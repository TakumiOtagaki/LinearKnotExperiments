"""Exercise budget packing and real tiny CLI tasks without a scheduler or build."""

from contextlib import redirect_stderr, redirect_stdout
import copy
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from aggregate_array import aggregate
from plan_array import exclusive, load_plan, make_tasks, prepare, validate_array
from run_array_task import run_task
from run_cli import digest, validate, write_json
from validate_outputs import validate_outputs
from task_storage import ArchiveReader, jsonl_rows


class ArrayTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        quota = patch("file_budget.query_quota", return_value={"used_files": 23000000, "limit_files": 24000000})
        quota.start()
        self.addCleanup(quota.stop)
        self.software = self.root / "software"
        self.build = self.software / "build"
        self.build.mkdir(parents=True)
        (self.software / "energy.par").write_text("tiny fixture parameters\n")
        (self.build / "CMakeCache.txt").write_text("CMAKE_BUILD_TYPE:STRING=Release\n")
        self.events = self.root / "events.jsonl"
        # Append event timestamps under flock to measure actual process overlap.
        code = f'''#!{sys.executable}
import fcntl, json, pathlib, sys, time
def event(kind):
    with open({str(self.events)!r}, 'a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        stream.write(json.dumps({{'kind': kind, 'time': time.monotonic()}})+'\\n')
event('start')
time.sleep(0.04)
if '--bpp' in sys.argv:
    record = pathlib.Path(sys.argv[sys.argv.index('-i')+1]).read_text().splitlines()[0]
    pathlib.Path(sys.argv[sys.argv.index('--bpp')+1]).write_text(record+'\\n')
event('finish')
sys.exit(7 if pathlib.Path(sys.argv[0]).name == 'failure' else 0)
'''
        for name in ("LinearCParty", "failure"):
            binary = self.build / name
            binary.write_text(code)
            binary.chmod(0o755)
        self.lock_path = self.root / "lock.json"
        self.lock = {"schema_version": 1, "source_commit": "1" * 40,
                     "build_revision_verified": True,
                     "cmake_cache_sha256": digest(self.build / "CMakeCache.txt"),
                     "energy_parameter_sha256": digest(self.software / "energy.par"),
                     "binaries": {name: digest(self.build / name) for name in ("LinearCParty", "failure")}}
        write_json(self.lock_path, self.lock)
        self.inputs = []
        for name, length in (("a", 5), ("b", 6), ("c", 11)):
            path = self.root / (name + ".txt")
            path.write_text(f'>{name}\n' + 'A' * length + '\n' + '.' * length + '\n')
            self.inputs.append({"name": name, "path": str(path), "length": length})
        self.config = {"schema_version": 1, "purpose": "development",
                       "energy_parameter_file": "energy.par", "software_lock": str(self.lock_path),
                       "memory_limit_mib": 64, "timeout_seconds": 5, "warmups": 1, "repetitions": 2,
                       "inputs": [{k: v for k, v in i.items() if k != "length"} for i in self.inputs],
                       "modes": [{"name": "exact"}, {"name": "beam"}],
                       "programs": [{"name": "LinearCParty", "outputs": [{"flag": "--bpp", "name": "bpp.tsv"}]},
                                    {"name": "failure"}]}
        tier = {"name": "short", "min_length": 1, "max_length": 10, "resource": "cpu_4",
                "h_rt": "00:01:00", "target_seconds": 15, "workers": 1, "memory_budget_mib": 512,
                "cost_model": {"basis": "tiny independent fixture estimates",
                               "seconds_per_nt_condition": 0.1, "memory_base_mib": 16, "mib_per_nt": 0.1}}
        self.settings = {"schema_version": 1, "measurement": "accuracy", "seed": 31,
                         "tiers": [tier, {**copy.deepcopy(tier), "name": "long", "min_length": 11,
                                          "max_length": 20, "target_seconds": 20}]}
        self.config_path, self.settings_path = self.root / "config.json", self.root / "array.json"
        self.output, self.jobs = self.root / "run", self.root / "jobs"
        # Mock Git only via a test executable. Runner, hashes, copied files, locks,
        # environment parsing, and CLI process/memory control are real.
        fake_bin = self.root / "test-bin"
        fake_bin.mkdir()
        fake_git = fake_bin / "git"
        fake_git.write_text(f"#!{sys.executable}\nimport sys\n"
                            "if 'rev-parse' in sys.argv: print('1'*40)\n")
        fake_git.chmod(0o755)
        self.env = {**os.environ, "PATH": str(fake_bin) + os.pathsep + os.environ["PATH"],
                    "JOB_ID": "fixture-job", "SGE_TASK_ID": "1", "TMPDIR": str(self.root)}

    def prepare(self, dry_run=False):
        write_json(self.config_path, self.config)
        write_json(self.settings_path, self.settings)
        with patch.dict(os.environ, self.env), redirect_stdout(io.StringIO()):
            self.plan = prepare(self.config_path, self.settings_path, self.output, self.jobs,
                                None if dry_run else self.software, None if dry_run else self.build, dry_run)
        self.plan_path = self.output / "array_plan.json"
        return self.plan

    def run_task(self, tier="short", resume=False):
        with patch.dict(os.environ, self.env), redirect_stdout(io.StringIO()):
            return run_task(self.plan_path, tier, software=self.software, build=self.build, resume=resume)

    def test_packing_estimates_and_complete_matched_matrices(self):
        tasks = make_tasks(self.config, self.settings, self.inputs)
        self.assertEqual([[i["name"] for i in t["inputs"]] for t in tasks], [["a", "b"], ["c"]])
        self.assertAlmostEqual(tasks[0]["estimated_seconds"], (5 + 6) * .1 * 2 * 2 * 3)
        self.assertAlmostEqual(tasks[0]["estimated_memory_mib"], 256 + 16.6)
        self.settings["tiers"][0]["target_seconds"] = 8
        tasks = make_tasks(self.config, self.settings, self.inputs)
        self.assertEqual([t["task_id"] for t in tasks], [1, 2, 1])
        # Nonadjacent scaffold records for one sequence stay in the same task.
        variants = [{**self.inputs[0], "sequence_sha256": "same-rna"}, self.inputs[1],
                    {**self.inputs[0], "name": "a_variant", "sequence_sha256": "same-rna"}]
        self.settings["tiers"][0]["target_seconds"] = 13
        tasks = make_tasks(self.config, self.settings, variants)
        self.assertEqual([[i["name"] for i in t["inputs"]] for t in tasks], [["a", "a_variant"], ["b"]])
        model = self.settings["tiers"][0]["cost_model"]
        model.update(length_exponent=2, mode_factors={"exact": 2}, safety_factor=2)
        with self.assertRaisesRegex(ValueError, "one complete input matrix"):
            make_tasks(self.config, self.settings, self.inputs)

    def test_resource_memory_and_wall_time_guards(self):
        cases = [("workers", 5), ("resource", "cpu_2"), ("memory_budget_mib", 100000),
                 ("memory_budget_mib", 64), ("h_rt", "00:00:00"), ("h_rt", "25:00:00"),
                 ("target_seconds", 55), ("max_wall_fraction", .9)]
        for key, value in cases:
            settings = copy.deepcopy(self.settings)
            settings["tiers"][0][key] = value
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                validate_array(self.config, settings)
        settings = copy.deepcopy(self.settings)
        settings["measurement"] = "timing"
        with self.assertRaisesRegex(ValueError, "node_f"):
            validate_array(self.config, settings)
        for t in settings["tiers"]:
            t["resource"] = "node_f"
        validate_array(self.config, settings)
        settings["tiers"][0]["workers"] = 2
        with self.assertRaisesRegex(ValueError, "workers=1"):
            validate_array(self.config, settings)
        self.settings["tiers"][1]["min_length"] = 10
        with self.assertRaisesRegex(ValueError, "nonoverlapping"):
            validate_array(self.config, self.settings)

    def test_dry_run_freezes_inputs_lock_and_plan_without_execution(self):
        plan = self.prepare(dry_run=True)
        self.assertFalse(plan["dependencies_verified_at_planning"])
        self.assertFalse(self.events.exists())
        self.assertEqual(load_plan(self.plan_path), plan)
        self.assertEqual(len(plan["tasks"][0]["planned_jobs"]), 24)
        self.assertEqual((self.output / "original_config.json").read_bytes(), self.config_path.read_bytes())
        self.assertEqual((self.output / "software_lock.json").read_bytes(), self.lock_path.read_bytes())
        script = self.jobs / "short.sh"
        self.assertEqual(subprocess.run(["bash", "-n", str(script)], capture_output=True).returncode, 0)
        self.assertIn("#$ -t 1-1", script.read_text())
        self.assertIn("#$ -l cpu_4=1", script.read_text())
        with self.assertRaisesRegex(ValueError, "immutable"):
            self.prepare(dry_run=True)
        copied = self.output / "inputs.tar"
        copied.write_text("changed\n")
        with self.assertRaisesRegex(ValueError, "frozen dependency"):
            load_plan(self.plan_path)

    def test_environment_task_execution_resume_and_partial_aggregation(self):
        plan = self.prepare()
        initial = aggregate(self.plan_path)
        self.assertEqual(initial["task_status_counts"], {"not_started": 2})
        self.assertEqual(initial["missing_measured_jobs"], 24)
        self.assertEqual(initial["array_config"]["seed"], 31)
        self.assertIsNone(initial["started_utc"])
        # Real CLI invocation selects the task exclusively via SGE_TASK_ID.
        result = subprocess.run([sys.executable, str(ROOT / "scripts/run_array_task.py"),
                                 "--plan", str(self.plan_path), "--tier", "short"],
                                env=self.env, text=True, capture_output=True, timeout=20)
        self.assertEqual(result.returncode, 1, result.stderr)
        directory = self.output / plan["tasks"][0]["directory"]
        metadata = json.loads((directory / "task.json").read_text())
        attempt = metadata["attempts"][0]
        self.assertEqual(attempt["JOB_ID"], "fixture-job")
        self.assertEqual(attempt["SGE_TASK_ID"], "1")
        self.assertTrue(attempt["hostname"] and attempt["started_utc"] and attempt["finished_utc"])
        self.assertEqual(attempt["requested_resource"], "cpu_4=1")
        run = directory / "run"
        manifest = json.loads((run / "run.json").read_text())
        self.assertEqual([a["id"] for a in manifest["execution_order"]],
                         [j["id"] for j in plan["tasks"][0]["planned_jobs"]])
        checkpoint = run / "task_results.jsonl"
        kept = {j["id"]: j for j in jsonl_rows(checkpoint)[:-1]}
        pending = manifest["jobs"][-1]
        lines = checkpoint.read_bytes().splitlines(keepends=True)
        checkpoint.write_bytes(b"".join(lines[:-1]) + b'{"torn":')
        with (run / "task_outputs.tar").open("ab") as stream:
            stream.write(b"uncommitted archive tail")
        with self.assertRaisesRegex(ValueError, "use --resume"):
            self.run_task()
        self.assertEqual(self.run_task(resume=True), 1)
        resumed = json.loads((run / "run.json").read_text())
        self.assertEqual(len(resumed["execution_order"]), 25)
        for identifier, data in kept.items():
            self.assertEqual(next(j for j in jsonl_rows(checkpoint) if j["id"] == identifier), data)
        self.assertFalse((run / "jobs").exists())
        self.assertLessEqual(len(list(run.iterdir())), 5)
        merged = aggregate(self.plan_path)
        self.assertEqual(merged["task_status_counts"], {"failed": 1, "not_started": 1})
        self.assertEqual((merged["planned_measured_jobs"], merged["completed_measured_jobs"],
                          merged["missing_measured_jobs"]), (24, 16, 8))
        validation = validate_outputs(self.output)
        self.assertEqual(sum(r["validation_status"] == "job_missing" for r in validation["observations"]), 8)
        self.assertEqual(sum(r["validation_status"] == "success" for r in validation["observations"]), 8)
        summary = self.root / "summary.tsv"
        result = subprocess.run([sys.executable, str(ROOT / "scripts/summarize.py"), str(self.output),
                                 "--out", str(summary)], env=self.env, capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(json.loads(summary.with_suffix(".json").read_text())["missing_measured_jobs"], 8)
        self.assertEqual(self.run_task(tier="long"), 1)
        merged = aggregate(self.plan_path)
        self.assertEqual(merged["status"], "failed")
        self.assertEqual(merged["missing_measured_jobs"], 0)
        self.assertEqual(merged["task_status_counts"], {"failed": 2})

    def test_parallel_accuracy_respects_worker_bound_order_and_resume(self):
        self.settings["tiers"][0]["workers"] = 2
        self.prepare()
        self.assertEqual(self.run_task(), 1)
        events = [json.loads(line) for line in self.events.read_text().splitlines()]
        active = maximum = 0
        for event in sorted(events, key=lambda e: e["time"]):
            active += 1 if event["kind"] == "start" else -1
            maximum = max(maximum, active)
        self.assertEqual((active, maximum), (0, 2))
        directory = self.output / self.plan["tasks"][0]["directory"] / "run"
        manifest = json.loads((directory / "run.json").read_text())
        self.assertEqual([a["id"] for a in manifest["execution_order"]],
                         [j["id"] for j in manifest["planned_jobs"]])
        attempts = {a["id"]: a for a in manifest["execution_order"]}
        for job in manifest["planned_jobs"]:
            prefix = job["id"].rsplit("-", 1)[0]
            if not job["warmup"] and job["repetition"] == 0:
                self.assertLessEqual(attempts[prefix + "-warmup0"]["finished_utc"],
                                     attempts[job["id"]]["started_utc"])
            if not job["warmup"] and job["repetition"] == 1:
                self.assertLessEqual(max(a["finished_utc"] for j in manifest["planned_jobs"]
                                        if not j["warmup"] and j["repetition"] == 0
                                        for a in [attempts[j["id"]]]), attempts[job["id"]]["started_utc"])
        before = self.events.read_bytes()
        self.assertEqual(self.run_task(resume=True), 1)
        self.assertEqual(self.events.read_bytes(), before)

    def test_task_index_locks_and_rejected_dependency_resume(self):
        self.prepare()
        task_dir = self.output / self.plan["tasks"][0]["directory"]
        with exclusive(task_dir / ".execution.lock"):
            with self.assertRaisesRegex(ValueError, "holds"):
                self.run_task()
            merged = aggregate(self.plan_path)
            self.assertEqual(merged["task_status_counts"]["running"], 1)
        for identifier in ("0", "undefined", "2", "-1"):
            with patch.dict(os.environ, {**self.env, "SGE_TASK_ID": identifier}), self.assertRaises(ValueError):
                run_task(self.plan_path, "short")
        self.assertEqual(self.run_task(), 1)
        run_file = task_dir / "run/run.json"
        before = run_file.read_bytes()
        copied = task_dir / "run/config.json"
        copied.write_text(copied.read_text() + "\n")
        self.assertEqual(self.run_task(resume=True), 2)
        self.assertEqual(run_file.read_bytes(), before)

    def test_formal_requires_array_allocation_and_verified_build(self):
        self.config["purpose"] = "formal"
        with self.assertRaisesRegex(ValueError, "array runner"):
            validate(self.config)
        self.config["model_settings"] = {"build_revision_verified": True}
        for flag in (False, 1, "true"):
            with self.subTest(flag=flag), self.assertRaises(ValueError):
                validate({**self.config, "model_settings": {"build_revision_verified": flag}}, array=True)
        self.prepare()
        env = {k: v for k, v in self.env.items() if k not in ("JOB_ID", "SGE_TASK_ID")}
        with patch.dict(os.environ, env, clear=True), self.assertRaisesRegex(ValueError, "allocation"):
            run_task(self.plan_path, "short", task_id=1)
        self.assertEqual(self.run_task(), 1)

    def test_startup_failure_is_a_failed_task_with_missing_observations(self):
        self.prepare(dry_run=True)
        env = {k: v for k, v in self.env.items() if k not in ("LCP_SOFTWARE_ROOT", "LCP_BUILD_DIR")}
        with patch.dict(os.environ, env, clear=True), redirect_stderr(io.StringIO()):
            self.assertEqual(run_task(self.plan_path, "short"), 2)
        merged = aggregate(self.plan_path)
        self.assertEqual(merged["task_status_counts"], {"failed": 1, "not_started": 1})
        self.assertEqual(merged["missing_measured_jobs"], 24)
        measured = [j for j in merged["jobs"] if not j["warmup"]]
        self.assertEqual(sum(j["status"] == "task_failed" for j in measured), 16)
        self.assertEqual(sum(j["status"] == "not_started" for j in measured), 8)

    def test_sigterm_records_interruption_and_stops_the_cli_group(self):
        self.config.update(inputs=self.config["inputs"][:1], programs=self.config["programs"][:1],
                           warmups=0, repetitions=1)
        binary = self.build / "LinearCParty"
        binary.write_text(binary.read_text().replace("time.sleep(0.04)", "time.sleep(20)"))
        self.lock["binaries"]["LinearCParty"] = digest(binary)
        write_json(self.lock_path, self.lock)
        self.prepare()
        process = subprocess.Popen([sys.executable, str(ROOT / "scripts/run_array_task.py"),
                                    "--plan", str(self.plan_path), "--tier", "short"],
                                   env=self.env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        try:
            deadline = time.monotonic() + 5
            while not self.events.exists() and time.monotonic() < deadline:
                time.sleep(.01)
            self.assertTrue(self.events.exists())
            process.terminate()
            stdout, stderr = process.communicate(timeout=5)
            self.assertEqual(process.returncode, 130, stderr)
        finally:
            if process.poll() is None:
                process.kill()
                process.communicate()
        metadata = json.loads((self.output / self.plan["tasks"][0]["directory"] / "task.json").read_text())
        self.assertEqual(metadata["status"], "interrupted")
        self.assertIsNotNone(metadata["attempts"][0]["finished_utc"])
        merged = aggregate(self.plan_path)
        self.assertEqual(merged["status"], "interrupted")
        self.assertEqual(merged["missing_measured_jobs"], 1)


if __name__ == "__main__":
    unittest.main()
