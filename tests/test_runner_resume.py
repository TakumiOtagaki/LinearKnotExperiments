"""Check repetition blocks, dependency checks and interrupted-run recovery."""

from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from run_cli import digest, main, plan_jobs, validate, write_json


class OrderTests(unittest.TestCase):
    def setUp(self):
        self.config = {"schema_version": 1, "purpose": "development",
                       "inputs": [{"name": name} for name in ("rna1", "rna2")],
                       "modes": [{"name": name} for name in ("exact", "beam")],
                       "programs": [{"name": name} for name in ("a", "b", "c")],
                       "warmups": 2, "repetitions": 3}

    def test_default_keeps_fixed_nested_order(self):
        plan = plan_jobs(self.config)
        expected = [f"{rna}-{mode}-{program}-{kind}{rep}"
                    for rna in ("rna1", "rna2") for mode in ("exact", "beam")
                    for program in ("a", "b", "c")
                    for kind, count in (("warmup", 2), ("rep", 3)) for rep in range(count)]
        self.assertEqual([job["id"] for job in plan], expected)
        self.config["order"] = {"kind": "fixed"}
        self.assertEqual(plan_jobs(self.config), plan)

    def test_seeded_repetition_blocks_and_first_measurement_warmups(self):
        self.config["order"] = {"kind": "blocked_random", "seed": 20261011}
        validate(self.config)
        plan = plan_jobs(self.config)
        self.assertEqual(plan_jobs(self.config), plan)
        combos = [(rna, mode, program) for rna in ("rna1", "rna2")
                  for mode in ("exact", "beam") for program in ("a", "b", "c")]
        rng = random.Random(20261011)
        for rep in range(3):
            expected = combos.copy()
            rng.shuffle(expected)
            block = [job for job in plan if not job["warmup"] and job["repetition"] == rep]
            self.assertEqual([(j["input"], j["mode"], j["program"]) for j in block], expected)
            for job in block:
                index = plan.index(job)
                if rep == 0:
                    self.assertEqual([j["repetition"] for j in plan[index-2:index]], [0, 1])
                    for warmup in plan[index-2:index]:
                        self.assertTrue(warmup["warmup"])
                        self.assertEqual(tuple(warmup[k] for k in ("input", "mode", "program")),
                                         tuple(job[k] for k in ("input", "mode", "program")))
        measured = [j["repetition"] for j in plan if not j["warmup"]]
        self.assertEqual(measured, [0]*12 + [1]*12 + [2]*12)
        self.assertEqual(len({j["id"] for j in plan}), 60)
        self.config["order"]["seed"] += 1
        self.assertNotEqual(plan_jobs(self.config), plan)

    def test_random_order_requires_an_integer_seed(self):
        for order in ({"kind": "unknown"}, {"kind": "blocked_random"},
                      {"kind": "blocked_random", "seed": True},
                      {"kind": "blocked_random", "seed": "1"}):
            with self.subTest(order=order), self.assertRaises(ValueError):
                validate({**self.config, "order": order})


class ResumeTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.software = self.root / "software"
        self.build = self.software / "build"
        self.build.mkdir(parents=True)
        self.parameter = self.software / "energy.par"
        self.parameter.write_text("fixture parameters\n")
        self.cache = self.build / "CMakeCache.txt"
        self.cache.write_text("CMAKE_BUILD_TYPE:STRING=Release\n")
        self.record = self.root / "input.txt"
        self.record.write_text(">rna\nAAAAA\n.....\n")
        for name, code in (("success", 0), ("failure", 7)):
            binary = self.build / name
            binary.write_text(f"#!/bin/sh\nprintf '{name}\\n'\nexit {code}\n")
            binary.chmod(0o755)
        self.lock_path = self.root / "lock.json"
        self.lock = {"schema_version": 1, "source_commit": "fixed-test-source",
                     "cmake_cache_sha256": digest(self.cache),
                     "energy_parameter_sha256": digest(self.parameter),
                     "binaries": {name: digest(self.build / name) for name in ("success", "failure")}}
        write_json(self.lock_path, self.lock)
        self.config_path = self.root / "config.json"
        self.config = {"schema_version": 1, "purpose": "development",
                       "energy_parameter_file": "energy.par",
                       "software_lock": str(self.lock_path), "timeout_seconds": 5,
                       "warmups": 1, "repetitions": 2,
                       "order": {"kind": "blocked_random", "seed": 20261011},
                       "inputs": [{"name": "rna", "path": str(self.record)}],
                       "modes": [{"name": "a"}, {"name": "b"}],
                       "programs": [{"name": "success"}, {"name": "failure"}]}
        write_json(self.config_path, self.config)
        self.output = self.root / "run"
        self.state = {"commit": "fixed-test-source", "dirty": False,
                      "status": "", "tracked_diff_sha256": "fixture"}

    def run_runner(self, resume=False):
        argv = ["run_cli.py", "--config", str(self.config_path),
                "--software-root", str(self.software), "--out", str(self.output)]
        if resume:
            argv.append("--resume")
        # Mock provenance only; dependency hashing, process execution and CLI parsing are real.
        with patch("sys.argv", argv), patch("run_cli.git_state", return_value=self.state), \
                redirect_stdout(io.StringIO()):
            return main()

    def manifest(self):
        return json.loads((self.output / "run.json").read_text())

    def test_resumes_only_missing_results_and_recovers_unlisted_results(self):
        self.assertEqual(self.run_runner(), 1)
        saved = self.manifest()
        self.assertEqual([a["id"] for a in saved["execution_order"]],
                         [j["id"] for j in saved["planned_jobs"]])
        pending = saved["planned_jobs"][-2:]
        pending_ids = {j["id"] for j in pending}
        kept = {j["id"]: (self.output / j["directory"] / "result.json").read_bytes()
                for j in saved["jobs"] if j["id"] not in pending_ids}
        for job in pending:
            directory = self.output / "jobs" / job["id"]
            (directory / "result.json").unlink()
            (directory / "partial.txt").write_text("unfinished output\n")
            (directory / "time.tsv").write_text("9999\t9999\n")
        # Emulate interruption after a result write but before its manifest update.
        saved["jobs"] = saved["jobs"][1:]
        saved["status"] = "running"
        for key in ("finished_utc", "software_changed_during_run", "binaries_changed_during_run"):
            saved.pop(key)
        write_json(self.output / "run.json", saved)
        self.assertEqual(self.run_runner(resume=True), 1)
        resumed = self.manifest()
        self.assertEqual(len(resumed["jobs"]), 12)
        self.assertEqual(len({j["id"] for j in resumed["jobs"]}), 12)
        self.assertEqual([a["id"] for a in resumed["execution_order"][-2:]],
                         [j["id"] for j in pending])
        self.assertEqual(len(resumed["execution_order"]), 14)
        self.assertEqual(resumed["resume_history"][0]["pending_jobs"], 2)
        self.assertEqual(resumed["resume_history"][0]["completed_jobs"], 10)
        self.assertEqual(resumed["resumed_utc"], resumed["resume_history"][0]["resumed_utc"])
        for job_id, data in kept.items():
            self.assertEqual((self.output / "jobs" / job_id / "result.json").read_bytes(), data)
        for job in pending:
            archived = self.output / "incomplete_attempts" / job["id"] / "0"
            self.assertEqual((archived / "partial.txt").read_text(), "unfinished output\n")
            result = json.loads((self.output / "jobs" / job["id"] / "result.json").read_text())
            self.assertLess(result["time_wall_seconds"], 5)
        self.assertEqual(resumed["status"], "failed")  # Completed failures remain observations.
        self.assertEqual(self.run_runner(resume=True), 1)
        self.assertEqual(len(self.manifest()["execution_order"]), 14)
        self.assertEqual(len(self.manifest()["resume_history"]), 2)

    def test_dependency_mismatches_refuse_without_modifying_the_run(self):
        self.run_runner()
        manifest_path = self.output / "run.json"
        before = manifest_path.read_bytes()
        cases = [(self.config_path, b" "), (self.lock_path, b" "),
                 (self.record, b"A"), (self.output / "inputs/rna.txt", b"A"),
                 (self.parameter, b"changed"), (self.output / "inputs/energy.par", b"changed"),
                 (self.build / "success", b"changed"), (self.cache, b"changed"),
                 (self.output / "config.json", b" ")]
        for path, extra in cases:
            with self.subTest(path=path):
                original = path.read_bytes()
                path.write_bytes(original + extra)
                try:
                    with self.assertRaises(ValueError):
                        self.run_runner(resume=True)
                    self.assertEqual(manifest_path.read_bytes(), before)
                finally:
                    path.write_bytes(original)
        self.state = {**self.state, "commit": "different-source"}
        with self.assertRaises(ValueError):
            self.run_runner(resume=True)
        self.assertEqual(manifest_path.read_bytes(), before)

    def test_legacy_fixed_manifest_resumes_without_order_or_finalization_fields(self):
        self.config.pop("order")
        write_json(self.config_path, self.config)
        self.run_runner()
        saved = self.manifest()
        pending = saved["jobs"].pop()
        (self.output / pending["directory"] / "result.json").unlink()
        for key in ("planned_jobs", "execution_order", "cmake_cache_sha256", "finished_utc",
                    "software_changed_during_run", "binaries_changed_during_run"):
            saved.pop(key)
        saved["status"] = "running"
        write_json(self.output / "run.json", saved)
        self.assertEqual(self.run_runner(resume=True), 1)
        resumed = self.manifest()
        self.assertEqual(len(resumed["jobs"]), 12)
        self.assertTrue(resumed["execution_order"][0]["recovered_legacy"])
        self.assertEqual(resumed["execution_order"][-1]["id"], pending["id"])

    def test_interrupted_summary_keeps_jobs_and_marks_unchecked_fields_unknown(self):
        self.run_runner()
        saved = self.manifest()
        saved["status"] = "running"
        for key in ("finished_utc", "software_changed_during_run", "binaries_changed_during_run"):
            saved.pop(key)
        write_json(self.output / "run.json", saved)
        output = self.root / "summary.tsv"
        result = subprocess.run([sys.executable, str(ROOT / "scripts/summarize.py"),
                                 str(self.output), "--out", str(output)],
                                capture_output=True, text=True, timeout=10)
        self.assertEqual(result.returncode, 0, result.stderr)
        summary = json.loads(output.with_suffix(".json").read_text())
        self.assertEqual(summary["run_status"], "interrupted")
        self.assertEqual(summary["measured_jobs"], 8)
        self.assertEqual(summary["status_counts"], {"success": 4, "failed": 4})
        self.assertIsNone(summary["software_changed_during_run"])
        self.assertIsNone(summary["binaries_changed_during_run"])


if __name__ == "__main__":
    unittest.main()
