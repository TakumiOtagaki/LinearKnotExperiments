"""Tiny packed/legacy parity fixtures, durable checkpoints, and inode guards."""

from contextlib import redirect_stdout
import copy
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from analyze_accuracy import analyze as accuracy
from analyze_e0 import POLICIES, analysis_settings, analyze as e0
from file_budget import check_budget, estimate_files, query_quota
from run_cli import digest, execute, write_json
from task_storage import ArchiveReader, RunReader, TaskStore, jsonl_rows, local_workspace
from plan_array import prepare
from aggregate_array import aggregate
import test_array


class StorageTests(unittest.TestCase):
    def test_durable_commit_torn_checkpoint_orphan_tar_and_hashes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            raw = root / "stdout.txt"; raw.write_text("preserved output\n")
            store = TaskStore(root)
            row = store.commit({"id": "a"}, {"jobs/a/stdout.txt": raw})
            first_end = row["archive_end"]
            with store.checkpoint.open("ab") as stream: stream.write(b'{"id":"torn')
            with store.archive.open("ab") as stream: stream.write(b"orphan")
            recovered = TaskStore(root, resume=True)
            self.assertEqual(recovered.rows, [row])
            self.assertEqual(recovered.archive.stat().st_size, first_end + 1024)
            self.assertEqual(ArchiveReader(store.archive).read("jobs/a/stdout.txt"), raw.read_bytes())
            recovered.commit({"id": "b"}, {"jobs/b/stdout.txt": raw})
            self.assertEqual(len(jsonl_rows(store.checkpoint)), 2)
            self.assertEqual(ArchiveReader(store.archive).read("jobs/b/stdout.txt"), raw.read_bytes())
            self.assertEqual(row["stdout_summary"]["bytes"], len(raw.read_bytes()))
            with self.assertRaisesRegex(ValueError, "already checkpointed"):
                recovered.commit({"id": "a"}, {"stdout.txt": raw})
            store.checkpoint.write_text('{"broken":}\n')
            with self.assertRaisesRegex(ValueError, "malformed checkpoint"):
                TaskStore(root, resume=True)

    def test_local_scratch_is_required_and_cleaned(self):
        with patch.dict(os.environ, {}, clear=True):
            with self.assertRaisesRegex(ValueError, "TMPDIR"):
                with local_workspace(): pass
            with self.assertRaisesRegex(ValueError, "node-local"):
                with local_workspace(ROOT): pass
            with tempfile.TemporaryDirectory() as temporary:
                with local_workspace(temporary) as workspace:
                    (workspace / "tiny.txt").write_text("fixture")
                self.assertFalse(workspace.exists())

    def test_quota_parser_fallback_fraction_and_conservative_cap(self):
        output = "Filesystem kbytes quota limit grace files quota limit grace\n/gs/bs 123 456 456 - 23117996 28000000 28000000 -\n"
        with patch("file_budget.subprocess.run", return_value=subprocess.CompletedProcess([], 0, output, "")):
            quota = query_quota(102571, "/gs/bs")
            self.assertEqual(quota["used_files"], 23117996)
            report = check_budget(estimate_files(841, 8), {})
            self.assertEqual(report["remaining_files"], 882004)
            self.assertEqual(report["estimated_shared_entries"], 10963)
        with patch("file_budget.query_quota", return_value=None):
            with self.assertRaisesRegex(ValueError, "fallback_remaining_files"):
                check_budget(estimate_files(1, 1), {})
            with self.assertRaisesRegex(ValueError, "budget refused"):
                check_budget(estimate_files(10, 1), {"file_budget": {"fallback_remaining_files": 100}})
            report = check_budget(estimate_files(1, 1), {"file_budget": {"fallback_remaining_files": 10000}})
            self.assertIn("config", report["remaining_source"])


class PackedParityTests(unittest.TestCase):
    def setUp(self):
        self.fixture = test_array.ArrayTests()
        self.fixture.setUp(); self.addCleanup(self.fixture.doCleanups)
        self.f = self.fixture
        self.f.config.update(inputs=self.f.config["inputs"][:1], warmups=0, repetitions=1)

    def programs(self, e0_mode=False):
        names = list(POLICIES) if e0_mode else ["LinearCParty", "LinearPKCapR", "LinearPKRaccess"]
        code = f'''#!{sys.executable}
import json,pathlib,sys
name,sequence,scaffold=pathlib.Path(sys.argv[sys.argv.index('-i')+1]).read_text().splitlines()
program=pathlib.Path(sys.argv[0]).name
policies={POLICIES!r}
if program in policies:
    z=.9 if program=='E0Default' else 1.
    print(json.dumps(dict(policy=policies[program],log_z=z,free_energy=-.6*z,rt=.6)))
elif program=='LinearCParty':
    pathlib.Path(sys.argv[sys.argv.index('--bpp')+1]).write_text(name+'\\n')
elif program=='LinearPKCapR':
    print(name)
    for c in ('Stem','Hairpin','Bulge','Internal','Multiloop','Exterior','PseudoLoop'):
        print(c+' '+' '.join(['1' if c=='Exterior' else '0']*len(sequence)))
else:
    print(name+'\\n# output_mode=probability\\n# access_lens=1\\n# scaffold='+scaffold+'\\n# folding_rt=.6\\n# ensemble_free_energy=-.6')
    for i in range(len(sequence)): print(str(i)+'\\t1,1;')
'''
        self.f.config["programs"] = []
        for name in names:
            binary = self.f.build / name; binary.write_text(code); binary.chmod(0o755)
            self.f.lock["binaries"][name] = digest(binary)
            program = {"name": name}
            if name == "LinearCParty": program["outputs"] = [{"flag": "--bpp", "name": "bpp.tsv"}]
            self.f.config["programs"].append(program)
        write_json(self.f.lock_path, self.f.lock)

    def compare_formats(self, analyze):
        self.f.prepare()
        self.assertEqual(self.f.run_task(), 0)
        packed_path = self.f.output / self.f.plan["tasks"][0]["directory"] / "run"
        old = self.f.root / "legacy"
        with patch.dict(os.environ, self.f.env), redirect_stdout(io.StringIO()):
            self.assertEqual(execute(self.f.config_path, self.f.software, self.f.build, old), 0)
        old_result, packed_result = analyze(old), analyze(packed_path)
        self.assertEqual(old_result["observations"], packed_result["observations"])
        archive = self.f.root / "legacy.tar.gz"
        with tarfile.open(archive, "w:gz") as stream: stream.add(old, arcname="legacy")
        shutil.rmtree(old)
        archived_result = analyze(archive)
        self.assertEqual(old_result["observations"], archived_result["observations"])
        for name, source in (("packed", packed_path), ("archive", archive)):
            result = subprocess.run([sys.executable, str(ROOT / "scripts/summarize.py"), str(source),
                                     "--out", str(self.f.root / (name + ".tsv"))],
                                    capture_output=True, text=True, timeout=10)
            self.assertEqual(result.returncode, 0, result.stderr)
        packed_summary = json.loads((self.f.root / "packed.json").read_text())
        archive_summary = json.loads((self.f.root / "archive.json").read_text())
        self.assertEqual(packed_summary["status_counts"], archive_summary["status_counts"])
        self.assertEqual(packed_summary["measured_jobs"], archive_summary["measured_jobs"])
        self.assertFalse(old.exists())
        self.assertFalse((packed_path / "jobs").exists())
        self.assertLessEqual(len(list(packed_path.iterdir())), 5)
        self.assertEqual(self.f.run_task(resume=True), 0)
        self.assertEqual(analyze(packed_path)["observations"], packed_result["observations"])

    def test_e0_values_counts_and_failed_denominators_match_legacy_and_archive(self):
        self.programs(e0_mode=True)
        self.f.config.update(modes=[{"name": "exact"}], analysis=analysis_settings())
        self.compare_formats(e0)

    def test_accuracy_bpp_context_accessibility_match_legacy_and_archive(self):
        self.programs()
        # Resource timings vary across actual processes; compare scientific outputs.
        def scientific(run):
            value = accuracy(run, "exact")
            for row in value["observations"]:
                row.pop("time_wall_seconds", None); row.pop("peak_rss_kib", None)
            value["observations"].sort(key=lambda row: row["id"])
            return value
        self.compare_formats(scientific)

    def test_plan_refusal_precedes_writes_and_estimate_creates_one_file(self):
        self.f.settings["file_budget"] = {"fallback_remaining_files": 1}
        write_json(self.f.config_path, self.f.config); write_json(self.f.settings_path, self.f.settings)
        with patch("file_budget.query_quota", return_value=None):
            with self.assertRaisesRegex(ValueError, "budget refused"):
                prepare(self.f.config_path, self.f.settings_path, self.f.output, self.f.jobs, dry_run=True)
        self.assertFalse(self.f.output.exists()); self.assertFalse(self.f.jobs.exists())
        self.f.settings.pop("file_budget"); write_json(self.f.settings_path, self.f.settings)
        report = prepare(self.f.config_path, self.f.settings_path, self.f.output, self.f.jobs,
                         dry_run=True, estimate_only=True)
        self.assertFalse(report["executable"])
        self.assertEqual(report["planned_jobs"], 4)
        self.assertEqual([p.name for p in self.f.output.iterdir()], ["array_estimate.json"])
        self.assertFalse(self.f.jobs.exists())

    def test_archived_array_recovers_markers_and_missing_jobs_without_extraction(self):
        self.f.prepare()
        old_root = self.f.root / "old-array"
        task = self.f.plan["tasks"][0]
        task_dir = old_root / task["directory"]
        task_dir.mkdir(parents=True)
        old_run = task_dir / "run"
        with patch.dict(os.environ, self.f.env), redirect_stdout(io.StringIO()):
            self.assertEqual(execute(self.f.config_path, self.f.software, self.f.build, old_run), 1)
        saved = json.loads((old_run / "run.json").read_text())
        missing = saved["planned_jobs"][-1]
        (old_run / "jobs" / missing["id"] / "result.json").unlink()
        (old_root / "inputs").mkdir()
        shutil.copyfile(old_run / "inputs/a.txt", old_root / "inputs/a.txt")
        metadata = {"status": "running", "attempts": [], "tier": task["tier"], "task_id": task["task_id"]}
        write_json(task_dir / "task.json", metadata)
        write_json(old_root / "array_plan.json", {**self.f.plan, "tasks": [task], "output_dir": str(old_root)})
        saved.update(run_kind="array", summary_only=True, jobs=[], inputs=[{
            "name": "a", "copied_path": str(old_root / "inputs/a.txt"), "sha256": digest(old_root / "inputs/a.txt")}],
            planned_tasks=1, task_status_counts={"running": 1})
        write_json(old_root / "run.json", saved)
        archive = self.f.root / "old-array.tar"
        with tarfile.open(archive, "w") as stream: stream.add(old_root, arcname="old-array")
        shutil.rmtree(old_root)
        value = aggregate(archive)
        self.assertEqual((value["completed_measured_jobs"], value["missing_measured_jobs"]), (3, 1))
        self.assertEqual(value["task_status_counts"], {"interrupted": 1})
        reader = RunReader(archive)
        completed = next(j for j in reader.manifest["jobs"] if j["result_present"])
        reader.artifact(completed, "stdout.txt")
        self.assertFalse(old_root.exists())
        selected = RunReader(str(archive) + "::" + task["directory"] + "/run")
        self.assertEqual(len(selected.manifest["jobs"]), 4)

    def test_shared_write_failure_stops_workers_and_keeps_jobs_pending(self):
        self.f.settings["tiers"][0]["workers"] = 2
        self.f.prepare()
        with patch("task_storage.TaskStore.commit", side_effect=OSError("Disk quota exceeded")):
            self.assertEqual(self.f.run_task(), 2)
        value = aggregate(self.f.plan_path)
        self.assertEqual(value["missing_measured_jobs"], 4)
        self.assertEqual(value["task_status_counts"], {"failed": 1})
        self.assertEqual(list(self.f.root.glob("lcp-task-*")), [])
        events = [json.loads(line) for line in self.f.events.read_text().splitlines()]
        self.assertLessEqual(sum(e["kind"] == "start" for e in events), 2)


if __name__ == "__main__": unittest.main()
