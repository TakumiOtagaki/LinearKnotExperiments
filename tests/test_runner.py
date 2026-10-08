"""Exercise failure retention, timeout cleanup and immutable run inputs."""

import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "scripts/run_cli.py"
SUMMARY = ROOT / "scripts/summarize.py"


class RunnerTests(unittest.TestCase):
    def test_failures_timeouts_and_warmups_remain_visible(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            software = root / "software"
            build = software / "build"
            build.mkdir(parents=True)
            (software / "params").mkdir()
            (software / "params/model.par").write_text("synthetic parameter fixture\n")
            record = root / "record.txt"
            record.write_text(">case\nAAAAA\n.....\n")
            success = build / "success"
            success.write_text("#!/bin/sh\nprintf 'output\\n'\n")
            success.chmod(0o755)
            failure = build / "failure"
            failure.write_text("#!/bin/sh\nprintf 'expected failure\\n' >&2\nexit 7\n")
            failure.chmod(0o755)
            timeout = build / "timeout"
            timeout.write_text("#!/bin/sh\nsleep 2\nprintf 'late output\\n'\n")
            timeout.chmod(0o755)
            config = {"schema_version": 1, "purpose": "development",
                      "energy_parameter_file": "params/model.par", "timeout_seconds": 0.1,
                      "warmups": 1, "repetitions": 1,
                      "inputs": [{"name": "fixture", "path": str(record)}],
                      "modes": [{"name": "test", "args": []}],
                      "programs": [{"name": name, "args": []}
                                   for name in ("failure", "timeout", "success")]}
            config_path = root / "config.json"
            config_path.write_text(json.dumps(config))
            output = root / "run"
            command = [sys.executable, str(RUNNER), "--config", str(config_path),
                       "--software-root", str(software), "--out", str(output)]
            result = subprocess.run(command, capture_output=True, text=True, timeout=20)
            self.assertEqual(result.returncode, 1, result.stderr)
            manifest = json.loads((output / "run.json").read_text())
            self.assertEqual(len(manifest["jobs"]), 6)
            measured = [job for job in manifest["jobs"] if not job["warmup"]]
            self.assertEqual([job["status"] for job in measured],
                             ["failed", "timeout", "success"])
            self.assertEqual(measured[0]["returncode"], 7)
            for job in manifest["jobs"]:
                directory = output / job["directory"]
                self.assertTrue((directory / "command.json").is_file())
                self.assertTrue((directory / "stdout.txt").is_file())
                self.assertTrue((directory / "stderr.txt").is_file())
                self.assertTrue((directory / "result.json").is_file())
            self.assertGreater(measured[2]["peak_rss_kib"], 0)
            record.write_text("changed after execution\n")
            frozen = Path(manifest["inputs"][0]["copied_path"])
            self.assertEqual(frozen.read_text(), ">case\nAAAAA\n.....\n")
            summary_path = root / "summary.tsv"
            summary = subprocess.run([sys.executable, str(SUMMARY), str(output),
                                      "--out", str(summary_path)], capture_output=True,
                                     text=True, timeout=10)
            self.assertEqual(summary.returncode, 0, summary.stderr)
            counts = json.loads(summary_path.with_suffix(".json").read_text())
            self.assertEqual(counts["measured_jobs"], 3)
            self.assertEqual(counts["status_counts"],
                             {"failed": 1, "timeout": 1, "success": 1})
            with summary_path.open() as stream:
                rows = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(len(rows), 3)
            before = (output / "run.json").read_bytes()
            repeat = subprocess.run(command, capture_output=True, text=True, timeout=10)
            self.assertEqual(repeat.returncode, 2)
            self.assertEqual((output / "run.json").read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
