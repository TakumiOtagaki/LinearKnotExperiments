"""PF-only parity, roundoff limits, closure attribution and failed denominators."""

import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from analyze_e0 import POLICIES, analyze, analysis_settings, compare, read_observation, within_tolerance
from run_cli import digest, write_json


class E0AnalysisTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.run = Path(self.temp.name); self.inputs = []; self.jobs = []

    def add(self, name, compatible=1.0, default=1.0, oracle=1.0, failed=None):
        record = self.run / (name + ".txt"); record.write_text(">" + name + "\nGAAAAAAC\nX.......\n")
        self.inputs.append({"name": name, "copied_path": str(record), "sha256": digest(record)})
        for program, z in zip(POLICIES, (compatible, default, oracle)):
            directory = self.run / name / program; directory.mkdir(parents=True)
            self.jobs.append({"input": name, "program": program, "status": "timeout" if program == failed else "success",
                              "warmup": False, "directory": str(directory.relative_to(self.run))})
            # Excluded historical fields are deliberately inconsistent/nonfinite.
            write_json(directory / "stdout.txt", {"policy": POLICIES[program], "log_z": z, "rt": .6,
                       "free_energy": -.6 * z, "mfe": float("nan"), "pairs": [[1, 8, 7.0]], "contexts": []})

    def save(self, status="success"):
        write_json(self.run / "run.json", {"status": status, "config": {"analysis": analysis_settings()},
                   "inputs": self.inputs, "jobs": self.jobs})

    def test_pf_only_parity_and_closure_counts_keep_failed_inputs(self):
        self.add("same"); self.add("closure", default=.9)
        self.add("bad_baseline", compatible=.8, default=.7)
        self.add("timeout", failed="E0CPartyOracle"); self.save()
        result = analyze(self.run)
        self.assertEqual(result["planned_inputs"], 4)
        self.assertEqual(result["numerical_parity_counts"], {"success": 2, "parity_failed": 1, "observation_unavailable": 1})
        self.assertEqual(result["closure_counts"], {"planned_inputs": 4, "attribution_valid": 2, "effect_detected": 1,
                         "within_roundoff": 1, "unassessable": 2, "restriction_direction_inconsistent": 0})
        for row in result["observations"]:
            if "metrics" in row: self.assertEqual(set(row["metrics"]), {"delta_log_z", "delta_free_energy"})

    def test_both_absolute_roundoff_limits_are_required(self):
        self.assertTrue(within_tolerance({"delta_log_z": 1e-9, "delta_free_energy": 1e-10}, analysis_settings()))
        self.assertFalse(within_tolerance({"delta_log_z": 4e-10, "delta_free_energy": 2.4e-10}, analysis_settings()))
        self.add("strict_energy", compatible=1 + 4e-10); self.save()
        self.assertEqual(analyze(self.run)["numerical_parity_counts"], {"parity_failed": 1})

    def test_bad_policy_and_malformed_metadata_are_artifact_errors(self):
        self.add("policy"); self.add("rt")
        path = self.run / "policy/E0LegacyCompatible/stdout.txt"
        v = json.loads(path.read_text()); v["policy"] = POLICIES["E0Default"]; write_json(path, v)
        path = self.run / "rt/E0CPartyOracle/stdout.txt"
        v = json.loads(path.read_text()); v["rt"] = .7; write_json(path, v)
        self.save(); result = analyze(self.run)
        self.assertEqual(result["numerical_parity_counts"], {"artifact_error": 2})
        self.assertEqual(result["closure_counts"]["unassessable"], 2)

    def test_zero_pf_values_and_nonfinite_rejection(self):
        path = self.run / "observation.json"
        write_json(path, {"log_z": 0, "free_energy": 0, "rt": .6})
        self.assertEqual(read_observation(path)["log_z"], 0)
        for key, value in (("log_z", float("nan")), ("free_energy", float("inf")), ("rt", 0)):
            record = {"log_z": 0, "free_energy": 0, "rt": .6}; record[key] = value; write_json(path, record)
            with self.assertRaises(ValueError): read_observation(path)
        with self.assertRaisesRegex(ValueError, "RT differs"):
            compare({"rt": .6}, {"rt": .7})

    def test_changed_input_running_and_duplicate_jobs_are_rejected(self):
        self.add("rna"); self.save("running")
        with self.assertRaisesRegex(ValueError, "finish"): analyze(self.run)
        self.save(); (self.run / "rna.txt").write_text(">rna\nAAAAAAAA\n........\n")
        with self.assertRaisesRegex(ValueError, "input changed"): analyze(self.run)
        self.inputs[0]["sha256"] = digest(self.run / "rna.txt"); self.jobs.append(self.jobs[0]); self.save()
        with self.assertRaisesRegex(ValueError, "one exact observation"): analyze(self.run)

    def test_restriction_direction_violation_is_counted(self):
        self.add("wrong_direction", default=1.1); self.save()
        result = analyze(self.run)
        self.assertEqual(result["numerical_parity_counts"], {"success": 1})
        self.assertEqual(result["closure_counts"]["restriction_direction_inconsistent"], 1)


if __name__ == "__main__": unittest.main()
