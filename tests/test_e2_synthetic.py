import hashlib
import io
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_e2_synthetic import (audit_packed_payload, check_revision, distribution, fit_log_model,
    grouped_validation, log_contributions, pair_comparison, predict,
    repetition_summary, selected_task_jobs, write_union_tsv)


class E2SyntheticTests(unittest.TestCase):
    def test_revision_pins_reject_other_releases_and_invalid_prefixes(self):
        commit = "febd2bc1fecfa6619a86b19a7382e3fc4a3746af"
        for expected in (None, "febd2bc1", commit):
            check_revision(commit, expected)
        for expected in ("e2cecbb4", "", "f", "febd2bc1*", commit + "0"):
            with self.subTest(expected=expected), self.assertRaisesRegex(ValueError, "source revision mismatch"):
                check_revision(commit, expected)

    def test_repetition_ranges_units_and_order(self):
        jobs = [{"repetition": i, "time_wall_seconds": t, "peak_rss_kib": r * 2**20,
                 "status": "success", "warmup": False} for i, t, r in ((2, 30, 3), (0, 10, 1), (1, 20, 2))]
        s = repetition_summary(jobs)
        self.assertEqual(s["wall_seconds_repetitions"], [10, 20, 30])
        self.assertEqual(s["wall_seconds_median"], 20)
        self.assertEqual(s["wall_seconds_relative_range"], 1)
        self.assertEqual(s["rss_gib_repetitions"], [1, 2, 3])
        self.assertEqual(s["rss_gib_median"], 2)

    def test_reject_invalid_repetitions(self):
        jobs = [{"repetition": i, "time_wall_seconds": i + 1, "peak_rss_kib": 100,
                 "status": "success", "warmup": False} for i in range(3)]
        for mutation in ({"repetition": 1}, {"status": "timeout"}, {"warmup": True},
                         {"time_wall_seconds": 0}, {"time_wall_seconds": None}, {"peak_rss_kib": float("nan")}):
            with self.subTest(mutation=mutation), self.assertRaises(ValueError):
                repetition_summary([jobs[0], jobs[1], {**jobs[2], **mutation}])

    def test_exact_power_law_and_contribution_identity(self):
        rows = [{"n": n, "wall_seconds_median": 7 * (1 + n)**1.5} for n in (3, 10, 40, 100, 500)]
        m = fit_log_model(rows, "wall_seconds", ["n"], alpha=0)
        self.assertAlmostEqual(m["coefficients"]["n"], 1.5)
        self.assertAlmostEqual(m["intercept"], np.log(7))
        np.testing.assert_allclose(predict(m, rows), [r["wall_seconds_median"] for r in rows])
        c = log_contributions(m, rows[1], rows[-1])
        self.assertAlmostEqual(sum(c.values()), np.log(predict(m, [rows[-1]])[0] / predict(m, [rows[1]])[0]))

    def test_ridge_constant_and_correlated_predictors(self):
        rows = [{"n": n, "duplicate_n": n, "constant": 3, "rss_gib_median": 2 * (1 + n)**.8}
                for n in (3, 10, 40, 100, 500)]
        m = fit_log_model(rows, "rss_gib", ["n", "duplicate_n", "constant"], 1)
        self.assertEqual(m["coefficients"]["constant"], 0)
        self.assertAlmostEqual(m["coefficients"]["n"], m["coefficients"]["duplicate_n"])
        self.assertTrue(np.all(np.isfinite(predict(m, rows))))

    def test_group_validation_refits_without_geometry_leakage(self):
        rows = [{"n": n, "family": f, "seed": s, "wall_seconds_median": (1 + n)**(1 + .1 * s)}
                for n in (10, 20, 80) for f in ("a", "b") for s in (1, 2, 3)]
        for group_key in ("n", "family", "seed"):
            calls = []
            original = fit_log_model
            def tracking(train, outcome, features, alpha):
                calls.append(train)
                return original(train, outcome, features, alpha)
            with patch("analyze_e2_synthetic.fit_log_model", side_effect=tracking):
                cv, predictions = grouped_validation(rows, "wall_seconds", ["n"], 0, group_key)
            self.assertEqual(cv["pooled"]["n"], len(rows))
            for fold, train in zip(cv["folds"], calls):
                held = fold["held_out"]
                self.assertTrue(all(r[group_key] != held for r in train))
                self.assertEqual(fold["n"], sum(r[group_key] == held for r in rows))
                indices = [i for i, r in enumerate(rows) if r[group_key] == held]
                manual = predict(original(train, "wall_seconds", ["n"], 0), [rows[i] for i in indices])
                np.testing.assert_allclose(predictions[indices], manual)

    def test_packed_scalar_audit_and_tampering(self):
        j = {"id": "job", "input": "rna", "program": "LinearCParty", "mode": "b100", "warmup": False,
             "repetition": 0, "status": "success", "returncode": 0, "time_wall_seconds": 2.3,
             "peak_rss_kib": 1234, "command": ["cli", "-i", "rna"], "directory": "jobs/job"}
        payloads = {"jobs/job/result.json": json.dumps(j).encode(), "jobs/job/time.tsv": b"2.30\t1234\n",
                    "jobs/job/stdout.txt": b"example raw output"}
        j["artifacts"] = {name: {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
                          for name, data in payloads.items()}
        with tempfile.TemporaryDirectory(dir="/tmp") as tmp:
            p = Path(tmp) / "outputs.tar"
            with tarfile.open(p, "w") as archive:
                for name, data in payloads.items():
                    member = tarfile.TarInfo(name)
                    member.size = len(data)
                    archive.addfile(member, io.BytesIO(data))
            with patch.object(tarfile.TarFile, "extractall", side_effect=AssertionError("no extraction")), \
                 patch.object(tarfile.TarFile, "extract", side_effect=AssertionError("no extraction")):
                self.assertEqual(audit_packed_payload(p, [j]), 3)
            self.assertEqual(list(Path(tmp).iterdir()), [p])
            with self.assertRaisesRegex(ValueError, "result/checkpoint mismatch"):
                audit_packed_payload(p, [{**j, "time_wall_seconds": 2.4}])
            broken = {**j, "artifacts": {k: {**v, "sha256": "0" * 64} for k, v in j["artifacts"].items()}}
            with self.assertRaisesRegex(ValueError, "hash/size mismatch"):
                audit_packed_payload(p, [broken])

    def test_union_tsv_preserves_missing_and_structured_values(self):
        with tempfile.TemporaryDirectory(dir="/tmp") as tmp:
            p = Path(tmp) / "summary.tsv"
            write_union_tsv(p, [{"name": "a", "prediction": None}, {"name": "b", "values": [1, 2]}])
            lines = p.read_text().splitlines()
            self.assertEqual(lines[0], "name\tprediction\tvalues")
            self.assertEqual(lines[1], "a\tNA\t")
            self.assertEqual(lines[2], "b\t\t[1,2]")

    def test_distribution(self):
        d = distribution([1, 3, 2])
        self.assertEqual((d["count"], d["min"], d["median"], d["max"]), (3, 1, 2, 3))
        with self.assertRaises(ValueError):
            distribution([])

    def test_unexecuted_unselected_planner_placeholders(self):
        jobs = [{"input": "small", "directory": "tasks/small/1/run/jobs/a", "status": "success"},
                {"input": "large", "status": "not_started"}]
        self.assertEqual(selected_task_jobs(jobs, {"small"}, "tasks/small/1"), jobs[:1])
        with self.assertRaises(KeyError):
            selected_task_jobs(jobs, {"small", "large"}, "tasks/small/1")


if __name__ == "__main__":
    unittest.main()
