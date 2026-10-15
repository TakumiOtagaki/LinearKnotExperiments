"""Check cohort matching, closure denominators, and finite-range slope units."""

import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_e0_e2_runs import FULL_FEATURES, empirical_slopes, historical_comparisons, summarize_e0


class E0E2RunTests(unittest.TestCase):
    def test_slopes_keep_each_seed_and_series_separate(self):
        rows = [{"program": "LinearCParty", "series": series, "seed": seed, "n": n,
                 "wall_seconds_median": 7 * n**power, "rss_gib_median": .03 * n**.5}
                for series, power in (("empty", 1.2), ("stem", .8)) for seed in (1, 2, 3)
                for n in (250, 500, 1000)]
        slopes = empirical_slopes(rows)
        self.assertEqual(len(slopes), 12)
        for row in slopes:
            expected = .5 if row["outcome"] == "rss_gib" else 1.2 if row["series"] == "empty" else .8
            self.assertAlmostEqual(row["slope_ln_y_vs_ln_n"], expected)
            self.assertAlmostEqual(row["log_r2"], 1)

    def test_history_matches_by_identity_and_rejects_changed_scaffold(self):
        a = {**dict.fromkeys(FULL_FEATURES, 0), "run": "old", "name": "rna", "program": "LinearCParty", "family": "rnacentral",
             "series": "rnacentral", "seed": None, "n": 100, "sequence_sha256": "s", "scaffold_sha256": "g",
             "beam": 100, "wall_seconds_median": 10, "rss_gib_median": 2}
        b = {**a, "run": "new", "wall_seconds_median": 5, "rss_gib_median": .2}
        pairs = historical_comparisons([a], [b])
        self.assertEqual(pairs[0]["wall_seconds_ratio"], .5)
        self.assertEqual(pairs[0]["rss_gib_ratio"], .1)
        with self.assertRaisesRegex(ValueError, "S/G/beam mismatch"):
            historical_comparisons([a], [{**b, "scaffold_sha256": "changed"}])
        with self.assertRaisesRegex(ValueError, "cohort mismatch"):
            historical_comparisons([], [b])

    def test_e0_cohorts_and_either_tolerance_counts(self):
        observations = []
        for name, dz, df in (("bpRNA_example", 4e-10, 2.4e-10), ("witness", -2e-7, 1e-7)):
            for kind in ("numerical_parity", "closure_model_difference"):
                observations.append({"input": name, "status": "success", "comparison_kind": kind,
                    "metrics": {"delta_log_z": 0., "delta_free_energy": 0.},
                    **({"default_minus_compatible": {"delta_log_z": dz, "delta_free_energy": df},
                        "closure_attribution_valid": True} if kind == "closure_model_difference" else {})})
        result = {"observations": observations, "run_manifest_sha256": "hash",
            "tolerances": {"log_z_absolute_tolerance": 1e-9, "free_energy_absolute_tolerance": 1e-10},
            "numerical_parity_counts": {"success": 2}, "closure_counts": {"effect_detected": 2}}
        with tempfile.TemporaryDirectory(dir="/tmp") as tmp:
            p = Path(tmp) / "analysis.json"
            p.write_text(json.dumps(result))
            with patch("analyze_e0_e2_runs.RunReader") as reader:
                reader.return_value.manifest_sha256 = "hash"
                reader.return_value.manifest = {"inputs": [{"name": "bpRNA_example"}, {"name": "witness"}]}
                summarize_e0(p, "unused")
            summary = json.loads((p.parent / "summary.json").read_text())
            self.assertEqual(summary["input_counts"], {"bpRNA": 1, "witness": 1})
            rows = list(csv.DictReader((p.parent / "distributions.tsv").open(), delimiter="\t"))
            row = next(r for r in rows if r["cohort"] == "all" and r["comparison"] == "weak_closure_vs_compatible"
                       and r["quantity"] == "delta_log_z" and r["absolute"] == "False")
            self.assertEqual(int(row["count"]), 2)
            self.assertEqual(int(row["exceeds_either_tolerance"]), 2)


if __name__ == "__main__":
    unittest.main()
