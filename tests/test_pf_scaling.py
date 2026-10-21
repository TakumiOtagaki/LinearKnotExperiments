"""Scientific comparison and censoring checks for the new PF summary."""
import copy
import sys
from pathlib import Path
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_pf_scaling import crossings, matched_cparty_modes, summarize, value


class PFScalingTests(unittest.TestCase):
    def test_cparty_reuse_requires_matching_inputs_model_and_parameters(self):
        current = {"inputs": [{"name": "one", "sha256": "bytes", "length": 250, "family": "empty", "seed": 1}],
                   "config": {"model_settings": {"wmb_closure_policy": "LegacyCompatible"},
                              "warmups": 1, "repetitions": 3,
                              "modes": [{"name": "exact"}, {"name": "beam100"}]},
                   "software_lock": {"energy_parameter_sha256": "parameters"}}
        reference = copy.deepcopy(current)
        reference["config"]["modes"].insert(0, {"name": "cparty"})
        self.assertEqual([m["name"] for m in matched_cparty_modes(current, reference)],
                         ["cparty", "exact", "beam100"])
        for change in ("input", "closure", "repetitions", "parameters", "ambiguous"):
            modified = copy.deepcopy(reference)
            if change == "input": modified["inputs"][0]["sha256"] = "different"
            if change == "closure": modified["config"]["model_settings"]["wmb_closure_policy"] = "RequireWeakClosure"
            if change == "repetitions": modified["config"]["repetitions"] = 1
            if change == "parameters": modified["software_lock"]["energy_parameter_sha256"] = "different"
            if change == "ambiguous": modified["config"]["modes"].append({"name": "cparty"})
            with self.subTest(change=change), self.assertRaises(ValueError):
                matched_cparty_modes(current, modified)

    def fixture(self):
        items = [{"name": f"n{n}", "length": n, "family": "empty", "seed": 1} for n in (250, 500)]
        modes = [{"name": m} for m in ("cparty", "exact", "beam100")]
        jobs = []
        for item in items:
            for mode in modes:
                for rep in range(3):
                    seconds = {250: {"cparty": 1, "exact": 2, "beam100": 2},
                               500: {"cparty": 8, "exact": 4, "beam100": 4}}[item["length"]][mode["name"]]
                    jobs.append({"input": item["name"], "mode": mode["name"], "warmup": False,
                                 "repetition": rep, "status": "success", "elapsed_seconds": seconds,
                                 "time_wall_seconds": seconds, "peak_rss_kib": seconds * 1024,
                                 "value": {"log_z": 10, "free_energy": -6, "rt": .6}})
        return items, modes, jobs

    def test_failure_keeps_denominator_and_removes_whole_paired_input(self):
        items, modes, jobs = self.fixture()
        jobs[0].update(status="timeout", elapsed_seconds=120, peak_rss_kib=None)
        warmup = copy.deepcopy(jobs[-1]); warmup.update(warmup=True, elapsed_seconds=10000)
        rows, parity, paired = summarize(jobs + [warmup], items, modes)
        self.assertEqual(len(rows), 6)
        self.assertEqual(rows[0]["successful_repetitions"], 2)
        self.assertEqual({r["input"] for r in paired}, {"n500"})
        self.assertEqual(parity[0]["status"], "unavailable")
        self.assertEqual(next(r for r in rows if r["input"] == "n500" and r["method"] == "beam100")["elapsed_seconds"], 4)

    def test_parity_gate_and_missing_repetition(self):
        items, modes, jobs = self.fixture()
        for r in jobs:
            if r["input"] == "n250" and r["mode"] == "exact":
                r["value"]["log_z"] += .1
        rows, parity, paired = summarize(jobs, items, modes)
        self.assertEqual(parity[0]["status"], "parity_failed")
        self.assertEqual({r["input"] for r in paired}, {"n500"})
        with self.assertRaisesRegex(ValueError, "missing or duplicated"):
            summarize(jobs[:-1], items, modes)

    def test_crossing_requires_measured_bracket(self):
        items, modes, jobs = self.fixture()
        _, _, paired = summarize(jobs, items, modes)
        report = crossings(paired)
        self.assertEqual(report[0]["brackets"][0]["lower_nt"], 250)
        self.assertAlmostEqual(report[0]["brackets"][0]["interpolated_nt"], 2**.5 * 250)
        self.assertFalse(any(r["brackets"] for r in crossings([r for r in paired if r["length"] == 250])))
        with self.assertRaises(ValueError):
            value('{"policy":"legacy-compatible-wmb","log_z":Infinity,"free_energy":-6,"rt":0.6}')

    def test_exact_exclusion_preserves_cparty_beam_pair(self):
        items, modes, jobs = self.fixture()
        modes[1]["length_range"] = {"min_length": 1, "max_length": 250}
        jobs = [j for j in jobs if not (j["input"] == "n500" and j["mode"] == "exact")]
        rows, parity, paired = summarize(jobs, items, modes)
        self.assertEqual(next(r for r in rows if r["input"] == "n500" and r["method"] == "exact")["status"], "excluded_by_plan")
        self.assertEqual({r["method"] for r in paired if r["input"] == "n500"}, {"cparty", "beam100"})
        self.assertEqual(parity[1]["status"], "unavailable")


if __name__ == "__main__":
    unittest.main()
