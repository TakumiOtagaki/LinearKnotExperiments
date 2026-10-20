"""Check concurrency bounds against complete seeded, beam-ordered schedules."""
import copy
import heapq
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from plan_array import make_tasks
from run_cli import plan_jobs
from replan_e1_long_cpu16_89baea52 import calibrated_settings


class ParallelPackingTests(unittest.TestCase):
    def setUp(self):
        self.inputs = [{"name": n, "length": 100, "sha256": n * 64} for n in ("a", "b", "c")]
        self.config = {"schema_version": 1, "purpose": "development", "software_lock": "unused.json",
                       "memory_limit_mib": 4096, "timeout_seconds": 7200,
                       "warmups": 0, "repetitions": 1, "inputs": self.inputs,
                       "modes": [{"name": "b20", "args": ["-b", "20"]},
                                 {"name": "b400", "args": ["-b", "400"]}],
                       "programs": [{"name": p} for p in ("LinearCParty", "LinearPKCapR", "LinearPKRaccess")]}
        observed = {i["name"]: {"sha256": i["sha256"], "memory_mib": 3000,
                                "seconds_by_mode": {"b20": 82, "b400": 82},
                                "max_seconds_by_mode": {"b20": 80, "b400": 80}} for i in self.inputs}
        self.tier = {"name": "long", "min_length": 1, "max_length": 200,
                     "resource": "cpu_16", "h_rt": "00:10:00", "target_seconds": 420,
                     "workers": 8, "memory_budget_mib": 35095, "task_overhead_mib": 1024,
                     "order_kind": "beam_ascending",
                     "cost_model": {"basis": "independent skewed durations", "seconds_per_nt_condition": 1,
                                    "mib_per_nt": 0, "observed_inputs": observed},
                     "parallel_packing": {"kind": "beam_list_bound", "efficiency": .7,
                                          "overhead_seconds": 10, "basis": "explicit scheduling allowance"}}
        self.settings = {"schema_version": 1, "measurement": "accuracy", "seed": 31, "tiers": [self.tier]}

    def test_bound_covers_seeded_schedules_and_keeps_whole_matrices(self):
        task = make_tasks(self.config, self.settings, self.inputs)[0]
        self.assertEqual(len(task["inputs"]), 3)
        self.assertEqual(task["memory_cap_envelope_mib"], 33792)
        self.assertEqual(task["estimated_sequential_seconds"], 492)
        # Heavy jobs and beam barriers prevent the naive 492/8 estimate.
        self.assertGreaterEqual(task["estimated_seconds"], 160)
        for seed in range(30):
            config = {**self.config, "order": {"kind": "beam_ascending", "seed": seed}}
            schedule = plan_jobs(config)
            wall = 0
            for mode in ("b20", "b400"):
                workers = [0] * 8
                for job in (j for j in schedule if j["mode"] == mode):
                    duration = 80 if job["program"] == "LinearCParty" else 1
                    heapq.heappush(workers, heapq.heappop(workers) + duration)
                wall += max(workers)
            self.assertLessEqual(wall / .7 + 10, task["estimated_seconds"])
        self.tier["target_seconds"] = 245
        tasks = make_tasks(self.config, self.settings, self.inputs)
        self.assertEqual([len(t["inputs"]) for t in tasks], [2, 1])

    def test_refuses_changed_inputs_incomplete_costs_and_unsupported_orders(self):
        for change in ("hash", "mode", "memory", "repetitions", "order", "efficiency"):
            config, settings = copy.deepcopy(self.config), copy.deepcopy(self.settings)
            tier = settings["tiers"][0]
            observed = tier["cost_model"]["observed_inputs"]["a"]
            if change == "hash": observed["sha256"] = "different"
            if change == "mode": observed["seconds_by_mode"].pop("b400")
            if change == "memory": observed["memory_mib"] = 4097
            if change == "repetitions": config["repetitions"] = 2
            if change == "order": tier["order_kind"] = "blocked_random"
            if change == "efficiency": tier["parallel_packing"]["efficiency"] = 1.1
            with self.subTest(change=change), self.assertRaises(ValueError):
                make_tasks(config, settings, self.inputs)

    def test_calibration_uses_cli_elapsed_and_limits_transfer_from_two_probes(self):
        config = {**self.config, "inputs": self.inputs[:2]}
        schedule = plan_jobs(config)
        rows = [{**j, "status": "success", "elapsed_seconds": 100, "time_wall_seconds": 1}
                for j in schedule]
        base = {"tasks": [{"planned_jobs": schedule}], "inputs": config["inputs"], "config": config,
                "array_config": self.settings, "software_lock": {"binaries": {"LinearCParty": "binary"}}}
        pilots = [{"input": i["name"], "mode": "b400", "program": "LinearCParty", "status": "success",
                   "memory_limit_mib": 4096, "input_sha256": i["sha256"], "binary_sha256": "binary",
                   "elapsed_seconds": 10, "peak_rss_kib": 1000} for i in config["inputs"]]
        settings, receipt = calibrated_settings(base, {"rows": rows, "source_sha256": {}}, pilots)
        self.assertEqual(receipt["release_factor"], .5)
        self.assertEqual(receipt["historical_cli_hours"], len(rows) * 100 / 3600)
        item = settings["tiers"][0]["cost_model"]["observed_inputs"]["a"]
        self.assertEqual(item["seconds_by_mode"]["b400"], 3 * (100 * .5 + .5))
        _, one_probe = calibrated_settings(base, {"rows": rows, "source_sha256": {}}, pilots[:1])
        self.assertEqual(one_probe["completed_probe_inputs"], 1)
        with self.assertRaisesRegex(ValueError, "exact frozen matrix"):
            calibrated_settings(base, {"rows": rows[:-1], "source_sha256": {}}, pilots)
        rows[0]["status"] = "timeout"
        with self.assertRaisesRegex(ValueError, "historical failure"):
            calibrated_settings(base, {"rows": rows, "source_sha256": {}}, pilots)


if __name__ == "__main__":
    unittest.main()
