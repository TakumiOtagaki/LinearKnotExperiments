"""Fixed E0 provenance, historical WIP regression and resource limits."""

import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from analyze_e0 import analysis_settings
import build_e0_tools as builder
from check_e0_fixed_witnesses import wip_comparison
from plan_array import make_tasks
from rna_records import read_record
from run_cli import plan_jobs


class FixedOracleTests(unittest.TestCase):
    def test_old_missing_reordered_and_wrong_base_oracles_are_rejected(self):
        manifest = {"upstream_commit": builder.ORACLE_COMMIT,
                    "patches": [{"name": n + "-fixture.patch"} for n in builder.PATCH_NUMBERS]}
        with tempfile.TemporaryDirectory() as temporary:
            owner = Path(temporary)
            config = owner / "test/oracle/legacy_cparty/oracle_config.py"
            config.parent.mkdir(parents=True)

            def verify(value):
                config.write_text("def verify_checkout(path):\n    return " + repr(value) + "\n")
                with patch.object(builder, "git_state", return_value={"commit": builder.ORACLE_OWNER_COMMIT, "dirty": False}):
                    return builder.verify_oracle(owner, owner / "oracle")

            self.assertEqual(verify(manifest), manifest)
            old = copy.deepcopy(manifest); old["patches"].pop()
            reordered = copy.deepcopy(manifest); reordered["patches"].reverse()
            wrong_base = copy.deepcopy(manifest); wrong_base["upstream_commit"] = "0" * 40
            for value in (old, reordered, wrong_base):
                with self.assertRaisesRegex(ValueError, "five-patch"):
                    verify(value)

    def test_wrong_or_dirty_source_is_rejected_before_build_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); output = root / "output"
            for state in ({"commit": "0" * 40, "dirty": False},
                          {"commit": builder.LINEAR_COMMIT, "dirty": True}):
                with patch.object(builder, "git_state", return_value=state):
                    with self.assertRaisesRegex(ValueError, "clean fixed"):
                        builder.build(root, root, root, output)
                    with self.assertRaisesRegex(ValueError, "clean fixed"):
                        builder.verify_oracle(root, root)
                self.assertFalse(output.exists())


class WipRegressionTests(unittest.TestCase):
    def values(self, energy):
        return {program: {"log_z": -energy / .6163207755, "free_energy": energy, "rt": .6163207755}
                for program in ("E0LegacyCompatible", "E0CPartyOracle")}

    def test_both_new_values_change_and_each_version_keeps_parity(self):
        result = wip_comparison(self.values(13.082603367763157), self.values(13.1), analysis_settings())
        self.assertTrue(result["passed"])
        self.assertGreater(result["new_minus_old"]["E0CPartyOracle"]["delta_log_z"], .028)

    def test_unchanged_one_side_wrong_direction_or_parity_failure_is_rejected(self):
        old = self.values(13.1)
        one_side = self.values(13.082603367763157)
        one_side["E0CPartyOracle"] = old["E0CPartyOracle"]
        bad_old = self.values(13.1); bad_old["E0CPartyOracle"] = self.values(13.2)["E0CPartyOracle"]
        for new, previous in ((old, old), (one_side, old), (self.values(13.2), old),
                              (self.values(13.082603367763157), bad_old)):
            self.assertFalse(wip_comparison(new, previous, analysis_settings())["passed"])


class FixedE0PlanTests(unittest.TestCase):
    def setUp(self):
        self.config = json.loads((ROOT / "configs/e0_bprna_legacy_d21618c5_20261014.json").read_text())
        self.settings = json.loads((ROOT / "configs/arrays/e0_bprna_legacy_d21618c5_20261014.json").read_text())
        self.previous = json.loads((ROOT / "configs/e0_bprna_legacy_20261011.json").read_text())
        self.cohort = json.loads((ROOT / self.config["input_cohort"]).read_text())
        self.inputs = []
        for entry in self.config["inputs"]:
            _, sequence, _ = read_record(ROOT / entry["path"])
            self.inputs.append({**entry, "length": len(sequence),
                                "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest()})
        self.inputs += [{"name": r["id"], "length": r["length"], "sequence_sha256": r["sequence_sha256"]}
                        for r in self.cohort["records"]]

    def test_same_304_inputs_plus_completed_three_base_hairpin_witness(self):
        self.assertEqual(self.config["input_cohort"], self.previous["input_cohort"])
        self.assertEqual(self.config["inputs"][:-1], self.previous["inputs"])
        self.assertEqual((len(self.inputs), len({i["name"] for i in self.inputs})), (305, 305))
        self.assertEqual(self.config["inputs"][-1]["name"], "WIP29")
        _, sequence, scaffold = read_record(ROOT / self.config["inputs"][-1]["path"])
        self.assertEqual((len(sequence), sequence[2:7], scaffold[2:7]), (29, "GAAAC", "(...)"))
        self.assertEqual(self.config["analysis"], analysis_settings())

    def test_complete_matrix_packs_with_cpu16_cpu40_and_inode_gate(self):
        tasks = make_tasks(self.config, self.settings, self.inputs)
        self.assertEqual(len(tasks), 21)
        self.assertEqual(len(plan_jobs({**self.config, "inputs": self.inputs})), 915)
        for task in tasks:
            length = max(i["length"] for i in task["inputs"])
            self.assertEqual(task["requested_resource"], "cpu_16=1" if length < 200 else "cpu_40=1")
            self.assertEqual(task["workers"], 1)
            tier = next(t for t in self.settings["tiers"] if t["name"] == task["tier"])
            self.assertLessEqual(task["memory_cap_envelope_mib"], tier["memory_budget_mib"])
            self.assertLessEqual(task["estimated_memory_mib"], tier["memory_budget_mib"])
        self.assertTrue(self.settings["compact_plan"])
        self.assertLessEqual(self.settings["file_budget"]["max_remaining_fraction"], .1)

    def test_former_128gib_cap_is_rejected_on_cpu40(self):
        settings = copy.deepcopy(self.settings)
        settings["tiers"][-1]["job_limits"]["memory_limit_mib"] = 131072
        with self.assertRaisesRegex(ValueError, "cap plus overhead"):
            make_tasks(self.config, settings, self.inputs)


if __name__ == "__main__":
    unittest.main()
