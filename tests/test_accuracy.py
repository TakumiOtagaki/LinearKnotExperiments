"""Check analysis against an independent two-structure distribution."""

import hashlib
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import analyze_accuracy as accuracy
from prepare_pilot import prepare
from rna_records import read_record, scaffold_pairs
from validate_outputs import validate_outputs


class AccuracyTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.record = ("toy", "ACAAAAGU", "(......)")
        self.rt = 0.6163207755
        # Independent structures: fixed (1,8), with/without generated (2,7).
        self.structures = [(0.6, {(1, 8)}), (0.4, {(1, 8), (2, 7)})]

    def artifacts(self, prefix, structures, free_energy):
        pairs = {}
        contexts = {c: [0.0] * 8 for c in accuracy.CONTEXTS}
        for weight, structure in structures:
            paired = {i for pair in structure for i in pair}
            for pair in structure:
                pairs[pair] = pairs.get(pair, 0) + weight
            for i in range(1, 9):
                contexts["Stem" if i in paired else "Hairpin"][i - 1] += weight
        bpp = self.root / (prefix + ".bpp")
        bpp.write_text(">toy\n" + "".join(f"{i} {j} {p:.17g}\n" for (i, j), p in pairs.items()))
        profile = self.root / (prefix + ".contexts")
        profile.write_text(">toy\n" + "".join(
            c + " " + " ".join(format(p, ".17g") for p in contexts[c]) + "\n"
            for c in accuracy.CONTEXTS))
        access = self.root / (prefix + ".access")
        lines = ["# output_mode=probability", "# access_lens=1,4,6", ">toy",
                 "# scaffold=(......)", f"# folding_rt={self.rt}",
                 f"# ensemble_free_energy={free_energy:.17g}"]
        for start in range(1, 9):
            tokens = []
            for width in (1, 4, 6):
                if start + width - 1 > 8:
                    continue
                interval = set(range(start, start + width))
                p = sum(weight for weight, structure in structures
                        if not interval.intersection(i for pair in structure for i in pair))
                tokens.append(f"{width},{p:.17g};")
            lines.append(str(start - 1) + "\t" + "".join(tokens))
        access.write_text("\n".join(lines) + "\n")
        return (accuracy.read_bpp(bpp, self.record),
                accuracy.read_contexts(profile, self.record),
                accuracy.read_accessibility(access, self.record)), access

    def test_free_base_metrics_and_joint_windows(self):
        reference, _ = self.artifacts("reference", self.structures, 0)
        candidate, _ = self.artifacts("candidate", [(1.0, {(1, 8)})], -self.rt * math.log(0.6))
        pair_metrics = accuracy.compare_bpp(candidate[0], reference[0], self.record)
        self.assertEqual(pair_metrics["free_bases"], 6)
        self.assertAlmostEqual(pair_metrics["generated_pair_l1_sum"], 0.4)
        self.assertAlmostEqual(pair_metrics["base_state_tv_mean_Gfree"], 0.8 / 6)
        context_metrics = accuracy.compare_contexts(candidate[1], reference[1], self.record)
        self.assertAlmostEqual(context_metrics["context_tv_mean_Gfree"], 0.8 / 6)
        metrics = accuracy.compare_accessibility(candidate[2], reference[2], self.record)
        width = metrics["accessibility"]["6"]
        self.assertEqual((width["all_windows"], width["eligible_windows"]), (3, 1))
        self.assertAlmostEqual(width["probability_mae_eligible"], 0.4)
        self.assertAlmostEqual(metrics["delta_log_z"], math.log(0.6))
        self.assertAlmostEqual(reference[2]["windows"][2, 6], 0.6)
        # The two correlated endpoint marginals would instead give 0.36.
        self.assertAlmostEqual(reference[0]["unpaired"][1] * reference[0]["unpaired"][6], 0.36)

    def test_zero_probability_and_incomplete_artifacts_are_distinct(self):
        reference, path = self.artifacts("reference", self.structures, 0)
        candidate = {**reference[2], "windows": dict(reference[2]["windows"])}
        candidate["windows"][2, 6] = 0.0
        metrics = accuracy.compare_accessibility(candidate, reference[2], self.record)
        width = metrics["accessibility"]["6"]
        self.assertEqual(width["reference_positive_candidate_zero"], 1)
        self.assertIsNone(width["opening_mae_both_positive_kcal_mol"])
        self.assertAlmostEqual(width["probability_mae_eligible"], 0.6)
        path.write_text(path.read_text().replace("6,0.59999999999999998;", ""))
        with self.assertRaisesRegex(ValueError, "incomplete"):
            accuracy.read_accessibility(path, self.record)
        broken = self.root / "broken.bpp"
        broken.write_text(">toy\n1 8 1\n2 7 0.6\n2 6 0.6\n")
        with self.assertRaisesRegex(ValueError, "incident mass"):
            accuracy.read_bpp(broken, self.record)

    def test_failed_jobs_and_missing_reference_do_not_become_zero_error(self):
        record_path = self.root / "input.txt"
        record_path.write_text(">toy\nACAAAAGU\n(......)\n")
        job = {"id": "failed", "input": "toy", "mode": "beam20", "program": "LinearCParty",
               "repetition": 0, "warmup": False, "status": "timeout"}
        manifest = {"status": "failed", "purpose": "development", "config_sha256": "fixture",
                    "inputs": [{"name": "toy", "copied_path": str(record_path),
                                "sha256": hashlib.sha256(record_path.read_bytes()).hexdigest()}],
                    "jobs": [job], "software_before": {"commit": None},
                    "experiment_before": {"commit": None},
                    "software_changed_during_run": False, "binaries_changed_during_run": False}
        (self.root / "run.json").write_text(json.dumps(manifest))
        result = accuracy.analyze(self.root, "exact")
        self.assertEqual(result["measured_jobs"], 1)
        row = result["observations"][0]
        self.assertEqual(row["analysis_status"], "job_failed")
        self.assertNotIn("metrics", row)
        self.assertEqual(result["job_status_counts"], {"timeout": 1})
        # A successful job with no exact counterpart remains unpaired.
        self.artifacts("success", [(1.0, {(1, 8)})], 0)
        job.update(status="success", command=["binary", "--bpp", str(self.root / "success.bpp")],
                   directory=".")
        (self.root / "run.json").write_text(json.dumps(manifest))
        row = accuracy.analyze(self.root, "exact")["observations"][0]
        self.assertEqual(row["analysis_status"], "reference_unavailable")
        self.assertNotIn("metrics", row)

    def test_consistency_validation_needs_no_exact_reference(self):
        parsed, _ = self.artifacts("reference", self.structures, 0)
        group = dict(zip(("LinearCParty", "LinearPKCapR", "LinearPKRaccess"), parsed))
        self.assertEqual(accuracy.check_consistency("toy", "beam100", 0, group)["status"], "success")
        group["LinearPKRaccess"]["windows"][2, 1] = 0.1
        self.assertEqual(accuracy.check_consistency("toy", "beam100", 0, group)["status"], "failed")
        record_path = self.root / "input.txt"
        record_path.write_text(">toy\nACAAAAGU\n(......)\n")
        job = {"id": "no_reference", "input": "toy", "mode": "beam100", "program": "LinearCParty",
               "repetition": 0, "warmup": False, "status": "success", "directory": ".",
               "command": ["binary", "--bpp", str(self.root / "reference.bpp")]}
        manifest = {"status": "success", "jobs": [job],
                    "inputs": [{"name": "toy", "copied_path": str(record_path),
                                "sha256": hashlib.sha256(record_path.read_bytes()).hexdigest()}]}
        (self.root / "run.json").write_text(json.dumps(manifest))
        result = validate_outputs(self.root)
        self.assertEqual(result["observations"][0]["validation_status"], "success")
        self.assertEqual(result["complete_groups_checked"], 0)
        self.assertNotIn("metrics", result["observations"][0])
        job["status"] = "timeout"
        (self.root / "run.json").write_text(json.dumps(manifest))
        self.assertEqual(validate_outputs(self.root)["observations"][0]["validation_status"], "job_failed")


class PilotTests(unittest.TestCase):
    def test_reproducible_matched_sequences_and_immutable_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            count, jobs = prepare(directory / "first", directory / "first.json", 20261008)
            prepare(directory / "second", directory / "second.json", 20261008)
            self.assertEqual((count, jobs), (18, 324))
            groups = {}
            for path in sorted((directory / "first").glob("*.txt")):
                self.assertEqual(path.read_bytes(), (directory / "second" / path.name).read_bytes())
                name, sequence, scaffold = read_record(path)
                groups.setdefault(name.rsplit("_", 1)[0], []).append((sequence, scaffold_pairs(scaffold)))
            self.assertEqual(len(groups), 6)
            for records in groups.values():
                self.assertEqual(len({record[0] for record in records}), 1)
                sets = sorted((record[1] for record in records), key=len)
                self.assertEqual([len(pairs) for pairs in sets], [0, 3, 9])
                self.assertTrue(sets[0] <= sets[1] <= sets[2])
            before = (directory / "first/dataset.json").read_bytes()
            with self.assertRaisesRegex(ValueError, "immutable"):
                prepare(directory / "first", directory / "first.json", 20261008)
            self.assertEqual(before, (directory / "first/dataset.json").read_bytes())
