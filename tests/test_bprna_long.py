"""Independent packed long-E1 fixture and checks of reference semantics."""
import csv
import json
import math
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import analyze_accuracy as accuracy
import analyze_bprna_long as long
from run_cli import digest
from scaffold_features import scaffold_features
from task_storage import FORMAT, RunReader, TaskStore


class LongAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.run = self.root / "array"
        self.run.mkdir()
        record = ("toy", "ACAAAAGU" + "A" * 392, "(......)" + "." * 392)
        self.record = record
        source = self.root / "toy.txt"
        source.write_text(">toy\n" + record[1] + "\n" + record[2] + "\n")
        with tarfile.open(self.run / "inputs.tar", "w") as archive:
            archive.add(source, arcname="inputs/toy.txt")
        entry = dict(name="toy", length=400, sha256=digest(source),
                     input_archive=str(self.run / "inputs.tar"), input_member="inputs/toy.txt")
        source.unlink()  # Only the archive is available to the analysis.
        config = dict(model_settings={"temperature_celsius": 37},
                      modes=[dict(name=m, args=["-b", str(b), "--max-hairpin", str(max(b,30)),
                             "--max-case3-depth", str(max(b,30))]) for m,b in zip(long.MODES,long.BEAMS)],
                      programs=[dict(name="LinearCParty", args=["--bpp-cutoff", "0"]),
                                dict(name="LinearPKCapR", args=[]),
                                dict(name="LinearPKRaccess", args=["--access-len", "1,4,8,16", "--probability"])])
        self.base = dict(output_format=FORMAT, config=config, software_before={"commit": "fixture"},
                         software_lock={"energy_parameter_sha256": "params", "binaries": {"binary": "hash"}},
                         inputs=[entry], software_changed_during_run=False, binaries_changed_during_run=False,
                         software_lock_valid_after_run=True)
        self.task = self.run / "tasks/nt400_799/000001/run"
        self.task.mkdir(parents=True)
        store = TaskStore(self.task)
        schedule = []
        for mode, beam, weight, free_energy in zip(long.MODES, long.BEAMS, (0,.2,.1,.3), (-9,-9.7,-9.6,-9.9)):
            q = [0, 1-weight, 1, 1, 1, 1, 1-weight, 0] + [1] * 392
            context = {c: [0.] * 400 for c in accuracy.CONTEXTS}
            context["Stem"] = [1-p for p in q]
            context["Exterior"] = q
            for program in sorted(accuracy.PROGRAMS):
                mp = next(m for m in config["modes"] if m["name"] == mode)
                pp = next(p for p in config["programs"] if p["name"] == program)
                command = ["binary", "-i", "local-input"] + mp["args"] + pp["args"]
                if program == "LinearCParty":
                    command += ["--bpp", "local-output"]
                job = dict(id=f"toy-{mode}-{program}", input="toy", mode=mode, program=program,
                           warmup=False, repetition=0, command=command)
                schedule.append(job)
                row = {**job, "status": "success", "returncode": 0, "error": None,
                       "directory": "jobs/" + job["id"], "elapsed_seconds": 2.,
                       "time_wall_seconds": 1., "peak_rss_kib": 100}
                if program == "LinearPKCapR":
                    stdout = ">toy\n" + "".join(c + " " + " ".join(map(str,context[c])) + "\n" for c in accuracy.CONTEXTS)
                elif program == "LinearPKRaccess":
                    stdout = (">toy\n# output_mode=probability\n# access_lens=1,4,8,16\n# scaffold=" + record[2]
                              + f"\n# folding_rt=.6\n# ensemble_free_energy={free_energy}\n")
                    for s in range(1,401):
                        tokens = []
                        for w in (1,4,8,16):
                            if s+w-1 > 400:
                                continue
                            interval = set(range(s,s+w))
                            p = 0 if interval & {1,8} else 1-weight if interval & {2,7} else 1
                            tokens.append(f"{w},{p};")
                        stdout += str(s-1) + "\t" + "".join(tokens) + "\n"
                else:
                    stdout = ">toy\n"
                contents = {"stdout.txt": stdout, "command.json": json.dumps({"argv": command}),
                            "result.json": json.dumps(row), "time.tsv": "1\t100\n"}
                if program == "LinearCParty":
                    contents["bpp.tsv"] = ">toy\n1 8 1\n" + (f"2 7 {weight}\n" if weight else "")
                files = {}
                for name, content in contents.items():
                    path = self.root / ("raw-" + name)
                    path.write_text(content)
                    files[row["directory"] + "/" + name] = path
                store.commit(row, files)
        (self.task / "run.json").write_text(json.dumps({**self.base, "status": "success", "planned_jobs": schedule,
            "jobs": [], "host": {"hostname": "fixture"}, "workers": 1, "resume_history": []}))
        tasks = [dict(tier="nt400_799", task_id=1, directory="tasks/nt400_799/000001", inputs=[entry]),
                 dict(tier="nt800_1599", task_id=2, directory="does-not-exist", inputs=[])]
        (self.run / "array_plan.json").write_text(json.dumps({**self.base, "tasks": tasks,
            "inputs_archive_sha256": digest(self.run / "inputs.tar")}))
        self.feature_path = self.root / "features.tsv"
        feature = {"cohort": "bprna", "name": "toy", **scaffold_features(record[1],record[2])}
        with self.feature_path.open("w") as stream:
            writer = csv.DictWriter(stream, list(feature), delimiter="\t")
            writer.writeheader()
            writer.writerow(feature)

    def test_packed_reference_and_nonmonotonicity_without_extraction(self):
        output = self.root / "analysis"
        with patch("tarfile.TarFile.extract", side_effect=AssertionError("extraction")), \
             patch("tarfile.TarFile.extractall", side_effect=AssertionError("extraction")):
            provenance = long.write_metrics([(self.run,("nt400_799",))], self.feature_path, output)
            summary = long.summarize(output, make_plots=False)
        self.assertEqual(provenance["jobs"], {"success":12})
        self.assertEqual(provenance["consistency_checks"], {"success":4})
        self.assertEqual(provenance["runs"][0]["excluded_planned_tasks"], {"nt800_1599":1})
        geometry = summary["bins"][0]["local_cap_geometry"]
        self.assertEqual(geometry["maximum_case3_branching_ancestors"], 0)
        self.assertEqual(geometry["sequences_with_free_run_gt_h_plus2"]["200"], 1)
        with (output / "metrics.tsv").open() as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        self.assertNotIn("exact", " ".join(rows[0]))
        self.assertAlmostEqual(float(rows[2]["delta_F_vs_b200_kcal_mol"]), .3)
        self.assertAlmostEqual(float(rows[2]["Z_ratio_vs_b200"]), math.exp(-.3/.6))
        self.assertAlmostEqual(float(rows[2]["base_state_tv_mean_Gfree"]), 2*.2/398)
        self.assertAlmostEqual(float(rows[2]["context_tv_mean_Gfree"]), 2*.2/398)
        self.assertEqual(float(rows[3]["delta_F_vs_b200_kcal_mol"]), 0)
        self.assertEqual(float(rows[3]["Z_ratio_vs_b200"]), 1)
        for metric in ("free_energy_kcal_mol", "log_Z_ratio_vs_b200", "base_state_tv_mean_Gfree"):
            m = next(m for m in summary["monotonicity"] if m["bin"] == "all" and m["metric"] == metric)
            self.assertEqual(m["monotone_sequences"], 0)
        self.assertEqual(len(list(output.iterdir())), 4)

    def test_negative_difference_ratio_above_one_and_no_eligible_windows(self):
        record = ("paired", "A"*16, "()"*8)
        pair = {"pairs": {(i,i+1):1. for i in range(1,17,2)}, "unpaired": [0.]*16}
        context = {c: [float(c == "Stem")]*16 for c in accuracy.CONTEXTS}
        access = {"windows": {(s,w):0. for w in (1,4,8,16) for s in range(1,18-w)},
                  "widths": [1,4,8,16], "rt": .6, "free_energy": -10}
        reference = {"LinearCParty":pair, "LinearPKCapR":context, "LinearPKRaccess":access}
        candidate = {**reference, "LinearPKRaccess":{**access, "free_energy":-11}}
        row = long.compare_condition(candidate,reference,record)
        self.assertEqual(row["delta_F_vs_b200_kcal_mol"], -1)
        self.assertGreater(row["Z_ratio_vs_b200"], 1)
        self.assertIsNone(row["base_state_tv_mean_Gfree"])
        self.assertEqual(row["access_w4_eligible_windows"], 0)
        self.assertIsNone(row["access_w4_probability_mae_eligible"])

    def test_reject_modified_raw_data_and_feature_join(self):
        reader = RunReader(self.task)
        job = reader.manifest["jobs"][0]
        archive = reader.archive(self.task / "task_outputs.tar")
        offset, size = archive.members[job["directory"] + "/time.tsv"]
        with (self.task / "task_outputs.tar").open("r+b") as stream:
            stream.seek(offset)
            stream.write(b"9")
        with self.assertRaisesRegex(ValueError,"hash changed"):
            long.audit_job(reader,job)
        features = long.load_features(self.feature_path)
        features["toy"]["scaffold_sha256"] = "bad"
        with self.assertRaisesRegex(ValueError,"hash or length mismatch"):
            long.join_features(self.record,features)

    def test_thresholds_signed_absolute_strict_and_empty_bins(self):
        counts = long.threshold_counts([-.11,-.1,0,.1,.11,float("nan")])["0.1"]
        self.assertEqual(counts["eligible_sequences"], 5)
        self.assertEqual(counts["positive"]["count"], 1)
        self.assertEqual(counts["negative"]["count"], 1)
        self.assertEqual(counts["absolute"]["count"], 2)
        self.assertEqual(long.statistics([])["n"],0)
        self.assertIsNone(long.statistics([])["p99"])
        self.assertIsNone(long.threshold_counts([])["1"]["absolute"]["fraction"])

    def test_partial_rank_correlation_handles_length_and_constants(self):
        n = np.arange(10.)
        result = long.association(n, n*n, n)
        self.assertAlmostEqual(result["spearman"], 1)
        self.assertIsNone(result["partial_spearman_length"])
        self.assertIsNone(long.association(np.ones(10),n,n)["spearman"])

    def test_pilot_selects_distinct_challenges_and_near_length_control(self):
        names = [f"rna{i}" for i in range(6)]
        metadata = [{k: str(i) for k in ("sequence_sha256", "scaffold_sha256", "run", "task", "input_sha256")}
                    for i in range(6)]
        data = {k: np.ones((6,4)) for k in long.FEATURES}
        data["length"] = np.tile(np.array([400,401,402,403,500,600])[:,None], (1,4))
        data[long.CORE[0]] = np.tile(np.array([10,9,0,.1,.3,8])[:,None], (1,4))
        data["base_state_tv_mean_Gfree"] = np.tile(np.array([.5,.6,0,.01,.02,.2])[:,None], (1,4))
        data["context_tv_mean_Gfree"] = data["base_state_tv_mean_Gfree"]
        data["free_energy_kcal_mol"] = np.tile(np.array([0,-1,-2,-3]), (6,1))
        data["free_energy_kcal_mol"][4] = [0,1,-2,-3]
        for k in ("Z_ratio_vs_b200", "LinearCParty_elapsed_seconds", "LinearCParty_peak_rss_kib"):
            data[k] = np.ones((6,4))
        rows = long.candidate_rows(names,metadata,data)
        self.assertEqual([r["input"] for r in rows], ["rna0","rna1","rna4","rna2"])
        self.assertEqual(rows[-1]["reason"], "low_BPP_sensitivity_length_matched_control")
        # If the only nonmonotone RNA is already selected for a large F gap,
        # the third candidate must be labelled for its F gap, not monotonicity.
        data["free_energy_kcal_mol"][4] = [0,-1,-2,-3]
        data["free_energy_kcal_mol"][0] = [0,1,-2,-3]
        rows = long.candidate_rows(names,metadata,data)
        self.assertEqual(rows[2]["input"], "rna5")
        self.assertEqual(rows[2]["reason"], "second_large_abs_F_remaining")

    def test_reject_wrong_cap_and_missing_condition(self):
        reader = RunReader(self.task)
        job = {**reader.manifest["jobs"][0]}
        job["command"] = job["command"][:]
        job["command"][job["command"].index("--max-hairpin")+1] = "31"
        with self.assertRaisesRegex(ValueError,"planned condition"):
            long.check_command(job,self.base["config"])
        rows = (self.task / "task_results.jsonl").read_text().splitlines()
        (self.task / "task_results.jsonl").write_text("\n".join(rows[:-1])+"\n")
        with self.assertRaisesRegex(ValueError,"incomplete or failed"):
            long.write_metrics([(self.run,("nt400_799",))],self.feature_path,self.root/"partial")


if __name__ == "__main__":
    unittest.main()
