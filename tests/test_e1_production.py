"""Tiny packed arrays exercise the split exact join and widest-beam options."""
import copy
import csv
import json
from pathlib import Path
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import analyze_accuracy as accuracy
import analyze_e1_production as production
import report_e1_production as report
from run_cli import digest
from task_storage import FORMAT, TaskStore


class ProductionTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def pack(self, label, n, modes, bad_old_exact=False, task_number=1):
        run = self.root / label
        run.mkdir()
        record = ("toy", "ACAAAAGU" + "A"*(n-8), "(......)" + "."*(n-8))
        source = self.root / "input.txt"
        source.write_text(">toy\n" + record[1] + "\n" + record[2] + "\n")
        with tarfile.open(run / "inputs.tar", "w") as archive:
            archive.add(source, arcname="inputs/toy.txt")
        entry = dict(name="toy", length=n, sha256=digest(source),
                     input_archive=str(run / "inputs.tar"), input_member="inputs/toy.txt")
        source.unlink()
        config = dict(model_settings={"temperature_celsius": 37}, modes=[], programs=[
            dict(name="LinearCParty", args=["--bpp-cutoff", "0"]),
            dict(name="LinearPKCapR", args=[]),
            dict(name="LinearPKRaccess", args=["--access-len", "1,4,8,16", "--probability"])])
        for mode in modes:
            b = int(mode[4:].split("_")[0]) if mode != "exact" else None
            args = ["--exact"] if b is None else ["-b", str(b), "--max-hairpin", str(max(b,30)),
                "--max-case3-depth", str(max(b,30)), "--ranking", "outer"]
            config["modes"].append(dict(name=mode, args=args))
        tier = "nt50_99" if n < 400 else "nt400_799"
        task = dict(tier=tier, task_id=task_number, directory=f"tasks/{tier}/{task_number:06d}", inputs=[entry])
        path = run / task["directory"] / "run"
        path.mkdir(parents=True)
        store = TaskStore(path)
        schedule = []
        for mode in modes:
            weight = .4 if mode == "exact" or mode.startswith("beam400") else 1.
            energy = -10. if mode == "exact" or mode.startswith("beam400") else -10.1
            q = [0, 1-weight, 1, 1, 1, 1, 1-weight, 0] + [1]*(n-8)
            for p in config["programs"]:
                program = p["name"]
                m = next(m for m in config["modes"] if m["name"] == mode)
                command = ["binary", "-i", "local-input"] + m["args"] + p["args"]
                if program == "LinearCParty": command += ["--bpp", "local-bpp"]
                job = dict(id=f"toy-{mode}-{program}", input="toy", mode=mode, program=program,
                           command=command, warmup=False, repetition=0)
                schedule.append(job)
                row = dict(job, directory="jobs/"+job["id"], status="failed" if bad_old_exact and mode == "exact" else "success",
                           elapsed_seconds=2., time_wall_seconds=1., peak_rss_kib=100)
                if program == "LinearCParty":
                    stdout = ">toy\n"
                elif program == "LinearPKCapR":
                    contexts = {c: [0.]*n for c in accuracy.CONTEXTS}
                    contexts.update(Stem=[1-p for p in q], Exterior=q)
                    stdout = ">toy\n" + "".join(c + " " + " ".join(map(str,v)) + "\n" for c,v in contexts.items())
                else:
                    stdout = (">toy\n# output_mode=probability\n# access_lens=1,4,8,16\n# scaffold=" + record[2]
                              + f"\n# folding_rt=.6\n# ensemble_free_energy={energy}\n")
                    for s in range(1,n+1):
                        values = []
                        for w in production.WIDTHS:
                            if s+w-1 > n: continue
                            interval = set(range(s,s+w))
                            value = 0 if interval & {1,8} else 1-weight if interval & {2,7} else 1
                            values.append(f"{w},{value};")
                        stdout += str(s-1) + "\t" + "".join(values) + "\n"
                if bad_old_exact and mode == "exact": stdout = "do not parse old exact"
                contents = {"stdout.txt": stdout, "command.json": json.dumps({"argv": command}),
                            "result.json": json.dumps(row), "time.tsv": "1\t100\n"}
                if program == "LinearCParty": contents["bpp.tsv"] = f">toy\n1 8 1\n2 7 {weight}\n"
                files = {}
                for filename, content in contents.items():
                    raw = self.root / ("raw-" + filename)
                    raw.write_text(content)
                    files[row["directory"] + "/" + filename] = raw
                store.commit(row, files)
        base = dict(output_format=FORMAT, config=config, software_before={"commit": "fixture"},
                    software_lock={"energy_parameter_sha256": "params", "binaries": {"binary": "hash"}},
                    software_changed_during_run=False, binaries_changed_during_run=False,
                    software_lock_valid_after_run=True, inputs=[entry])
        (path / "run.json").write_text(json.dumps(dict(base, status="failed" if bad_old_exact else "success",
            planned_jobs=schedule, jobs=[], host={"hostname": "fixture"}, workers=2)))
        (run / "array_plan.json").write_text(json.dumps(dict(base, tasks=[task], inputs_archive_sha256=digest(run / "inputs.tar"))))
        return run

    def test_split_exact_ignores_old_exact_and_keeps_signed_and_zero_counts(self):
        beam = self.pack("beam", 50, ["exact"]+[f"beam{b}_outer" for b in (20,50,100,200)], True)
        exact = self.pack("exact", 50, ["exact"], task_number=7)
        out = self.root / "out"
        with patch("tarfile.TarFile.extract", side_effect=AssertionError("extraction")), \
             patch("tarfile.TarFile.extractall", side_effect=AssertionError("extraction")):
            p = production.write_metrics(beam, exact, "unused", out)
            summary = production.summarize(out, make_plots=False)
        self.assertEqual(p["excluded_job_status_counts"], {"failed": 3})
        self.assertEqual(p["selected_job_status_counts"], {"success": 15})
        self.assertEqual(p["consistency_checks"], {"success": 5})
        self.assertEqual(summary["comparison_kind"], "exact_error")
        group = summary["groups"][0]
        self.assertAlmostEqual(group["metrics"][production.CORE[0]]["median"], -.1)
        self.assertAlmostEqual(group["metrics"]["base_state_tv_mean_Gfree"]["median"], 1.2/48)
        self.assertEqual(group["window_counts"]["access_w1_reference_positive_candidate_zero"]["total_windows"], 2)
        empty = next(g for g in summary["groups"] if g["bin"] == "nt100_199")
        self.assertIsNone(empty["metrics"][production.CORE[0]]["median"])
        self.assertEqual(len(list(out.iterdir())), 6)

    def test_long_uses_400_and_rejects_200_reference(self):
        run = self.pack("long", 400, [f"beam{b}_outer" for b in (20,50,100,200,300,400)])
        with self.assertRaisesRegex(ValueError, "maximum measured"):
            production.write_metrics(run, None, "beam200_outer", self.root / "bad")
        out = self.root / "out"
        production.write_metrics(run, None, "beam400_outer", out)
        summary = production.summarize(out, make_plots=False)
        self.assertEqual(summary["comparison_kind"], "beam_width_sensitivity")
        b300, b400 = summary["groups"][4:6]
        self.assertAlmostEqual(b300["metrics"][production.CORE[0]]["median"], -.1)
        self.assertEqual(b400["metrics"]["base_state_tv_mean_Gfree"]["max"], 0)
        with (out / "metrics.tsv").open() as stream:
            first = next(csv.DictReader(stream, delimiter="\t"))
        report.validate_row(first)
        broken = dict(first, access_w1_reference_positive_candidate_zero="999")
        with self.assertRaises(AssertionError):
            report.validate_row(broken)

    def test_combined_exact_uses_the_same_frozen_plan_and_rejects_another_matrix(self):
        run = self.pack("combined", 50, ["exact"]+[f"beam{b}_outer" for b in (20,50,100,200)])
        plan_path = run / "array_plan.json"
        plan = json.loads(plan_path.read_text())
        plan["output_dir"] = str(run)
        plan_path.write_text(json.dumps(plan))
        result = production.write_metrics(run, run, "unused", self.root / "combined-out")
        self.assertEqual(result["selected_job_status_counts"], {"success": 15})
        other = copy.deepcopy(plan)
        other["output_dir"] = str(self.root / "other")
        with self.assertRaisesRegex(ValueError, "identical combined plan"):
            production.check_exact_modes(plan, other)
        other["output_dir"] = plan["output_dir"]
        other["config"]["modes"].pop()
        with self.assertRaisesRegex(ValueError, "identical combined plan"):
            production.check_exact_modes(plan, other)

    def test_options_and_cross_revision_guard(self):
        a = production.arguments(["short", "--beam-run", "b", "--exact-run", "e", "--out", "o"])
        self.assertEqual(a.exact_run, Path("e"))
        a = production.arguments(["long", "--run", "r", "--reference", "beam400_outer", "--out", "o"])
        self.assertEqual(a.reference, "beam400_outer")
        beam = {"software_before": {"commit": "unverified"}, "config": {}, "software_lock": {}}
        exact = copy.deepcopy(beam); exact["software_before"]["commit"] = "different"
        with self.assertRaisesRegex(ValueError, "unverified"):
            production.check_reference_plans(beam, exact)

    def test_summary_only_regenerates_tables_and_preserves_empty_domains(self):
        run = self.pack("long", 400, [f"beam{b}_outer" for b in (20,50,100,200,300,400)])
        out = self.root / "out"
        production.write_metrics(run, None, "beam400_outer", out)
        # Model a valid RNA with no eligible w16 window, as in existing accuracy tests.
        with (out / "metrics.tsv").open() as stream:
            saved = list(csv.DictReader(stream, delimiter="\t"))
        for row in saved:
            for key in production.NUMERIC:
                if key.startswith("access_w16_") and ("mae" in key or "max_error" in key):
                    row[key] = "NA"
                elif key.startswith("access_w16_") and key != "access_w16_all_windows":
                    row[key] = "0"
        with (out / "metrics.tsv").open("w") as stream:
            writer = csv.DictWriter(stream, production.FIELDS, delimiter="\t")
            writer.writeheader(); writer.writerows(saved)
        p = json.loads((out / "provenance.json").read_text())
        p["metrics_sha256"] = digest(out / "metrics.tsv")
        (out / "provenance.json").write_text(json.dumps(p))
        with patch.object(production, "plot_summary"):
            production.main(["long", "--summary-only", "--out", str(out)])
        summary = json.loads((out / "summary.json").read_text())
        metric = summary["groups"][0]["metrics"]["access_w16_probability_mae_eligible"]
        self.assertEqual(metric["n"], 0)
        self.assertEqual(metric["missing_or_out_of_domain_rna"], 1)
        self.assertIsNone(metric["median"])
        self.assertTrue((out / "quantiles.tsv").is_file())
        self.assertTrue((out / "window_counts.tsv").is_file())


if __name__ == "__main__":
    unittest.main()
