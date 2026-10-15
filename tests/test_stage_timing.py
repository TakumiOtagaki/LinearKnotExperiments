"""Small public-API/CLI parity checks; never a performance benchmark."""

import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from build_stage_timing_driver import build
from plan_array import prepare
from run_array_task import run_task
from run_cli import write_json

SOURCE = ROOT / "deps/LinearCParty-d21618c5"
LIB_BUILD = ROOT / "build/linear-cparty-d21618c5"
DRIVER_BUILD = ROOT / "build/stage-timing-d21618c5"
DRIVER_LOCK = "manifests/software-stage-timing-d21618c5-20261016.json"
STAGES = ["inside", "inside_repeat", "mfe", "outside_bpp", "loop_profile", "accessibility", "sampling", "mea"]
CONTEXTS = ["Stem", "Hairpin", "Bulge", "Internal", "Multiloop", "Exterior", "PseudoLoop"]
PK_S = "AAAGGGAAAGGGCCCAAACCCAAA"
PK_G = "...(((......)))........."
# Same 24 nt record as data/smoke.
CASES = [(PK_S, PK_G, 3, 12, 1),
         ("GCAUAUUAAUAUUCUGCUAUAAAUAUAGGUGCACUUAGGGCAUGC",
          "((.(((....))(.(.....))).((....((.(....)))))))", 20, 45, 45),
         (PK_S * 6, PK_G * 6, 5, 30, 2)]


def strict_json(raw):
    return json.loads(raw, parse_constant=lambda token: (_ for _ in ()).throw(ValueError(token)))


class StageTimingTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        needed = [DRIVER_BUILD / "stage_timing_driver", LIB_BUILD / "LinearCParty",
                  LIB_BUILD / "LinearPKCapR", LIB_BUILD / "LinearPKRaccess",
                  LIB_BUILD / "LinearPKSampling", LIB_BUILD / "LinearPKMEA"]
        if not all(p.is_file() for p in needed):
            raise unittest.SkipTest("build_stage_timing_driver.py and the pinned CLI build are required")
        cls.env = {**os.environ, "CPARTY_FIXED_ENERGY_PARAM_FILE": str(SOURCE / "params/rna_Turner04.par"),
                   "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "LC_ALL": "C"}

    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="stage-timing-test-", dir="/tmp")
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)

    def launch(self, binary, args, expected=0):
        p = subprocess.run([str(binary), *args], text=True, capture_output=True,
                           env=self.env, timeout=30)
        self.assertEqual(p.returncode, expected, p.stderr + p.stdout[:1000])
        return p

    def driver(self, args, expected=0):
        p = self.launch(DRIVER_BUILD / "stage_timing_driver", args, expected)
        self.assertEqual(len(p.stdout.splitlines()), 1)
        return strict_json(p.stdout)

    def test_all_stages_match_clis_and_shared_session(self):
        for sequence, scaffold, beam, hairpin, depth in CASES:
            with self.subTest(n=len(sequence), b=beam, h=hairpin, d=depth):
                args = ["--sequence", sequence, "--scaffold", scaffold, "-b", str(beam),
                        "--max-hairpin", str(hairpin), "--max-case3-depth", str(depth)]
                record = self.driver([*args, "--session-mode", "both", "--include-values",
                                      "--access-len", "1,2,4,7", "-k", "8", "--seed", "20261013", "--gamma", "2"])
                self.assertTrue(record["ok"])
                self.assertEqual(record["config"]["beam"], beam)
                self.assertEqual(record["config"]["max_hairpin"], hairpin)
                self.assertEqual(record["config"]["max_case3_depth"], depth)
                fresh, shared = record["runs"]
                self.assertEqual([s["stage"] for s in fresh["stages"]], STAGES)
                for left, right in zip(fresh["stages"], shared["stages"]):
                    self.assertEqual(left["summary"], right["summary"])
                    self.assertEqual(left["ok"], right["ok"])
                    for stage in (left, right):
                        self.assertGreaterEqual(stage["elapsed_seconds"], 0)
                        self.assertFalse(stage["graph_reused"])
                        self.assertEqual(stage["error"], "")
                        for point in ("rss_before", "rss_after"):
                            self.assertGreater(stage[point]["vmrss_kib"], 0)
                            self.assertGreaterEqual(stage[point]["vmhwm_kib"], stage[point]["vmrss_kib"])
                        self.assertGreaterEqual(stage["rss_after"]["vmhwm_kib"], stage["rss_before"]["vmhwm_kib"])
                self.assertFalse(fresh["stages"][1]["inside_result_cached"])
                self.assertTrue(shared["stages"][1]["inside_result_cached"])
                self.assertIsNone(shared["stages"][1]["api_elapsed_seconds"])
                self.assertEqual(fresh["stages"][0]["summary"], fresh["stages"][1]["summary"])
                s = {stage["stage"]: stage["summary"] for stage in fresh["stages"]}
                for name in ("outside_bpp", "loop_profile", "accessibility", "sampling", "mea"):
                    self.assertAlmostEqual(s[name]["log_z"], s["inside"]["log_z"], places=10)

                bpp_path = self.directory / "bpp.tsv"
                folding = self.launch(LIB_BUILD / "LinearCParty", [*args, "--bpp", str(bpp_path)]).stdout
                energy = float(re.search(r"\((-?\d+(?:\.\d+)?)\)", folding)[1])
                free_energy = float(re.search(r"ensemble free energy:\s*(-?[0-9.]+)", folding)[1])
                self.assertLessEqual(abs(energy - s["mfe"]["energy"]), 0.00501)
                self.assertLessEqual(abs(free_energy - s["inside"]["free_energy"]), 0.00005001)
                self.assertIn(s["mfe"]["union_structure"], folding)
                cli_pairs = {(int(i), int(j)): float(p) for i, j, p in
                             (line.split() for line in bpp_path.read_text().splitlines() if line and not line.startswith(">"))}
                driver_pairs = {(p[0], p[1]): p[4] for p in s["outside_bpp"]["pairs"]}
                self.assertEqual(cli_pairs.keys(), driver_pairs.keys())
                for pair in cli_pairs:
                    self.assertAlmostEqual(cli_pairs[pair], driver_pairs[pair], places=12)
                self.assertAlmostEqual(sum(cli_pairs.values()), s["outside_bpp"]["bpp_sum"], places=12)
                self.assertAlmostEqual(2 * s["outside_bpp"]["bpp_sum"] + s["outside_bpp"]["unpaired_sum"], len(sequence), places=10)

                capr = self.launch(LIB_BUILD / "LinearPKCapR", args).stdout
                cli_context = {parts[0]: list(map(float, parts[1:])) for parts in
                               (line.split() for line in capr.splitlines() if line and not line.startswith(">"))}
                for c, name in enumerate(CONTEXTS):
                    values = [row[c] for row in s["loop_profile"]["contexts"]]
                    self.assertEqual(len(values), len(sequence))
                    for a, b in zip(cli_context[name], values): self.assertAlmostEqual(a, b, places=12)
                    self.assertAlmostEqual(sum(cli_context[name]), s["loop_profile"]["context_sums"][c], places=12)
                self.assertEqual(s["loop_profile"]["contexts"], s["accessibility"]["contexts"])
                raccess = self.launch(LIB_BUILD / "LinearPKRaccess", [*args, "--access-len", "1,2,4,7", "--probability"]).stdout
                cli_access = {w: [] for w in (1, 2, 4, 7)}
                for line in raccess.splitlines():
                    if not line or line[0] in ">#": continue
                    for term in line.split("\t")[1].split(";"):
                        if term:
                            w, p = term.split(","); cli_access[int(w)].append(float(p))
                for window in s["accessibility"]["accessibility"]:
                    expected = cli_access[window["length"]]
                    self.assertEqual(len(expected), window["window_count"])
                    for a, b in zip(expected, window["probabilities"]): self.assertAlmostEqual(a, b, places=12)
                    self.assertAlmostEqual(sum(expected), window["probability_sum"], places=12)
                for a, b in zip(s["accessibility"]["accessibility"][0]["probabilities"], s["outside_bpp"]["unpaired"]):
                    self.assertAlmostEqual(a, b, places=11)
                metadata = dict(re.findall(r"^# (\w+)=(.+)$", raccess, re.M))
                self.assertAlmostEqual(float(metadata["ensemble_free_energy"]), s["accessibility"]["free_energy"], places=12)

                sampled = self.launch(LIB_BUILD / "LinearPKSampling", [*args, "-k", "8", "--seed", "20261013"]).stdout
                self.assertEqual([line for line in sampled.splitlines() if line and not line.startswith(">")][1:],
                                 s["sampling"]["union_structures"])
                self.assertEqual(s["sampling"]["sample_count"], 8)
                mea = self.launch(LIB_BUILD / "LinearPKMEA", [*args, "--gamma", "2", "--stats"])
                self.assertIn(s["mea"]["union_structure"], mea.stdout)
                self.assertAlmostEqual(float(re.search(r"EA=([^;]+)", mea.stdout)[1]), s["mea"]["expected_accuracy"], places=12)
                self.assertAlmostEqual(float(re.search(r"log_z=(\S+)", mea.stderr)[1]), s["inside"]["log_z"], places=10)
                mea_stage = fresh["stages"][-1]
                self.assertGreaterEqual(mea_stage["api_decoding_seconds"], 0)
                self.assertLessEqual(mea_stage["api_decoding_seconds"], mea_stage["elapsed_seconds"])

    def test_one_record_normalization_json_escaping_and_isolated_stage(self):
        path = self.directory / "record.txt"
        path.write_text('# comment\n> a"b\\c\td\n\ngggtttccc\n.........\n')
        record = self.driver(["-i", str(path), "--stage", "outside_bpp"])
        self.assertEqual(record["input"]["name"], 'a"b\\c\td')
        self.assertEqual(record["input"]["sequence"], "GGGUUUCCC")
        self.assertEqual(len(record["runs"][0]["stages"]), 1)
        self.assertIsNone(record["runs"][0]["stages"][0]["summary"]["pairs"])
        path.write_text(path.read_text() + "AAAA\n....\n")
        self.assertFalse(self.driver(["-i", str(path)], 2)["ok"])

    def test_isolated_stages_ranking_and_complete_work_match_cli(self):
        # Separate processes, with the same explicit b/h/d as the E2 matrix.
        for name in ("empty_scaffold.txt", "pk_scaffold.txt"):
            for ranking in (["--ranking", "outer"], ["--ranking", "baseline"],
                            ["--ranking", "outer", "--outer-pilot-exponent", "0.5"]):
                with self.subTest(input=name, ranking=ranking):
                    args = ["-i", str(ROOT / "data/smoke" / name), "-b", "100",
                            "--max-hairpin", "100", "--max-case3-depth", "100", *ranking]
                    summaries = {}
                    for stage in ("inside", "mfe", "outside_bpp"):
                        record = self.driver([*args, "--stage", stage, "--skip-inside-repeat", "--include-values"])
                        stages = record["runs"][0]["stages"]
                        self.assertEqual([s["stage"] for s in stages], [stage])
                        self.assertEqual(record["config"]["ranking"], ranking[1])
                        self.assertFalse(record["config"]["inside_repeat"])
                        summaries[stage] = stages[0]["summary"]
                        work = summaries[stage]["work"]
                        fields = {"item_visits", "positive_items", "frontier_candidates", "frontier_selected",
                                  "vp_candidates", "vp_selected", "wmbp_candidates", "wmbp_selected",
                                  "case3_candidates", "case3_selected", "binary_join_probes", "extension_probes",
                                  "v_parent_probes", "vp_loop_parent_probes", "case3_truncated_endpoints", "families"}
                        self.assertEqual(set(work), fields)
                        families = work["families"]
                        self.assertTrue(families)
                        self.assertEqual(len({f["nonterminal"] for f in families}), len(families))
                        self.assertIn("V", {f["nonterminal"] for f in families})
                        for key in ("item_visits", "positive_items", "frontier_candidates", "frontier_selected"):
                            self.assertEqual(sum(f[key] for f in families), work[key])
                        self.assertGreater(work["positive_items"], 0)
                    self.assertEqual(summaries["inside"]["work"], summaries["outside_bpp"]["work"])
                    self.assertAlmostEqual(summaries["inside"]["log_z"], summaries["outside_bpp"]["log_z"], places=12)
                    bpp = self.directory / "isolated-bpp.tsv"
                    folding = self.launch(LIB_BUILD / "LinearCParty", [*args, "--bpp", str(bpp), "--bpp-cutoff", "0"]).stdout
                    mfe = float(re.search(r"\((-?\d+(?:\.\d+)?)\)", folding)[1])
                    pf = float(re.search(r"ensemble free energy:\s*(-?[0-9.]+)", folding)[1])
                    self.assertLessEqual(abs(mfe - summaries["mfe"]["energy"]), 0.00501)
                    self.assertLessEqual(abs(pf - summaries["inside"]["free_energy"]), 0.00005001)
                    self.assertIn(summaries["mfe"]["union_structure"], folding)
                    expected = {(int(i), int(j)): float(p) for i, j, p in
                                (line.split() for line in bpp.read_text().splitlines() if line and not line.startswith(">"))}
                    observed = {(p[0], p[1]): p[4] for p in summaries["outside_bpp"]["pairs"]}
                    self.assertEqual(expected.keys(), observed.keys())
                    for pair in expected:
                        self.assertAlmostEqual(expected[pair], observed[pair], places=12)

    def test_invalid_options_and_api_errors_are_strict_json(self):
        args = ["--sequence", "GGGAAACCC", "--scaffold", "........."]
        for extra in (["-b", "0"], ["-k", "-1"], ["--seed", "1x"], ["--gamma", "nan"],
                      ["--gamma", "-1"], ["--access-len", "0"], ["--access-len", "1,"],
                      ["--access-len", "1,1"], ["--session-mode", "graph"], ["--stage", "beta"], ["--seed"],
                      ["--ranking", "invalid"], ["--outer-pilot-exponent", "0.7"],
                      ["--ranking", "baseline", "--outer-pilot-exponent", "0.5"]):
            with self.subTest(extra=extra): self.assertFalse(self.driver([*args, *extra], 2)["ok"])
        record = self.driver(["--sequence", "ACGU", "--scaffold", "(", "--stage", "inside", "--session-mode", "both"], 1)
        self.assertFalse(record["ok"])
        for run in record["runs"]:
            for stage in run["stages"]:
                self.assertFalse(stage["ok"]); self.assertTrue(stage["error"])
                self.assertIsNone(stage["summary"]["log_z"])

    def test_zero_samples_and_windows_longer_than_sequence(self):
        args = ["--sequence", "GGGAAACCC", "--scaffold", ".........", "--include-values"]
        sample = self.driver([*args, "--stage", "sampling", "-k", "0"])
        self.assertEqual(sample["runs"][0]["stages"][0]["summary"]["sample_count"], 0)
        access = self.driver([*args, "--stage", "accessibility", "--access-len", "1,50"])
        window = access["runs"][0]["stages"][0]["summary"]["accessibility"][1]
        self.assertEqual(window, {"length": 50, "window_count": 0, "probability_sum": 0, "probabilities": []})

    def test_unreachable_root_keeps_all_stage_timings_and_api_errors(self):
        # Balanced/canonical scaffold, but its hairpin is too short for TURN.
        record = self.driver(["--sequence", "GC", "--scaffold", "()", "-b", "1",
                              "--session-mode", "both", "--include-values"], 1)
        self.assertFalse(record["ok"])
        for run in record["runs"]:
            self.assertEqual([s["stage"] for s in run["stages"]], STAGES)
            stages = {s["stage"]: s for s in run["stages"]}
            self.assertTrue(stages["inside"]["ok"])
            self.assertFalse(stages["inside"]["reachable"])
            self.assertIsNone(stages["inside"]["summary"]["log_z"])
            for name in ("loop_profile", "accessibility", "sampling"):
                self.assertFalse(stages[name]["ok"])
                self.assertTrue(stages[name]["error"])
                self.assertGreaterEqual(stages[name]["elapsed_seconds"], 0)
            for window in stages["accessibility"]["summary"]["accessibility"]:
                self.assertIsNone(window["probability_sum"])
                self.assertIsNone(window["probabilities"])

    def test_build_rejects_unverified_dependency_before_creating_files(self):
        lock = json.loads((ROOT / "manifests/software-pilot-d21618c5.json").read_text())
        lock["build_revision_verified"] = False
        path = self.directory / "unverified.json"; write_json(path, lock)
        target = self.directory / "build"
        with self.assertRaisesRegex(ValueError, "certify"):
            build(SOURCE, LIB_BUILD, target, path)
        self.assertFalse(target.exists())

    def test_build_rejects_mismatched_e0_archive_before_creating_files(self):
        archive = json.loads((ROOT / "manifests/software-e0-d21618c5-20261014.json").read_text())
        for entry in archive["extra_files"]:
            if entry["path"].endswith("libRNA.a"):
                entry["sha256"] = "0" * 64
        path = self.directory / "bad-archive.json"; write_json(path, archive)
        target = self.directory / "build"
        with self.assertRaisesRegex(ValueError, "differs from E0 lock"):
            build(SOURCE, LIB_BUILD, target, ROOT / "manifests/software-pilot-d21618c5.json", path)
        self.assertFalse(target.exists())

    def test_existing_array_runner_packs_json_and_resumes_without_jobs(self):
        config = json.loads((ROOT / "configs/stage_timing_smoke_15aa7742.json").read_text())
        # The historical lock deliberately freezes the old observer source.
        config["software_lock"] = DRIVER_LOCK
        config_path = self.directory / "config.json"; write_json(config_path, config)
        plan = prepare(config_path, ROOT / "configs/arrays/stage_timing_smoke.json",
                       self.directory / "plan", self.directory / "scripts", SOURCE, DRIVER_BUILD)
        self.assertEqual(len(plan["tasks"]), 1)
        task = plan["tasks"][0]
        scratch = self.directory / "scratch"; scratch.mkdir()
        with patch.dict(os.environ, {"SGE_TASK_ID": "1"}):
            self.assertEqual(run_task(self.directory / "plan/array_plan.json", "smoke", software=SOURCE,
                                      build=DRIVER_BUILD, scratch_dir=scratch), 0)
        output = self.directory / "plan" / task["directory"] / "run"
        self.assertFalse((output / "jobs").exists())
        results = output / "task_results.jsonl"
        rows = [strict_json(line) for line in results.read_text().splitlines()]
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(row["status"] == "success" for row in rows))
        archive = output / "task_outputs.tar"
        with tarfile.open(archive) as tar:
            members = [m for m in tar.getmembers() if m.name.endswith("stdout.txt")]
            self.assertEqual(len(members), 4)
            for member in members:
                self.assertTrue(strict_json(tar.extractfile(member).read())["ok"])
        before_results = results.read_bytes()
        with tarfile.open(archive) as tar:
            before_members = {m.name: tar.extractfile(m).read() for m in tar.getmembers() if m.isfile()}
        with patch.dict(os.environ, {"SGE_TASK_ID": "1"}):
            self.assertEqual(run_task(self.directory / "plan/array_plan.json", "smoke", software=SOURCE,
                                      build=DRIVER_BUILD, scratch_dir=scratch, resume=True), 0)
        self.assertEqual(results.read_bytes(), before_results)
        # Checkpoint recovery can trim tar end padding, preserving all raw data.
        with tarfile.open(archive) as tar:
            self.assertEqual({m.name: tar.extractfile(m).read() for m in tar.getmembers() if m.isfile()}, before_members)
        self.assertEqual(list(scratch.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
