"""Release tooling checks using tiny fake CLIs and node-local temporary files."""

import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from compare_release_outputs import PROGRAMS, compare, compare_bytes, specification
from make_rerun_configs import generate, write_configs
from prepare_release import inspect_build, prepare, test_build_link as build_link


class ReleaseComparisonTests(unittest.TestCase):
    def test_numeric_differences_missing_lines_and_only_elapsed_is_removed(self):
        result = compare_bytes(b"# elapsed_seconds 1\nx 1.0 -2e-4\n", b"# elapsed_seconds 3\nx 1.5 -3e-4\n")
        self.assertFalse(result["identical"])
        self.assertEqual(result["max_absolute_numeric_difference"], .5)
        self.assertTrue(result["all_differences_numeric"])
        self.assertEqual(result["different_line_examples"][0]["line"], 1)
        self.assertFalse(compare_bytes(b"x 1\n", b"y 1\n")["all_differences_numeric"])
        self.assertFalse(compare_bytes(b"x 1\n", b"")["all_differences_numeric"])
        self.assertFalse(compare_bytes(b"x 1\n", b"x 1")["identical"])
        self.assertTrue(compare_bytes(b"# elapsed_seconds 1\nx 1\n", b"# elapsed_seconds 4\nx 1\n")["identical"])

    def test_three_clis_split_exact_reference_failure_and_scratch_cleanup(self):
        with tempfile.TemporaryDirectory(dir="/tmp") as directory:
            root = Path(directory)
            builds = [root / n for n in ("beam", "exact")]
            for build in builds:
                build.mkdir()
                for name in PROGRAMS:
                    binary = build / name
                    binary.write_text(f"#!{sys.executable}\nimport pathlib,sys,time\n"
                                      "print('# elapsed_seconds',time.monotonic())\nprint('value 1.25')\n"
                                      "if '--bpp' in sys.argv: pathlib.Path(sys.argv[sys.argv.index('--bpp')+1]).write_text('1 2 0.125\\n')\n")
                    binary.chmod(0o755)
            source, parameter = root / "input.txt", root / "energy.par"
            source.write_text(">tiny\nAAAA\n....\n")
            parameter.write_text("test parameter")
            spec = specification(builds[0], builds[0], builds[1], [source], parameter, [20, 200], builds[1])
            self.assertEqual(compare(spec, dry_run=True)["planned_processes"], 18)
            report = compare(spec, scratch_dir=root)
            self.assertTrue(report["identical"])
            self.assertEqual(report["identical_comparisons"], 9)
            for row in report["rows"]:
                self.assertGreater(row["candidate"]["peak_rss_kib"], 0)
                if row["mode"] == "exact":
                    self.assertEqual(Path(row["candidate"]["argv"][0]).parent, builds[1])
            self.assertFalse(list(root.glob("lcp-task-*")))
            bad = builds[0] / "LinearPKCapR"
            bad.write_text(f"#!{sys.executable}\nimport sys\nsys.exit(7)\n")
            spec = specification(builds[0], builds[0], builds[1], [source], parameter, [20], builds[1])
            report = compare(spec, scratch_dir=root)
            self.assertFalse(report["successful"])
            self.assertFalse(report["identical"])


class ReleasePreparationTests(unittest.TestCase):
    def test_new_release_stage_order_suite_environment_and_idempotence(self):
        import argparse
        from types import SimpleNamespace
        with tempfile.TemporaryDirectory(dir="/tmp") as directory:
            root = Path(directory)
            commit = "a" * 40
            source = root / "deps/LinearCParty-aaaaaaaa"
            build = root / "build/linear-cparty-aaaaaaaa"
            lock = root / "manifests/software-pilot-aaaaaaaa.json"
            lock.parent.mkdir()
            input_path = root / "input.txt"
            input_path.write_text(">tiny\nAAAA\n....\n")
            args = argparse.Namespace(clone_from=root, commit=commit, software_root=None, build_dir=None, lock=None,
                out=root / "receipt.json", dry_run=False, estimate_only=False, jobs=8, estimated_build_files=100,
                fallback_remaining_files=100000, scratch_dir=root, beam_reference_build=build, exact_reference_build=build,
                candidate_exact_build=None, input=[input_path], beams="20,200", timeout=2)
            calls = []

            def run(argv, cwd=None, env=None, stdout=None, **kwargs):
                calls.append(argv)
                if argv[:2] == ["git", "clone"]:
                    (source / "params").mkdir(parents=True)
                    (source / "params/rna_Turner04.par").write_text("parameter")
                    (source / "test/tools").mkdir(parents=True)
                    (source / "test/tools/run_required_scfg2_tests.sh").write_text("fake required suite")
                elif argv[:2] == ["cmake", "-S"]:
                    build.mkdir(parents=True)
                    (build / "CMakeCache.txt").write_text(f"CMAKE_BUILD_TYPE:STRING=Release\nCMAKE_HOME_DIRECTORY:INTERNAL={source}\n")
                elif argv[:2] == ["cmake", "--build"]:
                    for name in PROGRAMS:
                        (build / name).write_text("fake binary")
                elif argv[0] == "bash":
                    self.assertEqual(env["SCFG2_BUILD_DIR"], str(build))
                    self.assertEqual((source / "build").resolve(), build)
                    self.assertTrue(Path(env["TMPDIR"]).is_relative_to(root))
                    stdout.write(b"[scfg2-required] required set ok\n")
                elif argv[0] == sys.executable:
                    self.assertFalse((source / "build").exists())
                    lock.write_text('{}')
                return SimpleNamespace(returncode=0)

            with patch("prepare_release.ROOT", root), \
                    patch("prepare_release.git", side_effect=lambda _, command, *a: commit if command == "rev-parse" else "file.cc"), \
                    patch("prepare_release.inspect_source", return_value={"commit": commit, "dirty": False}), \
                    patch("prepare_release.verify_lock", return_value={"build_revision_verified": True, "required_suite": "required set ok"}), \
                    patch("file_budget.query_quota", return_value=None), \
                    patch("prepare_release.subprocess.run", side_effect=run), \
                    patch("prepare_release.compare", side_effect=lambda spec, *a: {"specification": spec, "successful": True, "identical": True}) as comparison:
                report = prepare(args)
                self.assertEqual([c[0] for c in calls], ["git", "git", "cmake", "cmake", "bash", sys.executable])
                self.assertEqual(report["stages"]["required_suite"]["status"], "passed")
                self.assertTrue(report["stages"]["comparison"]["identical"])
                calls.clear()
                prepare(args)
                self.assertEqual(calls, [])
                self.assertEqual(comparison.call_count, 1)

    def test_temporary_build_link_cleanup_and_conflicting_link_preserved(self):
        with tempfile.TemporaryDirectory(dir="/tmp") as directory:
            source, build = Path(directory) / "source", Path(directory) / "compiled"
            source.mkdir(); build.mkdir()
            with self.assertRaisesRegex(RuntimeError, "test failure"):
                with build_link(source, build):
                    self.assertEqual((source / "build").resolve(), build)
                    raise RuntimeError("test failure")
            self.assertFalse((source / "build").is_symlink())
            (source / "build").symlink_to(build)
            with build_link(source, build):
                pass
            self.assertTrue((source / "build").is_symlink())
            with self.assertRaisesRegex(ValueError, "elsewhere"):
                with build_link(source, Path(directory) / "other"):
                    pass
            self.assertEqual((source / "build").resolve(), build)

    def test_existing_release_is_read_only_and_stale_build_is_rejected(self):
        import argparse
        with tempfile.TemporaryDirectory(dir="/tmp") as directory:
            root = Path(directory)
            source, build = root / "source", root / "compiled"
            (source / "params").mkdir(parents=True); build.mkdir()
            (source / "params/rna_Turner04.par").write_text("test")
            cache = build / "CMakeCache.txt"
            cache.write_text(f"CMAKE_BUILD_TYPE:STRING=Release\nCMAKE_HOME_DIRECTORY:INTERNAL={source}\n")
            for name in PROGRAMS:
                (build / name).write_text("fake binary")
            lock, input_path = root / "lock.json", root / "input.txt"
            lock.write_text('{}')
            input_path.write_text(">tiny\nAAAA\n....\n")
            args = argparse.Namespace(clone_from=source, commit="1"*40, software_root=source, build_dir=build, lock=lock,
                out=root / "receipt.json", dry_run=True, estimate_only=True, jobs=8, estimated_build_files=20000,
                fallback_remaining_files=100000, scratch_dir=root, beam_reference_build=build, exact_reference_build=build,
                candidate_exact_build=None, input=[input_path], beams="20,200", timeout=2)
            before = {p: p.read_bytes() for p in [cache, lock, *(build / n for n in PROGRAMS)]}
            with patch("prepare_release.git", side_effect=lambda _, command, *a: "1"*40 if command == "rev-parse" else "file.cc"), \
                    patch("prepare_release.inspect_source", return_value={"commit": "1"*40, "dirty": False}), \
                    patch("prepare_release.verify_lock", return_value={"build_revision_verified": True, "required_suite": "required set ok"}), \
                    patch("file_budget.query_quota", return_value=None), \
                    patch("prepare_release.subprocess.run", side_effect=AssertionError("dry-run executed a command")):
                report = prepare(args)
                self.assertEqual(report["stages"]["required_suite"]["status"], "verified_existing_lock")
                self.assertEqual(report["stages"]["comparison"]["planned_processes"], 18)
                self.assertEqual(prepare(args)["release"], report["release"])
            self.assertEqual(before, {p: p.read_bytes() for p in before})
            cache.write_text("CMAKE_BUILD_TYPE:STRING=Release\nCMAKE_HOME_DIRECTORY:INTERNAL=/wrong\n")
            with self.assertRaisesRegex(ValueError, "another source"):
                inspect_build(source, build)


class RerunConfigTests(unittest.TestCase):
    def setUp(self):
        self.settings = dict(revision="rehearsal", date="20261016", software_root="deps/LinearCParty-d21618c5",
            build_dir="build/linear-cparty-d21618c5", software_lock="manifests/software-pilot-d21618c5.json",
            stage_build_dir="build/stage-timing-d21618c5", stage_software_lock="manifests/software-stage-timing-d21618c5-20261016.json",
            workers=8, memory_limit_mib=4096, task_overhead_mib=1024, parallel_slowdown_factor=1.5)

    def test_templates_keep_scientific_settings_and_separate_stage_lock(self):
        variants = generate(self.settings)
        self.assertEqual(len(variants), 7)
        for name, c, a in variants[:5]:
            original = json.loads((ROOT / c["rerun"]["template_config"]).read_text())
            self.assertEqual(c["rerun"]["historical_template_analysis"], original.get("analysis"))
            for key in ("modes", "programs", "inputs", "input_cohort", "warmups", "repetitions", "order"):
                self.assertEqual(c.get(key), original.get(key), (name, key))
            for t in a["tiers"]:
                self.assertLessEqual(t["workers"] * c["memory_limit_mib"] + t["task_overhead_mib"], t["memory_budget_mib"])
                old_model = c["rerun"]["original_cost_models"][t["name"]]
                self.assertEqual(t["cost_model"]["seconds_per_nt_condition"], old_model["seconds_per_nt_condition"])
                self.assertEqual(t["cost_model"]["safety_factor"], old_model.get("safety_factor", 1) * 1.5)
            if "memory_breakdown" in name:
                self.assertEqual(c["software_lock"], self.settings["stage_software_lock"])
        one, eight = variants[-2:]
        self.assertEqual(one[1]["inputs"], eight[1]["inputs"])
        self.assertEqual(one[1]["rerun"]["selected_lengths"], [4793, 5482, 6990])
        self.assertEqual([v[2]["tiers"][0]["workers"] for v in (one, eight)], [1, 8])
        self.assertTrue(eight[2]["parallel_timing"]["allowed"])

    def test_memory_oversubscription_and_invalid_settings_are_refused_before_writing(self):
        for changes in ({"workers": 10}, {"memory_limit_mib": 8192}, {"parallel_slowdown_factor": .5},
                        {"h_rt": "25:00:00"}, {"workers": 17}, {"resource": "cpu_8"}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                generate({**self.settings, **changes})
        cpu40 = generate({**self.settings, "resource": "cpu_40", "workers": 10, "memory_limit_mib": 8192})
        for _, config, array in cpu40:
            for tier in array["tiers"]:
                self.assertEqual(tier["resource"], "cpu_40")
                self.assertLessEqual(tier["workers"] * config["memory_limit_mib"] + tier["task_overhead_mib"],
                                     tier["memory_budget_mib"])
        with tempfile.TemporaryDirectory(dir="/tmp") as directory:
            out = Path(directory) / "new"
            with self.assertRaises(ValueError):
                write_configs({**self.settings, "workers": 10}, out)
            self.assertFalse(out.exists())
            report = write_configs(self.settings, out)
            self.assertEqual(len(report["variants"]), 7)
            with self.assertRaisesRegex(ValueError, "new config"):
                write_configs(self.settings, out)


if __name__ == "__main__":
    unittest.main()
