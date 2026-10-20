"""Check unified lock provenance and centroid gain against a small distribution."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import make_software_lock
from compare_release_outputs import compare_bytes


class LockLayoutTests(unittest.TestCase):
    def test_unified_and_legacy_locks_keep_actual_binary_and_configure_args(self):
        for unified in (False, True):
            with self.subTest(unified=unified), tempfile.TemporaryDirectory() as td:
                root = Path(td)
                (root / "params").mkdir()
                (root / "params/rna_Turner04.par").write_text("parameters")
                (root / "CMakeCache.txt").write_text("CMAKE_BUILD_TYPE:STRING=Release\n")
                programs = ("linearknot",) if unified else make_software_lock.PROGRAMS
                for name in programs:
                    (root / name).write_bytes(name.encode())
                out = root / "lock.json"
                argv = ["lock", "--software-root", str(root), "--build-dir", str(root), "--out", str(out),
                        "--configure-argument=-DCMAKE_BUILD_TYPE=Release", "--configure-argument=-DLINEARKNOT_ENABLE_LTO=OFF"]
                state = {"commit": "a"*40, "dirty": False}
                def run(command, **kwargs):
                    text = "linearknot 0.1" if command[-1] == "--version" else "ok"
                    return subprocess.CompletedProcess(command, 0, text, "")
                with patch.object(sys, "argv", argv), patch.object(make_software_lock, "git_state", return_value=state), \
                        patch.object(make_software_lock.subprocess, "run", side_effect=run), contextlib.redirect_stdout(io.StringIO()):
                    make_software_lock.main()
                lock = json.loads(out.read_text())
                self.assertEqual(set(lock["binaries"]), set(programs))
                self.assertEqual(lock["source_repository"].split("/")[-1], "LinearKnot" if unified else "LinearCParty")
                self.assertIn("-DLINEARKNOT_ENABLE_LTO=OFF", lock["configure_arguments"])
                self.assertEqual(len(lock.get("subcommands", [])), 7 if unified else 0)

    def test_release_comparison_retains_metadata_and_numerical_differences(self):
        self.assertTrue(compare_bytes(b"1 2 0.5\n# elapsed_seconds 1\n", b"1 2 0.5\n# elapsed_seconds 2\n")["identical"])
        self.assertFalse(compare_bytes(b"# engine=linearcparty-scfg2\n", b"# engine=linearknot-scfg2\n")["identical"])
        self.assertFalse(compare_bytes(b"1 2 0.5\n", b"1 2 0.6\n")["identical"])


class CentroidGainTests(unittest.TestCase):
    def test_full_gain_matches_expectation_over_truth_structures(self):
        tool = ROOT / "tools/validation_centroid_89baea52_20261020/decoding"
        sys.path.insert(0, str(tool))
        spec = importlib.util.spec_from_file_location("centroid_gain_analysis", tool / "analyze.py")
        module = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(module)
        finally:
            sys.path.remove(str(tool))
        truths = [(set(), .25), ({(1, 2)}, .4), ({(1, 3)}, .35)]
        universe = {(1, 2), (1, 3), (2, 3)}
        d = {"q": [.25, .6, .65], "bpp": {(1, 2): .4, (1, 3): .35}}
        for gamma in (.5, 1, 2):
            for prediction, _ in truths:
                expected = sum(prob * (gamma*len(prediction & truth) + len(universe - (prediction | truth)))
                               for truth, prob in truths)
                self.assertAlmostEqual(module.centroid_score(prediction, d, gamma), expected, places=13)
        # gamma=1 selects the empty prediction; gamma=2 selects the 0.4 pair.
        self.assertGreater(module.centroid_score(set(), d, 1), module.centroid_score({(1, 2)}, d, 1))
        self.assertLess(module.centroid_score(set(), d, 2), module.centroid_score({(1, 2)}, d, 2))


if __name__ == "__main__":
    unittest.main()
