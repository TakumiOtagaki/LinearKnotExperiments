"""Validate raw-output interpretation without running RNA software."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from analyze_small_rerun_89baea52 import mfe_square


class RawRecoveryTests(unittest.TestCase):
    def test_scaffold_echo_is_not_the_mfe_structure(self):
        text = ">toy\nGAAAACAA\n(....)..\n[....].. (-1.23)\nensemble free energy: -1.24 kcal/mol\n"
        self.assertEqual(mfe_square(text, 8), {(1, 6)})

    def test_absent_prediction_is_not_an_empty_pair_set(self):
        with self.assertRaisesRegex(ValueError, "MFE structure missing"):
            mfe_square(">toy\nGAAAACAA\n(....)..\n", 8)


if __name__ == "__main__":
    unittest.main()
