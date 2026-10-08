"""Check input selection and incompatible threshold-pair rejection."""

from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from prepare_biological_inputs import select_records, threshold_scaffold


class BiologicalInputTests(unittest.TestCase):
    def test_selection_ignores_file_order_and_rejects_ambiguous_records(self):
        records = [(f"n{n}_{k}", "A" * (n+k)) for n in (50,100,200) for k in range(6)]
        records += [("ambiguous", "A"*60+"N")]
        a, excluded, bins = select_records(records, "bprna", 12)
        b, _, _ = select_records(list(reversed(records)), "bprna", 12)
        self.assertEqual(a,b)
        self.assertEqual(len(a),12)
        self.assertEqual([r["candidate_count"] for r in bins],[6,6,6])
        self.assertIn({"id":"ambiguous","length":61,"reason":"non_ACGU"},excluded)

    def test_scaffold_keeps_stacked_pairs_and_rejects_crossing_candidates(self):
        with tempfile.TemporaryDirectory() as temporary:
            path=Path(temporary)/"bpp.tsv"
            path.write_text("1 9 0.9\n2 8 0.8\n3 7 0.6\n")
            self.assertEqual(threshold_scaffold("GGGAAACCC",path,0.7,2),"((.....))")
            path.write_text("1 8 0.9\n2 9 0.8\n")
            with self.assertRaisesRegex(ValueError,"cross"):
                threshold_scaffold("GGGAAACCC",path,0.7,2)


if __name__ == "__main__":
    unittest.main()
