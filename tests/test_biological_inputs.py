"""Check input selection and incompatible threshold-pair rejection."""

from pathlib import Path
import hashlib
import json
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from prepare_biological_inputs import prepare, select_records, threshold_scaffold
from run_cli import digest


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

    def test_subset_retry_keeps_scaffold_settings_and_records_extended_budget(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            sequence = "GG" + "A"*3996 + "CC"
            source = root / "panel.fa"
            source.write_text(">first\n" + sequence + "\n>second\n" + sequence + "A\n")
            acquisition = root / "source.json"
            acquisition.write_text(json.dumps({"downloaded_path": str(source), "sha256": digest(source)}))
            lp = root / "linearpartition"
            lp.write_text(f"#!{sys.executable}\n"
                          "import os, pathlib, sys\n"
                          "sequence = sys.stdin.read().strip()\n"
                          "output = pathlib.Path(sys.argv[sys.argv.index('-r')+1])\n"
                          "output.write_text(f'1 {len(sequence)} 0.9\\n2 {len(sequence)-1} 0.8\\n')\n"
                          "print(','.join(os.environ[k] for k in "
                          "('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')))\n")
            lp.chmod(0o755)
            with patch("prepare_biological_inputs.git_state", return_value={"commit": "fixed", "dirty": False}):
                result = prepare(acquisition, root / "retry", lp, "rnacentral", 20261008,
                                 timeout_seconds=3600, record_ids=["first"])
            self.assertEqual(result["status"], "success")
            self.assertEqual(result["cohort_ids_before_subset"], ["first", "second"])
            self.assertEqual(result["selected_ids_before_scaffold"], ["first"])
            self.assertEqual(result["scaffold_settings"],
                             {"beam": 200, "dangles": 2, "pair_probability_threshold": 0.7, "min_helix": 2})
            record = result["records"][0]
            self.assertEqual(record["sequence_sha256"], hashlib.sha256(sequence.encode()).hexdigest())
            self.assertEqual(record["input_sha256"], digest(root / "retry/first.txt"))
            log = root / "retry/scaffold_logs/first"
            command = json.loads((log / "command.json").read_text())
            self.assertEqual(command["timeout_seconds"], 3600)
            self.assertEqual(command["memory_limit_mib"], 8192)
            self.assertEqual(command["argv"][1:8], ["-V", "-b", "200", "--dangles", "2", "-c", "0.7"])
            self.assertEqual((log / "stdout.txt").read_text(), "1,1,1\n")
            self.assertGreater(record["scaffold_elapsed_seconds"], 0)


if __name__ == "__main__":
    unittest.main()
