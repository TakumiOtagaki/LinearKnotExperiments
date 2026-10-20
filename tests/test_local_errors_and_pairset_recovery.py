"""Scientific definitions used by the 2026-10-19 reanalysis."""
import itertools
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))

from analyze_e1_local_errors_89baea52 import tail_counts
from analyze_e1_pk_stratified import crossing_pairs
from analyze_pk_membership import CANONICAL, parse_structure
from analyze_pk_recovery_pairset_89baea52 import crossing_set, mfe_pairs


class LocalErrorDefinitions(unittest.TestCase):
    def test_strict_low_and_inclusive_high_boundaries_and_sparse_zero(self):
        exact = {(1,6): .5, (2,7): .9, (3,8): .099, (4,9): .1}
        beam = {(1,6): .1, (2,7): .099, (3,8): .5, (4,9): .8, (5,10): .7}
        counts, missed, created = tail_counts(exact, beam)
        self.assertEqual(counts['missed_events'], 1)
        self.assertEqual(counts['created_events'], 2)
        self.assertEqual({r[0] for r in missed}, {(2,7)})
        self.assertEqual({r[0] for r in created}, {(3,8), (5,10)})
        self.assertEqual(counts['exact_ge_half'], 2)
        self.assertEqual(counts['beam_ge_half'], 3)
        self.assertAlmostEqual(counts['max_absolute_error'], .801)

    def test_empty_domain_has_no_events(self):
        counts, missed, created = tail_counts({}, {})
        self.assertFalse(missed or created)
        self.assertTrue(all(v == 0 for v in counts.values()))

    def test_crossing_masks_against_endpoint_geometry(self):
        for text in ('((..))..', '(..)(..)', '((.)(.))', '.(....).'):
            record = ('test', 'A'*len(text), text)
            g = {(i,j) for i,j,_ in parse_structure(text)}
            pairs = {p: .7 for p in itertools.combinations(range(1,len(text)+1), 2)
                     if all(text[i-1] == '.' for i in p)}
            self.assertEqual(set(crossing_pairs(record, {'pairs': pairs})), crossing_set(set(pairs), g))


class RecoveryDefinitions(unittest.TestCase):
    def test_bracket_color_swap_preserves_coordinate_metric(self):
        expected = {(1,5), (2,8)}
        for text in ('([..)..]', '[(..]..)'):
            triples = mfe_pairs(text+' (-1.0)',8)
            self.assertEqual({(i,j) for i,j,_ in triples}, expected)
            self.assertEqual(crossing_set(expected,{(1,5)}), {(2,8)})

    def test_curly_reference_pair_is_selected_by_crossing(self):
        ref = '({..)..}'
        pairs = {(i,j) for i,j,_ in parse_structure(ref)}
        self.assertEqual(crossing_set(pairs,{(1,5)}), {(2,8)})

    def test_removal_of_noncanonical_scaffold_can_remove_crossing(self):
        # The square pair crosses only a UU pair, which is absent from G.
        seq, ref = 'UCAAUAAG', '([..)..]'
        triples = parse_structure(ref)
        canonical = {(i,j) for i,j,_ in triples if seq[i-1]+seq[j-1] in CANONICAL}
        self.assertEqual(canonical, {(2,8)})
        # Use coordinates independently of the symbol assigned to the pair.
        self.assertEqual(crossing_set(canonical, {(1,5)}), {(2,8)})
        self.assertEqual(crossing_set(canonical, set()), set())


if __name__ == '__main__':
    unittest.main()
