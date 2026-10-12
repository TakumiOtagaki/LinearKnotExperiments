import random
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from scaffold_features import scaffold_features
from make_synthetic_scaffolds import embed_sequence, generate, make_scaffold, save
from rna_records import scaffold_pairs, read_record
from task_storage import ArchiveReader, MemoryFile


class ScaffoldFeatureTests(unittest.TestCase):
    def check_oracle(self, scaffold, seed=17, d=2, h=8):
        pairs = scaffold_pairs(scaffold)
        rng = random.Random(seed)
        seq = [rng.choice("ACGU") for _ in scaffold]
        for i, j in pairs:
            seq[i - 1], seq[j - 1] = rng.choice(("GC", "CG", "AU", "UA", "GU", "UG"))
        seq = "".join(seq)
        r = scaffold_features(seq, scaffold, hairpin=h, case3_depth=d)
        free = [i for i, s in enumerate(scaffold, 1) if s == "."]
        parents = {g: next(iter(sorted((p for p in pairs if p[0] < g[0] < g[1] < p[1]),
                                      key=lambda p: p[1] - p[0])), None) for g in pairs}
        children = {p: [g for g in pairs if parents[g] == p] for p in [None] + list(pairs)}
        owners = {i: next(iter(sorted((p for p in pairs if p[0] < i < p[1]),
                                     key=lambda p: p[1] - p[0])), None) for i in free}
        paths = {}
        for i in free:
            p, count = owners[i], 0
            while p is not None:
                p = parents[p]
                count += len(children[p]) >= 2
            paths[i] = count
        cross = [(i, j) for i in free for j in free if i < j and
                 any(a < i < b < j or i < a < j < b for a, b in pairs)]
        from scaffold_features import CANONICAL
        self.assertEqual(r["cross_loop_pairs"], len(cross))
        self.assertEqual(r["cross_loop_pairs_turn"], sum(j - i >= 4 for i, j in cross))
        self.assertEqual(r["cross_loop_canonical_pairs_turn"],
                         sum(j - i >= 4 and seq[i - 1] + seq[j - 1] in CANONICAL for i, j in cross))
        self.assertEqual(r["B_G"], sum(1 + sum(a < i < j < b for a, b in pairs) for i, j in pairs))
        self.assertEqual(r["case3_ancestors_max"], max(paths.values(), default=0))
        self.assertEqual(r["case3_path_entries"], sum(min(d, x) for x in paths.values()))
        self.assertEqual(r["case3_truncated_endpoints"], sum(x > d for x in paths.values()))
        hp = [(i, j) for i in free for j in free if 4 <= j - i <= h + 1 and
              set(scaffold[i - 1:j]) == {"."}]
        self.assertEqual(r["hairpin_free_geometries"], len(hp))
        self.assertEqual(r["same_loop_pairs_turn"], sum(i < j - 3 and owners[i] == owners[j]
                                                       for i in free for j in free))
        sibling = 0
        for i in free:
            for j in free:
                if i >= j:
                    continue
                common = [p for p in pairs if p[0] < i < j < p[1]]
                lca = min(common, key=lambda p: p[1] - p[0]) if common else None
                ci = [p for p in children[lca] if p[0] < i < p[1]]
                cj = [p for p in children[lca] if p[0] < j < p[1]]
                sibling += bool(ci and cj and ci != cj)
        self.assertEqual(r["case3_sibling_endpoint_pairs"], sibling)

    def test_independent_small_pair_oracles(self):
        fixtures = ["...........", "(((...)))", "(...)..(...)", "..((...)(...))..",
                    "(.(...)(.(...)(...)).)..", "(.((...).).(...))", "()", "((()))"]
        for scaffold in fixtures:
            for d in (0, 1, 2, 100):
                with self.subTest(scaffold=scaffold, d=d):
                    self.check_oracle(scaffold, d=d)
        rng = random.Random(12)
        def tree(depth):
            if depth == 0 or rng.random() < .4:
                return "." * rng.randrange(1, 5)
            return "(" + tree(depth - 1) + tree(depth - 1) + ")"
        for seed in range(30):
            self.check_oracle(tree(3), seed)

    def test_direct_free_in_branching_loop_is_not_path_entry(self):
        g = "(.(...)(...).)"
        r = scaffold_features(embed_sequence(g, 1), g)
        # Direct loop dots have 0; six hairpin dots have 1.
        self.assertEqual(r["case3_ancestor_hist"], {0: 2, 1: 6})

    def test_errors(self):
        for s, g in [("", ""), ("AAA", "...x"), ("AC", ")("), ("AA", "()"),
                     ("AX", ".."), ("AU", "xx")]:
            with self.subTest(s=s, g=g), self.assertRaises(ValueError):
                scaffold_features(s, g)

    def test_synthetic_matrix_and_feature_controls(self):
        rows = list(generate())
        self.assertEqual(len(rows), 156)
        self.assertEqual(len({r["name"] for r in rows}), len(rows))
        for row in rows:
            f = row["features"]
            self.assertEqual(len(row["sequence"]), len(row["scaffold"]))
            self.assertEqual(f["scaffold_pairs"], 0 if row["family"] == "empty" else 3 * f["n"] // 20)
            self.assertEqual(f["gc_fraction"], .5)
            self.assertEqual(f["B_G"], f["pair_depth_sum"])
        keyed = {(r["features"]["n"], r["seed"], r["family"]): r["features"] for r in rows}
        for n in (250, 500, 1000, 2000, 4000):
            for seed in (2026101201, 2026101202, 2026101203):
                a, b = [keyed[n, seed, k] for k in ("deep_branching", "balanced_branching")]
                for k in ("scaffold_pairs", "free_bases", "free_run_length_hist", "stack_helix_pairs_hist",
                          "branch_points_including_root", "max_pair_children", "cross_loop_pairs",
                          "case3_sibling_endpoint_pairs"):
                    self.assertEqual(a[k], b[k], (n, k))
                a, b = [keyed[n, seed, k] for k in ("separated_helices", "separated_concentrated")]
                for k in ("B_G", "pair_depth_max", "free_runs", "case3_ancestor_hist", "cross_loop_pairs",
                          "stack_helix_pairs_hist", "loop_free_length_hist"):
                    self.assertEqual(a[k], b[k], (n, k))
                self.assertNotEqual(a["free_run_length_hist"], b["free_run_length_hist"])

    def test_cap_boundaries_and_packed_reproducibility(self):
        for depth in (30, 31, 99, 100, 101):
            g = make_scaffold(2000, "branch_depth", depth)
            r = scaffold_features(embed_sequence(g, 5), g, case3_depth=100)
            self.assertEqual(r["case3_ancestors_max"], depth)
            self.assertEqual(r["case3_truncated_endpoints"], 6 if depth == 101 else 0)
        with tempfile.TemporaryDirectory(dir="/tmp") as tmp:
            a, b = Path(tmp) / "a", Path(tmp) / "b"
            first, second = [save(p, (250,), (19,)) for p in (a, b)]
            self.assertEqual(first["files"], second["files"])
            self.assertEqual(len(list(a.iterdir())), 4)
            archive = ArchiveReader(a / "inputs.tar")
            for entry in first["inputs"]:
                name, s, g = read_record(MemoryFile(archive.read(entry["input_member"])))
                self.assertEqual(name, entry["name"])
                scaffold_features(s, g)


if __name__ == "__main__":
    unittest.main()
