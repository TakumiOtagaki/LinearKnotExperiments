# Synthetic inputs used in Supplementary Figure S4

The saved set has 156 inputs in
[`data/e2-synthetic-shapes-20261012/manifest.json`](../data/e2-synthetic-shapes-20261012/manifest.json)
and `inputs.jsonl`; its generator is
[`scripts/make_synthetic_scaffolds.py`](../scripts/make_synthetic_scaffolds.py).
The adopted run configuration is
[`configs/linearknot_4f97fd3_20261021/e2-synthetic-shapes-4f97fd3-20261021.json`](../configs/linearknot_4f97fd3_20261021/e2-synthetic-shapes-4f97fd3-20261021.json).
It records all inputs and their hashes, one warm-up, three measurements,
model settings, and the requested outputs. This note documents existing
inputs; it does not change or regenerate them.

Lengths are 250, 500, 1,000, 2,000 and 4,000 nt. At each length,
seeds 2026101201, 2026101202 and 2026101203 each generate seven main shapes:
empty, long_stem, separated_helices, separated_concentrated, separated_12,
deep_branching (nested multiloops), and balanced_branching.
These 105 inputs are plotted in Figure S4. Nonempty shapes have
`m = 3*n//20 = floor(0.15*n)` scaffold pairs. Three- and twelve-helix
shapes have three-base hairpins, pair counts divided as evenly as possible,
and at least one free base in each exterior gap. Remaining free bases
are divided evenly (earlier gaps receive the remainder), or placed in the
first gap for separated_concentrated. The long stem has a three-base
hairpin and evenly divided exterior flanks.

For the two multiloop shapes, a full binary tree has `k` branching nodes,
`k+1` hairpins and `2*k+1` pairs. The main choice is
`k = min(101, (m-1)//2)`. Nested trees put all remaining branching in
the right child; balanced trees divide it between the children.
The remaining `m-(2*k+1)` pairs form a stacked outer trunk.
Free positions outside that tree are divided between its two exterior flanks.
Thus the matched trees preserve pair count, branch count, leaves, trunk
length and the multiset of free-run lengths.

The additional 51 inputs use the nested tree with each of
`k = 15, 30, 31, 99, 100, 101` that satisfies `k <= min(101,(m-1)//2)`,
excluding a duplicate of the main nested tree. Branch counts by length are:
250: 15; 500 and 1,000: 15, 30, 31; 2,000 and 4,000: 15, 30, 31, 99, 100.
Together with three seeds this is `(1+3+3+5+5)*3 = 51` inputs.
The main nested shape already covers 101 levels where feasible.
These inputs probe the search-depth limits 30 and 100 and are excluded
from Figure S4's seven main series.

For each shape and length, the actual sequence seed is `seed+n`.
`embed_sequence()` assigns canonical AU, UA, GC or CG scaffold pairs,
then fills free positions with A, C, G or U to reach exactly `ceil(n/2)`
GC bases (50% at these lengths). It shuffles pair-type and free-base
composition labels with a local Python `random.Random(sequence_seed)`.
Because scaffold endpoints change with shape, the embedded sequence also
changes; these are joint sequence/scaffold controls, not a same-sequence
causal comparison of scaffold shape.
