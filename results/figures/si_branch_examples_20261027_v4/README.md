# SI Figure S1: common arc geometry (2026-10-28)

Panels A, B and C now use the same 2.0 mm spacing per base and equal horizontal
and vertical coordinate scales. Every main-panel arc has a height of 0.45 times
its endpoint distance, close to a semicircle. The panel widths follow the sequence
lengths: 42, 52 and 70 mm. This widens A and B from their previous 40 mm panels
while giving all three examples the same scale.

A and B are stacked on the left, with C in the middle and its existing upper
and lower comparisons on the right. The canvas remains 165 × 65 mm, so the
paper's 165 mm inclusion width preserves all 6.6–7.6 pt text. Pair maps, colors,
line styles, shaded regions, labels and legend are unchanged. The v4 PNG was
opened and checked for arc shape, overlapping labels and clipping.

Reproduce from the experiments repository:

```bash
MPLCONFIGDIR=/tmp/si-branch-v4-mpl PYTHONDONTWRITEBYTECODE=1 \
  python3 scripts/make_si_branch_examples_20261027.py \
  --out results/figures/si_branch_examples_20261027_v4 \
  --context-reference ../../LinearKnot/test/tools/pseudoloop_context_reference.py
```

No RNA calculations were run. `examples.json` is byte-identical to v3, and the
read-only reference gives identical contexts for all five examples. Extracted
PDF text and font sizes also match v3. `validation.json` checks the pair maps,
marked contexts, label bounds, minimum font size and rendered coordinate scales;
it records 2.0 mm horizontal and vertical units and a 0.45 arc height/distance
ratio in A, B and C. `manifest.json` records the script and output hashes.

Only `si_branch_examples.pdf` is exported to the paper. The unchanged
inside/outside schematic is also rendered here; its existing paper copy is kept.
The renderer base commit is `1ba98bf4fc3d48c35da00bcba0423871b2575801`, with this
revision left uncommitted.

The paper's `bash build.sh` and `bash make_diff.sh 05ca09b` succeeded. Both the
normal SI and the SI diff against the preceding paper source have 18 pages;
Figure S1 on page 1 was opened and checked in both. SI logs have no errors,
undefined references/citations or overfull boxes. The SI text and caption are
unchanged, and v3 retains the preceding figure for visual comparison.

The shared project quota reported 20,747,353 of 28,000,000 files used before this
small revision. This directory contains eight files; rendering caches and PDF
review images stay under `/tmp`.
