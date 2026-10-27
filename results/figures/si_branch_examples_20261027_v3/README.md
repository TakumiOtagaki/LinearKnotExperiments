# SI Figure S1: shaded loop and comparisons (2026-10-27)

Panel C shades the loop face bounded by the outer and continuing inner pairs
of one band and an ordinary side branch. The side branch and inner region
remain white. Purple dots mark the unpaired segments of the shaded loop.
The upper inset removes the side branch while retaining the same two closing
pairs: the unpaired segments on both sides then belong to an internal loop.
The lower inset exchanges the scaffold and crossing-band roles, so the same
multiloop is closed by extension pairs. A side branch may also be a whole
pseudoknot branch, as the paper caption explains.

A and B retain their pair maps and branch labels. The main examples have
comparable arc heights. The canvas is 165 × 65 mm, matching the 165 mm paper
inclusion width. Text is 6.6–7.6 pt. The PNG and the normal and difference SI
pages were opened and checked for overlapping labels and clipping.

Reproduce from the experiments repository:

```bash
MPLCONFIGDIR=/tmp/si-branch-v3-mpl python3 scripts/make_si_branch_examples_20261027.py \
  --out results/figures/si_branch_examples_20261027_v3 \
  --context-reference ../../LinearKnot/test/tools/pseudoloop_context_reference.py
```

No RNA calculations were run. The script checks disjoint pair endpoints,
noncrossing scaffold and extension components, weakly closed branch
intervals, and that pairs of the other layer cross both closing pairs.
The read-only context reference classifies every marked unpaired base as
Multiloop in the main example and lower inset, and Internal in the upper
inset. `validation.json` also records minimum font sizes and checks that all
labels fit on the canvas. `examples.json` retains all five pair maps;
`manifest.json` records the script and output hashes.

Only `si_branch_examples.pdf` was exported to the paper. The script also
generates the inside/outside schematic; the paper's existing copy was kept.
The renderer base commit is `840b72fd70164c444cafa47107203fe2ec392b4d`, with
this revision left uncommitted.

Before this small revision, the shared project quota reported 20,586,911 of
28,000,000 files used. This directory contains eight small files; rendering
caches and PDF review images were kept under `/tmp`.
