# SI branch examples and inside/outside (2026-10-27)

Two new schematic figures for the author's SI comments on `supplementary_1027_1110.pdf` (d2302bf). No RNA sequence, energy, probability, or timing was calculated.

- `si_branch_examples.pdf` (165 × 55 mm): ordinary V branch, whole WMB pseudoknot branch, and a multiloop between two pairs of a crossing band. Grey solid scaffold pairs, blue dotted new noncrossing pairs, orange dashed new crossing pairs, and purple Multiloop bases match `scripts/make_figure1.py`.
- `si_inside_outside.pdf` (152 × 54 mm): original binary decomposition example, with upward inside arrows and downward outside arrows. The outside contribution to one child uses the other child's inside weight. Squares denote rule applications; rounded boxes denote states. This illustrates the general interpretation in Eisner (2016), [official paper](https://aclanthology.org/W16-5901.pdf), especially Section 6. The layout and example were drawn independently; no figure was copied.

Reproduce from the experiments repository:

```bash
MPLCONFIGDIR=/tmp/si-branch-mpl python3 scripts/make_si_branch_examples_20261027.py \
  --out results/figures/si_branch_examples_20261027 \
  --context-reference ../../LinearKnot/test/tools/pseudoloop_context_reference.py
```

The context-reference path is configurable and read only. The reference file at LinearKnot v0.2.1 has SHA-256 `30547e8701d3c28b7e56e97384b207d405157e0c6df4b7cf71d605892740a48f`. Its context labels are recorded in `validation.json`. All purple bases in C are Multiloop. Each component is noncrossing, pair endpoints are disjoint, each pair has at least three intervening bases, and the branch intervals are weakly closed in G. `examples.json` records all coordinates; `manifest.json` records renderer and artifact hashes. The displayed examples are fictional, rather than predictions for biological RNAs.

The shared quota was checked before drawing: 21,919,432 of 28,000,000 files used. This task creates two PDFs, two previews, and a few small provenance files, with temporary rendering under `/tmp`. Existing figures and run archives were preserved.

日本語：通常枝は V の区間、pseudoknot 枝は隙間の構造を含む WMB の全区間である。Band をまたぐ multiloop は一つの band の外側・内側の対の間にあり、側枝を持つ。紫の非対合塩基は全て Multiloop に分類される。Inside/outside 図は独自の分解例で、前者は子から親へ、後者は根から子へ伝わり、子の outside には他の子の inside が必要になる。
