# Figure 4: NEAT1 (2026-10-16)

[PDF](figure4.pdf) / [PNG](figure4.png) / [caption](caption.md) / [provenance](provenance.json)

The PDF is exactly 183 × 90 mm. The PNG is rendered at 400 dpi.
The script is [make_figure4.py](../../../scripts/make_figure4.py).
From the paper repository run:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python experiments/scripts/make_figure4.py
```

From the experiment repository, omit the `experiments/` prefix.
The script reads the fixed E3 analysis snapshot and streams only its recorded archive members.
It does not discover newly completed jobs, run RNA software, or extract archives.

## Panel sources

All table links below refer to [the E3 analysis snapshot](../../local/e3-neat1-20261016/README.md).
`provenance.json` records input and table hashes plus every consumed archive path, member, hash, and saved command.

| Panel | Direct data | Saved run and conditions |
|---|---|---|
| A | [NEAT1_2.txt](../../../data/neat1-20261016/NEAT1_2.txt), [input manifest](../../../data/neat1-20261016/manifest.json) | NR_131012.1, 22,743 nt; coordinates are one-based and inclusive. The whole-transcript track is schematic. |
| B | [bpp_pairs.tsv.gz](../../local/e3-neat1-20261016/bpp_pairs.tsv.gz), cross-checked against streamed `bpp.tsv` | `e3-neat1-bpp-e2cecbb4-20261016` and `e3-neat1-variants-bpp-e2cecbb4-20261016`, `neat1_long/000001`; NEAT1_2, threshold-0.7 G versus empty G, beam 100. |
| C | [position_profiles.tsv.gz](../../local/e3-neat1-20261016/position_profiles.tsv.gz), [accessibility.tsv.gz](../../local/e3-neat1-20261016/accessibility.tsv.gz), cross-checked against streamed CapR/Raccess stdout | Corresponding `e3-neat1-[variants-]capr/raccess-e2cecbb4-20261016` runs, `neat1_long/000001`; the same two conditions as B. |
| D | [sampling_histograms.tsv](../../local/e3-neat1-20261016/sampling_histograms.tsv), [sampling_summary.tsv](../../local/e3-neat1-20261016/sampling_summary.tsv), [sampling_reference_distances.tsv.gz](../../local/e3-neat1-20261016/sampling_reference_distances.tsv.gz), [representative_distances.tsv](../../local/e3-neat1-20261016/representative_distances.tsv) | `e3-neat1-sampling/bpp/mea-e2cecbb4-20261016` runs, `neat1_short/000001`; NEAT1_1 (NR_028272.1), threshold-0.7 G, beam 100/200, 1,000 samples per beam, seed 20261016, MEA accuracy weight 1. Every distance was recomputed from streamed structures. |

## Scientific scope and display choices

- B and C are local views of calculations on the full 22,743-nt transcript. They do not use a refolded fragment.
- Both bases of a displayed arc must lie in 20,350–20,600. B draws probabilities ≥0.01; the original BPP export cutoff is 0.001. For G, 17 scaffold pairs and 130 new pairs are drawn; for empty G, 217 pairs are drawn. 55 / 24 pairs at the display cutoff touch the interval but have their other base outside and are omitted. Contexts and accessibility still include full-transcript contributions.
- Gray, blue, and orange follow Figure 1: scaffold pairs, new pairs without G crossings, and new pairs crossing G. No orange arc passes the display cutoff here. The seven-context colors match Figure 1 exactly. No small PseudoLoop value is inflated for visibility.
- C shows the seven probabilities as per-base steps with no smoothing. Their sum is one. Hatch marks identify 8-nt windows containing at least one scaffold-paired base. Their exact accessibility is zero; microscopic raw numerical residuals are checked before assigning this exact condition.
- Every highlighted window starting at 20,451–20,500 contains only bases unpaired in G. The 50-window means are 0.526765463082 with G and 0.013738773627 with empty G. The 57 covered base positions run from 20,451 to 20,507. B/context boundaries cover these bases; the accessibility highlight covers the 50 starts.
- The focus was selected exploratorily for its large accessibility difference in the earlier E3 analysis; it is not a pre-specified biological region. Its coordinates were explicitly requested for this figure. No experimental validation or causal pseudoknot effect is inferred.
- D is a whole-transcript analysis of the shorter NEAT1_1 isoform, separate from the NEAT1_2 interval. Removing G prevents its 727 shared pairs from reducing the apparent diversity. Each beam has 1,000 distinct sampled structures and 499,500 unordered sample pairs. Pairs sharing a sample are not independent replicates.
- D uses 100 bins over [0,1], width 0.01, with density normalized separately to area one. The displayed distance axis is 0.30–0.95; all observed distances lie within it. Means come from exact distances, not bin centers. Beam 100 is solid and beam 200 dashed. MFE and MEA are separate reference structures for each beam; finite-beam MFE may use a different search set from the sampled ensemble.

## Development checks

Input hashes, archived member hashes and byte counts, saved beam/ranking/local-limit/seed/MEA settings, transcript identity, all per-base context probabilities, all width-8 accessibility values, BPP tables, and scaffold constraints are checked before plotting.
All 499,500 pairwise distances and 2,000 sample-to-reference distances per beam are independently recomputed with integer bitsets and compared with the saved tables.
Histogram counts and exact means must agree. The source analysis used sparse pair-incidence matrix products.
The script checks text bounds and overlapping text. Visual PNG review and the final PDF/PNG dimensions are recorded in the paper repository's `docs/progress.md`.
No manuscript file, bibliography, external implementation, RNA CLI, scheduler job, or Git commit is part of this script.

## 著者に選んでほしい点

1. 20,350–20,600 の周辺も含めた表示区間を本文に採用するか。この範囲は注目 window に加えて近傍の対と固定塩基を含む window を同時に示す。
2. D に beam 100/200 の両方を残すか、beam 100 のみとしてさらに簡潔にするか。現在は幅による分布の違いも小さく示せる両方を採用した。
