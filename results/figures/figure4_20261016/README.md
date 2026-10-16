# Figure 4: NEAT1 (2026-10-16)

[PDF](figure4.pdf) / [PNG](figure4.png) / [caption](caption.md) / [short caption](caption_short.md) / [provenance](provenance.json)

The PDF is exactly 183 × 90 mm. The PNG is rendered at 400 dpi.
The script is [make_figure4.py](../../../scripts/make_figure4.py).
From the paper repository run:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python experiments/scripts/make_figure4.py --beam 200 --region 21000 21300 --focus 21134 21143
```

From the experiment repository, omit the `experiments/` prefix.
The `--beam` option selects the main width for B, C, and D; its default is 200.
The saved sampling results allow 100 or 200. Use `--out` for an alternative figure directory.
The `--region START END` option sets the display interval (default: 21,000–21,300).
The `--focus START END` option sets the highlighted window starts (default: 21,134–21,143).
Both ranges are one-based and inclusive; all highlighted window bases must fit within the display interval and be unpaired in G.
The author selected the new interval and beam 200 on 2026-10-16 after checking the accessibility changes at other widths.
The pale-blue G curve in C always uses b=400; empty-G b=100 is checked numerically and is not drawn.
The script reads the fixed E3 analysis snapshot and streams recorded archive members.
That snapshot predates the empty-G b=200 BPP result.
For this one missing result, the script reads only the named `NEAT1_2_emptyG-beam200_outer-LinearCParty-rep0` job from the saved BPP checkpoint.
It verifies success, archive boundaries, byte counts, and hashes before using its `stdout.txt` and `bpp.tsv`.
It does not scan other runs, run RNA software, extract archives, or update the E3 snapshot.

## Panel sources

All table links below refer to [the E3 analysis snapshot](../../local/e3-neat1-20261016/README.md).
`provenance.json` records the selected coordinates and beam, input and table hashes, every consumed archive path, member, hash, and saved command, and the supplemental checkpoint snapshot and selected job.
It also records the four focused accessibility series and means, the per-window width comparisons, and the additional b=400 curve.

| Panel | Direct data | Saved run and conditions |
|---|---|---|
| A | [NEAT1_2.txt](../../../data/neat1-20261016/NEAT1_2.txt), [input manifest](../../../data/neat1-20261016/manifest.json) | NR_131012.1, 22,743 nt; coordinates are one-based and inclusive. The whole-transcript track is schematic. |
| B | Streamed `bpp.tsv`; [bpp_pairs.tsv.gz](../../local/e3-neat1-20261016/bpp_pairs.tsv.gz) provides an independent table cross-check except for empty G at b=200, which uses the supplemental checkpoint above. | `e3-neat1-bpp-e2cecbb4-20261016` and `e3-neat1-variants-bpp-e2cecbb4-20261016`, `neat1_long/000001`; NEAT1_2, threshold-0.7 G versus empty G, beam 200. |
| C | [position_profiles.tsv.gz](../../local/e3-neat1-20261016/position_profiles.tsv.gz), [accessibility.tsv.gz](../../local/e3-neat1-20261016/accessibility.tsv.gz), cross-checked against streamed LinearPKCapR/LinearPKRaccess stdout | Corresponding `e3-neat1-[variants-]capr/raccess-e2cecbb4-20261016` runs, `neat1_long/000001`; the same two conditions and main beam as B, plus G b=400 accessibility (drawn) and empty-G b=100 accessibility (numerical check). |
| D | [sampling_histograms.tsv](../../local/e3-neat1-20261016/sampling_histograms.tsv), [sampling_summary.tsv](../../local/e3-neat1-20261016/sampling_summary.tsv), [sampling_reference_distances.tsv.gz](../../local/e3-neat1-20261016/sampling_reference_distances.tsv.gz), [representative_distances.tsv](../../local/e3-neat1-20261016/representative_distances.tsv) | `e3-neat1-sampling/bpp/mea-e2cecbb4-20261016` runs, `neat1_short/000001`; NEAT1_1 (NR_028272.1), threshold-0.7 G, beam 200 only, 1,000 samples, seed 20261016, MEA accuracy weight 1. Every distance was recomputed from streamed structures. |

## Scientific scope and display choices

- B and C are local views of calculations on the full 22,743-nt transcript. They do not use a refolded fragment.
- Both bases of a displayed arc must lie in 21,000–21,300. B draws probabilities ≥0.01; the original BPP export cutoff is 0.001. For G, 44 scaffold pairs and 139 new pairs are drawn; for empty G, 165 pairs are drawn. 53 / 204 pairs at the display cutoff touch the interval but have their other base outside and are omitted. Contexts and accessibility still include full-transcript contributions.
- Gray, blue, and orange follow Figure 1: scaffold pairs, new pairs without G crossings, and new pairs crossing G. No orange arc passes the display cutoff here. The seven-context colors match Figure 1 exactly. No small PseudoLoop value is inflated for visibility.
- C shows the seven LinearPKCapR probabilities as per-base steps with no smoothing. Their sum is one. LinearPKRaccess supplies the 8-nt accessibility values. Hatch marks identify windows containing at least one scaffold-paired base. Their exact accessibility is zero; microscopic raw numerical residuals are checked before assigning this exact condition.
- Every highlighted window starting at 21,134–21,143 contains only bases unpaired in G. The 10-window means at b=200 are 0.785137595592 with G and 0.011004753320 with empty G. The 17 covered base positions run from 21,134 to 21,150. B/context boundaries cover these bases; the accessibility highlight covers the 10 starts.
- The pale-blue G b=400 curve is drawn behind the black main accessibility curve with lower opacity and a slightly wider line. The line key distinguishes the widths. Agreement is checked per highlighted window, not just between means; changes elsewhere in the displayed interval are not claimed to be small.
- The focus was selected exploratorily after comparing accessibility and checking its sensitivity to the saved widths; it is not a pre-specified biological region. The previous focus at 20,451–20,500 had nearly equal means at b=200 (0.0136 / 0.0137); its earlier b=100 contrast reflected width sensitivity with G. The author explicitly requested the replacement. No experimental validation or causal pseudoknot effect is inferred.
- D is a whole-transcript analysis of the shorter NEAT1_1 isoform, separate from the NEAT1_2 interval. Removing G prevents its 727 shared pairs from reducing the apparent diversity. The selected beam has 1,000 distinct sampled structures and 499,500 unordered sample pairs. Pairs sharing a sample are not independent replicates.
- D uses 100 bins over [0,1], width 0.01, with density normalized separately to area one. The displayed distance axis is 0.30–0.95; all observed distances lie within it. Means come from exact distances, not bin centers. Only beam 200 is plotted, with a solid outline and gray fill. The mean distances are 0.685568009553 (sample–sample), 0.604910056350 (sample–MFE), and 0.577992608021 (sample–MEA). Finite-beam MFE may use a different search set from the sampled ensemble.

## Focused accessibility and beam-width sensitivity

| Scaffold | Beam width | Mean over the 10 highlighted windows |
|---|---:|---:|
| G, threshold 0.7 | 200 | 0.785137595592449 |
| G, threshold 0.7 | 400 | 0.7851375955939182 |
| Empty G | 100 | 3.103196239844372e-05 |
| Empty G | 200 | 0.01100475331979723 |

For G, the largest absolute change among the highlighted windows between b=200 and b=400 is 1.555977569012157e-12; the absolute change of the mean is 1.46915812848647e-12.
For empty G, the corresponding changes between b=100 and b=200 are 0.03340806386802516 and 0.01097372135739879.
These comparisons describe sensitivity to the saved beam/local-limit settings (`--max-hairpin` and `--max-case3-depth` each equal b), rather than agreement with an exact calculation.

## Development checks

Input hashes, archived member hashes and byte counts, saved beam/ranking/local-limit/seed/MEA settings, transcript identity, all per-base context probabilities, all width-8 accessibility values, available BPP tables, and scaffold constraints are checked before plotting.
The supplemental empty-G BPP also passes probability bounds and incident-mass checks against the full LinearPKCapR stem probabilities; no absent table comparison is claimed.
All 499,500 pairwise distances and 2,000 sample-to-reference distances at the selected beam are independently recomputed with integer bitsets and compared with the saved tables.
Histogram counts and exact means must agree. The source analysis used sparse pair-incidence matrix products.
The script checks text bounds and overlapping text. Visual PNG review and the final PDF/PNG dimensions are recorded in the paper repository's `docs/progress.md`.
No manuscript file, bibliography, external implementation, RNA CLI, scheduler job, or Git commit is part of this script.

The author's 2026-10-16 choices settle the display interval and main beam selection for this figure.
