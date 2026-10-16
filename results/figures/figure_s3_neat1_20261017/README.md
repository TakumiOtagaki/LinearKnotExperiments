# Supplementary Figure S3: NEAT1 (2026-10-17)

[PDF](figure_s3.pdf) / [PNG](figure_s3.png) / [caption](caption.md) / [provenance](provenance.json) / [validation](validation.json) / [focus values](focus_accessibility.tsv)

The PDF is exactly 165 × 145 mm, for inclusion at the SI width of 165 mm.
The PNG is 300 dpi.
Ordinary text is 7 pt; independent bold A–D labels are 9 pt.
No mathematical superscripts or reduced-size logarithmic ticks are used.
The wider two-column comparison and three bottom histograms replace the 183 × 90 mm Figure 4 layout without overwriting it.
C uses seven named probability tracks with colored areas and black outlines instead of a color-only stack.
Every track uses its original 0–1 scale; no small component is enlarged.

## Reproduce

From the paper repository:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python experiments/scripts/make_figure4.py --beam 200 --region 21000 21300 --focus 21134 21143 --out experiments/results/figures/figure_s3_neat1_20261017
```

From the experiment repository, omit the `experiments/` prefixes.
Requirements: Python, NumPy, Matplotlib, Pillow, pypdf, and Ghostscript (`gs`).
The default output is the new S3 directory.
The saved sampling results support main beam widths 100 and 200; this figure uses 200.
`--region` and `--focus` use one-based inclusive positions.
All focus windows must fit inside the view and contain only bases unpaired in G.
The script streams saved tar members without extracting them, calling RNA software, or changing the analysis snapshot.
Temporary font caches and PDF inspection files use node-local temporary storage.

## Data and version

The scientific outputs remain those of LinearCParty `e2cecbb4`, not the ongoing resource rerun.
The [fixed E3 analysis](../../local/e3-neat1-20261016/README.md) supplies all table values.
`provenance.json` records the script hash, selected beam and positions, input/table hashes, saved commands, and every consumed archive/member hash.
The fixed E3 snapshot predates empty-G b=200 BPP.
Only the successful `NEAT1_2_emptyG-beam200_outer-LinearCParty-rep0` record and its committed `stdout.txt` and `bpp.tsv` members are read from the saved variants-BPP checkpoint.
Their archive boundaries, byte counts, and SHA-256 hashes are verified.

| Panel | Input and saved source | Scope |
|---|---|---|
| A | [NEAT1_2 input](../../../data/neat1-20261016/NEAT1_2.txt), [manifest](../../../data/neat1-20261016/manifest.json) | NR_131012.1, 22,743 nt; schematic full-transcript position track. |
| B | BPP archives `e3-neat1-[variants-]bpp-e2cecbb4-20261016`, `neat1_long/000001`; [BPP table](../../local/e3-neat1-20261016/bpp_pairs.tsv.gz) | G at threshold 0.7 versus empty G, b=200; the independent BPP table lacks empty-G b=200, covered by the named checkpoint above. |
| C | Corresponding `capr`/`raccess` archives; [contexts](../../local/e3-neat1-20261016/position_profiles.tsv.gz), [accessibility](../../local/e3-neat1-20261016/accessibility.tsv.gz) | Whole-transcript probabilities displayed locally; 8-nt total accessibility, not accessibility by loop type. Saved beam/local-limit comparisons are numerical, with no hidden overlay. |
| D | `sampling`/`bpp`/`mea` archives, `neat1_short/000001`; [histograms](../../local/e3-neat1-20261016/sampling_histograms.tsv), [summaries](../../local/e3-neat1-20261016/sampling_summary.tsv), [reference distances](../../local/e3-neat1-20261016/sampling_reference_distances.tsv.gz), [representative distances](../../local/e3-neat1-20261016/representative_distances.tsv) | NEAT1_1, NR_028272.1, 3,756 nt; 1,000 samples, b=200, seed 20261016, MEA weight 1. |

## Scientific interpretation

B and C display local results of full-transcript calculations, not refolded fragments.
Only arcs with probability ≥0.01 and both endpoints in 21,000–21,300 are drawn; the export cutoff was 0.001.
G has 44 displayed scaffold pairs and 139 new non-crossing pairs; empty G has 165 displayed pairs.
Pairs with only one endpoint in view are counted once and omitted: 53 with G, 204 with empty G.
Their contributions remain in C.
The numbers above are recomputed from the verified BPP exports.
Solid gray (#8a8a85) identifies G, dashed blue (#2a78d6) identifies new pairs without G crossings, and dotted orange (#eb6834) identifies pairs crossing G.
There are 0 displayed crossing pairs with G.
The seven context colors follow the shared palette: Stem #2a78d6, Hairpin #1baf7a, Bulge #eda100, Internal #e87ba4, Multiloop #6250d6, Exterior #b9b8b1, PseudoLoop #eb6834.
Direct row labels and black probability outlines allow identification without color.
A `<0.001 everywhere` label is used only when every saved value in that track is below 0.001.

The 10 highlighted starts are 21,134–21,143; their windows cover 17 bases through 21,150.
Dashed boundaries in B/context rows enclose the bases; the pale accessibility shading encloses the starts.
Every highlighted window contains only bases unpaired in G.
Hatching in the accessibility graph and legend uses the same `////` pattern for windows containing any G-paired base.
Such windows have exact accessibility zero after checking their raw numerical residuals.

The focus was selected exploratorily after inspecting the saved predictions and their width sensitivity; it is not a pre-specified biological region.
The author's accepted focus and b=200 setting are unchanged.
Changing G also changes the conditional ensemble; this comparison does not isolate a causal pseudoknot effect or experimentally validate the predictions.

## Focused accessibility and sensitivity

| Scaffold | Beam width | Mean over the 10 highlighted windows |
|---|---:|---:|
| G, threshold 0.7 | 200 | 0.785137595592449 |
| G, threshold 0.7 | 400 | 0.7851375955939182 |
| Empty G | 100 | 3.103196239844372e-05 |
| Empty G | 200 | 0.01100475331979723 |

The maximum per-window change is 1.555977569012157e-12 for G (b=200/400) and 0.03340806386802516 for empty G (b=100/200).
All four series are saved in `focus_accessibility.tsv`.
Only the highlighted windows are claimed to be stable within 0.05 across these saved settings.
The local hairpin and case-3 depth limits each equal `max(b, 30)` in these runs; the comparison measures beam/local-limit sensitivity, not agreement with exact calculations.

D removes 727 common G pairs before computing Jaccard distances.
It uses all 499,500 unordered sample pairs and 1,000 distances to each reference structure.
Pairs sharing a sample are not independent replicates.
Each histogram uses 100 bins of width 0.01 over [0, 1] and is independently normalized to area one.
All observed distances lie within the displayed 0.30–0.95 range.
Dashed vertical lines and labels use exact means, not bin centers: 0.685568009553, 0.604910056350, and 0.577992608021.
Finite-beam MFE may search a different structure set from the sampled ensemble.

## Checks and review coverage

- S3-1 / G1: 165-mm native width; two 67-mm data columns; D spans the bottom; ordinary text 7 pt.
- S3-2: omitted outside-endpoint pair counts in B; caption explains full-transcript contributions in C.
- S3-3: means explicitly refer to highlighted starts; base and start boundaries are defined separately.
- S3-4: largest per-window changes replace the obscured b=400 overlay; all saved-width means remain in the source table.
- S3-5: exact mean lines/labels, distance counts, independent density normalization, bin width, and MEA weight are explicit.
- S3-6: graph and legend use identical hatching.
- G2: named context rows with black outlines; pair line styles distinguish G and new pairs.
- G3: A–D use independent 9-pt bold letters at the same left position.

Before plotting, the script verifies input and artifact hashes, byte counts, archive boundaries, saved beam/ranking/local-limit/seed/MEA settings, and scaffold constraints.
Full-transcript context and accessibility arrays must match the independently saved tables.
Available BPP tables must match; the supplemental empty-G BPP must also satisfy incident-mass and context bounds.
All 499,500 sample-pair and 2,000 sample-reference distances are independently recomputed with integer bitsets and checked against the saved counts and exact means.
Text bounds, pairwise text overlaps, embedded PDF font sizes, PDF page dimensions, and PNG dimensions/DPI are checked in `validation.json`.
PNG visual inspection must additionally be recorded in the task log.
No manuscript, bibliography, external implementation, scheduler job, or Git commit is modified by this script.
The final placement and caption insertion in the SI belong to the manuscript editor.

## Final rendering review (2026-10-17)

The final 300-dpi PNG was opened at its original resolution and checked for text overlap, clipping, and readability.
The final PDF has exactly 165 × 145 mm pages, with embedded text only at 7.0 and 9.0 pt (100 text spans); no glyph is below 6 pt.
A node-local `article` 11-pt / 1-inch-margin wrapper included this PDF at 165 mm with the complete proposed English caption.
It compiled to one page without oversized-float or overfull-box warnings; its rendered page was also opened.
The minimum embedded text size after inclusion remained 7.0 pt.
This checks the isolated figure and caption, rather than the concurrently edited complete SI.

Thirteen scientific provenance fields exactly match the previous figure, including all focus values, context/accessibility diagnostics, pair counts, source tables, consumed artifacts, and recomputed distance summaries.
The old figure directory's six files retain their hashes.
`validation.json` records the final review and the tested PDF/caption hashes.
The complete SI must still be rebuilt by its editor after inserting the figure path and caption.
The shared `docs/progress.md` was left to the coordinating editor because this task restricts edits to the assigned files; the task log records the same changes and checks.
