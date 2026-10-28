# Figure 4: completion counts by length group (2026-10-28)

Canvas 178 × 64 mm (previously 178 × 58 mm). Plot rectangles and legends keep
their physical sizes; all plotted resource values and ranges are unchanged.
plot_data.json is byte-identical to ../figure4_review_stage1_20261027/plot_data.json.
Six millimeters of extra height leave room between the legends and count table.
Text is 7–9 pt; count rows use method colors and 7.5 pt.

The unfinished-input rugs and labels are replaced by a shared numeric table
for C/D. Columns are the five bpRNA length-group centers (250, 500, 1,000,
2,000, 4,000 nt), followed by the seven individual RNAcentral lengths.
CParty: 20/20, 20/20, 19/20, 0/20, 0/20, then seven 0/1.
Spark: 20/20, 19/20, 20/20, 20/20, 20/20, then four 1/1 and three 0/1.
LinearKnot: all 107 inputs completed at b=100 for each output group; all 80
planned exact inputs completed. Longer exact inputs were not planned.
Counts require all three measured runs to succeed, as for the plotted values.
See completion_counts.tsv and docs/comparison_resources_20261025.md.

The PNG and the vector PDF were opened and inspected for clipping, overlap, colors and counts.
No new RNA calculations. The unused synthetic-ratio companion is not exported.

Reproduce from experiments/ with a fresh output directory:

```sh
python3 -B scripts/make_figure4_resources_20261026.py \
  --out results/figures/figure4_resources_20261028_counts
```
