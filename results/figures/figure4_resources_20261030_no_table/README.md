# Figure 4: remove completion-count table (2026-10-30)

Author feedback in the 2026-10-28 main manuscript requests removal of the table below panels A–D.
The completion counts remain in `completion_counts.tsv` for auditing; no completion-count text is drawn in the figure.
The canvas is 178 × 54 mm (previously 178 × 64 mm). Plot rectangles, legends, font sizes, curves, symbols and measured values retain their physical sizes.
The removed lower 10 mm contained the table and completion-count line.

Reproduce from the experiments repository:

```bash
python3 -B scripts/make_figure4_resources_20261026.py \
  --out results/figures/NEW_DIRECTORY
```

`plot_data.json` and `completion_counts.tsv` are byte-identical to `wording_regression_20261028/resources`.
The generated PNG was opened and checked: all four panels and both legends are visible without clipping or overlap.
The script checks canvas size, font sizes, text bounds, completion totals and source validation.
No RNA computation, input extraction or scheduler submission was performed.
`figure4.pdf` is exported to the paper as `figures/figure3_main.pdf`; the unused regenerated SI ratio figure is not exported.
The manuscript caption and Results describe CParty's failures; their source is `docs/comparison_resources_20261025.md`.
