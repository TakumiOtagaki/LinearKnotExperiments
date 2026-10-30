# Figure 4: complete final legend (2026-10-30 21:20 feedback)

The last line of the right legend was clipped in the previous PDF.
Add 2 mm of bottom margin: the canvas is 178 × 56 mm instead of 178 × 54 mm.
Plot rectangles, legend positions relative to the plots, font sizes, colors,
markers and measured values keep their physical sizes. Completion counts
remain an audit artifact and are not drawn in the figure.
The script checks text bounds at the PNG export resolution (300 dpi).

Reproduce from the experiments repository:

```bash
python3 -B scripts/make_figure4_resources_20261026.py \
  --out results/figures/NEW_DIRECTORY
```

`plot_data.json` and `completion_counts.tsv` are byte-identical to
`figure4_resources_20261030_no_table`. The PNG and standalone PDF were opened
and checked: the complete two-line final legend is visible with bottom margin.
All PDF glyph bounds lie inside the canvas. Figure 4 on page 7 of both the
built main manuscript and the diff against `cb1e1db` was opened and checked;
the final legend is complete and clear of the caption. See `validation.json`.
`figure4.pdf` is exported as the paper's `figures/figure3_main.pdf`.
The regenerated SI ratio figure is not exported.
No RNA computation or scheduler submission was performed.
