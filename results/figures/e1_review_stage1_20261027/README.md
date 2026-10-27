# Review stage 1: Figure 2 and Figure S3 (2026-10-27)

Long length-group labels use two lines; Figure S3 shows length groups only on its bottom row.
Figure 2 has one shared horizontal-axis title.
Text sizes, canvas dimensions, measurements, cohorts, and metric definitions are unchanged.
`plot_data.tsv` is byte-identical to `../e1-context-0.2.0-20261023/plot_data.tsv`.
Generated PNGs were opened and checked; the renderer also checks text overlap and page bounds.

Reproduce from experiments/ with absolute paths for both data arguments:

```sh
python3 -B scripts/render_e1_89baea52.py \
  --short-data "$PWD/results/local/e1-context-0.2.0-20261023/e1-bprna-short-context-0.2.0-20261023" \
  --long-data "$PWD/results/local/e1-context-0.2.0-20261023/e1-bprna-long-context-0.2.0-20261023" \
  --out results/figures/NEW_REVIEW_E1
```

Source and output hashes are in `manifest.json`. No RNA calculation was run.
