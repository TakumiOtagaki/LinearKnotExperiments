# SI Figure S3: ratio wording (2026-10-26)

The axis now reads “Cost at 4,000 nt divided by cost at 2,000 nt”.
The data are byte-identical to ../figure4_resources_20261026_si_labels/plot_data.json.
Dimensions (165 × 42 mm), font sizes (at least 7 pt), points, scales, and legend are unchanged. The PNG was visually checked. No RNA calculations were rerun.

Reproduce from experiments into a fresh directory:

```sh
python -B scripts/make_figure4_resources_20261026.py --out /tmp/fb7_si_ratio_render
```

The renderer also makes Figure 4, which is not exported by this change. Only SI Figure S3 and its data and checks are kept here.
