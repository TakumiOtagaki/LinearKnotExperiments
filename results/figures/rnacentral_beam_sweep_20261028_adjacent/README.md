# Figure S5 without a reference beam (2026-10-28)

165 × 122 mm, unchanged canvas. Base text 7 pt, ticks 6.5 pt, legend 6.7 pt.
A: ensemble free energy per base, -R_gas T log(Z_hat_b)/n, 49 points.
B: mean pairing total variation between adjacent widths, 42 points;
the x coordinate is the width after each increase (50, 75, 100, 200, 300, 400).
C/D: original wall time and peak memory, unchanged 49 observations per panel.
All seven RNAs appear, b=100 retains its guide and filled symbols.
No abbreviations for metrics. PNG opened and checked for labels and clipping.
No beam serves as the sweep's reference; h=d=max(b,30) varies with b.

Reproduce from the experiments repository:

```sh
python3 scripts/plot_beam_sweep_20261028.py \
  --source results/analysis/beam-sweep-rnacentral-20261028 \
  --summary-out results/analysis/beam-sweep-rnacentral-20261028-adjacent \
  --figure-out results/figures/rnacentral_beam_sweep_20261028_adjacent
```

Use fresh output directories. Source summaries and original raw calculations
are unchanged. No new RNA calculations.
