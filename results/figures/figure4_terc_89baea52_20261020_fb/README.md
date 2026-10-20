# Figure feedback — 2026-10-20

New manuscript figure: `figure4.pdf` (178 mm wide); preview: `figure4.png` (300 dpi).
Proposed paired Japanese/English caption: [caption.md](caption.md).
The manuscript editor should replace the old image path with `experiments/results/figures/figure4_terc_89baea52_20261020_fb/figure4.pdf` and use that caption.

Feedback addressed: F27, F40, F44–F52.
Source data: `results/local/figure-feedback-20261020/data.json`.
Raw CLI output, frozen inputs and saved Terc output are packed in the adjacent `raw_outputs.tar`.
Commands, parameter/binary/source hashes and checks are in `calls.json` and `provenance.json` there.
Figure 1 was computed with LinearKnot 9c7d6e1 (exact, Turner04, 37°C); Terc probabilities and 1,000 samples are saved 89baea52 outputs.
The newly computed γ-centroid uses `linearknot centroid --gamma 1`; Terc uses beam/hairpin/depth 200 and exact for comparison.
These small local calculations are not performance measurements.
All pair-probability arc plots use p ≥ 0.05; complete probability data remain in the source archive.
When inserting the figure, use the caption in this directory and the pair-probability threshold 0.05.
Internal is a plain pink fill; no hatching is used.
`validation.json` checks actual PDF text sizes, physical dimensions, text overlap and source-data identity.
The old D pair table and MFE/MEA counts are removed. The figure is now 178 × 104 mm (previously 178 × 100 mm). Exact comparisons near 40–160 are in data.json: context max difference 0.00283443, 4-nt total accessibility 0.00191013. The stored centroid output supplies the inset. At both exact and b200, it contains all six outer reference pairs and none of the three inner reference pairs; these replace the old MEA counts.

Regenerate from experiments/: `python3 -B scripts/make_figures_feedback_20261020.py --examples-only` (use `--output-root NEW_DIRECTORY` to preserve this rendering).
To regenerate the raw small calculations in a new location: `python3 -B scripts/prepare_feedback_figure_data.py --out NEW_DIRECTORY`.
No LaTeX, bibliography or external implementation tree was edited.

## Final visual and font review — 2026-10-20

Opened `figure4.png` and checked panel labels, legends, line styles, shading, and text placement.
All PDF text is at least 7 pt at the stated physical widths; PNGs are 300 dpi. Text bounds/overlap checks passed.
