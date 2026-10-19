# E1 accuracy figures — 89baea52, 2026-10-19

Finite-beam and exact outputs both use 89baea52. Long RNAs use the completed
`e1-bprna-long-89baea52-cpu16-20261018` run and the finite b=400 reference.
RNA cohorts, metric definitions, linear-interpolation quantiles and the
`main_20261017/figure2_main` layout are unchanged.

Run from the experiments repository:
`python3 -B scripts/render_e1_89baea52.py`
The output directory must not exist. The renderer does not change old figures.

`figure2_main` is 178 mm wide. The SI companions `figure2_accessibility_si`
(all four window lengths) and `figure2_d_si` (positive delta F threshold)
are 165 mm wide. Each has a PDF and 300-dpi PNG.
`caption.md` contains paired Japanese/English captions.
`plot_data.tsv` records all plotted values and `evaluation_counts.tsv` the
eligible RNA/window denominators. `validation.json` records checks against
the source quantiles, geometry, overlap and actual PDF font sizes.
`manifest.json` records sources, scripts and output hashes.
PNG inspection and manuscript integration are recorded in `review_report.json`.
