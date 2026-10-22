# E1 accuracy figures — 89baea52, 2026-10-19

Finite-beam and exact outputs both use 89baea52. Long RNAs use the completed
`e1-bprna-long-89baea52-cpu16-20261018` run and the finite b=400 reference.
RNA cohorts, metric definitions and linear-interpolation quantiles are unchanged.
Round four makes small nonzero TV values visible on a zero-preserving scale.

Run from the experiments repository:
`python3 -B scripts/render_e1_89baea52.py --out results/figures/e1_89baea52_20261019_r4`
The output directory must not exist. The renderer does not change old figures.

`figure2_main` is 178 mm wide. The SI companions `figure2_accessibility_si`
(all four window lengths) and `figure2_d_si` (positive delta F threshold)
are 165 mm wide. Each has a PDF and 300-dpi PNG.
`caption.md` contains paired Japanese/English captions.
`plot_data.tsv` records all plotted values and `evaluation_counts.tsv` the
eligible RNA/window denominators. `validation.json` records checks against
the source quantiles, geometry, overlap and actual PDF font sizes.
`manifest.json` records sources, scripts and output hashes.
Round-four changes: B uses the same zero-preserving scale as C, with separate
filled/open metric symbols; the two symbols have wider horizontal separation.
The figure labels the linear regions in A and B/C. The distribution key combines
the range, percentile interval and a median symbol, with the text “symbol: median”.
S1 gives its own linear-region key, references and local limits in the caption.
S2 labels each group's denominator N using saved RNA counts.
Figure 2 remains 178 × 92 mm; S1 and S2 remain 165 × 137 and 165 × 55 mm.
The captions use the adopted Figure 2 / S1 / S2 numbering.
PNG inspection is recorded below after generation; manuscript insertion belongs
to the manuscript editor. No LaTeX or existing figure is changed.

## 2026-10-23 author feedback G51, G74–G76

Figure 2C replaces MAE with “absolute difference”; the main-text caption defines the per-RNA mean.
The plotted values and eligibility counts are byte-identical to the adopted 2026-10-20 figures.
Opened the new PNG and checked the expanded label and unchanged axes.
This is a label-only rendering from saved aggregates; no RNA calculation was run.
Exported figure2_main.pdf to the paper repository's figures/.
