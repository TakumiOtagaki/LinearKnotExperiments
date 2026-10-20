# Accuracy figures — author feedback, 2026-10-20

The annotated Figure 2 (main p5) and SI accuracy plot (p20) contain no direct figure comments.
The new paths are under `experiments/results/figures/e1_89baea52_20261020_fb/`.
Use `figure2_main.pdf` and `figure2_accessibility_si.pdf` (S1), with the paired caption proposals in caption.md.
`figure2_d_si.pdf` is the existing extra threshold companion; the manuscript already omits this duplicate plot. The synthetic scaffold plot is now S2 after NEAT1 figure removal.
The caption names beam widths instead of unspecified settings, defines h as the new-hairpin length limit, and defines d as the number of enclosing scaffold multiloop levels searched outward.
No old CLI names occur in the image labels or caption. Numerical data and plot layout are unchanged.
`plot_data.tsv` and `evaluation_counts.tsv` are byte-identical to the previous r4 figure.
Finite-beam and exact calculations both use 89baea52. Long RNAs use the completed `e1-bprna-long-89baea52-cpu16-20261018` run and the finite b=400 reference.
The RNA cohorts, error definitions and linear-interpolation quantiles are unchanged.
Figure 2 is 178 × 92 mm; S1 is 165 × 137 mm; the extra threshold plot is 165 × 55 mm. PNGs are 300 dpi.
`validation.json` records 1,750 checked quantile scalars, 350 checked summaries, plot geometry, text overlap and embedded PDF font checks.
`manifest.json` records source/output/script hashes; `evaluation_counts.tsv` gives the eligible RNA and window counts.
Regenerate: `python3 -B scripts/make_figures_feedback_20261020.py --accuracy-only --output-root NEW_DIRECTORY`.
The E1 renderer refuses an existing output directory.

## Final visual and font review — 2026-10-20

Opened `figure2_accessibility_si.png`, `figure2_d_si.png`, `figure2_main.png` and checked panel labels, legends, line styles, shading, and text placement.
All PDF text is at least 7 pt at the stated physical widths; PNGs are 300 dpi. Text bounds/overlap checks passed.
