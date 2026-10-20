# Time and memory — author feedback, 2026-10-20

Figure 4: `figure3_main.pdf` / `.png`, 178 × 55 mm.
Figure S2: `figure3_synthetic_all_cli_si.pdf` / `.png`, 165 × 178 mm.
Use the new paths under `experiments/results/figures/main_89baea52_20261020_fb/`.
Proposed paired captions are in [caption.md](caption.md).
Source: `results/local/e2-89baea52-20261018/`; source revision is 89baea52.
The updated names describe the matching LinearKnot subcommands; these are the saved 89baea52 timings, not new timings for the renamed CLI.
The 9c7d6e1 output comparison is documented in `docs/linearknot_9c7d6e1_20261020/README.md`.

F57–F64, F74, F88–F89: replace CLI names by output roles, explain b/h/d and the range of three measurements, define seeds as sequence-generation seeds, and name the synthetic scaffold shapes.
The two multiloop designs are labelled Nested multiloops and Balanced multiloops.
No measured value was changed; `plot_data.tsv` matches the previous r4 figure exactly.
The old crossed-out caption sentences about bounds and causal isolation are omitted.
`validation.json` checks source quantiles, plot geometry, overlapping text and actual PDF text size.
Run `python3 -B scripts/make_figures_feedback_20261020.py --resources-only` from experiments/.

## Final visual and font review — 2026-10-20

Opened `figure3_main.png`, `figure3_synthetic_all_cli_si.png` and checked panel labels, legends, line styles, shading, and text placement.
All PDF text is at least 7 pt at the stated physical widths; PNGs are 300 dpi. Text bounds/overlap checks passed.
