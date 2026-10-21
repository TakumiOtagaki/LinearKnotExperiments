# LinearKnot 4f97fd3 time and memory

The four PDF/PNG pairs replace the resource figure and its SI companions.
`figure3_main` retains the author-approved 178 x 55 mm layout.
`figure3_synthetic_all_cli_si` shows all five lengths for the seven main shapes;
`figure3_mode_memory_si` and `figure3_parallel_impact_si` show the stage and concurrency checks.
Only time/RSS measurements change. Accuracy and structure figures retain their adopted sources.
All measurement points use `results/local/e2-4f97fd3-20261021/`.
Paired Japanese/English captions are in `caption.md`; plot scalars and validation are retained.

Reproduce: `python3 -B scripts/plot_timing_4f97fd3.py` from experiments/.
