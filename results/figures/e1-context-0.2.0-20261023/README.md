# E1 context update — LinearKnot 0.2.0, 2026-10-23

Profiles use f93b55e (v0.2.0); all other E1 metrics use the certified saved
89baea52 outputs. Cohorts, definitions, reference modes and rendering functions
are unchanged. See manifest.json for the exact source and script hashes.
The SI companions use total accessibility and free energy; neither uses contexts.

Reproduce from the experiments repository:
```sh
python3 -B scripts/render_e1_89baea52.py \
  --short-data results/local/e1-context-0.2.0-20261023/e1-bprna-short-context-0.2.0-20261023 \
  --long-data results/local/e1-context-0.2.0-20261023/e1-bprna-long-context-0.2.0-20261023 \
  --out results/figures/e1-context-0.2.0-20261023
```

PNG inspection (2026-10-23): Figure 2 axes, legends, open context symbols and labels are clear, with no overlap or clipping. Non-context plot rows and evaluation counts match the adopted figure; both SI PNGs are byte-identical. Figure 2 was exported to the paper figures directory.
