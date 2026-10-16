# Revised main Figures 2/3 and SI companions — 2026-10-17

E1: beam source d21618c5; exact reference d9804df6.
E2: **d02d261c3a07e7308a3a3d37543d040adfd0f86e** (one revision for RNAcentral, synthetic controls and mode memory).
The resource figure keeps its task filename `figure3_main`; plan B may number it Figure 4 in the manuscript.
No manuscript, Japanese draft, bibliography or external implementation is edited.

## Files and validation

Every figure has a PDF and 300-dpi PNG.
Main figures are exactly 178 mm wide; SI figures are exactly 165 mm wide.
Normal text is 7 pt and independent bold panel letters are 9 pt.
Ghostscript checks the actual PDF text spans, including mathematical text.
All logarithmic ticks use ordinary decimal digits.

| PDF | Size (mm) | Minimum PDF text (pt) |
|---|---|---:|
| figure3_main.pdf | 178 × 55 | 7 |
| figure3_mode_memory_si.pdf | 165 × 65 | 7 |
| figure3_synthetic_all_cli_si.pdf | 165 × 178 | 7 |
| figure3_parallel_impact_si.pdf | 165 × 72 | 7 |

`caption.md` contains paired Japanese/English caption drafts, one sentence per line.
`plot_data.tsv` retains all checked E1 quantiles (including quartiles), threshold counts, E2 summaries and the exact numerator/denominator of every ratio.
When E1 is rendered, `evaluation_counts.tsv` records each accessibility domain, including excluded RNAs.
`validation.json` records source checks, geometry checks and PDF font sizes.
`manifest.json` records source/script/output hashes.
PNG visual checks and manuscript integration are recorded in `review_report.json` after inspection.

## Review mapping and choices

- F2-1: colors, symbols and within-group beam positions stay fixed in A–C and S1; length groups form the horizontal categories.
- F2-2: main C separates all length groups at 16 nt; the additional SI figure shows all four windows, with a declared linear band and no artificial positive floor.
- F2-3: A uses full range, 10th–90th percentiles and median; quartiles remain in the audit table.
- F2-4: pairing is filled and contexts are open, with a local key; caption defines TV and RNA averaging before the percentile.
- F2-5: beam values are settings within categorical length groups, not numeric line slopes; no lines connect the settings.
- F2-6: b=400 self-comparisons are omitted in all revised accuracy panels and explained in the captions.
- F3-1/F3-3: C shows the observed 4,000/2,000-nt ratio for each named series; all full curves and saved scaffold sketches move to a 165-mm SI companion.
- F3-2: larger open square and smaller filled triangle expose nearly coincident programs without moving the RNA coordinates; longest-input annotations are computed from aggregate medians.
- F3-4: program legend lists requested outputs; caption states whole-program measurement and distinct outputs.
- S1-1/S1-2: full threshold in the axis label; A/B, explicit denominators below, one-sided criterion and omitted self-comparison.
- S2-1/S2-2: explicit independent-process names, preliminary ranking, one measurement, points without connecting lines, and data-derived y limits.
- G1: final physical widths, 7-pt normal labels and decimal logarithmic ticks; embedded PDF font verification.
- G2: beam and program shapes, filled/open metric keys and direct synthetic-series labels supplement color; the saved scaffold sketches use G grey #8a8a85.
- G3: independent bold 9-pt panel letters at fixed offsets; Figure 2 uses one A/B/C per row and explicit reference columns.

The author/manuscript editor must choose final supplementary numbering for the two added detail figures and update figure references after plan B.
The E1 total-accessibility benchmark and historical E2 outputs do not validate the new per-loop accessibility components.
Resource values for the new per-loop option are not claimed.

## Numerical sources

| ID | Path relative to this directory | SHA-256 |
|---|---|---|
| e2_metrics.tsv | `../../local/e2-d02d261c-20261017/metrics.tsv` | ec5353edfaadbba1ccfb676b217f5afd09e5b2d1fc2d022fac90f4e6ccd410df |
| e2_comparisons.tsv | `../../local/e2-d02d261c-20261017/comparisons.tsv` | aee3aec54a95b834a92e32dfae94a0b2d976d07ba1b0754b18c673c368aad462 |
| e2_memory_metrics.tsv | `../../local/e2-d02d261c-20261017/memory_metrics.tsv` | adea6a83ddd8415dec8f938e5aac1664945abd005f216d0b7e90980d4eae58d8 |
| e2_summary.json | `../../local/e2-d02d261c-20261017/summary.json` | c5e67403c8b78bdda11e4f290db5f0128bebe1bad03c066e5a4e8c88231da6cd |
| g_features | `../../local/e2-e2cecbb4-g-features-20261016/rnacentral_g_features.tsv` | 856c1f56b79ef524038d1885540203ecdf9f037b5946592f0179e12032978c82 |
| synthetic_designs | `../../../data/e2-synthetic-shapes-20261012/inputs.jsonl` | acdab9208952bd2a53ce6a1cfce84f50e024dfcc707158ef1d949071027c8fbe |

E1 uses Turner 2004, 37°C, preliminary outer ranking and h=d=max(b,30).
All quantiles weight each evaluable RNA equally and use NumPy's linear quantile interpolation.
E2 uses b=h=d=100, outer ranking, Turner 2004, 37°C, one thread/process, at most 10 concurrent processes on cpu_40, and RLIMIT_AS=8 GiB.
RSS GiB = GNU time KiB / 2^20.
Synthetic ratios divide the medians already shown in the full curves; they are not medians of seed-wise ratios.

| Synthetic series | Time ratio | RSS ratio |
|---|---:|---:|
| Empty G | 2.05662 | 2.26857 |
| Long stem | 1.96016 | 2.70748 |
| 3 helices | 2.02297 | 2.56133 |
| Free bases together | 2.08745 | 2.50568 |
| 12 helices | 2.43318 | 2.53010 |
| Deep branching | 1.74273 | 3.52519 |
| Balanced branching | 1.65274 | 3.56919 |

**Integration check:** the time-ratio range from the plotted medians is 1.65–2.43.
If the manuscript retains a different range from medians of seed-wise ratios, update the text or explicitly distinguish the two statistics.

## Accessibility denominators

| Reference cohort | Length group (nt) | Window (nt) | Evaluated RNA | Excluded RNA | Eligible windows |
|---|---|---:|---:|---:|---:|


## Regenerate

Run from the paper repository:

```bash
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR="$TMPDIR/e2-figure-mpl" python3 experiments/scripts/make_figures.py --review-layout --figure3-only --e2-version d02d261c --out experiments/results/figures/main_d02d261c_20261017
```

`--e2-version febd2bc1` selects `results/local/e2-febd2bc1-20261016`.
`--e2-version d02d261c` selects `results/local/e2-d02d261c-20261017`.
For arbitrary aggregates, omit the preset and supply `--e2 PATH --e2-revision COMMIT`.
Historical defaults, `--main-compact`, `--figure3-only` and `--figure3-panels` remain available.
Use a new output directory for a different E2 revision, for example `main_20261017/d02d261c/`, to preserve this febd2bc1 rendering.

## Finish d02d261c after every run completes

Do not analyze a partial run or mix source revisions.
The required runs are `e2-rnacentral`, `e2-synthetic-shapes`, `e2-memory-breakdown`, `e2-parallel-impact-w1` and `e2-parallel-impact-w10`, all suffixed `-d02d261c-20261017`.
Check that every planned task's `task.json` has status `success`, then refresh the root manifests with `aggregate_array.py`.
Root `run.json` can be stale until this aggregation; `.execution.lock` alone is not completion evidence.
The analyzer validates artifacts, input identities, repetitions, settings and software provenance without extraction to shared storage.
The command below produces the same aggregate schema as the febd2bc1 input:

```bash
for cohort in rnacentral synthetic-shapes memory-breakdown parallel-impact-w1 parallel-impact-w10
do
  PYTHONDONTWRITEBYTECODE=1 python3 experiments/scripts/aggregate_array.py \
    --plan "$PWD/experiments/results/runs/e2-${cohort}-d02d261c-20261017/array_plan.json"
done
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR="$TMPDIR/e2-figure-mpl" python3 experiments/scripts/analyze_e0_e2_runs.py \
  --synthetic-run experiments/results/runs/e2-synthetic-shapes-d02d261c-20261017 \
  --rna-run experiments/results/runs/e2-rnacentral-d02d261c-20261017 \
  --previous-synthetic-run experiments/results/runs/e2-synthetic-shapes-febd2bc1-20261016 \
  --previous-rna-run experiments/results/runs/e2-rnacentral-febd2bc1-20261016 \
  --memory-run experiments/results/runs/e2-memory-breakdown-d02d261c-20261017 \
  --previous-memory-run experiments/results/runs/e2-memory-breakdown-febd2bc1-20261016 \
  --parallel-w1-run experiments/results/runs/e2-parallel-impact-w1-d02d261c-20261017 \
  --parallel-w10-run experiments/results/runs/e2-parallel-impact-w10-d02d261c-20261017 \
  --workers 10 --revision d02d261c --previous-workers 10 --previous-revision febd2bc1 \
  --previous-resource cpu_40=1 --skip-cap30-history \
  --out experiments/results/local/e2-d02d261c-20261017
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR="$TMPDIR/e2-figure-mpl" python3 experiments/scripts/make_figures.py --review-layout --figure3-only --e2-version d02d261c --out experiments/results/figures/main_d02d261c_20261017
```

Remove the temporary Matplotlib cache after analysis.
No RNA program, scheduler submission or implementation checkout is needed for aggregation or plotting.
