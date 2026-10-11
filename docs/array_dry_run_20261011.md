# TSUBAME array preparation: 2026-10-11 dry-runs

These examples prepare the two existing development configs. No CLI measurement
or scheduler submission was performed. The original configs and the external
LinearCParty working tree were unchanged. Both pinned source/build pairs passed
the existing clean-source, build-cache, binary and parameter software-lock checks.
All cost coefficients below are uncalibrated examples, not measured scaling.

Run from the experiment repository, choosing new output/job directory names
when repeating preparation:

```bash
python3 scripts/plan_array.py \
  --config configs/e2_rnacentral_15aa7742_20261011.json \
  --array-config configs/arrays/e2_rnacentral_20261011.json \
  --software-root deps/LinearCParty-15aa7742 \
  --build-dir build/linear-cparty-15aa7742 \
  --out results/runs/array-example-e2-20261011 --dry-run

python3 scripts/plan_array.py \
  --config configs/e1_bprna_pilot_20261008.json \
  --array-config configs/arrays/e1_bprna_20261011.json \
  --software-root deps/LinearCParty-20261008 \
  --build-dir build/linear-cparty-20261008 \
  --out results/runs/array-example-e1-20261011 --dry-run
```

The local prepared directories contain input/config/lock copies, task configs,
task manifests, a top-level manifest and the checksummed plan. Generated job
scripts are under `jobs/array-example-e{1,2}-20261011/`. These local artifacts
are ignored by Git. Exact copies of both plans, checksums, five scripts and
unstarted summaries are saved below for review. Their absolute paths refer to
the preparation workspace; regenerate for another installation.

## RNAcentral: 7 sequences, 21 CLI jobs, 3 arrays, 5 tasks

One measurement per program, three programs, beam 100, no warm-up,
3,600-second timeout and 65,536 MiB address-space cap per CLI job are unchanged.
All tiers request `node_f=1` and use one CLI worker. The tier memory budget is
73,728 MiB; the conservative cap envelope is 66,560 MiB including task overhead.
Every task uses blocked random order with seed 20261011.

| Tier / SGE_TASK_ID | Sequence IDs | Lengths (nt) | CLI jobs | Estimated seconds | Estimated peak MiB | h_rt |
|---|---|---|---:|---:|---:|---|
| rna_4_8kb / 1 | URS000233CDDE, URS00019626CE, URS000239D8A4 | 4793, 5482, 6990 | 9 | 2330.775 | 12110 | 03:00:00 |
| rna_8_20kb / 1 | URS00023C9544 | 11832 | 3 | 1597.320 | 21048 | 02:00:00 |
| rna_8_20kb / 2 | URS000241043C | 13652 | 3 | 1843.020 | 22868 | 02:00:00 |
| rna_20_30kb / 1 | URS000238C72C | 22543 | 3 | 7101.045 | 43030.5 | 05:00:00 |
| rna_20_30kb / 2 | URS00023A4DB9 | 29106 | 3 | 9168.390 | 52875 | 05:00:00 |

The medium tier's 3,000-second target separates its two inputs. The long tier's
10,000-second target keeps each 20–30 kb sequence alone. Short inputs share one
task under the 3,600-second target. Counts are consequences of estimated budgets,
not fixed sequences-per-task settings.

- [Complete plan](array_examples/20261011/e2/array_plan.json) and [checksum](array_examples/20261011/e2/array_plan.sha256).
- Generated [4–8 kb script](array_examples/20261011/e2/rna_4_8kb.sh), [8–20 kb script](array_examples/20261011/e2/rna_8_20kb.sh), [20–30 kb script](array_examples/20261011/e2/rna_20_30kb.sh).
- Unstarted [summary JSON](array_examples/20261011/e2/summary.json) and [TSV](array_examples/20261011/e2/summary.tsv).

## bpRNA: 12 sequences, 72 CLI jobs, 2 arrays, 2 tasks

Exact and beam 100, three programs, one measurement, no warm-up,
180-second timeout and 8,192 MiB per-job address-space cap are unchanged.
The array overlay supplies seed 20261011. Both tiers use `cpu_16=1`, two CLI
workers and a 24,576 MiB memory budget. The cap envelope is 16,896 MiB.
The time estimate still sums all work rather than assuming a twofold speedup.
These concurrent elapsed times would be accuracy-run diagnostics only.

| Tier / SGE_TASK_ID | Inputs | Length range (nt) | CLI jobs | Estimated seconds | Estimated peak MiB | h_rt |
|---|---:|---|---:|---:|---:|---|
| bprna_under200 / 1 | 8 | 52–139 | 48 | 14.094 | 2648 | 01:00:00 |
| bprna_200_800 / 1 | 4 | 201–367 | 24 | 91.152 | 8432 | 02:00:00 |

- [Complete plan](array_examples/20261011/e1/array_plan.json) and [checksum](array_examples/20261011/e1/array_plan.sha256).
- Generated [under-200 nt script](array_examples/20261011/e1/bprna_under200.sh) and [200–800 nt script](array_examples/20261011/e1/bprna_200_800.sh).
- Unstarted [summary JSON](array_examples/20261011/e1/summary.json) and [TSV](array_examples/20261011/e1/summary.tsv).

## Aggregation and remaining checks

Both prepared runs were aggregated and summarized without starting tasks.
RNAcentral reports five not-started tasks, 21 planned and 21 missing measured
jobs; bpRNA reports two not-started tasks, 72 planned and 72 missing measured
jobs. Each summary has zero completed measurements, null end-of-run dependency
checks and run status `interrupted` for the incomplete matrix. Unstarted jobs
are explicit rows, not omitted. `validate_outputs.py` and the accuracy reader
can use the merged artifact paths; the small integration tests check successful
and missing outputs, failures and resume preservation.

The [README resource table](../README.md#tsubame4-resources-and-remaining-checks)
records the official source and the failed live `qconf` check. Before submission,
the accounting group, current scheduler limits, concurrency and required runtime
modules are **要確認**. Cost coefficients need pilot calibration on the adopted
build. The 20–30 kb CLI feasibility and RSS estimates are **要確認**. This work
creates a formal-experiment mechanism; these historical configs retain their
development purpose and do not establish scientific results.
