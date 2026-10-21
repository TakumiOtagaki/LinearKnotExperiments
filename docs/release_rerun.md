# Release preparation and reruns

The [4f97fd3 preparation](linearknot_4f97fd3_20261021/README.md) records strict
378/378 output agreement with 9c7d6e1 and unsubmitted resource-only reruns.
Their configs use `executable` and `subcommand` to call the unified CLI while
preserving distinct analysis names in the packed runner.

For the renamed unified CLI, see [LinearKnot 9c7d6e1 preparation and comparison](linearknot_9c7d6e1_20261020/README.md).
`make_software_lock.py` accepts the single `linearknot` binary and records its subcommands; the historical `prepare_release.py` naming and production runners below still target the separate CLIs.
The [centroid replacement plans](centroid_89baea52_20261020/README.md) use the adopted 89baea52 build and include local commands and unsubmitted array scripts.

Run from the experiment repository with `PYTHONDONTWRITEBYTECODE=1` and node-local `TMPDIR` or `--scratch-dir /tmp`.
These tools never submit scheduler jobs.

`prepare_release.py COMMIT` performs a local clone, HEAD/clean checks, CMake Release build, required suite, software lock and output comparison in that order.
It uses `deps/LinearCParty-<first-eight-hex>` and `build/linear-cparty-<first-eight-hex>`.
Existing stages are inspected without rebuilding; an existing build needs its matching certified lock or preparation receipt.
The suite receives `SCFG2_BUILD_DIR` and a temporary `source/build` symlink, which the wrapper removes on success or failure.
Logs and test output stay in scratch; one JSON receipt retains results, hashes and bounded log tails.
New clones/builds pass an inode check using a tracked-file allowance plus `--estimated-build-files` (default 20,000).

```bash
python3 scripts/prepare_release.py e2cecbb4 \
  --clone-from deps/LinearCParty-e2cecbb4 \
  --input data/smoke/pk_scaffold.txt --input data/smoke/empty_scaffold.txt \
  --scratch-dir /tmp --out manifests/release-preparation-e2cecbb4.json \
  --dry-run --estimate-only
```

Remove the final two flags and choose a fresh receipt path to perform preparation and comparison.
The donor clone is read only, and dry-run never creates clones, builds or locks or runs tests/RNA CLIs.
Repeated real preparation reuses a successful comparison only while its inputs, parameters, arguments and binary fingerprints match.

`compare_release_outputs.py` compares stdout, stderr and BPP exports, removing only lines starting with `# elapsed_seconds`.
It runs LinearCParty with `--bpp-cutoff 0 --bpp`, LinearPKCapR, and LinearPKRaccess with `--access-len 1,4,8,16 --probability`.
Default beams are 20/200 with `h=d=max(b,30)` and `--ranking outer`; exact uses `--exact` without local caps.
Beam and exact references default to d21618c5 and d9804df6 and can be changed independently with `--beam-reference-build` and `--exact-reference-build`.
The JSON records byte agreement, aligned numeric maxima, differing-line examples, commands, hashes, individual RSS and wall times; failed or timed-out processes cannot count as agreement.
Numeric maxima cover lines with matching nonnumeric text, so inspect `all_differences_numeric` and examples for unmatched text/rows.
Output differences inform the author's rerun decision; the tool imposes no scientific acceptance threshold.

```bash
python3 scripts/compare_release_outputs.py \
  --candidate-build build/linear-cparty-d21618c5 \
  --candidate-exact-build build/linear-cparty-d9804df6 \
  --parameter-file deps/LinearCParty-d21618c5/params/rna_Turner04.par \
  --input data/smoke/pk_scaffold.txt --input data/smoke/empty_scaffold.txt \
  --scratch-dir /tmp --out /tmp/release-self-comparison.json
```

The separate exact candidate is for the rehearsal; omit it for a single new release build.
All input copies and per-process outputs are removed from scratch after comparison.

`make_rerun_configs.py` takes the small [release/resource settings](../configs/rerun_e2cecbb4_20261016.json) and preserves the five production templates' modes, inputs, programs, metrics, warm-ups and repetitions.
It adds paired workers=1/8 configs with the same 4,793/5,482/6,990 nt RNAcentral records and all nine input/CLI combinations on one task.
Warm-up dependencies and repetition barriers can use fewer active slots than the configured maximum.
The [generated configs](../configs/rerun_e2cecbb4_20261016/rerun_manifest.json) use cpu_16, single-thread processes, explicit timing-concurrency reasons and a memory envelope of 8 × 4,096 + 1,024 = 33,792 MiB within 35,095 MiB.
Ten workers need a smaller cap, such as 3,072 MiB; oversubscription is rejected before writing configs.
Generated scripts use `-tc 30`; coordinate multiple arrays so their combined weekday concurrency also stays within 30.

```bash
python3 scripts/make_rerun_configs.py \
  --settings configs/rerun_e2cecbb4_20261016.json \
  --out-dir /tmp/e2cecbb4-configs --dry-run --estimate-only
```

Without `--dry-run`, the generator writes configs and a manifest only.
Offline estimates require the referenced locks; pass `--software-root` and `--build-dir` to `plan_array.py` when verifying source/build hashes.
The memory-breakdown config needs its own same-release stage driver and lock, prepared separately with the existing builder.

```bash
python3 scripts/build_stage_timing_driver.py \
  --software-root deps/LinearCParty-e2cecbb4 \
  --build-dir build/linear-cparty-e2cecbb4 \
  --base-lock manifests/software-pilot-e2cecbb4.json \
  --out build/stage-timing-e2cecbb4 \
  --lock manifests/software-stage-timing-e2cecbb4-20261016.json --purpose formal
```

Time coefficients come from the templates, multiplied by `parallel_slowdown_factor` (1.5 here), with no assumed parallel speedup.
The example uses 24-hour h_rt, a 70% packing target and at most 1,024 inputs per task; the paired check uses 12 hours.
Original memory fits remain in provenance; the explicit RLIMIT_AS cap becomes the planning envelope, whose feasibility still needs calibration.

The [rehearsal report](release_rerun_rehearsal_20261016.json) records d21618c5 beam/d9804df6 exact, 18/18 identical comparisons over two small inputs, and all offline estimates.
The five production configs need 441 tasks and the paired checks add two, totaling 443 tasks, 397,322 CLI jobs including warm-ups, and 5,899 planned shared entries including directories.
Config generation adds 17 permanent entries; estimate-only generation adds 14 scratch entries.
Sequential-work envelopes total 5,782.36 hours, the largest task is 16.79 hours, and h_rt requests total 10,608 task-hours.
E2 plus the paired check alone needs 31 tasks, 2,063 CLI jobs and 499 planned shared entries.
The actual e2cecbb4 output comparison, stage-driver build/lock and memory/time calibration remain outstanding; no new build or scheduler submission occurred in this rehearsal.
If the author retains E0/E1 after output comparison, submit only the selected E2 reruns.

## Completed E2 e2cecbb4 rerun analysis (2026-10-16)

The actual reruns use cpu_40, at most 10 single-thread processes per allocation,
and an 8 GiB RLIMIT_AS per process, superseding the cpu_16 rehearsal above.
RNAcentral and synthetic runs retain one warm-up and three measurements;
the memory driver uses one fresh process per input/stage without warm-ups.
The w1/w10 check has nine input/CLI combinations, so w10 reaches nine active
processes. Hosts differ between w1, w10 and production; report paired ratios
without attributing the entire difference to concurrency or applying a correction.

`analyze_e0_e2_runs.py` keeps its historical run defaults. Omit `--e0-analysis`
for E2-only work. Use `--workers 10 --previous-resource cpu_40=1
--skip-cap30-history`, `--previous-synthetic-run`, `--memory-run` with
`--previous-memory-run`, and `--parallel-w1-run` with `--parallel-w10-run`
to compare the completed reruns. Specify the new RNAcentral and synthetic
runs with the existing `--rna-run` and `--synthetic-run` options and d21618c5
with `--previous-rna-run`. All paths are relative to the caller for analysis;
`aggregate_array.py --plan` resolves relative paths from the experiment root,
so use an absolute plan path when calling from the paper root.

The output `results/local/e2-e2cecbb4-20261016/` contains 14 files: CLI metrics,
comparisons, slopes and model/provenance JSON; four CLI figures; stage metrics,
family counters, stage comparisons and provenance JSON; one stage figure; and
validation JSON. Packed artifacts are streamed and never extracted.
The public work counters are cumulative main-forward counts, not simultaneous
live retention or allocation bytes. Profiles do not expose work; these fields
remain NA. Stage RSS differences are not additive memory components.

`validate_e0_e2_runs.py OUTPUT --test-log LOG` accepts an E2-only output;
the historical `--e0-output` option still validates E0 as before.
Exclude `test_stage_timing.py` when collecting tests without RNA CLI execution.

## E2 febd2bc1 replacement (2026-10-16)

The five febd2bc1 runs use the same cpu_40, workers=10, RLIMIT_AS=8 GiB,
beam/local caps and requested outputs as e2cecbb4 (the paired w1 run uses
workers=1). Refresh their root manifests with `aggregate_array.py --plan`
before analysis; task JSONL and raw tar files are unchanged.
The release receipt is `manifests/release-preparation-febd2bc1.json`.

Analysis accepts `--revision` and `--previous-revision` as full commits or
hexadecimal prefixes. These pins check CLI and memory-driver task revisions;
current parallel checks must also match the current production signature.
`--previous-workers` defaults to 1 for historical reproduction; specify 10
when comparing with e2cecbb4. Existing run defaults and E0 options remain
unchanged. Omitting revision pins retains historical behavior.

From the paper root:

```bash
PYTHONDONTWRITEBYTECODE=1 TMPDIR=/tmp \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
MPLCONFIGDIR=/tmp/e2-febd2bc1-mpl \
python3 experiments/scripts/analyze_e0_e2_runs.py \
  --revision febd2bc1 --previous-revision e2cecbb4 \
  --workers 10 --previous-workers 10 --previous-resource cpu_40=1 \
  --skip-cap30-history \
  --synthetic-run experiments/results/runs/e2-synthetic-shapes-febd2bc1-20261016 \
  --rna-run experiments/results/runs/e2-rnacentral-febd2bc1-20261016 \
  --previous-rna-run experiments/results/runs/e2-rnacentral-e2cecbb4-20261016 \
  --previous-synthetic-run experiments/results/runs/e2-synthetic-shapes-e2cecbb4-20261016 \
  --memory-run experiments/results/runs/e2-memory-breakdown-febd2bc1-20261016 \
  --previous-memory-run experiments/results/runs/e2-memory-breakdown-e2cecbb4-20261016 \
  --parallel-w1-run experiments/results/runs/e2-parallel-impact-w1-febd2bc1-20261016 \
  --parallel-w10-run experiments/results/runs/e2-parallel-impact-w10-febd2bc1-20261016 \
  --out experiments/results/local/e2-febd2bc1-20261016
python3 experiments/scripts/validate_e0_e2_runs.py \
  experiments/results/local/e2-febd2bc1-20261016 \
  --revision febd2bc1 --previous-revision e2cecbb4 --test-log /tmp/e2_febd2bc1/tests.log
PYTHONDONTWRITEBYTECODE=1 python3 experiments/scripts/make_figures.py \
  --figure3-only --figure3-panels abc --e2-revision febd2bc1 \
  --e2 experiments/results/local/e2-febd2bc1-20261016 \
  --out experiments/results/figures/draft_febd2bc1_20261016
```

Use a fresh analysis output directory for reproduction.
The validator defaults to the pins saved by the analysis and checks them
against task metadata as well as the exported provenance.
The figure builder checks both CLI and stage revisions against `--e2-revision`;
its historical default is e2cecbb4. `--figure3-panels abc` separates stage RSS
into `figure3_mode_memory_{183,89}mm.pdf/png` for the Supplementary material.
The default `abcd` layout remains available for historical figures.
No output in `results/figures/draft_20261016/` is overwritten by this command.
