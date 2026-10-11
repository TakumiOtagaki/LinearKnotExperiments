# Experiment preparation progress

## 2026-10-11

### Runner and verification

- Added optional `order.kind: blocked_random` with an integer seed. Each
  repetition shuffles all input × mode × program combinations; warm-ups precede
  each combination's first measurement. The default nested order is unchanged.
  `run.json` keeps the planned schedule and actual attempts, including starts
  that do not finish.
- Added `--resume`. Config bytes/snapshot, software lock, source/binary state,
  build settings/cache, energy parameters and source/copied input hashes must
  match before the existing run is changed. Every `result.json` is retained,
  including failed jobs and results not yet appended to the manifest. Missing
  results restart with old logs preserved under `incomplete_attempts/`.
  Resume times, environment, node and previous checks are recorded. Legacy
  fixed-order manifests are supported.
- Fixed interrupted-run summaries: missing finalization fields give
  `run_status: interrupted` and null end-of-run change checks.
- `python3 -m unittest discover -s tests -v`: all 18 tests passed, including
  seven new runner/order/resume/summary checks and one scaffold subset check.
  `python3 -m pytest tests` was unavailable (`No module named pytest`); no
  package was installed. `git diff --check` passed.
- Summarized the original `e2-rnacentral-20261008` directory into ignored
  `results/local/e2-rnacentral-20261008-interrupted-prep-20261011.{tsv,json}`.
  It correctly reports four completed measurements (one success, two failures,
  one timeout) and an interrupted run. The original manifest stayed unchanged
  (SHA-256 `2d1617f01fc60fce7d684a5f8d10a671b06c992c7c477e1cbf963062ee56ae7e`).

### Inputs and the next feasibility measurement

- Extended `prepare_biological_inputs.py` with `--timeout-seconds` and `--ids`,
  preserving its default 120-second budget and original scaffold construction.
- Reused the saved G for URS00023C9544 (11,832 nt, original 64.98 s) and
  URS000241043C (13,652 nt, original 81.82 s), with unchanged input/G hashes.
- Regenerated URS000238C72C (22,543 nt) successfully in 137.592357 s, retaining
  4,465 pairs, and URS00023A4DB9 (29,106 nt) successfully in 182.176338 s,
  retaining 6,539 pairs. They ran sequentially with thread settings at 1 on
  node r8n10, scheduler job 8970385, using 3,600 seconds and the same 8 GiB
  address-space cap. LinearPartition-V b450fb3e, wrapper/binary hashes, beam
  200, dangles 2, posterior threshold 0.7 and minimum helix length 2 matched the
  original manifest and stayed unchanged during preparation.
- Saved the two new inputs and dataset under `data/rnacentral-retry-20261011/`,
  retaining the initial timeout records. Exact commands/results and raw-log
  hashes are also in `manifests/rnacentral-scaffold-retry-20261011.json`.
- Added `configs/e2_rnacentral_long_NEXT_20261011.json`: seven original records,
  three CLIs, beam 100 and local limits 30/30, matching output flags,
  3,600 seconds/job, 64 GiB address space, warm-ups 0, repetitions 1 and seed
  20261011. Removed the former 8,191 nt exclusion in this new matrix only.
  Old configs and measurements remain historical records.
- Updated README, protocol and the biological pilot notes, including the
  requirement to avoid concurrent calculations on the measurement node.

### Scientific scope and remaining work

The new matrix is prepared; its conditional CLI jobs have not been run.
Scaffold preparation success establishes saved G availability, not conditional
CLI feasibility, approximation accuracy or timing rankings. The author must
replace `manifests/software-pilot-NEXT.json` with a verified lock for the next
fixed commit/build and set `model_settings.build_revision_verified` to true
before running the matrix. The selected software must support all seven lengths.

All edits are inside the experiment repository. The external LinearCParty
working tree was not changed or built. No git commit or push was performed.

## 2026-10-11 TSUBAME array-job preparation

### Changed files and behavior

- Added `scripts/plan_array.py`, `scripts/run_array_task.py` and
  `scripts/aggregate_array.py`. Plans freeze input/config/lock copies, assign
  inclusive length tiers, pack complete sequence matrices by time/memory
  estimates, and generate one Grid Engine array script per populated tier.
  Scaffold variants of a shared sequence remain together. Planning never
  submits jobs and refuses existing output/script directories.
- Refactored `scripts/run_cli.py` to expose its existing provenance/resume path
  to array tasks and support bounded process concurrency for accuracy arrays.
  Seeded dispatch order, warm-up dependencies and repetition barriers are
  retained. Timing arrays require `node_f=1`, one worker and single-thread CLIs.
  Formal configs require verified build flags in both config and software lock;
  direct CLI execution remains development-only. SIGTERM is handled by the task
  wrapper so available interruption records are finalized.
- Updated `scripts/summarize.py` and `scripts/validate_outputs.py` for the merged
  run, including task counts and planned/completed/missing measurement counts.
  Aggregation reads result markers, recovers missed manifest updates, preserves
  failures and emits explicit rows for missing jobs. Task manifests retain
  scheduler IDs, host, resources, timestamps and attempts. Kernel-released task
  locks reject duplicate controllers. Shared manifests are updated only by
  explicit aggregation.
- Added `tests/test_array.py`, `configs/arrays/e2_rnacentral_20261011.json` and
  `configs/arrays/e1_bprna_20261011.json`. Updated `.gitignore` for generated job
  scripts. Updated README and `docs/experiment_protocol.md`; added
  `docs/array_dry_run_20261011.md` and exact plans/scripts/summaries under
  `docs/array_examples/20261011/{e1,e2}/`. Original experiment configs are unchanged.

### Verification and dry-run artifacts

- `python3 -m unittest discover -s tests`: 27 tests passed, including all 18
  existing tests. The new tests use small dummy CLIs and cover packing and
  estimates, scaffold grouping, resource/wall/memory guards, SGE environment
  selection, immutable preparation, matching and rejected resumes, task locks,
  worker bounds, warm-up/block barriers, SIGTERM, startup failure and partial
  aggregation/summary/output validation. No biological CLI matrix was run.
- Python compilation and `git diff --check` passed. Both saved plan checksums,
  all five saved script checksums and `bash -n` for those scripts passed.
- Prepared `results/runs/array-example-e2-20261011` from the existing seven-input
  `e2_rnacentral_15aa7742_20261011.json`, using the read-only dedicated
  `deps/LinearCParty-15aa7742` and `build/linear-cparty-15aa7742`. The software lock
  verified. There are three arrays and five tasks: three 4.8–7.0 kb inputs share
  a task; 11.8 kb and 13.7 kb inputs each have a task; 22.5 kb and 29.1 kb inputs
  each have a task. There are 21 planned CLI jobs, all not started.
- Prepared `results/runs/array-example-e1-20261011` from the existing 12-input
  `e1_bprna_pilot_20261008.json`, using read-only `deps/LinearCParty-20261008`
  and `build/linear-cparty-20261008`. The software lock verified. Two arrays
  contain one task each (8 inputs below 200 nt and 4 inputs at 200–800 nt), with
  two accuracy workers and 72 planned CLI jobs, all not started.
- Aggregated and summarized both unstarted matrices: 5/2 not-started tasks,
  21/72 missing measured jobs and zero completed measurements. Generated scripts
  live under `jobs/array-example-e{1,2}-20261011/`; review copies are in docs.

### Scientific scope and remaining checks

The example cost coefficients are explicitly uncalibrated scheduling assumptions.
Before formal collection, calibrate them on the adopted fixed-build pilot and
confirm longest-input CLI feasibility, memory and an appropriate failure budget.
Historical input/config purposes remain development; these dry-runs establish
preparation and aggregation behavior, not biological accuracy or runtime scaling.
If an interrupted timing task resumes on another node, examine its attempt hosts
and CLI resume history before using paired comparisons.

The resource table was checked against the official TSUBAME4 handbook, section
5.1.1, linked in README. Live `qconf -sc`, `-spl` and `-sql` checks failed with
`getgrgid(3215) failed: gid resolve error`. Current accounting group, scheduler
limits, concurrency and runtime modules are marked **要確認** in README.

All changes remain inside `experiments/`. No qsub, biological measurement,
external source modification/build, manuscript/textbook edit, git commit or push
was performed.
