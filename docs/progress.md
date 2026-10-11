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
