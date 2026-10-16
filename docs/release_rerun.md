# Release preparation and reruns

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
