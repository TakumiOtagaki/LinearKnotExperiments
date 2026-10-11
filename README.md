# LinearCPartyExperiments

Experiments for the LinearCParty paper: exact-reference validation, beam-width
accuracy, runtime and memory scaling, RNA examples and scaffold sensitivity.
The experiment repository can run independently of the manuscript repository.

The initial runner uses Python's standard library and GNU `/usr/bin/time` on
Linux. It runs the existing LinearCParty CLI binaries without modifying or
building their source tree.

## Accuracy pilot preparation

`configs/e1_pilot.json` defines 18 saved synthetic inputs and 324 development
jobs at exact and finite beam settings. Inputs, seeds, hashes and scaffold
geometry are in `data/pilot-v1/`. The matrix is prepared but has not been run.
See [accuracy analysis](docs/accuracy_analysis.md) for commands, numerical
definitions, cap controls and limits of the measurement.

`scripts/analyze_accuracy.py` compares BPP, seven structural contexts and
joint window accessibility. It keeps failed jobs and missing references,
uses full-precision Raccess thermodynamic metadata, and checks consistency
across the three programs. It does not yet score decoder structures or
sampling frequencies.

## Biological input candidates

The next accuracy pilot will use short bpRNA sequences. The author approved
using the existing LinearRaccess-paper collection, containing 26,063 distinct
ACGU sequences. A byte-identical copy is now stored locally, with its source and
hash in [the acquisition manifest](manifests/data-sources/bprna-local-20261008.json).
This collection is different from the LinearCapR paper's evaluated subset.
Read the [bpRNA data audit](docs/data_audit_bprna_20261008.md) for the filtering
history and the original-distribution checks that remain unresolved.

For initial long-input timing and RSS measurements, seven ACGU-only sequences
of 4,793–29,106 nt are available in the published LinearCapR RNAcentral panel.
The complete 20-record panel has been downloaded and its expected SHA-256
verified. The [acquisition manifest](manifests/data-sources/rnacentral-panel-20261008.json)
records the fixed commit and the seven eligible IDs and hashes. Read the
[RNAcentral data audit](docs/data_audit_rnacentral_20261008.md).
The bpRNA pilot and three shorter RNAcentral records have been measured on
separate pinned development versions. See [the biological pilot](docs/biological_pilot_20261008.md)
for their scope and the next seven-record feasibility matrix.

## Formal experiments on TSUBAME arrays

The [full bpRNA preparation](docs/full_bprna_array_plan_20261011.md) covers
26,063 scaffold-generation jobs, the 699,387-job E1 matrix, staged b300/b400
pilots, and a 304-input CParty comparison. Configs, measured calibration,
software locks, tier assignments and script copies are saved there.
No scheduler submission was made. Finish the G arrays and regenerate verified
E1/E0 plans before execution; their current dry-runs explicitly reject pending G.

`scripts/scaffold_array.py plan/run/aggregate` generates and collects the G
arrays using the existing posterior/helix procedure. `scripts/build_e0_tools.py`
builds inside-only exact observers for both closure policies and the single
four-patch CParty oracle (upstream `800e95a`, without a closure patch).
`scripts/analyze_e0.py` compares log Z / PF free energy only, with absolute
roundoff limits of 1e-9 and 1e-10 kcal/mol from the implementation's audit.
It counts closure differences after compatible-versus-oracle parity passes.
BPP and MFE equality are outside E0; the historical BE41 MFE difference is
preserved in the earlier records. All four PF witnesses pass. The adopted
LinearCParty remains `15aa7742`; 28 small comparisons with `2456f9b2` have
identical PF values. See [the update evidence](docs/e0_oracle_update_20261011/).

Run formal experiments as tier-specific TSUBAME Grid Engine array jobs.
Interactive-node execution is limited to small development pilots.
`scripts/plan_array.py` freezes inputs, checks a software lock, saves a plan
and generates `jobs/<run>/<tier>.sh`. It never calls `qsub`.
`scripts/run_array_task.py` selects a tier's task using `SGE_TASK_ID`.
`scripts/aggregate_array.py` merges all tasks into the run's `run.json`.

One task contains a batch of sequences and **all** their programs, exact/beam
settings, beam widths, warm-ups and repetitions. Supplied scaffold variants
with the same sequence also stay together. Task execution uses the existing
`blocked_random` repetition order with the array config's seed. Comparisons for
a sequence therefore share the allocated node within that attempt.
Timing arrays require `node_f=1`, one CLI process at a time and single-thread
settings. Accuracy arrays may choose `workers` per tier; their concurrent
elapsed times must not be used as isolated runtime measurements.

Keep the original matrix config and put scheduler settings in its `array`
object or a separate `--array-config`. Examples are
[`configs/arrays/e2_rnacentral_20261011.json`](configs/arrays/e2_rnacentral_20261011.json)
and [`configs/arrays/e1_bprna_20261011.json`](configs/arrays/e1_bprna_20261011.json).
Each tier has inclusive `min_length`/`max_length`, `resource`, `h_rt`,
`target_seconds`, `memory_budget_mib`, `task_overhead_mib`, `workers` and
`cost_model`. Intervals must be sorted and nonoverlapping; every input must fit.
Packing uses estimated work rather than a fixed sequence count. A complete
sequence matrix that cannot fit is refused rather than split across nodes.

For length `n`, per-condition time is
`seconds_constant_per_condition + seconds_per_nt_condition * n**length_exponent`.
Multiply by programs, warm-ups plus repetitions, the sum of `mode_factors`
(default 1 per mode), and `safety_factor` (at least 1). Sum all records' work
in a task; accuracy planning does not assume ideal parallel speedup.
Estimated memory is `task_overhead_mib + workers * max(memory_base_mib +
mib_per_nt * n)`. Both this estimate and the stricter envelope
`task_overhead_mib + workers * memory_limit_mib` must fit the tier's memory
budget. The original job timeout and address-space cap stay unchanged by default.
An explicit tier `job_limits` overlay may set `timeout_seconds` and
`memory_limit_mib`; these values are validated and saved in each task config.
Modes may have an inclusive `length_range`, used in packing and execution
so exact is absent at n≥400 in full bpRNA E1. A `cost_model.per_mode` map can
give distinct constant/coefficient/exponent time and memory envelopes per mode.
See the full bpRNA configs and calibration manifest for the adopted formulas.
`target_seconds` must fit `max_wall_fraction` of `h_rt` (default 0.7, maximum 0.8).
Tasks save their estimates, coefficients and the stated `basis` of the model.
The example coefficients are uncalibrated scheduling estimates, not scientific
complexity claims. Refit them with the adopted build's pilot, including failures.
Many individual timeouts can exceed an estimated task wall budget.

Full bpRNA uses `compact_plan: true`: each task's complete schedule is a
checksummed `planned_jobs.json`, and the initial aggregate stores counts only.
Use `aggregate_array.py --summary-only` for status summaries; omit the flag
when observations are required. Accuracy analysis rejects a counts-only
aggregate. Analyze large matrices per task and combine RNA-level observations.
A `manifest_checkpoint_jobs` setting reduces repeated full manifest writes;
an append journal records every start/finish and recovers execution order on
resume, while each job's `result.json` remains the marker.
Long tiers support `--through-beam 20`, then `--resume --through-beam 50`,
100 and 200 as separate invocations after reviewing each stage.
Larger widths remain missing until explicitly resumed.

Formal configs require `purpose: formal`, a real `software_lock`, and
`model_settings.build_revision_verified: true`. The lock must also record a
verified build revision. Supply the clean pinned source and dedicated build:

```bash
python3 scripts/plan_array.py --config configs/formal.json \
  --array-config configs/arrays/formal.json \
  --software-root /path/to/pinned/LinearCParty \
  --build-dir /path/to/dedicated/build --out results/runs/formal-001
```

The planner saves byte copies of the original config and lock, input copies
and hashes, task configs, full seeded schedules, and `array_plan.json` with a
checksum. Existing output or script directories are refused. Task execution
checks frozen files and the live lock again. An offline `--dry-run` may omit
software/build paths: it checks lock structure without claiming to verify live
dependencies. At execution, supply `LCP_SOFTWARE_ROOT` and `LCP_BUILD_DIR`, or
the wrapper's path options. A dry-run with both paths supplied also checks
the actual dependencies.

Submit each tier separately from this repository with the authorized accounting
group and concurrency limit. The script already contains its range, resource
and wall time. These commands are examples for later submission; they were
**not executed** during preparation:

```bash
qsub -g YOUR_TSUBAME_GROUP -tc 2 jobs/formal-001/short.sh
# Resume only task 3 after its earlier processes have stopped:
qsub -g YOUR_TSUBAME_GROUP -t 3-3 jobs/formal-001/short.sh --resume
```

Do not override the saved resource or `h_rt`. Outputs live under
`tasks/<tier>/<six-digit-task-id>/run/`. The adjacent `task.json` records task
and global numbers, `JOB_ID`, `SGE_TASK_ID`, hostname, requested resource,
wall time, workers, start/end times and all attempts. A per-task `flock` rejects
concurrent execution; the kernel releases it even if the controller is killed.
Resume keeps every `result.json`, including failures/timeouts, and runs only
missing results. Partial jobs are archived by the existing runner.
A killed controller retains a null end time and is counted as interrupted.
Stop surviving processes before resuming. If a resume changes nodes, inspect
task attempts and CLI `resume_history` before pooling timing segments.

Refresh the aggregate after execution or resume, including interrupted runs:

```bash
python3 scripts/aggregate_array.py --plan results/runs/formal-001/array_plan.json
python3 scripts/summarize.py results/runs/formal-001 --out results/local/formal-001.tsv
python3 scripts/validate_outputs.py results/runs/formal-001 \
  --out results/local/formal-001-validation.json
```

The aggregate recovers result markers even after a missed manifest update,
keeps missing jobs with null measurements, and counts not-started, running,
interrupted and failed tasks. Summary JSON separates planned, completed and
missing measured jobs. Validation reports missing observations separately and
refuses an actively running aggregate. Tasks do not update the shared manifest;
re-run aggregation when their state changes.
See [the two dry-run plans and scripts](docs/array_dry_run_20261011.md).

### TSUBAME4 resources and remaining checks

The [official handbook, section 5.1.1](https://helpdesk.t3.gsic.titech.ac.jp/t4docs_preview/master/handbook.ja/jobs/)
lists these resources per allocation (checked 2026-10-11). The planner
conservatively converts the documented GB as decimal bytes to MiB.

| Resource | Physical CPU cores | Memory (GB, as documented) |
|---|---:|---:|
| `node_f` | 192 | 768 |
| `node_h` | 96 | 384 |
| `node_q` | 48 | 192 |
| `node_o` | 24 | 96 |
| `gpu_1` | 8 | 96 |
| `gpu_h` | 4 | 48 |
| `cpu_160` | 160 | 368 |
| `cpu_80` | 80 | 184 |
| `cpu_40` | 40 | 92 |
| `cpu_16` | 16 | 36.8 |
| `cpu_8` | 8 | 18.4 |
| `cpu_4` | 4 | 9.2 |

Fractional node and CPU allocations share a physical node, so timing uses
`node_f=1`. `NSLOTS` is recorded, not used to choose workers.
Live `qconf -sc`, `-spl` and `-sql` failed here with
`getgrgid(3215) failed: gid resolve error`.
Current scheduler/group limits, the authorized accounting group, concurrency
and required runtime modules remain **要確認** before submission.
Example cost coefficients and longest-input feasibility are also **要確認**;
the dry-runs perform no CLI measurements.

## Development smoke

Supply the implementation and build locations explicitly while LinearCParty is
under active development:

```bash
python3 scripts/run_cli.py --config configs/smoke.json \
  --software-root /path/to/LinearCParty \
  --build-dir /path/to/LinearCParty/build \
  --out results/runs/smoke-001
python3 scripts/summarize.py results/runs/smoke-001 \
  --out results/local/smoke-001.tsv
```

Use a new output directory for each run. To continue an interrupted run, repeat
the same command with `--resume`. The runner requires the original config bytes,
saved config, software lock, source state, binaries, build settings, parameter
file and original/copied input hashes to match. It keeps every job with a
`result.json`, including failures and timeouts, and runs only missing results.
Unfinished job files are moved to `incomplete_attempts/` before a fresh attempt;
resume times and execution attempts are recorded in `run.json`.
Stop any surviving child process before resuming. A run directory still needs
`--resume` even when all its jobs are complete.
Config paths, inputs and output paths supplied on the command line
are interpreted relative to the experiment repository. The supplied software
and build paths are interpreted relative to the shell's working directory.
Once added, `deps/LinearCParty` will be the default software location; the
current checkout does not pin the changing implementation yet.

The smoke has two small synthetic inputs, two modes (exact and beam width 100)
and five programs: LinearCParty, LinearPKCapR, LinearPKRaccess,
LinearPKSampling and LinearPKMEA. It checks the execution and logging path.
It is not an exactness benchmark or a biological experiment.

Each run keeps a config snapshot, copied input and energy-parameter files,
software and experiment Git state, executable hashes, selected CMake build
settings, the command and output of every job, elapsed time, peak RSS in KiB,
and return status. The runner returns nonzero if any measured job fails but
continues the remaining jobs. The summary includes failed jobs in its counts.
An unfinished manifest can be summarized: its run status is `interrupted`, and
unavailable end-of-run software/binary checks are JSON null.

The default execution order remains input → mode → program, with warm-ups and
all repetitions for each combination. Set
`"order": {"kind": "blocked_random", "seed": 20261011}` to shuffle all
input × mode × program combinations within each repetition block. Warm-ups
occur just before each combination's first measurement. `planned_jobs` saves
the full schedule; `execution_order` saves actual starts, including retries.

The next long-input config is
[`configs/e2_rnacentral_long_NEXT_20261011.json`](configs/e2_rnacentral_long_NEXT_20261011.json):
seven sequences (4,793–29,106 nt), three CLIs, default beam 100, 3,600 seconds/job,
64 GiB address space, no warm-up, one measurement and seed 20261011.
Its `manifests/software-pilot-NEXT.json` lock is an explicit placeholder.
Replace it after the author fixes the next source commit and verifies its
dedicated build, then update `model_settings.build_revision_verified` to true.
The template does not run until that real software lock exists.

An executable hash identifies the binary actually run. A source commit alone
does not establish that an existing build corresponds to that commit. The
initial smoke therefore records `build_revision_verified: false` and
`purpose: development`. Formal measurements need a separate pinned build.

## Repository contents

| Path | Contents |
|---|---|
| `configs/` | Versioned experiment settings; no personal software paths. |
| `data/smoke/` | Small synthetic sequence/scaffold records. |
| `data/pilot-v1/` | Frozen synthetic accuracy cohort and generation metadata. |
| `scripts/` | CLI runner and aggregation. |
| `docs/` | Execution protocol and experiment stages. |
| `manifests/` | Dependency-freeze guidance; formal locks will be added later. |
| `results/` | Reviewed summaries; local raw runs are ignored by Git. |
| `deps/` | Reserved for pinned dependencies. |

See [the protocol](docs/experiment_protocol.md) for interpretation and
[the experiment stages](docs/experiment_stages.md) for the next work.
The implementation is [LinearCParty](https://github.com/TakumiOtagaki/LinearCParty).
The manuscript is [LinearCPartyPaper](https://github.com/TakumiOtagaki/LinearCPartyPaper).

Run the runner's integration checks with:

```bash
python3 -m unittest discover -s tests -v
```

They exercise failures, process-group timeouts, warm-up exclusion, copied
inputs, fixed/seeded measurement order, matching and rejected resumes,
legacy manifests, interrupted summaries, scaffold subset preparation, array
packing, SGE task selection, task locks, parallel accuracy bounds and aggregation.
