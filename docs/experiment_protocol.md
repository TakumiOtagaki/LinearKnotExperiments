# Experiment protocol

Use the same sequence, supplied scaffold, thermodynamic parameter file,
temperature, local limits and requested output for a matched comparison.
Save the actual scaffold, rather than only the parameters used to create it.
For cap/beam controls and numerical domains, read
[accuracy analysis](accuracy_analysis.md). Comparisons against full exact
with capped finite beams measure the complete approximation, rather than
beam pruning alone. Keep uncapped-width observations separate.

## Development and formal measurements

The initial config is a development smoke. It permits a changing source tree
and an existing build, records their state and does not establish agreement
with upstream CParty. Its timings only verify measurement capture.

Before formal measurements, pin the implementation and reference commits,
build in dedicated directories, record the compiler and build flags, and
verify that the binaries come from those sources. Write the software lock in
`manifests/`. Source state and executable hashes are recorded before and after
a run; a changed source or binary must be investigated before using the run.
The runner does not certify a build merely because the source tree is clean.
Formal experiments use TSUBAME Grid Engine array jobs prepared by
`plan_array.py` and executed by `run_array_task.py`.
Direct interactive-node execution is limited to small development pilots.
Formal configs require `purpose: formal`, a software lock, and
`model_settings.build_revision_verified: true`; the lock must also record
the verified build revision. The direct `run_cli.py` command accepts development
configs only. The array wrapper checks the pinned lock at planning (unless
an offline dry-run), before every task and after execution, using the existing
source, build-cache, parameter and binary hash checks.
Development configs may already specify `software_lock`: the runner checks
the pinned source, cache, parameter and binary hashes before and after running.

Upstream CParty issues 3, 4 and 5 were reported by the author as accepted and
incorporated on 2026-10-08. Pin and inspect the reference actually used before
E0. Keep the completion-domain ambiguity follow-up in implementation notes
until its resolution is known; do not change either implementation here.

The [full bpRNA/E0 plan](full_bprna_array_plan_20261011.md) pins the
implementation's five-patch `800e95a` oracle (ordered 0001,0002,0006,0007,0008)
and LinearCParty `d21618c5` for the new E0 plan.
That single oracle omits completed-WMB
weak closure. Compare `LegacyCompatible` for numerical parity and
`RequireWeakClosure` for model differences, preserving both results.
E0 compares exact inside log Z and PF free energy only, using absolute
roundoff limits 1e-9 and 1e-10 kcal/mol respectively from the implementation's
full-precision audit. Both limits must pass. Original oracle objects provide
the reference; E0 performs no AD, outside, MFE or profile calculation.
The earlier BPP observer remains a historical diagnostic and is not used in
E0 reports. BPP/MFE equality, including the earlier BE41 MFE discrepancy,
is outside this comparison. Default-versus-compatible metrics count closure
effects only when compatible-versus-oracle PF parity passes; failures and
unassessable inputs keep their original denominator.

## Stored observations

`run.json` contains the config and environment metadata, Git provenance,
parameter and input hashes, executable hashes and one entry per job.
Development and historical runs retain job directories. New array jobs create
those directories only on node-local `$TMPDIR` (or explicit `--scratch-dir`,
outside `/gs/bs`). All completed results go into one `task_results.jsonl`, and
commands, stdout/stderr, GNU time, result JSON and BPP exports go into one
`task_outputs.tar` per task. Local job files are removed after their durable
checkpoint; the entire temporary workspace is removed at normal exit.
Commands use local copies of frozen inputs and parameters. The plan freezes all
inputs in one `inputs.tar`, instead of one shared file per RNA.

Each JSONL row preserves argv, status/returncode/error, start/finish times,
elapsed and GNU wall time, peak RSS, stdout/stderr sizes/hashes/bounded tails,
and raw member hashes. The tar is flushed and fsynced before the row is
appended, flushed and fsynced. The row's `archive_end` gives the last committed
tar offset. Resume retains successful, failed and timed-out completed jobs,
repairs an unterminated JSONL tail, truncates uncommitted tar bytes, and reruns
only missing jobs. Hard-killed local work can be lost; a start without a durable
completion is pending. Start/finish events separately recover attempt order.

`elapsed_seconds` measures elapsed wall time around the command, including
the GNU time wrapper. `time_wall_seconds` is GNU time's wall measurement;
`peak_rss_kib` is the peak resident set size of that individual process.
Timing tasks run sequentially, in separate processes, with `node_f=1` and
thread counts at 1. Accuracy-only arrays may choose concurrent CLI processes
through their tier's `workers` setting; do not use these concurrent timings as
isolated runtime measurements. Warm-ups complete before the matching first
measurement, and each seeded repetition block finishes before the next begins.
Configured warm-ups
are logged separately and excluded from the summary's measured-job count.
Do not run measurements alongside other calculations on the same node,
including scaffold generation, implementation checks or another CLI matrix.
Finish input preparation and checks before measuring; use a scheduler allocation
that avoids concurrent calculations on that node.

The default order remains input → mode → program, followed by that combination's
warm-ups and repetitions. An optional config entry
`"order": {"kind": "blocked_random", "seed": 20261011}` makes each repetition
a block: shuffle all input × mode × program combinations using Python
`random.Random(seed)`, continuing the generator between blocks.
Place all warm-ups for a combination directly before its first measured job.
Save the seed/config and complete `planned_jobs` schedule, and record each actual
job start in `execution_order` before launching it. Actual attempt records gain
a finish time and status when their results are saved.

Resume with the original command plus `--resume` after stopping any surviving
child processes. The runner first compares the original config bytes and saved
config with `run.json`, then verifies software paths/state, lock content/hash,
binary hashes, build settings/cache, energy parameters and every source/copied
input hash. A mismatch refuses the resume without changing the run.
The experiment repository may have changed; its current state is recorded for
each resume. Host/environment and resume time are also recorded in
`resume_history`; avoid pooling timing segments from different nodes.
Legacy runs without a cache hash compare their saved selected build settings
and software lock, if present. Their known execution order is recovered from
job command files, with unknown start times left null.

In historical directory runs, `result.json` is the completion marker regardless of success, failure or timeout.
Completed results remain unchanged, including a result saved just before the
manifest update was interrupted. Missing results are rerun in the saved schedule
order; their old job files move to `incomplete_attempts/` so stale outputs or
GNU time files cannot enter the new result. This preserves completed warm-ups
and avoids repeating them on resume. A manifest without finalization fields can
be summarized as `interrupted`; unavailable software/binary change checks are
null, which does not assert that those dependencies stayed unchanged.

Failures and timeouts remain in the manifest and summary. A timeout may lack
a completed GNU time measurement; that value is null rather than zero.
An unavailable binary is recorded as a failed launch through GNU time.
Never drop failed inputs when reporting an experiment's denominator.

The manifest records only selected environment variables relevant to the
calculation and scheduler. Do not dump credentials or a complete environment.
Configs can explicitly set OMP/BLAS thread counts and LC_ALL. The manifest
records the effective values passed to the process, including these overrides.
The explicit copied parameter file overrides any inherited
`CPARTY_FIXED_ENERGY_PARAM_FILE` for executed commands.

CLI time includes every requested calculation and output. LinearCParty also
computes MFE; Raccess calculates common structural profiles. Isolated inside,
outside and additional-consumer costs require an external API driver and must
not be inferred by relabeling these whole-process measurements.

## Array planning, recovery and aggregation

Keep all programs, exact/beam settings, widths, warm-ups and repetitions for
each sequence on the same task. Supplied scaffold variants of that sequence
also stay together. Separate inclusive length tiers into distinct arrays with
their own resources and `h_rt`. Pack by the cost model's estimated work and
memory, leaving headroom below `h_rt`; do not pack by a fixed number of inputs.
Save model coefficients, their basis, input assignments, estimates and full
seeded schedules in `array_plan.json` or checksummed per-task schedules.
The per-job timeout and memory cap remain unchanged unless an explicit tier
`job_limits` overlay is saved. That overlay is validated and copied into the
task config. Length-restricted modes are applied identically in packing and
the schedule. A sequence matrix exceeding its tier budget is rejected.
Cost coefficients are planning assumptions until calibrated with the adopted
build's pilot; they are not empirical or theoretical complexity results.

Each task writes a separate CLI run and `task.json`, including task numbers,
`JOB_ID`, `SGE_TASK_ID`, hostname, requested resources and start/end times for
every attempt. A per-task lock rejects duplicate controllers. Restart only the
interrupted task range with the generated script's `--resume`, after stopping
surviving CLI processes. Resume checks the frozen plan and dependency hashes,
then the JSONL checkpoint; completed failures and timeouts remain observations.
Old directory-layout task plans must be replaced for execution; their results
remain readable.
Check attempt hosts before combining timing segments from different nodes.

After execution or resume, use `aggregate_array.py --plan .../array_plan.json`
to regenerate the run-level manifest. Aggregation reads all planned tasks and
JSONL checkpoints (or historical result markers), including results saved before an interrupted manifest
update. Missing jobs receive explicit statuses and null measurements.
Task counts and planned/completed/missing measurement counts accompany the
ordinary `summarize.py` and `validate_outputs.py` outputs. An active task is
reported as running; an inactive controller without finalization is interrupted,
with unknown end time. Aggregation is explicit, so refresh it after task changes.

### File-count gate before submission

The original full E1 run exhausted the laboratory project file quota and
prevented other users from writing. Never create per-CLI or per-RNA output
files on the shared filesystem. Both CLI and G planners estimate all new shared
entries, including directories, plan files, task metadata/locks, scheduler logs,
and simultaneous atomic-write temporaries, before creating any plan output.
They run `lfs quota -p <project> <filesystem>` and refuse a plan above 10% of
the remaining entries (or a smaller configured fraction). If the query fails,
an explicit `file_budget.fallback_remaining_files` is required.

Use `array.file_budget` in G/embedded configs, or `file_budget` in a separate
CLI array config. Defaults are project 102571, filesystem `/gs/bs`, hard cap
24,000,000 and fraction 0.1. The effective cap is the smaller of the configured
cap and positive live soft/hard file limits. A later live limit is 28,000,000;
the 24,000,000 incident cap remains the conservative planning boundary.
The full configs' fallback is 800,000 entries. A fallback is a recorded estimate,
not a reservation; recheck it and the live quota immediately before submission.
Repeated submissions add scheduler log files beyond the first-attempt estimate.

New CLI tasks have five files under `run/` and no shared per-job directories.
G tasks have three files under `run/`; G publication has `dataset.json` and
one `inputs.tar`. Cohort preparation consumes members without extracting a full
dataset. Readers also support historical directory output and legacy
`results/archives/*.tar[.gz]`. Tar headers are scanned and member data are read
from streams without creating extracted files. Use `archive.tar::tasks/<tier>/<id>/run`
to analyze one legacy task. Whole-array aggregation retains missing observations
but requires memory proportional to the number of planned jobs.

`plan_array.py --dry-run --estimate-only` saves only `array_estimate.json` and
its directory. It checks real input/dependency hashes when source/build paths
are supplied, saves task assignments and the file budget, and is non-executable.
The [full E1 report](packed_array_20261011/array_estimate.json) estimates
10,963 entries for 841 tasks / 699,387 jobs, including temporary headroom.
Local scratch capacity, shared byte quota and task wall time for output transfer
must also fit the requested outputs. Packing does not reduce their byte volume.

The [README](../README.md#formal-experiments-on-tsubame-arrays) specifies the
config schema, memory checks, commands and confirmed resource table.
The [2026-10-11 dry-run record](array_dry_run_20261011.md) contains generated
plans/scripts for the existing bpRNA and RNAcentral configs. Those configs
remain development pilots and are not relabeled as formal measurements.

## Release boundary

Commit scripts, configs, small redistributable inputs and reviewed summaries.
Keep raw runs under ignored `results/runs/`, then archive the runs used in the
paper with checksums and a stable download reference. Figure generation belongs
here; export adopted figures to the manuscript repository with their run IDs.
The experiment release should reproduce its figures without access to the
private manuscript checkout.
