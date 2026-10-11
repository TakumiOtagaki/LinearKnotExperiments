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
The current runner accepts development configs only; a formal mode that
checks the source/build lock must be added before manuscript collection.
Development configs may already specify `software_lock`: the runner checks
the pinned source, cache, parameter and binary hashes before and after running.

Upstream CParty issues 3, 4 and 5 were reported by the author as accepted and
incorporated on 2026-10-08. Pin and inspect the reference actually used before
E0. Keep the completion-domain ambiguity follow-up in implementation notes
until its resolution is known; do not change either implementation here.

## Stored observations

`run.json` contains the config and environment metadata, Git provenance,
parameter and input hashes, executable hashes and one entry per job.
Every job directory contains its command, stdout, stderr, GNU time output and
result JSON; LinearCParty additionally exports `bpp.tsv` in the smoke.
Commands use copied inputs and parameters to isolate them from later edits.

`elapsed_seconds` measures elapsed wall time around the command, including
the GNU time wrapper. `time_wall_seconds` is GNU time's wall measurement;
`peak_rss_kib` is the peak resident set size of that individual process.
Measured jobs run sequentially, in separate processes. Configured warm-ups
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

`result.json` is the completion marker regardless of success, failure or timeout.
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

## Release boundary

Commit scripts, configs, small redistributable inputs and reviewed summaries.
Keep raw runs under ignored `results/runs/`, then archive the runs used in the
paper with checksums and a stable download reference. Figure generation belongs
here; export adopted figures to the manuscript repository with their run IDs.
The experiment release should reproduce its figures without access to the
private manuscript checkout.
