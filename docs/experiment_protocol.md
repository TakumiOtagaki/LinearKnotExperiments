# Experiment protocol

Use the same sequence, supplied scaffold, thermodynamic parameter file,
temperature, local limits and requested output for a matched comparison.
Save the actual scaffold, rather than only the parameters used to create it.

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

Failures and timeouts remain in the manifest and summary. A timeout may lack
a completed GNU time measurement; that value is null rather than zero.
An unavailable binary is recorded as a failed launch through GNU time.
Never drop failed inputs when reporting an experiment's denominator.

The manifest records only selected environment variables relevant to the
calculation and scheduler. Do not dump credentials or a complete environment.
The explicit copied parameter file overrides any inherited
`CPARTY_FIXED_ENERGY_PARAM_FILE` for executed commands.

## Release boundary

Commit scripts, configs, small redistributable inputs and reviewed summaries.
Keep raw runs under ignored `results/runs/`, then archive the runs used in the
paper with checksums and a stable download reference. Figure generation belongs
here; export adopted figures to the manuscript repository with their run IDs.
The experiment release should reproduce its figures without access to the
private manuscript checkout.
