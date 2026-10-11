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
legacy manifests, interrupted summaries and scaffold subset preparation.
