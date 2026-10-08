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

The next accuracy pilot will use short bpRNA sequences. The existing local
filtered FASTA contains 26,063 distinct ACGU sequences and enough candidates
in each proposed length bin. Its source-file hash and filtering procedure still
need to be reconciled with the original distribution before formal adoption.
Read the [bpRNA data audit](docs/data_audit_bprna_20261008.md).

For initial long-input timing and RSS measurements, seven ACGU-only sequences
of 4,793–29,106 nt are available in the published LinearCapR RNAcentral panel.
The panel has a fixed public commit and file hash; the audit records sequence
IDs, individual hashes and exclusions. Read the
[RNAcentral data audit](docs/data_audit_rnacentral_20261008.md).
Neither biological cohort has been selected or run in this repository yet.

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

Use a new output directory for each run; existing directories are never
overwritten. Config paths, inputs and output paths supplied on the command line
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

They exercise an expected nonzero exit, a process-group timeout, continuation,
warm-up exclusion, copied inputs and refusal to overwrite an existing run.
