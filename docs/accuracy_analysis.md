# Accuracy analysis

The pilot is a development experiment. Source/build and reference semantics
must be frozen before these observations become manuscript measurements.
The analyzer processes one named RNA record per input and uses only the
standard library. Each measured job remains recorded, including failed
jobs; missing reference output does not become zero error.

## Run and analyze

The committed cohort has six synthetic sequences, each with three scaffold
choices, at 40, 80 and 160 nt. Background bases are uniform ACGU with nine
canonical pairs imposed. These are geometric controls, not biological data.
`data/pilot-v1/dataset.json` records seeds, hashes and scaffold descriptors.
The CASE3 depth cap of 30 is inactive on this cohort; it is not a depth stress test.

```bash
python3 scripts/run_cli.py --config configs/e1_pilot.json \
  --software-root /path/to/LinearCParty \
  --build-dir /path/to/LinearCParty/build \
  --out results/runs/e1-pilot-001
python3 scripts/analyze_accuracy.py results/runs/e1-pilot-001 \
  --reference-mode exact --out results/local/e1-pilot-001-accuracy.json
python3 scripts/summarize.py results/runs/e1-pilot-001 \
  --out results/local/e1-pilot-001-time.tsv
```

There are 324 jobs: 18 inputs, six modes and three CLIs. The runner has not
executed this matrix as part of preparation. Use a scheduler for the full run.
Each job has a 60-second timeout, with no warm-up and one measured repetition.
These settings check feasibility and logging; they do not establish speed rankings.
The config applies and records single-thread environment values and `LC_ALL=C`.

The four `nocap` modes set both CLI geometry caps to 1,000,000, above every
possible length/depth in this bounded cohort. They isolate beam-width effects
relative to `--exact`, conditional on matching grammar/weights.
The default-cap b=100 mode is an additional setting-sensitivity observation.
Changing caps can change the retained beam as well; the finite-beam difference
is not a cap-only partition-mass decomposition. `--exact` currently ignores
the beam geometry caps.

To compare two finite settings, choose the intended reference explicitly:

```bash
python3 scripts/analyze_accuracy.py results/runs/e1-pilot-001 \
  --reference-mode beam100_nocap --out results/local/e1-pilot-001-settings.json
```

Inputs/configs can be reproduced with `scripts/prepare_pilot.py` using fresh
paths. Existing paths are refused. Stored input hashes are authoritative.

## Numerical definitions

Let U be the positions unpaired in G. For each i in U, its base-state
distribution has categories unpaired and paired to each other position.
The primary BPP error is its total variation, averaged over U:

```text
e_BPP = (sum_{i<j; i,j in U} |p_ij - p*_ij|
         + 0.5 sum_{i in U} |q_i - q*_i|) / |U|.
```

The starred quantities belong to the chosen reference. This uses the union
of sparse coordinates and does not materialize an n-by-n matrix. BPP export
must have cutoff zero; an omitted sparse pair then means zero. Duplicate
coordinates, invalid probabilities, incident mass exceeding one and missing
fixed pairs are rejected. The coordinates are one-based.

Context error is mean per-base total variation over the seven contexts on U.
Context-specific MAEs are retained. All seven rows and all positions must
be present. The probabilities must sum to one at each position, and fixed
scaffold positions must have Stem probability one.

For each requested window length, the primary accessibility domain contains
only windows entirely in U. Other windows necessarily have probability zero.
The analyzer records total and eligible window counts, probability MAE/max
error, and opening-energy MAE on windows where both probabilities are positive.
Reference-positive/candidate-zero and reference-zero/candidate-positive
windows are counted separately. Empty denominators produce null, not zero.
All requested coordinates must be present; missing accessibility is an artifact
error, not an implicit zero. Raccess zero-based starts are converted to one-based.
The calculation uses the emitted joint probability, not a product of q_i.

Thermodynamic errors use the full-precision Raccess metadata:

```text
delta_F = F_candidate - F_reference
delta_log_Z = -delta_F / RT
```

Both RT values must match. LinearCParty's printed free energy has only four
decimal places and its MFE has two; those displays are not used for strict
numerical verification. BPP output itself has full precision. A future API
driver must supply full-precision MFE/decoder values and isolated stage timings.

The analyzer also compares Stem to BPP paired mass, accessibility of length
one to q_i, and all window probabilities to the minimum unpaired marginal.
The default diagnostic tolerance is 1e-9. These are consistency checks; they
do not prove unique structure decomposition. Public profile/accessibility
CLIs clip values before printing, so output checks cannot detect every
invalid intermediate value. Independent enumeration/API checks belong to E0.

## Reporting boundary

Accuracy references use measured repetition zero, without silently choosing
another successful reference. Warm-ups are excluded. Unanalyzed programs such
as sampling and MEA remain listed alongside the total measured-job count.
No across-RNA aggregate or confidence interval is produced yet: summarize
per RNA first, and preserve the paired scaffold/beam blocks when resampling.
The analysis records its own Git state, script hashes and the source manifest
hash, as well as the output hashes observed during analysis. The experiment
commit of the original run is kept separately; reanalyzing an old run does
not mean the new analyzer existed at that original commit.

Whole CLI time includes startup, requested calculations and output.
LinearCParty also computes MFE; Raccess uses the common profiles computation.
These timings are not isolated inside or marginal aggregation timings.
Do not add repeated CLI wall times and describe the sum as a shared-graph run.

The next analysis additions are independent sampling-frequency checks,
decoder structures rescored on the same reference posterior, an external API
driver for stage timings/work counters, and a locked biological input cohort.
