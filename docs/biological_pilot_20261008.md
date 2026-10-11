# Biological pilot, 2026-10-08

The author approved the local bpRNA collection identified in the
[acquisition manifest](../manifests/data-sources/bprna-local-20261008.json).
The pilot selects four records in each of [50,100), [100,200) and [200,400) nt,
from ID-sorted pools with Python `random.Random(20261008)`.
The twelve selected IDs and sequence/scaffold hashes are saved in
[the dataset manifest](../data/bprna-pilot-20261008/dataset.json).
No input is replaced in response to a calculation failure or error size.

Both cohorts use LinearPartition-V commit
`b450fb3e63189073b68d385589035f992080aa3a`, beam 200 and dangles 2.
We retain noncrossing, disjoint canonical pairs with posterior probability at
least 0.7, then remove maximal stacked runs shorter than two pairs.
The actual G is saved as parentheses in each input.
We do not use LinearPartition MEA in this pilot: its implementation allocates
quadratic tables, which would burden long-input scaffold construction.
This replaces the earlier MEA/gamma=3 proposal and is a scaffold preparation
choice, rather than a new structure prediction method.
Scaffold wall times are saved separately from the conditional CLI measurements.

The pinned dedicated LinearCParty build is identified by
[the software lock](../manifests/software-pilot-66eaa2f1.json).
The runner checks source state, build cache, energy parameters and binary hashes
before and after a locked run. Three maintained implementation checks passed:
exact rule weights, an exact inside oracle and beam API behavior.
These checks do not substitute for an independent ensemble validation against
the updated upstream CParty. The fixed source still records a completed-WMB
closure follow-up; this run therefore has `purpose: development`.

[The short-input config](../configs/e1_bprna_pilot_20261008.json) specifies
12 inputs × exact/default beam 100 × three CLI programs = 72 measured jobs.
It requests complete sparse BPP output, seven contexts, and accessibility for
window lengths 1, 4, 8 and 16. Default finite-beam local limits are hairpin 30
and case-3 depth 30; the exact comparison measures their combined approximation.
Each job has 180 seconds and an 8 GiB address-space cap. Jobs run sequentially
with thread environment variables set to 1, no warm-up and one measurement.
This is a feasibility pilot; it does not establish timing rankings.

```bash
python3 scripts/run_cli.py --config configs/e1_bprna_pilot_20261008.json \
  --software-root deps/LinearCParty-20261008 \
  --build-dir build/linear-cparty-20261008 \
  --out results/runs/e1-bprna-20261008
python3 scripts/summarize.py results/runs/e1-bprna-20261008 \
  --out results/local/e1-bprna-20261008.tsv
python3 scripts/analyze_accuracy.py results/runs/e1-bprna-20261008 \
  --out results/local/e1-bprna-20261008-accuracy.json
```

CLI times include process startup, all requested calculations and output.
LinearCParty also computes MFE; Raccess includes profile computation.
We do not label these measurements as isolated inside time or add their
separate process times to estimate a reused graph's total cost.

## Long-input range fixed before running LinearCParty

Scaffold construction was initially attempted for all seven preselected RNAcentral
records. Five succeeded; the 22,543 and 29,106 nt records reached the 120-second
scaffold timeout. The candidate manifest retains both failures and the
unchanged source/binary hashes. No replacement record was selected.

Code inspection found that the 2026-10-08 pinned implementation's `sparse_tree` has
14 RMQ columns and a power-of-two array of the same length. Its Euler walk
contains `2*n+1` entries. That pilot conservatively required `2*n+1 < 2^14` for general
interval queries, which gives n <= 8,191. The four longer candidates
were deferred pending the implementation team's capacity adjustment, independently
of their scaffold status. We do not interpret their omitted CLI measurements
as runtime failures or successful runs.

[The long-input config](../configs/e2_rnacentral_pilot_20261008.json) includes
the 4,793, 5,482 and 6,990 nt records, default beam 100 and three CLI programs:
nine measured jobs with 300 seconds/job and a 16 GiB address-space cap.
It retains the four deferred IDs and reasons. Thread settings, output scope
and repetitions match the pilot above. This is a feasibility measurement;
output consistency has no long-input exact reference.

```bash
python3 scripts/run_cli.py --config configs/e2_rnacentral_pilot_20261008.json \
  --software-root deps/LinearCParty-20261008 \
  --build-dir build/linear-cparty-20261008 \
  --out results/runs/e2-rnacentral-20261008
python3 scripts/summarize.py results/runs/e2-rnacentral-20261008 \
  --out results/local/e2-rnacentral-20261008.tsv
python3 scripts/validate_outputs.py results/runs/e2-rnacentral-20261008 \
  --out results/local/e2-rnacentral-20261008-validation.json
```

The run above stopped at 22:41:59 JST when the controlling session ended.
Four jobs have results; the URS00019626CE LinearPKCapR job had started but left
no result, and four jobs had not started. The interrupted directory is kept
unchanged. Two resume configs measure the remaining jobs with the same settings:
[URS00019626CE, LinearPKCapR and LinearPKRaccess](../configs/e2_rnacentral_resume_a_20261008.json)
and [URS000239D8A4, all three programs](../configs/e2_rnacentral_resume_b_20261008.json).
They run one after the other on one node.
The runner now also supports `--resume` for matching interrupted runs; the
historical directory and separate resume runs above remain unchanged.

```bash
for part in a b; do
  python3 scripts/run_cli.py --config configs/e2_rnacentral_resume_${part}_20261008.json \
    --software-root deps/LinearCParty-20261008 \
    --build-dir build/linear-cparty-20261008 \
    --out results/runs/e2-rnacentral-resume-${part}-20261008
done
```

## Rerun on LinearCParty df274e87 (2026-10-09)

The implementation added constant-time scaffold geometry indexes and bounded
WMBP prefix heaps. A dedicated clean checkout of commit `df274e87` was built
with the same Release settings, and the same three maintained checks passed
([software lock](../manifests/software-pilot-df274e87.json)). The planned O(nb)
memory change is documented upstream but not implemented in this commit.
[This config](../configs/e2_rnacentral_df274e87_20261009.json) repeats the nine
RNAcentral jobs with unchanged inputs and limits. Its results are recorded as a
separate software version and do not replace the `66eaa2f1` results.

```bash
python3 scripts/run_cli.py --config configs/e2_rnacentral_df274e87_20261009.json \
  --software-root deps/LinearCParty-20261009 \
  --build-dir build/linear-cparty-20261009 \
  --out results/runs/e2-rnacentral-df274e87-20261009
```

## Rerun on LinearCParty 58a2b796 (2026-10-10)

This version adds the retained-state memory reduction (through `f676c5ea`),
full VP split ranges (`509b5de6`), unscaled log weights, weak closure for
completed WMB by default (`e1bab899`) and the long BE loop correction. The
maintained required suite (68 C++ targets, 7 Python checks) passed on the
dedicated build ([software lock](../manifests/software-pilot-58a2b796.json)).
[This config](../configs/e2_rnacentral_58a2b796_20261010.json) repeats the same
nine jobs and limits. The model change (closure) means outputs are not
expected to equal earlier versions.

```bash
python3 scripts/run_cli.py --config configs/e2_rnacentral_58a2b796_20261010.json \
  --software-root deps/LinearCParty-58a2b796 \
  --build-dir build/linear-cparty-58a2b796 \
  --out results/runs/e2-rnacentral-58a2b796-20261010
```

## Seven-record feasibility preparation (2026-10-11)

The next feasibility matrix includes all seven original RNAcentral candidates.
It removes the previous n <= 8,191 exclusion from `runtime_eligibility` and
retains the initial three-record measurements as separate historical runs.
The 11,832 nt `URS00023C9544` and 13,652 nt `URS000241043C` inputs reuse their
saved G, whose original generation times were 64.98 and 81.82 seconds.
Their input and scaffold hashes remain unchanged.

The two longer records were prepared again using the same LinearPartition-V
commit, wrapper/binary hashes, beam 200, dangles 2, probability threshold 0.7,
minimum stacked helix length 2 and 8 GiB address-space budget.
Only the time budget changed from 120 to 3,600 seconds per sequence.
They ran one at a time, with OMP/BLAS thread settings at 1, on node `r8n10`
(scheduler job `8970385`).

| Record | Length | New G generation | Retained pairs |
|---|---:|---:|---:|
| URS000238C72C | 22,543 nt | success, 137.59 s | 4,465 |
| URS00023A4DB9 | 29,106 nt | success, 182.18 s | 6,539 |

The [retry dataset](../data/rnacentral-retry-20261011/dataset.json) contains
the new inputs, sequence/G/input hashes, descriptors and generation times.
The [retry provenance manifest](../manifests/rnacentral-scaffold-retry-20261011.json)
also preserves the exact preparation/per-record commands, result records and
raw-log hashes. Source, wrapper and binary were unchanged before and after.
The earlier 120-second failures remain in the initial dataset; neither record
was replaced. G generation times are separate from conditional CLI measurements.

```bash
python3 scripts/prepare_biological_inputs.py --kind rnacentral \
  --source-manifest manifests/data-sources/rnacentral-panel-20261008.json \
  --linearpartition deps/LinearPartition-20261008/linearpartition \
  --out data/rnacentral-retry-20261011 \
  --timeout-seconds 3600 --ids URS000238C72C URS00023A4DB9
```

The [NEXT config](../configs/e2_rnacentral_long_NEXT_20261011.json) is a template
for seven inputs × default beam 100 × LinearCParty, LinearPKCapR and
LinearPKRaccess: 21 measured jobs. Hairpin and case-3 depth limits remain 30/30.
Requested outputs match the `58a2b796` config: complete sparse BPP, seven
contexts and accessibility for window lengths 1, 4, 8 and 16.
The budget is 3,600 seconds/job and 64 GiB address space, with zero warm-ups,
one measurement and `blocked_random` order with seed 20261011.
This is a feasibility measurement, with no long-input exact reference or
timing ranking from a single repetition.

`software_lock: manifests/software-pilot-NEXT.json` is deliberately unresolved.
Replace this value once the author fixes the next implementation commit and
verifies its dedicated build, including support for the full selected range;
then set `model_settings.build_revision_verified` to true.
No conditional CLI measurement from this template has been launched.
Finish all scaffold generation and implementation checks first, and run the
matrix without other calculations on the same node, as required by
[the protocol](experiment_protocol.md).
