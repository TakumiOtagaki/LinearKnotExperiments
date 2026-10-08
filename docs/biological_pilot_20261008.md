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

Scaffold construction was attempted for all seven preselected RNAcentral
records. Five succeeded; the 25,543 and 29,106 nt records reached the 120-second
scaffold timeout. The candidate manifest retains both failures and the
unchanged source/binary hashes. No replacement record was selected.

Code inspection found that the pinned implementation's `sparse_tree` has
14 RMQ columns and a power-of-two array of the same length. Its Euler walk
contains `2*n+1` entries. We conservatively require `2*n+1 < 2^14` for general
interval queries, which gives n <= 8,191. The four longer candidates are
deferred pending the implementation team's capacity adjustment, independently
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
Five jobs have results; the URS00019626CE LinearPKCapR job had started but left
no result, and four jobs had not started. The interrupted directory is kept
unchanged. Two resume configs measure the remaining jobs with the same settings:
[URS00019626CE, LinearPKCapR and LinearPKRaccess](../configs/e2_rnacentral_resume_a_20261008.json)
and [URS000239D8A4, all three programs](../configs/e2_rnacentral_resume_b_20261008.json).
They run one after the other on one node.

```bash
for part in a b; do
  python3 scripts/run_cli.py --config configs/e2_rnacentral_resume_${part}_20261008.json \
    --software-root deps/LinearCParty-20261008 \
    --build-dir build/linear-cparty-20261008 \
    --out results/runs/e2-rnacentral-resume-${part}-20261008
done
```
