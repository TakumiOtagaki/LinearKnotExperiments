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
