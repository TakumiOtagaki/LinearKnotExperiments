# E0 oracle update evidence, 2026-10-11

This is a small development check and an unstarted plan, not the bpRNA result.
The [updated plan](../full_bprna_array_plan_20261011.md) defines the adopted
single CParty oracle and PF-only comparison.

- [Version check](e0_version_parity.json): 14 saved inputs, both closure policies,
  28 comparisons; all PF values identical between 15aa7742 and 2456f9b2.
- [Four-witness analysis](e0_witness_analysis.json): four compatible/oracle
  PF matches; one detectable closure difference (nested45).
- [Witness observations](e0_witness_observations.json): all 12 full-precision
  observations and raw artifact hashes.
- [Plan summary](E0/summary.json), [assignments](E0/tasks.tsv) and the five
  saved scripts: 304 inputs, 21 tasks, 912 unstarted jobs. Dependencies verified,
  pending G prevents execution. Historical resource envelopes remain provisional.
- [Verification](verification.json): all 42 experiment tests passed.

Raw development observations live in ignored `results/runs/` and
`results/local/`. The earlier seven-patch lock and older E0 reports remain
historical records. BPP and MFE agreement are outside this E0 comparison.

The initial observer build uses `python3 scripts/build_e0_tools.py` after the
new oracle is built. To repeat the checks with the completed pinned clones,
observer build and latest exact library, choose fresh result directories:

```bash
python3 scripts/check_e0_version_parity.py --out results/local/e0-version-parity-REPEAT
python3 scripts/run_cli.py --config configs/e0_witnesses_800e95a_4patches_20261011.json \
  --software-root deps/LinearCParty-15aa7742 \
  --build-dir build/e0-observers-15aa7742-800e95a-4patches \
  --out results/runs/e0-witnesses-REPEAT
python3 scripts/analyze_e0.py results/runs/e0-witnesses-REPEAT \
  --out results/local/e0-witnesses-REPEAT-analysis.json
python3 scripts/plan_array.py --config configs/e0_bprna_legacy_20261011.json \
  --array-config configs/arrays/e0_bprna_legacy_20261011.json \
  --software-root deps/LinearCParty-15aa7742 \
  --build-dir build/e0-observers-15aa7742-800e95a-4patches \
  --out results/runs/e0-plan-REPEAT --dry-run
```

The build tool refuses an existing output directory. Reuse the completed
observer build to repeat the checks; rebuilding into a new directory requires
updating the lock and observer paths before making new plans. No qsub is needed
for these preparation checks.
