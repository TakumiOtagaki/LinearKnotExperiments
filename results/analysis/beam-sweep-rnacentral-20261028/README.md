# RNAcentral beam-width and local-limit sensitivity (2026-10-28)

All 49 frozen formal settings completed and passed output validation; all 91 comparisons are available.
The seven sequence/scaffold inputs match the original E2 hashes.
See ../../../docs/beam_sweep_rnacentral_20261027.md for definitions, full provenance and commands.

- observations.tsv: all per-RNA/beam full-precision scalars, signed free-energy differences, TV summaries and resources, including limits and status.
- comparisons.tsv: b400 and adjacent-width comparisons; mean/p90/maximum TV over scaffold-free bases, including the unpaired state, without a pair-probability cutoff.
- by_beam.tsv: ranges across the seven RNAs; not percentiles across RNAs.
- b100.tsv: the seven b100 observations.
- summary.json: analyzer hash, frozen plan hashes, counts and interpretation.

The b400 reference is finite, not exact. Beam width and both local limits vary together.
Each time is one single-thread whole-fold observation on a shared node, without warm-up.
No RNA calculations were run during this aggregation. Large raw outputs remain in ignored results/runs/.
