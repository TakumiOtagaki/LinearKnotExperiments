# Synthetic pilot cohort

Generated with `python3 scripts/prepare_pilot.py`, master seed 20261008.
Lengths are 40, 80 and 160 nt with two sequence seeds each. Every sequence
has empty, single-helix and branched scaffold versions. The single helix is a
subset of the branched scaffold; sequences are identical within each group.

The generator imposes nine canonical pairs on a uniform ACGU background.
Scaffolds are supplied geometry controls, not biological reference structures
or predictions from an RNA folding program. They do not exercise depth cap30.
`dataset.json` records provenance, per-record hashes and scaffold descriptors.

The committed files are immutable inputs. Generate a new cohort/config under
fresh paths when changing the seed, geometry or selection rule. Accuracy
measurements use the same saved S and G for every mode.
