# NEAT1: fresh 89baea52 analysis

All 51 main/variant conditions and the six additional by-loop conditions were recomputed.
`inventory.tsv` preserves the condition count; `resources.tsv` gives all main/by-loop resource values.
Accessibility, profiles, BPPs, sampling summaries and distances use only the new saved outputs.
`validation.json` gives probability and by-loop checks.
`provenance.json` records input/checkpoint/artifact hashes and commands; its separate `resource_audit` inventories the additional resource checks.
Raw tar members are streamed without extraction.

Reproduce with `python3 -B scripts/analyze_neat1_89baea52.py --out NEW_DIRECTORY`.
Render with `python3 -B scripts/render_non_e1_89baea52.py`.
