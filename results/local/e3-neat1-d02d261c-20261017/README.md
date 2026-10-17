# NEAT1 fixed rerun audit (2026-10-17)

Regenerate with `python experiments/scripts/analyze_e3_neat1_rerun.py`.
The script streams archives in memory; it does not extract tar files or run RNA software.
It checks all planned jobs against complete JSONL records, all artifact sizes/SHA256 and
archive boundaries, GNU time values, unchanged software/binaries, and valid final locks.
Matched commands differ only in executable and temporary paths; inputs have identical hashes.

`resources.tsv` contains old e2cecbb4 and new d02d261c whole-process wall seconds and
peak RSS (KiB / 2^20 = GiB). Each condition is one run. Missing sampling settings are
not imputed. By-loop runs request lengths 1,4,8,16 and six joint probability components.
`output_comparison.tsv` treats absent BPP pairs as zero and compares raw probabilities
before numerical cleanup. MFE and MEA stdout include their structures and scores.
Sampling compares all 1,000 structures in order with seed 20261016.
All Raccess stdout lines except the measured elapsed_seconds metadata agree exactly.

`byloop_components.tsv` gives full-transcript means and maxima by transcript, beam,
window length and loop type; maximum starts are one-based and inclusive.
`byloop_regions.tsv` summarizes the ten highlighted eight-base windows in Figure S6.
`validation.json` records component-sum, total-accessibility and length-one context checks.
`provenance.json` records consumed artifact hashes and checkpoint/software provenance.
The older e3-neat1-20261016 analysis and its variant/empty-scaffold outputs are retained.
