# Supplementary Figure S2: plain output labels (2026-10-26)

Memory axes read Peak memory (GiB). Panel headings name the partition function, pairing probabilities and MFE; loop profiles; and loop profiles with window accessibility. Empty G reads Empty scaffold.

All numerical plot data are byte-identical to ../main_4f97fd3_20261021/plot_data.tsv.
RNA calculations were not rerun. The PDF dimensions, data, curves, error bars, scales and axis limits are unchanged.
The PNG was opened and checked for clipping and overlap. Normal text is at least 7 pt.

Reproduce from experiments/ into a fresh output directory:

~~~sh
python3 -B scripts/make_figures.py --review-layout --figure3-only --e2 results/local/e2-4f97fd3-20261021 --e2-revision 4f97fd3 --out results/figures/NEW_DIRECTORY
~~~

The renderer also creates companion figures. This directory retains only Figure S2, its unchanged plot data, captions and checks; temporary companion figures were removed.
The paper exports figure3_synthetic_all_cli_si.pdf to figures/.
The manuscript caption is maintained in English and Japanese in the paper repository.
manifest.json records source and output hashes; validation.json records checks and PNG review.
