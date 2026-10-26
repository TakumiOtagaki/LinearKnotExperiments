# Supplementary Figure S3: plain output labels (2026-10-26)

The memory legend reads Peak memory and sits in the lower-left margin. Its diamond is open, matching the plotted symbols and caption.

All numerical plot data are byte-identical to ../figure4_resources_20261026_legend/plot_data.json.
RNA calculations were not rerun. The PDF dimensions, data, curves, error bars, scales and axis limits are unchanged.
The PNG was opened and checked for clipping and overlap. Normal text is at least 7 pt.

Reproduce from experiments/ into a fresh output directory:

~~~sh
python3 -B scripts/make_figure4_resources_20261026.py --out results/figures/NEW_DIRECTORY
~~~

The renderer also creates companion figures. This directory retains only Figure S3, its unchanged plot data, captions and checks; temporary companion figures were removed.
The paper exports figure4_synthetic_ratio_si.pdf to figures/.
The manuscript caption is maintained in English and Japanese in the paper repository.
manifest.json records source and output hashes; validation.json records checks and PNG review.
