# Figure 4: readable legends and completion marks (2026-10-26)

`figure4.pdf` / `figure4.png`: 178 × 58 mm, PNG at 300 dpi.
Axes and headings use 7 pt; legends and completion-row labels use 6.5 pt.
Panel letters use 9 pt. No displayed text is below 6.5 pt.

A/B show the seven RNAcentral RNAs. The legend names the calculated outputs:
`partition function + pairing probabilities + MFE`, `loop profiles`, and
`loop profiles + window accessibility`. The manuscript caption maps these to
`fold`, `profile`, and `access`, respectively.
C/D are headed `bpRNA and RNAcentral RNAs` and use the method labels
`CParty, partition function`, `LinearKnot exact`, `LinearKnot, beam width 100`,
`Spark, MFE`, and `LinearKnot, partition function + MFE`.
Memory axes read `Peak memory (GiB)` and `Peak memory (MiB)`; the caption
identifies this as maximum resident set size (RSS).

Below C/D, one row per method marks input lengths with a single vertical
symbol, colored like that method's resource curve. Labels give
`CParty, not completed (48 of 107)` and `Spark, not completed (4 of 107)`.
The rows share labels across C/D; the marks align with each panel's length axis.
An input appears once per panel regardless of how many causes were recorded
across repetitions. The marks carry no time or memory values.
Cause-specific symbols, their legend, and the old explanatory note are removed.
The existing SI tables S21 (`tab:si_cparty_real`) and S22
(`tab:si_spark_resources`) retain the cause breakdown and counts.

All saved plot data are identical to `../figure4_resources_20261026_final58/plot_data.json`.
This includes all resource measurements, group medians, and synthetic ratios.
The generated SI ratio PNG also has identical pixels; its adopted export is
unchanged. No RNA calculations were run.

Reproduce from `experiments/` into a new directory:

```sh
python3 -B scripts/make_figure4_resources_20261026.py --out results/figures/NEW_DIRECTORY
```

`manifest.json` records input, renderer, and output SHA-256 hashes.
`validation.json` records dimensions, font sizes, counts, and visual review.
The PNG was opened to check all headings, axes, legends, marks, and labels for
overlap, clipping, and readability. The manuscript exports this figure to
`figures/figure3_main.pdf`; its Figure 4 caption is maintained in English and
Japanese in the paper repository.
