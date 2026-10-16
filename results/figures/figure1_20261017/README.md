# Figure 1 revision — 2026-10-17

`figure1.pdf` and `figure1.png` are 178 × 145 mm at final main-text width.
The PNG is 300 dpi.
The prior 183 × 150 mm figure becomes 145.90 mm tall at 178 mm width, so this revision does not increase its published height.
`caption.md` is a proposed caption with one sentence per line.
`validation.json` records input SHA-256 hashes, values derived from saved outputs, the script hash, and PDF/PNG checks.

## Reproduce

From the experiments repository, use a node-local temporary directory for Matplotlib:

```sh
python3 - <<'PYTHON'
import os, subprocess, tempfile
with tempfile.TemporaryDirectory(prefix="figure1-", dir="/tmp") as scratch:
    subprocess.run([
        "python3", "scripts/make_figure1.py",
        "--data", "results/figures/figure1_20261016/data",
        "--out", "results/figures/figure1_20261017",
    ], env=dict(os.environ, MPLCONFIGDIR=scratch, PYTHONDONTWRITEBYTECODE="1"), check=True)
PYTHON
```


Dependencies: Python, Matplotlib, Pillow, pypdf and Ghostscript (`gs`).
No folding CLI, network access or new experimental run is needed.
The script writes only five final files to the requested output directory.

## Inputs and implementation version

All scientific values come from `results/figures/figure1_20261016/data/`, saved for the 2026-10-16 Figure 1.
The source README records `build/linear-cparty-febd2bc1`, exact mode, one thread, Turner 2004 parameters (`params/rna_Turner04.par`).
It states that those outputs are identical to the versions adopted at that time; this figure revision does not independently test cross-version equivalence.
The source record originally came from `results/local/pk-recovery-bprna-20261016/pk_sel.json`.
No saved input or old figure has been modified.

- `in.txt` supplies the RNA sequence and actual scaffold; `rec.json` supplies the bpRNA identifier and reference extension in A.
- `lc.out` supplies the exact MFE structure; `bpp.tsv` supplies pair probabilities; `capr.out` supplies seven context probabilities.
- `racc.out` supplies the exact-mode/scaffold metadata and probabilities for 1, 4, 8 nt windows; starts are converted from 0-based output to 1-based coordinates.
- `samp.out` supplies all unique structures and counts from `LinearPKSampling --exact -k 1000 --seed 20261016 --unique`.
- `mea.out` supplies the exact MEA structure and γ; the source README records `LinearPKMEA --exact --gamma 1`.
- Original ensemble commands: `LinearCParty --exact --bpp-cutoff 0 --bpp bpp.tsv`; `LinearPKCapR --exact`; `LinearPKRaccess --exact --access-len 1,4,8 --probability`.

The input has 63 bases and 11 scaffold pairs.
The caption says canonical nested reference pairs because the original reference also has a noncanonical nested pair that is absent from the saved scaffold.
A illustrates the saved reference extension; C contains model predictions, not experimentally established structures.
B is a schematic, not a measured beam run: its example width is 3.
The two spans in A are computed from the two pair sets, and the dotted line intersects 2 bands in their connected crossing group.
The drawing is an example of density-2, not a complete graphical definition of all allowed structures.

## Derived display values

The pair plot uses the strict display cutoff p > 0.01 (18 visible pairs).
The two samples have counts 975 and 6, leaving 19 of 1000 draws unshown.
The second sample is selected to contain a new pair that does not cross G, rather than being the second most frequent structure overall.
PseudoLoop's direct label marks bases 39–51, the longest contiguous region where it is the largest context probability.
Accessibility uses saved windows of 1, 4, 8 nt; no new loop-specific accessibility is claimed.

## Review changes

- F1-1 / G1: remove sequences, dot-bracket strings, suboptimal lists and the in-figure methods paragraph; use 7–8 pt text and 9 pt panel letters at 178 mm width.
- F1-2: keep only the five ensemble-derived outputs in the blue frame; draw a separate MFE branch from the rearranged decomposition.
- F1-3 / L16 / L17: mark both band spans and their overlap in A; put the definition and scope of density-2 in the caption.
- F1-4: show rank increasing upward; define the sum of log inside weight and estimated surrounding log weight in the caption.
- F1-5: compute shown/other sample counts and define the strict probability cutoff, window lengths, line styles and coordinates in the caption.
- G2 / L28: label PseudoLoop directly in its region; use dashed crossing arcs and discarded-cross symbols as cues beyond color.
- G3: give A, B and C independent bold 9 pt labels with consistent offsets.
- L21: remove suboptimal outputs from the main figure.

## Validation and handoff

Saved-data consistency checks passed: actual scaffold equals the recorded nested subset, all its pairs have BPP 1, every sampled/decoded structure retains G, window row counts match the RNA length, and context probabilities sum to one (maximum error 1.11e-16).
Ghostscript extracted embedded PDF sizes [7.0, 7.3, 7.5, 8.0, 9.0] pt; minimum 7 pt.
PDF media box: 178.000 × 145.000 mm.
PNG: 2102 × 1712 pixels, 299.999 dpi.
The rendered-text bounds check found no text outside the page.
The final 300 dpi PNG was opened with `view_image` and inspected after the last render.
No text overlaps, clipping or out-of-frame labels were found.
The lower frame was moved clear of the x-axis label, the crowded 60 tick was removed, and the stage arrow was separated from the vertical rank label.
The band brackets, direct PseudoLoop label, sample counts and separate MFE branch are legible.
Reproducing the figure resets the visual-check field to pending; open each newly rendered PNG before use.

The manuscript and its Japanese sources have not been edited.
The manuscript owner must insert this PDF and adapt the caption, then check pagination and the final typeset text sizes.
The caption must replace the old statement that the sample diagrams show the top three draws.
A third-band exclusion example is left out to keep A limited to G, G′ and R as the author requested; the caption gives the general condition.
No scientific claim or reference is added and no author decision is required to reproduce the figure.
The author may choose to redraw the same layout in PowerPoint.
Because this task is restricted to its assigned script/output/log, the manuscript owner should copy the substantive changes and validation into `docs/progress.md`.
