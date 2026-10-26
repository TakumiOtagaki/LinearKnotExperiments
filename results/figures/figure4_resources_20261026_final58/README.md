# Figure 4: resource comparisons (2026-10-26)

`figure4.pdf` / `figure4.png`: 178 × 58 mm, 300 dpi PNG, ordinary text 7 pt, method legends 6.5 pt, minimum 6 pt.
A/B retain the adopted E2 RNAcentral PF+BPP+MFE, profile, and access measurements.
C/D add the newly verified 107-real-RNA comparison. Each faint point is an input's three-run median; lines join medians across successful inputs within each bpRNA group or individual RNAcentral inputs. Methods can have different successful sets. Exact stops at 2,158 nt by design. Failure strips use length positions only, without elapsed or RSS values.
CParty failed on 48/107 inputs (43 nonfinite, one timeout-only, three mixed timeout/signal 11, one allocation failure); Spark failed on four. Signal 11 is not identified as memory exhaustion.

**Output scope:** CParty, exact and beam100 compare PF only with matched compatibility closure. Spark computes MFE. The new fold100 command includes PF+MFE without BPP; the existing A/B fold is a separate PF+BPP+MFE measurement. No structure/energy parity or speedup from censored runs is claimed.

`figure4_synthetic_ratio_si.pdf` / `.png` is the former Figure 4C, now SI Figure S3 (165 × 42 mm, minimum 7 pt); its values are unchanged. Figure S2 retains the full synthetic curves.

再現は experiments/ で `python3 -B scripts/make_figure4_resources_20261026.py --out results/figures/NEW_DIRECTORY`。
入力・集計・検証のhashはmanifest.json、全測定の検証はdocs/comparison_resources_20261026_completed_validation.json。
英日の正式captionは論文repoのmanuscript/latex/sections/results.texとmanuscript/results_ja.md。新比較のfoldはPFとMFEでありBPPを含めない。従来のBPP込み測定と区別する。

Opened both final PNGs and normal/diff PDF pages. No clipping or text/legend overlap was seen. Main figure typography was checked from actual PDF text sizes. No RNA program was run while aggregating/rendering.
