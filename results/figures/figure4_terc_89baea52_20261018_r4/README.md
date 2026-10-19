# figure4_terc: 89baea52 (20261018)

All probabilities, structures and sample counts were recomputed with 89baea52.
Raw archive: `results/local/small-rerun-89baea52-20261018/terc_data.tar`.
The archive contains the frozen input and new packed-task outputs.
Only node-local TMPDIR is used for extraction.
The caption and JSON report record the derived values and figure checks.

Reproduce from the experiment repository:

```bash
python3 -B scripts/render_non_e1_89baea52.py --examples-only --suffix _r4
```

## Round-four changes

This is manuscript Figure 3 (Terc), retaining filename figure4.pdf and size 178 × 100 mm. A is a positional pair-probability overview. O / I leaders identify both arms of the outer / inner reference groups, with the same letters in D. Pair styles match Figure 1 and S4.

The proposed bilingual caption is in `caption.md`. `validation.json` records actual PDF font sizes and layout checks. PNG inspection is recorded below after generation.

## PNG inspection and final PDF checks — 2026-10-19

- Figure 3: `figure4.png` を実際に開いて確認。O/Iの四本の指示線が両腕の外側／内側の参照群を指し、Dの同じ文字へ対応する。

Ghostscript `txtwrite`, `TextFormat=0` でPDF内の実際の文字サイズを検査した。
すべて最小7 ptで、6 pt未満の文字はない。PDFの掲載幅は本文178 mm、SI 165 mmで、各図の高さは旧版と同一。
PNGは300 dpi。原稿への挿入後のページ割り確認は原稿担当が行う。
