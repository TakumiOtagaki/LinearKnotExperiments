# figure1: 89baea52 (20261018)

All probabilities, structures and sample counts were recomputed with 89baea52.
Raw archive: `results/local/small-rerun-89baea52-20261018/figure1_data.tar`.
The archive contains the frozen input and new packed-task outputs.
Only node-local TMPDIR is used for extraction.
The caption and JSON report record the derived values and figure checks.

Reproduce from the experiment repository:

```bash
python3 -B scripts/render_non_e1_89baea52.py --examples-only --suffix _r4
```

## Round-four changes

Figure 1 keeps its 178 × 145 mm size. G is solid grey, new non-crossing pairs are dotted blue, and pairs crossing G are dashed orange. Internal uses sparse / hatching in both the context track and key. The sample labels explain the 975 / 6 draws; an arrow marks a new non-crossing pair in the second sample.

The proposed bilingual caption is in `caption.md`. `validation.json` records actual PDF font sizes and layout checks. PNG inspection is recorded below after generation.

## PNG inspection and final PDF checks — 2026-10-19

- Figure 1: `figure1.png` を実際に開いて確認。975回／6回の選択ラベル、青点線の追加対への矢印、Internalの疎な斜線と凡例を確認した。

Ghostscript `txtwrite`, `TextFormat=0` でPDF内の実際の文字サイズを検査した。
すべて最小7 ptで、6 pt未満の文字はない。PDFの掲載幅は本文178 mm、SI 165 mmで、各図の高さは旧版と同一。
PNGは300 dpi。原稿への挿入後のページ割り確認は原稿担当が行う。
