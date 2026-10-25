# Figure 1 — seven structural contexts, 2026-10-25

[PDF](figure1.pdf)・[PNG](figure1.png)（300 dpi）。
178 × 145 mm、高さの増分は 0 mm。最小文字 7 pt。
A 下段に一つの架空の構造を追加した。
A 上部の scaffold と extension を横に並べ、完成構造と交差 band の下に配置した。
本文からの参照は **Figure 1A (lower sketch)**（**Figure 1A 下段**）。

[PDF](figure1.pdf) and [PNG](figure1.png) (300 dpi).
The canvas remains 178 × 145 mm; the height increase is 0 mm.
The minimum embedded text size is 7 pt.
The new fictional structure is in the lower part of A.
The scaffold and extension above it are side by side, followed by their union and the band spans.
Refer to it as **Figure 1A (lower sketch)**.

## Caption addition / caption への追加案

A下段の模式図は架空の一つの構造を示す。
各円は一つの塩基で、その色は構造文脈を示す（Cの凡例）。
橙色の非対合塩基は、交差するbandの間の隙間であるPseudoLoopに属する。

The lower sketch in A shows one fictional structure.
Each circle is one base, colored by its structural context (key in C).
Orange unpaired bases belong to PseudoLoop, the gaps between crossing bands.

追加位置は既存の A の説明の最後。[caption.md](caption.md) も同じ案。
LaTeX は編集していない。高さを変えていないため、図の占有寸法は変わらない。
本文 8 ページの確認と参照文の反映は本文担当に引き渡す。

Add these sentences at the end of the existing description of A.
The same proposal is in [caption.md](caption.md).
LaTeX was not edited. The figure has the same manuscript footprint.
The manuscript editor can add the reference and check the eight-page manuscript.

## Structure and classification / 構造と分類

これは形だけの例であり、測定値や構造予測ではない。
配列・エネルギー・確率は設定していない。
灰実線は G、青点線は G と交差しない追加対、橙破線は G と交差する追加対。
各円は一つの塩基で、色は C と図 3 の既存の七色に合わせた。
G と G′ はそれぞれ交差対を持たず、対合塩基を共有しない。
H 型 pseudoknot は各二対の二つの band からなる。

This example illustrates shape; it is not a measurement or a structure prediction.
It has no sequence, energy or probability values.
Solid grey arcs show G; dotted blue arcs show new pairs that do not cross G;
dashed orange arcs show new pairs crossing G.
Each circle is one base; its color uses the existing seven-color key in C and Figure 3.
G and G′ each have no crossing pairs and share no paired bases.
The H-type pseudoknot has two bands of two pairs each.

- Length / 長さ: 37 nt.
- G: `.(.(.(.(...).......).)).((.....))....`
- R = G ∪ G′: `.(.(.(.(...).[...].).)).((.[[..)).]].`
- G pairs / G の対: `[(2, 23), (4, 22), (6, 20), (8, 12), (25, 33), (26, 32)]`
- G′ pairs / G′ の対: `[(14, 18), (28, 36), (29, 35)]`
- H-type bands / H 型の band: `[[(25, 33), (26, 32)], [(28, 36), (29, 35)]]`
- New noncrossing pair / 交差しない追加対: `[(14, 18)]`

| Context | Positions (1-based) / 塩基位置 | Color / 色 |
|---|---|---|
| Stem | 2, 4, 6, 8, 12, 14, 18, 20, 22, 23, 25, 26, 28, 29, 32, 33, 35, 36 | `#2a78d6` |
| Hairpin | 9, 10, 11, 15, 16, 17 | `#1baf7a` |
| Bulge | 3 | `#eda100` |
| Internal | 5, 21 | `#e87ba4` |
| Multiloop | 7, 13, 19 | `#6250d6` |
| Exterior | 1, 24, 37 | `#b9b8b1` |
| PseudoLoop | 27, 30, 31, 34 | `#eb6834` |

分類は手で割り当てず、LinearKnot **v0.2.0** の形状分類を実行した。
`git archive v0.2.0` で参照コードと fixture だけを node-local に取り出した。
版の commit は `f93b55ed9749dd1c198cd1f1809aece50f12f7d8`。
`test/tools/pseudoloop_context_reference.py` の `classify(37, completed_pairs)` を使用した。
参照コードの SHA-256 は `30547e8701d3c28b7e56e97384b207d405157e0c6df4b7cf71d605892740a48f`。
作図スクリプトはこの hash を確認してから分類し、その出力で各円を塗る。
参照コードの 14 fixture・545 非対合窓の検証も成功した。
全塩基の分類と判定根拠は [structural_contexts.json](structural_contexts.json) にある。

The labels were computed with the unmodified shape classifier from LinearKnot **v0.2.0**.
Only the reference code and its fixtures were extracted to node-local storage with `git archive v0.2.0`.
The commit is `f93b55ed9749dd1c198cd1f1809aece50f12f7d8`.
The call is `classify(37, completed_pairs)` in `test/tools/pseudoloop_context_reference.py`.
Its SHA-256 is `30547e8701d3c28b7e56e97384b207d405157e0c6df4b7cf71d605892740a48f`.
The renderer checks this hash and uses the returned labels to color the circles.
The reference also passed all 14 fixtures and 545 fully unpaired windows.
[structural_contexts.json](structural_contexts.json) saves every label and its geometric explanation.

PseudoLoop は塩基 27、30–31、34 の三つの隙間だけに現れる。
各隙間の両端は交差する別 band の対に属する。
ヘアピン 9–11・15–17、バルジ 3、内部ループ 5・21、multiloop 7・13・19、外部 1・24・37 も含む。
band の二対が閉じる multiloop の追加例は、図を小さく保つため入れていない。
この参照分類は形の分類であり、任意の構造の folding/admissibility の証明ではない。

PseudoLoop occurs only in the three gaps at bases 27, 30–31 and 34.
Each gap has boundaries from different crossing bands.
The example also includes hairpins at 9–11 and 15–17, a bulge at 3,
an internal loop at 5 and 21, a multiloop at 7, 13 and 19, and exterior bases at 1, 24 and 37.
The optional multiloop between two pairs of one crossing band is omitted to keep the sketch small.
This reference classifies shape; it does not certify folding or model admissibility for arbitrary pair maps.

## Reproduce / 再作図

experiments/ から実行する。参照コードは node-local のコピーを指定する。
六つの生成物を保存する。RNA の再計算や外部 repo の編集は行わない。

Run from experiments/ with the extracted reference code.
The renderer writes six files. It does not rerun RNA folding or modify an external repository.

~~~sh
MPLCONFIGDIR=/tmp/figure1-fb6-mpl python3 -B scripts/make_figure1.py \
    --data results/local/pseudoloop-0.2.0-20261023/examples/figure_data.json \
    --context-reference /tmp/figure1-reference/test/tools/pseudoloop_context_reference.py \
    --out results/figures/figure1_fb6_20261025
~~~

参照コードの取り出しと fixture 検証（論文 repo の root から実行し、既存の空ディレクトリを指定）:

Extract the reference and run its fixtures from the paper repository root
(use an existing empty node-local directory):

~~~sh
git -C ../LinearKnot archive v0.2.0 \
    test/tools/pseudoloop_context_reference.py \
    test/fixtures/pseudoloop_context_change.json | tar -x -C /tmp/figure1-reference
python3 -B /tmp/figure1-reference/test/tools/pseudoloop_context_reference.py
~~~

Dependencies: Matplotlib, NumPy, Pillow, pypdf, Ghostscript.
入力は既存の 0.2.0 集計。B と C の内容・配置・数値は保持する。
使用データ・作図コード・参照コードの hash は [validation.json](validation.json) に記録する。

The input is the adopted saved 0.2.0 data.
B and C keep their content, layout and values.
[validation.json](validation.json) records hashes for the data, renderer and classifier.

## Feedback review / 注釈の確認

`feedbacks/main_satoken_10252232.pdf` と s-1.png〜s-8.png を一枚ずつ開いて確認した。
pypdf で `/Annots` の本文注釈（Highlight / StrikeOut）と座標を取得した。
対象語句・コメント・コメントのないマークの判断根拠は [feedback_annotations.json](feedback_annotations.json) に保存する。
リンクと Popup は独立した本文注釈として数えない。
実際の markup は 32 件（Highlight 30、StrikeOut 2）。依頼文の 31 件とは一件差があるため全件を記録した。

Opened s-1.png through s-8.png individually and checked the annotated PDF.
pypdf supplied `/Annots` markup types, comments and coordinates.
[feedback_annotations.json](feedback_annotations.json) records target text and the grounds for interpreting uncommented marks.
Links and popups are not counted as independent text markup.
The PDF has 32 markup objects (30 highlights and two strikeouts); all were retained despite the request's count of 31.

- p1: Abstract の二か所の削除、CParty、時間・メモリの比較、pseudoknot のクラス、IPknot/ThreshKnot の重複、HFold の構造集合。
- p2: Spark、CapR/Raccess、time、交差の式、scaffold/with、新しい対、未定義の集合、Figure 1B の説明の細かさ。
- p3: Figure 1 全体を確認。本文 markup はない。
- p4: **“seven structural contexts” に “図”**。今回の A 下段で対応。他は TV、計算量、設定の説明に関する指摘。
- p5: GitHub bugs、centroid 検証、candidate selection。
- p6: ask、outer-pair probabilities と inner-pair probabilities。
- p7: binding/exposure の言い回し、complete-command、時間・メモリの手法比較。
- p8: 本文 markup はない。

今回の対象は p4 の文脈の図だけ。他の本文修正は本文担当に引き渡す。

This figure addresses the p4 request for a diagram of the seven contexts.
The other text comments remain with the manuscript editor.

## Validation / 検証

寸法・埋め込み文字サイズ・文字の重なり・ページ外へのはみ出しを作図時に確認する。
最終 PNG の目視結果、B/C の画素比較、export hash は validation.json に追記する。
export 先は論文 repo の figures/figure1.pdf。manifest は Figure 1 項目だけを更新する。
git commit は行わない。

The renderer checks dimensions, embedded font sizes, text intersections and page bounds.
Final PNG inspection, B/C pixel comparisons and export hashes are recorded in validation.json.
The export is paper figures/figure1.pdf; only the Figure 1 manifest entry is changed.
No Git commit is made.
