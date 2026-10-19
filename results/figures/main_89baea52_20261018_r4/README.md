# Time and memory: 89baea52

All current measurements use `results/local/e2-89baea52-20261018/`.
The renderer reads no E1 measurements in this invocation.
Two figures give manuscript Figure 4 and SI Figure S3.
`plot_data.tsv` gives plotted values and ratio numerators/denominators.
`manifest.json` gives source and script hashes; `validation.json` gives geometry and font checks.
Main width is 178 mm; SI width is 165 mm; PNGs are 300 dpi.
Each program requests the outputs specified in `caption.md`.
Measured growth over these lengths does not establish an asymptotic order.

Reproduce with `python3 -B scripts/render_non_e1_89baea52.py --suffix _r4 --published-only`.

## Round-four changes

Figure 4 remains 178 × 55 mm. The longest-input labels use short arrowheads ending at 29.1 kb, without crossing the 22.5-kb point. Time / RSS symbols in C use ±0.18 row offsets. Figure S3 retains its original geometry and values; its caption explicitly states base-2 length and base-10 cost axes.

## PNG inspection and final PDF checks — 2026-10-19

- Figure 4: `figure3_main.png` を実際に開いて確認。29.1 kbの短い矢印が22.5 kbの点を横切らない。12 helicesの丸と菱形を分けて確認できる。
- Figure S3: `figure3_synthetic_all_cli_si.png` を実際に開いて確認。六つの対数軸のパネルと七系列の凡例を確認した。変更はcaptionの両軸の説明であり、画像の配置・値は旧版と同じ。

Ghostscript `txtwrite`, `TextFormat=0` でPDF内の実際の文字サイズを検査した。
すべて最小7 ptで、6 pt未満の文字はない。PDFの掲載幅は本文178 mm、SI 165 mmで、各図の高さは旧版と同一。
PNGは300 dpi。原稿への挿入後のページ割り確認は原稿担当が行う。

## 第4回通読の図修正：原稿担当への報告

対象は「直すこと」の1、3、5、6、7、9、10、13、14、15の図の項目。
図、作図スクリプト、caption案、各図のREADMEと付随する生成検証ファイルを更新した。
LaTeX原稿・原稿の日本語ファイル・保存済み実験出力は編集していない。
新しいRNA計算やscheduler投入、外部実装の変更、git commitは行っていない。
この依頼の担当範囲に従い、作業記録は本READMEにまとめた。`docs/progress.md`への転記は原稿担当が行える。
図番号は本文Figure 3がTerc、Figure 4が資源使用量、SI Figure S3が合成系列、S4がNEAT1である。
旧ファイル名の`figure4.pdf`、`figure3_main.pdf`、`figure_s3.pdf`は互換性のため保持した。

### 新しい図とcaption案

パスはこのdirectoryからの相対パス。各caption.mdには対応する日本語と英語を置いた。

| 掲載図 | 新PDF | 新PNG | Caption案 | 掲載サイズ mm |
|---|---|---|---|---|
| Figure 1 | [figure1.pdf](../figure1_89baea52_20261018_r4/figure1.pdf) | [figure1.png](../figure1_89baea52_20261018_r4/figure1.png) | [caption.md](../figure1_89baea52_20261018_r4/caption.md) | 178 × 145 |
| Figure 2 | [figure2_main.pdf](../e1_89baea52_20261019_r4/figure2_main.pdf) | [figure2_main.png](../e1_89baea52_20261019_r4/figure2_main.png) | [caption.md](../e1_89baea52_20261019_r4/caption.md) | 178 × 92 |
| Figure 3 | [figure4.pdf](../figure4_terc_89baea52_20261018_r4/figure4.pdf) | [figure4.png](../figure4_terc_89baea52_20261018_r4/figure4.png) | [caption.md](../figure4_terc_89baea52_20261018_r4/caption.md) | 178 × 100 |
| Figure 4 | [figure3_main.pdf](figure3_main.pdf) | [figure3_main.png](figure3_main.png) | [caption.md](caption.md) | 178 × 55 |
| Figure S1 | [figure2_accessibility_si.pdf](../e1_89baea52_20261019_r4/figure2_accessibility_si.pdf) | [figure2_accessibility_si.png](../e1_89baea52_20261019_r4/figure2_accessibility_si.png) | [caption.md](../e1_89baea52_20261019_r4/caption.md) | 165 × 137 |
| Figure S2 | [figure2_d_si.pdf](../e1_89baea52_20261019_r4/figure2_d_si.pdf) | [figure2_d_si.png](../e1_89baea52_20261019_r4/figure2_d_si.png) | [caption.md](../e1_89baea52_20261019_r4/caption.md) | 165 × 55 |
| Figure S3 | [figure3_synthetic_all_cli_si.pdf](figure3_synthetic_all_cli_si.pdf) | [figure3_synthetic_all_cli_si.png](figure3_synthetic_all_cli_si.png) | [caption.md](caption.md) | 165 × 178 |
| Figure S4 | [figure_s3.pdf](../figure_s3_neat1_89baea52_20261018_r4/figure_s3.pdf) | [figure_s3.png](../figure_s3_neat1_89baea52_20261018_r4/figure_s3.png) | [caption.md](../figure_s3_neat1_89baea52_20261018_r4/caption.md) | 165 × 145 |

### 直したこと

| 指摘番号 | 修正内容 |
|---|---|
| 1 / F2-1 | Figure 2Bをsymlogへ変更。0–1e−6は線形、その上は対数。正の下限値をデータへ代入せず、塗り／白抜きの横方向の間隔も広げた。300–399 nt・b=200のp90 0.0014963 / 0.0017751はゼロから8.990 / 9.169 mm離れる。captionはp90とRNA内平均TVの最大値0.373 / 0.206を区別する。局所最大誤差はこの軸へ追加していない。 |
| 3 / G-1, G-2 | Figure 1 / 3 / S4でGを灰実線、新しい非交差対を青点線、Gと交差する新しい対を橙破線へ統一。七つのcontextの色は共通。Figure 1のInternalの領域と凡例にTercと同じ疎な`/`を付けた。 |
| 5 / S4-1 | Figure S4CにGを使った注目10窓の六成分平均を小表として追加。`byloop_regions.tsv`のNEAT1_2、b=200、w=8、starts 21134–21143を読み、成分和と保存total、既存total曲線の注目平均の一致を確認した。単一塩基contextから成分を推定していない。空のGには未保存の成分を描かない。 |
| 6 / F2-2, S1-1 | Figure 2Aには−0.001–0.001 kcal/mol、B/CとS1には0–1e−6の線形域のキーを追加。S1のcaptionは右列400–4,381 nt、有限参照b=400、局所上限h=d=max(b,30)を単独で説明する。 |
| 7 / F4-1 | Figure 4A/Bの29.1 kb注記を短い矢印へ変更。22.5 kbの点を横切る長い線を除いた。 |
| 9 / F1-1 | Figure 1CにMost frequent（975 draws）とWith a new non-crossing pair（6 draws）を記載。後者の追加対を短い青矢印で示し、captionは二つ目の選択理由と残り19 drawsを説明する。 |
| 10 / F3-1 | TercのAで外側・内側参照群の両腕へO/Iの指示線を追加。DにもO: Outer / I: Innerを付けた。 |
| 13 / S3-1 | Figure S3のcaptionにRNA長は底2、時間とpeak RSSは底10という両軸の説明を追加。画像の配置、目盛り、集計値は旧版と同一。 |
| 14 / S2-1, S4-2 | Figure S2の各群名直下へ保存済みのNを追加。S4のキーをOpacity = pair probability、微小contextの表示を<0.001 in this viewへ変更し、captionに弧の重なりを説明した。 |
| 15 / F2-3, F3-2, F4-2 | Figure 2Aの薄線・太線・記号を同じ位置へ重ねた分布キーとsymbol: medianの文言を採用。Terc Aの見出しをPair-probability overviewへ変更。Figure 4Cの時間/RSSを±0.18行ずらし、12 helicesの記号を離した。表に関する項目は原稿担当の範囲。 |

### Captionへの追加・差し替え案の要点

以下は変更箇所の短い英語案。定義・条件を含む全文は上表のcaption.mdを使う。
日本語はそれぞれのcaption.mdで対応させた。

- 共通：Solid grey arcs show G, dotted blue arcs show new pairs that do not cross G, and dashed orange arcs show new pairs crossing G.
- Figure 1：The diagrams show the most frequent sampled structure (975 draws) and the most frequent structure containing a new non-crossing pair (6 draws). The remaining 19 of 1,000 draws are not shown.
- Figure 2：B and C use a linear scale in the shaded band from 0 to 0.000001 and a logarithmic scale above it. For 300–399 nt at b=200, p90 is 0.0015 for pairing and 0.0018 for contexts.
- Figure 3：O and I point to both arms of the outer and inner reference groups. D lists all nine reference pairs and their probability ranges.
- Figure S1：Right panels use b=400 as a finite reference for RNAs of 400–4,381 nt. The local search limits h and d are set to max(b,30).
- Figure S2：N gives the number of RNAs in each length group, the denominator of the percentage.
- Figure S3：Both axes are logarithmic: RNA length uses base 2, and time and peak RSS use base 10.
- Figure S4：The small table gives mean joint accessibility by six loop types for the ten highlighted windows with G. Mean total accessibility is 0.785; Hairpin contributes 0.749, and the other five types sum to 0.036. Arc opacity represents pair probability; overlapping arcs may appear darker.

### 検証結果

全8 PNGを実際に開き、最終版を確認した。
Ghostscriptで全8 PDFの文字を検査し、最小はすべて7 pt（パネル記号9 pt）。
PDF media boxは本文178 mm、SI 165 mmで、全8図の高さは旧版と同一。
PNGは300 dpi。E1、Terc、NEAT1の文字領域検査、E1・Terc・資源図・S1–S4の文字重なり検査も通った。
Figure 1のInternal凡例は2 mm内でも斜線を読み取れるよう確認した。

E1では保存された1,750 quantile scalar、350 summaries、短RNA98,500行・長RNA8,628行を照合した。
`plot_data.tsv`、`evaluation_counts.tsv`は旧図版とバイト単位で同じ。
Tercの`reference_pairs.tsv`と`accessibility_region.tsv`、NEAT1の`focus_accessibility.tsv`も旧版とバイト単位で同じ。
資源図のaudit 1,244行はすべて旧版と一致し、未掲載の並列実行図にだけ対応する36行を新directoryには出していない。
E1/E2のsource・script・出力・派生ファイルのhash、NEAT1のsource table/inputとscriptのhashを照合した。
旧directory内の47ファイルは作業前のSHA-256とすべて一致した。

NEAT1の六成分はHairpin 0.7492128263072918、Bulge 0.02899050444391345、Internal 0.004971175847930974、Multiloop 0.001963088991727549、Exterior 0、PseudoLoop 2.5908891996426885e−23。
合計0.785137595591022。図は三桁または通常文字のe表記で示し、READMEとprovenanceには元の精度を残した。

作図前のinode使用量は21,544,482 / 28,000,000で、残り約645万。
共有領域の新出力を約65ファイルと見積もり、tarの展開・Matplotlibキャッシュ・検査用一時ファイルはnode-localの一時領域だけを使った。
旧図を消す作業は行っていない。

### 再生成

experiments repositoryから実行する。E1の出力directoryは存在しない状態で始める。
非E1は指定した_r4 directoryだけへ書き込む。新しいRNA計算は行わない。

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B scripts/render_e1_89baea52.py --out results/figures/e1_89baea52_20261019_r4
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B scripts/render_non_e1_89baea52.py --suffix _r4 --published-only
```

原稿担当が新パスとcaption案をLaTeXへ入れた後に、本文・SIの通常PDFと直前提示版との差分PDFでページ割りを確認する。
