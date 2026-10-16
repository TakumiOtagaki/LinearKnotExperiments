# Figure 2・3 の caption 案（2026-10-16）

JA：既存の集計 TSV だけから作成した下書きである。
原稿の図番号・caption は変更していない。
Figure 2 の D1/D2 と Figure 3 の注記あり/なしを著者が比較して選ぶ。

EN: These drafts use only existing aggregate TSV files.
The manuscript figure numbers and captions have not been edited.
The author can compare D1/D2 for Figure 2 and the annotated/plain versions of Figure 3.

## ファイルと体裁 / Files and layout

- `figure2_d1_183mm.pdf/png`：長さ別の ΔF>0.01 の割合を D に置く案。
- `figure2_d2_183mm.pdf/png`：参照ごとに RNA をまとめ、正負の差と絶対差の割合を D に置く案。
- `figure3_plain_183mm.pdf/png`：注記なし、D を含む案。
- `figure3_annotated_183mm.pdf/png`：11.8/13.7 kb の点と比・G 特徴量を注記した案。
- `figure3_synthetic_all_cli_183mm.pdf/png`：C の詳細。3 CLI を同じ軸範囲で比較する。
- 各ファイルに `_89mm` 版もある / Each file also has an `_89mm` version.

JA：183 mm 版は2段幅、89 mm 版は各軸を縦に並べた1段幅のレビュー用下書きである。
89 mm 版はページ高を超えるため、採用時にはパネルを別々に配置する。
PDF のページ幅は指定値を保ち、PNG は300 dpiとした。
本文・軸名は8.5 pt、目盛り・凡例・注記は8 pt、図の表題は9 ptである。
色は Okabe–Ito 系を用い、白黒では点の形、線種、塗りとハッチでも区別する。
PDF は埋め込み TrueType フォントを使う。

EN: The 183 mm versions use two-column width; the 89 mm versions stack axes vertically for review at one-column width.
The 89 mm drafts exceed a page in height, so their panels should be placed separately if selected.
PDFs retain the stated page width; PNGs use 300 dpi.
Axis labels use 8.5 pt, ticks/legends/annotations use 8 pt, and figure headings use 9 pt.
Okabe–Ito colors are supplemented by marker shapes, line styles, fill and hatching for grayscale reading.
PDFs embed TrueType fonts.

## Figure 2 caption

JA：**Figure 2. Beam 幅に対する条件付き ensemble の解析量の感度。**
400 nt 未満の24,625 RNAは exact を参照し、400–4,381 ntの1,438 RNAは有限幅 b=400を参照した。
左右のパネル（89 mm版では上下）を参照ごとに分け、長い区分の値を exact に対する誤差とは解釈しない。
有限幅の版は d21618c5、短い区分の exact は d9804df6で、Turner 2004、37°C、outer順位、h=d=max(b,30)を用いた。
幅と局所上限を同時に変える設定全体の感度であり、pruningだけの効果や各RNAにおける単調減少は仮定しない。

EN: **Figure 2. Sensitivity of conditional-ensemble outputs to beam width.**
The reference is exact for 24,625 RNAs shorter than 400 nt and finite beam b=400 for 1,438 RNAs of 400–4,381 nt.
The left and right panels (upper and lower axes at 89 mm) separate the references; differences for long RNAs are not errors relative to exact values.
Finite-beam calculations use d21618c5 and short-RNA exact references use d9804df6, with Turner 2004, 37°C, outer ranking and h=d=max(b,30).
The plots describe sensitivity to the combined settings, which vary beam width and local caps together; they do not isolate pruning or assume monotonic changes for individual RNAs.

JA：**A**：符号付きの ΔF=F_b−F_ref（kcal/mol）を長さ区分・幅ごとに示す。
細線は最小–最大、太線はp10–p90、箱はp25–p75、記号は中央値である。
縦軸は対称対数軸で、−0.001から+0.001の間だけ線形とし、0と負値を切り捨てない。
微小な差はこの中央部に圧縮される。
各RNAを等重みにした分位点を、昇順値 x_0,…,x_(N−1) の位置 h=(N−1)p、j=floor(h)に対する (1−h+j)x_j+(h−j)x_(j+1)（端点は最小・最大）で計算した。
既存quantiles.tsvの中央値・p90・p99・最小・最大と照合し、p10・p25・p75はmetrics.tsvから同じ線形補間で追加した。

EN: **A**, Signed ΔF=F_b−F_ref (kcal/mol) by length bin and beam width.
Thin lines show the minimum–maximum, thick lines p10–p90, boxes p25–p75 and symbols the median.
The symmetric logarithmic y-axis is linear only between −0.001 and +0.001; zeros and negative values are retained.
Very small differences are compressed in this central interval.
Quantiles give each RNA equal weight: for sorted values x_0,…,x_(N−1), h=(N−1)p and j=floor(h), the quantile is (1−h+j)x_j+(h−j)x_(j+1), with minimum/maximum at the endpoints.
Medians, p90, p99 and extrema were checked against quantiles.tsv; p10, p25 and p75 were added from metrics.tsv using the same linear interpolation.

JA：**B**：RNAごとの平均BPP TVと平均loop-profile TVのp90を長さ区分別に示す。
実線・塗った記号はBPP、破線・白抜き記号はloop profileを表す。
BPPはGで自由な塩基について相手塩基と非対合を含む分布のTV、loop profileは同じ塩基の七成分のTVを平均した値である。
TVはそれぞれの分布で 1/2 Σ|p_b−p_ref| とする。
平均TVは完成構造全体の分布のTVではない。

EN: **B**, p90 across RNAs of the mean BPP TV and mean loop-profile TV, separated by length bin.
Solid lines and filled symbols show BPP; dashed lines and open symbols show loop profiles.
BPP TV compares distributions over partner bases and the unpaired state at each base free in G; loop-profile TV compares the seven context probabilities at the same bases.
TV is 1/2 Σ|p_b−p_ref| for each distribution and is averaged over these bases within each RNA.
These mean marginal TVs do not measure TV over completed structures.

JA：**C**：窓長1、4、8、16 ntについて、RNAごとのaccessibility確率MAEのp90を示す。
各RNA内で、Gの対合塩基を含まないwindowの |p_b−p_ref| を平均した後、参照ごとに全長さ区分のRNAを等重みでまとめた。
長さ区分ごとのp90の中央値やwindowを混ぜたp90ではない。
評価可能なwindowを持たないRNAはNAとして分母から除き、0に置き換えない。
確率0のwindowは評価域内なら含める。
評価RNA数とwindow数は下表に示し、窓長間で評価域が異なる点に注意する。
Opening energyはこのパネルに含めない。

EN: **C**, p90 across RNAs of accessibility probability MAE for window lengths 1, 4, 8 and 16 nt.
Within each RNA, |p_b−p_ref| is averaged over windows containing no base paired in G; RNAs from all length bins are then pooled separately for each reference with equal RNA weights.
These are neither medians of bin-level p90 values nor quantiles pooled over windows.
RNAs without an eligible window are NA and excluded from the denominator, without replacement by zero.
Zero-probability windows are included when eligible.
The table below gives RNA and window counts; the evaluation domain differs among window lengths.
Opening energies are not included in this panel.

JA：**D1案**：各長さ区分・幅で、ΔF>0.01 kcal/molのRNA数を、その区分・幅でΔFを評価できる全RNA数で割り、100を掛けた割合を示す。
0.01は表示のための閾値であり、統計的有意性や全解析に共通の許容誤差を定義しない。

EN: **D1 option**, The percentage of RNAs with ΔF>0.01 kcal/mol in each length-bin/beam group: 100 times the number above the threshold divided by all RNAs with evaluable ΔF in that group.
The threshold is a display choice, not a significance threshold or a shared tolerance for all outputs.

JA：**D2案（D1の代替）**：参照ごとに全長さ区分のRNAをまとめ、ΔF>0.01、ΔF<−0.01、|ΔF|>1 kcal/molの割合を示す。
各割合の分母はshort全24,625 RNAまたはlong全1,438 RNAである。
曲線は重なる事象を含み、足して100%にはならない。
shortでは短いRNAの件数が多いため、そのRNAが集計を強く左右する。
b=400の自己比較の0は、exactとの一致を示さない。
E1の実行時間・RSSはこの図の全案で比較に使わない。

EN: **D2 option (alternative to D1)**, Percentages with ΔF>0.01, ΔF<−0.01 and |ΔF|>1 kcal/mol after pooling all length bins separately for each reference.
Each denominator is all 24,625 short RNAs or all 1,438 long RNAs.
The events overlap, so the curves do not sum to 100%.
The short-RNA summary is strongly influenced by the more numerous shorter RNAs.
Zero differences in the b=400 self-comparison do not establish agreement with exact values.
E1 wall time and RSS are excluded from every version of this figure.

## Figure 3 caption

JA：**Figure 3. 長いRNAと合成対照における時間とpeak RSS。**
版e2cecbb4、b=h=d=100、outer順位、Turner 2004、37°Cを用いた。
cpu_40割り当て内の同時実行数上限は10プロセス、各プロセス1 thread、RLIMIT_ASは8 GiBである。
LinearCPartyはMFEと分配関数・BPP、LinearPKCapRはloop profile、LinearPKRaccessはprofileと窓長1/4/8/16のaccessibilityを要求した。
Wall timeは起動・予備計算・本計算・復元・出力を含み、scaffold構築とbuildを含まない。
Peak RSSはプロセス全体の値で、GNU timeのKiBを2^20で割ったGiBを使う。
出力の異なるCLI全体の比較から、個別解析の追加コストは分離しない。

EN: **Figure 3. Time and peak RSS for long RNAs and synthetic controls.**
Calculations use e2cecbb4, b=h=d=100, outer ranking, Turner 2004 and 37°C.
Each cpu_40 allocation permits at most 10 concurrent processes, each with one thread and an 8 GiB RLIMIT_AS limit.
LinearCParty requests MFE and the partition function/BPP; LinearPKCapR requests loop profiles; LinearPKRaccess requests profiles and accessibility at window lengths 1/4/8/16.
Wall time includes startup, preliminary and main computation, traceback and output, but excludes scaffold construction and builds.
Peak RSS is for the entire process, in GiB obtained by dividing GNU time KiB by 2^20.
The CLIs request different outputs, so their total costs do not isolate the cost of an individual analysis.

JA：**A, B**：RNAcentralの7配列（4,793–29,106 nt）の長さとwall time（A）、peak RSS（B）を示す。
点はwarm-up1回を除く測定3回の中央値、縦線は最小–最大であり、信頼区間ではない。
3回の中央値は昇順の中央の測定値である。
ほぼ重なるCapR/Raccessの点も省略しない。
各長さは異なる配列Sとscaffold Gの組であるため、点を結ぶ回帰線や漸近次数は示さない。

EN: **A, B**, Sequence length against wall time (A) and peak RSS (B) for seven RNAcentral inputs of 4,793–29,106 nt.
Points show the median of three measured repetitions after one excluded warm-up; vertical lines show the minimum–maximum, not confidence intervals.
The median of three repetitions is the middle ordered measurement.
Nearly coincident CapR/Raccess points are retained.
Each length represents a different sequence S and scaffold G; no regression line or asymptotic order is shown.

JA：**C**：合成対照の主要7系列について、250/500/1,000/2,000/4,000 ntにおけるLinearCPartyの時間とRSSを両対数軸で示す。
各系列・長さの記号は、各seedの測定3回の中央値を求めた後、その3 seedの中央値を取った値である。
縦線は3 seed×3回の計9測定の最小–最大を示し、反復とseedの両方の幅を含む。
線は同じ系列の表示点を結ぶだけで、補間値や理論上界ではない。
全156入力のうち、5長さを共通に持つ7系列×5長さ×3 seed=105入力を図示した。
残る51入力は追加の分岐深さ対照であり、この図の対象には含めない。
非emptyのGはfloor(0.15n)対、GC比率は50%である。
形を変える際に配列Sも変わるため、形や特徴量の独立な因果効果とは解釈しない。
3 CLIを同じ軸範囲で示す詳細図は `figure3_synthetic_all_cli_*` にある。

EN: **C**, LinearCParty wall time and RSS for seven main synthetic series at 250/500/1,000/2,000/4,000 nt, on logarithmic axes.
For each series/length, symbols show the median across three seeds of the within-seed medians of three repetitions.
Vertical lines span the minimum–maximum of all nine measurements (three seeds × three repetitions), combining repeat and seed variation.
Lines only connect plotted values in the same series; they are not interpolated measurements or theoretical bounds.
The plots show 105 of the 156 synthetic inputs: seven series × five shared lengths × three seeds.
The remaining 51 inputs are additional branch-depth controls outside this plotted subset.
Nonempty G has floor(0.15n) pairs and GC content is 50%.
Sequence S also changes with scaffold shape, so these comparisons do not identify independent causal effects of shape or features.
The `figure3_synthetic_all_cli_*` detail compares all three CLIs on shared axis ranges.

JA：**D案**：同じ7 RNAに対して、inside（outer）、MFE、outside+BPP、loop profileの独立したプロセスのpeak RSSを示す。
各モードはwarm-upなしの1回測定であり、反復の範囲はない。
後三者は必要なforwardと準備も含む別計算で、insideへの追加分ではない。
Peakを積み上げたり差し引いたりしてbyte内訳にはしない。
予備計算の一部を省くinside_baselineと旧版の測定は図に含めない。
この案は本文のDまたはSupplementaryへの配置を選べる。

EN: **D candidate**, Peak RSS of independent processes for inside (outer), MFE, outside+BPP and loop profiles on the same seven RNAs.
Each mode has one measurement without warm-up and therefore has no repeat range.
The latter three modes include their required forward computation and preparation; they are separate computations rather than additions to an existing inside calculation.
Peaks are neither stacked nor subtracted to infer byte components.
The inside_baseline mode that omits part of the preliminary work, and measurements of the earlier version, are not plotted.
This panel can be placed in the main figure or the Supplementary material.

JA：**注記あり案**：11,832 nt（URS00023C9544）と13,652 nt（URS000241043C）のLinearCPartyの点を矢印で示す。
注記の時間比・RSS比の範囲は3 CLIの入力別中央値の比の最小–最大であり、反復範囲や信頼区間ではない。
Gで一つのloopに直接属する自由塩基数の最大と、同じloop内の自由端点対数J_Vを既存特徴量表から添えた。
J_Vの比は13,652 ntの値を11,832 ntの値で割ったものである。
異なるS/Gの観測比較であり、時間差の原因をこの図だけで確定しない。

EN: **Annotated option**, Arrows mark the LinearCParty points for 11,832 nt (URS00023C9544) and 13,652 nt (URS000241043C).
The annotated time/RSS ratio ranges span ratios of input medians across the three CLIs, not repetitions or confidence intervals.
The largest number of free bases directly belonging to one G loop and the number J_V of free endpoint pairs in the same loop come from the existing feature table.
The J_V ratio divides the 13,652 nt value by the 11,832 nt value.
These are observations for different S/G inputs; the figure alone does not establish the cause of the time difference.

## 出典 / Sources

| Source ID | Aggregate file |
|---|---|
| `e1_short_metrics.tsv` | [metrics.tsv](../../local/e1-bprna-short-d21618c5-20261016/metrics.tsv) |
| `e1_short_quantiles.tsv` | [quantiles.tsv](../../local/e1-bprna-short-d21618c5-20261016/quantiles.tsv) |
| `e1_long_metrics.tsv` | [metrics.tsv](../../local/e1-bprna-long-d21618c5-20261016/metrics.tsv) |
| `e1_long_quantiles.tsv` | [quantiles.tsv](../../local/e1-bprna-long-d21618c5-20261016/quantiles.tsv) |
| `e2_metrics.tsv` | [metrics.tsv](../../local/e2-e2cecbb4-20261016/metrics.tsv) |
| `e2_comparisons.tsv` | [comparisons.tsv](../../local/e2-e2cecbb4-20261016/comparisons.tsv) |
| `e2_memory_metrics.tsv` | [memory_metrics.tsv](../../local/e2-e2cecbb4-20261016/memory_metrics.tsv) |
| `g_features` | [rnacentral_g_features.tsv](../../local/e2-e2cecbb4-g-features-20261016/rnacentral_g_features.tsv) |

| Panel | Rows / columns and processing |
|---|---|
| 2A | E1 short/long `metrics.tsv`: `delta_free_energy_kcal_mol`, grouped by `tier,beam`; existing `quantiles.tsv` checked; p10/p25/p75 added by linear interpolation. |
| 2B | E1 short/long `quantiles.tsv`: `base_state_tv_mean_Gfree`, `context_tv_mean_Gfree`, p90 in each length bin/beam; recomputed from metrics. |
| 2C | E1 short/long `quantiles.tsv`: `bin=all`, `access_w{1,4,8,16}_probability_mae_eligible`, p90; RNA/window counts summed from metrics; no cross-reference pooling. |
| 2D1/2D2 | E1 short/long `metrics.tsv`: signed `delta_free_energy_kcal_mol`; direct strict-threshold counts and denominators, then ×100. |
| 3A/3B | E2 `metrics.tsv`, `run=rnacentral_current`, all 3 programs: `n`, `wall_seconds_{median,min,max}`, `rss_gib_{median,min,max}`; checked against the saved repetition arrays. |
| 3C/detail | E2 `metrics.tsv`, `run=synthetic_current`, seven named `series`, all five lengths; main C uses `program=LinearCParty`, detail uses all 3; median of three seed medians, extrema over all nine repeats. |
| 3D | E2 `memory_metrics.tsv`, `run=memory_current`, `mode=inside_outer/mfe_outer/outside_bpp_outer/loop_profile_outer`, `n,rss_gib`; independent one-shot measurements. |
| 3 annotation | E2 `comparisons.tsv`, `kind=rna_13652_over_11832`, three programs, `wall_seconds_ratio,rss_gib_ratio`; G features `loop_free_max,same_loop_pairs_turn` for URS00023C9544/URS000241043C. |

JA：元データのSHA-256は `manifest.json`、実際に描いた集約値・閾値件数は `plot_data.tsv` にある。
RNA-levelの分位点は元のquantiles.tsvと照合し、E2の中央値と範囲はmetrics.tsv内の反復配列と照合した。
検証件数と出力の物理サイズは `validation.json` に記録する。

EN: `manifest.json` records source SHA-256 hashes; `plot_data.tsv` records plotted summaries and threshold counts.
RNA-level quantiles were checked against the original quantiles.tsv, and E2 medians/ranges against repetition arrays in metrics.tsv.
`validation.json` records check counts and physical output sizes.

### Figure 2A/B/D1 RNA counts

| Reference | Length bin (nt) | RNA |
|---|---:|---:|
| <400 nt: exact reference | 50–99 | 12,428 |
| <400 nt: exact reference | 100–199 | 9,656 |
| <400 nt: exact reference | 200–299 | 1,446 |
| <400 nt: exact reference | 300–399 | 1,095 |
| 400+ nt: b=400 reference | 400–799 | 584 |
| 400+ nt: b=400 reference | 800–1599 | 423 |
| 400+ nt: b=400 reference | 1600–3199 | 340 |
| 400+ nt: b=400 reference | 3200–4381 | 91 |

### Figure 2C evaluation domain (same for every beam within each row)

| Cohort | Window (nt) | Evaluated RNA | NA RNA | Eligible windows |
|---|---:|---:|---:|---:|
| short | 1 | 24,625 | 0 | 1,684,864 |
| short | 4 | 24,592 | 33 | 1,088,368 |
| short | 8 | 22,322 | 2,303 | 707,799 |
| short | 16 | 14,060 | 10,565 | 418,465 |
| long | 1 | 1,438 | 0 | 1,176,310 |
| long | 4 | 1,438 | 0 | 790,413 |
| long | 8 | 1,438 | 0 | 551,798 |
| long | 16 | 1,426 | 12 | 360,700 |

JA：各参照・窓長のRNA数とwindow数が全幅で同じであることを検査した。
評価域内の参照正・有限幅0のwindow数も `plot_data.tsv` に保持する。

EN: RNA/window counts for each reference/window length were checked to be the same at every beam width.
`plot_data.tsv` also retains counts of eligible reference-positive/candidate-zero windows.

## 著者が選ぶ点 / Choices for the author

JA：Figure 2 D1は長さ別の尾部の残り方を読みやすく、A/Bの長さ区分と対応する。
D2は符号の逆転と大きな絶対差をまとめて示すが、全長さをまとめるため長さ依存を隠し、負側の割合は同じ0–100%軸では小さく見える。
推奨はD1を本文、D2をSupplementary候補とする。

EN: Figure 2 D1 shows length-dependent tails and matches the bins in A/B.
D2 summarizes negative differences and large absolute differences, but pooling hides length dependence and small negative fractions are less visible on the shared 0–100% axis.
The suggested choice is D1 for the main figure and D2 as a Supplementary candidate.

JA：Figure 3は注記あり案が11.8/13.7 kbの比較を本文と対応させやすい。
図の情報量を抑える場合は注記なし案を選び、G特徴量をSupplementaryに置く。
DはMFEの別計算のメモリを示せる一方、一回測定なので、本文に置く場合もA/Bの3回測定と区別する必要がある。
Cを3 CLIすべてに広げる場合は詳細図を使い、主図の高さ・凡例量との兼ね合いで決める。

EN: The annotated Figure 3 makes the 11.8/13.7 kb comparison easy to connect to the text.
The plain version reduces visual load and leaves G features to the Supplementary material.
D shows memory use of the separate MFE computation, but its one-shot measurements must remain distinct from the three-repetition summaries in A/B.
The all-CLI detail can replace/extend C if its extra height and legend space are acceptable.

## 再生成 / Regeneration

```bash
PYTHONDONTWRITEBYTECODE=1 python3 experiments/scripts/make_figures.py
```

JA：入力はすべてTSVで、出力先以外を変更しない。
`--e1-short`, `--e1-long`, `--e2`, `--g-features`, `--out` で場所を変更できる。
既存の入力ディレクトリを出力先として使うことは拒否する。
RNA CLI、scheduler、外部実装の実行・変更は行わない。

EN: All inputs are TSV files; only the output directory is written.
Locations can be changed with `--e1-short`, `--e1-long`, `--e2`, `--g-features` and `--out`.
An output directory overlapping input directories is rejected.
The script neither runs RNA CLIs/scheduler jobs nor executes or changes the external implementation.
