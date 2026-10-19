# Caption drafts / caption 案

Figure numbers match the adopted E1 manuscript figures.
各英語段落の直前に、対応する日本語段落を置く。

## Figure 2 — figure2_main.pdf

**Beam-width sensitivity of conditional-ensemble outputs.**
左列は400 nt未満の24,625 RNAをexactと比べ、右列は400–4,381 ntの1,438 RNAを有限幅b=400と比べる。
色と記号は全行でbeam幅bを示す。
幅400の自己比較は表示しない。

The left column compares 24,625 RNAs shorter than 400 nt with exact calculations.
The right column compares 1,438 RNAs of 400–4,381 nt with finite beam width b=400.
Colors and symbol shapes identify beam width b in every row.
The b=400 self-comparison is omitted because its differences are zero by definition.

Aは自由エネルギーの符号付き差ΔF=F_b−F_refを示す。
F_bとF_refは有限幅と参照のensemble free energyである。
細線は全範囲、太線は10–90パーセンタイル、記号は中央値を示す。
灰色帯の−0.001から0.001 kcal/molの間だけ線形とした対称対数軸を使う。

A shows the signed free-energy difference ΔF=F_b−F_ref, where F_b and F_ref are the ensemble free energies at finite beam width and in the reference calculation.
Thin lines show the full range, thick lines the 10th–90th percentiles, and symbols the median across RNAs.
The symmetric logarithmic axis is linear within the shaded band from −0.001 to 0.001 kcal/mol.

Bでは、各RNA内でGの非対合塩基についてtotal variation distance（TV）を平均した後、RNA間の90パーセンタイルを取る。
塗った記号は相手塩基と非対合状態の分布、白抜き記号は七つのcontextの分布を比べる。
TVは各分布の確率差の絶対値の和の半分である。

B first averages total variation distance (TV) over bases unpaired in scaffold G within each RNA, then takes the 90th percentile across RNAs in each length group.
Filled symbols compare distributions over pairing partners and the unpaired state; open symbols compare the seven structural-context probabilities.
TV is half the sum of absolute probability differences over the states in each distribution.

Cは、16 nt窓の全塩基が同時に非対合となる確率の平均絶対誤差を示す。
Gの対合塩基を含まない窓についてRNA内で平均し、長さ群ごとにRNA間の90パーセンタイルを取る。
評価可能な窓がないRNAは除き、確率ゼロの窓は含める。
対称対数軸は灰色帯の0–0.000001の間で線形であり、差に正の下限を代入していない。
他の窓長と評価件数は補足図とevaluation_counts.tsvに示す。

C shows the 90th percentile across RNAs of the mean absolute error (MAE) in 16-nt window accessibility, the probability that all bases in the window are unpaired together.
Within each RNA, the mean includes only windows with no base paired in G.
RNAs without an eligible window are excluded; eligible windows with zero probability are included.
The symmetric logarithmic axis is linear in the shaded band from 0 to 0.000001, and no positive floor is applied to the data.
The supplementary accessibility figure shows all four tested window lengths.

横軸は長さ群であり、群の中では設定を等間隔に並べる。
hとdは局所探索の上限であり、h=d=max(b,30)とした。
そのため、これはbeam幅と局所上限を合わせた設定の感度である。

The horizontal positions group RNAs by length and space the tested settings equally within each group.
The local search limits h and d are set to max(b,30), so the figure describes sensitivity to the combined settings.

## Figure S1 — figure2_accessibility_si.pdf

**Accessibility differences by window length and RNA length.**
A–Dは1/4/8/16 nt窓について、Figure 2Cと同じRNA内平均とRNA間90パーセンタイルを使う。
参照、色、記号、評価域、対称対数軸はFigure 2と同じである。
長さ群・窓長ごとの評価数と除外数はevaluation_counts.tsvに保存する。

A–D show 1-, 4-, 8- and 16-nt windows using the same within-RNA means and between-RNA 90th percentiles as Figure 2C.
References, colors, symbols, eligibility and symmetric logarithmic axes match Figure 2.
Evaluation counts and excluded RNAs for every length group and window length are provided in evaluation_counts.tsv.

## Figure S2 — figure2_d_si.pdf

**RNAs above a positive free-energy-difference threshold.**
A/BはFigure 2と同じexact/有限参照の群を使う。
ΔF=F_b−F_refが0.01 kcal/molを超えるRNA数を、各長さ群のRNA数で割った割合を示す。
この片側の条件は負の差を数えない。
幅400の自己比較はゼロになるため表示しない。
色と記号、群内の等間隔の設定はFigure 2と同じである。

A and B use the exact and finite references from Figure 2, respectively.
Each percentage divides the number of RNAs with ΔF=F_b−F_ref greater than 0.01 kcal/mol by all RNAs in that length group.
This one-sided criterion excludes negative differences.
The b=400 self-comparison is omitted because its percentage is zero by definition.
Colors, symbols and equally spaced settings within each length group match Figure 2.
RNA counts by length group are 50–99 nt: 12,428; 100–199 nt: 9,656; 200–299 nt: 1,446; 300–399 nt: 1,095; 400–799 nt: 584; 800–1599 nt: 423; 1600–3199 nt: 340; 3200–4381 nt: 91.
