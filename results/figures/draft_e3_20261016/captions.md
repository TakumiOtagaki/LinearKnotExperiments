# E3 figure draft captions / 図の caption 案

All panels show predictions from saved finite-beam calculations, not experimentally determined structures.
All RNA and window coordinates below are one-based.
Color is supplemented by labels or line styles; the palette follows Okabe-Ito colors.
Scaffold constraints are imposed exactly after checking raw numerical residuals: q=0 and Stem=1 at paired scaffold bases.
For LinearPartition, q is clipped to [0,1] after subtracting the exported incident BPP sum from one; raw-sum sensitivity is recorded in the analysis tables.
全図は保存済みの有限 beam 計算による予測を示す。
RNA と window の座標は1起点とする。
色に加えてラベルまたは線種で条件を識別する。
生値の残差を検査した後、scaffold の固定塩基は q=0、Stem=1 の厳密な条件を代入した。
LinearPartition の q は保存 BPP の incident sum を1から引いて [0,1] に補正し、生の和に対する感度は集計結果に残した。

## a_shape_profiles_bins

JA：NEAT1_1 の SHAPE 反応性と非対合確率を同じ座標に示す。
上段は SHAPE の生値と中心51 nt の移動平均、中段は b=200 の LinearPartition、空の G、G=0.7 の非対合確率の移動平均であり、G=0.7 の生値を薄く重ねた。
下段は五条件の非対合確率を幅0.1の区間に分けた SHAPE 平均であり、空の区間は表示しない。
SHAPE の位置1を転写産物の nt 1と仮定し、測定のある3,595塩基を用いた。
LinearPartition の非対合確率は閾値1e-4の BPP 出力から求め、LinearCParty は窓長1の Raccess を用いた。

EN: SHAPE reactivity and predicted unpaired probability along NEAT1_1.
The top panel shows raw SHAPE values and a centered 51-nt moving mean.
The middle panel shows moving means of unpaired probabilities for LinearPartition, empty G, and G=0.7 at b=200, with raw G=0.7 values in the background.
The bottom panel shows mean SHAPE reactivity in probability bins of width 0.1 for five conditions; empty bins are omitted.
We assume that SHAPE position 1 corresponds to transcript nucleotide 1 and use all 3,595 measured positions.
Unpaired probabilities come from length-1 Raccess for LinearCParty and exported BPPs with cutoff 1e-4 for LinearPartition.

## b_shape_metrics

JA：NEAT1_1 の SHAPE に対する非対合確率の Spearman 相関と AUROC を、scaffold と beam 幅ごとに比べる。
左列は3,595塩基、右列はそのうち G=0.7 で自由な共通の位置集合を示す。
AUROC の正例は全測定値の75百分位以上の SHAPE とし、全条件で同じ閾値を用いた。
線は測定点を結ぶもので、幅に対する単調性を仮定しない。

EN: Spearman correlation and AUROC of unpaired probability against NEAT1_1 SHAPE reactivity across scaffolds and beam widths.
The left column uses all 3,595 measured positions; the right column uses the common subset free in G=0.7.
Positive labels for AUROC are SHAPE values at or above the global 75th percentile, with the same threshold in every condition.
Lines connect evaluated settings without assuming monotonic beam-width sensitivity.

## c1_bpp_arcs_neat1_1

JA：G=0.7 の NEAT1_1 全体の BPP を b=100、200、400 で示す。
灰色は scaffold の対、青は G と交差しない新しい対、朱色は G の少なくとも一対と交差する新しい対であり、濃さを確率に対応させた。
読みやすさのため確率0.01以上の対だけを描いた。
BPP の保存形式は i、j、確率の3列であり、交差は座標から判定したもので square/round の出力列ではない。

EN: BPP arcs across NEAT1_1 for G=0.7 at b=100, 200, and 400.
Gray arcs denote scaffold pairs, blue arcs denote new pairs that do not cross G, and vermilion arcs denote new pairs crossing at least one scaffold pair.
Opacity increases with probability; only pairs with probability at least 0.01 are drawn.
The saved BPP format contains only i, j, and probability, so crossing is classified geometrically rather than from square/round output columns.

## c2_bpp_arcs_long_range

JA：NEAT1_2 の nt 282–546 または600–840と nt 20,761–21,120 の間にある保存 BPP を示す。
二つの表示領域の間の座標を圧縮し、切れ目を記した。
線がない条件は記載の閾値以上の対が出力されなかったことを意味し、真の確率が0であるとはしない。
未完了の条件は欠測と表示した。
この図は報告された区間での計算結果を示し、実験で確認された個別の塩基対を示さない。

EN: Exported NEAT1_2 BPPs joining nt 282–546 or 600–840 to nt 20,761–21,120.
The intervening sequence is compressed and the coordinate gap is marked.
An empty panel means that no pair reached the stated export cutoff, which does not establish zero pairing probability.
Unfinished conditions are labeled missing.
These panels show predictions in the specified intervals, not experimentally established individual base pairs.

## d_context_profiles

JA：NEAT1_1 の nt 3020–3520 における七つの context 確率、SHAPE、非対合確率、窓長8の accessibility を対応づける。
七つの context は各位置で和が1となる。
表示区間は b=200 の G=0.7 と空の G の context 分布の total variation の501 nt 平均が最大となる区間として探索的に選んだ。
SHAPE の黒線は中心21 nt の移動平均である。
この区間選択は予測を見やすくするためのもので、生物学的な domain の指定ではない。

EN: Seven structural-context probabilities, SHAPE reactivity, unpaired probability, and 8-nt window accessibility at NEAT1_1 nt 3020–3520.
The seven context probabilities sum to one at each position.
This exploratory interval maximizes the 501-nt mean total variation between the context distributions for G=0.7 and empty G at b=200.
The black SHAPE curve is a centered 21-nt moving mean.
The interval is selected to display model output and is not assigned a biological domain.

## e_accessibility_differences

JA：b=100 における G=0.7 と空の G の window accessibility の差を、窓長1、4、8、16で示す。
Accessibility は window の全塩基が同時に非対合となる確率である。
灰色は生の差、青は NEAT1_1 で中心51 nt、NEAT1_2 で中心101 nt の移動平均を示す。
橙色の帯は指定された NEAT1_2 の長距離相互作用の区間である。
G の変更は条件とする scaffold と ensemble を同時に変えるため、この差を pseudoknot だけの効果とは解釈しない。

EN: Differences in window accessibility between G=0.7 and empty G at b=100 for window lengths 1, 4, 8, and 16.
Accessibility is the joint probability that all bases in a window are unpaired.
Gray curves show raw differences; blue curves show centered moving means of 51 nt for NEAT1_1 and 101 nt for NEAT1_2.
Orange bands mark the specified NEAT1_2 long-range intervals.
Changing G changes both the scaffold condition and the allowed ensemble, so this contrast does not isolate the effect of pseudoknots alone.

## f_sampling_diversity

JA：G=0.7 の NEAT1_1 から b=100 と200で各1,000構造を抽出した。
左と中央は全499,500標本対についての塩基対集合の Jaccard 距離であり、完成構造と G を除いた extension を分けて示す。
右は標本と MFE・MEA 構造の extension の距離の中央値および5–95百分位を示す。
乱数 seed は20261016、MEA の γ は1である。
MFE の有限 beam 探索集合は分配関数の保持 ensemble と一致するとは限らない。

EN: Diversity of 1,000 NEAT1_1 structures sampled with G=0.7 at each of b=100 and b=200.
The left and middle panels show pair-set Jaccard distances for all 499,500 unordered sample pairs, using completed structures and extensions after removing G, respectively.
The right panel shows the median and 5th–95th percentiles of extension distances from samples to the MFE and MEA structures.
The sampling seed is 20261016 and MEA uses gamma=1.
The finite-beam MFE search set need not equal the retained partition-function ensemble.
