<!-- 日本語：NEAT1 の予測を対応づける。G は全予測構造で保持する対の集合で、LinearPartition の確率0.7以上かつ連続2対以上の部分から作る。beam 幅は計算中に保持する候補を制限し、B–D の主曲線は 200。B・C は全長計算の局所表示である。 -->
**Figure 4. Predicted pairing, local structure, and window accessibility in NEAT1.**
The scaffold G is a set of base pairs kept in every predicted structure.
The threshold-0.7 scaffold uses LinearPartition pairs with probability at least 0.7 in runs of two or more stacked pairs.
Beam width limits the candidates retained during calculation.
The main results in panels B–D use beam width 200.
Panels B and C show local views of full-transcript calculations.

<!-- 日本語：A は全長22,743塩基に対する表示区間と開始21,134–21,143の10個の8塩基窓。B は両端が21,000–21,300内にある確率0.01以上の対を示し、左右は閾値0.7のGありと空の条件。灰色はGの対、青色はGと交差しない新しい対で、濃さは確率を表す。 -->
**(A)** Location of the interval in NEAT1_2 (22,743 nucleotides, nt).
The highlight marks 10 eight-base windows starting at 21,134–21,143.
**(B)** Base-pairing probabilities with G (left) and an empty scaffold (right).
Gray marks pairs in G; blue marks new pairs without crossings of G.
Darker arcs have higher probabilities.
Arcs show only pairs with both bases in 21,000–21,300 and probability at least 0.01.

<!-- 日本語：C の LinearPKCapR は各塩基がstemまたは六つの非対合contextに属する確率、LinearPKRaccess は窓内8塩基が同時に非対合である確率を示す。網掛けはGの対合塩基を含みaccessibilityがゼロの窓。注目10窓はすべてGで自由な塩基のみを含み、b=200の平均はGあり0.7851、空のGで0.0110。b=200と比べ、各注目窓の変化はGありのb=400（薄青）、空のGのb=100でいずれも0.05未満。覆う塩基は21,150まで、開始位置は21,143までである。 -->
**(C)** LinearPKCapR gives the probability that each base lies in a paired stem or one of six unpaired structural contexts (defined in Methods).
LinearPKRaccess gives window accessibility, the probability that all eight bases are unpaired together.
Hatching marks windows containing a base paired in G; their accessibility is zero.
All 10 highlighted windows contain only bases unpaired in G.
Their means at b=200 are 0.7851 with G and 0.0110 with an empty scaffold.
Each highlighted value changes by <0.05 at b=400 with G (pale blue) or b=100 with empty G, relative to b=200.
The highlighted bases extend to 21,150; window starts end at 21,143.

<!-- 日本語：D は別のアイソフォーム NEAT1_1 全長3,756塩基の1,000標本を、標本同士および最小自由エネルギーMFE・期待正確度最大化MEA構造と比較する。MEAは対合確率との期待一致度を最大化する。共通Gの727対を除き、Jaccard距離は1から共有対数／和集合対数を引く。標本間499,500距離、各参照への1,000距離を幅0.01・面積1の分布とし、数値は正確な平均を示す。MEA重みは1。有限beamのMFE探索集合は異なりうる。モデル予測の比較である。 -->
**(D)** NEAT1_1 (3,756 nt): 1,000 sampled structures, compared with each other and with minimum-free-energy (MFE) and maximum expected accuracy (MEA) structures.
MEA selects a structure that maximizes expected agreement with pairing probabilities.
Distances exclude the 727 common pairs in G.
Jaccard distance is one minus the shared-pair count divided by the union-pair count.
Each distribution uses 499,500 distinct unordered sample pairs or 1,000 sample-to-reference distances.
Histograms have bin width 0.01 and area one; numbers give exact means.
MEA uses accuracy weight 1.
Finite-beam MFE may search a different structure set.
This figure compares model predictions.
