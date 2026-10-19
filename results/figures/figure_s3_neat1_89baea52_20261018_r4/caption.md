<!-- 日本語：主な予測は beam 幅 b=200 を使う。
beam 幅は保持する候補を制限する。
scaffold G は全構造で保持する対である。 -->
**Figure S4. Predicted probabilities and structure diversity in NEAT1.**
The main predictions use beam width b=200.
Beam width limits retained candidates.
Scaffold G contains pairs retained in every structure.

<!-- 日本語：A は NEAT1_2 の表示区間を示す。
注目する10個の8塩基窓は21,134–21,143から始まり、塩基21,134–21,150を覆う。 -->
**(A)** NEAT1_2 (22,743 nucleotides, nt).
The 10 highlighted eight-base windows start at 21,134–21,143 and span bases 21,134–21,150.

<!-- 日本語：B の灰色実線は G、青色点線は G と交差しない新しい対、橙色破線はGと交差する新しい対を示す。
弧の不透明度は対確率であり、重なった弧は濃く見えることがある。
弧は確率0.01以上で両端が表示区間内にある対だけを示す。
区間外へ伸びる対は G ありで53、空の G で204対あり、省略する。
C はそれらの寄与も含む全長計算である。 -->
**(B)** Solid grey arcs show G, dotted blue arcs show new pairs that do not cross G, and dashed orange arcs show new pairs crossing G.
Arc opacity represents pair probability; overlapping arcs may appear darker.
Arcs require probability ≥0.01 and both endpoints in view.
B omits pairs with a base outside the view: 53 with G and 204 with empty G.
C includes their full-transcript contributions.

<!-- 日本語：C の各 context は同じ0–1の尺度で示す。
0.001未満というラベルは表示区間内のすべての保存値がこの閾値未満であることを示す。
accessibility は窓内8塩基が同時に非対合である確率である。
網掛けは G の対合塩基を含み accessibility が0の窓を示す。
注目窓はすべて G で非対合の塩基だけを含む。
図内の平均と beam 幅・局所制限を変えた最大差は注目10窓について示す。 -->
**(C)** Context tracks share a 0–1 scale (definitions in Methods).
Context values labelled <0.001 are below this threshold throughout the displayed interval.
Accessibility is the probability that all eight window bases are unpaired together.
Hatched windows contain G-paired bases and have zero accessibility.
All highlighted windows contain only bases unpaired in G.
Means and changes between settings refer to the 10 highlighted windows.

<!-- 日本語：左下の小表はGを使った注目10窓のjoint accessibilityを六つのloopの種類に分けた平均である。
合計は0.785、Hairpinの寄与は0.749、他の五成分の和は0.036である。 -->
The small table gives mean joint accessibility by six loop types for the 10 highlighted windows with G.
Mean total accessibility is 0.785; Hairpin contributes 0.749, and the other five types sum to 0.036.


<!-- 日本語：D は NEAT1_1 全長の1,000標本で、共通 G の727対を距離から除く。
Jaccard 距離は1から共有対数を和集合の対数で割った値を引く。
MFE は最小自由エネルギー、MEA は最大期待正確度を意味し、MEA の重みは1である。
標本間の距離は499,500個、各参照との距離は1,000個を使う。
幅0.01のヒストグラムは各々の面積を1とし、破線は正確な平均、括弧内は距離数を示す。
有限 beam の MFE は異なる構造集合を探索しうる。 -->
**(D)** NEAT1_1 (3,756 nt), 1,000 samples; distances exclude 727 shared G pairs.
Jaccard distance is one minus shared-pair count divided by union-pair count.
MFE denotes minimum free energy.
MEA denotes maximum expected accuracy and uses weight 1.
Sample–sample uses all 499,500 unordered pairs; each sample-to-reference comparison uses 1,000 distances.
Histograms use bin width 0.01 and are normalized separately to area one.
Dashed lines mark exact means.
Finite-beam MFE may search a different structure set.
