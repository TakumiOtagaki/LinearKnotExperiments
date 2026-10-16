<!-- 日本語：NEAT1の主な予測をbeam幅200で示す。Gは全構造で保持する対。beam幅は計算中に保持する候補を制限する。AはNEAT1_2の21,000–21,300、BはGあり／空の対合確率で区間内の0.01以上の対を示す。 -->
**Figure 4. NEAT1 predictions at beam width 200.**
Scaffold G contains pairs kept in every structure.
Beam width limits retained candidates.
**(A)** NEAT1_2 interval 21,000–21,300.
**(B)** Pairing probabilities with G or empty G; arcs require probability ≥0.01.

<!-- 日本語：CのLinearPKCapRは各塩基の七つの局所構造の確率、LinearPKRaccessは8塩基が同時に非対合の確率を示す。注目窓の開始は21,134–21,143で全塩基がGで自由。b=200の平均はGあり0.7851／空0.0110。b=200と比べ、各注目窓の変化はGありのb=400（薄青）、空のGのb=100でいずれも0.05未満。網掛けはGの対合塩基を含む窓。 -->
**(C)** LinearPKCapR gives seven structural-context probabilities per base.
LinearPKRaccess gives accessibility: the probability that eight window bases are unpaired together.
Highlighted windows start at 21,134–21,143; all bases are unpaired in G.
Their means are 0.7851 with G and 0.0110 with empty G.
Each highlighted value changes by <0.05 at b=400 with G (pale blue) or b=100 with empty G, relative to b=200.
Hatching marks windows containing scaffold-paired bases.

<!-- 日本語：DはNEAT1_1の1,000標本。Gを除き、標本間およびMFE・MEA構造との距離を示す。Jaccard距離は1から共有対数／和集合対数を引く。 -->
**(D)** NEAT1_1: 1,000 samples, compared pairwise and with minimum-free-energy (MFE) or maximum expected accuracy (MEA) structures.
Distances exclude G.
Jaccard distance is one minus shared pairs divided by all distinct pairs.
