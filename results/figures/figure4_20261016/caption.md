<!-- 日本語：同じ長い RNA の同じ区間で複数の予測量を読み合わせる。G は保持する固定対の集合である。 -->
**Figure 4. Reading predicted pairing, local structure, and window exposure together in NEAT1.**
The scaffold G is a set of base pairs kept in every predicted structure.
The threshold-0.7 scaffold uses LinearPartition pairs with probability at least 0.7 in runs of two or more stacked pairs.
Beam width limits the candidates retained during calculation.
Panels B and C show local views of full-transcript calculations.

<!-- 日本語：A は全長に対する表示区間。B は同じ表示区間内に両端をもつ対の確率。左右は G ありと空の条件。 -->
**(A)** Location of the interval in NEAT1_2 (22,743 nucleotides, nt).
The highlight marks 50 eight-base windows starting at 20,451–20,500.
**(B)** Base-pairing probabilities with G (left) and an empty scaffold (right), at beam width 100.
Gray marks pairs in G; blue marks new pairs without crossings of G.
Darker arcs have higher probabilities.
Arcs show only pairs with both bases in 20,350–20,600 and probability at least 0.01.

<!-- 日本語：C は同じ座標の七つの context と accessibility。網掛けは固定塩基を含む window の定義上のゼロを示す。 -->
**(C)** Seven structural-context probabilities (defined in Methods) and window accessibility, the probability that all eight bases are unpaired together.
Hatching marks windows containing a base paired in G; their accessibility is zero.
All 50 highlighted windows contain only bases unpaired in G.
Their means are 0.527 with G and 0.014 with an empty scaffold.
The highlighted bases extend to 20,507; window starts end at 20,500.

<!-- 日本語：D は別のアイソフォーム NEAT1_1 の全長標本。共通の G を除いて構造多様性を測る。 -->
**(D)** NEAT1_1 (3,756 nt): 1,000 samples at each beam width, compared with each other and with minimum-free-energy (MFE) and maximum expected accuracy (MEA) structures.
Distances exclude the 727 common pairs in G.
Jaccard distance is one minus the shared-pair count divided by the union-pair count.
Each distribution uses 499,500 distinct unordered sample pairs or 1,000 sample-to-reference distances.
Histograms have bin width 0.01 and area one; numbers give exact means.
MEA uses accuracy weight 1.
Finite-beam MFE may search a different structure set.
The panels compare predictions, rather than validate them experimentally.
