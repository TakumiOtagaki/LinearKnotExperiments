**Figure 1. Conditional ensemble analysis with LinearKnot.**

Aでは、与えたscaffold Gに、対合塩基を共有しないextension G′を加えてRを作る。
GとG′はそれぞれ交差対を持たない。
灰色と橙色の下線は、互いに交差する入れ子対の二つの群（band）の範囲を示す。
縦の点線は二つの範囲を通る。

**(A)** Adding extension G′ to the supplied scaffold G gives structure R.
G and G′ share no paired bases, and each has no crossing pairs.
The grey and orange underlines show the spans of two groups of nested pairs, called bands, that cross each other.
The dotted vertical line passes through both spans.

BではCPartyの分解を二つずつの小さい組合せに分け、右端ごとに上位b候補を保持する（模式図はb=3）。
丸は保持、×は破棄を示す。
同じ分解でBoltzmann重みの和を取ると分配関数、エネルギーの最小を取ると最小自由エネルギー（MFE）構造を得る。

**(B)** LinearKnot divides CParty's decomposition into smaller combinations of two parts and keeps the top b candidates at each right endpoint (b=3 in this sketch).
Circles mark retained candidates; crosses mark discarded candidates.
The same decomposition gives the partition function by summing Boltzmann weights and a minimum-free-energy (MFE) structure by taking the minimum energy.

CはbpRNA_RFAM_31397（61 nt）のexact計算を示す。
Gは参照の入れ子のcanonical対（AU/UA、CG/GC、GU/UG）の8対である。
左上の対確率図のみ、確率0.05以上を表示し、弧の不透明度を対確率とする（数値キー参照）。
他の構造図の弧は不透明である。
灰色実線はG、青点線はGと交差しない追加対、橙破線はGと交差する追加対を示す。

**(C)** Exact calculations use bpRNA_RFAM_31397 (61 nt), with eight canonical nested reference pairs (AU/UA, CG/GC or GU/UG) as G.
Only the pair-probability plot at upper left uses arc opacity to show probability, as illustrated by its numeric key; it displays pairs with probability at least 0.05.
Arcs in the other structure diagrams have full opacity.
Solid grey arcs show G, dotted blue arcs show new pairs that do not cross G, and dashed orange arcs show new pairs crossing G.

各位置の七つのcontext確率の和は1であり、PseudoLoopはpseudoknotに伴う非対合loopを示す（定義はMethods）。
accessibilityは窓内の全塩基が同時に非対合である確率であり、実線・破線・点線は1・4・8 ntの窓を示す。
座標は1始まりで、accessibilityでは窓の開始位置を表す。
右上は10,000回のBoltzmann標本化で最多と2番目の構造（3,431回と1,426回）を示す。
残る5,143回の標本は示さない。
右下はγ=1のγ-centroid構造である。

The seven context probabilities sum to one at each position.
PseudoLoop denotes unpaired loops associated with a pseudoknot, as defined in Methods.
Accessibility is the probability that all bases in a window are unpaired together; solid, dashed and dotted lines show 1-, 4- and 8-nt windows.
Coordinates are 1-based positions or window starts.
The upper-right diagrams show the two most frequent structures in 10,000 Boltzmann draws (3,431 and 1,426 draws); the other 5,143 draws are not shown.
The lower-right diagram shows the γ-centroid structure at γ=1.
