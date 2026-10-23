**Figure 1. Structure prediction and conditional ensemble analysis with LinearKnot.**

Aでは、与えたscaffold Gにextension G′を加えて完成構造Rを作る。
GとG′は対合塩基を共有せず、それぞれは交差対を持たない。
灰色と橙色の下線は、互いに交差する入れ子対の群（band）の範囲を示す。
density-2の制限では、交差でつながったband群のうち一つの位置を覆うbandは最大二つである。

**(A)** Adding extension G′ to the supplied scaffold G gives the completed structure R.
G and G′ share no paired bases, and each has no crossing pairs.
The grey and orange underlines mark the spans of two groups of nested pairs, called bands, that cross each other.
The density-2 restriction allows at most two bands in a connected crossing group to span one position.

Bはscaffoldの一つのmultiloopを示す。
multiloopとは、枝と呼ぶ複数の閉じた部分構造をつなぐloopである。
CPartyは内側のscaffoldの部分と左右の枝を一度につなぐ。
LinearKnotはまず内側と右側をつなぎ、破線の枠で示す途中の状態を作る。
次に左側を加えてmultiloopを閉じる。
beam pruningでは、右端ごとに同じ種類の候補を順位づけ、上位b個を残す（模式図はb=3）。
丸は保持候補、×は破棄候補を示す。
同じ分解でBoltzmann重みの和を取ると分配関数、エネルギーの最小を取ると最小自由エネルギー（MFE）構造を得る。

**(B)** The arcs show one multiloop in the scaffold.
A multiloop connects several closed substructures, called branches.
CParty joins the inner scaffold part and the branches on its left and right in one step.
LinearKnot first joins the inner part to the right branch, forming the intermediate state in the dashed box.
It then adds the left branch and closes the multiloop.
Beam pruning ranks candidates of the same type at each right endpoint and keeps the top b (b=3 in this sketch).
Circles mark retained candidates; crosses mark discarded candidates.
The same decomposition gives the partition function by summing Boltzmann weights and an MFE structure by minimizing energy.

CはbpRNA_RFAM_31397（61 nt）のexact計算を示す。
Gは、pseudoknotを持つ参照構造の入れ子部分にある8対のAU、GC、GU対からなる。
左上の対確率図だけは、0.05以上の確率を弧の濃さで示す（数値キー参照）。
右の二つの標本図とγ-centroid図は、それぞれ一つの完成構造を示す。
灰色実線はG、青点線はGと交差しない追加対、橙破線はGと交差する追加対を示す。

**(C)** Panel C shows exact calculations for bpRNA_RFAM_31397 (61 nt).
G contains eight AU, GC, or GU pairs from the nested part of the reference structure, which also contains a pseudoknot.
Only the upper-left pair-probability diagram uses arc opacity to show probability (numeric key; pairs with probability at least 0.05).
The sample diagrams and the γ-centroid diagram each show one completed structure.
Solid grey arcs show G; dotted blue arcs show new pairs that do not cross G; dashed orange arcs show new pairs crossing G.

七つのcontext確率の和は各塩基で1となる。
PseudoLoopは、Methodsで定義した、交差するbandの間の非対合の隙間を示す。
accessibilityは窓内の全塩基が同時に非対合となる確率である。
実線・破線・点線は1・4・8 ntの窓を示す。
accessibilityの横軸は窓の先頭の塩基の位置を1始まりで示す。
例えば、8 ntの曲線の位置20は塩基20–27の窓を表す。
標本図は10,000回のBoltzmann標本化で最多と2番目に多い構造（3,431回と1,426回）を示す。
残る5,143回の標本は図示しない。
代表構造はγ=1のγ-centroid構造である。

The seven context probabilities sum to one at each base.
PseudoLoop denotes unpaired gaps between crossing bands, as defined in Methods.
Accessibility is the probability that all bases in a window are unpaired together.
Solid, dashed, and dotted curves show 1-, 4-, and 8-nt windows.
The accessibility x-axis gives the 1-based position of the first base in each window.
For example, position 20 on the 8-nt curve represents bases 20–27.
The sample diagrams show the two most frequent structures in 10,000 Boltzmann draws (3,431 and 1,426 draws).
The other 5,143 draws are not shown.
The representative structure is the γ-centroid structure at γ=1.
