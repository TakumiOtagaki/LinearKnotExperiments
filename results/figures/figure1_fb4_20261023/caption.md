**Figure 1. Structure prediction and conditional ensemble analysis with LinearKnot.**

Aでは、与えたscaffold Gにextension G′を加えて完成構造Rを作る。
GとG′は対合塩基を共有せず、それぞれは交差対を持たない。
灰色と橙色の下線は、互いに交差する入れ子対の群（band）の範囲を示す。
density-2の制限では、交差でつながったband群のうち一つの位置を覆うbandは最大二つである。

**(A)** Adding extension G′ to the supplied scaffold G gives the completed structure R.
G and G′ share no paired bases, and each has no crossing pairs.
The grey and orange underlines mark the spans of two groups of nested pairs, called bands, that cross each other.
The density-2 restriction allows at most two bands in a connected crossing group to span one position.

Bは、二つのscaffold対の間に左右の枝を一つずつ加えたmultiloopを示す。
multiloopは、枝と呼ぶ閉じた部分構造をつなぐloopである。
各枠は状態、すなわち部分構造の集まりを表し、弧はその一例を示す。
黒い点で合流する線は、子を親につなぐ分解の一段階を表す。
CPartyは左の枝・内側のscaffold部分・右の枝を一度につなぐ。
LinearKnotは内側と右の枝を途中の状態にまとめ、次に左の枝を加える。
Eは局所エネルギーを表す。
最初の段階では0、次の段階でEを加える。
beam pruningは、部分構造の最後の塩基である右端ごとに同じ種類の状態を順位づけ、上位b個を残す（模式図はb=3）。
途中の状態も枝刈りの対象となる。
丸は保持、×は破棄を示す。
Cは枝刈りを行わない計算例である。

**(B)** The example has one new branch on each side of a multiloop between two scaffold pairs.
A multiloop connects closed substructures, called branches.
Each box is a state: a set of partial structures.
The arcs show one example from that set.
Lines that meet at a black dot form one decomposition step, joining the children to their parent.
CParty joins the left branch, inner scaffold part and right branch at once.
LinearKnot first joins the inner part and right branch into an intermediate state, then adds the left branch.
E denotes the local energy term.
The first step adds zero; the second adds E.
Beam pruning ranks states of the same type at each right endpoint, the last base in the partial structure, and keeps the top b (b=3 in this sketch).
It also prunes intermediate states.
Circles remain; crosses are discarded.
Panel C uses no pruning.

CはbpRNA_RFAM_31397（61 nt）のexact計算を示す。
Gは、pseudoknotを持つ参照構造の入れ子部分にある8対のAU、GC、GU対からなる。
左上の対確率図だけは、0.05以上の確率を弧の濃さで示す。
右の二つの標本図とγ-centroid図は、それぞれ一つの完成構造を示す。
灰色実線はG、青点線はGと交差しない追加対、橙破線はGと交差する追加対を示す。

**(C)** Panel C shows exact calculations for bpRNA_RFAM_31397 (61 nt).
G contains eight AU, GC, or GU pairs from the nested part of the reference structure, which also contains a pseudoknot.
Only the upper-left pair-probability diagram uses arc opacity to show probability (pairs with probability at least 0.05).
The sample diagrams and the γ-centroid diagram each show one completed structure.
Solid grey arcs show G; dotted blue arcs show new pairs that do not cross G; dashed orange arcs show new pairs crossing G.

七つのcontext確率の和は各塩基で1となる。
PseudoLoopは、交差するbandの間の非対合の隙間を示す。
accessibilityは窓内の全塩基が同時に非対合となる確率である。
実線・破線・点線は1・4・8 ntの窓を示す。
accessibilityの横軸は窓の先頭の塩基の位置を1始まりで示す。
例えば、8 ntの曲線の位置20は塩基20–27の窓を表す。
標本図は10,000回のBoltzmann標本化で最多と2番目に多い構造（3,431回と1,426回）を示す。
残る5,143回の標本は図示しない。
代表構造はγ=1のγ-centroid構造である。

The seven context probabilities sum to one at each base.
PseudoLoop denotes unpaired gaps between crossing bands.
Accessibility is the probability that all bases in a window are unpaired together.
Solid, dashed, and dotted curves show 1-, 4-, and 8-nt windows.
The accessibility x-axis gives the 1-based position of the first base in each window.
For example, position 20 on the 8-nt curve represents bases 20–27.
The sample diagrams show the two most frequent structures in 10,000 Boltzmann draws (3,431 and 1,426 draws).
The other 5,143 draws are not shown.
The representative structure is the γ-centroid structure at γ=1.
