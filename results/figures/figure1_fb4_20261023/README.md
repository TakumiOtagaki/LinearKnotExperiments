# Figure 1 — states, partial structures and decomposition trees, 2026-10-23

[figure1.pdf](figure1.pdf): 178 × 145 mm.
[figure1.png](figure1.png): 300 dpi preview.
Minimum embedded text: 7 pt.
The bilingual caption proposal is below and in [caption.md](caption.md).

## Reproduce

From experiments/:

~~~sh
MPLCONFIGDIR=/tmp/figure1-fb4-mpl python3 -B scripts/make_figure1.py \
    --data results/local/pseudoloop-0.2.0-20261023/examples/figure_data.json \
    --out results/figures/figure1_fb4_20261023
~~~

Dependencies: Matplotlib, NumPy, Pillow, pypdf, Ghostscript.
Five output files; no RNA folding or experimental recomputation.
The input is the adopted LinearKnot 0.2.0 drawing data, including updated contexts.
Source, panel-C data and renderer hashes are in validation.json.

## Scientific drawing

B shows one multiloop between two scaffold pairs, with a new branch on each side.
Grey arcs are scaffold G; dotted blue arcs are added pairs that do not cross G.
The inner grey boundary delimits this band's calculation; its enclosed region is omitted.
Each outlined node is a state representing a set of partial structures.
Its thumbnail is one example, not the complete set.
Lines meeting at a small black dot are one edge with several children (a hyperedge).
The dot is an edge junction, not another state.
The arrows point from children toward the parent.
CParty has one three-child edge; LinearKnot has two two-child edges.
The intermediate contains the inner part and right branch, without the outer closing pair.
Its first edge adds zero local energy; the completion edge adds the original local term E once.
The two parent thumbnails and the larger lower-left multiloop have the same pairs.
The intermediate's dashed box points to a beam of intermediate states, for one state type.
At each right endpoint (last base of a partial structure), the sketch keeps 3 top-ranked states.
Black circles remain; grey crosses are discarded. These ranks are illustrative, not measured.
Coordinates, scaffold anchors, unpaired-base factors and admissibility details are omitted.

A and all C content keep their adopted geometry, labels and saved values.
The misleading B-to-C arrow is removed: C is an exact example without pruning.
C retains sample frequencies 3,431, 1,426, and 5,143 other draws out of 10,000.
No inference about these exact outputs is made from the illustrative beam grid.

## Implementation evidence (not figure labels)

The example is BE_WIP_WIP: left WIP × inner BE × right WIP.
Ordinary V/VM multiloop recursion is already binary in CParty; this drawing makes no general claim that CParty multiloops all have three children.
LinearKnot first combines inner BE × right WIP, then left WIP and the local outer factor.
The split coordinates and scaffold anchors are retained when grouping the products.
Contracting the intermediate restores the original edge and its weight before pruning.
The diagram does not claim that other local energy terms vanish.

Read-only evidence:
- results/figures/figure1_fb2_20261023/README.md.
- LinearKnot v0.2.0 documents/model.md, “Derivation correspondence” (BE wings).
- LinearKnot v0.2.0 src/scfg2/deductions/expand_basic.hh (BE_WIP_WIP).
- LinearKnot v0.2.0 src/scfg2/beam/be_frame_cache.hh: the intermediate is a product of the two child factors; the local factor enters on completion.
- LinearKnot src/scfg2/beam/be_right_tail_selection.hh for the separate intermediate selection.

## What was taken from the earlier versions

Viewed the earlier exported PDF at paper commit 867ccd7 and the adopted 0.2.0 PNG individually.
Read the earlier renderer history (experiments 2076dd8 / 350ba2f) and the current renderer.
The earlier boxes and branching lines made the decomposition tree clear; those are restored.
The current multiloop arcs, the inner-plus-right order, and the pair colors are retained.
Every tree node now has a structural thumbnail; the intermediate is explicitly a state.
A shared edge junction groups all three CParty children into one step.
Both binary steps and their local-energy assignments are visible in LinearKnot.

## Author feedback review

Opened m-1.png through m-9.png individually and checked the marked sentences against the PDF annotation coordinates.
Author StrikeOut annotations were distinguished from red latexdiff deletions.
The subtype-free Annots entries on pages 2, 4 and 5 resolve to XRef streams, not extra author comments.
Marks outside Figure 1 are recorded here for scope awareness and left to the manuscript editor.

- p1: highlight on the Lyngsø and Pedersen (2000) citation in the NP-hard sentence.
- p2: highlights on “No scaffold pair”, “may connect a base inside a completed pseudoknot region”, and “internal-loop length limits”; the margin asks what the undefined pseudoknot region means and asks for concrete limit values (the latter is visible in the image, not a typed annotation).
- p3: highlight on Intermediate and comment asking for the state–structure–tree relationship; comment contrasting the beam sketch with exact C; author StrikeOut on “numeric key;” and “with opaque arcs”. The new B addresses the tree comment, the B-to-C arrow is removed, and the caption proposal omits the struck phrases.
- p4: highlights on “calculation with a small beam” / “for ranking” and on the outside-calculation sentence; margin comments request the reason for the pilot beam and a short outside explanation with an example. An image-only red comment next to the TV formula asks for wording such as “TV is defined by”.
- p5: StrikeOut on the end of the context paragraph (other loops retaining labels, and the band multiloop example); margin suggests deleting from “other loops” or “For example”; highlights and comments on sampling “its parts”, computational-cost “dependence on the local limits” (image-only note: is “local limits” understandable?), and BPP “alternatives”.
- p6: StrikeOut on the sentence excluding RNAs without eligible windows; highlight and comment on “four window lengths”, asking to list them.
- p7: no author annotation or handwritten comment found.
- p8: highlight on the partition-function ratio in the accessibility explanation, with a comment questioning its placement and “no-crossing calculation”; StrikeOut and comment on the Discussion probability-output / centroid sentence, suggesting sampling or subopt instead; bottom comment asks for the default scaffold construction in Methods.
- p9: two adjacent StrikeOut annotations remove the statements that retaining every scaffold pair excludes replacement structures and that a crossing pair can have low probability when other structures have lower energies.

## Validation and handoff

Page size, embedded text size, text intersections and saved inputs are checked during rendering.
PNG visual inspection and unchanged A/C comparisons are recorded in validation.json.
Export: paper figures/figure1.pdf; only its export-manifest entry is updated.
LaTeX, bibliography and external implementation repositories are unchanged.
The canvas remains 178 × 145 mm, so the figure's manuscript footprint is unchanged.
No Git commit is made.

## Caption proposal (Japanese / English)

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
