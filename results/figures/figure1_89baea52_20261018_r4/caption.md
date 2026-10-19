**Figure 1. Conditional ensemble analysis with LinearCParty.**

<!-- 日本語：Aのscaffold Gとextension G′は対合塩基を共有せず、それぞれは交差対を持たない。
和集合がRを作る。
括弧はband、すなわち外部の同じ対群と交差する入れ子対の極大群の範囲を示す。
点線は二つのbandの範囲を通る。
density-2制限では、交差で連結したband群のうち同じ位置を覆うbandは最大二つである。
残りのloopと境界の制限はMethodsに示す。 -->
**(A)** Scaffold G and extension G′ share no paired bases, and each has no crossing pairs.
Their union forms R.
The brackets mark bands: maximal groups of nested pairs that cross the same outside pairs.
The dotted line passes through both band spans.
The density-2 restriction allows at most two bands from any connected group of crossing bands to span one position.
Methods give the remaining loop and boundary restrictions.

<!-- 日本語：BはCPartyの分解を小さな組合せへ並べ替え、右端の順に候補を処理する。
丸は保持候補、×は破棄候補を示す。
高い点ほど順位が高い（図の例では右端ごとに3候補を保持）。
順位はinside weightの対数に予備計算の周囲の重みの対数の推定値を足して決める。
inside weightは部分構造のBoltzmann重みの和である。 -->
**(B)** LinearCParty rearranges CParty's decomposition into smaller combinations and processes candidates by right endpoint.
Circles mark retained candidates and crosses mark discarded candidates.
Higher points have higher ranks (3 retained per endpoint in this sketch).
Ranking adds the log inside weight to an estimate of surrounding log weight from a preliminary run.
Inside weight sums the Boltzmann weights of partial structures.

<!-- 日本語：Cの青枠は一つの保持ensembleから得る対確率、七つのcontext確率、窓accessibility、Boltzmann標本、最大期待正確度（MEA）復元の五出力を含む。
AとCはbpRNA_PDB_345（63 nt）を使い、canonicalな入れ子の参照11対をGとする。
Cはモデルのexact予測である。
灰色実線はG、青色点線は新しい非交差対、橙色破線はGと交差する新しい対を示す。
対の弧は確率0.01より大きいものを示し、高確率ほど濃くする。 -->
**(C)** The blue frame contains five outputs from one retained ensemble: pair probabilities, seven context probabilities, window accessibility, Boltzmann samples, and maximum expected accuracy (MEA) decoding.
A and C use bpRNA_PDB_345 (63 nt), with 11 canonical nested reference pairs as G.
C shows exact model predictions.
Solid grey, dotted blue, and dashed orange arcs mark G, new non-crossing pairs, and new pairs crossing G, respectively.
Pair arcs show probabilities greater than 0.01, with darker arcs for higher probabilities.

<!-- 日本語：context確率の和は1である。
Internalにはこの図とTercの例で共通の疎な斜線を付ける。
括弧は塩基39–51のPseudoLoopを示す。
PseudoLoopは交差bandをつなぐ、またはbandをまたぐ分岐内の非対合領域を表す。
accessibilityは窓の全塩基が同時に非対合となる確率である（1, 4, 8 nt、実線・破線・点線）。
座標は1始まりの塩基位置または窓の開始位置である。 -->
Context probabilities sum to one.
Sparse hatching identifies Internal in this figure and the Terc example.
The bracket labels PseudoLoop at bases 39–51.
PseudoLoop labels unpaired regions connecting crossing bands or in branches spanning a band.
Accessibility is the probability that all window bases are unpaired (1, 4, 8 nt; solid, dashed, dotted lines).
Coordinates give 1-based positions or window starts.

<!-- 日本語：標本図は最多の構造（975回）と、新しい非交差対を含む構造の中で最多のもの（6回）を示す。
全1000回のうち残り19回の標本は示さない。
MEAの対の重みはγ=1とする。
別の最小自由エネルギー（MFE）探索は同じ分解を使い、独自のbeam候補を持つ。
この模式図の構造もexactの結果を示す。 -->
The diagrams show the most frequent sampled structure (975 draws) and the most frequent structure containing a new non-crossing pair (6 draws).
The remaining 19 of 1,000 draws are not shown.
MEA uses pair weight γ = 1.
The separate minimum-free-energy (MFE) search uses the same decomposition with its own beam candidates.
Its diagram also shows an exact result.
