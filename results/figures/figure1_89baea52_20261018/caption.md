Conditional ensemble analysis with LinearCParty.
**(A)** Scaffold G and extension G′ share no paired bases, and each has no crossing pairs.
Their union forms R.
The brackets mark bands: maximal groups of nested pairs that cross the same outside pairs.
The dotted line passes through both band spans.
The density-2 restriction allows at most two bands from any connected group of crossing bands to span one position.
Methods give the remaining loop and boundary restrictions.
**(B)** LinearCParty rearranges CParty's decomposition into smaller combinations and processes candidates by right endpoint.
Circles mark retained candidates and crosses mark discarded candidates.
Higher points have higher ranks (3 retained per endpoint in this sketch).
Ranking adds the log inside weight to an estimate of surrounding log weight from a preliminary run.
Inside weight sums the Boltzmann weights of partial structures.
**(C)** The blue frame contains five outputs from one retained ensemble: pair probabilities, seven context probabilities, window accessibility, Boltzmann samples, and maximum expected accuracy (MEA) decoding.
A and C use bpRNA_PDB_345 (63 nt), with 11 canonical nested reference pairs as G.
C shows exact model predictions.
Solid grey, solid blue, and dashed orange arcs mark G, new non-crossing pairs, and new pairs crossing G, respectively.
Pair arcs show probabilities greater than 0.01, with darker arcs for higher probabilities.
Context probabilities sum to one.
The bracket labels PseudoLoop at bases 39–51.
PseudoLoop labels unpaired regions connecting crossing bands or in branches spanning a band.
Accessibility is the probability that all window bases are unpaired (1, 4, 8 nt; solid, dashed, dotted lines).
Coordinates give 1-based positions or window starts.
The samples show the most frequent structure and the most frequent structure with a new non-crossing pair (981 of 1,000 draws; 19 other draws).
MEA uses pair weight γ = 1.
The separate minimum-free-energy (MFE) search uses the same decomposition with its own beam candidates.
Its diagram also shows an exact result.
