# Compact manuscript Figures 2 and 3

本文用の図は183 mm幅で、Figure 2は92 mm高、Figure 3は55 mm高とする。
Figure 2のDは55 mm高の独立したSI図に移す。
文字は原寸で7 pt以上とし、凡例はプロットの外に置く。

The manuscript figures are 183 mm wide: Figure 2 is 92 mm high and Figure 3 is 55 mm high.
The separate Figure 2 D for the Supplementary Information is 55 mm high.
All text is at least 7 pt at this width; legends sit outside the plots.

## Figure 2

左は400 nt未満のexact参照、右は400 nt以上の幅400参照を使う。
Aは符号付き自由エネルギー差、BはBPPとloop profileのTV、Cはaccessibilityの平均絶対誤差を示す。
分布・分位点・評価域は既存の図と同じである。
色分けの対象はA/Cではbeam幅、Bでは長さ群である。

Left panels use exact references below 400 nt; right panels use width 400 for longer RNAs.
A shows signed free-energy differences; B shows total variation distances for BPPs and loop profiles; C shows mean absolute accessibility differences.
The distributions, quantiles, and evaluation domains match the existing figure.
Colors identify beam widths in A/C and length groups in B.

## Figure 3

A/Bは実配列の時間とpeak RSS、Cの2枚はLinearCPartyの合成対照の時間とpeak RSSを示す。
集計はfebd2bc1を使い、既存の図と同じ中央値と全範囲を描く。

A/B show time and peak RSS for real sequences; the two plots in C show time and peak RSS for the LinearCParty synthetic controls.
The source revision is febd2bc1; medians and full ranges match the existing figure.

## Supplementary Figure: fraction above the free-energy threshold

符号付き差ΔFが0.01 kcal/molを超えるRNAの割合を示す。
左は400 nt未満のexact参照、右は400 nt以上の幅400参照を使い、色と記号は長さ群を示す。
閾値は表示のための選択である。

Fraction of RNAs with a signed free-energy difference ΔF above 0.01 kcal/mol.
The reference is exact below 400 nt and width 400 for longer inputs; colors and symbols identify length groups.
The threshold is a display choice.

## Sources

| Input | Path relative to this directory | SHA-256 |
|---|---|---|
| e1_short_metrics.tsv | `../../local/e1-bprna-short-d21618c5-20261016/metrics.tsv` | dbc68f96c9ebb670b49e4bf091e7a3a818fb6057b8a5fd84d944f38c6614078a |
| e1_short_quantiles.tsv | `../../local/e1-bprna-short-d21618c5-20261016/quantiles.tsv` | d19d2ac65d3ac140203de9e0f96e044e7afcf76e115e26328f870f02488e45d8 |
| e1_long_metrics.tsv | `../../local/e1-bprna-long-d21618c5-20261016/metrics.tsv` | 5201ba0fe4a543d27a92fa809fc717810324f5a1bc1c4fafb200339d6c7c2751 |
| e1_long_quantiles.tsv | `../../local/e1-bprna-long-d21618c5-20261016/quantiles.tsv` | 10e9bf08bbd0ca3c3295a09a48142778dcaea218409966731f8e084c8dd0d367 |
| e2_metrics.tsv | `../../local/e2-febd2bc1-20261016/metrics.tsv` | 6c7646276892b3f81933416afb1b0f8afd04ef0d1688e16c47a283f43278cfd1 |
| e2_comparisons.tsv | `../../local/e2-febd2bc1-20261016/comparisons.tsv` | 8ff532359967045c583697f3759595a1f799a442c53d5d162a6154ff97dda00a |
| e2_memory_metrics.tsv | `../../local/e2-febd2bc1-20261016/memory_metrics.tsv` | 87334e9c6cf0560213cca4f61465d362b872d491dd10b0b31a9001e6c8329201 |
| e2_summary.json | `../../local/e2-febd2bc1-20261016/summary.json` | 5059e8e1b6b4010805b0227e1d4a92a4a38fb0ece2535482257cc61b53250ad7 |
| g_features | `../../local/e2-e2cecbb4-g-features-20261016/rnacentral_g_features.tsv` | 856c1f56b79ef524038d1885540203ecdf9f037b5946592f0179e12032978c82 |

## Regeneration

From the paper repository (the default layouts and existing output directories are unchanged):

```bash
PYTHONDONTWRITEBYTECODE=1 python3 experiments/scripts/make_figures.py --main-compact --e2 experiments/results/local/e2-febd2bc1-20261016 --e2-revision febd2bc1 --out experiments/results/figures/main_20261016
```
