# エネルギー独立照合の結果（2026-10-17）

**C(1,n) と片端に接する multiloop の closing assignment は、標準 ViennaRNA の dangles=2 と一致しなかった。**
SI 表 S3 は source と一致している。
表の転記の問題ではなく、固定版の局所スコアの境界条件による違いである。
意図した独自のモデルかどうかは source だけでは決められない。

依頼の Major 2・3 のうち、通常 loop の独立評価と境界の原因特定に対応した。
第2回 global review の C17、S6、T4、T10 にも根拠を追加する。
loop 別 accessibility、pseudoknot 固有のエネルギー、一般の完全性・非曖昧性を解決したとは扱わない。
原稿への反映は担当外なので行っていない。

## 全件の比較

ランダム150本、bpRNA 150本、境界を調べる設計例37本、計337本を比較した。
各入力について固定版と観測用対照を独立 ViennaRNA 2.7.2 と照合し、**674 / 674 条件が成功**した。
失敗・timeout は0だった。
全 i<j を比較し、各実装につき2,226,060座標を調べた。
表示 cutoff で小さな確率を捨てていない。

同じ Turner 2004 parameter file、37°C、dangles=2、betaScale=1、pf_smooth=1、最小 hairpin 3、内部 loop の非対合塩基合計上限30を用いた。
主比較の G は空である。
設定、入力の選び方、hash、丸めへの対処は [README](README.md) と [provenance](provenance.json) に示す。

「内側隣接対照」は node-local の一時 header で closing の二つの guard だけを除いた観測器である。
修正版の採用・実装変更を行ったという意味ではない。

| 入力 | n | 本数 | 重み | 最大絶対 log Z差 | 最大絶対 PF free-energy差 (kcal/mol) | 最大絶対 BPP差 |
|---|---:|---:|---|---:|---:|---:|
| ランダム | 20–200 | 150 | 固定版 | 8.65928e-3 | 5.33689e-3 | 4.26569e-3 |
| 同上 | 同上 | 150 | 内側隣接対照 | 1.44951e-12 | 8.95284e-13 | 1.24534e-12 |
| bpRNA | 53–199 | 150 | 固定版 | 1.16957e-2 | 7.20828e-3 | 5.70433e-3 |
| 同上 | 同上 | 150 | 内側隣接対照 | 1.23634e-12 | 7.60281e-13 | 1.22280e-12 |
| 設計例 | 5–27 | 37 | 固定版 | 1.31976e-3 | 8.13394e-4 | 1.25321e-3 |
| 同上 | 同上 | 37 | 内側隣接対照 | 7.10543e-15 | 4.44089e-15 | 7.21645e-15 |

出典：[summary.tsv](summary.tsv)、[ensemble_metrics.tsv](ensemble_metrics.tsv)。
全群の固定版の符号付き log Z 差は −0.01169566594518301 から +0.001030188719099101 だった。
上表は正側の最大値ではなく、絶対差の最大値を使う。

固定版で log Z と BPP の最大差を示したのは `bpRNA_RFAM_9283`（187 nt）だった。
log Z は固定版126.43769845235263、ViennaRNA 126.44939411829782だった。
対 (36,176) の BPP は、固定版0.4421584353983448、ViennaRNA 0.4478627693154214だった。
同じ配列の MFE は両者で一致した。
これは「MFE が同じなら ensemble も同じ」とは言えない例である。

## 固定構造のエネルギーと項別会計

548条件を評価し、未到達は0だった。
これには MFE、ViennaRNA の36標本から重複をまとめた33構造、設計例、dangles=0 の73条件、非空 G の48条件を含む。
同じ配列・構造を G の割り当てや dangles を変えて再評価した条件を、独立配列数として数えていない。
非空 G では、設計した24個の非交差 multiloop 構造について「全対が G」と「outer pair だけが G」を比較した。

ViennaRNA の全対を強制し、他の位置を非対合に固定した PF と、通常の `eval_structure` を両方使った。
全ての指定対の確率が1になることを確認し、その最大誤差は2.11e-15だった。
さらに、インストールされた ViennaRNA 自身の局所関数で全項を足し直した。
この項別計算では SCFG2 の局所重みを利用していない。

| 検査 | 最大絶対差 (kcal/mol) |
|---|---:|
| ViennaRNA 局所項の和 − 独立な一構造 PF | 5.68434e-14 |
| その和に closing guard の差を足した値 − 固定版 | 5.68434e-14 |
| ViennaRNA 局所項の和 − 内側隣接対照 | 5.68434e-14 |
| 固定版 − 独立な一構造 PF | **1.50** |
| 内側隣接対照 − 独立な一構造 PF | 2.84217e-14 |

45条件で guard による非ゼロのエネルギー差を検出し、全てをこの局所項の差で説明できた。
最大1.50の例は `ml_adjacent_GC_pad00` である。
内側に非対合塩基を置く25 ntの例 `ml_gapped_GC_pad00` では、固定版 −0.60、ViennaRNA −1.70だった。
差の1.10は、closing の内側 mismatch 項だけである。
MLclosing、各 MLintern、hairpin、stack、exterior の項は一致した。
この構造は固定した検査用構造であり、この配列の MFE ではない。
具体的な全項は [README の表](README.md)、[term_breakdown.tsv](term_breakdown.tsv)、raw trace に保存した。

差は常に同じ向きではない。
右端に接する GU closing の `ml_gapped_GU_pad10` は、固定版 −0.70、ViennaRNA −0.20であった。
片方だけの dangle は、標準の二つの内側隣接塩基による mismatch より有利にも不利にもなる。
両端を padding して closing を配列内部に移した例、および設計した dangles=0 の対照では、この差が消えた。

出典：[structure_metrics.tsv](structure_metrics.tsv)、[structure_accounting.tsv](structure_accounting.tsv)、[accounting_summary.json](accounting_summary.json)。

## MFE と平滑化

ViennaRNA の MFE と同じ pair set を出したのは310 / 337本だった。
残る27本は、両候補を共通の PF スコアで再評価すると同点の別解だった（最大スコア差2.14e-14以内）。
固定版と内側隣接対照の MFE 出力は337 / 337本で同じ構造・同じエネルギーだった。
本集合では closing guard の変更による MFE の変化を検出しなかった。
他の配列や finite beam まで MFE が不変だと一般化しない。

全337本で、LinearCParty が出した MFE を ViennaRNA の整数評価に揃えると、ViennaRNA の最適エネルギーと一致した。
一方、LinearCParty の PF 用スコアと整数 MFE 値の最大差は0.01698166 kcal/molだった。
これは、PF の平滑化による差であって、closing guard による差ではない。
例 `random_010_n20` は同じ MFE 構造について、整数評価 −0.30、両実装の PF スコア −0.3169816578337869となった。
exterior の一つの neighboring-base 項の平滑化でこの差を説明できた。
直接の `eval_structure` の C float 丸めも別に保存した。

出典：[mfe_comparison.tsv](mfe_comparison.tsv)、[accounting_witnesses.json](accounting_witnesses.json)。

## Source の原因と引き継ぎ

固定版 `src/scfg2/replay/w_final_exact_inside.hh:430–455` は closing の pair type を逆にし、次の条件で内側隣接塩基を渡す。

```cpp
si = j < owner_.length() ? owner_.encoded_base(j-1) : -1;
sj = i > 1 ? owner_.encoded_base(i+1) : -1;
```

`C(1,n)` では両方が −1となる。
一方、exterior と multiloop branch は、外側隣接塩基に対して i>1、j<n を検査し、その条件は今回の ViennaRNA と合っている。
`src/scfg2/weights/log_loops.hh:78–86` は、二つの隣接塩基なら mismatch、一つなら dangle、なしなら隣接項なしとし、terminal AU/GU と MLintern を加える。
したがって closing の二つの guard を変えるだけで、今回の構造差と ensemble 差が消える。
parameter の読み込みの違いを原因とする証拠は見つからなかった。

著者・実装担当の判断が必要なのは次の二点である。

1. 標準 Turner/ViennaRNA の closing assignment へ修正するか、現行の端に依存する局所スコアをモデルとして明示するか。
2. 修正する場合、実装側の検証と新しい固定版を用意し、どの結果を再計算するか。

finite beam も同じ重み context を使う。
今回の短い空 G の差から、長鎖・非空 G・pseudoknot の解析例の影響量は決められない。
観測用対照をそのまま採用版にせず、実装担当へ渡す。
指定された新規文書 `../LinearCParty/documents/scfg2/energy_boundary_check_20261017.md` に source と witness をまとめた。
実装本体、既存文書、原稿、bib は変更していない。
commit も行っていない。

非空 G の ensemble は、交差する extension を禁止する公開設定を確認できなかったため省いた。
ViennaRNA で G の対を強制するだけでは、両実装の構造集合が一致しない。
非空 G での固定構造の照合だけを報告する。
band-spanning multiloop と crossing stack/internal loop の独立な検証は未実施である。
また MLbase=0 のため、非対合塩基数や loop 別出力の正しさをエネルギー一致だけから主張しない。

## 検証・保存

本比較は2026-10-17 02:17:09–02:42:52 JSTに逐次実行した。
観測器の ensemble 計算の最大 wall time は9.27秒で、設定した45秒 timeout以内だった。
性能測定としては利用しない。
記録した source、parameter、library、FASTA、実行した主 script の前後 hash は全て一致した。
provenance が保存した結果ファイルの hash も全て一致した。

既存の固定 `LinearCParty --exact` と観測器を2個の短い入力で照合し、全 BPP の最大差は0、MFE の構造と表示エネルギーも2 / 2で一致した。
その出力と executable hash は [observer_cli_parity.json](observer_cli_parity.json) に保存した。
全ての設計例の SCFG2 trace で、局所エネルギーの和と root エネルギーが一致することも検査した。

生出力は約32.2 MBの `raw_outputs.tar.gz` 一つにまとめた。
README と小さな集計を通常の Git 履歴へ保存できるように分け、生出力を追加するための git 操作はしていない。
担当ファイルの制限を優先し、共有の `docs/progress.md` には追記していない。
