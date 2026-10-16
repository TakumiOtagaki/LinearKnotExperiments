# Energy assignment versus ViennaRNA, 2026-10-17

担当範囲：エネルギー割り当ての独立照合（依頼の Major 2・3）。
第2回 global review の C17、S6、T4、T10 にも対応する。
原稿、日本語原稿、bib、実装 checkout は変更していない。
git commit はしていない。
最終集計、検証結果、著者への引き継ぎは [REPORT.md](REPORT.md) に示す。
以下では分配関数を PF、塩基対確率を BPP、最小自由エネルギー構造を MFE と略記する。

## 分かったこと

SI 表 S3 の closing factor `C(i,j)` は固定版 source を正しく記述している。
ただし、その境界条件は標準 ViennaRNA の dangles=2 と異なる。
`C(1,n)` は内側の隣接塩基を両方省く。
片端だけに接する closing では一方だけを使い、通常の mismatch の代わりに dangle を選ぶ。
同じ完成構造で、この条件だけによる差を確認した。
外側の枝の条件や、paired な隣接塩基も使うこと自体を誤りとはしていない。

比較の全数、成功数、最大差は [summary.tsv](summary.tsv) に記録する。
`pinned` は採用固定版の重みをそのまま使う観測プログラムである。
`inner_neighbors_diagnostic` は、node-local の一時 header で closing の二つの guard だけを取り除いた観測用対照である。
実装の修正版・新しい採用版として扱わない。

## 比較条件と入力

| 項目 | 条件 |
|---|---|
| LinearCParty | d02d261c3a07e7308a3a3d37543d040adfd0f86e、既存の固定 build の library をリンク |
| 独立参照 | Python ViennaRNA 2.7.2、別にインストールされた native library |
| パラメータ | 固定 source の `params/rna_Turner04.par` を両方でロード |
| parameter SHA-256 | `0c671a7455f5686bc160fc655194b11888594ab7d506626cdad1470bcadd17d9` |
| 共通設定 | 37°C、dangles=2、betaScale=1、pf_smooth=1、special hairpin 有効 |
| 構造制約 | pseudoknot-free、最小 hairpin 3、内部 loop の非対合塩基合計 ≤30、AU/CG/GU 対を許す |
| その他 | noLP=noGU=noGUclosure=logML=circ=gquad=0、pair span と exact hairpin 長の追加上限なし |
| Vienna 固有の設定 | salt=1.021（基準濃度、補正なし）、uniq_ML=1（sampling 用） |
| 主比較の G | 空 |
| 指標 | exact log Z、ensemble free energy、全 i<j の BPP、MFE 構造とエネルギー |
| 資源制限 | 逐次実行、観測プロセスごとに address space 2 GiB、CPU 40秒、wall timeout 45秒 |

provenance の `limits` は C++ 観測器に対する記録である。
ViennaRNA native 部分の thread 数は計測していない。
保存した elapsed time は実行記録であり、性能比較には使わない。

ランダム 150 本は長さ 20、40、80、120、200 nt 各30本である。
bpRNA 150 本は保存済みの clean FASTA から、50–99、100–149、150–200 nt を各50本、seed 20261017 で選んだ。
保存 FASTA の SHA-256 は `c76025522f1f9e174c25207595bc079e454deb656acbfd3a51de71c60c167bc5` である。
新しいデータ取得は行っていない。
設計した 37 入力（5–27 nt）は、GC/AU/GU closing、両端・片端・配列内部、非対合塩基のある/ない multiloop、hairpin、隣接する exterior stems を含む。
入力の配列・ID・hash は [inputs.jsonl](inputs.jsonl) に保存した。
この配列集合は、生物学的な予測精度を評価するための benchmark ではない。

## 丸めと平滑化を分けた独立評価

ViennaRNA の `mfe()` と `eval_structure()` は整数パラメータのエネルギーを使う。
固定版 LinearCParty の MFE は PF 用の平滑化された重みを最適化する。
このため、同じ構造でも両者のスコアに小さな差が出ることがある。
これを closing の境界差と混ぜないため、次の三つを保存した。

1. `RNA.fold_compound.eval_structure(R)` の直接の評価値。
2. R の全対を強制し、全非対合位置を `x` で固定した、独立 ViennaRNA の一構造だけの PF スコア。
3. ViennaRNA 自身の `exp_E_Hairpin`、`exp_E_IntLoop`、`exp_E_MLstem`、`exp_E_ExtLoop` を使った項別の和。

一構造だけの PF では、指定対の BPP が1になることも検査した。
`pf()` が Python に返す free energy は C float に丸められるため、log Z は double precision の `exp_matrices.q1k[n]` と `pf_scale` から復元した。
`log Z = log(q1k[n]) + n log(pf_scale)` とし、RT=0.6163207755 kcal/mol を両者で照合した。
通常の整数 MFE との構造の不一致は、同点の別解と、スコアの違いに分けて [mfe_comparison.tsv](mfe_comparison.tsv) に記録する。

## 構造単位の具体例

`ml_gapped_GC_pad00` は25 ntで、(1,25) が multiloop を閉じる。

```text
S  GAGGGAAAACCCAGGGAAAACCCAC
G  .........................
R  (.(((....))).(((....))).)
```

| 項（kcal/mol） | 固定版 | ViennaRNA |
|---|---:|---:|
| MLclosing | 9.30 | 9.30 |
| closing の MLintern | −0.90 | −0.90 |
| closing の内側 mismatch | 0.00 | −1.10 |
| 二つの枝の MLintern | −1.80 | −1.80 |
| 二つの枝の隣接塩基項 | −3.00 | −3.00 |
| multiloop の非対合塩基3個 | 0.00 | 0.00 |
| hairpin 2個 | 9.00 | 9.00 |
| stack 4個 | −13.20 | −13.20 |
| exterior / terminal AU/GU | 0.00 | 0.00 |
| 合計 | **−0.60** | **−1.70** |

右に A を一つ足すと、同じ closing が左端だけに接し、固定版 −2.80、ViennaRNA −3.40 となる。
両側に A を足して closing を配列内部へ動かすと、両者とも −3.20 となる。
各比較では exterior の項の変化も含めている。
GU closing が右端に接する `ml_gapped_GU_pad10` では、固定版 −0.70、ViennaRNA −0.20 である。
差の方向は一方向ではない。
paired な隣接塩基を持つ `ml_adjacent_GC_pad00` の差は1.50 kcal/molである。
ここでも dangles=2 による paired neighbor の使用は両者で共通で、closing の guard が差を生む。

固定した構造の追跡は [term_breakdown.tsv](term_breakdown.tsv) と raw archive に保存する。
設計例では dangles=0 の対照も評価した。
同じ非交差の完成構造について、全対を G にした場合、closing pair だけを G にした場合も評価した。
MFE と、最初の12ランダム配列から ViennaRNA が生成した各3標本も構造評価の対象にした。
構造評価の条件数と未到達数は [summary.json](summary.json) に示す。
同じ (S,R) の dangles や G の割り当てを変えた評価を、別の独立配列数として数えない。

## Source との対応

行番号は固定版に対するもの。

| 場所 | 確認したこと |
|---|---|
| `src/scfg2/replay/w_final_exact_inside.hh:430–455` | pair type の向き、E/M の外側隣接位置、C の内側隣接位置と外側の端による guard |
| 同 `:83–94` | VM の closing に MLclosing と closing-stem 項を一度足す |
| 同 `:109–145` | multiloop/exterior の枝の局所項 |
| `src/scfg2/weights/log_loops.hh:78–86` | terminal AU/GU、MLintern、mismatch/dangle の選択 |
| `src/scfg2/weights/runtime_parameters.hh` | 両者に渡した parameter file を固定版でロード |
| `src/ViennaRNA/params/params.c:575–778` | PF 用の連続値・平滑化・log weights |
| `src/scfg2/deductions/generate_exact_basic.hh:621–640` | 通常の内部 loop の非対合塩基合計 ≤30 |
| `src/scfg2/api/exact.cc` | empty G の inside/outside、Viterbi、固定構造評価の経路 |
| `src/scfg2/api/beam.cc` | 同じ局所重み context を finite beam も使うこと |

外側の枝 E(i,j) は、i>1 なら s[i−1]、j<n なら s[j+1] を使う。
M(i,j) は同じ位置を使い、MLintern を加える。
C(i,j) は pair type を逆にし、j<n の場合だけ s[j−1]、i>1 の場合だけ s[i+1] を使う。
全てで、GC/CG 以外の許容対には terminal AU/GU を加える。
隣接塩基を二つ渡せば mismatch、一つなら dangle、なければその項なしとなる。

## 残る判断・調べていない範囲

標準 Turner/ViennaRNA の closing assignment を採るか、現行の端に依存した局所スコアをモデルとして明記して維持するかは著者・実装担当の判断が要る。
source の「guard を保つ」というコメントだけでは、物理モデルとして意図した規約かどうかを確定できない。
修正する場合は、実装側の検証、新しい固定版、影響を受ける結果の再測定が必要になる。
観測用の一時 header を採用版に流用しない。
この比較から、既存の長鎖・pseudoknot の結果への影響量を推定していない。

G が空でない ensemble では、交差する extension を禁止する公開 API/CLI 設定を確認できなかった。
ViennaRNA で G の対を強制しただけの ensemble とは構造集合が違うため、その log Z/BPP/MFE 比較は省いた。
代わりに、非空 G を使う固定構造の評価を行った。
band-spanning multiloop、crossing stack/internal loop、pseudoknot 固有のエネルギー、任意入力の完全性・非曖昧性は独立に検証していない。
Turner 2004 の MLbase は0なので、multiloop の非対合塩基数の誤りをエネルギー一致だけで検出することはできない。
loop profile と loop 別 accessibility の正しさを今回の比較で確認したとは扱わない。
今回の結果で Major 2・3 全体を解決済みとは扱わない。

実装側の引き継ぎ文書は、指定された一つの新規ファイル
`../LinearCParty/documents/scfg2/energy_boundary_check_20261017.md` だけに書いた。
担当ファイルの制限を優先し、共有の `docs/progress.md` には追記していない。

## 再現と保存

実験 repo の root から実行する。
既存の出力を上書きしないため、再実行では新しい output 名を指定する。

```bash
python scripts/check_energy_vienna_20261017.py --output results/local/energy-vienna-check-rerun
python scripts/summarize_energy_vienna_20261017.py --input results/local/energy-vienna-check-rerun
```

source、build、FASTA はそれぞれ `--source`、`--build`、`--fasta` で指定できる。
既定の seed と入力選択条件はスクリプトと provenance に保存している。
最初に inode を確認し、使用数は21,361,065 / 28,000,000だった。
コンパイル、shadow header、CLI の確認用 BPP は node-local の一時領域で処理し、終了時に削除する。
実行ごとの source/parameter/library hash、argv、戻り値、stderr、全 BPP と trace は、一つの `raw_outputs.tar.gz` にまとめた。

Git に保存する小さな成果物は、三つの新規 script、README、REPORT.md、inputs.jsonl、provenance.json、summary.*、ensemble_metrics.tsv、structure_metrics.tsv、structure_accounting.tsv、term_breakdown.tsv、mfe_comparison.tsv、accounting_summary.json、accounting_witnesses.json、observer_cli_parity.json である。
`results/local/` は既存の ignore 対象なので、保存するときは対象ファイルを個別に `git add -f` する。
raw archive は通常の Git 履歴に追加しない。
