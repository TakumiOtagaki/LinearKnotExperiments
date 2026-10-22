# Scaffold の分岐祖先段数（2026-10-23）

保存済みの入力と scaffold G だけを読み、`d` が制限する CASE3 の分岐祖先段数を数えた。
RNA の計算、scaffold の再構築、scheduler への投入は行っていない。

## 結果

| データ | RNA 本数 | 長さ（nt） | 最大段数 | 30 を超える RNA 数 | 30 を超える端点数（d=30 で打ち切り） |
|---|---:|---:|---:|---:|---:|
| bpRNA 400 nt 以上（E1 long） | 1,438 | 400–4,381 | 16 | 0 | 0 |
| NEAT1_1、scaffold 閾値 0.7 | 1 | 3,756 | 8 | 0 | 0 |
| NEAT1_2、scaffold 閾値 0.7 | 1 | 22,743 | 26 | 0 | 0 |

対象は計 1,440 本で、`d=30` によって祖先経路が打ち切られる端点はない。
bpRNA の最大値 16 は `bpRNA_RFAM_43116`（3,665 nt）で得られた。
分母となる G の自由塩基数は、bpRNA long が 1,176,310、NEAT1_1 が 2,302、NEAT1_2 が 14,421。
同じ自由塩基を、複数の相手や beam 条件ごとに重複して数えていない。

- [summary.tsv](summary.tsv)：データごとの本数、最大、30 を超える RNA 数・端点数。
- [per_rna.tsv](per_rna.tsv)：1,440 本について、最大、自由塩基ごとの段数分布、打ち切り端点数、入力・S・G の SHA-256。
- [validation.json](validation.json)：入力・コードの hash、独立した数え方との照合、既存 RNAcentral 7 本との照合。

## 使用した定義

既存の [scaffold_features.py](../../../scripts/scaffold_features.py) の `scaffold_features(sequence, scaffold, case3_depth=30)` をそのまま呼び出した。
各自由塩基の owner は、その塩基を最も内側で囲む scaffold 対合（外側なら仮想 root）である。
owner 自身を除く祖先のうち、直下に二つ以上の対合の子を持つ節点を一段と数える。
仮想 root もこの条件を満たして経路上にあれば数える。
外側の自由塩基の段数は 0。

- `case3_ancestors_max`：全自由塩基について、上限を適用する前の段数の最大。
- `case3_truncated_endpoints`：段数が **30 より大きい**自由塩基の数。30 ちょうどは打ち切られない。
- `rnas_over_30`：この端点を一つ以上持つ RNA の数。

E1 short の既存列 `maximum_branching_ancestors_of_free_base` は、[rna_records.py](../../../scripts/rna_records.py) の `describe_scaffold` が計算しており、owner 自身が分岐ならそれも数える。
外側の自由塩基についても、root が分岐ならそれを数えるため、E2 の CASE3 の定義と完全には同じではない。
今回の TSV にはこの E1 の列も併記した。
対象三群では、群ごとの最大はこの定義でも順に 16、8、26 だった。

既知の bpRNA 400 nt 未満の「24,625 本、最大 5」は E1 の列による値であり、今回は再集計していない。
RNAcentral 7 本の既存 TSV は **d=100** で作られていた。
今回、保存入力から同じ E2 関数で段数分布を再計算し、全 7 本の分布・最大・S/G hash を照合したうえで、d=30 の打ち切り端点数を数え直した。
最大は既知の値と同じ 28、d=30 での打ち切り端点数も全て 0 だった。

## 入力と再現方法

bpRNA long は [cohort manifest](../../../manifests/cohorts/bprna-400plus-20261012.json) の全 1,438 本を対象にし、採用した E1 long run `results/runs/e1-bprna-long-89baea52-cpu16-20261018/inputs.tar` の凍結入力を展開せずに読んだ。
run の `array_plan.json` が持つ archive hash・各入力 hash と、cohort・scaffold dataset の hash を確認した。
配列の名前・長さ・hash、scaffold hash、既存 E1 の分岐段数も scaffold dataset と一致した。

NEAT1 は [data/neat1-20261016/manifest.json](../../../data/neat1-20261016/manifest.json) と `NEAT1_1.txt`、`NEAT1_2.txt` を使用した。
manifest の配列長・S/G hash・対合数が一致し、scaffold 設定は LinearPartition-V、beam 200、dangles 2、対合確率閾値 0.7、最小 helix 長 2 である。
SHAPE を用いた別 scaffold や閾値の異なる scaffold は含めていない。

集計用の [analyze_scaffold_branching_depth.py](../../../scripts/analyze_scaffold_branching_depth.py) は既存の E1/E2 関数を変更せずに呼ぶ。
論文 repo の `experiments/` から、次のコマンドで結果の三ファイルを作成した。
再実行する場合は `--out` に新しいディレクトリを指定する。

```bash
python3 -B scripts/analyze_scaffold_branching_depth.py \
  --out results/local/scaffold-branching-depth-20261023
```

## 検証記録

対象 1,440 本と参照用 RNAcentral 7 本の全自由塩基について、owner から親を明示的にたどる独立の集計と、既存 E2 関数の段数分布が一致した。
各分布の合計は自由塩基数と一致し、段数が 30 より大きい部分の和は打ち切り端点数と一致した。

```bash
python3 -B -m unittest discover -s tests -p test_scaffold_features.py -v
```

既存の 5 tests が成功した。
さらに既存の合成 scaffold 生成関数で段数 30 と 31 の境界を確認し、d=30 の打ち切り端点数がそれぞれ 0 と 6 になることを確認した。

共有領域の quota は開始前に確認し、使用 inode は 21,846,212、live 上限は 28,000,000、保守的な 24,000,000 の境界まで 2,153,788 の余裕があった。
新規共有 entry は集計スクリプト 1、結果ディレクトリ 1、結果ファイル 4 の計 6。
入力の展開、一実行ごとの出力、追加の生ログは作成していない。
原稿・外部実装は変更せず、git commit は行っていない。
