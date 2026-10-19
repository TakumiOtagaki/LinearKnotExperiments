# 解析 D: 座標集合と G との交差による pseudoknot 回復率

2026-10-19。原稿は変更せず、保存済み 89baea52 出力だけを再集計した。
新しい RNA CLI の実行、tar の展開、外部実装の変更、git commit は行っていない。

## 定義と出典

塩基対を 1-based の座標 (i,j), i<j の集合として読む。
(i,j) と (k,l) は i<k<j<l または k<i<l<j のとき交差する。
固定 scaffold の対集合を G、canonical な参照対の集合を A、全 MFE 対の集合を M とする。
Canonical は AU/UA/CG/GC/GU/UG。

- 新参照集合 T = {p ∈ A \ G : G の少なくとも一対と交差する}。
- 新予測集合 P = {p ∈ M \ G : G の少なくとも一対と交差する}。
- sensitivity = Σ|T ∩ M| / Σ|T|、PPV = Σ|T ∩ P| / Σ|P|。
  全出力で T∩M = T∩P を検査した。集計は RNA 間の分子・分母を合計する micro average。
- 確率 >0.5 は厳密な不等号。RNA 内平均は新参照集合 T 上で取り直す。
  RNA 間平均は等重み。距離群は RNA 内の参照対の j−i の中央値で分ける。
- 括弧種は座標を読むためだけに使い、参照・予測の選択には使わない。
  全追加対を pseudoknot とみなす定義とは区別する。

入力・G・出力は `results/runs/pk-recovery-{200-turner04,200-dp09,80,rrna}-89baea52-20261018/`。
200 本は exact、80 本と rRNA は b=h=d=400 / outer ranking。
既存値は `results/local/small-rerun-89baea52-20261018/{pk_results.json,mid300_700_results.json,recovery_summary.json,rrna_artifact_index.json}`。
この 89baea52 集計ディレクトリには README はなかったため、`scripts/analyze_small_rerun_89baea52.py` を確認した。
参照は `results/local/pk-recovery-bprna-20261016/pk_sel.json` と `results/local/pk-recovery-long-rrna-20261016/mid300_700_selection.json`。
rRNA は同じ旧ディレクトリの README と、89baea52 の `rrna_reference_metadata.json` を確認した。
旧 README の「全24対がゼロ」は旧実装の記録であり、新版の比較値に流用していない。

## 旧定義との比較／SI の回復表の差し替え案

`crossing_G` の行を採用する。200 本の2行で既存の回復表を差し替え、80 本と rRNA は別 panel に置く案。
表 caption には上記 T/P と micro average、exact と有限幅の違いを明記する。

| RNA 群 | モデル | 定義 | 参照対 | 予測対 | 一致対 | Sensitivity | PPV | 参照 p>0.5 | RNA 内平均 p>0.5 |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| 200 | turner04 | old_symbols | 983 | 850 | 210 | 21.36% | 24.71% | 179/983 (18.21%) | 34/200 |
| 200 | turner04 | crossing_G | 1040 | 850 | 210 | 20.19% | 24.71% | 179/1040 (17.21%) | 34/200 |
| 200 | dp09 | old_symbols | 983 | 413 | 110 | 11.19% | 26.63% | 80/983 (8.14%) | 13/200 |
| 200 | dp09 | crossing_G | 1040 | 413 | 110 | 10.58% | 26.63% | 80/1040 (7.69%) | 13/200 |
| 80 | turner04 | old_symbols | 811 | 301 | 40 | 4.93% | 13.29% | 24/811 (2.96%) | 3/80 |
| 80 | turner04 | crossing_G | 809 | 301 | 40 | 4.94% | 13.29% | 24/809 (2.97%) | 3/80 |
| rrna | turner04 | old_symbols | 24 | 14 | 0 | 0.00% | 0.00% | 0/24 (0.00%) | 0/3 |
| rrna | turner04 | crossing_G | 24 | 14 | 0 | 0.00% | 0.00% | 0/24 (0.00%) | 0/3 |

全483出力で予測 `[]` 集合と新予測集合 P は**座標ごとに完全一致**した。
予測側の符号化に由来する差は0で、PPV と一致対数は全条件で変わらない。
参照側には次の差がある。

- 200 本では、23 RNA の canonical な `{}` 対57対が G と交差するため加わる（983→1040）。
  Turner04/DP09 とも、この57対には MFE 一致も p>0.5 もない。
  そのため sensitivity は 21.36→20.19%、11.19→10.58%、p>0.5 の割合は18.21→17.21%、8.14→7.69%となる。
  RNA 内平均 p>0.5 の件数34/13は変わらない。
- 80 本では `bpRNA_CRW_4958` (133,166) と `bpRNA_CRW_4764` (242,275) が除かれる（811→809）。
  元注釈でこれらと交差する `()` 対はそれぞれ (135,167), (244,276) の UU 対だけで、G の構築時に除かれている。
  これら二対は新しい G に対する crossing pair ではない。
  両対とも MFE 一致でも p>0.5 でもないため、分子40/24は変わらない。
  sensitivity は4.932→4.944%、p>0.5 の割合は2.959→2.967%。
- rRNA の保存済み canonical 参照24対はすべて G と交差し、集合も指標も変わらない。
  MFE 一致0/24、予測14対、p>0.5 は0/24、報告値ゼロ9/24、最大確率0.00045326084596248895。

### 参照から除かれる対の内訳

| 群 | canonical 参照対 | G の対 | 非canonical除外 | G外で非交差 | 旧集合に追加 | 旧集合から除外 |
|---|---:|---:|---:|---:|---:|---:|
| 200 | 6100 | 5042 | 725 | 18 | 57 | 0 |
| 80 | 7082 | 6261 | 544 | 12 | 0 | 2 |
| rrna | 1499 | 1475 | 98 | 0 | 0 | 0 |

200/80 本の canonical 参照は元の全括弧種を読む。元注釈の対数は canonical と非canonical の和。
G に残す条件は従来どおり canonical な `()` 対かつ j−i≥4。
Canonical でも G 外で非交差の対は、新しい pseudoknot 参照集合には含めない。

rRNA の利用可能な保存物には元の完全な DBN がなく、G と canonical な `[]`/`{}` 参照座標24対がある。
したがって rRNA の新集合は、この凍結された参照座標と G に同じ交差判定を適用した結果である。
上表の rRNA の非canonical除外数は保存済み metadata の合計で、元DBNから再分類した数ではない。
元DBNの未保存の括弧種を含めた完全な再監査とは主張しない。
`[]`/`{}` を参照選定に使い分ける方針は廃し、保存された候補座標について同じ幾何学的判定を使う。

### 確率・距離の他の指標

200 本 Turner04 の距離表も分母と平均を更新する案（全条件は `span_summary.tsv`）。

| j−i (nt) | 参照対 | p>0.5 | 対ごとの平均確率 |
|---|---:|---:|---:|
| le50 | 706 | 150 (21.25%) | 0.209387 |
| 51to100 | 183 | 25 (13.66%) | 0.134448 |
| gt100 | 151 | 4 (2.65%) | 0.016233 |

RNA 内の参照距離中央値 ≤100 / >100 の群数は、200 本では174/26、80 本では37/43のまま。
群内の等重み RNA 平均確率は、200 本 Turner04 で0.18995349→0.18549277 / 0.02350847（不変）、
DP09 で0.08843006→0.08478907 / 0.00015225（不変）、
80 本で0.08624786→0.08333152 / 0.00488461（不変）。
RNA平均に使う対集合を変更したためで、各対の確率を再計算した差ではない。
今回の1040対という新しい分母を、旧983対に対する別の到達可能性検査の分母へ自動で移さない。

## 英文案（SI、3文）と対応する日本語

Reference pairs are canonical annotation pairs that cross at least one scaffold pair in G, and predicted pairs are new MFE pairs that cross G.
We matched pairs by nucleotide coordinates and obtained sensitivities of 210/1,040 (20.19%) and 110/1,040 (10.58%) with Turner 2004 and DP09, respectively; PPVs were 210/850 (24.71%) and 110/413 (26.63%).
The fractions of reference pairs with probability above 0.5 were 179/1,040 (17.21%) and 80/1,040 (7.69%), respectively.

参照集合は G の少なくとも一対と交差する canonical な注釈対、予測集合は G と交差する新規 MFE 対とする。
座標で照合すると、Turner04/DP09 の sensitivity は210/1040（20.19%）/110/1040（10.58%）、PPV は210/850（24.71%）/110/413（26.63%）である。
参照対のうち確率0.5を超える割合は179/1040（17.21%）/80/1040（7.69%）である。

## 再現と検証

```bash
OPENBLAS_NUM_THREADS=1 python3 -B scripts/analyze_pk_recovery_pairset_89baea52.py \
  --out results/local/pk-recovery-pairset-89baea52-20261019
python3 -B scripts/report_local_errors_pairset_89baea52.py recovery
OPENBLAS_NUM_THREADS=1 python3 -B -m unittest discover -s tests -p test_local_errors_and_pairset_recovery.py -v
```

解析の出力先は既存ディレクトリを上書きしないため、再集計時は別の出力先を使う。
`rna_metrics.tsv` は483条件、`pair_metrics.tsv` は新旧参照の和集合、`definition_changes.tsv` は全変更座標、
`reference_audit.tsv` は参照除外数、`summary.{tsv,json}` は集計、`span_summary.tsv` は距離別の値。
`provenance.json` は入力・plan・task manifest/checkpoint・解析コード・出力のSHA-256を記録する。
全480旧MFE集計と全483旧確率ベクトルを保存済み89baea52集計に完全一致で照合した。
raw artifact の checksum、保存コマンド、入力 S/G、MFE の G 保持、全参照について T∩M=T∩P を検査した。
境界0.1/0.5、疎なBPPの欠落ゼロ、括弧色の交換、交差判定、非canonical除去を6 testsで検査した。
科学的な集合・指標の変更と、これらの開発上の検査を分けて扱う。
