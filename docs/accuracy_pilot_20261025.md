# 構造予測精度の小規模 pilot（2026-10-25 依頼）

実行は 2026-10-26 00:03:51–00:08:28 JST。
**43 本で計測した。各群 20・20・3・0 本であり、長い群の比較は未完了である。**
保存済みの完全な参照注釈を使った予備調査として保持する。現時点では論文本文への精度比較表の追加を見送るのが妥当と考える。

50–199 nt の交差対 F 値は LinearKnot MFE 30.66%、centroid 17.72%、IPknot 10.38%、ThreshKnot 5.08% だった。
200–499 nt では LinearKnot 両設定が 0%、IPknot が 10.49%、ThreshKnot が 2.31% だった。
500–999 nt は 3 本のみで、LinearKnot 両設定 0%、IPknot 0.90%、ThreshKnot 3.12% だった。
全対 F 値では centroid が各群で 38.51%、48.62%、45.87% だった。少数かつ保存済み cohort に条件づけた結果なので、手法一般の順位を確定する用途には使わない。

## データと抽出

E1 と同じ FASTA 26,063 本（SHA-256 `c76025522f1f9e174c25207595bc079e454deb656acbfd3a51de71c60c167bc5`）に、保存済みの完全な DBN 注釈を ID と配列の両方で照合した。
参照は `results/local/pk-recovery-bprna-20261016/pk_sel.json` の 200 本と、`results/local/pk-recovery-long-rrna-20261016/mid300_700_selection.json` の 80 本に由来する。
**抽出の母集団は、この保存済み SI cohort 280 本である。** 元の E1 全集合に注釈を付けて選び直す作業は、公式 zip を再取得できる時点に残った。
元 zip の保存先は既に存在しない旧 node-local `/tmp/9028807.1.all.q/...` だった。今回の公式 URL への `curl` は名前解決に失敗した。

参照構造内に幾何学的に交差する対が一つ以上ある RNA を候補とした。
Python `random.Random(20261025)` を一つ使い、ID でソートした候補から、短い群から順に `sample(pool, min(20, len(pool)))` で選んだ。
表示・実行順は長さ、ID の順。予測結果による選び直しは行っていない。
群の端点が重ならないよう、200、500、1,000 nt を次の群に割り当てた。

| 長さ群 (nt) | E1 の配列数 | 完全参照がある候補 | 抽出数 | 実際の長さ | 参照全対 | 参照交差対 | 非 canonical 参照対（うち交差） |
|---|---:|---:|---:|---|---:|---:|---:|
| 50–199 | 22084 | 166 | 20 | 58–199 | 514 | 376 | 87 (46) |
| 200–499 | 2888 | 111 | 20 | 225–454 | 1819 | 485 | 109 (33) |
| 500–999 | 327 | 3 | 3 | 572–644 | 497 | 44 | 27 (1) |
| 1,000–3,000 | 628 | 0 | 0 | — | 0 | 0 | 0 (0) |

50–199 nt の 20 本はすべて RFAM 由来。200–499 nt は RFAM 11、RNP 5、CRW 3、tmRNA 1 本。500–999 nt の 3 本は CRW 由来である。
Canonical は AU、UA、CG、GC、GU、UG。参照の非 canonical 対も主集計に含めた。
bpRNA 注釈を参照構造として評価し、全レコードを実験だけで決まった構造とは解釈しない（[公式 help の Method / Validation Technique](https://bprna.cgrb.oregonstate.edu/help.php)）。

ArchiveII は EternaFold のローカル `archiveII_original.fasta`（3,992 レコード）が見つかったが、参照構造は含まれていなかった。
併せて手元の EternaFold `input_data.zip` の StructureData 2,793 注釈を確認したところ、交差を持つものは 0 本だった。
利用できる pseudoknot 参照が揃わなかったため、ArchiveII の精度計算は省略した。この調査だけで ArchiveII の全配布版の性質を断定しない。
候補・選択 ID・seed・各参照の出所・配列と参照対は [selection.json](../results/local/accuracy-pilot-20261025/selection.json) と [cohort.jsonl](../results/local/accuracy-pilot-20261025/cohort.jsonl) に保存した。

## 道具と設定

| 設定 | 使用したもの |
|---|---|
| LinearKnot centroid / MFE | `v0.2.0` → `f93b55ed9749dd1c198cd1f1809aece50f12f7d8`。`../../LinearKnot` から `git archive v0.2.0` で node-local に取り出し、CMake Release、g++ 11.5.0、`-j8 --target linearknot` で build。`b=100`、`h=d=100`、ranking `outer`、centroid は `--gamma 1`。 |
| Scaffold | LinearPartition-V `b=200`、dangles 2、確率 `p≥0.7`、連続 2 対以上。既存 `scripts/prepare_biological_inputs.py:threshold_scaffold` をそのまま呼んだ。参照注釈は scaffold 作成に使っていない。 |
| LinearPartition MEA | `deps/LinearPartition-20261008`、commit `b450fb3e63189073b68d385589035f992080aa3a`。`-V -b 200 --dangles 2 -M -g 1`。γ はこの pilot の固定設定。 |
| ThreshKnot | 同じ LinearPartition、`-V -b 200 --dangles 2 -T --threshold 0.3`。 |
| IPknot | ローカル Docker archive 内の 1.1.0。`-e LinearPartition-V --beam-size 200 -n 1`、既定の refinement 1 と auto/auto thresholds。 |
| Spark † | `deps/Spark-upstream` commit `3105c311ca6243a31788e74ac33640435e53b62c`、既存の未変更 Release binary。既存 software lock の binary hash と一致。同じ scaffold、`-P rna_Turner04.par -d2 -r G S`。 |
| Knotty | GitHub source の取得を試したが DNS 解決失敗。binary なしとして全 RNA を skip。入手できた場合は 200 nt 以下、30 s 上限で実行する設定を runner に保存した。 |

IPknot の LinearPartition-V 指定は [公式 README](https://github.com/satoken/ipknot) にあるモデル設定を使用した。
新規 clone は DNS 解決に失敗し、既存 SIF の Apptainer 起動も `getsockopt: operation not permitted` で失敗した。
既存 `PseudoknotInverseFolding/submodules/ipknot/ipknot_docker.tar` から必要な binary と libraries を node-local に取り出し、musl loader で通常プロセスとして起動できた。
image の SHA-256 と実行 binary / library の hash を保存した。image の正確な build 元 source は今回独立には認証していない。近くの source checkout は `dcaa3695041552df05f4665f3e4227cd617b37f6` だが、image の build revision と同一とは主張しない。
LinearPartition の MEA / ThreshKnot オプションは [公式 README](https://github.com/LinearFold/LinearPartition) にも記載されている。

LinearKnot と Spark の通常ループにはバイト一致の Turner 2004 parameter file（SHA-256 `0c671a7455f5686bc160fc655194b11888594ab7d506626cdad1470bcadd17d9`）を使い、37°C、dangles 2、GU を許した。
Spark との既存の比較では、PK 局所項の丸め、BE の扱い、WI→WMB 項などの意味の差も記録されている（[2026-10-12 の監査](progress.md#2026-10-12-offline-spark-mfe-comparison)）。
今回そろえたのは parameter bytes と nominal な設定であり、全エネルギー意味論の一致は前提にしない。
**Spark の精度は traceback の既知の不具合を含む参考値（†）とする。** 診断用の修正版は使っていない。

## 指標と計測範囲

全対は 1-based の座標集合を比較した。括弧の種類、1 塩基ずれの許容は使わない。
交差対の集合には、その構造の中で別の対と `i<k<j<l` または `k<i<l<j` を満たす全ての対を含めた。
参照・予測それぞれで交差集合を求め、その集合どうしの共通部分を TP とした。
参照で交差する対と同じ座標を予測しても、予測構造で交差しなければ交差 TP にはならない。
このため pseudoknot-free な MEA の交差 F は 0 になる。参照交差座標が全予測対のどこかに見つかる個数も JSONL に別保存した。

主表は群内の RNA を合計する micro 集計：sensitivity = TP / 参照対数、PPV = TP / 予測対数、F = 2 TP / (参照対数 + 予測対数)。
予測交差対が 0 の場合、PPV は未定義（—）、参照交差対があるので sensitivity と F は 0。
RNA 単位の平均 F も `summary.json` に保持した。
異常終了・不正出力・skip は数と理由を残し、精度と成功時間の分母には成功 RNA だけを使った。Spark の欠測による cohort 差は共通成功集合の表で確認した。

既存の 16-slot allocation（host `r11n8`、JOB_ID `9139876`、NSLOTS `16`）で、single-thread の CLI を常時一つずつ実行した。
測定は各入力・設定 1 回、計画的な warm-up は 0 回。別に実施した 3 RNA・21 process の動作 smoke は主集計から除外した。
LP、LinearKnot、IPknot の各 process は 120 s、Spark と Knotty は 30 s、AS limit は 4 GiB。
同一手法・群で 2 timeout 後は残りを skip、Spark 累計 180 s、pilot 全体 7,200 s の上限を設けた。実際の timeout は 0 件だった。
主計測は 301 process（scaffold 43、予測 258）で、成功 299、Spark の不正出力・異常終了が各 1。Knotty は別に 43 skip を記録した。
計測本体は 277.25 s（4.62 分）。大きな job の投入や `qsub` はしていない。

時間は `/usr/bin/time` で包んだ command の外側の wall-clock elapsed time、メモリは GNU time の peak RSS。
下の資源表の「LP+CLI」は、LinearKnot / Spark には scaffold の LinearPartition 時間を加え、他の道具には単独 command の時間を使った。
controller 側の threshold 処理・ファイル準備はこの合計に含めない。RSS は別々に実行した LP と decoder の最大値で、加算しない。
**LinearKnot MFE の `fold` CLI は Viterbi に加え inside / ensemble free energy も計算する。Spark は MFE 計算なので、この CLI 時間から MFE kernel 単独の速度比は求めない。**
繰り返しとノード占有を揃えた正式な性能測定は今後の計画に残る。

## 精度（%）

S = sensitivity、P = PPV。各手法で成功した RNA を集計。† は Spark の参考値。

### 50–199 nt

| 設定 | 成功/選択 | 全対 S | 全対 P | 全対 F | 交差 S | 交差 P | 交差 F |
|---|---:|---:|---:|---:|---:|---:|---:|
| LinearKnot centroid γ=1 | 20/20 | 38.33 | 38.70 | 38.51 | 11.17 | 42.86 | 17.72 |
| LinearKnot MFE | 20/20 | 40.08 | 34.33 | 36.98 | 23.40 | 44.44 | 30.66 |
| LinearPartition MEA γ=1 | 20/20 | 37.55 | 36.83 | 37.19 | 0.00 | — | 0.00 |
| ThreshKnot | 20/20 | 37.74 | 33.92 | 35.73 | 2.93 | 19.30 | 5.08 |
| IPknot LP-V | 20/20 | 37.74 | 32.50 | 34.92 | 7.71 | 15.85 | 10.38 |
| Spark † | 19/20 | 40.28 | 35.41 | 37.69 | 21.91 | 43.58 | 29.16 |

### 200–499 nt

| 設定 | 成功/選択 | 全対 S | 全対 P | 全対 F | 交差 S | 交差 P | 交差 F |
|---|---:|---:|---:|---:|---:|---:|---:|
| LinearKnot centroid γ=1 | 20/20 | 48.49 | 48.76 | 48.62 | 0.00 | 0.00 | 0.00 |
| LinearKnot MFE | 20/20 | 50.47 | 42.56 | 46.18 | 0.00 | 0.00 | 0.00 |
| LinearPartition MEA γ=1 | 20/20 | 49.97 | 45.96 | 47.88 | 0.00 | — | 0.00 |
| ThreshKnot | 20/20 | 50.63 | 44.34 | 47.28 | 1.65 | 3.83 | 2.31 |
| IPknot LP-V | 20/20 | 51.07 | 42.71 | 46.52 | 12.37 | 9.10 | 10.49 |
| Spark † | 19/20 | 50.56 | 44.53 | 47.35 | 0.00 | 0.00 | 0.00 |

### 500–999 nt

| 設定 | 成功/選択 | 全対 S | 全対 P | 全対 F | 交差 S | 交差 P | 交差 F |
|---|---:|---:|---:|---:|---:|---:|---:|
| LinearKnot centroid γ=1 | 3/3 | 47.48 | 44.36 | 45.87 | 0.00 | — | 0.00 |
| LinearKnot MFE | 3/3 | 49.70 | 41.79 | 45.40 | 0.00 | 0.00 | 0.00 |
| LinearPartition MEA γ=1 | 3/3 | 47.28 | 42.11 | 44.55 | 0.00 | — | 0.00 |
| ThreshKnot | 3/3 | 48.49 | 42.58 | 45.34 | 2.27 | 5.00 | 3.12 |
| IPknot LP-V | 3/3 | 50.10 | 41.16 | 45.19 | 2.27 | 0.56 | 0.90 |
| Spark † | 3/3 | 48.49 | 40.92 | 44.38 | 0.00 | 0.00 | 0.00 |

1,000–3,000 nt は完全参照がないため全手法で未実施。500–999 nt の 3 本は 572、612、644 nt であり、群全体の精度推定には小さすぎる。

全六設定が成功した共通集合の交差 F 値を次に示す。Spark の一部の有効出力にも scaffold 脱落があるため、† の位置づけは変わらない。

| 長さ (nt) | 共通成功数 | LK centroid | LK MFE | LP MEA | ThreshKnot | IPknot | Spark † |
|---|---:|---:|---:|---:|---:|---:|---:|
| 50–199 | 19 | 18.50 | 30.71 | 0.00 | 5.33 | 10.76 | 29.16 |
| 200–499 | 19 | 0.00 | 0.00 | 0.00 | 2.36 | 10.97 | 0.00 |
| 500–999 | 3 | 0.00 | 0.00 | 0.00 | 3.12 | 0.90 | 0.00 |

## 時間・メモリ（成功 RNA の中央値）

「CLI」は与えた scaffold に対する command の時間、または配列単独からの command の時間。
「LP+CLI」と peak RSS は前節の pipeline の範囲。異なる機能範囲を含む CLI の観測値である。

| 長さ (nt) | 設定 | n | CLI (s) | LP+CLI (s) | LP+CLI peak RSS (MiB) |
|---|---|---:|---:|---:|---:|
| 50–199 | LinearKnot centroid γ=1 | 20 | 0.051 | 0.103 | 9.00 |
| 50–199 | LinearKnot MFE | 20 | 0.044 | 0.099 | 9.00 |
| 50–199 | LinearPartition MEA γ=1 | 20 | 0.053 | 0.053 | 9.00 |
| 50–199 | ThreshKnot | 20 | 0.054 | 0.054 | 9.00 |
| 50–199 | IPknot LP-V | 20 | 0.173 | 0.173 | 6.00 |
| 50–199 | Spark † | 19 | 0.011 | 0.064 | 9.00 |
| 200–499 | LinearKnot centroid γ=1 | 20 | 2.641 | 2.866 | 43.50 |
| 200–499 | LinearKnot MFE | 20 | 1.136 | 1.438 | 40.50 |
| 200–499 | LinearPartition MEA γ=1 | 20 | 0.383 | 0.383 | 10.50 |
| 200–499 | ThreshKnot | 20 | 0.379 | 0.379 | 9.00 |
| 200–499 | IPknot LP-V | 20 | 2.288 | 2.288 | 15.00 |
| 200–499 | Spark † | 19 | 0.026 | 0.390 | 9.00 |
| 500–999 | LinearKnot centroid γ=1 | 3 | 5.996 | 7.224 | 72.00 |
| 500–999 | LinearKnot MFE | 3 | 3.157 | 4.385 | 75.12 |
| 500–999 | LinearPartition MEA γ=1 | 3 | 1.234 | 1.234 | 27.00 |
| 500–999 | ThreshKnot | 3 | 1.229 | 1.229 | 24.00 |
| 500–999 | IPknot LP-V | 3 | 4.390 | 4.390 | 27.00 |
| 500–999 | Spark † | 3 | 0.084 | 1.312 | 24.00 |

IPknot の最長 1 process は 43.07 s（`bpRNA_CRW_2185`, 572 nt）だった。実行した process の timeout は 0 件だった。
GNU time の wall time、user/system time、個別 CLI の RSS、失敗時の出力も JSONL / tar に保存した。

## 失敗・観察・論文に入れる価値

- Spark `bpRNA_RFAM_7554`（113 nt）：return code 0 だが `[` が 3、`]` が 5 の不正な構造。精度は欠測とした。
- Spark `bpRNA_RFAM_35962`（349 nt）：signal 11、GNU time return code 139。精度は欠測とした。原因の箇所をこの結果から特定しない。
- Spark `bpRNA_RFAM_4615`（299 nt）：構造自体は parse できたが scaffold の 18 対を脱落。返された構造の精度は参考値として保持し、脱落座標を `analysis.json` に保存した。
- Knotty：GitHub source の取得が DNS エラー。ArchiveII：使用できる pseudoknot 参照が揃わず省略。IPknot：container 起動は失敗したが native loader による実行は全 43 本で成功。

短い群では LinearKnot MFE の交差 F が centroid より高く、全対 F は centroid の方が高かった。
200–499 nt では IPknot が交差参照対を 60 対回復し、LinearKnot 両設定は 0 対だった。
この群には canonical な交差参照対が 452 対あり、500–999 nt にも 43 対あった。これらについても LinearKnot の交差 TP は 0 対だった。
この条件での crossing recovery の弱さは、参照の非 canonical 対の数だけでは説明できない。
Scaffold、beam / 局所制限、エネルギーモデルの寄与を分けるには、同じ入力で scaffold と幅を替える診断が必要になる。

本 pilot の価値は、配列から scaffold を予測する実用設定で、全対と交差対の成績を分けて観察できたことにある。
現在の結果は、短い RFAM への偏り、保存済み SI cohort への条件づけ、長い群の不足、Spark の出力不具合、単回の CLI 計測を含む。
**当面は作業メモに保持する。精度の比較を論文に入れる判断は、完全注釈から群を揃えた追試の後がよい。**
その場合も、既存 SI の「参照から作った scaffold に対する回復」と、今回の「予測 scaffold からの構造予測」を明示して分ける。

## 本格的に行う場合

1. 公式 bpRNA-1m(90) zip を node-local に取得する。前回記録された zip hash は `4756eb63d705ec32ba838660846c65ba36ad2f66a1fefecab1010ad59ed46320`。配布版が違えばその差を調査し、E1 FASTA と ID・配列を照合する。幾何学的交差で選別してから各群約 20 本を固定する。群の候補が足りない場合は全例と欠測理由を明示する。
2. 元データ由来・RNA family・非 canonical 対の数を監査し、全注釈に対する主評価と、canonical 対だけの補助評価を事前に定義する。Macro / micro の双方、RNA 単位の paired 差と区間推定、全手法の共通成功集合、全選択入力に対する成功率を保存する。
3. IPknot の build 元 source と solver を固定する。今回の LP-V 設定を主条件にするならその理由を明記し、既定 LP-C を補助条件にするか先に決める。Knotty は取得できた場合に短い群だけで試す。Spark の不正出力は修復して混ぜず、版と意味の差を分けて管理する。
4. 1,000–3,000 nt は先に各 2–3 本で時間・メモリを校正し、小さい allocation に収まる上限を決める。single-thread / single-worker、必要なら warm-up 1 + 3 repetitions。MFE 単独の資源比較には v0.2.0 API の Viterbi だけを呼ぶ observer を使い、`fold` の inside 計算を混ぜない。Scaffold LP と decoder の時間を別に測り、pipeline 合計も示す。
5. 少数の診断では threshold scaffold と参照由来 scaffold、b=100/200、短い RNA の exact を比較する。参照対がモデル内で到達可能か、エネルギー差、scaffold により禁止される対を分けて調べる。主精度比較の入力選択と設定は結果を見る前に固定する。

## 保存物・検証・再実行

実行スクリプトは [accuracy_pilot_20261025.py](../scripts/accuracy_pilot_20261025.py)。
結果は [results/local/accuracy-pilot-20261025/](../results/local/accuracy-pilot-20261025/) に保存した。

| ファイル | 内容 |
|---|---|
| `cohort.jsonl`, `selection.json` | 完全な配列・参照 DBN / 対集合・出所、候補全 ID、seed、群、source hashes |
| `results.jsonl` | 344 行。301 実行（299 成功、2 Spark 不正/失敗）と 43 Knotty skip。command、状態、対集合、全指標、時間・RSS、raw member hashes |
| `raw_outputs.tar` | 1 実行ごとの入力、BPP、実 scaffold、stdout、stderr、time、command / result JSON と script / input metadata の snapshot |
| `run.json` | source / binary / parameter provenance、実際の設定・順序・環境・計測範囲、前後の binary hashes |
| `summary.json`, `analysis.json` | 群別の全指標・資源中央値、macro F、共通成功集合の再集計、参照内訳、Spark 欠測・scaffold 脱落、source 前後照合 |
| `validation.json` | 独立な座標再集計、scaffold 規則、parse、hash の検査結果 |
| `setup.tar.gz` | node-local source archive と build log / binary / cache、IPknot runtime、取得失敗の生ログ、初期 smoke、集計・検証の補助コード |
| `artifact_hashes.json` | この結果 directory の保存物の SHA-256 |

500 組の無作為 pair maps を、別実装 `scfg2_pair_accuracy.py` と比較し、全対・交差対の TP、参照数、予測数、F が一致した。
参照では交差する座標が予測では交差しない小例と、不完全・非対称 BPSEQ を別に検査した。
主計測の全 256 有効予測を tar の stdout から再 parse し、参照と予測それぞれの交差判定と全指標を独立に再計算した。
全 43 scaffold は BPP `p≥0.7` と隣接する対の有無から別に再構成し、一致した。
raw member 2,193 個の hash、入力 / archive の hash、binary の前後 hash、外部 source の commit と status、Spark の既存 lock との一致を確認した。
`validation.json` は `passed`、検査エラー 0 件である。

既存結果の検証と再集計は次で再実行できる（再予測なし）。

```bash
PYTHONDONTWRITEBYTECODE=1 python3 scripts/accuracy_pilot_20261025.py validate --out results/local/accuracy-pilot-20261025
PYTHONDONTWRITEBYTECODE=1 python3 scripts/accuracy_pilot_20261025.py summarize --out results/local/accuracy-pilot-20261025
```

再集計は JSON の作成日時を更新するので、後で保存物の hash を取り直す。保存済み artifact hashes は今回提示した版の値である。
完全注釈を入手した後の新しい抽出は、既存 directory を上書きせず次のように行える。

```bash
PYTHONDONTWRITEBYTECODE=1 python3 scripts/accuracy_pilot_20261025.py prepare   --bprna-zip /node-local/bpRNA_1m_90.zip --per-bin 20   --out results/local/accuracy-pilot-full-new
```

`run` の実際の呼び出しは `run.json:invocation`、個々の CLI は `results.jsonl:argv` に保持した。
Binary を `setup.tar.gz` から node-local に展開して各 path を指定する。`run` は測定済み directory の上書きを拒否する。中断後の自動 resume は実装していない。
共有領域の quota は開始時 21,864,090 / 28,000,000 inode。新規共有 entry の予算を 16 個とし、各 process のファイルは全て node-local で作成して tar / JSONL にまとめた。
Build と実行用の node-local 中間物は保存後に除去した。外部実装 checkout の変更・checkout・commit、原稿の変更、git commit はしていない。
