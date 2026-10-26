# 実 RNA による CParty PF / Spark MFE 資源比較（2026-10-25 計画、10-26 集計）

**2026-10-26 集計済み。全63 taskの結果と原稿への反映は末尾を参照。以下の「準備」記録は投入前の状態であり、未開始という記述は当時のもの。**

準備作業は 2026-10-25 に開始し、10-26 JST に検証を完了した。準備担当は **qsub・git commit を実行していない。**
正式計画は 107 入力、15 arrays / 63 tasks、warm-up 508 回＋測定 1,524 回＝2,032 CLI。
全 63 task は `not_started`、attempt は空である。原稿・SI の数値と既存の図は編集していない。

## 入力の固定

E1 の bpRNA-1m(90) 26,063 本から、各中心の整数 ±10%（両端を含む）の候補を ID 昇順に並べる。
`random.Random(20261025)` を一つ作り、中心 250→500→1000→2000→4000 の順で、各群から重複なしで 20 本を `sample` する。
選んだ群内は長さ・ID 順に保存する。群が足りない場合は全候補とする規則だが、今回は全群に 20 本以上あった。
family、scaffold の複雑さ、計算の成功・失敗では選別していない。

| 中心 nt | 選択範囲 nt | E1 内の候補 | 固定本数 | 選択された実際の長さ nt |
|---:|---:|---:|---:|---:|
| 250 | 225–275 | 479 | 20 | 225–275 |
| 500 | 450–550 | 145 | 20 | 451–550 |
| 1,000 | 900–1,100 | 95 | 20 | 905–1,077 |
| 2,000 | 1,800–2,200 | 88 | 20 | 1,805–2,158 |
| 4,000 | 3,600–4,400 | 41 | 20 | 3,607–4,065 |
| RNAcentral E2 | 4,793–29,106 | 7 | 7 | 4,793–29,106 |
| 合計 | | | **107** | |

RNAcentral は URS000233CDDE (4,793)、URS00019626CE (5,482)、URS000239D8A4 (6,990)、
URS00023C9544 (11,832)、URS000241043C (13,652)、URS000238C72C (22,543)、URS00023A4DB9 (29,106)。
採用 E2 の入力・scaffold をそのまま使う。

bpRNA の scaffold も作り直さず、採用 E1 の保存 bytes を使った。
LinearPartition-V `b450fb3e`、beam 200、dangles=2、posterior 閾値 0.7、連続 helix 最小 2 対、MEA なしという同じ手順である。
cohort、scaffold dataset、採用 E1 short/long plan の配列 hash と S/G hash、E2 の入力 hash を照合した。
全 ID、群、長さ、family、配列・scaffold・三行入力の SHA-256、元 archive/member は
[data/comparison-resources-20261025/manifest.json](../data/comparison-resources-20261025/manifest.json) にある。
共有側の入力は `inputs.tar` と manifest の 2 ファイルであり、107 個の入力ファイルを展開しない。

- 選択 manifest SHA-256: `2a9535f75d5e1cc69ed3de380125f8000eacbf02c827d9fff0af910a709c2140`
- 入力 tar SHA-256: `57325191ccab2e879090f0e4768fef4b020915e7f93231dce40bc27b066bd2c5`

選定を再現する script は `scripts/prepare_comparison_inputs_20261025.py`。既存の固定出力を上書きしない。

## 固定したソフトウェアと測る範囲

全実行は Turner 2004（SHA-256 `0c671a7455f5686bc160fc655194b11888594ab7d506626cdad1470bcadd17d9`）、37 °C、dangles=2。
1 プロセス 1 thread とし、OMP/BLAS/MKL と locale を固定する。G の作成と build は測定に含めない。

| 方式名 | 版・実行内容 | 設定 |
|---|---|---|
| `cparty` | CParty `800e95a`＋論文と同じ 6 patch、既存 E0 scalar observer | PF のみ、MFE/outside なし、sampling=0、scaling seed=0。zero-sample 後処理は含む |
| `exact` | LinearKnot v0.2.0 `f93b55ed9749dd1c198cd1f1809aece50f12f7d8` の inside API | PF のみ、beam/h/d 上限なし、内部 loop 制限は保持 |
| `beam100` | 同 v0.2.0 の inside API | PF のみ、b=h=d=100、outer ranking、内部 loop 制限は保持 |
| `spark` | Spark `3105c311ca6243a31788e74ac33640435e53b62c` の公式 CLI | `-d 2 -P Turner04 -r G S`、GU と PK を許可、`-k` なし、ソース修正なし |
| `fold100` | 同 v0.2.0 `linearknot fold` 全体 | b=h=d=100、outer ranking、BPP 出力なし。**MFE と PF を含む** |

PF の 3 方式は `LegacyCompatible`、すなわち完成 WMB の weak closure なしにそろえた。
本文の既定の `RequireWeakClosure` と区別する。`fold100` はその既定を使う。
`fold --help` と `src/scfg2/cli/fold.cc` を確認したが、MFE だけを求めるオプションはなく、BPP を要求しなくても viterbi と inside の両方を呼ぶ。
Spark 対 fold はこの要求出力の範囲を明記した資源比較とし、energy・構造の一致は主張しない。
10-12 の整数丸め、GU 分岐、nested PK の局所モデル差、scaffold 対が抜ける traceback の確認事項は残す。Spark の修復は行わない。

`build_comparison_resources_20261025.py` で各ソースを node-local に取り出した。
LinearKnot は指定どおり `git archive v0.2.0`、Spark は固定 commit の archive、CParty は固定 upstream の archive に既存 patch を順番に適用した。
Spark の全 tracked file と CParty の patch・修正後ソースの hash を従来の software manifest と照合した。
既存の外部 repo を変更・checkout していない。build は TSUBAME `r11n8`、GCC 11.5、Release、build worker 2、LinearKnot LTO link worker 1。
Spark の既存 `-march=native` を含む flags、全 build command、compiler、CPU、各 library と executable の hash は保存済み。

- lock: [manifests/software-comparison-resources-v0.2.0-20261025.json](../manifests/software-comparison-resources-v0.2.0-20261025.json)
- node-local build の原ソース・ログ・CMakeCache・compile_commands: `results/local/comparison-resources-build-20261025/build_raw.tar.gz`
- node-local 展開用の固定ソース・実行ファイル: 同 `software.tar.gz`（約 5.1 MiB）
- build コマンドと環境: 同 `build_provenance.json`

各 task は既に build した同じ hash の bundle を node-local `/tmp` に展開する。測定時の再 build はない。
照合用の Git index だけを node-local に作り、元実装の objects を読み取り参照して source commit と clean 状態を検証する。
外部実装の作業 tree は変更しない。adapter が呼ぶ native executable も `build_files` で実行前後に検証する。
resume のため node-local パスは task ごとに固定し、同一パスへの展開と task の両方に排他 lock を使う。

実際の native コマンドは次のとおり。全方式の計測は `Comparison -i INPUT --method METHOD` という同じ Python adapter を経由する。
adapter は `exec` で native tool に置き換わる。wall-clock は起動、adapter、tool とその出力を含むため、特に短い値の解釈でこの範囲を保つ。

```bash
E0CPartyOracle -i INPUT
PFScaling -i INPUT --method exact
PFScaling -i INPUT --method beam100
Spark -d 2 -P TURNER04 -r SCAFFOLD SEQUENCE
linearknot fold -i INPUT -b 100 --max-hairpin 100 --max-case3-depth 100 --ranking outer
```

## 試走と失敗・打ち切りの扱い

6 実 RNA × 5 方式 × warm-up 1＋測定 3＝120 回を行い、116 回成功、同じ 1,805 nt RNA の CParty 4 回が `nonfinite oracle inside` だった。
下は各入力の測定 3 回の中央値（秒 / MiB）。node 内最大 8 single-thread processes の開発試走であり、本実験の群全体の性能値ではない。

| 入力 | nt | CParty PF | exact PF | b=100 PF | Spark | fold b=100 |
|---|---:|---:|---:|---:|---:|---:|
| bpRNA_RFAM_38321 | 225 | 0.066 / 6.0 | 0.166 / 14.0 | 0.367 / 18.0 | 0.066 / 6.0 | 0.367 / 18.0 |
| bpRNA_RFAM_9076 | 226 | 0.066 / 6.0 | 0.116 / 15.0 | 0.217 / 15.0 | 0.066 / 6.0 | 0.167 / 12.0 |
| bpRNA_RFAM_10156 | 451 | 0.317 / 21.0 | 1.469 / 45.4 | 2.770 / 63.0 | 0.117 / 6.0 | 3.572 / 75.0 |
| bpRNA_RFAM_2435 | 454 | 0.317 / 21.0 | 1.218 / 45.4 | 2.321 / 60.0 | 0.116 / 6.0 | 3.172 / 73.8 |
| bpRNA_CRW_54685 | 905 | 2.120 / 70.9 | 10.387 / 166.5 | 8.687 / 146.5 | 0.317 / 6.0 | 12.040 / 162.7 |
| bpRNA_RFAM_37115 | 1,805 | nonfinite (3/3) | 42.755 / 649.2 | 14.902 / 221.2 | 1.569 / 6.0 | 20.717 / 264.3 |

有限な CParty/exact の 5 RNA × 3 測定＝15 組は、log Z の絶対差 1e-9、PF energy の絶対差 1e-10 kcal/mol を全て通った。
最大差はそれぞれ 4.04e-12、2.5e-12。
これは PF の照合であり、Spark の energy/構造の照合はしていない。

同じ最終 staging runner でも 2 RNA・40 CLI が成功した。
resume 後も完了済み 40 行の JSONL hash が同じで、200 個の raw member の hash を再照合した。CLI の再実行はない。
tar の末尾 padding は checkpoint recovery で正規化されるため、tar 全体の hash が変わっても member の保存内容で検証する。
別の打ち切り対照では、native Spark の 0.1 秒 timeout、adapter の 0.001 秒 timeout、exact の 64 MiB AS cap による `std::bad_alloc` を保存した。
timeout で GNU time が完了しなかったときの RSS は null のままである。

- 試走 config・JSONL・tar・数値 summary: `results/local/comparison-resources-pilot-20261025/`
- staging/resume の開発 plan: `results/runs/comparison-resources-staging-smoke-20261025/`
- 全 dry-run の検証結果: [comparison_resources_20261025_validation.json](comparison_resources_20261025_validation.json)

CParty の nonfinite は失敗として保持し、scaling seed や実装を変更して回避しない。
timeout とメモリ上限はその条件での打ち切りとして保持する。AS cap は仮想アドレス空間の上限であり、RSS の上限と同じではない。
raw の elapsed/RSS は失敗診断として残るが、失敗や打ち切りの elapsed を成功時の資源値に置き換えない。
`analyze_comparison_resources_20261025.py` は 3 測定の成功と有限出力を必要とし、欠測・nonfinite・timeout・memory_limit を分母に残す。
Spark の既知の traceback 問題があるため、出力を受理しても構造の正しさや G 保持を保証したとは扱わない。

## Array と実行上限

各入力・方式に warm-up 1、測定 3、入力×方式の seeded blocked_random（seed 20261025）。warm-up 完了後に最初の測定を行い、測定 block 間に barrier を置く。
時間は既存 runner の `elapsed_seconds` と GNU time の wall、メモリはプロセスごとの最大 RSS を保存する。
集計する値は各入力・方式の測定 3 回の中央値である。

全 task は `cpu_16=1`、計画メモリ 35,095 MiB、task overhead 1,024 MiB。
2,200 nt までは pool 最大 8、3,600 nt 以上は最大 4。32 GiB/process を与える長い CParty/Spark は memory envelope を守るため pool 1 とする。
実際の同時 CLI 数は pool 上限と実行可能な combination 数の小さい方となる。異なる pool 数の結果では同時実行の影響も条件として報告する。
短い 3 群の全 5 方式は同じ task に置く。2,000/4,000 nt と RNAcentral は CParty、Spark、exact を別 array に分けた。
この長い範囲の異なる array の方式は同一 node の対応比較ではないため、host と workers を併記する。

exact は 2,200 nt まで（80 本）、それ以上の 27 本は `exact_exclusions` に明記した。
1,805 nt で約 43 秒の測定中央値を確認したことから 2,000 nt 群を加えた。
4,000 nt 以上は今回の exact 対象に含めず、実行不能と断定せず、資源値も外挿しない。
長い CParty/Spark の 1 CLI 上限は 1 時間である（24 時間は採用しない）。これにより 1+3 の各実行の打ち切りを task の h_rt 内で記録できる。

| run 末尾 / tier | 入力 | tasks | pool | AS cap MiB/process | timeout 秒/CLI | h_rt/task | h_rt 全使用 point |
|---|---:|---:|---:|---:|---:|---|---:|
| short / n250 | 20 | 5 | 8 | 1,024 | 10 | 00:10:00 | 0.0400 |
| short / n500 | 20 | 5 | 8 | 1,024 | 30 | 00:25:00 | 0.1000 |
| short / n1000 | 20 | 5 | 8 | 4,096 | 120 | 01:20:00 | 0.3200 |
| bounded / n2000 | 20 | 5 | 8 | 4,096 | 180 | 00:45:00 | 0.1800 |
| bounded / n4000 | 20 | 10 | 4 | 8,192 | 300 | 00:50:00 | 0.4000 |
| bounded / rnacentral | 7 | 4 | 4 | 8,192 | 900 | 02:15:00 | 0.4320 |
| cparty / n2000 | 20 | 3 | 8 | 4,096 | 180 | 00:45:00 | 0.1080 |
| cparty / n4000 | 20 | 5 | 4 | 8,192 | 600 | 01:30:00 | 0.3600 |
| cparty / rna_4793_6990 | 3 | 1 | 4 | 8,192 | 900 | 02:15:00 | 0.1080 |
| cparty / rna_11832_29106 | 4 | 4 | 1 | 32,768 | 3,600 | 05:05:00 | 0.9760 |
| spark / n2000 | 20 | 3 | 8 | 4,096 | 180 | 00:45:00 | 0.1080 |
| spark / n4000 | 20 | 5 | 4 | 8,192 | 600 | 01:30:00 | 0.3600 |
| spark / rna_4793_6990 | 3 | 1 | 4 | 8,192 | 900 | 02:15:00 | 0.1080 |
| spark / rna_11832_29106 | 4 | 4 | 1 | 32,768 | 3,600 | 05:05:00 | 0.9760 |
| exact / n2000 | 20 | 3 | 8 | 4,096 | 900 | 03:30:00 | 0.5040 |
| 合計 | 107 unique | **63** | | | | 105.833 allocation h | **5.0800** |

各表の run 接頭辞は `comparison-resources-v0.2.0-20261025-`。
config は [configs/comparison-resources-v0.2.0-20261025/](../configs/comparison-resources-v0.2.0-20261025/)、
script は `jobs/<run>/<tier>.sh`、凍結 plan は `results/runs/<run>/array_plan.json`。
機械可読の一覧・除外一覧は [comparison_resources_20261025_plan.json](comparison_resources_20261025_plan.json)。

h_rt は過去の実測と今回の小規模完走値を参照して短い群の CLI 上限を選び、その上限が全て発生しても各 1+3 matrix を収めるよう設定した。
長い未測定入力の完了時間は推定しない。packing は長さのべき乗の当てはめを使わず、実際の CLI timeout から計算する。
combination 数 C、pool w、各 CLI 上限 T について、既存 runner が warm-up 直後の最初の測定を待つ起動順を保った。
その待ちのため次の combination の起動が遅れるので、warm-up の仕事量 `CT` を直列として全て加算する。
測定 1 block の上限を、C≤w なら `T`、それ以外なら `T(C/w + 1 − 1/w)` とし、その 3 block 分と staging/checkpoint の 120 秒を加える。
すなわち `120 + T[C + 3 max_block_waves]` を task 上限とし、h_rt はこの上限が 80% 以下になる 5 分刻みの値とした。
同時実行が完全に満たされるという仮定を置かず、20 combination・8 workers の先頭待ちを独立に再現した検査でも上限が実行順を覆うことを確認した。
従来の accuracy packing は維持し、この timing 用の timeout-bound option と adapter の依存実行ファイル検証を追加した。

この timeout 予算に基づく課金見積もりは **4.02417 point**、全 task が h_rt を使う場合は **5.0800 point**。

通常優先度 0、`0.060 × (0.7 × max(task wall 秒,300) + 0.1 × h_rt 秒) / 3600`。
[2026 年度の公式料金表](https://www.t4.cii.isct.ac.jp/tsubame-kyodo/tariff2026) を 2026-10-26 JST に確認した。
0.0001 point 未満の切捨て、ストレージ、再投入は計算に含めない。これらの上限は実消費時間の予測ではない。

## ファイル数と dry-run

共有 entry は準備・初回 task 実行・scheduler log・atomic temporary を含め **1,019 entries** の枠。CLI ごとの共有ファイルは作らず、payload は 63 task × JSONL/tar の **126 ファイル**。
live quota は 21,864,590/28,000,000、保守的 24,000,000 まで残り 2,135,410。残りの 10% の gate を全 plan と全体予算が通った。

これは予約ではないため投入直前にも確認する。再投入の scheduler log は別途増える。
node-local build と各 CLI の入力・出力・一時ファイルは node-local に作り、保存する raw は tar にまとめた。検証後の再生成可能な node-local build は削除した。
通常 plan は依存パスを node-local に固定しない offline dry-run で作り、runtime bundle と全実行ファイル hash を staging 時に再照合する。
最終 checker は全 plan・63 task config/schedule・107 入力 hash・15 shell script・全 dependency hash・memory/wall/inode・全 task 未開始を検証した。
未開始 plan の集計でも全 1,524 測定が missing として残り、成功値は 0 件であることを確認した。
新規 4 tests、既存 array 11、software lock 2、packed storage 8 の計 25 tests が成功した。

```bash
cd /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments
python3 -B scripts/check_comparison_resources_20261025.py
# 単独 task の実行なしの確認例
python3 -B scripts/comparison_resources_array_20261025.py \
  --plan results/runs/comparison-resources-v0.2.0-20261025-short/array_plan.json \
  --tier n250 --task-id 1 --dry-run
```

## 後で投入するコマンド（この作業では実行していない）

既存 array と合わせて平日 30 task 以下にする。下は自分の比較を最大 4 task とし、`-sync y` で array を一つずつ進める。
まず short の 3 arrays の結果と失敗率を確認する。nonfinite や timeout を含む array の終了コードは非 0 になり得る。
失敗を成功に書き換えず、記録を見て以降の投入を進める。指定 resource/h_rt を qsub 側で変更しない。

```bash
cd /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments
: "${TSUBAME_GROUP:?Set the authorized accounting group}"
COMPARISON_TASK_BUDGET=4  # 同時に動く他の全 array を差し引き、合計 30 以下
lfs quota -p 102571 /gs/bs
qsub -g "$TSUBAME_GROUP" -sync y -tc "$COMPARISON_TASK_BUDGET" jobs/comparison-resources-v0.2.0-20261025-short/n250.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$COMPARISON_TASK_BUDGET" jobs/comparison-resources-v0.2.0-20261025-short/n500.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$COMPARISON_TASK_BUDGET" jobs/comparison-resources-v0.2.0-20261025-short/n1000.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$COMPARISON_TASK_BUDGET" jobs/comparison-resources-v0.2.0-20261025-bounded/n2000.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$COMPARISON_TASK_BUDGET" jobs/comparison-resources-v0.2.0-20261025-bounded/n4000.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$COMPARISON_TASK_BUDGET" jobs/comparison-resources-v0.2.0-20261025-bounded/rnacentral.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$COMPARISON_TASK_BUDGET" jobs/comparison-resources-v0.2.0-20261025-cparty/n2000.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$COMPARISON_TASK_BUDGET" jobs/comparison-resources-v0.2.0-20261025-cparty/n4000.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$COMPARISON_TASK_BUDGET" jobs/comparison-resources-v0.2.0-20261025-cparty/rna_4793_6990.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$COMPARISON_TASK_BUDGET" jobs/comparison-resources-v0.2.0-20261025-cparty/rna_11832_29106.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$COMPARISON_TASK_BUDGET" jobs/comparison-resources-v0.2.0-20261025-spark/n2000.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$COMPARISON_TASK_BUDGET" jobs/comparison-resources-v0.2.0-20261025-spark/n4000.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$COMPARISON_TASK_BUDGET" jobs/comparison-resources-v0.2.0-20261025-spark/rna_4793_6990.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$COMPARISON_TASK_BUDGET" jobs/comparison-resources-v0.2.0-20261025-spark/rna_11832_29106.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$COMPARISON_TASK_BUDGET" jobs/comparison-resources-v0.2.0-20261025-exact/n2000.sh
```

停止後に同じ task を再開する例（古い process が止まってから）：

```bash
qsub -g "$TSUBAME_GROUP" -sync y -tc 1 -t 2-2 \
  jobs/comparison-resources-v0.2.0-20261025-cparty/rna_11832_29106.sh --resume
```

成功・失敗・timeout の完了レコードは再実行しない。別 node に移った場合は task attempts の host を確認し、タイミングの取り合わせを記録する。
上限や設定を変える再試行には新しい run を作る。

## 後の集計

全 task が停止した後、既存 `aggregate_array.py` で run ごとの集計を更新できる。
今回の全方式・欠測・打ち切りを共通の表にするには次を使う。出力先は毎回新しくする。

```bash
python3 -B scripts/analyze_comparison_resources_20261025.py \
  --out results/local/comparison-resources-completed-20261025
```

CParty/exact の有限な対応測定だけで従来の PF 許容差を確認し、Spark と fold の energy/構造一致は計算しない。
共通成功入力に限った資源比較では採用本数と元の分母を併記する。全て失敗・打ち切りだった群は成功値なしとする。
今回の開発試走を正式な比較結果として原稿に移さず、正式 array 後に群別の完了率と対応関係を確認する。


## 2026-10-26 全 task の検証・集計と原稿反映

実行済みの15 arrays / 63 tasks、2,032 calls（warm-up 508、測定1,524）を検証した。
Task status は success 45 / failed 18だが、**failed task も全予定実行を保存済みで、欠測は0**。
全107 S/G、凍結plan/config/schedule、software lock/bundle/native binary、前後のsource/binary状態、JSONLの一意性・全1+3反復、生出力10,160 memberのbyte数とSHA-256を照合した。
`signal_11`（return code 139 / -11）は `bad_alloc` と分け、メモリ不足とは同定しなかった。

- [完了後検証JSON](comparison_resources_20261026_completed_validation.json)：taskごとのhost、workers、上限、status、全planとtask/run/JSONL/tarのhash。
- [入力ごとの中央値と群集計](../results/analysis/comparison-resources-20261026-final/summary.json)、[resources.tsv](../results/analysis/comparison-resources-20261026-final/resources.tsv)、[length_groups.tsv](../results/analysis/comparison-resources-20261026-final/length_groups.tsv)、[対応比](../results/analysis/comparison-resources-20261026-final/paired_comparisons.tsv)。
- `check_comparison_resources_20261025.py --completed` を追加。従来の投入前/pilot検証の経路も保持。
- 既存analyzerを長さ群と対応入力比の集計に拡張。中央値は入力ごとに3測定、次に成功入力間で中央値。三測定すべてが成功して有限であることを必要とし、失敗の時間/RSSから成功値を作らず、外挿しない。

### 長さごとの資源

各cellは **秒 / MiB**。群の長さ範囲は上の入力表どおりで、RNAcentralの行は各1本。
各手法の成功集合を使い、CPartyは1,000 nt群で19本、Sparkは500 nt群で19本、他の有限なbpRNA群は20本。
LinearKnot exactは計画80/80、PF b=100とfoldは各107/107成功。

| nt / 群 | CParty PF | exact PF | b=100 PF | Spark MFE | fold PF+MFE | CParty成功 | Spark成功 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 250群 | 0.116 / 9 | 0.217 / 16.6 | 0.417 / 21 | 0.0665 / 6 | 0.468 / 21 | 20/20 | 20/20 |
| 500群 | 0.367 / 24 | 1.19 / 53.2 | 2.35 / 61.5 | 0.116 / 6 | 3.15 / 70.5 | 20/20 | 19/20 |
| 1000群 | 2.52 / 78 | 10.5 / 186 | 7.38 / 140 | 0.367 / 6 | 9.84 / 160 | 19/20 | 20/20 |
| 2000群 | --（nonfinite） | 52.8 / 743 | 12.8 / 216 | 1.8 / 6 | 17 / 246 | 0/20 | 20/20 |
| 4000群 | --（nonfinite） | --（未測定） | 26.7 / 444 | 14 / 12 | 35.1 / 506 | 0/20 | 20/20 |
| 4,793 | --（nonfinite） | --（未測定） | 32.2 / 549 | 27.3 / 18 | 45.7 / 579 | 0/1 | 1/1 |
| 5,482 | --（nonfinite） | --（未測定） | 55.7 / 751 | 46.4 / 24.8 | 82.7 / 875 | 0/1 | 1/1 |
| 6,990 | --（timeout） | --（未測定） | 65.8 / 939 | 87.4 / 35.1 | 95.1 / 980 | 0/1 | 1/1 |
| 11,832 | --（signal_11,timeout） | --（未測定） | 87.2 / 1.31e+03 | 581 / 90.8 | 126 / 1.45e+03 | 0/1 | 1/1 |
| 13,652 | --（timeout,signal_11） | --（未測定） | 157 / 1.94e+03 | --（signal_11） | 247 / 2.19e+03 | 0/1 | 0/1 |
| 22,543 | --（signal_11,timeout） | --（未測定） | 240 / 2.86e+03 | --（timeout） | 339 / 3.12e+03 | 0/1 | 0/1 |
| 29,106 | --（memory_limit） | --（未測定） | 219 / 3.19e+03 | --（timeout） | 306 / 3.3e+03 | 0/1 | 0/1 |

### 失敗と比較できる範囲

CPartyは59/107本で成功、48/107本で失敗した。入力単位の排他的な内訳は、43 nonfinite、1 timeoutのみ（6,990 nt）、3 timeout/signal 11混在（11,832、13,652、22,543 nt）、1 memory allocation failure（29,106 nt、AS cap 32 GiB）。
Nonfiniteの43本は994 ntの1本、2,000/4,000 nt群40本、4,793/5,482 ntの2本。
測定単位では成功177、nonfinite129、timeout7、signal 11が5、memory_limit3。
Warm-upは成功59、nonfinite43、timeout2、signal 11が2、memory_limit1。
Nonfiniteは実装の数値の範囲の問題であり、アルゴリズムの限界とはしない。
AS capは仮想アドレス空間の上限で、RSSの上限ではない。signal 11の原因は保存ログだけでは同定しない。

共通成功の実RNA59本（225–1,077 nt）は全てCPartyがLinearKnot exactおよびb=100より速かった。
**1,805 nt以上のCPartyには有限な成功測定がないため、実RNAでLinearKnotが速くなる長さは決められない。**
合成入力の2,000 ntでの交差は既存の別実験（SI表S20）であり、実RNAへ外挿しなかった。
CPartyの成功群のRSS中央値は9→24→78 MiB。長さごとの成功集合が違うことも明記した。
CParty/exactの有限177測定組（59入力）は全てPF許容差内：最大log Z差1.1823431123048067e-11、energy差7.30437932361383e-12 kcal/mol。

Sparkは103/107本で成功。550 ntと13,652 ntの2本が全4回signal 11、22,543 ntと29,106 ntの2本が全4回3,600秒timeout。
測定単位では成功309、signal 11が6、timeout6。構造・energy・予測精度の一致は評価していない。
下の比は**共通成功入力ごとのSpark/fold比の中央値**であり、群中央値を割った値ではない。1より大きい時間比はfoldが速い。

| 群 / RNA ID | 共通成功本数 | Spark/fold 時間比 | Spark/fold RSS比 |
|---|---:|---:|---:|
| bprna_250 | 20 | 0.1414 | 0.2857 |
| bprna_500 | 19 | 0.03984 | 0.08696 |
| bprna_1000 | 20 | 0.03725 | 0.03741 |
| bprna_2000 | 20 | 0.1049 | 0.0244 |
| bprna_4000 | 20 | 0.4239 | 0.02734 |
| URS000233CDDE | 1 | 0.5978 | 0.03107 |
| URS00019626CE | 1 | 0.5616 | 0.0284 |
| URS000239D8A4 | 1 | 0.9182 | 0.03578 |
| URS00023C9544 | 1 | 4.625 | 0.0626 |

11,832 ntはfold 125.703 s / 1,449.922 MiB、Spark 581.365 s / 90.770 MiBで、foldが4.6249倍速い一方、RSSは15.9737倍。
29,106 ntのfoldは305.978 s（5.09964 min）、3,300.289 MiB（3.22294 GiB）。Sparkはtimeout、CPartyはbad_allocであり、これらから速度比を作らなかった。

**依頼文のBPPとの相違**：この比較の`fold100`コマンドは`--bpp`を持たない。
v0.2.0の`src/scfg2/cli/fold.cc`ではBPP要求なしの場合にviterbiとinsideを実行し、outsideは呼ばない。
したがって新比較は「fold全体＝PF+MFE、BPP出力なし」と記載した。既存E2のBPP込み測定（29,106 ntで約7 min / 3.22 GiB）はFigure 4A/B・SI表S16に別測定として保持した。
BPP込みの追加測定は行っていない。

### 図・英日原稿

- [新Figure 4と移動した合成比のSI図](../results/figures/figure4_resources_20261026_final58/)：Figure 4は178×58 mm、主要7 pt・凡例6.5 pt・最小6 pt、PNG 300 dpi。旧BPP込みA/Bを保持し、新C/DにPF比較3方式、Spark MFE、fold PF+MFEを重ねた。失敗は資源値のない別の帯で色と形を区別し、48/107と4/107を注記。
- 旧Figure 4Cの合成比はSupplementary Figure S3。既存Figure S2の全曲線と役割を分け、数値は変更しない。
- Results 3.3、caption、Methods 2.5、Abstract、SI S3.6の英日を更新。合成PF表S20を保持し、実RNA PFをS21、SparkをS22に整理。表S16などの英日参照を照合。
- Abstract比較文：`On real RNAs, CParty was faster through 1,077 nt, but numerical or resource failures prevented longer comparisons.`
  続く29 kb文：`At beam width 100, LinearKnot computed the partition function and MFE for a 29-kilobase RNA in 5.1 minutes and 3.22 GiB.`
- Sato feedback S04/S32を済に更新。外部実装repoの変更・checkout、RNA再計算、qsub、git commitは行っていない。

再生成（出力先は未使用directoryを指定）：

```bash
python3 -B scripts/check_comparison_resources_20261025.py --completed
python3 -B scripts/analyze_comparison_resources_20261025.py --out results/analysis/NEW_COMPARISON
python3 -B scripts/make_figure4_resources_20261026.py --summary results/analysis/NEW_COMPARISON/summary.json --out results/figures/NEW_FIGURE4
python3 -B -m unittest discover -s tests -p test_comparison_resources.py
```

新たな重い実験は投入していない。rawを展開せず、task単位のtar/JSONLを読み、小さい集計と図だけを追加した。20 MB超のファイルをGitへ追加していない。


最終検証：`cd manuscript/latex && bash build.sh` と `bash make_diff.sh 5e3b503` 成功。
通常本文8頁・SI23頁、差分本文9頁・SI23頁。4 logのエラー・未定義参照/引用0。
図の最終PNG、通常本文全頁概観とAbstract/Figure 4、通常/差分SIの比較表、差分本文全頁概観と追加・削除を画像で確認した。
英日追加24行・65 label・42種の参照、既存SI結果表の数値不変、全9 export hashを照合した。
通常/差分SIにoverfull・大きすぎるfloatなし。通常本文の既存class警告と差分Figure 1の既存vbox警告は残るが、画像上の欠け・重なりはない。
原稿側の照合とPDF hashは [resources_validation.json](../../docs/feedback_20261025_sato/resources_validation.json) に保存した。
