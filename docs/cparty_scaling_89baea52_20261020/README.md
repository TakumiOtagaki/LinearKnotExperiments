# CParty と LinearKnot の分配関数・時間・メモリ比較の準備（2026-10-20）

**未投入。qsub・git commit は行っていない。** 5 array / 5 task、105入力、warm-up込み1,092 CLI（測定819、warm-up273）を新しい run に準備した。
既存の E2 run と外部の `../LinearKnot`・`../LinearCParty` は変更していない。

## 比較する量と設定

元の `build/cparty-oracle-800e95a-6patches/cparty-oracle` は通常の CParty CLI を起動する。
その CLI は MFE、PF、既定1,000回の sampling とそこから推定する BPP、MEA、centroid 等も計算する。
これを LinearKnot の PF だけの時間と並べない。
今回の `PFScaling --method cparty` は、同じ6 patch oracleの検証済み object に link された既存 `E0CPartyOracle` を exec し、**分配関数だけを観測**する（sampling=0、スケーリング seed=0、outsideなし）。
`hfold_pf` 内の upstream の zero-sample 後処理は残り、その時間も含む。未使用の BPP/frequency 値は確率結果として扱わない。
BPP は CParty の sampling 推定と LinearKnot の outside 計算で量と精度が異なるため、今回の計画から外した。

LinearKnot exact / b=100 は、新しい外部 observer `tools/pf_scaling_89baea52.cc` が採用版 `build/linear-cparty-89baea52` の固定 library を呼ぶ。
MFE・outside・BPP は呼ばない。exact は h/d 上限なし、b=100 は本文と同じ outer ranking、h=d=100、内部ループ L=30。
モデルは Turner04、37 °C、dangles=2。パラメータ hash と RT は記録・検証する。

**互換性：全3方式を E0 と同じ `LegacyCompatible`（completed-WMB weak closure なし）に揃えた。**
6 patch oracle は `800e95a` + `0001,0002,0006,0007,0008,0009`。closure patch は含まない。
本文の LinearKnot/E2 は `RequireWeakClosure` なので、この比較は本文の時間と同じ版89baea52だが、closure と要求出力が異なる。本文の全 CLI 時間に混ぜない。
新しい lock は `manifests/software-pf-scaling-89baea52-20261020.json`。元の oracle/E0/stage lock、全 object・library・observer の hash を含み、実行前後に検証する。

## 入力

採用 E2 の `results/runs/e2-synthetic-shapes-89baea52-20261018/inputs.tar` から hash が同じ byte を固定した。
250/500/1,000/2,000/4,000 nt × 3 seed（2026101201–03）× 7形状：empty、long_stem、separated_helices、deep_branching、separated_concentrated、separated_12、balanced_branching。
branch_depth 感度対照は含めない。形状を作る際に canonical pair を埋めるので、S と G の合成対照であり、同一配列で G だけを変える実験ではない。

## 開発 pilot（論文の正式な性能値ではない）

seed 2026101201 の empty / long_stem / balanced_branching を各長さで試した。
`cparty_pilot/pilot.json` と `raw.tar` は CParty 9入力の warm-up 1・測定3をすべて保存する。
下の範囲は成功した入力ごとの3測定の elapsed median の最小～最大。peak RSS も測定反復から集計した。

| 長さ nt | 成功した入力 / 3 | elapsed 秒 | peak RSS MiB |
|---:|---:|---:|---:|
| 250 | 3/3 | 0.067–0.117 | 9.0–9.0 |
| 500 | 3/3 | 0.367–0.467 | 24.0–24.0 |
| 1,000 | 2/3 | 2.722–3.423 | 83.6–84.1 |

1,000 nt long_stem は warm-up と3測定すべて `nonfinite oracle inside`（約2.47秒、約84–86 MiB）。失敗は保存した。
別の小さい seed probe で、スケーリング seed を変えると有限値でも PF が変わった。例えば250 nt empty の log Z は seed=0 で121.51212325053449、seed=-100で122.63471329919057。
したがって数値を有限にするだけの seed 変更は採用していない。probe の source/command/output は `scaling_seed_probe.tar`。これは実装担当へ渡す開発上の確認事項で、性能や物理モデルの結論ではない。
長い CParty 入力も非有限値になる可能性がある。完走・数値一致・交差点が得られるとはまだ主張しない。

LinearKnot の追加 feasibility check：250 nt exact は約2.1–3.9秒、500 nt exact は約26–45秒。
1,000 nt exact の balanced_branching warm-up は120秒で打ち切った（RSSは未取得）。
広い1+3 pilot はこの時点で中断し、CParty の短い完全 pilot と追加の単発較正に切り替えた。
中断前の13観測は `comparison_pilot.json` / `comparison_pilot_raw.tar`、中断時の未完了ファイルは `interrupted_pilot_tail.tar` に残した。
1,000 nt b=100 の3形状は16.85–26.87秒、peak RSS 180–278 MiB。
`calibration_pilot.json` はこの観測をまとめたもの。全て単一 worker のローカル開発値で、正式 array の同時実行時の時間を保証しない。
250/500 nt の6入力で新 observer の exact と CParty の PF を照合し、全件が E0 の許容差を通った（最大絶対差 log Z 2.05e-12、PF energy 1.26e-12 kcal/mol）。`pilot_parity.json` に保存した。

## 外挿と採否

CParty は全成功 pilot の最大 `elapsed/n³` と最大 `RSS/n²` を使い、それぞれ1.5倍の余裕を付けた。
短いプロセスの起動時間も係数に含むため、単に1,000 ntの値を8倍・64倍するより保守的である。
LinearKnot exact は三次式と250→500 ntで観測した最大指数の大きい方（約3.647）を採用。b=100 は長さに比例する式。
これは資源を計画する仮定であり、理論計算量・実測の長さ依存の結論ではない。120秒の timeout は完走値に変換していない。

| 方式 | 2,000 nt 秒 / RSS GiB | 4,000 nt 秒 / RSS GiB | 採用長 |
|---|---:|---:|---|
| CParty | 89.8 / 0.844 | 718.3 / 3.375 | 250, 500, 1000, 2000, 4000 nt |
| LinearKnot exact | 11576.7 / 3.673 | 145040.7 / 14.692 | 250, 500, 1000 nt |
| LinearKnot b=100 | 80.6 / 0.814 | 161.2 / 1.628 | 250, 500, 1000, 2000, 4000 nt |

1 CLIあたり6時間・80 GiBを上限とした。さらに warm-up 1＋測定3の一式を同じtaskに置き、見積もりの約2倍の h_rt を24時間以内に取れることを条件にした。
exact 2,000 nt は1 CLI約3.22時間で6時間以内だが、一式は約12.9時間、h_rt約25.8時間となるので除外した。exact 4,000 nt は約40.3時間/CLIで除外した。
CParty / b=100 の2,000・4,000 nt は資源上は含める。CParty の数値失敗は別の制約として記録し、時間・メモリの成功値には流用しない。

## Array と確保するポイント

全 task は cpu_16=1（35,095 MiBを計画上の予算）、1 thread/process。各長さの21入力を同じtaskに置く。
node 内の固定 worker pool、seed付き blocked_random、warm-up完了後に最初の測定、反復間 barrier は E2 と同じ runner を使う。
peak RSS は各 CLI の値。`RLIMIT_AS` は仮想アドレス空間の上限であり、RSS 上限とは異なる。
各taskは node-local `$TMPDIR` で実行し、共有側は task ごとの JSONL と tar にまとめる。checkpoint/resume は完了した失敗も保存し、未完了だけを再実行する。

| tier | tasks | workers | AS cap / process MiB | h_rt | task wall見積もり h | h_rt枠 h | 見積もり point | h_rt上限 point |
|---|---:|---:|---:|---|---:|---:|---:|---:|
| n250 | 1 | 8 | 1024 | 00:15:00 | 0.096 | 0.250 | 0.0055 | 0.0120 |
| n500 | 1 | 8 | 1024 | 00:50:00 | 0.404 | 0.833 | 0.0220 | 0.0400 |
| n1000 | 1 | 8 | 2048 | 07:40:00 | 3.827 | 7.667 | 0.2067 | 0.3680 |
| n2000 | 1 | 8 | 2048 | 01:25:00 | 0.696 | 1.417 | 0.0377 | 0.0680 |
| n4000 | 1 | 4 | 8192 | 11:50:00 | 5.897 | 11.833 | 0.3187 | 0.5680 |
| 合計 | 5 | | | | 10.920 | 22.000 | 0.5906 | **1.0560** |

課金は優先度0、`0.060 × (0.7×max(task wall秒,300) + 0.1×h_rt秒) / 3600`。上限列は wall=h_rt。
0.0001 point未満の切捨て、ストレージ、再投入は含めない。h_rtは使い切る予定時間ではなく、確保する上限の計算に使った。
[公式課金規則・別表2/3](https://www.somuka.titech.ac.jp/reiki_int/reiki_honbun/x385RG00001693.html) を2026-10-20に確認した。
task wall は、4回の実行の仕事量/worker数＋最長処理の余裕を方式ごとに足し、120秒/taskを加えた保守的な値。
固定した予測式と実際の各batchの独立した上限確認は `preparation.json` にある。並列効率の新しい実測値としては扱わない。

## 後で投入するコマンド（この作業では実行していない）

同時 array の合計が平日30 task以下になるよう、他の実行枠を差し引いた `PF_TASK_BUDGET` を使う。
以下は `-sync y` で一つずつ待つ。既知の数値失敗があるarrayの終了コードは非0になり得る。失敗を成功に書き換えない。
まず n250/n500/n1000 を確認する。長いtierの外挿が外れた場合は元のplanを変更せず、新しいrun名で再計画する。

```bash
cd /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments
: "${TSUBAME_GROUP:?Set the authorized accounting group}"
PF_TASK_BUDGET=1  # この比較は各tier 1 task。別arrayとの合計は30以下。
lfs quota -p 102571 /gs/bs
qsub -g "$TSUBAME_GROUP" -sync y -tc "$PF_TASK_BUDGET" jobs/cparty-linearknot-pf-89baea52-20261020/n250.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$PF_TASK_BUDGET" jobs/cparty-linearknot-pf-89baea52-20261020/n500.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$PF_TASK_BUDGET" jobs/cparty-linearknot-pf-89baea52-20261020/n1000.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$PF_TASK_BUDGET" jobs/cparty-linearknot-pf-89baea52-20261020/n2000.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$PF_TASK_BUDGET" jobs/cparty-linearknot-pf-89baea52-20261020/n4000.sh
```

停止後の再開は対応scriptに `--resume` を付け、古いprocessが停止済みであることを確認して同じtaskだけを指定する。
成功・失敗・timeout の完了結果は再実行しない。失敗を別設定で測り直す場合は新しいrunを作る。

## 集計・SI 作図

全taskが終了してから実行する。既存の集計・図ディレクトリは上書きしない。

```bash
python3 -B scripts/aggregate_array.py --plan results/runs/cparty-linearknot-pf-89baea52-20261020/array_plan.json
python3 -B scripts/analyze_pf_scaling.py \
  --plan results/runs/cparty-linearknot-pf-89baea52-20261020/array_plan.json \
  --out results/local/cparty-linearknot-pf-89baea52-20261020
python3 -B scripts/plot_pf_scaling.py \
  --summary results/local/cparty-linearknot-pf-89baea52-20261020/summary.json \
  --out results/figures/cparty_linearknot_pf_89baea52_20261020
```

`resources.tsv` は全入力の全方式、成功反復数、失敗、計画からの除外、時間/RSSのmedian・min/maxを残す。
`parity.tsv` は E0 と同じ絶対許容差（log Z 1e-9、PF energy 1e-10 kcal/mol）で CParty/exact を照合する。
図は165×66 mm、時間・RSSを配列長に対して両対数、PDF/SVG/PNG。各入力の3測定medianを取ってから入力間medianを描き、棒は入力・反復を通したmin/max。
各長さで同じ成功入力を全ての採用方式に使い、使用数/予定21本を図中に明記する。exact がある長さではPFの一致も要求する。exact を除いた長さではCParty/beamの有限値を比較し、exactとの一致は未評価と明記する。
交差点は形状×seedごとの隣接する実測長の時間/RSS比が1を挟む場合だけ、両対数補間値と区間をJSONに記録する。
範囲外への外挿や、失敗入力を使う交差点は作らない。区間が得られなければ `not_observed_in_valid_measured_range`。
特に CParty の数値失敗が長い長さに残る場合、現在のoracleから交差点を示せない可能性がある。

## 検証とファイル数

新しい通常arrayは未実行（5 taskとも `not_started`、attempt空、`run/`なし）。入力105本は採用E2とhash一致。
新 observer は採用libraryと既存E0 oracleのhashを検証してlinkした。外部の実装working treeは変更していない。
24/48/96 ntの小さい開発array（2 workers、warm-up 1・測定3、36 CLI）で packed runner → 集計 → PF一致 → 165 mm作図を通した。
その生出力は `pipeline_smoke.tar`、数値の検証は `pipeline_verification.json`。この開発図・時間は論文の結果ではない。
失敗/欠測の分母、PF一致のgate、交差区間、exact除外の4 testsと、既存arrayの10 testsも成功。

quotaは `quota_before_plan.txt` と plannerのlive gateに保存。保守的な24,000,000までの残りは2,415,389 entries。
arrayの計画～初回実行・scheduler log・atomic temporaryは89 entries。新 observer・config・文書・検証artifactを含めても全作業の枠は300 entries以内を見込む。
28,000,000のlive上限までの残りは約642万。これは予約ではないので投入直前にも再確認する。
共有側にCLIごとのファイルは作らず、生出力はtar。作成したnode-local中間ファイルは検証後に削除する。
原稿・既存の図は変更しておらず、この時点の試行値を本文・SIの実験結果には追加していない。
