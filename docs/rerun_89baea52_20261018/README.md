# 89baea52 で論文の全数値を取り直す準備

multiloop closing の修正を含む `89baea52c9435b18829ef0867a66a24be56d70f9` に対し、25個の入力・依存関係を検証した plan、57個の job script、399,000 CLI job を準備した。
array 用は54 script / 980 task、Terc と Figure 1 の local 用は3 script / 3 task。
5個の単独検証 task は980に含み、399,000という数は25個の通常 matrix の CLI job 数である。
**この担当は実験を実行せず、scheduler への投入と git commit も行っていない。**
通常matrixの全978 taskは `not_started`、attemptは空で、taskの `run/` は存在しない。
単独検証5件も未実行である。

## 数値の出典と作成物

- [全21項目の出典表](inventory.md)：実験、現在の版、入力、出力、保存先、取り直し方法、array/local、見積もり。
- [全57 script の表](jobs.md)：実在するパス、task数、資源、workers、プロセスごとの上限、h_rt、見積もり。
- [config manifest](../../configs/rerun_89baea52_20261018/matrix_manifest.json)、[job manifest](job_scripts.json)、[検証記録](verification.json)。
- [全 experiment inventory の JSON](experiment_inventory.json)。表中の保存先は原則 experiments/ からの相対パス。

| 区分 | task数 | 主な条件 |
|---|---:|---|
| E0 CParty 比較 | 21 | 同じ305入力、互換・既定・6パッチ oracle |
| E1 short / long | 106 / 694 | 24,625 / 1,438本、前回本番と同じ科学条件・tier・資源 |
| E2 全5 matrix | 65 | cpu_40、最大10 workers、8 GiB RLIMIT_AS、並列対照は w1/w10 |
| E3 NEAT1 全9 matrix | 18 | 主要6出力、空 scaffold、閾値0.5/0.9も再測定 |
| 回復200本×2 parameter、80本、rRNA3本、所属検査 | 71 | 凍結した同じ入力、所属検査も全200本をfreshに計算 |
| 列挙、decoding、loop別検証、ViennaRNA、MEA/centroid再構成 | 5 | 各1 task、cpu_16、h_rt 2時間 |
| Terc、Terc by-loop、Figure 1 | 3 | local の packed task、同じ入力・seed・要求出力 |
| 合計 | 983 | array 980、local 3 |

scaffold の生成器と保存済み G、配列選択、parameter file は保持する。
配列長、入力件数、G の対数、定数表などはその凍結資料から確認する。
旧 pilot、旧 cap/ranking の選択実験、Spark 開発診断、旧 sampling 検証は、現行の数値の出典でないものまで再実行しない。
現行のsampling数値は729入力を使った e2cecbb4 の検証から取り直す。
削除済みの SHAPE・ウイルス難例は対象に戻さない。
既存 LinearPartition scaffold と補助的な LinearPartition 出力は、この LinearCParty の重み修正の影響を受けない固定資料である。

## 完了したビルドと oracle

Claude による [release receipt](../../manifests/release-preparation-89baea52.json) の suite・lock・comparison 完了を確認してから提供 build を使用した。
比較は14入力 / 126条件 / 252 processが正常終了し、byte一致は34条件だった。
この比較の最大差は種類の異なる出力行を含むため、確率や自由エネルギーの最大誤差とは解釈しない。
著者の「全実験を取り直す」という判断を維持する。

| 作成したもの | パス | lock |
|---|---|---|
| 6パッチ oracle | `build/cparty-oracle-800e95a-6patches` | `manifests/software-cparty-oracle-800e95a-6patches-20261018.json` |
| inside-only E0 observer | `build/e0-observers-89baea52-800e95a-6patches` | `manifests/software-e0-89baea52-20261018.json` |
| stage driver | `build/stage-timing-89baea52` | `manifests/software-stage-timing-89baea52-20261018.json` |
| 到達可能性 observer | `build/membership-89baea52` | `manifests/software-membership-89baea52-20261018.json` |
| DP09 用の同じ CLI build | `build/linear-cparty-89baea52` | `manifests/software-pilot-89baea52-dp09.json` |
| Suboptを含む例用の同じ CLI build | `build/linear-cparty-89baea52` | `manifests/software-examples-89baea52-20261018.json` |

oracle は upstream `800e95a` に `0001,0002,0006,0007,0008,0009` の順でpatchを適用したもの。
completed-WMB closure patch は含まない。
E0 の互換モードと比較し、BPP/MFE の oracle 一致は主張しない。
[第6パッチのコピー](0009-Multiloop-closing-use-inner-neighbors.patch) は固定版の同名ファイルと同じSHA256を持つ。
旧5パッチ版は変更していない。

再構築時は新しい出力先を使い、次の手順を再現する（この作業では実行済み）。
実装 clone は prepare_release が作成したものを使用する。

```bash
cd /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments
export PYTHONDONTWRITEBYTECODE=1
: "${TMPDIR:?node-local TMPDIR is required}"

# 先に file budget を確認し、fresh な oracle 出力先を用意する。
git clone -q --no-hardlinks build/cparty-oracle-800e95a-5patches/src \
  build/cparty-oracle-800e95a-6patches/src
bash deps/LinearCParty-89baea52/test/oracle/legacy_cparty/build.sh \
  "$PWD/build/cparty-oracle-800e95a-6patches" \
  > "$TMPDIR/oracle-build.log" 2>&1
cp "$TMPDIR/oracle-build.log" build/cparty-oracle-800e95a-6patches/build.log
python3 -B scripts/build_e0_tools_89baea52.py
python3 -B scripts/build_stage_timing_driver.py \
  --software-root deps/LinearCParty-89baea52 --build-dir build/linear-cparty-89baea52 \
  --base-lock manifests/software-pilot-89baea52.json \
  --out build/stage-timing-89baea52 \
  --lock manifests/software-stage-timing-89baea52-20261018.json --purpose formal
python3 -B scripts/build_membership_rerun_89baea52.py
python3 -B scripts/prepare_paper_rerun_89baea52.py
python3 -B scripts/prepare_paper_rerun_89baea52.py --plan-existing
```

E2 config は既存 `make_rerun_configs.py` と d02d261c の cpu_40 settings を使って生成した。
E1 は通常のコピーで前回 tier と条件を完全に保ち、全 matrix の plan は `plan_array.py` の `prepare()` で作成した。
通常 script は同じ node-local 実行・task JSONL/tar・checkpoint/resume を使用する。
local の例3個は `purpose: development` として、scheduler ID を捏造せず実行できる。
版と入力は正式な実験と同様に固定し、その local elapsed time は論文の性能測定には使わない。

E0 の旧道具は削除済みAPIを使っていたため、新版用コピーだけでRT loader、runtime owner、削除済みの root-dead VP alias を更新した。
列挙器の入力・整数演算・context割当表・許容差、decodingとsamplingの検定条件は維持した。
旧729入力と整数比較に使う4つの集約資料は、`tools/validation_89baea52/baseline_15aa7742/` に同じbyteのコピーを固定した。
新しい数値は計算し直し、この資料から旧エネルギー値を流用しない。
最終確認で報告書writerのpath変数の再利用によるrunnerの上書きを検出し、復元後に両runnerのimportと全依存hashを確認した。
[差分の根拠とhash](../../tools/validation_89baea52/source_adaptations.json) に記録した。
最初のE0 link失敗と2回のplanner設定エラーも [準備ログのtar](preparation_logs.tar) に残し、完了したplanは再検証して再利用した。

## 投入コマンド

**下のコマンドは実行していない。Claudeが投入するときに使う。**
会計groupはこの依頼と既存READMEで特定されていないため、認可済みの値を `TSUBAME_GROUP` に設定する。
他のarrayが走っている場合、この一式に使える枠を `LCP_ARRAY_BUDGET=30−他の実行枠` とする。

```bash
cd /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper
# TSUBAME_GROUP は Claude が使用している認可済みgroupを設定する。
: "${TSUBAME_GROUP:?Set the authorized accounting group}"
export LCP_ARRAY_BUDGET=30  # 他のarrayがある場合は減らす。
bash experiments/docs/rerun_89baea52_20261018/submit_arrays.sh
```

[submit_arrays.sh](submit_arrays.sh) は54 scriptを `qsub -sync y` で投入する。
標準設定は最大3つの同期qsub、各 `-tc 10` なので、この一式のactive taskは30以下。
1つが終了したら次を投入する。
資源とh_rtは保存scriptの指定を使い、上書きしない。
個別のqsub行を見たい場合は [逐次投入コマンド一覧](qsub_commands_sequential.sh) を参照する。
こちらは各arrayを `-sync y -tc 30` で順番に待つ代替手順であり、両方を実行しない。
各scriptの `-tc 30` だけでは複数arrayの合計を制限できない。

小さい例と単独検証のlocalコマンドも用意した。

```bash
# Tercの12条件、by-loopの2条件、Figure 1の6条件。逐次、同じpacked runner。
bash experiments/docs/rerun_89baea52_20261018/run_local_examples.sh

# 単独検証をarrayではなくlocalで行う場合は、対応arrayを投入しない。
# それぞれ前回18–26分程度なので、推奨は準備済みの1-task array。
python3 -B experiments/scripts/run_paper_validation_89baea52.py enumeration
python3 -B experiments/scripts/run_paper_validation_89baea52.py decoding
python3 -B experiments/scripts/run_paper_validation_89baea52.py byloop
python3 -B experiments/scripts/run_paper_validation_89baea52.py vienna
python3 -B experiments/scripts/run_paper_validation_89baea52.py linearspace
```

停止・timeout の再開は既存の task checkpoint を使う。
生きている旧processがないことを確認し、該当taskのscriptを `-t N-N ... --resume` で再投入する。
E1 long は既定でb=20/50/100/200/300/400を計画順に実行する。
途中段階だけ実行する場合は既存 `--through-beam` と `--resume` を使い、未実行幅を欠測として残す。

## 集計と原稿への反映

全条件の完了後、[analyze.sh](analyze.sh) に集約・集計コマンドをまとめた。
E1のrootはcounts-onlyで集約し、解析はtaskのtarとJSONLを直接読む。
shortのbeam参照とexact参照には同じ新planを指定し、exactを二重に測定しない。
層別の `--previous` は新しいE1集計であり、旧版の誤差表ではない。

```bash
bash experiments/docs/rerun_89baea52_20261018/analyze.sh
```

到達可能性は旧版BPPへのfallbackを使わず、全200本の新observer結果を要求する。
900秒上限に達した入力は欠測のまま残す。
通常arrayは失敗もcheckpoint済みの観測なので、`--resume` が失敗済み条件を自動再計算するわけではない。
失敗を再測定するときは、保存した失敗を維持した新config/planを作る。

ViennaRNA照合は保存済み337入力と548構造条件を再利用し、数値を新たに計算する。
dangles=2は実際の `ExactSession`、dangles=0の73条件は既存内部controlである。
全2,226,060 BPP座標を比較し、MFE候補も独立なsingle-structure PFで再採点する。
integer MFE/eval energyとの等値は要求しない。

TercとFigure 1の出力は単一の `*_data.tar` にexportし、`terc_summary.json` / `figure1_summary.json` に数値を再集計する。
sampling件数、表示cutoffを通る対数、PseudoLoop表示区間も再計算する。
既存の図を描き直す際は、描画scriptとcaptionに残る旧版名・旧source READMEの表記も更新する。
E2/NEAT1の新版と旧版のbyte一致を前提にする旧差替え手順は使わず、ここに指定した集計から本文・SI・図を更新する。
モデル定義のclosing-neighbor表、oracle patch数・版表、図captionの数値も更新対象である。

## 時間と共有ファイル数

時間は新版のpilot値ではなく、前回の費用モデルと小規模解析の保守的なallowanceである。
全983 taskの逐次仕事量は約4,211時間、plan上のh_rt枠合計は9,252 task-hours（local例の12時間枠を除く投入arrayは9,240 task-hours）、最長の計画taskは10.98時間。
大半はE1 longである。
4,211時間を30枠に均等配分した約140時間は、計画量の比較用の値であり、queue待ち・array境界・task内並列・timeoutを含む完了予測ではない。
細かい資源・時間は [jobs.md](jobs.md) にある。

前回の同じtierについて、controllerの全attemptの開始・終了時刻から足した実績も併記する。
これは既にtask内の並列を含む割当時間であり、上の逐次CLI仕事量とは異なる。
途中段階・再開・失敗attemptも含め、終点が不明のattemptはない。

| 前回の範囲 | task数 | 記録されたtask wall合計 h |
|---|---:|---:|
| E0 | 21 | 0.39 |
| E1 short | 106 | 26.31 |
| E1 long | 694 | 536.80 |
| E2 | 65 | 23.18 |
| E3（変種を含む） | 18 | 11.62 |
| 合計 | 904 | 598.30 |

[task別の元metadata hashとtier別実績](historical_wall_estimates.json)を保存した。
598.3 task-hoursを30枠へ均等配分した約20時間は、前回実績からみた参考規模である。
新しくarray化した回復・所属・単独検証、queue、端数task、3 array×10枠の配分による待ち、新しい重みによる実行時間の変化は別に必要になる。
E1 shortの旧exactはd21618c5で、新版のexact経路への速度補正は行っていない。


E2は `10 × 8,192 + 1,024 = 82,944 MiB` のcap envelopeを `cpu_40` の87,738 MiB内に収める。
E1は前回の2 workers/cpu_16と1 worker/cpu_40のtier別capを保持した。
NEAT1は通常3×20 GiB、sampling2×40 GiBである。
RLIMIT_ASはaddress-space上限でありRSS上限ではない。

| inode見積もり | entries |
|---|---:|
| 25 matrixの計画・実行・初回scheduler log・atomic temporary | 13,168 |
| このうち現在すでにあるplan/jobの実entries | 4,316 |
| 単独検証5個のcontroller・結果・logの余裕 | 200 |
| config・文書・集計用の余裕 | 300 |
| oracle・observerのビルド用の事前allowance | 1,500 |
| 全期間の保守的見積もり | 15,168 |
| 現在から実行・集計で加わる保守枠 | 9,352（約9,400） |

ビルド前はoracle既存版の431 entriesを数え、1,500 entries枠を [live quotaで確認](build_file_budget.json) した。
完成した6パッチoracleは433 entries。
matrix生成前には全体20,000 entries枠を [別に確認](matrix_file_budget.json) した。
liveは28,000,000 files上限、plannerは事故時の24,000,000を保守的上限として使用する。
[最終quota確認](final_file_budget.json)では保守的上限まで2,448,339 filesが残り、今回の全期間見積もりはその1%未満。
このcheckは予約ではないので、投入直前の `lfs quota -p 102571 /gs/bs` も再確認する。
再投入のscheduler logは1 attempt/taskごとに追加される。

小規模実験用の既存285入力は1つの `data/rerun-89baea52-20261018/inputs.tar` にまとめた。
全CLIの一時ファイルはnode-local TMPDIR、通常taskの共有側出力は既存のJSONL/tar方式である。
原稿用の新しいCLIごとの共有ディレクトリは作らない。
旧採用runのtar量からbyte量は数十GiB規模と見込み、80–100 GiB程度の共有byte枠を想定するが、新版の出力サイズの上限を保証する値ではない。
各実行nodeのTMPDIR容量も出力量と同時worker数に合わせる。

## 開発上の検証

31 tests成功、実際のRNA CLIを実行するstage tests 9件は明示的にskipした（[log](preparation_tests.log)）。
E0 observer、stage driver、所属observerを実際にlinkした。
列挙・decoding・ViennaRNA・loop別probeのC++構文確認と、57 job scriptとrunbook shellのbash構文確認も成功した。
提供suiteのlogには依存gateの `rg` が見つからない行があったため、そのgateを利用可能なrgで単独再実行し、成功を保存した。
E1の科学条件とtier資源の一致、全planの入力hash・依存lock・未実行状態を確認した。
この検証は新しい科学的な一致・精度・計算時間の測定を代替しない。
原稿・図の数値と `../LinearCParty` は変更していない。


## E1 long の cpu_16 再計画（2026-10-18、未投入）

E1 long は h_rt 全体に対するポイント確保を減らすため、新run `e1-bprna-long-89baea52-cpu16-20261018` に再計画した。
旧runと旧job scriptは変更せず、入力1,438本・全25,884条件・mode/program/要求出力・順序方針・seedを維持した。
全区分を cpu_16、8 workers、4,096 MiB/process、task余裕1,024 MiBとし、cap合計33,792 MiBを35,095 MiBの予算内に収めた。
4,381 nt・b=400のRaccessは4 GiB制限で完走（595.17秒、GNU time peak RSS 1.60 GiB）。
BPPの180/600秒試行と3,995 ntのRaccessの900秒試行は時間で打ち切った。打ち切りはメモリ不足を示さず、BPPと3,995 nt入力の完走時の最大メモリは未確認である。
試行はこの2入力に限り、全attemptの生出力を [pilots.tar](cpu16/pilots.tar) に保存した。

[比較表と試行一覧](jobs.md#e1-long-cpu_16-再計画2026-10-18)、[較正](cpu16/calibration.json)、[計画と不変性検証](cpu16/summary.json) を参照。
旧CLI実測536.03時間と完走試行の新旧比 0.516349 を使い、各job 0.5秒追加・仕事量1.5倍、並列効率80%の仮定、幅ごとの待ち合わせ、120秒/taskを含めた。
並列効率は8 workerの実測値ではない。
h_rt は1/2/3/4時間、全taskの見積もりは各h_rtの70%以内。
**694 → 60 task、予約枠は旧cpu_40の7,088時間から新cpu_16の157時間。資源が異なるので、この時間比をポイント比とはしない。**
新planと付帯資料の全期間inode枠は852 entriesで、[live quota gate](cpu16/file_budget_before.json) を通過した。

新しいjob scriptは次の4個。

- `jobs/e1-bprna-long-89baea52-cpu16-20261018/nt400_799.sh`
- `jobs/e1-bprna-long-89baea52-cpu16-20261018/nt800_1599.sh`
- `jobs/e1-bprna-long-89baea52-cpu16-20261018/nt1600_3199.sh`
- `jobs/e1-bprna-long-89baea52-cpu16-20261018/nt3200_4381.sh`

**この作業ではqsubもgit commitもしていない。以下は後で投入するためのコマンド。**
既存の `submit_arrays.sh` と `qsub_commands_sequential.sh` は旧E1 longを参照するため、E1 longにそのまま使わない。
新helperは4 arrayを `qsub -sync y` で順番に待ち、他arrayを除いた `LCP_ARRAY_BUDGET` を各arrayの `-tc` に渡す。
既存 `analyze.sh` で新runを集計するときも、E1 longのplan指定をこの新runに置き換える。

```bash
cd /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper
: "${TSUBAME_GROUP:?Set the authorized accounting group}"
export TSUBAME_GROUP
lfs quota -p 102571 /gs/bs
export LCP_ARRAY_BUDGET=30  # 他のarrayに使う枠があれば減らす。
bash experiments/docs/rerun_89baea52_20261018/cpu16/submit_e1_long_cpu16.sh
```

## E1 を除く集計と図の差し替え（2026-10-18）

`analyze_non_e1.sh` は上の `analyze.sh` から E1 short/long と E1 に依存する層別精度を除いたもの。
完了した E0・E2・回復・到達可能性・小例を先に集計し、最後に NEAT1 の未完了 task を待つ。
各 plan の全 task が success になってから集約する。
新しい集計の出力先は元の `analyze.sh` と同じで、旧版の科学的測定値を代用しない。

```bash
bash docs/rerun_89baea52_20261018/analyze_non_e1.sh
python3 -B scripts/render_non_e1_89baea52.py
python3 -B scripts/audit_non_e1_89baea52.py
```

描画は `results/figures/` の `figure1_89baea52_20261018/`、
`figure4_terc_89baea52_20261018/`、`main_89baea52_20261018/`、
`figure_s3_neat1_89baea52_20261018/` に保存する。
Terc と Figure 1 の tar は TMPDIR だけに展開し、図の再現用 manifest に tar と各入力の hash を記録する。
E2 の図だけを先に描く場合も、新しい E2 集計と `--e2-revision 89baea52 --figure3-only` を指定する。

ViennaRNA の初回照合では、本番の dangles=2 は完了したが、内部 dangles=0 対照 73 条件が観測プログラムの終了時に失敗した。
原因は、別名の関数として呼び出した旧 `main` に明示的な `return 0` がなかったことだった。
`energy_vienna_internal_89baea52.cc` にその return を加え、公開側の観測プログラムから参照するようにした。
実装の固定ソースは変更していない。
元の失敗した controller と `energy-vienna-check-89baea52-20261018/` を残し、次の別出力先で全条件を計算し直した。

```bash
python3 -B scripts/validate_vienna_89baea52.py \
  --out results/local/energy-vienna-check-89baea52-repaired-20261018
```

修正後は337配列、2,226,060 BPP座標、548構造条件（公開dangles=2が475、内部dangles=0が73）すべて成功。
最大絶対差は log Z が1.4495071809506044e-12、BPPが1.245337166722038e-12、
単一構造のPFエネルギーとMFE再採点が2.842170943040401e-14 kcal/molだった。
元の失敗を成功として書き換えていない。
