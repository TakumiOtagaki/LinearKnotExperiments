# LinearKnot 4f97fd3 の固定・出力照合・時間とメモリの再測定準備

**固定と照合を完了し、12 matrix・82 task・30 job scriptを準備した。qsub と git commit は行っていない。**
全82 taskは `not_started`、attemptは空、taskの `run/` は存在しない。
外部の実装repoと既存runは変更していない。
精度・検証・Terc・回復などの採用済み89baea52の結果は保持し、原稿・図の数値は今回更新しない。

## 固定と再開

前回セッションの専用cloneと途中のbuildを引き継いだ。
commitは `4f97fd332e793a50fe157e47767a1a75a306fd01`、sourceは `deps/LinearKnot-4f97fd3`、buildは `build/linearknot-4f97fd3`。
originは `https://github.com/TakumiOtagaki/LinearKnot.git`。
再開時にHEAD、clean state、CMakeのsource directoryと設定を確認した。
configureは前回と同じ `-DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DLINEARKNOT_ENABLE_LTO=OFF`、buildは `-j8`。
不足targetはmaintained suiteのbuild段階で補った。
source rootから `test/tools/run_required_scfg2_tests.sh` を実行し、**96項目とdependency gateが成功**した。
lock作成時のmaintained checks 3件も成功した。
一時的な `source/build` symlinkは削除済みで、sourceはcleanのまま。
`linearknot --version` は `LinearKnot 0.1.1` を返す。

- [release receipt](../../manifests/release-preparation-linearknot-4f97fd3.json)：引き継いだclone/configure記録、実際のcompile flags/link command、suite・照合の結果とhash。
- [software lock](../../manifests/software-linearknot-4f97fd3.json)：commit/tree、compiler、cache・CLI・parameter hash。
- [suite・lock作成ログ](preparation_logs.tar.gz)、[再開時の確認](resume_inspection.json)。

## 出力照合

referenceは `build/linearknot-9c7d6e1`、candidateは上の4f97fd3 build。
前回と同じ14入力、50–390 ntを凍結tarからTMPDIRに取り出し、各SHA256を検査した。
Turner04 parameter bytesも両sourceで一致した。
b=20/200は `h=d=max(b,30), ranking=outer`、exactは `--exact` のみ。
各processは1 thread、timeout180秒、RLIMIT_AS 16 GiB。
stdout・stderr・BPP fileを比較し、**`# elapsed_seconds` で始まる行だけ**を除いた。
名称や数値の置換、許容差による受理は行っていない。

| 出力 | 条件 | 完全一致 |
|---|---|---:|
| fold / BPP file | `--bpp-cutoff 0 --bpp FILE` | 42/42 |
| profile | 通常出力 | 42/42 |
| access | 窓1/4/8/16、`--probability` | 42/42 |
| access-byloop | 上記＋`--by-loop` | 42/42 |
| sample | count100、seed20261016 | 42/42 |
| sample-unique | 上記＋`--unique` | 42/42 |
| mea | gamma1 | 42/42 |
| centroid | gamma1 | 42/42 |
| subopt（前回と同じ追加対照） | top-k8 | 42/42 |
| 合計 | 14入力×3 mode×9出力 | **378/378** |

**756 processすべて正常終了。失敗・timeout・再試行は0。**
source・CLI・parameter・CMake cacheの前後一致も確認した。
[summary](comparison_summary.json) の `byte_identical=378` がstrict比較で、metadata renameの適用件数は0。
全コマンド・差の検査・hashは `results/local/linearknot-4f97fd3-comparison-20261021/comparison.json`、実行記録は `processes.jsonl`、生出力は `raw_outputs.tar`。
この14入力・指定optionでの回帰確認を受け、著者の方針どおり時間とメモリだけを取り直す。
照合中の時間・RSSは正式な性能値に使わない。

```bash
# experiments/ から。再現時は別の --out を選ぶ。
PYTHONDONTWRITEBYTECODE=1 python3 -B scripts/compare_linearknot_outputs.py \
  --candidate-source deps/LinearKnot-4f97fd3 --candidate-build build/linearknot-4f97fd3 \
  --reference-source deps/LinearKnot-9c7d6e1 --reference-build build/linearknot-9c7d6e1 \
  --reference-cli-layout unified --scratch-dir "$TMPDIR" \
  --out results/local/linearknot-4f97fd3-comparison-REPRODUCTION
```

## 再測定する条件

[全config・array config・plan・scriptの表](jobs.md)、[準備manifest](preparation.json)、[30 scriptのmanifest](job_scripts.json)、[前回task実績と出典hash](calibration.json) に全パスと根拠を保存した。
configは `configs/linearknot_4f97fd3_20261021/`、run・scriptはそれぞれ新しい4f97fd3の名前。
全matrixで入力件数・名前・hashが前回と一致し、幅、局所制限、要求出力、seed、順序、warm-up・反復数を保持した。

| 範囲 | tasks | CLI jobs（warm-up含む） | 計画point | h_rt上限point |
|---|---:|---:|---:|---:|
| E2 | 65 | 2063 | 4.9408 | 7.8300 |
| NEAT1 | 12 | 33 | 0.9198 | 1.4200 |
| PF（LinearKnot側） | 5 | 672 | 0.5204 | 0.7840 |
| 合計 | **82** | **2768** | **6.3811** | **10.0340** |

E2はRNAcentral 7本、合成対照156本、memory breakdown、w1/w10の5 matrix。
`cpu_40=1`、最大10 workers、各process 1 thread、RLIMIT_AS 8192 MiBを保持した。
cap envelopeは `10×8192+1024=82944 MiB` で87738 MiBの予算内。
RNAcentral・合成・並列対照はwarm-up1＋測定3。
memory breakdownは前回と同じ例外で、fresh processによるstageごとのwarm-up0・測定1。
inside/MFE/outside-BPP/loop-profileのouterとbaseline inside、計35 processを測る。
w1/w10は同じ3 RNA×3 CLIの9組で、w10の実際の同時実行は最大9本。
block・warm-up待ち合わせを保ち、旧版の並列比から新しい時間を補正しない。

CLIは `linearknot fold/profile/access`。
BPP cutoff0、窓1/4/8/16、probabilityなど旧CLIと同じ出力を要求する。
runnerの解析名 `program.name` を保ち、`executable: linearknot` と `subcommand` で起動先を指定する。
lockと実行前後のhash検査は実際の `linearknot` binaryを検査する。

NEAT1は本文・SI資源表と同じ閾値0.7のG、両isoform。
fold/BPP（cutoff0.001）、profile、access total/by-loop、gamma1 centroidはb=h=d=100/200/400、outer。
samplingは1,000構造・seed20261016・uniqueで、NEAT1_1のb100/200とNEAT1_2のb100だけ。
前回と同じwarm-up0・測定1、通常はcpu_40/3 workers/20480 MiB、samplingは2 workers/40960 MiB。
centroidは `--gamma 1 --stats` を保持する。
閾値0.5/0.9や空Gの科学的出力を追加で測り直す範囲には含めない。

## ObserverとCParty側の再利用

[stage lock](../../manifests/software-stage-timing-4f97fd3-20261021.json) は `build/stage-timing-4f97fd3/stage_timing_driver` と新版libraryを固定する。
[PF lock](../../manifests/software-pf-scaling-4f97fd3-20261021.json) は `build/pf-scaling-4f97fd3-20261021/PFScaling` を固定する。
PF sourceは `tools/pf_scaling_linearknot.cc`。
`tools/pf_scaling_89baea52.cc` と同じinside-only API・出力形式で、CPartyを起動する分岐は持たない。
MFE・outside・BPP・samplingを呼ばない。

PFは前回と同じ7形状×3 seed×250/500/1000/2000/4000 ntの105入力。
exactは前回の採用範囲250/500/1000 ntの63入力、b100は全105入力。
exactの2,000/4,000 ntは前回同様に除外し、開発用の高速化率だけでは範囲を広げていない。
exactは局所上限なし、b100はh=d=100・outer・L=30。
両方式ともTurner04、37°C、dangles2、`LegacyCompatible`。
本文E2の `RequireWeakClosure` とは異なる条件である。
資源は前回と同じcpu_16、250/500/1000/2000 ntが8 workers、4000 ntが4 workers、warm-up1＋測定3。

**CPartyは `results/runs/cparty-linearknot-pf-89baea52-20261020/` の観測をそのまま使い、再投入・再計算しない。**
新planに `cparty` modeはなく、旧LinearKnotの時間も流用しない。
CPartyの失敗・非有限値は残す。
`--cparty-plan` は入力hash・モデル・parameter hash・反復数を検査し、CPartyの観測だけを読む。

今回は次を実行済み。再構築するときはfreshな出力先を使う。

```bash
python3 -B scripts/build_stage_timing_driver.py \
  --software-root deps/LinearKnot-4f97fd3 --build-dir build/linearknot-4f97fd3 \
  --base-lock manifests/software-linearknot-4f97fd3.json \
  --out build/stage-timing-4f97fd3 \
  --lock manifests/software-stage-timing-4f97fd3-20261021.json --purpose formal
python3 -B scripts/build_pf_scaling_linearknot.py
python3 -B scripts/prepare_linearknot_timing_4f97fd3.py
```

## h_rtとポイント

資源・workers・AS capは同じ前回runから引き継いだ。
前回の完了taskの全attemptの開始・終了時刻を加え、同じtier・workersでの割当時間を校正に使った。
task wallを予定job数で割った最大値を割当時間のcreditとし、25%の余裕を加えた。
batchは前回の最大入力数以下。
h_rtは、校正した最大batch・前回の最長task・最長processを覆う値の2倍に240秒を加え、5分単位で切り上げた。
転送用120秒を加えた各task見積もりがh_rtの70%内に収まることを検査した。
これは同じ並列度での資源計画の仮定で、逐次CLI仕事量や新版の実測時間ではない。
開発測定の5–10倍・2.3倍という高速化率は使っていない。
process timeoutは旧上限とh_rt−60秒の小さい方で、AS capは保持した。

課金は優先度0、`資源係数×(0.7×max(task wall秒,300)+0.1×h_rt秒)/3600`。
cpu_40係数0.150、cpu_16係数0.060。
h_rt上限列はtask wall=h_rtとした確保枠である。
0.0001 point未満の切捨て、storage、再投入は含めない。
[公式課金規則の別表2–4](https://www.somuka.titech.ac.jp/reiki_int/reiki_honbun/x385RG00001693.html) を2026-10-21に確認した。
全30 tierのh_rt・前回最大task・cap・pointsは [jobs.md](jobs.md) に示す。

## 投入コマンド（未実行）

Claudeが後で投入するときのコマンド。
他arrayを除いた残りを `LCP_ARRAY_BUDGET` に設定し、合計30 task以内にする。
helperは各arrayを `-sync y` で順番に待つ。
保存した資源・h_rtは上書きしない。

```bash
cd /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments
: "${TSUBAME_GROUP:?Set the authorized accounting group}"
export TSUBAME_GROUP
export LCP_ARRAY_BUDGET=30  # 他のarrayがあれば差し引く。
lfs quota -p 102571 /gs/bs
python3 -B scripts/prepare_linearknot_timing_4f97fd3.py --verify
bash docs/linearknot_4f97fd3_20261021/submit_sequential.sh
```

[個別qsub行](qsub_commands_sequential.sh) は代替手順で、helperと両方を実行しない。
中断後は古いprocessの停止を確認し、対応scriptを `-t N-N ... --resume` で再開する。
例：`qsub -g "$TSUBAME_GROUP" -t 1-1 jobs/e2-rnacentral-4f97fd3-20261021/rna_4_8kb.sh --resume`。
完了した失敗・timeoutはcheckpointに残り、`--resume` は再計算しない。
失敗の再測定は別runで行う。

全PF task終了後の集計例も未実行。

```bash
python3 -B scripts/aggregate_array.py \
  --plan results/runs/linearknot-pf-4f97fd3-20261021/array_plan.json
python3 -B scripts/analyze_pf_scaling.py \
  --plan results/runs/linearknot-pf-4f97fd3-20261021/array_plan.json \
  --cparty-plan results/runs/cparty-linearknot-pf-89baea52-20261020/array_plan.json \
  --out results/local/cparty-linearknot-pf-4f97fd3-20261021
python3 -B scripts/plot_pf_scaling.py \
  --summary results/local/cparty-linearknot-pf-4f97fd3-20261021/summary.json \
  --out results/figures/cparty_linearknot_pf_4f97fd3_20261021
```

E2・NEAT1も終了後に各新planの `aggregate_array.py --plan ...` で更新する。
旧採用値をelapsed timeの比だけで変更せず、新しいJSONLから時間・RSSを集計する。

## 検証とinode

runner・lock・checkpoint/resume・packed storage・PF集計などの37 testsが成功した。
unified CLIの解析名ごとの結果、subcommandの位置、BPPのtar保存、実binary改変の拒否を検査した。
CParty再利用では入力・closure・parameter・反復数の不一致と重複modeを拒否する。
小入力2本でPFの旧新版4比較が完全一致し、新stage driverの10 processが成功した。
[smoke記録](observer_smoke.json) と [raw tar](observer_smoke.tar.gz) は機能確認で、新しい正式な性能測定ではない。
[最終検証記録](verification.json) にsource/build/lock・入力・未開始状態・shell構文・inode検査をまとめる。
準備中に修正した親directory作成と新規test fixtureの初期失敗も、最終成功logとともに `development_checks.tar.gz` に保持する。

再開時のlive inodeは21,596,002 / 28,000,000、保守的な24,000,000上限まで2,403,998 entriesだった。
12 planの計画・実行・初回scheduler log・atomic temporaryは合計 **1,294 entries**。
observer・config・文書・検証の余裕を含む追加5,000 entries枠は [timing gate](timing_file_budget.json) を通過した。
clone/CLI buildは前回の25,000 entries枠を引き継いだ。
[最終inode検査](final_file_budget.json) に実測entriesと全期間の保守枠を示す。
cloneは1,232、CLI buildは2,446、observer buildは各5 entriesで、新plan・jobの現在の実体は532 entries。
これらと初回実行・文書・検証を含む全期間の保守枠は **8,000 entries**、今後の実行・集計の追加余裕は2,000 entriesとした。
最終live使用は21,601,064 / 28,000,000、保守的上限まで2,398,936 entriesが残り、gateを通過した。
quota検査は予約ではないので投入直前も再確認する。
再投入のscheduler logは別に増える。
CLIごとの一時ファイルはnode-local TMPDIRで作り、共有側にはtaskごとのJSONLとtarを残す。
既存runの生ファイルを展開した共有directoryは作らない。
