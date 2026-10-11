# bpRNA 全体の G・E1 と CParty 比較 E0 の array 準備

2026-10-11。qsub、commit、push は実行していない。
実装の作業 tree と論文 repo のファイルは変更していない。
計画・job script・係数合わせ・テストと、52 nt の G 生成および4証人の小さな確認だけを行った。
本実験の結果ではない。

## 入力と条件

[取得 manifest](../manifests/data-sources/bprna-local-20261008.json) の clean 集合を全件使用する。
26,063 本、50–4,381 nt、ACGU-only、配列重複なし、FASTA SHA-256 は
`c76025522f1f9e174c25207595bc079e454deb656acbfd3a51de71c60c167bc5`。
保存した [全件 cohort](../manifests/cohorts/bprna-full-20261011.json) は長さ・ID順で、ID・長さ・配列 hash を持つ。
成功しやすい入力への差し替えや、G の失敗による分母の削減を行わない。

G は LinearPartition-V `b450fb3e`、beam200、dangles2、posterior ≥0.7、stacked helix 長 ≥2。
26,063 本すべてで生成し、同じ保存済み G を E1 と E0 の bpRNA 部分集合へ渡す。
[G config](../configs/scaffold_bprna_full_20261011.json) と
[LP・生成コードの lock](../manifests/software-scaffold-b450fb3e-20261011.json) を用意した。
`scaffold_array.py` が tier 別の計画、SGE_TASK_ID による実行、明示的な集計を行う。
各 task の `dataset.json` は既存の生成 manifest の形を引き継ぎ、コマンド、elapsed、GNU time の wall/RSS、
入力・配列・G・wrapper・binary・生成コードの hash、終了 status を保存する。
中断 task は `--resume` で再開でき、完了した failure/timeout と completion marker を再実行しない。
G が全件成功するまで、完成 dataset の公開と E1/E0 の通常 planning を拒否する。

E1 は固定した LinearCParty `15aa7742` と3 CLI を使う。
`--exact` は n<400 の24,625 本のみ（73,875 CLI jobs）。
全26,063 本に b=20/50/100/200 × 局所上限30/30・100000を適用する（625,512 CLI jobs）。
LinearCParty は cutoff0 の BPP と MFE、CapR は7 context、Raccess は probability の w=1/4/8/16。
Warm-upなし、測定1回。短い tier の2–4並列は精度測定用で、時間を単独性能測定へ転用しない。
[E1 config](../configs/e1_bprna_full_15aa7742_20261011.json) と
[array config](../configs/arrays/e1_bprna_full_20261011.json) に全条件を保存した。

## G と E1 の tier

「時間」は task に割り当てた仕事の逐次和で、安全係数を含む。並列数で割っていない。
全 task を1個ずつ実行した場合の合計であり、queue 待ちを含む暦時間ではない。

| 長さ nt | RNA 数 | G tasks / jobs | G 時間 h | E1 tasks | E1 CLI jobs | E1 時間 h | E1 resource / workers / h_rt |
|---|---:|---:|---:|---:|---:|---:|---|
| 50–99 | 12,428 | 25 / 12,428 | 0.74 | 25 | 335,556 | 14.46 | cpu_4 / 4 / 2 h |
| 100–199 | 9,656 | 20 / 9,656 | 0.85 | 35 | 260,712 | 46.96 | cpu_8 / 4 / 3 h |
| 200–299 | 1,446 | 3 / 1,446 | 0.20 | 10 | 39,042 | 11.03 | cpu_16 / 2 / 3 h |
| 300–399 | 1,095 | 3 / 1,095 | 0.21 | 38 | 29,565 | 90.98 | cpu_40 / 2 / 5 h |
| 400–799 | 584 | 6 / 584 | 0.87 | 182 | 14,016 | 312.14 | cpu_16 / 1 / 3 h |
| 800–1,599 | 423 | 7 / 423 | 1.64 | 186 | 10,152 | 590.97 | cpu_40 / 1 / 6 h |
| 1,600–3,199 | 340 | 11 / 340 | 2.51 | 274 | 8,160 | 904.96 | cpu_40 / 1 / 8 h |
| 3,200–4,381 | 91 | 4 / 91 | 0.96 | 91 | 2,184 | 346.51 | cpu_40 / 1 / 10 h |
| 合計 | 26,063 | 79 / 26,063 | **7.99** | **841** | **699,387** | **2,318.01** | 8 arraysずつ |

G は全 tier で cpu_4、worker1、h_rt=1 h、120秒/job、8,192 MiB address-space cap。
Controller 込みの cap envelope 8,704 MiB は、公式の9.2 GBをdecimalとして換算した資源枠内。
短い G tasks は最大500本、長いものは費用に応じて20–100本程度。
E1 は50–99 ntで428–500本/task、100–199 ntで90–300本/task。
200–299 ntは96–150本、300–399 ntは9–30本、長い tier は1–4本/task。
同じ RNA の全 mode と全 CLI は一つの task に保つ。
E1 の per-job memory cap は順に1/2/8/16/16/32/64/64 GiB。
Timeout は n<400 が600秒、それ以外7,200秒。
どちらも保存済み tier の `job_limits` から task config に明示的に展開する。

上記の長配列の時間は非常に保守的な外挿である。
短い pilot の速度を全件へ直接引き延ばした値ではない。
特に上限なしの資源量はまだ測れておらず、長い入力で小幅から段階的に確認し、実測で見積りを更新する。

## 長配列の b300/b400 予備計画

[選定 manifest](../manifests/cohorts/bprna-long-beam-pilot-20261011.json) は各長さ区間から3本、計12本。
Seed20261011、ID順の候補から抽出し、結果を見る前に固定した。
主 E1 の同じ12本で b20→50→100→200 の完了率、時間、RSSを確認する。
次に [b300](../configs/e1_bprna_long_beam300_pilot_20261011.json)、その結果を確認してから
[b400](../configs/e1_bprna_long_beam400_pilot_20261011.json) を実行する別計画にした。
両方とも同じ局所上限2条件・3 CLI・測定1回で、1 RNA/task。
全長配列への b300/b400 拡大を自動実行する config は作っていない。

| 長さ nt | RNA | b300 tasks / jobs / 時間 h | b400 tasks / jobs / 時間 h |
|---|---:|---:|---:|
| 400–799 | 3 | 3 / 18 / 2.47 | 3 / 18 / 4.38 |
| 800–1,599 | 3 | 3 / 18 / 8.11 | 3 / 18 / 14.41 |
| 1,600–3,199 | 3 | 3 / 18 / 11.27 | 3 / 18 / 20.03 |
| 3,200–4,381 | 3 | 3 / 18 / 18.55 | 3 / 18 / 32.98 |
| 合計 | 12 | **12 / 72 / 40.39** | **12 / 72 / 71.81** |

b300 は node_q、128 GiB cap、h_rt16 h、10,800秒/job。
b400 は node_h、192 GiB cap、h_rt24 h、18,000秒/job。
これらの大きな枠は b100 の RSS を b² で外挿した暫定 envelope による。
各 stage の実測で減らせる可能性がある。予備測定でこれらの条件を再検討してから全体への拡大を決める。
b1000 は計画に含めない。

`--through-beam 20` は保存した task の b20 までを実行して止める。
集計・点検後に `--resume --through-beam 50`、100、200へ進められる。
この中間停止は `interrupted` として残り、まだ実行していない大幅の jobs は missing のまま。
長い task の通常 schedule も幅の昇順で、大幅を小幅の終了前に開始しない。
各予備 RNA が属する主 E1 task は [task割当](full_bprna_arrays_20261011/E1/tasks.tsv) のID列から確認する。
同じ task 内の RNA を事後的に取り除かない。

長配列の config analysis は **最大の測定幅との差（beam-width sensitivity）** と定めた。
Exact誤差と呼ばず、default と nocap の参照を分け、単調な改善を仮定しない。
`analyze_accuracy.py --compare-runs ...` は同一版・モデル・S・Gを確認して幅を結合し、RNAごと・上限条件ごとに
最大の計画幅を参照にする。最大幅が失敗した場合は `reference_unavailable` にして、小幅へ置き換えない。
主 E1 だけなら長配列の参照は b200。
全件の解析は task ごとに行い、RNA単位の結果を結合して tier 分布を作る。
既存の全出力をメモリへ読む解析を26,063本へ一括適用しない。

## E0 の oracle・入力・比較

固定 clean clone `deps/LinearCParty-15aa7742` にある `test/oracle/legacy_cparty/build.sh <out_dir>` を使用し、
`build/legacy-cparty-59593650-7patches/` に upstream `59593650` + 7 patches の oracle を build した。
ローカルにある upstream の Git objectsを読むだけでcloneし、再利用するmirrorを `deps/CParty-upstream-59593650.git` に置いた。
出力の `src` がそのmirrorをoriginに持つので、以後のbuildも調整中の実装treeを変更しない。
[Oracle lock](../manifests/software-cparty-oracle-59593650-7patches-20261011.json) にsource、patches、cache、parameter、CParty・wrapperのhashとbuild結果を保存した。
Buildは成功。完成WMBのweak closureはこのoracleへ追加していない。

互換CLIフラグがないため、[外部 API observer](../tools/e0_linear.cc) を作り、固定buildのライブラリにリンクした。
`E0LegacyCompatible` と `E0Default` がそれぞれ `LegacyCompatible` と `RequireWeakClosure` を選ぶ。
Log Z、ensemble free energy、BPP、MFE、7 contextをfull precisionで出力する。
上流CLIのBPPはsampling由来なので、その出力をexact BPPとみなさない。
[Oracle observer](../tools/e0_legacy.cc) は外部の観測用sourceコピーでPFのscalarをreverse ADへ置き換え、
V/VPのpair-emissionにlog multiplierを付け、その微分からBPPを得る。
Oracle本体のsource・7 patches・生成する導出は変更していない。
元のoracle objectsへリンクした [非微分 observer](../tools/e0_legacy_scalar.cc) のroot値とも照合した。
[Build手順](../scripts/build_e0_tools.py) と [observer lock](../manifests/software-e0-15aa7742-20261011.json) は
ライブラリ・道具・観測コピー・oracle sourceとbinaryのhashも保持する。

bpRNAは50–99/100–199/200–299/300–399 ntから各75本、計300本をseed20261011で抽出した。
[選定一覧](../manifests/cohorts/bprna-e0-subset-20261011.json) に全ID・hashを保存した。
[証人4本](../data/e0-witnesses-20261011/dataset.json) はclosure14、nested45、BE41、BE65。
論文 repo の保存例を読むだけでコピーし、BE65のforced-unpaired `X` も両実装で保存した。
同じparameter file、37°C、dangles2、exactで3 observersを比較する。
独立したanalytic profileのoracle APIはないので、profileの外部exact一致は認定しない。

| E0 tier nt | 入力数 | tasks | CLI jobs | 時間 h | resource / h_rt |
|---|---:|---:|---:|---:|---|
| 1–49（証人3本） | 3 | 1 | 9 | 0.005 | cpu_4 / 8 h |
| 50–99（bpRNA75 + BE65） | 76 | 2 | 228 | 0.12 | cpu_4 / 8 h |
| 100–199 | 75 | 2 | 225 | 0.83 | cpu_8 / 8 h |
| 200–299 | 75 | 8 | 225 | 3.16 | cpu_40 / 8 h |
| 300–399 | 75 | 8 | 225 | 19.99 | node_q / 8 h |
| 合計 | **304** | **21** | **912** | **24.11** | 5 arrays |

E0のAD道具は12本E1のexact費用に時間×10、RSS係数×6を仮置きした。
小証人は数値・出力の確認に限り、300 nt級の資源モデルを認定する測定ではない。
短いtierから時間とメモリを確かめ、必要なら新しい計画版で改定する。

小さな確認は12/12 CLI jobs成功。
互換設定対oracleの最大 |Δlog Z| は2.14×10⁻¹⁴、BPP最大差は1.25×10⁻¹⁴。
非微分observerとAD observerのlog Z・free energy・MFE・RTは4本とも同値。
45 ntではdefault対互換のΔlog Zは約−2.33×10⁻⁷、BPP最大差は約2.33×10⁻⁷で、closureの効果を観測した。
**BE41のMFEは互換設定でもoracleと0.0001706856 kcal/mol異なる。**
事前の厳しい許容差1×10⁻⁷ではparity failureとして保存し、許容差を事後的に広げていない。
整数エネルギーとPF重み由来の長loop値の丸め差が候補だが、原因は実装担当に要確認。
`analyze_e0.py` はdefault−oracleだけでなく、default−互換と互換−oracleを分けて保存し、
互換一致が確認できない量の差をclosureだけへ帰属しない。
[確認結果](full_bprna_arrays_20261011/e0_witness_analysis.json) と
[元oracleとのroot照合](full_bprna_arrays_20261011/e0_uninstrumented_checks.json) を参照する。

## 係数合わせとポイント

[Calibration manifest](../manifests/calibration/bprna-array-20261011.json) はE1の324ジョブのwall/RSSと、
RNAcentral `58a2b796` の9ジョブ（8成功・1 timeout）を、元run hashとともに保存する。
短いtier・modeごとに、3 CLIを合わせた `max((wall−0.1)/n³)` の正部分を時間係数にした。
観測のないtierは12本の同mode全体へ戻る。
RSSは `128 + 1.5*max((RSS_MiB−128)/n³)*n³`。
Task時間は全CLI・全modeの逐次和に安全係数1.5を掛ける。
n³はpackingの便宜的な式で、exactの理論計算量の主張ではない。

長配列は成功例の `max(wall/n)` =0.0531941秒/nt（b100）を使い、固定capの O(nb²) を参考に `(b/100)²` で外挿する。
RSSはb100の最大1.82086 MiB/ntに1.5の余裕を掛け、同じb²のplanning envelopeを使った。
RSSの理論量を認定する式ではない。
上限なしの時間は暫定×2、task時間にはさらに1.5を掛けた。
上限なしの探索費用を固定capの O(nb²) が保証するとは扱わない。
RNAcentralは異なる版・node・一部併行測定なので、採用版の性能結果として混ぜず、保守的なpriorに限る。
Timeout時間を必要時間の観測にせず、失敗の情報を保存した。
Gの時間係数は既存生成の成功elapsedに基づく短長2式で、RSSは未校正の余裕とcapによる。
個々のtimeoutが多ければtask wallを超えるため、欠測と中断を保存してtaskを再開する。

課金は [2026年度の公式規則・別表2/3](https://www.somuka.titech.ac.jp/reiki_int/reiki_honbun/x385RG00001693.html) を2026-10-11に確認した。
通常優先度・1 allocationでtaskごとに
`resource係数 * (0.7*max(実行秒,300) + 0.1*h_rt秒) / 3600` を合計する。
CLI jobs数ではなくscheduler tasks数が課金単位。
cpu_4/8/16/40は0.015/0.03/0.06/0.15、node_q/hは0.25/0.5。
実測がないので、上記task時間を実行秒へ代入した概算と、全taskがh_rtまで使う枠の概算を併記する。

| 計画 | task時間 h | 概算ポイント | h_rtまで使う場合のポイント |
|---|---:|---:|---:|
| G全体 | 7.99 | 0.23 | 0.95 |
| E1全体 | 2,318.01 | 287.69 | 559.73 |
| E0 | 24.11 | 6.49 | 21.15 |
| b300予備 | 40.39 | 11.87 | 38.40 |
| b400予備 | 71.81 | 39.53 | 115.20 |
| 合計 | 2,462.31 | **345.81** | **735.43** |

Queue待ち、再実行、優先度追加、storageを含まない。料金の0.0001ポイント単位切捨ては概算では省いた。
利用group・残額・同時task数・queueのresource/h_rt制限・module環境は投入前に要確認。
E1は標準のraw配置で約373万ファイル＋約70万jobディレクトリになるため、容量とi-node quotaも確認する。
短いCLIを一つずつarrayへ投げず、数百本のmatrixをbatchすることでscheduler task数と最低課金の影響を抑える。

## Dry-run artifacts と実行前の再計画

[機械可読 summary](full_bprna_arrays_20261011/summary.json) に全tier・task/job数・時間・ポイント・plan/script hashを保存した。
各サブディレクトリにtask割当TSVと29本のscriptのcopyがある。
実際の凍結plan、入力copy、task config、seed付き全scheduleは以下のignored領域。

```text
results/runs/g-bprna-full-plan-20261011-v3/
results/runs/e1-bprna-full-dry-20261011-v2/
results/runs/e0-bprna-legacy-dry-20261011-v2/
results/runs/e1-long-beam300-pilot-dry-20261011/
results/runs/e1-long-beam400-pilot-dry-20261011/
jobs/<対応するrun名>/<tier>.sh
```

G以外は `pending_scaffolds=true` のdry-run。仮のdotsは長さによるpacking専用で、計算への使用をrunnerが拒否する。
既存12本のGを全体datasetの代わりに流用していない。
E1の大きなscheduleはtaskごとに分割保存し、dry-runの集計は `summary_only=true`。
699,387件すべてがmissingで、完了jobは0。
生成scriptは全て `bash -n` を通し、qsubは実行していない。

Gの各arrayを将来実行した後、experiment repoから次を行う。下記も今回は実行していない。

```bash
python3 scripts/scaffold_array.py aggregate \
  --plan results/runs/g-bprna-full-plan-20261011-v3/array_plan.json \
  --dataset-out results/local/bprna-full-scaffolds-20261011

python3 scripts/plan_array.py \
  --config configs/e1_bprna_full_15aa7742_20261011.json \
  --array-config configs/arrays/e1_bprna_full_20261011.json \
  --software-root deps/LinearCParty-15aa7742 --build-dir build/linear-cparty-15aa7742 \
  --out results/runs/e1-bprna-full-ready-20261011

python3 scripts/plan_array.py \
  --config configs/e0_bprna_legacy_20261011.json \
  --array-config configs/arrays/e0_bprna_legacy_20261011.json \
  --software-root deps/LinearCParty-15aa7742 --build-dir build/e0-observers-15aa7742-v2 \
  --out results/runs/e0-bprna-legacy-ready-20261011
```

同じ方法でb300/b400のconfigとarray configを使い、新しいready planを作る。
G全件のstatus・ID・S・G hash・LP版・設定が一致しない場合は再計画を拒否する。
長配列を段階的に進める際はgenerated scriptへ `--through-beam` と、再開時の `--resume` を渡す。
E0のMFE差を確認し、上限なしの長配列の資源見積りを実測で更新してから結果を解釈する。

## 検証

`python3 -m unittest discover -s tests`：36テスト成功。
399/400のexact境界、modeごとの係数、tierのcap、分割schedule、G未生成の実行拒否、
Gのcompletion marker回復、変更されたGの拒否、幅ごとの停止・再開・失敗保持・journal回復、
上限条件別の参照と最大幅失敗時の欠測、ADの独立した2構造ensemble、forced-unpaired入力を確認した。
52 ntのGは既存pilotとinput/G hashが一致した
（[小確認](full_bprna_arrays_20261011/g_one_record_smoke.json)）。
