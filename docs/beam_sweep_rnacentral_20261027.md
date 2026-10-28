# RNAcentral 7 本の beam 幅・局所上限感度：準備と試走（2026-10-27）

著者の 2026-10-27 の判断に従い、E2 と同じ 7 本の S/G で、b=25,50,75,100,200,300,400、h=d=max(b,30) を各 1 回計算する準備をした。**qsub、git commit、本番 49 条件の実行は行っていない。** [機械可読の計画](beam_sweep_rnacentral_20261027_plan.json)、[全 49 条件の推定 TSV](beam_sweep_rnacentral_20261027_estimates.tsv)、[最終 dry-run 検証](beam_sweep_rnacentral_20261027_validation.json) を併せて参照する。

これは **beam 幅と局所上限を一緒に変えた設定感度**である。b=400 は有限幅の比較基準であり、exact ではない。符号・単調性・収束を仮定せず、失敗・上限超えを残し、未完了条件の F/BPP を外挿しない。原稿・bib・実装 repo は編集していない。

## 入力と計算条件

- 入力は `results/runs/e2-rnacentral-4f97fd3-20261021/array_plan.json` の 7 本の凍結 S/G。各 plan の `inputs.tar` へ再凍結し、全入力 hash を元の E2 と照合した。scaffold は作り直さない。
- LinearKnot v0.2.1、commit `2101e5d86959661ba0a3da6d598f8f202e87186b`。外部 `../../LinearKnot` は読取りだけに使い、node-local に `git archive v0.2.1` で取り出して Release build。g++ と build flags、元 archive・binary・parameter の hash は lock に保存した。
- Turner2004、37°C、dangles=2、通常の `RequireWeakClosure`、`--ranking outer`。MFE は native fold と同じ別の枝刈りによる探索である。
- 出力は `fold` の MFE 構造・energy、分配関数の log Z・ensemble free energy F、`--bpp-cutoff 0` の BPP。保持した ensemble で CLI が返す全対を保存し、未出力の対のみ 0 として比較する。cutoff による切捨てを加えない。
- warm-up 0、測定 1、single thread、1 process/task。時間は MFE+PF+outside+BPP 書出しを含む fold 全体の GNU wall time、最大メモリは同一 process の GNU peak RSS。build・gzip・checkpoint・集計を CLI 時間から除く。CPU 時間の計測ではない。
- `cpu_16=1`（36.8 decimal GB）、計画メモリ 35,095 MiB、task overhead 1,024 MiB。1 task に 1 RNA × 1 幅。大きい幅を別 array にして幅ごとの timeout を設定するため、既存の「同じ RNA の全設定を同じ task に置く」方式を今回だけ分けた。allocation/host は task metadata に残る。cpu_16 は共有 node で、各条件 1 回の時間は厳密な単独性能評価ではない。

### log Z を丸めずに記録するための計測

v0.2.1 の通常の `linearknot fold` は F を小数 4 桁、MFE を 2 桁で表示し、log Z を表示しない。その値から log Z を逆算して高精度の値として扱わない。

実験側の `beam_sweep_metadata_20261027.py` が archived `fold.cc` のコピーに記録だけを加え、補助 binary `FoldMetadata` を作る。native fold と同じ viterbi/outside/inside 呼出しから `log_partition_weight`、丸め前 F、MFE を stderr の `BEAM_SWEEP_SCALARS` JSON に 17 桁で保存する。追加の DP は実行しない。元 source は変更せず、native `linearknot`、補助 binary、生成した observer source と compile command を保存した。試走 5 条件で native CLI と stdout/BPP の byte 一致を確認した。この計測用 binary を使用したことを結果の provenance に明記する。

## 推定の根拠

基準は E2 `4f97fd3` の `metrics.tsv`（`rnacentral_current`、`LinearCParty`）の b=100 の測定 3 回の中央値。時間の和は 1,597.98 秒、最大 RSS は 3.221 GiB。既存 E2 は workers=10 の併走測定なので、v0.2.1 の単独 process の実測に読み替えない。

幅の倍率は `results/local/timing-4f97fd3-20261021/neat1_resources.tsv` の `analysis=bpp`、NEAT1_1/2 の b=100/200/400 を用いた。時間比は b200/b100=2.53–2.69、b400/b100=5.88–7.13、RSS 比は 1.63–1.76、2.82–3.17。NEAT1 の BPP は cutoff=0.001 なので、cutoff=0 の RNAcentral への推定には書出し量の不確実性がある。

各 NEAT1 の 100–200、200–400 間を log-log で補間し、2 本のうち大きい倍率を採用した。b<100 は 100–200 の傾きから下方へ延長した未校正の見積もりである。RNA ごとに `T_b=5+(T_100−5)r_T(b)` 秒、`M_b=0.1+(M_100−0.1)r_M(b)` GiB とし、小さい幅での起動・固定メモリを残した。下表はこの大きい側の推定。2 本の倍率による小さい側も TSV に残した。

O(n b² log b) から、b100 の時間に `(b/100)² log(b)/log(100)` を掛けた保守的な形の比較も各条件に記録した。同じ n では b200=4.60 倍、b400=20.82 倍になる。これは漸近式の形を用いた**予算感度の別シナリオ**で、実測に校正した上界でも、h/d を同時に変える今回の runtime 保証でもない。採用した timeout は過去の実測倍率の約 2 倍から決め、これを超える遅延は観測された「上限超え」とする。

### 条件ごとの推定 wall time（秒）

| RNAcentral ID | nt | b25 | b50 | b75 | b100 | b200 | b300 | b400 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| URS000233CDDE | 4793 | 13.7 | 26.9 | 42.6 | 60.2 | 153.4 | 267.6 | 398.7 |
| URS00019626CE | 5482 | 21.9 | 47.6 | 78.2 | 112.5 | 294.0 | 516.5 | 771.8 |
| URS000239D8A4 | 6990 | 24.0 | 52.9 | 87.4 | 126.1 | 330.4 | 580.9 | 868.3 |
| URS00023C9544 | 11832 | 29.3 | 66.3 | 110.3 | 159.7 | 420.8 | 740.7 | 1108.0 |
| URS000241043C | 13652 | 51.5 | 122.4 | 206.9 | 301.5 | 802.0 | 1415.4 | 2119.4 |
| URS000238C72C | 22543 | 70.7 | 170.9 | 290.2 | 424.0 | 1131.1 | 1997.6 | 2992.4 |
| URS00023A4DB9 | 29106 | 69.1 | 166.9 | 283.4 | 413.9 | 1104.1 | 1950.0 | 2920.9 |

### 条件ごとの推定 peak RSS（GiB）

| RNAcentral ID | b25 | b50 | b75 | b100 | b200 | b300 | b400 |
|---|---:|---:|---:|---:|---:|---:|---:|
| URS000233CDDE | 0.284 | 0.401 | 0.500 | 0.590 | 0.961 | 1.316 | 1.652 |
| URS00019626CE | 0.392 | 0.575 | 0.733 | 0.875 | 1.462 | 2.023 | 2.555 |
| URS000239D8A4 | 0.426 | 0.632 | 0.808 | 0.967 | 1.623 | 2.250 | 2.846 |
| URS00023C9544 | 0.597 | 0.910 | 1.179 | 1.421 | 2.422 | 3.377 | 4.285 |
| URS000241043C | 0.906 | 1.414 | 1.849 | 2.242 | 3.864 | 5.413 | 6.886 |
| URS000238C72C | 1.208 | 1.906 | 2.504 | 3.045 | 5.275 | 7.405 | 9.430 |
| URS00023A4DB9 | 1.274 | 2.014 | 2.648 | 3.221 | 5.584 | 7.841 | 9.986 |

全 49 CLI の wall time の和は **6.46–7.37 process 時間**。実測 CPU 時間でも array の経過時間でもない。O(n b² log b) の形を用いた別シナリオでは **17.02 process 時間**になる。最大の採用 RSS 推定は 9.99 GiB。メモリの直線幅倍率を補助的に置くと b400 は最大 12.88 GiB で、下記 21 GiB AS cap に収まるという仮定で準備した。いずれも完走やメモリ使用量を保証しない。

## array、上限、ポイント

7 個の幅別 plan、各 3 長さ tier、**21 arrays / 49 tasks / 49 CLI**。4–8 kb は 3 tasks、8–20 kb と 20–30 kb は各 2 tasks/array。全て 1 worker、1 task に 1 条件。

CLI timeout は群内の最大推定時間の約 2 倍を分単位に切り上げ、最短 180 秒。h_rt はそれに staging/圧縮/保存用の 180 秒を加え、5 分単位に切り上げた。AS cap は群内最大 RSS 推定の 2 倍+0.5 GiB を GiB 単位で切り上げる。**AS cap は仮想アドレス空間の制限、推定 RSS とは別の量**である。

| b | tier | tasks | CLI timeout | h_rt/task | AS cap/process (MiB) |
|---:|---|---:|---:|---:|---:|
| 25 | rna_4_8kb | 3 | 180 s | 00:10:00 | 2,048 |
| 25 | rna_8_20kb | 2 | 180 s | 00:10:00 | 3,072 |
| 25 | rna_20_30kb | 2 | 180 s | 00:10:00 | 4,096 |
| 50 | rna_4_8kb | 3 | 180 s | 00:10:00 | 2,048 |
| 50 | rna_8_20kb | 2 | 300 s | 00:10:00 | 4,096 |
| 50 | rna_20_30kb | 2 | 360 s | 00:10:00 | 5,120 |
| 75 | rna_4_8kb | 3 | 180 s | 00:10:00 | 3,072 |
| 75 | rna_8_20kb | 2 | 420 s | 00:10:00 | 5,120 |
| 75 | rna_20_30kb | 2 | 600 s | 00:15:00 | 6,144 |
| 100 | rna_4_8kb | 3 | 300 s | 00:10:00 | 3,072 |
| 100 | rna_8_20kb | 2 | 660 s | 00:15:00 | 5,120 |
| 100 | rna_20_30kb | 2 | 900 s | 00:20:00 | 7,168 |
| 200 | rna_4_8kb | 3 | 720 s | 00:15:00 | 4,096 |
| 200 | rna_8_20kb | 2 | 1620 s | 00:30:00 | 9,216 |
| 200 | rna_20_30kb | 2 | 2280 s | 00:45:00 | 12,288 |
| 300 | rna_4_8kb | 3 | 1200 s | 00:25:00 | 5,120 |
| 300 | rna_8_20kb | 2 | 2880 s | 00:55:00 | 12,288 |
| 300 | rna_20_30kb | 2 | 4020 s | 01:10:00 | 17,408 |
| 400 | rna_4_8kb | 3 | 1740 s | 00:35:00 | 7,168 |
| 400 | rna_8_20kb | 2 | 4260 s | 01:15:00 | 15,360 |
| 400 | rna_20_30kb | 2 | 6000 s | 01:45:00 | 21,504 |

最長群の b400 は **100 分/CLI、h_rt=1:45、AS cap=21,504 MiB（21 GiB）**。時間に達したら `status=timeout` として raw/stdout/stderr を保存する。得られなかった F/TV/RSS は欠測。GNU time が完了せず RSS が取れない場合も null のままにする。allocation 全体が強制終了した場合は未 checkpoint の条件を interrupted/unfinished とし、同じ上限での再開は未完了条件だけを実行する。完了した timeout は再開でも再実行しない。

通常優先度 0、cpu_16 の係数 0.060、ポイントは `0.060 × (0.7 × max(task 実行秒,300) + 0.1 × h_rt 秒) / 3600`。**[2026 年度の公式料金表](https://www.t4.cii.isct.ac.jp/tsubame-kyodo/tariff2026)を 2026-10-27 に再確認した。** task ごとに最低 5 分を適用し、実行時間の推定には CLI+90 秒の staging/gzip/checkpoint allowance を置く。

| 予算の意味 | ポイント |
|---|---:|
| 採用した時間推定で終了（task overhead 90 秒込み） | **0.525926** |
| 全 CLI が timeout まで使用（task overhead 90 秒込み） | **0.919400** |
| 全 49 task が h_rt を使う場合 | **1.044000** |

h_rt の合計は **21.75 allocation 時間**、CLI timeout の合計は 17.45 時間。待ち時間・ストレージ・再投入・丸めによる各課金の切捨ては見積もりに含めない。少数の短い試走と build は既存 interactive allocation 内で行い、この新規 qsub 予算に含めない。

## ファイル数、保存と復元

既存の E2 の `plan_array.py`、`run_array_task.py`、`TaskStore` を流用した。正式 task で CLI ごとの共有ファイルを作らない。BPP は node-local の `bpp.tsv` に生成し、CLI の計測終了後に gzip level 6 で圧縮し、`bpp.tsv.gz` として task の `task_outputs.tar` に保存する。元 BPP と圧縮 member の両方の hash・サイズ、圧縮秒数を JSONL に記録する。tar の fsync 後に JSONL を fsync する既存 checkpoint/recovery を維持する。

- 主要 payload は **49 task × (JSONL+tar) = 98 ファイル**。
- planner の共有 entry 見積もりは 7×20+49×13=777。dirs、task config/schedule、metadata、lock、5 run files、scheduler log、atomic temporary を含む。
- build bundle/source/log、config、docs、集計出力、試走・検証・再投入 log の余裕を加え、**全体 1,050 shared entries** を予算とした。inode gate は live quota と 24,000,000 の保守境界、残数の 10% を使い、fallback は設定しない。
- 最終確認時：使用 20,735,133 / live limit 28,000,000、保守境界まで 3,264,867、10% gate 326,486。1,050 entries は gate 内。投入直前に再検査する。

既存 b100 の 7 本の BPP は計 124,238,728 bytes。サイズが幅に比例する仮定では今回 11.5 倍で、**非圧縮約 1.33 GiB**。これは上限ではなく、低確率対の数で変わる。試走の gzip は非圧縮の 37–38% だったが、全幅には保証しない。共有には非圧縮を送らず、node-local の scratch は source/runtime と生 BPP・gzip の同時保持を見込んで少なくとも数 GiB の空きを確保する。再生成可能な build tree は検証後に削除し、source/runtime/build metadata は tar にまとめて保存する。

source は node-local build で一度コンパイルし、その source+runtime を 1 個の hash 付き `software.tar.gz` に保存した。各正式 task は bundle を node-local `/tmp` に展開し、clean source と binary/cache/parameter/runner の hash を照合してから計算する。各 task での再コンパイルは不要。実験 repo の lock に記録した read-only source repository の object database を `attach_git` に用いるため、投入時も `../../LinearKnot` が参照可能であることが必要である。

## dry-run と試走

- 7 plan の全 49 入力・schedule・条件・memory/wall envelope を確認し、21 script の `bash -n` が成功。全正式 task は `not_started`、attempt は空。
- b25/4–8 kb、b200/8–20 kb、b400/20–30 kb の 3 代表 task を node-local に展開し、native と observer の binary、source、parameter、cache、凍結した runner の hash を検証。RNA CLI は実行していない。
- `data/smoke/pk_scaffold.txt` は **24 nt**。この入力に全 7 幅、最短 RNAcentral 4,793 nt に b25/50 の計 **9 条件が成功**。これは development 実行で、正式 49 条件に混ぜない。

| 入力 | b | h=d | wall 秒 | peak RSS GiB | log Z | F (kcal/mol) | MFE (kcal/mol) | BPP 生 bytes | gzip bytes |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| URS000233CDDE | 25 | 30 | 10.17 | 0.211 | 2123.238574275 | -1308.596044669 | -1232.700 | 1,870,861 | 711,931 |
| URS000233CDDE | 50 | 50 | 24.56 | 0.346 | 2174.842022941 | -1340.400322169 | -1259.600 | 3,139,211 | 1,150,227 |

- 実 RNA の b50−b25 は ΔF=−31.8042775 kcal/mol、自由 2,991 塩基の TV 平均 0.141866、p90 0.565138、最大 0.999642。小さい幅間では posterior の局所差が大きい位置もある。b100 以上の安定性はまだ評価していない。
- native CLI と比べたのは 24 nt の b25/50/400 と実 RNA の b25/50 の **5 条件**。stdout と BPP は全て byte 一致。
- 2 task の `--resume` で JSONL 全 9 行が byte 一致し、再計算されなかった。
- 実 RNA b400 を 0.05 秒で強制 timeout にした別の診断では、`timeout` と elapsed 約 0.055 秒を JSONL/tar に保存し、再開でも timeout を保持した。GNU wall/RSS は null。これは runtime の予測ではなく失敗保存経路の検査。
- 既知の小分布による TV、非対合状態、固定塩基の除外、座標和集合、重複座標、確率和、raw/gzip hash、p90 の **新規 7 tests** と既存 packed storage の **8 tests** が成功。quota 超過を模擬した既存 test の `Disk quota exceeded` は期待した検査出力。
- 試走の 2 図を開いて確認し、7 本の曲線を使う配置でも PNG/PDF の軸ラベル・凡例・title が canvas 内に収まることを検査した（配置検査の値は合成で、科学的な測定ではない）。
- 未実行の正式 49 条件を集計し、49 `not_started`、91 unavailable comparison を保持することも確認。欠測を 0 や b400 で代用しない。
- 最初の staging dry-run で source の `.github`/`.gitignore` も `.git` と一緒に除いてしまう packing filter の問題を検出し、path component が `.git` の場合だけ除くよう修正した。最終の 3 staging checks は成功。interactive allocation の `SGE_TASK_ID=undefined` は development 試走の child env から除き、架空の array ID は与えていない。

raw、native parity、timeout の JSONL/tar、`validation.json`、試走の表と 4 パネル/隣接差図は `results/local/beam-sweep-rnacentral-v0.2.1-20261027/smoke/`。図は集計・描画の動作確認用で、論文用の結果ではない。

## 投入前の検証と qsub コマンド（未実行）

実験 repo の root から実行する。`TSUBAME_GROUP` は著者が使う accounting group を設定する。下記は **各 array が完了してから次を投入**するので、この sweep 自体の同時 task は最大 3（最も短い群の task 数）。他の実行と合わせて平日の 30 tasks 上限内に収める。指定済み h_rt/resource を上書きしない。

```bash
cd /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments
export PYTHONDONTWRITEBYTECODE=1
python3 scripts/verify_beam_sweep_20261027.py
lfs quota -p 102571 /gs/bs
: "${TSUBAME_GROUP:?Set the authorized accounting group}"
for b in 25 50 75 100 200 300 400; do
  for tier in rna_4_8kb rna_8_20kb rna_20_30kb; do
    qsub -g "$TSUBAME_GROUP" -sync y -tc 3 \
      "jobs/beam-sweep-rnacentral-v0.2.1-20261027-b${b}/${tier}.sh"
  done
done
```

21 個の実際の script 名・task range・h_rt は機械可読計画の `arrays` と `jobs/beam-sweep-rnacentral-v0.2.1-20261027-b*/` にある。長い群だけを後から投げるなら、先に `tier` を 4–8/8–20 kb に限定し、20–30 kb の各幅を同じコマンドで別途実行する。複数 arrays を同時投入する場合は `-tc` の合計と他の arrays を含めた全体上限を著者側で割り当てる。

例：b400 の最長群だけを投入するコマンド（この task も未実行）。

```bash
qsub -g "$TSUBAME_GROUP" -tc 2 \
  jobs/beam-sweep-rnacentral-v0.2.1-20261027-b400/rna_20_30kb.sh
# controller が中断した task 2 だけ再開する場合：残存 process を確認した後
qsub -g "$TSUBAME_GROUP" -t 2-2 -tc 1 \
  jobs/beam-sweep-rnacentral-v0.2.1-20261027-b400/rna_20_30kb.sh --resume
```

## 実行後の集計と図

```bash
python3 scripts/analyze_beam_sweep_20261027.py \
  --plan-report docs/beam_sweep_rnacentral_20261027_plan.json \
  --scratch-dir /tmp \
  --out results/local/beam-sweep-rnacentral-v0.2.1-20261027/analysis-complete
```

`observations.tsv` は F、log Z、MFE、F_b−F_400、F_b−F_previous、各差/nt、時間、RSS、上限、状態を保存する。`comparisons.tsv` は b400 および隣接幅との TV の平均・p90・最大、最大 pair 確率差と座標、比較可否を保存する。隣接差は幅を昇順に並べた**新しい幅−直前の幅**。p90 は自由塩基の TV を昇順に並べ、index 0.9(N−1) の線形補間を使う。

各位置の非対合確率は `q_i=1−Σ_j p_ij`、TV は `(abs(q_i(b)−q_i(c))+Σ_j abs(p_ij(b)−p_ij(c)))/2`。平均・p90・最大は scaffold が `.` の塩基だけで計算する。入力/出力/hash、固定対の確率 1、確率範囲と incident mass≤1（許容差 1e−9）を検査し、丸めの許容差以内だけ [0,1] に戻す。分布の再正規化はしない。dense n×n 行列を作らず、node-local SQLite と対の sorted union の stream を用いる。

`beam_sweep.png/pdf` は F_b−F_400、b400 に対する平均 TV、wall time、RSS の 4 パネル、`beam_sweep_adjacent.png/pdf` は隣接幅の F/平均 TV と b400 に対する p90/最大 TV。全て横軸 b、線は RNA ごと、b100 に縦点線。b400 が欠測の RNA の reference 差は欠測のまま、隣接幅の成功比較は残す。図は valid な完了観測だけを結び、上限超えは TSV の status/timeout/elapsed に残す。

集計は array の task JSONL/tar を直接読むので、`aggregate_array.py` の実行は必須でない。必要なら既存 aggregate も plan ごとに再生成できる。script は source/model/binary が違う run を混ぜない。CLI が成功していても、task の終了時の依存検査が未確定なら `provenance_unverified` として比較に使わない。

## 再生成に使うスクリプト

- `scripts/build_beam_sweep_20261027.py --scratch /tmp/NEW_BUILD_PATH`：fresh な node-local build と bundle/lock。既存 bundle は上書きしない。
- `scripts/prepare_beam_sweep_20261027.py`：fresh な 7 plan/21 script と推定表。既存 config/run/script は上書きせず、qsub は呼ばない。
- `scripts/beam_sweep_array_20261027.py`：bundle stage、同じ standard runner、gzip BPP packing。
- `scripts/verify_beam_sweep_20261027.py`：49 未開始 task の整合性、script 構文、3 staging checks、live inode gate。
- `scripts/smoke_beam_sweep_20261027.py`：fresh な development 試走。`--collect-only` は保存済み native raw と JSONL/tar の集計を再検証する（既存 `smoke/analysis` は再生成先を空にして使う）。
- `scripts/analyze_beam_sweep_20261027.py`：実行後の集計/作図。出力先は新規 dir を指定する。

検証を再実行するコマンド：

```bash
python3 -m unittest discover -s tests -p 'test_beam_sweep_20261027.py' -v
python3 -m unittest discover -s tests -p 'test_task_storage.py' -v
```

## 2026-10-28：全 49 条件の集計と論文用図

保存済みの正式 49 tasks を集計した。全 49 tasks / 49 CLI が成功し、未完了・timeout・メモリ上限による終了は 0。全 49 出力と 91 比較（対 b400 の 49 組、隣接幅の 42 組）が有効だった。新しい RNA 計算や scheduler 投入はしていない。

入力 7 組の hash は元の E2 と全て一致し、全幅で同じ凍結 S/G、モデル、source `2101e5d86959661ba0a3da6d598f8f202e87186b`、計測用 FoldMetadata binary `8f505927b05bbaa8d160cdf76b2190b01ac6ca158d0104ab15bf565cc099a27e` を確認した。各 task の終了時の依存検査が有効で、stdout と丸め前 scalar の一致、圧縮 member と非圧縮 BPP の hash、固定対、確率範囲・確率和・座標重複を検査した。

[全 49 条件の値](../results/analysis/beam-sweep-rnacentral-20261028/observations.tsv)に各 RNA の F_b、log Z、F_b−F_400、隣接幅差、各差/nt、TV の平均・p90・最大、時間、最大メモリ、上限、状態を保存した。[全 91 比較](../results/analysis/beam-sweep-rnacentral-20261028/comparisons.tsv)に隣接幅の TV と最大の対確率差・座標も保存した。[幅別の範囲](../results/analysis/beam-sweep-rnacentral-20261028/by_beam.tsv)、[b100 の 7 本](../results/analysis/beam-sweep-rnacentral-20261028/b100.tsv)、[集計来歴](../results/analysis/beam-sweep-rnacentral-20261028/summary.json)も参照する。大きい生データは既存の ignored `results/runs/` に保ち、SQLite は `/tmp` で作って集計終了時に削除した。

### 自由エネルギーと塩基ごとの距離

以下の範囲は 7 RNA 間の最小–最大。TV の平均・p90・最大は、**各 RNA 内の scaffold で自由な塩基**に対して求めた値で、RNA 間の平均 TV の p90 ではない。非対合状態を含む。自由エネルギーの単位は kcal/mol、隣接差は大きい b の F から直前の幅の F を引いた値。p90 は index 0.9(N−1) の線形補間。

| b | F_b−F_400 | F_b−F_previous | 平均 TV 対 b400 | TV p90 対 b400 | 最大 TV 対 b400 |
|---:|---:|---:|---:|---:|---:|
| 25 | 39.57–574.81 | — | 0.1305–0.3119 | 0.4706–0.8809 | 0.994904–0.999949 |
| 50 | 7.76–332.50 | -270.97–-31.80 | 0.05667–0.2145 | 0.2014–0.7808 | 0.847529–0.999977 |
| 75 | 3.85–154.20 | -178.30–-3.91 | 0.0247–0.2059 | 0.09254–0.7872 | 0.79729–0.999949 |
| 100 | 2.22–107.93 | -52.32–-1.63 | 0.0128–0.1511 | 0.02992–0.664 | 0.53763–0.999926 |
| 200 | -0.01762–28.45 | -81.72–-2.24 | 0.0002129–0.06167 | 1.061e-08–0.244 | 0.0234073–0.999111 |
| 300 | 5.702e-08–11.43 | -18.81–0.02333 | 3.653e-10–0.04676 | 1.71e-10–0.1166 | 2.24562e-08–0.999111 |
| 400 | 0–0 | -11.43–-5.702e-08 | 0–0 | 0–0 | 0–0 |

| RNAcentral ID | F_100−F_400 | 平均 TV | TV p90 | 最大 TV |
|---|---:|---:|---:|---:|
| URS000233CDDE | 2.22 | 0.0128 | 0.02992 | 0.53763 |
| URS00019626CE | 32.42 | 0.1511 | 0.664 | 0.991672 |
| URS000239D8A4 | 5.28 | 0.04955 | 0.2346 | 0.855594 |
| URS00023C9544 | 22.79 | 0.06787 | 0.2951 | 0.982848 |
| URS000241043C | 71.64 | 0.1317 | 0.5245 | 0.999926 |
| URS000238C72C | 107.93 | 0.12 | 0.5255 | 0.999121 |
| URS00023A4DB9 | 97.68 | 0.07794 | 0.2814 | 0.999922 |

- b100 の最大 energy 差は URS000238C72C（22,543 nt）の **107.928882 kcal/mol**。BPP の平均・p90 が最大なのは URS00019626CE（5,482 nt）で、平均 **0.151128**、p90 **0.663988**。最大の塩基ごとの TV は URS000241043C の **0.999926**。
- b300 でも URS000238C72C の差は **11.429979 kcal/mol**、TV 平均 **0.0361146**、p90 **0.00913625**、最大 **0.999111**。最大の対確率差は座標 (6658,9368) の **0.990611**。平均や p90 だけから全位置の安定性を主張できない。
- URS000233CDDE は F_200−F_400 = **−0.0176198 kcal/mol**、F_300−F_200 = **+0.0233333 kcal/mol**。符号・単調性・厳密値への収束を仮定しない。b400 は exact ではなく、この差から保持確率質量や厳密誤差を推定しない。profile/accessibility の感度は今回測定していない。

### 時間と最大メモリ

各条件 1 回、予備実行 0、single thread、1 process/task。cpu_16 の共有ノード上の GNU wall time と peak resident memory。時間は MFE・分配関数・outside・BPP 出力を含み、build・scaffold 構築・gzip・集計を除く。既存 E2 の三測定中央値と混ぜない。全 CLI の wall time の和は **28,453.30 秒（7.904 process 時間）**で、CPU 時間・scheduler allocation 時間ではない。

| b | 時間の範囲（秒） | 最大メモリの範囲（GiB） | 7 本の時間の和（秒） |
|---:|---:|---:|---:|
| 25 | 10.0–74.3 | 0.207–1.204 | 265.79 |
| 50 | 24.3–157.3 | 0.346–1.942 | 593.41 |
| 75 | 40.7–277.2 | 0.490–2.690 | 1,026.24 |
| 100 | 55.7–423.5 | 0.588–3.222 | 1,580.51 |
| 200 | 131.5–1,211.7 | 0.980–5.681 | 4,368.26 |
| 300 | 213.2–2,328.1 | 1.361–8.691 | 8,234.49 |
| 400 | 295.0–3,592.6 | 1.611–10.615 | 12,384.60 |

全条件の時間は設定した CLI timeout 未満。メモリ上限（仮想アドレス空間の制限）による失敗はない。最大メモリは b400 の URS00023A4DB9（29,106 nt）の **10.615105 GiB**、最長の時間は b400 の URS000238C72C の **3,592.63 秒**。

### 採用図と再生成

[図 PDF](../results/figures/rnacentral_beam_sweep_20261028/rnacentral_beam_sweep_si.pdf)、[PNG](../results/figures/rnacentral_beam_sweep_20261028/rnacentral_beam_sweep_si.png)。165 × 122 mm、基本 7 pt・最小 6.5 pt、4 パネル（F_b−F_400、平均 total variation、時間、最大メモリ）。全 7 RNA・全 49 点を各パネルに示し、b100 を縦点線と大きい塗りつぶし記号で強調。PF、RSS、TV の略語を図内で使わない。PNG を開いて軸・凡例・点・文字の欠けと重なりを確認し、隣接する 75/100 の目盛りは回転させた。論文 repo の `figures/rnacentral_beam_sweep_si.pdf` と PNG に export し、manifest を更新した。

```bash
python3 -B scripts/analyze_beam_sweep_20261027.py \
  --plan-report docs/beam_sweep_rnacentral_20261027_plan.json \
  --scratch-dir /tmp --no-plot \
  --out results/local/beam-sweep-rnacentral-v0.2.1-20261027/NEW_ANALYSIS
python3 -B scripts/plot_beam_sweep_20261028.py \
  --source results/local/beam-sweep-rnacentral-v0.2.1-20261027/NEW_ANALYSIS \
  --summary-out results/analysis/NEW_SUMMARY \
  --figure-out results/figures/NEW_FIGURE
```

集計の既存 7 tests が成功。論文では SI S3.7・図 S5 と本文の近似精度の段落に反映し、b400 が近似であることと大きい局所差を残した。原稿・PDF の検証は論文 repo の `docs/beam_sweep_rnacentral_20261028_validation.json` に記録する。commit なし。
