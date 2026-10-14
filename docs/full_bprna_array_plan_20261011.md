# bpRNA 全体の G・E1 と CParty 比較 E0 の array 準備

2026-10-11。qsub、commit、push は実行していない。
実装の作業 tree と論文 repo のファイルは変更していない。
計画・job script・係数合わせ・テストと、52 nt の G 生成、4証人の小確認、およびE0更新時の14入力・両policyの版間確認だけを行った。
本実験の結果ではない。

## 共有 project のファイル数事故後の更新（2026-10-11）

旧 E1 の投入後、研究室共有 project 102571 のファイル数上限に達し、
研究室全員の書き込みが `Disk quota exceeded` になった。
著者が array を停止し、旧 run を `results/archives/` の tar にまとめて約85万ファイルを空けた。
本更新では qsub、commit、push、実装 tree や論文 repo の他のファイルの変更は行っていない。

新形式 `task-jsonl-tar-v1` は各 CLI を node-local `$TMPDIR` で実行する。
`TMPDIR` がなければ `--scratch-dir` を明示し、`/gs/bs` 以下は拒否する。
結果・command・stdout/stderr のサイズ/hash/末尾要約・time・peak RSS を
タスクの `task_results.jsonl` に保存し、全生出力を `task_outputs.tar` にまとめる。
tar を flush/fsync してから JSONL に完了を追記し、flush/fsync する。
再開は成功・failure・timeout を保ち、JSONL の途中行と tar の未確定末尾だけを除いて、
未完了の CLI を保存済み順序で実行する。ローカルの途中出力は allocation 終了時に失われ得る。
共有領域の `run/` は5ファイルで、ジョブ別 directory を作らない。

入力凍結も全件で1つの `inputs.tar` にした。
G は1 task 1 tar と JSONL、metadata の計3ファイル。
G 集約後の公開 dataset も `dataset.json` と `inputs.tar` の2ファイルだけ。
旧 directory 形式と `results/archives/*.tar[.gz]` は member のストリーム読み込みで解析し、
共有領域へ展開しない。大きな旧 E1 archive は `archive.tar::tasks/<tier>/<id>/run` で task を指定できる。
旧計画は読み取り用に保持し、実行には新しい計画を作る。

両 planner は出力作成前に file-count budget を検査する。
`lfs quota -p 102571 /gs/bs` の残りを取得し、取得失敗時は config の
`file_budget.fallback_remaining_files` を使う（明示値がなければ拒否）。
残りの10%を超える計画は拒否し、見積もりと quota 原文・判定を plan に保存する。
事故時の上限2,400万を config の保守的上限に残した。
今回の live 表示は2,800万で異なるため、有効上限は両者の小さい値を使う。
fallback は80万とし、投入直前に実際の残りを再確認する。

[新しい E1 全体の dry-run 見積もり](packed_array_20261011/array_estimate.json) は、
実依存 lock と既存 G/input の hash 検証を通り、26,063 RNA・841 tasks・699,387 CLI jobs。
`--dry-run --estimate-only` で1つの JSON とその directory だけを作り、task directory、
699,387件の schedule、生出力、入力コピーは作っていない。実行用の計画ではない。

| 対象 | 共有領域の見積もり（directory と一時ファイル余裕込み） |
|---|---:|
| E1 全体 | **10,963 entries** = 841 × 13 + 30 |
| G 全体（79 tasks） | 743 entries = 79 × 9 + 32 |
| 新 E1 見積もり artifact 自体 | 2 entries（JSON + directory） |

E1 判定時の使用数は23,118,679、有効上限24,000,000、残り881,321。
10%の許容枠88,132に対して10,963なので通過（残りの約1.24%）。
初回 scheduler log、plan/config/schedule、task metadata/lock、tar/JSONL、directory、
同時 atomic-write temporary を含む。再投入の log と別の解析結果は別途数える。
検査は quota を予約しない。ローカル scratch 容量、共有 byte quota、tar 転送時間を含む
task wall-time、同時 task 数、group/queue/module は投入前に要確認。

G の旧 lock は保持し、改訂 script の hash を持つ別の packed lock を新 G config に指定する。
既存 G dataset を読むことと、新形式で G を再実行することは別である。
G の新公開先と cohort の `scaffold_dataset` を対応させてから ready plan を作る。

開発上の検証は一時領域の数ジョブで行い、後で削除する。
新旧形式の BPP/context/accessibility/E0 値、失敗の分母、停止・再開、JSONL の途中行、
tar の未確定末尾、集計、旧 archive と G archive、予算超過時の作成前拒否を確認する。
検証記録は [progress](progress.md) に追記する。性能や科学的な全件一致の主張は追加していない。

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

## E0 の oracle・入力・比較（2026-10-11 更新）

E0 の LinearCParty は E1/E2 と同じ固定版 `15aa7742` を維持する。
Oracle の定義と build 手順だけを clean clone `deps/LinearCParty-2456f9b2/` から採用した。
開発 tree は `git clone file:///gs/bs/tga-satolab-gtex/otgk/LinearCParty` の参照元として使い、
clone 側を `2456f9b2` に checkout した。開発 tree の変更・checkout・build は行っていない。

その clone の `test/oracle/legacy_cparty/build.sh <out_dir>` を変更せずに実行し、
`build/cparty-oracle-800e95a-4patches/` への Release build は終了値0で成功した。
Upstream base は `800e95af9d1a97555233e2c7b37d73dbaa1a92b4`、パッチは
**0001 → 0002 → 0006 → 0007** の4本。Closure のパッチはなく、公式 oracle は1種類だけ。
ネットワークでの再取得を避けるため、開発 tree 内に既にある oracle checkout の Git objects を
`file://` で新しい `out/src` へ clone した。以後の reset・patch 適用・build は新しい clone 内で行った。
`manifest.json` の base commit、patch 順序・hash、patched diff・source・binary/object hash を
固定 clone の `oracle_config.verify_checkout` で検証した。
[新しい oracle lock](../manifests/software-cparty-oracle-800e95a-4patches-20261011.json) は
従来と同じ lock 項目に加え、検証した manifest 全体とその hash を保存する。
[旧7パッチ lock](../manifests/software-cparty-oracle-59593650-7patches-20261011.json) は履歴として変更せず残した。

[外部 API observer](../tools/e0_linear.cc) は `15aa7742` の既存 exact/RNA libraries にリンクする。
`E0LegacyCompatible` は `LegacyCompatible`、`E0Default` は既定の `RequireWeakClosure` を選び、
inside の log Z・PF free energy と RT のみを full precision で出力する。
[Oracle observer](../tools/e0_legacy_scalar.cc) は元の4パッチ oracle objects へ直接リンクし、
MFE や sampling を求めず、スケーリング seed 0 の exact inside だけを観測する。
AD の source copy、outside、BPP、MFE、profile は E0 の実行・報告経路から外した。
既存の [AD 道具](../tools/e0_legacy.cc) と [AD header](../tools/legacy_ad.hh) は歴史的な診断用であり、
新しい E0 の報告には使わない。
[Build script](../scripts/build_e0_tools.py) と [E0 lock](../manifests/software-e0-15aa7742-20261011.json) は
元の library hash を照合し、4パッチ manifest・元 objects・道具・binary の hash を保持する。
新しい observer build は `build/e0-observers-15aa7742-800e95a-4patches/` に保存した。

比較量は **log Z と PF free energy のみ**。
[Config](../configs/e0_bprna_legacy_20261011.json) の合格条件は
`|Δlog Z| ≤ 1e-9` **かつ** `|ΔPF free energy| ≤ 1e-10 kcal/mol`。
`2456f9b2/test/tools/audit_legacy_inside_parity.py` の既定値を採用した、full-precision 値の絶対丸め誤差許容であり、
相対誤差許容や public CLI の表示桁による許容を追加しない。
実装担当の550件一致は実装側の証拠として参照する。今回このrepoで550件を再実行したとは扱わない。
既定モード対 oracle は参考値として記録する。
[Analysis](../scripts/analyze_e0.py) は default−oracle、default−compatible、compatible−oracle を分け、
compatible の両許容条件が通った入力だけで、default−compatible がどちらかの許容を超えた件数を数える。
失敗・timeout・artifact error・帰属できない入力を分母から落とさない。
Closure による制限方向は Δlog Z ≤ 許容、Δfree energy ≥ −許容でも点検する。

`15aa7742..2456f9b2` の src 変更（parameter-file validation、suboptimal、factorized forest）を確認し、
固定 Turner04 parameter file、37°C、dangles2、両 closure policy で小入力を比較した。
4証人＋長さ20/30/50/80の乱数8入力（x/Xを含む）＋VP split25＋BE境界41の計14入力、
**28/28比較が同値、最大 |Δlog Z| = 0、最大 |ΔPF free energy| = 0**。
最新側は clone の exact library target のみを別 build で作った。
[再現 script](../scripts/check_e0_version_parity.py) と
[値・入力・source diff・build/hash 記録](e0_oracle_update_20261011/e0_version_parity.json) を保存した。
有効な固定 parameter file に対する小入力確認であり、beam/suboptimal の一致や一般的な同値性は主張しない。

bpRNA は50–99/100–199/200–299/300–399 ntから各75本、計300本をseed20261011で抽出した既存 cohort を維持する。
[選定一覧](../manifests/cohorts/bprna-e0-subset-20261011.json) と [証人4本](../data/e0-witnesses-20261011/dataset.json) は変更していない。
同じ保存済み G を用い、BE65 の forced-unpaired X も保存する。

| E0 tier nt | 入力数 | tasks | CLI jobs | 時間 h | resource / h_rt |
|---|---:|---:|---:|---:|---|
| 1–49（証人3本） | 3 | 1 | 9 | 0.005 | cpu_4 / 8 h |
| 50–99（bpRNA75 + BE65） | 76 | 2 | 228 | 0.12 | cpu_4 / 8 h |
| 100–199 | 75 | 2 | 225 | 0.83 | cpu_8 / 8 h |
| 200–299 | 75 | 8 | 225 | 3.16 | cpu_40 / 8 h |
| 300–399 | 75 | 8 | 225 | 19.99 | node_q / 8 h |
| 合計 | **304** | **21** | **912** | **24.11** | 5 arrays |

[Array config](../configs/arrays/e0_bprna_legacy_20261011.json) は inside 専用 observer の計画として更新した。
費用式は旧E0の保守的な envelope（E1 exact時間×10・memory係数×6）を当面維持する。
AD を計算する費用とは説明しない。今回の小確認だけで、300 nt級のscalar inside費用を校正したとは扱わない。
短い tier の実測で時間・RSS・capを確認し、必要なら新しい計画版で下げる。
見積りは逐次仕事量24.11275 h、通常優先度の保存済み課金式で6.49496ポイント、
全taskがh_rtまで使う枠で21.152ポイント。Queue待ちや再実行は含まない。

[4証人 config](../configs/e0_witnesses_800e95a_4patches_20261011.json) の小確認は **12/12 jobs成功、4/4 PF parity成功**。
最大 |Δlog Z| は `2.1316282072803006e-14`、最大 |ΔPF free energy| は `1.4210854715202004e-14 kcal/mol`。

| 証人 | compatible−oracle Δlog Z | compatible−oracle ΔPF free energy (kcal/mol) | closure差を検出 |
|---|---:|---:|---|
| closure14 | −6.07e-16 | 3.71e-16 | なし（許容内） |
| nested45 | −2.13e-14 | 1.42e-14 | あり |
| BE41 | −1.78e-15 | 0 | なし（許容内） |
| BE65 | −7.11e-15 | 1.78e-15 | なし（許容内） |

Closure差は **1/4入力**（nested45）で検出した。
Default−compatible は Δlog Z = `−2.3319792319398402e-7`、ΔPF free energy = `1.437247227897842e-7 kcal/mol`。
Closure差の方向も全4入力で整合する。
以前のBE41のMFE差 `0.0001706856 kcal/mol` は履歴として残すが、今回のE0では **比較対象外**。
MFEやBPPが一致したとの主張へ広げない。
[解析結果](e0_oracle_update_20261011/e0_witness_analysis.json) と
[観測値・生ログhash](e0_oracle_update_20261011/e0_witness_observations.json) を保存した。

新しい dry-run は実依存 lock の検証を通り、304入力・21tasks・912jobsすべて未開始。
`pending_scaffolds=true` なので実行は拒否される。
[更新版の summary・task割当・5scripts](e0_oracle_update_20261011/E0/summary.json) を正本とし、
旧 `full_bprna_arrays_20261011/E0/` の artifact は履歴として残す。
G全件を生成・検証した後に新しい ready plan を作る。

### E0 の取り直し準備：d21618c5 と5パッチ oracle（2026-10-14）

新しい E0 は clean clone `deps/LinearCParty-d21618c5/`、固定 commit
`d21618c5115453174204329fcc40534b1cb9aa9c` を使う。
その clone の変更していない `test/oracle/legacy_cparty/build.sh <out_dir>` で
`build/cparty-oracle-800e95a-5patches/` を Release build し、終了値0で成功した。
既存4パッチ oracle の Git objects を `file://` で新しい `out/src` に clone してから実行した。
元の checkout と開発 tree は変更していない。
パッチ順序は config.json どおり **0001 → 0002 → 0006 → 0007 → 0008**。
0008 は WIP が3塩基の最小 hairpin の V 枝を数える修正であり、closure のパッチは追加していない。
[Oracle lock](../manifests/software-cparty-oracle-800e95a-5patches-20261014.json) に
manifest、patch/source/binary/object hash、build 設定・ログ hash を保存した。

[Build script](../scripts/build_e0_tools.py) は
[既存の d21618c5 pilot lock](../manifests/software-pilot-d21618c5.json) と
`build/linear-cparty-d21618c5/` の対応を検証し、exact/RNA library をリンク前後に hash 検査する。
新 oracle の元 objects に直接リンクした inside 専用道具は
`build/e0-observers-d21618c5-800e95a-5patches/` に build した（終了値0）。
API の変更対応は不要だった。
[E0 lock](../manifests/software-e0-d21618c5-20261014.json) に実依存を記録した。
比較は従来と同じ **exact inside の log Z / PF free energy のみ**で、
Turner04、37°C、dangles2、絶対許容 `1e-9` / `1e-10 kcal/mol` の両方を必要とする。
`LegacyCompatible` 対 oracle が本比較、既定 closure 対 oracle は参考値。
旧 lock と旧 config は履歴としてそのまま残した。

新しい [WIP29 入力](../data/e0-witnesses-20261014/WIP29.txt) は
固定 clone の `scfg2_wip_hairpin_boundary_test.cc::check_root_consumers` の29 nt証人を保存したもの。
G の (3,7) の間に3塩基 (4,5,6) があり、WIP 内の最小 hairpin 枝を含む。
実装側で用いる extension (13,29) の completed witness も dataset metadata に記録した。
ここでは構造のエネルギーを別に採点せず、同じ Turner04 条件で ensemble の PF だけを観測する。
元の4証人と合わせて **15/15 jobs成功、5/5 PF parity成功**。
最大 `|Δlog Z| = 2.1316282072803006e-14`、
最大 `|ΔPF free energy| = 1.4210854715202004e-14 kcal/mol`。
Closure 差は nested45 の1/5入力で検出し、制限方向は全5入力で整合した。

| WIP29 の観測 | log Z | PF free energy (kcal/mol) |
|---|---:|---:|
| d21618c5 LegacyCompatible | −21.226938775753077 | 13.082603367763157 |
| 5パッチ oracle | −21.226938775753080 | 13.082603367763157 |
| 15aa7742 LegacyCompatible | −21.255165363154305 | 13.1 |
| 4パッチ oracle | −21.255165363154305 | 13.1 |

旧2ジョブも成功し、旧版同士の PF は一致した。
新−旧の差は log Z 約 **+0.02822658740123**、PF free energy 約 **−0.01739663223684 kcal/mol**。
この例で両実装が最小 hairpin を含める変更を反映したことを確認した。
旧 binary・library・oracle source/object・manifest・parameter の hash を旧 lock で確認した。
編集済みの実験 script の hash を合わせるために旧 lock を書き換えてはいない。
[再現 script](../scripts/check_e0_fixed_witnesses.py) はジョブ別ファイルを `/tmp` に置き、
生出力を1つの tar にまとめる。保存 tar の読み込みでも解析値が一致した。
[観測値](e0_oracle_update_20261014/e0_witness_observations.json) と
[解析](e0_oracle_update_20261014/e0_witness_analysis.json) を保存した。

[新 config](../configs/e0_bprna_legacy_d21618c5_20261014.json) は
従来と同じ300本の bpRNA cohort と4証人、追加の WIP29 の **305入力**。
[新 array config](../configs/arrays/e0_bprna_legacy_d21618c5_20261014.json) は
199 nt以下を cpu_16、200–399 ntを cpu_40、全 tier を workers=1、h_rt=8 h とする。
300–399 nt の旧128 GiB cap は cpu_40 に収まらないため、cap を80,000 MiB、
task memory budget を81,920 MiBにした。既存の費用式と時間条件は維持した。

| E0 tier nt | 入力数 | tasks | CLI jobs | 見積時間 h | resource |
|---|---:|---:|---:|---:|---|
| 1–49 | 4 | 1 | 12 | 0.006 | cpu_16 |
| 50–99 | 76 | 2 | 228 | 0.124 | cpu_16 |
| 100–199 | 75 | 2 | 225 | 0.828 | cpu_16 |
| 200–299 | 75 | 8 | 225 | 3.163 | cpu_40 |
| 300–399 | 75 | 8 | 225 | 19.994 | cpu_40 |
| 合計 | **305** | **21** | **915** | **24.114** | 5 arrays |

新しい dry-run は実依存 lock と既存の full G dataset を検証済みで、
`pending_scaffolds=false`。915 jobs はすべて未開始。
凍結 plan は `results/runs/e0-bprna-legacy-dry-d21618c5-5patches-20261014/array_plan.json`、
job script は `jobs/e0-bprna-legacy-dry-d21618c5-5patches-20261014/` の
`nt1_49.sh`、`nt50_99.sh`、`nt100_199.sh`、`nt200_299.sh`、`nt300_399.sh`。
[Review 用の summary・割当・script copy](e0_oracle_update_20261014/E0/summary.json) を保存し、
全5 script が `bash -n` を通った。
共有結果は task ごとの JSONL/tar、compact plan と inode gate を有効にした。
新計画の見積り297 entries は、live quota と保守的24,000,000上限での
残り2,951,789 entriesの10%以内に収まる。

既存の費用式は scalar inside の300 nt級で未校正であり、短い tier の実測による
時間・RSS・cap・timeoutの点検、scheduler resource/h_rt制限、利用group・同時task数・module、
投入直前の inode/byte quota と node-local scratch容量は要確認。
94テストが成功し、旧 oracle の拒否、WIP の新旧差の両側検査、305入力の維持、
cpu_40 に収まらない旧capの拒否を追加した。
小確認は本実験の全件一致や性能の証拠ではない。
qsub、git commit、開発 tree の変更、重い計算実験は行っていない。
[再現手順と検証記録](e0_oracle_update_20261014/README.md) を参照する。

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
旧E1のraw配置は約373万ファイル＋約70万jobディレクトリとなり、共有projectでは投入禁止。
新形式のファイル数の事前検査は上の更新節を参照する。byte容量の確認も必要。
短いCLIを一つずつarrayへ投げず、数百本のmatrixをbatchすることでscheduler task数と最低課金の影響を抑える。

## Dry-run artifacts と実行前の再計画

[当初の機械可読 summary](full_bprna_arrays_20261011/summary.json) に全tier・task/job数・時間・ポイント・plan/script hashを保存した。
E0の依存とscriptは [更新版summary](e0_oracle_update_20261011/summary.json) を使う（見積り値は同じ）。
各サブディレクトリにtask割当TSVと29本のscriptのcopyがある。
実際の凍結plan、入力copy、task config、seed付き全scheduleは以下のignored領域。

```text
results/runs/g-bprna-full-plan-20261011-v3/
results/runs/e1-bprna-full-dry-20261011-v2/
results/runs/e0-bprna-legacy-dry-20261011-800e95a-4patches/
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
  --software-root deps/LinearCParty-15aa7742 --build-dir build/e0-observers-15aa7742-800e95a-4patches \
  --out results/runs/e0-bprna-legacy-ready-20261011
```

同じ方法でb300/b400のconfigとarray configを使い、新しいready planを作る。
G全件のstatus・ID・S・G hash・LP版・設定が一致しない場合は再計画を拒否する。
長配列を段階的に進める際はgenerated scriptへ `--through-beam` と、再開時の `--resume` を渡す。
E0の短いtierでinside専用の資源見積りを校正し、上限なしの長配列の資源見積りも実測で更新してから結果を解釈する。

## 検証

`python3 -m unittest discover -s tests`：更新後42テスト成功。
399/400のexact境界、modeごとの係数、tierのcap、分割schedule、G未生成の実行拒否、
Gのcompletion marker回復、変更されたGの拒否、幅ごとの停止・再開・失敗保持・journal回復、
上限条件別の参照と最大幅失敗時の欠測、歴史的AD診断の独立した2構造ensemble、forced-unpaired入力を確認した。
PFだけの比較、両絶対許容の適用、zero/nonfinite/RT metadata、closure帰属と分母、
policy誤り、変更input、重複観測、制限方向違反を追加確認した。
[検証記録](e0_oracle_update_20261011/verification.json) を保存した。
52 ntのGは既存pilotとinput/G hashが一致した
（[小確認](full_bprna_arrays_20261011/g_one_record_smoke.json)）。
