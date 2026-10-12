# Experiment preparation progress

## 2026-10-11

### Runner and verification

- Added optional `order.kind: blocked_random` with an integer seed. Each
  repetition shuffles all input × mode × program combinations; warm-ups precede
  each combination's first measurement. The default nested order is unchanged.
  `run.json` keeps the planned schedule and actual attempts, including starts
  that do not finish.
- Added `--resume`. Config bytes/snapshot, software lock, source/binary state,
  build settings/cache, energy parameters and source/copied input hashes must
  match before the existing run is changed. Every `result.json` is retained,
  including failed jobs and results not yet appended to the manifest. Missing
  results restart with old logs preserved under `incomplete_attempts/`.
  Resume times, environment, node and previous checks are recorded. Legacy
  fixed-order manifests are supported.
- Fixed interrupted-run summaries: missing finalization fields give
  `run_status: interrupted` and null end-of-run change checks.
- `python3 -m unittest discover -s tests -v`: all 18 tests passed, including
  seven new runner/order/resume/summary checks and one scaffold subset check.
  `python3 -m pytest tests` was unavailable (`No module named pytest`); no
  package was installed. `git diff --check` passed.
- Summarized the original `e2-rnacentral-20261008` directory into ignored
  `results/local/e2-rnacentral-20261008-interrupted-prep-20261011.{tsv,json}`.
  It correctly reports four completed measurements (one success, two failures,
  one timeout) and an interrupted run. The original manifest stayed unchanged
  (SHA-256 `2d1617f01fc60fce7d684a5f8d10a671b06c992c7c477e1cbf963062ee56ae7e`).

### Inputs and the next feasibility measurement

- Extended `prepare_biological_inputs.py` with `--timeout-seconds` and `--ids`,
  preserving its default 120-second budget and original scaffold construction.
- Reused the saved G for URS00023C9544 (11,832 nt, original 64.98 s) and
  URS000241043C (13,652 nt, original 81.82 s), with unchanged input/G hashes.
- Regenerated URS000238C72C (22,543 nt) successfully in 137.592357 s, retaining
  4,465 pairs, and URS00023A4DB9 (29,106 nt) successfully in 182.176338 s,
  retaining 6,539 pairs. They ran sequentially with thread settings at 1 on
  node r8n10, scheduler job 8970385, using 3,600 seconds and the same 8 GiB
  address-space cap. LinearPartition-V b450fb3e, wrapper/binary hashes, beam
  200, dangles 2, posterior threshold 0.7 and minimum helix length 2 matched the
  original manifest and stayed unchanged during preparation.
- Saved the two new inputs and dataset under `data/rnacentral-retry-20261011/`,
  retaining the initial timeout records. Exact commands/results and raw-log
  hashes are also in `manifests/rnacentral-scaffold-retry-20261011.json`.
- Added `configs/e2_rnacentral_long_NEXT_20261011.json`: seven original records,
  three CLIs, beam 100 and local limits 30/30, matching output flags,
  3,600 seconds/job, 64 GiB address space, warm-ups 0, repetitions 1 and seed
  20261011. Removed the former 8,191 nt exclusion in this new matrix only.
  Old configs and measurements remain historical records.
- Updated README, protocol and the biological pilot notes, including the
  requirement to avoid concurrent calculations on the measurement node.

### Scientific scope and remaining work

The new matrix is prepared; its conditional CLI jobs have not been run.
Scaffold preparation success establishes saved G availability, not conditional
CLI feasibility, approximation accuracy or timing rankings. The author must
replace `manifests/software-pilot-NEXT.json` with a verified lock for the next
fixed commit/build and set `model_settings.build_revision_verified` to true
before running the matrix. The selected software must support all seven lengths.

All edits are inside the experiment repository. The external LinearCParty
working tree was not changed or built. No git commit or push was performed.

## 2026-10-11 TSUBAME array-job preparation

### Changed files and behavior

- Added `scripts/plan_array.py`, `scripts/run_array_task.py` and
  `scripts/aggregate_array.py`. Plans freeze input/config/lock copies, assign
  inclusive length tiers, pack complete sequence matrices by time/memory
  estimates, and generate one Grid Engine array script per populated tier.
  Scaffold variants of a shared sequence remain together. Planning never
  submits jobs and refuses existing output/script directories.
- Refactored `scripts/run_cli.py` to expose its existing provenance/resume path
  to array tasks and support bounded process concurrency for accuracy arrays.
  Seeded dispatch order, warm-up dependencies and repetition barriers are
  retained. Timing arrays require `node_f=1`, one worker and single-thread CLIs.
  Formal configs require verified build flags in both config and software lock;
  direct CLI execution remains development-only. SIGTERM is handled by the task
  wrapper so available interruption records are finalized.
- Updated `scripts/summarize.py` and `scripts/validate_outputs.py` for the merged
  run, including task counts and planned/completed/missing measurement counts.
  Aggregation reads result markers, recovers missed manifest updates, preserves
  failures and emits explicit rows for missing jobs. Task manifests retain
  scheduler IDs, host, resources, timestamps and attempts. Kernel-released task
  locks reject duplicate controllers. Shared manifests are updated only by
  explicit aggregation.
- Added `tests/test_array.py`, `configs/arrays/e2_rnacentral_20261011.json` and
  `configs/arrays/e1_bprna_20261011.json`. Updated `.gitignore` for generated job
  scripts. Updated README and `docs/experiment_protocol.md`; added
  `docs/array_dry_run_20261011.md` and exact plans/scripts/summaries under
  `docs/array_examples/20261011/{e1,e2}/`. Original experiment configs are unchanged.

### Verification and dry-run artifacts

- `python3 -m unittest discover -s tests`: 27 tests passed, including all 18
  existing tests. The new tests use small dummy CLIs and cover packing and
  estimates, scaffold grouping, resource/wall/memory guards, SGE environment
  selection, immutable preparation, matching and rejected resumes, task locks,
  worker bounds, warm-up/block barriers, SIGTERM, startup failure and partial
  aggregation/summary/output validation. No biological CLI matrix was run.
- Python compilation and `git diff --check` passed. Both saved plan checksums,
  all five saved script checksums and `bash -n` for those scripts passed.
- Prepared `results/runs/array-example-e2-20261011` from the existing seven-input
  `e2_rnacentral_15aa7742_20261011.json`, using the read-only dedicated
  `deps/LinearCParty-15aa7742` and `build/linear-cparty-15aa7742`. The software lock
  verified. There are three arrays and five tasks: three 4.8–7.0 kb inputs share
  a task; 11.8 kb and 13.7 kb inputs each have a task; 22.5 kb and 29.1 kb inputs
  each have a task. There are 21 planned CLI jobs, all not started.
- Prepared `results/runs/array-example-e1-20261011` from the existing 12-input
  `e1_bprna_pilot_20261008.json`, using read-only `deps/LinearCParty-20261008`
  and `build/linear-cparty-20261008`. The software lock verified. Two arrays
  contain one task each (8 inputs below 200 nt and 4 inputs at 200–800 nt), with
  two accuracy workers and 72 planned CLI jobs, all not started.
- Aggregated and summarized both unstarted matrices: 5/2 not-started tasks,
  21/72 missing measured jobs and zero completed measurements. Generated scripts
  live under `jobs/array-example-e{1,2}-20261011/`; review copies are in docs.

### Scientific scope and remaining checks

The example cost coefficients are explicitly uncalibrated scheduling assumptions.
Before formal collection, calibrate them on the adopted fixed-build pilot and
confirm longest-input CLI feasibility, memory and an appropriate failure budget.
Historical input/config purposes remain development; these dry-runs establish
preparation and aggregation behavior, not biological accuracy or runtime scaling.
If an interrupted timing task resumes on another node, examine its attempt hosts
and CLI resume history before using paired comparisons.

The resource table was checked against the official TSUBAME4 handbook, section
5.1.1, linked in README. Live `qconf -sc`, `-spl` and `-sql` checks failed with
`getgrgid(3215) failed: gid resolve error`. Current accounting group, scheduler
limits, concurrency and runtime modules are marked **要確認** in README.

All changes remain inside `experiments/`. No qsub, biological measurement,
external source modification/build, manuscript/textbook edit, git commit or push
was performed.

## 2026-10-11 bpRNA 全体と CParty 比較の array 準備

- 全clean集合26,063本のcohort、G生成array、E1全体、長配列b300/b400の段階的な予備計画、E0の300本subset＋4証人を保存した。G79 tasks / 26,063 jobs、E1 841 tasks / 699,387 jobs、E0 21 tasks / 912 jobs。予備幅は各12 tasks / 72 jobs。全29 scriptsはbash構文検査を通った。
- 12本E1の全324 jobsのwall/RSS、RNAcentralの成功とtimeout、元run hashをcalibration manifestに保存し、tierとmodeごとの費用式を設定した。長配列・上限なし・AD oracleの係数は暫定で、実測による見直しを計画に残した。通常優先度の概算は全計画345.81ポイント、全taskがh_rtまで使う枠は735.43ポイント。公式規則を確認し、group・同時実行数・queue・残額・storage/i-nodeを要確認とした。
- 固定clean cloneのbuild.shを使用し、このrepoのbuild下にupstream59593650＋7patchesのoracleを作った。Depsのローカルmirrorとsource/binary/parameter/cache/patch hashを保存した。調整中の実装repoは読むだけで、変更していない。
- 固定ライブラリにリンクする互換/default API observerと、oracle PFの外部コピーを微分するBPP observerを追加した。独立した2構造ensembleで微分を検査し、元oracle objectsのfull-precision probeともroot値を照合した。Profileは両Linear設定で保存するが、独立したanalytic oracleの一致とは扱わない。
- 開発上の検査：Python全36 tests成功。exactの399/400境界、per-mode費用、tier cap、分割schedule、G未生成の拒否、G completion markerとCLI journalからの回復、幅stageの停止/再開、失敗保持、上限条件別の最大幅参照と欠測を確認した。52 ntだけのG生成で既存pilotとinput/G hashが一致した。
- 科学的な小確認：4証人の12 CLI jobsは成功し、互換対oracleのlogZ/BPPは最大約2.14e-14/1.25e-14。45 ntではclosureによる差を観測した。BE41のMFEには互換でも0.0001706856 kcal/molの差があり、事前toleranceを変えずparity failureとして保存した。丸め・エネルギー規約が候補だが原因は実装担当に要確認。全bpRNAについて一致したとの主張はしていない。
- G全体は未実行。G未生成のE1/E0 dry-runは実行できず、全件Gを検証して新しいready planを作る。qsub、commit、pushは実行していない。論文repoの他のファイルは変更していない。

計画、再生成手順、各tierの時間とポイント、review用artifactsは
[計画の要約](full_bprna_array_plan_20261011.md) にまとめた。

## 2026-10-11 E0 を更新された単一 oracle と PF 比較へ変更

- 開発 tree を file:// で clean clone し、`deps/LinearCParty-2456f9b2/` を固定した。clone の変更していない build.sh で upstream `800e95a` + 0001/0002/0006/0007 の4パッチ oracle を `build/cparty-oracle-800e95a-4patches/` に build（終了値0）。manifest の base・順序・patch/source/binary/object hash を検証し、新しい oracle lock に保存した。旧7パッチ lock は変更せず残した。Closureパッチはなく、oracleは1種類。
- E0 の LinearCParty は `15aa7742` のまま。E1/E2の元buildとexact/RNA library hashを検証し、inside専用observerをリンクした。Oracleも元objectsへ直接リンクする。AD/outside/BPP/MFE/profileをE0実行と報告から外し、既存AD道具は歴史的診断用と明記した。
- 科学的な小確認：版間14入力・両policyの28/28比較で log Z / PF free energy が同値（最大差0）。新oracle対LegacyCompatibleは4証人すべて一致、12/12 jobs成功、最大絶対差はlog Z 2.1316282072803006e-14、free energy 1.4210854715202004e-14 kcal/mol。既定closureの差はnested45の1/4入力で検出した。BE41の過去のMFE差は今回の比較対象外であり、MFE/BPP一致は主張しない。
- Config・array config・生成元script・解析をPFだけへ変更。実装側の絶対丸め誤差許容（log Z 1e-9、free energy 1e-10 kcal/mol、両方必要）を採用した。Closure差はcompatible parityが通った入力で数え、失敗と帰属できない入力も元の分母へ残す。
- 新dry-runは実依存lock検証済み、304入力・21tasks・912jobsすべて未開始、5scriptsはbash構文検査成功。Gはpendingのため実行不可。旧envelopeを保守的に残した見込みは逐次24.11275 h、保存済み課金式で6.49496ポイント、h_rt envelopeで21.152ポイント。Inside専用の短tier実測による再校正、全件G、scheduler/group/module確認は未完了。
- 開発上の検査：`python3 -m unittest discover -s tests` は全42テスト成功。PF比較の範囲、両許容、zero/nonfinite/RT、policy誤り、closure帰属と分母、重複・改変入力・方向違反を追加検証した。元のAD診断テストも残るがE0一致の証拠には使わない。
- すべての変更は実験repo内。開発treeの変更・checkout・build、qsub、重い計算実験、commit、push、論文repoの他のファイルの変更は行っていない。

[更新版計画](full_bprna_array_plan_20261011.md) と [再現・証拠一覧](e0_oracle_update_20261011/README.md) を参照。

## 2026-10-11 共有 project のファイル数事故への対応

- CLI array と G array のジョブ別 directory/生出力を node-local `$TMPDIR`（なければ明示した `--scratch-dir`、`/gs/bs` 以下は拒否）へ移した。共有結果は task ごとの `task_results.jsonl` と `task_outputs.tar`。CLI の `run/` は5ファイル、G は3ファイル。完了ごとに tar を fsync してから JSONL を追記・flush/fsync し、再開は JSONL を正本として成功・失敗・timeout を保つ。途中行と未確定 tar 末尾から回復し、共有書き込み失敗時は残りの worker を停止する。
- 計画の入力凍結を全 RNA で1つの `inputs.tar` に変更。G 公開 dataset も `dataset.json` と `inputs.tar` の2ファイルにした。集計・summary・accuracy・E0・domain validation・cohort resolution は member をストリームで読み、共有領域へ展開しない。旧 directory と旧 tar/tar.gz の互換読み込み、旧 manifest 更新前の result marker 回復、archive 内の task 指定にも対応した。
- 両 planner に共有 entry 数（directory、scheduler log、一時ファイル余裕込み）の見積もりと作成前拒否を追加。project 102571 の live quota と、config の保守的上限2,400万の小さい値を使い、残りの10%以内に制限する。quota 取得失敗時は明示 fallback が必要。全件 config の fallback は80万。旧 G lock を保持し、改訂 code の8 hash を持つ新しい packed G lock を作成・検証した。LP source/wrapper/binary は変更していない。
- E1 全体の `--dry-run --estimate-only` は実依存 lock と既存 G/input hash の検証を通過。26,063 RNA・841 tasks・699,387 jobs、共有領域は **10,963 entries = 841 × 13 + 30**。live 使用数23,118,679、有効上限24,000,000、残り881,321、10%枠88,132に対し約1.24%で通過した。live の表示上限28,000,000と事故時の24,000,000の差を明記した。dry-run は1つの JSON と directory のみを作り、task directory・全 schedule・入力コピー・CLI 出力は作っていない。
- 開発上の検証：`python3 -m unittest discover -s tests` は **52 tests、17.295秒、全成功**。小さな実 CLI fixture の新旧 accuracy/E0 値、summary、失敗と欠測の分母、段階停止/再開、完了 JSONL と stale manifest、途中 JSONL と orphan tar、quota 拒否前の無書き込み、G の新旧形式・公開 tar、SIGTERM と共有書き込み失敗時の worker 停止を確認した。fixture は一時領域で作成し削除した。既存の synthetic pilot テストも共有 data directory を使わないよう修正した。
- 既存の実 archive も読み取りだけで確認：E1 `nt50_99/000001` の13,500計画ジョブから一件の stdout（1,939 bytes）と入力（119 bytes）を読み、G tar の79 tasks / 26,063 jobsを集計して全26,063件の完了を確認した。展開ファイルは0。G 全件を今回再生成した、または E1 全件の科学的な結果を検証したとの主張はしない。
- README、protocol、全件計画を更新した。すべての変更は実験repo内。qsub、git commit/push、外部の実装 tree の変更、論文repoの他のファイルの変更、大量のファイルを作るテストは行っていない。旧計画からの実行再開は新計画への置き換えが必要。投入直前の quota/fallback、ローカル scratch 容量、共有 byte quota、tar 転送を含む wall-time、scheduler/group/module、G 公開先と cohort の対応は要確認。

[E1 見積もり](packed_array_20261011/array_estimate.json)、
[検証記録](packed_array_20261011/verification.json)、
[改訂計画](full_bprna_array_plan_20261011.md) を参照。

## 2026-10-12 E1 short cohort に capb と採用条件を追加

- `analyze_bprna_short.py --extend <validated-output> --capb-run <packed-run>` を追加した。新helper `bprna_capb.py` は同じS/G cohort・model・binary・計算引数を検証し、以前の全task/checkpointとexact再実行のhashを照合してから旧scalar指標を再利用する。tarはmember streamで読み、単一SQLiteを `/tmp` で一時使用して削除する。共有領域へraw出力を展開しない。
- 固定版 `15aa7742` の24,625 RNAを、default30/30・capb=max(b,30)・nocap100000で比較した。capb b20と採用条件b20は既存defaultの同じ観測を共用し、b50/100/200は新capbを使う。TSVは実測11条件・270,875行、summaryは3上限条件と採用条件の4幅を表示する。既存197,000行の指標は全列そのまま保持した。
- b200 capbとnocapのΔFは全配列で差0。旧defaultより小さいΔFは3,743 RNA、大きいΔFは0（許容10⁻⁸）。採用条件では全体99.27%がΔF≤0.01 kcal/molだが、300–399 ntは14.98%が0.01、2.56%が1を超え、最大4.198 kcal/mol・最小保持率0.1101%。b50 capbとnocapの差には538 RNAでcapb大、1 RNAでnocap大があり、設定の一般的な単調改善は仮定しない。
- 採用条件の非増加はΔF 99.45%、BPP TV 97.48%、context TV 98.18%。旧defaultやnocapとの差は同じRNAのexact誤差指標の差であり、2設定のposterior同士のTVではない。LinearCPartyおよび他のCLIの時間・RSS比も保存し、並列accuracy runの参考値と明記した。0の分母と欠測を除き、分母を保持する。
- 59/59 unittestと24/24成果物検証が成功。exact失敗349件の指定再実行置換と元履歴を保持し、差し替え後886,500成功CLI観測、295,500組の整合を確認した。新capb 73,875組とexact再確認24,625組も成功。既存pilotの96条件・4,416 scalarは最大差0。alias・配列内差・資源比・採用条件の単調性をTSVから照合した。
- 成果物は `results/local/e1-bprna-short-capb-20261012/` のTSV 1、JSON 3、PNG 2（6ファイル、約247 MiB）。解析はsingle thread、全体1313.85秒・peak RSS 795,172 KiB。読み取り・集計だけで、qsub・git commit・外部実装変更は行っていない。論文側ノート `textbook/notes/e1_bprna_short_capb_20261012.md` に全表・主張の範囲・400 nt以上への含意・再生成コマンドをまとめた。

## 2026-10-12 Offline Spark MFE comparison

- Built unmodified Spark `3105c311ca6243a31788e74ac33640435e53b62c` out of tree in `build/spark-3105c311/`, CMake Release, g++ 11.5.0, C++17, one build process. Saved source commit/tree/byte hashes, build flags/cache/link command, binary and Turner04/DP09/PK hashes in `manifests/software-spark-3105c311-20261012.json`. No network, qsub, commit, or external LinearCParty changes.
- `scripts/compare_spark_mfe.py` runs the prepared 32 inputs with common Turner04 bytes, 37°C, dangles2 and GU allowed; it rechecks pinned LinearCParty `15aa7742` MFE under both closure policies. `spark_mfe_linear_probe.cc` observes output-constrained local energies, `spark_mfe_observer.py` generates separate Spark observers in `/tmp`, and `spark_mfe_domain_probe.cc` checks unconstrained short-span membership. Official Spark remains the reference; masked runs count as reverse evaluations only when the exact target pair set is returned. Trace-repair binaries are labeled diagnostics.
- Both closure policies yield the same comparison counts: 14 same structure/energy, 13 same structure/different energy, 4 different structure/same energy, 1 different structure/energy. Twelve energy differences come from local PK rounding. Two GU-stack differences of −0.047 kcal/mol come from Spark allowing its BE internal-loop alternative on an adjacent positive stack; intended model semantics need confirmation. Three structural differences are ties. The `fresh002` extension is excluded by LinearCParty's short partial-WMBP domain even though its VP is reachable.
- Source audit also found a nested-PK weight difference: LinearCParty's WI→WMB uses Psp+Pps=15.96 kcal/mol, while Spark uses Psm+Pps=11.03; Spark declares Psp but does not reference it in `Spark.cc`. This branch occurs zero times in the scored traces and does not explain the 32-input observations. Its 4.93 kcal/mol difference remains a model question, not a measured nested-PK discrepancy.
- `pk009` drops scaffold pair (12,17) in Spark traceback; observing the recomputed right WIP energy restores the pair without changing the MFE. A separate 6-nt boundary witness `CAAAGU` / `(...).` drops the added pair because Spark rejects valid `j=n` in `trace_VP`; changing that boundary in a diagnostic translation unit restores `([..)]` at the same 3.66 kcal/mol. These are traceback bug witnesses, distinct from the energy-model and state-domain differences. Unreliable terminal trace residuals involving INF are retained but not interpreted as physical loop energies.
- `scripts/validate_spark_mfe.py` passed packed validation of 273 command logs, 1,315 archive file members, input and artifact hashes, unchanged sources/libraries, all earlier 64 MFE values, exact observer/official stdout identity, and local sums. Maximum local-sum discrepancy is 4.44×10⁻¹⁵ kcal/mol; reconstructing Spark local rounding and positive-BE alternatives explains all 32 MFE values to 3.55×10⁻¹⁵ kcal/mol. This is a finite development diagnosis, not a complete enumeration of Spark's domain or a timing benchmark.
- Saved `results/local/spark-mfe-comparison-20261012/{comparisons.jsonl,raw_outputs.tar.gz,config.json,summary.json,validation.json}` (5 files, about 3.1 MiB). Local run files are packed; 61 build files were hash-checked in the archive before removing regenerable intermediates. The shared file budget was below 200 entries, with initial project usage 22,166,785/28,000,000. Keep manuscript adoption pending until model/domain questions are resolved; the paper note contains the full condition table and per-case diagnosis.
