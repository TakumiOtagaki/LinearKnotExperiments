# Experiment preparation progress

## 2026-10-20 CParty / LinearKnot PF scaling preparation (not submitted)

- Prepared [the runbook](cparty_scaling_89baea52_20261020/README.md), configs, a verified observer lock and five new `cpu_16` array scripts under `jobs/cparty-linearknot-pf-89baea52-20261020/`. No qsub or git commit. Existing runs, E2 input bytes and external implementation working trees are unchanged.
- The original six-patch CParty CLI includes MFE, PF, sampling-derived BPP and decoding. The external `PFScaling` driver execs the existing E0 CParty observer (sampling=0) or calls adopted 89baea52 exact/beam inside APIs. All methods use compatible closure; beam h=d=100 and outer ranking. BPP is excluded. The main-text default weak-closure model and whole-CLI timings remain separate.
- Completed a local CParty pilot: three shapes at 250/500/1,000 nt, warm-up 1 plus measurements 3, 4 GiB AS cap. Successful input medians span 0.067–0.117 s, 0.367–0.467 s, and 2.722–3.423 s. The 1,000-nt long-stem warm-up and all measurements return a nonfinite PF; raw failures are retained. A separate seed probe changes finite PF values, so no scaling workaround was adopted. These numerical limitations may prevent measuring a crossover with this oracle.
- Preserved a broader interrupted pilot after a 120 s exact warm-up timeout and added one-shot calibration. Exact at 500 nt takes 26–45 s; beam100 at 1,000 nt takes 16.85–26.87 s. CParty/exact PF agree on all six 250/500-nt inputs (maximum absolute log Z difference 2.05e-12). These development observations are not manuscript timing results.
- Planning uses CParty cubic time/quadratic RSS, exact max(cubic, observed exponent), beam linear envelopes, a 1.5 factor, 6 h / 80 GiB process limits, and approximately doubled h_rt for complete 1+3 matrices. Exact is admitted through 1,000 nt, CParty/beam through 4,000 nt. Exact 2,000 nt cannot fit the doubled matrix allowance in 24 h; exact 4,000 nt also exceeds the process-time limit.
- The new plan has 105 inputs, five tasks, 1,092 CLI calls (819 measured, 273 warm-ups), 8/8/8/8/4 workers, 1/1/2/2/8 GiB AS caps, and h_rt 00:15/00:50/07:40/01:25/11:50. All tasks remain not_started with empty attempts. Priority-0 estimated points are 0.59065, with h_rt envelope 1.056. Official billing rules were checked and linked in the runbook.
- Added summary and 165×66 mm log-log figure scripts. They retain failures, exclusions and denominators, use a common complete input cohort at each length, and check CParty/exact PF parity where exact is measured. Crossings require adjacent valid measurements bracketing a ratio of one; missing results are not extrapolated.
- Validation: four new comparison tests and ten existing array tests pass. A 24/48/96-nt development array ran 36 CLI calls with two workers, passed PF parity and rendered PDF/SVG/PNG. Frozen-input, live-lock, memory-envelope, shell-syntax and unstarted-plan checks pass. Raw evidence is packed in tar files. The inode gate estimates 89 array lifecycle entries against 2,415,389 remaining below the conservative 24-million boundary; the full preparation has a 300-entry allowance. Scratch uses `$TMPDIR` and local intermediates are removed after inspection.

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

## 2026-10-15 E0 再実行の準備（d21618c5・20261014 固定版 artifact）

- 新しい oracle は固定 clone の build.sh で upstream `800e95a` + config.json の順序 0001/0002/0006/0007/0008 を適用し、`build/cparty-oracle-800e95a-5patches/` へ build（終了値0）。新しい oracle lock に manifest 全体と patch/source/binary/object・build log の hash を保存した。既存4パッチ oracle から新しい出力先へ Git objects を file:// clone した。旧 checkout、旧 lock と開発 tree は変更していない。
- Inside 専用の外部 observer を、pilot lock で検証した `build/linear-cparty-d21618c5/` の exact/RNA libraries と新 oracle の元 objects に対して build（終了値0）。出力は `build/e0-observers-d21618c5-800e95a-5patches/`、lock は `manifests/software-e0-d21618c5-20261014.json`。API 対応は不要だった。Library はリンク前後の hash と lock の実依存を検査した。
- 科学的な小確認：既存4証人と新しい WIP29（scaffold の (3,7) に3塩基の hairpin）の全5入力で `LegacyCompatible` と5パッチ oracle の exact inside log Z / PF free energy が一致。新15ジョブと旧WIPの2ジョブが成功。最大絶対差は log Z `2.1316282072803006e-14`、PF free energy `1.4210854715202004e-14 kcal/mol`。旧 `15aa7742` と4パッチ oracle のWIP値も一致し、新−旧は log Z 約 `+0.02822658740123`、PF free energy 約 `−0.01739663223684 kcal/mol`。Turner04、37°C、dangles2 と従来の両絶対許容を維持した。既定 closure の差はnested45の1/5入力で検出。BPP/MFEや全bpRNAの一致・性能の主張は追加していない。
- 新configは同じ300本のbpRNA cohort＋元4証人＋WIP29の305入力。新array configは199 nt以下cpu_16、200–399 nt cpu_40、workers=1。300–399 ntの旧128 GiB capを80,000 MiBへ変更し、task budget 81,920 MiBに収めた。既存の時間・memory envelopeは未校正の計画仮定として維持した。
- Dry-run は実依存lockと既存full G datasetを検証済み（`pending_scaffolds=false`）。21tasks・915jobsはすべて未開始。5scriptsは`bash -n`成功。共有出力はtaskごとにJSONL/tar、compact planとinode gateを有効にし、新計画297 entriesはlive quotaと保守的上限による残り2,951,789 entriesの10%以内。各CLIの小確認は`/tmp`で行い、1つのtarにまとめて一時ファイルと作り直した初回確認のarchiveを削除した。
- 開発上の検査：追加7テストを含む `python3 -m unittest discover -s tests` は全94テスト成功。旧/順序違いoracle、汚れたsource、WIP両実装の新旧差、同じ305入力、cpu_40に収まらない旧capの拒否を確認した。保存tarからの再解析も同値。実依存lock・plan/script hash・旧lock/config/cohort不変・`git diff --check`を確認した。実装側maintained testsは今回再実行していない。
- 要確認：scalar insideの短いtier実測による時間・RSS/cap/timeout校正（特に300–399 nt）、schedulerのresource/h_rt制限・利用group・同時task数・module、投入直前のinode/byte quotaとnode-local scratch容量。qsub、git commit、重い計算実験、開発treeの変更は行っていない。

[E0追記](full_bprna_array_plan_20261011.md) と [再現・証拠](e0_oracle_update_20261014/README.md) を参照。

## 2026-10-16 E2 メモリ内訳測定の準備（d21618c5）

- 外部 stage driver の既存 JSON を保ち、`positive_items`、VP/WMBP/CASE3 の候補・選択数、4種の probe 数、全族の visits/positive/candidate/selected を追加した。`--ranking outer|baseline`、`--outer-pilot-exponent 0.5|2/3`、一度だけ inside を呼ぶ `--skip-inside-repeat` を追加し、既定の二度呼び出す診断は維持した。
- 公開 `BeamConfig` と固定版 api/beam/replay を調査した。baseline は outer pilot/guide を省く。b=100 の pilot 幅は p=22 または10で、任意幅の指定はない。scaffold-only LogMax 前走は baseline でも残る。内部 replay に forcing の flag はあるが公開 config に経路がなく、前走・強制保持の停止や前走幅 max(1000,b) の変更を public API driver では指定できない。CASE3 祖先登録の byte 数も公開統計にはない。
- `software-pilot-d21618c5.json` と E0 lock の exact/RNA archive hash を照合して `build/stage-timing-d21618c5/` に buildし、新しい `manifests/software-stage-timing-d21618c5-20261016.json` に source/builder/両lock/archives/CLI/provenance を固定した。依存 source を変更せず、従来15aa7742 APIへの syntax-only compileも成功した。
- 元E2の入力tar全体のhashと7入力すべての元byteをストリーム照合した。同じS/G、b=h=d=100、outerの inside/mfe/outside_bpp/loop_profile＋baseline insideを各別プロセスで1回、warm-up 0・repeat 1とするconfigとarray configを作成した。cpu_40・1 worker/thread、49,152 MiBのプロセス address-space cap、task budget 51,200 MiB。1 RNAの5条件を1 taskで逐次実行する。
- `plan_array.py --dry-run --estimate-only` のみ実行した。7 tasks・35 jobs・共有111 entries（7×13＋20）。h_rt は3/2/2 tasksに2/4/6時間、timeoutは各job 900/2100/2700秒。既存tierの最大whole-CLI時間を切り上げ、0.001n秒と2倍の余裕を加えた見積もりは全体56,743.98秒（15.76時間）、要求h_rt合計26 task-hours。stage/baselineの校正値ではない。live使用21,069,547、保守的2,400万上限まで残り2,930,453、10%枠293,045に対しinode gate成功。estimateはJSON1個とdirectoryのみで、本番plan・script・入力コピーは作成していない。
- 開発上の検証：全103 testsが25.240秒で成功。2つの既存smoke入力×outer/baseline/outer p10×inside/mfe/outsideを別プロセスで計算し、CLIのMFE構造・丸め表示のenergy/PF・full-precision BPP、全族と総数、一度だけのinside呼出しを照合した。既存全stageのCLI値とpacked runner/resumeも成功。E0 archive改変はbuild出力作成前に拒否した。smoke出力は`/tmp`のみ。両固定dependency treeはcleanを維持した。
- 科学的な範囲：workはmain forwardの活動数であり、pilot・scaffold前走・CASE3登録のbyte数を含まない。profile APIにはwork統計がない。別processのRSS peakやbaselineとの差は加算可能な内訳ではなく、baselineは保持集合も変える。長入力の測定・メモリ増大原因の確定は未実施。本番plan、qsub、git commit、外部実装・原稿・文献の変更は行っていない。

[driverの追記](stage_timing_driver.md#d21618c5-memory-breakdown-preparation-2026-10-16)、
[測定config](../configs/e2_memory_breakdown_d21618c5_20261016.json)、
[array config](../configs/arrays/e2_memory_breakdown_d21618c5_20261016.json)、
[見積もり](e2_memory_breakdown_d21618c5_20261016/array_estimate.json) を参照。


## 2026-10-17 Major 2: fixed-scaffold pseudoknot membership audit

- Added observation-only exact consumers and aggregation/report scripts for the saved 200-RNA / 983-pair cohort. The pinned d02d261c implementation and build remain read-only. All 983 pairs have current-build membership evidence: 973 positive exact marginals and 10 pairs in one successfully evaluated complete reference. The >0.5 count remains 179/983 (18.21%). Singleton additions: 980 reachable, 0 rejected, 3 unfinished.
- Complete canonical annotations: 164/200 evaluable; 36 rejected, with 22 extension-crossing, 14 minimum-span and 1 completed-WMB closure violations (nonexclusive). The historical G plus saved [] pairs is evaluable for 198/200 RNA; one closure rejection and one unresolved combination remain. Eref−EMFE for complete evaluable references has median 11.0395 kcal/mol (Q1–Q3: 4.80025–20.212); MFE overlap 210/983 pairs. Energy components and distance strata are saved separately from implementation validation.
- One 297-nt RNA hit the 60-second limit in both the original run and an energy-only retry; computation was stopped. Complete JSON values were recovered without inventing unfinished values. Three singleton tests remain unmeasured; a cpu_4, one-task, 30-minute array proposal is saved and was not submitted. This RNA reuses previous-build BPP, while 199 RNA / 973 recalculated pair probabilities match the previous exact output exactly.
- Validation: 200/200 public evaluator MFE rechecks; 1579 reachable energy audits; maximum local-sum/root error 1.35e-13 kcal/mol and component-sum error 1.42e-14. Both linked archive hashes and CMakeCache match the verified software-stage-timing-d02d261c lock. Source remains clean; manuscript and external development checkout were not edited. Temporary files were removed; raw attempts are in one JSONL and two packed archives. No commit.
- Raw results:results/local/pk-membership-20261017/ (ignored). Versionable summaries, README with five-sentence SI proposal, failure inventory, provenance and array proposal: [results/pk-membership-20261017](../results/pk-membership-20261017/README.md).

## 2026-10-18 Full paper rerun preparation for 89baea52

- Prepared the [21-row numerical-source inventory and runbook](rerun_89baea52_20261018/README.md), 25 live-verified plans and 57 scripts. The standard matrices contain 399,000 unstarted CLI jobs. Including five standalone validation tasks, there are 980 array tasks and three local example tasks. No scheduler submission or paper measurement was made.
- Built a separate upstream 800e95a oracle with ordered patches 0001/0002/0006/0007/0008/0009, plus same-release inside-only E0, stage and membership observers. Kept the five-patch oracle and external development tree unchanged. Added parameter-specific DP09 and Subopt-inclusive example locks derived from the certified release build.
- Preserved E1 science and tier resources; the combined short matrix supplies its own exact reference. E2 retains cpu_40, up to ten workers and an 8 GiB address-space cap. NEAT1 major outputs, by-loop and all existing scaffold variants are included. Fresh membership analysis requires all 200 new observations and does not substitute old-build BPP.
- Prepared submission commands with at most three synchronous arrays, each capped at ten tasks by default, and separate local/analysis commands. The shared-file estimate is 15,168 over preparation/execution and a conservative 9,352 additional entries after planning. Live quota gates passed; all CLI temporaries remain node-local and raw observations use task JSONL/tar.
- Validation: 31 tests passed and nine actual-RNA stage tests were explicitly skipped; four C++ probe syntax checks, all job shell checks, E1 template/resource equality, live dependency hashes and unstarted task states passed. Adapted only the new validation copies to removed runtime APIs/dead VP aliases. Preserved the initial link failure and two planning setup failures in one tar. No git commit or manuscript edit.


## 2026-10-18 E1 long cpu_16 replan for 89baea52 (not submitted)

- Prepared `e1-bprna-long-89baea52-cpu16-20261018`: 60 tasks (7/19/24/10), cpu_16, eight workers, 4,096 MiB/process and 1,024 MiB controller allowance; cap envelope 33,792 MiB <= 35,095 MiB. h_rt is 1/2/3/4 h, totaling 157 cpu_16 allocation hours versus the retained 694-task / 7,088 cpu_40-hour plan. All 25,884 conditions, 1,438 inputs, frozen input tar/lock bytes, modes/programs/outputs and order policy/seed are preserved. Previous run and script hashes are unchanged.
- Development resource probes used only two RNA inputs at b400. Raccess on 4,381 nt completed in 595.17 s with GNU peak RSS 1.598 GiB under a 4 GiB address-space cap. BPP on this input timed out at 180/600 s; Raccess on the 3,995 nt historical RSS outlier timed out at 900 s. The latter two sampled RSS maxima were 1.626/1.721 GiB, with no observed memory failure; their completed-run peaks remain unknown. Retained every attempt, command and raw output in one pilot tar. No manuscript performance/accuracy result was added.
- Replaced length-only allowance with all 694 historical task_results.jsonl files: 25,884 successful CLI elapsed measurements, 536.03 h. Transfer factor 0.516349 comes from the completed paired probe; timeouts remain censored lower bounds. Add 0.5 s/job, multiply work by 1.5, sum beam-wise list bounds, divide by assumed worker efficiency 0.8 and add 120 s/task. Efficiency is not an eight-worker throughput measurement. Every predicted task fits 70% of h_rt; the maximum fraction is 0.6999059.
- Added optional observed-input / parallel packing to the planner while retaining default sequential packing. 24 relevant tests passed; all 60 saved schedules also passed an independent modeled-duration replay. Four job scripts and the sequential submission helper pass bash syntax checks. All formal tasks remain not_started with empty attempts and no run directory. The first preflight rejected the unexpanded cohort config before creating outputs; it now validates the frozen expanded config and saves the original cohort config unchanged apart from name/cap.
- File budget: 802 entries for the new plan/run lifecycle plus 50 for ancillary artifacts, against 2,445,015 remaining under the conservative 24-million ceiling at the final live check. Temporary CLI files are node-local; retained raw output is one tar. [Runbook](rerun_89baea52_20261018/README.md), [comparison/probes](rerun_89baea52_20261018/jobs.md), [receipts](rerun_89baea52_20261018/cpu16/summary.json). qsub, git commit, existing-run removal/overwrite, manuscript edits and external implementation edits were not performed.

## 2026-10-19 Analyses C/D: both local error tails and coordinate-based recovery

- Analysis C: read the complete saved E1 short cohort (24,625 RNA, exact and b=20/50/100/200) with checksum/command/domain and cross-program checks. Counted exact >=0.5 / finite <0.1 and exact <0.1 / finite >=0.5 events for all BPP, BPP crossing G, PseudoLoop, and 4-nt total accessibility. Saved 394,000 RNA/width/metric rows, 80 summary rows, witnesses, validation and an SI table/three-sentence English proposal in [e1-local-errors-89baea52-20261019](../results/local/e1-local-errors-89baea52-20261019/README.md).
- Scientific result C: at b=200 all such events were in the 1,095 RNA of 300–399 nt. G-crossing BPP had 27 missed pairs in 7 RNA (0.64%) and 2 created pairs in 2 RNA (0.18%); PseudoLoop had 80 missed bases in 15 RNA (1.37%) and no created events; 4-nt accessibility had 10 missed windows in 7 RNA (0.64%) and 8 created windows in 4 RNA (0.37%). All-BPP counts were 117/8 events in 10/5 RNA. The all-window maximum accessibility error was 0.966789044 (bpRNA_RFAM_35940, window 318–321, exact 0.004942600 / finite 0.971731644), exceeding the old exact>=0.1-selected maximum 0.878913843.
- Analysis D: defined reference canonical pairs and new MFE pairs by whether their coordinates cross at least one G pair; matched coordinates independently of bracket symbols. In 200 RNA, 57 canonical curly-bracket reference pairs in 23 RNA enter the target set (983 -> 1,040). Turner04/DP09 sensitivities become 210/1,040=20.19% / 110/1,040=10.58% (old 21.36% / 11.19%); fractions with p>0.5 become 179/1,040=17.21% / 80/1,040=7.69%. All 483 predicted crossing sets equal the old square-bracket coordinate sets, so PPV remains 210/850=24.71% / 110/413=26.63%. Saved comparison tables, per-pair changes, distance summaries and the SI replacement/English proposal in [pk-recovery-pairset-89baea52-20261019](../results/local/pk-recovery-pairset-89baea52-20261019/README.md).
- Scientific result D, longer RNA: two old square-bracket reference pairs in the 80-RNA cohort cross only noncanonical UU pairs removed from G, so the reference count becomes 809 (old 811); sensitivity is 40/809=4.944% and PPV remains 40/301=13.29%. The three rRNA retain their saved 24 reference targets, with 0 MFE matches, 14 predicted crossing pairs, no p>0.5, and nine reported zeros. The available rRNA archive contains G and canonical square/curly reference coordinates, but no full original DBN; the recheck applies the same geometry to these frozen targets and does not claim a complete original-annotation audit.
- Development verification: six focused definition tests passed. Independently reproduced all C counts/fractions from the compressed RNA table and checked all 98,500 accessibility maxima against the existing E1 scalar table (maximum absolute difference 1.11e-16). All 369,375 saved E1 jobs succeeded and 123,125 cross-program consistency checks passed. D reproduced all 480 old MFE rows and all 483 old probability vectors exactly; a second crossing implementation agreed on all 280 full short/mid references, and all eight recovery summaries were reproduced from the pair table. Both result directories contain provenance and validation JSON. `git diff --check` passed.
- Operational scope: only saved outputs were analyzed; no RNA CLI, qsub, tar extraction, manuscript edit, external implementation edit, index operation or git commit. Checked project quota first (21,544,854 / 28,000,000 inodes); fewer than 30 new small files/directories were needed, with about 6 MiB of aggregate output. The shared raw archives were read in memory and analysis logs stayed under TMPDIR until copied as one log per result directory. Existing unrelated working-tree changes were preserved.


## 2026-10-19 Analysis A: same-G Terc crossing/no-crossing comparison

- Added `scripts/analyze_terc_nopk_89baea52.py` and `scripts/report_terc_nopk_89baea52.py`; results, settings, four proposed manuscript sentences and a Figure 3 panel proposal are in `results/local/terc-nopk-89baea52-20261019/README.md`. No manuscript or implementation edit, commit or scheduler submission was made.
- Scientific results: exact window 123–126 accessibility is 0.8605129487361832 with crossings and 0.12286507219770805 without crossings; dominant joint components are PseudoLoop 0.8415475077311951 and Hairpin 0.11695022375854795, respectively. Exact ensemble free energies are −121.40250658126598 and −119.6652734677487 kcal/mol. All nine reference crossing pairs are forbidden without crossings, while their endpoints can form different pairs or enter ordinary loop contexts; per-base 40–160 and four endpoint-region summaries are saved.
- Method: frozen release exact outputs versus ViennaRNA 2.7.2 with every G pair enforced, same Turner04 parameter bytes, temperature, dangling ends and ordinary loop settings. Exact E/H/M/B+I context masses use constrained PF ratios; only B/I separation uses sampling, with 100,000 unconditional and two 100,000 conditional draws and recorded seeds/standard errors. Read `README.md` for software/energy-model correspondence and scope.
- Verification: two constraint interfaces give identical log Z; forced-G BPP error <=8.882e−16; context sum error <=2.354e−14; nine fixed-target public-API energy checks differ by <=5.684e−14 kcal/mol; 3,000 sampled structures agree with ViennaRNA native context labels. Source, parameter, binaries and libraries were unchanged across measurement. Raw samples and commands are in one tar.gz, all disposable files were in TMPDIR, and output checksums pass. The exploratory get_Z binding error is preserved and its value was not used.


## 2026-10-20 Example search: Figure 1C and NEAT1 (analysis only)

- Added `scripts/search_examples_20261020.py` and `scripts/report_example_search_20261020.py`. Reports and scientific PNG/PDF previews are in [figure1-candidates-20261020](../results/local/figure1-candidates-20261020/README.md) and [neat1-pk-search-20261020](../results/local/neat1-pk-search-20261020/README.md). No manuscript, existing figure, external implementation, scheduler or Git-index/commit operation was performed.
- Figure 1C: screened all 105 reference-G inputs of length 50–100 in the 200-RNA Turner04 exact recovery cohort. Reused saved cutoff=0 BPP; ran only exact CapR and 10,000 exact samples each for 67 inputs with crossing BPP in [0.2,0.8]. All 67 also have intermediate PseudoLoop positions. The 134 new CLI calls took 16.14 seconds in summed process wall time, with the pinned 89baea52 build and matching parameter/binary hashes. Recommended bpRNA_RFAM_31397 (61 nt): PseudoLoop max 0.608737, 15 intermediate positions, three crossing pairs 0.495965–0.517433, sample ranks 1/2 = 34.31%/14.26%. Explicitly recorded that its top frequency does not meet the illustrative 50–80% range; PDB_379 is the one screened alternative in that range (50.13%) but retains strong 0/1 contrast. The report gives five candidates, all requested metrics, mixed-context counts and visual notes.
- NEAT1: read 28 BPP/context raw artifacts for 14 saved conditions and scanned 25,901 overlapping 300-nt windows at stride 1. No crossing pairs have p>=0.05 or p>=0.1. The full-length maximum is NEAT1_1 G0.7/b400, pair 6–3732, p=0.0122941016222; nine remote pairs span 3,722–3,730 nt. The two terminal windows 1–300 and 3457–3756 are only limited display candidates, not a stable PK example: b200 has no exported crossing pair and the terminal PseudoLoop peak changes strongly. NEAT1_2 crossing probabilities are all below the 0.001 export cutoff at b100/200/400. G0.5/G0.9 variants exist only for NEAT1_1 at b100/200 and show no increase in exported crossing pairs. Missing variants are not treated as measured zeros. Local limits h=d=b change along with b.
- Validation: checked reference-G definitions and all 105 frozen inputs; checked frozen NEAT1 inputs and raw-member hashes; independently classified coordinate crossings and matched the prior full-length NEAT1 crossing set and probabilities exactly. Checked context sums, all sampled structures retaining G, and sample/BPP errors (selected maximum 0.009370). Raw new calculations are one packed tar.gz; all disposable files were under TMPDIR, with no extraction on shared storage. Start quota was 21,567,719 / 28,000,000 inodes; fewer than 50 new shared entries and under 10 MiB were needed. Output checksums and both repository diff whitespace checks passed; PNG previews were inspected.


## 2026-10-20 LinearKnot pin, release parity and centroid preparation

- Pinned LinearKnot `9c7d6e170b2b21a94a2e626d9939817e62898253` in `deps/LinearKnot-9c7d6e1` and built Release in `build/linearknot-9c7d6e1` with `-j8`, tests enabled and LTO disabled per the maintained test guide. The maintained suite passed 87 checks plus its dependency gate; three selected lock checks and five Python tests passed. Source/tree, CLI, parameters, effective flags and build provenance are in `manifests/software-linearknot-9c7d6e1.json` and `manifests/release-preparation-linearknot-9c7d6e1.json`.
- Regression: the same 14 frozen release inputs, b20/b200/exact, fold/BPP, profile, access1/4/8/16 total/by-loop, sample normal/unique with identical seed, MEA/centroid gamma1 and Subopt top8. All 756 effective processes succeeded after retrying two exact sampling comparisons at 16 GiB (the old CLI hit the initial 8 GiB address-space cap). All 760 process attempts, including two initial failures, are retained; strict agreement after removing elapsed-time lines is 294/378. Agreement after the two documented access metadata renames is 378/378. See [the full report and CLI mapping](linearknot_9c7d6e1_20261020/README.md). Raw outputs remain one tar plus one process JSONL under `results/local/linearknot-9c7d6e1-comparison-20261020`; source/build fingerprints were unchanged.
- Prepared [centroid replacements](centroid_89baea52_20261020/README.md) on the adopted 89baea52 build: NEAT1 gamma1 at b100/200/400 (cpu_40, three workers, 20 GiB/process; h_rt 22 minutes and 2 h 50 minutes), full local Terc exact/b200 and Figure 1 exact matrices with MEA replaced, and a separate E0 decoding suite adding 6561 centroid objective checks and 6561 CLI records against enumeration. Gamma1 is ordinary centroid; its objective differs from MEA. Existing inputs, scaffold, other outputs, seeds and original plans are preserved.
- Verification: all four new matrix tasks and the E0 full suite remain unstarted; dependencies, frozen plans and seven shell scripts passed checks. A small E0 preparation smoke compiled the new probe and passed 27 API optimum/score checks plus 27 CLI checks, maximum error 2.842170943040401e-14. Full E0 is estimated at about 40 minutes in a 2-hour cpu_16 allocation. Initial suite cwd/parameter failure and the successful full rerun are retained in one log tar.
- Operations: checked quota before clone (21,567,719 / 28,000,000 files), kept planner inode checks enabled and all transient per-process files in TMPDIR. The clone/build use 833/2369 descendant entries. No scheduler submission, git commit, manuscript edit or modification to either external implementation checkout was made. Existing unrelated edits were preserved.

## 2026-10-20 Author figure annotations

- Added `prepare_feedback_figure_data.py` and `make_figures_feedback_20261020.py`; updated output-role labels and scaffold-shape names in `make_figures_review.py`. New figure directories are `figure1_89baea52_20261020_fb`, `figure4_terc_89baea52_20261020_fb`, `main_89baea52_20261020_fb`, and `e1_89baea52_20261020_fb` under `results/figures/`, each with paired caption proposals, README and validation.
- Scientific outputs: seven local 9c7d6e1 processes computed the 61-nt RFAM31397 exact fold/BPP, profile, accessibility, 10,000 samples and gamma1 centroid, plus Terc exact/b200 centroid. Both Terc centroids contain all six outer reference pairs and none of the three inner pairs. Terc probabilities/samples and all E1/E2 values remain saved 89baea52 results. Small local data, commands and one raw tar are in `results/local/figure-feedback-20261020/`.
- Display changes: remove the separate MFE box and Terc D table; threshold pair arcs at p≥0.05 with a numeric opacity key; remove Internal hatching; mark 123–126 in Terc A/B; give small sampling co-occurrence bars and a centroid inset; define generation seeds and measurement ranges; label fold/profile/access by output role and the two multiloop shapes by name. Retain the annotated deletions. No new NEAT1 figure is produced.
- Verification: all seven final RNA processes succeeded; binary/parameter hashes match the release lock and sources are unchanged. Initial centroid parsing expected the wrong gamma delimiter; it was corrected to the actual semicolon and the calculations rerun. Context sums, scaffold retention, window counts and aligned coordinates passed. All seven PNGs were opened; actual PDF text is at least 7 pt at the stated 178/165-mm widths. E1/E2 plotting data are byte-identical to r4; 1,750 E1 quantile scalars and 350 summaries were checked. No scheduler submission, commit, manuscript edit or external implementation write.
- Handoff: use each new caption.md, including p=0.05 in place of the parallel manuscript draft's p=0.1. Terc is now 104 mm high. With NEAT1 removed, accessibility is S1 and synthetic shapes are S2; the extra threshold accuracy plot remains outside the manuscript.

## 2026-10-21 LinearKnot 4f97fd3 pin and resource rerun preparation

- Resumed the existing clean `deps/LinearKnot-4f97fd3` checkout and partial Release build. The maintained required suite passed 96 checks plus the dependency gate, followed by three lock checks. Saved actual compile/link settings, source/tree/cache/CLI/parameter hashes and the suite log archive. The temporary source/build link was removed.
- Compared all 14 frozen release inputs at b20/b200/exact, including fold/BPP, profile, access1/4/8/16 total/by-loop, seeded samples normal/unique, MEA, centroid and Subopt. All 378 comparisons were byte-identical after excluding only elapsed-time lines; all 756 processes succeeded without retries. No metadata rename was needed. Source/build fingerprints were unchanged.
- Prepared [12 resource matrices, 82 unstarted tasks and 30 scripts](linearknot_4f97fd3_20261021/README.md), totaling 2,768 CLI jobs including warm-ups. E2 retains cpu_40, up to ten workers and 8 GiB AS caps; memory stages retain their 0+1 policy. NEAT1 retains the table's outputs, widths, seeds and sample length restrictions. Added same-release stage/PF observers; PF retains exact through 1,000 nt and b100 through 4,000 nt, compatible closure and the same synthetic inputs. The existing CParty observations, including failures, are reused without execution.
- Added executable/subcommand mapping to the packed runner while preserving distinct analysis names and actual-binary lock checks. PF analysis can reuse only CParty rows from a separately specified matched plan. Thirty-seven workflow tests, 32 shell syntax checks, four PF observer comparisons and ten stage smoke processes passed. All dependencies and frozen input hashes were verified; previous observations and both external implementation trees were unchanged.
- Chose h_rt from prior completed allocation times at unchanged concurrency, with explicit padding and no assumed new-release speedup. Priority-0 planning cost is 6.3811 points; the full h_rt envelope is 10.0340 points. The 12-plan lifecycle gate covers 1,294 entries; the whole preparation lifecycle allowance is 8,000 entries. Final conservative remaining inode count was 2,398,936. Per-process files stayed in TMPDIR and retained raw outputs/logs are packed. No qsub, commit, manuscript/figure numerical update or existing-run overwrite occurred.
