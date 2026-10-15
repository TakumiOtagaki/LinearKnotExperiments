# Public API stage timing observer

This E2 development tool links the fixed `15aa7742` beam, exact and RNA
archives without editing LinearCParty. Only small functional checks were run.
Formal measurements wait for the pinned prefix-ranking release.

## Build and use

From the experiment repository:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 scripts/build_stage_timing_driver.py
export CPARTY_FIXED_ENERGY_PARAM_FILE="$PWD/deps/LinearCParty-15aa7742/params/rna_Turner04.par"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 LC_ALL=C
build/stage-timing-15aa7742/stage_timing_driver \
  -i data/smoke/pk_scaffold.txt -b 20 --max-hairpin 30 --max-case3-depth 30 \
  --access-len 1,2,4,7 -k 10 --seed 20261013 --gamma 1 --session-mode both
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover -s tests -v
```

The builder checks the clean source, dedicated CLI build and parameter against
`software-pilot-15aa7742.json`, and exact/RNA archives against the E0 lock.
It records the beam archive hash from that verified build. Its generated lock
freezes archives, base lock, CLI binaries, driver source, build script and
provenance. Existing build directories and locks are refused.
For a new release, supply its `--software-root`, `--build-dir`, `--base-lock`,
`--out` and `--lock`; rebuild rather than editing a lock.

Supply `-i ONE_RECORD` or `--sequence S --scaffold G`. The file reader accepts
comments, blank lines and an optional name, and normalizes lowercase/T as the
CLIs do; multiple records are refused. Set b/h/d explicitly: API defaults here
are 100/30/30, while the planned paper setting is h=d=max(b,30).
The model uses the thermodynamic default 37 C and `RequireWeakClosure`.

## Measurement boundaries and reuse

| Stage | Public API | Work included |
|---|---|---|
| `inside` | `beam_inside` / `BeamSession::inside` | Preparation and numeric forward |
| `inside_repeat` | Same API again | Fresh forward or cached result copy |
| `mfe` | `beam_viterbi` / `BeamSession::viterbi` | Independent LogMax forward and traceback |
| `outside_bpp` | `beam_outside` / `BeamSession::outside` | Forward, outside and BPP collection |
| `loop_profile` | `beam_structural_profiles(...,{})` | Forward, outside and context aggregation |
| `accessibility` | `beam_structural_profiles(...,lengths)` | Forward, outside, contexts and joint windows |
| `sampling` | `beam_sample` | Factorized forward, sampler setup and draws; no outside |
| `mea` | `beam_mea` | Factorized forward, outside, marginals and decoding |

`fresh` uses free functions that construct new sessions. `shared` uses one
`BeamSession` for inside/MFE/outside; other APIs remain free functions.
`both` puts both runs in one JSON line for value checks. Only the inside result
and work counters are cached. **There is no public graph handoff.** Other calls
rebuild even with the same session (CLI と同じく毎回作り直す).
The cached inside repeat measures a result copy; its `api_elapsed_seconds` is
null because the stored API timer belongs to the original forward pass.
MFE ranks its own LogMax beam, distinct from the inside-ranked ensemble.

This follows the [API README](../deps/LinearCParty-15aa7742/src/scfg2/api/README.md),
[session declaration](../deps/LinearCParty-15aa7742/src/scfg2/api/beam.hh) and
[beam](../deps/LinearCParty-15aa7742/src/scfg2/api/beam.cc),
[profile](../deps/LinearCParty-15aa7742/src/scfg2/api/profiles.cc),
[sampling](../deps/LinearCParty-15aa7742/src/scfg2/api/sampling.cc) and
[MEA](../deps/LinearCParty-15aa7742/src/scfg2/api/mea.cc) implementations.
The session stores `inside_cache_`; other calls allocate their own charts.

`steady_clock` surrounds each whole public call, including returned-result
assembly and API-local destruction. It excludes RSS reads, driver summaries
and stdout serialization. MEA's `api_decoding_seconds` is its internal timer
for rewards, maximization and output after marginals are available.
Outside-only, aggregation-only and draw-only timers are unavailable; differences
between whole calls do not establish those times. More public timing fields or
a prepared-graph API would be needed for that decomposition.

Before/after `/proc/self/status` readings give current `vmrss_kib` and
**process-lifetime** `vmhwm_kib`. Returned results are alive at the after point;
API-local graphs have been destroyed. Allocators, prior calls, cached results
and driver summaries affect later readings. RSS differences are not additive
stage costs. Use separate processes with `--stage outside_bpp` (or another
stage) for an isolated call's lifetime high-water mark. `--stage inside` also
runs the separate cache diagnostic. Use separate fresh/shared invocations for
resource comparisons; `both` carries memory history into the second run.

## JSON and array execution

One strict JSON line contains revision, input/config, reuse contract and runs.
Each stage has scope, `ok/reachable/error`, external/API timers, RSS points,
`graph_reused=false`, a cache flag and numerical summaries: log Z/free energy,
MFE structure/energy, BPP/unpaired sums, context sums/normalization residual,
window sums/counts, sample counts, and MEA score/forest counts. CASE3-truncated
endpoint counters are included where exposed. Nonfinite numbers become null;
`reachable` is null for APIs without that field. Exit codes are 0 for success,
1 for failed/unreachable API results, and 2 for invalid arguments.

`--include-values` exports complete BPP, contexts, joint accessibility and
sampled structures for small checks. It increases memory held between stages.
BPP rows are `[i,j,round,square,union]` with 1-based coordinates; arrays omit
API index zero. Context columns are Stem/Hairpin/Bulge/Internal/Multi/Exterior/
PseudoLoop. Window arrays cover starts 1..n-w+1; w>n has zero windows.

The [matrix](../configs/stage_timing_smoke_15aa7742.json) and
[array template](../configs/arrays/stage_timing_smoke.json) use the existing
runner: two tiny inputs, fresh/shared modes, workers=1, no warm-up and one repeat.
The template uses `cpu_4` and an uncalibrated cost model.

```bash
python3 scripts/plan_array.py --config configs/stage_timing_smoke_15aa7742.json \
  --array-config configs/arrays/stage_timing_smoke.json \
  --software-root deps/LinearCParty-15aa7742 --build-dir build/stage-timing-15aa7742 \
  --out /tmp/stage-timing-plan --jobs-dir /tmp/stage-timing-jobs
python3 scripts/run_array_task.py --plan /tmp/stage-timing-plan/array_plan.json \
  --tier smoke --task-id 1 --scratch-dir /tmp
```

The planner checks project inode quota and creates a wrapper without submitting.
Per-job files stay in node-local scratch. Each task stores one JSONL and one tar,
plus constant-size runner metadata; driver JSON is each archived `stdout.txt`.
Read tar members without extracting them. Existing checkpoint/resume behavior
preserves failures and skips completed jobs. The saved four-condition smoke is
in `results/local/stage-timing-smoke-20261013/` (JSONL, tar and a small report).
No qsub or commit was performed.

For E2, pin/rebuild/check the prefix release, create a new config/lock, set adopted
b/h/d and output sizes, warm-up 1/repeats 3 with seeded block order, and calibrate
resource tiers. Isolated calls can be separate modes with `--stage`. The smoke
memory/time budget is not evidence of long-input feasibility.

## Verification

All **87 tests passed**, including seven new observer integration tests. Small
24/45/144 nt inputs match LinearCParty, LinearPKCapR, LinearPKRaccess,
LinearPKSampling and LinearPKMEA. BPP/profile/window values and MEA scores use
full-precision comparisons; MFE/PF free energy respect CLI display rounding
(0.00501/0.00005001 kcal/mol). MEA stats also check full-precision log Z.
Fresh/shared summaries and cached inside values agree. JSON/RSS, invalid inputs,
zero samples, overlong windows, unreachable roots, build rejection, packed task execution and resume
are checked. Build first: integration tests skip if dependencies are absent.
These checks do not establish ensemble semantics or performance scaling.

## d21618c5 memory breakdown preparation (2026-10-16)

The sections above describe the historical 15aa7742 observer. The new observer
links the clean `d21618c5115453174204329fcc40534b1cb9aa9c` dedicated build.
The new [formal config](../configs/e2_memory_breakdown_d21618c5_20261016.json)
collects isolated whole-call observations; no long-input measurement was run
during this preparation.

```bash
PYTHONDONTWRITEBYTECODE=1 python3 scripts/build_stage_timing_driver.py \
  --software-root deps/LinearCParty-d21618c5 \
  --build-dir build/linear-cparty-d21618c5 \
  --base-lock manifests/software-pilot-d21618c5.json \
  --archive-lock manifests/software-e0-d21618c5-20261014.json \
  --out build/stage-timing-d21618c5 \
  --lock manifests/software-stage-timing-d21618c5-20261016.json --purpose formal
export CPARTY_FIXED_ENERGY_PARAM_FILE="$PWD/deps/LinearCParty-d21618c5/params/rna_Turner04.par"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 LC_ALL=C
build/stage-timing-d21618c5/stage_timing_driver \
  -i data/smoke/pk_scaffold.txt -b 100 --max-hairpin 100 --max-case3-depth 100 \
  --ranking outer --session-mode fresh --stage inside --skip-inside-repeat
```

The builder verifies the base lock before and after linking, checks exact/RNA
archive hashes against the same-revision E0 lock, and freezes both locks,
archives, CLI binaries, driver source, builder and provenance in the
[new software lock](../manifests/software-stage-timing-d21618c5-20261016.json).
`--archive-lock` is inferred for 15aa7742 and d21618c5 when omitted.
`--purpose` defaults to `development`; this build explicitly uses `formal`.
Existing output directories and locks remain protected from overwriting.
The pre-ranking 15aa7742 API still passes a syntax-only compilation check.
Historical locks were not rewritten to accept the edited observer source.

### Added JSON fields and options

All previous JSON keys and the four previous `work` entries remain.
For inside, MFE and outside/BPP, `summary.work` now exports every public
`BeamWorkStats` field. The existing total aliases are `item_visits`,
`frontier_candidates` and `frontier_selected`; the added total alias is
`positive_items`. Additional fields are:

- `vp_candidates`, `vp_selected`, `wmbp_candidates`, `wmbp_selected`,
  `case3_candidates`, `case3_selected`;
- `binary_join_probes`, `extension_probes`, `v_parent_probes`,
  `vp_loop_parent_probes`;
- `families`, an array of objects with `nonterminal`, `item_visits`,
  `positive_items`, `frontier_candidates` and `frontier_selected`, including
  families whose counts are zero.

`case3_truncated_endpoints` keeps its previous meaning. These counters summarize
the **main forward DP**; they exclude scaffold-only prepass and pilot/guide work.
They count activity, not live objects, container capacity or allocated bytes.
In particular, CASE3 candidates are not branching-ancestor registration counts.
`StructuralProfileResult` exposes no `BeamWorkStats`, so loop-profile JSON
continues to contain numerical profile summaries and RSS/timing, without work
counters. Counter absence there does not mean zero work.

New driver options are:

| Option | Behavior |
|---|---|
| `--ranking outer\|baseline` | Sets public `BeamConfig::ranking`; default outer on d21618c5 |
| `--outer-pilot-exponent 0.5\|2/3` | Selects the pilot width rule; default 2/3, explicit option requires outer |
| `--skip-inside-repeat` | Runs one inside call; omits the second diagnostic call |

`config` adds `ranking`, `outer_pilot_exponent` (null for baseline) and
`inside_repeat`. Omitting `--skip-inside-repeat` preserves the old diagnostic
behavior: fresh executes two forwards, while shared copies the cached inside
result on its second call. The formal config uses fresh plus this option,
so each process performs exactly one selected API call. Full-value export is
disabled. On the old API the driver defaults to baseline and rejects outer
options rather than silently ignoring them.

### Preparation controls and limits of the comparison

Inspection covered the fixed [public config](../deps/LinearCParty-d21618c5/src/scfg2/api/beam.hh),
[API adapter](../deps/LinearCParty-d21618c5/src/scfg2/api/beam_detail.hh),
[ranking preparation](../deps/LinearCParty-d21618c5/src/scfg2/beam/prepare_ranking.hh),
[pilot width calculation](../deps/LinearCParty-d21618c5/src/scfg2/beam/pilot_outer.hh)
and [online replay](../deps/LinearCParty-d21618c5/src/scfg2/replay/w_final_vp_wmbp_online.hh).

| Preparation | Public control on d21618c5 | Consequence |
|---|---|---|
| (d) LogSumExp pilot and bounded outer guide | `ranking=Baseline` skips both; `outer_pilot_exponent` chooses p=ceil(sqrt(b)) or ceil(b^(2/3)) | At b=100, p=10 or 22; guide frontier is 4p. No arbitrary pilot-width field exists. |
| (e) Scaffold-only LogMax prepass and forced retention | No control in public `BeamConfig` or CLI settings | Baseline still runs the prepass. Its width is internally max(1000,b) for these finite uniform beams. |

Replay's internal widths type has `force_scaffold_derivation` and `forced_items`.
The public adapter does not expose them, and the scaffold-only function computes
its own width. No public-API driver option was added for disabling the prepass,
changing its width or dropping forced retention. Replacing the public calls
with internal replay calls would change this observer's measurement contract.
Outer prepares the forced set before its pilot; baseline obtains it from the
online runner. Thus baseline disables (d) while retaining (e), but it also
changes main-beam selection. RSS subtraction does not isolate guide bytes.

Unlike 15aa7742, a shared d21618c5 `BeamSession` caches immutable guide scores
and forced keys as well as the scalar inside result. It still does not hand
off the main graph. Fresh free-function calls reconstruct their own preparation,
and standalone profiles prepare their own ranking. All formal conditions use
separate fresh processes, so this guide reuse does not cross observations.

The main replay constructs CASE3 child indexes using the branching-ancestor
depth cap, including
[endpoint paths](../deps/LinearCParty-d21618c5/src/scfg2/beam/vp_case3_border_event_plan.hh)
and the
[nonleft index](../deps/LinearCParty-d21618c5/src/scfg2/beam/vp_case3_nonleft_last_child_index.hh).
Their internal path/registration sizes and byte allocations are not returned
by `BeamWorkStats`. This preparation can help distinguish (a), (b) and whole
(c) by process peak RSS and main-family activity, but cannot assign bytes to
(e), (f), individual containers, outside alone or BPP aggregation alone.
Lifetime RSS peaks and stage differences are not additive memory costs.

### Formal matrix and estimate only

The seven inputs match both the saved E2 plan hashes and the bytes streamed
from its hash-verified `inputs.tar`; no input was regenerated or extracted.
Each input has four outer conditions (`inside`, `mfe`, `outside_bpp`,
`loop_profile`) and baseline `inside`: **35 isolated processes in 7 tasks**.
All use b=h=d=100, 37 C, RequireWeakClosure, warm-ups 0, repetition 1,
cpu_40, one worker and one thread, with a blocked order seed of 20261016.

The [array config](../configs/arrays/e2_memory_breakdown_d21618c5_20261016.json)
keeps one RNA's five calls sequential in one task. Each process has a 49,152 MiB
(48 GiB) `RLIMIT_AS` address-space cap; this is not an RSS limit requested from
the scheduler. Task memory budget is 51,200 MiB including 1,024 MiB overhead.
The planning memory envelope is 40,960 MiB per process, above the historical
35.7734 GiB peak RSS; address-space use and baseline retention remain unmeasured
on long inputs.

| Length tier | Inputs/tasks | Job timeout | h_rt per task | Estimated five-call task time |
|---|---:|---:|---:|---:|
| 4–8 kb | 3 | 900 s | 02:00:00 | 3,647.93–3,669.90 s |
| 8–20 kb | 2 | 2,100 s | 04:00:00 | 9,618.32–9,636.52 s |
| 20–30 kb | 2 | 2,700 s | 06:00:00 | 13,225.43–13,291.06 s |

The budget uses the
[2026-10-16 E2 note](../../textbook/notes/e0_e2_d21618c5_20261016.md):
the largest LinearCParty wall times by tier were 350.41, 937.37 and 1,293.19 s.
Every condition receives a rounded whole-CLI maximum (360/950/1300 s), plus
0.001n seconds, then a factor of two. This is a conservative scheduling envelope,
not a calibrated prediction for isolated stages or baseline. The sequential
sum is 56,743.98 s (15.76 h); summed requested h_rt is 26 task-hours. Even five
job timeouts fit each task's h_rt with runner headroom.

Only the following estimate was created, with actual dependency verification:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 scripts/plan_array.py \
  --config configs/e2_memory_breakdown_d21618c5_20261016.json \
  --array-config configs/arrays/e2_memory_breakdown_d21618c5_20261016.json \
  --software-root deps/LinearCParty-d21618c5 --build-dir build/stage-timing-d21618c5 \
  --out docs/e2_memory_breakdown_d21618c5_20261016 --dry-run --estimate-only
```

The [estimate](e2_memory_breakdown_d21618c5_20261016/array_estimate.json) reports
no pending scaffolds and **111 planned shared entries = 7 × 13 + 20**, including
directories, scheduler logs and temporary-file headroom. Live quota was
21,069,547 / 28,000,000 files; the planner's conservative 24,000,000 cap leaves
2,930,453, with a 10% allowance of 293,045. The inode gate passed. Per-task
results remain one JSONL and one tar, plus constant-size metadata. Estimate-only
created one JSON and its directory; it created no executable plan, task
directories, schedules, input copies or scheduler scripts.

### Development verification

All **103 tests passed** (`PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover
-s tests -v`, 25.240 s). The observer tests now use the d21618c5 build and new
lock while retaining the historical smoke config as a fixture. New checks run
inside, MFE and outside/BPP as separate processes on both existing smoke inputs,
with outer, baseline and outer p=10. They verify complete work fields, family
totals, one-call selection, CLI MFE structure/rounded energy, PF rounded free
energy and full-precision BPP. Existing all-stage CLI/profile/accessibility/
sampling/MEA parity and packed runner/resume checks passed. A changed E0 archive
hash is rejected before an output directory is created.

Smoke outputs and test logs stayed in `/tmp`; temporary test files were cleaned.
The retained log is `/tmp/e2-memory-breakdown-tests-20261016.log` for this session.
Both fixed dependency trees remain clean. No implementation source, manuscript
or bibliography was edited, and no production plan, qsub or git commit was made.
