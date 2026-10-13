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
