#!/usr/bin/env bash
set -euo pipefail
cd /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 LC_ALL=C
: "${TMPDIR:?node-local TMPDIR is required}"

wait_plan() {
  python3 -B - "$1" <<'WAITPY'
import json,pathlib,sys,time
plan_path=pathlib.Path(sys.argv[1]);plan=json.loads(plan_path.read_text())
while True:
 states=[json.loads((plan_path.parent/t['directory']/'task.json').read_text()).get('status') for t in plan['tasks']]
 if all(s=='success' for s in states):break
 if any(s in ('failed','interrupted') for s in states):raise SystemExit(f'failed task in {plan_path}: {states}')
 print(f'Waiting for {plan_path.parent.name}: {dict((s,states.count(s)) for s in set(states))}',flush=True)
 time.sleep(30)
WAITPY
}
# Analyze completed independent cohorts before waiting for NEAT1 variants.
wait_plan results/runs/e0-bprna-legacy-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/e0-bprna-legacy-89baea52-20261018/array_plan.json
wait_plan results/runs/e2-rnacentral-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/e2-rnacentral-89baea52-20261018/array_plan.json
wait_plan results/runs/e2-synthetic-shapes-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/e2-synthetic-shapes-89baea52-20261018/array_plan.json
wait_plan results/runs/e2-memory-breakdown-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/e2-memory-breakdown-89baea52-20261018/array_plan.json
wait_plan results/runs/e2-parallel-impact-w1-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/e2-parallel-impact-w1-89baea52-20261018/array_plan.json
wait_plan results/runs/e2-parallel-impact-w10-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/e2-parallel-impact-w10-89baea52-20261018/array_plan.json
wait_plan results/runs/pk-recovery-200-turner04-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/pk-recovery-200-turner04-89baea52-20261018/array_plan.json
wait_plan results/runs/pk-recovery-200-dp09-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/pk-recovery-200-dp09-89baea52-20261018/array_plan.json
wait_plan results/runs/pk-recovery-80-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/pk-recovery-80-89baea52-20261018/array_plan.json
wait_plan results/runs/pk-recovery-rrna-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/pk-recovery-rrna-89baea52-20261018/array_plan.json
wait_plan results/runs/terc-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/terc-89baea52-20261018/array_plan.json
wait_plan results/runs/terc-byloop-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/terc-byloop-89baea52-20261018/array_plan.json
wait_plan results/runs/pk-membership-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/pk-membership-89baea52-20261018/array_plan.json
wait_plan results/runs/figure1-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/figure1-89baea52-20261018/array_plan.json
python3 -B scripts/analyze_e0.py results/runs/e0-bprna-legacy-89baea52-20261018 --out results/local/e0-bprna-legacy-89baea52-20261018.json
python3 -B scripts/analyze_e0_e2_runs.py --revision 89baea52 --previous-revision d02d261c --workers 10 --previous-workers 10 --previous-resource cpu_40=1 --skip-cap30-history --rna-run results/runs/e2-rnacentral-89baea52-20261018 --synthetic-run results/runs/e2-synthetic-shapes-89baea52-20261018 --memory-run results/runs/e2-memory-breakdown-89baea52-20261018 --parallel-w1-run results/runs/e2-parallel-impact-w1-89baea52-20261018 --parallel-w10-run results/runs/e2-parallel-impact-w10-89baea52-20261018 --previous-rna-run results/runs/e2-rnacentral-d02d261c-20261017 --previous-synthetic-run results/runs/e2-synthetic-shapes-d02d261c-20261017 --previous-memory-run results/runs/e2-memory-breakdown-d02d261c-20261017 --out results/local/e2-89baea52-20261018
python3 -B scripts/validate_e0_e2_runs.py results/local/e2-89baea52-20261018 --revision 89baea52 --previous-revision d02d261c --test-log docs/rerun_89baea52_20261018/preparation_tests.log
python3 -B scripts/analyze_small_rerun_89baea52.py recovery
python3 -B scripts/analyze_small_rerun_89baea52.py membership
python3 -B scripts/analyze_small_rerun_89baea52.py examples
wait_plan results/runs/e3-neat1-bpp-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/e3-neat1-bpp-89baea52-20261018/array_plan.json
wait_plan results/runs/e3-neat1-capr-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/e3-neat1-capr-89baea52-20261018/array_plan.json
wait_plan results/runs/e3-neat1-mea-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/e3-neat1-mea-89baea52-20261018/array_plan.json
wait_plan results/runs/e3-neat1-raccess-byloop-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/e3-neat1-raccess-byloop-89baea52-20261018/array_plan.json
wait_plan results/runs/e3-neat1-raccess-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/e3-neat1-raccess-89baea52-20261018/array_plan.json
wait_plan results/runs/e3-neat1-sampling-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/e3-neat1-sampling-89baea52-20261018/array_plan.json
wait_plan results/runs/e3-neat1-variants-bpp-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/e3-neat1-variants-bpp-89baea52-20261018/array_plan.json
wait_plan results/runs/e3-neat1-variants-capr-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/e3-neat1-variants-capr-89baea52-20261018/array_plan.json
wait_plan results/runs/e3-neat1-variants-raccess-89baea52-20261018/array_plan.json
python3 -B scripts/aggregate_array.py --plan results/runs/e3-neat1-variants-raccess-89baea52-20261018/array_plan.json
python3 -B scripts/analyze_neat1_89baea52.py
