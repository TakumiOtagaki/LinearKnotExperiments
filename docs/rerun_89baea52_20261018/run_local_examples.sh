#!/usr/bin/env bash
set -euo pipefail
cd /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments
unset SGE_TASK_ID
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 LC_ALL=C
: "${TMPDIR:?node-local TMPDIR is required}"

python3 -B scripts/run_array_task.py --plan results/runs/terc-89baea52-20261018/array_plan.json --tier terc --task-id 1 --software-root "$PWD/deps/LinearCParty-89baea52" --build-dir "$PWD/build/linear-cparty-89baea52" --scratch-dir "$TMPDIR"
python3 -B scripts/run_array_task.py --plan results/runs/terc-byloop-89baea52-20261018/array_plan.json --tier terc --task-id 1 --software-root "$PWD/deps/LinearCParty-89baea52" --build-dir "$PWD/build/linear-cparty-89baea52" --scratch-dir "$TMPDIR"
python3 -B scripts/run_array_task.py --plan results/runs/figure1-89baea52-20261018/array_plan.json --tier figure1 --task-id 1 --software-root "$PWD/deps/LinearCParty-89baea52" --build-dir "$PWD/build/linear-cparty-89baea52" --scratch-dir "$TMPDIR"
