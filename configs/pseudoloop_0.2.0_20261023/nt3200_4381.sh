#!/usr/bin/env bash
#$ -cwd
#$ -S /bin/bash
#$ -j y
#$ -N pl020_nt3200_4381
#$ -l cpu_16=1
#$ -l h_rt=02:20:00
#$ -t 1-12
#$ -tc 30
#$ -o /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments/results/runs/pseudoloop-context-0.2.0-20261023/scheduler_logs

set -euo pipefail
cd /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 LC_ALL=C
: "${TMPDIR:?node-local TMPDIR is required}"
exec "${LCP_PYTHON:-python3}" -B scripts/pseudoloop_context_array_20261023.py run --plan /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments/results/runs/pseudoloop-context-0.2.0-20261023/array_plan.json --tier nt3200_4381 --task-id "${SGE_TASK_ID}" "$@"
