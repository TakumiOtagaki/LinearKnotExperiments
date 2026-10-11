#!/usr/bin/env bash
#$ -cwd
#$ -S /bin/bash
#$ -j y
#$ -N g_nt3200_4381
#$ -l cpu_4=1
#$ -l h_rt=01:00:00
#$ -t 1-4
#$ -o /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments/results/runs/g-bprna-full-plan-20261011-v3/scheduler_logs
set -euo pipefail
: "${SGE_TASK_ID:?Run as an array task}"
: "${JOB_ID:?Missing allocation}"
cd /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec "${LCP_PYTHON:-python3}" scripts/scaffold_array.py run --plan /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments/results/runs/g-bprna-full-plan-20261011-v3/array_plan.json --tier nt3200_4381 "$@"
