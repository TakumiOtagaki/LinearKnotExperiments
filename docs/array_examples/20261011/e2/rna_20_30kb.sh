#!/usr/bin/env bash
#$ -cwd
#$ -S /bin/bash
#$ -j y
#$ -N lcp_rna_20_30kb
#$ -l node_f=1
#$ -l h_rt=05:00:00
#$ -t 1-2
#$ -o /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments/results/runs/array-example-e2-20261011/scheduler_logs

set -euo pipefail
: "${SGE_TASK_ID:?Run as an SGE array task}"
: "${JOB_ID:?Missing scheduler allocation}"
cd /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# Optional first argument: --resume (submit only the interrupted task range).
exec "${LCP_PYTHON:-python3}" scripts/run_array_task.py \
  --plan /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments/results/runs/array-example-e2-20261011/array_plan.json \
  --tier rna_20_30kb "$@"
