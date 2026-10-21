#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
: "${TSUBAME_GROUP:?Set the authorized accounting group}"
export PYTHONDONTWRITEBYTECODE=1
LCP_ARRAY_BUDGET="${LCP_ARRAY_BUDGET:-30}"
[[ "$LCP_ARRAY_BUDGET" =~ ^[0-9]+$ ]] && (( LCP_ARRAY_BUDGET >= 1 && LCP_ARRAY_BUDGET <= 30 )) || { echo "LCP_ARRAY_BUDGET must be 1..30" >&2; exit 2; }
python3 -B scripts/prepare_linearknot_timing_4f97fd3.py --verify
lfs quota -p 102571 /gs/bs
jobs=(
  'jobs/e2-rnacentral-4f97fd3-20261021/rna_4_8kb.sh'
  'jobs/e2-rnacentral-4f97fd3-20261021/rna_8_20kb.sh'
  'jobs/e2-rnacentral-4f97fd3-20261021/rna_20_30kb.sh'
  'jobs/e2-synthetic-shapes-4f97fd3-20261021/synthetic_n250.sh'
  'jobs/e2-synthetic-shapes-4f97fd3-20261021/synthetic_n500.sh'
  'jobs/e2-synthetic-shapes-4f97fd3-20261021/synthetic_n1000.sh'
  'jobs/e2-synthetic-shapes-4f97fd3-20261021/synthetic_n2000.sh'
  'jobs/e2-synthetic-shapes-4f97fd3-20261021/synthetic_n4000.sh'
  'jobs/e2-memory-breakdown-4f97fd3-20261021/rna_4_8kb.sh'
  'jobs/e2-memory-breakdown-4f97fd3-20261021/rna_8_20kb.sh'
  'jobs/e2-memory-breakdown-4f97fd3-20261021/rna_20_30kb.sh'
  'jobs/e2-parallel-impact-w1-4f97fd3-20261021/parallel_impact.sh'
  'jobs/e2-parallel-impact-w10-4f97fd3-20261021/parallel_impact.sh'
  'jobs/e3-neat1-bpp-4f97fd3-20261021/neat1_short.sh'
  'jobs/e3-neat1-bpp-4f97fd3-20261021/neat1_long.sh'
  'jobs/e3-neat1-capr-4f97fd3-20261021/neat1_short.sh'
  'jobs/e3-neat1-capr-4f97fd3-20261021/neat1_long.sh'
  'jobs/e3-neat1-raccess-4f97fd3-20261021/neat1_short.sh'
  'jobs/e3-neat1-raccess-4f97fd3-20261021/neat1_long.sh'
  'jobs/e3-neat1-raccess-byloop-4f97fd3-20261021/neat1_short.sh'
  'jobs/e3-neat1-raccess-byloop-4f97fd3-20261021/neat1_long.sh'
  'jobs/e3-neat1-sampling-4f97fd3-20261021/neat1_short.sh'
  'jobs/e3-neat1-sampling-4f97fd3-20261021/neat1_long.sh'
  'jobs/e3-neat1-centroid-4f97fd3-20261021/neat1_short.sh'
  'jobs/e3-neat1-centroid-4f97fd3-20261021/neat1_long.sh'
  'jobs/linearknot-pf-4f97fd3-20261021/n250.sh'
  'jobs/linearknot-pf-4f97fd3-20261021/n500.sh'
  'jobs/linearknot-pf-4f97fd3-20261021/n1000.sh'
  'jobs/linearknot-pf-4f97fd3-20261021/n2000.sh'
  'jobs/linearknot-pf-4f97fd3-20261021/n4000.sh'
)
failed=0
for job in "${jobs[@]}"; do
  qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" "$job" || failed=1
done
exit "$failed"
