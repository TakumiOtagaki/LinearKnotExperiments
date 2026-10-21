#!/usr/bin/env bash
# Alternative explicit commands, not executed. Use from experiments/.
set -euo pipefail
: "${TSUBAME_GROUP:?Set the authorized accounting group}"
LCP_ARRAY_BUDGET="${LCP_ARRAY_BUDGET:-30}"
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e2-rnacentral-4f97fd3-20261021/rna_4_8kb.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e2-rnacentral-4f97fd3-20261021/rna_8_20kb.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e2-rnacentral-4f97fd3-20261021/rna_20_30kb.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e2-synthetic-shapes-4f97fd3-20261021/synthetic_n250.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e2-synthetic-shapes-4f97fd3-20261021/synthetic_n500.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e2-synthetic-shapes-4f97fd3-20261021/synthetic_n1000.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e2-synthetic-shapes-4f97fd3-20261021/synthetic_n2000.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e2-synthetic-shapes-4f97fd3-20261021/synthetic_n4000.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e2-memory-breakdown-4f97fd3-20261021/rna_4_8kb.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e2-memory-breakdown-4f97fd3-20261021/rna_8_20kb.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e2-memory-breakdown-4f97fd3-20261021/rna_20_30kb.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e2-parallel-impact-w1-4f97fd3-20261021/parallel_impact.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e2-parallel-impact-w10-4f97fd3-20261021/parallel_impact.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e3-neat1-bpp-4f97fd3-20261021/neat1_short.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e3-neat1-bpp-4f97fd3-20261021/neat1_long.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e3-neat1-capr-4f97fd3-20261021/neat1_short.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e3-neat1-capr-4f97fd3-20261021/neat1_long.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e3-neat1-raccess-4f97fd3-20261021/neat1_short.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e3-neat1-raccess-4f97fd3-20261021/neat1_long.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e3-neat1-raccess-byloop-4f97fd3-20261021/neat1_short.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e3-neat1-raccess-byloop-4f97fd3-20261021/neat1_long.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e3-neat1-sampling-4f97fd3-20261021/neat1_short.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e3-neat1-sampling-4f97fd3-20261021/neat1_long.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e3-neat1-centroid-4f97fd3-20261021/neat1_short.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/e3-neat1-centroid-4f97fd3-20261021/neat1_long.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/linearknot-pf-4f97fd3-20261021/n250.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/linearknot-pf-4f97fd3-20261021/n500.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/linearknot-pf-4f97fd3-20261021/n1000.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/linearknot-pf-4f97fd3-20261021/n2000.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc "$LCP_ARRAY_BUDGET" jobs/linearknot-pf-4f97fd3-20261021/n4000.sh
