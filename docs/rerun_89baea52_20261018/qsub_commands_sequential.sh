#!/usr/bin/env bash
set -euo pipefail
cd /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments
: "${TSUBAME_GROUP:?Set the authorized accounting group}"
export PYTHONDONTWRITEBYTECODE=1

# Alternative to submit_arrays.sh. Do not execute both submission files.
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e0-bprna-legacy-89baea52-20261018/nt1_49.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e0-bprna-legacy-89baea52-20261018/nt50_99.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e0-bprna-legacy-89baea52-20261018/nt100_199.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e0-bprna-legacy-89baea52-20261018/nt200_299.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e0-bprna-legacy-89baea52-20261018/nt300_399.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e1-bprna-short-89baea52-20261018/nt50_99.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e1-bprna-short-89baea52-20261018/nt100_199.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e1-bprna-short-89baea52-20261018/nt200_299.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e1-bprna-short-89baea52-20261018/nt300_399.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e1-bprna-long-89baea52-20261018/nt400_799.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e1-bprna-long-89baea52-20261018/nt800_1599.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e1-bprna-long-89baea52-20261018/nt1600_3199.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e1-bprna-long-89baea52-20261018/nt3200_4381.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e2-rnacentral-89baea52-20261018/rna_4_8kb.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e2-rnacentral-89baea52-20261018/rna_8_20kb.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e2-rnacentral-89baea52-20261018/rna_20_30kb.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e2-synthetic-shapes-89baea52-20261018/synthetic_n250.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e2-synthetic-shapes-89baea52-20261018/synthetic_n500.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e2-synthetic-shapes-89baea52-20261018/synthetic_n1000.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e2-synthetic-shapes-89baea52-20261018/synthetic_n2000.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e2-synthetic-shapes-89baea52-20261018/synthetic_n4000.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e2-memory-breakdown-89baea52-20261018/rna_4_8kb.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e2-memory-breakdown-89baea52-20261018/rna_8_20kb.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e2-memory-breakdown-89baea52-20261018/rna_20_30kb.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e2-parallel-impact-w1-89baea52-20261018/parallel_impact.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e2-parallel-impact-w10-89baea52-20261018/parallel_impact.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e3-neat1-bpp-89baea52-20261018/neat1_short.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e3-neat1-bpp-89baea52-20261018/neat1_long.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e3-neat1-capr-89baea52-20261018/neat1_short.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e3-neat1-capr-89baea52-20261018/neat1_long.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e3-neat1-mea-89baea52-20261018/neat1_short.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e3-neat1-mea-89baea52-20261018/neat1_long.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e3-neat1-raccess-byloop-89baea52-20261018/neat1_short.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e3-neat1-raccess-byloop-89baea52-20261018/neat1_long.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e3-neat1-raccess-89baea52-20261018/neat1_short.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e3-neat1-raccess-89baea52-20261018/neat1_long.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e3-neat1-sampling-89baea52-20261018/neat1_short.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e3-neat1-sampling-89baea52-20261018/neat1_long.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e3-neat1-variants-bpp-89baea52-20261018/neat1_short.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e3-neat1-variants-bpp-89baea52-20261018/neat1_long.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e3-neat1-variants-capr-89baea52-20261018/neat1_short.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e3-neat1-variants-capr-89baea52-20261018/neat1_long.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e3-neat1-variants-raccess-89baea52-20261018/neat1_short.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e3-neat1-variants-raccess-89baea52-20261018/neat1_long.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/pk-recovery-200-turner04-89baea52-20261018/pk200.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/pk-recovery-200-dp09-89baea52-20261018/pk200.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/pk-recovery-80-89baea52-20261018/pk80.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/pk-recovery-rrna-89baea52-20261018/rrna.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/pk-membership-89baea52-20261018/membership.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e0-validation-89baea52-20261018/enumeration.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e0-validation-89baea52-20261018/decoding.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e0-validation-89baea52-20261018/byloop.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e0-validation-89baea52-20261018/vienna.sh
qsub -g "$TSUBAME_GROUP" -sync y -tc 30 jobs/e0-validation-89baea52-20261018/linearspace.sh
