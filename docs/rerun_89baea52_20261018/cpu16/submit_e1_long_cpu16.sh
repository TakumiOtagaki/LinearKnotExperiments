#!/usr/bin/env bash
# Prepared commands only. Execute explicitly when submission is authorized.
set -euo pipefail
cd /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments
: "${TSUBAME_GROUP:?Set the authorized accounting group}"
concurrency=${LCP_ARRAY_BUDGET:-30}
if ! [[ "$concurrency" =~ ^[0-9]+$ ]] || (( concurrency < 1 || concurrency > 30 )); then
  printf '%s\n' 'LCP_ARRAY_BUDGET must be in 1..30' >&2
  exit 2
fi
for tier in nt400_799 nt800_1599 nt1600_3199 nt3200_4381; do
  qsub -g "$TSUBAME_GROUP" -sync y -tc "$concurrency" \
    "jobs/e1-bprna-long-89baea52-cpu16-20261018/$tier.sh"
done
