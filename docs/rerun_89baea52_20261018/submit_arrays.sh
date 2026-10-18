#!/usr/bin/env bash
set -euo pipefail
cd /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments
: "${TSUBAME_GROUP:?Set the authorized accounting group}"
export PYTHONDONTWRITEBYTECODE=1

# Count this bundle together with any already running arrays; total cap is 30.
: "${TMPDIR:?node-local TMPDIR is required}"
budget="${LCP_ARRAY_BUDGET:-30}"
[[ "$budget" =~ ^[0-9]+$ ]] && (( budget >= 1 && budget <= 30 )) || { echo 'LCP_ARRAY_BUDGET must be 1..30' >&2; exit 2; }
max_arrays=3
(( budget >= max_arrays )) || max_arrays="$budget"
tc=$(( budget / max_arrays ))
submission_scratch="$(mktemp -d "$TMPDIR/lcp-submit-89baea52.XXXXXX")"
trap 'rm -rf "$submission_scratch"' EXIT
scripts=(
  jobs/e0-bprna-legacy-89baea52-20261018/nt1_49.sh
  jobs/e0-bprna-legacy-89baea52-20261018/nt50_99.sh
  jobs/e0-bprna-legacy-89baea52-20261018/nt100_199.sh
  jobs/e0-bprna-legacy-89baea52-20261018/nt200_299.sh
  jobs/e0-bprna-legacy-89baea52-20261018/nt300_399.sh
  jobs/e1-bprna-short-89baea52-20261018/nt50_99.sh
  jobs/e1-bprna-short-89baea52-20261018/nt100_199.sh
  jobs/e1-bprna-short-89baea52-20261018/nt200_299.sh
  jobs/e1-bprna-short-89baea52-20261018/nt300_399.sh
  jobs/e1-bprna-long-89baea52-20261018/nt400_799.sh
  jobs/e1-bprna-long-89baea52-20261018/nt800_1599.sh
  jobs/e1-bprna-long-89baea52-20261018/nt1600_3199.sh
  jobs/e1-bprna-long-89baea52-20261018/nt3200_4381.sh
  jobs/e2-rnacentral-89baea52-20261018/rna_4_8kb.sh
  jobs/e2-rnacentral-89baea52-20261018/rna_8_20kb.sh
  jobs/e2-rnacentral-89baea52-20261018/rna_20_30kb.sh
  jobs/e2-synthetic-shapes-89baea52-20261018/synthetic_n250.sh
  jobs/e2-synthetic-shapes-89baea52-20261018/synthetic_n500.sh
  jobs/e2-synthetic-shapes-89baea52-20261018/synthetic_n1000.sh
  jobs/e2-synthetic-shapes-89baea52-20261018/synthetic_n2000.sh
  jobs/e2-synthetic-shapes-89baea52-20261018/synthetic_n4000.sh
  jobs/e2-memory-breakdown-89baea52-20261018/rna_4_8kb.sh
  jobs/e2-memory-breakdown-89baea52-20261018/rna_8_20kb.sh
  jobs/e2-memory-breakdown-89baea52-20261018/rna_20_30kb.sh
  jobs/e2-parallel-impact-w1-89baea52-20261018/parallel_impact.sh
  jobs/e2-parallel-impact-w10-89baea52-20261018/parallel_impact.sh
  jobs/e3-neat1-bpp-89baea52-20261018/neat1_short.sh
  jobs/e3-neat1-bpp-89baea52-20261018/neat1_long.sh
  jobs/e3-neat1-capr-89baea52-20261018/neat1_short.sh
  jobs/e3-neat1-capr-89baea52-20261018/neat1_long.sh
  jobs/e3-neat1-mea-89baea52-20261018/neat1_short.sh
  jobs/e3-neat1-mea-89baea52-20261018/neat1_long.sh
  jobs/e3-neat1-raccess-byloop-89baea52-20261018/neat1_short.sh
  jobs/e3-neat1-raccess-byloop-89baea52-20261018/neat1_long.sh
  jobs/e3-neat1-raccess-89baea52-20261018/neat1_short.sh
  jobs/e3-neat1-raccess-89baea52-20261018/neat1_long.sh
  jobs/e3-neat1-sampling-89baea52-20261018/neat1_short.sh
  jobs/e3-neat1-sampling-89baea52-20261018/neat1_long.sh
  jobs/e3-neat1-variants-bpp-89baea52-20261018/neat1_short.sh
  jobs/e3-neat1-variants-bpp-89baea52-20261018/neat1_long.sh
  jobs/e3-neat1-variants-capr-89baea52-20261018/neat1_short.sh
  jobs/e3-neat1-variants-capr-89baea52-20261018/neat1_long.sh
  jobs/e3-neat1-variants-raccess-89baea52-20261018/neat1_short.sh
  jobs/e3-neat1-variants-raccess-89baea52-20261018/neat1_long.sh
  jobs/pk-recovery-200-turner04-89baea52-20261018/pk200.sh
  jobs/pk-recovery-200-dp09-89baea52-20261018/pk200.sh
  jobs/pk-recovery-80-89baea52-20261018/pk80.sh
  jobs/pk-recovery-rrna-89baea52-20261018/rrna.sh
  jobs/pk-membership-89baea52-20261018/membership.sh
  jobs/e0-validation-89baea52-20261018/enumeration.sh
  jobs/e0-validation-89baea52-20261018/decoding.sh
  jobs/e0-validation-89baea52-20261018/byloop.sh
  jobs/e0-validation-89baea52-20261018/vienna.sh
  jobs/e0-validation-89baea52-20261018/linearspace.sh
)
# qsub -sync y remains active until the allocation's tasks have finished.
# At most max_arrays synchronous submissions run, each with at most tc tasks.
index=0
for script in "${scripts[@]}"; do
  while (( $(jobs -pr | wc -l) >= max_arrays )); do
    wait -n || true
  done
  index=$(( index + 1 ))
  (
    set +e
    qsub -g "$TSUBAME_GROUP" -sync y -tc "$tc" "$script"
    result=$?
    printf '%s\n' "$result" > "$submission_scratch/$index"
    exit 0
  ) &
done
wait
failed=0
for (( i=1; i<=index; i++ )); do
  if [[ ! -f "$submission_scratch/$i" ]] || [[ "$(cat "$submission_scratch/$i")" != 0 ]]; then
    failed=$(( failed + 1 ))
  fi
done
printf 'Completed %s synchronous submissions; %s reported a failure. Aggregate and inspect every planned result.\n' "$index" "$failed"
(( failed == 0 ))
