# E0 preparation on d21618c5 and the five-patch oracle

The new plan has 305 inputs, 21 tasks and 915 unstarted jobs. It uses the same
300-RNA cohort and four original witnesses, plus WIP29. The full G dataset and
live dependencies passed verification; `pending_scaffolds` is false.
See the [E0 plan update](../full_bprna_array_plan_20261011.md).

The oracle build succeeded using the unchanged pinned owner's build.sh, with
ordered patches 0001/0002/0006/0007/0008. The new inside-only observers linked
against the verified dedicated d21618c5 exact/RNA libraries and original oracle
objects. No API adaptation was needed. The
[oracle lock](../../manifests/software-cparty-oracle-800e95a-5patches-20261014.json)
and [observer lock](../../manifests/software-e0-d21618c5-20261014.json) freeze hashes.
Historical locks and binaries were retained.

The [five-witness analysis](e0_witness_analysis.json) passes both absolute PF
tolerances on all five inputs. The [observations](e0_witness_observations.json)
also retain the two historical WIP29 observations. Both implementations change
log Z by about +0.02822658740123 and PF free energy by about −0.01739663223684
kcal/mol relative to their old versions. This is a small development check with
Turner04, 37°C and dangles2; it establishes no full-cohort or BPP/MFE result.

The [plan summary](E0/summary.json), [assignments](E0/tasks.tsv), five script
copies and [verification](verification.json) are review artifacts. The runnable
scripts live in `jobs/e0-bprna-legacy-dry-d21618c5-5patches-20261014/`.
Raw witnesses are retained as one ignored tar, with its hash in the observations.
Individual job files were created under `/tmp` and removed after archiving.

The initial build commands were:

```bash
mkdir build/cparty-oracle-800e95a-5patches
git clone -q file://"$PWD/build/cparty-oracle-800e95a-4patches/src" \
  build/cparty-oracle-800e95a-5patches/src
bash deps/LinearCParty-d21618c5/test/oracle/legacy_cparty/build.sh \
  "$PWD/build/cparty-oracle-800e95a-5patches" \
  > build/cparty-oracle-800e95a-5patches/build.log 2>&1
python3 scripts/build_e0_tools.py
```

Those paths now exist. The builder refuses existing output/lock paths; a rebuild
must use new paths and a matching new config. To repeat the small checks and
planning with the retained builds, use fresh output names:

```bash
python3 scripts/check_e0_fixed_witnesses.py --out results/local/e0-witnesses-REPEAT
python3 scripts/plan_array.py \
  --config configs/e0_bprna_legacy_d21618c5_20261014.json \
  --array-config configs/arrays/e0_bprna_legacy_d21618c5_20261014.json \
  --software-root deps/LinearCParty-d21618c5 \
  --build-dir build/e0-observers-d21618c5-800e95a-5patches \
  --out results/runs/e0-plan-REPEAT --dry-run
python3 -m unittest discover -s tests
```

The resource envelope remains provisional. Check short-tier wall time/RSS/caps,
scheduler limits/group/modules/concurrency, fresh shared inode/byte quotas and
node-local scratch capacity before submission. No qsub, git commit or heavy
experiment was run during this preparation; the external development tree was
left unchanged.
