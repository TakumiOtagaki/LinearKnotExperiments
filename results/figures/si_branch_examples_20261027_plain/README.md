# SI branch wording revision (2026-10-27)

Figure S1 uses “one ordinary branch” and “one pseudoknot branch” on the lines under A and B. These labels replace the grammar names “V interval” and “WMB interval”. Pair coordinates, colors, line styles, and the 165 × 55 mm canvas are unchanged from `../si_branch_examples_20261027/`. The PNG was opened and checked for placement and legibility.

Reproduce from the experiments repository:

```bash
MPLCONFIGDIR=/tmp/si-branch-fix9-mpl python3 scripts/make_si_branch_examples_20261027.py \
  --out results/figures/si_branch_examples_20261027_plain \
  --context-reference ../../LinearKnot/test/tools/pseudoloop_context_reference.py
```

No RNA calculations were run. The script validates disjoint pair endpoints, noncrossing scaffold and extension components, and weakly closed branch intervals. The external context reference is read only; all purple wing bases remain Multiloop. `examples.json`, `validation.json`, and `manifest.json` retain coordinates, validation results, and hashes.

The renderer also writes the unchanged inside/outside schematic in this directory; only `si_branch_examples.pdf` was exported to the paper repository. The paper's `si_inside_outside.pdf` was preserved.

Before rendering, the shared project quota reported 21,923,932 of 28,000,000 files used. This revision creates eight small files in one new directory, with rendering caches under `/tmp`.
