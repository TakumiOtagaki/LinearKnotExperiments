# Wording regression: count labels (2026-10-28)

Figure 1: `3,431 / 10,000`, `1,426 / 10,000`, `5,143 / 10,000` become counts `of 10,000`.
Figure 3 (Terc): `Samples / 1,000 draws` becomes `1,000 samples`.
Figure 4 (resources): completion counts use `of`; the LinearKnot line is `LinearKnot completed 107 of 107 (b=100); 80 of 80 (exact)`. The exact denominator is the 80 planned inputs. `C/D` becomes `C and D`.
These changes follow the author's general instruction in docs/feedback_20261026_si/response.md SI44 (paper repo).

元の図・保存データは保持し、数値・対座標・線種・色・軸・フォント・寸法は変えていない。
PNG 全3点を開いて確認した。旧 PNG との画素差は件数ラベルの矩形内だけ。
資源値の `resources/plot_data.json` と完了件数の `resources/completion_counts.tsv` は元図とバイト一致する。
`validation.json` に出典・hash・画素差の範囲を記録した。RNA の再計算はしていない。

## Reproduce from experiments/

Use a fresh output directory for each renderer. Figure 1 uses the unmodified v0.2.0 geometry reference; extract it with `git show v0.2.0:test/tools/pseudoloop_context_reference.py` from the read-only implementation repository into node-local scratch, as documented in figure1_fb6_20261025/README.md.

```sh
python3 -B scripts/make_figure1.py \
  --data results/local/pseudoloop-0.2.0-20261023/examples/figure_data.json \
  --context-reference /tmp/figure1-reference/test/tools/pseudoloop_context_reference.py \
  --out results/figures/wording_regression_20261028/figure1
python3 -B scripts/make_figure4_resources_20261026.py \
  --out /tmp/wording-regression-resources
```

For Terc, use only the renderer and adopted JSON; do not rerun render_pseudoloop_20261023.py, which also rewrites that JSON.

```python
import json, sys
from pathlib import Path
sys.path.insert(0, 'scripts')
import make_figures_feedback_20261020 as renderer
renderer.setup()
data = json.loads(Path('results/local/pseudoloop-0.2.0-20261023/examples/figure_data.json').read_text())
fig, details = renderer.terc_figure(data['terc'])
renderer.save(fig, Path('/tmp/wording-regression-terc'), 'figure4', details)
```

Only the three adopted PDFs were exported to the paper figures directory. Renderer side outputs are not new paper results. No git commits were made.

The resource renderer also checks that this completion line does not overlap any legend text.
