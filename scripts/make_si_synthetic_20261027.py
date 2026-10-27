#!/usr/bin/env python3
"""Render Figure S4's time/memory panels from the adopted saved plot data.

Run from experiments/:
  python3 -B scripts/make_si_synthetic_20261027.py --out results/figures/NEW
No RNA calculations; existing six-panel outputs are left intact.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
from types import SimpleNamespace

import make_figures as base
import make_figures_review as review
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter, NullLocator
import matplotlib.pyplot as plt
import numpy as np


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=base.ROOT / 'results/figures/main_4f97fd3_20261026_si_labels/plot_data.tsv')
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    with args.source.open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    summaries = [r for r in rows if r['panel'] == '3C']
    selected = [r for r in summaries if r['program'] == 'LinearCParty']
    assert len(selected) == 70
    data = SimpleNamespace(synthetic_summary={
        (r['program'], r['series'], int(r['length']), r['metric']):
        {k: float(r[k]) for k in ('median', 'min', 'max')} for r in selected})
    saved = [json.loads(line) for line in (base.ROOT / 'data/e2-synthetic-shapes-20261012/inputs.jsonl').read_text().splitlines()]
    designs = {r['family']: r['scaffold'] for r in saved if len(r['scaffold']) == 250 and r['family'] in base.SERIES}
    assert len(designs) == 7
    review.style(base)
    fig = plt.figure(figsize=(165 * base.MM, 90 * base.MM))
    for col, metric in enumerate(('wall_seconds', 'rss_gib')):
        ax = base.mm_axes(fig, (19, 102)[col], 47, 60, 36)
        base.plot_synthetic(ax, data, 'LinearCParty', metric)
        # Preserve the original A/B limits as well as every point and error bar.
        all_values = [r for r in summaries if r['metric'] == metric]
        lo = min(float(r['min']) for r in all_values) * .65
        hi = max(float(r['max']) for r in all_values) * 1.4
        ax.set_ylim(lo, hi)
        ticks = [10.**k for k in range(math.floor(math.log10(lo)), math.ceil(math.log10(hi))+1)]
        ax.yaxis.set_major_locator(FixedLocator(ticks))
        ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: np.format_float_positional(x, trim='-')))
        ax.yaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.xaxis.set_minor_locator(NullLocator())
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.set_ylabel('Wall time (s)' if col == 0 else 'Peak memory (GiB)')
        review.panel(fig, 'AB'[col], (1, 84)[col], 85)
    for i, series in enumerate(base.SERIES):
        col, row = i % 4, i // 4
        left, bottom = 3 + col * 41, 17 - row * 14
        review.draw_design(base.mm_axes(fig, left, bottom, 17, 8), designs[series])
        label = review.SHORT_SERIES[i].replace('Empty G', 'Empty\nscaffold').replace('Free bases together', 'Free bases\ntogether').replace('Balanced multiloops', 'Balanced\nmultiloops').replace('Nested multiloops', 'Nested\nmultiloops')
        review.text_mm(fig, left + 20, bottom + 6, label, va='center', color='.15')
        key = base.mm_axes(fig, left + 20, bottom - .5, 15, 2)
        key.plot([0, 1], [.5, .5], color=base.COLORS[i], ls=base.LINES[i], marker=base.MARKERS[i], markevery=[0], markersize=3, markerfacecolor='white')
        key.set_axis_off()
    records = []
    base.save(fig, 'figure3_synthetic_all_cli_si', args.out, 165, records, check_text_overlap=True)
    plt.close(fig)
    fields = list(selected[0])
    with (args.out / 'plot_data.tsv').open('w') as stream:
        writer = csv.DictWriter(stream, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(selected)
    assert list(csv.DictReader((args.out / 'plot_data.tsv').open(), delimiter='\t')) == selected
    manifest = {'source': str(args.source.relative_to(base.ROOT)), 'source_sha256': digest(args.source),
                'script_sha256': digest(Path(__file__)),
                'helper_sha256': {name: digest(base.ROOT / 'scripts' / name) for name in ('make_figures.py', 'make_figures_review.py')}, 'source_base_commit': '955807764f058e356af8614d655e07ff74bb9f0e',
                'data_rows': len(selected), 'series': 7, 'lengths': [250,500,1000,2000,4000],
                'panels': ['wall_seconds','rss_gib'], 'pdf_dimensions_mm': [165,90],
                'preserved': 'All 70 A/B point, min and max records, axis scales and limits; original series styles and scaffold sketches.',
                'files': {p.name:digest(p) for p in args.out.iterdir() if p.is_file()}}
    (args.out / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    (args.out / 'README.md').write_text('''# Figure S4: two panels (2026-10-27)

Time and peak memory for partition function, BPPs and MFE. All 70 saved
A/B records and their nine-measurement ranges are preserved, along with
axis scales, limits, series styles and scaffold sketches. The other four
panels remain in the earlier six-panel directory. No RNA calculation.

Reproduce from experiments/ into a fresh directory:

```sh
python3 -B scripts/make_si_synthetic_20261027.py --out results/figures/NEW
```

The figure is 165 x 90 mm. Source and output hashes are in manifest.json.
''')
    print(json.dumps({'output':str(args.out),'rows':len(selected),'files':len(list(args.out.iterdir()))}))


if __name__ == '__main__':
    main()
