#!/usr/bin/env python3
"""Publish small reviewed summaries and a 165-mm RNAcentral sensitivity figure.

Read analyze_beam_sweep_20261027.py outputs; do not run RNA calculations.
Use fresh --summary-out and --figure-out directories.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(path, rows):
    with path.open('w') as stream:
        writer = csv.DictWriter(stream, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def publish(source, summary_out, figure_out):
    with (source / 'observations.tsv').open() as stream:
        raw = list(csv.DictReader(stream, delimiter='\t'))
    report = json.loads((source / 'summary.json').read_text())
    # This publication is for the complete, validated 7-by-7 sweep.
    assert report['planned'] == report['valid'] == 49
    assert report['status_counts'] == {'success': 49}
    assert report['comparisons_complete'] == 91 and not report['comparisons_unavailable']
    beams = [25, 50, 75, 100, 200, 300, 400]
    names = sorted({r['name'] for r in raw}, key=lambda name: next(int(r['n']) for r in raw if r['name'] == name))
    assert len(names) == 7 and len({(r['name'], int(r['beam'])) for r in raw}) == 49
    for row in raw:
        assert row['status'] == 'success' and row['valid'] == 'True'
        assert float(row['time_wall_s']) < float(row['timeout_s'])
        assert float(row['peak_rss_gib']) * 1024 < float(row['memory_limit_mib'])
        row['F_per_base'] = float(row['F_kcal_mol']) / int(row['n'])
        row['absolute_F_change_per_base'] = (abs(float(row['F_minus_previous_per_nt']))
                                              if row['F_minus_previous_per_nt'] != 'NA' else 'NA')
    summary_out.mkdir(parents=True, exist_ok=False)
    figure_out.mkdir(parents=True, exist_ok=False)
    for name in ('observations.tsv', 'comparisons.tsv', 'summary.json'):
        shutil.copy2(source / name, summary_out / name)
    metrics = ['F_per_base', 'F_minus_previous_per_nt', 'absolute_F_change_per_base',
               'tv_previous_mean', 'tv_previous_p90', 'tv_previous_max', 'time_wall_s', 'peak_rss_gib']
    rollup = []
    for b in beams:
        group = [r for r in raw if int(r['beam']) == b]
        assert len(group) == 7
        row = {'beam': b, 'h': max(b, 30), 'd': max(b, 30), 'completed': 7, 'planned': 7}
        for key in metrics:
            values = [float(r[key]) for r in group if r[key] != 'NA']
            row[key + '_min'] = min(values) if values else 'NA'
            row[key + '_max'] = max(values) if values else 'NA'
        row['time_wall_s_sum'] = sum(float(r['time_wall_s']) for r in group)
        rollup.append(row)
    write_tsv(summary_out / 'by_beam.tsv', rollup)
    write_tsv(summary_out / 'b100.tsv', [r for r in raw if int(r['beam']) == 100])
    adjacent = [{'name': r['name'], 'n': r['n'], 'beam_before': beams[beams.index(int(r['beam'])) - 1],
                 'beam_after': r['beam'], 'F_change_per_base': r['F_minus_previous_per_nt'],
                 'absolute_F_change_per_base': r['absolute_F_change_per_base'],
                 'mean_pairing_total_variation': r['tv_previous_mean'],
                 'maximum_pairing_total_variation': r['tv_previous_max']}
                for r in raw if int(r['beam']) != beams[0]]
    assert len(adjacent) == 42
    write_tsv(summary_out / 'adjacent_widths.tsv', adjacent)

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullLocator, ScalarFormatter
    matplotlib.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7,
                               'axes.labelsize': 7, 'axes.titlesize': 7,
                               'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
                               'legend.fontsize': 6.7, 'pdf.fonttype': 42,
                               'axes.linewidth': .5, 'lines.linewidth': .85})
    fig, axes = plt.subplots(2, 2, figsize=(165 / 25.4, 122 / 25.4))
    fig.subplots_adjust(left=.105, right=.98, bottom=.27, top=.94, hspace=.48, wspace=.36)
    colors = ['#0072B2', '#E69F00', '#009E73', '#CC79A7', '#D55E00', '#56B4E9', '#333333']
    markers = ['o', 's', '^', 'D', 'v', 'P', 'X']
    keys = ['F_per_base', 'tv_previous_mean', 'time_wall_s', 'peak_rss_gib']
    titles = ['Ensemble free energy per base', 'Change in pairing distribution',
              'Calculation time', 'Peak memory']
    labels = ['Ensemble free energy\n(kcal/mol per base)', 'Mean per-base\ntotal variation',
              'Wall-clock time (seconds)', 'Peak memory (GiB)']
    for name, color, marker in zip(names, colors, markers):
        group = sorted((r for r in raw if r['name'] == name), key=lambda r: int(r['beam']))
        assert [int(r['beam']) for r in group] == beams
        label = f'{name} ({int(group[0]["n"]):,} bases)'
        for ax, key in zip(axes.flat, keys):
            points = [r for r in group if r[key] != 'NA']
            x = [int(r['beam']) for r in points]
            y = [float(r[key]) for r in points]
            ax.plot(x, y, marker=marker, markersize=3, markeredgewidth=.55,
                    color=color, label=label, markerfacecolor='white')
            ax.plot([100], [y[x.index(100)]], marker=marker, markersize=4.4, color=color,
                    markeredgewidth=.6)
    for letter, ax, title, label in zip('ABCD', axes.flat, titles, labels):
        ax.set_title(title, loc='left', pad=6)
        ax.text(-.19, 1.07, letter, transform=ax.transAxes, fontweight='bold', fontsize=8)
        ax.axvline(100, color='.45', linewidth=.7, linestyle=':', zorder=0)
        ax.set(xlabel='Beam width b', ylabel=label, xlim=(15, 410), xticks=beams)
        if ax is axes[0, 1]:
            ax.set(xlabel='Beam width after increase', xticks=beams[1:])
        ax.tick_params(length=2.5, width=.5, pad=2)
        ax.grid(axis='y', alpha=.18, linewidth=.5)
        ax.spines[['top', 'right']].set_visible(False)
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)
            tick.set_ha('right')
            tick.set_rotation_mode('anchor')
            if tick.get_text() == '100':
                tick.set_fontweight('bold')
    axes[0, 1].set_ylim(bottom=0)
    for ax in axes[1]:
        ax.set_yscale('log')
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.yaxis.set_minor_locator(NullLocator())
    handles, legend_labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, legend_labels, loc='lower center', bbox_to_anchor=(.52, .006),
               ncol=2, frameon=False, columnspacing=1.6, handlelength=2, labelspacing=.65)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text()]
    for t in texts:
        assert t.get_fontsize() >= 6
        box = t.get_window_extent(renderer)
        assert box.x0 >= 0 and box.y0 >= 0 and box.x1 <= fig.bbox.width and box.y1 <= fig.bbox.height, t.get_text()
    # Keep the physical size fixed: no tight bounding-box rescaling.
    stem = 'rnacentral_beam_sweep_si'
    for extension in ('pdf', 'png'):
        fig.savefig(figure_out / f'{stem}.{extension}', dpi=300)
    plt.close(fig)
    manifest = {
        'source': str(source), 'source_sha256': {f: digest(source / f) for f in ('observations.tsv', 'comparisons.tsv', 'summary.json')},
        'script_sha256': digest(Path(__file__)), 'pdf_dimensions_mm': [165, 122],
        'minimum_font_pt': 6.5, 'body_font_pt': 7, 'series': names, 'beams': beams,
        'points_per_panel': {'A': 49, 'B': 42, 'C': 49, 'D': 49},
        'panels': keys, 'highlighted_beam': 100,
        'adjacent_beam_pairs': list(zip(beams[:-1], beams[1:])),
        'interpretation': 'Beam width and local limits vary together: h=d=max(b,30). A shows absolute ensemble free energy per base; B compares adjacent widths over scaffold-free bases, including the unpaired state. No width is a reference for the sweep.',
        'timing': 'One single-thread process per task; one observation, no warm-up. Whole fold includes partition function, minimum free energy, outside and zero-cutoff pairing-probability output; excludes build, scaffold construction and compression.',
        'files': {p.name: digest(p) for p in figure_out.iterdir() if p.is_file()}}
    (figure_out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return rollup


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, type=Path)
    parser.add_argument('--summary-out', required=True, type=Path)
    parser.add_argument('--figure-out', required=True, type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='beam-sweep-plot-', dir='/tmp') as cache:
        os.environ['MPLCONFIGDIR'] = cache
        publish(args.source, args.summary_out, args.figure_out)
    print(args.figure_out / 'rnacentral_beam_sweep_si.pdf')


if __name__ == '__main__':
    main()
