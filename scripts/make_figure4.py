#!/usr/bin/env python3
"""Make the 183 x 90 mm NEAT1 manuscript figure from saved E3 outputs.

Run from either repository root. No RNA CLI, tar extraction, or new experiment
is used. Archived members are streamed, checked against recorded hashes,
and compared with the existing analysis tables where available.
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys
import tarfile
import tempfile

import numpy as np

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'figure4-mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import MaxNLocator, StrMethodFormatter

from analyze_e3_neat1 import crossing, pairs_of, parse_access, parse_bpp, parse_context, structure_from_stdout

ROOT = Path(__file__).resolve().parents[1]
MM = 1 / 25.4
WIDTH, HEIGHT = 183, 90
DEFAULT_REGION = (21000, 21300)
DEFAULT_FOCUS = (21134, 21143)
WINDOW = 8
PAIR_CUTOFF = .01
INK, MUTED, GRID = '#1f1f1e', '#6b6b66', '#d9d8d3'
G_COLOR, NEW_COLOR, CROSS_COLOR = '#8a8a85', '#2a78d6', '#eb6834'
FOCUS_COLOR = '#f4efcf'
CONTEXTS = ('Stem', 'Hairpin', 'Bulge', 'Internal', 'Multiloop', 'Exterior', 'PseudoLoop')
CONTEXT_COLORS = ('#2a78d6', '#1baf7a', '#eda100', '#e87ba4', '#6250d6', '#b9b8b1', '#eb6834')


def sha(data):
    return hashlib.sha256(data).hexdigest()


def rows(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as f:
        yield from csv.DictReader(f, delimiter='\t')


def supplemental_bpp(provenance, beam):
    """Read only the named empty-G BPP job omitted by the earlier snapshot."""
    directory = f'jobs/NEAT1_2_emptyG-beam{beam}_outer-LinearCParty-rep0'
    if any(a['member'] == directory + '/bpp.tsv' for a in provenance['artifacts']):
        return [], []
    task = ROOT / 'results/runs/e3-neat1-variants-bpp-e2cecbb4-20261016/tasks/neat1_long/000001/run'
    checkpoint = task / 'task_results.jsonl'
    blob = checkpoint.read_bytes()
    complete = b''.join(line for line in blob.splitlines(keepends=True) if line.endswith(b'\n'))
    matches = [line for line in complete.splitlines() if line.strip()
               and json.loads(line)['id'] == Path(directory).name]
    assert len(matches) == 1, 'Required empty-G BPP checkpoint is missing or duplicated'
    record = json.loads(matches[0])
    assert record['status'] == 'success' and record['returncode'] == 0 and not record['warmup']
    assert record['input'] == 'NEAT1_2_emptyG' and record['directory'] == directory
    assert record['output_archive'] == 'task_outputs.tar'
    archive = str((task / record['output_archive']).relative_to(ROOT))
    extra = []
    for basename in ('stdout.txt', 'bpp.tsv'):
        member = directory + '/' + basename
        if any(a['archive'] == archive and a['member'] == member for a in provenance['artifacts']):
            continue
        extra.append(dict(archive=archive, member=member, **record['artifacts'][member],
                          command=record['command'], finished_utc=record['finished_utc'],
                          archive_end=record['archive_end'], source='supplemental_checkpoint'))
    return extra, [dict(path=str(checkpoint.relative_to(ROOT)), snapshot_bytes=len(blob),
                        complete_bytes=len(complete), sha256_complete_lines=sha(complete),
                        ignored_partial_tail=len(blob) - len(complete), selected_job=record,
                        selected_record_sha256=sha(matches[0]))]


def read_data(data, beam, region, focus_starts):
    """Read the E3 snapshot and the explicitly named supplemental BPP job."""
    provenance = json.loads((data / 'provenance.json').read_text())
    table_names = ('provenance.json', 'position_profiles.tsv.gz', 'accessibility.tsv.gz',
                   'bpp_pairs.tsv.gz', 'sampling_histograms.tsv', 'sampling_summary.tsv',
                   'sampling_reference_distances.tsv.gz', 'representative_distances.tsv')
    sources = [{'path': str((data / name).relative_to(ROOT)),
                'sha256': sha((data / name).read_bytes())} for name in table_names]
    inputs = {}
    for name, path in [('NEAT1_2', 'data/neat1-20261016/NEAT1_2.txt'),
                       ('NEAT1_2_emptyG', 'data/neat1-variants-20261016/NEAT1_2_emptyG.txt'),
                       ('NEAT1_1', 'data/neat1-20261016/NEAT1_1.txt')]:
        raw = (ROOT / path).read_bytes()
        expected = next(r for r in provenance['inputs'] if r['path'] == path)
        assert sha(raw) == expected['sha256']
        header, seq, db = [s.strip() for s in raw.decode().splitlines() if s.strip()]
        assert header == '>' + name and len(seq) == len(db)
        inputs[name] = dict(name=name, sequence=seq, scaffold=db, n=len(seq), fixed=pairs_of(db))
        sources.append(expected)
    assert inputs['NEAT1_2']['n'] == 22743 and inputs['NEAT1_1']['n'] == 3756
    assert len(inputs['NEAT1_2']['fixed']) == 4161
    assert len(inputs['NEAT1_1']['fixed']) == 727
    assert not inputs['NEAT1_2_emptyG']['fixed']
    assert inputs['NEAT1_2']['sequence'] == inputs['NEAT1_2_emptyG']['sequence']

    # The earlier snapshot predates empty-G BPP at b=200. Read that one saved
    # checkpoint if needed; do not scan other runs or regenerate the E3 tables.
    extra, checkpoints = supplemental_bpp(provenance, beam)
    access_beams = {'NEAT1_2': {beam, 200, 400}, 'NEAT1_2_emptyG': {beam, 100, 200}}
    # Only the saved artifacts needed by panels B--D are consumed.
    wanted = {}
    for a in provenance['artifacts'] + extra:
        name = Path(a['member']).parent.name
        selected = any(name.startswith(f'{v}-beam{beam}_outer-') for v in ['NEAT1_2', 'NEAT1_2_emptyG'])
        kind = a['archive'].split('/')[2].replace('e3-neat1-', '').replace('variants-', '').split('-e2')[0]
        selected = selected and kind in ('bpp', 'capr', 'raccess')
        selected |= kind == 'raccess' and any(
            name.startswith(f'{v}-beam{b}_outer-') for v, beams in access_beams.items() for b in beams)
        selected |= name.startswith(f'NEAT1_1-beam{beam}_outer-') and kind in ('bpp', 'sampling', 'mea') and a['member'].endswith('stdout.txt')
        if selected:
            wanted.setdefault(a['archive'], {})[a['member']] = (kind, a)
    saved, consumed = {}, []
    for archive, members in sorted(wanted.items()):
        with tarfile.open(ROOT / archive, 'r|') as stream:
            for member in stream:
                if member.name not in members:
                    continue
                kind, expected = members.pop(member.name)
                if 'archive_end' in expected:
                    assert member.offset_data + member.size <= expected['archive_end']
                raw = stream.extractfile(member).read()
                assert len(raw) == expected['bytes'] and sha(raw) == expected['sha256']
                cmd = expected['command']
                saved_beam = int(cmd[cmd.index('-b') + 1])
                assert f'-beam{saved_beam}_outer-' in member.name
                assert saved_beam == beam or (kind == 'raccess' and saved_beam in access_beams[
                    Path(member.name).parent.name.split('-beam')[0]])
                assert cmd[cmd.index('--ranking') + 1] == 'outer'
                for flag in ('--max-hairpin', '--max-case3-depth'):
                    assert int(cmd[cmd.index(flag) + 1]) == max(saved_beam, 30)
                if kind == 'bpp':
                    assert float(cmd[cmd.index('--bpp-cutoff') + 1]) == .001
                name = Path(member.name).parent.name.split('-beam')[0]
                saved[name, saved_beam, kind, Path(member.name).name] = raw.decode()
                consumed.append(expected)
                if not members:
                    break
        assert not members

    paired = np.array([c != '.' for c in inputs['NEAT1_2']['scaffold']])
    fixed_window = np.convolve(paired.astype(int), np.ones(WINDOW, int), 'valid') > 0
    focal = slice(focus_starts[0] - 1, focus_starts[1])
    focal_bases = slice(focus_starts[0] - 1, focus_starts[1] + WINDOW - 1)
    assert not fixed_window[focal].any(), 'Focus windows contain scaffold-paired bases'
    access_by_beam, access_checks = {}, {}
    for name, beams in access_beams.items():
        method = 'G0.7' if name == 'NEAT1_2' else 'emptyG'
        access_by_beam[method], access_checks[method] = {}, {}
        for b in sorted(beams):
            acc = parse_access(saved[name, b, 'raccess', 'stdout.txt'], inputs[name])[WINDOW]
            residual = float(np.max(acc[fixed_window])) if method == 'G0.7' else 0.
            assert residual < 1e-8
            if method == 'G0.7':
                acc[fixed_window] = 0
            table = [r for r in rows(data / 'accessibility.tsv.gz') if r['iso'] == 'NEAT1_2'
                     and r['method'] == method and r['beam'] == str(b) and r['width'] == str(WINDOW)]
            assert [int(r['start']) for r in table] == list(range(1, inputs[name]['n'] - WINDOW + 2))
            assert np.allclose(acc, [float(r['accessibility']) for r in table], atol=1e-12, rtol=1e-12)
            access_by_beam[method][b] = acc
            access_checks[method][b] = dict(focus_mean=float(acc[focal].mean()),
                                            focus_values=acc[focal].tolist(),
                                            fixed_window_raw_residual=residual,
                                            full_transcript_table_crosschecked=True)
    stability = {}
    for method, widths in [('G0.7', (200, 400)), ('emptyG', (100, 200))]:
        a, b = (access_by_beam[method][w][focal] for w in widths)
        change = float(np.max(np.abs(a - b)))
        stability[method] = dict(beam_widths=widths, focus_windows=len(a),
                                  mean_absolute_change=float(np.abs(a - b).mean()),
                                  absolute_mean_change=float(abs(a.mean() - b.mean())),
                                  max_absolute_window_change=change,
                                  every_window_change_below_0_05=change < .05)
    profile, access, pairs, diagnostics = {}, {}, {}, {}
    for name, method in [('NEAT1_2', 'G0.7'), ('NEAT1_2_emptyG', 'emptyG')]:
        v = inputs[name]
        ctx = parse_context(saved[name, beam, 'capr', 'stdout.txt'], v)
        acc = access_by_beam[method][beam]
        own_paired = np.array([c != '.' for c in v['scaffold']])
        raw_fixed_residual = access_checks[method][beam]['fixed_window_raw_residual']
        for c in CONTEXTS:
            assert np.allclose(ctx[c][own_paired], 1 if c == 'Stem' else 0, atol=1e-8)
            ctx[c][own_paired] = 1 if c == 'Stem' else 0
        assert np.max(np.abs(sum(ctx.values()) - 1)) < 1e-8
        # Compare full profiles to the independently saved table coordinates.
        table = [r for r in rows(data / 'position_profiles.tsv.gz') if r['iso'] == 'NEAT1_2' and r['method'] == method and r['beam'] == str(beam)]
        assert [int(r['position']) for r in table] == list(range(1, v['n'] + 1))
        for c in CONTEXTS:
            assert np.allclose(ctx[c], [float(r[c]) for r in table], atol=1e-12, rtol=1e-12)
        fixed = inputs['NEAT1_2']['fixed'] if method == 'G0.7' else set()
        bp, _ = parse_bpp(saved[name, beam, 'bpp', 'bpp.tsv'], dict(v, mask=own_paired))
        dbs = structure_from_stdout(saved[name, beam, 'bpp', 'stdout.txt'], v['n'])
        assert len(dbs) == 2 and dbs[0] == v['scaffold']
        bptable = [r for r in rows(data / 'bpp_pairs.tsv.gz') if r['iso'] == 'NEAT1_2' and r['method'] == method and r['beam'] == str(beam)]
        table_available = bool(bptable)
        if table_available:
            assert len(bp) == len(bptable)
            for r in bptable:
                pair = int(r['i']), int(r['j'])
                assert bp[pair] == float(r['probability'])
                assert int(r['is_scaffold']) == int(pair in fixed)
                assert int(r['crosses_G']) == int(crossing(pair, fixed))
        else:
            assert method == 'emptyG' and beam == 200 and extra
            bptable = [dict(i=i, j=j, probability=p, is_scaffold='0', crosses_G='0')
                       for (i, j), p in sorted(bp.items())]
        incident = np.zeros(v['n'])
        for (i, j), p in bp.items():
            incident[[i - 1, j - 1]] += p
            if (i, j) not in fixed:
                assert not own_paired[i - 1] and not own_paired[j - 1]
        assert incident.max() <= 1 + 1e-8
        assert np.min(ctx['Stem'] - incident) > -1e-8
        profile[method], access[method], pairs[method] = ctx, acc, bptable
        diagnostics[method] = dict(focus_mean=float(acc[focal].mean()),
                                   fixed_window_raw_residual=raw_fixed_residual,
                                   bpp_table_crosschecked=table_available,
                                   bpp_pairs=len(bp),
                                   focus_windows=focus_starts[1] - focus_starts[0] + 1,
                                   focus_scaffold_paired_bases=int(own_paired[focal_bases].sum()),
                                   focus_context_max={c: float(ctx[c][focal_bases].max()) for c in CONTEXTS})

    # Recompute every pair-set distance from streamed samples with integer bitsets.
    # This checks the plotted histograms, references, and exact means independently.
    hist = [r for r in rows(data / 'sampling_histograms.tsv') if r['beam'] == str(beam)]
    summary = [r for r in rows(data / 'sampling_summary.tsv') if r['beam'] == str(beam)]
    refdist = list(rows(data / 'sampling_reference_distances.tsv.gz'))
    representative = list(rows(data / 'representative_distances.tsv'))
    sample_checks = []
    fixed = inputs['NEAT1_1']['fixed']
    sample_cmd = next(a['command'] for a in consumed if f'NEAT1_1-beam{beam}_' in a['member'] and '-sampling-' in a['archive'])
    assert sample_cmd[sample_cmd.index('--seed') + 1] == '20261016'
    seqs = structure_from_stdout(saved['NEAT1_1', beam, 'sampling', 'stdout.txt'], 3756, sampled=True)
    assert len(seqs) == len(set(seqs)) == 1000
    ss = [pairs_of(db) for db in seqs]
    assert all(fixed <= s for s in ss)
    refs = {}
    for ref, kind in [('MFE', 'bpp'), ('MEA', 'mea')]:
        if kind == 'mea':
            cmd = next(a['command'] for a in consumed if f'NEAT1_1-beam{beam}_' in a['member'] and '-mea-' in a['archive'])
            assert float(cmd[cmd.index('--gamma') + 1]) == 1.
        dbs = structure_from_stdout(saved['NEAT1_1', beam, kind, 'stdout.txt'], 3756)
        assert len(dbs) == 2 and dbs[0] == inputs['NEAT1_1']['scaffold']
        refs[ref] = pairs_of(dbs[-1])
        assert fixed <= refs[ref]
    sets = [s - fixed for s in ss] + [refs['MFE'] - fixed, refs['MEA'] - fixed]
    universe = {p: k for k, p in enumerate(sorted(set().union(*sets)))}
    masks = [sum(1 << universe[p] for p in s) for s in sets]
    popcount = getattr(int, 'bit_count', lambda x: bin(x).count('1'))
    def distance(a, b):
        u = popcount(a | b)
        return 1 - popcount(a & b) / u if u else 0.
    distances = {'sample_sample': np.fromiter((distance(masks[i], masks[j]) for i in range(1000) for j in range(i + 1, 1000)), float, count=499500)}
    for ref, m in zip(('MFE', 'MEA'), masks[-2:]):
        distances['sample_' + ref] = np.array([distance(s, m) for s in masks[:1000]])
        expected = [float(r['jaccard_distance']) for r in refdist if r['beam'] == str(beam) and r['space'] == 'extension' and r['reference'] == ref]
        assert np.allclose(distances['sample_' + ref], expected, atol=1e-14, rtol=0)
    for comparison, values in distances.items():
        assert values.min() >= .3 and values.max() <= .95
        h = [r for r in hist if r['beam'] == str(beam) and r['space'] == 'extension' and r['comparison'] == comparison]
        assert len(h) == 100
        counts, _ = np.histogram(values, bins=np.linspace(0, 1, 101))
        assert np.array_equal(counts, [int(r['count']) for r in h])
        s = next(r for r in summary if r['beam'] == str(beam) and r['space'] == 'extension' and r['comparison'] == comparison)
        assert len(values) == int(s['n']) and abs(values.mean() - float(s['mean'])) < 1e-12
    r = next(r for r in representative if r['beam'] == str(beam) and r['space'] == 'extension')
    assert abs(distance(*masks[-2:]) - float(r['jaccard_distance'])) < 1e-14
    sample_checks.append(dict(beam=beam, samples=1000, sample_pairs=499500, histogram_counts_exact=True,
                              means_verified=True, reference_distances_verified=True,
                              exact_means={k: float(v.mean()) for k, v in distances.items()}))
    return dict(beam=beam, region=region, focus=focus_starts, inputs=inputs,
                profile=profile, access=access, access_by_beam=access_by_beam,
                access_checks=access_checks, stability=stability, pairs=pairs, fixed_window=fixed_window,
                hist=hist, summary=summary, diagnostics=diagnostics, sample_checks=sample_checks,
                sources=sources, artifacts=consumed, supplemental_checkpoints=checkpoints)


def axis_mm(fig, x, y, w, h):
    return fig.add_axes([x / WIDTH, y / HEIGHT, w / WIDTH, h / HEIGHT])


def text_mm(fig, x, y, s, **kwargs):
    return fig.text(x / WIDTH, y / HEIGHT, s, **kwargs)


def panel(fig, x, y, letter, title):
    text_mm(fig, x, y, letter, fontsize=9, weight='bold', va='top')
    text_mm(fig, x + 5, y - .1, title, fontsize=7, weight='bold', va='top')


def clean_axis(ax):
    ax.spines[['top', 'right']].set_visible(False)
    ax.tick_params(pad=1.5)


def focus(ax, window_starts, starts=False):
    ax.axvspan(window_starts[0] - .5, window_starts[1] + (.5 if starts else WINDOW - .5),
               facecolor=FOCUS_COLOR, lw=0, zorder=-1)
    for x in (window_starts[0] - .5, window_starts[1] + (.5 if starts else WINDOW - .5)):
        ax.axvline(x, color='#c3b879', lw=.4, zorder=0)


def draw(fig, d):
    region, focal = d['region'], d['focus']
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 6,
                         'text.color': INK, 'axes.labelcolor': INK,
                         'axes.edgecolor': MUTED, 'axes.linewidth': .5,
                         'xtick.color': MUTED, 'ytick.color': MUTED,
                         'xtick.major.size': 2, 'ytick.major.size': 2,
                         'xtick.major.width': .5, 'ytick.major.width': .5,
                         'pdf.fonttype': 42, 'hatch.linewidth': .35})
    panel(fig, 1.5, 88.8, 'A', 'NEAT1_2  ·  22,743 nt')
    ax = axis_mm(fig, 13, 77.7, 109, 6)
    n = d['inputs']['NEAT1_2']['n']
    ax.plot([1, n], [0, 0], color=G_COLOR, lw=2, solid_capstyle='butt')
    ax.plot(region, [0, 0], color=INK, lw=3, solid_capstyle='butt')
    ax.plot([np.mean(focal)], [0], marker='v', color=INK, ms=3)
    ax.set(xlim=(1, n), ylim=(-1, 1), xticks=[1, 10000, n],
           xticklabels=['1', '10,000', '22,743'], yticks=[])
    ax.spines[:].set_visible(False)
    ax.tick_params(axis='x', length=0, pad=-1)
    ax.annotate(f'Shown below: {region[0]:,}–{region[1]:,}', xy=(np.mean(region), .15), xycoords='data',
                xytext=(.64, .98), textcoords='axes fraction', ha='center', va='bottom',
                fontsize=5.9, arrowprops=dict(arrowstyle='-', color=MUTED, lw=.5))
    text_mm(fig, 13, 75.5, f'Highlighted: {focal[1] - focal[0] + 1} window starts, {focal[0]:,}–{focal[1]:,}', fontsize=6, va='top',
            bbox=dict(facecolor=FOCUS_COLOR, edgecolor='none', pad=1.2))

    panel(fig, 1.5, 70.3, 'B', f'Base-pairing probabilities  ·  beam width {d["beam"]}')
    # B and C have identical coordinate limits and aligned comparison columns.
    columns = [('G0.7', 13, 'Scaffold G: threshold 0.7'), ('emptyG', 72, 'Empty scaffold')]
    arc_counts = {}
    arc_height = max(((int(r['j']) - int(r['i'])) / 2 for method, _, _ in columns
                     for r in d['pairs'][method] if float(r['probability']) >= PAIR_CUTOFF
                     and region[0] <= int(r['i']) < int(r['j']) <= region[1]), default=1) * 1.05
    for method, x, label in columns:
        text_mm(fig, x, 65.6, label, fontsize=6.5, weight='bold', va='top')
        ax = axis_mm(fig, x, 48.3, 50, 13.5)
        focus(ax, focal)
        paths, colors = [], []
        counts = {'scaffold': 0, 'new_nested': 0, 'crossing': 0, 'external_omitted': 0}
        selected = sorted(d['pairs'][method], key=lambda r: float(r['probability']))
        for r in selected:
            i, j, p = int(r['i']), int(r['j']), float(r['probability'])
            if p < PAIR_CUTOFF:
                continue
            if (region[0] <= i <= region[1]) != (region[0] <= j <= region[1]):
                counts['external_omitted'] += 1
            if not region[0] <= i < j <= region[1]:
                continue
            category = 'scaffold' if r['is_scaffold'] == '1' else ('crossing' if r['crosses_G'] == '1' else 'new_nested')
            color = {'scaffold': G_COLOR, 'new_nested': NEW_COLOR, 'crossing': CROSS_COLOR}[category]
            t = np.linspace(0, np.pi, 70)
            y = (j - i) / 2 * np.sin(t)
            paths.append(np.column_stack(((i + j) / 2 + (j - i) / 2 * np.cos(t), y)))
            colors.append(matplotlib.colors.to_rgba(color, p))
            counts[category] += 1
        ax.add_collection(LineCollection(paths, colors=colors, linewidths=.65))
        ax.axhline(0, color=MUTED, lw=.4)
        ax.set(xlim=(region[0] - .5, region[1] + .5), ylim=(-.02 * arc_height, arc_height), xticks=[], yticks=[])
        ax.spines[:].set_visible(False)
        arc_counts[method] = counts
    handles = [Line2D([], [], color=G_COLOR, lw=1, label='Pair in G'),
               Line2D([], [], color=NEW_COLOR, lw=1, label='New pair')]
    fig.legend(handles=handles, loc='center left', bbox_to_anchor=(13 / WIDTH, 45.8 / HEIGHT),
               frameon=False, ncol=2, fontsize=5.7, borderaxespad=0, handlelength=1.5,
               columnspacing=1.2, handletextpad=.45)
    text_mm(fig, 122, 45.8, '≥0.01; darker = higher', fontsize=5.7, ha='right', va='center')

    panel(fig, 1.5, 42.8, 'C', 'Seven structural-context probabilities (LinearPKCapR)')
    fig.legend(handles=[Line2D([], [], color=INK, lw=.8, label=f'Accessibility: b={d["beam"]}'),
                        Line2D([], [], color=NEW_COLOR, lw=1.4, alpha=.4, label='G: b=400')],
               loc='center left', bbox_to_anchor=(13 / WIDTH, 39 / HEIGHT), ncol=2,
               fontsize=5.1, frameon=False, borderaxespad=0, handlelength=1.8,
               handletextpad=.4, columnspacing=1.2)
    for method, x, label in columns:
        ax = axis_mm(fig, x, 22.8, 50, 10.8)
        pos = np.arange(region[0], region[1] + 1)
        ctx = np.array([d['profile'][method][c][region[0] - 1:region[1]] for c in CONTEXTS])
        # Per-base steps, no smoothing or interpolation between base probabilities.
        ax.stackplot(pos, ctx, colors=CONTEXT_COLORS, linewidth=0, step='mid')
        for xx in (focal[0] - .5, focal[1] + WINDOW - .5):
            ax.axvline(xx, color=INK, lw=.45, ls=(0, (2, 2)))
        ax.set(xlim=(region[0] - .5, region[1] + .5), ylim=(0, 1), yticks=[0, 1], xticks=[])
        clean_axis(ax)
        ax = axis_mm(fig, x, 7.7, 50, 11.7)
        focus(ax, focal, starts=True)
        acc = d['access'][method][region[0] - 1:region[1]]
        if method == 'G0.7':
            ax.plot(pos, d['access_by_beam'][method][400][region[0] - 1:region[1]],
                    color=NEW_COLOR, lw=1.4, alpha=.4, drawstyle='steps-mid', zorder=2.5)
        ax.plot(pos, acc, color=INK, lw=.8, drawstyle='steps-mid', zorder=3)
        if method == 'G0.7':
            mask = d['fixed_window'][region[0] - 1:region[1]]
            for a, b in boolean_runs(mask):
                ax.axvspan(region[0] + a - .5, region[0] + b - .5,
                           facecolor='#e8e8e4', edgecolor=G_COLOR, hatch='////', lw=0, zorder=1)
        # Mark the focal segment on the zero baseline even when the curve is tiny.
        ax.plot([focal[0] - .5, focal[1] + .5], [0, 0], color='#9e8b36', lw=1.1, zorder=2)
        mean = d['diagnostics'][method]['focus_mean']
        text_mm(fig, x, 21.1, f'8-nt accessibility (LinearPKRaccess) · mean {mean:.4f}',
                fontsize=5.4, va='top')
        ax.set(xlim=(region[0] - .5, region[1] + .5), ylim=(-.015, 1.05), yticks=[0, 1])
        ticks = MaxNLocator(nbins=3, integer=True).tick_values(*region)
        ax.set_xticks(ticks[(ticks >= region[0]) & (ticks <= region[1])])
        ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
        clean_axis(ax)
        ax.tick_params(labelsize=5.5)
    text_mm(fig, 67.5, 1.3, 'Nucleotide position (B, contexts); window start (accessibility)', fontsize=5.7, ha='center', va='bottom')

    # One compact seven-context key is shared by both comparison columns.
    # The title band above C carries it, so no legend overlaps the data.
    legend = [Patch(facecolor=c, edgecolor='none', label=n) for n, c in zip(CONTEXTS, CONTEXT_COLORS)]
    legend.append(Patch(facecolor='#e8e8e4', edgecolor=G_COLOR, hatch='////', label='Fixed-base window'))
    fig.legend(handles=legend, loc='lower left', bbox_to_anchor=(13 / WIDTH, 35.1 / HEIGHT),
               ncol=8, fontsize=5.1, frameon=False, borderaxespad=0, handlelength=.8,
               handletextpad=.3, columnspacing=.65)

    # Panel D is a separate isoform; it is explicitly labelled, not coordinate-linked.
    fig.add_artist(Line2D([130 / WIDTH] * 2, [3 / HEIGHT, 89 / HEIGHT],
                         transform=fig.transFigure, color=GRID, lw=.5))
    panel(fig, 133.5, 88.8, 'D', 'Sampled structure diversity')
    text_mm(fig, 138.5, 83.7, f'NEAT1_1  ·  3,756 nt\n1,000 samples; beam width {d["beam"]}',
            fontsize=6, va='top', linespacing=1.4)
    text_mm(fig, 138.5, 74.5, 'Scaffold G: threshold 0.7', fontsize=5.7, va='top')
    for k, (comparison, title) in enumerate([('sample_sample', 'Sample–sample'),
                                            ('sample_MFE', 'Sample–MFE'),
                                            ('sample_MEA', 'Sample–MEA')]):
        y = 50.8 - k * 20
        mean = float(next(r['mean'] for r in d['summary'] if r['space'] == 'extension'
                          and r['beam'] == str(d['beam']) and r['comparison'] == comparison))
        text_mm(fig, 138.5, y + 16.4, title, fontsize=6.3, weight='bold', va='top')
        text_mm(fig, 181.5, y + 16.4, f'{mean:.3f}',
                fontsize=5.5, ha='right', va='top')
        ax = axis_mm(fig, 139.5, y, 41, 12.8)
        h = [r for r in d['hist'] if r['space'] == 'extension' and r['beam'] == str(d['beam']) and r['comparison'] == comparison]
        edges = np.array([float(r['left']) for r in h] + [float(h[-1]['right'])])
        counts = np.array([int(r['count']) for r in h])
        density = counts / counts.sum() / np.diff(edges)
        ax.stairs(density, edges, fill=True, color=GRID, lw=0)
        ax.stairs(density, edges, color=INK, lw=.75)
        ax.set(xlim=(.3, .95), ylim=(0, 8), yticks=[0, 5], xticks=[.4, .6, .8])
        clean_axis(ax)
        if k < 2:
            ax.set_xticklabels([])
        if k == 1:
            ax.set_ylabel('Density', fontsize=5.7, labelpad=1)
    text_mm(fig, 160, 7, 'Jaccard distance of pairs outside G', fontsize=5.6, ha='center', va='top')
    text_mm(fig, 138.5, 3.5, 'Numbers above: mean distances', fontsize=5.4, va='top')
    return arc_counts


def boolean_runs(mask):
    starts = np.flatnonzero(np.diff(np.r_[False, mask, False].astype(int)) == 1)
    ends = np.flatnonzero(np.diff(np.r_[False, mask, False].astype(int)) == -1)
    return zip(starts, ends)


def stability_sentence(d):
    if all(s['every_window_change_below_0_05'] for s in d['stability'].values()):
        return ('Each highlighted value changes by <0.05 at b=400 with G (pale blue) or b=100 with empty G, relative to b=200.',
                'b=200と比べ、各注目窓の変化はGありのb=400（薄青）、空のGのb=100でいずれも0.05未満。')
    g = d['stability']['G0.7']['max_absolute_window_change']
    e = d['stability']['emptyG']['max_absolute_window_change']
    return (f'Relative to b=200, the largest window changes are {g:.4g} at b=400 with G (pale blue) and {e:.4g} at b=100 with empty G.',
            f'b=200に対する窓ごとの最大変化はGありのb=400（薄青）で{g:.4g}、空のGのb=100で{e:.4g}。')


def caption(d):
    g = d['diagnostics']['G0.7']['focus_mean']
    e = d['diagnostics']['emptyG']['focus_mean']
    region, focal = d['region'], d['focus']
    windows = focal[1] - focal[0] + 1
    stability_en, stability_ja = stability_sentence(d)
    return f'''<!-- 日本語：NEAT1 の予測を対応づける。G は全予測構造で保持する対の集合で、LinearPartition の確率0.7以上かつ連続2対以上の部分から作る。beam 幅は計算中に保持する候補を制限し、B–D の主曲線は {d['beam']}。B・C は全長計算の局所表示である。 -->
**Figure 4. Predicted pairing, local structure, and window accessibility in NEAT1.**
The scaffold G is a set of base pairs kept in every predicted structure.
The threshold-0.7 scaffold uses LinearPartition pairs with probability at least 0.7 in runs of two or more stacked pairs.
Beam width limits the candidates retained during calculation.
The main results in panels B–D use beam width {d['beam']}.
Panels B and C show local views of full-transcript calculations.

<!-- 日本語：A は全長22,743塩基に対する表示区間と開始{focal[0]:,}–{focal[1]:,}の{windows}個の8塩基窓。B は両端が{region[0]:,}–{region[1]:,}内にある確率0.01以上の対を示し、左右は閾値0.7のGありと空の条件。灰色はGの対、青色はGと交差しない新しい対で、濃さは確率を表す。 -->
**(A)** Location of the interval in NEAT1_2 (22,743 nucleotides, nt).
The highlight marks {windows} eight-base windows starting at {focal[0]:,}–{focal[1]:,}.
**(B)** Base-pairing probabilities with G (left) and an empty scaffold (right).
Gray marks pairs in G; blue marks new pairs without crossings of G.
Darker arcs have higher probabilities.
Arcs show only pairs with both bases in {region[0]:,}–{region[1]:,} and probability at least 0.01.

<!-- 日本語：C の LinearPKCapR は各塩基がstemまたは六つの非対合contextに属する確率、LinearPKRaccess は窓内8塩基が同時に非対合である確率を示す。網掛けはGの対合塩基を含みaccessibilityがゼロの窓。注目{windows}窓はすべてGで自由な塩基のみを含み、b={d['beam']}の平均はGあり{g:.4f}、空のGで{e:.4f}。{stability_ja}覆う塩基は{focal[1] + WINDOW - 1:,}まで、開始位置は{focal[1]:,}までである。 -->
**(C)** LinearPKCapR gives the probability that each base lies in a paired stem or one of six unpaired structural contexts (defined in Methods).
LinearPKRaccess gives window accessibility, the probability that all eight bases are unpaired together.
Hatching marks windows containing a base paired in G; their accessibility is zero.
All {windows} highlighted windows contain only bases unpaired in G.
Their means at b={d['beam']} are {g:.4f} with G and {e:.4f} with an empty scaffold.
{stability_en}
The highlighted bases extend to {focal[1] + WINDOW - 1:,}; window starts end at {focal[1]:,}.

<!-- 日本語：D は別のアイソフォーム NEAT1_1 全長3,756塩基の1,000標本を、標本同士および最小自由エネルギーMFE・期待正確度最大化MEA構造と比較する。MEAは対合確率との期待一致度を最大化する。共通Gの727対を除き、Jaccard距離は1から共有対数／和集合対数を引く。標本間499,500距離、各参照への1,000距離を幅0.01・面積1の分布とし、数値は正確な平均を示す。MEA重みは1。有限beamのMFE探索集合は異なりうる。モデル予測の比較である。 -->
**(D)** NEAT1_1 (3,756 nt): 1,000 sampled structures, compared with each other and with minimum-free-energy (MFE) and maximum expected accuracy (MEA) structures.
MEA selects a structure that maximizes expected agreement with pairing probabilities.
Distances exclude the 727 common pairs in G.
Jaccard distance is one minus the shared-pair count divided by the union-pair count.
Each distribution uses 499,500 distinct unordered sample pairs or 1,000 sample-to-reference distances.
Histograms have bin width 0.01 and area one; numbers give exact means.
MEA uses accuracy weight 1.
Finite-beam MFE may search a different structure set.
This figure compares model predictions.
'''


def short_caption(d):
    g = d['diagnostics']['G0.7']['focus_mean']
    e = d['diagnostics']['emptyG']['focus_mean']
    region, focal = d['region'], d['focus']
    stability_en, stability_ja = stability_sentence(d)
    return f'''<!-- 日本語：NEAT1の主な予測をbeam幅{d['beam']}で示す。Gは全構造で保持する対。beam幅は計算中に保持する候補を制限する。AはNEAT1_2の{region[0]:,}–{region[1]:,}、BはGあり／空の対合確率で区間内の0.01以上の対を示す。 -->
**Figure 4. NEAT1 predictions at beam width {d['beam']}.**
Scaffold G contains pairs kept in every structure.
Beam width limits retained candidates.
**(A)** NEAT1_2 interval {region[0]:,}–{region[1]:,}.
**(B)** Pairing probabilities with G or empty G; arcs require probability ≥0.01.

<!-- 日本語：CのLinearPKCapRは各塩基の七つの局所構造の確率、LinearPKRaccessは8塩基が同時に非対合の確率を示す。注目窓の開始は{focal[0]:,}–{focal[1]:,}で全塩基がGで自由。b={d['beam']}の平均はGあり{g:.4f}／空{e:.4f}。{stability_ja}網掛けはGの対合塩基を含む窓。 -->
**(C)** LinearPKCapR gives seven structural-context probabilities per base.
LinearPKRaccess gives accessibility: the probability that eight window bases are unpaired together.
Highlighted windows start at {focal[0]:,}–{focal[1]:,}; all bases are unpaired in G.
Their means are {g:.4f} with G and {e:.4f} with empty G.
{stability_en}
Hatching marks windows containing scaffold-paired bases.

<!-- 日本語：DはNEAT1_1の1,000標本。Gを除き、標本間およびMFE・MEA構造との距離を示す。Jaccard距離は1から共有対数／和集合対数を引く。 -->
**(D)** NEAT1_1: 1,000 samples, compared pairwise and with minimum-free-energy (MFE) or maximum expected accuracy (MEA) structures.
Distances exclude G.
Jaccard distance is one minus shared pairs divided by all distinct pairs.
'''


def source_note(d, counts):
    g = d['diagnostics']['G0.7']['focus_mean']
    e = d['diagnostics']['emptyG']['focus_mean']
    means = d['sample_checks'][0]['exact_means']
    region, focal = d['region'], d['focus']
    windows = focal[1] - focal[0] + 1
    gs, es = d['stability']['G0.7'], d['stability']['emptyG']
    checks = d['access_checks']
    crossing_note = ('No orange arc passes the display cutoff here.' if not counts['G0.7']['crossing']
                     else f"{counts['G0.7']['crossing']} orange arcs pass the display cutoff.")
    checkpoint_note = ('''That snapshot predates the empty-G b=200 BPP result.
For this one missing result, the script reads only the named `NEAT1_2_emptyG-beam200_outer-LinearCParty-rep0` job from the saved BPP checkpoint.
It verifies success, archive boundaries, byte counts, and hashes before using its `stdout.txt` and `bpp.tsv`.'''
                       if d['supplemental_checkpoints'] else
                       'All selected BPP results are included in the fixed E3 snapshot.')
    bpp_note = ('except for empty G at b=200, which uses the supplemental checkpoint above'
                if d['supplemental_checkpoints'] else 'for both selected conditions')
    return f'''# Figure 4: NEAT1 (2026-10-16)

[PDF](figure4.pdf) / [PNG](figure4.png) / [caption](caption.md) / [short caption](caption_short.md) / [provenance](provenance.json)

The PDF is exactly 183 × 90 mm. The PNG is rendered at 400 dpi.
The script is [make_figure4.py](../../../scripts/make_figure4.py).
From the paper repository run:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python experiments/scripts/make_figure4.py --beam {d['beam']} --region {region[0]} {region[1]} --focus {focal[0]} {focal[1]}
```

From the experiment repository, omit the `experiments/` prefix.
The `--beam` option selects the main width for B, C, and D; its default is 200.
The saved sampling results allow 100 or 200. Use `--out` for an alternative figure directory.
The `--region START END` option sets the display interval (default: 21,000–21,300).
The `--focus START END` option sets the highlighted window starts (default: 21,134–21,143).
Both ranges are one-based and inclusive; all highlighted window bases must fit within the display interval and be unpaired in G.
The author selected the new interval and beam 200 on 2026-10-16 after checking the accessibility changes at other widths.
The pale-blue G curve in C always uses b=400; empty-G b=100 is checked numerically and is not drawn.
The script reads the fixed E3 analysis snapshot and streams recorded archive members.
{checkpoint_note}
It does not scan other runs, run RNA software, extract archives, or update the E3 snapshot.

## Panel sources

All table links below refer to [the E3 analysis snapshot](../../local/e3-neat1-20261016/README.md).
`provenance.json` records the selected coordinates and beam, input and table hashes, every consumed archive path, member, hash, and saved command, and the supplemental checkpoint snapshot and selected job.
It also records the four focused accessibility series and means, the per-window width comparisons, and the additional b=400 curve.

| Panel | Direct data | Saved run and conditions |
|---|---|---|
| A | [NEAT1_2.txt](../../../data/neat1-20261016/NEAT1_2.txt), [input manifest](../../../data/neat1-20261016/manifest.json) | NR_131012.1, 22,743 nt; coordinates are one-based and inclusive. The whole-transcript track is schematic. |
| B | Streamed `bpp.tsv`; [bpp_pairs.tsv.gz](../../local/e3-neat1-20261016/bpp_pairs.tsv.gz) provides an independent table cross-check {bpp_note}. | `e3-neat1-bpp-e2cecbb4-20261016` and `e3-neat1-variants-bpp-e2cecbb4-20261016`, `neat1_long/000001`; NEAT1_2, threshold-0.7 G versus empty G, beam {d['beam']}. |
| C | [position_profiles.tsv.gz](../../local/e3-neat1-20261016/position_profiles.tsv.gz), [accessibility.tsv.gz](../../local/e3-neat1-20261016/accessibility.tsv.gz), cross-checked against streamed LinearPKCapR/LinearPKRaccess stdout | Corresponding `e3-neat1-[variants-]capr/raccess-e2cecbb4-20261016` runs, `neat1_long/000001`; the same two conditions and main beam as B, plus G b=400 accessibility (drawn) and empty-G b=100 accessibility (numerical check). |
| D | [sampling_histograms.tsv](../../local/e3-neat1-20261016/sampling_histograms.tsv), [sampling_summary.tsv](../../local/e3-neat1-20261016/sampling_summary.tsv), [sampling_reference_distances.tsv.gz](../../local/e3-neat1-20261016/sampling_reference_distances.tsv.gz), [representative_distances.tsv](../../local/e3-neat1-20261016/representative_distances.tsv) | `e3-neat1-sampling/bpp/mea-e2cecbb4-20261016` runs, `neat1_short/000001`; NEAT1_1 (NR_028272.1), threshold-0.7 G, beam {d['beam']} only, 1,000 samples, seed 20261016, MEA accuracy weight 1. Every distance was recomputed from streamed structures. |

## Scientific scope and display choices

- B and C are local views of calculations on the full 22,743-nt transcript. They do not use a refolded fragment.
- Both bases of a displayed arc must lie in {region[0]:,}–{region[1]:,}. B draws probabilities ≥0.01; the original BPP export cutoff is 0.001. For G, {counts['G0.7']['scaffold']} scaffold pairs and {counts['G0.7']['new_nested']} new pairs are drawn; for empty G, {counts['emptyG']['new_nested']} pairs are drawn. {counts['G0.7']['external_omitted']} / {counts['emptyG']['external_omitted']} pairs at the display cutoff touch the interval but have their other base outside and are omitted. Contexts and accessibility still include full-transcript contributions.
- Gray, blue, and orange follow Figure 1: scaffold pairs, new pairs without G crossings, and new pairs crossing G. {crossing_note} The seven-context colors match Figure 1 exactly. No small PseudoLoop value is inflated for visibility.
- C shows the seven LinearPKCapR probabilities as per-base steps with no smoothing. Their sum is one. LinearPKRaccess supplies the 8-nt accessibility values. Hatch marks identify windows containing at least one scaffold-paired base. Their exact accessibility is zero; microscopic raw numerical residuals are checked before assigning this exact condition.
- Every highlighted window starting at {focal[0]:,}–{focal[1]:,} contains only bases unpaired in G. The {windows}-window means at b={d['beam']} are {g:.12f} with G and {e:.12f} with empty G. The {windows + WINDOW - 1} covered base positions run from {focal[0]:,} to {focal[1] + WINDOW - 1:,}. B/context boundaries cover these bases; the accessibility highlight covers the {windows} starts.
- The pale-blue G b=400 curve is drawn behind the black main accessibility curve with lower opacity and a slightly wider line. The line key distinguishes the widths. Agreement is checked per highlighted window, not just between means; changes elsewhere in the displayed interval are not claimed to be small.
- The focus was selected exploratorily after comparing accessibility and checking its sensitivity to the saved widths; it is not a pre-specified biological region. The previous focus at 20,451–20,500 had nearly equal means at b=200 (0.0136 / 0.0137); its earlier b=100 contrast reflected width sensitivity with G. The author explicitly requested the replacement. No experimental validation or causal pseudoknot effect is inferred.
- D is a whole-transcript analysis of the shorter NEAT1_1 isoform, separate from the NEAT1_2 interval. Removing G prevents its 727 shared pairs from reducing the apparent diversity. The selected beam has 1,000 distinct sampled structures and 499,500 unordered sample pairs. Pairs sharing a sample are not independent replicates.
- D uses 100 bins over [0,1], width 0.01, with density normalized separately to area one. The displayed distance axis is 0.30–0.95; all observed distances lie within it. Means come from exact distances, not bin centers. Only beam {d['beam']} is plotted, with a solid outline and gray fill. The mean distances are {means['sample_sample']:.12f} (sample–sample), {means['sample_MFE']:.12f} (sample–MFE), and {means['sample_MEA']:.12f} (sample–MEA). Finite-beam MFE may use a different search set from the sampled ensemble.

## Focused accessibility and beam-width sensitivity

| Scaffold | Beam width | Mean over the {windows} highlighted windows |
|---|---:|---:|
| G, threshold 0.7 | 200 | {checks['G0.7'][200]['focus_mean']:.16g} |
| G, threshold 0.7 | 400 | {checks['G0.7'][400]['focus_mean']:.16g} |
| Empty G | 100 | {checks['emptyG'][100]['focus_mean']:.16g} |
| Empty G | 200 | {checks['emptyG'][200]['focus_mean']:.16g} |

For G, the largest absolute change among the highlighted windows between b=200 and b=400 is {gs['max_absolute_window_change']:.16g}; the absolute change of the mean is {gs['absolute_mean_change']:.16g}.
For empty G, the corresponding changes between b=100 and b=200 are {es['max_absolute_window_change']:.16g} and {es['absolute_mean_change']:.16g}.
These comparisons describe sensitivity to the saved beam/local-limit settings (`--max-hairpin` and `--max-case3-depth` each equal b), rather than agreement with an exact calculation.

## Development checks

Input hashes, archived member hashes and byte counts, saved beam/ranking/local-limit/seed/MEA settings, transcript identity, all per-base context probabilities, all width-8 accessibility values, available BPP tables, and scaffold constraints are checked before plotting.
The supplemental empty-G BPP also passes probability bounds and incident-mass checks against the full LinearPKCapR stem probabilities; no absent table comparison is claimed.
All 499,500 pairwise distances and 2,000 sample-to-reference distances at the selected beam are independently recomputed with integer bitsets and compared with the saved tables.
Histogram counts and exact means must agree. The source analysis used sparse pair-incidence matrix products.
The script checks text bounds and overlapping text. Visual PNG review and the final PDF/PNG dimensions are recorded in the paper repository's `docs/progress.md`.
No manuscript file, bibliography, external implementation, RNA CLI, scheduler job, or Git commit is part of this script.

The author's 2026-10-16 choices settle the display interval and main beam selection for this figure.
'''


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data', type=Path, default=ROOT / 'results/local/e3-neat1-20261016')
    ap.add_argument('--out', type=Path, default=ROOT / 'results/figures/figure4_20261016')
    ap.add_argument('--beam', type=int, choices=(100, 200), default=200,
                    help='main beam width for panels B, C, and D (default: 200)')
    ap.add_argument('--region', type=int, nargs=2, metavar=('START', 'END'), default=DEFAULT_REGION,
                    help='one-based inclusive display interval (default: 21000 21300)')
    ap.add_argument('--focus', type=int, nargs=2, metavar=('START', 'END'), default=DEFAULT_FOCUS,
                    help='one-based inclusive highlighted window starts (default: 21134 21143)')
    args = ap.parse_args()
    if not 1 <= args.region[0] < args.region[1] <= 22743 - WINDOW + 1:
        ap.error('--region must be increasing and contain valid width-8 window starts in NEAT1_2')
    if not args.region[0] <= args.focus[0] <= args.focus[1] <= args.region[1] - WINDOW + 1:
        ap.error('--focus windows must fit entirely within --region')
    d = read_data(args.data.resolve(), args.beam, tuple(args.region), tuple(args.focus))
    args.out.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(WIDTH * MM, HEIGHT * MM))
    counts = draw(fig, d)
    # No tight bounding box: the PDF remains exactly 183 x 90 mm.
    fig.savefig(args.out / 'figure4.pdf', facecolor='white')
    fig.savefig(args.out / 'figure4.png', dpi=400, facecolor='white')
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    outside = []
    texts = []
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_visible() and t.get_text():
            bb = t.get_window_extent(renderer)
            texts.append((t.get_text(), bb))
            if not fig.bbox.contains(bb.x0, bb.y0) or not fig.bbox.contains(bb.x1, bb.y1):
                outside.append(t.get_text())
    assert not outside, 'Text outside figure: ' + repr(outside)
    overlap = []
    for k, (a, bb) in enumerate(texts):
        for b, other in texts[k + 1:]:
            if min(bb.x1, other.x1) - max(bb.x0, other.x0) > .5 and min(bb.y1, other.y1) - max(bb.y0, other.y0) > .5:
                overlap.append([a, b])
    assert not overlap, 'Overlapping text: ' + repr(overlap)
    plt.close(fig)
    (args.out / 'caption.md').write_text(caption(d))
    (args.out / 'caption_short.md').write_text(short_caption(d))
    (args.out / 'README.md').write_text(source_note(d, counts))
    metadata = dict(figure_size_mm=[WIDTH, HEIGHT], beam_width=args.beam,
                    panel_beam_widths={p: args.beam for p in ('B', 'C', 'D')},
                    accessibility_overlay=dict(panel='C', scaffold='G0.7', beam_width=400,
                                               color=NEW_COLOR, alpha=.4, linewidth=1.4),
                    context_program='LinearPKCapR', accessibility_program='LinearPKRaccess',
                    region_one_based_inclusive=d['region'],
                    focus_window_starts_one_based_inclusive=d['focus'], window_length=WINDOW,
                    focus_bases_one_based_inclusive=[d['focus'][0], d['focus'][1] + WINDOW - 1],
                    bpp_display_cutoff=PAIR_CUTOFF, bpp_export_cutoff=.001,
                    arc_counts=counts, diagnostics=d['diagnostics'], sampling_validation=d['sample_checks'],
                    accessibility_validation=d['access_checks'], beam_width_sensitivity=d['stability'],
                    source_tables_and_inputs=d['sources'], streamed_artifacts=d['artifacts'],
                    supplemental_checkpoints=d['supplemental_checkpoints'],
                    script_sha256=sha(Path(__file__).read_bytes()), command=sys.argv,
                    software_revision='e2cecbb4', coordinate_convention='one-based inclusive',
                    validation=dict(text_inside_canvas=True, text_overlap=False, smoothing=False, raw_archives_extracted=False,
                                    rna_cli_called=False, all_pairwise_distances_recomputed=True,
                                    focus_windows_contain_only_scaffold_free_bases=True))
    (args.out / 'provenance.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps(dict(output=str(args.out), figure_size_mm=[WIDTH, HEIGHT], beam_width=args.beam, arc_counts=counts,
                          focus_mean={k: v['focus_mean'] for k, v in d['diagnostics'].items()},
                          accessibility_validation=d['access_checks'], beam_width_sensitivity=d['stability'],
                          sampling_validation=d['sample_checks']), indent=2))


if __name__ == '__main__':
    main()
