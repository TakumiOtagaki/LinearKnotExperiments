#!/usr/bin/env python3
"""Make the 165-mm-wide NEAT1 Supplementary Figure S3 from saved E3 outputs.

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
import re
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import xml.etree.ElementTree as ET

import numpy as np

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'neat1-s3-mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import StrMethodFormatter

from analyze_e3_neat1 import crossing, pairs_of, parse_access, parse_bpp, parse_context, structure_from_stdout

ROOT = Path(__file__).resolve().parents[1]
MM = 1 / 25.4
WIDTH, HEIGHT = 165, 145
# Vertical coordinates use the 165-mm layout grid; fonts stay at their native size.
LAYOUT_HEIGHT = 165
PNG_DPI = 300
FONT_SIZE = 7
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
    # These are saved joint window probabilities, not single-base contexts.
    byloop_table = data / 'byloop_regions.tsv'
    byloop = None
    if byloop_table.is_file():
        selected = [r for r in rows(byloop_table) if r['iso'] == 'NEAT1_2'
                    and int(r['beam']) == beam and int(r['window']) == WINDOW
                    and (int(r['start_min']), int(r['start_max'])) == focus_starts]
        assert len(selected) == 1, 'Missing or ambiguous saved loop-type focus means'
        byloop = {c: float(selected[0][c + '_mean']) for c in CONTEXTS if c != 'Stem'}
        byloop['total'] = float(selected[0]['total_mean'])
        assert all(v >= 0 for v in byloop.values())
        assert abs(sum(v for c, v in byloop.items() if c != 'total') - byloop['total']) < 1e-12
        sources.append(dict(path=str(byloop_table.relative_to(ROOT)), sha256=sha(byloop_table.read_bytes())))
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
        kind = re.sub(r'-[0-9a-f]{8}-\d{8}$', '', a['archive'].split('/')[2]
                      .removeprefix('e3-neat1-').removeprefix('variants-'))
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

    if byloop is not None:
        assert abs(byloop['total'] - diagnostics['G0.7']['focus_mean']) < 1e-12

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
                sources=sources, artifacts=consumed, supplemental_checkpoints=checkpoints,
                focus_byloop=byloop)


def axis_mm(fig, x, y, w, h):
    return fig.add_axes([x / WIDTH, y / LAYOUT_HEIGHT, w / WIDTH, h / LAYOUT_HEIGHT])


def text_mm(fig, x, y, s, **kwargs):
    return fig.text(x / WIDTH, y / LAYOUT_HEIGHT, s, **kwargs)


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
    windows = focal[1] - focal[0] + 1
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': FONT_SIZE,
                         'text.color': INK, 'axes.labelcolor': INK,
                         'axes.edgecolor': MUTED, 'axes.linewidth': .5,
                         'xtick.color': MUTED, 'ytick.color': MUTED,
                         'xtick.labelsize': FONT_SIZE, 'ytick.labelsize': FONT_SIZE,
                         'xtick.major.size': 2, 'ytick.major.size': 2,
                         'xtick.major.width': .5, 'ytick.major.width': .5,
                         'pdf.fonttype': 42, 'hatch.linewidth': .35})
    n = d['inputs']['NEAT1_2']['n']
    panel(fig, 1.5, 163.5, 'A', f'NEAT1_2 · {n:,} nt · local predictions at b={d["beam"]}')
    text_mm(fig, 163, 163.5, f'View: {region[0]:,}–{region[1]:,}', ha='right', va='top')
    ax = axis_mm(fig, 20, 157, 143, 2.5)
    ax.plot([1, n], [0, 0], color=G_COLOR, lw=2, solid_capstyle='butt')
    ax.plot(region, [0, 0], color=INK, lw=3, solid_capstyle='butt')
    ax.plot([np.mean(focal)], [0], marker='v', color=INK, ms=3)
    ax.set(xlim=(1, n), ylim=(-1, 1), xticks=[1, 10000, n],
           xticklabels=['1', '10,000', f'{n:,}'], yticks=[])
    ax.spines[:].set_visible(False)
    ax.tick_params(axis='x', length=0, pad=1)
    ax.get_xticklabels()[-1].set_ha('right')
    text_mm(fig, 20, 150.5,
            f'{windows} highlighted starts: {focal[0]:,}–{focal[1]:,}; '
            f'window bases: {focal[0]:,}–{focal[1] + WINDOW - 1:,}', va='top',
            bbox=dict(facecolor=FOCUS_COLOR, edgecolor='none', pad=1.2))

    panel(fig, 1.5, 147.5, 'B', 'Base-pairing probabilities')
    columns = [('G0.7', 20, 'Scaffold G: threshold 0.7'), ('emptyG', 96, 'Empty scaffold')]
    column_width = 67
    arc_counts = {}
    arc_height = max(((int(r['j']) - int(r['i'])) / 2 for method, _, _ in columns
                     for r in d['pairs'][method] if float(r['probability']) >= PAIR_CUTOFF
                     and region[0] <= int(r['i']) < int(r['j']) <= region[1]), default=1) * 1.05
    styles = {'scaffold': 'solid', 'new_nested': (0, (1, 1)), 'crossing': (0, (4, 2))}
    for method, x, label in columns:
        text_mm(fig, x, 143.5, label, weight='bold', va='top')
        ax = axis_mm(fig, x, 127.5, column_width, 13.5)
        focus(ax, focal)
        paths, colors, linestyles = [], [], []
        counts = {'scaffold': 0, 'new_nested': 0, 'crossing': 0, 'external_omitted': 0}
        for r in sorted(d['pairs'][method], key=lambda r: float(r['probability'])):
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
            paths.append(np.column_stack(((i + j) / 2 + (j - i) / 2 * np.cos(t),
                                          (j - i) / 2 * np.sin(t))))
            colors.append(matplotlib.colors.to_rgba(color, p))
            linestyles.append(styles[category])
            counts[category] += 1
        ax.add_collection(LineCollection(paths, colors=colors, linewidths=.65, linestyles=linestyles))
        ax.axhline(0, color=MUTED, lw=.4)
        ax.set(xlim=(region[0] - .5, region[1] + .5), ylim=(-.02 * arc_height, arc_height), xticks=[], yticks=[])
        ax.spines[:].set_visible(False)
        arc_counts[method] = counts
        text_mm(fig, x, 124.5, f'{counts["external_omitted"]:,} pairs leave the view (omitted)', va='top')
    handles = [Line2D([], [], color=G_COLOR, lw=1, label='Pair in G'),
               Line2D([], [], color=NEW_COLOR, lw=1, ls=styles['new_nested'], label='New pair, no G crossing')]
    if any(c['crossing'] for c in arc_counts.values()):
        handles.append(Line2D([], [], color=CROSS_COLOR, lw=1, ls=styles['crossing'], label='Pair crossing G'))
    fig.legend(handles=handles, loc='center left', bbox_to_anchor=(20 / WIDTH, 120.2 / LAYOUT_HEIGHT),
               frameon=False, ncol=len(handles), fontsize=FONT_SIZE, borderaxespad=0,
               handlelength=1.7, columnspacing=1, handletextpad=.5)
    text_mm(fig, 163, 120.2, 'Opacity = pair probability', ha='right', va='center')

    # Separate, directly named tracks preserve every context in grayscale.
    # All tracks retain the same 0--1 scale; no tiny component is enlarged.
    panel(fig, 1.5, 115.5, 'C', 'Structural contexts and 8-nt window accessibility')
    text_mm(fig, 163, 115.5, 'Each context row: 0–1 (bottom–top)', ha='right', va='top')
    pos = np.arange(region[0], region[1] + 1)
    for row, (context, color) in enumerate(zip(CONTEXTS, CONTEXT_COLORS)):
        y = 107.5 - row * 3.8
        text_mm(fig, 18.2, y + 1.65, context, ha='right', va='center')
        for method, x, _ in columns:
            ax = axis_mm(fig, x, y, column_width, 3.3)
            values = d['profile'][method][context][region[0] - 1:region[1]]
            ax.axhline(1, color=GRID, lw=.35, ls=(0, (1, 2)))
            ax.fill_between(pos, 0, values, step='mid', facecolor=color, linewidth=0)
            ax.plot(pos, values, color=INK, lw=.45, drawstyle='steps-mid')
            for xx in (focal[0] - .5, focal[1] + WINDOW - .5):
                ax.axvline(xx, color=INK, lw=.45, ls=(0, (2, 2)))
            ax.set(xlim=(region[0] - .5, region[1] + .5), ylim=(0, 1), xticks=[], yticks=[])
            ax.spines[:].set_visible(False)
            ax.spines['bottom'].set_visible(True)
            ax.spines['bottom'].set_color(GRID)
            if np.max(values) < .001:
                text_mm(fig, x + 1, y + 1.65, '<0.001 in this view', va='center')

    for method, x, _ in columns:
        mean = d['diagnostics'][method]['focus_mean']
        text_mm(fig, x, 81.3, f'8-nt accessibility · highlighted mean {mean:.4f}', va='top')
        ax = axis_mm(fig, x, 66, column_width, 10)
        focus(ax, focal, starts=True)
        if method == 'G0.7':
            mask = d['fixed_window'][region[0] - 1:region[1]]
            for a, b in boolean_runs(mask):
                ax.axvspan(region[0] + a - .5, region[0] + b - .5,
                           facecolor='#e8e8e4', edgecolor=G_COLOR, hatch='////', lw=0, zorder=1)
        ax.plot(pos, d['access'][method][region[0] - 1:region[1]],
                color=INK, lw=.85, drawstyle='steps-mid', zorder=3)
        ax.plot([focal[0] - .5, focal[1] + .5], [0, 0], color=INK, lw=1.1, zorder=2)
        ax.set(xlim=(region[0] - .5, region[1] + .5), ylim=(-.015, 1.05), yticks=[0, 1],
               xticks=np.linspace(*region, 4))
        ax.xaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))
        ax.get_xticklabels()[0].set_ha('left')
        ax.get_xticklabels()[-1].set_ha('right')
        clean_axis(ax)
    text_mm(fig, 91.5, 60, 'Nucleotide position / window start', ha='center', va='center')
    fig.legend(handles=[Patch(facecolor='#e8e8e4', edgecolor=G_COLOR, hatch='////',
                              label='Window with a G-paired base: 0')],
               loc='center left', bbox_to_anchor=(96 / WIDTH, 45.8 / LAYOUT_HEIGHT),
               fontsize=FONT_SIZE, frameon=False, borderaxespad=0, handlelength=1.7)
    if d['focus_byloop'] is not None:
        text_mm(fig, 20, 56.6, 'Highlighted means by loop type (G)', weight='bold', va='center')
        for i, context in enumerate(('Hairpin', 'Bulge', 'Internal', 'Multiloop', 'Exterior', 'PseudoLoop')):
            value = d['focus_byloop'][context]
            formatted = (f'{value:.2e}' if 0 < value < .001 else f'{value:.3f}').replace('e-', 'e−')
            text_mm(fig, 20 + 34 * (i % 2), 52.8 - 3.5 * (i // 2),
                    f'{context} {formatted}', va='center')
    text_mm(fig, 96, 56.6, 'Highlighted-window sensitivity', weight='bold', va='center')
    for i, (method, label) in enumerate((('G0.7', 'G'), ('emptyG', 'Empty G'))):
        widths = d['stability'][method]['beam_widths']
        delta = d['stability'][method]['max_absolute_window_change']
        text_mm(fig, 96, 52.8 - 3.5 * i,
                f'{label}, b={widths[0]} / {widths[1]}: max change {delta:.3g}'.replace('e-', 'e−'), va='center')

    fig.add_artist(Line2D([1.5 / WIDTH, 163 / WIDTH], [44 / LAYOUT_HEIGHT] * 2,
                         transform=fig.transFigure, color=GRID, lw=.5))
    panel(fig, 1.5, 41.5, 'D', 'NEAT1_1: structure diversity')
    sample_check = d['sample_checks'][0]
    text_mm(fig, 20, 37.3, f'{d["inputs"]["NEAT1_1"]["n"]:,} nt; {sample_check["samples"]:,} samples; '
            f'scaffold G at threshold 0.7; b={d["beam"]}', va='top')
    for k, (comparison, title) in enumerate([('sample_sample', 'Sample–sample'),
                                            ('sample_MFE', 'Sample–MFE'),
                                            ('sample_MEA', 'Sample–MEA')]):
        x = 20 + k * 51
        mean = sample_check['exact_means'][comparison]
        h = [r for r in d['hist'] if r['space'] == 'extension'
             and r['beam'] == str(d['beam']) and r['comparison'] == comparison]
        edges = np.array([float(r['left']) for r in h] + [float(h[-1]['right'])])
        counts = np.array([int(r['count']) for r in h])
        density = counts / counts.sum() / np.diff(edges)
        assert np.isclose(np.sum(density * np.diff(edges)), 1)
        text_mm(fig, x, 32.8, f'{title} ({counts.sum():,})', weight='bold', va='top')
        ax = axis_mm(fig, x, 10, 41, 20)
        ax.stairs(density, edges, fill=True, color=GRID, lw=0)
        ax.stairs(density, edges, color=INK, lw=.75)
        ax.axvline(mean, color=INK, lw=.75, ls=(0, (3, 2)))
        ax.text(.97, .95, f'Mean {mean:.3f}', transform=ax.transAxes,
                fontsize=FONT_SIZE, ha='right', va='top')
        ax.set(xlim=(.3, .95), ylim=(0, 8), yticks=[0, 5], xticks=[.4, .6, .8])
        clean_axis(ax)
        if k == 0:
            ax.set_ylabel('Density', fontsize=FONT_SIZE, labelpad=3)
    text_mm(fig, 91.5, 3.2, 'Jaccard distance of pairs outside G', ha='center', va='center')
    return arc_counts

def boolean_runs(mask):
    starts = np.flatnonzero(np.diff(np.r_[False, mask, False].astype(int)) == 1)
    ends = np.flatnonzero(np.diff(np.r_[False, mask, False].astype(int)) == -1)
    return zip(starts, ends)


def caption(d, counts):
    region, focal = d['region'], d['focus']
    windows = focal[1] - focal[0] + 1
    sample = d['sample_checks'][0]
    fixed = len(d['inputs']['NEAT1_1']['fixed'])
    byloop_note = ''
    if d['focus_byloop'] is not None:
        q = d['focus_byloop']
        others = sum(v for c, v in q.items() if c not in ('total', 'Hairpin'))
        byloop_note = f'''\n<!-- 日本語：左下の小表はGを使った注目{windows}窓のjoint accessibilityを六つのloopの種類に分けた平均である。
合計は{q['total']:.3f}、Hairpinの寄与は{q['Hairpin']:.3f}、他の五成分の和は{others:.3f}である。 -->
The small table gives mean joint accessibility by six loop types for the {windows} highlighted windows with G.
Mean total accessibility is {q['total']:.3f}; Hairpin contributes {q['Hairpin']:.3f}, and the other five types sum to {others:.3f}.
'''
    return f'''<!-- 日本語：主な予測は beam 幅 b={d['beam']} を使う。
beam 幅は保持する候補を制限する。
scaffold G は全構造で保持する対である。 -->
**Figure S4. Predicted probabilities and structure diversity in NEAT1.**
The main predictions use beam width b={d['beam']}.
Beam width limits retained candidates.
Scaffold G contains pairs retained in every structure.

<!-- 日本語：A は NEAT1_2 の表示区間を示す。
注目する{windows}個の8塩基窓は{focal[0]:,}–{focal[1]:,}から始まり、塩基{focal[0]:,}–{focal[1] + WINDOW - 1:,}を覆う。 -->
**(A)** NEAT1_2 ({d['inputs']['NEAT1_2']['n']:,} nucleotides, nt).
The {windows} highlighted eight-base windows start at {focal[0]:,}–{focal[1]:,} and span bases {focal[0]:,}–{focal[1] + WINDOW - 1:,}.

<!-- 日本語：B の灰色実線は G、青色点線は G と交差しない新しい対、橙色破線はGと交差する新しい対を示す。
弧の不透明度は対確率であり、重なった弧は濃く見えることがある。
弧は確率{PAIR_CUTOFF}以上で両端が表示区間内にある対だけを示す。
区間外へ伸びる対は G ありで{counts['G0.7']['external_omitted']}、空の G で{counts['emptyG']['external_omitted']}対あり、省略する。
C はそれらの寄与も含む全長計算である。 -->
**(B)** Solid grey arcs show G, dotted blue arcs show new pairs that do not cross G, and dashed orange arcs show new pairs crossing G.
Arc opacity represents pair probability; overlapping arcs may appear darker.
Arcs require probability ≥{PAIR_CUTOFF} and both endpoints in view.
B omits pairs with a base outside the view: {counts['G0.7']['external_omitted']} with G and {counts['emptyG']['external_omitted']} with empty G.
C includes their full-transcript contributions.

<!-- 日本語：C の各 context は同じ0–1の尺度で示す。
0.001未満というラベルは表示区間内のすべての保存値がこの閾値未満であることを示す。
accessibility は窓内8塩基が同時に非対合である確率である。
網掛けは G の対合塩基を含み accessibility が0の窓を示す。
注目窓はすべて G で非対合の塩基だけを含む。
図内の平均と beam 幅・局所制限を変えた最大差は注目{windows}窓について示す。 -->
**(C)** Context tracks share a 0–1 scale (definitions in Methods).
Context values labelled <0.001 are below this threshold throughout the displayed interval.
Accessibility is the probability that all eight window bases are unpaired together.
Hatched windows contain G-paired bases and have zero accessibility.
All highlighted windows contain only bases unpaired in G.
Means and changes between settings refer to the {windows} highlighted windows.
{byloop_note}

<!-- 日本語：D は NEAT1_1 全長の{sample['samples']:,}標本で、共通 G の{fixed}対を距離から除く。
Jaccard 距離は1から共有対数を和集合の対数で割った値を引く。
MFE は最小自由エネルギー、MEA は最大期待正確度を意味し、MEA の重みは1である。
標本間の距離は{sample['sample_pairs']:,}個、各参照との距離は{sample['samples']:,}個を使う。
幅0.01のヒストグラムは各々の面積を1とし、破線は正確な平均、括弧内は距離数を示す。
有限 beam の MFE は異なる構造集合を探索しうる。 -->
**(D)** NEAT1_1 ({d['inputs']['NEAT1_1']['n']:,} nt), {sample['samples']:,} samples; distances exclude {fixed} shared G pairs.
Jaccard distance is one minus shared-pair count divided by union-pair count.
MFE denotes minimum free energy.
MEA denotes maximum expected accuracy and uses weight 1.
Sample–sample uses all {sample['sample_pairs']:,} unordered pairs; each sample-to-reference comparison uses {sample['samples']:,} distances.
Histograms use bin width 0.01 and are normalized separately to area one.
Dashed lines mark exact means.
Finite-beam MFE may search a different structure set.
'''

def source_note(d, counts):
    region, focal = d['region'], d['focus']
    windows = focal[1] - focal[0] + 1
    gs, es = (d['stability'][m] for m in ('G0.7', 'emptyG'))
    checks = d['access_checks']
    sample = d['sample_checks'][0]
    means = sample['exact_means']
    checkpoint_note = ('''The fixed E3 snapshot predates empty-G b=200 BPP.
Only the successful `NEAT1_2_emptyG-beam200_outer-LinearCParty-rep0` record and its committed `stdout.txt` and `bpp.tsv` members are read from the saved variants-BPP checkpoint.
Their archive boundaries, byte counts, and SHA-256 hashes are verified.'''
                       if d['supplemental_checkpoints'] else
                       'All selected BPP results are included in the fixed E3 snapshot.')
    byloop_note = ''
    if d['focus_byloop'] is not None:
        byloop_note = ('\n## Joint accessibility by loop type\n\n'
                      'The small table in C uses the saved `byloop_regions.tsv` row for NEAT1_2, '
                      f'b={d["beam"]}, w=8, starts {focal[0]}–{focal[1]}. '
                      'These are joint window probabilities; no component is inferred from single-base contexts. '
                      'The component sum is checked against the saved total and the independently validated focus curve.\n\n'
                      '| Loop type | Mean contribution |\n|---|---:|\n' +
                      ''.join(f'| {c} | {v:.17g} |\n' for c, v in d['focus_byloop'].items()))
    return f'''# Supplementary Figure S4: NEAT1 (2026-10-19, round four)

[PDF](figure_s3.pdf) / [PNG](figure_s3.png) / [caption](caption.md) / [provenance](provenance.json) / [validation](validation.json) / [focus values](focus_accessibility.tsv)

The PDF is exactly {WIDTH} × {HEIGHT} mm, for inclusion at the SI width of {WIDTH} mm.
The PNG is {PNG_DPI} dpi.
Ordinary text is {FONT_SIZE} pt; independent bold A–D labels are 9 pt.
No mathematical superscripts or reduced-size logarithmic ticks are used.
The wider two-column comparison and three bottom histograms replace the 183 × 90 mm Figure 4 layout without overwriting it.
C uses seven named probability tracks with colored areas and black outlines instead of a color-only stack.
Every track uses its original 0–1 scale; no small component is enlarged.

## Reproduce

From the paper repository:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python experiments/scripts/make_figure4.py --beam {d['beam']} --region {region[0]} {region[1]} --focus {focal[0]} {focal[1]} --out experiments/results/figures/figure_s3_neat1_20261017
```

From the experiment repository, omit the `experiments/` prefixes.
Requirements: Python, NumPy, Matplotlib, Pillow, pypdf, and Ghostscript (`gs`).
The default output is the new S3 directory.
The saved sampling results support main beam widths 100 and 200; this figure uses {d['beam']}.
`--region` and `--focus` use one-based inclusive positions.
All focus windows must fit inside the view and contain only bases unpaired in G.
The script streams saved tar members without extracting them, calling RNA software, or changing the analysis snapshot.
Temporary font caches and PDF inspection files use node-local temporary storage.

## Data and version

The scientific outputs remain those of LinearCParty `e2cecbb4`, not the ongoing resource rerun.
The [fixed E3 analysis](../../local/e3-neat1-20261016/README.md) supplies all table values.
`provenance.json` records the script hash, selected beam and positions, input/table hashes, saved commands, and every consumed archive/member hash.
{checkpoint_note}

| Panel | Input and saved source | Scope |
|---|---|---|
| A | [NEAT1_2 input](../../../data/neat1-20261016/NEAT1_2.txt), [manifest](../../../data/neat1-20261016/manifest.json) | NR_131012.1, {d['inputs']['NEAT1_2']['n']:,} nt; schematic full-transcript position track. |
| B | BPP archives `e3-neat1-[variants-]bpp-e2cecbb4-20261016`, `neat1_long/000001`; [BPP table](../../local/e3-neat1-20261016/bpp_pairs.tsv.gz) | G at threshold 0.7 versus empty G, b={d['beam']}; the independent BPP table lacks empty-G b=200, covered by the named checkpoint above. |
| C | Corresponding `capr`/`raccess` archives; [contexts](../../local/e3-neat1-20261016/position_profiles.tsv.gz), [accessibility](../../local/e3-neat1-20261016/accessibility.tsv.gz), [loop-type means](../../local/e3-neat1-20261016/byloop_regions.tsv) | Whole-transcript probabilities displayed locally; 8-nt total accessibility curves and the six saved mean loop-type contributions for the highlighted windows with G. Saved beam/local-limit comparisons are numerical. |
| D | `sampling`/`bpp`/`mea` archives, `neat1_short/000001`; [histograms](../../local/e3-neat1-20261016/sampling_histograms.tsv), [summaries](../../local/e3-neat1-20261016/sampling_summary.tsv), [reference distances](../../local/e3-neat1-20261016/sampling_reference_distances.tsv.gz), [representative distances](../../local/e3-neat1-20261016/representative_distances.tsv) | NEAT1_1, NR_028272.1, {d['inputs']['NEAT1_1']['n']:,} nt; {sample['samples']:,} samples, b={d['beam']}, seed 20261016, MEA weight 1. |

## Scientific interpretation

B and C display local results of full-transcript calculations, not refolded fragments.
Only arcs with probability ≥{PAIR_CUTOFF} and both endpoints in {region[0]:,}–{region[1]:,} are drawn; the export cutoff was 0.001.
G has {counts['G0.7']['scaffold']} displayed scaffold pairs and {counts['G0.7']['new_nested']} new non-crossing pairs; empty G has {counts['emptyG']['new_nested']} displayed pairs.
Pairs with only one endpoint in view are counted once and omitted: {counts['G0.7']['external_omitted']} with G, {counts['emptyG']['external_omitted']} with empty G.
Their contributions remain in C.
The numbers above are recomputed from the verified BPP exports.
Solid gray (#8a8a85) identifies G, dotted blue (#2a78d6) identifies new pairs without G crossings, and dashed orange (#eb6834) identifies pairs crossing G.
There are {counts['G0.7']['crossing']} displayed crossing pairs with G.
The seven context colors follow the shared palette: Stem #2a78d6, Hairpin #1baf7a, Bulge #eda100, Internal #e87ba4, Multiloop #6250d6, Exterior #b9b8b1, PseudoLoop #eb6834.
Direct row labels and black probability outlines allow identification without color.
A `<0.001 in this view` label is used only when every saved value in the displayed interval is below 0.001.

The {windows} highlighted starts are {focal[0]:,}–{focal[1]:,}; their windows cover {windows + WINDOW - 1} bases through {focal[1] + WINDOW - 1:,}.
Dashed boundaries in B/context rows enclose the bases; the pale accessibility shading encloses the starts.
Every highlighted window contains only bases unpaired in G.
Hatching in the accessibility graph and legend uses the same `////` pattern for windows containing any G-paired base.
Such windows have exact accessibility zero after checking their raw numerical residuals.

The focus was selected exploratorily after inspecting the saved predictions and their width sensitivity; it is not a pre-specified biological region.
The author's accepted focus and b=200 setting are unchanged.
Changing G also changes the conditional ensemble; this comparison does not isolate a causal pseudoknot effect or experimentally validate the predictions.

{byloop_note}
## Focused accessibility and sensitivity

| Scaffold | Beam width | Mean over the {windows} highlighted windows |
|---|---:|---:|
| G, threshold 0.7 | 200 | {checks['G0.7'][200]['focus_mean']:.16g} |
| G, threshold 0.7 | 400 | {checks['G0.7'][400]['focus_mean']:.16g} |
| Empty G | 100 | {checks['emptyG'][100]['focus_mean']:.16g} |
| Empty G | 200 | {checks['emptyG'][200]['focus_mean']:.16g} |

The maximum per-window change is {gs['max_absolute_window_change']:.16g} for G (b=200/400) and {es['max_absolute_window_change']:.16g} for empty G (b=100/200).
All four series are saved in `focus_accessibility.tsv`.
Only the highlighted windows are claimed to be stable within 0.05 across these saved settings.
The local hairpin and case-3 depth limits each equal `max(b, 30)` in these runs; the comparison measures beam/local-limit sensitivity, not agreement with exact calculations.

D removes {len(d['inputs']['NEAT1_1']['fixed'])} common G pairs before computing Jaccard distances.
It uses all {sample['sample_pairs']:,} unordered sample pairs and {sample['samples']:,} distances to each reference structure.
Pairs sharing a sample are not independent replicates.
Each histogram uses 100 bins of width 0.01 over [0, 1] and is independently normalized to area one.
All observed distances lie within the displayed 0.30–0.95 range.
Dashed vertical lines and labels use exact means, not bin centers: {means['sample_sample']:.12f}, {means['sample_MFE']:.12f}, and {means['sample_MEA']:.12f}.
Finite-beam MFE may search a different structure set from the sampled ensemble.

## Checks and review coverage

- S3-1 / G1: {WIDTH}-mm native width; two {67}-mm data columns; D spans the bottom; ordinary text {FONT_SIZE} pt.
- S3-2: omitted outside-endpoint pair counts in B; caption explains full-transcript contributions in C.
- S3-3: means explicitly refer to highlighted starts; base and start boundaries are defined separately.
- S3-4: largest per-window changes replace the obscured b=400 overlay; all saved-width means remain in the source table.
- S3-5: exact mean lines/labels, distance counts, independent density normalization, bin width, and MEA weight are explicit.
- S3-6: graph and legend use identical hatching.
- G2: named context rows with black outlines; pair line styles distinguish G and new pairs.
- G3: A–D use independent 9-pt bold letters at the same left position.

Before plotting, the script verifies input and artifact hashes, byte counts, archive boundaries, saved beam/ranking/local-limit/seed/MEA settings, and scaffold constraints.
Full-transcript context and accessibility arrays must match the independently saved tables.
Available BPP tables must match; the supplemental empty-G BPP must also satisfy incident-mass and context bounds.
All {sample['sample_pairs']:,} sample-pair and {2 * sample['samples']:,} sample-reference distances are independently recomputed with integer bitsets and checked against the saved counts and exact means.
Text bounds, pairwise text overlaps, embedded PDF font sizes, PDF page dimensions, and PNG dimensions/DPI are checked in `validation.json`.
PNG visual inspection must additionally be recorded in the task log.
No manuscript, bibliography, external implementation, scheduler job, or Git commit is modified by this script.
The final placement and caption insertion in the SI belong to the manuscript editor.
'''

def inspect_outputs(pdf, png):
    """Inspect actual PDF glyph sizes, not just Matplotlib Text properties."""
    from PIL import Image
    from pypdf import PdfReader

    with tempfile.TemporaryDirectory(prefix='neat1-s3-pdf-') as tmp:
        xml = Path(tmp) / 'text.xml'
        command = ['gs', '-q', '-dNOPAUSE', '-dBATCH', '-sDEVICE=txtwrite',
                   '-dTextFormat=0', '-sOutputFile=' + str(xml), str(pdf)]
        subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        spans = ET.fromstring('<document>' + xml.read_text() + '</document>').findall('.//span')
        sizes = [float(s.attrib['size']) for s in spans if s.findall('char')]
        assert sizes and min(sizes) >= FONT_SIZE - .001, 'Embedded PDF glyph below font-size target'
        small = [dict(size_pt=float(s.attrib['size']),
                      text=''.join(c.attrib['c'] for c in s.findall('char')))
                 for s in spans if float(s.attrib['size']) < 6 - .001]
        assert not small, 'Embedded PDF glyph below 6 pt: ' + repr(small)
    pages = PdfReader(pdf).pages
    assert len(pages) == 1
    box = pages[0].mediabox
    dimensions = [float(box.width) / 72 * 25.4, float(box.height) / 72 * 25.4]
    assert np.allclose(dimensions, [WIDTH, HEIGHT], atol=.001, rtol=0)
    with Image.open(png) as im:
        pixels = list(im.size)
        dpi = list(im.info['dpi'])
    assert np.allclose(pixels, [WIDTH * MM * PNG_DPI, HEIGHT * MM * PNG_DPI], atol=1, rtol=0)
    assert np.allclose(dpi, [PNG_DPI, PNG_DPI], atol=.01, rtol=0)
    return dict(pdf_size_mm=dimensions, png_size_px=pixels, png_dpi=dpi,
                embedded_pdf_minimum_pt=min(sizes), embedded_pdf_maximum_pt=max(sizes),
                embedded_pdf_distinct_sizes_pt=sorted(set(sizes)),
                embedded_pdf_spans_checked=len(sizes), embedded_pdf_glyphs_below_6pt=small,
                pdf_font_inspection='Ghostscript txtwrite TextFormat=0; no NODISPLAY flag')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data', type=Path, default=ROOT / 'results/local/e3-neat1-20261016')
    ap.add_argument('--revision', default='e2cecbb4')
    ap.add_argument('--out', type=Path, default=ROOT / 'results/figures/figure_s3_neat1_20261017')
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
    source_provenance = json.loads((args.data / 'provenance.json').read_text())
    if args.revision != 'e2cecbb4':
        assert source_provenance['release'].startswith(args.revision)
        assert not d['supplemental_checkpoints'], 'New figures must use only the new analysis'
    args.out.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(WIDTH * MM, HEIGHT * MM))
    counts = draw(fig, d)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    outside = []
    texts = []
    for t in fig.findobj(matplotlib.text.Text):
        if t.get_visible() and t.get_text():
            assert t.get_fontsize() >= FONT_SIZE, 'Small text: ' + t.get_text()
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
    # Do not crop: native SI width and minimum font size must survive inclusion.
    pdf = args.out / 'figure_s3.pdf'
    png = args.out / 'figure_s3.png'
    fig.savefig(pdf, facecolor='white')
    fig.savefig(png, dpi=PNG_DPI, facecolor='white')
    plt.close(fig)
    validation = inspect_outputs(pdf, png)
    validation.update(text_inside_canvas=True, text_overlap=False,
                      minimum_matplotlib_text_pt=FONT_SIZE,
                      histogram_densities_integrate_to_one=True,
                      all_pairwise_distances_recomputed=True,
                      focus_windows_contain_only_scaffold_free_bases=True)
    (args.out / 'validation.json').write_text(json.dumps(validation, indent=2) + '\n')
    (args.out / 'caption.md').write_text(caption(d, counts))
    (args.out / 'README.md').write_text(source_note(d, counts))
    if args.revision != 'e2cecbb4':
        note = (args.out / 'README.md').read_text()
        note = note.replace('The scientific outputs remain those of LinearCParty `e2cecbb4`, not the ongoing resource rerun.',
                            f'All scientific outputs were recomputed with LinearCParty `{args.revision}`.')
        note = note.replace('e3-neat1-20261016', args.data.name).replace('e2cecbb4-20261016', args.data.name.removeprefix('e3-neat1-'))
        note = note.replace('figure_s3_neat1_20261017', args.out.name)
        note = note.replace('python experiments/scripts/make_figure4.py --beam',
                            f'python experiments/scripts/make_figure4.py --data experiments/{args.data.resolve().relative_to(ROOT)} --revision {args.revision} --beam')
        note = note.replace('# Supplementary Figure S3: NEAT1 (2026-10-17)',
                            f'# Supplementary Figure S3: NEAT1 ({args.revision})')
        note = note.replace('the independent BPP table lacks empty-G b=200, covered by the named checkpoint above',
                            'the new BPP table includes the empty-G b=200 output')
        (args.out / 'README.md').write_text(note)
    with (args.out / 'focus_accessibility.tsv').open('w') as f:
        writer = csv.writer(f, delimiter='\t', lineterminator='\n')
        series = [('G0.7', 200), ('G0.7', 400), ('emptyG', 100), ('emptyG', 200)]
        writer.writerow(['start', 'end'] + [f'{m}_b{b}' for m, b in series])
        for start in range(d['focus'][0], d['focus'][1] + 1):
            writer.writerow([start, start + WINDOW - 1] +
                            [format(d['access_by_beam'][m][b][start - 1], '.17g') for m, b in series])
    metadata = dict(figure_size_mm=[WIDTH, HEIGHT], beam_width=args.beam,
                    panel_beam_widths={p: args.beam for p in ('B', 'C', 'D')},
                    png_dpi=PNG_DPI, ordinary_text_pt=FONT_SIZE, panel_label_pt=9,
                    accessibility_overlay=None,
                    sensitivity_display='numerical main-width means and maximum per-window changes',
                    context_display='seven directly named probability tracks, each on 0--1 scale',
                    context_colors=dict(zip(CONTEXTS, CONTEXT_COLORS)),
                    pair_line_styles=dict(scaffold='solid', new_without_G_crossing='dotted', crossing_G='dashed'),
                    focus_joint_accessibility_by_loop_type=d['focus_byloop'],
                    region_context_max={m: {c: float(d['profile'][m][c][d['region'][0]-1:d['region'][1]].max())
                                            for c in CONTEXTS} for m in ('G0.7', 'emptyG')},
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
                    software_revision=args.revision, coordinate_convention='one-based inclusive',
                    validation=dict(validation, smoothing=False, raw_archives_extracted=False,
                                    rna_cli_called=False))
    (args.out / 'provenance.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps(dict(output=str(args.out), figure_size_mm=[WIDTH, HEIGHT], beam_width=args.beam, arc_counts=counts,
                          focus_mean={k: v['focus_mean'] for k, v in d['diagnostics'].items()},
                          accessibility_validation=d['access_checks'], beam_width_sensitivity=d['stability'],
                          sampling_validation=d['sample_checks']), indent=2))


if __name__ == '__main__':
    main()
