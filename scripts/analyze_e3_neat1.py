#!/usr/bin/env python3
"""Summarize saved E3 outputs; never run RNA software or extract a tar archive.

Run from any directory: python experiments/scripts/analyze_e3_neat1.py
Requires numpy and scipy. Plotting is a separate, rerunnable script.
"""
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
import sys
import tarfile

import numpy as np
from scipy import sparse
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
CONTEXTS = ('Stem', 'Hairpin', 'Bulge', 'Internal', 'Multiloop', 'Exterior', 'PseudoLoop')
WIDTHS = (1, 4, 8, 16)
REGIONS = {'5p_A': (282, 546), '5p_B': (600, 840), '3p': (20761, 21120)}
METHODS = {'NEAT1_1': 'G0.7', 'NEAT1_2': 'G0.7',
           'NEAT1_1_emptyG': 'emptyG', 'NEAT1_2_emptyG': 'emptyG',
           'NEAT1_1_p05': 'G0.5', 'NEAT1_1_p09': 'G0.9'}
PROGRAMS = {'bpp': 'LinearCParty', 'capr': 'LinearPKCapR',
            'raccess': 'LinearPKRaccess', 'mea': 'LinearPKMEA', 'sampling': 'LinearPKSampling'}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.generic):
        x = x.item()
    if isinstance(x, float) and not np.isfinite(x):
        return None
    return x


def write_json(path, value):
    path.write_text(json.dumps(clean(value), indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def write_table(path, rows, fields=None):
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row))
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'wt', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='raise')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: ('NA' if clean(v) is None else clean(v)) for k, v in row.items()})


def pairs_of(db):
    stacks = {c: [] for c in '([{<'}
    closing = dict(zip(')]}>', '([{<'))
    pairs = set()
    for i, c in enumerate(db, 1):
        if c in stacks:
            stacks[c].append(i)
        elif c in closing:
            stack = stacks[closing[c]]
            if not stack:
                raise ValueError('unbalanced structure')
            pairs.add((stack.pop(), i))
        elif c != '.':
            raise ValueError('unknown structure symbol: ' + c)
    if any(stacks.values()):
        raise ValueError('unbalanced structure')
    return pairs


def read_inputs(provenance):
    inputs = {}
    for folder in ('neat1-20261016', 'neat1-variants-20261016'):
        for p in sorted((ROOT / 'data' / folder).glob('*.txt')):
            data = p.read_bytes()
            lines = [l.strip() for l in data.decode().splitlines() if l.strip()]
            name, sequence, scaffold = lines[0][1:], lines[1], lines[2]
            assert name == p.stem and len(sequence) == len(scaffold)
            fixed = pairs_of(scaffold)
            mask = np.zeros(len(sequence), dtype=bool)
            for i, j in fixed:
                mask[[i - 1, j - 1]] = True
            inputs[name] = dict(name=name, iso=name[:7], method=METHODS[name], sequence=sequence,
                                scaffold=scaffold, fixed=fixed, mask=mask, n=len(sequence))
            provenance['inputs'].append({'path': str(p.relative_to(ROOT)), 'sha256': sha(data),
                                         'length': len(sequence), 'scaffold_pairs': len(fixed)})
    for name, v in inputs.items():
        assert v['sequence'] == inputs[v['iso']]['sequence']
    return inputs


def read_runs(inputs, provenance, inventory, run_suffix='e2cecbb4-20261016'):
    """Snapshot complete JSONL lines, then stream only their committed tar members."""
    saved = {}
    for run in sorted((ROOT / 'results/runs').glob('e3-neat1-*-'+run_suffix)):
        kind = run.name.removeprefix('e3-neat1-').removeprefix('variants-').removesuffix('-'+run_suffix)
        if kind not in PROGRAMS:
            continue
        for task in sorted(run.glob('tasks/*/*')):
            plan = json.loads((task / 'planned_jobs.json').read_text())
            jsonl = task / 'run/task_results.jsonl'
            blob = jsonl.read_bytes() if jsonl.exists() else b''
            lines = blob.splitlines(keepends=True)
            complete = b''.join(l for l in lines if l.endswith(b'\n'))
            rows = [json.loads(l) for l in complete.splitlines() if l.strip()]
            by_id = {r['id']: r for r in rows}
            assert len(by_id) == len(rows), 'duplicate checkpoint job'
            provenance['checkpoints'].append({'path': str(jsonl.relative_to(ROOT)),
                'snapshot_bytes': len(blob), 'complete_bytes': len(complete),
                'sha256_complete_lines': sha(complete), 'ignored_partial_tail': len(blob) - len(complete)})
            wanted = {}
            for job in plan:
                r = by_id.get(job['id'])
                beam = int(re.search(r'beam(\d+)', job['mode'])[1])
                v = inputs[job['input']]
                inventory.append(dict(iso=v['iso'], method=v['method'], beam=beam, kind=kind,
                    input=job['input'], status=r['status'] if r else 'missing_checkpoint',
                    job_id=job['id'], checkpoint=str(jsonl.relative_to(ROOT))))
                if not r or r['status'] != 'success' or r.get('warmup'):
                    continue
                command = r['command']
                assert command[command.index('-b') + 1] == str(beam)
                assert command[command.index('--ranking') + 1] == 'outer'
                for flag in ('--max-hairpin', '--max-case3-depth'):
                    assert int(command[command.index(flag) + 1]) == max(beam, 30)
                key = (v['iso'], v['method'], beam)
                group = saved.setdefault(key, {'input': v, 'raw': {}, 'commands': {}})
                assert kind not in group['raw']
                group['commands'][kind] = command
                artifacts = ['stdout.txt'] + (['bpp.tsv'] if kind == 'bpp' else [])
                for basename in artifacts:
                    member = r['directory'] + '/' + basename
                    wanted[member] = (group, kind, basename, r)
            tar = task / 'run/task_outputs.tar'
            if not wanted:
                continue
            with tarfile.open(tar, 'r|') as archive:
                for member in archive:
                    if member.name not in wanted:
                        continue
                    group, kind, basename, row = wanted.pop(member.name)
                    assert member.offset_data + member.size <= row['archive_end']
                    data = archive.extractfile(member).read()
                    expected = row['artifacts'][member.name]
                    assert len(data) == expected['bytes'] and sha(data) == expected['sha256']
                    group['raw'].setdefault(kind, {})[basename] = data.decode()
                    provenance['artifacts'].append({'archive': str(tar.relative_to(ROOT)),
                        'member': member.name, 'sha256': sha(data), 'bytes': len(data),
                        'command': row['command'], 'finished_utc': row.get('finished_utc')})
                    if not wanted:
                        break
            assert not wanted, 'committed artifact missing from archive'
    return saved


def probability(a):
    a = np.asarray(a, dtype=float)
    assert np.all(np.isfinite(a)) and a.min() >= -1e-8 and a.max() <= 1 + 1e-8
    return np.clip(a, 0, 1)


def parse_bpp(text, v):
    pairs, classes = {}, {}
    for line in text.splitlines():
        if not line.strip() or line.startswith(('>', '#')):
            continue
        t = line.split()
        # The saved exports currently contain exactly i, j, probability.
        # Reject schema changes rather than guessing a square/round column.
        assert len(t) == 3, 'BPP schema changed: inspect square/round columns before reuse'
        i, j, p = int(t[0]), int(t[1]), float(t[2])
        assert 1 <= i < j <= v['n'] and (i, j) not in pairs and 0 <= p <= 1 + 1e-8
        pairs[i, j] = p
        if (i, j) not in v['fixed']:
            assert not v['mask'][i - 1] and not v['mask'][j - 1]
        classes[i, j] = 'not_exported'
    assert all(abs(pairs.get(pair, 0) - 1) < 1e-8 for pair in v['fixed'])
    return pairs, classes


def parse_context(text, v):
    rows = {}
    assert [l for l in text.splitlines() if l.startswith('>')] == ['>' + v['name']]
    for line in text.splitlines():
        t = line.split()
        if t and t[0] in CONTEXTS:
            assert t[0] not in rows and len(t) == v['n'] + 1
            rows[t[0]] = probability(t[1:])
    assert set(rows) == set(CONTEXTS)
    return rows


def parse_access(text, v):
    a = {w: np.full(v['n'] - w + 1, np.nan) for w in WIDTHS}
    metadata = {}
    for line in text.splitlines():
        if line.startswith('#'):
            if '=' in line:
                k, value = line[1:].strip().split('=', 1)
                metadata[k] = value
        elif line and not line.startswith('>'):
            start, tokens = line.split('\t')
            start = int(start)  # Raccess is zero-based; output tables below are one-based.
            for token in tokens.strip(';').split(';'):
                w, value = token.split(',')
                w = int(w)
                assert 0 <= start < len(a[w]) and np.isnan(a[w][start])
                a[w][start] = float(value)
    assert metadata['scaffold'] == v['scaffold'] and metadata['output_mode'] == 'probability'
    return {w: probability(values) for w, values in a.items()}


def crossing(pair, fixed):
    i, j = pair
    return any(i < k < j < l or k < i < l < j for k, l in fixed)


def parse_saved(saved, provenance, validation):
    for key, g in saved.items():
        v = g['input']
        if 'bpp' in g['raw']:
            g['pairs'], g['classes'] = parse_bpp(g['raw']['bpp']['bpp.tsv'], v)
            cmd = g['commands']['bpp']
            g['cutoff'] = float(cmd[cmd.index('--bpp-cutoff') + 1])
            total = np.zeros(v['n'])
            cross = np.zeros(v['n'])
            g['cross_pairs'] = {}
            for (i, j), p in g['pairs'].items():
                total[[i - 1, j - 1]] += p
                if (i, j) not in v['fixed'] and crossing((i, j), v['fixed']):
                    g['cross_pairs'][i, j] = p
                    cross[[i - 1, j - 1]] += p
            g['bpp_sum'] = total
            g['cross_mass'] = cross
            assert total.max() <= 1 + 1e-8
        if 'capr' in g['raw']:
            g['context'] = parse_context(g['raw']['capr']['stdout.txt'], v)
        if 'raccess' in g['raw']:
            g['access'] = parse_access(g['raw']['raccess']['stdout.txt'], v)
            g['q'] = g['access'][1]
            g['q_source'] = 'Raccess_w1'
        elif 'context' in g:
            g['q'] = 1 - g['context']['Stem']
            g['q_source'] = '1-CapR_Stem'
        elif 'bpp_sum' in g:
            g['q'] = np.clip(1 - g['bpp_sum'], 0, 1)
            g['q_source'] = '1-truncated_BPP_sum'
        row = dict(iso=key[0], method=key[1], beam=key[2])
        if 'q' in g:
            g['q_before_scaffold_constraints'] = g['q'].copy()
            row['scaffold_q_raw_max'] = float(g['q'][v['mask']].max()) if v['mask'].any() else 0.
            row['scaffold_q_raw_nonzero_n'] = int(np.sum(g['q'][v['mask']] != 0))
        if 'context' in g:
            row['scaffold_stem_raw_max_deviation'] = float(np.max(np.abs(
                g['context']['Stem'][v['mask']] - 1))) if v['mask'].any() else 0.
            err = np.max(np.abs(sum(g['context'].values()) - 1))
            row['context_sum_max_error'] = err
            assert err < 1e-8
            if 'q' in g:
                err = np.max(np.abs(g['q'] - (1 - g['context']['Stem'])))
                row['q_vs_context_max_error'] = err
                assert err < 1e-8
        if 'bpp_sum' in g and 'q' in g:
            omitted = 1 - g['q'] - g['bpp_sum']
            row.update(bpp_omitted_incident_mean=float(omitted.mean()),
                       bpp_omitted_incident_max=float(omitted.max()),
                       bpp_omitted_incident_min=float(omitted.min()))
            assert omitted.min() > -1e-8
        if 'access' in g:
            for w, values in g['access'].items():
                fixed_windows = np.convolve(v['mask'].astype(int), np.ones(w, int), 'valid') > 0
                assert not np.any(values[fixed_windows] > 1e-8)
                upper = np.minimum.reduce([g['q'][j:j + len(values)] for j in range(w)])
                assert np.max(values - upper) < 1e-8
                values[fixed_windows] = 0.
        # Exact conditioning constraints prevent roundoff from inventing ranks
        # among scaffold-paired bases in Spearman/AUROC calculations.
        if 'context' in g:
            for c, values in g['context'].items():
                values[v['mask']] = 1. if c == 'Stem' else 0.
        if 'q' in g:
            g['q'][v['mask']] = 0.
        validation['conditions'].append(row)
    lp = ROOT / 'results/local/e3-neat1-linearpartition-20261016'
    for iso in ('NEAT1_1', 'NEAT1_2'):
        for beam in (100, 200):
            p = lp / f'{iso}_b{beam}_bpp.tsv'
            if not p.exists():
                continue
            data = p.read_bytes()
            v = dict(saved[next(k for k in saved if k[:2] == (iso, 'G0.7'))]['input'])
            v.update(name=iso + '_LP', method='LinearPartition', fixed=set(), mask=np.zeros(v['n'], bool))
            pairs, classes = parse_bpp(data.decode(), v)
            total = np.zeros(v['n'])
            for (i, j), value in pairs.items():
                total[[i - 1, j - 1]] += value
            # LP writes rounded probabilities and omits p < 1e-4.
            # Saved LP output can exceed one beyond the scale of decimal rounding.
            # Keep the raw sums and record the correction rather than renormalizing pairs.
            assert total.max() <= 1.02, 'LP incident mass anomaly exceeds inspected range'
            saved[iso, 'LinearPartition', beam] = dict(input=v, pairs=pairs, classes=classes,
                cutoff=1e-4, bpp_sum=total, cross_mass=np.zeros(v['n']), cross_pairs={},
                q=np.clip(1 - total, 0, 1), q_source='1-truncated_rounded_BPP_sum', raw={}, commands={})
            provenance['inputs'].append(dict(path=str(p.relative_to(ROOT)), sha256=sha(data)))
            validation['conditions'].append(dict(iso=iso, method='LinearPartition', beam=beam,
                exported_BPP_incident_max=float(total.max()), clipped_q_count=int(np.sum(total > 1)),
                incident_excess_max=float(max(0, total.max() - 1)),
                interpretation='incident sums exceed one beyond decimal rounding; origin not resolved here'))


def correlations(x, y):
    if len(x) < 3 or np.ptp(x) < 1e-10 or np.ptp(y) < 1e-10:
        return np.nan, np.nan
    return float(np.corrcoef(rankdata(x), rankdata(y))[0, 1]), float(np.corrcoef(x, y)[0, 1])


def auc_ap(labels, scores):
    """Tie-aware Mann-Whitney AUROC and step-integrated precision-recall (AP)."""
    labels = np.asarray(labels, bool)
    npos, nneg = int(labels.sum()), int((~labels).sum())
    if not npos or not nneg:
        return np.nan, np.nan
    if np.ptp(scores) < 1e-10:
        return .5, npos / len(labels)
    auc = (rankdata(scores)[labels].sum() - npos * (npos + 1) / 2) / (npos * nneg)
    order = np.argsort(-scores, kind='stable')
    ends = np.r_[np.flatnonzero(np.diff(scores[order]) != 0), len(order) - 1]
    tp = np.cumsum(labels[order])[ends]
    precision = tp / (ends + 1)
    ap = np.sum(np.diff(np.r_[0, tp]) / npos * precision)
    return float(auc), float(ap)


def features(g):
    f = {'q': g['q']}
    if 'bpp_sum' in g:
        f['bpp_sum_exported'] = g['bpp_sum']
    if g['q_source'] == 'Raccess_w1':
        f['bpp_from_q'] = 1 - g['q']
    f.update(g.get('context', {}))
    return f


def strata(v, ref):
    return {'all': np.ones(v['n'], bool), 'ownG_paired': v['mask'], 'ownG_free': ~v['mask'],
            'G07_paired': ref['mask'], 'G07_free': ~ref['mask']}


def regression(x, y):
    active = np.ptp(x, axis=0) > 1e-10
    design = np.column_stack([np.ones(len(y)), x[:, active]])
    beta, _, rank, _ = np.linalg.lstsq(design, y, rcond=None)
    pred = design @ beta
    denom = np.sum((y - y.mean()) ** 2)
    r2 = 1 - np.sum((y - pred) ** 2) / denom if denom else np.nan
    cv_pred = np.full(len(y), np.nan)
    # Contiguous held-out blocks preserve a modest separation from neighboring positions.
    for test in np.array_split(np.arange(len(y)), 5):
        train = np.ones(len(y), bool)
        train[test] = False
        cv_pred[test] = design[test] @ np.linalg.lstsq(design[train], y[train], rcond=None)[0]
    cv_r2 = 1 - np.sum((y - cv_pred) ** 2) / denom if denom else np.nan
    coefficients = np.full(x.shape[1], np.nan)
    if rank == design.shape[1]:
        coefficients[active] = beta[1:]
    return dict(n=len(y), rank=int(rank), active_terms=int(active.sum()),
                condition_number=float(np.linalg.cond(design)),
                r2=float(r2), blocked5_r2=float(cv_r2), intercept=float(beta[0]),
                status='full_rank_active_terms' if rank == design.shape[1] else 'rank_deficient'), coefficients


def analyze_shape(saved, inputs, shape, threshold):
    metrics, bins, offsets, regressions = [], [], [], []
    pos = np.flatnonzero(np.isfinite(shape))
    for key, g in sorted(saved.items()):
        if key[0] != 'NEAT1_1' or 'q' not in g:
            continue
        common = dict(iso=key[0], method=key[1], beam=key[2], q_source=g['q_source'])
        v = g['input']
        for stratum, mask in strata(v, inputs[key[0]]).items():
            use = pos[mask[pos]]
            y = shape[use]
            labels = y >= threshold
            for name, values in features(g).items():
                x = values[use]
                rho, pearson = correlations(x, y)
                sign = -1 if name in ('Stem', 'bpp_sum_exported', 'bpp_from_q') else 1
                auc_raw, ap_raw = auc_ap(labels, x)
                auc, ap = auc_ap(labels, sign * x)
                metrics.append(dict(**common, stratum=stratum, feature=name, n=len(use),
                    positives=int(labels.sum()), prevalence=float(labels.mean()) if len(y) else np.nan,
                    shape_threshold=threshold, spearman=rho, pearson=pearson,
                    high_shape_score='1-p' if sign < 0 else 'p', auroc=auc, auprc_ap=ap,
                    auroc_raw_probability=auc_raw, auprc_raw_probability=ap_raw,
                    status='constant_predictor' if len(x) and np.ptp(x) < 1e-10 else ('empty' if not len(x) else 'ok')))
                ids = np.minimum((np.clip(x, 0, 1) * 10).astype(int), 9)
                for b in range(10):
                    select = ids == b
                    yy = y[select]
                    bins.append(dict(**common, stratum=stratum, feature=name, bin_left=b / 10,
                        bin_right=(b + 1) / 10, n=len(yy),
                        probability_mean=float(x[select].mean()) if len(yy) else np.nan,
                        shape_mean=float(yy.mean()) if len(yy) else np.nan,
                        shape_sd=float(yy.std(ddof=1)) if len(yy) > 1 else np.nan))
            if len(use) > 20 and 'context' in g:
                names = list(CONTEXTS[1:])  # Stem is the omitted compositional reference.
                x = np.column_stack([g['context'][c][use] for c in names])
                result, coef = regression(x, y)
                q_result, _ = regression(g['q'][use, None], y)
                for term, value in zip(['Intercept'] + names, [result['intercept']] + list(coef)):
                    regressions.append(dict(**common, stratum=stratum, dropped_reference='Stem',
                        term=term, coefficient=value, **result,
                        q_only_r2=q_result['r2'], q_only_blocked5_r2=q_result['blocked5_r2']))
        # The same SHAPE positions are used for every shift; +d maps SHAPE p to RNA p+d.
        fixed_support = pos[(pos >= 5) & (pos < min(len(shape), 3595) - 5)]
        for shift in range(-5, 6):
            for name, values in features(g).items():
                rho, r = correlations(values[fixed_support + shift], shape[fixed_support])
                offsets.append(dict(**common, feature=name, shift_nt=shift, n=len(fixed_support),
                                    spearman=rho, pearson=r))
    return metrics, bins, regressions, offsets


def desc(a):
    a = np.asarray(a)
    return dict(n=len(a), mean=float(a.mean()) if len(a) else np.nan,
        median=float(np.median(a)) if len(a) else np.nan,
        p05=float(np.quantile(a, .05)) if len(a) else np.nan,
        p95=float(np.quantile(a, .95)) if len(a) else np.nan,
        min=float(a.min()) if len(a) else np.nan, max=float(a.max()) if len(a) else np.nan)


def access_analysis(saved):
    summary, differences, tiles = [], [], []
    for key, g in sorted(saved.items()):
        regions = {'all': (1, g['input']['n'])}
        if key[0] == 'NEAT1_2':
            regions.update(REGIONS)
        available = g.get('access', {1: g['q']} if 'q' in g else {})
        for w, a in available.items():
            for region, (lo, hi) in regions.items():
                values = a[lo - 1:hi - w + 1]  # windows wholly contained in region
                summary.append(dict(iso=key[0], method=key[1], beam=key[2], region=region,
                    start=lo, end=hi, width=w,
                    source='Raccess' if 'access' in g else g['q_source'], **desc(values)))
    for key, g in sorted(saved.items()):
        if 'access' not in g:
            continue
        iso, method, beam = key
        comparison_keys = []
        if method != 'emptyG':
            comparison_keys.append((iso, 'emptyG', beam))
        if method in ('G0.5', 'G0.9'):
            comparison_keys.append((iso, 'G0.7', beam))
        previous = {200: 100, 400: 200}.get(beam)
        if previous:
            comparison_keys.append((iso, method, previous))
        for refkey in comparison_keys:
            ref = saved.get(refkey)
            if not ref or 'access' not in ref:
                continue
            for w in WIDTHS:
                a, b = g['access'][w], ref['access'][w]
                delta = a - b
                unionmask = g['input']['mask'] | ref['input']['mask']
                common_free = np.convolve(unionmask.astype(int), np.ones(w, int), 'valid') == 0
                for label, mask in (('all', np.ones(len(a), bool)), ('bothG_free', common_free)):
                    d = delta[mask]
                    rho, _ = correlations(a[mask], b[mask])
                    differences.append(dict(iso=iso, method=method, beam=beam, reference=refkey[1],
                        reference_beam=refkey[2], width=w, stratum=label, **desc(d),
                        mean_abs_difference=float(np.mean(np.abs(d))) if len(d) else np.nan,
                        n_abs_gt_0_2=int(np.sum(np.abs(d) > .2)), spearman=rho))
                # Fixed 50-nt start-coordinate tiles, including a short terminal tile.
                for start in range(0, len(a), 50):
                    d = delta[start:start + 50]
                    tiles.append(dict(iso=iso, method=method, beam=beam, reference=refkey[1],
                        reference_beam=refkey[2], width=w, start=start + 1,
                        end=min(start + 50, len(a)), n=len(d), mean_delta=float(d.mean()),
                        mean_abs_delta=float(np.abs(d).mean()),
                        bothG_free_n=int(common_free[start:start + 50].sum()),
                        candidate_mean=float(a[start:start + 50].mean()), reference_mean=float(b[start:start + 50].mean())))
    return summary, differences, tiles


def long_range(saved):
    rows = []
    for method, beams in [('LinearPartition', (100, 200)), ('emptyG', (100, 200)), ('G0.7', (100, 200, 400))]:
        for beam in beams:
            g = saved.get(('NEAT1_2', method, beam))
            for region in ('5p_A', '5p_B', '5p_union'):
                ranges = [REGIONS[region]] if region != '5p_union' else [REGIONS['5p_A'], REGIONS['5p_B']]
                pairs = {} if not g else {p: x for p, x in g.get('pairs', {}).items()
                    if any(lo <= p[0] <= hi for lo, hi in ranges) and REGIONS['3p'][0] <= p[1] <= REGIONS['3p'][1]}
                fixed = set() if not g else {p for p in g['input']['fixed']
                    if any(lo <= p[0] <= hi for lo, hi in ranges) and REGIONS['3p'][0] <= p[1] <= REGIONS['3p'][1]}
                assert not fixed
                present = g is not None and 'pairs' in g
                maxpair = max(pairs, key=pairs.get) if pairs else (None, None)
                cutoff = g.get('cutoff') if g else (1e-4 if method == 'LinearPartition' else .001)
                possible = sum(hi - lo + 1 for lo, hi in ranges) * (REGIONS['3p'][1] - REGIONS['3p'][0] + 1)
                standardized = {p: x for p, x in pairs.items() if x >= .001}
                rows.append(dict(iso='NEAT1_2', method=method, beam=beam, region=region,
                    partner='3p', status='available' if present else 'missing_BPP', cutoff=cutoff,
                    retained_pairs=len(pairs) if present else None,
                    probability_sum_retained=sum(pairs.values()) if present else None,
                    probability_max_retained=max(pairs.values(), default=0) if present else None,
                    max_i=maxpair[0], max_j=maxpair[1],
                    common_cutoff=.001, common_cutoff_pairs=len(standardized) if present else None,
                    common_cutoff_sum=sum(standardized.values()) if present else None,
                    common_cutoff_max=max(standardized.values(), default=0) if present else None,
                    omitted_sum_bound_loose=(possible - len(pairs)) * cutoff if present else None,
                    interpretation='retained sum is a lower bound; absent pairs are below cutoff, not measured zero'))
    return rows


def crossing_analysis(saved, shape):
    summary, positions = [], []
    for key, g in sorted(saved.items()):
        if key[1] != 'G0.7' or 'cross_pairs' not in g:
            continue
        pairs = g['cross_pairs']
        summary.append(dict(iso=key[0], method=key[1], beam=key[2], cutoff=g['cutoff'],
            n_crossing_pairs=len(pairs), sum_pair_probability=sum(pairs.values()),
            max_pair_probability=max(pairs.values(), default=0),
            n_pairs_ge_0_5=sum(p >= .5 for p in pairs.values()),
            n_positions_mass_ge_0_5=int(np.sum(g['cross_mass'] >= .5))))
        positive = np.flatnonzero(g['cross_mass'] > 0)
        ranked = sorted(positive, key=lambda p: (-g['cross_mass'][p], p))[:30]
        for rank, p in enumerate(ranked, 1):
            partners = [(j if i == p + 1 else i, value) for (i, j), value in pairs.items() if p + 1 in (i, j)]
            partner, prob = max(partners, key=lambda t: t[1])
            positions.append(dict(iso=key[0], method=key[1], beam=key[2], rank=rank, position=p + 1,
                crossing_incident_mass=g['cross_mass'][p], max_partner=partner, max_pair_probability=prob,
                q=g['q'][p], shape=shape[p] if key[0] == 'NEAT1_1' else None))
        if key[0] == 'NEAT1_1':
            for label, mask in [('mass_ge_0_5', g['cross_mass'] >= .5),
                                ('mass_lt_0_5_Gfree', (g['cross_mass'] < .5) & ~g['input']['mask'])]:
                use = mask & np.isfinite(shape)
                summary.append(dict(iso=key[0], method=key[1], beam=key[2], group=label,
                    shape_n=int(use.sum()), shape_mean=float(shape[use].mean()) if use.any() else None,
                    shape_median=float(np.median(shape[use])) if use.any() else None))
    return summary, positions


def structure_from_stdout(text, n, sampled=False):
    found = []
    for line in text.splitlines():
        tokens = line.split()
        if sampled:
            tokens = tokens[1:] if tokens and tokens[0].isdigit() else []
        for t in tokens[:1]:
            if len(t) == n and set(t) <= set('.()[]{}<>'):
                found.append(t)
    return found


def jaccard(a, b):
    union = len(a | b)
    return 1 - len(a & b) / union if union else 0.


def sampling_analysis(saved):
    summary, histograms, distances, representatives, diagnostics = [], [], [], [], []
    for key, g in sorted(saved.items()):
        if key[0] != 'NEAT1_1' or 'sampling' not in g['raw']:
            continue
        v = g['input']
        dbs = structure_from_stdout(g['raw']['sampling']['stdout.txt'], v['n'], sampled=True)
        assert len(dbs) == 1000
        structures = [pairs_of(db) for db in dbs]
        assert all(v['fixed'] <= s for s in structures)
        mfe = structure_from_stdout(g['raw']['bpp']['stdout.txt'], v['n'])
        mea = structure_from_stdout(g['raw']['mea']['stdout.txt'], v['n'])
        assert len(mfe) == len(mea) == 2 and mfe[0] == mea[0] == v['scaffold']
        refs = {'MFE': pairs_of(mfe[-1]), 'MEA': pairs_of(mea[-1])}
        cmd = g['commands']['sampling']
        seed = cmd[cmd.index('--seed') + 1] if '--seed' in cmd else 'not_recorded'
        for space in ('completed', 'extension'):
            ss = [s if space == 'completed' else s - v['fixed'] for s in structures]
            union = sorted(set().union(*ss))
            index = {p: i for i, p in enumerate(union)}
            rr = np.repeat(np.arange(len(ss)), [len(s) for s in ss])
            cc = np.array([index[p] for s in ss for p in sorted(s)])
            mat = sparse.csr_matrix((np.ones(len(cc), np.int32), (rr, cc)), shape=(len(ss), len(union)))
            intersection = (mat @ mat.T).toarray()
            counts = np.array([len(s) for s in ss])
            unions = counts[:, None] + counts[None, :] - intersection
            d = 1 - np.divide(intersection, unions, out=np.ones_like(unions, dtype=float), where=unions != 0)
            values = d[np.triu_indices(len(ss), 1)]
            # Independent set arithmetic checks sparse intersection and pair encoding.
            for i, j in [(0, 1), (0, 999), (123, 456)]:
                assert abs(d[i, j] - jaccard(ss[i], ss[j])) < 1e-12
            common = dict(iso=key[0], beam=key[2], space=space, seed=seed)
            summary.append(dict(**common, comparison='sample_sample', unique_structures=len(set(dbs)),
                                samples=len(dbs), **desc(values)))
            count, edges = np.histogram(values, bins=np.linspace(0, 1, 101))
            for i, n in enumerate(count):
                histograms.append(dict(**common, comparison='sample_sample', left=edges[i], right=edges[i + 1], count=n))
            for refname, ref in refs.items():
                target = ref if space == 'completed' else ref - v['fixed']
                values = np.array([jaccard(s, target) for s in ss])
                summary.append(dict(**common, comparison='sample_' + refname, **desc(values)))
                for i, value in enumerate(values, 1):
                    distances.append(dict(**common, sample=i, reference=refname, jaccard_distance=value))
                count, edges = np.histogram(values, bins=np.linspace(0, 1, 101))
                for i, n in enumerate(count):
                    histograms.append(dict(**common, comparison='sample_' + refname, left=edges[i], right=edges[i + 1], count=n))
            a, b = refs['MFE'], refs['MEA']
            if space == 'extension':
                a, b = a - v['fixed'], b - v['fixed']
            representatives.append(dict(**common, comparison='MFE_MEA', jaccard_distance=jaccard(a, b),
                mfe_pairs=len(a), mea_pairs=len(b)))
        empirical = Counter(p for s in structures for p in s)
        errors = [abs(empirical[p] / len(structures) - prob) for p, prob in g['pairs'].items()]
        diagnostics.append(dict(beam=key[2], comparison='sampling_frequency_vs_retained_BPP',
            n_pairs=len(errors), mean_abs_error=float(np.mean(errors)), max_abs_error=max(errors),
            interpretation='descriptive Monte Carlo comparison; pairs below BPP cutoff are omitted'))
    return summary, histograms, distances, representatives, diagnostics


def profile_rows(saved, shape):
    for key, g in sorted(saved.items()):
        v = g['input']
        if 'q' not in g:
            continue
        missing = np.full(v['n'], np.nan)
        bpp = g.get('bpp_sum', missing)
        cross = g.get('cross_mass', missing)
        contexts = {c: g.get('context', {}).get(c, missing) for c in CONTEXTS}
        for i in range(v['n']):
            yield dict(iso=key[0], method=key[1], beam=key[2], position=i + 1,
                ownG_paired=int(v['mask'][i]), q=g['q'][i], q_source=g['q_source'],
                bpp_sum_exported=bpp[i],
                crossing_incident_mass=cross[i],
                shape=shape[i] if key[0] == 'NEAT1_1' else None,
                **{c: contexts[c][i] for c in CONTEXTS})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'results/local/e3-neat1-20261016')
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    provenance = dict(created_utc=datetime.now(timezone.utc).isoformat(),
        argv=sys.argv, script_sha256=sha(Path(__file__).read_bytes()), inputs=[], checkpoints=[], artifacts=[])
    validation = {'conditions': [], 'checks': [
        'Consumed tar members match committed byte lengths and SHA256; no extraction',
        'Input variant sequences identical to their base isoform',
        'Raccess zero-based starts converted to one-based; all requested windows present',
        'Contexts normalized; Raccess w1 agrees with 1-Stem',
        'Windows touching scaffold bases have zero probability; joint probability <= each marginal',
        'Scaffold pairs retained with probability one; new pairs use only free bases',
        'No scaffold pairs join the specified long-range regions',
        'Sampling includes scaffold; sparse Jaccard checked against independent set arithmetic']}
    inputs = read_inputs(provenance)
    provenance['metadata'] = []
    for path in [ROOT / 'data/neat1-20261016/manifest.json',
                 ROOT / 'data/neat1-variants-20261016/manifest.json',
                 ROOT / 'data/neat1-20261016/shape_lin2018/manifest.json',
                 ROOT / 'manifests/software-e3-e2cecbb4-20261016.json']:
        provenance['metadata'].append(dict(path=str(path.relative_to(ROOT)), sha256=sha(path.read_bytes())))
    inventory = []
    saved = read_runs(inputs, provenance, inventory)
    print('Read saved outputs:', len(inventory), 'planned jobs', flush=True)
    parse_saved(saved, provenance, validation)
    print('Probability and coordinate checks passed; analyzing SHAPE and sampling', flush=True)
    sp = ROOT / 'data/neat1-20261016/shape_lin2018/hNEAT1_fl_shape.txt'
    shape_data = np.loadtxt(sp)
    assert np.array_equal(shape_data[:, 0], np.arange(1, 3596))
    shape = np.full(inputs['NEAT1_1']['n'], np.nan)
    shape[:3595] = shape_data[:, 1]  # Keep negative observations; no tail winsorization.
    provenance['inputs'].append(dict(path=str(sp.relative_to(ROOT)), sha256=sha(sp.read_bytes())))
    threshold = float(np.quantile(shape_data[:, 1], .75))
    validation['scaffold_constraint_sensitivity'] = []
    for key, g in sorted(saved.items()):
        if key[0] != 'NEAT1_1' or key[1] == 'LinearPartition':
            continue
        for treatment, x in [('raw_Raccess', g['q_before_scaffold_constraints'][:3595]),
                             ('exact_scaffold_constraints', g['q'][:3595])]:
            rho, r = correlations(x, shape[:3595])
            auc, ap = auc_ap(shape[:3595] >= threshold, x)
            validation['scaffold_constraint_sensitivity'].append(dict(method=key[1], beam=key[2],
                treatment=treatment, spearman=rho, pearson=r, auroc=auc, auprc_ap=ap))
    validation['linearpartition_q_sensitivity'] = []
    for beam in (100, 200):
        g = saved['NEAT1_1', 'LinearPartition', beam]
        x, y = 1 - g['bpp_sum'][:3595], shape[:3595]
        for treatment, mask in [('raw_unclipped', np.ones(len(y), bool)),
                                ('drop_incident_sum_gt_1', x >= 0)]:
            rho, r = correlations(x[mask], y[mask])
            auc, ap = auc_ap(y[mask] >= threshold, x[mask])
            validation['linearpartition_q_sensitivity'].append(dict(beam=beam, treatment=treatment,
                n=int(mask.sum()), spearman=rho, pearson=r, auroc=auc, auprc_ap=ap))
    metrics, bins, regressions, offsets = analyze_shape(saved, inputs, shape, threshold)
    access, differences, tiles = access_analysis(saved)
    cross, crosspos = crossing_analysis(saved, shape)
    sampling, hist, distances, reps, samp_checks = sampling_analysis(saved)
    validation['sampling'] = samp_checks
    validation['status'] = 'all_checks_passed'
    write_table(out / 'inventory.tsv', inventory)
    write_table(out / 'shape_metrics.tsv', metrics)
    write_table(out / 'shape_probability_bins.tsv', bins)
    write_table(out / 'shape_regression.tsv', regressions)
    write_table(out / 'shape_offsets.tsv', offsets)
    write_table(out / 'accessibility_summary.tsv', access)
    write_table(out / 'accessibility_differences.tsv', differences)
    write_table(out / 'accessibility_regions.tsv', tiles)
    write_table(out / 'crossing_summary.tsv', cross)
    write_table(out / 'crossing_positions.tsv', crosspos,
        ['iso', 'method', 'beam', 'rank', 'position', 'crossing_incident_mass', 'max_partner', 'max_pair_probability', 'q', 'shape'])
    write_table(out / 'long_range.tsv', long_range(saved))
    write_table(out / 'sampling_summary.tsv', sampling)
    write_table(out / 'sampling_histograms.tsv', hist)
    write_table(out / 'sampling_reference_distances.tsv.gz', distances)
    write_table(out / 'representative_distances.tsv', reps)
    write_table(out / 'position_profiles.tsv.gz', profile_rows(saved, shape),
        ['iso', 'method', 'beam', 'position', 'ownG_paired', 'q', 'q_source', 'bpp_sum_exported',
         'crossing_incident_mass', 'shape'] + list(CONTEXTS))
    write_table(out / 'accessibility.tsv.gz', (dict(iso=k[0], method=k[1], beam=k[2],
        start=i + 1, width=w, accessibility=value)
        for k, g in sorted(saved.items()) for w, a in g.get('access', {}).items() for i, value in enumerate(a)),
        ['iso', 'method', 'beam', 'start', 'width', 'accessibility'])
    write_table(out / 'bpp_pairs.tsv.gz', (dict(iso=k[0], method=k[1], beam=k[2], i=i, j=j,
        probability=p, is_scaffold=int((i, j) in g['input']['fixed']),
        crosses_G=int((i, j) in g.get('cross_pairs', {})), output_pair_class=g['classes'][i, j],
        shape_i=shape[i - 1] if k[0] == 'NEAT1_1' else None,
        shape_j=shape[j - 1] if k[0] == 'NEAT1_1' else None, cutoff=g['cutoff'])
        for k, g in sorted(saved.items()) for (i, j), p in sorted(g.get('pairs', {}).items())),
        ['iso', 'method', 'beam', 'i', 'j', 'probability', 'is_scaffold', 'crosses_G',
         'output_pair_class', 'shape_i', 'shape_j', 'cutoff'])
    summary = dict(shape_n=3595, shape_negative_n=int(np.sum(shape_data[:, 1] < 0)),
        high_shape_threshold=threshold, high_shape_n=int(np.sum(shape_data[:, 1] >= threshold)),
        inventory_status=dict(Counter(r['status'] for r in inventory)),
        condition_count=len(saved), validation='all_checks_passed',
        note='Coordinates assume SHAPE position 1 = transcript nt 1; positive shift maps SHAPE p to RNA p+shift.')
    write_json(out / 'summary.json', summary)
    write_json(out / 'validation.json', validation)
    write_json(out / 'provenance.json', provenance)
    (out / 'README.md').write_text('''# E3 saved-output analysis (2026-10-16)

Regenerate from either repository directory:
`python experiments/scripts/analyze_e3_neat1.py` (paper repository), or
`python scripts/analyze_e3_neat1.py` (experiment repository).
Then run the corresponding `plot_e3_neat1.py` command.
No RNA CLI, scheduler submission, or tar extraction is performed.
Re-running discovers newly checkpointed jobs in the same eight fixed run directories.
Only newline-terminated JSONL records and SHA256-verified committed tar members are used.
Missing conditions remain explicit in inventory.tsv and long_range.tsv.
Unscheduled variant b400 conditions are not imputed.

All output positions and window starts are one-based and inclusive.
Raccess raw starts are zero-based and are converted once.
Window summaries include only windows wholly inside the named interval.
SHAPE values, including 36 negative values, are used without clipping, winsorization or end trimming.
High SHAPE means at least the global 75th percentile, computed once on the 3,595 observations.
The label threshold is unchanged across conditions and strata.
ownG_* uses each condition's scaffold; G07_* uses the common G=0.7 scaffold for fair subsets.
LP and emptyG have no ownG_paired positions.
Constant predictors have undefined correlations and tie-based AUROC=0.5/AP=prevalence.
Probability bins are [0,.1), ...,[.9,1], with raw SHAPE mean and sample SD; empty bins are NA.
Paired-probability predictors (Stem/bpp) use 1-p for high-SHAPE ranking; raw-score AUROC/AP are also saved.
AUPRC is average precision (step integral), not trapezoidal PR area.

q uses Raccess length 1 when available, then 1-CapR Stem, then 1-exported BPP sum.
After recording raw numerical residuals, scaffold bases are assigned exact q=0,
Stem=1 and other contexts=0 to prevent artificial ranks from roundoff.
All windows touching G are assigned exact accessibility zero after raw checks.
Raw-versus-canonical SHAPE metrics are in validation.json because tiny residuals
can substantially change ranks among the many scaffold-paired positions.
LinearCParty BPP cutoff is .001; LP cutoff is .0001 and probabilities are rounded in its export.
bpp_sum_exported sums the retained incident pair probabilities, not an exact paired marginal.
bpp_from_q uses the full Raccess marginal and is reported separately.
Some LP exported incident sums exceed one beyond the scale of decimal rounding.
Their origin is unresolved here; LP q is clipped to [0,1] and raw sums are retained.
Counts, excess sizes, and unclipped/drop-excess SHAPE sensitivity are in validation.json.
No longer-window accessibility or context probabilities are inferred from LP BPP.
Long-range sums are retained-mass lower bounds: a zero means no exported pair at/above cutoff.
The loose omission bound is number of absent position pairs times cutoff and is not an error estimate.
No conclusion about exact zero long-range pairing is justified.
long_range.tsv also applies a common probability cutoff .001 to both tools for comparison.

BPP exports have only i,j,p, with no square/round column.
output_pair_class=not_exported preserves that absence; crosses_G is a separate geometric test.
Crossing probabilities are sums of incident retained pairs crossing at least one G pair.
Given mutually exclusive incident pairs, this sum is a lower bound on a base participating in such a pair.
Comparing G0.7 with emptyG changes both the conditioning scaffold and allowed ensemble.
It does not isolate the causal effect of allowing pseudoknots under a fixed G.

Regression is unweighted OLS with intercept and six contexts, dropping Stem as compositional reference.
Predictors constant within a stratum are additionally excluded and their coefficients are NA.
Rank and active-term count are explicit; rank-deficient coefficients are NA.
R2 is descriptive in-sample fit; blocked5_r2 holds out each of five contiguous position blocks.
No independence-based p-values or confidence intervals are reported for neighboring bases.
Coordinate sensitivity uses a common set of SHAPE positions 6..3590 for offsets -5..5.
Positive offset means matching SHAPE position p to RNA nucleotide p+offset.

Sampling uses all 1,000 structures, all 499,500 distinct unordered sample pairs,
and exact pair-set Jaccard distance 1-|A intersect B|/|A union B|.
Both completed pair sets and extension-only pair sets (remove G) are reported.
Distances have no pseudoknot bracket-type dependence.
Histograms have 100 equal-width bins; all sample-to-MFE/MEA distances are in one gzip TSV.
MFE is the structure in saved BPP stdout; MEA gamma=1 is in saved MEA stdout.
Finite-beam MFE can use a different retained search set from the partition-function ensemble.

position_profiles.tsv.gz, accessibility.tsv.gz and bpp_pairs.tsv.gz pool all positions in three files.
provenance.json records checkpoint snapshots, consumed-member hashes, inputs, commands and script hash.
validation.json records numerical checks and the BPP mass omitted relative to full marginals.
Figures and bilingual draft captions are in results/figures/draft_e3_20261016/.
''')
    print(json.dumps(summary, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
