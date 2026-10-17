#!/usr/bin/env python3
"""Audit the fixed NEAT1 rerun and compare saved outputs without extracting tar files."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re
import tarfile

import numpy as np

from analyze_e3_neat1 import ROOT, WIDTHS, read_inputs, parse_access, parse_bpp, parse_context, write_json, write_table

OLD = 'e2cecbb4-20261016'
NEW = 'd02d261c-20261017'
KINDS = ('bpp', 'capr', 'raccess', 'mea', 'sampling', 'raccess-byloop')
LOOPS = dict(E='Exterior', H='Hairpin', B='Bulge', I='Internal', M='Multiloop', P='PseudoLoop')


def sha(data):
    return hashlib.sha256(data).hexdigest()


def normalized_command(command):
    command = list(command)
    command[0] = Path(command[0]).name
    for flag in ('-i', '--bpp'):
        if flag in command:
            command[command.index(flag) + 1] = Path(command[command.index(flag) + 1]).name
    return command


def load(version, inputs, provenance):
    saved = {}
    for kind in KINDS:
        if version == OLD and kind == 'raccess-byloop':
            continue
        run = ROOT / 'results/runs' / f'e3-neat1-{kind}-{version}'
        manifest = json.loads((run / 'run.json').read_text())
        assert manifest['software_before']['commit'].startswith(version.split('-')[0])
        for v in manifest['inputs']:
            assert sha((ROOT / 'data/neat1-20261016' / (v['name'] + '.txt')).read_bytes()) == v['sha256']
        for task in sorted(run.glob('tasks/*/*')):
            meta = json.loads((task / 'run/run.json').read_text())
            assert meta['status'] == 'success'
            assert meta['software_before'] == meta['software_after']
            assert meta['binaries_before'] == meta['binaries_after']
            assert meta['software_lock_valid_after_run'] is True
            assert meta['software_changed_during_run'] is False
            assert meta['binaries_changed_during_run'] is False
            blob = (task / 'run/task_results.jsonl').read_bytes()
            assert blob.endswith(b'\n'), 'incomplete checkpoint'
            rows = [json.loads(line) for line in blob.splitlines()]
            plan = json.loads((task / 'planned_jobs.json').read_text())
            assert len(rows) == len(plan) == len({r['id'] for r in rows})
            assert {r['id'] for r in rows} == {r['id'] for r in plan}
            wanted = {}
            for row in rows:
                assert row['status'] == 'success' and row['returncode'] == 0 and not row['error']
                assert not row['warmup'] and row['repetition'] == 0
                beam = int(re.fullmatch(r'beam(\d+)_outer', row['mode'])[1])
                command = row['command']
                assert command[command.index('-b') + 1] == str(beam)
                for flag in ('--max-hairpin', '--max-case3-depth'):
                    assert int(command[command.index(flag) + 1]) == max(beam, 30)
                assert command[command.index('--ranking') + 1] == 'outer'
                key = row['input'], beam, kind
                assert key not in saved
                saved[key] = dict(row=row, raw={})
                for member, info in row['artifacts'].items():
                    assert member not in wanted
                    wanted[member] = key, row, info
            with tarfile.open(task / 'run/task_outputs.tar', 'r|') as tar:
                for member in tar:
                    if member.name not in wanted:
                        continue
                    key, row, info = wanted.pop(member.name)
                    assert member.offset_data + member.size <= row['archive_end']
                    data = tar.extractfile(member).read()
                    assert len(data) == info['bytes'] and sha(data) == info['sha256']
                    basename = Path(member.name).name
                    if basename in ('stdout.txt', 'bpp.tsv'):
                        saved[key]['raw'][basename] = data.decode()
                    elif basename == 'time.tsv':
                        wall, rss = data.decode().split()
                        assert float(wall) == row['time_wall_seconds'] and int(rss) == row['peak_rss_kib']
                    provenance['artifacts'].append(dict(archive=str((task / 'run/task_outputs.tar').relative_to(ROOT)),
                        member=member.name, **info))
            assert not wanted
            provenance['checkpoints'].append(dict(path=str((task / 'run/task_results.jsonl').relative_to(ROOT)),
                sha256=sha(blob), jobs=len(rows), status=meta['status'], workers=meta['workers'],
                host=meta['host'], software_before=meta['software_before'],
                energy_parameter_sha256=meta['energy_parameter_sha256'],
                software_lock_sha256=meta['software_lock_sha256']))
    expected = {(iso, b, k) for iso in ('NEAT1_1', 'NEAT1_2') for b in (100, 200, 400)
                for k in KINDS if not (version == OLD and k == 'raccess-byloop')
                and (k != 'sampling' or (iso, b) in (('NEAT1_1', 100), ('NEAT1_1', 200), ('NEAT1_2', 100)))}
    assert set(saved) == expected
    return saved


def parse_byloop(text, n):
    total = {w: np.full(n-w+1, np.nan) for w in WIDTHS}
    components = {w: {c: np.full(n-w+1, np.nan) for c in LOOPS} for w in WIDTHS}
    assert '# byloop_value=joint_probability' in text
    for line in text.splitlines():
        if not line or line.startswith(('#', '>')):
            continue
        pos, totals, values = line.split('\t')
        pos = int(pos)
        assert values.startswith('byloop:')
        for item in totals.strip(';').split(';'):
            w, p = item.split(',')
            w = int(w)
            assert np.isnan(total[w][pos])
            total[w][pos] = float(p)
        for item in values[len('byloop:'):].strip(';').split(';'):
            w, *entries = item.split(',')
            w = int(w)
            assert {entry.split('=')[0] for entry in entries} == set(LOOPS)
            for entry in entries:
                c, p = entry.split('=')
                assert np.isnan(components[w][c][pos])
                components[w][c][pos] = float(p)
    for w in WIDTHS:
        for a in (total[w], *components[w].values()):
            assert np.all(np.isfinite(a)) and a.min() >= 0 and a.max() <= 1+1e-8
    return total, components


def structures(text, n):
    return [token for line in text.splitlines() for token in line.split()
            if len(token) == n and set(token) <= set('.()[]{}<>')]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT / 'results/local/e3-neat1-d02d261c-20261017')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    provenance = dict(inputs=[], checkpoints=[], artifacts=[], script_sha256=sha(Path(__file__).read_bytes()))
    inputs = read_inputs(provenance)
    old, new = (load(version, inputs, provenance) for version in (OLD, NEW))
    comparisons, resources, summaries, checks, regions = [], [], [], [], []
    for (iso, beam, kind), g in sorted(new.items()):
        row = g['row']
        prior = old.get((iso, beam, kind))
        resources.append(dict(iso=iso, beam=beam, analysis=kind,
            old_wall_s=prior['row']['time_wall_seconds'] if prior else None,
            new_wall_s=row['time_wall_seconds'],
            old_rss_gib=prior['row']['peak_rss_kib']/2**20 if prior else None,
            new_rss_gib=row['peak_rss_kib']/2**20,
            new_peak_rss_kib=row['peak_rss_kib'], command=json.dumps(row['command'])))
        if prior is None:
            continue
        assert normalized_command(row['command']) == normalized_command(prior['row']['command'])
        v = inputs[iso]
        raw, previous = g['raw'], prior['raw']
        result = dict(iso=iso, beam=beam, analysis=kind,
            stdout_byte_equal=raw['stdout.txt'] == previous['stdout.txt'], max_abs_difference=0.)
        if kind == 'bpp':
            a, _ = parse_bpp(raw['bpp.tsv'], v)
            b, _ = parse_bpp(previous['bpp.tsv'], v)
            result.update(old_rows=len(b), new_rows=len(a),
                bpp_byte_equal=raw['bpp.tsv'] == previous['bpp.tsv'],
                omitted_zero_rows=sum(p == 0 and pair not in a for pair, p in b.items()),
                max_abs_difference=max(abs(a.get(p, 0)-b.get(p, 0)) for p in a.keys() | b.keys()),
                bpp_numeric_equal={p:x for p,x in a.items() if x != 0} == {p:x for p,x in b.items() if x != 0})
            assert result['stdout_byte_equal']  # includes MFE structure and energy
        elif kind in ('capr', 'raccess'):
            # Raccess embeds measured elapsed time in an otherwise identical stdout.
            omit_time = lambda text: [line for line in text.splitlines()
                                      if not line.startswith('# elapsed_seconds=')]
            result['stdout_equal_except_elapsed'] = omit_time(raw['stdout.txt']) == omit_time(previous['stdout.txt'])
            assert result['stdout_equal_except_elapsed']
            parse = parse_context if kind == 'capr' else parse_access
            a, b = parse(raw['stdout.txt'], v), parse(previous['stdout.txt'], v)
            result['max_abs_difference'] = max(float(np.max(np.abs(a[k]-b[k]))) for k in a)
        elif kind == 'mea':
            result['structure_equal'] = structures(raw['stdout.txt'], v['n']) == structures(previous['stdout.txt'], v['n'])
            assert result['structure_equal'] and result['stdout_byte_equal']  # includes EA score
        elif kind == 'sampling':
            a, b = structures(raw['stdout.txt'], v['n']), structures(previous['stdout.txt'], v['n'])
            assert len(a) == len(b) == 1000
            result.update(sample_count=1000, structures_equal_in_order=a == b,
                          seed=row['command'][row['command'].index('--seed')+1])
            assert result['structures_equal_in_order'] and result['stdout_byte_equal']
        comparisons.append(result)
    for iso in ('NEAT1_1', 'NEAT1_2'):
        for beam in (100, 200, 400):
            row = new[iso, beam, 'raccess-byloop']['row']
            command = row['command']
            assert '--by-loop' in command and command[command.index('--access-len')+1] == '1,4,8,16'
            assert normalized_command([c for c in command if c != '--by-loop']) == normalized_command(new[iso, beam, 'raccess']['row']['command'])
            total, parts = parse_byloop(new[iso, beam, 'raccess-byloop']['raw']['stdout.txt'], inputs[iso]['n'])
            plain = parse_access(new[iso, beam, 'raccess']['raw']['stdout.txt'], inputs[iso])
            contexts = parse_context(new[iso, beam, 'capr']['raw']['stdout.txt'], inputs[iso])
            residual = max(float(np.max(np.abs(sum(parts[w].values())-total[w]))) for w in WIDTHS)
            diff = max(float(np.max(np.abs(total[w]-plain[w]))) for w in WIDTHS)
            w1diff = max(float(np.max(np.abs(parts[1][c]-contexts[name]))) for c, name in LOOPS.items())
            checks.append(dict(iso=iso, beam=beam, sum_vs_total_max_abs=residual,
                               byloop_total_vs_plain_max_abs=diff, w1_vs_context_max_abs=w1diff))
            assert max(residual, diff, w1diff) < 1e-8
            for w in WIDTHS:
                for c, name in LOOPS.items():
                    values = parts[w][c]
                    i = int(np.argmax(values))
                    summaries.append(dict(iso=iso, beam=beam, window=w, component=name,
                        windows=len(values), mean=float(np.mean(values)), max=float(values[i]),
                        max_start_1based=i+1, count_ge_0_001=int(np.sum(values >= .001))))
            if iso == 'NEAT1_2':
                chosen = slice(21134-1, 21143)
                regions.append(dict(iso=iso, beam=beam, window=8, start_min=21134, start_max=21143,
                    total_mean=float(total[8][chosen].mean()),
                    **{name+'_mean': float(parts[8][c][chosen].mean()) for c, name in LOOPS.items()}))
    equal = all(r['max_abs_difference'] == 0 for r in comparisons)
    assert len({r['energy_parameter_sha256'] for r in provenance['checkpoints']}) == 1
    write_table(args.out / 'resources.tsv', resources)
    write_table(args.out / 'output_comparison.tsv', comparisons)
    write_table(args.out / 'byloop_components.tsv', summaries)
    write_table(args.out / 'byloop_regions.tsv', regions)
    write_json(args.out / 'provenance.json', provenance)
    validation = dict(old=OLD, new=NEW, old_jobs=len(old), new_jobs=len(new),
        new_status_counts=dict(Counter(g['row']['status'] for g in new.values())),
        compared_conditions=len(comparisons), outputs_numerically_identical=equal,
        byloop_checks=checks, verified_artifacts=len(provenance['artifacts']))
    write_json(args.out / 'validation.json', validation)
    (args.out / 'README.md').write_text('''# NEAT1 fixed rerun audit (2026-10-17)

Regenerate with `python experiments/scripts/analyze_e3_neat1_rerun.py`.
The script streams archives in memory; it does not extract tar files or run RNA software.
It checks all planned jobs against complete JSONL records, all artifact sizes/SHA256 and
archive boundaries, GNU time values, unchanged software/binaries, and valid final locks.
Matched commands differ only in executable and temporary paths; inputs have identical hashes.

`resources.tsv` contains old e2cecbb4 and new d02d261c whole-process wall seconds and
peak RSS (KiB / 2^20 = GiB). Each condition is one run. Missing sampling settings are
not imputed. By-loop runs request lengths 1,4,8,16 and six joint probability components.
`output_comparison.tsv` treats absent BPP pairs as zero and compares raw probabilities
before numerical cleanup. MFE and MEA stdout include their structures and scores.
Sampling compares all 1,000 structures in order with seed 20261016.
All Raccess stdout lines except the measured elapsed_seconds metadata agree exactly.

`byloop_components.tsv` gives full-transcript means and maxima by transcript, beam,
window length and loop type; maximum starts are one-based and inclusive.
`byloop_regions.tsv` summarizes the ten highlighted eight-base windows in Figure S6.
`validation.json` records component-sum, total-accessibility and length-one context checks.
`provenance.json` records consumed artifact hashes and checkpoint/software provenance.
The older e3-neat1-20261016 analysis and its variant/empty-scaffold outputs are retained.
''')
    print(json.dumps(validation, indent=2))
    assert equal, 'changed outputs: rerun downstream numerical summaries before using the manuscript values'


if __name__ == '__main__':
    main()
