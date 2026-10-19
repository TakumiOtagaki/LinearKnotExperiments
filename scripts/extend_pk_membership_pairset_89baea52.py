#!/usr/bin/env python3
"""Check only new crossing-G targets, then join the saved 983-pair audit."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path

from analyze_pk_membership import execute, render

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    old = ROOT / 'results/local/pk-membership-89baea52-20261018'
    pairset = ROOT / 'results/local/pk-recovery-pairset-89baea52-20261019'
    cohort_path = ROOT / 'results/local/pk-recovery-bprna-20261016/pk_sel.json'
    cohort = {r['id']: r for r in json.loads(cohort_path.read_text())}
    plan_path = ROOT / 'results/runs/pk-membership-89baea52-20261018/array_plan.json'
    lock = json.loads(plan_path.read_text())['software_lock']
    driver = ROOT / 'build/membership-89baea52/membership_driver'
    parameter = ROOT / 'deps/LinearCParty-89baea52/params/rna_Turner04.par'
    expected_driver = next(x['sha256'] for x in lock['extra_files']
                           if x['path'] == str(driver.relative_to(ROOT)))
    assert sha(driver) == expected_driver
    assert sha(parameter) == lock['energy_parameter_sha256']
    inputs = [plan_path, cohort_path, old/'pairs.tsv', old/'summary.json',
              pairset/'pair_metrics.tsv', driver, parameter]
    before = {str(p.relative_to(ROOT)): sha(p) for p in inputs}
    original = {(r['id'], int(r['i']), int(r['j'])): r for r in
                csv.DictReader((old/'pairs.tsv').open(), delimiter='\t')}
    selected = [r for r in csv.DictReader((pairset/'pair_metrics.tsv').open(), delimiter='\t')
                if r['cohort'] == '200' and r['parameter'] == 'turner04'
                and r['new_reference'] == '1']
    additions = defaultdict(list)
    for r in selected:
        key = r['input'], int(r['i']), int(r['j'])
        if key not in original:
            additions[r['input']].append(key[1:])
    assert len(selected) == 1040 and sum(map(len, additions.values())) == 57
    fresh = {}
    # One process per RNA; commands and full stdout results share one JSONL.
    with (args.out/'raw.jsonl').open('w') as raw:
        for name, pairs in sorted(additions.items()):
            record = cohort[name]
            targets = {f'added_{i}_{j}': render(len(record['seq']), [(i, j)], '[')
                       for i, j in pairs}
            prepared = dict(record=record, scaffold=render(len(record['seq']), record['nested']),
                            unique=targets, aliases={k: k for k in targets})
            result = execute(prepared, driver, parameter, timeout=180, energy_only=True)
            raw.write(json.dumps(result, allow_nan=False)+'\n')
            raw.flush()
            if result.get('returncode') != 0:
                raise RuntimeError(f'{name}: observer failed; raw attempt preserved')
            for i, j in pairs:
                evaluation = result['result']['evaluations'][f'added_{i}_{j}']
                if evaluation['reachable']:
                    assert {tuple(p) for p in evaluation['pairs']} == {(i, j)}
                    assert abs(evaluation['trace_sum_error']) < 1e-8
                fresh[name, i, j] = evaluation
            print(name, len(pairs), 'checked', flush=True)
    rows = []
    for r in selected:
        name, i, j = r['input'], int(r['i']), int(r['j'])
        prior = original.get((name, i, j))
        if prior:
            assert float(prior['probability']) == float(r['probability'])
            reachable = prior['single_status'] == 'reachable'
            source = 'saved_89baea52_single_target'
        else:
            reachable = fresh[name, i, j]['reachable']
            source = 'added_89baea52_single_target'
        rows.append(dict(input=name, i=i, j=j, span=j-i,
                         probability=float(r['probability']), reachable=int(reachable),
                         evidence=source))
    with (args.out/'pairs.tsv').open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader(); w.writerows(rows)
    spans = []
    for name, lo, hi in [('le50', 0, 50), ('51to100', 51, 100), ('gt100', 101, 10**9)]:
        subset = [r for r in rows if lo <= r['span'] <= hi]
        spans.append(dict(span=name, pairs=len(subset), reachable=sum(r['reachable'] for r in subset),
                          probability_gt_half=sum(r['probability'] > .5 for r in subset),
                          mean_probability=sum(r['probability'] for r in subset)/len(subset)))
    summary = dict(revision=lock['source_commit'], rnas=len({r['input'] for r in rows}),
                   reference_pairs=len(rows), reused_pairs=len(original), added_pairs=len(fresh),
                   added_rnas=len(additions), reachable=sum(r['reachable'] for r in rows),
                   positive_marginals=sum(r['probability'] > 0 for r in rows), distance_groups=spans,
                   full_reference=json.loads((old/'summary.json').read_text())['membership']['full'],
                   full_reference_energy_gap=json.loads((old/'summary.json').read_text())['energy_gaps']['full'])
    (args.out/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    after = {str(p.relative_to(ROOT)): sha(p) for p in inputs}
    assert before == after
    provenance = dict(inputs_sha256=before, inputs_unchanged=True, script_sha256=sha(Path(__file__)),
                      outputs_sha256={p.name: sha(p) for p in args.out.iterdir() if p.is_file()})
    (args.out/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
