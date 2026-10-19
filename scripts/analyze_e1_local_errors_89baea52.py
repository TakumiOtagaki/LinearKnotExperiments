#!/usr/bin/env python3
"""Count both tails of local probability errors from the complete packed E1 run.

No RNA programs are executed and no tar members are extracted. Workers process
whole saved tasks; all probabilities pass the existing checksum/domain checks.
"""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import csv
import gzip
import json
import math
from pathlib import Path
import time

import numpy as np

import analyze_e1_production as production
from analyze_bprna_short import CHECKS
from analyze_e1_pk_stratified import crossing_pairs, write_tsv
from run_cli import digest, now, write_json

BEAMS = (20, 50, 100, 200)
METRICS = ('bpp_all', 'bpp_crossing_G', 'pseudoloop', 'access_w4_total')
FIELDS = ('input', 'input_sha256', 'tier', 'length', 'beam', 'metric',
          'exact_ge_half', 'beam_ge_half', 'missed_events', 'created_events',
          'max_absolute_error')
PLAN = RUN = None


def tail_counts(exact, beam):
    """Sparse BPP omissions are zero; dense domains are checked by the reader."""
    missed = [(k, p, beam.get(k, 0.)) for k, p in exact.items()
              if p >= .5 and beam.get(k, 0.) < .1]
    created = [(k, exact.get(k, 0.), p) for k, p in beam.items()
               if p >= .5 and exact.get(k, 0.) < .1]
    maximum = max((abs(exact.get(k, 0.) - beam.get(k, 0.))
                   for k in exact.keys() | beam.keys()), default=0.)
    return dict(exact_ge_half=sum(p >= .5 for p in exact.values()),
                beam_ge_half=sum(p >= .5 for p in beam.values()),
                missed_events=len(missed), created_events=len(created),
                max_absolute_error=maximum), missed, created


def domains(record, parsed):
    access = parsed['LinearPKRaccess']['windows']
    w4 = {start: p for (start, width), p in access.items() if width == 4}
    # Export includes zero windows containing scaffold-paired bases.
    if set(w4) != set(range(1, len(record[1]) - 2)):
        raise ValueError('incomplete 4-nt domain')
    return (parsed['LinearCParty']['pairs'],
            crossing_pairs(record, parsed['LinearCParty']),
            dict(enumerate(parsed['LinearPKCapR']['PseudoLoop'], 1)), w4)


def initialize(run):
    global RUN, PLAN
    RUN = Path(run)
    PLAN = production.load_plan(RUN)


def process_task(task):
    modes = ['exact'] + [f'beam{b}_outer' for b in BEAMS]
    checks = dict(tasks=[], selected_job_status_counts=Counter(),
                  excluded_job_status_counts=Counter(), consistency_checks=Counter(),
                  consistency_maxima=dict.fromkeys(CHECKS, 0.))
    reader, hashes, jobs, records = production.open_task(RUN, PLAN, task, modes, checks)
    rows, witnesses = [], {}
    for name, record in records.items():
        reference, _ = production.read_group(reader, jobs, name, 'exact', record, checks)
        reference_domains = domains(record, reference)
        for beam in BEAMS:
            parsed, _ = production.read_group(reader, jobs, name, f'beam{beam}_outer', record, checks)
            for metric, exact, finite in zip(METRICS, reference_domains, domains(record, parsed)):
                counts, missed, created = tail_counts(exact, finite)
                row = dict(input=name, input_sha256=hashes[name], tier=task['tier'],
                           length=len(record[1]), beam=beam, metric=metric, **counts)
                rows.append(row)
                for direction, events in (('missed', missed), ('created', created)):
                    if not events:
                        continue
                    key, p, q = max(events, key=lambda e: abs(e[1] - e[2]))
                    wkey = (task['tier'], beam, metric, direction)
                    witness = dict(tier=task['tier'], beam=beam, metric=metric,
                                   direction=direction, input=name, length=len(record[1]),
                                   coordinate=json.dumps(key), exact_probability=p,
                                   beam_probability=q, absolute_error=abs(p-q))
                    if wkey not in witnesses or abs(p-q) > witnesses[wkey]['absolute_error']:
                        witnesses[wkey] = witness
    return rows, list(witnesses.values()), checks


def analyze(args):
    started = time.monotonic()
    plan = production.load_plan(args.run)
    production.check_reference_plans(plan, plan)
    production.check_exact_modes(plan, plan)
    tasks, hashes = production.cohort(plan, production.SHORT)
    if plan['software_before']['commit'] != '89baea52c9435b18829ef0867a66a24be56d70f9':
        raise ValueError('wrong software revision')
    if {m['name'] for m in plan['config']['modes']} != {'exact', *(f'beam{b}_outer' for b in BEAMS)}:
        raise ValueError('unexpected modes')
    args.out.mkdir(parents=True, exist_ok=False)
    provenance = dict(generated_utc=now(), run=str(args.run),
                      array_plan_sha256=digest(args.run/'array_plan.json'),
                      inputs_archive_sha256=plan['inputs_archive_sha256'],
                      source_revision=plan['software_before']['commit'],
                      tasks=[], selected_job_status_counts=Counter(),
                      excluded_job_status_counts=Counter(), consistency_checks=Counter(),
                      consistency_maxima=dict.fromkeys(CHECKS, 0.),
                      script_sha256={p: digest(Path(__file__).with_name(p)) for p in
                         (Path(__file__).name, 'analyze_e1_production.py', 'analyze_accuracy.py',
                          'analyze_e1_pk_stratified.py', 'task_storage.py', 'bprna_capb.py')},
                      thresholds='missed: exact>=0.5 and beam<0.1; created: exact<0.1 and beam>=0.5',
                      sparse_bpp='cutoff 0; absent pairs are zero; compare union of both exports',
                      domains='all coordinate pairs; G-crossing extension pairs; all bases; all 4-nt windows',
                      rna_denominator='all RNA in length bin, including empty eligible domains',
                      setting='exact uncapped; finite outer ranking, h=d=max(b,30)',
                      storage='read checksum-verified tar members in memory; no extraction',
                      workers=args.workers)
    grouped, maxima, witnesses, seen = defaultdict(Counter), defaultdict(list), {}, set()
    with gzip.open(args.out/'rna_counts.tsv.gz', 'wt', newline='') as stream:
        writer = csv.DictWriter(stream, FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        with ProcessPoolExecutor(args.workers, initializer=initialize, initargs=(str(args.run),)) as pool:
            for number, (rows, examples, checks) in enumerate(pool.map(process_task, tasks), 1):
                writer.writerows(rows)
                for row in rows:
                    rkey = row['input'], row['beam'], row['metric']
                    if rkey in seen:
                        raise ValueError('duplicate comparison')
                    seen.add(rkey)
                    for tier in (row['tier'], 'all'):
                        key = tier, row['beam'], row['metric']
                        group = grouped[key]
                        group['rnas'] += 1
                        for field in ('exact_ge_half', 'beam_ge_half', 'missed_events', 'created_events'):
                            group[field] += row[field]
                        for direction in ('missed', 'created'):
                            group[direction+'_rnas'] += row[direction+'_events'] > 0
                        group['either_rnas'] += bool(row['missed_events'] or row['created_events'])
                        group['max_error_gt_0.1_rnas'] += row['max_absolute_error'] > .1
                        group['max_error_gt_0.25_rnas'] += row['max_absolute_error'] > .25
                        maxima[key].append(row['max_absolute_error'])
                for example in examples:
                    for tier in (example['tier'], 'all'):
                        key = tier, example['beam'], example['metric'], example['direction']
                        if key not in witnesses or example['absolute_error'] > witnesses[key]['absolute_error']:
                            witnesses[key] = dict(example, tier=tier)
                provenance['tasks'].extend(checks['tasks'])
                for key in ('selected_job_status_counts', 'excluded_job_status_counts', 'consistency_checks'):
                    provenance[key].update(checks[key])
                for key in CHECKS:
                    provenance['consistency_maxima'][key] = max(provenance['consistency_maxima'][key], checks['consistency_maxima'][key])
                print(f'task {number}/{len(tasks)}; {time.monotonic()-started:.1f}s', flush=True)
    if len(seen) != len(hashes)*len(BEAMS)*len(METRICS):
        raise ValueError('incomplete cohort')
    summary = []
    for (tier, beam, metric), counts in sorted(grouped.items()):
        row = dict(tier=tier, beam=beam, metric=metric, **counts)
        for key in ('missed_rnas', 'created_rnas', 'either_rnas', 'max_error_gt_0.1_rnas', 'max_error_gt_0.25_rnas'):
            row[key+'_fraction'] = counts[key]/counts['rnas']
        values = maxima[tier, beam, metric]
        row.update(max_error_p95=float(np.quantile(values, .95)),
                   max_error_p99=float(np.quantile(values, .99)), max_error=max(values))
        summary.append(row)
    write_tsv(args.out/'summary.tsv', list(summary[0]), summary)
    examples = [witnesses[k] for k in sorted(witnesses)]
    write_tsv(args.out/'examples.tsv', list(examples[0]), examples)
    provenance.update(rnas=len(hashes), comparisons=len(seen), elapsed_seconds=time.monotonic()-started,
                      outputs_sha256={p.name: digest(p) for p in args.out.iterdir()})
    write_json(args.out/'provenance.json', provenance)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=4)
    analyze(parser.parse_args())
