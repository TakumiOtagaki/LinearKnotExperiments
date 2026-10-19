#!/usr/bin/env python3
"""Recompute saved pseudoknot recovery with coordinate sets and G crossings."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
import statistics as st

import analyze_accuracy as accuracy
from analyze_bprna_long import audit_job
from analyze_e1_pk_stratified import write_tsv
from analyze_pk_membership import CANONICAL, crosses, parse_structure
from analyze_small_rerun_89baea52 import groups
from rna_records import scaffold_pairs
from run_cli import digest, now, write_json
from task_storage import MemoryFile, input_records

ROOT = Path(__file__).resolve().parents[1]
REVISION = '89baea52c9435b18829ef0867a66a24be56d70f9'
SOURCES = [('200', 'turner04', 'pk-recovery-200-turner04'),
           ('200', 'dp09', 'pk-recovery-200-dp09'),
           ('80', 'turner04', 'pk-recovery-80'),
           ('rrna', 'turner04', 'pk-recovery-rrna')]


def crossing_set(pairs, scaffold):
    return {p for p in pairs-scaffold if any(crosses(p, g) for g in scaffold)}


def mfe_pairs(text, n):
    found = []
    for line in text.splitlines():
        cols = line.split()
        if len(cols) >= 2 and len(cols[0]) == n and set(cols[0]) <= set('.()[]{}<>'):
            found.append(parse_structure(cols[0]))
    if len(found) != 1:
        raise ValueError('expected exactly one MFE structure')
    return found[0]


def fraction(a, b):
    return a/b if b else None


def summarize(rows, pairs):
    summaries = []
    for cohort, parameter, _ in SOURCES:
        selected = [r for r in rows if (r['cohort'], r['parameter']) == (cohort, parameter)]
        for definition in ('old_symbols', 'crossing_G'):
            prefix = 'old_' if definition == 'old_symbols' else 'new_'
            row = dict(cohort=cohort, parameter=parameter, definition=definition, rnas=len(selected))
            for field in ('reference_pairs', 'mfe_predicted_pairs', 'mfe_matched_pairs',
                          'reference_pairs_gt_half', 'rnas_mean_probability_gt_half',
                          'reference_reported_zero'):
                row[field] = sum(r[prefix+field] for r in selected)
            row['reference_pairs_gt_half_fraction'] = fraction(row['reference_pairs_gt_half'], row['reference_pairs'])
            row['sensitivity'] = fraction(row['mfe_matched_pairs'], row['reference_pairs'])
            row['ppv'] = fraction(row['mfe_matched_pairs'], row['mfe_predicted_pairs'])
            row['mean_rna_probability'] = st.mean(r[prefix+'mean_probability'] for r in selected)
            row['max_reference_probability'] = max(r[prefix+'max_probability'] for r in selected)
            for label, short in (('median_span_le_100', True), ('median_span_gt_100', False)):
                subset = [r for r in selected if (r[prefix+'median_span'] <= 100) == short]
                row[label+'_rnas'] = len(subset)
                row[label+'_mean_rna_probability'] = st.mean(r[prefix+'mean_probability'] for r in subset) if subset else None
            summaries.append(row)
    distances = []
    for cohort, parameter, _ in SOURCES:
        for definition, flag in (('old_symbols', 'old_reference'), ('crossing_G', 'new_reference')):
            for label, lo, hi in (('le50', 0, 50), ('51to100', 51, 100), ('gt100', 101, 10**9)):
                subset = [p for p in pairs if (p['cohort'], p['parameter']) == (cohort, parameter)
                          and p[flag] and lo <= p['j']-p['i'] <= hi]
                distances.append(dict(cohort=cohort, parameter=parameter, definition=definition,
                    span=label, reference_pairs=len(subset), reference_pairs_gt_half=sum(p['probability'] > .5 for p in subset),
                    mean_probability=st.mean(p['probability'] for p in subset) if subset else None))
    return summaries, distances


def analyze(args):
    old_dir = ROOT/'results/local/small-rerun-89baea52-20261018'
    sources = {'200': ROOT/'results/local/pk-recovery-bprna-20261016/pk_sel.json',
               '80': ROOT/'results/local/pk-recovery-long-rrna-20261016/mid300_700_selection.json'}
    cohorts = {k: {r['id']: r for r in json.loads(p.read_text())} for k, p in sources.items()}
    rrna = json.loads((old_dir/'rrna_reference_metadata.json').read_text())
    old_values = {(r['cohort'], r['par'], r['id']): r for filename in ('pk_results.json', 'mid300_700_results.json')
                  for r in json.loads((old_dir/filename).read_text())}
    rrna_values = {r['input']: r for r in json.loads((old_dir/'rrna_artifact_index.json').read_text())
                   if r['program'] == 'LinearCParty'}
    rows, pairs, changes, audits, task_provenance = [], [], [], [], {}
    for cohort, parameter, group in SOURCES:
        run = ROOT/'results/runs'/(group+'-89baea52-20261018')
        plan = json.loads((run/'array_plan.json').read_text())
        if digest(run/'inputs.tar') != plan['inputs_archive_sha256']:
            raise ValueError('input archive checksum mismatch')
        records, last_reader, seen = {}, None, set()
        for reader, job in groups(group):
            if reader is not last_reader:
                records = input_records(reader)
                last_reader = reader
                task_provenance[str(reader.path.relative_to(ROOT))] = dict(
                    run_manifest_sha256=reader.manifest_sha256,
                    checkpoint_sha256=digest(reader.path/'task_results.jsonl'))
            if job['program'] != 'LinearCParty':
                continue
            audit_job(reader, job)
            if job['input'] in seen:
                raise ValueError('duplicate RNA')
            seen.add(job['input'])
            record = records[job['input']]
            name, seq, scaffold_text = record
            g = scaffold_pairs(scaffold_text)
            source_record = cohorts[cohort][name] if cohort != 'rrna' else None
            if source_record:
                if seq != source_record['seq'] or g != set(map(tuple, source_record['nested'])):
                    raise ValueError('reference S/G differs from frozen run')
                annotated = parse_structure(source_record['ref'])
                canonical = {(i,j) for i,j,_ in annotated if seq[i-1]+seq[j-1] in CANONICAL}
                old_ref = set(map(tuple, source_record['pk']))
                if old_ref != {(i,j) for i,j,o in annotated if o == '[' and (i,j) in canonical}:
                    raise ValueError('old reference rule differs')
                tags = {(i,j): o for i,j,o in annotated}
                new_ref = crossing_set(canonical, g)
                excluded_noncanonical = len(annotated)-len(canonical)
                reference_scope = 'all canonical pairs in saved full annotation'
            else:
                # This archive stores canonical nested G and all saved canonical
                # square/curly targets as coordinates, without a full raw DBN.
                old_ref = set(map(tuple, rrna[name]['pk']))
                canonical = g | old_ref
                new_ref = crossing_set(canonical, g)
                tags = {p: 'saved_canonical_target' for p in old_ref}
                excluded_noncanonical = rrna[name]['noncanonical_dropped']
                reference_scope = 'frozen G plus saved canonical square/curly coordinate targets; raw DBN unavailable'
                if len(seq) != rrna[name]['n'] or len(g) != rrna[name]['nested']:
                    raise ValueError('rRNA reference metadata differs')
            raw_bpp = reader.artifact(job, 'bpp.tsv')
            bpp = accuracy.read_bpp(MemoryFile(raw_bpp), record)['pairs']
            triples = mfe_pairs(reader.artifact(job, 'stdout.txt').decode(), len(seq))
            mfe = {(i,j) for i,j,_ in triples}
            square = {(i,j) for i,j,o in triples if o == '['}
            if not g <= mfe:
                raise ValueError('MFE lost a scaffold pair')
            predicted = crossing_set(mfe, g)
            row = dict(cohort=cohort, parameter=parameter, input=name, length=len(seq),
                       mode=job['mode'], scaffold_pairs=len(g), mfe_total_pairs=len(mfe),
                       mfe_all_added_pairs=len(mfe-g), reference_scope=reference_scope)
            for prefix, ref, pred in (('old_', old_ref, square), ('new_', new_ref, predicted)):
                probs = [bpp.get(p, 0.) for p in ref]
                row.update({prefix+'reference_pairs': len(ref),
                            prefix+'mfe_predicted_pairs': len(pred),
                            prefix+'mfe_matched_pairs': len(ref & pred),
                            prefix+'reference_pairs_gt_half': sum(p > .5 for p in probs),
                            prefix+'rnas_mean_probability_gt_half': int(st.mean(probs) > .5) if probs else 0,
                            prefix+'mean_probability': st.mean(probs) if probs else None,
                            prefix+'max_probability': max(probs) if probs else None,
                            prefix+'reference_reported_zero': sum(p == 0 for p in probs),
                            prefix+'median_span': st.median(j-i for i,j in ref) if ref else None,
                            prefix+'sensitivity': fraction(len(ref & pred),len(ref)),
                            prefix+'ppv': fraction(len(ref & pred),len(pred))})
            if cohort != 'rrna':
                old = old_values[cohort, parameter, name]
                if (row['old_reference_pairs'], row['old_mfe_predicted_pairs'], row['old_mfe_matched_pairs']) != (
                        old['n_pk'], old['mfe_crossing'], old['mfe_crossing_correct']):
                    raise ValueError('old MFE counts not reproduced')
                if [bpp.get(tuple(p), 0.) for p in source_record['pk']] != old['pk_probs']:
                    raise ValueError('old reference probabilities changed')
            else:
                old = rrna_values[name]
                if [bpp.get(tuple(p),0.) for p in old['reference_pk_pairs']] != old['reference_probabilities']:
                    raise ValueError('old rRNA probabilities changed')
            if new_ref & mfe != new_ref & predicted:
                raise ValueError('coordinate sensitivity differs between full MFE and crossing subset')
            row['reference_pairs_added'] = len(new_ref-old_ref)
            row['reference_pairs_removed'] = len(old_ref-new_ref)
            row['mfe_crossing_not_square'] = len(predicted-square)
            row['mfe_square_not_crossing'] = len(square-predicted)
            rows.append(row)
            audits.append(dict(cohort=cohort, parameter=parameter, input=name,
                               canonical_reference_pairs=len(canonical), scaffold_pairs=len(g),
                               noncanonical_reference_pairs_removed=excluded_noncanonical,
                               canonical_nonscaffold_not_crossing=len(canonical-g-new_ref),
                               reference_pairs_added=len(new_ref-old_ref), reference_pairs_removed=len(old_ref-new_ref),
                               added_reference_symbols=json.dumps(dict(Counter(tags[p] for p in new_ref-old_ref))),
                               reference_scope=reference_scope))
            for i,j in sorted(old_ref|new_ref):
                pairs.append(dict(cohort=cohort, parameter=parameter, input=name, i=i, j=j,
                                  reference_symbol=tags[i,j], old_reference=int((i,j) in old_ref),
                                  new_reference=int((i,j) in new_ref), probability=bpp.get((i,j),0.),
                                  in_mfe=int((i,j) in mfe), in_mfe_square=int((i,j) in square)))
            for kind, items in (('reference_added', new_ref-old_ref), ('reference_removed', old_ref-new_ref),
                                ('mfe_crossing_not_square', predicted-square), ('mfe_square_not_crossing', square-predicted)):
                for i,j in sorted(items):
                    changes.append(dict(cohort=cohort, parameter=parameter, input=name,
                                        change=kind, i=i, j=j, reference_symbol=tags.get((i,j), ''),
                                        probability=bpp.get((i,j),0.), in_mfe=int((i,j) in mfe)))
        if seen != (cohorts[cohort].keys() if cohort != 'rrna' else rrna.keys()):
            raise ValueError('incomplete cohort')
        print(group, 'completed', len(seen), flush=True)
    summaries, distances = summarize(rows, pairs)
    args.out.mkdir(parents=True, exist_ok=False)
    for filename, output in (('rna_metrics.tsv', rows), ('pair_metrics.tsv', pairs),
                             ('definition_changes.tsv', changes), ('reference_audit.tsv', audits),
                             ('summary.tsv', summaries), ('span_summary.tsv', distances)):
        write_tsv(args.out/filename, list(output[0]), output)
    write_json(args.out/'summary.json', summaries)
    write_json(args.out/'provenance.json', dict(generated_utc=now(), source_revision=REVISION,
        scripts_sha256={p: digest(Path(__file__).with_name(p)) for p in
                        (Path(__file__).name, 'analyze_small_rerun_89baea52.py', 'analyze_pk_membership.py',
                         'analyze_accuracy.py', 'task_storage.py')},
        reference_sources={str(p.relative_to(ROOT)): digest(p) for p in
                           [*sources.values(), old_dir/'rrna_reference_metadata.json',
                            old_dir/'recovery_summary.json', old_dir/'pk_results.json',
                            old_dir/'mid300_700_results.json', old_dir/'rrna_artifact_index.json']},
        plans={str((ROOT/'results/runs'/(group+'-89baea52-20261018')/'array_plan.json').relative_to(ROOT)):
               digest(ROOT/'results/runs'/(group+'-89baea52-20261018')/'array_plan.json') for _,_,group in SOURCES},
        tasks=task_provenance, input_sets=4, observations=len(rows),
        old_value_checks='All 480 old short/mid MFE counts and all 483 saved probability vectors match exactly',
        definition='crosses((i,j),(k,l)) iff i<k<j<l or k<i<l<j; reference canonical and new relative to G; predicted new relative to G',
        sensitivity='sum |reference_crossing intersect all_MFE_pairs| / sum |reference_crossing|',
        ppv='sum |reference_crossing intersect MFE_crossing| / sum |MFE_crossing|',
        probabilities='strict >0.5; BPP cutoff 0; absent exports counted as zero; no probability clipping',
        rrna_reference_limit='Frozen canonical square/curly coordinate targets plus G; full original DBN was not retained in the available archive',
        storage='checksum-verified in-memory tar member reads; no extraction; no CLI execution',
        outputs_sha256={p.name: digest(p) for p in args.out.iterdir()}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    analyze(parser.parse_args())
