"""Incrementally add capb to a validated packed short-array analysis.

Raw artifacts are read through streams. A single node-local SQLite file holds
new scalar rows and old resource observations until the paired TSV is written.
"""
from collections import Counter, OrderedDict
import csv
import gc
import json
from pathlib import Path
import sqlite3
import tempfile
import time

import analyze_accuracy as accuracy
import analyze_bprna_short as short
from run_cli import digest, git_state, now
from task_storage import FORMAT, RunReader, input_records


def check_task(reader, plan, task):
    saved = reader.manifest
    if saved['status'] not in ('success', 'failed'):
        raise ValueError('task is incomplete: ' + str(reader.path))
    if accuracy.model_signature(saved) != accuracy.model_signature(plan):
        raise ValueError('task models or binaries differ from plan')
    if (saved.get('software_changed_during_run') is not False or
            saved.get('binaries_changed_during_run') is not False or
            saved.get('software_lock_valid_after_run') is not True):
        raise ValueError('task dependencies were not verified unchanged')
    hashes = {i['name']: i['sha256'] for i in saved['inputs']}
    if hashes != {i['name']: i['sha256'] for i in task['inputs']}:
        raise ValueError('task input hashes differ from plan')
    jobs = {}
    for job in saved['jobs']:
        key = (job['input'], job['mode'], job['program'])
        if job['warmup'] or job['repetition'] != 0 or key in jobs:
            raise ValueError('duplicate observation, warmup or repetition')
        jobs[key] = job
    return hashes, jobs


def without_caps(command):
    args = accuracy.algorithm_arguments(command)
    for flag in ('--max-hairpin', '--max-case3-depth'):
        if flag in args:
            i = args.index(flag)
            del args[i:i+2]
    return args


def check_capb_command(job, default_job, beam):
    args = accuracy.algorithm_arguments(job['command'])
    for flag, expected in (('-b', beam), ('--max-hairpin', max(beam, 30)),
                           ('--max-case3-depth', max(beam, 30))):
        if args.count(flag) != 1 or args[args.index(flag)+1] != str(expected):
            raise ValueError('capb algorithm arguments differ from cap rule')
    if without_caps(job['command']) != without_caps(default_job['command']):
        raise ValueError('capb and default non-cap arguments differ')


def resource_fields(jobs, name, mode):
    return {f'{p}_{k}': jobs[name, mode, p].get(k) for p in accuracy.PROGRAMS
            for k in ('elapsed_seconds', 'time_wall_seconds', 'peak_rss_kib')}


def compare_row(name, mode, task, parsed, reference, record, rerun_flags):
    check = accuracy.check_consistency(name, mode, 0, parsed)
    pair = accuracy.compare_bpp(parsed['LinearCParty'], reference['LinearCParty'], record)
    context = accuracy.compare_contexts(parsed['LinearPKCapR'], reference['LinearPKCapR'], record)
    access = accuracy.compare_accessibility(parsed['LinearPKRaccess'], reference['LinearPKRaccess'], record)
    row = dict(input=name, tier=task['tier'], task=task['task_id'], mode=mode,
               cap='capb', beam=int(mode[4:].split('_')[0]), analysis_status='success')
    row.update({k: pair[k] for k in short.CORE if k in pair})
    row.update({k: access[k] for k in short.CORE if k in access})
    row['context_tv_mean_Gfree'] = context['context_tv_mean_Gfree']
    row.update({'context_'+c+'_mae_Gfree': v for c, v in context['context_mae_Gfree'].items()})
    for w, metrics in access['accessibility'].items():
        row.update({f'access_w{w}_{k}': v for k, v in metrics.items() if f'access_w{w}_{k}' in short.ACCESS})
    row.update({k: check[k] for k in short.CHECKS})
    row.update(rt_kcal_mol=access['rt_kcal_mol'], exact_free_energy_kcal_mol=reference['LinearPKRaccess']['free_energy'])
    row.update({'exact_'+p+'_rerun': flag for p, flag in rerun_flags.items()})
    return row, check


def extend_capb(previous, capb_run, output):
    start = time.monotonic()
    old = json.loads((previous/'provenance.json').read_text())
    validation = json.loads((previous/'validation.json').read_text())
    if not validation['all_passed'] or digest(previous/'metrics.tsv') != old['metrics_sha256']:
        raise ValueError('previous analysis must be validated and its TSV unchanged')
    if old.get('caps', list(short.CAPS)) != list(short.CAPS):
        raise ValueError('extend requires the original default/nocap analysis')
    main_run = Path(old['main_run']['path'])
    main = json.loads((main_run/'array_plan.json').read_text())
    capb = json.loads((capb_run/'array_plan.json').read_text())
    if any(p.get('output_format') != FORMAT for p in (main, capb)):
        raise ValueError('requires packed task-jsonl-tar-v1 arrays')
    if digest(main_run/'array_plan.json') != old['main_run']['array_plan_sha256']:
        raise ValueError('previous array plan changed')
    if accuracy.model_signature(main) != accuracy.model_signature(capb):
        raise ValueError('capb and main models or binaries differ')
    main_tasks = [t for t in main['tasks'] if t['tier'] in short.TIERS]
    main_inputs = {i['name']: i['sha256'] for t in main_tasks for i in t['inputs']}
    capb_inputs = {i['name']: i['sha256'] for t in capb['tasks'] for i in t['inputs']}
    if capb_inputs != main_inputs or len(capb_inputs) != old['inputs']:
        raise ValueError('capb requires the identical complete S/G cohort')
    assignments = {i['name']: n for n, t in enumerate(main_tasks) for i in t['inputs']}
    reruns = accuracy.ExactReruns([RunReader(p['path']) for p in old['exact_rerun_runs']], main)
    if reruns.provenance != old['exact_rerun_runs']:
        raise ValueError('previous exact rerun manifests changed')
    provenance = {**old, 'schema_version': 2, 'generated_utc': now(),
                  'caps': ['default', 'capb', 'nocap'], 'modes': short.measured_modes(['default','capb','nocap']),
                  'previous_analysis': {'path': str(previous), 'metrics_sha256': old['metrics_sha256'],
                      'provenance_sha256': digest(previous/'provenance.json'),
                      'validation_sha256': digest(previous/'validation.json'),
                      'script_sha256': old['script_sha256']},
                  'capb_run': short.plan_provenance(capb_run, capb), 'capb_tasks': [],
                  'analysis_git_state': git_state(accuracy.ROOT),
                  'script_sha256': {name: digest(Path(__file__).parent/name) for name in
                      ('analyze_bprna_short.py','bprna_capb.py','analyze_accuracy.py','task_storage.py','rna_records.py')},
                  'adopted_rule': 'h=d=max(b,30); b20 aliases beam20_default; b50/100/200 use measured capb',
                  'storage': 'tar member streams; task-bounded raw reads; one temporary SQLite file outside /gs/bs; one output TSV',
                  'resource_caveat': 'Reference values only: concurrent accuracy runs on different allocations/dates, one repetition.'}
    statuses, checks, reference_checks, seen = Counter(), Counter(), Counter(), set()
    output.mkdir(parents=True, exist_ok=False)
    with tempfile.TemporaryDirectory(prefix='e1-capb-', dir='/tmp') as scratch:
        db = sqlite3.connect(str(Path(scratch)/'scalars.sqlite'))
        db.execute('CREATE TABLE scalars (name TEXT, mode TEXT, payload TEXT, PRIMARY KEY(name,mode))')
        db.execute('CREATE TABLE resources (name TEXT, mode TEXT, payload TEXT, PRIMARY KEY(name,mode))')
        original_statuses = Counter()
        old_tasks = {t['path']: t for t in old['tasks']}
        # Audit the reused scalar data against unchanged task checkpoints; retain only resources.
        for task in main_tasks:
            path = main_run/task['directory']/'run'
            reader = RunReader(path)
            _, jobs = check_task(reader, main, task)
            prior = old_tasks[str(path)]
            if reader.manifest_sha256 != prior['run_manifest_sha256'] or digest(path/'task_results.jsonl') != prior['task_results_sha256']:
                raise ValueError('previous task manifest or checkpoint changed')
            original_statuses.update(j['status'] for j in jobs.values())
            for name in (i['name'] for i in task['inputs']):
                for mode in short.measured_modes(short.CAPS):
                    db.execute('INSERT INTO resources VALUES (?,?,?)', (name, mode, json.dumps(resource_fields(jobs,name,mode))))
            del reader, jobs
        db.commit()
        if dict(original_statuses) != old['original_job_status_counts']:
            raise ValueError('previous status counts changed')
        cache = OrderedDict()
        def main_reader(name):
            n = assignments[name]
            if n not in cache:
                if len(cache) == 2:
                    cache.popitem(last=False)
                task = main_tasks[n]
                reader = RunReader(main_run/task['directory']/'run')
                _, jobs = check_task(reader, main, task)
                cache[n] = (reader, jobs, input_records(reader))
            cache.move_to_end(n)
            return cache[n]
        for ti, task in enumerate(capb['tasks'], 1):
            if task['tier'] not in short.TIERS:
                raise ValueError('capb task outside short tiers')
            path = capb_run/task['directory']/'run'
            reader = RunReader(path)
            hashes, jobs = check_task(reader, capb, task)
            records = input_records(reader)
            expected = {(name,f'beam{b}_capb',p) for name in records for b in short.BEAMS[1:] for p in accuracy.PROGRAMS}
            if set(jobs) != expected:
                raise ValueError('incomplete capb matrix')
            statuses.update(j['status'] for j in jobs.values())
            provenance['capb_tasks'].append({'path': str(path), 'run_manifest_sha256': reader.manifest_sha256,
                'task_results_sha256': digest(path/'task_results.jsonl'), 'status': reader.manifest['status'],
                'inputs': len(records), 'jobs': len(jobs)})
            for name, record in records.items():
                if name in seen or not 50 <= len(record[1]) < 400:
                    raise ValueError('duplicate RNA or length outside short cohort')
                seen.add(name)
                source, main_jobs, main_records = main_reader(name)
                if main_records[name] != record:
                    raise ValueError('capb requires identical S and G')
                reference, flags = {}, {}
                for p in accuracy.PROGRAMS:
                    selected_source, job = reruns.replace(source, main_jobs[name,'exact',p], record, hashes[name])
                    if job['status'] != 'success':
                        raise ValueError('exact reference unavailable')
                    flags[p] = int(selected_source is not source)
                    reference[p] = accuracy.read_job(selected_source, job, record)
                reference_check = accuracy.check_consistency(name,'exact',0,reference)
                reference_checks[reference_check['status']] += 1
                for b in short.BEAMS[1:]:
                    mode = f'beam{b}_capb'
                    parsed = {}
                    errors = []
                    for p in accuracy.PROGRAMS:
                        job = jobs[name,mode,p]
                        check_capb_command(job,main_jobs[name,f'beam{b}_default',p],b)
                        if job['status'] != 'success':
                            errors.append({'program':p,'status':job['status'],'error':job.get('error')})
                            continue
                        try:
                            parsed[p] = accuracy.read_job(reader,job,record)
                        except (ValueError, KeyError, OSError) as exc:
                            errors.append({'program':p,'error':str(exc)})
                    if errors:
                        row = dict(input=name,tier=task['tier'],task=task['task_id'],mode=mode,cap='capb',beam=b,analysis_status='unavailable')
                        provenance.setdefault('analysis_errors',[]).append({'input':name,'mode':mode,'errors':errors})
                    else:
                        row, check = compare_row(name,mode,task,parsed,reference,record,flags)
                        checks[check['status']] += 1
                        if check['status'] != 'success':
                            provenance.setdefault('consistency_failures',[]).append(check)
                    row.update(resource_fields(jobs,name,mode))
                    db.execute('INSERT INTO scalars VALUES (?,?,?)',(name,mode,json.dumps(row)))
            db.commit()
            print(f'capb task {ti}/{len(capb["tasks"])} {task["tier"]}/{task["task_id"]}: {len(seen)} RNA, {time.monotonic()-start:.1f}s',flush=True)
            gc.collect()
        if seen != set(main_inputs) or reruns.used != set(reruns.jobs):
            raise ValueError('missing RNA or unused exact reruns')
        if sorted(reruns.replacements,key=lambda r:r['id']) != sorted(old['exact_replacements'],key=lambda r:r['id']):
            raise ValueError('exact replacements differ from previous analysis')
        counts = Counter()
        with (previous/'metrics.tsv').open() as prior, (output/'metrics.tsv').open('w',newline='') as stream:
            rows = csv.DictReader(prior,delimiter='\t')
            writer = csv.DictWriter(stream,short.FIELDS,delimiter='\t',lineterminator='\n',restval='NA')
            writer.writeheader()
            for _ in range(old['inputs']):
                block = [next(rows) for _ in range(8)]
                name = block[0]['input']
                if [r['mode'] for r in block] != short.measured_modes(short.CAPS) or any(r['input'] != name for r in block):
                    raise ValueError('previous paired block is invalid')
                capb_rows = [json.loads(db.execute('SELECT payload FROM scalars WHERE name=? AND mode=?',(name,f'beam{b}_capb')).fetchone()[0]) for b in short.BEAMS[1:]]
                for row in block[:4] + capb_rows + block[4:]:
                    if row['cap'] != 'capb':
                        row.update(json.loads(db.execute('SELECT payload FROM resources WHERE name=? AND mode=?',(name,row['mode'])).fetchone()[0]))
                    else:
                        row.update({k:block[0][k] for k in short.DESCRIPTORS})
                    counts[row['analysis_status']] += 1
                    writer.writerow({k:'NA' if v is None else v for k,v in row.items()})
            if next(rows,None) is not None:
                raise ValueError('previous TSV has excess rows')
        db.close()
    combined = Counter(old['original_job_status_counts']) + statuses
    combined_checks = Counter(old['consistency_check_counts']) + checks
    provenance.update(inputs=len(seen),capb_job_status_counts=dict(statuses),
        original_job_status_counts=dict(combined),
        effective_job_status_counts={'success':combined['success']+len(reruns.used),
            'failed':sum(v for k,v in combined.items() if k!='success')-len(reruns.used)},
        analysis_status_counts=dict(counts),consistency_check_counts=dict(combined_checks),
        capb_consistency_check_counts=dict(checks),rechecked_exact_consistency_counts=dict(reference_checks),
        elapsed_seconds=time.monotonic()-start,metrics_sha256=digest(output/'metrics.tsv'))
    (output/'provenance.json').write_text(json.dumps(provenance,indent=2,allow_nan=False)+'\n')
    if set(counts) != {'success'} or set(checks) != {'success'} or set(reference_checks) != {'success'}:
        raise ValueError('analysis preserved failed or inconsistent observations; review provenance')
    return provenance
