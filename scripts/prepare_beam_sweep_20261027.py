#!/usr/bin/env python3
"""Estimate and freeze 49 RNAcentral beam/cap conditions; never submit jobs."""
import argparse
import csv
import json
import math
from pathlib import Path

from file_budget import check_budget, estimate_files
from plan_array import prepare, wall_seconds
from run_cli import ROOT, digest, write_json

BASE = 'beam-sweep-rnacentral-v0.2.1-20261027'
BEAMS = (25, 50, 75, 100, 200, 300, 400)
LOCK = ROOT/'manifests'/('software-'+BASE+'.json')
E2 = ROOT/'results/runs/e2-rnacentral-4f97fd3-20261021/array_plan.json'
GROUPS = (('rna_4_8kb', 1, 8191), ('rna_8_20kb', 8192, 19999), ('rna_20_30kb', 20000, 30000))


def ratios(rows, field, beam):
    """Log-linear between measured beams; extrapolate 25--75 from 100--200."""
    values = []
    for iso in sorted({r['iso'] for r in rows}):
        data = {int(r['beam']): float(r[field]) for r in rows if r['iso'] == iso}
        low, high = (100, 200) if beam <= 200 else (200, 400)
        exponent = math.log(data[high]/data[low])/math.log(high/low)
        values.append(data[low]/data[100]*(beam/low)**exponent)
    return min(values), max(values)


def config(name, inputs, beam, purpose='formal'):
    cap = max(beam, 30)
    return dict(schema_version=1, name=name, purpose=purpose,
                energy_parameter_file='params/rna_Turner04.par', software_lock=str(LOCK.relative_to(ROOT)),
                model_settings=dict(temperature_celsius=37, dangles=2, build_revision_verified=True,
                                    closure='RequireWeakClosure', ranking='outer', h=cap, d=cap),
                environment=dict(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', LC_ALL='C'),
                timeout_seconds=300, memory_limit_mib=1024, warmups=0, repetitions=1,
                inputs=inputs, modes=[dict(name=f'beam{beam}_capb_outer', args=[
                    '-b', str(beam), '--max-hairpin', str(cap), '--max-case3-depth', str(cap), '--ranking', 'outer'])],
                programs=[dict(name='FoldMetadata', executable='FoldMetadata', subcommand='fold',
                               args=['--bpp-cutoff', '0'], outputs=[dict(name='bpp.tsv', flag='--bpp')])],
                analysis=dict(scope='beam width and local caps varied together; finite-width sensitivity, b400 is not exact',
                              timing='one single-thread fold process; MFE+PF+outside+BPP emission; gzip/build excluded',
                              scalars='unchanged v0.2.1 fold computation; full-precision stderr metadata observer',
                              bpp='cutoff zero, all sparse union pairs; absent entries zero',
                              tv='partner and unpaired distribution, scaffold-free positions only',
                              preserve_failures=True, extrapolate_missing=False))


def freeze_lock():
    lock = json.loads(LOCK.read_text())
    extras = {r['path']: r['sha256'] for r in lock['extra_files']}
    for name in ('run_cli.py', 'plan_array.py', 'run_array_task.py', 'task_storage.py', 'file_budget.py',
                 'build_comparison_resources_20261025.py', 'comparison_resources_array_20261025.py',
                 'build_beam_sweep_20261027.py', 'beam_sweep_metadata_20261027.py',
                 'beam_sweep_array_20261027.py', 'prepare_beam_sweep_20261027.py'):
        path = ROOT/'scripts'/name
        extras[str(path.relative_to(ROOT))] = digest(path)
    lock['extra_files'] = [dict(path=k, sha256=v) for k, v in sorted(extras.items())]
    write_json(LOCK, lock)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, default=ROOT/'docs/beam_sweep_rnacentral_20261027_plan.json')
    a = p.parse_args()
    cpdir = ROOT/'configs'/BASE
    for b in BEAMS:
        for path in (ROOT/'results/runs'/(BASE+f'-b{b}'), ROOT/'jobs'/(BASE+f'-b{b}'),
                     cpdir/f'b{b}.json', cpdir/f'b{b}-array.json'):
            if path.exists():
                raise ValueError('refuse to overwrite a prepared sweep: '+str(path))
    old = json.loads(E2.read_text())
    metrics_path = ROOT/'results/local/e2-4f97fd3-20261021/metrics.tsv'
    neat_path = ROOT/'results/local/timing-4f97fd3-20261021/neat1_resources.tsv'
    with metrics_path.open() as stream:
        metrics = {r['name']: r for r in csv.DictReader(stream, delimiter='\t')
                   if r['run'] == 'rnacentral_current' and r['program'] == 'LinearCParty'}
    with neat_path.open() as stream:
        neat = [r for r in csv.DictReader(stream, delimiter='\t') if r['analysis'] == 'bpp']
    # Use the actual adopted plan's input archive, not paths in an older template.
    inputs = []
    for r in old['inputs']:
        inputs.append({k: r[k] for k in ('name', 'path', 'length', 'input_archive', 'input_member', 'sha256')})
    estimates = []
    for item in inputs:
        m = metrics[item['name']]
        t, rss = float(m['wall_seconds_median']), float(m['rss_gib_median'])
        for b in BEAMS:
            tl, th = ratios(neat, 'wall_s', b)
            ml, mh = ratios(neat, 'rss_gib', b)
            estimates.append(dict(name=item['name'], n=item['length'], beam=b, h=max(b,30), d=max(b,30),
                                  time_low_s=5+(t-5)*tl, time_s=5+(t-5)*th,
                                  rss_low_gib=.1+(rss-.1)*ml, rss_gib=.1+(rss-.1)*mh,
                                  complexity_proxy_s=t*(b/100)**2*math.log(b)/math.log(100),
                                  free_bases=int(m['free_bases']), input_sha256=item['sha256']))
    # Include all seven plans and their config/build/analysis files in a combined gate.
    budget = check_budget(dict(estimate_files(49,21), estimated_shared_entries=1050), {})
    freeze_lock()
    cpdir.mkdir(exist_ok=True)
    report = dict(schema_version=1, base=BASE, beams=BEAMS, source_plan=str(E2.relative_to(ROOT)),
                  source_plan_sha256=digest(E2), metrics_sha256=digest(metrics_path), neat1_sha256=digest(neat_path),
                  estimates=estimates, arrays=[], file_budget=budget,
                  tariff_url='https://www.t4.cii.isct.ac.jp/tsubame-kyodo/tariff2026', tariff_checked_date='2026-10-27',
                  point_formula='.06*(.7*max(task_wall_s,300)+.1*h_rt_s)/3600; priority=0',
                  estimate_basis='Per-RNA E2 median times/RSS at b100; min/max NEAT1 log-interpolated ratios; 5s/0.1GiB fixed floors. b<100 extrapolation is uncalibrated. nb^2 log b is a scheduling stress scenario, not a measured bound.',
                  scope='49 fresh conditions, no warmups; one process per task. Separate beam arrays explicitly split the matched RNA matrix across allocations for width-specific time/memory limits.')
    for b in BEAMS:
        name = BASE+f'-b{b}'
        cfg = config(name, inputs, b)
        tiers = []
        for label, low, high in GROUPS:
            rows = [e for e in estimates if e['beam']==b and low <= e['n'] <= high]
            timeout = max(180, math.ceil(2*max(e['time_s'] for e in rows)/60)*60)
            wall = math.ceil((timeout+180)/300)*300
            memory = max(1024, math.ceil((2*max(e['rss_gib'] for e in rows)+.5))*1024)
            if memory+1024 > 35095:
                raise ValueError('process memory cap cannot fit cpu_16')
            observed = {e['name']: dict(sha256=e['input_sha256'],
                        seconds_by_mode={cfg['modes'][0]['name']:e['time_s']},
                        max_seconds_by_mode={cfg['modes'][0]['name']:e['time_s']}, memory_mib=e['rss_gib']*1024) for e in rows}
            tiers.append(dict(name=label,min_length=low,max_length=high,resource='cpu_16',
                              h_rt=f'{wall//3600:02d}:{wall%3600//60:02d}:00', target_seconds=.7*wall,
                              workers=1,max_inputs_per_task=1,memory_budget_mib=35095,task_overhead_mib=1024,
                              job_limits=dict(timeout_seconds=timeout,memory_limit_mib=memory),
                              cost_model=dict(basis=report['estimate_basis'], seconds_per_nt_condition=1e-12,
                                              mib_per_nt=0, safety_factor=1, observed_inputs=observed)))
            for e in rows:
                e.update(timeout_s=timeout, h_rt_s=wall, memory_limit_mib=memory,
                         expected_points=.06*(.7*max(e['time_s']+90,300)+.1*wall)/3600)
        settings = dict(schema_version=1,measurement='accuracy',seed=20261027,compact_plan=True,
                        max_concurrent_tasks=10, task_runner_script='scripts/beam_sweep_array_20261027.py',
                        file_budget=dict(project=102571,filesystem='/gs/bs',hard_file_limit=24000000,max_remaining_fraction=.1),
                        tiers=tiers)
        cfgpath, apath = cpdir/f'b{b}.json', cpdir/f'b{b}-array.json'
        write_json(cfgpath,cfg); write_json(apath,settings)
        plan = prepare(cfgpath,apath,ROOT/'results/runs'/name,ROOT/'jobs'/name,dry_run=True)
        for t in tiers:
            tasks = [x for x in plan['tasks'] if x['tier']==t['name']]
            report['arrays'].append(dict(run=name,tier=t['name'],tasks=len(tasks),beam=b,
                         inputs=[i['name'] for x in tasks for i in x['inputs']],
                         h_rt=t['h_rt'],timeout_s=t['job_limits']['timeout_seconds'],
                         memory_limit_mib=t['job_limits']['memory_limit_mib'],workers=1,
                         script=f'jobs/{name}/{t["name"]}.sh',plan=f'results/runs/{name}/array_plan.json'))
    report.update(process_hours_low=sum(e['time_low_s'] for e in estimates)/3600,
                  process_hours=sum(e['time_s'] for e in estimates)/3600,
                  complexity_proxy_hours=sum(e['complexity_proxy_s'] for e in estimates)/3600,
                  expected_points=sum(e['expected_points'] for e in estimates),
                  timeout_process_hours=sum(e['timeout_s'] for e in estimates)/3600,
                  timeout_points=sum(.06*(.7*max(e['timeout_s']+90,300)+.1*e['h_rt_s'])/3600 for e in estimates),
                  h_rt_allocation_hours=sum(e['h_rt_s'] for e in estimates)/3600,
                  h_rt_points=sum(.06*.8*e['h_rt_s']/3600 for e in estimates))
    write_json(a.out,report)
    fields=list(estimates[0])
    with (ROOT/'docs/beam_sweep_rnacentral_20261027_estimates.tsv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(estimates)
    print(json.dumps({k:report[k] for k in ('process_hours_low','process_hours','complexity_proxy_hours','expected_points','timeout_points','h_rt_points','h_rt_allocation_hours')},indent=2))


if __name__ == '__main__':
    main()
