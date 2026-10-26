#!/usr/bin/env python3
"""Prepare capped real-RNA resource arrays; never calls qsub."""
import json
import math
from pathlib import Path

from file_budget import check_budget, estimate_files
from plan_array import prepare, wall_seconds
from run_cli import ROOT, digest, write_json

BASE='comparison-resources-v0.2.0-20261025'
LOCK=ROOT/'manifests/software-comparison-resources-v0.2.0-20261025.json'


def tier(name,low,high,methods,timeout,memory,workers,max_inputs):
    # Fit every input's complete 1+3 matrix under real process timeouts.
    pairs=max_inputs*len(methods)
    waves=1 if pairs<=workers else pairs/workers+1-1/workers
    bound=120+timeout*(pairs+3*waves)
    wall=math.ceil(bound/.8/300)*300
    if wall>86400:
        raise ValueError('complete capped matrix exceeds scheduler wall limit')
    return dict(name=name,min_length=low,max_length=high,resource='cpu_16',
                h_rt=f'{wall//3600:02d}:{wall%3600//60:02d}:00',
                target_seconds=bound,max_wall_fraction=.8,max_inputs_per_task=max_inputs,
                workers=workers,memory_budget_mib=35095,task_overhead_mib=1024,
                job_limits=dict(timeout_seconds=timeout,memory_limit_mib=memory),
                cost_model=dict(basis='Explicit process timeout ceiling, no length extrapolation; parallel packing uses blocked_random_timeout_bound',
                                seconds_constant_per_condition=timeout,seconds_per_nt_condition=1e-12,
                                mib_per_nt=0,memory_base_mib=memory,safety_factor=1),
                parallel_packing=dict(kind='blocked_random_timeout_bound',timeout_seconds=timeout,
                                      combinations_per_input=len(methods),overhead_seconds=120,
                                      basis='Preserve first-pending launch order: charge all warmups serially plus three measured list bounds work/w+(1-1/w)*largest_job; 120 s staging/checkpoint allowance'))


def main():
    inputs=json.loads((ROOT/'data/comparison-resources-20261025/manifest.json').read_text())
    # Freeze orchestration as well as the adapter's executable dependencies.
    lock=json.loads(LOCK.read_text())
    extras={r['path']:r['sha256'] for r in lock['extra_files']}
    for name in ('run_cli.py','plan_array.py','run_array_task.py','task_storage.py','file_budget.py',
                 'comparison_resources_array_20261025.py','prepare_comparison_resources_20261025.py',
                 'prepare_comparison_inputs_20261025.py'):
        path=ROOT/'scripts'/name
        extras[str(path.relative_to(ROOT))]=digest(path)
    lock['extra_files']=[dict(path=p,sha256=h) for p,h in sorted(extras.items())]
    write_json(LOCK,lock)
    specs=[]
    short=['cparty','exact','beam100','spark','fold100']
    specs.append(('short',short,[r for r in inputs['records'] if r['length']<=1100],
                  [tier('n250',225,275,short,10,1024,8,4),
                   tier('n500',450,550,short,30,1024,8,4),
                   tier('n1000',900,1100,short,120,4096,8,4)]))
    bounded=['beam100','fold100']
    specs.append(('bounded',bounded,[r for r in inputs['records'] if r['length']>1100],
                  [tier('n2000',1800,2200,bounded,180,4096,8,4),
                   tier('n4000',3600,4400,bounded,300,8192,4,2),
                   tier('rnacentral',4793,29106,bounded,900,8192,4,2)]))
    for method in ('cparty','spark'):
        specs.append((method,[method],[r for r in inputs['records'] if r['length']>1100],
                      [tier('n2000',1800,2200,[method],180,4096,8,8),
                       tier('n4000',3600,4400,[method],600,8192,4,4),
                       tier('rna_4793_6990',4793,6990,[method],900,8192,4,4),
                       tier('rna_11832_29106',11832,29106,[method],3600,32768,1,1)]))
    specs.append(('exact',['exact'],[r for r in inputs['records'] if r['group']=='bprna_2000'],
                  [tier('n2000',1800,2200,['exact'],900,4096,8,8)]))
    ntasks=sum(sum(math.ceil(sum(t['min_length']<=r['length']<=t['max_length'] for r in records)/t['max_inputs_per_task'])
                   for t in tiers) for _,_,records,tiers in specs)
    ntiers=sum(len(tiers) for _,_,_,tiers in specs)
    report=dict(schema_version=1,input_manifest_sha256=digest(ROOT/'data/comparison-resources-20261025/manifest.json'),
                software_lock_sha256=digest(LOCK),input_count=len(inputs['records']),
                point_formula='0.06*(0.7*max(task_wall_seconds,300)+0.1*h_rt_seconds)/3600; priority 0',
                tariff_url='https://www.t4.cii.isct.ac.jp/tsubame-kyodo/tariff2026',tariff_checked_date='2026-10-26',
                estimates='Timeout-budget ceilings, not predicted completion times or fitted scaling; storage/retries excluded',
                file_budget=check_budget(estimate_files(ntasks,ntiers),{}),arrays=[],
                exact_exclusions=[dict(input=r['name'],length=r['length'],reason='not planned above 2200 nt; practical exact range; no extrapolated resource value')
                                  for r in inputs['records'] if r['length']>2200])
    config_dir=ROOT/'configs'/BASE
    config_dir.mkdir(exist_ok=False)
    for label,methods,records,tiers in specs:
        name=BASE+'-'+label
        config=dict(schema_version=1,name=name,purpose='formal',energy_parameter_file='params/rna_Turner04.par',
                    model_settings=dict(temperature_celsius=37,dangles=2,build_revision_verified=True,
                                        pf_wmb_closure_policy='LegacyCompatible',fold_wmb_closure_policy='RequireWeakClosure'),
                    environment=dict(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',LC_ALL='C'),
                    timeout_seconds=tiers[0]['job_limits']['timeout_seconds'],memory_limit_mib=tiers[0]['job_limits']['memory_limit_mib'],
                    warmups=1,repetitions=3,order=dict(kind='blocked_random',seed=20261025),software_lock=str(LOCK.relative_to(ROOT)),
                    inputs=records,modes=[dict(name=m,args=['--method',m]) for m in methods],programs=[dict(name='Comparison',args=[])],
                    analysis=dict(scope='PF-only: CParty/exact/beam100; MFE resources: unmodified Spark versus whole fold100 (MFE+PF; no BPP); no energy/structure equality claim',
                                  preserve_failures=True,resource_extrapolation=False))
        settings=dict(schema_version=1,measurement='timing',seed=20261025,compact_plan=True,max_concurrent_tasks=30,
                      task_runner_script='scripts/comparison_resources_array_20261025.py',
                      file_budget=dict(project=102571,filesystem='/gs/bs',hard_file_limit=24000000,max_remaining_fraction=.1),
                      parallel_timing=dict(allowed=True,reason='Author-requested fixed single-thread pool; 8 through 2200 nt, 4 from 3600 nt; one worker when 32 GiB/process is required. Split long methods to preserve explicit censoring budgets.'),
                      tiers=tiers)
        cp,ap=config_dir/(label+'.json'),config_dir/(label+'-array.json')
        write_json(cp,config);write_json(ap,settings)
        plan=prepare(cp,ap,ROOT/'results/runs'/name,ROOT/'jobs'/name,dry_run=True)
        for t in tiers:
            tasks=[x for x in plan['tasks'] if x['tier']==t['name']]
            wall=wall_seconds(t['h_rt'])
            count=len(tasks)
            jobs=sum(len(x['inputs'])*len(methods)*4 for x in tasks)
            report['arrays'].append(dict(run=name,tier=t['name'],methods=methods,inputs=sum(len(x['inputs']) for x in tasks),
                                         tasks=count,cli_calls=jobs,workers=t['workers'],h_rt=t['h_rt'],
                                         timeout_seconds=t['job_limits']['timeout_seconds'],memory_limit_mib=t['job_limits']['memory_limit_mib'],
                                         timeout_budget_hours=sum(x['estimated_seconds'] for x in tasks)/3600,
                                         timeout_budget_points=sum(.06*(.7*max(x['estimated_seconds'],300)+.1*wall)/3600 for x in tasks),
                                         h_rt_points=count*.06*.8*wall/3600,
                                         script='jobs/'+name+'/'+t['name']+'.sh',plan='results/runs/'+name+'/array_plan.json'))
    report.update(tasks=sum(r['tasks'] for r in report['arrays']),cli_calls=sum(r['cli_calls'] for r in report['arrays']),
                  timeout_budget_points=sum(r['timeout_budget_points'] for r in report['arrays']),
                  h_rt_points=sum(r['h_rt_points'] for r in report['arrays']),
                  requested_allocation_hours=sum(r['tasks']*wall_seconds(r['h_rt'])/3600 for r in report['arrays']))
    # Five separate plans also need their own fixed entries, configs and scripts.
    report['total_shared_entry_budget']=report['file_budget']['estimated_shared_entries']+4*14+100
    write_json(ROOT/'docs/comparison_resources_20261025_plan.json',report)
    print(json.dumps({k:report[k] for k in ('tasks','cli_calls','timeout_budget_points','h_rt_points','total_shared_entry_budget')},indent=2))


if __name__=='__main__':
    main()
