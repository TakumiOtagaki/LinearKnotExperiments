#!/usr/bin/env python3
"""Read-only pre-submission checks, including live quota; never calls qsub."""
import argparse
import json
from pathlib import Path
import subprocess

from file_budget import check_budget
from plan_array import load_plan, task_jobs, wall_seconds
from prepare_beam_sweep_20261027 import BASE, BEAMS, E2, LOCK
from run_cli import ROOT, digest, now, write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=Path,default=ROOT/'docs/beam_sweep_rnacentral_20261027_validation.json')
    p.add_argument('--skip-staging',action='store_true')
    a=p.parse_args()
    report=json.loads((ROOT/'docs/beam_sweep_rnacentral_20261027_plan.json').read_text())
    old=json.loads(E2.read_text())
    original={r['name']:r['sha256'] for r in old['inputs']}
    seen=set();scripts=[];plans=[]
    for relative in dict.fromkeys(r['plan'] for r in report['arrays']):
        path=ROOT/relative;plan=load_plan(path)
        assert plan['software_lock_sha256']==digest(LOCK)
        assert plan['config']['warmups']==0 and plan['config']['repetitions']==1
        assert plan['config']['programs'][0]['args']==['--bpp-cutoff','0']
        assert plan['pending_scaffolds'] is False
        assert {r['name']:r['sha256'] for r in plan['inputs']}==original
        b=int(plan['config']['modes'][0]['args'][1]);assert b in BEAMS
        for task in plan['tasks']:
            schedule=task_jobs(plan,task)
            assert len(schedule)==1 and task['workers']==1 and task['requested_resource']=='cpu_16=1'
            assert task['memory_cap_envelope_mib']<=35095
            assert task['estimated_seconds']<=.7*wall_seconds(task['h_rt'])
            assert task['job_limits']['timeout_seconds']+180<=wall_seconds(task['h_rt'])
            assert len(task['inputs'])==1
            item=task['inputs'][0]
            assert (item['name'],b) not in seen;seen.add((item['name'],b))
            metadata_path=Path(plan['output_dir'])/task['directory']/'task.json'
            metadata=json.loads(metadata_path.read_text())
            assert metadata['status']=='not_started' and not metadata['attempts']
            assert not (metadata_path.parent/'run').exists()
        for script in plan['job_scripts']:
            assert digest(Path(script['path']))==script['sha256']
            subprocess.run(['bash','-n',script['path']],check=True)
            scripts.append(script)
        plans.append(dict(path=relative,sha256=digest(path),tasks=len(plan['tasks'])))
    assert seen=={(name,b) for name in original for b in BEAMS}
    staged=[]
    if not a.skip_staging:
        for b,tier,task in ((25,'rna_4_8kb',1),(200,'rna_8_20kb',1),(400,'rna_20_30kb',2)):
            command=['python3','scripts/beam_sweep_array_20261027.py','--plan',
                     f'results/runs/{BASE}-b{b}/array_plan.json','--tier',tier,'--task-id',str(task),'--dry-run']
            result=subprocess.run(command,check=True,text=True,capture_output=True)
            staged.append(json.loads(result.stdout))
    quota=check_budget(dict(report['file_budget'],estimated_shared_entries=1050),{})
    result=dict(created_utc=now(),status='passed',conditions=len(seen),arrays=len(scripts),tasks=len(seen),
                formal_tasks_not_started=True,plans=plans,script_syntax_checks=len(scripts),
                staged_dependency_checks=staged,file_budget=quota,
                software_lock_sha256=digest(LOCK),cli_executed=False,qsub_executed=False,
                notes='Staging checks verify native+metadata binaries, clean archive source, cache, params, and frozen runner dependencies without executing RNA calculations.')
    write_json(a.out,result)
    print(json.dumps({k:result[k] for k in ('status','conditions','arrays','tasks','formal_tasks_not_started')},indent=2))


if __name__=='__main__':
    main()
