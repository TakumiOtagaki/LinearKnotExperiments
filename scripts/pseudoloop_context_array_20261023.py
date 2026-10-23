#!/usr/bin/env python3
"""Plan/check/run profile-only E1 arrays; never submits scheduler jobs.

Use the old certified E1 plans and measured profile jobs, with the same inputs,
modes, caps, rank, seed, parameter and per-process limits. The archived release
is built on node-local storage once; tasks unpack its hashed runtime locally.
TaskStore supplies the repository's tar-before-JSONL durability and recovery.
"""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import fcntl
import hashlib
import io
import json
import math
import os
from pathlib import Path
import random
import shlex
import shutil
import signal
import subprocess
import tarfile
import time

from file_budget import check_budget, estimate_files
from task_storage import TaskStore, local_workspace

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT/'results/runs/pseudoloop-context-0.2.0-20261023'
CONTEXTS = ['Stem','Hairpin','Bulge','Internal','Multiloop','Exterior','PseudoLoop']
OLD = {'short':'e1-bprna-short-89baea52-20261018',
       'long':'e1-bprna-long-89baea52-cpu16-20261018'}


def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat()


def dump(p,d):
    Path(p).write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False)+'\n')


def upper_bound(observations, workers):
    # List-scheduling bound, inflated 2x; 0.5 s per command and 80% efficiency.
    durations=[2*(r['elapsed_seconds']+.5) for r in observations]
    return 60+(sum(durations)/workers+(1-1/workers)*max(durations,default=0))/.8


def plan(args):
    out=args.out.resolve()
    if out.exists(): raise ValueError('choose a new plan output directory')
    source,build=args.source.resolve(),args.build.resolve()
    if any(str(p).startswith('/gs/bs/') for p in (source,build)):
        raise ValueError('source and build must be node-local archives')
    assert (source/'VERSION').read_text().strip()=='0.2.0'
    report={'schema_version':1,'kind':'pseudoloop-profile-only-v1','created_utc':now(),
        'source_commit':subprocess.check_output(['git','-C',str(ROOT/'../../LinearKnot'),
            'rev-parse','f93b55e'],text=True).strip(),'source_acquisition':'git archive; node-local build; no worktree',
        'binary_sha256':sha(build/'linearknot'),'parameter_sha256':sha(source/'params/rna_Turner04.par'),
        'cmake_cache_sha256':sha(build/'CMakeCache.txt'),'script_sha256':sha(__file__),
        'task_storage_sha256':sha(ROOT/'scripts/task_storage.py'),
        'measurement':'context accuracy only; concurrent wall/RSS are not new performance results',
        'seed':20261011,'warmups':0,'repetitions':1,'programs':[{'subcommand':'profile','args':[]}],
        'old_plans':{},'input_archives':{},'modes':{},'tiers':{},'inputs':{},'tasks':[],
        'calibration':'matched old profile elapsed_seconds; factor 2 plus 0.5 s/job, list bound / 0.8 + 60 s/task',
        'point_formula':'0.06*(0.7*max(elapsed_seconds,300)+0.1*h_rt_seconds)/3600; priority 0',
        'tariff_url':'https://www.t4.cii.isct.ac.jp/tsubame-kyodo/tariff2026'}
    allobs={}; references=[]
    for cohort,run in OLD.items():
        oldroot=ROOT/'results/runs'/run; pp=oldroot/'array_plan.json'
        p=json.loads(pp.read_text()); references.append(pp)
        assert p['config']['warmups']==0 and p['config']['repetitions']==1
        assert p['software_lock']['energy_parameter_sha256']==report['parameter_sha256']
        report['old_plans'][cohort]={'path':str(pp.relative_to(ROOT)),'sha256':sha(pp),
            'cohort_size':len(p['inputs']),'source_commit':p['software_before']['commit'],
            'config':p['config'],'input_hash_set_sha256':hashlib.sha256(json.dumps(
                sorted((r['name'],r['sha256']) for r in p['inputs'])).encode()).hexdigest()}
        inputsfile=oldroot/'inputs.tar'
        assert sha(inputsfile)==p['inputs_archive_sha256']
        report['input_archives'][cohort]={'path':str(inputsfile.relative_to(ROOT)),'sha256':sha(inputsfile)}
        report['modes'][cohort]=p['config']['modes']
        tierinputs=defaultdict(list)
        for task in p['tasks']:
            tier=task['tier']
            if tier not in report['tiers']:
                cap=task['job_limits']['memory_limit_mib']
                workers=min(8,(35045-1024)//cap)
                wall=900 if tier in ('nt50_99','nt100_199') else 1800
                report['tiers'][tier]={'cohort':cohort,'resource':'cpu_16','workers':workers,
                    'h_rt_seconds':wall,'memory_budget_mib':35045,'task_overhead_mib':1024,
                    'memory_limit_mib':cap,'timeout_seconds':task['job_limits']['timeout_seconds'],
                    'target_wall_fraction':.65,'max_inputs_per_task':500}
            for r in task['inputs']: tierinputs[tier].append(r['name'])
            checkpoint=oldroot/task['directory']/'run/task_results.jsonl'
            checkpoint_hash=sha(checkpoint)
            with checkpoint.open() as f:
                for line in f:
                    r=json.loads(line)
                    if r['program']!='LinearPKCapR': continue
                    assert not r['warmup'] and r['repetition']==0
                    key=(r['input'],r['mode']); assert key not in allobs
                    if r['status']!='success': raise ValueError('prior profile failed: '+str(key))
                    allobs[key]={'input':r['input'],'mode':r['mode'],'elapsed_seconds':r['elapsed_seconds'],
                        'peak_rss_kib':r['peak_rss_kib'],'stdout_bytes':r['stdout_summary']['bytes'],
                        'stdout_sha256':r['stdout_summary']['sha256'],
                        'archive':str((checkpoint.parent/'task_outputs.tar').relative_to(ROOT)),
                        'member':r['directory']+'/stdout.txt','checkpoint_sha256':checkpoint_hash}
        for r in p['inputs']:
            name=r['name']; assert name not in report['inputs']
            obs=[allobs[name,m['name']] for m in p['config']['modes']]
            report['inputs'][name]={'name':name,'length':r['length'],'sha256':r['sha256'],
                'sequence_sha256':r['sequence_sha256'],'cohort':cohort,'member':r['input_member'],
                'observations':obs}
        for tier,names in tierinputs.items():
            settings=report['tiers'][tier]; workers=settings['workers']; wall=settings['h_rt_seconds']
            # A single sequence's full matrix stays together, even for slow exact.
            largest=max(upper_bound(report['inputs'][name]['observations'],workers) for name in names)
            if largest > .65*wall:
                wall=math.ceil(largest/.65/300)*300; settings['h_rt_seconds']=wall
            current=[]; obs=[]; tid=0
            def save():
                nonlocal current,obs,tid
                if not current:return
                tid+=1; estimate=upper_bound(obs,workers)
                report['tasks'].append({'tier':tier,'task_id':tid,'inputs':current,
                    'profile_jobs':len(obs),'estimated_wall_seconds':estimate,
                    'observed_sequential_seconds':sum(r['elapsed_seconds'] for r in obs)})
                current=[];obs=[]
            for name in names:
                more=report['inputs'][name]['observations']
                if current and (len(current)>=500 or upper_bound(obs+more,workers)>.65*wall):save()
                current.append(name);obs+=more
            save()
    assert len(report['inputs'])==26063 and len(allobs)==24625*5+1438*6
    count=len(report['tasks']); ntiers=len(report['tiers'])
    report['file_budget']=check_budget(estimate_files(count,ntiers),{'file_budget':{'hard_file_limit':24000000}})
    report['profile_jobs']=len(allobs)
    report['predicted_raw_profile_stdout_bytes']=sum(r['stdout_bytes'] for r in allobs.values())
    report['estimated_allocation_hours']=sum(t['estimated_wall_seconds']/3600 for t in report['tasks'])
    report['requested_allocation_hours']=sum(report['tiers'][t['tier']]['h_rt_seconds']/3600 for t in report['tasks'])
    report['estimated_points']=sum(.06*(.7*max(t['estimated_wall_seconds'],300)+.1*report['tiers'][t['tier']]['h_rt_seconds'])/3600 for t in report['tasks'])
    report['points_at_h_rt']=.06*.8*report['requested_allocation_hours']
    report['old_observed_profile_seconds']=sum(r['elapsed_seconds'] for r in allobs.values())
    out.mkdir(parents=True)
    # Runtime/source/build logs are a few packed files, never a shared clone/build.
    with tarfile.open(out/'runtime.tar.gz','w:gz') as tf:
        tf.add(build/'linearknot',arcname='linearknot')
        tf.add(source/'params/rna_Turner04.par',arcname='rna_Turner04.par')
        tf.add(build/'CMakeCache.txt',arcname='CMakeCache.txt')
        for filename in ['configure.log','build.log','build_verified.log']:
            p=source.parent/filename
            if p.exists():tf.add(p,arcname=filename)
    report['runtime_sha256']=sha(out/'runtime.tar.gz')
    with (out/'source.tar').open('wb') as f:
        subprocess.run(['git','-C',str(ROOT/'../../LinearKnot'),'archive',report['source_commit']],stdout=f,check=True)
    report['source_archive_sha256']=sha(out/'source.tar')
    scripts=ROOT/'jobs/pseudoloop-context-0.2.0-20261023'; scripts.mkdir(parents=True,exist_ok=False)
    (out/'scheduler_logs').mkdir()
    for tier,settings in report['tiers'].items():
        tasks=[t for t in report['tasks'] if t['tier']==tier]
        sec=settings['h_rt_seconds'];wall=f'{sec//3600:02d}:{sec%3600//60:02d}:{sec%60:02d}'
        script=scripts/(tier+'.sh')
        script.write_text('#!/usr/bin/env bash\n#$ -cwd\n#$ -S /bin/bash\n#$ -j y\n'
            f'#$ -N pl020_{tier}\n#$ -l cpu_16=1\n#$ -l h_rt={wall}\n'
            f'#$ -t 1-{len(tasks)}\n#$ -tc 30\n#$ -o {out}/scheduler_logs\n\n'
            'set -euo pipefail\n'+f'cd {shlex.quote(str(ROOT))}\n'
            'export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 LC_ALL=C\n'
            ': "${TMPDIR:?node-local TMPDIR is required}"\n'
            'exec "${LCP_PYTHON:-python3}" -B scripts/pseudoloop_context_array_20261023.py run '
            f'--plan {shlex.quote(str(out/"array_plan.json"))} --tier {tier} --task-id "${{SGE_TASK_ID}}" "$@"\n')
        settings['script']=str(script.relative_to(ROOT));settings['script_sha256']=sha(script)
    dump(out/'array_plan.json',report)
    (out/'array_plan.sha256').write_text(sha(out/'array_plan.json')+'\n')
    print(json.dumps({k:report[k] for k in ['profile_jobs','estimated_allocation_hours',
        'requested_allocation_hours','estimated_points','points_at_h_rt','predicted_raw_profile_stdout_bytes']},indent=2))
    check(argparse.Namespace(plan=out/'array_plan.json'))


def load(path):
    path=path.resolve()
    assert sha(path)==path.with_suffix('.sha256').read_text().strip(), 'plan changed'
    p=json.loads(path.read_text())
    assert sha(__file__)==p['script_sha256'], 'runner changed: replan before execution'
    assert sha(ROOT/'scripts/task_storage.py')==p['task_storage_sha256']
    return p


def check(args):
    path=args.plan.resolve();p=load(path);out=path.parent
    assert sha(out/'runtime.tar.gz')==p['runtime_sha256']
    assert sha(out/'source.tar')==p['source_archive_sha256']
    for r in p['old_plans'].values(): assert sha(ROOT/r['path'])==r['sha256']
    for r in p['input_archives'].values(): assert sha(ROOT/r['path'])==r['sha256']
    assigned=[n for t in p['tasks'] for n in t['inputs']]
    assert len(assigned)==len(set(assigned))==len(p['inputs'])
    for t in p['tasks']:
        s=p['tiers'][t['tier']]
        assert s['workers']*s['memory_limit_mib']+s['task_overhead_mib']<=s['memory_budget_mib']
        assert t['estimated_wall_seconds']<=.65*s['h_rt_seconds']
        assert sha(ROOT/s['script'])==s['script_sha256']
    # Stream frozen inputs and verify every byte hash, sequence hash and length.
    for cohort,r in p['input_archives'].items():
        expected={v['member']:v for v in p['inputs'].values() if v['cohort']==cohort};seen=set()
        with tarfile.open(ROOT/r['path'],'r|') as tf:
            for m in tf:
                if m.name not in expected: continue
                data=tf.extractfile(m).read();v=expected[m.name];seq=data.decode().splitlines()[1]
                assert hashlib.sha256(data).hexdigest()==v['sha256']
                assert len(seq)==v['length'] and hashlib.sha256(seq.encode()).hexdigest()==v['sequence_sha256']
                seen.add(m.name)
        assert seen==set(expected)
    quota=check_budget(estimate_files(len(p['tasks']),len(p['tiers'])),{'file_budget':{'hard_file_limit':24000000}})
    dump(out/'dry_run.json',{'checked_utc':now(),'status':'passed','cli_executed':False,
        'inputs':len(assigned),'profile_jobs':p['profile_jobs'],'tasks':len(p['tasks']),
        'checks':'plan/runtime/source/scripts/storage/old-plan hashes, every input hash/length/sequence, complete disjoint cohort, memory/wall budgets',
        'file_budget':quota})
    print(f'dry-run passed: {len(assigned)} inputs, {p["profile_jobs"]} profiles, {len(p["tasks"])} tasks')


def run(args):
    path=args.plan.resolve();p=load(path);home=path.parent
    task=next(t for t in p['tasks'] if t['tier']==args.tier and t['task_id']==args.task_id)
    settings=p['tiers'][args.tier]; names=list(task['inputs'])
    if args.limit_inputs:
        if not args.out: raise ValueError('a local pilot requires a distinct --out')
        names=names[:args.limit_inputs]
    out=args.out.resolve() if args.out else home/'tasks'/args.tier/f'{args.task_id:06d}'
    out.mkdir(parents=True,exist_ok=True)
    lock=(out/'task.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    if (out/'task.json').exists() and not args.resume: raise ValueError('use --resume for existing output')
    for item in ['runtime.tar.gz']:
        assert sha(home/item)==p['runtime_sha256']
    meta={'started_utc':now(),'status':'running','tier':args.tier,'task_id':args.task_id,
        'inputs':names,'workers':settings['workers'],'requested_resource':'cpu_16=1',
        'h_rt_seconds':settings['h_rt_seconds'],'plan_sha256':sha(path),'runner_sha256':sha(__file__),
        'job_id':os.environ.get('JOB_ID'),'sge_task_id':os.environ.get('SGE_TASK_ID'),
        'hostname':os.uname().nodename,'pilot':bool(args.limit_inputs),'previous_attempts':[]}
    if (out/'task.json').exists():
        prev=json.loads((out/'task.json').read_text())
        assert prev['plan_sha256']==meta['plan_sha256'] and prev['inputs']==names
        meta['previous_attempts']=prev['previous_attempts']+[{k:v for k,v in prev.items() if k!='previous_attempts'}]
    dump(out/'task.json',meta)
    store=TaskStore(out,resume=args.resume)
    with local_workspace() as scratch:
        with tarfile.open(home/'runtime.tar.gz') as tf:
            for m in tf:
                assert m.isfile() and Path(m.name).name==m.name
                dest=scratch/m.name;dest.write_bytes(tf.extractfile(m).read())
                dest.chmod(m.mode)
        assert sha(scratch/'linearknot')==p['binary_sha256']
        assert sha(scratch/'rna_Turner04.par')==p['parameter_sha256']
        cohort=settings['cohort'];entry=p['input_archives'][cohort]
        assert sha(ROOT/entry['path'])==entry['sha256']
        members={p['inputs'][name]['member']:name for name in names}; found=set()
        with tarfile.open(ROOT/entry['path'],'r|') as tf:
            for m in tf:
                if m.name not in members:continue
                name=members[m.name];data=tf.extractfile(m).read()
                assert hashlib.sha256(data).hexdigest()==p['inputs'][name]['sha256']
                (scratch/(name+'.txt')).write_bytes(data);found.add(name)
        assert found==set(names)
        modes=p['modes'][cohort]
        jobs=[{'id':name+'-'+mode['name']+'-LinearPKCapR-rep0','input':name,'mode':mode['name'],
               'args':mode['args']} for name in names for mode in modes]
        random.Random(p['seed']).shuffle(jobs)
        assert store.ids <= {j['id'] for j in jobs}
        env={k:v for k,v in os.environ.items() if not k.startswith('CPARTY_')}
        env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',LC_ALL='C',
            CPARTY_FIXED_ENERGY_PARAM_FILE=str(scratch/'rna_Turner04.par'))
        def execute(job):
            folder=scratch/job['id'];folder.mkdir()
            argv=[str(scratch/'linearknot'),'profile','-i',str(scratch/(job['input']+'.txt')),*job['args']]
            command=['/usr/bin/time','-f','%e\t%M','-o',str(folder/'time.tsv'),*argv]
            # Bash ulimit avoids preexec_fn in a threaded controller.
            limit=f'ulimit -v {settings["memory_limit_mib"]*1024}; exec "$@"'
            t=time.monotonic();started=now();status='success';error=None
            with (folder/'stdout.txt').open('wb') as fo,(folder/'stderr.txt').open('wb') as fe:
                proc=subprocess.Popen(['bash','-c',limit,'profile-limit',*command],stdout=fo,stderr=fe,
                    cwd=scratch,env=env,start_new_session=True)
                try:proc.wait(timeout=settings['timeout_seconds'])
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid,signal.SIGKILL);proc.wait();status='timeout';error='per-process timeout'
            if proc.returncode and status=='success':status='failed';error='nonzero exit'
            numeric={}
            if status=='success':
                for line in (folder/'stdout.txt').read_text().splitlines():
                    fields=line.split()
                    if fields and fields[0] in CONTEXTS:numeric[fields[0]]=list(map(float,fields[1:]))
                n=p['inputs'][job['input']]['length']
                if set(numeric)!=set(CONTEXTS) or any(len(v)!=n for v in numeric.values()):
                    status='invalid_output';error='missing or wrong-length context rows'
                elif any(not math.isfinite(x) or x < -1e-9 or x > 1+1e-9 for v in numeric.values() for x in v):
                    status='invalid_output';error='context probability outside [0,1]'
                elif max(abs(sum(numeric[c][i] for c in CONTEXTS)-1) for i in range(n))>1e-9:
                    status='invalid_output';error='context sum differs from one'
            row={k:v for k,v in job.items() if k!='args'}
            row.update(program='LinearPKCapR',subcommand='profile',warmup=False,repetition=0,
                command=argv,started_utc=started,finished_utc=now(),returncode=proc.returncode,
                status=status,error=error,elapsed_seconds=time.monotonic()-t,
                input_sha256=p['inputs'][job['input']]['sha256'],binary_sha256=p['binary_sha256'],
                parameter_sha256=p['parameter_sha256'],time_wall_seconds=None,peak_rss_kib=None)
            if (folder/'time.tsv').exists():
                fields=(folder/'time.tsv').read_text().splitlines()[-1].split('\t')
                if len(fields)==2:row.update(time_wall_seconds=float(fields[0]),peak_rss_kib=int(fields[1]))
            dump(folder/'command.json',argv)
            return row,folder
        with ThreadPoolExecutor(max_workers=settings['workers']) as pool:
            futures=[pool.submit(execute,j) for j in jobs if j['id'] not in store.ids]
            for f in as_completed(futures):
                row,folder=f.result();files={row['id']+'/'+q.name:q for q in folder.iterdir()}
                store.commit(row,files);shutil.rmtree(folder)
        meta.update(finished_utc=now(),completed_jobs=len(store.rows),
            status='success' if all(r['status']=='success' for r in store.rows) else 'completed_with_failures',
            status_counts=dict(Counter(r['status'] for r in store.rows)))
        assert sha(scratch/'linearknot')==p['binary_sha256']
        dump(out/'task.json',meta)
    print(json.dumps(meta,indent=2))
    if meta['status'] != 'success':
        raise SystemExit(1)


def main():
    ap=argparse.ArgumentParser(description=__doc__);sub=ap.add_subparsers(dest='action',required=True)
    pp=sub.add_parser('plan');pp.add_argument('--source',type=Path,required=True);pp.add_argument('--build',type=Path,required=True)
    pp.add_argument('--out',type=Path,default=DEFAULT)
    cp=sub.add_parser('check');cp.add_argument('--plan',type=Path,default=DEFAULT/'array_plan.json')
    rp=sub.add_parser('run');rp.add_argument('--plan',type=Path,default=DEFAULT/'array_plan.json')
    rp.add_argument('--tier',required=True);rp.add_argument('--task-id',type=int,required=True)
    rp.add_argument('--out',type=Path);rp.add_argument('--limit-inputs',type=int);rp.add_argument('--resume',action='store_true')
    args=ap.parse_args();globals()[args.action](args)


if __name__=='__main__':main()
