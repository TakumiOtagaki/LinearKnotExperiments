#!/usr/bin/env python3
"""Small development folds, native-CLI parity, timeout packing and resume checks."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile

from analyze_beam_sweep_20261027 import analyze
from beam_sweep_array_20261027 import GzipBppStore
from plan_array import prepare
from prepare_beam_sweep_20261027 import BASE, BEAMS, E2, config
from run_cli import ROOT, digest, run_one, write_json
from task_storage import ArchiveReader, MemoryFile, jsonl_rows
from rna_records import read_record


def finish_validation(out):
    pp=out/'run/array_plan.json'
    plan=json.loads(pp.read_text())
    reader=ArchiveReader(out/'native_parity.tar.gz')
    parity=[]
    for task in plan['tasks']:
        for row in jsonl_rows(out/'run'/task['directory']/'run/task_results.jsonl'):
            prefix=row['id']+'/'
            if prefix+'bpp.tsv' not in reader.members:
                continue
            stdout=next(v for n,v in row['artifacts'].items() if n.endswith('/stdout.txt'))
            assert hashlib.sha256(reader.read(prefix+'stdout.txt')).hexdigest()==stdout['sha256']
            assert hashlib.sha256(reader.read(prefix+'bpp.tsv')).hexdigest()==row['bpp_uncompressed_sha256']
            timing=reader.read(prefix+'time.tsv').decode().split()
            parity.append(dict(input=row['input'],mode=row['mode'],native_wall_s=float(timing[0]),
                               stdout_identical=True,bpp_identical=True))
    assert len(parity)==5
    timeout=jsonl_rows(out/'timeout/task_results.jsonl')[0]
    assert timeout['status']=='timeout'
    summary=analyze([pp],out/'analysis',Path('/tmp'))
    validation=dict(parity=parity,resume_rows_unchanged=True,
                    forced_timeout={k:timeout[k] for k in ('status','elapsed_seconds','time_wall_seconds','peak_rss_kib')},summary=summary)
    write_json(out/'validation.json',validation)
    print(json.dumps(validation,indent=2))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--collect-only',action='store_true')
    args=parser.parse_args()
    out=ROOT/'results/local'/BASE/'smoke'
    if args.collect_only:
        finish_validation(out)
        return
    out.mkdir(exist_ok=False)
    cfgdir=ROOT/'configs'/BASE
    old=json.loads(E2.read_text())
    item=next(r for r in old['inputs'] if r['name']=='URS000233CDDE')
    real={k:item[k] for k in ('name','path','length','input_archive','input_member','sha256')}
    small=dict(name='pk_scaffold',path='data/smoke/pk_scaffold.txt',
               length=len(read_record(ROOT/'data/smoke/pk_scaffold.txt')[1]),sha256=digest(ROOT/'data/smoke/pk_scaffold.txt'))
    cfg=config(BASE+'-smoke',[small,real],25,purpose='development')
    cfg['modes']=[]
    for b in BEAMS:
        mode=config('x',[small],b)['modes'][0]
        if b>50:
            mode['length_range']=dict(min_length=1,max_length=100)
        cfg['modes'].append(mode)
    cfg.update(timeout_seconds=180,memory_limit_mib=4096)
    tiers=[]
    for name,lo,hi in (('tiny',1,100),('short_rna',4793,4793)):
        tiers.append(dict(name=name,min_length=lo,max_length=hi,resource='cpu_16',h_rt='00:10:00',
                          target_seconds=420,max_inputs_per_task=1,workers=1,memory_budget_mib=35095,task_overhead_mib=1024,
                          cost_model=dict(basis='Development smoke only, 2 s per tiny or .006 s per nt condition; not formal timing',
                                          seconds_constant_per_condition=2,seconds_per_nt_condition=.006,memory_base_mib=4096,mib_per_nt=0,safety_factor=1)))
    settings=dict(schema_version=1,measurement='accuracy',seed=20261027,compact_plan=True,tiers=tiers,
                  task_runner_script='scripts/beam_sweep_array_20261027.py')
    cp,ap=cfgdir/'smoke.json',cfgdir/'smoke-array.json'
    write_json(cp,cfg);write_json(ap,settings)
    plan=prepare(cp,ap,out/'run',ROOT/'jobs'/(BASE+'-smoke'),dry_run=True)
    pp=out/'run/array_plan.json'
    commands=[]
    # This interactive allocation is not an array; SGE exports "undefined".
    child_env=dict(os.environ)
    child_env.pop('SGE_TASK_ID',None)
    for task in plan['tasks']:
        command=[sys.executable,'scripts/beam_sweep_array_20261027.py','--plan',str(pp),
                 '--tier',task['tier'],'--task-id',str(task['task_id'])]
        commands.append(command)
        subprocess.run(command,check=True,env=child_env)
    before={t['tier']:(out/'run'/t['directory']/'run/task_results.jsonl').read_bytes() for t in plan['tasks']}
    for command in commands:
        subprocess.run(command+['--resume'],check=True,env=child_env)
    assert all(before[t['tier']]==(out/'run'/t['directory']/'run/task_results.jsonl').read_bytes() for t in plan['tasks'])
    parity=[]
    with tempfile.TemporaryDirectory(prefix='beam-sweep-parity-',dir='/tmp') as temporary:
        work=Path(temporary)
        lock=json.loads((ROOT/'manifests'/('software-'+BASE+'.json')).read_text())
        with tarfile.open(ROOT/lock['software_bundle']) as tar:
            tar.extractall(work)
        parameter=work/'source/params/rna_Turner04.par'
        env=dict(os.environ,CPARTY_FIXED_ENERGY_PARAM_FILE=str(parameter),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',LC_ALL='C')
        for task in plan['tasks']:
            run=out/'run'/task['directory']/'run'
            reader=ArchiveReader(run/'task_outputs.tar')
            for row in jsonl_rows(run/'task_results.jsonl'):
                if row['input']=='pk_scaffold' and not any(x in row['mode'] for x in ('beam25_','beam50_','beam400_')):
                    continue
                directory=work/row['id'];directory.mkdir()
                if row['input']=='pk_scaffold':
                    data=(ROOT/small['path']).read_bytes()
                else:
                    data=ArchiveReader(Path(real['input_archive'])).read(real['input_member'])
                input_file=directory/'input.txt';input_file.write_bytes(data)
                command=list(row['command']);command[0]=str(work/'runtime/linearknot')
                command[command.index('-i')+1]=str(input_file)
                command[command.index('--bpp')+1]=str(directory/'bpp.tsv')
                result=run_one(command,directory,env,180,4096)
                assert result['status']=='success'
                stdout=next(v for n,v in row['artifacts'].items() if n.endswith('/stdout.txt'))
                assert digest(directory/'stdout.txt')==stdout['sha256']
                assert digest(directory/'bpp.tsv')==row['bpp_uncompressed_sha256']
                parity.append(dict(input=row['input'],mode=row['mode'],native_wall_s=result['time_wall_seconds'],
                                   stdout_identical=True,bpp_identical=True))
        timeout_dir=work/'forced-timeout';timeout_dir.mkdir()
        input_file=timeout_dir/'input.txt'
        input_file.write_bytes(ArchiveReader(Path(real['input_archive'])).read(real['input_member']))
        command=[str(work/'runtime/FoldMetadata'),'fold','-i',str(input_file),'-b','400',
                 '--max-hairpin','400','--max-case3-depth','400','--ranking','outer','--bpp-cutoff','0','--bpp',str(timeout_dir/'bpp.tsv')]
        timeout=run_one(command,timeout_dir,env,.05,4096)
        assert timeout['status']=='timeout'
        (out/'timeout').mkdir()
        store=GzipBppStore(out/'timeout')
        store.commit(dict(timeout,id='forced-timeout'),{'forced-timeout/'+p.name:p for p in timeout_dir.iterdir()})
        assert GzipBppStore(out/'timeout',resume=True).rows[0]['status']=='timeout'
        with tarfile.open(out/'native_parity.tar.gz','w:gz') as tar:
            for task in plan['tasks']:
                for row in jsonl_rows(out/'run'/task['directory']/'run/task_results.jsonl'):
                    path=work/row['id']
                    if path.exists():
                        tar.add(path,arcname=row['id'])
    finish_validation(out)


if __name__=='__main__':
    main()
