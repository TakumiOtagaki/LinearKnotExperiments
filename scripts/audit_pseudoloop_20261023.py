#!/usr/bin/env python3
"""Audit pilots, the published error witness, and saved exact Terc samples."""
import argparse
from collections import Counter
import hashlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import time

from rerun_pseudoloop_small_20261023 import access, dump, invariants, sha, CONTEXTS
from task_storage import ArchiveReader

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/local/pseudoloop-0.2.0-20261023/audit'


def profile(text):
    d={}
    for line in text.splitlines():
        f=line.split()
        if f and f[0] in CONTEXTS:d[f[0]]=list(map(float,f[1:]))
    assert set(d)==set(CONTEXTS)
    return d


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source',type=Path,required=True);ap.add_argument('--build',type=Path,required=True)
    args=ap.parse_args();OUT.mkdir(parents=True,exist_ok=False)
    planpath=ROOT/'results/runs/pseudoloop-context-0.2.0-20261023/array_plan.json'
    plan=json.loads(planpath.read_text());archives={};reports={};calls=[]
    def read(path,member,expected):
        path=ROOT/path
        if path not in archives:archives[path]=ArchiveReader(path)
        data=archives[path].read(member)
        assert hashlib.sha256(data).hexdigest()==expected
        return data
    for pilot in sorted((OUT.parent/'pilots').iterdir()):
        rows=[json.loads(s) for s in (pilot/'task_results.jsonl').read_text().splitlines()]
        local=ArchiveReader(pilot/'task_outputs.tar')
        tests=[]
        for row in rows:
            assert row['status']=='success'
            old=next(r for r in plan['inputs'][row['input']]['observations'] if r['mode']==row['mode'])
            olddata=read(old['archive'],old['member'],old['stdout_sha256'])
            member=row['id']+'/stdout.txt';newdata=local.read(member)
            assert hashlib.sha256(newdata).hexdigest()==row['artifacts'][member]['sha256']
            test=invariants(profile(olddata.decode()),profile(newdata.decode()),CONTEXTS,'Multiloop','PseudoLoop')
            test.update(input=row['input'],mode=row['mode'],old_seconds=old['elapsed_seconds'],
                new_seconds=row['elapsed_seconds'],ratio=row['elapsed_seconds']/old['elapsed_seconds'],peak_rss_kib=row['peak_rss_kib'])
            tests.append(test)
        reports[pilot.name]=tests
    env={k:v for k,v in os.environ.items() if not k.startswith('CPARTY_')}
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',LC_ALL='C',
               CPARTY_FIXED_ENERGY_PARAM_FILE=str(args.source/'params/rna_Turner04.par'))
    scratch=args.source.parent/'small-audit';scratch.mkdir(exist_ok=False)
    with tarfile.open(OUT/'raw_outputs.tar.gz','w:gz') as archive:
        def run(argv,name):
            t=time.monotonic();p=subprocess.run(list(map(str,argv)),env=env,capture_output=True,timeout=600)
            for channel,data in [('stdout',p.stdout),('stderr',p.stderr)]:
                tmp=scratch/(name+'.'+channel);tmp.write_bytes(data);archive.add(tmp,arcname=tmp.name)
            calls.append(dict(command=list(map(str,argv)),returncode=p.returncode,
                elapsed_seconds=time.monotonic()-t,stdout_sha256=hashlib.sha256(p.stdout).hexdigest(),
                stderr_sha256=hashlib.sha256(p.stderr).hexdigest()))
            assert p.returncode==0,p.stderr.decode()
            return p.stdout.decode()
        # Recompute the named 335-nt SI example now: one full small matrix.
        name='bpRNA_RFAM_35940';item=plan['inputs'][name]
        with tarfile.open(ROOT/plan['input_archives']['short']['path']) as tf:
            data=tf.extractfile(item['member']).read()
        assert hashlib.sha256(data).hexdigest()==item['sha256']
        inp=scratch/(name+'.txt');inp.write_bytes(data);archive.add(inp,arcname=inp.name)
        witness={};newprofiles={}
        for m in plan['modes']['short']:
            text=run([args.build/'linearknot','profile','-i',inp,*m['args']],m['name'])
            new=profile(text);newprofiles[m['name']]=new
            old=next(r for r in item['observations'] if r['mode']==m['name'])
            previous=profile(read(old['archive'],old['member'],old['stdout_sha256']).decode())
            witness[m['name']]=dict(invariants(previous,new,CONTEXTS,'Multiloop','PseudoLoop'),
                old_base325={c:previous[c][324] for c in CONTEXTS},new_base325={c:new[c][324] for c in CONTEXTS})
        for mode in witness:
            witness[mode]['new_base325_pseudoloop_abs_error']=abs(newprofiles[mode]['PseudoLoop'][324]-newprofiles['exact']['PseudoLoop'][324])
        reports['bpRNA_RFAM_35940']=witness
        # The 61-nt saved figure had total-only accessibility: obtain its old
        # by-loop output from the already certified, unchanged old executable.
        fig=json.loads((OUT.parent/'examples/data.json').read_text())['figure1']
        inp=scratch/'figure1.txt';inp.write_text('>'+fig['record']['id']+'\n'+fig['record']['seq']+'\n'+fig['scaffold']+'\n')
        oldbinary=ROOT/'build/linearknot-9c7d6e1/linearknot'
        oldtext=run([oldbinary,'access','-i',inp,'--exact','--access-len','1,4,8','--probability','--by-loop'],'figure1-old-byloop')
        newtext=run([args.build/'linearknot','access','-i',inp,'--exact','--access-len','1,4,8','--probability','--by-loop'],'figure1-new-byloop')
        oldp=scratch/'old';newp=scratch/'new';oldp.write_text(oldtext);newp.write_text(newtext)
        oldacc,newacc=access(oldp),access(newp)
        reports['figure1_byloop']=invariants({c:[v[c] for _,v in sorted(oldacc.items())] for c in 'HBIMEP'},
            {c:[v[c] for _,v in sorted(newacc.items())] for c in 'HBIMEP'},'HBIMEP','M','P')
        reports['figure1_byloop']['old_binary_sha256']=sha(oldbinary)
        # Reclassify the original sample and the adopted saved rerun separately.
        sys.path.insert(0,str(args.source/'test/tools'))
        from pseudoloop_context_reference import classify,pairs_from_db
        sample_sources={}
        early=ROOT/'results/local/figure4-terc-20261017/exact.samp'
        sample_sources['original_20261017']=early.read_bytes()
        with tarfile.open(OUT.parent/'examples/terc_data.tar') as tf:
            sample_sources['adopted_89baea52']=tf.extractfile('exact.samp').read()
        samples={}
        for label,data in sample_sources.items():
            counts=Counter();total=0
            for line in data.decode().splitlines():
                f=line.split('\t')
                if len(f)!=2 or not f[0].isdigit():continue
                count,db=int(f[0]),f[1];total+=count
                if db[122:126]!='....':continue
                labels,_=classify(len(db),pairs_from_db(db));window=labels[122:126]
                assert len(set(window))==1,window
                counts[window[0]]+=count
            assert total==1000
            samples[label]={'sha256':hashlib.sha256(data).hexdigest(),'draws':total,
                'unpaired_window_draws':sum(counts.values()),'context_counts':dict(counts)}
        reports['terc_exact_samples']=samples
    dump(OUT/'summary.json',reports)
    dump(OUT/'provenance.json',dict(plan_sha256=sha(planpath),binary_sha256=sha(args.build/'linearknot'),
        parameter_sha256=sha(args.source/'params/rna_Turner04.par'),script_sha256=sha(__file__),
        calls=calls,raw_archive_sha256=sha(OUT/'raw_outputs.tar.gz')))
    print(json.dumps(reports['terc_exact_samples'],indent=2))


if __name__=='__main__':main()
