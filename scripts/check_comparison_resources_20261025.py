#!/usr/bin/env python3
"""Read-only full dry-run and retained-pilot verification; never submits jobs."""
from collections import Counter,defaultdict
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile

from analyze_comparison_resources_20261025 import classify,collect,summarize
from comparison_resources_array_20261025 import check_bundle
from file_budget import check_budget
from plan_array import load_plan,wall_seconds
from run_cli import ROOT,digest,now,write_json
from task_storage import ArchiveReader,jsonl_rows


def completed(out):
    """Validate terminal tasks against frozen schedules and every retained member."""
    report=json.loads((ROOT/'docs/comparison_resources_20261025_plan.json').read_text())
    inputs_path=ROOT/'data/comparison-resources-20261025/manifest.json'
    assert digest(inputs_path)==report['input_manifest_sha256']
    inputs=json.loads(inputs_path.read_text())
    names={r['name']:r for r in inputs['records']}
    reader=ArchiveReader(ROOT/'data/comparison-resources-20261025/inputs.tar')
    for r in names.values():
        raw=reader.read(r['input_member'])
        assert hashlib.sha256(raw).hexdigest()==r['sha256']
        lines=[s.strip() for s in raw.decode().splitlines() if s.strip()]
        assert len(lines[1])==len(lines[2])==r['length']
        for key,value in [('sequence',lines[1]),('scaffold',lines[2])]:
            assert hashlib.sha256(value.encode()).hexdigest()==r[key+'_sha256']
    tasks=[];allrows=[];snapshots=[];members=0;assigned=defaultdict(list)
    for run in sorted({r['run'] for r in report['arrays']}):
        plan_path=ROOT/'results/runs'/run/'array_plan.json'
        plan=load_plan(plan_path);bundle=check_bundle(plan)
        snapshots.append(dict(path=str(plan_path.relative_to(ROOT)),sha256=digest(plan_path)))
        for script in plan['job_scripts']:
            assert digest(Path(script['path']))==script['sha256']
            subprocess.run(['bash','-n',script['path']],check=True)
        for task in plan['tasks']:
            d=Path(plan['output_dir'])/task['directory'];rd=d/'run'
            meta=json.loads((d/'task.json').read_text());saved=json.loads((rd/'run.json').read_text())
            assert meta['status'] in ('success','failed') and meta['finished_utc']
            assert saved['status']==meta['status'] and saved['finished_utc']
            assert saved['software_lock_valid_after_run'] is True
            assert saved['software_changed_during_run'] is False and saved['binaries_changed_during_run'] is False
            assert saved['software_before']['commit']==plan['software_lock']['source_commit']
            assert not saved['software_before']['dirty'] and not saved['software_after']['dirty']
            assert saved['binaries_before']==saved['binaries_after']
            assert saved['software_lock_sha256']==plan['software_lock_sha256']
            assert saved['workers']==task['workers']
            schedule=json.loads(Path(task['planned_jobs_path']).read_text())
            checkpoint=rd/'task_results.jsonl';blob=checkpoint.read_bytes()
            assert blob.endswith(b'\n')
            rows=list(jsonl_rows(checkpoint));assert len(rows)==len(schedule)
            assert len({r['id'] for r in rows})==len(rows)
            expected={r['id']:r for r in schedule};archive=ArchiveReader(rd/'task_outputs.tar')
            for row in rows:
                assert {k:row[k] for k in expected[row['id']]}==expected[row['id']]
                for member,record in row['artifacts'].items():
                    raw=archive.read(member)
                    assert len(raw)==record['bytes'] and hashlib.sha256(raw).hexdigest()==record['sha256']
                    members+=1
                for stream in ('stdout','stderr'):
                    assert row[stream+'_summary']['sha256']==row['artifacts'][row['directory']+'/'+stream+'.txt']['sha256']
                row['resource_status']=classify(row,archive.read(row['directory']+'/stdout.txt').decode(errors='replace'))
                allrows.append(row)
            for source in saved['inputs']:
                assert source['sha256']==names[source['name']]['sha256']
            grouped=defaultdict(list)
            for row in rows:grouped[row['input'],row['mode']].append(row)
            for (name,method),rs in grouped.items():
                assert len(rs)==4 and sum(r['warmup'] for r in rs)==1
                assert sorted(r['repetition'] for r in rs if not r['warmup'])==[0,1,2]
                assigned[method].append(name)
            tasks.append(dict(run=run,tier=task['tier'],task_id=task['task_id'],status=meta['status'],
                              host=saved['host']['hostname'],workers=task['workers'],limits=task['job_limits'],
                              calls=len(rows),statuses=dict(Counter(r['resource_status'] for r in rows))))
            for path in (d/'task.json',rd/'run.json',checkpoint,rd/'task_outputs.tar'):
                snapshots.append(dict(path=str(path.relative_to(ROOT)),bytes=path.stat().st_size,sha256=digest(path)))
    assert len(tasks)==63 and len(allrows)==2032
    for method,ns in assigned.items():
        expected={n for n,r in names.items() if method!='exact' or r['length']<=2200}
        assert len(ns)==len(set(ns)) and set(ns)==expected
    with tarfile.open(bundle) as archive:
        for name,sha in plan['software_lock']['build_files'].items():
            assert hashlib.sha256(archive.extractfile('runtime/'+name).read()).hexdigest()==sha
    result=dict(status='passed',checked_utc=now(),tasks=len(tasks),arrays=len(report['arrays']),
                cli_calls=len(allrows),measurements=sum(not r['warmup'] for r in allrows),raw_member_hashes_verified=members,
                task_statuses=dict(Counter(t['status'] for t in tasks)),
                measurement_statuses=dict(Counter(r['resource_status'] for r in allrows if not r['warmup'])),
                warmup_statuses=dict(Counter(r['resource_status'] for r in allrows if r['warmup'])),
                method_input_counts={m:len(ns) for m,ns in assigned.items()},task_details=tasks,snapshots=snapshots,
                checks='frozen plans/configs/schedules, inputs S/G, scripts, software lock/bundle/runtime hashes; terminal tasks; complete unique 1+3 schedules; every raw member size/hash; unchanged source/binaries',
                extrapolation=False)
    write_json(out,result)
    print(json.dumps({k:v for k,v in result.items() if k not in ('task_details','snapshots')},indent=2))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--completed',action='store_true')
    parser.add_argument('--out',type=Path,default=ROOT/'docs/comparison_resources_20261026_completed_validation.json')
    args=parser.parse_args()
    if args.completed:
        completed(args.out)
        return
    report=json.loads((ROOT/'docs/comparison_resources_20261025_plan.json').read_text())
    inputs=json.loads((ROOT/'data/comparison-resources-20261025/manifest.json').read_text())
    assert digest(ROOT/'data/comparison-resources-20261025/manifest.json')==report['input_manifest_sha256']
    names={r['name']:r for r in inputs['records']}
    assigned=defaultdict(list);ntasks=njobs=0;scripts=[]
    for run in sorted({r['run'] for r in report['arrays']}):
        plan=load_plan(ROOT/'results/runs'/run/'array_plan.json')
        bundle=check_bundle(plan)
        for script in plan['job_scripts']:
            assert digest(Path(script['path']))==script['sha256']
            subprocess.run(['bash','-n',script['path']],check=True)
            scripts.append(script['path'])
        for task in plan['tasks']:
            ntasks+=1
            d=Path(plan['output_dir'])/task['directory']
            meta=json.loads((d/'task.json').read_text())
            assert meta['status']=='not_started' and not meta['attempts'] and not (d/'run').exists()
            assert task['memory_cap_envelope_mib']<=35095
            assert task['estimated_seconds']<=.8*wall_seconds(task['h_rt'])+1e-6
            jobs=json.loads(Path(task['planned_jobs_path']).read_text());njobs+=len(jobs)
            grouped=defaultdict(list)
            for j in jobs:grouped[j['input'],j['mode']].append(j)
            for (name,method),rs in grouped.items():
                assert len(rs)==4 and sum(r['warmup'] for r in rs)==1
                assert sorted(r['repetition'] for r in rs if not r['warmup'])==[0,1,2]
                assigned[method].append(name)
        lock=plan['software_lock']
        assert digest(ROOT/'manifests/software-comparison-resources-v0.2.0-20261025.json')==plan['software_lock_sha256']
    assert ntasks==report['tasks']==63 and njobs==report['cli_calls']==2032 and len(scripts)==15
    for method,ns in assigned.items():
        expected={n for n,r in names.items() if method!='exact' or r['length']<=2200}
        assert len(ns)==len(set(ns)) and set(ns)==expected
    input_reader=ArchiveReader(ROOT/'data/comparison-resources-20261025/inputs.tar')
    for r in names.values():
        data=input_reader.read(r['input_member']);assert hashlib.sha256(data).hexdigest()==r['sha256']
    runtime_files=[]
    with tarfile.open(bundle) as archive:
        for name,sha in lock['build_files'].items():
            data=archive.extractfile('runtime/'+name).read();assert hashlib.sha256(data).hexdigest()==sha
            runtime_files.append(name)
    pilot=ROOT/'results/local/comparison-resources-pilot-20261025'
    allrows=[];parity=[];raw_count=0
    for phase in ('short','medium'):
        directory=pilot/phase;rs=collect([directory]);allrows.extend(rs)
        reader=ArchiveReader(directory/'task_outputs.tar');values=defaultdict(dict)
        for r in rs:
            for name,meta in r['artifacts'].items():
                assert hashlib.sha256(reader.read(name)).hexdigest()==meta['sha256'];raw_count+=1
            if r['resource_status']=='success' and not r['warmup'] and r['mode'] in ('cparty','exact'):
                out=reader.read(r['directory']+'/stdout.txt').decode()
                values[r['input'],r['repetition']][r['mode']]=next(json.loads(s) for s in out.splitlines() if s.startswith('{'))
        for (name,rep),v in values.items():
            if set(v)=={'cparty','exact'}:
                dz=abs(v['cparty']['log_z']-v['exact']['log_z']);df=abs(v['cparty']['free_energy']-v['exact']['free_energy'])
                assert dz<=1e-9 and df<=1e-10
                parity.append(dict(input=name,repetition=rep,log_z_absolute_error=dz,free_energy_absolute_error=df))
    assert len(allrows)==120 and Counter(r['resource_status'] for r in allrows)=={'success':116,'nonfinite':4}
    assert len(parity)==15
    resources=summarize(allrows,[(r['input'],r['mode']) for r in allrows])
    for r in resources:r.update(length=names[r['input']]['length'],group=names[r['input']]['group'])
    staging=ROOT/'results/runs/comparison-resources-staging-smoke-20261025/tasks/pilot/000001/run'
    srows=collect([staging]);assert len(srows)==40 and all(r['resource_status']=='success' for r in srows)
    resume=json.loads((pilot/'resume_verification.json').read_text());assert resume['status']=='passed'
    controls=collect([pilot/'controls'])
    assert Counter(r['resource_status'] for r in controls)=={'timeout':2,'memory_limit':1}
    quota=check_budget(dict(estimated_shared_entries=report['total_shared_entry_budget'],basis='all preparation, 63 tasks, first scheduler attempt and atomic-write headroom'),{})
    verification=dict(status='passed',checked_utc=now(),formal_tasks=ntasks,formal_tasks_started=0,
                      arrays=len(scripts),cli_calls=njobs,measurements=njobs*3//4,
                      input_count=len(names),method_input_counts={k:len(v) for k,v in assigned.items()},
                      frozen_runtime_files_verified=runtime_files,file_budget=quota,
                      checks='frozen plans, scripts, locks, bundle and runtime hashes; every input hash; complete/disjoint method schedules; process memory, timeout and wall budgets; bash -n; all formal tasks unstarted',
                      short_medium_calls=120,short_medium_statuses=dict(Counter(r['resource_status'] for r in allrows)),
                      raw_pilot_member_hashes_verified=raw_count,staging_calls=40,resume=resume,
                      limit_controls=[dict(id=r['id'],status=r['resource_status'],elapsed_seconds=r['elapsed_seconds'],
                                           peak_rss_kib=r['peak_rss_kib']) for r in controls],
                      pf_parity_checks=len(parity),pf_log_z_max_absolute_error=max(r['log_z_absolute_error'] for r in parity),
                      pf_free_energy_max_absolute_error=max(r['free_energy_absolute_error'] for r in parity),
                      qsub_executed=False)
    write_json(ROOT/'docs/comparison_resources_20261025_validation.json',verification)
    write_json(pilot/'summary.json',dict(resources=resources,pf_parity=parity,verification=verification))
    print(json.dumps(verification,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
