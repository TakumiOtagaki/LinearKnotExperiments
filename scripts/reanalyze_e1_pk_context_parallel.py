"""Task-parallel adapter for analyze_e1_pk_stratified's combined exact/beam run.

All metric, selection, quantile and witness definitions come from the adopted
scripts. Reproduce every old stratified TSV cell as an additional check before
publishing the new-context tables. Read each packed task only once.
"""
from array import array
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
import csv
import json
import math
from pathlib import Path
import pickle
import sqlite3
import time
import zlib

import analyze_accuracy as accuracy
import analyze_e1_pk_stratified as pk
import analyze_e1_production as production
from analyze_bprna_short import CHECKS
from run_cli import digest, now

RUN = PLAN = DB = RAW = None


def initialize(run, db):
    global RUN, PLAN, DB
    RUN = Path(run)
    PLAN = production.load_plan(RUN)
    DB = sqlite3.connect(f'file:{db}?mode=ro', uri=True)


def with_new_profile(parsed, name, mode):
    result = dict(parsed)
    blob = DB.execute('SELECT value FROM profiles WHERE name=? AND mode=?',(name,mode)).fetchone()[0]
    result['LinearPKCapR'] = pickle.loads(zlib.decompress(blob))
    assert accuracy.check_consistency(name,mode,0,result)['status'] == 'success'
    return result


def process(task):
    checks = dict(tasks=[],selected_job_status_counts=Counter(),excluded_job_status_counts=Counter(),
                  consistency_checks=Counter(),consistency_maxima=dict.fromkeys(CHECKS,0.))
    modes = ['exact']+[f'beam{b}_outer' for b in pk.BEAMS]
    reader, hashes, jobs, records = production.open_task(RUN,PLAN,task,modes,checks)
    result = []
    for name,record in records.items():
        ref_old,_ = RAW(reader,jobs,name,'exact',record,checks)
        ref_new = with_new_profile(ref_old,name,'exact')
        descriptors_old,selected_old = pk.select_reference(record,ref_old)
        descriptors_new,selected_new = pk.select_reference(record,ref_new)
        beams = []
        for beam in pk.BEAMS:
            mode = f'beam{beam}_outer'
            old,_ = RAW(reader,jobs,name,mode,record,checks)
            new = with_new_profile(old,name,mode)
            row_old = production.compare(old,ref_old,record)
            row_new = dict(row_old)
            contexts = accuracy.compare_contexts(new['LinearPKCapR'],ref_new['LinearPKCapR'],record)
            row_new['context_tv_mean_Gfree'] = contexts['context_tv_mean_Gfree']
            local_old = pk.local_errors(old,selected_old)
            local_new = dict(local_old)
            local_new[pk.LOCAL[0]] = {i:abs(new['LinearPKCapR']['PseudoLoop'][i]-p)
                                    for i,p in selected_new[pk.LOCAL[0]].items()}
            for row,local in ((row_old,local_old),(row_new,local_new)):
                for metric,errors in local.items():
                    row[metric+'_mae'] = accuracy.mean(list(errors.values()))
                    row[metric+'_max'] = max(errors.values()) if errors else None
            beams.append((beam,{k:row_old[k] for k in pk.RNA_METRICS},
                          {k:row_new[k] for k in pk.RNA_METRICS},local_old,local_new))
        result.append((name,hashes[name],len(record[1]),descriptors_old,descriptors_new,
                       selected_old,selected_new,beams))
    return task['tier'],result,checks


class Aggregate:
    def __init__(self):
        self.counts = {g:Counter() for g in pk.GROUPS}
        self.rna = defaultdict(lambda:array('d'))
        self.local = defaultdict(lambda:array('d'))
        self.witnesses = {}

    def add(self,name,length,tier,descriptors,selected,beams,new):
        groups = [('all','all')]+[(a,pk.stratum(descriptors[a])) for a in pk.AXES]
        for group in groups:
            c = self.counts[group]
            c.update(rna=1,**{tier:1})
            for key in pk.LOCAL:
                c[key+'_rna'] += bool(selected[key])
                c[key+'_events'] += len(selected[key])
            c['crossing_pair_sum'] += descriptors[pk.AXES[0]]
            c['pseudoloop_sum'] += descriptors[pk.AXES[1]]
        for b,oldrow,newrow,oldlocal,newlocal in beams:
            row,local = (newrow,newlocal) if new else (oldrow,oldlocal)
            for a,s in groups:
                for metric in pk.RNA_METRICS:
                    self.rna[a,s,b,metric].append(float('nan') if row[metric] is None else row[metric])
                for metric,errors in local.items():
                    key = a,s,b,metric
                    self.local[key].extend(errors.values())
                    if errors:
                        coord,maximum = max(errors.items(),key=lambda x:x[1])
                        if key not in self.witnesses or maximum > self.witnesses[key]['absolute_error']:
                            location = str(coord+1) if metric == pk.LOCAL[0] else (f'{coord[0]}-{coord[1]}'
                                if metric == pk.LOCAL[1] else f'{coord[0]}-{coord[0]+coord[1]-1}')
                            self.witnesses[key] = dict(axis=a,stratum=s,beam=b,metric=metric,input=name,
                                coordinate_1based=location,absolute_error=maximum,exact_probability=selected[metric][coord],
                                length=length,**descriptors)

    def tables(self):
        fields = ['axis','stratum','rna']+[t[0] for t in production.SHORT]+[
            'crossing_pair_sum','pseudoloop_sum']+[m+s for m in pk.LOCAL for s in ('_rna','_events')]
        counts = [dict(axis=a,stratum=s,**{k:self.counts[a,s][k] for k in fields[2:]}) for a,s in pk.GROUPS]
        rna,local = [],[]
        for a,s in pk.GROUPS:
            for b in pk.BEAMS:
                c = self.counts[a,s]
                for metric in pk.RNA_METRICS:
                    assert len(self.rna[a,s,b,metric]) == c['rna']
                    q = pk.quantiles(self.rna[a,s,b,metric])
                    rna.append(dict(axis=a,stratum=s,beam=b,metric=metric,rna_total=c['rna'],
                                    missing_rna=c['rna']-q['n'],**q))
                for metric in pk.LOCAL:
                    q = pk.quantiles(self.local[a,s,b,metric])
                    assert q['n'] == c[metric+'_events']
                    local.append(dict(axis=a,stratum=s,beam=b,metric=metric,rna_total=c['rna'],
                                      rna_selected=c[metric+'_rna'],**q))
        return {'stratum_counts.tsv':counts,'rna_error_quantiles.tsv':rna,
                'local_error_quantiles.tsv':local,'max_error_examples.tsv':[
                    self.witnesses[k] for k in sorted(self.witnesses)]}


def run(run,db,previous,old_previous,old_stratified,out,workers,raw_reader):
    global RAW
    RAW = raw_reader
    started = time.monotonic()
    plan = production.load_plan(run)
    tasks,hashes = production.cohort(plan,production.SHORT)
    assert not out.exists()
    old_metrics,_ = pk.read_previous(old_previous,[(run,plan)])
    new_metrics,_ = pk.read_previous(previous,[(run,plan)])
    aggregate_old,aggregate_new = Aggregate(),Aggregate()
    differences_old,differences_new = Counter(),Counter()
    provenance = dict(generated_utc=now(),workers=workers,tasks=[],selected_job_status_counts=Counter(),
        consistency_checks=Counter(),consistency_maxima=dict.fromkeys(CHECKS,0.),
        definition_scripts_sha256={n:digest(Path(__file__).with_name(n)) for n in
            ('analyze_e1_pk_stratified.py','analyze_e1_production.py','analyze_accuracy.py')},
        adapter_sha256=digest(Path(__file__)),previous_metrics_sha256=digest(previous/'metrics.tsv'),
        quantile_method='numpy linear; NA excluded; zeros retained',
        local_selection='exact probability >=0.1, inclusive',
        storage='checksum-verified in-memory task tar; node-local profile database',validation='passed')
    seen = set()
    with ProcessPoolExecutor(workers,initializer=initialize,initargs=(str(run),str(db))) as pool:
        for number,(tier,records,checks) in enumerate(pool.map(process,tasks),1):
            for name,checksum,length,do,dn,so,sn,beams in records:
                assert name not in seen and checksum == hashes[name]
                seen.add(name)
                for b,ro,rn,_,_ in beams:
                    pk.check_previous(old_metrics,name,b,checksum,tier,ro,differences_old)
                    pk.check_previous(new_metrics,name,b,checksum,tier,rn,differences_new)
                aggregate_old.add(name,length,tier,do,so,beams,False)
                aggregate_new.add(name,length,tier,dn,sn,beams,True)
            provenance['tasks'].extend(checks['tasks'])
            for key in ('selected_job_status_counts','consistency_checks'):
                provenance[key].update(checks[key])
            for key in CHECKS:
                provenance['consistency_maxima'][key] = max(provenance['consistency_maxima'][key],checks['consistency_maxima'][key])
            if number % 10 == 0 or number == len(tasks):
                print(f'combined stratified {number}/{len(tasks)} tasks; {len(seen)} RNA; {time.monotonic()-started:.1f}s',flush=True)
    assert seen == set(hashes) and not old_metrics and not new_metrics
    assert provenance['consistency_checks'] == {'success':len(seen)*5}
    reproduced = aggregate_old.tables()
    replication = {}
    for name,actual in reproduced.items():
        with (old_stratified/name).open() as f:
            saved = list(csv.DictReader(f,delimiter='\t'))
        assert len(saved) == len(actual)
        maximum = 0.
        for a,b in zip(saved,actual):
            assert set(a) == set(b)
            for k,v in a.items():
                if v == 'NA':
                    assert b[k] is None
                elif v != str(b[k]):
                    delta = abs(float(v)-float(b[k]))
                    assert math.isclose(float(v),float(b[k]),rel_tol=1e-12,abs_tol=1e-12),(name,k,v,b[k])
                    maximum = max(maximum,delta)
        replication[name] = dict(rows=len(saved),max_absolute_difference=maximum,source_sha256=digest(old_stratified/name))
    out.mkdir()
    for name,rows in aggregate_new.tables().items():
        pk.write_tsv(out/name,list(rows[0]),rows)
    provenance.update(inputs=len(seen),rna_comparisons=len(seen)*4,elapsed_seconds=time.monotonic()-started,
        old_tables_reproduced=replication,old_metric_max_differences=dict(differences_old),
        new_metric_max_differences=dict(differences_new),output_sha256={p.name:digest(p) for p in out.glob('*.tsv')})
    production.write_json(out/'provenance.json',provenance)
