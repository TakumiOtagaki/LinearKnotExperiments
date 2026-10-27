#!/usr/bin/env python3
"""Stream packed beam sweeps through node-local SQLite; keep missing comparisons."""
import argparse
from contextlib import contextmanager
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import sqlite3
import tarfile
import tempfile

from plan_array import load_plan, task_jobs
from prepare_beam_sweep_20261027 import BEAMS
from rna_records import read_record, scaffold_pairs
from run_cli import ROOT, digest, now, write_json
from task_storage import ArchiveReader, MemoryFile, jsonl_rows

TOL = 1e-9


def quantile(values, p):
    if not values:
        return None
    values = sorted(values)
    at = (len(values)-1)*p
    lo = int(at)
    return values[lo]+(values[min(lo+1,len(values)-1)]-values[lo])*(at-lo)


@contextmanager
def member_stream(row, suffix):
    matches = [n for n in row['_result']['artifacts'] if n.endswith('/'+suffix)]
    if len(matches) != 1:
        raise ValueError('missing or ambiguous artifact: '+suffix)
    name = matches[0]
    expected = row['_result']['artifacts'][name]
    with tarfile.open(row['_run']/'task_outputs.tar', 'r:') as tar:
        member = tar.getmember(name)
        if member.size != expected['bytes']:
            raise ValueError('artifact size mismatch')
        with tar.extractfile(member) as raw:
            h = hashlib.sha256()
            for block in iter(lambda: raw.read(1024*1024), b''):
                h.update(block)
            if h.hexdigest() != expected['sha256']:
                raise ValueError('artifact hash mismatch')
        with tar.extractfile(member) as raw:
            yield raw


def scalar_values(row):
    with member_stream(row, 'stderr.txt') as stream:
        lines = stream.read().decode().splitlines()
    metadata = [json.loads(s[len('BEAM_SWEEP_SCALARS '):]) for s in lines if s.startswith('BEAM_SWEEP_SCALARS ')]
    if len(metadata) != 1:
        raise ValueError('missing full precision scalar metadata')
    out = metadata[0]
    if set(out) != {'log_z','ensemble_free_energy','mfe'} or not all(math.isfinite(x) for x in out.values()):
        raise ValueError('invalid scalar metadata')
    with member_stream(row, 'stdout.txt') as stream:
        text = stream.read().decode()
    mfe = re.search(r'^(\S+) \((-?[\d.]+)\)$', text, re.M)
    energy = re.search(r'^ensemble free energy: (-?[\d.]+) kcal/mol$', text, re.M)
    if (not mfe or not energy or len(mfe[1]) != row['n'] or text.splitlines()[0] != '>'+row['name']
            or abs(float(mfe[2])-out['mfe']) > .005001
            or abs(float(energy[1])-out['ensemble_free_energy']) > .000050001):
        raise ValueError('native stdout/scalar metadata disagreement')
    if out['log_z'] and not .61 < -out['ensemble_free_energy']/out['log_z'] < .62:
        raise ValueError('unexpected thermodynamic RT')
    return dict(log_z=out['log_z'], F_kcal_mol=out['ensemble_free_energy'],
                MFE_kcal_mol=out['mfe'], mfe_structure_sha256=hashlib.sha256(mfe[1].encode()).hexdigest())


def import_bpp(db, row, record):
    name, seq, scaffold = record
    fixed = scaffold_pairs(scaffold)
    occupied = {x for pair in fixed for x in pair}
    seen_fixed, headers, incident, batch = set(), [], [0.]*(len(seq)+1), []
    count = 0
    raw_hash = hashlib.sha256()
    with member_stream(row, 'bpp.tsv.gz') as compressed:
        with gzip.GzipFile(fileobj=compressed) as stream:
            for line in stream:
                raw_hash.update(line)
                raw = line.decode().strip()
                if not raw:
                    continue
                if raw.startswith('>'):
                    headers.append(raw[1:]); continue
                fields = raw.split()
                if len(fields) != 3:
                    raise ValueError('BPP requires three columns')
                i, j, p = int(fields[0]), int(fields[1]), float(fields[2])
                if not 1 <= i < j <= len(seq) or not math.isfinite(p) or not -TOL <= p <= 1+TOL:
                    raise ValueError('invalid BPP coordinate/probability')
                p = max(0., min(1., p))
                if (i,j) in fixed:
                    if abs(p-1) > TOL:
                        raise ValueError('scaffold pair mass is not one')
                    seen_fixed.add((i,j))
                elif (i in occupied or j in occupied) and p > TOL:
                    raise ValueError('generated pair touches fixed base')
                incident[i] += p; incident[j] += p
                batch.append((row['beam'],i,j,p)); count += 1
                if len(batch) == 8192:
                    db.executemany('INSERT INTO pairs VALUES(?,?,?,?)',batch); batch=[]
    db.executemany('INSERT INTO pairs VALUES(?,?,?,?)',batch)
    if headers != [name] or seen_fixed != fixed or max(incident) > 1+TOL:
        raise ValueError('BPP header/scaffold/incident mass check failed')
    if raw_hash.hexdigest() != row['_result']['bpp_uncompressed_sha256']:
        raise ValueError('uncompressed BPP hash mismatch')
    db.commit()
    row['bpp_pairs'] = count
    row['max_incident_mass'] = max(incident)
    return [max(0.,1-p) for p in incident]


def tv_compare(db, a, b, qa, qb, scaffold):
    """Merge sparse coordinate union; no dense n by n allocation."""
    first = iter(db.execute('SELECT i,j,p FROM pairs WHERE b=? ORDER BY i,j',(a,)))
    second = iter(db.execute('SELECT i,j,p FROM pairs WHERE b=? ORDER BY i,j',(b,)))
    x, y = next(first,None), next(second,None)
    delta, maximum, max_pair = [0.]*(len(scaffold)+1), 0., None
    while x is not None or y is not None:
        if y is None or (x is not None and x[:2] < y[:2]):
            i,j,p = x; d=abs(p); x=next(first,None)
        elif x is None or y[:2] < x[:2]:
            i,j,p = y; d=abs(p); y=next(second,None)
        else:
            i,j,p = x; d=abs(p-y[2]); x=next(first,None); y=next(second,None)
        delta[i] += d; delta[j] += d
        if d > maximum:
            maximum,max_pair = d,[i,j]
    values = [.5*(delta[i]+abs(qa[i]-qb[i])) for i,c in enumerate(scaffold,1) if c=='.']
    if values and max(values) > 1+TOL:
        raise ValueError('TV exceeds one')
    return dict(tv_mean=sum(values)/len(values) if values else None,
                tv_p90=quantile(values,.9), tv_max=max(values) if values else None,
                max_pair_probability_difference=maximum, max_pair_coordinates=max_pair,
                free_bases=len(values))


def collect(plans):
    rows, records, signatures = {}, {}, set()
    for path in plans:
        plan=load_plan(path.resolve())
        lock=plan['software_lock']
        signatures.add((lock['source_commit'],lock['energy_parameter_sha256'],lock['binaries']['FoldMetadata']))
        for item in plan['inputs']:
            data=ArchiveReader(Path(item['input_archive'])).read(item['input_member'])
            if hashlib.sha256(data).hexdigest()!=item['sha256']:
                raise ValueError('input hash mismatch')
            record=read_record(MemoryFile(data))
            if item['name'] in records and records[item['name']]!=record:
                raise ValueError('different inputs with same name')
            records[item['name']]=record
        modes={m['name']:m for m in plan['config']['modes']}
        for task in plan['tasks']:
            run=Path(plan['output_dir'])/task['directory']/'run'
            saved=jsonl_rows(run/'task_results.jsonl')
            byid={r['id']:r for r in saved}
            schedule=task_jobs(plan,task)
            if len(byid)!=len(saved) or set(byid)-{j['id'] for j in schedule}:
                raise ValueError('duplicate or unplanned observations')
            manifest=json.loads((run/'run.json').read_text()) if (run/'run.json').exists() else {}
            task_metadata=json.loads((run.parent/'task.json').read_text())
            started={r['id'] for r in jsonl_rows(run/'execution_events.jsonl') if r.get('event')=='started'}
            for job in schedule:
                mode=modes[job['mode']]
                args=mode['args']; beam=int(args[args.index('-b')+1])
                cap=max(beam,30)
                if (int(args[args.index('--max-hairpin')+1])!=cap or int(args[args.index('--max-case3-depth')+1])!=cap
                        or args[args.index('--ranking')+1]!='outer'):
                    raise ValueError('unexpected beam/cap/ranking policy')
                key=(job['input'],beam)
                if key in rows or job['warmup'] or job['repetition']!=0:
                    raise ValueError('expected one measurement per condition')
                result=byid.get(job['id'])
                if result and any(result[k]!=v for k,v in job.items()):
                    raise ValueError('checkpoint differs from frozen schedule')
                missing_status='unfinished' if job['id'] in started else 'task_failed' if task_metadata['status']=='failed' else 'not_started'
                row=dict(name=key[0],n=len(records[key[0]][1]),beam=beam,h=cap,d=cap,
                         status=result['status'] if result else missing_status, valid=False,
                         time_wall_s=result.get('time_wall_seconds') if result else None,
                         elapsed_s=result.get('elapsed_seconds') if result else None,
                         peak_rss_gib=result['peak_rss_kib']/2**20 if result and result.get('peak_rss_kib') is not None else None,
                         timeout_s=task['job_limits'].get('timeout_seconds',plan['config']['timeout_seconds']),
                         memory_limit_mib=task['job_limits'].get('memory_limit_mib',plan['config']['memory_limit_mib']),
                         plan=str(path),tier=task['tier'],task_id=task['task_id'],
                         free_bases=records[key[0]][2].count('.'),
                         _result=result,_run=run)
                for field in ('log_z','F_kcal_mol','MFE_kcal_mol','mfe_structure_sha256','bpp_pairs','max_incident_mass',
                              'F_minus_400','F_minus_400_per_nt','F_minus_previous','F_minus_previous_per_nt',
                              'tv_400_mean','tv_400_p90','tv_400_max','tv_previous_mean','tv_previous_p90','tv_previous_max'):
                    row[field]=None
                if result and result['status']!='success':
                    row['failure_detail']=result.get('error') or (result.get('stderr_summary') or {}).get('tail')
                if result and row['status']=='success':
                    command=result['command']
                    if command[command.index('--bpp-cutoff')+1]!='0' or '--bpp' not in command:
                        raise ValueError('TV requires cutoff zero BPP')
                    if (manifest.get('software_lock_valid_after_run') is not True
                            or manifest.get('software_changed_during_run') is not False
                            or manifest.get('binaries_changed_during_run') is not False):
                        row.update(status='provenance_unverified',analysis_error='task must finalize with valid dependency checks')
                    else:
                        try:
                            row.update(scalar_values(row)); row['valid']=True
                        except (ValueError,KeyError,OSError) as exc:
                            row.update(status='invalid_output',analysis_error=str(exc))
                rows[key]=row
    if len(signatures)!=1:
        raise ValueError('cannot mix software/model/builds')
    return rows,records


def write_tsv(path,rows):
    fields=list(dict.fromkeys(k for r in rows for k in r if not k.startswith('_')))
    with path.open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields,delimiter='\t');writer.writeheader()
        for row in rows:
            writer.writerow({k:'NA' if row.get(k) is None else json.dumps(row[k]) if isinstance(row[k],(list,dict)) else row[k] for k in fields})


def plot(rows, comparisons, out):
    config_dir=tempfile.TemporaryDirectory(prefix='beam-sweep-matplotlib-',dir='/tmp')
    previous_config=os.environ.get('MPLCONFIGDIR')
    os.environ['MPLCONFIGDIR']=config_dir.name
    try:
        draw_plots(rows,comparisons,out)
    finally:
        if previous_config is None:
            os.environ.pop('MPLCONFIGDIR',None)
        else:
            os.environ['MPLCONFIGDIR']=previous_config
        config_dir.cleanup()


def draw_plots(rows, comparisons, out):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    names=sorted({r['name'] for r in rows},key=lambda n:next(r['n'] for r in rows if r['name']==n))
    fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    adjacent,axadj=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    colors=plt.get_cmap('tab10')
    for index,name in enumerate(names):
        group=sorted((r for r in rows if r['name']==name),key=lambda r:r['beam'])
        x=[r['beam'] for r in group]; label=f'{name} ({group[0]["n"]:,} nt)'; color=colors(index%10)
        for ax,key in zip(axes.flat,('F_minus_400','tv_400_mean','time_wall_s','peak_rss_gib')):
            y=[r.get(key) if r['valid'] and r.get(key) is not None else np.nan for r in group]
            ax.plot(x,y,'o-',markersize=4,label=label,color=color)
        for ax,key in zip(axadj.flat,('F_minus_previous','tv_previous_mean','tv_400_p90','tv_400_max')):
            ax.plot(x,[r.get(key) if r['valid'] and r.get(key) is not None else np.nan for r in group],
                    'o-',markersize=4,label=label,color=color)
    for ax,title in zip(axes.flat,('F_b − F_400 (kcal/mol)','Mean per-base TV versus b=400','Fold wall time (s)','Peak RSS (GiB)')):
        ax.set(ylabel=title,xlabel='Beam width b');ax.grid(alpha=.2)
    for ax,title in zip(axadj.flat,('F_b − F_previous (kcal/mol)','Mean TV versus previous width','TV p90 versus b=400','Maximum TV versus b=400')):
        ax.set(ylabel=title,xlabel='Beam width b');ax.grid(alpha=.2)
    for figure,axs in ((fig,axes),(adjacent,axadj)):
        for ax in axs.flat:
            ax.axvline(100,color='.5',lw=.7,ls=':')
            ax.set_xticks(BEAMS)
        figure.legend(*axs[0,0].get_legend_handles_labels(),loc='outside lower center',ncol=2,fontsize=8)
        figure.suptitle('Beam width and local-cap sensitivity: h = d = max(b, 30)')
    for ax in (axes[0,1],axes[1,0],axes[1,1],axadj[0,1],axadj[1,0],axadj[1,1]):
        ax.set_ylim(bottom=0)
    for figure,name in ((fig,'beam_sweep'),(adjacent,'beam_sweep_adjacent')):
        for suffix in ('png','pdf'):
            figure.savefig(out/(name+'.'+suffix),dpi=180)
        plt.close(figure)


def analyze(plans,out,scratch,make_plot=True):
    if str(scratch.resolve()).startswith('/gs/bs/'):
        raise ValueError('SQLite scratch must be node-local')
    rows,records=collect(plans)
    comparisons=[]
    with tempfile.TemporaryDirectory(prefix='beam-sweep-analysis-',dir=scratch) as temporary:
        db=sqlite3.connect(str(Path(temporary)/'bpp.sqlite'))
        db.execute('PRAGMA journal_mode=OFF');db.execute('PRAGMA temp_store=FILE');db.execute('PRAGMA cache_size=-32768')
        db.execute('CREATE TABLE pairs(b INTEGER,i INTEGER,j INTEGER,p REAL,PRIMARY KEY(b,i,j)) WITHOUT ROWID')
        for name,record in records.items():
            db.execute('DELETE FROM pairs');db.commit()
            q={}
            for b in BEAMS:
                row=rows.get((name,b))
                if row and row['valid']:
                    try:
                        q[b]=import_bpp(db,row,record)
                    except (ValueError,KeyError,OSError,EOFError,sqlite3.IntegrityError) as exc:
                        db.rollback();db.execute('DELETE FROM pairs WHERE b=?',(b,));db.commit()
                        row.update(status='invalid_bpp',valid=False,analysis_error=str(exc))
            for k,b in enumerate(BEAMS):
                row=rows.get((name,b))
                if row is None:
                    continue
                for label,ref in (('400',400),('previous',BEAMS[k-1] if k else None)):
                    if ref is None:
                        continue
                    comparison=dict(name=name,beam=b,reference_beam=ref,comparison=label,
                                    status='complete' if b in q and ref in q else 'unavailable',
                                    candidate_status=row['status'],reference_status=rows.get((name,ref),{}).get('status','not_planned'),
                                    free_bases=record[2].count('.'),tv_mean=None,tv_p90=None,tv_max=None,
                                    max_pair_probability_difference=None,max_pair_coordinates=None,F_difference_kcal_mol=None)
                    if comparison['status']=='complete':
                        comparison.update(tv_compare(db,b,ref,q[b],q[ref],record[2]))
                        difference=row['F_kcal_mol']-rows[name,ref]['F_kcal_mol']
                        comparison['F_difference_kcal_mol']=difference
                        row['F_minus_'+label]=difference;row['F_minus_'+label+'_per_nt']=difference/len(record[1])
                        for stat in ('mean','p90','max'):
                            row['tv_'+label+'_'+stat]=comparison['tv_'+stat]
                    comparisons.append(comparison)
        db.close()
    out.mkdir(parents=True,exist_ok=False)
    observations=sorted(rows.values(),key=lambda r:(r['n'],r['name'],r['beam']))
    write_tsv(out/'observations.tsv',observations);write_tsv(out/'comparisons.tsv',comparisons)
    report=dict(created_utc=now(),plans=[dict(path=str(p),sha256=digest(p)) for p in plans],
                planned=len(observations),valid=sum(r['valid'] for r in observations),
                status_counts=dict(Counter(r['status'] for r in observations)),
                comparisons_complete=sum(r['status']=='complete' for r in comparisons),
                comparisons_unavailable=sum(r['status']!='complete' for r in comparisons),
                script_sha256=digest(Path(__file__)),
                interpretation='Finite-width and local-cap sensitivity; b400 is not exact; missing/timeout values not extrapolated. Single observations are not timing medians.',
                p90='linear interpolation between sorted per-free-base TV values at index .9*(N-1)',
                sign='F_b-F_400; adjacent F_b-F_previous (increasing beam); no assumed monotonicity',
                probability_tolerance=TOL,
                figures='Only valid completed folds; unavailable reference leaves gaps; timeouts stay in TSV with their cap and elapsed time')
    write_json(out/'summary.json',report)
    if make_plot:
        plot(observations,comparisons,out)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plans',nargs='+',type=Path)
    p.add_argument('--plan-report',type=Path,default=ROOT/'docs/beam_sweep_rnacentral_20261027_plan.json')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--scratch-dir',type=Path,default=Path('/tmp'))
    p.add_argument('--no-plot',action='store_true')
    a=p.parse_args()
    plans=a.plans or [ROOT/s for s in dict.fromkeys(r['plan'] for r in json.loads(a.plan_report.read_text())['arrays'])]
    print(json.dumps(analyze(plans,a.out,a.scratch_dir,not a.no_plot),indent=2))


if __name__=='__main__':
    main()
