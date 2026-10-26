#!/usr/bin/env python3
"""Keep denominators, numerical failures and missing/censored resource results."""
import argparse
from collections import Counter,defaultdict
import csv
import json
import math
from pathlib import Path
import re
from statistics import median

from run_cli import ROOT,write_json
from task_storage import ArchiveReader,jsonl_rows


def classify(row,stdout=None):
    if row['status']!='success':
        error=row.get('stderr_summary',{}).get('tail','')+' '+str(row.get('error',''))
        if row['status']=='timeout':return 'timeout'
        if 'nonfinite' in error.lower():return 'nonfinite'
        if any(s in error for s in ('bad_alloc','Cannot allocate memory','cannot allocate memory')):return 'memory_limit'
        if row.get('returncode') in (-11,139):return 'signal_11'
        return row['status']
    if stdout is None:return 'success'
    try:
        if row['mode'] in ('cparty','exact','beam100'):
            values=[json.loads(s) for s in stdout.splitlines() if s.startswith('{')]
            if len(values)!=1 or not all(math.isfinite(values[0][k]) for k in ('log_z','free_energy','rt')):
                return 'invalid_output'
            if values[0]['rt']<=0:return 'invalid_output'
        else:
            found=re.findall(r'^([.()\[\]{}<>]+) \(([-\d.]+)\)$',stdout,re.M)
            if len(found)!=1 or not math.isfinite(float(found[0][1])):return 'invalid_output'
            if row['mode']=='fold100':
                energy=re.search(r'ensemble free energy: (\S+) kcal/mol',stdout)
                if not energy or not math.isfinite(float(energy[1])):return 'invalid_output'
    except (ValueError,KeyError,TypeError):return 'invalid_output'
    return 'success'


def summarize(rows,planned):
    groups=defaultdict(list)
    for row in rows:groups[row['input'],row['mode']].append(row)
    output=[]
    for name,method in sorted(set(planned)):
        rs=groups[name,method];measured=[r for r in rs if not r['warmup']]
        statuses=Counter(r.get('resource_status',classify(r)) for r in measured)
        warm=Counter(r.get('resource_status',classify(r)) for r in rs if r['warmup'])
        eligible=len(measured)==3 and statuses=={'success':3}
        rss=[r['peak_rss_kib'] for r in measured]
        output.append(dict(input=name,method=method,planned_measurements=3,completed_measurements=len(measured),
                           missing_measurements=3-len(measured),measurement_statuses=dict(statuses),warmup_statuses=dict(warm),
                           resource_eligible=eligible,
                           elapsed_median_seconds=median(r['elapsed_seconds'] for r in measured) if eligible else None,
                           peak_rss_median_mib=median(rss)/1024 if eligible and all(v is not None for v in rss) else None))
    return output


def collect(directories):
    rows=[]
    for directory in directories:
        checkpoint=directory/'task_results.jsonl'
        if not checkpoint.exists():continue
        reader=ArchiveReader(directory/'task_outputs.tar')
        for row in jsonl_rows(checkpoint):
            stdout=reader.read(row['directory']+'/stdout.txt').decode(errors='replace')
            row['resource_status']=classify(row,stdout)
            rows.append(row)
    return rows


def length_summaries(resources,inputs):
    """Summarize RNA medians; retain method and paired-comparison denominators."""
    lookup={(r['input'],r['method']):r for r in resources}
    groups=defaultdict(list)
    for name,record in inputs.items():
        group=record['group'] if record['group']!='rnacentral' else name
        groups[group].append(name)
    output=[];paired=[]
    for group,names in sorted(groups.items(),key=lambda item:min(inputs[n]['length'] for n in item[1])):
        for method in ('cparty','exact','beam100','spark','fold100'):
            rs=[lookup[n,method] for n in names if (n,method) in lookup]
            if not rs:continue
            good=[r for r in rs if r['resource_eligible']]
            statuses=Counter()
            for r in rs:statuses.update(r['measurement_statuses'])
            output.append(dict(group=group,method=method,planned_inputs=len(rs),eligible_inputs=len(good),
                length_min=min(inputs[n]['length'] for n in names),length_max=max(inputs[n]['length'] for n in names),
                length_median=median(inputs[n]['length'] for n in names),measurement_statuses=dict(statuses),
                seconds=median(r['elapsed_median_seconds'] for r in good) if good else None,
                mib=median(r['peak_rss_median_mib'] for r in good) if good else None))
        for reference,candidate in [('cparty','beam100'),('cparty','exact'),('spark','fold100')]:
            common=[n for n in names if all((n,m) in lookup and lookup[n,m]['resource_eligible'] for m in (reference,candidate))]
            if not common:continue
            paired.append(dict(group=group,reference=reference,candidate=candidate,common_inputs=len(common),inputs=common,
                reference_seconds=median(lookup[n,reference]['elapsed_median_seconds'] for n in common),
                candidate_seconds=median(lookup[n,candidate]['elapsed_median_seconds'] for n in common),
                reference_mib=median(lookup[n,reference]['peak_rss_median_mib'] for n in common),
                candidate_mib=median(lookup[n,candidate]['peak_rss_median_mib'] for n in common),
                reference_over_candidate_time=median(lookup[n,reference]['elapsed_median_seconds']/lookup[n,candidate]['elapsed_median_seconds'] for n in common),
                reference_over_candidate_rss=median(lookup[n,reference]['peak_rss_median_mib']/lookup[n,candidate]['peak_rss_median_mib'] for n in common),
                candidate_faster_inputs=sum(lookup[n,candidate]['elapsed_median_seconds']<lookup[n,reference]['elapsed_median_seconds'] for n in common)))
    return output,paired


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan-report',type=Path,default=ROOT/'docs/comparison_resources_20261025_plan.json')
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    report=json.loads(args.plan_report.read_text());directories=[];planned=[]
    for name in sorted({r['run'] for r in report['arrays']}):
        plan=json.loads((ROOT/'results/runs'/name/'array_plan.json').read_text())
        for task in plan['tasks']:
            directory=Path(plan['output_dir'])/task['directory']
            metadata=json.loads((directory/'task.json').read_text())
            if metadata['status']=='running':raise ValueError('wait for all tasks to stop before aggregation')
            directories.append(directory/'run')
            planned.extend((r['name'],m['name']) for r in task['inputs'] for m in plan['config']['modes'])
    rows=collect(directories)
    result=summarize(rows,planned)
    inputs={r['name']:r for r in json.loads((ROOT/'data/comparison-resources-20261025/manifest.json').read_text())['records']}
    for r in result:r.update(length=inputs[r['input']]['length'],group=inputs[r['input']]['group'])
    grouped,paired=length_summaries(result,inputs)
    # PF equality is checked only for CParty/exact, never for Spark/fold.
    parity=[]
    values=defaultdict(dict)
    for directory in directories:
        checkpoint=directory/'task_results.jsonl'
        if not checkpoint.exists():continue
        reader=ArchiveReader(directory/'task_outputs.tar')
        for row in jsonl_rows(checkpoint):
            if row['status']=='success' and not row['warmup'] and row['mode'] in ('cparty','exact'):
                text=reader.read(row['directory']+'/stdout.txt').decode()
                if classify(row,text)=='success':
                    values[row['input'],row['repetition']][row['mode']]=next(json.loads(s) for s in text.splitlines() if s.startswith('{'))
    for (name,rep),v in values.items():
        if set(v)=={'cparty','exact'}:
            dz=abs(v['cparty']['log_z']-v['exact']['log_z']);df=abs(v['cparty']['free_energy']-v['exact']['free_energy'])
            parity.append(dict(input=name,repetition=rep,log_z_absolute_error=dz,free_energy_absolute_error=df,
                               passed=dz<=1e-9 and df<=1e-10))
    args.out.mkdir(parents=True,exist_ok=False)
    write_json(args.out/'summary.json',dict(resources=result,exact_exclusions=report['exact_exclusions'],
                                           pf_parity=parity,planned_measurements=len(planned)*3,
                                           completed_measurements=sum(not r['warmup'] for r in rows),
                                           length_groups=grouped,paired_comparisons=paired,
                                           extrapolation=False,mfe_energy_or_structure_parity_claim=False))
    with (args.out/'resources.tsv').open('w') as output:
        writer=csv.DictWriter(output,fieldnames=list(result[0]),delimiter='\t');writer.writeheader()
        for row in result:
            writer.writerow({k:json.dumps(v,sort_keys=True) if isinstance(v,dict) else v for k,v in row.items()})
    for filename,records in [('length_groups.tsv',grouped),('paired_comparisons.tsv',paired)]:
        with (args.out/filename).open('w') as output:
            writer=csv.DictWriter(output,fieldnames=list(records[0]),delimiter='\t');writer.writeheader()
            for row in records:writer.writerow({k:json.dumps(v,sort_keys=True) if isinstance(v,(dict,list)) else v for k,v in row.items()})
    print(json.dumps(dict(planned=len(planned)*3,completed=sum(not r['warmup'] for r in rows),eligible=sum(r['resource_eligible'] for r in result))))


if __name__=='__main__':
    main()
