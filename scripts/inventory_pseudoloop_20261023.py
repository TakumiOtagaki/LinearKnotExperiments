#!/usr/bin/env python3
"""Save cell-level context dependencies and a compact reviewable array config."""
import csv
import json
from pathlib import Path
import re

from rerun_pseudoloop_small_20261023 import dump,sha

ROOT=Path(__file__).resolve().parents[1]
LOCAL=ROOT/'results/local/pseudoloop-0.2.0-20261023'


def main():
    rows=[]
    def add(location,condition,metric,old,new,source):
        rows.append(dict(location=location,condition=condition,metric=metric,
            old_value=old,new_value=new,source=source))
    plot=ROOT/'results/figures/e1_fb2_final_20261023/plot_data.tsv'
    for r in csv.DictReader(plot.open(),delimiter='\t'):
        if r['metric']!='context_tv_mean_Gfree':continue
        for q in ('median','p10','p90','p99','min','max','q25','q75'):
            if r.get(q):add('Figure 2B',f"{r['cohort']}/{r['bin']}/b{r['beam']}",q,r[q],'pending array',str(plot.relative_to(ROOT)))
    # The rendered SI cells are preserved verbatim, including all length strata.
    path=ROOT/'../manuscript/latex/supplementary/results.tex'
    lines=path.read_text().splitlines();section=None;group=''
    for lineno,line in enumerate(lines,1):
        if '\\label{tab:si_accuracy_quantiles}' in line:section='S13';group=''
        if '\\label{tab:si_pk_stratified}' in line:section='S14A';group=''
        if section and '\\bottomrule' in line:section=None
        if not section or '&' not in line:continue
        fields=[s.strip() for s in line.split('&')]
        idx=6 if section=='S13' else 5
        if len(fields)<=idx or fields[idx] in ('Context TV','') or not re.fullmatch(r'[0-9.eE+\-]+',fields[idx]):continue
        if fields[0]:group=fields[0]
        add(section,group+' / '+fields[3 if section=='S13' else 2],
            'Context TV',fields[idx],'pending array',f'../manuscript/latex/supplementary/results.tex:{lineno}')
    for name in ['rna_error_quantiles.tsv','local_error_quantiles.tsv']:
        path=ROOT/'results/local/e1-pk-stratified-89baea52-20261018'/name
        for r in csv.DictReader(path.open(),delimiter='\t'):
            if not ('context' in r['metric'] or 'pseudoloop' in r['metric'] or 'pseudoloop' in r['axis']):continue
            for q in ('rna_selected','n','p50','p90','max','min'):
                if q in r:add('E1 stratified derived data',f"{r['axis']}/{r['stratum']}/b{r['beam']}",
                    r['metric']+'/'+q,r[q],'pending array',str(path.relative_to(ROOT)))
    path=ROOT/'results/local/e1-local-errors-89baea52-20261019/summary.tsv'
    for r in csv.DictReader(path.open(),delimiter='\t'):
        if r['metric']!='pseudoloop':continue
        for q in ('exact_ge_half','beam_ge_half','missed_events','created_events','missed_rnas','created_rnas','either_rnas',
                  'max_error_gt_0.1_rnas','max_error_gt_0.25_rnas','max_error_p95','max_error_p99','max_error'):
            add('S15 / E1 local derived data',f"{r['tier']}/b{r['beam']}",q,r[q],'pending array',str(path.relative_to(ROOT)))
    s=json.loads((LOCAL/'examples/summary.json').read_text())
    f=s['figure1-exact']
    add('Figure 1C / SI validation','61 nt exact','max PseudoLoop',f['old_pseudoloop_max'],f['new_pseudoloop_max'],'examples/summary.json')
    add('Figure 1C / SI validation','61 nt exact','positions P in [0.2,0.8]',f['old_pseudoloop_between_02_08'],f['new_pseudoloop_between_02_08'],'examples/summary.json')
    for mode in ['exact','b200']:
        w=s['terc-'+mode]['window_123_126']
        for c in ['total','H','B','I','M','E','P']:
            add('Figure 3 / Terc no-crossing comparison',mode+' 123-126','by-loop '+c,w['old'][c],w['new'][c],'examples/summary.json')
    old=json.loads((ROOT/'results/local/byloop-validation-89baea52-20261018/summary.json').read_text())
    new=json.loads((LOCAL/'validation/summary.json').read_text())
    for i,c in enumerate(new['contexts']):
        add('SI by-loop validation','21 exact conditions',c,old['comparisons']['exact_vs_enumeration']['component_max_abs_error'][i],
            new['comparisons']['exact_vs_enumeration']['component_max_abs_error'][i],'validation/summary.json')
    for k in new['comparisons']:
        add('SI by-loop validation',k,'maximum over components',max(old['comparisons'][k]['component_max_abs_error']),
            max(new['comparisons'][k]['component_max_abs_error']),'validation/summary.json')
    a=json.loads((LOCAL/'audit/summary.json').read_text())['bpRNA_RFAM_35940']
    add('SI example base 325','bpRNA_RFAM_35940 b200','PseudoLoop absolute error',
        abs(a['exact']['old_base325']['PseudoLoop']-a['beam200_outer']['old_base325']['PseudoLoop']),
        a['beam200_outer']['new_base325_pseudoloop_abs_error'],'audit/summary.json')
    enum=LOCAL/'enumeration/summary.json'
    if enum.exists():
        new=json.loads(enum.read_text());old=json.loads((ROOT/'results/local/e0-enumeration-89baea52-20261018/summary.json').read_text())
        for k in ['context_error','context_error_M','context_error_P']:
            add('SI enumeration validation','729 cases',k,old['maxima'][k],new['maxima'][k],'enumeration/summary.json')
    target=ROOT/'docs/pseudoloop_0.2.0_20261023_inventory.tsv'
    with target.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    # Record the exact LaTeX snapshot reviewed, including nonnumeric definitions.
    mentions=[];sources={}
    for p in sorted((ROOT/'../manuscript/latex').rglob('*.tex')):
        if any(x in p.parts for x in ['build','latexdiff']):continue
        sources[str(p.relative_to(ROOT.parent))]=sha(p)
        for i,line in enumerate(p.read_text().splitlines(),1):
            if re.search(r'PseudoLoop|Multiloop|context|loop.type',line,re.I):
                mentions.append(dict(file=str(p.relative_to(ROOT.parent)),line=i,text=line))
    with (LOCAL/'manuscript_context_mentions.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=['file','line','text'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(mentions)
    dump(LOCAL/'manuscript_review_hashes.json',sources)
    planpath=ROOT/'results/runs/pseudoloop-context-0.2.0-20261023/array_plan.json'
    p=json.loads(planpath.read_text());config=ROOT/'configs/pseudoloop_0.2.0_20261023'
    config.mkdir(exist_ok=True)
    dump(config/'e1_profiles.json',{'schema_version':1,'purpose':'profile-only context recalculation; no submission',
        'plan':str(planpath.relative_to(ROOT)),'plan_sha256':sha(planpath),
        'source_commit':p['source_commit'],'binary_sha256':p['binary_sha256'],
        'parameter_sha256':p['parameter_sha256'],'runtime_sha256':p['runtime_sha256'],
        'programs':p['programs'],'modes':p['modes'],'tiers':p['tiers'],
        'seed':p['seed'],'warmups':0,'repetitions':1,'old_plans':{k:{a:b for a,b in v.items() if a!='config'} for k,v in p['old_plans'].items()},
        'calibration':p['calibration'],'file_budget':p['file_budget'],
        'resource_policy':'cpu_16 for all; preserve original memory caps (4/8 GiB beam, 8/16 GiB exact); choose concurrency to fit 36.8 decimal GB',
        'commands':{'plan':'python3 -B scripts/pseudoloop_context_array_20261023.py plan --source NODE_LOCAL_SOURCE --build NODE_LOCAL_BUILD',
                    'check':'python3 -B scripts/pseudoloop_context_array_20261023.py check'}})
    # Job copies remain reviewable even though generated jobs/ is gitignored.
    for s in p['tiers'].values():
        script=ROOT/s['script'];(config/script.name).write_bytes(script.read_bytes())
    print(f'{len(rows)} inventoried scalar entries; {len(mentions)} LaTeX mentions; eight job scripts copied beside config')


if __name__=='__main__':main()
