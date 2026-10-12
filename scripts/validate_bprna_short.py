#!/usr/bin/env python3
"""Check completeness and regress the short-array TSV against the saved 12-RNA grid."""
import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
import re

from run_cli import digest, now


def validate(output,pilot,tests_log,resource_log):
    p=json.loads((output/'provenance.json').read_text())
    s=json.loads((output/'summary.json').read_text())
    old=json.loads(pilot.read_text())
    pilot_manifest=pilot.parent.parent/'runs/e1-bprna-15aa7742-20261011/run.json'
    old_inputs={i['name']:i['sha256'] for i in json.loads(pilot_manifest.read_text())['inputs']}
    plan=json.loads((Path(p['main_run']['path'])/'array_plan.json').read_text())
    new_inputs={i['name']:i['sha256'] for i in plan['inputs']}
    matches={name for name,h in old_inputs.items() if new_inputs.get(name)==h}
    reference={}
    for row in old['observations']:
        if row['mode']=='exact' or row['input'] not in matches:continue
        metrics=row['metrics'];key=(row['input'],row['mode'])
        target=reference.setdefault(key,{})
        for field in ('generated_pair_l1_sum','base_state_tv_mean_Gfree','unpaired_mae_Gfree',
                      'context_tv_mean_Gfree','delta_free_energy_kcal_mol','delta_free_energy_per_nt','rt_kcal_mol'):
            if field in metrics:target[field]=metrics[field]
        for c,value in metrics.get('context_mae_Gfree',{}).items():target['context_'+c+'_mae_Gfree']=value
        for w,values in metrics.get('accessibility',{}).items():
            target.update({f'access_w{w}_{k}':v for k,v in values.items() if k!='all_windows'})
    checks=[]
    extended='capb_run' in p
    expected_jobs=886500 if extended else 664875
    expected_rows=270875 if extended else 197000
    expected_groups=295500 if extended else 221625
    def check(name,passed,**evidence):checks.append({'name':name,'passed':bool(passed),**evidence})
    text=tests_log.read_text()
    count=re.search(r'Ran (\d+) tests',text)
    check('unit_suite',count is not None and text.rstrip().endswith('OK'),log=text)
    check('cohort',p['inputs']==24625 and s['inputs']==24625,inputs=p['inputs'])
    expected={'nt50_99':12428,'nt100_199':9656,'nt200_299':1446,'nt300_399':1095}
    observed={b['bin']:b['sequences'] for b in s['bins'] if b['bin'] in expected}
    check('tier_denominators',observed==expected,observed=observed)
    check('original_failures_preserved',p['original_job_status_counts']=={'success':expected_jobs-349,'failed':349},counts=p['original_job_status_counts'])
    check('exact_replacements',len(p['exact_replacements'])==349 and p['exact_replacement_program_counts']=={'LinearPKCapR':172,'LinearPKRaccess':172,'LinearCParty':5},counts=p['exact_replacement_program_counts'])
    by_program={program:{r['input'] for r in p['exact_replacements'] if r['program']==program}
                for program in ('LinearCParty','LinearPKCapR','LinearPKRaccess')}
    check('unique_rerun_sequences',len(set.union(*by_program.values()))==172 and
          by_program['LinearCParty']<=by_program['LinearPKCapR']==by_program['LinearPKRaccess'],
          unique_sequences=len(set.union(*by_program.values())),linear_overlap=len(by_program['LinearCParty'] & by_program['LinearPKCapR']))
    check('effective_jobs_complete',p['effective_job_status_counts']=={'success':expected_jobs,'failed':0},counts=p['effective_job_status_counts'])
    check('all_metrics_complete',p['analysis_status_counts']=={'success':expected_rows},counts=p['analysis_status_counts'])
    check('all_cli_groups_consistent',p['consistency_check_counts']=={'success':expected_groups},counts=p['consistency_check_counts'])
    check('tsv_hash',digest(output/'metrics.tsv')==p['metrics_sha256'])
    hashes={name:digest(Path(__file__).parent/name) for name in p['script_sha256']}
    check('analysis_script_hashes',hashes==p['script_sha256'],hashes=hashes)
    check('pilot_source',old['software_commit']==p['main_run']['software_commit'])
    check('pilot_identical_inputs',len(matches)==12,matched=len(matches),total=len(old_inputs))
    found,comparisons,max_error,mismatches=set(),0,0.,[]
    status,rows=Counter(),0
    with (output/'metrics.tsv').open() as stream:
        for row in csv.DictReader(stream,delimiter='\t'):
            rows+=1;status[row['analysis_status']]+=1
            key=(row['input'],row['mode'])
            if key not in reference:continue
            found.add(key)
            for field,value in reference[key].items():
                actual=None if row[field]=='NA' else float(row[field]);comparisons+=1
                if value is None or actual is None:
                    passed=value is None and actual is None;error=0.
                else:
                    error=abs(actual-value);passed=error<=1e-9;max_error=max(error,max_error)
                if not passed:mismatches.append({'input':key[0],'mode':key[1],'field':field,'expected':value,'actual':actual})
    check('tsv_rows',rows==expected_rows and status=={'success':expected_rows},rows=rows,status=dict(status))
    check('pilot_metric_regression',len(found)==96 and not mismatches,
          matched_conditions=len(found),scalar_comparisons=comparisons,max_absolute_error=max_error,
          tolerance=1e-9,mismatches=mismatches)
    if extended:
        import analyze_bprna_short as short
        import numpy as np
        previous=Path(p['previous_analysis']['path'])
        check('previous_analysis_hashes',all(digest(previous/(name+'.json'))==p['previous_analysis'][name+'_sha256']
              for name in ('provenance','validation')) and digest(previous/'metrics.tsv')==p['previous_analysis']['metrics_sha256'])
        check('capb_jobs_complete',p['capb_job_status_counts']=={'success':221625},counts=p['capb_job_status_counts'])
        check('capb_groups_consistent',p['capb_consistency_check_counts']=={'success':73875})
        check('exact_groups_rechecked',p['rechecked_exact_consistency_counts']=={'success':24625})
        prior_rows=0
        with (previous/'metrics.tsv').open() as before,(output/'metrics.tsv').open() as after:
            original=iter(csv.DictReader(before,delimiter='\t'))
            preserved=True
            for row in csv.DictReader(after,delimiter='\t'):
                if row['cap']=='capb':continue
                prior=next(original,None)
                if prior is None or any(row[k]!=value for k,value in prior.items()):
                    preserved=False;break
                prior_rows+=1
            preserved=preserved and next(original,None) is None and prior_rows==197000
        check('all_previous_scalar_rows_unchanged',preserved,rows=prior_rows)
        aliases=True
        for group in s['bins']:
            conditions={(c['cap'],c['beam']):c for c in group['conditions']}
            for b in short.BEAMS:
                source='default' if b==20 else 'capb'
                adopted=conditions['adopted',b]
                measured=conditions[source,b]
                aliases=aliases and all(adopted[k]==measured[k] for k in adopted if k!='cap')
        check('adopted_matches_measured_sources',aliases)
        _,data,_,_=short.load_metrics(output)
        np.testing.assert_array_equal(data[short.CORE[0]][:,0,0],data[short.CORE[0]][:,1,0])
        check('capb_b20_alias',True)
        # Independently reconstruct paired differences and resource ratios for every tier/width.
        numeric_ok=True
        for group in s['bins']:
            chosen=(data['length'][:,0,0]>=group['min_length']) & (data['length'][:,0,0]<=group['max_length'])
            for comparison in group['paired_cap_differences']:
                bi=short.BEAMS.index(comparison['beam'])
                ci=p['caps'].index(comparison['right'])
                for key in short.METRICS:
                    diff=data[key][chosen,1,bi]-data[key][chosen,ci,bi]
                    stat=short.statistics(diff)
                    saved=comparison['metrics'][key]
                    finite=diff[np.isfinite(diff)]
                    tolerance=saved['comparison_tolerance']
                    numeric_ok=numeric_ok and all(saved[k]==value for k,value in stat.items())
                    numeric_ok=numeric_ok and saved['left_larger_count']==int(np.sum(finite>tolerance))
                    numeric_ok=numeric_ok and saved['right_larger_count']==int(np.sum(finite < -tolerance))
            for ratio in group['resource_ratios_capb_over_default']:
                bi=short.BEAMS.index(ratio['beam'])
                for key in short.RESOURCES:
                    a=data[key][chosen,1,bi];b=data[key][chosen,0,bi]
                    valid=np.isfinite(a)&np.isfinite(b)&(b>0)&(a>=0)
                    calculated=short.statistics(a[valid]/b[valid])
                    saved=ratio['metrics'][key]
                    numeric_ok=numeric_ok and all(saved[k]==value for k,value in calculated.items())
                    numeric_ok=numeric_ok and saved['invalid_or_zero_denominator']==int((~valid).sum())
        check('paired_differences_and_resource_ratios',numeric_ok)
        monotone_ok=True
        for saved in s['monotonicity']:
            if saved['cap']!='adopted':continue
            group=next(g for g in s['bins'] if g['bin']==saved['bin'])
            chosen=(data['length'][:,0,0]>=group['min_length']) & (data['length'][:,0,0]<=group['max_length'])
            matrix=data[saved['metric']][chosen,1,:]
            valid=np.all(np.isfinite(matrix),axis=1)
            changes=np.diff(matrix,axis=1)
            direction=-1 if saved['metric']=='partition_function_retention' else 1
            success=valid & np.all(direction*changes<=saved['tolerance'],axis=1)
            monotone_ok=monotone_ok and saved['eligible_sequences']==int(valid.sum()) and saved['monotone_sequences']==int(success.sum())
        check('adopted_monotonicity',monotone_ok)
    result={'schema_version':1,'generated_utc':now(),'all_passed':all(x['passed'] for x in checks),
            'unit_tests_passed':int(count[1]) if count else None,'checks':checks,
            'analysis_resources':resource_log.read_text().strip(),
            'pilot_accuracy_sha256':digest(pilot),'validation_script_sha256':digest(Path(__file__)),
            'report_script_sha256':digest(Path(__file__).with_name('report_bprna_short.py'))}
    result['plot_artifact_sha256']={name:digest(output/name) for name in
                                   ('delta_F_distribution.png','retention_distribution.png')}
    result['output_files_expected']=6
    (output/'validation.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ('all_passed','unit_tests_passed','analysis_resources')}))
    print(json.dumps({'pilot_scalar_comparisons':comparisons,'pilot_max_absolute_error':max_error,'checks_passed':sum(c['passed'] for c in checks),'checks_total':len(checks)}))
    return int(not result['all_passed'])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--pilot',type=Path,required=True)
    parser.add_argument('--tests-log',type=Path,required=True)
    parser.add_argument('--resource-log',type=Path,required=True)
    args=parser.parse_args()
    raise SystemExit(validate(args.output,args.pilot,args.tests_log,args.resource_log))
