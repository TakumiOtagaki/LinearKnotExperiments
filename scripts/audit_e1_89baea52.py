"""Audit E1 summaries and export manuscript values; no RNA calculations."""
from pathlib import Path
import csv, json, sys
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'scripts'))
from report_e1_production import validate_row
from run_cli import digest
import numpy as np

suffix={'short':'20261018','long':'cpu16-20261018'}
report={'revision':'89baea52','runs':{},'cohorts':{},'comparisons':{}}
for cohort in suffix:
 run=root/f'results/runs/e1-bprna-{cohort}-89baea52-{suffix[cohort]}'
 j=json.loads((run/'run.json').read_text())
 assert j['status']=='success' and j['missing_measured_jobs']==0
 assert j['software_before']['commit'].startswith('89baea52')
 assert j['software_changed_during_run'] is False and j['binaries_changed_during_run'] is False
 assert j['software_lock_valid_after_run'] is True
 report['runs'][cohort]={k:j[k] for k in ['task_status_counts','planned_measured_jobs','completed_measured_jobs','missing_measured_jobs']}
 directory=root/f'results/local/e1-bprna-{cohort}-89baea52-{suffix[cohort]}'
 provenance=json.loads((directory/'provenance.json').read_text())
 assert digest(directory/'metrics.tsv')==provenance['metrics_sha256']
 assert all(r['software_commit'].startswith('89baea52') for r in provenance['runs'])
 summary=json.loads((directory/'summary.json').read_text())
 old_dir=root/f'results/local/e1-bprna-{cohort}-d21618c5-20261016'
 old=json.loads((old_dir/'summary.json').read_text())
 for key in ['quantile_method','unit','delta_F_definition','opening_domain','missing_definition','reference_mode','comparison_kind']:
  assert summary[key]==old[key],key
 assert [(g['bin'],g['beam'],g['rna']) for g in summary['groups']]==[(g['bin'],g['beam'],g['rna']) for g in old['groups']]
 ids=set(); selected={}; old_ids=set(); id_hashes={}
 for row in csv.DictReader((directory/'metrics.tsv').open(),delimiter='\t'):
  validate_row(row)
  key=(row['input'],row['beam']);assert key not in ids;ids.add(key);id_hashes[key]=row['input_sha256']
  for field in ['delta_free_energy_kcal_mol','base_state_tv_mean_Gfree','context_tv_mean_Gfree']:
   assert row[field]!='NA',field
  if cohort=='short' and row['beam']=='200':
   selected[row['input']]=row
 for row in csv.DictReader((old_dir/'metrics.tsv').open(),delimiter='\t'):
  key=(row['input'],row['beam']);old_ids.add(key)
  assert id_hashes[key]==row['input_sha256'], 'S/G input changed from the previous cohort'
 assert ids==old_ids
 assert len(ids)==provenance['metric_rows']
 report['cohorts'][cohort]={'rows_validated':len(ids),'same_cohort_and_definitions_as_previous':True,
  'new_metrics_sha256':provenance['metrics_sha256'],'new_groups':summary['groups'], 'old_groups':old['groups']}
 if cohort=='short':
  below200=[float(r['abs_delta_free_energy_kcal_mol']) for r in selected.values() if int(r['length'])<200]
  mid=[float(r['abs_delta_free_energy_kcal_mol']) for r in selected.values() if 200<=int(r['length'])<400]
  report['short_abstract']={'below200_n':len(below200),'below200_max_absolute':max(below200),
    'nt200_399_n':len(mid),'nt200_399_p90_absolute':float(np.quantile(mid,.9,method='linear'))}
strata=root/'results/local/e1-pk-stratified-89baea52-20261018'
p=json.loads((strata/'provenance.json').read_text());assert p['validation']=='passed'
assert all(v<=1e-12 for v in p['previous_metric_max_absolute_differences'].values())
for name,h in p['output_sha256'].items():assert digest(strata/name)==h
report['stratification']={'provenance':p,'counts':list(csv.DictReader((strata/'stratum_counts.tsv').open(),delimiter='\t')),
 'rna':list(csv.DictReader((strata/'rna_error_quantiles.tsv').open(),delimiter='\t')),
 'local':list(csv.DictReader((strata/'local_error_quantiles.tsv').open(),delimiter='\t'))}
specs=[
 ('below200_max_abs_F','short',None,200,'abs_delta_free_energy_kcal_mol','max','kcal/mol'),
 ('short300_p90_F_b20','short','nt300_399',20,'delta_free_energy_kcal_mol','p90','kcal/mol'),
 ('short300_p90_F_b200','short','nt300_399',200,'delta_free_energy_kcal_mol','p90','kcal/mol'),
 ('short300_max_F_b200','short','nt300_399',200,'delta_free_energy_kcal_mol','max','kcal/mol'),
 ('short300_p90_BPP_b200','short','nt300_399',200,'base_state_tv_mean_Gfree','p90','TV'),
 ('short300_p90_context_b200','short','nt300_399',200,'context_tv_mean_Gfree','p90','TV'),
 ('short300_p90_w16_b200','short','nt300_399',200,'access_w16_probability_mae_eligible','p90','probability'),
 ('long_p90_F_b100','long','all',100,'delta_free_energy_kcal_mol','p90','kcal/mol'),
 ('long_p90_F_b300','long','all',300,'delta_free_energy_kcal_mol','p90','kcal/mol'),
 ('long_p90_BPP_b300','long','all',300,'base_state_tv_mean_Gfree','p90','TV'),
 ('long_p90_context_b300','long','all',300,'context_tv_mean_Gfree','p90','TV'),
 ('long_p90_w16_b300','long','all',300,'access_w16_probability_mae_eligible','p90','probability')]
for name,c,t,b,m,q,unit in specs:
 vals={}
 for kind in ('old','new'):
  gs=report['cohorts'][c][kind+'_groups']
  if t is None:
   vals[kind]=max(g['metrics'][m][q] for g in gs if g['beam']==b and g['bin'] in ('nt50_99','nt100_199'))
  else:
   vals[kind]=next(g['metrics'][m][q] for g in gs if g['beam']==b and g['bin']==t)
 report['comparisons'][name]=dict(**vals,unit=unit)
report['source_hashes']={str(path.relative_to(root)):digest(path) for path in [
 *[root/f'results/local/e1-bprna-{c}-89baea52-{s}/summary.json' for c,s in suffix.items()],strata/'provenance.json']}
out=root/'results/local/e1-paper-89baea52-20261019';out.mkdir(exist_ok=False)
report['all_passed']=True
report['audit_script_sha256']=digest(Path(__file__))
(out/'paper_values.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'output':str(out),'abstract':report['short_abstract'],'all_passed':True}))
