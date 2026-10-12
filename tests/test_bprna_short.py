"""Independent packed fixture: exact replacement, paired aggregation and NA."""
import copy
import csv
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
import analyze_accuracy as accuracy
import analyze_bprna_short as short
import bprna_capb as capb
from task_storage import FORMAT, RunReader, TaskStore
from run_cli import digest


class ShortAnalysisTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.record = ('toy', 'ACAAAAGU' + 'A'*42, '(......)' + '.'*42)
        self.input = self.root/'input.txt'
        self.input.write_text('>toy\n'+self.record[1]+'\n'+self.record[2]+'\n')
        self.entry = dict(name='toy',length=50,sha256=digest(self.input),copied_path=str(self.input))
        self.base = dict(status='success',purpose='fixture',config_sha256='fixture',
                         config=dict(model_settings={'temperature_celsius':37}),
                         software_before=dict(commit='fixture'),experiment_before=dict(commit='fixture'),
                         software_lock=dict(energy_parameter_sha256='params',binaries={'fixture':'hash'}),
                         software_changed_during_run=False,binaries_changed_during_run=False,
                         software_lock_valid_after_run=True,inputs=[self.entry],output_format=FORMAT)

    def pack(self, directory, modes, failed=False):
        directory.mkdir(parents=True)
        store = TaskStore(directory)
        planned = []
        for mode in modes:
            weight = .4 if mode=='exact' else {20:0,50:.2,100:.1,200:.3}[int(mode.split('_')[0][4:])]
            delta = 0 if mode=='exact' else {20:1,50:.3,100:.4,200:.1}[int(mode.split('_')[0][4:])]
            if mode.endswith('_capb'):
                beam=int(mode[4:].split('_')[0])
                weight={50:.25,100:.05,200:.35}[beam]
                delta={50:.2,100:.45,200:.05}[beam]
            for program in accuracy.PROGRAMS:
                command=['binary','-i',str(self.input)] + (['--exact'] if mode=='exact' else ['-b',mode[4:].split('_')[0]])
                if mode.endswith('_capb'):
                    beam = mode[4:].split('_')[0]
                    command += ['--max-hairpin',beam,'--max-case3-depth',beam]
                if program=='LinearCParty':command+=['--bpp-cutoff','0','--bpp','local/bpp.tsv']
                job=dict(id=f'toy-{mode}-{program}',input='toy',mode=mode,program=program,warmup=False,repetition=0,command=command)
                planned.append(job)
                broken = failed and mode=='exact' and program=='LinearPKRaccess'
                row={**job,'status':'failed' if broken else 'success','error':'bad_alloc' if broken else None,
                     'directory':'jobs/'+job['id'], 'elapsed_seconds':4 if mode.endswith('_capb') else 2,
                     'time_wall_seconds':0 if mode.startswith('beam50_') else 1,
                     'peak_rss_kib':150 if mode.endswith('_capb') else 100}
                bpp='>toy\n1 8 1\n'+(f'2 7 {weight}\n' if weight else '')
                q=[0,1-weight,1,1,1,1,1-weight,0]+[1]*42
                contexts={c:[0.]*50 for c in accuracy.CONTEXTS}
                for i in range(50):contexts['Stem'][i]=1-q[i];contexts['Exterior'][i]=q[i]
                if program=='LinearPKCapR':
                    stdout='>toy\n'+''.join(c+' '+' '.join(map(str,values))+'\n' for c,values in contexts.items())
                elif program=='LinearPKRaccess':
                    stdout='>toy\n# output_mode=probability\n# access_lens=1,4,8,16\n# scaffold='+self.record[2]+f'\n# folding_rt=.6\n# ensemble_free_energy={-1+delta}\n'
                    for start in range(1,51):
                        tokens=[]
                        for w in (1,4,8,16):
                            if start+w-1>50:continue
                            interval=set(range(start,start+w))
                            p=0 if interval & {1,8} else (1-weight if interval & {2,7} else 1)
                            tokens.append(f'{w},{p};')
                        stdout+=str(start-1)+'\t'+''.join(tokens)+'\n'
                else:stdout='>toy\n'
                raw=self.root/'raw.txt';raw.write_text(stdout)
                files={row['directory']+'/stdout.txt':raw}
                if program=='LinearCParty':
                    raw_bpp=self.root/'raw.bpp';raw_bpp.write_text(bpp);files[row['directory']+'/bpp.tsv']=raw_bpp
                store.commit(row,files)
        manifest={**self.base,'status':'failed' if failed else 'success','planned_jobs':planned,'jobs':[]}
        (directory/'run.json').write_text(json.dumps(manifest))
        return directory

    def test_packed_rerun_and_summary_keep_original_failure_and_na(self):
        run=self.root/'array';task=run/'tasks/nt50_99/000001/run'
        modes=['exact']+[f'beam{b}_{c}' for c in short.CAPS for b in short.BEAMS]
        self.pack(task,modes,failed=True)
        rerun=self.pack(self.root/'rerun',['exact'])
        # Supply only the one job which actually failed, not successful exact duplicates.
        checkpoint=rerun/'task_results.jsonl'
        rows=[json.loads(line) for line in checkpoint.read_text().splitlines()]
        rows=[r for r in rows if r['program']=='LinearPKRaccess']
        checkpoint.write_text(''.join(json.dumps(r)+'\n' for r in rows))
        m=json.loads((rerun/'run.json').read_text());m['planned_jobs']=[{k:r[k] for k in ('id','input','mode','program','warmup','repetition','command')} for r in rows]
        (rerun/'run.json').write_text(json.dumps(m))
        plan={**self.base,'tasks':[dict(tier='nt50_99',task_id=1,directory='tasks/nt50_99/000001',inputs=[self.entry])]}
        (run/'array_plan.json').write_text(json.dumps(plan))
        ordinary=accuracy.analyze(task,'exact',exact_reruns=[rerun])
        self.assertEqual(ordinary['job_status_counts'],{'success':27})
        self.assertTrue(all(r['analysis_status']=='success' for r in ordinary['observations']))
        self.assertEqual(ordinary['analysis_provenance']['original_job_status_counts'],{'success':26,'failed':1})
        self.assertEqual(len(ordinary['analysis_provenance']['exact_replacements']),1)
        out=self.root/'out'
        with patch('tarfile.TarFile.extractall',side_effect=AssertionError('must not extract')):
            provenance=short.write_metrics(run,[rerun],out)
        self.assertEqual(provenance['original_job_status_counts'],{'success':26,'failed':1})
        self.assertEqual(len(provenance['exact_replacements']),1)
        self.assertEqual(provenance['effective_job_status_counts'],{'success':27,'failed':0})
        self.assertEqual(provenance['consistency_check_counts'],{'success':9})
        with (out/'metrics.tsv').open() as stream:
            rows=list(csv.DictReader(stream,delimiter='\t'))
        self.assertEqual(len(rows),8)
        self.assertAlmostEqual(float(rows[0]['generated_pair_l1_sum']),.4)
        self.assertAlmostEqual(float(rows[0]['base_state_tv_mean_Gfree']),.8/48)
        self.assertAlmostEqual(float(rows[0]['partition_function_retention']),math.exp(-1/.6))
        self.assertEqual(rows[0]['exact_LinearPKRaccess_rerun'],'1')
        summary=short.summarize(out,make_plots=False)
        nonmonotone=next(x for x in summary['monotonicity'] if x['bin']=='all' and x['cap']=='default' and x['metric']==short.CORE[0])
        self.assertEqual(nonmonotone['monotone_sequences'],0)
        self.assertAlmostEqual(summary['nonmonotone_examples']['default_'+short.CORE[0]][0]['largest_violation'],.1)
        fine_empty=next(x for x in summary['bins'] if x['bin']=='nt75_99')
        self.assertEqual(fine_empty['conditions'][0]['metrics'][short.CORE[0]]['n'],0)
        self.assertIsNone(fine_empty['conditions'][0]['metrics'][short.CORE[0]]['median'])
        self.assertFalse((task/'jobs').exists())

    def test_exact_replacement_rejects_success_input_and_binary_changes(self):
        run=self.pack(self.root/'base',['exact'])
        rerun=self.pack(self.root/'rerun',['exact'])
        reader=RunReader(run);extra=RunReader(rerun)
        saved=copy.deepcopy(extra.manifest)
        extra.manifest['software_lock']={**saved['software_lock'],'binaries':{'fixture':'changed'}}
        with self.assertRaisesRegex(ValueError,'models or binaries'):
            accuracy.ExactReruns([extra],reader.manifest)
        extra.manifest=saved
        replacements=accuracy.ExactReruns([extra],reader.manifest)
        job=reader.manifest['jobs'][0]
        with self.assertRaisesRegex(ValueError,'successful original'):
            replacements.replace(reader,job,self.record,self.entry['sha256'])
        with self.assertRaisesRegex(ValueError,'identical S and G'):
            replacements.replace(reader,{**job,'status':'failed'},self.record,'changed')
        with self.assertRaisesRegex(ValueError,'algorithm arguments'):
            replacements.replace(reader,{**job,'status':'failed','command':job['command']+['--other']},self.record,self.entry['sha256'])

    def test_crossing_geometry_and_monotonic_tolerance(self):
        import numpy as np
        record=('x','A'*8,'..(..)..')
        bpp={'pairs':{(1,5):.4,(1,8):.2,(4,5):.3}}
        contexts={'PseudoLoop':[0.]*8}
        result=short.crossing_descriptors(record,bpp,contexts)
        self.assertEqual(result['exact_crossing_pair_count_gt_1e_minus12'],1)
        self.assertAlmostEqual(result['exact_crossing_pair_probability_sum'],.4)
        stat,_,_=short.monotonicity(np.array([[1,.5,.5+1e-10,0],[1,0,.1,0],[np.nan,1,0,0]]),1,1e-8)
        self.assertEqual((stat['eligible_sequences'],stat['monotone_sequences']),(2,1))

    def test_rna_quantiles_exclude_na_without_excluding_zero(self):
        import numpy as np
        result=short.statistics(np.array([0.,1.,10.,np.nan]))
        self.assertEqual(result['n'],3)
        self.assertEqual(result['median'],1)
        self.assertAlmostEqual(result['mean'],11/3)
        self.assertAlmostEqual(result['p90'],8.2)
        self.assertAlmostEqual(result['p99'],9.82)
        self.assertEqual((result['min'],result['max']),(0,10))

    def capb_fixture(self, failed=False):
        run=self.root/'array'
        task=run/'tasks/nt50_99/000001/run'
        modes=['exact']+short.measured_modes(short.CAPS)
        self.pack(task,modes,failed=failed)
        reruns=[]
        if failed:
            rerun=self.pack(self.root/'rerun',['exact'])
            checkpoint=rerun/'task_results.jsonl'
            rows=[json.loads(line) for line in checkpoint.read_text().splitlines()]
            rows=[r for r in rows if r['program']=='LinearPKRaccess']
            checkpoint.write_text(''.join(json.dumps(r)+'\n' for r in rows))
            m=json.loads((rerun/'run.json').read_text())
            m['planned_jobs']=[{k:r[k] for k in ('id','input','mode','program','warmup','repetition','command')} for r in rows]
            (rerun/'run.json').write_text(json.dumps(m))
            reruns=[rerun]
        task_entry=dict(tier='nt50_99',task_id=1,directory='tasks/nt50_99/000001',inputs=[self.entry])
        plan={**self.base,'inputs':[self.entry],'tasks':[task_entry]}
        (run/'array_plan.json').write_text(json.dumps(plan))
        previous=self.root/'previous'
        short.write_metrics(run,reruns,previous)
        (previous/'validation.json').write_text(json.dumps({'all_passed':True}))
        extra=self.root/'capb'
        self.pack(extra/'tasks/nt50_99/000001/run',[f'beam{b}_capb' for b in short.BEAMS[1:]])
        (extra/'array_plan.json').write_text(json.dumps(plan))
        return previous,extra

    def test_capb_extend_alias_pairing_resources_and_adopted_monotonicity(self):
        previous,extra=self.capb_fixture(failed=True)
        with patch('tarfile.TarFile.extractall',side_effect=AssertionError('must not extract')):
            p=capb.extend_capb(previous,extra,self.root/'out')
        self.assertEqual(p['original_job_status_counts'],{'success':35,'failed':1})
        self.assertEqual(p['effective_job_status_counts'],{'success':36,'failed':0})
        self.assertEqual(p['analysis_status_counts'],{'success':11})
        self.assertEqual(p['consistency_check_counts'],{'success':12})
        self.assertEqual(p['rechecked_exact_consistency_counts'],{'success':1})
        _,data,_,_=short.load_metrics(self.root/'out')
        import numpy as np
        np.testing.assert_array_equal(data[short.CORE[0]][:,0,0],data[short.CORE[0]][:,1,0])
        s=short.summarize(self.root/'out',make_plots=False)
        overall=next(g for g in s['bins'] if g['bin']=='all')
        adopted=[c for c in overall['conditions'] if c['cap']=='adopted']
        np.testing.assert_allclose([c['metrics'][short.CORE[0]]['mean'] for c in adopted],[1,.2,.45,.05],atol=1e-14)
        m=next(m for m in s['monotonicity'] if m['bin']=='all' and m['cap']=='adopted' and m['metric']==short.CORE[0])
        self.assertEqual(m['nonmonotone_sequences'],1)
        r=overall['resource_ratios_capb_over_default']
        self.assertEqual(r[0]['metrics']['LinearCParty_elapsed_seconds']['median'],1)
        self.assertEqual(r[1]['metrics']['LinearCParty_elapsed_seconds']['median'],2)
        self.assertEqual(r[1]['metrics']['LinearCParty_peak_rss_kib']['median'],1.5)
        self.assertEqual(r[1]['metrics']['LinearCParty_time_wall_seconds']['n'],0)
        self.assertEqual(r[1]['metrics']['LinearCParty_time_wall_seconds']['invalid_or_zero_denominator'],1)
        difference=overall['paired_cap_differences'][0]['metrics'][short.CORE[0]]
        self.assertEqual(difference['equal_within_tolerance_count'],1)
        difference=next(c for c in overall['paired_cap_differences'] if c['right']=='default' and c['beam']==50)['metrics'][short.CORE[0]]
        self.assertAlmostEqual(difference['mean'],-.1)
        self.assertEqual((difference['left_larger_count'],difference['right_larger_count']),(0,1))
        probability=next(c for c in overall['paired_cap_differences'] if c['right']=='default' and c['beam']==50)['metrics'][short.CORE[4]]
        self.assertAlmostEqual(probability['mean'],-.1/48)
        difference=next(c for c in overall['paired_cap_differences'] if c['right']=='default' and c['beam']==100)['metrics'][short.CORE[0]]
        self.assertAlmostEqual(difference['mean'],.05)
        self.assertEqual((difference['left_larger_count'],difference['right_larger_count']),(1,0))
        self.assertEqual(s['condition_sources']['adopted']['20'],'beam20_default')
        self.assertEqual(s['condition_sources']['adopted']['50'],'beam50_capb')
        self.assertEqual(len(list((self.root/'out').iterdir())),3)
        import report_bprna_short as report
        destination=self.root/'note.md'
        report.write_capb_note(self.root/'out',destination,s,p,
                              {'unit_tests_passed':59,'checks':[{'passed':True}],'analysis_resources':'fixture'})
        rendered=destination.read_text()
        self.assertIn('採用条件 max(b,30)',rendered)
        self.assertIn('参考値、LinearCParty',rendered)
        self.assertIn('400 nt 以上の区分に向けた含意',rendered)
        self.assertIn('e1-bprna-15aa7742-20261011-accuracy.json',rendered)

    def test_capb_rejects_changed_cohort_and_model(self):
        previous,extra=self.capb_fixture()
        plan_path=extra/'array_plan.json'
        saved=json.loads(plan_path.read_text())
        changed=copy.deepcopy(saved)
        changed['tasks'][0]['inputs'][0]['sha256']='changed'
        plan_path.write_text(json.dumps(changed))
        with self.assertRaisesRegex(ValueError,'identical complete S/G'):
            capb.extend_capb(previous,extra,self.root/'out')
        changed=copy.deepcopy(saved)
        changed['software_lock']['binaries']={'fixture':'changed'}
        plan_path.write_text(json.dumps(changed))
        with self.assertRaisesRegex(ValueError,'models or binaries'):
            capb.extend_capb(previous,extra,self.root/'out')

    def test_capb_rejects_wrong_cap_arguments(self):
        job={'command':['binary','-b','50','--max-hairpin','30','--max-case3-depth','50']}
        with self.assertRaisesRegex(ValueError,'cap rule'):
            capb.check_capb_command(job,{'command':['binary','-b','50']},50)


if __name__=='__main__':unittest.main()
