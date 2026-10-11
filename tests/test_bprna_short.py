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
            for program in accuracy.PROGRAMS:
                command=['binary','-i',str(self.input)] + (['--exact'] if mode=='exact' else ['-b',mode[4:].split('_')[0]])
                if program=='LinearCParty':command+=['--bpp-cutoff','0','--bpp','local/bpp.tsv']
                job=dict(id=f'toy-{mode}-{program}',input='toy',mode=mode,program=program,warmup=False,repetition=0,command=command)
                planned.append(job)
                broken = failed and mode=='exact' and program=='LinearPKRaccess'
                row={**job,'status':'failed' if broken else 'success','error':'bad_alloc' if broken else None,
                     'directory':'jobs/'+job['id']}
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


if __name__=='__main__':unittest.main()
