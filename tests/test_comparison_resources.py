"""Guard capped timing matrices and executable dependencies behind adapters."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from plan_array import batch_wall_seconds, validate_array
from prepare_comparison_resources_20261025 import tier
from run_cli import digest, verify_lock
from analyze_comparison_resources_20261025 import classify, summarize, length_summaries


class ComparisonResourcesTests(unittest.TestCase):
    def test_failures_and_missing_repetitions_never_become_resource_medians(self):
        rows=[dict(input='rna',mode='cparty',warmup=False,status='success',
                   elapsed_seconds=1,peak_rss_kib=1024)]
        result=summarize(rows,[('rna','cparty'),('rna','spark')])
        self.assertEqual(result[0]['missing_measurements'],2)
        self.assertIsNone(result[0]['elapsed_median_seconds'])
        self.assertEqual(result[1]['missing_measurements'],3)
        self.assertEqual(classify(dict(status='failed',stderr_summary={'tail':'nonfinite oracle inside'})),'nonfinite')
        self.assertEqual(classify(dict(status='timeout')),'timeout')
        self.assertEqual(classify(dict(status='failed',returncode=139)),'signal_11')
        self.assertEqual(classify(dict(status='success',mode='cparty'),'{"log_z":Infinity,"free_energy":0,"rt":1}'),'invalid_output')

    def test_length_groups_keep_failures_and_use_paired_input_ratios(self):
        inputs={n:dict(group='bprna_250',length=250) for n in ('a','b','c')}
        records=[]
        for name,reference,candidate in [('a',2,1),('b',9,3),('c',None,100)]:
            for method,value in [('spark',reference),('fold100',candidate)]:
                records.append(dict(input=name,method=method,resource_eligible=value is not None,
                    elapsed_median_seconds=value,peak_rss_median_mib=value,
                    measurement_statuses={'success':3} if value is not None else {'timeout':3}))
        groups,pairs=length_summaries(records,inputs)
        spark=next(r for r in groups if r['method']=='spark')
        self.assertEqual((spark['planned_inputs'],spark['eligible_inputs']),(3,2))
        self.assertEqual(spark['measurement_statuses'],{'success':6,'timeout':3})
        self.assertEqual(spark['seconds'],5.5)
        self.assertEqual(pairs[0]['common_inputs'],2)
        self.assertEqual(pairs[0]['reference_over_candidate_time'],2.5)

    def test_adapter_dependency_is_checked(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for name in ('adapter','actual-tool','CMakeCache.txt','energy.par'):
                (root/name).write_text(name)
            lock=dict(schema_version=1,source_commit='fixed',cmake_cache_sha256=digest(root/'CMakeCache.txt'),
                      energy_parameter_sha256=digest(root/'energy.par'),binaries={'adapter':digest(root/'adapter')},
                      build_files={'actual-tool':digest(root/'actual-tool')})
            path=root/'lock.json';path.write_text(json.dumps(lock))
            with patch('run_cli.git_state',return_value=dict(commit='fixed',dirty=False)):
                verify_lock(path,root,root,[dict(name='adapter')],root/'energy.par')
                (root/'actual-tool').write_text('tampered')
                with self.assertRaisesRegex(ValueError,'build dependency'):
                    verify_lock(path,root,root,[dict(name='adapter')],root/'energy.par')

    def test_timeout_bound_and_required_matrix(self):
        methods=['cparty','exact','beam100','spark','fold100']
        t=tier('short',225,275,methods,10,1024,8,4)
        # Charge 20 serial warmups and three parallel measurement blocks.
        self.assertEqual(batch_wall_seconds([{'estimated_seconds':200}]*4,t),421.25)
        self.assertEqual(batch_wall_seconds([{'estimated_seconds':200}],t),200)
        config=dict(schema_version=1,purpose='development',inputs=[dict(name='a',path='unused')],
                    modes=[dict(name=m) for m in methods],programs=[dict(name='Comparison')],
                    software_lock='unused',energy_parameter_file='unused',warmups=1,repetitions=3,
                    timeout_seconds=10,memory_limit_mib=1024)
        settings=dict(schema_version=1,measurement='timing',seed=1,tiers=[t],
                      parallel_timing=dict(allowed=True,reason='bounded fixture'))
        validate_array(config,settings)
        changed=copy.deepcopy(settings)
        changed['tiers'][0]['parallel_packing']['timeout_seconds']=11
        with self.assertRaisesRegex(ValueError,'timeout bound'):
            validate_array(config,changed)
        changed=copy.deepcopy(config);changed['modes'][0]['length_range']=dict(min_length=1,max_length=500)
        with self.assertRaisesRegex(ValueError,'timeout bound'):
            validate_array(changed,settings)

    def test_bound_covers_waiting_on_the_first_pending_warmup(self):
        # Independent discrete-event replay of the existing first-pending rule.
        # Equal capped jobs make the warmup queue wait even with eight workers.
        pending=[(i,0) for i in range(20) for _ in range(2)]
        pending += [(i,block) for block in (1,2) for i in range(20)]
        active=[];clock=0
        while pending or active:
            if pending and len(active)<8:
                combo,block=pending[0]
                if not any(c==combo or b!=block for c,b,end in active):
                    pending.pop(0);active.append((combo,block,clock+10));continue
            clock=min(end for c,b,end in active)
            active=[job for job in active if job[2]>clock]
        self.assertEqual(clock,270)
        t=tier('short',225,275,['a','b','c','d','e'],10,1024,8,4)
        self.assertGreaterEqual(batch_wall_seconds([{'estimated_seconds':200}]*4,t),120+clock)


if __name__=='__main__':
    unittest.main()
