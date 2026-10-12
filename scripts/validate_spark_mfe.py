#!/usr/bin/env python3
"""Validate packed Spark comparison evidence without extracting the archive."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import tarfile

from compare_spark_mfe import EXP, TOL, digest, pairs, parse_spark, validate, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('result', type=Path)
    args = parser.parse_args()
    result = args.result
    rows = [json.loads(x) for x in (result/'comparisons.jsonl').read_text().splitlines()]
    summary = json.loads((result/'summary.json').read_text())
    config = json.loads((result/'config.json').read_text())
    lock = json.loads(Path(config['software_lock_path']).read_text())
    assert digest(Path(config['software_lock_path'])) == config['software_lock_sha256']
    assert digest(result/'raw_outputs.tar.gz') == summary['raw_archive_sha256']
    assert len(rows) == len({r['case']['id'] for r in rows}) == 32
    for name, hash_value in lock['tracked_files_sha256'].items():
        assert digest(Path(lock['source_path'])/name) == hash_value
    assert digest(Path(lock['binary_path'])) == lock['binary_sha256']
    assert digest(Path(lock['build_path'])/'CMakeCache.txt') == lock['cmake_cache_sha256']
    assert digest(Path(lock['build_path'])/'compile_commands.json') == lock['compile_commands_sha256']
    for path, hash_value in lock['parameter_files'].items():
        assert digest(Path(path)) == hash_value
    for path, hash_value in config['linear_artifacts_before'].items():
        assert digest(Path(path)) == hash_value
    plan = json.loads((EXP/'results/local/spark-mfe-preparation-20261012/comparison_plan.json').read_text())
    planned = {r['id']:r for r in plan['cases']}
    max_local_sum_error = max_prediction_error = 0.
    nested_wi_wmb_trace_count = 0
    raw_command_count = 0
    raw_member_count = 0
    anomalies = []
    with tarfile.open(result/'raw_outputs.tar.gz', 'r:gz') as archive:
        members = {m.name:m for m in archive if m.isfile()}
        raw_member_count = len(members)
        for name in members:
            if name.endswith('.command.json') and not name.startswith('preflight/'):
                meta = json.load(archive.extractfile(members[name]))
                if 'stdout_sha256' not in meta:
                    continue
                prefix = name[:-len('.command.json')]
                assert meta['returncode'] == 0
                for stream in ['stdout','stderr']:
                    payload = archive.extractfile(members[prefix+'.'+stream]).read()
                    assert hashlib.sha256(payload).hexdigest() == meta[stream+'_sha256']
                assert meta['environment_overrides']['OMP_NUM_THREADS'] == '1'
                raw_command_count += 1
        for r in rows:
            case, name = r['case'], r['case']['id']
            assert planned[name]['sequence'] == case['sequence'] and planned[name]['scaffold'] == case['scaffold']
            def read(member):
                return archive.extractfile(members[member]).read().decode()
            official = read(name+'/official.stdout')
            assert set(parse_spark(official)['pairs']) == {tuple(p) for p in r['spark']['pairs']}
            assert official == read(name+'/observer.stdout')
            assert json.loads(json.dumps(validate(case,r['spark']))) == r['spark']['validation']
            assert r['spark']['validation']['extension_disjoint']
            assert r['spark']['validation']['extension_noncrossing'] and r['spark']['validation']['canonical']
            for policy in ['default','legacy']:
                linear = r['linear'][policy]
                expected = linear['predicted_spark_energy_of_linear_mfe']['energy']
                max_prediction_error = max(max_prediction_error,abs(expected-r['spark']['energy']))
                assert abs(expected-r['spark']['energy']) < TOL
                for target in ['linear_mfe','spark_extension']:
                    evaluated = linear[target]
                    if evaluated['evaluation']['reachable']:
                        nested_wi_wmb_trace_count += sum(x['rule']=='WI_SPLIT_WMB' for x in evaluated['loops'])
                        error = abs(sum(x['local_energy'] for x in evaluated['loops'])-evaluated['evaluation']['energy'])
                        max_local_sum_error = max(max_local_sum_error,error)
                        assert error < TOL
                assert linear['same_structure'] == (pairs(r['spark']['structure']) == pairs(linear['linear_mfe']['mfe']['completed_structure']))
                assert linear['same_energy'] == (abs(r['spark']['energy']-linear['linear_mfe']['mfe']['energy']) <= TOL)
            assert r['linear']['default']['linear_mfe']['mfe']['energy'] == r['linear']['legacy']['linear_mfe']['mfe']['energy']
            forward = r['linear_mfe_in_spark']
            if forward['exact_target_returned']:
                assert forward['adopted_evaluation_energy'] == forward['energy']
            else:
                assert forward['adopted_evaluation_energy'] is None
                assert r['reverse_trace_diagnostic']['exact_target_returned']
            nodes = [x.split('\t') for x in read(name+'/observer.stderr').splitlines() if x.startswith('loop\t')]
            if any(abs(int(x[6])) > 100000 for x in nodes):
                anomalies.append(name)
        domain = read('followup/domain.stdout')
        assert 'VP\t13\t17\t-0.19470380485332189\t0.12' in domain
        assert 'WMBP\t13\t17\t-inf\tinf' in domain and 'partial_WMBP_domain_active\t0' in domain
    validation = dict(status='passed', count=32, raw_command_count=raw_command_count,
                      raw_member_count=raw_member_count, max_loop_sum_error_kcal_mol=max_local_sum_error,
                      max_local_rounding_and_BE_coefficient_reconstruction_residual_kcal_mol=max_prediction_error,
                      official_vs_unmasked_observer_identical=32,
                      nested_wi_wmb_trace_count=nested_wi_wmb_trace_count,
                      source_parameter_binary_hashes_unchanged=True,
                      raw_trace_with_nonphysical_INF_residual_cases=anomalies,
                      note='INF residuals are traceback failures, not physical loop energies; use forward candidates and LinearCParty trace.',
                      artifact_hashes={p.name:digest(p) for p in result.iterdir() if p.is_file() and p.name!='validation.json'})
    write_json(result/'validation.json', validation)
    print(json.dumps({k:v for k,v in validation.items() if k not in ['artifact_hashes','raw_trace_with_nonphysical_INF_residual_cases']},indent=2))


if __name__ == '__main__':
    main()
