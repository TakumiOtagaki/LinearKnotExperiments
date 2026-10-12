#!/usr/bin/env python3
"""Small offline Spark/LinearCParty comparison; one process and one thread at a time.

Build the unmodified Spark checkout out of tree; compile observers in /tmp.
Keep stdout/stderr/commands in one tar and comparisons in one JSONL. Masked
Spark runs exclude extra pairs but do not force targets; adopt their energies
only after checking exact pair-set equality. Trace-repair runs are diagnostic.
"""
import argparse
import collections
import hashlib
import io
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tarfile
import tempfile
import time

from spark_mfe_observer import generate

EXP = Path(__file__).resolve().parents[1]
COMMIT = '3105c311ca6243a31788e74ac33640435e53b62c'
TOL = 1e-9


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def provenance(source):
    def git(*args):
        return subprocess.check_output(['git', '-C', str(source), *args], text=True).strip()
    files = subprocess.check_output(['git', '-C', str(source), 'ls-files', '-z']).split(b'\0')
    hashes = {os.fsdecode(f): digest(source / os.fsdecode(f)) for f in files if f and (source/os.fsdecode(f)).is_file()}
    gitlinks = {line.split('\t', 1)[1]: line.split()[1] for line in git('ls-files', '-s').splitlines()
                if line.startswith('160000 ')}
    return dict(commit=git('rev-parse', 'HEAD'), tree=git('rev-parse', 'HEAD^{tree}'),
                status=git('status', '--porcelain'), gitlinks=gitlinks, tracked_files_sha256=hashes,
                tracked_tree_sha256=hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest())


def pairs(structure):
    stacks = {c: [] for c in '([{<'}
    out = set()
    for i, c in enumerate(structure, 1):
        if c in stacks:
            stacks[c].append(i)
        elif c in ')]}>':
            left = '([{<'[')]}>'.index(c)]
            if not stacks[left]:
                raise ValueError('unbalanced structure')
            out.add((stacks[left].pop(), i))
        elif c != '.':
            raise ValueError('unexpected structure character')
    if any(stacks.values()):
        raise ValueError('unbalanced structure')
    return out


def render(extension, n):
    out = ['.'] * n
    for i, j in extension:
        out[i-1], out[j-1] = '(', ')'
    return ''.join(out)


def parse_spark(stdout):
    matches = re.findall(r'^([.()\[\]{}<>]+) \(([-\d.]+)\)$', stdout, re.M)
    if len(matches) != 1:
        raise ValueError('Spark did not emit one parseable result')
    structure, energy = matches[0]
    return dict(structure=structure, energy=float(energy), pairs=sorted(pairs(structure)),
                stdout_error_lines=[x for x in stdout.splitlines() if x.startswith('Error:')])


def validate(case, result):
    g, r = pairs(case['scaffold']), pairs(result['structure'])
    extension = r - g
    used_g = {x for p in g for x in p}
    used = [x for p in extension for x in p]
    canonical = {'AU', 'UA', 'CG', 'GC', 'GU', 'UG'}
    return dict(missing_scaffold_pairs=sorted(g-r),
                retains_scaffold=g <= r,
                extension_disjoint=not (set(used) & used_g) and len(used) == len(set(used)),
                extension_noncrossing=not any(a < c < b < d or c < a < d < b
                                               for a, b in extension for c, d in extension),
                canonical=all(case['sequence'][a-1]+case['sequence'][b-1] in canonical for a, b in r),
                minimum_pair_distance=min((b-a for a, b in r), default=None))


def parse_linear(stdout):
    result, loops = {}, []
    for line in stdout.splitlines():
        c = line.split('\t')
        if c[0] in ('mfe', 'evaluation'):
            energy = float(c[3])
            result[c[0]] = dict(ok=c[1] == '1', reachable=c[2] == '1',
                               energy=energy if math.isfinite(energy) else None,
                               extension=c[4], completed_structure=c[5])
        elif c[0] == 'loop':
            loops.append(dict(depth=int(c[1]), state=c[2], i=int(c[3]), j=int(c[4]),
                              ip=int(c[5]), jp=int(c[6]), rule=c[7],
                              local_energy=float(c[8]), subtree_energy=float(c[9])))
        elif c[0] == 'rt':
            result['rt_kcal_mol'] = float(c[1])
        elif c[0] == 'parameters':
            result['parameters'] = dict(temperature_celsius=float(c[1]), dangles=int(c[2]), pf_smooth=int(c[3]),
                                        MLclosing=float(c[4]), MLintern_CG=float(c[5]), MLbase=float(c[6]))
    result['loops'] = loops
    if result['evaluation']['reachable']:
        result['local_sum_error'] = sum(x['local_energy'] for x in loops)-result['evaluation']['energy']
        assert abs(result['local_sum_error']) <= TOL
    return result


def predicted_integer_energy(linear):
    """Reconstruct observed Spark alternatives; not an independent Spark evaluator."""
    if not linear['evaluation']['reachable']:
        return None
    changes = []
    for loop in linear['loops']:
        rule, energy = loop['rule'], loop['local_energy']
        if rule in ('BE_STACK', 'VP_STACK', 'BE_INTERNAL_LOOP', 'VP_INTERNAL_LOOP'):
            candidate = energy
            # Spark compute_BE also permits the internal alternative on a stack.
            if rule == 'BE_STACK' and energy > 0:
                candidate = energy / .89 * .74
            integer_energy = round(candidate * 100) / 100
            changes.append(dict(rule=rule, i=loop['i'], j=loop['j'],
                                linear_local_energy=energy, spark_candidate_before_rounding=candidate,
                                spark_integer_local_energy=integer_energy,
                                coefficient_difference=candidate-energy,
                                rounding_difference=integer_energy-candidate))
    return dict(energy=linear['evaluation']['energy'] + sum(
        c['spark_integer_local_energy']-c['linear_local_energy'] for c in changes), changes=changes)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=EXP/'deps/Spark-upstream')
    parser.add_argument('--build', type=Path, default=EXP/'build/spark-3105c311')
    parser.add_argument('--linear-source', type=Path, default=EXP/'deps/LinearCParty-15aa7742')
    parser.add_argument('--linear-build', type=Path, default=EXP/'build/linear-cparty-15aa7742')
    parser.add_argument('--preparation', type=Path, default=EXP/'results/local/spark-mfe-preparation-20261012')
    parser.add_argument('--out', type=Path, default=EXP/'results/local/spark-mfe-comparison-20261012')
    parser.add_argument('--lock', type=Path, default=EXP/'manifests/software-spark-3105c311-20261012.json')
    parser.add_argument('--preflight', type=Path, help='Archive prior local smoke/build logs')
    args = parser.parse_args()
    source, build = args.source.resolve(), args.build.resolve()
    ls, lb = args.linear_source.resolve(), args.linear_build.resolve()
    parameter = ls/'params/rna_Turner04.par'
    args.out.mkdir(parents=True, exist_ok=False)
    before = provenance(source)
    linear_before = provenance(ls)
    assert before['commit'] == COMMIT and not before['status']
    assert linear_before['commit'].startswith('15aa7742') and not linear_before['status']
    frozen_linear = {str(p): digest(p) for p in [parameter, lb/'libcparty_scfg2_exact.a', lb/'libRNA.a']}
    env = {k: v for k, v in os.environ.items() if not k.startswith(('CPARTY_', 'SPARK_OBSERVER_'))}
    overrides = dict(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
                     NUMEXPR_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1', LC_ALL='C',
                     GIT_OPTIONAL_LOCKS='0', CPARTY_FIXED_ENERGY_PARAM_FILE=str(parameter),
                     CPARTY_PF_OVERRIDE_SCALE='1')
    env.update(overrides)
    rows = [json.loads(x) for x in (args.preparation/'linear_mfe.jsonl').read_text().splitlines()]
    assert len(rows) == 32
    comparisons = []
    with tempfile.TemporaryDirectory(prefix='spark-mfe-', dir='/tmp') as tmp, \
            tarfile.open(args.out/'raw_outputs.tar.gz', 'w:gz') as archive:
        scratch = Path(tmp)

        def add(name, content):
            if isinstance(content, str):
                content = content.encode()
            info = tarfile.TarInfo(name); info.size = len(content)
            archive.addfile(info, io.BytesIO(content))

        def run(name, command, extra=None, timeout=120):
            effective = dict(env, **(extra or {}))
            start = time.time()
            try:
                p = subprocess.run([str(x) for x in command], capture_output=True, env=effective, timeout=timeout)
                stdout, stderr, returncode, status = p.stdout, p.stderr, p.returncode, 'completed'
            except subprocess.TimeoutExpired as error:
                stdout, stderr, returncode, status = error.stdout or b'', error.stderr or b'', None, 'timeout'
            metadata = dict(command=[str(x) for x in command], environment_overrides=dict(overrides, **(extra or {})),
                            start_unix=start, end_unix=time.time(), returncode=returncode, status=status,
                            stdout_sha256=hashlib.sha256(stdout).hexdigest(), stderr_sha256=hashlib.sha256(stderr).hexdigest())
            for suffix, payload in [('stdout', stdout), ('stderr', stderr), ('command.json', json.dumps(metadata, indent=2))]:
                add(name+'.'+suffix, payload)
            if returncode != 0:
                raise RuntimeError(f'{name}: {status}, returncode {returncode}; raw output retained')
            return stdout.decode(), stderr.decode(), metadata

        if args.preflight:
            for p in sorted(args.preflight.rglob('*')):
                if p.is_file():
                    archive.add(p, arcname='preflight/'+str(p.relative_to(args.preflight)), recursive=False)
        configure = ['cmake', '-S', source, '-B', build, '-DCMAKE_BUILD_TYPE=Release',
                     '-DCMAKE_CXX_COMPILER=g++', '-DCMAKE_EXPORT_COMPILE_COMMANDS=ON']
        build_command = ['cmake', '--build', build, '--parallel', '1']
        run('build/configure', configure)
        run('build/build', build_command)
        generator_files = [EXP/'scripts/spark_mfe_observer.py', EXP/'scripts/spark_mfe_linear_probe.cc', Path(__file__)]
        for p in generator_files:
            archive.add(p, arcname='scripts/'+p.name, recursive=False)
        for variant, repaired in [('spark_observer', False), ('spark_trace_diagnostic', True)]:
            file = scratch/(variant+'.cc'); generate(source, file, repaired)
            command = ['g++', '-std=c++17', '-O1', '-DNDEBUG', '-I'+str(source/'src'),
                       '-DPARAMS_DIR="'+str(source/'params')+'"', file,
                       build/'CMakeFiles/Spark.dir/src/cmdline.cc.o', build/'CMakeFiles/Spark.dir/src/trace_arrow.cc.o',
                       build/'libRNA.a', '-o', scratch/variant]
            run('build/'+variant, command)
            archive.add(file, arcname='observer/'+file.name, recursive=False)
            archive.add(scratch/variant, arcname='observer/'+variant, recursive=False)
        command = ['g++', '-std=c++20', '-O1', '-I'+str(ls/'src'), EXP/'scripts/spark_mfe_linear_probe.cc',
                   lb/'libcparty_scfg2_exact.a', lb/'libRNA.a', '-o', scratch/'linear_probe']
        run('build/linear_probe', command)
        archive.add(scratch/'linear_probe', arcname='observer/linear_probe', recursive=False)
        domain_source = EXP/'scripts/spark_mfe_domain_probe.cc'
        command = ['g++', '-std=c++20', '-O1', '-I'+str(ls/'src'), domain_source,
                   lb/'libcparty_scfg2_exact.a', lb/'libRNA.a', '-o', scratch/'domain_probe']
        run('build/domain_probe', command)
        run('followup/domain', [scratch/'domain_probe'])
        archive.add(domain_source, arcname='scripts/'+domain_source.name, recursive=False)
        archive.add(scratch/'domain_probe', arcname='observer/domain_probe', recursive=False)
        for p in [parameter, source/'params/rna_DirksPierce09.par', source/'src/PK_globals.hh']:
            archive.add(p, arcname='parameters/'+p.name, recursive=False)
        for p in sorted(build.rglob('*')):
            if p.is_file():
                archive.add(p, arcname='official_build/'+str(p.relative_to(build)), recursive=False)
        lock = dict(status='built_verified', usable_software_lock=True, purpose='small development comparison',
                    repository='git@github.com:TheCOBRALab/Spark.git', source_path=str(source),
                    source_status=before['status'], **{k: v for k, v in before.items() if k != 'status'},
                    build_revision_verified=True, build_path=str(build), configure_command=list(map(str, configure)),
                    build_command=list(map(str, build_command)), build_type='Release', generator='Unix Makefiles',
                    compiler=subprocess.check_output(['g++', '--version'], text=True).splitlines()[0],
                    cmake_version=subprocess.check_output(['cmake', '--version'], text=True).splitlines()[0],
                    flags=(build/'CMakeFiles/Spark.dir/flags.make').read_text(),
                    linker_command=(build/'CMakeFiles/Spark.dir/link.txt').read_text(),
                    binary_path=str(build/'Spark'), binary_sha256=digest(build/'Spark'),
                    cmake_cache_sha256=digest(build/'CMakeCache.txt'),
                    compile_commands_sha256=digest(build/'compile_commands.json'),
                    parameter_files={str(p): digest(p) for p in [parameter, source/'params/rna_Turner04.par',
                                                               source/'params/rna_DirksPierce09.par']},
                    pk_globals_sha256=digest(source/'src/PK_globals.hh'),
                    temperature_celsius=37, dangles=2, gu_allowed=True, pseudoknot_prediction=True, pk_only=False,
                    environment_overrides=overrides,
                    observer_binary_sha256={n: digest(scratch/n) for n in ['linear_probe','spark_observer','spark_trace_diagnostic']},
                    observer_note='Separate generated diagnostics; official binary and checkout are unmodified.',
                    retained_build='Spark, CMakeCache.txt, compile_commands.json; other build files packed in raw tar',
                    observed_on='2026-10-12')
        write_json(args.lock, lock)
        config = dict(purpose='small offline development comparison', count=32, source_before=before,
                      linear_source_before=linear_before, linear_artifacts_before=frozen_linear,
                      preparation_hashes={p.name: digest(p) for p in [args.preparation/'comparison_plan.json',
                                                                    args.preparation/'linear_mfe.jsonl']},
                      software_lock_path=str(args.lock), software_lock_sha256=digest(args.lock), tolerance_kcal_mol=TOL,
                      file_budget_estimate=200, network_used=False, qsub_used=False, git_commit_used=False)
        write_json(args.out/'config.json', config)

        def spark_command(binary, case):
            return [binary, '-d', '2', '-P', parameter, '-r', case['scaffold'], case['sequence']]

        for old in rows:
            case = dict(id=old['case'], sequence=old['sequence'], scaffold=old['scaffold'])
            name = case['id']
            out, _, meta = run(name+'/official', spark_command(build/'Spark', case))
            spark = parse_spark(out); spark['validation'] = validate(case, spark)
            obs_out, obs_err, _ = run(name+'/observer', spark_command(scratch/'spark_observer', case))
            assert obs_out == out, 'unmasked observer changed official output'
            metadata_line = next(x for x in obs_err.splitlines() if x.startswith('metadata\t'))
            assert metadata_line == 'metadata\t37\t2\t930\t-90\t0'
            item = dict(case=case, spark=spark, official_run=meta, effective_metadata=metadata_line,
                        spark_trace_raw_member=name+'/observer.stderr', linear={})
            g = pairs(case['scaffold']); gp = pairs(spark['structure'])-g
            for policy in ['default', 'legacy']:
                linear = {}
                for target, structure in [('spark_extension', render(gp, len(case['sequence']))),
                                          ('linear_mfe', old['linear_'+policy]['extension'])]:
                    stdout, _, _ = run(name+'/'+policy+'/'+target,
                                       [scratch/'linear_probe', case['sequence'], case['scaffold'], structure, policy])
                    linear[target] = parse_linear(stdout)
                mfe = linear['linear_mfe']['mfe']
                assert abs(mfe['energy']-old['linear_'+policy]['energy']) < TOL
                assert pairs(mfe['completed_structure']) == pairs(old['linear_'+policy]['completed_structure'])
                linear['same_energy'] = abs(mfe['energy']-spark['energy']) <= TOL
                linear['same_structure'] = pairs(mfe['completed_structure']) == pairs(spark['structure'])
                linear['spark_evaluation_restores_missing_scaffold'] = not spark['validation']['retains_scaffold']
                linear['predicted_spark_energy_of_linear_mfe'] = predicted_integer_energy(linear['linear_mfe'])
                linear['predicted_spark_energy_of_spark_extension_with_G'] = predicted_integer_energy(linear['spark_extension'])
                item['linear'][policy] = linear
            target_mask = render(pairs(item['linear']['default']['linear_mfe']['mfe']['completed_structure'])-g,
                                 len(case['sequence']))
            stdout, _, _ = run(name+'/masked_observer', spark_command(scratch/'spark_observer', case),
                               dict(SPARK_OBSERVER_TARGET=target_mask))
            reverse = parse_spark(stdout)
            reverse['exact_target_returned'] = pairs(reverse['structure']) == pairs(
                item['linear']['default']['linear_mfe']['mfe']['completed_structure'])
            reverse['adopted_evaluation_energy'] = reverse['energy'] if reverse['exact_target_returned'] else None
            item['linear_mfe_in_spark'] = reverse
            stdout, _, _ = run(name+'/trace_diagnostic', spark_command(scratch/'spark_trace_diagnostic', case))
            repaired = parse_spark(stdout); repaired['validation'] = validate(case, repaired)
            assert abs(repaired['energy']-spark['energy']) < TOL
            item['trace_diagnostic'] = repaired
            if not reverse['exact_target_returned']:
                stdout, _, _ = run(name+'/masked_trace_diagnostic', spark_command(scratch/'spark_trace_diagnostic', case),
                                   dict(SPARK_OBSERVER_TARGET=target_mask))
                diagnostic = parse_spark(stdout)
                diagnostic['exact_target_returned'] = pairs(diagnostic['structure']) == pairs(
                    item['linear']['default']['linear_mfe']['mfe']['completed_structure'])
                item['reverse_trace_diagnostic'] = diagnostic
            comparisons.append(item)
            with (args.out/'comparisons.jsonl').open('a') as stream:
                stream.write(json.dumps(item, allow_nan=False)+'\n')

        controls = []
        for case in [dict(id='terminal6', sequence='CAAAGU', scaffold='(...).'),
                     dict(id='gu_stack9', sequence='GUGUGUGUG', scaffold='.((....))')]:
            out, _, _ = run('controls/'+case['id']+'/official', spark_command(build/'Spark', case))
            obs, err, _ = run('controls/'+case['id']+'/observer', spark_command(scratch/'spark_observer', case))
            assert out == obs
            diag, _, _ = run('controls/'+case['id']+'/trace_diagnostic', spark_command(scratch/'spark_trace_diagnostic', case))
            target = '([..)]' if case['id'] == 'terminal6' else '[((..].))'
            lin, _, _ = run('controls/'+case['id']+'/linear', [scratch/'linear_probe',case['sequence'],case['scaffold'],target,'default'])
            controls.append(dict(case=case, spark=parse_spark(out), diagnostic=parse_spark(diag),
                                 linear=parse_linear(lin), candidate_evidence=[x for x in err.splitlines()
                                                                         if x.startswith('BE_candidates')]))
        summary = dict(count=32, controls=controls, categories={}, source_after=provenance(source),
                       linear_source_after=provenance(ls), linear_artifacts_after={p: digest(Path(p)) for p in frozen_linear})
        for policy in ['default', 'legacy']:
            counts = collections.Counter()
            for r in comparisons:
                comparison = r['linear'][policy]
                counts[('same_structure' if comparison['same_structure'] else 'different_structure')+'_'+
                       ('same_energy' if comparison['same_energy'] else 'different_energy')] += 1
            summary['categories'][policy] = dict(counts)
        summary['retains_scaffold'] = sum(r['spark']['validation']['retains_scaffold'] for r in comparisons)
        summary['spark_extension_reachable_in_linear'] = sum(
            r['linear']['default']['spark_extension']['evaluation']['reachable'] for r in comparisons)
        summary['reverse_exact_target_returned'] = sum(r['linear_mfe_in_spark']['exact_target_returned'] for r in comparisons)
        summary['reverse_exact_after_trace_diagnostic'] = sum(r['linear_mfe_in_spark']['exact_target_returned'] or
            r.get('reverse_trace_diagnostic',{}).get('exact_target_returned',False) for r in comparisons)
        assert summary['source_after'] == before and summary['linear_source_after'] == linear_before
        assert summary['linear_artifacts_after'] == frozen_linear
        assert digest(build/'Spark') == lock['binary_sha256']
        assert digest(parameter) == lock['parameter_files'][str(parameter)]
    summary['raw_archive_sha256'] = digest(args.out/'raw_outputs.tar.gz')
    write_json(args.out/'summary.json', summary)
    print(json.dumps({k: summary[k] for k in ['categories','retains_scaffold','spark_extension_reachable_in_linear',
                                             'reverse_exact_target_returned','reverse_exact_after_trace_diagnostic']}, indent=2))


if __name__ == '__main__':
    main()
