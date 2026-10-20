#!/usr/bin/env python3
"""Serial pinned-version validation; transient files stay in /tmp."""
import collections
import datetime
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import tarfile
import tempfile
import time

TOOL = Path(__file__).resolve().parent
ROOT = TOOL.parents[2]
HERE = ROOT / "results/local/e0-decoding-centroid-89baea52-20261020"
SRC = ROOT / 'deps/LinearCParty-89baea52'
BUILD = ROOT / 'build/linear-cparty-89baea52'
LOCK = ROOT / 'manifests/software-pilot-89baea52.json'
CASES = ROOT / 'tools/validation_89baea52/baseline_15aa7742/inputs/cases.json'
CONFIG = ROOT / 'configs/centroid_89baea52_20261020/e0-decoding-centroid-89baea52-20261020.json'
ORACLE = ROOT / 'tools/validation_89baea52/enumeration/enumerate.cc'
SAMPLE_IDS = ['closure14', 'random000', 'random001', 'random004', 'random031', 'random263']


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def dump(path, obj):
    Path(path).write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def fingerprint():
    paths = [CONFIG, LOCK, CASES, ORACLE, SRC / 'params/rna_Turner04.par', BUILD / 'CMakeCache.txt']
    paths += list(BUILD.glob('*.a'))
    paths += [BUILD / n for n in ['LinearCParty', 'LinearPKMEA', 'LinearPKCentroid', 'LinearPKSampling']]
    paths += [TOOL / n for n in ['probe.cc', 'run.py', 'analyze.py', 'validate.py']]
    return dict(commit=subprocess.check_output(['git', '-C', str(SRC), 'rev-parse', 'HEAD'], text=True).strip(),
                status=subprocess.check_output(['git', '-C', str(SRC), 'status', '--porcelain'], text=True),
                sha256={str(p.relative_to(ROOT)): sha(p) for p in paths})


def main():
    HERE.mkdir(parents=True, exist_ok=False)
    before = fingerprint()
    lock = json.loads(LOCK.read_text())
    assert before['commit'] == lock['source_commit'] == '89baea52c9435b18829ef0867a66a24be56d70f9'
    assert not before['status']
    assert sha(SRC / 'params/rna_Turner04.par') == lock['energy_parameter_sha256']
    assert sha(BUILD / 'CMakeCache.txt') == lock['cmake_cache_sha256']
    assert sha(BUILD / 'LinearCParty') == lock['binaries']['LinearCParty']
    cases = json.loads(CASES.read_text())
    assert len(cases) == 729
    sample_cases = [next(c for c in cases if c['id'] == cid) for cid in SAMPLE_IDS]
    config = json.loads(CONFIG.read_text())
    assert config['version'] == before['commit']
    assert config['sample_cases'] == sample_cases
    assert config['gamma'] == [0.5, 1, 2]
    dump(HERE / 'config.json', config)
    manifest = dict(before=before, config_sha256=sha(HERE / 'config.json'), node=platform.node(),
                    started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    execution='serial, thread 1', job_id=os.environ.get('JOB_ID'),
                    sge_task_id=os.environ.get('SGE_TASK_ID'), estimated_shared_files=12)
    dump(HERE / 'manifest.json', manifest)
    with tempfile.TemporaryDirectory(prefix='verify-centroid-89baea52-', dir=os.environ['TMPDIR']) as tmp:
        scratch = Path(tmp)
        manifest['scratch'] = tmp
        env = {k: v for k, v in os.environ.items() if not k.startswith('CPARTY_')}
        env.update(CPARTY_FIXED_ENERGY_PARAM_FILE=str(SRC / 'params/rna_Turner04.par'),
                   CPARTY_PF_OVERRIDE_SCALE='1', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                   MKL_NUM_THREADS='1', TMPDIR=tmp)
        manifest['environment_overrides'] = {k: env[k] for k in env if k.startswith('CPARTY_') or k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'TMPDIR']}
        members = {}
        records = []
        with tarfile.open(HERE / 'raw.tar.gz', 'w:gz') as archive, (HERE / 'runs.jsonl').open('w') as log:
            def add_bytes(name, data):
                info = tarfile.TarInfo(name); info.size = len(data)
                archive.addfile(info, io.BytesIO(data))
                members[name] = dict(sha256=hashlib.sha256(data).hexdigest(), bytes=len(data))
            def run(label, command, timeout=180):
                start = time.monotonic()
                out, err = scratch / 'stdout', scratch / 'stderr'
                timed_out = False
                with out.open('wb') as fo, err.open('wb') as fe:
                    p = subprocess.Popen([str(x) for x in command], cwd=ROOT, env=env, stdout=fo, stderr=fe, start_new_session=True)
                    try:
                        p.wait(timeout=timeout)
                    except subprocess.TimeoutExpired:
                        timed_out = True; os.killpg(p.pid, signal.SIGKILL); p.wait()
                data, errors = out.read_bytes(), err.read_bytes()
                add_bytes(label + '.stdout', data); add_bytes(label + '.stderr', errors)
                rec = dict(label=label, command=[str(x) for x in command], returncode=p.returncode,
                           timed_out=timed_out, seconds=time.monotonic()-start,
                           stdout_member=label+'.stdout', stderr_member=label+'.stderr')
                records.append(rec); log.write(json.dumps(rec)+'\n'); log.flush()
                if p.returncode or timed_out:
                    raise RuntimeError(label + ': ' + errors.decode(errors='replace'))
                if label != 'build/compile' and errors:
                    raise RuntimeError(label + ': unexpected stderr ' + errors.decode())
                return data.decode()
            for name in ['LinearCParty', 'LinearPKMEA', 'LinearPKCentroid', 'LinearPKSampling']:
                run('help/' + name, [BUILD / name, '--help'])
            binary = scratch / 'probe'
            cmd = ['g++', '-O1', '-std=c++17', '-I'+str(SRC / 'src'),
                   '-DCPARTY_VALIDATION_ORACLE="'+str(ORACLE)+'"', TOOL / 'probe.cc',
                   BUILD / 'libcparty_scfg2_beam.a', BUILD / 'libcparty_scfg2_exact.a',
                   BUILD / 'libRNA.a', '-o', binary]
            run('build/compile', cmd, 600)
            add_bytes('build/probe', binary.read_bytes())
            manifest['binary_sha256'] = sha(binary)
            dump(HERE / 'manifest.json', manifest)
            for i, case in enumerate(cases):
                run('enumeration/' + case['id'], [binary, case['sequence'], case['scaffold'], 'full'])
                if i % 25 == 0 or i == 728:
                    print(f'enumeration {i+1}/729', flush=True)
            records_file = scratch / 'cases.txt'
            records_file.write_text(''.join('>'+c['id']+'\n'+c['sequence']+'\n'+c['scaffold']+'\n' for c in cases))
            add_bytes('cli/input.txt', records_file.read_bytes())
            for mode in ['exact', 'beam20', 'beam100']:
                opts = ['--exact'] if mode == 'exact' else ['-b', mode[4:], '--ranking', 'outer', '--max-hairpin', '30', '--max-case3-depth', '30']
                run('cli/mfe_' + mode, [BUILD / 'LinearCParty', '-i', records_file, *opts], 600)
                for gamma in [0.5, 1, 2]:
                    run('cli/mea_' + mode + '_' + str(gamma), [BUILD / 'LinearPKMEA', '-i', records_file, '--gamma', str(gamma), *opts], 600)
                    run('cli/centroid_' + mode + '_' + str(gamma), [BUILD / 'LinearPKCentroid', '-i', records_file, '--gamma', str(gamma), *opts], 600)
                print('CLI ' + mode + ' complete', flush=True)
            for mi, mode in enumerate(config['sample_modes']):
                for ci, case in enumerate(sample_cases):
                    if mode != 'exact':
                        run('reference/'+case['id']+'_'+mode, [binary, case['sequence'], case['scaffold'], mode[4:]])
                    seed = config['sample_base_seed']+1000*mi+ci
                    opts = ['--exact'] if mode == 'exact' else ['-b', mode[4:], '--ranking', 'outer', '--max-hairpin', '30', '--max-case3-depth', '30']
                    run('sampling/'+case['id']+'_'+mode, [BUILD / 'LinearPKSampling', '--sequence', case['sequence'], '--scaffold', case['scaffold'], '--count', str(config['draws_per_condition']), '--seed', str(seed), *opts])
                    print('sampling ' + case['id'] + ' ' + mode + ' complete', flush=True)
        after = fingerprint()
        manifest.update(after=after, integrity_ok=before == after,
                        finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                        process_count=len(records), members=members,
                        artifacts_sha256={n: sha(HERE / n) for n in ['raw.tar.gz', 'runs.jsonl', 'config.json']})
        dump(HERE / 'manifest.json', manifest)
        assert before == after
    print('all processes complete; protected fingerprints unchanged', flush=True)


if __name__ == '__main__':
    main()
