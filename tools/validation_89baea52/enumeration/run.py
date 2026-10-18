#!/usr/bin/env python3
"""Serial enumeration only; per-run files live in node-local scratch."""
import argparse
import atexit
import shutil
import collections
import datetime
import hashlib
import json
import math
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
HERE = ROOT / "results/local/e0-enumeration-89baea52-20261018"
OLD = TOOL.parent / 'baseline_15aa7742'
SRC = ROOT / 'deps/LinearCParty-89baea52'
BUILD = ROOT / 'build/linear-cparty-89baea52'
LOCK = ROOT / 'manifests/software-pilot-89baea52.json'
COMMIT = '89baea52c9435b18829ef0867a66a24be56d70f9'
ERRORS = ['logZ_error', 'grouping_logZ_error', 'energy_error', 'bpp_error',
          'unpaired_error', 'context_error', 'access_error']
ERRORS += [f'context_error_{c}' for c in 'SHBIMEP']
ERRORS += [f'access_error_{w}' for w in [1, 2, 3, 5]]
ZEROS = ['ambiguous', 'context_ambiguous', 'ownership_bad', 'ownership_duplicate',
         'ownership_missing', 'ownership_on_paired', 'allowed_absent', 'overflow']
TOTALS = ['structures', 'derivations', 'crossing_structures', 'candidate_superset',
          'candidate_rejected', *ZEROS]


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def git(*args):
    return subprocess.check_output(['git', '-C', str(SRC), *args], text=True).strip()


def fingerprint():
    paths = [LOCK, SRC / 'params/rna_Turner04.par', BUILD / 'CMakeCache.txt',
             TOOL / 'enumerate.cc', TOOL / 'run.py', OLD / 'enumerate.cc',
             OLD / 'inputs/cases.json', OLD / 'report.json', OLD / 'summary.json']
    paths += sorted(BUILD.glob('*.a'))
    paths += [BUILD / name for name in ['LinearCParty', 'LinearPKCapR', 'LinearPKRaccess']]
    return dict(commit=git('rev-parse', 'HEAD'), status=git('status', '--porcelain'),
                sha256={str(p.relative_to(ROOT)): sha(p) for p in paths})


def parse(text):
    values, structures, ablated = {}, {}, {}
    for line in text.splitlines():
        c = line.split('\t')
        if c[0] == 'value':
            values[c[1]] = float(c[2]) if any(x in c[2] for x in '.eE') else int(c[2])
        elif c[0] == 'structure':
            assert c[1] not in structures
            structures[c[1]] = dict(count=int(c[2]), energy=float(c[6]), labels=c[7],
                                     best_log_weight=float(c[4]), worst_log_weight=float(c[5]))
        elif c[0] == 'pre_wip_structure':
            assert c[1] not in ablated
            ablated[c[1]] = int(c[2])
    assert all(math.isfinite(x) for x in values.values())
    return values, structures, ablated


def main():
    HERE.mkdir(parents=True, exist_ok=False)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scratch', type=Path)
    parser.add_argument('--reuse-binary', action='store_true')
    args = parser.parse_args()
    scratch = args.scratch or Path(tempfile.mkdtemp(prefix='verify-d9804df6-'))
    scratch.mkdir(exist_ok=True)
    if str(scratch.resolve()).startswith("/gs/bs/"):
        raise ValueError("scratch must be node-local")
    if args.scratch is None:
        atexit.register(shutil.rmtree, scratch)
    env = {k: v for k, v in os.environ.items() if not k.startswith('CPARTY_')}
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               CPARTY_PF_OVERRIDE_SCALE='1', TMPDIR=str(scratch),
               CPARTY_FIXED_ENERGY_PARAM_FILE=str(SRC / 'params/rna_Turner04.par'))
    before = fingerprint()
    lock = json.loads(LOCK.read_text())
    assert before['commit'] == lock['source_commit'] == COMMIT and not before['status']
    assert sha(SRC / 'params/rna_Turner04.par') == lock['energy_parameter_sha256']
    assert sha(BUILD / 'CMakeCache.txt') == lock['cmake_cache_sha256']
    for name, expected in lock['binaries'].items():
        assert sha(BUILD / name) == expected
    cases = json.loads((OLD / 'inputs/cases.json').read_text())
    assert len(cases) == 729 and len({c['id'] for c in cases}) == 729
    old_rows = {r['id']: r for r in json.loads((OLD / 'report.json').read_text())['enumeration']}
    assert set(old_rows) == {c['id'] for c in cases}
    assert all(all(old_rows[c['id']][k] == c[k] for k in ['sequence', 'scaffold']) for c in cases)
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    manifest = dict(started_utc=started, node=platform.node(), scratch=str(scratch),
                    before=before, environment_overrides={k: env[k] for k in env if k.startswith('CPARTY_') or k in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'TMPDIR']},
                    execution='serial', timeout_s=180, scope='enumeration and WIP-boundary ablation only')
    dump(HERE / 'manifest.json', manifest)
    binary = scratch / 'enumerate'
    if not args.reuse_binary:
        command = ['g++', '-O1', '-std=c++17', '-I' + str(SRC / 'src'), str(TOOL / 'enumerate.cc'),
                   str(BUILD / 'libcparty_scfg2_beam.a'), str(BUILD / 'libcparty_scfg2_exact.a'),
                   str(BUILD / 'libRNA.a'), '-o', str(binary)]
        p = subprocess.run(command, cwd=ROOT, env=env, text=True, capture_output=True, timeout=300)
        (scratch / 'compile.log').write_text(p.stdout + p.stderr)
        dump(scratch / 'compile.json', dict(command=command, returncode=p.returncode))
        if p.returncode:
            raise RuntimeError('compile failed: ' + p.stderr)
    compile_record = json.loads((scratch / 'compile.json').read_text())
    assert compile_record['returncode'] == 0
    manifest.update(compile=compile_record, binary_sha256=sha(binary))
    dump(HERE / 'manifest.json', manifest)
    failures, rows, deltas, witnesses = [], [], [], {}
    process_count = 0
    with tarfile.open(HERE / 'raw.tar.gz', 'w:gz') as archive, (HERE / 'runs.jsonl').open('w') as runlog:
        for name in ['enumerate', 'compile.log', 'compile.json']:
            archive.add(scratch / name, arcname='build/' + name)
        for name in ['initial_link_attempt.log', 'initial_link_attempt.json']:
            if (scratch / name).exists():
                archive.add(scratch / name, arcname='build/' + name)

        def run(case, mode='default', pre_wip=False):
            nonlocal process_count
            label = case['id'] + ('_pre_wip' if pre_wip else '_' + mode)
            command = [str(binary), case['sequence'], case['scaffold'], mode]
            if pre_wip:
                command.append('pre-wip')
            out, err = scratch / 'stdout.tsv', scratch / 'stderr.txt'
            start = time.monotonic()
            timed_out = False
            with out.open('w') as fo, err.open('w') as fe:
                proc = subprocess.Popen(command, cwd=ROOT, env=env, stdout=fo, stderr=fe, start_new_session=True)
                try:
                    proc.wait(timeout=180)
                except subprocess.TimeoutExpired:
                    timed_out = True
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
            record = dict(label=label, command=command, returncode=proc.returncode,
                          timed_out=timed_out, elapsed_s=time.monotonic() - start,
                          stdout_member=label + '.tsv', stderr_member=label + '.stderr',
                          stdout_sha256=sha(out), stderr_sha256=sha(err))
            for path, member in [(out, record['stdout_member']), (err, record['stderr_member'])]:
                archive.add(path, arcname=member)
            runlog.write(json.dumps(record, allow_nan=False) + '\n'); runlog.flush()
            process_count += 1
            if proc.returncode or timed_out or err.stat().st_size:
                failures.append(dict(reason='process failure, timeout, or stderr', **record))
                raise RuntimeError(f'failed {label}: {err.read_text()}')
            return parse(out.read_text())

        for idx, case in enumerate(cases):
            values, structures, _ = run(case)
            if any(values[k] for k in ZEROS) or values['max_count'] != 1 or any(values[k] > 1e-9 for k in ERRORS):
                failures.append(dict(label=case['id'], reason='enumeration or numerical check'))
            assert len(structures) == values['structures']
            assert sum(x['count'] for x in structures.values()) == values['derivations']
            row = dict(case, **values)
            rows.append(row)
            previous = old_rows[case['id']]
            pre, _, pre_structures = run(case, pre_wip=True)
            lost = sorted(pre_structures.keys() - structures.keys())
            added = sorted(structures.keys() - pre_structures.keys())
            restored = all(pre[k] == previous[k] for k in ['structures', 'derivations'])
            if not restored or lost or any(structures[gp]['count'] != c for gp, c in pre_structures.items() if gp in structures):
                failures.append(dict(label=case['id'], reason='WIP ablation does not explain prior counts'))
            delta = dict(id=case['id'], sequence=case['sequence'], scaffold=case['scaffold'],
                         family=case['family'], old_structures=int(previous['structures']),
                         new_structures=values['structures'], delta_structures=values['structures']-int(previous['structures']),
                         old_crossing_structures=int(previous['crossing_structures']),
                         new_crossing_structures=values['crossing_structures'],
                         pre_wip_structures=pre['structures'], pre_wip_derivations=pre['derivations'],
                         pre_wip_matches_old_counts=restored, lost_by_current=len(lost),
                         added_by_wip=len(added), added_structures=added)
            deltas.append(delta)
            if case['id'] == 'closure14':
                witnesses['default'] = dict(values=values, target=structures['.[....][....].'])
            if idx % 25 == 0 or idx == len(cases)-1:
                print(f'{idx+1}/729 structures={sum(r["structures"] for r in rows)} processes={process_count} failures={len(failures)}', flush=True)
        case = next(c for c in cases if c['id'] == 'closure14')
        values, structures, _ = run(case, mode='legacy')
        witnesses['legacy'] = dict(values=values, target=structures['.[....][....].'])
        if witnesses['default']['target']['count'] != 1 or witnesses['legacy']['target']['count'] != 2:
            failures.append(dict(label='closure14', reason='witness counts'))
    after = fingerprint()
    if before != after:
        failures.append(dict(label='integrity', reason='protected inputs, source, build, or tools changed'))
    summary = dict(cases=len(cases), unique_inputs=len({(c['sequence'], c['scaffold']) for c in cases}),
                   min_length=min(len(c['sequence']) for c in cases), max_length=max(len(c['sequence']) for c in cases),
                   families=dict(collections.Counter(c['family'] for c in cases)),
                   totals={k: sum(r[k] for r in rows) for k in TOTALS},
                   max_structures=max(r['structures'] for r in rows),
                   max_derivation_count=max(r['max_count'] for r in rows),
                   maxima={k: max(r[k] for r in rows) for k in ERRORS},
                   maxima_cases={k: [r['id'] for r in rows if r[k] == max(x[k] for x in rows)] for k in ERRORS},
                   previous_totals={k: sum(int(r[k]) for r in old_rows.values()) for k in ['structures', 'derivations', 'crossing_structures']},
                   changed_inputs=sum(d['delta_structures'] != 0 for d in deltas),
                   changed_unique_inputs=len({(d['sequence'], d['scaffold']) for d in deltas if d['delta_structures']}),
                   all_pre_wip_counts_match_old=all(d['pre_wip_matches_old_counts'] for d in deltas),
                   added_by_wip=sum(d['added_by_wip'] for d in deltas),
                   lost_by_current=sum(d['lost_by_current'] for d in deltas),
                   witnesses=witnesses, processes=process_count, failures=failures,
                   verification_complete=not failures, integrity_ok=before == after)
    dump(HERE / 'report.json', dict(enumeration=rows, comparison=deltas))
    dump(HERE / 'summary.json', summary)
    with (HERE / 'counts_vs_15aa7742.tsv').open('w') as f:
        columns = [k for k in deltas[0] if k != 'added_structures']
        f.write('\t'.join(columns) + '\n')
        for d in deltas:
            f.write('\t'.join(str(d[k]) for k in columns) + '\n')
    manifest.update(after=after, finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    processes=process_count, artifacts_sha256={p.name: sha(p) for p in [HERE / name for name in ['raw.tar.gz', 'runs.jsonl', 'report.json', 'summary.json', 'counts_vs_15aa7742.tsv']]})
    dump(HERE / 'manifest.json', manifest)
    print(json.dumps({k: summary[k] for k in ['totals', 'maxima', 'changed_inputs', 'verification_complete']}, ensure_ascii=False), flush=True)
    for name in ['stdout.tsv', 'stderr.txt']:
        (scratch / name).unlink(missing_ok=True)
    return 1 if failures else 0


if __name__ == '__main__':
    raise SystemExit(main())
