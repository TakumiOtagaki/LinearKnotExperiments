#!/usr/bin/env python3
"""Serial small-input validation; four shared result files, scratch per call.

Run from the paper root:
  python3 -B experiments/scripts/validate_sampling.py --out experiments/results/local/sampling-validation-20261012
The historical oracle is an explicit external input; --oracle and --cases allow
the experiment repository to use copies at other locations.
"""
import argparse
from collections import Counter
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import tarfile
import tempfile
import time

from scipy.stats import chi2

EXP = Path(__file__).resolve().parents[1]
PAPER = EXP.parent
COMMIT = '15aa7742a5497164573294e37972ff7ccf9ddfe1'


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def git(source):
    def call(*args):
        return subprocess.check_output(['git', '-C', str(source), *args], text=True).strip()
    return dict(commit=call('rev-parse', 'HEAD'), status=call('status', '--porcelain'))


def pairs(db):
    stack, result = [], []
    for i, ch in enumerate(db, 1):
        if ch in '([':
            stack.append(i)
        elif ch in ')]':
            result.append((stack.pop(), i))
    if stack:
        raise ValueError('unbalanced extension')
    return result


def choose_cases(path):
    historical = json.loads(path.read_text())
    indexed = {x['id']: x for x in historical}
    selected = [
        dict(id='seven_empty', sequence='GGAAACC', scaffold='.......', family='control'),
        dict(id='nine_empty', sequence='GGGAAACCC', scaffold='.........', family='control'),
        dict(id='nine_fixed', sequence='GGGAAACCC', scaffold='(((...)))', family='control'),
        dict(id='crossing20', sequence='GGGGAGGGAACCCCCCCAAA', scaffold='((((......))))......', family='crossing_witness'),
        indexed['closure14'],
        dict(id='twenty_empty', sequence='GCGCGCGCGCGCGCGCGCGC', scaffold='....................', family='control'),
    ]
    selected += [indexed[f'pk{i:03d}'] for i in range(14)]
    for family in ['random', 'fresh_random_17_20']:
        selected += [x for x in historical if x['family'] == family and '(' in x['scaffold']][:4]
    selected += [x for x in historical if x['family'] == 'single_pair_grid' and len(x['sequence']) == 12][:4]
    if len(selected) != 32 or len({(x['sequence'], x['scaffold']) for x in selected}) != 32:
        raise ValueError('cohort must contain 32 distinct S,G inputs')
    return selected


def read_reference(path):
    dist, values, bpp, mfe = {}, {}, {}, {}
    for line in path.read_text().splitlines():
        cells = line.split('\t')
        if cells[0] == 'structure':
            dist[cells[1]] = dict(derivations=int(cells[2]), log_weight=float(cells[3]), energy=float(cells[4]))
        elif cells[0] == 'value':
            values[cells[1]] = float(cells[2])
        elif cells[0] == 'bpp':
            bpp[(int(cells[1]), int(cells[2]))] = float(cells[3])
        elif cells[0] == 'mfe':
            mfe[cells[1]] = dict(energy=float(cells[2]), extension=cells[3], completed_structure=cells[4])
    return dist, values, bpp, mfe


def goodness(counts, probabilities, n):
    # Pool all categories with E<5 using only the reference, before looking at O.
    bins = [[n * p, counts.get(k, 0), [k]] for k, p in sorted(probabilities.items()) if n * p >= 5]
    rare = [[n * p, counts.get(k, 0), k] for k, p in sorted(probabilities.items()) if n * p < 5]
    if rare:
        pooled = [sum(x[0] for x in rare), sum(x[1] for x in rare), [x[2] for x in rare]]
        if pooled[0] >= 5 or not bins:
            bins.append(pooled)
        else:
            ix = min(range(len(bins)), key=lambda i: bins[i][0])
            bins[ix][0] += pooled[0]
            bins[ix][1] += pooled[1]
            bins[ix][2] += pooled[2]
    statistic = sum(2 * o * math.log(o / e) for e, o, _ in bins if o)
    df = len(bins) - 1
    return dict(statistic=max(0., statistic), df=df,
                p_value=float(chi2.sf(max(0., statistic), df)) if df > 0 else None,
                bins=len(bins), min_expected=min(x[0] for x in bins),
                rare_structures=len(rare), pooling='E<5 pooled; pool below 5 merged into smallest common bin',
                asymptotic=True)


def analyze(case, dist, values, api_bpp, counts, n):
    probabilities = {gp: math.exp(v['log_weight'] - values['enumerated_logZ']) for gp, v in dist.items()}
    unexpected = sorted(set(counts) - set(probabilities))
    fixed = pairs(case['scaffold'])
    theory_bpp, empirical_bpp = Counter(), Counter()
    theory_windows, empirical_windows = Counter(), Counter()
    pk_probability = pk_frequency = 0.
    for gp, p in probabilities.items():
        added = pairs(gp)
        empirical = counts.get(gp, 0) / n
        for pair in fixed + added:
            theory_bpp[pair] += p
            empirical_bpp[pair] += empirical
        pk = any(i < k < j < l or k < i < l < j for i, j in added for k, l in fixed)
        if pk:
            pk_probability += p
            pk_frequency += empirical
        paired = {i for bp in fixed + added for i in bp}
        for width in [1, 3, 5]:
            for start in range(1, len(gp) - width + 2):
                if not any(i in paired for i in range(start, start + width)):
                    theory_windows[(start, width)] += p
                    empirical_windows[(start, width)] += empirical
    bpp_keys = {(i, j) for i in range(1, len(case['sequence'])) for j in range(i + 1, len(case['sequence']) + 1)}
    window_keys = {(start, width) for width in [1, 3, 5] for start in range(1, len(case['sequence']) - width + 2)}
    structure_error = max(abs(counts.get(gp, 0) / n - p) for gp, p in probabilities.items())
    max_bpp_error = max(abs(empirical_bpp[k] - theory_bpp[k]) for k in bpp_keys)
    api_error = max(abs(api_bpp.get(k, 0.) - theory_bpp[k]) for k in bpp_keys)
    max_window_error = max(abs(empirical_windows[k] - theory_windows[k]) for k in window_keys)
    result = dict(structures=len(dist), derivations=sum(v['derivations'] for v in dist.values()),
                  max_derivations_per_structure=max(v['derivations'] for v in dist.values()),
                  draws=sum(counts.values()), unexpected_sample_structures=unexpected,
                  total_variation=0.5 * (sum(abs(counts.get(gp, 0) / n - p) for gp, p in probabilities.items())
                                         + sum(counts[gp] / n for gp in unexpected)),
                  max_structure_frequency_error=structure_error,
                  max_bpp_sample_minus_enumeration=max_bpp_error,
                  max_bpp_sample_minus_api=max(abs(empirical_bpp[k] - api_bpp.get(k, 0.)) for k in bpp_keys),
                  max_bpp_enumeration_minus_api=api_error,
                  max_window_frequency_error=max_window_error,
                  pk_probability=pk_probability, pk_frequency=pk_frequency,
                  events=len(dist) + len(bpp_keys) + len(window_keys),
                  g_test=goodness(counts, probabilities, n), reference=values,
                  frequencies=[dict(extension=k, probability=probabilities[k], count=counts.get(k, 0),
                                    derivations=dist[k]['derivations'], log_weight=dist[k]['log_weight'],
                                    energy=dist[k]['energy'] if math.isfinite(dist[k]['energy']) else None)
                               for k in sorted(dist)],
                  bpp=[dict(i=i, j=j, enumerated=theory_bpp[(i,j)], api=api_bpp.get((i,j), 0.),
                            sample=empirical_bpp[(i,j)]) for i, j in sorted(bpp_keys)],
                  windows=[dict(start=i, width=w, enumerated=theory_windows[(i,w)], sample=empirical_windows[(i,w)])
                           for i,w in sorted(window_keys)])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--source', type=Path, default=EXP / 'deps/LinearCParty-15aa7742')
    parser.add_argument('--build', type=Path, default=EXP / 'build/linear-cparty-15aa7742')
    parser.add_argument('--oracle', type=Path, default=PAPER / 'textbook/examples/verify_15aa7742/enumerate.cc')
    parser.add_argument('--cases', type=Path, default=PAPER / 'textbook/examples/verify_15aa7742/inputs/cases.json')
    parser.add_argument('--probe', type=Path, help='Use an already compiled probe (recorded, not rebuilt)')
    parser.add_argument('--preflight-dir', type=Path, help='Archive earlier tool-build or smoke logs')
    args = parser.parse_args()
    out = args.out.resolve()
    if out.exists():
        raise FileExistsError(out)
    source, build = args.source.resolve(), args.build.resolve()
    before = git(source)
    if before != dict(commit=COMMIT, status=''):
        raise RuntimeError(f'fixed clone mismatch: {before}')
    cohort = choose_cases(args.cases)
    parameters = source / 'params/rna_Turner04.par'
    files = [parameters, build / 'CMakeCache.txt', build / 'LinearPKSampling',
             *[build / f for f in ['libcparty_scfg2_beam.a', 'libcparty_scfg2_exact.a', 'libRNA.a']],
             args.oracle, args.cases, Path(__file__), EXP / 'scripts/sampling_validation_probe.cc']
    hashes = {str(p.resolve()): digest(p) for p in files}
    env = {k: v for k, v in os.environ.items() if not k.startswith('CPARTY_')}
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               CPARTY_FIXED_ENERGY_PARAM_FILE=str(parameters), CPARTY_PF_OVERRIDE_SCALE='1', LC_ALL='C')
    quota = subprocess.run(['lfs', 'quota', '-p', '102571', '/gs/bs'], capture_output=True, text=True)
    if quota.returncode:
        raise RuntimeError('quota unavailable; no shared output created')
    config = dict(purpose='development finite validation', fixed_source=before, source=str(source), build=str(build),
                  hashes=hashes, cases=cohort, count=100000, alpha=0.01, seed_base=202610120000,
                  modes=['exact', 'beam1', 'beam3'], beam_caps='n (nonbinding for this cohort)',
                  closure='RequireWeakClosure', temperature_celsius=37, dangles=2,
                  event_rule='every reference extension, every i<j BPP, every width 1/3/5 window',
                  windows=[1, 3, 5], g_test_rule='pool E<5 before examining frequencies; Holm family alpha=.01',
                  thread_environment={k: env[k] for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']},
                  quota_stdout=quota.stdout, new_shared_entries_estimate=5,
                  hostname=platform.node(), python=platform.python_version(),
                  compiler=subprocess.check_output(['g++','--version'], text=True).splitlines()[0],
                  started=time.strftime('%Y-%m-%dT%H:%M:%S%z'))
    out.mkdir(parents=True)
    (out / 'config.json').write_text(json.dumps(config, ensure_ascii=False, indent=2) + '\n')
    rows = []
    with tempfile.TemporaryDirectory(prefix='sampling-validation-', dir='/tmp') as tmp, \
         tarfile.open(out / 'raw.tar.gz', 'w:gz') as archive, \
         open(out / 'results.jsonl', 'w') as jsonl:
        scratch = Path(tmp)
        if args.preflight_dir:
            for path in sorted(args.preflight_dir.iterdir()):
                if path.is_file() and path.name != 'probe':
                    archive.add(path, arcname='preflight/'+path.name, recursive=False)
        if args.probe:
            probe = args.probe.resolve()
        else:
            probe = scratch / 'probe'
            command = ['g++', '-O1', '-std=c++17', f'-DCPARTY_VALIDATION_ORACLE="{args.oracle.resolve()}"',
                       '-I'+str(source / 'src'), str(EXP / 'scripts/sampling_validation_probe.cc'),
                       *[str(build / f) for f in ['libcparty_scfg2_beam.a','libcparty_scfg2_exact.a','libRNA.a']], '-o', str(probe)]
            p = subprocess.run(command, env=env, capture_output=True, timeout=180)
            for name, content in [('compile.stdout',p.stdout),('compile.stderr',p.stderr),
                                  ('compile.command.json', json.dumps(command).encode())]:
                info = tarfile.TarInfo(name);info.size=len(content);archive.addfile(info,io.BytesIO(content))
            if p.returncode:
                raise RuntimeError('probe compile failed; raw log retained')
        config['probe_sha256'] = digest(probe)
        archive.add(probe, arcname='probe', recursive=False)
        archive.add(parameters, arcname='rna_Turner04.par', recursive=False)
        archive.add(args.oracle, arcname='historical_enumerate.cc', recursive=False)
        for i, case in enumerate(cohort):
            for mode_index, mode in enumerate(config['modes']):
                width = 0 if mode == 'exact' else int(mode[4:])
                seed = config['seed_base'] + mode_index * 1000 + i
                draws_path = scratch / 'draws.txt'
                command = [str(probe), case['sequence'], case['scaffold'], 'exact' if width == 0 else 'beam',
                           str(width), str(config['count']), str(seed), str(draws_path)]
                start = time.monotonic()
                with open(scratch / 'reference.tsv', 'wb') as stdout, open(scratch / 'stderr', 'wb') as stderr:
                    try:
                        p = subprocess.run(command, env=env, cwd=scratch, stdout=stdout, stderr=stderr, timeout=180)
                        code = p.returncode
                    except subprocess.TimeoutExpired:
                        code = 124
                prefix = case['id'] + '/' + mode
                raw = {}
                for name in ['reference.tsv', 'stderr', 'draws.txt']:
                    path = scratch / name
                    if path.exists():
                        raw[name] = dict(bytes=path.stat().st_size, sha256=digest(path), archive_member=prefix+'/'+name)
                        archive.add(path, arcname=prefix+'/'+name, recursive=False)
                row = dict(case=case['id'], mode=mode, sequence=case['sequence'], scaffold=case['scaffold'],
                           n=len(case['sequence']), seed=seed, sample_count=config['count'],
                           command=command, exit_code=code, elapsed_seconds=time.monotonic()-start, raw=raw)
                if code == 0:
                    dist, values, bpp, mfe = read_reference(scratch / 'reference.tsv')
                    with open(draws_path) as f:
                        counts = Counter(line.strip() for line in f)
                    row['metrics'] = analyze(case, dist, values, bpp, counts, config['count'])
                    row['mfe_preparation'] = mfe
                    if width:
                        exact_row = next((r for r in rows if r['case']==case['id'] and r['mode']=='exact' and 'metrics' in r), None)
                        if exact_row:
                            full = {x['extension']: x for x in exact_row['metrics']['frequencies']}
                            missing = sorted(set(dist)-set(full))
                            row['metrics']['retained_support_absent_from_exact'] = missing
                            row['metrics']['max_retained_weight_minus_exact_weight'] = max(
                                (abs(v['log_weight']-full[k]['log_weight']) for k,v in dist.items() if k in full), default=0.)
                            row['metrics']['retained_support_structures_removed'] = len(full)-len(dist)
                else:
                    row['error'] = (scratch / 'stderr').read_text()[-4000:]
                rows.append(row)
                jsonl.write(json.dumps(row, ensure_ascii=False, allow_nan=False) + '\n');jsonl.flush()
                for name in ['reference.tsv','stderr','draws.txt']:
                    (scratch / name).unlink(missing_ok=True)
            print(f"{i+1}/32 {case['id']}: " + ', '.join(f"{r['mode']} rc={r['exit_code']}" for r in rows[-3:]), flush=True)
    config['finished'] = time.strftime('%Y-%m-%dT%H:%M:%S%z')
    config['source_after'] = git(source)
    config['hashes_unchanged'] = all(digest(Path(p)) == h for p, h in hashes.items())
    config['archive_sha256'] = digest(out / 'raw.tar.gz')
    (out / 'config.json').write_text(json.dumps(config, ensure_ascii=False, indent=2) + '\n')
    good = [r for r in rows if r['exit_code'] == 0]
    events = sum(r['metrics']['events'] for r in good)
    epsilon = math.sqrt(math.log(2 * events / config['alpha']) / (2 * config['count']))
    # Fixed 96-test family; degenerate one-category tests have p=1.
    def test_p(row):
        p=row.get('metrics', {}).get('g_test', {}).get('p_value')
        return p if p is not None else 1.
    ranked = sorted(rows, key=test_p)
    rejecting = True
    for rank, r in enumerate(ranked):
        p = r.get('metrics', {}).get('g_test', {}).get('p_value')
        reject = rejecting and p is not None and p <= config['alpha'] / (len(rows) - rank)
        rejecting &= reject
        if 'metrics' in r:
            r['metrics']['g_test']['holm_reject'] = reject
    summary = dict(planned=96, completed=len(rows), successful=len(good), failures=[r['case']+'/'+r['mode'] for r in rows if r['exit_code']],
                   count_per_case_mode=config['count'], total_draws=sum(r['metrics']['draws'] for r in good),
                   events=events, alpha=config['alpha'], simultaneous_hoeffding_epsilon=epsilon,
                   hashes_unchanged=config['hashes_unchanged'], source_unchanged=config['source_after']==before,
                   archive_sha256=config['archive_sha256'], modes={})
    numerical_failures = []
    for r in good:
        m = r['metrics'];v = m['reference']
        checks = [abs(v['enumerated_logZ']-v['api_logZ']) <= 1e-9,
                  abs(v['enumerated_logZ']-v['sampling_logZ']) <= 1e-9,
                  m['max_bpp_enumeration_minus_api'] <= 1e-9,
                  m['max_derivations_per_structure'] == 1, m['draws']==config['count'],
                  not m['unexpected_sample_structures']]
        if r['mode']!='exact':
            checks += [not m.get('retained_support_absent_from_exact', ['exact reference missing']),
                       m.get('max_retained_weight_minus_exact_weight', math.inf)<=1e-9]
        if not all(checks):
            numerical_failures.append(r['case']+'/'+r['mode'])
    for mode in config['modes']:
        subset = [r for r in good if r['mode']==mode]
        if not subset:
            summary['modes'][mode] = dict(successful=0);continue
        def extreme(key):
            r=max(subset,key=lambda x:x['metrics'][key]);return dict(value=r['metrics'][key],case=r['case'])
        summary['modes'][mode] = dict(successful=len(subset), total_structures=sum(r['metrics']['structures'] for r in subset),
                                      pk_positive_inputs=sum(r['metrics']['pk_probability']>0 for r in subset),
                                      max_tv=extreme('total_variation'), max_bpp_error=extreme('max_bpp_sample_minus_api'),
                                      max_structure_error=extreme('max_structure_frequency_error'),
                                      max_window_error=extreme('max_window_frequency_error'),
                                      max_api_bpp_reference_error=extreme('max_bpp_enumeration_minus_api'),
                                      holm_rejections=[r['case'] for r in subset if r['metrics']['g_test']['holm_reject']],
                                      min_raw_p=min((r['metrics']['g_test']['p_value'] for r in subset if r['metrics']['g_test']['p_value'] is not None),default=None))
    summary['numerical_failures'] = numerical_failures
    summary['hoeffding_failures'] = [r['case']+'/'+r['mode'] for r in good
                                    if max(r['metrics']['max_structure_frequency_error'],r['metrics']['max_bpp_sample_minus_enumeration'],
                                           r['metrics']['max_window_frequency_error']) > epsilon]
    summary['g_test_holm_rejections'] = [r['case']+'/'+r['mode'] for r in good if r['metrics']['g_test']['holm_reject']]
    summary['all_checks_pass'] = (len(good)==96 and not numerical_failures and not summary['hoeffding_failures']
                                  and not summary['g_test_holm_rejections'] and config['hashes_unchanged'] and summary['source_unchanged'])
    # Preserve adjusted tests as well as unadjusted per-run output.
    with open(out / 'results.jsonl', 'w') as f:
        for r in rows:f.write(json.dumps(r,ensure_ascii=False,allow_nan=False)+'\n')
    (out / 'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(summary,ensure_ascii=False,indent=2))
    return 0 if summary['all_checks_pass'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
