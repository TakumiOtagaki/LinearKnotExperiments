#!/usr/bin/env python3
"""Recompute objectives and sampling diagnostics directly from archived outputs."""
import collections
import json
import math
import os
from pathlib import Path
import re
import tarfile
os.environ['OMP_NUM_THREADS'] = os.environ['OPENBLAS_NUM_THREADS'] = os.environ['MKL_NUM_THREADS'] = '1'
import scipy
from scipy.special import logsumexp
from scipy.stats import chi2
from run import HERE, ROOT, CASES, dump


def pairs(s):
    stacks = {x: [] for x in '([{<'}
    closing = dict(zip(')]}>', '([{<'))
    result = set()
    for i, ch in enumerate(s, 1):
        if ch in stacks:
            stacks[ch].append(i)
        elif ch in closing:
            a = stacks[closing[ch]].pop(); result.add((a, i))
        else:
            assert ch == '.'
    assert not any(stacks.values())
    return result


def parse(text):
    sections, mfes = {}, {}
    sec = None
    for line in text.splitlines():
        c = line.split('\t')
        if c[0] == 'section':
            sec = dict(structures={}, values={}, bpp={}, mea={}, centroid={})
            sections[c[1]] = sec
        elif c[0] == 'structure':
            assert c[1] not in sec['structures']
            sec['structures'][c[1]] = dict(count=int(c[2]), log_weight=float(c[3]),
                                           energy=None if c[4] == 'nan' else float(c[4]))
        elif c[0] == 'value':
            sec['values'][c[1]] = float(c[2])
        elif c[0] == 'bpp':
            sec['bpp'][int(c[1]), int(c[2])] = float(c[3])
        elif c[0] == 'mfe':
            mfes[c[1]] = dict(energy=float(c[2]), generated=canonical_extension(c[3]), generated_raw=c[3], union=c[4])
        elif c[0] in ('mea', 'centroid'):
            sec[c[0]][float(c[1])] = dict(score=float(c[2]), generated=canonical_extension(c[3]), generated_raw=c[3], union=c[4], logZ=float(c[5]))
        else:
            raise ValueError(line)
    return sections, mfes


def canonical_extension(s):
    result = ['.'] * len(s)
    for i, j in pairs(s):
        result[i-1], result[j-1] = '[', ']'
    return ''.join(result)


def distribution(sec, scaffold):
    structs = sec['structures']
    z = float(logsumexp([v['log_weight'] for v in structs.values()]))
    p = {gp: math.exp(v['log_weight'] - z) for gp, v in structs.items()}
    fixed = pairs(scaffold)
    union = {gp: fixed | pairs(gp) for gp in structs}
    bpp = collections.defaultdict(float)
    q = [0.] * len(scaffold)
    for gp, prob in p.items():
        paired = set()
        assert not (set(t for pair in pairs(gp) for t in pair) & set(t for pair in fixed for t in pair))
        for pair in union[gp]:
            bpp[pair] += prob; paired.update(pair)
        for i in range(1, len(scaffold)+1):
            if i not in paired:
                q[i-1] += prob
    bpp_error = max((abs(bpp[k]-sec['bpp'].get(k, 0.)) for k in bpp.keys() | sec['bpp'].keys()), default=0.)
    logz_error = max(abs(z-sec['values'][k]) for k in ['api_logZ', 'enumerated_logZ'] + (['outside_logZ'] if 'outside_logZ' in sec['values'] else []))
    return dict(p=p, union=union, bpp=bpp, q=q, z=z, bpp_error=bpp_error, logz_error=logz_error)


def score(pairset, d, gamma):
    paired = {i for pair in pairset for i in pair}
    return 2*gamma*sum(d['bpp'][pair] for pair in pairset) + sum(q for i, q in enumerate(d['q'], 1) if i not in paired)


def centroid_score(pairset, d, gamma):
    # Full expected gain, including TN for every unselected i<j.
    n = len(d['q'])
    return sum(gamma*d['bpp'].get((i, j), 0.) if (i, j) in pairset
               else 1.-d['bpp'].get((i, j), 0.)
               for i in range(1, n+1) for j in range(i+1, n+1))


def chisquare(p, observed, n):
    bins, pooled = [], [0., 0]
    for gp, prob in p.items():
        expected = n*prob
        if expected < 5:
            pooled[0] += expected; pooled[1] += observed[gp]
        else:
            bins.append([expected, observed[gp]])
    if pooled[0] > 0:
        if pooled[0] >= 5 or not bins:
            bins.append(pooled)
        else:
            target = min(bins, key=lambda x: x[0]); target[0] += pooled[0]; target[1] += pooled[1]
    df = len(bins)-1
    statistic = sum((o-e)**2/e for e, o in bins)
    return dict(chi_square=statistic, df=df, p_value=float(chi2.sf(statistic, df)) if df else None,
                pooled_structures=sum(n*x < 5 for x in p.values()),
                minimum_bin_expected=min(x[0] for x in bins), bins=bins)


def collect():
    config = json.loads((HERE / 'config.json').read_text())
    tol = config['tolerance']
    cases = json.loads(CASES.read_text())
    selected = set(c['id'] for c in config['sample_cases'])
    rows, lookup, references, cli_raw, samples_raw = [], {}, {}, {}, {}
    with tarfile.open(HERE / 'raw.tar.gz', 'r:gz') as tar:
        for member in tar:
            if not member.isfile():
                continue
            if member.name.startswith('enumeration/') and member.name.endswith('.stdout'):
                cid = Path(member.name).stem
                case = next(c for c in cases if c['id'] == cid)
                sections, mfes = parse(tar.extractfile(member).read().decode())
                exact = sections['exact']['structures']
                minimum = min(x['energy'] for x in exact.values())
                ties = [gp for gp, x in exact.items() if abs(x['energy']-minimum) <= tol]
                row = dict(case, structures=len(exact), mfe_minimum=minimum, mfe_ties=len(ties), mfe={}, mea={}, centroid={}, ensembles={})
                for mode, m in mfes.items():
                    assert m['generated'] in exact
                    assert pairs(m['union']) == pairs(case['scaffold']) | pairs(m['generated'])
                    gap = m['energy'] - minimum
                    evaluated_error = abs(m['energy'] - exact[m['generated']]['energy'])
                    assert gap >= -tol and evaluated_error <= tol
                    if mode == 'exact':
                        assert abs(gap) <= tol and m['generated'] in ties
                    row['mfe'][mode] = dict(**m, gap=gap, evaluated_energy_error=evaluated_error,
                                            is_global_optimum=abs(gap)<=tol)
                for mode, sec in sections.items():
                    d = distribution(sec, case['scaffold'])
                    assert d['logz_error'] <= tol and d['bpp_error'] <= tol
                    assert all(x['count'] == 1 for x in sec['structures'].values())
                    assert sec['structures'].keys() <= exact.keys()
                    weight_error = max(abs(v['log_weight']-exact[gp]['log_weight']) for gp, v in sec['structures'].items())
                    assert weight_error <= tol
                    if mode == 'exact':
                        rt = sec['values']['rt_kcal_mol']
                        assert all(abs(v['log_weight'] + v['energy']/rt) <= tol for v in exact.values())
                    row['ensembles'][mode] = dict(structures=len(sec['structures']), logZ=d['z'],
                                                  logz_error=d['logz_error'], bpp_error=d['bpp_error'],
                                                  weight_error=weight_error, equals_exact=sec['structures'].keys() == exact.keys())
                    row['mea'][mode] = {}
                    for gamma, m in sec['mea'].items():
                        scores = {gp: score(ps, d, gamma) for gp, ps in d['union'].items()}
                        best = max(scores.values())
                        gp = m['generated']
                        assert gp in scores and pairs(m['union']) == d['union'][gp]
                        # Independently check the Methods incremental form.
                        increment_score = 2*gamma*sum(d['bpp'][pair] for pair in pairs(case['scaffold'])) + sum(d['q'])
                        increment_score += sum(2*gamma*d['bpp'][i,j]-d['q'][i-1]-d['q'][j-1] for i,j in pairs(gp))
                        errors = dict(optimum_error=abs(best-scores[gp]), reported_score_error=abs(m['score']-scores[gp]),
                                      reported_optimum_error=abs(m['score']-best), incremental_error=abs(increment_score-scores[gp]))
                        assert max(errors.values()) <= tol and abs(m['logZ']-d['z']) <= tol
                        row['mea'][mode][str(gamma)] = dict(**m, maximum=best, rescored=scores[gp],
                                                           ties=sum(abs(x-best)<=tol for x in scores.values()), **errors)
                    row['centroid'][mode] = {}
                    assert set(sec['centroid']) == set(config['gamma'])
                    for gamma, m in sec['centroid'].items():
                        scores = {gp: centroid_score(ps, d, gamma) for gp, ps in d['union'].items()}
                        best = max(scores.values())
                        gp = m['generated']
                        assert gp in scores and pairs(m['union']) == d['union'][gp]
                        # Pair-universe constant plus predicted-pair increments.
                        increment_score = len(case['sequence'])*(len(case['sequence'])-1)/2 - sum(d['bpp'].values())
                        increment_score += sum((gamma+1)*d['bpp'][pair]-1 for pair in pairs(m['union']))
                        errors = dict(optimum_error=abs(best-scores[gp]), reported_score_error=abs(m['score']-scores[gp]),
                                      reported_optimum_error=abs(m['score']-best), incremental_error=abs(increment_score-scores[gp]))
                        assert max(errors.values()) <= tol and abs(m['logZ']-d['z']) <= tol
                        row['centroid'][mode][str(gamma)] = dict(**m, maximum=best, rescored=scores[gp],
                            ties=sum(abs(x-best)<=tol for x in scores.values()), **errors)
                    if cid in selected and mode == 'exact':
                        references[cid, mode] = (sec, d)
                rows.append(row); lookup[cid] = row
            elif member.name.startswith('reference/') and member.name.endswith('.stdout'):
                cid, mode = Path(member.name).stem.rsplit('_', 1)
                case = next(c for c in cases if c['id'] == cid)
                sections, _ = parse(tar.extractfile(member).read().decode())
                sec = sections[mode]
                d = distribution(sec, case['scaffold'])
                assert d['logz_error'] <= tol and d['bpp_error'] <= tol
                assert all(v['count'] == 1 for v in sec['structures'].values())
                full = references[cid, 'exact'][0]['structures']
                assert sec['structures'].keys() <= full.keys()
                assert max(abs(v['log_weight']-full[gp]['log_weight']) for gp,v in sec['structures'].items()) <= tol
                references[cid,mode] = (sec,d)
            elif member.name.startswith('sampling/') and member.name.endswith('.stdout'):
                samples_raw[Path(member.name).stem] = collections.Counter(x.strip() for x in tar.extractfile(member).read().decode().splitlines()[2:] if x.strip())
            elif member.name.startswith('cli/') and member.name.endswith('.stdout'):
                cli_raw[Path(member.name).stem] = tar.extractfile(member).read().decode()
    assert len(rows) == len(cases) == 729 and set(lookup) == set(c['id'] for c in cases)
    cli_checks, cli_mfe_max_error, cli_mea_max_error = 0, 0., 0.
    cli_centroid_max_error = 0.
    for label, text in cli_raw.items():
        kind, mode, *gamma_raw = label.split('_')
        chunks = text.strip().split('>')
        count = 0
        for chunk in chunks:
            if not chunk:
                continue
            lines = chunk.strip().splitlines()
            row = lookup[lines[0]]
            assert lines[1:3] == [row['sequence'], row['scaffold']]
            rendered, value = lines[3].split(' ', 1)
            if kind == 'mfe':
                m = row['mfe'][mode]
                reported = float(value.strip('()'))
                err = abs(reported-m['energy'])
                cli_mfe_max_error = max(cli_mfe_max_error, err)
                assert err <= 0.0050000001
            else:
                assert kind in ('mea', 'centroid')
                gamma = str(float(gamma_raw[0])); m = row[kind][mode][gamma]
                reported = float(re.search(r'E[AG]=([^;]+)', value).group(1))
                err = abs(reported-m['score'])
                if kind == 'mea':
                    cli_mea_max_error = max(cli_mea_max_error, err)
                else:
                    cli_centroid_max_error = max(cli_centroid_max_error, err)
                assert err <= tol
            assert pairs(rendered) == pairs(m['union'])
            count += 1
        assert count == 729
        cli_checks += count
    assert len(cli_raw) == 21
    sampling = []
    for mi, mode in enumerate(config['sample_modes']):
        for ci, case in enumerate(config['sample_cases']):
            sec,d = references[case['id'],mode]
            observed_union = samples_raw[case['id']+'_'+mode]
            by_pairset = {frozenset(ps):gp for gp,ps in d['union'].items()}
            observed = collections.Counter()
            for union, count in observed_union.items():
                assert frozenset(pairs(union)) in by_pairset
                observed[by_pairset[frozenset(pairs(union))]] += count
            n = sum(observed.values()); assert n == config['draws_per_condition']
            errors = {gp:abs(observed[gp]/n-p) for gp,p in d['p'].items()}
            gp = max(errors, key=errors.get); p = d['p'][gp]
            stderr = math.sqrt(p*(1-p)/n)
            z_common = [errors[k]/math.sqrt(v*(1-v)/n) for k,v in d['p'].items() if n*v >= 5 and n*(1-v) >= 5]
            crossing = {k for k,ps in d['union'].items() if any(a<c<b<e or c<a<e<b for a,b in ps for c,e in ps)}
            sampling.append(dict(id=case['id'], mode=mode, n=n, seed=config['sample_base_seed']+1000*mi+ci,
                                 structures=len(d['p']), logZ=d['z'], logz_error=d['logz_error'], bpp_error=d['bpp_error'],
                                 total_variation=0.5*sum(errors.values()), max_abs_difference=errors[gp],
                                 max_difference_structure=gp, max_difference_probability=p,
                                 max_difference_standard_error=stderr, max_difference_z=errors[gp]/stderr if stderr else 0.,
                                 max_common_structure_z=max(z_common,default=0.),
                                 pruned=len(sec['structures']) < len(references[case['id'],'exact'][0]['structures']),
                                 weight_error=max(abs(v['log_weight']-references[case['id'],'exact'][0]['structures'][k]['log_weight']) for k,v in sec['structures'].items()),
                                 crossing_structures=len(crossing), crossing_draws=sum(observed[k] for k in crossing),
                                 frequencies={k:dict(probability=v, observed=observed[k]) for k,v in d['p'].items()},
                                 **chisquare(d['p'], observed,n)))
    tested = sorted([r for r in sampling if r['p_value'] is not None], key=lambda r:r['p_value'])
    adjusted = 0.
    for i,r in enumerate(tested):
        adjusted = max(adjusted, min(1., (len(tested)-i)*r['p_value']))
        r['holm_p_value'] = adjusted; r['holm_rejected'] = adjusted <= config['family_alpha']
    for r in sampling:
        if r['p_value'] is None:
            r.update(holm_p_value=None, holm_rejected=False)
    events = sum(r['structures'] for r in sampling)
    eps = math.sqrt(math.log(2*events/config['family_alpha'])/(2*config['draws_per_condition']))
    assert all(r['max_abs_difference'] <= eps for r in sampling)
    def maxfield(items, field):
        v = max(r[field] for r in items)
        return dict(value=v, cases=[r.get('id',r.get('mode')) for r in items if r[field]==v])
    mfe = {}
    for mode in config['mfe_modes']:
        values = [r['mfe'][mode] for r in rows]
        mfe[mode] = dict(cases=len(values), matched=sum(v['is_global_optimum'] for v in values),
                         min_gap=min(v['gap'] for v in values), max_gap=max(v['gap'] for v in values),
                         max_abs_gap=max(abs(v['gap']) for v in values),
                         max_evaluated_energy_error=max(v['evaluated_energy_error'] for v in values),
                         nonoptimal_ids=[r['id'] for r in rows if not r['mfe'][mode]['is_global_optimum']])
    mea = {}
    for mode in config['mea_modes']:
        values = [dict(id=r['id'],**m) for r in rows for m in r['mea'][mode].values()]
        mea[mode] = dict(comparisons=len(values),
                         maxima={k:maxfield(values,k) for k in ['optimum_error','reported_score_error','reported_optimum_error','incremental_error']},
                         tied_comparisons=sum(v['ties']>1 for v in values),
                         max_ties=max(v['ties'] for v in values),
                         pruned_cases=sum(not r['ensembles'][mode]['equals_exact'] for r in rows))
    centroid = {}
    for mode in config['centroid_modes']:
        values = [dict(id=r['id'], **m) for r in rows for m in r['centroid'][mode].values()]
        centroid[mode] = dict(comparisons=len(values),
            maxima={k:maxfield(values,k) for k in ['optimum_error','reported_score_error','reported_optimum_error','incremental_error']},
            tied_comparisons=sum(v['ties']>1 for v in values), max_ties=max(v['ties'] for v in values),
            pruned_cases=sum(not r['ensembles'][mode]['equals_exact'] for r in rows))
    sample_summary = {}
    for mode in config['sample_modes']:
        values = [r for r in sampling if r['mode']==mode]
        sample_summary[mode] = dict(conditions=len(values), total_draws=sum(r['n'] for r in values),
                                    total_structures=sum(r['structures'] for r in values),
                                    pruned_conditions=sum(r['pruned'] for r in values),
                                    min_p_value=min((r['p_value'] for r in values if r['p_value'] is not None),default=None),
                                    holm_rejections=sum(r['holm_rejected'] for r in values),
                                    maxima={k:maxfield(values,k) for k in ['total_variation','max_abs_difference','max_difference_z','max_common_structure_z','logz_error','bpp_error','weight_error']})
    summary = dict(version=config['version'], cases=len(rows), unique_inputs=len({(r['sequence'],r['scaffold']) for r in rows}),
                   min_length=min(len(r['sequence']) for r in rows), max_length=max(len(r['sequence']) for r in rows),
                   total_structures=sum(r['structures'] for r in rows), tolerance=tol,
                   mfe=mfe, mfe_tied_inputs=sum(r['mfe_ties']>1 for r in rows),
                   max_mfe_ties=max(r['mfe_ties'] for r in rows), mea=mea, centroid=centroid,
                   cli=dict(comparisons=cli_checks, mfe_max_reported_error=cli_mfe_max_error, mea_max_reported_error=cli_mea_max_error, centroid_max_reported_error=cli_centroid_max_error),
                   ensemble_max_logz_error=max(e['logz_error'] for r in rows for e in r['ensembles'].values()),
                   ensemble_max_bpp_error=max(e['bpp_error'] for r in rows for e in r['ensembles'].values()),
                   sampling=dict(conditions=len(sampling), draws=sum(r['n'] for r in sampling), modes=sample_summary,
                                 family_alpha=config['family_alpha'], frequency_events=events, hoeffding_epsilon=eps,
                                 tested_conditions=len(tested), holm_rejections=sum(r['holm_rejected'] for r in sampling)),
                   scientific_failures=[]) # Assertions above stop output on any deterministic failure.
    return dict(decoding=rows, sampling=sampling), summary


if __name__ == '__main__':
    report, summary = collect()
    dump(HERE / 'report.json', report)
    dump(HERE / 'summary.json', summary)
    print(json.dumps(summary,ensure_ascii=False,indent=2))
