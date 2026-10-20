#!/usr/bin/env python3
"""Explore frozen 89baea52 examples. No manuscript or external source writes."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile
import time

import numpy as np

from analyze_e3_neat1 import CONTEXTS, pairs_of, parse_context, read_inputs
from analyze_pk_membership import CANONICAL, parse_structure, render
from analyze_small_rerun_89baea52 import groups
from task_storage import RunReader

ROOT = Path(__file__).resolve().parents[1]
REV = "89baea52c9435b18829ef0867a66a24be56d70f9"
SEED = 20261020


def sha(data):
    return hashlib.sha256(data).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n")


def table(path, rows):
    rows = list(rows)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'wt', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def cross_vector(pairs, fixed):
    """Independent coordinate test, vectorized over G, without bracket labels."""
    fixed = np.array(sorted(fixed), dtype=int).reshape(-1, 2)
    if not len(fixed):
        return set()
    k, l = fixed.T
    return {(i, j) for i, j in pairs
            if np.any(((i < k) & (k < j) & (j < l)) |
                      ((k < i) & (i < l) & (l < j)))}


def bpp_parse(data):
    out = {}
    for line in data.decode().splitlines():
        c = line.split()
        if not c or line.startswith(('#', '>')):
            continue
        assert len(c) == 3
        i, j, p = int(c[0]), int(c[1]), float(c[2])
        assert i < j and 0 <= p <= 1 + 1e-8 and (i, j) not in out
        out[i, j] = p
    return out


def provenance_artifact(reader, job, basename, provenance):
    data = reader.artifact(job, basename)
    provenance.append(dict(run=str(reader.path.relative_to(ROOT)),
                           manifest_sha256=reader.manifest_sha256,
                           job=job['id'], member=job['directory']+'/'+basename,
                           sha256=sha(data), bytes=len(data), command=job['command']))
    return data


def figure_search(out):
    out.mkdir(exist_ok=False)
    cohort_path = ROOT/'results/local/pk-recovery-bprna-20261016/pk_sel.json'
    cohort = {r['id']: r for r in json.loads(cohort_path.read_text())
              if 50 <= len(r['seq']) <= 100}
    source, data, rows = [], {}, []
    for reader, job in groups('pk-recovery-200-turner04'):
        if job['input'] not in cohort:
            continue
        assert '--exact' in job['command']
        assert job['command'][job['command'].index('--bpp-cutoff')+1] == '0'
        r = cohort[job['input']]
        seq, n = r['seq'], len(r['seq'])
        fixed = set(map(tuple, r['nested']))
        expected = {(i, j) for i, j, o in parse_structure(r['ref'])
                    if o == '(' and j-i >= 4 and seq[i-1]+seq[j-1] in CANONICAL}
        assert fixed == expected and fixed
        scaffold = render(n, fixed)
        entry = next(e for e in reader.manifest['inputs'] if e['name'] == r['id'])
        frozen = reader.input(entry)
        assert sha(frozen) == entry['sha256']
        assert frozen.decode().splitlines() == ['>'+r['id'], seq, scaffold]
        raw = provenance_artifact(reader, job, 'bpp.tsv', source)
        pairs = bpp_parse(raw)
        crossing = cross_vector(pairs, fixed)
        reference_crossing = cross_vector(
            {(i, j) for i, j, _ in parse_structure(r['ref'])
             if seq[i-1]+seq[j-1] in CANONICAL}, fixed)
        assert reference_crossing
        assert all(abs(pairs[p]-1) < 1e-8 for p in fixed)
        mid = [p for p in crossing if .2 <= pairs[p] <= .8]
        best = min(mid, key=lambda p: (abs(pairs[p]-.4), p)) if mid else None
        row = dict(id=r['id'], n=n, G_pairs=len(fixed),
                   reference_crossing_pairs=len(reference_crossing),
                   crossing_max=max((pairs[p] for p in crossing), default=0.),
                   crossing_mid_n=len(mid),
                   crossing_near_04_i=best[0] if best else None,
                   crossing_near_04_j=best[1] if best else None,
                   crossing_near_04_p=pairs[best] if best else None)
        rows.append(row)
        data[r['id']] = dict(record=r, scaffold=scaffold, pairs=pairs,
                              crossing=crossing, bpp_raw=raw)
    assert len(rows) == len(cohort) == 105
    table(out/'bpp_screen.tsv', sorted(rows, key=lambda r:r['id']))
    build = ROOT/'build/linear-cparty-89baea52'
    lockpath = ROOT/'results/runs/pk-recovery-200-turner04-89baea52-20261018/software_lock.json'
    lock = json.loads(lockpath.read_text())
    parameter = ROOT/'deps/LinearCParty-89baea52/params/rna_Turner04.par'
    hashes = {p: sha((build/p).read_bytes()) for p in ('LinearPKCapR', 'LinearPKSampling')}
    assert lock['source_commit'] == REV and lock['energy_parameter_sha256'] == sha(parameter.read_bytes())
    assert all(lock['binaries'][k] == v for k, v in hashes.items())
    env = dict(os.environ, CPARTY_FIXED_ENERGY_PARAM_FILE=str(parameter),
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', LC_ALL='C')
    scratch_root = Path(os.environ['TMPDIR']).resolve()
    assert not str(scratch_root).startswith('/gs/bs/')
    start = time.monotonic()
    calls, profiles, samples = [], [], []
    with tempfile.TemporaryDirectory(prefix='example-search-', dir=scratch_root) as td:
        scratch = Path(td)
        with tarfile.open(out/'raw_small_calculations.tar.gz', 'w:gz') as archive:
            def add(name, raw):
                info = tarfile.TarInfo(name); info.size = len(raw)
                archive.addfile(info, io.BytesIO(raw))

            def call(name, program, args):
                elapsed = time.monotonic()-start
                if elapsed > 180:
                    raise RuntimeError('Small-computation budget exceeded; completed raw calls retained')
                cmd = [str(build/program), '-i', str(scratch/'input.txt'), '--exact', *args]
                st = time.monotonic()
                p = subprocess.run(cmd, capture_output=True, env=env, timeout=min(30, 185-elapsed))
                result = dict(id=name, program=program, command=cmd,
                              elapsed_seconds=time.monotonic()-st, returncode=p.returncode,
                              stdout_sha256=sha(p.stdout), stderr_sha256=sha(p.stderr))
                calls.append(result)
                add(name+'/'+program+'.stdout', p.stdout)
                add(name+'/'+program+'.stderr', p.stderr)
                assert p.returncode == 0, result
                return p.stdout.decode()

            for row in rows:
                if not row['crossing_mid_n']:
                    continue
                name = row['id']; d = data[name]; r = d['record']
                raw_input = f">{name}\n{r['seq']}\n{d['scaffold']}\n".encode()
                (scratch/'input.txt').write_bytes(raw_input)
                add(name+'/input.txt', raw_input)
                context = parse_context(call(name, 'LinearPKCapR', []), dict(name=name, n=row['n']))
                d['context'] = context
                assert np.max(np.abs(sum(context.values())-1)) < 1e-8
                for p in range(row['n']):
                    profiles.append(dict(id=name, position=p+1,
                                         **{c: float(context[c][p]) for c in CONTEXTS}))
                pseudo = context['PseudoLoop']
                fixedmask = np.zeros(row['n'], dtype=bool)
                for i,j in r['nested']:
                    fixedmask[[i-1,j-1]] = True
                mat = np.array([context[c] for c in CONTEXTS])
                entropy = -np.sum(mat*np.log2(np.maximum(mat, 1e-300)), axis=0)
                row.update(pseudoloop_max=float(max(pseudo)),
                           pseudoloop_max_position=int(np.argmax(pseudo))+1,
                           pseudoloop_mid_n=int(np.sum((pseudo>=.2)&(pseudo<=.8))),
                           context_mixed_positions=int(np.sum(np.max(mat, axis=0)<.95)),
                           free_positions=int(np.sum(~fixedmask)),
                           free_context_entropy_bits=float(np.mean(entropy[~fixedmask])))
                if not row['pseudoloop_mid_n']:
                    continue
                text = call(name, 'LinearPKSampling', ['--count','10000','--seed',str(SEED),'--unique'])
                counts = []
                for line in text.splitlines():
                    c = line.split()
                    if len(c) == 2 and c[0].isdigit() and len(c[1]) == row['n']:
                        counts.append((int(c[0]), c[1]))
                counts.sort(key=lambda x:(-x[0],x[1]))
                assert sum(c for c, _ in counts) == 10000
                G = pairs_of(d['scaffold'])
                assert all(G <= pairs_of(db) for _,db in counts)
                d['counts'] = counts
                for rank,(c,db) in enumerate(counts,1):
                    samples.append(dict(id=name, rank=rank, count=c, fraction=c/10000, structure=db))
                row.update(sample_n=10000, sample_unique=len(counts),
                           sample_top1=counts[0][0]/10000,
                           sample_top2=counts[1][0]/10000 if len(counts)>1 else 0.,
                           sample_top3=counts[2][0]/10000 if len(counts)>2 else 0.)
                # A meaningful sampling/BPP check includes every reported crossing pair.
                empirical = Counter()
                for c,db in counts:
                    for pair in pairs_of(db):
                        empirical[pair] += c
                row['sample_crossing_bpp_max_abs_error'] = max(
                    (abs(empirical[p]/10000-d['pairs'][p]) for p in d['crossing']), default=0.)
                print(name, row['pseudoloop_max'], row['pseudoloop_mid_n'],
                      row['sample_top1'], row['sample_top2'], flush=True)
            add('calls.json', json.dumps(calls, indent=2).encode())
        dump(out/'calls.json', calls)
    assert hashes == {p:sha((build/p).read_bytes()) for p in hashes}
    assert sha(parameter.read_bytes()) == lock['energy_parameter_sha256']
    # Unevaluated columns are explicit blanks, never silently treated as zeros.
    fields = list(dict.fromkeys(k for r in rows for k in r))
    table(out/'screen.tsv', [{k:r.get(k) for k in fields} for r in sorted(rows,key=lambda x:x['id'])])
    table(out/'position_profiles.tsv.gz', profiles)
    table(out/'sampling_structures.tsv.gz', samples)
    dump(out/'provenance.json', dict(revision=REV, seed=SEED,
         cohort_path=str(cohort_path.relative_to(ROOT)), cohort_sha256=sha(cohort_path.read_bytes()),
         source_artifacts=source, binary_sha256=hashes,
         parameter_sha256=sha(parameter.read_bytes()), new_cli_seconds=sum(c['elapsed_seconds'] for c in calls),
         new_cli_count=len(calls), script_sha256=sha(Path(__file__).read_bytes())))
    usable = [r for r in rows if r.get('sample_n')]
    dump(out/'candidate_data.json', {r['id']: dict(record=data[r['id']]['record'],
        scaffold=data[r['id']]['scaffold'],
        crossing_pairs=[dict(i=i,j=j,probability=data[r['id']]['pairs'][i,j])
                        for i,j in sorted(data[r['id']]['crossing'])],
        contexts={c:v.tolist() for c,v in data[r['id']]['context'].items()},
        top_samples=data[r['id']]['counts'][:20]) for r in usable})
    print('FIGURE COMPLETE', len(usable), 'eligible; elapsed', time.monotonic()-start, flush=True)


def neat_search(out):
    out.mkdir(exist_ok=False)
    provenance = dict(inputs=[], source_artifacts=[], revision=REV,
                      script_sha256=sha(Path(__file__).read_bytes()))
    inputs = read_inputs(provenance)
    saved = {}
    for base in ('e3-neat1-bpp', 'e3-neat1-capr', 'e3-neat1-variants-bpp', 'e3-neat1-variants-capr'):
        run = ROOT/'results/runs'/(base+'-89baea52-20261018')
        plan = json.loads((run/'array_plan.json').read_text())
        assert plan['software_before']['commit'] == REV
        for task in plan['tasks']:
            reader = RunReader(run/task['directory']/'run')
            assert reader.manifest['software_before']['commit'] == REV
            assert reader.manifest['software_lock_valid_after_run'] is True
            for job in reader.manifest['jobs']:
                assert job['status'] == 'success'
                if job['warmup']:
                    continue
                v = inputs[job['input']]
                entry = next(e for e in reader.manifest['inputs'] if e['name'] == job['input'])
                frozen = reader.input(entry)
                assert sha(frozen) == entry['sha256']
                assert frozen.decode().splitlines() == ['>'+v['name'],v['sequence'],v['scaffold']]
                cmd = job['command']; beam = int(cmd[cmd.index('-b')+1])
                for flag in ('--max-hairpin','--max-case3-depth'):
                    assert int(cmd[cmd.index(flag)+1]) == beam
                assert cmd[cmd.index('--ranking')+1] == 'outer'
                key = (v['iso'],v['method'],beam)
                g = saved.setdefault(key, dict(input=v))
                if job['program'] == 'LinearCParty':
                    g['cutoff'] = float(cmd[cmd.index('--bpp-cutoff')+1]); assert g['cutoff'] == .001
                    raw = provenance_artifact(reader, job, 'bpp.tsv', provenance['source_artifacts'])
                    g['pairs'] = bpp_parse(raw)
                    crossing = cross_vector(g['pairs'],v['fixed'])
                    g['crossing'] = {p:g['pairs'][p] for p in crossing}
                else:
                    raw = provenance_artifact(reader,job,'stdout.txt',provenance['source_artifacts'])
                    context = parse_context(raw.decode(),v)
                    assert np.max(np.abs(sum(context.values())-1)) < 1e-8
                    g['pseudo'] = context['PseudoLoop']
    assert len(saved) == 14
    # Audit the existing derived tables against freshly parsed and hashed raw artifacts.
    prior = ROOT/'results/local/e3-neat1-89baea52-20261018'
    derived_crossing = {}
    with gzip.open(prior/'bpp_pairs.tsv.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['method'] == 'LinearPartition':
                continue
            k = (r['iso'],r['method'],int(r['beam']))
            if r['crosses_G'] == '1':
                derived_crossing.setdefault(k,{})[int(r['i']),int(r['j'])] = float(r['probability'])
    assert all(g['crossing'] == derived_crossing.get(k,{}) for k,g in saved.items())
    full, crossings, positions = [], [], []
    for (iso,method,b),g in sorted(saved.items()):
        x = g['crossing']; pseudo=g['pseudo']; v=g['input']
        maximum = max(x.values()) if x else None
        full.append(dict(iso=iso,n=v['n'],G=method,G_pairs=len(v['fixed']),beam=b,
             crossing_max_exported=maximum, unexported_pair_upper_bound=g['cutoff'],
             crossing_exported_n=len(x),crossing_ge_005=sum(p>=.05 for p in x.values()),
             crossing_ge_01=sum(p>=.1 for p in x.values()),
             pseudoloop_max=float(max(pseudo)),pseudoloop_max_position=int(np.argmax(pseudo))+1))
        for (i,j),p in sorted(x.items()):
            crossings.append(dict(iso=iso,G=method,beam=b,i=i,j=j,span=j-i,probability=p))
        for i,p in enumerate(pseudo,1):
            positions.append(dict(iso=iso,G=method,beam=b,position=i,PseudoLoop=float(p)))
    table(out/'full_length_summary.tsv',full)
    table(out/'crossing_pairs.tsv',crossings)
    table(out/'pseudoloop_positions.tsv.gz',positions)
    # Exhaust all 300-nt windows (stride 1). Count pairs incident to the window
    # separately from pairs with both endpoints inside it; remote arcs matter here.
    all_windows, selected = [], []
    for iso in ('NEAT1_1','NEAT1_2'):
        n=inputs[iso]['n']; base=saved[iso,'G0.7',400]
        for a in range(1,n-300+2):
            z=a+299
            rec=dict(iso=iso,start=a,end=z)
            for b in (100,200,400):
                g=saved[iso,'G0.7',b]
                incident={p:q for p,q in g['crossing'].items() if a<=p[0]<=z or a<=p[1]<=z}
                internal={p:q for p,q in incident.items() if a<=p[0] and p[1]<=z}
                rec.update({f'b{b}_cross_max':max(incident.values()) if incident else None,
                    f'b{b}_cross_n':len(incident),
                    f'b{b}_cross_ge005':sum(q>=.05 for q in incident.values()),
                    f'b{b}_cross_ge01':sum(q>=.1 for q in incident.values()),
                    f'b{b}_cross_both_inside_n':len(internal),
                    f'b{b}_pseudo_max':float(max(g['pseudo'][a-1:z]))})
            all_windows.append(rec)
        # Rank every position by its incident crossing-pair maximum and then
        # PseudoLoop. Center the window on the anchor, clipping at sequence ends.
        # This avoids cutting off a peak to win arbitrary sliding-window ties.
        incident_peak = np.zeros(n)
        for (i,j),q in base['crossing'].items():
            incident_peak[i-1] = max(incident_peak[i-1],q)
            incident_peak[j-1] = max(incident_peak[j-1],q)
        anchors=sorted(range(1,n+1), key=lambda i:(-incident_peak[i-1],-base['pseudo'][i-1],i))
        windows={r['start']:r for r in all_windows if r['iso']==iso}
        chosen=[]
        for pos in anchors:
            a=max(1,min(n-299,int(pos)-149)); r=windows[a]
            if any(not(r['end']<s['start'] or s['end']<r['start']) for s in chosen):
                continue
            chosen.append(r)
            if len(chosen)==5:
                break
        selected.extend(chosen)
    selected.sort(key=lambda r:(-(r['b400_cross_max'] or 0),-r['b400_pseudo_max'],r['iso'],r['start']))
    for rank,r in enumerate(selected,1):
        r['rank']=rank
        r['pseudo_b400_minus_b200']=r['b400_pseudo_max']-r['b200_pseudo_max']
        r['pseudo_b400_over_b200']=r['b400_pseudo_max']/r['b200_pseudo_max'] if r['b200_pseudo_max'] else None
    # all_windows shares selected dict objects; write a uniform schema.
    allfields=list(dict.fromkeys(k for r in all_windows for k in r))
    table(out/'all_300nt_windows.tsv.gz',({k:r.get(k) for k in allfields} for r in all_windows))
    table(out/'ranked_windows.tsv',selected)
    dump(out/'provenance.json',provenance)
    dump(out/'validation.json',dict(status='passed',conditions=len(saved),
         raw_vs_prior_crossing_pairs_identical=True,
         source_member_hashes_checked=len(provenance['source_artifacts']),
         all_windows=len(all_windows),stride=1,width=300,
         crossing_definition='i<k<j<l or k<i<l<j for at least one (k,l) in G',
         absent_exported_pairs='unknown below 0.001; never report zero maximum'))
    print('NEAT COMPLETE',len(saved),'conditions;',len(all_windows),'windows',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('which',choices=['figure','neat'])
    p.add_argument('--out',type=Path)
    args=p.parse_args()
    default='figure1-candidates-20261020' if args.which=='figure' else 'neat1-pk-search-20261020'
    out=args.out or ROOT/'results/local'/default
    (figure_search if args.which=='figure' else neat_search)(out)


if __name__=='__main__':
    main()
