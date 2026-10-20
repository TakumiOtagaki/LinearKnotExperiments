#!/usr/bin/env python3
"""Small local calculations for the 2026-10-19 figure feedback.

Per-process files live in node-local scratch; retain one tar and JSON summaries.
Existing Terc probabilities and samples remain those measured with 89baea52.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile
import time

import make_figure1 as f1
import make_figure4_terc as terc
from analyze_pk_membership import parse_structure, CANONICAL

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / 'results/local/figure-feedback-20261020'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def decode(path, n, g):
    lines = [s for s in path.read_text().splitlines() if 'EG=' in s]
    assert len(lines) == 1 and 'gamma=1;' in lines[0]
    db = lines[0].split()[0]
    assert len(db) == n and g <= terc.pairs_of(db)
    return db


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', type=Path, default=DEFAULT)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    assert not list(args.out.iterdir()), 'Choose an empty output directory'
    binary = ROOT / 'build/linearknot-9c7d6e1/linearknot'
    parameter = ROOT / 'deps/LinearKnot-9c7d6e1/params/rna_Turner04.par'
    lock = ROOT / 'manifests/software-linearknot-9c7d6e1.json'
    selected = ROOT / 'results/local/figure1-candidates-20261020/selected_inputs.json'
    source_tar = ROOT / 'results/local/small-rerun-89baea52-20261018/terc_data.tar'
    sources = {str(p.relative_to(ROOT)): sha(p) for p in (binary, parameter, lock, selected, source_tar)}
    scratch_root = Path(os.environ.get('TMPDIR', '/tmp')).resolve()
    assert not str(scratch_root).startswith('/gs/bs/')
    env = dict(os.environ, CPARTY_FIXED_ENERGY_PARAM_FILE=str(parameter),
               OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', LC_ALL='C')
    calls = []
    with tempfile.TemporaryDirectory(prefix='figure-feedback-', dir=scratch_root) as td:
        scratch = Path(td)
        figdir, tercdir = scratch / 'figure1', scratch / 'terc'
        figdir.mkdir(); tercdir.mkdir()
        def run(subcommand, folder, output, *options):
            argv = [str(binary), subcommand, '-i', str(folder / 'in.txt'), *options]
            start = time.monotonic()
            result = subprocess.run(argv, env=env, cwd=folder, capture_output=True, timeout=1800)
            (folder / output).write_bytes(result.stdout)
            (folder / (output + '.stderr')).write_bytes(result.stderr)
            calls.append(dict(command=argv, returncode=result.returncode,
                              elapsed_seconds=time.monotonic() - start,
                              stdout_sha256=sha(folder / output),
                              stderr_sha256=sha(folder / (output + '.stderr'))))
            assert result.returncode == 0, result.stderr.decode()
            print(f'{folder.name}: {subcommand} complete', flush=True)
        item = json.loads(selected.read_text())['bpRNA_RFAM_31397']
        rec, scaffold = item['record'], item['scaffold']
        seq, n = rec['seq'], len(rec['seq'])
        g = terc.pairs_of(scaffold)
        expected = {(i,j) for i,j,o in parse_structure(rec['ref'])
                    if o == '(' and j-i >= 4 and seq[i-1]+seq[j-1] in CANONICAL}
        assert g == expected == set(map(tuple, rec['nested'])) and n == 61 and len(g) == 8
        dump(figdir / 'rec.json', rec)
        (figdir / 'in.txt').write_text(f">{rec['id']}\n{seq}\n{scaffold}\n")
        run('fold', figdir, 'lc.out', '--exact', '--bpp', str(figdir/'bpp.tsv'), '--bpp-cutoff', '0')
        run('profile', figdir, 'capr.out', '--exact')
        run('access', figdir, 'racc.out', '--exact', '--access-len', '1,4,8', '--probability')
        run('sample', figdir, 'samp.out', '--exact', '-k', '10000', '--seed', '20261020', '--unique')
        run('centroid', figdir, 'centroid.out', '--exact', '--gamma', '1')
        ctx = {c: v.tolist() for c,v in terc.read_context(figdir/'capr.out', n).items()}
        bpp = terc.read_bpp(figdir/'bpp.tsv', n)
        acc = {}
        for line in (figdir/'racc.out').read_text().splitlines():
            if not line.strip() or line.startswith(('#','>')): continue
            pos, rest = line.split('\t')
            for entry in rest.strip(';').split(';'):
                w,p = entry.split(',')
                acc.setdefault(w, {})[str(int(pos)+1)] = float(p)
        samples = sorted([(int(s.split('\t')[0]), s.split('\t')[1])
                          for s in (figdir/'samp.out').read_text().splitlines()
                          if '\t' in s and s.split('\t')[0].isdigit()], reverse=True)
        assert sum(c for c,_ in samples) == 10000
        assert all(g <= terc.pairs_of(db) and len(db) == n for _,db in samples)
        assert all(abs(bpp[p]-1) < 1e-10 for p in g)
        sum_error = max(abs(sum(ctx[c][i] for c in f1.CONTEXTS)-1) for i in range(n))
        assert sum_error < 1e-10
        assert all(len(v) == n-int(w)+1 for w,v in acc.items())
        sample_error = max(abs(sum(c for c,db in samples if pair in terc.pairs_of(db))/10000-p)
                           for pair,p in bpp.items())
        figure = dict(record=rec, scaffold=scaffold, context=ctx,
                      bpp=[[i,j,p] for (i,j),p in sorted(bpp.items())], access=acc,
                      samples=samples, centroid=decode(figdir/'centroid.out', n, g),
                      validation=dict(max_context_sum_error=sum_error,
                                      max_sample_bpp_difference=sample_error))
        with tarfile.open(source_tar) as t:
            for member in t:
                assert member.isfile() and Path(member.name).name == member.name
                (tercdir/member.name).write_bytes(t.extractfile(member).read())
        trec = json.loads((tercdir/'rec.json').read_text())
        tscaffold = (tercdir/'in.txt').read_text().splitlines()[-1]
        tg = terc.pairs_of(tscaffold)
        centroid = {}
        for mode, opts in (('b200', ['-b','200','--max-hairpin','200','--max-case3-depth','200']),
                           ('exact', ['--exact'])):
            run('centroid', tercdir, mode+'.centroid', *opts, '--gamma', '1')
            centroid[mode] = decode(tercdir/(mode+'.centroid'), 388, tg)
        pk = set(map(tuple, trec['pk']))
        stems = dict(outer={p for p in pk if p[0] <= 68}, inner={p for p in pk if p[0] >= 69}, both=pk)
        tdata = {}
        for mode in ('b200', 'exact'):
            setting = terc.read_setting(tercdir, mode, 388, tscaffold)
            stats = terc.summarize(setting, stems, 388)
            for name,pairs in stems.items():
                stats['stems'][name].pop('mea_pairs_present')
                stats['stems'][name]['centroid_pairs_present'] = sorted(pairs & terc.pairs_of(centroid[mode]))
            tdata[mode] = dict(context={c:v.tolist() for c,v in setting['context'].items()},
                              bpp=[[i,j,p] for (i,j),p in sorted(setting['bpp'].items())],
                              access=[[i,w,v] for (i,w),v in sorted(setting['access'].items())],
                              centroid=centroid[mode], stats=stats, validation=setting['validation'])
        differences = terc.compare(*(terc.read_setting(tercdir,m,388,tscaffold) for m in ('b200','exact')))
        dump(args.out/'data.json', dict(figure1=figure, terc=dict(record=trec, scaffold=tscaffold,
                                                               settings=tdata, differences=differences)))
        with tarfile.open(args.out/'raw_outputs.tar', 'w') as t:
            for p in sorted(scratch.rglob('*')):
                if p.is_file(): t.add(p, arcname=str(p.relative_to(scratch)))
    assert sources == {p: sha(ROOT/p) for p in sources}
    dump(args.out/'calls.json', calls)
    dump(args.out/'provenance.json', dict(sources=sources, calls=7, rna_computation='local; no timing benchmark',
        new_revision='9c7d6e170b2b21a94a2e626d9939817e62898253', terc_saved_revision='89baea52',
        source_files_unchanged=True, raw_archive_sha256=sha(args.out/'raw_outputs.tar'),
        data_sha256=sha(args.out/'data.json'), script_sha256=sha(Path(__file__))))


if __name__ == '__main__':
    main()
