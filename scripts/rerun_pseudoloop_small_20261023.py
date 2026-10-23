#!/usr/bin/env python3
"""Recompute the saved small context examples with an archived 0.2.0 build.

The source/build arguments must be node-local. Raw files are packed, and all
scientific comparisons use the saved manuscript inputs rather than new folds.
"""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tarfile
import time

import make_figure4_terc as terc
import validate_byloop_89baea52 as byloop

ROOT = Path(__file__).resolve().parents[1]
CONTEXTS = ['Stem', 'Hairpin', 'Bulge', 'Internal', 'Multiloop', 'Exterior', 'PseudoLoop']


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def dump(p, d):
    Path(p).write_text(json.dumps(d, indent=2, ensure_ascii=False, allow_nan=False)+'\n')


def access(p):
    rows = {}
    for line in Path(p).read_text().splitlines():
        f = line.split('\t')
        if len(f) < 2 or not f[0].isdigit():
            continue
        start = int(f[0])+1
        for entry in f[1].strip(';').split(';'):
            w, v = entry.split(','); rows[start, int(w)] = {'total': float(v)}
        if len(f) > 2 and f[2].startswith('byloop:'):
            for entry in f[2][7:].strip(';').split(';'):
                w, *parts = entry.split(',')
                rows[start, int(w)].update({k: float(v) for k, v in (x.split('=') for x in parts)})
    return rows


def invariants(old, new, channels, multi, pseudo):
    other = max(abs(old[c][i]-new[c][i]) for c in channels if c not in (multi, pseudo)
                for i in range(len(old[c])))
    mass = max(abs(old[multi][i]+old[pseudo][i]-new[multi][i]-new[pseudo][i])
               for i in range(len(old[multi])))
    gain = [new[multi][i]-old[multi][i] for i in range(len(old[multi]))]
    loss = [old[pseudo][i]-new[pseudo][i] for i in range(len(old[multi]))]
    result = dict(other_channels_max_abs_difference=other, multi_plus_pseudo_max_abs_difference=mass,
                  min_multiloop_gain=min(gain), min_pseudoloop_loss=min(loss),
                  max_transferred_probability=max(gain))
    assert other < 1e-9 and mass < 1e-9 and min(gain) > -1e-9 and min(loss) > -1e-9, result
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source', type=Path, required=True)
    ap.add_argument('--build', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--section', choices=['examples', 'validation', 'enumeration'], required=True)
    args = ap.parse_args()
    source, build, out = args.source.resolve(), args.build.resolve(), args.out.resolve()
    assert not str(source).startswith('/gs/bs/') and not str(build).startswith('/gs/bs/')
    out.mkdir(parents=True, exist_ok=False)
    scratch = source.parent / ('small-'+args.section)
    scratch.mkdir(exist_ok=False)
    env = {k:v for k,v in os.environ.items() if not k.startswith('CPARTY_')}
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', LC_ALL='C',
               PYTHONDONTWRITEBYTECODE='1', CPARTY_PF_OVERRIDE_SCALE='1',
               CPARTY_FIXED_ENERGY_PARAM_FILE=str(source/'params/rna_Turner04.par'))
    provenance = dict(source_commit=subprocess.check_output(['git','-C',str(ROOT/'../../LinearKnot'),
                      'rev-parse','f93b55e'],text=True).strip(), version=(source/'VERSION').read_text().strip(),
                      source_acquisition='git archive f93b55e | tar -x; node-local; no worktree',
                      binary_sha256=sha(build/'linearknot'), parameter_sha256=sha(source/'params/rna_Turner04.par'),
                      cmake_cache_sha256=sha(build/'CMakeCache.txt'), script_sha256=sha(__file__), calls=[])
    archive = tarfile.open(out/'raw_outputs.tar.gz','w:gz')

    def run(argv, label, cwd=source, timeout=600, stdin=None):
        t = time.monotonic()
        p = subprocess.run(list(map(str,argv)),input=stdin,capture_output=True,cwd=cwd,env=env,timeout=timeout)
        so, se = scratch/(label+'.stdout'), scratch/(label+'.stderr')
        so.write_bytes(p.stdout); se.write_bytes(p.stderr)
        archive.add(so,arcname=label+'.stdout'); archive.add(se,arcname=label+'.stderr')
        provenance['calls'].append(dict(argv=list(map(str,argv)),cwd=str(cwd),returncode=p.returncode,
            elapsed_seconds=time.monotonic()-t,stdout_sha256=sha(so),stderr_sha256=sha(se)))
        dump(out/'provenance.json',provenance)
        assert p.returncode == 0, (label,p.stderr.decode()[-3000:])
        print(label+' complete',flush=True)
        return so

    if args.section == 'examples':
        saved = ROOT/'results/local/figure-feedback-20261020/data.json'
        d = json.loads(saved.read_text()); provenance['saved_data_sha256'] = sha(saved)
        report = {}
        for group in ['figure1','terc']:
            folder = scratch/group; folder.mkdir()
            if group == 'figure1':
                item = d[group]
                (folder/'in.txt').write_text('>'+item['record']['id']+'\n'+item['record']['seq']+'\n'+item['scaffold']+'\n')
            else:
                oldtar = ROOT/'results/local/small-rerun-89baea52-20261018/terc_data.tar'
                provenance['terc_saved_archive_sha256'] = sha(oldtar)
                with tarfile.open(oldtar) as tf:
                    for m in tf:
                        assert m.isfile() and Path(m.name).name == m.name
                        (folder/m.name).write_bytes(tf.extractfile(m).read())
            provenance[group+'_input_sha256'] = sha(folder/'in.txt')
            n = len((folder/'in.txt').read_text().splitlines()[1])
            settings = [('exact',['--exact'])] if group == 'figure1' else [
                ('exact',['--exact']),('b200',['-b','200','--max-hairpin','200','--max-case3-depth','200','--ranking','outer'])]
            for mode,opts in settings:
                prefix = group+'-'+mode
                oldctx = d['figure1']['context'] if group == 'figure1' else d['terc']['settings'][mode]['context']
                prof = run([build/'linearknot','profile','-i',folder/'in.txt',*opts],prefix+'-profile')
                ctx = {c: v.tolist() for c,v in terc.read_context(prof,n).items()}
                acc = run([build/'linearknot','access','-i',folder/'in.txt',*opts,'--access-len','1,4,8','--probability','--by-loop'],prefix+'-access')
                newacc = access(acc)
                inv = invariants(oldctx,ctx,CONTEXTS,'Multiloop','PseudoLoop')
                inv['context_sum_error'] = max(abs(sum(ctx[c][i] for c in CONTEXTS)-1) for i in range(n))
                inv['byloop_sum_error'] = max(abs(sum(v[c] for c in 'HBIMEP')-v['total']) for v in newacc.values())
                if group == 'figure1':
                    oldacc = {(int(i),int(w)): {'total':p} for w,rows in d[group]['access'].items() for i,p in rows.items()}
                    d[group]['context'] = ctx
                    inv.update(old_pseudoloop_max=max(oldctx['PseudoLoop']),new_pseudoloop_max=max(ctx['PseudoLoop']),
                        old_pseudoloop_between_02_08=sum(.2<=x<=.8 for x in oldctx['PseudoLoop']),
                        new_pseudoloop_between_02_08=sum(.2<=x<=.8 for x in ctx['PseudoLoop']))
                else:
                    oldacc = access(folder/(mode+'.racc'))
                    oldch = {c:[v[c] for _,v in sorted(oldacc.items())] for c in 'HBIMEP'}
                    newch = {c:[newacc[k][c] for k in sorted(oldacc)] for c in 'HBIMEP'}
                    inv['byloop_old_new'] = invariants(oldch,newch,'HBIMEP','M','P')
                    inv['window_123_126'] = dict(old=oldacc[123,4],new=newacc[123,4])
                    inv['window_123_126']['new_pseudoloop_fraction'] = newacc[123,4]['P']/newacc[123,4]['total']
                    d['terc']['settings'][mode]['context'] = ctx
                    (folder/(mode+'.capr')).write_bytes(prof.read_bytes())
                    (folder/(mode+'.racc')).write_bytes(acc.read_bytes())
                inv['total_access_max_abs_difference'] = max(abs(v['total']-newacc[k]['total']) for k,v in oldacc.items())
                assert inv['total_access_max_abs_difference'] < 1e-9 and inv['byloop_sum_error'] < 1e-9
                report[prefix] = inv
                archive.add(folder/'in.txt',arcname=prefix+'/in.txt')
            if group == 'terc':
                with tarfile.open(out/'terc_data.tar','w') as tf:
                    for p in folder.iterdir(): tf.add(p,arcname=p.name)
        bctx, ectx = [d['terc']['settings'][m]['context'] for m in ['b200','exact']]
        report['terc_displayed_40_160_context_max_abs_error'] = max(abs(bctx[c][i]-ectx[c][i]) for c in CONTEXTS for i in range(39,160))
        dump(out/'data.json',d); dump(out/'summary.json',report)
    elif args.section == 'validation':
        run([build/'scfg2_pk_profiles_oracle_test'],'maintained_oracle')
        run([build/'scfg2_pk_profiles_beam_test'],'maintained_beam')
        geometry = run(['python3','-B',source/'test/tools/test_pseudoloop_geometry.py',build],'geometry')
        native = (source/'test/scfg2_pk_profiles_oracle_test.cc').read_text()
        # Retain exactly the original SI table's 21 strict-policy conditions.
        old = '  return check_case_policy(label, sequence, scaffold, WmbClosurePolicy::LegacyCompatible) && ok;'
        assert native.count(old) == 1
        native = native.replace(old,'  return ok;')
        probe = scratch/'byloop.cc'; probe.write_text(byloop.instrument(native))
        archive.add(probe,arcname='byloop.cc')
        run(['c++','-O2','-std=gnu++17','-I'+str(source/'src'),'-I'+str(source/'test'),
             '-I'+str(build/'generated'),probe,build/'libcparty_scfg2_exact.a',build/'libRNA.a','-o',scratch/'byloop'],'compile_byloop')
        provenance['byloop_binary_sha256'] = sha(scratch/'byloop')
        output = run([scratch/'byloop'],'byloop_metrics')
        rows = [json.loads(s[7:]) for s in output.read_text().splitlines() if s.startswith('METRIC ')]
        assert len(rows) == 294
        groups = {'exact_vs_enumeration':[r for r in rows if r['mode']=='exact_vs_enumeration'],
                  'unpruned_vs_exact':[r for r in rows if r['mode']=='unpruned_vs_exact'],
                  'finite_vs_taped':[r for r in rows if r['mode'].startswith('finite_')]}
        summary = dict(contexts=['Hairpin','Bulge','Internal','Multi','Exterior','PseudoLoop'],
            geometry=json.loads(geometry.read_text()),comparisons={})
        for name,rs in groups.items():
            summary['comparisons'][name] = dict(conditions=len(rs),
                component_max_abs_error=[max(r['component_max_abs_error'][c] for r in rs) for c in range(6)],
                sum_max_abs_error=max(r['sum_max_abs_error'] for r in rs),unit_max_abs_error=max(r['unit_max_abs_error'] for r in rs))
        summary.update(distinct_inputs=len({(r['sequence'],r['scaffold']) for r in groups['exact_vs_enumeration']}),
            enumerated_groups=sum(r['groups'] for r in groups['exact_vs_enumeration']),
            positive_longer_windows=[sum(r['positive_longer_windows'][c] for r in groups['exact_vs_enumeration']) for c in range(6)])
        dump(out/'summary.json',summary)
        (out/'case_metrics.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
    else:
        # Independent table transcription: only seven wing rule IDs change.
        oldsource = ROOT/'tools/validation_89baea52/enumeration/enumerate.cc'
        code = oldsource.read_text()
        changes = {
            'case RI::WI_BASE_SINGLE:case RI::WI_EXTEND_UNPAIRED:case RI::WIP_EXTEND_UNPAIRED:\n      emit(j,j,\'P\');break;':
            'case RI::WI_BASE_SINGLE:case RI::WI_EXTEND_UNPAIRED:\n      emit(j,j,\'P\');break;\n    case RI::WIP_EXTEND_UNPAIRED:emit(j,j,\'M\');break;',
            "emit(i,k-1,'P');break;":"emit(i,k-1,'M');break;",
            "case RI::VPR_SPLIT_VP_BASEPAIR:emit(k+1,j,'P');break;":"case RI::VPR_SPLIT_VP_BASEPAIR:emit(k+1,j,'M');break;",
            "case RI::BE_WIP_BASEPAIR:emit(l+1,j-1,'P');break;":"case RI::BE_WIP_BASEPAIR:emit(l+1,j-1,'M');break;",
            "case RI::BE_BASEPAIR_WIP:emit(i+1,k-1,'P');break;":"case RI::BE_BASEPAIR_WIP:emit(i+1,k-1,'M');break;"}
        for a,b in changes.items():
            assert code.count(a)==1,a
            code=code.replace(a,b)
        # The archived observer predates removal of a redundant option field.
        assert code.count('weights::RuntimeOptions{2,false}') == 1
        code = code.replace('weights::RuntimeOptions{2,false}', 'weights::RuntimeOptions{2}')
        probe = scratch/'enumerate.cc'; probe.write_text(code); archive.add(probe,arcname='enumerate.cc')
        run(['c++','-O1','-std=c++17','-I'+str(source/'src'),probe,build/'libcparty_scfg2_beam.a',
             build/'libcparty_scfg2_exact.a',build/'libRNA.a','-o',scratch/'enumerate'],'compile_enumeration')
        provenance['enumeration_binary_sha256'] = sha(scratch/'enumerate')
        casespath = ROOT/'tools/validation_89baea52/baseline_15aa7742/inputs/cases.json'
        cases = json.loads(casespath.read_text()); provenance['cases_sha256']=sha(casespath)
        from importlib.util import spec_from_file_location, module_from_spec
        spec=spec_from_file_location('old_enum',ROOT/'tools/validation_89baea52/enumeration/run.py')
        mod=module_from_spec(spec); spec.loader.exec_module(mod)
        rows=[]
        for case in cases:
            output=run([scratch/'enumerate',case['sequence'],case['scaffold'],'default'],case['id'])
            values,structures,_=mod.parse(output.read_text())
            assert all(values[k]==0 for k in mod.ZEROS) and values['max_count']==1
            assert all(values[k]<1e-9 for k in mod.ERRORS)
            rows.append(dict(case,**values))
        summary=dict(cases=len(rows),unique_inputs=len({(r['sequence'],r['scaffold']) for r in rows}),
            totals={k:sum(r[k] for r in rows) for k in mod.TOTALS},
            maxima={k:max(r[k] for r in rows) for k in mod.ERRORS},failures=[])
        dump(out/'summary.json',summary); dump(out/'report.json',dict(enumeration=rows))
    archive.close()
    provenance['raw_archive_sha256']=sha(out/'raw_outputs.tar.gz')
    dump(out/'provenance.json',provenance)


if __name__ == '__main__':
    main()
