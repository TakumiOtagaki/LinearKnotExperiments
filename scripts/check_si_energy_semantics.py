#!/usr/bin/env python3
"""Small structural accounting audit; compile and run only in node-local scratch.

The SI columns change one disputed printed term at a time, keeping BE base=1.
They do not pretend that the literal SI (whose BE base is zero) has finite PK weights.
"""
import argparse
import difflib
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def example(name, n, g, gp, leaves=(), ml=None, pk=None, band=None):
    s = ['A'] * n
    for i, j in [*g, *gp]:
        assert s[i-1] == s[j-1] == 'A'
        s[i-1], s[j-1] = 'G', 'C'
    def dots(pairs):
        text = ['.'] * n
        for i, j in pairs:
            text[i-1], text[j-1] = '(', ')'
        return ''.join(text)
    return dict(name=name, n=n, S=''.join(s), G=dots(g), Gprime=dots(gp),
                G_pairs=g, Gprime_pairs=gp, leaves=leaves, ml=ml, pk=pk, band=band)


def examples():
    return [
        example('ordinary_ml', 24, [(2,23),(5,9),(14,18)], [], [(5,9),(14,18)],
                ml=dict(closing=(2,23), components=[(5,9),(14,18)], pk_branches=0)),
        example('be_single_pair', 20, [(1,12)], [(5,20)],
                pk=dict(span=(1,20), ext=(5,20), outer=(1,12), inner=(1,12), nested=False)),
        example('wi_empty_wing', 7, [(1,6)], [(2,7)],
                pk=dict(span=(1,7), ext=(2,7), outer=(1,6), inner=(1,6), nested=False)),
        example('ml_single_pk', 22, [(1,22),(4,13)], [(8,18)],
                ml=dict(closing=(1,22), components=[(4,18)], pk_branches=1),
                pk=dict(span=(4,18), ext=(8,18), outer=(4,13), inner=(4,13), nested=True)),
        example('ml_terminal_pk', 32, [(1,32),(3,7),(11,20)], [(15,25)], [(3,7)],
                ml=dict(closing=(1,32), components=[(3,7),(11,25)], pk_branches=1),
                pk=dict(span=(11,25), ext=(15,25), outer=(11,20), inner=(11,20), nested=True)),
        example('ml_nonterminal_pk', 32, [(1,32),(4,13),(24,28)], [(8,18)], [(24,28)],
                ml=dict(closing=(1,32), components=[(4,18),(24,28)], pk_branches=1),
                pk=dict(span=(4,18), ext=(8,18), outer=(4,13), inner=(4,13), nested=True)),
        example('be_both_wings', 38, [(1,33),(10,24),(3,8),(26,31)], [(13,38)], [(3,8),(26,31)],
                pk=dict(span=(1,38), ext=(13,38), outer=(1,33), inner=(10,24), nested=False),
                band=dict(branches=[(3,8),(26,31)], disputed='both')),
        example('be_right_unpaired_zero', 29, [(1,24),(10,23),(3,8)], [(13,29)], [(3,8)],
                pk=dict(span=(1,29), ext=(13,29), outer=(1,24), inner=(10,23), nested=False),
                band=dict(branches=[(3,8)], disputed='right_unpaired')),
        example('be_right_unpaired_three', 32, [(1,27),(10,23),(3,8)], [(13,32)], [(3,8)],
                pk=dict(span=(1,32), ext=(13,32), outer=(1,27), inner=(10,23), nested=False),
                band=dict(branches=[(3,8)], disputed='right_unpaired')),
        example('be_left_unpaired_three', 32, [(1,27),(5,18),(20,25)], [(8,32)], [(20,25)],
                pk=dict(span=(1,32), ext=(8,32), outer=(1,27), inner=(5,18), nested=False),
                band=dict(branches=[(20,25)], disputed='left_unpaired')),
        example('wip_shortest_branch_root', 29, [(1,24),(10,23),(3,7)], [(13,29)], [(3,7)],
                pk=dict(span=(1,29), ext=(13,29), outer=(1,24), inner=(10,23), nested=False),
                band=dict(branches=[(3,7)], disputed='missing_wip')),
        example('four_band_chain', 30, [(4,17),(20,30)], [(2,9),(12,25)]),
        example('vp_right_lower_bound', 21, [(1,15),(5,10)], [(2,21),(11,16)], [(5,10)]),
    ]


def interval(i, j):
    return set(range(i, j+1))


def hand_energy(e, hairpins, affine):
    parts = dict(hairpins=sum(hairpins[tuple(p)] for p in e['leaves']))
    counts = {}
    if e['name']=='four_band_chain':
        free = interval(2,30)-{p for pair in e['G_pairs']+e['Gprime_pairs'] for p in pair}
        counts.update(pk_unpaired=sorted(free), bands=4)
        parts['bands_and_pseudoloop'] = 4*2.46+.06*len(free)-1.38
    if e['name']=='vp_right_lower_bound':
        free = (interval(3,10) | interval(17,20))-interval(5,10)
        counts.update(band_unpaired=sorted(free),band_degree=3,pk_unpaired=[12,13,14])
        parts['band_multiloop'] = 3.41+3*.56+.12*len(free)
        parts['bands_and_pseudoloop'] = 2*2.46+3*.06-1.38
    if e['ml']:
        ml = e['ml']
        i, j = ml['closing']
        free = interval(i+1,j-1)
        for a,b in ml['components']:
            free -= interval(a,b)
        counts.update(ml_unpaired=sorted(free), ml_degree=1+len(ml['components']))
        # Existing Turner04 code uses CParty b for a PK branch, Vienna b elsewhere.
        a,b,c = affine
        parts['ordinary_ml'] = a+b*(counts['ml_degree']-ml['pk_branches']) + .03*ml['pk_branches'] + c*len(free)
        parts['PSM'] = 10.07*ml['pk_branches']
    if e['pk']:
        pk = e['pk']
        i,j = pk['outer']; l,lp = pk['inner']
        free = interval(*pk['span']) - interval(i,l) - interval(lp,j) - set(pk['ext'])
        counts['pk_unpaired'] = sorted(free)
        parts['bands_and_pseudoloop'] = 2*2.46+.06*len(free)+(0 if pk['nested'] else -1.38)
        if e['band']:
            free_band = interval(i+1,l-1) | interval(lp+1,j-1)
            for a,b in e['band']['branches']:
                free_band -= interval(a,b)
            counts.update(band_unpaired=sorted(free_band), band_degree=2+len(e['band']['branches']),
                          right_unpaired=j-lp-1, left_span=l-i-1)
            parts['band_multiloop'] = 3.41+.56*counts['band_degree']+.12*len(free_band)
    return parts, counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--software-root', type=Path, default=Path('experiments/deps/LinearCParty-15aa7742'))
    parser.add_argument('--build-dir', type=Path, default=Path('experiments/build/linear-cparty-15aa7742'))
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--oracle-root', type=Path, default=Path('experiments/build/cparty-oracle-800e95a-4patches/src'))
    args = parser.parse_args()
    src, build = args.software_root.resolve(), args.build_dir.resolve()
    probe = Path(__file__).with_name('si_energy_semantics_probe.cc').resolve()
    legacy_probe = probe.with_name('si_energy_legacy_probe.cc')
    oracle = args.oracle_root.resolve()
    legacy_objects = [oracle/'build/CMakeFiles/CParty.dir/src'/f'{name}.cc.o'
                      for name in ['part_func','sparse_tree','SHAPE','dot_plot','mea','centroid']]
    tracked = [probe, Path(__file__).resolve(), build/'libcparty_scfg2_exact.a', build/'libRNA.a',
               src/'src/scfg2/replay/w_final_exact_inside.hh', src/'src/scfg2/deductions/generate_exact_basic.hh',
               src/'params/rna_Turner04.par', src/'params/rna_DirksPierce09.par']
    tracked += [legacy_probe,oracle/'src/part_func.cc',oracle/'src/part_func.hh',
                oracle/'src/pseudo_loop.cc',oracle/'build/libRNA.a',*legacy_objects]
    before = {str(p):digest(p) for p in tracked}
    records, raw, legacy = [], [], []
    env = os.environ.copy()
    env['CPARTY_PF_OVERRIDE_SCALE'] = '1'
    with tempfile.TemporaryDirectory(prefix='si-energy-', dir='/tmp') as scratch:
        binary = Path(scratch)/'probe'
        command = ['g++','-O0','-std=gnu++17','-I'+str(src/'src'),str(probe),'-o',str(binary),
                   str(build/'libcparty_scfg2_exact.a'),str(build/'libRNA.a')]
        subprocess.run(command, check=True, capture_output=True, text=True)
        for model in ['rna_Turner04.par','rna_DirksPierce09.par']:
            for e in examples():
                command = [str(binary),e['S'],e['G'],e['Gprime'],str(src/'params'/model),'0']
                p = subprocess.run(command, env=env, text=True, capture_output=True, check=True, timeout=30)
                raw.append(dict(case=e['name'],model=model,command=command,stdout=p.stdout,stderr=p.stderr))
                lines = [x.split('\t') for x in p.stdout.splitlines()]
                root = next(x for x in lines if x[0]=='ROOT')
                param = next(x for x in lines if x[0]=='PARAM')
                hairpins = {(int(x[1]),int(x[2])):float(x[3]) for x in lines if x[0]=='HAIRPIN'}
                affine = tuple(float(x) for x in param[2:5])
                parts, counts = hand_energy(e, hairpins, (3.39,.03,.02))
                code_parts, _ = hand_energy(e, hairpins, affine)
                hand, code_hand = sum(parts.values()), sum(code_parts.values())
                reached = root[1:3] == ['1','1']
                observed = float(root[3]) if reached else None
                if e['name']=='wip_shortest_branch_root':
                    assert not reached, (e['name'],root)
                else:
                    assert reached and math.isclose(observed,code_hand,abs_tol=1e-10), (e['name'],model,observed,code_hand)
                api = next(x for x in lines if x[0]=='API_D2')
                rules = [dict(rule=x[1],i=int(x[2]),j=int(x[3]),ip=int(x[4]),jp=int(x[5]),
                              k=int(x[6]),l=int(x[7]),energy=float(x[8])) for x in lines if x[0]=='RULE']
                isolated = hand
                if e['name'] in ['ml_single_pk','ml_terminal_pk']: isolated += 10.07
                if e['band'] and e['band']['disputed']=='both': isolated += .56
                if e['band'] and e['band']['disputed']=='right_unpaired': isolated += .24
                if e['name']=='wi_empty_wing': isolated = None  # Literal WI(empty)=0 kills this product.
                records.append(dict(**e,model=model,dangles=0,affine=affine,counts=counts,manual_parts=parts,
                                    manual_energy=hand,code_manual_energy=code_hand,observed_energy=observed,
                                    api_dangles2_energy=float(api[3]) if api[1:3]==['1','1'] else None,
                                    si_isolated_energy=isolated,rules=rules))
        original = (oracle/'src/part_func.cc').read_text()
        begin = original.index('void W_final_pf::compute_WIP(')
        end = original.index('void W_final_pf::compute_VPL(',begin)
        block = original[begin:end]
        assert block.count('k < j - TURN - 1') == 1
        patched = original[:begin]+block.replace('k < j - TURN - 1','k <= j - TURN - 1')+original[end:]
        patch_text = ''.join(difflib.unified_diff(original.splitlines(True),patched.splitlines(True),
                                                fromfile='a/src/part_func.cc',tofile='b/src/part_func.cc'))
        changed_source = Path(scratch)/'part_func.cc'
        changed_source.write_text(patched)
        changed_object = Path(scratch)/'part_func.cc.o'
        subprocess.run(['g++','-O0','-std=gnu++17','-I'+str(oracle/'src'),'-c',str(changed_source),
                        '-o',str(changed_object)],check=True,capture_output=True,text=True)
        for variant, objects in [('original',legacy_objects),('boundary_only', [changed_object,*legacy_objects[1:]])]:
            binary = Path(scratch)/('legacy_'+variant)
            subprocess.run(['g++','-O0','-std=gnu++17','-I'+str(oracle/'src'),str(legacy_probe),'-o',str(binary),
                            *map(str,objects),str(oracle/'build/libRNA.a')],check=True,capture_output=True,text=True)
            for model, h3, h4 in [('rna_Turner04.par',5.4,4.5),('rna_DirksPierce09.par',3.69,2.88)]:
                for name,s,i,j,expected in [('missing_prefix','AAGAAAC',2,7,h3+.56+.12),
                                           ('base_shortest','AAGAAAC',3,7,h3+.56),
                                           ('prefix_four','AAGAAAAC',2,8,h4+.56+.12)]:
                    command = [str(binary),s,'.'*len(s),str(src/'params'/model),str(i),str(j)]
                    p = subprocess.run(command,env=env,text=True,capture_output=True,check=True,timeout=30)
                    row = json.loads(p.stdout)
                    assert row['pf_scale']==row['WI_empty']==1
                    expect_zero = variant=='original' and name=='missing_prefix'
                    if expect_zero: assert row['WIP']==0
                    else: assert math.isclose(-row['rt']*math.log(row['WIP']),expected,abs_tol=1e-10)
                    legacy.append(dict(case=name,variant=variant,model=model,S=s,G='.'*len(s),i=i,j=j,
                                       expected_energy=expected,**row))
                    raw.append(dict(case='legacy_'+name,variant=variant,model=model,command=command,
                                    stdout=p.stdout,stderr=p.stderr))
    assert before == {str(p):digest(p) for p in tracked}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir/'cases.jsonl').write_text(''.join(json.dumps(r,ensure_ascii=False,allow_nan=False)+'\n' for r in records))
    (args.output_dir/'raw.jsonl').write_text(''.join(json.dumps(r,ensure_ascii=False)+'\n' for r in raw))
    (args.output_dir/'legacy.jsonl').write_text(''.join(json.dumps(r,ensure_ascii=False)+'\n' for r in legacy))
    (args.output_dir/'boundary_pf_only.patch').write_text(patch_text)
    (args.output_dir/'manifest.json').write_text(json.dumps(dict(source_commit=subprocess.check_output(
        ['git','-C',str(src),'rev-parse','HEAD'],text=True).strip(),files_sha256=before,temperature_c=37,
        dangles=0,api_dangles=2,observations=len(records),successful_accounting_checks=len(records)-2,
        expected_missing_structures=2,legacy_local_checks=len(legacy),
        scope='development; unchanged fixed libraries; scratch-only legacy PF boundary experiment; no installed model change'),indent=2)+'\n')
    print(f'PASS: {len(records)-2} finite structure accounting checks; 2 expected WIP omissions; source/library hashes unchanged')
    for r in records:
        print(r['model'],r['name'],r['observed_energy'],r['manual_energy'],r['si_isolated_energy'])
    print(f'PASS: {len(legacy)} legacy local checks; scratch-only one-character WIP boundary correction restores expected weights')


if __name__ == '__main__':
    main()
