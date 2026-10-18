#!/usr/bin/env python3
"""Recount archived structures and recompute distribution checks without API calls."""
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tarfile

TOOL = Path(__file__).resolve().parent
ROOT = TOOL.parents[2]
HERE = ROOT / "results/local/e0-enumeration-89baea52-20261018"
WIDTHS = [1, 2, 3, 5]
CODES = 'SHBIMEP'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def pairs(text):
    stack, result = [], []
    for i, char in enumerate(text, 1):
        if char in '([':
            stack.append(i)
        elif char in ')]':
            assert stack
            result.append((stack.pop(), i))
        else:
            assert char == '.'
    assert not stack
    return result


def values(lines):
    return {c[1]: float(c[2]) for c in lines if c[0] == 'value'}


def check_distribution(case, lines, legacy=False):
    n = len(case['sequence'])
    vals = values(lines)
    structures = [c for c in lines if c[0] == 'structure']
    assert len({c[1] for c in structures}) == len(structures) == vals['structures']
    assert sum(int(c[2]) for c in structures) == vals['derivations']
    assert max(int(c[2]) for c in structures) == vals['max_count']
    assert sum(int(c[2]) > 1 for c in structures) == vals['ambiguous']
    if legacy:
        # The control tests integer ambiguity; it need not define the same
        # structure distribution or satisfy default-policy ownership checks.
        structure_map = {c[1]: int(c[2]) for c in structures}
        assert structure_map['.[....][....].'] == 2 and vals['overflow'] == 0
        assert any(c[0] == 'witness' and c[1] == '.[....][....].' for c in lines)
        return vals, structure_map
    fixed_pairs = pairs(case['scaffold'])
    fixed = {i for pair in fixed_pairs for i in pair}
    logs = [-float(c[6]) / vals['rt_kcal_mol'] for c in structures]
    zmax = max(logs)
    z = zmax + math.log(math.fsum(math.exp(w-zmax) for w in logs))
    assert abs(z-vals['structure_logZ']) < 1e-12
    assert abs(abs(vals['structure_logZ']-vals['exact_logZ'])-vals['logZ_error']) < 1e-25
    deriv_logs = [float(c[3]) for c in structures]
    deriv_max = max(deriv_logs)
    deriv_z = deriv_max + math.log(math.fsum(math.exp(w-deriv_max) for w in deriv_logs))
    assert abs(deriv_z-vals['derivation_logZ']) < 1e-12
    assert abs(vals['derivation_logZ']-vals['exact_logZ']) == vals['grouping_logZ_error']
    bpp, q, cp, acc = {}, [0.]*(n+1), {}, {}
    crossing = duplicate = missing = on_paired = 0
    energy_error = 0.
    for c, lw in zip(structures, logs):
        gp, labels = c[1], c[7]
        generated = pairs(gp)
        free_paired = {i for pair in generated for i in pair}
        assert not (fixed & free_paired)
        assert len(free_paired) == 2*len(generated)
        assert all(case['sequence'][i-1]+case['sequence'][j-1] in ['AU','UA','GC','CG','GU','UG'] for i,j in generated)
        assert all(j-i > 3 for i,j in generated)
        paired = fixed | free_paired
        crossing += any(i<k<j<l or k<i<l<j for i,j in generated for k,l in fixed_pairs)
        p = math.exp(lw-vals['structure_logZ'])
        energy_error = max(energy_error, abs(lw-float(c[4])))
        for ij in [*fixed_pairs, *generated]:
            bpp[ij] = bpp.get(ij, 0.) + p
        for i in range(1, n+1):
            duplicate += labels[i-1] == '!'
            missing += i not in paired and labels[i-1] == '-'
            on_paired += i in paired and labels[i-1] != '-'
            if i not in paired:
                q[i] += p
            code = 'S' if i in paired else labels[i-1]
            assert code in CODES
            cp[i, code] = cp.get((i, code), 0.) + p
        for width in WIDTHS:
            for i in range(1, n-width+2):
                if all(t not in paired for t in range(i, i+width)):
                    acc[width, i] = acc.get((width, i), 0.) + p
    assert crossing == vals['crossing_structures']
    assert duplicate == vals['ownership_duplicate'] == 0
    assert missing == vals['ownership_missing'] == 0
    assert on_paired == vals['ownership_on_paired'] == 0
    assert energy_error == vals['energy_error']
    assert vals['overflow'] == 0
    errors = {key: 0. for key in ['bpp_error','unpaired_error','context_error','access_error']}
    errors.update({f'context_error_{c}': 0. for c in CODES})
    errors.update({f'access_error_{w}': 0. for w in WIDTHS})
    seen = {key: set() for key in ['pair','unpaired','context','access']}
    for c in lines:
        typ = c[0].removesuffix('_compare')
        if not c[0].endswith('_compare'):
            continue
        if typ == 'pair':
            key = int(c[1]), int(c[2]); expected = bpp.get(key, 0.); metric = 'bpp_error'
        elif typ == 'unpaired':
            key = int(c[1]); expected = q[key]; metric = 'unpaired_error'
        elif typ == 'context':
            key = int(c[1]), c[2]; expected = cp.get(key, 0.); metric = 'context_error'
        elif typ == 'access':
            key = int(c[1]), int(c[2]); expected = acc.get(key, 0.); metric = 'access_error'
        else:
            raise AssertionError(typ)
        assert key not in seen[typ]
        seen[typ].add(key)
        assert abs(expected-float(c[-2])) < 1e-13
        error = abs(float(c[-2])-float(c[-1]))
        errors[metric] = max(errors[metric], error)
        if typ in ['context', 'access']:
            submetric = metric + '_' + str(c[2] if typ == 'context' else c[1])
            errors[submetric] = max(errors[submetric], error)
    assert bpp.keys() <= seen['pair']
    assert seen['unpaired'] == set(range(1,n+1))
    assert seen['context'] == {(i,c) for i in range(1,n+1) for c in CODES}
    assert seen['access'] == {(w,i) for w in WIDTHS for i in range(1,n-w+2)}
    assert all(error == vals[key] for key,error in errors.items())
    if not legacy:
        assert vals['structures'] == vals['derivations'] and vals['max_count'] == 1
        assert vals['context_ambiguous'] == vals['allowed_absent'] == 0
        assert max(errors.values()) < 1e-9
    return vals, {c[1]: int(c[2]) for c in structures}


def main():
    manifest = json.loads((HERE / 'manifest.json').read_text())
    summary = json.loads((HERE / 'summary.json').read_text())
    report = json.loads((HERE / 'report.json').read_text())
    cases = {c['id']: c for c in json.loads((TOOL.parent / 'baseline_15aa7742/inputs/cases.json').read_text())}
    previous = {r['id']: r for r in json.loads((TOOL.parent / 'baseline_15aa7742/report.json').read_text())['enumeration']}
    rows = {r['id']: r for r in report['enumeration']}
    deltas = {r['id']: r for r in report['comparison']}
    assert manifest['before'] == manifest['after']
    for name, digest in manifest['after']['sha256'].items():
        assert sha(ROOT / name) == digest, name
    for name, digest in manifest['artifacts_sha256'].items():
        assert sha(HERE / name) == digest, name
    source = ROOT / 'deps/LinearCParty-89baea52'
    assert subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip() == manifest['after']['commit']
    assert not subprocess.check_output(['git','-C',str(source),'status','--porcelain'],text=True).strip()
    runs = [json.loads(line) for line in (HERE / 'runs.jsonl').read_text().splitlines()]
    assert len(runs) == summary['processes'] == 1459
    assert all(r['returncode'] == 0 and not r['timed_out'] for r in runs)
    members = {}
    for r in runs:
        members[r['stdout_member']] = r['stdout_sha256']
        members[r['stderr_member']] = r['stderr_sha256']
    checked, structure_count, added_count, seen_members = {}, 0, 0, set()
    with tarfile.open(HERE / 'raw.tar.gz','r|gz') as archive:
        for member in archive:
            assert member.isfile() and member.name not in seen_members
            seen_members.add(member.name)
            data = archive.extractfile(member).read()
            if member.name.startswith('build/'):
                if member.name == 'build/enumerate':
                    assert hashlib.sha256(data).hexdigest() == manifest['binary_sha256']
                continue
            assert hashlib.sha256(data).hexdigest() == members[member.name]
            if member.name.endswith('.stderr'):
                assert not data
                continue
            lines = [line.split('\t') for line in data.decode().splitlines()]
            label = member.name.removesuffix('.tsv')
            cid, mode = label.rsplit('_',1)
            if label.endswith('_pre_wip'):
                cid = label.removesuffix('_pre_wip')
                vals = values(lines)
                pre = {c[1]:int(c[2]) for c in lines if c[0]=='pre_wip_structure'}
                assert len(pre) == vals['structures'] == previous[cid]['structures']
                assert sum(pre.values()) == vals['derivations'] == previous[cid]['derivations']
                assert all(c == checked[cid][gp] for gp,c in pre.items())
                added = sorted(checked[cid].keys()-pre.keys())
                assert added == deltas[cid]['added_structures']
                assert len(added) == deltas[cid]['delta_structures'] == deltas[cid]['added_by_wip']
                added_count += len(added)
            else:
                vals, structure_map = check_distribution(cases[cid], lines, legacy=mode=='legacy')
                if mode == 'default':
                    assert all(vals[k] == rows[cid][k] for k in vals)
                    checked[cid] = structure_map
                    structure_count += len(structure_map)
                else:
                    assert mode == 'legacy' and cid == 'closure14'
                    assert structure_map['.[....][....].'] == 2
    assert members.keys() <= seen_members
    assert set(checked) == cases.keys() == rows.keys() == deltas.keys()
    assert checked['closure14']['.[....][....].'] == 1
    assert structure_count == summary['totals']['structures']
    assert added_count == summary['added_by_wip']
    for k, value in summary['totals'].items():
        assert sum(r[k] for r in rows.values()) == value
    for k, value in summary['maxima'].items():
        assert max(r[k] for r in rows.values()) == value
    assert sum(d['delta_structures']!=0 for d in deltas.values()) == summary['changed_inputs']
    assert summary['verification_complete'] and summary['integrity_ok'] and not summary['failures']
    result = dict(ok=True, cases=len(checked), processes=len(runs), structures=structure_count,
                  added_structures=added_count, archive_members=len(seen_members),
                  independent_distribution_reaggregation=True, protected_hashes_match=True,
                  validator_sha256=sha(Path(__file__).resolve()))
    (HERE / 'validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))


if __name__ == '__main__':
    main()
