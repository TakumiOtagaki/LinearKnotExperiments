#!/usr/bin/env python3
"""Recount ordered raw draws and check saved validation results independently."""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import tarfile


def parse_pairs(text):
    stacks = {'(': [], '[': []}
    result = set()
    for i, c in enumerate(text, 1):
        if c in stacks:
            stacks[c].append(i)
        elif c in ')]':
            opening = '(' if c == ')' else '['
            result.add((stacks[opening].pop(), i))
        elif c != '.':
            raise ValueError('unexpected structure character')
    assert not any(stacks.values())
    assert len({i for pair in result for i in pair}) == 2 * len(result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('run', type=Path)
    args = p.parse_args()
    run = args.run
    config = json.loads((run/'config.json').read_text())
    summary = json.loads((run/'summary.json').read_text())
    rows = [json.loads(line) for line in (run/'results.jsonl').read_text().splitlines()]
    failures = []
    checked_draws = checked_members = 0
    max_mfe_error = 0.
    computed_events = 0
    expected_conditions = {(c['id'],m) for c in config['cases'] for m in config['modes']}
    assert {(r['case'],r['mode']) for r in rows} == expected_conditions
    assert len(rows)==len(expected_conditions)==96
    with tarfile.open(run/'raw.tar.gz', 'r|gz') as archive:
        references = {}
        by_prefix = {r['case']+'/'+r['mode']:r for r in rows}
        for member in archive:
            prefix, _, name = member.name.rpartition('/')
            if prefix not in by_prefix or name not in by_prefix[prefix]['raw']:
                continue
            r=by_prefix[prefix]
            metadata=r['raw'][name]
            content=archive.extractfile(member).read()
            assert hashlib.sha256(content).hexdigest()==metadata['sha256']
            assert len(content)==metadata['bytes']
            checked_members+=1
            if name=='reference.tsv':
                mass={}
                for line in content.decode().splitlines():
                    fields=line.split('\t')
                    if fields[0]=='structure':mass[fields[1]]=(int(fields[2]),float(fields[3]))
                references[prefix]=mass
            if name!='draws.txt':continue
            counts=Counter(content.decode().splitlines())
            expected={f['extension']:f for f in r['metrics']['frequencies']}
            assert set(counts)<=set(expected)
            assert sum(counts.values())==r['sample_count']==100000
            checked_draws+=sum(counts.values())
            for key,f in expected.items():
                assert f['count']==counts.get(key,0)
                assert references[prefix][key][0]==f['derivations']==1
                assert abs(references[prefix][key][1]-f['log_weight'])<1e-12
            probabilities={key:math.exp(weight-r['metrics']['reference']['enumerated_logZ'])
                           for key,(_,weight) in references[prefix].items()}
            assert abs(sum(probabilities.values())-1.)<1e-10
            tv=sum(abs(counts.get(k,0)/100000-p) for k,p in probabilities.items())/2
            assert abs(tv-r['metrics']['total_variation'])<1e-12
            observed_bpp=Counter();reference_bpp=Counter()
            fixed=parse_pairs(r['scaffold'])
            for gp in probabilities:
                added=parse_pairs(gp)
                assert not ({i for bp in fixed for i in bp}&{i for bp in added for i in bp})
                for pair in fixed|added:
                    observed_bpp[pair]+=counts.get(gp,0)/100000
                    reference_bpp[pair]+=probabilities[gp]
            assert len(r['metrics']['bpp'])==r['n']*(r['n']-1)//2
            for f in r['metrics']['bpp']:
                key=(f['i'],f['j'])
                assert abs(observed_bpp[key]-f['sample'])<1e-12
                assert abs(reference_bpp[key]-f['enumerated'])<1e-12
            computed_events+=len(expected)+r['n']*(r['n']-1)//2+sum(r['n']-w+1 for w in [1,3,5])
            if r['mode']=='exact':
                energies={f['extension']:f['energy'] for f in r['metrics']['frequencies']}
                v=r['mfe_preparation']['default']
                canonical=v['extension'].translate(str.maketrans('()','[]'))
                error=max(abs(v['energy']-min(energies.values())),abs(v['energy']-energies[canonical]))
                max_mfe_error=max(max_mfe_error,error)
                assert error<1e-9
    assert checked_members==288 and checked_draws==9600000
    assert computed_events==summary['events']
    epsilon=math.sqrt(math.log(2*computed_events/.01)/(2*100000))
    assert abs(epsilon-summary['simultaneous_hoeffding_epsilon'])<1e-15
    assert summary['all_checks_pass']
    report=dict(valid=True,conditions=len(rows),raw_members_checked=checked_members,
                ordered_draws_recounted=checked_draws,event_count=computed_events,
                simultaneous_epsilon=epsilon,default_mfe_vs_enumerated_min_max_error=max_mfe_error,
                source_and_binary_hashes_unchanged=summary['source_unchanged'] and summary['hashes_unchanged'])
    (run/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
