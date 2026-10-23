#!/usr/bin/env python3
"""Classify every saved enumerated pair map without grammar or rule IDs."""
import argparse
from collections import Counter
from functools import lru_cache
import json
from pathlib import Path
import re
import sys
import tarfile
import time

from rerun_pseudoloop_small_20261023 import dump,sha

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/local/pseudoloop-0.2.0-20261023/enumeration'


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--source',type=Path,required=True)
    args=ap.parse_args();sys.path.insert(0,str(args.source/'test/tools'))
    from pseudoloop_context_reference import classify,pairs_from_db
    @lru_cache(maxsize=32768)
    def reference(n,pairs):return classify(n,pairs)[0]
    rows=json.loads((OUT/'report.json').read_text())['enumeration']
    cases={r['id']:r for r in rows};seen=set();counts=Counter();bad=[];t=time.monotonic()
    with tarfile.open(OUT/'raw_outputs.tar.gz','r|gz') as tf:
        for m in tf:
            name=m.name.removesuffix('.stdout')
            if not m.name.endswith('.stdout') or name not in cases:continue
            case=cases[name];seen.add(name);fixed=pairs_from_db(case['scaffold'])
            for line in tf.extractfile(m).read().decode().splitlines():
                fields=line.split('\t')
                if fields[0]!='structure':continue
                pairs=tuple(sorted(fixed | pairs_from_db(fields[1])))
                got=fields[7].replace('-','S');want=reference(len(got),pairs)
                counts['structure_instances']+=1;counts['base_instances']+=len(got)
                counts['derivations']+=int(fields[2])
                if got!=want:
                    bad.append({'case':name,'sequence':case['sequence'],'scaffold':case['scaffold'],
                        'pairs':pairs,'emitted':got,'geometry':want})
                for run in filter(None,re.split('S+',got)):
                    counts['fully_unpaired_window_instances']+=len(run)*(len(run)+1)//2
                    if len(set(run))!=1:counts['mixed_unpaired_runs']+=1
    assert seen==set(cases)
    result={'cases':len(seen),'counts':dict(counts),'counterexamples':bad,'elapsed_seconds':time.monotonic()-t,
        'enumeration_archive_sha256':sha(OUT/'raw_outputs.tar.gz'),
        'reference_sha256':sha(args.source/'test/tools/pseudoloop_context_reference.py'),
        'script_sha256':sha(__file__),'status':'passed' if not bad and not counts['mixed_unpaired_runs'] else 'failed'}
    dump(OUT/'geometry_validation.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='counterexamples'},indent=2))
    assert result['status']=='passed',f'{len(bad)} geometric counterexamples; see saved JSON'


if __name__=='__main__':main()
