#!/usr/bin/env python3
"""Save reference context mass/means from the audited node-local overlay DB.

Run after the profile audit and before reanalyze_e1_context removes its scratch
database. The short reference is exact; the long reference remains finite b400.
"""
import argparse
import csv
import json
from pathlib import Path
import pickle
import sqlite3
import tarfile
import zlib

import analyze_accuracy as accuracy
from pseudoloop_context_array_20261023 import ROOT
from run_cli import digest


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plan', type=Path, required=True)
    ap.add_argument('--db', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    plan = json.loads(args.plan.read_text())
    audit = json.loads((args.out/'task_output_validation.json').read_text())
    assert audit['all_passed'] and audit['plan_sha256'] == digest(args.plan)
    db = sqlite3.connect(f'file:{args.db}?mode=ro', uri=True)
    for cohort, archive in plan['input_archives'].items():
        names = {v['member']: n for n,v in plan['inputs'].items() if v['cohort'] == cohort}
        ref = 'exact' if cohort == 'short' else 'beam400_outer'
        path = args.out/f'e1-bprna-{cohort}-context-0.2.0-20261023/reference_context_descriptors.tsv'
        with path.open('w', newline='') as f, tarfile.open(ROOT/archive['path'], 'r|') as tf:
            writer = csv.DictWriter(f, ('input','input_sha256','reference_mode','free_bases',
                'reference_pseudoloop_sum','reference_pseudoloop_mean_Gfree',
                'exact_pseudoloop_mean_Gfree'), delimiter='\t', lineterminator='\n')
            writer.writeheader()
            seen = set()
            for member in tf:
                if member.name not in names:
                    continue
                name = names[member.name]
                scaffold = tf.extractfile(member).read().decode().splitlines()[2]
                row = db.execute('SELECT value FROM profiles WHERE name=? AND mode=?',(name,ref)).fetchone()
                profile = pickle.loads(zlib.decompress(row[0]))
                free = [i for i,c in enumerate(scaffold) if c == '.']
                mean = accuracy.mean([profile['PseudoLoop'][i] for i in free])
                writer.writerow(dict(input=name,input_sha256=plan['inputs'][name]['sha256'],
                    reference_mode=ref,free_bases=len(free),reference_pseudoloop_sum=sum(profile['PseudoLoop']),
                    reference_pseudoloop_mean_Gfree=mean if mean is not None else 'NA',
                    exact_pseudoloop_mean_Gfree=mean if cohort == 'short' and mean is not None else 'NA'))
                seen.add(name)
            assert seen == set(names.values())
    db.close()


if __name__ == '__main__':
    main()
