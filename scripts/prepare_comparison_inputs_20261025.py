#!/usr/bin/env python3
"""Select seeded real-RNA cohorts and freeze existing E1/E2 scaffold bytes."""
import hashlib
import io
import json
from pathlib import Path
import random
import tarfile

from run_cli import ROOT, digest, write_json
from task_storage import ArchiveReader


def prepare(output):
    if output.exists():
        raise ValueError('choose a fresh input directory')
    cohort_path = ROOT/'manifests/cohorts/bprna-full-20261011.json'
    cohort = json.loads(cohort_path.read_text())
    dataset_path = ROOT/cohort['scaffold_dataset']
    dataset = json.loads(dataset_path.read_text())
    if len(cohort['records']) != 26063 or not dataset['software_unchanged']:
        raise ValueError('unexpected E1 cohort/scaffold provenance')
    by_id = {r['id']:r for r in dataset['records']}
    e1_plans = [ROOT/'results/runs/e1-bprna-short-89baea52-20261018/array_plan.json',
                ROOT/'results/runs/e1-bprna-long-89baea52-cpu16-20261018/array_plan.json']
    # The adopted E1 plans independently certify both sequence and supplied G.
    adopted = {}
    for path in e1_plans:
        plan = json.loads(path.read_text())
        adopted.update({r['name']:r for r in plan['inputs']})
    rng = random.Random(20261025)
    selected, bins = [], []
    for center in (250, 500, 1000, 2000, 4000):
        low, high = center*9//10, center*11//10
        pool = sorted((r for r in cohort['records'] if low<=r['length']<=high), key=lambda r:r['id'])
        chosen = sorted(rng.sample(pool, min(20,len(pool))), key=lambda r:(r['length'],r['id']))
        bins.append(dict(center=center, min_length=low, max_length=high, candidates=len(pool), selected=len(chosen)))
        for r in chosen:
            g, e1 = by_id[r['id']], adopted[r['id']]
            if (g['scaffold_status']!='success' or g['sequence_sha256']!=r['sequence_sha256']
                    or g['input_sha256']!=e1['sha256']):
                raise ValueError('E1 input/scaffold mismatch: '+r['id'])
            selected.append(dict(name=r['id'], group='bprna_'+str(center), length=r['length'],
                                 source_archive=e1['input_archive'], source_member=e1['input_member'],
                                 sha256=e1['sha256'], sequence_sha256=r['sequence_sha256'],
                                 scaffold_sha256=g['scaffold_sha256'], family=r['id'].split('_')[1]))
    e2_path = ROOT/'configs/linearknot_4f97fd3_20261021/e2-rnacentral-4f97fd3-20261021.json'
    for r in json.loads(e2_path.read_text())['inputs']:
        selected.append(dict(name=r['name'],group='rnacentral',length=r['length'],
                             source_archive=r['input_archive'],source_member=r['input_member'],sha256=r['sha256']))
    output.mkdir(parents=True)
    readers = {}
    with tarfile.open(output/'inputs.tar','w') as archive:
        for r in selected:
            path = ROOT/r['source_archive']
            if path not in readers:
                readers[path]=ArchiveReader(path)
            data = readers[path].read(r['source_member'])
            if hashlib.sha256(data).hexdigest()!=r['sha256']:
                raise ValueError('frozen input changed: '+r['name'])
            name, sequence, scaffold = data.decode().splitlines()[:3]
            if len(sequence)!=r['length'] or len(scaffold)!=r['length'] or set(sequence)-set('ACGU'):
                raise ValueError('invalid RNA: '+r['name'])
            for key,value in [('sequence_sha256',sequence),('scaffold_sha256',scaffold)]:
                sha=hashlib.sha256(value.encode()).hexdigest()
                if r.get(key,sha)!=sha:
                    raise ValueError('input component changed: '+r['name'])
                r[key]=sha
            r.update(path=str((output/'inputs'/ (r['name']+'.txt')).relative_to(ROOT)),
                     input_archive=str((output/'inputs.tar').relative_to(ROOT)),input_member='inputs/'+r['name']+'.txt')
            info=tarfile.TarInfo(r['input_member']);info.size=len(data)
            archive.addfile(info,io.BytesIO(data))
    manifest=dict(schema_version=1,seed=20261025,
                  selection='inclusive integer ±10%; ID-sort each pool; random.Random(seed), sample without replacement in ascending center order; then length/ID-sort selected',
                  requested_per_bin=20,bins=bins,records=selected,total=len(selected),
                  scaffold_method=dataset['scaffold_method'],scaffold_settings=dataset['scaffold_settings'],
                  source_files_sha256={str(p.relative_to(ROOT)):digest(p) for p in [cohort_path,dataset_path,e2_path,*e1_plans]},
                  inputs_tar_sha256=digest(output/'inputs.tar'))
    write_json(output/'manifest.json',manifest)
    print(json.dumps(dict(total=len(selected),bins=bins,manifest_sha256=digest(output/'manifest.json')),ensure_ascii=False))


if __name__=='__main__':
    prepare(ROOT/'data/comparison-resources-20261025')
