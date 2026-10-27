#!/usr/bin/env python3
"""Build archived LinearKnot v0.2.1 locally, preserving a relocatable bundle."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile

from build_comparison_resources_20261025 import attach_git, extract_archive
from run_cli import ROOT, build_settings, digest, now, write_json

COMMIT = '2101e5d86959661ba0a3da6d598f8f202e87186b'
BASE = 'beam-sweep-rnacentral-v0.2.1-20261027'


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scratch', type=Path, required=True)
    p.add_argument('--source-repo', type=Path, default=ROOT.parent.parent/'LinearKnot')
    p.add_argument('--build-workers', type=int, default=2)
    a = p.parse_args()
    work, repo = a.scratch.resolve(), a.source_repo.resolve()
    if work.exists() or str(work).startswith('/gs/bs/'):
        raise ValueError('choose a fresh node-local scratch path')
    env = dict(os.environ, GIT_OPTIONAL_LOCKS='0', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    if subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'v0.2.1^{commit}'], env=env, text=True).strip() != COMMIT:
        raise ValueError('v0.2.1 tag changed')
    work.mkdir(parents=True)
    commands = []

    def run(command):
        command = list(map(str, command))
        commands.append(command)
        with (work/'build.log').open('a') as log:
            log.write(json.dumps(command)+'\n'); log.flush()
            subprocess.run(command, env=env, stdout=log, stderr=log, check=True)

    archive_sha = extract_archive(repo, 'v0.2.1', work/'source', run)
    run(['cmake', '-S', work/'source', '-B', work/'build', '-DCMAKE_BUILD_TYPE=Release',
         '-DCMAKE_EXPORT_COMPILE_COMMANDS=ON', '-DLINEARKNOT_LTO_JOBS=1'])
    run(['cmake', '--build', work/'build', '--target', 'linearknot', '--parallel', a.build_workers])
    runtime = work/'runtime'
    runtime.mkdir()
    for name in ('linearknot', 'CMakeCache.txt'):
        shutil.copy2(work/'build'/name, runtime/name)
    version = subprocess.check_output([str(runtime/'linearknot'), '--version'], text=True).strip()
    if '0.2.1' not in version:
        raise ValueError('unexpected binary version: '+version)
    out = ROOT/'results/local'/BASE/'build'
    out.mkdir(parents=True, exist_ok=False)
    provenance = dict(created_utc=now(), source_commit=COMMIT, source_tag='v0.2.1',
                      source_archive_sha256=archive_sha, commands=commands, version=version,
                      compiler=subprocess.check_output(['g++', '--version'], text=True).splitlines()[0],
                      host=os.uname().nodename, build_settings=build_settings(runtime),
                      build_workers=a.build_workers)
    write_json(out/'provenance.json', provenance)
    with tarfile.open(out/'build_raw.tar.gz', 'w:gz') as tar:
        for name in ('build.log', 'build/CMakeCache.txt', 'build/compile_commands.json'):
            tar.add(work/name, arcname=name)
    with tarfile.open(out/'software.tar.gz', 'w:gz') as tar:
        tar.add(work/'source', arcname='source')
        tar.add(runtime, arcname='runtime')
    lock = dict(schema_version=1, created_utc=now(), source_commit=COMMIT,
                source_git_state=dict(commit=COMMIT, dirty=False), build_revision_verified=True,
                verification_method='git archive v0.2.1 and Release build in node-local scratch; source archive indexed against read-only objects',
                source_repositories=dict(linearknot=str(repo)), software_bundle=str((out/'software.tar.gz').relative_to(ROOT)),
                energy_parameter_file='params/rna_Turner04.par', energy_parameter_sha256=digest(work/'source/params/rna_Turner04.par'),
                cmake_cache_sha256=digest(runtime/'CMakeCache.txt'), build_settings=build_settings(runtime),
                binaries=dict(linearknot=digest(runtime/'linearknot')), provenance=provenance,
                extra_files=[dict(path=str((out/'software.tar.gz').relative_to(ROOT)), sha256=digest(out/'software.tar.gz'))])
    write_json(ROOT/'manifests'/('software-'+BASE+'.json'), lock)
    from beam_sweep_metadata_20261027 import add_metadata
    add_metadata(work, out, ROOT/'manifests'/('software-'+BASE+'.json'))
    attach_git(work/'source', repo, COMMIT)
    print(json.dumps(dict(scratch=str(work), version=version, bundle=str(out/'software.tar.gz'))))


if __name__ == '__main__':
    main()
