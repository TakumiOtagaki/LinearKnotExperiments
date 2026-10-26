#!/usr/bin/env python3
"""Build three pinned tools in node-local scratch; save a portable tar and lock."""
import argparse
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile

from run_cli import ROOT, build_settings, digest, now, write_json

LK = 'f93b55ed9749dd1c198cd1f1809aece50f12f7d8'
SPARK = '3105c311ca6243a31788e74ac33640435e53b62c'
CPARTY = '800e95af9d1a97555233e2c7b37d73dbaa1a92b4'
OBJECTS = ('part_func', 'dot_plot', 'mea', 'centroid', 'W_final', 'pseudo_loop',
           's_energy_matrix', 'Hotspot', 'sparse_tree', 'SHAPE')


def attach_git(source, repository, commit):
    """Index archive bytes against read-only local objects, without a checkout."""
    env = dict(os.environ, GIT_OPTIONAL_LOCKS='0')
    subprocess.run(['git', 'init', '-q', str(source)], check=True, env=env)
    objects = subprocess.check_output(['git', '-C', str(repository), 'rev-parse',
                                       '--path-format=absolute', '--git-path', 'objects'], text=True).strip()
    (source/'.git/objects/info/alternates').write_text(objects+'\n')
    subprocess.run(['git', '-C', str(source), 'update-ref', 'HEAD', commit], check=True, env=env)
    subprocess.run(['git', '-C', str(source), 'read-tree', commit], check=True, env=env)
    status = subprocess.check_output(['git', '-C', str(source), 'status', '--porcelain'], text=True, env=env)
    if status:
        raise ValueError('archive does not match the pinned Git tree: '+status[:500])


def extract_archive(repository, revision, destination, run):
    destination.mkdir()
    archive = destination.parent/(destination.name+'.tar')
    run(['git', '-C', str(repository), 'archive', '--format=tar', '-o', str(archive), revision])
    sha = digest(archive)
    with tarfile.open(archive) as stream:
        stream.extractall(destination)
    archive.unlink()
    return sha


def build(args):
    work = args.scratch.resolve()
    if str(work).startswith('/gs/bs/') or work.exists():
        raise ValueError('choose a fresh node-local scratch directory')
    work.mkdir(parents=True)
    args.out.mkdir(parents=True, exist_ok=False)
    commands = []
    log = work/'build.log'

    def run(command):
        command = list(map(str, command))
        commands.append(command)
        with log.open('a') as output:
            output.write(json.dumps(command)+'\n'); output.flush()
            subprocess.run(command, check=True, stdout=output, stderr=output,
                           env=dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                                    MKL_NUM_THREADS='1', GIT_OPTIONAL_LOCKS='0'))

    oracle_lock = json.loads((ROOT/'manifests/software-cparty-oracle-800e95a-6patches-20261018.json').read_text())
    spark_lock = json.loads((ROOT/'manifests/software-spark-3105c311-20261012.json').read_text())
    lk = args.linearknot.resolve()
    if subprocess.check_output(['git', '-C', str(lk), 'rev-parse', 'v0.2.0^{commit}'], text=True).strip() != LK:
        raise ValueError('v0.2.0 tag differs')
    sources = {}
    sources['linearknot'] = extract_archive(lk, 'v0.2.0', work/'source', run)
    sources['spark'] = extract_archive(ROOT/'deps/Spark-upstream', SPARK, work/'spark-source', run)
    for name, sha in spark_lock['tracked_files_sha256'].items():
        if digest(work/'spark-source'/name) != sha:
            raise ValueError('Spark differs from 2026-10-12 source lock: '+name)
    sources['cparty_upstream'] = extract_archive(ROOT/'build/cparty-oracle-800e95a-6patches/src',
                                                CPARTY, work/'cparty-source', run)
    patches = []
    for patch in oracle_lock['patches']:
        path = ROOT/patch['path']
        if digest(path) != patch['sha256']:
            raise ValueError('CParty patch changed: '+str(path))
        patches.append(patch)
        run(['git', 'apply', '--unsafe-paths', '--directory='+str(work/'cparty-source'), path])
    for name, sha in oracle_lock['source_files_sha256'].items():
        if digest(work/'cparty-source/src'/name) != sha:
            raise ValueError('patched CParty differs from manuscript oracle: '+name)
    for name, source in [('lk', 'source'), ('spark', 'spark-source'), ('cparty', 'cparty-source')]:
        options = ['-DCMAKE_BUILD_TYPE=Release', '-DCMAKE_EXPORT_COMPILE_COMMANDS=ON']
        if name == 'lk':
            options.append('-DLINEARKNOT_LTO_JOBS=1')
        run(['cmake', '-S', work/source, '-B', work/(name+'-build'), *options])
        target = ['--target', 'linearknot'] if name == 'lk' else []
        run(['cmake', '--build', work/(name+'-build'), *target, '--parallel', str(args.build_workers)])
    runtime = work/'runtime'
    runtime.mkdir()
    archives = [work/'lk-build'/n for n in ('libcparty_scfg2_beam.a', 'libcparty_scfg2_exact.a', 'libRNA.a')]
    run(['g++', '-O2', '-std=c++17', '-I'+str(work/'source/src'), ROOT/'tools/pf_scaling_linearknot.cc',
         *archives, '-o', runtime/'PFScaling'])
    objects = work/'cparty-build/CMakeFiles/CParty.dir/src'
    run(['g++', '-O2', '-std=gnu++17', '-I'+str(work/'cparty-source/src'), ROOT/'tools/e0_legacy_scalar.cc',
         *[objects/(n+'.cc.o') for n in OBJECTS], work/'cparty-build/libRNA.a', '-o', runtime/'E0CPartyOracle'])
    for src, name in [(work/'lk-build/linearknot', 'linearknot'), (work/'spark-build/Spark', 'Spark'),
                      (work/'lk-build/CMakeCache.txt', 'CMakeCache.txt'),
                      (ROOT/'scripts/comparison_resources_command.py', 'Comparison')]:
        shutil.copy2(src, runtime/name)
    (runtime/'Comparison').chmod(0o755)
    provenance = dict(created_utc=now(), source_commits=dict(linearknot=LK, spark=SPARK, cparty=CPARTY),
                      source_archive_sha256=sources, patches=patches, commands=commands,
                      compiler=subprocess.check_output(['g++', '--version'], text=True).splitlines()[0],
                      hostname=os.uname().nodename, build_workers=args.build_workers,
                      cpu_model=subprocess.check_output(['lscpu'], text=True),
                      build_settings={name: build_settings(work/(name+'-build')) for name in ('lk','spark','cparty')},
                      library_sha256={p.name:digest(p) for p in archives},
                      scope='PF inside-only compatible; Spark official MFE; linearknot fold includes MFE+PF, no BPP')
    write_json(work/'build_provenance.json', provenance)
    # Keep complete source and build records in one archive; no per-file shared build.
    with tarfile.open(args.out/'build_raw.tar.gz', 'w:gz') as archive:
        for name in ('source', 'spark-source', 'cparty-source', 'build.log', 'build_provenance.json'):
            archive.add(work/name, arcname=name)
        for name in ('lk', 'spark', 'cparty'):
            for filename in ('CMakeCache.txt', 'compile_commands.json'):
                archive.add(work/(name+'-build')/filename, arcname=name+'-build/'+filename)
    with tarfile.open(args.out/'software.tar.gz', 'w:gz') as archive:
        archive.add(work/'source', arcname='source')
        archive.add(runtime, arcname='runtime')
    extra = [ROOT/'tools/pf_scaling_linearknot.cc', ROOT/'tools/e0_legacy_scalar.cc',
             ROOT/'scripts/comparison_resources_command.py', Path(__file__).resolve(),
             args.out/'software.tar.gz']
    lock = dict(schema_version=1, created_utc=now(), purpose='formal', source_commit=LK,
                source_git_state=dict(commit=LK, dirty=False), build_revision_verified=True,
                verification_method='git archive v0.2.0; rebuild on TSUBAME node-local; old CParty/Spark source hashes checked',
                energy_parameter_file='params/rna_Turner04.par',
                energy_parameter_sha256=digest(work/'source/params/rna_Turner04.par'),
                cmake_cache_sha256=digest(runtime/'CMakeCache.txt'),
                build_settings=build_settings(runtime),
                binaries={p.name:digest(p) for p in runtime.iterdir() if p.name!='CMakeCache.txt'},
                build_files={p.name:digest(p) for p in runtime.iterdir()},
                extra_files=[dict(path=str(p.relative_to(ROOT)),sha256=digest(p)) for p in extra],
                software_bundle=str((args.out/'software.tar.gz').relative_to(ROOT)),
                source_repositories=dict(linearknot=str(lk)), provenance=provenance)
    write_json(args.lock, lock)
    shutil.copy2(work/'build_provenance.json', args.out/'build_provenance.json')
    attach_git(work/'source', lk, LK)
    print(json.dumps(dict(scratch=str(work), lock=str(args.lock), bundle_sha256=digest(args.out/'software.tar.gz'))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scratch', type=Path, required=True)
    parser.add_argument('--linearknot', type=Path, default=ROOT.parent.parent/'LinearKnot')
    parser.add_argument('--build-workers', type=int, default=2)
    parser.add_argument('--out', type=Path, default=ROOT/'results/local/comparison-resources-build-20261025')
    parser.add_argument('--lock', type=Path, default=ROOT/'manifests/software-comparison-resources-v0.2.0-20261025.json')
    build(parser.parse_args())
