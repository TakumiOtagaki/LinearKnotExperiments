"""Add scalar logging to a copy of fold.cc; leave archived sources untouched."""
import json
from pathlib import Path
import subprocess
import tarfile

from run_cli import ROOT, digest, write_json


def add_metadata(work, out, lock_path):
    source = (work/'source/src/scfg2/cli/fold.cc').read_text()
    replacements = [
        ('  double mfe =', '  double log_z = 0.0;\n  double mfe ='),
        ('      out.ensemble_free_energy = outside.ensemble_free_energy;',
         '      out.log_z = outside.log_partition_weight;\n      out.ensemble_free_energy = outside.ensemble_free_energy;'),
        ('      out.ensemble_free_energy = inside.ensemble_free_energy;',
         '      out.log_z = inside.log_partition_weight;\n      out.ensemble_free_energy = inside.ensemble_free_energy;'),
        ('    std::cout << ">" << record.name',
         '    std::cerr << "BEAM_SWEEP_SCALARS " << std::setprecision(17)\n'
         '              << "{\\"log_z\\":" << result.log_z\n'
         '              << ",\\"ensemble_free_energy\\":" << result.ensemble_free_energy\n'
         '              << ",\\"mfe\\":" << result.mfe << "}\\n";\n'
         '    std::cout << ">" << record.name'),
    ]
    for (old, new), count in zip(replacements, (1, 2, 2, 1)):
        if source.count(old) != count:
            raise ValueError('unexpected v0.2.1 fold source')
        source = source.replace(old, new)
    source += '\nint main(int argc, char **argv) { return cparty::scfg2::cli::run_fold(argc-1, argv+1); }\n'
    observer = work/'fold_metadata.cc'
    observer.write_text(source)
    command = ['g++', '-O3', '-flto=1', '-std=c++17', '-I'+str(work/'source/src'),
               str(observer), *[str(work/'build'/n) for n in
                               ('libcparty_scfg2_beam.a', 'libcparty_scfg2_exact.a', 'libRNA.a')],
               '-o', str(work/'runtime/FoldMetadata')]
    with (work/'metadata_build.log').open('w') as log:
        subprocess.run(command, check=True, stdout=log, stderr=log)
    with tarfile.open(out/'metadata_build.tar.gz', 'w:gz') as tar:
        for name in ('fold_metadata.cc', 'metadata_build.log'):
            tar.add(work/name, arcname=name)
    with tarfile.open(out/'software.tar.gz', 'w:gz') as tar:
        # .git is node-local provenance, reconstructed at staging time.
        tar.add(work/'source', arcname='source', filter=lambda t: None if '.git' in Path(t.name).parts else t)
        tar.add(work/'runtime', arcname='runtime')
    lock = json.loads(lock_path.read_text())
    lock['binaries']['FoldMetadata'] = digest(work/'runtime/FoldMetadata')
    lock['build_files'] = {n: digest(work/'runtime'/n) for n in ('linearknot', 'FoldMetadata', 'CMakeCache.txt')}
    lock['provenance']['metadata_observer'] = dict(
        generator='scripts/beam_sweep_metadata_20261027.py', command=command,
        generated_source_sha256=digest(observer),
        scope='Unchanged fold.cc calculation and stdout; append full precision log Z/F/MFE to stderr. No additional DP run.')
    for entry in lock['extra_files']:
        if entry['path'] == str((out/'software.tar.gz').relative_to(ROOT)):
            entry['sha256'] = digest(out/'software.tar.gz')
    write_json(lock_path, lock)
    write_json(out/'provenance.json', lock['provenance'])
