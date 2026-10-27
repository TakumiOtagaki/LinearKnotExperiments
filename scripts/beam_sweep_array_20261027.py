#!/usr/bin/env python3
"""Run standard E2 packed tasks with node-local staging and BPP gzip members."""
import argparse
import gzip
import json
import os
from pathlib import Path
import shutil
import signal
import tarfile
import time

from build_comparison_resources_20261025 import attach_git
from comparison_resources_array_20261025 import check_bundle
from plan_array import exclusive, load_plan
import run_cli
from run_array_task import run_task
from task_storage import TaskStore


class GzipBppStore(TaskStore):
    def commit(self, row, files):
        files = dict(files)
        started = time.perf_counter()
        for name, path in list(files.items()):
            if name.endswith('/bpp.tsv'):
                compressed = path.with_suffix('.tsv.gz')
                with path.open('rb') as src, compressed.open('wb') as dst:
                    with gzip.GzipFile(filename='', mode='wb', fileobj=dst, mtime=0, compresslevel=6) as out:
                        shutil.copyfileobj(src, out, length=1024*1024)
                row['bpp_uncompressed_bytes'] = path.stat().st_size
                row['bpp_uncompressed_sha256'] = run_cli.digest(path)
                row['bpp_gzip_bytes'] = compressed.stat().st_size
                del files[name]
                files[name+'.gz'] = compressed
        row['bpp_compress_seconds'] = time.perf_counter()-started
        # Timed CLI already exited. Compression/checkpoint time is task overhead.
        return super().commit(row, files)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan', type=Path, required=True)
    p.add_argument('--tier', required=True)
    p.add_argument('--task-id', type=int)
    p.add_argument('--resume', action='store_true')
    p.add_argument('--dry-run', action='store_true')
    a = p.parse_args()
    task_id = a.task_id or int(os.environ['SGE_TASK_ID'])
    plan = load_plan(a.plan.resolve(), task_selector=(a.tier, task_id))
    bundle = check_bundle(plan)
    parent = Path('/tmp')/('beam-sweep-'+str(os.getuid()))/Path(plan['output_dir']).name/a.tier
    parent.mkdir(parents=True, exist_ok=True)
    work = parent/str(task_id)
    with exclusive(parent/(str(task_id)+'.stage.lock')):
        if work.exists():
            shutil.rmtree(work)
        work.mkdir()
        try:
            with tarfile.open(bundle) as tar:
                tar.extractall(work)
            attach_git(work/'source', Path(plan['software_lock']['source_repositories']['linearknot']),
                       plan['software_lock']['source_commit'])
            cfg = json.loads(next(Path(t['config_path']).read_text() for t in plan['tasks']
                                  if (t['tier'], t['task_id']) == (a.tier, task_id)))
            run_cli.verify_lock(Path(plan['output_dir'])/'software_lock.json', work/'source', work/'runtime',
                                cfg['programs'], work/'source'/cfg['energy_parameter_file'])
            if a.dry_run:
                print(json.dumps(dict(status='passed', tier=a.tier, task_id=task_id,
                                      staged_build_verified=True, cli_executed=False)))
                return 0
            run_cli.TaskStore = GzipBppStore
            def terminate(signum, frame):
                raise KeyboardInterrupt
            signal.signal(signal.SIGTERM, terminate)
            return run_task(a.plan.resolve(), a.tier, task_id, work/'source', work/'runtime',
                            resume=a.resume, scratch_dir=work)
        finally:
            shutil.rmtree(work)


if __name__ == '__main__':
    raise SystemExit(main())
