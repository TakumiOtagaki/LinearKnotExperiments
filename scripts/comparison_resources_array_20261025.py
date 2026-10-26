#!/usr/bin/env python3
"""Stage the hashed node-local build and run a standard packed array task."""
import argparse
import json
import os
from pathlib import Path
import shutil
import signal
import tarfile

from build_comparison_resources_20261025 import attach_git
from plan_array import exclusive, load_plan
from run_array_task import run_task
from run_cli import ROOT, digest


def check_bundle(plan):
    lock = plan['software_lock']
    bundle = ROOT/lock['software_bundle']
    expected = next(r['sha256'] for r in lock['extra_files'] if ROOT/r['path']==bundle)
    if digest(bundle)!=expected:
        raise ValueError('software bundle changed')
    for r in lock['extra_files']:
        if digest(ROOT/r['path'])!=r['sha256']:
            raise ValueError('frozen dependency changed: '+r['path'])
    return bundle


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True)
    parser.add_argument('--tier',required=True)
    parser.add_argument('--task-id',type=int)
    parser.add_argument('--resume',action='store_true')
    parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args()
    task_id=args.task_id or int(os.environ['SGE_TASK_ID'])
    plan=load_plan(args.plan.resolve(),task_selector=(args.tier,task_id))
    bundle=check_bundle(plan)
    if args.dry_run:
        print(json.dumps(dict(status='passed',task_id=task_id,tier=args.tier,cli_executed=False)))
        return 0
    # Stable node-local paths permit the existing hash-checked resume mechanism
    # to work across allocations; an outer local lock protects staging itself.
    parent=Path('/tmp')/('comparison-resources-'+str(os.getuid()))/Path(plan['output_dir']).name/args.tier
    parent.mkdir(parents=True,exist_ok=True)
    work=parent/str(task_id)
    with exclusive(parent/(str(task_id)+'.stage.lock')):
        if work.exists():
            shutil.rmtree(work)
        work.mkdir()
        try:
            with tarfile.open(bundle) as archive:
                archive.extractall(work)
            attach_git(work/'source',Path(plan['software_lock']['source_repositories']['linearknot']),
                       plan['software_lock']['source_commit'])
            def terminate(signum,frame):
                raise KeyboardInterrupt
            signal.signal(signal.SIGTERM,terminate)
            return run_task(args.plan.resolve(),args.tier,task_id,work/'source',work/'runtime',
                            resume=args.resume,scratch_dir=work)
        finally:
            shutil.rmtree(work)


if __name__=='__main__':
    raise SystemExit(main())
