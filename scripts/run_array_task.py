#!/usr/bin/env python3
"""Execute the SGE_TASK_ID-th matched batch in one tier; never submit jobs."""

import argparse
import json
import os
from pathlib import Path
import platform
import re
import signal
import sys

from plan_array import exclusive, load_plan
from run_cli import digest, execute, experiment_path, now, verify_lock, write_json


def run_task(plan_path, tier, task_id=None, software=None, build=None, resume=False, through_beam=None, scratch_dir=None):
    identifier = os.environ.get("SGE_TASK_ID")
    if identifier is None and task_id is not None:
        identifier = str(task_id)
    if identifier is None or not re.fullmatch(r"[1-9][0-9]*", identifier):
        raise ValueError("SGE_TASK_ID (or development --task-id) must be a positive integer")
    if task_id is not None and int(identifier) != task_id:
        raise ValueError("task-id differs from SGE_TASK_ID")
    plan = load_plan(plan_path, task_selector=(tier,int(identifier)))
    if plan.get("pending_scaffolds"):
        raise ValueError("pending scaffolds: finish G array, aggregate it, then regenerate a verified plan")
    if plan["config"]["purpose"] == "formal" and not (os.environ.get("JOB_ID") and os.environ.get("SGE_TASK_ID")):
        raise ValueError("formal tasks require JOB_ID and SGE_TASK_ID from an array allocation")
    matches = [t for t in plan["tasks"] if t["tier"] == tier and t["task_id"] == int(identifier)]
    if len(matches) != 1:
        raise ValueError("task is outside the saved tier array")
    task = matches[0]
    directory = Path(plan["output_dir"]) / task["directory"]
    with exclusive(directory / ".execution.lock"):
        metadata_path = directory / "task.json"
        metadata = json.loads(metadata_path.read_text())
        if (metadata["attempts"] or (directory / "run").exists()) and not resume:
            raise ValueError("task already attempted; use --resume")
        software_value = software or plan["software_root"] or os.environ.get("LCP_SOFTWARE_ROOT")
        build_value = build or plan["build_dir"] or os.environ.get("LCP_BUILD_DIR")
        software = Path(software_value).resolve() if software_value else None
        build = Path(build_value).resolve() if build_value else None
        attempt = {"attempt_number": len(metadata["attempts"]) + 1, "started_utc": now(),
                   "finished_utc": None, "status": "running", "resume": resume,
                   "global_task_number": task["global_task_number"], "task_id": task["task_id"],
                   "tier": tier, "JOB_ID": os.environ.get("JOB_ID"),
                   "SGE_TASK_ID": os.environ.get("SGE_TASK_ID"), "hostname": platform.node(),
                   "requested_resource": task["requested_resource"], "h_rt": task["h_rt"],
                   "workers": task["workers"], "software_root": str(software) if software else None,
                   "build_dir": str(build) if build else None, "through_beam": through_beam}
        metadata["attempts"].append(attempt)
        metadata.update(status="running", started_utc=attempt["started_utc"], finished_utc=None,
                        global_task_number=task["global_task_number"], requested_resource=task["requested_resource"],
                        h_rt=task["h_rt"], workers=task["workers"], array_plan_sha256=digest(plan_path))
        write_json(metadata_path, metadata)
        context = {key: metadata[key] for key in ("tier", "task_id", "global_task_number",
                                                 "requested_resource", "h_rt", "workers", "array_plan_sha256")}
        output = directory / "run"
        try:
            if software is None or build is None:
                raise ValueError("supply software-root/build-dir or LCP_SOFTWARE_ROOT/LCP_BUILD_DIR")
            for value, saved in ((software, plan["software_root"]), (build, plan["build_dir"])):
                if saved and str(value) != saved:
                    raise ValueError("software/build paths differ from verified plan")
            config = json.loads(Path(task["config_path"]).read_text())
            verify_lock(Path(plan["output_dir"]) / "software_lock.json", software, build,
                        config["programs"], software / config["energy_parameter_file"])
            if output.exists() and not (output / "run.json").exists():
                files = list(output.iterdir())
                if not resume or any(not p.is_file() or p.name not in
                                     {"config.json", "config.json.tmp", "run.json.tmp"} for p in files):
                    raise ValueError("incomplete task startup has no run manifest; inspect before restarting")
                # No CLI can start before the initial manifest is published.
                for file in files:
                    file.unlink()
                output.rmdir()
            code = execute(Path(task["config_path"]), software, build, output,
                           resume=resume and (output / "run.json").is_file(),
                           array_context=context, workers=task["workers"], through_beam=through_beam,
                           scratch_dir=scratch_dir)
            status = json.loads((output / "run.json").read_text())["status"]
        except KeyboardInterrupt:
            code, status = 130, "interrupted"
        except (ValueError, KeyError, OSError) as exc:
            code, status = 2, "failed"
            attempt["error"] = str(exc)
            print("error: " + str(exc), file=sys.stderr)
        attempt.update(finished_utc=now(), status=status, returncode=code)
        metadata.update(finished_utc=attempt["finished_utc"], status=status)
        write_json(metadata_path, metadata)
        return code


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True)
    parser.add_argument("--tier", required=True)
    parser.add_argument("--task-id", type=int, help="development test fallback if SGE_TASK_ID is absent")
    parser.add_argument("--software-root", type=Path)
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--scratch-dir", type=Path, help="explicit node-local directory if TMPDIR is unavailable")
    parser.add_argument("--through-beam", type=int, choices=(20,50,100,200,300,400),
                        help="complete only widths up to this stage; review, then --resume the next stage")
    args = parser.parse_args()

    def terminate(signum, frame):
        raise KeyboardInterrupt

    signal.signal(signal.SIGTERM, terminate)
    return run_task(experiment_path(args.plan), args.tier, args.task_id,
                    args.software_root, args.build_dir, args.resume, args.through_beam, args.scratch_dir)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, OSError) as exc:
        print("error: " + str(exc), file=sys.stderr)
        sys.exit(2)
