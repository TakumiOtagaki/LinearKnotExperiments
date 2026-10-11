#!/usr/bin/env python3
"""Merge task result markers into a run, retaining every missing observation."""

import argparse
from collections import Counter
import fcntl
import json
from pathlib import Path
import sys

from plan_array import exclusive, load_plan
from run_cli import digest, experiment_path, now, write_json


def active_task(directory):
    path = directory / ".execution.lock"
    if not path.exists():
        return False
    with path.open("r") as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return True
        fcntl.flock(stream, fcntl.LOCK_UN)
    return False


def aggregate(plan_path):
    plan = load_plan(plan_path)
    output = Path(plan["output_dir"])
    with exclusive(output / ".aggregate.lock"):
        jobs, tasks, checks = [], [], []
        for task in plan["tasks"]:
            directory = output / task["directory"]
            metadata = json.loads((directory / "task.json").read_text())
            run_file = directory / "run/run.json"
            manifest = json.loads(run_file.read_text()) if run_file.is_file() else {}
            active = active_task(directory)
            status = "running" if active else metadata["status"]
            if status == "running":
                status = "running" if active else "interrupted"
            missing = 0
            for planned in task["planned_jobs"]:
                relative = Path(task["directory"]) / "run/jobs" / planned["id"]
                result_file = output / relative / "result.json"
                if result_file.is_file():
                    result = json.loads(result_file.read_text())
                    if any(result.get(k) != v for k, v in planned.items()):
                        raise ValueError("result identity differs from plan: " + str(result_file))
                    row = {**result, "directory": str(relative), "result_present": True}
                else:
                    missing += 1
                    row = {**planned, "directory": str(relative), "result_present": False,
                           "status": "task_failed" if status == "failed" else
                                     "not_started" if status == "not_started" else "unfinished",
                           "returncode": None, "elapsed_seconds": None,
                           "time_wall_seconds": None, "peak_rss_kib": None}
                row.update(tier=task["tier"], task_id=task["task_id"],
                           global_task_number=task["global_task_number"])
                jobs.append(row)
            if not active and missing and status == "success":
                status = "interrupted"
            if not active and status == "success" and any(j["status"] != "success"
                    for j in jobs if j["global_task_number"] == task["global_task_number"]):
                status = "failed"
            tasks.append({**metadata, "status": status, "directory": task["directory"],
                          "missing_jobs": missing, "planned_jobs": len(task["planned_jobs"]),
                          "estimated_seconds": task["estimated_seconds"],
                          "estimated_memory_mib": task["estimated_memory_mib"],
                          "run_manifest_sha256": digest(run_file)})
            checks.append({k: manifest.get(k) for k in (
                "software_changed_during_run", "binaries_changed_during_run", "software_lock_valid_after_run")})
        counts = dict(Counter(t["status"] for t in tasks))
        finished = all(t["status"] in ("success", "failed") for t in tasks)
        status = "running" if "running" in counts else "failed" if finished and "failed" in counts else \
                 "success" if finished else "interrupted"

        def combined(key):
            values = [c[key] for c in checks]
            if key == "software_lock_valid_after_run":
                return False if False in values else True if all(v is True for v in values) else None
            return True if True in values else False if all(v is False for v in values) else None

        manifest = {"schema_version": 1, "run_kind": "array", "status": status,
                    "purpose": plan["config"]["purpose"], "config": plan["config"],
                    "array_config": plan["array_config"], "created_utc": plan["created_utc"],
                    "started_utc": min((a["started_utc"] for t in tasks for a in t["attempts"]), default=None),
                    "dependencies_verified_at_planning": plan["dependencies_verified_at_planning"],
                    "config_sha256": plan["original_config_sha256"], "inputs": plan["inputs"],
                    "software_before": plan["software_before"], "experiment_before": plan["experiment_before"],
                    "software_lock": plan["software_lock"], "array_plan_sha256": digest(plan_path),
                    "aggregated_utc": now(), "tasks": tasks, "task_status_counts": counts,
                    "planned_tasks": len(tasks), "jobs": jobs,
                    "planned_measured_jobs": sum(not j["warmup"] for j in jobs),
                    "completed_measured_jobs": sum(not j["warmup"] and j["result_present"] for j in jobs),
                    "missing_measured_jobs": sum(not j["warmup"] and not j["result_present"] for j in jobs),
                    **{k: combined(k) for k in checks[0]}}
        if finished:
            manifest["finished_utc"] = max(t.get("finished_utc") or "" for t in tasks)
        write_json(output / "run.json", manifest)
        return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True)
    args = parser.parse_args()
    manifest = aggregate(experiment_path(args.plan))
    print(json.dumps({k: manifest[k] for k in ("status", "planned_tasks", "task_status_counts",
                      "planned_measured_jobs", "completed_measured_jobs", "missing_measured_jobs")}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, OSError) as exc:
        print("error: " + str(exc), file=sys.stderr)
        sys.exit(2)
