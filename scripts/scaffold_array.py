#!/usr/bin/env python3
"""Plan, execute, or aggregate full bpRNA G arrays. Never submit jobs."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shlex
import shutil
import signal
import sys

from biological_cohort import LP_COMMIT, SETTINGS, acquired_records
from plan_array import exclusive, make_tasks
from prepare_biological_inputs import prepare as prepare_inputs
from run_cli import ROOT, digest, experiment_path, git_state, now, write_json


def verify_lp(lp, lock):
    state = git_state(lp.parent)
    if state["commit"] != lock["source_commit"] or state["dirty"]:
        raise ValueError("LinearPartition differs from clean pinned source")
    for key, path in (("wrapper_sha256", lp), ("vienna_binary_sha256", lp.parent / "bin/linearpartition_v")):
        if digest(path) != lock[key]:
            raise ValueError("LinearPartition " + key + " changed")
    for entry in lock.get("extra_files", []):
        if digest(experiment_path(entry["path"])) != entry["sha256"]:
            raise ValueError("scaffold preparation code changed: " + entry["path"])


def load(path):
    plan = json.loads(path.read_text())
    output = Path(plan["output_dir"])
    if path.resolve() != output / "array_plan.json" or digest(path) != (output / "array_plan.sha256").read_text().strip():
        raise ValueError("scaffold plan changed")
    for file, sha in ((Path(plan["source_path"]), plan["source_sha256"]),
                      (Path(plan["source_manifest"]), plan["source_manifest_sha256"]),
                      (output / "config.json", plan["config_sha256"]),
                      (output / "software_lock.json", plan["software_lock_sha256"])):
        if digest(file) != sha:
            raise ValueError("frozen scaffold dependency changed: " + str(file))
    return plan


def plan_array(config_path, output, jobs):
    config = json.loads(config_path.read_text())
    acquisition = experiment_path(config["source_manifest"])
    info, source, sequences = acquired_records(acquisition)
    lock = json.loads(experiment_path(config["software_lock"]).read_text())
    if lock["source_commit"] != LP_COMMIT:
        raise ValueError("unexpected LinearPartition revision")
    lp = experiment_path(config["linearpartition"])
    verify_lp(lp, lock)
    records = [{"name": name, "length": len(seq), "sequence_sha256": hashlib.sha256(seq.encode()).hexdigest()}
               for name, seq in sorted(sequences, key=lambda r: (len(r[1]), r[0]))]
    matrix = {"schema_version": 1, "purpose": "development", "inputs": records,
              "software_lock": config["software_lock"], "programs": [{"name": "LinearPartitionV"}],
              "modes": [{"name": "scaffold"}], "repetitions": 1, "warmups": 0,
              "memory_limit_mib": 8192, "timeout_seconds": config["timeout_seconds"]}
    tasks = make_tasks(matrix, config["array"], records)
    if output.exists() or jobs.exists():
        raise ValueError("use new immutable scaffold plan directories")
    output.mkdir(parents=True)
    jobs.mkdir(parents=True)
    (output / "scheduler_logs").mkdir()
    shutil.copyfile(config_path, output / "config.json")
    write_json(output / "software_lock.json", lock)
    plan = {"schema_version": 1, "run_kind": "scaffold_array", "created_utc": now(),
            "output_dir": str(output), "config": config, "config_sha256": digest(output / "config.json"),
            "source_manifest": str(acquisition), "source_manifest_sha256": digest(acquisition),
            "source_path": str(source), "source_sha256": info["sha256"],
            "software_lock": lock, "software_lock_sha256": digest(output / "software_lock.json"),
            "linearpartition": str(lp), "inputs": records, "tasks": tasks, "job_scripts": [],
            "dependencies_verified_at_planning": True, "scheduler_submitted": False}
    for task in tasks:
        task["directory"] = str(Path("tasks") / task["tier"] / f'{task["task_id"]:06d}')
        directory = output / task["directory"]
        directory.mkdir(parents=True)
        write_json(directory / "task.json", {"schema_version": 1, "status": "not_started", "attempts": []})
    for tier in config["array"]["tiers"]:
        count = sum(t["tier"] == tier["name"] for t in tasks)
        if not count:
            continue
        script = jobs / (tier["name"] + ".sh")
        script.write_text("\n".join(["#!/usr/bin/env bash", "#$ -cwd", "#$ -S /bin/bash", "#$ -j y",
            "#$ -N g_" + tier["name"], "#$ -l " + tier["resource"] + "=1", "#$ -l h_rt=" + tier["h_rt"],
            "#$ -t 1-" + str(count), "#$ -o " + str(output / "scheduler_logs"), "set -euo pipefail",
            ': "${SGE_TASK_ID:?Run as an array task}"', ': "${JOB_ID:?Missing allocation}"',
            "cd " + shlex.quote(str(ROOT)), "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1",
            'exec "${LCP_PYTHON:-python3}" scripts/scaffold_array.py run --plan ' +
            shlex.quote(str(output / "array_plan.json")) + " --tier " + tier["name"] + ' "$@"', ""]))
        script.chmod(0o755)
        plan["job_scripts"].append({"tier": tier["name"], "tasks": count, "path": str(script), "sha256": digest(script)})
    write_json(output / "array_plan.json", plan)
    (output / "array_plan.sha256").write_text(digest(output / "array_plan.json") + "\n")
    aggregate(output / "array_plan.json")
    return plan


def run(path, tier, task_id=None, resume=False):
    plan = load(path)
    identifier = os.environ.get("SGE_TASK_ID", str(task_id) if task_id is not None else "")
    if not re.fullmatch(r"[1-9][0-9]*", identifier) or (task_id is not None and int(identifier) != task_id):
        raise ValueError("invalid or conflicting SGE_TASK_ID")
    if plan["config"]["purpose"] == "formal" and not (os.environ.get("JOB_ID") and os.environ.get("SGE_TASK_ID")):
        raise ValueError("formal scaffold tasks require array allocation")
    task = next((t for t in plan["tasks"] if t["tier"] == tier and t["task_id"] == int(identifier)), None)
    if task is None:
        raise ValueError("task outside saved scaffold array")
    output = Path(plan["output_dir"]) / task["directory"]
    with exclusive(output / ".execution.lock"):
        metadata = json.loads((output / "task.json").read_text())
        if metadata["attempts"] and not resume:
            raise ValueError("task attempted: use --resume")
        lp = Path(plan["linearpartition"])
        verify_lp(lp, plan["software_lock"])
        attempt = {"started_utc": now(), "finished_utc": None, "hostname": platform.node(),
                   "JOB_ID": os.environ.get("JOB_ID"), "SGE_TASK_ID": os.environ.get("SGE_TASK_ID"),
                   "task_id": task["task_id"], "tier": tier, "resume": resume,
                   "requested_resource": task["requested_resource"], "h_rt": task["h_rt"]}
        metadata["attempts"].append(attempt)
        metadata["status"] = "running"
        write_json(output / "task.json", metadata)
        try:
            result = prepare_inputs(plan["source_manifest"], output / "run", lp, "bprna_full", 20261011,
                                    plan["config"]["timeout_seconds"], [i["name"] for i in task["inputs"]],
                                    resume=resume and (output / "run/dataset.json").exists())
            verify_lp(lp, plan["software_lock"])
            status = result["status"]
        except KeyboardInterrupt:
            status = "interrupted"
        except (ValueError, OSError) as exc:
            status = "failed"
            attempt["error"] = str(exc)
        attempt.update(status=status, finished_utc=now())
        metadata["status"] = status
        write_json(output / "task.json", metadata)
        return 0 if status == "success" else 130 if status == "interrupted" else 1


def aggregate(path, dataset_output=None):
    from aggregate_array import active_task
    plan = load(path)
    output = Path(plan["output_dir"])
    with exclusive(output / ".aggregate.lock"):
        records, statuses, datasets = [], [], []
        for task in plan["tasks"]:
            directory = output / task["directory"]
            meta = json.loads((directory / "task.json").read_text())
            if active_task(directory):
                raise ValueError("wait for active scaffold tasks before aggregation")
            saved = directory / "run/dataset.json"
            dataset = json.loads(saved.read_text()) if saved.exists() else None
            found = {r["id"]: r for r in dataset["records"]} if dataset else {}
            for item in task["inputs"]:
                record = found.get(item["name"], {"id": item["name"], "length": item["length"],
                                                 "scaffold_status": "not_started" if meta["status"] == "not_started" else "unfinished"})
                if record.get("sequence_sha256", item["sequence_sha256"]) != item["sequence_sha256"]:
                    raise ValueError("scaffold sequence differs from plan")
                if record["scaffold_status"] == "success":
                    file = directory / "run" / record["input_file"]
                    if digest(file) != record["input_sha256"]:
                        raise ValueError("generated G changed")
                    record = {**record, "task_input_path": str(file), "task_dataset_sha256": digest(saved)}
                records.append(record)
            statuses.append(meta["status"])
            datasets.append(dataset)
        success = all(r["scaffold_status"] == "success" for r in records) and all(s == "success" for s in statuses)
        summary = {"schema_version": 1, "status": "success" if success else "incomplete",
                   "planned_tasks": len(plan["tasks"]), "planned_jobs": len(records),
                   "completed_jobs": sum(r["scaffold_status"] in ("success", "failed") for r in records),
                   "array_plan_sha256": digest(path), "records": records}
        write_json(output / "summary.json", summary)
        if dataset_output is not None:
            if not success:
                raise ValueError("G incomplete: retain every missing/failed input; dataset not published")
            if any(d.get("software_unchanged") is not True or d["source_sha256"] != plan["source_sha256"] or
                   d["scaffold_settings"] != SETTINGS or d["linearpartition_source"]["commit"] != LP_COMMIT
                   for d in datasets):
                raise ValueError("task G provenance changed")
            verify_lp(Path(plan["linearpartition"]), plan["software_lock"])
            dataset_output.mkdir(parents=True, exist_ok=False)
            for record in records:
                shutil.copyfile(record["task_input_path"], dataset_output / record["input_file"])
            manifest = {**datasets[0], "kind": "bprna_full", "records": records, "status": "success",
                        "software_unchanged": True, "source_sha256": plan["source_sha256"],
                        "acquisition_manifest_sha256": plan["source_manifest_sha256"],
                        "selected_ids_before_scaffold": [r["id"] for r in records],
                        "cohort_ids_before_subset": [r["id"] for r in records], "requested_record_ids": None,
                        "array_plan_sha256": digest(path), "aggregated_utc": now()}
            write_json(dataset_output / "dataset.json", manifest)
        return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action", required=True)
    p = sub.add_parser("plan")
    p.add_argument("--config", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--jobs-dir")
    p = sub.add_parser("run")
    p.add_argument("--plan", required=True)
    p.add_argument("--tier", required=True)
    p.add_argument("--task-id", type=int)
    p.add_argument("--resume", action="store_true")
    p = sub.add_parser("aggregate")
    p.add_argument("--plan", required=True)
    p.add_argument("--dataset-out")
    args = parser.parse_args()
    if args.action == "plan":
        output = experiment_path(args.out)
        plan = plan_array(experiment_path(args.config), output,
                          experiment_path(args.jobs_dir) if args.jobs_dir else ROOT / "jobs" / output.name)
        print(json.dumps({"inputs": len(plan["inputs"]), "tasks": len(plan["tasks"]), "scripts": plan["job_scripts"]}))
    elif args.action == "run":
        return run(experiment_path(args.plan), args.tier, args.task_id, args.resume)
    else:
        summary = aggregate(experiment_path(args.plan), experiment_path(args.dataset_out) if args.dataset_out else None)
        print(json.dumps({k: summary[k] for k in ("status", "planned_tasks", "planned_jobs", "completed_jobs")}))
    return 0


if __name__ == "__main__":
    signal.signal(signal.SIGTERM, lambda *args: (_ for _ in ()).throw(KeyboardInterrupt()))
    try:
        sys.exit(main())
    except (ValueError, OSError, KeyError) as exc:
        print("error: " + str(exc), file=sys.stderr)
        sys.exit(2)
