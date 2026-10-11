#!/usr/bin/env python3
"""Freeze matched sequence batches and write one TSUBAME array script per tier.

This command never submits jobs. --dry-run permits offline software-lock checks;
without it, the source/build/parameters/binaries must match the saved lock.
"""

import argparse
from contextlib import contextmanager
import copy
import fcntl
import hashlib
import json
import math
import re
import shlex
import shutil
import sys
import io
import tarfile
from pathlib import Path

from rna_records import read_record
from task_storage import FORMAT, ArchiveReader, MemoryFile
from file_budget import check_budget, estimate_files
from run_cli import (ROOT, digest, experiment_path, git_state, now,
                     mode_applies, plan_jobs, validate, verify_lock, write_json)


# Official handbook, section 5.1.1; checked 2026-10-11. GB is decimal here,
# so the planner conservatively converts GB to MiB, not GB to GiB.
RESOURCE_SOURCE = "https://helpdesk.t3.gsic.titech.ac.jp/t4docs_preview/master/handbook.ja/jobs/"
RESOURCES = {
    "node_f": (192, 768), "node_h": (96, 384), "node_q": (48, 192),
    "node_o": (24, 96), "gpu_1": (8, 96), "gpu_h": (4, 48),
    "cpu_160": (160, 368), "cpu_80": (80, 184), "cpu_40": (40, 92),
    "cpu_16": (16, 36.8), "cpu_8": (8, 18.4), "cpu_4": (4, 9.2),
}


@contextmanager
def exclusive(path):
    """Kernel-released lock: a killed task leaves no stale ownership marker."""
    with path.open("a+") as stream:
        try:
            fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError("another process holds " + str(path)) from None
        try:
            yield
        finally:
            fcntl.flock(stream, fcntl.LOCK_UN)


def positive(value, label, zero=False):
    if type(value) not in (float, int) or not math.isfinite(value) or value < 0 or (not zero and value == 0):
        raise ValueError(label + " must be finite and " + ("nonnegative" if zero else "positive"))
    return value


def wall_seconds(value):
    if not isinstance(value, str) or not re.fullmatch(r"\d{1,2}:\d{2}:\d{2}", value):
        raise ValueError("h_rt must be HH:MM:SS")
    hours, minutes, seconds = map(int, value.split(":"))
    total = hours * 3600 + minutes * 60 + seconds
    if minutes >= 60 or seconds >= 60 or not 0 < total <= 86400:
        raise ValueError("h_rt must be positive and at most 24 hours")
    return total


def validate_array(config, settings):
    validate(config, array=True)
    if settings.get("schema_version") != 1 or settings.get("measurement") not in ("timing", "accuracy"):
        raise ValueError("array schema_version=1 and measurement=timing/accuracy are required")
    if type(settings.get("seed")) is not int:
        raise ValueError("array seed must be an integer")
    if not config.get("software_lock"):
        raise ValueError("array runs require a software lock, including development arrays")
    if "memory_limit_mib" not in config:
        raise ValueError("array runs require the existing per-job memory_limit_mib")
    for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        if config.get("environment", {}).get(key, "1") != "1":
            raise ValueError("array CLI processes require single-thread settings")
    tiers = settings.get("tiers", [])
    names, previous_max = set(), 0
    if not tiers:
        raise ValueError("array needs tiers")
    for tier in tiers:
        effective_config = {**config, **tier.get("job_limits", {})}
        if set(tier.get("job_limits", {})) - {"timeout_seconds", "memory_limit_mib"}:
            raise ValueError("job_limits supports timeout_seconds and memory_limit_mib only")
        validate(effective_config, array=True)
        name = tier.get("name", "")
        if not re.fullmatch(r"[A-Za-z0-9_-]+", name) or name in names:
            raise ValueError("tiers need unique simple names")
        names.add(name)
        low, high = tier.get("min_length"), tier.get("max_length")
        if type(low) is not int or type(high) is not int or low < 1 or high < low or low <= previous_max:
            raise ValueError("tiers must be sorted nonoverlapping inclusive length intervals")
        previous_max = high
        workers = tier.get("workers", 1)
        resource = tier.get("resource")
        if resource not in RESOURCES or type(workers) is not int or not 1 <= workers <= RESOURCES[resource][0]:
            raise ValueError("unknown resource or workers exceed allocated cores")
        if settings["measurement"] == "timing" and (resource != "node_f" or workers != 1):
            raise ValueError("timing arrays require node_f=1 and workers=1")
        wall = wall_seconds(tier.get("h_rt"))
        target = positive(tier.get("target_seconds"), "target_seconds")
        fraction = positive(tier.get("max_wall_fraction", 0.7), "max_wall_fraction")
        if fraction > 0.8 or target > wall * fraction:
            raise ValueError("target_seconds must fit max_wall_fraction (at most 0.8) of h_rt")
        memory = positive(tier.get("memory_budget_mib"), "memory_budget_mib")
        if memory > RESOURCES[resource][1] * 1e9 / 2**20:
            raise ValueError("memory_budget_mib exceeds the documented resource memory")
        overhead = positive(tier.get("task_overhead_mib", 256), "task_overhead_mib", zero=True)
        if workers * effective_config["memory_limit_mib"] + overhead > memory:
            raise ValueError("workers times per-job memory cap plus overhead exceeds memory budget")
        model = tier.get("cost_model", {})
        if not model.get("basis") or not isinstance(model["basis"], str):
            raise ValueError("cost_model.basis must explain the estimates")
        for key in ("seconds_per_nt_condition", "mib_per_nt"):
            positive(model.get(key), "cost_model." + key, zero=key == "mib_per_nt")
        for key, default in (("seconds_constant_per_condition", 0), ("memory_base_mib", 0)):
            positive(model.get(key, default), "cost_model." + key, zero=True)
        positive(model.get("length_exponent", 1), "length_exponent")
        if positive(model.get("safety_factor", 1), "safety_factor") < 1:
            raise ValueError("cost_model.safety_factor must be at least 1")
        for mode, factor in model.get("mode_factors", {}).items():
            if mode not in {m["name"] for m in config["modes"]}:
                raise ValueError("unknown mode in cost model")
            positive(factor, "mode factor")
        for mode, fit in model.get("per_mode", {}).items():
            if mode not in {m["name"] for m in config["modes"]}:
                raise ValueError("unknown per_mode cost")
            for key in ("constant_seconds", "seconds_coefficient", "memory_base_mib", "memory_coefficient"):
                positive(fit.get(key, 0), key, zero=True)
            positive(fit.get("length_exponent", 1), "length_exponent")
            positive(fit.get("memory_exponent", 1), "memory_exponent")


def estimate(length, config, tier):
    model = tier["cost_model"]
    conditions = len(config["programs"]) * (config.get("warmups", 0) + config.get("repetitions", 1))
    modes = [m for m in config["modes"] if mode_applies(m, {"length": length})]
    weighted = sum(model.get("mode_factors", {}).get(mode["name"], 1) for mode in modes)
    single = (model.get("seconds_constant_per_condition", 0)
              + model["seconds_per_nt_condition"] * length**model.get("length_exponent", 1))
    seconds = model.get("safety_factor", 1) * single * conditions * weighted
    memory = model.get("memory_base_mib", 0) + model["mib_per_nt"] * length
    if "per_mode" in model:
        seconds, memory = 0, 0
        for mode in modes:
            fit = model["per_mode"][mode["name"]]
            seconds += (fit.get("constant_seconds", 0) + fit["seconds_coefficient"] *
                        length ** fit.get("length_exponent", 1)) * conditions
            memory = max(memory, fit.get("memory_base_mib", 0) + fit["memory_coefficient"] *
                         length ** fit.get("memory_exponent", 1))
        seconds *= model.get("safety_factor", 1)
    if memory > config["memory_limit_mib"]:
        raise ValueError("estimated job memory exceeds its unchanged per-job cap")
    return {"estimated_seconds": seconds, "estimated_job_memory_mib": memory}


def make_tasks(config, settings, inputs):
    """Greedy batches by estimated work, without splitting any RNA's matrix.

    The summed sequential work is a conservative wall estimate even for accuracy
    workers. It deliberately does not assume ideal parallel speedup.
    """
    validate_array(config, settings)
    tasks = []
    for item in inputs:
        matches = [t for t in settings["tiers"] if t["min_length"] <= item["length"] <= t["max_length"]]
        if len(matches) != 1:
            raise ValueError("no unique tier for input " + item["name"])
    for tier in settings["tiers"]:
        effective_config = {**config, **tier.get("job_limits", {})}
        batches, batch = [], []

        def cost(items):
            return (sum(i["estimated_seconds"] for i in items),
                    tier.get("task_overhead_mib", 256) + max(
                        (i["estimated_job_memory_mib"] for i in items), default=0) * tier.get("workers", 1))

        groups = {}
        for item in inputs:
            if not tier["min_length"] <= item["length"] <= tier["max_length"]:
                continue
            # Multiple supplied scaffolds for the same RNA must stay together.
            groups.setdefault(item.get("sequence_sha256", item["name"]), []).append(
                {**item, **estimate(item["length"], effective_config, tier)})
        for group in groups.values():
            seconds, memory = cost(group)
            if seconds > tier["target_seconds"] or memory > tier["memory_budget_mib"]:
                raise ValueError("one complete input matrix exceeds tier budget: " + group[0]["name"])
            seconds, memory = cost(batch + group)
            if batch and (seconds > tier["target_seconds"] or memory > tier["memory_budget_mib"] or
                          len(batch + group) > tier.get("max_inputs_per_task", float("inf"))):
                batches.append(batch)
                batch = []
            batch.extend(group)
        if batch:
            batches.append(batch)
        for index, batch in enumerate(batches, 1):
            seconds, memory = cost(batch)
            tasks.append({"tier": tier["name"], "task_id": index,
                          "global_task_number": len(tasks) + 1, "inputs": batch,
                          "estimated_seconds": seconds, "estimated_memory_mib": memory,
                          "memory_cap_envelope_mib": tier.get("workers", 1) * effective_config["memory_limit_mib"]
                                                     + tier.get("task_overhead_mib", 256),
                          "job_limits": tier.get("job_limits", {}),
                          "requested_resource": tier["resource"] + "=1", "h_rt": tier["h_rt"],
                          "workers": tier.get("workers", 1), "cost_model": tier["cost_model"]})
    return tasks


def check_lock(config):
    path = experiment_path(config["software_lock"])
    lock = json.loads(path.read_text())
    hashes = [lock.get("cmake_cache_sha256"), lock.get("energy_parameter_sha256")]
    hashes += [lock.get("binaries", {}).get(p["name"]) for p in config["programs"]]
    if (lock.get("schema_version") != 1 or not re.fullmatch(r"[0-9a-f]{40}", str(lock.get("source_commit")))
            or any(not isinstance(h, str) or not re.fullmatch(r"[0-9a-f]{64}", h) for h in hashes)):
        raise ValueError("software lock needs a source commit, cache, parameter and binary hashes")
    if config["purpose"] == "formal" and lock.get("build_revision_verified") is not True:
        raise ValueError("formal lock must record build_revision_verified=true")
    return path, lock


def script_text(plan, tier, count):
    output = Path(plan["output_dir"])
    script = ["#!/usr/bin/env bash", "#$ -cwd", "#$ -S /bin/bash", "#$ -j y",
              "#$ -N lcp_" + tier["name"], "#$ -l " + tier["resource"] + "=1",
              "#$ -l h_rt=" + tier["h_rt"], "#$ -t 1-" + str(count),
              "#$ -o " + str(output / "scheduler_logs"), "", "set -euo pipefail",
              ': "${SGE_TASK_ID:?Run as an SGE array task}"',
              ': "${JOB_ID:?Missing scheduler allocation}"',
              "cd " + shlex.quote(str(ROOT)),
              "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1",
              '# Optional first argument: --resume (submit only the interrupted task range).',
              'exec "${LCP_PYTHON:-python3}" scripts/run_array_task.py \\',
              "  --plan " + shlex.quote(str(output / "array_plan.json")) + " \\",
              "  --tier " + shlex.quote(tier["name"]) + ' "$@"', ""]
    return "\n".join(script)


def prepare(config_path, settings_path, output, jobs_dir, software=None, build=None, dry_run=False, estimate_only=False):
    config = json.loads(config_path.read_text())
    from biological_cohort import expand_cohort
    config, cohort_provenance, pending = expand_cohort(config, dry_run)
    pending |= any(i.get("scaffold_pending", False) for i in config["inputs"])
    if pending and not dry_run:
        raise ValueError("pending scaffolds permit only a dry-run")
    settings = json.loads(settings_path.read_text()) if settings_path else config.get("array", {})
    validate_array(config, settings)
    lock_path, lock = check_lock(config)
    if bool(software) != bool(build):
        raise ValueError("supply both software-root and build-dir")
    if not dry_run and software is None:
        raise ValueError("software-root/build-dir are required unless --dry-run")
    state, verified = lock.get("source_git_state", {"commit": lock["source_commit"]}), False
    if software is not None:
        verify_lock(lock_path, software, build, config["programs"], software / config["energy_parameter_file"])
        state, verified = git_state(software), True
    inputs, existing_archives = [], {}
    for entry in config["inputs"]:
        source = experiment_path(entry["path"])
        if entry.get("scaffold_pending"):
            sequence = entry["planning_sequence"]
            payload = ">" + entry["name"] + "\n" + sequence + "\n" + "." * len(sequence) + "\n"
            sha = hashlib.sha256(payload.encode()).hexdigest()
        else:
            if "input_archive" in entry:
                archive_path = experiment_path(entry["input_archive"])
                if archive_path not in existing_archives:
                    existing_archives[archive_path] = ArchiveReader(archive_path)
                data = existing_archives[archive_path].read(entry["input_member"])
                _, sequence, _ = read_record(MemoryFile(data))
                sha = hashlib.sha256(data).hexdigest()
            else:
                _, sequence, _ = read_record(source)
                sha = digest(source)
        inputs.append({**entry, "source_path": str(source), "length": len(sequence), "sha256": sha,
                       "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest()})
    tasks = make_tasks(config, settings, inputs)
    nonempty = len({t["tier"] for t in tasks})
    budget = check_budget(estimate_files(len(tasks), nonempty, settings.get("compact_plan", False)), settings)
    if estimate_only:
        if not dry_run:
            raise ValueError("--estimate-only requires --dry-run")
        if output.exists():
            raise ValueError("use new immutable estimate output directory")
        planned_jobs = sum(len(config["programs"]) * (config.get("warmups", 0) + config.get("repetitions", 1)) *
                           sum(mode_applies(m, i) for m in config["modes"]) for i in inputs)
        report = {"schema_version": 1, "run_kind": "array_estimate", "output_format": FORMAT,
                  "dry_run": True, "executable": False, "created_utc": now(),
                  "config_sha256": digest(config_path), "array_config_sha256": digest(settings_path) if settings_path else None,
                  "dependencies_verified_at_planning": verified, "pending_scaffolds": pending,
                  "inputs": len(inputs), "tasks": len(tasks), "planned_jobs": planned_jobs,
                  "file_budget": budget, "estimate_artifact_entries": 2,
                  "tiers": [{"tier": name, "tasks": sum(t["tier"] == name for t in tasks),
                             "inputs": sum(len(t["inputs"]) for t in tasks if t["tier"] == name)}
                            for name in sorted({t["tier"] for t in tasks})],
                  "task_assignments": [{"tier": t["tier"], "task_id": t["task_id"],
                                        "inputs": [i["name"] for i in t["inputs"]]} for t in tasks]}
        output.mkdir(parents=True)
        write_json(output / "array_estimate.json", report)
        return report
    if output.exists() or jobs_dir.exists():
        raise ValueError("use new immutable output and jobs directories")
    if any(c.isspace() for c in str(output)):
        raise ValueError("scheduler output directory cannot contain whitespace")
    output.mkdir(parents=True)
    jobs_dir.mkdir(parents=True)
    (output / "scheduler_logs").mkdir()
    shutil.copyfile(config_path, output / "original_config.json")
    shutil.copyfile(lock_path, output / "software_lock.json")
    write_json(output / "array_config.json", settings)
    source_archives = {}
    input_tar = output / "inputs.tar"
    with tarfile.open(input_tar, "w") as archive:
        for entry in inputs:
            if entry.get("scaffold_pending"):
                sequence = entry.pop("planning_sequence")
                data = (">" + entry["name"] + "\n" + sequence + "\n" + "." * len(sequence) + "\n").encode()
            elif "input_archive" in entry:
                path = experiment_path(entry["input_archive"])
                if path not in source_archives:
                    source_archives[path] = ArchiveReader(path)
                data = source_archives[path].read(entry["input_member"])
            else:
                data = Path(entry["source_path"]).read_bytes()
            if hashlib.sha256(data).hexdigest() != entry["sha256"]:
                raise ValueError("input changed during planning")
            member = "inputs/" + entry["name"] + ".txt"
            info = tarfile.TarInfo(member)
            info.size = len(data)
            archive.addfile(info, io.BytesIO(data))
            entry.update(input_archive=str(input_tar), input_member=member,
                         copied_path=str(output / member))
    for entry in config["inputs"]:
        entry.pop("planning_sequence", None)
    frozen = {i["name"]: i for i in inputs}
    config_entries = {i["name"]: i for i in config["inputs"]}
    for task in tasks:
        for item in task["inputs"]:
            item.pop("planning_sequence", None)
        task["inputs"] = [{**i, "copied_path": frozen[i["name"]]["copied_path"]} for i in task["inputs"]]
        directory = Path("tasks") / task["tier"] / f'{task["task_id"]:06d}'
        (output / directory).mkdir(parents=True)
        task_config = copy.deepcopy(config)
        task_config.update(task["job_limits"])
        task_config.pop("array", None)
        task_config.pop("input_cohort", None)
        task_config["inputs"] = [{**config_entries[i["name"]], "length": i["length"],
                                  "path": frozen[i["name"]]["copied_path"],
                                  **{k: frozen[i["name"]][k] for k in ("input_archive", "input_member", "sha256")}}
                                 for i in task["inputs"]]
        task_config["software_lock"] = str(output / "software_lock.json")
        tier = next(t for t in settings["tiers"] if t["name"] == task["tier"])
        task_config["order"] = {"kind": tier.get("order_kind", "blocked_random"), "seed": settings["seed"]}
        task_config.setdefault("environment", {}).update(
            OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        config_file = output / directory / "task_config.json"
        write_json(config_file, task_config)
        task.update(directory=str(directory), config_path=str(config_file),
                    config_sha256=digest(config_file), planned_jobs=plan_jobs(task_config))
        if settings.get("compact_plan"):
            schedule = output / directory / "planned_jobs.json"
            write_json(schedule, task.pop("planned_jobs"))
            task.update(planned_jobs_path=str(schedule), planned_jobs_sha256=digest(schedule))
        write_json(output / directory / "task.json", {"schema_version": 1, "status": "not_started",
                   "tier": task["tier"], "task_id": task["task_id"], "attempts": []})
    plan = {"schema_version": 1, "run_kind": "array", "created_utc": now(),
            "output_format": FORMAT, "file_budget": budget, "inputs_archive_sha256": digest(input_tar),
            "output_dir": str(output), "config": config, "array_config": settings,
            "original_config_sha256": digest(output / "original_config.json"),
            "array_config_sha256": digest(output / "array_config.json"),
            "software_lock": lock, "software_lock_sha256": digest(output / "software_lock.json"),
            "software_root": str(software) if software else None, "build_dir": str(build) if build else None,
            "software_before": state, "dependencies_verified_at_planning": verified,
            "dry_run": dry_run, "experiment_before": git_state(ROOT), "inputs": inputs, "tasks": tasks,
            "pending_scaffolds": pending, "cohort_provenance": cohort_provenance,
            "resource_reference": {"url": RESOURCE_SOURCE, "checked_date": "2026-10-11",
                                   "memory_unit": "decimal GB; conservative conversion to MiB"}}
    scripts = []
    for tier in settings["tiers"]:
        count = sum(t["tier"] == tier["name"] for t in tasks)
        if not count:
            continue
        path = jobs_dir / (tier["name"] + ".sh")
        path.write_text(script_text(plan, tier, count))
        path.chmod(0o755)
        scripts.append({"tier": tier["name"], "tasks": count, "path": str(path), "sha256": digest(path)})
    plan["job_scripts"] = scripts
    write_json(output / "array_plan.json", plan)
    (output / "array_plan.sha256").write_text(digest(output / "array_plan.json") + "\n")
    from aggregate_array import aggregate
    aggregate(output / "array_plan.json", summary_only=settings.get("compact_plan", False))
    return plan


def load_plan(path, task_selector=None, verify_frozen_inputs=True):
    plan = json.loads(path.read_text())
    output = Path(plan["output_dir"])
    if path.resolve() != output / "array_plan.json" or plan.get("schema_version") != 1:
        raise ValueError("unexpected plan path or schema")
    if digest(path) != (output / "array_plan.sha256").read_text().strip():
        raise ValueError("array plan changed")
    checks = [(output / "original_config.json", plan["original_config_sha256"]),
              (output / "array_config.json", plan["array_config_sha256"]),
              (output / "software_lock.json", plan["software_lock_sha256"])]
    selected = [t for t in plan["tasks"] if task_selector is None or (t["tier"],t["task_id"]) == task_selector]
    if task_selector is not None and len(selected) != 1:
        raise ValueError("task is outside the saved tier array")
    inputs = plan["inputs"] if task_selector is None else selected[0]["inputs"]
    if verify_frozen_inputs:
        if "inputs_archive_sha256" in plan:
            checks += [(output / "inputs.tar", plan["inputs_archive_sha256"])]
        else:
            checks += [(Path(i["copied_path"]), i["sha256"]) for i in inputs]
    checks += [(Path(t["config_path"]), t["config_sha256"]) for t in selected]
    checks += [(Path(t["planned_jobs_path"]), t["planned_jobs_sha256"]) for t in selected if "planned_jobs_path" in t]
    for file, expected in checks:
        if digest(file) != expected:
            raise ValueError("frozen dependency changed: " + str(file))
    return plan


def task_jobs(plan, task):
    if "planned_jobs" in task:
        return task["planned_jobs"]
    path = Path(task["planned_jobs_path"])
    if digest(path) != task["planned_jobs_sha256"]:
        raise ValueError("task schedule changed")
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--array-config", help="array settings overlay; otherwise read config.array")
    parser.add_argument("--out", required=True)
    parser.add_argument("--jobs-dir", help="default: jobs/<output basename>")
    parser.add_argument("--software-root", type=Path)
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--estimate-only", action="store_true", help="dry-run counts and assignments in one JSON; no task files or input copies")
    args = parser.parse_args()
    output = experiment_path(args.out)
    plan = prepare(experiment_path(args.config), experiment_path(args.array_config) if args.array_config else None,
                   output, experiment_path(args.jobs_dir) if args.jobs_dir else ROOT / "jobs" / output.name,
                   args.software_root.resolve() if args.software_root else None,
                   args.build_dir.resolve() if args.build_dir else None, args.dry_run, args.estimate_only)
    if args.estimate_only:
        print(json.dumps({k: plan[k] for k in ("inputs", "tasks", "planned_jobs", "file_budget", "pending_scaffolds")}, indent=2))
        return
    print(json.dumps({"plan": str(output / "array_plan.json"), "dry_run": plan["dry_run"],
                      "dependencies_verified": plan["dependencies_verified_at_planning"],
                      "inputs": len(plan["inputs"]), "tasks": len(plan["tasks"]),
                      "job_scripts": plan["job_scripts"], "file_budget": plan["file_budget"]}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, OSError) as exc:
        print("error: " + str(exc), file=sys.stderr)
        sys.exit(2)
