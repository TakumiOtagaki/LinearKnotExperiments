#!/usr/bin/env python3
"""Prepare the E1 long cpu_16 replacement from frozen inputs and measured CLI times.

Never submits jobs. Requires completed development probes; keeps unsuccessful
attempts in the pilot receipt/archive. Existing output/config paths are refused.
"""
import argparse
from collections import Counter
import copy
import hashlib
import json
import os
from pathlib import Path
import tarfile

from file_budget import check_budget, estimate_files
from plan_array import load_plan, make_tasks, prepare, task_jobs, wall_seconds
from run_cli import ROOT, digest, now, write_json

BASE = ROOT / "results/runs/e1-bprna-long-89baea52-20261018"
HISTORY = ROOT / "results/runs/e1-bprna-long-d21618c5-20261014"
DOCS = ROOT / "docs/rerun_89baea52_20261018"
NAME = "e1-bprna-long-89baea52-cpu16-20261018"


def snapshot(paths):
    files = sorted({p for root in paths for p in ([root] if root.is_file() else root.rglob("*")) if p.is_file()})
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in files}
    return {"files": len(files), "manifest_sha256": hashlib.sha256(
        json.dumps(hashes, sort_keys=True).encode()).hexdigest()}


def read_history():
    rows, sources = [], {}
    for path in sorted(HISTORY.glob("tasks/*/*/run/task_results.jsonl")):
        data = path.read_bytes()
        sources[str(path.relative_to(ROOT))] = hashlib.sha256(data).hexdigest()
        for line in data.splitlines():
            record = json.loads(line)
            rows.append({k: record.get(k) for k in (
                "id", "input", "mode", "program", "warmup", "repetition", "status",
                "elapsed_seconds", "time_wall_seconds", "peak_rss_kib")})
    return {"rows": rows, "source_sha256": sources}


def calibrated_settings(base, history, pilots):
    rows = history["rows"]
    expected = {job["id"] for task in base["tasks"] for job in task_jobs(base, task)}
    if len(rows) != len(expected) or {r["id"] for r in rows} != expected:
        raise ValueError("historical records must cover the exact frozen matrix once")
    if any(r["status"] != "success" or not r["elapsed_seconds"] or r["elapsed_seconds"] <= 0 for r in rows):
        raise ValueError("historical failure/missing time needs an explicit conservative model")
    by_key = {(r["input"], r["mode"], r["program"]): r for r in rows}
    if len(by_key) != len(rows):
        raise ValueError("historical matrix must have one measurement per condition")
    completed = [p for p in pilots if p["status"] == "success"]
    if not 1 <= len({p["input"] for p in completed}) <= 2 or len({p["input"] for p in pilots}) > 2:
        raise ValueError("require a successful probe of at least one of at most two long inputs")
    inputs = {i["name"]: i for i in base["inputs"]}
    ratios = []
    for pilot in completed:
        if (pilot["memory_limit_mib"] != 4096 or pilot["input_sha256"] != inputs[pilot["input"]]["sha256"]
                or pilot["binary_sha256"] != base["software_lock"]["binaries"][pilot["program"]]
                or not pilot.get("peak_rss_kib")):
            raise ValueError("pilot inputs/build/cap differ from the frozen plan")
        old = by_key[pilot["input"], pilot["mode"], pilot["program"]]
        ratios.append({"input": pilot["input"], "program": pilot["program"], "mode": pilot["mode"],
                       "old_elapsed_seconds": old["elapsed_seconds"],
                       "new_elapsed_seconds": pilot["elapsed_seconds"],
                       "new_over_old": pilot["elapsed_seconds"] / old["elapsed_seconds"]})
    # One or two examples cannot establish a cohort-wide speedup. Credit at most 2x,
    # transfer the slowest observed ratio, then apply 1.5 work safety below.
    release_factor = max(.5, max(r["new_over_old"] for r in ratios))
    censored = [p for p in pilots if p["status"] == "timeout"]
    # A timeout is a lower bound, not a completed speed measurement. Never let
    # the final per-job estimate fall below that observed lower bound.
    for pilot in censored:
        old = by_key[pilot["input"], pilot["mode"], pilot["program"]]
        release_factor = max(release_factor, pilot["elapsed_seconds"] / old["elapsed_seconds"] / 1.5)
    settings = copy.deepcopy(base["array_config"])
    for tier, hours in zip(settings["tiers"], (1, 2, 3, 4)):
        observed = {}
        for item in base["inputs"]:
            if not tier["min_length"] <= item["length"] <= tier["max_length"]:
                continue
            sums, maxima = {}, {}
            for mode in base["config"]["modes"]:
                # Includes launch, calculation and output; 0.5 s/job additionally
                # allows controller/checkpoint work before worker-efficiency loss.
                values = [by_key[item["name"], mode["name"], p["name"]]["elapsed_seconds"] * release_factor + .5
                          for p in base["config"]["programs"]]
                sums[mode["name"]], maxima[mode["name"]] = sum(values), max(values)
            observed[item["name"]] = {"sha256": item["sha256"], "seconds_by_mode": sums,
                                      "max_seconds_by_mode": maxima, "memory_mib": 4096}
        tier.update(resource="cpu_16", workers=8, memory_budget_mib=35095, task_overhead_mib=1024,
                    h_rt=f"{hours:02d}:00:00", target_seconds=hours * 3600 * .7, max_wall_fraction=.7,
                    job_limits={"memory_limit_mib": 4096, "timeout_seconds": base["config"]["timeout_seconds"]},
                    cost_model={"basis": "Exact d21618c5 task_results.jsonl elapsed_seconds matched by input/mode/program; "
                                "release factor and paired 89baea52 probes in docs/rerun_89baea52_20261018/cpu16/calibration.json; "
                                "0.5 s/job overhead, safety 1.5; scheduling estimate, not isolated performance evidence",
                                "seconds_per_nt_condition": 1, "mib_per_nt": 0,
                                "safety_factor": 1.5, "observed_inputs": observed},
                    parallel_packing={"kind": "beam_list_bound", "efficiency": .8, "overhead_seconds": 120,
                                      "basis": "Per beam W/8 + 7L/8 upper bound for work-conserving seeded lists; "
                                      "use L if at most 8 jobs in that beam; "
                                      "sum all six beam barriers, divide by assumed 0.8 worker efficiency, add 120 s/task"})
    settings["note"] = "E1 long 89baea52 cpu_16 replacement; unchanged scientific matrix, seed and ascending-beam policy. Not submitted."
    calibration = {"created_utc": now(), "historical_run": str(HISTORY.relative_to(ROOT)),
                   "historical_status_counts": dict(Counter(r["status"] for r in rows)),
                   "historical_cli_jobs": len(rows), "historical_cli_hours": sum(r["elapsed_seconds"] for r in rows)/3600,
                   "historical_jsonl_sha256": history["source_sha256"], "paired_probes": ratios,
                   "completed_probe_inputs": len({p["input"] for p in completed}),
                   "pilot_attempts": pilots, "release_factor": release_factor,
                   "release_factor_rule": "max(0.5, max(completed_probe_new_over_old)); credit at most a 2x speedup, then apply work safety 1.5",
                   "censored_probe_rule": "Also raise the factor if a timed-out attempt's elapsed lower bound exceeds its predicted job time",
                   "job_extra_seconds": .5, "work_safety_factor": 1.5, "worker_efficiency_assumption": .8,
                   "task_overhead_seconds": 120, "target_wall_fraction": .7,
                   "scope": "Planning from whole CLI times and one or two completed development probes, retaining timeouts. At most two probe CLIs overlapped. "
                            "No cpu_16 throughput benchmark; no claim of universal memory or runtime bounds."}
    return settings, calibration


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-results", nargs="+", required=True, type=Path)
    parser.add_argument("--history-cache", type=Path, help="optional node-local cache produced by the same historical scan")
    parser.add_argument("--name", default=NAME)
    args = parser.parse_args()
    if not os.environ.get("TMPDIR") or not Path(os.environ["TMPDIR"]).is_dir():
        raise ValueError("node-local TMPDIR is required")
    base = load_plan(BASE / "array_plan.json")
    old_jobs = ROOT / "jobs" / BASE.name
    before = snapshot([BASE, old_jobs])
    history = json.loads(args.history_cache.read_text()) if args.history_cache else read_history()
    pilots = [json.loads(p.read_text()) for p in args.pilot_results]
    settings, calibration = calibrated_settings(base, history, pilots)
    config = json.loads((BASE / "original_config.json").read_text())
    config.update(name=args.name, memory_limit_mib=4096)
    preflight_config = {**base["config"], "name": args.name, "memory_limit_mib": 4096}
    tasks = make_tasks(preflight_config, settings, base["inputs"])
    budget = check_budget(estimate_files(len(tasks), 4, True), settings)
    # Extra entries cover configs, receipts, pilot tar, docs and validation log.
    extra_budget = check_budget({**budget, "estimated_shared_entries": budget["estimated_shared_entries"] + 50}, settings)
    destination = DOCS / "cpu16"
    config_path = ROOT / "configs/rerun_89baea52_20261018" / (args.name + ".json")
    array_path = config_path.parent / "arrays" / config_path.name
    output = ROOT / "results/runs" / args.name
    jobs = ROOT / "jobs" / args.name
    for path in (destination, config_path, array_path, output, jobs):
        if path.exists():
            raise ValueError("refusing existing output: " + str(path))
    destination.mkdir()
    write_json(destination / "calibration.json", calibration)
    write_json(destination / "file_budget_before.json", extra_budget)
    write_json(config_path, config)
    write_json(array_path, settings)
    with tarfile.open(destination / "pilots.tar", "w") as archive:
        for index, pilot in enumerate(pilots, 1):
            archive.add(pilot["workspace"], arcname=f"attempt{index}")
        probe = Path(os.environ["TMPDIR"]) / "e1_long_cpu16_probe.py"
        if probe.is_file():
            archive.add(probe, arcname=probe.name)
    plan = prepare(config_path, array_path, output, jobs,
                   ROOT / "deps/LinearCParty-89baea52", ROOT / "build/linear-cparty-89baea52")
    summary = []
    historical_walls = json.loads((DOCS / "historical_wall_estimates.json").read_text())["tiers"]
    lengths = {i["name"]: i["length"] for i in base["inputs"]}
    for tier in settings["tiers"]:
        selected = [t for t in plan["tasks"] if t["tier"] == tier["name"]]
        old_selected = [t for t in base["tasks"] if t["tier"] == tier["name"]]
        observed = next(t for t in historical_walls if t["old_run"] == HISTORY.name and t["tier"] == tier["name"])
        summary.append({"tier": tier["name"], "inputs": sum(len(t["inputs"]) for t in selected),
                        "tasks": len(selected), "workers": 8, "h_rt": tier["h_rt"],
                        "reserved_cpu16_hours": len(selected) * wall_seconds(tier["h_rt"])/3600,
                        "max_predicted_task_hours": max(t["estimated_seconds"] for t in selected)/3600,
                        "predicted_task_hours_sum": sum(t["estimated_seconds"] for t in selected)/3600,
                        "predicted_sequential_cli_hours": sum(t["estimated_sequential_seconds"] for t in selected)/3600,
                        "old_cpu40_tasks": len(old_selected), "old_cpu40_h_rt": old_selected[0]["h_rt"],
                        "old_reserved_cpu40_hours": sum(wall_seconds(t["h_rt"]) for t in old_selected)/3600,
                        "historical_cli_hours": sum(r["elapsed_seconds"] for r in history["rows"] if
                            tier["min_length"] <= lengths[r["input"]] <= tier["max_length"])/3600,
                        "historical_task_wall_hours": observed["known_attempt_wall_hours"]})
    after = snapshot([BASE, old_jobs])
    if before != after:
        raise ValueError("previous run or scripts changed during preparation")
    original = json.loads((BASE / "original_config.json").read_text())
    if {k:v for k,v in config.items() if k not in ("name", "memory_limit_mib")} != {
            k:v for k,v in original.items() if k not in ("name", "memory_limit_mib")}:
        raise ValueError("scientific config changed")
    if [(i["name"],i["sha256"]) for i in base["inputs"]] != [(i["name"],i["sha256"]) for i in plan["inputs"]]:
        raise ValueError("frozen input order or bytes changed")
    if any(t["estimated_seconds"] > .7 * wall_seconds(t["h_rt"]) or wall_seconds(t["h_rt"]) > 14400 for t in plan["tasks"]):
        raise ValueError("task exceeds requested wall limits")
    expected = [j["id"] for t in base["tasks"] for j in task_jobs(base, t)]
    planned = [j["id"] for t in plan["tasks"] for j in task_jobs(plan, t)]
    if len(planned) != len(expected) or set(planned) != set(expected):
        raise ValueError("repacked matrix changed or duplicated a CLI condition")
    for task in plan["tasks"]:
        saved = json.loads(Path(task["config_path"]).read_text())
        old_tier = next(t for t in base["array_config"]["tiers"] if t["name"] == task["tier"])
        if saved["order"] != {"kind": old_tier["order_kind"], "seed": base["array_config"]["seed"]}:
            raise ValueError("task ordering policy or seed changed")
        metadata = json.loads((output/task["directory"]/"task.json").read_text())
        if metadata["status"] != "not_started" or metadata["attempts"] or (output/task["directory"]/"run").exists():
            raise ValueError("new plan has already executed a task")
    result = {"name": args.name, "tiers": summary, "old_run_before": before, "old_run_after": after,
              "old_run_and_scripts_unchanged": True, "scientific_config_unchanged_except_name_and_memory_cap": True,
              "input_roster_order_and_bytes_identical": True, "plan_sha256": digest(output/"array_plan.json"),
              "cli_matrix_identity_verified": True, "planned_cli_jobs": len(planned),
              "task_order_policy_and_seed_identical": True, "all_tasks_unstarted": True,
              "job_scripts": plan["job_scripts"], "file_budget": plan["file_budget"],
              "submitted": False, "git_commit_created": False}
    write_json(destination / "summary.json", result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
