#!/usr/bin/env python3
"""Prepare and verify resource reruns after strict release parity; never submit."""
import argparse
from collections import Counter
import copy
from datetime import datetime
import json
import math
from pathlib import Path

from file_budget import check_budget
from plan_array import load_plan, prepare, task_jobs, wall_seconds
from run_cli import ROOT, digest, mode_applies, now, plan_jobs, verify_lock, write_json

DOC = ROOT / "docs/linearknot_4f97fd3_20261021"
CONFIG = ROOT / "configs/linearknot_4f97fd3_20261021"
SOURCE = ROOT / "deps/LinearKnot-4f97fd3"
BUILD = ROOT / "build/linearknot-4f97fd3"
COMPARISON = ROOT / "results/local/linearknot-4f97fd3-comparison-20261021/comparison.json"
MAP = {"LinearCParty": "fold", "LinearPKCapR": "profile", "LinearPKRaccess": "access",
       "LinearPKSampling": "sample", "LinearPKCentroid": "centroid"}
RUNS = ["e2-rnacentral-89baea52-20261018", "e2-synthetic-shapes-89baea52-20261018",
        "e2-memory-breakdown-89baea52-20261018", "e2-parallel-impact-w1-89baea52-20261018",
        "e2-parallel-impact-w10-89baea52-20261018",
        *[f"e3-neat1-{x}-89baea52-20261018" for x in ("bpp", "capr", "raccess", "raccess-byloop", "sampling")],
        "e3-neat1-centroid-89baea52-20261020", "cparty-linearknot-pf-89baea52-20261020"]


def relative(p):
    p = Path(p)
    return str(p.relative_to(ROOT)) if p.is_absolute() and p.is_relative_to(ROOT) else str(p)


def release_gate():
    r = json.loads(COMPARISON.read_text())
    if (r.get("status") != "complete" or not r.get("successful") or not r.get("identical")
            or r.get("comparisons") != 378 or r.get("identical_comparisons") != 378):
        raise ValueError("all 378 strict release comparisons must pass first")
    if not r["fingerprints"]["candidate"]["source"]["commit"].startswith("4f97fd3"):
        raise ValueError("comparison candidate differs")
    if not r["fingerprints"]["reference"]["source"]["commit"].startswith("9c7d6e1"):
        raise ValueError("comparison reference differs")
    for n, h in r["artifacts_sha256"].items():
        if digest(COMPARISON.parent / n) != h:
            raise ValueError("comparison artifact changed")
    for key, file in (("cmake_cache_sha256", BUILD / "CMakeCache.txt"),
                      ("parameter_sha256", SOURCE / "params/rna_Turner04.par")):
        if digest(file) != r["fingerprints"]["candidate"][key]:
            raise ValueError("compared dependency changed")
    if digest(BUILD / "linearknot") != r["fingerprints"]["candidate"]["binaries"]["linearknot"]:
        raise ValueError("compared CLI changed")
    return r


def frozen_inputs(plan):
    keys = ("name", "path", "input_archive", "input_member", "sha256", "family", "seed", "length")
    return [{k: relative(x[k]) if k in ("path", "input_archive") else x[k]
             for k in keys if k in x} for x in plan["inputs"]]


def calibration(old, config, settings, source_path):
    hashes = {relative(source_path): digest(source_path)}
    evidence = []
    for tier in settings["tiers"]:
        tasks = [t for t in old["tasks"] if t["tier"] == tier["name"]]
        rates, walls, max_cli, max_batch = [], [], 0, 0
        for task in tasks:
            directory = Path(old["output_dir"]) / task["directory"]
            meta_path = directory / "task.json"
            rows_path = directory / "run/task_results.jsonl"
            meta = json.loads(meta_path.read_text())
            if meta["status"] not in ("success", "failed") or not meta["attempts"]:
                raise ValueError("historical task is unfinished")
            attempts = meta["attempts"]
            if any(not a.get("finished_utc") for a in attempts):
                raise ValueError("historical attempt lacks an end time")
            wall = sum((datetime.fromisoformat(a["finished_utc"]) -
                        datetime.fromisoformat(a["started_utc"])).total_seconds() for a in attempts)
            modes = {m["name"] for m in config["modes"]}
            expected = [j for j in task_jobs(old, task) if j["mode"] in modes]
            observed = [json.loads(line) for line in rows_path.read_text().splitlines()]
            selected = [j for j in observed if j["mode"] in modes]
            if ({j["id"] for j in selected} != {j["id"] for j in expected} or
                    len(selected) != len(expected) or any(j["status"] != "success" for j in selected)):
                raise ValueError("retained historical matrix has missing or failed jobs")
            max_cli = max(max_cli, max(j["elapsed_seconds"] for j in selected))
            rates.append(wall / len(selected)); walls.append(wall)
            max_batch = max(max_batch, len(task["inputs"]))
            for p in (meta_path, rows_path):
                hashes[relative(p)] = digest(p)
            evidence.append({"tier": tier["name"], "task_id": task["task_id"],
                             "inputs": len(task["inputs"]), "retained_jobs": len(selected),
                             "previous_task_wall_seconds": wall,
                             "statuses_all_methods": dict(Counter(j["status"] for j in observed))})
        # Fit allocation wall credits at unchanged concurrency. These are not
        # estimates of sequential CLI work, nor a measured new-release speedup.
        rate = max(rates)
        jobs_per_input = max(sum(mode_applies(m, x) for m in config["modes"]) *
                             len(config["programs"]) * (config["warmups"] + config["repetitions"])
                             for x in old["inputs"] if tier["min_length"] <= x["length"] <= tier["max_length"])
        base = max(max(walls), max_batch * jobs_per_input * rate, max_cli)
        seconds = max(600, math.ceil((2 * base + 240) / 300) * 300)
        if seconds > 86400:
            raise ValueError("historical bound needs an explicit longer-tier decision")
        tier.update(h_rt=f"{seconds//3600:02d}:{seconds%3600//60:02d}:00",
                    target_seconds=.7 * seconds - 1, max_wall_fraction=.7,
                    max_inputs_per_task=max_batch)
        tier.pop("parallel_packing", None)
        tier.pop("order_kind", None)
        cap = tier.get("job_limits", {}).get("memory_limit_mib", config["memory_limit_mib"])
        timeout = tier.get("job_limits", {}).get("timeout_seconds", config["timeout_seconds"])
        tier["job_limits"] = {"memory_limit_mib": cap, "timeout_seconds": min(timeout, seconds - 60)}
        tier["cost_model"] = {
            "basis": "Allocation-wall credits: maximum prior complete task wall / retained scheduled jobs at the same worker count and batch cap; 1.25x allowance, no new-release speed discount. Not sequential CLI work or a complexity claim. h_rt >= twice the fitted full batch and prior maximum task plus 240 s, rounded up to 5 min.",
            "seconds_per_nt_condition": 1e-12, "seconds_constant_per_condition": 1.25 * rate,
            "safety_factor": 1, "memory_base_mib": cap, "mib_per_nt": 0,
            "calibration": {"source_plan": relative(source_path), "max_wall_credit_per_job": rate,
                            "max_previous_task_seconds": max(walls), "max_cli_seconds": max_cli,
                            "max_previous_batch_inputs": max_batch}}
    return {"source_hashes": hashes, "tasks": evidence}


def write_once(path, data):
    if path.exists():
        if json.loads(path.read_text()) != data:
            raise ValueError("existing preparation file differs: " + str(path))
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        write_json(path, data)


def verify_preparation():
    release_gate()
    report = json.loads((DOC / "preparation.json").read_text())
    for path, h in report["files_sha256"].items():
        if digest(ROOT / path) != h:
            raise ValueError("preparation artifact changed: " + path)
    for entry in report["matrices"]:
        p = ROOT / entry["plan"]
        plan = load_plan(p)
        verify_lock(p.parent / "software_lock.json", SOURCE, ROOT / entry["build"],
                    plan["config"]["programs"], SOURCE / plan["config"]["energy_parameter_file"])
        for t in plan["tasks"]:
            directory = p.parent / t["directory"]
            meta = json.loads((directory / "task.json").read_text())
            if meta["status"] != "not_started" or meta["attempts"] or (directory / "run").exists():
                raise ValueError("resource task already attempted: " + str(directory))
    for path, h in report["historical_source_hashes"].items():
        if digest(ROOT / path) != h:
            raise ValueError("historical observation changed: " + path)
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--verify", action="store_true")
    args = p.parse_args()
    if args.verify or (DOC / "preparation.json").exists():
        r = verify_preparation()
        print(json.dumps({"verified": True, "matrices": len(r["matrices"]), "tasks": r["tasks"]}))
        return
    release_gate()
    budget = check_budget({"estimated_shared_entries": 5000,
                           "basis": "All new timing plans, configs, observer builds, tests, reports, first execution and scheduler logs; clone/CLI build already certified separately"}, {})
    write_json(DOC / "timing_file_budget.json", budget)
    matrices, scripts, history, calibrations = [], [], {}, {}
    for old_name in RUNS:
        old_path = ROOT / "results/runs" / old_name / "array_plan.json"
        old = load_plan(old_path)
        config = copy.deepcopy(old["config"])
        settings = copy.deepcopy(old["array_config"])
        pf = old_name.startswith("cparty-")
        memory = "memory-breakdown" in old_name
        new_name = ("linearknot-pf-4f97fd3-20261021" if pf else
                    old_name.replace("89baea52-20261018", "4f97fd3-20261021").replace("89baea52-20261020", "4f97fd3-20261021"))
        config.update(name=new_name, inputs=frozen_inputs(old))
        if pf:
            config["modes"] = [m for m in config["modes"] if m["name"] != "cparty"]
            config["analysis"]["cparty_reference_plan"] = relative(old_path)
            build = ROOT / "build/pf-scaling-4f97fd3-20261021"
            lock = "manifests/software-pf-scaling-4f97fd3-20261021.json"
        elif memory:
            build = ROOT / "build/stage-timing-4f97fd3"
            lock = "manifests/software-stage-timing-4f97fd3-20261021.json"
        else:
            build, lock = BUILD, "manifests/software-linearknot-4f97fd3.json"
            config["programs"] = [{**x, "executable": "linearknot", "subcommand": MAP[x["name"]]}
                                  for x in config["programs"]]
        config["software_lock"] = lock
        config["analysis"]["scope"] = "Resource rerun on LinearKnot 4f97fd3; retain adopted scientific results; same inputs, model, outputs, seeds and repetition policy."
        config["rerun_preparation"] = {"release": json.loads((ROOT / lock).read_text())["source_commit"],
                                       "source_plan": relative(old_path), "source_plan_sha256": digest(old_path),
                                       "strict_comparison_sha256": digest(COMPARISON)}
        settings.update(measurement="timing", compact_plan=True, max_concurrent_tasks=30)
        settings["parallel_timing"] = {"allowed": True,
            "reason": "Retain the prior fixed single-thread worker count for matched resource reruns; allocation wall calibration keeps concurrency and host effects, without a new-release speed discount."}
        evidence = calibration(old, config, settings, old_path)
        history.update(evidence["source_hashes"]); calibrations[new_name] = evidence
        config_path, array_path = CONFIG / (new_name + ".json"), CONFIG / "arrays" / (new_name + ".json")
        write_once(config_path, config); write_once(array_path, settings)
        out, jobs = ROOT / "results/runs" / new_name, ROOT / "jobs" / new_name
        if (out / "array_plan.json").exists():
            plan = load_plan(out / "array_plan.json")
            if digest(config_path) != plan["original_config_sha256"] or digest(array_path) != plan["array_config_sha256"]:
                raise ValueError("existing plan differs")
        else:
            plan = prepare(config_path, array_path, out, jobs, SOURCE, build)
        if {x["name"]: x["sha256"] for x in plan["inputs"]} != {x["name"]: x["sha256"] for x in old["inputs"]}:
            raise ValueError("frozen input set differs")
        if config["modes"] != [m for m in old["config"]["modes"] if not pf or m["name"] != "cparty"]:
            raise ValueError("scientific modes differ")
        if [{k: v for k, v in x.items() if k not in ("executable", "subcommand")} for x in config["programs"]] != old["config"]["programs"]:
            raise ValueError("requested program outputs differ")
        matrices.append({"name": new_name, "source_plan": relative(old_path), "config": relative(config_path),
                         "array_config": relative(array_path), "plan": relative(out / "array_plan.json"),
                         "build": relative(build), "tasks": len(plan["tasks"]),
                         "cli_jobs": len(plan_jobs(config)), "inputs": len(plan["inputs"]),
                         "file_budget_entries": plan["file_budget"]["estimated_shared_entries"]})
        for script in plan["job_scripts"]:
            tier = next(t for t in settings["tiers"] if t["name"] == script["tier"])
            assigned = [t for t in plan["tasks"] if t["tier"] == script["tier"]]
            h = wall_seconds(tier["h_rt"]); coefficient = .15 if tier["resource"] == "cpu_40" else .06
            estimates = [t["estimated_seconds"] + 120 for t in assigned]
            if any(s > .7 * h for s in estimates):
                raise ValueError("transfer allowance exceeds planning fraction")
            scripts.append({"matrix": new_name, **script, "path": relative(script["path"]),
                "resource": tier["resource"] + "=1", "h_rt": tier["h_rt"], "workers": tier["workers"],
                "cap_mib": tier["job_limits"]["memory_limit_mib"],
                "timeout_seconds": tier["job_limits"]["timeout_seconds"],
                "previous_max_task_seconds": tier["cost_model"]["calibration"]["max_previous_task_seconds"],
                "estimated_task_hours": sum(estimates) / 3600,
                "h_rt_task_hours": len(assigned) * h / 3600,
                "estimated_points": sum(coefficient * (.7 * max(s, 300) + .1 * h) / 3600 for s in estimates),
                "h_rt_envelope_points": len(assigned) * coefficient * .8 * h / 3600})
    write_json(DOC / "calibration.json", calibrations)
    write_json(DOC / "job_scripts.json", scripts)
    files = [COMPARISON, DOC / "calibration.json", DOC / "job_scripts.json", Path(__file__).resolve(),
             ROOT / "scripts/run_cli.py", ROOT / "scripts/plan_array.py"]
    files += [ROOT / m[k] for m in matrices for k in ("config", "array_config", "plan")]
    files += [ROOT / s["path"] for s in scripts]
    report = {"created_utc": now(), "scheduler_submission": False, "git_commit": False,
              "matrices": matrices, "scripts": len(scripts), "tasks": sum(m["tasks"] for m in matrices),
              "cli_jobs": sum(m["cli_jobs"] for m in matrices), "historical_source_hashes": history,
              "files_sha256": {relative(f): digest(f) for f in files},
              "total_h_rt_envelope_points": sum(s["h_rt_envelope_points"] for s in scripts),
              "total_estimated_points": sum(s["estimated_points"] for s in scripts),
              "total_plan_execution_entries": sum(m["file_budget_entries"] for m in matrices),
              "billing_source": "https://www.somuka.titech.ac.jp/reiki_int/reiki_honbun/x385RG00001693.html",
              "billing_checked_date": "2026-10-21", "billing_priority": 0,
              "billing_note": "resource coefficient*(0.7*max(actual seconds,300)+0.1*h_rt seconds)/3600; truncation, storage, retries excluded",
              "timing_file_budget": budget}
    write_json(DOC / "preparation.json", report)
    verify_preparation()
    print(json.dumps({k: report[k] for k in ("tasks", "scripts", "cli_jobs", "total_h_rt_envelope_points")}))


if __name__ == "__main__":
    main()
