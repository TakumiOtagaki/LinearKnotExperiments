#!/usr/bin/env python3
"""Prepare (never submit) the CParty/LinearKnot PF scaling arrays."""
import argparse
from collections import defaultdict
import json
import math
from pathlib import Path

from plan_array import prepare
from run_cli import ROOT, digest, mode_applies, now, plan_jobs, write_json

DOC = ROOT / "docs/cparty_scaling_89baea52_20261020"
RUN = "cparty-linearknot-pf-89baea52-20261020"
LENGTHS = (250, 500, 1000, 2000, 4000)
METHODS = ("cparty", "exact", "beam100")


def calibrate(pilot):
    rows = pilot["observations"]
    fits = {}
    for method in METHODS:
        successful = [r for r in rows if r["method"] == method and r["status"] == "success"]
        if not successful:
            raise ValueError("no successful calibration for " + method)
        exponent = 3 if method != "beam100" else 1
        # Exact need not share CParty's cubic empirical behavior. Use the larger
        # observed 250->500 exponent if necessary; keep timeouts as censored data.
        if method == "exact":
            by_shape = defaultdict(dict)
            for r in successful:
                by_shape[r["family"]].setdefault(r["length"], []).append(r["elapsed_seconds"])
            for lengths in by_shape.values():
                if 250 in lengths and 500 in lengths:
                    exponent = max(exponent, math.log(max(lengths[500]) / max(lengths[250]), 2))
        time_coefficient = max(r["elapsed_seconds"] / r["length"]**exponent for r in successful)
        # RSS and AS differ. Use a 50% RSS envelope for extrapolation, then a
        # separate explicit per-process AS cap for the allocation feasibility.
        memory_exponent = 1 if method == "beam100" else 2
        memory_coefficient = max(r["peak_rss_kib"] / 1024 / r["length"]**memory_exponent for r in successful)
        predictions = {}
        for n in LENGTHS:
            seconds = 1.5 * time_coefficient * n**exponent
            mib = 1.5 * memory_coefficient * n**memory_exponent
            predictions[n] = {"seconds": seconds, "rss_mib": mib,
                              "within_per_process_limit": seconds <= 21600 and mib <= 81920}
        fits[method] = {"time_exponent": exponent, "time_coefficient": time_coefficient,
                        "rss_exponent": memory_exponent, "rss_mib_coefficient": memory_coefficient,
                        "planning_factor": 1.5, "predictions": predictions,
                        "successful_calibration_jobs": len(successful),
                        "censored_jobs": [r for r in rows if r["method"] == method and r["status"] == "timeout"],
                        "numerical_failures": [r for r in rows if r["method"] == method and r["status"] == "failed"]}
    return fits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", type=Path, default=DOC / "calibration_pilot.json")
    args = parser.parse_args()
    pilot = json.loads(args.pilot.read_text())
    if not pilot.get("finished_utc"):
        raise ValueError("finish the development pilot before freezing the plan")
    config_path = ROOT / "configs" / (RUN + ".json")
    array_path = ROOT / "configs/arrays" / (RUN + ".json")
    if any(p.exists() for p in (config_path, array_path, ROOT / "results/runs" / RUN, ROOT / "jobs" / RUN)):
        raise ValueError("refuse to overwrite existing config/plan/jobs")
    fits = calibrate(pilot)
    original = ROOT / "results/runs/e2-synthetic-shapes-89baea52-20261018/array_plan.json"
    e2 = json.loads(original.read_text())
    selected = [{k: x[k] for k in ("name", "path", "input_archive", "input_member", "sha256", "family", "seed", "length")}
                for x in e2["inputs"] if x["family"] != "branch_depth"]
    if len(selected) != 105 or any(sum(x["length"] == n for x in selected) != 21 for n in LENGTHS):
        raise ValueError("expected seven shapes x three seeds x five lengths")
    modes, decisions = [], []
    for method in METHODS:
        admitted = []
        for n in LENGTHS:
            prediction = fits[method]["predictions"][n]
            # The absolute 6h/80GiB bound is checked first. Also require a
            # single complete 1+3 matrix to fit half a 24h scheduler limit.
            accepted = prediction["within_per_process_limit"] and 4 * prediction["seconds"] + 120 <= 43200
            decisions.append({"method": method, "length": n, "included": accepted,
                              **prediction, "exclusion_reason": None if accepted else "6h/job or 80GiB/process, or 1+3 matrix cannot fit half of 24h"})
            if accepted:
                admitted.append(n)
        if not admitted or admitted != list(LENGTHS[:len(admitted)]):
            raise ValueError("expected a nonempty contiguous length prefix for each method")
        modes.append({"name": method, "args": ["--method", method],
                      "length_range": {"min_length": min(admitted), "max_length": max(admitted)}})
    config = {"schema_version": 1, "name": RUN, "purpose": "formal",
              "energy_parameter_file": "params/rna_Turner04.par",
              "model_settings": {"temperature_celsius": 37, "dangles": 2,
                                 "wmb_closure_policy": "LegacyCompatible", "build_revision_verified": True},
              "environment": {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "LC_ALL": "C"},
              "timeout_seconds": 21600, "memory_limit_mib": 4096, "warmups": 1, "repetitions": 3,
              "order": {"kind": "blocked_random", "seed": 20261020},
              "software_lock": "manifests/software-pf-scaling-89baea52-20261020.json",
              "inputs": selected, "modes": modes, "programs": [{"name": "PFScaling", "args": []}],
              "analysis": {"scope": "Inside-only conditional PF; no MFE, outside, sampling or BPP; compatible weak-closure setting differs from adopted main-text default",
                           "log_z_absolute_tolerance": 1e-9, "free_energy_absolute_tolerance": 1e-10,
                           "input_source_plan_sha256": digest(original), "development_pilot_sha256": digest(args.pilot)}}
    tiers, envelopes = [], []
    for n in LENGTHS:
        active = [m for m in modes if mode_applies(m, {"length": n})]
        if not active:
            continue
        predicted = {m["name"]: fits[m["name"]]["predictions"][n] for m in active}
        max_rss = max(p["rss_mib"] for p in predicted.values())
        # AS headroom beyond the already padded RSS model; keep all caps small
        # enough for cpu_16, reducing worker count if a larger cap is needed.
        cap = 1024
        while cap < max(1024, 1.3 * max_rss):
            cap *= 2
        workers = min(8, (35095 - 1024) // cap)
        if workers < 1:
            raise ValueError("prediction needs a larger resource; prepare a new explicit tier")
        # Fixed-size batch envelope. Work/W + longest chain bounds seeded
        # work-conserving dispatch. Sum the conservative bound across methods
        # and four executions, then add 120 s/task. Encode wall credits in the
        # existing planner, whose sum remains a conservative task wall estimate.
        total = sum(p["seconds"] for p in predicted.values())
        batch = 21
        while batch > 1:
            bound = 4 * (batch / workers + (1 - 1 / workers)) * total + 120
            if 2 * bound <= 86400:
                break
            batch -= 1
        credit = 1 / workers + (1 - 1 / workers) / batch
        full_bound = 4 * batch * credit * total + 120
        if 2 * full_bound > 86400:
            raise ValueError("complete matrix cannot fit 24h with twofold h_rt allowance")
        # Round upward in 5-minute increments, retaining about twofold padding.
        h_seconds = max(600, math.ceil(2 * full_bound / 300) * 300)
        if h_seconds > 86400:
            raise ValueError("rounded h_rt exceeds 24h")
        h_rt = f"{h_seconds//3600:02d}:{h_seconds%3600//60:02d}:00"
        # A partial batch has a larger chain fraction. Independently check its
        # physical wall bound below; full-batch h_rt also covers this shorter tail.
        per_mode = {m["name"]: {"constant_seconds": 120 / batch / 4 / len(active),
                                 "seconds_coefficient": fits[m["name"]]["time_coefficient"] * 1.5 * credit,
                                 "length_exponent": fits[m["name"]]["time_exponent"],
                                 "memory_base_mib": cap, "memory_coefficient": 0}
                    for m in active}
        tier = {"name": f"n{n}", "min_length": n, "max_length": n, "resource": "cpu_16",
                "h_rt": h_rt, "target_seconds": full_bound + 1, "max_wall_fraction": .7,
                "max_inputs_per_task": batch, "workers": workers, "memory_budget_mib": 35095,
                "task_overhead_mib": 1024,
                "job_limits": {"memory_limit_mib": cap, "timeout_seconds": 21600},
                "cost_model": {"basis": "Development 250/500/1000 pilot; 1.5x time envelope, CParty cubic/RSS quadratic, exact max(cubic, observed exponent), beam linear. Credits include fixed batch list bound and 120s overhead; every actual batch is independently checked. Not an asymptotic complexity claim.",
                               "seconds_per_nt_condition": 1e-12, "mib_per_nt": 0,
                               "memory_base_mib": cap, "safety_factor": 1, "per_mode": per_mode}}
        tiers.append(tier)
        envelopes.append({"length": n, "batch_size": batch, "workers": workers, "cap_mib": cap,
                          "per_method_prediction": predicted, "full_batch_bound_seconds": full_bound,
                          "h_rt_seconds": h_seconds, "credit": credit})
    settings = {"schema_version": 1, "measurement": "timing", "seed": 20261020,
                "compact_plan": True, "max_concurrent_tasks": 30,
                "file_budget": {"project": 102571, "filesystem": "/gs/bs", "hard_file_limit": 24000000, "max_remaining_fraction": .1},
                "parallel_timing": {"allowed": True, "reason": "Author-requested E2 procedure: bounded node-local single-thread processes, one warmup and three measured repetition blocks; save workers/host and report concurrent wall shifts."},
                "tiers": tiers}
    # Preview packing before mutating a run: check each actual batch's physical
    # wall envelope, including a final short batch's critical path and overhead.
    from plan_array import make_tasks
    from run_cli import digest as sha
    preview = [{**x, "sequence_sha256": next(y["sequence_sha256"] for y in e2["inputs"] if y["name"] == x["name"])} for x in selected]
    tasks = make_tasks(config, settings, preview)
    for t in tasks:
        e = next(e for e in envelopes if t["tier"] == f"n{e['length']}")
        k = len(t["inputs"])
        actual_bound = 4 * (k / t["workers"] + 1 - 1 / t["workers"]) * sum(p["seconds"] for p in e["per_method_prediction"].values()) + 120
        if actual_bound > .5 * e["h_rt_seconds"] + 1:
            raise ValueError("partial batch needs a fresh longer h_rt")
    write_json(config_path, config)
    write_json(array_path, settings)
    plan = prepare(config_path, array_path, ROOT / "results/runs" / RUN, ROOT / "jobs" / RUN,
                   ROOT / "deps/LinearCParty-89baea52", ROOT / "build/pf-scaling-89baea52-20261020")
    expected_hashes = {x["name"]: x["sha256"] for x in selected}
    if any(x["sha256"] != expected_hashes[x["name"]] for x in plan["inputs"]):
        raise ValueError("new frozen input differs from adopted E2")
    budgets = []
    actual_bounds = []
    for task in plan["tasks"]:
        envelope = next(e for e in envelopes if task["tier"] == f"n{e['length']}")
        seconds = 4 * (len(task["inputs"]) / task["workers"] + 1 - 1 / task["workers"]) * sum(p["seconds"] for p in envelope["per_method_prediction"].values()) + 120
        actual_bounds.append({"tier": task["tier"], "task_id": task["task_id"], "wall_bound_seconds": seconds})
    for tier in tiers:
        assigned = [t for t in plan["tasks"] if t["tier"] == tier["name"]]
        wall = next(e["h_rt_seconds"] for e in envelopes if tier["name"] == f"n{e['length']}")
        bounds = [b["wall_bound_seconds"] for b in actual_bounds if b["tier"] == tier["name"]]
        predicted_seconds = sum(bounds)
        budgets.append({"tier": tier["name"], "tasks": len(assigned), "resource": "cpu_16=1",
                        "workers": tier["workers"], "cap_mib": tier["job_limits"]["memory_limit_mib"],
                        "h_rt": tier["h_rt"], "reserved_task_hours": len(assigned) * wall / 3600,
                        "estimated_task_hours": predicted_seconds / 3600,
                        "estimated_points": sum(.06 * (.7 * max(seconds, 300) + .1 * wall) / 3600 for seconds in bounds),
                        "h_rt_envelope_points": len(assigned) * .06 * (.7 * max(wall, 300) + .1 * wall) / 3600})
    write_json(DOC / "preparation.json", {"created_utc": now(), "no_qsub": True, "no_commit": True,
               "source_e2_plan_sha256": sha(original), "pilot_sha256": sha(args.pilot), "fits": fits,
               "method_length_decisions": decisions, "batch_envelopes": envelopes, "actual_task_wall_bounds": actual_bounds, "budgets": budgets,
               "total_h_rt_envelope_points": sum(b["h_rt_envelope_points"] for b in budgets),
               "total_estimated_points": sum(b["estimated_points"] for b in budgets),
               "billing_source": "https://www.somuka.titech.ac.jp/reiki_int/reiki_honbun/x385RG00001693.html",
               "billing_checked_date": "2026-10-20", "billing_priority": 0,
               "billing_note": "cpu_16 coefficient .060; .7*max(wall,300)+.1*h_rt; /3600. Envelope uses h_rt for wall. Before 0.0001-point truncation; storage and retries excluded.",
               "plan_sha256": sha(ROOT / "results/runs" / RUN / "array_plan.json"),
               "file_budget": plan["file_budget"], "all_tasks_not_started": all(not json.loads((Path(plan["output_dir"]) / t["directory"] / "task.json").read_text())["attempts"] for t in plan["tasks"]),
               "planned_jobs": len(plan_jobs(config))})
    print(json.dumps({"tasks": len(plan["tasks"]), "budgets": budgets}, indent=2))


if __name__ == "__main__":
    main()
