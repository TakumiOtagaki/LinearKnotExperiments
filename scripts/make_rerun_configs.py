#!/usr/bin/env python3
"""Copy adopted science settings and generate cpu_16 release/parallelism configs."""

import argparse
import copy
import json
import math
from pathlib import Path
import re
import sys

from plan_array import RESOURCES, positive, prepare, validate_array, wall_seconds
from rna_records import read_record
from run_cli import ROOT, digest, experiment_path, now, write_json

TEMPLATES = (
    "e1_bprna_short_d21618c5_20261014", "e1_bprna_long_d21618c5_20261014",
    "e2_rnacentral_d21618c5_20261014", "e2_synthetic_shapes_d21618c5_20261014",
    "e2_memory_breakdown_d21618c5_20261016",
)
PARALLEL_REASON = ("Author policy 2026-10-16: fixed single-thread CLI concurrency on cpu_16; "
                   "record workers and describe concurrent wall-time shifts in Supplementary.")


def saved_path(path):
    path = path.resolve()
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def generate(settings):
    required = ("revision", "date", "software_root", "build_dir", "software_lock",
                "stage_build_dir", "stage_software_lock", "workers", "memory_limit_mib")
    if any(k not in settings for k in required):
        raise ValueError("settings require " + ", ".join(required))
    if not re.fullmatch(r"[A-Za-z0-9_-]+", settings["revision"]) or not re.fullmatch(r"\d{8}", settings["date"]):
        raise ValueError("revision must be a simple name and date must be YYYYMMDD")
    if settings.get("resource", "cpu_16") != "cpu_16":
        raise ValueError("release reruns use cpu_16")
    slowdown = positive(settings.get("parallel_slowdown_factor", 1), "parallel_slowdown_factor")
    if slowdown < 1:
        raise ValueError("parallel_slowdown_factor must be at least 1")
    variants = []
    for template in TEMPLATES:
        cp, ap = ROOT / "configs" / (template + ".json"), ROOT / "configs/arrays" / (template + ".json")
        config, array = json.loads(cp.read_text()), json.loads(ap.read_text())
        name = template.split("_d21618c5_")[0] + "_" + settings["revision"] + "_" + settings["date"]
        stage = "memory_breakdown" in template
        config["name"] = name.replace("_", "-")
        config["software_lock"] = settings["stage_software_lock" if stage else "software_lock"]
        config["memory_limit_mib"] = settings["memory_limit_mib"]
        if "timeout_seconds" in settings:
            config["timeout_seconds"] = settings["timeout_seconds"]
        config.setdefault("environment", {}).update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
        config["rerun"] = {"template_config": saved_path(cp), "template_config_sha256": digest(cp),
                           "template_array": saved_path(ap), "template_array_sha256": digest(ap),
                           "software_root": settings["software_root"],
                           "build_dir": settings["stage_build_dir" if stage else "build_dir"],
                           "workers": settings["workers"], "parallel_slowdown_factor": slowdown,
                           "memory_basis": "Explicit per-process RLIMIT_AS cap; allocation feasibility is not a measured RSS prediction.",
                           "historical_template_analysis": copy.deepcopy(config.get("analysis")),
                           "original_cost_models": {t["name"]: copy.deepcopy(t["cost_model"]) for t in array["tiers"]}}
        # Preserve metric definitions; replace resource descriptions inherited from old runs.
        analysis = config.setdefault("analysis", {})
        if "scope" in analysis:
            analysis["scope"] = re.sub(r"node_f, one thread", "cpu_16, one thread per process", analysis["scope"])
        if "memory_limit_kind" in analysis:
            analysis["memory_limit_kind"] = f"Per-process RLIMIT_AS address-space cap, {settings['memory_limit_mib']} MiB; not an RSS limit."
        array["compact_plan"] = True
        array["max_concurrent_tasks"] = 30
        array["note"] = "Release rerun; cpu_16; time coefficients inherited from the recorded template; sequential-work wall estimate."
        if array["measurement"] == "timing":
            array["parallel_timing"] = {"allowed": True, "reason": PARALLEL_REASON}
        for tier in array["tiers"]:
            tier.update(resource="cpu_16", workers=settings["workers"],
                        memory_budget_mib=settings.get("memory_budget_mib", math.floor(RESOURCES["cpu_16"][1] * 1e9 / 2**20)),
                        task_overhead_mib=settings.get("task_overhead_mib", 1024))
            limits = tier.setdefault("job_limits", {})
            limits["memory_limit_mib"] = settings["memory_limit_mib"]
            if "timeout_seconds" in settings:
                limits["timeout_seconds"] = settings["timeout_seconds"]
            if "h_rt" in settings:
                tier["h_rt"] = settings["h_rt"]
                tier["target_seconds"] = math.floor(wall_seconds(tier["h_rt"]) * min(tier.get("max_wall_fraction", .7), .7))
            if "max_inputs_per_task" in settings:
                if type(settings["max_inputs_per_task"]) is not int or settings["max_inputs_per_task"] < 1:
                    raise ValueError("max_inputs_per_task must be a positive integer")
                tier["max_inputs_per_task"] = settings["max_inputs_per_task"]
            model = tier["cost_model"]
            model["safety_factor"] = model.get("safety_factor", 1) * slowdown
            model["basis"] += f"; rerun time slowdown factor {slowdown}; memory envelope replaced by explicit {settings['memory_limit_mib']} MiB RLIMIT_AS cap (not a calibrated RSS prediction)"
            model.update(memory_base_mib=settings["memory_limit_mib"], mib_per_nt=0)
            for fit in model.get("per_mode", {}).values():
                fit.update(memory_base_mib=settings["memory_limit_mib"], memory_coefficient=0)
        variants.append((name, config, array))

    # Paired Supplementary configs use exactly the same three S/G, modes and repeats.
    _, base_config, base_array = next(v for v in variants if v[0].startswith("e2_rnacentral_"))
    selected = [i for i in base_config["inputs"] if 4500 <= len(read_record(experiment_path(i["path"]))[1]) <= 7100]
    lengths = [len(read_record(experiment_path(i["path"]))[1]) for i in selected]
    if len(selected) != 3:
        raise ValueError("RNAcentral parallel check requires the frozen 4.8/5.5/7.0 kb inputs")
    for workers in (1, 8):
        config, array = copy.deepcopy(base_config), copy.deepcopy(base_array)
        name = f"e2_parallel_impact_w{workers}_{settings['revision']}_{settings['date']}"
        config["name"] = name.replace("_", "-")
        config["inputs"] = selected
        if "runtime_eligibility" in config:
            config["runtime_eligibility"]["selected_ids"] = [i["name"] for i in selected]
        config["rerun"].update(workers=workers, paired_workers=[1, 8], selected_lengths=lengths,
                               scope="Supplementary concurrency check; compare whole-CLI wall time at fixed cpu_16 and thread count.")
        array["tiers"] = [copy.deepcopy(base_array["tiers"][0])]
        # Keep all 9 combinations on one task, making 8 CLI slots available.
        tier = array["tiers"][0]
        tier.update(name="parallel_impact", workers=workers, max_inputs_per_task=3,
                    h_rt=settings.get("parallel_check_h_rt", "12:00:00"))
        tier["target_seconds"] = math.floor(wall_seconds(tier["h_rt"]) * .7)
        variants.append((name, config, array))

    for _, config, array in variants:
        # Full input/provenance checks belong to the planner. Validate resource
        # settings here without rereading every frozen RNA during generation.
        validation_config = copy.deepcopy(config)
        if "input_cohort" in config and not config.get("inputs"):
            cohort = json.loads(experiment_path(config["input_cohort"]).read_text())
            validation_config["inputs"] = [{"name": i["id"], "length": i["length"],
                                             "path": "resolved-by-array-planner"} for i in cohort["records"]]
        validate_array(validation_config, array)
    return variants


def write_configs(settings, output, dry_run=False, estimate_only=False):
    if estimate_only and not dry_run:
        raise ValueError("--estimate-only requires --dry-run")
    variants = generate(settings)
    if output.exists():
        raise ValueError("use a new config output directory")
    output.mkdir(parents=True)
    (output / "arrays").mkdir()
    report = {"schema_version": 1, "created_utc": now(), "settings": settings,
              "dry_run": dry_run, "estimate_only": estimate_only, "variants": []}
    for name, config, array in variants:
        cp, ap = output / (name + ".json"), output / "arrays" / (name + ".json")
        write_json(cp, config)
        write_json(ap, array)
        row = {"name": name, "config": saved_path(cp), "array_config": saved_path(ap),
               "software_root": settings["software_root"], "build_dir": config["rerun"]["build_dir"],
               "memory_cap_envelope_mib": max(t["workers"] * config["memory_limit_mib"] + t["task_overhead_mib"] for t in array["tiers"])}
        if dry_run:
            estimate = prepare(cp, ap, output / ("estimate_" + name), output / ("jobs_" + name),
                               dry_run=True, estimate_only=True)
            row["estimate"] = {k: v for k, v in estimate.items() if k != "task_assignments"}
        report["variants"].append(row)
        write_json(output / "rerun_manifest.json", report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--settings", type=Path, required=True, help="small release/resource JSON")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true", help="write configs and offline planner estimates; never execute or submit")
    parser.add_argument("--estimate-only", action="store_true")
    args = parser.parse_args()
    report = write_configs(json.loads(args.settings.read_text()), args.out_dir.resolve(), args.dry_run, args.estimate_only)
    print(json.dumps({"out_dir": str(args.out_dir), "variants": [
        {"name": v["name"], **{k: v.get("estimate", {}).get(k) for k in ("inputs", "tasks", "planned_jobs")}} for v in report["variants"]]}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, OSError) as exc:
        print("error: " + str(exc), file=sys.stderr)
        sys.exit(2)
