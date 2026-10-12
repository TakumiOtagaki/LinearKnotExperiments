#!/usr/bin/env python3
"""Prepare dry-run E2 shape and RNAcentral capb arrays; never submit jobs."""

import copy
import json
from pathlib import Path

from plan_array import prepare
from run_cli import ROOT, digest, write_json


def configurations():
    original = json.loads((ROOT / "configs/e2_rnacentral_15aa7742_timing_20261011.json").read_text())
    common = copy.deepcopy(original)
    common.update(purpose="formal", name="e2-rnacentral-15aa7742-capb-timing-20261012")
    common["model_settings"]["build_revision_verified"] = True
    common["modes"] = [{"name": "beam100_capb", "args": ["-b", "100", "--max-hairpin", "100", "--max-case3-depth", "100"]}]
    common["order"] = {"kind": "blocked_random", "seed": 20261012}
    common["analysis"] = {"scope": "Whole-CLI timing/RSS; node_f, one thread, warmup 1, measured 3. Explicit adopted h=d=max(b,30)=100; no exact long-RNA accuracy claim."}
    common["runtime_eligibility"]["reason"] = "Seven frozen, previously successful inputs; same S/G hashes, fixed 15aa7742. Re-measure only local caps at 100/100."
    synthetic = copy.deepcopy(common)
    synthetic.update(name="e2-synthetic-shapes-15aa7742-capb-timing-20261012")
    synthetic.pop("runtime_eligibility")
    manifest_path = ROOT / "data/e2-synthetic-shapes-20261012/manifest.json"
    manifest = json.loads(manifest_path.read_text())
    synthetic["inputs"] = manifest["inputs"]
    synthetic["input_generation"] = {"manifest": str(manifest_path), "sha256": digest(manifest_path),
                                     "same_sequence_shape_only": False, "main_inputs": 60, "total_inputs": 156}
    synthetic["analysis"]["scope"] += " Joint (S,G) synthetic controls; embedded canonical pairs change sequence. Run/calibrate n250/500 tiers before n1000/2000/4000."
    configs = {"e2_rnacentral_15aa7742_capb_timing_20261012.json": common,
               "e2_synthetic_shapes_15aa7742_capb_timing_20261012.json": synthetic}
    settings = {"schema_version": 1, "measurement": "timing", "seed": 20261012,
                "file_budget": {"project": 102571, "filesystem": "/gs/bs", "hard_file_limit": 24000000,
                                "max_remaining_fraction": .1}, "tiers": []}
    for n in (250, 500, 1000, 2000, 4000):
        settings["tiers"].append({"name": "synthetic_n" + str(n), "min_length": n, "max_length": n,
            "resource": "node_f", "h_rt": "12:00:00", "target_seconds": 24000,
            "max_inputs_per_task": 4, "workers": 1, "memory_budget_mib": 140000, "task_overhead_mib": 1024,
            "cost_model": {"basis": "Planning envelope only: historical isolated h=d=30 E2 worst wall/n is below 0.07 s/nt/CLI. Factor 3 allows cap100 and synthetic uncertainty; no calibrated shape model or asymptotic law. Start n250/500, then refit immutable larger-tier plan if needed.",
                           "seconds_per_nt_condition": .07, "seconds_constant_per_condition": 5,
                           "safety_factor": 3, "memory_base_mib": 2048, "mib_per_nt": 7}})
    rna = copy.deepcopy(settings)
    rna["tiers"] = []
    for low, high, name, wall, target in (
            (1, 8191, "rna_4_8kb", "20:00:00", 48000),
            (8192, 19999, "rna_8_20kb", "24:00:00", 60000),
            (20000, 30000, "rna_20_30kb", "24:00:00", 60000)):
        # RNAcentral maximum observed wall/n = <0.07; factor 2 reserves cap change.
        rna["tiers"].append({"name": name, "min_length": low, "max_length": high,
            "resource": "node_f", "h_rt": wall, "target_seconds": target,
            "max_inputs_per_task": 1, "workers": 1, "memory_budget_mib": 300000, "task_overhead_mib": 1024,
            "job_limits": {"timeout_seconds": 7200, "memory_limit_mib": 262144},
            "cost_model": {"basis": "Historical isolated 30/30 E2 wall/n <0.07 s/nt/CLI, peak <2.3 MiB/nt. Factor 2 time and 7 MiB/nt memory allow cap100 uncertainty; envelope only, not a measured cap100 prediction. One RNA per task.",
                           "seconds_per_nt_condition": .07, "seconds_constant_per_condition": 5,
                           "safety_factor": 2, "memory_base_mib": 2048, "mib_per_nt": 7}})
    arrays = {"e2_synthetic_shapes_capb_timing_20261012.json": settings,
              "e2_rnacentral_capb_timing_20261012.json": rna}
    return configs, arrays


def main():
    configs, arrays = configurations()
    for name, config in configs.items():
        path = ROOT / "configs" / name
        if path.exists():
            raise ValueError("immutable configuration already exists: " + str(path))
        write_json(path, config)
    for name, settings in arrays.items():
        path = ROOT / "configs/arrays" / name
        if path.exists():
            raise ValueError("immutable array configuration already exists: " + str(path))
        write_json(path, settings)
    for kind, config_name, array_name in (
            ("synthetic-shapes", "e2_synthetic_shapes_15aa7742_capb_timing_20261012.json", "e2_synthetic_shapes_capb_timing_20261012.json"),
            ("rnacentral", "e2_rnacentral_15aa7742_capb_timing_20261012.json", "e2_rnacentral_capb_timing_20261012.json")):
        name = "e2-" + kind + "-15aa7742-capb-timing-20261012"
        plan = prepare(ROOT / "configs" / config_name, ROOT / "configs/arrays" / array_name,
                       ROOT / "results/runs" / name, ROOT / "jobs" / name,
                       ROOT / "deps/LinearCParty-15aa7742", ROOT / "build/linear-cparty-15aa7742", dry_run=True)
        print(json.dumps({"kind": kind, "inputs": len(plan["inputs"]), "tasks": len(plan["tasks"]),
                          "estimated_node_hours": sum(t["estimated_seconds"] for t in plan["tasks"]) / 3600,
                          "estimated_entries": plan["file_budget"]["estimated_shared_entries"],
                          "dependencies_verified": plan["dependencies_verified_at_planning"]}))


if __name__ == "__main__":
    main()
