#!/usr/bin/env python3
"""Compare exact inside log Z / PF free energy and count closure differences."""

import argparse
from collections import Counter
import json
import math
from pathlib import Path

from rna_records import read_record
from run_cli import digest, experiment_path, now, write_json
from task_storage import RunReader, MemoryFile, input_records

QUANTITIES = ("log_z", "free_energy")
POLICIES = {"E0LegacyCompatible": "legacy-compatible-wmb", "E0Default": "weakly-closed-wmb",
            "E0CPartyOracle": "legacy-compatible-wmb"}


def analysis_settings():
    return {
        "quantities": list(QUANTITIES), "excluded_quantities": ["bpp", "mfe", "contexts"],
        "comparisons": {
            "compatible_vs_oracle": "Exact inside log Z / PF free energy only; single upstream 800e95a + patches 0001,0002,0006,0007 oracle without closure",
            "default_vs_oracle": "Reference model difference; also report default minus compatible, attributing closure only when compatible parity passes"},
        "log_z_absolute_tolerance": 1e-9, "free_energy_absolute_tolerance": 1e-10,
        "tolerance_source": "LinearCParty 2456f9b2 test/tools/audit_legacy_inside_parity.py defaults",
        "tolerance_scope": "Absolute roundoff allowance for full-precision inside values; no relative tolerance or printed-CLI tolerance",
        "closure_count_rule": "Count differences exceeding either PF tolerance only when compatible-versus-oracle passes both tolerances"}


def read_observation(path, record=None):
    value = json.loads(path.read_text())
    for key in (*QUANTITIES, "rt"):
        if isinstance(value[key], bool) or not math.isfinite(value[key]):
            raise ValueError("nonfinite or invalid E0 observation")
    if value["rt"] <= 0 or abs(value["log_z"] + value["free_energy"] / value["rt"]) > 1e-9:
        raise ValueError("thermodynamic metadata inconsistent")
    # Ignore historical BPP/MFE/profile fields even when an old artifact has them.
    return {k: value[k] for k in (*QUANTITIES, "rt", "policy") if k in value}


def compare(candidate, reference):
    if abs(candidate["rt"] - reference["rt"]) > 1e-12:
        raise ValueError("RT differs")
    return {"delta_log_z": candidate["log_z"] - reference["log_z"],
            "delta_free_energy": candidate["free_energy"] - reference["free_energy"]}


def within_tolerance(metrics, config):
    return (abs(metrics["delta_log_z"]) <= config["log_z_absolute_tolerance"] and
            abs(metrics["delta_free_energy"]) <= config["free_energy_absolute_tolerance"])


def analyze(run):
    reader = RunReader(run)
    manifest = reader.manifest
    if manifest["status"] == "running" or manifest.get("summary_only"):
        raise ValueError("finish and aggregate observations before E0 analysis")
    records = input_records(reader)
    seen_inputs = set()
    for entry in manifest["inputs"]:
        if entry["name"] in seen_inputs:
            raise ValueError("duplicate E0 input")
        seen_inputs.add(entry["name"])
    config = manifest["config"]["analysis"]
    for key in ("log_z_absolute_tolerance", "free_energy_absolute_tolerance"):
        if not math.isfinite(config[key]) or config[key] < 0:
            raise ValueError("finite nonnegative tolerances required")
    parsed, errors, seen = {}, {}, set()
    jobs = [j for j in manifest["jobs"] if not j["warmup"]]
    for job in jobs:
        key = (job["input"], job["program"])
        if key in seen:
            raise ValueError("E0 expects one exact observation per input/program")
        seen.add(key)
        if job["status"] != "success":
            continue
        try:
            value = read_observation(MemoryFile(reader.artifact(job, "stdout.txt")), records[job["input"]])
            if value.get("policy") != POLICIES[job["program"]]:
                raise ValueError("E0 closure policy differs")
            parsed[key] = value
        except (ValueError, KeyError, TypeError, OSError) as exc:
            errors[key] = str(exc)
    rows = []
    for name in records:
        for program, kind in (("E0LegacyCompatible", "numerical_parity"), ("E0Default", "closure_model_difference")):
            keys = [(name, program), (name, "E0CPartyOracle")]
            row = {"input": name, "comparison_kind": kind, "program": program, "reference_program": "E0CPartyOracle"}
            if any(k in errors for k in keys):
                row.update(status="artifact_error", errors=[errors[k] for k in keys if k in errors])
            elif any(k not in parsed for k in keys):
                row["status"] = "observation_unavailable"
            else:
                try:
                    metrics = compare(parsed[keys[0]], parsed[keys[1]])
                    row.update(metrics=metrics, status="success" if kind == "closure_model_difference" or
                               within_tolerance(metrics, config) else "parity_failed")
                    if kind == "closure_model_difference":
                        row["restriction_direction_consistent"] = (metrics["delta_log_z"] <= config["log_z_absolute_tolerance"] and
                            metrics["delta_free_energy"] >= -config["free_energy_absolute_tolerance"])
                        row["closure_attribution_valid"] = False
                        compatible = parsed.get((name, "E0LegacyCompatible"))
                        if compatible is not None:
                            baseline = compare(compatible, parsed[keys[1]])
                            effect = compare(parsed[keys[0]], compatible)
                            row.update(default_minus_compatible=effect, compatible_minus_oracle=baseline,
                                       closure_attribution_valid=within_tolerance(baseline, config))
                            row["closure_effect_detected"] = row["closure_attribution_valid"] and not within_tolerance(effect, config)
                except ValueError as exc:
                    row.update(status="artifact_error", errors=[str(exc)])
            rows.append(row)
    parity = [r for r in rows if r["comparison_kind"] == "numerical_parity"]
    closure = [r for r in rows if r["comparison_kind"] == "closure_model_difference"]
    eligible = [r for r in closure if r["status"] == "success" and r.get("closure_attribution_valid")]
    return {"schema_version": 2, "created_utc": now(), "run_manifest_sha256": reader.manifest_sha256,
        "analysis_script_sha256": digest(Path(__file__)), "planned_inputs": len(records),
        "comparison_quantities": list(QUANTITIES), "excluded_quantities": ["bpp", "mfe", "contexts"],
        "tolerances": {k: config[k] for k in ("log_z_absolute_tolerance", "free_energy_absolute_tolerance")},
        "tolerance_scope": "Absolute roundoff allowance from implementation audit_legacy_inside_parity.py; no relative tolerance",
        "job_status_counts": dict(Counter(j["status"] for j in jobs)),
        "numerical_parity_counts": dict(Counter(r["status"] for r in parity)),
        "closure_counts": {"planned_inputs": len(records), "attribution_valid": len(eligible),
            "effect_detected": sum(r["closure_effect_detected"] for r in eligible),
            "within_roundoff": sum(not r["closure_effect_detected"] for r in eligible),
            "unassessable": len(records) - len(eligible),
            "restriction_direction_inconsistent": sum(not r["restriction_direction_consistent"] for r in eligible)},
        "observations": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run"); parser.add_argument("--out", required=True)
    args = parser.parse_args(); out = experiment_path(args.out)
    if out.exists(): parser.error("use new analysis output")
    result = analyze(experiment_path(args.run)); out.parent.mkdir(parents=True, exist_ok=True); write_json(out, result)
    print(json.dumps({"parity": result["numerical_parity_counts"], "closure": result["closure_counts"]}))
    return int(any(r["status"] != "success" for r in result["observations"]) or
               result["closure_counts"]["restriction_direction_inconsistent"] > 0)


if __name__ == "__main__": raise SystemExit(main())
