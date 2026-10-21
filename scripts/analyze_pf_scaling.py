#!/usr/bin/env python3
"""Summarize the matched PF array, retaining failures and measured denominators."""
import argparse
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
from statistics import median

from plan_array import load_plan, task_jobs
from run_cli import digest, mode_applies, now, write_json
from task_storage import RunReader

METHODS = ("cparty", "exact", "beam100")


def value(text):
    z = json.loads(text)
    if z.get("policy") != "legacy-compatible-wmb":
        raise ValueError("closure policy differs")
    for key in ("log_z", "free_energy", "rt"):
        if isinstance(z[key], bool) or not isinstance(z[key], (int, float)) or not math.isfinite(z[key]):
            raise ValueError("nonfinite or missing PF observation")
    if z["rt"] <= 0 or abs(z["log_z"] + z["free_energy"] / z["rt"]) > 1e-9:
        raise ValueError("inconsistent PF thermodynamic metadata")
    return z


def summarize(observations, inputs, modes, repetitions=3):
    groups = defaultdict(list)
    for row in observations:
        if not row["warmup"]:
            groups[row["input"], row["mode"]].append(row)
    rows, parity = [], []
    for item in inputs:
        per_method = {}
        for mode in modes:
            n = item["length"]
            if not mode_applies(mode, item):
                rows.append({"input": item["name"], "length": n, "family": item["family"],
                             "seed": item["seed"], "method": mode["name"], "status": "excluded_by_plan",
                             "planned_repetitions": 0, "successful_repetitions": 0})
                continue
            jobs = groups[item["name"], mode["name"]]
            if len(jobs) != repetitions or len({j["repetition"] for j in jobs}) != repetitions:
                raise ValueError("missing or duplicated planned repetition: " + item["name"])
            good = [j for j in jobs if j["status"] == "success" and j.get("value") is not None]
            row = {"input": item["name"], "length": n, "family": item["family"], "seed": item["seed"],
                   "method": mode["name"], "planned_repetitions": repetitions,
                   "successful_repetitions": len(good), "statuses": dict(Counter(j["status"] for j in jobs)),
                   "status": "success" if len(good) == repetitions else "unavailable"}
            if len(good) == repetitions:
                scalars = [j["value"] for j in good]
                if any(abs(z["log_z"] - scalars[0]["log_z"]) > 1e-9 for z in scalars):
                    raise ValueError("PF changed between repetitions")
                for field in ("elapsed_seconds", "time_wall_seconds", "peak_rss_kib"):
                    vals = [j[field] for j in good]
                    if any(v is None or v < 0 or (v == 0 and field != "time_wall_seconds") for v in vals):
                        raise ValueError("missing/invalid successful resource measurement")
                    row[field] = median(vals)
                    row[field + "_min"] = min(vals)
                    row[field + "_max"] = max(vals)
                row.update({k: scalars[0][k] for k in ("log_z", "free_energy", "rt")})
                per_method[mode["name"]] = row
            rows.append(row)
        pair = {"input": item["name"], "length": item["length"], "status": "unavailable"}
        if all(m in per_method for m in ("cparty", "exact")):
            a, b = per_method["cparty"], per_method["exact"]
            pair.update(delta_log_z=b["log_z"] - a["log_z"],
                        delta_free_energy=b["free_energy"] - a["free_energy"])
            pair["status"] = "pass" if (abs(pair["delta_log_z"]) <= 1e-9 and
                                        abs(pair["delta_free_energy"]) <= 1e-10 and
                                        abs(b["rt"] - a["rt"]) <= 1e-12) else "parity_failed"
        parity.append(pair)
    # Plot only a common cohort of successful, PF-compatible input records at
    # each length. Do not give each method a different success-only denominator.
    good_keys = {(r["input"], r["method"]) for r in rows if r["status"] == "success"}
    passing = {r["input"] for r in parity if r["status"] == "pass"}
    paired = []
    for item in inputs:
        applicable = [m["name"] for m in modes if mode_applies(m, item)]
        # Long lengths may have only CParty and beam100. Their finite outputs
        # remain comparable, but exact parity has not been assessed there.
        parity_ok = item["name"] in passing if "exact" in applicable else True
        if parity_ok and all((item["name"], m) in good_keys for m in applicable):
            paired.extend(r for r in rows if r["input"] == item["name"] and r["status"] == "success")
    return rows, parity, paired


def crossings(rows):
    lookup = {(r["input"], r["method"]): r for r in rows}
    cohorts = defaultdict(list)
    for r in rows:
        if r["method"] != "cparty":
            continue
        for other in ("exact", "beam100"):
            b = lookup.get((r["input"], other))
            if b:
                cohorts[r["family"], r["seed"], other].append((r, b))
    reports = []
    for (family, seed, other), points in sorted(cohorts.items()):
        points.sort(key=lambda pair: pair[0]["length"])
        for metric in ("elapsed_seconds", "peak_rss_kib"):
            brackets = []
            for (a, b), (c, d) in zip(points, points[1:]):
                # Only neighboring planned doubling lengths can bracket a crossing.
                if c["length"] != 2 * a["length"]:
                    continue
                y0, y1 = math.log(a[metric] / b[metric]), math.log(c[metric] / d[metric])
                if y0 * y1 <= 0 and y0 != y1:
                    x = math.exp(math.log(a["length"]) - y0 * math.log(c["length"] / a["length"]) / (y1 - y0))
                    brackets.append({"lower_nt": a["length"], "upper_nt": c["length"], "interpolated_nt": x})
            reports.append({"family": family, "seed": seed, "other": other, "metric": metric,
                            "status": "bracketed" if brackets else "not_observed_in_valid_measured_range",
                            "brackets": brackets, "valid_lengths": [a["length"] for a, _ in points]})
    return reports


def tsv(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, keys, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow({k: json.dumps(v, sort_keys=True) if isinstance(v, (list, dict)) else v for k, v in r.items()})


def matched_cparty_modes(plan, reference):
    """Reuse only CParty observations, with the same frozen inputs and model."""
    fields = ("name", "sha256", "length", "family", "seed")
    inputs = lambda p: sorted(tuple(x[k] for k in fields) for x in p["inputs"])
    if inputs(plan) != inputs(reference):
        raise ValueError("CParty reference inputs differ")
    for key in ("model_settings", "warmups", "repetitions"):
        if plan["config"][key] != reference["config"][key]:
            raise ValueError("CParty reference model or repetitions differ")
    if plan["software_lock"]["energy_parameter_sha256"] != reference["software_lock"]["energy_parameter_sha256"]:
        raise ValueError("CParty reference parameters differ")
    if {m["name"] for m in plan["config"]["modes"]} != {"exact", "beam100"}:
        raise ValueError("new plan must contain only exact and beam100")
    modes = [m for m in reference["config"]["modes"] if m["name"] == "cparty"]
    if len(modes) != 1:
        raise ValueError("reference must contain a single CParty mode")
    return modes + plan["config"]["modes"]


def read_observations(plan, selected_modes=None):
    observations, attempts = [], []
    for task in plan["tasks"]:
        directory = Path(plan["output_dir"]) / task["directory"]
        metadata = json.loads((directory / "task.json").read_text())
        if metadata["status"] not in ("success", "failed"):
            raise ValueError("task unfinished; preserve missing results and analyze after completion")
        attempts.extend(metadata["attempts"])
        if len({a["hostname"] for a in metadata["attempts"]}) > 1:
            raise ValueError("task moved between nodes; review timing attempts before pooling")
        reader = RunReader(directory / "run")
        manifest = reader.manifest
        if (manifest.get("software_changed_during_run") is not False or
                manifest.get("binaries_changed_during_run") is not False or
                manifest.get("software_lock_valid_after_run") is not True):
            raise ValueError("post-run dependency verification failed or unavailable")
        expected = {j["id"] for j in task_jobs(plan, task)}
        if {j["id"] for j in manifest["jobs"]} != expected:
            raise ValueError("task schedule differs from frozen plan")
        for job in manifest["jobs"]:
            if selected_modes is not None and job["mode"] not in selected_modes:
                continue
            row = dict(job)
            if row["status"] == "success":
                try:
                    row["value"] = value(reader.artifact(job, "stdout.txt"))
                except (ValueError, KeyError, TypeError) as exc:
                    row.update(status="invalid_pf", parse_error=str(exc))
            observations.append(row)
    return observations, attempts


def analyze(plan_path, output, cparty_plan_path=None):
    if output.exists():
        raise ValueError("use a new analysis output directory")
    plan = load_plan(plan_path)
    observations, attempts = read_observations(plan)
    modes = plan["config"]["modes"]
    reference_info = None
    if cparty_plan_path is not None:
        reference = load_plan(cparty_plan_path)
        modes = matched_cparty_modes(plan, reference)
        previous, previous_attempts = read_observations(reference, {"cparty"})
        observations.extend(previous)
        reference_info = {"plan": str(cparty_plan_path), "plan_sha256": digest(cparty_plan_path),
                          "source_commit": reference["software_lock"]["source_commit"],
                          "reused_modes": ["cparty"], "observations": len(previous),
                          "attempts": previous_attempts}
    rows, parity, paired = summarize(observations, plan["inputs"], modes, plan["config"]["repetitions"])
    output.mkdir(parents=True)
    tsv(output / "resources.tsv", rows)
    tsv(output / "parity.tsv", parity)
    report = {"created_utc": now(), "plan_sha256": digest(plan_path), "script_sha256": digest(Path(__file__)),
              "planned_inputs": len(plan["inputs"]), "observations": observations, "resources": rows,
              "parity": parity, "paired_resources": paired, "crossings": crossings(paired),
              "status_counts": dict(Counter(r["status"] for r in rows)), "attempts": attempts,
              "source_commit": plan["software_lock"]["source_commit"], "cparty_reference": reference_info,
              "pairing_rule": "Same input and all applicable methods, three successful repetitions; CParty/exact parity required where exact is measured. Exact-excluded lengths use paired CParty/beam100 finite outputs, with exact parity unassessed.",
              "time_definition": "Median of three process elapsed_seconds including GNU time wrapper; warmups excluded",
              "memory_definition": "Median peak RSS of three processes; min/max also saved; RLIMIT_AS is separate"}
    write_json(output / "summary.json", report)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--cparty-plan", type=Path, help="reuse only CParty records from a matched completed plan")
    args = parser.parse_args()
    report = analyze(args.plan.resolve(), args.out.resolve(), args.cparty_plan.resolve() if args.cparty_plan else None)
    print(json.dumps({"status_counts": report["status_counts"], "paired_rows": len(report["paired_resources"])}))
