#!/usr/bin/env python3
"""Recompute E0/E2 exported scalars without RNA calculations or extraction."""

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import statistics

from task_storage import RunReader


def near(a, b):
    if not math.isclose(float(a), float(b), rel_tol=1e-10, abs_tol=1e-12):
        raise ValueError(f"scalar mismatch: {a} versus {b}")


def table(path):
    with Path(path).open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(2**20), b""):
            h.update(chunk)
    return h.hexdigest()


def validate_memory(output):
    summary = json.loads((output / "memory_summary.json").read_text())
    rows = table(output / "memory_metrics.tsv")
    families = table(output / "memory_families.tsv")
    lookup = {(r["run"], r["name"], r["mode"]): r for r in rows}
    if len(lookup) != len(rows):
        raise ValueError("duplicate memory rows")
    checks = 0
    for kind, provenance in summary["provenance"].items():
        reader = RunReader(provenance["path"])
        if reader.manifest_sha256 != provenance["manifest_sha256"]:
            raise ValueError("memory manifest changed")
        label = "memory_" + kind
        if len([r for r in rows if r["run"] == label]) != len(reader.manifest["jobs"]):
            raise ValueError("memory denominator mismatch")
        for job in reader.manifest["jobs"]:
            row = lookup[label, job["input"], job["mode"]]
            raw = json.loads(reader.artifact(job, "stdout.txt"))
            stage = raw["runs"][0]["stages"][0]
            near(row["rss_gib"], job["peak_rss_kib"] / 2**20)
            near(row["wall_seconds"], job["time_wall_seconds"])
            near(row["driver_vmhwm_kib"], stage["rss_after"]["vmhwm_kib"])
            checks += 3
            work = stage["summary"].get("work", {})
            if work:
                near(row["rss_bytes_per_item_visit"], job["peak_rss_kib"] * 1024 / work["item_visits"])
                checks += 1
            elif row["rss_bytes_per_item_visit"] != "NA":
                raise ValueError("missing profile work must stay NA")
            for k, value in work.items():
                if k == "families":
                    selected = [f for f in families if all(f[key] == row[key] for key in ("run", "name", "mode"))]
                    if len(selected) != len(value):
                        raise ValueError("family denominator mismatch")
                    for expected in value:
                        f = next(f for f in selected if f["nonterminal"] == expected["nonterminal"])
                        for field in expected.keys() - {"nonterminal"}:
                            near(f[field], expected[field])
                            checks += 1
                else:
                    near(row[k], value)
                    checks += 1
    for p in table(output / "memory_comparisons.tsv"):
        a, b = [lookup[p["run_" + side], p["name_" + side], p["mode"]] for side in ("a", "b")]
        for key, value in p.items():
            if not key.endswith("_ratio") or not value:
                continue
            k = key[:-6]
            if value == "NA":
                if a.get(k) not in (None, "", "NA", "0", "0.0"):
                    raise ValueError("unexpected NA ratio")
            else:
                near(value, float(b[k]) / float(a[k]))
                checks += 1
    return checks


def validate(output, e0_output, test_log):
    output = Path(output)
    e0_output = Path(e0_output) if e0_output else None
    summary = json.loads((output / "summary.json").read_text())
    rows = table(output / "metrics.tsv")
    lookup = {(r["run"], r["name"], r["program"]): r for r in rows}
    if len(lookup) != len(rows):
        raise ValueError("duplicate exported metrics")
    for label, provenance in summary["provenance"]["runs"].items():
        reader = RunReader(provenance["path"])
        if reader.manifest_sha256 != provenance["manifest_sha256"]:
            raise ValueError("run manifest changed since analysis")
        jobs = defaultdict(list)
        for job in reader.manifest["jobs"]:
            if not job["warmup"]:
                jobs[job["input"], job["program"]].append(job)
        selected = [r for r in rows if r["run"] == label]
        if len(selected) != len(jobs):
            raise ValueError("raw/exported input/CLI denominator mismatch")
        for r in selected:
            original = sorted(jobs[r["name"], r["program"]], key=lambda j: j["repetition"])
            if len(original) != 3 or any(j["status"] != "success" for j in original):
                raise ValueError("three successful measurements required")
            for outcome, key, scale in (("wall_seconds", "time_wall_seconds", 1), ("rss_gib", "peak_rss_kib", 2**20)):
                values = [j[key] / scale for j in original]
                for a, b in zip(json.loads(r[outcome + "_repetitions"]), values):
                    near(a, b)
                mid = statistics.median(values)
                for suffix, value in (("median", mid), ("min", min(values)), ("max", max(values)),
                                      ("relative_range", (max(values) - min(values)) / mid)):
                    near(r[outcome + "_" + suffix], value)
    pairs = table(output / "comparisons.tsv")
    for p in pairs:
        a, b = [lookup[p["run_" + side], p["name_" + side], p["program"]] for side in ("a", "b")]
        near(p["n_ratio"], float(b["n"]) / float(a["n"]))
        for outcome in ("wall_seconds", "rss_gib"):
            va, vb = float(a[outcome + "_median"]), float(b[outcome + "_median"])
            near(p[outcome + "_ratio"], vb / va)
            near(p[outcome + "_difference"], vb - va)
            if p[outcome + "_model_ratio"]:
                model = summary["models"][p["program"]][outcome]["geometry_full"]
                delta = sum(beta * (math.log1p(float(b[k])) - math.log1p(float(a[k])))
                            for k, beta in model["coefficients"].items())
                near(p[outcome + "_model_ratio"], math.exp(delta))
                near(sum(json.loads(p[outcome + "_model_log_contributions"]).values()), delta)
                if float(p["n_ratio"]) != 1:
                    near(p[outcome + "_growth_exponent"], math.log(vb / va) / math.log(float(p["n_ratio"])))
    for item in summary.get("comparison_distributions", []):
        values = [float(p[item["outcome"] + "_ratio"]) for p in pairs
                  if p["kind"] == item["kind"] and p["program"] == item["program"]]
        if len(values) != item["count"]:
            raise ValueError("comparison distribution denominator mismatch")
        for key, value in (("min", min(values)), ("median", statistics.median(values)), ("max", max(values))):
            near(item[key], value)
    slopes = table(output / "slopes.tsv")
    for slope in slopes:
        selected = [r for r in rows if r["run"] == "synthetic_current" and all(
            r[k] == slope[k] for k in ("program", "series", "seed"))]
        x, y = [math.log(float(r["n"])) for r in selected], [math.log(float(r[slope["outcome"] + "_median"])) for r in selected]
        mx, my = statistics.mean(x), statistics.mean(y)
        beta = sum((a-mx)*(b-my) for a, b in zip(x, y)) / sum((a-mx)**2 for a in x)
        near(slope["slope_ln_y_vs_ln_n"], beta)
    cv_checks = prediction_checks = 0
    for program, outcomes in summary["models"].items():
        group = [r for r in rows if r["program"] == program and r["run"] == "synthetic_current"]
        for outcome, models in outcomes.items():
            for name, model in models.items():
                for row in [r for r in rows if r["program"] == program and r["run"] in ("synthetic_current", "rnacentral_current")]:
                    prediction = math.exp(model["intercept"] + sum(beta * math.log1p(float(row[k]))
                                          for k, beta in model["coefficients"].items()))
                    near(row[name + "_" + outcome + "_prediction"], prediction)
                    prediction_checks += 1
            for key in ("n", "family", "seed"):
                residuals = [math.log(float(r["geometry_full_" + outcome + "_cv_" + key]) / float(r[outcome + "_median"])) for r in group]
                cv = models["geometry_full"]["validation"]["leave_one_" + key + "_out"]
                near(cv["pooled"]["log_rmse"], math.sqrt(statistics.mean(v*v for v in residuals)))
                if sum(f["n"] for f in cv["folds"]) != len(group):
                    raise ValueError("CV denominator mismatch")
                cv_checks += 1
    e0 = []
    if e0_output is not None:
        e0 = table(e0_output / "comparisons.tsv")
        e0_analysis = json.loads((e0_output / "analysis.json").read_text())
        expected = {}
        for row in e0_analysis["observations"]:
            expected[row["input"], row["comparison_kind"]] = row["metrics"]
            if "default_minus_compatible" in row:
                expected[row["input"], "weak_closure_vs_compatible"] = row["default_minus_compatible"]
        if {(r["input"], r["comparison"]) for r in e0} != set(expected) or len(e0) != len(expected):
            raise ValueError("E0 exported comparison denominator mismatch")
        for row in e0:
            for key, value in expected[row["input"], row["comparison"]].items():
                near(row[key], value)
        for row in table(e0_output / "distributions.tsv"):
            values = sorted((abs(float(r[row["quantity"]])) if row["absolute"] == "True" else float(r[row["quantity"]]))
                            for r in e0 if r["comparison"] == row["comparison"] and
                            (row["cohort"] == "all" or row["cohort"] == r["cohort"]))
            if len(values) != int(row["count"]):
                raise ValueError("E0 distribution denominator mismatch")
            exceeded = sum(abs(float(r["delta_log_z"])) > e0_analysis["tolerances"]["log_z_absolute_tolerance"] or
                           abs(float(r["delta_free_energy"])) > e0_analysis["tolerances"]["free_energy_absolute_tolerance"]
                           for r in e0 if r["comparison"] == row["comparison"] and
                           (row["cohort"] == "all" or row["cohort"] == r["cohort"]))
            if int(row["exceeds_either_tolerance"]) != exceeded:
                raise ValueError("E0 tolerance count mismatch")
            for key, q in (("min", 0), ("p10", .1), ("median", .5), ("p90", .9), ("p95", .95), ("max", 1)):
                pos = (len(values) - 1) * q
                lo, hi = math.floor(pos), math.ceil(pos)
                near(row[key], values[lo] + (values[hi] - values[lo]) * (pos - lo))
    memory_checks = validate_memory(output) if (output / "memory_summary.json").exists() else 0
    log = Path(test_log).read_text()
    match = re.search(r"Ran (\d+) tests in ([\d.]+)s\n\nOK \(skipped=(\d+)\)\s*$", log)
    if not match:
        raise ValueError("successful test log with explicit RNA CLI skips required")
    files = [p for directory in (e0_output, output) if directory is not None for p in directory.iterdir() if p.is_file() and p.name != "validation.json"]
    report = {"metrics_rows_recomputed": len(rows), "comparisons_recomputed": len(pairs), "slopes_recomputed": len(slopes),
        "e0_rows_recomputed": len(e0), "memory_scalar_checks": memory_checks, "prediction_checks": prediction_checks, "grouped_cv_scalar_checks": cv_checks,
        "test_cases": int(match[1]), "tests_passed": int(match[1]) - int(match[3]), "tests_skipped": int(match[3]),
        "test_log": log, "test_log_sha256": digest(test_log), "validator_sha256": digest(__file__),
        "metrics_by_run": dict(Counter(r["run"] for r in rows)),
        "outputs": {str(p): {"bytes": p.stat().st_size, "sha256": digest(p)} for p in files},
        "tar_extractions_to_filesystem": 0}
    (output / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k.endswith(("recomputed", "checks", "passed", "skipped"))}))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output")
    parser.add_argument("--e0-output")
    parser.add_argument("--test-log", required=True)
    args = parser.parse_args()
    validate(args.output, args.e0_output, args.test_log)
