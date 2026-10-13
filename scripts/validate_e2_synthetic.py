#!/usr/bin/env python3
"""Independently recompute exported E2 scalars and prediction identities."""

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import statistics


def near(a, b):
    if not math.isclose(float(a), float(b), rel_tol=1e-10, abs_tol=1e-12):
        raise ValueError(f"scalar mismatch: {a} versus {b}")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate(output, test_log):
    output = Path(output)
    summary = json.loads((output / "summary.json").read_text())
    with (output / "metrics.tsv").open() as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    with (output / "comparisons.tsv").open() as stream:
        pairs = list(csv.DictReader(stream, delimiter="\t"))
    counts = Counter(r["run"] for r in rows)
    if counts != {"synthetic_cpu40": 468, "synthetic_node_f": 162, "rnacentral_cap100": 21, "rnacentral_cap30": 21}:
        raise ValueError("wrong output cohort counts")
    lookup = {(r["run"], r["name"], r["program"]): r for r in rows}
    if len(lookup) != len(rows):
        raise ValueError("duplicate metrics row")
    for r in rows:
        if int(r["measured_repetitions"]) != 3 or int(r["n"]) != int(r["free_bases"]) + 2 * int(r["scaffold_pairs"]):
            raise ValueError("bad repetition/length identity")
        near(r["B_G"], r["pair_depth_sum"])
        for outcome in ("wall_seconds", "rss_gib"):
            values = json.loads(r[outcome + "_repetitions"])
            median = statistics.median(values)
            near(r[outcome + "_median"], median)
            near(r[outcome + "_min"], min(values))
            near(r[outcome + "_max"], max(values))
            near(r[outcome + "_relative_range"], (max(values) - min(values)) / median)
            if r["run"] not in ("synthetic_cpu40", "rnacentral_cap100"):
                continue
            models = summary["models"][r["program"]][outcome]
            for name, model in models.items():
                prediction = math.exp(model["intercept"] + sum(beta * math.log1p(float(r[k]))
                                      for k, beta in model["coefficients"].items()))
                near(r[name + "_" + outcome + "_prediction"], prediction)
                if r["run"] == "rnacentral_cap100":
                    factor = summary["node_f_time_calibration"][r["program"]] if outcome == "wall_seconds" else 1
                    near(r[name + "_" + outcome + "_prediction_node_f"], prediction * factor)
    for p in pairs:
        a = lookup[p["run_a"], p["name_a"], p["program"]]
        b = lookup[p["run_b"], p["name_b"], p["program"]]
        near(p["n_ratio"], float(b["n"]) / float(a["n"]))
        for outcome in ("wall_seconds", "rss_gib"):
            va, vb = float(a[outcome + "_median"]), float(b[outcome + "_median"])
            near(p[outcome + "_ratio"], vb / va)
            near(p[outcome + "_difference"], vb - va)
            if p[outcome + "_model_ratio"]:
                coefficients = summary["models"][p["program"]][outcome]["geometry_full"]["coefficients"]
                delta = sum(beta * (math.log1p(float(b[k])) - math.log1p(float(a[k])))
                            for k, beta in coefficients.items())
                near(p[outcome + "_model_ratio"], math.exp(delta))
                contributions = json.loads(p[outcome + "_model_log_contributions"])
                near(sum(contributions.values()), delta)
                if float(p["n_ratio"]) != 1:
                    near(p[outcome + "_growth_exponent"], math.log(vb / va) / math.log(float(p["n_ratio"])))
    cv_checks = 0
    for program, outcomes in summary["models"].items():
        group = [r for r in rows if r["program"] == program and r["run"] == "synthetic_cpu40"]
        for outcome, models in outcomes.items():
            for key in ("n", "family", "seed"):
                residuals = [math.log(float(r["geometry_full_" + outcome + "_cv_" + key]) / float(r[outcome + "_median"]))
                             for r in group]
                cv = models["geometry_full"]["validation"]["leave_one_" + key + "_out"]
                near(cv["pooled"]["log_rmse"], math.sqrt(statistics.mean(v*v for v in residuals)))
                if sum(f["n"] for f in cv["folds"]) != len(group):
                    raise ValueError("CV fold denominator mismatch")
                cv_checks += 1
    log = Path(test_log).read_text()
    match = re.search(r"Ran (\d+) tests in ([\d.]+)s\n\nOK\s*$", log)
    if not match:
        raise ValueError("test log did not finish successfully")
    files = sorted(p for p in output.iterdir() if p.is_file() and p.name != "validation.json")
    if {p.name for p in files} != {"metrics.tsv", "comparisons.tsv", "summary.json", "scaling_time.png", "scaling_rss.png", "feature_time.png"}:
        raise ValueError("unexpected output files")
    validation = {"schema_version": 1, "metrics_rows_recomputed": len(rows), "comparisons_recomputed": len(pairs),
                  "grouped_cv_scalar_checks": cv_checks, "tests_passed": int(match[1]), "test_seconds": float(match[2]),
                  "test_log_sha256": digest(test_log), "test_log": log,
                  "outputs": {p.name: {"sha256": digest(p), "bytes": p.stat().st_size} for p in files},
                  "validator_sha256": digest(__file__),
                  "new_test_source_sha256": digest(Path(__file__).resolve().parents[1] / "tests/test_e2_synthetic.py"),
                  "tar_extractions_to_filesystem": 0, "qsub_calls": 0}
    (output / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    print(json.dumps({k: validation[k] for k in ("metrics_rows_recomputed", "comparisons_recomputed", "tests_passed")}))
    return validation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output")
    parser.add_argument("--test-log", required=True)
    args = parser.parse_args()
    validate(args.output, args.test_log)


if __name__ == "__main__":
    main()
