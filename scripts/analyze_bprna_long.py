#!/usr/bin/env python3
"""Packed E1 beam-width sensitivity, using b=200 (never exact) as reference.

Read one RNA at a time through tar member streams, without extraction or CLI
execution. Write two TSVs, two JSONs and at most two PNGs. RNA is the unit of
observation; preserve signed F differences, missing windows and zero events.
"""

import argparse
from collections import Counter, defaultdict
import csv
import gc
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

import analyze_accuracy as accuracy
from analyze_bprna_short import ACCESS, CHECKS, CONTEXT, monotonicity, plan_provenance
from analyze_e2_shape import ranks
from bprna_capb import check_task, resource_fields
from run_cli import digest, git_state, now
from task_storage import FORMAT, RunReader, input_records

BEAMS = (20, 50, 100, 200)
MODES = ("beam20_default", "beam50_capb", "beam100_capb", "beam200_capb")
TIERS = (("nt400_799", 400, 799), ("nt800_1599", 800, 1599),
         ("nt1600_3199", 1600, 3199), ("nt3200_4381", 3200, 4381))
FEATURES = ("scaffold_pairs", "free_bases", "paired_fraction", "gc_fraction", "B_G",
            "pair_depth_max", "free_runs", "free_run_max", "loop_free_max",
            "case3_ancestors_max", "case3_ancestors_mean", "max_pair_children",
            "cross_loop_canonical_pairs_turn", "same_loop_pairs_turn")
CORE = ("delta_F_vs_b200_kcal_mol", "abs_delta_F_vs_b200_kcal_mol",
        "delta_F_vs_b200_per_nt", "log_Z_ratio_vs_b200", "Z_ratio_vs_b200",
        "generated_pair_l1_sum", "base_state_tv_mean_Gfree", "unpaired_mae_Gfree",
        "context_tv_mean_Gfree")
RESOURCES = tuple(f"{p}_{k}" for p in sorted(accuracy.PROGRAMS) for k in
                  ("elapsed_seconds", "time_wall_seconds", "peak_rss_kib"))
IDENTITY = ("input", "run", "tier", "task", "mode", "beam", "local_cap", "length",
            "input_sha256", "sequence_sha256", "scaffold_sha256")
THERMO = ("rt_kcal_mol", "free_energy_kcal_mol", "reference_b200_free_energy_kcal_mol")
FIELDS = list(IDENTITY + FEATURES + CORE) + CONTEXT + ACCESS + CHECKS + list(THERMO + RESOURCES)
METRICS = list(CORE) + CONTEXT + [k for k in ACCESS if "mae" in k or "max_error" in k] + list(RESOURCES)


def load_features(path):
    result = {}
    with path.open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["cohort"] != "bprna":
                continue
            if row["name"] in result:
                raise ValueError("duplicate bpRNA feature row")
            result[row["name"]] = row
    return result


def join_features(record, features):
    name, sequence, scaffold = record
    row = features[name]
    hashes = {"sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
              "scaffold_sha256": hashlib.sha256(scaffold.encode()).hexdigest()}
    if any(row[k] != v for k, v in hashes.items()) or int(row["n"]) != len(sequence):
        raise ValueError("feature S/G hash or length mismatch: " + name)
    if int(row["free_bases"]) != scaffold.count("."):
        raise ValueError("feature free-base count mismatch")
    return {**hashes, **{k: float(row[k]) for k in FEATURES}}


def check_command(job, config):
    """Require the complete planned algorithm arguments, including output flags."""
    mode = next(m for m in config["modes"] if m["name"] == job["mode"])
    program = next(p for p in config["programs"] if p["name"] == job["program"])
    expected = mode["args"] + program.get("args", [])
    if accuracy.algorithm_arguments(job["command"]) != expected:
        raise ValueError("command differs from planned condition: " + job["id"])
    beam = BEAMS[MODES.index(job["mode"])]
    for flag, value in (("-b", beam), ("--max-hairpin", max(beam, 30)),
                        ("--max-case3-depth", max(beam, 30))):
        args = job["command"]
        if args.count(flag) != 1 or args[args.index(flag) + 1] != str(value):
            raise ValueError("condition violates max(b,30) cap rule")
    if job["program"] == "LinearCParty" and job["command"].count("--bpp") != 1:
        raise ValueError("BPP export is required")


def audit_job(reader, job):
    for name in ("command.json", "result.json", "time.tsv"):
        if not job.get("artifacts", {}).get(job["directory"] + "/" + name, {}).get("sha256"):
            raise ValueError("packed artifact has no checksum")
    if json.loads(reader.artifact(job, "command.json"))["argv"] != job["command"]:
        raise ValueError("raw argv differs from checkpoint")
    saved = json.loads(reader.artifact(job, "result.json"))
    if any(saved[k] != job[k] for k in saved if k != "directory"):
        raise ValueError("raw result differs from checkpoint")
    wall, rss = reader.artifact(job, "time.tsv").decode().strip().split("\t")
    if float(wall) != job["time_wall_seconds"] or int(rss) != job["peak_rss_kib"]:
        raise ValueError("raw resource observation differs from checkpoint")
    artifact = "bpp.tsv" if job["program"] == "LinearCParty" else "stdout.txt"
    if not job["artifacts"].get(job["directory"] + "/" + artifact, {}).get("sha256"):
        raise ValueError("probability artifact has no checksum")


def compare_condition(parsed, reference, record):
    pair = accuracy.compare_bpp(parsed["LinearCParty"], reference["LinearCParty"], record)
    context = accuracy.compare_contexts(parsed["LinearPKCapR"], reference["LinearPKCapR"], record)
    access = accuracy.compare_accessibility(parsed["LinearPKRaccess"], reference["LinearPKRaccess"], record)
    delta = access["delta_free_energy_kcal_mol"]
    row = {"delta_F_vs_b200_kcal_mol": delta, "abs_delta_F_vs_b200_kcal_mol": abs(delta),
           "delta_F_vs_b200_per_nt": access["delta_free_energy_per_nt"],
           "log_Z_ratio_vs_b200": access["delta_log_z"],
           "Z_ratio_vs_b200": access["partition_function_retention"],
           **{k: pair[k] for k in ("generated_pair_l1_sum", "base_state_tv_mean_Gfree", "unpaired_mae_Gfree")},
           "context_tv_mean_Gfree": context["context_tv_mean_Gfree"],
           **{"context_" + c + "_mae_Gfree": v for c, v in context["context_mae_Gfree"].items()},
           "rt_kcal_mol": access["rt_kcal_mol"],
           "free_energy_kcal_mol": parsed["LinearPKRaccess"]["free_energy"],
           "reference_b200_free_energy_kcal_mol": reference["LinearPKRaccess"]["free_energy"]}
    for w, values in access["accessibility"].items():
        row.update({f"access_w{w}_{k}": v for k, v in values.items() if f"access_w{w}_{k}" in ACCESS})
    return row


def write_metrics(run_specs, features_path, output):
    start = time.monotonic()
    features = load_features(features_path)
    plans = [(run, tiers, json.loads((run / "array_plan.json").read_text())) for run, tiers in run_specs]
    baseline = plans[0][2]
    for run, tiers, plan in plans:
        if plan.get("output_format") != FORMAT or accuracy.model_signature(plan) != accuracy.model_signature(baseline):
            raise ValueError("run format, source, model or binary mismatch")
        if tuple(m["name"] for m in plan["config"]["modes"]) != MODES:
            raise ValueError("requires the four adopted beam/cap conditions, without exact")
        if digest(run / "inputs.tar") != plan["inputs_archive_sha256"]:
            raise ValueError("frozen input archive checksum mismatch")
    output.mkdir(parents=True, exist_ok=False)
    provenance = {"schema_version": 1, "generated_utc": now(), "reference": "beam200_capb; not exact",
                  "scope": "400 <= n <= 4381; one RNA is one observation",
                  "cap_rule": "h=d=max(b,30)", "feature_path": str(features_path),
                  "feature_sha256": digest(features_path), "analysis_git_state": git_state(accuracy.ROOT),
                  "script_sha256": {p.name: digest(p) for p in (Path(__file__), Path(accuracy.__file__),
                      Path(__file__).with_name("task_storage.py"), Path(__file__).with_name("bprna_capb.py"),
                      Path(__file__).with_name("analyze_bprna_short.py"), Path(__file__).with_name("analyze_e2_shape.py"))},
                  "runs": [], "tasks": [], "storage": "tar member streams; no extraction; one RNA in memory",
                  "resource_scope": "accuracy arrays; workers=1 within task, concurrent array tasks; one repetition; not isolated timing"}
    seen, statuses, checks, hosts, workers = set(), Counter(), Counter(), Counter(), Counter()
    check_maxima = dict.fromkeys(CHECKS, 0.)
    with (output / "metrics.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for run, selected_tiers, plan in plans:
            tasks = [t for t in plan["tasks"] if t["tier"] in selected_tiers]
            info = {**plan_provenance(run, plan), "selected_tiers": list(selected_tiers),
                    "excluded_planned_tasks": dict(Counter(t["tier"] for t in plan["tasks"] if t not in tasks)),
                    "selected_tasks": len(tasks)}
            provenance["runs"].append(info)
            for number, task in enumerate(tasks, 1):
                reader = RunReader(run / task["directory"] / "run")
                hashes, jobs = check_task(reader, plan, task)
                if reader.manifest["status"] != "success" or any(j["status"] != "success" for j in jobs.values()):
                    raise ValueError("selected task has incomplete or failed jobs")
                records = input_records(reader)
                if set(jobs) != {(name, mode, program) for name in records for mode in MODES for program in accuracy.PROGRAMS}:
                    raise ValueError("selected task has incomplete or extra observations")
                statuses.update(j["status"] for j in jobs.values())
                host = reader.manifest["host"]["hostname"]
                hosts[host] += 1
                workers[reader.manifest["workers"]] += 1
                provenance["tasks"].append({"run": run.name, "tier": task["tier"], "task": task["task_id"],
                    "run_manifest_sha256": reader.manifest_sha256,
                    "task_results_sha256": digest(reader.path / "task_results.jsonl"),
                    "inputs": len(records), "jobs": len(jobs), "host": host,
                    "resume_hosts": [h["host"]["hostname"] for h in reader.manifest.get("resume_history", [])]})
                for name, record in records.items():
                    lo, hi = next((lo, hi) for tier, lo, hi in TIERS if tier == task["tier"])
                    if name in seen or not lo <= len(record[1]) <= hi:
                        raise ValueError("duplicate RNA or length outside selected tier")
                    seen.add(name)
                    desc = join_features(record, features)
                    # Read widest first, then one candidate at a time.
                    reference = None
                    rows = {}
                    for mode in reversed(MODES):
                        parsed = {}
                        for program in sorted(accuracy.PROGRAMS):
                            job = jobs[name, mode, program]
                            check_command(job, plan["config"])
                            audit_job(reader, job)
                            parsed[program] = accuracy.read_job(reader, job, record)
                        if parsed["LinearPKRaccess"]["widths"] != [1, 4, 8, 16]:
                            raise ValueError("unexpected accessibility window lengths")
                        check = accuracy.check_consistency(name, mode, 0, parsed)
                        checks[check["status"]] += 1
                        for k in CHECKS:
                            check_maxima[k] = max(check_maxima[k], check[k])
                        if check["status"] != "success":
                            raise ValueError("CLI consistency check failed: " + str(check))
                        if reference is None:
                            reference = parsed
                        beam = BEAMS[MODES.index(mode)]
                        rows[mode] = {"input": name, "run": run.name, "tier": task["tier"], "task": task["task_id"],
                            "mode": mode, "beam": beam, "local_cap": max(beam, 30), "length": len(record[1]),
                            "input_sha256": hashes[name], **desc, **compare_condition(parsed, reference, record),
                            **{k: check[k] for k in CHECKS}, **resource_fields(jobs, name, mode)}
                    for mode in MODES:
                        writer.writerow({k: "NA" if v is None else v for k, v in rows[mode].items()})
                    del parsed, reference, rows
                if number % 10 == 0 or number == len(tasks):
                    print(f"{run.name} {number}/{len(tasks)} tasks; {len(seen)} RNA; {time.monotonic()-start:.1f}s", flush=True)
                del reader, jobs, records
                gc.collect()
    provenance.update(inputs=len(seen), jobs=dict(statuses), consistency_checks=dict(checks),
                      consistency_maxima=check_maxima, task_hosts=dict(hosts), task_workers=dict(workers),
                      metrics_sha256=digest(output / "metrics.tsv"), elapsed_seconds=time.monotonic()-start,
                      validation="all selected jobs, input/feature hashes, used artifact hashes, raw commands/results/resources, model and dependency checks passed")
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")
    return provenance


def statistics(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return {"n": 0, **dict.fromkeys(("min", "p10", "p25", "median", "p75", "p90", "p95", "p99", "max", "mean"))}
    q = np.quantile(values, [.1, .25, .5, .75, .9, .95, .99], method="linear")
    return {"n": len(values), "min": float(values.min()), "max": float(values.max()), "mean": float(values.mean()),
            **{k: float(v) for k, v in zip(("p10", "p25", "median", "p75", "p90", "p95", "p99"), q)}}


def threshold_counts(values):
    values = np.asarray(values, dtype=float)
    valid = values[np.isfinite(values)]
    return {str(t): {"eligible_sequences": len(valid),
                    **{label: {"count": int(test.sum()), "fraction": float(test.mean()) if len(test) else None}
                       for label, test in (("positive", valid > t), ("negative", valid < -t), ("absolute", abs(valid) > t))}}
            for t in (.01, .1, 1)}


def load_metrics(output):
    provenance = json.loads((output / "provenance.json").read_text())
    if digest(output / "metrics.tsv") != provenance["metrics_sha256"]:
        raise ValueError("metrics checksum mismatch")
    blocks = defaultdict(list)
    with (output / "metrics.tsv").open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            blocks[row["input"]].append(row)
    if len(blocks) != provenance["inputs"]:
        raise ValueError("metrics RNA count mismatch")
    names, metadata = [], []
    numeric = [k for k in FIELDS if k not in IDENTITY or k in ("length", "beam", "local_cap")]
    data = {k: [] for k in numeric}
    for name, rows in blocks.items():
        if tuple(r["mode"] for r in rows) != MODES:
            raise ValueError("incomplete or unordered four-condition RNA block")
        if any(r[k] != rows[0][k] for r in rows for k in FEATURES + ("length", "sequence_sha256", "scaffold_sha256")):
            raise ValueError("descriptor differs between beam conditions")
        names.append(name)
        metadata.append({k: rows[0][k] for k in IDENTITY if k not in ("beam", "local_cap", "mode")})
        for k in numeric:
            data[k].append([np.nan if r[k] == "NA" else float(r[k]) for r in rows])
    return provenance, names, metadata, {k: np.asarray(v) for k, v in data.items()}


def association(x, y, length):
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(length)
    x, y, length = x[valid], y[valid], length[valid]
    result = {"n": len(x), "spearman": None, "partial_spearman_length": None}
    if len(x) < 6:
        return result
    rx, ry, rn = [np.asarray(ranks(v.tolist())) for v in (x, y, length)]
    if np.std(rx) == 0 or np.std(ry) == 0:
        return result
    result["spearman"] = float(np.corrcoef(rx, ry)[0, 1])
    controls = np.column_stack((np.ones(len(rn)), rn))
    ex = rx - controls @ np.linalg.lstsq(controls, rx, rcond=None)[0]
    ey = ry - controls @ np.linalg.lstsq(controls, ry, rcond=None)[0]
    if np.linalg.norm(ex) > 1e-10 and np.linalg.norm(ey) > 1e-10:
        result["partial_spearman_length"] = float(np.corrcoef(ex, ey)[0, 1])
    return result


def candidate_rows(names, metadata, data):
    rows, chosen_hashes = [], set()
    delta, tv = data[CORE[0]][:, 2], data["base_state_tv_mean_Gfree"][:, 2]
    length = data["length"][:, 0]
    rise = np.maximum(0, np.diff(data["free_energy_kcal_mol"], axis=1)).max(axis=1)
    for tier, lo, hi in TIERS:
        ids = np.flatnonzero((length >= lo) & (length <= hi)).tolist()
        if not ids:
            continue
        remaining = set(ids)
        def select(order, reason):
            for i in order:
                hashes = tuple(metadata[i][k] for k in ("sequence_sha256", "scaffold_sha256"))
                if i not in remaining or hashes in chosen_hashes:
                    continue
                remaining.remove(i)
                chosen_hashes.add(hashes)
                row = {"tier": tier, "input": names[i], "reason": reason,
                       "proposed_beams": "300,400", "proposed_local_caps": "300,400",
                       **{k: metadata[i][k] for k in ("run", "task", "input_sha256", "sequence_sha256", "scaffold_sha256")},
                       "length": int(length[i]), **{k: data[k][i, 0] for k in FEATURES},
                       "F100_minus_F200_kcal_mol": delta[i], "Z100_over_Z200": data["Z_ratio_vs_b200"][i, 2],
                       "BPP_TV100_vs200": tv[i], "context_TV100_vs200": data["context_tv_mean_Gfree"][i, 2],
                       "largest_F_rise_kcal_mol": rise[i],
                       "b200_elapsed_seconds": data["LinearCParty_elapsed_seconds"][i, 3],
                       "b200_peak_rss_mib": data["LinearCParty_peak_rss_kib"][i, 3] / 1024}
                rows.append(row)
                return i
            return None
        worst = select(sorted(ids, key=lambda i: (-abs(delta[i]), names[i])), "largest_abs_F100_minus_F200")
        select(sorted(ids, key=lambda i: (-tv[i], names[i])), "largest_BPP_TV100_vs200_remaining")
        nonmonotone = [i for i in remaining if rise[i] > 1e-8 and
                      tuple(metadata[i][k] for k in ("sequence_sha256", "scaffold_sha256")) not in chosen_hashes]
        if nonmonotone:
            select(sorted(nonmonotone, key=lambda i: (-rise[i], names[i])), "largest_nonmonotone_F_rise_remaining")
        else:
            select(sorted(ids, key=lambda i: (-abs(delta[i]), names[i])), "second_large_abs_F_remaining")
        if worst is not None:
            # A control near the difficult input's length, from the lower third
            # of BPP sensitivity within this tier; do not imply exact accuracy.
            limit = float(np.quantile(tv[ids], 1/3))
            low = [i for i in remaining if tv[i] <= limit]
            select(sorted(low, key=lambda i: (abs(length[i]-length[worst]), abs(delta[i]), names[i])),
                   "low_BPP_sensitivity_length_matched_control")
    return rows


def summarize(output, make_plots=True):
    provenance, names, metadata, data = load_metrics(output)
    length = data["length"][:, 0]
    bins = [("all", 400, 4381)] + list(TIERS)
    bins += [(f"fine_nt{lo}_{hi}", lo, hi) for lo, hi in
             [(lo, lo+99) for lo in range(400,1600,100)] + [(lo,min(lo+199,4381)) for lo in range(1600,4382,200)]]
    summary = {"schema_version": 1, "generated_utc": now(), "inputs": len(names),
               "reference": "b=200, h=d=200, largest measured width; not exact",
               "quantiles": "RNA-weighted; numpy linear interpolation; finite values only, denominator retained",
               "condition_sources": dict(zip(map(str, BEAMS), MODES)), "metrics_sha256": provenance["metrics_sha256"],
               "bins": [], "monotonicity": [], "associations": [], "paired_resource_ratios": [],
               "association_scope": "descriptive Spearman and correlation of rank residuals after regressing both ranks on length rank; no causal or theoretical complexity inference"}
    for label, lo, hi in bins:
        selected = (length >= lo) & (length <= hi)
        group = {"bin": label, "min_length": lo, "max_length": hi, "sequences": int(selected.sum()), "conditions": []}
        group["local_cap_geometry"] = {
            "maximum_case3_branching_ancestors": float(data["case3_ancestors_max"][selected, 0].max()) if selected.any() else None,
            "sequences_with_free_run_gt_h_plus2": {str(h): int((data["free_run_max"][selected, 0] > h+2).sum()) for h in (30,50,100,200)},
            "interpretation": "free-run length permits geometries excluded by h; canonical pairing, retained states and posterior mass not assessed"}
        for bi, beam in enumerate(BEAMS):
            delta = data[CORE[0]][selected, bi]
            group["conditions"].append({"beam": beam,
                "metrics": {k: statistics(data[k][selected, bi]) for k in METRICS},
                "delta_F_thresholds": threshold_counts(delta),
                "negative_delta_F_below_minus_1e8": int((delta < -1e-8).sum()),
                "Z_ratio_underflow_zero_sequences": int((data["Z_ratio_vs_b200"][selected, bi] == 0).sum()),
                "window_counts": {k: {"NA_sequences": int(np.isnan(data[k][selected, bi]).sum()),
                                      "sum": float(np.nansum(data[k][selected, bi]))} for k in ACCESS if k not in METRICS}})
        summary["bins"].append(group)
        for key, direction, tolerance in (("free_energy_kcal_mol", 1, 1e-8), ("log_Z_ratio_vs_b200", -1, 1e-8 / data["rt_kcal_mol"][0, 0]),
                                          ("base_state_tv_mean_Gfree", 1, 1e-9), ("context_tv_mean_Gfree", 1, 1e-9)):
            matrix = data[key][selected, :]
            stat, differences, valid = monotonicity(matrix, direction, tolerance)
            steps = [{"from_beam": BEAMS[j], "to_beam": BEAMS[j+1],
                      "violations": int(((differences[:,j] > tolerance) & valid).sum()),
                      "directional_change": statistics(differences[valid,j])} for j in range(3)]
            strict, _, _ = monotonicity(matrix, direction, 0)
            summary["monotonicity"].append({"bin": label, "metric": key, **stat, "steps": steps,
                                            "zero_tolerance_monotone_sequences": strict["monotone_sequences"]})
        for bi, beam in enumerate(BEAMS[:-1]):
            ratios = {}
            for suffix in ("elapsed_seconds", "time_wall_seconds", "peak_rss_kib"):
                a, b = data["LinearCParty_" + suffix][selected, 3], data["LinearCParty_" + suffix][selected, bi]
                valid = np.isfinite(a) & np.isfinite(b) & (a >= 0) & (b > 0)
                ratios[suffix] = {**statistics(a[valid] / b[valid]), "excluded_sequences": int((~valid).sum())}
            summary["paired_resource_ratios"].append({"bin": label, "numerator_beam": 200, "denominator_beam": beam, "ratios": ratios})
        if label == "all" or label in [t[0] for t in TIERS]:
            outcomes = {"signed_F100_minus_F200": data[CORE[0]][:, 2], "abs_F100_minus_F200": data[CORE[1]][:, 2],
                        "BPP_TV100_vs200": data["base_state_tv_mean_Gfree"][:, 2],
                        "context_TV100_vs200": data["context_tv_mean_Gfree"][:, 2],
                        "b200_elapsed_seconds": data["LinearCParty_elapsed_seconds"][:, 3],
                        "b200_peak_rss_kib": data["LinearCParty_peak_rss_kib"][:, 3]}
            for outcome, values in outcomes.items():
                for feature in ("length",) + FEATURES:
                    summary["associations"].append({"bin": label, "outcome": outcome, "feature": feature,
                        **association(data[feature][selected, 0], values[selected], length[selected])})
            group["feature_distributions_by_abs_F100_minus_F200_gt1"] = {}
            for label2, test in (("above_1", data[CORE[1]][:,2] > 1), ("at_most_1", data[CORE[1]][:,2] <= 1)):
                group["feature_distributions_by_abs_F100_minus_F200_gt1"][label2] = {
                    k: statistics(data[k][selected & test, 0]) for k in ("length",) + FEATURES}
    summary["worst_examples"] = {}
    for key in (CORE[1], "base_state_tv_mean_Gfree", "context_tv_mean_Gfree"):
        ids = np.argsort(-data[key][:, 2], kind="stable")[:15]
        summary["worst_examples"][key] = [{"input": names[i], "length": int(length[i]),
            **{k: float(data[k][i, 2]) for k in (CORE[0], "Z_ratio_vs_b200", "base_state_tv_mean_Gfree", "context_tv_mean_Gfree")},
            **{k: float(data[k][i, 0]) for k in FEATURES}} for i in ids]
    summary["nonmonotone_examples"] = {}
    for key, tol in (("free_energy_kcal_mol", 1e-8), ("base_state_tv_mean_Gfree", 1e-9)):
        rise = np.diff(data[key], axis=1)
        ids = np.flatnonzero((rise > tol).any(axis=1))
        ids = sorted(ids, key=lambda i: -max(rise[i]))
        summary["nonmonotone_examples"][key] = [{"input": names[i], "length": int(length[i]),
            "values_by_beam": data[key][i].tolist(), "largest_rise": float(max(rise[i]))} for i in ids]
    candidates = candidate_rows(names, metadata, data)
    summary["pilot_candidates"] = candidates
    if candidates:
        with (output / "pilot_candidates.tsv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, list(candidates[0]), delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(candidates)
    summary["candidate_sha256"] = digest(output / "pilot_candidates.tsv")
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    if make_plots:
        plot(output, data)
    return summary


def plot(output, data):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colors = ("#d55e00", "#e69f00", "#0072b2", "#009e73")
    length = data["length"][:, 0]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), constrained_layout=True)
    panels = ((CORE[0], "F(b) - F(200) [kcal/mol]"), ("log_Z_ratio_vs_b200", "log[Z(b) / Z(200)]"),
              ("base_state_tv_mean_Gfree", "BPP TV on free bases"), ("context_tv_mean_Gfree", "Context TV on free bases"),
              ("access_w4_probability_mae_eligible", "Accessibility MAE, w=4"),
              ("access_w16_probability_mae_eligible", "Accessibility MAE, w=16"))
    for ax, (key, title) in zip(axes.flat, panels):
        for bi, beam in enumerate(BEAMS[:-1]):
            ax.scatter(length, data[key][:, bi], s=5, alpha=.3, color=colors[bi], linewidths=0, label=f"b={beam}")
        if key in CORE[:4]:
            ax.set_yscale("symlog", linthresh=.01)
            ax.axhline(0, color="grey", linewidth=.6)
        ax.set(title=title, xlabel="Length [nt]")
        ax.grid(alpha=.15)
    axes[0,0].legend(markerscale=3)
    fig.suptitle("bpRNA beam-width sensitivity: reference b=200, h=d=max(b,30); no exact reference", fontsize=12)
    fig.savefig(output / "beam_sensitivity.png", dpi=170)
    plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), constrained_layout=True)
    for column, (key, title, divisor) in enumerate((("LinearCParty_elapsed_seconds", "Elapsed time [s]", 1),
                                                  ("LinearCParty_peak_rss_kib", "Peak RSS [MiB]", 1024))):
        ax = axes[0,column]
        for bi, beam in enumerate(BEAMS):
            ax.scatter(length, data[key][:,bi] / divisor, s=5, alpha=.35, linewidths=0, color=colors[bi], label=f"b={beam}")
        ax.set(xscale="log", yscale="log", xlabel="Length [nt]", title=title)
        ax.grid(alpha=.15)
        boxes = []
        for _, lo, hi in TIERS:
            select = (length >= lo) & (length <= hi) & (data[key][:,2] > 0)
            boxes.append(data[key][select,3] / data[key][select,2])
        axes[1,column].boxplot(boxes, tick_labels=[f"{lo}–{hi}" for _,lo,hi in TIERS], showfliers=True,
                              flierprops={"markersize": 2, "alpha": .3})
        axes[1,column].set(yscale="log", ylabel="Paired ratio b=200 / b=100", xlabel="Length bin [nt]")
        axes[1,column].axhline(1, color="grey", linewidth=.7)
        axes[1,column].grid(axis="y", alpha=.15)
    axes[0,0].legend(markerscale=3)
    fig.suptitle("LinearCParty resources: accuracy arrays, one observation per condition; descriptive only", fontsize=12)
    fig.savefig(output / "time_rss.png", dpi=170)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run400", type=Path, required=True)
    parser.add_argument("--run800", type=Path, required=True)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary-only", action="store_true")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()
    if not args.summary_only:
        write_metrics([(args.run400, (TIERS[0][0],)), (args.run800, tuple(t[0] for t in TIERS[1:]))], args.features, args.output)
    summary = summarize(args.output, not args.no_plots)
    print(json.dumps({"inputs": summary["inputs"], "pilot_candidates": len(summary["pilot_candidates"])}))


if __name__ == "__main__":
    main()
