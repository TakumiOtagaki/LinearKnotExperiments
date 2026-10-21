#!/usr/bin/env python3
"""Analyze saved E2 repetitions, scaffold contrasts, and exploratory log models.

No RNA calculation or scheduler calls. Tar members are read as streams; no
extraction. Fit only cpu_40 synthetic input medians, separately for each CLI.
All three seeds of a geometry stay together in length/family validation folds.
"""

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import tarfile

import numpy as np

from scaffold_features import scaffold_features
from task_storage import RunReader, input_records

PROGRAMS = ("LinearCParty", "LinearPKCapR", "LinearPKRaccess")
OUTCOMES = ("wall_seconds", "rss_gib")
FEATURE_BLOCKS = {
    "length_free": ("n", "free_bases"),
    "free_loop": ("same_loop_pairs_turn", "loop_free_max", "v_geometry_beam_cap_sum"),
    "local_runs": ("hairpin_free_geometries", "free_run_max", "free_runs"),
    "BE_depth": ("B_G", "be_geometry_beam_cap_sum", "pair_depth_max"),
    "CASE3": ("case3_path_entries", "case3_ancestors_max", "case3_sibling_endpoint_pairs"),
    "canonical_crossing": ("cross_loop_canonical_pairs_turn",),
}
FULL_FEATURES = tuple(k for block in FEATURE_BLOCKS.values() for k in block)
MODEL_SPECS = {
    "n_only": (("n",), 0.),
    "free_loop": (("n", "same_loop_pairs_turn", "loop_free_max"), 1.),
    "geometry_core": (("n", "same_loop_pairs_turn", "loop_free_max", "B_G",
                       "case3_path_entries", "case3_ancestors_max"), 1.),
    "geometry_full": (FULL_FEATURES, 1.),
}


def distribution(values):
    values = np.asarray(values, dtype=float)
    if not len(values) or not np.all(np.isfinite(values)):
        raise ValueError("nonempty finite observations required")
    return {"count": len(values), **{k: float(v) for k, v in zip(
        ("min", "p10", "median", "p90", "p95", "max"), np.quantile(values, (0, .1, .5, .9, .95, 1)))}}


def repetition_summary(jobs):
    if len(jobs) != 3 or sorted(j["repetition"] for j in jobs) != [0, 1, 2]:
        raise ValueError("expected exactly measured repetitions 0,1,2")
    ordered = sorted(jobs, key=lambda j: j["repetition"])
    result = {"measured_repetitions": 3}
    for key, prefix, scale in (("time_wall_seconds", "wall_seconds", 1), ("peak_rss_kib", "rss_gib", 2**20)):
        if any(j["status"] != "success" or j.get(key) is None or
               not math.isfinite(j[key]) or j[key] <= 0 or j["warmup"] for j in ordered):
            raise ValueError("unsuccessful, warmup, or nonpositive measured observation")
        values = [j[key] / scale for j in ordered]
        median = statistics.median(values)
        result.update({prefix + "_repetitions": values, prefix + "_median": median,
                       prefix + "_min": min(values), prefix + "_max": max(values),
                       prefix + "_relative_range": (max(values) - min(values)) / median})
    return result


def model_errors(actual, predicted):
    actual, predicted = np.asarray(actual, float), np.asarray(predicted, float)
    residual = np.log(predicted) - np.log(actual)
    total = float(np.sum((np.log(actual) - np.log(actual).mean())**2))
    return {"n": len(actual), "log_rmse": float(np.sqrt(np.mean(residual**2))),
            "multiplicative_rmse": float(np.exp(np.sqrt(np.mean(residual**2)))),
            "log_r2": 1 - float(np.sum(residual**2)) / total if total else None,
            "median_absolute_percent_error": float(np.median(np.abs(predicted / actual - 1)) * 100),
            "predicted_over_observed": distribution(predicted / actual)}


def fit_log_model(rows, outcome, features, alpha=1.):
    x = np.log1p(np.asarray([[r[k] for k in features] for r in rows], float))
    y = np.log([r[outcome + "_median"] for r in rows])
    center, scale = x.mean(0), x.std(0)
    scale[scale == 0] = 1
    z = (x - center) / scale
    if alpha:
        beta = np.linalg.solve(z.T @ z + alpha * np.eye(len(features)), z.T @ (y - y.mean()))
    else:
        beta = np.linalg.lstsq(z, y - y.mean(), rcond=None)[0]
    unscaled = beta / scale
    intercept = float(y.mean() - center @ unscaled)
    model = {"features": list(features), "alpha": alpha, "intercept": intercept,
             "coefficients": {k: float(v) for k, v in zip(features, unscaled)},
             "standardized_coefficients": {k: float(v) for k, v in zip(features, beta)},
             "log1p_centers": dict(zip(features, center.tolist())),
             "log1p_scales": dict(zip(features, scale.tolist())),
             "design_rank": int(np.linalg.matrix_rank(z)),
             "design_condition_number": float(np.linalg.cond(z)),
             "training_range": {k: [min(r[k] for r in rows), max(r[k] for r in rows)] for k in features}}
    model["training"] = model_errors([r[outcome + "_median"] for r in rows], predict(model, rows))
    return model


def predict(model, rows):
    return np.exp([model["intercept"] + sum(beta * math.log1p(r[k])
                  for k, beta in model["coefficients"].items()) for r in rows])


def grouped_validation(rows, outcome, features, alpha, group_key):
    predictions = np.empty(len(rows))
    folds = []
    for group in sorted({r[group_key] for r in rows}):
        train = [r for r in rows if r[group_key] != group]
        indices = [i for i, r in enumerate(rows) if r[group_key] == group]
        test = [rows[i] for i in indices]
        model = fit_log_model(train, outcome, features, alpha)
        predictions[indices] = predict(model, test)
        folds.append({"held_out": group, "train_inputs": len(train),
                      **model_errors([r[outcome + "_median"] for r in test], predictions[indices])})
    return {"pooled": model_errors([r[outcome + "_median"] for r in rows], predictions),
            "folds": folds}, predictions


def log_contributions(model, a, b):
    return {k: beta * (math.log1p(b[k]) - math.log1p(a[k]))
            for k, beta in model["coefficients"].items()}


def file_digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(2**20), b""):
            h.update(chunk)
    return h.hexdigest()


def audit_packed_payload(path, jobs):
    """One sequential scan, hash every raw member, retain only small scalars."""
    expected = {name: entry for j in jobs for name, entry in j["artifacts"].items()}
    actual, scalars = set(), {}
    with tarfile.open(path, "r|*", bufsize=1024 * 1024) as archive:
        for member in archive:
            if not member.isfile():
                continue
            if member.name not in expected or member.name in actual:
                raise ValueError("unplanned or duplicate packed artifact: " + member.name)
            h, chunks = hashlib.sha256(), []
            keep = member.name.endswith(("/time.tsv", "/result.json"))
            if keep and member.size > 65536:
                raise ValueError("unexpected scalar artifact size")
            stream = archive.extractfile(member)
            for chunk in iter(lambda: stream.read(2**20), b""):
                h.update(chunk)
                if keep:
                    chunks.append(chunk)
            entry = expected[member.name]
            if member.size != entry["bytes"] or h.hexdigest() != entry["sha256"]:
                raise ValueError("packed artifact hash/size mismatch: " + member.name)
            actual.add(member.name)
            if keep:
                scalars[member.name] = b"".join(chunks)
            archive.members.clear()
    if actual != set(expected):
        raise ValueError("missing packed artifacts")
    for j in jobs:
        prefix = j["directory"] + "/"
        raw = json.loads(scalars[prefix + "result.json"])
        if any(raw.get(k) != j.get(k) for k in ("id", "input", "program", "mode", "warmup",
                "repetition", "status", "returncode", "time_wall_seconds", "peak_rss_kib", "command")):
            raise ValueError("packed result/checkpoint mismatch")
        time, rss = scalars[prefix + "time.tsv"].decode().strip().split("\t")
        if float(time) != j["time_wall_seconds"] or int(rss) != j["peak_rss_kib"]:
            raise ValueError("GNU time/checkpoint mismatch")
    return len(actual)


def load_dataset(directory):
    directory = Path(directory)
    records = {}
    stored = {r["name"]: r for r in csv.DictReader((directory / "features.tsv").open(), delimiter="\t")}
    manifest = json.loads((directory / "manifest.json").read_text())
    for line in (directory / "inputs.jsonl").open():
        r = json.loads(line)
        if r["name"] in records:
            raise ValueError("duplicate synthetic input")
        f = json.loads(json.dumps(scaffold_features(r["sequence"], r["scaffold"])))
        if f != r["features"] or any(json.loads(stored[r["name"]][k]) != v
                for k, v in f.items() if not k.endswith("sha256")) or any(
                stored[r["name"]][k] != f[k] for k in ("sequence_sha256", "scaffold_sha256")):
            raise ValueError("synthetic saved/recomputed feature mismatch")
        records[r["name"]] = r
    if len(records) != manifest["records"] or set(records) != set(stored):
        raise ValueError("synthetic feature cohort mismatch")
    return records, {p.name: file_digest(p) for p in directory.iterdir() if p.is_file()}


def selected_task_jobs(jobs, selected, task_directory):
    # Unexecuted planner placeholders need not carry execution task metadata.
    return [j for j in jobs if j["input"] in selected and j["directory"].startswith(task_directory + "/run/")]


def check_revision(actual, expected):
    """Allow a pinned full commit or hexadecimal prefix; omission preserves history."""
    if expected is not None and (not re.fullmatch(r"[0-9a-f]{7,40}", expected)
                                 or not actual.startswith(expected)):
        raise ValueError(f"source revision mismatch: expected {expected}, found {actual}")


def load_run(path, label, dataset, caps, resource, selected_lengths=None, expected_workers=1,
             expected_revision=None):
    reader = RunReader(path)
    m, groups = reader.manifest, defaultdict(list)
    c = m["config"]
    if c["warmups"] != 1 or c["repetitions"] != 3 or len(c["modes"]) != 1:
        raise ValueError("expected one mode with warmup1/measured3")
    args = c["modes"][0]["args"]
    actual_caps = tuple(int(args[args.index(flag) + 1]) for flag in ("-b", "--max-hairpin", "--max-case3-depth"))
    if actual_caps != caps or tuple(p["name"] for p in c["programs"]) != PROGRAMS:
        raise ValueError("unexpected CLI/cap configuration")
    records = input_records(reader)
    selected = {name for name, (_, s, _) in records.items() if selected_lengths is None or len(s) in selected_lengths}
    if any(j["input"] not in records for j in m["jobs"]):
        raise ValueError("unplanned job input")
    for j in m["jobs"]:
        if j["input"] in selected:
            groups[j["input"], j["program"]].append(j)
    tasks, task_by_number, signature, artifacts = [], {}, None, 0
    audit_keys = ("id", "input", "program", "mode", "warmup", "repetition", "status", "returncode",
                  "time_wall_seconds", "peak_rss_kib", "command")
    for t in m["tasks"]:
        tj = selected_task_jobs(m["jobs"], selected, t["directory"])
        if not tj:
            continue
        if t["status"] != "success" or t["workers"] != expected_workers or t["requested_resource"] != resource:
            raise ValueError("selected timing task has wrong status/allocation/workers")
        if expected_workers > 1 and not (t.get("parallel_timing", {}).get("allowed") is True
                and t["parallel_timing"].get("reason", "").strip()):
            raise ValueError("parallel timing lacks explicit authorization metadata")
        if len({a["hostname"] for a in t["attempts"]}) != 1:
            raise ValueError("timing repetitions pooled across hosts")
        base = Path(path) / t["directory"] / "run"
        metadata = json.loads((base / "run.json").read_text())
        if file_digest(base / "run.json") != t["run_manifest_sha256"]:
            raise ValueError("task manifest hash changed")
        if metadata["binaries_before"] != metadata["binaries_after"] or any(metadata[k] for k in
                ("software_changed_during_run", "binaries_changed_during_run")) or not metadata["software_lock_valid_after_run"]:
            raise ValueError("task dependencies changed")
        if any(metadata["environment"].get(k) != "1" for k in
               ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")):
            raise ValueError("task not single threaded")
        sig = {"binaries": metadata["binaries_before"], "energy_parameter_sha256": metadata["energy_parameter_sha256"],
               "commit": metadata["software_before"]["commit"], "build_settings": metadata["build_settings"],
               "programs": c["programs"], "model_settings_temperature": c["model_settings"]["temperature_celsius"]}
        check_revision(sig["commit"], expected_revision)
        if signature is not None and sig != signature:
            raise ValueError("inconsistent task software/model")
        signature = sig
        checkpoint = base / "task_results.jsonl"
        if checkpoint.exists():
            durable = [json.loads(line) for line in checkpoint.open()]
            if len({j["id"] for j in durable}) != len(durable):
                raise ValueError("duplicate checkpoint")
            selected_checkpoint = {j["id"]: j for j in durable if j["input"] in selected}
            if set(selected_checkpoint) != {j["id"] for j in tj} or any(
                    any(selected_checkpoint[j["id"]].get(k) != j.get(k) for k in audit_keys) for j in tj):
                raise ValueError("aggregate/checkpoint mismatch")
            artifacts += audit_packed_payload(base / "task_outputs.tar", durable)
        else:
            for j in tj:
                raw = json.loads(reader.artifact(j, "result.json"))
                if any(raw.get(k) != j.get(k) for k in audit_keys):
                    raise ValueError("legacy raw/aggregate mismatch")
                time, rss = reader.artifact(j, "time.tsv").decode().strip().split("\t")
                if float(time) != j["time_wall_seconds"] or int(rss) != j["peak_rss_kib"]:
                    raise ValueError("legacy GNU time mismatch")
                artifacts += 2
        info = {"global_task_number": t["global_task_number"], "directory": t["directory"],
                "hostname": t["attempts"][0]["hostname"], "task_manifest_sha256": file_digest(base / "run.json"),
                "checkpoint_sha256": file_digest(checkpoint) if checkpoint.exists() else None,
                "raw_artifact_hashes_verified": checkpoint.exists()}
        if all(j.get("started_utc") and j.get("finished_utc") for j in tj):
            events = sorted([(j["started_utc"], 1) for j in tj] + [(j["finished_utc"], -1) for j in tj])
            active = peak = 0
            for _, change in events:
                active += change
                peak = max(peak, active)
            info["observed_max_active_processes"] = peak
        tasks.append(info)
        task_by_number[t["global_task_number"]] = info
    rows = []
    for name in sorted(selected):
        header, s, g = records[name]
        if header != name:
            raise ValueError("record name mismatch")
        f = scaffold_features(s, g, *caps)
        saved = dataset.get(name)
        if saved and (saved["sequence"] != s or saved["scaffold"] != g):
            raise ValueError("timing/dataset S/G mismatch")
        family = saved["family"] if saved else "rnacentral"
        series = family + ("_d" + str(saved["branch_budget"]) if family == "branch_depth" else "")
        for program in PROGRAMS:
            jobs = groups[name, program]
            if len(jobs) != 4 or sum(j["warmup"] for j in jobs) != 1 or any(j["status"] != "success" for j in jobs):
                raise ValueError("missing/failed warmup or measurement")
            if len({j["global_task_number"] for j in jobs}) != 1:
                raise ValueError("combination spread across tasks")
            measured = [j for j in jobs if not j["warmup"]]
            task = task_by_number[jobs[0]["global_task_number"]]
            rows.append({"run": label, "name": name, "program": program, "family": family, "series": series,
                         "seed": saved["seed"] if saved else None, "branch_budget": saved["branch_budget"] if saved else None,
                         "resource": resource, "workers": expected_workers, "task": task["directory"], "hostname": task["hostname"],
                         **repetition_summary(measured), **f})
    provenance = {"path": str(path), "manifest_sha256": reader.manifest_sha256, "overall_status": m["status"],
                  "selected_inputs": len(selected), "selected_tasks": len(tasks), "measured_jobs": len(rows) * 3,
                  "warmups": len(rows), "unselected_planned_jobs": sum(j["input"] not in selected for j in m["jobs"]),
                  "unselected_statuses": dict(Counter(j["status"] for j in m["jobs"] if j["input"] not in selected)),
                  "caps": dict(zip(("b", "h", "d"), caps)), "resource": resource, "workers": expected_workers, "signature": signature,
                  "memory_limit_mib": c.get("memory_limit_mib"), "mode": c["modes"][0],
                  "raw_artifacts_verified": artifacts, "tasks": tasks}
    return rows, provenance


def pair_comparison(kind, a, b, models=None):
    if a["program"] != b["program"]:
        raise ValueError("paired comparison must use same CLI")
    row = {"kind": kind, "program": a["program"], "run_a": a["run"], "run_b": b["run"],
           "name_a": a["name"], "name_b": b["name"], "family_a": a["family"], "family_b": b["family"],
           "series_a": a["series"], "series_b": b["series"], "seed": a["seed"], "n_a": a["n"], "n_b": b["n"],
           "n_ratio": b["n"] / a["n"], "same_S_G": all(a[k] == b[k] for k in ("sequence_sha256", "scaffold_sha256"))}
    for outcome in OUTCOMES:
        row[outcome + "_a"] = a[outcome + "_median"]
        row[outcome + "_b"] = b[outcome + "_median"]
        row[outcome + "_ratio"] = b[outcome + "_median"] / a[outcome + "_median"]
        row[outcome + "_difference"] = b[outcome + "_median"] - a[outcome + "_median"]
        if models:
            model = models[a["program"]][outcome]["geometry_full"]
            contributions = log_contributions(model, a, b)
            row[outcome + "_model_log_contributions"] = contributions
            row[outcome + "_model_ratio"] = math.exp(sum(contributions.values()))
            row[outcome + "_growth_exponent"] = (math.log(row[outcome + "_ratio"]) / math.log(row["n_ratio"])
                                                    if row["n_ratio"] != 1 else None)
    row["feature_values_a"] = {k: a[k] for k in FULL_FEATURES}
    row["feature_values_b"] = {k: b[k] for k in FULL_FEATURES}
    return row


def build_comparisons(cpu, node, rna, historical, models):
    pairs = []
    lookup = {(r["name"], r["program"]): r for r in cpu}
    for r in node:
        a = lookup[r["name"], r["program"]]
        if any(a[k] != r[k] for k in ("sequence_sha256", "scaffold_sha256", "beam", "hairpin_cap", "case3_depth_cap")):
            raise ValueError("unmatched resource comparison")
        pairs.append(pair_comparison("node_f_over_cpu40", a, r))
    old = {(r["name"], r["program"]): r for r in historical}
    for r in rna:
        a = old[r["name"], r["program"]]
        if any(a[k] != r[k] for k in ("sequence_sha256", "scaffold_sha256", "beam")):
            raise ValueError("unmatched cap comparison")
        pairs.append(pair_comparison("rna_cap100_over_cap30", a, r))
    indexed = {(r["series"], r["n"], r["seed"], r["program"]): r for r in cpu}
    for r in cpu:
        b = indexed.get((r["series"], r["n"] * 2, r["seed"], r["program"]))
        if b:
            pairs.append(pair_comparison("n_doubling", r, b, models))
    contrasts = (("separated_helices", "separated_concentrated"), ("balanced_branching", "deep_branching"),
                 ("separated_helices", "separated_12"), ("separated_helices", "long_stem"),
                 ("branch_depth_d30", "branch_depth_d31"), ("branch_depth_d99", "branch_depth_d100"),
                 ("branch_depth_d100", "deep_branching"))
    for sa, sb in contrasts:
        for a in cpu:
            if a["series"] != sa:
                continue
            b = indexed.get((sb, a["n"], a["seed"], a["program"]))
            if b and (sa != "branch_depth_d100" or b["case3_ancestors_max"] == 101):
                pairs.append(pair_comparison("shape_contrast", a, b, models))
    real = {(r["name"], r["program"]): r for r in rna}
    for program in PROGRAMS:
        pairs.append(pair_comparison("rna_13652_over_11832", real["URS00023C9544", program],
                                     real["URS000241043C", program], models))
    return pairs


def write_union_tsv(path, rows):
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "NA" if v is None else json.dumps(v, separators=(",", ":"))
                             if isinstance(v, (list, dict)) else v for k, v in row.items()})


def figures(output, cpu, rna, models):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    series = sorted({r["series"] for r in cpu})
    colors = dict(zip(series, plt.get_cmap("tab20").colors))
    for outcome, ylabel in (("wall_seconds", "CLI wall time (s)"), ("rss_gib", "Peak RSS (GiB)")):
        fig, axes = plt.subplots(1, 3, figsize=(15, 5.7), constrained_layout=True)
        for ax, program in zip(axes, PROGRAMS):
            for shape in series:
                selected = [r for r in cpu if r["program"] == program and r["series"] == shape]
                ns = sorted({r["n"] for r in selected})
                ys, lo, hi = [], [], []
                for n in ns:
                    group = [r for r in selected if r["n"] == n]
                    ys.append(statistics.median(r[outcome + "_median"] for r in group))
                    lo.append(min(r[outcome + "_min"] for r in group))
                    hi.append(max(r[outcome + "_max"] for r in group))
                ax.plot(ns, ys, "o-", color=colors[shape], lw=1, ms=3, label=shape)
                ax.fill_between(ns, lo, hi, color=colors[shape], alpha=.09)
            ax.set(xscale="log", yscale="log", xlabel="Sequence length n (nt)", ylabel=ylabel, title=program)
            ax.set_xticks([250, 500, 1000, 2000, 4000], labels=["250", "500", "1000", "2000", "4000"])
            ax.grid(alpha=.2)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="outside lower center", ncol=4, fontsize=8)
        fig.suptitle("cpu_40; b=h=d=100. Lines: median across 3 seeds; shading: all 9 measured values.", fontsize=11)
        fig.savefig(output / ("scaling_" + ("time" if outcome == "wall_seconds" else "rss") + ".png"), dpi=170)
        plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5), constrained_layout=True)
    for col, program in enumerate(PROGRAMS):
        group = [r for r in cpu if r["program"] == program]
        model = models[program]["wall_seconds"]["geometry_full"]
        for shape in series:
            rows = [r for r in group if r["series"] == shape]
            axes[0, col].scatter([r["same_loop_pairs_turn"] for r in rows],
                                 [r["wall_seconds_median"] for r in rows], s=12, color=colors[shape])
            axes[1, col].scatter(predict(model, rows), [r["wall_seconds_median"] for r in rows],
                                 s=12, color=colors[shape])
        real = [r for r in rna if r["program"] == program]
        axes[1, col].scatter([r["geometry_full_wall_seconds_prediction_node_f"] for r in real],
                             [r["wall_seconds_median"] for r in real], s=45, marker="x", c="black", label="RNAcentral, external")
        target = next(r for r in real if r["name"] == "URS000241043C")
        axes[1, col].annotate("13.7 kb", (target["geometry_full_wall_seconds_prediction_node_f"], target["wall_seconds_median"]),
                              xytext=(5, -13), textcoords="offset points", fontsize=8)
        limit = max(max(r["wall_seconds_median"] for r in real), max(predict(model, real))) * 1.2
        axes[1, col].plot([1, limit], [1, limit], "k--", lw=.7)
        for ax in axes[:, col]:
            ax.set(xscale="log", yscale="log", ylabel="Observed CLI wall time (s)", title=program)
            ax.grid(alpha=.2)
        axes[0, col].set_xlabel("Same-loop free endpoint pairs J_V (TURN=3)")
        axes[1, col].set_xlabel("Full G model prediction (s); RNAcentral time calibrated to node_f")
        axes[1, col].legend(fontsize=8)
    fig.suptitle("Scaffold supply and time. Colored points: synthetic training inputs; black crosses: external predictions.", fontsize=11)
    fig.savefig(output / "feature_time.png", dpi=170)
    plt.close(fig)


def analyze(dataset_path, cpu_path, node_path, rna_path, historical_path, output):
    output = Path(output)
    if output.exists():
        raise ValueError("use a new analysis output directory")
    dataset, dataset_hashes = load_dataset(dataset_path)
    specifications = ((cpu_path, "synthetic_cpu40", (100, 100, 100), "cpu_40=1", None),
                      (node_path, "synthetic_node_f", (100, 100, 100), "node_f=1", {250, 500}),
                      (rna_path, "rnacentral_cap100", (100, 100, 100), "node_f=1", None),
                      (historical_path, "rnacentral_cap30", (100, 30, 30), "node_f=1", None))
    cohorts, provenance = [], {}
    for path, label, caps, resource, lengths in specifications:
        rows, p = load_run(path, label, dataset, caps, resource, lengths)
        cohorts.append(rows)
        provenance[label] = p
    cpu, node, rna, historical = cohorts
    if len(cpu) != 468 or len(node) != 162 or len(rna) != 21 or len(historical) != 21:
        raise ValueError("unexpected E2 cohort counts")
    if any(p["signature"] != provenance["synthetic_cpu40"]["signature"] for p in provenance.values()):
        raise ValueError("cross-run software/model/output mismatch")
    models, ablation = {}, []
    for program in PROGRAMS:
        group = [r for r in cpu if r["program"] == program]
        models[program] = {}
        for outcome in OUTCOMES:
            models[program][outcome] = {}
            for name, (features, alpha) in MODEL_SPECS.items():
                model = fit_log_model(group, outcome, features, alpha)
                model["validation"] = {}
                for group_key in ("n", "family", "seed"):
                    cv, predictions = grouped_validation(group, outcome, features, alpha, group_key)
                    model["validation"]["leave_one_" + group_key + "_out"] = cv
                    if name == "geometry_full":
                        for row, value in zip(group, predictions):
                            row[name + "_" + outcome + "_cv_" + group_key] = float(value)
                models[program][outcome][name] = model
                for row, value in zip(group, predict(model, group)):
                    row[name + "_" + outcome + "_prediction"] = float(value)
            for block, excluded in FEATURE_BLOCKS.items():
                features = [k for k in FULL_FEATURES if k not in excluded]
                entry = {"program": program, "outcome": outcome, "removed_block": block, "removed_features": list(excluded)}
                for key in ("n", "family"):
                    cv, _ = grouped_validation(group, outcome, features, 1., key)
                    baseline = models[program][outcome]["geometry_full"]["validation"]["leave_one_" + key + "_out"]["pooled"]["log_rmse"]
                    entry["leave_one_" + key + "_out"] = {**cv["pooled"], "log_rmse_increase_vs_full": cv["pooled"]["log_rmse"] - baseline}
                ablation.append(entry)
    comparisons = build_comparisons(cpu, node, rna, historical, models)
    calibration, variability, pair_summaries = {}, [], []
    for label, rows in zip(provenance, cohorts):
        for program in PROGRAMS:
            for outcome in OUTCOMES:
                group = [r for r in rows if r["program"] == program]
                variability.append({"run": label, "program": program, "outcome": outcome,
                                    "relative_range": distribution([r[outcome + "_relative_range"] for r in group]),
                                    "over_5_percent": sum(r[outcome + "_relative_range"] > .05 for r in group),
                                    "over_10_percent": sum(r[outcome + "_relative_range"] > .1 for r in group),
                                    "largest_range": sorted(({"name": r["name"], "n": r["n"], "hostname": r["hostname"],
                                                              "values": r[outcome + "_repetitions"],
                                                              "relative_range": r[outcome + "_relative_range"]} for r in group),
                                                             key=lambda x: x["relative_range"], reverse=True)[:3]})
                if label == "synthetic_cpu40":
                    subset = [r for r in group if r["n"] in (250, 500)]
                    variability.append({"run": label + "_matched_n250_500", "program": program, "outcome": outcome,
                                        "relative_range": distribution([r[outcome + "_relative_range"] for r in subset])})
    for kind in ("node_f_over_cpu40", "rna_cap100_over_cap30"):
        for program in PROGRAMS:
            for n in ([None, 250, 500] if kind == "node_f_over_cpu40" else [None]):
                group = [r for r in comparisons if r["kind"] == kind and r["program"] == program and (n is None or r["n_a"] == n)]
                entry = {"kind": kind, "program": program, "n": n,
                         **{outcome + "_ratio": distribution([r[outcome + "_ratio"] for r in group]) for outcome in OUTCOMES}}
                pair_summaries.append(entry)
                if kind == "node_f_over_cpu40" and n is None:
                    calibration[program] = entry["wall_seconds_ratio"]["median"]
    for program in PROGRAMS:
        group = [r for r in rna if r["program"] == program]
        for outcome in OUTCOMES:
            for name, model in models[program][outcome].items():
                raw = predict(model, group)
                calibrated = raw * (calibration[program] if outcome == "wall_seconds" else 1)
                model["external_rnacentral"] = {"raw": model_errors([r[outcome + "_median"] for r in group], raw),
                    "node_f_calibrated": model_errors([r[outcome + "_median"] for r in group], calibrated),
                    "time_calibration_factor": calibration[program] if outcome == "wall_seconds" else 1,
                    "calibration_caveat": "Median node_f/cpu_40 ratio at n250/500; constancy for longer RNAs is untested."}
                for row, value, corrected in zip(group, raw, calibrated):
                    row[name + "_" + outcome + "_prediction"] = float(value)
                    row[name + "_" + outcome + "_prediction_node_f"] = float(corrected)
                    row[name + "_outside_training_features"] = [k for k, (lo, hi) in model["training_range"].items()
                                                                 if not lo <= row[k] <= hi]
    shape_summary = []
    for program in PROGRAMS:
        for series in sorted({r["series"] for r in cpu}):
            group = [r for r in cpu if r["program"] == program and r["series"] == series]
            for n in sorted({r["n"] for r in group}):
                subset = [r for r in group if r["n"] == n]
                entry = {"program": program, "series": series, "n": n, "inputs": len(subset),
                         "features": {k: distribution([r[k] for r in subset]) for k in FULL_FEATURES},
                         **{o: distribution([r[o + "_median"] for r in subset]) for o in OUTCOMES}}
                for outcome in OUTCOMES:
                    model = models[program][outcome]["geometry_full"]
                    entry[outcome + "_mean_centered_log_contributions"] = {
                        k: float(np.mean([beta * (math.log1p(r[k]) - model["log1p_centers"][k]) for r in subset]))
                        for k, beta in model["coefficients"].items()}
                shape_summary.append(entry)
    contrast_summary = []
    keys = sorted({(p["kind"], p["program"], p["series_a"], p["series_b"], p["n_a"], p["n_b"])
                   for p in comparisons if p["kind"] in ("shape_contrast", "n_doubling", "rna_13652_over_11832")})
    for kind, program, sa, sb, na, nb in keys:
        group = [p for p in comparisons if (p["kind"], p["program"], p["series_a"], p["series_b"], p["n_a"], p["n_b"]) ==
                 (kind, program, sa, sb, na, nb)]
        contrast_summary.append({"kind": kind, "program": program, "series_a": sa, "series_b": sb, "n_a": na, "n_b": nb,
                                 **{o + "_ratio": distribution([r[o + "_ratio"] for r in group]) for o in OUTCOMES},
                                 **{o + "_model_ratio": distribution([r[o + "_model_ratio"] for r in group]) for o in OUTCOMES}})
    shape_features = [r for r in cpu if r["program"] == PROGRAMS[0]]
    x = np.log1p([[r[k] for k in FULL_FEATURES] for r in shape_features])
    correlation = np.corrcoef(x.T)
    correlated = [{"a": a, "b": b, "log1p_pearson": float(correlation[i, j])} for i, a in enumerate(FULL_FEATURES)
                  for j, b in enumerate(FULL_FEATURES) if i < j and abs(correlation[i, j]) >= .9]
    summary = {"schema_version": 1, "scope": "Exploratory whole-CLI empirical description; no complexity proof or causal effect.",
               "definitions": {"relative_range": "(max-min)/median of 3 measurements; not a confidence interval",
                   "fit": "ln(y)=intercept+sum beta_k ln(1+x_k). Unpenalized intercept; standardized predictors.",
                   "ridge_objective": "sum squared log residuals + alpha*sum standardized_beta^2; alpha=1 fixed, not tuned",
                   "units": "GNU time wall seconds and peak RSS KiB/2^20 (GiB)",
                   "validation": "Refit centers/scales/coefficients within each training fold. Leave n/family/seed out; input medians only.",
                   "prediction": "exp(fitted log y), no bias correction; conditional geometric prediction, not arithmetic mean",
                   "contrasts": "Same seed/n where possible; shape embedding changes S. Doubling changes S and multiple G descriptors."},
               "input_counts": {label: p["selected_inputs"] for label, p in provenance.items()},
               "variability": variability, "resource_and_cap_comparisons": pair_summaries,
               "node_f_time_calibration": calibration, "models": models, "feature_block_ablation": ablation,
               "strongly_correlated_predictors": correlated, "shape_summary": shape_summary,
               "contrast_summary": contrast_summary,
               "limitations": ["Only 3 repetitions: ranges are descriptive, not precision guarantees.",
                   "Resource comparison covers 54 identical inputs at n250/500, different hosts/dates; no randomized allocation experiment.",
                   "No n-only linear scaling claim: input shape, free run lengths, retained work and startup change with n.",
                   "Correlated G descriptors and changed S prevent identifying individual causal effects; ridge coefficients may be negative.",
                   "RNAcentral n4793-29106 lies beyond training n250-4000 and has different G shapes/compositions; predictions are extrapolations.",
                   "Peak RSS is whole process memory. G supplies no allocator byte identity; no work/live/forcing counters collected.",
                   "Full forward and whole CLI asymptotic upper bounds remain unresolved in the implementation note."],
               "provenance": {"runs": provenance, "synthetic_dataset_sha256": dataset_hashes,
                   "scripts_sha256": {p.name: file_digest(p) for p in (Path(__file__), Path(__file__).with_name("scaffold_features.py"),
                                                                      Path(__file__).with_name("task_storage.py"))},
                   "tar_extractions_to_filesystem": 0, "qsub_calls": 0}}
    output.mkdir(parents=True)
    write_union_tsv(output / "metrics.tsv", sum(cohorts, []))
    write_union_tsv(output / "comparisons.tsv", comparisons)
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    figures(output, cpu, rna, models)
    print(json.dumps({"metrics_rows": sum(map(len, cohorts)), "comparison_rows": len(comparisons),
                      "successful_measured_jobs": sum(len(rows) * 3 for rows in cohorts),
                      "raw_artifacts_verified": sum(p["raw_artifacts_verified"] for p in provenance.values()),
                      "output": str(output)}))
    return summary


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", default=str(root / "data/e2-synthetic-shapes-20261012"))
    parser.add_argument("--cpu-run", default=str(root / "results/runs/e2-synthetic-shapes-cpu40-20261013"))
    parser.add_argument("--node-run", default=str(root / "results/runs/e2-synthetic-shapes-capb-timing-20261012"))
    parser.add_argument("--rna-run", default=str(root / "results/runs/e2-rnacentral-capb-timing-20261012"))
    parser.add_argument("--historical-run", default=str(root / "results/runs/e2-rnacentral-15aa7742-timing-20261011"))
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    analyze(args.dataset, args.cpu_run, args.node_run, args.rna_run, args.historical_run, args.out)


if __name__ == "__main__":
    main()
