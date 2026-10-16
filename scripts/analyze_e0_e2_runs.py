#!/usr/bin/env python3
"""Summarize completed E0/E2 runs using the existing parity and E2 estimators.

Read packed artifacts without extraction. Keep historical analyzers unchanged.
E0 input is the JSON from analyze_e0.py; E2 comparisons allow different versions
and allocations and are descriptive, without a speed calibration.
"""

import argparse
from collections import Counter
import json
import math
from pathlib import Path
import re

import numpy as np

from analyze_e0 import within_tolerance
from analyze_e2_shape import FEATURES, spearman
from analyze_e2_synthetic import (FEATURE_BLOCKS, FULL_FEATURES, MODEL_SPECS,
    OUTCOMES, PROGRAMS, distribution, file_digest, fit_log_model,
    grouped_validation, load_dataset, load_run, model_errors, pair_comparison,
    predict, write_union_tsv)
from task_storage import RunReader


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def summarize_e0(analysis_path, run_path):
    analysis_path = Path(analysis_path)
    result = json.loads(analysis_path.read_text())
    reader = RunReader(run_path)
    if result["run_manifest_sha256"] != reader.manifest_sha256:
        raise ValueError("E0 analysis/run hash mismatch")
    names = {r["name"] for r in reader.manifest["inputs"]}
    rows, summaries = [], []
    for r in result["observations"]:
        if r["input"] not in names or r["status"] != "success":
            raise ValueError("E0 requires complete successful comparisons")
        cohort = "bpRNA" if r["input"].startswith("bpRNA_") else "witness"
        for key, kind in (("metrics", r["comparison_kind"]),
                          ("default_minus_compatible", "weak_closure_vs_compatible")):
            if key in r:
                rows.append({"input": r["input"], "cohort": cohort, "comparison": kind,
                    **r[key], "exceeds_either_tolerance": not within_tolerance(r[key], result["tolerances"]),
                    "closure_attribution_valid": r.get("closure_attribution_valid")})
    for cohort in ("all", "bpRNA", "witness"):
        for kind in sorted({r["comparison"] for r in rows}):
            group = [r for r in rows if r["comparison"] == kind and (cohort == "all" or r["cohort"] == cohort)]
            for quantity in ("delta_log_z", "delta_free_energy"):
                for absolute in (False, True):
                    values = [abs(r[quantity]) if absolute else r[quantity] for r in group]
                    summaries.append({"cohort": cohort, "comparison": kind, "quantity": quantity,
                        "absolute": absolute, **distribution(values),
                        "exceeds_either_tolerance": sum(r["exceeds_either_tolerance"] for r in group)})
    out = analysis_path.parent
    write_union_tsv(out / "comparisons.tsv", rows)
    write_union_tsv(out / "distributions.tsv", summaries)
    save_json(out / "summary.json", {"input_counts": dict(Counter(r["cohort"] for r in rows
        if r["comparison"] == "numerical_parity")), "tolerances": result["tolerances"],
        "numerical_parity_counts": result["numerical_parity_counts"], "closure_counts": result["closure_counts"],
        "distributions": summaries, "provenance": {"analysis_sha256": file_digest(analysis_path),
        "run_manifest_sha256": reader.manifest_sha256, "script_sha256": file_digest(__file__)}})


def audit_outcomes(path):
    """Check exit status and root evidence, independently of timing eligibility."""
    reader = RunReader(path)
    groups = {}
    for warmup in (False, True):
        for program in PROGRAMS:
            jobs = [j for j in reader.manifest["jobs"] if j["warmup"] == warmup and j["program"] == program]
            counts = Counter(planned=len(jobs))
            for job in jobs:
                counts["status_" + job["status"]] += 1
                stdout = reader.artifact(job, "stdout.txt").decode()
                stderr = reader.artifact(job, "stderr.txt").decode()
                unreachable = bool(re.search(r"unreachable|retained no structure|no reachable", stderr, re.I))
                if job["status"] != "success" or job["returncode"] != 0:
                    counts["failed"] += 1
                elif program == "LinearPKCapR":
                    # This CLI fails when profiles::beam_profiles finds a
                    # nonfinite root; it prints no root scalar on success.
                    labels = [line.split()[0] for line in stdout.splitlines() if line and not line.startswith(">")]
                    if labels != ["Bulge", "Exterior", "Hairpin", "Internal", "Multiloop", "PseudoLoop", "Stem"]:
                        raise ValueError("incomplete CapR profile")
                    counts["reachable_via_successful_profile_cli"] += 1
                else:
                    pattern = (r"^ensemble free energy:\s*(\S+)" if program == "LinearCParty" else
                               r"^# ensemble_free_energy=(\S+)")
                    matches = re.findall(pattern, stdout, re.M)
                    if len(matches) != 1:
                        raise ValueError("missing or duplicate ensemble free energy")
                    if math.isfinite(float(matches[0])):
                        counts["finite_ensemble_free_energy"] += 1
                    else:
                        unreachable = True
                counts["unreachable"] += unreachable
            counts.setdefault("failed", 0)
            groups[("warmup" if warmup else "measured") + "/" + program] = dict(counts)
    return groups


def fit_models(cpu):
    models, ablation = {}, []
    for program in PROGRAMS:
        rows = [r for r in cpu if r["program"] == program]
        models[program] = {}
        for outcome in OUTCOMES:
            fitted = models[program][outcome] = {}
            for name, (features, alpha) in MODEL_SPECS.items():
                model = fitted[name] = fit_log_model(rows, outcome, features, alpha)
                model["validation"] = {}
                for key in ("n", "family", "seed"):
                    cv, predictions = grouped_validation(rows, outcome, features, alpha, key)
                    model["validation"]["leave_one_" + key + "_out"] = cv
                    if name == "geometry_full":
                        for row, value in zip(rows, predictions):
                            row[name + "_" + outcome + "_cv_" + key] = float(value)
                for row, value in zip(rows, predict(model, rows)):
                    row[name + "_" + outcome + "_prediction"] = float(value)
            for block, excluded in FEATURE_BLOCKS.items():
                features = [k for k in FULL_FEATURES if k not in excluded]
                entry = {"program": program, "outcome": outcome, "removed_block": block}
                for key in ("n", "family"):
                    cv, _ = grouped_validation(rows, outcome, features, 1., key)
                    baseline = fitted["geometry_full"]["validation"]["leave_one_" + key + "_out"]["pooled"]["log_rmse"]
                    entry["leave_one_" + key + "_out"] = {**cv["pooled"],
                        "log_rmse_increase_vs_full": cv["pooled"]["log_rmse"] - baseline}
                ablation.append(entry)
    return models, ablation


def geometry_comparisons(cpu, rna, models):
    pairs = []
    indexed = {(r["series"], r["n"], r["seed"], r["program"]): r for r in cpu}
    for a in cpu:
        b = indexed.get((a["series"], a["n"] * 2, a["seed"], a["program"]))
        if b:
            pairs.append(pair_comparison("n_doubling", a, b, models))
    contrasts = (("separated_helices", "separated_concentrated"), ("balanced_branching", "deep_branching"),
        ("separated_helices", "separated_12"), ("separated_helices", "long_stem"),
        ("branch_depth_d30", "branch_depth_d31"), ("branch_depth_d99", "branch_depth_d100"),
        ("branch_depth_d100", "deep_branching"))
    for sa, sb in contrasts:
        for a in cpu:
            b = indexed.get((sb, a["n"], a["seed"], a["program"]))
            if a["series"] == sa and b and (sa != "branch_depth_d100" or b["case3_ancestors_max"] == 101):
                pairs.append(pair_comparison("shape_contrast", a, b, models))
    for program in PROGRAMS:
        group = {r["name"]: r for r in rna if r["program"] == program}
        pairs.append(pair_comparison("rna_13652_over_11832", group["URS00023C9544"], group["URS000241043C"], models))
    return pairs


def historical_comparisons(previous, current):
    lookup = {(r["name"], r["program"]): r for r in previous}
    if set(lookup) != {(r["name"], r["program"]) for r in current}:
        raise ValueError("historical cohort mismatch")
    pairs = []
    for b in current:
        a = lookup[b["name"], b["program"]]
        if any(a[k] != b[k] for k in ("sequence_sha256", "scaffold_sha256", "beam")):
            raise ValueError("historical S/G/beam mismatch")
        pairs.append(pair_comparison("current_over_" + a["run"], a, b))
    return pairs


def empirical_slopes(cpu):
    rows = []
    for program in PROGRAMS:
        for series in sorted({r["series"] for r in cpu}):
            for seed in sorted({r["seed"] for r in cpu}):
                group = sorted((r for r in cpu if (r["program"], r["series"], r["seed"]) ==
                                (program, series, seed)), key=lambda r: r["n"])
                if len(group) < 2:
                    continue
                for outcome in OUTCOMES:
                    x, y = np.log([r["n"] for r in group]), np.log([r[outcome + "_median"] for r in group])
                    beta, intercept = np.polyfit(x, y, 1)
                    rows.append({"program": program, "series": series, "seed": seed, "outcome": outcome,
                        "lengths": [r["n"] for r in group], "slope_ln_y_vs_ln_n": float(beta),
                        "intercept": float(intercept), "log_r2": 1 - float(np.sum((y - intercept - beta*x)**2) /
                        np.sum((y - y.mean())**2))})
    return rows


def figures(output, cpu, rna, revision="d21618c5", workers=1):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    series = sorted({r["series"] for r in cpu})
    colors = dict(zip(series, plt.get_cmap("tab20").colors))
    for outcome, ylabel in (("wall_seconds", "CLI wall time (s)"), ("rss_gib", "Peak RSS (GiB)")):
        fig, axes = plt.subplots(1, 3, figsize=(15, 5.7), constrained_layout=True)
        for ax, program in zip(axes, PROGRAMS):
            for shape in series:
                group = [r for r in cpu if r["program"] == program and r["series"] == shape]
                ns = sorted({r["n"] for r in group})
                mid, lo, hi = [], [], []
                for n in ns:
                    selected = [r for r in group if r["n"] == n]
                    mid.append(float(np.median([r[outcome + "_median"] for r in selected])))
                    lo.append(min(r[outcome + "_min"] for r in selected))
                    hi.append(max(r[outcome + "_max"] for r in selected))
                ax.plot(ns, mid, "o-", color=colors[shape], lw=1, ms=3, label=shape)
                ax.fill_between(ns, lo, hi, color=colors[shape], alpha=.09)
            ax.set(xscale="log", yscale="log", xlabel="Sequence length n (nt)", ylabel=ylabel, title=program)
            ax.set_xticks([250, 500, 1000, 2000, 4000], labels=["250", "500", "1000", "2000", "4000"])
            ax.grid(alpha=.2)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="outside lower center", ncol=4, fontsize=8)
        fig.suptitle(f"{revision}; cpu_40, workers={workers}, outer, b=h=d=100. Median across 3 seeds; band spans 9 measurements.", fontsize=10)
        fig.savefig(output / ("scaling_time.png" if outcome == "wall_seconds" else "scaling_rss.png"), dpi=170)
        plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=True)
    for col, program in enumerate(PROGRAMS):
        for ax, outcome, ylabel in zip(axes[:, col], OUTCOMES, ("CLI wall time (s)", "Peak RSS (GiB)")):
            for shape in series:
                group = [r for r in cpu if r["program"] == program and r["series"] == shape]
                ax.scatter([r["same_loop_pairs_turn"] for r in group], [r[outcome + "_median"] for r in group],
                           s=12, c=[colors[shape]], label=shape)
            ax.set(xscale="log", yscale="log", xlabel="Same-loop free endpoint pairs J_V (TURN=3)", ylabel=ylabel, title=program)
            ax.grid(alpha=.2)
    fig.suptitle(f"Synthetic inputs: scaffold descriptor and observed whole-CLI resource use ({revision})")
    fig.savefig(output / "feature_resources.png", dpi=170)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for program in PROGRAMS:
        group = sorted((r for r in rna if r["program"] == program), key=lambda r: r["n"])
        for ax, outcome, ylabel in zip(axes, OUTCOMES, ("CLI wall time (s)", "Peak RSS (GiB)")):
            mid = [r[outcome + "_median"] for r in group]
            ax.errorbar([r["n"] / 1000 for r in group], mid,
                yerr=[[v - r[outcome + "_min"] for v, r in zip(mid, group)],
                      [r[outcome + "_max"] - v for v, r in zip(mid, group)]], fmt="o-", lw=1, ms=4, label=program)
            ax.set(xlabel="Sequence length (kb)", ylabel=ylabel)
            ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    fig.suptitle(f"RNAcentral {revision}; cpu_40, workers={workers}, outer, b=h=d=100; median and range of 3 measurements", fontsize=10)
    fig.savefig(output / "rnacentral_time_rss.png", dpi=180)
    plt.close(fig)


def analyze(args):
    output = Path(args.out)
    if output.exists():
        raise ValueError("use a new E2 output directory")
    dataset, hashes = load_dataset(args.dataset)
    specifications = [(args.synthetic_run, "synthetic_current", dataset, (100, 100, 100), "cpu_40=1", args.workers),
        (args.rna_run, "rnacentral_current", {}, (100, 100, 100), "cpu_40=1", args.workers),
        (args.previous_rna_run, "rnacentral_previous" if args.skip_cap30_history else "rnacentral_15aa_cap100",
         {}, (100, 100, 100), args.previous_resource, 1)]
    if not args.skip_cap30_history:
        specifications.append((args.previous_rna_cap30_run, "rnacentral_15aa_cap30", {}, (100, 30, 30), "node_f=1", 1))
    if args.previous_synthetic_run:
        specifications.append((args.previous_synthetic_run, "synthetic_previous", dataset, (100, 100, 100), args.previous_resource, 1))
    cohorts, provenance, outcomes = [], {}, {}
    for path, label, data, caps, resource, workers in specifications:
        rows, p = load_run(path, label, data, caps, resource, expected_workers=workers)
        provenance[label] = {**p, "mode": RunReader(path).manifest["config"]["modes"][0]}
        cohorts.append(rows)
        if label.endswith("current"):
            if provenance[label]["mode"]["args"][-2:] != ["--ranking", "outer"]:
                raise ValueError("current run must explicitly select outer ranking")
            outcomes[label] = audit_outcomes(path)
    cpu, rna = cohorts[:2]
    if [len(c) for c in cohorts] != [468 if s[1].startswith("synthetic") else 21 for s in specifications]:
        raise ValueError("unexpected cohort counts")
    if provenance["synthetic_current"]["signature"] != provenance["rnacentral_current"]["signature"]:
        raise ValueError("current E2 software/model/output mismatch")
    models, ablation = fit_models(cpu)
    for program in PROGRAMS:
        group = [r for r in rna if r["program"] == program]
        for outcome in OUTCOMES:
            for name, model in models[program][outcome].items():
                prediction = predict(model, group)
                model["external_rnacentral"] = model_errors([r[outcome + "_median"] for r in group], prediction)
                for row, value in zip(group, prediction):
                    row[name + "_" + outcome + "_prediction"] = float(value)
                    row[name + "_outside_training_features"] = [k for k, (lo, hi) in model["training_range"].items()
                                                                 if not lo <= row[k] <= hi]
    pairs = geometry_comparisons(cpu, rna, models)
    for old, specification in zip(cohorts[2:], specifications[2:]):
        pairs.extend(historical_comparisons(old, cpu if specification[1].startswith("synthetic") else rna))
    if args.parallel_w1_run or args.parallel_w10_run:
        if not (args.parallel_w1_run and args.parallel_w10_run):
            raise ValueError("supply both parallel comparison runs")
        for path, label, workers in ((args.parallel_w1_run, "parallel_w1", 1), (args.parallel_w10_run, "parallel_w10", 10)):
            rows, p = load_run(path, label, {}, (100, 100, 100), "cpu_40=1", expected_workers=workers)
            if (len(rows) != 9 or p["signature"] != provenance["rnacentral_current"]["signature"]
                    or p["mode"] != provenance["rnacentral_current"]["mode"]):
                raise ValueError("parallel comparison cohort/software mismatch")
            outcomes[label] = audit_outcomes(path)
            cohorts.append(rows)
            provenance[label] = p
        w1, w10 = cohorts[-2:]
        names = {r["name"] for r in w1}
        production = [r for r in rna if r["name"] in names]
        for group, kind in ((w10, "parallel_w10_over_w1"), (production, "production_over_w1")):
            for pair in historical_comparisons(w1, group):
                pair["kind"] = kind
                pairs.append(pair)
    variability = [{"run": label, "program": program, "outcome": outcome,
        "relative_range": distribution([r[outcome + "_relative_range"] for r in rows if r["program"] == program]),
        "over_5_percent": sum(r[outcome + "_relative_range"] > .05 for r in rows if r["program"] == program)}
        for label, rows in zip(provenance, cohorts) for program in PROGRAMS for outcome in OUTCOMES]
    correlations = [{"run": label, "program": program, "feature": f, "outcome": outcome,
        "spearman": spearman([r[f] for r in rows if r["program"] == program],
                             [r[outcome + "_median"] for r in rows if r["program"] == program])}
        for label, rows in zip(list(provenance)[:2], cohorts[:2]) for program in PROGRAMS for f in FEATURES for outcome in OUTCOMES]
    summary = {"schema_version": 1, "scope": "Observed whole-CLI resources, not theoretical complexity or causal effects.",
        "definitions": {"rss_gib": "GNU time peak RSS KiB / 2^20", "relative_range": "(max-min)/median; not a confidence interval",
        "models": "ln(y)=c+sum beta ln(1+x); n_only OLS, other models standardized ridge alpha=1; medians of 3 repetitions",
        "slopes": "OLS ln(y) versus ln(n), per series and seed; finite-range observations, not asymptotic orders",
        "validation": "Refit within each training fold, leave one n/family/seed out; synthetic inputs only"},
        "outcome_counts": outcomes, "input_counts": {k: p["selected_inputs"] for k, p in provenance.items()},
        "variability": variability, "correlations": correlations, "models": models, "feature_block_ablation": ablation,
        "historical_comparison_caveat": (
            "Same S/G, b/h/d, outer ranking and requested outputs; version, workers, host and date differ. No implementation-only speedup claim."
            if args.skip_cap30_history else
            "Same S/G and b; version, ranking, allocation and date differ; cap30 also changes h,d. No isolated speedup claim."),
        "limitations": ["Three repetitions, fixed input cohort; S also changes between synthetic shapes.",
            "Correlated geometry descriptors do not identify causal work or RSS bytes.",
            "RNAcentral predictions extrapolate beyond synthetic n=250..4000, without refitting or resource calibration.",
            "CapR reports no numeric root; successful complete profiles use the CLI's finite-root check."],
        "provenance": {"runs": provenance, "dataset_sha256": hashes,
            "scripts_sha256": {name: file_digest(Path(__file__).with_name(name)) for name in
                (Path(__file__).name, "analyze_e0.py", "analyze_e2_shape.py", "analyze_e2_synthetic.py", "analyze_e2_diagnostics.py", "scaffold_features.py", "task_storage.py")},
            "tar_extractions_to_filesystem": 0}}
    summary["comparison_distributions"] = [{"kind": kind, "program": program, "outcome": outcome,
        **distribution([p[outcome + "_ratio"] for p in pairs if p["kind"] == kind and p["program"] == program])}
        for kind in sorted({p["kind"] for p in pairs}) for program in PROGRAMS for outcome in OUTCOMES]
    if args.e0_analysis:
        summarize_e0(args.e0_analysis, args.e0_run)
    output.mkdir(parents=True)
    write_union_tsv(output / "metrics.tsv", sum(cohorts, []))
    write_union_tsv(output / "comparisons.tsv", pairs)
    write_union_tsv(output / "slopes.tsv", empirical_slopes(cpu))
    save_json(output / "summary.json", summary)
    revision = provenance["rnacentral_current"]["signature"]["commit"][:8]
    figures(output, cpu, rna, revision, args.workers)
    if args.memory_run or args.previous_memory_run:
        if not (args.memory_run and args.previous_memory_run):
            raise ValueError("supply both memory breakdown runs")
        from analyze_e2_diagnostics import analyze_memory
        analyze_memory(output, args.memory_run, args.previous_memory_run, args.workers)
    print(json.dumps({"metrics_rows": sum(map(len, cohorts)), "comparison_rows": len(pairs),
        "raw_artifacts_verified": sum(p["raw_artifacts_verified"] for p in provenance.values()), "output": str(output)}))


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", default=str(root / "data/e2-synthetic-shapes-20261012"))
    parser.add_argument("--e0-run", default=str(root / "results/runs/e0-bprna-legacy-d21618c5-20261014"))
    parser.add_argument("--e0-analysis", help="Omit for E2-only analysis; existing E0 output stays untouched")
    parser.add_argument("--synthetic-run", default=str(root / "results/runs/e2-synthetic-shapes-d21618c5-20261014"))
    parser.add_argument("--rna-run", default=str(root / "results/runs/e2-rnacentral-d21618c5-20261014"))
    parser.add_argument("--previous-rna-run", default=str(root / "results/runs/e2-rnacentral-capb-timing-20261012"))
    parser.add_argument("--previous-rna-cap30-run", default=str(root / "results/runs/e2-rnacentral-15aa7742-timing-20261011"))
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--previous-resource", default="node_f=1")
    parser.add_argument("--skip-cap30-history", action="store_true")
    parser.add_argument("--previous-synthetic-run")
    parser.add_argument("--memory-run")
    parser.add_argument("--previous-memory-run")
    parser.add_argument("--parallel-w1-run")
    parser.add_argument("--parallel-w10-run")
    parser.add_argument("--out", required=True)
    analyze(parser.parse_args())


if __name__ == "__main__":
    main()
