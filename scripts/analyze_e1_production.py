#!/usr/bin/env python3
"""E1 packed production arrays: separate exact run, or widest-beam sensitivity.

Use the unchanged analyze_accuracy definitions. Never extract tar members or
execute RNA programs. One temporary SQLite file in /tmp joins differently packed
exact/beam tasks; shared outputs are three TSVs, three JSONs and two PNGs per cohort.
"""

import argparse
from collections import Counter
import csv
import gc
import json
import math
from pathlib import Path
import pickle
import re
import sqlite3
import tempfile
import time

import numpy as np

import analyze_accuracy as accuracy
from analyze_bprna_long import audit_job, threshold_counts
from analyze_bprna_short import ACCESS, CHECKS, CONTEXT, CORE, plan_provenance, statistics
from bprna_capb import check_task, resource_fields
from run_cli import digest, git_state, now
from task_storage import FORMAT, RunReader, input_records

SHORT = (("nt50_99", 50, 99), ("nt100_199", 100, 199),
         ("nt200_299", 200, 299), ("nt300_399", 300, 399))
LONG = (("nt400_799", 400, 799), ("nt800_1599", 800, 1599),
        ("nt1600_3199", 1600, 3199), ("nt3200_4381", 3200, 4381))
WIDTHS = (1, 4, 8, 16)
RESOURCES = sorted(f"{p}_{k}" for p in accuracy.PROGRAMS for k in
                   ("elapsed_seconds", "time_wall_seconds", "peak_rss_kib"))
METRICS = CORE + ["abs_delta_free_energy_kcal_mol", "delta_log_z"] + CONTEXT + [
    k for k in ACCESS if "mae" in k or "max_error" in k]
COUNTS = [k for k in ACCESS if k not in METRICS] + [f"access_w{w}_all_windows" for w in WIDTHS]
NUMERIC = ["length", "free_bases", "rt_kcal_mol", "free_energy_kcal_mol",
           "reference_free_energy_kcal_mol"] + METRICS + COUNTS + CHECKS + RESOURCES
FIELDS = ["input", "input_sha256", "tier", "task", "mode", "beam", "reference_mode",
          "comparison_kind", "analysis_status"] + NUMERIC


def load_plan(run):
    plan = json.loads((run / "array_plan.json").read_text())
    if plan.get("output_format") != FORMAT:
        raise ValueError("requires task-jsonl-tar-v1")
    if digest(run / "inputs.tar") != plan["inputs_archive_sha256"]:
        raise ValueError("input archive checksum mismatch")
    return plan


def cohort(plan, tiers):
    allowed = {t[0] for t in tiers}
    tasks = [t for t in plan["tasks"] if t["tier"] in allowed]
    entries = [i for t in tasks for i in t["inputs"]]
    hashes = {i["name"]: i["sha256"] for i in entries}
    if len(hashes) != len(entries) or hashes != {i["name"]: i["sha256"] for i in plan["inputs"]}:
        raise ValueError("duplicate, omitted or out-of-cohort input")
    return tasks, hashes


def check_reference_plans(beam, exact):
    # The author verified this specific outside-only revision on 14 inputs.
    sources = tuple(p["software_before"]["commit"] for p in (beam, exact))
    if sources[0] != sources[1] and sources != (
            "d21618c5115453174204329fcc40534b1cb9aa9c",
            "d9804df6a6a20590d6762475878d24a14c4e8239"):
        raise ValueError("unverified beam/exact source combination")
    for key in ("model_settings", "programs"):
        if beam["config"][key] != exact["config"][key]:
            raise ValueError("beam/exact model or requested outputs differ")
    if beam["software_lock"]["energy_parameter_sha256"] != exact["software_lock"]["energy_parameter_sha256"]:
        raise ValueError("beam/exact energy parameters differ")
    if sources[0] == sources[1] and accuracy.model_signature(beam) != accuracy.model_signature(exact):
        raise ValueError("same-source beam/exact binaries differ")


def check_command(job, config):
    mode = next(m for m in config["modes"] if m["name"] == job["mode"])
    program = next(p for p in config["programs"] if p["name"] == job["program"])
    args = accuracy.algorithm_arguments(job["command"])
    if args != mode["args"] + program.get("args", []):
        raise ValueError("raw command differs from planned condition")
    if job["mode"] == "exact":
        if mode["args"] != ["--exact"]:
            raise ValueError("exact reference must have no local limits")
    else:
        match = re.fullmatch(r"beam(\d+)_outer", job["mode"])
        if not match:
            raise ValueError("requires outer beam modes")
        beam = int(match[1])
        for flag, value in (("-b", str(beam)), ("--max-hairpin", str(max(beam, 30))),
                            ("--max-case3-depth", str(max(beam, 30))), ("--ranking", "outer")):
            if args.count(flag) != 1 or args[args.index(flag) + 1] != value:
                raise ValueError("condition violates outer/max(b,30) rule")
    if job["program"] == "LinearCParty":
        if args.count("--bpp-cutoff") != 1 or args[args.index("--bpp-cutoff") + 1] != "0" or job["command"].count("--bpp") != 1:
            raise ValueError("requires unthresholded BPP export")


def read_group(reader, jobs, name, mode, record, provenance):
    parsed = {}
    for program in sorted(accuracy.PROGRAMS):
        job = jobs[name, mode, program]
        if job["status"] != "success":
            raise ValueError("selected job did not succeed: " + job["id"])
        check_command(job, reader.manifest["config"])
        audit_job(reader, job)
        parsed[program] = accuracy.read_job(reader, job, record)
    if parsed["LinearPKRaccess"]["widths"] != list(WIDTHS):
        raise ValueError("unexpected accessibility window lengths")
    check = accuracy.check_consistency(name, mode, 0, parsed)
    provenance["consistency_checks"][check["status"]] += 1
    for key in CHECKS:
        provenance["consistency_maxima"][key] = max(provenance["consistency_maxima"][key], check[key])
    if check["status"] != "success":
        raise ValueError("CLI consistency failed: " + str(check))
    return parsed, {k: check[k] for k in CHECKS}


def open_task(run, plan, task, modes, provenance):
    reader = RunReader(run / task["directory"] / "run")
    hashes, all_jobs = check_task(reader, plan, task)
    jobs = {k: v for k, v in all_jobs.items() if k[1] in modes}
    if set(jobs) != {(n, m, p) for n in hashes for m in modes for p in accuracy.PROGRAMS}:
        raise ValueError("incomplete selected task matrix")
    # Old exact jobs are deliberately excluded, including their statuses/artifacts.
    provenance["selected_job_status_counts"].update(j["status"] for j in jobs.values())
    provenance["excluded_job_status_counts"].update(j["status"] for k, j in all_jobs.items() if k not in jobs)
    records = input_records(reader)
    lo, hi = next((lo, hi) for tier, lo, hi in SHORT + LONG if tier == task["tier"])
    if any(not lo <= len(r[1]) <= hi for r in records.values()):
        raise ValueError("RNA length outside tier")
    provenance["tasks"].append({"run": str(run), "tier": task["tier"], "task": task["task_id"],
        "run_manifest_sha256": reader.manifest_sha256,
        "checkpoint_sha256": digest(reader.path / "task_results.jsonl"), "inputs": len(records),
        "selected_jobs": len(jobs), "excluded_jobs": len(all_jobs) - len(jobs),
        "hostname": reader.manifest["host"]["hostname"], "workers": reader.manifest["workers"]})
    return reader, hashes, jobs, records


def compare(parsed, reference, record):
    pair = accuracy.compare_bpp(parsed["LinearCParty"], reference["LinearCParty"], record)
    context = accuracy.compare_contexts(parsed["LinearPKCapR"], reference["LinearPKCapR"], record)
    access = accuracy.compare_accessibility(parsed["LinearPKRaccess"], reference["LinearPKRaccess"], record)
    row = {k: pair[k] for k in CORE if k in pair}
    row.update({k: access[k] for k in CORE if k in access})
    row.update(free_bases=pair["free_bases"], context_tv_mean_Gfree=context["context_tv_mean_Gfree"],
               abs_delta_free_energy_kcal_mol=abs(access["delta_free_energy_kcal_mol"]),
               delta_log_z=access["delta_log_z"], rt_kcal_mol=access["rt_kcal_mol"],
               free_energy_kcal_mol=parsed["LinearPKRaccess"]["free_energy"],
               reference_free_energy_kcal_mol=reference["LinearPKRaccess"]["free_energy"])
    row.update({"context_" + c + "_mae_Gfree": v for c, v in context["context_mae_Gfree"].items()})
    for w, values in access["accessibility"].items():
        row.update({f"access_w{w}_{k}": v for k, v in values.items()})
    return row


def write_metrics(run, exact_run, reference_mode, output):
    start = time.monotonic()
    plan = load_plan(run)
    tiers = SHORT if exact_run else LONG
    beams = (20, 50, 100, 200) if exact_run else (20, 50, 100, 200, 300, 400)
    modes = [f"beam{b}_outer" for b in beams]
    tasks, hashes = cohort(plan, tiers)
    if set(m["name"] for m in plan["config"]["modes"]) != set(modes + (["exact"] if exact_run else [])):
        raise ValueError("unexpected beam-run modes")
    reference_mode = "exact" if exact_run else reference_mode
    if not exact_run and reference_mode != modes[-1]:
        raise ValueError("long reference must be the maximum measured width: beam400_outer")
    exact_plan = None
    if exact_run:
        exact_plan = load_plan(exact_run)
        check_reference_plans(plan, exact_plan)
        exact_tasks, exact_hashes = cohort(exact_plan, tiers)
        if hashes != exact_hashes or [m["name"] for m in exact_plan["config"]["modes"]] != ["exact"]:
            raise ValueError("exact run must contain the identical S/G cohort and only exact")
    kind = "exact_error" if exact_run else "beam_width_sensitivity"
    provenance = {"schema_version": 1, "generated_utc": now(), "comparison_kind": kind,
        "reference_mode": reference_mode, "beams": beams, "bins": tiers,
        "runs": [plan_provenance(run, plan)] + ([plan_provenance(exact_run, exact_plan)] if exact_run else []),
        "analysis_git_state": git_state(accuracy.ROOT), "script_sha256": {n: digest(Path(__file__).with_name(n)) for n in
            (Path(__file__).name, "analyze_accuracy.py", "analyze_bprna_short.py", "analyze_bprna_long.py", "bprna_capb.py", "task_storage.py", "rna_records.py")},
        "storage": "tar member streams; no extraction; one temporary SQLite file under /tmp for exact join",
        "resource_scope": "reference distributions only; concurrent accuracy jobs; one repetition; no performance comparison",
        "tasks": [], "selected_job_status_counts": Counter(), "excluded_job_status_counts": Counter(),
        "consistency_checks": Counter(), "consistency_maxima": dict.fromkeys(CHECKS, 0.)}
    if exact_plan and plan["software_before"]["commit"] != exact_plan["software_before"]["commit"]:
        evidence = accuracy.ROOT / "results/local/exact-d9804df6-vs-d21618c5-20261015/summary.txt"
        provenance["cross_revision_evidence"] = {"path": str(evidence), "sha256": digest(evidence),
            "scope": "author-verified 14 RNA; beam b20/b200 identical; exact numeric difference <=2.2e-13"}
    output.mkdir(parents=True, exist_ok=False)
    seen = set()
    with tempfile.TemporaryDirectory(prefix="e1-production-", dir="/tmp") as scratch:
        db = sqlite3.connect(str(Path(scratch) / "exact.sqlite"))
        try:
            db.execute("CREATE TABLE refs (name TEXT PRIMARY KEY, hash TEXT, record BLOB, parsed BLOB)")
            if exact_run:
                for number, task in enumerate(exact_tasks, 1):
                    reader, task_hashes, jobs, records = open_task(exact_run, exact_plan, task, ["exact"], provenance)
                    for name, record in records.items():
                        parsed, _ = read_group(reader, jobs, name, "exact", record, provenance)
                        db.execute("INSERT INTO refs VALUES (?,?,?,?)", (name, task_hashes[name],
                            pickle.dumps(record), pickle.dumps(parsed)))
                    db.commit()
                    if number % 5 == 0 or number == len(exact_tasks):
                        print(f"exact {number}/{len(exact_tasks)} tasks; {time.monotonic()-start:.1f}s", flush=True)
                    del reader, jobs, records, parsed
                    gc.collect()
            with (output / "metrics.tsv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, FIELDS, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                for number, task in enumerate(tasks, 1):
                    reader, task_hashes, jobs, records = open_task(run, plan, task, modes, provenance)
                    for name, record in records.items():
                        if name in seen:
                            raise ValueError("duplicate beam RNA")
                        seen.add(name)
                        if exact_run:
                            saved_hash, saved_record, saved_parsed = db.execute("SELECT hash,record,parsed FROM refs WHERE name=?", (name,)).fetchone()
                            if saved_hash != task_hashes[name] or pickle.loads(saved_record) != record:
                                raise ValueError("beam/exact S/G mismatch")
                            reference = pickle.loads(saved_parsed)
                        else:
                            reference, reference_check = read_group(reader, jobs, name, reference_mode, record, provenance)
                        for mode, beam in zip(modes, beams):
                            if mode == reference_mode:
                                parsed, check = reference, reference_check
                            else:
                                parsed, check = read_group(reader, jobs, name, mode, record, provenance)
                            row = {"input": name, "input_sha256": task_hashes[name], "tier": task["tier"],
                                "task": task["task_id"], "mode": mode, "beam": beam, "length": len(record[1]),
                                "reference_mode": reference_mode, "comparison_kind": kind, "analysis_status": "success",
                                **compare(parsed, reference, record), **check, **resource_fields(jobs, name, mode)}
                            writer.writerow({k: "NA" if v is None else v for k, v in row.items()})
                        del parsed, reference
                    if number % 5 == 0 or number == len(tasks):
                        print(f"beam {number}/{len(tasks)} tasks; {len(seen)} RNA; {time.monotonic()-start:.1f}s", flush=True)
                    del reader, jobs, records
                    gc.collect()
        finally:
            db.close()
    if seen != set(hashes):
        raise ValueError("beam cohort incomplete")
    provenance.update(inputs=len(seen), metric_rows=len(seen)*len(beams),
                      metrics_sha256=digest(output / "metrics.tsv"), elapsed_seconds=time.monotonic()-start,
                      validation="selected matrix, unchanged task dependencies, S/G hashes, raw commands/results/resources, used artifact hashes and CLI consistency passed")
    write_json(output / "provenance.json", provenance)
    return provenance


def write_json(path, payload):
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


def summarize(output, make_plots=True):
    provenance = json.loads((output / "provenance.json").read_text())
    if digest(output / "metrics.tsv") != provenance["metrics_sha256"]:
        raise ValueError("metrics checksum mismatch")
    data = np.empty((provenance["metric_rows"], len(NUMERIC)))
    bins, beams, names = [], [], []
    with (output / "metrics.tsv").open() as stream:
        for i, row in enumerate(csv.DictReader(stream, delimiter="\t")):
            data[i] = [float(row[k]) if row[k] != "NA" else np.nan for k in NUMERIC]
            bins.append(row["tier"]); beams.append(int(row["beam"])); names.append(row["input"])
    if len(names) != provenance["metric_rows"]:
        raise ValueError("metrics row count differs")
    bins, beams = np.asarray(bins), np.asarray(beams)
    groups = []
    for tier in ["all"] + [t[0] for t in provenance["bins"]]:
        for beam in provenance["beams"]:
            mask = (beams == beam) & ((bins == tier) if tier != "all" else True)
            block = data[mask]
            n = int(mask.sum())
            metrics = {k: {**statistics(block[:, NUMERIC.index(k)]),
                           "missing_or_out_of_domain_rna": int(n - np.isfinite(block[:, NUMERIC.index(k)]).sum())}
                       for k in METRICS}
            counts = {k: {"total_windows": int(np.nansum(block[:, NUMERIC.index(k)])),
                           "rna_with_positive_count": int((block[:, NUMERIC.index(k)] > 0).sum()),
                           "missing_rna": int(np.isnan(block[:, NUMERIC.index(k)]).sum())} for k in COUNTS}
            groups.append({"bin": tier, "beam": beam, "rna": n, "metrics": metrics, "window_counts": counts,
                "delta_F_thresholds": threshold_counts(block[:, NUMERIC.index(CORE[0])]),
                "no_free_bases_rna": int((block[:, NUMERIC.index("free_bases")] == 0).sum())})
    summary = {"schema_version": 1, "generated_utc": now(),
        "metrics_sha256": provenance["metrics_sha256"], "summary_script_sha256": digest(Path(__file__)),
        "reference_mode": provenance["reference_mode"],
        "comparison_kind": provenance["comparison_kind"], "inputs": provenance["inputs"],
        "unit": "one RNA; equal weight; per-RNA means before RNA quantiles",
        "quantile_method": "linear interpolation; NA excluded; zeros retained",
        "delta_F_definition": "signed F_candidate - F_reference, kcal/mol; never clipped",
        "opening_domain": "both probabilities positive; reference-positive/candidate-zero windows counted separately",
        "missing_definition": "missing_or_out_of_domain_rna includes empty eligible/both-positive domains; not zero-filled",
        "groups": groups}
    write_json(output / "summary.json", summary)
    with (output / "quantiles.tsv").open("w", newline="") as stream:
        fields = ["bin", "beam", "metric", "rna", "n", "missing_or_out_of_domain_rna",
                  "median", "p90", "p99", "max", "min", "mean"]
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for group in groups:
            for metric, values in group["metrics"].items():
                writer.writerow({k: "NA" if v is None else v for k, v in
                    dict(bin=group["bin"], beam=group["beam"], rna=group["rna"], metric=metric, **values).items()})
    with (output / "window_counts.tsv").open("w", newline="") as stream:
        fields = ["bin", "beam", "count", "rna", "total_windows", "rna_with_positive_count", "missing_rna"]
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for group in groups:
            for count, values in group["window_counts"].items():
                writer.writerow(dict(bin=group["bin"], beam=group["beam"], rna=group["rna"], count=count, **values))
    resources = {"scope": provenance["resource_scope"], "groups": [
        {"bin": tier, "beam": beam, "distributions": {k: statistics(data[mask, NUMERIC.index(k)]) for k in RESOURCES}}
        for tier in ["all"] + [t[0] for t in provenance["bins"]] for beam in provenance["beams"]
        for mask in [(beams == beam) & ((bins == tier) if tier != "all" else True)]]}
    write_json(output / "resource_reference.json", resources)
    if make_plots:
        plot_summary(summary, provenance, output)
    return summary


def plot_summary(summary, provenance, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    selections = [("core_quantiles.png", [CORE[0], "base_state_tv_mean_Gfree", "context_tv_mean_Gfree"],
                   ["Signed delta F (kcal/mol)", "BPP TV on free bases", "Loop-profile TV on free bases"]),
                  ("accessibility_quantiles.png", [f"access_w{w}_probability_mae_eligible" for w in WIDTHS] +
                   [f"access_w{w}_opening_mae_both_positive_kcal_mol" for w in WIDTHS],
                   [f"Accessibility MAE\nw={w}" for w in WIDTHS] + [f"Opening MAE (kcal/mol)\nw={w}" for w in WIDTHS])]
    for filename, metrics, labels in selections:
        fig, axes = plt.subplots(len(metrics), 4, figsize=(15, 2.4*len(metrics)), squeeze=False)
        for j, (tier, lo, hi) in enumerate(provenance["bins"]):
            rows = [r for r in summary["groups"] if r["bin"] == tier]
            for i, (metric, label) in enumerate(zip(metrics, labels)):
                ax = axes[i,j]
                for q, style in (("median", "-o"), ("p90", "-s"), ("p99", "-^"), ("max", ":")):
                    ax.plot([r["beam"] for r in rows], [np.nan if r["metrics"][metric][q] is None else r["metrics"][metric][q] for r in rows], style, label=q, markersize=3)
                ax.set_yscale("symlog", linthresh=1e-10)
                ax.yaxis.get_major_locator().set_params(numticks=6)
                ax.grid(alpha=.2)
                if i == 0: ax.set_title(f"{lo}-{hi} nt")
                if j == 0: ax.set_ylabel(label)
                if i == len(metrics)-1: ax.set_xlabel("Beam width")
        axes[0,-1].legend(fontsize=8)
        reference = "error relative to exact" if summary["comparison_kind"] == "exact_error" else "differences relative to b=400"
        fig.suptitle(f"E1: {reference}\nouter; h=d=max(b,30); linear y within |value| <= 1e-10", fontsize=12)
        fig.tight_layout(rect=(0,0,1,.94 if len(metrics) == 3 else .975))
        fig.savefig(output / filename, dpi=150)
        plt.close(fig)


def arguments(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cohort", choices=("short", "long"))
    parser.add_argument("--beam-run", type=Path)
    parser.add_argument("--exact-run", type=Path)
    parser.add_argument("--run", type=Path)
    parser.add_argument("--reference", default="beam400_outer")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--summary-only", action="store_true", help="recreate summaries/figures from an unchanged metrics TSV")
    args = parser.parse_args(argv)
    if not args.summary_only:
        if args.cohort == "short" and (not args.beam_run or not args.exact_run or args.run):
            parser.error("short requires --beam-run and --exact-run")
        if args.cohort == "long" and (not args.run or args.beam_run or args.exact_run or args.reference != "beam400_outer"):
            parser.error("long requires --run and --reference beam400_outer")
    return args


def main(argv=None):
    args = arguments(argv)
    if not args.summary_only:
        write_metrics(args.beam_run or args.run, args.exact_run, args.reference, args.out)
    summarize(args.out)


if __name__ == "__main__":
    main()
