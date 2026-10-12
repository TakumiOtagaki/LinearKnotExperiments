#!/usr/bin/env python3
"""Task-bounded bpRNA E1 analysis using analyze_accuracy's validated metrics.

No tar extraction, per-job output files, CLI measurements, or worker threads.
One row per RNA and finite-beam condition; every exact replacement is audited.
"""

import csv
from collections import Counter, defaultdict
import gc
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import analyze_accuracy as accuracy
from rna_records import describe_scaffold, scaffold_pairs
from run_cli import digest, git_state, now
from task_storage import FORMAT, RunReader, input_records

TIERS = ("nt50_99", "nt100_199", "nt200_299", "nt300_399")
BEAMS = (20, 50, 100, 200)
CAPS = ("default", "nocap")
IDENTITY = ["input", "tier", "task", "mode", "cap", "beam", "analysis_status"]
DESCRIPTORS = ["length", "scaffold_pairs", "free_bases", "paired_fraction", "maximum_pair_depth",
               "maximum_branching_ancestors_of_free_base", "branching_nodes_including_exterior",
               "gc_fraction", "exact_crossing_pair_count_gt_1e_minus12", "exact_crossing_pair_probability_sum",
               "exact_pseudoloop_mean_Gfree"]
CORE = ["delta_free_energy_kcal_mol", "partition_function_retention", "delta_free_energy_per_nt",
        "generated_pair_l1_sum", "base_state_tv_mean_Gfree", "unpaired_mae_Gfree", "context_tv_mean_Gfree"]
CONTEXT = ["context_" + c + "_mae_Gfree" for c in accuracy.CONTEXTS]
ACCESS = [f"access_w{w}_{k}" for w in (1, 4, 8, 16) for k in
          ("eligible_windows", "probability_mae_eligible", "probability_max_error_eligible",
           "opening_mae_both_positive_kcal_mol", "both_positive", "both_zero",
           "reference_positive_candidate_zero", "reference_zero_candidate_positive")]
CHECKS = ["stem_vs_bpp_max_error", "accessibility_w1_vs_unpaired_max_error",
          "accessibility_above_marginal_max_violation"]
REFERENCES = ["exact_LinearCParty_rerun", "exact_LinearPKCapR_rerun", "exact_LinearPKRaccess_rerun",
              "rt_kcal_mol", "exact_free_energy_kcal_mol"]
RESOURCES = [f"{p}_{k}" for p in accuracy.PROGRAMS for k in
             ("elapsed_seconds", "time_wall_seconds", "peak_rss_kib")]
FIELDS = IDENTITY + DESCRIPTORS + CORE + CONTEXT + ACCESS + CHECKS + REFERENCES + RESOURCES
METRICS = CORE + CONTEXT + [key for key in ACCESS if "mae" in key or "max_error" in key]


def measured_modes(caps):
    return [f"beam{b}_{c}" for c in caps for b in BEAMS if c != "capb" or b != 20]


def crossing_descriptors(record, bpp, contexts):
    # A generated pair crosses G iff its endpoints have different enclosing-G masks.
    masks, mask, stack = [], 0, []
    for i, symbol in enumerate(record[2]):
        if symbol == "(":
            stack.append(1 << i)
            mask |= stack[-1]
        elif symbol == ")":
            mask ^= stack.pop()
        masks.append(mask)
    occupied = {i for pair in scaffold_pairs(record[2]) for i in pair}
    cross = [p for (i, j), p in bpp["pairs"].items()
             if i not in occupied and j not in occupied and masks[i - 1] != masks[j - 1]]
    free = [i for i, c in enumerate(record[2]) if c == "."]
    return {"exact_crossing_pair_count_gt_1e_minus12": sum(p > 1e-12 for p in cross),
            "exact_crossing_pair_probability_sum": sum(cross),
            "exact_pseudoloop_mean_Gfree": accuracy.mean([contexts["PseudoLoop"][i] for i in free])}


def plan_provenance(path, plan):
    return {"path": str(path), "array_plan_sha256": digest(path / "array_plan.json"),
            "inputs_archive_sha256": plan.get("inputs_archive_sha256"),
            "software_commit": plan["software_before"]["commit"],
            "software_lock": plan["software_lock"], "config": {k: v for k, v in plan["config"].items() if k != "inputs"}}


def write_metrics(run, rerun_paths, output):
    start = time.monotonic()
    plan = json.loads((run / "array_plan.json").read_text())
    if plan.get("output_format") != FORMAT:
        raise ValueError("packed short-array analysis requires task-jsonl-tar-v1")
    reruns = accuracy.ExactReruns([RunReader(p) for p in rerun_paths], plan)
    tasks = [t for t in plan["tasks"] if t["tier"] in TIERS]
    seen, status_counts, analysis_counts, checks = set(), Counter(), Counter(), Counter()
    provenance = {"schema_version": 1, "generated_utc": now(), "main_run": plan_provenance(run, plan),
                  "analysis_git_state": git_state(accuracy.ROOT),
                  "script_sha256": {p.name: digest(p) for p in
                      (Path(__file__), Path(accuracy.__file__), Path(__file__).with_name("task_storage.py"),
                       Path(__file__).with_name("rna_records.py"))},
                  "exact_rerun_runs": reruns.provenance, "tasks": [],
                  "scope": "50 <= n < 400; >=400 planned but excluded; one RNA is one observation",
                  "storage": "task-bounded reads of tar member streams; no extraction; one metrics TSV",
                  "threads": 1}
    output.mkdir(parents=True, exist_ok=False)
    with (output / "metrics.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, FIELDS, delimiter="\t", lineterminator="\n", restval="NA")
        writer.writeheader()
        for task_number, task in enumerate(tasks, 1):
            path = run / task["directory"] / "run"
            reader = RunReader(path)
            saved = reader.manifest
            if saved["status"] not in ("success", "failed"):
                raise ValueError("short task is incomplete: " + str(path))
            if accuracy.model_signature(saved) != accuracy.model_signature(plan):
                raise ValueError("task source/model/binaries differ from plan")
            if (saved.get("software_changed_during_run") is not False or
                    saved.get("binaries_changed_during_run") is not False or
                    saved.get("software_lock_valid_after_run") is not True):
                raise ValueError("task provenance does not establish unchanged dependencies")
            records = input_records(reader)
            expected_inputs = {i["name"]: i["sha256"] for i in task["inputs"]}
            actual_inputs = {i["name"]: i["sha256"] for i in saved["inputs"]}
            if actual_inputs != expected_inputs:
                raise ValueError("task input hashes differ from array plan")
            by_input = defaultdict(dict)
            for job in saved["jobs"]:
                if job["warmup"] or job["repetition"] != 0:
                    raise ValueError("short E1 requires no warmups and one repetition")
                key = (job["mode"], job["program"])
                if key in by_input[job["input"]]:
                    raise ValueError("duplicate E1 observation")
                by_input[job["input"]][key] = job
                status_counts[job["status"]] += 1
            task_info = {"path": str(path), "run_manifest_sha256": reader.manifest_sha256,
                         "task_results_sha256": digest(path / "task_results.jsonl"),
                         "status": saved["status"], "inputs": len(records), "jobs": len(saved["jobs"])}
            provenance["tasks"].append(task_info)
            for name in records:
                if name in seen:
                    raise ValueError("RNA assigned to multiple tasks")
                seen.add(name)
                record = records[name]
                if not 50 <= len(record[1]) < 400:
                    raise ValueError("length outside requested short cohort")
                jobs = by_input[name]
                modes = ["exact"] + [f"beam{b}_{c}" for c in CAPS for b in BEAMS]
                if set(jobs) != {(mode, program) for mode in modes for program in accuracy.PROGRAMS}:
                    raise ValueError("incomplete planned E1 matrix for " + name)
                parsed, replaced, errors = {}, {}, {}
                for key, job in jobs.items():
                    source, selected = reruns.replace(reader, job, record, actual_inputs[name])
                    replaced[key] = int(source is not reader)
                    if selected["status"] != "success":
                        errors[key] = "job_failed: " + selected["status"]
                        continue
                    try:
                        parsed[key] = accuracy.read_job(source, selected, record)
                    except (ValueError, KeyError, OSError) as exc:
                        errors[key] = "artifact_error: " + str(exc)
                desc = describe_scaffold(record[1], record[2])
                if all(("exact", p) in parsed for p in accuracy.PROGRAMS):
                    desc.update(crossing_descriptors(record, parsed["exact", "LinearCParty"],
                                                     parsed["exact", "LinearPKCapR"]))
                for mode in modes:
                    group = {p: parsed[mode, p] for p in accuracy.PROGRAMS if (mode, p) in parsed}
                    if len(group) == 3:
                        check = accuracy.check_consistency(name, mode, 0, group)
                        checks[check["status"]] += 1
                        if check["status"] != "success":
                            provenance.setdefault("consistency_failures", []).append(check)
                    if mode == "exact":
                        continue
                    match = accuracy.re.fullmatch(r"beam(\d+)_(default|nocap)", mode)
                    row = {"input": name, "tier": task["tier"], "task": task["task_id"],
                           "mode": mode, "cap": match[2], "beam": int(match[1]), **desc}
                    for p in accuracy.PROGRAMS:
                        row["exact_" + p + "_rerun"] = replaced["exact", p]
                    needed = [(m, p) for m in ("exact", mode) for p in accuracy.PROGRAMS]
                    if any(key not in parsed for key in needed):
                        row["analysis_status"] = "unavailable"
                        provenance.setdefault("analysis_errors", []).append(
                            {"input": name, "mode": mode, "errors": {str(k): errors[k] for k in needed if k in errors}})
                    else:
                        pair = accuracy.compare_bpp(parsed[mode, "LinearCParty"], parsed["exact", "LinearCParty"], record)
                        context = accuracy.compare_contexts(parsed[mode, "LinearPKCapR"], parsed["exact", "LinearPKCapR"], record)
                        access = accuracy.compare_accessibility(parsed[mode, "LinearPKRaccess"], parsed["exact", "LinearPKRaccess"], record)
                        row.update({k: pair[k] for k in CORE if k in pair})
                        row.update({k: access[k] for k in CORE if k in access})
                        row["context_tv_mean_Gfree"] = context["context_tv_mean_Gfree"]
                        row.update({"context_" + c + "_mae_Gfree": value for c, value in context["context_mae_Gfree"].items()})
                        for w, metrics in access["accessibility"].items():
                            row.update({f"access_w{w}_{k}": v for k, v in metrics.items() if f"access_w{w}_{k}" in ACCESS})
                        row.update({k: check[k] for k in CHECKS})
                        row.update(rt_kcal_mol=access["rt_kcal_mol"],
                                   exact_free_energy_kcal_mol=parsed["exact", "LinearPKRaccess"]["free_energy"],
                                   analysis_status="success")
                    analysis_counts[row["analysis_status"]] += 1
                    writer.writerow({k: "NA" if v is None else v for k, v in row.items()})
            print(f"task {task_number}/{len(tasks)} {task['tier']}/{task['task_id']}: {len(seen)} RNA, {time.monotonic()-start:.1f}s", flush=True)
            del reader, saved, records, by_input, parsed, jobs, group
            gc.collect()
    if reruns.used != set(reruns.jobs):
        raise ValueError("unused exact rerun jobs: " + str(sorted(set(reruns.jobs) - reruns.used)))
    provenance.update(inputs=len(seen), original_job_status_counts=dict(status_counts),
                      exact_replacements=reruns.replacements,
                      exact_replacement_program_counts=dict(Counter(r["program"] for r in reruns.replacements)),
                      effective_job_status_counts={"success": status_counts["success"] + len(reruns.used),
                                                   "failed": sum(v for k, v in status_counts.items() if k != "success") - len(reruns.used)},
                      analysis_status_counts=dict(analysis_counts), consistency_check_counts=dict(checks),
                      elapsed_seconds=time.monotonic()-start, metrics_sha256=digest(output / "metrics.tsv"))
    (output / "provenance.json").write_text(json.dumps(provenance, indent=2, allow_nan=False) + "\n")
    return provenance


def statistics(values):
    import numpy as np
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return {"n": 0, "median": None, "mean": None, "p90": None, "p99": None, "max": None, "min": None}
    q = np.quantile(values, [.5, .9, .99], method="linear")
    return {"n": len(values), "median": float(q[0]), "mean": float(np.mean(values)),
            "p90": float(q[1]), "p99": float(q[2]), "max": float(np.max(values)), "min": float(np.min(values))}


def monotonicity(matrix, direction, tolerance):
    import numpy as np
    valid = np.all(np.isfinite(matrix), axis=1)
    increase = direction * np.diff(matrix, axis=1)
    success = valid & np.all(increase <= tolerance, axis=1)
    return {"total_sequences": len(matrix), "eligible_sequences": int(valid.sum()),
            "monotone_sequences": int(success.sum()), "fraction_eligible": float(success.sum()/valid.sum()) if valid.any() else None,
            "nonmonotone_sequences": int((valid & ~success).sum()), "tolerance": tolerance}, increase, valid


def load_metrics(output):
    import numpy as np
    provenance = json.loads((output / "provenance.json").read_text())
    numeric = [key for key in FIELDS if key not in IDENTITY]
    caps = provenance.get("caps", CAPS)
    modes = measured_modes(caps)
    numbers = np.full((provenance["inputs"], len(caps), len(BEAMS), len(numeric)), np.nan)
    metadata, names = [], []
    with (output / "metrics.tsv").open() as stream:
        for i, row in enumerate(csv.DictReader(stream, delimiter="\t")):
            if i >= provenance["inputs"] * len(modes):
                raise ValueError("TSV row count differs from provenance")
            ci, bi = list(caps).index(row["cap"]), BEAMS.index(int(row["beam"]))
            numbers[i // len(modes), ci, bi] = [float(row[k]) if row.get(k, "NA") != "NA" else np.nan for k in numeric]
            if i % len(modes) == 0:
                names.append(row["input"])
                metadata.append({k: (None if row[k] == "NA" else float(row[k])) for k in DESCRIPTORS})
            expected = modes[i % len(modes)]
            if row["input"] != names[-1] or row["mode"] != expected:
                raise ValueError("TSV must keep each RNA's complete paired beam/cap block")
    if len(names) * len(modes) != provenance["inputs"] * len(modes) or i + 1 != provenance["inputs"] * len(modes):
        raise ValueError("TSV row count differs from provenance")
    if "capb" in caps:
        numbers[:, list(caps).index("capb"), 0] = numbers[:, list(caps).index("default"), 0]
    return provenance, {key: numbers[..., j] for j, key in enumerate(numeric)}, names, metadata


def summarize(output, make_plots=True):
    import numpy as np
    provenance, data, names, metadata = load_metrics(output)
    caps = list(provenance.get("caps", CAPS))
    if "capb" in caps:
        # The adopted policy uses the measured default b20 and measured capb b>=50.
        for key in data:
            data[key] = np.concatenate((data[key], data[key][:, caps.index("capb"):caps.index("capb")+1, :]), axis=1)
        caps.append("adopted")
    length = data["length"][:, 0, 0]
    bins = [("all", 50, 399)] + [(t, lo, hi) for t, lo, hi in zip(TIERS, (50,100,200,300), (99,199,299,399))]
    bins += [(f"nt{lo}_{lo+24}", lo, lo+24) for lo in range(50,400,25)]
    summary = {"schema_version": 1, "generated_utc": now(), "inputs": len(names),
               "quantiles": "RNA-weighted; numpy linear interpolation; NA excluded per metric with denominator retained",
               "metrics_sha256": provenance["metrics_sha256"], "caps": caps, "bins": [], "monotonicity": [],
               "condition_sources": {c: {str(b): f"beam{b}_{'default' if b == 20 and c in ('capb','adopted') else 'capb' if c == 'adopted' else c}" for b in BEAMS} for c in caps},
               "worst_examples": {}, "nonmonotone_examples": {}, "diagnostic_inputs": {}}
    for label, lo, hi in bins:
        selected = (length >= lo) & (length <= hi)
        group = {"bin": label, "min_length": lo, "max_length": hi, "sequences": int(selected.sum()),
                 "conditions": [], "cap_differences_default_minus_nocap": []}
        for ci, cap in enumerate(caps):
            for bi, beam in enumerate(BEAMS):
                metrics = {k: statistics(data[k][selected, ci, bi]) for k in METRICS}
                delta = data[CORE[0]][selected, ci, bi]
                condition = {"cap": cap, "beam": beam, "metrics": metrics,
                             "delta_F_above": {str(threshold): {"count": int(np.sum(delta > threshold)),
                                 "eligible_sequences": int(np.isfinite(delta).sum()), "total_sequences": int(selected.sum()),
                                 "fraction_eligible": float(np.mean(delta[np.isfinite(delta)] > threshold)) if np.isfinite(delta).any() else None}
                                 for threshold in (.01, .1, 1)},
                             "delta_F_below_minus_1e8": int(np.sum(delta < -1e-8)),
                             "window_counts": {k: {"NA_sequences": int(np.isnan(data[k][selected,ci,bi]).sum()),
                                                     "sum": float(np.nansum(data[k][selected,ci,bi]))} for k in ACCESS if k not in METRICS}}
                group["conditions"].append(condition)
            for key, direction, tolerance in [(CORE[0], 1, 1e-8), (CORE[1], -1, 1e-8),
                                               (CORE[4], 1, 1e-9), (CORE[6], 1, 1e-9)]:
                stat, _, _ = monotonicity(data[key][selected,ci,:], direction, tolerance)
                summary["monotonicity"].append({"bin": label, "cap": cap, "metric": key, **stat})
        for bi, beam in enumerate(BEAMS):
            values = {}
            for key in METRICS:
                diff = data[key][selected,caps.index("default"),bi] - data[key][selected,caps.index("nocap"),bi]
                tolerance = 1e-8 if key.startswith("delta_free_energy") else 1e-9
                finite = diff[np.isfinite(diff)]
                values[key] = {**statistics(diff), "default_larger_count": int(np.sum(finite > tolerance)),
                               "nocap_larger_count": int(np.sum(finite < -tolerance)), "comparison_tolerance": tolerance}
            group["cap_differences_default_minus_nocap"].append({"beam": beam, "metrics": values})
        if "capb" in caps:
            group["paired_cap_differences"] = []
            for right in ("nocap", "default"):
                for bi, beam in enumerate(BEAMS):
                    values = {}
                    for key in METRICS:
                        diff = data[key][selected,caps.index("capb"),bi] - data[key][selected,caps.index(right),bi]
                        tolerance = 1e-8 if key.startswith("delta_free_energy") or "kcal_mol" in key else 1e-9
                        finite = diff[np.isfinite(diff)]
                        values[key] = {**statistics(diff), "left_larger_count": int(np.sum(finite > tolerance)),
                                       "right_larger_count": int(np.sum(finite < -tolerance)),
                                       "equal_within_tolerance_count": int(np.sum(np.abs(finite) <= tolerance)),
                                       "comparison_tolerance": tolerance}
                    group["paired_cap_differences"].append({"left": "capb", "right": right, "beam": beam, "metrics": values})
            group["resource_ratios_capb_over_default"] = []
            for bi, beam in enumerate(BEAMS):
                metrics = {}
                for key in RESOURCES:
                    numerator = data[key][selected,caps.index("capb"),bi]
                    denominator = data[key][selected,caps.index("default"),bi]
                    valid = np.isfinite(numerator) & np.isfinite(denominator) & (denominator > 0) & (numerator >= 0)
                    ratios = np.full(len(numerator), np.nan)
                    np.divide(numerator, denominator, out=ratios, where=valid)
                    metrics[key] = {**statistics(ratios), "total_sequences": int(selected.sum()),
                                    "invalid_or_zero_denominator": int((~valid).sum())}
                group["resource_ratios_capb_over_default"].append({"beam": beam, "metrics": metrics})
        summary["bins"].append(group)
    def example(index, ci, bi):
        return {"input": names[index], "cap": caps[ci], "beam": BEAMS[bi], **metadata[index],
                **{k: (float(data[k][index,ci,bi]) if np.isfinite(data[k][index,ci,bi]) else None) for k in CORE}}
    for ci, cap in enumerate(caps):
        for key in (CORE[0], CORE[4], CORE[6]):
            for bi, beam in enumerate(BEAMS):
                values = data[key][:,ci,bi]
                order = np.argsort(np.where(np.isfinite(values), values, -np.inf))[::-1][:10]
                summary["worst_examples"][f"{cap}_b{beam}_{key}"] = [example(i,ci,bi) for i in order]
        for key, direction, tolerance in [(CORE[0],1,1e-8),(CORE[1],-1,1e-8),(CORE[4],1,1e-9),(CORE[6],1,1e-9)]:
            _, increases, valid = monotonicity(data[key][:,ci,:], direction, tolerance)
            amount = np.max(increases, axis=1)
            eligible = np.flatnonzero(valid & (amount > tolerance))
            order = eligible[np.argsort(amount[eligible])[::-1][:10]]
            summary["nonmonotone_examples"][cap+"_"+key] = [
                {**example(i,ci,int(np.argmax(increases[i]))+1), "values_b20_50_100_200": data[key][i,ci,:].tolist(),
                 "largest_violation": float(amount[i]), "transition": [BEAMS[int(np.argmax(increases[i]))], BEAMS[int(np.argmax(increases[i]))+1]]}
                for i in order]
    for name in ("bpRNA_RFAM_5936", "bpRNA_RFAM_2449", "bpRNA_RFAM_2886", "bpRNA_RFAM_27342"):
        if name not in names:
            continue
        i = names.index(name)
        summary["diagnostic_inputs"][name] = [example(i,ci,bi) for ci in range(len(caps)) for bi in range(4)]
    if "capb" in caps:
        summary["resource_ratio_caveat"] = "Reference values only: concurrent accuracy runs, one measurement, different allocations/dates; elapsed includes GNU time wrapper; GNU wall rounds to 0.01s and zero denominators are excluded. b20 reuses the identical default observation."
    summary["maximum_consistency_errors"] = {k: float(np.nanmax(data[k])) for k in CHECKS}
    summary["rt_range"] = [float(np.nanmin(data["rt_kcal_mol"])), float(np.nanmax(data["rt_kcal_mol"]))]
    summary["eligible_zero_window_sequences"] = {str(w): int(np.sum(data[f"access_w{w}_eligible_windows"][:,0,0] == 0)) for w in (1,4,8,16)}
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    if make_plots:
        plot_distributions(output, data, length, caps)
    return summary


def plot_distributions(output, data, length, caps=CAPS):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    intervals = [(50,99),(100,199),(200,299),(300,399)]
    colors = ("#416a9c", "#299768", "#cd7b35", "#a1559b")
    plotted_caps = [c for c in caps if c != "adopted"]
    for metric, filename in [(CORE[0], "delta_F_distribution.png"), (CORE[1], "retention_distribution.png")]:
        fig, axes = plt.subplots(1,4,figsize=(14,4.2),sharey=True)
        for ax, (lo,hi) in zip(axes, intervals):
            selected = (length >= lo) & (length <= hi)
            for ci, cap in enumerate(plotted_caps):
                source_ci = list(caps).index(cap)
                samples = [data[metric][selected,source_ci,bi] for bi in range(4)]
                samples = [x[np.isfinite(x)] for x in samples]
                if metric == CORE[0]:
                    # Only display clips roundoff; signed values remain unchanged in TSV/JSON.
                    samples = [np.maximum(x,0) for x in samples]
                positions = np.arange(4)+1+(ci-(len(plotted_caps)-1)/2)*.24
                bp = ax.boxplot(samples,positions=positions,widths=.2,whis=(10,90),showfliers=False,
                                patch_artist=True,manage_ticks=False)
                for box in bp["boxes"]: box.set(facecolor=colors[ci],alpha=.65)
                for bi, sample in enumerate(samples):
                    if len(sample):
                        ax.scatter(positions[bi],np.quantile(sample,.99),marker="^",s=23,color=colors[ci],zorder=3)
                        extremum = max(sample) if metric == CORE[0] else min(sample)
                        ax.scatter(positions[bi],extremum,marker="x",s=25,color=colors[ci],zorder=3)
            ax.set(title=f"{lo}–{hi} nt (N={selected.sum():,})",xticks=range(1,5),xticklabels=BEAMS,xlabel="Beam width b")
            ax.grid(axis="y",alpha=.25)
        if metric == CORE[0]:
            axes[0].set_yscale("symlog",linthresh=.001)
            axes[0].set_ylabel("ΔF (kcal/mol; symlog)")
        else:
            axes[0].set(ylim=(-.03,1.04),ylabel="Partition-function ratio Zb / Zexact")
        from matplotlib.patches import Patch
        fig.legend(handles=[Patch(facecolor=colors[ci],alpha=.65,label=cap) for ci,cap in enumerate(plotted_caps)],
                   loc="upper center",ncol=len(plotted_caps),frameon=False)
        fig.text(.5,.01,"Box: 25–75%, line: median, whiskers: 10–90%, triangle: 99%, x: extreme",ha="center",fontsize=9)
        fig.tight_layout(rect=(0,.055,1,.91))
        fig.savefig(output/filename,dpi=170)
        plt.close(fig)


def analyze_packed_array(run, rerun_paths, output):
    provenance = write_metrics(run, rerun_paths, output)
    summary = summarize(output)
    print(json.dumps({"inputs": summary["inputs"], "exact_replacements": len(provenance["exact_replacements"]),
                      "original_job_status_counts": provenance["original_job_status_counts"],
                      "analysis_status_counts": provenance["analysis_status_counts"],
                      "consistency_check_counts": provenance["consistency_check_counts"]}), flush=True)
    return int(set(provenance["analysis_status_counts"]) != {"success"} or
               set(provenance["consistency_check_counts"]) != {"success"})


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--extend", type=Path, help="reuse the validated earlier short-array metrics")
    parser.add_argument("--capb-run", type=Path)
    args = parser.parse_args()
    if bool(args.extend) != bool(args.capb_run):
        parser.error("--extend and --capb-run must be supplied together")
    if args.extend:
        from bprna_capb import extend_capb
        extend_capb(args.extend, args.capb_run, args.output)
    summarize(args.output)
