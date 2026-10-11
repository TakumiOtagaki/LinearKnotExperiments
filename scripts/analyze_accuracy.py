#!/usr/bin/env python3
"""Compare sparse BPP, seven contexts and joint accessibility within a CLI run."""

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re

from rna_records import paired_positions, read_record, scaffold_pairs
from run_cli import digest, git_state, now
from task_storage import RunReader, MemoryFile, input_records, artifact_digest


ROOT = Path(__file__).resolve().parents[1]
TOLERANCE = 1e-9
CONTEXTS = ("Stem", "Hairpin", "Bulge", "Internal", "Multiloop", "Exterior", "PseudoLoop")
PROGRAMS = {"LinearCParty", "LinearPKCapR", "LinearPKRaccess"}


def probability(raw):
    value = float(raw)
    if not math.isfinite(value) or not -TOLERANCE <= value <= 1 + TOLERANCE:
        raise ValueError("invalid probability")
    return value


def output_lines(path, name):
    lines = (path.read_text() if hasattr(path, "read_text") else Path(path).read_text()).splitlines()
    headers = [line.strip()[1:] for line in lines if line.strip().startswith(">")]
    if headers != [name]:
        raise ValueError("output must contain exactly the matching named record")
    return lines


def unpaired(pairs, length):
    values = [1.0] * length
    for (i, j), value in pairs.items():
        values[i - 1] -= value
        values[j - 1] -= value
    if min(values) < -TOLERANCE:
        raise ValueError("BPP incident mass exceeds one")
    return values


def read_bpp(path, record):
    name, sequence, scaffold = record
    fixed = scaffold_pairs(scaffold)
    occupied = paired_positions(fixed)
    pairs = {}
    for raw in output_lines(path, name):
        if not raw.strip() or raw.lstrip().startswith((">", "#")):
            continue
        fields = raw.split()
        if len(fields) != 3:
            raise ValueError("BPP needs i j probability")
        i, j = map(int, fields[:2])
        value = probability(fields[2])
        if not 1 <= i < j <= len(sequence) or (i, j) in pairs:
            raise ValueError("invalid or duplicate BPP coordinate")
        if (i, j) not in fixed and (i in occupied or j in occupied) and value > TOLERANCE:
            raise ValueError("generated pair touches a fixed scaffold base")
        pairs[i, j] = value
    if any(abs(pairs.get(pair, 0) - 1) > TOLERANCE for pair in fixed):
        raise ValueError("fixed scaffold pair does not have probability one")
    return {"pairs": pairs, "unpaired": unpaired(pairs, len(sequence))}


def read_contexts(path, record):
    name, sequence, scaffold = record
    contexts = {}
    for raw in output_lines(path, name):
        if not raw.strip() or raw.lstrip().startswith((">", "#")):
            continue
        fields = raw.split()
        if fields[0] not in CONTEXTS or fields[0] in contexts or len(fields) != len(sequence) + 1:
            raise ValueError("missing, duplicate or malformed structural context")
        contexts[fields[0]] = [probability(value) for value in fields[1:]]
    if set(contexts) != set(CONTEXTS):
        raise ValueError("all seven contexts are required")
    for i, symbol in enumerate(scaffold):
        if abs(sum(contexts[c][i] for c in CONTEXTS) - 1) > TOLERANCE:
            raise ValueError("context probabilities do not sum to one")
        if symbol != "." and abs(contexts["Stem"][i] - 1) > TOLERANCE:
            raise ValueError("scaffold base is not always Stem")
    return contexts


def read_accessibility(path, record):
    name, sequence, scaffold = record
    metadata, windows, positions = {}, {}, set()
    for raw in output_lines(path, name):
        line = raw.strip()
        if line.startswith("#"):
            if "=" in line:
                key, value = line[1:].strip().split("=", 1)
                metadata[key] = value
        elif line and not line.startswith(">"):
            fields = line.split("\t")
            start = int(fields[0]) + 1  # Raccess output is zero-based.
            if not 1 <= start <= len(sequence) or start in positions:
                raise ValueError("invalid or duplicate accessibility start")
            positions.add(start)
            for token in "".join(fields[1:]).split(";"):
                if not token:
                    continue
                width, value = token.split(",")
                width = int(width)
                if width < 1 or start + width - 1 > len(sequence) or (start, width) in windows:
                    raise ValueError("invalid or duplicate accessibility coordinate")
                windows[start, width] = probability(value)
    if metadata.get("output_mode") != "probability" or metadata.get("scaffold") != scaffold:
        raise ValueError("probability output and matching scaffold metadata are required")
    widths = [int(w) for w in metadata["access_lens"].split(",")]
    if len(widths) != len(set(widths)) or any(w < 1 for w in widths):
        raise ValueError("invalid requested window lengths")
    expected = {(s, w) for w in widths for s in range(1, len(sequence) - w + 2)}
    if positions != set(range(1, len(sequence) + 1)) or set(windows) != expected:
        raise ValueError("accessibility coordinates are incomplete; do not fill with zero")
    rt, free_energy = float(metadata["folding_rt"]), float(metadata["ensemble_free_energy"])
    if not math.isfinite(rt) or rt <= 0 or not math.isfinite(free_energy):
        raise ValueError("finite high-precision thermodynamic metadata are required")
    return {"windows": windows, "widths": widths, "rt": rt, "free_energy": free_energy}


def mean(values):
    return sum(values) / len(values) if values else None


def compare_bpp(candidate, reference, record):
    free = {i for i, symbol in enumerate(record[2], 1) if symbol == "."}
    keys = set(candidate["pairs"]) | set(reference["pairs"])
    error = sum(abs(candidate["pairs"].get(pair, 0) - reference["pairs"].get(pair, 0))
                for pair in keys if set(pair) <= free)
    q_error = sum(abs(candidate["unpaired"][i - 1] - reference["unpaired"][i - 1]) for i in free)
    return {"free_bases": len(free), "generated_pair_l1_sum": error,
            "base_state_tv_mean_Gfree": (error + q_error / 2) / len(free) if free else None,
            "unpaired_mae_Gfree": q_error / len(free) if free else None,
            "sparse_coordinate_union_count": len(keys)}


def compare_contexts(candidate, reference, record):
    free = [i for i, symbol in enumerate(record[2]) if symbol == "."]
    return {"free_bases": len(free),
            "context_tv_mean_Gfree": mean([
                sum(abs(candidate[c][i] - reference[c][i]) for c in CONTEXTS) / 2 for i in free]),
            "context_mae_Gfree": {c: mean([abs(candidate[c][i] - reference[c][i]) for i in free])
                                  for c in CONTEXTS}}


def compare_accessibility(candidate, reference, record):
    if set(candidate["windows"]) != set(reference["windows"]):
        raise ValueError("window sets differ")
    if abs(candidate["rt"] - reference["rt"]) > 1e-12:
        raise ValueError("RT differs between candidate and reference")
    rt = reference["rt"]
    delta = candidate["free_energy"] - reference["free_energy"]
    paired_prefix = [0]
    for symbol in record[2]:
        paired_prefix.append(paired_prefix[-1] + (symbol != "."))
    per_width = {}
    for width in reference["widths"]:
        keys = [(s, w) for s, w in reference["windows"] if w == width]
        eligible = [(s, w) for s, w in keys if
                    paired_prefix[s + w - 1] == paired_prefix[s - 1]]
        errors, opening_errors = [], []
        lost = gained = both_zero = both_positive = 0
        for key in eligible:
            p, q = candidate["windows"][key], reference["windows"][key]
            errors.append(abs(p - q))
            if p > 0 and q > 0:
                opening_errors.append(rt * abs(math.log(p) - math.log(q)))
                both_positive += 1
            elif q > 0:
                lost += 1
            elif p > 0:
                gained += 1
            else:
                both_zero += 1
        per_width[str(width)] = {"all_windows": len(keys), "eligible_windows": len(eligible),
                                 "probability_mae_eligible": mean(errors),
                                 "probability_max_error_eligible": max(errors) if errors else None,
                                 "opening_mae_both_positive_kcal_mol": mean(opening_errors),
                                 "both_positive": both_positive, "both_zero": both_zero,
                                 "reference_positive_candidate_zero": lost,
                                 "reference_zero_candidate_positive": gained}
    return {"delta_free_energy_kcal_mol": delta,
            "delta_free_energy_per_nt": delta / len(record[1]),
            "delta_log_z": -delta / rt, "rt_kcal_mol": rt,
            "free_energy_source": "Raccess full-precision metadata", "accessibility": per_width}


def read_job(run, job, record):
    reader = run if isinstance(run, RunReader) else RunReader(run)
    if job["program"] == "LinearCParty":
        command = job["command"]
        if "--bpp-cutoff" in command and float(command[command.index("--bpp-cutoff") + 1]) != 0:
            raise ValueError("accuracy requires unthresholded BPP output")
        if "--bpp" not in command:
            raise ValueError("BPP export is required")
        return read_bpp(MemoryFile(reader.artifact(job, "bpp.tsv")), record)
    if job["program"] == "LinearPKCapR":
        return read_contexts(MemoryFile(reader.artifact(job, "stdout.txt")), record)
    return read_accessibility(MemoryFile(reader.artifact(job, "stdout.txt")), record)


def check_consistency(input_name, mode, repetition, group):
    q = group["LinearCParty"]["unpaired"]
    stem = group["LinearPKCapR"]["Stem"]
    windows = group["LinearPKRaccess"]["windows"]
    stem_error = max(abs(1 - value - stem[i]) for i, value in enumerate(q))
    single_errors = [abs(value - windows[i + 1, 1]) for i, value in enumerate(q)
                     if (i + 1, 1) in windows]
    singleton_error = max(single_errors) if single_errors else None
    bound_violation = max([0.0] + [p - min(q[s - 1:s + w - 1])
                                   for (s, w), p in windows.items()])
    return {"input": input_name, "mode": mode, "repetition": repetition,
            "stem_vs_bpp_max_error": stem_error,
            "accessibility_w1_vs_unpaired_max_error": singleton_error,
            "accessibility_above_marginal_max_violation": bound_violation,
            "status": "success" if max(stem_error, singleton_error or 0,
                                         bound_violation) <= TOLERANCE else "failed"}


def analyze(run, reference_mode, extra_runs=()):
    reader = RunReader(run)
    analysis_provenance = {"generated_utc": now(), "experiment_git_state": git_state(ROOT),
                           "analysis_script_sha256": digest(Path(__file__)),
                           "record_helpers_sha256": digest(Path(__file__).with_name("rna_records.py")),
                           "run_manifest_sha256": reader.manifest_sha256}
    manifest = reader.manifest
    if manifest["status"] == "running" or manifest.get("summary_only"):
        raise ValueError("wait for the run to finish")
    manifest = {**manifest, "inputs": manifest["inputs"].copy(), "jobs": manifest["jobs"].copy()}
    analysis_provenance["additional_runs"] = []
    records = input_records(reader)
    input_hashes = {i["name"]: i["sha256"] for i in manifest["inputs"]}
    job_ids = {j["id"] for j in manifest["jobs"]}
    for extra in extra_runs:
        extra_reader = RunReader(extra)
        saved = extra_reader.manifest
        if saved["status"] == "running" or saved.get("summary_only"):
            raise ValueError("finish and aggregate every sensitivity run")
        if (saved["software_before"]["commit"] != manifest["software_before"]["commit"] or
                saved["config"]["model_settings"] != manifest["config"]["model_settings"] or
                saved["software_lock"]["energy_parameter_sha256"] != manifest["software_lock"]["energy_parameter_sha256"]):
            raise ValueError("sensitivity run models differ")
        extra_records = input_records(extra_reader)
        for item in saved["inputs"]:
            if item["name"] in input_hashes:
                if item["sha256"] != input_hashes[item["name"]]:
                    raise ValueError("sensitivity runs require identical S and G")
            else:
                input_hashes[item["name"]] = item["sha256"]
                manifest["inputs"].append(item)
                records[item["name"]] = extra_records[item["name"]]
        for job in saved["jobs"]:
            if job["id"] in job_ids:
                raise ValueError("duplicate observations in sensitivity runs")
            job_ids.add(job["id"])
            manifest["jobs"].append({**job, "artifact_run": str(extra)})
        analysis_provenance["additional_runs"].append({"path": str(extra), "run_manifest_sha256": extra_reader.manifest_sha256})
    all_jobs = [job for job in manifest["jobs"] if not job["warmup"]]
    jobs = [job for job in all_jobs if job["program"] in PROGRAMS]
    parsed, errors = {}, {}
    for job in jobs:
        if job["status"] != "success":
            continue
        try:
            parsed[job["id"]] = read_job(reader, job, records[job["input"]])
        except (ValueError, KeyError, OSError) as exc:
            errors[job["id"]] = str(exc)
    references = {(job["input"], job["program"], job["mode"]): job["id"] for job in jobs
                  if job["repetition"] == 0}
    largest_modes = {}
    if extra_runs:
        for name, program, mode in references:
            match = re.fullmatch(r"beam(\d+)_(default|nocap)", mode)
            if match:
                key = (name, program, match[2])
                if key not in largest_modes or int(match[1]) > largest_modes[key][0]:
                    largest_modes[key] = (int(match[1]), mode)
    policy = manifest.get("config", {}).get("analysis", {})
    compare = {"LinearCParty": compare_bpp, "LinearPKCapR": compare_contexts,
               "LinearPKRaccess": compare_accessibility}
    observations = []
    for job in jobs:
        row = {key: job[key] for key in ("id", "input", "mode", "program", "repetition", "status")}
        row["measurement_scope"] = "whole CLI process, including requested output"
        artifacts = {"stdout.txt": artifact_digest(reader, job, "stdout.txt") if "directory" in job else None}
        command = job.get("command", [])
        if "--bpp" in command:
            artifacts["bpp_export"] = artifact_digest(reader, job, "bpp.tsv")
        row["observed_artifact_sha256"] = artifacts
        row.update({key: job.get(key) for key in ("time_wall_seconds", "peak_rss_kib")})
        selected_mode, kind = reference_mode, "exact_error"
        if len(records[job["input"]][1]) >= policy.get("exact_max_length_exclusive", float("inf")):
            group = "nocap" if job["mode"].endswith("_nocap") else "default"
            selected_mode = policy["long_reference_modes"][group]
            if extra_runs:
                selected_mode = largest_modes.get((job["input"], job["program"], group), (0,selected_mode))[1]
            kind = "beam_width_sensitivity"
        reference = references.get((job["input"], job["program"], selected_mode))
        row.update(reference_mode=selected_mode, comparison_kind=kind)
        row["reference_job"] = reference
        if job["status"] != "success":
            row["analysis_status"] = "job_failed"
        elif job["id"] in errors:
            row.update(analysis_status="artifact_error", error=errors[job["id"]])
        elif reference not in parsed:
            row["analysis_status"] = "reference_unavailable"
        else:
            try:
                row["metrics"] = compare[job["program"]](
                    parsed[job["id"]], parsed[reference], records[job["input"]])
                row["analysis_status"] = "success"
            except ValueError as exc:
                row.update(analysis_status="comparison_error", error=str(exc))
        observations.append(row)
    by_group = {}
    for job in jobs:
        if job["id"] in parsed:
            by_group.setdefault((job["input"], job["mode"], job["repetition"]), {})[
                job["program"]] = parsed[job["id"]]
    checks = []
    for (input_name, mode, repetition), group in by_group.items():
        if not PROGRAMS <= set(group):
            continue
        checks.append(check_consistency(input_name, mode, repetition, group))
    return {"schema_version": 1, "purpose": manifest["purpose"],
            "analysis_provenance": analysis_provenance,
            "reference_mode": reference_mode, "reference_repetition": 0,
            "run_status": manifest["status"], "measured_jobs": len(all_jobs),
            "job_status_counts": dict(Counter(job["status"] for job in all_jobs)),
            "software_commit": manifest["software_before"]["commit"],
            "experiment_commit": manifest["experiment_before"]["commit"],
            "config_sha256": manifest["config_sha256"],
            "software_changed_during_run": manifest["software_changed_during_run"],
            "binaries_changed_during_run": manifest["binaries_changed_during_run"],
            "observations": observations, "consistency_checks": checks,
            "not_analyzed": [{"id": j["id"], "status": j["status"], "program": j["program"]}
                             for j in all_jobs if j["program"] not in PROGRAMS]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("--reference-mode", default="exact")
    parser.add_argument("--compare-runs", nargs="*", type=Path, default=[],
                        help="merge staged long-width runs; use each RNA's largest planned width per cap regime")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run = args.run if args.run.is_absolute() else ROOT / args.run
    output = args.out if args.out.is_absolute() else ROOT / args.out
    if output.exists():
        parser.error("use a new analysis output path")
    extra = [p if p.is_absolute() else ROOT / p for p in args.compare_runs]
    result = analyze(run, args.reference_mode, extra)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    failures = sum(row["analysis_status"] != "success" for row in result["observations"])
    failures += sum(check["status"] != "success" for check in result["consistency_checks"])
    print(json.dumps({"measured_jobs": result["measured_jobs"],
                      "analyzed_observations": len(result["observations"]),
                      "consistency_checks": len(result["consistency_checks"]),
                      "analysis_failures": failures}))
    return int(failures > 0)


if __name__ == "__main__":
    raise SystemExit(main())
