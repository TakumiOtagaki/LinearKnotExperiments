#!/usr/bin/env python3
"""Stratify adopted E1 exact errors by crossing-pair/PseudoLoop mass.

Read the entire short cohort through the existing validated tar readers.
Recompute and cross-check the published RNA metrics; retain no extracted files.
One node-local SQLite database joins the differently packed exact/beam tasks.
"""

import argparse
from array import array
from collections import Counter, defaultdict
import csv
import gc
import json
import math
from pathlib import Path
import pickle
import shutil
import sqlite3
import tempfile
import time

import numpy as np

import analyze_accuracy as accuracy
import analyze_e1_production as production
from analyze_bprna_short import CHECKS, crossing_descriptors
from rna_records import scaffold_pairs
from run_cli import digest, now
from task_storage import local_workspace


BEAMS = (20, 50, 100, 200)
AXES = ("exact_crossing_pair_sum", "exact_pseudoloop_sum")
STRATA = ("lt0.1", "0.1_to_lt1", "ge1")
GROUPS = [("all", "all")] + [(a, s) for a in AXES for s in STRATA]
BASE_METRICS = ("delta_free_energy_kcal_mol", "abs_delta_free_energy_kcal_mol",
                "base_state_tv_mean_Gfree", "context_tv_mean_Gfree",
                "access_w16_probability_mae_eligible")
LOCAL = ("pseudoloop_ge0.1", "crossing_bpp_ge0.1", "access_w4_ge0.1")
RNA_METRICS = BASE_METRICS + tuple(f"{m}_{s}" for m in LOCAL for s in ("mae", "max"))
DEPENDENCIES = (Path(__file__).name, "analyze_e1_production.py", "analyze_accuracy.py",
                "analyze_bprna_short.py", "analyze_bprna_long.py", "bprna_capb.py",
                "task_storage.py", "rna_records.py", "run_cli.py")


def stratum(value):
    return STRATA[0] if value < 0.1 else STRATA[1] if value < 1 else STRATA[2]


def crossing_pairs(record, bpp):
    """Endpoints in different enclosing-G masks cross at least one G pair.

    Scaffold endpoints cannot belong to an extension pair. Count each pair
    once, even if it crosses several scaffold pairs.
    """
    masks, mask, stack = [], 0, []
    for i, symbol in enumerate(record[2]):
        if symbol == "(":
            stack.append(1 << i)
            mask |= stack[-1]
        elif symbol == ")":
            mask ^= stack.pop()
        masks.append(mask)
    return {pair: p for pair, p in bpp["pairs"].items()
            if all(record[2][i - 1] == "." for i in pair)
            and masks[pair[0] - 1] != masks[pair[1] - 1]}


def select_reference(record, parsed):
    pairs = crossing_pairs(record, parsed["LinearCParty"])
    contexts = parsed["LinearPKCapR"]
    descriptors = dict(zip(AXES, (sum(pairs.values()), sum(contexts["PseudoLoop"]))))
    previous = crossing_descriptors(record, parsed["LinearCParty"], contexts)
    if descriptors[AXES[0]] != previous["exact_crossing_pair_probability_sum"]:
        raise ValueError("crossing descriptor differs from existing definition")
    # Verify the geometric meaning independently for every selected reference pair.
    fixed = scaffold_pairs(record[2])
    selected_pairs = {pair: p for pair, p in pairs.items() if p >= 0.1}
    for i, j in selected_pairs:
        if not any(i < k < j < l or k < i < l < j for k, l in fixed):
            raise ValueError("selected crossing pair does not cross G")
    selected = {
        LOCAL[0]: {i: p for i, p in enumerate(contexts["PseudoLoop"]) if p >= 0.1},
        LOCAL[1]: selected_pairs,
        LOCAL[2]: {key: p for key, p in parsed["LinearPKRaccess"]["windows"].items()
                   if key[1] == 4 and p >= 0.1},
    }
    for start, width in selected[LOCAL[2]]:
        if any(c != "." for c in record[2][start - 1:start + width - 1]):
            raise ValueError("positive accessibility window contains a scaffold-paired base")
    return descriptors, selected


def local_errors(parsed, selected):
    values = (
        {i: abs(parsed["LinearPKCapR"]["PseudoLoop"][i] - q)
         for i, q in selected[LOCAL[0]].items()},
        {pair: abs(parsed["LinearCParty"]["pairs"].get(pair, 0.0) - q)
         for pair, q in selected[LOCAL[1]].items()},
        {key: abs(parsed["LinearPKRaccess"]["windows"][key] - q)
         for key, q in selected[LOCAL[2]].items()},
    )
    return dict(zip(LOCAL, values))


def quantiles(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return dict(n=0, p50=None, p90=None, max=None, min=None)
    p50, p90 = np.quantile(values, [0.5, 0.9], method="linear")
    return dict(n=len(values), p50=float(p50), p90=float(p90),
                max=float(values.max()), min=float(values.min()))


def write_tsv(path, fields, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: "NA" if v is None else v for k, v in row.items()})


def read_previous(path, plans):
    provenance = json.loads((path / "provenance.json").read_text())
    validation = json.loads((path / "validation.json").read_text())
    checksum = digest(path / "metrics.tsv")
    if (not validation["all_passed"] or checksum != provenance["metrics_sha256"]
            or checksum != validation["metrics_sha256"]):
        raise ValueError("previous E1 metrics are not unchanged and validated")
    for run, _ in plans:
        prior = next(r for r in provenance["runs"] if Path(r["path"]).name == run.name)
        if prior["array_plan_sha256"] != digest(run / "array_plan.json"):
            raise ValueError("previous array plan changed")
    rows = {}
    with (path / "metrics.tsv").open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            key = row["input"], int(row["beam"])
            if key in rows or row["analysis_status"] != "success" or row["reference_mode"] != "exact":
                raise ValueError("duplicate or invalid previous E1 row")
            values = {k: None if row[k] == "NA" else float(row[k]) for k in BASE_METRICS}
            rows[key] = row["input_sha256"], row["tier"], values
    return rows, {"directory": str(path), "metrics_sha256": checksum,
                  "validation_sha256": digest(path / "validation.json"),
                  "provenance_sha256": digest(path / "provenance.json")}


def check_previous(old, name, beam, input_hash, tier, row, maxima):
    saved_hash, saved_tier, saved = old.pop((name, beam))
    if saved_hash != input_hash or saved_tier != tier:
        raise ValueError("previous input hash or tier differs")
    for key in BASE_METRICS:
        if row[key] is None or saved[key] is None:
            if row[key] != saved[key]:
                raise ValueError("previous empty metric domain differs")
        else:
            error = abs(row[key] - saved[key])
            maxima[key] = max(maxima[key], error)
            if not math.isclose(row[key], saved[key], rel_tol=1e-12, abs_tol=1e-12):
                raise ValueError("recomputed RNA metric differs: " + key)


def analyze(args):
    started = time.monotonic()
    beam_plan, exact_plan = [production.load_plan(p) for p in (args.beam_run, args.exact_run)]
    production.check_reference_plans(beam_plan, exact_plan)
    production.check_exact_modes(beam_plan, exact_plan)
    beam_tasks, hashes = production.cohort(beam_plan, production.SHORT)
    exact_tasks, exact_hashes = production.cohort(exact_plan, production.SHORT)
    modes = [f"beam{b}_outer" for b in BEAMS]
    if hashes != exact_hashes:
        raise ValueError("exact cohort or reference mode mismatch")
    if {m["name"] for m in beam_plan["config"]["modes"]} != set(modes + ["exact"]):
        raise ValueError("unexpected beam conditions")
    plans = [(args.beam_run, beam_plan), (args.exact_run, exact_plan)]
    old, prior = read_previous(args.previous, plans)
    if len(old) != len(hashes) * len(BEAMS):
        raise ValueError("previous cohort is incomplete")
    if args.out.exists():
        raise ValueError("output already exists; do not overwrite a completed analysis")
    provenance = {
        "generated_utc": now(), "schema_version": 1, "sampling": "none; full short E1 cohort",
        "previous_analysis": prior,
        "runs": [{"name": p.name, "array_plan_sha256": digest(p / "array_plan.json"),
                  "inputs_archive_sha256": plan["inputs_archive_sha256"],
                  "software_commit": plan["software_before"]["commit"]} for p, plan in plans],
        "script_sha256": {n: digest(Path(__file__).with_name(n)) for n in DEPENDENCIES},
        "tasks": [], "selected_job_status_counts": Counter(), "excluded_job_status_counts": Counter(),
        "consistency_checks": Counter(), "consistency_maxima": dict.fromkeys(CHECKS, 0.0),
        "previous_metric_max_absolute_differences": dict.fromkeys(BASE_METRICS, 0.0),
        "quantile_method": "numpy linear interpolation; exact exported probabilities; no error clipping",
        "strata": {a: "<0.1; 0.1 <= sum <1; >=1" for a in AXES},
        "local_selection": "exact probability >=0.1, inclusive; no selection by beam output",
        "rna_weighting": "equal RNA weight; means/maxima computed inside each selected RNA first",
        "local_weighting": "pooled bases, pairs, or windows; overlapping windows are descriptive, not independent replicates",
        "missing": "empty eligible domains are NA and omitted from quantiles; zeros in nonempty domains retained",
        "storage": "tar member streams; no extraction; temporary SQLite only in node-local scratch",
        "scope": "same S/G/model; exact uncapped, beam b=20/50/100/200 with h=d=max(b,30); total accessibility only",
    }
    counts = {g: Counter() for g in GROUPS}
    rna_values, pooled_values, witnesses = defaultdict(lambda: array("d")), defaultdict(lambda: array("d")), {}
    seen = set()
    with local_workspace(args.scratch_dir) as scratch:
        db = sqlite3.connect(str(scratch / "exact.sqlite"))
        try:
            db.execute("CREATE TABLE refs (name TEXT PRIMARY KEY, hash TEXT, record BLOB, parsed BLOB, descriptors BLOB, selected BLOB)")
            for number, task in enumerate(exact_tasks, 1):
                reader, task_hashes, jobs, records = production.open_task(args.exact_run, exact_plan, task, ["exact"], provenance)
                for name, record in records.items():
                    parsed, _ = production.read_group(reader, jobs, name, "exact", record, provenance)
                    descriptors, selected = select_reference(record, parsed)
                    groups = [("all", "all")] + [(a, stratum(descriptors[a])) for a in AXES]
                    for group in groups:
                        counts[group].update(rna=1, **{task["tier"]: 1})
                        for key in LOCAL:
                            counts[group][key + "_rna"] += bool(selected[key])
                            counts[group][key + "_events"] += len(selected[key])
                        counts[group]["crossing_pair_sum"] += descriptors[AXES[0]]
                        counts[group]["pseudoloop_sum"] += descriptors[AXES[1]]
                    db.execute("INSERT INTO refs VALUES (?,?,?,?,?,?)", (name, task_hashes[name],
                        *(pickle.dumps(v, protocol=pickle.HIGHEST_PROTOCOL) for v in (record, parsed, descriptors, selected))))
                db.commit()
                if number % 5 == 0 or number == len(exact_tasks):
                    print(f"exact {number}/{len(exact_tasks)} tasks; {time.monotonic()-started:.1f}s", flush=True)
                del reader, jobs, records, parsed
                gc.collect()
            for number, task in enumerate(beam_tasks, 1):
                reader, task_hashes, jobs, records = production.open_task(args.beam_run, beam_plan, task, modes, provenance)
                for name, record in records.items():
                    if name in seen:
                        raise ValueError("duplicate beam RNA")
                    seen.add(name)
                    saved = db.execute("SELECT hash,record,parsed,descriptors,selected FROM refs WHERE name=?", (name,)).fetchone()
                    if saved is None or saved[0] != task_hashes[name] or pickle.loads(saved[1]) != record:
                        raise ValueError("beam/exact sequence or scaffold differs")
                    reference, descriptors, selected = (pickle.loads(v) for v in saved[2:])
                    groups = [("all", "all")] + [(a, stratum(descriptors[a])) for a in AXES]
                    for beam, mode in zip(BEAMS, modes):
                        parsed, _ = production.read_group(reader, jobs, name, mode, record, provenance)
                        row = production.compare(parsed, reference, record)
                        check_previous(old, name, beam, task_hashes[name], task["tier"], row,
                                       provenance["previous_metric_max_absolute_differences"])
                        local = local_errors(parsed, selected)
                        for metric, errors in local.items():
                            values = list(errors.values())
                            row[metric + "_mae"] = accuracy.mean(values)
                            row[metric + "_max"] = max(values) if values else None
                        for axis, layer in groups:
                            for metric in RNA_METRICS:
                                rna_values[axis, layer, beam, metric].append(
                                    float("nan") if row[metric] is None else row[metric])
                            for metric, errors in local.items():
                                key = axis, layer, beam, metric
                                pooled_values[key].extend(errors.values())
                                if errors:
                                    coordinate, maximum = max(errors.items(), key=lambda item: item[1])
                                    if key not in witnesses or maximum > witnesses[key]["absolute_error"]:
                                        # Human-readable coordinates are always 1-based, inclusive.
                                        if metric == LOCAL[0]:
                                            location = str(coordinate + 1)
                                        elif metric == LOCAL[1]:
                                            location = f"{coordinate[0]}-{coordinate[1]}"
                                        else:
                                            location = f"{coordinate[0]}-{coordinate[0]+coordinate[1]-1}"
                                        witnesses[key] = dict(axis=axis, stratum=layer, beam=beam, metric=metric,
                                            input=name, coordinate_1based=location, absolute_error=maximum,
                                            exact_probability=selected[metric][coordinate],
                                            length=len(record[1]), **descriptors)
                        del parsed
                    del reference
                if number % 5 == 0 or number == len(beam_tasks):
                    print(f"beam {number}/{len(beam_tasks)} tasks; {len(seen)} RNA; {time.monotonic()-started:.1f}s", flush=True)
                del reader, jobs, records
                gc.collect()
        finally:
            db.close()
        if seen != set(hashes) or old:
            raise ValueError("incomplete new/previous cohort join")
        if provenance["consistency_checks"] != {"success": len(hashes) * 5}:
            raise ValueError("missing CLI consistency checks")
        for axis in AXES:
            if sum(counts[axis, s]["rna"] for s in STRATA) != len(hashes):
                raise ValueError("strata do not partition the RNA cohort")
        output = scratch / "summary"
        output.mkdir()
        count_fields = ["axis", "stratum", "rna"] + [t[0] for t in production.SHORT] + [
            "crossing_pair_sum", "pseudoloop_sum"] + [m + s for m in LOCAL for s in ("_rna", "_events")]
        write_tsv(output / "stratum_counts.tsv", count_fields,
                  [dict(axis=a, stratum=s, **{k: counts[a, s][k] for k in count_fields[2:]}) for a, s in GROUPS])
        rna_rows, local_rows = [], []
        for axis, layer in GROUPS:
            for beam in BEAMS:
                for metric in RNA_METRICS:
                    q = quantiles(rna_values[axis, layer, beam, metric])
                    if len(rna_values[axis, layer, beam, metric]) != counts[axis, layer]["rna"]:
                        raise ValueError("missing RNA-level observation")
                    rna_rows.append(dict(axis=axis, stratum=layer, beam=beam, metric=metric,
                        rna_total=counts[axis, layer]["rna"], missing_rna=counts[axis, layer]["rna"]-q["n"], **q))
                for metric in LOCAL:
                    q = quantiles(pooled_values[axis, layer, beam, metric])
                    if q["n"] != counts[axis, layer][metric + "_events"]:
                        raise ValueError("local event counts differ between beams")
                    local_rows.append(dict(axis=axis, stratum=layer, beam=beam, metric=metric,
                        rna_total=counts[axis, layer]["rna"], rna_selected=counts[axis, layer][metric+"_rna"], **q))
        write_tsv(output / "rna_error_quantiles.tsv", ["axis", "stratum", "beam", "metric", "rna_total", "missing_rna", "n", "p50", "p90", "max", "min"], rna_rows)
        write_tsv(output / "local_error_quantiles.tsv", ["axis", "stratum", "beam", "metric", "rna_total", "rna_selected", "n", "p50", "p90", "max", "min"], local_rows)
        write_tsv(output / "max_error_examples.tsv", ["axis", "stratum", "beam", "metric", "input", "coordinate_1based", "absolute_error", "exact_probability", "length"] + list(AXES),
                  [witnesses[key] for key in sorted(witnesses)])
        provenance.update(inputs=len(seen), rna_comparisons=len(seen)*len(BEAMS),
                          elapsed_seconds=time.monotonic()-started, validation="passed",
                          output_sha256={p.name: digest(p) for p in output.glob("*.tsv")})
        production.write_json(output / "provenance.json", provenance)
        # Publish only the small completed summaries, after all audits have passed.
        shutil.copytree(output, args.out)
    print(f"completed: {args.out}; {len(seen)} RNA, no subsampling", flush=True)


def self_test():
    assert [stratum(v) for v in (0, 0.099999, 0.1, 0.999999, 1)] == [STRATA[0], STRATA[0], STRATA[1], STRATA[1], STRATA[2]]
    record = ("toy", "ACGUACGUACGU", "..((....))..")
    pairs = {(1, 5): .4, (5, 12): .1, (1, 12): .8, (5, 8): .2, (3, 10): 1.}
    assert crossing_pairs(record, {"pairs": pairs}) == {(1, 5): .4, (5, 12): .1}
    for pair in pairs:
        i, j = pair
        direct = all(record[2][x-1] == "." for x in pair) and any(
            i < k < j < l or k < i < l < j for k, l in scaffold_pairs(record[2]))
        assert (pair in crossing_pairs(record, {"pairs": pairs})) == direct
    selected = {LOCAL[0]: {0: .1}, LOCAL[1]: {(1, 5): .3}, LOCAL[2]: {(2, 4): .4}}
    parsed = {"LinearPKCapR": {"PseudoLoop": [.25]}, "LinearCParty": {"pairs": {}},
              "LinearPKRaccess": {"windows": {(2, 4): .1}}}
    errors = local_errors(parsed, selected)
    assert math.isclose(errors[LOCAL[0]][0], .15)
    assert errors[LOCAL[1]][1, 5] == .3  # absent sparse BPP is zero
    assert math.isclose(errors[LOCAL[2]][2, 4], .3)
    assert quantiles([0, np.nan, 1])["n"] == 2
    assert quantiles([0, np.nan, 1])["p50"] == .5
    assert quantiles([])["max"] is None  # empty domain is not zero error
    assert accuracy.mean([0, 0, 1]) == 1/3  # RNA means precede RNA quantiles
    print("self-test: crossing geometry, inclusive thresholds, sparse-zero convention, local errors, quantiles and empty domains passed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--beam-run", type=Path)
    parser.add_argument("--exact-run", type=Path)
    parser.add_argument("--previous", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--scratch-dir", type=Path)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    else:
        if any(getattr(args, key) is None for key in ("beam_run", "exact_run", "previous", "out")):
            parser.error("--beam-run, --exact-run, --previous and --out are required")
        analyze(args)


if __name__ == "__main__":
    main()
