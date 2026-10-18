#!/usr/bin/env python3
"""Audit fixed-scaffold reference membership using the pinned exact API kernel.

All builds and temporary files live under TMPDIR. Raw process results are kept
in one JSONL; cohort, build log and observer sources are kept in one tar.gz.
The implementation checkout and build directory are read-only inputs.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
import os
from pathlib import Path
import statistics as st
import subprocess
import tarfile
import tempfile
import time

CANONICAL = {"AU", "UA", "CG", "GC", "GU", "UG"}
ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def command(cmd, **kwargs):
    return subprocess.run(list(map(str, cmd)), capture_output=True, text=True,
                          check=True, **kwargs).stdout.strip()


def parse_structure(text):
    stacks = {x: [] for x in "([{<"}
    out = []
    for j, c in enumerate(text, 1):
        if c in stacks:
            stacks[c].append(j)
        elif c in ")]}>" :
            opening = "([{<"[")]}>".index(c)]
            out.append((stacks[opening].pop(), j, opening))
        else:
            assert c in ". _xX"
    assert not any(stacks.values())
    return out


def render(n, pairs, opening="("):
    text = ["."] * n
    closing = ")]}>"["([{<".index(opening)]
    for i, j in pairs:
        assert text[i - 1] == text[j - 1] == "."
        text[i - 1], text[j - 1] = opening, closing
    return "".join(text)


def crosses(a, b):
    i, j = a
    k, l = b
    return i < k < j < l or k < i < l < j


def prepare(r):
    seq, ref = r["seq"], r["ref"]
    assert len(seq) == len(ref) and set(seq) <= set("ACGU")
    n = len(seq)
    annotated = parse_structure(ref)
    canonical = [(i, j, o) for i, j, o in annotated if seq[i-1] + seq[j-1] in CANONICAL]
    g, pk = set(map(tuple, r["nested"])), list(map(tuple, r["pk"]))
    assert set(pk) == {(i, j) for i, j, o in canonical if o == "["}
    assert g == {(i, j) for i, j, o in canonical if o == "(" and j-i >= 4}
    full = ["."] * n
    filtered = ["."] * n
    for i, j, o in canonical:
        if (i, j) in g:
            continue
        full[i-1], full[j-1] = o, ")]}>"["([{<".index(o)]
        if j-i >= 4:
            filtered[i-1], filtered[j-1] = full[i-1], full[j-1]
    targets = {"full": "".join(full), "filtered": "".join(filtered),
               "reported": render(n, pk, "["), "scaffold": "." * n}
    targets.update({f"single_{k}": render(n, [p], "[") for k, p in enumerate(pk)})
    # Identical targets share one DP call; aliases are retained in each record.
    unique, aliases, seen = {}, {}, {}
    for label, target in targets.items():
        if target not in seen:
            unique[label] = target
            seen[target] = label
        aliases[label] = seen[target]
    return {"record": r, "scaffold": render(n, g), "targets": targets,
            "unique": unique, "aliases": aliases,
            "noncanonical_removed": len(annotated)-len(canonical),
            "short_canonical": [(i, j) for i, j, o in canonical if j-i < 4],
            "other_bracket_pairs": [(i, j) for i, j, o in canonical if o in "{<"]}


def execute(prepared, driver, parameter, timeout, energy_only=False):
    r = prepared["record"]
    cmd = [str(driver), r["seq"], prepared["scaffold"], str(parameter)]
    if energy_only:
        cmd.append("--energy-only")
    target_input = "".join(k + "\t" + v + "\n" for k, v in prepared["unique"].items())
    start = time.monotonic()
    record = {"id": r["id"], "argv": cmd, "stdin": target_input,
              "aliases": prepared["aliases"], "energy_only": energy_only}
    try:
        p = subprocess.run(cmd, input=target_input, text=True, capture_output=True,
                           timeout=timeout,
                           env=dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                                    MKL_NUM_THREADS="1", LC_ALL="C"))
        record.update(returncode=p.returncode, stderr=p.stderr,
                      stdout_sha256=hashlib.sha256(p.stdout.encode()).hexdigest())
        if p.returncode == 0:
            record["result"] = json.loads(p.stdout)
        else:
            record["failed_stdout"] = p.stdout
    except subprocess.TimeoutExpired as e:
        record.update(returncode=None, error="timeout", stderr=str(e.stderr or ""))
        record["failed_stdout"] = (e.stdout or b"").decode(errors="replace") if isinstance(e.stdout, bytes) else e.stdout
    record["elapsed_seconds"] = time.monotonic() - start
    return record


def get_eval(raw, label):
    return raw["result"]["evaluations"].get(raw["aliases"][label], {"reachable": None, "missing": True})


def salvage_partial_record(raw):
    """Read only complete JSON values emitted before a timed-out observer died.

    Keep the failed attempt and its original stdout unchanged in raw.jsonl.
    An incomplete JSON value is never inferred or repaired into a measurement.
    """
    if raw.get("returncode") == 0 or not raw.get("energy_only"):
        return raw
    text = raw.get("failed_stdout") or ""
    decoder = json.JSONDecoder()
    try:
        mfe, _ = decoder.raw_decode(text, text.index('"mfe":')+len('"mfe":'))
        pos = text.index('"evaluations":{')+len('"evaluations":{')
    except (ValueError, json.JSONDecodeError):
        return raw
    evaluations = {}
    while pos < len(text):
        while pos < len(text) and text[pos] in " ,\n\t":
            pos += 1
        try:
            label, pos = decoder.raw_decode(text, pos)
            assert text[pos] == ":"
            value, pos = decoder.raw_decode(text, pos+1)
        except (ValueError, json.JSONDecodeError, AssertionError):
            break
        evaluations[label] = value
    if not mfe.get("reachable") or not all(raw["aliases"][x] in evaluations for x in ("full", "reported", "scaffold")):
        return raw
    if not mfe.get("log_weight"):
        return raw
    return {**raw, "partial_salvaged": True, "result": {
        "mfe": mfe, "evaluations": evaluations, "log_z": None, "probabilities": [],
        "rt_kcal": -mfe["energy"]/mfe["log_weight"],
        "rt_source": "derived from complete emitted MFE energy and log_weight"}}


def diagnostics(prepared, label, raw):
    if get_eval(raw, label).get("missing"):
        return "unmeasured", ["light_run_time_limit"]
    pairs = [(i, j) for i, j, _ in parse_structure(prepared["targets"][label])]
    g = set(map(tuple, prepared["record"]["nested"]))
    seq = prepared["record"]["seq"]
    reasons = []
    if any(j-i < 4 for i, j in pairs):
        reasons.append("minimum_pair_span")
    if any(seq[i-1]+seq[j-1] not in CANONICAL for i, j in pairs):
        reasons.append("noncanonical")
    endpoints = [x for p in list(g) + pairs for x in p]
    if len(set(endpoints)) != len(endpoints):
        reasons.append("endpoint_reuse")
    if any(crosses(a, b) for k, a in enumerate(pairs) for b in pairs[k+1:]):
        reasons.append("extension_crossing")
    # A fixed round region with no generated endpoint inside must be parsed by
    # V. For its single direct child, the two unpaired wings sum to <= MAXLOOP.
    # Do not apply this V restriction to a pseudoknot band (BE has other bounds).
    over_limit = []
    for i, j in g:
        if any(i < x < j for pair in pairs for x in pair):
            continue
        children = [(k, l) for k, l in g if i < k < l < j and
                    not any(i < a < k < l < b < j for a, b in g)]
        if len(children) == 1:
            k, l = children[0]
            if k-i-1+j-l-1 > 30:
                over_limit.append((i, j, k, l))
    if over_limit:
        reasons.append("V_internal_loop_MAXLOOP_30")
    evaluated = get_eval(raw, label)
    if evaluated["reachable"]:
        # Topology audit and grammar admission are recorded separately. A
        # successful positive evaluator result is never discarded by a proxy.
        return "reachable", reasons
    legacy = raw["result"]["evaluations"].get(raw["aliases"][label]+"_legacy", {})
    if legacy.get("reachable"):
        reasons.append("completed_WMB_weak_closure")
    if not reasons:
        reasons.append("constrained_DP_rejection_unresolved")
    return "rejected", reasons


def distribution(values):
    values = sorted(x for x in values if x is not None and math.isfinite(x))
    if not values:
        return {"n": 0}
    def q(f):
        pos = f * (len(values)-1)
        k = int(pos)
        return values[k] + (values[min(k+1, len(values)-1)]-values[k]) * (pos-k)
    return {"n": len(values), "min": values[0], "q25": q(.25),
            "median": q(.5), "q75": q(.75), "max": values[-1], "mean": st.mean(values)}


def ranks(values):
    result = [0.] * len(values)
    order = sorted(range(len(values)), key=values.__getitem__)
    a = 0
    while a < len(order):
        b = a+1
        while b < len(order) and values[order[b]] == values[order[a]]:
            b += 1
        for k in order[a:b]:
            result[k] = (a+b-1)/2
        a = b
    return result


def spearman(a, b):
    if len(a) < 2 or len(set(a)) < 2 or len(set(b)) < 2:
        return None
    ra, rb = ranks(a), ranks(b)
    ma, mb = st.mean(ra), st.mean(rb)
    return sum((x-ma)*(y-mb) for x, y in zip(ra, rb)) / math.sqrt(
        sum((x-ma)**2 for x in ra)*sum((y-mb)**2 for y in rb))


def write_tsv(path, rows):
    with Path(path).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows({k: ("NA" if v is None else json.dumps(v) if isinstance(v, (dict, list))
                             else int(v) if isinstance(v, bool) else v)
                         for k, v in row.items()} for row in rows)


def aggregate(preps, records, old_results, out, revision="d02d261c",
              previous_probability_source="previous_febd2bc1_exact"):
    old = {r["id"]: r for r in old_results if r["par"] == "turner04"}
    by_id = {r["id"]: salvage_partial_record(r) for r in records}
    rna_rows, pair_rows, term_rows = [], [], []
    max_probability_delta = 0.
    for prepared in preps:
        r = prepared["record"]
        raw = by_id[r["id"]]
        if raw.get("returncode") != 0 and not raw.get("partial_salvaged"):
            raise RuntimeError("cannot summarize incomplete cohort: " + r["id"])
        result = raw["result"]
        probs = {(i, j): p for i, j, p in result["probabilities"]}
        recovered = raw.get("energy_only", False)
        if recovered:
            probs = dict(zip(map(tuple, r["pk"]), old[r["id"]]["pk_probs"]))
        mfe = result["mfe"]
        assert mfe["reachable"]
        mfe_pairs = set(map(tuple, mfe["pairs"]))
        mfe_square = {(i, j) for i, j, o in parse_structure(mfe["generated"]) if o == "["}
        pk = list(map(tuple, r["pk"]))
        row = {"id": r["id"], "n": len(r["seq"]), "n_scaffold": len(r["nested"]),
               "n_reference_pk": len(pk), "noncanonical_removed": prepared["noncanonical_removed"],
               "n_short_canonical": len(prepared["short_canonical"]),
               "n_other_bracket_pairs": len(prepared["other_bracket_pairs"]),
               "median_pk_span": st.median(j-i for i, j in pk),
               "mean_pk_probability": st.mean(probs.get(p, 0.) for p in pk),
               "mfe_energy": mfe["energy"], "mfe_generated_pairs": len(mfe_pairs),
               "mfe_square_pairs": len(mfe_square), "mfe_pk_overlap": len(mfe_pairs & set(pk)),
               "mfe_square_pk_overlap": len(mfe_square & set(pk)),
               "log_z": result["log_z"], "rt_kcal": result["rt_kcal"],
               "probability_source": previous_probability_source if recovered else revision+"_exact"}
        for label in ("full", "filtered", "reported", "scaffold"):
            ev = get_eval(raw, label)
            status, reasons = diagnostics(prepared, label, raw)
            row[label+"_status"] = status
            row[label+"_reasons"] = reasons
            row[label+"_energy"] = ev.get("energy")
            row[label+"_delta_mfe"] = ev.get("energy", mfe["energy"])-mfe["energy"] if ev["reachable"] else None
            row[label+"_log_best_derivation_probability"] = ev["log_weight"]-result["log_z"] if ev["reachable"] and result["log_z"] is not None else None
            if ev["reachable"]:
                wanted = {(i, j) for i, j, _ in parse_structure(prepared["targets"][label])}
                assert set(map(tuple, ev["pairs"])) == wanted
                assert ev["energy"] >= mfe["energy"] - 1e-8
                term_rows.append({"id": r["id"], "target": label, "energy": ev["energy"], **ev["terms"]})
        for label in ("full", "reported"):
            row[label+"_delta_scaffold"] = (row[label+"_energy"]-row["scaffold_energy"]
                if row[label+"_energy"] is not None and row["scaffold_energy"] is not None else None)
        term_rows.append({"id": r["id"], "target": "mfe", "energy": mfe["energy"], **mfe["terms"]})
        # Maximal stacks in the saved, canonical [] target; no gap bridging.
        components = []
        todo = set(pk)
        while todo:
            seed = min(todo)
            component = {seed}
            pending = [seed]
            todo.remove(seed)
            while pending:
                i, j = pending.pop()
                for neighbor in ((i-1, j+1), (i+1, j-1)):
                    if neighbor in todo:
                        todo.remove(neighbor); component.add(neighbor); pending.append(neighbor)
            components.append(component)
        for k, p in enumerate(pk):
            prob = probs.get(p, 0.)
            delta = abs(prob-old[r["id"]]["pk_probs"][k])
            max_probability_delta = max(max_probability_delta, delta)
            ev = get_eval(raw, f"single_{k}")
            status, reasons = diagnostics(prepared, f"single_{k}", raw)
            stem = next(c for c in components if p in c)
            full_ev = get_eval(raw, "full")
            full_contains_pair = full_ev["reachable"] and p in set(map(tuple, full_ev["pairs"]))
            evidence = (revision+"_positive_marginal" if not recovered and prob > 0 else
                        revision+"_full_reference_trace" if full_contains_pair else
                        revision+"_singleton_trace" if ev["reachable"] else "unresolved")
            pair_rows.append({"id": r["id"], "pair_index": k, "i": p[0], "j": p[1],
                              "span": p[1]-p[0], "stack_pairs": len(stem),
                              "stack_median_span": st.median(j-i for i, j in stem),
                              "crosses_scaffold": any(crosses(p, g) for g in map(tuple, r["nested"])),
                              "probability": prob, "probability_source": row["probability_source"],
                              "previous_probability": old[r["id"]]["pk_probs"][k],
                              "probability_delta": delta, "p_gt_half": prob > .5,
                              "positive_marginal": prob > 0,
                              "any_completion_supported": evidence != "unresolved", "membership_evidence": evidence,
                              "single_status": status, "single_reasons": reasons,
                              "single_energy": ev.get("energy"),
                              "full_status": row["full_status"], "reported_status": row["reported_status"],
                              "in_mfe": p in mfe_pairs})
        rna_rows.append(row)
    def fraction(rows):
        n = len(rows)
        return {"pairs": n, "p_gt_half": sum(x["p_gt_half"] for x in rows),
                "fraction_gt_half": sum(x["p_gt_half"] for x in rows)/n if n else None,
                "positive_marginal": sum(x["positive_marginal"] for x in rows),
                "mean_probability": st.mean(x["probability"] for x in rows) if n else None}
    summary = {"rna_count": len(rna_rows), "reference_pair_count": len(pair_rows),
               "pair_groups": {"all": fraction(pair_rows)}, "membership": {},
               "energy_gaps": {}, "distance_groups": {}, "validation": {
                   "public_evaluator_recovers_mfe": len(rna_rows),
                   "max_probability_delta_vs_previous_build": max_probability_delta,
                   "mfe_overlap": sum(x["mfe_pk_overlap"] for x in rna_rows),
                   "mfe_square_total": sum(x["mfe_square_pairs"] for x in rna_rows),
                   "mfe_square_overlap": sum(x["mfe_square_pk_overlap"] for x in rna_rows)}}
    summary["validation"]["fresh_exact_bpp_rnas"] = sum(x["probability_source"] == revision+"_exact" for x in rna_rows)
    summary["validation"]["previous_bpp_rnas"] = sum(x["probability_source"] != revision+"_exact" for x in rna_rows)
    summary["validation"]["fresh_exact_bpp_pairs"] = sum(x["probability_source"] == revision+"_exact" for x in pair_rows)
    summary["validation"]["aggregation_python_sha256"] = sha(__file__)
    summary["validation"]["partial_salvaged_rnas"] = sum(bool(r.get("partial_salvaged")) for r in by_id.values())
    summary["validation"]["timeout_attempts_preserved"] = sum(r.get("error") == "timeout" for r in records)
    summary["incomplete_singleton_tests"] = [
        {k: x[k] for k in ("id", "pair_index", "i", "j", "probability", "single_reasons")}
        for x in pair_rows if x["single_status"] == "unmeasured"]
    checked = []
    max_trace_error = 0.
    max_component_error = 0.
    for raw in by_id.values():
        evaluations = [raw["result"]["mfe"]] + list(raw["result"]["evaluations"].values())
        for ev in evaluations:
            if not ev["reachable"]:
                continue
            checked.append(ev)
            assert abs(ev["energy"] + raw["result"]["rt_kcal"] * ev["log_weight"]) < 1e-8
            max_trace_error = max(max_trace_error, abs(ev["trace_sum_error"]))
            if "terms" in ev:
                component_sum = sum(v for k, v in ev["terms"].items() if not k.endswith("_count"))
                max_component_error = max(max_component_error, abs(component_sum-ev["energy"]))
    summary["validation"].update(reachable_evaluations_audited=len(checked),
        max_trace_sum_error=max_trace_error, max_energy_component_sum_error=max_component_error,
        single_reachable_but_zero_marginal=sum(x["single_status"] == "reachable" and not x["positive_marginal"] for x in pair_rows),
        fresh_bpp_threshold_changes_vs_previous=sum((x["probability"] > .5) != (x["previous_probability"] > .5) for x in pair_rows))
    assert max_trace_error < 1e-8 and max_component_error < 1e-8
    for label in ("full", "filtered", "reported", "scaffold"):
        rejected = [x for x in rna_rows if x[label+"_status"] == "rejected"]
        reason_names = sorted({y for x in rejected for y in x[label+"_reasons"]})
        summary["membership"][label] = {"reachable": len(rna_rows)-len(rejected), "rejected": len(rejected),
            "reasons_nonexclusive": {y: sum(y in x[label+"_reasons"] for x in rejected) for y in reason_names}}
        summary["energy_gaps"][label] = distribution(x[label+"_delta_mfe"] for x in rna_rows)
        if label in ("full", "reported"):
            summary["energy_gaps"][label+"_minus_scaffold"] = distribution(x[label+"_delta_scaffold"] for x in rna_rows)
            summary["energy_gaps"][label]["above_mfe_by_gt_1e_8"] = sum(x[label+"_delta_mfe"] is not None and x[label+"_delta_mfe"] > 1e-8 for x in rna_rows)
    for key in ("single", "full", "reported"):
        for status in ("reachable", "rejected", "unmeasured"):
            group = [x for x in pair_rows if x[key+"_status"] == status]
            summary["pair_groups"][key+"_"+status] = fraction(group)
    summary["pair_groups"]["any_completion_supported"] = fraction([x for x in pair_rows if x["any_completion_supported"]])
    summary["pair_groups"]["any_completion_unresolved"] = fraction([x for x in pair_rows if not x["any_completion_supported"]])
    summary["any_completion_evidence"] = {e: sum(x["membership_evidence"] == e for x in pair_rows)
        for e in sorted({x["membership_evidence"] for x in pair_rows})}
    summary["single_rejection_reasons_nonexclusive"] = {
        reason: sum(reason in x["single_reasons"] for x in pair_rows if x["single_status"] == "rejected")
        for reason in sorted({y for x in pair_rows if x["single_status"] == "rejected" for y in x["single_reasons"]})}
    summary["input_audit"] = {"short_canonical_pairs": sum(x["n_short_canonical"] for x in rna_rows),
        "rna_with_short_canonical": sum(x["n_short_canonical"] > 0 for x in rna_rows),
        "other_bracket_pairs": sum(x["n_other_bracket_pairs"] for x in rna_rows),
        "rna_with_other_brackets": sum(x["n_other_bracket_pairs"] > 0 for x in rna_rows),
        "noncanonical_pairs_removed": sum(x["noncanonical_removed"] for x in rna_rows)}
    for low, high, name in ((0, 50, "span_le_50"), (50, 100, "span_51_100"), (100, 10000, "span_gt_100")):
        group = [x for x in pair_rows if low < x["span"] <= high]
        rnas = [x for x in rna_rows if low < x["median_pk_span"] <= high]
        summary["distance_groups"][name] = {**fraction(group), "rna_count": len(rnas),
            "rna_mean_probability": st.mean(x["mean_pk_probability"] for x in rnas) if rnas else None,
            "full_energy_gap": distribution(x["full_delta_mfe"] for x in rnas),
            "reported_energy_gap": distribution(x["reported_delta_mfe"] for x in rnas),
            "full_reachable_rnas": sum(x["full_status"] == "reachable" for x in rnas),
            "single_reachable_pairs": sum(x["single_status"] == "reachable" for x in group)}
    summary["spearman"] = {
        "rna_median_span_vs_mean_probability": spearman([x["median_pk_span"] for x in rna_rows], [x["mean_pk_probability"] for x in rna_rows]),
        "pair_span_vs_probability": spearman([x["span"] for x in pair_rows], [x["probability"] for x in pair_rows])}
    for label in ("full", "reported"):
        valid = [x for x in rna_rows if x[label+"_delta_mfe"] is not None]
        summary["spearman"][label+"_rna_span_vs_energy_gap"] = spearman([x["median_pk_span"] for x in valid], [x[label+"_delta_mfe"] for x in valid])
    by_term = {(r["id"], r["target"]): r for r in term_rows}
    summary["paired_energy_terms"] = {}
    for label in ("full", "reported"):
        ids = [r["id"] for r in rna_rows if (r["id"], label) in by_term]
        summary["paired_energy_terms"][label] = {term: {
            "reference": distribution(by_term[(rid, label)][term] for rid in ids),
            "mfe": distribution(by_term[(rid, "mfe")][term] for rid in ids),
            "difference": distribution(by_term[(rid, label)][term]-by_term[(rid, "mfe")][term] for rid in ids)}
            for term in ("exterior", "in_multiloop", "in_pseudoloop", "band", "unpaired",
                         "nested_region", "band_stack", "band_internal",
                         "band_multiloop_closing", "band_multiloop_branch", "band_multiloop_unpaired",
                         "multiloop_closing", "multiloop_branch", "multiloop_unpaired",
                         "other_thermodynamic")}
    write_tsv(out/"rna.tsv", rna_rows)
    write_tsv(out/"pairs.tsv", pair_rows)
    write_tsv(out/"energy_terms.tsv", term_rows)
    dump(out/"summary.json", summary)
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cohort", type=Path, default=ROOT/"results/local/pk-recovery-bprna-20261016/pk_sel.json")
    p.add_argument("--previous", type=Path, default=ROOT/"results/local/pk-recovery-bprna-20261016/pk_results.json")
    p.add_argument("--source", type=Path, default=ROOT/"deps/LinearCParty-d02d261c")
    p.add_argument("--build", type=Path, default=ROOT/"build/linear-cparty-d02d261c")
    p.add_argument("--out", type=Path, default=ROOT/"results/local/pk-membership-20261017")
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--timeout", type=float, default=60)
    p.add_argument("--budget-seconds", type=float, default=600)
    p.add_argument("--aggregate-only", action="store_true")
    p.add_argument("--only-id", action="append", default=[], help="restrict a fresh run to selected cohort IDs")
    p.add_argument("--recover-energy-only", action="store_true",
                   help="retry incomplete RNA without recomputing outside; label reused previous BPP")
    args = p.parse_args()
    cohort = json.loads(args.cohort.read_text())
    if args.only_id:
        available = {r["id"] for r in cohort}
        if set(args.only_id)-available:
            p.error("unknown cohort ID")
        cohort = [r for r in cohort if r["id"] in args.only_id]
    preps = [prepare(r) for r in cohort]
    if args.aggregate_only:
        records = [json.loads(line) for line in (args.out/"raw.jsonl").read_text().splitlines()]
        summary = aggregate(preps, records, json.loads(args.previous.read_text()), args.out)
        print(json.dumps(summary, indent=2))
        return
    previous_attempts = []
    if args.recover_energy_only:
        previous_attempts = [json.loads(line) for line in (args.out/"raw.jsonl").read_text().splitlines()]
        latest = {r["id"]: r for r in previous_attempts}
        preps = [x for x in preps if latest.get(x["record"]["id"], {}).get("returncode") != 0]
        if not preps:
            raise RuntimeError("no incomplete RNA to recover")
    if not 1 <= args.workers <= 16:
        p.error("workers must be in 1..16")
    scratch_root = Path(os.environ["TMPDIR"]).resolve()
    if str(scratch_root).startswith("/gs/bs/"):
        p.error("TMPDIR must be node-local")
    if not args.recover_energy_only:
        args.out.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    source_before = command(["git", "-C", args.source, "status", "--porcelain"])
    if source_before:
        raise RuntimeError("pinned source is not clean")
    parameter = (args.source/"params/rna_Turner04.par").resolve()
    cpp = ROOT/"scripts/pk_membership_driver.cc"
    provenance = {"started_utc": datetime.now(timezone.utc).isoformat(), "argv": os.sys.argv,
                  "implementation_revision": command(["git", "-C", args.source, "rev-parse", "HEAD"]),
                  "experiment_revision": command(["git", "-C", ROOT, "rev-parse", "HEAD"]),
                  "experiment_status": command(["git", "-C", ROOT, "status", "--porcelain"]),
                  "source": str(args.source.resolve()), "build": str(args.build.resolve()),
                  "parameter_file": str(parameter), "parameter_sha256": sha(parameter),
                  "cohort_sha256": sha(args.cohort), "previous_results_sha256": sha(args.previous),
                  "python_sha256": sha(__file__), "cpp_sha256": sha(cpp),
                  "library_sha256": {name: sha(args.build/name) for name in ("libcparty_scfg2_exact.a", "libRNA.a")},
                  "policy": "RequireWeakClosure", "mode": "exact", "dangles": 2,
                  "workers": args.workers, "per_rna_timeout_seconds": args.timeout,
                  "wall_budget_seconds": args.budget_seconds, "local_caps": "exact grammar; no beam caps"}
    records = []
    with tempfile.TemporaryDirectory(prefix="pk-membership-", dir=scratch_root) as scratch:
        driver = Path(scratch)/"driver"
        build_cmd = ["/bin/c++", "-std=c++17", "-O2", "-I"+str(args.source.resolve()/"src"),
                     "-I"+str(args.build.resolve()/"generated"), str(cpp),
                     str(args.build.resolve()/"libcparty_scfg2_exact.a"),
                     str(args.build.resolve()/"libRNA.a"), "-o", str(driver)]
        build = subprocess.run(build_cmd, capture_output=True, text=True, timeout=120)
        provenance.update(build_argv=build_cmd, build_returncode=build.returncode,
                          compiler=command(["/bin/c++", "--version"]))
        prefix = "recovery_" if args.recover_energy_only else ""
        dump(args.out/(prefix+"provenance.json"), provenance)
        with tarfile.open(args.out/(prefix+"inputs_and_build.tar.gz"), "w:gz") as tar:
            for path, name in ((args.cohort, "pk_sel.json"), (args.previous, "previous_pk_results.json"),
                               (cpp, "pk_membership_driver.cc"), (Path(__file__), "analyze_pk_membership.py")):
                tar.add(path, arcname=name)
            log = (build.stdout+build.stderr).encode()
            info = tarfile.TarInfo("build.log"); info.size = len(log)
            tar.addfile(info, io.BytesIO(log))
        if build.returncode:
            raise RuntimeError("observer build failed: " + build.stderr[-5000:])
        provenance["driver_sha256"] = sha(driver)
        with (args.out/"raw.jsonl").open("a" if args.recover_energy_only else "w") as f, ThreadPoolExecutor(max_workers=args.workers) as pool:
            # Bound in-flight work; stop submission if the light-run budget is exceeded.
            todo = iter(preps)
            pending = {pool.submit(execute, x, driver, parameter, args.timeout, args.recover_energy_only)
                       for x in [next(todo) for _ in range(min(args.workers, len(preps)))]}
            while pending:
                future = next(as_completed(pending))
                pending.remove(future)
                record = future.result()
                f.write(json.dumps(record, allow_nan=False)+"\n"); f.flush(); os.fsync(f.fileno())
                records.append(record)
                if len(records) % 20 == 0:
                    print(f"completed {len(records)}/{len(preps)}; elapsed {time.monotonic()-start:.1f} s", flush=True)
                if time.monotonic()-start <= args.budget_seconds:
                    nxt = next(todo, None)
                    if nxt is not None:
                        pending.add(pool.submit(execute, nxt, driver, parameter, args.timeout, args.recover_energy_only))
        provenance.update(elapsed_seconds=time.monotonic()-start, completed=len(records),
                          failed=sum(r.get("returncode") != 0 for r in records),
                          source_unchanged=command(["git", "-C", args.source, "status", "--porcelain"]) == source_before,
                          libraries_unchanged=all(sha(args.build/name) == value for name, value in provenance["library_sha256"].items()),
                          finished_utc=datetime.now(timezone.utc).isoformat())
        dump(args.out/(prefix+"provenance.json"), provenance)
    if len(records) != len(preps) or provenance["failed"]:
        raise RuntimeError("stopped/incomplete light run: inspect raw.jsonl and use an array for remaining RNA")
    aggregate([prepare(r) for r in cohort], previous_attempts+records,
              json.loads(args.previous.read_text()), args.out)
    print("complete:", args.out)


if __name__ == "__main__":
    main()
