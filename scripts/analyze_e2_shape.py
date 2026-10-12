#!/usr/bin/env python3
"""Join isolated E2 timings to G features; E1 is descriptive parallel-run evidence.

Reads existing summaries only; does not unpack tar files or invoke RNA CLIs.
Writes one joined TSV, one JSON summary, and one standalone figure.
"""

import argparse
from collections import defaultdict
import csv
import json
import math
from pathlib import Path
import statistics
import tarfile

from scaffold_features import scaffold_features, sha256, write_tsv
from task_storage import MemoryFile, RunReader, input_records
from rna_records import read_record

FEATURES = ("n", "free_bases", "B_G", "pair_depth_max", "free_run_max", "free_runs",
            "loop_free_max", "same_loop_pairs_turn", "v_geometry_beam_cap_sum",
            "hairpin_free_geometries", "cross_loop_canonical_pairs_turn",
            "case3_sibling_endpoint_pairs", "case3_path_entries",
            "case3_ancestors_max", "max_pair_children")


def ranks(values):
    order = sorted(range(len(values)), key=values.__getitem__)
    result = [0.] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        for i in order[start:end]:
            result[i] = (start + end - 1) / 2
        start = end
    return result


def spearman(xs, ys):
    if len(xs) < 3 or len(xs) != len(ys):
        return None
    x, y = ranks(xs), ranks(ys)
    mx, my = statistics.mean(x), statistics.mean(y)
    xx, yy = sum((v - mx)**2 for v in x), sum((v - my)**2 for v in y)
    if not xx or not yy:
        return None
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(xx * yy)


def summarize(values):
    return {"median": statistics.median(values), "min": min(values), "max": max(values)}


def analyze(features_path, timing_run, output, e1_metrics=None):
    output = Path(output)
    if output.exists():
        raise ValueError("use a new analysis output directory")
    with Path(features_path).open() as stream:
        features = {r["name"]: r for r in csv.DictReader(stream, delimiter="\t")}
    reader = RunReader(timing_run)
    manifest = reader.manifest
    if manifest["config"].get("warmups") != 1 or manifest["config"].get("repetitions") != 3:
        raise ValueError("expected warmup 1 and repetitions 3")
    if any(t.get("workers") != 1 or t.get("requested_resource") != "node_f=1" for t in manifest["tasks"]):
        raise ValueError("E2 timing must have isolated node_f and one worker")
    records = input_records(reader)
    joined = []
    excluded = []
    for name, (_, sequence, scaffold) in records.items():
        modes = manifest["config"]["modes"]
        if len(modes) != 1:
            raise ValueError("expected one measured mode")
        args = modes[0]["args"]
        b, h, d = [int(args[args.index(flag) + 1]) for flag in ("-b", "--max-hairpin", "--max-case3-depth")]
        # Recompute under actual historical 30/30 caps, not adopted 100/100.
        f = scaffold_features(sequence, scaffold, b, h, d)
        stored = features[name]
        if any(stored[k] != f[k] for k in ("sequence_sha256", "scaffold_sha256")):
            raise ValueError("feature/timing input hash mismatch: " + name)
        for program in manifest["config"]["programs"]:
            jobs = [j for j in manifest["jobs"] if j["input"] == name and
                    j["program"] == program["name"] and not j["warmup"]]
            if len(jobs) != 3 or sorted(j["repetition"] for j in jobs) != [0, 1, 2]:
                raise ValueError("missing or duplicate measured repetitions")
            ok = [j for j in jobs if j["status"] == "success" and
                  j.get("time_wall_seconds") is not None and j.get("peak_rss_kib") is not None]
            excluded.extend({"input": name, "program": program["name"], "repetition": j["repetition"],
                             "status": j["status"]} for j in jobs if j not in ok)
            if not ok:
                continue
            time = summarize([j["time_wall_seconds"] for j in ok])
            rss = summarize([j["peak_rss_kib"] / 2**20 for j in ok])
            joined.append({"cohort": "rnacentral", "name": name, "program": program["name"],
                           "measured_repetitions": len(ok),
                           **{"wall_seconds_" + k: v for k, v in time.items()},
                           **{"rss_gib_" + k: v for k, v in rss.items()}, **f})
    correlations = []
    for program in manifest["config"]["programs"]:
        group = [r for r in joined if r["program"] == program["name"]]
        for feature in FEATURES:
            for outcome in ("wall_seconds_median", "rss_gib_median"):
                correlations.append({"program": program["name"], "feature": feature, "outcome": outcome,
                                     "n_inputs": len(group),
                                     "spearman": spearman([r[feature] for r in group], [r[outcome] for r in group])})
    comparisons = {}
    for program in manifest["config"]["programs"]:
        by_id = {r["name"]: r for r in joined if r["program"] == program["name"]}
        a, b = by_id["URS00023C9544"], by_id["URS000241043C"]
        comparisons[program["name"]] = {k: b[k] / a[k] for k in FEATURES + ("wall_seconds_median", "rss_gib_median")
                                            if a[k] != 0}
    e1_summary = {}
    e1_validation = {}
    if e1_metrics:
        provenance_path = Path(e1_metrics).with_name("provenance.json")
        provenance = json.loads(provenance_path.read_text())
        if sha256(Path(e1_metrics).read_bytes()) != provenance["metrics_sha256"]:
            raise ValueError("E1 metrics/provenance hash mismatch")
        if provenance["main_run"]["software_commit"] != manifest["software_before"]["commit"]:
            raise ValueError("E1/E2 fixed version mismatch")
        rows = []
        seen = set()
        with Path(e1_metrics).open() as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                if row["mode"] != "beam100_default":
                    continue
                name = row["input"]
                if name in seen:
                    raise ValueError("duplicate E1 input/mode")
                seen.add(name)
                f = features[name]
                if (any(int(row[a]) != int(f[b]) for a, b in
                        (("length", "n"), ("scaffold_pairs", "scaffold_pairs"), ("free_bases", "free_bases")))
                        or row["analysis_status"] != "success"):
                    raise ValueError("unexpected E1 identity/status: " + name)
                f = dict(f)
                # Evaluate cap-dependent descriptors under the measured 30/30.
                hp = 0
                for length, count in json.loads(f["free_run_length_hist"]).items():
                    length, high = int(length), min(31, int(length) - 1)
                    if high >= 4:
                        hp += count * ((high - 3) * length - (4 + high) * (high - 3) // 2)
                f["hairpin_free_geometries"] = hp
                f["case3_path_entries"] = sum(min(30, int(k)) * v for k, v in json.loads(f["case3_ancestor_hist"]).items())
                rows.append((row, f))
        # Verify S/G identities from one packed E1 input archive, in memory.
        archive_path = Path(provenance["main_run"]["path"]) / "inputs.tar"
        if sha256(archive_path.read_bytes()) != provenance["main_run"]["inputs_archive_sha256"]:
            raise ValueError("E1 frozen input archive hash mismatch")
        verified = set()
        with tarfile.open(archive_path, "r|") as archive:
            for member in archive:
                if not member.isfile() or Path(member.name).stem not in seen:
                    continue
                data = archive.extractfile(member).read()
                name, sequence, scaffold = read_record(MemoryFile(data))
                f = features[name]
                if name in verified or (sha256(sequence.encode()) != f["sequence_sha256"] or
                                        sha256(scaffold.encode()) != f["scaffold_sha256"]):
                    raise ValueError("E1 feature/input identity mismatch: " + name)
                verified.add(name)
        if verified != seen:
            raise ValueError("E1 frozen archive missing selected inputs")
        e1_validation = {"verified_sequence_scaffold_hashes": len(verified), "archive_sha256": provenance["main_run"]["inputs_archive_sha256"],
                         "metrics_sha256": provenance["metrics_sha256"], "provenance_sha256": sha256(provenance_path.read_bytes()),
                         "descriptor_caps": {"h": 30, "d": 30}, "tar_extractions_to_filesystem": 0}
        # Within 25-nt bins subtract each variable's mean rank, then correlate
        # the residual ranks. No speed ranking or fitted asymptotic law.
        for program in ("LinearCParty", "LinearPKCapR", "LinearPKRaccess"):
            for measure in ("time_wall_seconds", "peak_rss_kib"):
                key = program + "_" + measure
                selected = [(r, f) for r, f in rows if r[key] not in ("", "NA") and float(r[key]) > 0]
                tests = []
                for feature in FEATURES:
                    xs, ys = [], []
                    bins = defaultdict(list)
                    for row, f in selected:
                        bins[int(f["n"]) // 25].append((float(f[feature]), float(row[key])))
                    for group in bins.values():
                        x, y = ranks([p[0] for p in group]), ranks([p[1] for p in group])
                        if len(group) >= 3:
                            avg = (len(group) - 1) / 2
                            # Normalize rank range so large bins do not dominate.
                            xs.extend((v - avg) / len(group) for v in x)
                            ys.extend((v - avg) / len(group) for v in y)
                    tests.append({"feature": feature, "rank_residual_spearman": spearman(xs, ys), "n": len(xs)})
                e1_summary[key] = {"input_rows": len(rows), "positive_observations": len(selected),
                                   "excluded_zero_or_missing": len(rows) - len(selected),
                                   "within_25nt_bins": tests}
    output.mkdir(parents=True)
    write_tsv(output / "rnacentral_joined.tsv", joined)
    summary = {"schema_version": 1, "scope": "Descriptive association only; seven fixed real RNAs, S and G both differ.",
               "historical_caps": {"b": b, "h": h, "d": d}, "successful_measured_jobs": sum(r["measured_repetitions"] for r in joined),
               "excluded_jobs": excluded, "rnacentral_correlations": correlations,
               "13652_vs_11832_ratios": comparisons, "e1_parallel_short_trends": e1_summary,
               "e1_validation": e1_validation,
               "e1_caveat": "50-399 nt subset, old beam100_default h=d=30, parallel accuracy runs. Timing granularity and shared-node contention matter; no isolated speed claim.",
               "limitations": ["No internal live-memory, candidate or forced-prepass counters in the historical CLI outputs.",
                               "G-only geometry supplies bounds and hypotheses, not allocator byte counts or retained numeric work.",
                               "Raw seven-point correlations are strongly confounded by n; no regression/p-value or causal claim."],
               "provenance": {"features_sha256": sha256(Path(features_path).read_bytes()),
                              "timing_run_sha256": reader.manifest_sha256,
                              "e1_metrics_sha256": sha256(Path(e1_metrics).read_bytes()) if e1_metrics else None,
                              "script_sha256": sha256(Path(__file__).read_bytes())}}
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for program in manifest["config"]["programs"]:
        group = [r for r in joined if r["program"] == program["name"]]
        for ax, prefix, ylabel in zip(axes, ("wall_seconds", "rss_gib"), ("CLI wall time (s)", "Peak RSS (GiB)")):
            x, y = [r["n"] / 1000 for r in group], [r[prefix + "_median"] for r in group]
            error = [[r[prefix + "_median"] - r[prefix + "_min"] for r in group],
                     [r[prefix + "_max"] - r[prefix + "_median"] for r in group]]
            ax.errorbar(x, y, yerr=error, fmt="o-", ms=4, lw=1, label=program["name"])
            ax.set(xlabel="Sequence length (kb)", ylabel=ylabel)
            ax.grid(alpha=.2)
    axes[0].legend(fontsize=8)
    axes[1].annotate("URS000241043C", (13.652, 28.55), xytext=(9, 37),
                     arrowprops={"arrowstyle": "->"}, fontsize=8)
    fig.suptitle("Fixed 15aa7742, b=100, h=d=30; median and range of 3 isolated runs", fontsize=10)
    fig.savefig(output / "rnacentral_time_rss.png", dpi=180)
    plt.close(fig)
    print(json.dumps({"joined_rows": len(joined), "e1_rows": len(rows) if e1_metrics else 0,
                      "successful_measured_jobs": summary["successful_measured_jobs"]}))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", required=True)
    parser.add_argument("--timing-run", required=True)
    parser.add_argument("--e1-metrics")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    analyze(args.features, args.timing_run, args.out, args.e1_metrics)


if __name__ == "__main__":
    main()
