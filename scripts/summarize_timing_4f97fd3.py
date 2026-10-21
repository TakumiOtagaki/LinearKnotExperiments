#!/usr/bin/env python3
"""Audit the 4f97fd3 timing rerun and export paired resource tables.

Run after aggregate_array.py and the E2/PF analyses. Read packed artifacts
without extracting them; retain the old scientific results separately.
"""
from collections import Counter
import argparse
import csv
import json
from pathlib import Path
from statistics import median

from analyze_e2_synthetic import audit_packed_payload
from run_cli import digest

ROOT = Path(__file__).resolve().parents[1]
NEW = "4f97fd3-20261021"
OLD = "89baea52-20261018"


def table(path, rows):
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def read_table(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter="\t"))


def neat_resources(kind, version, verify=False):
    run = ROOT / "results/runs" / f"e3-neat1-{kind}-{version}"
    m = json.loads((run / "run.json").read_text())
    rows = []
    for task in m["tasks"]:
        base = run / task["directory"] / "run"
        meta = json.loads((base / "run.json").read_text())
        assert meta["status"] == "success"
        assert meta["software_before"] == meta["software_after"]
        assert meta["binaries_before"] == meta["binaries_after"]
        assert meta["software_lock_valid_after_run"] is True
        assert not meta["software_changed_during_run"] and not meta["binaries_changed_during_run"]
        durable = [json.loads(s) for s in (base / "task_results.jsonl").read_text().splitlines()]
        planned = json.loads((run / task["directory"] / "planned_jobs.json").read_text())
        assert len(durable) == len(planned) == len({r["id"] for r in durable})
        assert {r["id"] for r in durable} == {r["id"] for r in planned}
        if verify:
            audit_packed_payload(base / "task_outputs.tar", durable)
        for r in durable:
            assert r["status"] == "success" and r["returncode"] == 0
            assert not r["warmup"] and r["repetition"] == 0
            assert r["time_wall_seconds"] > 0 and r["peak_rss_kib"] > 0
            b = int(r["mode"].split("_")[0].removeprefix("beam"))
            rows.append(dict(iso=r["input"], beam=b, analysis=kind,
                             wall_s=r["time_wall_seconds"], rss_gib=r["peak_rss_kib"] / 2**20,
                             input_sha256=next(i["sha256"] for i in m["inputs"] if i["name"] == r["input"]),
                             parameter_sha256=meta["energy_parameter_sha256"],
                             checkpoint_sha256=digest(base / "task_results.jsonl"),
                             command=json.dumps(r["command"])))
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "results/local/timing-4f97fd3-20261021")
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(exist_ok=False)
    audit = []
    for run in sorted((ROOT / "results/runs").glob(f"*-{NEW}")):
        m = json.loads((run / "run.json").read_text())
        assert m["status"] == "success" and m["missing_measured_jobs"] == 0
        assert all(j["status"] == "success" for j in m["jobs"])
        assert m["planned_measured_jobs"] == m["completed_measured_jobs"]
        audit.append(dict(run=run.name, tasks=len(m["tasks"]), processes=len(m["jobs"]),
                          measured=m["completed_measured_jobs"], missing=0, failed=0,
                          manifest_sha256=digest(run / "run.json")))
    assert len(audit) == 12 and sum(r["tasks"] for r in audit) == 82
    assert sum(r["processes"] for r in audit) == 2768
    table(out / "run_audit.tsv", audit)
    neat, comparisons = [], []
    for kind in ("bpp", "capr", "raccess", "raccess-byloop", "sampling", "centroid"):
        a = neat_resources(kind, "89baea52-20261020" if kind == "centroid" else OLD)
        b = neat_resources(kind, NEW, verify=True)
        keys = lambda r: (r["iso"], r["beam"], r["analysis"])
        prior = {keys(r): r for r in a}
        assert set(prior) == {keys(r) for r in b}
        neat.extend(b)
        for r in b:
            p = prior[keys(r)]
            assert p["input_sha256"] == r["input_sha256"]
            assert p["parameter_sha256"] == r["parameter_sha256"]
            comparisons.append(dict(scope="NEAT1", input=r["iso"], mode=r["beam"], analysis=kind,
                                    old_wall_s=p["wall_s"], new_wall_s=r["wall_s"],
                                    old_rss_gib=p["rss_gib"], new_rss_gib=r["rss_gib"]))
    assert len(neat) == 33
    table(out / "neat1_resources.tsv", neat)
    e2 = ROOT / "results/local/e2-4f97fd3-20261021"
    metrics = read_table(e2 / "metrics.tsv")
    for scope in ("rnacentral", "synthetic"):
        lookup = {(r["name"], r["program"]): r for r in metrics if r["run"] == scope + "_previous"}
        for r in metrics:
            if r["run"] != scope + "_current":
                continue
            p = lookup[r["name"], r["program"]]
            assert all(p[k] == r[k] for k in ("sequence_sha256", "scaffold_sha256"))
            comparisons.append(dict(scope=scope, input=r["name"], mode=100, analysis=r["program"],
                                    old_wall_s=float(p["wall_seconds_median"]), new_wall_s=float(r["wall_seconds_median"]),
                                    old_rss_gib=float(p["rss_gib_median"]), new_rss_gib=float(r["rss_gib_median"])))
    memory = read_table(e2 / "memory_metrics.tsv")
    prior = {(r["name"], r["mode"]): r for r in memory if r["run"] == "memory_previous"}
    for r in memory:
        if r["run"] == "memory_current":
            p = prior[r["name"], r["mode"]]
            comparisons.append(dict(scope="memory", input=r["name"], mode=100, analysis=r["mode"],
                                    old_wall_s=float(p["wall_seconds"]), new_wall_s=float(r["wall_seconds"]),
                                    old_rss_gib=float(p["rss_gib"]), new_rss_gib=float(r["rss_gib"])))
    for workers in (1, 10):
        run = ROOT / "results/runs" / f"e2-parallel-impact-w{workers}-{OLD}"
        previous = json.loads((run / "run.json").read_text())
        for r in metrics:
            if r["run"] != f"parallel_w{workers}":
                continue
            jobs = [j for j in previous["jobs"] if j["input"] == r["name"]
                    and j["program"] == r["program"] and not j["warmup"]]
            assert len(jobs) == 3 and all(j["status"] == "success" for j in jobs)
            comparisons.append(dict(scope=f"parallel_w{workers}", input=r["name"], mode=100,
                                    analysis=r["program"], old_wall_s=median(j["time_wall_seconds"] for j in jobs),
                                    new_wall_s=float(r["wall_seconds_median"]),
                                    old_rss_gib=median(j["peak_rss_kib"] for j in jobs) / 2**20,
                                    new_rss_gib=float(r["rss_gib_median"])))
    assert len(comparisons) == 575
    for r in comparisons:
        r.update(new_over_old_wall=r["new_wall_s"] / r["old_wall_s"],
                 new_over_old_rss=r["new_rss_gib"] / r["old_rss_gib"])
    table(out / "version_comparison.tsv", comparisons)
    pf = json.loads((ROOT / "results/local/cparty-linearknot-pf-4f97fd3-20261021/summary.json").read_text())
    failed = [r for r in pf["observations"] if r["mode"] == "cparty" and r["status"] != "success"]
    assert len(failed) == 192
    assert all(r["stderr_summary"]["tail"].strip() == "nonfinite oracle inside" for r in failed)
    summaries = []
    for n in (250, 500, 1000, 2000, 4000):
        for method in ("cparty", "exact", "beam100"):
            all_rows = [r for r in pf["resources"] if r["length"] == n and r["method"] == method]
            paired = [r for r in pf["paired_resources"] if r["length"] == n and r["method"] == method]
            full = [r for r in all_rows if r["status"] == "success"]
            r = dict(length=n, method=method, planned_inputs=len(all_rows), successful_inputs=len(full),
                     nonfinite_inputs=sum(v["status"] == "unavailable" for v in all_rows), paired_inputs=len(paired))
            for label, group in (("paired", paired), ("all_successful", full)):
                if group:
                    r.update({label + "_wall_s": median(v["elapsed_seconds"] for v in group),
                              label + "_rss_gib": median(v["peak_rss_kib"] for v in group) / 2**20})
            summaries.append(r)
    table(out / "pf_by_length.tsv", summaries)
    # PF version comparison uses matched inputs and the old PF-only observations.
    from analyze_pf_scaling import summarize, read_observations
    from plan_array import load_plan
    oldplan = load_plan(ROOT / "results/runs/cparty-linearknot-pf-89baea52-20261020/array_plan.json")
    obs, _ = read_observations(oldplan)
    oldrows, _, _ = summarize(obs, oldplan["inputs"], oldplan["config"]["modes"])
    prior = {(r["input"], r["method"]): r for r in oldrows}
    pf_comparisons = []
    for r in pf["resources"]:
        if r["method"] == "cparty" or r["status"] != "success":
            continue
        p = prior[r["input"], r["method"]]
        assert p["status"] == "success"
        pf_comparisons.append(dict(input=r["input"], length=r["length"], method=r["method"],
                                   old_wall_s=p["elapsed_seconds"], new_wall_s=r["elapsed_seconds"],
                                   old_rss_gib=p["peak_rss_kib"] / 2**20, new_rss_gib=r["peak_rss_kib"] / 2**20))
    assert len(pf_comparisons) == 168
    table(out / "pf_version_comparison.tsv", pf_comparisons)
    (out / "validation.json").write_text(json.dumps(dict(status="passed", matrices=12, tasks=82,
        processes=2768, failures=0, missing=0, neat1_conditions=33, pf_conditions=168,
        version_comparison_rows=len(comparisons),
        cparty_nonfinite_inputs=48, cparty_nonfinite_processes=192,
        cparty_exact_parity=dict(Counter(r["status"] for r in pf["parity"])),
        script_sha256=digest(Path(__file__)), tar_extractions_to_filesystem=0), indent=2) + "\n")
    print(out)


if __name__ == "__main__":
    main()
