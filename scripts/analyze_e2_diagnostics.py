"""Read paired E2 stage processes; counters describe work, not allocated bytes."""

import hashlib
import json
import math
from pathlib import Path

from analyze_e2_synthetic import audit_packed_payload, file_digest, write_union_tsv
from task_storage import RunReader, input_records


MODES = {"inside_outer": ("inside", "outer"), "inside_baseline": ("inside", "baseline"),
         "mfe_outer": ("mfe", "outer"), "outside_bpp_outer": ("outside_bpp", "outer"),
         "loop_profile_outer": ("loop_profile", "outer")}


def read_memory(path, label, workers):
    reader = RunReader(path)
    manifest = reader.manifest
    records = input_records(reader)
    if (manifest["status"] != "success" or manifest["config"]["warmups"] != 0
            or manifest["config"]["repetitions"] != 1 or len(records) != 7):
        raise ValueError("expected seven complete one-shot memory inputs")
    tasks, artifacts = {}, 0
    for task in manifest["tasks"]:
        if (task["status"] != "success" or task["workers"] != workers
                or task["requested_resource"] != "cpu_40=1"
                or len({a["hostname"] for a in task["attempts"]}) != 1):
            raise ValueError("memory task allocation/status mismatch")
        base = Path(path) / task["directory"] / "run"
        metadata = json.loads((base / "run.json").read_text())
        if (file_digest(base / "run.json") != task["run_manifest_sha256"]
                or metadata["binaries_before"] != metadata["binaries_after"]
                or metadata["software_changed_during_run"] or metadata["binaries_changed_during_run"]
                or not metadata["software_lock_valid_after_run"]):
            raise ValueError("memory task dependencies changed")
        if any(metadata["environment"].get(k) != "1" for k in
               ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")):
            raise ValueError("memory task not single threaded")
        durable = [json.loads(line) for line in (base / "task_results.jsonl").open()]
        artifacts += audit_packed_payload(base / "task_outputs.tar", durable)
        aggregate = [j for j in manifest["jobs"] if j["global_task_number"] == task["global_task_number"]]
        if {j["id"] for j in durable} != {j["id"] for j in aggregate}:
            raise ValueError("memory checkpoint denominator mismatch")
        for a, b in ((a, next(b for b in durable if b["id"] == a["id"])) for a in aggregate):
            if any(a[k] != b[k] for k in ("command", "status", "returncode", "peak_rss_kib", "time_wall_seconds")):
                raise ValueError("memory checkpoint mismatch")
        tasks[task["global_task_number"]] = {"task": task["directory"], "hostname": task["attempts"][0]["hostname"],
            "task_manifest_sha256": file_digest(base / "run.json"),
            "checkpoint_sha256": file_digest(base / "task_results.jsonl"),
            "source_revision": metadata["software_before"]["commit"],
            "energy_parameter_sha256": metadata["energy_parameter_sha256"]}
    rows, families, seen = [], [], set()
    for job in manifest["jobs"]:
        key = job["input"], job["mode"]
        if key in seen or job["warmup"] or job["status"] != "success" or job["returncode"] != 0:
            raise ValueError("duplicate/failed memory observation")
        seen.add(key)
        raw = json.loads(reader.artifact(job, "stdout.txt"))
        _, sequence, scaffold = records[job["input"]]
        task = tasks[job["global_task_number"]]
        stage, ranking = MODES[job["mode"]]
        if (not raw["ok"] or raw["source_revision"] != task["source_revision"]
                or raw["input"]["sequence"] != sequence or raw["input"]["scaffold"] != scaffold
                or raw["config"]["ranking"] != ranking or raw["config"]["inside_repeat"]
                or any(raw["config"][k] != 100 for k in ("beam", "max_hairpin", "max_case3_depth"))
                or len(raw["runs"]) != 1 or len(raw["runs"][0]["stages"]) != 1):
            raise ValueError("memory driver input/config mismatch")
        result = raw["runs"][0]["stages"][0]
        if (result["stage"] != stage or not result["ok"] or result["reachable"] is False
                or (result["reachable"] is None and not math.isfinite(result["summary"]["log_z"]))):
            raise ValueError("unreachable memory stage")
        work = result["summary"].get("work", {})
        row = {"run": label, "name": job["input"], "mode": job["mode"], "n": len(sequence),
            "workers": workers, "resource": "cpu_40=1", **task,
            "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
            "scaffold_sha256": hashlib.sha256(scaffold.encode()).hexdigest(),
            "wall_seconds": job["time_wall_seconds"], "rss_kib": job["peak_rss_kib"],
            "rss_gib": job["peak_rss_kib"] / 2**20,
            "driver_vmhwm_kib": result["rss_after"]["vmhwm_kib"],
            "driver_elapsed_seconds": result["elapsed_seconds"],
            **{k: v for k, v in work.items() if k != "families"}}
        row["rss_bytes_per_item_visit"] = row["rss_kib"] * 1024 / work["item_visits"] if work.get("item_visits") else None
        rows.append(row)
        for family in work.get("families", []):
            families.append({"run": label, "name": job["input"], "mode": job["mode"], **family})
        if work and any(sum(f[k] for f in work["families"]) != work[k] for k in
               ("item_visits", "positive_items", "frontier_candidates", "frontier_selected")):
            raise ValueError("family totals do not equal work counters")
    if seen != {(name, mode) for name in records for mode in MODES}:
        raise ValueError("memory input/stage denominator mismatch")
    return rows, families, {"path": str(path), "manifest_sha256": reader.manifest_sha256,
        "processes": len(rows), "raw_artifacts_verified": artifacts, "tasks": list(tasks.values())}


def memory_pairs(old, new):
    lookup = {(r["name"], r["mode"]): r for r in old}
    if set(lookup) != {(r["name"], r["mode"]) for r in new}:
        raise ValueError("memory cohorts differ")
    pairs = []
    for b in new:
        a = lookup[b["name"], b["mode"]]
        if any(a[k] != b[k] for k in ("sequence_sha256", "scaffold_sha256", "energy_parameter_sha256")):
            raise ValueError("memory S/G/parameter mismatch")
        pairs.append({"kind": "new_over_previous", "name_a": a["name"], "name_b": b["name"],
            "mode": b["mode"], "run_a": a["run"], "run_b": b["run"],
            **{k + "_ratio": b[k] / a[k] if a.get(k) and b.get(k) is not None else None for k in
               ("wall_seconds", "rss_gib", "rss_bytes_per_item_visit", "item_visits", "positive_items",
                "binary_join_probes", "extension_probes", "v_parent_probes", "vp_loop_parent_probes")}})
    for label, rows in (("previous", old), ("current", new)):
        for mode in MODES:
            a, b = [next(r for r in rows if r["name"] == name and r["mode"] == mode)
                    for name in ("URS00023C9544", "URS000241043C")]
            pairs.append({"kind": "13652_over_11832", "run_a": a["run"], "run_b": b["run"],
                "name_a": a["name"], "name_b": b["name"], "mode": mode,
                **{k + "_ratio": b[k] / a[k] if a.get(k) and b.get(k) is not None else None for k in
                   ("wall_seconds", "rss_gib", "item_visits", "positive_items", "binary_join_probes",
                    "extension_probes", "v_parent_probes", "vp_loop_parent_probes")}})
    return pairs


def analyze_memory(output, current_path, previous_path, workers):
    new, nf, np = read_memory(current_path, "memory_current", workers)
    old, of, op = read_memory(previous_path, "memory_previous", 1)
    write_union_tsv(output / "memory_metrics.tsv", old + new)
    write_union_tsv(output / "memory_families.tsv", of + nf)
    write_union_tsv(output / "memory_comparisons.tsv", memory_pairs(old, new))
    (output / "memory_summary.json").write_text(json.dumps({"provenance": {"current": np, "previous": op},
        "rss_source": "GNU time peak RSS; driver VmHWM retained separately",
        "counter_scope": "Cumulative main-forward work, including family counts; excludes forcing/pilot/outer. Not allocation bytes or simultaneous live counts. Profile work is unavailable, recorded as NA.",
        "script_sha256": file_digest(__file__),
        "caveat": "Fresh separate processes, one measurement/stage; version, worker count, host and date differ. Stage RSS differences are not additive memory components.",
        "tar_extractions_to_filesystem": 0}, indent=2) + "\n")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    for ax, rows, title in zip(axes, (old, new), ("d21618c5; workers=1", "e2cecbb4; workers=10")):
        for mode in MODES:
            group = sorted((r for r in rows if r["mode"] == mode), key=lambda r: r["n"])
            ax.plot([r["n"] / 1000 for r in group], [r["rss_gib"] for r in group], "o-", ms=3, label=mode)
        ax.set(xlabel="Sequence length (kb)", ylabel="Process peak RSS (GiB)", title=title)
        ax.grid(alpha=.2)
    axes[1].legend(fontsize=7)
    fig.suptitle("RNAcentral stage processes: separate one-shot measurements; panel y scales differ", fontsize=10)
    fig.savefig(output / "memory_stages.png", dpi=170)
    plt.close(fig)
