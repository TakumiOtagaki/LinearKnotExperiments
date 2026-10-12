#!/usr/bin/env python3
"""Small single-thread checks, pinned API path oracle, and immutable dry-run audit.

Compilation and short CLI probes use /tmp only. Nothing is submitted.
"""

import csv
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from make_synthetic_scaffolds import embed_sequence
from plan_array import load_plan, task_jobs
from run_cli import ROOT, digest, verify_lock
from scaffold_features import scaffold_features


def main():
    source = ROOT / "deps/LinearCParty-15aa7742"
    build = ROOT / "build/linear-cparty-15aa7742"
    lock = ROOT / "manifests/software-pilot-15aa7742.json"
    programs = [{"name": n} for n in ("LinearCParty", "LinearPKCapR", "LinearPKRaccess")]
    verify_lock(lock, source, build, programs, source / "params/rna_Turner04.par")
    env = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1", "OMP_NUM_THREADS": "1",
           "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "LC_ALL": "C",
           "CPARTY_FIXED_ENERGY_PARAM_FILE": str(source / "params/rna_Turner04.par")}
    test = subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", str(ROOT / "tests"),
                           "-p", "test_scaffold_features.py", "-v"], env=env, capture_output=True, text=True, check=True)
    probes = []
    with tempfile.TemporaryDirectory(prefix="lcp-e2-shape-", dir="/tmp") as tmp:
        binary = Path(tmp) / "e2-shape-work"
        libraries = [build / n for n in ("libcparty_scfg2_beam.a", "libcparty_scfg2_exact.a", "libRNA.a")]
        command = ["c++", "-O2", "-std=c++17", "-I" + str(source / "src"), str(ROOT / "tools/e2_shape_work.cc"),
                   *map(str, libraries), "-o", str(binary)]
        compile_result = subprocess.run(command, env=env, capture_output=True, text=True, check=True, timeout=60)
        for index, g in enumerate(("...............", "(((...)))", "(.(...)(...).)", "((...)(.(...)(...)))")):
            s = embed_sequence(g, 23)
            record = Path(tmp) / "record.txt"
            record.write_text(f'>probe{index}\n{s}\n{g}\n')
            for d in (0, 1, 2):
                run = subprocess.run([str(binary), "-i", str(record), "3", "10", str(d)], env=env,
                                     capture_output=True, text=True, check=True, timeout=30)
                row = json.loads(run.stdout)
                expected = scaffold_features(s, g, 3, 10, d)["case3_truncated_endpoints"]
                if row["case3_truncated_endpoints"] != expected:
                    raise ValueError("pinned API CASE3 endpoint counter disagrees with feature extractor")
                row["expected_truncated_endpoints"] = expected
                probes.append(row)
        observer_binary_sha256 = digest(binary)
    plans = []
    for kind in ("synthetic-shapes", "rnacentral"):
        path = ROOT / "results/runs" / f"e2-{kind}-15aa7742-capb-timing-20261012" / "array_plan.json"
        plan = load_plan(path)
        if not plan["dry_run"] or not plan["dependencies_verified_at_planning"] or plan["pending_scaffolds"]:
            raise ValueError("unexpected dry-run provenance")
        jobs = []
        for task in plan["tasks"]:
            if task["workers"] != 1 or task["requested_resource"] != "node_f=1":
                raise ValueError("timing allocation mismatch")
            config = json.loads(Path(task["config_path"]).read_text())
            if config["modes"] != [{"name": "beam100_capb", "args": ["-b", "100", "--max-hairpin", "100", "--max-case3-depth", "100"]}]:
                raise ValueError("adopted cap policy mismatch")
            schedule = task_jobs(plan, task)
            warmed = set()
            previous_block = 0
            for job in schedule:
                combination = (job["input"], job["mode"], job["program"])
                if job["warmup"]:
                    warmed.add(combination)
                elif combination not in warmed or job["repetition"] < previous_block:
                    raise ValueError("warmup/repetition block order violated")
                else:
                    previous_block = job["repetition"]
            if len(schedule) != 12 * len(task["inputs"]):
                raise ValueError("expected 3 CLIs x (warmup1 + measured3)")
            if json.loads((path.parent / task["directory"] / "task.json").read_text())["attempts"]:
                raise ValueError("task unexpectedly executed")
            jobs.extend(schedule)
        for script in plan["job_scripts"]:
            if digest(Path(script["path"])) != script["sha256"]:
                raise ValueError("job script changed")
            subprocess.run(["bash", "-n", script["path"]], check=True)
        plans.append({"kind": kind, "plan_sha256": digest(path), "tasks": len(plan["tasks"]),
                      "jobs_including_warmups": len(jobs), "measured_jobs": sum(not j["warmup"] for j in jobs),
                      "estimated_node_hours": sum(t["estimated_seconds"] for t in plan["tasks"]) / 3600,
                      "file_budget": plan["file_budget"]})
    with (ROOT / "results/local/e2-scaffold-features-20261012.tsv").open() as stream:
        features = list(csv.DictReader(stream, delimiter="\t"))
    if len(features) != 26070 or len({(r["cohort"], r["name"]) for r in features}) != 26070:
        raise ValueError("feature cohort count/uniqueness mismatch")
    for row in features:
        if int(row["B_G"]) != int(row["pair_depth_sum"]) or int(row["n"]) != int(row["free_bases"]) + 2 * int(row["scaffold_pairs"]):
            raise ValueError("feature identity mismatch")
    verify_lock(lock, source, build, programs, source / "params/rna_Turner04.par")
    summary = {"schema_version": 1, "unit_tests": test.stderr, "observer_compile_command": command,
               "observer_source_sha256": digest(ROOT / "tools/e2_shape_work.cc"),
               "observer_binary_sha256": observer_binary_sha256,
               "observer_libraries": {p.name: digest(p) for p in libraries},
               "observer_compile_stderr": compile_result.stderr, "api_oracle_probes": probes,
               "features_records": len(features), "plans": plans,
               "software_lock_verified_before_and_after": True, "qsub_calls": 0,
               "tar_extractions_to_filesystem": 0, "scratch": "/tmp, cleaned at exit"}
    destination = ROOT / "results/local/e2-shape-validation-20261012.json"
    destination.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"unit_tests": 5, "api_oracle_probes": len(probes), "features": len(features),
                      "plans": [{k: p[k] for k in ("kind", "tasks", "measured_jobs")} for p in plans]}))


if __name__ == "__main__":
    main()
