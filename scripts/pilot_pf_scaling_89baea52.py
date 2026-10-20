#!/usr/bin/env python3
"""Small local calibration only; never submit a scheduler job."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import random
import socket
import tarfile

from run_cli import ROOT, now, run_one, verify_lock, write_json
from task_storage import local_workspace


def pilot(output):
    if output.exists():
        raise ValueError("use a new pilot directory")
    source = ROOT / "deps/LinearCParty-89baea52"
    build = ROOT / "build/pf-scaling-89baea52-20261020"
    lock = ROOT / "manifests/software-pf-scaling-89baea52-20261020.json"
    verify_lock(lock, source, build, [{"name": "PFScaling"}], source / "params/rna_Turner04.par")
    e2 = ROOT / "results/runs/e2-synthetic-shapes-89baea52-20261018"
    plan = json.loads((e2 / "array_plan.json").read_text())
    inputs = [x for x in plan["inputs"] if x["length"] <= 1000 and x["seed"] == 2026101201
              and x["family"] in ("empty", "long_stem", "balanced_branching")]
    schedule = []
    rng = random.Random(20261020)
    for rep in range(3):
        block = inputs.copy(); rng.shuffle(block)
        for item in block:
            if rep == 0:
                schedule.append((item, "cparty", True, 0))
            schedule.append((item, "cparty", False, rep))
    # Additional one-shot LinearKnot feasibility checks, deliberately avoiding
    # repetition of the known 120s exact timeout at n=1000.
    for item in inputs:
        if item["length"] == 1000:
            schedule.append((item, "beam100", False, 0))
        if item["length"] == 250 and item["family"] in ("empty", "long_stem"):
            schedule.append((item, "exact", False, 0))
        if item["length"] == 500 and item["family"] == "empty":
            schedule.append((item, "exact", False, 0))
    result = {"purpose": "development calibration only", "workers": 1, "host": socket.gethostname(),
              "started_utc": now(), "timeout_seconds": 120, "memory_limit_mib": 4096,
              "cparty_warmups": 1, "cparty_repetitions": 3,
              "linear_checks": "one-shot feasibility controls, not formal measurements", "observations": []}
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1",
           "LC_ALL": "C", "CPARTY_FIXED_ENERGY_PARAM_FILE": str(source / "params/rna_Turner04.par")}
    output.mkdir(parents=True)
    with local_workspace() as work, tarfile.open(e2 / "inputs.tar") as ar, tarfile.open(output / "raw.tar", "x") as raw:
        for item in inputs:
            data = ar.extractfile(item["input_member"]).read()
            if hashlib.sha256(data).hexdigest() != item["sha256"]:
                raise ValueError("frozen E2 input hash differs")
            (work / (item["name"] + ".txt")).write_bytes(data)
        for index, (item, method, warmup, repetition) in enumerate(schedule):
            directory = work / f"job{index:04d}"; directory.mkdir()
            command = [str(build / "PFScaling"), "-i", str(work / (item["name"] + ".txt")), "--method", method]
            row = run_one(command, directory, env, 120, 4096)
            row.update({k: item[k] for k in ("name", "length", "family", "seed", "sha256")})
            row.update(method=method, warmup=warmup, repetition=repetition,
                       stdout=(directory / "stdout.txt").read_text(), stderr=(directory / "stderr.txt").read_text())
            result["observations"].append(row)
            raw.add(directory, arcname=f"job{index:04d}")
            raw.fileobj.flush()
            write_json(output / "pilot.json", result)
            print(item["length"], item["family"], method, "warmup" if warmup else repetition,
                  row["status"], round(row["elapsed_seconds"], 3), row["peak_rss_kib"], flush=True)
    result["finished_utc"] = now()
    verify_lock(lock, source, build, [{"name": "PFScaling"}], source / "params/rna_Turner04.par")
    result["software_lock_valid_after_run"] = True
    write_json(output / "pilot.json", result)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    pilot(args.out.resolve())
