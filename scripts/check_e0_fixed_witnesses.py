#!/usr/bin/env python3
"""Check five PF witnesses and the WIP change; keep per-job files under /tmp."""

import argparse
import json
import os
from pathlib import Path
import tarfile
import tempfile

from analyze_e0 import analyze, compare, read_observation, within_tolerance
from run_cli import ROOT, digest, execute, experiment_path, git_state, now, run_one, write_json
from task_storage import MemoryFile, RunReader

CONFIG = ROOT / "configs/e0_witnesses_d21618c5_5patches_20261014.json"
SOURCE = ROOT / "deps/LinearCParty-d21618c5"
BUILD = ROOT / "build/e0-observers-d21618c5-800e95a-5patches"
OLD_SOURCE = ROOT / "deps/LinearCParty-15aa7742"
OLD_BUILD = ROOT / "build/e0-observers-15aa7742-800e95a-4patches"
OLD_LOCK = ROOT / "manifests/software-e0-15aa7742-20261011.json"
OLD_ORACLE_LOCK = ROOT / "manifests/software-cparty-oracle-800e95a-4patches-20261011.json"


def verify_historical(config):
    """Verify retained artifacts, without requiring edited experiment scripts.

    Old locks remain unchanged. Their builder/observer source files are historical
    provenance, not the source of a rebuild in this check. Verify every retained
    build/dependency entry, including the scalar binary, libraries and objects.
    """
    lock = json.loads(OLD_LOCK.read_text())
    oracle = json.loads(OLD_ORACLE_LOCK.read_text())
    if git_state(OLD_SOURCE) != lock["source_git_state"]:
        raise ValueError("historical LinearCParty source changed")
    oracle_source = ROOT / "build/cparty-oracle-800e95a-4patches/src"
    if git_state(oracle_source) != oracle["source_git_state"]:
        raise ValueError("historical oracle source changed")
    required = [(OLD_BUILD / "CMakeCache.txt", lock["cmake_cache_sha256"])]
    required += [(OLD_BUILD / name, sha) for name, sha in lock["binaries"].items()]
    required += [(ROOT / entry["path"], entry["sha256"]) for entry in lock["extra_files"]
                 if entry["path"].startswith(("build/", "deps/"))]
    parameter = SOURCE / config["energy_parameter_file"]
    required += [(parameter, lock["energy_parameter_sha256"]),
                 (OLD_SOURCE / config["energy_parameter_file"], lock["energy_parameter_sha256"])]
    for path, sha in required:
        if sha is None or digest(path) != sha:
            raise ValueError("historical artifact changed: " + str(path))
    return {"source_commit": lock["source_commit"], "oracle": lock["oracle"],
            "lock_sha256": digest(OLD_LOCK), "oracle_lock_sha256": digest(OLD_ORACLE_LOCK),
            "verified_artifacts": [{"path": str(path.relative_to(ROOT)), "sha256": sha}
                                   for path, sha in required],
            "scope": "Retained historical binaries, cache, libraries, oracle manifest/sources/objects, owner config/patches and parameters; edited experiment scripts are not rebuilt"}


def wip_comparison(new, old, settings):
    baseline = compare(old["E0LegacyCompatible"], old["E0CPartyOracle"])
    current = compare(new["E0LegacyCompatible"], new["E0CPartyOracle"])
    changes = {program: compare(new[program], old[program]) for program in old}
    return {"old_compatible_minus_oracle": baseline, "new_compatible_minus_oracle": current,
            "new_minus_old": changes,
            "passed": (within_tolerance(baseline, settings) and within_tolerance(current, settings)
                       and all(not within_tolerance(value, settings)
                               and value["delta_log_z"] > 0 and value["delta_free_energy"] < 0
                               for value in changes.values()))}


def check(output):
    if output.exists():
        raise ValueError("use a new witness output directory")
    config = json.loads(CONFIG.read_text())
    historical = verify_historical(config)
    output.mkdir(parents=True)
    with tempfile.TemporaryDirectory(prefix="e0-fixed-witnesses-", dir="/tmp") as temporary:
        workspace = Path(temporary)
        run = workspace / "new"
        code = execute(CONFIG, SOURCE, BUILD, run)
        reader = RunReader(run)
        analysis = analyze(run)
        observations = []
        new_wip = {}
        for job in reader.manifest["jobs"]:
            row = {"input": job["input"], "program": job["program"], "status": job["status"],
                   "command": job["command"], "result": job}
            if job["status"] == "success":
                value = read_observation(MemoryFile(reader.artifact(job, "stdout.txt")))
                row["value"] = value
                if job["input"] == "WIP29":
                    new_wip[job["program"]] = value
            row["raw_sha256"] = {name: digest(run / job["directory"] / name)
                                 for name in ("stdout.txt", "stderr.txt", "command.json", "result.json", "time.tsv")}
            observations.append(row)
        old_wip, old_jobs = {}, []
        record = ROOT / next(entry["path"] for entry in config["inputs"] if entry["name"] == "WIP29")
        env = {**os.environ, **config["environment"],
               "CPARTY_FIXED_ENERGY_PARAM_FILE": str(SOURCE / config["energy_parameter_file"])}
        for program in ("E0LegacyCompatible", "E0CPartyOracle"):
            directory = workspace / "old" / program
            directory.mkdir(parents=True)
            argv = [str(OLD_BUILD / program), "-i", str(record)]
            result = run_one(argv, directory, env, 60, 4096)
            write_json(directory / "result.json", result)
            row = {"program": program, "command": argv, "result": result,
                   "raw_sha256": {name: digest(directory / name) for name in
                                  ("stdout.txt", "stderr.txt", "command.json", "result.json", "time.tsv")}}
            if result["status"] == "success":
                old_wip[program] = read_observation(directory / "stdout.txt")
                row["value"] = old_wip[program]
            old_jobs.append(row)
        regression = (wip_comparison(new_wip, old_wip, config["analysis"])
                      if all(p in new_wip and p in old_wip for p in ("E0LegacyCompatible", "E0CPartyOracle"))
                      else {"passed": False, "status": "observation_unavailable"})
        verify_historical(config)
        archive = output / "witness_run.tar"
        with tarfile.open(archive, "w") as stream:
            stream.add(workspace, arcname="e0-fixed-witnesses")
        # Check that the retained tar can be analyzed without extracted files.
        archived = analyze(str(archive) + "::e0-fixed-witnesses/new")
        for key in ("numerical_parity_counts", "closure_counts", "observations"):
            if analysis[key] != archived[key]:
                raise ValueError("archived witness analysis differs: " + key)
        parity = [r for r in analysis["observations"] if r["comparison_kind"] == "numerical_parity"]
        passed = (code == 0 and len(parity) == 5 and all(r["status"] == "success" for r in parity)
                  and regression["passed"] and analysis["closure_counts"]["restriction_direction_inconsistent"] == 0)
        summary = {"schema_version": 1, "created_utc": now(), "purpose": "development", "passed": passed,
                   "script_sha256": digest(Path(__file__)), "config_sha256": digest(CONFIG),
                   "new_lock_sha256": digest(ROOT / config["software_lock"]), "historical": historical,
                   "raw_archive_path": str(archive), "raw_archive_sha256": digest(archive),
                   "new_jobs": len(observations), "old_jobs": len(old_jobs),
                   "parameter_sha256": digest(SOURCE / config["energy_parameter_file"]),
                   "max_absolute_delta_log_z": max((abs(r["metrics"]["delta_log_z"]) for r in parity if "metrics" in r), default=None),
                   "max_absolute_delta_free_energy": max((abs(r["metrics"]["delta_free_energy"]) for r in parity if "metrics" in r), default=None),
                   "wip_regression": regression, "observations": observations, "old_wip_observations": old_jobs,
                   "scope": "Five small exact-inside PF witnesses and historical WIP change only; no full bpRNA, BPP/MFE or performance claim"}
        write_json(output / "e0_witness_observations.json", summary)
        write_json(output / "e0_witness_analysis.json", archived)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    result = check(experiment_path(args.out))
    print(json.dumps({key: result[key] for key in
                     ("passed", "new_jobs", "old_jobs", "max_absolute_delta_log_z", "max_absolute_delta_free_energy", "wip_regression")}))
    return int(not result["passed"])


if __name__ == "__main__":
    raise SystemExit(main())
