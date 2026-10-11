#!/usr/bin/env python3
"""Check output domains and cross-program identities without an exact reference."""

import argparse
from collections import Counter
import json
from pathlib import Path

from analyze_accuracy import PROGRAMS, check_consistency, read_job
from rna_records import read_record
from run_cli import ROOT, digest, experiment_path, git_state, now, write_json


def validate_outputs(run):
    manifest = json.loads((run / "run.json").read_text())
    if manifest["status"] == "running":
        raise ValueError("wait for the run to finish")
    records = {}
    for entry in manifest["inputs"]:
        path = Path(entry["copied_path"])
        if digest(path) != entry["sha256"]:
            raise ValueError("copied input hash changed")
        records[entry["name"]] = read_record(path)
    groups, observations = {}, []
    jobs = [job for job in manifest["jobs"] if not job["warmup"]]
    for job in jobs:
        row = {key: job[key] for key in ("id", "input", "mode", "program", "status", "repetition")}
        row["stdout_sha256"] = digest(run / job["directory"] / "stdout.txt")
        if job["status"] != "success":
            row["validation_status"] = "job_missing" if job.get("result_present") is False else "job_failed"
        elif job["program"] not in PROGRAMS:
            row["validation_status"] = "unsupported_program"
        else:
            try:
                parsed = read_job(run, job, records[job["input"]])
                groups.setdefault((job["input"], job["mode"], job["repetition"]), {})[job["program"]] = parsed
                row["validation_status"] = "success"
                if job["program"] == "LinearCParty":
                    row["bpp_sha256"] = digest(Path(job["command"][job["command"].index("--bpp") + 1]))
                if job["program"] == "LinearPKRaccess":
                    row.update(ensemble_free_energy_kcal_mol=parsed["free_energy"], rt_kcal_mol=parsed["rt"])
            except (ValueError, KeyError, OSError) as exc:
                row.update(validation_status="artifact_error", error=str(exc))
        observations.append(row)
    checks = [check_consistency(*key, group) for key, group in groups.items() if PROGRAMS <= set(group)]
    result = {"schema_version": 1, "generated_utc": now(),
            "scope": "output domains and cross-program consistency; no exact accuracy reference",
            "experiment_git_state": git_state(ROOT), "run_manifest_sha256": digest(run / "run.json"),
            "validator_sha256": digest(Path(__file__)),
            "analysis_helpers_sha256": digest(Path(__file__).with_name("analyze_accuracy.py")),
            "record_helpers_sha256": digest(Path(__file__).with_name("rna_records.py")),
            "run_status": manifest["status"], "measured_jobs": len(jobs),
            "job_status_counts": dict(Counter(job["status"] for job in jobs)),
            "observations": observations, "consistency_checks": checks,
            "complete_groups_checked": len(checks)}
    if manifest.get("run_kind") == "array":
        result.update({key: manifest[key] for key in ("planned_tasks", "task_status_counts",
                       "planned_measured_jobs", "completed_measured_jobs", "missing_measured_jobs")})
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    output = experiment_path(args.out)
    if output.exists():
        parser.error("use a new output path")
    result = validate_outputs(experiment_path(args.run))
    output.parent.mkdir(parents=True, exist_ok=True)
    write_json(output, result)
    failures = sum(r["validation_status"] != "success" for r in result["observations"])
    failures += sum(r["status"] != "success" for r in result["consistency_checks"])
    print(json.dumps({"measured_jobs": result["measured_jobs"], "failures": failures,
                      "complete_groups_checked": result["complete_groups_checked"]}))
    return int(failures > 0)


if __name__ == "__main__":
    raise SystemExit(main())
