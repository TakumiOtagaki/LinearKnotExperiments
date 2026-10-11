#!/usr/bin/env python3
"""Export measured jobs, including failures, from a CLI run manifest."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path
from task_storage import RunReader


ROOT = Path(__file__).resolve().parents[1]
FIELDS = ("id", "input", "mode", "program", "repetition", "status", "returncode",
          "elapsed_seconds", "time_wall_seconds", "peak_rss_kib", "directory",
          "tier", "task_id", "global_task_number", "result_present", "output_archive",
          "archive_member_directory")


def resolve(value):
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    run = resolve(args.run)
    manifest = RunReader(run).manifest
    if manifest.get("summary_only"):
        raise ValueError("aggregate without --summary-only before exporting observations")
    jobs = [job for job in manifest["jobs"] if not job["warmup"]]
    output = resolve(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(jobs)
    counts = dict(Counter(job["status"] for job in jobs))
    finalized = all(key in manifest for key in (
        "finished_utc", "software_changed_during_run", "binaries_changed_during_run"))
    summary = {"run_status": manifest["status"] if finalized or manifest.get("run_kind") == "array" else "interrupted",
               "measured_jobs": len(jobs),
               "status_counts": counts, "purpose": manifest["purpose"],
               "software_commit": manifest["software_before"]["commit"],
               "experiment_commit": manifest["experiment_before"]["commit"],
               "software_changed_during_run": manifest.get("software_changed_during_run"),
               "binaries_changed_during_run": manifest.get("binaries_changed_during_run")}
    if manifest.get("run_kind") == "array":
        summary.update({key: manifest[key] for key in ("planned_tasks", "task_status_counts",
                       "planned_measured_jobs", "completed_measured_jobs", "missing_measured_jobs")})
        summary["software_lock_valid_after_run"] = manifest["software_lock_valid_after_run"]
    output.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
