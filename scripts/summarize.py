#!/usr/bin/env python3
"""Export measured jobs, including failures, from a CLI run manifest."""

import argparse
from collections import Counter
import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIELDS = ("id", "input", "mode", "program", "repetition", "status", "returncode",
          "elapsed_seconds", "time_wall_seconds", "peak_rss_kib", "directory")


def resolve(value):
    path = Path(value)
    return path if path.is_absolute() else ROOT / path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    run = resolve(args.run)
    manifest = json.loads((run / "run.json").read_text())
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
    summary = {"run_status": manifest["status"] if finalized else "interrupted",
               "measured_jobs": len(jobs),
               "status_counts": counts, "purpose": manifest["purpose"],
               "software_commit": manifest["software_before"]["commit"],
               "experiment_commit": manifest["experiment_before"]["commit"],
               "software_changed_during_run": manifest.get("software_changed_during_run"),
               "binaries_changed_during_run": manifest.get("binaries_changed_during_run")}
    output.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
