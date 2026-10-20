#!/usr/bin/env python3
"""Audit strict agreement and the two explicit access metadata renames."""
import argparse
from collections import Counter
import json
from pathlib import Path
import tarfile

from compare_release_outputs import normalize
from run_cli import digest, now, write_json

RENAMES = {b"# linearpkraccess_metadata=1\n": b"# linearknot_access_metadata=1\n",
           b"# engine=linearcparty-scfg2\n": b"# engine=linearknot-scfg2\n"}


def audit(directory):
    report = json.loads((directory / "comparison.json").read_text())
    assert report["status"] == "complete" and report["successful"]
    for name, expected in report["artifacts_sha256"].items():
        assert digest(directory / name) == expected
    counts, differences, other = {}, Counter(), []
    with tarfile.open(directory / "raw_outputs.tar") as archive:
        for row in report["rows"]:
            item = counts.setdefault(row["command"], dict(comparisons=0, successful=0, byte_identical=0, identical_after_metadata_rename=0))
            item["comparisons"] += 1
            item["successful"] += row["successful"]
            item["byte_identical"] += row["identical"]
            same = True
            for name, comparison in row["files"].items():
                prefix = "/".join((row["input"], row["mode"], row["command"]))
                left = normalize(archive.extractfile(prefix + "/reference/" + name).read())
                right = normalize(archive.extractfile(prefix + "/candidate/" + name).read())
                changed = b"".join(RENAMES.get(line, line) for line in left.splitlines(keepends=True))
                if left != right:
                    for line in left.splitlines(keepends=True):
                        if line in RENAMES:
                            differences[line.decode().strip()] += 1
                if changed != right:
                    same = False
                    other.append({"input": row["input"], "mode": row["mode"], "command": row["command"], "file": name})
            item["identical_after_metadata_rename"] += same and row["successful"]
    return dict(created_utc=now(), comparison_report_sha256=digest(directory / "comparison.json"),
        comparison_directory=str(directory), inputs=len(report["inputs"]), comparisons=report["comparisons"],
        processes=2*report["comparisons"], byte_identical=report["identical_comparisons"],
        total_process_attempts=report.get("total_process_attempts", 2*report["comparisons"]),
        retried_comparisons=report.get("retried_comparisons", 0),
        initial_report=report.get("initial_report"),
        by_command=counts, metadata_rename_occurrences=dict(differences), other_differences=other,
        all_identical_after_metadata_rename=not other, ignored_elapsed_line_prefix="# elapsed_seconds",
        metadata_renames={k.decode().strip(): v.decode().strip() for k, v in RENAMES.items()},
        fingerprints=report["fingerprints"], artifacts_sha256=report["artifacts_sha256"])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    result = audit(args.directory.resolve())
    write_json(args.out, result)
    print(json.dumps({k: result[k] for k in ("comparisons", "byte_identical", "all_identical_after_metadata_rename")}))


if __name__ == "__main__":
    main()
