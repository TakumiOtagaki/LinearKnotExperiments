#!/usr/bin/env python3
"""Validate RNA-level E1 metrics before using them for stratified analysis.

This reads only the new aggregate and its provenance. It does not open raw
archives or execute RNA software. Metric identities use the existing E1 audit.
"""

import argparse
from collections import Counter
import csv
import json

from report_e1_production import validate_row
from run_cli import digest


def main():
    from pathlib import Path
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    directory = args.directory
    provenance = json.loads((directory / "provenance.json").read_text())
    checksum = digest(directory / "metrics.tsv")
    assert checksum == provenance["metrics_sha256"], "metrics checksum changed"
    assert set(provenance["selected_job_status_counts"]) == {"success"}
    seen, counts = set(), Counter()
    with (directory / "metrics.tsv").open() as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            validate_row(row)
            key = row["input"], int(row["beam"])
            assert key not in seen, "duplicate input/beam"
            seen.add(key)
            counts[row["tier"], int(row["beam"])] += 1
            assert row["reference_mode"] == provenance["reference_mode"]
            for field in ("delta_free_energy_kcal_mol", "base_state_tv_mean_Gfree", "context_tv_mean_Gfree"):
                assert row[field] != "NA", "missing core metric"
    assert len(seen) == provenance["metric_rows"]
    assert len({key[0] for key in seen}) == provenance["inputs"]
    assert {beam for _, beam in seen} == set(provenance["beams"])
    for tier, _, _ in provenance["bins"]:
        assert len({counts[tier, beam] for beam in provenance["beams"]}) == 1
    summary = json.loads((directory / "summary.json").read_text())
    assert summary["metrics_sha256"] == checksum
    report = dict(all_passed=True, metrics_sha256=checksum, rows_validated=len(seen),
                  provenance_sha256=digest(directory / "provenance.json"),
                  script_sha256=digest(Path(__file__)),
                  checks="E1 scalar identities, finite core metrics, unique complete input/beam matrix, window-domain partitions, NA domains and hashes")
    (directory / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report))


if __name__ == "__main__":
    main()
