#!/usr/bin/env python3
"""Adapt the existing membership observer to the packed task CLI runner."""
import argparse
import json
import os
from pathlib import Path
import sys

from analyze_pk_membership import execute, prepare
from rna_records import read_record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("-i", type=Path, required=True)
    p.add_argument("--driver", type=Path, required=True)
    p.add_argument("--cohort", type=Path, required=True)
    args = p.parse_args()
    name, sequence, scaffold = read_record(args.i)
    cohort = {r["id"]: r for r in json.loads(args.cohort.read_text())}
    prepared = prepare(cohort[name])
    if sequence != cohort[name]["seq"] or scaffold != prepared["scaffold"]:
        raise ValueError("frozen membership input differs from the original cohort")
    # The enclosing run_cli process-group timeout/cap governs the observer too.
    result = execute(prepared, args.driver,
                     Path(os.environ["CPARTY_FIXED_ENERGY_PARAM_FILE"]), 86400)
    result["probability_source_revision"] = "89baea52c9435b18829ef0867a66a24be56d70f9"
    print(json.dumps(result, allow_nan=False))
    return 0 if result.get("returncode") == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
