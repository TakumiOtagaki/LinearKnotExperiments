#!/usr/bin/env python3
"""Compare three release CLIs using local scratch and one durable JSON report."""

import argparse
from itertools import zip_longest
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import sys

from run_cli import digest, experiment_path, now, run_one, write_json
from rna_records import read_record
from task_storage import local_workspace

PROGRAMS = ("LinearCParty", "LinearPKCapR", "LinearPKRaccess")
NUMBER = re.compile(rb"(?<![A-Za-z_])[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


def normalize(data):
    # Keep every other byte, including whitespace and the final newline.
    return b"".join(line for line in data.splitlines(keepends=True)
                    if not line.startswith(b"# elapsed_seconds"))


def compare_bytes(left, right, max_examples=3):
    left, right = normalize(left), normalize(right)
    maximum, comparable, examples = 0.0, True, []
    different = 0
    for number, (a, b) in enumerate(zip_longest(left.splitlines(keepends=True),
                                              right.splitlines(keepends=True)), 1):
        if a == b:
            continue
        different += 1
        if len(examples) < max_examples:
            examples.append({"line": number, "reference": a.decode("utf-8", "replace")[:500] if a is not None else None,
                             "candidate": b.decode("utf-8", "replace")[:500] if b is not None else None})
        if a is None or b is None or NUMBER.sub(b"<number>", a) != NUMBER.sub(b"<number>", b):
            comparable = False
            continue
        for x, y in zip(NUMBER.findall(a), NUMBER.findall(b)):
            delta = abs(float(x) - float(y))
            if not math.isfinite(delta):
                comparable = False
            else:
                maximum = max(maximum, delta)
    return {"identical": left == right, "max_absolute_numeric_difference": maximum,
            "all_differences_numeric": comparable, "different_lines": different,
            "different_line_examples": examples,
            "reference_normalized_sha256": hashlib.sha256(left).hexdigest(),
            "candidate_normalized_sha256": hashlib.sha256(right).hexdigest()}


def specification(candidate, beam_reference, exact_reference, inputs, parameter, beams, candidate_exact=None):
    builds = {"candidate_beam": candidate, "candidate_exact": candidate_exact or candidate,
              "reference_beam": beam_reference, "reference_exact": exact_reference}
    if not beams or any(type(b) is not int or b < 1 for b in beams) or len(set(beams)) != len(beams):
        raise ValueError("beams must be unique positive integers")
    if not inputs or len(set(map(str, inputs))) != len(inputs):
        raise ValueError("supply distinct input files")
    records = []
    for path in inputs:
        name, sequence, _ = read_record(path)
        records.append({"path": str(path), "name": name, "length": len(sequence), "sha256": digest(path)})
    if not parameter.is_file():
        raise ValueError("energy parameter file is missing")
    return {"builds": {k: {"path": str(v), "binaries": {n: digest(v / n) for n in PROGRAMS}}
                        for k, v in builds.items()}, "inputs": records,
            "parameter": {"path": str(parameter), "sha256": digest(parameter)}, "beams": beams,
            "beam_settings": "h=d=max(b,30), ranking=outer", "exact_args": ["--exact"],
            "ignored_lines": ["# elapsed_seconds"], "threads": 1,
            "numeric_difference_scope": "Matched lines with the same nonnumeric tokens; inspect all_differences_numeric and line examples for unmatched text/rows."}


def compare(spec, timeout=120, scratch_dir=None, dry_run=False):
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError("timeout must be finite and positive")
    report = {"schema_version": 1, "created_utc": now(), "specification": spec,
              "timeout_seconds": timeout, "dry_run": dry_run,
              "planned_processes": len(spec["inputs"]) * (len(spec["beams"]) + 1) * len(PROGRAMS) * 2}
    if dry_run:
        return report
    if any(h is None for b in spec["builds"].values() for h in b["binaries"].values()):
        raise ValueError("a comparison binary is missing")
    env = {**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "LC_ALL": "C"}
    rows = []
    with local_workspace(scratch_dir or os.environ.get("TMPDIR") or "/tmp") as work:
        parameter = work / "energy.par"
        parameter.write_bytes(Path(spec["parameter"]["path"]).read_bytes())
        if digest(parameter) != spec["parameter"]["sha256"]:
            raise ValueError("parameter changed during comparison")
        env["CPARTY_FIXED_ENERGY_PARAM_FILE"] = str(parameter)
        for item in spec["inputs"]:
            source = work / "input.txt"
            source.write_bytes(Path(item["path"]).read_bytes())
            if digest(source) != item["sha256"]:
                raise ValueError("input changed during comparison")
            for beam in [*spec["beams"], None]:
                kind = "exact" if beam is None else "beam"
                mode_args = ["--exact"] if beam is None else ["-b", str(beam), "--max-hairpin", str(max(beam, 30)),
                              "--max-case3-depth", str(max(beam, 30)), "--ranking", "outer"]
                for program in PROGRAMS:
                    row = {"input": item["name"], "input_sha256": item["sha256"], "mode": kind, "beam": beam, "program": program}
                    outputs = {}
                    for side in ("reference", "candidate"):
                        directory = work / side
                        directory.mkdir()
                        binary = Path(spec["builds"][side + "_" + kind]["path"]) / program
                        args = [str(binary), "-i", str(source), *mode_args]
                        if program == "LinearCParty":
                            args += ["--bpp-cutoff", "0", "--bpp", str(directory / "bpp.tsv")]
                        elif program == "LinearPKRaccess":
                            args += ["--access-len", "1,4,8,16", "--probability"]
                        result = run_one(args, directory, env, timeout)
                        row[side] = {**result, "argv": args}
                        names = ["stdout.txt", "stderr.txt"] + (["bpp.tsv"] if program == "LinearCParty" else [])
                        outputs[side] = {n: (directory / n).read_bytes() if (directory / n).is_file() else None for n in names}
                        shutil.rmtree(directory)
                    files = {}
                    for name in outputs["reference"]:
                        a, b = outputs["reference"][name], outputs["candidate"][name]
                        files[name] = compare_bytes(a, b) if a is not None and b is not None else {
                            "identical": False, "missing_reference": a is None, "missing_candidate": b is None,
                            "max_absolute_numeric_difference": 0, "all_differences_numeric": False}
                    row["files"] = files
                    row["successful"] = all(row[s]["status"] == "success" for s in outputs)
                    row["identical"] = row["successful"] and all(f["identical"] for f in files.values())
                    rows.append(row)
    # A changing binary must not be accepted as a comparison of pinned releases.
    for build in spec["builds"].values():
        if build["binaries"] != {n: digest(Path(build["path"]) / n) for n in PROGRAMS}:
            raise ValueError("comparison binaries changed during execution")
    report.update(rows=rows, successful=all(r["successful"] for r in rows),
                  identical=all(r["identical"] for r in rows),
                  identical_comparisons=sum(r["identical"] for r in rows), comparisons=len(rows),
                  max_absolute_numeric_difference=max(f["max_absolute_numeric_difference"] for r in rows for f in r["files"].values()),
                  all_differences_numeric=all(f["all_differences_numeric"] for r in rows for f in r["files"].values()))
    return report


def add_arguments(parser):
    parser.add_argument("--beam-reference-build", type=Path, default=experiment_path("build/linear-cparty-d21618c5"))
    parser.add_argument("--exact-reference-build", type=Path, default=experiment_path("build/linear-cparty-d9804df6"))
    parser.add_argument("--candidate-exact-build", type=Path, help="optional exact candidate (rehearsal with split adopted releases)")
    parser.add_argument("--input", type=Path, action="append", required=True)
    parser.add_argument("--beams", default="20,200")
    parser.add_argument("--timeout", type=float, default=120)
    parser.add_argument("--scratch-dir", type=Path)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--estimate-only", action="store_true")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-build", type=Path, required=True)
    parser.add_argument("--parameter-file", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    add_arguments(parser)
    args = parser.parse_args()
    if args.estimate_only and not args.dry_run:
        parser.error("--estimate-only requires --dry-run")
    spec = specification(args.candidate_build.resolve(), args.beam_reference_build.resolve(), args.exact_reference_build.resolve(),
                         [p.resolve() for p in args.input], args.parameter_file.resolve(), list(map(int, args.beams.split(","))),
                         args.candidate_exact_build.resolve() if args.candidate_exact_build else None)
    if args.out.exists():
        parser.error("comparison report already exists; choose a new --out")
    report = compare(spec, args.timeout, args.scratch_dir, args.dry_run)
    write_json(args.out, report)
    print(json.dumps({k: v for k, v in report.items() if k not in ("rows", "specification")}))
    return 0 if report.get("successful", True) else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, OSError) as exc:
        print("error: " + str(exc), file=sys.stderr)
        sys.exit(2)
