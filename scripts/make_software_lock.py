#!/usr/bin/env python3
"""Write a software lock for a dedicated clean LinearCParty/LinearKnot build.

The lock records the source commit and tree, build settings, CMake cache hash,
CLI binary hashes, energy parameter hash and three maintained checks. The
checks run from the source root, as the maintained test scripts expect.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_cli import build_settings, digest, git_state, now, write_json  # noqa: E402

PROGRAMS = ("LinearCParty", "LinearPKCapR", "LinearPKRaccess", "LinearPKSampling", "LinearPKMEA", "LinearPKCentroid")
CHECKS = ("scfg2_rules_core_exact_weight_test", "scfg2_exact_inside_oracle_test",
          "scfg2_beam_api_test")
PARAMETER = "params/rna_Turner04.par"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--software-root", type=Path, required=True)
    parser.add_argument("--build-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--note", default="", help="known follow-ups for this version")
    parser.add_argument("--required-suite", default="",
                        help="summary of the maintained required suite run on this build")
    parser.add_argument("--cli-layout", choices=("auto", "legacy", "unified"), default="auto")
    parser.add_argument("--configure-argument", action="append",
                        help="actual CMake argument, repeatable (use =-D...)")
    args = parser.parse_args()
    software, build = args.software_root.resolve(), args.build_dir.resolve()
    unified = args.cli_layout == "unified" or (args.cli_layout == "auto" and (build / "linearknot").is_file())
    programs = ("linearknot",) if unified else PROGRAMS
    state = git_state(software)
    if state["dirty"] or state["commit"] is None:
        parser.error("software root must be a clean git checkout")
    checks = []
    for name in CHECKS:
        result = subprocess.run([str(build / name)], cwd=software, capture_output=True, text=True)
        checks.append({"name": name, "returncode": result.returncode,
                       "stdout": result.stdout[-2000:], "stderr": result.stderr[-2000:]})
        if result.returncode != 0:
            parser.error("maintained check failed: " + name)
    tree = subprocess.run(["git", "-C", str(software), "rev-parse", "HEAD^{tree}"],
                          capture_output=True, text=True, check=True).stdout.strip()
    compiler = subprocess.run(["c++", "--version"], capture_output=True, text=True).stdout
    lock = {"schema_version": 1, "created_utc": now(), "purpose": "development",
            "source_repository": "https://github.com/TakumiOtagaki/" + ("LinearKnot" if unified else "LinearCParty"),
            "source_commit": state["commit"], "source_git_state": state, "source_tree": tree,
            "build_revision_verified": True,
            "verification_method": "dedicated clean checkout; CMake configure/build; selected maintained checks",
            "compiler": compiler.splitlines()[0] if compiler else None,
            "build_settings": build_settings(build),
            "configure_arguments": args.configure_argument or ["-DBUILD_TESTING=ON", "-DCMAKE_BUILD_TYPE=Release"],
            "cmake_cache_sha256": digest(build / "CMakeCache.txt"),
            "binaries": {name: digest(build / name) for name in programs},
            "energy_parameter_file": PARAMETER,
            "energy_parameter_sha256": digest(software / PARAMETER),
            "verification_results": checks,
            "required_suite": args.required_suite,
            "scientific_scope": "fixed implementation pilot, not independent proof of ensemble semantics",
            "known_semantic_followup": args.note}
    if unified:
        lock["cli_layout"] = "unified"
        lock["subcommands"] = ["fold", "profile", "access", "sample", "subopt", "mea", "centroid"]
        lock["version"] = subprocess.run([str(build / "linearknot"), "--version"],
                                         check=True, capture_output=True, text=True).stdout.strip()
    if any(value is None for value in lock["binaries"].values()):
        parser.error("a CLI binary is missing from the build")
    write_json(args.out, lock)
    print(args.out, state["commit"])


if __name__ == "__main__":
    main()
