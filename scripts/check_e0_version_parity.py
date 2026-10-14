#!/usr/bin/env python3
"""Small development check of exact inside on 15aa7742 and 2456f9b2."""

import argparse
import importlib.util
import json
import os
from pathlib import Path
import random
import subprocess

from analyze_e0 import analysis_settings, compare, read_observation, within_tolerance
from run_cli import ROOT, digest, experiment_path, git_state, now, run_one, verify_lock, write_json

# This script reproduces the historical 2026-10-11 comparison, not the latest
# E0 builder's version pair.
LINEAR_COMMIT = "15aa7742a5497164573294e37972ff7ccf9ddfe1"
ORACLE_OWNER_COMMIT = "2456f9b21691a14f63511df29068628a9b2e4193"


def check(output):
    old = ROOT / "deps/LinearCParty-15aa7742"
    latest = ROOT / "deps/LinearCParty-2456f9b2"
    oldbuild = ROOT / "build/e0-observers-15aa7742-800e95a-4patches"
    newbuild = ROOT / "build/linear-cparty-2456f9b2-e0"
    config = json.loads((ROOT / "configs/e0_bprna_legacy_20261011.json").read_text())
    before = [git_state(p) for p in (old, latest)]
    if any(s["dirty"] or s["commit"] != c for s, c in zip(before, (LINEAR_COMMIT, ORACLE_OWNER_COMMIT))):
        raise ValueError("clean pinned source trees required")
    verify_lock(ROOT / config["software_lock"], old, oldbuild, config["programs"], old / config["energy_parameter_file"])
    output.mkdir(parents=True, exist_ok=False)
    # The latest exact library was built separately with BUILD_TESTING=OFF.
    command = ["c++", "-O2", "-std=c++17", "-I" + str(latest / "src"), str(ROOT / "tools/e0_linear.cc"),
               str(newbuild / "libcparty_scfg2_exact.a"), str(newbuild / "libRNA.a"),
               "-o", str(output / "e0-linear-2456f9b2")]
    result = subprocess.run(command, capture_output=True, text=True)
    (output / "observer-build.stdout").write_text(result.stdout)
    (output / "observer-build.stderr").write_text(result.stderr)
    result.check_returncode()
    generator = latest / "test/tools/gen_random_seq_struct.py"
    spec = importlib.util.spec_from_file_location("e0_random_input", generator)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    cases = [(entry["name"], (ROOT / entry["path"]).read_text()) for entry in config["inputs"]]
    rng = random.Random(20261011)
    for n in (20, 30, 50, 80):
        for index in range(2):
            g = module.generate_structure(rng, n, 3); s = module.assign_sequence(rng, g)
            if index == 1:
                g = "".join(("x" if i % 2 else "X") if c == "." and i % 3 == 0 else c for i, c in enumerate(g))
            name = f"random{n}_{index}"
            cases.append((name, f">{name}\n{s}\n{g}\n"))
    cases += [("vp_split25", ">vp_split25\nGGGAAAAACGGAAAAACAAAACCAC\n(.........(.....)....)...\n"),
              ("be_boundary41", ">be_boundary41\nG" + "A" * 16 + "G" + "A" * 7 + "C" + "A" * 14 + "C\n(" + "." * 16 + "(" + "." * 7 + ")" + "." * 14 + ")\n")]
    env = {**os.environ, **config["environment"], "CPARTY_FIXED_ENERGY_PARAM_FILE": str(old / config["energy_parameter_file"])}
    rows = []
    for name, payload in cases:
        record = output / (name + ".txt"); record.write_text(payload)
        for policy in ("legacy", "default"):
            row = {"input": name, "policy": policy, "input_sha256": digest(record), "passed": False, "jobs": []}
            values = []
            for version, binary in (("15aa7742", oldbuild / "e0-linear"), ("2456f9b2", output / "e0-linear-2456f9b2")):
                directory = output / name / policy / version; directory.mkdir(parents=True)
                argv = [str(binary), "-i", str(record), "--policy", policy]
                observation = run_one(argv, directory, env, 60, 4096)
                row["jobs"].append({"version": version, "directory": str(directory.relative_to(output)), **observation})
                if observation["status"] == "success":
                    values.append(read_observation(directory / "stdout.txt"))
            if len(values) == 2:
                row.update(old=values[0], latest=values[1], metrics=compare(values[1], values[0]))
                row["passed"] = within_tolerance(row["metrics"], config["analysis"])
            rows.append(row)
    after = [git_state(p) for p in (old, latest)]
    if after != before: raise ValueError("source trees changed during check")
    artifacts = [generator, ROOT / "tools/e0_linear.cc", newbuild / "CMakeCache.txt",
                 newbuild / "libcparty_scfg2_exact.a", newbuild / "libRNA.a", oldbuild / "e0-linear",
                 output / "e0-linear-2456f9b2", newbuild / "configure.log", newbuild / "build.log"]
    changed = subprocess.check_output(["git", "-C", str(latest), "diff", "--name-only", LINEAR_COMMIT, ORACLE_OWNER_COMMIT, "--", "src"], text=True).splitlines()
    summary = {"schema_version": 1, "created_utc": now(), "purpose": "development", "inputs": len(cases),
        "comparisons": len(rows), "passed": sum(r["passed"] for r in rows), "source_before": before, "source_after": after,
        "src_changes": changed, "observer_compile_command": command,
        "latest_library_build_commands": [["cmake", "-S", str(latest), "-B", str(newbuild), "-DCMAKE_BUILD_TYPE=Release", "-DBUILD_TESTING=OFF"],
                                          ["cmake", "--build", str(newbuild), "-j2", "--target", "cparty_scfg2_exact"]],
        "parameter_sha256": digest(old / config["energy_parameter_file"]), "generator_seed": 20261011,
        "script_sha256": digest(Path(__file__)), "tolerances": analysis_settings(),
        "artifacts": [{"path": str(p.relative_to(ROOT)), "sha256": digest(p)} for p in artifacts],
        "max_absolute_delta_log_z": max((abs(r["metrics"]["delta_log_z"]) for r in rows if "metrics" in r), default=None),
        "max_absolute_delta_free_energy": max((abs(r["metrics"]["delta_free_energy"]) for r in rows if "metrics" in r), default=None),
        "scope": "Finite small-input check with the valid pinned Turner04 file; no claim about beam/suboptimal/factorized outputs or invalid parameter inputs",
        "observations": rows}
    write_json(output / "summary.json", summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--out", required=True)
    args = parser.parse_args(); result = check(experiment_path(args.out))
    print(json.dumps({k: result[k] for k in ("inputs", "comparisons", "passed", "max_absolute_delta_log_z", "max_absolute_delta_free_energy")}))
    return int(result["passed"] != result["comparisons"])


if __name__ == "__main__": raise SystemExit(main())
