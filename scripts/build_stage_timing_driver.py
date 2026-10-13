#!/usr/bin/env python3
"""Link the external public-API timing observer to a verified dedicated build."""

import argparse
import json
import shlex
import subprocess
from pathlib import Path

from run_cli import ROOT, build_settings, digest, experiment_path, git_state, now, verify_lock, write_json

LIBRARIES = ("libcparty_scfg2_beam.a", "libcparty_scfg2_exact.a", "libRNA.a")


def saved_path(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def build(source, library_build, output, base_lock_path):
    base = json.loads(base_lock_path.read_text())
    if base.get("build_revision_verified") is not True:
        raise ValueError("base lock must certify its dedicated build revision")
    verify_lock(base_lock_path, source, library_build,
                [{"name": n} for n in base["binaries"]], source / base["energy_parameter_file"])
    paths = [library_build / n for n in LIBRARIES]
    hashes = {p.name: digest(p) for p in paths}
    # The E0 lock independently froze exact and RNA archives for this checkpoint.
    # The beam archive is recorded here from the original verified dedicated build.
    if base["source_commit"].startswith("15aa7742"):
        e0 = json.loads((ROOT / "manifests/software-e0-15aa7742-20261011.json").read_text())
        for p in paths[1:]:
            expected = next(e["sha256"] for e in e0["extra_files"] if Path(e["path"]).name == p.name)
            if hashes[p.name] != expected:
                raise ValueError("fixed library differs from E0 lock: " + p.name)
    driver_source = ROOT / "tools/stage_timing_driver.cc"
    output.mkdir(parents=True, exist_ok=False)
    command = ["c++", "-O2", "-std=c++17", "-Wall", "-Wextra", "-pedantic",
               '-DLCP_SOURCE_REVISION="' + base["source_commit"] + '"',
               "-I" + str(source / "src"), str(driver_source), *map(str, paths),
               "-o", str(output / "stage_timing_driver")]
    with (output / "build.log").open("w") as log:
        log.write(shlex.join(command) + "\n"); log.flush()
        subprocess.run(command, stdout=log, stderr=log, check=True)
    # run_cli/plan_array use this cache and the companion software lock.
    (output / "CMakeCache.txt").write_text(
        "CMAKE_BUILD_TYPE:STRING=Release\nCMAKE_CXX_COMPILER:FILEPATH=/bin/c++\n"
        "CMAKE_CXX_FLAGS:STRING=-O2 -std=c++17 -Wall -Wextra -pedantic\n")
    compiler = subprocess.check_output(["c++", "--version"], text=True).splitlines()[0]
    write_json(output / "build_provenance.json", {
        "created_utc": now(), "command": command, "compiler": compiler,
        "source_git_state": git_state(source), "library_build": str(library_build),
        "libraries_sha256": hashes, "base_lock_sha256": digest(base_lock_path),
        "base_verification_method": base.get("verification_method"),
        "observation_scope": "Public API calls; only scalar inside-result caching; development checks only"})
    extras = [driver_source, Path(__file__), base_lock_path, *paths,
              library_build / "CMakeCache.txt", output / "build_provenance.json"]
    extras += [library_build / n for n in base["binaries"]]
    if base["source_commit"].startswith("15aa7742"):
        extras.append(ROOT / "manifests/software-e0-15aa7742-20261011.json")
    # Guard against a concurrent source/build change while linking.
    verify_lock(base_lock_path, source, library_build,
                [{"name": n} for n in base["binaries"]], source / base["energy_parameter_file"])
    if hashes != {p.name: digest(p) for p in paths}:
        raise ValueError("libraries changed during build")
    return {
        "schema_version": 1, "created_utc": now(), "purpose": "development",
        "source_commit": base["source_commit"], "source_git_state": git_state(source),
        "build_revision_verified": True,
        "verification_method": "External driver linked to base-lock verified dedicated build; archive hashes saved; exact/RNA also checked against E0 for 15aa7742",
        "build_settings": build_settings(output), "cmake_cache_sha256": digest(output / "CMakeCache.txt"),
        "compiler": compiler, "energy_parameter_file": base["energy_parameter_file"],
        "energy_parameter_sha256": base["energy_parameter_sha256"],
        "binaries": {"stage_timing_driver": digest(output / "stage_timing_driver")},
        "extra_files": [{"path": saved_path(p), "sha256": digest(p)} for p in extras],
        "scientific_scope": "Development observer, not E2 measurements; recollect on the prefix-ranking release"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--software-root", default="deps/LinearCParty-15aa7742")
    parser.add_argument("--build-dir", default="build/linear-cparty-15aa7742")
    parser.add_argument("--base-lock", default="manifests/software-pilot-15aa7742.json")
    parser.add_argument("--out", default="build/stage-timing-15aa7742")
    parser.add_argument("--lock", default="manifests/software-stage-timing-15aa7742.json")
    args = parser.parse_args()
    lock = experiment_path(args.lock)
    if lock.exists():
        parser.error("refuse to overwrite existing software lock; choose a new --lock")
    result = build(*[experiment_path(v) for v in
                     (args.software_root, args.build_dir, args.out, args.base_lock)])
    write_json(lock, result)
    print("Built stage_timing_driver; lock: " + str(lock))


if __name__ == "__main__":
    main()
