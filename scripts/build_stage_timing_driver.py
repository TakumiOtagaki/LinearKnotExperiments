#!/usr/bin/env python3
"""Link the external public-API timing observer to a verified dedicated build."""

import argparse
import json
import shlex
import subprocess
from pathlib import Path

from run_cli import ROOT, build_settings, digest, experiment_path, git_state, now, verify_lock, write_json

LIBRARIES = ("libcparty_scfg2_beam.a", "libcparty_scfg2_exact.a", "libRNA.a")
ARCHIVE_LOCKS = {
    "15aa7742": "manifests/software-e0-15aa7742-20261011.json",
    "d21618c5": "manifests/software-e0-d21618c5-20261014.json",
}


def saved_path(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def build(source, library_build, output, base_lock_path, archive_lock_path=None, purpose="development"):
    base = json.loads(base_lock_path.read_text())
    if base.get("build_revision_verified") is not True:
        raise ValueError("base lock must certify its dedicated build revision")
    verify_lock(base_lock_path, source, library_build,
                [{"name": n} for n in base["binaries"]], source / base["energy_parameter_file"])
    paths = [library_build / n for n in LIBRARIES]
    hashes = {p.name: digest(p) for p in paths}
    # The E0 lock independently froze exact and RNA archives for this checkpoint.
    # The beam archive is recorded here from the original verified dedicated build.
    if archive_lock_path is None:
        archive_lock_path = next((ROOT / path for revision, path in ARCHIVE_LOCKS.items()
                                  if base["source_commit"].startswith(revision)), None)
    if archive_lock_path is not None:
        e0 = json.loads(archive_lock_path.read_text())
        if e0["source_commit"] != base["source_commit"] or e0.get("build_revision_verified") is not True:
            raise ValueError("archive lock must certify the same source revision")
        for p in paths[1:]:
            expected = next(e["sha256"] for e in e0["extra_files"] if Path(e["path"]).name == p.name)
            if hashes[p.name] != expected:
                raise ValueError("fixed library differs from E0 lock: " + p.name)
    frozen_files = [ROOT / "tools/stage_timing_driver.cc", Path(__file__), base_lock_path]
    if archive_lock_path is not None:
        frozen_files.append(archive_lock_path)
    frozen_hashes = {p: digest(p) for p in frozen_files}
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
        "archive_lock_sha256": digest(archive_lock_path) if archive_lock_path else None,
        "observation_scope": "Whole public API calls; inside result and (where supported) ranking guide cached; no graph handoff"})
    extras = [*frozen_files, *paths,
              library_build / "CMakeCache.txt", output / "build_provenance.json"]
    extras += [library_build / n for n in base["binaries"]]
    # Guard against a concurrent source/build change while linking.
    verify_lock(base_lock_path, source, library_build,
                [{"name": n} for n in base["binaries"]], source / base["energy_parameter_file"])
    if hashes != {p.name: digest(p) for p in paths}:
        raise ValueError("libraries changed during build")
    if frozen_hashes != {p: digest(p) for p in frozen_files}:
        raise ValueError("driver source, build script or locks changed during build")
    return {
        "schema_version": 1, "created_utc": now(), "purpose": purpose,
        "source_commit": base["source_commit"], "source_git_state": git_state(source),
        "build_revision_verified": True,
        "verification_method": "External driver linked to base-lock verified dedicated build; archive hashes saved" +
                               ("; exact/RNA checked against same-revision E0 lock" if archive_lock_path else ""),
        "build_settings": build_settings(output), "cmake_cache_sha256": digest(output / "CMakeCache.txt"),
        "compiler": compiler, "energy_parameter_file": base["energy_parameter_file"],
        "energy_parameter_sha256": base["energy_parameter_sha256"],
        "binaries": {"stage_timing_driver": digest(output / "stage_timing_driver")},
        "extra_files": [{"path": saved_path(p), "sha256": digest(p)} for p in extras],
        "scientific_scope": "Isolated whole public API calls; work counters describe main DP activity, not allocated bytes or a complexity proof"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--software-root", default="deps/LinearCParty-15aa7742")
    parser.add_argument("--build-dir", default="build/linear-cparty-15aa7742")
    parser.add_argument("--base-lock", default="manifests/software-pilot-15aa7742.json")
    parser.add_argument("--out", default="build/stage-timing-15aa7742")
    parser.add_argument("--lock", default="manifests/software-stage-timing-15aa7742.json")
    parser.add_argument("--archive-lock", help="E0 exact/RNA archive lock; inferred for 15aa7742/d21618c5")
    parser.add_argument("--purpose", choices=("development", "formal"), default="development")
    args = parser.parse_args()
    lock = experiment_path(args.lock)
    if lock.exists():
        parser.error("refuse to overwrite existing software lock; choose a new --lock")
    result = build(*[experiment_path(v) for v in
                     (args.software_root, args.build_dir, args.out, args.base_lock)],
                   archive_lock_path=experiment_path(args.archive_lock) if args.archive_lock else None,
                   purpose=args.purpose)
    write_json(lock, result)
    print("Built stage_timing_driver; lock: " + str(lock))


if __name__ == "__main__":
    main()
