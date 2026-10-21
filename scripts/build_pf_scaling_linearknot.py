#!/usr/bin/env python3
"""Link an inside-only observer to a certified LinearKnot release; never submit."""
import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile

from run_cli import ROOT, build_settings, digest, git_state, now, verify_lock, write_json


def build(source, library, base_path, stage_build, output, lock_path):
    if output.exists() or lock_path.exists():
        raise ValueError("choose fresh observer build and lock paths")
    base = json.loads(base_path.read_text())
    if not base.get("build_revision_verified"):
        raise ValueError("base lock must certify the dedicated release")
    verify_lock(base_path, source, stage_build, [{"name": "stage_timing_driver"}],
                source / base["energy_parameter_file"])
    archives = [library / n for n in ("libcparty_scfg2_beam.a", "libcparty_scfg2_exact.a", "libRNA.a")]
    frozen = {p: digest(p) for p in archives}
    for p in archives:
        expected = next(x["sha256"] for x in base["extra_files"]
                        if (ROOT / x["path"]).resolve() == p.resolve())
        if frozen[p] != expected:
            raise ValueError("library differs from stage lock: " + str(p))
    code = ROOT / "tools/pf_scaling_linearknot.cc"
    output.mkdir(parents=True)
    with tempfile.TemporaryDirectory(prefix="pf-linearknot-build-", dir=os.environ["TMPDIR"]) as tmp:
        binary = Path(tmp) / "PFScaling"
        command = ["c++", "-O2", "-std=c++17", "-Wall", "-Wextra", "-pedantic",
                   "-I" + str(source / "src"), str(code), *map(str, archives), "-o", str(binary)]
        process = subprocess.run(command, capture_output=True, text=True)
        (output / "build.log").write_text(shlex.join(command) + "\n" + process.stdout + process.stderr)
        process.check_returncode()
        shutil.copy2(binary, output / "PFScaling")
    (output / "CMakeCache.txt").write_text(
        "CMAKE_BUILD_TYPE:STRING=Release\nCMAKE_CXX_COMPILER:FILEPATH=/bin/c++\n"
        "CMAKE_CXX_FLAGS:STRING=-O2 -std=c++17 -Wall -Wextra -pedantic\n")
    write_json(output / "build_provenance.json", {
        "created_utc": now(), "command": command,
        "libraries_sha256": {str(p.relative_to(ROOT)): h for p, h in frozen.items()},
        "compiler": subprocess.check_output(["c++", "--version"], text=True).splitlines()[0],
        "source_commit": base["source_commit"], "base_lock_sha256": digest(base_path),
        "scope": "PF only; exact and beam100; LegacyCompatible; no MFE/outside/sampling or CParty execution"})
    extras = {x["path"]: x["sha256"] for x in base["extra_files"]}
    for p in (code, Path(__file__).resolve(), base_path, output / "build_provenance.json"):
        extras[str(p.relative_to(ROOT))] = digest(p)
    lock = {"schema_version": 1, "created_utc": now(), "purpose": "formal",
            "source_commit": base["source_commit"], "source_git_state": git_state(source),
            "build_revision_verified": True,
            "verification_method": "Inside-only observer linked to same-release stage-lock certified archives",
            "build_settings": build_settings(output), "cmake_cache_sha256": digest(output / "CMakeCache.txt"),
            "energy_parameter_file": base["energy_parameter_file"],
            "energy_parameter_sha256": base["energy_parameter_sha256"],
            "binaries": {"PFScaling": digest(output / "PFScaling")},
            "extra_files": [{"path": p, "sha256": h} for p, h in sorted(extras.items())],
            "scientific_scope": "Inside-only PF; exact uncapped and beam h=d=100, outer ranking; compatible closure"}
    verify_lock(base_path, source, stage_build, [{"name": "stage_timing_driver"}],
                source / base["energy_parameter_file"])
    if frozen != {p: digest(p) for p in archives}:
        raise ValueError("libraries changed during link")
    write_json(lock_path, lock)
    verify_lock(lock_path, source, output, [{"name": "PFScaling"}], source / base["energy_parameter_file"])
    return lock


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--software-root", type=Path, default=ROOT / "deps/LinearKnot-4f97fd3")
    p.add_argument("--library-build", type=Path, default=ROOT / "build/linearknot-4f97fd3")
    p.add_argument("--base-lock", type=Path, default=ROOT / "manifests/software-stage-timing-4f97fd3-20261021.json")
    p.add_argument("--stage-build", type=Path, default=ROOT / "build/stage-timing-4f97fd3")
    p.add_argument("--out", type=Path, default=ROOT / "build/pf-scaling-4f97fd3-20261021")
    p.add_argument("--lock", type=Path, default=ROOT / "manifests/software-pf-scaling-4f97fd3-20261021.json")
    a = p.parse_args()
    build(*[x.resolve() for x in (a.software_root, a.library_build, a.base_lock, a.stage_build, a.out, a.lock)])
    print(a.lock)
