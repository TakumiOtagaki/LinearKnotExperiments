#!/usr/bin/env python3
"""Build an inside-only benchmark against the existing locked builds."""
import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile

from build_e0_tools_89baea52 import verify_oracle
from run_cli import ROOT, build_settings, digest, git_state, now, verify_lock, write_json


def build(output, lock_path):
    if output.exists() or lock_path.exists():
        raise ValueError("use new observer build and lock paths")
    source = ROOT / "deps/LinearCParty-89baea52"
    library = ROOT / "build/linear-cparty-89baea52"
    oracle = ROOT / "build/cparty-oracle-800e95a-6patches"
    e0build = ROOT / "build/e0-observers-89baea52-800e95a-6patches"
    base_path = ROOT / "manifests/software-stage-timing-89baea52-20261018.json"
    e0_path = ROOT / "manifests/software-e0-89baea52-20261018.json"
    base = json.loads(base_path.read_text())
    e0 = json.loads(e0_path.read_text())
    verify_lock(base_path, source, ROOT / "build/stage-timing-89baea52",
                [{"name": "stage_timing_driver"}], source / base["energy_parameter_file"])
    verify_lock(e0_path, source, e0build, [{"name": "E0CPartyOracle"}],
                source / e0["energy_parameter_file"])
    verify_oracle(source, oracle)
    archives = [library / n for n in ("libcparty_scfg2_beam.a", "libcparty_scfg2_exact.a", "libRNA.a")]
    frozen = {str(p): digest(p) for p in archives}
    for p in archives:
        expected = next(x["sha256"] for x in base["extra_files"] if x["path"] == str(p.relative_to(ROOT)))
        if frozen[str(p)] != expected:
            raise ValueError("archive changed: " + str(p))
    code = ROOT / "tools/pf_scaling_89baea52.cc"
    output.mkdir(parents=True)
    with tempfile.TemporaryDirectory(prefix="pf-scaling-build-", dir=os.environ["TMPDIR"]) as tmp:
        binary = Path(tmp) / "PFScaling"
        command = ["c++", "-O2", "-std=c++17", "-Wall", "-Wextra", "-pedantic",
                   '-DCPARTY_ORACLE_PATH="' + str(e0build / "E0CPartyOracle") + '"',
                   "-I" + str(source / "src"), str(code), *map(str, archives), "-o", str(binary)]
        process = subprocess.run(command, capture_output=True, text=True)
        (output / "build.log").write_text(shlex.join(command) + "\n" + process.stdout + process.stderr)
        process.check_returncode()
        shutil.copy2(binary, output / "PFScaling")
    (output / "CMakeCache.txt").write_text("CMAKE_BUILD_TYPE:STRING=Release\nCMAKE_CXX_COMPILER:FILEPATH=/bin/c++\nCMAKE_CXX_FLAGS:STRING=-O2 -std=c++17 -Wall -Wextra -pedantic\n")
    write_json(output / "build_provenance.json", {
        "created_utc": now(), "command": command, "library_sha256": frozen,
        "compiler": subprocess.check_output(["c++", "--version"], text=True).splitlines()[0],
        "scope": "PF only; CParty execs unchanged E0 zero-scaling observer; no sampling/outside/MFE"})
    extras = {x["path"]: x["sha256"] for x in [*base["extra_files"], *e0["extra_files"]]}
    for p in [code, Path(__file__).resolve(), base_path, e0_path, e0build / "E0CPartyOracle",
              ROOT / "manifests/software-cparty-oracle-800e95a-6patches-20261018.json",
              output / "build_provenance.json"]:
        extras[str(p.relative_to(ROOT))] = digest(p)
    lock = {"schema_version": 1, "created_utc": now(), "purpose": "formal",
            "source_commit": e0["source_commit"], "source_git_state": git_state(source),
            "build_revision_verified": True, "verification_method": "Locked E0 observer plus stage-lock verified adopted libraries; new external inside-only driver",
            "build_settings": build_settings(output), "cmake_cache_sha256": digest(output / "CMakeCache.txt"),
            "energy_parameter_file": e0["energy_parameter_file"],
            "energy_parameter_sha256": e0["energy_parameter_sha256"],
            "binaries": {"PFScaling": digest(output / "PFScaling")},
            "extra_files": [{"path": p, "sha256": h} for p, h in sorted(extras.items())],
            "scientific_scope": "Inside-only PF; legacy-compatible closure for all three methods; beam h=d=100; exact uncapped; numerical CParty failures retained"}
    if frozen != {str(p): digest(p) for p in archives}:
        raise ValueError("libraries changed during linking")
    verify_oracle(source, oracle)
    verify_lock(e0_path, source, e0build, [{"name": "E0CPartyOracle"}], source / e0["energy_parameter_file"])
    write_json(lock_path, lock)
    verify_lock(lock_path, source, output, [{"name": "PFScaling"}], source / e0["energy_parameter_file"])
    return lock


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "build/pf-scaling-89baea52-20261020")
    parser.add_argument("--lock", type=Path, default=ROOT / "manifests/software-pf-scaling-89baea52-20261020.json")
    args = parser.parse_args()
    build(args.out.resolve(), args.lock.resolve())
