#!/usr/bin/env python3
"""Build inside-only E0 observers without modifying either implementation."""

import argparse
import importlib.util
import json
import shlex
import subprocess
from pathlib import Path

from run_cli import ROOT, build_settings, digest, experiment_path, git_state, now, verify_lock, write_json

LINEAR_COMMIT = "89baea52c9435b18829ef0867a66a24be56d70f9"
ORACLE_OWNER_COMMIT = LINEAR_COMMIT
ORACLE_COMMIT = "800e95af9d1a97555233e2c7b37d73dbaa1a92b4"
PATCH_NUMBERS = ["0001", "0002", "0006", "0007", "0008", "0009"]
OBJECT_NAMES = ("part_func", "dot_plot", "mea", "centroid", "W_final", "pseudo_loop",
                "s_energy_matrix", "Hotspot", "sparse_tree", "SHAPE")


def verify_oracle(owner, oracle):
    state = git_state(owner)
    if state["dirty"] or state["commit"] != ORACLE_OWNER_COMMIT:
        raise ValueError("use the clean fixed LinearCParty-89baea52 oracle owner")
    spec = importlib.util.spec_from_file_location("e0_oracle_config", owner / "test/oracle/legacy_cparty/oracle_config.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    manifest = module.verify_checkout(oracle / "src")
    if (manifest["upstream_commit"] != ORACLE_COMMIT or
            [p["name"].split("-", 1)[0] for p in manifest["patches"]] != PATCH_NUMBERS):
        raise ValueError("E0 requires the single six-patch oracle without closure")
    return manifest


def oracle_lock(owner, oracle):
    manifest = verify_oracle(owner, oracle)
    patchdir = owner / "test/oracle/legacy_cparty"
    return {
        "schema_version": 1, "created_utc": now(), "purpose": "oracle",
        "source_repository": manifest["upstream_url"], "source_commit": ORACLE_COMMIT,
        "source_git_state": git_state(oracle / "src"),
        "source_modifications": "Exactly ordered patches 0001,0002,0006,0007,0008,0009; completed WMB weak closure absent",
        "build_revision_verified": True, "oracle_owner_commit": ORACLE_OWNER_COMMIT,
        "verification_method": "Clean pinned owner build.sh; oracle_config.verify_checkout verifies commit, ordered patches, patched sources and binary/object hashes",
        "build_command": ["bash", str((patchdir / "build.sh").relative_to(ROOT)), str(oracle)],
        "build_script_sha256": digest(patchdir / "build.sh"), "build_log_sha256": digest(oracle / "build.log"),
        "build_settings": build_settings(oracle / "src/build"),
        "cmake_cache_sha256": digest(oracle / "src/build/CMakeCache.txt"),
        "binaries": {"CParty": digest(oracle / "src/build/CParty"), "cparty-oracle": digest(oracle / "cparty-oracle")},
        "energy_parameter_file": "params/rna_Turner04.par",
        "energy_parameter_sha256": digest(owner / "params/rna_Turner04.par"),
        "patches": [{"path": str((patchdir / "patches" / p["name"]).relative_to(ROOT)), "sha256": p["sha256"]}
                    for p in manifest["patches"]],
        "source_files_sha256": {p.name: digest(p) for p in sorted((oracle / "src/src").glob("*")) if p.is_file()},
        "manifest_path": str((oracle / "manifest.json").relative_to(ROOT)),
        "manifest_sha256": digest(oracle / "manifest.json"), "oracle_manifest": manifest,
        "build_exit_status": 0,
        "scientific_scope": "Exact inside log Z / PF free energy only; BPP and MFE parity excluded"}


def build(source, library_build, oracle, output, oracle_owner=None):
    oracle_owner = oracle_owner or ROOT / "deps/LinearCParty-89baea52"
    state = git_state(source)
    if state["dirty"] or state["commit"] != LINEAR_COMMIT:
        raise ValueError("use the clean fixed LinearCParty-89baea52 checkout")
    manifest = verify_oracle(oracle_owner, oracle)
    base_lock_path = ROOT / "manifests/software-pilot-89baea52.json"
    base_lock = json.loads(base_lock_path.read_text())
    if base_lock.get("build_revision_verified") is not True:
        raise ValueError("the fixed library build must have verified source provenance")
    verify_lock(base_lock_path, source, library_build, [{"name": n} for n in base_lock["binaries"]],
                source / "params/rna_Turner04.par")
    # The supplied pilot lock certifies the dedicated build. Freeze its exact/RNA
    # archives here and check them again after linking; do not reuse old hashes.
    libraries = {name: digest(library_build / name)
                 for name in ("libcparty_scfg2_exact.a", "libRNA.a")}
    output.mkdir(parents=True, exist_ok=False)
    commands = []
    compiler = subprocess.check_output(["c++", "--version"], text=True).splitlines()[0]

    def run(argv):
        commands.append(argv)
        with (output / "build.log").open("a") as stream:
            stream.write(shlex.join(argv) + "\n"); stream.flush()
            subprocess.run(argv, stdout=stream, stderr=stream, check=True)

    run(["c++", "-O2", "-std=c++17", "-I" + str(source / "src"), str(ROOT / "tools/e0_linear_89baea52.cc"),
         str(library_build / "libcparty_scfg2_exact.a"), str(library_build / "libRNA.a"), "-o", str(output / "e0-linear")])
    objects = oracle / "src/build/CMakeFiles/CParty.dir/src"
    run(["c++", "-O2", "-std=gnu++17", "-I" + str(oracle / "src/src"), str(ROOT / "tools/e0_legacy_scalar.cc"),
         *[str(objects / (n + ".cc.o")) for n in OBJECT_NAMES], str(oracle / "src/build/libRNA.a"),
         "-o", str(output / "E0CPartyOracle")])
    for name, policy in (("E0LegacyCompatible", "legacy"), ("E0Default", "default")):
        wrapper = output / name
        wrapper.write_text('#!/usr/bin/env bash\nset -euo pipefail\nexec ' +
                           shlex.quote(str(output / "e0-linear")) + ' "$@" --policy ' + policy + "\n")
        wrapper.chmod(0o755)
    (output / "CMakeCache.txt").write_text("CMAKE_BUILD_TYPE:STRING=Release\nCMAKE_CXX_COMPILER:FILEPATH=/bin/c++\n"
                                         "CMAKE_CXX_FLAGS:STRING=-O2 -std=c++17 (linear); -O2 -std=gnu++17 (oracle)\n")
    write_json(output / "build_provenance.json", {"created_utc": now(), "compiler": compiler, "commands": commands,
        "source_commit": state["commit"], "oracle_owner_commit": ORACLE_OWNER_COMMIT,
        "library_sha256": libraries, "library_source_lock_sha256": digest(base_lock_path),
        "oracle_manifest_sha256": digest(oracle / "manifest.json"), "library_build": str(library_build),
        "oracle_build": str(oracle), "observation_scope": "inside only; original oracle objects; no AD/outside/MFE/profiles"})
    extras = [ROOT / "tools/e0_linear_89baea52.cc", ROOT / "tools/e0_legacy_scalar.cc", Path(__file__), base_lock_path,
              library_build / "CMakeCache.txt", library_build / "libcparty_scfg2_exact.a", library_build / "libRNA.a",
              output / "e0-linear", output / "build_provenance.json", oracle / "manifest.json", oracle / "cparty-oracle",
              oracle / "src/build/CParty", oracle / "src/build/CMakeCache.txt", oracle / "src/build/libRNA.a"]
    extras += [objects / (n + ".cc.o") for n in OBJECT_NAMES]
    configdir = oracle_owner / "test/oracle/legacy_cparty"
    extras += [configdir / n for n in ("config.json", "build.sh", "oracle_config.py")]
    extras += [configdir / "patches" / p["name"] for p in manifest["patches"]]
    extras += [p for p in sorted((oracle / "src/src").glob("*")) if p.is_file()]
    for name, expected in libraries.items():
        if digest(library_build / name) != expected:
            raise ValueError("fixed library changed during observer build: " + name)
    verify_lock(base_lock_path, source, library_build, [{"name": n} for n in base_lock["binaries"]],
                source / "params/rna_Turner04.par")
    verify_oracle(oracle_owner, oracle)
    return {"schema_version": 1, "created_utc": now(), "purpose": "formal", "source_commit": state["commit"],
        "source_git_state": state, "build_revision_verified": True,
        "verification_method": "Inside-only observer linked to dedicated pilot-lock-verified 89baea52 exact/RNA libraries, hashed before/after linking, and manifest-verified six-patch oracle objects",
        "compiler": compiler, "build_settings": build_settings(output), "cmake_cache_sha256": digest(output / "CMakeCache.txt"),
        "energy_parameter_file": "params/rna_Turner04.par", "energy_parameter_sha256": digest(source / "params/rna_Turner04.par"),
        "binaries": {n: digest(output / n) for n in ("E0LegacyCompatible", "E0Default", "E0CPartyOracle")},
        "extra_files": [{"path": str(p.relative_to(ROOT)), "sha256": digest(p)} for p in extras],
        "oracle": {"source_commit": ORACLE_COMMIT, "patches": 6, "patch_order": PATCH_NUMBERS,
                   "completed_wmb_weak_closure": False, "manifest_sha256": digest(oracle / "manifest.json")},
        "observation_scope": "log_z, free_energy, rt only; BPP and MFE parity excluded"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--software-root", default="deps/LinearCParty-89baea52")
    parser.add_argument("--build-dir", default="build/linear-cparty-89baea52")
    parser.add_argument("--oracle-owner", default="deps/LinearCParty-89baea52")
    parser.add_argument("--oracle-dir", default="build/cparty-oracle-800e95a-6patches")
    parser.add_argument("--out", default="build/e0-observers-89baea52-800e95a-6patches")
    parser.add_argument("--lock", default="manifests/software-e0-89baea52-20261018.json")
    parser.add_argument("--oracle-lock", default="manifests/software-cparty-oracle-800e95a-6patches-20261018.json")
    args = parser.parse_args()
    for destination in (args.out, args.lock, args.oracle_lock):
        if experiment_path(destination).exists():
            parser.error("use a new output/lock path: " + destination)
    owner, oracle = experiment_path(args.oracle_owner), experiment_path(args.oracle_dir)
    write_json(experiment_path(args.oracle_lock), oracle_lock(owner, oracle))
    result = build(*[experiment_path(v) for v in (args.software_root, args.build_dir, args.oracle_dir, args.out)], owner)
    write_json(experiment_path(args.lock), result)


if __name__ == "__main__":
    main()
