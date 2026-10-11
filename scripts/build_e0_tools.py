#!/usr/bin/env python3
"""Build external E0 observers; never change the pinned implementation/oracle."""

import argparse
import json
from pathlib import Path
import shlex
import shutil
import subprocess

from run_cli import ROOT, build_settings, digest, experiment_path, git_state, now, write_json


def build_scalar(oracle, output):
    objects = oracle / "src/build/CMakeFiles/CParty.dir/src"
    names = ("part_func", "dot_plot", "mea", "centroid", "W_final", "pseudo_loop", "s_energy_matrix", "Hotspot", "sparse_tree", "SHAPE")
    command = ["c++", "-O2", "-std=gnu++17", "-I" + str(oracle / "src/src"),
               str(ROOT / "tools/e0_legacy_scalar.cc"), *[str(objects / (n + ".cc.o")) for n in names],
               str(oracle / "src/build/libRNA.a"), "-o", str(output / "e0-legacy-uninstrumented")]
    with (output / "build.log").open("a") as stream:
        stream.write(shlex.join(command) + "\n");stream.flush()
        subprocess.run(command, stdout=stream, stderr=stream, check=True)
    return command


def build(source, library_build, oracle, output):
    if git_state(source)["dirty"] or git_state(source)["commit"] != "15aa7742a5497164573294e37972ff7ccf9ddfe1":
        raise ValueError("use the clean fixed LinearCParty-15aa7742 checkout")
    if git_state(oracle / "src")["commit"] != "59593650ad5c5b2ba54bf92ee843c6b28403181e":
        raise ValueError("oracle source revision differs")
    output.mkdir(parents=True, exist_ok=False)
    observation = output / "legacy-observation"
    shutil.copytree(oracle / "src/src", observation)
    types = observation / "base_types.hh"
    original = types.read_text()
    if original.count("typedef double pf_t;") != 1:
        raise ValueError("unexpected legacy scalar definition")
    types.write_text(original.replace("typedef double pf_t;", '#include "legacy_ad.hh"\ntypedef oracle_ad::Real pf_t;'))
    pf = observation / "part_func.cc"
    original = pf.read_text()
    replacements = {"V[ij] = contributions;": "V[ij] = unpaired ? oracle_ad::tag(contributions, i, j) : contributions;",
                    "VP[ij] = contributions;": "VP[ij] = oracle_ad::tag(contributions, i, j);"}
    for before, after in replacements.items():
        if original.count(before) != 1:
            raise ValueError("unexpected oracle pair emission: " + before)
        original = original.replace(before, after)
    # Make the two scalar conditional literals explicit for the observer type.
    if original.count("(k > 1) ? W[k - 1] : 1;") != 2:
        raise ValueError("unexpected oracle scalar conditional")
    original = original.replace("(k > 1) ? W[k - 1] : 1;", "(k > 1) ? W[k - 1] : pf_t(1);")
    pf.write_text(original)
    commands = []
    compiler = subprocess.check_output(["c++", "--version"], text=True).splitlines()[0]

    def run(argv):
        commands.append(argv)
        with (output / "build.log").open("a") as stream:
            stream.write(shlex.join(argv) + "\n")
            stream.flush()
            subprocess.run(argv, stdout=stream, stderr=stream, check=True)

    run(["c++", "-O2", "-std=c++17", "-I" + str(source / "src"), str(ROOT / "tools/e0_linear.cc"),
         str(library_build / "libcparty_scfg2_beam.a"), str(library_build / "libcparty_scfg2_exact.a"),
         str(library_build / "libRNA.a"), "-o", str(output / "e0-linear")])
    objects = []
    for name in ("part_func", "dot_plot", "mea", "centroid"):
        obj = output / (name + ".o")
        run(["c++", "-O1", "-std=gnu++17", "-I" + str(observation), "-I" + str(ROOT / "tools"),
             "-c", str(observation / (name + ".cc")), "-o", str(obj)])
        objects.append(str(obj))
    original_objects = oracle / "src/build/CMakeFiles/CParty.dir/src"
    for name in ("W_final", "pseudo_loop", "s_energy_matrix", "Hotspot", "sparse_tree", "SHAPE"):
        objects.append(str(original_objects / (name + ".cc.o")))
    run(["c++", "-O1", "-std=gnu++17", "-I" + str(observation), "-I" + str(ROOT / "tools"),
         str(ROOT / "tools/e0_legacy.cc"), *objects, str(oracle / "src/build/libRNA.a"),
         "-o", str(output / "E0CPartyOracle")])
    commands.append(build_scalar(oracle, output))
    for name, policy in (("E0LegacyCompatible", "legacy"), ("E0Default", "default")):
        wrapper = output / name
        wrapper.write_text('#!/usr/bin/env bash\nset -euo pipefail\nexec ' +
                           shlex.quote(str(output / "e0-linear")) + ' "$@" --policy ' + policy + "\n")
        wrapper.chmod(0o755)
    # A real external build record, rather than relabeling the implementation cache.
    (output / "CMakeCache.txt").write_text("CMAKE_BUILD_TYPE:STRING=Release\nCMAKE_CXX_COMPILER:FILEPATH=/bin/c++\n"
                                           "CMAKE_CXX_FLAGS:STRING=-O2 -std=c++17 (linear); -O1 -std=gnu++17 (AD oracle)\n")
    write_json(output / "build_provenance.json", {"created_utc": now(), "compiler": compiler,
               "commands": commands, "source_commit": git_state(source)["commit"],
               "instrumentation": replacements, "oracle_scalar": "external copy uses reverse AD pf_t",
               "library_build": str(library_build), "oracle_build": str(oracle)})
    extras = [ROOT / "tools/e0_linear.cc", ROOT / "tools/e0_legacy.cc", ROOT / "tools/legacy_ad.hh",
              ROOT / "tools/e0_legacy_scalar.cc", output / "e0-legacy-uninstrumented",
              Path(__file__), library_build / "CMakeCache.txt", library_build / "libcparty_scfg2_beam.a",
              library_build / "libcparty_scfg2_exact.a", library_build / "libRNA.a", output / "e0-linear",
              oracle / "cparty-oracle", oracle / "src/build/CParty", oracle / "src/build/CMakeCache.txt",
              oracle / "src/build/libRNA.a", output / "build_provenance.json"]
    extras += list((source / "test/oracle/legacy_cparty/patches").glob("*.patch"))
    extras += sorted(observation.glob("*"))
    extras += sorted((oracle / "src/src").glob("*"))
    extras = [p for p in extras if p.is_file()]
    lock = {"schema_version": 1, "created_utc": now(), "purpose": "formal",
            "source_commit": git_state(source)["commit"], "source_git_state": git_state(source),
            "build_revision_verified": True, "verification_method": "external observer linked to locked libraries; seven-patch oracle observation build",
            "compiler": compiler, "build_settings": build_settings(output),
            "cmake_cache_sha256": digest(output / "CMakeCache.txt"),
            "energy_parameter_file": "params/rna_Turner04.par",
            "energy_parameter_sha256": digest(source / "params/rna_Turner04.par"),
            "binaries": {name: digest(output / name) for name in ("E0LegacyCompatible", "E0Default", "E0CPartyOracle")},
            "extra_files": [{"path": str(p.relative_to(ROOT)), "sha256": digest(p)} for p in extras],
            "oracle": {"source_commit": "59593650ad5c5b2ba54bf92ee843c6b28403181e", "patches": 7,
                       "completed_wmb_weak_closure": False},
            "profile_scope": "Linear exact profiles in both policies; upstream has sampling profiles only, no independent analytic profile oracle"}
    return lock


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--software-root", default="deps/LinearCParty-15aa7742")
    parser.add_argument("--build-dir", default="build/linear-cparty-15aa7742")
    parser.add_argument("--oracle-dir", default="build/legacy-cparty-59593650-7patches")
    parser.add_argument("--out", default="build/e0-observers-15aa7742-v2")
    parser.add_argument("--lock", default="manifests/software-e0-15aa7742-20261011.json")
    args = parser.parse_args()
    result = build(*[experiment_path(v) for v in (args.software_root, args.build_dir, args.oracle_dir, args.out)])
    write_json(experiment_path(args.lock), result)


if __name__ == "__main__":
    main()
