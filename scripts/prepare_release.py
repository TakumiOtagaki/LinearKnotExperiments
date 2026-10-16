#!/usr/bin/env python3
"""Prepare a dedicated release; inspect existing stages without rebuilding them."""

import argparse
from contextlib import contextmanager
import json
import os
from pathlib import Path
import subprocess
import sys

from compare_release_outputs import add_arguments, compare, specification
from file_budget import check_budget
from make_software_lock import CHECKS, PARAMETER, PROGRAMS
from run_cli import ROOT, build_settings, digest, git_state, now, verify_lock, write_json
from task_storage import local_workspace


def git(source, *args):
    return subprocess.check_output(["git", "-C", str(source), *args], text=True).strip()


def inspect_source(source, commit):
    state = git_state(source)
    if state["commit"] != commit or state["dirty"] is not False:
        raise ValueError("release checkout must have the requested HEAD and be clean: " + str(source))
    return state


def inspect_build(source, build):
    cache = build / "CMakeCache.txt"
    if build_settings(build).get("CMAKE_BUILD_TYPE") != "Release":
        raise ValueError("existing build is not CMake Release")
    home = next((l.split("=", 1)[1] for l in cache.read_text().splitlines()
                 if l.startswith("CMAKE_HOME_DIRECTORY:INTERNAL=")), None)
    if home is None or Path(home).resolve() != source.resolve():
        raise ValueError("existing build belongs to another source directory")
    binaries = {name: digest(build / name) for name in PROGRAMS}
    if any(v is None for v in binaries.values()):
        raise ValueError("existing build is incomplete; use a fresh build path")
    return {"cmake_cache_sha256": digest(cache), "binaries": binaries}


@contextmanager
def test_build_link(source, build):
    """Some maintained Python tests still resolve <source>/build directly."""
    link = source / "build"
    created = False
    if link.exists() or link.is_symlink():
        if not link.is_symlink() or link.resolve() != build.resolve():
            raise ValueError("required-suite build link already points elsewhere")
    else:
        link.symlink_to(build.resolve(), target_is_directory=True)
        created = True
    try:
        yield
    finally:
        if created:
            link.unlink()


def prepare(args):
    donor = args.clone_from.resolve()
    commit = git(donor, "rev-parse", "--verify", args.commit + "^{commit}")
    short = commit[:8]
    source = args.software_root.resolve() if args.software_root else ROOT / "deps" / ("LinearCParty-" + short)
    build = args.build_dir.resolve() if args.build_dir else ROOT / "build" / ("linear-cparty-" + short)
    lock = args.lock.resolve() if args.lock else ROOT / "manifests" / ("software-pilot-" + short + ".json")
    out = args.out.resolve()
    if args.estimate_only and not args.dry_run:
        raise ValueError("--estimate-only requires --dry-run")
    if args.jobs < 1 or args.estimated_build_files < 1:
        raise ValueError("jobs and estimated-build-files must be positive")
    # Source/build overrides support existing releases in a read-only rehearsal.
    if not source.exists() and source != ROOT / "deps" / ("LinearCParty-" + short):
        raise ValueError("new clones must use deps/LinearCParty-<short>")
    if not build.exists() and build != ROOT / "build" / ("linear-cparty-" + short):
        raise ValueError("new builds must use build/linear-cparty-<short>")
    previous = json.loads(out.read_text()) if out.exists() else {}
    identity = {"commit": commit, "source": str(source), "build": str(build), "lock": str(lock)}
    if previous and (previous.get("release") != identity or previous.get("dry_run") != args.dry_run):
        raise ValueError("existing preparation report has another release or dry-run mode")
    report = {"schema_version": 1, "created_utc": now(), "dry_run": args.dry_run,
              "release": identity, "stages": {}}

    def save(stage, result):
        report["stages"][stage] = result
        write_json(out, report)

    tracked = len(git(donor, "ls-tree", "-r", "--name-only", commit).splitlines())
    entries = (tracked * 2 + 100 if not source.exists() else 0) + (args.estimated_build_files if not build.exists() else 0) + 2
    report["file_budget"] = check_budget({"estimated_shared_entries": entries,
        "basis": "clone tracked files plus directory allowance (2*tracked+100), explicit build allowance, one report plus atomic temporary; actual build may vary"},
        {"file_budget": {"fallback_remaining_files": args.fallback_remaining_files}})
    suite_receipt = None
    created_build = False
    with local_workspace(args.scratch_dir or os.environ.get("TMPDIR") or "/tmp") as work:
        env = {**os.environ, "TMPDIR": str(work), "TMP": str(work), "TEMP": str(work),
               "PYTHONDONTWRITEBYTECODE": "1", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "LC_ALL": "C",
               "CPARTY_FIXED_ENERGY_PARAM_FILE": str(source / PARAMETER)}

        def command(argv, cwd=None):
            log = work / "command.log"
            with log.open("wb") as stream:
                result = subprocess.run(list(map(str, argv)), cwd=cwd, env=env, stdout=stream, stderr=stream)
            with log.open("rb") as stream:
                stream.seek(max(0, log.stat().st_size - 4000))
                tail = stream.read().decode("utf-8", "replace")
            record = {"argv": list(map(str, argv)), "returncode": result.returncode,
                      "log_sha256": digest(log), "log_tail": tail}
            if result.returncode:
                save("failed_command", record)
                raise ValueError("release preparation command failed: " + str(argv[0]) + " (see JSON log tail)")
            return record

        if source.exists():
            save("clone", {"status": "verified_existing", "state": inspect_source(source, commit)})
        elif args.dry_run:
            save("clone", {"status": "planned", "donor": str(donor)})
        else:
            command(["git", "clone", "--no-hardlinks", "--no-checkout", donor, source])
            command(["git", "-C", source, "checkout", "--detach", commit])
            save("clone", {"status": "created", "state": inspect_source(source, commit)})

        if build.exists():
            snapshot = inspect_build(source, build)
            prior_build = previous.get("stages", {}).get("build", {})
            created_build = (prior_build.get("created_by_prepare") is True and
                             all(prior_build.get(k) == v for k, v in snapshot.items()))
            save("build", {"status": "verified_existing", "created_by_prepare": created_build, **snapshot})
        elif args.dry_run:
            save("build", {"status": "planned", "build_type": "Release", "parallel_build_jobs": args.jobs})
            snapshot = None
        else:
            command(["cmake", "-S", source, "-B", build, "-DCMAKE_BUILD_TYPE=Release", "-DBUILD_TESTING=ON"])
            # Do not build unmaintained optional targets through the default 'all'.
            command(["cmake", "--build", build, "-j", args.jobs, "--target", *PROGRAMS, *CHECKS])
            snapshot = inspect_build(source, build)
            created_build = True
            save("build", {"status": "created", "created_by_prepare": True, **snapshot})

        if lock.exists():
            if not source.exists() or not build.exists():
                raise ValueError("an existing lock requires its source and build")
            locked = verify_lock(lock, source, build, [{"name": n} for n in PROGRAMS], source / PARAMETER)
            if (locked.get("build_revision_verified") is not True or
                    "required set ok" not in str(locked.get("required_suite", ""))):
                raise ValueError("existing lock must certify its revision and required suite")
            suite_receipt = {"status": "verified_existing_lock", "summary": locked["required_suite"], "lock_sha256": digest(lock)}
        else:
            prior_suite = previous.get("stages", {}).get("required_suite", {})
            prior_build = previous.get("stages", {}).get("build", {})
            if (prior_suite.get("status") == "passed" and snapshot is not None and
                    all(prior_build.get(k) == v for k, v in snapshot.items()) and
                    prior_suite.get("script_sha256") == digest(source / "test/tools/run_required_scfg2_tests.sh")):
                suite_receipt = prior_suite
            elif args.dry_run:
                suite_receipt = {"status": "planned", "env": "SCFG2_BUILD_DIR=<build>; temporary source/build symlink; node-local TMPDIR"}
            elif not created_build:
                raise ValueError("existing build has no verified required-suite receipt; refusing to modify it")
            else:
                with test_build_link(source, build):
                    env["SCFG2_BUILD_DIR"] = str(build)
                    receipt = command(["bash", source / "test/tools/run_required_scfg2_tests.sh"], cwd=source)
                if "[scfg2-required] required set ok" not in receipt["log_tail"]:
                    raise ValueError("required suite returned without its success marker")
                suite_receipt = {"status": "passed", **receipt,
                                 "script_sha256": digest(source / "test/tools/run_required_scfg2_tests.sh")}
                inspect_source(source, commit)
                # The suite may build additional maintained targets.
                snapshot = inspect_build(source, build)
                save("build", {"status": "created", "created_by_prepare": True, **snapshot})
        save("required_suite", suite_receipt)
        if lock.exists():
            save("lock", {"status": "verified_existing", "sha256": digest(lock)})
        elif args.dry_run:
            save("lock", {"status": "planned"})
        else:
            command([sys.executable, ROOT / "scripts/make_software_lock.py", "--software-root", source,
                     "--build-dir", build, "--out", lock, "--required-suite", "test/tools/run_required_scfg2_tests.sh: required set ok; receipt " + str(out)])
            verify_lock(lock, source, build, [{"name": n} for n in PROGRAMS], source / PARAMETER)
            save("lock", {"status": "created", "sha256": digest(lock)})

        parameter = source / PARAMETER
        if not source.exists():
            parameter = work / "requested-release-energy.par"
            parameter.write_bytes(subprocess.check_output(["git", "-C", str(donor), "show", commit + ":" + PARAMETER]))
        spec = specification(build, args.beam_reference_build.resolve(), args.exact_reference_build.resolve(),
                             [p.resolve() for p in args.input], parameter, list(map(int, args.beams.split(","))),
                             args.candidate_exact_build.resolve() if args.candidate_exact_build else None)
        spec["parameter"]["path"] = str(source / PARAMETER)
        prior = previous.get("stages", {}).get("comparison", {})
        if not args.dry_run and prior.get("specification") == spec and prior.get("successful") is True:
            save("comparison", prior)
        else:
            save("comparison", compare(spec, args.timeout, args.scratch_dir, args.dry_run))
        if source.exists():
            inspect_source(source, commit)
        if lock.exists():
            verify_lock(lock, source, build, [{"name": n} for n in PROGRAMS], source / PARAMETER)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("commit")
    parser.add_argument("--clone-from", type=Path, default=ROOT.parent.parent / "LinearCParty")
    parser.add_argument("--software-root", type=Path)
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--lock", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=8)
    parser.add_argument("--estimated-build-files", type=int, default=20000)
    parser.add_argument("--fallback-remaining-files", type=int, help="explicit conservative fallback if lfs quota is unavailable")
    add_arguments(parser)
    report = prepare(parser.parse_args())
    print(json.dumps({"release": report["release"], "dry_run": report["dry_run"], "file_budget": report["file_budget"],
                      "stages": {k: v.get("status", "planned" if v.get("dry_run") else "compared") for k, v in report["stages"].items()}}))
    return 0 if report["stages"]["comparison"].get("successful", True) else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, OSError, subprocess.CalledProcessError) as exc:
        print("error: " + str(exc), file=sys.stderr)
        sys.exit(2)
