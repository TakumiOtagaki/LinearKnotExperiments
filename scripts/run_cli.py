#!/usr/bin/env python3
"""Run a sequential CLI matrix and retain provenance and failed jobs."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import re
import resource
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime, timezone


ROOT = Path(__file__).resolve().parents[1]
TIME = Path("/usr/bin/time")
ENV_KEYS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "LC_ALL", "LANG", "JOB_ID", "SGE_TASK_ID", "NSLOTS")
BUILD_KEYS = ("CMAKE_BUILD_TYPE", "CMAKE_CXX_COMPILER", "CMAKE_CXX_FLAGS",
              "CMAKE_CXX_FLAGS_RELEASE", "CMAKE_GENERATOR")


def now():
    return datetime.now(timezone.utc).isoformat()


def digest(path):
    if not path.is_file():
        return None
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(path)


def git_state(path):
    def git(*args):
        result = subprocess.run(["git", "-C", str(path), *args],
                                capture_output=True, text=True)
        return result.stdout.strip() if result.returncode == 0 else None
    changes = git("status", "--porcelain", "--untracked-files=normal")
    diff = git("diff", "--binary", "HEAD")
    return {"commit": git("rev-parse", "--verify", "HEAD"),
            "dirty": bool(changes) if changes is not None else None,
            "status": changes,
            "tracked_diff_sha256": hashlib.sha256(diff.encode()).hexdigest()
            if diff is not None else None}


def build_settings(build):
    cache = build / "CMakeCache.txt"
    settings = {}
    if cache.is_file():
        for line in cache.read_text().splitlines():
            if ":" not in line or "=" not in line or line.startswith("//"):
                continue
            key = line.split(":", 1)[0]
            if key in BUILD_KEYS:
                settings[key] = line.split("=", 1)[1]
    return settings


def stop_group(process):
    try:
        os.killpg(process.pid, signal.SIGTERM)
        process.wait(timeout=1)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait()
    except ProcessLookupError:
        process.wait()


def memory_cap(mib):
    def apply():
        limit = mib * 1024 * 1024
        resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    return apply


def verify_lock(path, software, build, programs, parameter):
    lock = json.loads(path.read_text())
    state = git_state(software)
    if lock.get("schema_version") != 1 or state["dirty"] or state["commit"] != lock["source_commit"]:
        raise ValueError("source differs from clean software lock")
    if digest(build / "CMakeCache.txt") != lock["cmake_cache_sha256"]:
        raise ValueError("build cache differs from software lock")
    if digest(parameter) != lock["energy_parameter_sha256"]:
        raise ValueError("energy parameters differ from software lock")
    for program in programs:
        name = program["name"]
        if digest(build / name) != lock["binaries"].get(name) or lock["binaries"].get(name) is None:
            raise ValueError("binary differs from software lock: " + name)
    return lock


def run_one(command, directory, env, timeout, memory_limit_mib=None):
    write_json(directory / "command.json", {"argv": command, "cwd": str(directory)})
    started = time.perf_counter()
    status, error = "success", None
    timed_command = [str(TIME), "-f", "%e\t%M", "-o",
                     str(directory / "time.tsv"), "--", *command]
    with (directory / "stdout.txt").open("wb") as stdout, \
            (directory / "stderr.txt").open("wb") as stderr:
        try:
            process = subprocess.Popen(timed_command, cwd=directory, env=env,
                                       stdout=stdout, stderr=stderr,
                                       preexec_fn=memory_cap(memory_limit_mib) if memory_limit_mib else None,
                                       start_new_session=True)
        except OSError as exc:
            status, error, returncode = "launch_error", str(exc), None
        else:
            try:
                returncode = process.wait(timeout=timeout)
                if returncode != 0:
                    status = "failed"
            except subprocess.TimeoutExpired:
                stop_group(process)
                status, returncode = "timeout", process.returncode
            except KeyboardInterrupt:
                stop_group(process)
                status, returncode = "interrupted", process.returncode
    elapsed = time.perf_counter() - started
    wall, rss = None, None
    stamp = directory / "time.tsv"
    if stamp.is_file():
        for line in reversed(stamp.read_text().splitlines()):
            fields = line.split("\t")
            if len(fields) == 2:
                try:
                    wall, rss = float(fields[0]), int(fields[1])
                except ValueError:
                    continue
                break
    return {"status": status, "returncode": returncode, "error": error,
            "elapsed_seconds": elapsed, "time_wall_seconds": wall,
            "peak_rss_kib": rss}


def validate(config):
    if config.get("schema_version") != 1:
        raise ValueError("config schema_version must be 1")
    if config.get("purpose") != "development":
        raise ValueError("initial runner supports development runs; freeze protocol first")
    for group in ("inputs", "modes", "programs"):
        entries = config.get(group, [])
        if not entries:
            raise ValueError("config needs " + group)
        names = [entry["name"] for entry in entries]
        if len(set(names)) != len(names) or any(
                not re.fullmatch(r"[A-Za-z0-9_-]+", name) for name in names):
            raise ValueError(group + " needs unique simple names")
        for entry in entries:
            args = entry.get("args", [])
            if not isinstance(args, list) or any(not isinstance(x, str) for x in args):
                raise ValueError("args must be a list of strings")
    for program in config["programs"]:
        for product in program.get("outputs", []):
            name = product["name"]
            if Path(name).name != name or name in (".", ".."):
                raise ValueError("output must be a file name")
    for key, minimum in (("warmups", 0), ("repetitions", 1)):
        value = config.get(key, minimum)
        if type(value) is not int or value < minimum:
            raise ValueError(key + " is out of range")
    order = config.get("order", {"kind": "fixed"})
    if not isinstance(order, dict) or order.get("kind") not in ("fixed", "blocked_random"):
        raise ValueError("order.kind must be fixed or blocked_random")
    if order["kind"] == "blocked_random" and type(order.get("seed")) is not int:
        raise ValueError("blocked_random order needs an integer seed")
    if config.get("timeout_seconds", 30) <= 0:
        raise ValueError("timeout_seconds must be positive")
    if "memory_limit_mib" in config and (type(config["memory_limit_mib"]) is not int or config["memory_limit_mib"] < 64):
        raise ValueError("memory_limit_mib must be an integer of at least 64")
    environment = config.get("environment", {})
    if not isinstance(environment, dict) or any(
            key not in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "LC_ALL")
            or not isinstance(value, str) for key, value in environment.items()):
        raise ValueError("environment supports string thread-count and LC_ALL overrides")


def experiment_path(value):
    path = Path(value)
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def plan_jobs(config):
    """Keep the old nested order by default; shuffle each repetition otherwise."""
    combinations = [(item["name"], mode["name"], program["name"])
                    for item in config["inputs"] for mode in config["modes"]
                    for program in config["programs"]]
    jobs = []

    def add(combination, warmup, repetition):
        item, mode, program = combination
        jobs.append({"id": "{}-{}-{}-{}{}".format(
            item, mode, program, "warmup" if warmup else "rep", repetition),
            "input": item, "mode": mode, "program": program,
            "warmup": warmup, "repetition": repetition})

    warmups, repetitions = config.get("warmups", 0), config.get("repetitions", 1)
    if config.get("order", {}).get("kind", "fixed") == "fixed":
        for combination in combinations:
            for repetition in range(warmups):
                add(combination, True, repetition)
            for repetition in range(repetitions):
                add(combination, False, repetition)
    else:
        rng = random.Random(config["order"]["seed"])
        for repetition in range(repetitions):
            block = combinations.copy()
            rng.shuffle(block)
            for combination in block:
                if repetition == 0:
                    for index in range(warmups):
                        add(combination, True, index)
                add(combination, False, repetition)
    return jobs


def resume_manifest(output, config_path, config, software, build, parameter, sources,
                    lock_path, lock, state, binaries, plan):
    """Validate all dependencies before changing any existing run files."""
    manifest = json.loads((output / "run.json").read_text())

    def require(condition, label):
        if not condition:
            raise ValueError("cannot resume: " + label + " differs from run.json")

    require(manifest.get("schema_version") == 1, "schema")
    require(manifest["config"] == config
            and manifest["original_config_sha256"] == digest(config_path)
            and manifest["config_sha256"] == digest(output / "config.json")
            and json.loads((output / "config.json").read_text()) == config, "config")
    require(manifest["software_root"] == str(software)
            and manifest["build_dir"] == str(build), "software/build paths")
    require(manifest["software_lock"] == lock
            and manifest["software_lock_sha256"] == (digest(lock_path) if lock_path else None),
            "software lock")
    require(manifest["software_before"] == state, "software state")
    require(manifest["binaries_before"] == binaries, "binary hashes")
    require(manifest["build_settings"] == build_settings(build), "build settings")
    if "cmake_cache_sha256" in manifest:
        require(manifest["cmake_cache_sha256"] == digest(build / "CMakeCache.txt"), "build cache")
    require(manifest["energy_parameter_sha256"] == digest(parameter)
            == digest(output / "inputs/energy.par"), "energy parameter hash")
    require(len(manifest["inputs"]) == len(sources), "input count")
    for saved, (entry, source) in zip(manifest["inputs"], sources):
        copied = output / "inputs" / (entry["name"] + ".txt")
        require(all(saved.get(key) == value for key, value in entry.items())
                and saved["copied_path"] == str(copied)
                and saved["sha256"] == digest(source) == digest(copied), "input " + entry["name"])
    require(manifest.get("planned_jobs", plan) == plan, "job order")
    completed = []
    for job in plan:
        directory = Path("jobs") / job["id"]
        result_path = output / directory / "result.json"
        if result_path.is_file():
            result = json.loads(result_path.read_text())
            require(all(result.get(key) == value for key, value in job.items())
                    and result.get("directory") == str(directory), "job result " + job["id"])
            completed.append(result)
    manifest["jobs"] = completed
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/smoke.json")
    parser.add_argument("--software-root", type=Path, default=ROOT / "deps/LinearCParty")
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--out", required=True)
    parser.add_argument("--resume", action="store_true",
                        help="continue an existing matching run; keep all jobs with result.json")
    args = parser.parse_args()
    config_path = experiment_path(args.config)
    config = json.loads(config_path.read_text())
    validate(config)
    software = args.software_root.resolve()
    build = (args.build_dir or software / "build").resolve()
    if not software.is_dir() or not build.is_dir():
        parser.error("supply existing --software-root and --build-dir paths")
    if not TIME.is_file():
        parser.error("GNU /usr/bin/time is required")
    parameter = software / config["energy_parameter_file"]
    if not parameter.is_file():
        parser.error("energy parameter file is missing: " + str(parameter))
    lock_path = experiment_path(config["software_lock"]) if config.get("software_lock") else None
    lock = verify_lock(lock_path, software, build, config["programs"], parameter) if lock_path else None
    sources = [(entry, experiment_path(entry["path"])) for entry in config["inputs"]]
    if any(not path.is_file() for _, path in sources):
        parser.error("an input file is missing")
    output = experiment_path(args.out)
    plan = plan_jobs(config)
    state = git_state(software)
    binaries = {p["name"]: digest(build / p["name"]) for p in config["programs"]}
    copied_parameter = output / "inputs/energy.par"
    env = os.environ.copy()
    env.update(config.get("environment", {}))
    env["CPARTY_FIXED_ENERGY_PARAM_FILE"] = str(copied_parameter)
    host = {"hostname": platform.node(), "platform": platform.platform(),
            "python": platform.python_version(), "cpu_count": os.cpu_count()}
    environment = {key: env[key] for key in ENV_KEYS if key in env}
    if args.resume:
        manifest = resume_manifest(output, config_path, config, software, build, parameter,
                                   sources, lock_path, lock, state, binaries, plan)
        resumed = now()
        manifest.setdefault("resume_history", []).append({
            "resumed_utc": resumed, "previous_status": manifest.get("status"),
            "previous_finished_utc": manifest.get("finished_utc"),
            "previous_checks": {key: manifest[key] for key in (
                "software_after", "binaries_after", "experiment_after",
                "software_changed_during_run", "binaries_changed_during_run",
                "software_lock_valid_after_run", "software_lock_error_after_run")
                if key in manifest},
            "completed_jobs": len(manifest["jobs"]),
            "pending_jobs": len(plan) - len(manifest["jobs"]),
            "host": host, "environment": environment, "experiment": git_state(ROOT),
            "runner_sha256": digest(Path(__file__))})
        manifest["resumed_utc"] = resumed
        # Earlier runners recorded results only, so recover the known legacy order.
        if "execution_order" not in manifest:
            manifest["execution_order"] = [
                {"id": job["id"], "started_utc": None, "recovered_legacy": True}
                for job in plan if (output / "jobs" / job["id"] / "command.json").is_file()]
        for key in ("finished_utc", "software_after", "experiment_after", "binaries_after",
                    "software_changed_during_run", "binaries_changed_during_run",
                    "software_lock_valid_after_run", "software_lock_error_after_run"):
            manifest.pop(key, None)
        manifest["status"] = "running"
    else:
        output.mkdir(parents=True, exist_ok=False)
        copied_parameter.parent.mkdir()
        shutil.copyfile(parameter, copied_parameter)
        frozen_inputs = []
        for entry, source in sources:
            copied = copied_parameter.parent / (entry["name"] + ".txt")
            shutil.copyfile(source, copied)
            frozen_inputs.append({**entry, "copied_path": str(copied), "sha256": digest(copied)})
        write_json(output / "config.json", config)
        manifest = {"schema_version": 1, "status": "running", "started_utc": now(),
                "purpose": config["purpose"], "config": config,
                "config_sha256": digest(output / "config.json"),
                "original_config_sha256": digest(config_path),
                "software_root": str(software), "build_dir": str(build),
                "software_before": state,
                "experiment_before": git_state(ROOT), "binaries_before": binaries,
                "build_settings": build_settings(build), "inputs": frozen_inputs,
                "cmake_cache_sha256": digest(build / "CMakeCache.txt"),
                "software_lock": lock, "software_lock_sha256": digest(lock_path) if lock_path else None,
                "energy_parameter_sha256": digest(copied_parameter),
                "host": host, "environment": environment,
                "execution_order": [], "jobs": []}
    manifest["planned_jobs"] = plan
    write_json(output / "run.json", manifest)
    completed = {job["id"] for job in manifest["jobs"]}
    failed = any(job["status"] != "success" for job in manifest["jobs"])
    interrupted = False
    items = {item["name"]: item for item in manifest["inputs"]}
    modes = {mode["name"]: mode for mode in config["modes"]}
    programs = {program["name"]: program for program in config["programs"]}
    for job in plan:
        if job["id"] in completed:
            continue
        directory = output / "jobs" / job["id"]
        if directory.exists():
            archive = output / "incomplete_attempts" / job["id"]
            archive.mkdir(parents=True, exist_ok=True)
            index = 0
            while (archive / str(index)).exists():
                index += 1
            attempt = archive / str(index)
            shutil.move(str(directory), str(attempt))
        directory.mkdir(parents=True)
        item, mode, program = items[job["input"]], modes[job["mode"]], programs[job["program"]]
        command = [str(build / program["name"]), "-i", item["copied_path"],
                   *mode.get("args", []), *program.get("args", [])]
        for product in program.get("outputs", []):
            command.extend([product["flag"], str(directory / product["name"])])
        attempt = {"id": job["id"], "started_utc": now()}
        manifest["execution_order"].append(attempt)
        write_json(output / "run.json", manifest)
        result = run_one(command, directory, env, config.get("timeout_seconds", 30),
                         config.get("memory_limit_mib"))
        result.update({**job, "command": command, "directory": str(directory.relative_to(output))})
        write_json(directory / "result.json", result)
        manifest["jobs"].append(result)
        attempt.update({"finished_utc": now(), "status": result["status"]})
        write_json(output / "run.json", manifest)
        print(job["id"] + ": " + result["status"], flush=True)
        failed |= result["status"] != "success"
        interrupted = result["status"] == "interrupted"
        if interrupted:
            break
    manifest.update({"finished_utc": now(), "status": "interrupted" if interrupted else
                     "failed" if failed else "success", "software_after": git_state(software),
                     "experiment_after": git_state(ROOT),
                     "binaries_after": {p["name"]: digest(build / p["name"])
                                        for p in config["programs"]}})
    manifest["software_changed_during_run"] = manifest["software_before"] != manifest["software_after"]
    manifest["binaries_changed_during_run"] = manifest["binaries_before"] != manifest["binaries_after"]
    if lock_path:
        try:
            verify_lock(lock_path, software, build, config["programs"], parameter)
            manifest["software_lock_valid_after_run"] = digest(lock_path) == manifest["software_lock_sha256"]
        except ValueError as exc:
            manifest["software_lock_valid_after_run"] = False
            manifest["software_lock_error_after_run"] = str(exc)
        if not manifest["software_lock_valid_after_run"]:
            failed = True
            manifest["status"] = "failed"
    write_json(output / "run.json", manifest)
    return 130 if interrupted else 1 if failed else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, OSError) as exc:
        print("error: " + str(exc), file=sys.stderr)
        sys.exit(2)
