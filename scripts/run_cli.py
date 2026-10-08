#!/usr/bin/env python3
"""Run a sequential CLI matrix and retain provenance and failed jobs."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
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


def run_one(command, directory, env, timeout):
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
    if config.get("timeout_seconds", 30) <= 0:
        raise ValueError("timeout_seconds must be positive")


def experiment_path(value):
    path = Path(value)
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/smoke.json")
    parser.add_argument("--software-root", type=Path, default=ROOT / "deps/LinearCParty")
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--out", required=True)
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
    sources = [(entry, experiment_path(entry["path"])) for entry in config["inputs"]]
    if any(not path.is_file() for _, path in sources):
        parser.error("an input file is missing")
    output = experiment_path(args.out)
    output.mkdir(parents=True, exist_ok=False)
    inputs_dir = output / "inputs"
    inputs_dir.mkdir()
    copied_parameter = inputs_dir / "energy.par"
    shutil.copyfile(parameter, copied_parameter)
    frozen_inputs = []
    for entry, source in sources:
        copied = inputs_dir / (entry["name"] + ".txt")
        shutil.copyfile(source, copied)
        frozen_inputs.append({**entry, "copied_path": str(copied), "sha256": digest(copied)})
    write_json(output / "config.json", config)
    binaries = {p["name"]: digest(build / p["name"]) for p in config["programs"]}
    manifest = {"schema_version": 1, "status": "running", "started_utc": now(),
                "purpose": config["purpose"], "config": config,
                "config_sha256": digest(output / "config.json"),
                "original_config_sha256": digest(config_path),
                "software_root": str(software), "build_dir": str(build),
                "software_before": git_state(software),
                "experiment_before": git_state(ROOT), "binaries_before": binaries,
                "build_settings": build_settings(build), "inputs": frozen_inputs,
                "energy_parameter_sha256": digest(copied_parameter),
                "host": {"hostname": platform.node(), "platform": platform.platform(),
                         "python": platform.python_version(), "cpu_count": os.cpu_count()},
                "environment": {key: os.environ[key] for key in ENV_KEYS if key in os.environ},
                "jobs": []}
    env = os.environ.copy()
    env["CPARTY_FIXED_ENERGY_PARAM_FILE"] = str(copied_parameter)
    write_json(output / "run.json", manifest)
    failed = interrupted = False
    for item in frozen_inputs:
        for mode in config["modes"]:
            for program in config["programs"]:
                for index in range(config.get("warmups", 0) + config.get("repetitions", 1)):
                    warmup = index < config.get("warmups", 0)
                    repetition = index if warmup else index - config.get("warmups", 0)
                    job_id = "{}-{}-{}-{}{}".format(item["name"], mode["name"],
                               program["name"], "warmup" if warmup else "rep", repetition)
                    directory = output / "jobs" / job_id
                    directory.mkdir(parents=True)
                    command = [str(build / program["name"]), "-i", item["copied_path"],
                               *mode.get("args", []), *program.get("args", [])]
                    for product in program.get("outputs", []):
                        command.extend([product["flag"], str(directory / product["name"])])
                    result = run_one(command, directory, env, config.get("timeout_seconds", 30))
                    result.update({"id": job_id, "input": item["name"], "mode": mode["name"],
                                   "program": program["name"], "warmup": warmup,
                                   "repetition": repetition, "command": command,
                                   "directory": str(directory.relative_to(output))})
                    write_json(directory / "result.json", result)
                    manifest["jobs"].append(result)
                    write_json(output / "run.json", manifest)
                    print(job_id + ": " + result["status"], flush=True)
                    failed |= result["status"] != "success"
                    interrupted |= result["status"] == "interrupted"
                    if interrupted:
                        break
                if interrupted:
                    break
            if interrupted:
                break
        if interrupted:
            break
    manifest.update({"finished_utc": now(), "status": "interrupted" if interrupted else
                     "failed" if failed else "success", "software_after": git_state(software),
                     "experiment_after": git_state(ROOT),
                     "binaries_after": {p["name"]: digest(build / p["name"])
                                        for p in config["programs"]}})
    manifest["software_changed_during_run"] = manifest["software_before"] != manifest["software_after"]
    manifest["binaries_changed_during_run"] = manifest["binaries_before"] != manifest["binaries_after"]
    write_json(output / "run.json", manifest)
    return 130 if interrupted else 1 if failed else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, OSError) as exc:
        print("error: " + str(exc), file=sys.stderr)
        sys.exit(2)
