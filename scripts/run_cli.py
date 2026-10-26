#!/usr/bin/env python3
"""Run a CLI matrix and retain provenance and failed jobs."""

import argparse
import hashlib
import json
import multiprocessing
from multiprocessing.connection import wait
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
from task_storage import FORMAT, ArchiveReader, TaskStore, local_workspace, jsonl_rows


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
        name = program.get("executable", program["name"])
        if digest(build / name) != lock["binaries"].get(name) or lock["binaries"].get(name) is None:
            raise ValueError("binary differs from software lock: " + name)
    # Relocatable adapters exec other binaries from the same staged build.
    # Pin those dependencies too, even when only the adapter is requested.
    for name, expected in lock.get("build_files", {}).items():
        if Path(name).name != name or expected is None or digest(build / name) != expected:
            raise ValueError("build dependency differs from software lock: " + name)
    for entry in lock.get("extra_files", []):
        if digest(experiment_path(entry["path"])) != entry["sha256"]:
            raise ValueError("extra dependency differs from software lock: " + entry["path"])
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


def validate(config, array=False):
    if config.get("schema_version") != 1:
        raise ValueError("config schema_version must be 1")
    if config.get("purpose") not in ("development", "formal"):
        raise ValueError("purpose must be development or formal")
    if config["purpose"] == "formal":
        if not array:
            raise ValueError("formal measurements require the array runner")
        if not config.get("software_lock") or config.get("model_settings", {}).get("build_revision_verified") is not True:
            raise ValueError("formal measurements require a software lock and verified build revision")
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
            if group == "modes" and "length_range" in entry:
                bounds = entry["length_range"]
                if (not isinstance(bounds, dict) or
                        any(type(bounds.get(k)) is not int for k in ("min_length", "max_length")) or
                        not 1 <= bounds["min_length"] <= bounds["max_length"]):
                    raise ValueError("mode length_range needs positive inclusive bounds")
    for program in config["programs"]:
        for key in ("executable", "subcommand"):
            if key in program and (not isinstance(program[key], str) or
                                   not re.fullmatch(r"[A-Za-z0-9_-]+", program[key])):
                raise ValueError(key + " must be a simple executable or subcommand name")
        for product in program.get("outputs", []):
            name = product["name"]
            if Path(name).name != name or name in (".", ".."):
                raise ValueError("output must be a file name")
    for key, minimum in (("warmups", 0), ("repetitions", 1), ("manifest_checkpoint_jobs", 1)):
        value = config.get(key, minimum)
        if type(value) is not int or value < minimum:
            raise ValueError(key + " is out of range")
    order = config.get("order", {"kind": "fixed"})
    if not isinstance(order, dict) or order.get("kind") not in ("fixed", "blocked_random", "beam_ascending"):
        raise ValueError("order.kind must be fixed, blocked_random or beam_ascending")
    if order["kind"] in ("blocked_random", "beam_ascending") and type(order.get("seed")) is not int:
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
                    if mode_applies(mode, item) for program in config["programs"]]
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
            if config["order"]["kind"] == "beam_ascending":
                widths = {m["name"]: beam_width(m) for m in config["modes"]}
                block.sort(key=lambda c: widths[c[1]])
            for combination in block:
                if repetition == 0:
                    for index in range(warmups):
                        add(combination, True, index)
                add(combination, False, repetition)
    return jobs


def beam_width(mode):
    args = mode.get("args", [])
    return int(args[args.index("-b") + 1]) if "-b" in args else 0


def mode_applies(mode, item):
    bounds = mode.get("length_range")
    if bounds is None:
        return True
    length = item.get("length")
    if type(length) is not int:
        from rna_records import read_record
        _, sequence, _ = read_record(experiment_path(item["path"]))
        length = len(sequence)
    return bounds["min_length"] <= length <= bounds["max_length"]


def resume_manifest(output, config_path, config, software, build, parameter, sources,
                    lock_path, lock, state, binaries, plan, store=None):
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
    require(manifest["energy_parameter_sha256"] == digest(parameter), "energy parameter hash")
    if store is None:
        require(manifest["energy_parameter_sha256"] == digest(output / "inputs/energy.par"), "copied energy parameter hash")
    require(len(manifest["inputs"]) == len(sources), "input count")
    for saved, (entry, source) in zip(manifest["inputs"], sources):
        copied = output / "inputs" / (entry["name"] + ".txt")
        require(all(saved.get(key) == value for key, value in entry.items())
                and (store is not None or saved["copied_path"] == str(copied))
                and saved["sha256"] == digest(source)
                and (store is not None or saved["sha256"] == digest(copied)), "input " + entry["name"])
    require(manifest.get("planned_jobs", plan) == plan, "job order")
    completed = []
    if store is not None:
        require(manifest.get("output_format") == FORMAT, "output format (old arrays must be replanned)")
        planned = {j["id"]: j for j in plan}
        for result in store.rows:
            require(result["id"] in planned and all(result.get(k) == v for k, v in planned[result["id"]].items()),
                    "checkpoint job " + result["id"])
        completed = store.rows.copy()
    for job in plan:
        if store is not None:
            break
        directory = Path("jobs") / job["id"]
        result_path = output / directory / "result.json"
        if result_path.is_file():
            result = json.loads(result_path.read_text())
            require(all(result.get(key) == value for key, value in job.items())
                    and result.get("directory") == str(directory), "job result " + job["id"])
            completed.append(result)
    manifest["jobs"] = completed
    journal = output / "execution_events.jsonl"
    if journal.exists():
        attempts = []
        for event in jsonl_rows(journal):
            kind = event.pop("event")
            if kind == "started":
                attempts.append(event)
            else:
                match = next((a for a in reversed(attempts) if a["id"] == event["id"] and "finished_utc" not in a), None)
                if match is None:
                    raise ValueError("cannot resume: execution journal has no matching start")
                match.update(event)
        manifest["execution_order"] = attempts
    return manifest


def job_process(connection, command, directory, env, timeout, memory):
    """A separate controller keeps RLIMIT_AS setup out of threaded Python."""
    try:
        connection.send(run_one(command, directory, env, timeout, memory))
    except BaseException as exc:
        connection.send({"status": "launch_error", "returncode": None, "error": str(exc),
                         "elapsed_seconds": None, "time_wall_seconds": None, "peak_rss_kib": None})
    finally:
        connection.close()


def execute_jobs(manifest, output, build, env, workers, through_beam=None, workspace=None, store=None):
    """Launch in saved order; serialize each combination and repetition block."""
    config = manifest["config"]
    completed = {job["id"] for job in manifest["jobs"]}
    pending = [job for job in manifest["planned_jobs"] if job["id"] not in completed]
    items = {item["name"]: item for item in manifest["inputs"]}
    modes = {mode["name"]: mode for mode in config["modes"]}
    if through_beam is not None:
        pending = [j for j in pending if beam_width(modes[j["mode"]]) <= through_beam]
    programs = {program["name"]: program for program in config["programs"]}
    active = {}
    interrupted = False
    checkpoint = config.get("manifest_checkpoint_jobs", 1)

    def save_progress(event):
        if checkpoint > 1 or store is not None:
            with (output / "execution_events.jsonl").open("a") as stream:
                stream.write(json.dumps(event) + "\n")
                stream.flush()
                if store is not None:
                    os.fsync(stream.fileno())
        if checkpoint == 1 or len(manifest["jobs"]) % checkpoint == 0:
            write_json(output / "run.json", manifest)

    def finish(job, directory, attempt, result):
        nonlocal interrupted
        result.update({**job, "command": attempt.pop("command"),
                       "started_utc": attempt["started_utc"], "finished_utc": now(),
                       "directory": str(directory.relative_to(workspace or output))})
        write_json(directory / "result.json", result)
        if store is not None:
            result = store.commit(result, {str(p.relative_to(workspace)): p
                                          for p in directory.rglob("*") if p.is_file()})
        manifest["jobs"].append(result)
        attempt.update({"finished_utc": now(), "status": result["status"]})
        save_progress({"event": "finished", **attempt})
        print(job["id"] + ": " + result["status"], flush=True)
        interrupted |= result["status"] == "interrupted"
        if store is not None:
            shutil.rmtree(directory)

    def collect():
        for connection in wait(list(active), timeout=0.2):
            process, job, directory, attempt = active.pop(connection)
            try:
                result = connection.recv()
            except EOFError:
                result = {"status": "launch_error", "returncode": process.exitcode,
                          "error": "job controller exited without a result",
                          "elapsed_seconds": None, "time_wall_seconds": None, "peak_rss_kib": None}
            process.join()
            connection.close()
            finish(job, directory, attempt, result)

    try:
        while pending or active:
            if pending and not interrupted and len(active) < workers:
                job = pending[0]
                combination = tuple(job[k] for k in ("input", "mode", "program"))
                block = 0 if job["warmup"] else job["repetition"]
                blocked = any(tuple(a[1][k] for k in ("input", "mode", "program")) == combination
                              or (config.get("order", {}).get("kind") in ("blocked_random", "beam_ascending")
                                  and (0 if a[1]["warmup"] else a[1]["repetition"]) != block)
                              or (config.get("order", {}).get("kind") == "beam_ascending" and
                                  beam_width(modes[a[1]["mode"]]) != beam_width(modes[job["mode"]]))
                              for a in active.values())
                if not blocked:
                    pending.pop(0)
                    directory = (workspace or output) / "jobs" / job["id"]
                    if directory.exists():
                        archive = (workspace or output) / "incomplete_attempts" / job["id"]
                        archive.mkdir(parents=True, exist_ok=True)
                        index = 0
                        while (archive / str(index)).exists():
                            index += 1
                        shutil.move(str(directory), str(archive / str(index)))
                    directory.mkdir(parents=True)
                    item, mode, program = items[job["input"]], modes[job["mode"]], programs[job["program"]]
                    command = [str(build / program.get("executable", program["name"])),
                               *([program["subcommand"]] if "subcommand" in program else []),
                               "-i", item["copied_path"],
                               *mode.get("args", []), *program.get("args", [])]
                    for product in program.get("outputs", []):
                        command.extend([product["flag"], str(directory / product["name"])])
                    attempt = {"id": job["id"], "started_utc": now()}
                    manifest["execution_order"].append(attempt)
                    save_progress({"event": "started", **attempt})
                    attempt["command"] = command
                    if workers == 1:
                        finish(job, directory, attempt,
                               run_one(command, directory, env, config.get("timeout_seconds", 30),
                                       config.get("memory_limit_mib")))
                    else:
                        context = multiprocessing.get_context("fork")
                        receive, send = context.Pipe(duplex=False)
                        process = context.Process(target=job_process, args=(
                            send, command, directory, env, config.get("timeout_seconds", 30),
                            config.get("memory_limit_mib")))
                        process.start()
                        send.close()
                        active[receive] = (process, job, directory, attempt)
                    continue
            if active:
                collect()
            elif interrupted:
                break
    except KeyboardInterrupt:
        interrupted = True
        # Controllers catch SIGINT inside run_one and stop each CLI process group.
        for process, *_ in active.values():
            if process.is_alive():
                os.kill(process.pid, signal.SIGINT)
        while active:
            collect()
    except BaseException:
        # If durable shared storage fails, stop outstanding local controllers.
        # Their uncommitted jobs remain pending for the next allocation.
        for process, *_ in active.values():
            if process.is_alive():
                os.kill(process.pid, signal.SIGINT)
        for connection, (process, *_) in active.items():
            process.join(timeout=5)
            if process.is_alive():
                process.terminate()
                process.join()
            connection.close()
        raise
    return interrupted


def execute(config_path, software, build, output, resume=False, array_context=None, workers=1,
            through_beam=None, scratch_dir=None, _workspace=None):
    if array_context is not None and _workspace is None:
        with local_workspace(scratch_dir) as workspace:
            return execute(config_path, software, build, output, resume, array_context,
                           workers, through_beam, scratch_dir, workspace)
    packed = _workspace is not None
    config = json.loads(config_path.read_text())
    if any(i.get("scaffold_pending") for i in config.get("inputs", [])):
        raise ValueError("pending scaffolds cannot be executed")
    validate(config, array=array_context is not None)
    if type(workers) is not int or workers < 1 or (workers != 1 and array_context is None):
        raise ValueError("parallel workers are available only through an array plan")
    if not software.is_dir() or not build.is_dir():
        raise ValueError("supply existing --software-root and --build-dir paths")
    if not TIME.is_file():
        raise ValueError("GNU /usr/bin/time is required")
    parameter = software / config["energy_parameter_file"]
    if not parameter.is_file():
        raise ValueError("energy parameter file is missing: " + str(parameter))
    lock_path = experiment_path(config["software_lock"]) if config.get("software_lock") else None
    lock = verify_lock(lock_path, software, build, config["programs"], parameter) if lock_path else None
    if config["purpose"] == "formal" and lock.get("build_revision_verified") is not True:
        raise ValueError("formal software lock must verify the build revision")
    sources = []
    archives = {}
    for entry in config["inputs"]:
        if "input_archive" in entry:
            archive_path = experiment_path(entry["input_archive"])
            if archive_path not in archives:
                archives[archive_path] = ArchiveReader(archive_path)
            source = _workspace / (entry["name"] + "-source.txt")
            data = archives[archive_path].read(entry["input_member"])
            if hashlib.sha256(data).hexdigest() != entry["sha256"]:
                raise ValueError("frozen archive input changed")
            source.write_bytes(data)
        else:
            source = experiment_path(entry["path"])
        sources.append((entry, source))
    if any(not path.is_file() for _, path in sources):
        raise ValueError("an input file is missing")
    plan = plan_jobs(config)
    state = git_state(software)
    binaries = {p["name"]: digest(build / p.get("executable", p["name"])) for p in config["programs"]}
    copied_parameter = (_workspace or output) / "inputs/energy.par"
    env = os.environ.copy()
    env.update(config.get("environment", {}))
    env["CPARTY_FIXED_ENERGY_PARAM_FILE"] = str(copied_parameter)
    host = {"hostname": platform.node(), "platform": platform.platform(),
            "python": platform.python_version(), "cpu_count": os.cpu_count()}
    environment = {key: env[key] for key in ENV_KEYS if key in env}
    if packed:
        output.mkdir(parents=True, exist_ok=resume)
    store = TaskStore(output) if packed else None
    if resume:
        manifest = resume_manifest(output, config_path, config, software, build, parameter,
                                   sources, lock_path, lock, state, binaries, plan, store)
        if array_context is not None and manifest.get("array_task") != array_context:
            raise ValueError("cannot resume: array task differs from run.json")
        if manifest.get("workers", 1) != workers:
            raise ValueError("cannot resume: workers differ from run.json")
        if store is not None:
            store = TaskStore(output, resume=True)
            jsonl_rows(output / "execution_events.jsonl", repair=True)
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
        if not packed:
            output.mkdir(parents=True, exist_ok=False)
        copied_parameter.parent.mkdir()
        shutil.copyfile(parameter, copied_parameter)
        frozen_inputs = []
        for entry, source in sources:
            copied = copied_parameter.parent / (entry["name"] + ".txt")
            shutil.copyfile(source, copied)
            frozen_inputs.append({**entry, "copied_path": str(copied), "sha256": digest(copied),
                                  **({"source_input_path": str(source)} if packed and "input_archive" not in entry else {})})
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
        if array_context is not None:
            manifest.update(array_task=array_context, workers=workers, output_format=FORMAT)
    if packed and resume:
        copied_parameter.parent.mkdir()
        shutil.copyfile(parameter, copied_parameter)
        for saved, (entry, source) in zip(manifest["inputs"], sources):
            copied = copied_parameter.parent / (entry["name"] + ".txt")
            shutil.copyfile(source, copied)
            saved["copied_path"] = str(copied)
    manifest["planned_jobs"] = plan
    write_json(output / "run.json", manifest)
    interrupted = execute_jobs(manifest, output, build, env, workers, through_beam, _workspace, store)
    interrupted |= len(manifest["jobs"]) < len(plan)
    failed = any(job["status"] != "success" for job in manifest["jobs"])
    manifest.update({"finished_utc": now(), "status": "interrupted" if interrupted else
                     "failed" if failed else "success", "software_after": git_state(software),
                     "experiment_after": git_state(ROOT),
                     "binaries_after": {p["name"]: digest(build / p.get("executable", p["name"]))
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/smoke.json")
    parser.add_argument("--software-root", type=Path, default=ROOT / "deps/LinearCParty")
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--out", required=True)
    parser.add_argument("--resume", action="store_true",
                        help="continue an existing matching run; keep all jobs with result.json")
    args = parser.parse_args()
    software = args.software_root.resolve()
    build = (args.build_dir or software / "build").resolve()
    return execute(experiment_path(args.config), software, build,
                   experiment_path(args.out), resume=args.resume)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, OSError) as exc:
        print("error: " + str(exc), file=sys.stderr)
        sys.exit(2)
