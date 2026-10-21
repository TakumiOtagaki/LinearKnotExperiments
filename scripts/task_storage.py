"""Durable task-sized archives and readers that never unpack shared results."""

from contextlib import contextmanager
import gzip
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import tarfile
import tempfile


FORMAT = "task-jsonl-tar-v1"


@contextmanager
def local_workspace(explicit=None):
    root = explicit or os.environ.get("TMPDIR")
    if not root:
        raise ValueError("node-local TMPDIR or explicit --scratch-dir is required")
    root = Path(root).resolve()
    if not root.is_dir():
        raise ValueError("scratch directory must already exist: " + str(root))
    # Prevent the common accidental fallback onto this project's shared storage.
    if str(root) == "/gs/bs" or str(root).startswith("/gs/bs/"):
        raise ValueError("scratch must be node-local, outside /gs/bs")
    with tempfile.TemporaryDirectory(prefix="lcp-task-", dir=root) as value:
        yield Path(value)


def jsonl_rows(path, repair=False):
    """Ignore only an unterminated last record; malformed complete lines fail."""
    rows, end = [], 0
    if not path.exists():
        return rows
    with path.open("rb") as stream:
        for line in stream:
            if not line.endswith(b"\n"):
                break
            try:
                rows.append(json.loads(line))
            except (ValueError, UnicodeError) as exc:
                raise ValueError("malformed checkpoint: " + str(path)) from exc
            end = stream.tell()
    if repair and path.stat().st_size != end:
        with path.open("r+b") as stream:
            stream.truncate(end)
    return rows


def file_summary(path):
    h, size, tail = hashlib.sha256(), 0, b""
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
            size += len(chunk)
            tail = (tail + chunk)[-2048:]
    return {"bytes": size, "sha256": h.hexdigest(), "tail": tail.decode("utf-8", "replace")}


class TaskStore:
    """Tar is fsynced before JSONL completion. Resume discards uncommitted tail."""

    def __init__(self, directory, resume=False):
        self.directory = Path(directory)
        self.checkpoint = self.directory / "task_results.jsonl"
        self.archive = self.directory / "task_outputs.tar"
        self.rows = jsonl_rows(self.checkpoint, repair=resume)
        if len({r["id"] for r in self.rows}) != len(self.rows):
            raise ValueError("duplicate checkpoint job")
        self.ids = {r["id"] for r in self.rows}
        self.end = 0
        for row in self.rows:
            end = row["archive_end"]
            if type(end) is not int or end <= self.end or end % 512:
                raise ValueError("invalid archive checkpoint offset")
            self.end = end
        if self.rows and (not self.archive.is_file() or self.archive.stat().st_size < self.end):
            raise ValueError("checkpoint raw archive is missing or truncated")
        if resume and self.archive.exists():
            with self.archive.open("r+b") as stream:
                stream.seek(self.end)
                stream.write(b"\0" * 1024)
                stream.truncate()
                stream.flush()
                os.fsync(stream.fileno())

    def commit(self, row, files):
        if row["id"] in self.ids:
            raise ValueError("job already checkpointed: " + row["id"])
        artifacts, channels = {}, {}
        with self.archive.open("r+b" if self.archive.exists() else "w+b") as stream:
            stream.seek(self.end)
            with tarfile.open(fileobj=stream, mode="w", format=tarfile.PAX_FORMAT) as archive:
                for name, path in files.items():
                    path = Path(path)
                    summary = file_summary(path)
                    artifacts[name] = {k: summary[k] for k in ("bytes", "sha256")}
                    for channel in ("stdout", "stderr"):
                        if name.endswith("/" + channel + ".txt"):
                            channels[channel] = summary
                    archive.add(path, arcname=name, recursive=False)
                end = stream.tell()
            stream.truncate()
            stream.flush()
            os.fsync(stream.fileno())
        saved = {**row, "output_archive": "task_outputs.tar", "archive_end": end,
                 "artifacts": artifacts}
        for channel in ("stdout", "stderr"):
            saved[channel + "_summary"] = channels.get(channel)
        with self.checkpoint.open("ab") as stream:
            stream.write((json.dumps(saved, ensure_ascii=False) + "\n").encode())
            stream.flush()
            os.fsync(stream.fileno())
        self.end = end
        self.rows.append(saved)
        self.ids.add(saved["id"])
        return saved


class ArchiveReader:
    """Scan headers once; read member bytes through file/gzip streams, no extraction."""

    def __init__(self, path):
        self.path = Path(path)
        self.members = {}
        try:
            with tarfile.open(self.path, "r|*", bufsize=1024 * 1024) as archive:
                for member in archive:
                    if member.isfile():
                        name = str(PurePosixPath(member.name))
                        self.members[name] = (member.offset_data, member.size)
                    archive.members.clear()
        except tarfile.TarError as exc:
            raise ValueError("unreadable archive: " + str(self.path)) from exc

    def read(self, name):
        name = str(PurePosixPath(name))
        if name not in self.members:
            raise FileNotFoundError("archive member missing: " + name)
        offset, size = self.members[name]
        opener = gzip.open if self.path.name.endswith((".gz", ".tgz")) else open
        with opener(self.path, "rb") as stream:
            stream.seek(offset)
            data = stream.read(size)
        if len(data) != size:
            raise ValueError("archive member truncated: " + name)
        return data


class MemoryFile:
    def __init__(self, data):
        self.data = data

    def read_text(self):
        return self.data.decode()

    def read_bytes(self):
        return self.data


class RunReader:
    """Directory, new task archive, or an archived legacy run (tar[.gz])."""

    def __init__(self, path):
        parts = str(path).split("::", 1)
        self.path = Path(parts[0])
        requested = parts[1].strip("/") if len(parts) == 2 else None
        self.archives = {}
        self.other_runs = {}
        self.root = None
        if self.path.is_file():
            archive = self.archive(self.path)
            candidates = [n for n in archive.members if n.endswith("/run.json") or n == "run.json"]
            if not candidates:
                raise ValueError("archive contains no run.json")
            shortest = min(len(PurePosixPath(n).parts) for n in candidates)
            candidates = [n for n in candidates if len(PurePosixPath(n).parts) == shortest]
            if len(candidates) != 1:
                raise ValueError("archive contains multiple runs; choose an individual run archive")
            self.root = str(PurePosixPath(candidates[0]).parent)
            if requested:
                candidate = str(PurePosixPath(requested) / "run.json")
                if candidate not in archive.members:
                    candidate = str(PurePosixPath(self.root) / requested / "run.json")
                self.root = str(PurePosixPath(candidate).parent)
            else:
                candidate = candidates[0]
            self.manifest = json.loads(archive.read(candidate))
            self.manifest_sha256 = hashlib.sha256(archive.read(candidate)).hexdigest()
            self._recover_legacy_archive()
        else:
            data = (self.path / "run.json").read_bytes()
            self.manifest = json.loads(data)
            self.manifest_sha256 = hashlib.sha256(data).hexdigest()
            if self.manifest.get("output_format") == FORMAT and self.manifest.get("run_kind") != "array":
                rows = jsonl_rows(self.path / "task_results.jsonl")
                completed = {r["id"]: r for r in rows}
                if len(rows) != len(completed):
                    raise ValueError("duplicate checkpoint job")
                schedule = self.manifest.get("planned_jobs", rows)
                if set(completed) - {j["id"] for j in schedule}:
                    raise ValueError("unplanned checkpoint job")
                self.manifest["jobs"] = [completed.get(j["id"], {**j, "directory": "jobs/" + j["id"],
                    "result_present": False, "status": "unfinished", "returncode": None,
                    "elapsed_seconds": None, "time_wall_seconds": None, "peak_rss_kib": None}) for j in schedule]

    def _recover_legacy_archive(self):
        archive = self.archive(self.path)
        root = PurePosixPath(self.root)
        plan_path = str(root / "array_plan.json")
        if self.manifest.get("run_kind") == "array" and plan_path in archive.members:
            plan = json.loads(archive.read(plan_path))
            jobs, statuses, tasks = [], [], []
            for task in plan["tasks"]:
                task_root = root / task["directory"]
                metadata = json.loads(archive.read(str(task_root / "task.json")))
                status = metadata["status"]
                if status == "running":
                    status = "interrupted"
                schedule = task.get("planned_jobs")
                if schedule is None:
                    schedule = json.loads(archive.read(str(task_root / "planned_jobs.json")))
                missing = 0
                for planned in schedule:
                    directory = str(PurePosixPath(task["directory"]) / "run/jobs" / planned["id"])
                    marker = str(root / directory / "result.json")
                    if marker in archive.members:
                        result = json.loads(archive.read(marker))
                        if any(result.get(k) != v for k, v in planned.items()):
                            raise ValueError("archived result differs from plan")
                        row = {**result, "result_present": True}
                    else:
                        missing += 1
                        row = {**planned, "result_present": False, "status": "task_failed" if status == "failed" else
                               "not_started" if status == "not_started" else "unfinished",
                               "returncode": None, "elapsed_seconds": None,
                               "time_wall_seconds": None, "peak_rss_kib": None}
                    row.update(directory=directory, tier=task["tier"], task_id=task["task_id"],
                               global_task_number=task["global_task_number"])
                    jobs.append(row)
                if missing and status == "success":
                    status = "interrupted"
                tasks.append({**metadata, "status": status, "missing_jobs": missing})
                statuses.append(status)
            from collections import Counter
            measured = [j for j in jobs if not j["warmup"]]
            self.manifest.update(jobs=jobs, summary_only=False, tasks=tasks,
                                 task_status_counts=dict(Counter(statuses)), planned_measured_jobs=len(measured),
                                 completed_measured_jobs=sum(j["result_present"] for j in measured),
                                 missing_measured_jobs=sum(not j["result_present"] for j in measured))
            self.manifest["status"] = "interrupted" if any(s not in ("success", "failed") for s in statuses) else \
                                      "failed" if "failed" in statuses else "success"
        elif self.manifest.get("output_format") != FORMAT:
            schedule = self.manifest.get("planned_jobs", self.manifest["jobs"])
            found = []
            for planned in schedule:
                directory = planned.get("directory") or "jobs/" + planned["id"]
                marker = str(root / directory / "result.json")
                if marker in archive.members:
                    result = json.loads(archive.read(marker))
                    found.append({**result, "result_present": True})
                elif "planned_jobs" in self.manifest:
                    found.append({**planned, "directory": directory, "result_present": False,
                                  "status": "unfinished", "returncode": None, "elapsed_seconds": None,
                                  "time_wall_seconds": None, "peak_rss_kib": None})
                else:
                    found.append(planned)
            self.manifest["jobs"] = found
            if self.manifest["status"] == "running":
                self.manifest["status"] = "interrupted"

    def archive(self, path):
        path = Path(path)
        if path not in self.archives:
            self.archives[path] = ArchiveReader(path)
        return self.archives[path]

    def read(self, relative):
        if self.root is not None:
            return self.archive(self.path).read(str(PurePosixPath(self.root) / relative))
        return (self.path / relative).read_bytes()

    def input(self, entry):
        if "input_archive" in entry:
            return self.archive(self.path / entry["input_archive"]).read(entry["input_member"])
        if self.root is not None:
            archive = self.archive(self.path)
            name = str(PurePosixPath(self.root) / "inputs" / (entry["name"] + ".txt"))
            if name not in archive.members:
                # Array archives store plan-level inputs; prefer them to vanished absolute paths.
                original = Path(entry["copied_path"])
                suffix = "/inputs/" + original.name
                matches = [n for n in archive.members if n.endswith(suffix)]
                if not matches:
                    raise FileNotFoundError("archived input missing: " + entry["name"])
                name = min(matches, key=lambda n: len(PurePosixPath(n).parts))
            return archive.read(name)
        return Path(entry.get("source_input_path", entry["copied_path"])).read_bytes()

    def artifact(self, job, name, verify=True):
        if "artifact_run" in job:
            path = job["artifact_run"]
            if path not in self.other_runs:
                self.other_runs[path] = RunReader(path)
            return self.other_runs[path].artifact({k: v for k, v in job.items() if k != "artifact_run"}, name, verify)
        member = str(PurePosixPath(job.get("archive_member_directory", job["directory"])) / name)
        if "output_archive" in job:
            data = self.archive(self.path / job["output_archive"]).read(member)
            expected = job.get("artifacts", {}).get(member, {}).get("sha256")
            # Aggregation changes directory but retains the member names from the task.
            if expected is None:
                expected = job.get("artifact_hashes", {}).get(name)
            if verify and expected and hashlib.sha256(data).hexdigest() != expected:
                raise ValueError("output artifact hash changed: " + member)
            return data
        if name == "bpp.tsv" and "--bpp" in job.get("command", []):
            exported = Path(job["command"][job["command"].index("--bpp") + 1])
            if self.root is None:
                return exported.read_bytes()
            try:
                return self.read(member)
            except FileNotFoundError:
                return self.read(str(PurePosixPath(job["directory"]) / exported.name))
        return self.read(member)


def input_records(reader):
    from rna_records import read_record
    records = {}
    for entry in reader.manifest["inputs"]:
        data = reader.input(entry)
        if hashlib.sha256(data).hexdigest() != entry["sha256"]:
            raise ValueError("copied input changed: " + entry["name"])
        records[entry["name"]] = read_record(MemoryFile(data))
    return records


def artifact_digest(reader, job, name):
    try:
        return hashlib.sha256(reader.artifact(job, name, verify=False)).hexdigest()
    except FileNotFoundError:
        return None
