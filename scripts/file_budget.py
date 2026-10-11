"""Reject plans that consume too much of a shared project's remaining inodes."""

import math
import subprocess


def query_quota(project, filesystem):
    try:
        result = subprocess.run(["lfs", "quota", "-p", str(project), filesystem],
                                capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.TimeoutExpired):
        return None
    if result.returncode:
        return None
    lines = result.stdout.splitlines()
    for index, line in enumerate(lines):
        fields = line.split()
        if fields and fields[0] == filesystem:
            if len(fields) == 1 and index + 1 < len(lines):
                fields += lines[index + 1].split()
            try:
                # Lustre: filesystem kbytes quota limit grace files quota limit grace.
                used = int(fields[5].rstrip("*"))
                limits = [int(fields[i].rstrip("*")) for i in (6, 7)]
                limit = min(v for v in limits if v > 0)
            except (ValueError, IndexError):
                return None
            return {"used_files": used, "limit_files": limit, "raw": result.stdout.strip()}
    return None


def estimate_files(tasks, tiers, compact_plan=True, scaffold=False):
    # Include directories, scheduler stdout per task, and one simultaneous atomic
    # JSON temporary per task. No entry depends on the number of CLI jobs/RNAs.
    per_task = 9 if scaffold else 12 + int(compact_plan)
    fixed = (16 if scaffold else 14) + 2 * tiers
    return {"format": "task-jsonl-tar-v1", "tasks": tasks, "nonempty_tiers": tiers,
            "entries_per_task": per_task, "fixed_entries": fixed,
            "estimated_shared_entries": fixed + tasks * per_task,
            "counts_directories": True, "scheduler_log_attempts": 1,
            "basis": "plan files, directories, task lock/metadata, JSONL, tar, run metadata, one atomic-write temporary, one scheduler log per task; retries add scheduler logs"}


def check_budget(estimate, settings):
    config = settings.get("file_budget", {})
    project, filesystem = config.get("project", 102571), config.get("filesystem", "/gs/bs")
    fraction = config.get("max_remaining_fraction", .1)
    cap = config.get("hard_file_limit", 24000000)
    if (type(fraction) not in (int, float) or not math.isfinite(fraction) or not 0 < fraction <= .1
            or type(cap) is not int or cap <= 0):
        raise ValueError("file_budget needs a positive hard limit and fraction in (0, 0.1]")
    quota = query_quota(project, filesystem)
    if quota is not None:
        remaining = max(0, min(cap, quota["limit_files"]) - quota["used_files"])
        source = "lfs quota (with configured conservative hard limit)"
    else:
        remaining = config.get("fallback_remaining_files")
        if type(remaining) is not int or remaining < 0:
            raise ValueError("cannot read project quota: configure file_budget.fallback_remaining_files")
        source = "config fallback_remaining_files; recheck live quota before submission"
    report = {**estimate, "project": project, "filesystem": filesystem, "quota": quota,
              "hard_file_limit": cap, "remaining_files": remaining, "remaining_source": source,
              "max_remaining_fraction": fraction, "allowed_entries": math.floor(remaining * fraction)}
    report["accepted"] = estimate["estimated_shared_entries"] <= report["allowed_entries"]
    if not report["accepted"]:
        raise ValueError("file-count budget refused: {} entries > {} allowed (remaining {}, fraction {})".format(
            estimate["estimated_shared_entries"], report["allowed_entries"], remaining, fraction))
    return report
