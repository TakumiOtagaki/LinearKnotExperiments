#!/usr/bin/env python3
"""Run the prepared E0 decoding suite with centroid; never submits jobs."""
import fcntl
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import time

from prepare_paper_rerun_89baea52 import ROOT, SOURCE, release_ready
from run_cli import digest, now, write_json


def main():
    release_ready()
    preparation = ROOT / "docs/centroid_89baea52_20261020/preparation.json"
    frozen = json.loads(preparation.read_text())
    def verify():
        for name, expected in frozen["files_sha256"].items():
            if digest(ROOT / name) != expected:
                raise ValueError("prepared dependency changed: " + name)
    verify()
    scratch = Path(os.environ["TMPDIR"]).resolve()
    if scratch == Path("/gs/bs") or str(scratch).startswith("/gs/bs/"):
        raise ValueError("TMPDIR must be node-local")
    for name in list(os.environ):
        if name.startswith("CPARTY_"):
            del os.environ[name]
    os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", LC_ALL="C",
        PYTHONDONTWRITEBYTECODE="1", CPARTY_PF_OVERRIDE_SCALE="1",
        CPARTY_FIXED_ENERGY_PARAM_FILE=str(SOURCE / "params/rna_Turner04.par"))
    tool = ROOT / "tools/validation_centroid_89baea52_20261020/decoding"
    out = ROOT / "results/local/validation-controller-decoding-centroid-89baea52-20261020"
    lock_path = out.parent / ".validation-centroid-89baea52.lock"
    with lock_path.open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        out.mkdir(exist_ok=False)
        commands = [[sys.executable, "-B", str(tool / n)] for n in ("run.py", "analyze.py", "validate.py")]
        meta = dict(started_utc=now(), status="running", commands=commands,
                    preparation_sha256=digest(preparation), job_id=os.environ.get("JOB_ID"), sge_task_id=os.environ.get("SGE_TASK_ID"))
        write_json(out / "provenance.json", meta)
        results = []
        with tempfile.TemporaryDirectory(prefix="centroid-e0-", dir=scratch) as td, \
                tarfile.open(out / "task_outputs.tar", "w") as archive, (out / "task_results.jsonl").open("w") as log:
            for i, command in enumerate(commands):
                paths = [Path(td) / n for n in ("stdout", "stderr")]
                start = time.monotonic()
                with paths[0].open("w") as stdout, paths[1].open("w") as stderr:
                    result = subprocess.run(command, stdout=stdout, stderr=stderr, cwd=ROOT)
                row = dict(argv=command, returncode=result.returncode, elapsed_seconds=time.monotonic()-start,
                           stdout_sha256=digest(paths[0]), stderr_sha256=digest(paths[1]))
                for path in paths:
                    archive.add(path, arcname=str(i) + "." + path.name)
                log.write(json.dumps(row) + "\n"); log.flush(); os.fsync(log.fileno())
                results.append(row)
                if result.returncode:
                    meta.update(status="failed", finished_utc=now(), results=results)
                    write_json(out / "provenance.json", meta)
                    raise RuntimeError("E0 validation failed; inspect packed stdout/stderr")
        release_ready(); verify()
        meta.update(status="success", finished_utc=now(), results=results)
        write_json(out / "provenance.json", meta)


if __name__ == "__main__":
    main()
