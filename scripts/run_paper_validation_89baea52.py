#!/usr/bin/env python3
"""Execute one prepared small validation, locally or as a one-task array."""
import argparse
import fcntl
import io
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tarfile
import tempfile
import time

from prepare_paper_rerun_89baea52 import BUILD, COMMIT, ROOT, SOURCE, release_ready
from run_cli import digest, now, write_json

KINDS = ("enumeration", "decoding", "byloop", "vienna", "linearspace")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("kind", choices=KINDS)
    args = p.parse_args()
    release_ready()
    scratch_root = Path(os.environ["TMPDIR"]).resolve()
    if str(scratch_root).startswith("/gs/bs/"):
        p.error("TMPDIR must be node-local")
    for key in list(os.environ):
        if key.startswith("CPARTY_"):
            del os.environ[key]
    os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", LC_ALL="C",
                      PYTHONDONTWRITEBYTECODE="1", CPARTY_PF_OVERRIDE_SCALE="1",
                      CPARTY_FIXED_ENERGY_PARAM_FILE=str(SOURCE / "params/rna_Turner04.par"))
    local = ROOT / "results/local"
    lock = (local / (".validation-89baea52-"+args.kind+".lock")).open("w")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    tool = ROOT / "tools/validation_89baea52" / args.kind
    if args.kind == "enumeration":
        commands = [[sys.executable, "-B", str(tool / name)] for name in ("run.py", "validate.py")]
    elif args.kind == "decoding":
        commands = [[sys.executable, "-B", str(tool / name)] for name in ("run.py", "analyze.py", "validate.py")]
    elif args.kind in ("byloop", "vienna"):
        script = "validate_byloop_89baea52.py" if args.kind == "byloop" else "validate_vienna_89baea52.py"
        commands = [[sys.executable, "-B", str(ROOT / "scripts" / script)]]
    else:
        commands = [[str(BUILD / "scfg2_inside_pair_max_test")]]
    report = local / ("validation-controller-"+args.kind+"-89baea52-20261018")
    report.mkdir(exist_ok=False)
    meta = {"started_utc": now(), "release": COMMIT, "commands": commands,
            "controller_sha256": digest(Path(__file__)), "job_id": os.environ.get("JOB_ID"),
            "sge_task_id": os.environ.get("SGE_TASK_ID"), "status": "running"}
    write_json(report / "provenance.json", meta)
    results = []
    with tempfile.TemporaryDirectory(prefix="validation-89baea52-", dir=scratch_root) as td, \
            tarfile.open(report / "task_outputs.tar", "w") as archive, \
            (report / "task_results.jsonl").open("w") as checkpoint:
        for i, command in enumerate(commands):
            out, err = Path(td)/"stdout", Path(td)/"stderr"
            print(shlex.join(command), flush=True)
            start = time.monotonic()
            with out.open("w") as fo, err.open("w") as fe:
                run = subprocess.run(command, stdout=fo, stderr=fe, cwd=ROOT.parent)
            result = {"argv": command, "returncode": run.returncode, "elapsed_seconds": time.monotonic()-start,
                      "stdout_sha256": digest(out), "stderr_sha256": digest(err)}
            for path, suffix in ((out, "stdout"), (err, "stderr")):
                archive.add(path, arcname=f"{i}.{suffix}")
            checkpoint.write(json.dumps(result)+"\n"); checkpoint.flush(); os.fsync(checkpoint.fileno())
            results.append(result)
            if run.returncode:
                meta.update(status="failed", finished_utc=now(), results=results)
                write_json(report / "provenance.json", meta)
                raise RuntimeError("validation failed; inspect its packed stdout/stderr")
    release_ready()
    meta.update(status="success", finished_utc=now(), results=results)
    write_json(report / "provenance.json", meta)


if __name__ == "__main__":
    main()
