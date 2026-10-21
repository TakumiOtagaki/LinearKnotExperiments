#!/usr/bin/env python3
"""Compare the renamed CLI with 89baea52 on the frozen 14 release inputs.

Only elapsed-time lines are removed. Raw outputs are kept in one tar; transient
inputs and per-process files stay in node-local scratch. Never submits jobs.
"""
import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import tarfile

from compare_release_outputs import compare_bytes
from run_cli import ROOT, digest, git_state, now, run_one, write_json
from task_storage import local_workspace

COMMANDS = (
    ("fold", "LinearCParty", ["--bpp-cutoff", "0"]),
    ("profile", "LinearPKCapR", []),
    ("access", "LinearPKRaccess", ["--access-len", "1,4,8,16", "--probability"]),
    ("access-byloop", "LinearPKRaccess", ["--access-len", "1,4,8,16", "--probability", "--by-loop"]),
    ("sample", "LinearPKSampling", ["--count", "100", "--seed", "20261016"]),
    ("sample-unique", "LinearPKSampling", ["--count", "100", "--seed", "20261016", "--unique"]),
    ("mea", "LinearPKMEA", ["--gamma", "1"]),
    ("centroid", "LinearPKCentroid", ["--gamma", "1"]),
    ("subopt", "LinearPKSubopt", ["--top-k", "8"]),
)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidate-build", type=Path, default=ROOT / "build/linearknot-9c7d6e1")
    p.add_argument("--reference-build", type=Path, default=ROOT / "build/linear-cparty-89baea52")
    p.add_argument("--candidate-source", type=Path, default=ROOT / "deps/LinearKnot-9c7d6e1")
    p.add_argument("--reference-source", type=Path, default=ROOT / "deps/LinearCParty-89baea52")
    p.add_argument("--reference-cli-layout", choices=("legacy", "unified"), default="legacy")
    p.add_argument("--scratch-dir", type=Path)
    p.add_argument("--timeout", type=float, default=180)
    p.add_argument("--memory-limit-mib", type=int, default=16384)
    p.add_argument("--out", type=Path, required=True, help="new report directory")
    a = p.parse_args()
    if a.memory_limit_mib < 64:
        p.error("memory-limit-mib must be at least 64")
    if a.out.exists():
        p.error("choose a fresh --out directory")
    receipt = ROOT / "manifests/release-preparation-89baea52.json"
    frozen = json.loads(receipt.read_text())["stages"]["comparison"]["specification"]["inputs"]
    assert len(frozen) == 14
    sources = {"reference": a.reference_source.resolve(), "candidate": a.candidate_source.resolve()}
    builds = {"reference": a.reference_build.resolve(), "candidate": a.candidate_build.resolve()}
    binaries = {"reference": ["linearknot"] if a.reference_cli_layout == "unified" else [v[1] for v in COMMANDS], "candidate": ["linearknot"]}
    fingerprints = {s: {"source": git_state(sources[s]), "binaries": {n: digest(builds[s] / n) for n in set(binaries[s])},
                         "cmake_cache_sha256": digest(builds[s] / "CMakeCache.txt"),
                         "parameter_sha256": digest(sources[s] / "params/rna_Turner04.par")} for s in sources}
    assert all(not f["source"]["dirty"] and all(f["binaries"].values()) for f in fingerprints.values())
    assert fingerprints["reference"]["parameter_sha256"] == fingerprints["candidate"]["parameter_sha256"]
    a.out.mkdir(parents=True)
    env = {k: v for k, v in os.environ.items() if not k.startswith("CPARTY_")}
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", LC_ALL="C")
    report = {"schema_version": 1, "started_utc": now(), "fingerprints": fingerprints,
              "inputs": frozen, "input_archive": "results/runs/e1-bprna-short-d21618c5-20261014/inputs.tar",
              "release_receipt_sha256": digest(receipt), "commands": COMMANDS,
              "ignored_lines": ["# elapsed_seconds"], "timeout_seconds": a.timeout,
              "memory_limit_mib": a.memory_limit_mib, "script_sha256": digest(Path(__file__)),
              "reference_cli_layout": a.reference_cli_layout,
              "planned_comparisons": 14 * 3 * len(COMMANDS), "rows": [], "status": "running"}
    write_json(a.out / "comparison.json", report)
    with local_workspace(a.scratch_dir) as work, tarfile.open(ROOT / report["input_archive"]) as inputs, \
            tarfile.open(a.out / "raw_outputs.tar", "w") as archive, (a.out / "processes.jsonl").open("w") as log:
        env["CPARTY_FIXED_ENERGY_PARAM_FILE"] = str(work / "energy.par")
        (work / "energy.par").write_bytes((sources["candidate"] / "params/rna_Turner04.par").read_bytes())
        members = {Path(m.name).name: m for m in inputs if m.isfile()}
        for item in frozen:
            payload = inputs.extractfile(members[item["name"] + ".txt"]).read()
            assert hashlib.sha256(payload).hexdigest() == item["sha256"]
            source = work / "input.txt"
            source.write_bytes(payload)
            for beam in (20, 200, None):
                mode = "exact" if beam is None else "beam" + str(beam)
                options = ["--exact"] if beam is None else ["-b", str(beam), "--max-hairpin", str(max(beam, 30)),
                            "--max-case3-depth", str(max(beam, 30)), "--ranking", "outer"]
                for label, legacy, extra in COMMANDS:
                    row = {"input": item["name"], "input_sha256": item["sha256"], "mode": mode, "command": label}
                    outputs = {}
                    for side in ("reference", "candidate"):
                        directory = work / side
                        directory.mkdir()
                        prefix = [str(builds[side] / legacy)] if side == "reference" and a.reference_cli_layout == "legacy" else [str(builds[side] / "linearknot"), label.split("-")[0]]
                        argv = [*prefix, "-i", str(source), *options, *extra]
                        if label == "fold":
                            argv += ["--bpp", str(directory / "bpp.tsv")]
                        result = run_one(argv, directory, env, a.timeout, memory_limit_mib=a.memory_limit_mib)
                        row[side] = {**result, "argv": argv}
                        log.write(json.dumps({"input": item["name"], "mode": mode, "command": label, "side": side, **row[side]}) + "\n")
                        log.flush()
                        outputs[side] = {}
                        for filename in ["stdout.txt", "stderr.txt"] + (["bpp.tsv"] if label == "fold" else []):
                            path = directory / filename
                            data = path.read_bytes() if path.is_file() else None
                            outputs[side][filename] = data
                            if data is not None:
                                info = tarfile.TarInfo("/".join((item["name"], mode, label, side, filename)))
                                info.size = len(data)
                                archive.addfile(info, io.BytesIO(data))
                        shutil.rmtree(directory)
                    row["files"] = {n: compare_bytes(outputs["reference"][n], outputs["candidate"][n])
                        if outputs["reference"][n] is not None and outputs["candidate"][n] is not None else {"identical": False, "missing": True}
                        for n in outputs["reference"]}
                    row["successful"] = all(row[s]["status"] == "success" for s in outputs)
                    row["identical"] = row["successful"] and all(f["identical"] for f in row["files"].values())
                    report["rows"].append(row)
            write_json(a.out / "comparison.json", report)
            print(item["name"], len(report["rows"]), "/", report["planned_comparisons"], flush=True)
    for side, saved in fingerprints.items():
        assert saved["source"] == git_state(sources[side])
        assert saved["binaries"] == {n: digest(builds[side] / n) for n in saved["binaries"]}
        assert saved["cmake_cache_sha256"] == digest(builds[side] / "CMakeCache.txt")
        assert saved["parameter_sha256"] == digest(sources[side] / "params/rna_Turner04.par")
    report.update(finished_utc=now(), status="complete", comparisons=len(report["rows"]),
                  successful=all(r["successful"] for r in report["rows"]),
                  identical=all(r["identical"] for r in report["rows"]),
                  identical_comparisons=sum(r["identical"] for r in report["rows"]),
                  artifacts_sha256={n: digest(a.out / n) for n in ("raw_outputs.tar", "processes.jsonl")})
    write_json(a.out / "comparison.json", report)
    print(json.dumps({k: report[k] for k in ("successful", "identical", "comparisons", "identical_comparisons")}))
    return 0 if report["successful"] and (report["identical"] or a.reference_cli_layout == "legacy") else 1


if __name__ == "__main__":
    raise SystemExit(main())
