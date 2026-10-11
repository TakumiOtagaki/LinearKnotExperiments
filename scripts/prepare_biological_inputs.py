#!/usr/bin/env python3
"""Freeze biological cohorts and construct thresholded LinearPartition scaffolds."""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import random
import resource
import subprocess
import sys
import time

from rna_records import describe_scaffold, read_record, scaffold_pairs
from run_cli import digest, git_state, now, stop_group, write_json

ROOT = Path(__file__).resolve().parents[1]


def read_fasta(path):
    records, name, chunks = [], None, []
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            header = line[1:].split()
            if not header:
                raise ValueError("empty FASTA header")
            name, chunks = header[0], []
        else:
            if name is None:
                raise ValueError("FASTA sequence before header")
            chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    if not records or len({name for name, _ in records}) != len(records):
        raise ValueError("empty FASTA or duplicate IDs")
    return records


def select_records(records, kind, seed):
    eligible, excluded, seen = [], [], set()
    for name, sequence in sorted(records):
        reason = None
        if not sequence or set(sequence) - set("ACGU"):
            reason = "non_ACGU"
        elif sequence in seen:
            reason = "duplicate_sequence"
        elif kind == "rnacentral" and not 4000 <= len(sequence) <= 30000:
            reason = "outside_initial_length_range"
        if reason:
            excluded.append({"id": name, "length": len(sequence), "reason": reason})
        else:
            seen.add(sequence)
            eligible.append((name, sequence))
    if kind == "rnacentral":
        return sorted(eligible, key=lambda r: (len(r[1]), r[0])), excluded, None
    rng, selected, bins = random.Random(seed), [], []
    for lo, hi in ((50, 100), (100, 200), (200, 400)):
        pool = [r for r in eligible if lo <= len(r[1]) < hi]
        if len(pool) < 4:
            raise ValueError("insufficient bpRNA candidates")
        chosen = sorted(rng.sample(pool, 4), key=lambda r: (len(r[1]), r[0]))
        selected.extend(chosen)
        bins.append({"min_inclusive": lo, "max_exclusive": hi,
                     "candidate_count": len(pool), "selected_ids": [r[0] for r in chosen]})
    return selected, excluded, bins


def threshold_scaffold(sequence, bpp_path, threshold, min_helix):
    pairs = {}
    for line in Path(bpp_path).read_text().splitlines():
        fields = line.split()
        if not fields or fields[0].startswith(("#", ">")):
            continue
        if len(fields) != 3:
            raise ValueError("malformed LinearPartition BPP")
        i, j, p = int(fields[0]), int(fields[1]), float(fields[2])
        if not 1 <= i < j <= len(sequence) or not math.isfinite(p) or not 0 <= p <= 1.00001:
            raise ValueError("invalid LinearPartition probability or coordinate")
        if p >= threshold:
            if sequence[i-1] + sequence[j-1] not in {"AU", "UA", "CG", "GC", "GU", "UG"}:
                raise ValueError("noncanonical threshold pair")
            if i in pairs:
                raise ValueError("duplicate threshold pair opening")
            pairs[i] = j
    positions = [x for pair in pairs.items() for x in pair]
    if len(set(positions)) != len(positions):
        raise ValueError("threshold pairs share a base")
    raw = ["."] * len(sequence)
    for i, j in pairs.items():
        raw[i-1], raw[j-1] = "(", ")"
    if scaffold_pairs("".join(raw)) != set(pairs.items()):
        raise ValueError("threshold pairs cross")
    retained = {}
    for i, j in pairs.items():
        if pairs.get(i-1) == j+1:
            continue
        run = [(i, j)]
        while pairs.get(run[-1][0]+1) == run[-1][1]-1:
            run.append((run[-1][0]+1, run[-1][1]-1))
        if len(run) >= min_helix:
            retained.update(run)
    result = ["."] * len(sequence)
    for i, j in retained.items():
        result[i-1], result[j-1] = "(", ")"
    return "".join(result)


def memory_cap(mib):
    def apply():
        limit = mib * 1024 * 1024
        resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    return apply


def prepare(source_manifest, output, lp, kind, seed, timeout_seconds=120, record_ids=None):
    source_manifest, output, lp = map(Path, (source_manifest, output, lp))
    source_info = json.loads(source_manifest.read_text())
    source = ROOT / source_info.get("copied_path", source_info.get("downloaded_path"))
    if digest(source) != source_info["sha256"]:
        raise ValueError("source checksum differs from acquisition manifest")
    state = git_state(lp.parent)
    if not state["commit"] or state["dirty"]:
        raise ValueError("LinearPartition needs a clean pinned checkout")
    selected, excluded, bins = select_records(read_fasta(source), kind, seed)
    cohort_ids = [r[0] for r in selected]
    if timeout_seconds <= 0 or not math.isfinite(timeout_seconds):
        raise ValueError("timeout_seconds must be finite and positive")
    if record_ids is not None:
        if len(set(record_ids)) != len(record_ids) or set(record_ids) - set(cohort_ids):
            raise ValueError("record IDs must be unique members of the selected cohort")
        selected = [r for r in selected if r[0] in record_ids]
    if not selected:
        raise ValueError("empty cohort")
    output.mkdir(parents=True, exist_ok=False)
    logs = output / "scaffold_logs"
    logs.mkdir()
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", LC_ALL="C")
    manifest = {"schema_version": 1, "created_utc": now(), "status": "running",
                "kind": kind, "source_sha256": digest(source),
                "acquisition_manifest_sha256": digest(source_manifest),
                "selection_seed": seed if kind == "bprna" else None,
                "selection_order": "ID-sorted pool; length/ID-sorted selected records",
                "candidate_bins": bins, "excluded_records": excluded,
                "selected_ids_before_scaffold": [r[0] for r in selected],
                "cohort_ids_before_subset": cohort_ids,
                "requested_record_ids": record_ids,
                "timeout_seconds": timeout_seconds, "memory_limit_mib": 8192,
                "environment": {key: env[key] for key in
                                ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "LC_ALL")},
                "scaffold_method": "LinearPartition-V posterior threshold, no MEA",
                "scaffold_settings": {"beam": 200, "dangles": 2,
                                      "pair_probability_threshold": 0.7, "min_helix": 2},
                "linearpartition_source": state, "wrapper_sha256": digest(lp),
                "vienna_binary_sha256": digest(lp.parent / "bin/linearpartition_v"),
                "script_sha256": digest(Path(__file__)), "records": []}
    write_json(output / "dataset.json", manifest)
    for name, sequence in selected:
        job = logs / name
        job.mkdir()
        bpp = job / "pairs.tsv"
        command = [str(lp.resolve()), "-V", "-b", "200", "--dangles", "2",
                   "-c", "0.7", "-r", str(bpp.resolve())]
        write_json(job / "command.json", {"argv": command, "started_utc": now(),
                                         "timeout_seconds": timeout_seconds, "memory_limit_mib": 8192})
        started = time.perf_counter()
        status, returncode, error = "success", None, None
        try:
            process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE, text=True, env=env,
                                       preexec_fn=memory_cap(8192), start_new_session=True)
            try:
                stdout, stderr = process.communicate(sequence+"\n", timeout=timeout_seconds)
            except subprocess.TimeoutExpired:
                stop_group(process)
                stdout, stderr = process.communicate()
                status = "timeout"
            (job / "stdout.txt").write_text(stdout)
            (job / "stderr.txt").write_text(stderr)
            returncode = process.returncode
            if status == "timeout" or returncode:
                raise ValueError("LinearPartition " + status + "; returncode=" + str(returncode))
            scaffold = threshold_scaffold(sequence, bpp, 0.7, 2)
            path = output / (name + ".txt")
            path.write_text(">"+name+"\n"+sequence+"\n"+scaffold+"\n")
            read_record(path)
            manifest["records"].append({"id": name, "input_file": path.name,
                "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                "input_sha256": digest(path), "scaffold_sha256": hashlib.sha256(scaffold.encode()).hexdigest(),
                **describe_scaffold(sequence, scaffold), "scaffold_status": "success"})
        except (ValueError, OSError) as exc:
            status, error = "failed" if status == "success" else status, str(exc)
            manifest["records"].append({"id": name, "length": len(sequence),
                                        "scaffold_status": "failed", "error": str(exc)})
        elapsed = time.perf_counter() - started
        manifest["records"][-1]["scaffold_elapsed_seconds"] = elapsed
        write_json(job / "result.json", {"status": status, "returncode": returncode,
                                         "error": error, "elapsed_seconds": elapsed})
        write_json(output / "dataset.json", manifest)
        print(name + ": " + manifest["records"][-1]["scaffold_status"], flush=True)
    manifest["linearpartition_source_after"] = git_state(lp.parent)
    manifest["vienna_binary_sha256_after"] = digest(lp.parent / "bin/linearpartition_v")
    manifest["wrapper_sha256_after"] = digest(lp)
    manifest["software_unchanged"] = (
        manifest["linearpartition_source_after"] == state
        and manifest["vienna_binary_sha256_after"] == manifest["vienna_binary_sha256"]
        and manifest["wrapper_sha256_after"] == manifest["wrapper_sha256"])
    manifest["status"] = "success" if manifest["software_unchanged"] and all(
        r["scaffold_status"] == "success" for r in manifest["records"]) else "failed"
    manifest["finished_utc"] = now()
    write_json(output / "dataset.json", manifest)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kind", choices=("bprna", "rnacentral"), required=True)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--linearpartition", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20261008)
    parser.add_argument("--timeout-seconds", type=float, default=120)
    parser.add_argument("--ids", nargs="+", help="prepare only these IDs from the unchanged cohort")
    args = parser.parse_args()
    result = prepare(args.source_manifest, args.out, args.linearpartition, args.kind, args.seed,
                     args.timeout_seconds, args.ids)
    return int(result["status"] != "success")


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ValueError, OSError) as exc:
        print("error: " + str(exc), file=sys.stderr)
        raise SystemExit(2)
