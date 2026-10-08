#!/usr/bin/env python3
"""Freeze a small synthetic cohort and development accuracy configuration."""

import argparse
import hashlib
import json
from pathlib import Path
import random

from rna_records import describe_scaffold, read_record


ROOT = Path(__file__).resolve().parents[1]


def dot_bracket(length, pairs):
    symbols = ["."] * length
    for i, j in pairs:
        symbols[i - 1], symbols[j - 1] = "(", ")"
    return "".join(symbols)


def prepare(destination, config_path, seed):
    if destination.exists() or config_path.exists():
        raise ValueError("use new dataset and config paths; outputs are immutable")
    destination.relative_to(ROOT)
    config_path.relative_to(ROOT)
    destination.mkdir(parents=True)
    config_path.parent.mkdir(parents=True, exist_ok=True)
    records = []
    for length in (40, 80, 160):
        outer = [(i, length + 1 - i) for i in (1, 2, 3)]
        left = [(7 + k, length // 2 - 3 - k) for k in range(3)]
        right = [(length // 2 + 4 + k, length - 6 - k) for k in range(3)]
        full = outer + left + right
        for replicate in range(2):
            sequence_seed = seed + 1000 * length + replicate
            rng = random.Random(sequence_seed)
            sequence = [rng.choice("ACGU") for _ in range(length)]
            for i, j in full:
                sequence[i - 1], sequence[j - 1] = rng.choice(
                    ("AU", "UA", "CG", "GC", "GU", "UG"))
            sequence = "".join(sequence)
            for kind, pairs in (("empty", []), ("single", left), ("branched", full)):
                name = f"n{length}_r{replicate}_{kind}"
                scaffold = dot_bracket(length, pairs)
                path = destination / (name + ".txt")
                path.write_text(f">{name}\n{sequence}\n{scaffold}\n")
                read_record(path)
                records.append({"name": name, "path": str(path.relative_to(ROOT)),
                                "sequence_group": f"n{length}_r{replicate}",
                                "sequence_seed": sequence_seed, "scaffold_kind": kind,
                                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                                "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                                **describe_scaffold(sequence, scaffold)})
    dataset = {"schema_version": 1, "purpose": "synthetic diagnostic",
               "master_seed": seed, "lengths": [40, 80, 160],
               "sequence_replicates": 2, "records": records,
               "generation": "Uniform ACGU backgrounds with nine canonical pairs imposed; "
               "three nested scaffold choices share each resulting sequence. "
               "These are supplied geometric controls, not predicted biological scaffolds."}
    (destination / "dataset.json").write_text(json.dumps(dataset, indent=2) + "\n")
    modes = [{"name": "exact", "args": ["--exact"]}]
    # Both caps exceed every possible length/depth in this bounded cohort.
    for beam in (20, 50, 100, 200):
        modes.append({"name": f"beam{beam}_nocap", "args":
                      ["-b", str(beam), "--max-hairpin", "1000000",
                       "--max-case3-depth", "1000000"]})
    modes.append({"name": "beam100_default", "args":
                  ["-b", "100", "--max-hairpin", "30", "--max-case3-depth", "30"]})
    config = {"schema_version": 1, "name": "e1-synthetic-pilot",
              "purpose": "development", "energy_parameter_file": "params/rna_Turner04.par",
              "model_settings": {"temperature_celsius": 37,
                                 "temperature_source": "CLI/API default; verify RT metadata",
                                 "build_revision_verified": False},
              "environment": {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                              "MKL_NUM_THREADS": "1", "LC_ALL": "C"},
              "timeout_seconds": 60, "warmups": 0, "repetitions": 1,
              "inputs": [{"name": r["name"], "path": r["path"]} for r in records],
              "modes": modes,
              "programs": [
                  {"name": "LinearCParty", "args": ["--bpp-cutoff", "0"],
                   "outputs": [{"name": "bpp.tsv", "flag": "--bpp"}]},
                  {"name": "LinearPKCapR", "args": []},
                  {"name": "LinearPKRaccess", "args":
                   ["--access-len", "1,4,8,16", "--probability"]}],
              "analysis": {"reference_mode": "exact",
                           "default_cap_comparison": "Finite-beam setting sensitivity; "
                           "not an isolated cap-only error decomposition."}}
    config_path.write_text(json.dumps(config, indent=2) + "\n")
    return len(records), len(records) * len(modes) * len(config["programs"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="data/pilot-v1")
    parser.add_argument("--config-out", default="configs/e1_pilot.json")
    parser.add_argument("--seed", type=int, default=20261008)
    args = parser.parse_args()
    count, jobs = prepare((ROOT / args.out).resolve(), (ROOT / args.config_out).resolve(), args.seed)
    print(f"Prepared {count} records and {jobs} jobs; no CLI jobs were executed.")


if __name__ == "__main__":
    main()
