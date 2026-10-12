#!/usr/bin/env python3
"""Generate E2 shape controls in four files, without per-record shared files.

All nonempty inputs have floor(0.15*n) pairs. Main families: empty,
long_stem, separated_helices, deep_branching. Additional matched contrasts
change gap allocation, helix count, branching topology or branching budget.
Embedding compatible pairs changes S: this is a joint (S,G) control.
"""

import argparse
import hashlib
import io
import json
from pathlib import Path
import random
import tarfile

from scaffold_features import scaffold_features, write_tsv

DEFAULT_LENGTHS = (250, 500, 1000, 2000, 4000)
DEFAULT_SEEDS = (2026101201, 2026101202, 2026101203)
DEPTHS = (15, 30, 31, 99, 100, 101)


def allocate(total, count, minimum=0, concentrated=False):
    if count < 1 or total < minimum * count:
        raise ValueError("insufficient positions")
    residual = total - minimum * count
    if concentrated:
        return [minimum + residual] + [minimum] * (count - 1)
    q, r = divmod(residual, count)
    return [minimum + q + int(i < r) for i in range(count)]


def binary_scaffold(branches, balanced=False):
    """Full binary tree: branches internal nodes and branches+1 hairpin leaves."""
    if branches == 0:
        return "(...)"
    if balanced:
        left = (branches - 1) // 2
        right = branches - 1 - left
    else:
        left, right = 0, branches - 1
    return "(" + binary_scaffold(left, balanced) + binary_scaffold(right, balanced) + ")"


def make_scaffold(n, family, branches=None):
    if n < 50:
        raise ValueError("synthetic lengths must be at least 50")
    m = 3 * n // 20
    if family == "empty":
        return "." * n
    if family == "long_stem":
        # Allocate free bases to exterior and hairpin, avoiding a giant fixed hairpin.
        free = n - 2 * m
        left, right = allocate(free - 3, 2)
        return "." * left + "(" * m + "..." + ")" * m + "." * right
    if family.startswith("separated"):
        count = 12 if family == "separated_12" else 3
        if m < count or n - 2 * m < 4 * count + 1:
            raise ValueError("insufficient pair/free budget for helix count")
        stems = allocate(m, count, 1)
        gaps = allocate(n - 2 * m - 3 * count, count + 1, 1,
                        family == "separated_concentrated")
        pieces = ["." * gaps[0]]
        for stem, gap in zip(stems, gaps[1:]):
            pieces.extend(["(" * stem + "..." + ")" * stem, "." * gap])
        return "".join(pieces)
    if family in ("deep_branching", "balanced_branching", "branch_depth"):
        k = min(101, (m - 1) // 2) if branches is None else branches
        if k < 1 or 2 * k + 1 > m:
            raise ValueError("branching budget exceeds fixed pair count")
        # Unary trunk pads pair count; topology comparison preserves its length,
        # branch count, leaves, and free-run histogram exactly.
        trunk = m - (2 * k + 1)
        body = "(" * trunk + binary_scaffold(k, family == "balanced_branching") + ")" * trunk
        slack = n - len(body)
        if slack < 0:
            raise ValueError("branching hairpins exceed length")
        left, right = allocate(slack, 2)
        return "." * left + body + "." * right
    raise ValueError("unknown family: " + family)


def embed_sequence(scaffold, seed):
    """Exactly round(n/2) GC bases; canonical AU/GC scaffold pairs only."""
    rng = random.Random(seed)
    stack, pairs, free = [], [], []
    for i, s in enumerate(scaffold):
        if s == "(":
            stack.append(i)
        elif s == ")":
            pairs.append((stack.pop(), i))
        else:
            free.append(i)
    target_gc = (len(scaffold) + 1) // 2
    gc_pairs = min(len(pairs), target_gc // 4)
    labels = [True] * gc_pairs + [False] * (len(pairs) - gc_pairs)
    rng.shuffle(labels)
    seq = [""] * len(scaffold)
    for (i, j), gc in zip(pairs, labels):
        a, b = rng.choice(("GC", "CG") if gc else ("AU", "UA"))
        seq[i], seq[j] = a, b
    free_gc = target_gc - 2 * gc_pairs
    if not 0 <= free_gc <= len(free):
        raise ValueError("GC target infeasible")
    labels = [True] * free_gc + [False] * (len(free) - free_gc)
    rng.shuffle(labels)
    for i, gc in zip(free, labels):
        seq[i] = rng.choice("GC" if gc else "AU")
    return "".join(seq)


def generate(lengths=DEFAULT_LENGTHS, seeds=DEFAULT_SEEDS):
    for n in lengths:
        max_depth = min(101, (3 * n // 20 - 1) // 2)
        variants = [(family, None) for family in (
            "empty", "long_stem", "separated_helices", "deep_branching",
            "separated_concentrated", "separated_12", "balanced_branching")]
        variants += [("branch_depth", d) for d in DEPTHS if d <= max_depth and d != max_depth]
        for seed in seeds:
            for family, branches in variants:
                suffix = f"_d{branches}" if branches is not None else ""
                name = f"syn_n{n}_s{seed}_{family}{suffix}"
                scaffold = make_scaffold(n, family, branches)
                # Same seed within n gives reproducible paired backgrounds;
                # different pair endpoints still change the sequence.
                sequence = embed_sequence(scaffold, seed + n)
                features = scaffold_features(sequence, scaffold)
                yield {"name": name, "family": family, "branch_budget": branches,
                       "seed": seed, "sequence_seed": seed + n,
                       "sequence": sequence, "scaffold": scaffold,
                       "features": features, "same_sequence_shape_only": False}


def save(output, lengths=DEFAULT_LENGTHS, seeds=DEFAULT_SEEDS):
    output = Path(output)
    if output.exists():
        raise ValueError("use a new immutable output directory")
    output.mkdir(parents=True)
    records = list(generate(lengths, seeds))
    jsonl, archive_path = output / "inputs.jsonl", output / "inputs.tar"
    entries = []
    with jsonl.open("w") as stream, tarfile.open(archive_path, "w") as archive:
        for row in records:
            stream.write(json.dumps(row, separators=(",", ":"), sort_keys=True) + "\n")
            data = f'>{row["name"]}\n{row["sequence"]}\n{row["scaffold"]}\n'.encode()
            member = "inputs/" + row["name"] + ".txt"
            info = tarfile.TarInfo(member)
            info.size, info.mtime, info.mode = len(data), 0, 0o644
            archive.addfile(info, io.BytesIO(data))
            entries.append({"name": row["name"], "path": str(output / member),
                            "input_archive": str(archive_path), "input_member": member,
                            "sha256": hashlib.sha256(data).hexdigest(),
                            "family": row["family"], "seed": row["seed"]})
    write_tsv(output / "features.tsv", ({"cohort": "synthetic", "name": r["name"], **r["features"]}
                                        for r in records))
    manifest = {"schema_version": 1, "generator_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "feature_script_sha256": hashlib.sha256(Path(__file__).with_name("scaffold_features.py").read_bytes()).hexdigest(),
                "lengths": list(lengths), "seeds": list(seeds), "records": len(records),
                "paired_fraction_policy": "2*floor(0.15*n)/n; empty has zero pairs",
                "gc_policy": "ceil(n/2) GC positions; AU/GC fixed pairs; record free-base composition",
                "sequence_caveat": "Compatible-pair embedding changes S between shapes; no same-S causal claim.",
                "inputs": entries,
                "files": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in (jsonl, archive_path, output / "features.tsv")}}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--lengths", default=",".join(map(str, DEFAULT_LENGTHS)))
    parser.add_argument("--seeds", default=",".join(map(str, DEFAULT_SEEDS)))
    args = parser.parse_args()
    manifest = save(args.out.resolve(), tuple(map(int, args.lengths.split(","))),
                    tuple(map(int, args.seeds.split(","))))
    print(json.dumps({"records": manifest["records"], "files": len(manifest["files"]) + 1}))


if __name__ == "__main__":
    main()
