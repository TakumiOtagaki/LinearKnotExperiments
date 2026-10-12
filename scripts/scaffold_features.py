#!/usr/bin/env python3
"""Extract scaffold geometry in O(n) passes; never enumerate endpoint products.

Only ACGU and ordinary .() scaffolds are supported. A loop owner is the
innermost enclosing pair (0 for exterior). CASE3 counts *proper* branching
ancestors of that owner, matching scaffold_branching_paths.hh, not all
branching nodes enclosing the free position. Chains are G-only proxies,
not CParty bands in an unknown completed structure R.
"""

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import sys

from rna_records import read_record
from task_storage import MemoryFile, RunReader, input_records

CANONICAL = {"AU", "UA", "GC", "CG", "GU", "UG"}


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def histogram(values):
    return dict(sorted(Counter(values).items()))


def quantile(values, fraction):
    if not values:
        return 0
    values = sorted(values)
    x = (len(values) - 1) * fraction
    low = int(x)
    return values[low] + (values[min(low + 1, len(values) - 1)] - values[low]) * (x - low)


def canonical_count(counts):
    return counts["A"] * counts["U"] + counts["G"] * counts["C"] + counts["G"] * counts["U"]


def scaffold_features(sequence, scaffold, beam=100, hairpin=100, case3_depth=100, turn=3):
    if not sequence or set(sequence) - set("ACGU") or len(sequence) != len(scaffold):
        raise ValueError("equal-length nonempty ACGU sequence/scaffold required")
    if set(scaffold) - set(".()"):
        raise ValueError("features require ordinary .() scaffold; x/X/_ are not free bases")
    if any(type(x) is not int or x < 0 for x in (beam, hairpin, case3_depth, turn)) or beam == 0:
        raise ValueError("finite positive beam and nonnegative integer local caps required")
    n = len(sequence)
    children, parent, close, depth = {0: []}, {}, {}, {0: 0}
    owners, stack = [0] * n, [0]
    direct_free, loop_bases, runs, run = Counter(), {}, [], 0
    for pos, (base, symbol) in enumerate(zip(sequence, scaffold), 1):
        if symbol != "." and run:
            runs.append(run)
            run = 0
        if symbol == "(":
            p = stack[-1]
            children[p].append(pos)
            children[pos] = []
            parent[pos], depth[pos] = p, depth[p] + 1
            stack.append(pos)
        elif symbol == ")":
            if len(stack) == 1:
                raise ValueError("unbalanced scaffold")
            opening = stack.pop()
            close[opening] = pos
            if sequence[opening - 1] + base not in CANONICAL:
                raise ValueError("scaffold contains a noncanonical pair")
        else:
            owner = stack[-1]
            owners[pos - 1] = owner
            direct_free[owner] += 1
            loop_bases.setdefault(owner, Counter())[base] += 1
            run += 1
    if len(stack) != 1:
        raise ValueError("unbalanced scaffold")
    if run:
        runs.append(run)
    pair_nodes = list(parent)
    proper_branch = {0: 0}
    for node in pair_nodes:
        p = parent[node]
        proper_branch[node] = proper_branch[p] + int(len(children[p]) >= 2)
    # CASE3 uses paths of the containing pair to its proper ancestors.
    ancestors = [proper_branch[owners[i]] for i, s in enumerate(scaffold) if s == "."]
    subtree_pairs, subtree_free = {}, dict(direct_free)
    sibling_cross = 0
    for node in reversed([0] + pair_nodes):
        subtotal = 0
        for child in children[node]:
            weight = subtree_free[child]
            sibling_cross += subtotal * weight
            subtotal += weight
        subtree_free[node] = direct_free[node] + subtotal
        subtree_pairs[node] = int(node != 0) + sum(subtree_pairs[c] for c in children[node])
    f, m = sum(direct_free.values()), len(pair_nodes)
    cross = f * (f - 1) // 2 - sum(x * (x - 1) // 2 for x in direct_free.values())
    total_bases = sum(loop_bases.values(), Counter())
    cross_canonical = canonical_count(total_bases) - sum(canonical_count(c) for c in loop_bases.values())
    cross_turn, cross_canonical_turn = cross, cross_canonical
    # Only O(n*turn) nearby pairs, with fixed turn=3; no quadratic pair list.
    for j, symbol in enumerate(scaffold):
        if symbol != ".":
            continue
        for i in range(max(0, j - turn), j):
            if scaffold[i] == "." and owners[i] != owners[j]:
                cross_turn -= 1
                cross_canonical_turn -= sequence[i] + sequence[j] in CANONICAL
    # Count hairpin endpoint geometries within each contiguous free run.
    hp = 0
    for length in runs:
        high = min(hairpin + 1, length - 1)
        low = turn + 1
        if high >= low:
            count = high - low + 1
            hp += count * length - (low + high) * count // 2
    # Same-loop V endpoint supply per column; excluded close neighbours.
    seen, v_cap, same_loop_turn = Counter(), 0, 0
    for j, symbol in enumerate(scaffold):
        if j > turn and scaffold[j - turn - 1] == ".":
            seen[owners[j - turn - 1]] += 1
        if symbol == ".":
            eligible = seen[owners[j]]
            same_loop_turn += eligible
            v_cap += min(beam, eligible)
        elif symbol == ")":
            v_cap += 1
    unary_chains, stack_helices, left_arms, right_arms = [], [], [], []
    for node in pair_nodes:
        p = parent[node]
        if p == 0 or len(children[p]) != 1:
            end, count = node, 1
            while len(children[end]) == 1:
                end = children[end][0]
                count += 1
            unary_chains.append(count)
            left_arms.append(end - node + 1)
            right_arms.append(close[node] - close[end] + 1)
        if p == 0 or not (node == p + 1 and close[node] == close[p] - 1):
            end, count = node, 1
            while len(children[end]) == 1:
                c = children[end][0]
                if c != end + 1 or close[c] != close[end] - 1:
                    break
                end, count = c, count + 1
            stack_helices.append(count)
    branches = sum(len(c) >= 2 for c in children.values())
    return {
        "n": n, "scaffold_pairs": m, "free_bases": f,
        "paired_fraction": 2 * m / n, "gc_fraction": sum(b in "GC" for b in sequence) / n,
        "free_runs": len(runs), "free_run_max": max(runs, default=0),
        "free_run_mean": f / len(runs) if runs else 0,
        "free_run_median": quantile(runs, .5), "free_run_p90": quantile(runs, .9),
        "free_run_length_hist": histogram(runs), "nonempty_free_loops": len(direct_free),
        "loop_free_max": max(direct_free.values(), default=0),
        "loop_free_length_hist": histogram(direct_free.values()),
        "pair_depth_max": max(depth.values()), "pair_depth_sum": sum(depth[x] for x in pair_nodes),
        "B_G": sum(subtree_pairs[x] for x in pair_nodes),
        "branch_points_including_root": branches,
        "branch_points_excluding_root": branches - int(len(children[0]) >= 2),
        "max_pair_children": max(map(len, children.values())),
        "case3_ancestors_max": max(ancestors, default=0),
        "case3_ancestors_mean": sum(ancestors) / f if f else 0,
        "case3_ancestor_hist": histogram(ancestors),
        "case3_path_entries": sum(min(case3_depth, x) for x in ancestors),
        "case3_truncated_endpoints": sum(x > case3_depth for x in ancestors),
        "case3_sibling_endpoint_pairs": sibling_cross,
        "cross_loop_pairs": cross, "cross_loop_pairs_turn": cross_turn,
        "cross_loop_canonical_pairs": cross_canonical,
        "cross_loop_canonical_pairs_turn": cross_canonical_turn,
        "same_loop_pairs_turn": same_loop_turn,
        "hairpin_free_geometries": hp,
        "v_geometry_beam_cap_sum": v_cap,
        "be_geometry_beam_cap_sum": sum(min(beam, subtree_pairs[x]) for x in pair_nodes),
        "unary_chains": len(unary_chains), "unary_chain_pairs_hist": histogram(unary_chains),
        "unary_chain_left_arm_hist": histogram(left_arms), "unary_chain_right_arm_hist": histogram(right_arms),
        "unary_chain_arm_max": max(left_arms + right_arms, default=0),
        "stack_helices": len(stack_helices), "stack_helix_pairs_hist": histogram(stack_helices),
        "stack_helix_pairs_max": max(stack_helices, default=0),
        "beam": beam, "hairpin_cap": hairpin, "case3_depth_cap": case3_depth, "turn": turn,
        "sequence_sha256": sha256(sequence.encode()), "scaffold_sha256": sha256(scaffold.encode()),
    }


def write_tsv(path, records):
    with Path(path).open("w", newline="") as stream:
        writer = None
        for row in records:
            if writer is None:
                writer = csv.DictWriter(stream, fieldnames=list(row), delimiter="\t")
                writer.writeheader()
            writer.writerow({k: json.dumps(v, separators=(",", ":")) if isinstance(v, (dict, list)) else v
                             for k, v in row.items()})


def records_from_sources(rnacentral_run=None, bprna_dataset=None, jsonl=None):
    if rnacentral_run:
        reader = RunReader(rnacentral_run)
        for name, (_, sequence, scaffold) in input_records(reader).items():
            yield "rnacentral", name, sequence, scaffold
    if bprna_dataset:
        dataset = Path(bprna_dataset)
        manifest = json.loads(dataset.read_text())
        for entry in manifest["records"]:
            if entry["scaffold_status"] != "success":
                raise ValueError("unavailable scaffold: " + entry["id"])
            data = (dataset.parent / entry["input_file"]).read_bytes()
            if sha256(data) != entry["input_sha256"]:
                raise ValueError("input hash mismatch: " + entry["id"])
            name, sequence, scaffold = read_record(MemoryFile(data))
            if name != entry["id"] or sha256(sequence.encode()) != entry["sequence_sha256"]:
                raise ValueError("sequence identity mismatch: " + entry["id"])
            yield "bprna", name, sequence, scaffold
    if jsonl:
        with Path(jsonl).open() as stream:
            for line in stream:
                entry = json.loads(line)
                yield "synthetic", entry["name"], entry["sequence"], entry["scaffold"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rnacentral-run")
    parser.add_argument("--bprna-dataset")
    parser.add_argument("--jsonl")
    parser.add_argument("--out", required=True)
    parser.add_argument("--beam", type=int, default=100)
    parser.add_argument("--hairpin", type=int, default=100)
    parser.add_argument("--case3-depth", type=int, default=100)
    args = parser.parse_args()
    if not any((args.rnacentral_run, args.bprna_dataset, args.jsonl)):
        parser.error("at least one input source required")
    count = 0

    def rows():
        nonlocal count
        for cohort, name, sequence, scaffold in records_from_sources(
                args.rnacentral_run, args.bprna_dataset, args.jsonl):
            count += 1
            yield {"cohort": cohort, "name": name, **scaffold_features(
                sequence, scaffold, args.beam, args.hairpin, args.case3_depth)}
    write_tsv(args.out, rows())
    print(json.dumps({"records": count, "out": args.out}))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError) as exc:
        print("error: " + str(exc), file=sys.stderr)
        sys.exit(2)
