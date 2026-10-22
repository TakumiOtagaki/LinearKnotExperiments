#!/usr/bin/env python3
"""Count d-limited scaffold branching paths using the existing E2 definition.

Read frozen E1 long inputs and saved NEAT1 inputs without running RNA software.
Also reproduce the seven saved RNAcentral feature rows as a reference check.
"""

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess

from rna_records import describe_scaffold, read_record
from scaffold_features import scaffold_features
from task_storage import ArchiveReader, MemoryFile

ROOT = Path(__file__).resolve().parents[1]
LONG_RUN = Path("results/runs/e1-bprna-long-89baea52-cpu16-20261018")
RNA_RUN = Path("results/runs/e2-rnacentral-e2cecbb4-20261016")
COHORT = Path("manifests/cohorts/bprna-400plus-20261012.json")
NEAT1 = Path("data/neat1-20261016/manifest.json")
RNA_FEATURES = Path("results/local/e2-e2cecbb4-g-features-20261016/rnacentral_g_features.tsv")
CAP = 30


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def frozen_records(run, sources):
    plan_path, archive_path = run / "array_plan.json", run / "inputs.tar"
    plan = json.loads(plan_path.read_text())
    archive_hash = digest(archive_path)
    require(archive_hash == plan["inputs_archive_sha256"], "frozen archive changed")
    sources[str(run.relative_to(ROOT))] = {
        "array_plan_sha256": digest(plan_path), "inputs_archive_sha256": archive_hash,
        "inputs": len(plan["inputs"]),
    }
    archive = ArchiveReader(archive_path)
    seen = set()
    for entry in plan["inputs"]:
        data = archive.read(entry["input_member"])
        name, sequence, scaffold = read_record(MemoryFile(data))
        require(name == entry["name"] and name not in seen, "input identity/duplicate")
        require(sha256(data) == entry["sha256"], "frozen input hash mismatch")
        require(sha256(sequence.encode()) == entry["sequence_sha256"], "sequence hash mismatch")
        require(len(sequence) == entry["length"], "frozen input length mismatch")
        seen.add(name)
        yield name, sequence, scaffold, sha256(data)


def ancestor_walk_histogram(scaffold):
    """Independent check: explicitly walk parents from each free base's owner."""
    parent, children, owners, stack = {}, Counter(), [], [0]
    for pos, symbol in enumerate(scaffold, 1):
        if symbol == "(":
            parent[pos] = stack[-1]
            children[stack[-1]] += 1
            stack.append(pos)
        elif symbol == ")":
            stack.pop()
        else:
            owners.append(stack[-1])
    histogram = Counter()
    for owner in owners:
        count = 0
        while owner:
            owner = parent[owner]
            count += children[owner] >= 2
        histogram[count] += 1
    return dict(histogram)


def measure(dataset, name, sequence, scaffold, input_hash):
    features = scaffold_features(sequence, scaffold, case3_depth=CAP)
    histogram = ancestor_walk_histogram(scaffold)
    require(histogram == features["case3_ancestor_hist"], "independent ancestor histogram mismatch")
    require(sum(histogram.values()) == features["free_bases"], "free-base denominator mismatch")
    require(sum(count for depth, count in histogram.items() if depth > CAP)
            == features["case3_truncated_endpoints"], "truncation count mismatch")
    old = describe_scaffold(sequence, scaffold)
    return {
        "dataset": dataset, "input": name,
        **{key: features[key] for key in (
            "n", "scaffold_pairs", "free_bases", "case3_ancestors_max",
            "case3_depth_cap", "case3_truncated_endpoints")},
        "maximum_branching_ancestors_of_free_base": old["maximum_branching_ancestors_of_free_base"],
        "case3_ancestor_hist": features["case3_ancestor_hist"],
        "input_sha256": input_hash, "sequence_sha256": features["sequence_sha256"],
        "scaffold_sha256": features["scaffold_sha256"],
    }


def write_tsv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, separators=(",", ":"))
                             if isinstance(value, dict) else value for key, value in row.items()})


def analyze(output):
    sources, rows = {}, []
    cohort_path = ROOT / COHORT
    cohort = json.loads(cohort_path.read_text())
    expected = {entry["id"]: entry for entry in cohort["records"]}
    require(len(expected) == len(cohort["records"]) == 1438, "expected 1,438 distinct E1 long RNAs")
    dataset_path = ROOT / cohort["scaffold_dataset"]
    dataset = json.loads(dataset_path.read_text())
    scaffolds = {entry["id"]: entry for entry in dataset["records"]}
    plan = json.loads((ROOT / LONG_RUN / "array_plan.json").read_text())
    require(digest(cohort_path) == plan["cohort_provenance"]["cohort_sha256"], "cohort changed")
    require(digest(dataset_path) == plan["cohort_provenance"]["scaffold_dataset_sha256"], "dataset changed")
    for name, sequence, scaffold, input_hash in frozen_records(ROOT / LONG_RUN, sources):
        entry, saved = expected[name], scaffolds[name]
        require(400 <= len(sequence) == entry["length"], "not an E1 long input")
        require(sha256(sequence.encode()) == entry["sequence_sha256"], "cohort sequence mismatch")
        require(saved["scaffold_status"] == "success" and input_hash == saved["input_sha256"],
                "frozen input differs from published scaffold dataset")
        row = measure("bpRNA_ge400", name, sequence, scaffold, input_hash)
        require(row["scaffold_sha256"] == saved["scaffold_sha256"], "dataset scaffold mismatch")
        require(row["maximum_branching_ancestors_of_free_base"]
                == saved["maximum_branching_ancestors_of_free_base"], "E1 feature mismatch")
        rows.append(row)
    require({row["input"] for row in rows} == set(expected), "E1 long cohort incomplete")

    neat_manifest = json.loads((ROOT / NEAT1).read_text())
    require(neat_manifest["scaffold_settings"]["pair_probability_threshold"] == 0.7,
            "expected threshold 0.7")
    require({entry["name"] for entry in neat_manifest["records"]} == {"NEAT1_1", "NEAT1_2"}
            and len(neat_manifest["records"]) == 2, "expected two NEAT1 isoforms")
    for entry in neat_manifest["records"]:
        path = ROOT / NEAT1.parent / entry["path"]
        name, sequence, scaffold = read_record(path)
        row = measure(name, name, sequence, scaffold, digest(path))
        require(name == entry["name"] and len(sequence) == entry["length"], "NEAT1 identity mismatch")
        for key in ("sequence_sha256", "scaffold_sha256", "scaffold_pairs"):
            require(row[key] == entry[key], "NEAT1 manifest mismatch: " + key)
        rows.append(row)

    with (ROOT / RNA_FEATURES).open() as stream:
        saved_rows = list(csv.DictReader(stream, delimiter="\t"))
    reference = {row["name"]: row for row in saved_rows}
    require(len(reference) == len(saved_rows) == 7, "expected seven distinct RNAcentral RNAs")
    reference_checks = []
    for name, sequence, scaffold, input_hash in frozen_records(ROOT / RNA_RUN, sources):
        row = measure("rnacentral_reference", name, sequence, scaffold, input_hash)
        saved = reference[name]
        for key in ("n", "free_bases", "case3_ancestors_max"):
            require(row[key] == int(saved[key]), "RNAcentral feature mismatch: " + key)
        for key in ("sequence_sha256", "scaffold_sha256"):
            require(row[key] == saved[key], "RNAcentral S/G hash mismatch")
        require({int(k): v for k, v in json.loads(saved["case3_ancestor_hist"]).items()}
                == row["case3_ancestor_hist"], "RNAcentral ancestor histogram changed")
        # The historical TSV used d=100; explicitly evaluate its histogram at d=30.
        require(sum(v for k, v in row["case3_ancestor_hist"].items() if k > int(saved["case3_depth_cap"]))
                == int(saved["case3_truncated_endpoints"]), "historical truncation mismatch")
        reference_checks.append({key: row[key] for key in (
            "input", "case3_ancestors_max", "case3_depth_cap", "case3_truncated_endpoints")})
    require({row["input"] for row in reference_checks} == set(reference), "RNAcentral reference incomplete")

    summary = []
    for dataset_name in ("bpRNA_ge400", "NEAT1_1", "NEAT1_2"):
        group = [row for row in rows if row["dataset"] == dataset_name]
        maximum = max(row["case3_ancestors_max"] for row in group)
        summary.append({
            "dataset": dataset_name, "rna_count": len(group),
            "n_min": min(row["n"] for row in group), "n_max": max(row["n"] for row in group),
            "case3_ancestors_max": maximum, "case3_depth_cap": CAP,
            "rnas_over_30": sum(row["case3_ancestors_max"] > CAP for row in group),
            "case3_truncated_endpoints": sum(row["case3_truncated_endpoints"] for row in group),
            "free_bases": sum(row["free_bases"] for row in group),
            "maximum_branching_ancestors_of_free_base":
                max(row["maximum_branching_ancestors_of_free_base"] for row in group),
            "inputs_at_maximum": ";".join(row["input"] for row in group if row["case3_ancestors_max"] == maximum),
        })

    file_sources = [COHORT, dataset_path.relative_to(ROOT), NEAT1, RNA_FEATURES,
                    Path(__file__).relative_to(ROOT)] + [Path("scripts") / name for name in
                    ("scaffold_features.py", "rna_records.py", "task_storage.py")]
    validation = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "analysis_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "analysis_kind": "saved scaffold geometry only; no RNA CLI or scaffold generation",
        "definition": "scaffold_features.scaffold_features; proper branching ancestors of free-base owner; branching root included",
        "case3_depth_cap": CAP, "target_rnas_checked": len(rows), "independent_histograms_checked": len(rows) + len(reference_checks),
        "input_hashes_verified": True, "cohort_and_dataset_hashes_verified": True,
        "neat1_scaffold_settings": neat_manifest["scaffold_settings"],
        "frozen_sources": sources, "source_sha256": {str(path): digest(ROOT / path) for path in file_sources},
        "rnacentral_reference": reference_checks,
        "rnacentral_saved_depth_cap": sorted({int(row["case3_depth_cap"]) for row in saved_rows}),
        "rnacentral_maximum": max(row["case3_ancestors_max"] for row in reference_checks),
        "rnacentral_truncated_endpoints_at_30": sum(row["case3_truncated_endpoints"] for row in reference_checks),
    }
    output.mkdir(parents=True, exist_ok=True)
    for name in ("per_rna.tsv", "summary.tsv", "validation.json"):
        require(not (output / name).exists(), "refusing to replace existing output: " + name)
    write_tsv(output / "per_rna.tsv", rows)
    write_tsv(output / "summary.tsv", summary)
    validation["output_sha256"] = {name: digest(output / name) for name in ("per_rna.tsv", "summary.tsv")}
    (output / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    analyze(parser.parse_args().out)
