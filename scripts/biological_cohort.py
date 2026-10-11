"""Resolve a saved sequence selection and require verified generated scaffolds."""

import copy
import hashlib
import json
from pathlib import Path

from prepare_biological_inputs import read_fasta
from rna_records import read_record
from run_cli import digest, experiment_path

SETTINGS = {"beam": 200, "dangles": 2, "pair_probability_threshold": 0.7, "min_helix": 2}
LP_COMMIT = "b450fb3e63189073b68d385589035f992080aa3a"


def acquired_records(manifest_path):
    info = json.loads(manifest_path.read_text())
    source = experiment_path(info.get("copied_path", info.get("downloaded_path")))
    if digest(source) != info["sha256"]:
        raise ValueError("source checksum differs from acquisition manifest")
    records = read_fasta(source)
    if any(not s or set(s) - set("ACGU") for _, s in records):
        raise ValueError("the full clean cohort must contain only ACGU")
    if len({s for _, s in records}) != len(records) or len(records) != info["records"]:
        raise ValueError("clean cohort count or uniqueness changed")
    return info, source, records


def expand_cohort(config, dry_run=False):
    if "input_cohort" not in config:
        return config, None, False
    config = copy.deepcopy(config)
    cohort_path = experiment_path(config["input_cohort"])
    cohort = json.loads(cohort_path.read_text())
    acquisition = experiment_path(cohort["source_manifest"])
    if digest(acquisition) != cohort["source_manifest_sha256"]:
        raise ValueError("cohort acquisition manifest changed")
    info, source, records = acquired_records(acquisition)
    sequences = dict(records)
    dataset_path = experiment_path(cohort["scaffold_dataset"])
    dataset = json.loads(dataset_path.read_text()) if dataset_path.exists() else None
    if dataset is not None:
        if (dataset.get("status") != "success" or dataset.get("software_unchanged") is not True or
                dataset["source_sha256"] != info["sha256"] or
                dataset["acquisition_manifest_sha256"] != digest(acquisition) or
                dataset["scaffold_settings"] != SETTINGS or
                dataset["linearpartition_source"]["commit"] != LP_COMMIT):
            raise ValueError("scaffold dataset incomplete or provenance differs")
        generated = {r["id"]: r for r in dataset["records"]}
        if len(generated) != len(records) or set(generated) != set(sequences):
            raise ValueError("full scaffold cohort differs; do not discard failures")
    elif not dry_run:
        raise ValueError("pending scaffolds: finish and aggregate G before verified planning")
    config.setdefault("inputs", [])
    for entry in cohort["records"]:
        name = entry["id"]
        sequence = sequences[name]
        sha = hashlib.sha256(sequence.encode()).hexdigest()
        if sha != entry["sequence_sha256"] or len(sequence) != entry["length"]:
            raise ValueError("selected sequence changed: " + name)
        path = dataset_path.parent / (name + ".txt")
        item = {"name": name, "path": str(path), "length": len(sequence)}
        if dataset is None:
            item.update(scaffold_pending=True, planning_sequence=sequence)
        else:
            saved = generated[name]
            record = read_record(path)
            if (saved["scaffold_status"] != "success" or digest(path) != saved["input_sha256"] or
                    record[:2] != (name, sequence) or saved["sequence_sha256"] != sha or
                    hashlib.sha256(record[2].encode()).hexdigest() != saved["scaffold_sha256"]):
                raise ValueError("generated scaffold input changed: " + name)
        config["inputs"].append(item)
    provenance = {"cohort_sha256": digest(cohort_path), "source_sha256": digest(source),
                  "acquisition_manifest_sha256": digest(acquisition),
                  "scaffold_dataset_sha256": digest(dataset_path), "pending_scaffolds": dataset is None}
    return config, provenance, dataset is None
