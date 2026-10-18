#!/usr/bin/env python3
"""Collect the fresh paper values and completion evidence, excluding E1."""
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import tarfile

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT / "results/local"
SUFFIX = "89baea52-20261018"


def read(path):
    return json.loads(path.read_text())


def table(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main():
    runs = []
    for path in sorted((ROOT / "results/runs").glob("*" + SUFFIX)):
        if path.name.startswith("e1-") or not (path / "array_plan.json").exists():
            continue
        manifest = read(path / "run.json")
        assert manifest["status"] == "success", path
        assert manifest["software_before"]["commit"].startswith("89baea52")
        counts = dict(Counter(job["status"] for job in manifest["jobs"]))
        assert set(counts) == {"success"}, (path, counts)
        task_count = 0
        for metadata_path in path.glob("tasks/*/*/run/run.json"):
            metadata = read(metadata_path)
            assert metadata["status"] == "success", metadata_path
            assert metadata["software_before"]["commit"].startswith("89baea52")
            assert metadata["software_before"] == metadata["software_after"]
            assert metadata["binaries_before"] == metadata["binaries_after"]
            assert metadata["software_lock_valid_after_run"] is True
            assert not metadata["software_changed_during_run"]
            assert not metadata["binaries_changed_during_run"]
            task_count += 1
        runs.append(dict(run=path.name, statuses=counts,
            tasks_with_unchanged_dependencies=task_count,
            measured=sum(not job["warmup"] for job in manifest["jobs"]),
            manifest_sha256=hashlib.sha256((path / "run.json").read_bytes()).hexdigest()))
    assert len(runs) == 23, len(runs)
    e0 = read(LOCAL / f"e0-bprna-legacy-{SUFFIX}.json")
    parity = [r["metrics"] for r in e0["observations"] if r["comparison_kind"] == "numerical_parity"]
    sources = dict(
        enumeration=LOCAL / f"e0-enumeration-{SUFFIX}/summary.json",
        decoding=LOCAL / f"e0-decoding-{SUFFIX}/summary.json",
        byloop=LOCAL / f"byloop-validation-{SUFFIX}/summary.json",
        vienna=LOCAL / "energy-vienna-check-89baea52-repaired-20261018/summary.json",
        linearspace=LOCAL / f"validation-controller-linearspace-{SUFFIX}/provenance.json",
        e2=LOCAL / f"e2-{SUFFIX}/validation.json",
        neat1=LOCAL / f"e3-neat1-{SUFFIX}/validation.json",
        recovery=LOCAL / f"small-rerun-{SUFFIX}/recovery_summary.json",
        membership=LOCAL / f"pk-membership-{SUFFIX}/summary.json",
        terc=LOCAL / f"small-rerun-{SUFFIX}/terc_summary.json",
        figure1=LOCAL / f"small-rerun-{SUFFIX}/figure1_summary.json")
    values = {key: read(path) for key, path in sources.items()}
    assert values["vienna"]["passed"] and values["neat1"]["status"] == "passed"
    assert values["enumeration"]["integrity_ok"] and not values["decoding"]["scientific_failures"]
    assert values["linearspace"]["status"] == "success"
    assert all(r["returncode"] == 0 for r in values["linearspace"]["results"])
    with tarfile.open(sources["linearspace"].parent / "task_outputs.tar") as stream:
        linearspace_result = stream.extractfile("0.stdout").read().decode().splitlines()[-1]
    assert linearspace_result == "scfg2_inside_pair_max=ok comparisons=1120"
    resources = table(LOCAL / f"e2-{SUFFIX}/metrics.tsv")
    plot = table(ROOT / "results/figures/main_89baea52_20261018/plot_data.tsv")
    rrna = read(LOCAL / f"small-rerun-{SUFFIX}/rrna_artifact_index.json")
    rrna_probs = [p for row in rrna for p in row.get("reference_probabilities", [])]
    assert len(rrna_probs) == 24
    report = dict(revision="89baea52", runs=runs,
        cli_observations=sum(sum(r["statuses"].values()) for r in runs),
        e0=dict(max_abs_logz=max(abs(r["delta_log_z"]) for r in parity),
                max_abs_free_energy=max(abs(r["delta_free_energy"]) for r in parity),
                parity=e0["numerical_parity_counts"], closure=e0["closure_counts"]),
        validation_and_examples=values,
        linearspace_result=linearspace_result,
        rnacentral=[r for r in resources if r["run"] == "rnacentral_current"],
        resource_ratios=[r for r in plot if r["panel"] == "3C_ratio"],
        neat1_resources=table(LOCAL / f"e3-neat1-{SUFFIX}/resources.tsv"),
        neat1_region=table(LOCAL / f"e3-neat1-{SUFFIX}/byloop_regions.tsv"),
        rrna=dict(reference_pairs=24, reported_zeros=sum(p == 0 for p in rrna_probs),
                  max_probability=max(rrna_probs)),
        source_hashes={str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
                       for path in sources.values()})
    out = LOCAL / f"small-rerun-{SUFFIX}/paper_values.json"
    out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(dict(output=str(out), runs=len(runs), observations=report["cli_observations"], rrna=report["rrna"])))


if __name__ == "__main__":
    main()
