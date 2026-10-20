#!/usr/bin/env python3
"""Prepare centroid replacements on the adopted build; never execute/submit them."""
import copy
import json
from pathlib import Path
import shlex

from file_budget import check_budget
from plan_array import load_plan, prepare
from prepare_paper_rerun_89baea52 import BUILD, COMMIT, ROOT, SOURCE, release_ready
from run_cli import digest, now, verify_lock, write_json

TAG = "89baea52-20261020"
CONFIG = ROOT / "configs/centroid_89baea52_20261020"
DOC = ROOT / "docs/centroid_89baea52_20261020"
OLD = ROOT / "configs/rerun_89baea52_20261018"
TOOL = ROOT / "tools/validation_centroid_89baea52_20261020/decoding"


def rel(p):
    return str(p.relative_to(ROOT))


def save(path, value):
    if path.exists():
        if json.loads(path.read_text()) != value:
            raise ValueError("existing configuration differs: " + str(path))
    else:
        write_json(path, value)


def main():
    release_ready()
    budget = check_budget({"estimated_shared_entries": 500,
        "basis": "Three matrix plans (four tasks), E0 wrapper/config, source copies, reports and execution allowance"}, {})
    CONFIG.mkdir(exist_ok=True)
    (CONFIG / "arrays").mkdir(exist_ok=True)
    DOC.mkdir(exist_ok=True)
    write_json(DOC / "file_budget.json", budget)
    variants = []
    observed_files = []
    old_run = ROOT / "results/runs/e3-neat1-mea-d02d261c-20261017"
    for stem in ("e3-neat1", "terc", "figure1"):
        original = OLD / (("e3-neat1-mea-89baea52-20261018" if stem == "e3-neat1" else stem + "-89baea52-20261018") + ".json")
        config = json.loads(original.read_text())
        array = json.loads((OLD / "arrays" / original.name).read_text())
        name = stem + "-centroid-" + TAG
        config["name"] = name
        for program in config["programs"]:
            if program["name"] == "LinearPKMEA":
                program["name"] = "LinearPKCentroid"
                assert program["args"][:2] == ["--gamma", "1"]
        config["analysis"] = {"scope": "Replace MEA with gamma-centroid, gamma=1; preserve all other scientific settings.",
                              "gamma_meaning": "gamma=1 is the ordinary centroid; gamma is not the MEA objective's pair/unpaired tradeoff."}
        config["rerun_preparation"] = {"release": COMMIT, "template": rel(original),
                                       "template_sha256": digest(original), "prepared_for": "centroid replacement"}
        if stem == "e3-neat1":
            array["note"] = "Same cpu_40/3 workers/20 GiB per process as prior MEA; h_rt about twice the observed slowest width. No centroid timing measured."
            for tier in array["tiers"]:
                path = old_run / "tasks" / tier["name"] / "000001/run/task_results.jsonl"
                observed_files.append(path)
                records = [json.loads(line) for line in path.read_text().splitlines()]
                assert len(records) == 3 and all(r["status"] == "success" for r in records)
                assert len({r["input"] for r in records}) == 1
                input_name = records[0]["input"]
                input_path = ROOT / next(i["path"] for i in config["inputs"] if i["name"] == input_name)
                seconds = {r["mode"]: r["elapsed_seconds"] for r in records}
                short = tier["name"] == "neat1_short"
                tier.update(h_rt="00:22:00" if short else "02:50:00",
                            max_wall_fraction=.8, target_seconds=1056 if short else 8160)
                tier["job_limits"]["timeout_seconds"] = 1260 if short else 10080
                tier["cost_model"] = {"basis": "Same-input d02d261c MEA measured elapsed, 5% planning allowance; summed sequential work, no parallel speedup assumed. Centroid uncalibrated.",
                    "seconds_per_nt_condition": 1e-9, "mib_per_nt": 0, "safety_factor": 1.05,
                    "observed_inputs": {input_name: {"sha256": digest(input_path), "seconds_by_mode": seconds,
                        "max_seconds_by_mode": seconds, "memory_mib": 20480,
                        "source": rel(path), "source_sha256": digest(path),
                        "measured_peak_rss_kib": {r["mode"]: r["peak_rss_kib"] for r in records}}}}
        cp, ap = CONFIG / (name + ".json"), CONFIG / "arrays" / (name + ".json")
        save(cp, config); save(ap, array)
        out, jobs = ROOT / "results/runs" / name, ROOT / "jobs" / name
        if out.exists():
            plan = load_plan(out / "array_plan.json")
            assert plan["original_config_sha256"] == digest(cp) and plan["array_config_sha256"] == digest(ap)
            assert plan["dependencies_verified_at_planning"] and not plan["dry_run"]
            verify_lock(out / "software_lock.json", SOURCE, BUILD, config["programs"], SOURCE / config["energy_parameter_file"])
        else:
            plan = prepare(cp, ap, out, jobs, SOURCE, BUILD)
        variants.append({"name": name, "config": rel(cp), "array_config": rel(ap),
            "plan": rel(out / "array_plan.json"), "execution": "array" if stem == "e3-neat1" else "local",
            "tasks": [{k: t[k] for k in ("tier", "task_id", "estimated_seconds", "memory_cap_envelope_mib", "h_rt", "workers", "requested_resource")} for t in plan["tasks"]],
            "planned_jobs": sum(len(config["programs"])*len(config["modes"]) for _ in config["inputs"]),
            "scripts": [rel(Path(s["path"])) for s in plan["job_scripts"]]})
        if stem != "e3-neat1":
            local = jobs / "local.sh"
            content = "\n".join(["#!/usr/bin/env bash", "set -euo pipefail",
                "cd " + shlex.quote(str(ROOT)), ': "${TMPDIR:?node-local TMPDIR is required}"',
                "export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 LC_ALL=C",
                'exec "${LCP_PYTHON:-python3}" -B scripts/run_array_task.py --plan ' + shlex.quote(str(out / "array_plan.json"))
                + " --tier " + stem + ' --task-id 1 "$@"', ""])
            if local.exists() and local.read_text() != content:
                raise ValueError("existing local wrapper differs")
            local.write_text(content); local.chmod(0o755)
            variants[-1]["local_script"] = rel(local)
    cases_path = ROOT / "tools/validation_89baea52/baseline_15aa7742/inputs/cases.json"
    cases = json.loads(cases_path.read_text())
    ids = ["closure14", "random000", "random001", "random004", "random031", "random263"]
    config = dict(version=COMMIT, cases_file=rel(cases_path),
        sample_cases=[next(c for c in cases if c["id"] == cid) for cid in ids],
        sample_modes=["exact", "beam1", "beam3"], draws_per_condition=200000,
        sample_base_seed=202610160000, seed_rule="base + 1000*mode_index + case_index",
        gamma=[.5, 1, 2], main_gamma=1, mfe_modes=["exact", "beam1", "beam3", "beam20", "beam100"],
        mea_modes=["exact", "beam20", "beam100"], centroid_modes=["exact", "beam20", "beam100"],
        ranking="outer", outer_pilot_exponent="2/3", max_hairpin=30, max_case3_depth=30,
        temperature_C=37, dangles=2, tolerance=1e-9, family_alpha=.01,
        chi_square="Pearson; theory-only pooling expected<5; Holm across 18 conditions",
        frequency_bound="Hoeffding union bound over all theoretical structures in all conditions",
        centroid_objective="E[gamma*TP_pair + TN_pair] over all i<j; retain fixed G; gamma=1 is centroid")
    save(CONFIG / ("e0-decoding-centroid-" + TAG + ".json"), config)
    job = ROOT / "jobs" / ("e0-decoding-centroid-" + TAG) / "decoding.sh"
    logs = ROOT / "results/runs" / ("e0-decoding-centroid-" + TAG) / "scheduler_logs"
    job.parent.mkdir(exist_ok=True); logs.mkdir(parents=True, exist_ok=True)
    content = "\n".join(["#!/usr/bin/env bash", "#$ -cwd", "#$ -S /bin/bash", "#$ -j y", "#$ -N lcp_centroid_e0",
        "#$ -l cpu_16=1", "#$ -l h_rt=02:00:00", "#$ -t 1-1", "#$ -tc 30", "#$ -o " + str(logs),
        "", "set -euo pipefail", "cd " + shlex.quote(str(ROOT)),
        "export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 LC_ALL=C",
        ': "${TMPDIR:?node-local TMPDIR is required}"', "ulimit -v 16777216",
        'exec "${LCP_PYTHON:-python3}" -B scripts/run_centroid_validation_89baea52.py', ""])
    if job.exists() and job.read_text() != content:
        raise ValueError("existing E0 wrapper differs")
    job.write_text(content); job.chmod(0o755)
    paths = list(CONFIG.rglob("*.json")) + list(TOOL.glob("*")) + observed_files + [cases_path,
        ROOT / "tools/validation_89baea52/enumeration/enumerate.cc",
        ROOT / "scripts/prepare_centroid_89baea52.py", ROOT / "scripts/run_centroid_validation_89baea52.py", job]
    for row in variants:
        paths += [ROOT / row["plan"], *[ROOT / n for n in row["scripts"]]]
        if "local_script" in row:
            paths.append(ROOT / row["local_script"])
    paths += [ROOT / "manifests/software-pilot-89baea52.json", ROOT / "manifests/software-examples-89baea52-20261018.json"]
    manifest = {"created_utc": now(), "release": COMMIT, "software_root": rel(SOURCE), "build_dir": rel(BUILD),
        "variants": variants, "e0": {"script": rel(job), "cases": 729, "centroid_objective_checks": 6561,
            "centroid_cli_checks": 6561, "gamma": [.5, 1, 2], "estimated_minutes": 40, "h_rt": "02:00:00"},
        "files_sha256": {rel(p): digest(p) for p in paths if p.is_file()}, "submitted": False}
    write_json(DOC / "preparation.json", manifest)
    print(json.dumps({"variants": len(variants), "tasks": sum(len(r["tasks"]) for r in variants) + 1,
                      "planned_matrix_jobs": sum(r["planned_jobs"] for r in variants), "submitted": False}))


if __name__ == "__main__":
    main()
