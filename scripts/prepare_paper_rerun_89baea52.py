#!/usr/bin/env python3
"""Freeze the complete corrected-release paper matrix. Never submit jobs."""
import argparse
import copy
import io
import json
import os
from pathlib import Path
import shlex
import subprocess
import tarfile

from file_budget import check_budget
from make_rerun_configs import write_configs
from plan_array import load_plan, prepare, wall_seconds
from run_cli import ROOT, digest, now, verify_lock, write_json

REV = "89baea52"
COMMIT = "89baea52c9435b18829ef0867a66a24be56d70f9"
DATE = "20261018"
SOURCE = ROOT / ("deps/LinearCParty-" + REV)
BUILD = ROOT / ("build/linear-cparty-" + REV)
LOCK = "manifests/software-pilot-89baea52.json"
DOC = ROOT / "docs/rerun_89baea52_20261018"
CONFIG = ROOT / "configs/rerun_89baea52_20261018"
DATA = ROOT / "data/rerun-89baea52-20261018"


def rel(path):
    return str(Path(path).relative_to(ROOT))


def release_ready():
    receipt = json.loads((ROOT / "manifests/release-preparation-89baea52.json").read_text())
    stages = receipt["stages"]
    if (receipt.get("dry_run") or receipt["release"]["commit"] != COMMIT or
            stages.get("required_suite", {}).get("status") not in ("passed", "verified_existing_lock") or
            stages.get("lock", {}).get("status") not in ("created", "verified_existing") or
            stages.get("comparison", {}).get("successful") is not True):
        raise ValueError("Wait for prepare_release.py to finish its suite, lock and comparison")
    lock = json.loads((ROOT / LOCK).read_text())
    if lock["source_commit"] != COMMIT:
        raise ValueError("wrong release lock")
    verify_lock(ROOT / LOCK, SOURCE, BUILD, [{"name": n} for n in lock["binaries"]],
                SOURCE / "params/rna_Turner04.par")


def base_config(name, inputs, modes, programs, cap=8192, timeout=1800):
    return {"schema_version": 1, "name": name, "purpose": "formal", "software_lock": LOCK,
            "energy_parameter_file": "params/rna_Turner04.par",
            "model_settings": {"temperature_celsius": 37, "build_revision_verified": True},
            "environment": {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
                            "MKL_NUM_THREADS": "1", "LC_ALL": "C"},
            "timeout_seconds": timeout, "memory_limit_mib": cap, "warmups": 0, "repetitions": 1,
            "order": {"kind": "blocked_random", "seed": 20261017},
            "inputs": inputs, "modes": modes, "programs": programs}


def small_array(name, hi, seconds, per_task, workers=1, cap=8192, timeout=1800):
    return {"schema_version": 1, "measurement": "accuracy", "seed": 20261017,
            "compact_plan": True, "max_concurrent_tasks": 30,
            "file_budget": {"project": 102571, "filesystem": "/gs/bs", "hard_file_limit": 24000000,
                            "max_remaining_fraction": 0.1},
            "tiers": [{"name": name, "min_length": 1, "max_length": hi,
                       "resource": "cpu_16", "h_rt": "04:00:00", "target_seconds": 10000,
                       "max_inputs_per_task": per_task, "workers": workers,
                       "memory_budget_mib": 35045, "task_overhead_mib": 1024,
                       "job_limits": {"memory_limit_mib": cap, "timeout_seconds": timeout},
                       "cost_model": {"seconds_constant_per_condition": seconds,
                                      "seconds_per_nt_condition": 1e-9,
                                      "memory_base_mib": cap, "mib_per_nt": 0, "safety_factor": 1,
                                      "basis": "Conservative scheduling allowance; no new-release pilot or assumed parallel speedup."}}]}


def cli(name, args=None, bpp=False):
    p = {"name": name, "args": args or []}
    if bpp:
        p["outputs"] = [{"name": "bpp.tsv", "flag": "--bpp"}]
    return p


def beam(b):
    return {"name": f"beam{b}_outer", "args": ["-b", str(b), "--max-hairpin", str(max(b, 30)),
            "--max-case3-depth", str(max(b, 30)), "--ranking", "outer"]}


def write_inputs():
    DATA.mkdir(exist_ok=False)
    groups = {}
    origins = {}
    with tarfile.open(DATA / "inputs.tar", "w") as archive:
        def add(group, name, payload):
            member = group + "/" + name + ".txt"
            info = tarfile.TarInfo(member); info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))
            groups.setdefault(group, []).append({"name": name, "path": rel(DATA / member),
                "input_archive": rel(DATA / "inputs.tar"), "input_member": member})

        for group, path in (("pk200", ROOT / "results/local/pk-recovery-bprna-20261016/pk_sel.json"),
                            ("pk80", ROOT / "results/local/pk-recovery-long-rrna-20261016/mid300_700_selection.json")):
            origins[rel(path)] = digest(path)
            for r in json.loads(path.read_text()):
                g = ["."] * len(r["seq"])
                for i, j in r["nested"]:
                    g[i-1], g[j-1] = "(", ")"
                add(group, r["id"], (f'>{r["id"]}\n{r["seq"]}\n' + "".join(g) + "\n").encode())
        path = ROOT / "results/local/pk-recovery-long-rrna-20261016/outputs.tar.gz"
        origins[rel(path)] = digest(path)
        with tarfile.open(path) as old:
            for rid in ("bpRNA_CRW_20565", "bpRNA_CRW_54735", "bpRNA_CRW_55337"):
                add("rrna", rid, old.extractfile(rid + ".txt").read())
        for group, path in (("terc", ROOT / "results/local/figure4-terc-20261017/in.txt"),
                            ("figure1", ROOT / "results/figures/figure1_20261016/data/in.txt")):
            origins[rel(path)] = digest(path)
            payload = path.read_bytes()
            add(group, payload.decode().splitlines()[0][1:], payload)
    write_json(DATA / "provenance.json", {"created_utc": now(), "origins_sha256": origins,
               "inputs_tar_sha256": digest(DATA / "inputs.tar"), "groups": groups})
    return groups


def generate():
    if CONFIG.exists() or DATA.exists():
        raise ValueError("config/data output already exists; preserve immutable preparation")
    # Joint bound includes all plans, task execution, custom validation wrappers,
    # configs and documentation; it is checked before any matrix is written.
    budget = check_budget({"estimated_shared_entries": 20000,
                           "basis": "Whole paper upper allowance, including all array and local plans, config/data and validation files"}, {})
    DOC.mkdir(exist_ok=True)
    write_json(DOC / "matrix_file_budget.json", budget)
    CONFIG.mkdir(); (CONFIG / "arrays").mkdir()
    groups = write_inputs()
    rows = []

    def save(name, config, array, build=BUILD, local=False, template=None):
        config["name"] = name
        if local:
            config["purpose"] = "development"
            config["analysis"] = {"scope": "Small example predictions on the certified release; local execution is allowed; no formal performance measurement."}
        config["rerun_preparation"] = {"release": COMMIT, "template": str(template) if template else None}
        array.update(compact_plan=True, max_concurrent_tasks=30)
        cp, ap = CONFIG / (name + ".json"), CONFIG / "arrays" / (name + ".json")
        write_json(cp, config); write_json(ap, array)
        rows.append({"name": name, "config": rel(cp), "array_config": rel(ap),
                     "software_root": rel(SOURCE), "build_dir": rel(build), "execution": "local" if local else "array"})

    for stem in ("e0_bprna_legacy_d21618c5_20261014", "e1_bprna_short_d21618c5_20261014",
                 "e1_bprna_long_d21618c5_20261014"):
        cp, ap = ROOT / "configs" / (stem + ".json"), ROOT / "configs/arrays" / (stem + ".json")
        c, a = json.loads(cp.read_text()), json.loads(ap.read_text())
        c["software_lock"] = LOCK
        build = BUILD
        if stem.startswith("e0"):
            c["software_lock"] = "manifests/software-e0-89baea52-20261018.json"
            c["analysis"]["comparisons"]["compatible_vs_oracle"] += "; additionally patch 0009, inner multiloop closing neighbors"
            build = ROOT / "build/e0-observers-89baea52-800e95a-6patches"
        name = stem.replace("d21618c5_20261014", "89baea52_20261018").replace("_", "-")
        save(name, c, a, build, template=rel(cp))

    settings = json.loads((ROOT / "configs/rerun_d02d261c_cpu40_20261017.json").read_text())
    settings = json.loads(json.dumps(settings).replace("d02d261c", REV).replace("20261017", DATE))
    settings["note"] = "All E2 measurements rerun after the ordinary multiloop closing correction; same cpu_40/10 workers/8192 MiB as d02d261c."
    settings_path = CONFIG / "e2_settings.json"
    write_json(settings_path, settings)
    generated = write_configs(settings, CONFIG / "e2")
    for row in generated["variants"]:
        row["name"] = json.loads((ROOT / row["config"]).read_text())["name"]
        row["execution"] = "array"
        rows.append(row)

    for cp in sorted((ROOT / "configs/e3_neat1_d02d261c_20261017").glob("e3*.json")):
        ap = cp.with_name("array_" + cp.name)
        c, a = json.loads(cp.read_text()), json.loads(ap.read_text())
        c["software_lock"] = LOCK
        name = c["name"].replace("d02d261c-20261017", "89baea52-20261018")
        save(name, c, a, template=rel(cp))
    for kind in ("bpp", "capr", "raccess"):
        cp = ROOT / f"configs/e3_neat1_variants_{kind}_e2cecbb4_20261016.json"
        c = json.loads(cp.read_text()); a = json.loads((cp.parent / "arrays" / cp.name).read_text())
        c["software_lock"] = LOCK
        save(f"e3-neat1-variants-{kind}-89baea52-20261018", c, a, template=rel(cp))

    exact = {"name": "exact", "args": ["--exact"]}
    for param in ("turner04", "dp09"):
        name = f"pk-recovery-200-{param}-89baea52-20261018"
        c = base_config(name, groups["pk200"], [exact], [cli("LinearCParty", ["--bpp-cutoff", "0"], True)])
        if param == "dp09":
            c["energy_parameter_file"] = "params/rna_DirksPierce09.par"
            c["software_lock"] = "manifests/software-pilot-89baea52-dp09.json"
        save(name, c, small_array("pk200", 400, 600, 10))
    name = "pk-recovery-80-89baea52-20261018"
    save(name, base_config(name, groups["pk80"], [beam(400)], [cli("LinearCParty", ["--bpp-cutoff", "0"], True)]),
         small_array("pk80", 700, 600, 10))
    name = "pk-recovery-rrna-89baea52-20261018"
    save(name, base_config(name, groups["rrna"], [beam(400)],
                          [cli("LinearCParty", ["--bpp-cutoff", "0"], True), cli("LinearPKCapR")], cap=16384, timeout=7200),
         small_array("rrna", 5000, 3600, 1, cap=16384, timeout=7200))
    for group, seed, modes in (("terc", "20261017", [exact, beam(200)]),
                               ("figure1", "20261016", [exact])):
        name = f"{group}-89baea52-20261018"
        programs = [cli("LinearCParty", ["--bpp-cutoff", "0"], True), cli("LinearPKCapR"),
                    cli("LinearPKRaccess", ["--access-len", "1,4,8", "--probability"]),
                    cli("LinearPKSampling", ["-k", "1000", "--seed", seed, "--unique"]),
                    cli("LinearPKMEA", ["--gamma", "1"]), cli("LinearPKSubopt", ["-k", "5"])]
        c = base_config(name, groups[group], modes, programs)
        c["software_lock"] = "manifests/software-examples-89baea52-20261018.json"
        save(name, c,
             small_array(group, 400, 120, 1), local=True)
        if group == "terc":
            name = "terc-byloop-89baea52-20261018"
            save(name, base_config(name, groups[group], modes,
                 [cli("LinearPKRaccess", ["--access-len", "1,4,8", "--probability", "--by-loop"])]),
                 small_array(group, 400, 120, 1), local=True)

    name = "pk-membership-89baea52-20261018"
    c = base_config(name, groups["pk200"], [{"name": "exact", "args": []}], [cli("E0Membership")], timeout=900)
    c["software_lock"] = "manifests/software-membership-89baea52-20261018.json"
    save(name, c, small_array("membership", 400, 900, 10, timeout=900),
         ROOT / "build/membership-89baea52")
    write_json(CONFIG / "matrix_manifest.json", {"release": COMMIT, "created_utc": now(), "variants": rows})
    return rows


def plan_all(rows):
    example_lock = ROOT / "manifests/software-examples-89baea52-20261018.json"
    if not example_lock.exists():
        release_ready()
        lock = copy.deepcopy(json.loads((ROOT / LOCK).read_text()))
        lock["binaries"]["LinearPKSubopt"] = digest(BUILD / "LinearPKSubopt")
        lock.update(created_utc=now(), base_lock_sha256=digest(ROOT / LOCK),
            verification_method=lock["verification_method"]+"; Subopt from the same full-target dedicated build, exercised by the certified required suite's linear_pk_subopt_cli check")
        lock.setdefault("extra_files", []).append({"path": LOCK, "sha256": digest(ROOT / LOCK)})
        write_json(example_lock, lock)
    dp_lock = ROOT / "manifests/software-pilot-89baea52-dp09.json"
    if not dp_lock.exists():
        lock = copy.deepcopy(json.loads((ROOT / LOCK).read_text()))
        lock.update(created_utc=now(), energy_parameter_file="params/rna_DirksPierce09.par",
                    energy_parameter_sha256=digest(SOURCE / "params/rna_DirksPierce09.par"),
                    parameter_variant_from_lock_sha256=digest(ROOT / LOCK))
        write_json(dp_lock, lock)
    summary = []
    for row in rows:
        out = ROOT / "results/runs" / row["name"]
        jobs = ROOT / "jobs" / row["name"]
        if out.exists() or jobs.exists():
            plan = load_plan(out / "array_plan.json")
            if (plan["original_config_sha256"] != digest(ROOT / row["config"]) or
                    plan["array_config_sha256"] != digest(ROOT / row["array_config"]) or
                    plan["software_root"] != str(ROOT / row["software_root"]) or
                    plan["build_dir"] != str(ROOT / row["build_dir"]) or
                    not plan["dependencies_verified_at_planning"] or plan["dry_run"]):
                raise ValueError("existing frozen plan differs: "+row["name"])
            verify_lock(out / "software_lock.json", ROOT / row["software_root"], ROOT / row["build_dir"],
                        plan["config"]["programs"], SOURCE / plan["config"]["energy_parameter_file"])
        else:
            plan = prepare(ROOT / row["config"], ROOT / row["array_config"], out, jobs,
                           ROOT / row["software_root"], ROOT / row["build_dir"])
        for script in plan["job_scripts"]:
            tasks = [t for t in plan["tasks"] if t["tier"] == script["tier"]]
            summary.append({"name": row["name"], "execution": row["execution"],
                "path": rel(script["path"]), "tier": script["tier"], "tasks": len(tasks),
                "resource": tasks[0]["requested_resource"], "workers": tasks[0]["workers"],
                "h_rt": tasks[0]["h_rt"], "sequential_hours": sum(t["estimated_seconds"] for t in tasks)/3600,
                "max_task_hours": max(t["estimated_seconds"] for t in tasks)/3600,
                "requested_task_hours": len(tasks)*wall_seconds(tasks[0]["h_rt"])/3600,
                "memory_limit_mib": tasks[0]["job_limits"]["memory_limit_mib"],
                "plan": rel(out / "array_plan.json"), "script_sha256": script["sha256"]})
    directory = ROOT / "jobs/e0-validation-89baea52-20261018"
    logs = ROOT / "results/runs/e0-validation-89baea52-20261018/scheduler_logs"
    directory.mkdir(exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)
    for kind, hours in (("enumeration", .5), ("decoding", .5), ("byloop", .25), ("vienna", .75), ("linearspace", .1)):
        path = directory / (kind+".sh")
        content = "\n".join(["#!/usr/bin/env bash", "#$ -cwd", "#$ -S /bin/bash", "#$ -j y",
            "#$ -N lcp_"+kind, "#$ -l cpu_16=1", "#$ -l h_rt=02:00:00", "#$ -t 1-1", "#$ -tc 30",
            "#$ -o "+str(logs), "", "set -euo pipefail", "cd "+shlex.quote(str(ROOT)),
            "export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 LC_ALL=C",
            ': "${TMPDIR:?node-local TMPDIR is required}"', "ulimit -v 16777216",
            'exec "${LCP_PYTHON:-python3}" -B scripts/run_paper_validation_89baea52.py '+kind, ""])
        if path.exists():
            if path.read_text() != content:
                raise ValueError("existing validation wrapper differs: "+str(path))
        else:
            path.write_text(content)
            path.chmod(0o755)
        summary.append(dict(name="e0-validation-89baea52-20261018", execution="array", path=rel(path),
            tier=kind, tasks=1, resource="cpu_16=1", workers=1, h_rt="02:00:00", sequential_hours=hours,
            max_task_hours=hours, requested_task_hours=2, memory_limit_mib=16384,
            plan="scripts/run_paper_validation_89baea52.py", script_sha256=digest(path)))
    write_json(DOC / "job_scripts.json", summary)
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--plan-existing", action="store_true", help="plan previously generated configs after all observer locks exist")
    args = p.parse_args()
    if args.plan_existing:
        release_ready()
        rows = json.loads((CONFIG / "matrix_manifest.json").read_text())["variants"]
        print(json.dumps(plan_all(rows), indent=2))
    else:
        print(json.dumps(generate(), indent=2))


if __name__ == "__main__":
    main()
