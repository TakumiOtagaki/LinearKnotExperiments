#!/usr/bin/env python3
"""Remeasure 337 frozen inputs and 548 frozen targets against ViennaRNA."""
import csv
import json
import os
from pathlib import Path
import subprocess
import tarfile
import tempfile

import RNA
import check_energy_vienna_20261017 as old
from prepare_paper_rerun_89baea52 import BUILD, DOC, ROOT, SOURCE, release_ready
from run_cli import digest, now, write_json


def main():
    release_ready()
    for key in list(os.environ):
        if key.startswith("CPARTY_"):
            del os.environ[key]
    os.environ.update(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", LC_ALL="C")
    output = ROOT / "results/local/energy-vienna-check-89baea52-20261018"
    output.mkdir(parents=True, exist_ok=False)
    original = ROOT / "results/local/energy-vienna-check-20261017"
    inputs = [json.loads(line) for line in (original / "inputs.jsonl").read_text().splitlines()]
    with (original / "structure_metrics.tsv").open() as stream:
        fixed = list(csv.DictReader(stream, delimiter="\t"))
    assert len(inputs) == 337 and len(fixed) == 548
    records = {r["id"]: r for r in inputs}
    param = SOURCE / "params/rna_Turner04.par"
    RNA.params_load(str(param))
    RNA.init_rand(20261017)
    deps = [Path(__file__), Path(old.__file__), ROOT / "scripts/energy_vienna_public_89baea52.cc",
            ROOT / "scripts/energy_vienna_probe_20261017.cc", original / "inputs.jsonl",
            original / "structure_metrics.tsv", param, BUILD / "CMakeCache.txt",
            BUILD / "libcparty_scfg2_exact.a", BUILD / "libRNA.a"]
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in deps}
    provenance = {"started_utc": now(), "source": str(SOURCE), "build": str(BUILD),
        "vienna_version": RNA.__version__, "dependencies_sha256_before": hashes,
        "scope": "dangles=2 public ExactSession; dangles=0 historical internal control; PF-smoothed fixed-target energy, all i<j BPP coordinates",
        "historical_numbers_reused": False, "frozen_inputs": 337, "frozen_structure_conditions": 548}
    write_json(output / "provenance.json", provenance)
    failures, ensemble_rows, structure_rows, mfe_rows = [], [], [], []
    with tempfile.TemporaryDirectory(prefix="vienna-89baea52-", dir=os.environ["TMPDIR"]) as td, \
            tarfile.open(output / "raw_outputs.tar.gz", "w:gz") as archive, \
            (output / "runs.jsonl").open("w") as log:
        scratch = Path(td)
        binary = scratch / "public_probe"
        command = ["c++", "-O2", "-std=c++17", "-I"+str(SOURCE / "src"),
                   "-I"+str(ROOT / "scripts"), str(ROOT / "scripts/energy_vienna_public_89baea52.cc"),
                   str(BUILD / "libcparty_scfg2_exact.a"), str(BUILD / "libRNA.a"), "-o", str(binary)]
        compiled = subprocess.run(command, capture_output=True, text=True, timeout=180)
        old.archive_bytes(archive, "build.json", json.dumps({"argv": command, "returncode": compiled.returncode,
            "stdout": compiled.stdout, "stderr": compiled.stderr}))
        compiled.check_returncode()
        for path in deps[:6]:
            archive.add(path, arcname="sources/"+path.name)

        def invoke(record, mode, dangles, targets, name):
            r = old.invoke(binary, param, record, mode, dangles, targets, archive, name)
            log.write(json.dumps({k: v for k, v in r.items() if k not in ("bpp", "rules", "targets")})+"\n")
            log.flush()
            if r["status"] != "success":
                failures.append({"condition": name, "status": r["status"]})
                return None
            return r

        for index, r in enumerate(inputs):
            fc, reference = old.vienna_ensemble(r["sequence"])
            bpp = fc.bpp()
            observed = invoke(r, "ensemble", 2, [], r["id"]+"/public")
            if observed is None:
                continue
            all_pairs = [(i, j) for i in range(1, r["length"]+1) for j in range(i+1, r["length"]+1)]
            row = {"id": r["id"], "length": r["length"], "coordinates": len(all_pairs),
                "abs_logz_error": abs(observed["logz"]-reference["vienna_logz"]),
                "max_abs_bpp_error": max(abs(observed["bpp"].get(p, 0)-bpp[p[0]][p[1]]) for p in all_pairs)}
            ensemble_rows.append(row)
            old.archive_bytes(archive, r["id"]+"/vienna_reference.json", json.dumps(reference))
            old.archive_bytes(archive, r["id"]+"/vienna_bpp.tsv", "".join(f"{i}\t{j}\t{bpp[i][j]:.17g}\n" for i,j in all_pairs))
            independent = old.structure_reference(r["sequence"], observed["mfe_structure"], 2)
            candidate = old.structure_reference(r["sequence"], reference["vienna_mfe_structure"], 2)
            mfe_rows.append({"id": r["id"], "public_energy": observed["mfe_energy"],
                "independent_pf_energy": independent["vienna_structure_pf_energy"],
                "abs_rescore_error": abs(observed["mfe_energy"]-independent["vienna_structure_pf_energy"]),
                "common_pf_gap_to_vienna_candidate": independent["vienna_structure_pf_energy"]-candidate["vienna_structure_pf_energy"]})
            if (index+1) % 25 == 0:
                print(f"Vienna public ensemble {index+1}/{len(inputs)}", flush=True)
        grouped = {}
        for index, row in enumerate(fixed):
            grouped.setdefault((row["id"], row["scaffold"], int(row["dangles"])), []).append((index, row))
        for gi, ((rid, scaffold, dangles), rows) in enumerate(grouped.items()):
            record = dict(records[rid], scaffold=scaffold)
            targets = [old.db_from_pairs(record["length"], set(old.pairs(row["structure"]))-set(old.pairs(scaffold))) for _, row in rows]
            observed = invoke(record, "eval", dangles, targets, f"fixed/{gi:04d}")
            if observed is None:
                continue
            if len(observed["targets"]) != len(rows):
                raise ValueError("fixed-target observation count differs")
            for (index, row), evaluated in zip(rows, observed["targets"]):
                reference = old.structure_reference(record["sequence"], row["structure"], dangles)
                ok = evaluated["same_pairs"] and evaluated["reachable"]
                structure_rows.append({"condition": index, "id": rid, "dangles": dangles,
                    "reachable": ok, "energy": evaluated["energy"],
                    "abs_pf_energy_error": abs(evaluated["energy"]-reference["vienna_structure_pf_energy"]) if ok else None})
                old.archive_bytes(archive, f"fixed/{index:04d}_reference.json", json.dumps(reference))
                if not ok:
                    failures.append({"condition": index, "id": rid, "status": "unreachable"})
    old.write_tsv(output / "ensemble_metrics.tsv", ensemble_rows)
    old.write_tsv(output / "structure_metrics.tsv", structure_rows)
    old.write_tsv(output / "mfe_comparison.tsv", mfe_rows)
    summary = {"frozen_inputs": len(inputs), "completed_ensembles": len(ensemble_rows),
        "bpp_coordinates": sum(r["coordinates"] for r in ensemble_rows),
        "frozen_structure_conditions": len(fixed), "completed_structure_conditions": len(structure_rows),
        "public_dangles2_conditions": sum(r["dangles"] == 2 for r in structure_rows),
        "internal_dangles0_controls": sum(r["dangles"] == 0 for r in structure_rows),
        "max_abs_logz_error": max((r["abs_logz_error"] for r in ensemble_rows), default=None),
        "max_abs_bpp_error": max((r["max_abs_bpp_error"] for r in ensemble_rows), default=None),
        "max_abs_fixed_pf_energy_error": max((r["abs_pf_energy_error"] for r in structure_rows if r["reachable"]), default=None),
        "max_abs_mfe_rescore_error": max((r["abs_rescore_error"] for r in mfe_rows), default=None),
        "max_common_pf_mfe_gap": max((r["common_pf_gap_to_vienna_candidate"] for r in mfe_rows), default=None),
        "failures": failures, "tolerance": 1e-9}
    summary["passed"] = (not failures and len(ensemble_rows) == 337 and len(structure_rows) == 548 and
        all(summary[k] is not None and summary[k] <= 1e-9 for k in
            ("max_abs_logz_error", "max_abs_bpp_error", "max_abs_fixed_pf_energy_error", "max_abs_mfe_rescore_error", "max_common_pf_mfe_gap")))
    write_json(output / "summary.json", summary)
    provenance.update(finished_utc=now(), dependencies_unchanged=hashes == {str(p.relative_to(ROOT)): digest(p) for p in deps},
                      outputs_sha256={p.name: digest(p) for p in output.iterdir() if p.name != "provenance.json"})
    write_json(output / "provenance.json", provenance)
    if not summary["passed"] or not provenance["dependencies_unchanged"]:
        raise RuntimeError("Vienna comparison failed; retained failures are part of the denominator")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
