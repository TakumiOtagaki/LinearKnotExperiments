#!/usr/bin/env python3
"""Account for saved fixed-structure differences using ViennaRNA local helpers.

No SCFG2 code is imported here. ViennaRNA supplies every numerical local factor.
"""
import argparse
import csv
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import tarfile
import tempfile

import RNA
import RNA._RNA as native_rna_extension

sys.dont_write_bytecode = True
from check_energy_vienna_20261017 import ROOT, digest, loop_tree, make_md, write_tsv


def native_terms(sequence, structure, dangles):
    n = len(sequence)
    encoded = [0]+["ACGU".index(x)+1 for x in sequence]
    md = make_md(dangles)
    p = RNA.exp_param(md)
    rt = p.kT/1000.
    energy = lambda q: -rt*math.log(q)
    pp, children = loop_tree(structure)
    ends = dict(pp)
    terms = []

    def add(kind, i, j, value, **extras):
        terms.append(dict(term=kind, i=i, j=j, native_pf_energy=value,
                          guarded_pf_energy=value, **extras))

    def stem(i, j, multiloop=False, closing=False):
        t = md.pair[encoded[j]][encoded[i]] if closing else md.pair[encoded[i]][encoded[j]]
        au = energy(p.expTermAU) if t > 2 else 0.
        if closing:
            si, sj = encoded[j-1], encoded[i+1]
            gsi = si if j < n else -1
            gsj = sj if i > 1 else -1
        else:
            si = encoded[i-1] if i > 1 else -1
            sj = encoded[j+1] if j < n else -1
            gsi, gsj = si, sj
        if dangles == 0:
            si = sj = gsi = gsj = -1
        helper = RNA.exp_E_MLstem if multiloop else RNA.exp_E_ExtLoop
        baseline = energy(helper(t, -1, -1, p))
        full = energy(helper(t, si, sj, p))
        guarded = energy(helper(t, gsi, gsj, p))
        prefix = "closing" if closing else "multiloop_branch" if multiloop else "exterior_branch"
        add(prefix+"_terminal_AU_GU", i, j, au, pair_type=t)
        if multiloop:
            add(prefix+"_MLintern", i, j, baseline-au, pair_type=t)
        add(prefix+"_neighbor", i, j, full-baseline, pair_type=t,
            neighbor5=si, neighbor3=sj, guarded_neighbor5=gsi, guarded_neighbor3=gsj)
        terms[-1]["guarded_pf_energy"] = guarded-baseline

    for i in children[0]:
        stem(i, ends[i])
    for i, j in pp:
        branch = children[i]
        t = md.pair[encoded[i]][encoded[j]]
        if not branch:
            add("hairpin", i, j, energy(RNA.exp_E_Hairpin(j-i-1, t, encoded[i+1], encoded[j-1], sequence[i-1:j], p)))
        elif len(branch) == 1:
            k, l = branch[0], ends[branch[0]]
            inner_reverse = md.pair[encoded[l]][encoded[k]]
            value = energy(RNA.exp_E_IntLoop(k-i-1, j-l-1, t, inner_reverse,
                                           encoded[i+1], encoded[j-1], encoded[k-1], encoded[l+1], p))
            add("stack" if k == i+1 and l == j-1 else "bulge_or_internal", i, j, value,
                inner_i=k, inner_j=l)
        else:
            add("MLclosing", i, j, energy(p.expMLclosing))
            stem(i, j, multiloop=True, closing=True)
            unpaired = j-i-1-sum(ends[k]-k+1 for k in branch)
            add("MLbase", i, j, unpaired*energy(p.expMLbase), unpaired=unpaired)
            for k in branch:
                stem(k, ends[k], multiloop=True)
    return terms


def verify_cli_observer(inputs, provenance, input_dir):
    """Two small checks that the observation harness matches the adopted CLI."""
    executable = Path(provenance["build"])/"LinearCParty"
    parameter = Path(provenance["source"])/"params/rna_Turner04.par"
    before = digest(executable)
    results = []
    with tarfile.open(input_dir/"raw_outputs.tar.gz") as archive, \
            tempfile.TemporaryDirectory(prefix="energy-vienna-parity-", dir=os.environ.get("TMPDIR", "/tmp")) as temp:
        if str(Path(temp).resolve()).startswith("/gs/bs/"):
            raise ValueError("CLI scratch must be node-local")
        for identifier in ["random_000_n20", "exterior_adjacent_stems"]:
            if identifier not in inputs:
                continue
            record = inputs[identifier]
            stdout = archive.extractfile(identifier+"/pinned/stdout.tsv").read().decode()
            observed_bpp = {}
            for line in stdout.splitlines():
                f = line.split("\t")
                if f[0] == "BPP":
                    observed_bpp[int(f[1]), int(f[2])] = float(f[3])
                elif f[0] == "MFE":
                    observed_energy, observed_structure = float(f[1]), f[2]
            bpp_path = Path(temp)/"bpp.tsv"
            command = [str(executable), "--exact", "--sequence", record["sequence"],
                       "--scaffold", "."*record["length"], "--bpp", str(bpp_path), "--bpp-cutoff", "0"]
            run = subprocess.run(command, env=dict(os.environ, CPARTY_FIXED_ENERGY_PARAM_FILE=str(parameter)),
                                 text=True, capture_output=True, timeout=45)
            run.check_returncode()
            cli_bpp = {}
            for line in bpp_path.read_text().splitlines():
                if not line.strip() or line.startswith(">"):
                    continue
                i, j, p = line.split()
                cli_bpp[int(i), int(j)] = float(p)
            max_difference = max([abs(cli_bpp.get(k, 0.)-observed_bpp.get(k, 0.))
                                  for k in cli_bpp.keys() | observed_bpp.keys()]+[0.])
            cli_mfe = run.stdout.splitlines()[3]
            expected_mfe = f"{observed_structure} ({observed_energy:.2f})"
            row = dict(id=identifier, command=command, returncode=run.returncode,
                       stdout=run.stdout, stderr=run.stderr, bpp= bpp_path.read_text(),
                       max_abs_bpp_difference=max_difference, mfe_text_agrees=cli_mfe == expected_mfe)
            results.append(row)
            if max_difference > 1e-12 or not row["mfe_text_agrees"]:
                raise ValueError("observation harness disagrees with pinned CLI")
    return dict(executable=str(executable), executable_sha256=before,
                executable_unchanged=digest(executable) == before, checks=results)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT/"results/local/energy-vienna-check-20261017")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    out = args.output or args.input
    out.mkdir(parents=True, exist_ok=True)
    inputs = {r["id"]: r for r in map(json.loads, (args.input/"inputs.jsonl").read_text().splitlines())}
    provenance = json.loads((args.input/"provenance.json").read_text())
    cli_parity = verify_cli_observer(inputs, provenance, args.input)
    (out/"observer_cli_parity.json").write_text(json.dumps(cli_parity, indent=2)+"\n")
    param = Path(provenance["source"])/"params/rna_Turner04.par"
    if digest(param) != provenance["dependencies_sha256_before"][str(param)]:
        raise ValueError("parameter hash changed")
    RNA.params_load(str(param))
    rows = list(csv.DictReader((args.input/"structure_metrics.tsv").open(), delimiter="\t"))
    result, breakdown = [], []
    for index, row in enumerate(rows):
        if row["reachable"] != "True":
            continue
        sequence = inputs[row["id"]]["sequence"]
        terms = native_terms(sequence, row["structure"], int(row["dangles"]))
        native = sum(x["native_pf_energy"] for x in terms)
        guarded = sum(x["guarded_pf_energy"] for x in terms)
        result.append(dict(row=index, id=row["id"], group=row["group"], roles=row["roles"],
                           dangles=row["dangles"], length=row["length"], structure=row["structure"],
                           scaffold=row["scaffold"], native_local_sum=native, guarded_local_sum=guarded,
                           native_sum_vs_single_structure_pf=native-float(row["vienna_structure_pf_energy"]),
                           guarded_sum_vs_pinned=guarded-float(row["pinned_energy"]),
                           native_sum_vs_diagnostic=native-float(row["diagnostic_energy"]),
                           expected_boundary_difference=guarded-native,
                           observed_pinned_difference=float(row["delta_pinned_vs_vienna_pf"])))
        if row["roles"] == "constructed_target" and row["scaffold"] == "."*len(sequence):
            breakdown.extend(dict(row=index, id=row["id"], dangles=row["dangles"], **term) for term in terms)
    write_tsv(out/"structure_accounting.tsv", result)
    write_tsv(out/"term_breakdown.tsv", breakdown)
    ensemble = list(csv.DictReader((args.input/"ensemble_metrics.tsv").open(), delimiter="\t"))
    by_structure = {(r["id"], r["structure"]): r for r in rows
                    if r["dangles"] == "2" and r["scaffold"] == "."*int(r["length"])}
    mfe_comparisons = []
    for row in ensemble:
        if row["status"] != "success":
            continue
        reference = by_structure[row["id"], row["vienna_mfe_structure"]]
        candidate = by_structure[row["id"], row["mfe_structure"]]
        key = "pinned_energy" if row["variant"] == "pinned" else "diagnostic_energy"
        gap = float(reference[key])-float(candidate[key])
        pf_gap = float(reference["vienna_structure_pf_energy"])-float(candidate["vienna_structure_pf_energy"])
        mfe_comparisons.append(dict(id=row["id"], group=row["group"], variant=row["variant"], length=row["length"],
            same_pair_set=row["mfe_same_pair_set"],
            lcp_score_gap_vienna_mfe_minus_lcp_mfe=gap,
            independent_pf_gap_vienna_mfe_minus_lcp_mfe=pf_gap,
            integer_gap_lcp_mfe_minus_vienna_mfe=round(float(candidate["vienna_eval_energy"]), 2)-float(row["vienna_mfe_energy"]),
            lcp_mfe_vs_its_fixed_structure_evaluation=float(row["mfe_energy"])-float(candidate[key]),
            classification="same_structure" if row["mfe_same_pair_set"] == "True" else
                           "same_score_tie" if abs(gap) < 1e-9 and abs(pf_gap) < 1e-9 else
                           "different_scores_inspect_pf_and_boundary_terms"))
    write_tsv(out/"mfe_comparison.tsv", mfe_comparisons)
    witnesses = []
    special_rows = {max(range(len(rows)), key=lambda k: abs(
        round(float(rows[k]["vienna_eval_energy"]), 2)-float(rows[k]["vienna_structure_pf_energy"]))),
        max(range(len(rows)), key=lambda k: abs(float(rows[k]["delta_pinned_vs_vienna_pf"])))}
    for index in sorted(special_rows):
        row = rows[index]
        sequence = inputs[row["id"]]["sequence"]
        witnesses.append(dict(row=index, sequence=sequence, **row,
                              terms=native_terms(sequence, row["structure"], int(row["dangles"]))))
    (out/"accounting_witnesses.json").write_text(json.dumps(witnesses, indent=2)+"\n")
    summary = dict(evaluations=len(result),
                   max_native_local_sum_vs_independent_pf=max(abs(r["native_sum_vs_single_structure_pf"]) for r in result),
                   max_guarded_local_sum_vs_pinned=max(abs(r["guarded_sum_vs_pinned"]) for r in result),
                   max_native_local_sum_vs_diagnostic=max(abs(r["native_sum_vs_diagnostic"]) for r in result),
                   changed_by_boundary_guard=sum(abs(r["expected_boundary_difference"]) > 1e-9 for r in result),
                   max_abs_boundary_difference=max(abs(r["expected_boundary_difference"]) for r in result),
                   native_eval_vs_pf_max=max(abs(float(r["vienna_eval_energy"])-float(r["vienna_structure_pf_energy"])) for r in rows),
                   native_integer_eval_vs_pf_max=max(abs(round(float(r["vienna_eval_energy"]), 2)-float(r["vienna_structure_pf_energy"])) for r in rows),
                   mfe_classification={variant: {key: sum(r["classification"] == key and r["variant"] == variant for r in mfe_comparisons)
                        for key in sorted({r["classification"] for r in mfe_comparisons})}
                        for variant in ["pinned", "inner_neighbors_diagnostic"]},
                   max_mfe_vs_its_fixed_structure_evaluation=max(abs(r["lcp_mfe_vs_its_fixed_structure_evaluation"]) for r in mfe_comparisons),
                   source_parameter_sha256=digest(param),
                   environment=dict(python_version=sys.version, platform=platform.platform(),
                                    vienna_version=RNA.__version__,
                                    vienna_native_extension=str(native_rna_extension.__file__),
                                    vienna_native_extension_sha256=digest(Path(native_rna_extension.__file__))),
                   input_hashes={name: digest(args.input/name) for name in
                       ["inputs.jsonl", "ensemble_metrics.tsv", "structure_metrics.tsv"]},
                   script_sha256=digest(Path(__file__)))
    (out/"accounting_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))
    if max(summary[k] for k in ["max_native_local_sum_vs_independent_pf", "max_guarded_local_sum_vs_pinned",
                                "max_native_local_sum_vs_diagnostic"]) > 1e-9:
        raise SystemExit("unexplained fixed-structure residual; inspect structure_accounting.tsv")


if __name__ == "__main__":
    main()
