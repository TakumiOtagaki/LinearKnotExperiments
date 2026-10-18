#!/usr/bin/env python3
"""Recompute NEAT1 predictions and by-loop checks; impose no old-output equality."""
import argparse
import json
from pathlib import Path

import numpy as np
import analyze_e3_neat1 as core
import analyze_e3_neat1_rerun as audit
from prepare_paper_rerun_89baea52 import ROOT, COMMIT
from run_cli import digest, now


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=ROOT / "results/local/e3-neat1-89baea52-20261018")
    out = parser.parse_args().out.resolve()
    out.mkdir(exist_ok=False)
    provenance = dict(created_utc=now(), release=COMMIT, inputs=[], checkpoints=[], artifacts=[],
                      script_sha256=digest(Path(__file__)))
    inputs = core.read_inputs(provenance)
    inventory = []
    saved = core.read_runs(inputs, provenance, inventory, "89baea52-20261018")
    assert len(inventory) == 51 and all(r["status"] == "success" for r in inventory)
    validation = {"conditions": []}
    core.parse_saved(saved, provenance, validation)
    access, differences, regions = core.access_analysis(saved)
    sampling, hist, distances, reps, checks = core.sampling_analysis(saved)
    validation["sampling"] = checks
    for name, rows in (("inventory", inventory), ("accessibility_summary", access),
                       ("accessibility_differences", differences), ("accessibility_regions", regions),
                       ("sampling_summary", sampling), ("sampling_histograms", hist),
                       ("representative_distances", reps)):
        core.write_table(out / (name+".tsv"), rows)
    core.write_table(out / "sampling_reference_distances.tsv.gz", distances)
    shape_placeholder = np.full(inputs["NEAT1_1"]["n"], np.nan)
    profiles = []
    for row in core.profile_rows(saved, shape_placeholder):
        row.pop("shape", None)
        profiles.append(row)
    core.write_table(out / "position_profiles.tsv.gz", profiles)
    core.write_table(out / "accessibility.tsv.gz", (dict(iso=k[0], method=k[1], beam=k[2],
        start=i+1, width=w, accessibility=p) for k,g in sorted(saved.items())
        for w,values in g.get("access", {}).items() for i,p in enumerate(values)),
        ["iso", "method", "beam", "start", "width", "accessibility"])
    core.write_table(out / "bpp_pairs.tsv.gz", (dict(iso=k[0], method=k[1], beam=k[2],
        i=i, j=j, probability=p, is_scaffold=int((i,j) in g["input"]["fixed"]),
        crosses_G=int((i,j) in g.get("cross_pairs", {})))
        for k,g in sorted(saved.items()) for (i,j),p in sorted(g.get("pairs", {}).items())),
        ["iso", "method", "beam", "i", "j", "probability", "is_scaffold", "crosses_G"])
    # Keep the plot's stdout/BPP records (which include saved commands)
    # separate from the full resource audit, which also inventories other members.
    resource_provenance = {"checkpoints": [], "artifacts": []}
    new = audit.load("89baea52-20261018", inputs, resource_provenance)
    provenance["resource_audit"] = resource_provenance
    resource_rows, components, focus, residuals = [], [], [], []
    for (iso, beam, kind), g in new.items():
        row = g["row"]
        resource_rows.append(dict(iso=iso, beam=beam, kind=kind,
            elapsed_seconds=row["elapsed_seconds"], time_wall_seconds=row["time_wall_seconds"],
            peak_rss_kib=row["peak_rss_kib"], peak_rss_gib=row["peak_rss_kib"]/2**20))
        if kind != "raccess-byloop":
            continue
        total, parts = audit.parse_byloop(g["raw"]["stdout.txt"], inputs[iso]["n"])
        plain = core.parse_access(new[iso, beam, "raccess"]["raw"]["stdout.txt"], inputs[iso])
        contexts = core.parse_context(new[iso, beam, "capr"]["raw"]["stdout.txt"], inputs[iso])
        residual = max(float(np.max(np.abs(sum(parts[w].values())-total[w]))) for w in core.WIDTHS)
        diff = max(float(np.max(np.abs(total[w]-plain[w]))) for w in core.WIDTHS)
        unit = max(float(np.max(np.abs(parts[1][c]-contexts[name]))) for c,name in audit.LOOPS.items())
        assert max(residual, diff, unit) < 1e-8
        residuals.append(dict(iso=iso, beam=beam, sum_vs_total=residual, total_vs_plain=diff, unit_vs_context=unit))
        for w in core.WIDTHS:
            for c,name in audit.LOOPS.items():
                values = parts[w][c]; i = int(np.argmax(values))
                components.append(dict(iso=iso, beam=beam, window=w, component=name,
                    mean=float(values.mean()), max=float(values[i]), max_start_1based=i+1))
        if iso == "NEAT1_2":
            chosen = slice(21134-1, 21143)
            focus.append(dict(iso=iso, beam=beam, window=8, start_min=21134, start_max=21143,
                total_mean=float(total[8][chosen].mean()),
                **{name+"_mean": float(parts[8][c][chosen].mean()) for c,name in audit.LOOPS.items()}))
    for name, rows in (("resources", resource_rows), ("byloop_components", components), ("byloop_regions", focus)):
        core.write_table(out / (name+".tsv"), rows)
    validation.update(byloop=residuals, status="passed", probability_outputs_recomputed=True)
    core.write_json(out / "validation.json", validation)
    core.write_json(out / "provenance.json", provenance)
    (out / "README.md").write_text(
        "# NEAT1: fresh 89baea52 analysis\n\n"
        "All 51 main/variant conditions and the six additional by-loop conditions were recomputed.\n"
        "`inventory.tsv` preserves the condition count; `resources.tsv` gives all main/by-loop resource values.\n"
        "Accessibility, profiles, BPPs, sampling summaries and distances use only the new saved outputs.\n"
        "`validation.json` gives probability and by-loop checks.\n"
        "`provenance.json` records input/checkpoint/artifact hashes and commands; "
        "its separate `resource_audit` inventories the additional resource checks.\n"
        "Raw tar members are streamed without extraction.\n\n"
        "Reproduce with `python3 -B scripts/analyze_neat1_89baea52.py --out NEW_DIRECTORY`.\n"
        "Render with `python3 -B scripts/render_non_e1_89baea52.py`.\n")


if __name__ == "__main__":
    main()
