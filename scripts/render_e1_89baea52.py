#!/usr/bin/env python3
"""Render only the adopted E1 figures with the unchanged 20261017 layout.

No RNA programs or archives are opened. Validate quantiles against the new
RNA-level metrics, and preserve source and output hashes beside the figures.
"""

from collections import defaultdict
import json
import os
from pathlib import Path

import make_figures as base
import make_figures_review as review


def main():
    out = base.ROOT / "results/figures/e1_89baea52_20261019"
    out.mkdir(parents=True, exist_ok=False)
    data = base.Data.__new__(base.Data)
    data.paths = {}
    data.audit = []
    data.checks = defaultdict(int)
    data.q = {}
    data.thresholds = {}
    data.figure3_only = False
    revisions = {}
    for cohort, suffix in (("short", "20261018"), ("long", "cpu16-20261018")):
        directory = base.ROOT / f"results/local/e1-bprna-{cohort}-89baea52-{suffix}"
        provenance = json.loads((directory / "provenance.json").read_text())
        revisions[cohort] = [run["software_commit"] for run in provenance["runs"]]
        base.require(all(r.startswith("89baea52") for r in revisions[cohort]),
                     "E1 sources must all be 89baea52")
        for name in ("metrics.tsv", "quantiles.tsv", "provenance.json", "window_counts.tsv"):
            data.paths[f"e1_{cohort}_{name}"] = directory / name
        base.require(base.digest(directory / "metrics.tsv") == provenance["metrics_sha256"],
                     "E1 metrics differ from audited provenance")
    data.hashes = {key: base.digest(path) for key, path in data.paths.items()}
    data.load_e1()
    review.style(base)
    counts = review.evaluation_counts(data, base, out)
    records = []
    for draw, name, width in ((review.accuracy_figure, "figure2_main", 178),
                              (review.threshold_figure, "figure2_d_si", 165),
                              (review.accessibility_detail, "figure2_accessibility_si", 165)):
        base.save(draw(data, base), name, out, width, records, check_text_overlap=True)
        print(f"Rendered {name}", flush=True)
    pdf_checks = {r["file"]: review.pdf_text_check(out / r["file"])
                  for r in records if r["file"].endswith(".pdf")}
    base.write_audit(data, out)
    # Reuse the paired E1 caption definitions, excluding unrelated E2 captions.
    review.captions(data, base, out, counts)
    blocks = (out / "caption.md").read_text().split("\n## ")
    selected = [block for block in blocks[1:]
                if block.startswith(("Figure 2", "Figure S1", "Additional SI — figure2_"))]
    selected.sort(key=lambda block: (0 if block.startswith("Figure 2") else
                                     1 if "figure2_accessibility_si.pdf" in block.splitlines()[0] else 2))
    caption = blocks[0] + "".join("\n## " + block for block in selected)
    caption = caption.replace("Figure S1 — figure2_d_si.pdf", "Figure S2 — figure2_d_si.pdf")
    caption = caption.replace("Additional SI — figure2_accessibility_si.pdf", "Figure S1 — figure2_accessibility_si.pdf")
    caption = caption.replace("Figure numbers follow the task labels; the manuscript editor can renumber the resource figure after Terc.",
                              "Figure numbers match the adopted E1 manuscript figures.")
    (out / "caption.md").write_text(caption)
    for key, path in data.paths.items():
        base.require(base.digest(path) == data.hashes[key], f"Source changed: {path}")
    validation = dict(all_passed=True, checks=dict(data.checks),
                      source_files_unchanged=True, text_within_pages=True,
                      text_without_overlap=True, plotted_points_within_axis_limits=True,
                      minimum_normal_font_pt=7, pdf_text_checks=pdf_checks, outputs=records)
    (out / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    (out / "README.md").write_text("""# E1 accuracy figures — 89baea52, 2026-10-19

Finite-beam and exact outputs both use 89baea52. Long RNAs use the completed
`e1-bprna-long-89baea52-cpu16-20261018` run and the finite b=400 reference.
RNA cohorts, metric definitions, linear-interpolation quantiles and the
`main_20261017/figure2_main` layout are unchanged.

Run from the experiments repository:
`python3 -B scripts/render_e1_89baea52.py`
The output directory must not exist. The renderer does not change old figures.

`figure2_main` is 178 mm wide. The SI companions `figure2_accessibility_si`
(all four window lengths) and `figure2_d_si` (positive delta F threshold)
are 165 mm wide. Each has a PDF and 300-dpi PNG.
`caption.md` contains paired Japanese/English captions.
`plot_data.tsv` records all plotted values and `evaluation_counts.tsv` the
eligible RNA/window denominators. `validation.json` records checks against
the source quantiles, geometry, overlap and actual PDF font sizes.
`manifest.json` records sources, scripts and output hashes.
PNG inspection and manuscript integration are recorded in `review_report.json`.
""")
    derived = ("caption.md", "README.md", "plot_data.tsv", "evaluation_counts.tsv", "validation.json")
    manifest = dict(layout="review-20261017", e1_revisions=revisions,
                    scripts_sha256={p.name: base.digest(p) for p in
                                    (Path(__file__), Path(base.__file__), Path(review.__file__))},
                    numpy=base.np.__version__, matplotlib=base.plt.matplotlib.__version__,
                    sources={key: dict(path=os.path.relpath(path, out), sha256=data.hashes[key])
                             for key, path in data.paths.items()}, outputs=records,
                    derived_files={name: base.digest(out / name) for name in derived})
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(validation["checks"], indent=2))


if __name__ == "__main__":
    main()
