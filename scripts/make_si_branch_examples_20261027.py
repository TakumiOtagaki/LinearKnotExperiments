#!/usr/bin/env python3
"""Draw SI branch examples and inside/outside flow, without RNA calculations.

Run from the experiment repository:
  MPLCONFIGDIR=/tmp/si-branch-mpl python3 scripts/make_si_branch_examples_20261027.py \
    --out results/figures/si_branch_examples_20261027_v4 \
    --context-reference ../../LinearKnot/test/tools/pseudoloop_context_reference.py

The optional context reference is read only. Pair colors and line styles match
make_figure1.py. The examples are fictional structures, not measured results.
"""
import argparse
import hashlib
import json
from pathlib import Path
import runpy

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Arc, Polygon, Rectangle

MM = 1 / 25.4
INK, GREY, BLUE, ORANGE, PURPLE = (
    "#1f1f1e", "#8a8a85", "#2a78d6", "#eb6834", "#6250d6")
STYLES = ("-", (0, (1, 1)), (0, (4, 2)))
EXAMPLES = [
    dict(name="Ordinary branch", n=20, G=[], extension=[(3, 17), (4, 16), (7, 12)],
         branch=(3, 17), branch_label="one ordinary branch"),
    dict(name="Pseudoknot branch", n=25, G=[(2, 15), (3, 14)],
         extension=[(5, 9), (10, 23), (11, 22)], branch=(2, 23),
         branch_label="one pseudoknot branch"),
    dict(name="Band-spanning multiloop", n=34, G=[(2, 27), (14, 23)],
         extension=[(5, 10), (18, 32), (19, 31)],
         outer=(2, 27), inner=(14, 23), side=(5, 10),
         unpaired_segments=[3, 4, 11, 12, 13, 24, 25, 26], context="M"),
    dict(name="Internal loop without a side branch", n=34,
         G=[(2, 27), (14, 23)], extension=[(18, 32), (19, 31)],
         outer=(2, 27), inner=(14, 23),
         unpaired_segments=list(range(3, 14)) + [24, 25, 26], context="I"),
    dict(name="Band-spanning multiloop closed by extension pairs", n=34,
         G=[(18, 32), (19, 31)], extension=[(2, 27), (14, 23), (5, 10)],
         outer=(2, 27), inner=(14, 23), side=(5, 10),
         unpaired_segments=[3, 4, 11, 12, 13, 24, 25, 26], context="M"),
]


def crosses(p, q):
    i, j = p
    k, ell = q
    return i < k < j < ell or k < i < ell < j


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate(reference):
    result = {"fictional_structures": True, "folding_runs": 0,
              "pair_colors_match_figure1": True, "examples": []}
    for ex in EXAMPLES:
        pairs = ex["G"] + ex["extension"]
        endpoints = [x for pair in pairs for x in pair]
        assert len(endpoints) == len(set(endpoints))
        assert all(1 <= i < j <= ex["n"] and j - i > 3 for i, j in pairs)
        for group in (ex["G"], ex["extension"]):
            assert not any(crosses(p, q) for p in group for q in group)
        if "branch" in ex:
            lo, hi = ex["branch"]
            assert all(not ((lo <= i <= hi) != (lo <= j <= hi)) for i, j in ex["G"])
        if "outer" in ex:
            closing_in_scaffold = ex["outer"] in ex["G"]
            closing_layer = ex["G"] if closing_in_scaffold else ex["extension"]
            other_layer = ex["extension"] if closing_in_scaffold else ex["G"]
            assert ex["inner"] in closing_layer
            assert any(crosses(ex["outer"], p) and crosses(ex["inner"], p)
                       for p in other_layer)
            assert ("side" in ex) == (ex["context"] == "M")
        result["examples"].append({"name": ex["name"], "valid_pair_map": True})
    if reference:
        classify = runpy.run_path(str(reference))["classify"]
        for ex, rec in zip(EXAMPLES, result["examples"]):
            labels, _ = classify(ex["n"], set(ex["G"] + ex["extension"]))
            rec["context_labels"] = labels
            if "unpaired_segments" in ex:
                assert all(labels[i - 1] == ex["context"] for i in ex["unpaired_segments"])
                rec["all_marked_unpaired_bases_match_context"] = ex["context"]
        result["context_reference_sha256"] = digest(reference)
    return result


def arc_points(pair, reverse=False):
    """Points along the same upper half-ellipse used for the pair arcs."""
    i, j = pair
    theta = np.linspace(np.pi, 0, 121)
    points = np.column_stack(((i + j) / 2 + (j - i) / 2 * np.cos(theta),
                              .45 * (j - i) * np.sin(theta)))
    return points[::-1] if reverse else points


def diagram(fig, ex, box):
    width, height = fig.get_size_inches() / MM
    x, y, w, h = box
    ax = fig.add_axes([x / width, y / height, w / width, h / height])
    max_arc_height = .45 * max(j - i for i, j in ex["G"] + ex["extension"])
    ax.set(xlim=(0, ex["n"] + 1), ylim=(-1, 1.5 * max_arc_height))
    ax.axis("off")
    if "outer" in ex:
        # The shaded face excludes both the continuing inner region and the
        # side branch's interior. Crossing arcs belong to the other layer.
        lower_pairs = sorted([ex["inner"]] + ([ex["side"]] if "side" in ex else []),
                             reverse=True)
        boundary = np.vstack([arc_points(ex["outer"])] +
                             [arc_points(pair, reverse=True) for pair in lower_pairs])
        ax.add_patch(Polygon(boundary, closed=True, fc=PURPLE, alpha=.12, ec="none", zorder=0))
    ax.plot([1, ex["n"]], [0, 0], lw=.6, color=GREY)
    for i, j in sorted(ex["G"] + ex["extension"], key=lambda p: p[1] - p[0], reverse=True):
        kind = 0 if (i, j) in ex["G"] else 2 if any(crosses((i, j), q) for q in ex["G"]) else 1
        ax.add_patch(Arc(((i + j) / 2, 0), j - i, .9 * (j - i),
                         theta1=0, theta2=180, color=(GREY, BLUE, ORANGE)[kind],
                         lw=1.15, ls=STYLES[kind]))
    if "unpaired_segments" in ex:
        ax.scatter(ex["unpaired_segments"], [0] * len(ex["unpaired_segments"]),
                   s=7 if w > 30 else 4, color=PURPLE, zorder=4)
    return ax


def branches():
    width, height = 165, 65
    fig = plt.figure(figsize=(width * MM, height * MM))
    # Give each base the same physical spacing and preserve the arc geometry.
    # Stack A/B to leave room for C and its existing upper/lower comparisons.
    base_spacing = 2  # mm per base at the 165 mm paper inclusion width
    for k, ex in enumerate(EXAMPLES[:3]):
        x, baseline, title_y = ((2, 42, 59), (2, 15, 30), (60, 29, 59))[k]
        fig.text(x / width, title_y / height, f"{'ABC'[k]}  {ex['name']}", size=7.6, weight="bold")
        panel_width = base_spacing * (ex["n"] + 1)
        ymax = 1.5 * .45 * max(j - i for i, j in ex["G"] + ex["extension"])
        diagram_height = base_spacing * (ymax + 1)
        ax = diagram(fig, ex, (x, baseline - base_spacing, panel_width, diagram_height))
        ax.set_aspect("equal", adjustable="box", anchor="S")
        if "branch" in ex:
            lo, hi = ex["branch"]
            # Use axes coordinates to keep the bracket below every arc.
            left, right = lo / (ex["n"] + 1), hi / (ex["n"] + 1)
            ax.plot([left, left, right, right], [-.04, -.12, -.12, -.04],
                    transform=ax.transAxes, color=INK, lw=.75, clip_on=False)
            ax.text((left + right) / 2, -.16, ex["branch_label"],
                    transform=ax.transAxes, ha="center", va="top", size=7)
        else:
            for label, target, pos in [
                    ("Outer pair", (13, 11.1), (7, 13)),
                    ("Inner pair", (22, 2.55), (27, 7.5)),
                    ("Side branch", (7.5, 2.25), (3, 11.5))]:
                ax.annotate(label, xy=target, xytext=pos, fontsize=6.6, ha="center",
                            arrowprops=dict(arrowstyle="-", lw=.55, color=INK))
            ax.annotate("Band-spanning multiloop\n(this shaded loop)", xy=(12.5, 5.8),
                        xytext=(.47, -.03), textcoords="axes fraction", ha="center", va="top",
                        fontsize=7, arrowprops=dict(arrowstyle="-", lw=.6, color=INK))
    fig.text(148 / width, 52 / height, "No side branch", size=7, ha="center")
    diagram(fig, EXAMPLES[3], (134, 34, 29, 18))
    fig.text(148 / width, 31.5 / height, "Internal loop", size=7, ha="center")
    fig.text(148 / width, 25 / height, "Closing pairs in $G'$", size=7, ha="center")
    diagram(fig, EXAMPLES[4], (134, 7.5, 29, 18))
    fig.text(148 / width, 5 / height, "Multiloop", size=7, ha="center")
    fig.text(60 / width, 7 / height, "Purple: unpaired segments\nof this loop", color=PURPLE,
             size=7, va="bottom")
    handles = [Line2D([], [], color=c, ls=ls, lw=1.2) for c, ls in zip((GREY, BLUE, ORANGE), STYLES)]
    fig.legend(handles, ["Scaffold pair", "New noncrossing pair", "New crossing pair"],
               loc="lower center", bbox_to_anchor=(.5, .01), ncol=3, frameon=False,
               handlelength=2.6, columnspacing=2, fontsize=7)
    return fig


def tree():
    fig, axes = plt.subplots(1, 2, figsize=(152 * MM, 54 * MM))
    fig.subplots_adjust(left=.03, right=.99, top=.99, bottom=.03, wspace=.14)
    nodes = {"root": (.5, .86), "v": (.32, .56), "x": (.79, .56),
             "u1": (.13, .16), "u2": (.51, .16)}
    rules = [("root", "r0", (.5, .71), ("v", "x")),
             ("v", "r", (.32, .37), ("u1", "u2"))]
    for k, ax in enumerate(axes):
        ax.set(xlim=(0, 1), ylim=(0, 1))
        ax.axis("off")
        ax.text(.01, .98, ("A  Inside: children to root", "B  Outside: root to children")[k],
                size=8, weight="bold", va="top")
        for parent, rule, pos, children in rules:
            for start, end in [(nodes[parent], pos)] + [(pos, nodes[c]) for c in children]:
                if k == 0:
                    start, end = end, start
                ax.annotate("", xy=end, xytext=start,
                            arrowprops=dict(arrowstyle="-|>", lw=.85, color=INK,
                                            linestyle="-" if k == 0 else "--",
                                            shrinkA=6, shrinkB=6, mutation_scale=8))
            ax.add_patch(Rectangle((pos[0] - .028, pos[1] - .025), .056, .05,
                                   fc="white", ec=INK, lw=.7, zorder=4))
            ax.text(pos[0] + .045, pos[1], r"$\omega_0$" if rule == "r0" else r"$\omega_r$",
                    size=8, va="center")
        for name, pos in nodes.items():
            label = {"root": r"$v_0$", "v": r"$v$", "x": r"$x$",
                     "u1": r"$u_1$", "u2": r"$u_2$"}[name]
            ax.text(*pos, label, ha="center", va="center", size=10,
                    bbox=dict(boxstyle="round,pad=.25", fc="white", ec=GREY, lw=.6), zorder=5)
        if k == 0:
            ax.text(.73, .85, r"$\alpha(v_0)=\hat Z$", size=8, ha="center")
            ax.text(.32, .035, r"$\omega_r\alpha(u_1)\alpha(u_2)\ \to\ \alpha(v)$",
                    ha="center", size=8)
        else:
            ax.text(.75, .85, r"$\beta(v_0)=1$", size=8, ha="center")
            ax.text(.32, .035, r"$\beta(v)\omega_r\alpha(u_2)\ \to\ \beta(u_1)$",
                    ha="center", size=8)
            ax.text(.77, .31, "Use sibling\ninside weight", size=7, ha="center")
            ax.annotate("", xy=(.49, .19), xytext=(.72, .27),
                        arrowprops=dict(arrowstyle="-", lw=.55, color=GREY))
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--context-reference", type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7,
                         "text.color": INK, "pdf.fonttype": 42, "mathtext.fontset": "dejavusans"})
    validation = validate(args.context_reference)
    validation["layout"] = {}
    for name, create in [("si_branch_examples", branches), ("si_inside_outside", tree)]:
        fig = create()
        fig.canvas.draw()
        labels = (fig.texts + [t for ax in fig.axes for t in ax.texts] +
                  [t for legend in fig.legends for t in legend.get_texts()])
        renderer = fig.canvas.get_renderer()
        for label in labels:
            box = label.get_window_extent(renderer)
            assert 0 <= box.x0 <= box.x1 <= fig.bbox.width, label.get_text()
            assert 0 <= box.y0 <= box.y1 <= fig.bbox.height, label.get_text()
        minimum_font = min(label.get_fontsize() for label in labels)
        assert minimum_font >= 6
        validation["layout"][name] = {
            "size_mm": [round(float(v), 3) for v in fig.get_size_inches() / MM],
            "minimum_font_pt": minimum_font, "all_labels_inside_canvas": True}
        if name == "si_branch_examples":
            scales = []
            for panel, ax in zip("ABC", fig.axes[:3]):
                origin, dx, dy = ax.transData.transform([(0, 0), (1, 0), (0, 1)])
                xscale = (dx[0] - origin[0]) / fig.dpi / MM
                yscale = (dy[1] - origin[1]) / fig.dpi / MM
                scales.append((xscale, yscale))
                validation["layout"][name].setdefault("main_panels", {})[panel] = {
                    "base_spacing_mm": round(float(xscale), 6),
                    "vertical_unit_mm": round(float(yscale), 6),
                    "arc_height_over_endpoint_distance": round(float(.45 * yscale / xscale), 6)}
            assert np.allclose(scales, scales[0][0]), "Main panel scales must match"
        fig.savefig(args.out / (name + ".pdf"))
        fig.savefig(args.out / (name + ".png"), dpi=220)
        plt.close(fig)
    (args.out / "examples.json").write_text(json.dumps(EXAMPLES, indent=2) + "\n")
    (args.out / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    manifest = {"script": Path(__file__).name, "script_sha256": digest(Path(__file__)),
                "figure1_pair_colors": {"scaffold": GREY, "new_noncrossing": BLUE,
                                         "new_crossing": ORANGE, "Multiloop": PURPLE},
                "files": {p.name: digest(p) for p in sorted(args.out.glob("*"))
                          if p.suffix in {".pdf", ".png", ".json"} and p.name != "manifest.json"}}
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
