#!/usr/bin/env python3
"""Draft Figure 1: what LinearCParty conditions on, approximates and computes.

Panel A shows the scaffold G, an extension G' and the completed structure R for
a short pseudoknotted RNA (bpRNA_PDB_345; G = nested pairs of its reference).
Panel B is a schematic of the rearranged decomposition and beam pruning.
Panel C shows real exact outputs for the same RNA, computed by the CLIs listed
in --data (BPP, seven structural contexts, window accessibility, Boltzmann
samples, MEA, MFE and suboptimal structures).

  python3 scripts/make_figure1.py --data <dir with lc.out bpp.tsv capr.out
      racc.out samp.out mea.out subopt.out rec.json> --out <output dir>
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, FancyArrowPatch, FancyBboxPatch, Rectangle
import matplotlib.lines

MM = 1 / 25.4
INK, MUTED, GRID = "#1f1f1e", "#6b6b66", "#d9d8d3"
G_COLOR, NEW_COLOR, CROSS_COLOR = "#8a8a85", "#2a78d6", "#eb6834"
CONTEXTS = ["Stem", "Hairpin", "Bulge", "Internal", "Multiloop", "Exterior", "PseudoLoop"]
CONTEXT_COLORS = {"Stem": "#2a78d6", "Hairpin": "#1baf7a", "Bulge": "#eda100", "Internal": "#e87ba4",
                  "Multiloop": "#6250d6", "Exterior": "#b9b8b1", "PseudoLoop": "#eb6834"}

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.5, "axes.linewidth": 0.6,
                     "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED,
                     "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2,
                     "ytick.major.size": 2, "pdf.fonttype": 42})


def pairs_from_db(db):
    stacks = {"(": [], "[": [], "{": []}
    close = {")": "(", "]": "[", "}": "{"}
    out = []
    for i, c in enumerate(db, 1):
        if c in stacks:
            stacks[c].append(i)
        elif c in close:
            out.append((stacks[close[c]].pop(), i, c))
    return out


def crosses(a, b):
    (i, j), (k, l) = a, b
    return i < k < j < l or k < i < l < j


def draw_arcs(ax, pairs, n, color, alpha=1.0, lw=0.9, above=True):
    for p in pairs:
        i, j = p[0], p[1]
        a = alpha[(i, j)] if isinstance(alpha, dict) else alpha
        if a <= 0.01:
            continue
        w = j - i
        ax.add_patch(Arc(((i + j) / 2, 0), w, w * 0.9, theta1=0 if above else 180,
                         theta2=180 if above else 360, color=color, lw=lw, alpha=min(1, a)))


def setup_arc_axis(ax, n, height, label=None):
    ax.set_xlim(0, n + 1)
    ax.set_ylim(-0.04 * height, height)
    ax.axhline(0, color=MUTED, lw=0.6)
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.set_xticks([])
    if label:
        ax.text(0.0, 1.0, label, transform=ax.transAxes, ha="left", va="top", fontsize=6.5, color=INK)


def panel_label(fig, x, y, s):
    fig.text(x, y, s, fontsize=9, fontweight="bold", color=INK, ha="left", va="top")


def read_outputs(d):
    rec = json.loads((d / "rec.json").read_text())
    seq = rec["seq"]
    n = len(seq)
    G = [(i, j) for i, j in rec["nested"]]
    Gp = [(i, j) for i, j in rec["pk"]]
    bpp = {}
    for line in (d / "bpp.tsv").read_text().splitlines():
        f = line.split()
        if len(f) == 3 and f[0].isdigit():
            bpp[(int(f[0]), int(f[1]))] = float(f[2])
    ctx = {}
    for line in (d / "capr.out").read_text().splitlines():
        f = line.split()
        if f and f[0] in CONTEXTS:
            ctx[f[0]] = [float(x) for x in f[1:]]
    acc = {}
    for line in (d / "racc.out").read_text().splitlines():
        if line.startswith(("#", ">")) or not line.strip():
            continue
        pos, rest = line.split("\t")
        for item in rest.strip(";\n").split(";"):
            w, p = item.split(",")
            acc.setdefault(int(w), {})[int(pos) + 1] = float(p)
    samples = []
    for line in (d / "samp.out").read_text().splitlines():
        f = line.split("\t")
        if len(f) == 2 and f[0].isdigit():
            samples.append((int(f[0]), f[1]))
    mea = [l.split()[0] for l in (d / "mea.out").read_text().splitlines() if "EA=" in l][0]
    subopt = []
    for line in (d / "subopt.out").read_text().splitlines():
        f = line.split()
        if len(f) == 2 and f[1].startswith("(") and set(f[0]) <= set(".()[]{}"):
            subopt.append((f[0], float(f[1].strip("()"))))
    return seq, n, G, Gp, bpp, ctx, acc, samples, mea, subopt


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    seq, n, G, Gp, bpp, ctx, acc, samples, mea, subopt = read_outputs(args.data)
    args.out.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(183 * MM, 150 * MM))
    T = fig.transFigure

    def box(ax, x, y, w, hh, text, fc="white", ec=MUTED, fs=6.0, lw=0.6, ls="-"):
        ax.add_patch(FancyBboxPatch((x, y), w, hh, boxstyle="round,pad=0.25,rounding_size=1.0",
                                    fc=fc, ec=ec, lw=lw, ls=ls))
        ax.text(x + w / 2, y + hh / 2, text, ha="center", va="center", fontsize=fs, color=INK)

    def arrow(ax, x0, y0, x1, y1, color=MUTED, lw=0.8, ls="-"):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=6,
                                     color=color, lw=lw, ls=ls, shrinkA=0, shrinkB=0))

    # ---------------- Panel A: what we condition on ----------------
    panel_label(fig, 0.005, 0.995, "A")
    fig.text(0.035, 0.993, "What is conditioned on", fontsize=7.5, fontweight="bold", va="top")
    rows = [("Scaffold $G$: input pairs without crossings", [G], [G_COLOR]),
            ("Extension $G'$: added pairs without crossings,\nsharing no base with $G$", [Gp], [CROSS_COLOR]),
            ("Completed structure $R = G \\cup G'$", [G, Gp], [G_COLOR, CROSS_COLOR])]
    ytop = 0.95
    for k, (label, sets, colors) in enumerate(rows):
        yl = ytop - k * 0.118
        fig.text(0.035, yl, label, fontsize=6.2, color=INK, va="top", linespacing=1.25)
        ax = fig.add_axes([0.035, yl - 0.108, 0.30, 0.07])
        setup_arc_axis(ax, n, 34)
        for pset, col in zip(sets, colors):
            draw_arcs(ax, pset, n, col, lw=0.9)
        ax.scatter(range(1, n + 1), [0] * n, s=1.0, color=MUTED, zorder=3, lw=0)
    fig.text(0.035, ytop - 3 * 0.118 + 0.002,
             "The ensemble holds every $R$ that keeps $G$:\n"
             "$P(R \\mid S, G) = w(R)\\,/\\,Z(S, G)$, where $w(R)$ is the Boltzmann\n"
             "weight and $Z(S, G)$ sums $w(R)$ over all such $R$.",
             fontsize=6.0, color=INK, va="top", linespacing=1.35)

    # ---------------- Panel B: what we approximate ----------------
    panel_label(fig, 0.375, 0.995, "B")
    fig.text(0.405, 0.993, "What is approximated", fontsize=7.5, fontweight="bold", va="top")
    axb = fig.add_axes([0.375, 0.575, 0.62, 0.385])
    axb.set_xlim(0, 100)
    axb.set_ylim(0, 60)
    axb.axis("off")
    # step 1: CParty
    axb.text(12, 59, "CParty", fontsize=6.4, fontweight="bold", ha="center", va="top")
    box(axb, 2, 46, 20, 5, "region $[i, j]$")
    for k, x in enumerate((1.5, 9.5, 17.5)):
        box(axb, x, 34, 5.5, 4.5, "part", fs=5.6)
        arrow(axb, 12, 46, x + 2.75, 38.7)
    axb.text(12, 30.5, "combines several\nparts in one step", ha="center", va="top", fontsize=5.6, color=MUTED)
    arrow(axb, 24.5, 42, 30.5, 42, color=INK, lw=1.0)
    axb.text(27.5, 44, "rearrange", ha="center", fontsize=5.6, color=MUTED)
    # step 2: LinearCParty
    axb.text(48, 59, "LinearCParty", fontsize=6.4, fontweight="bold", ha="center", va="top")
    for k, x in enumerate((32, 43.5, 55)):
        box(axb, x, 40, 9.5, 4.5, "partial" if k < 2 else "region", fs=5.6, fc="#eef4fc" if k < 2 else "white")
        if k:
            arrow(axb, x - 1.8, 42.25, x - 0.3, 42.25)
    axb.text(48, 36.5, "adds one part at a time;\nevery state is ordered by its right end $j$",
             ha="center", va="top", fontsize=5.6, color=MUTED)
    arrow(axb, 67, 42, 73, 42, color=INK, lw=1.0)
    axb.text(70, 44, "prune", ha="center", fontsize=5.6, color=MUTED)
    # step 3: beam
    axb.text(87, 59, "Beam pruning", fontsize=6.4, fontweight="bold", ha="center", va="top")
    rnd = __import__("random").Random(4)
    for c in (76, 81, 86, 91, 96):
        nc = rnd.randint(4, 6)
        for r in range(nc):
            y = 33 + r * 3.0
            kept = r >= nc - 3
            axb.scatter([c], [y], s=9, color=NEW_COLOR if kept else "white",
                        edgecolor=NEW_COLOR if kept else G_COLOR, lw=0.6, zorder=3)
    arrow(axb, 74, 30.5, 99, 30.5, color=MUTED)
    axb.text(99, 28, "right end $j$ (5$'$ to 3$'$)", ha="right", va="top", fontsize=5.4, color=MUTED)
    axb.scatter([75.5], [24.5], s=9, color=NEW_COLOR, lw=0)
    axb.text(77, 24.5, "kept (top $b$)", fontsize=5.4, va="center", color=INK)
    axb.scatter([89], [24.5], s=9, color="white", edgecolor=G_COLOR, lw=0.6)
    axb.text(90.5, 24.5, "removed", fontsize=5.4, va="center", color=INK)
    # notes spanning panel B
    axb.text(2, 19.5, "Ranking: a candidate gets its inside weight plus a rough estimate of its outside weight,\n"
                    "taken from a preliminary run with a small beam. The top $b$ candidates are kept\n"
                    "for each state type and right end $j$. Local limits $h$ (length of a new hairpin) and\n"
                    "$d$ (depth of scaffold branching) bound the candidates; here $h = d = \\max(b, 30)$.\n"
                    "Without the beam and these limits, the same calculation is exact.",
             fontsize=5.7, color=INK, va="top", linespacing=1.45)

    # ---------------- Panel C: what one retained ensemble gives ----------------
    panel_label(fig, 0.005, 0.535, "C")
    fig.text(0.035, 0.533, "What is computed (bpRNA_PDB_345, 63 nt; exact; $G$ = nested pairs of the reference structure)",
             fontsize=7.5, fontweight="bold", va="top")
    # pair legend
    for k, (lab, col) in enumerate((("pair in $G$", G_COLOR), ("new pair, no crossing", NEW_COLOR),
                                    ("new pair crossing $G$", CROSS_COLOR))):
        x = 0.035 + k * 0.16
        fig.add_artist(matplotlib.lines.Line2D([x, x + 0.022], [0.502, 0.502], transform=T, color=col, lw=1.4))
        fig.text(x + 0.027, 0.502, lab, fontsize=5.8, color=INK, va="center")
    # bracket "one retained ensemble"
    fig.add_artist(Rectangle((0.028, 0.022), 0.966, 0.461, transform=T, fill=False, ec=NEW_COLOR, lw=0.7))
    fig.text(0.036, 0.477, "From one retained ensemble (the same probabilities are shared by all of these)",
             fontsize=6.4, fontweight="bold", color=NEW_COLOR, va="top")
    arrow_fig = FancyArrowPatch((0.86, 0.575), (0.86, 0.485), transform=T, arrowstyle="-|>", mutation_scale=7,
                                color=NEW_COLOR, lw=1.0)
    fig.add_artist(arrow_fig)
    lw_ = 0.455
    # C1 BPP arcs
    fig.text(0.036, 0.452, "Base-pairing probabilities (darker = higher)", fontsize=6.1, color=INK, va="top")
    ax = fig.add_axes([0.036, 0.35, lw_, 0.085])
    setup_arc_axis(ax, n, 34)
    gset = set(G)
    newpairs = [p for p in bpp if p not in gset]
    crossing = [p for p in newpairs if any(crosses(p, g) for g in G)]
    nested_new = [p for p in newpairs if p not in set(crossing)]
    draw_arcs(ax, G, n, G_COLOR, lw=0.9)
    draw_arcs(ax, nested_new, n, NEW_COLOR, alpha=bpp, lw=0.8)
    draw_arcs(ax, crossing, n, CROSS_COLOR, alpha=bpp, lw=0.9)
    # C2 contexts
    fig.text(0.036, 0.33, "Seven structural-context probabilities", fontsize=6.1, color=INK, va="top")
    ax = fig.add_axes([0.036, 0.215, lw_, 0.1])
    bottom = [0.0] * n
    for c in CONTEXTS:
        vals = ctx[c]
        ax.bar(range(1, n + 1), vals, bottom=bottom, width=0.86, color=CONTEXT_COLORS[c], lw=0, label=c)
        bottom = [x + y for x, y in zip(bottom, vals)]
    ax.set_xlim(0, n + 1); ax.set_ylim(0, 1.0); ax.set_yticks([0, 1]); ax.set_xticks([1, 20, 40, 63])
    for s_ in ("right", "top"):
        ax.spines[s_].set_visible(False)
    ax.legend(ncol=7, fontsize=5.3, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.2),
              handlelength=0.9, columnspacing=0.7, handletextpad=0.3)
    # C3 window accessibility
    fig.text(0.036, 0.162, "Window accessibility (all bases unpaired)", fontsize=6.1, color=INK, va="top")
    ax = fig.add_axes([0.036, 0.075, lw_, 0.07])
    for w, ls, col in ((1, "-", INK), (4, "--", MUTED), (8, ":", MUTED)):
        xs_ = sorted(acc[w])
        ax.plot(xs_, [acc[w][x] for x in xs_], color=col, lw=0.8, ls=ls, label="%d nt" % w)
    ax.set_xlim(0, n + 1); ax.set_ylim(0, 1.05); ax.set_yticks([0, 1]); ax.set_xticks([1, 20, 40, 63])
    ax.set_xlabel("start of window", fontsize=5.8, color=MUTED, labelpad=1)
    for s_ in ("right", "top"):
        ax.spines[s_].set_visible(False)
    ax.legend(fontsize=5.3, frameon=False, loc="lower right", bbox_to_anchor=(1.0, 1.12), ncol=3,
              handlelength=1.6, borderaxespad=0)
    # right column
    xr, xe = 0.53, 0.985
    mono = dict(fontsize=4.95, family="DejaVu Sans Mono", va="top")
    fig.text(xr, 0.452, "Boltzmann sampling: 1,000 structures", fontsize=6.1, color=INK, va="top")
    fig.text(xe, 0.452, "count", fontsize=5.4, color=MUTED, va="top", ha="right")
    fig.text(xr, 0.428, seq, color=MUTED, **mono)
    for k, (cnt, s_) in enumerate(samples[:3]):
        fig.text(xr, 0.41 - k * 0.017, s_, color=INK, **mono)
        fig.text(xe, 0.41 - k * 0.017, "%d" % cnt, fontsize=5.4, color=MUTED, va="top", ha="right")
    fig.text(xr, 0.33, "Maximum expected accuracy structure ($\\gamma = 1$)", fontsize=6.1, color=INK, va="top")
    fig.text(xr, 0.308, mea, color=INK, **mono)
    # dashed box: separate maximization over the same decomposition
    fig.add_artist(Rectangle((xr - 0.007, 0.062), xe - xr + 0.012, 0.2, transform=T, fill=False,
                             ec=MUTED, lw=0.6, ls=(0, (3, 2))))
    fig.text(xr, 0.255, "Same decomposition, separate maximization", fontsize=6.1, color=MUTED,
             fontweight="bold", va="top")
    fig.text(xr, 0.233, "Minimum-free-energy and suboptimal structures", fontsize=6.1, color=INK, va="top")
    fig.text(xe, 0.233, "kcal/mol", fontsize=5.4, color=MUTED, va="top", ha="right")
    for k, (s_, e) in enumerate(subopt[:4]):
        fig.text(xr, 0.21 - k * 0.017, s_, color=INK, **mono)
        fig.text(xe, 0.21 - k * 0.017, "%.1f" % e, fontsize=5.4, color=MUTED, va="top", ha="right")
    fig.text(xr, 0.125, "The minimum-free-energy search uses the same pieces as the ensemble,\n"
                        "but takes the best one instead of the sum. With a finite beam it keeps\n"
                        "its own candidates.", fontsize=5.6, color=MUTED, va="top", linespacing=1.35)

    for ext in ("pdf", "png"):
        fig.savefig(args.out / f"figure1_draft.{ext}", dpi=300)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
