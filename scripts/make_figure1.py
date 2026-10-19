#!/usr/bin/env python3
"""Figure 1: scaffold, decomposition and outputs from one conditional ensemble.

Render the saved exact bpRNA_PDB_345 example at 178 x 145 mm.
No folding program is run and no experimental values are invented.
Example, from the experiments repository::

    MPLCONFIGDIR=/tmp/figure1-mpl python3 scripts/make_figure1.py \
        --data results/figures/figure1_20261016/data \
        --out results/figures/figure1_20261017

The output directory receives figure1.pdf, figure1.png (300 dpi), caption.md,
README.md and validation.json. The old figure and saved inputs are unchanged.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, FancyArrowPatch, FancyBboxPatch, Rectangle
from matplotlib.text import Text
from PIL import Image
from pypdf import PdfReader

MM = 1 / 25.4
WIDTH, HEIGHT = 178, 145
INK, MUTED, GRID = "#1f1f1e", "#6b6b66", "#d9d8d3"
G_COLOR, NEW_COLOR, CROSS_COLOR = "#8a8a85", "#2a78d6", "#eb6834"
CONTEXTS = ["Stem", "Hairpin", "Bulge", "Internal", "Multiloop", "Exterior", "PseudoLoop"]
CONTEXT_COLORS = dict(zip(CONTEXTS, [NEW_COLOR, "#1baf7a", "#eda100", "#e87ba4",
                                     "#6250d6", "#b9b8b1", CROSS_COLOR]))
PAIR_CUTOFF = 0.01
BEAM_DEMO = 3  # A drawing parameter only; panel C is exact.
INPUT_FILES = ["rec.json", "in.txt", "lc.out", "bpp.tsv", "capr.out", "racc.out", "samp.out", "mea.out"]

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 7, "text.color": INK,
    "axes.linewidth": 0.6, "axes.edgecolor": MUTED,
    "xtick.color": INK, "ytick.color": INK,
    "xtick.labelsize": 7, "ytick.labelsize": 7,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "pdf.fonttype": 42, "hatch.linewidth": 0.35,
})


def pairs_from_db(db):
    """Read 1-based pair coordinates, irrespective of the bracket alphabet."""
    stacks = {"(": [], "[": [], "{": [], "<": []}
    close = {")": "(", "]": "[", "}": "{", ">": "<"}
    pairs = []
    for pos, char in enumerate(db, 1):
        if char in stacks:
            stacks[char].append(pos)
        elif char in close:
            pairs.append((stacks[close[char]].pop(), pos))
        elif char != ".":
            raise ValueError(f"Unexpected structure character {char!r}")
    if any(stacks.values()):
        raise ValueError("Unbalanced structure")
    return sorted(pairs)


def crosses(a, b):
    (i, j), (k, ell) = a, b
    return i < k < j < ell or k < i < ell < j


def read_outputs(data):
    rec = json.loads((data / "rec.json").read_text())
    input_lines = (data / "in.txt").read_text().splitlines()
    seq, scaffold = input_lines[-2:]
    n = len(seq)
    G, Gp = pairs_from_db(scaffold), sorted(map(tuple, rec["pk"]))
    assert seq == rec["seq"] and len(scaffold) == n
    assert set(G) == set(map(tuple, rec["nested"]))
    assert set(G) <= set(pairs_from_db(rec["ref"]))
    assert all(seq[i - 1] + seq[j - 1] in {"AU", "UA", "CG", "GC", "GU", "UG"} for i, j in G)
    assert not ({x for p in G for x in p} & {x for p in Gp for x in p})
    assert not any(crosses(a, b) for k, a in enumerate(G) for b in G[k + 1:])
    assert not any(crosses(a, b) for k, a in enumerate(Gp) for b in Gp[k + 1:])
    bpp = {}
    for line in (data / "bpp.tsv").read_text().splitlines():
        fields = line.split()
        if len(fields) == 3 and fields[0].isdigit():
            bpp[int(fields[0]), int(fields[1])] = float(fields[2])
    ctx = {}
    for line in (data / "capr.out").read_text().splitlines():
        fields = line.split()
        if fields and fields[0] in CONTEXTS:
            ctx[fields[0]] = list(map(float, fields[1:]))
    acc = {}
    access_text = (data / "racc.out").read_text()
    for line in access_text.splitlines():
        if line.startswith(("#", ">")) or not line.strip():
            continue
        pos, rest = line.split("\t")
        for item in rest.strip(";\n").split(";"):
            width, probability = item.split(",")
            acc.setdefault(int(width), {})[int(pos) + 1] = float(probability)
    samples = []
    for line in (data / "samp.out").read_text().splitlines():
        fields = line.split("\t")
        if len(fields) == 2 and fields[0].isdigit():
            samples.append((int(fields[0]), fields[1]))
    samples.sort(key=lambda row: -row[0])  # Stable for ties in saved output.
    mea_text = (data / "mea.out").read_text()
    mea_line = next(line for line in mea_text.splitlines() if "EA=" in line)
    mea = mea_line.split()[0]
    gamma = float(re.search(r"gamma=([^)]*)", mea_line).group(1))
    mfe_line = next(line for line in (data / "lc.out").read_text().splitlines()
                    if re.fullmatch(r"[.()\[\]{}<>]+\s+\(-?[\d.]+\)", line))
    mfe = mfe_line.split()[0]
    assert "# beam=exact" in access_text
    assert f"# scaffold={scaffold}" in access_text
    assert set(ctx) == set(CONTEXTS) and all(len(row) == n for row in ctx.values())
    context_error = max(abs(sum(ctx[c][i] for c in CONTEXTS) - 1) for i in range(n))
    assert context_error < 1e-10
    assert all(abs(bpp[p] - 1) < 1e-10 for p in G)
    assert all(len(s) == n and set(G) <= set(pairs_from_db(s)) for _, s in samples)
    assert set(G) <= set(pairs_from_db(mea)) and set(G) <= set(pairs_from_db(mfe))
    assert all(0 <= p <= 1 for p in bpp.values())
    assert all(len(row) == n - width + 1 for width, row in acc.items())
    # Select the most common draw and the most common draw with a new pair
    # that does not cross G. This makes both added-pair classes visible.
    def has_noncrossing(db):
        return any(pair not in G and not any(crosses(pair, g) for g in G)
                   for pair in pairs_from_db(db))
    second = next((row for row in samples[1:] if has_noncrossing(row[1])), samples[1])
    shown = [samples[0], second]
    # A's two bands must have identical crossing partners within each band.
    assert all({p for p in Gp if crosses(g, p)} == set(Gp) for g in G)
    assert all({g for g in G if crosses(p, g)} == set(G) for p in Gp)
    spans = [(min(i for i, j in group), max(j for i, j in group)) for group in (G, Gp)]
    probe = (max(i for i, j in spans) + min(j for i, j in spans)) / 2
    overlap = sum(i < probe < j for i, j in spans)
    assert overlap == 2
    # Choose the longest region where PseudoLoop is the most likely context.
    runs, start = [], None
    for i in range(n + 1):
        dominant = i < n and ctx["PseudoLoop"][i] == max(ctx[c][i] for c in CONTEXTS)
        if dominant and start is None:
            start = i + 1
        if not dominant and start is not None:
            runs.append((start, i)); start = None
    pseudo_region = max(runs, key=lambda region: region[1] - region[0])
    return dict(rec=rec, seq=seq, n=n, G=G, Gp=Gp, bpp=bpp, ctx=ctx, acc=acc,
                samples=samples, shown=shown, total=sum(count for count, db in samples),
                mea=mea, gamma=gamma, mfe=mfe, spans=spans, probe=probe,
                overlap=overlap, pseudo_region=pseudo_region, context_error=context_error)


def make_figure(d):
    fig = plt.figure(figsize=(WIDTH * MM, HEIGHT * MM))
    canvas = fig.add_axes([0, 0, 1, 1], zorder=0)
    canvas.set(xlim=(0, WIDTH), ylim=(0, HEIGHT))
    canvas.axis("off")

    def text(x, y, label, size=7, **kwargs):
        kwargs.setdefault("va", "top")
        return canvas.text(x, y, label, fontsize=size, **kwargs)

    def arrow(x0, y0, x1, y1, color=MUTED, **kwargs):
        canvas.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>",
                                        mutation_scale=6, lw=0.8, color=color,
                                        shrinkA=0, shrinkB=0, **kwargs))

    def box(x, y, w, h, label=None, fill="white", edge=MUTED, **kwargs):
        canvas.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=0.7",
                                       fc=fill, ec=edge, lw=0.7, **kwargs))
        if label:
            canvas.text(x + w / 2, y + h / 2, label, fontsize=7, va="center", ha="center")

    def panel(x, y, letter, title):
        text(x, y, letter, 9, weight="bold")
        text(x + 6, y - 0.2, title, 8, weight="bold")

    def arc_axis(x, y, w, h):
        ax = fig.add_axes([x / WIDTH, y / HEIGHT, w / WIDTH, h / HEIGHT])
        ax.set(xlim=(0, d["n"] + 1), ylim=(-0.7, 30))
        ax.axhline(0, color=MUTED, lw=0.55)
        ax.axis("off")
        return ax

    def arcs(ax, pairs, color, probabilities=None, crossing=False, width=0.9):
        for i, j in pairs:
            probability = probabilities[i, j] if probabilities is not None else 1
            if probability <= PAIR_CUTOFF:
                continue
            ax.add_patch(Arc(((i + j) / 2, 0), j - i, (j - i) * 0.9,
                             theta1=0, theta2=180, color=color,
                             alpha=probability, lw=width,
                             ls=(0, (4, 2)) if crossing else
                                (0, (1, 1)) if color == NEW_COLOR else "-"))

    def completed(ax, pairs):
        new_pairs = [p for p in pairs if p not in d["G"]]
        crossing = [p for p in new_pairs if any(crosses(p, g) for g in d["G"])]
        arcs(ax, d["G"], G_COLOR)
        arcs(ax, [p for p in new_pairs if p not in crossing], NEW_COLOR)
        arcs(ax, crossing, CROSS_COLOR, crossing=True)

    def probability_axis(x, y, w, h, ticks=False):
        ax = fig.add_axes([x / WIDTH, y / HEIGHT, w / WIDTH, h / HEIGHT])
        ax.set(xlim=(0, d["n"] + 1), ylim=(0, 1.04), yticks=[0, 1])
        ax.spines[["right", "top"]].set_visible(False)
        if ticks:
            ax.set_xticks([1] + list(range(20, d["n"] - 8, 20)) + [d["n"]])
        else:
            ax.set_xticks([])
        ax.tick_params(pad=1.3)
        return ax

    # A: exactly the three components requested by the author.
    panel(2, 142, "A", "Keep G; add G′")
    for label, y, pairs, color, dashed in (
        ("Scaffold G", 134, d["G"], G_COLOR, False),
        ("Extension G′", 122, d["Gp"], CROSS_COLOR, True),
    ):
        text(7, y, label, 7.5)
        arcs(arc_axis(7, y - 10, 50, 7.5), pairs, color, crossing=dashed)
    text(7, 110, "R = G ∪ G′", 7.5)
    ax = arc_axis(7, 97, 50, 10)
    completed(ax, d["G"] + d["Gp"])
    # The brackets show the spans of the two bands, not individual pair arcs.
    map_x = lambda pos: 7 + pos * 50 / (d["n"] + 1)
    for (i, j), y, color, style in zip(d["spans"], (94, 91), (G_COLOR, CROSS_COLOR), ("-", "--")):
        canvas.plot([map_x(i), map_x(i), map_x(j), map_x(j)], [y + 0.8, y, y, y + 0.8],
                    color=color, lw=1, ls=style)
    probe_x = map_x(d["probe"])
    canvas.plot([probe_x, probe_x], [90, 106], color=INK, lw=0.7, ls=":")
    text(7, 87.5, "Two bands span this position", 7)
    arrow(probe_x, 87.8, probe_x, 90, color=INK)

    # B: decompose -> smaller combinations -> retained candidates.
    panel(63, 142, "B", "Rearrange CParty; retain candidates")
    text(79, 134, "CParty", 7.5, ha="center", weight="bold")
    text(121, 134, "Smaller steps", 7.5, ha="center", weight="bold")
    text(162, 134, "Beam pruning", 7.5, ha="center", weight="bold")
    box(69, 122, 20, 5, "region")
    for x in (65, 75, 85):
        box(x, 112, 8, 5, "part")
        arrow(79, 121.7, x + 4, 117.5)
    text(79, 109, "Combine several\nparts in one step", ha="center", linespacing=1.25)
    arrow(95, 121.5, 103, 121.5, color=INK)
    # A binary tree illustrates combining two pieces before the third.
    box(110, 123, 21, 5, "region")
    box(105, 114, 17, 5, "partial", fill="#eef4fc")
    box(126, 114, 8, 5, "part")
    arrow(113.5, 119.4, 118, 122.6)
    arrow(130, 119.4, 123, 122.6)
    for x in (106.5, 118.5):
        canvas.add_patch(Rectangle((x - 1.2, 109), 2.4, 2.4, facecolor="white", edgecolor=MUTED, lw=0.6))
        arrow(x, 111.6, 113.5, 113.6)
    text(120, 106.5, "Combine two, then add one", ha="center")
    arrow(136, 121.5, 143, 121.5, color=INK)
    # Top b filled circles versus discarded crosses, ordered by rank.
    for column in (155, 163, 171):
        for row in range(5):
            retained = row >= 5 - BEAM_DEMO
            canvas.plot(column, 112 + row * 3.2,
                        marker="o" if retained else "x", ms=3.3 if retained else 3.0,
                        color=NEW_COLOR if retained else G_COLOR, mew=0.7, ls="none")
    arrow(149, 111, 149, 127)
    text(146.7, 119, "Higher rank", rotation=90, ha="center", va="center")
    arrow(153, 109, 174, 109)
    text(163, 107.5, "Right endpoint", ha="center")
    # MFE has a separate arrow from the rearranged decomposition, outside C.
    canvas.plot([139.5, 139.5], [121.5, 102.8], color=MUTED, lw=0.8)
    arrow(139.5, 104.8, 139.5, 102.5)
    box(102, 87, 54, 15, edge=MUTED, ls=(0, (3, 2)))
    text(104, 100.5, "Separate MFE search", 7.3, weight="bold")
    text(104, 95.5, "Own beam\ncandidates", linespacing=1.2)
    completed(arc_axis(132, 89, 22, 7), pairs_from_db(d["mfe"]))
    # Ensemble branch descends from beam pruning directly into the blue frame.
    arrow(169, 103.5, 169, 77.3, color=NEW_COLOR)
    text(171, 94, "Retain", color=NEW_COLOR, rotation=90, ha="center", va="center")

    # C: only the five outputs derived from the same retained ensemble.
    panel(2, 82.5, "C", "Outputs from one retained ensemble")
    canvas.add_patch(Rectangle((5, 1.3), 170, 75.7, facecolor="none", edgecolor=NEW_COLOR, lw=0.8))
    text(12, 74, "Pair probabilities", 7.5, weight="bold")
    bpp_ax = arc_axis(12, 59.5, 85, 11)
    new_pairs = [p for p in d["bpp"] if p not in d["G"]]
    crossing = [p for p in new_pairs if any(crosses(p, g) for g in d["G"])]
    arcs(bpp_ax, d["G"], G_COLOR, probabilities=d["bpp"])
    arcs(bpp_ax, [p for p in new_pairs if p not in crossing], NEW_COLOR, probabilities=d["bpp"])
    arcs(bpp_ax, crossing, CROSS_COLOR, probabilities=d["bpp"], crossing=True)
    text(12, 56.8, "Seven context probabilities", 7.5, weight="bold")
    ax = probability_axis(12, 35.5, 85, 16.5)
    bottom = [0.0] * d["n"]
    # Continuous steps make an unpaired region legible without per-base gaps.
    x = [pos + 0.5 for pos in range(d["n"] + 1)]
    for context in CONTEXTS:
        top = [a + b for a, b in zip(bottom, d["ctx"][context])]
        ax.fill_between(x, bottom + [bottom[-1]], top + [top[-1]], step="post",
                        facecolor=CONTEXT_COLORS[context], edgecolor=MUTED,
                        hatch="/" if context == "Internal" else "", lw=0)
        bottom = top
    begin, end = d["pseudo_region"]
    center = (begin + end) / 2
    ax.text(center, 0.52, "PseudoLoop", fontsize=7, ha="center", va="center", color=INK)
    ax.plot([begin - 0.5, begin - 0.5, end + 0.5, end + 0.5], [1.015, 1.055, 1.055, 1.015],
            color=INK, lw=0.7, clip_on=False)
    # Named swatches, two rows; PseudoLoop also has a direct in-region label.
    for y, labels in ((32.4, CONTEXTS[:4]), (28.4, CONTEXTS[4:])):
        xpos = 12
        for label in labels:
            canvas.add_patch(Rectangle((xpos, y - 2.0), 2, 2,
                                       facecolor=CONTEXT_COLORS[label], edgecolor=MUTED,
                                       hatch="/" if label == "Internal" else "", lw=0))
            if label == "Internal":
                # One sparse slash remains visible in this 2-mm key swatch.
                canvas.plot([xpos + .2, xpos + 1.8], [y - 1.8, y - .2],
                            color=MUTED, lw=.35)
            text(xpos + 2.8, y, label, 7)
            xpos += {"Stem": 15, "Hairpin": 18, "Bulge": 16, "Internal": 20,
                     "Multiloop": 23, "Exterior": 21, "PseudoLoop": 23}[label]
    text(12, 23.6, "Window accessibility", 7.5, weight="bold")
    ax = probability_axis(12, 9, 85, 10.5, ticks=True)
    for width, style in zip(sorted(d["acc"]), ("-", "--", ":")):
        xs = sorted(d["acc"][width])
        ax.plot(xs, [d["acc"][width][pos] for pos in xs], color=INK, lw=0.9, ls=style)
    # Explicit line keys preserve the identities in greyscale.
    for k, (width, style) in enumerate(zip(sorted(d["acc"]), ("-", "--", ":"))):
        xx = 60 + 12 * k
        canvas.plot([xx, xx + 3.5], [22.2, 22.2], color=INK, lw=0.9, ls=style)
        text(xx + 4, 23.5, f"{width} nt", 7)
    text(54.5, 4.8, "Position / start of window", 7, ha="center")

    text(111, 74, "Boltzmann samples", 7.5, weight="bold")
    for k, ((count, db), base) in enumerate(zip(d["shown"], (58.5, 43))):
        text(111, base + 12, f"Most frequent: {count:,} draws" if k == 0 else
             f"With a new non-crossing pair: {count:,} draws", 7)
        sample_ax = arc_axis(111, base, 60, 9.7)
        sample_pairs = pairs_from_db(db)
        completed(sample_ax, sample_pairs)
        if k == 1:
            i, j = next(p for p in sample_pairs if p not in d["G"] and
                        not any(crosses(p, g) for g in d["G"]))
            sample_ax.annotate("", xy=((i + j) / 2, (j - i) * .45),
                               xytext=((i + j) / 2 + 5, (j - i) * .45 + 8),
                               arrowprops=dict(arrowstyle="->", color=NEW_COLOR, lw=.8))
    other = d["total"] - sum(count for count, db in d["shown"])
    text(111, 39.8, f"Other draws: {other:,} of {d['total']:,}", 7)
    text(111, 32.5, "MEA structure", 7.5, weight="bold")
    text(171, 32.3, f"γ = {d['gamma']:g}", 7, ha="right")
    completed(arc_axis(111, 19, 60, 9.7), pairs_from_db(d["mea"]))
    # A compact pair key uses line styles as well as the common colors.
    for yy, label, color, style in ((15.5, "Scaffold G", G_COLOR, "-"),
                                    (11.2, "New; no crossing", NEW_COLOR, (0, (1, 1))),
                                    (6.9, "New; crosses G", CROSS_COLOR, (0, (4, 2)))):
        canvas.plot([111, 117], [yy - 1.1, yy - 1.1], color=color, lw=1, ls=style)
        text(119, yy, label, 7)
    return fig


def validate_figure(fig, pdf, png):
    """Check actual PDF text sizes, figure dimensions and rendered text extents."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    text_artists = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    outside = []
    for artist in text_artists:
        bounds = artist.get_window_extent(renderer)
        # Tick labels from invisible axes can be stale; all used axes have
        # explicit limits and only active artists are included by their axes.
        if artist.axes is not None and not artist.axes.axison and artist not in artist.axes.texts:
            continue
        if not fig.bbox.contains(bounds.x0, bounds.y0) or not fig.bbox.contains(bounds.x1, bounds.y1):
            outside.append(artist.get_text())
    assert not outside, f"Text outside page: {outside}"
    command = ["gs", "-q", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0",
               "-sOutputFile=-", str(pdf)]
    extracted = subprocess.run(command, check=True, text=True, capture_output=True).stdout
    sizes = [float(size) for size in re.findall(r' size="([\d.]+)"', extracted)]
    assert sizes and min(sizes) >= 6, f"Embedded text below 6 pt: {min(sizes)}"
    page = PdfReader(str(pdf)).pages[0]
    page_mm = [float(page.mediabox.width) * 25.4 / 72, float(page.mediabox.height) * 25.4 / 72]
    assert all(abs(actual - target) < 0.01 for actual, target in zip(page_mm, (WIDTH, HEIGHT)))
    with Image.open(png) as rendered:
        pixels, dpi = list(rendered.size), list(rendered.info["dpi"])
    assert all(abs(value - 300) < 0.1 for value in dpi)
    return dict(pdf_page_mm=page_mm, png_pixels=pixels, png_dpi=dpi,
                pdf_text_sizes_pt=sorted(set(sizes)), minimum_embedded_text_pt=min(sizes),
                minimum_matplotlib_text_pt=min(t.get_fontsize() for t in text_artists),
                text_outside_page=outside, gs_command=command,
                visual_inspection="Required separately: open the generated PNG.")


def write_notes(data, out, d, report):
    prefix = Path(__file__).resolve().parent.parent
    source = data.resolve().relative_to(prefix).as_posix() if data.resolve().is_relative_to(prefix) else str(data.resolve())
    shown = sum(count for count, db in d["shown"])
    other = d["total"] - shown
    windows = ", ".join(str(w) for w in sorted(d["acc"]))
    lo, hi = d["pseudo_region"]
    caption = f"""**Figure 1. Conditional ensemble analysis with LinearCParty.**

<!-- 日本語：Aのscaffold Gとextension G′は対合塩基を共有せず、それぞれは交差対を持たない。
和集合がRを作る。
括弧はband、すなわち外部の同じ対群と交差する入れ子対の極大群の範囲を示す。
点線は二つのbandの範囲を通る。
density-2制限では、交差で連結したband群のうち同じ位置を覆うbandは最大二つである。
残りのloopと境界の制限はMethodsに示す。 -->
**(A)** Scaffold G and extension G′ share no paired bases, and each has no crossing pairs.
Their union forms R.
The brackets mark bands: maximal groups of nested pairs that cross the same outside pairs.
The dotted line passes through both band spans.
The density-2 restriction allows at most two bands from any connected group of crossing bands to span one position.
Methods give the remaining loop and boundary restrictions.

<!-- 日本語：BはCPartyの分解を小さな組合せへ並べ替え、右端の順に候補を処理する。
丸は保持候補、×は破棄候補を示す。
高い点ほど順位が高い（図の例では右端ごとに{BEAM_DEMO}候補を保持）。
順位はinside weightの対数に予備計算の周囲の重みの対数の推定値を足して決める。
inside weightは部分構造のBoltzmann重みの和である。 -->
**(B)** LinearCParty rearranges CParty's decomposition into smaller combinations and processes candidates by right endpoint.
Circles mark retained candidates and crosses mark discarded candidates.
Higher points have higher ranks ({BEAM_DEMO} retained per endpoint in this sketch).
Ranking adds the log inside weight to an estimate of surrounding log weight from a preliminary run.
Inside weight sums the Boltzmann weights of partial structures.

<!-- 日本語：Cの青枠は一つの保持ensembleから得る対確率、七つのcontext確率、窓accessibility、Boltzmann標本、最大期待正確度（MEA）復元の五出力を含む。
AとCは{d['rec']['id']}（{d['n']} nt）を使い、canonicalな入れ子の参照{len(d['G'])}対をGとする。
Cはモデルのexact予測である。
灰色実線はG、青色点線は新しい非交差対、橙色破線はGと交差する新しい対を示す。
対の弧は確率{PAIR_CUTOFF:g}より大きいものを示し、高確率ほど濃くする。 -->
**(C)** The blue frame contains five outputs from one retained ensemble: pair probabilities, seven context probabilities, window accessibility, Boltzmann samples, and maximum expected accuracy (MEA) decoding.
A and C use {d['rec']['id']} ({d['n']} nt), with {len(d['G'])} canonical nested reference pairs as G.
C shows exact model predictions.
Solid grey, dotted blue, and dashed orange arcs mark G, new non-crossing pairs, and new pairs crossing G, respectively.
Pair arcs show probabilities greater than {PAIR_CUTOFF:g}, with darker arcs for higher probabilities.

<!-- 日本語：context確率の和は1である。
Internalにはこの図とTercの例で共通の疎な斜線を付ける。
括弧は塩基{lo}–{hi}のPseudoLoopを示す。
PseudoLoopは交差bandをつなぐ、またはbandをまたぐ分岐内の非対合領域を表す。
accessibilityは窓の全塩基が同時に非対合となる確率である（{windows} nt、実線・破線・点線）。
座標は1始まりの塩基位置または窓の開始位置である。 -->
Context probabilities sum to one.
Sparse hatching identifies Internal in this figure and the Terc example.
The bracket labels PseudoLoop at bases {lo}–{hi}.
PseudoLoop labels unpaired regions connecting crossing bands or in branches spanning a band.
Accessibility is the probability that all window bases are unpaired ({windows} nt; solid, dashed, dotted lines).
Coordinates give 1-based positions or window starts.

<!-- 日本語：標本図は最多の構造（{d['shown'][0][0]}回）と、新しい非交差対を含む構造の中で最多のもの（{d['shown'][1][0]}回）を示す。
全{d['total']}回のうち残り{other}回の標本は示さない。
MEAの対の重みはγ={d['gamma']:g}とする。
別の最小自由エネルギー（MFE）探索は同じ分解を使い、独自のbeam候補を持つ。
この模式図の構造もexactの結果を示す。 -->
The diagrams show the most frequent sampled structure ({d['shown'][0][0]:,} draws) and the most frequent structure containing a new non-crossing pair ({d['shown'][1][0]:,} draws).
The remaining {other:,} of {d['total']:,} draws are not shown.
MEA uses pair weight γ = {d['gamma']:g}.
The separate minimum-free-energy (MFE) search uses the same decomposition with its own beam candidates.
Its diagram also shows an exact result.
"""
    (out / "caption.md").write_text(caption)
    hashes = {name: hashlib.sha256((data / name).read_bytes()).hexdigest() for name in INPUT_FILES}
    report.update(input_directory=source, input_sha256=hashes,
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  sequence_length=d["n"], scaffold_pair_count=len(d["G"]), extension_pair_count=len(d["Gp"]),
                  displayed_bpp_pair_count=sum(p > PAIR_CUTOFF for p in d["bpp"].values()),
                  bpp_display_cutoff_exclusive=PAIR_CUTOFF,
                  sample_draws_total=d["total"], sample_counts_shown=[c for c, db in d["shown"]],
                  sample_draws_other=other, sample_selection="Most frequent; most frequent with a new non-crossing pair",
                  accessibility_window_lengths=sorted(d["acc"]), mea_gamma=d["gamma"],
                  pseudoloop_label_region=[lo, hi], band_spans=d["spans"],
                  band_overlap_probe=d["probe"], bands_spanning_probe=d["overlap"],
                  max_context_sum_error=d["context_error"],
                  pair_linestyles=dict(scaffold="solid", new_no_crossing="dotted",
                                       new_crossing_G="dashed"),
                  context_hatches={c: "/" if c == "Internal" else "" for c in CONTEXTS},
                  illustrative_beam_width=BEAM_DEMO)
    (out / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    readme = f"""# Figure 1 revision — 2026-10-17

`figure1.pdf` and `figure1.png` are {WIDTH} × {HEIGHT} mm at final main-text width.
The PNG is 300 dpi.
The prior 183 × 150 mm figure becomes 145.90 mm tall at 178 mm width, so this revision does not increase its published height.
`caption.md` is a proposed caption with one sentence per line.
`validation.json` records input SHA-256 hashes, values derived from saved outputs, the script hash, and PDF/PNG checks.

## Reproduce

From the experiments repository, use a node-local temporary directory for Matplotlib:

```sh
python3 - <<'PYTHON'
import os, subprocess, tempfile
with tempfile.TemporaryDirectory(prefix="figure1-", dir="/tmp") as scratch:
    subprocess.run([
        "python3", "scripts/make_figure1.py",
        "--data", "{source}",
        "--out", "results/figures/figure1_20261017",
    ], env=dict(os.environ, MPLCONFIGDIR=scratch, PYTHONDONTWRITEBYTECODE="1"), check=True)
PYTHON
```


Dependencies: Python, Matplotlib, Pillow, pypdf and Ghostscript (`gs`).
No folding CLI, network access or new experimental run is needed.
The script writes only five final files to the requested output directory.

## Inputs and implementation version

All scientific values come from `{source}/`.
The source README records the implementation version and saved-output conditions:

{(data / 'README.md').read_text().strip()}

This rendering makes no cross-version equality assumption.
The source record originally came from `results/local/pk-recovery-bprna-20261016/pk_sel.json`.
No saved input or old figure has been modified.

- `in.txt` supplies the RNA sequence and actual scaffold; `rec.json` supplies the bpRNA identifier and reference extension in A.
- `lc.out` supplies the exact MFE structure; `bpp.tsv` supplies pair probabilities; `capr.out` supplies seven context probabilities.
- `racc.out` supplies the exact-mode/scaffold metadata and probabilities for {windows} nt windows; starts are converted from 0-based output to 1-based coordinates.
- `samp.out` supplies all unique structures and counts from `LinearPKSampling --exact -k 1000 --seed 20261016 --unique`.
- `mea.out` supplies the exact MEA structure and γ; the source README records `LinearPKMEA --exact --gamma 1`.
- Original ensemble commands: `LinearCParty --exact --bpp-cutoff 0 --bpp bpp.tsv`; `LinearPKCapR --exact`; `LinearPKRaccess --exact --access-len 1,4,8 --probability`.

The input has {d['n']} bases and {len(d['G'])} scaffold pairs.
The caption says canonical nested reference pairs because the original reference also has a noncanonical nested pair that is absent from the saved scaffold.
A illustrates the saved reference extension; C contains model predictions, not experimentally established structures.
B is a schematic, not a measured beam run: its example width is {BEAM_DEMO}.
The two spans in A are computed from the two pair sets, and the dotted line intersects {d['overlap']} bands in their connected crossing group.
The drawing is an example of density-2, not a complete graphical definition of all allowed structures.

## Derived display values

The pair plot uses the strict display cutoff p > {PAIR_CUTOFF:g} ({report['displayed_bpp_pair_count']} visible pairs).
The two samples have counts {d['shown'][0][0]} and {d['shown'][1][0]}, leaving {other} of {d['total']} draws unshown.
The second sample is selected to contain a new pair that does not cross G, rather than being the second most frequent structure overall.
PseudoLoop's direct label marks bases {lo}–{hi}, the longest contiguous region where it is the largest context probability.
Accessibility uses saved windows of {windows} nt; no new loop-specific accessibility is claimed.

## Review changes

- F1-1 / G1: remove sequences, dot-bracket strings, suboptimal lists and the in-figure methods paragraph; use 7–8 pt text and 9 pt panel letters at 178 mm width.
- F1-2: keep only the five ensemble-derived outputs in the blue frame; draw a separate MFE branch from the rearranged decomposition.
- F1-3 / L16 / L17: mark both band spans and their overlap in A; put the definition and scope of density-2 in the caption.
- F1-4: show rank increasing upward; define the sum of log inside weight and estimated surrounding log weight in the caption.
- F1-5: compute shown/other sample counts and define the strict probability cutoff, window lengths, line styles and coordinates in the caption.
- G2 / L28: label PseudoLoop directly in its region; use dashed crossing arcs and discarded-cross symbols as cues beyond color.
- G3: give A, B and C independent bold 9 pt labels with consistent offsets.
- L21: remove suboptimal outputs from the main figure.

## Validation and handoff

Saved-data consistency checks passed: actual scaffold equals the recorded nested subset, all its pairs have BPP 1, every sampled/decoded structure retains G, window row counts match the RNA length, and context probabilities sum to one (maximum error {d['context_error']:.3g}).
Ghostscript extracted embedded PDF sizes {report['pdf_text_sizes_pt']} pt; minimum {report['minimum_embedded_text_pt']:g} pt.
PDF media box: {report['pdf_page_mm'][0]:.3f} × {report['pdf_page_mm'][1]:.3f} mm.
PNG: {report['png_pixels'][0]} × {report['png_pixels'][1]} pixels, {report['png_dpi'][0]:.3f} dpi.
The rendered-text bounds check found no text outside the page.
A separate visual inspection of the PNG is required after generation; the final handoff log records that inspection.

The manuscript and its Japanese sources have not been edited.
The manuscript owner must insert this PDF and adapt the caption, then check pagination and the final typeset text sizes.
The caption must replace the old statement that the sample diagrams show the top three draws.
A third-band exclusion example is left out to keep A limited to G, G′ and R as the author requested; the caption gives the general condition.
No scientific claim or reference is added and no author decision is required to reproduce the figure.
The author may choose to redraw the same layout in PowerPoint.
Because this task is restricted to its assigned script/output/log, the manuscript owner should copy the substantive changes and validation into `docs/progress.md`.
"""
    (out / "README.md").write_text(readme)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.data.resolve().parent == args.out.resolve():
        parser.error("Choose a new output directory; do not overwrite the input figure.")
    data = read_outputs(args.data)
    args.out.mkdir(parents=True, exist_ok=True)
    fig = make_figure(data)
    pdf, png = args.out / "figure1.pdf", args.out / "figure1.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=300)
    report = validate_figure(fig, pdf, png)
    write_notes(args.data, args.out, data, report)
    plt.close(fig)
    print(json.dumps({"output": str(args.out), "minimum_embedded_text_pt": report["minimum_embedded_text_pt"],
                      "figure_mm": report["pdf_page_mm"], "samples_shown": [c for c, db in data["shown"]]}, indent=2))


if __name__ == "__main__":
    main()
