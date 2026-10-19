#!/usr/bin/env python3
"""Make Figure 4 from saved rat Terc outputs; never run an RNA executable.

Defaults are relative to this experiment repository, so the paper checkout
and its manuscript are not required. Coordinates in figures/tables are
1-based, inclusive; accessibility coordinates denote the window's first base.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "figure4-terc-mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Arc, Patch, Rectangle
from matplotlib.text import Text
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = ROOT / "results/local/figure4-terc-20261017/febedd98"
DEFAULT_OUT = ROOT / "results/figures/figure4_terc_20261017b"
MM = 1 / 25.4
WIDTH, HEIGHT, DPI = 178, 100, 300
REGION, WINDOW, CUTOFF = (40, 160), 4, 0.01
INK, MUTED, GRID = "#1f1f1e", "#6b6b66", "#d9d8d3"
G_COLOR, NEW_COLOR, CROSS_COLOR = "#8a8a85", "#2a78d6", "#eb6834"
CONTEXTS = ("Stem", "Hairpin", "Bulge", "Internal", "Multiloop", "Exterior", "PseudoLoop")
COLORS = dict(zip(CONTEXTS, ("#2a78d6", "#1baf7a", "#eda100", "#e87ba4",
                            "#6250d6", "#b9b8b1", "#eb6834")))
HATCHES = {c: "/" if c == "Internal" else "" for c in CONTEXTS}
LOOPS = dict(zip("EHBIMP", ("Exterior", "Hairpin", "Bulge", "Internal", "Multiloop", "PseudoLoop")))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative(path):
    return os.path.relpath(path.resolve(), ROOT)


def pairs_of(db):
    stacks = {c: [] for c in "([{<"}
    closing = dict(zip(")]}>", "([{<"))
    pairs = set()
    for i, c in enumerate(db, 1):
        if c in stacks:
            stacks[c].append(i)
        elif c in closing:
            opening = closing[c]
            require(bool(stacks[opening]), f"Unmatched bracket at {i}")
            pairs.add((stacks[opening].pop(), i))
        else:
            require(c == ".", f"Unknown structure character: {c}")
    require(not any(stacks.values()), "Unclosed structure bracket")
    return pairs


def crosses(a, b):
    i, j = a
    k, l = b
    return i < k < j < l or k < i < l < j


def read_bpp(path, n):
    out = {}
    for line in path.read_text().splitlines():
        if line.startswith(">") or not line.strip():
            continue
        i, j, p = line.split()
        key, p = (int(i), int(j)), float(p)
        require(1 <= key[0] < key[1] <= n and 0 <= p <= 1, "Invalid BPP")
        require(key not in out, "Repeated BPP pair")
        out[key] = p
    return out


def read_context(path, n):
    out = {}
    for line in path.read_text().splitlines():
        fields = line.split()
        if fields and fields[0] in CONTEXTS:
            require(fields[0] not in out, "Repeated context")
            out[fields[0]] = np.array([float(p) for p in fields[1:]])
    require(set(out) == set(CONTEXTS), "Missing context")
    require(all(len(v) == n for v in out.values()), "Context length mismatch")
    require(all(np.all((v >= 0) & (v <= 1)) for v in out.values()), "Invalid context probability")
    return out


def read_access(path, n):
    metadata, out = {}, {}
    for line in path.read_text().splitlines():
        if line.startswith("# ") and "=" in line:
            k, v = line[2:].split("=", 1)
            metadata[k] = v
        elif line.strip() and not line.startswith(("#", ">")):
            fields = line.split("\t")
            require(len(fields) == 3 and fields[2].startswith("byloop:"), "Missing --by-loop values")
            pos = int(fields[0]) + 1  # CLI start coordinates are zero-based.
            totals = dict((int(w), float(p)) for w, p in
                          (item.split(",") for item in fields[1].strip(";").split(";")))
            components = {}
            for item in fields[2][7:].strip(";").split(";"):
                w, *values = item.split(",")
                components[int(w)] = {k: float(v) for k, v in (q.split("=") for q in values)}
            require(set(totals) == set(components), "Mismatched window lengths")
            for w, total in totals.items():
                require((pos, w) not in out, "Repeated accessibility window")
                require(set(components[w]) == set(LOOPS), "Missing loop component")
                require(1 <= pos <= n - w + 1, "Invalid window coordinate")
                out[pos, w] = {"total": total, **components[w]}
    require(set(out) == {(i, w) for w in (1, 4, 8) for i in range(1, n - w + 2)},
            "Incomplete accessibility output")
    require(all(0 <= p <= 1 for v in out.values() for p in v.values()), "Invalid accessibility")
    require(metadata["byloop_value"] == "joint_probability", "Unexpected loop accessibility semantics")
    return metadata, out


def representative(path, n, marker):
    candidates = [line.split()[0] for line in path.read_text().splitlines() if marker in line]
    require(len(candidates) == 1 and len(candidates[0]) == n, f"Missing representative: {path.name}")
    return pairs_of(candidates[0])


def read_setting(data, setting, n, scaffold):
    bpp = read_bpp(data / f"{setting}.bpp.tsv", n)
    ctx = read_context(data / f"{setting}.capr", n)
    metadata, acc = read_access(data / f"{setting}.racc", n)
    samples = []
    for line in (data / f"{setting}.samp").read_text().splitlines():
        fields = line.split("\t")
        if len(fields) == 2 and fields[0].isdigit():
            count, db = int(fields[0]), fields[1]
            require(count > 0 and len(db) == n, "Invalid sample")
            samples.append((count, pairs_of(db)))
    require(sum(c for c, _ in samples) == 1000, "Expected 1,000 sample draws")
    require(len(samples) == len({frozenset(s) for _, s in samples}), "Repeated sample structure row")
    mfe = representative(data / f"{setting}.lc.out", n, " (-")
    mea = representative(data / f"{setting}.mea", n, "EA=")
    require("gamma=1)" in (data / f"{setting}.mea").read_text(), "Expected MEA gamma=1")
    require(metadata["scaffold"] == scaffold, "Scaffold metadata differs from input")
    g = pairs_of(scaffold)
    require(all(g <= s for _, s in samples) and g <= mfe and g <= mea, "Scaffold was not retained")
    context_sum_error = float(np.max(np.abs(sum(ctx.values()) - 1)))
    byloop_sum_error = max(abs(v["total"] - sum(v[k] for k in LOOPS)) for v in acc.values())
    single_base_error = max(abs(acc[i, 1][k] - ctx[c][i - 1])
                            for i in range(1, n + 1) for k, c in LOOPS.items())
    scaffold_error = max(abs(bpp.get(p, 0) - 1) for p in g)
    for error in (context_sum_error, byloop_sum_error, single_base_error, scaffold_error):
        require(error < 1e-10, "Saved outputs fail probability/scaffold consistency check")
    return dict(bpp=bpp, context=ctx, access=acc, samples=samples, mfe=mfe, mea=mea,
                metadata=metadata, validation=dict(context_sum_max_error=context_sum_error,
                byloop_sum_max_error=byloop_sum_error, window1_context_max_error=single_base_error,
                scaffold_pair_max_error=scaffold_error))


def maximum(values):
    value = max(values.values())
    pos = min(i for i, v in values.items() if v == value)
    return dict(position=pos, probability=value)


def summarize(setting, stems, n):
    ctx, acc, samples = setting["context"], setting["access"], setting["samples"]
    out = {"stems": {}, "sample_draws": sum(c for c, _ in samples),
           "unique_sample_structures": len(samples), "largest_structure_count": max(c for c, _ in samples)}
    for name, pairs in stems.items():
        probs = [setting["bpp"].get(p, 0) for p in sorted(pairs)]
        count = sum(c for c, s in samples if pairs <= s)
        out["stems"][name] = dict(pairs=[list(p) for p in sorted(pairs)],
            pair_probabilities=probs, pair_probability_min=min(probs), pair_probability_max=max(probs),
            mean_pair_probability=float(np.mean(probs)), complete_sample_count=count,
            complete_sample_fraction=count / out["sample_draws"],
            any_pair_sample_count=sum(c for c, s in samples if pairs & s),
            individual_pair_sample_counts=[sum(c for c, s in samples if p in s) for p in sorted(pairs)],
            mfe_pairs_present=[list(p) for p in sorted(pairs & setting["mfe"])],
            mea_pairs_present=[list(p) for p in sorted(pairs & setting["mea"])])
    out["pseudoloop"] = dict(
        global_peak=maximum({i: float(ctx["PseudoLoop"][i - 1]) for i in range(1, n + 1)}),
        region_peak=maximum({i: float(ctx["PseudoLoop"][i - 1]) for i in range(REGION[0], REGION[1] + 1)}),
        selected_positions={str(i): float(ctx["PseudoLoop"][i - 1]) for i in (55, 123, 308)})
    out["accessibility"] = dict(window_length=WINDOW,
        region_component_peaks={k: maximum({i: acc[i, WINDOW][k] for i in range(REGION[0], REGION[1] + 1)})
                                for k in ("total", *LOOPS)},
        selected_windows={str(i): {"end": i + WINDOW - 1, **acc[i, WINDOW]}
                          for i in (43, 123, 128, 130, 142, 152)})
    return out


def compare(a, b):
    lo, hi = REGION
    context_errors = []
    bpp_errors = []
    for c in CONTEXTS:
        for i, (p, q) in enumerate(zip(a["context"][c], b["context"][c]), 1):
            context_errors.append((abs(float(p - q)), c, i))
    for pair in a["bpp"].keys() | b["bpp"].keys():
        bpp_errors.append((abs(a["bpp"].get(pair, 0) - b["bpp"].get(pair, 0)), pair))
    out = dict(context_max_abs_difference=max(context_errors),
               context_region_max_abs_difference=max(e for e in context_errors if lo <= e[2] <= hi),
               bpp_max_abs_difference_in_stored_files=max(bpp_errors),
               bpp_region_max_abs_difference_in_stored_files=max(e for e in bpp_errors if lo <= e[1][0] <= e[1][1] <= hi),
               bpp_missing_row_handling="zero for comparison of stored files; low probabilities may be omitted",
               accessibility={})
    for w in (1, 4, 8):
        values = [(abs(v["total"] - b["access"][i, ww]["total"]), i)
                  for (i, ww), v in a["access"].items() if ww == w]
        components = [(abs(v[k] - b["access"][i, ww][k]), i, k)
                      for (i, ww), v in a["access"].items() if ww == w for k in LOOPS]
        out["accessibility"][str(w)] = dict(total_global_max=max(values),
            total_region_max=max(e for e in values if lo <= e[1] <= hi),
            component_region_max=max(e for e in components if lo <= e[1] <= hi))
    return out


def style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(pad=2)


def stack(ax, positions, values, labels):
    # Draw one continuous shape per context, so sparse hatching does not restart
    # in every narrow base column. Each base/window occupies a unit-width bin.
    edges = np.append(np.asarray(positions) - 0.5, positions[-1] + 0.5)
    bottom = np.zeros(len(positions) + 1)
    for label in labels:
        v = np.append(values[label], values[label][-1])
        ax.fill_between(edges, bottom, bottom + v, step="post", color=COLORS[label],
                        hatch=HATCHES[label], edgecolor=MUTED, linewidth=0, zorder=2)
        bottom += v
    ax.set_ylim(0, 1.01)
    ax.set_yticks((0, 0.5, 1), ("0", "0.5", "1"))
    style(ax)


def make_plot(out, rec, setting, stats):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7,
        "axes.linewidth": 0.6, "axes.edgecolor": MUTED, "text.color": INK,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.major.width": 0.5,
        "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2,
        "hatch.linewidth": 0.3, "pdf.fonttype": 42,
        "axes.labelsize": 7, "axes.labelcolor": MUTED})
    n, lo, hi = len(rec["seq"]), *REGION
    g = {tuple(p) for p in rec["nested"]}
    fig = plt.figure(figsize=(WIDTH * MM, HEIGHT * MM), dpi=DPI)

    def text(x, y, value, **kwargs):
        return fig.text(x / WIDTH, y / HEIGHT, value, va="top", **kwargs)

    def axes(x, y, width, height):
        return fig.add_axes((x / WIDTH, y / HEIGHT, width / WIDTH, height / HEIGHT))

    def title(x, y, letter, value):
        text(x, y, letter, fontsize=9, fontweight="bold")
        text(x + 7, y - 0.15, value, fontsize=7.5, fontweight="bold")

    text(8, 99, f"Rat Terc  |  {n} nt  |  G: {len(g)} canonical nested reference pairs  |  beam width 200",
         fontsize=8, fontweight="bold")
    title(1, 93, "A", "Pair-probability overview")
    pair_styles = (("G", G_COLOR, "solid"),
                   ("New, no crossing", NEW_COLOR, (0, (1, 1))),
                   ("Crossing G", CROSS_COLOR, (0, (4, 2))))
    fig.legend(handles=[Line2D([], [], color=c, lw=1.2, linestyle=ls, label=label)
                        for label, c, ls in pair_styles],
               loc="upper right", bbox_to_anchor=(176 / WIDTH, 93.5 / HEIGHT),
               ncol=3, frameon=False, fontsize=7, handlelength=1.8,
               handletextpad=0.4, columnspacing=1.1, borderaxespad=0)
    axa = axes(12, 77, 164, 12)
    axa.set_xlim(0.5, n + 0.5)
    axa.set_ylim(-10, 180)
    axa.axhline(0, color=MUTED, lw=0.5)
    displayed = {"scaffold": 0, "new_no_crossing": 0, "new_crossing_G": 0}
    groups = {k: [] for k in displayed}
    for pair, p in setting["bpp"].items():
        if p < CUTOFF:
            continue
        kind = "scaffold" if pair in g else "new_crossing_G" if any(crosses(pair, q) for q in g) else "new_no_crossing"
        groups[kind].append((pair, p))
    for kind, (_, color, linestyle) in zip(groups, pair_styles):
        for (i, j), p in sorted(groups[kind], key=lambda v: v[1]):
            axa.add_patch(Arc(((i + j) / 2, 0), j - i, (j - i) * 0.9,
                             theta1=0, theta2=180, color=color, lw=0.65,
                             linestyle=linestyle, alpha=p))
            displayed[kind] += 1
    # Separate solid/dotted marks distinguish the two reference groups even
    # without color; their complete pair identities are drawn in D.
    reference_marks = {}
    for name, level, linestyle in (("outer", -5, "solid"), ("inner", -9, (0, (1, 1)))):
        pairs = stats["stems"][name]["pairs"]
        arms = [(min(p[k] for p in pairs), max(p[k] for p in pairs)) for k in (0, 1)]
        reference_marks[name] = arms
        for start, end in arms:
            axa.plot((start, end), (level, level), color=INK, lw=1.5,
                     linestyle=linestyle, solid_capstyle="butt")
        labels = ((40, 120), (167, 120)) if name == "outer" else ((90, 155), (115, 155))
        for (start, end), (x, y) in zip(arms, labels):
            axa.annotate("O" if name == "outer" else "I", xy=((start + end) / 2, level),
                         xytext=(x, y), ha="center", va="center", fontsize=7,
                         bbox=dict(facecolor="white", edgecolor="none", pad=.3),
                         arrowprops=dict(arrowstyle="->", color=INK, lw=.55,
                                         shrinkA=2, shrinkB=1))
    axa.set_xticks((1, 50, 100, 150, 200, 250, 300, 350, n))
    axa.get_xticklabels()[-1].set_ha("right")
    axa.set_yticks([])
    for spine in axa.spines.values():
        spine.set_visible(False)

    title(1, 70, "B", "Structural contexts")
    fig.legend(handles=[Patch(facecolor=COLORS[c], hatch=HATCHES[c], edgecolor=MUTED,
                              linewidth=0, label=c) for c in CONTEXTS],
               loc="upper left", bbox_to_anchor=(12 / WIDTH, 66 / HEIGHT),
               ncol=7, frameon=False, fontsize=7, handlelength=1.0,
               handletextpad=0.35, columnspacing=1.1, borderaxespad=0)
    axfull = axes(12, 54, 164, 8)
    stack(axfull, np.arange(1, n + 1), setting["context"], CONTEXTS)
    axfull.set_xlim(0.5, n + 0.5)
    axfull.set_yticks((0, 1), ("0", "1"))
    axfull.set_xticks((1, 50, 100, 150, 200, 250, 300, 350, n))
    axfull.get_xticklabels()[-1].set_ha("right")
    axfull.add_patch(Rectangle((lo - 0.5, 0), hi - lo + 1, 1.01, fill=False,
                              edgecolor=INK, lw=0.8, linestyle=(0, (3, 2)),
                              zorder=4, clip_on=False))
    peak = stats["pseudoloop"]["region_peak"]
    text(12, 47.5, f"{lo}–{hi} enlarged")
    text(120, 47.5, f"Base {peak['position']}: PseudoLoop {peak['probability']:.3f}", ha="right")
    positions = np.arange(lo, hi + 1)
    axzoom = axes(12, 32, 108, 11)
    stack(axzoom, positions, {c: setting["context"][c][lo - 1:hi] for c in CONTEXTS}, CONTEXTS)
    axzoom.set_xlim(lo - 0.5, hi + 0.5)
    axzoom.set_xticks((40, 60, 80, 100, 120, 140, 160))
    axzoom.set_ylabel("Probability", labelpad=3)

    selected = stats["accessibility"]["selected_windows"]["123"]
    start, end = 123, selected["end"]
    # This bracket spans bases in B. C's guide instead marks a window START.
    axzoom.plot((start - 0.5, start - 0.5, end + 0.5, end + 0.5),
                (1.015, 1.08, 1.08, 1.015), transform=axzoom.get_xaxis_transform(),
                color=INK, lw=0.7, clip_on=False, zorder=4)
    axzoom.axvspan(start - 0.5, end + 0.5, fill=False, edgecolor=INK, lw=0.8, zorder=4)
    axzoom.axvline(start, color=INK, lw=0.6, linestyle=(0, (2, 2)), zorder=4)

    title(1, 25.5, "C", "4-nt window accessibility by loop type")
    text(12, 21.5, f"{start}–{end}: total {selected['total']:.3f}; PseudoLoop {selected['P']:.3f}")
    axc = axes(12, 7, 108, 11)
    values = {c: [setting["access"][i, WINDOW][k] for i in positions] for k, c in LOOPS.items()}
    stack(axc, positions, values, tuple(LOOPS.values()))
    axc.plot(positions, [setting["access"][i, WINDOW]["total"] for i in positions],
             color=INK, lw=0.9, zorder=3)
    axc.set_xlim(lo - 0.5, hi + 0.5)
    axc.set_xticks((40, 60, 80, 100, 120, 140, 160))
    axc.set_xlabel("Start of 4-nt window (nt)", labelpad=2)
    axc.set_ylabel("Probability", labelpad=3)
    axc.axvline(start, color=INK, lw=0.6, linestyle=(0, (2, 2)), zorder=4)
    axc.plot(start, selected["total"], "o", color=INK, ms=2.5, zorder=5)
    axc.annotate("", xy=(start, selected["total"]), xytext=(start - 2, 1.14),
                 arrowprops=dict(arrowstyle="->", lw=0.7, color=INK), annotation_clip=False)
    title(127, 47.5, "D", "Reference pairs and samples")
    # Pair ladder: both arms are explicit 1-based coordinates. The extra
    # right-arm base between the two groups is obtained from the saved pairs.
    column_positions = {"outer": np.linspace(128, 155, len(stats["stems"]["outer"]["pairs"])),
                        "inner": np.linspace(164.5, 175.3, len(stats["stems"]["inner"]["pairs"]))}
    for name, label in (("outer", "Outer"), ("inner", "Inner")):
        stem = stats["stems"][name]
        columns = column_positions[name]
        text(float(np.mean(columns)) if name == "outer" else 177, 43,
             f"{'O' if name == 'outer' else 'I'}: {label} ({len(stem['pairs'])})",
             ha="center" if name == "outer" else "right", fontweight="bold")
        text(float(np.mean(columns)) if name == "outer" else 177, 39.5,
             f"p {stem['pair_probability_min']:.3f}–{stem['pair_probability_max']:.3f}",
             ha="center" if name == "outer" else "right")
        for x, (i, j) in zip(columns, stem["pairs"]):
            text(x, 35.5, str(i), ha="center")
            text(x, 29.5, str(j), ha="center")
            fig.add_artist(Line2D([x / WIDTH] * 2, [30.5 / HEIGHT, 32.5 / HEIGHT],
                                 transform=fig.transFigure, color=CROSS_COLOR, lw=1,
                                 linestyle="solid" if name == "outer" else (0, (1, 1))))
    outer_j = min(p[1] for p in stats["stems"]["outer"]["pairs"])
    inner_j = max(p[1] for p in stats["stems"]["inner"]["pairs"])
    gap = list(range(inner_j + 1, outer_j))
    require(len(gap) == 1, "Expected one reference bulged base")
    text(159.75, 26.5, f"{gap[0]}: bulge", ha="center")
    text(129, 23, "Samples with every listed pair")
    # No per-pair average: each bar requires the complete specified set.
    # Fixed labels and percentage columns keep even the shorter bars readable.
    for y, name in ((19, "outer"), (15.5, "inner"), (12, "both")):
        fraction = stats["stems"][name]["complete_sample_fraction"]
        text(129, y + 1.1, f"{name.title()}: all {len(stats['stems'][name]['pairs'])}")
        fig.add_artist(Rectangle((148 / WIDTH, (y - 0.5) / HEIGHT),
                                16 * fraction / WIDTH, 1.7 / HEIGHT,
                                transform=fig.transFigure, facecolor=CROSS_COLOR,
                                edgecolor=INK, linewidth=0.4))
        text(176, y + 1.1, f"{fraction:.1%}", ha="right")
    text(158, 8.7, "Outer", ha="center")
    text(172, 8.7, "Inner", ha="center")
    for y, key, label in ((5.7, "mfe_pairs_present", "MFE"),
                           (2.7, "mea_pairs_present", "MEA")):
        text(129, y, label)
        for x, name in ((158, "outer"), (172, "inner")):
            stem = stats["stems"][name]
            text(x, y, f"{len(stem[key])}/{len(stem['pairs'])}", ha="center")

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    require(all(t.get_fontsize() >= 6 for t in texts), "Figure has text smaller than 6 pt")
    bounds = fig.bbox
    outside = [t.get_text() for t in texts if not bounds.contains(*t.get_window_extent(renderer).get_points()[0])
               or not bounds.contains(*t.get_window_extent(renderer).get_points()[1])]
    require(not outside, f"Text extends beyond figure: {outside}")
    overlaps = []
    for i, a in enumerate(texts):
        for b in texts[i + 1:]:
            if a.get_window_extent(renderer).overlaps(b.get_window_extent(renderer)):
                overlaps.append((a.get_text(), b.get_text()))
    require(not overlaps, f"Text overlap: {overlaps}")
    # Verify the actual coordinate transforms, not just the manually set limits.
    aligned = max(abs(axzoom.transData.transform((i, 0))[0] -
                      axc.transData.transform((i, 0))[0]) for i in positions)
    require(aligned < 1e-10, "B zoom and C coordinates are not aligned")
    for suffix in ("pdf", "png"):
        fig.savefig(out / f"figure4.{suffix}", dpi=DPI, facecolor="white")
    # NODISPLAY suppresses txtwrite output on the installed Ghostscript; use
    # its real text device and require a nonempty list of embedded text spans.
    gs_argv = ["gs", "-q", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0",
               "-sOutputFile=-", str(out / "figure4.pdf")]
    gs_result = subprocess.run(gs_argv, capture_output=True, text=True, check=True)
    pdf_sizes = [float(v) for v in re.findall(r'<span[^>]* size="([0-9.]+)"', gs_result.stdout)]
    require(pdf_sizes and min(pdf_sizes) >= 6, "PDF has absent text or embedded text below 6 pt")
    validation = dict(width_mm=WIDTH, height_mm=HEIGHT, png_dpi=DPI,
        min_font_pt=min(t.get_fontsize() for t in texts), visible_text_items=len(texts),
        pdf_embedded_font_sizes_pt=sorted(set(pdf_sizes)), pdf_text_span_count=len(pdf_sizes),
        pdf_font_inspection=gs_argv[:-1] + ["figure4.pdf"],
        text_outside_figure=outside, overlapping_text_pairs=overlaps,
        aligned_Bzoom_C_max_horizontal_difference_px=aligned,
        aligned_Bzoom_C_xlimits=list(axzoom.get_xlim()),
        aligned_Bzoom_C_width_mm=108, displayed_pair_counts=displayed,
        pair_display_cutoff=CUTOFF, arc_opacity="pair probability, without rescaling",
        reference_marks=reference_marks, reference_bulged_base=gap[0],
        context_hatches=HATCHES, direct_context_labels=["PseudoLoop in B and C callouts"],
        pair_linestyles={"scaffold": "solid", "new_no_crossing": "dotted", "new_crossing_G": "dashed"})
    plt.close(fig)
    return validation


def write_caption(out, rec, stats, differences):
    stems = stats["stems"]
    selected = stats["accessibility"]["selected_windows"]["123"]
    # Display upper bounds rounded upward, rather than claim an incorrectly
    # rounded maximum. Preserve full precision in the provenance and README.
    ctx_bound = math.ceil(differences["context_region_max_abs_difference"][0] * 1e5) / 1e5
    acc_bound = math.ceil(differences["accessibility"]["4"]["total_region_max"][0] * 1e5) / 1e5
    counts = [stems[name]["complete_sample_count"] for name in ("outer", "inner", "both")]
    n = len(rec["seq"])
    caption = f"""<!-- 日本語：ラットTerc（NR_001567.1、{n} nt）をbeam幅200で解析した。
GはbpRNA参照構造の入れ子の部分から選んだcanonicalな61対を含む。 -->
**Figure 3. Structural contexts and accessibility in rat telomerase RNA.**
Rat Terc (NR_001567.1; {n} nt) was analyzed at beam width 200.
Scaffold G contains the 61 canonical AU, GC, or GU pairs in the nested part of the bpRNA reference structure.

<!-- 日本語：Aは確率0.01以上の対を示す。
灰色の実線はG、青色の点線はGと交差しない新しい対、橙色の破線はGと交差する対を表す。
弧の不透明度は対確率そのものとし、重なった弧は濃く見えることがある。
OとIの指示線は外側と内側の参照対群の両腕を示す。
Dに全9対の座標と対確率pの範囲を示す。 -->
**(A)** The positional overview shows pairs with probability ≥0.01.
Grey solid arcs show G; blue dotted arcs show new pairs without crossings; orange dashed arcs show new pairs crossing G.
Arc opacity equals pair probability; overlapping arcs may appear darker.
O and I point to both arms of the outer and inner reference groups.
D lists all nine reference pairs and their probability ranges, where p denotes pair probability.

<!-- 日本語：Bは各塩基の七つのcontext確率を示す。
PseudoLoopはMethodsで定義したpseudoknotに対応する非対合contextである。
Internalには疎な斜線を付ける。
拡大図とCは同じ40–160の座標を共有する。
Bの括弧は塩基123–126、Cの補助線は窓の開始位置123を示す。 -->
**(B)** Stacks show seven context probabilities per base.
PseudoLoop is the pseudoknot-associated unpaired context defined in Methods.
Sparse hatching identifies Internal.
The enlarged track and C share coordinates 40–160.
The bracket in B spans bases 123–126; the guide in C marks the window start, 123.

<!-- 日本語：Cの色付き成分は、4塩基すべてが非対合で同じloopの種類に属する同時確率である。
黒線は六成分の和である。
窓123–126のaccessibilityは{selected['total']:.3f}で、そのうちPseudoLoopは{selected['P']:.3f}を占める。 -->
**(C)** Each colored component is the joint probability that all four bases are unpaired and assigned to that loop type.
The black line is their sum.
Window 123–126 has accessibility {selected['total']:.3f}, including a PseudoLoop contribution of {selected['P']:.3f}.

<!-- 日本語：Dは外側6対、内側3対、または全9対を含む標本の割合を示す。
1,000回のBoltzmann samplingで該当する標本は順に{counts[0]}、{counts[1]}、{counts[2]}個だった。
MFEは自由エネルギー最小構造、MEAは期待正確度最大構造で、表は含まれる参照対の数を示す。
MEAはgamma=1を使う。
40–160のcontextと4塩基accessibilityのexactからの最大絶対差は、それぞれ{ctx_bound:.5f}未満と{acc_bound:.5f}未満である。 -->
**(D)** Bars show the fractions of samples containing all six outer reference pairs, all three inner pairs, or all nine.
The corresponding counts are {counts[0]}, {counts[1]}, and {counts[2]} of {stats['sample_draws']:,} Boltzmann samples.
The table reports reference pairs present in the minimum-free-energy (MFE) and maximum expected accuracy (MEA) structures.
MEA uses γ=1.
Within coordinates 40–160, maximum absolute differences from exact calculations are <{ctx_bound:.5f} for contexts and <{acc_bound:.5f} for four-base accessibility (window start coordinates).
"""
    (out / "caption.md").write_text(caption)
    english = re.sub(r"<!--.*?-->", "", caption, flags=re.S)
    return len(english.split())


def write_tables(out, a, b, stats):
    rows = ["stem\ti\tj\tb200_probability\texact_probability\tb200_sample_pair_count"]
    for name in ("outer", "inner"):
        stem = stats["stems"][name]
        for pair, count in zip(stem["pairs"], stem["individual_pair_sample_counts"]):
            key = tuple(pair)
            rows.append(f"{name}\t{pair[0]}\t{pair[1]}\t{a['bpp'][key]:.17g}\t{b['bpp'][key]:.17g}\t{count}")
    (out / "reference_pairs.tsv").write_text("\n".join(rows) + "\n")
    rows = ["start_1based\tend_1based\ttotal\t" + "\t".join(LOOPS.values())]
    for i in range(REGION[0], REGION[1] + 1):
        v = a["access"][i, WINDOW]
        rows.append(f"{i}\t{i + WINDOW - 1}\t" + "\t".join(f"{v[k]:.17g}" for k in ("total", *LOOPS)))
    (out / "accessibility_region.tsv").write_text("\n".join(rows) + "\n")


def write_readme(out, data, stats, differences, validation, caption_words):
    stems = stats["stems"]
    windows = stats["accessibility"]["selected_windows"]
    lines = ["# Figure 4: ラット Terc（20261017b）", "",
        "[PDF](figure4.pdf) · [PNG、300 dpi](figure4.png) · [caption](caption.md) · [provenance](provenance.json)", "",
        "bpRNA_RFAM_4317、ラット Terc（RefSeq NR_001567.1、RNAcentral URS00001ED027）、388 nt。",
        "識別情報は保存済みデータの README に従う。外部データの取得と RNA CLI の再実行は行わない。",
        "G は bpRNA 参照構造の入れ子の73対のうち、AU/UA、GC/CG、GU/UG の61対（j − i > 3）。",
        "除いた12対は非canonical対であり、参照の入れ子の全対を保持する条件ではない。",
        "入力 in.txt と rec.json の配列・対集合を照合して確認する。表示する予測はすべて beam 幅200。", "",
        "## 入力・版・再生成", "",
        f"保存済みの入力・出力: `{relative(data)}`。",
        "入力・出力の版と保存条件は、データの README と provenance.json の software lock に従う。",
        (data / "README.md").read_text().strip(),
        "元図 results/figures/figure4_terc_20261017/ を変更せず、別ディレクトリへ出力する。", "",
        "実験 repo から次で再生成できる。", "", "```bash",
        f"python3 -B scripts/make_figure4_terc.py --data {relative(data)} --out {relative(out)}", "```", "",
        "依存: Python 3、NumPy、Matplotlib、Ghostscript（gs）。キャッシュと一時ファイルは /tmp を使う。",
        "共有領域の新規ファイルは PDF、PNG、caption、README、2 TSV、provenance の7個だけ。", "",
        "## 図の読み方", "",
        "- A: 全長の対確率の概観。0.01以上の対を表示し、弧の alpha を確率そのものとする。",
        "  G は灰色の実線、新しい非交差対は青の点線、G と交差する対は橙の破線。重なりの濃さを個別対の確率として読まない。",
        "  基線下の実線は外側6対の両腕、点線は内側3対の両腕。D が各対の座標と予測対確率 p の範囲を示す。",
        "- B: 七つの context の周辺確率。上段は全長、破線の40–160を下段で拡大する。Internal は疎な斜線を共有する。",
        "  括弧と枠は塩基123–126、破線は塩基123。PseudoLoop の値は拡大区間の最大値であり、全長の最大値ではない。",
        "- C: 開始位置40–160の4塩基窓。六成分は全4塩基が非対合かつ同じ loop の種類に属する同時確率。",
        "  六成分の和が全体の accessibility（黒線）。B の拡大図と C は同じ左端、108 mm幅、x範囲39.5–160.5を使う。",
        "  C の点・矢印は開始位置123の窓123–126であり、B の塩基123の周辺確率とは別の量。",
        "- D: 上部は bpRNA 由来の参照対の座標図、p の範囲は保存された beam200 の予測対確率。",
        "  外側6対は (63,139)–(68,134)、内側3対は (69,132)、(70,131)、(71,130)。133 の bulge を挟む。",
        "  棒は各群の全対を含む標本の割合であり、Both は全9対の同時出現。ほかの対も含んでよい。",
        "  周辺対確率の平均や積ではない。MFE/MEA 表は含まれる参照対の個数を示す。MEA は γ=1。", "",
        "## 科学的な集計", "",
        "| 参照対の群 | 対確率の範囲 | 全対を含む標本 | MFE | MEA (γ=1) |",
        "|---|---:|---:|---:|---:|"]
    for name in ("outer", "inner", "both"):
        v = stems[name]
        lines.append(f"| {name} ({len(v['pairs'])}対) | {v['pair_probability_min']:.6f}–{v['pair_probability_max']:.6f} | "
                     f"{v['complete_sample_count']}/{stats['sample_draws']} ({v['complete_sample_fraction']:.1%}) | "
                     f"{len(v['mfe_pairs_present'])}/{len(v['pairs'])} | {len(v['mea_pairs_present'])}/{len(v['pairs'])} |")
    p = stats["pseudoloop"]
    w = windows["123"]
    lines += ["", f"PseudoLoop の拡大区間内の最大は塩基{p['region_peak']['position']}で{p['region_peak']['probability']:.9f}。",
        f"全長の最大は塩基{p['global_peak']['position']}で{p['global_peak']['probability']:.9f}。",
        f"窓123–126の total は{w['total']:.9f}、PseudoLoop 成分は{w['P']:.9f}。",
        f"PseudoLoop は total の{w['P'] / w['total']:.3%}を占める。これは予測同時確率の内訳であり、物理的露出の実証ではない。", "",
        "| 4塩基窓 | Total | Exterior | Hairpin | Bulge | Internal | Multiloop | PseudoLoop |",
        "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for i, v in windows.items():
        lines.append(f"| {i}–{v['end']} | " + " | ".join(f"{v[k]:.6f}" for k in ("total", *LOOPS)) + " |")
    lines += ["", "各 loop 成分の40–160での最大値（開始位置）:", "",
        "| Loop | 最大同時確率 | 開始 |", "|---|---:|---:|"]
    for k, c in LOOPS.items():
        peak = stats["accessibility"]["region_component_peaks"][k]
        lines.append(f"| {c} | {peak['probability']:.6f} | {peak['position']} |")
    lines += ["", f"{stats['sample_draws']:,}標本は{stats['unique_sample_structures']}種類。同一構造の最多出現は{stats['largest_structure_count']}回。",
        "標本数は samp の各行の multiplicity で重みづけする。sampling seed は raw README によると20261017。", "",
        "### exact との比較", "",
        f"40–160の全 context の最大絶対差は{differences['context_region_max_abs_difference'][0]:.11f}。",
        f"開始位置40–160の4塩基 total accessibility の最大絶対差は{differences['accessibility']['4']['total_region_max'][0]:.11f}。",
        "caption の <0.00284 と <0.00191 は上の最大値を小数第5位で切り上げた上限（通常の四捨五入の最大値ではない）。",
        f"六つの4塩基成分の最大絶対差は{differences['accessibility']['4']['component_region_max'][0]:.11f}で、total と区別する。",
        f"全長では context の最大差{differences['context_max_abs_difference'][0]:.9f}、4塩基 total の最大差{differences['accessibility']['4']['total_global_max'][0]:.9f}。",
        "caption の比較を全長へ広げない。BPP の保存ファイルの欠行は比較に限り0とし、低確率の未保存値を復元しない。",
        "この条件付き解析例は、無条件の予測精度や生物学的機能を評価しない。", "",
        "## 開発上の検証", "",
        "G の61対と canonical フィルタ、参照9対、全標本と MFE/MEA での G 保持を照合する。",
        "context の和=1、loop 別 accessibility の和=total、1塩基窓と context の一致を beam200/exact の両方で検査する。",
        f"図は{WIDTH} × {HEIGHT} mm、PNG は{DPI} dpi。Matplotlib の可視文字は{validation['min_font_pt']:g} pt以上。",
        f"PDF 内の文字は Ghostscript txtwrite で{validation['pdf_text_span_count']} span を検出し、文字サイズは{validation['pdf_embedded_font_sizes_pt']} pt。",
        "-dNODISPLAY は使用環境で空の出力となったため付けず、空でない span を必須として検査する。",
        f"B 拡大と C の全塩基座標の横位置差は{validation['aligned_Bzoom_C_max_horizontal_difference_px']:g} px。",
        "文字同士の重複0、図外の文字0。PNG の最終目視結果と変更対応は下記の担当報告に記録する。",
        f"caption は英語{caption_words}語（日本語コメントを除く空白分割）。", "",
        "## レビュー対応", "",
        "- F4-1 / G1: 178 mm幅・100 mm高・通常7 pt。B 拡大/C の範囲と幅をそろえ、123の補助線を付ける。",
        "- F4-2 / C9 / L56: 全9対の座標図を D に追加し、外側6対・内側3対と133の bulge を明示。",
        "- F4-3: 弧の濃さは概観、9参照対の確率は D の数値範囲として示す。alpha=確率、閾値0.01は維持。",
        "- F4-4 / C10: 窓123–126の total/PseudoLoop を数値・矢印で示し、同時確率の内訳と定義する。",
        "- F4-5 / L61: 全対を含む標本の割合と個数を明記。MFE/MEA は対数表にする。",
        "- F4-6 / C2 / C11 / L60: G は canonical な61対。PseudoLoop は Methods の定義へ接続し、exact 比較を40–160に限定。",
        "- G2 / G3: 対の線種、Internal の疎な斜線、PseudoLoop の直接注記。独立した太字9 ptの A–D。", "",
        "## 未変更・著者への引継ぎ", "",
        "原稿、参考文献、他担当の図、raw、旧図、外部実装は変更しない。git commit もしない。",
        "各参照対の値は reference_pairs.tsv、C の全数表は accessibility_region.tsv、入力・出力 hash と環境は provenance.json。",
        "組み替え案 B に伴う図番号・caption の原稿への組込み・最終掲載 PDF の検査は統合担当へ引き継ぐ。",
        "新しい科学的主張・条件の変更はない。著者には新配置と caption の採用を確認いただく。",
        "同じ変更・未変更・判断事項・検証結果を論文 repo の textbook/notes/codex_logs/figure4_terc_round2.last.txt に記録する。",
        "docs/progress.md は担当外なので直接変更せず、この報告から統合担当が転記する。", ""]
    (out / "README.md").write_text("\n".join(lines))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, default=DEFAULT_DATA)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--software-lock", type=Path, default=ROOT / "manifests/software-pilot-febedd98.json")
    ap.add_argument("--identity-source", type=Path)
    args = ap.parse_args()
    data, out = args.data.resolve(), args.out.resolve()
    require(data != out and data not in out.parents, "Write figure artifacts outside the raw-data directory")
    rec = json.loads((data / "rec.json").read_text())
    header, seq, scaffold = (line for line in (data / "in.txt").read_text().splitlines() if line)
    n = len(seq)
    g = {tuple(p) for p in rec["nested"]}
    pk = {tuple(p) for p in rec["pk"]}
    require(rec["id"] == "bpRNA_RFAM_4317" and seq == rec["seq"] and n == len(scaffold) == 388, "Wrong RNA")
    require(g == pairs_of(scaffold) and len(g) == 61, "Unexpected scaffold")
    require(g <= pairs_of(rec["ref"]) and pk <= pairs_of(rec["ref"]), "Pairs absent from reference")
    canonical = {"AU", "UA", "CG", "GC", "GU", "UG"}
    reference_nested = pairs_of(rec["ref"]) - pk
    require(g == {p for p in reference_nested
                  if p[1] - p[0] > 3 and seq[p[0] - 1] + seq[p[1] - 1] in canonical},
            "Scaffold differs from canonical nested reference pairs")
    noncanonical = {p for p in reference_nested if seq[p[0] - 1] + seq[p[1] - 1] not in canonical}
    require(len(reference_nested) == 73 and len(noncanonical) == 12 and reference_nested - g == noncanonical,
            "Unexpected reference pair exclusions")
    require(len(pk) == 9 and all(any(crosses(p, q) for q in g) for p in pk), "Wrong reference pseudoknot")
    require(not any(crosses(p, q) for p in g for q in g), "Scaffold contains crossings")
    stems = dict(outer={p for p in pk if p[0] <= 68}, inner={p for p in pk if p[0] >= 69}, both=pk)
    require(len(stems["outer"]) == 6 and len(stems["inner"]) == 3, "Unexpected reference stems")
    settings = {name: read_setting(data, name, n, scaffold) for name in ("b200", "exact")}
    require(settings["b200"]["metadata"]["beam"] == "200", "Wrong beam width")
    require(settings["exact"]["metadata"]["beam"] == "exact", "Expected exact comparison")
    stats = {name: summarize(s, stems, n) for name, s in settings.items()}
    differences = compare(settings["b200"], settings["exact"])
    out.mkdir(parents=True, exist_ok=True)
    validation = make_plot(out, rec, settings["b200"], stats["b200"])
    caption_words = write_caption(out, rec, stats["b200"], differences)
    write_tables(out, settings["b200"], settings["exact"], stats["b200"])
    write_readme(out, data, stats["b200"], differences, validation, caption_words)
    sources = [data / name for name in ("in.txt", "rec.json", "README.md")]
    sources += [data / f"{s}.{ext}" for s in settings for ext in ("bpp.tsv", "capr", "racc", "samp", "lc.out", "mea")]
    identity_source = args.identity_source or data.parent / "README.md"
    sources += [identity_source, args.software_lock]
    lock = json.loads(sources[-1].read_text())
    provenance = dict(schema_version=1, generated_utc=datetime.now(timezone.utc).isoformat(),
        purpose="Figure 4 from existing rat Terc outputs", rna=dict(id=rec["id"], length=n,
        species="Rattus norvegicus", gene="Terc", refseq="NR_001567.1", rnacentral="URS00001ED027",
        identity_source=relative(identity_source)),
        software=dict(source_commit=lock["source_commit"], lock=relative(sources[-1]),
                      association_source=relative(data / "README.md")),
        environment=dict(python=platform.python_version(), numpy=np.__version__, matplotlib=matplotlib.__version__),
        regeneration=dict(argv=sys.argv, script=relative(Path(__file__)), script_sha256=sha(Path(__file__))),
        sources=[dict(path=relative(p), sha256=sha(p), bytes=p.stat().st_size) for p in sources],
        design=dict(beam=200, scaffold_pairs=len(g), reference_pk_pairs=len(pk), region=list(REGION),
            scaffold_definition="canonical AU/UA, GC/CG, GU/UG nested bpRNA reference pairs; j-i>3",
            reference_nested_pair_count=len(reference_nested),
            excluded_noncanonical_reference_pairs=[list(p) for p in sorted(noncanonical)],
            sampling_seed_reported_in_raw_readme=20261017, mea_gamma=1,
            coordinates="1-based inclusive; accessibility position is window start",
            accessibility_window=WINDOW, colors=COLORS,
            complete_stem_definition="every listed reference pair occurs; other pairs are allowed",
            sample_counts="multiplicity-weighted; do not treat each unique row as one draw",
            caption_exact_agreement_scope="probabilities near reference pseudoknot, positions 40–160",
            raw_command_line_receipt=("frozen packed-task JSONL and tar" if lock["source_commit"].startswith("89baea52")
                                      else "settings taken from raw README and Raccess metadata")),
        raw_metadata={k: v["metadata"] for k, v in settings.items()},
        summaries=stats, beam_exact_comparison=differences,
        validation=dict(saved_output_consistency={k: v["validation"] for k, v in settings.items()}, figure=validation,
                        caption_english_words=caption_words),
        outputs=[dict(path=p.name, sha256=sha(p), bytes=p.stat().st_size)
                 for p in sorted(out.iterdir()) if p.is_file() and p.name != "provenance.json"])
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(dict(output=relative(out), figure_validation=validation, caption_english_words=caption_words,
                          b200_summary=stats["b200"], beam_exact_comparison=differences), indent=2))


if __name__ == "__main__":
    main()
