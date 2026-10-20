#!/usr/bin/env python3
"""Publication-size Figures 2/3 and SI companions, selected by --review-layout.

Called by make_figures.py, with its validated aggregate data and helpers.
Historical drawing functions and their defaults remain untouched.
"""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import re
import subprocess

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Arc
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter, NullLocator
import numpy as np

MM = 1 / 25.4
BEAMS = (20, 50, 100, 200, 300)
BEAM_COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#222222")
BEAM_MARKERS = ("o", "s", "^", "D", "v")
PROGRAM_COLORS = ("#0072B2", "#D55E00", "#009E73")
OUTPUT_LABELS = (
    "fold: partition function + BPPs + MFE",
    "profile: loop profile",
    "access: accessibility",
)
MODE_LABELS = ("Inside + preliminary ranking", "MFE program",
               "BPP program", "Context program")
SHORT_SERIES = ("Empty G", "Long stem", "3 helices", "Free bases together",
                "12 helices", "Nested multiloops", "Balanced multiloops")


def style(base):
    base.style()
    plt.rcParams.update({
        "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "lines.markersize": 3, "lines.linewidth": .8, "axes.labelpad": 2,
        "xtick.major.pad": 2, "ytick.major.pad": 2, "axes.titlepad": 3,
        "axes.formatter.use_mathtext": False, "axes.unicode_minus": True,
        "xtick.minor.visible": False, "ytick.minor.visible": False,
    })


def text_mm(fig, x, y, text, **kwargs):
    return fig.text(x / (fig.get_figwidth() / MM), y / (fig.get_figheight() / MM),
                    text, **kwargs)


def panel(fig, letter, x, y):
    text_mm(fig, x, y, letter, fontsize=9, fontweight="bold", va="baseline")


def beam_handles(cohort=None):
    return [Line2D([], [], color=c, marker=m, ls="none", markersize=3.3,
                   label=str(b)) for b, c, m in zip(BEAMS, BEAM_COLORS, BEAM_MARKERS)
            if cohort != "short" or b != 300]


def legend_mm(fig, handles, x, y, ncol, **kwargs):
    kwargs.setdefault("columnspacing", .9)
    return fig.legend(handles=handles, loc="center",
                      bbox_to_anchor=(x / (fig.get_figwidth() / MM),
                                      y / (fig.get_figheight() / MM)),
                      ncol=ncol, frameon=False, handlelength=1.4,
                      handletextpad=.5, borderaxespad=0, **kwargs)


def grouped_points(ax, data, base, cohort, metric, statistic="p90", kind="single"):
    """The same beam position, color and shape in every accuracy panel."""
    beams = base.BEAMS[cohort]
    beams = tuple(b for b in beams if b != 400)
    for j, tier in enumerate(base.COHORTS[cohort]):
        for beam in beams:
            i = BEAMS.index(beam)
            x = j + (i - (len(beams) - 1) / 2) * .16
            q = data.q[(cohort, tier, beam, metric)]
            c, m = BEAM_COLORS[i], BEAM_MARKERS[i]
            if kind == "distribution":
                ax.plot([x, x], [q["min"], q["max"]], color=c, lw=.7, alpha=.75)
                ax.plot([x, x], [q["p10"], q["p90"]], color=c, lw=2)
                y = q["median"]
            else:
                y = q[statistic]
            dx = -.04 if kind == "pairing" else .04 if kind == "context" else 0
            ax.plot(x + dx, y, marker=m, ls="none", color=c, markersize=3,
                    markerfacecolor="white" if kind == "context" else c,
                    markeredgewidth=.65, zorder=4)
    ax.set_xlim(-.6, 3.6)
    ax.set_xticks(range(4), [base.bin_label(t) for t in base.COHORTS[cohort]])
    base.clean(ax)


def error_axis(ax):
    # No floor is applied to observations: zero remains zero on the linear band.
    ax.set_yscale("symlog", linthresh=.000001, linscale=.5)
    ax.set_ylim(-.00000012, .4)
    ax.set_yticks([0, .0001, .01, .3], ["0", "0.0001", "0.01", "0.3"])
    ax.yaxis.set_minor_locator(NullLocator())
    ax.axhspan(0, .000001, color=".92", zorder=0)


def accuracy_figure(data, base):
    fig = plt.figure(figsize=(178 * MM, 92 * MM))
    for col, cohort in enumerate(base.COHORTS):
        left = (18, 107)[col]
        heading = ("Exact reference · <400 nt" if col == 0 else
                   "Finite reference (b=400) · 400+ nt")
        text_mm(fig, left + 34, 89, heading, ha="center", fontweight="bold")
        ax = base.mm_axes(fig, left, 64, 68, 18)
        grouped_points(ax, data, base, cohort, base.DELTA, kind="distribution")
        ax.set_yscale("symlog", linthresh=.001, linscale=.5)
        ax.set_ylim(-4, 400)
        ax.set_yticks([-1, 0, 1, 100], ["−1", "0", "1", "100"])
        ax.yaxis.set_minor_locator(NullLocator())
        ax.axhspan(-.001, .001, color=".92", zorder=0)
        ax.set_ylabel("ΔF (kcal/mol)")
        ax.tick_params(axis="x", labelbottom=False)
        if col == 0:
            ax.set_title("Free-energy difference", loc="left")

        ax = base.mm_axes(fig, left, 34, 68, 15)
        for metric, kind in zip(base.TV, ("pairing", "context")):
            grouped_points(ax, data, base, cohort, metric, kind=kind)
        error_axis(ax)
        ax.tick_params(axis="x", labelbottom=False)
        ax.set_ylabel("p90 mean TV")
        ax.set_title("Pairing and context probabilities", loc="left")

        ax = base.mm_axes(fig, left, 10, 68, 13)
        grouped_points(ax, data, base, cohort, base.ACCESS[-1])
        error_axis(ax)
        ax.set_ylabel("90th percentile")
        ax.set_xlabel("RNA length group (nt)")
        ax.set_title("16-nt accessibility (MAE)", loc="left")
    for letter, y in zip("ABC", (84, 51, 25)):
        panel(fig, letter, 1, y)
    key = base.mm_axes(fig, 106, 83, 4, 4)
    key.plot([.5, .5], [0, 1], color=".3", lw=.7)
    key.plot([.5, .5], [.2, .8], color=".3", lw=2)
    key.plot(.5, .5, marker="D", color=".3", ms=3)
    key.set(xlim=(0, 1), ylim=(-.1, 1.1))
    key.axis("off")
    text_mm(fig, 111, 84.2, "Range; 10–90%; symbol: median")
    text_mm(fig, 18, 60.5, "Linear region: −0.001 to 0.001 kcal/mol")
    text_mm(fig, 28, 56, "Beam width b", ha="center")
    legend_mm(fig, beam_handles(), 105, 56.6, 5)
    legend_mm(fig, [Line2D([], [], color=".2", marker="o", ls="none", label="Pairing: filled"),
                    Line2D([], [], color=".2", marker="o", markerfacecolor="white", ls="none",
                           label="Contexts: open")], 130, 29, 2)
    text_mm(fig, 18, 28.2, "B/C linear region: 0 to 1e−6")
    return fig


def threshold_figure(data, base):
    fig = plt.figure(figsize=(165 * MM, 55 * MM))
    for col, cohort in enumerate(base.COHORTS):
        left = (18, 101)[col]
        ax = base.mm_axes(fig, left, 16, 61, 29)
        beams = tuple(b for b in base.BEAMS[cohort] if b != 400)
        for j, tier in enumerate(base.COHORTS[cohort]):
            for beam in beams:
                i = BEAMS.index(beam)
                x = j + (i - (len(beams) - 1) / 2) * .16
                y = data.thresholds[(cohort, tier, beam)]["positive_gt_0.01"]
                ax.plot(x, y, marker=BEAM_MARKERS[i], color=BEAM_COLORS[i], ls="none")
        ax.set_xlim(-.6, 3.6)
        labels = [f"{base.bin_label(t)}\nN={data.q[(cohort, t, beams[0], base.DELTA)]['n']:,}"
                  for t in base.COHORTS[cohort]]
        ax.set_xticks(range(4), labels)
        ax.set_ylim(-3, 105)
        ax.set_yticks([0, 25, 50, 75, 100])
        ax.set_ylabel("RNAs with ΔF > 0.01 (%)")
        ax.set_xlabel("RNA length group (nt)")
        ax.set_title("Exact reference" if col == 0 else "Finite reference (b=400)",
                     loc="left", fontweight="bold")
        base.clean(ax)
        panel(fig, "AB"[col], (2, 85)[col], 48)
    text_mm(fig, 25, 3, "Beam width b", ha="center")
    legend_mm(fig, beam_handles(), 101, 3.8, 5)
    return fig


def rna_points(ax, data, base, metric):
    for i, program in enumerate(base.PROGRAMS):
        rows = sorted((r for r in data.rna if r["program"] == program), key=lambda r: int(r["n"]))
        x = np.array([int(r["n"]) / 1000 for r in rows])
        y = np.array([base.value(r, metric + "_median") for r in rows])
        lo = np.array([base.value(r, metric + "_min") for r in rows])
        hi = np.array([base.value(r, metric + "_max") for r in rows])
        ax.errorbar(x, y, yerr=(y - lo, hi - y), fmt=("o", "s", "^")[i],
                    ls="none", markersize=(3.6, 5, 3)[i], markeredgewidth=.7,
                    markerfacecolor=PROGRAM_COLORS[i] if i == 2 else "white",
                    color=PROGRAM_COLORS[i], capsize=(2, 3, 1)[i], elinewidth=.8,
                    capthick=.7, zorder=3 + i)
    maximum = max(base.value(r, metric + "_max") for r in data.rna)
    ax.set_ylim(0, maximum * 1.24)
    ax.set_xlim(3, 31)
    ax.set_xticks([5, 15, 25])
    ax.set_xlabel("RNA length (kb)")
    ax.set_ylabel("Wall time (s)" if metric == "wall_seconds" else "Peak RSS (GiB)")
    if metric == "wall_seconds":
        ax.set_yticks([0, 500, 1000])
    else:
        ax.set_yticks([0, 2, 4])
    row = max((r for r in data.rna if r["program"] == "LinearCParty"), key=lambda r: int(r["n"]))
    cost = base.value(row, metric + "_median")
    amount = f"{cost / 60:.1f} min" if metric == "wall_seconds" else f"{cost:.2f} GiB"
    ax.annotate(f"{int(row['n']) / 1000:.1f} kb: {amount}",
                xy=(int(row["n"]) / 1000, cost),
                xytext=(int(row["n"]) / 1000, maximum * 1.12),
                textcoords="data", va="center", ha="right", fontsize=7,
                arrowprops=dict(arrowstyle="->", lw=.6, color=".3",
                                connectionstyle="angle,angleA=0,angleB=90,rad=0"))
    base.clean(ax)


def resource_figure(data, base):
    fig = plt.figure(figsize=(178 * MM, 55 * MM))
    for left, metric in ((13, "wall_seconds"), (64, "rss_gib")):
        ax = base.mm_axes(fig, left, 21, 35, 24)
        rna_points(ax, data, base, metric)
        ax.set_title("RNAcentral: time" if metric == "wall_seconds" else "RNAcentral: memory", loc="left")
    for letter, x in zip("ABC", (1, 52, 104)):
        panel(fig, letter, x, 49)
    text_mm(fig, 115, 49, "2,000 → 4,000 nt", va="baseline")
    ax = base.mm_axes(fig, 143, 21, 32, 24)
    ratios = []
    for j, series in enumerate(base.SERIES):
        for i, metric in enumerate(("wall_seconds", "rss_gib")):
            numerator = data.synthetic_summary[("LinearCParty", series, 4000, metric)]["median"]
            denominator = data.synthetic_summary[("LinearCParty", series, 2000, metric)]["median"]
            ratio = numerator / denominator
            ratios.append(ratio)
            ax.plot(ratio, 6 - j + (.18 if i == 0 else -.18), marker=("o", "D")[i],
                    color=("#0072B2", ".2")[i], ls="none", markersize=3.6,
                    markerfacecolor="#0072B2" if i == 0 else "white")
            data.audit.append(dict(panel="3C_ratio", program="LinearCParty", series=series,
                                   metric=metric, numerator=numerator, denominator=denominator,
                                   value=ratio, length_from=2000, length_to=4000,
                                   statistic="ratio of medians of three within-seed medians"))
    ax.set_xlim(1, max(4, math.ceil(max(ratios))))
    ax.set_xticks(range(1, int(ax.get_xlim()[1]) + 1))
    ax.set_ylim(-.6, 6.6)
    ax.set_yticks(range(7), SHORT_SERIES[::-1])
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Cost ratio")
    ax.axvline(2, color=".7", lw=.7, ls=":", zorder=0)
    ax.grid(axis="x", color=".9")
    legend_mm(fig, [Line2D([], [], marker=m, ls="none", color=c, label=label,
                           markersize=s, markerfacecolor=c if i == 2 else "white")
                    for i, (m, c, s, label) in enumerate(zip(("o", "s", "^"), PROGRAM_COLORS,
                                                            (3.6, 5, 3), OUTPUT_LABELS))],
              53, 6.4, 1, labelspacing=.35)
    legend_mm(fig, [Line2D([], [], marker="o", ls="none", color="#0072B2", label="Time"),
                    Line2D([], [], marker="D", ls="none", color=".2", markerfacecolor="white",
                           label="Memory")], 143, 6.4, 2)
    return fig


def memory_figure(data, base):
    fig = plt.figure(figsize=(165 * MM, 65 * MM))
    ax = base.mm_axes(fig, 18, 20, 143, 36)
    for i, stage in enumerate(base.STAGES):
        rows = sorted((r for r in data.memory if r["mode"] == stage), key=lambda r: int(r["n"]))
        ax.plot([int(r["n"]) / 1000 for r in rows], [base.value(r, "rss_gib") for r in rows],
                ls="none", marker=base.MARKERS[i], markersize=(5.5, 4.5, 3.5, 2.5)[i],
                markerfacecolor="white", markeredgewidth=.9, color=base.COLORS[i], zorder=2 + i)
    ax.set_xlim(3, 31)
    ax.set_ylim(0, max(base.value(r, "rss_gib") for r in data.memory) * 1.13)
    ax.set_xticks([5, 10, 15, 20, 25, 30])
    ax.set_xlabel("RNA length (kb)")
    ax.set_ylabel("Peak RSS (GiB)")
    ax.set_title("Peak memory of independent processes", loc="left", fontweight="bold")
    base.clean(ax)
    legend_mm(fig, [Line2D([], [], marker=base.MARKERS[i], ls="none", color=base.COLORS[i],
                           markerfacecolor="white", label=label)
                    for i, label in enumerate(MODE_LABELS)], 83, 6.5, 2, labelspacing=.5)
    return fig


def concurrency_figure(data, base):
    rows = base.read_tsv(data.paths["e2_metrics.tsv"])
    pairs = [r for r in base.read_tsv(data.paths["e2_comparisons.tsv"])
             if r["kind"] in ("parallel_w10_over_w1", "production_over_w1")]
    base.require(len(pairs) == 18, "Expected three RNAs, three programs and two concurrency comparisons")
    provenance = json.loads(data.paths["e2_summary.json"].read_text())["provenance"]["runs"]
    base.require(max(t["observed_max_active_processes"] for t in provenance["parallel_w10"]["tasks"]) == 9,
                 "Refresh the concurrency caption for the observed process count")
    for label, workers in (("parallel_w1", 1), ("parallel_w10", 10)):
        base.require(provenance[label]["signature"]["commit"] == data.e2_commit,
                     "Concurrency/production source revisions differ")
        selected = [r for r in rows if r["run"] == label]
        base.require(len(selected) == 9 and all(int(r["workers"]) == workers for r in selected),
                     "Concurrency cohort/worker mismatch")
    lookup = {(r["run"], r["name"], r["program"]): r for r in rows}
    base.require(len(lookup) == len(rows), "Duplicate concurrency source row")
    for pair in pairs:
        a, b = [lookup[(pair["run_" + side], pair["name_" + side], pair["program"])]
                for side in ("a", "b")]
        base.require(pair["name_a"] == pair["name_b"] and a["n"] == b["n"],
                     "Concurrency input mismatch")
        for metric in ("wall_seconds", "rss_gib"):
            for row in (a, b):
                repetitions = json.loads(row[metric + "_repetitions"])
                base.require(len(repetitions) == 3, "Concurrency repetition mismatch")
                base.close(np.median(repetitions), base.value(row, metric + "_median"),
                           "Concurrency replicate median")
            numerator, denominator = (base.value(r, metric + "_median") for r in (b, a))
            base.close(numerator / denominator, base.value(pair, metric + "_ratio"),
                       "Concurrency ratio")
            data.audit.append(dict(panel="SI_parallel", comparison=pair["kind"],
                                   program=pair["program"], name=a["name"], length=int(a["n"]),
                                   metric=metric, numerator=numerator, denominator=denominator,
                                   value=numerator / denominator))
    data.checks["e2_concurrency_ratios_checked"] = 2 * len(pairs)
    fig = plt.figure(figsize=(165 * MM, 72 * MM))
    for col, metric in enumerate(("wall_seconds", "rss_gib")):
        ax = base.mm_axes(fig, 17 + 82 * col, 23, 62, 37)
        for i, program in enumerate(base.PROGRAMS):
            for kind, filled in (("parallel_w10_over_w1", True), ("production_over_w1", False)):
                group = [r for r in pairs if r["program"] == program and r["kind"] == kind]
                ax.plot([float(r["n_a"]) / 1000 for r in group],
                        [base.value(r, metric + "_ratio") for r in group],
                        ls="none", marker=("o", "s", "^")[i],
                        markersize=(4.5, 5.5, 3.5)[i], color=PROGRAM_COLORS[i],
                        markerfacecolor=PROGRAM_COLORS[i] if filled else "white",
                        markeredgewidth=.8, zorder=2 + i)
        values = [base.value(r, metric + "_ratio") for r in pairs] + [1]
        span = max(values) - min(values)
        padding = max(.0005, .15 * span)
        ax.set_ylim(min(values) - padding, max(values) + padding)
        ax.set_xlim(4.5, 7.3)
        ax.set_xticks([4.793, 5.482, 6.990], ["4.793", "5.482", "6.990"])
        ax.set_xlabel("RNA length (kb)")
        ax.set_ylabel("Time ratio to w1" if col == 0 else "RSS ratio to w1")
        ax.axhline(1, color=".5", lw=.7, ls=":")
        base.clean(ax)
        panel(fig, "AB"[col], 1 + 82 * col, 66)
    legend_mm(fig, [Line2D([], [], marker=m, ls="none", color=c, label=p)
                    for m, c, p in zip(("o", "s", "^"), PROGRAM_COLORS, base.PROGRAMS)],
              83, 12, 3)
    legend_mm(fig, [Line2D([], [], marker="o", ls="none", color=".3", label=label,
                           markerfacecolor=".3" if filled else "white")
                    for label, filled in (("w10/w1", True), ("Main/w1", False))],
              83, 4, 2)
    return fig


def accessibility_detail(data, base):
    fig = plt.figure(figsize=(165 * MM, 137 * MM))
    for row, metric in enumerate(base.ACCESS):
        bottom = 107 - row * 30
        panel(fig, "ABCD"[row], 1, bottom + 20)
        for col, cohort in enumerate(base.COHORTS):
            left = (19, 102)[col]
            ax = base.mm_axes(fig, left, bottom, 60, 18)
            grouped_points(ax, data, base, cohort, metric)
            error_axis(ax)
            ax.set_title(f"{base.WINDOWS[row]}-nt accessibility", loc="left")
            ax.set_ylabel("90th percentile\nof mean error")
            if row == 3:
                ax.set_xlabel("RNA length group (nt)")
    text_mm(fig, 49, 133, "Exact reference", ha="center", fontweight="bold")
    text_mm(fig, 132, 133, "Finite reference (b=400)", ha="center", fontweight="bold")
    text_mm(fig, 19, 6.3, "Linear region: 0 to 1e−6; logarithmic above")
    text_mm(fig, 25, 3, "Beam width b", ha="center")
    legend_mm(fig, beam_handles(), 102, 3.8, 5)
    return fig


def load_designs(data, base):
    path = base.ROOT / "data/e2-synthetic-shapes-20261012/inputs.jsonl"
    designs = {}
    for line in path.read_text().splitlines():
        row = json.loads(line)
        if len(row["scaffold"]) == 250 and row["family"] in base.SERIES:
            designs.setdefault(row["family"], row["scaffold"])
    base.require(set(designs) == set(base.SERIES), "Missing saved scaffold designs")
    data.paths["synthetic_designs"] = path
    data.hashes["synthetic_designs"] = base.digest(path)
    expected = json.loads(data.paths["e2_summary.json"].read_text())["provenance"]["dataset_sha256"]["inputs.jsonl"]
    base.require(data.hashes["synthetic_designs"] == expected, "Scaffold designs differ from E2 dataset")
    return designs


def draw_design(ax, scaffold):
    # Saved n=250 scaffolds; all pairs are shown without changing their spacing.
    stack = []
    ax.plot([0, 1], [0, 0], color=".7", lw=.6)
    for j, char in enumerate(scaffold):
        if char == "(":
            stack.append(j)
        elif char == ")":
            i = stack.pop()
            span = (j - i) / len(scaffold)
            ax.add_patch(Arc(((i + j) / 2 / len(scaffold), 0), span, span,
                             theta1=0, theta2=180, color="#8a8a85", lw=.45))
    ax.set_xlim(-.02, 1.02)
    ax.set_ylim(-.02, .52)
    ax.set_axis_off()


def synthetic_detail(data, base, designs):
    fig = plt.figure(figsize=(165 * MM, 178 * MM))
    for row, program in enumerate(base.PROGRAMS):
        bottom = 135 - row * 44
        for col, metric in enumerate(("wall_seconds", "rss_gib")):
            left = (19, 102)[col]
            ax = base.mm_axes(fig, left, bottom, 60, 30)
            base.plot_synthetic(ax, data, program, metric)
            summaries = [q for (p, s, n, m), q in data.synthetic_summary.items() if m == metric]
            lo, hi = min(q["min"] for q in summaries) * .65, max(q["max"] for q in summaries) * 1.4
            ax.set_ylim(lo, hi)
            role = ("fold", "profile (loop profile)", "access (accessibility)")[row]
            ax.set_title(role + (": time" if col == 0 else ": memory"), loc="left")
            # Normal-sized decimal labels, including all logarithmic ticks.
            ticks = [10.**k for k in range(math.floor(math.log10(lo)), math.ceil(math.log10(hi)) + 1)]
            ax.yaxis.set_major_locator(FixedLocator(ticks))
            ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: np.format_float_positional(x, trim="-")))
            ax.yaxis.set_minor_locator(NullLocator())
            ax.yaxis.set_minor_formatter(NullFormatter())
            ax.xaxis.set_minor_locator(NullLocator())
            ax.xaxis.set_minor_formatter(NullFormatter())
            panel(fig, "ABCDEF"[row * 2 + col], (1, 84)[col], bottom + 33)
    for i, series in enumerate(base.SERIES):
        col, row = i % 4, i // 4
        left, bottom = 3 + col * 41, 17 - row * 14
        ax = base.mm_axes(fig, left, bottom, 17, 8)
        draw_design(ax, designs[series])
        text_mm(fig, left + 20, bottom + 6, SHORT_SERIES[i].replace("Free bases together", "Free bases\ntogether")
                .replace("Balanced multiloops", "Balanced\nmultiloops").replace("Nested multiloops", "Nested\nmultiloops"),
                va="center", color=".15")
        key = base.mm_axes(fig, left + 20, bottom - .5, 15, 2)
        key.plot([0, 1], [.5, .5], color=base.COLORS[i], ls=base.LINES[i],
                 marker=base.MARKERS[i], markevery=[0], markersize=3, markerfacecolor="white")
        key.set_axis_off()
    return fig


def pdf_text_check(path):
    result = subprocess.run(["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=txtwrite",
                             "-dTextFormat=0", "-sOutputFile=-", str(path)],
                            check=True, capture_output=True, text=True)
    # Ghostscript's XML contains size attributes on span elements.
    sizes = [float(x) for x in re.findall(r'\bsize="([0-9.]+)"', result.stdout)]
    if not sizes:
        raise ValueError(f"Ghostscript found no text sizes in {path}")
    if min(sizes) < 6:
        raise ValueError(f"PDF text below 6 pt: {path}: {min(sizes)}")
    return dict(minimum_pt=min(sizes), maximum_pt=max(sizes), spans=len(sizes),
                sizes_pt=sorted(set(sizes)), method="Ghostscript txtwrite TextFormat=0")


def evaluation_counts(data, base, out):
    rows = []
    for cohort, tiers in base.COHORTS.items():
        for tier in tiers:
            for window, metric in zip(base.WINDOWS, base.ACCESS):
                summaries = [data.q[(cohort, tier, b, metric)] for b in base.BEAMS[cohort]]
                domains = {(q["n"], q["missing"], q["eligible_windows"]) for q in summaries}
                base.require(len(domains) == 1, "Length-specific evaluation domains vary with beam")
                q = summaries[0]
                rows.append(dict(cohort=cohort, length_group=base.bin_label(tier), window_nt=window,
                                 total_rna=q["rna"], evaluated_rna=q["n"], excluded_rna=q["missing"],
                                 eligible_windows=q["eligible_windows"]))
    with (out / "evaluation_counts.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def captions(data, base, out, counts):
    group_counts = "; ".join(
        f"{base.bin_label(t)} nt: {data.q[(c,t,20,base.DELTA)]['rna']:,}"
        for c in base.COHORTS for t in base.COHORTS[c]) if counts else ""
    captions = """# Caption drafts / caption 案

Figure numbers follow the task labels; the manuscript editor can renumber the resource figure after Terc.
各英語段落の直前に、対応する日本語段落を置く。

## Figure 2 — figure2_main.pdf

**Beam-width sensitivity of conditional-ensemble outputs.**
左列は400 nt未満の24,625 RNAをexactと比べ、右列は400–4,381 ntの1,438 RNAを有限幅b=400と比べる。
色と記号は全行でbeam幅bを示す。
幅400の自己比較は表示しない。

The left column compares 24,625 RNAs shorter than 400 nt with exact calculations.
The right column compares 1,438 RNAs of 400–4,381 nt with finite beam width b=400.
Colors and symbol shapes identify beam width b in every row.
The b=400 self-comparison is omitted because its differences are zero by definition.

Aは自由エネルギーの符号付き差ΔF=F_b−F_refを示す。
F_bとF_refは有限幅と参照のensemble free energyである。
細線は全範囲、太線は10–90パーセンタイル、記号は中央値を示す。
灰色帯の−0.001から0.001 kcal/molの間だけ線形とした対称対数軸を使う。

A shows the signed free-energy difference ΔF=F_b−F_ref, where F_b and F_ref are the ensemble free energies at finite beam width and in the reference calculation.
Thin lines show the full range, thick lines the 10th–90th percentiles, and symbols the median across RNAs.
The symmetric logarithmic axis is linear within the shaded band from −0.001 to 0.001 kcal/mol.

Bでは、各RNA内でGの非対合塩基についてtotal variation distance（TV）を平均した後、RNA間の90パーセンタイル（p90）を取る。
塗った記号は相手塩基と非対合状態の分布、白抜き記号は七つのcontextの分布を比べる。
TVは各分布の確率差の絶対値の和の半分である。

B first averages total variation distance (TV) over bases unpaired in scaffold G within each RNA, then takes the 90th percentile (p90) across RNAs in each length group.
Filled symbols compare distributions over pairing partners and the unpaired state; open symbols compare the seven structural-context probabilities.
TV is half the sum of absolute probability differences over the states in each distribution.

BとCは灰色帯の0–0.000001で線形、その上で対数の尺度を使う。
ゼロを正の下限値へ置き換えていない。

B and C use a linear scale in the shaded band from 0 to 0.000001 and a logarithmic scale above it.
Zeros remain zero; no positive floor is applied to the observations.
FOCUS_TV_STATISTICS

Cは、16 nt窓の全塩基が同時に非対合となる確率の平均絶対誤差を示す。
Gの対合塩基を含まない窓についてRNA内で平均し、長さ群ごとにRNA間の90パーセンタイルを取る。
評価可能な窓がないRNAは除き、確率ゼロの窓は含める。
対称対数軸は灰色帯の0–0.000001の間で線形であり、差に正の下限を代入していない。
他の窓長と評価件数は補足図とevaluation_counts.tsvに示す。

C shows the 90th percentile across RNAs of the mean absolute error (MAE) in 16-nt window accessibility, the probability that all bases in the window are unpaired together.
Within each RNA, the mean includes only windows with no base paired in G.
RNAs without an eligible window are excluded; eligible windows with zero probability are included.
The symmetric logarithmic axis is linear in the shaded band from 0 to 0.000001, and no positive floor is applied to the data.
The supplementary accessibility figure shows all four tested window lengths.

横軸は長さ群であり、群の中では設定を等間隔に並べる。
hとdは局所探索の上限であり、h=d=max(b,30)とした。
そのため、これはbeam幅と局所上限を合わせた設定の感度である。

The horizontal positions group RNAs by length and space the tested settings equally within each group.
The local search limits h and d are set to max(b,30), so the figure describes sensitivity to the combined settings.

## Figure 3 — figure3_main.pdf

**Time and memory for long RNAs and synthetic controls.**
A/BはRNAcentralの七入力（4,793–29,106 nt）について異なる出力を要求したプログラム全体の費用を示す。
MFEはminimum free energy、BPPはbase-pairing probabilityを指す。
LinearCPartyはMFE・分配関数・BPP、LinearPKCapRはcontext確率、LinearPKRaccessはcontext確率と1/4/8/16 nt窓のtotal accessibilityを計算する。
個別の解析の追加コストは分離しない。

A and B report complete programs with different requested outputs on seven RNAcentral inputs of 4,793–29,106 nt.
MFE denotes minimum free energy and BPP denotes base-pairing probability.
LinearCParty computes MFE, the partition function and BPPs; LinearPKCapR computes context probabilities; LinearPKRaccess computes context probabilities and total accessibility for 1-, 4-, 8- and 16-nt windows.
These totals do not isolate the added cost of a single analysis.

点はwarm-up一回を除いた三測定の中央値、縦線は全範囲である。
一部の範囲は記号より小さい。
長さごとに異なる配列とGを使う。
時間は起動・予備計算・主計算・復元・出力を含み、scaffold構築とbuildを除く。
メモリはプロセス全体のpeak resident set size（RSS）である。

Symbols show medians of three measurements after one excluded warm-up; vertical bars show the full range.
Some ranges are smaller than the symbols.
Each length represents a different RNA and scaffold G.
Wall time includes startup, preliminary ranking, main computation, traceback and output, and excludes scaffold construction and builds.
Memory is the peak resident set size (RSS) of the entire process.

CはLinearCPartyの合成七系列における2,000→4,000 ntの時間とRSSの比を示す。
比は、各seedの三反復の中央値を三seed間でまとめた中央値同士から計算する。
点線は比2を示す。
比は測定区間の観測であり、漸近計算量の上界ではない。
系列ごとに配列とGの両方が変わるため、Gの形だけの因果効果とは解釈しない。
全長さと三プログラムの曲線は補足図に置く。

C shows the ratios of LinearCParty time and RSS between 4,000 and 2,000 nt for seven synthetic series.
At each length, we first take the median of three repetitions for each seed, then the median across three seeds.
Ratios divide the resulting 4,000-nt median by the 2,000-nt median; the dotted line marks a ratio of two.
These ratios describe the measured interval and are not asymptotic bounds.
Both sequence and scaffold change between series, so the comparisons do not identify a causal effect of scaffold shape alone.
The supplementary synthetic figure shows all five lengths and all three programs.

## Figure S1 — figure2_d_si.pdf

**RNAs above a positive free-energy-difference threshold.**
A/BはFigure 2と同じexact/有限参照の群を使う。
ΔF=F_b−F_refが0.01 kcal/molを超えるRNA数を、各長さ群のRNA数で割った割合を示す。
この片側の条件は負の差を数えない。
幅400の自己比較はゼロになるため表示しない。
色と記号、群内の等間隔の設定はFigure 2と同じである。
各長さ群のNは割合の分母となるRNA数を示す。

A and B use the exact and finite references from Figure 2, respectively.
Each percentage divides the number of RNAs with ΔF=F_b−F_ref greater than 0.01 kcal/mol by all RNAs in that length group.
This one-sided criterion excludes negative differences.
The b=400 self-comparison is omitted because its percentage is zero by definition.
Colors, symbols and equally spaced settings within each length group match Figure 2.
N gives the number of RNAs in each length group, the denominator of the percentage.
GROUP_COUNTS

## Figure S2 — figure3_mode_memory_si.pdf

**Peak memory of independent calculation modes.**
各点は別のRNAとGの組を示す。
insideは順位付けの予備計算を含み、他のモードは必要なforward計算と準備をそれぞれ含む。
各モードは独立したプロセスで、warm-upなしの一回測定である。
異なるプロセスのpeakを足したり引いたりしてメモリの内訳とはしない。

Each length represents a different RNA and scaffold G.
Inside includes preliminary ranking; the other modes each include their own required forward calculation and preparation.
Each mode runs in an independent process and has one measurement without warm-up.
The peaks cannot be added or subtracted to infer memory components.

## Additional SI — figure2_accessibility_si.pdf

**Accessibility differences by window length and RNA length.**
A–Dは1/4/8/16 nt窓の同時非対合確率について、RNA内の平均絶対誤差を求め、RNA間の90パーセンタイルを示す。
左列は400 nt未満をexactと比べ、右列は400–4,381 ntを有限幅b=400と比べる。
色と記号はbeam幅bを示し、b=400の自己比較は除く。
Gで対合している塩基を含む窓は除き、確率ゼロの窓は含める。
灰色帯の0–0.000001は線形、その上は対数とする。
hとdは局所探索の上限で、h=d=max(b,30)とした。
長さ群・窓長ごとの評価数と除外数はevaluation_counts.tsvに保存する。

A–D show the 90th percentile across RNAs of within-RNA mean absolute errors in 1-, 4-, 8- and 16-nt joint window accessibility.
Left panels compare RNAs shorter than 400 nt with exact calculations.
Right panels use b=400 as a finite reference for RNAs of 400–4,381 nt.
Colors and symbols identify beam width b; the b=400 self-comparison is omitted.
Eligible windows contain no base paired in G; eligible zero-probability windows are included.
The shaded band is linear from 0 to 0.000001; larger values use a logarithmic scale.
The local search limits h and d are set to max(b,30).
Evaluation counts and excluded RNAs for every length group and window length are provided in evaluation_counts.tsv.

## Additional SI — figure3_synthetic_all_cli_si.pdf

**Measured costs for seven synthetic sequence–scaffold series.**
両軸は対数であり、RNA長は底2、時間とpeak RSSは底10を使う。
点は各seed三反復の中央値を三seed間でまとめた中央値、縦線は九測定の全範囲である。
線は同じ系列の測定点を結ぶ。
全156入力のうち、共通の五長さを持つ七系列×五長さ×三seedの105入力を示す。
残る51入力は追加の分岐深さ対照である。
非空のGはfloor(0.15n)対を持ち、配列のGC含量は50%である。
下の弧図は保存済みの250 ntのGを同じ縮尺で示す。

Both axes are logarithmic: RNA length uses base 2, and time and peak RSS use base 10.
Symbols show medians across three seeds of within-seed medians of three repetitions; vertical bars span all nine measurements.
Lines connect measured values within each series.
The figure shows 105 of the 156 inputs: seven series, five shared lengths and three seeds.
The remaining 51 inputs are additional branch-depth controls.
Each nonempty G contains floor(0.15n) pairs, where n is sequence length; sequences have 50% GC content.
The arc sketches show the saved 250-nt scaffolds on a common coordinate scale.

## Additional SI — figure3_parallel_impact_si.pdf

**Concurrent and single-process measurements.**
A/Bは三つのRNAについて、各入力の三測定の中央値を単独実行w1の中央値で割った時間/RSS比を示す。
塗った記号はw10/w1、白抜き記号は本文の測定/w1を示す。
w10は最大十プロセスを許可し、観測された同時実行の最大数は九だった。
各設定は同じ資源型の異なるホストで実行し、点線は比1を示す。
これらの比で他の測定を補正しない。

A and B show time and RSS ratios for three RNAs, dividing each input's median of three measurements by its one-process (w1) median.
Filled symbols show w10/w1; open symbols show Main/w1.
The w10 setting allowed ten processes and reached nine concurrent processes.
Settings used different hosts of the same resource type; dotted lines mark a ratio of one.
These ratios do not correct the other measurements.
"""
    captions = captions.replace("GROUP_COUNTS", "RNA counts by length group are " + group_counts + ".")
    focus_note = ''
    if counts:
        pairing, context = (data.q[('short', 'nt300_399', 200, m)] for m in base.TV)
        focus_note = (f"\n300–399 nt・b=200のp90はpairing {pairing['p90']:.4f}、context {context['p90']:.4f}である。\n"
                      f"同群・同設定のRNA内平均TVの最大値はそれぞれ{pairing['max']:.3f}と{context['max']:.3f}であり、p90とは異なる。\n\n"
                      f"For 300–399 nt at b=200, p90 is {pairing['p90']:.4f} for pairing and {context['p90']:.4f} for contexts.\n"
                      f"The corresponding maxima of RNA mean TV are {pairing['max']:.3f} and {context['max']:.3f}; B plots p90, not these maxima.\n")
    captions = captions.replace('FOCUS_TV_STATISTICS', focus_note)
    if data.figure3_only:
        blocks = captions.split("\n## ")
        captions = blocks[0] + "".join("\n## " + block for block in blocks[1:]
                                     if not block.startswith(("Figure 2", "Figure S1",
                                                              "Additional SI — figure2_")))
    (out / "caption.md").write_text(captions)


def readme(data, base, out, records, pdf_checks, counts):
    sources = "\n".join(f"| {key} | `{os.path.relpath(path,out)}` | {data.hashes[key]} |"
                        for key, path in data.paths.items())
    output_table = "\n".join(f"| {r['file']} | {r['width_mm']:g} × {r['height_mm']:g} | "
                             f"{pdf_checks[r['file']]['minimum_pt']:g} |"
                             for r in records if r["file"].endswith(".pdf"))
    access_table = "\n".join(f"| {r['cohort']} | {r['length_group']} | {r['window_nt']} | "
                             f"{r['evaluated_rna']:,} | {r['excluded_rna']:,} | {r['eligible_windows']:,} |"
                             for r in counts)
    ratio_rows = [r for r in data.audit if r["panel"] == "3C_ratio"]
    ratio_table = "\n".join(f"| {label} | " + " | ".join(
        f"{next(r['value'] for r in ratio_rows if r['series'] == series and r['metric'] == metric):.5f}"
        for metric in ("wall_seconds", "rss_gib")) + " |"
        for series, label in zip(base.SERIES, SHORT_SERIES))
    time_ratios = [r["value"] for r in ratio_rows if r["metric"] == "wall_seconds"]
    text = f"""# Revised main Figures 2/3 and SI companions — 2026-10-17

E1: beam source d21618c5; exact reference d9804df6.
E2: **{data.e2_commit}** (one revision for RNAcentral, synthetic controls and mode memory).
The resource figure keeps its task filename `figure3_main`; plan B may number it Figure 4 in the manuscript.
No manuscript, Japanese draft, bibliography or external implementation is edited.

## Files and validation

Every figure has a PDF and 300-dpi PNG.
Main figures are exactly 178 mm wide; SI figures are exactly 165 mm wide.
Normal text is 7 pt and independent bold panel letters are 9 pt.
Ghostscript checks the actual PDF text spans, including mathematical text.
All logarithmic ticks use ordinary decimal digits.

| PDF | Size (mm) | Minimum PDF text (pt) |
|---|---|---:|
{output_table}

`caption.md` contains paired Japanese/English caption drafts, one sentence per line.
`plot_data.tsv` retains all checked E1 quantiles (including quartiles), threshold counts, E2 summaries and the exact numerator/denominator of every ratio.
When E1 is rendered, `evaluation_counts.tsv` records each accessibility domain, including excluded RNAs.
`validation.json` records source checks, geometry checks and PDF font sizes.
`manifest.json` records source/script/output hashes.
PNG visual checks and manuscript integration are recorded in `review_report.json` after inspection.

## Review mapping and choices

- F2-1: colors, symbols and within-group beam positions stay fixed in A–C and S1; length groups form the horizontal categories.
- F2-2: main C separates all length groups at 16 nt; the additional SI figure shows all four windows, with a declared linear band and no artificial positive floor.
- F2-3: A uses full range, 10th–90th percentiles and median; quartiles remain in the audit table.
- F2-4: pairing is filled and contexts are open, with a local key; caption defines TV and RNA averaging before the percentile.
- F2-5: beam values are settings within categorical length groups, not numeric line slopes; no lines connect the settings.
- F2-6: b=400 self-comparisons are omitted in all revised accuracy panels and explained in the captions.
- F3-1/F3-3: C shows the observed 4,000/2,000-nt ratio for each named series; all full curves and saved scaffold sketches move to a 165-mm SI companion.
- F3-2: larger open square and smaller filled triangle expose nearly coincident programs without moving the RNA coordinates; longest-input annotations are computed from aggregate medians.
- F3-4: program legend lists requested outputs; caption states whole-program measurement and distinct outputs.
- S1-1/S1-2: full threshold in the axis label; A/B, explicit denominators below, one-sided criterion and omitted self-comparison.
- S2-1/S2-2: explicit independent-process names, preliminary ranking, one measurement, points without connecting lines, and data-derived y limits.
- G1: final physical widths, 7-pt normal labels and decimal logarithmic ticks; embedded PDF font verification.
- G2: beam and program shapes, filled/open metric keys and direct synthetic-series labels supplement color; the saved scaffold sketches use G grey #8a8a85.
- G3: independent bold 9-pt panel letters at fixed offsets; Figure 2 uses one A/B/C per row and explicit reference columns.

The author/manuscript editor must choose final supplementary numbering for the two added detail figures and update figure references after plan B.
The E1 total-accessibility benchmark and historical E2 outputs do not validate the new per-loop accessibility components.
Resource values for the new per-loop option are not claimed.

## Numerical sources

| ID | Path relative to this directory | SHA-256 |
|---|---|---|
{sources}

E1 uses Turner 2004, 37°C, preliminary outer ranking and h=d=max(b,30).
All quantiles weight each evaluable RNA equally and use NumPy's linear quantile interpolation.
E2 uses b=h=d=100, outer ranking, Turner 2004, 37°C, one thread/process, at most 10 concurrent processes on cpu_40, and RLIMIT_AS=8 GiB.
RSS GiB = GNU time KiB / 2^20.
Synthetic ratios divide the medians already shown in the full curves; they are not medians of seed-wise ratios.

| Synthetic series | Time ratio | RSS ratio |
|---|---:|---:|
{ratio_table}

**Integration check:** the time-ratio range from the plotted medians is {min(time_ratios):.2f}–{max(time_ratios):.2f}.
If the manuscript retains a different range from medians of seed-wise ratios, update the text or explicitly distinguish the two statistics.

## Accessibility denominators

| Reference cohort | Length group (nt) | Window (nt) | Evaluated RNA | Excluded RNA | Eligible windows |
|---|---|---:|---:|---:|---:|
{access_table}

## Regenerate

Run from the paper repository:

```bash
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR="$TMPDIR/e2-figure-mpl" python3 experiments/scripts/make_figures.py --review-layout {"--figure3-only " if data.figure3_only else ""}--e2-version {data.e2_commit[:8]} --out {os.path.relpath(out,base.ROOT.parent)}
```

`--e2-version febd2bc1` selects `results/local/e2-febd2bc1-20261016`.
`--e2-version d02d261c` selects `results/local/e2-d02d261c-20261017`.
For arbitrary aggregates, omit the preset and supply `--e2 PATH --e2-revision COMMIT`.
Historical defaults, `--main-compact`, `--figure3-only` and `--figure3-panels` remain available.
Use a new output directory for a different E2 revision, for example `main_20261017/d02d261c/`, to preserve this febd2bc1 rendering.

## Finish d02d261c after every run completes

Do not analyze a partial run or mix source revisions.
The required runs are `e2-rnacentral`, `e2-synthetic-shapes`, `e2-memory-breakdown`, `e2-parallel-impact-w1` and `e2-parallel-impact-w10`, all suffixed `-d02d261c-20261017`.
Check that every planned task's `task.json` has status `success`, then refresh the root manifests with `aggregate_array.py`.
Root `run.json` can be stale until this aggregation; `.execution.lock` alone is not completion evidence.
The analyzer validates artifacts, input identities, repetitions, settings and software provenance without extraction to shared storage.
The command below produces the same aggregate schema as the febd2bc1 input:

```bash
for cohort in rnacentral synthetic-shapes memory-breakdown parallel-impact-w1 parallel-impact-w10
do
  PYTHONDONTWRITEBYTECODE=1 python3 experiments/scripts/aggregate_array.py \\
    --plan "$PWD/experiments/results/runs/e2-${{cohort}}-d02d261c-20261017/array_plan.json"
done
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR="$TMPDIR/e2-figure-mpl" python3 experiments/scripts/analyze_e0_e2_runs.py \\
  --synthetic-run experiments/results/runs/e2-synthetic-shapes-d02d261c-20261017 \\
  --rna-run experiments/results/runs/e2-rnacentral-d02d261c-20261017 \\
  --previous-synthetic-run experiments/results/runs/e2-synthetic-shapes-febd2bc1-20261016 \\
  --previous-rna-run experiments/results/runs/e2-rnacentral-febd2bc1-20261016 \\
  --memory-run experiments/results/runs/e2-memory-breakdown-d02d261c-20261017 \\
  --previous-memory-run experiments/results/runs/e2-memory-breakdown-febd2bc1-20261016 \\
  --parallel-w1-run experiments/results/runs/e2-parallel-impact-w1-d02d261c-20261017 \\
  --parallel-w10-run experiments/results/runs/e2-parallel-impact-w10-d02d261c-20261017 \\
  --workers 10 --revision d02d261c --previous-workers 10 --previous-revision febd2bc1 \\
  --previous-resource cpu_40=1 --skip-cap30-history \\
  --out experiments/results/local/e2-d02d261c-20261017
PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR="$TMPDIR/e2-figure-mpl" python3 experiments/scripts/make_figures.py --review-layout --figure3-only --e2-version d02d261c --out experiments/results/figures/main_d02d261c_20261017
```

Remove the temporary Matplotlib cache after analysis.
No RNA program, scheduler submission or implementation checkout is needed for aggregation or plotting.
"""
    (out / "README.md").write_text(text)


def render_review(args, data, base):
    out = args.out
    style(base)
    records = []
    designs = load_designs(data, base)
    counts = evaluation_counts(data, base, out) if not args.figure3_only else []
    figures = []
    if not args.figure3_only:
        figures.extend(((accuracy_figure, "figure2_main", 178),
                        (threshold_figure, "figure2_d_si", 165),
                        (accessibility_detail, "figure2_accessibility_si", 165)))
    figures.extend(((resource_figure, "figure3_main", 178),
                    (memory_figure, "figure3_mode_memory_si", 165),
                    (lambda d, b: synthetic_detail(d, b, designs), "figure3_synthetic_all_cli_si", 165),
                    (concurrency_figure, "figure3_parallel_impact_si", 165)))
    if args.published_only:
        figures = [f for f in figures if f[1] not in
                   ("figure3_mode_memory_si", "figure3_parallel_impact_si")]
    for draw, name, width in figures:
        base.save(draw(data, base), name, out, width, records, check_text_overlap=True)
        print(f"Rendered {name}", flush=True)
    pdf_checks = {r["file"]: pdf_text_check(out / r["file"])
                  for r in records if r["file"].endswith(".pdf")}
    base.write_audit(data, out)
    captions(data, base, out, counts)
    if args.published_only:
        blocks = (out / "caption.md").read_text().split("\n## ")
        selected = [b for b in blocks[1:] if 'figure3_mode_memory_si.pdf' not in b.splitlines()[0]
                    and 'figure3_parallel_impact_si.pdf' not in b.splitlines()[0]]
        caption = blocks[0] + ''.join('\n## ' + b for b in selected)
        caption = caption.replace('Figure 3 — figure3_main.pdf', 'Figure 4 — figure3_main.pdf')
        caption = caption.replace('Additional SI — figure3_synthetic_all_cli_si.pdf',
                                  'Figure S3 — figure3_synthetic_all_cli_si.pdf')
        caption = caption.replace('Figure numbers follow the task labels; the manuscript editor can renumber the resource figure after Terc.',
                                  'Figure numbers match the adopted manuscript: Figure 4 (resources) and Figure S3 (synthetic series).')
        (out / 'caption.md').write_text(caption)
    readme(data, base, out, records, pdf_checks, counts)
    for key, path in data.paths.items():
        base.require(base.digest(path) == data.hashes[key], f"Source changed: {path}")
    validation = dict(all_passed=True, checks=dict(data.checks), source_files_unchanged=True,
                      text_within_pages=True, text_without_overlap=True,
                      plotted_points_within_axis_limits=True, minimum_normal_font_pt=7,
                      pdf_text_checks=pdf_checks, outputs=records)
    (out / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    derived = ["caption.md", "README.md", "plot_data.tsv", "validation.json"]
    if counts:
        derived.append("evaluation_counts.tsv")
    manifest = dict(layout="review-20261019-r4" if out.name.endswith('_r4') else "review-20261017",
                    e2_revision=data.e2_commit,
                    scripts_sha256={p.name: base.digest(p) for p in (Path(base.__file__), Path(__file__))},
                    numpy=np.__version__, matplotlib=plt.matplotlib.__version__,
                    sources={key: dict(path=os.path.relpath(path,out), sha256=data.hashes[key])
                             for key, path in data.paths.items()}, outputs=records,
                    derived_files={name: base.digest(out / name) for name in derived})
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(dict(out=str(out), figures=len(records)//2, pdf_text_checks=pdf_checks), indent=2))
