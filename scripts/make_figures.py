#!/usr/bin/env python3
"""Draft Figures 2 and 3 from aggregate TSV files; never run RNA software.

Run from either repository with Python >=3.9, NumPy and Matplotlib:
    python3 experiments/scripts/make_figures.py
Only --out is written. Caches go to a temporary directory outside the repository.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
from collections import defaultdict

_CACHE = tempfile.TemporaryDirectory(prefix="linearcp-figures-", dir="/tmp")
os.environ["MPLCONFIGDIR"] = _CACHE.name
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.text import Text
import numpy as np
from analyze_e2_synthetic import check_revision

ROOT = Path(__file__).resolve().parents[1]
MM = 1 / 25.4
DELTA = "delta_free_energy_kcal_mol"
TV = ("base_state_tv_mean_Gfree", "context_tv_mean_Gfree")
WINDOWS = (1, 4, 8, 16)
ACCESS = tuple(f"access_w{w}_probability_mae_eligible" for w in WINDOWS)
METRICS = (DELTA, *TV, *ACCESS)
COHORTS = {
    "short": ("nt50_99", "nt100_199", "nt200_299", "nt300_399"),
    "long": ("nt400_799", "nt800_1599", "nt1600_3199", "nt3200_4381"),
}
REFERENCES = {"short": "exact", "long": "beam400_outer"}
REF_LABEL = {"short": "<400 nt: exact reference", "long": "400+ nt: b=400 reference"}
BEAMS = {"short": (20, 50, 100, 200), "long": (20, 50, 100, 200, 300, 400)}
# Okabe-Ito colors; redundant markers, line styles and hatches preserve identities.
COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#222222", "#56B4E9")
MARKERS = ("o", "s", "^", "D", "v", "P", "X")
HATCHES = ("", "///", "\\\\", "xx", "..", "++")
LINES = ("-", "--", "-.", ":", (0, (5, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1)), (0, (1, 2)))
PROGRAMS = ("LinearCParty", "LinearPKCapR", "LinearPKRaccess")
SERIES = (
    "empty", "long_stem", "separated_helices", "separated_concentrated",
    "separated_12", "deep_branching", "balanced_branching",
)
SERIES_LABELS = (
    "Empty G", "Long stem", "Separated helices", "Concentrated free bases",
    "12 helices", "Deep branching", "Balanced branching",
)
STAGES = ("inside_outer", "mfe_outer", "outside_bpp_outer", "loop_profile_outer")
STAGE_LABELS = ("Inside (outer)", "MFE", "Outside + BPP", "Loop profile")


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    checksum = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            checksum.update(chunk)
    return checksum.hexdigest()


def read_tsv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def value(row, key):
    return math.nan if row[key] in ("NA", "", "nan") else float(row[key])


def bin_label(tier):
    return tier.removeprefix("nt").replace("_", "–")


def close(a, b, context):
    require(math.isclose(a, b, rel_tol=2e-11, abs_tol=1e-12), f"Mismatch: {context}: {a} != {b}")


class Data:
    def __init__(self, args):
        self.paths = {}
        self.figure3_only = args.figure3_only
        self.figure3_panels = args.figure3_panels
        self.e2_revision = args.e2_revision
        for cohort, directory in (() if self.figure3_only else (("short", args.e1_short), ("long", args.e1_long))):
            for filename in ("metrics.tsv", "quantiles.tsv"):
                self.paths[f"e1_{cohort}_{filename}"] = directory / filename
        for filename in ("metrics.tsv", "comparisons.tsv", "memory_metrics.tsv", "summary.json"):
            self.paths[f"e2_{filename}"] = args.e2 / filename
        self.paths["g_features"] = args.g_features
        self.hashes = {key: digest(path) for key, path in self.paths.items()}
        self.audit = []
        self.checks = defaultdict(int)
        self.q = {}
        self.thresholds = {}
        if not self.figure3_only:
            self.load_e1()
        self.load_e2()

    def load_e1(self):
        for cohort, tiers in COHORTS.items():
            groups = defaultdict(list)
            seen = set()
            # Read only the requested accuracy and counting columns into memory.
            with self.paths[f"e1_{cohort}_metrics.tsv"].open() as stream:
                for row in csv.DictReader(stream, delimiter="\t"):
                    b = int(row["beam"])
                    identity = (row["input"], b)
                    require(identity not in seen, f"Duplicate E1 input: {identity}")
                    seen.add(identity)
                    require(row["tier"] in tiers and b in BEAMS[cohort], "Unexpected E1 bin/beam")
                    require(row["reference_mode"] == REFERENCES[cohort], "E1 reference changed")
                    require(row["analysis_status"] == "success", "Unsuccessful E1 metric row")
                    require((int(row["length"]) < 400) == (cohort == "short"), "E1 length/reference mismatch")
                    keep = {key: value(row, key) for key in METRICS}
                    for w in WINDOWS:
                        prefix = f"access_w{w}_"
                        for key in ("eligible_windows", "reference_positive_candidate_zero"):
                            keep[prefix + key] = int(row[prefix + key])
                        require(math.isnan(keep[prefix + "probability_mae_eligible"]) ==
                                (keep[prefix + "eligible_windows"] == 0), "E1 empty-domain NA mismatch")
                    groups[(row["tier"], b)].append(keep)
                    groups[("all", b)].append(keep)
            saved = {(r["bin"], int(r["beam"]), r["metric"]): r
                     for r in read_tsv(self.paths[f"e1_{cohort}_quantiles.tsv"])}
            for (tier, beam), rows in groups.items():
                deltas = np.array([r[DELTA] for r in rows])
                require(np.isfinite(deltas).all(), "Missing E1 free energy")
                counts = {
                    "positive_gt_0.01": int((deltas > 0.01).sum()),
                    "negative_lt_minus_0.01": int((deltas < -0.01).sum()),
                    "absolute_gt_1": int((np.abs(deltas) > 1).sum()),
                }
                self.thresholds[(cohort, tier, beam)] = {k: 100 * v / len(rows) for k, v in counts.items()}
                for key, count in counts.items():
                    self.audit.append(dict(panel="2D", cohort=cohort, bin=tier, beam=beam,
                                           metric=key, n=len(rows), count=count,
                                           percent=100 * count / len(rows)))
                for metric in METRICS:
                    raw = np.array([r[metric] for r in rows])
                    valid = raw[np.isfinite(raw)]
                    require(len(valid) > 0, "No evaluable RNA in E1 summary")
                    numbers = dict(zip(("min", "p10", "q25", "median", "q75", "p90", "p99", "max"),
                                       np.quantile(valid, (0, .1, .25, .5, .75, .9, .99, 1), method="linear")))
                    stored = saved[(tier, beam, metric)]
                    require(int(stored["n"]) == len(valid) and int(stored["rna"]) == len(rows),
                            "E1 quantile denominator mismatch")
                    require(int(stored["missing_or_out_of_domain_rna"]) == len(rows) - len(valid),
                            "E1 quantile NA mismatch")
                    for key in ("min", "median", "p90", "p99", "max"):
                        close(numbers[key], value(stored, key), f"E1 {cohort}/{tier}/{beam}/{metric}/{key}")
                        self.checks["e1_quantile_scalars_checked"] += 1
                    summary = dict(numbers, n=len(valid), rna=len(rows), missing=len(rows) - len(valid))
                    if metric in ACCESS:
                        w = WINDOWS[ACCESS.index(metric)]
                        summary["window_length"] = w
                        summary["eligible_windows"] = sum(r[f"access_w{w}_eligible_windows"] for r in rows)
                        summary["reference_positive_candidate_zero"] = sum(
                            r[f"access_w{w}_reference_positive_candidate_zero"] for r in rows)
                    self.q[(cohort, tier, beam, metric)] = summary
                    self.audit.append(dict(panel="2ABC", cohort=cohort, bin=tier, beam=beam,
                                           metric=metric, **summary))
                    self.checks["e1_summaries_checked"] += 1
                if cohort == "long" and beam == 400:
                    for metric in METRICS:
                        require(all(math.isnan(r[metric]) or r[metric] == 0 for r in rows),
                                "b=400 self-comparison is nonzero")
            self.checks[f"e1_{cohort}_rows"] = len(seen)
            for tier in tiers:
                ids = [len(groups[(tier, b)]) for b in BEAMS[cohort]]
                require(len(set(ids)) == 1, "E1 beam denominators differ")

    def load_e2(self):
        summary = json.loads(self.paths["e2_summary.json"].read_text())
        commits = {summary["provenance"]["runs"][label]["signature"]["commit"]
                   for label in ("rnacentral_current", "synthetic_current")}
        require(len(commits) == 1, "E2 CLI source revisions differ")
        self.e2_commit = commits.pop()
        check_revision(self.e2_commit, self.e2_revision)
        all_rows = read_tsv(self.paths["e2_metrics.tsv"])
        selected = [r for r in all_rows if r["run"] in ("rnacentral_current", "synthetic_current")]
        identities = set()
        for r in selected:
            identity = (r["run"], r["name"], r["program"])
            require(identity not in identities, f"Duplicate E2 row: {identity}")
            identities.add(identity)
            require(r["program"] in PROGRAMS, "Unknown CLI")
            require(int(r["measured_repetitions"]) == 3 and int(r["workers"]) == 10,
                    "E2 repetition/worker setting changed")
            require(all(int(r[k]) == 100 for k in ("beam", "hairpin_cap", "case3_depth_cap")),
                    "E2 beam/local caps changed")
            for metric in ("wall_seconds", "rss_gib"):
                repetitions = json.loads(r[metric + "_repetitions"])
                require(len(repetitions) == 3 and all(v > 0 for v in repetitions), "Invalid E2 repetitions")
                for statistic, measured in (("median", np.median(repetitions)),
                                             ("min", min(repetitions)), ("max", max(repetitions))):
                    close(float(measured), value(r, metric + "_" + statistic), "E2 replicate summary")
                    self.checks["e2_replicate_scalars_checked"] += 1
                self.audit.append(dict(panel="3AB" if r["run"] == "rnacentral_current" else "3C_detail",
                                       run=r["run"], name=r["name"], program=r["program"],
                                       series=r["series"], seed=r["seed"], length=int(r["n"]),
                                       metric=metric, median=value(r, metric + "_median"),
                                       min=value(r, metric + "_min"), max=value(r, metric + "_max"), n=3))
        self.rna = [r for r in selected if r["run"] == "rnacentral_current"]
        self.synthetic = [r for r in selected if r["run"] == "synthetic_current"]
        require(len(self.rna) == 21 and len(self.synthetic) == 468, "E2 cohort size changed")
        self.synthetic_summary = {}
        for program in PROGRAMS:
            for series in SERIES:
                for length in sorted({int(r["n"]) for r in self.synthetic}):
                    group = [r for r in self.synthetic if r["program"] == program and
                             r["series"] == series and int(r["n"]) == length]
                    require(len(group) == 3 and len({r["seed"] for r in group}) == 3,
                            "Synthetic group must have three seeds")
                    for metric in ("wall_seconds", "rss_gib"):
                        summary = dict(median=float(np.median([value(r, metric + "_median") for r in group])),
                                       min=min(value(r, metric + "_min") for r in group),
                                       max=max(value(r, metric + "_max") for r in group))
                        self.synthetic_summary[(program, series, length, metric)] = summary
                        self.audit.append(dict(panel="3C", run="synthetic_current", program=program,
                                               series=series, length=length, metric=metric,
                                               n=9, seeds=3, **summary))
        self.memory = [r for r in read_tsv(self.paths["e2_memory_metrics.tsv"])
                       if r["run"] == "memory_current" and r["mode"] in STAGES]
        require(len(self.memory) == 28, "Expected 7 RNA x 4 independent memory modes")
        require(len({(r["name"], r["mode"]) for r in self.memory}) == 28, "Duplicate memory row")
        for r in self.memory:
            require(r["source_revision"] == self.e2_commit, "Memory/CLI source revisions differ")
            close(value(r, "rss_gib"), value(r, "rss_kib") / 2**20, "Memory units")
            self.audit.append(dict(panel="SI_mode_memory" if self.figure3_panels == "abc" else "3D",
                                   name=r["name"], mode=r["mode"],
                                   length=int(r["n"]), metric="rss_gib", value=value(r, "rss_gib"), n=1))
        self.features = {r["name"]: r for r in read_tsv(self.paths["g_features"])}
        for r in self.rna + self.memory:
            f = self.features[r["name"]]
            for key in ("n", "sequence_sha256", "scaffold_sha256"):
                require(r[key] == f[key], f"Feature/RNA mismatch: {r['name']}/{key}")
                self.checks["e2_feature_identity_fields_checked"] += 1
        self.ratios = [r for r in read_tsv(self.paths["e2_comparisons.tsv"])
                       if r["kind"] == "rna_13652_over_11832"]
        require(len(self.ratios) == 3, "Missing RNA pair comparisons")
        for r in self.ratios:
            for metric in ("wall_seconds", "rss_gib"):
                a, b = [next(x for x in self.rna if x["name"] == r[f"name_{side}"] and
                             x["program"] == r["program"]) for side in ("a", "b")]
                close(value(b, metric + "_median") / value(a, metric + "_median"),
                      value(r, metric + "_ratio"), "RNA pair ratio")
                self.audit.append(dict(panel="3_annotation", program=r["program"],
                                       name=r["name_b"] + "/" + r["name_a"], metric=metric + "_ratio",
                                       value=value(r, metric + "_ratio")))
        for name in (self.ratios[0]["name_a"], self.ratios[0]["name_b"]):
            for metric in ("scaffold_pairs", "free_bases", "loop_free_max", "same_loop_pairs_turn"):
                self.audit.append(dict(panel="3_annotation", name=name, metric=metric,
                                       value=int(self.features[name][metric])))
        self.checks["e2_current_metric_rows"] = len(selected)
        self.checks["e2_current_stage_rows"] = len(self.memory)


def style():
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.5, "axes.titlesize": 8.5,
        "axes.labelsize": 8.5, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.fontsize": 8, "figure.titlesize": 9, "pdf.fonttype": 42,
        "ps.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": .6, "lines.linewidth": 1, "lines.markersize": 3.8,
        "grid.linewidth": .45, "grid.color": ".88", "savefig.facecolor": "white",
        "figure.constrained_layout.h_pad": .065, "figure.constrained_layout.w_pad": .065,
        "figure.constrained_layout.hspace": .06, "figure.constrained_layout.wspace": .06,
    })


def clean(ax):
    ax.set_axisbelow(True)
    ax.grid(axis="y")


def handles(labels, start=0):
    return [Line2D([], [], color=COLORS[i + start], marker=MARKERS[i + start],
                   linestyle=LINES[i + start], label=label, markerfacecolor="white")
            for i, label in enumerate(labels)]


def plot_distribution(ax, data, cohort):
    beams = BEAMS[cohort]
    width = .75 / len(beams)
    for j, tier in enumerate(COHORTS[cohort]):
        for i, beam in enumerate(beams):
            x = j + (i - (len(beams) - 1) / 2) * width
            q = data.q[(cohort, tier, beam, DELTA)]
            ax.plot([x, x], [q["min"], q["max"]], color=COLORS[i], lw=.5, alpha=.7)
            ax.plot([x, x], [q["p10"], q["p90"]], color=COLORS[i], lw=1.7)
            ax.add_patch(Rectangle((x - width * .36, q["q25"]), width * .72,
                                   q["q75"] - q["q25"], facecolor="white",
                                   edgecolor=COLORS[i], hatch=HATCHES[i], lw=.7, zorder=3))
            ax.plot(x, q["median"], marker=MARKERS[i], color=COLORS[i], markersize=3,
                    markerfacecolor="white", markeredgewidth=.8, zorder=4)
    labels = [bin_label(t) + "\nN=" + format(data.q[(cohort, t, beams[0], DELTA)]["rna"], ",")
              for t in COHORTS[cohort]]
    ax.set_xticks(range(4), labels)
    ax.set_xlim(-.55, 3.55)
    ax.set_yscale("symlog", linthresh=.001, linscale=.5)
    ax.set_ylim(-4, 400)
    ax.set_yticks([-1, 0, .01, 1, 100], ["−1", "0", "0.01", "1", "100"])
    ax.axhline(0, color=".4", lw=.5)
    ax.set_ylabel(r"$\Delta F$ (kcal/mol)")
    ax.set_xlabel("RNA length bin (nt)")
    clean(ax)


def plot_tv(ax, data, cohort):
    for i, tier in enumerate(COHORTS[cohort]):
        for k, metric in enumerate(TV):
            y = [data.q[(cohort, tier, b, metric)]["p90"] for b in BEAMS[cohort]]
            ax.plot(range(len(y)), y, color=COLORS[i], marker=MARKERS[i],
                    ls="-" if k == 0 else "--", markerfacecolor=COLORS[i] if k == 0 else "white",
                    markeredgewidth=.7, label=bin_label(tier) if k == 0 else None)
    ax.set_xticks(range(len(BEAMS[cohort])), BEAMS[cohort])
    ax.set_ylim(-.012, .42)
    ax.set_ylabel("p90 of mean TV")
    ax.set_xlabel("Beam width b")
    clean(ax)


def plot_access(ax, data, cohort):
    for i, beam in enumerate(BEAMS[cohort]):
        y = [data.q[(cohort, "all", beam, metric)]["p90"] for metric in ACCESS]
        ax.plot(range(4), y, color=COLORS[i], marker=MARKERS[i], linestyle=LINES[i],
                label=str(beam), markerfacecolor="white")
    ax.set_xticks(range(4), WINDOWS)
    maximum = max(data.q[(c, "all", b, m)]["p90"] for c in COHORTS
                  for b in BEAMS[c] for m in ACCESS)
    ax.set_ylim(-.008, math.ceil((maximum + .005) / .05) * .05)
    ax.set_xlabel("Window length (nt)")
    ax.set_ylabel("p90 of probability MAE")
    clean(ax)


def plot_threshold(ax, data, cohort, option):
    if option == "d1":
        for i, tier in enumerate(COHORTS[cohort]):
            y = [data.thresholds[(cohort, tier, b)]["positive_gt_0.01"] for b in BEAMS[cohort]]
            ax.plot(range(len(y)), y, color=COLORS[i], marker=MARKERS[i], linestyle=LINES[i],
                    markerfacecolor="white", label=bin_label(tier))
    else:
        for i, (key, label) in enumerate((
                ("positive_gt_0.01", r"$\Delta F>0.01$"),
                ("negative_lt_minus_0.01", r"$\Delta F<-0.01$"),
                ("absolute_gt_1", r"$|\Delta F|>1$"))):
            y = [data.thresholds[(cohort, "all", b)][key] for b in BEAMS[cohort]]
            ax.plot(range(len(y)), y, color=COLORS[i], marker=MARKERS[i], linestyle=LINES[i],
                    markerfacecolor="white", label=label)
    ax.set_xticks(range(len(BEAMS[cohort])), BEAMS[cohort])
    ax.set_ylim(-3, 105)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_ylabel("RNA (%)")
    ax.set_xlabel("Beam width b")
    clean(ax)


def figure2(data, width, option):
    narrow = width == 89
    fig = plt.figure(figsize=(width * MM, (508 if narrow else 274) * MM), layout="constrained")
    subfigs = fig.subfigures(4, 1, height_ratios=(1.08, 1.16, .94, 1 if option == "d2" else .92))
    titles = ("A  Free-energy difference", "B  BPP / loop-profile TV",
              "C  Window accessibility", "D1  Fraction with ΔF > 0.01" if option == "d1"
              else "D2  Signed differences and tails")
    for row, subfig in enumerate(subfigs):
        footer_needed = row in (0, 1) or (row == 3 and option == "d2")
        if narrow and row == 1:
            grid = subfig.add_gridspec(4, 1, height_ratios=(1, .20, 1, .20))
            positions = (grid[0, 0], grid[2, 0])
        else:
            ratios = (1, 1, .25) if narrow and footer_needed else (1, 1) if narrow else (
                (1, .20) if footer_needed else (1,))
            grid = subfig.add_gridspec(len(ratios), 1 if narrow else 2, height_ratios=ratios)
            positions = (grid[0, 0], grid[1, 0]) if narrow else (grid[0, 0], grid[0, 1])
        for col, cohort in enumerate(COHORTS):
            ax = subfig.add_subplot(positions[col])
            if row == 0:
                plot_distribution(ax, data, cohort)
            elif row == 1:
                plot_tv(ax, data, cohort)
            elif row == 2:
                plot_access(ax, data, cohort)
            else:
                plot_threshold(ax, data, cohort, option)
            detail = "Solid/filled: BPP; dashed/open: loop" if row == 1 else REF_LABEL[cohort]
            if row == 1:
                detail = REF_LABEL[cohort] + "\n" + detail
            if row == 2:
                detail += "\nBeam key as in A"
            if row == 3 and option == "d1":
                detail += "\nLength-bin key as in B"
            ax.set_title(titles[row] + "\n" + detail, loc="left", pad=7)
            if row == 1:
                footer = subfig.add_subplot(grid[1 + col * 2, 0] if narrow else grid[1, col])
                footer.set_axis_off()
                footer.legend(handles=handles([bin_label(t) for t in COHORTS[cohort]]),
                              loc="center", ncol=2, frameon=False, handlelength=1.4, columnspacing=1)
        if row == 0 or (row == 3 and option == "d2"):
            footer = subfig.add_subplot(grid[-1, :])
            footer.set_axis_off()
            if row == 0:
                legend_handles = handles([f"b={b}" for b in BEAMS["long"]])
            else:
                legend_handles = handles([r"$\Delta F>0.01$", r"$\Delta F<-0.01$", r"$|\Delta F|>1$"])
            footer.legend(handles=legend_handles, loc="center", ncol=3 if narrow else len(legend_handles),
                          frameon=False, columnspacing=.9, handlelength=1.6)
    fig.suptitle(("Figure 2 · Beam-width sensitivity\n" if narrow else
                  "Figure 2 · Beam-width sensitivity · ") + option.upper(), fontweight="bold")
    return fig


def plot_rna(ax, data, metric, annotated):
    for i, program in enumerate(PROGRAMS):
        rows = sorted([r for r in data.rna if r["program"] == program], key=lambda r: int(r["n"]))
        x = np.array([int(r["n"]) / 1000 for r in rows])
        y = np.array([value(r, metric + "_median") for r in rows])
        lo = np.array([value(r, metric + "_min") for r in rows])
        hi = np.array([value(r, metric + "_max") for r in rows])
        # Slightly different marker sizes reveal nearly coincident CapR/Raccess points.
        ax.errorbar(x, y, yerr=np.array([y - lo, hi - y]), fmt=MARKERS[i],
                    color=COLORS[i], markerfacecolor="white", markersize=(4.5, 5, 3.5)[i],
                    capsize=2 + i * .5, capthick=.7, elinewidth=.8, lw=0, zorder=2 + i)
    ax.set_xlabel("RNA length (kb)")
    ax.set_xticks([5, 10, 15, 20, 25, 30])
    ax.set_xlim(3, 31)
    if metric == "wall_seconds":
        ax.set_ylabel("Wall time (s)")
        ax.set_ylim(0, 1330 if annotated else 1160)
    else:
        ax.set_ylabel("Peak RSS (GiB)")
        ax.set_ylim(0, 5.3)
    if annotated:
        rows = {int(r["n"]): r for r in data.rna if r["program"] == "LinearCParty"}
        for length, xytext in ((11832, (-30, 28)), (13652, (20, 30))):
            row = rows[length]
            ax.annotate(f"{length / 1000:.1f} kb", (length / 1000, value(row, metric + "_median")),
                        xytext=xytext, textcoords="offset points", ha="center", fontsize=8,
                        arrowprops=dict(arrowstyle="-", lw=.6, color=".25"),
                        bbox=dict(facecolor="white", edgecolor="none", pad=.7, alpha=.9))
    clean(ax)


def plot_synthetic(ax, data, program, metric):
    for i, series in enumerate(SERIES):
        lengths = sorted(k[2] for k in data.synthetic_summary if k[0] == program and
                         k[1] == series and k[3] == metric)
        rows = [data.synthetic_summary[(program, series, length, metric)] for length in lengths]
        y = np.array([r["median"] for r in rows])
        lo, hi = np.array([r["min"] for r in rows]), np.array([r["max"] for r in rows])
        ax.plot(lengths, y, color=COLORS[i], marker=MARKERS[i], linestyle=LINES[i],
                markerfacecolor="white", markersize=3.5)
        ax.errorbar(lengths, y, yerr=np.array([y - lo, hi - y]), fmt="none", color=COLORS[i],
                    elinewidth=.6, capsize=1.5, alpha=.8)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks([250, 500, 1000, 2000, 4000], ["250", "500", "1,000", "2,000", "4,000"])
    ax.set_xlabel("RNA length (nt)")
    ax.set_ylabel("Wall time (s)" if metric == "wall_seconds" else "Peak RSS (GiB)")
    ax.set_ylim((.3, 700) if metric == "wall_seconds" else (.007, 3))
    clean(ax)


def plot_stages(ax, data):
    for i, stage in enumerate(STAGES):
        rows = sorted([r for r in data.memory if r["mode"] == stage], key=lambda r: int(r["n"]))
        ax.plot([int(r["n"]) / 1000 for r in rows], [value(r, "rss_gib") for r in rows],
                color=COLORS[i], marker=MARKERS[i], linestyle=LINES[i], markerfacecolor="white")
    ax.set_xlabel("RNA length (kb)")
    ax.set_ylabel("Peak RSS (GiB)")
    ax.set_xlim(3, 31)
    ax.set_ylim(0, 5.1)
    ax.set_xticks([5, 10, 15, 20, 25, 30])
    clean(ax)


def annotation_text(data, narrow):
    a, b = [data.features[data.ratios[0][f"name_{s}"]] for s in ("a", "b")]
    times = [value(r, "wall_seconds_ratio") for r in data.ratios]
    rss = [value(r, "rss_gib_ratio") for r in data.ratios]
    ratio = int(b["same_loop_pairs_turn"]) / int(a["same_loop_pairs_turn"])
    first = f"13.7 / 11.8 kb: time {min(times):.2f}–{max(times):.2f}×; RSS {min(rss):.2f}–{max(rss):.2f}×"
    second = f"Largest free loop: {a['loop_free_max']} → {b['loop_free_max']} nt; J_V: {ratio:.2f}× (different S/G)"
    if narrow:
        first = first.replace("; RSS", "\nRSS")
        second = second.replace("; J_V", "\nJ_V")
    return first + "\n" + second


def figure3(data, width, annotated):
    narrow = width == 89
    include_stages = data.figure3_panels == "abcd"
    height = ((463 if annotated else 438) if narrow else (256 if annotated else 241))
    if not include_stages:
        height -= 105 if narrow else 60
    fig = plt.figure(figsize=(width * MM, height * MM), layout="constrained")
    ratios = (1.12 if annotated else .96, 1.1) + ((.7,) if include_stages else ())
    subfigs = fig.subfigures(len(ratios), 1, height_ratios=ratios)
    top = subfigs[0]
    ratios = (1, 1, .75 if annotated else .35) if narrow else (1, .35 if annotated else .12)
    grid = top.add_gridspec(len(ratios), 1 if narrow else 2, height_ratios=ratios)
    axs = [top.add_subplot(grid[i, 0] if narrow else grid[0, i]) for i in range(2)]
    for ax, metric, title in zip(axs, ("wall_seconds", "rss_gib"),
                                  ("A  RNAcentral: wall time", "B  RNAcentral: peak RSS")):
        plot_rna(ax, data, metric, annotated)
        ax.set_title(title, loc="left", pad=7)
    footer = top.add_subplot(grid[-1, :])
    footer.set_axis_off()
    footer.legend(handles=handles(PROGRAMS), loc="upper center" if annotated else "center",
                  ncol=1 if narrow else 3, frameon=False, borderaxespad=0)
    if annotated:
        footer.text(.5, 0, annotation_text(data, narrow), ha="center", va="bottom", fontsize=8)
    middle = subfigs[1]
    ratios = (1, 1, .55) if narrow else (1, .30)
    grid = middle.add_gridspec(len(ratios), 1 if narrow else 2, height_ratios=ratios)
    axs = [middle.add_subplot(grid[i, 0] if narrow else grid[0, i]) for i in range(2)]
    for ax, metric in zip(axs, ("wall_seconds", "rss_gib")):
        plot_synthetic(ax, data, "LinearCParty", metric)
        ax.set_title("C  Synthetic controls: LinearCParty" if metric == "wall_seconds"
                     else "C  Synthetic controls: peak RSS", loc="left", pad=7)
    footer = middle.add_subplot(grid[-1, :])
    footer.set_axis_off()
    footer.legend(handles=handles(SERIES_LABELS), loc="center", ncol=2 if narrow else 3,
                  frameon=False, columnspacing=1, handlelength=2)
    if include_stages:
        stage_panel(subfigs[2], data, narrow, "D (candidate)  ")
    fig.suptitle("Figure 3 · Length and resource use" + (("\n" if narrow else " · ") + "annotated" if annotated else ""),
                 fontweight="bold")
    return fig


def stage_panel(bottom, data, narrow, prefix=""):
    grid = bottom.add_gridspec(2, 1, height_ratios=(1, .27 if narrow else .13))
    ax = bottom.add_subplot(grid[0, 0])
    plot_stages(ax, data)
    ax.set_title(prefix + ("Independent processes\nPeak RSS; one measurement per mode\nPeaks are not additive" if narrow else
                 "Independent processes: peak RSS\nOne measurement per mode; peaks are not additive"),
                 loc="left", pad=7)
    footer = bottom.add_subplot(grid[1, 0])
    footer.set_axis_off()
    footer.legend(handles=handles(STAGE_LABELS), loc="center", ncol=2 if narrow else 4,
                  frameon=False, columnspacing=1)


def stage_detail(data, width):
    fig = plt.figure(figsize=(width * MM, (115 if width == 89 else 85) * MM), layout="constrained")
    stage_panel(fig, data, width == 89)
    return fig


def synthetic_detail(data, width):
    narrow = width == 89
    fig = plt.figure(figsize=(width * MM, (388 if narrow else 171) * MM), layout="constrained")
    subfigs = fig.subfigures(4, 1, height_ratios=(1, 1, 1, .16)) if narrow else None
    if narrow:
        for subfig, program in zip(subfigs[:3], PROGRAMS):
            axs = subfig.subplots(2, 1)
            for ax, metric in zip(axs, ("wall_seconds", "rss_gib")):
                plot_synthetic(ax, data, program, metric)
                ax.set_title(program, loc="left")
        subfigs[-1].legend(handles=handles(SERIES_LABELS), loc="center", ncol=2, frameon=False)
    else:
        grid = fig.add_gridspec(3, 3, height_ratios=(1, 1, .35))
        for col, program in enumerate(PROGRAMS):
            for row, metric in enumerate(("wall_seconds", "rss_gib")):
                ax = fig.add_subplot(grid[row, col])
                plot_synthetic(ax, data, program, metric)
                ax.set_title(program, loc="left")
        footer = fig.add_subplot(grid[-1, :])
        footer.set_axis_off()
        footer.legend(handles=handles(SERIES_LABELS), loc="center", ncol=3, frameon=False)
    fig.suptitle("Figure 3C detail · Three CLIs" + ("\n" if narrow else " · ") +
                 "Seven scaffold series", fontweight="bold")
    return fig


def save(fig, name, out, width, records):
    # Resolve layout before checking visible text against the physical page.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    overflow = []
    unused_tick_labels = set()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            low, high = sorted(axis.get_view_interval())
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                if tick.get_loc() < low or tick.get_loc() > high:
                    unused_tick_labels.update((id(tick.label1), id(tick.label2)))
        for line in ax.lines:
            if line.get_transform() is ax.transData:
                x, y = np.asarray(line.get_xdata(), dtype=float), np.asarray(line.get_ydata(), dtype=float)
                for points, limits in ((x, ax.get_xlim()), (y, ax.get_ylim())):
                    low, high = sorted(limits)
                    finite = points[np.isfinite(points)]
                    require(np.all((finite >= low - 1e-9) & (finite <= high + 1e-9)),
                            f"Plotted data outside axis limits in {name}: {finite}, {limits}")
    for text in fig.findobj(match=Text):
        if not text.get_visible() or not text.get_text() or id(text) in unused_tick_labels:
            continue
        box = text.get_window_extent(renderer)
        if box.width and box.height and (box.x0 < -1 or box.y0 < -1 or
                                        box.x1 > fig.bbox.width + 1 or box.y1 > fig.bbox.height + 1):
            overflow.append(text.get_text())
    require(not overflow, f"Text outside page in {name}: {overflow}")
    for extension in ("pdf", "png"):
        path = out / f"{name}.{extension}"
        # Do not use bbox_inches='tight': preserve exactly 89/183 mm page width.
        fig.savefig(path, dpi=300, metadata={"Creator": "make_figures.py"} if extension == "pdf" else None)
        records.append(dict(file=path.name, width_mm=width, height_mm=fig.get_figheight() / MM,
                            bytes=path.stat().st_size, sha256=digest(path)))
    plt.close(fig)


def write_audit(data, out):
    keys = sorted({k for r in data.audit for k in r})
    with (out / "plot_data.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(data.audit)


def write_captions(data, out):
    sources = {k: os.path.relpath(v, out) for k, v in data.paths.items()}
    source_table = "\n".join(f"| `{k}` | [{v.name}]({sources[k]}) |" for k, v in data.paths.items())
    counts = []
    for cohort, tiers in (() if data.figure3_only else COHORTS.items()):
        for tier in tiers:
            first = data.q[(cohort, tier, BEAMS[cohort][0], DELTA)]
            counts.append(f"| {REF_LABEL[cohort]} | {bin_label(tier)} | {first['rna']:,} |")
    access_counts = []
    for cohort in (() if data.figure3_only else COHORTS):
        for w, metric in zip(WINDOWS, ACCESS):
            q = data.q[(cohort, "all", BEAMS[cohort][0], metric)]
            access_counts.append(f"| {cohort} | {w} | {q['n']:,} | {q['missing']:,} | {q['eligible_windows']:,} |")
    text = r"""# Figure 2・3 の caption 案（2026-10-16）

JA：既存の集計 TSV だけから作成した下書きである。
原稿の図番号・caption は変更していない。
Figure 2 の D1/D2 と Figure 3 の注記あり/なしを著者が比較して選ぶ。

EN: These drafts use only existing aggregate TSV files.
The manuscript figure numbers and captions have not been edited.
The author can compare D1/D2 for Figure 2 and the annotated/plain versions of Figure 3.

## ファイルと体裁 / Files and layout

- `figure2_d1_183mm.pdf/png`：長さ別の ΔF>0.01 の割合を D に置く案。
- `figure2_d2_183mm.pdf/png`：参照ごとに RNA をまとめ、正負の差と絶対差の割合を D に置く案。
- `figure3_plain_183mm.pdf/png`：注記なし、D を含む案。
- `figure3_annotated_183mm.pdf/png`：11.8/13.7 kb の点と比・G 特徴量を注記した案。
- `figure3_synthetic_all_cli_183mm.pdf/png`：C の詳細。3 CLI を同じ軸範囲で比較する。
- 各ファイルに `_89mm` 版もある / Each file also has an `_89mm` version.

JA：183 mm 版は2段幅、89 mm 版は各軸を縦に並べた1段幅のレビュー用下書きである。
89 mm 版はページ高を超えるため、採用時にはパネルを別々に配置する。
PDF のページ幅は指定値を保ち、PNG は300 dpiとした。
本文・軸名は8.5 pt、目盛り・凡例・注記は8 pt、図の表題は9 ptである。
色は Okabe–Ito 系を用い、白黒では点の形、線種、塗りとハッチでも区別する。
PDF は埋め込み TrueType フォントを使う。

EN: The 183 mm versions use two-column width; the 89 mm versions stack axes vertically for review at one-column width.
The 89 mm drafts exceed a page in height, so their panels should be placed separately if selected.
PDFs retain the stated page width; PNGs use 300 dpi.
Axis labels use 8.5 pt, ticks/legends/annotations use 8 pt, and figure headings use 9 pt.
Okabe–Ito colors are supplemented by marker shapes, line styles, fill and hatching for grayscale reading.
PDFs embed TrueType fonts.

## Figure 2 caption

JA：**Figure 2. Beam 幅に対する条件付き ensemble の解析量の感度。**
400 nt 未満の24,625 RNAは exact を参照し、400–4,381 ntの1,438 RNAは有限幅 b=400を参照した。
左右のパネル（89 mm版では上下）を参照ごとに分け、長い区分の値を exact に対する誤差とは解釈しない。
有限幅の版は d21618c5、短い区分の exact は d9804df6で、Turner 2004、37°C、outer順位、h=d=max(b,30)を用いた。
幅と局所上限を同時に変える設定全体の感度であり、pruningだけの効果や各RNAにおける単調減少は仮定しない。

EN: **Figure 2. Sensitivity of conditional-ensemble outputs to beam width.**
The reference is exact for 24,625 RNAs shorter than 400 nt and finite beam b=400 for 1,438 RNAs of 400–4,381 nt.
The left and right panels (upper and lower axes at 89 mm) separate the references; differences for long RNAs are not errors relative to exact values.
Finite-beam calculations use d21618c5 and short-RNA exact references use d9804df6, with Turner 2004, 37°C, outer ranking and h=d=max(b,30).
The plots describe sensitivity to the combined settings, which vary beam width and local caps together; they do not isolate pruning or assume monotonic changes for individual RNAs.

JA：**A**：符号付きの ΔF=F_b−F_ref（kcal/mol）を長さ区分・幅ごとに示す。
細線は最小–最大、太線はp10–p90、箱はp25–p75、記号は中央値である。
縦軸は対称対数軸で、−0.001から+0.001の間だけ線形とし、0と負値を切り捨てない。
微小な差はこの中央部に圧縮される。
各RNAを等重みにした分位点を、昇順値 x_0,…,x_(N−1) の位置 h=(N−1)p、j=floor(h)に対する (1−h+j)x_j+(h−j)x_(j+1)（端点は最小・最大）で計算した。
既存quantiles.tsvの中央値・p90・p99・最小・最大と照合し、p10・p25・p75はmetrics.tsvから同じ線形補間で追加した。

EN: **A**, Signed ΔF=F_b−F_ref (kcal/mol) by length bin and beam width.
Thin lines show the minimum–maximum, thick lines p10–p90, boxes p25–p75 and symbols the median.
The symmetric logarithmic y-axis is linear only between −0.001 and +0.001; zeros and negative values are retained.
Very small differences are compressed in this central interval.
Quantiles give each RNA equal weight: for sorted values x_0,…,x_(N−1), h=(N−1)p and j=floor(h), the quantile is (1−h+j)x_j+(h−j)x_(j+1), with minimum/maximum at the endpoints.
Medians, p90, p99 and extrema were checked against quantiles.tsv; p10, p25 and p75 were added from metrics.tsv using the same linear interpolation.

JA：**B**：RNAごとの平均BPP TVと平均loop-profile TVのp90を長さ区分別に示す。
実線・塗った記号はBPP、破線・白抜き記号はloop profileを表す。
BPPはGで自由な塩基について相手塩基と非対合を含む分布のTV、loop profileは同じ塩基の七成分のTVを平均した値である。
TVはそれぞれの分布で 1/2 Σ|p_b−p_ref| とする。
平均TVは完成構造全体の分布のTVではない。

EN: **B**, p90 across RNAs of the mean BPP TV and mean loop-profile TV, separated by length bin.
Solid lines and filled symbols show BPP; dashed lines and open symbols show loop profiles.
BPP TV compares distributions over partner bases and the unpaired state at each base free in G; loop-profile TV compares the seven context probabilities at the same bases.
TV is 1/2 Σ|p_b−p_ref| for each distribution and is averaged over these bases within each RNA.
These mean marginal TVs do not measure TV over completed structures.

JA：**C**：窓長1、4、8、16 ntについて、RNAごとのaccessibility確率MAEのp90を示す。
各RNA内で、Gの対合塩基を含まないwindowの |p_b−p_ref| を平均した後、参照ごとに全長さ区分のRNAを等重みでまとめた。
長さ区分ごとのp90の中央値やwindowを混ぜたp90ではない。
評価可能なwindowを持たないRNAはNAとして分母から除き、0に置き換えない。
確率0のwindowは評価域内なら含める。
評価RNA数とwindow数は下表に示し、窓長間で評価域が異なる点に注意する。
Opening energyはこのパネルに含めない。

EN: **C**, p90 across RNAs of accessibility probability MAE for window lengths 1, 4, 8 and 16 nt.
Within each RNA, |p_b−p_ref| is averaged over windows containing no base paired in G; RNAs from all length bins are then pooled separately for each reference with equal RNA weights.
These are neither medians of bin-level p90 values nor quantiles pooled over windows.
RNAs without an eligible window are NA and excluded from the denominator, without replacement by zero.
Zero-probability windows are included when eligible.
The table below gives RNA and window counts; the evaluation domain differs among window lengths.
Opening energies are not included in this panel.

JA：**D1案**：各長さ区分・幅で、ΔF>0.01 kcal/molのRNA数を、その区分・幅でΔFを評価できる全RNA数で割り、100を掛けた割合を示す。
0.01は表示のための閾値であり、統計的有意性や全解析に共通の許容誤差を定義しない。

EN: **D1 option**, The percentage of RNAs with ΔF>0.01 kcal/mol in each length-bin/beam group: 100 times the number above the threshold divided by all RNAs with evaluable ΔF in that group.
The threshold is a display choice, not a significance threshold or a shared tolerance for all outputs.

JA：**D2案（D1の代替）**：参照ごとに全長さ区分のRNAをまとめ、ΔF>0.01、ΔF<−0.01、|ΔF|>1 kcal/molの割合を示す。
各割合の分母はshort全24,625 RNAまたはlong全1,438 RNAである。
曲線は重なる事象を含み、足して100%にはならない。
shortでは短いRNAの件数が多いため、そのRNAが集計を強く左右する。
b=400の自己比較の0は、exactとの一致を示さない。
E1の実行時間・RSSはこの図の全案で比較に使わない。

EN: **D2 option (alternative to D1)**, Percentages with ΔF>0.01, ΔF<−0.01 and |ΔF|>1 kcal/mol after pooling all length bins separately for each reference.
Each denominator is all 24,625 short RNAs or all 1,438 long RNAs.
The events overlap, so the curves do not sum to 100%.
The short-RNA summary is strongly influenced by the more numerous shorter RNAs.
Zero differences in the b=400 self-comparison do not establish agreement with exact values.
E1 wall time and RSS are excluded from every version of this figure.

## Figure 3 caption

JA：**Figure 3. 長いRNAと合成対照における時間とpeak RSS。**
版e2cecbb4、b=h=d=100、outer順位、Turner 2004、37°Cを用いた。
cpu_40割り当て内の同時実行数上限は10プロセス、各プロセス1 thread、RLIMIT_ASは8 GiBである。
LinearCPartyはMFEと分配関数・BPP、LinearPKCapRはloop profile、LinearPKRaccessはprofileと窓長1/4/8/16のaccessibilityを要求した。
Wall timeは起動・予備計算・本計算・復元・出力を含み、scaffold構築とbuildを含まない。
Peak RSSはプロセス全体の値で、GNU timeのKiBを2^20で割ったGiBを使う。
出力の異なるCLI全体の比較から、個別解析の追加コストは分離しない。

EN: **Figure 3. Time and peak RSS for long RNAs and synthetic controls.**
Calculations use e2cecbb4, b=h=d=100, outer ranking, Turner 2004 and 37°C.
Each cpu_40 allocation permits at most 10 concurrent processes, each with one thread and an 8 GiB RLIMIT_AS limit.
LinearCParty requests MFE and the partition function/BPP; LinearPKCapR requests loop profiles; LinearPKRaccess requests profiles and accessibility at window lengths 1/4/8/16.
Wall time includes startup, preliminary and main computation, traceback and output, but excludes scaffold construction and builds.
Peak RSS is for the entire process, in GiB obtained by dividing GNU time KiB by 2^20.
The CLIs request different outputs, so their total costs do not isolate the cost of an individual analysis.

JA：**A, B**：RNAcentralの7配列（4,793–29,106 nt）の長さとwall time（A）、peak RSS（B）を示す。
点はwarm-up1回を除く測定3回の中央値、縦線は最小–最大であり、信頼区間ではない。
3回の中央値は昇順の中央の測定値である。
ほぼ重なるCapR/Raccessの点も省略しない。
各長さは異なる配列Sとscaffold Gの組であるため、点を結ぶ回帰線や漸近次数は示さない。

EN: **A, B**, Sequence length against wall time (A) and peak RSS (B) for seven RNAcentral inputs of 4,793–29,106 nt.
Points show the median of three measured repetitions after one excluded warm-up; vertical lines show the minimum–maximum, not confidence intervals.
The median of three repetitions is the middle ordered measurement.
Nearly coincident CapR/Raccess points are retained.
Each length represents a different sequence S and scaffold G; no regression line or asymptotic order is shown.

JA：**C**：合成対照の主要7系列について、250/500/1,000/2,000/4,000 ntにおけるLinearCPartyの時間とRSSを両対数軸で示す。
各系列・長さの記号は、各seedの測定3回の中央値を求めた後、その3 seedの中央値を取った値である。
縦線は3 seed×3回の計9測定の最小–最大を示し、反復とseedの両方の幅を含む。
線は同じ系列の表示点を結ぶだけで、補間値や理論上界ではない。
全156入力のうち、5長さを共通に持つ7系列×5長さ×3 seed=105入力を図示した。
残る51入力は追加の分岐深さ対照であり、この図の対象には含めない。
非emptyのGはfloor(0.15n)対、GC比率は50%である。
形を変える際に配列Sも変わるため、形や特徴量の独立な因果効果とは解釈しない。
3 CLIを同じ軸範囲で示す詳細図は `figure3_synthetic_all_cli_*` にある。

EN: **C**, LinearCParty wall time and RSS for seven main synthetic series at 250/500/1,000/2,000/4,000 nt, on logarithmic axes.
For each series/length, symbols show the median across three seeds of the within-seed medians of three repetitions.
Vertical lines span the minimum–maximum of all nine measurements (three seeds × three repetitions), combining repeat and seed variation.
Lines only connect plotted values in the same series; they are not interpolated measurements or theoretical bounds.
The plots show 105 of the 156 synthetic inputs: seven series × five shared lengths × three seeds.
The remaining 51 inputs are additional branch-depth controls outside this plotted subset.
Nonempty G has floor(0.15n) pairs and GC content is 50%.
Sequence S also changes with scaffold shape, so these comparisons do not identify independent causal effects of shape or features.
The `figure3_synthetic_all_cli_*` detail compares all three CLIs on shared axis ranges.

JA：**D案**：同じ7 RNAに対して、inside（outer）、MFE、outside+BPP、loop profileの独立したプロセスのpeak RSSを示す。
各モードはwarm-upなしの1回測定であり、反復の範囲はない。
後三者は必要なforwardと準備も含む別計算で、insideへの追加分ではない。
Peakを積み上げたり差し引いたりしてbyte内訳にはしない。
予備計算の一部を省くinside_baselineと旧版の測定は図に含めない。
この案は本文のDまたはSupplementaryへの配置を選べる。

EN: **D candidate**, Peak RSS of independent processes for inside (outer), MFE, outside+BPP and loop profiles on the same seven RNAs.
Each mode has one measurement without warm-up and therefore has no repeat range.
The latter three modes include their required forward computation and preparation; they are separate computations rather than additions to an existing inside calculation.
Peaks are neither stacked nor subtracted to infer byte components.
The inside_baseline mode that omits part of the preliminary work, and measurements of the earlier version, are not plotted.
This panel can be placed in the main figure or the Supplementary material.

JA：**注記あり案**：11,832 nt（URS00023C9544）と13,652 nt（URS000241043C）のLinearCPartyの点を矢印で示す。
注記の時間比・RSS比の範囲は3 CLIの入力別中央値の比の最小–最大であり、反復範囲や信頼区間ではない。
Gで一つのloopに直接属する自由塩基数の最大と、同じloop内の自由端点対数J_Vを既存特徴量表から添えた。
J_Vの比は13,652 ntの値を11,832 ntの値で割ったものである。
異なるS/Gの観測比較であり、時間差の原因をこの図だけで確定しない。

EN: **Annotated option**, Arrows mark the LinearCParty points for 11,832 nt (URS00023C9544) and 13,652 nt (URS000241043C).
The annotated time/RSS ratio ranges span ratios of input medians across the three CLIs, not repetitions or confidence intervals.
The largest number of free bases directly belonging to one G loop and the number J_V of free endpoint pairs in the same loop come from the existing feature table.
The J_V ratio divides the 13,652 nt value by the 11,832 nt value.
These are observations for different S/G inputs; the figure alone does not establish the cause of the time difference.

## 出典 / Sources

| Source ID | Aggregate file |
|---|---|
SOURCE_TABLE

| Panel | Rows / columns and processing |
|---|---|
| 2A | E1 short/long `metrics.tsv`: `delta_free_energy_kcal_mol`, grouped by `tier,beam`; existing `quantiles.tsv` checked; p10/p25/p75 added by linear interpolation. |
| 2B | E1 short/long `quantiles.tsv`: `base_state_tv_mean_Gfree`, `context_tv_mean_Gfree`, p90 in each length bin/beam; recomputed from metrics. |
| 2C | E1 short/long `quantiles.tsv`: `bin=all`, `access_w{1,4,8,16}_probability_mae_eligible`, p90; RNA/window counts summed from metrics; no cross-reference pooling. |
| 2D1/2D2 | E1 short/long `metrics.tsv`: signed `delta_free_energy_kcal_mol`; direct strict-threshold counts and denominators, then ×100. |
| 3A/3B | E2 `metrics.tsv`, `run=rnacentral_current`, all 3 programs: `n`, `wall_seconds_{median,min,max}`, `rss_gib_{median,min,max}`; checked against the saved repetition arrays. |
| 3C/detail | E2 `metrics.tsv`, `run=synthetic_current`, seven named `series`, all five lengths; main C uses `program=LinearCParty`, detail uses all 3; median of three seed medians, extrema over all nine repeats. |
| 3D | E2 `memory_metrics.tsv`, `run=memory_current`, `mode=inside_outer/mfe_outer/outside_bpp_outer/loop_profile_outer`, `n,rss_gib`; independent one-shot measurements. |
| 3 annotation | E2 `comparisons.tsv`, `kind=rna_13652_over_11832`, three programs, `wall_seconds_ratio,rss_gib_ratio`; G features `loop_free_max,same_loop_pairs_turn` for URS00023C9544/URS000241043C. |

JA：元データのSHA-256は `manifest.json`、実際に描いた集約値・閾値件数は `plot_data.tsv` にある。
RNA-levelの分位点は元のquantiles.tsvと照合し、E2の中央値と範囲はmetrics.tsv内の反復配列と照合した。
検証件数と出力の物理サイズは `validation.json` に記録する。

EN: `manifest.json` records source SHA-256 hashes; `plot_data.tsv` records plotted summaries and threshold counts.
RNA-level quantiles were checked against the original quantiles.tsv, and E2 medians/ranges against repetition arrays in metrics.tsv.
`validation.json` records check counts and physical output sizes.

### Figure 2A/B/D1 RNA counts

| Reference | Length bin (nt) | RNA |
|---|---:|---:|
RNA_COUNTS

### Figure 2C evaluation domain (same for every beam within each row)

| Cohort | Window (nt) | Evaluated RNA | NA RNA | Eligible windows |
|---|---:|---:|---:|---:|
ACCESS_COUNTS

JA：各参照・窓長のRNA数とwindow数が全幅で同じであることを検査した。
評価域内の参照正・有限幅0のwindow数も `plot_data.tsv` に保持する。

EN: RNA/window counts for each reference/window length were checked to be the same at every beam width.
`plot_data.tsv` also retains counts of eligible reference-positive/candidate-zero windows.

## 著者が選ぶ点 / Choices for the author

JA：Figure 2 D1は長さ別の尾部の残り方を読みやすく、A/Bの長さ区分と対応する。
D2は符号の逆転と大きな絶対差をまとめて示すが、全長さをまとめるため長さ依存を隠し、負側の割合は同じ0–100%軸では小さく見える。
推奨はD1を本文、D2をSupplementary候補とする。

EN: Figure 2 D1 shows length-dependent tails and matches the bins in A/B.
D2 summarizes negative differences and large absolute differences, but pooling hides length dependence and small negative fractions are less visible on the shared 0–100% axis.
The suggested choice is D1 for the main figure and D2 as a Supplementary candidate.

JA：Figure 3は注記あり案が11.8/13.7 kbの比較を本文と対応させやすい。
図の情報量を抑える場合は注記なし案を選び、G特徴量をSupplementaryに置く。
DはMFEの別計算のメモリを示せる一方、一回測定なので、本文に置く場合もA/Bの3回測定と区別する必要がある。
Cを3 CLIすべてに広げる場合は詳細図を使い、主図の高さ・凡例量との兼ね合いで決める。

EN: The annotated Figure 3 makes the 11.8/13.7 kb comparison easy to connect to the text.
The plain version reduces visual load and leaves G features to the Supplementary material.
D shows memory use of the separate MFE computation, but its one-shot measurements must remain distinct from the three-repetition summaries in A/B.
The all-CLI detail can replace/extend C if its extra height and legend space are acceptable.

## 再生成 / Regeneration

```bash
PYTHONDONTWRITEBYTECODE=1 python3 experiments/scripts/make_figures.py
```

JA：入力はすべてTSVで、出力先以外を変更しない。
`--e1-short`, `--e1-long`, `--e2`, `--g-features`, `--out` で場所を変更できる。
既存の入力ディレクトリを出力先として使うことは拒否する。
RNA CLI、scheduler、外部実装の実行・変更は行わない。

EN: All inputs are TSV files; only the output directory is written.
Locations can be changed with `--e1-short`, `--e1-long`, `--e2`, `--g-features` and `--out`.
An output directory overlapping input directories is rejected.
The script neither runs RNA CLIs/scheduler jobs nor executes or changes the external implementation.
"""
    text = text.replace("SOURCE_TABLE", source_table).replace("RNA_COUNTS", "\n".join(counts))
    text = text.replace("ACCESS_COUNTS", "\n".join(access_counts))
    text = text.replace("e2cecbb4", data.e2_commit[:8])
    if data.figure3_panels == "abc":
        text = text.replace("注記なし、D を含む案。", "注記なし、A–C の3パネル。")
        text = text.replace("**D案**", "**Supplementary 別図**").replace("**D candidate**", "**Separate Supplementary figure**")
        text = text.replace("この案は本文のDまたはSupplementaryへの配置を選べる。", "段階別 RSS は `figure3_mode_memory_*` に保存する。")
        text = text.replace("This panel can be placed in the main figure or the Supplementary material.",
                            "Stage RSS is saved separately as `figure3_mode_memory_*` for the Supplementary material.")
        text = text.replace("| 3D |", "| SI mode memory |")
    if data.figure3_only:
        text = "# Figure 3 と段階別 RSS（" + data.e2_commit[:8] + "）\n\n" + text.split("## Figure 3 caption\n", 1)[1].split("### Figure 2A/B/D1 RNA counts", 1)[0]
        text = "\n".join(line for line in text.splitlines() if not line.startswith("| 2")) + "\n"
        text += ("\n## 再生成 / Regeneration\n\n```bash\n"
                 f"PYTHONDONTWRITEBYTECODE=1 python3 experiments/scripts/make_figures.py --figure3-only --figure3-panels {data.figure3_panels} "
                 f"--e2 {data.paths['e2_metrics.tsv'].parent} --e2-revision {data.e2_commit[:8]} --out {out}\n```\n")
    (out / "captions.md").write_text(text, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--e1-short", type=Path, default=ROOT / "results/local/e1-bprna-short-d21618c5-20261016")
    parser.add_argument("--e1-long", type=Path, default=ROOT / "results/local/e1-bprna-long-d21618c5-20261016")
    parser.add_argument("--e2", type=Path, default=ROOT / "results/local/e2-e2cecbb4-20261016")
    parser.add_argument("--e2-revision", default="e2cecbb4", help="Required E2 commit (full hash or hex prefix)")
    parser.add_argument("--figure3-only", action="store_true", help="Read and render only E2 figures")
    parser.add_argument("--figure3-panels", choices=("abc", "abcd"), default="abcd",
                        help="Use abc for the manuscript; abcd preserves the historical layout")
    parser.add_argument("--g-features", type=Path,
                        default=ROOT / "results/local/e2-e2cecbb4-g-features-20261016/rnacentral_g_features.tsv")
    parser.add_argument("--out", type=Path, default=ROOT / "results/figures/draft_20261016")
    args = parser.parse_args()
    args.out = args.out.resolve()
    for directory in (args.e1_short, args.e1_long, args.e2, args.g_features.parent):
        directory = directory.resolve()
        require(not (args.out.is_relative_to(directory) or directory.is_relative_to(args.out)),
                "Output must not overlap an input directory")
    data = Data(args)
    for cohort in (() if args.figure3_only else COHORTS):
        for metric in ACCESS:
            counts = {(data.q[(cohort, "all", b, metric)]["n"],
                       data.q[(cohort, "all", b, metric)]["missing"],
                       data.q[(cohort, "all", b, metric)]["eligible_windows"]) for b in BEAMS[cohort]}
            require(len(counts) == 1, "E1 accessibility domains vary with beam")
    args.out.mkdir(parents=True, exist_ok=True)
    style()
    records = []
    for width in (183, 89):
        for option in (() if args.figure3_only else ("d1", "d2")):
            save(figure2(data, width, option), f"figure2_{option}_{width}mm", args.out, width, records)
        for annotated in (False, True):
            name = "annotated" if annotated else "plain"
            save(figure3(data, width, annotated), f"figure3_{name}_{width}mm", args.out, width, records)
        save(synthetic_detail(data, width), f"figure3_synthetic_all_cli_{width}mm", args.out, width, records)
        if args.figure3_panels == "abc":
            save(stage_detail(data, width), f"figure3_mode_memory_{width}mm", args.out, width, records)
        print(f"Rendered {width} mm drafts", flush=True)
    write_audit(data, args.out)
    write_captions(data, args.out)
    for key, path in data.paths.items():
        require(digest(path) == data.hashes[key], f"Source file changed: {path}")
    validation = dict(all_passed=True, checks=dict(data.checks), source_files_unchanged=True,
                      text_within_pages=True, plotted_points_within_axis_limits=True, outputs=records,
                      quantile_method="NumPy quantile(method='linear'), equal weight per evaluable RNA",
                      font_pt=dict(labels=8.5, ticks=8, legends=8, annotations=8, headings=9))
    (args.out / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
    manifest = dict(script_sha256=digest(Path(__file__)), numpy=np.__version__, matplotlib=matplotlib.__version__,
                    e2_revision=data.e2_commit, figure3_panels=args.figure3_panels,
                    sources={k: dict(path=os.path.relpath(p, args.out), sha256=data.hashes[k], bytes=p.stat().st_size)
                             for k, p in data.paths.items()}, outputs=records,
                    derived_files={name: digest(args.out / name) for name in
                                   ("captions.md", "plot_data.tsv", "validation.json")})
    (args.out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(dict(figures=len(records), checks=dict(data.checks), out=str(args.out)), indent=2))


if __name__ == "__main__":
    main()
