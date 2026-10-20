#!/usr/bin/env python3
"""Render the SI PF timing/RSS comparison at 165 mm width."""
import argparse
from collections import defaultdict
import json
import os
from pathlib import Path
from statistics import median


def plot(summary_path, output):
    if output.exists():
        raise ValueError("use a new figure directory")
    os.environ.setdefault("MPLCONFIGDIR", str(Path(os.environ["TMPDIR"]) / "pf-scaling-mpl"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import NullFormatter
    from run_cli import digest, write_json
    summary = json.loads(summary_path.read_text())
    rows = summary["paired_resources"]
    if not rows:
        raise ValueError("no complete matched inputs")
    plt.rcParams.update({"font.size": 7, "axes.labelsize": 8, "legend.fontsize": 7,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "lines.linewidth": 1.1})
    fig, axes = plt.subplots(1, 2, figsize=(165 / 25.4, 66 / 25.4))
    styles = {"cparty": ("CParty", "#444444", "o"), "exact": ("LinearKnot exact", "#0072B2", "s"),
              "beam100": ("LinearKnot b=100", "#D55E00", "^")}
    for method, (label, color, marker) in styles.items():
        groups = defaultdict(list)
        for row in rows:
            if row["method"] == method:
                groups[row["length"]].append(row)
        lengths = sorted(groups)
        for ax, field, scale in ((axes[0], "elapsed_seconds", 1), (axes[1], "peak_rss_kib", 1024**2)):
            centers = [median(r[field] for r in groups[n]) / scale for n in lengths]
            lows = [min(r[field + "_min"] for r in groups[n]) / scale for n in lengths]
            highs = [max(r[field + "_max"] for r in groups[n]) / scale for n in lengths]
            # Do not bridge an unavailable intermediate length.
            for i in range(len(lengths) - 1):
                if lengths[i + 1] == 2 * lengths[i]:
                    ax.plot(lengths[i:i+2], centers[i:i+2], color=color)
            ax.errorbar(lengths, centers, yerr=[[a-b for a,b in zip(centers,lows)],
                                              [b-a for a,b in zip(centers,highs)]],
                        fmt=marker, markersize=3, capsize=2, color=color, label=label, linestyle="none")
    planned = defaultdict(int)
    counts = defaultdict(set)
    for r in summary["resources"]:
        if r["method"] == "cparty":
            planned[r["length"]] += 1
    for r in rows:
        counts[r["length"]].add(r["input"])
    denominator = "; ".join(f"{n}: {len(counts[n])}/{planned[n]}" for n in sorted(planned))
    for ax, ylabel, letter in zip(axes, ("Elapsed time (s)", "Peak RSS (GiB)"), ("A", "B")):
        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("Sequence length (nt)"); ax.set_ylabel(ylabel)
        ax.set_xticks(sorted(planned), labels=[f"{n:,}" for n in sorted(planned)])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(axis="x", labelsize=6)
        ax.spines[["top", "right"]].set_visible(False)
        ax.text(-0.16, 1.04, letter, transform=ax.transAxes, fontweight="bold", fontsize=9)
    axes[0].legend(frameon=False, loc="upper left")
    fig.text(0.5, 0.015, "Complete matched inputs (nt: used/planned): " + denominator, ha="center", fontsize=6)
    fig.tight_layout(rect=(0, .09, 1, 1), w_pad=1.5)
    output.mkdir(parents=True)
    for suffix in ("pdf", "svg", "png"):
        fig.savefig(output / ("cparty_linearknot_pf_scaling." + suffix), dpi=300)
    plt.close(fig)
    write_json(output / "figure_manifest.json", {"width_mm": 165, "height_mm": 66,
               "summary_sha256": digest(summary_path), "script_sha256": digest(Path(__file__)),
               "complete_input_counts": {n: len(counts[n]) for n in sorted(planned)},
               "planned_input_counts": dict(planned), "crossings": summary["crossings"],
               "caption": "Conditional partition function only, compatible closure policy. Points: median across input-level medians of three measurements. Bars: minimum to maximum across inputs and repetitions. The same complete inputs enter all applicable curves at each length; denominators are printed. Exact/CParty PF parity is required where exact is measured; it is unassessed where exact is excluded. Missing or numerically invalid CParty values and excluded exact lengths are not extrapolated. See resources.tsv and parity.tsv for every planned input. Timing uses concurrent single-threaded processes; inspect the saved workers and attempts."})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    plot(args.summary.resolve(), args.out.resolve())
