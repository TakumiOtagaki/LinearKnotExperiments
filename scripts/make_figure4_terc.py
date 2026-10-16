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
import os
from pathlib import Path
import platform
import re
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
DEFAULT_OUT = ROOT / "results/figures/figure4_terc_20261017"
MM = 1 / 25.4
WIDTH, HEIGHT, DPI = 183, 100, 400
REGION, WINDOW, CUTOFF = (40, 160), 4, 0.01
INK, MUTED, GRID = "#1f1f1e", "#6b6b66", "#d9d8d3"
G_COLOR, NEW_COLOR, CROSS_COLOR = "#8a8a85", "#2a78d6", "#eb6834"
CONTEXTS = ("Stem", "Hairpin", "Bulge", "Internal", "Multiloop", "Exterior", "PseudoLoop")
COLORS = dict(zip(CONTEXTS, ("#2a78d6", "#1baf7a", "#eda100", "#e87ba4",
                            "#6250d6", "#b9b8b1", "#eb6834")))
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
    bottom = np.zeros(len(positions))
    for label in labels:
        v = np.array(values[label])
        ax.bar(positions, v, bottom=bottom, width=0.92, color=COLORS[label], linewidth=0,
               rasterized=False, zorder=2)
        bottom += v
    ax.set_ylim(0, 1.01)
    ax.set_yticks((0, 0.5, 1), ("0", "0.5", "1"))
    style(ax)


def make_plot(out, rec, setting, stats):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 6.3,
        "axes.linewidth": 0.6, "axes.edgecolor": MUTED, "text.color": INK,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.major.width": 0.5,
        "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2,
        "pdf.fonttype": 42, "axes.labelsize": 6.3, "axes.labelcolor": MUTED})
    n, lo, hi = len(rec["seq"]), *REGION
    g = {tuple(p) for p in rec["nested"]}
    fig = plt.figure(figsize=(WIDTH * MM, HEIGHT * MM), dpi=DPI)
    fig.text(0.04, 0.988, "Rat telomerase RNA  |  388 nt  |  G: reference nested pairs  |  beam width 200",
             fontsize=7.5, fontweight="bold", va="top")

    def title(x, y, letter, text):
        fig.text(x, y, letter, fontsize=9, fontweight="bold", va="top")
        fig.text(x + 0.033, y - 0.001, text, fontsize=7, fontweight="bold", va="top")

    title(0.007, 0.935, "A", "Base-pairing probabilities")
    legend = fig.legend(handles=[Line2D([], [], color=c, lw=1.3, label=l) for l, c in
                        (("G", G_COLOR), ("new, no crossing", NEW_COLOR), ("crossing G", CROSS_COLOR))],
               loc="upper left", bbox_to_anchor=(0.486, 0.941), ncol=3, frameon=False,
               fontsize=6.1, handlelength=1.5, handletextpad=0.4, columnspacing=1.0, borderaxespad=0)
    fig.text(0.98, 0.935, "Darker = higher", fontsize=6.1, va="top", ha="right")
    axa = fig.add_axes((0.065, 0.749, 0.915, 0.151))
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
    for kind, color in zip(groups, (G_COLOR, NEW_COLOR, CROSS_COLOR)):
        for (i, j), p in sorted(groups[kind], key=lambda v: v[1]):
            axa.add_patch(Arc(((i + j) / 2, 0), j - i, (j - i) * 0.9,
                             theta1=0, theta2=180, color=color, lw=0.65, alpha=p))
            displayed[kind] += 1
    # Small baseline marks locate both arms of the reference pseudoknot.
    for start, end in ((63, 71), (130, 139)):
        axa.plot((start, end), (-6, -6), color=CROSS_COLOR, lw=2, solid_capstyle="butt")
    axa.set_xticks((1, 50, 100, 150, 200, 250, 300, 350, n))
    axa.set_yticks([])
    for spine in axa.spines.values():
        spine.set_visible(False)

    title(0.007, 0.694, "B", "Structural contexts")
    fig.legend(handles=[Patch(facecolor=COLORS[c], label=c) for c in CONTEXTS],
               loc="upper left", bbox_to_anchor=(0.283, 0.699), ncol=7, frameon=False,
               fontsize=6, handlelength=0.9, handletextpad=0.35, columnspacing=0.9, borderaxespad=0)
    axfull = fig.add_axes((0.065, 0.58, 0.915, 0.08))
    stack(axfull, np.arange(1, n + 1), setting["context"], CONTEXTS)
    axfull.set_xlim(0.5, n + 0.5)
    axfull.set_yticks((0, 1), ("0", "1"))
    axfull.set_xticks((1, 50, 100, 150, 200, 250, 300, 350, n))
    axfull.add_patch(Rectangle((lo - 0.5, 0), hi - lo + 1, 1.01, fill=False,
                              edgecolor=INK, lw=0.65, linestyle=(0, (3, 2)), zorder=4, clip_on=False))
    fig.text(0.065, 0.542, "40–160 enlarged", fontsize=6.3, va="top")
    peak = stats["pseudoloop"]["region_peak"]
    fig.text(0.98, 0.542, f"PseudoLoop peak: {peak['probability']:.3f} at base {peak['position']}",
             fontsize=6.3, ha="right", va="top", color=CROSS_COLOR)
    axzoom = fig.add_axes((0.065, 0.405, 0.915, 0.11))
    positions = np.arange(lo, hi + 1)
    stack(axzoom, positions, {c: setting["context"][c][lo - 1:hi] for c in CONTEXTS}, CONTEXTS)
    axzoom.set_xlim(lo - 0.5, hi + 0.5)
    axzoom.set_xticks((40, 60, 80, 100, 120, 140, 160))
    fig.text(0.018, 0.535, "Probability", fontsize=6.3, rotation=90, ha="center", va="center", color=MUTED)
    axzoom.set_xlabel("Position (nt)", labelpad=2)

    title(0.007, 0.323, "C", "4-nt window accessibility, by loop type")
    fig.text(0.065, 0.288, "All four bases unpaired together; black line = total", fontsize=6.1, va="top")
    axc = fig.add_axes((0.065, 0.105, 0.605, 0.164))
    values = {c: [setting["access"][i, WINDOW][k] for i in positions] for k, c in LOOPS.items()}
    stack(axc, positions, values, tuple(LOOPS.values()))
    axc.plot(positions, [setting["access"][i, WINDOW]["total"] for i in positions], color=INK, lw=0.65, zorder=3)
    axc.set_xlim(lo - 0.5, hi + 0.5)
    axc.set_xticks((40, 60, 80, 100, 120, 140, 160))
    axc.set_xlabel("Start of 4-nt window (nt)", labelpad=2)
    axc.set_ylabel("Probability", labelpad=3)

    title(0.711, 0.323, "D", "Reference stems in samples")
    fig.text(0.75, 0.288, "All pairs present; 1,000 samples", fontsize=6.1, va="top")
    axd = fig.add_axes((0.82, 0.126, 0.157, 0.143))
    for y, name in zip((2, 1, 0), ("outer", "inner", "both")):
        value = stats["stems"][name]["complete_sample_fraction"] * 100
        axd.barh(y, value, height=0.53, color=CROSS_COLOR, lw=0)
        axd.text(value + 3, y, f"{value:.1f}%", fontsize=6.2, va="center")
    axd.set_yticks((2, 1, 0), ("Outer (6)", "Inner (3)", "Both (9)"))
    axd.set_xlim(0, 115)
    axd.set_ylim(-0.6, 2.6)
    axd.set_xticks((0, 50, 100))
    axd.set_xlabel("Samples (%)", labelpad=1)
    axd.spines["left"].set_visible(False)
    axd.tick_params(axis="y", length=0, pad=3)
    style(axd)
    for y, key, label in ((0.053, "mfe_pairs_present", "MFE"), (0.027, "mea_pairs_present", "MEA, γ = 1")):
        outer = len(stats["stems"]["outer"][key])
        inner = len(stats["stems"]["inner"][key])
        fig.text(0.75, y, f"{label}: outer {outer}/6, inner {inner}/3 pairs", fontsize=6.1, va="top")

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
            bb_a, bb_b = a.get_window_extent(renderer), b.get_window_extent(renderer)
            if bb_a.overlaps(bb_b):
                overlaps.append((a.get_text(), b.get_text()))
    require(not overlaps, f"Text overlap: {overlaps}")
    for suffix in ("pdf", "png"):
        fig.savefig(out / f"figure4.{suffix}", dpi=DPI, facecolor="white")
    validation = dict(width_mm=WIDTH, height_mm=HEIGHT, png_dpi=DPI,
        min_font_pt=min(t.get_fontsize() for t in texts), visible_text_items=len(texts),
        text_outside_figure=outside, overlapping_text_pairs=overlaps,
        displayed_pair_counts=displayed, pair_display_cutoff=CUTOFF,
        arc_opacity="pair probability, without rescaling", reference_marks=[[63, 71], [130, 139]])
    plt.close(fig)
    return validation


def write_caption(out):
    caption = """<!-- 日本語：ラットのテロメラーゼRNAの局所構造とaccessibility。beam幅200で参照構造の入れ子の対Gを保持する。Aの弧は対合確率で、濃いほど高い。確率0.01以上の対を示す。灰色はG、青色はGと交差しない新しい対、橙色はGと交差する対。橙色の目印は参照のpseudoknotの位置を示す。 -->
**Figure 4. Structural contexts and accessibility in rat telomerase RNA.**
Calculations retain scaffold G, the reference nested pairs, at beam width 200.
**(A)** Darker arcs indicate higher pairing probabilities.
Arcs require probability ≥0.01.
Grey pairs belong to G; new pairs are blue without crossings or orange when crossing G.
Orange marks locate the reference pseudoknot.

<!-- 日本語：BのLinearPKCapRは各塩基の七つの局所構造確率を示す。PseudoLoopはpseudoknotのloop内の非対合塩基を表す。拡大区間でのPseudoLoopの最大値は123の0.840。CのLinearPKRaccessは4塩基窓のaccessibilityをloopの種類ごとに分ける。六つの同時確率の和は4塩基全部が非対合になる確率（黒線）。 -->
**(B)** LinearPKCapR gives seven context probabilities per base.
PseudoLoop describes unpaired bases in pseudoknot loops.
PseudoLoop peaks at 0.840 at base 123 in the enlarged region.
**(C)** LinearPKRaccess partitions four-base window accessibility by loop type.
The six joint probabilities sum to the probability that all four bases are unpaired (black).

<!-- 日本語：DはBoltzmann samplingの1,000標本中、各参照stemの全対を含む割合を示す。MFEは自由エネルギー最小構造、MEAは期待正確度最大構造を表す。exact計算はpseudoknot付近でほぼ同じ確率を与える。 -->
**(D)** Bars count complete reference stems among 1,000 Boltzmann samples.
MFE and MEA denote minimum free energy and maximum expected accuracy structures.
Exact calculations give nearly identical probabilities near the pseudoknot.
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
    lines = ["# Figure 4: ラットのテロメラーゼ RNA", "",
        "[PDF](figure4.pdf) · [PNG](figure4.png) · [caption](caption.md) · [provenance](provenance.json)", "",
        "bpRNA_RFAM_4317、ラット Terc（RefSeq NR_001567.1、RNAcentral URS00001ED027）、388 nt。",
        "識別情報は保存済みデータの README に従い、外部データの取得は行っていない。",
        "G は参照構造の canonical な入れ子の61対（j − i > 3）。表示値はすべて beam 幅200。", "",
        "実験 repo から次で再生成できる。RNA CLI は実行しない。", "", "```bash",
        f"python3 scripts/make_figure4_terc.py --data {relative(data)} --out {relative(out)}", "```", "",
        "依存: Python 3、NumPy、Matplotlib。キャッシュは /tmp に置き、共有領域に塩基・窓・標本ごとのファイルを作らない。", "",
        "## 図の読み方", "",
        "- A: 全長のBPP。0.01以上の対を表示し、弧のalphaを確率そのものとする。Gは灰色、新しい非交差対は青、Gと交差する対は橙。基線下の橙の目印は参照pseudoknotの両腕63–71、130–139。",
        "- B: 七つのcontextをFigure 1と同じ色で積み上げる。上段は全長、破線の40–160を下段で拡大する。",
        "- C: 40–160を開始位置とする4塩基窓。各色は全塩基が非対合かつ同じloopの種類に属する同時確率。六色の和が全体のaccessibilityで、黒線は保存済みの全体値。Bと色を共有する。座標は1-based、窓の開始位置。",
        "- D: 各stemの指定した参照対がすべて存在する標本の割合。外側6対は(63,139)–(68,134)、内側3対は(69,132)、(70,131)、(71,130)。133のbulgeを挟む。Bothは9対すべての同時出現。棒は周辺対確率の平均や積ではない。MFE/MEAには存在する参照対の個数を示す。", "",
        "## 科学的な集計", "", "| 参照stem | 対ごとの確率の範囲 | 全対を含む標本 | MFE | MEA (γ=1) |",
        "|---|---:|---:|---:|---:|"]
    for name in ("outer", "inner", "both"):
        s = stems[name]
        lines.append(f"| {name} ({len(s['pairs'])}対) | {s['pair_probability_min']:.6f}–{s['pair_probability_max']:.6f} | "
                     f"{s['complete_sample_count']}/1000 ({s['complete_sample_fraction']:.1%}) | "
                     f"{len(s['mfe_pairs_present'])}/{len(s['pairs'])} | {len(s['mea_pairs_present'])}/{len(s['pairs'])} |")
    p = stats["pseudoloop"]
    lines += ["", f"PseudoLoop: 55で{p['selected_positions']['55']:.6f}、拡大区間内の最大は123で{p['region_peak']['probability']:.6f}。",
        f"全長の最大は308で{p['global_peak']['probability']:.6f}。参照pseudoknot付近の山と全長最大を区別する。", "",
        "| 4塩基窓 | Total | Exterior | Hairpin | Bulge | Internal | Multiloop | PseudoLoop |",
        "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for i, v in windows.items():
        lines.append(f"| {i}–{v['end']} | " + " | ".join(f"{v[k]:.6f}" for k in ("total", *LOOPS)) + " |")
    lines += ["", "各loop成分の拡大区間内の最大値（括弧内は窓の開始位置）:", "",
        "| Loop | 最大同時確率 | 開始 |", "|---|---:|---:|"]
    for k, c in LOOPS.items():
        peak = stats["accessibility"]["region_component_peaks"][k]
        lines.append(f"| {c} | {peak['probability']:.6f} | {peak['position']} |")
    lines += ["", f"1,000標本は{stats['unique_sample_structures']}種類で、同じ構造の最多出現は{stats['largest_structure_count']}回。sampling seedはraw READMEによれば20261017。", "",
        "### exactとの比較", "",
        f"40–160の全contextの最大絶対差は{differences['context_region_max_abs_difference'][0]:.8f}、4塩基accessibilityは{differences['accessibility']['4']['total_region_max'][0]:.8f}。",
        "参照9対の確率は [reference_pairs.tsv](reference_pairs.tsv) にbeam200とexactを並べる。",
        f"一方、全長のcontextでは最大{differences['context_max_abs_difference'][0]:.6f}、4塩基accessibilityでは最大{differences['accessibility']['4']['total_global_max'][0]:.6f}の差がある。",
        "captionの『ほぼ同じ』はpseudoknot付近に限り、全長の完全一致を主張しない。",
        "この条件付き計算例は、無条件の構造予測精度や生物学的機能を評価する実験ではない。", "",
        "## 開発上の検証とprovenance", "",
        "入力とrec.jsonの配列・G、参照9対、全標本と代表構造でのG保持を照合する。",
        "全contextの和=1、loop別accessibilityの和=total、1塩基窓とCapRの各loop確率の一致をbeam200/exactの両方で確認する。",
        "BPPは保存済みの数表を使う。省略された低確率の対を復元しない。BPP差の集計だけは欠行を0として扱う。",
        f"図は{WIDTH} × {HEIGHT} mm、{DPI} dpi、すべての可視文字は{validation['min_font_pt']:g} pt以上。",
        "作図時に文字同士の重複と図外へのはみ出しを検査する。視覚確認の結果は論文repoのdocs/progress.mdに記録する。",
        f"captionは英語{caption_words}語（日本語コメントを除き、空白で分割）。", "",
        "[accessibility_region.tsv](accessibility_region.tsv) はCの全数表、provenance.jsonは入力・スクリプト・出力hash、環境、設定、集計、整合検査結果を含む。",
        "既存のraw READMEとRaccessヘッダが実行条件の出典。元の全CLIコマンド行とstderrはこのデータディレクトリに保存されていないため復元していない。",
        "software lockは固定buildの出典として記録する。RNA CLI再実行、原稿・参考文献・外部実装repoの編集、git commitは行わない。", ""]
    (out / "README.md").write_text("\n".join(lines))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, default=DEFAULT_DATA)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
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
    require(g == {p for p in pairs_of(rec["ref"]) - pk
                  if p[1] - p[0] > 3 and seq[p[0] - 1] + seq[p[1] - 1] in canonical},
            "Scaffold differs from canonical nested reference pairs")
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
    caption_words = write_caption(out)
    write_tables(out, settings["b200"], settings["exact"], stats["b200"])
    write_readme(out, data, stats["b200"], differences, validation, caption_words)
    sources = [data / name for name in ("in.txt", "rec.json", "README.md")]
    sources += [data / f"{s}.{ext}" for s in settings for ext in ("bpp.tsv", "capr", "racc", "samp", "lc.out", "mea")]
    sources += [data.parent / "README.md", ROOT / "manifests/software-pilot-febedd98.json"]
    lock = json.loads(sources[-1].read_text())
    provenance = dict(schema_version=1, generated_utc=datetime.now(timezone.utc).isoformat(),
        purpose="Figure 4 from existing rat Terc outputs", rna=dict(id=rec["id"], length=n,
        species="Rattus norvegicus", gene="Terc", refseq="NR_001567.1", rnacentral="URS00001ED027",
        identity_source=relative(data.parent / "README.md")),
        software=dict(source_commit=lock["source_commit"], lock=relative(sources[-1]),
                      association_source=relative(data / "README.md")),
        environment=dict(python=platform.python_version(), numpy=np.__version__, matplotlib=matplotlib.__version__),
        regeneration=dict(argv=sys.argv, script=relative(Path(__file__)), script_sha256=sha(Path(__file__))),
        sources=[dict(path=relative(p), sha256=sha(p), bytes=p.stat().st_size) for p in sources],
        design=dict(beam=200, scaffold_pairs=len(g), reference_pk_pairs=len(pk), region=list(REGION),
            sampling_seed_reported_in_raw_readme=20261017, mea_gamma=1,
            coordinates="1-based inclusive; accessibility position is window start",
            accessibility_window=WINDOW, colors=COLORS,
            complete_stem_definition="every listed reference pair occurs; other pairs are allowed",
            sample_counts="multiplicity-weighted; do not treat each unique row as one draw",
            caption_exact_agreement_scope="probabilities near reference pseudoknot, positions 40–160",
            raw_command_line_receipt="not present; settings taken from raw README and Raccess metadata"),
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
