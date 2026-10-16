#!/usr/bin/env python3
"""Draw draft E3 figures from analyze_e3_neat1.py tables, without RNA software."""
import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import tempfile

import numpy as np

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'e3-matplotlib-cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MaxNLocator

ROOT = Path(__file__).resolve().parents[1]
CONTEXTS = ('Stem', 'Hairpin', 'Bulge', 'Internal', 'Multiloop', 'Exterior', 'PseudoLoop')
# Okabe-Ito palette, augmented by gray; labels/line styles also distinguish series.
COLORS = {'LinearPartition': '#0072B2', 'emptyG': '#E69F00', 'G0.5': '#CC79A7',
          'G0.7': '#009E73', 'G0.9': '#D55E00'}
CONTEXT_COLORS = ['#777777', '#56B4E9', '#E69F00', '#0072B2', '#CC79A7', '#F0E442', '#009E73']
METHODS = list(COLORS)
STYLES = ['-', '--', ':', '-.', (0, (3, 1, 1, 1))]


def table(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def number(value):
    return float(value) if value not in ('NA', '', None) else np.nan


def smooth(a, width=51):
    a = np.asarray(a, float)
    valid = np.isfinite(a)
    numerator = np.convolve(np.where(valid, a, 0), np.ones(width), 'same')
    denominator = np.convolve(valid.astype(float), np.ones(width), 'same')
    result = np.divide(numerator, denominator, out=np.full_like(a, np.nan), where=denominator > 0)
    result[~valid] = np.nan  # Do not extend SHAPE smoothing beyond measured positions.
    return result


def save(fig, out, name):
    fig.savefig(out / (name + '.pdf'), bbox_inches='tight')
    fig.savefig(out / (name + '.png'), dpi=180, bbox_inches='tight')
    plt.close(fig)


def profile_arrays(rows):
    result = {}
    for r in rows:
        k = (r['iso'], r['method'], int(r['beam']))
        result.setdefault(k, []).append(r)
    arrays = {}
    for k, group in result.items():
        group.sort(key=lambda r: int(r['position']))
        arrays[k] = {c: np.array([number(r[c]) for r in group])
                     for c in ('position', 'q', 'shape', 'crossing_incident_mass', *CONTEXTS)}
    return arrays


def shape_figure(data, bins, out):
    fig, axes = plt.subplots(3, 1, figsize=(11, 8.4), gridspec_kw={'height_ratios': [1, 1, 1.2]})
    g = data['NEAT1_1', 'G0.7', 200]
    x = g['position']
    axes[0].plot(x, g['shape'], color='.7', lw=.45, label='SHAPE (raw)')
    axes[0].plot(x, smooth(g['shape']), color='black', lw=1.3, label='SHAPE (51-nt mean)')
    axes[0].set(ylabel='SHAPE reactivity', xlim=(1, 3756), title='a  NEAT1_1: SHAPE and predicted unpaired probability')
    axes[0].legend(loc='upper right', fontsize=8, ncol=2)
    axes[1].plot(x, g['q'], color=COLORS['G0.7'], alpha=.2, lw=.4)
    for method, style in zip(('LinearPartition', 'emptyG', 'G0.7'), STYLES):
        y = data['NEAT1_1', method, 200]['q']
        axes[1].plot(x, smooth(y), color=COLORS[method], ls=style, label=f'{method}, b=200')
    axes[1].set(ylabel='Unpaired probability q', xlabel='Nucleotide position', xlim=(1, 3756), ylim=(0, 1))
    axes[1].legend(loc='upper right', fontsize=8, ncol=3)
    for method, style in zip(METHODS, STYLES):
        selected = [r for r in bins if r['method'] == method and r['beam'] == '200'
                    and r['feature'] == 'q' and r['stratum'] == 'all' and int(r['n']) > 0]
        axes[2].plot([number(r['probability_mean']) for r in selected],
                     [number(r['shape_mean']) for r in selected], marker='o', markersize=4,
                     color=COLORS[method], ls=style, label=f'{method}, b=200')
    axes[2].set(xlabel='Unpaired probability q (bin mean)', ylabel='Mean SHAPE reactivity', xlim=(0, 1))
    axes[2].legend(fontsize=8, ncol=3)
    fig.tight_layout()
    save(fig, out, 'a_shape_profiles_bins')


def metric_figure(metrics, out):
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    for row, stat in enumerate(('spearman', 'auroc')):
        for col, stratum in enumerate(('all', 'G07_free')):
            ax = axes[row, col]
            for method, style in zip(METHODS, STYLES):
                selected = [r for r in metrics if r['method'] == method and r['feature'] == 'q' and r['stratum'] == stratum]
                selected.sort(key=lambda r: int(r['beam']))
                ax.plot([int(r['beam']) for r in selected], [number(r[stat]) for r in selected],
                        marker='o', color=COLORS[method], ls=style, label=method)
            ax.set(title='All SHAPE positions' if col == 0 else 'Positions free in G=0.7',
                   ylabel='Spearman correlation' if row == 0 else 'AUROC (high SHAPE)', xticks=[100, 200, 400])
            if row:
                ax.axhline(.5, color='.7', lw=.6, ls=':')
                ax.set_xlabel('Beam width b')
    axes[0, 0].legend(fontsize=8, loc='best')
    fig.suptitle('b  SHAPE comparison across scaffold conditions')
    fig.tight_layout()
    save(fig, out, 'b_shape_metrics')


def arc(ax, rows, limits, threshold=.01, compress=None, scaffold=True):
    groups = {'Scaffold G': [], 'New, no G crossing': [], 'New, crosses G': []}
    alphas = {c: [] for c in groups}
    for r in rows:
        i, j, p = int(r['i']), int(r['j']), number(r['probability'])
        fixed = r['is_scaffold'] == '1'
        if p < threshold or (fixed and not scaffold):
            continue
        if not limits[0] <= i < j <= limits[1]:
            continue
        category = 'Scaffold G' if fixed else ('New, crosses G' if r['crosses_G'] == '1' else 'New, no G crossing')
        if compress:
            i, j = compress(i), compress(j)
        theta = np.linspace(0, np.pi, 35)
        xx = (i + j) / 2 + (j - i) / 2 * np.cos(theta)
        yy = (j - i) / 2 * np.sin(theta)
        groups[category].append(np.column_stack([xx, yy]))
        alphas[category].append(.15 + .85 * p)
    palette = {'Scaffold G': '#777777', 'New, no G crossing': '#0072B2', 'New, crosses G': '#D55E00'}
    for c, paths in groups.items():
        if paths:
            colors = [matplotlib.colors.to_rgba(palette[c], alpha=a) for a in alphas[c]]
            ax.add_collection(LineCollection(paths, colors=colors, linewidths=.55 if c == 'Scaffold G' else .8))
    ax.axhline(0, color='.3', lw=.6)
    ax.set_yticks([])
    ax.spines[['top', 'right', 'left']].set_visible(False)
    maxheight = max((path[:, 1].max() for paths in groups.values() for path in paths), default=1)
    ax.set_ylim(0, maxheight * 1.07)
    return {c: len(paths) for c, paths in groups.items()}


def arcs_figure(pairs, longrows, out):
    fig, axes = plt.subplots(3, 1, figsize=(12, 7.5), sharex=True)
    counts = []
    for ax, beam in zip(axes, (100, 200, 400)):
        rows = [r for r in pairs if r['iso'] == 'NEAT1_1' and r['method'] == 'G0.7' and int(r['beam']) == beam]
        counts.append(dict(iso='NEAT1_1', beam=beam, **arc(ax, rows, (1, 3756))))
        ax.set(xlim=(1, 3756), ylabel=f'b={beam}')
    handles = [Line2D([0], [0], color=c, label=l) for l, c in
               [('Scaffold G', '#777777'), ('New, no G crossing', '#0072B2'), ('New, crosses G', '#D55E00')]]
    axes[0].legend(handles=handles, fontsize=8, loc='upper right', ncol=3)
    axes[-1].set_xlabel('Nucleotide position')
    fig.suptitle('c1  NEAT1_1: retained pairs with probability at least 0.01')
    fig.tight_layout()
    save(fig, out, 'c1_bpp_arcs_neat1_1')
    # Two coordinate segments; a marked gap replaces the intervening ~20 kb.
    selected = [('LinearPartition', 200), ('emptyG', 100), ('emptyG', 200),
                ('G0.7', 100), ('G0.7', 200), ('G0.7', 400)]
    fig, axes = plt.subplots(len(selected), 1, figsize=(11, 8.7), sharex=True)
    def compress(x):
        return x - 282 if x < 1000 else x - 20761 + 680
    for ax, (method, beam) in zip(axes, selected):
        rows = [r for r in pairs if r['iso'] == 'NEAT1_2' and r['method'] == method and int(r['beam']) == beam
                and (282 <= int(r['i']) <= 546 or 600 <= int(r['i']) <= 840)
                and 20761 <= int(r['j']) <= 21120]
        counts.append(dict(iso='NEAT1_2', method=method, beam=beam,
                           **arc(ax, rows, (282, 21120), threshold=0, compress=compress)))
        status = next(r for r in longrows if r['method'] == method and int(r['beam']) == beam and r['region'] == '5p_union')
        ax.set_ylabel(f'{method}\nb={beam}', rotation=0, labelpad=55, va='center', fontsize=8)
        if not rows:
            label = 'BPP checkpoint missing' if status['status'] != 'available' else f'No exported pair (cutoff {status["cutoff"]})'
            ax.text(.5, .4, label, transform=ax.transAxes, ha='center', color='.35', fontsize=9, zorder=6)
        ax.axvspan(570, 665, color='white', zorder=4)
        ax.text(617, .02, '//', transform=ax.get_xaxis_transform(), ha='center', zorder=5)
        ax.set_xlim(-15, 1060)
    axes[-1].set_xticks([0, 264, 318, 558, 680, 1039], ['282', '546', '600', '840', '20,761', '21,120'])
    axes[-1].set_xlabel('Nucleotide position (gap compressed)')
    fig.suptitle('c2  NEAT1_2: pairs joining the specified 5′ and 3′ intervals')
    fig.tight_layout()
    save(fig, out, 'c2_bpp_arcs_long_range')
    return counts


def context_figure(data, accessrows, out):
    g = data['NEAT1_1', 'G0.7', 200]
    # Deterministic exploratory region: largest mean context total variation.
    # This is a figure-selection rule, not a pre-specified biological region.
    reference = data['NEAT1_1', 'emptyG', 200]
    tv = sum(np.abs(g[c] - reference[c]) for c in CONTEXTS) / 2
    scores = np.convolve(tv, np.ones(501) / 501, 'valid')
    start = int(np.argmax(scores))
    end = start + 501
    fig, axes = plt.subplots(6, 1, figsize=(11, 11), sharex=True,
                             gridspec_kw={'height_ratios': [1, 1, 1, .7, .8, .8]})
    for ax, (method, beam) in zip(axes[:3], [('emptyG', 200), ('G0.7', 100), ('G0.7', 200)]):
        a = data['NEAT1_1', method, beam]
        ax.stackplot(a['position'][start:end], [a[c][start:end] for c in CONTEXTS],
                     colors=CONTEXT_COLORS, labels=CONTEXTS, linewidth=0)
        ax.set(ylim=(0, 1), ylabel=f'{method}\nb={beam}', yticks=[0, .5, 1])
    axes[0].legend(fontsize=7.5, ncol=7, loc='upper center', bbox_to_anchor=(.5, 1.42))
    x = g['position'][start:end]
    axes[3].plot(x, g['shape'][start:end], color='.6', lw=.7)
    axes[3].plot(x, smooth(g['shape'], 21)[start:end], color='black', lw=1)
    axes[3].set_ylabel('SHAPE')
    for method, style in zip(('emptyG', 'G0.7'), STYLES):
        a = data['NEAT1_1', method, 200]
        axes[4].plot(x, a['q'][start:end], color=COLORS[method], ls=style, lw=1, label=method)
    axes[4].set(ylabel='Unpaired probability q', ylim=(0, 1), xlim=(x[0], x[-1]))
    axes[4].legend(fontsize=8, ncol=2)
    for method, style in zip(('emptyG', 'G0.7'), STYLES):
        rows = [r for r in accessrows if r['iso'] == 'NEAT1_1' and r['method'] == method
                and r['beam'] == '200' and r['width'] == '8' and start < int(r['start']) <= end]
        axes[5].plot([int(r['start']) for r in rows], [number(r['accessibility']) for r in rows],
                     color=COLORS[method], ls=style, lw=1)
    axes[5].set(ylabel='8-nt accessibility', xlabel='Nucleotide / window start', ylim=(0, 1))
    fig.suptitle(f'd  Structural contexts: nt {start + 1}–{end}\nExploratory region: largest 501-nt mean context difference (G=0.7 versus empty G, b=200)', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, .96))
    save(fig, out, 'd_context_profiles')
    return (start + 1, end)


def accessibility_figure(accessrows, out):
    a = {}
    for r in accessrows:
        if r['method'] in ('emptyG', 'G0.7') and r['beam'] == '100':
            k = (r['iso'], r['method'], int(r['width']))
            a.setdefault(k, []).append(number(r['accessibility']))
    fig, axes = plt.subplots(2, 4, figsize=(14, 6.4), sharey=True)
    for row, iso in enumerate(('NEAT1_1', 'NEAT1_2')):
        for col, width in enumerate((1, 4, 8, 16)):
            ax = axes[row, col]
            x = np.arange(1, len(a[iso, 'G0.7', width]) + 1)
            delta = np.array(a[iso, 'G0.7', width]) - np.array(a[iso, 'emptyG', width])
            ax.plot(x, delta, color='.75', lw=.35)
            ax.plot(x, smooth(delta, 51 if row == 0 else 101), color='#0072B2', lw=1)
            ax.axhline(0, color='black', lw=.5)
            ax.set(title=f'{iso}, window={width} nt', ylim=(-1.02, 1.02), xlabel='Window start (nt)')
            ax.xaxis.set_major_locator(MaxNLocator(4))
            if col == 0:
                ax.set_ylabel('Accessibility difference\nG=0.7 − empty G')
            if row:
                for lo, hi in ((282, 546), (600, 840), (20761, 21120)):
                    ax.axvspan(lo, hi, color='#E69F00', alpha=.2)
    fig.suptitle('e  Window accessibility changes, b=100 (raw and centered moving mean)')
    fig.tight_layout()
    save(fig, out, 'e_accessibility_differences')


def sampling_figure(rows, summary, out):
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    for ax, space in zip(axes[:2], ('completed', 'extension')):
        for beam, color, style in [(100, '#0072B2', '-'), (200, '#D55E00', '--')]:
            selected = [r for r in rows if r['space'] == space and r['comparison'] == 'sample_sample' and int(r['beam']) == beam]
            x = [(number(r['left']) + number(r['right'])) / 2 for r in selected]
            n = np.array([number(r['count']) for r in selected])
            ax.step(x, n / n.sum(), where='mid', color=color, ls=style, label=f'b={beam}')
        ax.set(xlabel='Pair-set Jaccard distance', ylabel='Fraction of sample pairs',
               title='Completed structures' if space == 'completed' else 'Extensions (remove G)', xlim=(0, 1))
        ax.legend(fontsize=8)
    for beam, color, shift in [(100, '#0072B2', -.06), (200, '#D55E00', .06)]:
        for j, reference in enumerate(('MFE', 'MEA')):
            r = next(r for r in summary if r['space'] == 'extension' and r['comparison'] == 'sample_' + reference and int(r['beam']) == beam)
            axes[2].errorbar(j + shift, number(r['median']),
                yerr=[[number(r['median']) - number(r['p05'])], [number(r['p95']) - number(r['median'])]],
                color=color, marker='o', capsize=4, label=f'b={beam}' if j == 0 else None)
    axes[2].set(xticks=[0, 1], xticklabels=['MFE', 'MEA'], ylabel='Extension Jaccard distance',
                title='Sample to representative\nmedian and 5th–95th percentile', ylim=(0, 1))
    axes[2].legend(fontsize=8)
    fig.suptitle('f  NEAT1_1: 1,000 samples per beam width')
    fig.tight_layout()
    save(fig, out, 'f_sampling_diversity')


def captions(out, region):
    lo, hi = region
    (out / 'captions.md').write_text(f'''# E3 figure draft captions / 図の caption 案

All panels show predictions from saved finite-beam calculations, not experimentally determined structures.
All RNA and window coordinates below are one-based.
Color is supplemented by labels or line styles; the palette follows Okabe-Ito colors.
Scaffold constraints are imposed exactly after checking raw numerical residuals: q=0 and Stem=1 at paired scaffold bases.
For LinearPartition, q is clipped to [0,1] after subtracting the exported incident BPP sum from one; raw-sum sensitivity is recorded in the analysis tables.
全図は保存済みの有限 beam 計算による予測を示す。
RNA と window の座標は1起点とする。
色に加えてラベルまたは線種で条件を識別する。
生値の残差を検査した後、scaffold の固定塩基は q=0、Stem=1 の厳密な条件を代入した。
LinearPartition の q は保存 BPP の incident sum を1から引いて [0,1] に補正し、生の和に対する感度は集計結果に残した。

## a_shape_profiles_bins

JA：NEAT1_1 の SHAPE 反応性と非対合確率を同じ座標に示す。
上段は SHAPE の生値と中心51 nt の移動平均、中段は b=200 の LinearPartition、空の G、G=0.7 の非対合確率の移動平均であり、G=0.7 の生値を薄く重ねた。
下段は五条件の非対合確率を幅0.1の区間に分けた SHAPE 平均であり、空の区間は表示しない。
SHAPE の位置1を転写産物の nt 1と仮定し、測定のある3,595塩基を用いた。
LinearPartition の非対合確率は閾値1e-4の BPP 出力から求め、LinearCParty は窓長1の Raccess を用いた。

EN: SHAPE reactivity and predicted unpaired probability along NEAT1_1.
The top panel shows raw SHAPE values and a centered 51-nt moving mean.
The middle panel shows moving means of unpaired probabilities for LinearPartition, empty G, and G=0.7 at b=200, with raw G=0.7 values in the background.
The bottom panel shows mean SHAPE reactivity in probability bins of width 0.1 for five conditions; empty bins are omitted.
We assume that SHAPE position 1 corresponds to transcript nucleotide 1 and use all 3,595 measured positions.
Unpaired probabilities come from length-1 Raccess for LinearCParty and exported BPPs with cutoff 1e-4 for LinearPartition.

## b_shape_metrics

JA：NEAT1_1 の SHAPE に対する非対合確率の Spearman 相関と AUROC を、scaffold と beam 幅ごとに比べる。
左列は3,595塩基、右列はそのうち G=0.7 で自由な共通の位置集合を示す。
AUROC の正例は全測定値の75百分位以上の SHAPE とし、全条件で同じ閾値を用いた。
線は測定点を結ぶもので、幅に対する単調性を仮定しない。

EN: Spearman correlation and AUROC of unpaired probability against NEAT1_1 SHAPE reactivity across scaffolds and beam widths.
The left column uses all 3,595 measured positions; the right column uses the common subset free in G=0.7.
Positive labels for AUROC are SHAPE values at or above the global 75th percentile, with the same threshold in every condition.
Lines connect evaluated settings without assuming monotonic beam-width sensitivity.

## c1_bpp_arcs_neat1_1

JA：G=0.7 の NEAT1_1 全体の BPP を b=100、200、400 で示す。
灰色は scaffold の対、青は G と交差しない新しい対、朱色は G の少なくとも一対と交差する新しい対であり、濃さを確率に対応させた。
読みやすさのため確率0.01以上の対だけを描いた。
BPP の保存形式は i、j、確率の3列であり、交差は座標から判定したもので square/round の出力列ではない。

EN: BPP arcs across NEAT1_1 for G=0.7 at b=100, 200, and 400.
Gray arcs denote scaffold pairs, blue arcs denote new pairs that do not cross G, and vermilion arcs denote new pairs crossing at least one scaffold pair.
Opacity increases with probability; only pairs with probability at least 0.01 are drawn.
The saved BPP format contains only i, j, and probability, so crossing is classified geometrically rather than from square/round output columns.

## c2_bpp_arcs_long_range

JA：NEAT1_2 の nt 282–546 または600–840と nt 20,761–21,120 の間にある保存 BPP を示す。
二つの表示領域の間の座標を圧縮し、切れ目を記した。
線がない条件は記載の閾値以上の対が出力されなかったことを意味し、真の確率が0であるとはしない。
未完了の条件は欠測と表示した。
この図は報告された区間での計算結果を示し、実験で確認された個別の塩基対を示さない。

EN: Exported NEAT1_2 BPPs joining nt 282–546 or 600–840 to nt 20,761–21,120.
The intervening sequence is compressed and the coordinate gap is marked.
An empty panel means that no pair reached the stated export cutoff, which does not establish zero pairing probability.
Unfinished conditions are labeled missing.
These panels show predictions in the specified intervals, not experimentally established individual base pairs.

## d_context_profiles

JA：NEAT1_1 の nt {lo}–{hi} における七つの context 確率、SHAPE、非対合確率、窓長8の accessibility を対応づける。
七つの context は各位置で和が1となる。
表示区間は b=200 の G=0.7 と空の G の context 分布の total variation の501 nt 平均が最大となる区間として探索的に選んだ。
SHAPE の黒線は中心21 nt の移動平均である。
この区間選択は予測を見やすくするためのもので、生物学的な domain の指定ではない。

EN: Seven structural-context probabilities, SHAPE reactivity, unpaired probability, and 8-nt window accessibility at NEAT1_1 nt {lo}–{hi}.
The seven context probabilities sum to one at each position.
This exploratory interval maximizes the 501-nt mean total variation between the context distributions for G=0.7 and empty G at b=200.
The black SHAPE curve is a centered 21-nt moving mean.
The interval is selected to display model output and is not assigned a biological domain.

## e_accessibility_differences

JA：b=100 における G=0.7 と空の G の window accessibility の差を、窓長1、4、8、16で示す。
Accessibility は window の全塩基が同時に非対合となる確率である。
灰色は生の差、青は NEAT1_1 で中心51 nt、NEAT1_2 で中心101 nt の移動平均を示す。
橙色の帯は指定された NEAT1_2 の長距離相互作用の区間である。
G の変更は条件とする scaffold と ensemble を同時に変えるため、この差を pseudoknot だけの効果とは解釈しない。

EN: Differences in window accessibility between G=0.7 and empty G at b=100 for window lengths 1, 4, 8, and 16.
Accessibility is the joint probability that all bases in a window are unpaired.
Gray curves show raw differences; blue curves show centered moving means of 51 nt for NEAT1_1 and 101 nt for NEAT1_2.
Orange bands mark the specified NEAT1_2 long-range intervals.
Changing G changes both the scaffold condition and the allowed ensemble, so this contrast does not isolate the effect of pseudoknots alone.

## f_sampling_diversity

JA：G=0.7 の NEAT1_1 から b=100 と200で各1,000構造を抽出した。
左と中央は全499,500標本対についての塩基対集合の Jaccard 距離であり、完成構造と G を除いた extension を分けて示す。
右は標本と MFE・MEA 構造の extension の距離の中央値および5–95百分位を示す。
乱数 seed は20261016、MEA の γ は1である。
MFE の有限 beam 探索集合は分配関数の保持 ensemble と一致するとは限らない。

EN: Diversity of 1,000 NEAT1_1 structures sampled with G=0.7 at each of b=100 and b=200.
The left and middle panels show pair-set Jaccard distances for all 499,500 unordered sample pairs, using completed structures and extensions after removing G, respectively.
The right panel shows the median and 5th–95th percentiles of extension distances from samples to the MFE and MEA structures.
The sampling seed is 20261016 and MEA uses gamma=1.
The finite-beam MFE search set need not equal the retained partition-function ensemble.
''')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=ROOT / 'results/local/e3-neat1-20261016')
    parser.add_argument('--output', type=Path, default=ROOT / 'results/figures/draft_e3_20261016')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'savefig.facecolor': 'white'})
    data = profile_arrays(table(args.input / 'position_profiles.tsv.gz'))
    accessrows = table(args.input / 'accessibility.tsv.gz')
    shape_figure(data, table(args.input / 'shape_probability_bins.tsv'), args.output)
    metric_figure(table(args.input / 'shape_metrics.tsv'), args.output)
    counts = arcs_figure(table(args.input / 'bpp_pairs.tsv.gz'), table(args.input / 'long_range.tsv'), args.output)
    region = context_figure(data, accessrows, args.output)
    accessibility_figure(accessrows, args.output)
    sampling_figure(table(args.input / 'sampling_histograms.tsv'), table(args.input / 'sampling_summary.tsv'), args.output)
    captions(args.output, region)
    (args.output / 'figure_metadata.json').write_text(json.dumps(dict(
        context_region=region, arc_threshold_neat1_1=.01, arc_counts=counts,
        smoothing='centered; shorter support at sequence ends; skip missing SHAPE',
        analysis_provenance_sha256=__import__('hashlib').sha256((args.input / 'provenance.json').read_bytes()).hexdigest(),
        plot_script_sha256=__import__('hashlib').sha256(Path(__file__).read_bytes()).hexdigest()), indent=2) + '\n')
    print('Saved seven PDF/PNG figure drafts and captions in', args.output)


if __name__ == '__main__':
    main()
