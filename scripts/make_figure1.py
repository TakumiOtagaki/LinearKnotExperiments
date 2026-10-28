#!/usr/bin/env python3
"""Figure 1: scaffold, seven contexts, decomposition and saved ensemble outputs.

Render the adopted exact bpRNA_RFAM_31397 example at 178 x 145 mm.
No folding program is run and no experimental values are invented.
Example, from the experiments repository::

    MPLCONFIGDIR=/tmp/figure1-mpl python3 scripts/make_figure1.py \
        --data results/local/pseudoloop-0.2.0-20261023/examples/figure_data.json \
        --context-reference /tmp/figure1-reference/test/tools/pseudoloop_context_reference.py \
        --out results/figures/figure1_fb6_20261025

The output directory receives figure1.pdf, figure1.png (300 dpi), caption.md,
README.md, structural_contexts.json and validation.json.
The old figure and saved inputs are unchanged.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import runpy
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, FancyArrowPatch, Rectangle
from matplotlib.lines import Line2D
import numpy as np
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
REFERENCE_SHA256 = "30547e8701d3c28b7e56e97384b207d405157e0c6df4b7cf71d605892740a48f"
REFERENCE_COMMIT = "f93b55ed9749dd1c198cd1f1809aece50f12f7d8"
CONTEXT_NAMES = dict(zip("SHBIMEP", CONTEXTS))


def structural_context_example(reference):
    """Classify a fictional pair map with the unmodified v0.2.0 reference."""
    digest = hashlib.sha256(reference.read_bytes()).hexdigest()
    assert digest == REFERENCE_SHA256, "Use the LinearKnot v0.2.0 geometry reference"
    classify = runpy.run_path(str(reference))["classify"]
    n = 37
    scaffold = {(2, 23), (4, 22), (6, 20), (8, 12), (25, 33), (26, 32)}
    extension = {(14, 18), (28, 36), (29, 35)}
    pairs = scaffold | extension
    assert not ({x for pair in scaffold for x in pair} &
                {x for pair in extension for x in pair})
    for group in (scaffold, extension):
        assert not any(crosses(p, q) for p in group for q in group)
    assert all(j - i >= 4 for i, j in pairs)
    labels, explanations = classify(n, pairs)
    assert set(labels) == set(CONTEXT_NAMES)
    bands = [{p for p in scaffold if any(crosses(p, q) for q in extension)},
             {q for q in extension if any(crosses(p, q) for p in scaffold)}]
    assert [len(band) for band in bands] == [2, 2]
    assert all(crosses(p, q) for p in bands[0] for q in bands[1])
    new_noncrossing = extension - bands[1]
    assert new_noncrossing
    db = ["."] * n
    g_db = db.copy()
    for i, j in scaffold:
        db[i - 1], db[j - 1] = "(", ")"
        g_db[i - 1], g_db[j - 1] = "(", ")"
    for i, j in extension:
        db[i - 1], db[j - 1] = "[", "]"
    return dict(length=n, scaffold_pairs=sorted(scaffold), extension_pairs=sorted(extension),
                completed_pairs=sorted(pairs), scaffold="".join(g_db),
                completed_structure="".join(db), labels=labels,
                contexts_by_position=[CONTEXT_NAMES[c] for c in labels],
                positions_by_context={name: [i for i, c in enumerate(labels, 1) if c == code]
                                      for code, name in CONTEXT_NAMES.items()},
                explanations_by_position=explanations, context_colors=CONTEXT_COLORS,
                reference_version="LinearKnot v0.2.0", reference_commit=REFERENCE_COMMIT,
                reference_file="test/tools/pseudoloop_context_reference.py",
                reference_sha256=digest, crossing_bands=[sorted(band) for band in bands],
                new_noncrossing_pairs=sorted(new_noncrossing),
                scope="Fictional shape illustration; no sequence, folding, energies or probabilities")

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



ROOT = Path(__file__).resolve().parents[1]
CUTOFF = .05
GREY, BLUE, ORANGE = G_COLOR, NEW_COLOR, CROSS_COLOR
COLORS = CONTEXT_COLORS
STYLES = ("-", (0, (1, 1)), (0, (4, 2)))


class Scene:
    def __init__(self, height):
        self.w, self.h = 178, height
        self.fig = plt.figure(figsize=(self.w*MM, height*MM))

    def text(self, x,y,s, **kw):
        return self.fig.text(x/self.w, y/self.h, s, va=kw.pop('va','top'), **kw)

    def panel(self,x,y,letter,title):
        self.text(x,y,letter,fontsize=9,weight='bold')
        self.text(x+6,y-.15,title,fontsize=8,weight='bold')

    def axis(self,x,y,w,h):
        return self.fig.add_axes((x/self.w,y/self.h,w/self.w,h/self.h))

    def line(self,x,y,xx,yy,**kw):
        self.fig.add_artist(Line2D([x/self.w,xx/self.w],[y/self.h,yy/self.h],
                                  transform=self.fig.transFigure,**kw))

    def rect(self,x,y,w,h,**kw):
        self.fig.add_artist(Rectangle((x/self.w,y/self.h),w/self.w,h/self.h,
                                     transform=self.fig.transFigure,**kw))

    def arrow(self,x,y,xx,yy,**kw):
        self.fig.add_artist(FancyArrowPatch((x/self.w,y/self.h),(xx/self.w,yy/self.h),
            transform=self.fig.transFigure,arrowstyle='-|>',mutation_scale=6,
            lw=.7,color=kw.pop('color',INK),**kw))

    def box(self,x,y,w,h,label):
        self.rect(x,y,w,h,fill=False,ec=MUTED,lw=.7)
        self.text(x+w/2,y+h/2,label,ha='center',va='center')


def setup():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.labelsize':7,
        'xtick.labelsize':7,'ytick.labelsize':7,'text.color':INK,'axes.linewidth':.6,
        'axes.edgecolor':MUTED,'axes.spines.top':False,'axes.spines.right':False,
        'xtick.major.size':2,'ytick.major.size':2,'xtick.major.pad':1.5,
        'ytick.major.pad':1.5,'pdf.fonttype':42})


def arc_axis(scene,x,y,w,h,lo,hi,max_height):
    ax = scene.axis(x,y,w,h)
    ax.set(xlim=(lo-.5,hi+.5),ylim=(-1,max_height))
    ax.axhline(0,color=MUTED,lw=.45)
    ax.axis('off')
    return ax


def pair_arcs(ax, pairs, g, probabilities=None):
    # Preserve probability itself as alpha; the nearby numeric key uses the
    # same colors, line styles, widths and alpha, on a white background.
    for i,j in sorted(pairs, key=lambda p: probabilities[p] if probabilities else 1):
        p = probabilities[i,j] if probabilities is not None else 1
        if probabilities is not None and p < CUTOFF: continue
        kind = 0 if (i,j) in g else 2 if any(crosses((i,j),q) for q in g) else 1
        ax.add_patch(Arc(((i+j)/2,0),j-i,(j-i)*.9,theta1=0,theta2=180,
                         lw=.9,color=(GREY,BLUE,ORANGE)[kind],ls=STYLES[kind],alpha=p))


def opacity_key(scene,x,y):
    scene.text(x,y,'Pair probability:')
    for xx,p in zip((x+26,x+39,x+52),(.1,.5,1)):
        scene.line(xx,y-1.3,xx+4,y-1.3,color=ORANGE,alpha=p,lw=.9,ls=STYLES[2])
        scene.text(xx+5,y,f'{p:g}')


def contexts(ax,lo,hi,data):
    positions = np.arange(lo,hi+1)
    edges = np.arange(lo-.5,hi+1.5)
    bottom = np.zeros(hi-lo+2)
    for c in CONTEXTS:
        v = np.array(data[c][lo-1:hi]); v = np.append(v,v[-1])
        ax.fill_between(edges,bottom,bottom+v,step='post',color=COLORS[c],lw=0)
        bottom += v
    ax.set(xlim=(lo-.5,hi+.5),ylim=(0,1.01),yticks=(0,1))


def context_key(scene,x,y):
    for row, names in enumerate((CONTEXTS[:4],CONTEXTS[4:])):
        xx=x
        for c in names:
            scene.rect(xx,y-row*4-2,2,2,fc=COLORS[c],ec='none')
            scene.text(xx+2.6,y-row*4,c)
            xx += {'Stem':15,'Hairpin':18,'Bulge':16,'Internal':20,
                   'Multiloop':23,'Exterior':21,'PseudoLoop':24}[c]


def decomposition_panel(s):
    """One BE-wing multiloop, drawn as states and grouped deduction edges.

    Each card represents a state (a set of partial structures), not a unique
    structure. Its arcs are one representative. The shared edge junction is
    not a state. In the binary expansion, the first edge has zero local
    energy; the original local factor is applied only on completion.
    """
    s.panel(63,142,'B','Partial structures, decomposition trees and pruning')
    for x,label in ((82,'CParty'),(121,'LinearKnot'),(161,'Beam pruning')):
        s.text(x,136,label,fontsize=7.5,ha='center',weight='bold')

    def arcs(x,y,w,parts=('Left','Inner','Right'),closed=False):
        # Fixed relative positions make every thumbnail the same example.
        scale=w/30
        if closed:
            s.line(x,y,x+w,y,color=MUTED,lw=.4)
            s.fig.add_artist(Arc(((x+w/2)/s.w,y/s.h),w*.97/s.w,w*.46/s.h,
                transform=s.fig.transFigure,theta1=0,theta2=180,color=GREY,lw=.85))
        for name,a,b,c,style in (('Left',3,9,BLUE,STYLES[1]),
                                 ('Inner',12,18,GREY,STYLES[0]),
                                 ('Right',21,27,BLUE,STYLES[1])):
            if name not in parts: continue
            if not closed:
                s.line(x+(a-1)*scale,y,x+(b+1)*scale,y,color=MUTED,lw=.4)
            for inset in (0,1):
                span=(b-a-2*inset)*scale
                s.fig.add_artist(Arc(((x+(a+b)/2*scale)/s.w,y/s.h),
                    span/s.w,span*.9/s.h,transform=s.fig.transFigure,
                    theta1=0,theta2=180,color=c,lw=.85,ls=style))

    def state(x,y,w,label,parts,closed=False,intermediate=False):
        s.rect(x,y,w,9,fc='white',ec=INK,lw=.65,
               ls=(0,(3,2)) if intermediate else '-')
        s.text(x+w/2,y+8.4,label,ha='center')
        # Cropped domains keep individual branches large enough to read.
        if len(parts)==1:
            full_w=(w-2)*30/8
            offset={'Left':2,'Inner':11,'Right':20}[parts[0]]
            arcs(x+1-offset*full_w/30,y+1.2,full_w,parts)
        elif intermediate:
            full_w=(w-2)*30/17
            arcs(x+1-11*full_w/30,y+1.2,full_w,parts)
        else:
            arcs(x+1,y+1.2,w-2,parts,closed)

    def join(children,junction,parent,label,label_xy):
        # All branches meet at one junction: one edge with several children.
        # A junction dot is deliberately much smaller than a state card.
        for x,y in children:
            s.line(x,y,*junction,color=INK,lw=.7)
        s.fig.add_artist(Line2D([junction[0]/s.w],[junction[1]/s.h],
            transform=s.fig.transFigure,marker='o',ms=1.5,color=INK,ls='none'))
        s.arrow(*junction,*parent)
        s.text(*label_xy,label)

    all_parts=('Left','Inner','Right')
    state(73,123,18,'Parent',all_parts,closed=True)
    for x,name in ((64,'Left'),(77,'Inner'),(90,'Right')):
        state(x,109,10.5,name,(name,))
    join([(69.25,118),(82.25,118),(95.25,118)],(82,120.5),(82,122.8),
         '+ E',(87,122.4))
    s.text(82,104,'Same multiloop',ha='center')
    arcs(68,94.5,28,closed=True)

    state(112,123,18,'Parent',all_parts,closed=True)
    state(103,109,10.5,'Left',('Left',))
    state(116,109,24,'Intermediate',('Inner','Right'),intermediate=True)
    join([(108.25,118),(128,118)],(121,120.5),(121,122.8),
         '+ E',(128,122.4))
    for x,name in ((116,'Inner'),(129.5,'Right')):
        state(x,95.5,10.5,name,(name,))
    join([(121.25,104.5),(134.75,104.5)],(128,106.5),(128,108.8),
         '+ 0',(132,108.3))

    s.text(161,132,'Keep top b (b = 3)',ha='center')
    s.text(161,128,'Intermediate states',ha='center')
    s.text(161,124.5,'Higher rank ↑',ha='center')
    for x in (151,160,169):
        for k in range(5):
            s.fig.add_artist(Line2D([x/s.w],[(109+k*3)/s.h],
                transform=s.fig.transFigure,marker='o' if k>=2 else 'x',
                ms=3.2,color=INK if k>=2 else GREY,ls='none'))
    s.line(148,113.5,173,113.5,color=MUTED,lw=.5,ls=(0,(2,2)))
    s.arrow(140.6,114,147.5,118.5,color=MUTED,ls=(0,(3,2)))
    s.arrow(148,105,174,105)
    s.text(161,103.5,'Right endpoint',ha='center')
    s.text(161,99.5,'× = discarded',ha='center')

    s.text(63,91,'Box: a state = a set of partial structures; arcs: one example.')
    s.text(63,87,'Joined lines: one step; E = local energy (added once).')


def context_schematic(s, example):
    """One colored circle per base, with the existing pair colors and styles."""
    s.text(7, 102, 'Seven contexts (schematic)', fontsize=7.5, weight='bold')
    ax = s.axis(7, 88.5, 50, 11)
    n = example['length']
    ax.set(xlim=(.3, n + .7), ylim=(-2, 11))
    ax.axis('off')
    ax.plot([1, n], [0, 0], color=MUTED, lw=.45, zorder=1)
    pair_arcs(ax, set(map(tuple, example['completed_pairs'])),
              set(map(tuple, example['scaffold_pairs'])))
    for patch in ax.patches:
        patch.set_zorder(2)
    ax.scatter(np.arange(1, n + 1), np.zeros(n), s=3**2,
               c=[COLORS[c] for c in example['contexts_by_position']],
               edgecolors=INK, linewidths=.22, zorder=3)
    s.text(7, 87, 'Base colors: context key in C')


def make_figure(d, example):
    s=Scene(145); n=len(d['record']['seq']); g=pair_set(d['scaffold'])
    gp=set(map(tuple,d['record']['pk']))
    s.panel(2,142,'A','Keep G; add G′')
    for x,label,pairs in ((7,'Scaffold G',g),(35,'Extension G′',gp)):
        s.text(x,134,label,fontsize=7.5)
        pair_arcs(arc_axis(s,x,125,23,6,1,n,26),pairs,g)
    s.text(7,122,'R = G ∪ G′',fontsize=7.5)
    pair_arcs(arc_axis(s,7,113,50,6,1,n,26),g|gp,g)
    # Only the scaffold band that crosses G′ is underlined; its enclosed
    # hairpin does not belong to that band in this new example.
    bandg={p for p in g if any(crosses(p,q) for q in gp)}
    spans=[(min(i for i,j in p),max(j for i,j in p)) for p in (bandg,gp)]
    for (i,j),y,c in zip(spans,(111,109),(GREY,ORANGE)):
        x=lambda p:7+(p-.5)*50/n
        s.line(x(i),y,x(j),y,color=c,lw=1.4)
    s.text(7,107,'Grey and orange bands cross')
    context_schematic(s, example)
    decomposition_panel(s)
    s.panel(2,82,'C','Outputs: exact example, bpRNA_RFAM_31397 (61 nt)')
    s.rect(5,1,170,76,fill=False,ec=BLUE,lw=.8)
    s.text(12,74,'Pair probabilities',fontsize=7.5,weight='bold')
    bpp={(i,j):p for i,j,p in d['bpp']}
    pair_arcs(arc_axis(s,12,61,85,10,1,n,26),bpp,g,bpp)
    opacity_key(s,12,58)
    s.text(12,52,'Seven context probabilities',fontsize=7.5,weight='bold')
    ax=s.axis(12,32,85,16); contexts(ax,1,n,d['context']); ax.set_xticks([])
    context_key(s,12,29)
    s.text(12,20.5,'Window accessibility',fontsize=7.5,weight='bold')
    ax=s.axis(12,9,85,8.5); ax.set(xlim=(.5,n+.5),ylim=(0,1.04),yticks=(0,1))
    ax.set_xticks((1,20,40,n))
    for k,(w,values) in enumerate(d['access'].items()):
        xs=sorted(map(int,values))
        ax.plot(xs,[values[str(i)] for i in xs],color=INK,ls=('-','--',':')[k],lw=.8)
        xx=61+12*k; s.line(xx,19,xx+3,19,color=INK,ls=('-','--',':')[k],lw=.8)
        s.text(xx+3.5,20.3,w+' nt')
    s.text(54.5,4.8,'Window start (nt, 1-based)',ha='center')
    s.text(111,74,'Sampled structures',fontsize=7.5,weight='bold')
    total=sum(c for c,db in d['samples'])
    for (count,db),y,label in zip(d['samples'][:2],(69,53),('Most frequent','Second most frequent')):
        s.text(111,y,f'{label}: {count:,} of {total:,}')
        pair_arcs(arc_axis(s,111,y-13,60,10,1,n,26),pair_set(db),g)
    other=total-sum(c for c,db in d['samples'][:2])
    s.text(111,37,f'Other draws: {other:,} of {total:,}')
    s.text(111,30,'γ-centroid structure (γ = 1)',fontsize=7.5,weight='bold')
    pair_arcs(arc_axis(s,111,17.5,60,10,1,n,26),pair_set(d['centroid']),g)
    for y,label,c,ls in zip((14,10,6),('Scaffold G','New; no crossing','New; crosses G'),(GREY,BLUE,ORANGE),STYLES):
        s.line(111,y-1,117,y-1,color=c,lw=.9,ls=ls); s.text(119,y,label)
    return s.fig,dict(sequence_length=n,scaffold_pairs=len(g),band_spans=spans,
        probability_display_cutoff_inclusive=CUTOFF,displayed_bpp_pairs=sum(p>=CUTOFF for p in bpp.values()),
        sample_draws=total,sample_counts_shown=[c for c,db in d['samples'][:2]],
        context_hatches='none',centroid_gamma=1,illustrative_beam_width=BEAM_DEMO,
        multiloop_rule='BE_WIP_WIP',intermediate_order=['Inner + Right','Left + Intermediate + local factor'],
        decomposition_child_counts={'CParty':[3],'LinearKnot':[2,2]},
        decomposition_local_energy={'CParty':['E'],'LinearKnot':['0','E']},
        state_thumbnails='One representative per set of partial structures',
        beam_grid_state_type='Intermediate states',
        panel_C_uses_illustrative_beam=False,
        **d['validation'])


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
    visible = [t for t in text_artists if not (
        t.axes is not None and not t.axes.axison and t not in t.axes.texts)]
    overlaps = [(t.get_text(), u.get_text()) for i, t in enumerate(visible)
                for u in visible[i+1:] if t.get_window_extent(renderer).overlaps(u.get_window_extent(renderer))]
    assert not overlaps, overlaps
    command = ["gs", "-q", "-dBATCH", "-dNOPAUSE", "-sDEVICE=txtwrite", "-dTextFormat=0",
               "-sOutputFile=-", str(pdf)]
    extracted = subprocess.run(command, check=True, text=True, capture_output=True).stdout
    sizes = [float(size) for size in re.findall(r' size="([\d.]+)"', extracted)]
    assert sizes and min(sizes) >= 7, f"Embedded text below 7 pt: {min(sizes)}"
    page = PdfReader(str(pdf)).pages[0]
    page_mm = [float(page.mediabox.width) * 25.4 / 72, float(page.mediabox.height) * 25.4 / 72]
    assert all(abs(actual - target) < 0.01 for actual, target in zip(page_mm, (WIDTH, HEIGHT)))
    with Image.open(png) as rendered:
        pixels, dpi = list(rendered.size), list(rendered.info["dpi"])
    assert all(abs(value - 300) < 0.1 for value in dpi)
    return dict(pdf_page_mm=page_mm, png_pixels=pixels, png_dpi=dpi,
                pdf_text_sizes_pt=sorted(set(sizes)), minimum_embedded_text_pt=min(sizes),
                minimum_matplotlib_text_pt=min(t.get_fontsize() for t in text_artists),
                text_outside_page=outside, text_overlaps=overlaps, gs_command=command,
                visual_inspection="Required separately: open the generated PNG.")


def pair_set(db):
    return set(pairs_from_db(db))


def read_feedback_data(source):
    """Validate the adopted saved example without computing any new RNA output."""
    d=json.loads(source.read_text())['figure1']
    assert d['record']['id']=='bpRNA_RFAM_31397'
    n=len(d['record']['seq']); g=pair_set(d['scaffold'])
    gp=set(map(tuple,d['record']['pk']))
    assert n==61 and len(g)==8 and len(d['scaffold'])==n
    assert g==set(map(tuple,d['record']['nested']))
    assert pair_set(d['record']['ref'])==g|gp
    assert not ({i for p in g for i in p}&{i for p in gp for i in p})
    for group in (g,gp):
        assert not any(crosses(p,q) for p in group for q in group)
    seq=d['record']['seq']
    assert all(seq[i-1]+seq[j-1] in {'AU','UA','CG','GC','GU','UG'} for i,j in g)
    bpp={(i,j):p for i,j,p in d['bpp']}
    assert all(0<=p<=1 for p in bpp.values())
    assert all(abs(bpp[p]-1)<1e-12 for p in g)
    assert set(d['context'])==set(CONTEXTS)
    assert all(len(d['context'][c])==n for c in CONTEXTS)
    assert max(abs(sum(d['context'][c][i] for c in CONTEXTS)-1) for i in range(n))<1e-12
    assert all(len(db)==n and g<=pair_set(db) for _,db in d['samples'])
    assert g<=pair_set(d['centroid']) and len(d['centroid'])==n
    assert d['samples']==sorted(d['samples'],key=lambda row:-row[0])
    for w,v in d['access'].items():
        assert set(map(int,v))==set(range(1,n-int(w)+2))
        assert all(0<=p<=1 for p in v.values())
    return d


def write_notes(source, out, d, report, example, reference):
    source_name = source.resolve().relative_to(ROOT).as_posix()
    out_name = out.resolve().relative_to(ROOT).as_posix()
    caption = """A下段の模式図は架空の一つの構造を示す。
各円は一つの塩基で、その色は構造文脈を示す（Cの凡例）。
橙色の非対合塩基は、交差するbandの間の隙間であるPseudoLoopに属する。

The lower sketch in A shows one fictional structure.
Each circle is one base, colored by its structural context (key in C).
Orange unpaired bases belong to PseudoLoop, the gaps between crossing bands.
"""
    (out / 'caption.md').write_text(caption)
    (out / 'structural_contexts.json').write_text(json.dumps(example, indent=2) + '\n')
    report.update(
        source_data=source_name,
        source_data_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        panel_C_data_sha256=hashlib.sha256(json.dumps(d, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        context_example_length=example['length'],
        context_reference_sha256=example['reference_sha256'],
        context_reference_commit=example['reference_commit'],
        context_reference_path_used=str(reference.resolve()),
        context_example_sha256=hashlib.sha256((out / 'structural_contexts.json').read_bytes()).hexdigest(),
        canvas_height_increase_mm=0,
        schematic_base_marker_diameter_pt=3,
        context_colors=CONTEXT_COLORS)
    (out / 'validation.json').write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
    table = '\n'.join(f"| {c} | {', '.join(map(str, example['positions_by_context'][c]))} | `{COLORS[c]}` |"
                      for c in CONTEXTS)
    (out / 'README.md').write_text(f"""# Figure 1 — seven structural contexts, 2026-10-25

[PDF](figure1.pdf)・[PNG](figure1.png)（300 dpi）。
178 × 145 mm、高さの増分は 0 mm。最小文字 {report['minimum_embedded_text_pt']:g} pt。
A 下段に一つの架空の構造を追加した。
A 上部の scaffold と extension を横に並べ、完成構造と交差 band の下に配置した。
本文からの参照は **Figure 1A (lower sketch)**（**Figure 1A 下段**）。

[PDF](figure1.pdf) and [PNG](figure1.png) (300 dpi).
The canvas remains 178 × 145 mm; the height increase is 0 mm.
The minimum embedded text size is {report['minimum_embedded_text_pt']:g} pt.
The new fictional structure is in the lower part of A.
The scaffold and extension above it are side by side, followed by their union and the band spans.
Refer to it as **Figure 1A (lower sketch)**.

## Caption addition / caption への追加案

{caption}
追加位置は既存の A の説明の最後。[caption.md](caption.md) も同じ案。
LaTeX は編集していない。高さを変えていないため、図の占有寸法は変わらない。
本文 8 ページの確認と参照文の反映は本文担当に引き渡す。

Add these sentences at the end of the existing description of A.
The same proposal is in [caption.md](caption.md).
LaTeX was not edited. The figure has the same manuscript footprint.
The manuscript editor can add the reference and check the eight-page manuscript.

## Structure and classification / 構造と分類

これは形だけの例であり、測定値や構造予測ではない。
配列・エネルギー・確率は設定していない。
灰実線は G、青点線は G と交差しない追加対、橙破線は G と交差する追加対。
各円は一つの塩基で、色は C と図 3 の既存の七色に合わせた。
G と G′ はそれぞれ交差対を持たず、対合塩基を共有しない。
H 型 pseudoknot は各二対の二つの band からなる。

This example illustrates shape; it is not a measurement or a structure prediction.
It has no sequence, energy or probability values.
Solid grey arcs show G; dotted blue arcs show new pairs that do not cross G;
dashed orange arcs show new pairs crossing G.
Each circle is one base; its color uses the existing seven-color key in C and Figure 3.
G and G′ each have no crossing pairs and share no paired bases.
The H-type pseudoknot has two bands of two pairs each.

- Length / 長さ: {example['length']} nt.
- G: `{example['scaffold']}`
- R = G ∪ G′: `{example['completed_structure']}`
- G pairs / G の対: `{example['scaffold_pairs']}`
- G′ pairs / G′ の対: `{example['extension_pairs']}`
- H-type bands / H 型の band: `{example['crossing_bands']}`
- New noncrossing pair / 交差しない追加対: `{example['new_noncrossing_pairs']}`

| Context | Positions (1-based) / 塩基位置 | Color / 色 |
|---|---|---|
{table}

分類は手で割り当てず、LinearKnot **v0.2.0** の形状分類を実行した。
`git archive v0.2.0` で参照コードと fixture だけを node-local に取り出した。
版の commit は `{REFERENCE_COMMIT}`。
`test/tools/pseudoloop_context_reference.py` の `classify(37, completed_pairs)` を使用した。
参照コードの SHA-256 は `{REFERENCE_SHA256}`。
作図スクリプトはこの hash を確認してから分類し、その出力で各円を塗る。
参照コードの 14 fixture・545 非対合窓の検証も成功した。
全塩基の分類と判定根拠は [structural_contexts.json](structural_contexts.json) にある。

The labels were computed with the unmodified shape classifier from LinearKnot **v0.2.0**.
Only the reference code and its fixtures were extracted to node-local storage with `git archive v0.2.0`.
The commit is `{REFERENCE_COMMIT}`.
The call is `classify(37, completed_pairs)` in `test/tools/pseudoloop_context_reference.py`.
Its SHA-256 is `{REFERENCE_SHA256}`.
The renderer checks this hash and uses the returned labels to color the circles.
The reference also passed all 14 fixtures and 545 fully unpaired windows.
[structural_contexts.json](structural_contexts.json) saves every label and its geometric explanation.

PseudoLoop は塩基 27、30–31、34 の三つの隙間だけに現れる。
各隙間の両端は交差する別 band の対に属する。
ヘアピン 9–11・15–17、バルジ 3、内部ループ 5・21、multiloop 7・13・19、外部 1・24・37 も含む。
band の二対が閉じる multiloop の追加例は、図を小さく保つため入れていない。
この参照分類は形の分類であり、任意の構造の folding/admissibility の証明ではない。

PseudoLoop occurs only in the three gaps at bases 27, 30–31 and 34.
Each gap has boundaries from different crossing bands.
The example also includes hairpins at 9–11 and 15–17, a bulge at 3,
an internal loop at 5 and 21, a multiloop at 7, 13 and 19, and exterior bases at 1, 24 and 37.
The optional multiloop between two pairs of one crossing band is omitted to keep the sketch small.
This reference classifies shape; it does not certify folding or model admissibility for arbitrary pair maps.

## Reproduce / 再作図

experiments/ から実行する。参照コードは node-local のコピーを指定する。
六つの生成物を保存する。RNA の再計算や外部 repo の編集は行わない。

Run from experiments/ with the extracted reference code.
The renderer writes six files. It does not rerun RNA folding or modify an external repository.

~~~sh
MPLCONFIGDIR=/tmp/figure1-fb6-mpl python3 -B scripts/make_figure1.py \\
    --data {source_name} \\
    --context-reference /tmp/figure1-reference/test/tools/pseudoloop_context_reference.py \\
    --out {out_name}
~~~

参照コードの取り出しと fixture 検証（論文 repo の root から実行し、既存の空ディレクトリを指定）:

Extract the reference and run its fixtures from the paper repository root
(use an existing empty node-local directory):

~~~sh
git -C ../LinearKnot archive v0.2.0 \\
    test/tools/pseudoloop_context_reference.py \\
    test/fixtures/pseudoloop_context_change.json | tar -x -C /tmp/figure1-reference
python3 -B /tmp/figure1-reference/test/tools/pseudoloop_context_reference.py
~~~

Dependencies: Matplotlib, NumPy, Pillow, pypdf, Ghostscript.
入力は既存の 0.2.0 集計。B と C の内容・配置・数値は保持する。
使用データ・作図コード・参照コードの hash は [validation.json](validation.json) に記録する。

The input is the adopted saved 0.2.0 data.
B and C keep their content, layout and values.
[validation.json](validation.json) records hashes for the data, renderer and classifier.

## Feedback review / 注釈の確認

`feedbacks/main_satoken_10252232.pdf` と s-1.png〜s-8.png を一枚ずつ開いて確認した。
pypdf で `/Annots` の本文注釈（Highlight / StrikeOut）と座標を取得した。
対象語句・コメント・コメントのないマークの判断根拠は [feedback_annotations.json](feedback_annotations.json) に保存する。
リンクと Popup は独立した本文注釈として数えない。
実際の markup は 32 件（Highlight 30、StrikeOut 2）。依頼文の 31 件とは一件差があるため全件を記録した。

Opened s-1.png through s-8.png individually and checked the annotated PDF.
pypdf supplied `/Annots` markup types, comments and coordinates.
[feedback_annotations.json](feedback_annotations.json) records target text and the grounds for interpreting uncommented marks.
Links and popups are not counted as independent text markup.
The PDF has 32 markup objects (30 highlights and two strikeouts); all were retained despite the request's count of 31.

- p1: Abstract の二か所の削除、CParty、時間・メモリの比較、pseudoknot のクラス、IPknot/ThreshKnot の重複、HFold の構造集合。
- p2: Spark、CapR/Raccess、time、交差の式、scaffold/with、新しい対、未定義の集合、Figure 1B の説明の細かさ。
- p3: Figure 1 全体を確認。本文 markup はない。
- p4: **“seven structural contexts” に “図”**。今回の A 下段で対応。他は TV、計算量、設定の説明に関する指摘。
- p5: GitHub bugs、centroid 検証、candidate selection。
- p6: ask、outer-pair probabilities と inner-pair probabilities。
- p7: binding/exposure の言い回し、complete-command、時間・メモリの手法比較。
- p8: 本文 markup はない。

今回の対象は p4 の文脈の図だけ。他の本文修正は本文担当に引き渡す。

This figure addresses the p4 request for a diagram of the seven contexts.
The other text comments remain with the manuscript editor.

## Validation / 検証

寸法・埋め込み文字サイズ・文字の重なり・ページ外へのはみ出しを作図時に確認する。
最終 PNG の目視結果、B/C の画素比較、export hash は validation.json に追記する。
export 先は論文 repo の figures/figure1.pdf。manifest は Figure 1 項目だけを更新する。
git commit は行わない。

The renderer checks dimensions, embedded font sizes, text intersections and page bounds.
Final PNG inspection, B/C pixel comparisons and export hashes are recorded in validation.json.
The export is paper figures/figure1.pdf; only the Figure 1 manifest entry is changed.
No Git commit is made.
""")


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data',type=Path,required=True,help='Adopted saved feedback data.json')
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--context-reference', type=Path, required=True,
                        help='Unmodified v0.2.0 pseudoloop_context_reference.py')
    args=parser.parse_args()
    d=read_feedback_data(args.data)
    setup()
    example=structural_context_example(args.context_reference)
    fig,details=make_figure(d,example)
    args.out.mkdir(parents=True,exist_ok=True)
    pdf,png=args.out/'figure1.pdf',args.out/'figure1.png'
    fig.savefig(pdf)
    fig.savefig(png,dpi=300)
    details.update(validate_figure(fig,pdf,png))
    write_notes(args.data,args.out,d,details,example,args.context_reference)
    plt.close(fig)
    print(json.dumps({'output':str(args.out),'figure_mm':details['pdf_page_mm'],
        'minimum_embedded_text_pt':details['minimum_embedded_text_pt'],
        'sample_counts_shown':details['sample_counts_shown']}))


if __name__=='__main__':
    main()
