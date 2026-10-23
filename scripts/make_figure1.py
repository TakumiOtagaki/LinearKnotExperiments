#!/usr/bin/env python3
"""Figure 1: scaffold, decomposition and outputs from one conditional ensemble.

Render the adopted exact bpRNA_RFAM_31397 example at 178 x 145 mm.
No folding program is run and no experimental values are invented.
Example, from the experiments repository::

    MPLCONFIGDIR=/tmp/figure1-mpl python3 scripts/make_figure1.py \
        --data results/local/pseudoloop-0.2.0-20261023/examples/figure_data.json \
        --out results/figures/figure1_fb4_20261023

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


def make_figure(d):
    s=Scene(145); n=len(d['record']['seq']); g=pair_set(d['scaffold'])
    gp=set(map(tuple,d['record']['pk']))
    s.panel(2,142,'A','Keep G; add G′')
    for y,label,pairs in ((134,'Scaffold G',g),(120,'Extension G′',gp),(106,'R = G ∪ G′',g|gp)):
        s.text(7,y,label,fontsize=7.5)
        pair_arcs(arc_axis(s,7,y-11,50,8,1,n,26),pairs,g)
    # Only the scaffold band that crosses G′ is underlined; its enclosed
    # hairpin does not belong to that band in this new example.
    bandg={p for p in g if any(crosses(p,q) for q in gp)}
    spans=[(min(i for i,j in p),max(j for i,j in p)) for p in (bandg,gp)]
    for (i,j),y,c in zip(spans,(92,89),(GREY,ORANGE)):
        x=lambda p:7+(p-.5)*50/n
        s.line(x(i),y,x(j),y,color=c,lw=1.4)
    s.text(7,86.5,'Grey and orange bands cross')
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
        s.text(111,y,f'{label}: {count:,} / {total:,}')
        pair_arcs(arc_axis(s,111,y-13,60,10,1,n,26),pair_set(db),g)
    other=total-sum(c for c,db in d['samples'][:2])
    s.text(111,37,f'Other draws: {other:,} / {total:,}')
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


def write_notes(source,out,d,report):
    source_name=source.resolve().relative_to(ROOT).as_posix()
    out_name=out.resolve().relative_to(ROOT).as_posix()
    total=report['sample_draws']; a,b=report['sample_counts_shown']; other=total-a-b
    caption=f"""**Figure 1. Structure prediction and conditional ensemble analysis with LinearKnot.**

Aでは、与えたscaffold Gにextension G′を加えて完成構造Rを作る。
GとG′は対合塩基を共有せず、それぞれは交差対を持たない。
灰色と橙色の下線は、互いに交差する入れ子対の群（band）の範囲を示す。
density-2の制限では、交差でつながったband群のうち一つの位置を覆うbandは最大二つである。

**(A)** Adding extension G′ to the supplied scaffold G gives the completed structure R.
G and G′ share no paired bases, and each has no crossing pairs.
The grey and orange underlines mark the spans of two groups of nested pairs, called bands, that cross each other.
The density-2 restriction allows at most two bands in a connected crossing group to span one position.

Bは、二つのscaffold対の間に左右の枝を一つずつ加えたmultiloopを示す。
multiloopは、枝と呼ぶ閉じた部分構造をつなぐloopである。
各枠は状態、すなわち部分構造の集まりを表し、弧はその一例を示す。
黒い点で合流する線は、子を親につなぐ分解の一段階を表す。
CPartyは左の枝・内側のscaffold部分・右の枝を一度につなぐ。
LinearKnotは内側と右の枝を途中の状態にまとめ、次に左の枝を加える。
Eは局所エネルギーを表す。
最初の段階では0、次の段階でEを加える。
beam pruningは、部分構造の最後の塩基である右端ごとに同じ種類の状態を順位づけ、上位b個を残す（模式図はb={BEAM_DEMO}）。
途中の状態も枝刈りの対象となる。
丸は保持、×は破棄を示す。
Cは枝刈りを行わない計算例である。

**(B)** The example has one new branch on each side of a multiloop between two scaffold pairs.
A multiloop connects closed substructures, called branches.
Each box is a state: a set of partial structures.
The arcs show one example from that set.
Lines that meet at a black dot form one decomposition step, joining the children to their parent.
CParty joins the left branch, inner scaffold part and right branch at once.
LinearKnot first joins the inner part and right branch into an intermediate state, then adds the left branch.
E denotes the local energy term.
The first step adds zero; the second adds E.
Beam pruning ranks states of the same type at each right endpoint, the last base in the partial structure, and keeps the top b (b={BEAM_DEMO} in this sketch).
It also prunes intermediate states.
Circles remain; crosses are discarded.
Panel C uses no pruning.

CはbpRNA_RFAM_31397（61 nt）のexact計算を示す。
Gは、pseudoknotを持つ参照構造の入れ子部分にある8対のAU、GC、GU対からなる。
左上の対確率図だけは、0.05以上の確率を弧の濃さで示す。
右の二つの標本図とγ-centroid図は、それぞれ一つの完成構造を示す。
灰色実線はG、青点線はGと交差しない追加対、橙破線はGと交差する追加対を示す。

**(C)** Panel C shows exact calculations for bpRNA_RFAM_31397 (61 nt).
G contains eight AU, GC, or GU pairs from the nested part of the reference structure, which also contains a pseudoknot.
Only the upper-left pair-probability diagram uses arc opacity to show probability (pairs with probability at least 0.05).
The sample diagrams and the γ-centroid diagram each show one completed structure.
Solid grey arcs show G; dotted blue arcs show new pairs that do not cross G; dashed orange arcs show new pairs crossing G.

七つのcontext確率の和は各塩基で1となる。
PseudoLoopは、交差するbandの間の非対合の隙間を示す。
accessibilityは窓内の全塩基が同時に非対合となる確率である。
実線・破線・点線は1・4・8 ntの窓を示す。
accessibilityの横軸は窓の先頭の塩基の位置を1始まりで示す。
例えば、8 ntの曲線の位置20は塩基20–27の窓を表す。
標本図は{total:,}回のBoltzmann標本化で最多と2番目に多い構造（{a:,}回と{b:,}回）を示す。
残る{other:,}回の標本は図示しない。
代表構造はγ=1のγ-centroid構造である。

The seven context probabilities sum to one at each base.
PseudoLoop denotes unpaired gaps between crossing bands.
Accessibility is the probability that all bases in a window are unpaired together.
Solid, dashed, and dotted curves show 1-, 4-, and 8-nt windows.
The accessibility x-axis gives the 1-based position of the first base in each window.
For example, position 20 on the 8-nt curve represents bases 20–27.
The sample diagrams show the two most frequent structures in {total:,} Boltzmann draws ({a:,} and {b:,} draws).
The other {other:,} draws are not shown.
The representative structure is the γ-centroid structure at γ=1.
"""
    (out/'caption.md').write_text(caption)
    report.update(source_data=source_name,source_data_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        panel_C_data_sha256=hashlib.sha256(json.dumps(d,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'validation.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    (out/'README.md').write_text(f"""# Figure 1 — states, partial structures and decomposition trees, 2026-10-23

[figure1.pdf](figure1.pdf): {WIDTH} × {HEIGHT} mm.
[figure1.png](figure1.png): 300 dpi preview.
Minimum embedded text: {report['minimum_embedded_text_pt']:g} pt.
The bilingual caption proposal is below and in [caption.md](caption.md).

## Reproduce

From experiments/:

~~~sh
MPLCONFIGDIR=/tmp/figure1-fb4-mpl python3 -B scripts/make_figure1.py \\
    --data {source_name} \\
    --out {out_name}
~~~

Dependencies: Matplotlib, NumPy, Pillow, pypdf, Ghostscript.
Five output files; no RNA folding or experimental recomputation.
The input is the adopted LinearKnot 0.2.0 drawing data, including updated contexts.
Source, panel-C data and renderer hashes are in validation.json.

## Scientific drawing

B shows one multiloop between two scaffold pairs, with a new branch on each side.
Grey arcs are scaffold G; dotted blue arcs are added pairs that do not cross G.
The inner grey boundary delimits this band's calculation; its enclosed region is omitted.
Each outlined node is a state representing a set of partial structures.
Its thumbnail is one example, not the complete set.
Lines meeting at a small black dot are one edge with several children (a hyperedge).
The dot is an edge junction, not another state.
The arrows point from children toward the parent.
CParty has one three-child edge; LinearKnot has two two-child edges.
The intermediate contains the inner part and right branch, without the outer closing pair.
Its first edge adds zero local energy; the completion edge adds the original local term E once.
The two parent thumbnails and the larger lower-left multiloop have the same pairs.
The intermediate's dashed box points to a beam of intermediate states, for one state type.
At each right endpoint (last base of a partial structure), the sketch keeps {BEAM_DEMO} top-ranked states.
Black circles remain; grey crosses are discarded. These ranks are illustrative, not measured.
Coordinates, scaffold anchors, unpaired-base factors and admissibility details are omitted.

A and all C content keep their adopted geometry, labels and saved values.
The misleading B-to-C arrow is removed: C is an exact example without pruning.
C retains sample frequencies {a:,}, {b:,}, and {other:,} other draws out of {total:,}.
No inference about these exact outputs is made from the illustrative beam grid.

## Implementation evidence (not figure labels)

The example is BE_WIP_WIP: left WIP × inner BE × right WIP.
Ordinary V/VM multiloop recursion is already binary in CParty; this drawing makes no general claim that CParty multiloops all have three children.
LinearKnot first combines inner BE × right WIP, then left WIP and the local outer factor.
The split coordinates and scaffold anchors are retained when grouping the products.
Contracting the intermediate restores the original edge and its weight before pruning.
The diagram does not claim that other local energy terms vanish.

Read-only evidence:
- results/figures/figure1_fb2_20261023/README.md.
- LinearKnot v0.2.0 documents/model.md, “Derivation correspondence” (BE wings).
- LinearKnot v0.2.0 src/scfg2/deductions/expand_basic.hh (BE_WIP_WIP).
- LinearKnot v0.2.0 src/scfg2/beam/be_frame_cache.hh: the intermediate is a product of the two child factors; the local factor enters on completion.
- LinearKnot src/scfg2/beam/be_right_tail_selection.hh for the separate intermediate selection.

## What was taken from the earlier versions

Viewed the earlier exported PDF at paper commit 867ccd7 and the adopted 0.2.0 PNG individually.
Read the earlier renderer history (experiments 2076dd8 / 350ba2f) and the current renderer.
The earlier boxes and branching lines made the decomposition tree clear; those are restored.
The current multiloop arcs, the inner-plus-right order, and the pair colors are retained.
Every tree node now has a structural thumbnail; the intermediate is explicitly a state.
A shared edge junction groups all three CParty children into one step.
Both binary steps and their local-energy assignments are visible in LinearKnot.

## Author feedback review

Opened m-1.png through m-9.png individually and checked the marked sentences against the PDF annotation coordinates.
Author StrikeOut annotations were distinguished from red latexdiff deletions.
The subtype-free Annots entries on pages 2, 4 and 5 resolve to XRef streams, not extra author comments.
Marks outside Figure 1 are recorded here for scope awareness and left to the manuscript editor.

- p1: highlight on the Lyngsø and Pedersen (2000) citation in the NP-hard sentence.
- p2: highlights on “No scaffold pair”, “may connect a base inside a completed pseudoknot region”, and “internal-loop length limits”; the margin asks what the undefined pseudoknot region means and asks for concrete limit values (the latter is visible in the image, not a typed annotation).
- p3: highlight on Intermediate and comment asking for the state–structure–tree relationship; comment contrasting the beam sketch with exact C; author StrikeOut on “numeric key;” and “with opaque arcs”. The new B addresses the tree comment, the B-to-C arrow is removed, and the caption proposal omits the struck phrases.
- p4: highlights on “calculation with a small beam” / “for ranking” and on the outside-calculation sentence; margin comments request the reason for the pilot beam and a short outside explanation with an example. An image-only red comment next to the TV formula asks for wording such as “TV is defined by”.
- p5: StrikeOut on the end of the context paragraph (other loops retaining labels, and the band multiloop example); margin suggests deleting from “other loops” or “For example”; highlights and comments on sampling “its parts”, computational-cost “dependence on the local limits” (image-only note: is “local limits” understandable?), and BPP “alternatives”.
- p6: StrikeOut on the sentence excluding RNAs without eligible windows; highlight and comment on “four window lengths”, asking to list them.
- p7: no author annotation or handwritten comment found.
- p8: highlight on the partition-function ratio in the accessibility explanation, with a comment questioning its placement and “no-crossing calculation”; StrikeOut and comment on the Discussion probability-output / centroid sentence, suggesting sampling or subopt instead; bottom comment asks for the default scaffold construction in Methods.
- p9: two adjacent StrikeOut annotations remove the statements that retaining every scaffold pair excludes replacement structures and that a crossing pair can have low probability when other structures have lower energies.

## Validation and handoff

Page size, embedded text size, text intersections and saved inputs are checked during rendering.
PNG visual inspection and unchanged A/C comparisons are recorded in validation.json.
Export: paper figures/figure1.pdf; only its export-manifest entry is updated.
LaTeX, bibliography and external implementation repositories are unchanged.
The canvas remains 178 × 145 mm, so the figure's manuscript footprint is unchanged.
No Git commit is made.

## Caption proposal (Japanese / English)

{caption}""")


def main():
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data',type=Path,required=True,help='Adopted saved feedback data.json')
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    d=read_feedback_data(args.data)
    setup()
    fig,details=make_figure(d)
    args.out.mkdir(parents=True,exist_ok=True)
    pdf,png=args.out/'figure1.pdf',args.out/'figure1.png'
    fig.savefig(pdf)
    fig.savefig(png,dpi=300)
    details.update(validate_figure(fig,pdf,png))
    write_notes(args.data,args.out,d,details)
    plt.close(fig)
    print(json.dumps({'output':str(args.out),'figure_mm':details['pdf_page_mm'],
        'minimum_embedded_text_pt':details['minimum_embedded_text_pt'],
        'sample_counts_shown':details['sample_counts_shown']}))


if __name__=='__main__':
    main()
