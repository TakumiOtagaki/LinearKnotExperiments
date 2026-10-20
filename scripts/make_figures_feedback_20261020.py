#!/usr/bin/env python3
"""Render author-feedback figures from the packed, validated local summaries.

Run prepare_feedback_figure_data.py first. No RNA program is run here.
All outputs use new directories ending in _fb; no manuscript is edited.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

os.environ.setdefault('MPLCONFIGDIR', '/tmp/figure-feedback-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Arc, Rectangle, FancyArrowPatch
from matplotlib.text import Text
import numpy as np
from pypdf import PdfReader

import make_figure1 as f1
import make_figure4_terc as terc
import make_figures_review as review

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'results/local/figure-feedback-20261020'
OUT = ROOT / 'results/figures'
MM = 1/25.4
CUTOFF = .05
CONTEXTS, COLORS = f1.CONTEXTS, f1.CONTEXT_COLORS
INK, GREY, BLUE, ORANGE = f1.INK, f1.G_COLOR, f1.NEW_COLOR, f1.CROSS_COLOR
STYLES = ('-', (0,(1,1)), (0,(4,2)))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False)+'\n')


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
        self.rect(x,y,w,h,fill=False,ec=f1.MUTED,lw=.7)
        self.text(x+w/2,y+h/2,label,ha='center',va='center')


def setup():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.labelsize':7,
        'xtick.labelsize':7,'ytick.labelsize':7,'text.color':INK,'axes.linewidth':.6,
        'axes.edgecolor':f1.MUTED,'axes.spines.top':False,'axes.spines.right':False,
        'xtick.major.size':2,'ytick.major.size':2,'xtick.major.pad':1.5,
        'ytick.major.pad':1.5,'pdf.fonttype':42})


def arc_axis(scene,x,y,w,h,lo,hi,max_height):
    ax = scene.axis(x,y,w,h)
    ax.set(xlim=(lo-.5,hi+.5),ylim=(-1,max_height))
    ax.axhline(0,color=f1.MUTED,lw=.45)
    ax.axis('off')
    return ax


def pair_arcs(ax, pairs, g, probabilities=None):
    # Preserve probability itself as alpha; the nearby numeric key uses the
    # same colors, line styles, widths and alpha, on a white background.
    for i,j in sorted(pairs, key=lambda p: probabilities[p] if probabilities else 1):
        p = probabilities[i,j] if probabilities is not None else 1
        if probabilities is not None and p < CUTOFF: continue
        kind = 0 if (i,j) in g else 2 if any(terc.crosses((i,j),q) for q in g) else 1
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


def figure1(d):
    s=Scene(145); n=len(d['record']['seq']); g=terc.pairs_of(d['scaffold'])
    gp=set(map(tuple,d['record']['pk']))
    s.panel(2,142,'A','Keep G; add G′')
    for y,label,pairs in ((134,'Scaffold G',g),(120,'Extension G′',gp),(106,'R = G ∪ G′',g|gp)):
        s.text(7,y,label,fontsize=7.5)
        pair_arcs(arc_axis(s,7,y-11,50,8,1,n,26),pairs,g)
    # Only the scaffold band that crosses G′ is underlined; its enclosed
    # hairpin does not belong to that band in this new example.
    bandg={p for p in g if any(terc.crosses(p,q) for q in gp)}
    spans=[(min(i for i,j in p),max(j for i,j in p)) for p in (bandg,gp)]
    for (i,j),y,c in zip(spans,(92,89),(GREY,ORANGE)):
        x=lambda p:7+(p-.5)*50/n
        s.line(x(i),y,x(j),y,color=c,lw=1.4)
    probe=(max(i for i,j in spans)+min(j for i,j in spans))/2
    xp=7+(probe-.5)*50/n
    s.line(xp,88,xp,104,color=INK,lw=.7,ls=':')
    s.text(7,86.5,'Two band spans overlap here')
    s.panel(63,142,'B','Smaller steps; beam pruning')
    s.text(79,134,'CParty',ha='center',weight='bold')
    s.text(119,134,'Smaller steps',ha='center',weight='bold')
    s.text(160,134,'Keep top b',ha='center',weight='bold')
    s.box(69,122,20,5,'region')
    for x in (65,75,85):
        s.box(x,112,8,5,'part'); s.arrow(79,122,x+4,117)
    s.text(79,109,'Combine several\nparts in one step',ha='center')
    s.arrow(95,124,102,124)
    s.box(110,123,21,5,'region'); s.box(104,114,18,5,'partial'); s.box(126,114,8,5,'part')
    s.arrow(113,119,118,123); s.arrow(130,119,125,123)
    for x in (108,119):
        s.box(x-1,108,2,2,''); s.arrow(x,110,113,114)
    s.text(120,106,'Combine two, then one',ha='center')
    s.arrow(136,124,142,124)
    for x in (153,161,169):
        for k in range(5):
            s.fig.add_artist(Line2D([x/s.w],[(112+k*3.2)/s.h],transform=s.fig.transFigure,
                marker='o' if k>=2 else 'x',ms=3,color=BLUE if k>=2 else GREY,ls='none'))
    s.arrow(147,111,147,128); s.text(144.5,120,'Higher rank',rotation=90,ha='center',va='center')
    s.arrow(151,109,174,109); s.text(162,107,'Right endpoint',ha='center')
    s.text(70,99,'Sum weights → partition function')
    s.text(70,93.5,'Minimum energy → MFE')
    s.arrow(170,99,170,78,color=BLUE)
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
    s.text(54.5,4.8,'Position / start of window',ha='center')
    s.text(111,74,'Boltzmann samples',fontsize=7.5,weight='bold')
    total=sum(c for c,db in d['samples'])
    for (count,db),y,label in zip(d['samples'][:2],(69,53),('Most frequent','Second most frequent')):
        s.text(111,y,f'{label}: {count:,} / {total:,}')
        pair_arcs(arc_axis(s,111,y-13,60,10,1,n,26),terc.pairs_of(db),g)
    other=total-sum(c for c,db in d['samples'][:2])
    s.text(111,37,f'Other draws: {other:,} / {total:,}')
    s.text(111,30,'γ-centroid structure (γ = 1)',fontsize=7.5,weight='bold')
    pair_arcs(arc_axis(s,111,17.5,60,10,1,n,26),terc.pairs_of(d['centroid']),g)
    for y,label,c,ls in zip((14,10,6),('Scaffold G','New; no crossing','New; crosses G'),(GREY,BLUE,ORANGE),STYLES):
        s.line(111,y-1,117,y-1,color=c,lw=.9,ls=ls); s.text(119,y,label)
    return s.fig,dict(sequence_length=n,scaffold_pairs=len(g),band_spans=spans,
        probability_display_cutoff_inclusive=CUTOFF,displayed_bpp_pairs=sum(p>=CUTOFF for p in bpp.values()),
        sample_draws=total,sample_counts_shown=[c for c,db in d['samples'][:2]],
        context_hatches='none',centroid_gamma=1,**d['validation'])


def terc_figure(d):
    s=Scene(104); n=388; v=d['settings']['b200']; stats=v['stats']
    g=terc.pairs_of(d['scaffold']); bpp={(i,j):p for i,j,p in v['bpp']}
    s.text(8,103,'Rat Terc  |  388 nt  |  G: 61 reference pairs  |  beam width 200',fontsize=8,weight='bold')
    s.panel(1,97,'A','Pair probabilities')
    for x,label,c,ls in zip((68,88,132),('G','New, no crossing','Crossing G'),(GREY,BLUE,ORANGE),STYLES):
        s.line(x,95.5,x+5,95.5,color=c,lw=.9,ls=ls); s.text(x+6,97,label)
    ax=arc_axis(s,12,80,164,13,1,n,180)
    pair_arcs(ax,bpp,g,bpp)
    ax.axis('on'); ax.spines[['left','bottom','top','right']].set_visible(False)
    ax.set_yticks([]); ax.set_xticks((1,50,100,150,200,250,300,350,388))
    ax.set_ylim(-10,180)
    refmarks={}
    for name,letter,level,locs in (('outer','O',-4,((43,108),(171,110))),
                                  ('inner','I',-8,((90,151),(115,151)))):
        pp=stats['stems'][name]['pairs']
        arms=[(min(p[k] for p in pp),max(p[k] for p in pp)) for k in (0,1)]
        refmarks[name]=arms
        for (lo,hi),xytext in zip(arms,locs):
            ax.plot([lo,hi],[level,level],color=INK,lw=1.4,ls='-' if name=='outer' else ':')
            ax.annotate(letter,xy=((lo+hi)/2,level),xytext=xytext,ha='center',va='center',
                bbox=dict(fc='white',ec='none',pad=.2),arrowprops=dict(arrowstyle='->',lw=.5,color=INK))
    ax.add_patch(Rectangle((122.5,-10),4,60,fill=False,ec=INK,lw=.8,zorder=6))
    opacity_key(s,12,74.8)
    s.text(102,74.8,'O: 6 pairs; I: 3 pairs; bulge at 133')
    s.panel(1,70,'B','Structural contexts')
    for x,c in zip((12,34,58,79,103,128,150),CONTEXTS):
        s.rect(x,64,2,2,fc=COLORS[c],ec='none'); s.text(x+2.6,66,c)
    axfull=s.axis(12,53,164,8); contexts(axfull,1,n,v['context'])
    axfull.set_xticks((1,50,100,150,200,250,300,350,388))
    axfull.add_patch(Rectangle((39.5,0),121,1.01,fill=False,ec=INK,lw=.7,ls='--',clip_on=False))
    peak=stats['pseudoloop']['region_peak']
    s.text(12,46.8,'40–160 enlarged')
    s.text(120,46.8,f"Base {peak['position']}: PseudoLoop {peak['probability']:.3f}",ha='right')
    axzoom=s.axis(12,31,108,11); contexts(axzoom,40,160,v['context'])
    axzoom.set_xticks((40,60,80,100,120,140,160)); axzoom.set_ylabel('Probability',labelpad=3)
    axzoom.set_yticks((0,.5,1),('0','0.5','1'))
    axzoom.add_patch(Rectangle((122.5,0),4,1.01,fill=False,ec=INK,lw=.8,zorder=4))
    s.panel(1,25,'C','4-nt window accessibility by loop type')
    acc={(i,w):p for i,w,p in v['access']}; selected=acc[123,4]
    s.text(12,21,f"123–126: total {selected['total']:.3f}; PseudoLoop {selected['P']:.3f}")
    axc=s.axis(12,7.2,108,11); xs=np.arange(40,161); bottom=np.zeros(122)
    for key,c in terc.LOOPS.items():
        vv=np.array([acc[int(i),4][key] for i in xs]); vv=np.append(vv,vv[-1])
        axc.fill_between(np.arange(39.5,161.5),bottom,bottom+vv,step='post',color=COLORS[c],lw=0)
        bottom+=vv
    axc.plot(xs,[acc[int(i),4]['total'] for i in xs],color=INK,lw=.8)
    axc.set(xlim=(39.5,160.5),ylim=(0,1.01),yticks=(0,.5,1),yticklabels=('0','0.5','1'),
            xticks=(40,60,80,100,120,140,160),xlabel='Start of 4-nt window (nt)',ylabel='Probability')
    axc.axvline(123,color=INK,lw=.6,ls=':')
    axc.xaxis.labelpad=2
    s.text(130,46.8,'Samples / 1,000 draws',weight='bold')
    for name,label,y in (('outer','All 6 O',41),('inner','All 3 I',35.5),('both','All 9 O + I',30)):
        fraction=stats['stems'][name]['complete_sample_fraction']
        s.text(130,y,label)
        s.rect(149,y-2.4,16,1.8,fc='#eeeeeb',ec='none')
        s.rect(149,y-2.4,16*fraction,1.8,fc=ORANGE,ec=INK,lw=.4)
        s.text(176,y,f'{fraction:.1%}',ha='right')
    s.text(130,22,'γ-centroid, γ = 1',weight='bold')
    cp=terc.pairs_of(v['centroid']); local={p for p in cp if 40<=p[0]<p[1]<=160}
    pair_arcs(arc_axis(s,130,10,46,9,40,160,55),local,g)
    s.text(130,7.8,'Region 40–160; arcs at full opacity')
    s.fig.canvas.draw()
    alignment=max(abs(axzoom.transData.transform((i,0))[0]-axc.transData.transform((i,0))[0]) for i in xs)
    assert alignment < 1e-10
    return s.fig,dict(probability_display_cutoff_inclusive=CUTOFF,context_hatches='none',
        highlighted_bases_A_B=[123,126],reference_bulged_base=133,reference_marks=refmarks,
        sample_fractions={k:stats['stems'][k]['complete_sample_fraction'] for k in ('outer','inner','both')},
        centroid_gamma=1,centroid_region=[40,160],
        centroid_pairs_outside_region=len(cp-local),Bzoom_C_alignment_px=alignment,
        centroid_reference_pairs={k:stats['stems'][k]['centroid_pairs_present'] for k in ('outer','inner')},
        **v['validation'])


def save(fig,out,name,details):
    out.mkdir(parents=True,exist_ok=True)
    fig.canvas.draw(); renderer=fig.canvas.get_renderer()
    texts=[]
    for t in fig.findobj(Text):
        if not t.get_visible() or not t.get_text().strip(): continue
        if t.axes is not None and not t.axes.axison and t not in t.axes.texts: continue
        box=t.get_window_extent(renderer)
        if not box.width or not box.height: continue
        texts.append((t,box))
    outside=[t.get_text() for t,b in texts if b.x0 < -1 or b.y0 < -1 or b.x1>fig.bbox.width+1 or b.y1>fig.bbox.height+1]
    overlaps=[(t.get_text(),u.get_text()) for i,(t,b) in enumerate(texts) for u,bb in texts[i+1:] if b.overlaps(bb)]
    assert not outside, outside
    assert not overlaps, overlaps
    assert min(t.get_fontsize() for t,b in texts)>=7
    for suffix in ('pdf','png'): fig.savefig(out/(name+'.'+suffix),dpi=300)
    check=review.pdf_text_check(out/(name+'.pdf'))
    page=PdfReader(out/(name+'.pdf')).pages[0]
    dims=[float(page.mediabox.width)*25.4/72,float(page.mediabox.height)*25.4/72]
    assert abs(dims[0]-178)<.01
    dump(out/'validation.json',dict(pdf_page_mm=dims,pdf_text_check=check,text_outside_page=outside,
        text_overlaps=overlaps,source_data_sha256=sha(SOURCE/'data.json'),
        script_sha256=sha(Path(__file__)),visual_inspection='Open the generated PNG separately.',**details))
    plt.close(fig)


def notes(out,name,caption,feedback,extra=''):
    (out/'caption.md').write_text(caption)
    figure_path = 'experiments/'+str(out.relative_to(ROOT)) if out.is_relative_to(ROOT) else str(out)
    (out/'README.md').write_text(f'''# Figure feedback — 2026-10-20

New manuscript figure: `{name}.pdf` (178 mm wide); preview: `{name}.png` (300 dpi).
Proposed paired Japanese/English caption: [caption.md](caption.md).
The manuscript editor should replace the old image path with `{figure_path}/{name}.pdf` and use that caption.

Feedback addressed: {feedback}.
Source data: `results/local/figure-feedback-20261020/data.json`.
Raw CLI output, frozen inputs and saved Terc output are packed in the adjacent `raw_outputs.tar`.
Commands, parameter/binary/source hashes and checks are in `calls.json` and `provenance.json` there.
Figure 1 was computed with LinearKnot 9c7d6e1 (exact, Turner04, 37°C); Terc probabilities and 1,000 samples are saved 89baea52 outputs.
The newly computed γ-centroid uses `linearknot centroid --gamma 1`; Terc uses beam/hairpin/depth 200 and exact for comparison.
These small local calculations are not performance measurements.
All pair-probability arc plots use p ≥ 0.05; complete probability data remain in the source archive.
When inserting the figure, use the caption in this directory and the pair-probability threshold 0.05.
Internal is a plain pink fill; no hatching is used.
`validation.json` checks actual PDF text sizes, physical dimensions, text overlap and source-data identity.
{extra}

Regenerate from experiments/: `python3 -B scripts/make_figures_feedback_20261020.py --examples-only` (use `--output-root NEW_DIRECTORY` to preserve this rendering).
To regenerate the raw small calculations in a new location: `python3 -B scripts/prepare_feedback_figure_data.py --out NEW_DIRECTORY`.
No LaTeX, bibliography or external implementation tree was edited.
''')


def example_captions(data):
    f=data['figure1']; n=sum(c for c,db in f['samples']); a,b=(c for c,db in f['samples'][:2])
    caption=f'''**Figure 1. Conditional ensemble analysis with LinearKnot.**

Aでは、与えたscaffold Gに、対合塩基を共有しないextension G′を加えてRを作る。
GとG′はそれぞれ交差対を持たない。
灰色と橙色の下線は、互いに交差する入れ子対の二つの群（band）の範囲を示す。
縦の点線は二つの範囲を通る。

**(A)** Adding extension G′ to the supplied scaffold G gives structure R.
G and G′ share no paired bases, and each has no crossing pairs.
The grey and orange underlines show the spans of two groups of nested pairs, called bands, that cross each other.
The dotted vertical line passes through both spans.

BではCPartyの分解を二つずつの小さい組合せに分け、右端ごとに上位b候補を保持する（模式図はb=3）。
丸は保持、×は破棄を示す。
同じ分解でBoltzmann重みの和を取ると分配関数、エネルギーの最小を取ると最小自由エネルギー（MFE）構造を得る。

**(B)** LinearKnot divides CParty's decomposition into smaller combinations of two parts and keeps the top b candidates at each right endpoint (b=3 in this sketch).
Circles mark retained candidates; crosses mark discarded candidates.
The same decomposition gives the partition function by summing Boltzmann weights and a minimum-free-energy (MFE) structure by taking the minimum energy.

CはbpRNA_RFAM_31397（61 nt）のexact計算を示す。
Gは参照の入れ子のcanonical対（AU/UA、CG/GC、GU/UG）の8対である。
左上の対確率図のみ、確率0.05以上を表示し、弧の不透明度を対確率とする（数値キー参照）。
他の構造図の弧は不透明である。
灰色実線はG、青点線はGと交差しない追加対、橙破線はGと交差する追加対を示す。

**(C)** Exact calculations use bpRNA_RFAM_31397 (61 nt), with eight canonical nested reference pairs (AU/UA, CG/GC or GU/UG) as G.
Only the pair-probability plot at upper left uses arc opacity to show probability, as illustrated by its numeric key; it displays pairs with probability at least 0.05.
Arcs in the other structure diagrams have full opacity.
Solid grey arcs show G, dotted blue arcs show new pairs that do not cross G, and dashed orange arcs show new pairs crossing G.

各位置の七つのcontext確率の和は1であり、PseudoLoopはpseudoknotに伴う非対合loopを示す（定義はMethods）。
accessibilityは窓内の全塩基が同時に非対合である確率であり、実線・破線・点線は1・4・8 ntの窓を示す。
座標は1始まりで、accessibilityでは窓の開始位置を表す。
右上は{n:,}回のBoltzmann標本化で最多と2番目の構造（{a:,}回と{b:,}回）を示す。
残る{n-a-b:,}回の標本は示さない。
右下はγ=1のγ-centroid構造である。

The seven context probabilities sum to one at each position.
PseudoLoop denotes unpaired loops associated with a pseudoknot, as defined in Methods.
Accessibility is the probability that all bases in a window are unpaired together; solid, dashed and dotted lines show 1-, 4- and 8-nt windows.
Coordinates are 1-based positions or window starts.
The upper-right diagrams show the two most frequent structures in {n:,} Boltzmann draws ({a:,} and {b:,} draws); the other {n-a-b:,} draws are not shown.
The lower-right diagram shows the γ-centroid structure at γ=1.
'''
    out=OUT/'figure1_89baea52_20261020_fb'
    notes(out,'figure1',caption,'F22–F27, F84',
          'A uses the reference extension as a schematic; C shows model predictions. The small extra scaffold hairpin is not part of the crossing band underlined in A.')
    t=data['terc']; v=t['settings']['b200']; selected=next(p for i,w,p in v['access'] if i==123 and w==4)
    stems=v['stats']['stems']; dif=t['differences']
    caption=f'''**Figure 3. Loop profiles, accessibility and samples for rat Terc.**

ラットTerc（NR_001567.1、388 nt）をbeam幅200、Turner04、37°Cで計算した。
GはbpRNA_RFAM_4317参照の入れ子のcanonical対61対である。
Aは確率0.05以上の対を示す。
灰色実線はG、青点線はGと交差しない追加対、橙破線はGと交差する追加対である。
弧の不透明度は対確率に等しく、重なった弧は濃く見えることがある。
OとIは参照pseudoknotの外側6対と内側3対であり、塩基133のbulgeで区切られる。
両腕のOとIの指示線がそれぞれの群を示す。

Rat Terc (NR_001567.1; 388 nt) was analyzed at beam width 200 with Turner04 parameters at 37°C.
G contains 61 canonical nested pairs from the bpRNA_RFAM_4317 reference.
**(A)** Arcs show pairs with probability at least 0.05.
Solid grey arcs show G, dotted blue arcs show new non-crossing pairs, and dashed orange arcs show new pairs crossing G.
Arc opacity equals pair probability, as illustrated by the numeric key; overlapping arcs may appear darker.
O and I mark the six outer and three inner reference pseudoknot pairs on both arms.
A bulge at base 133 separates these two groups.

Bの積み上げた色は、各位置の七つのcontext確率を示す。
破線の枠は40–160の拡大区間であり、拡大図とCの座標は共通である。
AとBの拡大図の実線の枠は塩基123–126を示す。
Internalは桃色で示し、斜線は用いない。

**(B)** The colored segments in each stacked bar show the seven context probabilities at one base.
The dashed rectangle marks the enlarged region 40–160; its enlarged view and C share coordinates.
Solid rectangles in A and the enlarged B view mark bases 123–126.
Internal uses a plain pink fill.

Cの各色は、窓の全4塩基が同時に非対合でそのloop型に属する確率を示す。
黒線はその和であるtotal accessibilityを示す。
縦の点線は窓の開始位置123を示す。
123–126のtotalは{selected['total']:.3f}、PseudoLoop成分は{selected['P']:.3f}である。

**(C)** Each colored component is the probability that all four window bases are unpaired and belong to that loop type.
The black line shows their sum, the total accessibility.
The dotted vertical line marks window start 123.
Window 123–126 has total accessibility {selected['total']:.3f}, including a PseudoLoop contribution of {selected['P']:.3f}.

右の小さい棒は1,000標本のうちOの6対すべて、Iの3対すべて、両群9対すべてを含む割合を示す。
対応する回数は{stems['outer']['complete_sample_count']}、{stems['inner']['complete_sample_count']}、{stems['both']['complete_sample_count']}である。
他の対の有無は問わず、割合には標本誤差がある。
その下の図はγ=1のγ-centroidの40–160を不透明な弧で示す。
全座標は1始まりで、Cは窓の開始位置を表す。

The small bars at right show the fractions of 1,000 samples containing all six O pairs, all three I pairs, or all nine pairs together.
The corresponding counts are {stems['outer']['complete_sample_count']}, {stems['inner']['complete_sample_count']} and {stems['both']['complete_sample_count']}.
Other pairs are allowed, and these fractions have sampling uncertainty.
Below the bars, the γ-centroid diagram (γ=1) shows region 40–160 with fully opaque arcs.
Coordinates are 1-based; C gives window starts.
'''
    out=OUT/'figure4_terc_89baea52_20261020_fb'
    notes(out,'figure4',caption,'F27, F40, F44–F52',
          f"The old D pair table and MFE/MEA counts are removed. The figure is now 178 × 104 mm (previously 178 × 100 mm). Exact comparisons near 40–160 are in data.json: context max difference {dif['context_region_max_abs_difference'][0]:.6g}, 4-nt total accessibility {dif['accessibility']['4']['total_region_max'][0]:.6g}. The stored centroid output supplies the inset. At both exact and b200, it contains all six outer reference pairs and none of the three inner reference pairs; these replace the old MEA counts.")


def resources():
    out=OUT/'main_89baea52_20261020_fb'
    subprocess.run([sys.executable,'-B',str(ROOT/'scripts/make_figures.py'),
        '--review-layout','--figure3-only','--published-only','--e2',
        str(ROOT/'results/local/e2-89baea52-20261018'),'--e2-revision','89baea52','--out',str(out)],check=True)
    caption='''**Figure 4. Time and memory for long RNAs and synthetic controls.**

AとBはRNAcentralの7入力（4,793–29,106 nt）で、三つの解析が要求する出力を含む実行全体の時間とメモリを示す。
LinearKnotのfoldは分配関数・塩基対確率（BPP）・最小自由エネルギー（MFE）構造、profileはloop profile、accessは1・4・8・16 nt窓のaccessibilityを求める。
計算にはbeam幅b=100、hairpin長の上限h=100、scaffoldのmultiloopを外側へ探す段数の上限d=100を使う。

**(A, B)** Whole-program time and memory are shown for seven RNAcentral inputs of 4,793–29,106 nt.
The legend gives the roles of the three LinearKnot subcommands: fold computes the partition function, base-pairing probabilities (BPPs) and a minimum-free-energy (MFE) structure; profile computes loop profiles; access computes accessibility for 1-, 4-, 8- and 16-nt windows.
Calculations use beam width b=100, hairpin-length limit h=100, and a limit d=100 on how many enclosing scaffold multiloop levels are searched outward.

記号はwarm-up一回を除いた3測定の中央値、縦の誤差線はその最小から最大の範囲である。
この範囲の多くは記号より小さく、見えにくい。
各長さは異なるRNAとGを表す。
メモリはプロセス全体のpeak resident set size（RSS）である。

Symbols show the median of three measurements after one excluded warm-up.
Vertical error bars extend from the smallest to the largest of those three measurements.
Many ranges are smaller than the symbols and are difficult to see.
Each length represents a different RNA and scaffold G.
Memory is the peak resident set size (RSS) of the whole process.

Cは7種類のscaffoldの形を持つ合成入力で、4,000 ntのfoldの時間とRSSを2,000 ntの値で割った比を示す。
各形・長さに対して、乱数seedを変えて3本の合成配列を生成した。
各配列の3測定の中央値を取り、次に3配列の中央値を取ってから、長さ間の比を求めた。
乱数は配列生成に用い、foldの計算は決定的である。
全長さの測定とscaffoldの図は補足図S2に示す。

**(C)** Ratios divide fold time and RSS at 4,000 nt by those at 2,000 nt for synthetic inputs with seven scaffold shapes.
For each shape and length, three random seeds generate three synthetic RNA sequences.
We take the median of three measurements for each sequence, then the median across the three sequences, before computing the ratio between lengths.
The random seeds control sequence generation; fold itself is deterministic.
Figure S2 shows all measured lengths and sketches of the scaffold shapes.

**Figure S2. Time and memory for seven synthetic scaffold shapes.**

A/Bはfold、C/Dはprofile、E/Fはaccessの時間とRSSである。
各解析の出力はFigure 4のcaptionに示す。
RNA長の対数軸は底2、時間とRSSは底10である。
各点の中央値の取り方はFigure 4Cと同じで、縦の誤差線は3配列×3測定の9測定の最小から最大を示す。
この範囲の一部は記号より小さい。
線は同じscaffoldの形の測定点を結ぶ。
全156入力のうち、7種類の形×5長さ×3配列の105入力を示す。
残る51入力は、scaffoldのmultiloopの入れ子の深さを変えた追加対照である。
空でないGはfloor(0.15n)対（nは配列長）を持ち、配列のGC含量は50%である。
下の弧図は保存済みの250 ntのGを同じ縮尺で示す。

**(A, B)** show fold, **(C, D)** show profile, and **(E, F)** show access; requested outputs are listed in the Figure 4 caption.
RNA length uses a base-2 logarithmic axis; time and RSS use base 10.
Symbols use the medians defined in Figure 4C; vertical error bars span the smallest to largest of the nine measurements (three generated sequences, three measurements each).
Some ranges are smaller than the symbols.
Lines connect measurements with the same scaffold shape.
The figure shows 105 of the 156 inputs: seven shapes, five lengths and three generated sequences.
The other 51 inputs are additional controls with different depths of nested scaffold multiloops.
Each nonempty G contains floor(0.15n) pairs, where n is sequence length; sequences have 50% GC content.
The arc sketches show the saved 250-nt scaffolds on a common coordinate scale.
'''
    (out/'caption.md').write_text(caption)
    (out/'README.md').write_text('''# Time and memory — author feedback, 2026-10-20

Figure 4: `figure3_main.pdf` / `.png`, 178 × 55 mm.
Figure S2: `figure3_synthetic_all_cli_si.pdf` / `.png`, 165 × 178 mm.
Use the new paths under `experiments/results/figures/main_89baea52_20261020_fb/`.
Proposed paired captions are in [caption.md](caption.md).
Source: `results/local/e2-89baea52-20261018/`; source revision is 89baea52.
The updated names describe the matching LinearKnot subcommands; these are the saved 89baea52 timings, not new timings for the renamed CLI.
The 9c7d6e1 output comparison is documented in `docs/linearknot_9c7d6e1_20261020/README.md`.

F57–F64, F74, F88–F89: replace CLI names by output roles, explain b/h/d and the range of three measurements, define seeds as sequence-generation seeds, and name the synthetic scaffold shapes.
The two multiloop designs are labelled Nested multiloops and Balanced multiloops.
No measured value was changed; `plot_data.tsv` matches the previous r4 figure exactly.
The old crossed-out caption sentences about bounds, causal isolation, time breakdown and the dotted ratio line are omitted.
`validation.json` checks source quantiles, plot geometry, overlapping text and actual PDF text size.
Run `python3 -B scripts/make_figures_feedback_20261020.py --resources-only` from experiments/.
''')
    old=OUT/'main_89baea52_20261018_r4/plot_data.tsv'
    assert (out/'plot_data.tsv').read_bytes()==old.read_bytes()
    manifest=json.loads((out/'manifest.json').read_text())
    manifest['feedback_script_sha256']=sha(Path(__file__))
    for name in ('README.md','caption.md'): manifest['derived_files'][name]=sha(out/name)
    dump(out/'manifest.json',manifest)


def accuracy():
    out=OUT/'e1_89baea52_20261020_fb'
    subprocess.run([sys.executable,'-B',str(ROOT/'scripts/render_e1_89baea52.py'),'--out',str(out)],check=True)
    p=out/'caption.md'; text=p.read_text()
    text=text.replace('The horizontal positions group RNAs by length and space the tested settings equally within each group.',
                      'The horizontal positions group RNAs by length and space beam widths equally within each group.')
    text=text.replace('横軸は長さ群であり、群の中では設定を等間隔に並べる。','横軸は長さ群であり、群の中ではbeam幅を等間隔に並べる。')
    text=text.replace('The local search limits h and d are set to max(b,30)',
                      'The limit h bounds the length of a new hairpin.\nThe limit d bounds how many enclosing scaffold multiloop levels are searched outward.\nBoth equal max(b,30)')
    text=text.replace('hとdは局所探索の上限であり、h=d=max(b,30)とした。','hは新しいhairpinの長さの上限、dは外側へ探すscaffold multiloopの段数の上限であり、両方をmax(b,30)とした。')
    text=text.replace('hとdは局所探索の上限で、h=d=max(b,30)とした。','hは新しいhairpinの長さの上限、dは外側へ探すscaffold multiloopの段数の上限であり、両方をmax(b,30)とした。')
    text=text.replace('## Figure S2 — figure2_d_si.pdf',
                      '## Additional accuracy companion — figure2_d_si.pdf (outside the manuscript)')
    p.write_text(text)
    (out/'README.md').write_text('''# Accuracy figures — author feedback, 2026-10-20

The annotated Figure 2 (main p5) and SI accuracy plot (p20) contain no direct figure comments.
The new paths are under `experiments/results/figures/e1_89baea52_20261020_fb/`.
Use `figure2_main.pdf` and `figure2_accessibility_si.pdf` (S1), with the paired caption proposals in caption.md.
`figure2_d_si.pdf` is the existing extra threshold companion; the manuscript already omits this duplicate plot. The synthetic scaffold plot is now S2 after NEAT1 figure removal.
The caption names beam widths instead of unspecified settings, defines h as the new-hairpin length limit, and defines d as the number of enclosing scaffold multiloop levels searched outward.
No old CLI names occur in the image labels or caption. Numerical data and plot layout are unchanged.
`plot_data.tsv` and `evaluation_counts.tsv` are byte-identical to the previous r4 figure.
Finite-beam and exact calculations both use 89baea52. Long RNAs use the completed `e1-bprna-long-89baea52-cpu16-20261018` run and the finite b=400 reference.
The RNA cohorts, error definitions and linear-interpolation quantiles are unchanged.
Figure 2 is 178 × 92 mm; S1 is 165 × 137 mm; the extra threshold plot is 165 × 55 mm. PNGs are 300 dpi.
`validation.json` records 1,750 checked quantile scalars, 350 checked summaries, plot geometry, text overlap and embedded PDF font checks.
`manifest.json` records source/output/script hashes; `evaluation_counts.tsv` gives the eligible RNA and window counts.
Regenerate: `python3 -B scripts/make_figures_feedback_20261020.py --accuracy-only --output-root NEW_DIRECTORY`.
The E1 renderer refuses an existing output directory.
''')
    for name in ('plot_data.tsv','evaluation_counts.tsv'):
        assert (out/name).read_bytes()==(OUT/'e1_89baea52_20261019_r4'/name).read_bytes()
    manifest=json.loads((out/'manifest.json').read_text())
    manifest['feedback_script_sha256']=sha(Path(__file__))
    for name in ('README.md','caption.md'): manifest['derived_files'][name]=sha(out/name)
    dump(out/'manifest.json',manifest)


def main():
    global OUT
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-root',type=Path,default=OUT,
                    help='Parent of the four new _fb directories')
    group=ap.add_mutually_exclusive_group()
    for name in ('examples','resources','accuracy'): group.add_argument('--'+name+'-only',action='store_true')
    args=ap.parse_args(); OUT=args.output_root.resolve(); setup()
    if not (args.resources_only or args.accuracy_only):
        provenance=json.loads((SOURCE/'provenance.json').read_text())
        assert sha(SOURCE/'data.json')==provenance['data_sha256']
        lock=json.loads((ROOT/'manifests/software-linearknot-9c7d6e1.json').read_text())
        assert provenance['sources']['build/linearknot-9c7d6e1/linearknot']==lock['binaries']['linearknot']
        assert provenance['sources']['deps/LinearKnot-9c7d6e1/params/rna_Turner04.par']==lock['energy_parameter_sha256']
        assert all(sha(ROOT/name)==value for name,value in provenance['sources'].items())
        data=json.loads((SOURCE/'data.json').read_text())
        for key,draw,directory,name in (('figure1',figure1,'figure1_89baea52_20261020_fb','figure1'),
                                       ('terc',terc_figure,'figure4_terc_89baea52_20261020_fb','figure4')):
            fig,details=draw(data[key]); save(fig,OUT/directory,name,details)
            print('Rendered '+name,flush=True)
        example_captions(data)
    if not (args.examples_only or args.accuracy_only): resources()
    if not (args.examples_only or args.resources_only): accuracy()


if __name__=='__main__': main()
