#!/usr/bin/env python3
"""Figure 4 and an SI ratio panel, from verified summaries without RNA execution."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/figure4-resources-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.text import Text
from matplotlib.ticker import NullLocator, ScalarFormatter
from pypdf import PdfReader

ROOT=Path(__file__).resolve().parents[1]
MM=1/25.4
HEIGHT=54
METHODS={'cparty':('#b55a10','o','CParty: partition function'),
 'exact':('#666666','s','LinearKnot: partition function, exact'),
 'beam100':('#0072b2','D','LinearKnot: partition function, b=100'),
 'spark':('#009e73','^','Spark: minimum free energy'),
 'fold100':('#a3459a','v','LinearKnot: partition function\n+ minimum free energy, b=100')}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,o):p.write_text(json.dumps(o,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def old_y(y):return (y*58/65+HEIGHT-58)/HEIGHT
def axmm(fig,x,y,w,h):return fig.add_axes([x/178,old_y(y),w/178,h*58/65/HEIGHT])
def textmm(fig,x,y,s,**kw):return fig.text(x/178,old_y(y),s,**kw)
def physical_text(fig,x,y,s,**kw):return fig.text(x/178,y/HEIGHT,s,**kw)
def key(items):return [Line2D([],[],color=c,marker=m,ms=3.5,lw=.8,label=label) for c,m,label in items]
def clean(ax):
 ax.grid(which='major',color='.90',lw=.4)
 ax.xaxis.set_minor_locator(NullLocator());ax.yaxis.set_minor_locator(NullLocator())
 ax.set_axisbelow(True)

def save(fig,out,name,width,height):
 fig.canvas.draw();renderer=fig.canvas.get_renderer()
 visible=[t for t in fig.findobj(Text) if t.get_visible() and t.get_text() and t.get_window_extent(renderer).width>0]
 assert min(t.get_fontsize() for t in visible)>=6
 bounds=fig.bbox
 outside=[t.get_text() for t in visible if not bounds.contains(*t.get_window_extent(renderer).p0) or not bounds.contains(*t.get_window_extent(renderer).p1)]
 assert not outside,outside
 completion=[t for t in visible if t.get_text().startswith('LinearKnot completed')]
 legend_texts=[t for legend in fig.legends for t in legend.get_texts()]
 collisions=[(t.get_text(),u.get_text()) for t in completion for u in legend_texts
             if t.get_window_extent(renderer).overlaps(u.get_window_extent(renderer))]
 assert not collisions,collisions
 for suffix in ('pdf','png'):fig.savefig(out/(name+'.'+suffix),dpi=300)
 pdf=PdfReader(out/(name+'.pdf'));box=pdf.pages[0].mediabox
 assert abs(float(box.width)*25.4/72-width)<.01
 assert abs(float(box.height)*25.4/72-height)<.01
 sizes=[]
 pdf.pages[0].extract_text(visitor_text=lambda text,cm,tm,font,size:sizes.append(size) if text.strip() else None)
 assert min(sizes)>=6
 return dict(width_mm=width,height_mm=height,png_dpi=300,min_font_pt=min(sizes),text_outside_canvas=outside)

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--summary',type=Path,default=ROOT/'results/analysis/comparison-resources-20261026-final/summary.json')
 parser.add_argument('--out',type=Path,default=ROOT/'results/figures/figure4_resources_20261026_legend')
 args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=False)
 validation=ROOT/'docs/comparison_resources_20261026_completed_validation.json'
 assert json.loads(validation.read_text())['status']=='passed'
 data=json.loads(args.summary.read_text());old=ROOT/'results/figures/main_4f97fd3_20261021/plot_data.tsv'
 rows=list(csv.DictReader(old.open(),delimiter='\t'))
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.labelsize':7,'axes.titlesize':7,
  'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.5,
  'axes.spines.top':False,'axes.spines.right':False,'xtick.major.size':2,'ytick.major.size':2,
  'xtick.major.pad':1.5,'ytick.major.pad':1.5,'pdf.fonttype':42})
 fig=plt.figure(figsize=(178*MM,HEIGHT*MM));audit=[]
 oldstyles=[('LinearCParty','#0072b2','o','partition function + pairing probabilities\n+ minimum free energy'),('LinearPKCapR','#e69f00','s','loop profiles'),('LinearPKRaccess','#009e73','^','loop profiles + window accessibility')]
 fig.legend(handles=key([(c,m,l) for _,c,m,l in oldstyles]),loc='upper left',bbox_to_anchor=(0,old_y(27)),frameon=False,ncol=1,handlelength=1.3,fontsize=7.5,labelspacing=.2)
 for col,metric in enumerate(('wall_seconds','rss_gib')):
  ax=axmm(fig,12+44*col,35,31,22)
  for program,c,m,label in oldstyles:
   rs=sorted([r for r in rows if r['panel']=='3AB' and r['metric']==metric and r['program']==program],key=lambda r:int(r['length']))
   x=[int(r['length']) for r in rs];y=[float(r['median']) for r in rs]
   ax.errorbar(x,y,yerr=[[float(r['median'])-float(r['min']) for r in rs],[float(r['max'])-float(r['median']) for r in rs]],color=c,marker=m,ms=3,lw=.8,capsize=1)
   audit.extend(rs)
  ax.set_xscale('log');ax.set_yscale('log');ax.set_xlim(4500,32000)
  ax.set_xticks([5000,30000],['5k','30k'])
  ax.set_yticks([50,100,500] if col==0 else [.5,1,3],['50','100','500'] if col==0 else ['0.5','1','3'])
  ax.set_ylabel('Time (s)' if col==0 else 'Peak memory (GiB)');ax.set_xlabel('RNA length (nt)')
  clean(ax);textmm(fig,1+44*col,60,'AB'[col],weight='bold',fontsize=9)
  textmm(fig,12+44*col,60,'RNAcentral RNAs',fontsize=7)
 fig.legend(handles=key(METHODS.values()),loc='upper left',bbox_to_anchor=(.49,old_y(27)),frameon=False,ncol=1,handlelength=1.0,labelspacing=.15,fontsize=7.5)
 incomplete={method:[r for r in data['resources'] if r['method']==method and not r['resource_eligible']]
             for method in ('cparty','spark')}
 for col,metric in enumerate(('seconds','mib')):
  ax=axmm(fig,100+44*col,35,31,22)
  for method,(c,m,label) in METHODS.items():
   raw=[r for r in data['resources'] if r['method']==method and r['resource_eligible']]
   ax.scatter([r['length'] for r in raw],[r['elapsed_median_seconds'] if col==0 else r['peak_rss_median_mib'] for r in raw],s=4,c=c,marker=m,alpha=.18,linewidths=0,zorder=2)
   rs=[r for r in data['length_groups'] if r['method']==method and r['eligible_inputs']]
   ax.plot([r['length_median'] for r in rs],[r[metric] for r in rs],color=c,marker=m,ms=3,lw=.8)
   audit.extend(dict(panel='CD'[col],**r) for r in rs)
  ax.set_xscale('log');ax.set_yscale('log');ax.set_xlim(200,35000)
  ax.set_xticks([250,4000,30000],['250','4k','30k'])
  ax.set_yticks([.1,10,1000] if col==0 else [10,100,1000],['0.1','10','1k'] if col==0 else ['10','100','1k'])
  ax.set_ylim((.04,1500) if col==0 else (4,5500));ax.set_ylabel('Time (s)' if col==0 else 'Peak memory (MiB)')
  ax.set_xlabel('RNA length (nt)');clean(ax)
  textmm(fig,89+44*col,60,'CD'[col],weight='bold',fontsize=9)
  textmm(fig,100+44*col,64,'bpRNA and\nRNAcentral RNAs',fontsize=7,va='top',linespacing=1.0)
 cp=incomplete['cparty'];sp=incomplete['spark']
 # Keep completion counts as an audit artifact, without a table in the figure.
 # Plot rectangles, legends, fonts and measured values keep their physical sizes.
 groups=list(dict.fromkeys(r['group'] for r in data['length_groups']))
 assert len(groups)==12
 counts=[]
 for index,group in enumerate(groups):
  group_rows=[r for r in data['length_groups'] if r['group']==group]
  length=int(group.split('_')[1]) if group.startswith('bprna_') else int(group_rows[0]['length_median'])
  for method,y in [('cparty',4.5),('spark',1.7)]:
   row=next(r for r in group_rows if r['method']==method)
   completed,planned=row['eligible_inputs'],row['planned_inputs']
   counts.append(dict(group=group,length=length,method=method,completed=completed,planned=planned))
 for method,total in [('cparty',59),('spark',103)]:
  assert sum(r['completed'] for r in counts if r['method']==method)==total
  assert sum(r['planned'] for r in counts if r['method']==method)==107
 for method,total in [('beam100',107),('fold100',107),('exact',80)]:
  planned=[r for r in data['resources'] if r['method']==method]
  assert len(planned)==total and all(r['resource_eligible'] for r in planned)
 with (args.out/'completion_counts.tsv').open('w') as stream:
  writer=csv.DictWriter(stream,list(counts[0]),delimiter='\t',lineterminator='\n')
  writer.writeheader();writer.writerows(counts)
 checks={'figure4':save(fig,args.out,'figure4',178,HEIGHT)};plt.close(fig)
 # Preserve the former Figure 4C as a separate SI figure; no numerical change.
 fig=plt.figure(figsize=(165*MM,42*MM));ax=fig.add_axes([.29,.25,.68,.68])
 ratios=[r for r in rows if r['panel']=='3C_ratio'];series=list(dict.fromkeys(r['series'] for r in ratios))
 labels=['Empty scaffold','One long stem','Three helices','Free bases together','12 helices','Nested multiloops','Balanced multiloops']
 for j,s in enumerate(series):
  for metric,c,m,dy in [('wall_seconds','#0072b2','o',.13),('rss_gib','#222222','D',-.13)]:
   row=next(r for r in ratios if r['series']==s and r['metric']==metric)
   ax.plot(float(row['value']),6-j+dy,marker=m,c=c,mfc=c if m=='o' else 'white',ms=3,ls='none')
 ax.set_yticks(range(7),labels[::-1]);ax.tick_params(axis='y',length=0);ax.set_xlim(1,3);ax.set_xticks([1,1.5,2,2.5,3]);ax.set_xlabel('Cost at 4,000 nt divided by cost at 2,000 nt')
 ax.axvline(2,c='.6',ls=':',lw=.6);ax.grid(axis='x',c='.9')
 ratio_handles=key([('#0072b2','o','Time'),('#222222','D','Peak memory')])
 ratio_handles[1].set_markerfacecolor('white')
 ax.legend(handles=ratio_handles,loc='lower left',bbox_to_anchor=(0,0),
           bbox_transform=fig.transFigure,borderaxespad=0,frameon=False)
 checks['synthetic_ratio_si']=save(fig,args.out,'figure4_synthetic_ratio_si',165,42);plt.close(fig)
 dump(args.out/'plot_data.json',dict(comparison=data,original_ensemble_and_ratio_rows=audit+ratios))
 sources=[args.summary,validation,old,ROOT/'docs/comparison_resources_20261025.md',Path(__file__)]
 dump(args.out/'manifest.json',dict(sources={str(p.relative_to(ROOT)):sha(p) for p in sources},outputs={p.name:sha(p) for p in args.out.iterdir() if p.is_file()},validation=checks))
 dump(args.out/'validation.json',dict(**checks,cparty_failed_inputs=len(cp),spark_failed_inputs=len(sp),completion_counts=counts,linearknot_completed={'beam100':107,'fold100':107,'exact':80},png_visually_checked=False))
 print(json.dumps(checks))
if __name__=='__main__':main()
