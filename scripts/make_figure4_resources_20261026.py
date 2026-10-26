#!/usr/bin/env python3
"""Figure 4 and an SI ratio panel, from verified summaries without RNA execution."""
import argparse
from collections import Counter
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
METHODS={'cparty':('#b55a10','o','CParty PF'), 'exact':('#666666','s','LinearKnot exact PF'),
 'beam100':('#0072b2','D','LinearKnot b=100 PF'), 'spark':('#009e73','^','Spark MFE'),
 'fold100':('#a3459a','v','LinearKnot fold (PF + MFE)')}
FAILURES={'nonfinite':('#d55e00','x','Nonfinite'), 'timeout':('#c52233','^','Time limit'),
 'memory_limit':('#7030a0','s','AS limit'), 'signal_11':('#222222','+','Signal 11')}

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,o):p.write_text(json.dumps(o,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def axmm(fig,x,y,w,h):return fig.add_axes([x/178,y/65,w/178,h/65])
def textmm(fig,x,y,s,**kw):return fig.text(x/178,y/65,s,**kw)
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
 parser.add_argument('--out',type=Path,default=ROOT/'results/figures/figure4_resources_20261026_final58')
 args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=False)
 validation=ROOT/'docs/comparison_resources_20261026_completed_validation.json'
 assert json.loads(validation.read_text())['status']=='passed'
 data=json.loads(args.summary.read_text());old=ROOT/'results/figures/main_4f97fd3_20261021/plot_data.tsv'
 rows=list(csv.DictReader(old.open(),delimiter='\t'))
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.labelsize':7,'axes.titlesize':7,
  'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'axes.linewidth':.5,
  'axes.spines.top':False,'axes.spines.right':False,'xtick.major.size':2,'ytick.major.size':2,
  'xtick.major.pad':1.5,'ytick.major.pad':1.5,'pdf.fonttype':42})
 fig=plt.figure(figsize=(178*MM,58*MM));audit=[]
 oldstyles=[('LinearCParty','#0072b2','o','fold: PF + BPP + MFE'),('LinearPKCapR','#e69f00','s','profile: contexts'),('LinearPKRaccess','#009e73','^','access: contexts + accessibility')]
 fig.legend(handles=key([(c,m,l) for _,c,m,l in oldstyles]),loc='upper left',bbox_to_anchor=(0,26/65),frameon=False,ncol=1,handlelength=1.3,fontsize=6.5,labelspacing=.3)
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
  ax.set_ylabel('Time (s)' if col==0 else 'RSS (GiB)');ax.set_xlabel('RNA length (nt)')
  clean(ax);textmm(fig,1+44*col,60,'AB'[col],weight='bold',fontsize=9)
  textmm(fig,12+44*col,60,'RNAcentral',fontsize=7)
 fig.legend(handles=key(METHODS.values()),loc='upper left',bbox_to_anchor=(.50,26/65),frameon=False,ncol=2,columnspacing=.7,handlelength=1.0,labelspacing=.3,fontsize=6.5)
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
  ax.set_ylim((.04,1500) if col==0 else (4,5500));ax.set_ylabel('Time (s)' if col==0 else 'RSS (MiB)')
  ax.set_xlabel('RNA length (nt)');clean(ax)
  textmm(fig,89+44*col,60,'CD'[col],weight='bold',fontsize=9)
  textmm(fig,100+44*col,60,'Real RNAs',fontsize=7)
  rug=axmm(fig,100+44*col,7,31,5);rug.set_xscale('log');rug.set_xlim(200,35000);rug.set_ylim(-.5,1.5)
  for method,y in [('cparty',1),('spark',0)]:
   for status,(c,m,label) in FAILURES.items():
    rs=[r for r in data['resources'] if r['method']==method and status in r['measurement_statuses']]
    # Small status-specific offsets distinguish mixed failure types, not resource values.
    shift={'nonfinite':0,'timeout':.16,'memory_limit':0,'signal_11':-.16}[status]
    rug.scatter([r['length'] for r in rs],[y+shift]*len(rs),s=12,c=c,marker=m,linewidths=.7)
  rug.set_yticks([1,0],['CP 48/107','Sp 4/107']);rug.tick_params(axis='y',length=0,pad=2,labelsize=6)
  rug.tick_params(axis='x',which='both',bottom=False,labelbottom=False)
  for spine in rug.spines.values():spine.set_visible(False)
  if col==0:textmm(fig,3,8,'Failure positions only; CP: CParty, Sp: Spark.\nMixed causes share input lengths.',fontsize=6)
 fig.legend(handles=key(FAILURES.values()),loc='upper center',bbox_to_anchor=(.5,6/65),frameon=False,ncol=4,columnspacing=1.5,handlelength=1)
 cp=[r for r in data['resources'] if r['method']=='cparty' and not r['resource_eligible']]
 sp=[r for r in data['resources'] if r['method']=='spark' and not r['resource_eligible']]
 
 checks={'figure4':save(fig,args.out,'figure4',178,58)};plt.close(fig)
 # Preserve the former Figure 4C as a separate SI figure; no numerical change.
 fig=plt.figure(figsize=(165*MM,42*MM));ax=fig.add_axes([.29,.25,.68,.68])
 ratios=[r for r in rows if r['panel']=='3C_ratio'];series=list(dict.fromkeys(r['series'] for r in ratios))
 labels=['Empty scaffold','One long stem','Three helices','Free bases together','12 helices','Nested multiloops','Balanced multiloops']
 for j,s in enumerate(series):
  for metric,c,m,dy in [('wall_seconds','#0072b2','o',.13),('rss_gib','#222222','D',-.13)]:
   row=next(r for r in ratios if r['series']==s and r['metric']==metric)
   ax.plot(float(row['value']),6-j+dy,marker=m,c=c,mfc=c if m=='o' else 'white',ms=3,ls='none')
 ax.set_yticks(range(7),labels[::-1]);ax.tick_params(axis='y',length=0);ax.set_xlim(1,3);ax.set_xticks([1,1.5,2,2.5,3]);ax.set_xlabel('4,000 / 2,000 nt cost ratio')
 ax.axvline(2,c='.6',ls=':',lw=.6);ax.grid(axis='x',c='.9')
 ax.legend(handles=key([('#0072b2','o','Time'),('#222222','D','Peak RSS')]),loc='upper right',frameon=False)
 checks['synthetic_ratio_si']=save(fig,args.out,'figure4_synthetic_ratio_si',165,42);plt.close(fig)
 dump(args.out/'plot_data.json',dict(comparison=data,original_ensemble_and_ratio_rows=audit+ratios))
 sources=[args.summary,validation,old,Path(__file__)]
 dump(args.out/'manifest.json',dict(sources={str(p.relative_to(ROOT)):sha(p) for p in sources},outputs={p.name:sha(p) for p in args.out.iterdir() if p.is_file()},validation=checks))
 dump(args.out/'validation.json',dict(**checks,cparty_failed_inputs=len(cp),spark_failed_inputs=len(sp),png_visually_checked=False))
 print(json.dumps(checks))
if __name__=='__main__':main()
