#!/usr/bin/env python3
"""Render the adopted Figure 1/3 layouts with the recomputed 0.2.0 contexts."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

import make_figure1 as f1
import make_figure4_terc as terc
import make_figures_feedback_20261020 as feedback
from rerun_pseudoloop_small_20261023 import dump, sha

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'results/local/pseudoloop-0.2.0-20261023/examples'
OUT=ROOT/'results/figures/pseudoloop-0.2.0-20261023'


def main():
    data=json.loads((DATA/'data.json').read_text())
    with tempfile.TemporaryDirectory(prefix='pl020-render-',dir=os.environ['TMPDIR']) as td:
        td=Path(td)
        with tarfile.open(DATA/'terc_data.tar') as tf:
            for m in tf:
                assert m.isfile() and Path(m.name).name==m.name
                (td/m.name).write_bytes(tf.extractfile(m).read())
        t=data['terc'];pk=set(map(tuple,t['record']['pk']))
        stems=dict(outer={p for p in pk if p[0]<=68},inner={p for p in pk if p[0]>=69},both=pk)
        settings={}
        for mode in ('b200','exact'):
            v=t['settings'][mode]
            settings[mode]=setting=terc.read_setting(td,mode,388,t['scaffold'])
            stats=terc.summarize(setting,stems,388)
            # Sample and representative structure outputs are the saved adopted data.
            for name,pairs in stems.items():
                stats['stems'][name].pop('mea_pairs_present')
                stats['stems'][name]['centroid_pairs_present']=sorted(pairs & terc.pairs_of(v['centroid']))
            v.update(context={c:x.tolist() for c,x in setting['context'].items()},
                access=[[i,w,p] for (i,w),p in sorted(setting['access'].items())],
                stats=stats,validation=setting['validation'])
        t['differences']=terc.compare(settings['b200'],settings['exact'])
    path=DATA/'figure_data.json';dump(path,data)
    OUT.mkdir(parents=True,exist_ok=False)
    subprocess.run([sys.executable,'-B',ROOT/'scripts/make_figure1.py','--data',path,'--out',OUT/'figure1'],check=True)
    feedback.SOURCE=DATA;feedback.OUT=OUT;feedback.setup()
    fig,details=feedback.terc_figure(data['terc'])
    feedback.save(fig,OUT/'figure3','figure4',details)
    # Reuse the paired captions, routing only the Figure 3 caption here.
    def notes(out,name,caption,*args):
        if name=='figure4':(OUT/'figure3/caption.md').write_text(caption)
    feedback.notes=notes
    feedback.example_captions(data)
    for dirname in ('figure1','figure3'):
        out=OUT/dirname
        p=out/'validation.json';v=json.loads(p.read_text())
        v.update(source_data_sha256=sha(path),context_release='0.2.0 / f93b55e',
            context_source=str(DATA.relative_to(ROOT)),render_wrapper_sha256=sha(__file__),
            retained_quantities='BPP, total accessibility, samples and centroid: adopted saved outputs; contexts and by-loop: freshly recomputed')
        dump(p,v)
        caption=out/'caption.md'
        text=caption.read_text().replace('pseudoknotに伴う非対合loop','交差するbandの間の非対合の隙間').replace(
            'unpaired loops associated with a pseudoknot','unpaired gaps between crossing bands')
        caption.write_text(text)
        (out/'README.md').write_text(f'# {dirname}: PseudoLoop 0.2.0, 2026-10-23\n\n'
            'Loop profiles and loop-type accessibility use the archived f93b55e source.\n'
            'BPP, total accessibility, samples and centroid retain the adopted saved outputs.\n'
            'Figure 3 retains filename figure4.pdf/figure4.png.\n'
            f'Raw calculations, input/executable hashes and commands: `{DATA.relative_to(ROOT)}`.\n'
            f'Complete plotted data: `{path.relative_to(ROOT)}`.\n'
            'Reproduce: `python3 -B scripts/render_pseudoloop_20261023.py` (choose a new OUT to retain this rendering).\n'
            'No export to the manuscript figures/ directory has been made.\n')
    dump(OUT/'provenance.json',{'figure_data_sha256':sha(path),'renderer_sha256':sha(__file__),
        'raw_provenance_sha256':sha(DATA/'provenance.json'),
        'outputs':{str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*') if p.is_file()}})


if __name__=='__main__':main()
