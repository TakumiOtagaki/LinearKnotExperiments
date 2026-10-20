#!/usr/bin/env python3
"""Write example-search reports and scientific preview artifacts from summaries."""
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import tarfile

os.environ.setdefault('MPLCONFIGDIR', str(Path(os.environ['TMPDIR'])/'example-search-mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
import numpy as np

from analyze_e3_neat1 import CONTEXTS, pairs_of
from analyze_small_rerun_89baea52 import groups
from task_storage import ArchiveReader

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT/'results/local/figure1-candidates-20261020'
NEAT = ROOT/'results/local/neat1-pk-search-20261020'
NAMES = ['bpRNA_RFAM_31397','bpRNA_RFAM_31352','bpRNA_PDB_598',
         'bpRNA_RFAM_29552','bpRNA_PDB_379']
COLORS = ['#808080','#f0c04a','#cfa18b','#c18bd0','#6aaf8b','#82bddb','#db7474']


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def dump(path,value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')


def read_table(path):
    with path.open() as f:
        return list(csv.DictReader(f,delimiter='\t'))


def table(path,rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t')
        w.writeheader();w.writerows(rows)


def figure_report():
    rows={r['id']:r for r in read_table(FIG/'screen.tsv')}
    d=json.loads((FIG/'candidate_data.json').read_text())
    # Verify the frozen reference-G input bytes after the exploratory calculation.
    plan=json.loads((ROOT/'results/runs/pk-recovery-200-turner04-89baea52-20261018/array_plan.json').read_text())
    archive=ArchiveReader(ROOT/'results/runs/pk-recovery-200-turner04-89baea52-20261018/inputs.tar')
    entries={e['name']:e for e in plan['inputs']}
    cohort={r['id']:r for r in json.loads((ROOT/'results/local/pk-recovery-bprna-20261016/pk_sel.json').read_text())}
    for name in rows:
        r=cohort[name];g=['.']*len(r['seq'])
        for i,j in r['nested']:g[i-1]='(';g[j-1]=')'
        raw=f">{name}\n{r['seq']}\n{''.join(g)}\n".encode()
        assert sha(raw)==entries[name]['sha256']
        assert raw==archive.read(entries[name]['input_member'])
    # Original example: use its existing profile and sampling, without new CLI.
    baseline={}
    baseline_sources=[]
    for reader,job in groups('figure1'):
        if job['program'] not in ('LinearPKCapR','LinearPKSampling'):
            continue
        raw=reader.artifact(job,'stdout.txt')
        baseline_sources.append(dict(run=str(reader.path.relative_to(ROOT)),
            job=job['id'],sha256=sha(raw),command=job['command']))
        if job['program']=='LinearPKCapR':
            context={c[0]:np.array(c[1:],dtype=float) for line in raw.decode().splitlines()
                     if (c:=line.split()) and c[0] in CONTEXTS}
            p=context['PseudoLoop']
            baseline.update(pseudoloop_max=float(max(p)),pseudoloop_mid_n=int(np.sum((p>=.2)&(p<=.8))))
        else:
            counts=sorted([int(c[0]) for line in raw.decode().splitlines()
                           if len(c:=line.split())==2 and c[0].isdigit()],reverse=True)
            baseline.update(sample_n=sum(counts),sample_top1=counts[0]/sum(counts))
    baseline.update(id='bpRNA_PDB_345',n=63,crossing_max=float(rows['bpRNA_PDB_345']['crossing_max']),
                    crossing_mid_n=0,sources=baseline_sources)
    dump(FIG/'baseline.json',baseline)
    descriptions={
      NAMES[0]:'G の2本の stem と、約0.5の3対の交差 stem が明瞭。11–13 と30–31、37–45に中間の PseudoLoop が続き、積み上げ context の色が広く混ざる。上位2標本は交差 stem の有無を見比べやすい。',
      NAMES[1]:'約0.55の3対の交差 stem と、36–42 の PseudoLoop のなだらかな帯が見える。18–22には異なる中間値が並ぶ。31397 と近い参照形状だが、標本の2位以下が細かく分散する。',
      NAMES[2]:'交差対40–48・41–47が約0.35、位置22の PseudoLoop が0.483。上位2標本は短い交差 stem の有無が中心で、28.5%対27.6%とほぼ拮抗する。中間 PseudoLoop は2位置だけで、細い棒が主役になる。',
      NAMES[3]:'17–20と28–31の4対が約0.63–0.66で揃う。位置15–16、26–27が中間 PseudoLoop になり、局所的な交差 stem と context の関係を示しやすい。上位2標本は交差 stem を共有し、別の小さな hairpin の有無で分かれる。',
      NAMES[4]:'交差 stem は7対で大半が0.8台。20と51にだけ中間 PseudoLoop がある。1位50.1%という標本比率を優先する代案だが、元の図の強いコントラストが残る。'}
    selected=[]
    for rank,name in enumerate(NAMES,1):
        r=rows[name]
        selected.append(dict(rank=rank,**r,visual_note=descriptions[name]))
    table(FIG/'candidates.tsv',selected)
    dump(FIG/'selected_inputs.json',{name:dict(record=d[name]['record'],scaffold=d[name]['scaffold']) for name in NAMES})
    lines=['# Figure 1C の例候補（2026-10-20）','',
      '**推薦：bpRNA_RFAM_31397（61 nt）。** PseudoLoop は最大0.608737、0.2–0.8の位置が15か所あり、交差 stem の3対は0.495965–0.517433。位置12は0.473670、11・13は約0.50で、著者の「0.4くらい」を広い領域で見せる目的に近い。参照由来の G は8対で、短い全長図でも読みやすい。',
      '',
      '**標本比率とのトレードオフ：** 推薦例の1位は34.31%、2位は14.26%で、例示された「1位50–80%」には届かない。調べた105本のうち、交差 BPP と PseudoLoop の両方に0.2–0.8の値がある67本を1万回ずつ標本化し、1位50–80%だったのは PDB_379 の1本だけだった。以下の5本すべてがその比率条件を満たすとはしない。50–80%を必須とする場合は PDB_379 を選ぶが、中間 context が広く見えることを優先して31397を推す。',
      '',
      '## 5候補','',
      '全確率は同じ参照 G に条件付けた exact、Turner04、37°C。座標は1始まり。標本は各10,000回、seed=20261020。',
      '',
      '| 推薦順・ID | nt / G対 | PseudoLoop 最大 | 0.2–0.8位置数 | 交差対 最大 / 中間対数 | 中間の交差対の例 | 標本1位 / 2位 / 3位 |',
      '|---|---:|---:|---:|---:|---|---|']
    for r in selected:
        lines.append(f"| {r['rank']}. {r['id']} | {r['n']} / {r['G_pairs']} | {float(r['pseudoloop_max']):.6f} | {r['pseudoloop_mid_n']} | {float(r['crossing_max']):.6f} / {r['crossing_mid_n']} | {r['crossing_near_04_i']}–{r['crossing_near_04_j']}: {float(r['crossing_near_04_p']):.6f} | {float(r['sample_top1']):.2%} / {float(r['sample_top2']):.2%} / {float(r['sample_top3']):.2%} |")
    lines.extend(['','| ID | 混合contextの位置 / Gで自由な位置 | 自由位置のcontext entropy（bit、平均） | 標本の異なる構造数 |','|---|---:|---:|---:|'])
    for r in selected:
        lines.append(f"| {r['id']} | {r['context_mixed_positions']} / {r['free_positions']} | {float(r['free_context_entropy_bits']):.3f} | {r['sample_unique']} |")
    lines.extend(['',
      '「混合context」は7 contextの最大値が0.95未満の位置。Gの固定対合位置は必ず Stem=1なので、自由位置を分母にする。entropy は各位置の7確率から −Σ p log₂ p を計算した。交差対は括弧の種類で数えず、Gの少なくとも1対(k,l)と i<k<j<l または k<i<l<j を満たす追加対を数えた。0.2–0.8は両端を含む。',
      '',
      '## 図にしたときの見え方',''])
    for name in NAMES:
        lines.extend([f'- **{name}**：{descriptions[name]}'])
    lines.extend(['',
      '[5候補のプレビュー](candidate_previews.png)／[PDF](candidate_previews.pdf)。左は G（灰）と確率≥0.02の交差対（赤、濃さと太さは確率）、中央は7 context、右は上位3標本と残りの頻度。低確率の交差対と非交差の追加対は左のプレビューには描かない。実データは全確率で保持している。',
      '',
      '## 探索範囲と検証','',
      f"- 200本の参照PK回復 cohort から50–100 ntの105本を全件確認。67本に中間の交差 BPP があり、67本すべてに中間 PseudoLoop があった。残り38本は BPP 段階で除外。Gは参照の `()` の canonical（AU/UA/CG/GC/GU/UG、j−i≥4）対。完全参照には canonical な交差対がある。元の図と同じ G 生成規則で、E1 の LinearPartition 由来 G は混ぜていない。",
      '- BPP は `results/runs/pk-recovery-200-turner04-89baea52-20261018/` の保存済み exact、cutoff=0を使用。105入力の配列・Gを frozen inputs.tar のbytesとhashに照合した。参照の `[]` に限らず、canonicalな全参照対とGの交差を確認した。',
      f"- 元の bpRNA_PDB_345（63 nt）は交差 BPP 最大{baseline['crossing_max']:.6f}、中間の交差対0本。既存Figure 1出力の PseudoLoop 最大{baseline['pseudoloop_max']:.6f}、中間位置{baseline['pseudoloop_mid_n']}か所、{baseline['sample_n']}標本の1位{baseline['sample_top1']:.2%}。値は `baseline.json`。",
      '- context と sampling だけを pinned `build/linear-cparty-89baea52/` で追加計算。binaryとTurner04 parameterのhashを既存lockに前後照合し、contextの総和、全標本のG保持、BPPとの標本誤差を確認した。',
      f"- 5候補の交差 BPP と標本頻度の最大絶対差は {max(float(r['sample_crossing_bpp_max_abs_error']) for r in selected):.6f}。1万回の比率は標本推定値で、母確率そのものではない（比率0.5の標準誤差は約0.5 percentage point）。",
      '- 一般的なbpRNA全体の最良例を証明する探索ではなく、この200本の既存cohort内の探索。67本を同じ標本数で比較した。',
      '',
      '## 出力と再現','',
      '- `candidates.tsv`：5候補の全指標と見え方。`selected_inputs.json`：配列・参照・G。',
      '- `bpp_screen.tsv` / `screen.tsv`：105本の探索記録。context未計算の欄は空欄で、ゼロではない。',
      '- `candidate_data.json`：67本の context、交差 BPP、上位20標本。`position_profiles.tsv.gz` と `sampling_structures.tsv.gz`：全67本の全位置・全標本集計。',
      '- `raw_small_calculations.tar.gz`：入力と追加CLIのstdout/stderr。`calls.json` と `provenance.json`：実行、hash、出典。',
      '- tarは共有領域へ展開していない。新しいCLI計算はTMPDIRで行い、生出力を1個のtarに集約。原稿・外部実装の変更、qsub、git commitは行っていない。',
      '',
      '実験repoから、既存ディレクトリを上書きせず再探索する場合：','',
      '```bash','python3 -B scripts/search_examples_20261020.py figure --out "$TMPDIR/figure-search-reproduce"','```','',
      'レポートとプレビューの再生成：`python3 -B scripts/report_example_search_20261020.py`。'])
    (FIG/'README.md').write_text('\n'.join(lines)+'\n')
    fig,axes=plt.subplots(5,3,figsize=(15,15),gridspec_kw={'width_ratios':[1.2,1.5,.8]})
    for row,name in enumerate(NAMES):
        v=d[name];r=rows[name];n=len(v['record']['seq'])
        a,b,c=axes[row]
        for i,j in sorted(pairs_of(v['scaffold'])):
            a.add_patch(Arc(((i+j)/2,0),j-i,(j-i)*.65,theta1=0,theta2=180,color='#888888',lw=1))
        for pair in sorted(v['crossing_pairs'],key=lambda x:x['probability']):
            i,j,p=pair['i'],pair['j'],pair['probability']
            if p>=.02:
                a.add_patch(Arc(((i+j)/2,0),j-i,(j-i)*.65,theta1=180,theta2=360,
                    color='#bf3b40',alpha=.15+.85*p,lw=.6+2.5*p))
        a.axhline(0,color='#aaaaaa',lw=.4)
        a.set_xlim(0,n+1);a.set_ylim(-n*.33,n*.33);a.set_yticks([])
        a.set_title(name+f" ({n} nt)",fontsize=11,loc='left')
        a.set_xlabel('Position; G above, crossing pairs below')
        b.stackplot(np.arange(1,n+1),*[v['contexts'][x] for x in CONTEXTS],colors=COLORS,step='mid')
        b.set_xlim(1,n);b.set_ylim(0,1);b.set_yticks([0,.5,1]);b.set_ylabel('Probability')
        b.set_title(f"PseudoLoop max {float(r['pseudoloop_max']):.3f}; {r['pseudoloop_mid_n']} intermediate positions",fontsize=10)
        b.set_xlabel('Position')
        ps=[float(r['sample_top1']),float(r['sample_top2']),float(r['sample_top3'])]
        ps.append(1-sum(ps))
        c.barh(['Rank 1','Rank 2','Rank 3','Others'],ps,color=['#406794','#6e8eb4','#a5bad2','#d1d1d1'])
        c.invert_yaxis();c.set_xlim(0,1)
        for idx,p in enumerate(ps):c.text(p+.02,idx,f'{p:.1%}',va='center',fontsize=9)
        c.set_title('10,000 exact samples',fontsize=10);c.set_xlabel('Fraction')
    fig.legend([plt.Rectangle((0,0),1,1,color=x) for x in COLORS],CONTEXTS,
               loc='lower center',ncol=7,frameon=False)
    fig.tight_layout(rect=(0,.035,1,1))
    fig.savefig(FIG/'candidate_previews.png',dpi=160)
    fig.savefig(FIG/'candidate_previews.pdf')
    plt.close(fig)
    maxerr=max(float(r['sample_crossing_bpp_max_abs_error']) for r in selected)
    assert maxerr<.03
    # Independently aggregate the packed raw context and unique-sample outputs.
    calls=json.loads((FIG/'calls.json').read_text())
    max_context_sum_error=0.
    raw_sample_totals={}
    checked=0
    with tarfile.open(FIG/'raw_small_calculations.tar.gz') as archive:
        for call in calls:
            name,program=call['id'],call['program']
            raw=archive.extractfile(name+'/'+program+'.stdout').read()
            stderr=archive.extractfile(name+'/'+program+'.stderr').read()
            assert sha(raw)==call['stdout_sha256'] and sha(stderr)==call['stderr_sha256']
            assert call['returncode']==0
            if program=='LinearPKCapR':
                context={c[0]:np.array(c[1:],dtype=float) for line in raw.decode().splitlines()
                         if (c:=line.split()) and c[0] in CONTEXTS}
                error=float(np.max(np.abs(sum(context.values())-1)))
                max_context_sum_error=max(max_context_sum_error,error)
                assert error<1e-8
                assert all(np.array_equal(context[c],d[name]['contexts'][c]) for c in CONTEXTS)
            else:
                raw_counts={c[1]:int(c[0]) for line in raw.decode().splitlines()
                    if len(c:=line.split())==2 and c[0].isdigit()}
                assert sum(raw_counts.values())==10000
                assert all(pairs_of(d[name]['scaffold'])<=pairs_of(db) for db in raw_counts)
                for count,db in d[name]['top_samples']:
                    assert raw_counts[db]==count
                raw_sample_totals[name]=sum(raw_counts.values())
            checked+=1
    dump(FIG/'validation.json',dict(status='passed',screened_inputs=105,
         frozen_input_hashes_checked=105,context_and_sampling_inputs=67,sample_n_per_input=10000,
         selected=NAMES,recommended=NAMES[0],selected_sample_cross_bpp_max_abs_error=maxerr,
         raw_calculation_outputs_checked=checked,context_sum_max_error=max_context_sum_error,
         sampled_structures_contain_G=True,raw_total_samples=sum(raw_sample_totals.values()),
         baseline=baseline,
         sample_top1_50_80_count=sum(bool(r['sample_top1']) and .5<=float(r['sample_top1'])<=.8 for r in rows.values()),
         sample_top1_criterion_met_by_recommendation=False,
         manuscript_modified=False,external_sources_modified=False,git_commit=False))


def neat_report():
    full=read_table(NEAT/'full_length_summary.tsv');windows=read_table(NEAT/'ranked_windows.tsv')
    table(NEAT/'top5_windows.tsv',windows[:5])
    lines=['# NEAT1 全長の交差対探索（2026-10-20）','',
      '**結論：確率≥0.05、≥0.1の交差対を示せる区間は、保存された全長出力にはない。** 主設定 G0.7 の全長最大は、NEAT1_1・b=400 の対6–3732の **0.0122941016**。NEAT1_2はb=100/200/400の全幅で、全長の交差対確率が **0.001未満**。空のG・G0.5・G0.9の保存出力にも0.05以上の交差対はない。',
      '',
      '**図に使うなら限定的な候補は NEAT1_1 の3457–3756（3′端）と1–300（5′端）の組。** b=400の0.012前後の遠距離対を、両端を対応させた図で示せる。ただしb=200には出力閾値を超える対応対がなく、安定した強いpseudoknotの例には推さない。SIの「交差対が見えない」問題を、表示窓の移動だけで解消できるという結果ではなかった。',
      '',
      '## 表の読み方','',
      '- 交差は追加対(i,j)がGの少なくとも1対(k,l)と i<k<j<l または k<i<l<j を満たすこと。括弧ラベルや、追加対どうしの交差は使わない。',
      '- 保存済み BPP の出力閾値は **0.001**。出力されない対の確率はゼロと断定せず **<0.001** と表示する。閾値≥0.05/≥0.1の本数は出力から漏れなく数えられる。',
      '- 窓の「交差対数」は**少なくとも片端が窓内にある出力対数**（確率≥0.001）。両端が窓内の本数は別に記載する。NEAT1_1の両端窓は9本の遠距離対を共有するので、窓の本数を足して全長本数にはしない。',
      '- 座標は1始まり、両端を含む300 nt。b比較では実行時の局所上限も h=d=b で変わる。したがって差をbeam pruningだけの効果と断定しない。',
      '',
      '## 上位5表示窓（主設定 G0.7）','',
      '| RNA・表示窓 | 交差対max b200 → b400 | 出力対数 b400（両端内） | ≥0.05 / ≥0.1 b400 | PseudoLoop max b200 → b400 | b200→400の評価 |',
      '|---|---|---:|---:|---|---|']
    notes=['交差対はb400だけ。PseudoLoopも約9.1×10^8倍で不安定。',
           '交差対はb400だけ。PseudoLoopは約1.36×10^5倍で不安定。',
           '両幅で交差対<0.001。PseudoLoopのピーク位置866は同じ、量は約1.44倍。',
           '両幅で交差対<0.001。PseudoLoopは約4.55倍、絶対値は小さい。',
           '両幅で交差対<0.001。PseudoLoopは約7.01倍、絶対値は小さい。']
    for r,note in zip(windows[:5],notes):
        maximum=float(r['b400_cross_max']) if r['b400_cross_max'] else None
        ma=f'{maximum:.6f}' if maximum is not None else '<0.001'
        lines.append(f"| {r['iso']} {r['start']}–{r['end']} | <0.001 → {ma} | {r['b400_cross_n']}（{r['b400_cross_both_inside_n']}） | {r['b400_cross_ge005']} / {r['b400_cross_ge01']} | {float(r['b200_pseudo_max']):.6g} → {float(r['b400_pseudo_max']):.6g} | {note} |")
    lines.extend(['',
      '全25,901窓（NEAT1_1:3,457、NEAT1_2:22,444、stride=1）を走査した。表示候補はb400の各位置の最大交差 BPP、次にPseudoLoopで位置を順位付けし、その位置を中心にした窓を作った（端では切り詰めず窓をRNA内へ移動）。同じRNA内で重ならない窓を各5個選び、窓の最大交差 BPP、PseudoLoopで全体を順位付けした。表はその上位5個。`ranked_windows.tsv` に10窓、`all_300nt_windows.tsv.gz` に全窓を保持した。',
      '',
      'NEAT1_1の遠距離対のspanは3722–3730 nt。単独の300 nt窓には両端を収められない。1–300内には対1–6（0.00222848）もあるが、強い交差 stem の代わりにはならない。3′端のPseudoLoopの全長最大は位置3724の0.0308042で、b200では同位置を含む端窓の最大が3.38×10^−11。',
      '',
      '[PseudoLoopの全長図](pseudoloop_full_length.png)／[PDF](pseudoloop_full_length.pdf) は比較用に対数目盛、下限10^−12を使う。小さい値が見えること自体を大きい確率と解釈しない。',
      '',
      '## 全長・全保存設定','',
      '| RNA | G | G対数 | b | 交差対 全長最大 | 出力交差対数 | ≥0.05 / ≥0.1 | PseudoLoop 全長最大（位置） |',
      '|---|---|---:|---:|---:|---:|---:|---|'])
    for r in full:
        ma=f"{float(r['crossing_max_exported']):.9g}" if r['crossing_max_exported'] else ('0（定義）' if r['G']=='emptyG' else '<0.001')
        lines.append(f"| {r['iso']} | {r['G']} | {r['G_pairs']} | {r['beam']} | {ma} | {r['crossing_exported_n']} | {r['crossing_ge_005']} / {r['crossing_ge_01']} | {float(r['pseudoloop_max']):.9g}（{r['pseudoloop_max_position']}） |")
    lines.extend(['',
      '**G閾値を変えても増えるという証拠はない。** NEAT1_1のG0.5・G0.9はb100/200でともに出力交差対0本（各対<0.001）。G0.5ではPseudoLoop最大が約5.85×10^−5とG0.7/b200の約1.05×10^−7より大きいが、表示閾値を超える交差対は増えない。G0.9も増えない。G0.5/0.9のb400、およびNEAT1_2のG0.5/0.9は指定の保存runに存在しないので比較できない。空Gには交差対象となるG対がないので交差対数は定義上0であり、Gを緩くした感度の代用にはしない。',
      '',
      '## 出典・検証・再現','',
      '- `results/runs/e3-neat1-{bpp,capr,variants-bpp,variants-capr}-89baea52-20261018/` の14条件×2種、28保存出力をtarからstreamで読んだ。NEAT1_1=3,756 nt、NEAT1_2=22,743 nt。主設定はG0.7、b=100/200/400。',
      '- 各入力の配列・Gを凍結inputs.tarと照合し、28 member のSHA-256をJSONL記録に照合した。revision、実行成功、run後のlock検証、BPP cutoff、局所上限とrankingを確認した。7 contextの各位置総和は誤差<10^−8。',
      '- 全長BPPをGとの座標交差で独立に再分類し、既存 `results/local/e3-neat1-89baea52-20261018/bpp_pairs.tsv.gz` の交差対座標・確率と完全一致した。',
      '- `full_length_summary.tsv`：全14条件。`crossing_pairs.tsv`：出力された全10交差対。`pseudoloop_positions.tsv.gz`：全長の全位置。`top5_windows.tsv`：上表。',
      '- `provenance.json`：入力と生出力のhash・実行コマンド。`validation.json`：条件数・member数・全窓数と検証。',
      '- NEAT1の新規RNA計算は0件。原稿・既存図・外部実装への変更、tarの共有領域への展開、qsub、git commitは行っていない。',
      '',
      '```bash','python3 -B scripts/search_examples_20261020.py neat --out "$TMPDIR/neat1-search-reproduce"',
      '```','',
      'レポートと図の再生成：`python3 -B scripts/report_example_search_20261020.py`。'])
    (NEAT/'README.md').write_text('\n'.join(lines)+'\n')
    profiles={}
    with gzip.open(NEAT/'pseudoloop_positions.tsv.gz','rt') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            if r['G']=='G0.7':profiles.setdefault((r['iso'],int(r['beam'])),[]).append(float(r['PseudoLoop']))
    fig,axes=plt.subplots(2,1,figsize=(13,6),sharey=True)
    for ax,iso in zip(axes,('NEAT1_1','NEAT1_2')):
        for b,c in zip((100,200,400),('#9d9d9d','#4e83ae','#c34b4b')):
            p=np.array(profiles[iso,b]);ax.plot(np.arange(1,len(p)+1),np.maximum(p,1e-12),color=c,lw=.8,label=f'b={b}')
        ax.set_yscale('log');ax.set_ylim(1e-12,1);ax.set_xlim(1,len(p));ax.set_ylabel('PseudoLoop probability')
        ax.set_title(iso+'; G0.7; values below 1e-12 clipped for display',fontsize=10)
        ax.legend(loc='upper right',ncol=3,frameon=False);ax.set_xlabel('Position')
    fig.tight_layout();fig.savefig(NEAT/'pseudoloop_full_length.png',dpi=160)
    fig.savefig(NEAT/'pseudoloop_full_length.pdf');plt.close(fig)


def receipts():
    for folder in (FIG,NEAT):
        dump(folder/'report_provenance.json',dict(
            report_script_sha256=sha(Path(__file__).read_bytes()),
            analysis_script_sha256=sha((ROOT/'scripts/search_examples_20261020.py').read_bytes())))
        (folder/'SHA256SUMS').write_text(''.join(f'{sha(p.read_bytes())}  {p.name}\n'
            for p in sorted(folder.iterdir()) if p.is_file() and p.name!='SHA256SUMS'))


if __name__=='__main__':
    figure_report();neat_report();receipts()
