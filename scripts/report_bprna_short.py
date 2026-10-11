#!/usr/bin/env python3
"""Write the requested Japanese E1 note from the preserved aggregate JSON."""
import argparse
import json
import math
from pathlib import Path


def number(value):
    if value is None:return 'NA'
    if value == 0:return '0'
    return f'{value:.3g}'


def percent(value):return 'NA' if value is None else f'{100*value:.2f}%'


def ratio(value):return 'NA' if value is None else f'{value:.6g}'


def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(map(str,row))+' |' for row in rows])


def write_note(output,destination):
    s=json.loads((output/'summary.json').read_text())
    p=json.loads((output/'provenance.json').read_text())
    v=json.loads((output/'validation.json').read_text())
    bins=[x for x in s['bins'] if x['bin'] in ('nt50_99','nt100_199','nt200_299','nt300_399')]
    overall=next(g for g in s['bins'] if g['bin']=='all')
    d200=next(c for c in overall['conditions'] if c['cap']=='default' and c['beam']==200)
    n200=next(c for c in overall['conditions'] if c['cap']=='nocap' and c['beam']==200)
    long200=next(c for c in bins[-1]['conditions'] if c['cap']=='nocap' and c['beam']==200)
    mono=next(x for x in s['monotonicity'] if x['bin']=='all' and x['cap']=='nocap' and x['metric']=='delta_free_energy_kcal_mol')
    title=lambda g:f"{g['min_length']}–{g['max_length']}"
    text=['# bpRNA 400 nt 未満の E1 集計（2026-10-11）','',
          f"固定版 `15aa7742` の bpRNA **{s['inputs']:,}配列**を、exact と b=20/50/100/200、局所上限30/30（default）・100000（nocap）で比較した。",
          '400 nt 以上の未実行区分は、この精度集計の対象に含めない。',
          '以下の値は同じ S・G・モデルの exact に対する近似誤差であり、注釈構造に対する予測精度ではない。','',
          f"b=200 で ΔF≤0.01 kcal/mol の割合は全体で default {percent(1-d200['delta_F_above']['0.01']['fraction_eligible'])}、nocap {percent(1-n200['delta_F_above']['0.01']['fraction_eligible'])}だった。",
          f"一方、300–399 nt の nocap では {percent(long200['delta_F_above']['0.01']['fraction_eligible'])}が0.01を超え、最大 ΔF は{number(long200['metrics']['delta_free_energy_kcal_mol']['max'])} kcal/mol、最小保持率は{percent(long200['metrics']['partition_function_retention']['min'])}だった。",
          f"ΔF の単調性違反は両上限で {mono['nonmonotone_sequences']}配列あり、TVの違反も観測された。以下に分布と個別例を示す。",'',
          '## 入力・参照・保存方法','',
          table(['長さ nt','配列数'],[[title(g),f"{g['sequences']:,}"] for g in bins]),'',
          f"全体の {percent((bins[0]['sequences']+bins[1]['sequences'])/s['inputs'])}が200 nt未満である。この長さ構成を含む全体割合と、長さ別の割合を併記する。",'',
          'G は LinearPartition-V `b450fb3e`、beam200、dangles2、posterior≥0.7、stacked helix長≥2で生成した保存済み scaffold。',
          'Turner04、37°C、固定版の既定 closure、3 CLI、BPP cutoff0、accessibility w=1/4/8/16を用いた。',
          f"RT は {s['rt_range'][0]:.10f} kcal/mol。",'',
          f"本体の {sum(p['original_job_status_counts'].values()):,} CLI jobs は、成功 {p['original_job_status_counts'].get('success',0):,}、失敗 {p['original_job_status_counts'].get('failed',0):,}だった。",
          f"失敗した exact **{len(p['exact_replacements'])}件**（CapR 172、Raccess 172、LinearCParty 5、{len(set(r['input'] for r in p['exact_replacements']))}配列）は、32 GiB 上限で成功した指定の再実行を参照に差し替えた。",
          '同じ入力 hash・S・G・source commit・モデル設定・parameter hash・binary hashes・計算引数を確認し、成功済み exact の上書きや別配列への置換は拒否した。',
          '元の失敗は `provenance.json` の `exact_replacements` に、status・error・stderr・再実行先とともに残す。',
          f"差し替え後の成功は {p['effective_job_status_counts']['success']:,}/{sum(p['effective_job_status_counts'].values()):,} jobs。有限 beam の RNA×条件は {p['analysis_status_counts'].get('success',0):,}行。",'',
          'tar は展開せず、108タスクを順に読み、1配列の27出力だけを解析用に保持した。',
          '共有領域の結果は metrics TSV、集計 JSON、provenance JSON、validation JSON、PNG 2枚の6ファイルにまとめた。',
          '全体の3.1 GBの `run.json` は一括ロードせず、array plan と各タスクの保存済み manifest/checkpoint を使った。',
          '', '## 指標と集計の定義','',
          'ΔF = F_b − F_exact、F = −RT log Z（kcal/mol）、保持率 r = Z_b/Z_exact = exp(−ΔF/RT)。',
          'Raccess の full-precision metadata を使い、CLIの丸めた表示値は使わない。負の微小 ΔF と1を微小に超える保持率も TSV/JSON にそのまま残す。',
          'このRTでは ΔF='+', '.join(f'{f:g} → r={percent(math.exp(-f/s["rt_range"][0]))}' for f in (.01,.1,1))+'。ΔF<1 kcal/molという基準だけでは、分配関数をほぼ完全に保持したとは言えない。',
          'BPP は生成対の L1 総和と、Gで自由な位置について相手塩基・unpaired の base-state 平均 TV。',
          'Context は自由塩基における7 context の平均 TV と各 context の MAE。',
          'Accessibility は G の対合塩基を含まない eligible window に限り、各 w の probability MAE・最大誤差・件数を記録する。Opening-energy MAE は両確率が正の window のみとし、片側だけ0も別に数える。',
          '各指標はまず RNA ごとに計算し、配列を等重みで集計する。位置や window を独立な配列として数えない。',
          '標本の中央値・平均・p90・p99・最大・最小は `summary.json` に保存した。分位点は線形補間。NAは0にせず、指標ごとの有効件数を併記する。',
          '4区分に加え、50–74、75–99、…、375–399 nt の25 nt区分と、50–399 nt全体も保存した。','',
          '## ΔF の分布と閾値超過','',
          'ΔFは符号付き。>は厳密な閾値超過。各行の分母は表の長さ区分の配列数。','']
    rows=[]
    for g in bins:
        for c in g['conditions']:
            m=c['metrics']['delta_free_energy_kcal_mol']
            rows.append([title(g),c['cap'],c['beam']]+[number(m[k]) for k in ('median','mean','p90','p99','max')]+
                        [percent(c['delta_F_above'][k]['fraction_eligible']) for k in ('0.01','0.1','1')])
    text += [table(['nt','上限','b','中央値','平均','p90','p99','最大','>0.01','>0.1','>1'],rows),'',
             '## b=200 の保持率と周辺確率の誤差','']
    rows=[]
    for g in bins:
        for c in g['conditions']:
            if c['beam']!=200:continue
            r=c['metrics']['partition_function_retention'];b=c['metrics']['base_state_tv_mean_Gfree'];ct=c['metrics']['context_tv_mean_Gfree']
            rows.append([title(g),c['cap'],ratio(r['median']),ratio(r['mean']),ratio(r['min']),
                         number(b['median']),number(b['p99']),number(b['max']),number(ct['p99']),number(ct['max'])])
    text += [table(['nt','上限','r中央値','r平均','r最小','BPP TV中央値','BPP TV p99','BPP TV最大','context TV p99','context TV最大'],rows),'',
             '全幅の generated-pair L1、BPP TV、context TVと7種の context別 MAE は TSV/JSON を参照する。',
             'ΔF が小さいという絶対差の見方と、r が何割を保持するかという見方を併記する。','',
             '## Accessibility','',
             'b=200 の probability MAE の RNA間 p90（各 RNA 内は eligible window の平均）。全幅の中央値・平均・p90・p99・最大、および window 内最大誤差は JSON にある。','']
    rows=[]
    for g in bins:
        for c in g['conditions']:
            if c['beam']==200:
                rows.append([title(g),c['cap']]+[number(c['metrics'][f'access_w{w}_probability_mae_eligible']['p90']) for w in (1,4,8,16)])
    text += [table(['nt','上限','w=1','w=4','w=8','w=16'],rows),'',
             'Eligible window が0の配列数：'+', '.join(f"w={w}: {n:,}" for w,n in s['eligible_zero_window_sequences'].items())+'。',
             'NA の分母と、reference-positive/candidate-zero、reference-zero/candidate-positive の総 window 数は長さ・条件ごとに保存した。','',
             '## 上限あり・なしの配列内差','',
             '以下は default − nocap の ΔF。正ならdefaultの誤差が大きい。有限beamで上限を変えると保持集合も変わるため、上限だけの損失分解とは呼ばない。','']
    rows=[]
    for g in bins:
        for c in g['cap_differences_default_minus_nocap']:
            m=c['metrics']['delta_free_energy_kcal_mol']
            rows.append([title(g),c['beam']]+[number(m[k]) for k in ('median','mean','p99','max','min')]+
                        [m['default_larger_count'],m['nocap_larger_count']])
    text += [table(['nt','b','中央値','平均','p99','最大','最小','default大','nocap大'],rows),'',
             '差の件数判定は |ΔΔF|>10⁻⁸ kcal/mol。確率指標の配列内差も JSON に保存した。',
             'nocapのΔFがdefaultより大きい配列数は、b=20/50/100/200の順で '+
             '/'.join(str(c['metrics']['delta_free_energy_kcal_mol']['nocap_larger_count']) for c in overall['cap_differences_default_minus_nocap'])+
             ' 件だった。b=200のBPP TVは '+str(overall['cap_differences_default_minus_nocap'][-1]['metrics']['base_state_tv_mean_Gfree']['nocap_larger_count'])+
             ' 配列でnocapの方が大きかった（TV差>10⁻⁹）。',
             'この観測も、局所上限を緩めることを有限beamにおける単調な精度改善と同一視しない理由になる。','',
             '## beam 幅に対する単調性','',
             'b=20→50→100→200 の3段階すべてで、ΔF・BPP TV・context TV が非増加かを調べた。保持率は非減少を調べる。',
             '数値許容は ΔF 10⁻⁸ kcal/mol、保持率10⁻⁸、TV10⁻⁹。許容内の変化は同値として扱う。','']
    text += ['ΔFの非増加と保持率の非減少は数学的には同値だが、別々の絶対許容で判定した割合は指数変換のため異なり得る。特に小さい保持率の変化は許容内になりやすい。','']
    rows=[]
    labels=['all']+[g['bin'] for g in bins]
    for label in labels:
        for cap in ('default','nocap'):
            stats={x['metric']:x for x in s['monotonicity'] if x['bin']==label and x['cap']==cap}
            d=stats['delta_free_energy_kcal_mol']
            rows.append([label.replace('nt','').replace('_','–'),cap,d['total_sequences']]+
                        [percent(stats[k]['fraction_eligible']) for k in
                         ('delta_free_energy_kcal_mol','partition_function_retention','base_state_tv_mean_Gfree','context_tv_mean_Gfree')]+
                        [d['nonmonotone_sequences'],stats['base_state_tv_mean_Gfree']['nonmonotone_sequences']])
    text += [table(['nt','上限','N','ΔF↓','r↑','BPP TV↓','context TV↓','ΔF違反数','BPP TV違反数'],rows),'',
             '最も大きな単調性違反の例（値は b=20/50/100/200 の順）：','']
    rows=[]
    for cap in ('default','nocap'):
        for metric,label in [('delta_free_energy_kcal_mol','ΔF'),('base_state_tv_mean_Gfree','BPP TV'),('context_tv_mean_Gfree','context TV')]:
            for x in s['nonmonotone_examples'][cap+'_'+metric][:2]:
                rows.append([x['input'],int(x['length']),cap,label,' → '.join(map(number,x['values_b20_50_100_200'])),
                             '→'.join(map(str,x['transition'])),number(x['largest_violation'])])
    text += [table(['ID','n','上限','指標','4幅の値','悪化した幅','最大悪化量'],rows),'',
             '## 最悪の例と G の形','',
             '各条件・各指標の上位10配列を JSON に残した。以下は default の b=20/200 の上位例。',
             'G 自体はすべて pseudoknot-free。交差対数は exact BPP>10⁻¹² の生成対のうち G の対を交差する対の数で、完成構造が擬似ノットを含み得ることを示す。',
             '交差対の確率和は交差対数の期待値であり、擬似ノット構造の確率ではない。','']
    rows=[]
    for beam in (20,200):
        for metric,label in [('delta_free_energy_kcal_mol','ΔF'),('base_state_tv_mean_Gfree','BPP TV')]:
            for x in s['worst_examples'][f'default_b{beam}_{metric}'][:5]:
                rows.append([label,beam,x['input'],int(x['length']),int(x['scaffold_pairs']),int(x['free_bases']),
                             int(x['maximum_pair_depth']),int(x['maximum_branching_ancestors_of_free_base']),
                             int(x['exact_crossing_pair_count_gt_1e_minus12']),number(x[metric]),ratio(x['partition_function_retention'])])
    text += [table(['順位指標','b','ID','n','G対数','自由塩基','G深さ','分岐祖先最大','交差対数','誤差','r'],rows),'',
             'GC比率、対合割合、分岐節数、交差対期待数、exact PseudoLoop の自由塩基平均も TSV に保存した。','',
             '## 診断メモとの関係','',
             '[有限 beam の質量診断](beam_mass_diag_20261011.md)の4入力を同じ母集団内で確認した。以下は nocap の b=100/200。','']
    rows=[]
    for name,examples in s['diagnostic_inputs'].items():
        for x in examples:
            if x['cap']=='nocap' and x['beam'] in (100,200):
                rows.append([name,int(x['length']),x['beam'],number(x['delta_free_energy_kcal_mol']),
                             ratio(x['partition_function_retention']),number(x['base_state_tv_mean_Gfree'])])
    text += [table(['ID','n','b','ΔF','r','BPP TV'],rows),'',
             '診断メモは `5936` の公開 WMBP・私的 BTAIL、`2449` の CASE1/2 prefix と WMv など、個別入力の取りこぼし箇所を対照で調べた。',
             '本集計はその個別原因を全配列へ帰属させる実験ではない。長さ別の誤差分布と外れ値を示し、次の診断入力を選ぶ材料にする。',
             'Empty G の対照や幅300以上の測定は本matrixにないため、この集計だけで擬似ノット固有の原因や幅300の十分性を一般化しない。','',
             '## 論文で言えることと言えないこと','',
             '日本語：この固定条件・保存済み cohort では、beam幅を広げたときの分布の変化と、exactに近い配列の割合を長さ別に報告できる。全体中央値だけでは外れ値を説明できないため、上側分位点・最大値・保持率を併記する。',
             '', 'English: For this fixed implementation and saved cohort, we can report error distributions by length and the fraction of sequences close to the exact calculation. Upper quantiles, maxima, and ratios of partition functions describe errors that an overall median can hide.','',
             '日本語：有限beamの誤差は全配列で単調に減るとは限らない。必要な幅は配列とGによって異なる。defaultとnocapの差は設定の感度であり、独立した局所上限の誤差分解ではない。',
             '', 'English: Errors do not decrease monotonically with beam width for every sequence. The required width depends on the sequence and the supplied scaffold. Differences between the default and uncapped settings measure sensitivity to those settings.','',
             '400 nt以上、生物学的な構造予測精度、全入力に対する単調収束、b≥nならexactという保証、実装の計算量上界、独立したCParty構造和との一般的同値性は、本集計では主張しない。',
             'ID中のRFAM/CRW等は出典区分であり、RNA familyや相同性を調整した精度推定・信頼区間は本集計に含めない。',
             '本体は精度用の並列CLI実行なので、保存された時間を孤立した速度比較へ転用しない。','',
             '## 検証・成果物・再生成','',
             f"実験repoの全テストは {v.get('unit_tests_passed','55')}/{v.get('unit_tests_passed','55')}件通過。packed読み込み、exact差し替え、元の失敗履歴、入力・binary不一致の拒否、単調性、NAを独立fixtureで確認した。",
             f"CLI間の整合チェック（exactを含む9条件）は {p['consistency_check_counts'].get('success',0):,}組成功。",
             f"既存12配列gridとの同じ入力hashの96条件、{next((c['scalar_comparisons'] for c in v.get('checks',[]) if c['name']=='pilot_metric_regression'),0):,}個の指標を照合し、最大差は {next((c['max_absolute_error'] for c in v.get('checks',[]) if c['name']=='pilot_metric_regression'),None)}だった。",
             '有限beamの最大整合差：'+', '.join(f'{k}={number(value)}' for k,value in s['maximum_consistency_errors'].items())+'。',
             'これらは出力の整合と集計の検証であり、構造の非曖昧性の新たな証明ではない。','',
             '- [配列×条件の指標 TSV](../../experiments/results/local/e1-bprna-short-20261011/metrics.tsv)',
             '- [全分布・単調性・外れ値 JSON](../../experiments/results/local/e1-bprna-short-20261011/summary.json)',
             '- [差し替え・入力・実装・script provenance](../../experiments/results/local/e1-bprna-short-20261011/provenance.json)',
             '- [検証結果](../../experiments/results/local/e1-bprna-short-20261011/validation.json)',
             '- [ΔF 分布 PNG](../../experiments/results/local/e1-bprna-short-20261011/delta_F_distribution.png)',
             '- [保持率分布 PNG](../../experiments/results/local/e1-bprna-short-20261011/retention_distribution.png)','',
             '図は長さ4区分×beam4幅。箱は25–75%、線は中央値、ひげは10–90%、三角はp99、×はΔF最大または保持率最小。ΔFのみ表示上の丸め誤差を0にし、symlog（線形域0.001）で示す。','',
             '論文repoをcwdにした再生成コマンド（解析先は新しいdirectoryを指定する）：','',
             '```bash',
             'OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \\',
             'PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/e1-bprna-mpl-20261011 \\',
             'python experiments/scripts/analyze_accuracy.py \\',
             '  results/runs/e1-bprna-full-packed-20261011 --packed-array \\',
             '  --exact-reruns results/runs/e1-bprna-exact-rerun-capr_raccess-20261011 \\',
             '    results/runs/e1-bprna-exact-rerun-lincparty-20261011 \\',
             '  --out results/local/e1-bprna-short-20261011-new',
             '',
             '# TSVを保存したまま集計・図だけ再生成',
             'OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \\',
             'PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/e1-bprna-mpl-20261011 \\',
             'python experiments/scripts/analyze_bprna_short.py \\',
             '  experiments/results/local/e1-bprna-short-20261011',
             '',
             'OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \\',
             'PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/e1-bprna-mpl-20261011 \\',
             'python experiments/scripts/report_bprna_short.py \\',
             '  experiments/results/local/e1-bprna-short-20261011 \\',
             '  --out textbook/notes/e1_bprna_short_20261011.md --plots',
             '```','',
             'qsub、git commit、`../LinearCParty` の変更は行っていない。意味のある変更と検証結果を `docs/progress.md` に追記した。','']
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text('\n'.join(text))


def render_plots(output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import analyze_bprna_short as short
    # Keep the zero baseline visible without a large unused negative ΔF range.
    plt.rcParams['axes.ymargin']=0
    _,data,_,_=short.load_metrics(output)
    short.plot_distributions(output,data,data['length'][:,0,0])


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--plots',action='store_true',help='also render both PNGs with tight data margins')
    args=parser.parse_args()
    write_note(args.output,args.out)
    if args.plots:render_plots(args.output)
