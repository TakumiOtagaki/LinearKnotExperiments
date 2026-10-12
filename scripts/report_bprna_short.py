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
    if 'capb_run' in p:
        return write_capb_note(output,destination,s,p,v)
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


def write_capb_note(output,destination,s,p,v):
    labels=('all','nt50_99','nt100_199','nt200_299','nt300_399')
    groups=[g for label in labels for g in s['bins'] if g['bin']==label]
    bins=groups[1:]
    title=lambda g:'全体' if g['bin']=='all' else f"{g['min_length']}–{g['max_length']}"
    condition=lambda g,cap,b:next(c for c in g['conditions'] if c['cap']==cap and c['beam']==b)
    adopted200=condition(groups[0],'adopted',200)
    long200=condition(bins[-1],'adopted',200)
    delta='delta_free_energy_kcal_mol';bpp='base_state_tv_mean_Gfree';context='context_tv_mean_Gfree'
    paired200={c['right']:c['metrics'][delta] for c in groups[0]['paired_cap_differences'] if c['beam']==200}
    nocap_differences={c['beam']:c['metrics'][delta] for c in groups[0]['paired_cap_differences'] if c['right']=='nocap'}
    note=['# bpRNA 400 nt 未満の E1：cap=max(b,30) の追加集計（2026-10-12）','',
          f"固定版 `15aa7742`、同じ **{s['inputs']:,}配列**・scaffold・モデル・exact 参照に、新 run の capb（b=50/100/200）を追加した。",
          '著者の決定（2026-10-12）に従い、局所上限 h=d=max(b,30) を採用条件とする。b=20 は既存 `beam20_default`、b≥50 は新しい `beam{b}_capb` を使う。',
          'ここで default は比較用の旧固定上限30/30、capb は max(b,30)、nocap は100000を意味する。defaultという mode 名を、新しく採用した既定と混同しない。','',
          f"採用条件の b=200 で ΔF≤0.01 kcal/mol は全体の {percent(1-adopted200['delta_F_above']['0.01']['fraction_eligible'])}。300–399 nt では {percent(long200['delta_F_above']['0.01']['fraction_eligible'])}が0.01を超え、{percent(long200['delta_F_above']['1']['fraction_eligible'])}が1を超える。",
          f"この区分の最大 ΔF は {number(long200['metrics'][delta]['max'])} kcal/mol、最小保持率は {percent(long200['metrics']['partition_function_retention']['min'])}。設定の変更だけでは有限 beam の大きな誤差は解消しない。",'',
          f"b200のΔFは旧defaultよりcapbで小さい配列が {paired200['default']['right_larger_count']:,}、大きい配列が {paired200['default']['left_larger_count']:,}。capbとnocapのΔFは {paired200['nocap']['equal_within_tolerance_count']:,}/{s['inputs']:,}配列で許容10⁻⁸ kcal/mol以内に一致した。詳細は配列内差の表に示す。",'',
          '## 入力・参照・保存方法','',
          table(['長さ nt','配列数'],[[title(g),f"{g['sequences']:,}"] for g in bins]),'',
          f"全体の {percent((bins[0]['sequences']+bins[1]['sequences'])/s['inputs'])}が200 nt未満。全体の割合と長さ別の上側分位点を併記する。",
          'G は保存済みの LinearPartition-V `b450fb3e`、beam200、dangles2、posterior≥0.7、stacked helix長≥2の scaffold。Turner04、37°C、既定 closure、BPP cutoff0、accessibility w=1/4/8/16、3 CLI を共通にした。',
          f"RT={s['rt_range'][0]:.10f} kcal/mol。新 run は {len(p['capb_tasks'])}タスク、{sum(p['capb_job_status_counts'].values()):,}/{sum(p['capb_job_status_counts'].values()):,} CLI 成功。元の短い区分は {len(p['tasks'])}タスク。",
          'exact は前回と同じ本体および指定の再実行。元の失敗349件（CapR/Raccess各172、LinearCParty 5、重複を除く172配列）を成功した32 GiB再実行で置換し、元の失敗履歴を保持した。',
          f"exact差し替え後の統合 CLI 観測は {p['effective_job_status_counts']['success']:,}件成功、指標行は {p['analysis_status_counts']['success']:,}行（RNA×実測11条件）。b20 capbや採用条件の再表示は新しい測定として数えない。",
          '入力hash・S/G・source commit・モデル設定・parameter hash・binary hashes・計算引数を照合した。前回のTSVを再利用する前に、全タスクのmanifest/checkpointのhashと前回の検証結果を確認した。',
          'tarは展開せず、member streamから読む。新しいscalar値と時間・RSSは /tmp の一つのSQLiteへ一時保存し、配列ごとの11条件を一つのTSVへ統合した。共有領域の出力はTSV 1、JSON 3、PNG 2の6ファイル。','',
          '## 指標と集計','',
          'ΔF=F_b−F_exact、F=−RT log Z。保持率r=Z_b/Z_exact=exp(−ΔF/RT)。Raccessのfull-precision metadataを使い、微小な負のΔFもそのまま保存した。',
          'BPP TVはGの自由塩基における相手塩基・unpaired状態の平均TV。Context TVは同じ自由塩基の7 contextの平均TV。Generated-pair L1、unpaired MAE、context別MAEも保存した。',
          'AccessibilityはGの対合塩基を含まないwindowを評価する。Probability MAEと最大誤差、両確率が正のwindowでのopening-energy MAE、片側だけ0の件数を分ける。',
          'まずRNAごとの指標を計算し、RNAを等重みで集計する。NAを0にせず、有効件数を保持する。分位点は線形補間。JSONには全体・4区分・25 nt区分、中央値・平均・p90・p99・最大・最小を保存した。',
          'このRTでは '+', '.join(f'ΔF={f:g} → r={percent(math.exp(-f/s["rt_range"][0]))}' for f in (.01,.1,1))+'。','',
          '## 3上限条件の ΔF 分布と閾値超過','',
          '単位はkcal/mol。閾値は厳密な >。b20 capbはdefaultの同じ観測である。','']
    rows=[]
    for g in bins:
        for c in g['conditions']:
            if c['cap']=='adopted':continue
            m=c['metrics'][delta]
            rows.append([title(g),c['cap'],c['beam']]+[number(m[k]) for k in ('median','mean','p90','p99','max')]+
                        [percent(c['delta_F_above'][k]['fraction_eligible']) for k in ('0.01','0.1','1')])
    note += [table(['nt','上限','b','中央値','平均','p90','p99','最大','>0.01','>0.1','>1'],rows),'',
             '## 採用条件を一つの列にまとめた比較','',
             '各セルは「ΔF p99 / ΔF>0.01 の割合」。採用条件はb20でdefault、b≥50でcapbを選ぶ。','']
    rows=[]
    for g in groups:
        for b in (20,50,100,200):
            rows.append([title(g),b]+[number(condition(g,cap,b)['metrics'][delta]['p99'])+' / '+percent(condition(g,cap,b)['delta_F_above']['0.01']['fraction_eligible']) for cap in ('default','nocap','adopted')])
    note += [table(['nt','b','旧default 30/30','nocap 100000','採用条件 max(b,30)'],rows),'',
             '## 採用条件の保持率・BPP・context・accessibility（全幅）','',
             'TVとaccessibility列はRNA間p99。保持率は平均と最小。Accessibilityは各RNA内のprobability MAEを集計した。','']
    rows=[]
    for g in bins:
        for c in g['conditions']:
            if c['cap']!='adopted':continue
            m=c['metrics'];r=m['partition_function_retention']
            rows.append([title(g),c['beam'],ratio(r['mean']),ratio(r['min']),number(m[bpp]['p99']),number(m[context]['p99'])]+
                        [number(m[f'access_w{w}_probability_mae_eligible']['p99']) for w in (1,4,8,16)])
    note += [table(['nt','b','r平均','r最小','BPP TV','context TV','w1 MAE','w4 MAE','w8 MAE','w16 MAE'],rows),'',
             '## b=200 の保持率と周辺確率（3条件）','']
    rows=[]
    for g in bins:
        for cap in ('default','capb','nocap'):
            m=condition(g,cap,200)['metrics'];r=m['partition_function_retention']
            rows.append([title(g),cap,ratio(r['median']),ratio(r['mean']),ratio(r['min']),number(m[bpp]['median']),number(m[bpp]['p99']),number(m[bpp]['max']),number(m[context]['p99']),number(m[context]['max'])])
    note += [table(['nt','上限','r中央値','r平均','r最小','BPP TV中央値','BPP TV p99','BPP TV最大','context TV p99','context TV最大'],rows),'',
             '## b=200 の accessibility（3条件）','',
             'Probability MAEのRNA間p90。全幅の分布・opening-energy誤差・window件数はJSON/TSVにある。','']
    rows=[]
    for g in bins:
        for cap in ('default','capb','nocap'):
            m=condition(g,cap,200)['metrics']
            rows.append([title(g),cap]+[number(m[f'access_w{w}_probability_mae_eligible']['p90']) for w in (1,4,8,16)])
    note += [table(['nt','上限','w1','w4','w8','w16'],rows),'',
             'Eligible windowが0の配列数：'+', '.join(f'w={w}: {n:,}' for w,n in s['eligible_zero_window_sequences'].items())+'。','',
             '## capb と nocap・旧default の配列内差','',
             '以下は同じRNA・同じbの誤差指標について capb − 比較条件。負はcapbの誤差が小さい。TV列もexactに対するTV誤差の差であり、2設定のposterior同士のTVではない。有限beamで上限を変えると保持集合も変わるため、上限だけの損失分解とは呼ばない。',
             '大小の判定許容はΔF 10⁻⁸ kcal/mol、TV 10⁻⁹。JSONは全指標の差・有効件数・両方向の件数・許容内同値の件数を保存する。','']
    rows=[]
    for g in bins:
        for c in g['paired_cap_differences']:
            m=c['metrics'][delta]
            rows.append([title(g),c['right'],c['beam']]+[number(m[k]) for k in ('mean','p99','max','min')]+
                        [m['left_larger_count'],m['right_larger_count'],number(c['metrics'][bpp]['p99']),number(c['metrics'][context]['p99'])])
    note += [table(['nt','比較条件','b','ΔΔF平均','p99','最大','最小','capb誤差大','比較条件誤差大','ΔBPP TV p99','Δcontext TV p99'],rows),'',
             '全体の配列内差：','']
    rows=[]
    for c in groups[0]['paired_cap_differences']:
        m=c['metrics'][delta]
        rows.append([c['right'],c['beam'],number(m['mean']),m['left_larger_count'],m['right_larger_count'],m['equal_within_tolerance_count']])
    note += [table(['比較条件','b','ΔΔF平均','capb誤差大','比較条件誤差大','許容内同値'],rows),'',
             f"b50ではcapbのΔFがnocapより大きい配列が {nocap_differences[50]['left_larger_count']}、小さい配列が {nocap_differences[50]['right_larger_count']}あり、最大差は {number(nocap_differences[50]['max'])} kcal/mol。b100で許容を超える差は {nocap_differences[100]['left_larger_count']}配列、最大差 {number(nocap_differences[100]['max'])} kcal/mol。b200では全配列のΔF差が0だった。幅ごとの観測として扱い、400 nt以上へ一致を外挿しない。",'',
             '## 時間・peak RSS の capb/default 比（参考値、LinearCParty）','',
             '同じRNA・同じbの個別比を計算してRNA間で集計した。時間はelapsed_seconds（GNU time wrapperを含む経過時間）。GNU wall timeの比もJSONと表に保存する。',
             'これらは精度用の並列実行、各条件1回、別日・別割当の値。短い配列では起動・出力・スケジューリングの影響が大きい。孤立した速度比較、速度向上、計算量の根拠には使わない。',
             '各比は有限の非負分子と正の分母を条件とし、分母0や欠測を除外して有効Nを示す。GNU wall timeは0.01秒に丸められる。RSSは個別プロセスのpeak。b20は同じdefault測定の共用なので比1は新測定の結果ではない。','']
    rows=[]
    for g in bins:
        for c in g['resource_ratios_capb_over_default']:
            m=c['metrics'];t=m['LinearCParty_elapsed_seconds'];r=m['LinearCParty_peak_rss_kib'];w=m['LinearCParty_time_wall_seconds']
            rows.append([title(g),c['beam']]+[number(t[k]) for k in ('median','p90','p99')]+[t['n']]+
                        [number(r[k]) for k in ('median','p90','p99')]+[r['n'],number(w['median']),w['invalid_or_zero_denominator']])
    note += [table(['nt','b','時間中央値','時間p90','時間p99','有効N','RSS中央値','RSSp90','RSSp99','有効N','GNU比中央値','GNU除外N'],rows),'',
             '## 採用条件の beam 幅に対する単調性','',
             'b=20→50→100→200の3段階すべてでΔF・BPP TV・context TVの非増加、保持率の非減少を確認した。ΔF許容10⁻⁸ kcal/mol、r許容10⁻⁸、TV許容10⁻⁹。',
             '採用条件では幅と局所上限の両方が変わる。この検査は採用する方法全体のbeam-width sensitivityであり、beamだけの効果とは呼ばない。ΔFとrの判定は数学的には同値でも、絶対許容の指数変換により割合が異なり得る。','']
    rows=[]
    for g in groups:
        stats={x['metric']:x for x in s['monotonicity'] if x['bin']==g['bin'] and x['cap']=='adopted'}
        rows.append([title(g),g['sequences']]+[percent(stats[k]['fraction_eligible']) for k in (delta,'partition_function_retention',bpp,context)]+
                    [stats[k]['nonmonotone_sequences'] for k in (delta,bpp,context)])
    note += [table(['nt','N','ΔF↓','r↑','BPP TV↓','context TV↓','ΔF違反N','BPP違反N','context違反N'],rows),'',
             '採用条件の大きな単調性違反（値はb20/50/100/200）：','']
    rows=[]
    for metric,label in ((delta,'ΔF'),(bpp,'BPP TV'),(context,'context TV')):
        for x in s['nonmonotone_examples']['adopted_'+metric][:3]:
            rows.append([x['input'],int(x['length']),label,' → '.join(map(number,x['values_b20_50_100_200'])),'→'.join(map(str,x['transition'])),number(x['largest_violation'])])
    note += [table(['ID','nt','指標','4幅の値','悪化した幅','悪化量'],rows),'',
             '## 採用条件 b200 の外れ値','']
    rows=[]
    for x in s['worst_examples']['adopted_b200_'+delta][:5]:
        rows.append([x['input'],int(x['length']),number(x[delta]),ratio(x['partition_function_retention']),number(x[bpp]),number(x[context]),int(x['scaffold_pairs']),int(x['free_bases'])])
    note += [table(['ID','nt','ΔF','r','BPP TV','context TV','G対数','自由塩基'],rows),'',
             '全条件の上位10例、単調性違反例、[質量診断](beam_mass_diag_20261011.md)の4入力はJSONに保存した。個別の取りこぼし箇所を全配列へ帰属させる実験は行っていない。','',
             '## 論文で言えること・言えないこと','',
             '日本語：保存済みの400 nt未満のcohortについて、新しい既定 max(b,30) の exact に対する誤差分布と保持率を、長さ別に報告できる。旧上限30/30と100000との配列内差から、固定条件下での設定の感度を報告できる。',
             '', 'English: For the saved cohort below 400 nt, we can report error distributions and partition-function ratios for the adopted local caps max(b,30). Paired differences from caps of 30/30 and 100000 describe sensitivity to these settings under the fixed experimental conditions.','',
             '日本語：採用条件で幅を広げても、すべての配列で誤差が単調に減るとは限らない。全体は短い配列が多いため、全体割合だけで長い配列の精度を説明せず、長さ別の上側分位点・外れ値を示す。',
             '', 'English: Errors under the adopted settings do not decrease monotonically with beam width for every sequence. Because short sequences dominate this cohort, length-specific upper quantiles and individual large errors are needed alongside the overall fractions.','',
             'この集計は注釈に対する構造予測精度、上限を緩めた場合の単調改善、全入力での収束、b≥nならexactという保証、実装の計算量の証明、独立したCParty構造和との一般的同値性を示さない。ID中の出典区分は相同性を調整した精度推定や信頼区間を意味しない。',
             '時間・RSS比は参考値で、正式な性能の主張には逐次・同一環境・反復測定が必要。著者の既定の決定は本集計の結果による新たな最適化選択とは区別する。','',
             '## 400 nt 以上の区分に向けた含意','',
             '300–399 ntに残る大きな誤差と単調性違反から、b200やmax(b,30)が400 nt以上にも十分とは外挿できない。長い区分では同じ採用規則を明記し、幅ごとの完了・timeout・メモリ不足を含めて保存する。',
             'exactを計算していない長い入力では、最大測定幅との差をexact誤差や保持率と呼ばない。幅を広げたときの出力変化として報告し、実行可能な境界付近のexact対照を別に用意する。',
             'bとともに局所上限も広がるため、精度だけでなく実行可能性・資源使用量を幅ごとに確認する。必要に応じて固定上限・nocapの対照を同じS/Gで残す。b200を十分と仮定せず、外れ値の診断とb300以上の対照を候補にする。',
             '本作業は既存runの集計のみ。新しい投入・幅300以上の測定・実装変更は行っていない。','',
             '## 開発上の検証・成果物・再生成','',
             f"全 {v['unit_tests_passed']}テスト成功。{sum(c['passed'] for c in v['checks'])}/{len(v['checks'])}成果物チェック成功。既存197,000行の指標は全列そのまま保持した。",
             f"統合した実測条件のCLI間整合は {p['consistency_check_counts']['success']:,}組、新capbは {p['capb_consistency_check_counts']['success']:,}組。exact {p['rechecked_exact_consistency_counts']['success']:,}組も再確認した（統合件数へ重複計上しない）。",
             '最大整合差：'+', '.join(f'{k}={number(value)}' for k,value in s['maximum_consistency_errors'].items())+'。',
             '既存12配列gridの96条件・4,416指標も回帰照合した。これらは集計と出力整合の検証であり、科学的なexactnessの新たな証明とは区別する。',
             '', '集計プロセスの資源使用量（CLI計測値ではない）：','', '```text',v['analysis_resources'],'```','']
    base='../../experiments/results/local/'+output.name+'/'
    note += ['- '+f'[{label}]({base}{filename})' for label,filename in (
        ('配列×実測条件 TSV','metrics.tsv'),('全分布・配列内差・資源比・単調性 JSON','summary.json'),
        ('参照と結合の provenance','provenance.json'),('検証結果','validation.json'),
        ('ΔF分布 PNG','delta_F_distribution.png'),('保持率分布 PNG','retention_distribution.png'))]
    note += ['', '図は3上限条件・4長さ区分・4幅。capb b20はdefaultと同じ値。箱25–75%、中央値、ひげ10–90%、三角p99、×最大ΔFまたは最小r。微小負ΔFは表示だけ0とし、TSV/JSONに符号を残す。','',
             '論文repoをcwdとする再生成コマンド（出力先は新しいdirectoryにする）：','', '```bash',
             'export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1',
             'export PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/e1-capb-mpl-20261012',
             '/usr/bin/time -v -o /tmp/e1-capb-resources-20261012.log \\',
             '  python experiments/scripts/analyze_bprna_short.py \\',
             '  experiments/results/local/e1-bprna-short-capb-20261012-new \\',
             '  --extend experiments/results/local/e1-bprna-short-20261011 \\',
             '  --capb-run experiments/results/runs/e1-bprna-capb-short-20261012',
             'python -m unittest discover -s experiments/tests -v > /tmp/e1-capb-tests-20261012.log 2>&1',
             'python experiments/scripts/validate_bprna_short.py \\',
             '  experiments/results/local/e1-bprna-short-capb-20261012-new \\',
             '  --pilot experiments/results/local/e1-bprna-15aa7742-20261011-accuracy.json \\',
             '  --tests-log /tmp/e1-capb-tests-20261012.log \\',
             '  --resource-log /tmp/e1-capb-resources-20261012.log',
             'python experiments/scripts/report_bprna_short.py \\',
             '  experiments/results/local/e1-bprna-short-capb-20261012-new \\',
             '  --out textbook/notes/e1_bprna_short_capb_20261012.md',
             '```','',
             'TSVを保存したまま集計・図だけ再生成する場合は `python experiments/scripts/analyze_bprna_short.py <output>` を使い、validationとnoteも再生成する。',
             'qsub・git commit・調整中の `../LinearCParty` の変更は行っていない。','']
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text('\n'.join(note))


def render_plots(output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import analyze_bprna_short as short
    # Keep the zero baseline visible without a large unused negative ΔF range.
    plt.rcParams['axes.ymargin']=0
    provenance,data,_,_=short.load_metrics(output)
    short.plot_distributions(output,data,data['length'][:,0,0],provenance.get('caps',short.CAPS))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--plots',action='store_true',help='also render both PNGs with tight data margins')
    args=parser.parse_args()
    write_note(args.output,args.out)
    if args.plots:render_plots(args.output)
