#!/usr/bin/env python3
"""Write reviewable SI proposals from the calculated C/D tables only."""
import argparse
import csv
import json
from pathlib import Path

from run_cli import digest, now, write_json

ROOT = Path(__file__).resolve().parents[1]
LOCAL = ROOT/'results/local/e1-local-errors-89baea52-20261019'
RECOVERY = ROOT/'results/local/pk-recovery-pairset-89baea52-20261019'
LABELS = {'bpp_all': 'All BPP', 'bpp_crossing_G': 'BPP crossing G',
          'pseudoloop': 'PseudoLoop', 'access_w4_total': '4-nt total accessibility'}


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def pct(a, b):
    return f'{100*int(a)/int(b):.2f}%'


def recovery():
    rows = read(RECOVERY/'summary.tsv')
    keys = {(r['cohort'], r['parameter'], r['definition']): r for r in rows}
    lines = []
    for cohort, parameter in [('200', 'turner04'), ('200', 'dp09'), ('80', 'turner04'), ('rrna', 'turner04')]:
        for definition in ['old_symbols', 'crossing_G']:
            r = keys[cohort, parameter, definition]
            n, predicted, matched = r['reference_pairs'], r['mfe_predicted_pairs'], r['mfe_matched_pairs']
            lines.append(f"| {cohort} | {parameter} | {definition} | {n} | {predicted} | {matched} | {pct(matched,n)} | {pct(matched,predicted)} | {r['reference_pairs_gt_half']}/{n} ({pct(r['reference_pairs_gt_half'],n)}) | {r['rnas_mean_probability_gt_half']}/{r['rnas']} |")
    audit = read(RECOVERY/'reference_audit.tsv')
    audits = []
    for cohort in ('200', '80', 'rrna'):
        selected = [r for r in audit if r['cohort'] == cohort and r['parameter'] == 'turner04']
        sums = {k:sum(int(r[k]) for r in selected) for k in
                ('canonical_reference_pairs','scaffold_pairs','noncanonical_reference_pairs_removed',
                 'canonical_nonscaffold_not_crossing','reference_pairs_added','reference_pairs_removed')}
        audits.append(f"| {cohort} | {sums['canonical_reference_pairs']} | {sums['scaffold_pairs']} | {sums['noncanonical_reference_pairs_removed']} | {sums['canonical_nonscaffold_not_crossing']} | {sums['reference_pairs_added']} | {sums['reference_pairs_removed']} |")
    distances = []
    for r in read(RECOVERY/'span_summary.tsv'):
        if r['cohort']=='200' and r['parameter']=='turner04' and r['definition']=='crossing_G':
            distances.append(f"| {r['span']} | {r['reference_pairs']} | {r['reference_pairs_gt_half']} ({pct(r['reference_pairs_gt_half'],r['reference_pairs'])}) | {float(r['mean_probability']):.6f} |")
    text = '''# 解析 D: 座標集合と G との交差による pseudoknot 回復率

2026-10-19。原稿は変更せず、保存済み 89baea52 出力だけを再集計した。
新しい RNA CLI の実行、tar の展開、外部実装の変更、git commit は行っていない。

## 定義と出典

塩基対を 1-based の座標 (i,j), i<j の集合として読む。
(i,j) と (k,l) は i<k<j<l または k<i<l<j のとき交差する。
固定 scaffold の対集合を G、canonical な参照対の集合を A、全 MFE 対の集合を M とする。
Canonical は AU/UA/CG/GC/GU/UG。

- 新参照集合 T = {p ∈ A \\ G : G の少なくとも一対と交差する}。
- 新予測集合 P = {p ∈ M \\ G : G の少なくとも一対と交差する}。
- sensitivity = Σ|T ∩ M| / Σ|T|、PPV = Σ|T ∩ P| / Σ|P|。
  全出力で T∩M = T∩P を検査した。集計は RNA 間の分子・分母を合計する micro average。
- 確率 >0.5 は厳密な不等号。RNA 内平均は新参照集合 T 上で取り直す。
  RNA 間平均は等重み。距離群は RNA 内の参照対の j−i の中央値で分ける。
- 括弧種は座標を読むためだけに使い、参照・予測の選択には使わない。
  全追加対を pseudoknot とみなす定義とは区別する。

入力・G・出力は `results/runs/pk-recovery-{200-turner04,200-dp09,80,rrna}-89baea52-20261018/`。
200 本は exact、80 本と rRNA は b=h=d=400 / outer ranking。
既存値は `results/local/small-rerun-89baea52-20261018/{pk_results.json,mid300_700_results.json,recovery_summary.json,rrna_artifact_index.json}`。
この 89baea52 集計ディレクトリには README はなかったため、`scripts/analyze_small_rerun_89baea52.py` を確認した。
参照は `results/local/pk-recovery-bprna-20261016/pk_sel.json` と `results/local/pk-recovery-long-rrna-20261016/mid300_700_selection.json`。
rRNA は同じ旧ディレクトリの README と、89baea52 の `rrna_reference_metadata.json` を確認した。
旧 README の「全24対がゼロ」は旧実装の記録であり、新版の比較値に流用していない。

## 旧定義との比較／SI の回復表の差し替え案

`crossing_G` の行を採用する。200 本の2行で既存の回復表を差し替え、80 本と rRNA は別 panel に置く案。
表 caption には上記 T/P と micro average、exact と有限幅の違いを明記する。

| RNA 群 | モデル | 定義 | 参照対 | 予測対 | 一致対 | Sensitivity | PPV | 参照 p>0.5 | RNA 内平均 p>0.5 |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
''' + '\n'.join(lines) + '''

全483出力で予測 `[]` 集合と新予測集合 P は**座標ごとに完全一致**した。
予測側の符号化に由来する差は0で、PPV と一致対数は全条件で変わらない。
参照側には次の差がある。

- 200 本では、23 RNA の canonical な `{}` 対57対が G と交差するため加わる（983→1040）。
  Turner04/DP09 とも、この57対には MFE 一致も p>0.5 もない。
  そのため sensitivity は 21.36→20.19%、11.19→10.58%、p>0.5 の割合は18.21→17.21%、8.14→7.69%となる。
  RNA 内平均 p>0.5 の件数34/13は変わらない。
- 80 本では `bpRNA_CRW_4958` (133,166) と `bpRNA_CRW_4764` (242,275) が除かれる（811→809）。
  元注釈でこれらと交差する `()` 対はそれぞれ (135,167), (244,276) の UU 対だけで、G の構築時に除かれている。
  これら二対は新しい G に対する crossing pair ではない。
  両対とも MFE 一致でも p>0.5 でもないため、分子40/24は変わらない。
  sensitivity は4.932→4.944%、p>0.5 の割合は2.959→2.967%。
- rRNA の保存済み canonical 参照24対はすべて G と交差し、集合も指標も変わらない。
  MFE 一致0/24、予測14対、p>0.5 は0/24、報告値ゼロ9/24、最大確率0.00045326084596248895。

### 参照から除かれる対の内訳

| 群 | canonical 参照対 | G の対 | 非canonical除外 | G外で非交差 | 旧集合に追加 | 旧集合から除外 |
|---|---:|---:|---:|---:|---:|---:|
'''+'\n'.join(audits)+'''

200/80 本の canonical 参照は元の全括弧種を読む。元注釈の対数は canonical と非canonical の和。
G に残す条件は従来どおり canonical な `()` 対かつ j−i≥4。
Canonical でも G 外で非交差の対は、新しい pseudoknot 参照集合には含めない。

rRNA の利用可能な保存物には元の完全な DBN がなく、G と canonical な `[]`/`{}` 参照座標24対がある。
したがって rRNA の新集合は、この凍結された参照座標と G に同じ交差判定を適用した結果である。
上表の rRNA の非canonical除外数は保存済み metadata の合計で、元DBNから再分類した数ではない。
元DBNの未保存の括弧種を含めた完全な再監査とは主張しない。
`[]`/`{}` を参照選定に使い分ける方針は廃し、保存された候補座標について同じ幾何学的判定を使う。

### 確率・距離の他の指標

200 本 Turner04 の距離表も分母と平均を更新する案（全条件は `span_summary.tsv`）。

| j−i (nt) | 参照対 | p>0.5 | 対ごとの平均確率 |
|---|---:|---:|---:|
'''+'\n'.join(distances)+'''

RNA 内の参照距離中央値 ≤100 / >100 の群数は、200 本では174/26、80 本では37/43のまま。
群内の等重み RNA 平均確率は、200 本 Turner04 で0.18995349→0.18549277 / 0.02350847（不変）、
DP09 で0.08843006→0.08478907 / 0.00015225（不変）、
80 本で0.08624786→0.08333152 / 0.00488461（不変）。
RNA平均に使う対集合を変更したためで、各対の確率を再計算した差ではない。
今回の1040対という新しい分母を、旧983対に対する別の到達可能性検査の分母へ自動で移さない。

## 英文案（SI、3文）と対応する日本語

Reference pairs are canonical annotation pairs that cross at least one scaffold pair in G, and predicted pairs are new MFE pairs that cross G.
We matched pairs by nucleotide coordinates and obtained sensitivities of 210/1,040 (20.19%) and 110/1,040 (10.58%) with Turner 2004 and DP09, respectively; PPVs were 210/850 (24.71%) and 110/413 (26.63%).
The fractions of reference pairs with probability above 0.5 were 179/1,040 (17.21%) and 80/1,040 (7.69%), respectively.

参照集合は G の少なくとも一対と交差する canonical な注釈対、予測集合は G と交差する新規 MFE 対とする。
座標で照合すると、Turner04/DP09 の sensitivity は210/1040（20.19%）/110/1040（10.58%）、PPV は210/850（24.71%）/110/413（26.63%）である。
参照対のうち確率0.5を超える割合は179/1040（17.21%）/80/1040（7.69%）である。

## 再現と検証

```bash
OPENBLAS_NUM_THREADS=1 python3 -B scripts/analyze_pk_recovery_pairset_89baea52.py \\
  --out results/local/pk-recovery-pairset-89baea52-20261019
python3 -B scripts/report_local_errors_pairset_89baea52.py recovery
OPENBLAS_NUM_THREADS=1 python3 -B -m unittest discover -s tests -p test_local_errors_and_pairset_recovery.py -v
```

解析の出力先は既存ディレクトリを上書きしないため、再集計時は別の出力先を使う。
`rna_metrics.tsv` は483条件、`pair_metrics.tsv` は新旧参照の和集合、`definition_changes.tsv` は全変更座標、
`reference_audit.tsv` は参照除外数、`summary.{tsv,json}` は集計、`span_summary.tsv` は距離別の値。
`provenance.json` は入力・plan・task manifest/checkpoint・解析コード・出力のSHA-256を記録する。
全480旧MFE集計と全483旧確率ベクトルを保存済み89baea52集計に完全一致で照合した。
raw artifact の checksum、保存コマンド、入力 S/G、MFE の G 保持、全参照について T∩M=T∩P を検査した。
境界0.1/0.5、疎なBPPの欠落ゼロ、括弧色の交換、交差判定、非canonical除去を6 testsで検査した。
科学的な集合・指標の変更と、これらの開発上の検査を分けて扱う。
'''
    (RECOVERY/'README.md').write_text(text)


def local_errors():
    rows = read(LOCAL/'summary.tsv')
    keys = {(r['tier'], int(r['beam']), r['metric']): r for r in rows}
    tables = []
    for metric, label in LABELS.items():
        lines = [f'### {label}', '', '| 長さ nt | b | RNA 数 | 取りこぼし対/塩基/窓 | 該当RNA n (%) | 作り出し対/塩基/窓 | 該当RNA n (%) |',
                 '|---|---:|---:|---:|---:|---:|---:|']
        for tier, length in [('nt50_99','50–99'),('nt100_199','100–199'),('nt200_299','200–299'),('nt300_399','300–399')]:
            for beam in (20,50,100,200):
                r = keys[tier,beam,metric]
                lines.append(f"| {length} | {beam} | {r['rnas']} | {r['missed_events']} | {r['missed_rnas']} ({pct(r['missed_rnas'],r['rnas'])}) | {r['created_events']} | {r['created_rnas']} ({pct(r['created_rnas'],r['rnas'])}) |")
        tables.append('\n'.join(lines))
    overall = []
    for beam in (20,50,100,200):
        for metric, label in LABELS.items():
            r = keys['all',beam,metric]
            overall.append(f"| {beam} | {label} | {r['missed_events']} | {r['missed_rnas']} ({pct(r['missed_rnas'],r['rnas'])}) | {r['created_events']} | {r['created_rnas']} ({pct(r['created_rnas'],r['rnas'])}) | {r['either_rnas']} ({pct(r['either_rnas'],r['rnas'])}) |")
    b200 = [keys['nt300_399',200,m] for m in LABELS]
    english = ("Small average errors can coexist with large local probability errors in either direction.\n"
        f"Among 1,095 RNAs of 300–399 nt at width 200, {b200[1]['missed_rnas']} ({pct(b200[1]['missed_rnas'],1095)}) had at least one G-crossing pair with exact probability at least 0.5 and finite-width probability below 0.1, while {b200[1]['created_rnas']} ({pct(b200[1]['created_rnas'],1095)}) had a pair with the reverse thresholds.\n"
        f"The corresponding RNA counts were {b200[2]['missed_rnas']} and {b200[2]['created_rnas']} for PseudoLoop probabilities, and {b200[3]['missed_rnas']} and {b200[3]['created_rnas']} for 4-nt total accessibility.")
    old = read(ROOT/'results/local/e1-pk-stratified-89baea52-20261018/local_error_quantiles.tsv')
    examples = read(LOCAL/'examples.tsv')
    created_window = next(r for r in examples if r['tier']=='all' and r['beam']=='200'
                          and r['metric']=='access_w4_total' and r['direction']=='created')
    comparisons=[]
    for metric, old_metric in [('bpp_crossing_G','crossing_bpp_ge0.1'),('pseudoloop','pseudoloop_ge0.1'),('access_w4_total','access_w4_ge0.1')]:
        r=next(r for r in old if r['axis']=='all' and r['stratum']=='all' and int(r['beam'])==200 and r['metric']==old_metric)
        new=keys['all',200,metric]
        comparisons.append(f"| {LABELS[metric]} | {float(r['max']):.9g} | {float(new['max_error']):.9g} |")
    text='''# 解析 C: 局所確率誤差の両側

2026-10-19。保存済み89baea52 E1 short の全24,625 RNA、exact と b=20/50/100/200を使った。
原稿・外部実装の変更、RNA CLI の実行、tar の展開、git commit は行っていない。

## 事象と分母

- 取りこぼし: exact ≥0.5 かつ有限幅 <0.1。
- 作り出し: exact <0.1 かつ有限幅 ≥0.5。
- BPP は全座標対と、G の少なくとも一対と交差する追加対を別に集計する。
  交差は i<k<j<l または k<i<l<j。複数のG対に交差しても一対一事象。
- PseudoLoop は塩基ごとの確率。Accessibility は4塩基が同時に非対合である total probability。
- 各表のRNA率の分母は長さ群の全RNA。評価域に事象を持たないRNAも分母に残す。
  50–99/100–199/200–299/300–399 nt のRNA数は12,428/9,656/1,446/1,095。
- BPP は cutoff 0 の両出力の座標和集合で比較し、欠落対だけをゼロとする。
  両出力にない対は両側ゼロなので、どちらの事象にもならない。
  PseudoLoop は全塩基、4-nt窓は全開始位置を読み、欠行をゼロで埋めない。
  Scaffold対、scaffold対合塩基、scaffold対合塩基を含む窓も全域に含まれるが、これらに閾値反転は生じない。
- 事象数は座標の延べ数。隣接塩基・重複窓・同じRNAの対を独立反復とはみなさない。
- exact は局所制限なし。有限幅は outer ranking / h=d=max(b,30) なので、beam と局所制限を合わせた設定の誤差。
  これは同じ S/G/モデルの数値近似誤差であり、参照注釈に対する構造予測の偽陽性・偽陰性ではない。

## SI 表案

長さ群×幅を行、取りこぼしと作り出しの事象数／RNA数・割合を並列の列とし、4種類の事象をpanelに分ける。
既存の平均誤差・exact ≥0.1に限定した表を補う表とする。

'''+ '\n\n'.join(tables)+'''

### 全24,625 RNA の集計

| b | 事象 | 取りこぼし数 | 該当RNA n (%) | 作り出し数 | 該当RNA n (%) | いずれかのRNA n (%) |
|---:|---|---:|---:|---:|---:|---:|
'''+ '\n'.join(overall)+'''

## 既存集計との違い

旧 `results/local/e1-pk-stratified-89baea52-20261018/` にREADMEはなかったため、
同版の `provenance.json`、`scripts/analyze_e1_pk_stratified.py` と、説明を持つ旧日付のREADMEを読んだ。
旧局所表は exact ≥0.1 だけを選択するので、今回の「作り出し」は定義上すべて評価域外だった。
今回の表はexact側の確率で座標を限定せず、両側を評価する。
既存の平均誤差を変更した結果ではない。

参考として b=200 の局所最大絶対誤差も比較する。
旧値はexact ≥0.1に限定、新値は全座標（対応する事象域）である。

| 事象 | 旧最大 | 全域での最大 |
|---|---:|---:|
'''+ '\n'.join(comparisons)+f'''

追加で見つかった最大差は {created_window['input']}（{created_window['length']} nt）の窓
{created_window['coordinate']}–{int(created_window['coordinate'])+3}で、exact は{float(created_window['exact_probability']):.9g}、
有限幅は{float(created_window['beam_probability']):.9g}だった。
この窓のexact確率は0.1未満なので、旧局所表の評価域には含まれなかった。
b=200の取りこぼし・作り出しは、4事象ともすべて300–399 nt群にあり、他の3長さ群では0件だった。

`summary.tsv` には局所最大誤差 >0.1 / >0.25 のRNA数と率、およびRNA内最大誤差のp95/p99/最大も収めた。
これらは上の強い閾値反転だけに限定しない指標。分位点はNumPyの線形補間。

## 英文案（3文）と対応する日本語

'''+english+f'''

RNA内の平均誤差が小さくても、局所の確率には両方向の大きな誤差が残る。
300–399 ntの1,095本をb=200で計算したとき、Gと交差する対の取りこぼしを持つRNAは{b200[1]['missed_rnas']}本（{pct(b200[1]['missed_rnas'],1095)}）、作り出しを持つRNAは{b200[1]['created_rnas']}本（{pct(b200[1]['created_rnas'],1095)}）だった。
同じRNA数はPseudoLoopでは{b200[2]['missed_rnas']}/{b200[2]['created_rnas']}本、4-nt total accessibilityでは{b200[3]['missed_rnas']}/{b200[3]['created_rnas']}本だった。

## 再現・出力・検証

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B -u scripts/analyze_e1_local_errors_89baea52.py \\
  --run results/runs/e1-bprna-short-89baea52-20261018 \\
  --out results/local/e1-local-errors-89baea52-20261019 --workers 4
python3 -B scripts/report_local_errors_pairset_89baea52.py local
OPENBLAS_NUM_THREADS=1 python3 -B -m unittest discover -s tests -p test_local_errors_and_pairset_recovery.py -v
```

再集計時は別の出力先を使う（解析scriptは既存directoryを上書きしない）。
`summary.tsv` は長さ群・全体×幅×事象の80行。
`rna_counts.tsv.gz` はRNA×幅×事象の394,000行、`examples.tsv` は各群・幅・方向・事象の最大差の例。
塩基・窓開始位置・対座標は1-based、4-nt窓の終端は開始+3。
`provenance.json` は全taskのhash、成功数、CLI間整合、コード・出力hashを記録する。
全369,375 CLI jobの保存成功と123,125組のCLI間整合を確認し、入力・コマンド・raw artifact checksumを検査した。
RNA CLIや既存の数値を生成し直す実装テストは実行していない。
閾値境界、欠落BPP、空評価域、交差の二実装の一致は追加の6 testsに含めた。
科学的な誤差の集計と、これらの開発上の検査を分けて扱う。
'''
    (LOCAL/'README.md').write_text(text)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('kind', choices=('local','recovery'))
    args=parser.parse_args()
    {'local':local_errors, 'recovery':recovery}[args.kind]()
    out=LOCAL if args.kind=='local' else RECOVERY
    write_json(out/'report_provenance.json', dict(generated_utc=now(),
        report_script_sha256=digest(Path(__file__)),
        test_script_sha256=digest(ROOT/'tests/test_local_errors_and_pairset_recovery.py'),
        source_provenance_sha256=digest(out/'provenance.json'),
        readme_sha256=digest(out/'README.md'),
        validation_sha256=digest(out/'validation.json') if (out/'validation.json').exists() else None))
