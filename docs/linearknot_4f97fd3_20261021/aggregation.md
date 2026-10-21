# 4f97fd3 の時間・メモリ集計（2026-10-21）

12 matrix・82 task・2,768 process はすべて成功し、測定の欠測は0。
精度・検証・Terc・回復の採用値は89baea52のまま保持した。

- [完了数と旧→新の主要値](../../results/local/timing-4f97fd3-20261021/README.md)
- [全575条件の旧新版比較](../../results/local/timing-4f97fd3-20261021/version_comparison.tsv)：E2、NEAT1、独立モード、w1/w10。
- [PFだけの168条件の旧新版比較](../../results/local/timing-4f97fd3-20261021/pf_version_comparison.tsv)
- [NEAT1全33条件](../../results/local/timing-4f97fd3-20261021/neat1_resources.tsv)
- [PFの長さ別集計](../../results/local/timing-4f97fd3-20261021/pf_by_length.tsv)：全成功入力と共通成功入力を区別する。
- [E2集計](../../results/local/e2-4f97fd3-20261021/summary.json)、[PF集計](../../results/local/cparty-linearknot-pf-4f97fd3-20261021/summary.json)、[検証](../../results/local/timing-4f97fd3-20261021/validation.json)
- [本文とSIの4図](../../results/figures/main_4f97fd3_20261021/README.md)：PNGを開いて軸・凡例・範囲・配置を確認した。最小文字7 pt、枠外・重なり・数値検査が成功した。

E2はwarm-up1＋測定3の中央値、独立モードはwarm-up0＋測定1、NEAT1もwarm-up0＋測定1。
PFはwarm-up1＋測定3の入力ごとの中央値を、共通成功入力間でさらに中央値にした。
CPartyは保存済みの6 patch oracleのPFのみを再利用し、LinearKnotの旧時間は再利用していない。
時間とRSSは新JSONLから集計した。
入力・パラメータ・版・反復数・checkpointとaggregate・raw出力のhashを確認した。
E2解析は旧新版・対照を合わせて21,248個のraw artifactを検証した。
全tarはstreamで読み、共有領域にもTMPDIRにも展開していない。
旧新版は日付・ホストも異なるため、観測値の差を実装だけの高速化率とは呼ばない。

## CPartyとの比較

250/500/1,000/2,000/4,000 ntの各21入力のうち、CPartyの非有限値は0/0/12/15/21入力。
失敗した192 processはすべて `nonfinite oracle inside`。
成功した共通入力では1,000 ntまでCPartyが速く、2,000 ntではb=100が23.40 s、CPartyが34.97 s（6入力）。
正確な交差長を補間して原稿に載せていない。
4,000 ntには共通成功入力がないので、SI表S21ではb=100全21入力の47.58 s・931.1 MiBを†で区別する。
CPartyの失敗値を成功時の時間・RSSに入れない。
成功したCParty/exactの51入力はすべてPFの許容差を通った。
Exactは検証用とし、1,000 ntより長い設定は計画どおり除外した。

## 再現コマンド

experiments/から。出力先が存在する解析・作図は上書きしない。
再現時は以下の出力先を別名に替える。

```bash
PYTHONPATH=scripts python3 -B - <<'PY'
from pathlib import Path
from aggregate_array import aggregate
for p in sorted(Path('results/runs').glob('*-4f97fd3-20261021/array_plan.json')):
    aggregate(p.resolve())
PY
```

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
MPLCONFIGDIR="$TMPDIR/lk4-mpl" python3 -B scripts/analyze_e0_e2_runs.py \
  --synthetic-run results/runs/e2-synthetic-shapes-4f97fd3-20261021 \
  --rna-run results/runs/e2-rnacentral-4f97fd3-20261021 \
  --previous-rna-run results/runs/e2-rnacentral-89baea52-20261018 \
  --previous-synthetic-run results/runs/e2-synthetic-shapes-89baea52-20261018 \
  --workers 10 --revision 4f97fd3 --previous-revision 89baea52 \
  --previous-workers 10 --previous-resource cpu_40=1 --skip-cap30-history \
  --memory-run results/runs/e2-memory-breakdown-4f97fd3-20261021 \
  --previous-memory-run results/runs/e2-memory-breakdown-89baea52-20261018 \
  --parallel-w1-run results/runs/e2-parallel-impact-w1-4f97fd3-20261021 \
  --parallel-w10-run results/runs/e2-parallel-impact-w10-4f97fd3-20261021 \
  --out results/local/e2-4f97fd3-20261021
python3 -B scripts/analyze_pf_scaling.py \
  --plan results/runs/linearknot-pf-4f97fd3-20261021/array_plan.json \
  --cparty-plan results/runs/cparty-linearknot-pf-89baea52-20261020/array_plan.json \
  --out results/local/cparty-linearknot-pf-4f97fd3-20261021
python3 -B scripts/summarize_timing_4f97fd3.py --out results/local/timing-4f97fd3-20261021
MPLCONFIGDIR="$TMPDIR/lk4-mpl" python3 -B scripts/plot_timing_4f97fd3.py \
  --out results/figures/main_4f97fd3_20261021
```

## 開発上の確認

共有tarのhash照合が小さいread単位で遅かったため、既存の逐次tar読み出しbufferを1 MiBにした。
読取内容とhash検査は変えず、入力・rawファイルは変更していない。
保存処理8件、E2解析10件、PF解析5件の既存testが成功した。
保存処理testの `Disk quota exceeded` は意図した失敗経路の出力で、実際のquota到達ではない。
LaTeXの新しい表に対するlatexdiffの罫線と、削除済みNEAT1図の旧参照は、論文repoの生成scriptだけで補正する。
新しいRNA計算、qsub、git commit、外部実装とbibの変更は行っていない。
