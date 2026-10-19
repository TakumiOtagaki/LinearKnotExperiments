# Reachability of the 1,040 crossing-G reference pairs

2026-10-19。解析 D の参照集合を、canonical な注釈対のうち G の少なくとも一対と交差する対に統一した。
`pk-recovery-pairset-89baea52-20261019/pair_metrics.tsv` の200本・Turner04・`new_reference=1`の1,040対を使う。

- 既存の89baea52到達可能性検査983対を座標で結合し、その対確率が保存済みの新集合の値と完全一致することを確認した。
- 追加の57対（23本）だけを、同版の固定済み `build/membership-89baea52/membership_driver` の制約付きexact評価で検査した。
  各ターゲットは追加対一つを指定し、Gは全体を保持する。57対すべてが到達可能で、生成された追加対集合も指定座標に一致した。
- 合計1,040/1,040対が個別に到達可能だった。参照対全体の同時成立を示す検査ではない。
  完全canonical参照の164/200本という結果とエネルギー差は、同じ200本の既存検査を維持した。

| j−i (nt) | Reference pairs | Individually reachable | Probability >0.5 | Mean probability |
|---|---:|---:|---:|---:|
| ≤50 | 706 | 706 | 150 | 0.20938669434073864 |
| 51–100 | 183 | 183 | 25 | 0.1344479323004064 |
| >100 | 151 | 151 | 4 | 0.016232824840285515 |

`pairs.tsv` は1,040対と出典、`summary.json` は距離別集計、`raw.jsonl` は23プロセスのコマンド・標準入力・評価結果・終了状態を保持する。
`provenance.json` は保存入力、固定binary、parameter、解析scriptと出力のhashを記録する。
Binaryとparameterは元のarray planのsoftware lockと照合し、全入力の計算前後のhashも一致した。
各RNAの出力を個別ファイルにせず、共有領域には一つのJSONLにまとめた。
外部実装の変更、build、scheduler投入、git commitは行っていない。

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B scripts/extend_pk_membership_pairset_89baea52.py \
  --out results/local/pk-membership-pairset-89baea52-20261019
```

再実行時は新しい出力先を指定する。既存directoryの上書きは拒否する。
