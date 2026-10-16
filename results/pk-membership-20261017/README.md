# 参照 pseudoknot の所属とエネルギー（Major 2、2026-10-17）

この README と隣接する小さな集計は `scripts/report_pk_membership.py` で再生成する。
原データ、RNA ごとの判定、対ごとの確率、エネルギー内訳は
`results/local/pk-membership-20261017/` に置く（git 管理外）。
原稿と実装は変更しない。

## 結果

| 評価対象 | 成功 / RNA 数 | 拒否 / RNA 数 |
| --- | --- | --- |
| full：canonical な注釈全体 | 164/200 | 36/200 |
| filtered：full から span < 4 の対を除外 | 176/200 | 24/200 |
| reported：既存 G ＋既存の [] 対 | 198/200 | 2/200 |
| scaffold：既存 G のみ | 199/200 | 1/200 |

full の拒否理由（重複可）は {"completed_WMB_weak_closure": 1, "extension_crossing": 22, "minimum_pair_span": 14}。
短い対を除いた `filtered` では、これらに加えて `bpRNA_PDB_598` の組合せが
`constrained_DP_rejection_unresolved` となる。短い対を含む full では同じ RNA に span 違反がある。
`bpRNA_RFAM_33327` の reported は LegacyCompatible では評価できるため、既定の completed-WMB weak closure による拒否と切り分けた。
`bpRNA_tmRNA_144` の G のみは V(51,232) → V(52,71) の unpaired な両側の合計が 160 nt となり、MAXLOOP=30 を超える。
この RNA は参照 pseudoknot を加えると評価できる。
拒否した構造のエネルギーは定義せず、TSV では NA とする。
`reported` は以前の回復率の集計と同じ G と 983 対を使う。
単独追加は各対 p について完成構造を G ∪ {p} に固定する検査であり、
同じ対を含む別の完成構造の有無とは区別する。
正の exact 周辺確率で、何らかの完成構造に現れることを確認できたのは
983/983 対だった。
最新固定版での独立した所属の証拠は、正の BPP 973 対と、
タイムアウト RNA の完了済み参照全体の traceback に含まれる
10 対である。
したがって最新固定版でも全 983 対が少なくとも一つの許容構造に含まれる。
個別の参照対を一律に禁止していることでは低回復率を説明できない。

| 対の分類 | 分母（対） | p > 0.5（対） | 割合 |
| --- | --- | --- | --- |
| 全参照 [] 対 | 983 | 179 | 18.21% |
| 何らかの完成構造で許容を確認 | 983 | 179 | 18.21% |
| 完成構造への所属の証拠なし | 0 | 0 | NA |
| 単独追加で評価成功 | 980 | 179 | 18.27% |
| 単独追加で拒否 | 0 | 0 | NA |
| 単独追加の検査未完了 | 3 | 0 | 0.00% |
| 注釈全体の評価成功 RNA の対 | 828 | 150 | 18.12% |
| 注釈全体が拒否された RNA の対 | 155 | 29 | 18.71% |

空の分母の割合は NA とする。
単独追加の検査未完了の 3 対は `bpRNA_RFAM_12431` の
(14,257), (13,258), (12,259)。
これらを拒否された対には数えない。
「注釈全体が拒否された RNA の対」は、その対自体が禁止されていることを意味しない。
確率は固定した G の条件付き ensemble に対する値である。

### エネルギー差

単位は kcal/mol。MFE も参照評価も `d02d261c` の exact 計算である。
分位点は並べた値の線形補間で求め、評価成功した RNA を分母とする。

| 対象 − exact MFE | 分母（RNA） | 中央値 | Q1–Q3 | 最小–最大 |
| --- | --- | --- | --- | --- |
| full | 164 | 11.039 | 4.800–20.212 | 0.000–82.976 |
| filtered | 176 | 11.513 | 4.827–22.515 | 0.000–82.976 |
| reported | 198 | 13.010 | 5.520–21.490 | 0.000–82.976 |
| scaffold | 199 | 11.336 | 6.957–19.042 | 2.100–66.888 |

full が評価できた 164 本のうち、差が 1e-8 kcal/mol を超えたのは
159 本だった。
MFE の生成対と参照 [] 対の重なりは 210/983 対。
MFE の Square lane は合計 850 対で、参照と一致したのは
210 対だった。
RNA ごとの MFE energy、重なり、参照 energy と G のみの energy の差も `rna.tsv` に保存する。

### stem の距離との関係

ここで距離は参照対の両端の隔たり j − i とする。
stem の長さは、保存済みの canonical [] 対で (i+1,j−1) と連続する最大 stack の対数とし、
欠けた対をつながない。対ごとの span、stack 長、stack 内の span の中央値を `pairs.tsv` に保存する。

| 対の span = j − i | 対数 | p > 0.5 | 平均 p | RNA 数¹ | RNA 平均 p¹ |
| --- | --- | --- | --- | --- | --- |
| ≤50 nt | 689 | 150/689 | 0.214554 | 144 | 0.195065 |
| 51–100 nt | 143 | 25/143 | 0.172022 | 30 | 0.165428 |
| >100 nt | 151 | 4/151 | 0.016236 | 26 | 0.023513 |

¹ RNA の分類は各 RNA の [] 対の span の中央値で行う。
RNA ごとの平均 p をさらに RNA 間で平均するため、対ごとの平均 p とは重みが異なる。
span と回復の Spearman 相関は RNA 単位で
-0.328、対単位で
-0.328。
対は同じ RNA 内で依存するため、対を独立な標本とする検定は行わない。
距離による違いは記述的な関連であり、エネルギー項の因果効果とはしない。

### 参照と MFE のエネルギー内訳

full が評価できた同じ 164 本を参照・MFE 双方の分母とする。
`exterior`、`in_multiloop`、`in_pseudoloop` は pseudoknot を置く場所に応じた開始項、
`band` は band の項、`unpaired` は pseudoloop 内の unpaired 項である。
`band_stack` と `band_internal` は band 内の stack / internal-loop のエネルギー。
`other_thermodynamic` は通常の loop・外側など残りの熱力学項を表す。

| 項（kcal/mol） | 参照平均 | MFE 平均 | 差の平均 |
| --- | --- | --- | --- |
| exterior | -0.934 | -1.279 | 0.345 |
| in_multiloop | 4.237 | 0.000 | 4.237 |
| in_pseudoloop | 0.366 | 0.000 | 0.366 |
| band | 6.195 | 4.815 | 1.380 |
| unpaired | 1.067 | 0.454 | 0.614 |
| nested_region | 0.779 | 0.293 | 0.486 |
| band_stack | -25.213 | -14.054 | -11.158 |
| band_internal | 2.554 | 1.369 | 1.186 |
| band_multiloop_closing | 0.956 | 0.229 | 0.728 |
| band_multiloop_branch | 0.471 | 0.123 | 0.348 |
| band_multiloop_unpaired | 0.504 | 0.078 | 0.427 |
| multiloop_branch | 0.013 | 0.000 | 0.013 |
| other_thermodynamic | -3.925 | -21.430 | 17.505 |

開始項は配置により符号も異なるため、全てを一律の正の penalty と解釈しない。
この比較は構造間の各項の差であり、パラメータを変えて再最適化する実験ではない。
全項と出現数、reported の対応比較は `energy_terms.tsv` と `energy_terms_summary.tsv` に残す。
reported と G のみの energy の差の中央値は 0.147 kcal/mol
（G のみも評価できた 197 本）。
MFE との差には、参照にない追加対を作る効果も含まれる。
MFE との差を pseudoknot の単一の energy penalty と同一視しない。

## 入力、手順と再現

1. 既存の `pk-recovery-bprna-20261016/pk_sel.json` と `pk_results.json` を読む。
   200 本と 983 個の canonical [] 対を維持する。
   `nested` が canonical かつ span ≥ 4 の () 対、`pk` が canonical [] 対と一致することを検査する。
2. G は既存 `nested` を固定する。full は生の注釈の全 bracket 種から canonical 対だけを残し、
   G 以外の対を指定する。元の注釈から除いた非 canonical 対は 725 対。
   既存 G が除いていた短い canonical 対は 18 対 / 14 本、
   既存の [] 集計に含めなかった {} / <> の canonical 対は 57 対 / 23 本。
3. read-only の固定版ソースと静的 library に外部 observer をリンクする。
   公開 `ExactSession::evaluate_energy` の実装と同じ output-constrained exact Viterbi を呼び、
   traceback の生成対が指定した対集合と一致することを確認する。
   MFE と outside も同じ固定版で求める。全て既定の `RequireWeakClosure`、dangles=2、Turner04。
4. 参照評価が失敗した場合、span、G′ 内の交差、endpoint の共有、canonical 条件を独立に確認する。
   G の対に囲まれ、G′ の endpoint を含まない単一-child の V internal loop には MAXLOOP=30 を確認する。
   さらに `LegacyCompatible` で再評価し、completed-WMB weak closure の変更だけで回復するか調べる。
   どの根拠でも説明できなければ `constrained_DP_rejection_unresolved` とする。
5. 選ばれた同じ traceback の local weight を足して root energy と照合する。
   pseudoknot の各線形パラメータを +1 kcal/mol とした snapshot で同じ trace だけを再採点し、
   差から出現数を求める。stack / internal の倍率を 0 とした再採点から band の熱力学項を分ける。
   各項の和、energy = −RT log(weight)、指定対の完全一致、参照 energy ≥ MFE を検査する。

一時 build は node-local `$TMPDIR` に置き、終了時に消す。
生ログは RNA ごとのファイルを作らず一つの `raw.jsonl` に保存する。
cohort、元の確率、observer source、build log は tar.gz にまとめる。
実行の最初に共有 inode quota を確認した。結果の runtime は並列解析の実行記録であり、性能評価ではない。

```bash
python scripts/analyze_pk_membership.py --workers 4
# 60 秒で打ち切った RNA は outside の再計算を省いて構造評価だけ回収する。
python scripts/analyze_pk_membership.py --recover-energy-only --workers 1
python scripts/analyze_pk_membership.py --aggregate-only
python scripts/report_pk_membership.py
```

この cohort では最初の計算と再試行は時間上限のため非ゼロ終了する。
その後の `--aggregate-only` が、保存済みの完全な値だけを集計し、未完了の検査を明示する。
失敗と再試行はどちらも JSONL に残し、RNA ごとに最後の試行を集計する。
`bpRNA_RFAM_12431` は初回と outside を省いた再試行の双方で 60 秒上限に達したため、追加計算を止めた。
再試行で出力済みの完全な JSON 値だけを読み、MFE・参照全体・G のみ・単独対 7 個を回収した。
未完了の JSON 値を修復して測定値として使うことはしない。
この部分出力の RT は、完了済みの MFE energy と log(weight) から計算した。
今回は 199 本 / 973 対で BPP を再計算し、
既存の `febd2bc1` exact 結果との差の最大は 0。
残り 1 本の BPP は以前の exact 出力を使い、`probability_source` 列で明示した。
全 RNA の所属・エネルギー・MFE は `d02d261c` で評価する。
MFE の再評価は 200/200 本で一致した。
local weight の和と root の最大 energy 差は 1.35e-13 kcal/mol、
項の和との差は 1.42e-14 kcal/mol。
入力・parameter・library・observer の SHA256、build command、失敗履歴は provenance と archive に保存する。
残る単独対 3 個と BPP 再確認用の未実行 array 案は `bpp_retry_array_plan.json` に置く。
cpu_4 × 1 task、1 worker、h_rt=30 分、RNA の上限600秒とする設定案であり、投入していない。

## 判定の根拠となる固定版ソース

全て `deps/LinearCParty-d02d261c`、commit `d02d261c3a07e7308a3a3d37543d040adfd0f86e`。

- `src/scfg2/api/exact.hh`：公開 `ExactSession::evaluate_energy`、`outside`、`viterbi`。
- `src/scfg2/api/exact.cc:180`：評価 API と内部 output-constrained Viterbi の接続。
- `src/scfg2/replay/w_final_exact_inside.hh:754`：target による deduction 制限と最終 pair-set の完全一致。
- `src/scfg2/engine/traceback_pairs.hh`：生成対の抽出と target の局所的な endpoint guard。
- `src/scfg2/deductions/generate_exact_basic.hh:120`：V / VP の endpoint、span、weak closure。
- 同ファイル `:208`：completed-WMB の weak closure、`:619` と `:1190`：V internal-loop の MAXLOOP 制限。
- `src/scfg2/weights/runtime_snapshot.cc:81` と `src/scfg2/weights/log_loops.hh`：canonical 判定と loop の重み。
- `src/ViennaRNA/params/constants.h:27`：TURN と MAXLOOP。
- `src/scfg2/replay/w_final_exact_inside.hh:145`、`:295`：開始・band・unpaired などの local weight の所有。
- `src/scfg2/weights/pseudoknot_parameters.hh` と `runtime_snapshot.cc:55`：pseudoknot モデル値と単位。
- `documents/scfg2/energy_evaluator_contract_20260716.md`：正の結果の soundness と、拒否結果の completeness の限界。
- `documents/eval_structure/scfg2_wmb_closure_implementation_20261010.md`：weak closure policy の意味。

## 限界

この所属は指定した G と実装文法に対する判定であり、別の G / G′ の割り当ては探索していない。
評価成功は指定した pair map の導出が存在する証拠になる。
失敗だけから一般に非所属とは断定できないため、拒否理由の独立した構造チェックも残す。
返す energy はその pair map の最良導出の energy であり、全長での非曖昧性の証明を仮定しない。
元の生物学的注釈そのものではなく canonical 化と既存 G の条件付き解析を行う。
低い周辺確率は競合構造とその重みの結果であり、特定のエネルギー項の誤りまでは特定しない。
exact にも文法固有の局所制限がある。beam pruning による近似誤差はこの解析には含まれない。

## SI 英文案（5 文、日本語と段落対応）

既存の G を固定すると、参照 pseudoknot の全 983 対は、
少なくとも一つの許容構造に含まれていた。
canonical な注釈全体は 164/200 本で評価できた。
全参照対の 179/983 対が確率 0.5 を超えた。
評価できた注釈全体の MFE との差の中央値は 11.04 kcal/mol
（四分位範囲 4.80–20.21）だった。
この結果はこの集合の低い回復を解釈する際に、個々の対の構造的な許容性と条件付き確率を区別する必要を示すが、
原因となるエネルギー項は特定しない。

With the supplied scaffolds held fixed, each of the 983 reference pseudoknot pairs occurred in at least one allowed structure.
The complete canonical reference annotations could be evaluated for 164 of the 200 RNAs.
Of all reference pseudoknot pairs, 179/983 had a conditional probability above 0.5.
For evaluable complete references, the median energy difference from the exact minimum-free-energy structure was 11.04 kcal/mol (interquartile range, 4.80–20.21 kcal/mol).
These results distinguish the structural admissibility of individual pairs from their conditional probabilities in this cohort, but do not identify a causal energy term.
