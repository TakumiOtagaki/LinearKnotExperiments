# Figure 3 と段階別 RSS（febd2bc1）


JA：**Figure 3. 長いRNAと合成対照における時間とpeak RSS。**
版febd2bc1、b=h=d=100、outer順位、Turner 2004、37°Cを用いた。
cpu_40割り当て内の同時実行数上限は10プロセス、各プロセス1 thread、RLIMIT_ASは8 GiBである。
LinearCPartyはMFEと分配関数・BPP、LinearPKCapRはloop profile、LinearPKRaccessはprofileと窓長1/4/8/16のaccessibilityを要求した。
Wall timeは起動・予備計算・本計算・復元・出力を含み、scaffold構築とbuildを含まない。
Peak RSSはプロセス全体の値で、GNU timeのKiBを2^20で割ったGiBを使う。
出力の異なるCLI全体の比較から、個別解析の追加コストは分離しない。

EN: **Figure 3. Time and peak RSS for long RNAs and synthetic controls.**
Calculations use febd2bc1, b=h=d=100, outer ranking, Turner 2004 and 37°C.
Each cpu_40 allocation permits at most 10 concurrent processes, each with one thread and an 8 GiB RLIMIT_AS limit.
LinearCParty requests MFE and the partition function/BPP; LinearPKCapR requests loop profiles; LinearPKRaccess requests profiles and accessibility at window lengths 1/4/8/16.
Wall time includes startup, preliminary and main computation, traceback and output, but excludes scaffold construction and builds.
Peak RSS is for the entire process, in GiB obtained by dividing GNU time KiB by 2^20.
The CLIs request different outputs, so their total costs do not isolate the cost of an individual analysis.

JA：**A, B**：RNAcentralの7配列（4,793–29,106 nt）の長さとwall time（A）、peak RSS（B）を示す。
点はwarm-up1回を除く測定3回の中央値、縦線は最小–最大であり、信頼区間ではない。
3回の中央値は昇順の中央の測定値である。
ほぼ重なるCapR/Raccessの点も省略しない。
各長さは異なる配列Sとscaffold Gの組であるため、点を結ぶ回帰線や漸近次数は示さない。

EN: **A, B**, Sequence length against wall time (A) and peak RSS (B) for seven RNAcentral inputs of 4,793–29,106 nt.
Points show the median of three measured repetitions after one excluded warm-up; vertical lines show the minimum–maximum, not confidence intervals.
The median of three repetitions is the middle ordered measurement.
Nearly coincident CapR/Raccess points are retained.
Each length represents a different sequence S and scaffold G; no regression line or asymptotic order is shown.

JA：**C**：合成対照の主要7系列について、250/500/1,000/2,000/4,000 ntにおけるLinearCPartyの時間とRSSを両対数軸で示す。
各系列・長さの記号は、各seedの測定3回の中央値を求めた後、その3 seedの中央値を取った値である。
縦線は3 seed×3回の計9測定の最小–最大を示し、反復とseedの両方の幅を含む。
線は同じ系列の表示点を結ぶだけで、補間値や理論上界ではない。
全156入力のうち、5長さを共通に持つ7系列×5長さ×3 seed=105入力を図示した。
残る51入力は追加の分岐深さ対照であり、この図の対象には含めない。
非emptyのGはfloor(0.15n)対、GC比率は50%である。
形を変える際に配列Sも変わるため、形や特徴量の独立な因果効果とは解釈しない。
3 CLIを同じ軸範囲で示す詳細図は `figure3_synthetic_all_cli_*` にある。

EN: **C**, LinearCParty wall time and RSS for seven main synthetic series at 250/500/1,000/2,000/4,000 nt, on logarithmic axes.
For each series/length, symbols show the median across three seeds of the within-seed medians of three repetitions.
Vertical lines span the minimum–maximum of all nine measurements (three seeds × three repetitions), combining repeat and seed variation.
Lines only connect plotted values in the same series; they are not interpolated measurements or theoretical bounds.
The plots show 105 of the 156 synthetic inputs: seven series × five shared lengths × three seeds.
The remaining 51 inputs are additional branch-depth controls outside this plotted subset.
Nonempty G has floor(0.15n) pairs and GC content is 50%.
Sequence S also changes with scaffold shape, so these comparisons do not identify independent causal effects of shape or features.
The `figure3_synthetic_all_cli_*` detail compares all three CLIs on shared axis ranges.

JA：**Supplementary 別図**：同じ7 RNAに対して、inside（outer）、MFE、outside+BPP、loop profileの独立したプロセスのpeak RSSを示す。
各モードはwarm-upなしの1回測定であり、反復の範囲はない。
後三者は必要なforwardと準備も含む別計算で、insideへの追加分ではない。
Peakを積み上げたり差し引いたりしてbyte内訳にはしない。
予備計算の一部を省くinside_baselineと旧版の測定は図に含めない。
段階別 RSS は `figure3_mode_memory_*` に保存する。

EN: **Separate Supplementary figure**, Peak RSS of independent processes for inside (outer), MFE, outside+BPP and loop profiles on the same seven RNAs.
Each mode has one measurement without warm-up and therefore has no repeat range.
The latter three modes include their required forward computation and preparation; they are separate computations rather than additions to an existing inside calculation.
Peaks are neither stacked nor subtracted to infer byte components.
The inside_baseline mode that omits part of the preliminary work, and measurements of the earlier version, are not plotted.
Stage RSS is saved separately as `figure3_mode_memory_*` for the Supplementary material.

JA：**注記あり案**：11,832 nt（URS00023C9544）と13,652 nt（URS000241043C）のLinearCPartyの点を矢印で示す。
注記の時間比・RSS比の範囲は3 CLIの入力別中央値の比の最小–最大であり、反復範囲や信頼区間ではない。
Gで一つのloopに直接属する自由塩基数の最大と、同じloop内の自由端点対数J_Vを既存特徴量表から添えた。
J_Vの比は13,652 ntの値を11,832 ntの値で割ったものである。
異なるS/Gの観測比較であり、時間差の原因をこの図だけで確定しない。

EN: **Annotated option**, Arrows mark the LinearCParty points for 11,832 nt (URS00023C9544) and 13,652 nt (URS000241043C).
The annotated time/RSS ratio ranges span ratios of input medians across the three CLIs, not repetitions or confidence intervals.
The largest number of free bases directly belonging to one G loop and the number J_V of free endpoint pairs in the same loop come from the existing feature table.
The J_V ratio divides the 13,652 nt value by the 11,832 nt value.
These are observations for different S/G inputs; the figure alone does not establish the cause of the time difference.

## 出典 / Sources

| Source ID | Aggregate file |
|---|---|
| `e2_metrics.tsv` | [metrics.tsv](../../local/e2-febd2bc1-20261016/metrics.tsv) |
| `e2_comparisons.tsv` | [comparisons.tsv](../../local/e2-febd2bc1-20261016/comparisons.tsv) |
| `e2_memory_metrics.tsv` | [memory_metrics.tsv](../../local/e2-febd2bc1-20261016/memory_metrics.tsv) |
| `e2_summary.json` | [summary.json](../../local/e2-febd2bc1-20261016/summary.json) |
| `g_features` | [rnacentral_g_features.tsv](../../local/e2-febd2bc1-g-features-20261016/rnacentral_g_features.tsv) |

| Panel | Rows / columns and processing |
|---|---|
| 3A/3B | E2 `metrics.tsv`, `run=rnacentral_current`, all 3 programs: `n`, `wall_seconds_{median,min,max}`, `rss_gib_{median,min,max}`; checked against the saved repetition arrays. |
| 3C/detail | E2 `metrics.tsv`, `run=synthetic_current`, seven named `series`, all five lengths; main C uses `program=LinearCParty`, detail uses all 3; median of three seed medians, extrema over all nine repeats. |
| SI mode memory | E2 `memory_metrics.tsv`, `run=memory_current`, `mode=inside_outer/mfe_outer/outside_bpp_outer/loop_profile_outer`, `n,rss_gib`; independent one-shot measurements. |
| 3 annotation | E2 `comparisons.tsv`, `kind=rna_13652_over_11832`, three programs, `wall_seconds_ratio,rss_gib_ratio`; G features `loop_free_max,same_loop_pairs_turn` for URS00023C9544/URS000241043C. |

JA：元データのSHA-256は `manifest.json`、実際に描いた集約値・閾値件数は `plot_data.tsv` にある。
RNA-levelの分位点は元のquantiles.tsvと照合し、E2の中央値と範囲はmetrics.tsv内の反復配列と照合した。
検証件数と出力の物理サイズは `validation.json` に記録する。

EN: `manifest.json` records source SHA-256 hashes; `plot_data.tsv` records plotted summaries and threshold counts.
RNA-level quantiles were checked against the original quantiles.tsv, and E2 medians/ranges against repetition arrays in metrics.tsv.
`validation.json` records check counts and physical output sizes.


## 再生成 / Regeneration

```bash
PYTHONDONTWRITEBYTECODE=1 python3 experiments/scripts/make_figures.py --figure3-only --figure3-panels abc --e2 experiments/results/local/e2-febd2bc1-20261016 --e2-revision febd2bc1 --out /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments/results/figures/draft_febd2bc1_20261016
```
