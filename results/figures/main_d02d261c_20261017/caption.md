# Caption drafts / caption 案

Figure numbers follow the task labels; the manuscript editor can renumber the resource figure after Terc.
各英語段落の直前に、対応する日本語段落を置く。

## Figure 3 — figure3_main.pdf

**Time and memory for long RNAs and synthetic controls.**
A/BはRNAcentralの七入力（4,793–29,106 nt）について異なる出力を要求したプログラム全体の費用を示す。
MFEはminimum free energy、BPPはbase-pairing probabilityを指す。
LinearCPartyはMFE・分配関数・BPP、LinearPKCapRはcontext確率、LinearPKRaccessはcontext確率と1/4/8/16 nt窓のtotal accessibilityを計算する。
個別の解析の追加コストは分離しない。

A and B report complete programs with different requested outputs on seven RNAcentral inputs of 4,793–29,106 nt.
MFE denotes minimum free energy and BPP denotes base-pairing probability.
LinearCParty computes MFE, the partition function and BPPs; LinearPKCapR computes context probabilities; LinearPKRaccess computes context probabilities and total accessibility for 1-, 4-, 8- and 16-nt windows.
These totals do not isolate the added cost of a single analysis.

点はwarm-up一回を除いた三測定の中央値、縦線は全範囲である。
一部の範囲は記号より小さい。
長さごとに異なる配列とGを使う。
時間は起動・予備計算・主計算・復元・出力を含み、scaffold構築とbuildを除く。
メモリはプロセス全体のpeak resident set size（RSS）である。

Symbols show medians of three measurements after one excluded warm-up; vertical bars show the full range.
Some ranges are smaller than the symbols.
Each length represents a different RNA and scaffold G.
Wall time includes startup, preliminary ranking, main computation, traceback and output, and excludes scaffold construction and builds.
Memory is the peak resident set size (RSS) of the entire process.

CはLinearCPartyの合成七系列における2,000→4,000 ntの時間とRSSの比を示す。
比は、各seedの三反復の中央値を三seed間でまとめた中央値同士から計算する。
点線は比2を示す。
比は測定区間の観測であり、漸近計算量の上界ではない。
系列ごとに配列とGの両方が変わるため、Gの形だけの因果効果とは解釈しない。
全長さと三プログラムの曲線は補足図に置く。

C shows the ratios of LinearCParty time and RSS between 4,000 and 2,000 nt for seven synthetic series.
At each length, we first take the median of three repetitions for each seed, then the median across three seeds.
Ratios divide the resulting 4,000-nt median by the 2,000-nt median; the dotted line marks a ratio of two.
These ratios describe the measured interval and are not asymptotic bounds.
Both sequence and scaffold change between series, so the comparisons do not identify a causal effect of scaffold shape alone.
The supplementary synthetic figure shows all five lengths and all three programs.

## Figure S2 — figure3_mode_memory_si.pdf

**Peak memory of independent calculation modes.**
各点は別のRNAとGの組を示す。
insideは順位付けの予備計算を含み、他のモードは必要なforward計算と準備をそれぞれ含む。
各モードは独立したプロセスで、warm-upなしの一回測定である。
異なるプロセスのpeakを足したり引いたりしてメモリの内訳とはしない。

Each length represents a different RNA and scaffold G.
Inside includes preliminary ranking; the other modes each include their own required forward calculation and preparation.
Each mode runs in an independent process and has one measurement without warm-up.
The peaks cannot be added or subtracted to infer memory components.

## Additional SI — figure3_synthetic_all_cli_si.pdf

**Measured costs for seven synthetic sequence–scaffold series.**
三プログラムの時間とRSSを対数軸で示す。
点は各seed三反復の中央値を三seed間でまとめた中央値、縦線は九測定の全範囲である。
線は同じ系列の測定点を結ぶ。
全156入力のうち、共通の五長さを持つ七系列×五長さ×三seedの105入力を示す。
残る51入力は追加の分岐深さ対照である。
非空のGはfloor(0.15n)対を持ち、配列のGC含量は50%である。
下の弧図は保存済みの250 ntのGを同じ縮尺で示す。

Time and RSS use logarithmic axes for all three programs.
Symbols show medians across three seeds of within-seed medians of three repetitions; vertical bars span all nine measurements.
Lines connect measured values within each series.
The figure shows 105 of the 156 inputs: seven series, five shared lengths and three seeds.
The remaining 51 inputs are additional branch-depth controls.
Each nonempty G contains floor(0.15n) pairs, where n is sequence length; sequences have 50% GC content.
The arc sketches show the saved 250-nt scaffolds on a common coordinate scale.

## Additional SI — figure3_parallel_impact_si.pdf

**Concurrent and single-process measurements.**
A/Bは三つのRNAについて、各入力の三測定の中央値を単独実行w1の中央値で割った時間/RSS比を示す。
塗った記号はw10/w1、白抜き記号は本文の測定/w1を示す。
w10は最大十プロセスを許可し、観測された同時実行の最大数は九だった。
各設定は同じ資源型の異なるホストで実行し、点線は比1を示す。
これらの比で他の測定を補正しない。

A and B show time and RSS ratios for three RNAs, dividing each input's median of three measurements by its one-process (w1) median.
Filled symbols show w10/w1; open symbols show Main/w1.
The w10 setting allowed ten processes and reached nine concurrent processes.
Settings used different hosts of the same resource type; dotted lines mark a ratio of one.
These ratios do not correct the other measurements.
