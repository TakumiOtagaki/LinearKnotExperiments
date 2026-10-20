**Figure 4. Time and memory for long RNAs and synthetic controls.**

AとBはRNAcentralの7入力（4,793–29,106 nt）で、三つの解析が要求する出力を含む実行全体の時間とメモリを示す。
LinearKnotのfoldは分配関数・塩基対確率（BPP）・最小自由エネルギー（MFE）構造、profileはloop profile、accessは1・4・8・16 nt窓のaccessibilityを求める。
計算にはbeam幅b=100、hairpin長の上限h=100、scaffoldのmultiloopを外側へ探す段数の上限d=100を使う。

**(A, B)** Whole-program time and memory are shown for seven RNAcentral inputs of 4,793–29,106 nt.
The legend gives the roles of the three LinearKnot subcommands: fold computes the partition function, base-pairing probabilities (BPPs) and a minimum-free-energy (MFE) structure; profile computes loop profiles; access computes accessibility for 1-, 4-, 8- and 16-nt windows.
Calculations use beam width b=100, hairpin-length limit h=100, and a limit d=100 on how many enclosing scaffold multiloop levels are searched outward.

記号はwarm-up一回を除いた3測定の中央値、縦の誤差線はその最小から最大の範囲である。
この範囲の多くは記号より小さく、見えにくい。
各長さは異なるRNAとGを表す。
メモリはプロセス全体のpeak resident set size（RSS）である。

Symbols show the median of three measurements after one excluded warm-up.
Vertical error bars extend from the smallest to the largest of those three measurements.
Many ranges are smaller than the symbols and are difficult to see.
Each length represents a different RNA and scaffold G.
Memory is the peak resident set size (RSS) of the whole process.

Cは7種類のscaffoldの形を持つ合成入力で、4,000 ntのfoldの時間とRSSを2,000 ntの値で割った比を示す。
各形・長さに対して、乱数seedを変えて3本の合成配列を生成した。
各配列の3測定の中央値を取り、次に3配列の中央値を取ってから、長さ間の比を求めた。
乱数は配列生成に用い、foldの計算は決定的である。
全長さの測定とscaffoldの図は補足図S2に示す。

**(C)** Ratios divide fold time and RSS at 4,000 nt by those at 2,000 nt for synthetic inputs with seven scaffold shapes.
For each shape and length, three random seeds generate three synthetic RNA sequences.
We take the median of three measurements for each sequence, then the median across the three sequences, before computing the ratio between lengths.
The random seeds control sequence generation; fold itself is deterministic.
Figure S2 shows all measured lengths and sketches of the scaffold shapes.

**Figure S2. Time and memory for seven synthetic scaffold shapes.**

A/Bはfold、C/Dはprofile、E/Fはaccessの時間とRSSである。
各解析の出力はFigure 4のcaptionに示す。
RNA長の対数軸は底2、時間とRSSは底10である。
各点の中央値の取り方はFigure 4Cと同じで、縦の誤差線は3配列×3測定の9測定の最小から最大を示す。
この範囲の一部は記号より小さい。
線は同じscaffoldの形の測定点を結ぶ。
全156入力のうち、7種類の形×5長さ×3配列の105入力を示す。
残る51入力は、scaffoldのmultiloopの入れ子の深さを変えた追加対照である。
空でないGはfloor(0.15n)対（nは配列長）を持ち、配列のGC含量は50%である。
下の弧図は保存済みの250 ntのGを同じ縮尺で示す。

**(A, B)** show fold, **(C, D)** show profile, and **(E, F)** show access; requested outputs are listed in the Figure 4 caption.
RNA length uses a base-2 logarithmic axis; time and RSS use base 10.
Symbols use the medians defined in Figure 4C; vertical error bars span the smallest to largest of the nine measurements (three generated sequences, three measurements each).
Some ranges are smaller than the symbols.
Lines connect measurements with the same scaffold shape.
The figure shows 105 of the 156 inputs: seven shapes, five lengths and three generated sequences.
The other 51 inputs are additional controls with different depths of nested scaffold multiloops.
Each nonempty G contains floor(0.15n) pairs, where n is sequence length; sequences have 50% GC content.
The arc sketches show the saved 250-nt scaffolds on a common coordinate scale.
