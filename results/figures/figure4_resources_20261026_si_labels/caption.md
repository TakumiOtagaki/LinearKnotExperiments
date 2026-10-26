# Figure S3

### 図 S3：合成入力の 2,000 から 4,000 塩基への分配関数・塩基対確率・MFE の資源使用量の比

図 S2 の七つの scaffold パターンを使う。
塗った丸は時間、白抜きのひし形は最大メモリを示す。
各長さで、各配列の三測定中央値を取り、次に三配列間の中央値を取り、長さ間で割る。
三つの seed は配列生成に使い、構造計算は決定的である。
時間比は 1.36–2.43、最大メモリ比は 2.05–2.67 であり、これらの測定は漸近的な計算量の次数を決めない。


English (LaTeX):

Synthetic partition-function, pairing-probability, and MFE cost ratios from 2,000 to 4,000 nt for the seven scaffold patterns of Figure~\ref{fig:si_synthetic_all}.
Filled circles show time; open diamonds show peak memory.
At each length, we take each sequence's three-measurement median, then the median across three generated sequences, before dividing between lengths.
The three seeds control sequence generation; folding is deterministic.
Time ratios range from 1.36 to 2.43 and peak memory ratios from 2.05 to 2.67.
These measurements do not determine an asymptotic order.
