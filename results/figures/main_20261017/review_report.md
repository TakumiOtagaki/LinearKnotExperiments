# Figure 2・3・S1・S2 の改訂記録（2026-10-17）

担当: `scripts/make_figures.py` と新しい補助スクリプト `scripts/make_figures_review.py`、本出力 directory。
本文・SI の LaTeX、日本語原稿、bib、他の作図担当のファイル、調整中の実装、docs/progress.md は変更していない。
担当範囲の指示に従い、意味のある変更と検証結果は本記録と自分の last.txt に残す。
git commit、RNA CLI の実行、scheduler 投入はしていない。

## 出力とレビューへの対応

- Figure 2 (`figure2_main.pdf/png`): 178 × 92 mm。
  F2-1 は長さ群を横軸のカテゴリ、beam 幅を全行で同じ色と記号に変更した。
  F2-2 は C を長さ群別16 nt窓にし、0–0.000001の線形域を持つ対称対数軸で小さい差を表示した。
  0や小さい正値への数値の置換はしていない。
  300–399 nt、b=200では評価RNA数 1,048、p90 2.47610098e-05、最大 0.0304740993 である。
  F2-3 は四分位箱とhatchを外し、全範囲・10–90パーセンタイル・中央値に絞った。
  四分位の数値はplot_data.tsvに残した。
  F2-4 は塗った記号と白抜き記号をpairing/contextのキーと対応させ、captionでTVとRNA内平均→RNA間分位点を定義した。
  F2-5 は数量軸の傾きに見える折れ線を使わず、群内の設定は等間隔と明記した。
  F2-6 はb=400の自己比較を全て省き、理由をcaptionに書いた。
- 時間とメモリの図 (`figure3_main.pdf/png`): 178 × 55 mm。
  F3-1/F3-3 はCを七系列の2,000→4,000 ntの比へ変更し、系列名を各行に直接表示した。
  比は図示する中央値同士の比で、seedごとの比の中央値ではない。
  全長さ・全CLIの曲線と保存済みGの模式図は追加SI図に移した。
  F3-2 は同じ座標で重なるsquare/triangleの大きさと塗りを変え、凡例から線を外し、最長入力の注記を集計値から生成した。
  F3-4 は凡例に要求出力を加え、captionで異なる出力を要求するプログラム全体の費用とした。
- S1 (`figure2_d_si.pdf/png`): 165 × 55 mm。
  S1-1/S1-2 は閾値を軸名に入れ、A/B、分母、片側の条件、自己比較の省略、カテゴリ軸を明記した。
  色と記号はFigure 2と同じbeam幅を表す。
- S2 (`figure3_mode_memory_si.pdf/png`): 165 × 65 mm。
  S2-1/S2-2 は予備順位計算と独立プログラムの名前を明記し、異なるRNAを結ぶ線を外した。
  縦軸上限は保存された最大RSSから計算した。
- 追加SI図は `figure2_accessibility_si.pdf/png` (165 × 137 mm) と `figure3_synthetic_all_cli_si.pdf/png` (165 × 178 mm)。
  前者は窓長1/4/8/16の全長さ群、後者は3 CLI・七系列・五長さを示す。
  Gの模式図は保存済みinputs.jsonlから読み、そのhashがE2集計のdataset hashと一致することも検査した。
- G1–G3 は掲載幅での作図、7 ptの通常文字、通常の数字による対数目盛り、色に加えた記号・塗り・線種・直接ラベル、独立した9 pt太字パネル記号で対応した。
  Figure 2は行A–Cを一度だけ表示し、exact/有限参照を列見出しで区別した。

## 未反映と理由、著者・統合担当の判断

E2は全て同じfebd2bc1の保存集計を使った。
`--review-layout --e2-version febd2bc1` / `d02d261c` で切り替えられる。
異なる版の図は新しい出力先で保存する。
確認時刻は 2026-10-16T16:23:34.593184+00:00（日本時間2026-10-17未明）で、d02d261cは全runがまだ完了していない。
メモリrunだけが全task成功であり、未完のRNAcentralや合成系列と混ぜない。

| Run | Planned tasks | task.jsonの状態 |
|---|---:|---|
| e2-memory-breakdown-d02d261c-20261017 | 4 | success: 4 |
| e2-parallel-impact-w1-d02d261c-20261017 | 1 | running: 1 |
| e2-parallel-impact-w10-d02d261c-20261017 | 1 | running: 1 |
| e2-rnacentral-d02d261c-20261017 | 4 | running: 4 |
| e2-synthetic-shapes-d02d261c-20261017 | 55 | not_started: 49, running: 6 |

READMEに全task成功の確認、aggregate_array.pyによるroot manifestの更新、analyze_e0_e2_runs.pyによる同じschemaの集計、別directoryへの再作図のコマンドを記した。
root run.jsonは集計前には古い状態を持つため、完了判定にそのまま使わない。
未完了版のpresetは集計ファイルがなければ拒否し、出力先を作らないことを確認した。

Cの時間比の範囲は **1.51–2.63**、RSS比は **2.27–3.58** である。
本文の1.51–2.54はseedごとの比を先にまとめた別統計なので、新図とそろえる場合は1.51–2.63へ更新する必要がある。
全14比の分子・分母・値はplot_data.tsvの3C_ratio行とREADMEに保存した。
この作業では数値を本文へ差し込んでいない。

案Bで時間とメモリの図を本文Figure 4にする場合も、図自体に番号を書いていないためファイル名figure3_mainのまま差し込める。
追加SI二図の採否と最終番号、本文・SIの参照とcaptionの取り込みは著者・統合担当が決める。
caption.mdは日本語と英語を段落ごとに対応させ、一文一行とした。
本文の8ページ確認はLaTeX担当に委ねる。
主図の高さは旧版の92/55 mmを保った。

## 検証

- 出力6 PNGを実際に開き、重なり・はみ出し・ラベル・記号の対応を確認した。
  合成詳細図の薄色系列名は黒い文字へ変更し、色と記号は直下のキーに残した。
- Matplotlibの文字の重なり、ページ外へのはみ出し、軸範囲外の点の検査を6図すべてで通した。
- PDF実体のGhostscript txtwrite/TextFormat=0の結果は以下。
  この環境では-dNODISPLAYを付けるとtext出力が空になるため、-dNOPAUSE -dBATCHを用いた。

| PDF | 最小pt | 最大pt |
|---|---:|---:|
| figure2_main.pdf | 7 | 9 |
| figure2_d_si.pdf | 7 | 9 |
| figure2_accessibility_si.pdf | 7 | 9 |
| figure3_main.pdf | 7 | 9 |
| figure3_mode_memory_si.pdf | 7 | 7 |
| figure3_synthetic_all_cli_si.pdf | 7 | 9 |

- E1の保存分位点1,750 scalar、E2の反復要約2,934 scalar、G特徴量の対応147項目を照合した。
  E1は350要約、98,500短配列行と8,628長配列行、E2は489 CLI行と28独立モード行を確認した。
  長さ群・窓長ごとに評価域が全beamで同じことも検査した。
- 6 PDFのMediaBoxの物理サイズ、6 PNGの画素数と300 dpiメタデータを確認した。
- 旧15レイアウトの描画は全画素で一致した（historical_regression.json）。
  main以外の従来33関数・classのASTも一致した。
  旧出力30ファイルのSHA-256も開始時と一致した。
- 担当Pythonの構文と末尾空白、担当分のgit diff --checkは成功した。
- source/output/script hashと検証結果はmanifest.jsonとvalidation.jsonに記録した。
  一時ファイルは/tmpだけに置き、作業後に自分の一時directoryを削除する。
  共有領域に追加する図関係の成果物は21ファイルにまとめた。
