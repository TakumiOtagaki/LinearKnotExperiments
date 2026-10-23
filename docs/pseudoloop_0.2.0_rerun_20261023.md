# PseudoLoop 0.2.0 に伴う文脈出力の再計算（2026-10-23）

対象実装は LinearKnot **0.2.0、`f93b55e`**（tag なし）。実装 repo は読み取りと `git archive` のみに使用した。小規模計算、作図、大規模 E1 の dry-run・local 試走を実施した。**qsub、git commit、論文 `figures/` へのコピーは実施していない。** 原稿の数値もこの作業では書き換えていない。

## 保存先と再現情報

- 小規模結果：[results/local/pseudoloop-0.2.0-20261023/](../results/local/pseudoloop-0.2.0-20261023/)。`examples/`、`validation/`、`enumeration/`、`audit/`、`pilots/` に分けた。
- 新しい図：[Figure 1](../results/figures/pseudoloop-0.2.0-20261023/figure1/figure1.pdf)、[Figure 3（ファイル名 figure4.pdf）](../results/figures/pseudoloop-0.2.0-20261023/figure3/figure4.pdf)。対応する PNG、caption、validation、描画データと hash も保存した。
- 大規模計画：[array_plan.json](../results/runs/pseudoloop-context-0.2.0-20261023/array_plan.json)、[dry_run.json](../results/runs/pseudoloop-context-0.2.0-20261023/dry_run.json)。入力ごとの旧 profile 実測、旧出力の member/hash、全入力 hash、task 割当、完全な mode 引数が入っている。
- レビュー用の小さい設定と job script：[configs/pseudoloop_0.2.0_20261023/](../configs/pseudoloop_0.2.0_20261023/)。実行する script は `jobs/pseudoloop-context-0.2.0-20261023/*.sh`。同じ内容のコピーを設定の隣にも保存した。
- 数値ごとの一覧：[pseudoloop_0.2.0_20261023_inventory.tsv](pseudoloop_0.2.0_20261023_inventory.tsv)。Figure 2B の分位点、SI の Context TV セル、PseudoLoop の層別・局所集計、小規模の旧値と新値を含む。`pending array` は未計算であり、旧定義の数値を新定義に読み替えたものではない。
- 読んだ LaTeX の文脈関連行と hash：結果 directory の `manuscript_context_mentions.tsv`、`manuscript_review_hashes.json`。文献に出る context の一般的説明も含め、本文・SI の全 `.tex` を検索した。

ビルドは node-local の `$TMPDIR/pseudoloop-0.2.0-20261023/` だけで行った。

```sh
git -C ../../LinearKnot archive f93b55e | tar -x -C "$TMPDIR/pseudoloop-0.2.0-20261023/source"
cmake -S "$TMPDIR/pseudoloop-0.2.0-20261023/source" \
  -B "$TMPDIR/pseudoloop-0.2.0-20261023/build" \
  -DCMAKE_BUILD_TYPE=Release -DLINEARKNOT_ENABLE_LTO=OFF
cmake --build "$TMPDIR/pseudoloop-0.2.0-20261023/build" \
  --target linearknot scfg2_pk_profiles_oracle_test scfg2_pk_profiles_beam_test -j 4
```

Turner04、37°C、1 thread、weak closure 有効という採用条件を使用した。Release の C++ 最適化は `-O3`、LTO は無効。build cache、configure/build log、実行ファイル、パラメータは `runtime.tar.gz`、固定 commit のソースは `source.tar` にまとめた。各小規模計算の provenance に入力・実行ファイル・パラメータ・observer の hash、実際の argv、戻り値、時間、stdout/stderr の hash がある。生出力は tar、実行ごとのファイルは node-local に置いた。

完全なcommitは `f93b55ed9749dd1c198cd1f1809aece50f12f7d8`。実行ファイルSHA-256は `90a5a39d1f1294898b64460480066ffa512cab3afeaa81e47a4e50d1b9b6b92d`、Turner04は `0c671a7455f5686bc160fc655194b11888594ab7d506626cdad1470bcadd17d9`。ソースarchiveの705通常ファイルをnode-localの展開内容とbyte比較し、全件一致した（`build_provenance.json`）。列挙とループ別検証の埋込み入力には、sequence改行scaffold改行という明示した形式の入力hash一覧 `input_manifest.json` を加え、Turner04とDP09両パラメータのhashを記録した。計算時点のscript snapshotも `calculation_scripts.tar.gz` に残した。

最初の build 呼出しに存在しない observer target が含まれていたため最後が非ゼロ終了した。必要な三 target は完成しており、上記の実在する target だけの build を再実行して成功した。旧列挙 observer の初回 compile も、削除済みの RuntimeOptions field を指定して失敗した。新 API の `RuntimeOptions{2}` に合わせて作り直した。初回の stderr と command は `enumeration_initial_compile_failure/` に残した。実装ソースの編集ではない。

## 原稿・図・集計の対応

表番号は変わり得るため LaTeX label も併記する。

| 原稿・図の対象 | 旧値・旧出力 | 元実験・集計・作図 | 今回の対応 |
|---|---|---|---|
| Figure 1C の七つの文脈の帯、SI の61 nt例 | 最大 P=0.608737、0.2–0.8 の15位置 | `results/local/figure-feedback-20261020/data.json`、同 `raw_outputs.tar`、`scripts/prepare_feedback_figure_data.py` → `scripts/make_figure1.py` | 実際の採用入力 bpRNA_RFAM_31397 を再計算・再描画。古い `figure1-89baea52` config の PDB345 は別の例なので使わない。 |
| Figure 2B の白抜き記号と分布、Methods の caption、Results の文脈誤差 | 300–399 nt・b200：p90 0.0018、最大0.206（本文0.21）。長鎖 b300：p90 0.043 | E1 short/long の `metrics.tsv`、`quantiles.tsv`、`summary.json` → `scripts/render_e1_89baea52.py`、`make_figures_review.py`。採用図は `results/figures/e1_fb2_final_20261023/plot_data.tsv` | profile-only array 待ち。全長さ群・全幅の context 分位点を更新する。 |
| SI S13 `tab:si_accuracy_quantiles` の Context TV 列 | 8長さ群 × p90/最大の16セル。例0.00178、0.206 | `results/local/e1-bprna-short-89baea52-20261018/summary.json` と `e1-bprna-long-89baea52-cpu16-20261018/summary.json` | 16セルを inventory に保存。array 待ち。 |
| SI S14A `tab:si_pk_stratified` の Context TV 列 | 交差 mass の3群 × p90/最大の6セル | `scripts/analyze_e1_pk_stratified.py` → `results/local/e1-pk-stratified-89baea52-20261018/rna_error_quantiles.tsv` | 交差 mass の群分け自体は同じ。各群の context TV は再集計する。 |
| SI S14B の PseudoLoop base | 12,740 RNA / 223,512 bases、p90 4.01e-14、最大0.982 | 同 `local_error_quantiles.tsv`、`max_error_examples.tsv`。exact P≥0.1 の塩基を選ぶ | 新 P で選択し直すため、分母も再集計する。 |
| SI の具体例 bpRNA_RFAM_35940、塩基325 | b200 の P 誤差0.98191 | 同 `max_error_examples.tsv`、旧 E1 の exact/beam profile | この入力の全5 mode を今ここで計算。**新誤差0.9819320729907435（丸め0.98193）**。全cohortの最大値の代用にはしない。 |
| SI S15 `tab:si_local_errors` の PseudoLoop bases 全4行と直前の文 | b20：694/85 RNA、59/12 RNA。b50：446/56、38/6。b100：307/42、11/4。b200：80/15、0/0 | `scripts/analyze_e1_local_errors_89baea52.py` → `results/local/e1-local-errors-89baea52-20261019/summary.tsv` | 閾値は exact≥0.5/beam<0.1 とその逆。新 P に基づく件数・割合を再集計する。 |
| SI列挙 `tab:validation` と文脈の一意性の記述 | 729 cases、706 inputs、156,534 structures、context最大差2.14e-14 | `tools/validation_89baea52/enumeration/{enumerate.cc,run.py}` → `results/local/e0-enumeration-89baea52-20261018/` | 独立な局所ラベル表の七 rule を更新した observer を node-local でコンパイルし、同じ729 casesで再検証。結果は下記。 |
| SI `tab:si_byloop_validation` | Multiloop 7.8e-16、PseudoLoop 3.4e-15、finite reverse最大2.6e-15 | `scripts/validate_byloop_89baea52.py`、`results/local/byloop-validation-89baea52-20261018/summary.json` | 同じ20入力・21条件・294比較で再測定。下表。 |
| Figure 3B/C、Results のTercの文・caption | beam200 P窓0.840、約98%、context局所差<0.00284 | `results/local/small-rerun-89baea52-20261018/terc_data.tar`、`figure-feedback-20261020/data.json` → `make_figures_feedback_20261020.py` の Terc 描画 | exact、outer b200、h=d=200 の profile と by-loop を再計算・再描画。丸め値は維持。 |
| Tercの交差なしとの比較の交差あり側 | exact P窓0.842 | `scripts/analyze_terc_nopk_89baea52.py` → `results/local/terc-nopk-89baea52-20261019/window_accessibility.tsv` | 交差あり exact 側を再計算。交差なし ViennaRNA 側は保存値を保持。 |

現在の Abstract に context 誤差の数値はない。progress の旧一覧にある Abstract の数値は、今回読んだ現行ソースには存在しなかった。Methods/Introduction/SI の定義そのものは著者による新定義への変更済み。数値の差し替えと混同しない。

E1 の `context_tv_mean_Gfree`、`context_Multiloop_mae_Gfree`、`context_PseudoLoop_mae_Gfree`、`exact_pseudoloop_mean_Gfree`、PseudoLoop mass での層別記述、選択塩基数と局所誤差の全派生集計も更新対象である。他の五つの文脈の MAE は数学的には同じだが、新 profile を通して整合性を検証する。原稿で採用している層別は交差 pair mass であり、その群の所属は変わらない。保存済みの PseudoLoop mass 層別は所属も変わり得る。

Figure 2A/C、SI の total accessibility 図、BPP、分配関数、標本・代表構造の結果は今回のラベル変更の対象外。NEAT1 は現行原稿では資源表に集約され、文脈の確率を示す採用図・表はない。時間/RSS の表は旧版での測定記録を保持する。今回の試走時間を新しい性能測定として差し替えたり、ラベル変更だけで実測時間が同一と断定したりはしない。

## 小規模の新しい値

### Figure 1C（61 nt、exact）

| 指標 | 旧 | 新 |
|---|---:|---:|
| 最大 PseudoLoop 確率 | 0.6087370439248755 | 0.6087370439248755 |
| PseudoLoop 確率0.2–0.8の位置数 | 15 | 15 |
| PseudoLoop→Multiloop の最大移動量 | — | 0.0015273386655017885 |

七文脈の帯は変わる。最大値・位置数が保たれることと、全ての位置が同じであることは別である。入力 hash、実際の profile/access command は `examples/provenance.json`。旧版の by-loop も同じ入力で呼び出し、profile と窓の Multi+Pseudo 不変性を確認した（`audit/summary.json`）。

### Terc の窓123–126（joint probability）

| 設定・成分 | 旧 | 新 |
|---|---:|---:|
| beam200 total | 0.8587446945808029 | 0.8587446945808029 |
| beam200 Hairpin | 0.015546938018613003 | 0.015546938018613003 |
| beam200 Bulge | 0.001248899184414132 | 0.001248899184414132 |
| beam200 Internal | 0.0005668095987658231 | 0.0005668095987658231 |
| beam200 Multiloop | 0.0009810601969629854 | **0.0009811445354800372** |
| beam200 Exterior | 0.000765860056895999 | 0.000765860056895999 |
| beam200 PseudoLoop | 0.839635127525148 | **0.839635043186631** |
| exact total | 0.8605129487361832 | 0.8605129487361832 |
| exact Hairpin | 0.015197152844990562 | 0.015197152844990562 |
| exact Bulge | 0.0012511993160502113 | 0.0012511993160502113 |
| exact Internal | 0.0005659810829103324 | 0.0005659810829103324 |
| exact Multiloop | 0.0009665318230559554 | **0.0009666164964465066** |
| exact Exterior | 0.0009845759379804393 | 0.0009845759379804393 |
| exact PseudoLoop | 0.8415475077311951 | **0.8415474230578046** |

新 PseudoLoop / total は beam200 **0.9777469933558072**、exact **0.9779602088426061**。原稿の0.840、0.842、約98%は引き続き妥当。交差なしの total **0.12286507219770805**、Hairpin **0.11695022375854795** は既存 ViennaRNA の別計算を保持する。交差あり側を exact の新値に置き換えればよい。

Terc 全塩基での最大移動量は exact 0.0005585582660358329、beam200 0.0005572695592495。表示領域40–160の新 context の exact/beam 最大絶対差は **0.00283444428800448**。caption の <0.00284 は維持される。区間全体の非対合確率の旧新版差は、全ての要求窓で **0**。

実装 examples の説明にある塩基123の PseudoLoop marginal 0.8397007149 は、本表の4塩基の joint probability とは別の量である。

### SI の難例の再計算

bpRNA_RFAM_35940（335 nt、旧入力と byte 一致）、塩基325：

| 設定 | 旧 PseudoLoop | 新 PseudoLoop |
|---|---:|---:|
| exact | 0.9898543696644012 | 0.9898541636624613 |
| b200、h=d=200、outer | 0.007941112929212707 | 0.00792209067171774 |
| 絶対差 | 0.9819132567351885 | **0.9819320729907435** |

丸めた記載は **0.98193** となる。全5 mode（exact、b20/50/100/200）を保存した。

### ループ別 accessibility 検証

旧表と同じ20 distinct inputs、21 input–parameter conditions、320 enumerated groups。Turner04 20条件とDP09追加1条件、長さ9–29 nt、要求窓1/2/3/5/n/n+3/2（重複）、有効な start を使用した。

| exact 対 列挙：最大絶対差 | 旧 | 新 |
|---|---:|---:|
| Hairpin | 4.3298697960381105e-15 | 4.3298697960381105e-15 |
| Bulge | 2.220446049250313e-16 | 2.220446049250313e-16 |
| Internal | 7.771561172376096e-16 | 7.771561172376096e-16 |
| Multiloop | 7.771561172376096e-16 | **6.661338147750939e-16** |
| Exterior | 1.887379141862766e-15 | 1.887379141862766e-15 |
| PseudoLoop | 3.4416913763379853e-15 | 3.4416913763379853e-15 |

| その他の比較 | 新しい最大絶対差 |
|---|---:|
| endpoint-ordered unpruned 対 exact（21条件） | 3.4416913763379853e-15 |
| finite 対 taped reverse（252条件、両ranking/closure、b2/5/20） | **3.552713678800501e-15**（旧2.55351295663786e-15） |
| 六成分和 対 total（全294計算） | 2.220446049250313e-16 |
| 1塩基窓 対 profile（全294計算） | 3.3306690738754696e-16 |

長さ>1の positive window instances は Hairpin/Bulge/Internal/Multi/Exterior/Pseudo の順で **1835/83/55/105/669/1175**（旧Multi101、Pseudo1179）。個々の probability のラベル移動と、誤差最大値の浮動小数点の変化を区別する。

独立な pair-map 形状参照でも検証した。14 fixtures、両 closure それぞれ461 groups・10,565 base instances・19,987 fully-unpaired window instancesで一致し、変更される七 rule と六つの target 構造を全て通った。maintained exact/beam tests も成功した。

同じ729 cases・706 inputsの独立な局所ラベル表による列挙も全件成功した。156,534 structure instances、156,534 derivations、60,222 crossing structuresを維持し、文脈の曖昧さ・欠落・重複は0。context最大絶対差 **2.142730437526552e-14** は旧値と同じ。Multiloop個別の最大差は1.6601845766184287e-19 → **2.1345230270808369e-19**、PseudoLoopは1.554312234475219e-14のまま。入力一覧hashと各caseの生出力・commandを `enumeration/` に保存した。

さらに全156,534 structure instancesの完成pair mapを、ruleや局所重みを読まない独立な形状参照で再分類した。**2,898,623 base instancesでラベルが一致し、4,328,947 fully-unpaired window instancesで混合文脈は0**。結果と参照コードhashは `enumeration/geometry_validation.json`。これは保存した有限入力集合の検証であり、全入力の証明としては扱わない。

### 保存済み Terc 標本の確認

元の exact.samp と採用89baea52の exact.samp は同じ SHA-256 `41ecd697f944ceae1aca785e7c6da927d92ffdc98ebf74dc5e3c0a325987c2b7`。1,000 draws のうち窓が全て非対合の871 draws を新しい形状参照で分類すると **PseudoLoop848、Hairpin18、Bulge3、Multiloop2**。以前の progress メモの850/21は再現しなかった。この標本の joint PseudoLoop 頻度は848/1000で、848/871は非対合条件付きの割合。新たに sampling を行った結果ではなく、保存標本の再分類である。

## 旧新版の不変性と作図の検証

Figure 1、Terc exact/b200、21 pilot jobs、難例の5 modeで、各塩基の Multi+Pseudo が保たれ、他の五文脈が一致することを絶対許容差 **1e-9** で確認した。小規模例の各窓でも六成分和と total、旧新版の Multi+Pseudo、他四成分の不変性を確認した。実測のずれは丸め誤差の範囲で、許容差を広げていない。Stem は非対合成分和の引算なので最後のbitの違いを許す。

Figure 1 は178×145 mm、Figure 3は178×104 mm、最小文字7 ptを維持し、領域外の文字・文字同士の重なりの検査に合格。両PNGを開き、profile帯、窓の積上げ、軸・凡例・文字の欠けがないことを確認した。新しい Figure 3 の描画データには新 profile だけでなく新しい by-loop 全窓も反映した。BPP・total・samples・centroid は採用済み保存値を用い、その出どころを provenance に明記した。

## 大規模 E1 の計画

対象は短鎖 **24,625入力×5 mode=123,125 profile**、長鎖 **1,438入力×6 mode=8,628 profile**、計 **131,753 profile**。短鎖は exact および b20/50/100/200、長鎖は b20/50/100/200/300/400。全有限幅で outer ranking、h=d=max(b,30)。長鎖の参照は b400 のまま、exact とみなさない。温度、パラメータ、scaffold、seed20261011、warmup0/repetition1も旧E1と同じ。

旧 E1 の accessibility は **total のみ**であり、loop別誤差を評価していない。従って E1 では `linearknot profile` だけを実行する。`fold`、BPP、sampling、decoder、total accessibility の再実行は不要。`access --by-loop` が必要な採用内容は今回済ませた Terc と小規模検証である。

`plan` は旧 `array_plan.json` と全 `task_results.jsonl` の profile 観測を読み、各 input/mode の実測を対応づけた。全て旧成功観測が存在することを確認し、不足や失敗を捨てて cohort を縮めない。入力tarをhash確認し、各入力のbyte hash・sequence hash・長さも全件検査した。元の per-process timeout とメモリ上限は保持する。8/16 GiBの上限を保持する短鎖exactでは4/2 workers、長鎖beamは4 GiB×8 workersとし、task overheadを含めcpu_16の36.8 decimal GB内に収めた。

時間見積もりは各旧profile観測tに対して `d=2*(t+0.5秒)`、workersをpとして `60秒+(Σd/p+(1−1/p)*max(d))/0.8`。最悪1jobを含むlist-schedulingの上界に、版・node差への2倍余裕、起動・転送余裕、80%効率を加えた。1配列の全modeを同一taskに置き、1taskの見積もりはh_rtの65%以下とした。1配列だけで超えるtierはh_rtを5分刻みで引き上げた。全profileの過去実測の逐次合計は144.37時間であり、これ自体を今回の並列wall timeとは扱わない。

| 長さ(nt) | 入力 | profile | task | workers | h_rt | 推定割当h | 旧h_rt |
|---|---:|---:|---:|---:|---|---:|---|
| 50–99 | 12,428 | 62,140 | 44 | 4 | 00:15:00 | 7.014 | 02:00:00 |
| 100–199 | 9,656 | 48,280 | 44 | 4 | 00:15:00 | 7.076 | 03:00:00 |
| 200–299 | 1,446 | 7,230 | 13 | 2 | 00:30:00 | 3.963 | 03:00:00 |
| 300–399 | 1,095 | 5,475 | 21 | 2 | 00:30:00 | 6.728 | 08:00:00 |
| 400–799 | 584 | 3,504 | 13 | 8 | 00:30:00 | 3.936 | 01:00:00 |
| 800–1599 | 423 | 2,538 | 41 | 8 | 00:55:00 | 21.988 | 02:00:00 |
| 1600–3199 | 340 | 2,040 | 42 | 8 | 01:30:00 | 36.989 | 03:00:00 |
| 3200–4381 | 91 | 546 | 12 | 8 | 02:20:00 | 16.242 | 04:00:00 |

合計 **8 arrays、230 tasks、見込み104時間弱のcpu_16割当時間、h_rt合計174.083時間**。node_f換算（資源係数0.06）は見込み約6.24時間、h_rt合計10.445時間。通常優先度の概算ポイントは **5.41197**、全taskがh_rtまで使う場合 **8.35600**。

計算式は `0.06*(0.7*max(実行秒,300)+0.1*h_rt秒)/3600` をtaskごとに合計した。[2026年度の公式料金表](https://www.t4.cii.isct.ac.jp/tsubame-kyodo/tariff2026)でcpu_16係数0.06、最低5分、実行時間と指定時間の係数を再確認した。優先度0、allocation1、課金単位の切捨ては無視した概算で、再実行・他job・storageは含まない。

inode gate は有効。作成見積もり **3,020 entries**（計画directory、task metadata/lock、JSONL/tar、scheduler log、atomic-write余裕を含む）。最初のlive quotaは21,850,968 files、live上限28,000,000、保守上限24,000,000に対する残り2,149,032。新規計画は残りの10%以下。gateはdry-run時にも再照会し、詳細を `dry_run.json` に保存した。per-input/per-commandの共有ファイルは作らず、task当たり一つのJSONLと一つのtarに集約する。再投入は追加scheduler logを生む。

旧profile stdout総量から新出力の生stdoutは約 **2.72 GB**。tar metadata/512-byte padding、command/stderr/GNU timeと計画を含む共有領域増加の余裕は約 **5 GB**（6 GB確保を推奨）。固定runtimeと入力の展開、計算中の出力はnode-localのみ。各taskの出力は完了直後にtarへ追加・fsyncし、その後JSONLへ完了記録・fsyncする。resumeは完了済み成功/失敗を再実行せず、未完のtar/JSONL末尾を復旧する。task lockで二重controllerを拒否する。

### dry-run と試走

`scripts/pseudoloop_context_array_20261023.py check` は全26,063入力、runtime/source/script/planのhash、完全で重複のないcohort、memory/wall budgetを検査して成功。8 script全て `bash -n` 成功。

同じrunnerを別出力directoryで試走：50 ntを2入力×5 mode、300 ntを1入力×5 mode、400 ntを1入力×6 mode、計 **21/21成功**。全ての新profileをhash確認済み旧tar memberと比較し、旧新版の不変性に合格。50/300/400の最大RSSは3.0/66.0/24.0 MiB。全21観測が上記のper-job計画余裕 `2*(旧実測+0.5秒)` 内に入った。小さい例では起動時間の比が大きくなるため、ratioだけでは判断しない。短鎖pilotを `--resume` しても完了件数10を維持し、再実行がないことを確認した。

投入前の見直しで、失敗観測をJSONLに保存した後にrunnerも非ゼロ終了するよう修正し、planのrunner hashを固定し直した。初回planの圧縮snapshotと初回pilot/auditを保持した上で、最終runnerでも21件を再試走した。未投入の大規模taskは0件のまま。最終planの `pre_submission_runner_revision` と `previous_plan_sha256` がこの変更を記録する。

### 著者確認後の qsub コマンド（未実行）

まず投入直前に再検査する。

```sh
cd /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments
python3 -B scripts/pseudoloop_context_array_20261023.py check
lfs quota -p 102571 /gs/bs
```

確認済みの会計groupを `TSUBAME_GROUP` に設定して、次を実行する。各arrayの同時実行を3にすると、8 arrays合計は最大24 task（weekday上限30以内）。同時に動く別計画がある場合は著者側で総数を調整する。

```sh
qsub -g "$TSUBAME_GROUP" -tc 3 jobs/pseudoloop-context-0.2.0-20261023/nt50_99.sh
qsub -g "$TSUBAME_GROUP" -tc 3 jobs/pseudoloop-context-0.2.0-20261023/nt100_199.sh
qsub -g "$TSUBAME_GROUP" -tc 3 jobs/pseudoloop-context-0.2.0-20261023/nt200_299.sh
qsub -g "$TSUBAME_GROUP" -tc 3 jobs/pseudoloop-context-0.2.0-20261023/nt300_399.sh
qsub -g "$TSUBAME_GROUP" -tc 3 jobs/pseudoloop-context-0.2.0-20261023/nt400_799.sh
qsub -g "$TSUBAME_GROUP" -tc 3 jobs/pseudoloop-context-0.2.0-20261023/nt800_1599.sh
qsub -g "$TSUBAME_GROUP" -tc 3 jobs/pseudoloop-context-0.2.0-20261023/nt1600_3199.sh
qsub -g "$TSUBAME_GROUP" -tc 3 jobs/pseudoloop-context-0.2.0-20261023/nt3200_4381.sh
```

中断taskの例：子processの停止を確認してから `qsub -g "$TSUBAME_GROUP" -tc 3 -t 2-2 jobs/pseudoloop-context-0.2.0-20261023/nt300_399.sh --resume`。plan、runtime、runnerを変更するとhash検査で停止する。h_rtとresourceを提出時に上書きしない。

array完了後は、新旧profileのMulti+Pseudo/他五文脈の不変性を全件確認してから、新profile同士（短鎖exact、長鎖b400参照）でcontext TV/各文脈MAE/PseudoLoopの選択・局所誤差を再集計する。BPP・total accessibility・free energyの旧値は入力hashで結合する。旧 `analyze_e1_production.py` は三CLIを要求するため、profile-only出力にそのまま実行しない。Figure 2BとSIセルはその集計を待つ。現時点で原稿のこれらの数値を更新してよいとは結論しない。


## 2026-10-23 E1全230 taskの検証・集計と英日原稿への反映

上の「tagなし」「未投入・未コピー」は準備段階の記録。この節はarray完了後の結果である。現在の `v0.2.0` tagのpeeled commitは、再計算に使った `f93b55ed9749dd1c198cd1f1809aece50f12f7d8` と一致する。

- [新しい全集計](../results/local/e1-context-0.2.0-20261023/)に短鎖・長鎖のmetrics/quantiles/summary/provenance/validation、参照側のPseudoLoop massと自由塩基上の平均、層別・局所誤差を保存した。旧集計は上書きしていない。長鎖の参照は有限b400のまま。
- [新しいFigure 2](../results/figures/e1-context-0.2.0-20261023/figure2_main.pdf)、300 dpi PNG、plot_data/evaluation_counts/caption/manifest/validationを保存し、論文repoのFigure 2とexport manifestを更新した。
- `check`が成功。plan/runtime/source/runner/storage/job-script/旧plan/入力tarと、全入力のbyte・配列・長さのhashを確認した。完了taskは **230/230**、profileは **131,753/131,753**、入力は **26,063/26,063**。欠落・重複・失敗・timeout・不正出力は0。task metadata、JSONLの終端・mode行列、終了コード、実際の引数、binary/parameter/input hash、tar内全4 artifactのbyte数とSHA-256、commandとstdoutの対応を検査した。task tar自体のhashも保存した。
- 全 **27,208,998塩基分**で旧profileと比較した。Multiloop+PseudoLoopの最大絶対差は **6.464965520642124e-12**、Stemは **3.019806626980426e-14**、Hairpin/Bulge/Internal/Exteriorは **0**。計画どおり絶対許容差1e-9を使った。
- 短鎖 **98,500行**・長鎖 **8,628行**の文脈以外の全列が、旧metricsと文字列単位で一致した。BPP・total accessibility・ΔF・窓数・保存資源値は旧値のまま。新metricsは既存 `validate_e1_production.py` に合格した。
- 文脈は既存 `analyze_accuracy.compare_contexts`、集約は `analyze_e1_production.summarize` を使用。層別adapterは `analyze_e1_pk_stratified.py` の選択・誤差・分位点関数をそのまま使い、旧profileを通した全4表も全セル再現して照合した。局所誤差は既存 `analyze_e1_local_errors_89baea52.py` に新profileを重ねて集計した。閾値、NA、等RNA重み、線形補間の分位点は変更していない。
- all/交差pair mass群の文脈以外の派生集計も旧値と一致した。PseudoLoop mass層別では新しいmassで群を選び直したため、群の所属と、その群で条件づけたBPP/accessibility等の集計は変わり得る。RNA単位の非文脈量の変化とは区別した。

### 原稿に使う主要値（旧 → 新）

| 対象 | 旧 | 新 |
|---|---:|---:|
| 300–399 nt、b200、Context TV p90 | 0.001775133296871157 | 0.001735394573176407 |
| 同、本文・caption | 0.0018 | 0.0017 |
| 同、S13 | 0.00178 | 0.00174 |
| 同、Context TV最大 | 0.20602449748510265 | 0.20602439304774917（丸め0.206を維持） |
| 全長鎖、b300、Context TV p90 | 0.04306175698759455 | 0.04306177015534478（丸め0.043を維持） |
| S14B、exact PseudoLoop≥0.1のRNA / 塩基数 | 12,740 / 223,512 | 12,687 / 219,028 |
| 同、PseudoLoop局所誤差p90 | 4.007905118896815e-14 | 4.030109579389318e-14 |
| 同、PseudoLoop局所誤差最大 | 0.9819132567351885 | 0.9819320729907435 |

S15のPseudoLoop行（300–399 nt、全1,095本。各欄は事象数 / RNA数）：

| b | 取りこぼし：旧 → 新 | 作り出し：旧 → 新 |
|---:|---|---|
| 20 | 694 / 85 → 687 / 85 | 59 / 12 → 59 / 12 |
| 50 | 446 / 56 → 449 / 57 | 38 / 6 → 41 / 6 |
| 100 | 307 / 42 → 306 / 42 | 11 / 4 → 14 / 4 |
| 200 | 80 / 15 → 80 / 15 | 0 / 0 → 0 / 0 |

塩基325の例はarrayでも小規模の新誤差0.9819320729907435を確認し、既に反映済みの丸め0.98193を維持した。小規模のFigure 1・Terc・列挙・by-loop検証は先の採用値を保持した。inventoryの **pending array 1,874件を全て新値と新しい出どころで埋めた**。S13のContext TV 16セル、S14Aの6セル、S14BのPseudoLoop選択数・誤差、S15の全4行、本文とFigure 2 captionを英日で対応させた。

本文の「少なくとも90%が…以下」は、丸めたp90を厳密な上限として読ませないため、「p90 mean TV distances were …」とp90を直接述べる文にした。日本語も「平均TVのp90は…」にした。
交差mass≥1の群でb20/50/100の四指標のp90が他の二群より大きいという文は、全12比較で引き続き成立した。大小関係の逆転による意味の変更はない。

Figure 2のPNGを開き、Bの白抜き記号、軸、凡例、文字の重なり・欠けがないことを確認した。178×92 mm、300 dpi、最小7 ptを維持し、1,750分位点・350集計と文字・座標検査に合格した。文脈以外のplot_data **450行**と評価窓数は旧図と完全一致。E1のSI二図はtotal accessibilityとΔFのみを使い、新しいPNGが両方とも旧図とbyte一致した。採用SI図の差し替えは不要だった。export先と新Figure 2のPDF hashも一致した。

### ビルド・再現・実行記録

`cd manuscript/latex && bash build.sh` と `bash make_diff.sh 37a0058` が成功。通常本文 **8頁**・SI **29頁**、差分本文 **10頁**・差分SI **30頁**。比較元は `37a0058c21897bac7cc55489e89435dff929caa8`。LaTeXエラー・未定義参照・未定義引用は0。変更頁と差分の赤青表示を画像で確認した。PDF・ソースhashと警告は論文repo `manuscript/latex/build/e1_context_0.2.0_20261023_validation.json` に保存した。

実験repoからの再現手順（再実行時は新しい `--out` と対応するrenderer入力pathを指定する）：

```sh
python3 -B scripts/pseudoloop_context_array_20261023.py check
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 -B scripts/reanalyze_e1_context_0_2_0_20261023.py --scratch-dir /tmp --workers 4
python3 -B scripts/render_e1_89baea52.py --short-data results/local/e1-context-0.2.0-20261023/e1-bprna-short-context-0.2.0-20261023 --long-data results/local/e1-context-0.2.0-20261023/e1-bprna-long-context-0.2.0-20261023 --out results/figures/e1-context-0.2.0-20261023
```

新集計rootの `task_output_validation.json`、`context_update_validation.json`、`noncontext_and_claim_validation.json`、`manuscript_update.json`、`output_manifest.json` に照合結果とhashを記録した。実行時codeと転記・検証scriptは `analysis_scripts.tar.gz`、ログは `execution_logs.tar.gz` にまとめた。初回auditのscriptは `initial_analysis_scripts.tar.gz` に保存した。

開発上の記録：参照descriptorの初回exportは長鎖出力directoryがまだなく失敗した。両集計のvalidationを待つ順に直して成功した。逐次の層別集計は読み直しが多かったため、検証済みprofile DBをnode-localに保持して停止・再開し、tarのmemory cacheとtask並列adapterへ切り替えた。旧4表全セルの再現、新旧RNA指標と非文脈集計の一致で検証した。これらは集計処理の実行記録であり、arrayのprofile失敗や科学的な性能測定ではない。

共有inodeを確認し、per-inputのファイルを共有領域に展開せず、profile DBと中間ファイルはnode-localに置いた。95 MBの短鎖metricsを含む `results/local/` と生出力はgitの対象外。20 MB超の新規ファイルをstageせず、git commitも行っていない。外部実装repoは読み取りだけで、RNA CLIの再実行・qsub・性能測定は行っていない。
