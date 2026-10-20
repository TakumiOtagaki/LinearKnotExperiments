# γ-centroid の計算準備（2026-10-20）

採用版 `89baea52c9435b18829ef0867a66a24be56d70f9` の既存の専用 build を使う。
NEAT1、ラット Terc、Figure 1 の MEA を γ-centroid に置き換える設定を作った。
E0 の列挙との比較には centroid を追加した。
本番の centroid matrix と E0 全件は未実行で、qsub・commit は行っていない。
実装 repo と既存の rerun plan は変更していない。

主設定は `--gamma 1` とする。
γ-centroid は、真陽性の対の数を TP、真陰性の対の数を TN としたときの期待値 `E[γ TP + TN]` を最大にする。
TN はすべての位置の組 `i<j` を対象とする。
対 `(i,j)` を予測に加える増分は、その対の確率を `p(i,j)` として `(γ+1)p(i,j)−1` である。
**γ=1 が通常の centroid**であり、期待される塩基対の対称差を最小にする。
固定した G の対は保持する。
MEA は対と非対合塩基を採点する別の目的関数なので、MEA の γ=1 と centroid の γ=1 が同じ構造を返すとは限らない。
定義は固定版の `src/scfg2/api/centroid.hh` と `src/scfg2/decoding/pair_posterior.hh` で確認した。

## 設定と計画

[準備 manifest](preparation.json) に config・array config・plan・script と依存ファイルの hash を保存した。
[file budget](file_budget.json) は作成・実行を合わせて500 entriesの余裕を見込み、既存 planner の inode 検査も有効にした。
すべての matrix plan は実際の source/build/lock/入力を照合済みである。
新しい plan は4 task・24 CLI jobで、そのうち centroid は9 jobである。
別に E0 の1 taskを準備した。

| 対象 | config（experiments/ からの相対パス） | 条件 | 実行 |
|---|---|---|---|
| NEAT1 | `configs/centroid_89baea52_20261020/e3-neat1-centroid-89baea52-20261020.json` | 両 isoform、b=h=d=100/200/400、outer、γ=1、stats | 2 tier、各1 array task、各 taskで3幅 |
| ラット Terc | `configs/centroid_89baea52_20261020/terc-centroid-89baea52-20261020.json` | 同じ388 nt/G、exact と b=h=d=200、γ=1 | local、1 task |
| Figure 1 | `configs/centroid_89baea52_20261020/figure1-centroid-89baea52-20261020.json` | 同じ PDB345/G、exact、γ=1 | local、1 task |
| E0 decoding | `configs/centroid_89baea52_20261020/e0-decoding-centroid-89baea52-20261020.json` | 同じ729件、exact/b20/b100、γ=.5/1/2 | 1 array task、localも可 |

Terc と Figure 1 は、既存の6出力のうち MEA だけを centroid に差し替えた。
BPP/MFE、profile、access、sampling と seed、Subopt、入力、scaffold、parameter は元の設定を引き継ぐ。
Terc の by-loop は既存の `terc-byloop-89baea52-20261018` を使う。
出力は各新 run の `tasks/<tier>/000001/run/task_results.jsonl` と `task_outputs.tar` に保存する。
共有領域に CLI ごとのファイルを作らず、実行中のファイルは node-local TMPDIR に置く。
この準備では図の描画・原稿の数値や MEA 表記を更新しない。
既存の `analyze_small_rerun_89baea52.py examples` は旧 run 名と MEA を前提とするので、この新出力の集計・描画には後続の対応が必要である。

## 資源と見積もり

NEAT1 は前回の MEA と同じ `cpu_40=1`、3 workers、各20 GiB RLIMIT_AS とした。
メモリの計画上限は `3×20480+1024=62464 MiB` で、割当の87738 MiB内に収まる。
時間の根拠は `results/runs/e3-neat1-mea-d02d261c-20261017/tasks/*/000001/run/task_results.jsonl` の同じ入力・同じ幅の実測で、hashを array config に記録した。
centroid 自体の実行時間を測った値ではない。

| 対象 | 前回 MEA b100/200/400（分） | 前回の逐次合計×1.05（分） | 新 h_rt | プロセス timeout |
|---|---:|---:|---|---:|
| NEAT1_1 | 1.72 / 4.46 / 10.48 | 17.50 | 00:22:00 | 1260秒 |
| NEAT1_2 | 12.07 / 33.24 / 82.69 | 134.40 | 02:50:00 | 10080秒 |

h_rt は、それぞれ前回の最も遅い幅の約2.10倍・2.06倍である。
planner は並列による速度向上を仮定せず、実測の逐次合計に5%を加える。
これが h_rt の80%以下に収まることを確認した。
3幅は既存の blocked-random order と3 workersで実行する。
centroid の時間や必要メモリは未校正なので、最初の実行の完了率と RSS を確認する。

Terc は保守的に24分、Figure 1 は12分の逐次仕事量を見込む。
local の各プロセスは8 GiB・1800秒で、local elapsed timeを正式な性能測定に使わない。
E0 は前回約19.7分に centroid API/CLI と独立再採点を加えて約40分を見込み、`cpu_16=1`・16 GiB・h_rt 2時間の1 taskとした。
NEAT1 の要求 task 時間は計3時間12分、E0 と合わせて5時間12分である。

## E0 に追加した照合

`tools/validation_centroid_89baea52_20261020/decoding/` は既存の検証の別コピーである。
列挙器、729入力、MEA/MFE、sampling の18条件・計360万 draw、seed、許容差1e−9と検定条件を引き継ぐ。
exact と保持した beam ensemble それぞれの列挙分布から、全 `i<j` に対する期待 gain を直接合計する。
その列挙最適値、予測構造の再採点値、API の score、定数と対の増分による別計算を照合する。
同じ目的関数で CLI の `EG` と構造も確認する。
tie は構造のbyte一致ではなく、列挙した集合への所属と目的値で判定する。

centroid は `729×3 mode×3 gamma=6561` 目的値照合と6561 CLI レコード照合を追加する。
MEA と MFE の CLI を含む合計は15309レコードである。
全件の検証結果はまだ存在しない。
準備確認として3つの小入力だけで27 API目的値・27 CLI照合を実行し、最大誤差は2.85e−14未満だった。
確認結果と初期失敗を含むログは [LinearKnot 固定の報告](../linearknot_9c7d6e1_20261020/README.md) から参照できる。

## 後で使う実行コマンド（未実行）

この一式は従来の MEA 版の対応 run の代わりに使う。
旧一式の `submit_arrays.sh` には自動で登録していない。
他の array と合わせた task の同時実行枠を30以下に保つ。
group は認可済みの値を設定する。

```bash
cd /gs/bs/tga-satolab-gtex/otgk/LinearCPartyPaper/experiments
: "${TMPDIR:?node-local TMPDIR is required}"
export PYTHONDONTWRITEBYTECODE=1

# Local examples (new plans; add --resume only after an interrupted attempt).
bash jobs/terc-centroid-89baea52-20261020/local.sh
bash jobs/figure1-centroid-89baea52-20261020/local.sh

# Scheduler submissions: examples only; none has been submitted.
: "${TSUBAME_GROUP:?Set the authorized accounting group}"
qsub -g "$TSUBAME_GROUP" -tc 1 jobs/e3-neat1-centroid-89baea52-20261020/neat1_short.sh
qsub -g "$TSUBAME_GROUP" -tc 1 jobs/e3-neat1-centroid-89baea52-20261020/neat1_long.sh
qsub -g "$TSUBAME_GROUP" -tc 1 jobs/e0-decoding-centroid-89baea52-20261020/decoding.sh

# Alternatively, E0 may run locally in a suitable allocation:
# python3 -B scripts/run_centroid_validation_89baea52.py
```

E0 controller は準備時の依存 hash と software lock を開始・終了時に確認し、ログをtarにまとめる。
NEAT1 が終了したら新 plan の `aggregate_array.py --plan` で root の集計を更新する。
local の例も同じ packed runner と checkpoint/resume を使う。
準備を再検証するコマンドは `python3 -B scripts/prepare_centroid_89baea52.py` である。
このコマンドは config と plan を検査・準備するだけで、RNA計算や qsub を実行しない。
