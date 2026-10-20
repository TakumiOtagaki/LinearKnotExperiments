# LinearKnot 9c7d6e1 の固定と γ-centroid の準備

正規の repo の固定版 `9c7d6e170b2b21a94a2e626d9939817e62898253` を
`deps/LinearKnot-9c7d6e1` に clone し、`build/linearknot-9c7d6e1` に Release build を作った。
専用 clone の origin は `https://github.com/TakumiOtagaki/LinearKnot.git` である。
提供元の `../LinearKnot` と `../LinearCParty` は変更せず、qsub・git commit は行っていない。

[release receipt](../../manifests/release-preparation-linearknot-9c7d6e1.json) と
[software lock](../../manifests/software-linearknot-9c7d6e1.json) に source commit/tree、clean state、parameter・CLI・cache hash、compiler、configure arguments を保存した。
`linearknot --version` は `LinearKnot 0.1.1` を返す。
configure は `-DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DLINEARKNOT_ENABLE_LTO=OFF`、build は `-j8` とした。
LTO OFF は固定版の `test/tools/README.md` の acceptance-suite 手順に合わせた。
receipt には cache の設定に加え、実際の CLI compile flags と link command も保存した。

維持された `test/tools/run_required_scfg2_tests.sh` は source root から実行し、87項目と dependency gate が成功した。
専用 build への一時的な `source/build` symlink は終了時に削除した。
初回は実行場所を experiments にしたため相対 parameter path が読めず停止したが、source root から suite 全体を再実行して成功した。
初回を含むログは [preparation_logs.tar.gz](preparation_logs.tar.gz) に保存した。
lock 作成時の3つの maintained checks、lock 回帰テスト2件、新しい準備テスト3件も成功した。

開始時の共有領域は21,567,719 / 28,000,000 inode、残り6,432,281だった。
固定 clone は833 entries、build は2369 entriesを使用した。
比較の入力とプロセスごとの出力、一時ビルド、テストログは TMPDIR に置き、共有領域に残す raw はtarにまとめた。

## CLI の対応

固定版の README と、旧7 CLI・新7 subcommand の `--help` をそれぞれ確認した。
Usage のコマンド名を除く help の全行は一致した。
**今回使う option 名に変更はない。**
[機械可読の対応表](cli_mapping.json) とログの `cli-help.json` に根拠を残した。

| 89baea52 の CLI | LinearKnot | 主な option（変更なし） |
|---|---|---|
| `LinearCParty` | `linearknot fold` | `--bpp FILE --bpp-cutoff 0` |
| `LinearPKCapR` | `linearknot profile` | `-o/--output` |
| `LinearPKRaccess` | `linearknot access` | `--access-len 1,4,8,16 --probability --by-loop`、別名 `--byloop/-byloop` |
| `LinearPKSampling` | `linearknot sample` | `-k/--count --seed --unique` |
| `LinearPKSubopt` | `linearknot subopt` | `-k/--top-k -e/--delta-energy --stats` |
| `LinearPKMEA` | `linearknot mea` | `--gamma --stats` |
| `LinearPKCentroid` | `linearknot centroid` | `--gamma --stats` |

共通の入力 option `-i/--input`、`--sequence/--scaffold` と、`-b/--beam`、`--exact`、
`--max-hairpin`、`--max-case3-depth`、`--ranking`、`--outer-pilot-exponent` も保持されている。
内部の互換環境変数 `CPARTY_FIXED_ENERGY_PARAM_FILE` も同じである。

## 旧版との比較方法と保存先

`manifests/release-preparation-89baea52.json` に記録された14入力、50–390 ntを使った。
元の `results/runs/e1-bprna-short-d21618c5-20261014/inputs.tar` から必要な14件だけを TMPDIR に取り出し、各入力の SHA256 を既存 receipt と照合した。
reference は全 mode で `build/linear-cparty-89baea52`、candidate は新 build である。
parameter bytes は両 source で一致した。

b=20 と b=200 は `h=d=max(b,30), ranking=outer`、exact は `--exact` のみとした。
fold は cutoff=0 の BPP file も比較した。
access は窓1/4/8/16の total と `--by-loop` の両方を実行した。
sample は count=100・seed=20261016を双方に与え、通常形式と `--unique` の両方を比較した。
MEA/centroid は γ=1、追加の Subopt は top-k=8である。
各プロセスは単一 thread、timeout180秒とした。初回の RLIMIT_AS は8 GiBで、旧版の exact sampling 2件が上限に達したため、その2条件だけを双方16 GiBで再実行した。再現スクリプトの既定上限は16 GiBに更新した。
`--stats` は比較では指定していない。

比較は `14 inputs × 3 modes × 9 outputs = 378条件、756 processes` である。
stdout・stderr・BPPの比較では **`# elapsed_seconds` で始まる行だけ**を除いた。
brand 名の違いを strict byte 比較から除いていない。
access の次の2行だけを対応させた別の照合も行う。

```text
# linearpkraccess_metadata=1  →  # linearknot_access_metadata=1
# engine=linearcparty-scfg2  →  # engine=linearknot-scfg2
```

初回の8 GiB上限では旧版の exact sampling 2件が `std::bad_alloc` で停止した。
その2条件を双方16 GiBで再実行すると4プロセスとも成功した。
各条件の最終実行では全378条件・756プロセスが正常終了した。
実際の実行回数は760で、初回の失敗2件も保存した。timeoutは0だった。
経過時間行を除く strict byte 一致は **294/378条件**だった。
数値・構造・BPP・sampling の出力は一致し、差は access の名称2行だけだった。
2行の名称だけを対応させた照合は **378/378条件**で一致した。

| 出力 | strict byte 一致 | 名称2行を対応させた一致 |
|---|---:|---:|
| fold | 42/42 | 42/42 |
| profile | 42/42 | 42/42 |
| access | 0/42 | 42/42 |
| access-byloop | 0/42 | 42/42 |
| sample | 42/42 | 42/42 |
| sample-unique | 42/42 | 42/42 |
| mea | 42/42 | 42/42 |
| centroid | 42/42 | 42/42 |
| subopt | 42/42 | 42/42 |

集計は [comparison_summary.json](comparison_summary.json) を参照する。
全条件のコマンド、戻り値、差の例と hash は `results/local/linearknot-9c7d6e1-comparison-20261020/comparison.json`、
全760実行の記録は同じ場所の `processes.jsonl`、raw stdout/stderr/BPP は `raw_outputs.tar` に保存する。
初回全体は隣の `linearknot-9c7d6e1-comparison-20261020-initial8g/` に保存した。最終tarでも差し替え前の sampling 出力を `initial8g/` 以下に残し、同名memberの重複を避けた。
source、CLI、CMake cache、parameter の hash が比較の前後で変わらないことも検査する。

この比較は14入力と上記条件の回帰確認であり、すべての配列・option に対する一致を証明するものではない。
build設定が異なり、比較中にテストも実行したため、保存した elapsed time を正式な性能測定や旧論文値の置き換えに使わない。
採用版89baea52の論文数値を取り直す既存の判断は変更しない。

## γ-centroid の計算準備

[専用 README](../centroid_89baea52_20261020/README.md) に全 config・plan・script・投入/localコマンドと見積もりを示した。
新しい CLI で既存の全 runner を変更する範囲は含めず、今回の centroid 計算には採用版89baea52の build を使う。
NEAT1 は γ=1、b100/200/400、同じ cpu_40/3 workers/20 GiB per processで準備した。
h_rt は前回 MEA の最長実績の約2倍、短 isoform 22分・長 isoform 2時間50分である。
Terc は exact/b200、Figure 1 は exact の local packed taskを用意した。
E0 は729入力・3 mode・γ=.5/1/2に centroid の列挙最適値とCLI照合を追加し、1 task・h_rt 2時間で準備した。

γ=1 は通常の centroidである。
MEA と目的関数が異なるので、同じ γ 数値が同じ構造を意味しない。
全件実行の前に、小入力3件で新しい E0 probe のコンパイルと27 API目的値・27 CLI照合を確認した。
[小入力の確認結果](centroid_smoke.json) の最大誤差は2.842170943040401e−14である。
本番 matrix と E0 全件の結果はまだ生成していない。

## 再現コマンド

experiments/ から、新しい出力先で実行する。
clone/checkout は専用の deps だけに行う。

```bash
export PYTHONDONTWRITEBYTECODE=1
: "${TMPDIR:?node-local TMPDIR is required}"
git clone --no-hardlinks --no-checkout /gs/bs/tga-satolab-gtex/otgk/LinearKnot deps/LinearKnot-9c7d6e1
git -C deps/LinearKnot-9c7d6e1 checkout --detach 9c7d6e170b2b21a94a2e626d9939817e62898253
cmake -S deps/LinearKnot-9c7d6e1 -B build/linearknot-9c7d6e1 \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON -DLINEARKNOT_ENABLE_LTO=OFF
cmake --build build/linearknot-9c7d6e1 -j8
# Run the maintained suite from source root, with SCFG2_BUILD_DIR set to
# the absolute dedicated build; temporarily link source/build to that build.
python3 -B scripts/make_software_lock.py \
  --software-root deps/LinearKnot-9c7d6e1 --build-dir build/linearknot-9c7d6e1 \
  --cli-layout unified --configure-argument=-DCMAKE_BUILD_TYPE=Release \
  --configure-argument=-DBUILD_TESTING=ON --configure-argument=-DLINEARKNOT_ENABLE_LTO=OFF \
  --required-suite 'Maintained suite passed from source root; retain its log hash' \
  --out manifests/software-linearknot-9c7d6e1.json
python3 -B scripts/compare_linearknot_outputs.py \
  --out results/local/linearknot-9c7d6e1-comparison-REPRODUCTION
python3 -B scripts/summarize_linearknot_comparison.py \
  results/local/linearknot-9c7d6e1-comparison-REPRODUCTION --out "$TMPDIR/linearknot-summary.json"
python3 -B scripts/prepare_centroid_89baea52.py
```
