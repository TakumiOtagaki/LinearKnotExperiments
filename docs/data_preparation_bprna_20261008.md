# bpRNA 入力の取得確認と local 集合の保存（2026-10-08）

著者の指示により、`../LinearRaccess-paper` の別集合ではなく、LinearCapR の公開実験 repo の固定版から bpRNA 入力を取得する作業を始めた。
今回は、公開元と保存先の確認を先に行い、配列選定・構造処理・アルゴリズム実行は行っていない。

その後、著者から `../LinearRaccess-paper` にある集合を使う方針が示された。
`bpRNA_1m_90_filtered_clean.fa` を実験 repo の local 領域へコピーし、
SHA-256 `c76025522f1f9e174c25207595bc079e454deb656acbfd3a51de71c60c167bc5` の一致を確認した。
26,063配列、ACGU-only、ID の重複なし、配列の完全重複なしを確認した。
取得元・コピー先・hash・候補区間の件数は
[取得 manifest](../manifests/data-sources/bprna-local-20261008.json) に保存した。
この集合を LinearCapR の評価 subset と同一とは記述しない。

## 指定した公開元

- Repository: [TakumiOtagaki/LinearCapR_ComputationalExperiments](https://github.com/TakumiOtagaki/LinearCapR_ComputationalExperiments)
- Commit: `26313637784231088fdbaecbf5dc9e67f8889711`
- [固定版 README](https://github.com/TakumiOtagaki/LinearCapR_ComputationalExperiments/blob/26313637784231088fdbaecbf5dc9e67f8889711/README.md)
- [固定版 .gitignore](https://github.com/TakumiOtagaki/LinearCapR_ComputationalExperiments/blob/26313637784231088fdbaecbf5dc9e67f8889711/.gitignore)
- [完全 tree API](https://api.github.com/repos/TakumiOtagaki/LinearCapR_ComputationalExperiments/git/trees/26313637784231088fdbaecbf5dc9e67f8889711?recursive=1)

Tree API は指定 SHA を返し、`truncated=false` だった。
README が記載する次のファイルは、この tree に存在しなかった。

| README の入力パス | 確認結果 |
|---|---|
| `data/processed/bprna/bpRNA_multifasta.fasta` | tree に不在。固定 commit の raw URL は HTTP 404 |
| `data/raw/bpRNA/bpRNA_1m_90.fasta` | tree に不在。固定 commit の raw URL は HTTP 404 |
| `data/raw/bprna/bpRNA_1m_90.fasta` | .gitignore の小文字表記のパス。tree に不在 |

`.gitignore` は `data/processed/bprna/` 全体、raw bpRNA の FASTA/ZIP と `*.dbn` を除外している。
README は processed datasets が既に存在することを前提としている。
したがって、**この commit から公開されている bpRNA 入力ファイルの直接取得は完了できていない**。
総 record 数、入力 ID 一覧、原本 FASTA hash は未取得であり、0 件と記録しない。
公開済みの計算結果や実行成功例から入力集合を選ぶ処理も行っていない。

## 取得済み metadata の保存先

`results/local/data-sources/linearcapr-bprna-2631363-20261008/` に新規保存した。
この領域は `.gitignore` の `results/local/` に含まれる。

- `github-tree.json`: 完全な公開 tree 応答。
- `README.md` と `upstream.gitignore`: 指定 commit の原本 bytes。
- `source_manifest.json`: repository、commit、URL、bytes、SHA-256、入力パスの不在と 404 の記録。
- 前処理 script の原本 5 ファイル。
- `preprocessing_manifest.json`: 各 script の URL・bytes・SHA-256 と、確認できた前処理・不足情報。

既存の成果物は上書きしていない。
取得した script は読み取りだけに使い、実行していない。

## 公開前処理から確認できたこと

- [`filter_bprna_no_pk.py`](https://github.com/TakumiOtagaki/LinearCapR_ComputationalExperiments/blob/26313637784231088fdbaecbf5dc9e67f8889711/scripts/data_preprocessing/filter_bprna_no_pk.py) は TSV の `PageNumber` が文字列 `1` の rows を抽出する。
- [`extract_pagenumber1_dbn.py`](https://github.com/TakumiOtagaki/LinearCapR_ComputationalExperiments/blob/26313637784231088fdbaecbf5dc9e67f8889711/scripts/data_preprocessing/extract_pagenumber1_dbn.py) は、DBN 中に `#PageNumber: 1` という行があるファイルを指定 directory へコピーする。
- [`bpRNA_multifasta.py`](https://github.com/TakumiOtagaki/LinearCapR_ComputationalExperiments/blob/26313637784231088fdbaecbf5dc9e67f8889711/scripts/analysis/bpRNAanalysis/bpRNA_multifasta.py) は、既定の `data/processed/bprna/dbn_pagenumber1` にある DBN を filename 順に読み、名前と配列を `data/processed/bprna/bpRNA_multifasta.fasta` に書く。
- `dataset_utils.py` は DBN の sort と書き出し処理を提供する。`dbn2fasta.py` も同じ directory の DBN から個別 FASTA を書く。

ただし、入力 directory の構成ファイルと ID 一覧は非公開であり、これらの scripts だけで論文の入力集合との一致は証明できない。
確認した scripts は、Rfam-only 条件を独立に確定するものでもない。
前処理を推測して別集合を作り、それを「LinearCapR と同じデータ」と記述しない。

README は bpRNA-1m(90) Version 1.0 と説明しながら、Source は全体版の `http://bprna.cgrb.oregonstate.edu/bpRNA_1m.zip` にリンクしている。
この不一致も、元データの再取得手順を確定する際に解消する必要がある。
約 4 GB の公式 ZIP は取得していない。

## 著者の local 作業領域と利用条件

追加指示に従い、`/gs/bs/tga-satolab-gtex/otgk` 直下の directory 名だけを確認した。
名称に `capr` を含む directory は見つからず、相当する local repo の指定 FASTA を確認する段階には進めなかった。
大きな results の検索は行っていない。

固定 commit の tree に LICENSE / LICENCE / COPYING に相当する名称のファイルは確認できなかった。
データの再配布条件は未確定であり、仮に配列本体の Git 保存条件が不明でも、local 実験と取得 recipe・ID・hash の保存は別に進められる。
今回は配列本体そのものが取得できていないため、cohort や取得スクリプトはまだ作成していない。

## 今後同じ subset を照合する場合

同じ集合を使うには、著者の元 `bpRNA_multifasta.fasta` または対応する DBN directory/入力 snapshot の保存場所を特定するか、公開 archive にその入力が残っているかを確認する。
原本を取得した時点で URL または元の保存先、版・commit、bytes、hash、全 ID を保存する。
現在の pilot では著者の選択に従って保存済みの local 26,063配列集合を使い、
長さ区間による独立選定と G の保存を進める。
