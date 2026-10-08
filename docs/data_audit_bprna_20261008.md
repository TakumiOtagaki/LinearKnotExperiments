# bpRNA データ候補の調査（2026-10-08）

この記録は、短い実配列 pilot と exact に対する近似誤差の評価に使う入力候補の調査である。
配列の選定、配列本体の取得・コピー、scaffold の生成、アルゴリズムの実行は行っていない。
他論文の計算結果や、実行に成功した選定済み配列は再利用しない。

## 判断

`../LinearRaccess-paper/benchmark/data/bpRNA_1m_90_filtered_clean.fa` は、pilot の母集団として再利用可能な候補である。
現存ファイルは全 26,063 配列が ACGU のみで、ID と配列の完全重複はない。
50–99、100–199、200–399 nt の各区間には各 4 配列を選ぶのに十分な件数がある。
ただし、原本と現在の公式配布物の一致、当初の取得 URL・日付、RNA family の metadata は未確認であり、正式採用済みの snapshot としては扱わない。

今回の主評価は、同じ配列・同じ scaffold G・同じエネルギー条件における exact と有限 beam の差である。
参照構造の pseudoknot-free subset に限定する必要は、現段階ではない。
生物学的な予測精度を評価する際は、参照構造の由来と、入力 G の作り方を別途固定する。

## 公式配布元で確認した事項

- [公式 About](https://bprna.cgrb.oregonstate.edu/about.php) は、bpRNA-1m Version 1.0 を 2018-03-25 の初版とし、7 出典から集めた 102,318 の single-molecule secondary structures と説明している。
- [公式 Download](https://bprna.cgrb.oregonstate.edu/download.php) は、bpRNA-1m(90) を 90% 未満の sequence similarity の subset と説明している。
- Download page の HTML の `window.open` に埋め込まれた公式 subset URL は `https://bprna.cgrb.oregonstate.edu/bpRNA_1m_90.zip` である。通常のページ解析では、この JavaScript 内の URL が表示されない。
- この URL の HTTP HEAD は 2026-10-08 に HTTP 200 を返した。`Last-Modified` は `Mon, 02 Apr 2018 20:37:25 GMT`、`Content-Length` は **4,023,074,135 bytes**、`Content-Type` は `application/zip` だった。ヘッダのみ取得し、ZIP 本体は取得していない。更新日時と ETag は hash の代わりにはならない。
- 公式 Download は FASTA、BPSEQ、dot-bracket、structure type 等を提供すると説明している。全 bpRNA-1m 用の別形式 URL は HTML から確認できるが、90% subset の小さい FASTA 専用 URL は今回確認していない。
- [原論文（Danaee et al., 2018）](https://academic.oup.com/nar/article/46/11/5381/4994207) と [著者の bpRNA repository](https://github.com/hendrixlab/bpRNA) が原資料である。後者の README は、同じ filebase の BPSEQ、dot-bracket、structure type を対応させる説明を含む。

## 現存する local FASTA の独立確認

2026-10-08 に Python 標準ライブラリで FASTA を読み、ファイル hash、全 ID、正規化後の全配列、長さ、文字集合を集計した。

| 項目 | 確認結果 |
|---|---|
| ファイル | `../LinearRaccess-paper/benchmark/data/bpRNA_1m_90_filtered_clean.fa` |
| bytes | 5,517,177 |
| SHA-256 | `c76025522f1f9e174c25207595bc079e454deb656acbfd3a51de71c60c167bc5` |
| 配列数 | 26,063 |
| 文字集合 | A, C, G, U のみ |
| 重複 ID | 0 |
| 完全一致する配列の余分な records | 0 |
| 最小・中央値・最大長 | 50・102・4,381 nt |
| 平均長 | 193.13640026090627 nt |

| 長さ区間 | 配列数 |
|---|---:|
| 50 ≤ n < 100 | 12,428 |
| 100 ≤ n < 200 | 9,656 |
| 200 ≤ n < 400 | 2,541 |
| 400 ≤ n < 500 | 347 |
| 500 ≤ n < 1,000 | 327 |
| 1,000 ≤ n < 2,000 | 440 |
| 2,000 ≤ n < 5,000 | 324 |

pilot の候補区間は端点を重ねず、上の最初の 3 区間を使う。
区間ごとに 4 配列を選ぶことは可能だが、今回は ID の選定をしていない。

## フィルタと元 snapshot の記録

現存する `../LinearRaccess-paper/benchmark/scripts/filter_bprna_clean.py` の処理は次のとおりである。

1. FASTA の各配列行の前後 whitespace を除き、配列行を連結する。
2. uppercase と T→U の置換を行う。
3. 長さ 50 nt 未満を除く。
4. ACGU 以外の文字を含む配列を除く。
5. 残った record の header と順序を保存する。

このスクリプトには、参照構造の pseudoknot、RNA family、予測結果、実行成功の有無によるフィルタはない。
ID は `bpRNA_<source>_<number>` の形を保つが、元データベースの accession や RNA family の ID は FASTA header に含まれない。
スクリプトの SHA-256 は `459e15331bfdc0c2f2c89ebf4e47ae644585af2527b12ff5d0f0fc4ab6d43b3f` である。

既存の入力監査 `../LinearRaccess-paper/benchmark/results/bprna_audit_8661385/manifest.json` は、2026-09-14 に次を記録している。

| 項目 | 既存 manifest に残る記録 |
|---|---|
| 原本名 | `benchmark/data/bpRNA_1m_90.fasta` |
| 原本 bytes | 6,535,356 |
| 原本 SHA-256 | `dbc65331d6051141d7679235a87ef97e1ef450162926b79094d897517ad93cb1` |
| 原本の配列数 | 28,370 |
| 長さ不足による除外 | 1,353 |
| 50 nt 以上で alphabet 不適合による除外 | 954 |
| retained | 26,063 |
| フィルタ再現 | 内容と順序の一致を当時確認したとの記録 |

現在の clean FASTA hash とフィルタスクリプト hash は、この manifest の hash と一致する。
監査 manifest 自体の SHA-256 は `9c411b85b871f8ffdf48c39c506d9e52a036b59e86f2b6975da2172aff2f6aa4` である。
原本 FASTA は、今回確認した現在の `benchmark/data/` には存在しない。
したがって、原本からのフィルタ再現は今回は行っておらず、上表の原本情報と除外件数は既存記録の確認である。
その manifest は retrospective audit であり、元の download URL、取得日、ZIP hash を確定する記録ではない。

## Source と RNA family、参照構造

local FASTA の ID に含まれる source の内訳は以下だった。

| source | records |
|---|---:|
| CRW | 3,420 |
| tmRNA | 235 |
| SRP | 333 |
| SPR | 207 |
| RNP | 246 |
| RFAM | 21,557 |
| PDB | 65 |

これらは **データの出典区分**であり、RNA family の分類ではない。
`RFAM` という文字列から特定の Rfam family accession を復元することもできない。
[公式 Help](https://bprna.cgrb.oregonstate.edu/help.php) は、bpRNA ID が出典を示すこと、RNA type が元の機能注釈に基づくことを説明している。
Family を単位とした選定上限を置くには、元 accession、Rfam family、RNA type 等を区別した metadata の join が必要である。

今回確認した local FASTA 自体には構造注釈がない。
公式は ID/filebase に対応する構造ファイルを提供すると説明しているが、local 26,063 件と公式 BPSEQ/DBN/ST の全件照合は未実施である。
現在の local 集合について、参照構造の pseudoknot 有無、構造の取得・検証方法、canonical pair の割合、モデルが許す構造クラスへの適合は未確認である。
既存の clean フィルタには pseudoknot を除く処理がない。

この集合は 7 出典を含むため、LinearCapR 論文の **pseudoknot-free Rfam subset と同じ集合とは扱えない**。
同じ bpRNA という名称だけで、入力数や性質の比較をしない。
その論文の subset の正確な選定条件・ID 一覧との一致は、今回の調査では確認していない。

また、公式 Help は構造の由来に comparative analysis、実験構造、予測・unknown の区別があると説明している。
bpRNA 全体を、実験で直接決まった構造だけの集合として記述しない。

## 利用・再配布条件

公式 About / Download / Help と著者 repository の表示を確認した範囲では、データ全体に適用される明示的な license と再配布条件を確定できなかった。
公開 download があることや原論文の open access 条件だけで、派生 FASTA の Git 再配布条件が確定したとは扱わない。
今後、ZIP 内 README 等と出典側の条件を確認する。
それまでの再現用成果物は、取得手順、必要な原本 hash、選定 ID、フィルタ設定、正規化後の配列 hash を中心に準備できる。

## 次に実行できる手順

1. 原本 FASTA の既存保存場所または元取得記録を確認する。公式 ZIP は約 4.02 GB あるため、小さい公式配列ファイルで照合できるかも先に調べる。
2. 取得物の URL・日時・bytes・hash を保存し、原本 FASTA を上記の既存 hash と照合する。一致しない場合は変更点を調べ、別 snapshot として記録する。
3. 固定した原本に明記したフィルタを適用し、現存 clean FASTA の hash と比較する。ここまでで local 入力の由来を確定する。
4. 主評価では配列長・alphabet・重複を事前条件とし、3 区間から seed と決定的な ID 順序を保存して各 4 配列を独立に選ぶ。既存論文の結果や成功例を選定条件にしない。
5. Family metadata を用いる場合は、その取得・欠測・join 条件を選定前に固定する。Source を family の代用にしない。
6. 選んだ配列の ID・長さ・GC・配列 hash と候補母集団・選定 script/seed を保存し、同じ G と energy 条件での exact 比較を準備する。構造注釈を用いる別評価は、必要な join と構造クラスの確認後に準備する。

この調査では、母集団の由来と容量、local ファイルの同一性、選定可能な長さ区間まで確認した。
正式な入力 manifest と pilot の選定 ID は、次のデータ準備作業で作成する。
