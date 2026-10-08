# RNAcentral データ候補の確認（2026-10-08）

本メモは長い実配列による時間・最大 RSS 評価の準備を目的とする。最終入力の選定、scaffold の生成、アルゴリズム実行は行っていない。公開済みの小規模 panel と少数 record の metadata のみを取得して確認した。

## 判断

LinearCapR の公開 RNAcentral panel を取得できる。元の 20 配列のうち、ACGU だけからなる配列は 11 本で、4,793–337,259 nt をカバーする。特に 4,793–29,106 nt の 7 本は、初期の長さに対する時間・RSS 評価の候補にできる。最長 853,910 nt の実行を初期段階の要件にする必要はない。

ただし、残り 9 本は N を含む。N の置換・削除は行わず、ACGU-only を採用条件にする場合は除外する。この 11 本やその部分集合を利用する場合は「LinearCapR の公開 panel の ACGU-only subset」と記載し、元の 20 配列と同一の panel と呼ばない。

論文では「RNAcentral 登録配列」と記述する。登録された RNA type が lncRNA であることは確認できるが、完全な mature transcript や実験的構造注釈の保証としては扱わない。時間・RSS 評価に生物学的機能の解明を先に要求する必要はない。

## 固定した取得元

- 論文: [LinearCapR, Bioinformatics 2026](https://doi.org/10.1093/bioinformatics/btag295)。Supplementary Text S2 と Table S1 も確認した。
- 公開実験 repo: [LinearCapR_ComputationalExperiments](https://github.com/TakumiOtagaki/LinearCapR_ComputationalExperiments)。2026-10-08 の main/HEAD は `26313637784231088fdbaecbf5dc9e67f8889711`（`git ls-remote` と GitHub commit API で確認）。
- [固定 commit の panel FASTA](https://raw.githubusercontent.com/TakumiOtagaki/LinearCapR_ComputationalExperiments/26313637784231088fdbaecbf5dc9e67f8889711/data/processed/rnacentral/rnacentral_active_20.fasta): 3,437,375 bytes、20 records。
- [固定 commit の ID list](https://raw.githubusercontent.com/TakumiOtagaki/LinearCapR_ComputationalExperiments/26313637784231088fdbaecbf5dc9e67f8889711/data/processed/rnacentral/rnacentral_active_20_ids.txt): 20 IDs。FASTA と一致する。
- FASTA 全体の SHA-256: `14b264bc91744d58b1262ec17e4140f8ae4cc7ecf0401092b31fe24c14535cc5`。

Repo README は RNAcentral Release 25（公開告知 2025-06-19）を由来としている。RNAcentral の [Release 25 告知](https://blog.rnacentral.org/2025/06/rnacentral-release-25.html) と [25.0 archive](https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/25.0/) が存在する。ただし今回、release 全体の FASTA と panel 全 20 records の対応を独立照合してはいない。したがって、release 25 は著者の公開 provenance、上記 commit/hash は今回確認した配列固定根拠として区別する。Archive 内部の release_notes は 2025-05-14 を記載しており、告知日とは異なる。

元の候補範囲は 4,000–985,945 nt、実際の panel は 4,793–853,910 nt である。Supplementary Text S2 は候補範囲と最終標本の違いを説明し、Table S1 は後者を報告している。

## 近いうちに使える 7 本

すべて原 FASTA の配列文字列が ACGU のみである。下記 SHA-256 は、改行を除き uppercase にした RNA 配列文字列だけに対する値（FASTA header は含めない）。元 FASTA に T は存在せず、今回 T→U 変換は必要なかった。

出典列は **2026-10-08 に取得した現在 API の xrefs** であり、release 25 時点の metadata の復元ではない。各 ID の `/api/v1/rna/<ID>/xrefs` を用いた。

| RNAcentral ID | 長さ (nt) | 現在 API の source / 外部 ID | Taxid / species |
|---|---:|---|---|
| URS000233CDDE | 4,793 | LncBook: HSALNT0322046; GeneCards: LOC101927259 | 9606 / Homo sapiens |
| URS00019626CE | 5,482 | Ensembl: ENSSSCT00025043577 | 9823 / Sus scrofa |
| URS000239D8A4 | 6,990 | PLncDB: MTRU_LNC008266.1 | 3880 / Medicago truncatula |
| URS00023C9544 | 11,832 | PLncDB: CARA_LNC001626.1 | 13443 / Coffea arabica |
| URS000241043C | 13,652 | PLncDB: CANN_LNC016619.5 | 4072 / Capsicum annuum |
| URS000238C72C | 22,543 | PLncDB: CARA_LNC002490.4 | 13443 / Coffea arabica |
| URS00023A4DB9 | 29,106 | PLncDB: ZMAY_LNC003640.4 | 4577 / Zea mays |

| RNAcentral ID | 配列 SHA-256 |
|---|---|
| URS000233CDDE | `3d757ed7666e716acc54ee142430b5536fe1f222bd1e7a7ba97b15d1f913738a` |
| URS00019626CE | `fabeaedff0618b12b5c1a045bae6d6fa5719264430510ec9f2779701059b5848` |
| URS000239D8A4 | `d4fc8ae23a15761c6f002dc199811356b74dc2f9f0174d7d7697910077b0c618` |
| URS00023C9544 | `b8acfc8040f3c0c7703c0e7ec459ddf5d16af323395c664b1ebdbb214b956b17` |
| URS000241043C | `ff1758c9958a530c7b9244d52062185f425427b28da0e324e681291d358d5627` |
| URS000238C72C | `2237659dd3b8ffc0ecb4bad5e75ef5e9a8fce8d637069f69a8243654a9d59ad7` |
| URS00023A4DB9 | `843bf7a3108552b5b767449d65c6b2d1e2bbbc8b77dd2935473fc9be9860f9a5` |

残る ACGU-only 4 本は URS00023A038B（60,668 nt）、URS00023B9FFB（79,797 nt）、URS00023D6C38（134,714 nt）、URS000241DEE3（337,259 nt）である。長さ区間 [4,000,10,000)、[10,000,20,000)、[20,000,50,000)、[50,000,100,000)、[100,000,500,000)、[500,000,1,000,000) の件数は順に 3、2、2、2、2、0。ACGU-only subset は元の 20 logarithmic bins を完全には保たない。

## N と record の意味

N を含む 9 本は以下の通り。N 以外の ACGU 外文字は原 panel で検出しなかった。

| ID | 長さ (nt) | N の数 |
|---|---:|---:|
| URS000241010F | 17,009 | 1,513 |
| URS000234B597 | 37,723 | 170 |
| URS0002383EB1 | 84,464 | 4,001 |
| URS0002411078 | 160,505 | 12,887 |
| URS000238BD54 | 192,460 | 7,532 |
| URS000240503C | 249,870 | 10,469 |
| URS00024138E3 | 460,939 | 31,090 |
| URS00023C899B | 630,484 | 5,301 |
| URS000234CCEC | 853,910 | 79,344 |

RNAcentral は 10% 以下の N を含む配列も保持するため、「RNAcentral active」だけでは ACGU-only を保証しない。[RNAcentral FAQ](https://rnacentral.org/help/) は 10% を超える N を除外対象としている。Active は少なくとも一つの expert database に存在するという意味であり、完全性・実験的検証の保証ではない。[FTP documentation](https://www.rnacentral.org/help/ftp)

最長 record の [URS000234CCEC API](https://rnacentral.org/api/v1/rna/URS000234CCEC) と [xrefs](https://rnacentral.org/api/v1/rna/URS000234CCEC/xrefs) は、PLncDB の barley lncRNA `HVUL_LNC011203.1`（taxid 4513）を指す。最大 ACGU-only record の [URS000241DEE3 xrefs](https://rnacentral.org/api/v1/rna/URS000241DEE3/xrefs) は tomato の `SLYC_LNC000716.2` を指す。これらは genome として登録されている record ではない。ただし record ごとの成熟転写物の完全性や exon のつながりまでは確認していない。PLncDB は RNA-seq からの lncRNA annotation を含むデータベースである。[PLncDB 原論文](https://doi.org/10.1093/nar/gkaa910)

現在 API から取得した URS000234CCEC、URS000241DEE3、URS00023A4DB9 の 3 配列は、固定公開 FASTA の同じ ID と配列 SHA-256 が一致した。RNAcentral が返す MD5 も DNA 表記（U→T）に対する再計算と一致した。全 panel 内で配列 SHA-256 の重複はなかった。

## 選定方法と取得方法

Supplementary Text S2 は、ID/長さ TSV を用い、20 logarithmic bins から seed 42 で各 1 本を選び、その ID で配列を取得したと説明する。公開 [sample_log_bins_from_tsv.py](https://github.com/TakumiOtagaki/LinearCapR_ComputationalExperiments/blob/26313637784231088fdbaecbf5dc9e67f8889711/scripts/data_preprocessing/sample_log_bins_from_tsv.py) と整合する。

Repo README が挙げる `download_rnacentral_bins.py` は冒頭で不使用と説明している。今回見た公開 tree には README が指す候補 pool TSV が含まれていなかったため、seed だけから元の 20 IDs を再選定できるとは判断しない。公開済み FASTA/ID list を固定して再利用する分には、この再選定の不確定性は障害にならない。

| 案 | 利点 | 今回の採用条件 |
|---|---|---|
| 既存 panel の subset | 数 MB で配列を入手でき、固定 commit/hash がある。過去研究との入力の対応が明確。 | ACGU-only と長さ範囲を事前に決め、除外 ID を残す。subset の範囲で時間・RSS を評価する。 |
| RNAcentral から新しい panel | Species/source/RNA type/長さごとの反復数を目的に合わせて選べる。 | 候補 pool と metadata を先に保存し、ページ順に依存しない選定を実装する。結果を見て差し替えない。 |

少数の入力なら `https://rnacentral.org/api/v1/rna/<URS>?format=fasta` で配列、同 ID の JSON と `/xrefs` で source/taxid/external accession を取得できる。API は version 1 でも **データ release の固定ではない**。[API documentation](https://www.rnacentral.org/api) により version は API schema の互換性を指す。API から新たに取得する場合は取得日時、応答本文、ID、MD5、正規化後 SHA-256 を保存して snapshot を作る。

Release を固定した新しい候補 pool を必要とする場合には `current_release` を避け、[25.0/sequences](https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/25.0/sequences/) のような版付き archive を使う。全 active FASTA は 8.4 GB（gzip）のため、今回 download していない。Release 別の ID mapping と Rfam mapping も archive にある。今回、当環境の curl では FTP host への HTTPS 接続が拒否されたが、web tool で archive の listing は確認できた。既存 panel と少数 API は取得できている。

## 本研究側で保存する項目

- 元 dataset/release、URL、公開 repo commit、取得日、元ファイル SHA-256。
- URS、必要に応じた taxid、外部 accession と source version、RNA type、gene/family metadata、metadata の取得版または日時。
- 元配列、uppercase/改行除去/T→U の処理を明記した実験入力、長さ、GC、配列 SHA-256。
- ACGU-only に合わない record の除外理由。N の置換・削除を行わない。
- 正規化後の同一配列は hash で重複除去する。URS の違いや taxid の違いだけで独立した配列として数えない。近縁・同じ gene/family の偏りを制限する場合は別の規則として事前に定める。
- Dataset 由来の再配布情報。[RNAcentral は配列・annotation を CC0 で公開](https://www.rnacentral.org/downloads)しており、小規模な固定入力の再配布が可能。出典引用は付す。

Zenodo DOI `10.5281/zenodo.19450645` は LinearCapR v1.0.4 の software archive であり、今回確認した record は software ZIP 一つを含む。実験 panel の archive とは区別する。Panel の固定根拠は実験 repo の commit と SHA-256 とする。
