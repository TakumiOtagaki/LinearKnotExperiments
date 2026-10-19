#!/usr/bin/env python3
"""Validate and document the saved Terc crossing/no-crossing comparison."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import runpy
import tarfile

import RNA

ROOT = Path(__file__).resolve().parents[1]
ALLOWED, NO = "crossings_allowed_exact", "no_crossings_exact"


def read_tsv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    out = args.directory
    energy = {r["ensemble"]: r for r in read_tsv(out / "ensemble.tsv")}
    window = {r["ensemble"]: r for r in read_tsv(out / "window_accessibility.tsv")}
    context = read_tsv(out / "context_region.tsv")
    base = {r["ensemble"]: r for r in context if r["position"] == "123"}
    regions = {(r["ensemble"], r["region"]): r for r in read_tsv(out / "region_summary.tsv")}
    pairs = read_tsv(out / "reference_pairs.tsv")
    validation = json.loads((out / "validation.json").read_text())
    provenance = json.loads((out / "provenance.json").read_text())
    assert provenance["dependencies_sha256_before"] == provenance["dependencies_sha256_after"]
    assert provenance["source_status_before"] == provenance["source_status_after"] == ""
    assert len(context) == 242 and len(pairs) == 18
    # A zero maximum has no supported partner. Render the identifier as blank.
    for r in pairs:
        for side in ("left", "right"):
            if float(r[side + "_top_partner_probability"]) == 0:
                r[side + "_top_partner"] = ""
    with (out / "reference_pairs.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, list(pairs[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(pairs)
    assert all(float(r["pair_probability"]) == 0 for r in pairs if r["ensemble"] == NO)
    for r in context:
        exact_mass = sum(float(r[c]) for c in
                         ("Stem", "Hairpin", "Multiloop", "Exterior", "PseudoLoop", "Bulge_or_Internal_exact"))
        assert abs(exact_mass - 1) < 1e-9
        if r["Bulge"] and r["Internal"]:
            assert abs(float(r["Bulge"]) + float(r["Internal"]) - float(r["Bulge_or_Internal_exact"])) < 1e-10
    with tarfile.open(out / "raw_outputs.tar.gz") as archive:
        seq, g = archive.extractfile("in.txt").read().decode().splitlines()[1:3]
        summary = json.loads(archive.extractfile("summary.json").read())
        counts = json.loads(archive.extractfile("conditional_counts.json").read())
        calls = json.loads(archive.extractfile("pf_calls.json").read())
        classify = runpy.run_path(str(ROOT / "scripts/analyze_terc_nopk_89baea52.py"))["classify"]
        native_context_checks = 0
        for filename in ("unconditional.samples.txt", "base123_BI.samples.txt", "window123_126_BI.samples.txt"):
            stream = archive.extractfile(filename)
            for _ in range(1000):
                db = stream.readline().decode().strip()
                assert len(db) == len(seq)
                ours = classify(db)[1:]
                native = RNA.db_to_element_string(db)
                for a, b in zip(ours, native):
                    expected = {"e": "Exterior", "h": "Hairpin", "m": "Multiloop"}.get(b)
                    if expected:
                        assert a == expected
                    elif b == "i":
                        assert a in ("Bulge", "Internal")
                    else:
                        assert a == "Stem"
                native_context_checks += 1
        pt = RNA.ptable(g)
        leaf_bounds = []
        multiple_child_bounds = []
        for i, j in enumerate(pt):
            if not i or j <= i:
                continue
            children, k = [], i + 1
            while k < j:
                if pt[k] > k:
                    children.append((k, pt[k]))
                    k = pt[k] + 1
                else:
                    k += 1
            if not children:
                leaf_bounds.append(max(0, j - i - 5))
            elif len(children) > 1:
                # An internal-loop child must enclose every mandatory child.
                # Its two endpoints consume at least one base on each wing.
                multiple_child_bounds.append(max(0, children[0][0] - i - 1 + j - children[-1][1] - 1 - 2))
        assert max(leaf_bounds) <= 30
        assert max(multiple_child_bounds, default=0) <= 30
    allowed_access = float(window[ALLOWED]["accessibility_exact"])
    no_access = float(window[NO]["accessibility_exact"])
    pkfe = float(energy[ALLOWED]["free_energy_kcal_mol"])
    nofe = float(energy[NO]["free_energy_kcal_mol"])
    context_table = "\n".join(
        f"| {c} | {float(base[ALLOWED][c]):.9g} | {float(base[NO][c]):.9g} |"
        for c in ("Stem", "Hairpin", "Bulge", "Internal", "Multiloop", "Exterior", "PseudoLoop"))
    window_table = "\n".join(
        f"| {c} | {float(window[ALLOWED][c]):.9g} | {float(window[NO][c]):.9g} |"
        for c in ("Exterior", "Hairpin", "Bulge", "Internal", "Multiloop", "PseudoLoop"))
    region_table = "\n".join(
        f"| {name} ({regions[NO, key]['start']}–{regions[NO, key]['end']}) | "
        f"{float(regions[ALLOWED, key]['paired_exact']):.6f} | {float(regions[NO, key]['paired_exact']):.6f} | "
        f"{max(('Hairpin', 'Bulge', 'Internal', 'Multiloop', 'Exterior'), key=lambda c: float(regions[NO,key][c] or 0))} |"
        for key, name in (("outer_left", "outer, left"), ("outer_right", "outer, right"),
                          ("inner_left", "inner, left"), ("inner_right", "inner, right")))
    english = (f"We compared two ensembles for rat Terc while retaining the same 61 scaffold pairs.\n"
        f"The probability that bases 123–126 were all unpaired was {allowed_access:.3f} when crossings were allowed and {no_access:.3f} when crossings were excluded.\n"
        f"Most of this accessibility came from PseudoLoop in the first ensemble and from a hairpin loop in the second.\n"
        f"This comparison shows how allowing crossings can change the predicted accessibility and structural context under a supplied scaffold.")
    readme = f"""# Terc: extensions with and without crossings of the same scaffold G

解析 A、2026-10-19。原稿・既存図は編集していない。採用版89baea52の交差を許す **exact** 出力と、同じ配列・Gを必ず保持する交差なし **exact** 分配関数を比較した。

## 主な結果

- 123–126の4塩基が**同時に**非対合となる確率：{allowed_access:.9f} → {no_access:.9f}（交差あり → 交差なし、差 {allowed_access-no_access:.9f}）。
- 交差ありでは窓のPseudoLoop成分が {float(window[ALLOWED]['PseudoLoop']):.9f}。交差なしではHairpin成分が {float(window[NO]['Hairpin']):.9f}。
- ensemble自由エネルギー：{pkfe:.12f} → {nofe:.12f} kcal/mol。差 F(交差あり)−F(交差なし) = **{pkfe-nofe:.9f} kcal/mol**。
- 参照pseudoknotの外側6対・内側3対は交差なしではすべて確率0。ただし、各塩基の代替対合と非対合contextは位置によって異なる。

## 方法と設定の対応

入力は既存の `results/local/small-rerun-89baea52-20261018/terc_data.tar::in.txt` と `rec.json`。
bpRNA_RFAM_4317、ラットTerc、{len(seq)} nt、Gはcanonicalな参照入れ子対{g.count('(')}対。
入力とexact出力のSHA-256は `results/figures/figure4_terc_89baea52_20261018/rerun_provenance.json` と全て照合した。
塩基・窓は1-based inclusive。生のRaccess出力の窓開始位置は0-basedなので122を123に変換した。

`build/linear-cparty-89baea52/LinearCParty`、`LinearPKCapR`、`LinearPKRaccess` の `--help`、sourceのCLI parser、`src/scfg2/api/exact.hh`、`api/beam.hh`、`api/README.md`、`src/scfg2/README.md` を確認した。
Gを保持したまま交差するextension対を禁止するCLI option・公開API設定は見当たらなかった。
`scfg2_restricted_partition_function_tool` の `restricted` は入力scaffoldを指し、交差禁止の設定ではない。
したがって、交差を許す側は上記tarの採用版exact出力を再利用し、交差なし側はViennaRNA Python **2.7.2** で新規計算した。
`RNA.fold_compound(seq, md)` に `hc_add_from_db(G, RNA.CONSTRAINT_DB_DEFAULT | RNA.CONSTRAINT_DB_ENFORCE_BP)` を適用した。
Gのdotは追加対を許す自由塩基であり、Gの対だけを強制する。
ENFORCE_BPを使う理由は[ViennaRNA hard-constraint documentation](https://viennarna.readthedocs.io/en/latest/grammar/constraints/hard.html)の説明に従う。
別の `hc_add_bp(i,j, ALL_LOOPS | ENFORCE)` の指定でも分配関数を独立に照合した。

| 設定 | 両条件でそろえた内容 |
|---|---|
| Parameter file | `deps/LinearCParty-89baea52/params/rna_Turner04.par`、SHA-256 `{provenance['dependencies_sha256_before']['deps/LinearCParty-89baea52/params/rna_Turner04.par']}` |
| Temperature / RT | 37 °C / {float(energy[NO]['RT_kcal_mol']):.15g} kcal/mol |
| Dangling ends / PF smoothing | dangles=2 / pf_smooth=1 |
| Canonical pairing | AU/UA、GC/CG、GU/UG。lonely pairsとGU closureを許す |
| Hairpin / ordinary internal loops | 最小hairpin 3 nt、通常のinternal/bulge探索は合計非対合数≤30 |
| Other model flags | special_hp=1、betaScale=1、logML=circ=gquad=0 |
| Salt | ViennaRNA salt=1.021 M（基準濃度、salt補正なし） |
| Beam / hairpin cap | 両側ともexact、beamなし、hairpin上限なし。有限beam用CASE3上限も使わない |

ViennaRNAではsampling用に `uniq_ML=1` を設定した（同じ構造の重複を避けるmultiloop分解）。
LinearCPartyはCPartyのpseudoknot固有エネルギー項を含む。交差なし側にはそのcontextが現れない。
ソフトウェアの実装は別であり、全設定を一つのCLI内で切り替えた比較ではない。
通常の入れ子構造の局所エネルギー規約は既存の `energy-vienna-check-89baea52-repaired-20261018` の照合を参照し、今回もViennaRNAのG制約MFE構造と8標本を採用版の公開 `ExactSession::evaluate_energy` で照合した。
全{validation['fixed_structure_energy_checks']}構造が到達可能で、単一構造PFエネルギーの最大差は {validation['fixed_structure_energy_max_abs_error']:.3g} kcal/molだった。
今回のエネルギー照合範囲はこの9構造であり、全ての入れ子構造の列挙照合は行っていない。

LinearCPartyの固定scaffoldに沿う一部のloop規則は30 ntを超えるloopも評価する。
ただし今回のGでは、直下のGの子対が一つの場合の非対合gap最大が{validation['max_single_child_unpaired_gap_in_G']} nt、Gの葉対の中に追加対を置く場合のinternal-loop上限が{max(leaf_bounds)} ntであり、いずれも30以下。
Gの子対が複数ある場合、それらをすべて囲む追加対とのinternal loopの上限も{max(multiple_child_bounds, default=0)} ntである。
Gを閉じるinternal loopでViennaRNAの30 nt制限が余分な構造を除くことはない。
追加対を閉じる通常のinternal loopは両実装で30 ntに制限される。

## accessibility と context の求め方

ZをG制約の分配関数、Z_Uをさらに123–126の全塩基を強制非対合にした分配関数とし、accessibilityを Z_U/Z とした。
位置ごとの非対合確率の積は使っていない。
分配関数のdouble精度の根 `exp_matrices.q1k[n]` と `pf_scale` からlog Zを復元した。
`pf()`の戻り値の自由エネルギーはこのbindingでfloatに丸められるため採用値には使わない。
BPPは `fc.bpp()`、非対合確率は1−BPP行和。

Exterior、Hairpin、Multiloop、Bulge/Internal合計は、対象の全塩基をそれぞれのloop contextで強制非対合にした分配関数の比からexactに求めた。
連続した非対合窓は一つのloopに属する。
BulgeとInternalの分割は[Boltzmann sampling](https://viennarna.readthedocs.io/en/latest/subopt/stochbt.html)で求めた。
一つの子対を含むloopで片側のgapが0ならBulge、両側が正ならInternalと数えた。
非対合contextの判定は、各塩基を直接囲む対の直下の子対の数とgapによる。

G制約だけのensembleから{provenance['draws_per_ensemble']:,} draws（seed {provenance['seed']}）を取り、40–160のBulge/Internal比を集計した。
123と窓123–126は、Bulge/Internalに属することも強制した別の条件付きensembleから各{provenance['draws_per_ensemble']:,} draws（seed {provenance['seed']+1}、{provenance['seed']+2}）を取り、稀な成分を分割した。
これらの標本頻度にexactなBulge/Internal合計を掛けた。
123のB/I分割の標準誤差は {float(base[NO]['BI_split_standard_error']):.4g}、窓では {float(window[NO]['BI_split_standard_error']):.4g}。
`BI_split_standard_error` は片方の成分の条件付き二項標準誤差で、両成分の誤差は逆符号に相関する。
一方の標本countが0の場合に標準誤差0となっても、成分が厳密に0であることは意味しない。
領域内でB/I標本が0なら分割成分は空欄とし、exactな合計だけを残した。
PseudoLoopは交差なしでは厳密に0。

### 窓123–126

各行は「4塩基が同時に非対合、かつそのloopに属する」joint probabilityであり、accessibilityを条件とした割合ではない。
Bulge/Internal分割以外はexact。

| Context | Crossings allowed | No crossings |
|---|---:|---:|
{window_table}
| Total | {allowed_access:.9g} | {no_access:.9g} |

### 塩基123と周辺

| Context | Crossings allowed | No crossings |
|---|---:|---:|
{context_table}

`context_region.tsv` に40–160の全塩基を保存した。
交差を許す場合は123付近のPseudoLoopが大きいが、交差なしの123ではHairpinが {float(base[NO]['Hairpin']):.6f}、対合が {float(base[NO]['Stem']):.6f} となる。

### 参照pseudoknotの9対とその塩基

外側は63–139、64–138、65–137、66–136、67–135、68–134、内側は69–132、70–131、71–130。
すべてGの少なくとも一対と交差し、交差なしでは許されない。
`reference_pairs.tsv` に各参照対のBPP、Gと交差する対の数・例、両端の総対合確率・最も確率の高い代替相手を記録した。
相手の確率が0のときは相手IDを空欄にした。

| Region | Mean paired, crossings allowed | Mean paired, no crossings | Largest unpaired context, no crossings |
|---|---:|---:|---|
{region_table}

左側の66と67は交差なしでも60と59に、それぞれ {next(float(r['left_top_partner_probability']) for r in pairs if r['ensemble']==NO and r['i']=='66'):.6f}、{next(float(r['left_top_partner_probability']) for r in pairs if r['ensemble']==NO and r['i']=='67'):.6f} の確率で対合する。
左側69–71はほぼInternal loopの非対合塩基となる。
右側の塩基には別の相手との対合や通常のloop contextが残る。
`region_summary.tsv` と `context_region.tsv` を使って、参照対が失われることと塩基が非対合になることを区別できる。

## 原稿に入れる英文案（4文、原稿には未反映）

同じ61対のscaffoldを保ち、ラットTercの二つのensembleを比較した。
123–126の全塩基が同時に非対合となる確率は、交差を許す場合{allowed_access:.3f}、交差を除く場合{no_access:.3f}だった。
このaccessibilityの大部分は、前者ではPseudoLoop、後者ではhairpin loopから生じた。
この比較は、与えたscaffoldの下で交差を許すことが、予測されるaccessibilityとstructural contextを変えうることを示す。

{english}

## Figure 3に加える小パネル案

窓123–126について、横軸を `Crossings allowed` / `No crossings`、縦軸を0–1のjoint probabilityとする2本の積み上げ棒を置く。
既存図と同じ6 contextの色を使い、棒の高さをtotal accessibilityとする。
1までの残りは白または薄灰で `At least one base paired` と示す。
棒上にtotal（{allowed_access:.3f}、{no_access:.3f}）、主な成分に `PseudoLoop` / `Hairpin` を直接書く。
短い共通ラベルで `Same sequence and 61-pair scaffold G; exact` と示す。
小パネルでは9対すべてや自由エネルギーを詰め込まず、隣の構造図の123–126に対応する印を付ける。
B/I分割を含む棒は推定内訳なので、captionに「totalとB+I合計はexact、B/Iの分割はsampling」と一言添える。
この比較だけから生物学的機能や実験的accessibilityの正確さを主張しない。

## 再現・検証・保存

```bash
python scripts/analyze_terc_nopk_89baea52.py --out results/local/terc-nopk-89baea52-20261019
python scripts/report_terc_nopk_89baea52.py results/local/terc-nopk-89baea52-20261019
```

既存出力先への計算の上書きは拒否する。再実行には新しい `--out` を使う。
一時領域は環境変数TMPDIRのnode-localパス。標本・CLI help・分配関数receipt・入力・観測用C++コード・失敗記録は `raw_outputs.tar.gz` にまとめ、共有領域へ展開しない。
小さなTSVは上記5ファイル。`provenance.json` は入力tar・parameter・binary・library・build cache・解析scriptの前後hashを持つ。
交差を許す側の実行コマンドは既存のfrozen packed-task JSONL/tarにあり、今回のarchiveはその参照と入力hashを保持する。
`receipt.jsonl` は今回の解析の完了記録、`validation.json` は数値検証。

- G全対のBPPの1からの最大差：{validation['forced_pair_max_abs_error']:.3g}。
- 二つのhard-constraint指定法のlog Z差：{validation['constraint_interfaces_logz_abs_difference']:.3g}。
- 40–160のexactなcontext合計の1からの最大差：{validation['context_mass_max_abs_error']:.3g}。
- 窓のexactなcontext合計とtotalの差：{validation['window_mass_abs_error']:.3g}。
- G違反の標本：{validation['sampled_structure_G_violations']}。loop classifierの手計算fixture：{validation['classifier_fixtures_passed']}件成功。
- 3,000標本の全塩基で独自classifierをViennaRNAの `db_to_element_string` と照合し、Stem・Hairpin・B/I合計・Multiloop・Exteriorが一致した（同APIはBulgeとInternalを区別しない）。
- 無条件標本の対合頻度とBPPの最大絶対差：{validation['sampling_max_stem_abs_error']:.6f}（40–160、参考確認）。
- 入力の幾何から、Gを閉じるinternal loopが30 ntを超えないことを追加確認した。
- 採用sourceはclean、parameter・binary・library・build cacheは計算前後で不変。原稿・調整中の外部実装の変更、git commitは行っていない。

開始時の共有inode確認：21,544,868 / 28,000,000 files。今回の永続追加は結果・小表・script・記録のみで、標本ごとのファイルを共有領域に作っていない。
試行中にこのPython bindingにない `exp_matrices.get_Z()` を呼んだ失敗はraw内に記録し、既存auditと同じ `q1k[n]` に修正した。失敗した値は使用していない。
"""
    (out / "README.md").write_text(readme)
    receipt = dict(status="success", analysis="A: Terc no-crossing comparison", source_revision=provenance["source_revision"],
                   outputs=str(out), conditional_draw_counts=counts, partition_function_calls=len(calls),
                   max_G_leaf_internal_loop_bound=max(leaf_bounds),
                   max_G_multiple_child_internal_loop_bound=max(multiple_child_bounds, default=0),
                   classifier_native_crosscheck_structures=native_context_checks,
                   scalar_validation=validation,
                   reporting_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                   normalization="Blank top_partner when its maximum probability is zero; numerical fields unchanged")
    (out / "receipt.jsonl").write_text(json.dumps(receipt, ensure_ascii=False) + "\n")
    checksums = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.iterdir())
                 if p.is_file() and p.name != "checksums.json"}
    (out / "checksums.json").write_text(json.dumps(checksums, indent=2) + "\n")
    print(english)
    print("\nValidated report:", out / "README.md")


if __name__ == "__main__":
    main()
