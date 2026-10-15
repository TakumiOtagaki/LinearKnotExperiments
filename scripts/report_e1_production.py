#!/usr/bin/env python3
"""Write the Japanese production E1 note and a scalar comparison with earlier E1.

Read aggregate TSV/JSON only. All historical comparisons retain their reference
and source, and match RNA IDs before comparing scalar metrics.
"""
import argparse
from collections import defaultdict
import csv
import json
import math
from pathlib import Path
import re

import numpy as np

from analyze_bprna_short import CORE, statistics
from run_cli import digest

ROOT = Path(__file__).resolve().parents[1]
PRIMARY = [CORE[0], "base_state_tv_mean_Gfree", "context_tv_mean_Gfree"]


def rows(path):
    with path.open() as stream:
        yield from csv.DictReader(stream, delimiter="\t")


def scalar(row, key):
    return np.nan if row[key] == "NA" else float(row[key])


def validate_row(row):
    """Check identities across independently computed per-RNA output fields."""
    assert row["analysis_status"] == "success", row["input"]
    delta = scalar(row, CORE[0])
    assert abs(delta - (scalar(row,"free_energy_kcal_mol") - scalar(row,"reference_free_energy_kcal_mol"))) < 1e-10
    assert abs(scalar(row,"delta_log_z") + delta/scalar(row,"rt_kcal_mol")) < 1e-9
    assert math.isclose(scalar(row,"partition_function_retention"), math.exp(-delta/scalar(row,"rt_kcal_mol")), rel_tol=1e-12, abs_tol=1e-300)
    if scalar(row,"free_bases") > 0:
        for metric in PRIMARY[1:]:
            assert -1e-9 <= scalar(row,metric) <= 1+1e-9
        component_sum = sum(scalar(row,"context_"+c+"_mae_Gfree") for c in
                            ("Stem","Hairpin","Bulge","Internal","Multiloop","Exterior","PseudoLoop"))
        assert abs(component_sum/2 - scalar(row,"context_tv_mean_Gfree")) < 1e-12
        assert abs(scalar(row,"unpaired_mae_Gfree") - scalar(row,"access_w1_probability_mae_eligible")) < 1e-9
    for w in (1,4,8,16):
        prefix = f"access_w{w}_"
        eligible = scalar(row,prefix+"eligible_windows")
        assert scalar(row,prefix+"all_windows") == scalar(row,"length")-w+1
        assert 0 <= eligible <= scalar(row,prefix+"all_windows")
        assert eligible == sum(scalar(row,prefix+k) for k in ("both_positive","both_zero",
            "reference_positive_candidate_zero","reference_zero_candidate_positive"))
        for metric in ("probability_mae_eligible","probability_max_error_eligible"):
            value = scalar(row,prefix+metric)
            assert np.isnan(value) == (eligible == 0)
            assert np.isnan(value) or -1e-9 <= value <= 1+1e-9
        opening = scalar(row,prefix+"opening_mae_both_positive_kcal_mol")
        assert np.isnan(opening) == (scalar(row,prefix+"both_positive") == 0)
        assert np.isnan(opening) or math.isfinite(opening) and opening >= 0


def paired(old, new):
    old, new = np.asarray(old, dtype=float), np.asarray(new, dtype=float)
    valid = np.isfinite(old) & np.isfinite(new)
    difference = new[valid] - old[valid]
    return {"old": statistics(old), "new": statistics(new), "new_minus_old": statistics(difference),
            "smaller": int((difference < -1e-8).sum()), "larger": int((difference > 1e-8).sum()),
            "equal_within_1e_minus8": int((abs(difference) <= 1e-8).sum())}


def comparisons(short, long, previous_short, previous_long, ranking_short, ranking_long):
    for previous in (previous_short, previous_long):
        saved = json.loads((previous / "provenance.json").read_text())
        if digest(previous / "metrics.tsv") != saved["metrics_sha256"]:
            raise ValueError("previous metrics checksum mismatch: " + str(previous))
    new = {}
    for directory in (short, long):
        provenance = json.loads((directory / "provenance.json").read_text())
        assert digest(directory / "metrics.tsv") == provenance["metrics_sha256"]
        count = 0
        for row in rows(directory / "metrics.tsv"):
            validate_row(row)
            key = row["input"], int(row["beam"])
            assert key not in new, key
            new[key] = {k: row[k] for k in ["tier", "free_energy_kcal_mol"] + PRIMARY}
            count += 1
        assert count == provenance["metric_rows"]
        (directory / "validation.json").write_text(json.dumps({"all_passed": True, "rows": count,
            "metrics_sha256": provenance["metrics_sha256"], "validator_sha256": digest(Path(__file__)),
            "checks": ["unique RNA/beam", "successful rows", "signed F difference and Z ratio",
                       "TV bounds and context-component identity", "w1 marginal identity",
                       "window-count partition", "empty-domain NA", "finite opening energies"]}, indent=2) + "\n")
    sources = [previous_short / "metrics.tsv", previous_short / "provenance.json",
               previous_long / "metrics.tsv", previous_long / "provenance.json", ranking_short, ranking_long]
    result = {"sources": [{"path": str(p), "sha256": digest(p)} for p in sources],
              "scope": "matched RNA IDs; earlier source versions differ; no attribution to ranking alone",
              "previous_short": [], "previous_long": [], "ranking_short": [], "ranking_long": []}
    groups = defaultdict(list)
    old_caps = {}
    for old in rows(previous_short / "metrics.tsv"):
        beam = int(old["beam"])
        if beam == 200:
            old_caps[old["input"],old["cap"]] = scalar(old,CORE[0])
        # b20 capb was an alias, and therefore has no independently measured row.
        if old["cap"] != ("default" if beam == 20 else "capb"): continue
        current = new[old["input"], beam]
        groups[old["tier"], beam].append(({k: old[k] for k in PRIMARY}, current))
    for (tier, beam), pairs in sorted(groups.items()):
        result["previous_short"].append({"bin": tier, "beam": beam, "rna": len(pairs),
            "reference": "each source version's exact", "metrics": {k: paired(
                [scalar(o,k) for o,n in pairs], [scalar(n,k) for o,n in pairs]) for k in PRIMARY}})
    cap_names = sorted(n for n,c in old_caps if c == "capb")
    result["previous_cap_control_b200"] = {
        "fixed30_vs_max_b30": paired([old_caps[n,"default"] for n in cap_names],
                                     [old_caps[n,"capb"] for n in cap_names]),
        "max_abs_max_b30_minus_nocap": max(abs(old_caps[n,"capb"]-old_caps[n,"nocap"]) for n in cap_names)}
    groups.clear()
    for old in rows(previous_long / "metrics.tsv"):
        beam = int(old["beam"])
        current = new[old["input"], beam]
        ref = new[old["input"], 200]
        groups[old["tier"], beam].append((float(old["delta_F_vs_b200_kcal_mol"]),
            float(current["free_energy_kcal_mol"]) - float(ref["free_energy_kcal_mol"])))
    for (tier, beam), pairs in sorted(groups.items()):
        result["previous_long"].append({"bin": tier, "beam": beam, "rna": len(pairs),
            "reference": "b200 of each source version", "delta_F": paired(*zip(*pairs))})
    metric_by_program = {"LinearCParty": "base_state_tv_mean_Gfree", "LinearPKCapR": "context_tv_mean_Gfree",
                         "LinearPKRaccess": CORE[0]}
    for label, path in (("ranking_short", ranking_short), ("ranking_long", ranking_long)):
        groups.clear()
        for old in json.loads(path.read_text())["observations"]:
            match = re.fullmatch(r"beam(\d+)_capb_(outer|baseline)", old["mode"])
            if not match or old["analysis_status"] != "success": continue
            beam = int(match[1]); current = new[old["input"], beam]
            metric = metric_by_program[old["program"]]
            groups[current["tier"], beam, match[2], metric].append((old["metrics"][metric], scalar(current,metric)))
        for (tier, beam, ranking, metric), pairs in sorted(groups.items()):
            result[label].append({"bin": tier, "beam": beam, "ranking": ranking, "metric": metric,
                "rna": len(pairs), "reference": "exact" if label.endswith("short") else "b400 outer; not exact",
                **paired(*zip(*pairs))})
    return result


def fmt(value):
    return "NA" if value is None else f"{value:.4g}"


def quantiles(group, metric):
    return " / ".join(fmt(group["metrics"][metric][q]) for q in ("median", "p90", "p99", "max"))


def write_note(short, long, destination, comparison):
    summaries = [json.loads((p / "summary.json").read_text()) for p in (short,long)]
    provenances = [json.loads((p / "provenance.json").read_text()) for p in (short,long)]
    def link(path, label=None):
        import os
        return f"[{label or path.name}]({os.path.relpath(path, destination.parent)})"
    lines = ["# E1 本番：d21618c5 の beam と d9804df6 の exact（2026-10-16）", "",
        "400 nt 未満の24,625配列を exact と比較し、400 nt 以上の1,438配列を b=400 と比較した。",
        "順位付けは outer、局所上限は h=d=max(b,30)、エネルギーパラメータは Turner 2004、温度は37°Cである。",
        "保存済みの同じ配列 S と scaffold G に対する3 CLI の出力を使った。", "",
        "## 入力と参照", "", "| 範囲 | 配列数 | 幅 b | 参照 |", "| --- | ---: | --- | --- |"]
    for summary, provenance in zip(summaries,provenances):
        for tier, lo, hi in provenance["bins"]:
            count = next(g["rna"] for g in summary["groups"] if g["bin"] == tier)
            lines.append(f"| {lo}–{hi} nt | {count:,} | {', '.join(map(str,provenance['beams']))} | {summary['reference_mode']} |")
    lines += ["", "短い区分の beam は `e1-bprna-short-d21618c5-20261014`、exact は全件 `e1-bprna-short-exact-d9804df6-20261015` から読んだ。",
        "beam run の古い exact 73,875件は、成功したものも含めて参照から除いた。",
        f"この除外には最終 checkpoint で記録された{provenances[0]['excluded_job_status_counts'].get('failed',0)}件の失敗が含まれ、依頼文と10月15日の記録の25件とは異なる。",
        "28件は14配列の CapR・Raccess 各14件の `std::bad_alloc` で、今回の参照には使っていない。",
        "d9804df6 の outside の変更については、14配列で beam b=20/200 の出力一致と exact の最大数値差2.2×10⁻¹³が確認済みである。",
        link(ROOT / "results/local/exact-d9804df6-vs-d21618c5-20261015/summary.txt", "版間の照合記録") + "を provenance に hash とともに記録した。",
        "長い区分は `e1-bprna-long-d21618c5-20261014` の全6幅を使い、最大幅 `beam400_outer` を参照とした。",
        "長い区分の値は有限幅同士の差であり、exact に対する誤差ではない。", "",
        "## 指標と評価できない値の扱い", "",
        "既存の `analyze_accuracy.py` の定義をそのまま使い、まず RNA ごとに計算してから RNA を等重みで集計した。",
        "ΔF=F_b−F_ref は符号付きの差で、Raccess の全桁の metadata を使い、負値も切り捨てずに残した。",
        "F=−RT log Z、RT=0.6163207755 kcal/mol とし、JSON と TSV に Z_b/Z_ref と長さで割った ΔF も残した。",
        "long の `partition_function_retention` 列も Z_b/Z_400 を表し、exact の質量保持率を意味しない。",
        "BPP TV は、G の自由塩基ごとに相手塩基と非対合を合わせた分布の TV を求めた平均である。",
        "loop-profile TV は、同じ自由塩基の7成分の分布の TV を求めた平均である。",
        "accessibility は、G の対合塩基を一つも含まない window で、確率差の絶対値の平均と最大を窓長1/4/8/16ごとに求めた。",
        "opening energy の差 RT|log p_b−log p_ref| は両確率が正の window でだけ平均し、参照で正・有限幅で0の window、逆向き、両方0を別に数えた。",
        "評価できる window がない RNA と、両側正の window がない RNA は NA とし、0で埋めなかった。",
        "分位点は線形補間で、以下の4数は中央値 / p90 / p99 / 最大を表す。",
        "全指標・全幅・全長さ区分の分位点と NA の RNA 数は `quantiles.tsv`、window 数は `window_counts.tsv` にある。", ""]
    for label, summary, width in zip(("400 nt 未満：exact に対する誤差", "400 nt 以上：b=400 との差"), summaries, (200,300)):
        lines += ["## " + label, "", "ΔF の単位は kcal/mol、TV は無次元である。", "",
            "| nt | b | ΔF：中央値 / p90 / p99 / 最大 | BPP TV p90 | loop TV p90 | ΔF>0.01 の RNA |",
            "| --- | ---: | --- | ---: | ---: | ---: |"]
        for g in summary["groups"]:
            if g["bin"] == "all": continue
            tier = g["bin"][2:].replace("_", "–")
            count = g["delta_F_thresholds"]["0.01"]["positive"]["count"]
            lines.append(f"| {tier} | {g['beam']} | {quantiles(g,CORE[0])} | {fmt(g['metrics'][PRIMARY[1]]['p90'])} | {fmt(g['metrics'][PRIMARY[2]]['p90'])} | {count}/{g['rna']} |")
        lines += ["", f"以下は b={width} の周辺確率の分位点である。", "",
            "| nt | BPP TV：中央値 / p90 / p99 / 最大 | loop TV：中央値 / p90 / p99 / 最大 |",
            "| --- | --- | --- |"]
        for g in summary["groups"]:
            if g["bin"] != "all" and g["beam"] == width:
                lines.append(f"| {g['bin'][2:].replace('_','–')} | {quantiles(g,PRIMARY[1])} | {quantiles(g,PRIMARY[2])} |")
        lines += ["", f"以下は b={width} の accessibility と opening energy の分位点である。",
            "NA は評価域のない RNA 数を確率 / opening の順に示し、失われた window 数は参照で正・有限幅で0の数である。", "",
            "| nt | w | 確率 MAE：中央値 / p90 / p99 / 最大 | opening MAE kcal/mol：中央値 / p90 / p99 / 最大 | NA RNA | 失われた window |",
            "| --- | ---: | --- | --- | ---: | ---: |"]
        for g in summary["groups"]:
            if g["bin"] == "all" or g["beam"] != width: continue
            for w in (1,4,8,16):
                p = f"access_w{w}_probability_mae_eligible"; e = f"access_w{w}_opening_mae_both_positive_kcal_mol"
                lost = g["window_counts"][f"access_w{w}_reference_positive_candidate_zero"]["total_windows"]
                lines.append(f"| {g['bin'][2:].replace('_','–')} | {w} | {quantiles(g,p)} | {quantiles(g,e)} | {g['metrics'][p]['missing_or_out_of_domain_rna']} / {g['metrics'][e]['missing_or_out_of_domain_rna']} | {lost} |")
        lines.append("")
    lines += ["## b=300 と b=400 の残る差", ""]
    for g in summaries[1]["groups"]:
        if g["beam"] != 300 or g["bin"] == "all": continue
        t = g["delta_F_thresholds"]["1"]["absolute"]["count"]
        negative = g["delta_F_thresholds"]["0.01"]["negative"]["count"]
        lines.append(f"{g['bin'][2:].replace('_','–')} nt の F300−F400 は中央値{fmt(g['metrics'][CORE[0]]['median'])}、p99 {fmt(g['metrics'][CORE[0]]['p99'])}、最大{fmt(g['metrics'][CORE[0]]['max'])}、最小{fmt(g['metrics'][CORE[0]]['min'])} kcal/mol で、絶対差1 kcal/mol超は{t}/{g['rna']}配列、−0.01 kcal/mol未満は{negative}配列だった。")
    lines += ["b=400 の自己比較は評価域のある指標で0となった。",
        "b=400 より広い幅との比較と、長い区分の exact との比較は行っていない。", "", "## 前回との比較", ""]
    def previous_short(beam):
        return next(g for g in comparison["previous_short"] if g["bin"] == "nt300_399" and g["beam"] == beam)["metrics"][CORE[0]]
    a, b = previous_short(100), previous_short(200)
    rank = next(g for g in comparison["ranking_short"] if g["bin"] == "nt300_399" and g["beam"] == 100 and g["ranking"] == "outer" and g["metric"] == CORE[0])
    cap = comparison["previous_cap_control_b200"]
    rank_long = [g for g in comparison["ranking_long"] if g["beam"] == 300 and g["ranking"] == "outer" and g["metric"] == CORE[0]]
    rank_long_difference = max(max(abs(g["new_minus_old"]["min"]),abs(g["new_minus_old"]["max"])) for g in rank_long)
    old_long = sorted((g for g in comparison["previous_long"] if g["beam"] == 100),
                      key=lambda g: int(g["bin"][2:].split("_")[0]))
    lines += [f"[前回の短い区分](e1_bprna_short_capb_20261012.md)と同じ max(b,30) 条件では、300–399 nt の b=100 の ΔF p90 は{fmt(a['old']['p90'])}から{fmt(a['new']['p90'])}、最大は{fmt(a['old']['max'])}から{fmt(a['new']['max'])} kcal/mol となり、b=200 の最大は{fmt(b['old']['max'])}から{fmt(b['new']['max'])} kcal/mol となった。",
        f"outer を使った[短い順位付け評価](../../experiments/results/local/e1-ranking-eval-short-20261014-accuracy.json)の同じ300–399 nt・b=100では、ΔF の最大は{fmt(rank['old']['max'])}、今回との差の絶対最大は{fmt(max(abs(rank['new_minus_old']['min']),abs(rank['new_minus_old']['max'])))} kcal/mol だった。",
        "max(b,30) は前回の採用条件から引き継いでおり、今回の版間比較だけでは上限と順位付けそれぞれの効果を分けられない。",
        f"前回の同じ15aa7742での b=200 の上限対照では、固定30に比べ max(b,30) の ΔF が小さい配列は{cap['fixed30_vs_max_b30']['smaller']:,}、大きい配列は{cap['fixed30_vs_max_b30']['larger']:,}で、max(b,30) と上限なしの ΔF 差の絶対最大は{fmt(cap['max_abs_max_b30_minus_nocap'])} kcal/mol だった。",
        "[前回の長い区分](e1_bprna_long_capb_20261013.md)の b=200 参照にそろえて F100−F200 を計算すると、4区分の中央値は " +
        "、".join(f"{fmt(g['delta_F']['old']['median'])}→{fmt(g['delta_F']['new']['median'])}" for g in old_long) + " kcal/mol だった。",
        "長い区分の BPP・loop profile の前回表は b=200 参照、今回表は b=400 参照なので、表の値の差をそのまま改善量とはしない。",
        f"16配列の[長い順位付け評価](../../experiments/results/local/e1-ranking-eval-long-20261014-vs-b400outer.json)についても同じ ID と b=400 参照で照合し、b=300・outer の ΔF の今回との差の絶対最大は{fmt(rank_long_difference)} kcal/mol だった。",
        "ΔF・BPP TV・loop TV の分位点と配列内差は `previous_comparison.json` に保存した。", "",
        "## 到達不能・欠測", ""]
    for label, p in zip(("short", "long"),provenances):
        lines.append(f"{label} の選択した CLI は{p['selected_job_status_counts'].get('success',0):,}件成功、指標行は{p['metric_rows']:,}行で、CLI 間の整合チェックは{p['consistency_checks'].get('success',0):,}組すべて通った。")
    lines += ["選択したジョブの失敗、到達不能な root、欠けた出力座標、読み取り失敗は0件だった。",
        "評価域が空である NA と確率0の window は、この0件に含めず別に数えた。", "",
        "次は全体の件数で、各欄の4数は w=1 / 4 / 8 / 16 の順である。",
        "長さ区分別の全件数は `window_counts.tsv` と `quantiles.tsv` に保存した。", "",
        "| 範囲 | b | 確率 NA RNA | opening NA RNA | 参照正・有限幅0 window | 参照0・有限幅正 window |",
        "| --- | ---: | --- | --- | --- | --- |"]
    for label, summary in zip(("short", "long"),summaries):
        for g in summary["groups"]:
            if g["bin"] != "all": continue
            probability_na, opening_na, lost, gained = [], [], [], []
            for w in (1,4,8,16):
                probability_na.append(g["metrics"][f"access_w{w}_probability_mae_eligible"]["missing_or_out_of_domain_rna"])
                opening_na.append(g["metrics"][f"access_w{w}_opening_mae_both_positive_kcal_mol"]["missing_or_out_of_domain_rna"])
                lost.append(g["window_counts"][f"access_w{w}_reference_positive_candidate_zero"]["total_windows"])
                gained.append(g["window_counts"][f"access_w{w}_reference_zero_candidate_positive"]["total_windows"])
            counts = [" / ".join(map(str,v)) for v in (probability_na,opening_na,lost,gained)]
            lines.append(f"| {label} | {g['beam']} | " + " | ".join(counts) + " |")
    lines += ["", "## 開発上の確認", "",
        "array plan と task の source・model・binary・入力 hash、変更検知フラグ、実際のコマンド・result・資源ログ・使った出力の hash を照合した。",
        "実装の一般的な正しさの証明と、これらの出力保存・集計の検証は区別する。",
        "追加した小規模テスト4件は、別 task の exact との照合、古い exact の除外、符号、0イベント、b=400 参照、表の再生成と空の評価域を確認した。",
        "最初の全テスト97件では、既存 stage-timing の RSS 単調性確認が1件失敗し、96件が成功した。",
        "この既存テストが実際の RNA CLI を呼ぶことを確認した後は除外し、追加テストを含む残る94件が成功した。",
        "集計本体は RNA CLI を呼ばず、scheduler にも投入していない。", "",
        "## 出力と再現", "",
        "tar は展開せず member の stream から読み、別に分割された exact task の照合には /tmp の一時 SQLite 1個を使った。",
        "共有領域には RNA ごとの入力・中間ファイルを作っていない。",
        "時間と RSS は各幅の分布だけを `resource_reference.json` に保存し、比較には使わない。", ""]
    for label, path in (("short",short),("long",long)):
        lines += [f"- {label}：" + "、".join(link(path / filename,filename) for filename in
            ("metrics.tsv","quantiles.tsv","window_counts.tsv","summary.json","provenance.json","validation.json","resource_reference.json")) + "。",
            f"- {label} 図：" + "、".join(link(path / filename,filename) for filename in
            ("core_quantiles.png","accessibility_quantiles.png")) + "。"]
    lines += ["- 前回との照合：" + link(short / "previous_comparison.json") + "。", "",
        "論文 repo の root から次を実行する。",
        "新規集計では未作成の出力先を使い、TSV から表と図だけを再生成するときは `--summary-only` を加える。", "", "```bash",
        "export PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/e1-production-matplotlib",
        "python3 experiments/scripts/analyze_e1_production.py short \\",
        "  --beam-run experiments/results/runs/e1-bprna-short-d21618c5-20261014 \\",
        "  --exact-run experiments/results/runs/e1-bprna-short-exact-d9804df6-20261015 \\",
        "  --out experiments/results/local/e1-bprna-short-d21618c5-20261016",
        "python3 experiments/scripts/analyze_e1_production.py long \\",
        "  --run experiments/results/runs/e1-bprna-long-d21618c5-20261014 \\",
        "  --reference beam400_outer \\",
        "  --out experiments/results/local/e1-bprna-long-d21618c5-20261016",
        "python3 experiments/scripts/report_e1_production.py",
        "python3 - <<'PY'",
        "from pathlib import Path",
        "import unittest",
        "suite, loader = unittest.TestSuite(), unittest.TestLoader()",
        "for path in sorted(Path('experiments/tests').glob('test_*.py')):",
        "    if path.name != 'test_stage_timing.py':",
        "        suite.addTests(loader.discover('experiments/tests', pattern=path.name))",
        "result = unittest.TextTestRunner().run(suite)",
        "raise SystemExit(not result.wasSuccessful())",
        "PY", "```", "", "## 限界", "",
        "対象は保存済みの bpRNA 配列、scaffold、モデル、測定した幅に限る。",
        "有限幅では beam と局所上限を同時に用いており、exact との誤差を beam pruning だけに帰さない。",
        "b=300 と b=400 の差が小さい配列でも、exact との差はこの集計からは決まらない。",
        "平均 TV は位置ごとの周辺分布の平均であり、完成構造全体の分布の TV を測っていない。",
        "opening energy の分位点は両側正の評価域に限り、参照で正・有限幅で0の window は含まれない。",
        "E1 は1回の実行と複数 worker を含むため、時間・RSS の比較や計算量の主張には用いない。",
        "この記録には生物学的な結論や、測定範囲を越える収束の結論を加えていない。", ""]
    destination.write_text("\n".join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--short", type=Path, default=ROOT / "results/local/e1-bprna-short-d21618c5-20261016")
    parser.add_argument("--long", type=Path, default=ROOT / "results/local/e1-bprna-long-d21618c5-20261016")
    parser.add_argument("--note", type=Path, default=ROOT.parent / "textbook/notes/e1_bprna_d21618c5_20261016.md")
    args = parser.parse_args()
    local = ROOT / "results/local"
    comparison = comparisons(args.short, args.long, local / "e1-bprna-short-capb-20261012",
        local / "e1-bprna-long-capb-20261012", local / "e1-ranking-eval-short-20261014-accuracy.json",
        local / "e1-ranking-eval-long-20261014-vs-b400outer.json")
    (args.short / "previous_comparison.json").write_text(json.dumps(comparison, indent=2, allow_nan=False) + "\n")
    write_note(args.short, args.long, args.note, comparison)


if __name__ == "__main__":
    main()
