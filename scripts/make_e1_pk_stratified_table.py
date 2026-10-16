#!/usr/bin/env python3
"""Render one SI table directly from the audited E1 pseudoknot summaries."""

import argparse
import csv
import json
from pathlib import Path

from run_cli import digest


def load_tsv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def number(value):
    if value == "NA":
        return "---"
    return format(float(value), ".3g")


def render(directory):
    provenance = json.loads((directory / "provenance.json").read_text())
    if provenance["validation"] != "passed":
        raise ValueError("analysis validation did not pass")
    for name, checksum in provenance["output_sha256"].items():
        if digest(directory / name) != checksum:
            raise ValueError("summary changed: " + name)
    rna = {(r["axis"], r["stratum"], int(r["beam"]), r["metric"]): r
           for r in load_tsv(directory / "rna_error_quantiles.tsv")}
    local = {(r["axis"], r["stratum"], int(r["beam"]), r["metric"]): r
             for r in load_tsv(directory / "local_error_quantiles.tsv")}
    counts = {(r["axis"], r["stratum"]): r for r in load_tsv(directory / "stratum_counts.tsv")}
    lines = [
        r"% exact の交差対の期待数で RNA を層別し、高確率の位置・対・窓の局所誤差も示す。",
        r"% A は RNA 内平均から求めた RNA 間 p90 と最大値、件数は全 RNA / 16 nt 適格 RNA。",
        r"% B は exact 確率 0.1 以上の事象の絶対誤差の中央値・p90・最大値と RNA 内最大値の RNA 間 p90。",
        r"% 窓は total accessibility であり、重なる窓を含む。",
        r"\begin{table}[!htbp]",
        r"\centering",
        r"\fontsize{9}{11}\selectfont",
        r"\caption{Exact-reference errors in RNAs and positions with pseudoknot probability mass.",
        r"A: crossing mass is the sum of exact probabilities of pairs that cross at least one pair in $G$; each pair is counted once.",
        r"$\Delta F=F_b-F_{\mathrm{exact}}$ is in kcal/mol; TV is total variation distance and MAE is mean absolute error.",
        r"The RNA counts give all RNAs / RNAs with eligible 16-nt windows.",
        r"Each error cell gives p90 / maximum across RNA-level values, using the definitions in Table~\ref{tab:si_accuracy_short}.",
        r"The groups are not matched by RNA length.",
        r"B: exact probability is at least 0.1 for the selected PseudoLoop bases, crossing pairs, or total-accessibility windows.",
        r"The counts give RNAs with selected events / selected bases, pairs, or windows.",
        r"The median, p90 and maximum summarize absolute errors pooled over those events; the last column gives the p90 of the maximum error within each selected RNA.",
        r"The 4-nt windows can overlap.",
        r"Accessibility in both panels is total accessibility; \texttt{e} denotes powers of ten.}",
        r"\label{tab:si_pk_stratified}",
        r"\setlength{\tabcolsep}{2.5pt}",
        r"\begin{tabular}{@{}lrrrrrr@{}}",
        r"\toprule",
        r"\multicolumn{7}{l}{A. RNA-level errors by exact crossing mass} \\",
        r"\midrule",
        r"Crossing mass & RNAs & $b$ & \shortstack[r]{$\Delta F$\\p90 / max} & \shortstack[r]{BPP TV\\p90 / max} & \shortstack[r]{Context TV\\p90 / max} & \shortstack[r]{16-nt MAE\\p90 / max} \\",
        r"\midrule",
    ]
    metrics = ("delta_free_energy_kcal_mol", "base_state_tv_mean_Gfree",
               "context_tv_mean_Gfree", "access_w16_probability_mae_eligible")
    for layer, label in (("lt0.1", r"$<0.1$"), ("0.1_to_lt1", r"$0.1$--$<1$"), ("ge1", r"$\geq1$")):
        count = counts["exact_crossing_pair_sum", layer]["rna"]
        n16 = rna["exact_crossing_pair_sum", layer, 20, metrics[-1]]["n"]
        for beam in (20, 50, 100, 200):
            values = [" / ".join(number(rna["exact_crossing_pair_sum", layer, beam, m][q])
                                for q in ("p90", "max")) for m in metrics]
            first = [label, f"{int(count):,} / {int(n16):,}"] if beam == 20 else ["", ""]
            lines.append(" & ".join(first + [str(beam)] + values) + r" \\")
        if layer != "ge1":
            lines.append(r"\addlinespace")
    lines.extend([
        r"\midrule",
        r"\multicolumn{7}{l}{B. Local absolute probability errors at exact probability $\geq0.1$} \\",
        r"\midrule",
        r"Selected event & RNAs / events & $b$ & Median & p90 & Maximum & \shortstack[r]{p90 RNA\\maximum} \\",
        r"\midrule",
    ])
    for metric, label in (("pseudoloop_ge0.1", "PseudoLoop base"), ("crossing_bpp_ge0.1", "Crossing pair"),
                          ("access_w4_ge0.1", "4-nt window")):
        for beam in (20, 50, 100, 200):
            row = local["all", "all", beam, metric]
            maximum = rna["all", "all", beam, metric + "_max"]
            first = [label, f'{int(row["rna_selected"]):,} / {int(row["n"]):,}'] if beam == 20 else ["", ""]
            lines.append(" & ".join(first + [str(beam)] + [number(row[k]) for k in ("p50", "p90", "max")]
                                    + [number(maximum["p90"])]) + r" \\")
        if metric != "access_w4_ge0.1":
            lines.append(r"\addlinespace")
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    table = render(args.directory)
    if args.out:
        args.out.write_text(table)
    else:
        print(table, end="")


if __name__ == "__main__":
    main()
