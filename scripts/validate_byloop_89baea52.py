#!/usr/bin/env python3
"""Measure existing by-loop oracle checks without editing the implementation.

All generated code, binaries and run logs live in node-local scratch.  The
maintained oracle's checks and inputs stay intact; additional observations
measure errors and compare finite-beam numeric/taped reverse calculations.
"""
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import tarfile
import tempfile


EXP = Path(__file__).resolve().parents[1]
OUT = EXP / "results/local/byloop-validation-89baea52-20261018"
SOURCE = EXP / "deps/LinearCParty-89baea52"
BUILD = EXP / "build/linear-cparty-89baea52"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_once(text, old, new):
    assert text.count(old) == 1, (text.count(old), old[:100])
    return text.replace(old, new)


def instrument(text):
    text = replace_once(text, '#include "scfg2/profiles/exact_profiles.hh"',
                        '#include "scfg2/profiles/exact_profiles.hh"\n'
                        '#include "scfg2/profiles/beam_profiles.hh"\n'
                        '#include <iomanip>\n#include <set>')
    text = replace_once(text, '  double weight = 0.0;',
                        '  double weight = 0.0;\n  std::vector<std::string> sources;')
    text = replace_once(text, '  const auto total_only =',
                        '  std::array<double, 6> enum_error{};\n'
                        '  std::array<std::size_t, 6> longer_window_count{};\n'
                        '  std::size_t crossing_groups = 0, boundary_windows = 0;\n'
                        '  double total_enum_error = 0.0;\n'
                        '  const auto total_only =')
    text = replace_once(text, '      initial.emitted.assign(static_cast<std::size_t>(n) + 1, kNone);',
                        '      initial.emitted.assign(static_cast<std::size_t>(n) + 1, kNone);\n'
                        '      initial.sources.resize(static_cast<std::size_t>(n) + 1);')
    text = replace_once(text, '          initial.emitted[static_cast<std::size_t>(pos)] = context_code(run.context);',
                        '          initial.emitted[static_cast<std::size_t>(pos)] = context_code(run.context);\n'
                        '          initial.sources[pos] = std::string(grammar::rule_name(deduction.rule)) +\n'
                        '              ":" + std::to_string(item.i) + ":" + std::to_string(item.j);')
    text = replace_once(text, '              combined.emitted[pos] = child.emitted[pos];',
                        '              combined.emitted[pos] = child.emitted[pos];\n'
                        '              combined.sources[pos] = child.sources[pos];')
    text = replace_once(text, '    const double p = summary.weight / z;', r'''
    const double p = summary.weight / z;
    auto complete_pairs = summary.pairs;
    std::vector<cand_pos_t> stack;
    for (cand_pos_t i = 1; i <= n; ++i) {
      if (scaffold[i-1] == '(') stack.push_back(i);
      if (scaffold[i-1] == ')') {
        engine::append_unique_traceback_pair(&complete_pairs, {stack.back(), i});
        stack.pop_back();
      }
    }
    bool crossing = false;
    for (const auto &a : complete_pairs) for (const auto &b : complete_pairs)
      crossing = crossing || (a.i < b.i && b.i < a.j && a.j < b.j);
    if (crossing) ++crossing_groups;
''')
    text = replace_once(text, '          if (w > 1 && p > 0.0) window_loop_coverage()[c] = true;', r'''
          if (w > 1 && p > 0.0) {
            window_loop_coverage()[c] = true;
            ++longer_window_count[c];
            bool crosses_emission = false;
            for (cand_pos_t t = s+1; t < s+w; ++t)
              crosses_emission = crosses_emission || summary.sources[t] != summary.sources[s];
            if (crossing && c == 5 && crosses_emission) {
              ++boundary_windows;
              if (boundary_windows == 1) {
                std::cout << "BOUNDARY {\"label\":" << std::quoted(label)
                          << ",\"start\":" << s << ",\"length\":" << w
                          << ",\"emitted\":" << std::quoted(summary.emitted.substr(1))
                          << ",\"emission_sources\":[";
                for (cand_pos_t t = s; t < s+w; ++t) {
                  if (t > s) std::cout << ',';
                  std::cout << std::quoted(summary.sources[t]);
                }
                std::cout << "]}\n";
              }
            }
          }
''')
    text = replace_once(text, '      const double want = access[w_idx][s];',
                        '      const double want = access[w_idx][s];\n'
                        '      total_enum_error = std::max(total_enum_error, std::fabs(got-want));')
    text = replace_once(text, '        const double value = profile.accessibility_by_loop[w_idx][s][c];',
                        '        const double value = profile.accessibility_by_loop[w_idx][s][c];\n'
                        '        enum_error[c] = std::max(enum_error[c], std::fabs(value-loop_access[w_idx][s][c]));')
    text = replace_once(text, '  std::size_t with_pseudoloop = 0;', r'''
  const auto emit = [&](const std::string &mode, const profiles::ProfileResult &got,
                        const profiles::ProfileResult &reference, bool enum_ref) {
    if (!got.ok || !reference.ok) {
      std::cerr << label << " " << mode << ": " << got.error << " / " << reference.error << "\n";
      std::exit(2);
    }
    std::array<double, 6> errors{};
    double sum_error = 0.0, unit_error = 0.0, min_value = 1.0, max_value = 0.0;
    for (std::size_t wi = 0; wi < windows.size(); ++wi) {
      for (cand_pos_t s = 1; s + windows[wi] - 1 <= n; ++s) {
        double sum = 0.0;
        for (std::size_t c = 0; c < 6; ++c) {
          const double value = got.accessibility_by_loop[wi][s][c];
          if (!std::isfinite(value)) std::exit(3);
          min_value = std::min(min_value, value);
          max_value = std::max(max_value, value);
          errors[c] = std::max(errors[c], std::fabs(value-reference.accessibility_by_loop[wi][s][c]));
          if (windows[wi] == 1)
            unit_error = std::max(unit_error, std::fabs(value-got.context[s][c+1]));
          sum += value;
        }
        sum_error = std::max(sum_error, std::fabs(sum-got.accessibility[wi][s]));
      }
    }
    if (enum_ref) errors = enum_error;
    if (*std::max_element(errors.begin(), errors.end()) > 1e-9 ||
        sum_error > 1e-9 || unit_error > 1e-9 || min_value < -1e-9 || max_value > 1+1e-9)
      std::exit(4);
    std::cout << std::setprecision(17) << "METRIC {\"label\":" << std::quoted(label)
              << ",\"sequence\":" << std::quoted(sequence) << ",\"scaffold\":" << std::quoted(scaffold)
              << ",\"parameter\":" << std::quoted(std::string(std::getenv("CPARTY_FIXED_ENERGY_PARAM_FILE")))
              << ",\"n\":" << n << ",\"mode\":" << std::quoted(mode)
              << ",\"groups\":" << derivations.size() << ",\"crossing_groups\":" << crossing_groups
              << ",\"pseudoloop_emission_boundary_windows\":" << boundary_windows
              << ",\"component_max_abs_error\":[";
    for (std::size_t c = 0; c < 6; ++c) std::cout << (c ? "," : "") << errors[c];
    std::cout << "],\"positive_longer_windows\":[";
    for (std::size_t c = 0; c < 6; ++c) std::cout << (c ? "," : "") << longer_window_count[c];
    std::cout << "],\"sum_max_abs_error\":" << sum_error << ",\"unit_max_abs_error\":" << unit_error
              << ",\"total_enumeration_max_abs_error\":" << total_enum_error
              << ",\"min_component\":" << min_value << ",\"max_component\":" << max_value << "}\n";
  };
  emit("exact_vs_enumeration", profile, profile, true);
  const auto keep_all = profiles::beam_profiles(sequence, scaffold,
      replay::WFinalVpWmbpOnlineWidths{}, windows, WmbClosurePolicy::RequireWeakClosure,
      beam::RankingMode::Baseline, beam::PilotExponent::TwoThirds, true);
  emit("unpruned_vs_exact", keep_all, profile, false);
  for (const auto policy : {WmbClosurePolicy::RequireWeakClosure, WmbClosurePolicy::LegacyCompatible}) {
    for (const auto rank : {beam::RankingMode::Baseline, beam::RankingMode::Outer}) {
      for (const std::size_t b : {2, 5, 20}) {
        const auto widths = replay::WFinalVpWmbpOnlineWidths::uniform(b);
        const auto numeric = profiles::beam_profiles(sequence, scaffold, widths, windows, policy,
            rank, beam::PilotExponent::TwoThirds, true);
        const auto taped = profiles::beam_profiles_taped(sequence, scaffold, widths, windows, policy,
            rank, beam::PilotExponent::TwoThirds, true);
        const std::string mode = "finite_b" + std::to_string(b) +
            (policy == WmbClosurePolicy::RequireWeakClosure ? "_strict" : "_compatible") +
            (rank == beam::RankingMode::Baseline ? "_baseline" : "_outer");
        emit(mode, numeric, taped, false);
      }
    }
  }
  std::size_t with_pseudoloop = 0;
''')
    return text


def main():
    OUT.mkdir(parents=True, exist_ok=False)
    native = SOURCE / "test/scfg2_pk_profiles_oracle_test.cc"
    provenance_files = [native, SOURCE / "test/scfg2_pk_profiles_beam_test.cc",
                        SOURCE / "test/tools/test_linear_pk_raccess.py",
                        SOURCE / "params/rna_Turner04.par", SOURCE / "params/rna_DirksPierce09.par",
                        BUILD / "libcparty_scfg2_exact.a", BUILD / "libRNA.a",
                        BUILD / "scfg2_pk_profiles_oracle_test", BUILD / "scfg2_pk_profiles_beam_test",
                        BUILD / "LinearPKRaccess", BUILD / "LinearPKCapR"]
    before = {str(p.relative_to(EXP)): sha(p) for p in provenance_files}
    records = []
    with tempfile.TemporaryDirectory(prefix="byloop-validation-", dir=os.environ.get("TMPDIR", "/tmp")) as temp:
        scratch = Path(temp)
        (scratch / "probe.cc").write_text(instrument(native.read_text()))
        commands = []

        def run(args, name, cwd=SOURCE, timeout=300):
            commands.append({"argv": list(map(str, args)), "cwd": str(cwd), "name": name})
            result = subprocess.run(list(map(str, args)), cwd=cwd, capture_output=True,
                                    text=True, timeout=timeout,
                                    env=dict(os.environ, SCFG2_BUILD_DIR=str(BUILD)))
            (scratch / (name + ".stdout")).write_text(result.stdout)
            (scratch / (name + ".stderr")).write_text(result.stderr)
            commands[-1]["returncode"] = result.returncode
            assert result.returncode == 0, (name, result.stderr[-2000:])
            return result.stdout

        try:
            run([BUILD / "scfg2_pk_profiles_oracle_test"], "maintained_oracle")
            run([BUILD / "scfg2_pk_profiles_beam_test"], "maintained_beam")
            run(["python3", SOURCE / "test/tools/test_linear_pk_raccess.py"], "maintained_cli")
            run(["/bin/c++", "-O2", "-std=gnu++17", "-I" + str(SOURCE / "src"),
                 "-I" + str(SOURCE / "test"), "-I" + str(BUILD / "generated"),
                 scratch / "probe.cc", BUILD / "libcparty_scfg2_exact.a", BUILD / "libRNA.a",
                 "-o", scratch / "probe"], "compile", timeout=300)
            stdout = run([scratch / "probe"], "measured_oracle", timeout=300)
            boundaries = []
            for line in stdout.splitlines():
                if line.startswith("METRIC "):
                    records.append(json.loads(line[7:]))
                elif line.startswith("BOUNDARY "):
                    boundaries.append(json.loads(line[9:]))
            assert len(records) == 21 * 14, len(records)
            exact = [r for r in records if r["mode"] == "exact_vs_enumeration"]
            groups = {"exact_vs_enumeration": exact,
                      "unpruned_vs_exact": [r for r in records if r["mode"] == "unpruned_vs_exact"],
                      "finite_vs_taped": [r for r in records if r["mode"].startswith("finite_")]}
            summary = {"contexts": ["Hairpin", "Bulge", "Internal", "Multi", "Exterior", "PseudoLoop"],
                       "distinct_sequence_scaffold_inputs": len({(r["sequence"], r["scaffold"]) for r in exact}),
                       "input_parameter_conditions": len(exact),
                       "length_range_nt": [min(r["n"] for r in exact), max(r["n"] for r in exact)],
                       "requested_windows": ["1", "2", "3", "5", "n", "n+3", "2 (duplicate)"],
                       "enumerated_groups_across_conditions": sum(r["groups"] for r in exact),
                       "crossing_groups_across_conditions": sum(r["crossing_groups"] for r in exact),
                       "positive_longer_window_counts": [sum(r["positive_longer_windows"][c] for r in exact) for c in range(6)],
                       "pseudoloop_emission_boundary_windows": sum(r["pseudoloop_emission_boundary_windows"] for r in exact),
                       "boundary_examples": boundaries, "comparisons": {}}
            for name, rows in groups.items():
                summary["comparisons"][name] = {
                    "conditions": len(rows),
                    "component_max_abs_error": [max(r["component_max_abs_error"][c] for r in rows) for c in range(6)],
                    "sum_max_abs_error": max(r["sum_max_abs_error"] for r in rows),
                    "unit_max_abs_error": max(r["unit_max_abs_error"] for r in rows),
                    "min_component": min(r["min_component"] for r in rows),
                    "max_component": max(r["max_component"] for r in rows)}
            assert all(c > 0 for c in summary["positive_longer_window_counts"])
            assert summary["pseudoloop_emission_boundary_windows"] > 0
            after = {str(p.relative_to(EXP)): sha(p) for p in provenance_files}
            assert before == after, "implementation or build changed during the run"
            (OUT / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
            (OUT / "case_metrics.jsonl").write_text("".join(json.dumps(r) + "\n" for r in records))
            print(json.dumps(summary, indent=2))
        finally:
            manifest = {"source_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=SOURCE, text=True).strip(),
                        "source_status": subprocess.check_output(["git", "status", "--porcelain"], cwd=SOURCE, text=True),
                        "script_sha256": sha(Path(__file__)), "file_sha256": before,
                        "compiler": subprocess.check_output(["/bin/c++", "--version"], text=True).splitlines()[0],
                        "python": platform.python_version(), "commands": commands,
                        "environment_override": {"SCFG2_BUILD_DIR": str(BUILD)},
                        "scope": "finite regression checks, not an independent grammar or energy-model validation"}
            (OUT / "provenance.json").write_text(json.dumps(manifest, indent=2) + "\n")
            with tarfile.open(OUT / "raw_outputs.tar.gz", "w:gz") as archive:
                for path in sorted(scratch.iterdir()):
                    if path.name != "probe":
                        archive.add(path, arcname=path.name)


if __name__ == "__main__":
    main()
