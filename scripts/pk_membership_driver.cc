// Observation-only consumer of the pinned LinearCParty exact implementation.
// Compile and run via analyze_pk_membership.py; never modify the implementation.
#include "scfg2/api/detail.hh"
#include "scfg2/api/exact.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace s = cparty::scfg2;
using Trace = s::engine::ViterbiTracebackResult<double>;

std::string quote(const std::string &v) {
  std::ostringstream o;
  o << '"';
  for (unsigned char c : v) {
    if (c == '"' || c == '\\') o << '\\' << c;
    else if (c < 32) o << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
    else o << c;
  }
  return o.str() + '"';
}

double trace_energy(const Trace &trace, s::weights::RuntimeSnapshot &owner) {
  s::replay::WFinalLogWeightContext ctx(owner);
  s::engine::AllLocalWeightModel<decltype(ctx)> weight(&ctx);
  double logw = 0;
  std::vector<const s::engine::DeductionTraceNode<double> *> pending{&trace.trace};
  while (!pending.empty()) {
    const auto *node = pending.back();
    pending.pop_back();
    if (node->chosen_deduction) logw += weight.log_weight_for(*node->chosen_deduction);
    for (const auto &child : node->children) pending.push_back(&child);
  }
  return owner.partition_free_energy_from_log_weight(logw);
}

// Extract linear parameter contributions by rescoring the SAME selected trace.
// No folding, target selection or thermodynamic table is changed by this audit.
std::map<std::string, double> terms(const Trace &trace, const std::string &seq,
                                    const std::string &g, double energy) {
  using Model = s::weights::PseudoknotEnergies;
  auto &defaults = s::weights::pseudoknot_energy_defaults();
  const Model saved = defaults;
  struct Restore { Model &ref; Model saved; ~Restore() { ref = saved; } } restore{defaults, saved};
  const std::vector<std::pair<std::string, int Model::*>> fields{
      {"exterior", &Model::exterior}, {"in_multiloop", &Model::in_multiloop},
      {"in_pseudoloop", &Model::in_pseudoloop}, {"band", &Model::band},
      {"unpaired", &Model::unpaired}, {"nested_region", &Model::nested_region},
      {"multiloop_closing", &Model::multiloop_closing},
      {"multiloop_branch", &Model::multiloop_branch},
      {"multiloop_unpaired", &Model::multiloop_unpaired},
      {"band_multiloop_closing", &Model::band_multiloop_closing},
      {"band_multiloop_branch", &Model::band_multiloop_branch},
      {"band_multiloop_unpaired", &Model::band_multiloop_unpaired}};
  std::map<std::string, double> out;
  double accounted = 0;
  for (const auto &[name, field] : fields) {
    defaults = saved;
    defaults.*field += 100; // +1 kcal/mol for each occurrence
    s::weights::RuntimeSnapshot altered(seq, g);
    const double count = trace_energy(trace, altered) - energy;
    if (std::abs(count - std::round(count)) > 1e-7)
      throw std::runtime_error("noninteger model coefficient: " + name);
    out[name + "_count"] = std::round(count);
    out[name] = std::round(count) * (saved.*field) / 100.;
    accounted += out[name];
  }
  for (const auto &[name, field] : std::vector<std::pair<std::string, double Model::*>>{
           {"band_stack", &Model::stack_multiplier}, {"band_internal", &Model::internal_multiplier}}) {
    defaults = saved;
    defaults.*field = 0;
    s::weights::RuntimeSnapshot altered(seq, g);
    out[name] = energy - trace_energy(trace, altered);
    accounted += out[name];
  }
  out["other_thermodynamic"] = energy - accounted;
  return out;
}

void emit_trace(const Trace &trace, s::weights::RuntimeSnapshot &owner,
                const std::string &seq, const std::string &g, bool breakdown) {
  const double logw = trace.run.chart.get(trace.root);
  const bool reachable = std::isfinite(logw);
  std::cout << "{\"reachable\":" << (reachable ? "true" : "false");
  if (reachable) {
    const double energy = owner.partition_free_energy_from_log_weight(logw);
    const double audited = trace_energy(trace, owner);
    if (std::abs(energy - audited) > 1e-8)
      throw std::runtime_error("trace local-weight sum disagrees with root");
    s::replay::WFinalLogWeightContext ctx(owner);
    const auto gp = s::api::detail::render_generated_structure(seq.size(), trace.pairs, ctx);
    std::cout << ",\"energy\":" << energy << ",\"log_weight\":" << logw
              << ",\"trace_sum_error\":" << audited - energy
              << ",\"generated\":" << quote(gp) << ",\"pairs\":[";
    bool first = true;
    for (const auto &pair : trace.generated_pairs) {
      if (!first) std::cout << ',';
      first = false;
      std::cout << '[' << pair.i << ',' << pair.j << ']';
    }
    std::cout << ']';
    if (breakdown) {
      std::cout << ",\"terms\":{";
      first = true;
      for (const auto &[name, value] : terms(trace, seq, g, energy)) {
        if (!first) std::cout << ',';
        first = false;
        std::cout << quote(name) << ':' << value;
      }
      std::cout << '}';
    }
  }
  std::cout << '}';
}

int main(int argc, char **argv) {
  if (argc != 4 && argc != 5) { std::cerr << "usage: DRIVER SEQUENCE SCAFFOLD PARAM_FILE [--energy-only]\n"; return 2; }
  try {
    const std::string seq = argv[1], g = argv[2];
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", argv[3], 1);
    std::cout << std::setprecision(17);
    s::api::ExactSession session(seq, g);
    const bool energy_only = argc == 5 && std::string(argv[4]) == "--energy-only";
    if (argc == 5 && !energy_only) throw std::runtime_error("unknown observer option");
    std::cout << "{\"log_z\":";
    bool first = true;
    if (energy_only) {
      std::cout << "null,\"probabilities\":[";
    } else {
      const auto outside = session.outside();
      if (!outside.ok || !outside.reachable) throw std::runtime_error("outside: " + outside.error);
      std::cout << outside.log_partition_weight << ",\"probabilities\":[";
      for (const auto &p : outside.pair_probabilities) {
        if (!first) std::cout << ',';
        first = false;
        std::cout << '[' << p.i << ',' << p.j << ',' << p.union_probability << ']';
      }
    }
    std::cout << "],\"mfe\":" << std::flush;
    {
      s::weights::RuntimeSnapshot owner(seq, g);
      const auto trace = s::replay::run_w_final_exact_viterbi_traceback<true>(owner);
      emit_trace(trace, owner, seq, g, true);
      std::cout << std::flush;
      const auto check = session.evaluate_energy(s::api::detail::to_base_pairs(trace.generated_pairs));
      const double e = owner.partition_free_energy_from_log_weight(trace.run.chart.get(trace.root));
      if (!check.ok || !check.reachable || std::abs(check.energy - e) > 1e-8)
        throw std::runtime_error("public evaluator failed to recover MFE target");
    }
    std::cout << ",\"evaluations\":{";
    first = true;
    std::string line;
    while (std::getline(std::cin, line)) {
      if (line.empty()) continue;
      const auto tab = line.find('\t');
      if (tab == std::string::npos) throw std::runtime_error("missing target tab");
      const auto label = line.substr(0, tab), target = line.substr(tab + 1);
      std::vector<s::api::BasePair> pairs;
      std::string error;
      if (!s::api::detail::parse_dot_bracket_pairs(target, &pairs, &error)) throw std::runtime_error(error);
      s::weights::RuntimeSnapshot owner(seq, g);
      const auto trace = s::replay::run_w_final_exact_output_constrained_viterbi_traceback<true>(
          owner, s::api::detail::to_traceback_pairs(pairs));
      if (!first) std::cout << ',';
      first = false;
      std::cout << quote(label) << ':';
      emit_trace(trace, owner, seq, g, label == "full" || label == "reported" || label == "filtered" || label == "scaffold");
      if (!std::isfinite(trace.run.chart.get(trace.root))) {
        const auto legacy = s::replay::run_w_final_exact_output_constrained_viterbi_traceback<true>(
            owner, s::api::detail::to_traceback_pairs(pairs), s::WmbClosurePolicy::LegacyCompatible);
        std::cout << ',' << quote(label + "_legacy") << ':';
        // Local factors are unchanged by the completed-WMB domain policy.
        emit_trace(legacy, owner, seq, g, false);
      }
    }
    std::cout << "},\"rt_kcal\":";
    s::weights::RuntimeSnapshot owner(seq, g);
    std::cout << owner.parameters().kT / 1000. << "}\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
