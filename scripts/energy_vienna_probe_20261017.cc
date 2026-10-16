// Read-only observer of the pinned SCFG2 runtime. Build only in local scratch.
#include "scfg2/api/detail.hh"
#include "scfg2/engine/outside_marginals.hh"
#include "scfg2/replay/w_final_exact_inside.hh"
#include "scfg2/replay/w_final_exact_outside.hh"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

using namespace cparty::scfg2;

static double trace(const engine::DeductionTraceNode<double>& node,
                    const replay::WFinalLogWeightContext& ctx, double rt,
                    std::size_t target) {
  double total = 0.;
  if (node.chosen_deduction) {
    const auto& d = *node.chosen_deduction;
    const double energy = -rt *
        engine::AllLocalWeightModel<replay::WFinalLogWeightContext>(&ctx).log_weight_for(d);
    std::cout << "RULE\t" << target << '\t' << grammar::rule_name(d.rule)
              << '\t' << d.parent.i << '\t' << d.parent.j
              << '\t' << d.parent.ip << '\t' << d.parent.jp
              << '\t' << d.split.k << '\t' << d.split.l << '\t' << energy << '\n';
    total += energy;
  }
  for (const auto& child : node.children) total += trace(child, ctx, rt, target);
  return total;
}

int main(int argc, char** argv) {
  if (argc != 7) {
    std::cerr << "usage: probe PARAM SEQUENCE SCAFFOLD DANGLES ensemble|eval TRACE\n";
    return 2;
  }
  try {
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", argv[1], 1);
    const std::string s = argv[2], g = argv[3], mode = argv[5];
    const bool emit_trace = std::atoi(argv[6]) != 0;
    weights::RuntimeSnapshot owner(s, g, weights::RuntimeOptions{std::atoi(argv[4]), false});
    replay::WFinalLogWeightContext ctx(owner);
    const double rt = owner.log_pf_parameters().kT / 1000.;
    const auto& p = owner.log_pf_parameters();
    std::cout << std::setprecision(17);
    std::cout << "META\t" << rt << '\t' << p.temperature << '\t'
              << p.model_details.dangles << '\t' << p.model_details.pf_smooth
              << '\t' << ctx.turn() << '\t' << MAXLOOP << '\n';
    std::cout << "MLPARAM\t" << -rt*p.expMLclosing << '\t' << -rt*p.expMLintern[1]
              << '\t' << -rt*p.expMLbase << '\t' << -rt*p.expTermAU << '\n';
    if (mode == "ensemble") {
      const auto n = owner.length();
      const auto schedule = engine::build_pkfree_basic_exact_schedule(n, ctx);
      using S = semirings::LogSumExpSemiring;
      const auto chart = replay::run_w_final_exact_dp_schedule<S>(schedule, ctx,
          [](const auto&) { return true; }, [](const auto&...) {});
      const deductions::ItemKey root{grammar::NonTerminal::W, 1, n, -1, -1};
      const double z = chart.get(root);
      std::cout << "Z\t" << z << '\t' << -rt*z << '\n';
      const auto beta = replay::run_w_final_exact_regenerated_outside(schedule, ctx, chart, root);
      const auto marginal = engine::collect_outside_marginals_by_item(chart,
          [&](const auto& item, double alpha) -> std::optional<double> {
            if (S::is_zero(alpha)) return std::nullopt;
            return std::exp(alpha + beta.get(item) - z);
          }, n, ctx);
      if (!marginal.ok) throw std::runtime_error(marginal.error);
      for (const auto& pair : marginal.pairs)
        std::cout << "BPP\t" << pair.i << '\t' << pair.j << '\t'
                  << pair.union_probability << '\n';
      const auto best = replay::run_w_final_exact_viterbi_traceback<true>(owner);
      const auto db = api::detail::render_generated_structure(n, best.pairs, ctx);
      std::cout << "MFE\t" << -rt*best.run.chart.get(best.root) << '\t' << db << '\n';
    } else if (mode != "eval") {
      throw std::runtime_error("unknown mode");
    }
    std::string target;
    std::size_t index = 0;
    while (std::getline(std::cin, target)) {
      if (target.empty()) continue;
      std::vector<api::BasePair> pairs;
      std::string error;
      if (target.size() != s.size() ||
          !api::detail::parse_dot_bracket_pairs(target, &pairs, &error))
        throw std::runtime_error("invalid target: " + error);
      const auto wanted = api::detail::to_traceback_pairs(pairs);
      const auto evaluated = replay::run_w_final_exact_output_constrained_viterbi_traceback<true>(owner, wanted);
      const double value = evaluated.run.chart.get(evaluated.root);
      const bool same = engine::same_traceback_pair_set(evaluated.generated_pairs, wanted);
      std::cout << "EVAL\t" << index << '\t' << same << '\t'
                << std::isfinite(value) << '\t' << -rt*value << '\n';
      if (same && std::isfinite(value) && emit_trace) {
        const double total = trace(evaluated.trace, ctx, rt, index);
        if (std::abs(total + rt*value) > 1e-9)
          throw std::runtime_error("trace energies do not sum to target energy");
      }
      ++index;
    }
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
