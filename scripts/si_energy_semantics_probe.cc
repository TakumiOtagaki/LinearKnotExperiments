// Small read-only observer for the pinned SCFG2 library; all build files go to /tmp.
#include "scfg2/api/exact.hh"
#include "scfg2/replay/w_final_exact_inside.hh"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
using namespace cparty::scfg2;

static double print_trace(const engine::DeductionTraceNode<double>& node,
                         const replay::WFinalLogWeightContext& ctx, double rt) {
  double sum = 0.;
  if (node.chosen_deduction) {
    const auto& d = *node.chosen_deduction;
    const double energy = -rt * engine::AllLocalWeightModel<replay::WFinalLogWeightContext>(&ctx).log_weight_for(d);
    std::cout << "RULE\t" << grammar::rule_name(d.rule) << '\t'
              << d.parent.i << '\t' << d.parent.j << '\t'
              << d.parent.ip << '\t' << d.parent.jp << '\t'
              << d.split.k << '\t' << d.split.l << '\t' << energy << '\n';
    sum += energy;
  }
  for (const auto& child : node.children) sum += print_trace(child, ctx, rt);
  return sum;
}

int main(int argc, char** argv) {
  if (argc != 6) return 2;
  std::string s = argv[1], g = argv[2], gp = argv[3];
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", argv[4], 1);
  const int dangles = std::atoi(argv[5]);
  W_final_pf pf(s, g, true, false, false, dangles, 0., 0,
                false, false, false, false, true);
  replay::WFinalLogWeightContext ctx(pf);
  const double rt = pf.log_pf_parameters().kT / 1000.;
  std::cout << std::setprecision(17);
  const auto& p = pf.log_pf_parameters();
  const auto& pk = pf.log_pk_parameters();
  std::cout << "PARAM\t" << rt << '\t' << -rt*p.expMLclosing << '\t'
            << -rt*p.expMLintern[1] << '\t' << -rt*p.expMLbase << '\t'
            << -rt*pk.a << '\t' << -rt*pk.b << '\t' << -rt*pk.c << '\n';
  std::vector<engine::TracebackPair> target;
  std::vector<int> stack;
  for (int i = 1; i <= static_cast<int>(gp.size()); ++i) {
    if (gp[i-1] == '(') stack.push_back(i);
    if (gp[i-1] == ')') {
      target.push_back({stack.back(), i, engine::TracebackPairLane::Square});
      stack.pop_back();
    }
  }
  for (const auto& text : {g, gp}) {
    stack.clear();
    for (int i = 1; i <= static_cast<int>(text.size()); ++i) {
      if (text[i-1] == '(') stack.push_back(i);
      if (text[i-1] == ')') {
        const int left = stack.back(); stack.pop_back();
        std::cout << "HAIRPIN\t" << left << '\t' << i << '\t'
                  << -rt*ctx.hairpin_weight(left, i) << '\n';
      }
    }
  }
  auto result = replay::run_w_final_exact_output_constrained_viterbi_traceback<true>(pf, target);
  const double root = result.run.chart.get(result.root);
  const bool exact = engine::same_traceback_pair_set(result.generated_pairs, target);
  std::cout << "ROOT\t" << exact << '\t' << std::isfinite(root) << '\t' << -rt*root << '\n';
  if (std::isfinite(root) && exact) {
    const double sum = print_trace(result.trace, ctx, rt);
    std::cout << "SUM\t" << sum << '\n';
    if (std::abs(sum + rt*root) > 1e-10) return 4;
  }
  const auto api_result = api::evaluate_energy(s, g, gp);
  std::cout << "API_D2\t" << api_result.ok << '\t' << api_result.reachable
            << '\t' << api_result.energy << '\n';
  return 0;
}
