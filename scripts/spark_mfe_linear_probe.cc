// Read-only observer for the pinned LinearCParty libraries; limited to small inputs.
#include "scfg2/api/exact.hh"
#include "scfg2/api/detail.hh"
#include "scfg2/replay/w_final_exact_inside.hh"
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace sc = cparty::scfg2;
void emit(const sc::engine::DeductionTraceNode<double> &node, double rt, int depth) {
  double local = node.chart_value;
  for (const auto &child : node.children) local -= child.chart_value;
  const auto &k = node.item;
  std::cout << "loop\t" << depth << '\t' << sc::grammar::nonterminal_name(k.nonterminal)
            << '\t' << k.i << '\t' << k.j << '\t' << k.ip << '\t' << k.jp << '\t'
            << (node.chosen_deduction ? sc::grammar::rule_name(node.chosen_deduction->rule) : "base")
            << '\t' << -rt * local << '\t' << -rt * node.chart_value << '\n';
  for (const auto &child : node.children) emit(child, rt, depth + 1);
}
int main(int argc, char **argv) {
  if (argc != 5) return 2;
  try {
    std::string s = argv[1], g = argv[2], target = argv[3];
    if (s.size() > 24) throw std::runtime_error("observer limited to 24 nt");
    auto policy = std::string(argv[4]) == "legacy" ? sc::WmbClosurePolicy::LegacyCompatible
                                                : sc::WmbClosurePolicy::RequireWeakClosure;
    sc::api::ExactSession session(s, g, policy);
    auto mfe = session.viterbi();
    auto eval = session.evaluate_energy(target);
    std::cout << std::setprecision(17);
    std::cout << "mfe\t" << mfe.ok << '\t' << mfe.reachable << '\t' << mfe.energy
              << '\t' << mfe.generated_structure << '\t' << mfe.union_structure << '\n';
    std::cout << "evaluation\t" << eval.ok << '\t' << eval.reachable << '\t' << eval.energy
              << '\t' << eval.generated_structure << '\t' << eval.union_structure << '\n';
    if (!eval.ok) std::cout << "error\t" << eval.error << '\n';
    if (!eval.ok || !eval.reachable) return 0;
    W_final_pf owner(s, g, true, false, false, 2, 0., 0, false, false, false, false, true);
    const double rt = -owner.partition_free_energy_from_log_weight(1.);
    const auto &parameters = owner.log_pf_parameters();
    std::cout << "parameters\t" << parameters.temperature << '\t' << parameters.model_details.dangles
              << '\t' << parameters.model_details.pf_smooth << '\t' << -rt * parameters.expMLclosing
              << '\t' << -rt * parameters.expMLintern[1] << '\t' << -rt * parameters.expMLbase << '\n';
    auto trace = sc::replay::run_w_final_exact_output_constrained_viterbi_traceback<true>(
        owner, sc::api::detail::to_traceback_pairs(eval.generated_pairs), policy);
    const double traced_energy = owner.partition_free_energy_from_log_weight(trace.run.chart.get(trace.root));
    if (std::abs(traced_energy - eval.energy) > 1e-10)
      throw std::runtime_error("observer disagrees with evaluate_energy");
    std::cout << "rt\t" << rt << '\n';
    emit(trace.trace, rt, 0);
    return 0;
  } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
