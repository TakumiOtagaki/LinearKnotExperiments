// Observe the short-span domain independently of output-constrained evaluation.
#include "scfg2/replay/w_final_exact_inside.hh"
#include <iomanip>
#include <iostream>

int main() {
  using namespace cparty::scfg2;
  std::string s="CACGUACAAGUAUCAAAGU", g=".............(...).";
  W_final_pf owner(s,g,true,false,false,2,0.,0,false,false,false,false,true);
  replay::WFinalLogWeightContext ctx(owner);
  const auto chart=replay::run_w_final_exact_dp_schedule<semirings::LogMaxSemiring>(
      engine::build_pkfree_basic_exact_schedule(s.size(),ctx),ctx,
      [](const auto &){return true;},[](const auto &,const auto &,const auto &,const auto &,auto){});
  std::cout<<std::setprecision(17);
  for(auto nt:{grammar::NonTerminal::VP,grammar::NonTerminal::WMBP}) {
    double value=chart.get({nt,13,17,-1,-1});
    std::cout<<grammar::nonterminal_name(nt)<<"\t13\t17\t"<<value<<'\t'
             <<owner.partition_free_energy_from_log_weight(value)<<'\n';
  }
  const auto active=deductions::exact_round_pseudoloop_region_active(ctx,13,17);
  std::cout<<"partial_WMBP_domain_active\t"<<active<<'\n';
  return active ? 1 : 0;
}
