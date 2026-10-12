// Paper-level finite validation; links the pinned library without changing it.
// Reuse the existing structure/count oracle, never its production sampler.
#define main historical_enumeration_main
#include CPARTY_VALIDATION_ORACLE
#undef main
#include "scfg2/api/sampling.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"
#include "scfg2/semirings/factorized_logsumexp.hh"
#include <fstream>
#include <functional>
#include <optional>

struct BeamMass { double log_mass=-INFINITY; uint64_t count=0; };
using BeamDist=std::map<std::string,BeamMass>; // byte partner table, one-based
void add_mass(BeamDist &out,const std::string &key,BeamMass mass) {
  auto &v=out[key]; v.log_mass=LS::plus(v.log_mass,mass.log_mass);
  v.count=add_count(v.count,mass.count);
}
std::string combine_mates(std::string a,const std::string &b) {
  for(size_t i=0;i<a.size();++i) if(b[i]) {
    if(a[i]) throw std::runtime_error("repeated generated endpoint in retained derivation");
    a[i]=b[i];
  }
  return a;
}
std::string extension(const std::string &mates) {
  std::string out(mates.size(),'.');
  std::vector<int> stack;
  for(size_t i=0;i<mates.size();++i) if(mates[i]) {
    int j=static_cast<unsigned char>(mates[i]);
    if(j<1||j>int(mates.size())||static_cast<unsigned char>(mates[j-1])!=i+1)
      throw std::runtime_error("invalid retained partner table");
    if(j>int(i+1)) {out[i]='[';stack.push_back(j);}
    else {out[i]=']';if(stack.empty()||stack.back()!=int(i+1))
      throw std::runtime_error("crossing extension in retained graph");stack.pop_back();}
  }
  return out;
}
void emit_bpp(const api::OutsideResult &out) {
  if(!out.ok||!out.reachable) throw std::runtime_error("outside failed: "+out.error);
  for(const auto &p:out.pair_probabilities)
    std::cout<<"bpp\t"<<p.i<<'\t'<<p.j<<'\t'<<p.union_probability<<'\n';
}
void exact_reference(std::string s,std::string g) {
  auto groups=enumerate(s,g,WmbClosurePolicy::RequireWeakClosure);
  api::ExactSession session(s,g);auto z=session.inside();
  if(!z.ok||!z.reachable) throw std::runtime_error("exact inside failed");
  W_final_pf owner(s,g,true,false,false,2,0.,0,false,false,false,false,true);
  const double rt=-owner.partition_free_energy_from_log_weight(1.0);
  std::map<std::string,BeamMass> dist;
  for(const auto &[_,v]:groups) add_mass(dist,v.pairs,{v.sum,v.count});
  double total=-INFINITY;
  for(const auto &[gp,v]:dist) {
    auto energy=session.evaluate_energy(gp);
    if(!energy.ok||!energy.reachable) throw std::runtime_error("oracle structure rejected");
    const double lw=-energy.energy/rt;
    // The structural probability uses evaluated energy only after count=1.
    if(v.count!=1||std::abs(v.log_mass-lw)>1e-9)
      throw std::runtime_error("ambiguous exact structure or energy disagreement");
    total=LS::plus(total,lw);
    std::cout<<"structure\t"<<gp<<'\t'<<v.count<<'\t'<<lw<<'\t'<<energy.energy<<'\n';
  }
  std::cout<<"value\tapi_logZ\t"<<z.log_partition_weight<<'\n'
           <<"value\tenumerated_logZ\t"<<total<<'\n'
           <<"value\trt_kcal_mol\t"<<rt<<'\n';
  emit_bpp(session.outside());
  for(auto policy:{WmbClosurePolicy::RequireWeakClosure,WmbClosurePolicy::LegacyCompatible}) {
    api::ExactSession mfe_session(s,g,policy);auto v=mfe_session.viterbi();
    if(!v.ok||!v.reachable) throw std::runtime_error("MFE preparation failed");
    std::cout<<"mfe\t"<<(policy==WmbClosurePolicy::RequireWeakClosure?"default":"legacy")
             <<'\t'<<v.energy<<'\t'<<v.generated_structure<<'\t'<<v.union_structure<<'\n';
  }
}
void beam_reference(std::string s,std::string g,api::BeamConfig cfg) {
  W_final_pf owner(s,g,true,false,false,2,0.,0,false,false,false,false,true);
  replay::WFinalLogWeightContext ctx(owner,cfg.wmb_closure_policy);
  semirings::FactorizedInsideGraph graph;
  using SR=semirings::FactorizedLogSumExpSemiring;
  auto chart=[&] {
    SR::Scope scope(graph);
    auto widths=replay::WFinalVpWmbpOnlineWidths::uniform(cfg.beam_size);
    widths.max_hairpin_unpaired=cfg.max_hairpin_unpaired;
    widths.max_case3_ancestor_depth=cfg.max_case3_ancestor_depth;
    return replay::run_w_final_vp_wmbp_online<SR>(ctx,s.size(),widths).chart;
  }();
  auto root=chart.get({NT::W,1,static_cast<int>(s.size()),-1,-1});
  if(!std::isfinite(root.log_value)) throw std::runtime_error("empty retained ensemble");
  // Attach semantic pair ownership from chart items, independently of sampler
  // helpers. Several aliases may own the same pair; conflicting owners fail.
  std::map<uint32_t,std::pair<int,int>> owners;
  for(const auto &[k,v]:chart.values()) {
    if(k.nonterminal!=NT::V&&k.nonterminal!=NT::VP&&k.nonterminal!=NT::VP_DIRECT&&k.nonterminal!=NT::VP_CLOSED) continue;
    if(ctx.pair_partner(k.i)==k.j||v.node>=graph.states().size()) continue;
    auto [it,inserted]=owners.emplace(v.node,std::make_pair(k.i,k.j));
    if(!inserted&&it->second!=std::make_pair(k.i,k.j))
      throw std::runtime_error("conflicting retained pair owners");
  }
  std::vector<std::optional<BeamDist>> memo(graph.states().size());
  size_t visited=0;
  std::function<const BeamDist &(uint32_t)> expand=[&](uint32_t id)->const BeamDist & {
    if(memo[id]) return *memo[id];
    BeamDist result;const auto &node=graph.states()[id];
    for(size_t t=node.first_term;t<size_t(node.first_term)+node.term_count;++t) {
      const auto &term=graph.terms()[t];
      double local=term.log_value;
      for(size_t f=term.first_factor;f<size_t(term.first_factor)+term.factor_count;++f) {
        auto child=graph.factors()[f];if(child>=id) throw std::runtime_error("non-topological retained graph");
        local-=graph.states()[child].log_value;
      }
      BeamDist partial;partial.emplace(std::string(s.size(),0),BeamMass{local,1});
      for(size_t f=term.first_factor;f<size_t(term.first_factor)+term.factor_count;++f) {
        const auto &child=expand(graph.factors()[f]);BeamDist next;
        for(const auto &[a,av]:partial) for(const auto &[b,bv]:child)
          add_mass(next,combine_mates(a,b),{av.log_mass+bv.log_mass,mul_count(av.count,bv.count)});
        partial=std::move(next);
      }
      for(const auto &[key,mass]:partial) {
        auto annotated=key;auto it=owners.find(id);
        if(it!=owners.end()) {
          auto [i,j]=it->second;
          if(annotated[i-1]||annotated[j-1]||g[i-1]!='.'||g[j-1]!='.')
            throw std::runtime_error("invalid retained pair emission");
          annotated[i-1]=static_cast<char>(j);annotated[j-1]=static_cast<char>(i);
        }
        add_mass(result,annotated,mass);
      }
    }
    memo[id]=std::move(result);++visited;return *memo[id];
  };
  double total=-INFINITY;
  for(const auto &[mates,mass]:expand(root.node)) {
    total=LS::plus(total,mass.log_mass);
    std::cout<<"structure\t"<<extension(mates)<<'\t'<<mass.count<<'\t'<<mass.log_mass<<"\tnan\n";
  }
  std::cout<<"value\tapi_logZ\t"<<root.log_value<<'\n'
           <<"value\tenumerated_logZ\t"<<total<<'\n'
           <<"value\tgraph_states\t"<<graph.states().size()<<'\n'
           <<"value\tvisited_states\t"<<visited<<'\n';
  auto out=api::beam_outside(s,g,cfg);
  if(!out.ok||!out.reachable) throw std::runtime_error("beam outside failed: "+out.error);
  for(const auto &p:out.pair_probabilities)
    std::cout<<"bpp\t"<<p.i<<'\t'<<p.j<<'\t'<<p.union_probability<<'\n';
  std::cout<<"value\toutside_logZ\t"<<out.log_partition_weight<<'\n';
}
int main(int argc,char **argv) {
  if(argc!=8) return 2;
  try {
    const std::string s=argv[1],g=argv[2],mode=argv[3];
    if(s.size()>24) throw std::runtime_error("probe is limited to 24 nt");
    const size_t width=std::stoull(argv[4]),count=std::stoull(argv[5]);
    const uint64_t seed=std::stoull(argv[6]);
    api::BeamConfig cfg;cfg.beam_size=width;cfg.max_hairpin_unpaired=s.size();cfg.max_case3_ancestor_depth=s.size();
    std::cout<<std::setprecision(17);
    if(mode=="exact") exact_reference(s,g);else if(mode=="beam") beam_reference(s,g,cfg);
    else throw std::runtime_error("unknown mode");
    auto draws=mode=="exact"?api::exact_sample(s,g,count,seed):api::beam_sample(s,g,count,seed,cfg);
    if(!draws.ok) throw std::runtime_error("sampling failed: "+draws.error);
    std::cout<<"value\tsampling_logZ\t"<<draws.log_partition_weight<<'\n';
    std::ofstream raw(argv[7]);
    for(const auto &pairs:draws.generated_pairs) {
      std::string mates(s.size(),0);
      for(const auto &p:pairs) {
        if(mates[p.i-1]||mates[p.j-1]||g[p.i-1]!='.'||g[p.j-1]!='.')
          throw std::runtime_error("invalid sample endpoint");
        mates[p.i-1]=static_cast<char>(p.j);mates[p.j-1]=static_cast<char>(p.i);
      }
      raw<<extension(mates)<<'\n';
    }
    if(!raw) throw std::runtime_error("raw draws write failed");
    return 0;
  } catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
}
