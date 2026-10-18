// Count every exact derivation and collect its generated pair map.
// No production profiles, pair-emission helper, or sum-product runner is used
// in this enumeration. The generator, schedule and log local weights are shared.
#include "scfg2/api/profiles.hh"
#include "scfg2/replay/w_final_exact_inside.hh"
#include <array>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_map>
using namespace cparty::scfg2;
using NT=grammar::NonTerminal; using RI=grammar::RuleId; using Key=deductions::ItemKey;
using LS=semirings::LogSumExpSemiring;
struct Group {
  std::string pairs,labels,trace;
  std::vector<int> mates; // Preserve actual partners, even for invalid crossings.
  uint64_t count=0;
  double sum=-INFINITY,best=-INFINITY,worst=INFINITY;
};
using Groups=std::map<std::string,Group>;
uint64_t add_count(uint64_t a,uint64_t b) {
  if(UINT64_MAX-a<b) throw std::runtime_error("count overflow");return a+b;
}
uint64_t mul_count(uint64_t a,uint64_t b) {
  if(b&&a>UINT64_MAX/b) throw std::runtime_error("count overflow");return a*b;
}
void merge(Groups &dest,Group x) {
  auto key=x.pairs+"|"+x.labels+"|";
  for(int partner:x.mates) key+=std::to_string(partner)+",";
  auto [it,added]=dest.emplace(key,x);
  if(!added) {
    auto &v=it->second; v.count=add_count(v.count,x.count);
    v.sum=LS::plus(v.sum,x.sum);v.best=std::max(v.best,x.best);v.worst=std::min(v.worst,x.worst);
    if(v.trace!=x.trace&&v.trace.find(" ALTERNATIVE ")==std::string::npos)
      v.trace+=" ALTERNATIVE "+x.trace;
  }
}
// Independent transcription of the context ownership table, deliberately
// avoiding profiles::unpaired_emission and all profile assembly routines.
void label_local(const deductions::Deduction &d,Group &g,const std::string &scaf) {
  auto emit=[&](int a,int b,char c) {for(int t=a;t<=b;++t) {
    if(t<1||t>int(scaf.size())) throw std::runtime_error("bad emitted position");
    auto &v=g.labels[t-1];v=v=='-'?c:'!';
  }};
  int i=d.parent.i,j=d.parent.j,k=d.split.k,l=d.split.l;
  switch(d.rule) {
    case RI::V_HAIRPIN:emit(i+1,j-1,'H');break;
    case RI::V_INTERNAL:case RI::VP_INTERNAL_LOOP:case RI::BE_INTERNAL_LOOP:
      {char c=(k>i+1&&j>l+1)?'I':'B';emit(i+1,k-1,c);emit(l+1,j-1,c);break;}
    case RI::W_EXTEND_UNPAIRED:
      if(d.child_count||scaf[j-1]=='.') emit(j,j,'E');break;
    case RI::VM_SPLIT_WMp_BASE:emit(i+1,k-1,'M');break;
    case RI::WM_START_V:case RI::WM_START_WMB:emit(i,k-1,'M');break;
    case RI::WM_EXTEND_UNPAIRED:case RI::WMv_EXTEND_UNPAIRED:case RI::WMp_EXTEND_UNPAIRED:
      emit(j,j,'M');break;
    case RI::WI_BASE_SINGLE:case RI::WI_EXTEND_UNPAIRED:case RI::WIP_EXTEND_UNPAIRED:
      emit(j,j,'P');break;
    case RI::WIP_BASEPAIR_V:case RI::WIP_BASEPAIR_WMB:case RI::VPL_SPLIT_VP:
      emit(i,k-1,'P');break;
    case RI::VPR_SPLIT_VP_BASEPAIR:emit(k+1,j,'P');break;
    case RI::BE_WIP_BASEPAIR:emit(l+1,j-1,'P');break;
    case RI::BE_BASEPAIR_WIP:emit(i+1,k-1,'P');break;
    default:break;
  }
}
Groups enumerate(std::string s,std::string g,WmbClosurePolicy policy,bool pre_wip=false) {
  weights::RuntimeSnapshot owner(s,g,weights::RuntimeOptions{2,false});
  replay::WFinalLogWeightContext ctx(owner,policy);
  engine::AllLocalWeightModel<replay::WFinalLogWeightContext> weights(&ctx);
  std::unordered_map<Key,Groups,deductions::ItemKeyHash,deductions::ItemKeyEq> chart;
  auto schedule=engine::build_pkfree_basic_exact_schedule(s.size(),ctx);
  for(const auto &key:schedule) {
    Groups result;
    for(const auto &d:deductions::generate_exact_basic_deductions_for_item(key,ctx)) {
      // Diagnostic ablation: restore exactly the four split upper bounds
      // changed by d9554a8e, without changing the implementation or its build.
      if(pre_wip && (d.rule==RI::WIP_SPLIT_V||d.rule==RI::WIP_SPLIT_WMB||
          d.rule==RI::WIP_BASEPAIR_V||d.rule==RI::WIP_BASEPAIR_WMB) &&
          d.split.k==d.parent.j-ctx.turn()-1) continue;
      if(!engine::legacy_span_schedule_supports(d)) continue;
      double w=weights.log_weight_for(d);if(!std::isfinite(w)) continue;
      Group start;start.pairs=std::string(s.size(),'.');start.labels=std::string(s.size(),'-');
      start.mates.assign(s.size(),0);
      start.count=1;start.sum=start.best=start.worst=w;
      start.trace=std::string(grammar::rule_name(d.rule))+"("+std::to_string(key.i)+","+std::to_string(key.j)+")";
      if((key.nonterminal==NT::V||key.nonterminal==NT::VP)
          &&ctx.pair_partner(key.i)!=key.j) {
        start.pairs[key.i-1]='[';start.pairs[key.j-1]=']';
        start.mates[key.i-1]=key.j;start.mates[key.j-1]=key.i;
      }
      label_local(d,start,g);
      Groups partial;merge(partial,std::move(start));
      for(std::size_t c=0;c<d.child_count&&!partial.empty();++c) {
        Key ck=d.children[c].item;
        auto it=chart.find(ck);if(it==chart.end()) {partial.clear();break;}
        Groups next;
        for(auto &[_,a]:partial) for(auto &[__,b]:it->second) {
          Group x=a;
          for(std::size_t p=0;p<s.size();++p) {
            if(b.pairs[p]!='.') {if(x.pairs[p]!='.') x.pairs[p]='!';else x.pairs[p]=b.pairs[p];}
            if(b.mates[p]) x.mates[p]=b.mates[p];
            if(b.labels[p]!='-') {if(x.labels[p]!='-') x.labels[p]='!';else x.labels[p]=b.labels[p];}
          }
          x.count=mul_count(x.count,b.count);x.sum+=b.sum;x.best+=b.best;x.worst+=b.worst;
          x.trace+="{"+b.trace+"}";merge(next,std::move(x));
        }
        partial=std::move(next);
      }
      for(auto &[_,x]:partial) merge(result,std::move(x));
    }
    chart[key]=std::move(result);
  }
  return chart[{NT::W,1,static_cast<int>(s.size()),-1,-1}];
}
std::vector<std::pair<int,int>> pairs_of(const std::string &g) {
  std::vector<int> st;std::vector<std::pair<int,int>> out;
  for(int i=0;i<int(g.size());++i) {
    if(g[i]=='['||g[i]=='(') st.push_back(i+1);
    if(g[i]==']'||g[i]==')') {if(st.empty()) throw std::runtime_error("bad pair map");out.push_back({st.back(),i+1});st.pop_back();}
  }
  if(!st.empty()||g.find('!')!=std::string::npos) throw std::runtime_error("bad pair ownership");
  return out;
}
// Enumerate all canonical, noncrossing G' disjoint from G. VP also requires
// j-i>3; the candidate superset does not impose the density-2 grammar.
std::vector<std::string> candidates(const std::string &s,const std::string &g,int lo,int hi) {
  if(lo>hi) return {""};
  auto tails=candidates(s,g,lo+1,hi);std::vector<std::string> out;
  for(auto &t:tails) out.push_back("."+t);
  if(g[lo]!='.') return out;
  for(int j=lo+4;j<=hi;++j) if(g[j]=='.') {
    std::string bp;bp+=s[lo];bp+=s[j];
    if(bp!="AU"&&bp!="UA"&&bp!="CG"&&bp!="GC"&&bp!="GU"&&bp!="UG") continue;
    auto inner=candidates(s,g,lo+1,j-1);auto right=candidates(s,g,j+1,hi);
    for(auto &a:inner) for(auto &b:right) out.push_back("["+a+"]"+b);
  }
  return out;
}
int main(int argc,char **argv) {
  if(argc<3) return 2;std::string s=argv[1],g=argv[2];
  auto policy=argc>3&&std::string(argv[3])=="legacy"?WmbClosurePolicy::LegacyCompatible:WmbClosurePolicy::RequireWeakClosure;
  std::cout<<std::setprecision(17);
  bool pre_wip=argc>4&&std::string(argv[4])=="pre-wip";
  auto root=enumerate(s,g,policy,pre_wip);int n=s.size();
  std::map<std::string,std::vector<Group>> structures;
  uint64_t count=0;double deriv_z=-INFINITY;
  for(auto &[_,x]:root) {
    // A single bracket lane cannot encode crossing G' pairs. Verify actual
    // partners before using the dot-bracket serialization as the root key.
    auto represented=pairs_of(x.pairs);std::vector<int> mates(n,0);
    for(auto [i,j]:represented) {mates[i-1]=j;mates[j-1]=i;}
    if(mates!=x.mates) throw std::runtime_error("crossing generated component or inconsistent pair map");
    for(int i=0;i<n;++i) if(mates[i]&&g[i]!='.') throw std::runtime_error("generated endpoint overlaps scaffold");
    count=add_count(count,x.count);deriv_z=LS::plus(deriv_z,x.sum);structures[x.pairs].push_back(x);
  }
  if(pre_wip) {
    for(auto &[gp,groups]:structures) {
      uint64_t c=0;for(auto &x:groups)c=add_count(c,x.count);
      std::cout<<"pre_wip_structure\t"<<gp<<'\t'<<c<<'\n';
    }
    std::cout<<"value\tstructures\t"<<structures.size()<<'\n'
             <<"value\tderivations\t"<<count<<'\n'
             <<"value\toverflow\t0\n";
    return 0;
  }
  api::ExactSession session(s,g,policy);auto z=session.inside();auto o=session.outside();
  auto pr=api::exact_structural_profiles(s,g,{1,2,3,5},policy);
  if(!z.ok||!o.ok||!pr.ok) throw std::runtime_error("public API failed");
  std::vector<bool> fixed(n+1,false);for(auto [i,j]:pairs_of(g)) fixed[i]=fixed[j]=true;
  double struct_z=-INFINITY,max_energy_error=0;uint64_t maxcount=0;
  int ambiguous=0,ownership_bad=0,context_ambiguous=0,crossing=0;
  int ownership_duplicate=0,ownership_missing=0,ownership_on_paired=0;
  std::vector<std::tuple<std::string,std::string,double>> dist;
  for(auto &[gp,groups]:structures) {
    uint64_t c=0;double sum=-INFINITY,best=-INFINITY,worst=INFINITY;
    for(auto &x:groups) {c=add_count(c,x.count);sum=LS::plus(sum,x.sum);best=std::max(best,x.best);worst=std::min(worst,x.worst);}
    maxcount=std::max(maxcount,c);ambiguous+=c>1;context_ambiguous+=groups.size()>1;
    auto ev=session.evaluate_energy(gp);if(!ev.ok||!ev.reachable) throw std::runtime_error("enumerated structure rejected by evaluator");
    double lw=-ev.energy/pr.rt_kcal_mol;struct_z=LS::plus(struct_z,lw);max_energy_error=std::max(max_energy_error,std::abs(lw-best));
    auto paired=fixed;for(auto [i,j]:pairs_of(gp)) paired[i]=paired[j]=true;
    bool pk=false;for(auto [i,j]:pairs_of(gp)) for(auto [k,l]:pairs_of(g)) pk|=(i<k&&k<j&&j<l)||(k<i&&i<l&&l<j);
    crossing+=pk;
    for(auto &x:groups) for(int i=1;i<=n;++i) {
      ownership_bad+=(x.labels[i-1]=='!'||(paired[i]!=(x.labels[i-1]=='-')));
      ownership_duplicate+=x.labels[i-1]=='!';
      ownership_missing+=!paired[i]&&x.labels[i-1]=='-';
      ownership_on_paired+=paired[i]&&x.labels[i-1]!='-';
    }
    dist.emplace_back(gp,groups[0].labels,lw);
    std::cout<<"structure\t"<<gp<<'\t'<<c<<'\t'<<sum<<'\t'<<best<<'\t'<<worst<<'\t'<<ev.energy<<'\t'<<groups[0].labels<<'\n';
    if(c>1) for(auto &x:groups) std::cout<<"witness\t"<<gp<<'\t'<<x.trace<<'\n';
  }
  auto all=candidates(s,g,0,n-1);int allowed_absent=0,rejected=0;
  for(auto &gp:all) if(!structures.count(gp)) {
    auto e=session.evaluate_energy(gp);if(!e.ok) throw std::runtime_error("candidate validation error: "+e.error);
    if(e.reachable) {++allowed_absent;std::cout<<"missing\t"<<gp<<'\n';}else ++rejected;
  }
  double bpp_err=0,qerr=0,ctxerr=0,accerr=0;
  std::array<double,7> ctx_errors{};
  std::array<double,4> acc_errors{};
  std::map<std::pair<int,int>,double> bpp;
  std::vector<double> q(n+1,0);std::vector<std::array<double,7>> cp(n+1);
  std::vector<std::vector<double>> ac(4,std::vector<double>(n+1,0));
  const std::string codes="SHBIMEP";int widths[]={1,2,3,5};
  for(auto &[gp,labels,lw]:dist) {
    double p=std::exp(lw-struct_z);auto paired=fixed;
    for(auto [i,j]:pairs_of(g)) bpp[{i,j}]+=p;
    for(auto [i,j]:pairs_of(gp)) {bpp[{i,j}]+=p;paired[i]=paired[j]=true;}
    for(int i=1;i<=n;++i) {if(!paired[i]) q[i]+=p;auto c=paired[i]?0:codes.find(labels[i-1]);if(c<7) cp[i][c]+=p;}
    for(int w=0;w<4;++w) for(int i=1;i+widths[w]-1<=n;++i) {
      bool free=true;for(int t=i;t<i+widths[w];++t) free&=!paired[t];if(free) ac[w][i]+=p;
    }
  }
  for(auto &p:o.pair_probabilities) {
    auto expected=bpp[{p.i,p.j}];
    bpp_err=std::max(bpp_err,std::abs(p.union_probability-expected));
    std::cout<<"pair_compare\t"<<p.i<<'\t'<<p.j<<'\t'<<expected<<'\t'<<p.union_probability<<'\n';
    bpp.erase({p.i,p.j});
  }
  for(auto &[ij,p]:bpp) {
    bpp_err=std::max(bpp_err,p);
    std::cout<<"pair_compare\t"<<ij.first<<'\t'<<ij.second<<'\t'<<p<<"\t0\n";
  }
  for(int i=1;i<=n;++i) {
    qerr=std::max(qerr,std::abs(q[i]-o.unpaired_probabilities[i]));
    std::cout<<"unpaired_compare\t"<<i<<'\t'<<q[i]<<'\t'<<o.unpaired_probabilities[i]<<'\n';
    for(int c=0;c<7;++c) {
      double e=std::abs(cp[i][c]-pr.contexts[i][c]);
      ctxerr=std::max(ctxerr,e);ctx_errors[c]=std::max(ctx_errors[c],e);
      std::cout<<"context_compare\t"<<i<<'\t'<<codes[c]<<'\t'<<cp[i][c]<<'\t'<<pr.contexts[i][c]<<'\n';
    }
    for(int w=0;w<4;++w) if(i+widths[w]-1<=n) {
      double e=std::abs(ac[w][i]-pr.accessibility[w][i]);
      accerr=std::max(accerr,e);acc_errors[w]=std::max(acc_errors[w],e);
      std::cout<<"access_compare\t"<<widths[w]<<'\t'<<i<<'\t'<<ac[w][i]<<'\t'<<pr.accessibility[w][i]<<'\n';
    }
  }
  auto value=[&](const char *k,auto v){std::cout<<"value\t"<<k<<'\t'<<v<<'\n';};
  value("n",n);value("structures",structures.size());value("derivations",count);value("max_count",maxcount);
  value("ambiguous",ambiguous);value("context_ambiguous",context_ambiguous);value("ownership_bad",ownership_bad);
  value("ownership_duplicate",ownership_duplicate);value("ownership_missing",ownership_missing);value("ownership_on_paired",ownership_on_paired);
  value("overflow",0);value("rt_kcal_mol",pr.rt_kcal_mol);
  value("crossing_structures",crossing);value("candidate_superset",all.size());value("candidate_rejected",rejected);value("allowed_absent",allowed_absent);
  value("exact_logZ",z.log_partition_weight);value("derivation_logZ",deriv_z);value("structure_logZ",struct_z);
  value("logZ_error",std::abs(struct_z-z.log_partition_weight));value("grouping_logZ_error",std::abs(deriv_z-z.log_partition_weight));
  value("energy_error",max_energy_error);value("bpp_error",bpp_err);value("unpaired_error",qerr);value("context_error",ctxerr);value("access_error",accerr);
  for(int c=0;c<7;++c) value(("context_error_"+codes.substr(c,1)).c_str(),ctx_errors[c]);
  for(int w=0;w<4;++w) value(("access_error_"+std::to_string(widths[w])).c_str(),acc_errors[w]);
  return 0;
}
