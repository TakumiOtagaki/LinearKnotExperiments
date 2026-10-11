// Historical BPP diagnostic only; not built or used in E0 PF parity reports.
// Linked against an observation copy of the seven-patch oracle. The only
// edits replace scalar arithmetic and mark V/VP pair emissions for reverse AD.
#include "h_globals.hh"
#include "part_func.hh"
#include "W_final.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>

int main(int argc, char **argv) {
  try {
    if (argc != 3 || std::string(argv[1]) != "-i") throw std::runtime_error("usage: e0-legacy -i record");
    std::ifstream input(argv[2]); std::string name, seq, g;
    std::getline(input,name); std::getline(input,seq); std::getline(input,g);
    const char *parameter=std::getenv("CPARTY_FIXED_ENERGY_PARAM_FILE");
    if (!parameter || !vrna_params_load(parameter, VRNA_PARAMETER_FORMAT_DEFAULT)) throw std::runtime_error("parameter load failed");
    sparse_tree tree(g,seq.size()); SHAPEData shape("",seq.size());
    W_final mfe(seq,g,shape,false,false,2);
    const double best=mfe.hfold(tree);
    // Zero MFE scaling seed matches the implementation's independent oracle probe.
    W_final_pf pf(seq,g,shape,false,false,2,0.0,0,false,false,1.0);
    oracle_ad::Real energy=pf.hfold_pf(tree);
    const double rt=pf.exp_params_->kT/1000.0;
    if (!std::isfinite(double(energy))) throw std::runtime_error("nonfinite legacy energy");
    oracle_ad::reverse(energy);
    std::map<std::pair<int,int>,double> pairs;
    for (const auto &t:oracle_ad::tags) pairs[{t.i,t.j}] += -oracle_ad::tape[t.node].adj*t.value/rt;
    for (int i=1;i<=int(seq.size());++i)
      if (tree.tree[i].pair > i) pairs[{i,tree.tree[i].pair}]=1;
    std::cout<<std::setprecision(17)<<"{\"policy\":\"legacy-seven-patch-oracle\",\"log_z\":"<<-double(energy)/rt
      <<",\"free_energy\":"<<double(energy)<<",\"rt\":"<<rt<<",\"mfe\":"<<best
      <<",\"bpp_method\":\"reverse derivative of oracle V/VP log multipliers\",\"pairs\":[";
    bool first=true;
    for (const auto &[ij,p]:pairs) { if(p==0)continue; if(!first)std::cout<<',';first=false;
      std::cout<<'['<<ij.first<<','<<ij.second<<','<<p<<']'; }
    std::cout<<"],\"ad_nodes\":"<<oracle_ad::tape.size()<<"}\n";
    return 0;
  } catch (const std::exception &e) { std::cerr<<e.what()<<'\n';return 1; }
}
