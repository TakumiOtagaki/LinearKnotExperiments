// Inside-only observation linked to the manifest-verified oracle objects.
#include "h_globals.hh"
#include "part_func.hh"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
int main(int argc,char **argv) {
  try {
  if(argc!=3||std::string(argv[1])!="-i")throw std::runtime_error("usage: e0-legacy-scalar -i record");
  std::ifstream input(argv[2]);std::string name,seq,g;
  std::getline(input,name);std::getline(input,seq);std::getline(input,g);
  if(!input||seq.empty()||seq.size()!=g.size())throw std::runtime_error("invalid input record");
  const char *parameter=std::getenv("CPARTY_FIXED_ENERGY_PARAM_FILE");
  if(!parameter||!vrna_params_load(parameter,VRNA_PARAMETER_FORMAT_DEFAULT))throw std::runtime_error("parameter load failed");
  sparse_tree tree(g,seq.size());SHAPEData shape("",seq.size());
  // Zero scaling seed and no sampling, as in the implementation's parity probe.
  W_final_pf pf(seq,g,shape,false,false,2,0.0,0,false,false,1.0);
  const double energy=pf.hfold_pf(tree),rt=pf.exp_params_->kT/1000.0;
  if(!std::isfinite(energy)||!std::isfinite(rt)||rt<=0)throw std::runtime_error("nonfinite oracle inside");
  std::cout<<std::setprecision(17)<<"{\"policy\":\"legacy-compatible-wmb\",\"free_energy\":"<<energy<<",\"rt\":"<<rt
      <<",\"log_z\":"<<-energy/rt<<"}\n";
  return 0;
  } catch(const std::exception &e) { std::cerr<<e.what()<<'\n';return 1; }
}
