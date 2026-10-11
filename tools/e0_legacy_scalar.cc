// Full-precision observation linked to the uninstrumented seven-patch objects.
#include "h_globals.hh"
#include "part_func.hh"
#include "W_final.hh"
#include <fstream>
#include <iomanip>
#include <iostream>
int main(int argc,char **argv) {
  if(argc!=3||std::string(argv[1])!="-i")return 2;
  std::ifstream input(argv[2]);std::string name,seq,g;
  std::getline(input,name);std::getline(input,seq);std::getline(input,g);
  const char *parameter=std::getenv("CPARTY_FIXED_ENERGY_PARAM_FILE");
  if(!parameter||!vrna_params_load(parameter,VRNA_PARAMETER_FORMAT_DEFAULT))return 2;
  sparse_tree tree(g,seq.size());SHAPEData shape("",seq.size());
  W_final mfe(seq,g,shape,false,false,2);
  const double best=mfe.hfold(tree);
  W_final_pf pf(seq,g,shape,false,false,2,0.0,0,false,false,1.0);
  const double energy=pf.hfold_pf(tree),rt=pf.exp_params_->kT/1000.0;
  std::cout<<std::setprecision(17)<<"{\"free_energy\":"<<energy<<",\"rt\":"<<rt
      <<",\"log_z\":"<<-energy/rt<<",\"mfe\":"<<best<<"}\n";
}
