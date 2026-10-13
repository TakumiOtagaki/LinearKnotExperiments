// Link exclusively against the recorded upstream 800e95a + four-patch objects.
#include "h_globals.hh"
#include "part_func.hh"
#include "sparse_tree.hh"
#include <iomanip>
#include <iostream>
#include <string>
int main(int argc, char** argv) {
  if (argc != 6 || !vrna_params_load(argv[3], VRNA_PARAMETER_FORMAT_DEFAULT)) return 2;
  std::string s=argv[1], g=argv[2];
  sparse_tree tree(g, static_cast<cand_pos_t>(s.size()));
  SHAPEData shape("", static_cast<cand_pos_t>(s.size()));
  W_final_pf pf(s,g,shape,false,false,0,0.,0,false,false,1.);
  pf.hfold_pf(tree);
  const int i=std::stoi(argv[4]), j=std::stoi(argv[5]);
  std::cout << std::setprecision(17) << "{\"WI\":" << pf.get_energy_WI(i,j)
            << ",\"WIP\":" << pf.get_energy_WIP(i,j)
            << ",\"V\":" << pf.get_energy(i,j)
            << ",\"WI_empty\":" << pf.get_energy_WI(i,i-1)
            << ",\"rt\":" << pf.exp_params_->kT/1000.
            << ",\"pf_scale\":" << pf.exp_params_->pf_scale << "}\n";
}
