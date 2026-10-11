// External observer linked to the pinned LinearCParty libraries.
#include "scfg2/api/exact.hh"
#include "scfg2/weights/runtime_parameters.hh"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>

int main(int argc, char **argv) {
  try {
    if (argc != 5 || std::string(argv[1]) != "-i" || std::string(argv[3]) != "--policy")
      throw std::runtime_error("usage: e0-linear -i record --policy legacy|default");
    std::ifstream input(argv[2]);
    std::string name, s, g;
    std::getline(input, name); std::getline(input, s); std::getline(input, g);
    if (!input || s.empty() || s.size() != g.size()) throw std::runtime_error("invalid input record");
    const std::string choice = argv[4];
    if (choice != "legacy" && choice != "default") throw std::runtime_error("invalid policy");
    using namespace cparty::scfg2;
    const auto policy = choice == "legacy" ? WmbClosurePolicy::LegacyCompatible : WmbClosurePolicy::RequireWeakClosure;
    api::ExactSession session(s, g, policy);
    auto z = session.inside();
    if (!z.ok || !z.reachable) throw std::runtime_error(z.error);
    // Read RT from the same parameter loader without computing another output.
    std::unique_ptr<vrna_exp_param_t, decltype(&std::free)> parameters(
        weights::runtime_pf_parameters(2, true), &std::free);
    const double rt = parameters->kT / 1000.0;
    if (!std::isfinite(rt) || rt <= 0) throw std::runtime_error("invalid RT");
    std::cout << std::setprecision(17) << "{\"policy\":\"" << wmb_closure_policy_name(policy)
      << "\",\"log_z\":" << z.log_partition_weight << ",\"free_energy\":" << z.ensemble_free_energy
      << ",\"rt\":" << rt << "}\n";
    return 0;
  } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
