// External observer linked to the pinned LinearCParty libraries.
#include "scfg2/api/exact.hh"
#include "scfg2/api/profiles.hh"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

int main(int argc, char **argv) {
  try {
    if (argc != 5 || std::string(argv[1]) != "-i" || std::string(argv[3]) != "--policy")
      throw std::runtime_error("usage: e0-linear -i record --policy legacy|default");
    std::ifstream input(argv[2]);
    std::string name, s, g;
    std::getline(input, name); std::getline(input, s); std::getline(input, g);
    const std::string choice = argv[4];
    if (choice != "legacy" && choice != "default") throw std::runtime_error("invalid policy");
    using namespace cparty::scfg2;
    const auto policy = choice == "legacy" ? WmbClosurePolicy::LegacyCompatible : WmbClosurePolicy::RequireWeakClosure;
    api::ExactSession session(s, g, policy);
    auto z = session.inside();
    auto b = session.outside();
    auto v = session.viterbi();
    auto p = api::exact_structural_profiles(s, g, {1}, policy);
    if (!z.ok || !z.reachable || !b.ok || !v.ok || !p.ok)
      throw std::runtime_error(z.error + b.error + v.error + p.error);
    std::cout << std::setprecision(17) << "{\"policy\":\"" << wmb_closure_policy_name(policy)
      << "\",\"log_z\":" << z.log_partition_weight << ",\"free_energy\":" << z.ensemble_free_energy
      << ",\"rt\":" << p.rt_kcal_mol << ",\"mfe\":" << v.energy << ",\"pairs\":[";
    bool first = true;
    for (const auto &pair : b.pair_probabilities) {
      if (!first) std::cout << ','; first = false;
      std::cout << '[' << pair.i << ',' << pair.j << ',' << pair.union_probability << ']';
    }
    std::cout << "],\"contexts\":[";
    for (std::size_t i = 1; i < p.contexts.size(); ++i) {
      if (i > 1) std::cout << ',';
      std::cout << '[';
      for (std::size_t c = 0; c < 7; ++c) { if (c) std::cout << ','; std::cout << p.contexts[i][c]; }
      std::cout << ']';
    }
    std::cout << "]}\n";
    return 0;
  } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
