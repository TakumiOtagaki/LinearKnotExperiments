// Inside-only comparison. Link the adopted libraries; leave their sources intact.
#include "scfg2/api/beam.hh"
#include "scfg2/weights/runtime_parameters.hh"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <unistd.h>

#ifndef CPARTY_ORACLE_PATH
#error CPARTY_ORACLE_PATH must name the locked E0 inside-only observer
#endif

int main(int argc, char **argv) {
  try {
    if (argc != 5 || std::string(argv[1]) != "-i" || std::string(argv[3]) != "--method")
      throw std::runtime_error("usage: PFScaling -i record --method cparty|exact|beam100");
    const std::string method = argv[4];
    if (method == "cparty") {
      execl(CPARTY_ORACLE_PATH, CPARTY_ORACLE_PATH, "-i", argv[2], nullptr);
      throw std::runtime_error("cannot execute locked CParty observer");
    }
    if (method != "exact" && method != "beam100") throw std::runtime_error("unknown method");
    std::ifstream input(argv[2]);
    std::string name, s, g;
    std::getline(input, name); std::getline(input, s); std::getline(input, g);
    if (!input || s.empty() || s.size() != g.size()) throw std::runtime_error("invalid input");
    using namespace cparty::scfg2;
    const auto policy = WmbClosurePolicy::LegacyCompatible;
    double log_z, energy;
    if (method == "exact") {
      api::ExactSession session(s, g, policy);
      const auto z = session.inside();
      if (!z.ok || !z.reachable) throw std::runtime_error(z.error);
      log_z = z.log_partition_weight; energy = z.ensemble_free_energy;
    } else {
      api::BeamConfig config;
      config.beam_size = 100;
      config.max_hairpin_unpaired = 100;
      config.max_case3_ancestor_depth = 100;
      config.wmb_closure_policy = policy;
      config.ranking = api::BeamRanking::Outer;
      api::BeamSession session(s, g, config);
      const auto z = session.inside();
      if (!z.ok || !z.reachable) throw std::runtime_error(z.error);
      log_z = z.log_partition_weight; energy = z.ensemble_free_energy;
    }
    std::unique_ptr<vrna_exp_param_t, decltype(&std::free)> parameters(
        weights::runtime_pf_parameters(2), &std::free);
    const double rt = parameters->kT / 1000.0;
    if (!std::isfinite(log_z) || !std::isfinite(energy) || !std::isfinite(rt) || rt <= 0)
      throw std::runtime_error("nonfinite partition function");
    std::cout << std::setprecision(17) << "{\"policy\":\"legacy-compatible-wmb\",\"method\":\""
              << method << "\",\"log_z\":" << log_z << ",\"free_energy\":" << energy
              << ",\"rt\":" << rt << "}\n";
    return 0;
  } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
