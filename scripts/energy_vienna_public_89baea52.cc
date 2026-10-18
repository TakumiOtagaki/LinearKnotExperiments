// Fresh public-API observation at dangles=2; dangles=0 uses the retained control.
#define main historical_internal_observer_main
#include "energy_vienna_internal_89baea52.cc"
#undef main
#include "scfg2/api/exact.hh"

int main(int argc, char **argv) {
  if (argc != 7) return 2;
  if (std::atoi(argv[4]) != 2) return historical_internal_observer_main(argc, argv);
  try {
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", argv[1], 1);
    const std::string sequence = argv[2], scaffold = argv[3], mode = argv[5];
    std::cout << std::setprecision(17);
    // Read thermodynamic metadata; all measured quantities below use ExactSession.
    {
      weights::RuntimeSnapshot owner(sequence, scaffold, weights::RuntimeOptions{2, false});
      const auto &p = owner.log_pf_parameters();
      std::cout << "META\t" << p.kT / 1000. << '\t' << p.temperature
                << "\t2\t" << p.model_details.pf_smooth << "\t3\t" << MAXLOOP << '\n';
    }
    api::ExactSession session(sequence, scaffold);
    if (mode == "ensemble") {
      const auto outside = session.outside();
      const auto mfe = session.viterbi();
      if (!outside.ok || !outside.reachable || !mfe.ok || !mfe.reachable)
        throw std::runtime_error("public exact calculation failed");
      std::cout << "Z\t" << outside.log_partition_weight << '\t' << outside.ensemble_free_energy << '\n';
      for (const auto &p : outside.pair_probabilities)
        std::cout << "BPP\t" << p.i << '\t' << p.j << '\t' << p.union_probability << '\n';
      std::cout << "MFE\t" << mfe.energy << '\t' << mfe.union_structure << '\n';
    } else if (mode != "eval") throw std::runtime_error("unknown observation mode");
    std::string target;
    std::size_t index = 0;
    while (std::getline(std::cin, target)) {
      if (target.empty()) continue;
      const auto result = session.evaluate_energy(target);
      std::vector<api::BasePair> wanted;
      std::string error;
      if (!api::detail::parse_dot_bracket_pairs(target, &wanted, &error))
        throw std::runtime_error(error);
      const auto a = api::detail::to_traceback_pairs(result.generated_pairs);
      const auto b = api::detail::to_traceback_pairs(wanted);
      std::cout << "EVAL\t" << index++ << '\t' << engine::same_traceback_pair_set(a, b)
                << '\t' << (result.ok && result.reachable) << '\t' << result.energy << '\n';
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
