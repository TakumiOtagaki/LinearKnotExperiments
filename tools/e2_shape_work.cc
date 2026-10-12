// Optional diagnostics, separate from the three whole-CLI timing measurements.
// Link only the pinned beam/exact/RNA static libraries. One short record per call.
#include "scfg2/api/beam.hh"
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char **argv) {
  try {
    if (argc != 6 || std::string(argv[1]) != "-i")
      throw std::runtime_error("usage: e2-shape-work -i RECORD b h d");
    std::ifstream input(argv[2]);
    std::string name, sequence, scaffold;
    std::getline(input, name);
    std::getline(input, sequence);
    std::getline(input, scaffold);
    if (!input || name.empty() || name[0] != '>' || sequence.empty() || sequence.size() != scaffold.size())
      throw std::runtime_error("invalid named record");
    name.erase(0, 1);
    for (const unsigned char ch : name)
      if (!std::isalnum(ch) && ch != '_' && ch != '-') throw std::runtime_error("invalid record name");
    const auto integer = [](const char *text) {
      const std::string s(text);
      if (s.empty() || s.find_first_not_of("0123456789") != std::string::npos)
        throw std::runtime_error("invalid integer cap");
      return std::stoull(s);
    };
    cparty::scfg2::api::BeamConfig config;
    config.beam_size = integer(argv[3]);
    config.max_hairpin_unpaired = integer(argv[4]);
    config.max_case3_ancestor_depth = integer(argv[5]);
    const auto result = cparty::scfg2::api::beam_inside(sequence, scaffold, config);
    if (!result.ok || !result.reachable) throw std::runtime_error(result.error);
    const auto &w = result.work;
    std::cout << std::setprecision(17) << "{\"name\":\"" << name << "\",\"n\":" << sequence.size()
      << ",\"b\":" << config.beam_size << ",\"h\":" << config.max_hairpin_unpaired
      << ",\"d\":" << config.max_case3_ancestor_depth << ",\"log_z\":" << result.log_partition_weight
      << ",\"case3_truncated_endpoints\":" << w.case3_truncated_endpoints
      << ",\"item_visits\":" << w.total_item_visits << ",\"positive_items\":" << w.total_positive_items
      << ",\"frontier_candidates\":" << w.total_frontier_candidates
      << ",\"frontier_selected\":" << w.total_frontier_selected
      << ",\"vp_candidates\":" << w.vp_candidates << ",\"wmbp_candidates\":" << w.wmbp_candidates
      << ",\"case3_candidates\":" << w.case3_candidates
      << ",\"binary_join_probes\":" << w.binary_join_probes
      << ",\"v_parent_probes\":" << w.v_parent_probes
      << ",\"vp_loop_parent_probes\":" << w.vp_loop_parent_probes << ",\"families\":[";
    bool first = true;
    for (const auto &f : w.families) {
      if (!first) std::cout << ',';
      first = false;
      std::cout << "{\"nt\":\"" << f.nonterminal << "\",\"visits\":" << f.item_visits
        << ",\"positive\":" << f.positive_items << ",\"candidates\":" << f.frontier_candidates
        << ",\"selected\":" << f.frontier_selected << '}';
    }
    std::cout << "]}\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
