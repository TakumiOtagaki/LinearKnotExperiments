// External observer: only public SCFG2 APIs; no implementation instrumentation.
#include "scfg2/api/beam.hh"
#include "scfg2/api/mea.hh"
#include "scfg2/api/profiles.hh"
#include "scfg2/api/sampling.hh"

// Keep the observer buildable against the historical pre-ranking API too.
#if __has_include("scfg2/beam/ranking_config.hh")
#define LCP_HAS_OUTER_RANKING 1
#else
#define LCP_HAS_OUTER_RANKING 0
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <locale>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <type_traits>

#ifndef LCP_SOURCE_REVISION
#define LCP_SOURCE_REVISION "unknown"
#endif

namespace api = cparty::scfg2::api;
namespace {
// Encode strict JSON without another dependency. Nonfinite numbers become null;
// ok/reachable/error preserve their meaning, including a disconnected root.
struct Json {
  std::string encoded;
  Json() : encoded("null") {}
  Json(bool x) : encoded(x ? "true" : "false") {}
  template <class T, std::enable_if_t<std::is_integral_v<T> &&
                                     !std::is_same_v<T, bool>, int> = 0>
  Json(T x) : encoded(std::to_string(x)) {}
  Json(double x) {
    if (!std::isfinite(x)) { encoded = "null"; return; }
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::setprecision(17) << x;
    encoded = out.str();
  }
  Json(const char *x) : Json(std::string(x)) {}
  Json(const std::string &x) {
    std::ostringstream out;
    out << '"';
    for (unsigned char c : x) {
      if (c == '"' || c == '\\') out << '\\' << c;
      else if (c < 32) out << "\\u" << std::hex << std::setw(4)
                           << std::setfill('0') << static_cast<int>(c);
      else out << c;
    }
    out << '"';
    encoded = out.str();
  }
  static Json array(const std::vector<Json> &xs) {
    Json out; out.encoded = "[";
    for (const auto &x : xs) {
      if (out.encoded.size() > 1) out.encoded += ',';
      out.encoded += x.encoded;
    }
    out.encoded += ']'; return out;
  }
  static Json object(std::initializer_list<std::pair<std::string, Json>> xs) {
    Json out; out.encoded = "{";
    for (const auto &x : xs) {
      if (out.encoded.size() > 1) out.encoded += ',';
      out.encoded += Json(x.first).encoded + ':' + x.second.encoded;
    }
    out.encoded += '}'; return out;
  }
};

template <class Range> Json numbers(const Range &xs, std::size_t begin = 0) {
  std::vector<Json> values;
  for (std::size_t i = begin; i < xs.size(); ++i) values.emplace_back(xs[i]);
  return Json::array(values);
}

Json rss() {
  std::ifstream in("/proc/self/status");
  Json resident, peak;
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream row(line);
    std::string key, unit; std::uint64_t value;
    if (row >> key >> value >> unit && unit == "kB") {
      if (key == "VmRSS:") resident = Json(value);
      if (key == "VmHWM:") peak = Json(value);
    }
  }
  return Json::object({{"vmrss_kib", resident}, {"vmhwm_kib", peak}});
}

struct Options {
  std::string input, sequence, scaffold, name = "seq1";
  std::string session_mode = "fresh", stage = "all";
  api::BeamConfig beam;
  std::vector<int> lengths{1, 2, 4, 7};
  std::size_t count = 10;
  std::uint64_t seed = 1;
  double gamma = 1;
  bool values = false;
  bool inside_repeat = true, pilot_explicit = false;
#if LCP_HAS_OUTER_RANKING
  std::string ranking = "outer";
#else
  std::string ranking = "baseline";
#endif
  std::string pilot_exponent = "2/3";
};

std::uint64_t integer(const std::string &s) {
  if (s.empty() || !std::all_of(s.begin(), s.end(), [](unsigned char c) {
        return std::isdigit(c); })) throw std::invalid_argument("expected unsigned integer: " + s);
  std::size_t end; const auto v = std::stoull(s, &end);
  if (end != s.size()) throw std::invalid_argument("invalid integer: " + s);
  return v;
}

std::string trim(const std::string &s) {
  const auto a = s.find_first_not_of(" \t\r\n");
  return a == std::string::npos ? "" : s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}

void read_one(Options &o) {
  std::ifstream in(o.input);
  if (!in) throw std::invalid_argument("cannot open input: " + o.input);
  std::vector<std::string> lines; std::string line;
  while (std::getline(in, line)) {
    line = trim(line);
    if (!line.empty() && line[0] != '#') lines.push_back(line);
  }
  if (!lines.empty() && lines[0][0] == '>') {
    o.name = trim(lines[0].substr(1)); lines.erase(lines.begin());
  }
  if (lines.size() != 2) throw std::invalid_argument("-i requires exactly one sequence/scaffold record");
  o.sequence = lines[0]; o.scaffold = lines[1];
  for (auto &c : o.sequence) {
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
    if (c == 'T') c = 'U';
  }
}

Options parse(int argc, char **argv) {
  Options o;
  for (int i = 1; i < argc; ++i) {
    const std::string key = argv[i];
    const auto value = [&]() -> std::string {
      if (++i == argc) throw std::invalid_argument("missing value for " + key);
      return argv[i];
    };
    if (key == "-i") o.input = value();
    else if (key == "--sequence") o.sequence = value();
    else if (key == "--scaffold") o.scaffold = value();
    else if (key == "-b" || key == "--beam") o.beam.beam_size = integer(value());
    else if (key == "--max-hairpin") o.beam.max_hairpin_unpaired = integer(value());
    else if (key == "--max-case3-depth") o.beam.max_case3_ancestor_depth = integer(value());
    else if (key == "-k" || key == "--count") o.count = integer(value());
    else if (key == "--seed") o.seed = integer(value());
    else if (key == "--session-mode") o.session_mode = value();
    else if (key == "--stage") o.stage = value();
    else if (key == "--include-values") o.values = true;
    else if (key == "--skip-inside-repeat") o.inside_repeat = false;
    else if (key == "--ranking") o.ranking = value();
    else if (key == "--outer-pilot-exponent") {
      o.pilot_exponent = value(); o.pilot_explicit = true;
    }
    else if (key == "--gamma") {
      const auto s = value(); std::size_t used;
      o.gamma = std::stod(s, &used);
      if (used != s.size() || !std::isfinite(o.gamma) || o.gamma < 0)
        throw std::invalid_argument("gamma must be finite and nonnegative");
    } else if (key == "--access-len") {
      const auto s = value(); o.lengths.clear();
      if (s.empty() || s.back() == ',') throw std::invalid_argument("invalid access lengths");
      std::istringstream parts(s); std::string token;
      while (std::getline(parts, token, ',')) {
        const auto w = integer(token);
        if (!w || w > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
          throw std::invalid_argument("access length out of range");
        if (std::find(o.lengths.begin(), o.lengths.end(), w) != o.lengths.end())
          throw std::invalid_argument("duplicate access length");
        o.lengths.push_back(static_cast<int>(w));
      }
    } else throw std::invalid_argument("unknown option: " + key);
  }
  if (!o.input.empty()) {
    if (!o.sequence.empty() || !o.scaffold.empty()) throw std::invalid_argument("use -i or sequence/scaffold");
    read_one(o);
  } else if (o.sequence.empty() || o.scaffold.empty()) throw std::invalid_argument("supply sequence and scaffold");
  if (!o.beam.beam_size) throw std::invalid_argument("beam must be positive");
  if (o.ranking != "outer" && o.ranking != "baseline")
    throw std::invalid_argument("ranking must be outer or baseline");
  if (o.pilot_exponent != "0.5" && o.pilot_exponent != "2/3")
    throw std::invalid_argument("outer pilot exponent must be 0.5 or 2/3");
  if (o.pilot_explicit && o.ranking != "outer")
    throw std::invalid_argument("outer pilot exponent requires outer ranking");
#if LCP_HAS_OUTER_RANKING
  o.beam.ranking = o.ranking == "outer" ? api::BeamRanking::Outer : api::BeamRanking::Baseline;
  o.beam.outer_pilot_exponent = o.pilot_exponent == "0.5" ?
      api::OuterPilotExponent::Half : api::OuterPilotExponent::TwoThirds;
#else
  if (o.ranking != "baseline" || o.pilot_explicit)
    throw std::invalid_argument("this source revision has no outer ranking API");
#endif
  if (o.session_mode != "fresh" && o.session_mode != "shared" && o.session_mode != "both")
    throw std::invalid_argument("session mode must be fresh, shared or both");
  const std::vector<std::string> stages{"all", "inside", "mfe", "outside_bpp", "loop_profile", "accessibility", "sampling", "mea"};
  if (std::find(stages.begin(), stages.end(), o.stage) == stages.end()) throw std::invalid_argument("unknown stage");
  return o;
}

Json work(const api::BeamWorkStats &w) {
  std::vector<Json> families;
  for (const auto &f : w.families)
    families.push_back(Json::object({{"nonterminal", f.nonterminal},
      {"item_visits", f.item_visits}, {"positive_items", f.positive_items},
      {"frontier_candidates", f.frontier_candidates}, {"frontier_selected", f.frontier_selected}}));
  return Json::object({{"item_visits", w.total_item_visits},
    {"frontier_candidates", w.total_frontier_candidates}, {"frontier_selected", w.total_frontier_selected},
    {"case3_truncated_endpoints", w.case3_truncated_endpoints},
    {"positive_items", w.total_positive_items},
    {"vp_candidates", w.vp_candidates}, {"vp_selected", w.vp_selected},
    {"wmbp_candidates", w.wmbp_candidates}, {"wmbp_selected", w.wmbp_selected},
    {"case3_candidates", w.case3_candidates}, {"case3_selected", w.case3_selected},
    {"binary_join_probes", w.binary_join_probes}, {"extension_probes", w.extension_probes},
    {"v_parent_probes", w.v_parent_probes}, {"vp_loop_parent_probes", w.vp_loop_parent_probes},
    {"families", Json::array(families)}});
}
Json summary(const api::BeamInsideResult &r, const Options &) {
  return Json::object({{"log_z", r.log_partition_weight}, {"free_energy", r.ensemble_free_energy}, {"work", work(r.work)}});
}
Json summary(const api::BeamViterbiResult &r, const Options &) {
  return Json::object({{"log_weight", r.log_weight}, {"energy", r.energy},
    {"generated_pair_count", r.generated_pairs.size()}, {"generated_structure", r.generated_structure},
    {"union_structure", r.union_structure}, {"work", work(r.work)}});
}
Json summary(const api::BeamOutsideResult &r, const Options &o) {
  double sum = 0, round = 0, square = 0;
  std::vector<Json> pairs;
  for (const auto &p : r.pair_probabilities) {
    sum += p.union_probability; round += p.round_probability; square += p.square_probability;
    if (o.values) pairs.push_back(Json::array({p.i, p.j, p.round_probability, p.square_probability, p.union_probability}));
  }
  const double unpaired = r.unpaired_probabilities.empty() ? 0 :
    std::accumulate(r.unpaired_probabilities.begin() + 1, r.unpaired_probabilities.end(), 0.0);
  return Json::object({{"log_z", r.log_partition_weight}, {"free_energy", r.ensemble_free_energy},
    {"pair_count", r.pair_probabilities.size()}, {"bpp_sum", sum}, {"round_bpp_sum", round}, {"square_bpp_sum", square},
    {"unpaired_sum", unpaired}, {"reverse_records", r.reverse_records}, {"work", work(r.work)},
    {"pairs", o.values ? Json::array(pairs) : Json()},
    {"unpaired", o.values ? numbers(r.unpaired_probabilities, 1) : Json()}});
}
Json summary(const api::StructuralProfileResult &r, const Options &o) {
  std::array<double, api::structural_context_count> sums{};
  double norm_error = 0; std::vector<Json> contexts, windows;
  for (std::size_t i = 1; i < r.contexts.size(); ++i) {
    double total = 0;
    for (std::size_t c = 0; c < sums.size(); ++c) { sums[c] += r.contexts[i][c]; total += r.contexts[i][c]; }
    norm_error = std::max(norm_error, std::abs(total - 1));
    if (o.values) contexts.push_back(numbers(r.contexts[i]));
  }
  for (std::size_t k = 0; k < r.access_lengths.size(); ++k) {
    const auto w = static_cast<std::size_t>(r.access_lengths[k]);
    const auto count = w > o.sequence.size() ? 0 : o.sequence.size() - w + 1;
    const bool available = r.ok && k < r.accessibility.size() &&
                           r.accessibility[k].size() > count;
    double sum = 0; std::vector<Json> values;
    for (std::size_t i = 1; available && i <= count; ++i) {
      sum += r.accessibility[k][i];
      if (o.values) values.emplace_back(r.accessibility[k][i]);
    }
    windows.push_back(Json::object({{"length", w}, {"window_count", count},
      {"probability_sum", available ? Json(sum) : Json()},
      {"probabilities", o.values && available ? Json::array(values) : Json()}}));
  }
  return Json::object({{"log_z", r.log_partition_weight}, {"free_energy", r.ensemble_free_energy},
    {"rt_kcal_mol", r.rt_kcal_mol}, {"position_count", r.contexts.empty() ? 0 : r.contexts.size() - 1},
    {"context_sums", numbers(sums)}, {"max_context_normalization_error", norm_error},
    {"contexts", o.values ? Json::array(contexts) : Json()}, {"accessibility", Json::array(windows)}});
}
Json summary(const api::SampleResult &r, const Options &o) {
  std::size_t pairs = 0; std::vector<Json> structures;
  for (const auto &sample : r.generated_pairs) pairs += sample.size();
  if (o.values) for (const auto &s : r.union_structures) structures.emplace_back(s);
  return Json::object({{"log_z", r.log_partition_weight}, {"free_energy", r.ensemble_free_energy},
    {"sample_count", r.union_structures.size()}, {"generated_pair_count_sum", pairs},
    {"union_structures", o.values ? Json::array(structures) : Json()}});
}
Json summary(const api::MeaResult &r, const Options &) {
  return Json::object({{"log_z", r.log_partition_weight}, {"expected_accuracy", r.expected_accuracy},
    {"generated_pair_count", r.generated_pairs.size()}, {"generated_structure", r.generated_structure},
    {"union_structure", r.union_structure}, {"decoding_states", r.decoding_states},
    {"decoding_terms", r.decoding_terms}, {"decoding_factors", r.decoding_factors}});
}

template <class T> Json reachable(const T &) { return Json(); }
Json reachable(const api::BeamInsideResult &r) { return r.reachable; }
Json reachable(const api::BeamOutsideResult &r) { return r.reachable; }
Json reachable(const api::BeamViterbiResult &r) { return r.reachable; }
Json reachable(const api::MeaResult &r) { return r.reachable; }
template <class T> Json internal_seconds(const T &r) { return r.elapsed_seconds; }
Json internal_seconds(const api::BeamViterbiResult &) { return Json(); }
template <class T> Json decoding_seconds(const T &) { return Json(); }
Json decoding_seconds(const api::MeaResult &r) { return r.decoding_seconds; }

template <class Fn>
Json measure(const std::string &name, const std::string &call, const std::string &scope,
             bool cached, const Options &o, bool &ok, Fn fn) {
  const auto before = rss();
  const auto start = std::chrono::steady_clock::now();
  const auto result = fn();
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  const auto after = rss();  // Returned values still alive; API temporaries were destroyed.
  ok = ok && result.ok && reachable(result).encoded != "false";
  return Json::object({{"stage", name}, {"api_call", call}, {"scope", scope}, {"ok", result.ok},
    {"reachable", reachable(result)}, {"error", result.error}, {"elapsed_seconds", seconds},
    {"rss_before", before}, {"rss_after", after}, {"graph_reused", false}, {"inside_result_cached", cached},
    {"api_elapsed_seconds", cached ? Json() : internal_seconds(result)},
    {"api_decoding_seconds", decoding_seconds(result)}, {"summary", summary(result, o)}});
}

Json run(const Options &o, bool shared, bool &ok) {
  api::BeamSession session(o.sequence, o.scaffold, o.beam);
  std::vector<Json> stages;
  const auto want = [&](const std::string &s) { return o.stage == "all" || o.stage == s; };
  if (want("inside")) {
    for (int repeat = 0; repeat < (o.inside_repeat ? 2 : 1); ++repeat)
      stages.push_back(measure(repeat ? "inside_repeat" : "inside",
        shared ? "BeamSession::inside" : "beam_inside",
        shared && repeat ? "cached result copy; no forward pass" : "numeric forward and preparation",
        shared && repeat, o, ok, [&] { return shared ? session.inside() : api::beam_inside(o.sequence, o.scaffold, o.beam); }));
  }
  if (want("mfe")) stages.push_back(measure("mfe", shared ? "BeamSession::viterbi" : "beam_viterbi",
    "independent LogMax forward, traceback and preparation", false, o, ok,
    [&] { return shared ? session.viterbi() : api::beam_viterbi(o.sequence, o.scaffold, o.beam); }));
  if (want("outside_bpp")) stages.push_back(measure("outside_bpp", shared ? "BeamSession::outside" : "beam_outside",
    "numeric forward, regenerated outside, BPP aggregation and preparation", false, o, ok,
    [&] { return shared ? session.outside() : api::beam_outside(o.sequence, o.scaffold, o.beam); }));
  if (want("loop_profile")) stages.push_back(measure("loop_profile", "beam_structural_profiles",
    "numeric forward, regenerated outside, context aggregation and preparation", false, o, ok,
    [&] { return api::beam_structural_profiles(o.sequence, o.scaffold, {}, o.beam); }));
  if (want("accessibility")) stages.push_back(measure("accessibility", "beam_structural_profiles",
    "numeric forward, regenerated outside, context and window aggregation and preparation", false, o, ok,
    [&] { return api::beam_structural_profiles(o.sequence, o.scaffold, o.lengths, o.beam); }));
  if (want("sampling")) stages.push_back(measure("sampling", "beam_sample",
    "factorized forward, sampler preparation and all draws; no outside", false, o, ok,
    [&] { return api::beam_sample(o.sequence, o.scaffold, o.count, o.seed, o.beam); }));
  if (want("mea")) stages.push_back(measure("mea", "beam_mea",
    "factorized forward, outside, BPP aggregation, rewards, maximization and traceback", false, o, ok,
    [&] { return api::beam_mea(o.sequence, o.scaffold, o.gamma, o.beam); }));
  return Json::object({{"session_mode", shared ? "shared" : "fresh"}, {"stages", Json::array(stages)}});
}
}  // namespace

int main(int argc, char **argv) {
  if (argc == 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
    std::cout << "stage_timing_driver -i ONE_RECORD | --sequence S --scaffold G\n"
      "  -b N --max-hairpin N --max-case3-depth N --access-len 1,2,4,7\n"
      "  -k N --seed N --gamma X --session-mode fresh|shared|both\n"
      "  --ranking outer|baseline --outer-pilot-exponent 0.5|2/3\n"
      "  --stage all|inside|mfe|outside_bpp|loop_profile|accessibility|sampling|mea\n"
      "  --skip-inside-repeat  (one inside call per run; omit cache diagnostic)\n"
      "  --include-values  (full numerical output for small validation inputs)\n";
    return 0;
  }
  try {
    const auto o = parse(argc, argv);
    bool ok = true; std::vector<Json> runs;
    if (o.session_mode != "shared") runs.push_back(run(o, false, ok));
    if (o.session_mode != "fresh") runs.push_back(run(o, true, ok));
    const char *parameter = std::getenv("CPARTY_FIXED_ENERGY_PARAM_FILE");
    std::cout << Json::object({{"schema_version", 1}, {"ok", ok}, {"source_revision", LCP_SOURCE_REVISION},
      {"clock", "std::chrono::steady_clock"}, {"rss_source", "/proc/self/status; process lifetime VmHWM"},
      {"input", Json::object({{"name", o.name}, {"sequence", o.sequence}, {"scaffold", o.scaffold}, {"length", o.sequence.size()}})},
      {"config", Json::object({{"beam", o.beam.beam_size}, {"max_hairpin", o.beam.max_hairpin_unpaired},
        {"max_case3_depth", o.beam.max_case3_ancestor_depth}, {"access_lengths", numbers(o.lengths)},
        {"sample_count", o.count}, {"seed", o.seed}, {"gamma", o.gamma}, {"include_values", o.values},
        {"ranking", o.ranking},
        {"outer_pilot_exponent", o.ranking == "outer" ? Json(o.pilot_exponent) : Json()},
        {"inside_repeat", o.inside_repeat},
        {"stage", o.stage}, {"wmb_closure_policy", "RequireWeakClosure"}, {"temperature_celsius", 37},
        {"energy_parameter_file", parameter ? Json(parameter) : Json()}})},
      {"reuse_contract", "Only inside result is cached. No public graph handoff; other calls rebuild, as in the CLIs."},
      {"runs", Json::array(runs)}}).encoded << '\n';
    return ok ? 0 : 1;
  } catch (const std::exception &e) {
    std::cout << Json::object({{"schema_version", 1}, {"ok", false}, {"error", e.what()}}).encoded << '\n';
    std::cerr << "stage_timing_driver: " << e.what() << '\n';
    return 2;
  }
}
