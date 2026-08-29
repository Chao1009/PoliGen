// Generates T0 spectator-tagged events -- 6Li(e, e' alpha)X on the embedded
// deuteron, 7Li(e, e' alpha)X on the quasi-free triton, or the d(e, e' p)X
// Cosyn-Weiss control -- through lipolgen::Pipeline, writes them as HepMC3
// Asciiv3, and prints the accepted cross section per spin category, the
// Roman-Pot tag fraction at every optics of the menu, and the throughput.
//
// Usage:
//   generate_tagged [options]
//     --channel 6Li-alpha|7Li-alpha|d-p   tagged channel      (default 6Li-alpha)
//     --config 0|1|2        energy point low/mid/top          (default 1)
//     --plan NAME           azz | apar | pure                 (default azz)
//     --optics NAME         high-acceptance | high-divergence |
//                           tagging | tagging-legacy          (default
//                                                              high-acceptance)
//     --events N            fixed event count                 (default 200000)
//     --lumi L              integrated luminosity [pb^-1] instead of --events
//     --seed S              RNG seed                          (default 20260713)
//     --threads N           worker threads (output is identical for any N)
//     --beta B              short-range scale of the radial waves (0.20-0.40)
//     --pd P                alpha-d D-state probability (6Li)
//     --inclusive-b1        put an inclusive b1 in the struck-cluster kernel
//                           (double-counts the wave-function tensor
//                            asymmetry -- for the RATE identity only)
//     --out FILE            HepMC3 output ("" = none)
//
// The route is RECOMPUTED per optics from the finished record (`route_of`),
// which is why one generated sample can be priced at the whole menu at once.
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "lipolgen/hepmc_writer.hpp"
#include "lipolgen/pipeline.hpp"

using namespace lipolgen;

namespace {

struct Args {
  std::string channel = "6Li-alpha";
  int config = 1;
  std::string plan = "azz";
  std::string optics = "high-acceptance";
  long long events = 200000;
  double lumi = -1.0;
  unsigned long long seed = 20260713;
  unsigned threads = 1;
  double beta = BETA_DEFAULT;
  double p_d = P_D_LI6;
  bool inclusive_b1 = false;
  std::string out = "generate_tagged.hepmc";
};

[[noreturn]] void usage(const char* argv0, const std::string& why) {
  std::cerr << argv0 << ": " << why << "\n"
            << "usage: " << argv0
            << " [--channel 6Li-alpha|7Li-alpha|d-p] [--config 0|1|2]\n"
               "       [--plan azz|apar|pure] [--optics high-acceptance|"
               "high-divergence|tagging|tagging-legacy]\n"
               "       [--events N | --lumi PB] [--seed S] [--threads N]"
               " [--beta B] [--pd P]\n"
               "       [--inclusive-b1] [--out FILE]\n";
  std::exit(2);
}

Args parse(int argc, char** argv) {
  Args a;
  for (int i = 1; i < argc; ++i) {
    const std::string k = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) usage(argv[0], "missing value for " + k);
      return argv[++i];
    };
    if (k == "--channel") a.channel = next();
    else if (k == "--config") a.config = std::stoi(next());
    else if (k == "--plan") a.plan = next();
    else if (k == "--optics") a.optics = next();
    else if (k == "--events") { a.events = std::stoll(next()); a.lumi = -1.0; }
    else if (k == "--lumi") { a.lumi = std::stod(next()); a.events = -1; }
    else if (k == "--seed") a.seed = std::stoull(next());
    else if (k == "--threads") a.threads = static_cast<unsigned>(std::stoul(next()));
    else if (k == "--beta") a.beta = std::stod(next());
    else if (k == "--pd") a.p_d = std::stod(next());
    else if (k == "--inclusive-b1") a.inclusive_b1 = true;
    else if (k == "--out") a.out = next();
    else if (k == "-h" || k == "--help") usage(argv[0], "help");
    else usage(argv[0], "unknown option " + k);
  }
  if (a.config < 0 || a.config > 2) usage(argv[0], "--config must be 0, 1 or 2");
  if (a.threads < 1) usage(argv[0], "--threads must be >= 1");
  return a;
}

PipelineChannel channel_of(const std::string& s, const char* argv0) {
  if (s == "6Li-alpha") return PipelineChannel::TaggedLi6Alpha;
  if (s == "7Li-alpha") return PipelineChannel::TaggedLi7Alpha;
  if (s == "d-p") return PipelineChannel::TaggedDeuteronP;
  usage(argv0, "unknown --channel " + s);
}

OpticsChoice optics_of(const std::string& s, const char* argv0) {
  if (s == "high-acceptance") return OpticsChoice::YellowReportHighAcceptance;
  if (s == "high-divergence") return OpticsChoice::YellowReportHighDivergence;
  if (s == "tagging") return OpticsChoice::Tagging;
  if (s == "tagging-legacy") return OpticsChoice::TaggingLegacyLevers;
  usage(argv0, "unknown --optics " + s);
}

RunPlan plan_of(const std::string& name, double j, const Scenario& sc,
                const char* argv0) {
  if (name == "apar") {
    return helicity_flip_plan(j, sc.pol_ion_vector, sc.pol_electron);
  }
  if (name == "pure") {
    std::vector<double> pops(static_cast<std::size_t>(2 * j + 1.5), 0.0);
    pops[0] = 1.0;   // the stretched state M = +J
    return RunPlan({SpinCategory("pure", j, pops, 0, 0.0, 0.0, 0.0, 1.0)}, 0.0,
                   1.0, j >= 1.0 ? 1.0 : 0.0);
  }
  if (name == "azz") {
    if (std::fabs(j - 1.0) > 1e-9) {
      usage(argv0, "--plan azz is spin-1 only; use apar or pure");
    }
    return tensor_thirds_plan(sc.pol_ion_vector, sc.pol_ion_tensor);
  }
  usage(argv0, "unknown --plan " + name);
}

}  // namespace

int main(int argc, char** argv) {
  const Args args = parse(argc, argv);
  try {
    PipelineConfig cfg;
    cfg.channel = channel_of(args.channel, argv[0]);
    cfg.isotope = channel_isotope(cfg.channel);
    cfg.beam_config = args.config;
    cfg.optics_choice = optics_of(args.optics, argv[0]);
    cfg.seed = args.seed;
    cfg.cluster_beta = args.beta;
    cfg.p_d = args.p_d;
    cfg.struck.inclusive_b1 = args.inclusive_b1;
    if (args.lumi > 0.0) cfg.lumi_pb = args.lumi;
    else cfg.n_events = static_cast<std::uint64_t>(args.events);

    const Ion& ion = ion_by_name(cfg.isotope);
    const RunPlan plan = plan_of(args.plan, ion.spin, cfg.scenario, argv[0]);

    const auto t_build0 = std::chrono::steady_clock::now();
    const Pipeline p(cfg, plan);
    const auto t_build1 = std::chrono::steady_clock::now();

    const BeamConfig& bc = p.beam_config();
    const TaggedChannel& ch = *p.tagged_channel();
    std::printf("LiPolGen tagged T0 generator -- %s\n", ch.label.c_str());
    std::printf("  beams        %s (sqrt(s_eN) = %.2f GeV), pot config %s\n",
                bc.label().c_str(), bc.sqrt_s_per_nucleon(),
                p.pot_config().c_str());
    std::printf("  cluster      spectator %s (A=%d, Z=%d), struck %s,"
                " kappa = %.4f GeV, beta = %.2f\n",
                ch.base.spectator.c_str(), ch.base.spectator_A,
                ch.base.spectator_Z, ch.dis_target.name.c_str(),
                ch.base.kappa(), args.beta);
    std::printf("  waves       ");
    for (const Wave& w : ch.waves) std::printf(" L=%d P=%.4f", w.l, w.prob);
    std::printf("   R(k=0) = %.5f, vector dilution %.4f\n",
                ch.base.r_at_k_zero(), p.tagged_model()->vector_dilution());
    std::printf("  window       Q2 >= %.2f, y in [%.3f, %.3f], W2 >= %.1f;"
                " %zu accepted (x, Q2) cells on the struck cluster\n",
                cfg.scenario.q2_min, cfg.scenario.y_min, cfg.scenario.y_max,
                cfg.scenario.w2_min, p.dis_sampler().n_cells());
    std::printf("  optics       %s: %.3f x %.3f mrad envelope, L/L_HA = %.4f\n",
                p.optics().name.c_str(), 1e3 * p.optics().envelope_x(),
                1e3 * p.optics().envelope_y(), p.optics().lumi_fraction);
    std::printf("  build        %.3f s\n",
                std::chrono::duration<double>(t_build1 - t_build0).count());

    std::printf("\n  %-10s %10s %6s %16s %14s\n", "category", "lumi_frac",
                "lam_e", "sigma [pb]", "events");
    for (std::size_t k = 0; k < plan.categories().size(); ++k) {
      const SpinCategory& c = plan.categories()[k];
      std::printf("  %-10s %10.5f %6d %16.6g %14llu\n", c.name.c_str(),
                  c.lumi_fraction, c.lam_e, p.sigma_per_category_pb()[k],
                  static_cast<unsigned long long>(p.counts()[k]));
    }
    std::printf("  %-10s %10.5f %6s %16.6g %14llu\n", "TOTAL", 1.0, "",
                p.sigma_pb(), static_cast<unsigned long long>(p.size()));

    // the whole optics menu, priced on the SAME sample
    struct Menu { const char* label; OpticsChoice choice; };
    const Menu menu[] = {
        {"YR high-acceptance", OpticsChoice::YellowReportHighAcceptance},
        {"YR high-divergence", OpticsChoice::YellowReportHighDivergence},
        {"tagging", OpticsChoice::Tagging},
        {"tagging (18x275 levers)", OpticsChoice::TaggingLegacyLevers}};
    std::vector<Optics> envelopes;
    std::vector<std::uint64_t> tagged_counts;
    for (const Menu& m : menu) {
      envelopes.push_back(optics_for(m.choice, cfg.isotope,
                                     bc.ion_momentum_per_nucleon));
      tagged_counts.push_back(0);
    }

    std::unique_ptr<HepMC3Writer> writer;
    if (!args.out.empty()) writer = std::make_unique<HepMC3Writer>(args.out);

    std::uint64_t written = 0;
    double worst_p = 0.0, worst_q = 0.0, sum_k = 0.0;
    const auto sink = [&](const Event& ev) {
      const Vec4 r = momentum_residual(ev);
      const Vec4 s = momentum_scale(ev);
      worst_p = std::max(worst_p,
                         std::max(std::fabs(r.e), std::fabs(r.pz)) / s.e);
      worst_q = std::max(worst_q, std::fabs(charge_residual(ev)));
      sum_k += ev.kin.k;
      for (std::size_t i = 0; i < envelopes.size(); ++i) {
        tagged_counts[i] += rp_tagged(ev, envelopes[i], p.pot_config());
      }
      if (writer) writer->write(ev);
      ++written;
    };

    const auto t0 = std::chrono::steady_clock::now();
    const std::uint64_t n = args.threads > 1 ? p.for_each(sink, args.threads)
                                             : p.for_each(sink);
    const auto t1 = std::chrono::steady_clock::now();
    if (writer) writer->close();
    const double secs = std::chrono::duration<double>(t1 - t0).count();

    std::printf("\n  %-26s %12s %12s %12s\n", "optics", "envelope [mrad]",
                "RP tag", "tag x L/L_HA");
    for (std::size_t i = 0; i < envelopes.size(); ++i) {
      const double f = n ? static_cast<double>(tagged_counts[i]) / n : 0.0;
      char env[32];
      std::snprintf(env, sizeof(env), "%.3fx%.3f", 1e3 * envelopes[i].envelope_x(),
                    1e3 * envelopes[i].envelope_y());
      std::printf("  %-26s %12s %12.4f %12.4f\n", menu[i].label, env, f,
                  f * envelopes[i].lumi_fraction);
    }

    std::printf("\n  generated    %llu events (seed %llu, %u thread%s)\n",
                static_cast<unsigned long long>(n),
                static_cast<unsigned long long>(args.seed), args.threads,
                args.threads == 1 ? "" : "s");
    std::printf("  time         %.3f s -> %.4g events/s%s\n", secs,
                secs > 0.0 ? static_cast<double>(n) / secs : 0.0,
                writer ? " (including HepMC3 output)" : "");
    std::printf("  <k>          %.4f GeV/c (spectator, ion rest frame)\n",
                n ? sum_k / static_cast<double>(n) : 0.0);
    std::printf("  worst |k + P_ion - k' - p_spec - X| / E_in = %.3g,"
                " |dQ| = %.3g\n", worst_p, worst_q);
    if (writer) {
      std::printf("  wrote        %llu events to %s\n",
                  static_cast<unsigned long long>(written), args.out.c_str());
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "generate_tagged: " << e.what() << std::endl;
    return 1;
  }
}
