// Generates T0 coherent diffractive events -- e + 6Li -> e' + X + 6Li(g.s.) --
// through lipolgen::Pipeline, writes them as HepMC3 Asciiv3, and prints the
// accepted cross section per spin category, the |t| spectrum against the
// exponential slope, the near-beam Roman-Pot tag fraction at every optics of
// the menu, and the throughput.
//
// EVERYTHING HERE IS A SCENARIO with explicit bands (coherent.hpp): the
// coherent fraction f0/(1 + (x/x_coh)^2), the |F(t)|^2 slope B, the flat
// gluon-transversity cos 2phi amplitude and the m = 0 slope modulation
// eps_b0.  No published calculation exists for any polarized A > 2 nucleus.
//
// Usage:
//   generate_coherent [options]
//     --config 0|1|2        energy point low/mid/top          (default 1)
//     --plan azz|cos2phi|flip                                 (default azz)
//     --optics NAME         high-acceptance | high-divergence |
//                           tagging | tagging-legacy          (default tagging
//                                                              -- the Yellow
//                             Report envelope sees nothing at all)
//     --events N            fixed event count                 (default 200000)
//     --lumi L              integrated luminosity [pb^-1] instead of --events
//     --f0 F                coherent fraction at x -> 0       (band 0.02-0.08)
//     --slope B             |F(t)|^2 slope [GeV^-2]           (band 40-60)
//     --amp A               flat cos 2phi amplitude at P_zz=1 (band 3e-3-1e-2)
//     --eps-b0 E            m = 0 relative slope modulation   (band -0.04..-0.13)
//     --tmax T              |t| truncation [GeV^2]            (default 0.2; see P4)
//     --weighted-azimuth    draw phi_t from the modulated density instead of
//                           carrying it as an event weight
//     --seed S              RNG seed                          (default 20260713)
//     --threads N           worker threads (output identical for any N)
//     --out FILE            HepMC3 output ("" = none)
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
  int config = 1;
  std::string plan = "azz";
  std::string optics = "tagging";
  long long events = 200000;
  double lumi = -1.0;
  double t_max = COHERENT_T_MAX_DEFAULT;   // P4: 0.5 makes the weight negative
  bool weighted_azimuth = false;
  unsigned long long seed = 20260713;
  unsigned threads = 1;
  std::string out = "generate_coherent.hepmc";
  CoherentScenario sc;
};

[[noreturn]] void usage(const char* argv0, const std::string& why) {
  std::cerr << argv0 << ": " << why << "\n"
            << "usage: " << argv0
            << " [--config 0|1|2] [--plan azz|cos2phi|flip]\n"
               "       [--optics high-acceptance|high-divergence|tagging|"
               "tagging-legacy]\n"
               "       [--events N | --lumi PB] [--f0 F] [--slope B] [--amp A]"
               " [--eps-b0 E]\n"
               "       [--tmax T] [--weighted-azimuth] [--seed S]"
               " [--threads N] [--out FILE]\n";
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
    if (k == "--config") a.config = std::stoi(next());
    else if (k == "--plan") a.plan = next();
    else if (k == "--optics") a.optics = next();
    else if (k == "--events") { a.events = std::stoll(next()); a.lumi = -1.0; }
    else if (k == "--lumi") { a.lumi = std::stod(next()); a.events = -1; }
    else if (k == "--f0") a.sc.f0 = std::stod(next());
    else if (k == "--slope") a.sc.slope_b = std::stod(next());
    else if (k == "--amp") a.sc.amp = std::stod(next());
    else if (k == "--eps-b0") a.sc.eps_b0 = std::stod(next());
    else if (k == "--tmax") a.t_max = std::stod(next());
    else if (k == "--weighted-azimuth") a.weighted_azimuth = true;
    else if (k == "--seed") a.seed = std::stoull(next());
    else if (k == "--threads") a.threads = static_cast<unsigned>(std::stoul(next()));
    else if (k == "--out") a.out = next();
    else if (k == "-h" || k == "--help") usage(argv[0], "help");
    else usage(argv[0], "unknown option " + k);
  }
  if (a.config < 0 || a.config > 2) usage(argv[0], "--config must be 0, 1 or 2");
  if (a.threads < 1) usage(argv[0], "--threads must be >= 1");
  return a;
}

OpticsChoice optics_of(const std::string& s, const char* argv0) {
  if (s == "high-acceptance") return OpticsChoice::YellowReportHighAcceptance;
  if (s == "high-divergence") return OpticsChoice::YellowReportHighDivergence;
  if (s == "tagging") return OpticsChoice::Tagging;
  if (s == "tagging-legacy") return OpticsChoice::TaggingLegacyLevers;
  usage(argv0, "unknown --optics " + s);
}

RunPlan plan_of(const std::string& name, const Scenario& sc, const char* argv0) {
  if (name == "azz") return tensor_thirds_plan(sc.pol_ion_vector, sc.pol_ion_tensor);
  if (name == "cos2phi") return transverse_tensor_plan(sc.pol_ion_tensor);
  if (name == "flip") return tensor_flip_plan(sc.pol_ion_tensor);
  usage(argv0, "unknown --plan " + name);
}

}  // namespace

int main(int argc, char** argv) {
  const Args args = parse(argc, argv);
  try {
    PipelineConfig cfg;
    cfg.channel = PipelineChannel::CoherentLi6;
    cfg.isotope = "6Li";
    cfg.beam_config = args.config;
    cfg.optics_choice = optics_of(args.optics, argv[0]);
    cfg.seed = args.seed;
    cfg.coherent = args.sc;
    cfg.coherent_t_max = args.t_max;
    cfg.coherent_weighted_azimuth = args.weighted_azimuth;
    if (args.lumi > 0.0) cfg.lumi_pb = args.lumi;
    else cfg.n_events = static_cast<std::uint64_t>(args.events);

    const RunPlan plan = plan_of(args.plan, cfg.scenario, argv[0]);
    const Pipeline p(cfg, plan);
    const BeamConfig& bc = p.beam_config();
    const double p_u = bc.ion_momentum_per_nucleon;

    std::printf("LiPolGen coherent T0 generator -- e + 6Li -> e' + X + 6Li(g.s.)"
                "  (SCENARIO inputs)\n");
    std::printf("  beams        %s (sqrt(s_eN) = %.2f GeV), pot config %s\n",
                bc.label().c_str(), bc.sqrt_s_per_nucleon(),
                p.pot_config().c_str());
    std::printf("  scenario     f0 = %.4g (band 0.02-0.08), x_coh = %.4g,"
                " B = %.1f GeV^-2 (band 40-60)\n",
                args.sc.f0, args.sc.x_coh, args.sc.slope_b);
    std::printf("               flat cos2phi amp = %.3g at P_zz = 1"
                " (band 3e-3 - 1e-2), eps_b0 = %+.3f\n",
                args.sc.amp, args.sc.eps_b0);
    std::printf("               |t| <= %.2f GeV^2, <|t|> = 1/B = %.4f GeV^2,"
                " R_rms(B) = %.3f fm\n",
                args.t_max, 1.0 / args.sc.slope_b,
                std::sqrt(3.0 * args.sc.slope_b) * GEV_PER_FM_INV);
    std::printf("  optics       %s: %.4f x %.4f mrad envelope, L/L_HA = %.4f\n",
                p.optics().name.c_str(), 1e3 * p.optics().envelope_x(),
                1e3 * p.optics().envelope_y(), p.optics().lumi_fraction);
    std::printf("  window       Q2 >= %.2f, y in [%.3f, %.3f], W2 >= %.1f;"
                " %zu accepted (x, Q2) cells\n",
                cfg.scenario.q2_min, cfg.scenario.y_min, cfg.scenario.y_max,
                cfg.scenario.w2_min, p.dis_sampler().n_cells());

    std::printf("\n  %-10s %10s %10s %16s %14s\n", "category", "lumi_frac",
                "P_zz", "sigma [pb]", "events");
    for (std::size_t k = 0; k < plan.categories().size(); ++k) {
      const SpinCategory& c = plan.categories()[k];
      std::printf("  %-10s %10.5f %10.4f %16.6g %14llu\n", c.name.c_str(),
                  c.lumi_fraction, c.moments().tensor,
                  p.sigma_per_category_pb()[k],
                  static_cast<unsigned long long>(p.counts()[k]));
    }
    std::printf("  %-10s %10.5f %10s %16.6g %14llu\n", "TOTAL", 1.0, "",
                p.sigma_pb(), static_cast<unsigned long long>(p.size()));

    struct Menu { const char* label; OpticsChoice choice; };
    const Menu menu[] = {
        {"YR high-acceptance", OpticsChoice::YellowReportHighAcceptance},
        {"YR high-divergence", OpticsChoice::YellowReportHighDivergence},
        {"tagging", OpticsChoice::Tagging},
        {"tagging (18x275 levers)", OpticsChoice::TaggingLegacyLevers}};
    std::vector<Optics> envelopes;
    std::vector<std::uint64_t> tag(4, 0);
    for (const Menu& m : menu) {
      envelopes.push_back(optics_for(m.choice, "6Li", p_u));
    }

    std::unique_ptr<HepMC3Writer> writer;
    if (!args.out.empty()) writer = std::make_unique<HepMC3Writer>(args.out);

    const int nt = 10;
    std::vector<double> hist(nt, 0.0);
    double sum_t = 0.0, sum_w = 0.0, sum_wc = 0.0, worst_p = 0.0, worst_q = 0.0;
    std::uint64_t written = 0;
    const auto sink = [&](const Event& ev) {
      const Vec4 r = momentum_residual(ev);
      const Vec4 s = momentum_scale(ev);
      worst_p = std::max(worst_p,
                         std::max(std::fabs(r.e), std::fabs(r.pz)) / s.e);
      worst_q = std::max(worst_q, std::fabs(charge_residual(ev)));
      sum_t += ev.kin.t;
      sum_w += ev.weight;
      // <w cos 2(phi_t - phi_S)>: with the azimuth drawn flat and the
      // modulation carried as a weight this estimates c2/2 directly
      const Particle* rec = ev.find(Role::IntactRecoil);
      const double phi_t = std::atan2(rec->p.py, rec->p.px);
      sum_wc += ev.weight * std::cos(2.0 * (phi_t - ev.spin.phi_s));
      const int b = static_cast<int>(ev.kin.t / args.t_max * nt);
      if (b >= 0 && b < nt) hist[static_cast<std::size_t>(b)] += 1.0;
      for (std::size_t i = 0; i < envelopes.size(); ++i) {
        tag[i] += rp_tagged(ev, envelopes[i], p.pot_config());
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

    const double bw = args.t_max / nt;
    const double norm = 1.0 - std::exp(-args.sc.slope_b * args.t_max);
    std::printf("\n  %-16s %12s %12s\n", "|t| bin [GeV^2]", "fraction",
                "B exp(-B|t|)");
    for (int b = 0; b < nt; ++b) {
      const double lo = b * bw, hi = lo + bw;
      const double want = (std::exp(-args.sc.slope_b * lo) -
                           std::exp(-args.sc.slope_b * hi)) / norm;
      if (want < 1e-5 && hist[static_cast<std::size_t>(b)] == 0.0) continue;
      char lbl[32];
      std::snprintf(lbl, sizeof(lbl), "%.3f-%.3f", lo, hi);
      std::printf("  %-16s %12.5f %12.5f\n", lbl,
                  n ? hist[static_cast<std::size_t>(b)] / n : 0.0, want);
    }

    std::printf("\n  %-26s %16s %12s %12s\n", "optics", "envelope [mrad]",
                "near-beam tag", "analytic");
    for (std::size_t i = 0; i < envelopes.size(); ++i) {
      const double f = n ? static_cast<double>(tag[i]) / n : 0.0;
      const double an = args.sc.tag_acceptance_angular(
          envelopes[i].sigma_theta, p_u, 6, envelopes[i].n_sigma);
      char env[40];
      std::snprintf(env, sizeof(env), "%.4fx%.4f", 1e3 * envelopes[i].envelope_x(),
                    1e3 * envelopes[i].envelope_y());
      std::printf("  %-26s %16s %12.5f %12.5f\n", menu[i].label, env, f, an);
    }
    std::printf("  (the analytic column is exp(-B pT_cut^2), the acceptance of"
                " the CIRCLE inscribed at\n   n_sigma sigma_h.  The pots are"
                " planar and the real cut is the RECTANGLE that circle is"
                "\n   inscribed in, so the measured tag is the smaller"
                " number.)\n");

    std::printf("\n  generated    %llu events (seed %llu, %u thread%s)\n",
                static_cast<unsigned long long>(n),
                static_cast<unsigned long long>(args.seed), args.threads,
                args.threads == 1 ? "" : "s");
    std::printf("  time         %.3f s -> %.4g events/s%s\n", secs,
                secs > 0.0 ? static_cast<double>(n) / secs : 0.0,
                writer ? " (including HepMC3 output)" : "");
    std::printf("  <|t|>        %.5f GeV^2 (truncated exponential %.5f)\n",
                n ? sum_t / n : 0.0,
                1.0 / args.sc.slope_b -
                    args.t_max * std::exp(-args.sc.slope_b * args.t_max) / norm);
    std::printf("  <w>          %.6f, 2<w cos2(phi_t - phi_S)> = %.3e"
                " (the ensemble cos 2phi coefficient)\n",
                n ? sum_w / n : 0.0, n ? 2.0 * sum_wc / n : 0.0);
    std::printf("  worst |k + P_ion - k' - P_recoil - X| / E_in = %.3g,"
                " |dQ| = %.3g\n", worst_p, worst_q);
    if (writer) {
      std::printf("  wrote        %llu events to %s\n",
                  static_cast<unsigned long long>(written), args.out.c_str());
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "generate_coherent: " << e.what() << std::endl;
    return 1;
  }
}
