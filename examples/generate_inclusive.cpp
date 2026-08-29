// Generates T0 inclusive e + 6Li / 7Li events for a standard run plan and
// writes them as HepMC3 Asciiv3 through lipolgen::HepMC3Writer, printing the
// accepted cross section per spin category and the generation throughput.
//
// Usage:
//   generate_inclusive [options]
//     --isotope 6Li|7Li     ion species                       (default 6Li)
//     --config 0|1|2        energy point low/mid/top          (default 1)
//     --plan NAME           apar | azz | cos2phi | flip       (default azz)
//     --events N            fixed event count
//     --lumi L              integrated luminosity [pb^-1]     (default 2)
//     --seed S              RNG seed                          (default 20260713)
//     --out FILE            HepMC3 output ("" = none)
//                                       (default generate_inclusive.hepmc)
//     --analysis-window     use the ANALYSIS cuts instead of the looser
//                           generator window
//
// --events and --lumi are exclusive; --events splits the count across the
// plan's categories in proportion to lumi_fraction * sigma.
//
// Verify the output externally with pyhepmc, e.g.:
//   python3 -c "
//   import pyhepmc
//   with pyhepmc.open('generate_inclusive.hepmc') as f:
//       ev = next(iter(f))
//   print(len(ev.particles), ev.attributes['spin_category'])
//   "
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "lipolgen/beams.hpp"
#include "lipolgen/bookkeeping.hpp"
#include "lipolgen/generator.hpp"
#include "lipolgen/hepmc_writer.hpp"
#include "lipolgen/sampler.hpp"
#include "lipolgen/sf.hpp"
#include "lipolgen/xsec.hpp"

using namespace lipolgen;

namespace {

struct Args {
  std::string isotope = "6Li";
  int config = 1;
  std::string plan = "azz";
  long long events = -1;
  double lumi = 2.0;
  unsigned long long seed = 20260713;
  std::string out = "generate_inclusive.hepmc";
  bool analysis_window = false;
};

[[noreturn]] void usage(const char* argv0, const std::string& why) {
  std::cerr << argv0 << ": " << why << "\n"
            << "usage: " << argv0
            << " [--isotope 6Li|7Li] [--config 0|1|2]"
               " [--plan apar|azz|cos2phi|flip]\n"
               "       [--events N | --lumi PB] [--seed S] [--out FILE]"
               " [--analysis-window]\n";
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
    if (k == "--isotope") a.isotope = next();
    else if (k == "--config") a.config = std::stoi(next());
    else if (k == "--plan") a.plan = next();
    else if (k == "--events") a.events = std::stoll(next());
    else if (k == "--lumi") a.lumi = std::stod(next());
    else if (k == "--seed") a.seed = std::stoull(next());
    else if (k == "--out") a.out = next();
    else if (k == "--analysis-window") a.analysis_window = true;
    else if (k == "-h" || k == "--help") usage(argv[0], "help");
    else usage(argv[0], "unknown option " + k);
  }
  if (a.config < 0 || a.config > 2) usage(argv[0], "--config must be 0, 1 or 2");
  if (a.isotope != "6Li" && a.isotope != "7Li") {
    usage(argv[0], "--isotope must be 6Li or 7Li");
  }
  return a;
}

/// The kernel each isotope is generated with: 6Li carries the Miller
/// deuteron b1 through the rank-2 embedded-deuteron transfer plus a
/// discovery-scale Delta; 7Li has no published rank-2 input, so its tensor
/// slots stay empty (xsec.hpp's documented default = identically zero).
std::shared_ptr<const InclusiveKernel> make_kernel(const Ion& ion) {
  InclusiveKernel::Options opt;
  if (std::fabs(ion.spin - 1.0) < 1e-9) {
    static const auto li6_b1 =
        std::make_shared<Li6B1>(std::make_shared<MillerB1>());
    opt.b1_func = li6_b1->b1_func();
    opt.delta_func = [](double x, double q2, double f1) {
      return toy_delta_gluon(x, q2, f1, 1e-2);
    };
  }
  return std::make_shared<InclusiveKernel>(ion, opt);
}

RunPlan make_plan(const std::string& name, const Ion& ion,
                  const Scenario& sc) {
  const bool spin1 = std::fabs(ion.spin - 1.0) < 1e-9;
  if (name == "apar") {
    return helicity_flip_plan(ion.spin, sc.pol_ion_vector, sc.pol_electron);
  }
  if (!spin1) {
    throw std::runtime_error("run plan '" + name +
                             "' is spin-1 only; " + ion.name +
                             " has spin " + std::to_string(ion.spin) +
                             " (use --plan apar)");
  }
  if (name == "azz") {
    return tensor_thirds_plan(sc.pol_ion_vector, sc.pol_ion_tensor);
  }
  if (name == "cos2phi") return transverse_tensor_plan(sc.pol_ion_tensor);
  if (name == "flip") return tensor_flip_plan(sc.pol_ion_tensor);
  throw std::runtime_error("unknown run plan '" + name + "'");
}

}  // namespace

int main(int argc, char** argv) {
  const Args args = parse(argc, argv);
  try {
    const Ion& ion = ion_by_name(args.isotope);
    const BeamConfig cfg =
        default_configs(args.isotope)[static_cast<std::size_t>(args.config)];
    const Scenario analysis;
    const Scenario scenario =
        args.analysis_window ? analysis : generator_scenario(analysis);

    const auto sampler = std::make_shared<InclusiveSampler>(
        make_kernel(ion), cfg, scenario, InclusiveSampler::GridSpec());
    const InclusiveGenerator gen(sampler);
    const RunPlan plan = make_plan(args.plan, ion, scenario);

    std::printf("LiPolGen inclusive T0 generator\n");
    std::printf("  beams        %s, %.1f GeV e- x %.1f GeV/u %s"
                " (sqrt(s_eN) = %.2f GeV)\n",
                cfg.label().c_str(), cfg.electron_energy,
                cfg.ion_momentum_per_nucleon, ion.name.c_str(),
                cfg.sqrt_s_per_nucleon());
    std::printf("  window       Q2 >= %.2f, y in [%.3f, %.3f], W2 >= %.1f,"
                " eta in [%.1f, %.1f], E' >= %.2f  (%s)\n",
                scenario.q2_min, scenario.y_min, scenario.y_max,
                scenario.w2_min, scenario.eta_min, scenario.eta_max,
                scenario.e_prime_min,
                args.analysis_window ? "analysis" : "generator");
    std::printf("  grid         %zu accepted (x, Q2) cells\n",
                sampler->n_cells());
    std::printf("  run plan     %s: P_e = %.3f, P_z = %.3f, P_zz = %.4f\n",
                args.plan.c_str(), plan.pe_true(), plan.pz_true(),
                plan.pzz_true());

    const std::vector<double> sigma = gen.sigma_per_category(plan);
    std::printf("\n  %-10s %8s %6s %14s %14s\n", "category", "lumi_frac",
                "lam_e", "sigma [pb]", "a2_eff");
    for (std::size_t k = 0; k < plan.categories().size(); ++k) {
      const SpinCategory& c = plan.categories()[k];
      const auto mod = sampler->effective_modulation(c);
      std::printf("  %-10s %8.5f %6d %14.6g %14.6g\n", c.name.c_str(),
                  c.lumi_fraction, c.lam_e, sigma[k], mod.a2);
    }

    std::unique_ptr<HepMC3Writer> writer;
    if (!args.out.empty()) writer = std::make_unique<HepMC3Writer>(args.out);

    std::uint64_t written = 0;
    Vec4 sum_e;   // running check that the per-nucleon balance holds
    double worst_residual = 0.0;
    const auto sink = [&](const Event& ev) {
      const Vec4 r = (ev.particles[0].p + ev.find(Role::StruckNucleon)->p) -
                     (ev.find(Role::ScatteredElectron)->p +
                      ev.find(Role::HadronicX)->p);
      worst_residual = std::max(worst_residual,
                                std::max(std::fabs(r.e), std::fabs(r.pz)));
      sum_e = sum_e + ev.find(Role::ScatteredElectron)->p;
      if (writer) writer->write(ev);
      ++written;
    };

    const auto t0 = std::chrono::steady_clock::now();
    const std::uint64_t n =
        args.events >= 0
            ? gen.run_n(plan, static_cast<std::uint64_t>(args.events),
                        args.seed, sink)
            : gen.run(plan, args.lumi, args.seed, sink);
    const auto t1 = std::chrono::steady_clock::now();
    if (writer) writer->close();

    const double secs = std::chrono::duration<double>(t1 - t0).count();
    std::printf("\n  generated    %llu events", static_cast<unsigned long long>(n));
    if (args.events < 0) std::printf(" at %g pb^-1", args.lumi);
    std::printf(" (seed %llu)\n", static_cast<unsigned long long>(args.seed));
    std::printf("  time         %.3f s -> %.3g events/s%s\n", secs,
                secs > 0.0 ? static_cast<double>(n) / secs : 0.0,
                writer ? " (including HepMC3 output)" : "");
    std::printf("  <E'>         %.4f GeV\n",
                n ? sum_e.e / static_cast<double>(n) : 0.0);
    std::printf("  worst |k + P_N - k' - X| = %.3g GeV\n", worst_residual);
    if (writer) std::printf("  wrote        %llu events to %s\n",
                            static_cast<unsigned long long>(written),
                            args.out.c_str());
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "generate_inclusive: " << e.what() << std::endl;
    return 1;
  }
}
