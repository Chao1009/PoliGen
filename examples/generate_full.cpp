// Generates T0+T2 events through the full chain -- lipolgen::Pipeline
// (inclusive | tagged 6Li-alpha | tagged 7Li-alpha | coherent 6Li), hadronized
// by lipolgen::PythiaBridge, written as HepMC3 Asciiv3 -- and prints the
// accepted cross section per spin category, the PYTHIA veto/retry/no-surrogate
// counters, the far-forward tag fraction, and the throughput.
//
// Usage:
//   generate_full [options]
//     --isotope 6Li|7Li     ion species (inclusive channel only; the tagged
//                           channels imply their own isotope, coherent is
//                           6Li only)                        (default 6Li)
//     --config 0|1|2        energy point low/mid/top          (default 1)
//     --channel NAME        inclusive | 6Li-alpha | 7Li-alpha | coherent
//                                                              (default
//                                                               inclusive)
//     --plan NAME           azz | apar | cos2phi | flip | pure
//                                                              (default azz)
//     --events N            fixed event count                 (default 20000)
//     --seed S              RNG seed                          (default
//                                                               20260713)
//     --out FILE             HepMC3 output ("" = none)
//                                       (default generate_full.hepmc)
//     --tier T0|T1           tagged final-state fidelity      (default T1:
//                            the struck cluster is resolved into a nucleon
//                            plus its partner spectator(s), breakup.hpp)
//     --no-hadronize         T0 only: PYTHIA is never linked into the run
//
// PYTHIA BINDING (docs/T2_CHAIN.md has the full story).  `PipelineConfig::
// hadronizer` is bound to `PythiaBridge::hadronize`, exactly the snippet
// docs/USAGE.md section 6 shows -- and, since the T1 tier landed, that is
// ALL a tagged channel needs.  The pipeline resolves the struck cluster into
// a struck NUCLEON plus its on-shell partner spectator(s) before the bridge
// sees the event (`PipelineConfig::tier`, default `Tier::T1`), so the bridge
// takes `Role::StruckNucleon` verbatim, the whole-record balance
// (spectator + partners + hadrons + e') closes to the numerical floor, and
// no `set_nucleon_in_cluster` hook is installed by this program any more.
// Run with `--tier T0` to see the old pseudo-cluster record instead -- the
// residual it prints below is then the honest ~19 % of the v0 fallback.
//
// COHERENT has no analogous fix available: a coherent event carries neither
// `Role::StruckNucleon` nor `Role::StruckCluster` (the diffractive system X
// is not a struck nucleon), so `hadronize()` always falls back to "a nucleon
// at rest in the ion frame" -- unrelated to the true (small, largely
// transverse) momentum transfer to the recoil.  The call still succeeds
// mechanically (no crash, no excess vetoes), but whole-record momentum and
// charge are NOT conserved with any binding available through the public
// PythiaBridge API.  This is reported here, not patched around: PythiaBridge
// v0 has no representation of a coherent-diffractive final state.
#include <algorithm>
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
#include "lipolgen/pythia_bridge.hpp"

using namespace lipolgen;

namespace {

struct Args {
  std::string isotope = "6Li";
  int config = 1;
  std::string channel = "inclusive";
  std::string plan = "azz";
  long long events = 20000;
  unsigned long long seed = 20260713;
  std::string out = "generate_full.hepmc";
  std::string tier = "T1";
  bool hadronize = true;
};

[[noreturn]] void usage(const char* argv0, const std::string& why) {
  std::cerr << argv0 << ": " << why << "\n"
            << "usage: " << argv0
            << " [--isotope 6Li|7Li] [--config 0|1|2]\n"
               "       [--channel inclusive|6Li-alpha|7Li-alpha|coherent]"
               " [--plan azz|apar|cos2phi|flip|pure]\n"
               "       [--events N] [--seed S] [--out FILE] [--tier T0|T1]"
               " [--no-hadronize]\n";
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
    else if (k == "--channel") a.channel = next();
    else if (k == "--plan") a.plan = next();
    else if (k == "--events") a.events = std::stoll(next());
    else if (k == "--seed") a.seed = std::stoull(next());
    else if (k == "--out") a.out = next();
    else if (k == "--tier") a.tier = next();
    else if (k == "--no-hadronize") a.hadronize = false;
    else if (k == "-h" || k == "--help") usage(argv[0], "help");
    else usage(argv[0], "unknown option " + k);
  }
  if (a.config < 0 || a.config > 2) usage(argv[0], "--config must be 0, 1 or 2");
  if (a.events <= 0) usage(argv[0], "--events must be positive");
  if (a.tier != "T0" && a.tier != "T1") usage(argv[0], "--tier must be T0 or T1");
  return a;
}

PipelineChannel channel_of(const std::string& s, const char* argv0) {
  if (s == "inclusive") return PipelineChannel::Inclusive;
  if (s == "6Li-alpha") return PipelineChannel::TaggedLi6Alpha;
  if (s == "7Li-alpha") return PipelineChannel::TaggedLi7Alpha;
  if (s == "coherent") return PipelineChannel::CoherentLi6;
  usage(argv0, "unknown --channel " + s);
}

RunPlan plan_of(const std::string& name, double j, const Scenario& sc,
                const char* argv0) {
  const bool spin1 = std::fabs(j - 1.0) < 1e-9;
  if (name == "apar") return helicity_flip_plan(j, sc.pol_ion_vector, sc.pol_electron);
  if (name == "pure") {
    std::vector<double> pops(static_cast<std::size_t>(2 * j + 1.5), 0.0);
    pops[0] = 1.0;  // stretched state M = +J
    // The stretched state has rank-2 moment T = +1 for ANY J >= 1 (7Li
    // |3/2,+3/2> included: <3Jz^2 - J(J+1)>/3 = (27/4 - 15/4)/3 = 1), exactly
    // as generate_tagged.cpp records; only a spin-1/2 fill has none.
    return RunPlan({SpinCategory("pure", j, pops, 0, 0.0, 0.0, 0.0, 1.0)}, 0.0,
                   1.0, j >= 1.0 ? 1.0 : 0.0);
  }
  if (!spin1) {
    usage(argv0, "--plan " + name + " is spin-1 only (6Li); 7Li needs "
                                    "--plan apar or --plan pure");
  }
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
    cfg.channel = channel_of(args.channel, argv[0]);
    const bool tagged = is_tagged(cfg.channel);
    const bool coherent = cfg.channel == PipelineChannel::CoherentLi6;
    cfg.isotope = tagged ? channel_isotope(cfg.channel)
                         : coherent ? std::string("6Li") : args.isotope;
    cfg.beam_config = args.config;
    cfg.optics_choice = tagged || coherent ? OpticsChoice::Tagging
                                           : OpticsChoice::YellowReportHighAcceptance;
    cfg.seed = args.seed;
    cfg.n_events = static_cast<std::uint64_t>(args.events);
    cfg.tier = (args.tier == "T0") ? Tier::T0 : Tier::T1;

    const Ion& ion = ion_by_name(cfg.isotope);
    const RunPlan plan = plan_of(args.plan, ion.spin, cfg.scenario, argv[0]);

    const BeamConfig bc = default_configs(cfg.isotope)[static_cast<std::size_t>(cfg.beam_config)];
    std::unique_ptr<PythiaBridge> bridge;
    PythiaBridgeStats snapshot;  // last stats() read (bridge dies before we print)
    if (args.hadronize) {
      bridge = std::make_unique<PythiaBridge>(bc);
      cfg.hadronizer = [&](Event& ev, Rng& rng) { bridge->hadronize(ev, rng); };
    }

    const auto t_build0 = std::chrono::steady_clock::now();
    const Pipeline p(cfg, plan);
    const auto t_build1 = std::chrono::steady_clock::now();

    std::printf("LiPolGen full T0+T2 chain -- %s, %s\n",
               pipeline_channel_name(cfg.channel), p.beam_config().label().c_str());
    std::printf("  beams        sqrt(s_eN) = %.2f GeV, pot config %s, optics %s\n",
               p.beam_config().sqrt_s_per_nucleon(), p.pot_config().c_str(),
               p.optics().name.c_str());
    std::printf("  hadronizer   %s\n",
               args.hadronize ? "PythiaBridge::hadronize (plain binding)"
                              : "none (--no-hadronize: T0 only)");
    if (tagged) {
      std::printf("  tier         %s\n",
                 p.tier() == Tier::T1
                     ? "T1 -- struck nucleon + partner spectators (breakup.hpp)"
                     : "T0 -- struck cluster as one pseudo-particle (does NOT "
                       "conserve through the bridge)");
    }
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

    std::unique_ptr<HepMC3Writer> writer;
    if (!args.out.empty()) writer = std::make_unique<HepMC3Writer>(args.out);

    std::uint64_t written = 0, n_tagged = 0, n_taggable = 0;
    double worst_p_rel = 0.0, worst_q = 0.0;
    const auto sink = [&](const Event& ev) {
      if (tagged || coherent) {
        ++n_taggable;
        n_tagged += rp_tagged(ev, p.optics(), p.pot_config());
      }
      if (args.hadronize) {
        // "beams in = Status::Final out", per-nucleon for inclusive
        // (pipeline.hpp: the (A-1) remnant is never written there),
        // whole-nucleus for tagged/coherent -- see docs/T2_CHAIN.md.
        const Particle* pn = ev.find(Role::StruckNucleon);
        const Vec4 want_p = (cfg.channel == PipelineChannel::Inclusive && pn)
                                 ? ev.particles[0].p + pn->p
                                 : ev.particles[0].p + ev.particles[1].p;
        const double want_q = (cfg.channel == PipelineChannel::Inclusive && pn)
                                   ? ev.particles[0].charge + pn->charge
                                   : ev.particles[0].charge + ev.particles[1].charge;
        const Vec4 r = want_p - ev.total_final();
        worst_p_rel = std::max(
            worst_p_rel, std::max({std::fabs(r.e), std::fabs(r.px),
                                   std::fabs(r.py), std::fabs(r.pz)}) /
                             want_p.e);
        worst_q = std::max(worst_q, std::fabs(want_q - ev.total_charge_final()));
      }
      if (writer) writer->write(ev);
      ++written;
    };

    const auto t0 = std::chrono::steady_clock::now();
    const std::uint64_t n = p.for_each(sink);
    const auto t1 = std::chrono::steady_clock::now();
    if (writer) writer->close();
    const double secs = std::chrono::duration<double>(t1 - t0).count();
    if (bridge) snapshot = bridge->stats();

    if (tagged || coherent) {
      std::printf("\n  far-forward tag (%s, this run's optics): %.4f"
                 " (%llu / %llu)\n",
                 coherent ? "intact recoil" : "spectator",
                 n_taggable ? double(n_tagged) / double(n_taggable) : 0.0,
                 static_cast<unsigned long long>(n_tagged),
                 static_cast<unsigned long long>(n_taggable));
    }

    if (args.hadronize) {
      std::printf("\n  PYTHIA        called %llu, ok %llu, failed %llu,"
                 " retries %llu, no-surrogate %llu, struck p/n %llu/%llu\n",
                 static_cast<unsigned long long>(snapshot.n_called),
                 static_cast<unsigned long long>(snapshot.n_ok),
                 static_cast<unsigned long long>(snapshot.n_failed),
                 static_cast<unsigned long long>(snapshot.n_retries),
                 static_cast<unsigned long long>(snapshot.n_no_surrogate),
                 static_cast<unsigned long long>(snapshot.n_proton),
                 static_cast<unsigned long long>(snapshot.n_neutron));
      std::printf("  conservation  worst |beams_in - Sum(Status::Final)| / E_in"
                 " = %.4g, worst |charge| = %.4g\n", worst_p_rel, worst_q);
      if (coherent) {
        std::printf("                NOT expected to be small: PythiaBridge v0"
                   " has no coherent-diffractive target (see the file header"
                   " and docs/T2_CHAIN.md).\n");
      } else if (worst_p_rel > 1e-6) {
        std::printf("                EXPECTED at --tier T0: the bridge's"
                   " DEPRECATED cluster branch feeds PYTHIA p_cluster/A_c"
                   " with nothing carrying the rest.  Tier::T1 (the default)"
                   " closes it; see docs/T2_CHAIN.md.\n");
      }
    }

    std::printf("\n  generated    %llu events (seed %llu)\n",
               static_cast<unsigned long long>(n),
               static_cast<unsigned long long>(args.seed));
    std::printf("  time         %.3f s -> %.4g events/s%s\n", secs,
               secs > 0.0 ? static_cast<double>(n) / secs : 0.0,
               writer ? " (including HepMC3 output)" : "");
    if (writer) {
      std::printf("  wrote        %llu events to %s\n",
                 static_cast<unsigned long long>(written), args.out.c_str());
    }
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "generate_full: " << e.what() << std::endl;
    return 1;
  }
}
