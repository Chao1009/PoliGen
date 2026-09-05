// SPDX-License-Identifier: GPL-3.0-or-later
// T0 event assembly: beams, the scattered-electron four-vector and its
// azimuth convention, the virtual photon, exact per-nucleon four-momentum
// conservation through the hadronic X, the spin labels, and the
// Fermi-smeared-target hook P4/P6 will supply.

#include <cmath>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "check_close.hpp"
#include "lipolgen/beams.hpp"
#include "lipolgen/bookkeeping.hpp"
#include "lipolgen/generator.hpp"
#include "lipolgen/sampler.hpp"
#include "lipolgen/sf.hpp"
#include "lipolgen/xsec.hpp"

using namespace lipolgen;

namespace {

std::shared_ptr<const InclusiveSampler> gen_sampler() {
  static std::shared_ptr<const InclusiveSampler> s = [] {
    InclusiveKernel::Options opt;
    opt.b1_func = [](double x, double q2, double f1) {
      return toy_b1(x, q2, f1);
    };
    opt.delta_func = [](double, double, double f1) { return -1e-2 * f1; };
    InclusiveSampler::GridSpec g;
    g.nx = 30; g.nq2 = 22;
    g.x_min = 1e-4; g.x_max = 1.0;
    g.q2_min = 0.7; g.q2_max = 2e3;
    return std::make_shared<InclusiveSampler>(
        std::make_shared<InclusiveKernel>(LI6(), opt), default_configs("6Li")[1],
        generator_scenario(Scenario()), g);
  }();
  return s;
}

double mdot(const Vec4& a, const Vec4& b) {
  return a.e * b.e - a.px * b.px - a.py * b.py - a.pz * b.pz;
}

double max_component(const Vec4& v) {
  return std::max(std::max(std::fabs(v.e), std::fabs(v.px)),
                  std::max(std::fabs(v.py), std::fabs(v.pz)));
}

}  // namespace

TEST_CASE("generator: the head-on beams are status-4 and carry the ion code") {
  const auto sampler = gen_sampler();
  const InclusiveGenerator gen(sampler);
  const BeamConfig& bc = sampler->config();
  CHECK(gen.ion_pdg() == 1000030060);          // 10LZZZAAAI, 6Li
  CHECK(nuclear_pdg(3, 7) == 1000030070);
  CHECK(nuclear_pdg(1, 2) == 1000010020);

  // massless electron along -z, whole-nucleus ion along +z
  CHECK(gen.beam_electron().e == bc.electron_energy);
  CHECK(gen.beam_electron().pz == -bc.electron_energy);
  CHECK(gen.beam_electron().px == 0.0);
  const double p_a = bc.ion.A * bc.ion_momentum_per_nucleon;
  CHECK_CLOSE(gen.beam_ion().pz, p_a, 1e-15);
  // m^2 = E^2 - p^2 at p/m ~ 100 loses ~4 digits to cancellation
  CHECK_CLOSE(std::sqrt(gen.beam_ion().m2()), bc.ion.mass(), 1e-9);
  // s per nucleon is the massless 4 E_e p_u the kernel was normalized with
  CHECK_CLOSE(sampler->s(), 4.0 * bc.electron_energy *
                                bc.ion_momentum_per_nucleon, 1e-14);

  const RunPlan plan = tensor_flip_plan(0.6);
  const auto cplan = sampler->make_plan(plan.categories()[0]);
  Rng rng(1, 1, 1, 1);
  const EventDraw d = sampler->draw_event(cplan, rng);
  const Event ev = gen.make_event(plan.categories()[0], plan, d, rng, 0, 1, 0);
  REQUIRE(ev.particles.size() == 6);
  CHECK(ev.particles[0].status == Status::Beam);
  CHECK(ev.particles[1].status == Status::Beam);
  CHECK(ev.particles[0].pdg == 11);
  CHECK(ev.particles[1].pdg == 1000030060);
  CHECK(ev.particles[1].charge == 3.0);
  CHECK(ev.find(Role::ScatteredElectron) != nullptr);
  CHECK(ev.find(Role::VirtualPhoton) != nullptr);
  CHECK(ev.find(Role::StruckNucleon) != nullptr);
  CHECK(ev.find(Role::HadronicX) != nullptr);
  CHECK(ev.channel == Channel::Inclusive);
}

TEST_CASE("generator: e' reproduces (x, Q2, y) and the phi convention") {
  const auto sampler = gen_sampler();
  const InclusiveGenerator gen(sampler);
  const BeamConfig& bc = sampler->config();
  const RunPlan plan = tensor_flip_plan(0.6);
  const SpinCategory& cat = plan.categories()[0];
  const auto cplan = sampler->make_plan(cat);

  for (std::uint64_t i = 0; i < 3000; ++i) {
    Rng rng(77, 0, 0, i);
    const EventDraw d = sampler->draw_event(cplan, rng);
    const Event ev = gen.make_event(cat, plan, d, rng, i, 1, 0);
    const Vec4& k = ev.particles[0].p;
    const Vec4& kp = ev.find(Role::ScatteredElectron)->p;
    const Vec4 q = k - kp;

    // Q2 = -q^2 exactly (the electron is built from Q2 by construction)
    CHECK_CLOSE(-mdot(q, q), d.q2, 1e-9);
    // y from the electron method: y = 1 - E'(1 - cos theta)/(2 E_e)
    const double ct = kp.pz / kp.e;
    CHECK_CLOSE(1.0 - kp.e * (1.0 - ct) / (2.0 * bc.electron_energy), d.y,
                1e-9);
    // x = Q2/(s y) with the massless per-nucleon s the kernel uses
    CHECK_CLOSE(d.q2 / (sampler->s() * d.y), d.x, 1e-12);
    // the electron is massless and goes backwards (theta ~ pi)
    CHECK_CLOSE_AT(kp.m2(), 0.0, 0.0, 1e-9 * kp.e * kp.e);
    // the electron lands inside the window the sampler accepted it in
    // (theta is measured from +z, the ION direction, so eta < 0 is the
    // backward/electron side -- but the generator window reaches eta = +3.8,
    // where a high-x electron does come out forward)
    CHECK(kp.e >= sampler->scenario().e_prime_min);
    const double eta = -std::log(std::tan(0.5 * std::acos(kp.pz / kp.e)));
    CHECK(eta >= sampler->scenario().eta_min - 1e-9);
    CHECK(eta <= sampler->scenario().eta_max + 1e-9);

    // THE azimuth convention: phi is atan2(py, px) of the SCATTERED ELECTRON
    // about the ion (+z) axis in the head-on frame, in [0, 2 pi).
    double phi_e = std::atan2(kp.py, kp.px);
    if (phi_e < 0.0) phi_e += 2.0 * kPi;
    CHECK_CLOSE_AT(phi_e, ev.kin.phi, 0.0, 1e-9);
    // the virtual photon shares the transverse direction of -kp
    const Vec4& g = ev.find(Role::VirtualPhoton)->p;
    CHECK_CLOSE_AT(std::atan2(g.py, g.px),
                   std::atan2(-kp.py, -kp.px), 0.0, 1e-12);
    CHECK(g.m2() < 0.0);   // spacelike
    CHECK_CLOSE(-g.m2(), d.q2, 1e-9);

    CHECK_CLOSE(ev.kin.w2, w2_from_xq2(d.x, d.q2), 1e-15);
    CHECK_CLOSE(ev.kin.nu, d.q2 / (2.0 * M_NUCLEON * d.x), 1e-15);
  }
}

TEST_CASE("generator: four-momentum is conserved on the per-nucleon target") {
  const auto sampler = gen_sampler();
  const InclusiveGenerator gen(sampler);
  const RunPlan plan = tensor_thirds_plan(0.7, 0.6);
  const SpinCategory& cat = plan.categories()[0];
  const auto cplan = sampler->make_plan(cat);
  for (std::uint64_t i = 0; i < 2000; ++i) {
    Rng rng(31, 0, 0, i);
    const EventDraw d = sampler->draw_event(cplan, rng);
    const Event ev = gen.make_event(cat, plan, d, rng, i, 1, 0);
    const Vec4& k = ev.particles[0].p;
    const Vec4& kp = ev.find(Role::ScatteredElectron)->p;
    const Vec4& pn = ev.find(Role::StruckNucleon)->p;
    const Vec4& xp = ev.find(Role::HadronicX)->p;
    const Vec4 residual = (k + pn) - (kp + xp);
    CHECK_CLOSE_AT(max_component(residual), 0.0, 0.0, 1e-9);
    // the default target is the on-shell nucleon moving with the beam
    CHECK_CLOSE(pn.pz, sampler->config().ion_momentum_per_nucleon, 1e-15);
    CHECK_CLOSE(std::sqrt(pn.m2()), M_NUCLEON, 1e-9);
    // X is a timelike hadronic system above the pion threshold
    CHECK(xp.m2() > 0.0);
    CHECK(std::sqrt(xp.m2()) > 1.0);
    // charge: only the struck nucleon's charge flows into X
    CHECK(ev.find(Role::HadronicX)->charge ==
          ev.find(Role::StruckNucleon)->charge);
    CHECK((ev.find(Role::StruckNucleon)->pdg == 2212 ||
           ev.find(Role::StruckNucleon)->pdg == 2112));
  }
}

TEST_CASE("generator: the target hook is what P4/P6 replace") {
  // Whatever the hook returns is used verbatim, and X follows it exactly --
  // that is the whole contract a Fermi-smeared nucleon needs.
  const auto sampler = gen_sampler();
  GeneratorConfig cfg;
  cfg.struck_nucleon_pdg = 2112;
  cfg.target = [](const EventDraw& d, Rng& rng) {
    // a deliberately silly "Fermi" kick, deterministic in the event stream
    const double kx = 0.2 * (rng.uniform() - 0.5);
    const double ky = 0.2 * (rng.uniform() - 0.5);
    const double pz = 0.9 * d.x + 90.0;
    return Vec4{std::sqrt(kx * kx + ky * ky + pz * pz +
                          M_NUCLEON * M_NUCLEON),
                kx, ky, pz};
  };
  const InclusiveGenerator gen(sampler, cfg);
  const RunPlan plan = transverse_tensor_plan(0.6);
  const SpinCategory& cat = plan.categories()[0];
  const auto cplan = sampler->make_plan(cat);
  for (std::uint64_t i = 0; i < 500; ++i) {
    Rng rng(4242, 0, 0, i);
    const EventDraw d = sampler->draw_event(cplan, rng);
    const Event ev = gen.make_event(cat, plan, d, rng, i, 1, 0);
    const Vec4& pn = ev.find(Role::StruckNucleon)->p;
    CHECK(std::fabs(pn.px) > 0.0);
    CHECK_CLOSE(pn.pz, 0.9 * d.x + 90.0, 1e-15);
    CHECK(ev.find(Role::StruckNucleon)->pdg == 2112);
    CHECK(ev.find(Role::HadronicX)->charge == 0.0);
    const Vec4 residual = (ev.particles[0].p + pn) -
                          (ev.find(Role::ScatteredElectron)->p +
                           ev.find(Role::HadronicX)->p);
    CHECK_CLOSE_AT(max_component(residual), 0.0, 0.0, 1e-9);
  }
}

TEST_CASE("generator: spin labels carry the run plan and the fill state") {
  const auto sampler = gen_sampler();
  const InclusiveGenerator gen(sampler);
  HelicityFlipOptions opt;
  opt.theta_s = 0.9;
  opt.phi_s = 0.4;
  const RunPlan plan = helicity_flip_plan(1.0, 0.7, 0.7, opt);
  const SpinCategory& cat = plan.categories()[1];   // lam_e = -1
  const auto cplan = sampler->make_plan(cat);
  Rng rng(8, 0, 0, 0);
  const EventDraw d = sampler->draw_event(cplan, rng);
  const Event ev = gen.make_event(cat, plan, d, rng, 17, 3, 1);
  CHECK(ev.number == 17u);
  CHECK(ev.spin.j == 1.0);
  CHECK(ev.spin.m_ion == d.m);
  CHECK(std::isnan(ev.spin.m_struck));    // inclusive: no struck cluster
  CHECK(ev.spin.lam_e == -1);
  CHECK(ev.spin.pe == 0.7);
  CHECK(ev.spin.theta_s == 0.9);
  CHECK(ev.spin.phi_s == 0.4);
  CHECK(ev.spin.pz == plan.pz_true());
  CHECK(ev.spin.pzz == plan.pzz_true());
  CHECK(ev.spin.category == "apar-");
  CHECK(ev.spin.run == 3);
  CHECK(ev.spin.bunch == 1);
  CHECK(ev.particles[0].pol == -1.0);           // SPINUP-style helicity label
  CHECK(ev.particles[1].pol == d.m);            // ion projection
  CHECK(ev.find(Role::ScatteredElectron)->pol == -1.0);
}

TEST_CASE("generator: run() streams a whole plan without storing it") {
  const auto sampler = gen_sampler();
  const InclusiveGenerator gen(sampler);
  const RunPlan plan = tensor_thirds_plan(0.7, 0.6, 1e-4);
  const double lumi = 4.0;   // pb^-1, a few thousand events

  std::map<std::string, std::size_t> counted;
  std::map<std::string, double> xsec_seen;
  std::uint64_t emitted = 0;
  const std::uint64_t total = gen.run(plan, lumi, 20260713,
                                      [&](const Event& ev) {
                                        counted[ev.spin.category] += 1;
                                        xsec_seen[ev.spin.category] = ev.xsec_pb;
                                        ++emitted;
                                      });
  CHECK(total == emitted);
  CHECK(total > 0u);
  CHECK(counted.size() == 3u);
  const std::vector<double> sigma = gen.sigma_per_category(plan);
  for (std::size_t k = 0; k < plan.categories().size(); ++k) {
    const std::string& name = plan.categories()[k].name;
    CHECK_CLOSE(xsec_seen[name], sigma[k], 1e-15);
    const double mu = lumi * plan.categories()[k].lumi_fraction * sigma[k];
    CHECK_CLOSE_AT(static_cast<double>(counted[name]), mu, 0.0,
                   6.0 * std::sqrt(mu));
  }
  // the same (seed, run) is the same run
  std::vector<double> a, b;
  gen.run(plan, lumi, 20260713, [&](const Event& ev) { a.push_back(ev.kin.x); });
  gen.run(plan, lumi, 20260713, [&](const Event& ev) { b.push_back(ev.kin.x); });
  CHECK(a == b);
  CHECK(a.size() == total);

  // fixed statistics: the split is exact and proportional to share * sigma
  std::map<std::string, std::size_t> n_fixed;
  const std::uint64_t nn = gen.run_n(plan, 5000, 99, [&](const Event& ev) {
    n_fixed[ev.spin.category] += 1;
  });
  CHECK(nn == 5000u);
  std::size_t sum = 0;
  for (const auto& kv : n_fixed) sum += kv.second;
  CHECK(sum == 5000u);
  CHECK_THROWS_AS(gen.run(plan, lumi, 1, nullptr), std::runtime_error);
}

// P1.  The struck nucleon used to be drawn flat Z : N, but the inclusive rate
// off a nucleus is Z F2p(x, Q2) + N F2n(x, Q2).  For an N = Z nucleus like 6Li
// the two agree only where F2n = F2p, i.e. nowhere: at x = 0.5 the true
// proton share is 0.616, not 0.5, and it keeps rising towards the valence edge.
TEST_CASE("generator: the struck nucleon is drawn Z F2p : N F2n") {
  const auto sampler = gen_sampler();
  const InclusiveGenerator gen(sampler);
  const Ion& ion = sampler->config().ion;
  const UnpolSF& sf = *sampler->kernel().nuclear_f2().base();

  // --- the analytic fraction, straight off the kernel's own backend
  for (double x : {0.01, 0.1, 0.3, 0.5, 0.7}) {
    const double q2 = 10.0;
    const double want = ion.Z * sf.f2p(x, q2)
                        / (ion.Z * sf.f2p(x, q2) + ion.N() * sf.f2n(x, q2));
    CHECK_CLOSE(gen.proton_fraction(x, q2), want, 1e-14);
    CHECK(want > 0.0);
    CHECK(want < 1.0);
  }
  // it is NOT Z/A, and the departure grows with x
  const double flat = static_cast<double>(ion.Z) / ion.A;
  CHECK(flat == 0.5);
  const double f05 = gen.proton_fraction(0.5, 10.0);
  const double f07 = gen.proton_fraction(0.7, 10.0);
  MESSAGE("6Li P(proton): Z/A = " << flat << ", Z F2p : N F2n = "
          << gen.proton_fraction(0.01, 10.0) << " at x = 0.01, " << f05
          << " at x = 0.5, " << f07 << " at x = 0.7");
  CHECK(f05 > flat + 0.05);
  CHECK(f07 > f05);

  // --- and the DRAW follows it.  make_event at a fixed (x, Q2) over many
  // streams: the only randomness left in the species is the one uniform.
  const RunPlan plan = tensor_flip_plan(0.6);
  const SpinCategory& cat = plan.categories()[0];
  for (double x : {0.05, 0.5}) {
    EventDraw d;
    d.x = x;
    d.q2 = 10.0;
    d.y = d.q2 / (sampler->s() * d.x);
    d.phi = 0.7;
    d.m = 0.0;
    const int n = 40000;
    int n_p = 0;
    for (int i = 0; i < n; ++i) {
      Rng rng(7, 0, 0, static_cast<std::uint64_t>(i));
      const Event ev = gen.make_event(cat, plan, d, rng, i, 1, 0);
      const Particle* sn = ev.find(Role::StruckNucleon);
      REQUIRE(sn != nullptr);
      if (sn->pdg == 2212) ++n_p;
      CHECK(sn->charge == (sn->pdg == 2212 ? 1.0 : 0.0));
      // P3: on shell at the FREE nucleon mass in both places it is built
      CHECK_CLOSE(sn->p.m2(), M_NUCLEON * M_NUCLEON, 1e-9);
    }
    const double got = static_cast<double>(n_p) / n;
    const double want = gen.proton_fraction(x, d.q2);
    const double err = std::sqrt(want * (1.0 - want) / n);
    MESSAGE("x = " << x << ": drew p in " << got << " +- " << err
            << " of events, Z F2p/(Z F2p + N F2n) = " << want);
    CHECK(std::fabs(got - want) < 4.0 * err);
  }
  // the pinned-species option still wins
  GeneratorConfig gc;
  gc.struck_nucleon_pdg = 2112;
  const InclusiveGenerator pinned(sampler, gc);
  EventDraw d;
  d.x = 0.5; d.q2 = 10.0; d.y = d.q2 / (sampler->s() * d.x); d.phi = 0.0;
  Rng rng(1, 0, 0, 0);
  CHECK(pinned.make_event(cat, plan, d, rng, 0, 1, 0)
            .find(Role::StruckNucleon)->pdg == 2112);
}
