// P6 -- the PYTHIA 8 hadronization tier.  Every event here is built by hand:
// the core generator (P3) is a concurrent work package, so the test owns its
// own head-on-frame kinematics and pins the convention it uses.
//
// FRAME AND CONVENTION (docs/CONVENTIONS.md, mirrored by
// lipolgen::dis_scattered_electron):
//   ion along +z at p_u GeV per nucleon, electron along -z at E_e,
//   k   = (E_e, 0, 0, -E_e),   p_N = (E_N, 0, 0, +p_u),
//   Q^2 = -(k-k')^2,  y = p_N.(k-k')/p_N.k,  x = Q^2/(2 p_N.(k-k')),
//   phi = atan2(k'_y, k'_x), the azimuth of e' about +z.
// From (x, Q^2) and the beams, y = Q^2/(2 x p_N.k) and
//   E'   = E_e (1-y) + p_u Q^2 / (2 E_e (E_N + p_u)),
//   k'_z = Q^2/(2 E_e) - E',   k'_T = sqrt(E'^2 - k'_z^2).
//
// Beams: 10 x 99.5 GeV/u 6Li -- the mid gamma-matched EIC configuration of
// beams.default_configs("6Li"), the one the standing PYTHIA production of
// PolarizedLithiumSim is made at.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "check_close.hpp"
#include "doctest.h"

#include "lipolgen/beams.hpp"
#include "lipolgen/event.hpp"
#include "lipolgen/pythia_bridge.hpp"
#include "lipolgen/rng.hpp"

#include "Pythia8/Pythia.h"

using namespace lipolgen;

namespace {

constexpr double kEe = 10.0;
constexpr double kPu = 99.5;
constexpr std::uint64_t kSeed = 20260829;

BeamConfig mid_config() {
  BeamConfig bc;
  bc.electron_energy = kEe;
  bc.ion = LI6();
  bc.ion_momentum_per_nucleon = kPu;
  return bc;
}

/// A nucleon at rest in the ion rest frame: P_ion / A.
/// P3.  The struck nucleon of the per-nucleon subsystem is ON SHELL AT THE
/// FREE NUCLEON MASS -- `InclusiveGenerator::target_nucleon` and the bridge's
/// own implicit target both use M_NUCLEON (0.9383), never the ion's mass per
/// nucleon (0.9338), because every per-nucleon kinematic label on the record
/// (W2, nu) is built on M_NUCLEON.  See docs/CONVENTIONS.md.
Vec4 nucleon_at_rest_in_ion(const BeamConfig& bc) {
  (void)bc;
  const double m = M_NUCLEON;
  return {std::sqrt(kPu * kPu + m * m), 0.0, 0.0, kPu};
}

struct Point {
  double x, q2, phi;
};
// A few (x, Q^2, phi) inside the generator window (Q^2 >= 0.7,
// y in [0.004, 0.985], W^2 >= 8).
const std::vector<Point>& points() {
  static const std::vector<Point> p = {
      {0.005, 2.0, 0.0},    {0.010, 5.0, 1.1},   {0.030, 10.0, 2.4},
      {0.080, 20.0, 3.9},   {0.200, 40.0, 5.2},  {0.400, 80.0, 0.7},
  };
  return p;
}

/// Build the T0 part of an event: beams, the struck nucleon, e'.
/// Returns false when (x, Q^2) is outside the physical region.
bool make_event(const BeamConfig& bc, const Point& pt, const Vec4& p_n,
                int id_n, std::uint64_t number, Event* out) {
  Vec4 kp;
  double y = 0.0;
  if (!dis_scattered_electron(bc.electron_energy, p_n, pt.x, pt.q2, pt.phi, &kp,
                              &y))
    return false;

  Event ev;
  ev.number = number;
  ev.channel = Channel::Inclusive;
  ev.spin.lam_e = +1;
  ev.kin.x = pt.x;
  ev.kin.q2 = pt.q2;
  ev.kin.y = y;
  ev.kin.phi = pt.phi;

  Particle be;
  be.pdg = 11;
  be.status = Status::Beam;
  be.role = Role::BeamElectron;
  be.p = {bc.electron_energy, 0.0, 0.0, -bc.electron_energy};
  be.charge = -1.0;
  be.pol = +1.0;
  ev.particles.push_back(be);

  Particle bi;
  bi.pdg = 1000030060;
  bi.status = Status::Beam;
  bi.role = Role::BeamIon;
  bi.p = {std::sqrt(kPu * kPu * 36.0 + bc.ion.mass() * bc.ion.mass()), 0.0, 0.0,
          6.0 * kPu};
  bi.charge = 3.0;
  ev.particles.push_back(bi);

  Particle sn;
  sn.pdg = id_n;
  sn.status = Status::Intermediate;
  sn.role = Role::StruckNucleon;
  sn.p = p_n;
  sn.mass = std::sqrt(std::max(0.0, p_n.m2()));
  sn.charge = (id_n == 2212) ? 1.0 : 0.0;
  ev.particles.push_back(sn);

  Particle se;
  se.pdg = 11;
  se.status = Status::Final;
  se.role = Role::ScatteredElectron;
  se.p = kp;
  se.charge = -1.0;
  ev.particles.push_back(se);

  Particle hx;   // the T0 pseudo-particle the bridge must demote
  hx.pdg = 0;
  hx.status = Status::Final;
  hx.role = Role::HadronicX;
  hx.p = p_n + (be.p - kp);
  ev.particles.push_back(hx);

  *out = ev;
  return true;
}

/// Sum of the hadrons plus the scattered electron, i.e. what must equal
/// k + p_N.
Vec4 hadrons_plus_electron(const Event& ev) {
  Vec4 s = ev.find(Role::ScatteredElectron)->p;
  for (const auto& p : ev.particles)
    if (p.role == Role::Hadron) s = s + p.p;
  return s;
}

double max_abs_dev(const Vec4& a, const Vec4& b) {
  return std::max(std::max(std::fabs(a.e - b.e), std::fabs(a.px - b.px)),
                  std::max(std::fabs(a.py - b.py), std::fabs(a.pz - b.pz)));
}

}  // namespace

// ---------------------------------------------------------------------------

TEST_CASE("pythia: the head-on (x, Q2, phi) -> e' construction is the "
          "convention it claims") {
  const BeamConfig bc = mid_config();
  const Vec4 p_n = nucleon_at_rest_in_ion(bc);
  const Vec4 k{kEe, 0.0, 0.0, -kEe};
  for (const Point& pt : points()) {
    Vec4 kp;
    double y = 0.0;
    REQUIRE(dis_scattered_electron(kEe, p_n, pt.x, pt.q2, pt.phi, &kp, &y));
    const Vec4 q = k - kp;
    // The three invariants come back exactly.
    CHECK_CLOSE(-q.m2(), pt.q2, 1e-12);
    CHECK_CLOSE(dot4(p_n, q) / dot4(p_n, k), y, 1e-12);
    CHECK_CLOSE(-q.m2() / (2.0 * dot4(p_n, q)), pt.x, 1e-12);
    // e' is massless and its azimuth is phi.
    CHECK_CLOSE_AT(kp.m2(), 0.0, 0.0, 1e-9);
    double az = std::atan2(kp.py, kp.px);
    if (az < 0.0) az += 2.0 * kPi;
    CHECK_CLOSE(az, pt.phi, 1e-12);
    // y = Q^2 / (x (s - m_N^2)) with s the per-nucleon (k + p_N)^2.
    const double s = (k + p_n).m2();
    CHECK_CLOSE(y, pt.q2 / (pt.x * (s - p_n.m2())), 1e-12);
  }
}

TEST_CASE("pythia: xi puts the outgoing quark on its mass shell, and reduces "
          "to x_Bj for a massless target") {
  const BeamConfig bc = mid_config();
  const Vec4 k{kEe, 0.0, 0.0, -kEe};

  SUBCASE("massive on-shell target") {
    const Vec4 p_n = nucleon_at_rest_in_ion(bc);
    for (const Point& pt : points()) {
      Vec4 kp;
      REQUIRE(dis_scattered_electron(kEe, p_n, pt.x, pt.q2, pt.phi, &kp));
      const Vec4 q = k - kp;
      double xi = 0.0;
      REQUIRE(dis_parton_fraction(q, p_n, &xi));
      const Vec4 out = q + scale4(xi, p_n);
      CHECK_CLOSE_AT(out.m2(), 0.0, 0.0, 1e-9 * pt.q2);
      // and with a charm mass
      double xic = 0.0;
      REQUIRE(dis_parton_fraction(q, p_n, &xic, 1.0, 1.5));
      const Vec4 outc = q + scale4(xic, p_n);
      CHECK_CLOSE(outc.m2(), 1.5 * 1.5, 1e-9);
      CHECK(xic > xi);
    }
  }

  SUBCASE("massless target: xi = x_Bj exactly") {
    const Vec4 p_n{kPu, 0.0, 0.0, kPu};
    for (const Point& pt : points()) {
      Vec4 kp;
      if (!dis_scattered_electron(kEe, p_n, pt.x, pt.q2, pt.phi, &kp)) continue;
      const Vec4 q = k - kp;
      double xi = 0.0;
      REQUIRE(dis_parton_fraction(q, p_n, &xi));
      CHECK_CLOSE(xi, -q.m2() / (2.0 * dot4(q, p_n)), 1e-12);
      CHECK_CLOSE(xi, pt.x, 1e-12);
    }
  }

  SUBCASE("off-shell, Fermi-moving target: still exactly massless out") {
    Vec4 p_n{0.0, 0.2, 0.0, kPu * 1.12};
    p_n.e = std::sqrt(p_n.px * p_n.px + p_n.py * p_n.py + p_n.pz * p_n.pz +
                      M_NUCLEON * M_NUCLEON - 0.05);
    // E ~ 111 GeV against m^2 ~ 0.83 GeV^2: the cancellation costs ~3e-12
    // of relative precision, which is the floor of any m^2 at collider energies.
    CHECK_CLOSE(p_n.m2(), M_NUCLEON * M_NUCLEON - 0.05, 1e-10);
    for (const Point& pt : points()) {
      Vec4 kp;
      if (!dis_scattered_electron(kEe, p_n, pt.x, pt.q2, pt.phi, &kp)) continue;
      const Vec4 q = k - kp;
      double xi = 0.0;
      REQUIRE(dis_parton_fraction(q, p_n, &xi));
      const Vec4 out = q + scale4(xi, p_n);
      CHECK_CLOSE_AT(out.m2(), 0.0, 0.0, 1e-9 * pt.q2);
    }
  }
}

// ---------------------------------------------------------------------------

TEST_CASE("pythia: 200 events on p and on n conserve four-momentum, charge, "
          "and satisfy both HFS truth identities") {
  const BeamConfig bc = mid_config();
  const Vec4 p_n = nucleon_at_rest_in_ion(bc);
  PythiaBridgeOptions opt;
  opt.seed = kSeed;
  PythiaBridge bridge(bc, opt);

  for (int id_n : {2212, 2112}) {
    CAPTURE(id_n);
    double worst_p = 0.0, worst_pt = 0.0, worst_sig = 0.0, worst_lambda = 0.0;
    int n_done = 0;
    long n_charged = 0;
    for (int i = 0; i < 200; ++i) {
      const Point& pt = points()[static_cast<std::size_t>(i) % points().size()];
      Point pti = pt;
      pti.phi = std::fmod(pt.phi + 0.017 * i, 2.0 * kPi);
      Event ev;
      REQUIRE(make_event(bc, pti, p_n, id_n, static_cast<std::uint64_t>(i),
                         &ev));
      Rng rng(kSeed, 0, static_cast<std::uint64_t>(id_n),
              static_cast<std::uint64_t>(i));
      REQUIRE(bridge.hadronize(ev, rng));
      ++n_done;

      // (a) four-momentum: hadrons + e' = k + p_N.
      const Vec4 want = ev.find(Role::BeamElectron)->p + p_n;
      worst_p = std::max(worst_p, max_abs_dev(hadrons_plus_electron(ev), want) /
                                      want.e);

      // (b) charge: the hadronic system carries the nucleon's charge.
      const HfsSummary h = hfs_summary(ev);
      CHECK_CLOSE_AT(h.charge, (id_n == 2212) ? 1.0 : 0.0, 0.0, 1e-9);

      // (c) Sum (E - p_z)_hadrons = 2 E_e y + m_N^2/(E_N + p_N).
      const double truth = hfs_sigma_empz_truth(kEe, ev.kin.y, p_n);
      worst_sig = std::max(worst_sig,
                           std::fabs(h.sigma_empz - truth) / std::fabs(truth));
      // and the exact form, to the numerical floor.
      CHECK_CLOSE(h.sigma_empz,
                  hfs_sigma_empz_exact(ev.find(Role::BeamElectron)->p, p_n,
                                       ev.find(Role::ScatteredElectron)->p),
                  1e-10);

      // (d) |Sum p_T,hadrons| = p_T,e' for a collinear target.
      worst_pt = std::max(worst_pt,
                          std::fabs(h.pt - ev.find(Role::ScatteredElectron)->p.pt()));

      // (e) exactly one scattered electron in the record; PYTHIA's copy of it
      //     is never carried over as a hadron.
      int n_scat = 0, n_hadron_e = 0;
      for (const auto& p : ev.particles) {
        if (p.role == Role::ScatteredElectron) ++n_scat;
        if (p.role == Role::Hadron && p.pdg == 11 &&
            std::fabs(p.p.e - ev.find(Role::ScatteredElectron)->p.e) < 1e-6)
          ++n_hadron_e;
      }
      CHECK(n_scat == 1);
      CHECK(n_hadron_e == 0);

      // (f) the T0 pseudo-particle was demoted.
      REQUIRE(ev.find(Role::HadronicX) != nullptr);
      CHECK(ev.find(Role::HadronicX)->status == Status::Intermediate);

      worst_lambda = std::max(worst_lambda, std::fabs(bridge.last_rescale() - 1.0));
      n_charged += h.n_charged;
    }
    CHECK(n_done == 200);
    MESSAGE("id_n = " << id_n << ": worst relative 4-momentum deviation "
                      << worst_p << ", worst |dp_T| " << worst_pt
                      << " GeV, worst relative Sigma(E-p_z) deviation "
                      << worst_sig << ", worst |lambda - 1| " << worst_lambda
                      << ", <n_charged> = "
                      << double(n_charged) / double(n_done));
    CHECK(worst_p < 1e-6);
    CHECK(worst_pt < 1e-8);
    CHECK(worst_sig < 1e-3);
    // The frame map is a *pure* Lorentz transformation whenever the surrogate
    // is realised exactly; lambda != 1 would mean PYTHIA moved our e'.
    CHECK(worst_lambda < 1e-6);
  }
  CHECK(bridge.stats().n_failed == 0);
  CHECK(bridge.stats().n_proton == 200);
  CHECK(bridge.stats().n_neutron == 200);
}

TEST_CASE("pythia: a Fermi-moving, off-shell target still conserves "
          "four-momentum and charge") {
  const BeamConfig bc = mid_config();
  PythiaBridgeOptions opt;
  opt.seed = kSeed + 1;
  PythiaBridge bridge(bc, opt);

  // 0.2 GeV of transverse momentum and an off-shellness of -0.05 GeV^2,
  // with the longitudinal momentum swung around the nominal by +-12 %.
  double worst_p = 0.0, worst_pt = 0.0;
  int n_done = 0;
  for (int i = 0; i < 60; ++i) {
    const double phi_n = 0.31 * i;
    const double frac = 1.0 + 0.12 * std::sin(0.7 * i);
    Vec4 p_n{0.0, 0.2 * std::cos(phi_n), 0.2 * std::sin(phi_n), kPu * frac};
    const double m2 = M_NUCLEON * M_NUCLEON - 0.05;
    p_n.e = std::sqrt(p_n.px * p_n.px + p_n.py * p_n.py + p_n.pz * p_n.pz + m2);
    REQUIRE(std::fabs(p_n.m2() - m2) < 1e-9);

    const Point& pt = points()[static_cast<std::size_t>(i) % points().size()];
    Event ev;
    // dis_scattered_electron is derived for a collinear p_N; the transverse
    // 0.2 GeV is a 0.2 % effect on the invariants and is irrelevant here --
    // the bridge takes (k, k', p_N) as given and never re-derives them.
    if (!make_event(bc, pt, p_n, (i % 2) ? 2212 : 2112,
                    static_cast<std::uint64_t>(i), &ev))
      continue;
    Rng rng(kSeed + 1, 0, 0, static_cast<std::uint64_t>(i));
    if (!bridge.hadronize(ev, rng)) continue;
    ++n_done;

    const Vec4 want = ev.find(Role::BeamElectron)->p + p_n;
    worst_p = std::max(worst_p, max_abs_dev(hadrons_plus_electron(ev), want) /
                                    want.e);
    const HfsSummary h = hfs_summary(ev);
    CHECK_CLOSE_AT(h.charge, (i % 2) ? 1.0 : 0.0, 0.0, 1e-9);
    // For a moving target the p_T identity picks up the target's own p_T:
    //   Sum p_T,hadrons = p_N,T - p_T,e'.
    const Vec4 kp = ev.find(Role::ScatteredElectron)->p;
    worst_pt = std::max(
        worst_pt, std::max(std::fabs(h.px - (p_n.px - kp.px)),
                           std::fabs(h.py - (p_n.py - kp.py))));
    // and Sigma(E - p_z) is exact in the general form.
    CHECK_CLOSE(h.sigma_empz,
                hfs_sigma_empz_exact(ev.find(Role::BeamElectron)->p, p_n, kp),
                1e-10);
  }
  MESSAGE("Fermi-moving: " << n_done << " events, worst relative 4-momentum "
                           << worst_p << ", worst |dp_T| " << worst_pt);
  CHECK(n_done > 50);
  CHECK(worst_p < 1e-6);
  CHECK(worst_pt < 1e-8);
}

TEST_CASE("pythia: the same (seed, run, bunch, event) gives the identical "
          "event") {
  const BeamConfig bc = mid_config();
  const Vec4 p_n = nucleon_at_rest_in_ion(bc);

  auto run_one = [&](std::uint64_t seed, int i, PythiaBridge& bridge) {
    Event ev;
    REQUIRE(make_event(bc, points()[static_cast<std::size_t>(i) % points().size()],
                       p_n, (i % 2) ? 2212 : 2112, static_cast<std::uint64_t>(i),
                       &ev));
    Rng rng(seed, 3, 7, static_cast<std::uint64_t>(i));
    REQUIRE(bridge.hadronize(ev, rng));
    return ev;
  };

  PythiaBridgeOptions opt;
  opt.seed = kSeed;
  PythiaBridge a(bc, opt);
  PythiaBridge b(bc, opt);

  // Deliberately different orders: b runs the events backwards, so only a
  // per-event counter-based stream can make the two agree.
  std::vector<Event> ea, eb(8);
  for (int i = 0; i < 8; ++i) ea.push_back(run_one(kSeed, i, a));
  for (int i = 7; i >= 0; --i) eb[static_cast<std::size_t>(i)] = run_one(kSeed, i, b);

  for (std::size_t i = 0; i < ea.size(); ++i) {
    REQUIRE(ea[i].particles.size() == eb[i].particles.size());
    for (std::size_t j = 0; j < ea[i].particles.size(); ++j) {
      CHECK(ea[i].particles[j].pdg == eb[i].particles[j].pdg);
      CHECK(ea[i].particles[j].p.e == doctest::Approx(eb[i].particles[j].p.e));
      CHECK(ea[i].particles[j].p.px == doctest::Approx(eb[i].particles[j].p.px));
      CHECK(ea[i].particles[j].p.py == doctest::Approx(eb[i].particles[j].p.py));
      CHECK(ea[i].particles[j].p.pz == doctest::Approx(eb[i].particles[j].p.pz));
    }
  }
}

TEST_CASE("pythia: an implicit target is drawn Z : N and a cluster target "
          "uses the hook") {
  BeamConfig bc = mid_config();
  PythiaBridgeOptions opt;
  opt.seed = kSeed + 2;
  PythiaBridge bridge(bc, opt);
  const Vec4 p_n = nucleon_at_rest_in_ion(bc);

  SUBCASE("inclusive fallback: no explicit target") {
    int n_p = 0, n_n = 0;
    for (int i = 0; i < 60; ++i) {
      Event ev;
      REQUIRE(make_event(bc, points()[static_cast<std::size_t>(i) % points().size()],
                         p_n, 2212, static_cast<std::uint64_t>(i), &ev));
      // drop the explicit struck nucleon
      std::vector<Particle> keep;
      for (const auto& p : ev.particles)
        if (p.role != Role::StruckNucleon) keep.push_back(p);
      ev.particles = keep;
      Rng rng(kSeed + 2, 0, 0, static_cast<std::uint64_t>(i));
      REQUIRE(bridge.hadronize(ev, rng));
      // the bridge documents its choice by appending the nucleon it used
      const Particle* sn = ev.find(Role::StruckNucleon);
      REQUIRE(sn != nullptr);
      CHECK(sn->status == Status::Intermediate);
      if (sn->pdg == 2212) ++n_p; else ++n_n;
      // p_N = P_ion / A, the nucleon at rest in the ion rest frame
      CHECK_CLOSE(sn->p.pz, kPu, 1e-12);
      // P3: the free nucleon mass, the same one InclusiveGenerator uses
      CHECK_CLOSE(sn->p.m2(), M_NUCLEON * M_NUCLEON, 1e-9);
      CHECK(std::fabs(sn->p.m2()
                      - bc.ion.mass_per_nucleon() * bc.ion.mass_per_nucleon())
            > 1e-3);
    }
    // 6Li is Z = N = 3, so the split must be compatible with 50:50.
    MESSAGE("inclusive fallback drew " << n_p << " p and " << n_n << " n");
    CHECK(n_p > 15);
    CHECK(n_n > 15);
  }

  SUBCASE("overridable choice") {
    PythiaBridgeOptions o2 = opt;
    o2.nucleon_choice = NucleonChoice::Neutron;
    PythiaBridge nb(bc, o2);
    Event ev;
    REQUIRE(make_event(bc, points()[2], p_n, 2212, 0, &ev));
    std::vector<Particle> keep;
    for (const auto& p : ev.particles)
      if (p.role != Role::StruckNucleon) keep.push_back(p);
    ev.particles = keep;
    Rng rng(1, 0, 0, 0);
    REQUIRE(nb.hadronize(ev, rng));
    CHECK(nb.last_struck_nucleon_pdg() == 2112);
  }

  SUBCASE("cluster target: v0 is p_cluster / A_c, hook overrides it") {
    // A deuteron cluster inside 6Li, carrying 2/6 of the ion momentum.
    Event ev;
    REQUIRE(make_event(bc, points()[2], p_n, 2212, 0, &ev));
    std::vector<Particle> keep;
    for (const auto& p : ev.particles)
      if (p.role != Role::StruckNucleon) keep.push_back(p);
    Particle cl;
    cl.pdg = 1000010020;   // 2H
    cl.status = Status::Intermediate;
    cl.role = Role::StruckCluster;
    cl.p = {2.0 * p_n.e, 0.0, 0.0, 2.0 * kPu};
    keep.push_back(cl);
    ev.particles = keep;
    Rng rng(5, 0, 0, 0);
    REQUIRE(bridge.hadronize(ev, rng));
    CHECK_CLOSE(bridge.last_struck_nucleon().pz, kPu, 1e-12);
    CHECK_CLOSE(bridge.last_struck_nucleon().m2(), M_NUCLEON * M_NUCLEON, 1e-9);

    bool hook_ran = false;
    bridge.set_nucleon_in_cluster(
        [&](const Vec4& pc, int a_c, int z_c, Rng&, int* pdg) {
          hook_ran = true;
          CHECK(a_c == 2);
          CHECK(z_c == 1);
          *pdg = 2112;
          Vec4 v{0.0, 0.0, 0.0, pc.pz / a_c};
          v.e = std::sqrt(v.pz * v.pz + M_NUCLEON * M_NUCLEON);
          return v;
        });
    Event ev2 = ev;
    ev2.particles.erase(
        std::remove_if(ev2.particles.begin(), ev2.particles.end(),
                       [](const Particle& p) {
                         return p.role == Role::Hadron ||
                                p.role == Role::StruckNucleon;
                       }),
        ev2.particles.end());
    for (auto& p : ev2.particles)
      if (p.role == Role::HadronicX) p.status = Status::Final;
    Rng rng2(6, 0, 0, 0);
    REQUIRE(bridge.hadronize(ev2, rng2));
    CHECK(hook_ran);
    CHECK(bridge.last_struck_nucleon_pdg() == 2112);
  }
}

// ---------------------------------------------------------------------------

TEST_CASE("pythia: the charged multiplicity agrees with a stock "
          "WeakBosonExchange run in the same (x, Q2) window") {
  // Stock reference: exactly the configuration of
  // PolarizedLithiumSim/tools/pythia8/gen_dis_hfs.py -- head-on beams,
  // ff2ff(t:gmZ), the two silent-cut fixes -- restricted after the fact to
  // the window the bridge is exercised in.  The bridge sees no PhaseSpace
  // cut at all (its hard process comes from LHAup), so this is a comparison
  // of the SHOWER + HADRONIZATION only.
  const double kQ2Lo = 4.0, kQ2Hi = 30.0;
  const double kXLo = 0.005, kXHi = 0.10;

  Pythia8::Pythia stock(std::string(LIPOLGEN_PYTHIA8_XMLDOC), false);
  for (const char* s : {"Beams:frameType = 2",
                        "Beams:idA = 2212",
                        "Beams:idB = 11",
                        "WeakBosonExchange:ff2ff(t:gmZ) = on",
                        "PhaseSpace:Q2Min = 4.0",
                        "PhaseSpace:pTHatMinDiverge = 0.5",
                        "PhaseSpace:mHatMin = 0.5",
                        "SpaceShower:dipoleRecoil = on",
                        "SpaceShower:pTmaxMatch = 2",
                        "PDF:lepton = off",
                        "TimeShower:QEDshowerByL = off",
                        "Random:setSeed = on",
                        "Random:seed = 101",
                        "Next:numberCount = 0",
                        "Next:numberShowInfo = 0",
                        "Next:numberShowProcess = 0",
                        "Next:numberShowEvent = 0",
                        "Print:quiet = on"})
    REQUIRE(stock.readString(s));
  REQUIRE(stock.readString("Beams:eA = " + std::to_string(kPu)));
  REQUIRE(stock.readString("Beams:eB = " + std::to_string(kEe)));
  REQUIRE(stock.init());

  const Vec4 k{kEe, 0.0, 0.0, -kEe};
  long stock_charged = 0;
  long stock_n = 0;
  double stock_x = 0.0, stock_q2 = 0.0;
  for (int i = 0; i < 2000 && stock_n < 2000; ++i) {
    if (!stock.next()) continue;
    const Pythia8::Event& pev = stock.event;
    int i_e = -1;
    double best = -1.0;
    for (int j = 0; j < pev.size(); ++j)
      if (pev[j].isFinal() && pev[j].id() == 11 && pev[j].e() > best) {
        best = pev[j].e();
        i_e = j;
      }
    if (i_e < 0) continue;
    const Vec4 p_n = {pev[1].e(), pev[1].px(), pev[1].py(), pev[1].pz()};
    const Vec4 kp = {pev[i_e].e(), pev[i_e].px(), pev[i_e].py(), pev[i_e].pz()};
    const Vec4 q = k - kp;
    const double q2 = -q.m2();
    const double x = q2 / (2.0 * dot4(p_n, q));
    if (q2 < kQ2Lo || q2 > kQ2Hi || x < kXLo || x > kXHi) continue;
    int nch = 0;
    for (int j = 0; j < pev.size(); ++j)
      if (pev[j].isFinal() && j != i_e && pev[j].isCharged()) ++nch;
    stock_charged += nch;
    stock_x += x;
    stock_q2 += q2;
    ++stock_n;
  }
  REQUIRE(stock_n > 200);
  const double stock_mean = double(stock_charged) / double(stock_n);
  const double xbar = stock_x / double(stock_n);
  const double q2bar = stock_q2 / double(stock_n);

  // Bridge: the same number of events, thrown at the stock run's own mean
  // (x, Q^2) so the two are compared at the same point of the window.
  BeamConfig bc = mid_config();
  PythiaBridgeOptions opt;
  opt.seed = kSeed + 3;
  PythiaBridge bridge(bc, opt);
  const Vec4 p_n_rest{std::sqrt(kPu * kPu + PROTON_MASS * PROTON_MASS), 0.0, 0.0,
                      kPu};
  long bridge_charged = 0, bridge_n = 0;
  for (int i = 0; i < 2000; ++i) {
    Point pt{xbar, q2bar, std::fmod(0.013 * i, 2.0 * kPi)};
    Event ev;
    if (!make_event(bc, pt, p_n_rest, 2212, static_cast<std::uint64_t>(i), &ev))
      continue;
    Rng rng(kSeed + 3, 0, 0, static_cast<std::uint64_t>(i));
    if (!bridge.hadronize(ev, rng)) continue;
    bridge_charged += hfs_summary(ev).n_charged;
    ++bridge_n;
  }
  REQUIRE(bridge_n > 1500);
  const double bridge_mean = double(bridge_charged) / double(bridge_n);

  std::printf(
      "\n[pythia] charged multiplicity, e + p at 10 x %.1f GeV/u,"
      " Q2 in [%.1f, %.1f], x in [%.3f, %.2f]\n"
      "         stock WeakBosonExchange : %.3f   (%ld events,"
      " <x> = %.4f, <Q2> = %.2f)\n"
      "         LiPolGen PythiaBridge   : %.3f   (%ld events at that"
      " exact (x, Q2))\n"
      "         ratio bridge/stock      : %.3f\n",
      kPu, kQ2Lo, kQ2Hi, kXLo, kXHi, stock_mean, stock_n, xbar, q2bar,
      bridge_mean, bridge_n, bridge_mean / stock_mean);

  // Loose on purpose: the stock run averages over the whole window while the
  // bridge sits at its mean point, and <n_ch> grows like ln W^2.
  CHECK(bridge_mean / stock_mean > 0.80);
  CHECK(bridge_mean / stock_mean < 1.20);
}

// P2.  The flavour pool used to be weighted at the MASSLESS zeta for every
// flavour, while c and b were then produced at their own, larger
// zeta_q = (Q^2 + m_q^2)/(P_A^+ qtil^-).  Two things went wrong: the pool was
// priced at one momentum fraction and drawn at another, and a heavy flavour
// whose zeta_q had already run past 1 -- which is exactly what happens at low
// Q^2, where zeta_light is small but m_c^2/(P_A^+ qtil^-) is not -- stayed in
// the pool and burned a pythia.next() retry every time it was picked.
TEST_CASE("pythia: heavy flavours out of phase space cost no retries") {
  BeamConfig bc = mid_config();
  PythiaBridgeOptions opt;
  opt.seed = kSeed + 31;
  opt.include_strange = true;
  opt.include_charm = true;
  opt.include_bottom = true;   // b is the flavour that runs out of zeta first
  PythiaBridge bridge(bc, opt);
  const Vec4 p_n = nucleon_at_rest_in_ion(bc);

  // A scan of the generator window: zeta_q = (Q^2 + m_q^2)/(P_A^+ qtil^-)
  // exceeds the light zeta by (1 + m_q^2/Q^2), so a heavy flavour runs out of
  // phase space at LARGE x and the SMALLEST Q^2 the W^2 >= 8 cut allows.
  std::vector<Point> scan;
  for (double x : {0.02, 0.1, 0.3, 0.5, 0.7, 0.8, 0.85}) {
    const double q2_min = std::max(0.8, 7.2 * x / (1.0 - x));
    for (double f : {1.0, 1.6, 4.0}) {
      scan.push_back(Point{x, q2_min * f, 0.3 + 1.7 * f});
    }
  }
  int n_ok = 0, n_built = 0;
  for (std::size_t i = 0; i < scan.size(); ++i) {
    for (int j = 0; j < 4; ++j) {
      Event ev;
      const std::uint64_t num = static_cast<std::uint64_t>(4 * i + j);
      if (!make_event(bc, scan[i], p_n, 2212, num, &ev)) continue;
      ++n_built;
      Rng rng(kSeed + 31, 0, 0, num);
      if (bridge.hadronize(ev, rng)) ++n_ok;
    }
  }
  const PythiaBridgeStats& st = bridge.stats();
  MESSAGE("flavour-pool scan: " << n_ok << "/" << st.n_called
          << " hadronized, " << st.n_flavour_dropped
          << " flavour offers dropped for zeta_q >= 1, " << st.n_retries
          << " retries, " << st.n_failed << " failures");
  CHECK(n_built > 50);
  CHECK(n_ok > n_built - 5);
  // The window really does drive heavy flavours out of phase space ...
  CHECK(st.n_flavour_dropped > 0);
  // ... and THE POINT: not one retry is spent on them.  Under the old
  // weighting they stayed in the pool, priced at the massless zeta, and every
  // pick of one burned a pythia.next().
  CHECK(st.n_retries == 0);
  CHECK(st.n_failed == 0);

  // Light flavours alone never leave the pool over the same scan.
  PythiaBridgeOptions light = opt;
  light.include_charm = false;
  light.include_bottom = false;
  PythiaBridge lb(bc, light);
  for (std::size_t i = 0; i < scan.size(); ++i) {
    for (int j = 0; j < 4; ++j) {
      Event ev;
      const std::uint64_t num = static_cast<std::uint64_t>(4 * i + j);
      if (!make_event(bc, scan[i], p_n, 2212, num, &ev)) continue;
      Rng rng(kSeed + 31, 0, 0, num);
      (void)lb.hadronize(ev, rng);
    }
  }
  CHECK(lb.stats().n_flavour_dropped == 0);
  CHECK(lb.stats().n_retries == 0);

  // THE MISPRICING, which is the other half of P2.  A heavy flavour is
  // produced at zeta_q = zeta_light (1 + m_q^2/Q^2) -- 1.23 at Q^2 = 10, 1.70
  // at Q^2 = 3.2 -- and its pool weight x f_q used to be read at zeta_light
  // instead, i.e. at a momentum fraction the event never uses, on a
  // distribution that falls steeply in between.  The zeta that reaches PYTHIA
  // is the one the flavour is now weighted at, and it carries the mass:
  PythiaBridgeOptions copt;
  copt.seed = kSeed + 51;
  copt.include_charm = true;
  copt.include_bottom = false;
  PythiaBridge cb(bc, copt);
  int n_c = 0, n_tot = 0;
  double worst_ratio_dev = 0.0;
  double zeta_light_ref = 0.0, zeta_c_ref = 0.0;
  const Point cpt{0.02, 12.0, 2.1};
  for (int i = 0; i < 400; ++i) {
    Event ev;
    if (!make_event(bc, cpt, p_n, 2212, static_cast<std::uint64_t>(i), &ev)) continue;
    Rng rng(kSeed + 51, 0, 0, static_cast<std::uint64_t>(i));
    if (!cb.hadronize(ev, rng)) continue;
    ++n_tot;
    const double z = cb.last_xi_pythia();
    if (std::abs(cb.last_quark_id()) == 4) {
      ++n_c;
      zeta_c_ref = z;
    } else {
      zeta_light_ref = z;   // every light flavour shares one zeta
    }
  }
  REQUIRE(n_tot > 350);
  REQUIRE(n_c > 0);
  REQUIRE(zeta_light_ref > 0.0);
  const double m_c = 1.5;   // the PYTHIA table charm mass, to ~1 %
  const double want_ratio = 1.0 + m_c * m_c / cpt.q2;
  const double got_ratio = zeta_c_ref / zeta_light_ref;
  worst_ratio_dev = std::fabs(got_ratio / want_ratio - 1.0);
  MESSAGE("charm picked in " << n_c << "/" << n_tot << " events at x = "
          << cpt.x << ", Q2 = " << cpt.q2 << "; zeta_c/zeta_light = "
          << got_ratio << " against 1 + m_c^2/Q^2 = " << want_ratio);
  // the produced zeta really does carry the quark mass ...
  CHECK(got_ratio > 1.0);
  CHECK(worst_ratio_dev < 0.05);
  // ... and it is a different momentum fraction from the light one, which is
  // the one the old code priced charm at
  CHECK(zeta_c_ref > zeta_light_ref * 1.1);
}

// P1.  The bridge's own fallback -- used only when the caller drives it
// without a `Role::StruckNucleon`, which a Pipeline event never does -- is
// weighted by the same structure functions as the generator's draw.
TEST_CASE("pythia: the implicit-target species follows Z F2p : N F2n") {
  BeamConfig bc = mid_config();
  const ToyF2 f2;
  const Ion& ion = bc.ion;
  const double x = 0.4, q2 = 40.0;
  const double want = ion.Z * f2.f2p(x, q2)
                      / (ion.Z * f2.f2p(x, q2) + ion.N() * f2.f2n(x, q2));
  CHECK(want > 0.55);             // 6Li is N = Z, and yet not 0.5
  PythiaBridgeOptions opt;
  opt.seed = kSeed + 41;
  CHECK(opt.nucleon_choice == NucleonChoice::ByStructureFunctions);
  PythiaBridge bridge(bc, opt);
  const Vec4 p_n = nucleon_at_rest_in_ion(bc);

  const Point pt{x, q2, 1.3};
  const int n = 400;
  int n_p = 0, n_tried = 0;
  for (int i = 0; i < n; ++i) {
    Event ev;
    if (!make_event(bc, pt, p_n, 2212, static_cast<std::uint64_t>(i), &ev)) continue;
    std::vector<Particle> keep;
    for (const auto& p : ev.particles)
      if (p.role != Role::StruckNucleon) keep.push_back(p);
    ev.particles = keep;
    Rng rng(kSeed + 41, 0, 0, static_cast<std::uint64_t>(i));
    if (!bridge.hadronize(ev, rng)) continue;
    ++n_tried;
    if (bridge.last_struck_nucleon_pdg() == 2212) ++n_p;
  }
  REQUIRE(n_tried > 300);
  const double got = static_cast<double>(n_p) / n_tried;
  const double err = std::sqrt(want * (1.0 - want) / n_tried);
  MESSAGE("bridge implicit target at x = " << x << ": p in " << got << " +- "
          << err << " of events, Z F2p/(Z F2p + N F2n) = " << want
          << " (flat Z/A would be " << 0.5 << ")");
  CHECK(std::fabs(got - want) < 4.0 * err);

  // the old flat rule is still reachable by name
  PythiaBridgeOptions zn = opt;
  zn.nucleon_choice = NucleonChoice::ByZN;
  PythiaBridge zb(bc, zn);
  int zn_p = 0, zn_tried = 0;
  for (int i = 0; i < n; ++i) {
    Event ev;
    if (!make_event(bc, pt, p_n, 2212, static_cast<std::uint64_t>(i), &ev)) continue;
    std::vector<Particle> keep;
    for (const auto& p : ev.particles)
      if (p.role != Role::StruckNucleon) keep.push_back(p);
    ev.particles = keep;
    Rng rng(kSeed + 41, 0, 0, static_cast<std::uint64_t>(i));
    if (!zb.hadronize(ev, rng)) continue;
    ++zn_tried;
    if (zb.last_struck_nucleon_pdg() == 2212) ++zn_p;
  }
  REQUIRE(zn_tried > 300);
  const double zn_got = static_cast<double>(zn_p) / zn_tried;
  CHECK(std::fabs(zn_got - 0.5) < 4.0 * std::sqrt(0.25 / zn_tried));
  CHECK(zn_got < got);
}
