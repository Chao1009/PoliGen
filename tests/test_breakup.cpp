// SPDX-License-Identifier: GPL-3.0-or-later
// Tier T1: the struck cluster resolved into a struck nucleon + partner
// spectator(s).  See breakup.hpp for every physics choice these tests pin.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <numeric>
#include <vector>

#include "check_close.hpp"
#include "doctest.h"

#include "lipolgen/breakup.hpp"
#include "lipolgen/pipeline.hpp"

using namespace lipolgen;

namespace {

constexpr std::uint64_t kSeed = 20260830;

double quantile(std::vector<double> v, double q) {
  std::sort(v.begin(), v.end());
  const double pos = q * static_cast<double>(v.size() - 1);
  const std::size_t i = static_cast<std::size_t>(pos);
  const double t = pos - static_cast<double>(i);
  return (i + 1 < v.size()) ? v[i] * (1.0 - t) + v[i + 1] * t : v.back();
}

double mean(const std::vector<double>& v) {
  return std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(v.size());
}

/// A struck cluster AT REST of invariant mass `m`: the frame in which the
/// drawn internal momentum is the partner's momentum, one to one.
Vec4 at_rest(double m) { return {m, 0.0, 0.0, 0.0}; }

PipelineConfig tagged_cfg(PipelineChannel ch, std::uint64_t n) {
  PipelineConfig cfg;
  cfg.channel = ch;
  cfg.isotope = channel_isotope(ch);
  cfg.beam_config = 1;
  cfg.n_events = n;
  cfg.seed = kSeed;
  // The lithium tagging working point is tabulated for 6Li and 7Li only; the
  // deuteron control keeps the Yellow Report envelope.
  if (ch != PipelineChannel::TaggedDeuteronP)
    cfg.optics_choice = OpticsChoice::Tagging;
  return cfg;
}

RunPlan plan_for_channel(PipelineChannel ch) {
  return (ch == PipelineChannel::TaggedLi7Alpha)
             ? helicity_flip_plan(1.5, 0.7, 0.7)
             : tensor_thirds_plan(0.7, 0.6);
}

}  // namespace

// ---------------------------------------------------------------------------

TEST_CASE("breakup: species table and the AME2020 separation energies") {
  CHECK(cluster_species(1, 1) == ClusterSpecies::Nucleon);
  CHECK(cluster_species(0, 1) == ClusterSpecies::Nucleon);
  CHECK(cluster_species(1, 2) == ClusterSpecies::Deuteron);
  CHECK(cluster_species(1, 3) == ClusterSpecies::Triton);
  CHECK_THROWS_AS(cluster_species(2, 4), std::runtime_error);

  // The three numbers breakup.hpp quotes, straight off the mass table that
  // spectator.cpp already carries -- so nothing is hard-coded twice.
  const double m_p = nuclear_mass(1, 1), m_n = nuclear_mass(0, 1);
  const double m_d = nuclear_mass(1, 2), m_t = nuclear_mass(1, 3);
  const double s_d = m_p + m_n - m_d;
  const double s_n_t = m_n + m_d - m_t;
  const double s_3body_t = m_p + 2.0 * m_n - m_t;
  MESSAGE("S_d = " << 1e3 * s_d << " MeV, S_n(3H) = " << 1e3 * s_n_t
                   << " MeV, S(p+n+n) = " << 1e3 * s_3body_t << " MeV");
  CHECK_CLOSE(1e3 * s_d, 2.2246, 1e-3);
  CHECK_CLOSE(1e3 * s_n_t, 6.2572, 1e-3);
  CHECK_CLOSE(1e3 * s_3body_t, 8.4818, 1e-3);
}

TEST_CASE("breakup: the deuteron partner reproduces the S-wave marginal at "
          "P_D = 0 (the spectator.py Hulthen sampler)") {
  BreakupOptions bo;
  bo.p_d = 0.0;                     // pure S wave -> isotropic, exact marginal
  ClusterBreakup br(bo);
  const double m_d = nuclear_mass(1, 2);

  std::vector<double> k;
  k.reserve(40000);
  Rng rng(kSeed, 1, 0, 0);
  BreakupInput in;
  in.species = ClusterSpecies::Deuteron;
  in.p_cluster = at_rest(m_d);
  in.z = 1;
  in.a = 2;
  in.m_s = 1.0;
  in.s_cluster = 1.0;
  in.x = 0.1;
  in.q2 = 10.0;
  BreakupResult out;
  for (int i = 0; i < 40000; ++i) {
    REQUIRE(br.resolve(in, rng, out));
    REQUIRE(out.partners.size() == 1);
    // At rest the partner's three-momentum IS the drawn relative momentum.
    CHECK_CLOSE(out.partners[0].p.p(), out.k, 1e-9);
    k.push_back(out.k);
  }

  // The independent reference: `spectator.MomentumSampler` on the deuteron
  // channel, a 30000-point inverse CDF of k^2 |psi_Hulthen|^2 -- a different
  // grid, a different code path, the same physics.
  MomentumSampler ref(DEUTERON_P_TAG(), BETA_DEFAULT);
  std::vector<double> want;
  want.reserve(40000);
  for (int i = 0; i < 40000; ++i)
    want.push_back(ref.k_of_u((i + 0.5) / 40000.0));

  const double q25 = quantile(k, 0.25), q50 = quantile(k, 0.5),
               q75 = quantile(k, 0.75), q95 = quantile(k, 0.95);
  const double r25 = quantile(want, 0.25), r50 = quantile(want, 0.5),
               r75 = quantile(want, 0.75), r95 = quantile(want, 0.95);
  MESSAGE("deuteron internal |k| [GeV] quantiles 25/50/75/95: "
          << q25 << "/" << q50 << "/" << q75 << "/" << q95
          << " against the Hulthen sampler's " << r25 << "/" << r50 << "/"
          << r75 << "/" << r95 << ", mean " << mean(k) << " vs " << mean(want));
  CHECK_CLOSE(q25, r25, 0.05);
  CHECK_CLOSE(q50, r50, 0.05);
  CHECK_CLOSE(q75, r75, 0.05);
  CHECK_CLOSE(q95, r95, 0.08);
}

TEST_CASE("breakup: the deuteron's struck nucleon carries the D-state-diluted "
          "polarization (1 - 3/2 P_D) m_S") {
  BreakupOptions bo;                       // P_D = P_D_DEUTERON = 0.045
  ClusterBreakup br(bo);
  const double m_d = nuclear_mass(1, 2);
  const double want = 1.0 - 1.5 * P_D_DEUTERON;

  for (double m_s : {1.0, 0.0, -1.0}) {
    Rng rng(kSeed, 2, 0, static_cast<std::uint64_t>(m_s + 2.0));
    BreakupInput in;
    in.species = ClusterSpecies::Deuteron;
    in.p_cluster = at_rest(m_d);
    in.z = 1; in.a = 2;
    in.m_s = m_s;
    in.s_cluster = 1.0;
    in.x = 0.1; in.q2 = 10.0;
    BreakupResult out;
    double sum = 0.0, sum_partner = 0.0;
    const int n = 200000;
    for (int i = 0; i < n; ++i) {
      REQUIRE(br.resolve(in, rng, out));
      CHECK((out.struck.pol == 1.0 || out.struck.pol == -1.0));
      sum += out.struck.pol;
      sum_partner += out.partners[0].pol;
    }
    const double got = sum / n;
    const double got_p = sum_partner / n;
    MESSAGE("m_S = " << m_s << ": <P_N,struck> = " << got << ", <P_N,partner> = "
                     << got_p << ", (1 - 3/2 P_D) m_S = " << want * m_s);
    CHECK(std::fabs(got - want * m_s) < 0.01);
    // The pair is CORRELATED, not independent: m_1 + m_2 = m_sc always.
    CHECK(std::fabs(got_p - want * m_s) < 0.01);
  }
}

TEST_CASE("breakup: the struck species follows F2p : F2n, not flat 1 : 1") {
  ClusterBreakup br;
  const double m_d = nuclear_mass(1, 2);
  const ToyF2 f2;
  for (double x : {0.05, 0.5}) {
    Rng rng(kSeed, 3, 0, static_cast<std::uint64_t>(1000 * x));
    BreakupInput in;
    in.species = ClusterSpecies::Deuteron;
    in.p_cluster = at_rest(m_d);
    in.z = 1; in.a = 2;
    in.m_s = 1.0; in.s_cluster = 1.0;
    in.x = x; in.q2 = 10.0;
    BreakupResult out;
    int n_p = 0;
    const int n = 100000;
    for (int i = 0; i < n; ++i) {
      REQUIRE(br.resolve(in, rng, out));
      if (out.struck.pdg == 2212) ++n_p;
    }
    const double got = static_cast<double>(n_p) / n;
    const double want = f2.f2p(x, 10.0) / (f2.f2p(x, 10.0) + f2.f2n(x, 10.0));
    MESSAGE("x = " << x << ": struck proton in " << got
                   << " of breakups, F2p/(F2p + F2n) = " << want);
    CHECK(std::fabs(got - want) < 0.006);
  }
}

TEST_CASE("breakup: the triton goes to n + d or p + (nn) with the right "
          "charges, baryon number and polarizations") {
  ClusterBreakup br;
  const double m_t = nuclear_mass(1, 3);
  const Ion& t = TRITON();

  for (double m_s : {0.5, -0.5}) {
    Rng rng(kSeed, 4, 0, m_s > 0 ? 0 : 1);
    BreakupInput in;
    in.species = ClusterSpecies::Triton;
    in.p_cluster = at_rest(m_t);
    in.z = 1; in.a = 3;
    in.m_s = m_s;
    in.s_cluster = 0.5;
    in.x = 0.1; in.q2 = 10.0;
    in.eff_pol_p = t.eff_pol_p;
    in.eff_pol_n = t.eff_pol_n;
    BreakupResult out;
    const int n = 200000;
    int n_p = 0;
    double pol_p = 0.0, pol_n = 0.0;
    int n_pol_p = 0, n_pol_n = 0;
    std::vector<double> k_d;
    for (int i = 0; i < n; ++i) {
      REQUIRE(br.resolve(in, rng, out));
      // charge and baryon number of the whole breakup
      double q = out.struck.charge;
      int a = 1;
      for (const Particle& f : out.partners) {
        q += f.charge;
        a += (std::abs(f.pdg) > 1000000000) ? (std::abs(f.pdg) / 10) % 1000 : 1;
      }
      REQUIRE(q == 1.0);
      REQUIRE(a == 3);
      // the four-momentum closes exactly
      Vec4 s = out.struck.p;
      for (const Particle& f : out.partners) s = s + f.p;
      REQUIRE(std::fabs(s.e - m_t) < 1e-9);
      REQUIRE(std::fabs(s.px) < 1e-9);
      REQUIRE(std::fabs(s.py) < 1e-9);
      REQUIRE(std::fabs(s.pz) < 1e-9);
      if (out.struck.pdg == 2212) {
        ++n_p;
        REQUIRE(out.partners.size() == 2);
        CHECK(out.partners[0].pdg == 2112);
        CHECK(out.partners[1].pdg == 2112);
        pol_p += out.struck.pol;
        ++n_pol_p;
      } else {
        REQUIRE(out.partners.size() == 1);
        CHECK(out.partners[0].pdg == nuclide_pdg(1, 2));
        pol_n += out.struck.pol;
        ++n_pol_n;
        k_d.push_back(out.k);
      }
    }
    const double f_p = static_cast<double>(n_p) / n;
    MESSAGE("triton m_S = " << m_s << ": struck p in " << f_p
                            << " of breakups, <pol_p> = " << pol_p / n_pol_p
                            << " (want " << 2 * m_s * t.eff_pol_p
                            << "), <pol_n> = " << pol_n / n_pol_n << " (want "
                            << 2 * m_s * t.eff_pol_n << "), <k(n+d)> = "
                            << mean(k_d) << " GeV");
    CHECK(std::fabs(pol_p / n_pol_p - 2 * m_s * t.eff_pol_p) < 0.01);
    CHECK(std::fabs(pol_n / n_pol_n - 2 * m_s * t.eff_pol_n) < 0.01);
  }
}

TEST_CASE("breakup: the nn remnant's invariant mass is 2 sqrt(m_n^2 + q^2) "
          "and it sits at the virtual-state scale") {
  ClusterBreakup br;
  const double m_t = nuclear_mass(1, 3), m_n = nuclear_mass(0, 1);
  Rng rng(kSeed, 5, 0, 0);
  BreakupInput in;
  in.species = ClusterSpecies::Triton;
  in.p_cluster = at_rest(m_t);
  in.z = 1; in.a = 3;
  in.m_s = 0.5; in.s_cluster = 0.5;
  in.x = 0.1; in.q2 = 10.0;
  BreakupResult out;
  std::vector<double> q;
  double worst = 0.0;
  for (int i = 0; i < 40000 && q.size() < 5000; ++i) {
    REQUIRE(br.resolve(in, rng, out));
    if (out.partners.size() != 2) continue;
    const Vec4 nn = out.partners[0].p + out.partners[1].p;
    worst = std::max(worst, std::fabs(std::sqrt(nn.m2()) - out.m_remnant));
    CHECK_CLOSE(out.m_remnant,
                2.0 * std::sqrt(m_n * m_n + out.q_nn * out.q_nn), 1e-12);
    q.push_back(out.q_nn);
  }
  MESSAGE("nn relative momentum: mean " << mean(q) << " GeV, median "
          << quantile(q, 0.5) << " GeV (kappa_nn = " << KAPPA_NN_VIRTUAL
          << "); worst |sqrt(M_nn^2) - m_remnant| = " << worst << " GeV");
  CHECK(worst < 1e-9);
  // The virtual state is very shallow, so the pair comes out nearly at rest
  // relative to itself: the median q is well below a typical Fermi momentum.
  CHECK(quantile(q, 0.5) < 0.1);
}

// ---------------------------------------------------------------------------
// The pipeline at T1.

TEST_CASE("T1 pipeline: whole-record 4-momentum and charge conservation on "
          "every tagged channel") {
  for (PipelineChannel ch : {PipelineChannel::TaggedLi6Alpha,
                             PipelineChannel::TaggedLi7Alpha,
                             PipelineChannel::TaggedDeuteronP}) {
    PipelineConfig cfg = tagged_cfg(ch, 20000);
    CHECK(cfg.tier == Tier::T1);          // the default
    Pipeline p(cfg, plan_for_channel(ch));
    REQUIRE(p.tier() == Tier::T1);

    double worst_p = 0.0, worst_q = 0.0, worst_rel = 0.0;
    double sum_virt = 0.0, min_virt = 0.0;
    std::size_t n_partner = 0, n_ev = 0;
    p.for_each([&](const Event& ev) {
      ++n_ev;
      const Vec4 r = momentum_residual(ev);
      const Vec4 sc = momentum_scale(ev);
      worst_p = std::max(worst_p, std::max({std::fabs(r.e), std::fabs(r.px),
                                            std::fabs(r.py), std::fabs(r.pz)}));
      worst_rel = std::max(worst_rel, std::fabs(r.e) / sc.e);
      worst_q = std::max(worst_q, std::fabs(charge_residual(ev)));

      // T1 names the struck nucleon and every partner spectator.
      const Particle* pn = ev.find(Role::StruckNucleon);
      REQUIRE(pn != nullptr);
      REQUIRE(pn->status == Status::Intermediate);
      CHECK((pn->pdg == 2212 || pn->pdg == 2112));
      CHECK((pn->pol == 1.0 || pn->pol == -1.0));
      const double virt = pn->p.m2() - M_NUCLEON * M_NUCLEON;
      sum_virt += virt;
      min_virt = std::min(min_virt, virt);

      // The struck cluster's own four-vector is exactly the struck nucleon
      // plus its partners -- the T1 tier ADDS to the T0 record, it does not
      // move anything already in it.
      const Particle* pc = ev.find(Role::StruckCluster);
      REQUIRE(pc != nullptr);
      Vec4 s = pn->p;
      std::size_t np = 0;
      for (const Particle& q : ev.particles) {
        if (q.role != Role::PartnerSpectator) continue;
        REQUIRE(q.status == Status::Final);
        s = s + q.p;
        ++np;
      }
      n_partner += np;
      CHECK(std::fabs(s.e - pc->p.e) < 1e-9);
      CHECK(std::fabs(s.pz - pc->p.pz) < 1e-9);
      CHECK(std::fabs(s.px - pc->p.px) < 1e-9);
      CHECK(std::fabs(s.py - pc->p.py) < 1e-9);
      // X is the PER-NUCLEON remainder now.
      const Particle* xx = ev.find(Role::HadronicX);
      REQUIRE(xx != nullptr);
      const Vec4 want = ev.particles[0].p + pn->p -
                        ev.find(Role::ScatteredElectron)->p;
      CHECK(std::fabs(xx->p.e - want.e) < 1e-9);
      CHECK(xx->p.m2() >= 0.0);
    });

    MESSAGE(std::string(pipeline_channel_name(ch))
            << " T1: worst |dp| = " << worst_p << " GeV (rel " << worst_rel
            << "), worst |dq| = " << worst_q << ", <virtuality> = "
            << sum_virt / n_ev << " GeV^2, min " << min_virt
            << " GeV^2, partners/event = "
            << static_cast<double>(n_partner) / n_ev);
    CHECK(worst_p < 1e-9);
    CHECK(worst_rel < 1e-12);
    CHECK(worst_q == 0.0);
    // Off shell, and by a Fermi-motion amount, not a nuclear-mass one.
    CHECK(sum_virt / n_ev < 0.0);
    CHECK(sum_virt / n_ev > -1.0);
  }
}

TEST_CASE("T1 pipeline: partner multiplicity per channel") {
  struct Row { PipelineChannel ch; std::size_t lo, hi; };
  const Row rows[] = {
      {PipelineChannel::TaggedLi6Alpha, 1, 1},      // d* -> N + N
      {PipelineChannel::TaggedLi7Alpha, 1, 2},      // t* -> n + d  /  p + nn
      {PipelineChannel::TaggedDeuteronP, 0, 0},     // already a nucleon
  };
  for (const Row& r : rows) {
    Pipeline p(tagged_cfg(r.ch, 3000), plan_for_channel(r.ch));
    std::size_t lo = 99, hi = 0;
    p.for_each([&](const Event& ev) {
      std::size_t n = 0;
      for (const Particle& q : ev.particles)
        if (q.role == Role::PartnerSpectator) ++n;
      lo = std::min(lo, n);
      hi = std::max(hi, n);
    });
    MESSAGE(std::string(pipeline_channel_name(r.ch)) << ": partners per event in [" << lo
                                        << ", " << hi << "]");
    CHECK(lo == r.lo);
    CHECK(hi == r.hi);
  }
}

TEST_CASE("T1 pipeline: T0 is unchanged by the tier switch, and T1 adds only "
          "the partners and the struck nucleon") {
  PipelineConfig c0 = tagged_cfg(PipelineChannel::TaggedLi6Alpha, 2000);
  c0.tier = Tier::T0;
  PipelineConfig c1 = tagged_cfg(PipelineChannel::TaggedLi6Alpha, 2000);
  Pipeline p0(c0, tensor_thirds_plan(0.7, 0.6));
  Pipeline p1(c1, tensor_thirds_plan(0.7, 0.6));
  REQUIRE(p0.tier() == Tier::T0);
  REQUIRE(p0.breakup() == nullptr);
  REQUIRE(p1.breakup() != nullptr);
  REQUIRE(p0.size() == p1.size());
  CHECK_CLOSE(p0.sigma_pb(), p1.sigma_pb(), 1e-15);

  // The breakup consumes the event's stream only AFTER every T0 quantity is
  // fixed, so the two tiers give the SAME T0 record -- except at an index
  // where T1's own timelike-X rejection fires and T0's does not.  X shrinks
  // from k + P_X - k' to the per-nucleon k + p_N - k', so a draw that was
  // fine at T0 can leave a spacelike remainder here; measured 0.024 % of
  // 6Li-alpha draws (0.022 % 7Li-alpha, 0.004 % d control) over 50 k events.
  std::size_t n_redrawn = 0;
  for (std::uint64_t i = 0; i < 2000; ++i) {
    const Event a = p0.event(i), b = p1.event(i);
    CHECK(a.find(Role::StruckNucleon) == nullptr);
    CHECK(b.find(Role::StruckNucleon) != nullptr);
    if (a.kin.x != b.kin.x) { ++n_redrawn; continue; }
    CHECK(a.kin.q2 == b.kin.q2);
    CHECK(a.kin.k == b.kin.k);
    CHECK(a.spin.m_struck == b.spin.m_struck);
    const Particle* sa = a.find(Role::Spectator);
    const Particle* sb = b.find(Role::Spectator);
    REQUIRE(sa != nullptr);
    REQUIRE(sb != nullptr);
    CHECK(sa->p.e == sb->p.e);
    CHECK(sa->p.pz == sb->p.pz);
    CHECK(a.find(Role::StruckCluster)->p.e == b.find(Role::StruckCluster)->p.e);
  }
  MESSAGE("T1 timelike-X redraws: " << n_redrawn << " of 2000 indices");
  CHECK(n_redrawn * 100 < 2000);          // < 1 %
}

TEST_CASE("T1 pipeline: determinism -- same index, same record, any order") {
  for (PipelineChannel ch : {PipelineChannel::TaggedLi6Alpha,
                             PipelineChannel::TaggedLi7Alpha}) {
    Pipeline p(tagged_cfg(ch, 500), plan_for_channel(ch));
    for (std::uint64_t i : {0ULL, 137ULL, 499ULL}) {
      const Event a = p.event(i);
      const Event b = p.event(i);
      REQUIRE(a.particles.size() == b.particles.size());
      for (std::size_t j = 0; j < a.particles.size(); ++j) {
        CHECK(a.particles[j].pdg == b.particles[j].pdg);
        CHECK(a.particles[j].pol == b.particles[j].pol);
        CHECK(a.particles[j].p.e == b.particles[j].p.e);
        CHECK(a.particles[j].p.px == b.particles[j].p.px);
        CHECK(a.particles[j].p.py == b.particles[j].p.py);
        CHECK(a.particles[j].p.pz == b.particles[j].p.pz);
      }
    }
    // and threaded generation is bit identical to serial
    std::vector<double> serial, threaded;
    p.for_each([&](const Event& ev) {
      const Particle* pn = ev.find(Role::StruckNucleon);
      serial.push_back(pn->p.pz + pn->pol);
    });
    p.for_each([&](const Event& ev) {
      const Particle* pn = ev.find(Role::StruckNucleon);
      threaded.push_back(pn->p.pz + pn->pol);
    }, 4, 64);
    REQUIRE(serial.size() == threaded.size());
    CHECK(serial == threaded);
  }
}

TEST_CASE("T1 pipeline: kin.cell is the sampler's own accepted-cell index") {
  // Inclusive: the ION-level sampler.
  PipelineConfig ci;
  ci.channel = PipelineChannel::Inclusive;
  ci.isotope = "6Li";
  ci.beam_config = 1;
  ci.n_events = 2000;
  ci.seed = kSeed;
  Pipeline pi(ci, tensor_thirds_plan(0.7, 0.6));
  const std::size_t n_cells_i = pi.dis_sampler().n_cells();
  std::map<int, int> seen;
  pi.for_each([&](const Event& ev) {
    REQUIRE(ev.kin.cell >= 0);
    REQUIRE(static_cast<std::size_t>(ev.kin.cell) < n_cells_i);
    ++seen[ev.kin.cell];
    // the cell the event was drawn in contains its (x, Q2)
    const std::size_t c = static_cast<std::size_t>(ev.kin.cell);
    CHECK(std::log(ev.kin.x) >= pi.dis_sampler().logx_lo()[c] - 1e-12);
    CHECK(std::log(ev.kin.x) <= pi.dis_sampler().logx_hi()[c] + 1e-12);
    CHECK(std::log(ev.kin.q2) >= pi.dis_sampler().logq2_lo()[c] - 1e-12);
    CHECK(std::log(ev.kin.q2) <= pi.dis_sampler().logq2_hi()[c] + 1e-12);
  });
  MESSAGE("inclusive: " << seen.size() << " distinct cells of " << n_cells_i
                        << " over 2000 events");
  CHECK(seen.size() > 20);

  // Tagged: the STRUCK-CLUSTER sampler.
  Pipeline pt(tagged_cfg(PipelineChannel::TaggedLi6Alpha, 2000),
              tensor_thirds_plan(0.7, 0.6));
  const std::size_t n_cells_t = pt.dis_sampler().n_cells();
  std::map<int, int> seen_t;
  pt.for_each([&](const Event& ev) {
    REQUIRE(ev.kin.cell >= 0);
    REQUIRE(static_cast<std::size_t>(ev.kin.cell) < n_cells_t);
    ++seen_t[ev.kin.cell];
    const std::size_t c = static_cast<std::size_t>(ev.kin.cell);
    CHECK(std::log(ev.kin.x) >= pt.dis_sampler().logx_lo()[c] - 1e-12);
    CHECK(std::log(ev.kin.x) <= pt.dis_sampler().logx_hi()[c] + 1e-12);
  });
  MESSAGE("tagged 6Li-alpha: " << seen_t.size() << " distinct cells of "
                               << n_cells_t << " over 2000 events");
  CHECK(seen_t.size() > 20);
}

TEST_CASE("T1 pipeline: the stored spectator-lab block IS boost_spectator's "
          "own answer") {
  for (PipelineChannel ch : {PipelineChannel::TaggedLi6Alpha,
                             PipelineChannel::TaggedLi7Alpha,
                             PipelineChannel::TaggedDeuteronP}) {
    Pipeline p(tagged_cfg(ch, 3000), plan_for_channel(ch));
    double worst = 0.0;
    p.for_each([&](const Event& ev) {
      const Particle* f = ev.find(Role::Spectator);
      const Particle* bi = ev.find(Role::BeamIon);
      REQUIRE(f != nullptr);
      REQUIRE(bi != nullptr);
      // re-derive exactly as `export.columns_from_events` does
      const double pt = f->p.pt(), plab = f->p.p();
      const double gamma = bi->p.e / bi->mass, gbeta = bi->p.pz / bi->mass;
      worst = std::max(worst, std::fabs(ev.kin.spec_pt - pt));
      worst = std::max(worst, std::fabs(ev.kin.spec_p_lab - plab));
      worst = std::max(worst,
                       std::fabs(ev.kin.spec_theta - std::atan2(pt, f->p.pz)));
      worst = std::max(worst, std::fabs(ev.kin.spec_kx - f->p.px));
      worst = std::max(worst, std::fabs(ev.kin.spec_ky - f->p.py));
      worst = std::max(worst, std::fabs(ev.kin.spec_kz -
                                        (gamma * f->p.pz - gbeta * f->p.e)));
      worst = std::max(worst, std::fabs(ev.kin.phi_spec -
                                        std::atan2(f->p.py, f->p.px)));
      if (f->charge > 0.0) {
        worst = std::max(worst,
                         std::fabs(ev.kin.spec_r - (plab / f->charge) /
                                                       (bi->p.pz / bi->charge)));
      }
    });
    MESSAGE(std::string(pipeline_channel_name(ch))
            << ": worst |stored - re-derived| spectator-lab component = "
            << worst);
    CHECK(worst < 1e-9);
  }
}
