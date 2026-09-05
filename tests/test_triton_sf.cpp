// SPDX-License-Identifier: GPL-3.0-or-later
// The triton spectral function (triton_sf.hpp): the Ciofi degli Atti-Simula
// transcription and its three-channel breakup.  Every number the header
// quotes as "MEASURED HERE" is measured here; the conservation and RNG
// discipline the model promises are checked through the T1 pipeline.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "check_close.hpp"
#include "doctest.h"

#include "lipolgen/pipeline.hpp"
#include "lipolgen/triton_sf.hpp"

using namespace lipolgen;

namespace {

constexpr std::uint64_t kSeed = 20260901;

}  // namespace

// ---------------------------------------------------------------------------

TEST_CASE("triton sf: the CS spectroscopic factors -- S_0 = 0.6525 for A = 3, "
          "and S_0 + S_1 closes on 1 untuned") {
  // A = 3: S_0 is the 3He(e,e'p)d spectroscopic factor, the ~2/3 the
  // literature quotes; S_0 + S_1 = 1 is CS Eq. (28) and nothing in the
  // transcription was fitted to make it come out.
  const double s0 = cs_norm(cs_n0_terms(3));
  const double s1 = cs_norm(cs_n1_terms(3));
  MESSAGE("A = 3: S_0 = " << s0 << ", S_1 = " << s1 << ", sum " << s0 + s1);
  CHECK_CLOSE(s0, 0.6525, 1e-3);           // the spec number
  CHECK_CLOSE(s0, 0.652548, 1e-4);         // as transcribed, tighter
  CHECK_CLOSE(s1, 0.347113, 1e-3);
  CHECK(std::fabs(s0 + s1 - 1.0) < 5e-4);  // header: 0.99966

  // A = 2: the deuteron's n_0 IS the whole distribution (S_0 = 1 up to the
  // parameterization's own 0.3 %) and it has no excitable remnant at all.
  CHECK_CLOSE(cs_norm(cs_n0_terms(2)), 1.0031, 1e-3);
  CHECK(cs_n1_terms(2).empty());
  CHECK(cs_norm(cs_n1_terms(2)) == 0.0);

  // A = 4, carried for the cross-check: S_0 ~ 0.8, sum closes to 2.2e-3.
  const double s0_4 = cs_norm(cs_n0_terms(4));
  const double s1_4 = cs_norm(cs_n1_terms(4));
  CHECK_CLOSE(s0_4, 0.799576, 1e-3);
  CHECK_CLOSE(s1_4, 0.198250, 1e-3);
  CHECK_CLOSE(s0_4 + s1_4, 0.99783, 1e-3);

  // CS give no other nucleus.
  CHECK_THROWS_AS(cs_n0_terms(5), std::runtime_error);
  CHECK_THROWS_AS(cs_n1_terms(1), std::runtime_error);
}

TEST_CASE("triton sf: the n_0 high-k tail table of the A = 3 set") {
  const std::vector<CsTerm>& n0 = cs_n0_terms(3);
  // P(|k| > k0) of k^2 n_0(k), k0 in GeV -- the four spec rows.
  CHECK_CLOSE(cs_p_above(n0, 0.10), 0.4266, 1e-3);
  CHECK_CLOSE(cs_p_above(n0, 0.20), 0.06276, 1e-3);
  CHECK_CLOSE(cs_p_above(n0, 0.30), 0.009621, 1e-3);
  CHECK_CLOSE(cs_p_above(n0, 0.45), 0.002365, 1e-3);
  // Sanity of the helper itself: P(> 0) = 1, P(> ceiling) = 0, monotone.
  CHECK_CLOSE(cs_p_above(n0, 0.0), 1.0, 1e-12);
  CHECK(cs_p_above(n0, 10.0 * HBARC_GEV_FM) == 0.0);
}

TEST_CASE("triton sf: mean k -- ~102 MeV on the bound-remnant channel, "
          "~126 MeV on the total (printed)") {
  const double mean_n0 = cs_mean_k(cs_n0_terms(3));
  // The total n_0 + n_1 is one CsTerm list too: same functional family.
  std::vector<CsTerm> total = cs_n0_terms(3);
  const std::vector<CsTerm>& n1 = cs_n1_terms(3);
  total.insert(total.end(), n1.begin(), n1.end());
  const double mean_tot = cs_mean_k(total);
  MESSAGE("A = 3 <k>: " << 1e3 * mean_n0 << " MeV on n_0 (the n + d channel), "
          << 1e3 * mean_tot << " MeV on n_0 + n_1 (all struck neutrons); the "
          "sequential model's were 133 (n + d) / 145 (p + nn)");
  CHECK_CLOSE(mean_n0, 0.1021, 2e-3);
  CHECK_CLOSE(mean_tot, 0.1260, 2e-3);
}

TEST_CASE("triton sf: the sampled branching integrates to S_0 = 0.653 +- "
          "0.005 at 1e5 draws, and p_two_body is the n_0/(n_0 + n_1) ratio") {
  const CiofiSimulaTriton sf;
  // The object's own S_0/S_1 (integrated over [0, k_max], where it samples).
  CHECK_CLOSE(sf.s0(), 0.6525, 1e-3);
  CHECK_CLOSE(sf.s0() + sf.s1(), 1.0, 5e-4);

  Rng rng(12345, 1, 0, 0);
  const int n = 100000;
  int n_bound = 0;
  double sum_k_bound = 0.0, sum_k = 0.0;
  for (int i = 0; i < n; ++i) {
    const TritonDraw d = sf.sample(0.0, rng);   // every draw a struck neutron
    REQUIRE(!d.struck_proton);
    sum_k += d.k;
    if (d.channel == TritonChannel::NeutronD) {
      ++n_bound;
      sum_k_bound += d.k;
      // The bound remnant: the deuteron, E = 0, no internal continuum.
      CHECK(d.e_rel == 0.0);
      CHECK(d.q_pair == 0.0);
      CHECK_CLOSE(d.m_remnant, nuclear_mass(1, 2), 1e-12);
    } else {
      REQUIRE(d.channel == TritonChannel::NeutronPnCont);
      // The continuum pn pair: M = sum sqrt(m^2 + q^2), E = M - m_1 - m_2.
      CHECK(d.e_rel > 0.0);
      CHECK(d.q_pair > 0.0);
      const double m1 = nuclear_mass(1, 1), m2 = nuclear_mass(0, 1);
      const double q2 = d.q_pair * d.q_pair;
      CHECK_CLOSE(d.m_remnant,
                  std::sqrt(m1 * m1 + q2) + std::sqrt(m2 * m2 + q2), 1e-12);
    }
  }
  const double frac = static_cast<double>(n_bound) / n;
  MESSAGE("bound (n + d) fraction " << frac << " of " << n
          << " struck neutrons (S_0/(S_0 + S_1) = "
          << sf.s0() / (sf.s0() + sf.s1()) << "); <k> "
          << 1e3 * sum_k_bound / n_bound << " MeV bound / " << 1e3 * sum_k / n
          << " MeV overall");
  CHECK(std::fabs(frac - 0.653) < 0.005);                  // the spec window
  CHECK(std::fabs(sum_k_bound / n_bound - 0.102) < 0.003); // ~102 MeV
  CHECK(std::fabs(sum_k / n - 0.126) < 0.003);             // ~126 MeV

  // The branching is a RATIO the parameterization already contains -- and a
  // struck proton NEVER leaves a bound remnant (no bound nn exists).
  for (double k : {0.0, 0.05, 0.10, 0.20, 0.30, 0.45, 0.80}) {
    CHECK(sf.p_two_body(k, 2212) == 0.0);
    const double n0 = sf.n0_cs(k), n1 = sf.n1_cs(k);
    CHECK(sf.p_two_body(k, 2112) == n0 / (n0 + n1));
  }
  // ... and it FALLS with k: the correlated continuum owns the tail.
  CHECK(sf.p_two_body(0.45, 2112) < sf.p_two_body(0.10, 2112));
  CHECK(sf.p_two_body(0.10, 2112) < sf.p_two_body(0.0, 2112));
}

TEST_CASE("triton sf: RNG stream discipline -- exactly 6 uniforms per "
          "sample() on every branch -- and determinism") {
  const CiofiSimulaTriton sf;

  // Two identical streams; one runs 1000 mixed-branch samples (alternating
  // all-proton / all-neutron, so every channel is exercised), the other burns
  // exactly 6 * 1000 raw uniforms.  If any branch consumed a different count
  // the two streams would land in different places.
  Rng a(kSeed, 1, 0, 0), b(kSeed, 1, 0, 0);
  for (int i = 0; i < 1000; ++i) {
    (void)sf.sample((i & 1) ? 1.0 : 0.0, a);
  }
  for (int i = 0; i < 6 * 1000; ++i) (void)b.uniform();
  CHECK(a.uniform() == b.uniform());
  CHECK(TritonSpectralFunction::kUniformsPerSample == 6);

  // Determinism: the same counter-based stream gives the identical
  // TritonDraw sequence, field for field.
  Rng c(987, 2, 3, 4), d(987, 2, 3, 4);
  for (int i = 0; i < 500; ++i) {
    const TritonDraw x = sf.sample(0.3, c);
    const TritonDraw y = sf.sample(0.3, d);
    CHECK(x.channel == y.channel);
    CHECK(x.struck_proton == y.struck_proton);
    CHECK(x.k == y.k);
    CHECK(x.cos_theta_k == y.cos_theta_k);
    CHECK(x.phi_k == y.phi_k);
    CHECK(x.e_rel == y.e_rel);
    CHECK(x.q_pair == y.q_pair);
    CHECK(x.m_remnant == y.m_remnant);
  }
}

TEST_CASE("triton sf: T1 pipeline with --triton-sf ciofi-simula -- "
          "whole-record conservation, the three channels, and determinism") {
  PipelineConfig cfg;
  cfg.channel = PipelineChannel::TaggedLi7Alpha;
  cfg.isotope = channel_isotope(cfg.channel);
  cfg.beam_config = 1;
  cfg.n_events = 20000;
  cfg.seed = kSeed;
  cfg.optics_choice = OpticsChoice::Tagging;
  cfg.triton_sf = TritonSfChoice::CiofiSimula;
  REQUIRE(cfg.tier == Tier::T1);
  Pipeline p(cfg, helicity_flip_plan(1.5, 0.7, 0.7));

  double worst_p = 0.0, worst_q = 0.0;
  // The channel, read off the record: struck n + ONE partner (the d) is the
  // bound channel; struck n + TWO partners is the NEW n + (pn) continuum,
  // recognizable by its proton partner; struck p always leaves two neutrons.
  std::size_t n_nd = 0, n_npn = 0, n_pnn = 0, n_ev = 0;
  std::vector<double> first_pz;
  p.for_each([&](const Event& ev) {
    ++n_ev;
    const Vec4 r = momentum_residual(ev);
    worst_p = std::max(worst_p, std::max({std::fabs(r.e), std::fabs(r.px),
                                          std::fabs(r.py), std::fabs(r.pz)}));
    worst_q = std::max(worst_q, std::fabs(charge_residual(ev)));

    const Particle* pn = ev.find(Role::StruckNucleon);
    REQUIRE(pn != nullptr);
    if (first_pz.size() < 64) first_pz.push_back(pn->p.pz);
    std::vector<int> partner_pdg;
    for (const Particle& q : ev.particles) {
      if (q.role == Role::PartnerSpectator) partner_pdg.push_back(q.pdg);
    }
    if (pn->pdg == 2112 && partner_pdg.size() == 1) {
      CHECK(partner_pdg[0] == nuclide_pdg(1, 2));   // the deuteron
      ++n_nd;
    } else if (pn->pdg == 2112) {
      REQUIRE(partner_pdg.size() == 2);             // the NEW third channel
      const int n_prot = static_cast<int>(
          std::count(partner_pdg.begin(), partner_pdg.end(), 2212));
      CHECK(n_prot == 1);                           // one p + one n
      ++n_npn;
    } else {
      REQUIRE(pn->pdg == 2212);
      REQUIRE(partner_pdg.size() == 2);             // always the nn continuum
      CHECK(partner_pdg[0] == 2112);
      CHECK(partner_pdg[1] == 2112);
      ++n_pnn;
    }
  });
  REQUIRE(n_ev == cfg.n_events);

  const double frac_bound =
      static_cast<double>(n_nd) / static_cast<double>(n_nd + n_npn);
  MESSAGE("ciofi-simula T1: worst |dp| " << worst_p << " GeV, worst |dq| "
          << worst_q << "; channels n+d " << n_nd << " / n+(pn) " << n_npn
          << " / p+(nn) " << n_pnn << " -> bound fraction among struck "
          "neutrons " << frac_bound);
  // The spec's conservation bound, with the option ON.
  CHECK(worst_p < 1e-9);
  CHECK(worst_q == 0.0);
  // All three channels populated, and the n+d : n+(pn) split is S_0's
  // (0.653; ~12k struck neutrons -> sigma ~ 0.004, gate at 3.5 sigma).
  CHECK(n_pnn > 0);
  CHECK(n_npn > 0);
  CHECK(std::fabs(frac_bound - 0.653) < 0.015);

  SUBCASE("same seed, same records; the Hulthen default differs") {
    Pipeline p2(cfg, helicity_flip_plan(1.5, 0.7, 0.7));
    std::vector<double> again;
    p2.for_each_range([&](const Event& ev) {
      again.push_back(ev.find(Role::StruckNucleon)->p.pz);
    }, 0, 64);
    REQUIRE(again.size() == first_pz.size());
    for (std::size_t i = 0; i < again.size(); ++i)
      CHECK(again[i] == first_pz[i]);

    PipelineConfig hcfg = cfg;
    hcfg.triton_sf = TritonSfChoice::Hulthen;
    Pipeline ph(hcfg, helicity_flip_plan(1.5, 0.7, 0.7));
    std::size_t n_same = 0;
    std::size_t i = 0;
    ph.for_each_range([&](const Event& ev) {
      if (ev.find(Role::StruckNucleon)->p.pz == first_pz[i++]) ++n_same;
    }, 0, 64);
    CHECK(n_same < 8);   // a different model gives different records
  }
}
