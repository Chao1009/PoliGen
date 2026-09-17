// SPDX-License-Identifier: GPL-3.0-or-later
// Cluster-spectator kinematics and far-forward routing, against
// validation/reference/spectator.json (rtol 1e-12, the boost is closed form)
// and against the analytic identities fastsim/tests/test_spectator.py pins.

#include <fstream>
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "check_close.hpp"
#include "doctest.h"
#include "json_min.hpp"
#include "lipolgen/beams.hpp"
#include "lipolgen/rng.hpp"
#include "lipolgen/spectator.hpp"

using namespace lipolgen;

#ifndef LIPOLGEN_REFERENCE_DIR
#define LIPOLGEN_REFERENCE_DIR "validation/reference"
#endif

namespace {

constexpr double kRtol = 1e-12;

bool load_spectator(jsonmin::Value& out) {
  return jsonmin::load_file(std::string(LIPOLGEN_REFERENCE_DIR) + "/spectator.json",
                            out);
}

/// Is the reference blob on disk?  Called at doctest REGISTRATION time by the
/// `doctest::skip()` decorators below, so that a checkout without
/// `validation/reference/` reports these cases as SKIPPED in the tally
/// instead of PASSED WITH ZERO ASSERTIONS.  The pattern is
/// `tests/test_coherent.cpp`'s (phase F, 2026-09-05), which introduced it for
/// exactly this reason and was applied to test_coherent / test_cluster /
/// test_b1_nuclear and to no other file -- these twelve rtol-1e-12 port gates
/// kept the in-case `MESSAGE(...); return;` and went on passing, silently and
/// vacuously, until 2026-09-16.  A blob that IS there and does not parse is a
/// FAILURE, not a skip: that is what the `REQUIRE` on the loader inside each
/// case is for.
bool spectator_present() {
  return std::ifstream(std::string(LIPOLGEN_REFERENCE_DIR) +
                       "/spectator.json")
      .good();
}


double median(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  const std::size_t n = v.size();
  return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

double quantile(std::vector<double> v, double q) {
  std::sort(v.begin(), v.end());
  const double pos = q * (v.size() - 1);
  const std::size_t i = static_cast<std::size_t>(pos);
  const double f = pos - i;
  return i + 1 < v.size() ? v[i] * (1.0 - f) + v[i + 1] * f : v.back();
}

}  // namespace

TEST_CASE("spectator: the AME2020 mass tables are the Python's" *
          doctest::skip(!spectator_present())) {
  jsonmin::Value ref;
  REQUIRE(load_spectator(ref));
  const jsonmin::Value& c = ref["constants"];
  CHECK_CLOSE(M_U, c["M_U"].num(), kRtol);
  for (const auto& kv : c["MASSES"].obj()) {
    CHECK_CLOSE(cluster_mass(kv.first), kv.second.num(), kRtol);
  }
  for (const auto& kv : c["NUCLEUS_MASS"].obj()) {
    const std::size_t sep = kv.first.find('_');
    const int z = std::stoi(kv.first.substr(0, sep));
    const int a = std::stoi(kv.first.substr(sep + 1));
    CHECK(nuclear_mass_known(z, a));
    CHECK_CLOSE(nuclear_mass(z, a), kv.second.num(), kRtol);
  }
}

TEST_CASE("spectator: every channel's masses, kappa and R(k=0)" *
          doctest::skip(!spectator_present())) {
  jsonmin::Value ref;
  REQUIRE(load_spectator(ref));
  for (const auto& kv : ref["channels"].obj()) {
    const ClusterChannel& ch = channel_by_name(kv.first);
    const jsonmin::Value& r = kv.second;
    CAPTURE(kv.first);
    CHECK(ch.beam_A == static_cast<int>(r["beam_A"].num()));
    CHECK(ch.beam_Z == static_cast<int>(r["beam_Z"].num()));
    CHECK(ch.spectator == r["spectator"].str());
    CHECK(ch.spectator_A == static_cast<int>(r["spectator_A"].num()));
    CHECK(ch.spectator_Z == static_cast<int>(r["spectator_Z"].num()));
    CHECK(ch.l_wave == static_cast<int>(r["l_wave"].num()));
    CHECK_CLOSE(ch.separation_energy, r["separation_energy"].num(), kRtol);
    CHECK_CLOSE(ch.m_spec(), r["m_spec"].num(), kRtol);
    CHECK_CLOSE(ch.m_beam(), r["m_beam"].num(), kRtol);
    CHECK_CLOSE(ch.m_partner(), r["m_partner"].num(), kRtol);
    CHECK_CLOSE(ch.kappa(), r["kappa"].num(), kRtol);
    if (ch.spectator_Z > 0) {
      CHECK_CLOSE(ch.r_at_k_zero(), r["R_at_k_zero"].num(), kRtol);
    } else {
      CHECK(std::isnan(ch.r_at_k_zero()));
    }
  }
}

TEST_CASE("spectator: the unnormalized momentum densities on the k grid" *
          doctest::skip(!spectator_present())) {
  jsonmin::Value ref;
  REQUIRE(load_spectator(ref));
  for (const auto& kv : ref["channels"].obj()) {
    const ClusterChannel& ch = channel_by_name(kv.first);
    const jsonmin::Value& md = kv.second["momentum_density"];
    const double beta = md["beta"].num();
    const std::vector<double> kg = md["k_grid"].flat();
    const std::vector<double> nk = md["n_of_k"].flat();
    CAPTURE(kv.first);
    for (std::size_t i = 0; i < kg.size(); ++i) {
      CHECK_CLOSE(momentum_density(kg[i], ch.kappa(), beta, ch.l_wave), nk[i],
                  kRtol);
    }
  }
}

TEST_CASE("spectator: _boost_fragment at fixed (kx, ky, kz)" *
          doctest::skip(!spectator_present())) {
  jsonmin::Value ref;
  REQUIRE(load_spectator(ref));
  for (const auto& kv : ref["channels"].obj()) {
    const ClusterChannel& ch = channel_by_name(kv.first);
    CAPTURE(kv.first);
    for (const jsonmin::Value& row : kv.second["boost_fragment"].arr()) {
      const FragmentLab lab = boost_spectator_fragment(
          ch, row["p_per_nucleon"].num(), row["kx"].num(), row["ky"].num(),
          row["kz"].num());
      CHECK_CLOSE(lab.pT, row["pT"].num(), kRtol);
      CHECK_CLOSE(lab.theta, row["theta"].num(), kRtol);
      CHECK_CLOSE(lab.phi, row["phi"].num(), kRtol);
      CHECK_CLOSE(lab.p_lab, row["p_lab"].num(), kRtol);
      CHECK_CLOSE(lab.xL, row["xL"].num(), kRtol);
      CHECK_CLOSE(lab.k, row["k"].num(), kRtol);
      if (ch.spectator_Z > 0) {
        CHECK_CLOSE(lab.R, row["R"].num(), kRtol);
      } else {
        CHECK(std::isnan(lab.R));
      }
    }
  }
}

TEST_CASE("spectator: kappa is a property of the beam, not of the tag") {
  // mu is the reduced mass of the two SEPARATED clusters and is symmetric
  // under swapping them, so the alpha-tag and the partner-tag of one beam
  // must share kappa.  With m_beam = A*M_U and m_partner = m_beam - m_spec
  // they did not (60.50 vs 60.62 MeV for 6Li) -- the symptom that fixed the
  // construction.
  CHECK_CLOSE(LI6_ALPHA_TAG().kappa(), LI6_D_TAG().kappa(), 1e-5);
  CHECK_CLOSE(LI7_ALPHA_TAG().kappa(), LI7_T_TAG().kappa(), 1e-5);
  CHECK_CLOSE_AT(LI6_ALPHA_TAG().m_partner(), cluster_mass("d"), 0.0, 5e-6);
  CHECK_CLOSE_AT(LI6_ALPHA_TAG().m_beam() - LI6_ALPHA_TAG().m_spec(),
                 cluster_mass("d") - 1.4743e-3, 0.0, 5e-6);
  // external anchor: the d -> p+n control is the textbook deuteron wave
  // number, 1/kappa = 4.318 fm (hbar c = 0.1973269804 GeV fm)
  CHECK_CLOSE_AT(DEUTERON_P_TAG().kappa(), 0.04570, 0.0, 5e-5);
  CHECK_CLOSE_AT(0.1973269804 / DEUTERON_P_TAG().kappa(), 4.318, 0.0, 5e-3);
}

TEST_CASE("spectator: the A*M_U fallback off the mass table") {
  // u is one twelfth of the mass of a neutral 12C ATOM, so at (Z, A) = (6, 12)
  // the fallback is HIGH of the nuclear mass -- the opposite sign to every
  // entry the table does carry.
  const double u = 931.49410242e-3, m_e = 0.51099895e-3, b_e6 = 1030.11e-9;
  const double m_c12 = 12 * u - 6 * m_e + b_e6;
  CHECK(!nuclear_mass_known(6, 12));
  CHECK_CLOSE(nuclear_mass(6, 12), 12 * M_U, kRtol);
  CHECK(nuclear_mass(6, 12) > m_c12);
  CHECK_CLOSE_AT(nuclear_mass(6, 12) - m_c12, 3.016e-3, 0.0, 2e-6);
  // the m_partner fallback keeps the DEFINITION S = m_spec + m_partner - m_beam
  const double s_alpha = 7.36659e-3;
  const ClusterChannel c12{"12C -> alpha + 8Be", 12, 6, "alpha", 4, 2, s_alpha, 0};
  CHECK(!nuclear_mass_known(4, 8));
  CHECK_CLOSE_AT(c12.m_spec() + c12.m_partner() - c12.m_beam(), s_alpha, 0.0, 1e-12);
}

TEST_CASE("spectator: R(k=0) is a ratio of mass-to-charge ratios") {
  struct Case { const ClusterChannel& ch; double p_u; double drop; };
  const Case cases[] = {{LI6_ALPHA_TAG(), 137.5, 1.87e-3},
                        {LI7_ALPHA_TAG(), 117.9, 1.67e-3}};
  for (const Case& c : cases) {
    const FragmentLab at_rest = boost_spectator_fragment(c.ch, c.p_u, 0, 0, 0);
    CHECK_CLOSE(at_rest.R, c.ch.r_at_k_zero(), kRtol);
    const double naive = (c.ch.spectator_A * static_cast<double>(c.ch.beam_Z))
                         / (c.ch.beam_A * static_cast<double>(c.ch.spectator_Z));
    CHECK(at_rest.R < naive);
    CHECK_CLOSE((naive - at_rest.R) / naive, c.drop, 0.02);
  }
  // the headline numbers
  CHECK_CLOSE_AT(LI6_ALPHA_TAG().r_at_k_zero(), 0.99813, 0.0, 5e-6);
  CHECK_CLOSE_AT(LI7_ALPHA_TAG().r_at_k_zero(), 0.85571, 0.0, 5e-6);
}

TEST_CASE("spectator: sampled distributions sit on the analytic centres") {
  // The Python's PCG64 stream cannot be matched, so this compares
  // DISTRIBUTIONS, not draws.
  const MomentumSampler s6(LI6_ALPHA_TAG());
  const MomentumSampler s7(LI7_ALPHA_TAG());
  std::vector<double> r6, r7, pt6, th6;
  Rng rng(12345, 0, 0, 0);
  const int n = 50000;
  for (int i = 0; i < n; ++i) {
    double kx, ky, kz;
    s6.sample(rng, kx, ky, kz);
    const FragmentLab a = boost_spectator_fragment(LI6_ALPHA_TAG(), 137.5, kx, ky, kz);
    r6.push_back(a.R);
    pt6.push_back(a.pT);
    th6.push_back(a.theta);
    s7.sample(rng, kx, ky, kz);
    const FragmentLab b = boost_spectator_fragment(LI7_ALPHA_TAG(), 117.9, kx, ky, kz);
    r7.push_back(b.R);
  }
  CHECK_CLOSE_AT(median(r6), 1.0, 0.0, 0.01);
  CHECK_CLOSE_AT(median(r7), 6.0 / 7.0, 0.0, 0.01);
  CHECK(median(r6) > LI6_ALPHA_TAG().r_at_k_zero());  // the smearing is one-sided
  CHECK(median(r6) < 1.0);                            // 0.9984, not the naive 1
  CHECK(median(pt6) > 0.02);
  CHECK(median(pt6) < 0.3);
  CHECK(quantile(th6, 0.99) < 5e-3);
}

TEST_CASE("spectator: the lab boost is by the beam velocity of the NUCLEAR mass") {
  // Rapidity is additive under a longitudinal boost and k is drawn
  // isotropically, so <y_spec> = y_beam with no Fermi-smearing bias at all.
  // A*M_U would put the mean 2.0e-3 to 6.8e-3 above it.
  struct Case { const ClusterChannel& ch; double p_u; };
  const Case cases[] = {{LI6_ALPHA_TAG(), 137.5}, {LI6_D_TAG(), 137.5},
                        {LI7_ALPHA_TAG(), 117.9}, {LI7_T_TAG(), 117.9},
                        {DEUTERON_P_TAG(), 130.0}, {HE3_P_TAG(), 166.0}};
  for (const Case& c : cases) {
    CAPTURE(c.ch.name);
    const MomentumSampler s(c.ch);
    const double p_beam = c.ch.beam_A * c.p_u;
    const double y_beam = std::asinh(p_beam / c.ch.m_beam());
    const double y_amu = std::asinh(p_beam / (c.ch.beam_A * M_U));
    Rng rng(777, 1, 0, 0);
    double sum = 0.0;
    const int n = 200000;
    for (int i = 0; i < n; ++i) {
      double kx, ky, kz;
      s.sample(rng, kx, ky, kz);
      const FragmentLab lab = boost_spectator_fragment(c.ch, c.p_u, kx, ky, kz);
      const double e = std::hypot(lab.p_lab, c.ch.m_spec());
      sum += 0.5 * std::log((e + lab.pz_lab) / (e - lab.pz_lab));
    }
    const double resid = sum / n - y_beam;
    CHECK(std::fabs(resid) < 0.25 * (y_amu - y_beam));
    CHECK(std::fabs(resid) < 5e-4);
  }
}

TEST_CASE("farforward: the routing windows") {
  const Optics& o = HIGH_ACCEPTANCE();
  const double env = o.n_sigma * o.sigma_theta;   // 0.727 mrad
  const double kNoPhi = std::nan("");
  // 7Li alpha: R = 0.857, off rigidity -> Roman Pots, no near-beam cut
  CHECK(route_charged(0.857, 1e-3, 0.01, o, kNoPhi) == kRouteRomanPots);
  // 6Li alpha: R = 1 -> the ANGULAR envelope decides, not pT
  CHECK(route_charged(1.0, 1e-3, 0.05, o, kNoPhi) == kRouteRPNearBeam);
  CHECK(route_charged(1.0, 0.5 * env, 0.30, o, kNoPhi) == kRouteLost);
  CHECK(route_charged(1.0, 1.5 * env, 0.01, o, kNoPhi) == kRouteRPNearBeam);
  CHECK(route_charged(0.50, 1e-3, 0.05, o, kNoPhi) == kRouteOMD);
  CHECK(route_charged(3.0 / 7.0, 1e-3, 0.05, o, kNoPhi) == kRouteLost);
  // over-rigid: the 7Li triton at R = 1.29 is on the pots' inner half,
  // the 6Li 3He+t triton at R = 1.5044 is past the last module
  CHECK(route_charged(9.0 / 7.0, 1e-3, 0.05, o, kNoPhi) == kRouteRPInner);
  CHECK(route_charged(1.5044, 1e-3, 0.05, o, kNoPhi) == kRouteLost);
  CHECK(route_charged(0.857, 8e-3, 0.05, o, kNoPhi) == kRouteB0);
  CHECK(route_neutral(1e-3) == kRouteZDC);
  CHECK(route_neutral(1e-2) == kRouteLost);
  CHECK(rp_accepted(kRouteRomanPots));
  CHECK(rp_accepted(kRouteRPNearBeam));
  CHECK(!rp_accepted(kRouteB0));
  CHECK(!rp_accepted(kRouteRPInner));
}

TEST_CASE("farforward: the near-beam cut is angular, not a momentum threshold") {
  const Optics& ha = HIGH_ACCEPTANCE();
  CHECK_CLOSE(ha.pt_cut_near_beam(), 0.20, 1e-12);
  CHECK_CLOSE(HIGH_DIVERGENCE().pt_cut_near_beam(), 0.45, 1e-12);
  const double p_alpha = 4 * 137.5;
  CHECK_CLOSE(ha.pt_cut_for(p_alpha), 0.40, 1e-6);
  CHECK_CLOSE(0.20 / (ha.sigma_theta * p_alpha), 5.0, 1e-6);
  CHECK_CLOSE(ha.pt_cut_for(6 * 20.5), 0.0895, 2e-3);
  CHECK_CLOSE(ha.pt_cut_for(6 * 50.0), 0.218, 2e-3);
  CHECK_CLOSE(ha.pt_cut_for(6 * 137.5), 0.600, 2e-3);
}

TEST_CASE("farforward: the envelope is a RECTANGLE once the azimuth is known") {
  const std::vector<BeamConfig> cfgs = default_configs("6Li");
  const Optics o = yr_optics("6Li", cfgs[0].ion_momentum_per_nucleon, true);
  CHECK_CLOSE_AT(o.envelope_x(), 2.20e-3, 0.0, 1e-5);
  CHECK_CLOSE_AT(o.envelope_y(), 3.80e-3, 0.0, 1e-5);
  CHECK(route_charged(1.0, 3.0e-3, 0.1, o, 0.0) == kRouteRPNearBeam);
  CHECK(route_charged(1.0, 3.0e-3, 0.1, o, 0.5 * 3.14159265358979323846)
        == kRouteLost);
  // dropping the azimuth degenerates to the inscribed circle and accepts the
  // vertical fragment too -- the 1.7x overstatement plans/09 B2 measures
  CHECK(route_charged(1.0, 3.0e-3, 0.1, o, std::nan("")) == kRouteRPNearBeam);
  // the configuration key and the species step
  CHECK(yr_config_key("6Li", cfgs[0].ion_momentum_per_nucleon) == "5x41");
  CHECK(yr_config_key("6Li", cfgs[1].ion_momentum_per_nucleon) == "10x100");
  CHECK(yr_config_key("6Li", cfgs[2].ion_momentum_per_nucleon) == "18x275");
  double h = 0, v = 0;
  sigma_theta_for("6Li", cfgs[2].ion_momentum_per_nucleon, true, h, v);
  // 6Li is rigidity-capped at the top configuration, so it picks up
  // sqrt(beta*gamma_p / beta*gamma_ion) on the 65 urad proton divergence --
  // sqrt(2) up to the 0.1 GeV/u rounding of the capped momentum (2.5e-3).
  CHECK_CLOSE(h / 65e-6, std::sqrt(2.0), 3e-3);
  CHECK_CLOSE(v, h, 1e-12);
}

// P5: the pot transport levers are `polli_fastsim.farforward.POT_LEVERS`
// (R12, R34, D) verbatim.  R34 at 5x41 was carried as -1 ("never measured")
// until 2026-08-29 -- it is 4.56 m, and this table is the only C++ copy.
TEST_CASE("spectator: pot levers pin farforward.POT_LEVERS") {
  struct Row { const char* key; double r12, r34, d; };
  const Row rows[3] = {{"5x41", 19.24, 4.56, 0.311},
                       {"10x100", 21.25, 3.35, 0.287},
                       {"18x275", 29.97, 2.93, 0.292}};
  for (const Row& r : rows) {
    const PotLevers& l = pot_levers(r.key);
    CHECK(l.r12 == r.r12);
    CHECK(l.r34 == r.r34);
    CHECK(l.dispersion == r.d);
    // every configuration's vertical lever is now MEASURED, so none of them
    // may carry the "never measured" sentinel
    CHECK(l.r34 > 0.0);
  }
  CHECK_THROWS(pot_levers("5x100"));
}
