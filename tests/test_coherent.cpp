// SPDX-License-Identifier: GPL-3.0-or-later
// Coherent (intact-6Li) channel: the scenario model, the analytic tagging
// identities, the m-state relation of the deformation mechanism, the recoil
// kinematics, and the breakup / veto table routing -- against
// validation/reference/coherent.json and evgen/tests/test_coherent.py.

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "check_close.hpp"
#include "doctest.h"
#include "json_min.hpp"
#include "lipolgen/beams.hpp"
#include "lipolgen/bookkeeping.hpp"
#include "lipolgen/cluster_config.hpp"
#include "lipolgen/coherent.hpp"
#include "lipolgen/constants.hpp"
#include "lipolgen/event.hpp"
#include "lipolgen/rc.hpp"
#include "lipolgen/rng.hpp"
#include "lipolgen/sampler.hpp"
#include "lipolgen/spectator.hpp"

using namespace lipolgen;

#ifndef LIPOLGEN_REFERENCE_DIR
#define LIPOLGEN_REFERENCE_DIR "validation/reference"
#endif

namespace {

constexpr double kRtol = 1e-12;

bool load_coherent(jsonmin::Value& out) {
  return jsonmin::load_file(std::string(LIPOLGEN_REFERENCE_DIR) + "/coherent.json",
                            out);
}

}  // namespace

TEST_CASE("coherent: constants and scenario defaults against polligen") {
  jsonmin::Value ref;
  if (!load_coherent(ref)) {
    MESSAGE("coherent.json not found -- skipping");
    return;
  }
  CHECK_CLOSE(GEV_PER_FM_INV, ref["constants"]["GEV_PER_FM_INV"].num(), kRtol);
  CHECK_CLOSE(M_LI6_DOC, ref["constants"]["M_LI6"].num(), kRtol);
  CHECK_CLOSE(RATE_WEIGHT_SYST, ref["constants"]["RATE_WEIGHT_SYST"].num(), kRtol);
  const CoherentScenario sc;
  CHECK_CLOSE(sc.f0, ref["scenario_defaults"]["f0"].num(), kRtol);
  CHECK_CLOSE(sc.x_coh, ref["scenario_defaults"]["x_coh"].num(), kRtol);
  CHECK_CLOSE(sc.slope_b, ref["scenario_defaults"]["slope_b"].num(), kRtol);
  CHECK_CLOSE(sc.amp, ref["scenario_defaults"]["amp"].num(), kRtol);
  CHECK_CLOSE(sc.eps_b0, ref["scenario_defaults"]["eps_b0"].num(), kRtol);
  // the digitized Mantysaari table
  const jsonmin::Object& tab = ref["constants"]["MANTYSAARI_A2_DEUTERON"].obj();
  CHECK(tab.size() == mantysaari_a2_deuteron().size());
  for (const MantysaariRow& r : mantysaari_a2_deuteron()) {
    bool found = false;
    for (const auto& kv : tab) {
      if (std::fabs(std::stod(kv.first) - r.t_abs) < 1e-12) {
        CHECK_CLOSE(r.a2_m0, kv.second[0].num(), kRtol);
        CHECK_CLOSE(r.a2_m1, kv.second[1].num(), kRtol);
        found = true;
      }
    }
    CHECK(found);
  }
}

TEST_CASE("coherent: the scenario formulas against polligen") {
  jsonmin::Value ref;
  if (!load_coherent(ref)) {
    MESSAGE("coherent.json not found -- skipping");
    return;
  }
  const CoherentScenario sc;
  for (const jsonmin::Value& r : ref["gaussian_slope"].arr()) {
    CHECK_CLOSE(gaussian_slope(r["r_rms_fm"].num()), r["B"].num(), kRtol);
  }
  for (const jsonmin::Value& r : ref["coherent_fraction"].arr()) {
    CHECK_CLOSE(sc.coherent_fraction(r["x"].num()), r["value"].num(), kRtol);
  }
  for (const jsonmin::Value& r : ref["tag_acceptance"].arr()) {
    CHECK_CLOSE(sc.tag_acceptance(r["pt_cut"].num()), r["value"].num(), kRtol);
    CHECK_CLOSE(sc.mean_t_tagged(r["pt_cut"].num()), r["mean_t_tagged"].num(),
                kRtol);
  }
  for (const jsonmin::Value& r : ref["tag_acceptance_angular"].arr()) {
    CHECK_CLOSE(sc.tag_acceptance_angular(r["sigma_theta"].num(),
                                          r["p_per_nucleon"].num(),
                                          static_cast<int>(r["a_beam"].num()),
                                          r["n_sigma"].num()),
                r["value"].num(), kRtol);
  }
  for (const jsonmin::Value& r : ref["a2_deformation"].arr()) {
    CHECK_CLOSE_AT(sc.a2_deformation(r["t_abs"].num(), r["pzz"].num()),
                   r["a2_deformation"].num(), kRtol, 1e-300);
    CHECK_CLOSE_AT(sc.cos2phi_coefficient_deformation(r["t_abs"].num(),
                                                      r["pzz"].num()),
                   r["cos2phi_coefficient_deformation"].num(), kRtol, 1e-300);
  }
  for (const jsonmin::Value& r : ref["a2_tagged"].arr()) {
    CHECK_CLOSE(sc.a2_tagged(r["pt_cut"].num(), r["pzz"].num()),
                r["value"].num(), kRtol);
  }
  for (const jsonmin::Value& r : ref["recoil_lab"].arr()) {
    const CoherentRecoil rec = recoil_lab(r["t_abs"].num(), r["phi_t"].num(),
                                          r["p_per_nucleon"].num());
    CHECK_CLOSE_AT(rec.pT, r["pT"].num(), kRtol, 1e-300);
    CHECK_CLOSE_AT(rec.theta, r["theta"].num(), kRtol, 1e-300);
    CHECK_CLOSE(rec.R, r["R"].num(), kRtol);
    CHECK_CLOSE(rec.xL, r["xL"].num(), kRtol);
  }
}

TEST_CASE("coherent: the coherent fraction and the slope have the right shape") {
  const CoherentScenario sc;
  CHECK_CLOSE(sc.coherent_fraction(1e-4), sc.f0, 1e-3);
  CHECK_CLOSE(sc.coherent_fraction(sc.x_coh), sc.f0 / 2.0, kRtol);
  CHECK(sc.coherent_fraction(0.1) < sc.f0 / 50.0);
  for (double x = 1e-4; x < 1.0; x *= 1.5) {
    CHECK(sc.coherent_fraction(x * 1.5) < sc.coherent_fraction(x));
  }
  CHECK(gaussian_slope(2.4) > 45.0);
  CHECK(gaussian_slope(2.4) < 55.0);
}

// T10a (open item O5).  The eSTARlight table is the SINGLE code home of the
// unpolarized coherent rates of estarlight_li6.md; the reach arithmetic built
// on it is `validation/o5_a2_reach.py` + `python/tests/test_o5_reach.py`, and
// the write-up is phase_C_numbers.md sec. C2.  What is pinned here is the
// table and the two closed forms the verdict turns on, so that the C++ side
// fails too if an input moves.
TEST_CASE("coherent: the eSTARlight 6Li baseline, and O5's reach arithmetic") {
  const std::vector<EstarlightLi6Row>& t = estarlight_li6_coherent();
  REQUIRE(t.size() == 3);
  CHECK(std::string(t[0].vm) == "jpsi");
  CHECK(std::string(t[1].vm) == "phi");
  CHECK(std::string(t[2].vm) == "rho");
  CHECK_CLOSE(t[0].sigma_nb, 1.773, 1e-12);
  CHECK_CLOSE(t[1].sigma_nb, 30.16, 1e-12);
  CHECK_CLOSE(t[2].sigma_nb, 506.4, 1e-12);
  CHECK_CLOSE(t[0].sigma_q7_nb, 0.605, 1e-12);
  CHECK_CLOSE(t[1].sigma_q7_nb, 1.897, 1e-12);
  CHECK_CLOSE(t[2].sigma_q7_nb, 15.32, 1e-12);
  CHECK_CLOSE(t[0].sigma_rmeas_nb, 1.255, 1e-12);
  CHECK_CLOSE(t[1].sigma_rmeas_nb, 22.08, 1e-12);
  CHECK_CLOSE(t[2].sigma_rmeas_nb, 379.5, 1e-12);
  CHECK_CLOSE(t[0].branching, 0.0597, 1e-12);   // J/psi -> e+e- ONLY
  CHECK_CLOSE(t[1].branching, 0.49, 1e-12);
  CHECK_CLOSE(t[2].branching, 1.0, 1e-12);
  // 2026-09-04: the sec. 2b slope column now has a code home (it used to be a
  // dict inside validation/o5_a2_reach.py, unreachable from a re-fit) ...
  CHECK_CLOSE(t[0].b_rmeas, 55.0, 1e-12);
  CHECK_CLOSE(t[1].b_rmeas, 54.8, 1e-12);
  CHECK_CLOSE(t[2].b_rmeas, 54.1, 1e-12);
  for (const EstarlightLi6Row& r : t) CHECK(r.b_rmeas > r.b_fit);
  // ... and so does the SECOND lepton channel, which the reach arithmetic
  // used to mention in prose and never use.  J/psi -> e+e- AND mu+mu-, from
  // eSTARlight's own JpsiBree + JpsiBrmumu = 0.05971 + 0.05961.
  CHECK_CLOSE(t[0].branching_all, 0.11932, 1e-12);
  CHECK_CLOSE(t[1].branching_all, 0.49, 1e-12);   // phi has no second channel
  CHECK_CLOSE(t[2].branching_all, 1.0, 1e-12);    // nor does rho
  CHECK_CLOSE(std::sqrt(t[0].branching_all / t[0].branching), 1.41374, 1e-4);
  // The efficiency is a NAMED CONSTANT now, not a literal in a script.
  CHECK_CLOSE(COHERENT_JPSI_EFF_IR8_LI7, 0.1775, 1e-12);

  // The slope is VM-independent -- (max - min)/mean = 1.54 %, the "1.5 %" of
  // estarlight_li6.md sec. 4 -- and every row UNDERSHOOTS the analytic
  // R_G^2/(3 hbar^2 c^2) = 40.70 by 3.7-5.2 %, one-sided.
  double lo = t[0].b_fit, hi = t[0].b_fit, sum = 0.0;
  for (const EstarlightLi6Row& r : t) {
    lo = std::min(lo, r.b_fit);
    hi = std::max(hi, r.b_fit);
    sum += r.b_fit;
  }
  CHECK_CLOSE((hi - lo) / (sum / 3.0), 0.0154, 2e-2);
  const double b_analytic = gaussian_slope(1.2 * std::cbrt(6.0));
  CHECK_CLOSE(b_analytic, 40.7026, 1e-4);
  CHECK_CLOSE(lo / b_analytic - 1.0, -0.0517, 1e-2);
  CHECK_CLOSE(hi / b_analytic - 1.0, -0.0369, 1e-2);
  for (const EstarlightLi6Row& r : t) {
    CHECK(r.b_fit < b_analytic);          // one-sided
    CHECK(r.sigma_rmeas_nb < r.sigma_nb);       // bigger R -> smaller sigma
  }

  // O5's two closed forms.  <t^2> of a truncated exponential, and the reach
  // delta(kappa) = 1/sqrt(2 <P_zz^2> N <t^2>) that follows from the same
  // sqrt(2/N) `estimators::cos2phi_fit_err` already carries.
  const CoherentScenario sc;                       // slope_b = 50
  const double tm = COHERENT_T_MAX_DEFAULT;        // 0.2
  const double u = sc.slope_b * tm;
  const double num = 1.0 - std::exp(-u) * (1.0 + u + 0.5 * u * u);
  const double t2 = (2.0 / (sc.slope_b * sc.slope_b)) * num
                    / (1.0 - std::exp(-u));
  CHECK_CLOSE(t2, 7.978207e-4, 1e-5);
  CHECK_CLOSE(std::sqrt(t2), 0.028246, 1e-4);
  // |t| < 0.2 keeps 99.72 % of the Fisher information and 99.995 % of the rate
  CHECK_CLOSE(t2 * (1.0 - std::exp(-u)) / (2.0 / (sc.slope_b * sc.slope_b)),
              0.99723, 1e-4);
  CHECK(1.0 - std::exp(-u) > 0.9999);

  // The tensor signal at the MEASURED quadrupole, and the fact the |t| slope
  // knocks it down by 10.6x from the a_2(0.3) that gets quoted.
  const double kappa = a2_from_quadrupole(2.0 * LI6_QUADRUPOLE_FM2, 6, 1.0, 1);
  CHECK_CLOSE(kappa, 0.08752978, 1e-6);
  CHECK(kappa > 0.0);                              // opposite to the deuteron
  CHECK_CLOSE(0.3 * kappa / (kappa * std::sqrt(t2)), 10.62, 1e-3);

  // tensor_flip_plan(0.6) -> P_zz = {+0.6, -1.2}, so <P_zz^2> = 2.5 pzz^2.
  const RunPlan plan = tensor_flip_plan(Scenario().pol_ion_tensor);
  REQUIRE(plan.categories().size() == 2);
  double pzz2 = 0.0;
  for (const SpinCategory& c : plan.categories()) {
    const double p = 1.0 - 3.0 * c.populations[1];
    pzz2 += c.lumi_fraction * p * p;
  }
  CHECK_CLOSE(pzz2, 0.90, 1e-12);
  CHECK_CLOSE(pzz2, 2.5 * 0.6 * 0.6, 1e-12);

  // The BACKGROUND-IMMUNE two-fill difference -- the estimator sec. C2.5's
  // separation argument actually requires -- costs 1/1.1111 of that:
  // <P_zz^2>_eff = 0.81.  The headline is quoted on THIS one.
  const double p1 = 1.0 - 3.0 * plan.categories()[0].populations[1];
  const double p2 = 1.0 - 3.0 * plan.categories()[1].populations[1];
  const double v_diff = (2.0 / 0.5 + 2.0 / 0.5) / ((2.0 * (p1 - p2))
                                                   * (2.0 * (p1 - p2)));
  const double v_opt = 1.0 / (2.0 * (0.5 * p1 * p1 + 0.5 * p2 * p2));
  CHECK_CLOSE(std::sqrt(v_diff / v_opt), 1.0541, 1e-3);
  const double pzz2_diff = pzz2 / (v_diff / v_opt);
  CHECK_CLOSE(pzz2_diff, 0.81, 1e-12);

  // THE RESTRICTED ROW, kept because it is what the item shipped with and it
  // is still arithmetically right FOR THAT ROW.  Coherent J/psi at one EIC
  // year (10 fb^-1/u per nucleon, hence 10/6 fb^-1 of e+6Li), J/psi -> e+e-
  // ONLY, over 0.1 < Q^2 < 100 GeV^2 ONLY, on the OPTIMAL <P_zz^2> = 0.90:
  // 3.13e4 events, 0.749 sigma, 3 sigma at 160 fb^-1/u.
  const double lumi_nucleus_fb = Scenario().lumi_fb_per_nucleon / 6.0;
  const double n_det = t[0].sigma_nb * 1e6 * lumi_nucleus_fb * t[0].branching
                       * (1.0 - std::exp(-t[0].b_fit * tm))
                       * COHERENT_JPSI_EFF_IR8_LI7;
  CHECK_CLOSE(n_det, 3.130e4, 2e-3);
  const double u0 = t[0].b_fit * tm;
  const double t2j = (2.0 / (t[0].b_fit * t[0].b_fit))
                     * (1.0 - std::exp(-u0) * (1.0 + u0 + 0.5 * u0 * u0))
                     / (1.0 - std::exp(-u0));
  const double dk = 1.0 / std::sqrt(2.0 * pzz2 * n_det * t2j);
  CHECK_CLOSE(0.3 * dk, 0.03505, 1e-3);
  const double s_jpsi = kappa / dk;
  CHECK_CLOSE(s_jpsi, 0.749, 2e-3);
  CHECK(s_jpsi < 1.0);
  CHECK_CLOSE(Scenario().lumi_fb_per_nucleon * (3.0 / s_jpsi) * (3.0 / s_jpsi),
              160.3, 1e-2);

  // THE VERDICT, which is NOT that row.  The same chain over the WHOLE Q^2
  // range (estarlight_li6_q2_floors(), no floor) with BOTH lepton channels
  // and the background-immune <P_zz^2>: MARGINAL, and 3 sigma INSIDE the
  // {1, 10, 100} fb^-1/u band.  This is T10c's subject; pinned here too so
  // the C++ side fails if the photoproduction rows move.
  const std::vector<EstarlightLi6Q2Row>& q = estarlight_li6_q2_floors();
  const EstarlightLi6Q2Row* full = nullptr;
  for (const EstarlightLi6Q2Row& r : q) {
    if (std::string(r.vm) == "jpsi" && r.q2_floor_gev2 == 0.0) full = &r;
  }
  REQUIRE(full != nullptr);
  const double n_full = full->sigma_nb * 1e6 * lumi_nucleus_fb
                        * t[0].branching_all
                        * (1.0 - std::exp(-full->b_fit * tm))
                        * COHERENT_JPSI_EFF_IR8_LI7;
  CHECK_CLOSE(n_full, 4.224e5, 2e-3);
  const double uf = full->b_fit * tm;
  const double t2f = (2.0 / (full->b_fit * full->b_fit))
                     * (1.0 - std::exp(-uf) * (1.0 + uf + 0.5 * uf * uf))
                     / (1.0 - std::exp(-uf));
  const double dkf = 1.0 / std::sqrt(2.0 * pzz2_diff * n_full * t2f);
  const double s_full = kappa / dkf;
  CHECK_CLOSE(s_full, 2.618, 3e-3);
  CHECK(s_full > 2.0);
  CHECK(s_full < 3.0);                              // MARGINAL, not YES
  const double lumi3 = Scenario().lumi_fb_per_nucleon
                       * (3.0 / s_full) * (3.0 / s_full);
  CHECK_CLOSE(lumi3, 13.13, 1e-2);
  CHECK(lumi3 < 100.0);   // INSIDE the band: the shipped "NO" was refuted
  CHECK(lumi3 > 1.0);
}

// T10c (open item O5, 2026-09-04).  The PHOTOPRODUCTION scan.  The 0.1 GeV^2
// floor `estarlight_li6_coherent()` carries is arXiv:2511.05638's
// acceptance-study kinematic range, not a physics window, and most of the
// coherent rate is below it.  These are the same eSTARlight runs -- same
// commit, beams and seed -- with MIN_GAMMA_Q2 lowered; the Q^2 > 0.1 rows
// reproduce sec. 2a / 2b exactly, which is what makes the rest comparable.
TEST_CASE("coherent: the eSTARlight 6Li Q^2-floor scan (photoproduction)") {
  const std::vector<EstarlightLi6Q2Row>& q = estarlight_li6_q2_floors();
  const std::vector<EstarlightLi6Row>& t = estarlight_li6_coherent();
  REQUIRE(q.size() == 9);                    // 3 vector mesons x 3 floors

  // The Q^2 > 0.1 rows ARE sec. 2a / 2b: same sigma, same slope, bit for bit.
  for (const EstarlightLi6Row& base : t) {
    const EstarlightLi6Q2Row* at01 = nullptr;
    for (const EstarlightLi6Q2Row& r : q) {
      if (std::string(r.vm) == base.vm && r.q2_floor_gev2 == 0.1) at01 = &r;
    }
    REQUIRE(at01 != nullptr);
    CHECK_CLOSE(at01->sigma_nb, base.sigma_nb, 1e-12);
    CHECK_CLOSE(at01->b_fit, base.b_fit, 1e-12);
    CHECK_CLOSE(at01->sigma_rmeas_nb, base.sigma_rmeas_nb, 1e-12);
    CHECK_CLOSE(at01->b_rmeas, base.b_rmeas, 1e-12);
  }

  // Removing the floor MULTIPLIES the rate -- 6.75x for J/psi, and more for
  // the lighter mesons, whose Q^2 suppression sets in earlier.
  const double want_full[3] = {6.7518, 21.696, 35.195};
  const double want_001[3] = {1.8883, 3.4213, 4.5119};
  const char* vms[3] = {"jpsi", "phi", "rho"};
  for (int i = 0; i < 3; ++i) {
    const EstarlightLi6Q2Row *a = nullptr, *b = nullptr, *c = nullptr;
    for (const EstarlightLi6Q2Row& r : q) {
      if (std::string(r.vm) != vms[i]) continue;
      if (r.q2_floor_gev2 == 0.1) a = &r;
      if (r.q2_floor_gev2 == 0.01) b = &r;
      if (r.q2_floor_gev2 == 0.0) c = &r;
    }
    REQUIRE(a != nullptr);
    REQUIRE(b != nullptr);
    REQUIRE(c != nullptr);
    CHECK_CLOSE(c->sigma_nb / a->sigma_nb, want_full[i], 1e-3);
    CHECK_CLOSE(b->sigma_nb / a->sigma_nb, want_001[i], 1e-3);
    // monotone in the floor, on BOTH densities
    CHECK(c->sigma_nb > b->sigma_nb);
    CHECK(b->sigma_nb > a->sigma_nb);
    CHECK(c->sigma_rmeas_nb > b->sigma_rmeas_nb);
    CHECK(b->sigma_rmeas_nb > a->sigma_rmeas_nb);
    // the measured-radius density is always the smaller rate, always the
    // steeper slope -- sigma and B are not independent knobs
    CHECK(c->sigma_rmeas_nb < c->sigma_nb);
    CHECK(c->b_rmeas > c->b_fit);
    // and the |t| SLOPE barely moves with the floor: the recoil p_T spectrum
    // the far-forward acceptance cuts on is the same sample.  < 1.5 %.
    CHECK(std::abs(c->b_fit / a->b_fit - 1.0) < 0.015);
    CHECK(std::abs(c->b_rmeas / a->b_rmeas - 1.0) < 0.015);
    // <W> FALLS as the floor is removed, on every meson.  In code since
    // 2026-09-04: both efficiency-transfer arguments -- the Q^2-floor leg and
    // the beam-energy leg of T10d -- run on it, and it lived only in prose.
    CHECK(c->w_mean_gev < a->w_mean_gev);
    CHECK(a->w_mean_gev > 0.0);
  }
  const EstarlightLi6Q2Row* jf = nullptr;
  for (const EstarlightLi6Q2Row& r : q) {
    if (std::string(r.vm) == "jpsi" && r.q2_floor_gev2 == 0.0) jf = &r;
  }
  REQUIRE(jf != nullptr);
  CHECK_CLOSE(jf->w_mean_gev, 30.20, 1e-12);   // estarlight_li6.md sec. 2f
}

// T10d (open item O5, 2026-09-04, third pass).  THE EFFICIENCY CHAIN.
// COHERENT_JPSI_EFF_IR8_LI7 = 0.1775 is a 7Li number at 7Li's own TOP energy
// and it is applied to a 6Li sample at 10 x 99.5 GeV/u.  Until this pass the
// tree recorded that as "none at 10 x 99.5" -- a gap with no direction and no
// size -- while the SAME paper measures the beam-energy dependence on 3He.
// These two tables are that measurement and the species list beside it, and
// what is pinned here is the DIRECTION of each lever and its SIZE.
TEST_CASE("coherent: the far-forward efficiency is not flat in beam energy") {
  const std::vector<Chang26EffEnergyRow>& e = chang26_he3_energy_scan();
  REQUIRE(e.size() == 3);
  // arXiv:2511.05638 sec. V.B, verbatim.
  CHECK_CLOSE(e[0].e_electron_gev, 18.0, 1e-12);
  CHECK_CLOSE(e[0].e_ion_gev, 183.0, 1e-12);
  CHECK_CLOSE(e[0].efficiency, 0.3223, 1e-12);
  CHECK_CLOSE(e[1].e_ion_gev, 100.0, 1e-12);
  CHECK_CLOSE(e[1].efficiency, 0.5438, 1e-12);
  CHECK_CLOSE(e[2].e_ion_gev, 41.0, 1e-12);
  CHECK_CLOSE(e[2].efficiency, 0.9977, 1e-12);
  // MONOTONE: lower energy, higher efficiency.  This is the direction the
  // whole correction rests on, and it is a measurement.
  for (std::size_t i = 1; i < e.size(); ++i) {
    CHECK(e[i].e_ion_gev < e[i - 1].e_ion_gev);
    CHECK(e[i].efficiency > e[i - 1].efficiency);
    CHECK(e[i].efficiency <= 1.0);
  }
  // ... and the two local power laws, which are what o5_a2_reach.py's
  // `beam_energy_scaling` uses.  They differ because eff is capped at 1.
  const double s1 = std::log(e[1].efficiency / e[0].efficiency)
                    / std::log(e[1].e_ion_gev / e[0].e_ion_gev);
  const double s2 = std::log(e[2].efficiency / e[1].efficiency)
                    / std::log(e[2].e_ion_gev / e[1].e_ion_gev);
  CHECK_CLOSE(s1, -0.8656, 1e-3);
  CHECK_CLOSE(s2, -0.6807, 1e-3);
  CHECK(s1 < 0.0);
  CHECK(s2 < 0.0);
  CHECK(s2 > s1);                       // flattening, as the cap requires

  // THE STEP THE O5 CHAIN NEEDS: 117.9 -> 99.5 GeV/u, and the factor is
  // 1.12-1.16, NOT the paper's own 183 -> 100 ratio of 1.687.  That ratio is
  // a step 3.56x larger in ln E; carrying it whole is the wrong lever arm.
  CHECK_CLOSE(COHERENT_JPSI_EFF_IR8_LI7_E_ION_GEV, 117.9, 1e-12);
  CHECK_CLOSE(COHERENT_JPSI_EFF_IR8_LI7_W_MEAN_GEV, 43.2, 1e-12);
  const double e_to = 99.5;
  const double f1 = std::exp(s1 * std::log(e_to
                                           / COHERENT_JPSI_EFF_IR8_LI7_E_ION_GEV));
  const double f2 = std::exp(s2 * std::log(e_to
                                           / COHERENT_JPSI_EFF_IR8_LI7_E_ION_GEV));
  CHECK_CLOSE(f1, 1.15821, 1e-4);
  CHECK_CLOSE(f2, 1.12243, 1e-4);
  CHECK(f1 > 1.0);                      // UP: the transfer is conservative
  CHECK(f2 > 1.0);
  CHECK(f1 < 1.4);                      // and it is NOT 1.687
  CHECK_CLOSE(e[1].efficiency / e[0].efficiency, 1.6872, 1e-3);
  CHECK(std::log(e[1].e_ion_gev / e[0].e_ion_gev)
            / std::log(e_to / COHERENT_JPSI_EFF_IR8_LI7_E_ION_GEV) > 3.0);
  // The sample this tree prices sits BELOW the <W> the efficiency was
  // measured at, which by the same paper's Fig. 2 points the same way.
  const EstarlightLi6Q2Row* jf = nullptr;
  for (const EstarlightLi6Q2Row& r : estarlight_li6_q2_floors()) {
    if (std::string(r.vm) == "jpsi" && r.q2_floor_gev2 == 0.0) jf = &r;
  }
  REQUIRE(jf != nullptr);
  CHECK(jf->w_mean_gev < COHERENT_JPSI_EFF_IR8_LI7_W_MEAN_GEV);
}

TEST_CASE("coherent: the species list is fixed RIGIDITY, NOT fixed E/u") {
  const std::vector<Chang26SpeciesEffRow>& t = chang26_species_efficiency();
  REQUIRE(t.size() == 7);
  // Every entry is at its own top energy Z/A x 275, so A/Z x E is the SAME
  // 275 GeV/e rigidity for all seven.  The fall from 2D to 16O is therefore
  // NOT a rigidity effect.  THAT IS THE WHOLE OF WHAT THE RIGIDITY BUYS.
  for (const Chang26SpeciesEffRow& r : t) {
    const double rigidity = static_cast<double>(r.a) / r.z * r.e_ion_gev;
    CHECK(rigidity > 273.0);
    CHECK(rigidity < 276.0);
    CHECK(r.efficiency > 0.0);
    CHECK(r.efficiency < 1.0);
  }
  // ... and it falls monotonically with A.
  for (std::size_t i = 1; i < t.size(); ++i) {
    CHECK(t[i].a > t[i - 1].a);
    CHECK(t[i].efficiency < t[i - 1].efficiency);
  }
  // The 7Li row IS COHERENT_JPSI_EFF_IR8_LI7 -- one number, one home.
  CHECK_CLOSE(t[3].efficiency, COHERENT_JPSI_EFF_IR8_LI7, 1e-12);
  CHECK(std::string(t[3].nucleus) == "7Li");
  CHECK_CLOSE(t[0].efficiency, 0.4712, 1e-12);    // 2D
  CHECK_CLOSE(t[1].efficiency, 0.3223, 1e-12);    // 3He, = the scan's top row
  CHECK_CLOSE(t[1].efficiency, chang26_he3_energy_scan()[0].efficiency, 1e-12);
  CHECK_CLOSE(t[2].efficiency, 0.2942, 1e-12);    // 4He
  CHECK_CLOSE(t[6].efficiency, 0.0159, 1e-12);    // 16O

  // THE RETRACTED SENTENCE.  Until 2026-09-04 this test was named "so it is a
  // lever" and five sites said "every entry is at the same rigidity, so its
  // A-ordering is a species lever on its own".  Eliminating rigidity does NOT
  // leave species alone: at fixed R = A E/(Z) the per-nucleon energy is
  // E/u = R Z/A and the total momentum is p_z = Z R, and BOTH vary down the
  // list.  This block is the retraction, as arithmetic.
  double eu_lo = 1e9, eu_hi = 0.0, pz_lo = 1e9, pz_hi = 0.0;
  int z_lo = 99, z_hi = 0;
  for (const Chang26SpeciesEffRow& r : t) {
    eu_lo = std::min(eu_lo, r.e_ion_gev);
    eu_hi = std::max(eu_hi, r.e_ion_gev);
    pz_lo = std::min(pz_lo, r.z * r.e_ion_gev * r.a / static_cast<double>(r.z));
    pz_hi = std::max(pz_hi, r.z * r.e_ion_gev * r.a / static_cast<double>(r.z));
    z_lo = std::min(z_lo, r.z);
    z_hi = std::max(z_hi, r.z);
  }
  CHECK_CLOSE(eu_lo, 118.0, 1e-12);   // 7Li -- the anchor the chain uses
  CHECK_CLOSE(eu_hi, 183.0, 1e-12);   // 3He
  CHECK(eu_hi / eu_lo > 1.5);         // a x1.55 spread in E/u, at fixed R
  CHECK(z_lo == 1);
  CHECK(z_hi == 8);
  CHECK_CLOSE(pz_lo, 274.0, 1e-9);    // = A x E, the total beam momentum
  CHECK_CLOSE(pz_hi, 2192.0, 1e-9);
  CHECK(pz_hi / pz_lo > 7.0);
  // AND THE CONFOUND IS THE SIZE OF THE EFFECT BEING READ OFF.  The chain's
  // OTHER leg measures d ln(eff)/d ln(E_ion) on 3He; carried across the
  // list's own E/u spread it is worth up to x1.46 -- as large as the whole
  // claimed "species gain" of x1.16-1.22.
  const std::vector<Chang26EffEnergyRow>& e = chang26_he3_energy_scan();
  const double slope_steep =
      std::log(e[1].efficiency / e[0].efficiency)
      / std::log(e[1].e_ion_gev / e[0].e_ion_gev);
  const double eu_worth = std::exp(slope_steep * std::log(eu_lo / eu_hi));
  CHECK_CLOSE(eu_worth, 1.4620, 1e-3);
  CHECK(eu_worth > 1.16);

  // WHAT THE TABLE DOES SUPPORT WITHOUT THE CONFOUND: the FOUR entries that
  // share a beam energy -- 2D, 4He, 12C, 16O, all at 137 GeV/u.  6Li's own
  // fixed-rigidity energy is Z/A x 275 = 137.5, so those four bracket A = 6
  // with NO energy step at all.
  int same_energy = 0;
  for (const Chang26SpeciesEffRow& r : t) {
    if (r.e_ion_gev == 137.0) ++same_energy;
  }
  CHECK(same_energy == 4);
  CHECK_CLOSE(0.5 * 275.0, 137.5, 1e-12);   // 6Li's own Z/A x 275

  // THE 7Li -> 6Li SUBSTITUTION, priced with this repository's OWN acceptance
  // model: eff = exp(-B pT_cut^2) = `CoherentScenario::tag_acceptance`, with
  // B = `gaussian_slope(1.2 A^(1/3))` and pT_cut inverted from the 7Li row.
  const auto r_g = [](double a) { return 1.2 * std::cbrt(a); };
  const auto b_of = [&](double a) { return gaussian_slope(r_g(a)); };
  const double pt2 = std::log(1.0 / COHERENT_JPSI_EFF_IR8_LI7) / b_of(7.0);
  CHECK_CLOSE(std::sqrt(pt2), 0.19577, 1e-4);
  const double eps6 = std::exp(-b_of(6.0) * pt2);
  CHECK_CLOSE(eps6, 0.21015, 1e-4);
  CHECK_CLOSE(eps6 / COHERENT_JPSI_EFF_IR8_LI7, 1.1839, 1e-3);
  // ... but that reading interpolates a 137 GeV/u anchor (4He) onto a
  // 118 GeV/u one (7Li), so it is the confounded one.  On the SAME-ENERGY
  // set -- 4He -> 12C, both at 137, bracketing A = 6 -- the same four closed
  // forms straddle 1 instead of clearing it.
  const double e4 = t[2].efficiency, e12 = t[5].efficiency;
  const double w_a = (6.0 - 4.0) / (12.0 - 4.0);
  const double w_r = (r_g(6.0) - r_g(4.0)) / (r_g(12.0) - r_g(4.0));
  const double se_log_a =
      std::exp(std::log(e4) + (std::log(e12) - std::log(e4)) * w_a);
  const double se_lin_a = e4 + (e12 - e4) * w_a;
  const double se_log_r =
      std::exp(std::log(e4) + (std::log(e12) - std::log(e4)) * w_r);
  const double se_pt4 =
      std::exp(-b_of(6.0) * (std::log(1.0 / e4) / b_of(4.0)));
  const double se_pt12 =
      std::exp(-b_of(6.0) * (std::log(1.0 / e12) / b_of(12.0)));
  CHECK_CLOSE(se_log_a, 0.20059, 1e-4);
  CHECK_CLOSE(se_lin_a, 0.23655, 1e-4);
  CHECK_CLOSE(se_log_r, 0.17823, 1e-4);
  CHECK_CLOSE(se_pt4, 0.20124, 1e-4);
  CHECK_CLOSE(se_pt12, 0.176289, 1e-4);
  CHECK(se_pt12 < COHERENT_JPSI_EFF_IR8_LI7);   // BELOW 1: the sign is open
  CHECK(se_lin_a > COHERENT_JPSI_EFF_IR8_LI7);
  CHECK_CLOSE(se_pt12 / COHERENT_JPSI_EFF_IR8_LI7, 0.9932, 1e-3);
  CHECK_CLOSE(se_lin_a / COHERENT_JPSI_EFF_IR8_LI7, 1.3327, 1e-3);

  // THE ONE PLACE THE pT-THRESHOLD FORM CAN BE TESTED -- 3He/4He, both Z = 2,
  // same rigidity -- AND THE Z THAT WAS WRONG IN IT UNTIL 2026-09-04.  The
  // criterion is "a safe distance FROM THE BEAM", so the cut is on the ANGLE:
  // pT_cut = theta p_z = theta Z R.  It carries the CHARGE.  7Li is Z = 3 and
  // the pair is Z = 2, so the pair's pT_cut^2 is (2/3)^2 of 7Li's.  With that
  // factor the model reproduces the measured ratio to 0.1 %; without it, it
  // was reported as "overstating by 12 %" and that was read as an instruction
  // to take the low end of the band.  BOTH numbers are kept here: the second
  // is what the wrong Z gives.
  const double meas = t[1].efficiency / t[2].efficiency;
  CHECK_CLOSE(meas, 1.0955, 1e-3);
  const double pred_z = std::exp(-(b_of(3.0) - b_of(4.0)) * pt2 * 4.0 / 9.0);
  CHECK_CLOSE(pred_z, 1.0967, 1e-3);
  CHECK_CLOSE(pred_z / meas, 1.0011, 1e-3);
  CHECK(std::abs(pred_z / meas - 1.0) < 0.01);       // 0.1 %, not 12 %
  const double pred_wrong_z = std::exp(-(b_of(3.0) - b_of(4.0)) * pt2);
  CHECK_CLOSE(pred_wrong_z, 1.2309, 1e-3);
  CHECK_CLOSE(pred_wrong_z / meas, 1.1236, 1e-3);

  // THE THIRD READING OF THE SAME TABLE, and the reason none of them settles
  // it: a Z-INDEPENDENT pT_cut fits the LIGHT end's ABSOLUTE values far
  // better than a Z-scaled one, and the one misfit there is 3He -- by 16 %,
  // which is what its E/u = 183 against the family's 137 is worth on the
  // chain's own energy scan (an 18-22 % shortfall).  IT IS A LOCAL FORM AND
  // NOT A DESCRIPTION OF THE LIST: it misses 12C by 32 % and 16O by x3.1,
  // and both are checked here so the claim cannot be overstated later.
  const auto abs_pred = [&](double a) { return std::exp(-b_of(a) * pt2); };
  CHECK_CLOSE(abs_pred(2.0) / t[0].efficiency, 1.0025, 1e-3);   // 2D
  CHECK_CLOSE(abs_pred(4.0) / t[2].efficiency, 1.0336, 1e-3);   // 4He
  CHECK_CLOSE(abs_pred(9.0) / t[4].efficiency, 1.0469, 1e-3);   // 9Be
  CHECK_CLOSE(abs_pred(3.0) / t[1].efficiency, 1.1613, 1e-3);   // 3He: 16 %
  CHECK_CLOSE(abs_pred(12.0) / t[5].efficiency, 1.3217, 1e-3);  // 12C: 32 %
  CHECK_CLOSE(abs_pred(16.0) / t[6].efficiency, 3.1320, 1e-3);  // 16O: x3.1
  CHECK(abs_pred(12.0) / t[5].efficiency > 1.20);   // the form breaks down
  CHECK(abs_pred(16.0) / t[6].efficiency > 2.0);    // ... badly, by 16O
  const double he3_energy_penalty =
      std::exp(slope_steep * std::log(183.0 / 137.0));
  CHECK_CLOSE(he3_energy_penalty, 0.7783, 1e-3);                // 22 % down
  CHECK(1.0 / he3_energy_penalty > 1.16);   // covers 3He's own 16 % residual
}

TEST_CASE("coherent: <|t|> = 1/B and the acceptance is exp(-B c^2)") {
  const CoherentScenario sc;   // slope_b = 50
  Rng rng(2026, 0, 0, 1);
  const std::size_t n = 400000;
  double sum = 0.0;
  std::size_t above = 0;
  const double cut = HIGH_ACCEPTANCE().pt_cut_near_beam();  // 0.20 GeV
  for (std::size_t i = 0; i < n; ++i) {
    const double t = sc.sample_t(rng);
    sum += t;
    if (t > cut * cut) ++above;
  }
  CHECK_CLOSE(sum / n, 1.0 / sc.slope_b, 0.02);
  CHECK_CLOSE(static_cast<double>(above) / n, sc.tag_acceptance(cut), 0.02);
}

TEST_CASE("coherent: the deformation a_2 scales, and a_0 = -2 a_1") {
  const CoherentScenario sc;
  const double a = sc.a2_deformation(0.1, 0.6);
  CHECK_CLOSE(a, 0.6 / 4.0 * 0.08 * 50.0 * 0.1, kRtol);
  CHECK(a > 0.0);   // eps_b0 < 0 for 6Li with pzz > 0 gives a positive a_2
  CHECK_CLOSE(sc.a2_deformation(0.2, 0.6), 2.0 * a, kRtol);
  CHECK(sc.a2_deformation(0.1, 0.0) == 0.0);
  const double cut = HIGH_ACCEPTANCE().pt_cut_near_beam();
  CHECK_CLOSE(sc.a2_tagged(cut, 0.6),
              sc.a2_deformation(sc.mean_t_tagged(cut), 0.6), kRtol);
  // THE m-STATE RELATION: a_2(m=0) = -2 a_2(m=+-1), the wave-function
  // symmetry the 6Li scaling inherits from the deuteron anchor
  for (double t : {0.02, 0.05, 0.1, 0.2}) {
    CHECK_CLOSE(sc.a2_m_state(t, 0), -2.0 * sc.a2_m_state(t, +1), kRtol);
    CHECK_CLOSE(sc.a2_m_state(t, +1), sc.a2_m_state(t, -1), kRtol);
    // and the population average over p_m reproduces a2_deformation
    for (double pzz : {0.6, -0.6, 1.0}) {
      // any fill with this P_zz will do: p+ = p- = (P_zz + 2)/6,
      // p0 = (1 - P_zz)/3, whose tensor moment p+ + p- - 2 p0 is P_zz
      const double pp = (pzz + 2.0) / 6.0, p0 = (1.0 - pzz) / 3.0;
      CHECK_CLOSE_AT(2.0 * pp * sc.a2_m_state(t, 1) + p0 * sc.a2_m_state(t, 0),
                     sc.a2_deformation(t, pzz), kRtol, 1e-300);
    }
  }
  // the digitized anchor obeys the same relation in the linear regime
  for (const MantysaariRow& r : mantysaari_a2_deuteron()) {
    if (r.t_abs <= 0.2) CHECK_CLOSE(r.a2_m0, -2.0 * r.a2_m1, 0.25);
    CHECK(r.a2_m0 * r.a2_m1 < 0.0);
  }
}

TEST_CASE("coherent: the total cos 2phi coefficient carries both mechanisms") {
  // 1 + c2 cos 2(phi - phi_S) with c2 = 2 a2_deformation + amp * P_zz --
  // `money_cos2phi_coherent.py` injects the flat piece as amp * pzz and the
  // reconstructed pseudo-experiments inject the deformation as 2 a_2.
  const CoherentScenario sc;
  const double pzz = 0.6, t = 0.06;
  CHECK_CLOSE(sc.cos2phi_coefficient(t, pzz),
              2.0 * sc.a2_deformation(t, pzz) + sc.amp * pzz, kRtol);
  const CoherentScenario flat{0.04, 0.01, 50.0, 0.01, 0.0};
  CHECK_CLOSE(flat.cos2phi_coefficient(t, pzz), 0.01 * pzz, kRtol);
  const CoherentScenario deform{0.04, 0.01, 50.0, 0.0, -0.08};
  CHECK_CLOSE(deform.cos2phi_coefficient(t, pzz),
              deform.cos2phi_coefficient_deformation(t, pzz), kRtol);
}

TEST_CASE("coherent: the recoil stays inside the near-beam band") {
  for (double t : {0.01, 0.05, 0.2}) {
    const CoherentRecoil r = recoil_lab(t, 0.0, 100.0, 0.01);
    CHECK(std::fabs(r.R - 1.0) < 0.05);
    CHECK(r.theta < 5e-3);
  }
}

TEST_CASE("coherent: fragment rigidities are mass-to-charge ratios") {
  jsonmin::Value ref;
  if (!load_coherent(ref)) {
    MESSAGE("coherent.json not found -- skipping");
    return;
  }
  for (const jsonmin::Value& r : ref["fragment_rigidity"].arr()) {
    const int a = static_cast<int>(r["a"].num());
    const int z = static_cast<int>(r["z"].num());
    const int ba = static_cast<int>(r["beam_a"].num());
    const int bz = static_cast<int>(r["beam_z"].num());
    const double got = fragment_rigidity(a, z, ba, bz);
    if (z == 0) {
      CHECK(std::isnan(got));
    } else {
      CHECK_CLOSE(got, r["rigidity"].num(), kRtol);
    }
  }
  // and the breakup tables themselves
  for (const auto& kv : ref["breakup_tables"].obj()) {
    const int ba = kv.first == "6Li" ? 6 : 7;
    const std::vector<BreakupChannel>& tab = breakup_table(ba, 3);
    CHECK(tab.size() == kv.second.size());
    for (std::size_t i = 0; i < tab.size(); ++i) {
      CHECK(tab[i].name == kv.second[i]["name"].str());
      CHECK_CLOSE(tab[i].threshold_mev, kv.second[i]["threshold_mev"].num(), kRtol);
      CHECK(tab[i].fragments.size() == kv.second[i]["fragments"].size());
      for (std::size_t f = 0; f < tab[i].fragments.size(); ++f) {
        CHECK(tab[i].fragments[f].name == kv.second[i]["fragments"][f]["name"].str());
        CHECK(tab[i].fragments[f].a == static_cast<int>(kv.second[i]["fragments"][f]["a"].num()));
        CHECK(tab[i].fragments[f].z == static_cast<int>(kv.second[i]["fragments"][f]["z"].num()));
      }
    }
  }
}

TEST_CASE("coherent: the veto table routing") {
  // The routing labels, and the rigidities as PHYSICAL mass-to-charge ratios:
  // the naive (A/Z)/(A/Z)_beam disagrees with the masses the spectator module
  // boosts the same fragments with.
  std::map<std::string, std::map<std::string, VetoRow>> frag;
  for (const VetoRow& row : veto_table()) frag[row.channel][row.fragment] = row;

  const VetoRow& alpha = frag["alpha+d"]["alpha"];
  CHECK_CLOSE_AT(alpha.rigidity, 0.99813, 0.0, 5e-6);
  CHECK(alpha.destination.find("beam-blind") != std::string::npos);
  const VetoRow& d = frag["alpha+d"]["d"];
  CHECK_CLOSE_AT(d.rigidity, 1.00452, 0.0, 5e-6);
  CHECK(d.destination.find("beam-blind") != std::string::npos);
  const VetoRow& he3 = frag["3He+t"]["3He"];
  CHECK_CLOSE_AT(he3.rigidity, 0.75204, 0.0, 5e-6);
  CHECK(he3.destination == "RomanPots");
  const VetoRow& t = frag["3He+t"]["t"];
  CHECK_CLOSE_AT(t.rigidity, 1.50437, 0.0, 5e-6);
  CHECK(t.destination.find("over-rigid") != std::string::npos);
  CHECK(frag["alpha+p+n"]["n"].destination == "ZDC");
  CHECK_CLOSE_AT(frag["alpha+p+n"]["p"].rigidity, 0.50251, 0.0, 5e-6);
  CHECK(frag["alpha+p+n"]["p"].destination == "OMD");

  // 7Li: the alpha the tag rides on, and the triton the far-forward scan
  // measured onto the pots' inner half
  std::map<std::string, std::map<std::string, VetoRow>> f7;
  for (const VetoRow& row : veto_table(7, 3)) f7[row.channel][row.fragment] = row;
  CHECK_CLOSE_AT(f7["alpha+t"]["alpha"].rigidity, 0.85571, 0.0, 5e-6);
  CHECK(f7["alpha+t"]["alpha"].destination == "RomanPots");
  CHECK_CLOSE_AT(f7["alpha+t"]["t"].rigidity, 1.28971, 0.0, 5e-6);
  CHECK(f7["alpha+t"]["t"].destination == "RP-inner (over-rigid)");
  // 6He is off the mass table, so its rigidity comes from the A*M_U fallback
  CHECK_CLOSE_AT(f7["6He+p"]["6He"].rigidity, 1.28308, 0.0, 5e-6);

  // the same alpha the fast simulation boosts: ONE rigidity, not two
  const MomentumSampler ms(LI6_ALPHA_TAG());
  Rng rng(3, 0, 0, 0);
  std::vector<double> r_at_rest;
  for (int i = 0; i < 20000; ++i) {
    double kx, ky, kz;
    ms.sample(rng, kx, ky, kz);
    if (std::sqrt(kx * kx + ky * ky + kz * kz) >= 0.01) continue;
    r_at_rest.push_back(
        boost_spectator_fragment(LI6_ALPHA_TAG(), 137.5, kx, ky, kz).R);
  }
  std::sort(r_at_rest.begin(), r_at_rest.end());
  CHECK(!r_at_rest.empty());
  CHECK_CLOSE_AT(r_at_rest[r_at_rest.size() / 2], alpha.rigidity, 0.0, 1e-3);
}

namespace {

/// A representative DIS side for the coherent sampler: one accepted (x, Q2)
/// of the generator window, with the exact head-on residual
/// R = k + P_ion - k' the recoil is solved against.
CoherentDis test_dis(const BeamConfig& cfg, double x, double q2) {
  CoherentDis d;
  d.x = x;
  d.q2 = q2;
  d.w2 = w2_from_xq2(x, q2);
  const double p_a = 6.0 * cfg.ion_momentum_per_nucleon;
  const double m_a = nuclear_mass(3, 6);
  const Vec4 k{cfg.electron_energy, 0.0, 0.0, -cfg.electron_energy};
  const Vec4 p_ion{std::sqrt(p_a * p_a + m_a * m_a), 0.0, 0.0, p_a};
  const double s4 = 4.0 * cfg.electron_energy * cfg.ion_momentum_per_nucleon;
  const double y = y_from_xq2(x, q2, s4);
  const ScatteredElectron e = scattered_electron(x, y, s4, cfg.electron_energy);
  const Vec4 kp{e.e_prime, e.e_prime * std::sin(e.theta), 0.0,
                e.e_prime * std::cos(e.theta)};
  d.residual = (k + p_ion) - kp;
  return d;
}

}  // namespace

TEST_CASE("coherent: the sampler's recoil, weights and event record") {
  const CoherentScenario sc;
  const BeamConfig cfg = default_configs("6Li")[1];
  const double p_u = cfg.ion_momentum_per_nucleon;
  const double pzz = 0.6, phi_s = 0.0;
  const CoherentDis dis = test_dis(cfg, 2e-3, 3.0);
  CoherentSampler s(sc, p_u, pzz, phi_s);
  Rng rng(7, 0, 0, 0);
  const std::size_t n = 200000;
  double tsum = 0.0, wsum = 0.0, wc2 = 0.0, c2sum = 0.0;
  const double m_beam = nuclear_mass(3, 6);
  for (std::size_t i = 0; i < n; ++i) {
    const CoherentEvent e = s.sample(rng, dis);
    tsum += e.t;
    wsum += e.weight;
    wc2 += e.weight * std::cos(2.0 * e.phi_t);
    c2sum += e.c2;
    if (i < 200) {
      // the recoil is on shell, at beam rigidity, with |t| = pT^2
      CHECK_CLOSE(e.p_recoil.m2(), m_beam * m_beam, 1e-9);
      CHECK_CLOSE(e.p_recoil.pt() * e.p_recoil.pt(), e.t, 1e-9);
      CHECK(std::fabs(e.recoil.R - 1.0) < NEAR_BEAM_BAND);
      CHECK_CLOSE(e.c2, sc.cos2phi_coefficient(e.t, pzz), kRtol);
    }
  }
  CHECK_CLOSE(tsum / n, 1.0 / sc.slope_b, 0.02);
  CHECK_CLOSE(c2sum / n, sc.cos2phi_coefficient(tsum / n, pzz), 0.02);
  CHECK_CLOSE(wsum / n, 1.0, 0.01);   // the weight is unbiased in the azimuth
  // E[2 w cos 2phi] = E[c2] exactly (the azimuth is uniform, so
  // E[cos^2 2phi] = 1/2); the MC error of the estimator is 2 sqrt(1/2n) =
  // 3.2e-3, so the band below is 4 sigma and not a tolerance on the physics.
  const double c2_bar = c2sum / n;
  CHECK_CLOSE(c2_bar, sc.cos2phi_coefficient(1.0 / sc.slope_b, pzz), 0.02);
  CHECK_CLOSE_AT(2.0 * wc2 / n, c2_bar, 0.0, 0.015);

  // weighted-azimuth mode: the same modulation, carried by the density
  CoherentSampler sw(sc, p_u, pzz, phi_s);
  sw.set_weighted_azimuth(true);
  double c2_hat = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const CoherentEvent e = sw.sample(rng, dis);
    CHECK(e.weight == 1.0);
    c2_hat += std::cos(2.0 * e.phi_t);
  }
  CHECK_CLOSE_AT(2.0 * c2_hat / n, c2_bar, 0.0, 0.015);

  // the event record
  Event ev;
  Particle beam_i;
  beam_i.pdg = 1000030060;
  beam_i.status = Status::Beam;
  beam_i.role = Role::BeamIon;
  const double p_beam = 6 * p_u;
  beam_i.p = Vec4{std::sqrt(p_beam * p_beam + m_beam * m_beam), 0, 0, p_beam};
  beam_i.mass = m_beam;
  beam_i.charge = 3;
  ev.particles.push_back(beam_i);
  const CoherentEvent ce = s.sample(rng, dis);
  s.fill_event(ev, ce);
  const Particle* rec = ev.find(Role::IntactRecoil);
  REQUIRE(rec != nullptr);
  CHECK(rec->pdg == 1000030060);
  CHECK(rec->charge == 3);
  CHECK(rec->status == Status::Final);
  CHECK(rec->mother1 == 0);
  CHECK_CLOSE(rec->mass, m_beam, kRtol);
  CHECK_CLOSE(ev.kin.t, ce.t, kRtol);
  CHECK(ev.kin.x_pom == ce.x_pom);
  CHECK(ev.kin.x_pom > 0.0);
  CHECK(ev.kin.m_x2 == ce.m_x2);
  CHECK(ev.kin.beta_pom == ce.beta);
  CHECK(ev.channel == Channel::CoherentLi6);
  CHECK_CLOSE(ev.weight, ce.weight, kRtol);
}

TEST_CASE("coherent: the angular near-beam cut kills the upper energies") {
  // The envelope is an ANGLE and the recoil momentum is A p_u, so the same
  // optics gives very different acceptances at the three configurations.
  const CoherentScenario sc;
  const double sig = HIGH_ACCEPTANCE().sigma_theta;
  CHECK_CLOSE_AT(sc.tag_acceptance_angular(sig, 20.5), 0.67, 0.0, 0.02);
  CHECK_CLOSE_AT(sc.tag_acceptance_angular(sig, 50.0), 0.093, 0.0, 0.005);
  CHECK(sc.tag_acceptance_angular(sig, 137.5) < 1e-7);
  const std::vector<BeamConfig> cfgs = default_configs("6Li");
  // at the machine's OWN per-configuration divergence it is dead everywhere
  double worst = 0.0;
  for (const BeamConfig& c : cfgs) {
    double h = 0.0, v = 0.0;
    sigma_theta_for("6Li", c.ion_momentum_per_nucleon, true, h, v);
    worst = std::fmax(worst,
                      sc.tag_acceptance_angular(h, c.ion_momentum_per_nucleon));
  }
  CHECK(worst < 1e-5);
}

// ---------------------------------------------------------------- open item O4
//
// T10b.  Delta B is DEFINED (coherent.hpp, `eps_b0`) and the definition is
// exercised, not merely written: Delta B_m is the cos 2(Phi - Phi_S)
// coefficient of the m state's |F|^2 slope, a_2(m) = -(Delta B_m/2)|t|, and
// eps_b0 = 2 Delta B_{+-1}/B = -Delta B_0/B.  And the consequence of the
// shipped default is pinned so it cannot be quietly forgotten:
// eps_b0 = -0.08 implies a 6Li charge quadrupole 11.4x the measured one.
// docs/open_items/run_2026-09-03/phase_C_numbers.md sec. C4.
TEST_CASE("T10b Delta B is defined, and eps_b0 = -0.08 is a deuteron number") {
  const CoherentScenario sc;

  // (1) THE DEFINITION.  a_2(m) = -(Delta B_m/2)|t|, exactly, at every m.
  for (double t : {0.01, 0.05, 0.1, 0.2, 0.3}) {
    for (int m : {0, +1, -1}) {
      CHECK_CLOSE(sc.a2_m_state(t, m), -0.5 * sc.delta_b_m(m) * t, 1e-15);
    }
  }
  // (2) THE TWO RELATIONS the old docstring got wrong by a sign.
  CHECK_CLOSE(sc.eps_b0, 2.0 * sc.delta_b_m(+1) / sc.slope_b, 1e-15);
  CHECK_CLOSE(sc.eps_b0, -sc.delta_b_m(0) / sc.slope_b, 1e-15);
  CHECK_CLOSE(sc.delta_b_m(0), -2.0 * sc.delta_b_m(+1), 1e-15);
  CHECK_CLOSE(sc.delta_b_m(+1), sc.delta_b_m(-1), 1e-15);
  // (3) `slope_at_azimuth` IS the definition: B + Delta B_m cos 2 phi, and
  // its phi average is B.  Expanding exp(-|t| B(phi)) to first order gives a
  // cos 2phi coefficient -Delta B_m |t| = 2 a_2(m), which is exactly the
  // anchor's 1 + 2 a_2 cos 2phi normalization (`cos2phi_coefficient_...`).
  const double kHalfPi = 2.0 * std::atan(1.0) * 1.0;   // pi/2
  for (int m : {0, +1}) {
    CHECK_CLOSE(sc.slope_at_azimuth(0.0, m), sc.slope_b + sc.delta_b_m(m),
                1e-15);
    CHECK_CLOSE(sc.slope_at_azimuth(kHalfPi, m), sc.slope_b - sc.delta_b_m(m),
                1e-13);
    CHECK_CLOSE_AT(sc.slope_at_azimuth(0.5 * kHalfPi, m), sc.slope_b, 1e-12,
                   1e-12);
    const double t = 0.05;
    const double num = (std::exp(-t * sc.slope_at_azimuth(0.0, m))
                        - std::exp(-t * sc.slope_at_azimuth(kHalfPi, m)))
                       / (std::exp(-t * sc.slope_at_azimuth(0.0, m))
                          + std::exp(-t * sc.slope_at_azimuth(kHalfPi, m)));
    // (F(0) - F(pi/2))/(F(0) + F(pi/2)) = -tanh(Delta B_m |t|) -> 2 a_2 at
    // small |t|; 2 % is the size of the cubic term at this |t| and m = 0.
    CHECK_CLOSE(num, 2.0 * sc.a2_m_state(t, m), 0.02);
  }
  // (4) THE SHIPPED DEFAULT IS NOT A 6Li NUMBER.  Read the quadrupole it
  // assumes straight back out of the map (`quadrupole_from_a2_slope`).
  const double q_matter = quadrupole_from_a2_slope(sc.a2_m_state(1.0, 1), 6);
  CHECK_CLOSE(q_matter, -1.8690781872, 1e-9);
  const double q_charge = 0.5 * q_matter;          // N = Z
  CHECK_CLOSE(q_charge, -0.9345390936, 1e-9);
  CHECK_CLOSE(q_charge / LI6_QUADRUPOLE_FM2, 11.4246832958, 1e-9);
  // ... and it is a bigger a_2 than the DEUTERON's own digitized one, for a
  // nucleus whose quadrupole is 3.5x smaller.
  CHECK(sc.a2_m_state(0.3, +1) > 0.29);
  for (const MantysaariRow& r : mantysaari_a2_deuteron()) {
    if (std::fabs(r.t_abs - 0.3) < 1e-9) {
      CHECK(std::fabs(r.a2_m1) < sc.a2_m_state(0.3, +1));
    }
  }
  // (5) THE MEASURED 6Li BAND, through the same map at this B.  The three
  // rows are `quadrupole_band_fm2()`; nothing is retyped.
  const ClusterConfigSampler s;
  const std::array<double, 3> band = s.quadrupole_band_fm2();
  const double eps_meas = a2_from_quadrupole(2.0 * band[0], 6, 1.0, 1) * -4.0
                          / sc.slope_b;
  const double eps_gfmc = a2_from_quadrupole(2.0 * band[1], 6, 1.0, 1) * -4.0
                          / sc.slope_b;
  const double eps_mod = a2_from_quadrupole(2.0 * band[2], 6, 1.0, 1) * -4.0
                         / sc.slope_b;
  // Measured 2026-09-05 at the defaults (B = 50): the three DERIVED eps_b0.
  // The docs' tables round these to -0.0070 / -0.0171 / -0.0527 and do their
  // arithmetic on the rounded values -- which is why "2.8000" is not a
  // publishable digit count and "2.80" is (`t_positivity_edge`).
  CHECK_CLOSE(eps_meas, -0.0070024, 1e-4);
  CHECK_CLOSE(eps_gfmc, -0.0171207, 1e-4);
  CHECK_CLOSE(eps_mod, -0.0526846, 1e-4);
  // the shipped default sits 1.52x ABOVE the top of that band ...
  CHECK(std::fabs(sc.eps_b0) > std::fabs(eps_mod));
  CHECK_CLOSE(sc.eps_b0 / eps_mod, 1.5185, 2e-3);
  // ... and the old documented band -(0.04 .. 0.13) overlaps the 6Li band
  // ONLY at its alpha+d model end.  Say it exactly: 0.0527 is inside
  // [0.04, 0.13], and the measured and GFMC rows are both BELOW the floor.
  CHECK(std::fabs(eps_mod) > 0.04);
  CHECK(std::fabs(eps_mod) < 0.13);
  CHECK(std::fabs(eps_gfmc) < 0.04);
  CHECK(std::fabs(eps_meas) < 0.04);
  // (6) eps_b0 AND slope_b ARE NOT INDEPENDENT: the physics is the product.
  for (double b : {40.0, 60.0}) {
    CoherentScenario v;
    v.slope_b = b;
    CHECK_CLOSE(0.5 * quadrupole_from_a2_slope(v.a2_m_state(1.0, 1), 6)
                    / q_charge,
                b / sc.slope_b, 1e-12);
  }
  // (7) The POSITIVITY EDGE is a consequence of (4), not a property of 6Li,
  // and since 2026-09-04 it is DERIVED by `t_positivity_edge` rather than
  // retyped as a closed form here (D5).  It is no longer the stated reason
  // for COHERENT_T_MAX_DEFAULT -- see (8) below and the constant's own
  // comment -- but it is still what `check_positivity` guards.
  CHECK_CLOSE(sc.t_positivity_edge(-2.0), 0.245, 1e-12);
  CHECK(COHERENT_T_MAX_DEFAULT < sc.t_positivity_edge(-2.0));
  CoherentScenario meas = sc;
  meas.eps_b0 = eps_meas;
  // At the MEASURED 6Li quadrupole the edge is 2.80 GeV^2, 9.3x outside the
  // |t| <= 0.30 the deformation input is digitized over: positivity stops
  // binding the moment eps_b0 is corrected.  The literal here is the DERIVED
  // value 2.7990 (eps_meas = -0.0070024 above), not the 2.8000 the docs'
  // tables get from the ROUNDED eps_b0 = -0.0070; both round to the same
  // "2.80 GeV^2", which is all three figures of the rounded input support.
  // (This read 2.7991 until 2026-09-05 -- inside the 1e-3 tolerance, but one
  // unit wrong in the digit it printed.)
  CHECK_CLOSE(meas.t_positivity_edge(-2.0), 2.7990, 1e-3);
  CoherentScenario mod = sc;
  mod.eps_b0 = eps_mod;
  CHECK_CLOSE(mod.t_positivity_edge(-2.0), 0.37202, 1e-4);
  // (8) ... and the reason that DOES survive a change of eps_b0 is the
  // anchor range: `mantysaari_a2_deuteron` is digitized to |t| <= 0.30, and
  // 0.2 is inside it while every one of the three quadrupole rows above puts
  // the positivity edge somewhere else entirely.
  double t_anchor_max = 0.0;
  for (const MantysaariRow& r : mantysaari_a2_deuteron())
    t_anchor_max = std::max(t_anchor_max, r.t_abs);
  CHECK(t_anchor_max == 0.30);
  CHECK(COHERENT_T_MAX_DEFAULT < t_anchor_max);
  for (double e : {sc.eps_b0, eps_mod, eps_gfmc, eps_meas}) {
    CoherentScenario v = sc;
    v.eps_b0 = e;
    CHECK(COHERENT_T_MAX_DEFAULT < v.t_positivity_edge(-2.0));
  }
  CHECK(meas.t_positivity_edge(-2.0) > 9.0 * t_anchor_max);
}

// D5.  `t_positivity_edge` in its own right: it must reproduce the zero of
// `positivity_margin` for every sign combination, not only the shipped one,
// and it must SAY what happens in the degenerate cases rather than return a
// plausible-looking number.
TEST_CASE("coherent: t_positivity_edge is the zero of positivity_margin") {
  const double inf = std::numeric_limits<double>::infinity();
  for (double eps : {-0.08, -0.0070, +0.05, +0.13}) {
    for (double b : {33.1, 50.0}) {
      for (double amp : {0.0, 0.01, -0.02}) {
        for (double pzz : {-2.0, -0.6, 0.6, 1.0}) {
          CoherentScenario v;
          v.eps_b0 = eps;
          v.slope_b = b;
          v.amp = amp;
          const double edge = v.t_positivity_edge(pzz);
          REQUIRE(std::isfinite(edge));
          REQUIRE(edge > 0.0);
          // AT the edge the margin is zero, and it is positive strictly
          // inside and negative strictly outside -- i.e. this is the FIRST
          // crossing, which is what a ceiling has to be.
          CHECK_CLOSE_AT(v.positivity_margin(edge, pzz), 0.0, 0.0, 1e-12);
          CHECK(v.positivity_margin(0.999 * edge, pzz) > 0.0);
          CHECK(v.positivity_margin(1.001 * edge, pzz) < 0.0);
        }
      }
    }
  }
  // Degenerate: c_2 constant in |t| (P_zz = 0, or eps_b0 = 0, or B = 0) --
  // the weight never stops being a density, so the edge is infinite.
  CoherentScenario z;
  CHECK(z.t_positivity_edge(0.0) == inf);
  CoherentScenario flat = z;
  flat.eps_b0 = 0.0;
  CHECK(flat.t_positivity_edge(-2.0) == inf);
  // ... unless the FLAT term alone already exceeds 1, in which case the
  // weight is not a density anywhere and the edge is 0, not infinity.
  flat.amp = 0.6;
  CHECK(flat.t_positivity_edge(-2.0) == 0.0);
  CoherentScenario over = z;
  over.amp = 0.75;                       // |amp P_zz| = 1.5 at P_zz = -2
  CHECK(over.t_positivity_edge(-2.0) == 0.0);
  CHECK(over.positivity_margin(0.0, -2.0) < 0.0);
}
