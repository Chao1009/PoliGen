// Coherent (intact-6Li) channel: the scenario model, the analytic tagging
// identities, the m-state relation of the deformation mechanism, the recoil
// kinematics, and the breakup / veto table routing -- against
// validation/reference/coherent.json and evgen/tests/test_coherent.py.

#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "check_close.hpp"
#include "doctest.h"
#include "json_min.hpp"
#include "lipolgen/beams.hpp"
#include "lipolgen/coherent.hpp"
#include "lipolgen/constants.hpp"
#include "lipolgen/event.hpp"
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
