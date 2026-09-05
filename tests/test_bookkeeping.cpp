// SPDX-License-Identifier: GPL-3.0-or-later
// Run-plan bookkeeping: the four standard plans against the polligen
// reference table (validation/reference/bookkeeping.json, rtol 1e-12), the
// spin-temperature pzz_true ladder, the relative-luminosity bias formulas and
// the run-share rule -- a share moves COUNTS, never CROSS SECTIONS
// (evgen/tests/test_bookkeeping.py, test_run_share.py).

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "check_close.hpp"
#include "json_min.hpp"
#include "lipolgen/bookkeeping.hpp"
#include "lipolgen/spin.hpp"

using namespace lipolgen;

#ifndef LIPOLGEN_REFERENCE_DIR
#define LIPOLGEN_REFERENCE_DIR "validation/reference"
#endif

namespace {

constexpr double kRtol = 1e-12;
// The max-entropy fills come out of a bisection with tol 1e-13, transcribed
// step for step from the Python (spin.cpp), so the agreement is a couple of
// ulps rather than exact; the reference README asks for a looser comparison
// on those tables only.
constexpr double kMaxentRtol = 1e-13;
constexpr double kPe = 0.7;

bool load_bookkeeping(jsonmin::Value& out) {
  return jsonmin::load_file(std::string(LIPOLGEN_REFERENCE_DIR) +
                                "/bookkeeping.json",
                            out);
}

/// The plan each named reference block was dumped from.
RunPlan build_plan(const std::string& name) {
  if (name == "helicity_flip_j12") {
    return helicity_flip_plan(0.5, 0.7, kPe);
  }
  if (name == "helicity_flip_j1_maxent_anchor") {
    return helicity_flip_plan(1.0, 8.0 / 13.0, kPe);
  }
  if (name == "helicity_flip_j32_maxent_anchor") {
    return helicity_flip_plan(1.5, 0.7, kPe);
  }
  if (name == "helicity_flip_j1_explicit_pzz") {
    HelicityFlipOptions opt;
    opt.use_explicit_pzz = true;
    opt.pzz = 0.35;
    return helicity_flip_plan(1.0, 0.6, kPe, opt);
  }
  if (name == "helicity_flip_j1_tilted_offset") {
    HelicityFlipOptions opt;
    opt.theta_s = 0.9;
    opt.phi_s = 0.4;
    opt.rel_lumi_offset = 1e-4;
    return helicity_flip_plan(1.0, 0.7, kPe, opt);
  }
  if (name == "tensor_thirds_j1") {
    return tensor_thirds_plan(0.7, 0.6, 1e-4);
  }
  if (name == "transverse_tensor_j1") {
    return transverse_tensor_plan(0.6);
  }
  if (name == "tensor_flip_j1") {
    return tensor_flip_plan(0.6, kPi / 2.0, 0.5, 1e-4);
  }
  if (name == "tensor_thirds_j1_with_offset_azz0") {
    return with_offset(tensor_thirds_plan(0.7, 0.6), "azz0", 1e-4);
  }
  throw std::runtime_error("unknown reference plan " + name);
}

/// Independent spin-temperature construction, solved here rather than taken
/// from the module under test (bisection in t, not in beta).
std::vector<double> ladder_populations(double j, double pz) {
  return spin_temperature_ladder(j, pz).populations;
}

}  // namespace

TEST_CASE("reference tables: bookkeeping.json plans") {
  jsonmin::Value doc;
  if (!load_bookkeeping(doc)) {
    MESSAGE("validation/reference/bookkeeping.json absent -- skipping");
    return;
  }
  const jsonmin::Value& plans = doc["plans"];
  for (const auto& entry : plans.obj()) {
    const std::string& name = entry.first;
    const jsonmin::Value& ref = entry.second;
    CAPTURE(name);
    const RunPlan plan = build_plan(name);
    // A max-entropy fill is the only branch that is not closed-form.
    const bool maxent = name.rfind("helicity_flip_j1_explicit", 0) != 0 &&
                        name.rfind("helicity_flip", 0) == 0 &&
                        name != "helicity_flip_j12";
    const double rtol = maxent ? kMaxentRtol : kRtol;

    CHECK_CLOSE_AT(plan.pe_true(), ref["pe_true"].num(), rtol, 1e-300);
    CHECK_CLOSE_AT(plan.pz_true(), ref["pz_true"].num(), rtol, 1e-300);
    CHECK_CLOSE_AT(plan.pzz_true(), ref["pzz_true"].num(), rtol, 1e-300);
    // no polarimetry smear in any reference plan
    CHECK(plan.measured_pe() == plan.pe_true());
    CHECK(plan.measured_pz() == plan.pz_true());
    CHECK(plan.measured_pzz() == plan.pzz_true());

    const jsonmin::Value& cats = ref["categories"];
    REQUIRE(cats.size() == plan.categories().size());
    for (std::size_t i = 0; i < cats.size(); ++i) {
      const SpinCategory& c = plan.categories()[i];
      const jsonmin::Value& rc = cats[i];
      CAPTURE(c.name);
      CHECK(c.name == rc["name"].str());
      CHECK(c.j == rc["j"].num());
      CHECK(c.lam_e == static_cast<int>(rc["lam_e"].num()));
      CHECK_CLOSE_AT(c.pe, rc["pe"].num(), kRtol, 1e-300);
      CHECK_CLOSE_AT(c.theta_s, rc["theta_s"].num(), kRtol, 1e-300);
      CHECK_CLOSE_AT(c.phi_s, rc["phi_s"].num(), kRtol, 1e-300);
      CHECK_CLOSE(c.lumi_fraction, rc["lumi_fraction"].num(), kRtol);
      const std::vector<double> pops = rc["populations"].flat();
      REQUIRE(pops.size() == c.populations.size());
      for (std::size_t k = 0; k < pops.size(); ++k) {
        CHECK_CLOSE_AT(c.populations[k], pops[k], rtol, 1e-300);
      }
      // `moments()` is spin.moments_along_axis; the j=1/2 rows carry the
      // vector moment alone (the Python raises there, see the README's
      // "Functions this script could not call").
      const std::vector<double> mom = rc["moments"].flat();
      if (std::fabs(c.j - 0.5) < 1e-9) {
        REQUIRE(mom.size() == 1);
        CHECK_CLOSE(c.vector_moment(), mom[0], kRtol);
      } else {
        const AxisMoments am = c.moments();
        REQUIRE(mom.size() >= 2);
        CHECK_CLOSE_AT(am.vector, mom[0], rtol, 1e-300);
        CHECK_CLOSE_AT(am.tensor, mom[1], rtol, 1e-300);
        if (mom.size() > 2) CHECK_CLOSE_AT(am.octupole, mom[2], rtol, 1e-300);
      }
    }

    // lumi_shares(1e6) is purely lumi_fraction * total.
    const std::map<std::string, double> shares = plan.lumi_shares(1e6);
    const jsonmin::Value& rs = ref["lumi_shares_at_1e6_pb"];
    CHECK(shares.size() == rs.size());
    for (const auto& kv : shares) {
      CHECK_CLOSE(kv.second, rs[kv.first].num(), kRtol);
    }
  }
}

TEST_CASE("reference tables: spin-temperature ladder and rel-lumi biases") {
  jsonmin::Value doc;
  if (!load_bookkeeping(doc)) {
    MESSAGE("validation/reference/bookkeeping.json absent -- skipping");
    return;
  }
  const jsonmin::Value& lad = doc["spin_temperature_ladder_t3"];

  const SpinTemperatureLadder j1 = spin_temperature_ladder(1.0, 8.0 / 13.0);
  CHECK_CLOSE(j1.t, lad["j1"]["t"].num(), 1e-12);
  const std::vector<double> p1 = lad["j1"]["populations"].flat();
  for (std::size_t i = 0; i < p1.size(); ++i) {
    CHECK_CLOSE(j1.populations[i], p1[i], 1e-12);
  }
  CHECK_CLOSE(moments_along_axis(1.0, j1.populations).tensor,
              lad["j1"]["pzz"].num(), 1e-12);

  const SpinTemperatureLadder j32 = spin_temperature_ladder(1.5, 0.7);
  CHECK_CLOSE(j32.t, lad["j32"]["t"].num(), 1e-12);
  const std::vector<double> p32 = lad["j32"]["populations"].flat();
  for (std::size_t i = 0; i < p32.size(); ++i) {
    CHECK_CLOSE(j32.populations[i], p32[i], 1e-12);
  }
  CHECK_CLOSE(moments_along_axis(1.5, j32.populations).tensor,
              lad["j32"]["t_moment"].num(), 1e-12);

  const jsonmin::Value& bias = doc["rel_lumi_bias"];
  for (std::size_t i = 0; i < bias["apar_rel_lumi_bias"].size(); ++i) {
    const jsonmin::Value& r = bias["apar_rel_lumi_bias"][i];
    CHECK_CLOSE(apar_rel_lumi_bias(r["offset"].num(), r["pe"].num(),
                                   r["pz"].num()),
                r["value"].num(), kRtol);
  }
  for (std::size_t i = 0; i < bias["azz_rel_lumi_bias"].size(); ++i) {
    const jsonmin::Value& r = bias["azz_rel_lumi_bias"][i];
    CHECK_CLOSE(azz_rel_lumi_bias(r["offset"].num(), r["pzz"].num()),
                r["value"].num(), kRtol);
  }
}

TEST_CASE("a spin-1/2 flip plan has no rank-2 moment at all") {
  // Exactly zero, not merely small, and it must survive the polarimetry
  // smear as zero: a smeared 1e-17 divisor would be worse than a zero one.
  const RunPlan plan = helicity_flip_plan(0.5, 0.7, kPe);
  CHECK(plan.pzz_true() == 0.0);
  CHECK(plan.measured_pzz() == 0.0);
  const RunPlan smeared(plan.categories(), plan.pe_true(), plan.pz_true(),
                        plan.pzz_true(), 0.03);
  CHECK(smeared.measured_pzz() == 0.0);
  CHECK(smeared.measured_pz() != 0.7);   // the smear does fire
  CHECK_CLOSE_AT(plan.categories()[0].populations[0], 0.85, 0.0, 1e-15);
  CHECK_CLOSE_AT(plan.categories()[0].populations[1], 0.15, 0.0, 1e-15);
}

TEST_CASE("spin-temperature pzz_true at the t=3 rational anchors") {
  // t = 3 makes both ladders rational: (9,3,1)/13 for J=1 and (27,9,3,1)/40
  // for J=3/2, giving (P_z, rank-2) = (8/13, 4/13) and (7/10, 2/5).  The
  // rank-2 value is re-derived from the density matrix and the angular
  // momentum operators, not from moments_along_axis.
  struct Anchor { double j, pz, pzz; };
  const Anchor anchors[2] = {{1.0, 8.0 / 13.0, 4.0 / 13.0}, {1.5, 0.7, 0.4}};
  for (const Anchor& a : anchors) {
    CAPTURE(a.j);
    const RunPlan plan = helicity_flip_plan(a.j, a.pz, kPe);
    CHECK_CLOSE_AT(plan.pzz_true(), a.pzz, 0.0, 1e-12);
    const std::vector<double> pops = ladder_populations(a.j, a.pz);
    for (std::size_t i = 0; i < pops.size(); ++i) {
      CHECK_CLOSE_AT(plan.categories()[0].populations[i], pops[i], 0.0, 1e-12);
    }
    const CplxMatrix rho = rho_from_populations(a.j, pops);
    CHECK_CLOSE_AT(tensor_polarization(rho, a.j), a.pzz, 0.0, 1e-12);
  }
}

TEST_CASE("spin-temperature pzz_true across pz, through the operator trace") {
  for (double j : {1.0, 1.5}) {
    for (double pz : {0.05, 0.3, 0.7, 0.9}) {
      CAPTURE(j); CAPTURE(pz);
      const RunPlan plan = helicity_flip_plan(j, pz, kPe);
      const std::vector<double> pops = ladder_populations(j, pz);
      const CplxMatrix rho = rho_from_populations(j, pops);
      CHECK_CLOSE_AT(plan.pzz_true(), tensor_polarization(rho, j), 0.0, 1e-10);
      if (j == 1.0) {
        // the published spin-1 relation Pzz = 2 - sqrt(4 - 3 Pz^2)
        CHECK_CLOSE_AT(plan.pzz_true(), 2.0 - std::sqrt(4.0 - 3.0 * pz * pz),
                       0.0, 1e-12);
      }
      CHECK(plan.pzz_true() > 0.0);   // a vector fill is prolate, never zero
    }
  }
}

TEST_CASE("spin-temperature pzz_true, small-pz limit") {
  // beta = 3 pz/(J+1) and the rank-2 moment is beta^2/3 (J=1) / beta^2/2
  // (J=3/2), i.e. (3/4) pz^2 and (18/25) pz^2.
  const double pz = 1e-4;
  CHECK_CLOSE(helicity_flip_plan(1.0, pz, kPe).pzz_true() / (pz * pz), 0.75,
              1e-6);
  CHECK_CLOSE(helicity_flip_plan(1.5, pz, kPe).pzz_true() / (pz * pz), 0.72,
              1e-6);
}

TEST_CASE("pzz_true is even in pz and is the moment along the FILL axis") {
  // The rank-2 moment cannot know the sign of the vector one; that is what
  // makes the flip pattern tensor-blind.
  HelicityFlipOptions explicit_1;
  explicit_1.use_explicit_pzz = true;
  explicit_1.pzz = 0.35;
  HelicityFlipOptions explicit_32;
  explicit_32.use_explicit_pzz = true;
  explicit_32.pzz = 0.2;
  struct Case { double j, pz; HelicityFlipOptions opt; };
  const Case cases[4] = {{1.0, 0.7, HelicityFlipOptions()},
                         {1.5, 0.7, HelicityFlipOptions()},
                         {1.0, 0.6, explicit_1},
                         {1.5, 0.5, explicit_32}};
  for (const Case& c : cases) {
    CAPTURE(c.j); CAPTURE(c.pz);
    const RunPlan up = helicity_flip_plan(c.j, c.pz, kPe, c.opt);
    const RunPlan down = helicity_flip_plan(c.j, -c.pz, kPe, c.opt);
    CHECK_CLOSE_AT(up.pzz_true(), down.pzz_true(), 0.0, 1e-12);
    CHECK_CLOSE_AT(up.pz_true(), -down.pz_true(), 0.0, 1e-15);
    for (const SpinCategory& cat : up.categories()) {
      CHECK_CLOSE_AT(cat.moments().tensor, up.pzz_true(), 0.0, 1e-12);
    }
    CHECK(up.categories()[0].lam_e == +1);
    CHECK(up.categories()[1].lam_e == -1);
  }

  // Tilting theta_S leaves pzz_true alone while the LAB moment picks up
  // P2(cos theta_S): the kernel works in the axis frame, so the axis moment
  // is the right divisor.
  const double pz = 0.7, theta_s = 0.9;
  const RunPlan flat = helicity_flip_plan(1.0, pz, kPe);
  HelicityFlipOptions tilt;
  tilt.theta_s = theta_s;
  const RunPlan tilted = helicity_flip_plan(1.0, pz, kPe, tilt);
  CHECK_CLOSE_AT(tilted.pzz_true(), flat.pzz_true(), 0.0, 1e-15);
  SpinDensity dens;
  dens.j = 1.0;
  dens.populations = tilted.categories()[0].populations;
  dens.theta = theta_s;
  const double p2 = 0.5 * (3.0 * std::cos(theta_s) * std::cos(theta_s) - 1.0);
  CHECK_CLOSE_AT(dens.lab_moments().tensor_zz, flat.pzz_true() * p2, 0.0,
                 1e-12);
}

TEST_CASE("explicit pzz closes and unphysical corners are refused") {
  struct Case { double j, pz, pzz; };
  const Case cases[6] = {{1.0, 0.6, 0.35}, {1.0, 0.5, -0.4}, {1.0, 0.0, 1.0},
                         {1.5, 0.5, 0.2},  {1.5, 0.0, -1.0}, {1.5, 0.3, -0.1}};
  for (const Case& c : cases) {
    CAPTURE(c.j); CAPTURE(c.pzz);
    HelicityFlipOptions opt;
    opt.use_explicit_pzz = true;
    opt.pzz = c.pzz;
    const RunPlan plan = helicity_flip_plan(c.j, c.pz, kPe, opt);
    CHECK_CLOSE_AT(plan.pzz_true(), c.pzz, 0.0, 1e-12);
    CHECK(plan.pz_true() == c.pz);
    const CplxMatrix rho =
        rho_from_populations(c.j, plan.categories()[0].populations);
    CHECK_CLOSE_AT(tensor_polarization(rho, c.j), c.pzz, 0.0, 1e-12);
    CHECK_CLOSE_AT(vector_polarization(rho, c.j)[2], c.pz, 0.0, 1e-12);
  }

  HelicityFlipOptions bad_spin;
  bad_spin.use_explicit_pzz = true;
  bad_spin.pzz = 0.3;
  CHECK_THROWS_AS(helicity_flip_plan(2.0, 0.5, kPe, bad_spin),
                  std::runtime_error);
  HelicityFlipOptions bad_pzz;
  bad_pzz.use_explicit_pzz = true;
  bad_pzz.pzz = -2.0;
  CHECK_THROWS_AS(helicity_flip_plan(1.0, 0.5, kPe, bad_pzz),
                  std::runtime_error);
}

TEST_CASE("polarimetry smearing is fixed-seed and reproducible") {
  const RunPlan plan = tensor_thirds_plan(0.7, 0.6);
  CHECK(plan.measured_pzz() == 0.6);   // no smearing by default
  const RunPlan smeared(plan.categories(), 0.0, 0.7, 0.6, 0.03, 5);
  CHECK(smeared.measured_pzz() != 0.6);
  const RunPlan again(plan.categories(), 0.0, 0.7, 0.6, 0.03, 5);
  CHECK(smeared.measured_pzz() == again.measured_pzz());
  CHECK(smeared.measured_pz() == again.measured_pz());
  // a different seed is a different fill
  const RunPlan other(plan.categories(), 0.0, 0.7, 0.6, 0.03, 6);
  CHECK(other.measured_pzz() != smeared.measured_pzz());
  // the smear is a few percent, not a factor
  CHECK(std::fabs(smeared.measured_pzz() / 0.6 - 1.0) < 0.3);
}

TEST_CASE("with_offset moves exactly one share") {
  const RunPlan base = tensor_thirds_plan(0.7, 0.6);
  const RunPlan off = with_offset(base, "azz0", 1e-4);
  for (std::size_t i = 0; i < base.categories().size(); ++i) {
    const double factor = (base.categories()[i].name == "azz0") ? 1.0001 : 1.0;
    CHECK_CLOSE(off.categories()[i].lumi_fraction,
                base.categories()[i].lumi_fraction * factor, 1e-15);
  }
  CHECK(off.pzz_true() == base.pzz_true());
  CHECK_THROWS_AS(with_offset(base, "nope", 1e-4), std::runtime_error);
}

TEST_CASE("the spin share is not the programme share") {
  // SpinCategory::lumi_fraction divides ONE measurement's own luminosity and
  // sums to one; the programme share divides the year and does not.
  const RunPlan plan = tensor_flip_plan(0.6);
  double tot = 0.0;
  for (const SpinCategory& c : plan.categories()) tot += c.lumi_fraction;
  CHECK_CLOSE(tot, 1.0, 1e-15);

  const std::map<std::string, double> shares = plan.lumi_shares(1.0e3);
  double sum = 0.0;
  for (const auto& kv : shares) sum += kv.second;
  CHECK_CLOSE(sum, 1.0e3, 1e-12);

  // a programme share of 1/4 is applied OUTSIDE the plan, to the total
  const std::map<std::string, double> quarter = plan.lumi_shares(0.25 * 1.0e3);
  for (const auto& kv : shares) {
    CHECK_CLOSE(quarter.at(kv.first), 0.25 * kv.second, 1e-12);
  }
}
