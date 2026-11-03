// SPDX-License-Identifier: GPL-3.0-or-later
// plans/10: the EIC energy menu, and that nothing has drifted off it.
// Ported from PolarizedLithiumSim/fastsim/tests/test_beams.py.
//
// Ion energies are GAMMA-matched to their proton configuration and then
// capped by the ring rigidity -- not rigidity-scaled.  Nothing below is a
// hard-coded rounded literal: every energy is derived from the physical
// (AME2020) nuclear mass and the two constraints.

#include <cmath>
#include <string>
#include <vector>

#include "check_close.hpp"
#include "lipolgen/beams.hpp"
#include "lipolgen/constants.hpp"

using namespace lipolgen;

TEST_CASE("a proton reproduces the machine configurations exactly") {
  // The strongest self-check available: run the ion machinery on a PROTON and
  // it must return the configuration energies themselves.
  const std::vector<BeamConfig> cfgs = default_configs("p");
  REQUIRE(cfgs.size() == 3);
  CHECK(cfgs[0].ion_momentum_per_nucleon == 41.0);
  CHECK(cfgs[1].ion_momentum_per_nucleon == 100.0);
  CHECK(cfgs[2].ion_momentum_per_nucleon == 275.0);
}

TEST_CASE("ions are gamma matched below the rigidity cap") {
  // Every species sits at ~41 GeV/u at the low configuration -- the same
  // SPEED as the 41 GeV proton, not the same rigidity.
  for (const char* name : {"d", "3He", "6Li", "7Li"}) {
    const Ion& ion = ion_by_name(name);
    const double lo = default_configs(name)[0].ion_momentum_per_nucleon;
    CHECK(lo >= 40.5);
    CHECK(lo <= 41.1);
    // and NOT the rigidity-scaled value, which for 6Li is 20.5
    const double rigidity_scaled = 41.0 * ion.Z / static_cast<double>(ion.A);
    CHECK((std::fabs(lo - rigidity_scaled) > 1.0 || ion.Z == ion.A));
  }
}

TEST_CASE("top energies are rigidity capped and match the published menu") {
  const std::vector<std::pair<std::string, double>> want = {
      {"d", 137.5}, {"3He", 183.3}, {"6Li", 137.5}, {"7Li", 117.9}};
  for (const auto& w : want) {
    const double got = default_configs(w.first)[2].ion_momentum_per_nucleon;
    CHECK_CLOSE_AT(got, w.second, 0.0, 0.15);
    // the top point is the rigidity cap, i.e. it is NOT gamma-matched
    CHECK_CLOSE_AT(got, ion_by_name(w.first).momentum_per_nucleon_max(), 0.0,
                   0.05);
  }
}

TEST_CASE("the 6Li and 7Li menus are 40.8/99.5/137.5 and 40.8/99.5/117.9") {
  const double want6[3] = {40.8, 99.5, 137.5};
  const double want7[3] = {40.8, 99.5, 117.9};
  const std::vector<BeamConfig> c6 = default_configs("6Li");
  const std::vector<BeamConfig> c7 = default_configs("7Li");
  for (int i = 0; i < 3; ++i) {
    CHECK(c6[static_cast<std::size_t>(i)].ion_momentum_per_nucleon == want6[i]);
    CHECK(c7[static_cast<std::size_t>(i)].ion_momentum_per_nucleon == want7[i]);
    CHECK(c6[static_cast<std::size_t>(i)].electron_energy
          == ELECTRON_ENERGIES[i]);
  }
  // ... and the mid configuration is the one every kernel test runs on
  CHECK_CLOSE(c6[1].s_per_nucleon(), 3980.0000000000005, 1e-15);
  CHECK_CLOSE(c6[1].sqrt_s_per_nucleon(), 63.08724118235002, 1e-15);
}

TEST_CASE("mass per nucleon uses the physical nuclear mass") {
  // Binding matters at the 0.5 % level, which is exactly what separates a
  // 41 GeV proton from a gamma-matched 6Li at 40.8 GeV/u.
  CHECK_CLOSE(nucleus_mass("6Li", 6, 3), 5.601518702, 1e-15);
  CHECK_CLOSE(LI6().mass_per_nucleon(), 5.601518702 / 6.0, 1e-15);
  CHECK_CLOSE_AT(PROTON_MASS - LI6().mass_per_nucleon(), 4.7e-3, 0.0, 1e-3);
  CHECK(LI6().mass_per_nucleon() < 0.93149 * 1.01);  // not A * M_U
  CHECK_THROWS(nucleus_mass("12C", 12, 6));
}

TEST_CASE("gamma_of and the EPIOS synchronisation windows") {
  // Yellow Report Table 10.2's gold at 41 GeV/u is what settles gamma
  // matching; the 41 and 275 GeV anchors sit inside the two EPIOS windows
  // and the 100 GeV point (gamma 106.6) sits in neither -- a KNOWN CONFLICT
  // the Python flags rather than resolves.
  CHECK_CLOSE(gamma_of(41.0), 43.70878676863391, 1e-14);
  CHECK_CLOSE(gamma_of(100.0), 106.58358375431271, 1e-14);
  CHECK_CLOSE(gamma_of(275.0), 293.0936603113653, 1e-14);
  CHECK(epios_window_of(41.0) == "bypass");
  CHECK(epios_window_of(275.0) == "radial-shift");
  CHECK(epios_window_of(100.0) == "");
}

TEST_CASE("the effective polarizations are per nucleon") {
  // 7Li: the verified whole-nucleus VMC sums P_p = +0.866, P_n = -0.037 are
  // stored divided by Z = 3 and N = 4, so Z*P_p and N*P_n return them.
  CHECK_CLOSE(LI7().Z * LI7().eff_pol_p, 0.866, 1e-14);
  CHECK_CLOSE(LI7().N() * LI7().eff_pol_n, -0.037, 1e-14);
  // 6Li: the slots hold LI6_CLUSTER_POLARIZATION/3 each, so Z*P_p = N*P_n is
  // the whole-nucleus factor.  A1 (2026-11-03): 0.81123 (the cluster product,
  // now LI6_CLUSTER_POLARIZATION_HULTHEN) -> the six-body VMC 0.848; the
  // Hulthen Ion is `li6_ion(LI6_CLUSTER_POLARIZATION_HULTHEN)`.
  CHECK_CLOSE(LI6().Z * LI6().eff_pol_p, LI6_CLUSTER_POLARIZATION, 1e-15);
  CHECK_CLOSE(LI6().N() * LI6().eff_pol_n, LI6_CLUSTER_POLARIZATION, 1e-15);
  CHECK_CLOSE(LI6_CLUSTER_POLARIZATION_HULTHEN, 0.81123, 5e-6);  // A1 (2026-11-03): the pre-A1 constant, named; the default is 0.848 (end of file)
  CHECK(LI6().spin == 1.0);
  CHECK(LI7().spin == 1.5);
}

TEST_CASE("the 6Li slot is the cluster wave function, nothing hard-coded") {
  // The 6Li slot is built from the wave function the tagged sector uses, not
  // from a transcribed 0.81 (author decision 2026-08-29, plans/04 #6).
  CHECK(P_D_LI6 == 0.0867);
  CHECK(P_D_DEUTERON_SCENARIO == 0.045);   // A1 (2026-11-03): the scenario, named
  CHECK(LI6_CLUSTER_POLARIZATION_HULTHEN   // A1: the pre-A1 constant, bit for bit
        == (1.0 - 1.5 * P_D_LI6) * (1.0 - 1.5 * P_D_DEUTERON_SCENARIO));
  // the SCENARIO deuteron carries the same expression, bit for bit, which is
  // what makes the Hulthen per-nucleon g1 ratio below exact (A1: re-pointed)
  CHECK(DEUTERON_SCENARIO().eff_pol_p == 1.0 - 1.5 * P_D_DEUTERON_SCENARIO);
  CHECK(DEUTERON().eff_pol_n == DEUTERON().eff_pol_p);
  // Hulthen per-nucleon g1(6Li)/g1(d) = (1 - 1.5 P_D_LI6)/3 = 0.290 (A1: the
  // explicit pre-A1 ions; the defaults' ratio is pinned at the end of the file)
  const double ratio = li6_ion(LI6_CLUSTER_POLARIZATION_HULTHEN).eff_pol_p / DEUTERON_SCENARIO().eff_pol_p;
  CHECK_CLOSE(ratio, (1.0 - 1.5 * P_D_LI6) / 3.0, 1e-15);
  CHECK_CLOSE(ratio, 0.290, 5e-4);
  // the retired Cloet convention, pinned: a whole-nucleus 1.0, 1.233 times
  // the cluster picture and the g1(6Li) every number published before
  // 2026-08-29 was computed with
  CHECK(LI6_NAIVE_ONE_THIRD == 1.0 / 3.0);
  CHECK_CLOSE(LI6_NAIVE_ONE_THIRD / (LI6_CLUSTER_POLARIZATION_HULTHEN / 3.0), 1.233,
              5e-4);   // A1 (2026-11-03): against the Hulthen reading; 1.179 against 0.848
}

TEST_CASE("beam config labels and sqrt(s)") {
  const BeamConfig c = default_configs("6Li")[1];
  CHECK(c.label() == "e(10) x 6Li(99.5/u)");
  CHECK_CLOSE(c.sqrt_s_per_nucleon(), std::sqrt(4.0 * 10.0 * 99.5), 1e-15);
}

// A1 (2026-11-03; docs/open_items/run_2026-11-03/DECISIONS.md, the
// supervisor's choice under the author's delegation): the three inclusive
// factors, their defaults, their named pre-A1 values and the parameters.
TEST_CASE("A1: the inclusive polarization factors and their parameters") {
  // (a) the vector factor IS the six-body constant, one home
  CHECK(LI6_CLUSTER_POLARIZATION == LI6_POLARIZATION_VMC_SIX_BODY);
  CHECK(LI6().eff_pol_p == LI6_CLUSTER_POLARIZATION / 3.0);
  // ... the Hulthen reading is still exactly the pre-A1 0.811228375
  CHECK(LI6_CLUSTER_POLARIZATION_HULTHEN
        == li6_cluster_polarization(P_D_LI6, P_D_DEUTERON_SCENARIO));
  CHECK(LI6_CLUSTER_POLARIZATION_HULTHEN == 0.811228375);   // beams.json's value
  // ... and `li6_ion` at it is the pre-A1 `LI6()` slot, bit for bit
  CHECK(li6_ion(LI6_CLUSTER_POLARIZATION_HULTHEN).eff_pol_p == 0.811228375 / 3.0);
  CHECK(li6_ion().eff_pol_p == LI6().eff_pol_p);
  // the band: the two cluster readings bracket the default
  CHECK(LI6_CLUSTER_POLARIZATION_HULTHEN < LI6_CLUSTER_POLARIZATION);
  // (c) the inclusive deuteron: AV18's P_D, 0.9136; the scenario under a name
  CHECK(P_D_DEUTERON == 0.0576);
  CHECK(DEUTERON_VECTOR_POLARIZATION == 1.0 - 1.5 * 0.0576);
  CHECK_CLOSE(DEUTERON_VECTOR_POLARIZATION, 0.9136, 1e-15);
  CHECK(DEUTERON().eff_pol_p == DEUTERON_VECTOR_POLARIZATION);
  CHECK(DEUTERON_SCENARIO().eff_pol_p == 0.9325);
  CHECK(DEUTERON_VECTOR_POLARIZATION_SCENARIO == 0.9325);
  CHECK(deuteron_ion(P_D_DEUTERON_SCENARIO).eff_pol_p == DEUTERON_SCENARIO().eff_pol_p);
  CHECK(deuteron_ion().eff_pol_p == DEUTERON().eff_pol_p);
  // the default ions' per-nucleon g1 ratio, (0.848/3)/0.9136
  CHECK_CLOSE(LI6().eff_pol_p / DEUTERON().eff_pol_p, 0.3093987157, 1e-9);
  // (b) the transfer (its tagged-model pin is in tests/test_tagged.cpp)
  CHECK(LI6_B1_RANK2_TRANSFER == LI6_B1_RANK2_TRANSFER_VMC);
  CHECK(LI6_B1_RANK2_TRANSFER_VMC == 0.982576);
  CHECK(LI6_B1_RANK2_TRANSFER_HULTHEN == 0.921947);
}
