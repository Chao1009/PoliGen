// SPDX-License-Identifier: GPL-3.0-or-later
// Benchmark group SPIN AND DEUTERON, the C++ side of BENCHMARK_PLAN.md sec. 4
// rows 1 and 9 (the Python harnesses are validation/benchmarks/
// t5_epios_source_modes.py and t5_est_identity.py; the record is
// docs/open_items/run_2026-09-23/phase_B1_spin_deuteron.md).
//
// Row 1: the eight EPIOS Table II (P_z, P_zz) deuteron source modes, read from
// the VENDORED table validation/benchmarks/data/
// epios2026_table2_deuteron_source_modes.json (Cartesian, per ion), through
// spin1_populations and back through the density matrix -- non-negative,
// normalised, moments to 1e-12.
//
// Row 9: spin_temperature_pzz(1, P_z) IS the polarized-target equal-spin-
// temperature closed form P_zz = 2 - sqrt(4 - 3 P_z^2) (COMPASS 6LiD common
// spin temperature, Koivuniemi et al. SPIN 2004; Keller, Crabb, Day, NIM A
// 981 (2020) 164504), to the bound populations_maxent's bisection guarantees:
// 0.5 * 1e-13 * max|dP_zz/dbeta| (= 0.42609) + rounding = 2.31e-14.  MEASURED
// 2026-09-23: 1.138e-14 over 19999 interior points -- the 1e-14 first asked
// for does NOT hold on a dense grid.  At |P_z| = 1 the tree throws.

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "doctest.h"
#include "json_min.hpp"
#include "lipolgen/bookkeeping.hpp"
#include "lipolgen/spin.hpp"

using namespace lipolgen;

#ifndef LIPOLGEN_REFERENCE_DIR
#define LIPOLGEN_REFERENCE_DIR "validation/reference"
#endif

namespace {

// "p/q" or "p" -> double, as the vendored table writes the ideal columns.
double rational(const std::string& s) {
  const std::size_t slash = s.find('/');
  if (slash == std::string::npos) return std::stod(s);
  return std::stod(s.substr(0, slash)) / std::stod(s.substr(slash + 1));
}

const char* kEpios =
    "/../benchmarks/data/epios2026_table2_deuteron_source_modes.json";

}  // namespace

TEST_CASE("bench row 1: EPIOS Table II modes round-trip through spin1_populations") {
  jsonmin::Value doc;
  REQUIRE(jsonmin::load_file(std::string(LIPOLGEN_REFERENCE_DIR) + kEpios, doc));
  const jsonmin::Value& rows = doc["rows"];
  REQUIRE(rows.size() == 8);
  for (std::size_t i = 0; i < rows.size(); ++i) {
    const double pz = rational(rows[i]["pz_ideal"].str());
    const double pzz = rational(rows[i]["pzz_ideal"].str());
    const std::array<double, 3> p = spin1_populations(pz, pzz);
    const std::vector<double> pv(p.begin(), p.end());
    CHECK(p[0] >= 0.0);
    CHECK(p[1] >= 0.0);
    CHECK(p[2] >= 0.0);
    CHECK(std::fabs(p[0] + p[1] + p[2] - 1.0) <= 1e-15);
    const AxisMoments m = moments_along_axis(1.0, pv);
    CHECK(std::fabs(m.vector - pz) <= 1e-12);
    CHECK(std::fabs(m.tensor - pzz) <= 1e-12);
    const CplxMatrix rho = rho_from_populations(1.0, pv);
    CHECK(std::fabs(vector_polarization(rho, 1.0)[2] - pz) <= 1e-12);
    CHECK(std::fabs(tensor_polarization(rho, 1.0) - pzz) <= 1e-12);
    // every polarized mode sits on the domain edge (min p_m = 0)
    if (rows[i]["mode"].num() != 0.0) {
      CHECK(std::min(p[0], std::min(p[1], p[2])) == 0.0);
    }
  }
}

TEST_CASE("bench row 9: spin_temperature_pzz(1, .) is the EST closed form") {
  const double tol = 0.5 * 1e-13 * 0.4260945789955075 + 8 * 2.2e-16;
  const int n = 20001;
  double worst = 0.0;
  for (int i = 1; i < n - 1; ++i) {
    const double pz = -1.0 + 2.0 * i / (n - 1);
    const double d = std::fabs(spin_temperature_pzz(1.0, pz) -
                               (2.0 - std::sqrt(4.0 - 3.0 * pz * pz)));
    worst = std::max(worst, d);
  }
  MESSAGE("EST identity: max |diff| = " << worst << " over " << n - 2
                                        << " interior points; tol " << tol);
  CHECK(worst <= tol);
  CHECK_THROWS(spin_temperature_pzz(1.0, 1.0));
  CHECK_THROWS(spin_temperature_pzz(1.0, -1.0));
  // J = 3/2 (7Li): the ladder is the Brillouin function B_{3/2}(J ln t)
  for (double pz : {-0.9, -0.3, 0.2, 0.6, 0.95}) {
    const SpinTemperatureLadder l = spin_temperature_ladder(1.5, pz);
    const double y = 1.5 * std::log(l.t);
    const double a = 4.0 / 3.0, b = 1.0 / 3.0;
    CHECK(std::fabs(a / std::tanh(a * y) - b / std::tanh(b * y) - pz) <= 1e-14);
  }
}
