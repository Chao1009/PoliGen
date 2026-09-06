// SPDX-License-Identifier: GPL-3.0-or-later
// Gate 1 (DEVELOPMENT_PLAN section 4.1): rho-matrix moments, all axes, exact.
// Ported one-for-one from PolarizedLithiumSim/evgen/tests/test_spin.py.

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "check_close.hpp"
#include "lipolgen/spin.hpp"

using namespace lipolgen;
using lipolgen_test::close;

namespace {

constexpr double kPi2 = 3.14159265358979323846;

CplxMatrix identity(std::size_t n) {
  CplxMatrix e(n);
  for (std::size_t i = 0; i < n; ++i) e(i, i) = 1.0;
  return e;
}

}  // namespace

TEST_CASE("wigner_d for spin 1 is the known 3x3 matrix") {
  const double beta = 0.7;
  const RealMatrix d = wigner_d(1.0, beta);
  const double c = std::cos(beta), s = std::sin(beta);
  const double r2 = std::sqrt(2.0);
  // rows/cols ordered m = +1, 0, -1
  const double expected[3][3] = {
      {(1 + c) / 2, -s / r2, (1 - c) / 2},
      {s / r2, c, -s / r2},
      {(1 - c) / 2, s / r2, (1 + c) / 2},
  };
  for (std::size_t a = 0; a < 3; ++a) {
    for (std::size_t b = 0; b < 3; ++b) {
      CHECK_CLOSE_AT(d(a, b), expected[a][b], 0.0, 1e-12);
    }
  }
}

TEST_CASE("wigner_d for spin 3/2 reproduces its diagonal entries") {
  const double beta = 1.1;
  const RealMatrix d = wigner_d(1.5, beta);
  const double ch = std::cos(beta / 2.0);
  CHECK_CLOSE_AT(d(0, 0), ch * ch * ch, 0.0, 1e-12);
  CHECK_CLOSE_AT(d(1, 1), (3 * std::cos(beta) - 1) / 2 * ch, 0.0, 1e-12);
}

TEST_CASE("the rotation matrix is unitary for every spin") {
  for (double j : {0.5, 1.0, 1.5}) {
    const CplxMatrix u = rotation_matrix(j, 0.9, 2.3);
    const CplxMatrix uu = matmul(u, adjoint(u));
    const CplxMatrix e = identity(u.size());
    for (std::size_t a = 0; a < u.size(); ++a) {
      for (std::size_t b = 0; b < u.size(); ++b) {
        CHECK_CLOSE_AT(uu(a, b).real(), e(a, b).real(), 0.0, 1e-12);
        CHECK_CLOSE_AT(uu(a, b).imag(), 0.0, 0.0, 1e-12);
      }
    }
  }
}

TEST_CASE("Clebsch-Gordan reference values 1/sqrt(3) and sqrt(2/3)") {
  CHECK_CLOSE(clebsch_gordan(1, 1, 1, -1, 0, 0), 1.0 / std::sqrt(3.0), 1e-14);
  CHECK_CLOSE(clebsch_gordan(0.5, 0.5, 0.5, -0.5, 1, 0), 1.0 / std::sqrt(2.0),
              1e-14);
  CHECK_CLOSE(clebsch_gordan(1, 0, 1, 0, 2, 0), std::sqrt(2.0 / 3.0), 1e-14);
  // selection rules
  CHECK(clebsch_gordan(1, 1, 1, 1, 0, 0) == 0.0);
  CHECK(clebsch_gordan(1, 1, 1, -1, 3, 0) == 0.0);
}

TEST_CASE("multipole operators are orthonormal under the trace product") {
  for (double j : {1.0, 1.5}) {
    const int kmax = static_cast<int>(std::lround(2.0 * j));
    std::vector<std::pair<std::pair<int, int>, CplxMatrix>> ops;
    for (int k = 0; k <= kmax; ++k) {
      for (int q = -k; q <= k; ++q) {
        ops.emplace_back(std::make_pair(k, q), multipole_operator(j, k, q));
      }
    }
    for (const auto& o1 : ops) {
      for (const auto& o2 : ops) {
        const std::complex<double> v = trace(matmul(adjoint(o1.second), o2.second));
        const double want = (o1.first == o2.first) ? 1.0 : 0.0;
        CHECK_CLOSE_AT(v.real(), want, 0.0, 1e-12);
        CHECK_CLOSE_AT(v.imag(), 0.0, 0.0, 1e-12);
      }
    }
  }
}

TEST_CASE("pure and unpolarized moments") {
  // pure m = +1 along z
  CplxMatrix rho = rho_from_populations(1.0, {1.0, 0.0, 0.0});
  const std::array<double, 3> v = vector_polarization(rho, 1.0);
  CHECK_CLOSE_AT(v[0], 0.0, 0.0, 1e-12);
  CHECK_CLOSE_AT(v[1], 0.0, 0.0, 1e-12);
  CHECK_CLOSE_AT(v[2], 1.0, 0.0, 1e-12);
  CHECK_CLOSE_AT(tensor_polarization(rho, 1.0), 1.0, 0.0, 1e-12);

  // unpolarized: every moment vanishes, on any axis
  for (auto jn : {std::make_pair(1.0, 3), std::make_pair(1.5, 4)}) {
    const std::vector<double> pops(static_cast<std::size_t>(jn.second),
                                   1.0 / jn.second);
    rho = rho_from_populations(jn.first, pops, 0.4, 1.2);
    const std::array<double, 3> vv = vector_polarization(rho, jn.first);
    for (double c : vv) CHECK_CLOSE_AT(c, 0.0, 0.0, 1e-12);
    CHECK_CLOSE_AT(tensor_polarization(rho, jn.first), 0.0, 0.0, 1e-12);
  }

  // HERMES-style m0-enriched fill: Pzz = -2
  rho = rho_from_populations(1.0, {0.0, 1.0, 0.0});
  CHECK_CLOSE_AT(tensor_polarization(rho, 1.0), -2.0, 0.0, 1e-12);
}

TEST_CASE("spin-1 populations round trip through the moments") {
  const double pz = 0.5, pzz = 0.3;
  const std::array<double, 3> pops = spin1_populations(pz, pzz);
  const AxisMoments m = moments_along_axis(1.0, {pops[0], pops[1], pops[2]});
  CHECK_CLOSE_AT(m.vector, pz, 0.0, 1e-12);
  CHECK_CLOSE_AT(m.tensor, pzz, 0.0, 1e-12);
  CHECK_THROWS(spin1_populations(1.0, -2.0));  // unphysical corner
}

TEST_CASE("spin-3/2 populations round trip through the moments") {
  const double pz = 0.6, t = 0.2, o = -0.1;
  const std::array<double, 4> pops = spin32_populations(pz, t, o);
  const AxisMoments m =
      moments_along_axis(1.5, {pops[0], pops[1], pops[2], pops[3]});
  CHECK_CLOSE_AT(m.vector, pz, 0.0, 1e-12);
  CHECK_CLOSE_AT(m.tensor, t, 0.0, 1e-12);
  CHECK(m.has_octupole);
  CHECK_CLOSE_AT(m.octupole, o, 0.0, 1e-12);
  // pure m = +3/2 has all normalized moments = +1
  const std::array<double, 4> pure = spin32_populations(1.0, 1.0, 1.0);
  const double want[4] = {1.0, 0.0, 0.0, 0.0};
  for (std::size_t i = 0; i < 4; ++i) CHECK_CLOSE_AT(pure[i], want[i], 0.0, 1e-12);
}

TEST_CASE("rotated moments follow the analytic n_hat / P2(cos theta) law") {
  const double theta = 0.8, phi = 2.1;
  const std::vector<std::pair<double, std::vector<double>>> cases = {
      {1.0, {0.6, 0.3, 0.1}},
      {1.5, {0.5, 0.25, 0.15, 0.1}},
  };
  for (const auto& c : cases) {
    const SpinDensity dens{c.first, c.second, theta, phi};
    const SpinDensity::LabMoments mom = dens.lab_moments();
    const AxisMoments axis = moments_along_axis(c.first, c.second);
    const double n_hat[3] = {std::sin(theta) * std::cos(phi),
                             std::sin(theta) * std::sin(phi), std::cos(theta)};
    for (std::size_t i = 0; i < 3; ++i) {
      CHECK_CLOSE_AT(mom.vector[i], axis.vector * n_hat[i], 0.0, 1e-12);
    }
    // the rank-2 moment transforms with P2(cos theta) along z
    const double p2 = 0.5 * (3 * std::cos(theta) * std::cos(theta) - 1);
    CHECK_CLOSE_AT(mom.tensor_zz, axis.tensor * p2, 0.0, 1e-12);
  }
}

TEST_CASE("unphysical populations are rejected") {
  CHECK_THROWS(rho_from_populations(1.0, {0.7, 0.6, -0.3}));
  CHECK_THROWS(rho_from_populations(1.0, {0.5, 0.5}));  // wrong length
  CHECK_THROWS(spin32_populations(1.2, 0.0, 0.0));
}

TEST_CASE("max-entropy populations obey the spin-temperature relation") {
  // p_m ~ exp(beta m) with <J_z>/J = pz.  For spin 1 the implied tensor
  // moment is P_zz = 2 - sqrt(4 - 3 pz^2), which is what
  // bookkeeping.helicity_flip_plan hands the money scripts as `pzz_true`.
  for (double pz : {0.2, 0.5, 0.7, 0.9}) {
    const std::vector<double> pops = populations_maxent(1.0, pz);
    double sum = 0.0;
    for (double p : pops) {
      CHECK(p >= 0.0);
      sum += p;
    }
    CHECK_CLOSE_AT(sum, 1.0, 0.0, 1e-12);
    const AxisMoments m = moments_along_axis(1.0, pops);
    CHECK_CLOSE_AT(m.vector, pz, 0.0, 1e-10);
    CHECK_CLOSE_AT(m.tensor, 2.0 - std::sqrt(4.0 - 3.0 * pz * pz), 0.0, 1e-12);
  }
  // spin temperature is a one-parameter family: the ratio is geometric
  const std::vector<double> p = populations_maxent(1.0, 0.6);
  CHECK_CLOSE(p[0] * p[2], p[1] * p[1], 1e-10);
  CHECK_THROWS(populations_maxent(1.0, 1.0));
}

TEST_CASE("max-entropy anchor (P_z, P_zz) = (8/13, 4/13) for J = 1") {
  // The rational t = 3 anchor of evgen/tests/test_bookkeeping.py: populations
  // (9, 3, 1)/13, whose moments are exactly 8/13 and 4/13.
  const std::vector<double> pops = populations_maxent(1.0, 8.0 / 13.0);
  const double want[3] = {9.0 / 13.0, 3.0 / 13.0, 1.0 / 13.0};
  for (std::size_t i = 0; i < 3; ++i) CHECK_CLOSE(pops[i], want[i], 1e-11);
  const AxisMoments m = moments_along_axis(1.0, pops);
  CHECK_CLOSE(m.vector, 8.0 / 13.0, 1e-11);
  CHECK_CLOSE(m.tensor, 4.0 / 13.0, 1e-11);
}

TEST_CASE("max-entropy populations for spin 3/2, anchor (0.7, 0.4)") {
  for (double pz : {0.3, 0.8}) {
    const std::vector<double> pops = populations_maxent(1.5, pz);
    double sum = 0.0;
    for (double p : pops) {
      CHECK(p >= 0.0);
      sum += p;
    }
    CHECK_CLOSE_AT(sum, 1.0, 0.0, 1e-12);
    CHECK_CLOSE_AT(moments_along_axis(1.5, pops).vector, pz, 0.0, 1e-10);
    // geometric, so consecutive ratios are equal
    const double r0 = pops[1] / pops[0], r1 = pops[2] / pops[1],
                 r2 = pops[3] / pops[2];
    CHECK_CLOSE(r0, r1, 1e-9);
    CHECK_CLOSE(r1, r2, 1e-9);
  }
  // the J = 3/2 anchor: pz = 0.7 gives tensor 0.4 and octupole 0.2 exactly,
  // from the geometric ratio 1/3 (populations 27, 9, 3, 1 over 40).
  const std::vector<double> pops = populations_maxent(1.5, 0.7);
  const double want[4] = {0.675, 0.225, 0.075, 0.025};
  for (std::size_t i = 0; i < 4; ++i) CHECK_CLOSE(pops[i], want[i], 1e-11);
  const AxisMoments m = moments_along_axis(1.5, pops);
  CHECK_CLOSE(m.vector, 0.7, 1e-11);
  CHECK_CLOSE(m.tensor, 0.4, 1e-11);
  CHECK_CLOSE(m.octupole, 0.2, 1e-10);
}

TEST_CASE("the density matrix of a pure state along a tilted axis") {
  // rho = |n><n| for a stretched state: rank one, unit trace, and the vector
  // moment is exactly n_hat.
  const double theta = kPi2 / 3.0, phi = 0.9;
  const CplxMatrix rho = rho_from_populations(1.0, {1.0, 0.0, 0.0}, theta, phi);
  CHECK_CLOSE_AT(trace(rho).real(), 1.0, 0.0, 1e-12);
  const std::array<double, 3> v = vector_polarization(rho, 1.0);
  CHECK_CLOSE_AT(v[0], std::sin(theta) * std::cos(phi), 0.0, 1e-12);
  CHECK_CLOSE_AT(v[1], std::sin(theta) * std::sin(phi), 0.0, 1e-12);
  CHECK_CLOSE_AT(v[2], std::cos(theta), 0.0, 1e-12);
}

TEST_CASE("D1/F3: the J = 3/2 population domain is stated exactly, and it is "
          "the function's own") {
  // The four populations are fixed UNIQUELY by (1, pz, t, o) -- the system is
  // square -- so an unphysical request is a DOMAIN statement and not a solver
  // failure.  Inverting by hand:
  //
  //   p(+3/2) + p(-3/2) = (1 + t)/2 ,  p(+3/2) - p(-3/2) = 0.9 pz + 0.1 o
  //   p(+1/2) + p(-1/2) = (1 - t)/2 ,  p(+1/2) - p(-1/2) = 0.3 (pz - o)
  //
  // so positivity is exactly the two inequalities below, and at o = 0
  // 1.8|pz| - 1 <= t <= 1 - 0.6|pz| (hence |pz| <= 5/6).  The refusal message
  // prints those bounds; this checks they ARE the boundary.
  auto reachable = [](double pz, double t, double o) {
    return std::fabs(0.9 * pz + 0.1 * o) <= 0.5 * (1.0 + t) + 1e-12
           && std::fabs(0.3 * (pz - o)) <= 0.5 * (1.0 - t) + 1e-12;
  };
  int checked = 0;
  for (int i = 0; i <= 20; ++i) {
    for (int j = 0; j <= 20; ++j) {
      for (int m = 0; m <= 4; ++m) {
        const double pz = -1.0 + 0.1 * i;
        const double t = -1.0 + 0.1 * j;
        const double o = -1.0 + 0.5 * m;
        bool ok = true;
        try {
          spin32_populations(pz, t, o);
        } catch (const std::runtime_error&) {
          ok = false;
        }
        CHECK(ok == reachable(pz, t, o));
        ++checked;
      }
    }
  }
  CHECK(checked == 21 * 21 * 5);
  // The inversion itself, on a point well inside: the four moments come back.
  const std::array<double, 4> p = spin32_populations(0.4, 0.2, 0.1);
  double s = 0.0, v = 0.0, tt = 0.0, oo = 0.0;
  const double ms[4] = {1.5, 0.5, -0.5, -1.5};
  for (int i = 0; i < 4; ++i) {
    s += p[i];
    v += p[i] * ms[i] / 1.5;
    tt += p[i] * (3.0 * ms[i] * ms[i] - 3.75) / 3.0;
    oo += p[i] * (ms[i] * ms[i] * ms[i] - (41.0 / 20.0) * ms[i]) / 0.3;
  }
  CHECK_CLOSE_AT(s, 1.0, 0.0, 1e-12);
  CHECK_CLOSE_AT(v, 0.4, 0.0, 1e-12);
  CHECK_CLOSE_AT(tt, 0.2, 0.0, 1e-12);
  CHECK_CLOSE_AT(oo, 0.1, 0.0, 1e-12);
  // ... and the message names the pair, the offending m and both edges.
  try {
    spin32_populations(0.7, 0.6, 0.0);
    CHECK(false);
  } catch (const std::runtime_error& e) {
    const std::string m(e.what());
    CHECK(m.find("unphysical (pz = 0.7, t = 0.6, o = 0)") != std::string::npos);
    CHECK(m.find("0.26 <= t <= 0.58") != std::string::npos);
    CHECK(m.find("|pz| <= 5/6") != std::string::npos);
  }
}

TEST_CASE("the spin-1 domain is a DOMAIN, and the refusal names its edges") {
  // The spin-1 counterpart of the J = 3/2 case above, and it exists because
  // `--pzz-mode typed` (2026-09-06) made this branch reachable from the
  // command line: `spin1_populations` used to throw "unphysical (pz, pzz):
  // negative population", which named neither the offending m nor the P_zz
  // this P_z admits.
  //
  // The three populations are fixed UNIQUELY by (1, pz, pzz), so positivity
  // is exactly 3|pz| - 2 <= pzz <= 1.  Checked against the function itself
  // on a 41 x 41 grid before it is quoted in a message.
  int disagreements = 0;
  for (int i = 0; i <= 40; ++i) {
    const double pz = -1.0 + 0.05 * i;
    for (int k = 0; k <= 40; ++k) {
      const double pzz = -2.0 + 0.075 * k;
      const bool inside = (pzz >= 3.0 * std::fabs(pz) - 2.0 - 1e-12) &&
                          (pzz <= 1.0 + 1e-12);
      bool threw = false;
      try {
        spin1_populations(pz, pzz);
      } catch (const std::runtime_error&) {
        threw = true;
      }
      if (threw == inside) ++disagreements;
    }
  }
  CHECK(disagreements == 0);

  // The edge is reached, not approached: at pzz = 3|pz| - 2 exactly, p(-1)
  // is 0 and the fill is still built.
  const std::array<double, 3> edge = spin1_populations(0.7, 0.1);
  CHECK_CLOSE_AT(edge[2], 0.0, 0.0, 1e-15);
  CHECK_CLOSE_AT(edge[0] + edge[1] + edge[2], 1.0, 0.0, 1e-15);

  // ... and the message names the offending m, its population and the edge.
  try {
    spin1_populations(0.7, -0.5);
    CHECK(false);
  } catch (const std::runtime_error& e) {
    const std::string m(e.what());
    CHECK(m.find("unphysical (pz = 0.7, pzz = -0.5)") != std::string::npos);
    CHECK(m.find("p(m = -1) = -0.1 < 0") != std::string::npos);
    CHECK(m.find("0.1 <= pzz <= 1") != std::string::npos);
    // and it says the J = 3/2 domain at the same pz is SMALLER, which is the
    // trap `--pzz-mode typed` walks into at the CLI's own default fill
    CHECK(m.find("0.26 <= T <= 0.58") != std::string::npos);
  }
}
