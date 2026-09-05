// SPDX-License-Identifier: GPL-3.0-or-later
// Cluster relative partial waves: the radial shapes and the azimuth-stripped
// spherical harmonics, against validation/reference/tagged.json (`waves[i]`)
// and against the analytic properties they must have.

#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <stdexcept>
#include <memory>
#include <string>
#include <vector>

#include "check_close.hpp"
#include "doctest.h"
#include "json_min.hpp"
#include "lipolgen/cluster.hpp"
#include "lipolgen/constants.hpp"
#include "lipolgen/numerics.hpp"
#include "lipolgen/spectator.hpp"

using namespace lipolgen;

#ifndef LIPOLGEN_REFERENCE_DIR
#define LIPOLGEN_REFERENCE_DIR "validation/reference"
#endif

namespace {

constexpr double kRtol = 1e-12;

bool load_tagged(jsonmin::Value& out) {
  return jsonmin::load_file(std::string(LIPOLGEN_REFERENCE_DIR) + "/tagged.json",
                            out);
}

/// Is the reference blob on disk?  Called at doctest REGISTRATION time by the
/// `doctest::skip()` decorators below, so that a checkout without
/// `validation/reference/` reports those cases as SKIPPED in the tally
/// instead of PASSED with zero assertions (the T1v pattern of
/// tests/test_b1_nuclear.cpp; phase F, 2026-09-05).  A blob that IS there and
/// does not parse is a failure, not a skip -- that is what the `REQUIRE` on
/// the loader inside each case is for.
bool tagged_present() {
  return std::ifstream(std::string(LIPOLGEN_REFERENCE_DIR) +
                       "/tagged.json")
      .good();
}

}  // namespace

TEST_CASE("cluster: Wave::radial reproduces the polligen radial tables" *
          doctest::skip(!tagged_present())) {
  jsonmin::Value ref;
  REQUIRE(load_tagged(ref));
  for (const auto& kv : ref["channels"].obj()) {
    CAPTURE(kv.first);
    for (const jsonmin::Value& w : kv.second["waves"].arr()) {
      const Wave wave{static_cast<int>(w["l"].num()), w["prob"].num(),
                      w["beta"].num()};
      const double kappa = w["kappa"].num();
      const std::vector<double> kg = w["k_grid"].flat();
      const std::vector<double> rad = w["radial"].flat();
      CAPTURE(wave.l);
      for (std::size_t i = 0; i < kg.size(); ++i) {
        CHECK_CLOSE(wave.radial(kg[i], kappa), rad[i], kRtol);
      }
    }
  }
}

TEST_CASE("cluster: the S wave IS spectator.momentum_density's psi") {
  // The L = 0 Hulthen radial and the fast simulation's S-wave amplitude are
  // the same function -- which is what makes the P_D = 0 limit of the tagged
  // sampler agree with the fast simulation's spectator sampler.
  const ClusterChannel& ch = LI6_ALPHA_TAG();
  const double kappa = ch.kappa();
  const Wave s{0, 1.0, BETA_DEFAULT};
  for (double k = 1e-4; k < 1.2; k += 0.037) {
    const double psi = s.radial(k, kappa);
    CHECK_CLOSE(psi * psi, momentum_density(k, kappa, BETA_DEFAULT, 0), kRtol);
  }
  const Wave p{1, 1.0, BETA_DEFAULT};
  const double kap7 = LI7_ALPHA_TAG().kappa();
  for (double k = 1e-4; k < 1.2; k += 0.037) {
    const double psi = p.radial(k, kap7);
    CHECK_CLOSE(psi * psi, momentum_density(k, kap7, BETA_DEFAULT, 1), kRtol);
  }
}

TEST_CASE("cluster: the P and D waves vanish at the origin") {
  const Wave p{1, 1.0, BETA_DEFAULT};
  const Wave d{2, 1.0, BETA_DEFAULT};
  const Wave s{0, 1.0, BETA_DEFAULT};
  const double kappa = 0.06;
  CHECK(p.radial(0.0, kappa) == 0.0);
  CHECK(d.radial(0.0, kappa) == 0.0);
  CHECK(s.radial(0.0, kappa) > 0.0);
  // and are strictly suppressed against the S wave at k -> 0
  CHECK(p.radial(1e-4, kappa) < 2e-3 * s.radial(1e-4, kappa));
  CHECK(d.radial(1e-4, kappa) < 2e-6 * s.radial(1e-4, kappa));
  CHECK_THROWS(Wave{3, 1.0, BETA_DEFAULT}.radial(0.1, kappa));
}

TEST_CASE("cluster: theta_lm is the Condon-Shortley Y_l^m with the azimuth off") {
  // |Y_l^m|^2 integrated over the sphere is 1, i.e.
  //   2 pi integral_{-1}^{1} Theta_lm(c)^2 dc = 1.
  for (int l = 0; l <= 2; ++l) {
    for (int m = -l; m <= l; ++m) {
      const std::size_t n = 20001;
      const std::vector<double> c = linspace(-1.0, 1.0, n);
      std::vector<double> f(n);
      for (std::size_t i = 0; i < n; ++i) {
        const double t = theta_lm(l, m, c[i]);
        f[i] = t * t;
      }
      CAPTURE(l);
      CAPTURE(m);
      CHECK_CLOSE(2.0 * kPi * trapezoid(f, c), 1.0, 1e-7);
    }
  }
  // orthogonality at equal m, different l
  const std::size_t n = 20001;
  const std::vector<double> c = linspace(-1.0, 1.0, n);
  std::vector<double> f(n);
  for (std::size_t i = 0; i < n; ++i) f[i] = theta_lm(0, 0, c[i]) * theta_lm(2, 0, c[i]);
  CHECK_CLOSE_AT(2.0 * kPi * trapezoid(f, c), 0.0, 0.0, 1e-7);
  // the explicit closed forms
  CHECK_CLOSE(theta_lm(0, 0, 0.3), std::sqrt(1.0 / (4.0 * kPi)), kRtol);
  CHECK_CLOSE(theta_lm(1, 0, 0.3), std::sqrt(3.0 / (4.0 * kPi)) * 0.3, kRtol);
  CHECK_CLOSE(theta_lm(1, 1, 0.0), -std::sqrt(3.0 / (8.0 * kPi)), kRtol);
  CHECK_CLOSE(theta_lm(1, -1, 0.0), +std::sqrt(3.0 / (8.0 * kPi)), kRtol);
  CHECK_CLOSE(theta_lm(2, 0, 1.0), std::sqrt(5.0 / (16.0 * kPi)) * 2.0, kRtol);
  CHECK_CLOSE(theta_lm(2, 2, 0.0), std::sqrt(15.0 / (32.0 * kPi)), kRtol);
  CHECK_CLOSE(theta_lm(2, -2, 0.0), std::sqrt(15.0 / (32.0 * kPi)), kRtol);
  CHECK_THROWS(theta_lm(3, 0, 0.0));
  CHECK_THROWS(theta_lm(2, 3, 0.0));
}

// ------------------------------------------------------- the VMC backend
//
// The ANL VMC tables (data/vmc/, provenance in data/vmc/README.md) reached
// through `VmcRadial`.  Every reference number here comes from
// validation/vmc_reconcile.py, which recomputes them independently in Python
// with one normalization; docs/open_items/vmc_reconciliation.md carries the
// table.

namespace {

const std::string kLi6Overlap = data_path("vmc/li6_alpha_d/li6.ad");
const std::string kLi7Overlap = data_path("vmc/li7_alpha_t/li7.at");
const std::string kLi6Momentum = data_path("vmc/momenta/li6_ad1.momentum");
const std::string kLi7Momentum = data_path("vmc/momenta/li7_at3.momentum");

bool have(const std::string& path) {
  std::ifstream f(path);
  return static_cast<bool>(f);
}

/// hbar c, the one conversion the VMC tables need.
constexpr double kHbarC = 0.1973269804;

}  // namespace

TEST_CASE("cluster: data_dir honours LIPOLGEN_DATA_DIR, else the built-in") {
  // The compiled-in default is ${CMAKE_SOURCE_DIR}/data; nothing here may
  // depend on the CWD.
  CHECK(!data_dir().empty());
  CHECK(data_path("x/y") == data_dir() + "/x/y");
}

TEST_CASE("cluster: the ANL overlap reader" *
          doctest::skip(!have(kLi6Overlap))) {
  REQUIRE(have(kLi6Overlap));   // the decorator already guaranteed it
  const std::vector<AnlTable> b = read_anl_overlap(kLi6Overlap);
  REQUIRE(b.size() == 2);
  // k block: 0.00 .. 5.00 fm^-1 in 0.10 steps, two signed columns + errors
  CHECK(b[0].x.size() == 51);
  CHECK_CLOSE_AT(b[0].x.front(), 0.0, 0.0, 1e-12);
  CHECK_CLOSE_AT(b[0].x.back(), 5.0, 0.0, 1e-12);
  CHECK_CLOSE_AT(b[0].col[0][0], -123.6, 0.0, 1e-9);      // Aad00(0)
  CHECK_CLOSE_AT(b[0].err[0][0], 5.129, 0.0, 1e-9);
  CHECK_CLOSE_AT(b[0].col[1][5], 3.987, 0.0, 1e-9);       // Aad22(0.5)
  // r block: 0.05 fm up, and the split_numbers() parser must survive the
  // `(.7251E-01)-0.1357E-01` run-together column the Fortran format writes
  CHECK_CLOSE_AT(b[1].x.front(), 0.05, 0.0, 1e-12);
  CHECK_CLOSE_AT(b[1].col[0][0], 0.5846, 0.0, 1e-9);
  CHECK_CLOSE_AT(b[1].col[1][0], -0.1357e-01, 0.0, 1e-12);
}

TEST_CASE("cluster: the ANL momentum reader and its printed normalizations" *
          doctest::skip(!have(kLi6Momentum))) {
  REQUIRE(have(kLi6Momentum));   // the decorator already guaranteed it
  const std::vector<AnlTable> b = read_anl_momentum(kLi6Momentum);
  REQUIRE(b.size() == 2);
  CHECK(b[0].col.size() == 1);            // total rho(K)
  CHECK(b[1].col.size() == 2);            // rho_0, rho_2
  CHECK_CLOSE_AT(b[0].col[0][0], 1041.5, 0.0, 1e-9);
  CHECK_CLOSE_AT(b[1].col[1][5], 1.0509, 0.0, 1e-9);
  // the S/D blocks add up to the total block, bin by bin -- to the 5
  // significant digits the file prints, which is all it carries
  for (std::size_t i = 1; i < b[1].x.size(); ++i) {
    CHECK_CLOSE_AT(b[1].col[0][i] + b[1].col[1][i], b[0].col[0][i], 1e-4, 1e-6);
  }
  const std::vector<double> n = read_anl_momentum_norms(kLi6Momentum);
  REQUIRE(n.size() == 3);
  CHECK_CLOSE_AT(n[0], 0.81971, 0.0, 1e-9);
  CHECK_CLOSE_AT(n[1], 0.80362, 0.0, 1e-9);
  CHECK_CLOSE_AT(n[2], 0.015861, 0.0, 1e-9);
  const std::vector<double> n7 = read_anl_momentum_norms(kLi7Momentum);
  REQUIRE(n7.size() == 1);
  CHECK_CLOSE_AT(n7[0], 1.0084, 0.0, 1e-9);
}

TEST_CASE("cluster: VmcRadial interpolates and is ZERO outside its table") {
  const VmcRadial v({0.0, 1.0, 2.0}, {0.0, 10.0, 20.0}, 0, "unit test");
  CHECK_CLOSE_AT(v(0.5), 5.0, 0.0, 1e-12);
  CHECK_CLOSE_AT(v(2.0), 20.0, 0.0, 1e-12);
  CHECK(v(2.0000001) == 0.0);
  CHECK(v(-0.1) == 0.0);
  CHECK(v.l() == 0);
  CHECK_CLOSE_AT(v.scaled(-2.0)(0.5), -10.0, 0.0, 1e-12);
  CHECK_THROWS(VmcRadial({1.0, 0.0}, {1.0, 2.0}, 0));       // not increasing
  CHECK_THROWS(VmcRadial({0.0}, {1.0}, 0));                 // one point
  // ... and a Wave with a table ignores kappa and beta entirely
  Wave w;
  w.l = 0;
  w.vmc = std::make_shared<const VmcRadial>(v);
  CHECK_CLOSE_AT(w.radial(0.5, 0.06), 5.0, 0.0, 1e-12);
  CHECK_CLOSE_AT(w.radial(0.5, 999.0), 5.0, 0.0, 1e-12);
}

TEST_CASE("cluster: the r-space Fourier-Bessel route reproduces the k block" *
          doctest::skip(!have(kLi6Overlap))) {
  REQUIRE(have(kLi6Overlap));   // the decorator already guaranteed it
  // The two blocks of li6.ad are INDEPENDENT Monte Carlo estimators of the
  // same overlap, one in r and one in k.  Transforming the r block with the
  // ANL cluster convention A(k) = 4pi int A(r) j_L(kr) r^2 dr must land on
  // the k block -- that is what validates the convention (and with it the
  // r-space-only tables, li6n_31.table / h2n.table, for future use).
  for (int col = 0; col < 2; ++col) {
    const int l = col == 0 ? 0 : 2;
    const VmcRadial ft = vmc_from_overlap_r(kLi6Overlap, col, l);
    const VmcRadial nat = vmc_from_overlap_k(kLi6Overlap, col, l);
    for (double k_fm = 0.2; k_fm <= 2.0; k_fm += 0.2) {
      const double k = k_fm * kHbarC;
      CHECK_CLOSE_AT(ft(k), nat(k), 5e-3, 5e-3);
    }
  }
}

TEST_CASE("cluster: the VMC S-D relative sign, and where it comes from" *
          doctest::skip(!have(kLi6Overlap) || !have(kLi6Momentum))) {
  // the decorator already guaranteed both
  REQUIRE(have(kLi6Overlap));
  REQUIRE(have(kLi6Momentum));
  // A momentum DENSITY carries no phase; the overlap amplitudes do.  Both
  // blocks of li6.ad, and the Fourier-Bessel transform of the r block, agree
  // that A22 and A00 have OPPOSITE sign wherever the D wave has strength.
  const VmcRadial a00 = vmc_from_overlap_k(kLi6Overlap, 0, 0);
  const VmcRadial a22 = vmc_from_overlap_k(kLi6Overlap, 1, 2);
  for (double k_fm : {0.1, 0.2, 0.4, 0.6}) {
    CHECK(a00(k_fm * kHbarC) * a22(k_fm * kHbarC) < 0.0);
  }
  // exactly ONE node each below 3 fm^-1, at the published positions
  const std::vector<double> k_fm00 = [&] {
    std::vector<double> x(a00.k().size());
    for (std::size_t i = 0; i < x.size(); ++i) x[i] = a00.k()[i] / kHbarC;
    return x;
  }();
  const std::vector<double> n00 = vmc_sign_steps(k_fm00, a00.psi(), 3.0);
  REQUIRE(n00.size() == 1);
  CHECK_CLOSE_AT(n00[0], 0.678, 0.0, 0.005);
  const std::vector<double> n22 = vmc_sign_steps(k_fm00, a22.psi(), 3.0);
  REQUIRE(n22.size() == 1);
  CHECK_CLOSE_AT(n22[0], 2.250, 0.0, 0.005);

  // The momentum-file tables inherit that structure through the anchor.
  const VmcRadial s = vmc_from_momentum(kLi6Momentum, 1, 0, 0, &a00);
  const VmcRadial d = vmc_from_momentum(kLi6Momentum, 1, 1, 2, &a22);
  for (double k_fm : {0.1, 0.2, 0.4, 0.6}) {
    CHECK(s(k_fm * kHbarC) * d(k_fm * kHbarC) < 0.0);
  }
  // and the magnitude is sqrt(rho), i.e. |psi|^2 = rho of the file
  CHECK_CLOSE_AT(s(0.5 * kHbarC) * s(0.5 * kHbarC), 43.685, 0.0, 1e-3);
  CHECK_CLOSE_AT(d(0.5 * kHbarC) * d(0.5 * kHbarC), 1.0509, 0.0, 1e-4);
  // no sign reference => positive definite (correct only for one wave)
  const VmcRadial plain = vmc_from_momentum(kLi6Momentum, 1, 0, 0, nullptr);
  for (double v : plain.psi()) CHECK(v >= 0.0);
}

// ------------------------------------------------ the CD-Bonn deuteron (A5)
//
// THE COEFFICIENT GATE.  `cluster.hpp` types twenty published numbers
// (Machleidt, PRC 63 (2001) 024001 = arXiv:nucl-th/0006014, Table XX, LaTeX
// label `tab_dwpar`) and computes four more from the r -> 0 boundary
// conditions.  Nothing else in the library would notice a mistyped digit: the
// wave function would simply be a slightly different one, the A = 2 gate
// would move, and the move would be read as physics.  So this test asks the
// coefficients for CD-Bonn's OWN published deuteron properties.
//
// It needs NO data file -- the parameterisation is analytic and the momentum
// moments are CLOSED FORM.  With
//   (2/pi) int_0^inf dp p^2/[(p^2+a^2)(p^2+b^2)] = 1/(a+b)
// the norms are the double sums sum_ij C_i C_j/(m_i+m_j) and the same in D,
// so what is tested here is the COEFFICIENTS, with no quadrature error
// anywhere to hide behind.
//
// HOW SHARP IS IT, measured rather than asserted.  Perturbing each of the 18
// published coefficients by one unit in its last printed digit and
// recomputing moves the normalisation by, in units of 1e-9 relative:
//   C1 25   C2 4.6  C3 0.21  C4 79   C5 42   C6 243  C7 143  C8 83
//   C9 45   C10 18  D1 16    D2 4.8  D3 0.96 D4 23   D5 70   D6 21
//   D7 55   D8 9.8
// and the D-state probability by 20 to 1400e-9 relative for the D_j.  The
// pins below sit at 1e-9, so a single wrong digit ANYWHERE in the table is
// caught, with ONE honest exception: the last digit of C_3 = -0.44114404e-01,
// whose eighth significant figure is an absolute 1e-9 on the smallest
// coefficient in the table, is worth 2.1e-10 and is NOT resolvable by any
// double-precision observable.  A slip in its seventh figure IS caught (2.1e-9).
TEST_CASE("cluster: the CD-Bonn deuteron coefficients (Machleidt Table XX)") {
  const CdBonnWave& w = cdbonn_wave();

  // ---- the four boundary conditions, Eqs. (D23)/(D24) in the sum-rule form.
  //      Machleidt: they "must be enforced by double precision (i.e., to
  //      about 15 decimal digits), otherwise the wave function is not
  //      reproduced correctly for r <= 0.5 fm".  The tolerances below are
  //      that demand made mechanical: each residual is compared to the SCALE
  //      of the sum it cancels (sum |term|), so they are relative statements
  //      and do not silently loosen if a coefficient changes.
  const std::array<double, 4> res = w.constraint_residuals();
  double s_c = 0.0, s_dm = 0.0, s_d = 0.0, s_dp = 0.0;
  for (std::size_t j = 0; j < w.m.size(); ++j) {
    const double x = w.m[j] * w.m[j];
    s_c += std::fabs(w.c[j]);
    s_dm += std::fabs(w.d[j]) / x;
    s_d += std::fabs(w.d[j]);
    s_dp += std::fabs(w.d[j]) * x;
  }
  MESSAGE("CD-Bonn constraint residuals: sum C = " << res[0] << ", sum D/m^2 = "
          << res[1] << ", sum D = " << res[2] << ", sum D m^2 = " << res[3]);
  CHECK(res[0] == 0.0);                       // C_11 = -sum_{j<11} C_j, exactly
  CHECK(std::fabs(res[1]) < 1e-14 * s_dm);
  CHECK(std::fabs(res[2]) < 1e-14 * s_d);
  CHECK(std::fabs(res[3]) < 1e-14 * s_dp);

  // ---- NORMALISATION, Eq. (D14): (2/pi) int dp p^2 (u^2 + w^2) = 1.
  //      TWO statements, and they are different tests.
  //      (a) the physics: the parameterisation is a FIT to Machleidt's
  //          numerical wave function (his quoted L2 quality is 2.2e-4 in u and
  //          1.1e-4 in w), so it reproduces unity to 1.7e-7, not to machine
  //          precision, and demanding better would be demanding the fit be
  //          something it is not;
  //      (b) the pin: 1e-9 relative on the computed value, which is what
  //          actually catches a mistyped digit (see the table above).  The
  //          value is route-dependent at the 5e-11 level -- the double sum
  //          cancels 8e+4 down to 1, so it carries ~1e-11 of rounding however
  //          it is summed -- and 1e-9 sits 20x above that.
  MESSAGE("CD-Bonn norm = " << w.norm() << " (S " << w.norm_s() << ", D "
                            << w.norm_d() << ")");
  CHECK(std::fabs(w.norm() - 1.0) < 5e-7);
  CHECK_CLOSE(w.norm(), 0.99999982615384, 1e-9);
  CHECK_CLOSE(w.norm_s(), 0.95143775250961, 1e-9);

  // ---- P_D.  CD-Bonn's Table XV publishes 4.85 %; the parameterisation
  //      gives 4.8562 %, which does NOT round to 4.85 -- it is a fit, and
  //      6.2e-5 absolute is the size of the fit residual, not of an error.
  //      That is stated rather than hidden by a loose tolerance: the
  //      published-value clause allows 1.5e-4 and the PIN is at 1e-9.
  MESSAGE("CD-Bonn P_D = " << 100.0 * w.norm_d() << " % (published 4.85 %)");
  CHECK(std::fabs(w.norm_d() - 0.0485) < 1.5e-4);
  CHECK_CLOSE(w.norm_d(), 0.048562073644234, 1e-9);

  // ---- the asymptotics, Eq. (D15): A_S = C_1 and A_D = D_1 exactly, so eta
  //      is a RATIO OF TWO PUBLISHED NUMBERS and needs no integral.  Table XV
  //      prints A_S = 0.8846(9) fm^-1/2 and eta = 0.0256(4) for CD-Bonn
  //      (the parentheses are the EMPIRICAL errors it is compared against).
  //      The parameterisation gives 0.884730 and 0.0255714, i.e. it misses
  //      Machleidt's own A_S by 1.3e-4 -- a seventh of the experimental
  //      uncertainty on that quantity, and the same kind of residual as P_D
  //      above: this is a FIT to his numerical wave function, not the wave
  //      function.  The tolerances say that and no more; the PINS are what
  //      catch a typo.
  MESSAGE("CD-Bonn A_S = " << w.a_s() << " (published 0.8846), eta = "
                           << w.eta() << " (published 0.0256)");
  CHECK(std::fabs(w.a_s() - 0.8846) < 2e-4);
  CHECK(std::fabs(w.eta() - 0.0256) < 5e-5);
  CHECK_CLOSE(w.a_s(), 0.88472985, 1e-15);
  CHECK_CLOSE(w.eta(), 0.0255713786530431, 1e-13);
  CHECK(w.a_d() == w.eta() * w.a_s());

  // ---- the SIGN, cheaply: w(p) = -psi_2^a must be POSITIVE at small p (the
  //      `fdeut.av18` / CDKS convention -- see the header, and
  //      `alpha_d_quadrupole_fm2` in tests/test_b1_nuclear.cpp for the gate
  //      that settles it against Q_d).  u(p) is positive there too, and both
  //      D-wave endpoints vanish: w(0) = 0 is the sum_j D_j/m_j^2 = 0
  //      constraint seen from momentum space.
  CHECK(w.psi_s(0.1) > 0.0);
  CHECK(w.psi_d(0.1) > 0.0);
  CHECK(std::fabs(w.psi_d(0.0)) < 1e-12 * w.psi_s(0.0));
  // The S wave has its node where the paper's own wave function does, and it
  // is 13 % above AV18's 2.0929 fm^-1 -- the single feature that drives most
  // of the difference the A = 2 gate sees (phase_A_cdbonn.md section 7).
  CHECK(w.psi_s(2.0) > 0.0);
  CHECK(w.psi_s(3.0) < 0.0);

  // ---- the `FdeutTable` adapter: `fdeut.av18`'s OWN grid by default, the
  //      same "both or neither" hbar c conversion `read_fdeut_k` applies, and
  //      CD-Bonn's own binding energy rather than AV18's.
  const FdeutTable t = cdbonn_fdeut_table();
  REQUIRE(t.k_gev.size() == 201);
  CHECK(t.k_gev.front() == 0.0);
  CHECK_CLOSE(t.k_gev.back(), 20.0 * HBARC_GEV_FM, 1e-15);
  CHECK_CLOSE(t.ebind_gev, CD_BONN_BINDING_MEV * 1e-3, 1e-15);
  const double scale = 1.0 / (HBARC_GEV_FM * std::sqrt(HBARC_GEV_FM));
  for (double p : {0.1, 0.7, 2.0, 5.0}) {
    const std::size_t i = static_cast<std::size_t>(std::llround(p / 0.1));
    CHECK_CLOSE(t.k_gev[i], p * HBARC_GEV_FM, 1e-14);
    CHECK_CLOSE(t.u[i], w.psi_s(p) * scale, 1e-15);
    CHECK_CLOSE(t.w[i], w.psi_d(p) * scale, 1e-15);
  }
  // "Both or neither": the GeV table must carry the SAME normalisation as the
  // fm one, i.e. hbar c^3 must not appear.  Trapezoid on the file's own 201
  // nodes, so this is 0.99998 rather than the closed form's 0.9999998 -- the
  // 2.2e-5 is the 0.1 fm^-1 spacing and the 20 fm^-1 truncation, exactly as
  // for `fdeut.av18` itself (0.9999764).
  double n = 0.0;
  for (std::size_t i = 1; i < t.k_gev.size(); ++i) {
    const double a = t.k_gev[i - 1] * t.k_gev[i - 1]
                     * (t.u[i - 1] * t.u[i - 1] + t.w[i - 1] * t.w[i - 1]);
    const double b = t.k_gev[i] * t.k_gev[i]
                     * (t.u[i] * t.u[i] + t.w[i] * t.w[i]);
    n += 0.5 * (a + b) * (t.k_gev[i] - t.k_gev[i - 1]);
  }
  MESSAGE("CD-Bonn on the fdeut grid: int k^2 (u^2 + w^2) dk = " << n);
  CHECK_CLOSE(n, 0.999977985, 1e-7);
  // A finer grid changes the TABLE and not the wave function.
  const FdeutTable f = cdbonn_fdeut_table(20.0, 0.02);
  CHECK(f.k_gev.size() == 1001);
  CHECK_CLOSE(f.u[50], t.u[10], 1e-15);   // both are p = 1.0 fm^-1
  CHECK_CLOSE(f.w[50], t.w[10], 1e-15);
  CHECK_THROWS_AS(cdbonn_fdeut_table(20.0, 0.0), std::runtime_error);
  CHECK_THROWS_AS(cdbonn_fdeut_table(0.01, 0.1), std::runtime_error);
}
