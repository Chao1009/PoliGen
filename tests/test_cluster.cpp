// Cluster relative partial waves: the radial shapes and the azimuth-stripped
// spherical harmonics, against validation/reference/tagged.json (`waves[i]`)
// and against the analytic properties they must have.

#include <cmath>
#include <fstream>
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

}  // namespace

TEST_CASE("cluster: Wave::radial reproduces the polligen radial tables") {
  jsonmin::Value ref;
  if (!load_tagged(ref)) {
    MESSAGE("tagged.json not found -- skipping");
    return;
  }
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

TEST_CASE("cluster: the ANL overlap reader") {
  if (!have(kLi6Overlap)) {
    MESSAGE("data/vmc not present -- skipping");
    return;
  }
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

TEST_CASE("cluster: the ANL momentum reader and its printed normalizations") {
  if (!have(kLi6Momentum)) {
    MESSAGE("data/vmc not present -- skipping");
    return;
  }
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

TEST_CASE("cluster: the r-space Fourier-Bessel route reproduces the k block") {
  if (!have(kLi6Overlap)) {
    MESSAGE("data/vmc not present -- skipping");
    return;
  }
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

TEST_CASE("cluster: the VMC S-D relative sign, and where it comes from") {
  if (!have(kLi6Overlap) || !have(kLi6Momentum)) {
    MESSAGE("data/vmc not present -- skipping");
    return;
  }
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
