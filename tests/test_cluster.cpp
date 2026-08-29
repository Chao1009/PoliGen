// Cluster relative partial waves: the radial shapes and the azimuth-stripped
// spherical harmonics, against validation/reference/tagged.json (`waves[i]`)
// and against the analytic properties they must have.

#include <cmath>
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
