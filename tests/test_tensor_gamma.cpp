// SPDX-License-Identifier: GPL-3.0-or-later
// The exact finite-gamma tensor kernel (plans/08 D2, Cosyn et al.).
//
// `InclusiveKernel::Options::tensor_gamma` replaces the massless
// Hoodbhoy-Jaffe-Manohar b-sector with Cosyn, Roldan Tomei, Sosa and Zec,
// EPJ A 61 (2025) 83 (arXiv:2410.12764) Eqs. (9), (10), (14), (16), (17) and
// (24) -- the four tensor structure functions with explicit geometry, the
// b3/b4 slots, and the rest-frame angle theta_q that tilts the spin axis away
// from the photon direction and so leaks the b-sector into cos 2phi.
//
// Ported one for one from PolarizedLithiumSim/evgen/tests/test_tensor_gamma.py.
//
// Two kinds of test are kept apart here, because only one of them can catch a
// transcription error.  `table1_row2_along_q` and `table1_row3_along_the_beam`
// are anchored on a DERIVED result of the paper -- the two finite-gamma rows
// of its Table 1, closed-form functions of (gamma^2, eps, theta_q) that the
// module has never seen -- and it was those two rows, and not the re-typing
// below, that caught the two transcription errors this file's history
// records.  Everything else is a second ROUTE through the same transcription:
// `sfs_paper` writes Eqs. (17) out again and `harmonics_numeric` projects the
// cross-section weight numerically instead of expanding it in closed form,
// which tests the projection and the geometry but shares the module's reading
// of the equations.

#include <cmath>
#include <vector>

#include "check_close.hpp"
#include "lipolgen/asymmetries.hpp"
#include "lipolgen/beams.hpp"
#include "lipolgen/sf.hpp"
#include "lipolgen/xsec.hpp"

using namespace lipolgen;

namespace {

struct Point {
  double x, q2, y;
};

/// (x, Q2, y): a low-Q2 sweet spot of the low configuration, a mid one, and a
/// deliberately large-gamma^2 point at the edge of the W^2 >= 10 cut.
const Point kPoints[3] = {{0.08913, 1.135, 0.01561},
                          {0.1413, 3.127, 0.02713},
                          {0.2, 1.0, 0.05}};

double s_mid_6li() { return default_configs("6Li")[1].s_per_nucleon(); }

/// A hand-made SF table: F1 = 1, F2 = 2x(1+R)F1, b2 = 2x b1.
SFTables hand_tables(double x, double b1 = 0.02, double b3 = 0.0,
                     double b4 = 0.0, double r = 0.3) {
  SFTables t;
  t.f1 = 1.0;
  t.f2 = 2.0 * x * (1.0 + r);
  t.b1 = b1;
  t.b2 = 2.0 * x * b1;
  t.b3 = b3;
  t.b4 = b4;
  t.delta = 0.0;
  t.g1 = 0.0;
  t.g2 = 0.0;
  t.has_g2 = true;
  return t;
}

// --- a second, literal transcription of the paper -------------------------

struct PaperSFs {
  double f_t, f_l, f_lt, f_tt;
};

/// Cosyn Eqs. (17a)-(17e), written out again from the PDF, with the library's
/// per-nucleon map x_d -> x.  The leading 2 of Eq. (17a) multiplies b1 alone:
/// the large bracket opens before it.
PaperSFs sfs_paper(double b1, double b2, double b3, double b4, double x,
                   double g2) {
  const double g = std::sqrt(g2);
  const double a = 1.0 + g2;
  PaperSFs o;
  o.f_t = -(2.0 * a * b1 - (g2 / x) * ((1.0 / 6.0) * b2 - (1.0 / 2.0) * b3));
  o.f_l = (1.0 / x) * (2.0 * a * x * b1
                       - a * a * ((1.0 / 3.0) * b2 + b3 + b4)
                       - a * ((1.0 / 3.0) * b2 - b4)
                       - ((1.0 / 3.0) * b2 - b3));
  o.f_lt = -(g / (2.0 * x)) * (a * ((1.0 / 3.0) * b2 - b4)
                               + ((2.0 / 3.0) * b2 - 2.0 * b3));
  o.f_tt = -(g2 / x) * ((1.0 / 6.0) * b2 - (1.0 / 2.0) * b3);
  return o;
}

/// (h0, h1, h2) by NUMERICAL projection of the tensor weight over the lepton
/// azimuth -- the same physics as the module's closed form, by a different
/// route.  Eqs. (9), (10), (14), (16), (17), (19), (24).
TensorHarmonics harmonics_numeric(const SFTables& t, double x, double q2,
                                  double y, double theta_s, double m = 1.0,
                                  double j = 1.0, int n = 8192) {
  const double g2 = 4.0 * M_NUCLEON * M_NUCLEON * x * x / q2;
  const double eps = ((1.0 - y - 0.25 * g2 * y * y)
                      / (1.0 - y + 0.5 * y * y + 0.25 * g2 * y * y));
  const double root = std::sqrt(1.0 + g2);
  const double c = (1.0 + 0.5 * g2 * y) / root;
  const double s = std::sqrt(g2 * (1.0 - y - 0.25 * g2 * y * y)) / root;
  const double ct = std::cos(theta_s), st = std::sin(theta_s);
  const double q_nn = (3.0 * m * m - j * (j + 1.0)) / 3.0;
  const double pref = 1.5 * q_nn;  // Eq. (9) for a pure state
  const PaperSFs f = sfs_paper(t.b1, t.b2, t.b3, t.b4, x, g2);
  const double den = 2.0 * t.f1 + eps * ((1.0 + g2) * t.f2 / x - 2.0 * t.f1);
  double s0 = 0.0, s1 = 0.0, s2 = 0.0;
  for (int i = 0; i < n; ++i) {
    const double phi = 2.0 * kPi * static_cast<double>(i) / n;
    const double cp = std::cos(phi);
    const double nx = c * st * cp + s * ct;
    const double ny = st * std::sin(phi);
    const double nz = c * ct - s * st * cp;
    const double t_ll = pref * (nz * nz - 1.0 / 3.0);       // Eq. (14a)
    const double t_lt = pref * nx * nz;                     // Eq. (14b)
    const double t_tt = pref * (nx * nx - ny * ny);         // Eq. (14c)
    const double num = (t_ll * (f.f_t + eps * f.f_l)
                        + t_lt * std::sqrt(2.0 * eps * (1.0 + eps)) * f.f_lt
                        + t_tt * eps * f.f_tt);             // Eq. (19)
    const double w = -TENSOR_LL_SIGN * num / den;
    s0 += w;
    s1 += w * cp;
    s2 += w * std::cos(2.0 * phi);
  }
  TensorHarmonics h;
  h.h0 = s0 / n;
  h.h1 = 2.0 * s1 / n;
  h.h2 = 2.0 * s2 / n;
  return h;
}

InclusiveKernel gamma_kernel(SFFunc3 b1 = nullptr, SFFunc3 b3 = nullptr,
                             SFFunc3 b4 = nullptr) {
  InclusiveKernel::Options o;
  o.b1_func = std::move(b1);
  o.b3_func = std::move(b3);
  o.b4_func = std::move(b4);
  o.tensor_gamma = true;
  return InclusiveKernel(LI6(), o);
}

SFFunc3 toy_b1_func() {
  return [](double x, double q2, double f1) { return toy_b1(x, q2, f1); };
}

}  // namespace

// --- the massless limit ---------------------------------------------------

TEST_CASE("tensor_gamma: the Cosyn SFs reduce to the HJM b-sector at gamma = 0") {
  // At gamma = 0, -(F_TLL_T + eps F_TLL_L)/(F_UU_T + eps F_UU_L) is EXACTLY
  // the massless (b1 + (1-y)/(x y^2) b2)/D_phi of the master formula -- for
  // any b2, and with b3, b4 cancelling identically.
  const double kRtol = 1e-12;
  const double b2_scales[3] = {2.0, 1.3, 2.0};
  const double b3s[3] = {0.0, 0.0, 0.011};
  const double b4s[3] = {0.0, 0.0, -0.007};
  for (const Point& p : kPoints) {
    for (int i = 0; i < 3; ++i) {
      SFTables t = hand_tables(p.x, 0.02, b3s[i], b4s[i]);
      t.b2 = b2_scales[i] * p.x * t.b1;
      const double eps = epsilon_gamma(p.y, 0.0);
      const CosynTensorSFs f =
          cosyn_tensor_sfs(t.b1, t.b2, t.b3, t.b4, p.x, 0.0);
      const std::pair<double, double> fu =
          cosyn_unpolarized_sfs(t.f1, t.f2, p.x, 0.0);
      const double got =
          -(f.f_t + eps * f.f_l) / (fu.first + eps * fu.second);
      const double kern = t.b1 + (1.0 - p.y) / (p.x * p.y * p.y) * t.b2;
      const double dphi = t.f1 + (1.0 - p.y) / (p.x * p.y * p.y) * t.f2;
      CHECK_CLOSE(got, kern / dphi, kRtol);
      CHECK(f.f_tt == 0.0);
    }
  }
}

TEST_CASE("tensor_gamma: theta_q is a unit vector and vanishes masslessly") {
  for (const Point& p : kPoints) {
    const double g2 = gamma_squared(p.x, p.q2);
    const std::pair<double, double> cs = theta_q_cos_sin(p.y, g2);
    CHECK_CLOSE(cs.first * cs.first + cs.second * cs.second, 1.0, 1e-14);
    CHECK(cs.second > 0.0);
  }
  const std::pair<double, double> c0 = theta_q_cos_sin(0.3, 0.0);
  CHECK(c0.first == 1.0);
  CHECK(c0.second == 0.0);
}

TEST_CASE("tensor_gamma: the exact path reduces to the massless one at large Q2") {
  // gamma^2 -> 0 by taking Q^2 -> infinity at fixed x: the two kernels must
  // agree, and at a physical Q^2 they must differ by O(gamma^2).
  const double x = 0.1;
  const double q2s[2] = {1.0e8, 1.0e4};
  const double tols[2] = {1e-9, 5e-5};
  for (int i = 0; i < 2; ++i) {
    // s is chosen so that y stays at 1e-3 while gamma^2 -> 0; the kernels see
    // s only through y
    const double s_big = q2s[i] / (1.0e-3 * x);
    InclusiveKernel::Options massless;
    massless.b1_func = toy_b1_func();
    const InclusiveKernel k0(LI6(), massless);
    const InclusiveKernel k1 = gamma_kernel(toy_b1_func());
    const EventSpinState st{0, 0.0, 1.0, 1.0, 0.0, 0.0};
    const double w0 =
        k0.amplitudes(k0.tables(x, q2s[i]), x, q2s[i], s_big, st).w_avg;
    const double w1 =
        k1.amplitudes(k1.tables(x, q2s[i]), x, q2s[i], s_big, st).w_avg;
    CHECK_CLOSE(w1, w0, tols[i]);
    // and the residual really is O(gamma^2), coefficient of order 3
    CHECK(std::fabs(w1 / w0 - 1.0) < 5.0 * gamma_squared(x, q2s[i]));
  }
}

TEST_CASE("tensor_gamma: the massless path is untouched by the new slots") {
  // The default kernel is bit-for-bit the pre-D2 one: b3/b4 tables are built
  // but no massless amplitude reads them.
  const double s = s_mid_6li();
  const double x = 0.056, q2 = 1.14;
  const double y = q2 / (s * x);
  InclusiveKernel::Options base;
  base.b1_func = toy_b1_func();
  InclusiveKernel::Options with_b34 = base;
  with_b34.b3_func = [](double, double, double f1) { return 0.3 * f1; };
  with_b34.b4_func = [](double, double, double f1) { return -0.2 * f1; };
  const InclusiveKernel kerns[2] = {InclusiveKernel(LI6(), base),
                                    InclusiveKernel(LI6(), with_b34)};
  const EventSpinState st{0, 0.0, 1.0, 1.0, 0.4, 0.0};
  for (const InclusiveKernel& k : kerns) {
    CHECK(!k.tensor_gamma());
    const SFTables t = k.tables(x, q2);
    const Amplitudes got = k.amplitudes(t, x, q2, s, st);
    const double expected =
        TENSOR_LL_SIGN * (1.0 / 3.0) * 0.5 * (3.0 * std::cos(0.4) * std::cos(0.4) - 1.0)
        * (t.b1 + (1.0 - y) / (x * y * y) * t.b2)
        / (t.f1 + (1.0 - y) / (x * y * y) * t.f2);
    CHECK(got.w_avg == expected);  // bit for bit
    CHECK(got.a1 == 0.0);
  }
}

// --- the harmonics against an independent projection ----------------------

TEST_CASE("tensor_gamma: the harmonics match a direct numerical projection") {
  const InclusiveKernel kern = gamma_kernel();
  const double thetas[4] = {0.0, 0.5 * kPi, 0.7, 2.3};
  const double ms[3] = {1.0, 0.0, -1.0};
  for (const Point& p : kPoints) {
    const SFTables t = hand_tables(p.x, 0.02, 0.004, -0.002);
    for (double theta_s : thetas) {
      for (double m : ms) {
        const EventSpinState st{0, 0.0, 1.0, m, theta_s, 0.0};
        const TensorHarmonics got =
            kern.tensor_harmonics_gamma(t, p.x, p.q2, p.y, st);
        const TensorHarmonics want =
            harmonics_numeric(t, p.x, p.q2, p.y, theta_s, m);
        CHECK_CLOSE_AT(got.h0, want.h0, 1e-10, 1e-15);
        CHECK_CLOSE_AT(got.h1, want.h1, 1e-10, 1e-15);
        CHECK_CLOSE_AT(got.h2, want.h2, 1e-10, 1e-15);
      }
    }
  }
}

TEST_CASE("tensor_gamma: the transverse fill has no cos phi harmonic") {
  // At theta_S = 90 deg every cos phi' coefficient carries sin theta_S
  // cos theta_S and vanishes: the leakage is cos 2phi' only, which is why it
  // lands on the observable.
  const InclusiveKernel kern = gamma_kernel();
  const SFTables t = hand_tables(0.1, 0.02, 0.004, -0.002);
  const EventSpinState st{0, 0.0, 1.0, 1.0, 0.5 * kPi, 0.0};
  const TensorHarmonics h = kern.tensor_harmonics_gamma(t, 0.1, 1.14, 0.03, st);
  CHECK(std::fabs(h.h1) < 1e-18);
  CHECK(std::fabs(h.h2) > 0.0);
}

// --- pinned numbers -------------------------------------------------------

TEST_CASE("tensor_gamma: the pinned values at three points") {
  // (h0, h1, h2) at theta_S = 90 deg, m = +1, for hand_tables(x) with
  // b1 = 0.02, F1 = 1, R = 0.3, b3 = b4 = 0 -- the hand-derived combination of
  // Eqs. (17d) and (17e) with the T_LL leakage, at the three points.  The
  // same numbers `evgen/tests/test_tensor_gamma.py` pins.
  const double want0[3] = {2.5543780132391e-03, 2.5558792618530e-03,
                           2.5180085726965e-03};
  const double want2[3] = {-3.0869622799178e-05, -2.7677909400843e-05,
                           -1.6754330602946e-04};
  const InclusiveKernel kern = gamma_kernel();
  const EventSpinState st{0, 0.0, 1.0, 1.0, 0.5 * kPi, 0.0};
  for (int i = 0; i < 3; ++i) {
    const Point& p = kPoints[i];
    const TensorHarmonics h = kern.tensor_harmonics_gamma(
        hand_tables(p.x), p.x, p.q2, p.y, st);
    CHECK_CLOSE(h.h0, want0[i], 1e-11);
    CHECK_CLOSE(h.h2, want2[i], 1e-11);
  }
}

TEST_CASE("tensor_gamma: the two leading channels cancel to the 17e term") {
  // The three channels of the cos 2phi' harmonic stand
  // T_LL : T_LT : T_TT = 3 : -3 : 1 as gamma^2 and y go to zero, so the
  // leading-twist rate leakage and the twist-3 Eq. (17d) one cancel almost
  // exactly and what survives is the twist-4 Eq. (17e) term alone.  That
  // cancellation is the whole size of the effect: plans/08 D2 had guessed a
  // factor 6.9 the other way, from the T_LT term ADDING.  It is what puts
  // Delta_fake/(gamma^2 b1) at 0.14-0.16 -- the 1/6 of Eq. (17e) -- in place
  // of the retired bound "gamma^2 b1 x 1.15"
  // (`evgen/scripts/tensor_gamma_leakage.py`).
  struct Channels { double r_ll, r_lt, r_tot; };
  const auto channels = [](double x, double q2, double y) {
    const double g2 = gamma_squared(x, q2);
    const double eps = epsilon_gamma(y, g2);
    const std::pair<double, double> cs = theta_q_cos_sin(y, g2);
    const SFTables t = hand_tables(x);
    const CosynTensorSFs f = cosyn_tensor_sfs(t.b1, t.b2, t.b3, t.b4, x, g2);
    const double ll = 0.5 * cs.second * cs.second * (f.f_t + eps * f.f_l);
    const double lt = -0.5 * cs.first * cs.second
                      * std::sqrt(2.0 * eps * (1.0 + eps)) * f.f_lt;
    const double tt = 0.5 * (cs.first * cs.first + 1.0) * eps * f.f_tt;
    return Channels{ll / tt, lt / tt, (ll + lt + tt) / tt};
  };
  for (const Point& p : kPoints) {
    const Channels c = channels(p.x, p.q2, p.y);
    CHECK(2.9 < c.r_ll);
    CHECK(c.r_ll < 3.2);
    CHECK(-3.0 < c.r_lt);
    CHECK(c.r_lt < -2.8);
    CHECK(0.9 < c.r_tot);
    CHECK(c.r_tot < 1.2);
  }
  const Channels c = channels(0.05, 5.0, 1.0e-5);
  CHECK_CLOSE_AT(c.r_ll, 3.0, 0.0, 5e-3);
  CHECK_CLOSE_AT(c.r_lt, -3.0, 0.0, 5e-3);
  CHECK_CLOSE_AT(c.r_tot, 1.0, 0.0, 1e-2);
}

TEST_CASE("tensor_gamma: the cos 2phi leakage coefficient replaces the 1.15 bound") {
  // What `evgen/scripts/tensor_gamma_leakage.py` measures and prints: the
  // equivalent fake Delta the tensor rate sector leaks into the cos 2phi
  // amplitude of a TRANSVERSELY tensor-polarized fill,
  //
  //   Delta_fake = -a2(tensor) y^2 D_phi / (1 - y),
  //
  // as a coefficient of gamma^2 b1.  Until 2026-08-29 that size was quoted as
  // a BOUND, Delta_fake = 1.15 gamma^2 b1 -- the Eq. (17e) term alone,
  // gamma^2 b1/6, times the 6.9 the full combination was guessed to be, with
  // the T_LT term ADDING.  It does not add: it cancels the T_LL one (the
  // 3 : -3 : 1 ratio above), so what survives is the twist-4 Eq. (17e) term
  // almost alone and the coefficient sits at the 1/6 of that term.  The
  // script quotes 0.14-0.16 over the twelve money-plot-5 sweet spots of the
  // three 6Li configurations; the three points of this file, which are three
  // of those spots' neighbourhood plus a deliberately large-gamma^2 one, give
  // 0.1628 / 0.1601 / 0.1548 on the same toy table.
  //
  // SIGN.  With the literature convention (TENSOR_LL_SIGN = -1) and a
  // positive b1 the leakage is NEGATIVE at every spot -- OPPOSITE in sign to
  // the cos 2phi amplitude of the moment-constrained (negative-Delta) models
  // -- so it cancels part of the measured amplitude instead of faking one.
  // Flipping the constant back flips it, which is why D2 was gated on D1.
  const double want[3] = {0.1628248344088873, 0.16006748610799001,
                          0.15477697871375315};
  const InclusiveKernel kern = gamma_kernel();
  const EventSpinState st{0, 0.0, 1.0, 1.0, 0.5 * kPi, 0.0};
  for (int i = 0; i < 3; ++i) {
    const Point& p = kPoints[i];
    const SFTables t = hand_tables(p.x);
    const TensorHarmonics h =
        kern.tensor_harmonics_gamma(t, p.x, p.q2, p.y, st);
    const double dphi = t.f1 + (1.0 - p.y) / (p.x * p.y * p.y) * t.f2;
    const double d_fake = -h.h2 * p.y * p.y * dphi / (1.0 - p.y);
    const double coef = d_fake / (gamma_squared(p.x, p.q2) * t.b1);
    CHECK_CLOSE(coef, want[i], 1e-11);
    CHECK(0.14 < coef);
    CHECK(coef < 0.17);          // the band, not the retired 1.15
    CHECK(h.h2 < 0.0);           // negative with b1 > 0 and the -1 convention
  }
}

// --- the paper's own Table 1: the anchor that is not a re-typing ----------

TEST_CASE("tensor_gamma: Table 1 row 2, a target polarized along q") {
  // For b3 = b4 = 0, b2 = 2x b1 and a target polarized along q (Eq. 20a:
  // T_LL = Q/3, T_LT = T_TT = 0),
  //
  //   b1/(F1 A_T) = -9(1 + eps g2)/(6 + 5 g2 + 2 eps g2^2),  g2 = gamma^2.
  //
  // That row is a DERIVED result the module has never seen, and it is what
  // fixes the bracketing of Eq. (17a): the leading 2 multiplies b1 alone.
  // Applying it to the whole bracket gives 4 g2 for the 5 g2 and fails here by
  // up to 2.5 % over the range the W^2 cut allows.
  for (const Point& p : kPoints) {
    const double g2 = gamma_squared(p.x, p.q2);
    const double eps = epsilon_gamma(p.y, g2);
    const SFTables t = hand_tables(p.x, 0.02, 0.0, 0.0, 0.0);  // F2 = 2x F1
    const CosynTensorSFs f = cosyn_tensor_sfs(t.b1, t.b2, t.b3, t.b4, p.x, g2);
    const std::pair<double, double> fu =
        cosyn_unpolarized_sfs(t.f1, t.f2, p.x, g2);
    const double a_t =
        (2.0 / 3.0) * (f.f_t + eps * f.f_l) / (fu.first + eps * fu.second);
    const double printed =
        -9.0 * (1.0 + eps * g2) / (6.0 + 5.0 * g2 + 2.0 * eps * g2 * g2);
    CHECK_CLOSE(t.b1 / (t.f1 * a_t), printed, 1e-10);
  }
}

TEST_CASE("tensor_gamma: Table 1 row 3, a target polarized along the beam") {
  // The same quantity for a target polarized along the BEAM, which the module
  // reaches at theta_S = 0.  This row is the stronger anchor of the two -- it
  // carries F[U T_LT] and F[U T_TT] and the whole photon-frame geometry, and
  // it is what fixes the frame: the incoming lepton sits at +sin theta_q, so
  // that T_LT cos phi_TL = +(Q/2) cos theta_q sin theta_q for N = N_e.
  // Eq. (22b) as printed has the opposite sign and misses this row by up to a
  // factor of three at the largest gamma^2 the W^2 cut allows.
  const InclusiveKernel kern = gamma_kernel();
  const EventSpinState st{0, 0.0, 1.0, 1.0, 0.0, 0.0};  // N = N_e, Q = 1
  for (const Point& p : kPoints) {
    const double g2 = gamma_squared(p.x, p.q2);
    const double g = std::sqrt(g2);
    const double eps = epsilon_gamma(p.y, g2);
    const std::pair<double, double> cs = theta_q_cos_sin(p.y, g2);
    const double c = cs.first, sn = cs.second;
    const SFTables t = hand_tables(p.x, 0.02, 0.0, 0.0, 0.0);
    const TensorHarmonics h =
        kern.tensor_harmonics_gamma(t, p.x, p.q2, p.y, st);
    // w = -TENSOR_LL_SIGN * (numerator/denominator) and A_T is twice that
    // numerator over the same denominator, at Q = 1
    const double a_t = -2.0 * h.h0 / TENSOR_LL_SIGN;
    const double printed =
        18.0 * (1.0 + eps * g2)
        / ((3.0 * c * c - 1.0) * (-6.0 - 5.0 * g2 - 2.0 * eps * g2 * g2)
           - c * sn * std::sqrt(2.0 * eps * (1.0 + eps))
                 * (9.0 * g + 3.0 * g * g * g)
           - sn * sn * eps * 3.0 * g2);
    CHECK_CLOSE(t.b1 / (t.f1 * a_t), printed, 1e-10);
    CHECK(h.h1 == 0.0);
    CHECK(h.h2 == 0.0);
  }
}

TEST_CASE("tensor_gamma: b3 and b4 default to zero and reach the kernel") {
  const double s = s_mid_6li();
  const double x = 0.1, q2 = 1.14;
  InclusiveKernel::Options plain;
  plain.b1_func = toy_b1_func();
  const SFTables t0 = InclusiveKernel(LI6(), plain).tables(x, q2);
  CHECK(t0.b3 == 0.0);
  CHECK(t0.b4 == 0.0);

  const InclusiveKernel kern = gamma_kernel(
      toy_b1_func(), [](double, double, double f1) { return 0.05 * f1; });
  const EventSpinState st{0, 0.0, 1.0, 1.0, 0.5 * kPi, 0.0};
  const SFTables t = kern.tables(x, q2);
  CHECK(t.b3 != 0.0);
  const InclusiveKernel base = gamma_kernel(toy_b1_func());
  const double a2 = kern.amplitudes(t, x, q2, s, st).a2;
  const double a2_0 =
      base.amplitudes(base.tables(x, q2), x, q2, s, st).a2;
  CHECK(a2 != a2_0);
}
