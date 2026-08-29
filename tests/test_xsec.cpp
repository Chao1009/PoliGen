// Gate 2 (DEVELOPMENT_PLAN section 4.2): the doubly polarized master formula
// must reproduce the analytic asymmetries sector by sector, and the finite-
// gamma option must reproduce E143's own lab-frame factors.
//
// Ported from PolarizedLithiumSim/evgen/tests/test_xsec_identity.py,
// test_tensor_convention.py and test_target_mass.py.

#include <cmath>
#include <memory>
#include <vector>

#include "check_close.hpp"
#include "lipolgen/asymmetries.hpp"
#include "lipolgen/beams.hpp"
#include "lipolgen/numerics.hpp"
#include "lipolgen/sf.hpp"
#include "lipolgen/xsec.hpp"

using namespace lipolgen;

namespace {

constexpr double kRtol = 1e-12;

struct Point {
  double x, q2, y;
};

/// The analysis grid of test_xsec_identity.py: x in logspace(-3, -0.15, 14),
/// Q2 in logspace(0.1, 2.3, 11), keeping 0.01 < y < 0.95.
std::vector<Point> analysis_grid(double s) {
  std::vector<Point> out;
  const std::vector<double> xs = logspace(-3.0, -0.15, 14);
  const std::vector<double> q2s = logspace(0.1, 2.3, 11);
  for (double xv : xs) {
    for (double qv : q2s) {
      const double yv = qv / (s * xv);
      if (yv > 0.01 && yv < 0.95) out.push_back(Point{xv, qv, yv});
    }
  }
  return out;
}

/// The (x, Q2) points test_tensor_convention.py walks across the window.
const Point kConventionPoints[6] = {
    {0.0224, 1.14, 0.0}, {0.056, 1.14, 0.0}, {0.141, 3.13, 0.0},
    {0.141, 14.3, 0.0},  {0.005, 2.0, 0.0},  {0.30, 30.0, 0.0}};

double s_mid_6li() { return default_configs("6Li")[1].s_per_nucleon(); }

double eps_of(double y) { return (1.0 - y) / (1.0 - y + 0.5 * y * y); }

/// The scenario b1 / Delta shapes the Python identity tests use.
SFFunc3 toy_b1_func() {
  return [](double x, double q2, double f1) { return toy_b1(x, q2, f1); };
}
SFFunc3 toy_delta_func(double scale = 1e-3) {
  return [scale](double x, double q2, double f1) {
    return toy_delta_gluon(x, q2, f1, scale);
  };
}

/// F2 scaled by a constant: F1 and the rate move, a_2 must not.
class ScaledF2 : public ToyF2 {
 public:
  explicit ScaledF2(double k) : k_(k) {}
  double f2p(double x, double q2) const override {
    return k_ * ToyF2::f2p(x, q2);
  }
  double f2n(double x, double q2) const override {
    return f2p(x, q2) * f2n_over_f2p(x);
  }

 private:
  double k_;
};

RFunc const_r(double value) {
  return [value](double, double) { return value; };
}

/// The MASSLESS kernel: `InclusiveKernel::Options::target_mass` defaults to
/// TRUE since 2026-08-29 (`xsec.py`'s own default), so every identity written
/// against A_par = D(y) g1/F1 has to ask for the massless kernel by name.
InclusiveKernel massless_kernel(const Ion& ion) {
  InclusiveKernel::Options opt;
  opt.target_mass = false;
  return InclusiveKernel(ion, opt);
}

}  // namespace

// --- the three sector identities ----------------------------------------

TEST_CASE("vector sector: (w+ - w-)/(2 + w+ + w-) = P_e A_par") {
  const double s = s_mid_6li();
  const InclusiveKernel kern = massless_kernel(LI6());  // pure vector sector
  const double pe = 0.7;
  for (const Point& p : analysis_grid(s)) {
    const SFTables t = kern.tables(p.x, p.q2);
    const Amplitudes wp =
        kern.amplitudes(t, p.x, p.q2, s, EventSpinState{+1, pe, 1.0, 1.0, 0.0, 0.0});
    const Amplitudes wm =
        kern.amplitudes(t, p.x, p.q2, s, EventSpinState{-1, pe, 1.0, 1.0, 0.0, 0.0});
    const double measured =
        (wp.w_avg - wm.w_avg) / (2.0 + wp.w_avg + wm.w_avg);
    const double expected = pe * a_parallel(t.g1, t.f1, p.y, p.x, p.q2);
    CHECK_CLOSE(measured, expected, kRtol);
  }
}

TEST_CASE("vector sector for spin 3/2, and rank-2 slots default to zero") {
  const double s = s_mid_6li();
  const InclusiveKernel kern = massless_kernel(LI7());
  for (const Point& p : analysis_grid(s)) {
    const SFTables t = kern.tables(p.x, p.q2);
    // m = +3/2: full vector polarization m/J = 1
    const Amplitudes wp =
        kern.amplitudes(t, p.x, p.q2, s, EventSpinState{+1, 1.0, 1.5, 1.5, 0.0, 0.0});
    const Amplitudes wm =
        kern.amplitudes(t, p.x, p.q2, s, EventSpinState{-1, 1.0, 1.5, 1.5, 0.0, 0.0});
    CHECK_CLOSE((wp.w_avg - wm.w_avg) / (2.0 + wp.w_avg + wm.w_avg),
                a_parallel(t.g1, t.f1, p.y, p.x, p.q2), kRtol);
    // unpolarized electron on a zero rank-2 kernel -> no modulation at all
    const Amplitudes w0 =
        kern.amplitudes(t, p.x, p.q2, s, EventSpinState{0, 0.0, 1.5, 1.5, 0.0, 0.0});
    CHECK(w0.w_avg == 0.0);
    CHECK(w0.a1 == 0.0);
    CHECK(w0.a2 == 0.0);
  }
}

TEST_CASE("tensor sector: the thirds combination is A_zz") {
  const double s = s_mid_6li();
  InclusiveKernel::Options opt;
  opt.b1_func = toy_b1_func();
  const InclusiveKernel kern(LI6(), opt);
  for (const Point& p : analysis_grid(s)) {
    const SFTables t = kern.tables(p.x, p.q2);
    double sigma[3];
    int i = 0;
    for (double m : {1.0, 0.0, -1.0}) {
      sigma[i++] = 1.0 + kern.amplitudes(t, p.x, p.q2, s,
                                         EventSpinState{0, 0.0, 1.0, m, 0.0, 0.0})
                             .w_avg;
    }
    const double measured =
        (sigma[0] + sigma[2] - 2 * sigma[1]) / (sigma[0] + sigma[2] + sigma[1]);
    CHECK_CLOSE(measured, azz(t.b1, t.f1, t.f2, p.x, p.y), kRtol);
  }
}

TEST_CASE("transverse tensor: a_2 is A_cos2phi") {
  const double s = s_mid_6li();
  InclusiveKernel::Options opt;
  opt.delta_func = toy_delta_func();
  const InclusiveKernel kern(LI6(), opt);
  for (const Point& p : analysis_grid(s)) {
    const SFTables t = kern.tables(p.x, p.q2);
    for (double m : {1.0, -1.0}) {
      const Amplitudes a = kern.amplitudes(
          t, p.x, p.q2, s, EventSpinState{0, 0.0, 1.0, m, kPi / 2.0, 0.0});
      CHECK(a.w_avg == 0.0);  // no b1 -> no phi-averaged tensor shift
      CHECK_CLOSE(a.a2, a_cos2phi(t.delta, t.f1, t.f2, p.x, p.y), kRtol);
    }
  }
}

TEST_CASE("the population-averaged cross section is the unpolarized one") {
  const double s = s_mid_6li();
  InclusiveKernel::Options opt;
  opt.b1_func = toy_b1_func();
  opt.delta_func = [](double x, double q2, double f1) {
    return toy_delta_gluon(x, q2, f1);
  };
  const InclusiveKernel kern(LI6(), opt);
  const std::vector<Point> grid = analysis_grid(s);
  for (std::size_t k = 0; k < grid.size(); k += 7) {
    const Point& p = grid[k];
    double tot = 0.0;
    for (double m : {1.0, 0.0, -1.0}) {
      const EventSpinState st{0, 0.0, 1.0, m, 0.7, 1.3};
      double mean = 0.0;
      for (int i = 0; i < 48; ++i) {
        const double phi = 2.0 * kPi * i / 48.0;
        mean += kern.dsigma(p.x, p.q2, phi, s, st);
      }
      tot += (mean / 48.0) / 3.0;
    }
    CHECK_CLOSE(tot * 2 * kPi, kern.dsigma_unpol(p.x, p.q2, s), 1e-9);
  }
}

TEST_CASE("a_perp is the leading-gamma transverse amplitude") {
  // Against an independent transcription of d(y) * gamma * ((y/2) g1 + g2)/F1
  // (E143 conventions).
  const double s = s_mid_6li();
  const InclusiveKernel kern(LI6());  // g2 = g2_WW default
  const double x = 0.2;
  for (double q2 : {5.0, 20.0, 80.0}) {
    const double y = q2 / (s * x);
    const SFTables t = kern.tables(x, q2, true);
    const Amplitudes a = kern.amplitudes(
        t, x, q2, s, EventSpinState{+1, 1.0, 1.0, 1.0, kPi / 2.0, 0.0}, true);
    const double eps = eps_of(y);
    const double d = depolarization_d(y, x, q2) * std::sqrt(2.0 * eps / (1.0 + eps));
    const double gamma = 2.0 * 0.9383 * x / std::sqrt(q2);
    const double expected = d * gamma * (0.5 * y * t.g1 + t.g2) / t.f1;
    CHECK_CLOSE(a.a1, expected, kRtol);
    CHECK(expected != 0.0);
    // and the shared free function agrees with the kernel's own
    CHECK_CLOSE(a.a1, a_perp(t.g1, t.g2, t.f1, y, x, q2), kRtol);
    // a longitudinal axis gives no cos(phi) modulation
    const Amplitudes al = kern.amplitudes(
        t, x, q2, s, EventSpinState{+1, 1.0, 1.0, 1.0, 0.0, 0.0}, true);
    CHECK(al.a1 == 0.0);
  }
}

// --- the tensor convention ------------------------------------------------

TEST_CASE("Cosyn Eq. 27: A_zz(theta_S=0)(1 + eps R) = sign (2/3) b1/F1") {
  // EXACTLY and at every y.  With b2 = 2x b1 and F2 = 2x(1+R) F1 the
  // (1-y)/(x y^2) terms of numerator and denominator combine into
  // 2(1-y+y^2/2) and 2[(1-y+y^2/2) + R(1-y)], whose ratio is 1/(1 + eps R).
  for (const Point& p : kConventionPoints) {
    const double f1 = 1.0;   // the identity is F1-free
    const double b1 = 0.037;  // any value; only the ratio enters
    const double r = r_sigma_lt(p.x, p.q2);
    const double f2 = 2.0 * p.x * (1.0 + r) * f1;
    for (double y : {0.01, 0.05, 0.2, 0.6, 0.9}) {
      const double got = azz(b1, f1, f2, p.x, y) * (1.0 + eps_of(y) * r);
      CHECK_CLOSE(got, TENSOR_LL_SIGN * (2.0 / 3.0) * b1 / f1, kRtol);
    }
  }
}

TEST_CASE("the program sign is opposite to the literature, deliberately") {
  // Stated as a test so that flipping TENSOR_LL_SIGN is a deliberate act with
  // a visible consequence (plans/08 D1).
  const double x = 0.056, q2 = 1.14, y = 0.05, b1 = 0.037;
  const double r = r_sigma_lt(x, q2);
  const double f2 = 2.0 * x * (1.0 + r);
  const double cosyn = -(2.0 / 3.0) * b1;  // Cosyn Eq. (27) with F1 = 1
  const double program = azz(b1, 1.0, f2, x, y) * (1.0 + eps_of(y) * r);
  CHECK(TENSOR_LL_SIGN == +1.0);
  CHECK_CLOSE(program, -cosyn, kRtol);
}

TEST_CASE("the kernel thirds combination carries the same sign as A_zz") {
  const double s = s_mid_6li();
  InclusiveKernel::Options opt;
  opt.b1_func = toy_b1_func();
  const InclusiveKernel kern(LI6(), opt);
  const double x = 0.056, q2 = 1.14;
  const double y = q2 / (s * x);
  const SFTables t = kern.tables(x, q2);
  double w[3];
  int i = 0;
  for (double m : {1.0, 0.0, -1.0}) {
    w[i++] = kern.amplitudes(t, x, q2, s, EventSpinState{0, 0.0, 1.0, m, 0.0, 0.0})
                 .w_avg;
  }
  const double expected = azz(t.b1, t.f1, t.f2, x, y, &t.b2);
  CHECK_CLOSE((w[0] + w[2] - 2.0 * w[1]) / (3.0 + w[0] + w[2] + w[1]), expected,
              kRtol);
  // A_zz is the thirds combination, not a single state: the m = +-1 shift
  // alone is half of it
  CHECK_CLOSE(w[0], 0.5 * expected, kRtol);
  CHECK(((w[0] > 0) == (TENSOR_LL_SIGN * t.b1 > 0)));
}

TEST_CASE("the Delta sector does not depend on the tensor convention") {
  const double s = s_mid_6li();
  InclusiveKernel::Options opt;
  opt.b1_func = toy_b1_func();
  opt.delta_func = [](double, double, double f1) { return -1e-2 * f1; };
  const InclusiveKernel kern(LI6(), opt);
  const double x = 0.056, q2 = 1.14;
  const double y = q2 / (s * x);
  const SFTables t = kern.tables(x, q2);
  const Amplitudes a = kern.amplitudes(
      t, x, q2, s, EventSpinState{0, 0.0, 1.0, 1.0, kPi / 2.0, 0.0});
  CHECK_CLOSE(a.a2, a_cos2phi(t.delta, t.f1, t.f2, x, y), kRtol);
  CHECK(a.a2 > 0.0);  // Delta < 0 with c_m = +1 -> a2 > 0
}

// --- one rank-2 geometry for both spins ----------------------------------

TEST_CASE("the rank-2 geometry is one formula for both spins") {
  // Cosyn Eq. (9): t_ij = (Q_NN/2)(3 n_i n_j - d_ij) for any J, so
  // T_LL = Q_NN P_2 and T_TT = (3/2) Q_NN sin^2.  The kernel returns
  // (Q_NN, 3 Q_NN) and builds both channels from it.
  const std::vector<std::pair<Ion, std::vector<double>>> cases = {
      {LI6(), {1.0, 0.0, -1.0}}, {LI7(), {1.5, 0.5, -0.5, -1.5}}};
  for (const auto& c : cases) {
    const InclusiveKernel kern(c.first);
    const double j = c.first.spin;
    for (double m : c.second) {
      const std::pair<double, double> qc = kern.tensor_moments(m);
      CHECK_CLOSE(qc.first, (3.0 * m * m - j * (j + 1.0)) / 3.0, kRtol);
      CHECK_CLOSE(qc.second, 3.0 * qc.first, kRtol);
    }
  }
  // spin 1/2 carries no rank-2 alignment
  const InclusiveKernel kp(PROTON());
  CHECK(kp.tensor_moments(0.5).first == 0.0);
  CHECK(kp.tensor_moments(0.5).second == 0.0);
}

TEST_CASE("the J = 1 geometry is the HJM transcription digit for digit") {
  const double s = s_mid_6li();
  InclusiveKernel::Options opt;
  opt.b1_func = toy_b1_func();
  const InclusiveKernel kern(LI6(), opt);
  for (double m : {1.0, 0.0, -1.0}) {
    const double c_m = 3.0 * m * m - 2.0;  // (1, -2, 1)
    const std::pair<double, double> qc = kern.tensor_moments(m);
    CHECK_CLOSE(qc.first, c_m / 3.0, 1e-13);
    CHECK_CLOSE(qc.second, c_m, 1e-13);
  }
  const double x = 0.056, q2 = 1.14;
  const SFTables t = kern.tables(x, q2);
  const Amplitudes a =
      kern.amplitudes(t, x, q2, s, EventSpinState{0, 0.0, 1.0, 1.0, 0.0, 0.0});
  CHECK_CLOSE(a.w_avg,
              0.5 * azz(t.b1, t.f1, t.f2, x, q2 / (s * x), &t.b2), 1e-13);
}

TEST_CASE("the spin-3/2 rate and cos-2phi channels are mutually consistent") {
  // CHARACTERIZATION, not a physics assertion: the spin-3/2 rank-2
  // normalization is plans/04 #14.  What is pinned is INTERNAL consistency --
  // the transverse/longitudinal ratio is Q_NN-free, so it must be
  // spin-independent.
  const double s = s_mid_6li();
  auto ratio = [s](const Ion& ion, double j) {
    InclusiveKernel::Options opt;
    opt.b1_func = [](double, double, double f1) { return 0.05 * f1; };
    opt.delta_func = [](double, double, double f1) { return -1e-2 * f1; };
    opt.b1_32_func = opt.b1_func;
    opt.delta_32_func = opt.delta_func;
    const InclusiveKernel kern(ion, opt);
    const double x = 0.056, q2 = 1.14;
    const SFTables t = kern.tables(x, q2);
    const double lon =
        kern.amplitudes(t, x, q2, s, EventSpinState{0, 0.0, j, j, 0.0, 0.0}).w_avg;
    const double tra =
        kern.amplitudes(t, x, q2, s, EventSpinState{0, 0.0, j, j, kPi / 2.0, 0.0})
            .a2;
    return tra / lon;
  };
  CHECK_CLOSE(ratio(LI7(), 1.5), ratio(LI6(), 1.0), 1e-10);
}

// --- b2, the axis angle and the magic angle -------------------------------

TEST_CASE("b2_func overrides the default b2 = 2 x b1") {
  const double s = s_mid_6li();
  const double xs[3] = {0.02, 0.056, 0.14};
  const double q2s[3] = {1.14, 1.14, 3.13};
  InclusiveKernel::Options def_opt, ovr_opt;
  def_opt.b1_func = [](double, double, double f1) { return 0.05 * f1; };
  ovr_opt.b1_func = def_opt.b1_func;
  ovr_opt.b2_func = [](double x, double, double f1) { return 5.0 * x * 0.05 * f1; };
  const InclusiveKernel def(LI6(), def_opt);
  const InclusiveKernel ovr(LI6(), ovr_opt);
  for (int i = 0; i < 3; ++i) {
    const SFTables td = def.tables(xs[i], q2s[i]);
    const SFTables to = ovr.tables(xs[i], q2s[i]);
    CHECK_CLOSE(td.b2, 2.0 * xs[i] * td.b1, kRtol);
    CHECK_CLOSE(to.b2, 5.0 * xs[i] * to.b1, kRtol);
    CHECK(td.b2 != to.b2);
    const double y = q2s[i] / (s * xs[i]);
    const EventSpinState st{0, 0.0, 1.0, 1.0, 0.0, 0.0};
    const double wo = ovr.amplitudes(to, xs[i], q2s[i], s, st).w_avg;
    CHECK_CLOSE(wo, 0.5 * azz(to.b1, to.f1, to.f2, xs[i], y, &to.b2), kRtol);
    CHECK(def.amplitudes(td, xs[i], q2s[i], s, st).w_avg != wo);
  }
}

TEST_CASE("the tensor rate follows P_2(cos theta_S), and the magic angle kills it") {
  const double s = s_mid_6li();
  InclusiveKernel::Options opt;
  opt.b1_func = [](double, double, double f1) { return 0.05 * f1; };
  const InclusiveKernel kern(LI6(), opt);
  const double x = 0.056, q2 = 1.14;
  const double y = q2 / (s * x);
  const SFTables t = kern.tables(x, q2);
  const double magic = std::acos(1.0 / std::sqrt(3.0));
  for (double th : {0.0, 0.4, 0.7, magic, kPi / 2.0}) {
    const Amplitudes a =
        kern.amplitudes(t, x, q2, s, EventSpinState{0, 0.0, 1.0, 1.0, th, 0.0});
    CHECK_CLOSE_AT(a.w_avg,
                   0.5 * azz(t.b1, t.f1, t.f2, x, y, &t.b2, th), kRtol, 1e-300);
  }
  const Amplitudes am =
      kern.amplitudes(t, x, q2, s, EventSpinState{0, 0.0, 1.0, 1.0, magic, 0.0});
  CHECK(std::fabs(am.w_avg) < 1e-15);
}

// --- R in the cos 2phi channel (plans/08 C2) ------------------------------

TEST_CASE("r_func = r_sigma_lt is bit-for-bit r_func = none") {
  const double s = s_mid_6li();
  InclusiveKernel::Options a_opt, b_opt;
  a_opt.b1_func = toy_b1_func();
  a_opt.delta_func = [](double x, double q2, double f1) {
    return toy_delta_gluon(x, q2, f1);
  };
  b_opt = a_opt;
  b_opt.r_func = [](double x, double q2) { return r_sigma_lt(x, q2); };
  const InclusiveKernel ka(LI6(), a_opt);
  const InclusiveKernel kb(LI6(), b_opt);
  for (const Point& p : analysis_grid(s)) {
    const SFTables ta = ka.tables(p.x, p.q2, true);
    const SFTables tb = kb.tables(p.x, p.q2, true);
    CHECK(ta.f1 == tb.f1);
    CHECK(ta.f2 == tb.f2);
    CHECK(ta.g1 == tb.g1);
    CHECK(ta.g2 == tb.g2);
    CHECK(ta.b1 == tb.b1);
    CHECK(ta.b2 == tb.b2);
    CHECK(ta.delta == tb.delta);
    const EventSpinState st{+1, 0.7, 1.0, 1.0, 1.0, 0.3};
    const Amplitudes aa = ka.amplitudes(ta, p.x, p.q2, s, st, true);
    const Amplitudes ab = kb.amplitudes(tb, p.x, p.q2, s, st, true);
    CHECK(aa.w_avg == ab.w_avg);
    CHECK(aa.a1 == ab.a1);
    CHECK(aa.a2 == ab.a2);
    CHECK(ka.dsigma_unpol(p.x, p.q2, s) == kb.dsigma_unpol(p.x, p.q2, s));
  }
}

TEST_CASE("the cos 2phi amplitude carries R and NOT F2") {
  // a_2 = -(eps/2) c_eff sin^2(theta_S) (Delta/F1)/(1 + eps R), written out by
  // hand from eps and F1 = F2/(2x(1+R)) -- no F2 anywhere in it.
  const double s = s_mid_6li();
  const double d_over_f1 = 1.3e-3;
  const EventSpinState st{0, 0.0, 1.0, 1.0, kPi / 2.0, 0.0};  // c_eff = 1
  for (double r_value : {0.0, 0.18, 0.40}) {
    InclusiveKernel::Options opt;
    opt.delta_func = [d_over_f1](double, double, double f1) {
      return d_over_f1 * f1;
    };
    opt.r_func = const_r(r_value);
    const InclusiveKernel kern(LI6(), opt);
    for (const Point& p : analysis_grid(s)) {
      const SFTables t = kern.tables(p.x, p.q2);
      const double expected =
          -(eps_of(p.y) / 2.0) * d_over_f1 / (1.0 + eps_of(p.y) * r_value);
      CHECK_CLOSE(kern.amplitudes(t, p.x, p.q2, s, st).a2, expected, kRtol);
    }
  }
}

TEST_CASE("the cos 2phi amplitude is independent of the F2 backend") {
  const double s = s_mid_6li();
  const EventSpinState st{0, 0.0, 1.0, 1.0, kPi / 2.0, 0.0};
  const std::vector<Point> grid = analysis_grid(s);
  for (std::size_t k = 0; k < grid.size(); k += 5) {
    const Point& p = grid[k];
    double f1_ref = 0.0, a2_ref = 0.0, sig_ref = 0.0;
    for (double kfac : {1.0, 1.7}) {
      InclusiveKernel::Options opt;
      opt.f2_source = std::make_shared<const ScaledF2>(kfac);
      opt.delta_func = [](double, double, double f1) { return 1.3e-3 * f1; };
      const InclusiveKernel kern(LI6(), opt);
      const SFTables t = kern.tables(p.x, p.q2);
      const double a2 = kern.amplitudes(t, p.x, p.q2, s, st).a2;
      const double sig = kern.dsigma_unpol(p.x, p.q2, s);
      if (kfac == 1.0) {
        f1_ref = t.f1;
        a2_ref = a2;
        sig_ref = sig;
      } else {
        CHECK_CLOSE(t.f1, 1.7 * f1_ref, kRtol);
        CHECK_CLOSE(sig, 1.7 * sig_ref, kRtol);
        CHECK_CLOSE(a2, a2_ref, 1e-14);
      }
    }
  }
}

TEST_CASE("R1998 rescales a_2 by (1 + eps R_toy)/(1 + eps R_1998) and nothing else") {
  const double s = s_mid_6li();
  const EventSpinState st{0, 0.0, 1.0, 1.0, kPi / 2.0, 0.0};
  InclusiveKernel::Options o_toy, o_98;
  o_toy.delta_func = [](double, double, double f1) { return 1.3e-3 * f1; };
  o_98 = o_toy;
  o_98.r_func = [](double x, double q2) { return r1998(x, q2); };
  const InclusiveKernel k_toy(LI6(), o_toy);
  const InclusiveKernel k_98(LI6(), o_98);
  bool saw_big = false;
  for (const Point& p : analysis_grid(s)) {
    const double a_toy =
        k_toy.amplitudes(k_toy.tables(p.x, p.q2), p.x, p.q2, s, st).a2;
    const double a_98 =
        k_98.amplitudes(k_98.tables(p.x, p.q2), p.x, p.q2, s, st).a2;
    const double eps = eps_of(p.y);
    const double expected = (1.0 + eps * r_sigma_lt(p.x, p.q2))
                            / (1.0 + eps * r1998(p.x, p.q2));
    CHECK_CLOSE(a_98 / a_toy, expected, kRtol);
    if (std::fabs(expected - 1.0) > 0.1) saw_big = true;
  }
  CHECK(saw_big);  // the swap is not cosmetic
}

// --- the finite-gamma (target-mass) option --------------------------------

TEST_CASE("the finite-gamma factors match E143's lab-frame definitions") {
  // eps = 1/[1 + 2(1 + nu^2/Q^2) tan^2(theta/2)],
  // D   = (1 - E' eps/E)/(1 + eps R),  eta = eps sqrt(Q^2)/(E - E' eps).
  // The two forms are the same object written in different variables, so
  // nothing but double-precision agreement is acceptable.
  const double r = 0.18;
  for (int i = 0; i < 64; ++i) {
    const double e_in = 5.0 + 45.0 * (i + 0.5) / 64.0;
    const double e_out = e_in * (0.05 + 0.9 * ((i * 7) % 64 + 0.5) / 64.0);
    const double theta = 0.05 + 1.15 * ((i * 13) % 64 + 0.5) / 64.0;
    const double nu = e_in - e_out;
    const double q2 = 4.0 * e_in * e_out * std::pow(std::sin(theta / 2.0), 2.0);
    const double x = q2 / (2.0 * M_NUCLEON * nu);
    const double y = nu / e_in;
    const double g2v = gamma_squared(x, q2);
    CHECK_CLOSE(g2v, q2 / (nu * nu), kRtol);
    const double eps_lab =
        1.0 / (1.0 + 2.0 * (1.0 + nu * nu / q2) * std::pow(std::tan(theta / 2.0), 2.0));
    const double d_lab = (1.0 - e_out * eps_lab / e_in) / (1.0 + eps_lab * r);
    const double eta_lab = eps_lab * std::sqrt(q2) / (e_in - e_out * eps_lab);
    CHECK_CLOSE(epsilon_gamma(y, g2v), eps_lab, kRtol);
    CHECK_CLOSE(depolarization_gamma(y, g2v, r), d_lab, kRtol);
    CHECK_CLOSE(eta_gamma(y, g2v), eta_lab, kRtol);
  }
}

TEST_CASE("at gamma = 0 the finite-gamma set IS the fast simulation's") {
  const double x = 0.1, q2 = 4.0;
  const double r = r_sigma_lt(x, q2);
  for (double y : {0.01, 0.05, 0.2, 0.5, 0.9}) {
    CHECK_CLOSE(epsilon_gamma(y, 0.0), eps_of(y), 1e-14);
    CHECK_CLOSE(depolarization_gamma(y, 0.0, r), depolarization_d(y, x, q2),
                1e-14);
    CHECK(eta_gamma(y, 0.0) == 0.0);
  }
}

TEST_CASE("target_mass off is bit-for-bit the published kernel") {
  const InclusiveKernel k0 = massless_kernel(LI6());
  InclusiveKernel::Options tm;
  tm.target_mass = true;
  const InclusiveKernel k1(LI6(), tm);
  const double s = s_mid_6li();
  for (const Point& p : analysis_grid(s)) {
    const SFTables t0 = k0.tables(p.x, p.q2);
    const SFTables t1 = k1.tables(p.x, p.q2);
    CHECK(!t0.has_g2);  // nothing extra is computed
    CHECK(t1.has_g2);   // ... and it IS when needed
    const double y = 0.2;
    CHECK(k0.a_parallel(t0, p.x, p.q2, y)
          == a_parallel(t0.g1, t0.f1, y, p.x, p.q2));
    // the two kernels share every unpolarized / tensor number
    CHECK(t0.f1 == t1.f1);
    CHECK(t0.f2 == t1.f2);
    CHECK(t0.g1 == t1.g1);
    CHECK(t0.b1 == t1.b1);
    CHECK(t0.b2 == t1.b2);
    CHECK(t0.delta == t1.delta);
  }
}

// C2b: `target_mass = true` with `G2Mode::kZero` is the twist-3 g2_scale = 0
// variation and is PERMITTED (xsec.py:126-129 zeroes g2 in the tables and
// proceeds); only handing the finite-gamma kernel a table with no g2 slot at
// all is an error.
TEST_CASE("target_mass with g2_mode = zero is the g2_scale = 0 variation") {
  InclusiveKernel::Options z;
  z.target_mass = true;
  z.g2_mode = G2Mode::kZero;
  const InclusiveKernel kz(LI6(), z);            // no longer refused
  InclusiveKernel::Options s0;
  s0.target_mass = true;
  s0.g2_scale = 0.0;
  const InclusiveKernel k0s(LI6(), s0);
  InclusiveKernel::Options tmw;
  tmw.target_mass = true;
  const InclusiveKernel k1(LI6(), tmw);          // g2 = g2_WW, scale 1
  InclusiveKernel::Options s15;
  s15.target_mass = true;
  s15.g2_scale = 1.5;
  const InclusiveKernel k15(LI6(), s15);

  bool moved = false;
  for (double q2 : {2.0, 5.0, 20.0}) {
    for (double x : {0.05, 0.2, 0.5}) {
      const SFTables tz = kz.tables(x, q2);
      const SFTables t0 = k0s.tables(x, q2);
      const SFTables t1 = k1.tables(x, q2);
      const SFTables t15 = k15.tables(x, q2);
      // g2_mode = zero and g2_scale = 0 are the SAME variation, bit for bit
      CHECK(tz.has_g2);
      CHECK(tz.g2 == 0.0);
      CHECK(t0.g2 == 0.0);
      // g2_scale multiplies the WW table exactly
      CHECK_CLOSE(t15.g2, 1.5 * t1.g2, kRtol);
      // ... and it is a real handle on A_par
      const double y = 0.3;
      const double a_zero = kz.a_parallel(tz, x, q2, y);
      const double a_ww = k1.a_parallel(t1, x, q2, y);
      CHECK(std::isfinite(a_zero));
      CHECK(k0s.a_parallel(t0, x, q2, y) == a_zero);
      if (std::fabs(a_ww / a_zero - 1.0) > 1e-6) moved = true;
    }
  }
  CHECK(moved);
  // the finite-gamma kernel still refuses a table with NO g2 slot
  CHECK_THROWS(k1.a_parallel(massless_kernel(LI6()).tables(0.2, 5.0), 0.2, 5.0,
                             0.3));
}

// C2: the DEFAULT-constructed kernel is the Python default -- target mass ON.
TEST_CASE("the default kernel is the Python default: target_mass on") {
  const InclusiveKernel kdef(LI6());
  CHECK(kdef.target_mass());
  CHECK(kdef.g2_scale() == 1.0);
  InclusiveKernel::Options tm2;
  tm2.target_mass = true;
  const InclusiveKernel ktm(LI6(), tm2);
  const double s = s_mid_6li();
  for (const Point& p : analysis_grid(s)) {
    const SFTables td = kdef.tables(p.x, p.q2);
    const SFTables tt = ktm.tables(p.x, p.q2);
    CHECK(td.has_g2);
    CHECK(td.g2 == tt.g2);
    CHECK(kdef.a_parallel(td, p.x, p.q2, 0.3)
          == ktm.a_parallel(tt, p.x, p.q2, 0.3));
  }
  CHECK(!massless_kernel(LI6()).target_mass());
}

// P8: a spin state whose J is not the kernel's ion spin is refused -- the
// rank-2 branch is gated on state.j while tensor_moments reads ion().spin.
TEST_CASE("amplitudes refuse a spin state of the wrong J") {
  const double s = s_mid_6li();
  const InclusiveKernel k6(LI6());   // J = 1
  const InclusiveKernel k7(LI7());   // J = 3/2
  const SFTables t6 = k6.tables(0.2, 5.0);
  const SFTables t7 = k7.tables(0.2, 5.0);
  CHECK_THROWS(k6.amplitudes(t6, 0.2, 5.0, s,
                             EventSpinState{+1, 1.0, 1.5, 1.5, 0.0, 0.0}));
  CHECK_THROWS(k7.amplitudes(t7, 0.2, 5.0, s,
                             EventSpinState{+1, 1.0, 1.0, 1.0, 0.0, 0.0}));
  CHECK_THROWS(k6.amplitudes(t6, 0.2, 5.0, s,
                             EventSpinState{+1, 1.0, 0.5, 0.5, 0.0, 0.0}));
  // the matching state is fine
  CHECK_NOTHROW(k6.amplitudes(t6, 0.2, 5.0, s,
                              EventSpinState{+1, 1.0, 1.0, 1.0, 0.0, 0.0}));
  CHECK_NOTHROW(k7.amplitudes(t7, 0.2, 5.0, s,
                              EventSpinState{+1, 1.0, 1.5, 0.5, 0.0, 0.0}));
  // ... and so is a J = 1 CHANNEL state on a spin-1/2 kernel: that is the
  // d(e,e'p) control, where m_S labels the S_c = 1 channel spin of p (x) n
  // while the DIS target is the struck NEUTRON.  The kernel has no rank-2
  // sector at all there, so the tensor term is identically zero.
  const Ion neutron{"n", 1, 0, 0.5, 0.0, 1.0};  // tagged.hpp NEUTRON_TARGET
  const InclusiveKernel kn(neutron);
  const SFTables tn = kn.tables(0.2, 5.0);
  for (double m : {1.0, 0.0, -1.0}) {
    Amplitudes a;
    CHECK_NOTHROW(a = kn.amplitudes(tn, 0.2, 5.0, s,
                                    EventSpinState{+1, 1.0, 1.0, m, 0.4, 0.0}));
    CHECK(a.a2 == 0.0);
    // the whole w_avg is the vector term m/J times the longitudinal asymmetry
    const double y = 5.0 / (s * 0.2);
    CHECK_CLOSE_AT(a.w_avg, m * std::cos(0.4)
                                * kn.a_parallel(tn, 0.2, 5.0, y), kRtol, 1e-300);
  }
}

TEST_CASE("the target-mass flag moves A_par by O(gamma^2) and by nothing else") {
  const InclusiveKernel k0 = massless_kernel(LI6());
  InclusiveKernel::Options tm;
  tm.target_mass = true;
  const InclusiveKernel k1(LI6(), tm);
  // the DIS region W^2 >= 10 GeV^2
  std::vector<Point> grid;
  for (double lx = -3.0; lx <= std::log10(0.7) + 1e-12; lx += 0.1) {
    for (double lq = 0.0; lq <= 2.0 + 1e-12; lq += 0.2) {
      const double x = std::pow(10.0, lx), q2 = std::pow(10.0, lq);
      if (q2 * (1.0 - x) / x + M_NUCLEON * M_NUCLEON >= 10.0) {
        grid.push_back(Point{x, q2, 0.0});
      }
    }
  }
  CHECK(grid.size() > 50);
  for (double yy : {0.02, 0.1, 0.3, 0.6, 0.9}) {
    for (const Point& p : grid) {
      const SFTables t = k1.tables(p.x, p.q2);
      if (std::fabs(t.g2 / t.g1) > 3.0) continue;  // away from the g1 zero
      const double shift = k1.a_parallel(t, p.x, p.q2, yy)
                               / k0.a_parallel(t, p.x, p.q2, yy) - 1.0;
      CHECK(std::fabs(shift) <= 3.0 * gamma_squared(p.x, p.q2));
    }
  }
  // it scales EXACTLY as 1/Q^2 at fixed (x, y)
  const double xy[3][2] = {{0.3, 0.1}, {0.1, 0.05}, {0.5, 0.3}};
  for (const auto& c : xy) {
    double sh[2];
    int i = 0;
    for (double qq : {2.0, 200.0}) {
      const SFTables t = k1.tables(c[0], qq);
      sh[i++] = k1.a_parallel(t, c[0], qq, c[1])
                    / k0.a_parallel(t, c[0], qq, c[1]) - 1.0;
    }
    CHECK_CLOSE_AT(sh[0] / sh[1], 100.0, 0.05, 0.0);
  }
}

TEST_CASE("at small y the target-mass shift collapses to (1 + gamma^2)") {
  const InclusiveKernel k0 = massless_kernel(LI6());
  InclusiveKernel::Options tm;
  tm.target_mass = true;
  const InclusiveKernel k1(LI6(), tm);
  std::vector<double> devs;
  for (double yy : {0.05, 0.01}) {
    double worst = 0.0;
    for (double lx = -3.0; lx <= std::log10(0.7); lx += 0.1) {
      for (double lq = 0.0; lq <= 2.0; lq += 0.2) {
        const double x = std::pow(10.0, lx), q2 = std::pow(10.0, lq);
        if (q2 * (1.0 - x) / x + M_NUCLEON * M_NUCLEON < 10.0) continue;
        const SFTables t = k1.tables(x, q2);
        if (std::fabs(t.g2 / t.g1) > 3.0) continue;
        const double g2v = gamma_squared(x, q2);
        const double shift =
            k1.a_parallel(t, x, q2, yy) / k0.a_parallel(t, x, q2, yy) - 1.0;
        worst = std::max(worst, std::fabs(shift / g2v - 1.0));
      }
    }
    devs.push_back(worst);
  }
  CHECK(devs[0] <= 0.08);  // y <= 0.05: the collapse holds to 8 %
  CHECK(devs[1] <= 0.02);  // ... and improves linearly in y
  CHECK(devs[0] / devs[1] > 4.0);
}

TEST_CASE("the target-mass flag leaves the tensor and unpolarized sectors alone") {
  const double s = s_mid_6li();
  InclusiveKernel::Options o0, o1;
  o0.b1_func = toy_b1_func();
  o0.target_mass = false;
  o1 = o0;
  o1.target_mass = true;
  const InclusiveKernel k0(LI6(), o0);
  const InclusiveKernel k1(LI6(), o1);
  const EventSpinState unpol{0, 0.0, 1.0, 1.0, 0.4, 0.0};
  const EventSpinState lam{+1, 0.7, 1.0, 1.0, 0.0, 0.0};
  bool moved = false;
  for (const Point& p : analysis_grid(s)) {
    const SFTables t0 = k0.tables(p.x, p.q2);
    const SFTables t1 = k1.tables(p.x, p.q2);
    const Amplitudes a0 = k0.amplitudes(t0, p.x, p.q2, s, unpol);
    const Amplitudes a1 = k1.amplitudes(t1, p.x, p.q2, s, unpol);
    CHECK(a0.w_avg == a1.w_avg);
    CHECK(a0.a1 == a1.a1);
    CHECK(a0.a2 == a1.a2);
    const double w0 = k0.amplitudes(t0, p.x, p.q2, s, lam).w_avg;
    const double w1 = k1.amplitudes(t1, p.x, p.q2, s, lam).w_avg;
    if (w0 != w1) moved = true;
    // the vector-L term moves by the a_parallel shift and nothing else
    CHECK_CLOSE_AT(w1 - w0,
                   0.7 * (k1.a_parallel(t1, p.x, p.q2, p.y)
                          - k0.a_parallel(t0, p.x, p.q2, p.y)),
                   1e-6, 1e-15);
  }
  CHECK(moved);
}

// --- the exact positivity minimum over phi --------------------------------

TEST_CASE("density_min is the exact minimum over phi, not the envelope") {
  // The guard uses the EXACT minimum of 1 + A cos phi + B cos 2phi, not the
  // accept-reject bound 1 + |A| + |B|; checked against a fine scan.
  unsigned seed = 3u;
  auto rnd = [&seed]() {
    seed = seed * 1103515245u + 12345u;
    return -2.0 + 4.0 * ((seed >> 8) % 100000u) / 100000.0;
  };
  for (int i = 0; i < 300; ++i) {
    const double a = rnd(), b = rnd();
    double brute = 1e30;
    for (int k = 0; k <= 8000; ++k) {
      const double phi = 2.0 * kPi * k / 8000.0;
      brute = std::min(brute, 1.0 + a * std::cos(phi) + b * std::cos(2.0 * phi));
    }
    CHECK_CLOSE_AT(density_min(a, b), brute, 0.0, 1e-6);
    CHECK(density_min(a, b) <= 1.0 + std::fabs(a) + std::fabs(b));
  }
}

TEST_CASE("production scenarios keep a healthy positivity margin") {
  const double s = s_mid_6li();
  InclusiveKernel::Options opt;
  opt.b1_func = toy_b1_func();
  opt.delta_func = [](double, double, double f1) { return -1e-2 * f1; };
  const InclusiveKernel kern(LI6(), opt);
  const EventSpinState st{0, 0.6, 1.0, 1.0, kPi / 2.0, 0.0};
  for (const Point& p : analysis_grid(s)) {
    const Amplitudes a = kern.amplitudes(kern.tables(p.x, p.q2), p.x, p.q2, s, st);
    CHECK(InclusiveKernel::positivity_margin(a) > 0.9);
  }
  // ... while an oversized Delta drives the density negative, which is what
  // the sampler's guard exists to catch
  InclusiveKernel::Options big_opt;
  big_opt.delta_func = [](double, double, double f1) { return 3.0 * f1; };
  const InclusiveKernel big(LI6(), big_opt);
  bool negative = false;
  for (const Point& p : analysis_grid(s)) {
    const Amplitudes a = big.amplitudes(big.tables(p.x, p.q2), p.x, p.q2, s, st);
    if (InclusiveKernel::positivity_margin(a) < 0.0) negative = true;
  }
  CHECK(negative);
}

TEST_CASE("dsigma is the unpolarized cross section times W(phi')/2pi") {
  const double s = s_mid_6li();
  InclusiveKernel::Options opt;
  opt.b1_func = toy_b1_func();
  opt.delta_func = toy_delta_func();
  const InclusiveKernel kern(LI6(), opt);
  const double x = 0.056, q2 = 1.14;
  const EventSpinState st{+1, 0.7, 1.0, 1.0, 0.9, 0.3};
  const SFTables t = kern.tables(x, q2, true);
  const Amplitudes a = kern.amplitudes(t, x, q2, s, st, true);
  for (int i = 0; i < 8; ++i) {
    const double phi = 2.0 * kPi * i / 8.0;
    const double w = InclusiveKernel::density(a, phi - st.phi_s);
    CHECK_CLOSE(kern.dsigma(x, q2, phi, s, st, true),
                kern.dsigma_unpol(x, q2, s) / (2.0 * kPi) * std::max(w, 0.0),
                kRtol);
  }
}

// --- the analytic estimator errors ----------------------------------------

TEST_CASE("the three statistical error formulas") {
  // These are the spreads the P3 estimator closure has to reproduce to 15 %:
  // err(A_par) = 1/(P_e P_z sqrt(N)) for two-state flips with equal
  // luminosity halves, err(A_zz) = sqrt(2/N)/P_zz for the thirds estimator,
  // and the same for a cos(2phi) fit amplitude.
  const double n = 1.0e4, pe = 0.7, pz = 0.6, pzz = 0.8;
  CHECK_CLOSE(err_a_parallel(n, pe, pz), 1.0 / (pe * pz * std::sqrt(n)), kRtol);
  CHECK_CLOSE(err_azz(n, pzz), std::sqrt(2.0 / n) / pzz, kRtol);
  CHECK_CLOSE(err_cos2phi_amplitude(n, pzz), std::sqrt(2.0 / n) / pzz, kRtol);
  // they all fall as 1/sqrt(N)
  CHECK_CLOSE(err_a_parallel(4.0 * n, pe, pz), 0.5 * err_a_parallel(n, pe, pz),
              kRtol);
  CHECK_CLOSE(err_azz(4.0 * n, pzz), 0.5 * err_azz(n, pzz), kRtol);
  // an empty bin is floored rather than dividing by zero
  CHECK(std::isfinite(err_azz(0.0, pzz)));
}

TEST_CASE("phi_averaged_density is the kernel's own D_phi") {
  const double s = s_mid_6li();
  const InclusiveKernel kern(LI6());
  for (const Point& p : analysis_grid(s)) {
    const SFTables t = kern.tables(p.x, p.q2);
    CHECK(InclusiveKernel::dphi(t, p.x, p.y)
          == phi_averaged_density(t.f1, t.f2, p.x, p.y));
    // ... and the tensor kernel K is the same combination on (b1, b2)
    CHECK(InclusiveKernel::tensor_kernel(t, p.x, p.y)
          == phi_averaged_density(t.b1, t.b2, p.x, p.y));
  }
}
