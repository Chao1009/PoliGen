#include "lipolgen/sf.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "lipolgen/numerics.hpp"

namespace lipolgen {
namespace {

// R1998: Table II of Abe et al., PLB 452 (1999) 194 (arXiv:hep-ex/9808028).
const double kR1998A[6] = {0.0485, 0.5470, 2.0621, -0.3804, 0.5090, -0.0285};
const double kR1998B[6] = {0.0481, 0.6114, -0.3509, -0.4611, 0.7172, -0.0317};
const double kR1998C[6] = {0.0577, 0.4644, 1.8288, 12.3708, -43.1043, 41.7415};

// The leading term of every form is 1/log(Q2/0.04), whose pole sits at
// Q2 = Lambda^2 = 0.04 GeV2; with clip=False the Q2 floor is still enforced.
constexpr double kR1998Q2PoleFloor = 0.05;

double r1998_theta(double x, double q2) {
  return 1.0 + 12.0 * (q2 / (q2 + 1.0))
                   * (0.125 * 0.125 / (0.125 * 0.125 + x * x));
}

void r1998_prepare(double& x, double& q2, bool do_clip) {
  if (do_clip) {
    x = clip(x, R1998_X_MIN, R1998_X_MAX);
    q2 = clip(q2, R1998_Q2_MIN, R1998_Q2_MAX);
  } else {
    q2 = std::max(q2, kR1998Q2PoleFloor);
    x = std::max(x, 1e-12);  // x^a6 with a6 < 0 diverges at x = 0
  }
}

// The 12-point canonical unpolarized EMC shape (SCENARIO).
const double kEmcX[12] = {1e-4, 0.01, 0.06, 0.10, 0.20, 0.30,
                          0.45, 0.60, 0.70, 0.80, 0.88, 0.95};
const double kEmcR[12] = {0.96, 0.98, 1.00, 1.01, 1.00, 0.98,
                          0.95, 0.91, 0.88, 0.90, 1.00, 1.15};

double mean_of(const std::vector<double>& v) {
  return pairwise_sum(v) / static_cast<double>(v.size());
}

}  // namespace

// ---------------------------------------------------------------------- R

double r_sigma_lt(double, double q2) { return 0.18 / (1.0 + q2 / 50.0); }

R1998Forms r1998_forms(double x, double q2, bool do_clip) {
  r1998_prepare(x, q2, do_clip);
  const double rlog = r1998_theta(x, q2) / std::log(q2 / 0.04);

  const double a1 = kR1998A[0], a2 = kR1998A[1], a3 = kR1998A[2];
  const double a4 = kR1998A[3], a5 = kR1998A[4], a6 = kR1998A[5];
  const double r_a = a1 * rlog + (a2 * std::pow(std::pow(q2, 4.0) + std::pow(a3, 4.0), -0.25)
                                  * (1.0 + a4 * x + a5 * x * x)
                                  * std::pow(x, a6));

  const double b1c = kR1998B[0], b2 = kR1998B[1], b3 = kR1998B[2];
  const double b4 = kR1998B[3], b5 = kR1998B[4], b6 = kR1998B[5];
  const double r_b = b1c * rlog + ((b2 / q2 + b3 / (q2 * q2 + 0.3 * 0.3))
                                   * (1.0 + b4 * x + b5 * x * x)
                                   * std::pow(x, b6));

  const double c1 = kR1998C[0], c2 = kR1998C[1], c3 = kR1998C[2];
  const double c4 = kR1998C[3], c5 = kR1998C[4], c6 = kR1998C[5];
  const double q2thr = c4 * x + c5 * x * x + c6 * std::pow(x, 3.0);
  const double r_c = c1 * rlog + c2 * std::pow((q2 - q2thr) * (q2 - q2thr)
                                               + c3 * c3, -0.5);

  return R1998Forms{r_a, r_b, r_c};
}

double r1998(double x, double q2, R1998Form form, bool do_clip) {
  const R1998Forms f = r1998_forms(x, q2, do_clip);
  switch (form) {
    case R1998Form::kAverage: return (f.a + f.b + f.c) / 3.0;
    case R1998Form::kA: return f.a;
    case R1998Form::kB: return f.b;
    case R1998Form::kC: return f.c;
  }
  throw std::runtime_error("unknown R1998 form");
}

double r1998_spread(double x, double q2, bool do_clip) {
  const R1998Forms f = r1998_forms(x, q2, do_clip);
  return std::max(std::max(f.a, f.b), f.c) - std::min(std::min(f.a, f.b), f.c);
}

double r1998_fit_error(double x, double q2) {
  return (0.0078 - 0.013 * x
          + (0.070 - 0.39 * x + 0.70 * x * x) / (1.7 + q2));
}

// ------------------------------------------------------------------- ToyF2

double ToyF2::f2p(double x, double q2) const {
  const double lq = std::log(std::max(q2, 1.1) / 0.04);  // Lambda = 0.2
  const double lam = 0.045 * lq;
  const double sea = 0.20 * std::pow(x, -lam) * std::pow(1.0 - x, 7.0);
  const double val = 1.05 * std::pow(x, 0.55) * std::pow(1.0 - x, 3.0);
  return sea + val;
}

double ToyF2::f2n_over_f2p(double x) const {
  return clip(1.0 - 0.75 * x, 0.25, 1.0);
}

double ToyF2::f2n(double x, double q2) const {
  return f2p(x, q2) * f2n_over_f2p(x);
}

// --------------------------------------------------------------- NuclearF2

NuclearF2::NuclearF2(Ion ion, std::shared_ptr<const UnpolSF> base,
                     std::function<double(double)> emc_ratio, RFunc r_func)
    : ion_(std::move(ion)),
      base_(base ? std::move(base) : std::make_shared<const ToyF2>()),
      emc_ratio_(std::move(emc_ratio)),
      r_func_(std::move(r_func)) {}

double NuclearF2::f2a(double x, double q2) const {
  double f2 = ion_.Z * base_->f2p(x, q2) + ion_.N() * base_->f2n(x, q2);
  if (emc_ratio_) f2 = f2 * emc_ratio_(x);
  return f2;
}

double NuclearF2::f1a(double x, double q2) const {
  const double r = resolve_r(r_func_, x, q2);
  return f2a(x, q2) / (2.0 * x * (1.0 + r));
}

double dsigma_dx_dq2(double x, double q2, double s, double f2, const double* fl,
                     const RFunc& r_func) {
  const double y = q2 / (s * x);
  double fl_val;
  if (fl == nullptr) {
    const double r = resolve_r(r_func, x, q2);
    fl_val = f2 * r / (1.0 + r);
  } else {
    fl_val = *fl;
  }
  const double bracket = (1.0 - y + 0.5 * y * y) * f2 - 0.5 * y * y * fl_val;
  const double xsec = 4.0 * kPi * ALPHA_EM * ALPHA_EM / (x * q2 * q2) * bracket;
  return std::max(xsec, 0.0) * GEV2_TO_PB;
}

// -------------------------------------------------------------------- g2 WW

double g2_ww(const std::function<double(double, double)>& g1, double x,
             double q2, int npts) {
  const std::vector<double> t = linspace(0.0, 1.0, static_cast<std::size_t>(npts));
  std::vector<double> g1u(t.size());
  for (std::size_t i = 0; i < t.size(); ++i) {
    g1u[i] = g1(std::pow(x, 1.0 - t[i]), q2);
  }
  const double integral = -std::log(std::max(x, 1e-12)) * trapezoid(g1u, t);
  return -g1(x, q2) + integral;
}

// -------------------------------------------------------------------- PolSF

double PolSF::g1_nucleus(const Ion& ion, double x, double q2,
                         const std::function<double(double)>& medium_ratio) const {
  double g1 = ion.Z * ion.eff_pol_p * g1p(x, q2)
              + ion.N() * ion.eff_pol_n * g1n(x, q2);
  if (medium_ratio) g1 = g1 * medium_ratio(x);
  return g1;
}

double PolSF::g2p(double x, double q2, int npts) const {
  return g2_ww([this](double xx, double qq) { return g1p(xx, qq); }, x, q2, npts);
}

double PolSF::g2n(double x, double q2, int npts) const {
  return g2_ww([this](double xx, double qq) { return g1n(xx, qq); }, x, q2, npts);
}

ToyG1::ToyG1(std::shared_ptr<const UnpolSF> base, RFunc r_func)
    : base_(base ? std::move(base) : std::make_shared<const ToyF2>()),
      r_func_(std::move(r_func)) {}

double ToyG1::a1p(double x) const {
  // x^0.7 tracks moderate/high-x world data; saturates below 1
  return clip(std::pow(x, 0.7), 0.0, 1.0);
}

double ToyG1::a1n(double x) const {
  // small negative at low/mid x, positive rise at high x
  return -0.07 * std::pow(1.0 - x, 2.0) + 0.8 * std::pow(x, 2.2);
}

double ToyG1::f1_of(double f2, double x, double q2) const {
  const double r = resolve_r(r_func_, x, q2);
  return f2 / (2.0 * x * (1.0 + r));
}

double ToyG1::g1p(double x, double q2) const {
  return a1p(x) * f1_of(base_->f2p(x, q2), x, q2);
}

double ToyG1::g1n(double x, double q2) const {
  return a1n(x) * f1_of(base_->f2n(x, q2), x, q2);
}

// --------------------------------------------------------- DigitizedTable

const double* DigitizedTable::column(const std::string& name) const {
  for (std::size_t i = 0; i < columns.size(); ++i) {
    if (columns[i] == name) return values[i];
  }
  throw std::runtime_error("digitized table has no column " + name);
}

double DigitizedTable::interp(const std::string& name, double xv) const {
  return np_interp(xv, x, column(name), n);
}

double DigitizedTable::interp_tapered(const std::string& name, double xv,
                                      double x_end, double power) const {
  const double* col = column(name);
  const double xmax = x_max();
  double val = np_interp(xv, x, col, n);
  if (xmax < x_end && xv > xmax) {
    const double frac = clip((x_end - xv) / (x_end - xmax), 0.0, 1.0);
    val = val * std::pow(frac, power);
  }
  return val;
}

// ------------------------------------------------------------------ EMC

double unpolarized_emc_ratio(double x) {
  return np_interp(x, kEmcX, kEmcR, 12);
}

double cbt_unpolarized_emc_ratio(double x) {
  return tables::kCbtPolemc7liQ5().interp("R_unpol", x);
}

double cbt_valence_scale() { return 1.0; }

double tmt_valence_scale() {
  static const double s = []() {
    const std::vector<double> g = linspace(POLEMC_VALENCE_WINDOW_LO,
                                           POLEMC_VALENCE_WINDOW_HI, 301);
    std::vector<double> d_ref(g.size()), d_mod(g.size());
    for (std::size_t i = 0; i < g.size(); ++i) {
      d_ref[i] = 1.0 - tables::kCbtPolemc7liQ5().interp("R_unpol", g[i]);
      d_mod[i] = 1.0 - tables::kTmtPolemcNmQ10().interp("R_unpol", g[i]);
    }
    return mean_of(d_ref) / mean_of(d_mod);
  }();
  return s;
}

double cbt_polarized_emc_ratio(double x, EmcMode mode, int eq) {
  if (mode == EmcMode::kConstant) {
    return 1.0 - 2.0 * (1.0 - unpolarized_emc_ratio(x));
  }
  if (eq != 23 && eq != 26) {
    throw std::runtime_error("eq must be 23 (R^{3/2 3/2}) or 26 (R^{(3/2 1)})");
  }
  return tables::kCbtPolemc7liQ5().interp(
      eq == 23 ? "R_pol_eq23" : "R_pol_eq26", x);
}

double tmt_published_emc_ratio(double x) {
  return tables::kTmtPolemcNmQ10().interp("R_pol", x);
}

double tmt_polarized_emc_ratio(double x, EmcMode mode) {
  if (mode == EmcMode::kConstant) return unpolarized_emc_ratio(x);
  return 1.0 - tmt_valence_scale() * (1.0 - tmt_published_emc_ratio(x));
}

double cbt_ratio_of_effects(double x, int eq) {
  const double d_unpol = 1.0 - cbt_unpolarized_emc_ratio(x);
  return (1.0 - cbt_polarized_emc_ratio(x, EmcMode::kDigitized, eq)) / d_unpol;
}

double tmt_ratio_of_effects(double x) {
  return (1.0 - tmt_published_emc_ratio(x))
         / (1.0 - tables::kTmtPolemcNmQ10().interp("R_unpol", x));
}

// ---------------------------------------------------------------- tensor

SFFunc3 TensorSF::b1_func() const {
  return [this](double x, double q2, double f1) { return b1(x, q2, f1); };
}
SFFunc3 TensorSF::b2_func() const {
  return [this](double x, double q2, double f1) { return b2(x, q2, f1); };
}
SFFunc3 TensorSF::delta_func() const {
  return [this](double x, double q2, double f1) { return delta(x, q2, f1); };
}

double toy_b1_shape(double x, double, double f1) {
  const double shape = 0.01 * std::pow(std::max(x, 1e-6), -0.2)
                       * (1.0 - x / 0.20);
  return shape * std::exp(-3.0 * x) * f1;
}

double toy_b1(double x, double q2, double f1, B1Mode mode) {
  if (mode == B1Mode::kToy) return toy_b1_shape(x, q2, f1);
  return B1_PER_DEUTERON_TO_PER_NUCLEON
         * tables::kB1Miller().interp_tapered("b1", x);
}

double b1_convolution(double x, double q2, double f1, B1Mode mode) {
  if (mode == B1Mode::kToy) return 0.1 * toy_b1_shape(x, q2, f1);
  const double xs = std::max(x, 1e-6);
  const double xb1 = tables::kB1CdksQ2p5().interp("xb1_theory1_sum", xs);
  return B1_PER_DEUTERON_TO_PER_NUCLEON * xb1 / xs;
}

double close_kumano_integral(bool cdks) {
  const DigitizedTable& t = cdks ? tables::kB1CdksQ2p5() : tables::kB1Miller();
  const char* col = cdks ? "xb1_theory1_sum" : "b1";
  const double* v = t.column(col);
  std::vector<double> xs(t.x, t.x + t.n), ys(t.n);
  for (std::size_t i = 0; i < t.n; ++i) ys[i] = cdks ? v[i] / t.x[i] : v[i];
  return trapezoid(ys, xs);
}

double b1_li6_from_deuteron(double b1_d, double transfer, double per_nucleon) {
  return transfer * per_nucleon * b1_d;
}

double toy_delta_gluon(double x, double, double f1, double scale) {
  return scale * f1 * std::pow(x, 0.3) * std::pow(1.0 - x, 4.0);
}

// ---------------------------------------------------------- Delta models

DeltaVariant delta_variant(const std::string& name) {
  if (name == "low_x") return DeltaVariant{0.3, 4.0};
  if (name == "mid_x") return DeltaVariant{0.7, 3.0};
  if (name == "high_x") return DeltaVariant{1.5, 2.0};
  throw std::runtime_error("unknown Delta shape variant " + name);
}

double xab_peak_value(double alpha, double beta) {
  const double xp = alpha / (alpha + beta);
  return std::pow(xp, alpha) * std::pow(1.0 - xp, beta);
}

double shape_normalized(double x, const std::string& variant) {
  const DeltaVariant v = delta_variant(variant);
  return std::pow(std::max(x, 1e-12), v.alpha)
         * std::pow(std::max(1.0 - x, 0.0), v.beta)
         / xab_peak_value(v.alpha, v.beta);
}

double alpha_s_lo(double q2, double lambda_qcd, double n_f) {
  const double lam2 = lambda_qcd * lambda_qcd;
  const double ln = std::log(std::max(q2, lam2 * 1.01) / lam2);
  return 12.0 * kPi / ((33.0 - 2.0 * n_f) * ln);
}

namespace {

std::vector<double> moment_x_grid() {
  std::vector<double> g = logspace(-5.0, std::log10(0.5), 220);
  const std::vector<double> hi = linspace(0.5, 0.9999, 120);
  g.insert(g.end(), hi.begin(), hi.end());
  std::sort(g.begin(), g.end());
  g.erase(std::unique(g.begin(), g.end()), g.end());  // np.unique
  return g;
}

}  // namespace

double solve_A_interp_a(const std::function<double(double, double)>& f1_func,
                        double q2_ref, const std::string& variant,
                        double c_moment) {
  const std::vector<double> xg = moment_x_grid();
  std::vector<double> integrand(xg.size());
  for (std::size_t i = 0; i < xg.size(); ++i) {
    integrand[i] = xg[i] * f1_func(xg[i], q2_ref) * shape_normalized(xg[i], variant);
  }
  const double integral = trapezoid(integrand, xg);
  if (std::fabs(integral) < 1e-30) {
    throw std::runtime_error("sum-rule integral vanished");
  }
  return c_moment / integral;
}

double solve_A_interp_b(const std::string& variant, double c_moment) {
  const DeltaVariant v = delta_variant(variant);
  const double beta_int = std::tgamma(v.alpha + 2.0) * std::tgamma(v.beta + 1.0)
                          / std::tgamma(v.alpha + v.beta + 3.0);
  return c_moment * xab_peak_value(v.alpha, v.beta) / beta_int;
}

namespace {

std::string info_string(const char* name, const std::string& variant, double c,
                        double a_solved, const double* q2_ref, double dilution) {
  char buf[256];
  if (q2_ref != nullptr) {
    std::snprintf(buf, sizeof(buf),
                  "%s(variant=%s, c=%.4g, A=%.4g, q2_ref=%.4g, dilution=%.4g)",
                  name, variant.c_str(), c, a_solved, *q2_ref, dilution);
  } else {
    std::snprintf(buf, sizeof(buf),
                  "%s(variant=%s, c=%.4g, A=%.4g, dilution=%.4g)", name,
                  variant.c_str(), c, a_solved, dilution);
  }
  return std::string(buf);
}

}  // namespace

DeltaModel make_delta_toy(double scale) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "toy(scale=%.4g)", scale);
  return DeltaModel("toy",
                    [scale](double x, double q2, double f1) {
                      return toy_delta_gluon(x, q2, f1, scale);
                    },
                    std::string(buf));
}

DeltaModel make_moment_a(const std::function<double(double, double)>& f1_func,
                         double q2_ref, const std::string& variant,
                         double c_moment, double dilution,
                         std::function<double(double)> alphas) {
  if (!alphas) alphas = [](double q2) { return alpha_s_lo(q2); };
  const double a_solved = solve_A_interp_a(f1_func, q2_ref, variant, c_moment);
  const std::string var = variant;
  return DeltaModel("moment_A",
                    [dilution, a_solved, alphas, var](double x, double q2,
                                                      double f1) {
                      return dilution * a_solved * alphas(q2) * f1
                             * shape_normalized(x, var);
                    },
                    info_string("moment_A", variant, c_moment, a_solved,
                                &q2_ref, dilution));
}

DeltaModel make_moment_b(const std::string& variant, double c_moment,
                         double dilution, std::function<double(double)> alphas) {
  if (!alphas) alphas = [](double q2) { return alpha_s_lo(q2); };
  const double a_solved = solve_A_interp_b(variant, c_moment);
  const std::string var = variant;
  return DeltaModel("moment_B",
                    [dilution, a_solved, alphas, var](double x, double q2,
                                                      double f1) {
                      (void)f1;
                      return dilution * a_solved * alphas(q2)
                             * shape_normalized(x, var);
                    },
                    info_string("moment_B", variant, c_moment, a_solved, nullptr,
                                dilution));
}

}  // namespace lipolgen
