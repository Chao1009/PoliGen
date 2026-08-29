#include "lipolgen/xsec.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace lipolgen {

double gamma_squared(double x, double q2, double m) {
  return 4.0 * m * m * x * x / q2;
}

double epsilon_gamma(double y, double gamma2) {
  return ((1.0 - y - 0.25 * gamma2 * y * y)
          / (1.0 - y + 0.5 * y * y + 0.25 * gamma2 * y * y));
}

double depolarization_gamma(double y, double gamma2, double r) {
  const double eps = epsilon_gamma(y, gamma2);
  return (1.0 - (1.0 - y) * eps) / (1.0 + eps * r);
}

double eta_gamma(double y, double gamma2) {
  const double eps = epsilon_gamma(y, gamma2);
  return eps * std::sqrt(gamma2) * y / (1.0 - (1.0 - y) * eps);
}

double density_min(double a1n, double a2n) {
  const double a = a1n;
  const double b = a2n;
  const double ends = 1.0 - std::fabs(a) + b;  // min(f(+1), f(-1))
  if (b > 0.0) {
    const double cstar = -a / (4.0 * b);
    if (std::fabs(cstar) <= 1.0) return (1.0 - b) - a * a / (8.0 * b);
  }
  return ends;
}

InclusiveKernel::InclusiveKernel(Ion ion)
    : InclusiveKernel(std::move(ion), Options{}) {}

InclusiveKernel::InclusiveKernel(Ion ion, Options options)
    : ion_(std::move(ion)),
      r_func_(options.r_func),
      nf2_(ion_,
           options.f2_source ? options.f2_source
                             : std::static_pointer_cast<const UnpolSF>(
                                   std::make_shared<const ToyF2>()),
           options.emc_ratio, options.r_func),
      g1_model_(options.g1_model),
      b1_func_(std::move(options.b1_func)),
      b2_func_(std::move(options.b2_func)),
      delta_func_(std::move(options.delta_func)),
      b1_32_func_(std::move(options.b1_32_func)),
      b2_32_func_(std::move(options.b2_32_func)),
      delta_32_func_(std::move(options.delta_32_func)),
      g2_mode_(options.g2_mode),
      g2_scale_(options.g2_scale),
      target_mass_(options.target_mass),
      g2_npts_(options.g2_npts) {
  if (!g1_model_) {
    g1_model_ = std::make_shared<const ToyG1>(nf2_.base(), options.r_func);
  }
  // `target_mass = true` with `G2Mode::kZero` is NOT refused (xsec.py:126-129):
  // g2 is filled with zeros and the finite-gamma kernel runs on them, which is
  // exactly the g2_scale = 0 twist-3 variation.
}

double InclusiveKernel::g1a(double x, double q2) const {
  return g1_model_->g1_nucleus(ion_, x, q2) / static_cast<double>(ion_.A);
}

SFTables InclusiveKernel::tables(double x, double q2, bool with_g2) const {
  SFTables t;
  t.f2 = nf2_.f2a(x, q2) / static_cast<double>(ion_.A);
  t.f1 = nf2_.f1a(x, q2) / static_cast<double>(ion_.A);
  t.g1 = g1a(x, q2);

  const double j = ion_.spin;
  const SFFunc3* b1f = nullptr;
  const SFFunc3* b2f = nullptr;
  const SFFunc3* df = nullptr;
  bool rank2 = false;
  if (std::fabs(j - 1.0) < 1e-9) {
    b1f = &b1_func_;
    b2f = &b2_func_;
    df = &delta_func_;
    rank2 = true;
  } else if (std::fabs(j - 1.5) < 1e-9) {
    b1f = &b1_32_func_;
    b2f = &b2_32_func_;
    df = &delta_32_func_;
    rank2 = true;
  }
  if (rank2) {
    t.b1 = (*b1f) ? (*b1f)(x, q2, t.f1) : 0.0;
    t.b2 = (*b2f) ? (*b2f)(x, q2, t.f1) : 2.0 * x * t.b1;
    t.delta = (*df) ? (*df)(x, q2, t.f1) : 0.0;
  }

  if (with_g2 || target_mass_) {
    t.has_g2 = true;
    t.g2 = (g2_mode_ == G2Mode::kWandzuraWilczek)
               ? g2_scale_ * g2_ww([this](double xx, double qq) {
                   return g1a(xx, qq);
                 }, x, q2, g2_npts_)
               : 0.0;
  }
  return t;
}

double InclusiveKernel::dphi(const SFTables& t, double x, double y) {
  return t.f1 + (1.0 - y) / (x * y * y) * t.f2;
}

double InclusiveKernel::tensor_kernel(const SFTables& t, double x, double y) {
  return t.b1 + (1.0 - y) / (x * y * y) * t.b2;
}

double InclusiveKernel::a_parallel(const SFTables& t, double x, double q2,
                                   double y) const {
  if (!target_mass_) {
    return depolarization_d(y, x, q2, r_func_) * t.g1 / std::max(t.f1, 1e-30);
  }
  if (!t.has_g2) {
    throw std::runtime_error(
        "target_mass=true needs g2 in tables(); build them with this kernel");
  }
  const double g2v = gamma_squared(x, q2);
  const double r = resolve_r(r_func_, x, q2);
  const double f1 = std::max(t.f1, 1e-30);
  const double a1 = (t.g1 - g2v * t.g2) / f1;
  const double a2 = std::sqrt(g2v) * (t.g1 + t.g2) / f1;
  return depolarization_gamma(y, g2v, r) * (a1 + eta_gamma(y, g2v) * a2);
}

double InclusiveKernel::a_perp(const SFTables& t, double x, double q2,
                               double y) const {
  if (!t.has_g2) {
    throw std::runtime_error("tables(..., with_g2=true) required for a_perp");
  }
  const double eps = (1.0 - y) / (1.0 - y + 0.5 * y * y);
  const double d = depolarization_d(y, x, q2, r_func_)
                   * std::sqrt(std::max(2.0 * eps / (1.0 + eps), 0.0));
  const double gamma = 2.0 * M_NUCLEON * x / std::sqrt(q2);
  const double amp = gamma * (0.5 * y * t.g1 + t.g2) / std::max(t.f1, 1e-30);
  return d * amp;
}

std::pair<double, double> InclusiveKernel::tensor_moments(double m) const {
  const double j = ion_.spin;
  if (j < 1.0 - 1e-9) return std::make_pair(0.0, 0.0);
  const double q_nn = (3.0 * m * m - j * (j + 1.0)) / 3.0;
  return std::make_pair(q_nn, 3.0 * q_nn);
}

Amplitudes InclusiveKernel::amplitudes(const SFTables& t, double x, double q2,
                                       double s, const EventSpinState& state,
                                       bool with_perp) const {
  // P8.  The rank-2 branch below used to be gated on `state.j` alone while
  // `tensor_moments(state.m)` reads the KERNEL's ion spin, so a mismatched
  // pair would compute the alignment Q_NN for one spin and apply it to a
  // population of another.  Two rules close that:
  //
  //   (a) a kernel that HAS a rank-2 sector (ion spin >= 1, the only case in
  //       which `tables()` fills b1/b2/Delta at all) refuses any spin state
  //       that is not its own J;
  //   (b) the branch is entered only when BOTH spins are >= 1, so the gate and
  //       the moments can no longer disagree.
  //
  // (b) changes no number: with ion spin < 1 `tensor_moments` returns (0, 0)
  // and `tables()` leaves b1 = b2 = Delta = 0, so the branch contributed
  // exactly zero.  That is the d(e,e'p) control, where the struck cluster is a
  // spin-1/2 NEUTRON but the projection m_S labels the S_c = 1 channel spin --
  // legitimate, and what the Python does.
  const double j_ion = ion_.spin;
  if (j_ion >= 1.0 - 1e-9 && std::fabs(state.j - j_ion) > 1e-9) {
    throw std::runtime_error(
        "InclusiveKernel::amplitudes: spin state J = " +
        std::to_string(state.j) + " is not the kernel's ion spin " +
        std::to_string(j_ion));
  }
  const double y = q2 / (s * x);
  const double den = dphi(t, x, y);
  const double j = state.j;
  const double ct = std::cos(state.theta_s);
  const double st = std::sin(state.theta_s);

  Amplitudes out;
  if (j >= 1.0 - 1e-9 && j_ion >= 1.0 - 1e-9) {
    const std::pair<double, double> qc = tensor_moments(state.m);
    // T_LL = Q_NN P_2(cos theta_S), one line for every spin
    const double t_geo = TENSOR_LL_SIGN * qc.first * 0.5 * (3.0 * ct * ct - 1.0);
    const double kern = tensor_kernel(t, x, y);
    out.w_avg = out.w_avg + t_geo * kern / std::max(den, 1e-30);
    out.a2 = -(1.0 - y) / (y * y) * qc.second * st * st * t.delta
             / std::max(den, 1e-30);
  }

  const double helicity = state.lam_e * state.pe;
  const double v = (j > 0.0) ? state.m / j : 0.0;
  if (helicity != 0.0 && v != 0.0) {
    out.w_avg = out.w_avg + helicity * v * ct * a_parallel(t, x, q2, y);
  }
  if (with_perp && helicity != 0.0 && v != 0.0 && std::fabs(st) > 1e-12) {
    out.a1 = helicity * v * st * a_perp(t, x, q2, y);
  }
  return out;
}

double InclusiveKernel::dsigma_unpol(double x, double q2, double s) const {
  const double f2 = nf2_.f2a(x, q2) / static_cast<double>(ion_.A);
  return dsigma_dx_dq2(x, q2, s, f2, nullptr, r_func_);
}

double InclusiveKernel::density(const Amplitudes& a, double phip) {
  return 1.0 + a.w_avg + a.a1 * std::cos(phip) + a.a2 * std::cos(2.0 * phip);
}

double InclusiveKernel::positivity_margin(const Amplitudes& a) {
  const double den = 1.0 + a.w_avg;
  return density_min(a.a1 / den, a.a2 / den);
}

double InclusiveKernel::dsigma(double x, double q2, double phi, double s,
                               const EventSpinState& state,
                               bool with_perp) const {
  const SFTables t = tables(x, q2, with_perp);
  const Amplitudes a = amplitudes(t, x, q2, s, state, with_perp);
  const double phip = phi - state.phi_s;
  const double w = density(a, phip);
  return dsigma_unpol(x, q2, s) / (2.0 * kPi) * std::max(w, 0.0);
}

}  // namespace lipolgen
