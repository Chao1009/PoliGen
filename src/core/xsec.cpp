#include "lipolgen/xsec.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace lipolgen {

std::pair<double, double> theta_q_cos_sin(double y, double gamma2) {
  const double root = std::sqrt(1.0 + gamma2);
  const double cos = (1.0 + 0.5 * gamma2 * y) / root;
  const double sin =
      std::sqrt(gamma2 * std::max(1.0 - y - 0.25 * gamma2 * y * y, 0.0)) / root;
  return std::make_pair(cos, sin);
}

CosynTensorSFs cosyn_tensor_sfs(double b1, double b2, double b3, double b4,
                                double x, double gamma2) {
  const double g = std::sqrt(gamma2);
  const double onep = 1.0 + gamma2;
  CosynTensorSFs o;
  // Eq. (17a): the leading 2 multiplies b1 ALONE -- the large bracket opens
  // before it.  (Applying it to the whole bracket gives 4 gamma^2 for the
  // 5 gamma^2 of Table 1's second row and misses it by up to 2.5 %.)
  o.f_t = -(2.0 * onep * b1 - (gamma2 / x) * (b2 / 6.0 - b3 / 2.0));
  o.f_l = (2.0 * onep * x * b1
           - onep * onep * (b2 / 3.0 + b3 + b4)
           - onep * (b2 / 3.0 - b4)
           - (b2 / 3.0 - b3)) / x;
  o.f_lt = -(g / (2.0 * x)) * (onep * (b2 / 3.0 - b4)
                               + (2.0 / 3.0 * b2 - 2.0 * b3));
  o.f_tt = -(gamma2 / x) * (b2 / 6.0 - b3 / 2.0);
  return o;
}

std::pair<double, double> cosyn_unpolarized_sfs(double f1, double f2, double x,
                                                double gamma2) {
  return std::make_pair(2.0 * f1, (1.0 + gamma2) * f2 / x - 2.0 * f1);
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
      b3_func_(std::move(options.b3_func)),
      b4_func_(std::move(options.b4_func)),
      g2_mode_(options.g2_mode),
      g2_scale_(options.g2_scale),
      target_mass_(options.target_mass),
      tensor_gamma_(options.tensor_gamma),
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
  // b3, b4 are filled for EVERY spin (xsec.py fills them outside the rank-2
  // branch) and read only by the finite-gamma tensor path.
  t.b3 = b3_func_ ? b3_func_(x, q2, t.f1) : 0.0;
  t.b4 = b4_func_ ? b4_func_(x, q2, t.f1) : 0.0;

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
  return a_parallel_exact(t.g1, t.g2, t.f1, y, x, q2, r_func_);
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

TensorHarmonics InclusiveKernel::tensor_harmonics_gamma(
    const SFTables& t, double x, double q2, double y,
    const EventSpinState& state) const {
  const std::pair<double, double> qc = tensor_moments(state.m);
  const double pref = 1.5 * qc.first;
  const double g2v = gamma_squared(x, q2);
  const double eps = epsilon_gamma(y, g2v);
  const std::pair<double, double> cs = theta_q_cos_sin(y, g2v);
  const double c = cs.first, sn = cs.second;
  const double ct = std::cos(state.theta_s), st = std::sin(state.theta_s);
  const CosynTensorSFs f = cosyn_tensor_sfs(t.b1, t.b2, t.b3, t.b4, x, g2v);
  const std::pair<double, double> fu = cosyn_unpolarized_sfs(t.f1, t.f2, x, g2v);
  const double den = fu.first + eps * fu.second;
  // the structure-function combination each alignment channel meets
  const double lam_ll = f.f_t + eps * f.f_l;
  const double lam_lt = std::sqrt(2.0 * eps * (1.0 + eps)) * f.f_lt;
  const double lam_tt = eps * f.f_tt;
  // harmonics of t_zz, t_xz, t_xx - t_yy in cos(n phi')
  const double zz[3] = {
      pref * (0.5 * sn * sn * st * st + c * c * ct * ct - 1.0 / 3.0),
      -pref * 2.0 * sn * c * st * ct,
      pref * 0.5 * sn * sn * st * st};
  const double xz[3] = {
      pref * c * sn * (ct * ct - 0.5 * st * st),
      pref * (c * c - sn * sn) * st * ct,
      -pref * 0.5 * c * sn * st * st};
  const double tt[3] = {
      pref * sn * sn * (ct * ct - 0.5 * st * st),
      pref * 2.0 * c * sn * st * ct,
      pref * 0.5 * st * st * (c * c + 1.0)};
  const double scale = -TENSOR_LL_SIGN / std::max(den, 1e-30);
  TensorHarmonics h;
  double* out[3] = {&h.h0, &h.h1, &h.h2};
  for (int n = 0; n < 3; ++n) {
    *out[n] = scale * (zz[n] * lam_ll + xz[n] * lam_lt + tt[n] * lam_tt);
  }
  return h;
}

Amplitudes InclusiveKernel::tensor_amplitudes(const SFTables& t, double x,
                                              double q2, double s,
                                              const EventSpinState& state,
                                              bool with_delta) const {
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
  //
  // The check lives HERE and not in `amplitudes()` because `amplitudes()` is
  // now defined as this function plus the vector terms and calls it first and
  // unconditionally, so both entry points refuse the same mismatched pair.
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
  double a1_tensor = 0.0;
  if (j >= 1.0 - 1e-9 && j_ion >= 1.0 - 1e-9) {
    const std::pair<double, double> qc = tensor_moments(state.m);
    if (tensor_gamma_) {
      const TensorHarmonics h = tensor_harmonics_gamma(t, x, q2, y, state);
      out.w_avg = out.w_avg + h.h0;
      a1_tensor = a1_tensor + h.h1;
      out.a2 = out.a2 + h.h2;
    } else {
      // T_LL = Q_NN P_2(cos theta_S), one line for every spin
      const double t_geo =
          TENSOR_LL_SIGN * qc.first * 0.5 * (3.0 * ct * ct - 1.0);
      const double kern = tensor_kernel(t, x, y);
      out.w_avg = out.w_avg + t_geo * kern / std::max(den, 1e-30);
    }
    // Delta (gluon transversity) is rank 2 as well but is NOT in POLRAD's
    // b1..b4 basis and nobody has computed RC for a phi-dependent tensor
    // observable, so `rc.hpp` leaves it out of the band by default
    // (RcScope::TensorRate) and `with_delta` is what buys it back
    // (RcScope::TensorAll, for pricing the omission).  `amplitudes()` always
    // wants it -- it is part of W.
    if (with_delta) {
      out.a2 = out.a2 + (-(1.0 - y) / (y * y) * qc.second * st * st * t.delta
                         / std::max(den, 1e-30));
    }
  }
  out.a1 = a1_tensor;
  return out;
}

Amplitudes InclusiveKernel::amplitudes(const SFTables& t, double x, double q2,
                                       double s, const EventSpinState& state,
                                       bool with_perp) const {
  // The whole b-sector, in ONE place (`tensor_amplitudes`, which `rc.hpp` and
  // `InclusiveSampler::StateTables::w_tensor` read as well), plus the vector
  // terms.  A pure extraction: it moves no number, which the rtol-1e-12
  // `validation/reference/*.json` gates of tests/test_reference.cpp prove.
  Amplitudes out = tensor_amplitudes(t, x, q2, s, state, /*with_delta=*/true);

  const double y = q2 / (s * x);
  const double j = state.j;
  const double ct = std::cos(state.theta_s);
  const double st = std::sin(state.theta_s);
  const double helicity = state.lam_e * state.pe;
  const double v = (j > 0.0) ? state.m / j : 0.0;
  if (helicity != 0.0 && v != 0.0) {
    out.w_avg = out.w_avg + helicity * v * ct * a_parallel(t, x, q2, y);
  }
  if (with_perp && helicity != 0.0 && v != 0.0 && std::fabs(st) > 1e-12) {
    out.a1 = out.a1 + helicity * v * st * a_perp(t, x, q2, y);
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
