// SPDX-License-Identifier: GPL-3.0-or-later
#include "lipolgen/rc.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#include "lipolgen/asymmetries.hpp"
#include "lipolgen/cluster.hpp"
#include "lipolgen/cluster_config.hpp"
#include "lipolgen/numerics.hpp"
#include "lipolgen/spin.hpp"

namespace lipolgen {
namespace {

// --------------------------------------------------------------- quadrature

/// 8-point Gauss-Legendre on [-1, 1].  The t-peak integrand lives in the few
/// e-folds of ln(eta_A) where the form factor is alive while the kinematic
/// range spans fifteen, so ONE high-order rule over the whole range would put
/// most of its nodes where the integrand is zero: `gl_log_quadrature` below
/// applies this rule PANEL BY PANEL instead.
constexpr int kGlOrder = 8;
constexpr double kGlNode[kGlOrder] = {
    -0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
    -0.1834346424956498,  0.1834346424956498,  0.5255324099163290,
     0.7966664774136267,  0.9602898564975363};
constexpr double kGlWeight[kGlOrder] = {
    0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
    0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
    0.2223810344533745, 0.1012285362903763};

/// INT_lo^hi f(v) dv/v  ==  INT_{ln lo}^{ln hi} f(e^u) du, by composite
/// Gauss-Legendre with `n` nodes in panels of 8.
template <class F>
double gl_log_quadrature(const F& f, double lo, double hi, int n) {
  if (!(hi > lo) || !(lo > 0.0)) return 0.0;
  const int n_panel = std::max(1, n / kGlOrder);
  const double a = std::log(lo);
  const double b = std::log(hi);
  const double h = (b - a) / static_cast<double>(n_panel);
  double acc = 0.0;
  for (int p = 0; p < n_panel; ++p) {
    const double c0 = a + h * static_cast<double>(p);
    const double mid = c0 + 0.5 * h;
    for (int i = 0; i < kGlOrder; ++i) {
      const double u = mid + 0.5 * h * kGlNode[i];
      acc += kGlWeight[i] * f(std::exp(u));
    }
  }
  return 0.5 * h * acc;
}

bool is_tagged_channel(Channel c) {
  switch (c) {
    case Channel::TaggedLi6Alpha:
    case Channel::TaggedLi6D:
    case Channel::TaggedLi7Alpha:
    case Channel::TaggedLi7T:
    case Channel::TaggedDeuteronP:
    case Channel::TaggedDeuteronN:
    case Channel::TaggedHe3P:
      return true;
    case Channel::Inclusive:
    case Channel::CoherentLi6:
      return false;
  }
  throw std::runtime_error("rc.cpp: unhandled Channel");
}

const char* channel_name(Channel c) {
  switch (c) {
    case Channel::Inclusive:       return "inclusive";
    case Channel::TaggedLi6Alpha:  return "tagged-6Li-alpha";
    case Channel::TaggedLi6D:      return "tagged-6Li-d";
    case Channel::TaggedLi7Alpha:  return "tagged-7Li-alpha";
    case Channel::TaggedLi7T:      return "tagged-7Li-t";
    case Channel::TaggedDeuteronP: return "tagged-d-p";
    case Channel::TaggedDeuteronN: return "tagged-d-n";
    case Channel::TaggedHe3P:      return "tagged-3He-p";
    case Channel::CoherentLi6:     return "coherent-6Li";
  }
  throw std::runtime_error("rc.cpp: unhandled Channel");
}

/// Bilinear index + fraction on a strictly increasing node vector, in ln v.
/// Throws rather than extrapolating: outside the table the tail is not a
/// small correction to anything the run generated.
void locate(const std::vector<double>& nodes, double v, const char* what,
            std::size_t* i0, double* frac) {
  if (nodes.size() < 2) {
    throw std::runtime_error(std::string("RcModel: the tail table has no ") +
                             what + " axis");
  }
  if (!(v >= nodes.front() * (1.0 - 1e-12)) ||
      !(v <= nodes.back() * (1.0 + 1e-12))) {
    throw std::runtime_error(
        "RcModel::tail_ratio_at: " + std::string(what) + " = " +
        std::to_string(v) + " is outside the tail table's support [" +
        std::to_string(nodes.front()) + ", " + std::to_string(nodes.back()) +
        "]; the tables refuse to extrapolate");
  }
  const std::size_t n = nodes.size();
  std::size_t k = static_cast<std::size_t>(
      std::upper_bound(nodes.begin(), nodes.end(), v) - nodes.begin());
  if (k == 0) k = 1;
  if (k > n - 1) k = n - 1;
  *i0 = k - 1;
  const double lo = std::log(nodes[k - 1]);
  const double hi = std::log(nodes[k]);
  *frac = (hi > lo) ? clip((std::log(v) - lo) / (hi - lo), 0.0, 1.0) : 0.0;
}

/// ln of a table entry, with the "this node is exactly zero" case pushed to
/// -inf-in-practice rather than to a NaN.
double safe_log(double v) { return std::log(std::max(v, 1e-300)); }

}  // namespace

// -------------------------------------------------------------- the band

double rc_delta(double x, double delta_high, double delta_low, double x_high,
                double x_low) {
  if (!(x_low > 0.0) || !(x_high > x_low)) {
    throw std::runtime_error(
        "rc_delta: the band anchors must satisfy 0 < x_low < x_high (got "
        "x_low = " + std::to_string(x_low) + ", x_high = " +
        std::to_string(x_high) + "); the interpolation divides by "
        "ln(x_high/x_low)");
  }
  // Both clamp branches return the anchor UNMODIFIED, which is what makes
  // `rc_delta(x_high) == delta_high` and `rc_delta(x_low) == delta_low` exact
  // (T3 compares them with ==, legitimately).
  if (x >= x_high) return delta_high;
  if (x <= x_low) return delta_low;
  const double f = std::log(x_high / x) / std::log(x_high / x_low);
  return delta_high + (delta_low - delta_high) * f;
}

bool rc_band_applies(Channel channel) {
  return channel != Channel::CoherentLi6;
}

bool rc_tail_applies(Channel channel, bool with_tail) {
  return rc_band_applies(channel) && !is_tagged_channel(channel) && with_tail;
}

const char* rc_mode_name(RcMode m) {
  switch (m) {
    case RcMode::Off:        return "off";
    case RcMode::TensorBand: return "tensor-band";
  }
  throw std::runtime_error("rc_mode_name: unhandled RcMode");
}

const char* rc_tail_model_name(RcTailModel m) {
  switch (m) {
    case RcTailModel::TPeak:       return "t-peak";
    case RcTailModel::PolradFull:  return "polrad-full";
    case RcTailModel::TPeakPlusLL: return "t-peak+ll";
  }
  throw std::runtime_error("rc_tail_model_name: unhandled RcTailModel");
}

const char* pipeline_rc_name(PipelineRc r) { return rc_mode_name(r); }

const char* c0_shape_name(C0Shape s) {
  switch (s) {
    case C0Shape::Ho:     return "ho";
    case C0Shape::VmcFt:  return "vmc-ft";
  }
  throw std::runtime_error("c0_shape_name: unknown C0Shape");
}

const char* rc_weight_name(std::size_t i) {
  switch (i) {
    case 0: return "rc_tensor_lo";
    case 1: return "rc_tensor_hi";
    case 2: return "rc_tail";
    default: break;
  }
  throw std::runtime_error("rc_weight_name: index " + std::to_string(i) +
                           " is outside the kRcWeightCount = 3 block");
}

std::string rc_weight_name(std::size_t i, std::size_t slot) {
  // Slot 0 is the event's OWN spin state and keeps the bare name, so an
  // unweighted run's HepMC3 weight names are exactly
  // {"nominal", "rc_tensor_lo", "rc_tensor_hi", "rc_tail"}.  Slot 1 + k is
  // spin category k and carries the "_k" suffix of "spin_weight_k".
  const std::string base = rc_weight_name(i);
  return slot == 0 ? base : base + "_" + std::to_string(slot);
}

// ---------------------------------------------------- nucleon form factors

NucleonFF nucleon_ff(double t_gev2) {
  const double t = std::max(t_gev2, 0.0);
  const double d = 1.0 + t / NUCLEON_FF_DIPOLE_GEV2;
  const double gd = 1.0 / (d * d);
  const double tau = t / (4.0 * PROTON_MASS * PROTON_MASS);
  NucleonFF f;
  f.ge_p = gd;
  f.gm_p = MU_PROTON * gd;
  f.gm_n = MU_NEUTRON * gd;
  // Galster: G_E^n = -mu_n tau G_D/(1 + a tau); mu_n < 0, so this is > 0.
  f.ge_n = -MU_NEUTRON * tau * gd / (1.0 + GALSTER_A * tau);
  return f;
}

// ------------------------------------------------ the eta_A quadratures

EtaLimits polrad_eta_limits(double x_a, double y, double s_a, double m_a) {
  EtaLimits lim;
  if (!(x_a > 0.0) || !(y > 0.0) || !(y < 1.0) || !(s_a > 0.0) ||
      !(m_a > 0.0)) {
    return lim;
  }
  const double q2 = x_a * y * s_a;
  const double sx = y * s_a;                      // S_x = S_A - X_A
  const double m2 = m_a * m_a;
  const double sqly = std::sqrt(sx * sx + 4.0 * m2 * q2);   // sqrt(lambda_Q)
  const double w2 = sx - q2 + m2;                 // W^2
  if (!(w2 > 0.0)) return lim;
  // adgh:8611, rewritten so that no digit cancels: the printed form has
  // (S_x - sqrt(lambda_Q)), a difference of two nearly equal large numbers,
  // and S_x - sqrt(lambda_Q) = -4 M_A^2 Q^2/(S_x + sqrt(lambda_Q)) exactly.
  lim.lo = q2 * (4.0 * m2 * q2 / (sx + sqly) + 2.0 * q2) /
           (4.0 * w2 * (sx + sqly));
  lim.hi = sx / (4.0 * m2);                       // adgh:8613
  if (!(lim.hi > lim.lo)) { lim.lo = 0.0; lim.hi = 0.0; }
  return lim;
}

double polrad_eta_min_ur(double x_a) {
  if (!(x_a > 0.0) || !(x_a < 1.0)) return 0.0;
  return x_a * x_a / (4.0 * (1.0 - x_a));
}

double polrad_tpeak_quadrature(
    const std::function<double(double, double, double)>& integrand,
    double x_a, double y, double s_a, EtaLimits lim, int n_eta) {
  if (!(x_a > 0.0) || !(y > 0.0) || !(y < 1.0) || !(s_a > 0.0)) return 0.0;
  if (!(lim.hi > lim.lo) || !(lim.lo > 0.0)) return 0.0;
  const double yp = (1.0 + (1.0 - y) * (1.0 - y)) / (1.0 - y);   // Y_+
  const double xa2 = x_a * x_a;
  const double acc = gl_log_quadrature(
      [&](double eta) {
        // POLRAD Eq. (39); adgh:8887 `xx1=xa**2 + 4.*xa*eta - 4.*eta`.
        const double x1 = xa2 + 4.0 * x_a * eta - 4.0 * eta;
        const double xt = x1 / (2.0 * eta * xa2);
        return integrand(eta, xt, x1);
      },
      lim.lo, lim.hi, n_eta);
  const double a3 = ALPHA_EM * ALPHA_EM * ALPHA_EM;
  // The LEADING MINUS of Eq. (18) (polrad2t.tex:580), which the
  // ultrarelativistic Eq. (38) does not print.  Without it every unpolarised
  // entry comes out NEGATIVE (Xt < 0 over essentially the whole range) and
  // `rc_tail` would ANTI-dilute A_zz.  T8(0) is the gate.
  return -(a3 / s_a) * yp * acc;
}

double polrad_sigma_el_u(const Spin1ElasticFF& ff, double x_a, double y,
                         double s_a, double m_a, int n_eta) {
  const double m2 = m_a * m_a;
  const std::function<double(double, double, double)> integ =
      [&](double eta, double xt, double) {
        const double t = 4.0 * m2 * eta;
        const double fc = ff.fc(t);
        const double fm = ff.fm(t);
        const double fq = ff.fq(t);
        // Xt Im_2|_{Q_N=0} - Im_1|_{Q_N=0}/eta   (adgh:8893-8894)
        return (fc * fc + (8.0 / 9.0) * eta * eta * fq * fq +
                (2.0 / 3.0) * eta * fm * fm) * xt
               - (2.0 / 3.0) * (1.0 + eta) * fm * fm;
      };
  return polrad_tpeak_quadrature(integ, x_a, y, s_a,
                                 polrad_eta_limits(x_a, y, s_a, m_a), n_eta);
}

double polrad_sigma_el_t(const Spin1ElasticFF& ff, double x_a, double y,
                         double s_a, double m_a, int n_eta) {
  const double m2 = m_a * m_a;
  const double xa2 = x_a * x_a;
  const std::function<double(double, double, double)> integ =
      [&](double eta, double xt, double x1) {
        const double t = 4.0 * m2 * eta;
        const double fc = ff.fc(t);
        const double fm = ff.fm(t);
        const double fq = ff.fq(t);
        // adgh:8985-8987, character for character:
        //   (1 + eta + (.75 xa^2 - eta) Xt) fm^2
        //   - Xt X1/(1+eta) fq (3 fc + 3 eta fm + eta fq)
        //   - 2 eta Xt fq (4 fc - 3 xa fm + 4/3 eta fq)
        return (1.0 + eta + (0.75 * xa2 - eta) * xt) * fm * fm
               - xt * x1 / (1.0 + eta) * fq *
                     (3.0 * fc + 3.0 * eta * fm + eta * fq)
               - 2.0 * eta * xt * fq *
                     (4.0 * fc - 3.0 * x_a * fm + (4.0 / 3.0) * eta * fq);
      };
  return polrad_tpeak_quadrature(integ, x_a, y, s_a,
                                 polrad_eta_limits(x_a, y, s_a, m_a), n_eta);
}

double pauli_suppression(double q_gev, double kf_gev) {
  // de Forest & Walecka, Adv. Phys. 15 (1966) 1, the Fermi-gas sum rule
  // POLRAD codes verbatim in `ffquas` (adgh:5605-5613):
  //     S(q) = (3/4)(q/k_F) - (1/16)(q/k_F)^3      q < 2 k_F
  //          = 1                                   q >= 2 k_F
  // (continuous at q = 2 k_F: 3/2 - 1/2 = 1).  `kf_gev <= 0` disables it,
  // which is the S_E = S_M = 1 of POLRAD Eq. (44) as v0 shipped it.
  if (!(kf_gev > 0.0) || !(q_gev > 0.0)) return 1.0;
  const double u = q_gev / kf_gev;
  if (u >= 2.0) return 1.0;
  return 0.75 * u - u * u * u / 16.0;
}

double polrad_sigma_qe_u(int z, int n, double x, double y, double s,
                         double m_n, int n_eta, double kf_gev) {
  const double m2 = m_n * m_n;
  const EtaLimits lim = polrad_eta_limits(x, y, s, m_n);
  auto one = [&](bool proton) {
    const std::function<double(double, double, double)> integ =
        [&, proton](double eta, double xt, double) {
          const double t = 4.0 * m2 * eta;
          const NucleonFF f = nucleon_ff(t);
          const double ge = proton ? f.ge_p : f.ge_n;
          const double gm = proton ? f.gm_p : f.gm_n;
          // POLRAD Eq. (44): F_e^2 -> S_E F_e^2, F_m^2 -> S_M F_m^2.  With
          // S_E = S_M = S(q) the whole integrand -- which is a combination of
          // G_E^2 and G_M^2 ALONE, no G_E G_M interference -- carries one
          // factor S.  `q` is the three-momentum transfer at the elastic
          // vertex in the target rest frame, q^2 = Q^2(1 + Q^2/4m^2) =
          // t(1 + eta), which is the argument `ffquas` suppresses on.
          const double sup =
              pauli_suppression(std::sqrt(t * (1.0 + eta)), kf_gev);
          // Eq. (40) + Eq. (A.5):  F1^2 + eta F2^2 = (G_E^2 + eta G_M^2)/(1+eta)
          //                        (F1 + F2)^2     = G_M^2
          return sup * ((ge * ge + eta * gm * gm) / (1.0 + eta) * xt -
                        gm * gm);
        };
    return polrad_tpeak_quadrature(integ, x, y, s, lim, n_eta);
  };
  return static_cast<double>(z) * one(true) +
         static_cast<double>(n) * one(false);
}

// -------------------------- POLRAD Eq. (18): the exact tau_A quadrature
//
// WHY `long double` INSIDE THIS BLOCK, and it is not a style choice.
// Eq. (B.3)'s a_ik for i = 5, 6 are {Q^2 - 3(eta q)^2, 6(eta q), -3}/M^2, i.e.
// the expansion of -M_A^2 k_n(q - k) = Q^2 - 3(eta(q-k))^2 in powers of the
// photon contraction, and in the collinear region the three terms cancel:
// MEASURED at 6Li config 1, x = 1e-4, y = 0.5, tau near tau_min,
// SUM|term| / |SUM term| = 2.9e4 inside one theta_5j and 6.6e4 at x = 0.1,
// i.e. ~4.8 decimal digits gone before the i-sum begins.  In `double` the
// TENSOR column of the tail then moves by 14 % between tanh-sinh steps
// h = 1/16 and h = 1/32 (-3.30e-4 -> -3.78e-4 at x = 1e-3, y = 0.5) -- the
// quadrature is converged and the ARITHMETIC is not.  In `long double`
// (64-bit mantissa on x86-64) the same scan is stable to 6e-5.  The
// unpolarised column never needed it (10 digits either way); it is carried
// along because it shares the theta evaluation.
//
// A PLATFORM WHERE `long double` IS `double` WOULD SILENTLY LOSE THAT, so
// `RcModel`'s constructor REFUSES `PolradFull` there rather than ship a
// tensor tail with two significant figures.  T14 and T19(e) are the gates.
namespace {
using ld = long double;

/// The eight Appendix-B F-functions at one "eta level": the bare set of
/// Eq. (B.12), or its image under one or two applications of the Eq. (B.8)
/// lift.  `FIR` is m^2 F_2+ - Q_m^2 F_d at the SAME level -- Eq. (B.8) does
/// not list F_IR^eta because it is not independent, and `adgh`'s `ffu`
/// (adgh:1125) rebuilds it from its arguments at every level:
/// `hi2`/`shi2`/`ehi2`/`ohi2` = `aml2*bis - ym*bi12` etc. at adgh:1134-1137,
/// with `ym = y + al2` = Q_m^2 (adgh:751).  (The "adgh:9134-9137" this line
/// carried until 2026-09-06 is inside `al2ll` and reads `refss=0d0`.)
struct FullFSet {
  ld F = 0, F2p = 0, F2m = 0, Fd = 0, F1p = 0, FIR = 0;
};

/// The invariants Eq. (18) is written in, all at the NUCLEAR map when the
/// caller is elastic and at the NUCLEON map when it is quasi-elastic -- the
/// formulae do not care which, only that S, X, M and the form factors agree.
struct FullKin {
  ld S, X, Sx, Sp, Q2, Qm2, M, M2, m2, lQ, sqlQ;
  ld apq, apn;    ///< (eta q) and (eta K), Eq. (B.16); 0 for a spin-1/2 target
  ld se, re0, re1;///< s_eta and r_eta = re0 + re1 tau, Eq. (B.10)
};

FullKin full_kin(ld x_a, ld y, ld s_a, ld m_a, ld m_lep) {
  FullKin k;
  k.S = s_a;
  k.X = (1.0L - y) * s_a;
  k.Sx = k.S - k.X;
  k.Sp = k.S + k.X;
  k.Q2 = x_a * y * s_a;
  k.m2 = m_lep * m_lep;
  k.Qm2 = k.Q2 + 2.0L * k.m2;
  k.M = m_a;
  k.M2 = m_a * m_a;
  k.lQ = k.Sx * k.Sx + 4.0L * k.M2 * k.Q2;
  k.sqlQ = std::sqrt(k.lQ);
  // LONGITUDINAL target polarisation, POLRAD `conkin` `+self,if=long`
  // (adgh:779-781): eta = 2(a_eta k_1 + c_eta p) with
  // a_eta = M/sqrt(lambda_s), c_eta = -S/(2 M sqrt(lambda_s)).  That is the
  // unique unit spacelike vector in the (k_1, p) plane with eta.p = 0:
  // eta.p = a_eta S + 2 c_eta M^2 = 0 and eta^2 = (4 M^2 m^2 - S^2)/lambda_s
  // = -1 exactly.  In the target rest frame it is the beam axis, which for
  // this library's collider variables is the ion's own momentum -- the same
  // axis Eq. (43) rotates off.
  const ld ls = k.S * k.S - 4.0L * k.m2 * k.M2;
  const ld sqls = std::sqrt(ls);
  const ld ae = k.M / sqls, be = 0.0L, ce = -k.S / (2.0L * k.M * sqls);
  k.apq = -k.Q2 * (ae - be) + ce * k.Sx;                 // adgh:788
  k.apn = (k.Q2 + 4.0L * k.m2) * (ae + be) + ce * k.Sp;  // adgh:789
  k.se = ae + be;                                        // s_eta, Eq. (B.10)
  k.re0 = 2.0L * ce;                        // r_eta = 2c + tau(a - b)
  k.re1 = ae - be;
  return k;
}

/// Eq. (B.12) on Eq. (B.13)'s B_{1,2}, C_{1,2}, with Eq. (B.14) for F_d.
/// adgh `tails`:1053-1065 (adgh:1066-1069 is `sps`/`spe`/`ccpe`/`ccps`, i.e.
/// s_eta and r_eta, which `full_kin` carries instead).
void full_level0(const FullKin& k, ld ta, FullFSet* out, ld* fi, ld* fii) {
  const ld b1 = -0.5L * (k.lQ * ta + k.Sp * (k.Sx * ta + 2.0L * k.Q2));
  const ld b2 = -0.5L * (k.lQ * ta - k.Sp * (k.Sx * ta + 2.0L * k.Q2));
  const ld com = 4.0L * k.m2 * (k.Q2 + ta * k.Sx - ta * ta * k.M2);
  const ld c1 = (k.S * ta + k.Q2) * (k.S * ta + k.Q2) + com;
  const ld c2 = (k.X * ta - k.Q2) * (k.X * ta - k.Q2) + com;
  const ld sc1 = std::sqrt(c1), sc2 = std::sqrt(c2);
  out->F = 1.0L / k.sqlQ;
  // Eq. (B.14), NOT the printed tau^{-1}(C_2^{-1/2} - C_1^{-1/2}): tau = 0 is
  // inside the range and the printed form is 0/0 there.
  out->Fd = k.Sp * (ta * k.Sx + 2.0L * k.Q2) / (sc1 * sc2 * (sc1 + sc2));
  out->F1p = 1.0L / sc2 + 1.0L / sc1;
  out->F2p = b2 / (sc2 * c2) - b1 / (sc1 * c1);
  out->F2m = b2 / (sc2 * c2) + b1 / (sc1 * c1);
  out->FIR = k.m2 * out->F2p - k.Qm2 * out->Fd;
  *fi = -b1 / (k.lQ * k.sqlQ);
  *fii = 0.5L * (3.0L * b1 * b1 - k.lQ * c1) / (k.lQ * k.lQ * k.sqlQ);
}

/// Eq. (B.8), FIRST lift (F -> F^eta).  adgh `tails`:1073-1075 (`eis`, `eir`,
/// `ei12` = F_2+^eta, F_2-^eta, F_d^eta), 1089 (`ei1pi2` = F_1+^eta) and 1092
/// (`eb` = F^eta); the same lift in the S direction is adgh:1070-1072.  (This
/// citation read "1075-1078, 1088, 1090" until 2026-09-06 -- 1076-1078 is
/// `ois`, the eta-s CROSS set, not this lift.)
FullFSet full_lift1(const FullKin& k, ld ta, const FullFSet& b, ld fi, ld se,
                    ld re) {
  FullFSet o;
  o.F2p = (2.0L * b.F1p * se + b.F2m * se * ta + b.F2p * re) / 2.0L;
  o.F2m = (2.0L * b.Fd * se * ta + b.F2m * re + b.F2p * se * ta) / 2.0L;
  o.Fd = (b.Fd * re + b.F1p * se) / 2.0L;
  o.F1p = (4.0L * b.F * se + b.Fd * se * ta * ta + b.F1p * re) / 2.0L;
  o.F = ((re - se * ta) * b.F + 2.0L * fi * se) / 2.0L;
  o.FIR = k.m2 * o.F2p - k.Qm2 * o.Fd;
  return o;
}

/// Eq. (B.8), SECOND lift (F -> F^{eta eta}), built from the BARE set exactly
/// as adgh `tails`:1083-1094 does (`eeis`, `eeir`, `eei12`, `eei1i2`, `eeb`).
/// TWO of these five differ from
/// polrad2t.tex's Appendix B and the FORTRAN is right in both:
///
///   * `4 F_{2-}^{eta eta}` carries `2(2F_d + F_2+) tau r_eta s_eta`; the
///     paper prints it with NO tau, which contradicts its own first-order
///     `2F_{2-}^eta = (2F_d + F_2+) tau s_eta + F_2- r_eta`.
///   * `4 F^{eta eta}` carries `4 F_i (r_eta - tau s_eta) s_eta`; the paper
///     prints `4 F_i (r_eta - tau s_eta)`, which is FIRST order in (r, s)
///     where every other term is second.
///
/// Both were caught by the x_A -> 0 gate, not by inspection: with the paper's
/// forms the tensor column does not reduce to Eq. (38)'s sigma_q^d.
FullFSet full_lift2(const FullKin& k, ld ta, const FullFSet& b, ld fi, ld fii,
                    ld se, ld re) {
  const ld rr = re * re + se * se * ta * ta;
  FullFSet o;
  o.F2p = (rr * b.F2p + 8.0L * b.F * se * se + 4.0L * b.Fd * se * se * ta * ta +
           4.0L * b.F1p * re * se + 2.0L * b.F2m * re * se * ta) / 4.0L;
  o.F2m = (rr * b.F2m + 4.0L * b.Fd * re * se * ta +
           4.0L * b.F1p * se * se * ta + 2.0L * b.F2p * re * se * ta) / 4.0L;
  o.Fd = (rr * b.Fd + 4.0L * b.F * se * se + 2.0L * b.F1p * re * se) / 4.0L;
  o.F1p = (rr * b.F1p + 4.0L * (2.0L * re - se * ta) * b.F * se +
           8.0L * fi * se * se + 2.0L * b.Fd * re * se * ta * ta) / 4.0L;
  o.F = ((re - se * ta) * (re - se * ta) * b.F +
         4.0L * (re - se * ta) * fi * se + 4.0L * fii * se * se) / 4.0L;
  o.FIR = k.m2 * o.F2p - k.Qm2 * o.Fd;
  return o;
}

/// Eq. (B.4) + Eq. (B.5) + the Eq. (B.6) lift, i.e. adgh's `ffu`, for the four
/// rows this file needs: `r1` = T_{1j n} (= T_{5j n} by Eq. (B.5)), `r2` =
/// T_{2j n} (= T_{6j n}), `r7` = T_{7j n}, `r8` = T_{8j n}.  T_3 and T_4 are
/// NOT here: every term of them carries P_L (Eq. (B.4)), and the whole A_zz
/// programme runs at an unpolarised beam.
///
/// `b` is the level-(n-1) set, `nx` the level-n one, `fdd` the level-(n+1)
/// F_d that T_{73} alone needs (adgh's `ffu` reads it from /bseo/, which is
/// correct because row 7 is only ever used at n = 1).
struct FullTm {
  ld r1[3], r2[3], r7[3], r8[3];
};
FullTm full_ffu(const FullKin& k, ld ta, const FullFSet& b, const FullFSet& nx,
                ld fdd) {
  FullTm o;
  o.r1[0] = -4.0L * (2.0L * k.m2 - k.Q2) * b.FIR;
  o.r1[1] = 4.0L * b.FIR * ta;
  o.r1[2] = -2.0L * (2.0L * b.F + b.Fd * ta * ta);
  o.r2[0] = ((k.Sp * k.Sp - k.Sx * k.Sx) - 4.0L * k.M2 * k.Q2) * b.FIR /
            (2.0L * k.M2);
  o.r2[1] = (2.0L * k.m2 * b.F2m * k.Sp - 4.0L * k.M2 * b.FIR * ta -
             b.Fd * k.Sp * k.Sp * ta + b.F1p * k.Sp * k.Sx +
             2.0L * b.FIR * k.Sx) / (2.0L * k.M2);
  o.r2[2] = (2.0L * (2.0L * b.F + b.Fd * ta * ta) * k.M2 + 4.0L * k.m2 * b.Fd -
             b.Fd * k.Sx * ta - b.F1p * k.Sp) / (2.0L * k.M2);
  o.r7[0] = -2.0L * (4.0L * k.m2 + 3.0L * k.apn * k.apn -
                     3.0L * k.apq * k.apq + k.Q2) * b.FIR;
  o.r7[1] = -2.0L * (6.0L * k.m2 * k.apn * nx.F2m -
                     3.0L * k.apn * k.apn * b.Fd * ta +
                     3.0L * k.apn * k.apq * b.F1p + 6.0L * k.apq * nx.FIR +
                     b.FIR * ta);
  o.r7[2] = -(24.0L * k.m2 * fdd - 6.0L * k.apn * nx.F1p -
              6.0L * k.apq * nx.Fd * ta - 2.0L * b.F - b.Fd * ta * ta);
  o.r8[0] = -3.0L * (k.apn * k.Sp - k.apq * k.Sx) * b.FIR / k.M;
  // THE PAIRING HERE IS `adgh`'s, NOT the paper's.  polrad2t.tex's T_{821}
  // prints `(S eta k_1 + X eta k_2) F_{1+}` = (apn S_p + apq S_x)F_1+/2;
  // adgh's `tm3(6,2,n)` (adgh:1161-1163, in `ffu`; the "adgh:9155-9157" this
  // comment carried until 2026-09-06 is `al2ll`'s own `write`/`end`) has
  // `apn F_1+ S_x + apq F_1+ S_p`, i.e.
  // (S eta k_1 - X eta k_2).  With the paper's the tensor column misses
  // Eq. (38)'s sigma_q^d by x14 in the x_A -> 0 limit; with adgh's it lands
  // on 1.003 at x_A = 0.003.  MEASURED, not argued.
  o.r8[1] = -3.0L * (2.0L * (k.apn * b.F2m + nx.F2m * k.Sp) * k.m2 -
                     2.0L * b.Fd * k.Sp * ta * k.apn +
                     k.apn * b.F1p * k.Sx + k.apq * b.F1p * k.Sp +
                     2.0L * b.FIR * k.apq + 2.0L * nx.FIR * k.Sx) /
            (2.0L * k.M);
  o.r8[2] = -3.0L * (8.0L * k.m2 * nx.Fd - k.apn * b.F1p -
                     k.apq * b.Fd * ta - nx.Fd * k.Sx * ta -
                     nx.F1p * k.Sp) / (2.0L * k.M);
  return o;
}

/// theta_ij(tau) of Eq. (B.1) for i in {1, 2, 5, 6, 7, 8}, laid out as
/// `th[i][j-1]` with j running to k_i = (3, 3, -, -, 5, 5, 3, 4).
struct FullTheta { ld th[9][5]; };

FullTheta full_thetas(const FullKin& k, ld ta, bool tensor) {
  FullFSet l0; ld fi = 0, fii = 0;
  full_level0(k, ta, &l0, &fi, &fii);
  FullTheta t;
  for (int i = 0; i < 9; ++i)
    for (int j = 0; j < 5; ++j) t.th[i][j] = 0.0L;
  if (!tensor) {
    // Spin 1/2, or the Q_N = 0 half of spin 1: only i = 1, 2 survive and the
    // eta four-vector is never needed.  This is the whole quasi-elastic path.
    const FullFSet z;
    const FullTm n1 = full_ffu(k, ta, l0, z, 0.0L);
    for (int j = 0; j < 3; ++j) { t.th[1][j] = n1.r1[j]; t.th[2][j] = n1.r2[j]; }
    return t;
  }
  const ld se = k.se, re = k.re0 + k.re1 * ta;
  const FullFSet l1 = full_lift1(k, ta, l0, fi, se, re);
  const FullFSet l2 = full_lift2(k, ta, l0, fi, fii, se, re);
  const FullFSet z;
  const FullTm n1 = full_ffu(k, ta, l0, l1, l2.Fd);
  const FullTm n2 = full_ffu(k, ta, l1, l2, 0.0L);
  const FullTm n3 = full_ffu(k, ta, l2, z, 0.0L);
  const FullTm* nn[3] = {&n1, &n2, &n3};
  for (int j = 0; j < 3; ++j) {
    t.th[1][j] = n1.r1[j];
    t.th[2][j] = n1.r2[j];
    t.th[7][j] = n1.r7[j];
  }
  // Eq. (B.3)'s a_ik, IN adgh's NORMALISATION -- with the 1/M^{l_i - 1} the
  // paper leaves implicit in Eq. (A.1)'s k_n = (3(q eta)^2 - Q^2)/M^2 and
  // Omega~/M^2.  `bornin` (adgh:817-825) settles it independently: the BORN
  // kernels obey tm(5) = -ek tm(1) with ek = (3 apq^2 - Q^2)/M^2 exactly, and
  // tm(8) = (apq/M) x (the i = 8 kernel).
  const ld a8[2] = {k.apq / k.M, -1.0L / k.M};
  const ld a5[3] = {(k.Q2 - 3.0L * k.apq * k.apq) / k.M2,
                    6.0L * k.apq / k.M2, -3.0L / k.M2};
  // Eq. (B.7)'s q_ik, IN adgh's FORM (adgh:1115-1116): the extra term is
  // tau/M^2 times T_{i,j-1,1} and is NOT multiplied by a_i2.  The paper's
  // "q_ik = tau/M inside T_{ij2}" would put a_i2 = 6(eta q)/M^2 in front of
  // it; MEASURED, that misses Eq. (38)'s sigma_q^d by five orders of
  // magnitude in the x_A -> 0 limit, while adgh's lands on 1.003.
  const ld qk = ta / k.M2;
  for (int kk = 1; kk <= 2; ++kk)
    for (int j = kk; j <= kk + 2; ++j)
      t.th[8][j - 1] += nn[kk - 1]->r8[j - kk] * a8[kk - 1];
  for (int kk = 1; kk <= 3; ++kk)
    for (int j = kk; j <= kk + 2; ++j) {
      t.th[5][j - 1] += nn[kk - 1]->r1[j - kk] * a5[kk - 1];
      t.th[6][j - 1] += nn[kk - 1]->r2[j - kk] * a5[kk - 1];
      if (kk == 2) {
        t.th[5][j - 1] += n1.r1[j - 2] * qk;
        t.th[6][j - 1] += n1.r2[j - 2] * qk;
      }
    }
  return t;
}

/// tanh-sinh (double-exponential) quadrature of `f` on [a, b] with `n` nodes.
///
/// WHY NOT THE Gauss-Legendre PANELS `gl_log_quadrature` USES.  The s- and
/// p-peaks sit AT tau_s and tau_p with a width set by 4 m_e^2 -- a relative
/// scale of 1e-11 on a range of order 1.  Splitting the range at the two peak
/// positions and putting a rule that clusters its nodes DOUBLE-exponentially
/// at every panel edge resolves both sides of both peaks with no adaptivity
/// and no peak-width estimate.  The abscissa is written as a + (b-a)*s with
/// s the logistic of pi/2 sinh(u), so the offset from the edge never cancels.
template <class F>
ld tanh_sinh_quadrature(const F& f, ld a, ld b, int n) {
  if (!(b > a)) return 0.0L;
  const ld umax = 3.5L;
  const ld h = 2.0L * umax / static_cast<ld>(std::max(2, n) - 1);
  const ld len = b - a;
  const int half = static_cast<int>(umax / h);
  ld acc = 0.0L;
  for (int i = -half; i <= half; ++i) {
    const ld u = static_cast<ld>(i) * h;
    const ld v = 0.5L * static_cast<ld>(kPi) * std::sinh(u);
    const ld s = 1.0L / (1.0L + std::exp(-2.0L * v));
    const ld sm = 1.0L / (1.0L + std::exp(2.0L * v));
    const ld w = len * static_cast<ld>(kPi) * s * sm * std::cosh(u);
    if (!(w > 0.0L)) continue;
    const ld fv = f(a + len * s);
    if (!std::isfinite(static_cast<double>(fv))) continue;
    acc += w * fv;
  }
  return acc * h;
}

/// The shared driver: `im` fills Im^el_i(t, eta) at one tau node, and the
/// caller says which of the eight it fills.
template <class IM>
ld polrad_full_quadrature(const FullKin& k, const IM& im, bool tensor,
                          int n_tau) {
  if (!(k.Q2 > 0.0L) || !(k.Sx > k.Q2)) return 0.0L;
  const ld tamax = (k.Sx + k.sqlQ) / (2.0L * k.M2);
  const ld tamin = -k.Q2 / (k.M2 * tamax);        // = (S_x - sqrt(lQ))/(2M^2)
  if (!(tamax > tamin)) return 0.0L;
  static const int kmax[9] = {0, 3, 3, 0, 0, 5, 5, 3, 4};
  auto integ = [&](ld ta) -> ld {
    const ld rel = (k.Sx - k.Q2) / (1.0L + ta);   // R_el, Eq. (17)
    const ld t = k.Q2 + rel * ta;
    if (!(t > 0.0L)) return 0.0L;
    const ld eta = t / (4.0L * k.M2);
    ld imv[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    im(static_cast<double>(t), static_cast<double>(eta), imv);
    const FullTheta th = full_thetas(k, ta, tensor);
    const ld pref = 2.0L * k.M2 / ((1.0L + ta) * t * t);
    ld acc = 0.0L;
    for (int i = 1; i <= 8; ++i) {
      if (imv[i] == 0.0L) continue;
      ld rp = 1.0L / rel;                          // R^{j-2} at j = 1
      ld sub = 0.0L;
      for (int j = 1; j <= kmax[i]; ++j) { sub += th.th[i][j - 1] * rp; rp *= rel; }
      acc += sub * imv[i];
    }
    return acc * pref;
  };
  // The peaks go on panel EDGES: tau_s = -Q^2/S (always inside, because
  // tau_min = -2Q^2/(S_x + sqrt(lQ)) < -Q^2/S_x <= -Q^2/S) and tau_p = Q^2/X.
  ld bp[5] = {tamin, -k.Q2 / k.S, 0.0L, k.Q2 / k.X, tamax};
  ld tot = 0.0L;
  ld lo = tamin;
  for (int i = 1; i < 5; ++i) {
    const ld hi = (i == 4) ? tamax : std::min(bp[i], tamax);
    if (!(hi > lo)) continue;
    tot += tanh_sinh_quadrature(integ, lo, hi, n_tau);
    lo = hi;
  }
  // Eq. (18)'s LEADING MINUS (polrad2t.tex:580) -- the same one
  // `polrad_tpeak_quadrature` applies, and for the same reason (T8(0)).
  const ld a3 = static_cast<ld>(ALPHA_EM) * static_cast<ld>(ALPHA_EM) *
                static_cast<ld>(ALPHA_EM);
  // y_Bj = S_x/S_A exactly (X_A = (1-y)S_A).
  return -a3 * (k.Sx / k.S) * tot;
}

}  // namespace

TauLimits polrad_tau_limits(double x_a, double y, double s_a, double m_a) {
  TauLimits lim;
  if (!(x_a > 0.0) || !(y > 0.0) || !(y < 1.0) || !(s_a > 0.0) ||
      !(m_a > 0.0)) {
    return lim;
  }
  const double q2 = x_a * y * s_a;
  const double sx = y * s_a;
  const double m2 = m_a * m_a;
  const double sqly = std::sqrt(sx * sx + 4.0 * m2 * q2);
  lim.hi = (sx + sqly) / (2.0 * m2);
  // tau_min tau_max = (S_x^2 - lambda_Q)/(4 M^4) = -Q^2/M^2 exactly, which is
  // how the printed (S_x - sqrt(lambda_Q))/(2M^2) is evaluated without
  // cancelling two nearly equal large numbers.
  lim.lo = (lim.hi > 0.0) ? -q2 / (m2 * lim.hi) : 0.0;
  return lim;
}

PolradIm polrad_im_el_spin1(const Spin1ElasticFF& ff, double t_gev2,
                            double m_a) {
  PolradIm o;
  if (!(m_a > 0.0)) return o;
  const double eta = t_gev2 / (4.0 * m_a * m_a);
  const double fc = ff.fc(t_gev2), fm = ff.fm(t_gev2), fq = ff.fq(t_gev2);
  // polrad2t.tex:2698-2709, Eq. (A.4), verbatim and split in Q_N.
  o.u[1] = (2.0 / 3.0) * eta * (1.0 + eta) * fm * fm;
  o.u[2] = fc * fc + (2.0 / 3.0) * eta * fm * fm +
           (8.0 / 9.0) * eta * eta * fq * fq;
  o.t[1] = eta * eta * fm * fm;
  o.t[2] = eta * fm * fm + (4.0 * eta * eta / (1.0 + eta)) *
                               ((eta / 3.0) * fq + fc - fm) * fq;
  o.t[5] = 0.25 * fm * fm;
  // 4/(1 + eta_A), NOT 4 eta_A/(1 + eta_A): design_C_tensor_rc.md sec. 1.4.6
  // carried the spurious eta_A until 2026-09-04 and
  // polrad_transcription_check.md sec. 7.1 is the correction.
  o.t[6] = 0.25 * (fm * fm + (4.0 / (1.0 + eta)) *
                                 ((eta / 3.0) * fq + fc + eta * fm) * fq);
  o.t[7] = eta * (1.0 + eta) * fm * fm;
  o.t[8] = -eta * fm * (fm + 2.0 * fq);
  return o;
}

PolradFullPair polrad_full_sigma_el(const Spin1ElasticFF& ff, double x_a,
                                    double y, double s_a, double m_a,
                                    int n_tau, double m_lepton) {
  PolradFullPair out;
  if (!(x_a > 0.0) || !(y > 0.0) || !(y < 1.0) || !(s_a > 0.0) ||
      !(m_a > 0.0) || !(m_lepton > 0.0)) {
    return out;
  }
  const FullKin k = full_kin(x_a, y, s_a, m_a, m_lepton);
  // Eq. (A.4) at Q_N = 0 -- Im_1|_0 = B/2, Im_2|_0 = A of the spin-1
  // Rosenbluth pair, the SAME algebra `rosenbluth_spin1` writes and not a
  // second definition of it (T9 pins the two against each other).
  auto im_u = [&](double t, double, ld* imv) {
    const PolradIm im = polrad_im_el_spin1(ff, t, m_a);
    for (int i = 1; i <= 8; ++i) imv[i] = static_cast<ld>(im.u[i]);
  };
  // ... and the Q_N part, Im_i = Im_i|_0 + (Q_N/6) Im_i^T, so what is
  // returned is Eq. (37)'s sigma_q^A and NOT (Q_N/6) sigma_q^A.
  auto im_t = [&](double t, double, ld* imv) {
    const PolradIm im = polrad_im_el_spin1(ff, t, m_a);
    for (int i = 1; i <= 8; ++i) imv[i] = static_cast<ld>(im.t[i]);
  };
  out.u = static_cast<double>(
      polrad_full_quadrature(k, im_u, /*tensor=*/false, n_tau));
  out.t = static_cast<double>(
      polrad_full_quadrature(k, im_t, /*tensor=*/true, n_tau));
  return out;
}

double polrad_full_sigma_qe_u(int z, int n, double x, double y, double s,
                              double m_n, int n_tau, double kf_gev,
                              double m_lepton) {
  if (!(x > 0.0) || !(y > 0.0) || !(y < 1.0) || !(s > 0.0) || !(m_n > 0.0) ||
      !(m_lepton > 0.0)) {
    return 0.0;
  }
  const FullKin k = full_kin(x, y, s, m_n, m_lepton);
  auto one = [&](bool proton) {
    auto im = [&](double t, double eta, ld* imv) {
      const NucleonFF f = nucleon_ff(t);
      const double ge = proton ? f.ge_p : f.ge_n;
      const double gm = proton ? f.gm_p : f.gm_n;
      // POLRAD Eq. (A.5) at Z = 1: Im_1 = eta G_M^2, Im_2 = (G_E^2 + eta
      // G_M^2)/(1 + eta).  Eq. (44)'s S_E = S_M = S(q) multiplies both, at
      // the SAME q^2 = t(1 + eta) the t-peak quasi-elastic tail uses.
      const double sup =
          pauli_suppression(std::sqrt(t * (1.0 + eta)), kf_gev);
      imv[1] = static_cast<ld>(sup * eta * gm * gm);
      imv[2] = static_cast<ld>(sup * (ge * ge + eta * gm * gm) / (1.0 + eta));
    };
    return static_cast<double>(
        polrad_full_quadrature(k, im, /*tensor=*/false, n_tau));
  };
  return static_cast<double>(z) * one(true) +
         static_cast<double>(n) * one(false);
}

// ----------------------------------------- the LEADING-LOG s- and p-peaks
//
// These six functions were written for T8(c) and lived in tests/test_rc.cpp
// until this run.  They are here now because the tail MODEL uses them
// (`RcTailModel::TPeakPlusLL`), and because a gate that measures a
// test-local construction certifies the construction and not the shipped
// code.  Nothing in them is new physics: same `nucleon_ff`, same
// `Spin1ElasticFF`, same ALPHA_EM / M_ELECTRON / PROTON_MASS / kPi.
//
// WHAT THEY ARE, stated once so no caller has to guess.  POLRAD sec. 2.1.3 B
// asserts that for a TAIL the s- and p-peaks are suppressed, and the whole
// `TPeak` tail rests on that one sentence.  This is the standard collinear
// (equivalent-radiator) estimate of the SAME observable, sharing no line with
// the eta_A quadrature above, so that the assertion can be MEASURED instead of
// believed.  It holds where this generator runs and fails elsewhere; T8(c)
// asserts both halves.

double ll_radiator(double z, double q2) {
  if (!(z > 0.0) || !(z < 1.0)) return 0.0;
  // ln(Q^2/m_e^2) < 0 below the electron mass: not a small radiator but a
  // NEGATIVE one, and a negative cross section is worse than a missing term.
  // DEFENSIVE, and measured to be so: the lowest tail-table node of the
  // shipped 6Li scenarios sits at Q^2 = 1.59e-3 GeV^2 (config 1) and
  // 3.96e-3 (config 2), i.e. 6097x and 1.5e4x above m_e^2.  What could reach
  // it is a USER scenario: `build_tail_tables` floors y at max(1e-6, y_min),
  // so at the sampler's x_min = 1e-4 and config 1's s any y_min below
  // 6.6e-4 crosses over -- the generator scenario's 0.004 clears it by 6.1.
  if (!(q2 > M_ELECTRON * M_ELECTRON)) return 0.0;
  return ALPHA_EM / kPi * std::log(q2 / (M_ELECTRON * M_ELECTRON)) *
         (1.0 + z * z) / (1.0 - z);
}

RosenbluthAB rosenbluth_spin1(const Spin1ElasticFF& ff, double q2, double m) {
  const double eta = q2 / (4.0 * m * m);
  const double fc = ff.fc(q2), fm = ff.fm(q2), fq = ff.fq(q2);
  RosenbluthAB ab;
  ab.a = fc * fc + (8.0 / 9.0) * eta * eta * fq * fq +
         (2.0 / 3.0) * eta * fm * fm;
  ab.b = (4.0 / 3.0) * eta * (1.0 + eta) * fm * fm;
  return ab;
}

RosenbluthAB rosenbluth_nucleon(bool proton, double q2) {
  const NucleonFF f = nucleon_ff(q2);
  const double ge = proton ? f.ge_p : f.ge_n;
  const double gm = proton ? f.gm_p : f.gm_n;
  const double tau = q2 / (4.0 * PROTON_MASS * PROTON_MASS);
  RosenbluthAB ab;
  ab.a = (ge * ge + tau * gm * gm) / (1.0 + tau);
  ab.b = 2.0 * tau * gm * gm;
  return ab;
}

double dsigma_el_dq2(const RosenbluthAB& ab, double q2p, double sp, double m) {
  if (!(q2p > 0.0) || !(sp > 0.0)) return 0.0;
  const double yp = q2p / sp;
  const double etap = q2p / (4.0 * m * m);
  const double kin = 1.0 - yp - m * m * yp * yp / q2p;
  if (!(kin > 0.0)) return 0.0;
  return 4.0 * kPi * ALPHA_EM * ALPHA_EM / (q2p * q2p) *
         (0.5 * yp * yp * ab.b / (2.0 * etap) + kin * ab.a);
}

LlPeaks ll_peaks_spin1(const Spin1ElasticFF& ff, double x_a, double y,
                       double s_a, double m) {
  const double q2 = x_a * y * s_a;
  const double zs = (1.0 - y) / (1.0 - x_a * y);
  const double zp = 1.0 - y * (1.0 - x_a);
  LlPeaks out;
  out.s = y * s_a * ll_radiator(zs, q2) *
          dsigma_el_dq2(rosenbluth_spin1(ff, zs * q2, m), zs * q2, zs * s_a,
                        m) *
          zs / (1.0 - x_a * y);
  out.p = y * s_a * ll_radiator(zp, q2) *
          dsigma_el_dq2(rosenbluth_spin1(ff, q2 / zp, m), q2 / zp, s_a, m) /
          zp;
  return out;
}

LlPeaks ll_peaks_qe(int z, int n, double x, double y, double s,
                    double kf_gev) {
  const double q2 = x * y * s;
  const double zs = (1.0 - y) / (1.0 - x * y);
  const double zp = 1.0 - y * (1.0 - x);
  // The SAME Fermi-gas S(q) `polrad_sigma_qe_u` applies inside its eta
  // integrand, here at each peak's OWN elastic-vertex q -- q'^2 = t'(1 + eta')
  // with t' the shifted Q'^2.  The Rosenbluth (A, B) is a combination of
  // G_E^2 and G_M^2 alone (no interference), so S multiplies it as a whole,
  // exactly as in Eq. (44).  Default kf_gev = 0 => S = 1, which is what T8(c)
  // compares against and what keeps this function's own gates unchanged.
  auto sup = [&](double t) {
    return pauli_suppression(
        std::sqrt(t * (1.0 + t / (4.0 * PROTON_MASS * PROTON_MASS))), kf_gev);
  };
  const double sup_s = sup(zs * q2), sup_p = sup(q2 / zp);
  LlPeaks out;
  for (int k = 0; k < 2; ++k) {
    const bool proton = (k == 0);
    const double mult = proton ? static_cast<double>(z)
                               : static_cast<double>(n);
    out.s += sup_s * mult * y * s * ll_radiator(zs, q2) *
             dsigma_el_dq2(rosenbluth_nucleon(proton, zs * q2), zs * q2,
                           zs * s, PROTON_MASS) *
             zs / (1.0 - x * y);
    out.p += sup_p * mult * y * s * ll_radiator(zp, q2) *
             dsigma_el_dq2(rosenbluth_nucleon(proton, q2 / zp), q2 / zp, s,
                           PROTON_MASS) /
             zp;
  }
  return out;
}

// ------------------------------------------- the 6Li spin-1 form factors

namespace {

/// The HO point-nucleon monopole, design sec. 2.1:
///   F_point(q) = [1 - (alpha/(2+3 alpha)) q^2 a^2/2] exp(-q^2 a^2/4)
double ho_point(double q2_fm2, double a_fm, double alpha) {
  const double qa2 = q2_fm2 * a_fm * a_fm;
  return (1.0 - alpha / (2.0 + 3.0 * alpha) * (0.5 * qa2)) *
         std::exp(-0.25 * qa2);
}

/// F_m's OWN shape: (1 - q^2/q_z^2) exp(-q^2 b^2/4).  It is NOT the monopole
/// above, because WS98's F_T has its first peak at q = 0.5 fm^-1, a zero, and
/// a second peak at q = 2 fm^-1 -- structure exactly where the tail lives
/// (q = 2 fm^-1 <=> t ~ 0.155 GeV^2) and a monopole has none there.
double ho_mag(double q2_fm2, double qz_fm, double b_fm) {
  return (1.0 - q2_fm2 / (qz_fm * qz_fm)) *
         std::exp(-0.25 * q2_fm2 * b_fm * b_fm);
}

double q2_fm2_of(double t_gev2) {
  const double t = std::max(t_gev2, 0.0);
  return t / (HBARC_GEV_FM * HBARC_GEV_FM);
}

/// j_0(x) = sin x / x, with the series where the quotient loses digits.  T11
/// differentiates F_point at q^2 = 1e-6 fm^-2, i.e. x <~ 0.02, and there
/// sin(x)/x cancels its two leading digits.
double j0_stable(double x) {
  const double ax = std::fabs(x);
  if (ax < 1e-2) {
    const double x2 = x * x;
    return 1.0 - x2 / 6.0 * (1.0 - x2 / 20.0 * (1.0 - x2 / 42.0));
  }
  return std::sin(x) / x;
}

/// `C0Shape::VmcFt`: the j_0 transform of `data/vmc/density/li6.density`,
/// r-rescaled to the MEASURED <r^2>_point and continued above
/// `LI6_VMC_C0_Q_CUT_FM`.  Built ONCE per process -- rc.hpp's block carries
/// the model argument, the provenance and both policy overrides.
class Li6VmcC0 {
 public:
  /// Function-local static: thread-safe under C++11, and a throw during
  /// construction (missing or malformed table) propagates to the caller and
  /// is retried on the next call rather than caching a broken object.
  static const Li6VmcC0& get() {
    static const Li6VmcC0 kOne;
    return kOne;
  }

  /// F_point(q) at the PHYSICAL q [fm^-1].  F(0) = 1 exactly.
  double f(double q_fm) const {
    const double q = std::fabs(q_fm);
    if (q <= LI6_VMC_C0_Q_CUT_FM) return sum(q);
    return f_cut_ * std::exp(-lambda_ * (q - LI6_VMC_C0_Q_CUT_FM));
  }

  double r_scale() const { return scale_; }
  double lambda_fm() const { return lambda_; }
  std::size_t n_rows() const { return r_fm_.size(); }

 private:
  Li6VmcC0();

  double sum(double q_fm) const {
    double acc = 0.0;
    const double qs = q_fm * scale_;
    for (std::size_t i = 0; i < r_fm_.size(); ++i) {
      acc += w_[i] * j0_stable(qs * r_fm_[i]);
    }
    return acc;
  }

  std::string path_;
  std::vector<double> r_fm_;   ///< the rows with rho > 0 (the rest add zero)
  std::vector<double> w_;      ///< rho r^2 / SUM rho r^2, so SUM w_ = 1
  double r2_table_ = 0.0;
  double scale_ = 1.0;
  double lambda_ = 0.0;
  double f_cut_ = 0.0;
};

Li6VmcC0::Li6VmcC0() {
  // The same `data_dir()` resolution every VMC table goes through
  // (docs/CONVENTIONS.md: "Nothing resolves a data path against the working
  // directory"), and the same `****`-ruled reader `cluster_config.cpp` uses
  // on the alpha core's he4.density -- NOT a second parser for the same
  // file format.
  path_ = data_path(VMC_LI6_DENSITY);
  const AnlTable t = read_anl_plain(path_, 2);
  if (t.x.size() < 2) {
    throw std::runtime_error("Li6VmcC0: " + path_ + " carries fewer than two "
                             "(R, RHORP, DRHORP) rows");
  }
  // The plain sum below is the MIDPOINT rule and its bin width cancels in the
  // ratio ONLY because the grid is uniform (0.05, 0.15, ... step 0.1).  Check
  // it rather than assume it: a re-fetched file on another grid would
  // silently mis-weight <r^2> and the transform together.
  const double dr = t.x[1] - t.x[0];
  if (!(dr > 0.0)) {
    throw std::runtime_error("Li6VmcC0: " + path_ + " is not increasing in R");
  }
  double s2 = 0.0, s4 = 0.0;
  for (std::size_t i = 0; i < t.x.size(); ++i) {
    if (i > 0 && std::fabs((t.x[i] - t.x[i - 1]) - dr) > 1e-6 * dr) {
      throw std::runtime_error(
          "Li6VmcC0: " + path_ + " is not on a UNIFORM R grid at row " +
          std::to_string(i) + ".  The transform is a midpoint sum whose bin "
          "width cancels in the ratio; a non-uniform grid needs per-row "
          "weights, which this shape deliberately does not carry");
    }
    const double r = t.x[i];
    const double rho = t.col[0][i];
    if (!(r > 0.0) || rho < 0.0) {
      throw std::runtime_error(
          "Li6VmcC0: " + path_ + " row " + std::to_string(i) +
          " has R <= 0 or a NEGATIVE density");
    }
    s2 += rho * r * r;
    s4 += rho * r * r * r * r;
  }
  if (!(s2 > 0.0)) {
    throw std::runtime_error("Li6VmcC0: " + path_ + " integrates to zero");
  }
  r2_table_ = s4 / s2;
  // <r^2> is put on the MEASURED point-nucleon moment, so BOTH edges of the
  // band carry the same second moment and T11 needs no edge-specific
  // tolerance.  Computed from THIS quadrature's own <r^2>, not from the
  // file's printed rms: the two differ by 1.4e-4 relative (2.44363 against
  // 2.4433) and T11 gates <r^2> at 1e-4.
  scale_ = std::sqrt(LI6_R2_POINT_FM2 / r2_table_);
  for (std::size_t i = 0; i < t.x.size(); ++i) {
    const double rho = t.col[0][i];
    if (rho == 0.0) continue;          // the r >~ 10 fm tail: exactly zero
    r_fm_.push_back(t.x[i]);
    w_.push_back(rho * t.x[i] * t.x[i] / s2);
  }
  // The continuation.  rc.hpp's block argues the choice; this only refuses a
  // table that cannot support it, rather than emitting a form factor that
  // GROWS with q.
  f_cut_ = sum(LI6_VMC_C0_Q_CUT_FM);
  const double f_match = sum(LI6_VMC_C0_Q_MATCH_FM);
  if (!(f_cut_ > 0.0) || !(f_match > f_cut_)) {
    throw std::runtime_error(
        "Li6VmcC0: the transform of " + path_ + " is not positive and falling "
        "over [" + std::to_string(LI6_VMC_C0_Q_MATCH_FM) + ", " +
        std::to_string(LI6_VMC_C0_Q_CUT_FM) + "] fm^-1 (F = " +
        std::to_string(f_match) + " -> " + std::to_string(f_cut_) +
        "), so the exponential continuation above the cut has no slope to "
        "take");
  }
  lambda_ = std::log(f_match / f_cut_) /
            (LI6_VMC_C0_Q_CUT_FM - LI6_VMC_C0_Q_MATCH_FM);
}

/// THE C0 SHAPE, shared by `F_c` and `F_q`.  The one place `C0Shape` is read.
double c0_point(const HoSpin1FFOptions& o, double q2_fm2) {
  if (o.c0_shape == C0Shape::VmcFt) {
    return Li6VmcC0::get().f(std::sqrt(std::max(q2_fm2, 0.0)));
  }
  return ho_point(q2_fm2, o.a_fm, o.alpha);
}

/// The ISOSCALAR nucleon folding G_E^p + G_E^n (N = Z = 3 for 6Li), which is
/// what turns a POINT-nucleon density into a charge form factor.  NOT G_E^p
/// alone: <r^2>_point was derived by subtracting BOTH <r^2>_p and
/// (N/Z)<r^2>_n, so folding with the proton one would be inconsistent with the
/// very number the shape is fitted to.  G_E^n(0) = 0 keeps F_c(0) = Z exact.
double nucleon_fold(const HoSpin1FF& ff, double t_gev2) {
  if (!ff.options().fold_nucleon) return 1.0;
  const NucleonFF n = nucleon_ff(t_gev2);
  return n.ge_p + n.ge_n;
}

}  // namespace

std::shared_ptr<HoSpin1FF> HoSpin1FF::for_ion(const Ion& ion,
                                              HoSpin1FFOptions opt) {
  // CONVENTIONS.md: no physics number in two places.  M_A and Z come from the
  // Ion (the beams.cpp AME mass table and LI6().Z), NEVER from a literal
  // here -- and an ion-specific DEFAULT would silently hand a 7Li or deuteron
  // run 6Li form factors, which is exactly what this factory exists to stop.
  if (opt.z == 0.0) opt.z = static_cast<double>(ion.Z);
  if (opt.m_a_gev == 0.0) opt.m_a_gev = ion.mass();
  if (ion.name == "6Li") {
    if (opt.a_fm == 0.0)     opt.a_fm     = LI6_FF_HO_A_FM;
    if (opt.alpha == 0.0)    opt.alpha    = LI6_FF_HO_ALPHA;
    if (opt.mu_n == 0.0)     opt.mu_n     = LI6_MU_N;
    if (opt.q_fm2 == 0.0)    opt.q_fm2    = LI6_QUADRUPOLE_FM2;
    if (opt.fm_qz_fm == 0.0) opt.fm_qz_fm = LI6_FF_FM_QZ_FM;
    if (opt.fm_b_fm == 0.0)  opt.fm_b_fm  = LI6_FF_FM_B_FM;
  }
  // Anything still at its 0 sentinel is a field this ion has no measured
  // block for.  THROW rather than default it: a silent 0 quadrupole moment
  // would zero the tensor tail and no unpolarised test would see it.
  std::string missing;
  auto need = [&missing](double v, const char* name) {
    if (v == 0.0) { if (!missing.empty()) missing += ", "; missing += name; }
  };
  need(opt.a_fm, "a_fm");
  need(opt.alpha, "alpha");
  need(opt.mu_n, "mu_n");
  need(opt.q_fm2, "q_fm2");
  need(opt.fm_qz_fm, "fm_qz_fm");
  need(opt.fm_b_fm, "fm_b_fm");
  if (!missing.empty()) {
    throw std::runtime_error(
        "HoSpin1FF::for_ion: ion '" + ion.name + "' has no measured-moment "
        "block in rc.hpp (today only 6Li does), and " + missing +
        " were left at the 0 sentinel.  Supply them on HoSpin1FFOptions "
        "explicitly -- design_C_tensor_rc.md sec. 2.1 -- or use an ion that "
        "has one.  (A field genuinely equal to zero cannot be expressed by "
        "the sentinel; use fq_scale = 0 / tail_tensor_scale = 0 for that.)");
  }
  if (std::fabs(ion.spin - 1.0) > 1e-9) {
    throw std::runtime_error(
        "HoSpin1FF::for_ion: '" + ion.name + "' has spin " +
        std::to_string(ion.spin) +
        "; POLRAD Eq. (A.4)'s (F_c, F_m, F_q) is the SPIN-1 line");
  }
  // `C0Shape::VmcFt` reads ONE file, 6Li's.  Every other field on this struct
  // is either ion-agnostic or filled from the ion's own block; this one is
  // not, so it gets the same refusal the measured-moment block gets rather
  // than silently handing a deuteron 6Li's charge shape.  (T10 builds a
  // deuteron HoSpin1FF with supplied moments and would be exactly that
  // caller.)
  if (opt.c0_shape == C0Shape::VmcFt && ion.name != "6Li") {
    throw std::runtime_error(
        "HoSpin1FF::for_ion: C0Shape::VmcFt is the j0 transform of "
        "6Li's OWN point-proton density (" + std::string(VMC_LI6_DENSITY) +
        "), and ion '" + ion.name + "' is not 6Li.  There is no equivalent "
        "table for it in data/vmc/density/ that this model knows how to "
        "normalise, so the shape is REFUSED rather than defaulted -- use "
        "C0Shape::Ho with that ion's own (a, alpha).");
  }
  return std::shared_ptr<HoSpin1FF>(new HoSpin1FF(opt));
}

double HoSpin1FF::fc(double t_gev2) const {
  const double q2 = q2_fm2_of(t_gev2);
  // F_c(0) = Z exactly, because G_E^n(0) = 0 and F_point(0) = 1.
  return opt_.z * c0_point(opt_, q2) * nucleon_fold(*this, t_gev2);
}

double HoSpin1FF::fm(double t_gev2) const {
  const double q2 = q2_fm2_of(t_gev2);
  // F_m(0) = (M_A/m_p) mu_A/mu_N -- the Rosenbluth G_M of POLRAD Eq. (A.4)'s
  // Q_N = 0 limit, validated on the deuteron at 1.71396 (T10).
  return (opt_.m_a_gev / PROTON_MASS) * opt_.mu_n *
         ho_mag(q2, opt_.fm_qz_fm, opt_.fm_b_fm) * opt_.tail_tensor_scale;
}

double HoSpin1FF::fq(double t_gev2) const {
  const double q2 = q2_fm2_of(t_gev2);
  const double m_fm = opt_.m_a_gev / HBARC_GEV_FM;   // M_A in fm^-1
  // F_q(0) = M_A^2 Q_A -- G_Q, validated on the deuteron at 25.83 (T10).
  return m_fm * m_fm * opt_.q_fm2 * c0_point(opt_, q2) *
         nucleon_fold(*this, t_gev2) * opt_.fq_scale;
}

std::string HoSpin1FF::provenance() const {
  std::ostringstream os;
  os.setf(std::ios::fmtflags(0), std::ios::floatfield);
  os.precision(6);
  os << "HoSpin1FF (design_C_tensor_rc.md sec. 2.1, PHENOMENOLOGICAL FIT, "
        "not a shell model): ";
  if (opt_.c0_shape == C0Shape::VmcFt) {
    const Li6VmcC0& v = Li6VmcC0::get();
    // The C0 EDGE has to be readable off the run banner, because it moves
    // F_c and F_q together and sigma^el_T is quadratic in the pair.
    // The RELATIVE name, not `v.path()`: `data_dir()` is machine-specific
    // and this string lands in the npz `meta` and the HepMC3 run info.
    os << "C0 shape = vmc-ft, the j0 transform of " << VMC_LI6_DENSITY
       << " (" << v.n_rows() << " non-zero rows), r-rescaled by "
       << v.r_scale() << " to <r^2>_point = " << LI6_R2_POINT_FM2
       << " fm^2, continued above q = " << LI6_VMC_C0_Q_CUT_FM
       << " fm^-1 as exp(-" << v.lambda_fm() << " (q - q_cut)) -- NO C0 ZERO "
          "at any q; (a = " << opt_.a_fm << " fm, alpha = " << opt_.alpha
       << ") are UNUSED on this edge";
  } else {
    os << "C0 shape = ho, F_point(a = " << opt_.a_fm << " fm, alpha = "
       << opt_.alpha << ")";
  }
  os << ", F_mag(q_z = " << opt_.fm_qz_fm << " fm^-1, b = "
     << opt_.fm_b_fm << " fm); MEASURED moments mu = " << opt_.mu_n
     << " mu_N, Q = " << opt_.q_fm2 << " fm^2 (NOT VMC: WS98's Q(6Li) = "
        "-0.23(9) fm^2 is 3x the measured one); Z = " << opt_.z
     << ", M_A = " << opt_.m_a_gev << " GeV (beams.cpp AME mass); "
     << (opt_.fold_nucleon ? "folded with G_E^p + G_E^n (N = Z isoscalar)"
                           : "POINT nucleons, no folding")
     << "; fq_scale = " << opt_.fq_scale
     << ", tail_tensor_scale = " << opt_.tail_tensor_scale
     << ".  UNFITTED STARTING VALUES: ";
  if (opt_.c0_shape == C0Shape::Ho) {
    os << "(a, alpha) reproduce <r^2>_point = 6.078 fm^2 and a first C0 zero "
          "at 3.10 fm^-1, and ";
  }
  os << "(q_z, b) are a two-parameter stand-in for WS98's THREE F_T "
        "landmarks -- design sec. 8.3 and Q10.  The C0 shape is a BAND: run "
        "c0_shape = ho AND vmc-ft (phase_B_numbers.md sec. B1), never one "
        "alone.";
  return os.str();
}

// ------------------------------------------------------ tabulated variant

std::shared_ptr<TabulatedSpin1FF> TabulatedSpin1FF::from_data_dir(
    const std::string& relative, HoSpin1FFOptions norm) {
  // CONVENTIONS.md: "Nothing resolves a data path against the working
  // directory."  `data_path` is $LIPOLGEN_DATA_DIR, else the compiled-in
  // prefix -- exactly what the VMC tables go through.
  const std::string path = data_path(relative);
  if (norm.m_a_gev == 0.0) {
    throw std::runtime_error(
        "TabulatedSpin1FF: `norm.m_a_gev` is 0.  F_C2 -> F_q needs M_A: the "
        "multipole convention is F_C2 = (2 sqrt(2)/3) eta_A F_q with "
        "eta_A = t/(4 M_A^2).  Fill it from Ion::mass(), never a literal.");
  }
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error(
        "TabulatedSpin1FF: cannot open '" + path + "' (relative = '" +
        relative + "' under data_dir()).  NOTE that this repository ships NO "
        "6Li elastic form-factor table: design_C_tensor_rc.md sec. 2.1 found "
        "WS98 (nucl-th/9807037) to be FIGURES ONLY, with no table anywhere in "
        "the paper, so there are no digitised numbers with provenance to "
        "ship.  Use HoSpin1FF::for_ion, or digitise WS98 Fig. 1 yourself and "
        "record the provenance in the file's 'source' column.");
  }
  auto out = std::shared_ptr<TabulatedSpin1FF>(new TabulatedSpin1FF());
  out->path_ = path;
  out->norm_ = norm;
  std::string line;
  std::size_t lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    // strip a comment and skip blanks / the header row
    const std::size_t hash = line.find('#');
    if (hash != std::string::npos) line.erase(hash);
    if (line.find_first_not_of(" \t\r\n,") == std::string::npos) continue;
    for (char& ch : line) { if (ch == ',') ch = ' '; }
    std::istringstream ls(line);
    double q = 0.0, c0 = 0.0, c2 = 0.0, m1 = 0.0;
    if (!(ls >> q >> c0 >> c2 >> m1)) continue;   // the header row lands here
    if (!(q > 0.0)) {
      throw std::runtime_error(
          "TabulatedSpin1FF: " + path + ":" + std::to_string(lineno) +
          " has q = " + std::to_string(q) +
          "; the interpolation is log-linear in q and needs q > 0");
    }
    if (!out->q_fm_.empty() && !(q > out->q_fm_.back())) {
      throw std::runtime_error(
          "TabulatedSpin1FF: " + path + ":" + std::to_string(lineno) +
          " is not strictly increasing in q");
    }
    out->q_fm_.push_back(q);
    out->c0_.push_back(c0);
    out->c2_.push_back(c2);
    out->m1_.push_back(m1);
  }
  if (out->q_fm_.size() < 2) {
    throw std::runtime_error("TabulatedSpin1FF: " + path +
                             " carries fewer than two (q, F_C0, F_C2, F_M1) "
                             "rows");
  }
  return out;
}

namespace {

/// Linear in ln q, REFUSING to extrapolate (design sec. 2.1: outside the
/// digitised range there is no measurement and a clamped end value would be a
/// silent model).
double log_interp_strict(const std::vector<double>& q, const std::vector<double>& f,
                         double qv, const std::string& path) {
  if (!(qv >= q.front()) || !(qv <= q.back())) {
    throw std::runtime_error(
        "TabulatedSpin1FF: q = " + std::to_string(qv) + " fm^-1 is outside " +
        path + "'s support [" + std::to_string(q.front()) + ", " +
        std::to_string(q.back()) + "]; the table refuses to extrapolate");
  }
  const std::size_t k = static_cast<std::size_t>(
      std::upper_bound(q.begin(), q.end(), qv) - q.begin());
  const std::size_t i = (k == 0) ? 0 : std::min(k, q.size() - 1) - 1;
  const double lo = std::log(q[i]), hi = std::log(q[i + 1]);
  const double t = (hi > lo) ? (std::log(qv) - lo) / (hi - lo) : 0.0;
  return f[i] + (f[i + 1] - f[i]) * t;
}

}  // namespace

double TabulatedSpin1FF::fc(double t_gev2) const {
  const double q = std::sqrt(std::max(t_gev2, 0.0)) / HBARC_GEV_FM;
  return log_interp_strict(q_fm_, c0_, q, path_);
}

double TabulatedSpin1FF::fm(double t_gev2) const {
  const double q = std::sqrt(std::max(t_gev2, 0.0)) / HBARC_GEV_FM;
  return log_interp_strict(q_fm_, m1_, q, path_) * norm_.tail_tensor_scale;
}

double TabulatedSpin1FF::fq(double t_gev2) const {
  // F_C2 = (2 sqrt(2)/3) eta_A F_q  =>  F_q = 3 F_C2 / (2 sqrt(2) eta_A),
  // which reproduces A(Q^2) = G_C^2 + (8/9) eta^2 G_Q^2 + (2/3) eta G_M^2
  // term by term (design sec. 2.1).
  const double q = std::sqrt(std::max(t_gev2, 0.0)) / HBARC_GEV_FM;
  const double c2 = log_interp_strict(q_fm_, c2_, q, path_);
  const double eta = t_gev2 / (4.0 * norm_.m_a_gev * norm_.m_a_gev);
  if (!(eta > 0.0)) return 0.0;
  return 3.0 * c2 / (2.0 * std::sqrt(2.0) * eta) * norm_.fq_scale;
}

std::string TabulatedSpin1FF::provenance() const {
  std::ostringstream os;
  os.precision(4);
  os << "TabulatedSpin1FF from " << path_ << " (" << q_fm_.size()
     << " rows, q = " << q_fm_.front() << " .. " << q_fm_.back()
     << " fm^-1, log-linear in q, NO extrapolation); F_C2 -> F_q through "
        "F_C2 = (2 sqrt(2)/3) eta_A F_q at M_A = " << norm_.m_a_gev
     << " GeV; fq_scale = " << norm_.fq_scale
     << ", tail_tensor_scale = " << norm_.tail_tensor_scale
     << ".  The file's own 'source' column carries the digitisation "
        "provenance -- read it before quoting anything.";
  return os.str();
}

// ------------------------------------------------------------- the model

namespace {

/// P_2(cos theta).
double p2_of(double ct) { return 0.5 * (3.0 * ct * ct - 1.0); }

/// `m_values(j)` without a per-event allocation, for the (never-built) plan
/// whose categories do not all share one J.  `RcModel::m_val_` covers the
/// normal case; this is the fallback that keeps `state_of` total.
const std::vector<double>& local_m_values(double j) {
  static thread_local double cached_j = -1.0;
  static thread_local std::vector<double> cached;
  if (cached_j != j) { cached = m_values(j); cached_j = j; }
  return cached;
}

/// The y-band a table node falls in: 0 = y < 0.5, 1 = 0.5-0.9, 2 = y > 0.9.
/// LOG ALL THREE: the tail grows like Y_+ ~ 1/(1-y), so a clipped edge hides
/// inside a small global number.
std::size_t y_band(double y) {
  if (y < 0.5) return 0;
  if (y <= 0.9) return 1;
  return 2;
}

}  // namespace

RcModel::RcModel(RcMode mode, RcOptions opt,
                 std::shared_ptr<const InclusiveSampler> dis, RunPlan plan,
                 Channel channel, Ion ion)
    : RcModel(mode, opt, std::move(dis), std::move(plan), channel,
              std::move(ion), nullptr) {}

RcModel::RcModel(RcMode mode, RcOptions opt,
                 std::shared_ptr<const InclusiveSampler> dis, RunPlan plan,
                 Channel channel, Ion ion,
                 std::shared_ptr<const TaggedModel> tagged)
    : mode_(mode),
      opt_(opt),
      dis_(std::move(dis)),
      tagged_(std::move(tagged)),
      plan_(std::move(plan)),
      channel_(channel),
      ion_(std::move(ion)) {
  if (mode_ == RcMode::Off) {
    exclusion_reason_ = "rc mode is off (the default): every weight is "
                        "exactly 1.0 and Event::rc_weights stays EMPTY";
    return;
  }
  if (opt_.tail_model != RcTailModel::TPeak &&
      opt_.tail_model != RcTailModel::TPeakPlusLL &&
      opt_.tail_model != RcTailModel::PolradFull) {
    throw std::runtime_error("RcModel: unhandled RcTailModel");
  }
  if (opt_.n_eta < 8) {
    throw std::runtime_error("RcModel: n_eta must be >= 8 (the eta_A "
                             "quadrature runs in Gauss-Legendre panels of 8)");
  }
  if (opt_.tail_model == RcTailModel::PolradFull) {
    // TWO REFUSALS THAT ARE NOT STYLE, both measured.
    //
    // (a) `long double` HAS TO BE WIDER THAN `double` HERE.  Eq. (B.3)'s
    //     a_ik for i = 5, 6 cancel to ~5 decimal digits in the collinear
    //     region (rc.cpp's `long double` note), and in plain `double` the
    //     TENSOR column of the tail moves by 14 % between two tanh-sinh step
    //     sizes that the unpolarised column agrees on to ten digits.  On a
    //     platform where the two types coincide the number would still be
    //     produced and would still be wrong, so it is refused instead.
    if (std::numeric_limits<long double>::digits <=
        std::numeric_limits<double>::digits) {
      throw std::runtime_error(
          "RcModel: RcTailModel::PolradFull needs a `long double` wider than "
          "`double`, and this platform's is not (digits = " +
          std::to_string(std::numeric_limits<long double>::digits) +
          ").  Eq. (B.3)'s a_ik cancel to ~5 decimal digits in the collinear "
          "region and the TENSOR column of the tail loses 14 % to it; the "
          "t-peak models are unaffected, so run --rc-tail-model t-peak or "
          "t-peak+ll");
    }
    // (b) `n_eta` IS THE tanh-sinh NODE COUNT PER PANEL HERE, and 8 of them
    //     do not resolve a peak of relative width 4 m_e^2/Q^2.  MEASURED
    //     (6Li config 1, x = 1e-3, y = 0.5): sigma^el_U = 0.2197428 at
    //     n_tau = 32 against 0.22086597 at 128, 256 and 512 -- 0.51 % low.
    //     At 64 it is 0.2208600, 2.7e-5 low.  64 is the floor.
    if (opt_.n_eta < 64) {
      throw std::runtime_error(
          "RcModel: RcTailModel::PolradFull needs n_eta >= 64 -- it is the "
          "tanh-sinh node count PER PANEL of the tau_A quadrature, and the "
          "s-/p-peaks it exists to resolve are 4 m_e^2/Q^2 wide.  MEASURED "
          "at 6Li config 1, x = 1e-3, y = 0.5: sigma^el_U is 0.51 % low at "
          "n_eta = 32 and 2.7e-5 low at 64, against a value stable to 1e-12 "
          "from 128 up (the default)");
    }
  }
  if (!(opt_.fq_scale >= 0.0) || !(opt_.tail_tensor_scale >= 0.0) ||
      !(opt_.qe_suppression >= 0.0) || !(opt_.qe_tensor_scale >= 0.0) ||
      !(opt_.sp_tensor_scale >= 0.0)) {
    throw std::runtime_error("RcModel: fq_scale / tail_tensor_scale / "
                             "qe_suppression / qe_tensor_scale / "
                             "sp_tensor_scale must be >= 0");
  }
  // A negative fraction would be silently squared away by the hypot in
  // `delta()` and recorded in meta as if it had meant something.
  if (!(opt_.a_transfer_frac >= 0.0)) {
    throw std::runtime_error(
        "RcModel: a_transfer_frac must be >= 0 -- it is a FRACTION of "
        "delta(x) added in quadrature (design_C_tensor_rc.md Q8), and the "
        "quadrature makes its sign meaningless");
  }
  // A KNOB THAT DID NOT RUN MAY NOT BE RECORDED AS IF IT HAD -- the rule
  // `m_lepton` is refused under.  `qe_tensor_scale` multiplies `sigma^q_U`,
  // which `with_qe_tail = false` never computes, so the pair would silently
  // record a polarised quasi-elastic price on a run that priced no
  // quasi-elastic tail at all.
  //
  // ALL THREE QUASI-ELASTIC KNOBS, since 2026-09-05, and not the polarised
  // one alone.  `qe_suppression` is a flat multiplier ON sigma^q_U and
  // `qe_kf_gev` the Fermi momentum of the Pauli factor S(q) INSIDE it, so
  // with the tail off neither is a factor of anything: measured, both are
  // bit-identical to `with_qe_tail = false` alone, and `meta` recorded them
  // as 0.5 and 0.25 with row status `read`.  That is the same defect the
  // clause was written for, one level up.  (What is NOT refused: either knob
  // at a non-default value with the tail ON -- `qe_suppression = 0`
  // degenerates the tail the same way `with_qe_tail = false` does, but it is
  // then a numeric edge of a term that RAN, which is a band row and not a
  // silent claim.)
  const RcOptions od;
  if (!opt_.with_qe_tail &&
      (opt_.qe_tensor_scale != 0.0 ||
       opt_.qe_suppression != od.qe_suppression ||
       opt_.qe_kf_gev != od.qe_kf_gev)) {
    throw std::runtime_error(
        "RcModel: qe_tensor_scale / qe_suppression / qe_kf_gev away from "
        "their defaults with with_qe_tail = false -- they scale, or sit "
        "inside, Eq. (44)'s sigma^q_U, which this run does not compute, so "
        "they would be recorded in meta as a price that was never paid.  "
        "Turn the quasi-elastic tail on, or leave all three at their "
        "defaults");
  }
  // THE SAME RULE, ONE LEVEL DOWN, FOR THE s-/p-PEAKS.  `sp_tensor_scale`
  // multiplies `TailTriple::u_sp`, the ELASTIC leading-log s+p column, which
  // `tail_sigma_at` fills ONLY under `RcTailModel::TPeakPlusLL`; under `TPeak`
  // that table is identically zero, so a non-zero scale is bit-identical to
  // the default and would be recorded in meta as a price that was never paid.
  // This is the `with_qe_tail` clause above with `u_sp` in place of
  // `sigma^q_U`, and it is refused rather than labelled for the reason stated
  // once at `KnobProvenance`: it names a VARIATION OF A PIECE THAT DID NOT
  // RUN.
  if (opt_.tail_model != RcTailModel::TPeakPlusLL &&
      opt_.sp_tensor_scale != 0.0) {
    throw std::runtime_error(
        std::string("RcModel: sp_tensor_scale != 0 needs tail_model = "
                    "TPeakPlusLL -- it is a fraction OF the leading-log "
                    "s-/p-peak column u_sp, and ") +
        (opt_.tail_model == RcTailModel::PolradFull
             // TWO REFUSALS, OPPOSITE REASONS, AND THE DISTINCTION MATTERS.
             ? "PolradFull does not compute that column because it does not "
               "NEED a stand-in: Eq. (18) carries the s- and p-peaks inside "
               "its own tau_A integral WITH their Eq. (A.4) tensor content, "
               "so the scale would double-count a term that ran"
             : "the t-peak-only tail does not compute it at all (the table "
               "is identically zero), so it would be recorded in meta as a "
               "price that was never paid") +
        " (rc.hpp, RcOptions::sp_tensor_scale)");
  }
  if (!(opt_.tail_max > 0.0)) {
    throw std::runtime_error("RcModel: tail_max must be > 0 (a Monte-Carlo "
                             "weight has to be bounded)");
  }
  // The whole A_zz programme runs at an UNPOLARISED beam (design sec. 1.1).
  // With lam_e * P_e != 0, tau = W_tensor/W would still be the rank-2
  // projection but `1 - rhobar/rho` would not, and the vector A_par term
  // would ride along in every closed form this header quotes.
  for (const SpinCategory& c : plan_.categories()) {
    if (c.lam_e != 0 && c.pe != 0.0) {
      throw std::runtime_error(
          "RcModel: --rc assumes an UNPOLARISED beam (lam_e * pe == 0); "
          "category '" + c.name + "' has lam_e = " + std::to_string(c.lam_e) +
          ", pe = " + std::to_string(c.pe) + ".  tau is the rank-2 projection "
          "and a polarised beam would price the VECTOR sector as if it were "
          "rank 2 (design_C_tensor_rc.md sec. 1.1)");
    }
  }

  // --- the per-Channel rule, TOTAL over `Channel` (design sec. 1.5) -------
  //
  // The two BOOLEANS come from `rc_band_applies` / `rc_tail_applies`, which
  // `PipelineConfig::validate()` also reads -- one definition, so the refusal
  // of a sub-knob whose piece did not run and the model that does not compute
  // that piece cannot disagree.  What stays here is the REASON, which only a
  // built model can state.
  applies_ = rc_band_applies(channel_);
  tail_applies_ = rc_tail_applies(channel_, opt_.with_tail);
  if (channel_ == Channel::CoherentLi6) {
    exclusion_reason_ =
        "coherent-6Li: its tensor dependence is entirely AZIMUTHAL "
        "(the cos 2phi coefficient of CoherentSampler), and nobody has "
        "computed RC for a phi-dependent tensor observable; its own "
        "M_X >= 1 GeV cut already excludes the elastic point.  Both weights "
        "are exactly 1.0";
  } else if (is_tagged_channel(channel_)) {
    // WHY THIS TEXT IS SO LONG.  `rc_tail == 1` here is HALF a kinematic fact
    // and half an OMISSION, and the two halves have opposite standing.  The
    // elastic half is a property of the route classification.  The
    // quasi-elastic half is a piece nobody computed, and the argument that
    // used to make it look safe -- "the tag vetoes the tail" -- is FALSE for
    // it.  Working that out is B4 (phase_B_numbers.md sec. B4); the run says
    // it because a reader of the npz would otherwise take `rc_tail == 1` for
    // a veto on the whole tail.
    exclusion_reason_ =
        std::string(channel_name(channel_)) +
        ": rc_tail == 1 exactly.  That is a FACT for the ELASTIC tail and an "
        "OMISSION for the QUASI-ELASTIC one.  FACT: e + A -> e' + gamma + A "
        "leaves the ion INTACT at x_L = 1, inside the 10-sigma "
        "beam-exclusion envelope, while the tag looks at "
        "x_L ~ A_spec/A_beam, so the tag itself vetoes it.  OMISSION: a "
        "QUASI-elastic knockout does not leave the ion intact.  It removes "
        "one nucleon; the A-1 remnant is UNBOUND (5Li -> alpha + p, "
        "5He -> alpha + n) and its alpha emerges at x_L ~ 2/3 -- the tag "
        "window itself -- and when the struck nucleon is one of the embedded "
        "deuteron's own the alpha is a TRUE SPECTATOR carrying the very "
        "n_M(k, c) this tag is built on.  So the tag neither vetoes it nor "
        "even suppresses it in the RATIO: the same spectator density and the "
        "same Roman-Pot acceptance multiply the tail and the tagged Born, "
        "and both select the 2 of 6 nucleons inside the deuteron.  "
        "rc_tail == 1 therefore drops a dilution of the same ORDER as the "
        "inclusive quasi-elastic one (22 % / 73 % / 99.9 % of the inclusive "
        "tail at x = 0.01 / 0.10 / 0.30), not a negligible one.  It is not "
        "implemented because it needs a TAGGED Born denominator and the tag "
        "acceptance folded into POLRAD Eq. (44), neither of which exists "
        "here or in the literature -- and a fourth piece, the cluster-"
        "elastic e + A -> e' + gamma + d + alpha, is not in Eq. (44)'s free-"
        "nucleon sum at all.  The BAND still applies, as an UNCITED "
        "extrapolation: no RC calculation exists for a tagged tensor "
        "asymmetry";
    if (!tagged_) {
      throw std::runtime_error(
          "RcModel: channel " + std::string(channel_name(channel_)) +
          " needs the TaggedModel constructor -- StruckClusterOptions::"
          "inclusive_b1 defaults to false, so the struck cluster's DIS kernel "
          "carries no b1 at all and the inclusive recipe would silently "
          "return tau = 0 and price nothing (design sec. 1.5.1)");
    }
  } else {
    if (!tail_applies_) {
      exclusion_reason_ =
          "with_tail = false: this is a BAND-ONLY diagnostic run and "
          "rc_tail == 1; it does NOT price the radiative-tail dilution";
    }
  }

  if (!dis_) {
    throw std::runtime_error("RcModel: the inclusive sampler is null");
  }
  s_ = dis_->s();

  // --- the (category, m) StateTables, resolved ONCE ----------------------
  // NEVER per event: `InclusiveSampler::state_tables` takes `cache_mutex_` on
  // every call and would serialise `Pipeline::for_each(sink, nthreads)`.
  // `m_val_` is hoisted out of the event loop: `m_values` allocates, and
  // `tensor_fraction` / `state_of` would otherwise call it per event.
  if (is_tagged_channel(channel_)) {
    m_val_ = m_values(tagged_->channel().j_ion);
  } else if (!plan_.categories().empty()) {
    m_val_ = m_values(plan_.categories().front().j);
  }
  if (applies_ && !is_tagged_channel(channel_)) {
    const std::vector<SpinCategory>& cats = plan_.categories();
    states_.resize(cats.size());
    for (std::size_t k = 0; k < cats.size(); ++k) {
      const std::vector<double> ms = m_values(cats[k].j);
      states_[k].reserve(ms.size());
      for (double m : ms) states_[k].push_back(&dis_->state_tables(cats[k], m));
    }
  }

  // --- the form factor ---------------------------------------------------
  if (tail_applies_) {
    ff_ = opt_.ff;
    if (!ff_) {
      HoSpin1FFOptions ho;
      ho.fq_scale = opt_.fq_scale;
      ho.tail_tensor_scale = opt_.tail_tensor_scale;
      ho.c0_shape = opt_.c0_shape;
      ff_ = HoSpin1FF::for_ion(ion_, ho);
    }
    build_tail_tables();
  }
}

const Spin1ElasticFF& RcModel::ff() const {
  if (!ff_) {
    throw std::runtime_error(
        "RcModel::ff: no elastic form factor was built -- the tail does not "
        "apply on this run (" + exclusion_reason_ + ")");
  }
  return *ff_;
}

std::string RcModel::ff_provenance() const {
  if (!ff_) {
    return "(no 6Li elastic form factor: the tail does not apply here -- " +
           exclusion_reason_ + ")";
  }
  return ff_->provenance();
}

double RcModel::delta(double x) const {
  const double d = rc_delta(x, opt_.delta_high_x, opt_.delta_low_x,
                            opt_.x_high, opt_.x_low);
  // The A = 2 -> A = 6 TRANSFER uncertainty, in quadrature (design Q8, and
  // `RcOptions::a_transfer_frac`).  Every number `rc_delta` interpolates
  // between is a DEUTERON number; nothing in the band is 6Li.  At the default
  // 0.0 this is `hypot(d, 0) == |d|`, and `rc_delta` is non-negative on every
  // configuration that can reach here (`PipelineConfig::validate()` refuses a
  // negative delta anchor, src/core/pipeline.cpp:394), so it returns EXACTLY
  // `d` -- which is what keeps the shipped band, and T3's `==` comparisons
  // against both anchors, bit for bit.
  //
  // ONE call site on purpose: this is the whole model's only use of
  // `rc_delta`, so the transfer term cannot be applied twice or skipped on a
  // path.  `rc_delta` itself stays the PUBLISHED shape and is what
  // `tests/test_rc.cpp` T3 and the Python `rc_delta` binding expose.
  return std::hypot(d, opt_.a_transfer_frac * d);
}

// ---------------------------------------------------------- the tail

double RcModel::born_pb_at(double x, double q2) const {
  // CONVENTIONS.md forbids defining a physics number twice, so the tail's
  // denominator is the LIBRARY's own Born and not a second transcription of
  // POLRAD Eq. (9): d^2 sigma/(dx dy) = x s d^2 sigma/(dx dQ^2).  T7 pins the
  // two against each other at 1e-12.
  return x * s_ * dis_->kernel().dsigma_unpol(x, q2, s_);
}

RcModel::TailTriple RcModel::tail_sigma_at(double x, double q2) const {
  TailTriple out;
  if (!ff_) return out;
  const double y = q2 / (x * s_);
  if (!(y > 0.0) || !(y < 1.0)) return out;

  // THE NUCLEAR MAP -- LiPolGen's own collider variables, NOT POLRAD's
  // fixed-target ones (review of 2026-09-03, minor "the nuclear map").
  //
  // POLRAD's `conkin` (adgh:747) defines the per-nucleon x with the FREE
  // NUCLEON mass against a target at rest, which makes its map
  // S_A = s m_A/m_p, x_A = x m_p/M_A.  This library does not: `beams.cpp`
  // carries the per-nucleon momentum p_u = p_A/A and the sampler's
  //     s = 4 E_e p_u  =>  S_A = 2 k1.p_A = A s   EXACTLY,
  //     x = Q^2/(y s)  =>  x_A = Q^2/(y S_A) = x/A   EXACTLY,
  // which is the design's own map (sec. 1.4.1) and is exact for the collider
  // kinematics the whole library is written in.  The two differ by
  // A m_p/M_A = 1.005 for 6Li; because the tail is steep in x_A through t_min
  // that is 0.26 % / 1.4 % / 6.7 % on sigma^el_U at x = 0.01 / 0.1 / 0.3.
  //
  //     t_min = 4 M_A^2 eta_min = M_A^2 x_A^2/(1 - x_A) = (x M_A/A)^2/(1 - x/A)
  //
  // and M_A/A = 0.9336 GeV for 6Li against m_p = 0.9383, so t_min is still
  // independent of A to 0.5 % -- which is all the design's "t_min ~ (x M_N)^2,
  // independent of A" ever claimed.  POLRAD's fixed-target form is kept only
  // as a documented alternative in the header.
  const double m_a = ion_.mass();
  const double a = static_cast<double>(ion_.A);
  const double s_a = a * s_;
  const double x_a = x / a;

  // PER-NUCLEON REDUCTION -- read this before changing it.
  //
  // Eq. (38) is the WHOLE-NUCLEUS d^2 sigma^el/(dx_A dy).  Two independent
  // derivations say so, and both were checked in the review of 2026-09-03:
  //
  //  (1) Weizsacker-Williams x Compton.  The t-peak IS the WW photon flux off
  //      the nucleus, dn = (Z^2 alpha/pi)(dz/z)(dt/t)(1 - t_min/t) F^2 at
  //      small z, folded with dsigma/dt^ = (2 pi alpha^2/s^^2) Y_+ and the
  //      Jacobian dz dQ^2 = x_A S_A dx_A dy.  That gives
  //        d^2 sigma/(dx_A dy) = 2 alpha^3 Z^2/(x_A^2 S_A) Y_+
  //                              INT (d eta/eta)(1 - eta_min/eta) F^2 ,
  //      which reproduces -Eq. (38) to a factor (1 - x_A) at eight (Z, x, y)
  //      points -- with NO 1/A in it.  The paper's "sigma_1^el = (1/A)
  //      d^2 sigma/dx_A dy" notation (polrad2t.tex:579) is loose.
  //  (2) POLRAD's own FORTRAN applies TWO factors to INT elu, not one:
  //      `ter = amh/amp` = m_p/M_A inside `apptai` (adgh:8607-8620) AND a
  //      division by `tara` = A in the main program (adgh:489, 497), against a
  //      Born that is itself PER NUCLEON (`f2sfun`, adgh:4330-4339).
  //
  // So, with this library's own map (dx_A/dx = 1/A):
  //
  //   want = (1/A) d^2 sigma/(dx dy)                 [ per nucleon ]
  //        = (1/A) (dx_A/dx) d^2 sigma/(dx_A dy)
  //        = (1/A) (1/A) [ Eq. (38) ]
  //        = Eq. (38) / A^2 .
  //
  // The earlier `per_nucleon = m_p/M_A` applied ONE of those two factors and
  // was 6.00x too large for 6Li (2x for the deuteron).  It rested on the
  // transcription check's sec. 4 "net difference 0.5 %" sentence, which was
  // itself wrong: the difference is A m_p/M_A = 5.97, not 0.5 %.  Both
  // documents are corrected.  The QRT below has always carried its /A, so the
  // old code's ERT : QRT ratio was internally inconsistent by a factor A.
  const double per_nucleon = 1.0 / (a * a);

  // POLRAD Eq. (18) + Appendix B + Eq. (A.4), the EXACT tau_A quadrature.
  //
  // IT TAKES THE SAME `per_nucleon` AS Eq. (38), and that is a MEASURED
  // statement rather than an assumed one: with a form factor dead at the
  // s-/p-peak vertex, so that only the t-peak survives, Eq. (18) as written
  // in `polrad_full_sigma_el` divided by Eq. (38) tends to 1 as x_A -> 0
  // (1.00230 at x_A = 0.003 unpolarised, 1.00313 tensor, on the F_m sector;
  // 1 + 1.64 x_A on a spin-0 Gaussian).  Both are therefore the WHOLE-NUCLEUS
  // d^2 sigma/(dx_A dy) in the same normalisation, and both get (1/A) x the
  // Jacobian dx_A/dx = 1/A.  The quasi-elastic partner is in NUCLEON
  // invariants and takes the same single 1/A its t-peak partner takes.
  //
  // AND IT HAS NO SEPARATE s-/p-PEAK COLUMNS.  `u_sp` and `qe_sp` exist
  // because `TPeakPlusLL` ADDS a leading-log estimate beside a t-peak that
  // does not contain one; Eq. (18) contains all three peaks in ONE tau_A
  // integral, with their tensor content, so splitting them out would be an
  // invented decomposition.  They stay identically zero here, which is also
  // why `RcOptions::sp_tensor_scale` is refused on this model -- not because
  // the s/p tensor part did not run, but because it DID.
  if (opt_.tail_model == RcTailModel::PolradFull) {
    const PolradFullPair el = polrad_full_sigma_el(*ff_, x_a, y, s_a, m_a,
                                                   opt_.n_eta, opt_.m_lepton);
    out.u = per_nucleon * el.u;
    out.t = per_nucleon * el.t;
    if (opt_.with_qe_tail) {
      out.qe = polrad_full_sigma_qe_u(ion_.Z, ion_.N(), x, y, s_, PROTON_MASS,
                                      opt_.n_eta, opt_.qe_kf_gev,
                                      opt_.m_lepton) /
               a;
    }
    return out;
  }

  out.u = per_nucleon *
          polrad_sigma_el_u(*ff_, x_a, y, s_a, m_a, opt_.n_eta);
  out.t = per_nucleon *
          polrad_sigma_el_t(*ff_, x_a, y, s_a, m_a, opt_.n_eta);
  if (opt_.with_qe_tail) {
    // The QRT is already in NUCLEON invariants (M -> m_p, x_A -> x, S_A -> s),
    // so it needs no Jacobian -- only the per-nucleon 1/A on the Z + N = A
    // nucleon sum of Eq. (44).
    out.qe = polrad_sigma_qe_u(ion_.Z, ion_.N(), x, y, s_, PROTON_MASS,
                               opt_.n_eta, opt_.qe_kf_gev) /
             a;
  }

  // --- the leading-log s- and p-peaks (RcTailModel::TPeakPlusLL) ----------
  //
  // SAME OBSERVABLE, SAME NORMALISATION, SAME REDUCTION.  `ll_peaks_spin1` is
  // the WHOLE-NUCLEUS d^2 sigma/(dx_A dy), exactly what Eq. (38) is, so it
  // takes the SAME `per_nucleon = 1/A^2` -- both the 1/A that makes it per
  // nucleon and the Jacobian dx_A/dx = 1/A.  `ll_peaks_qe` is already in
  // nucleon invariants and carries the Z + N multiplicity, exactly like
  // `polrad_sigma_qe_u`, so it takes the same single 1/A.  Getting these two
  // reductions crossed is the one way this block can be wrong by a factor 6
  // and still look plausible; T8(c') and T8(e) gate the sum against the
  // pieces.
  //
  // NO TENSOR PARTNER IS COMPUTED HERE.  See the `RcTailModel::TPeakPlusLL`
  // comment in rc.hpp: `ll_peaks_spin1` returns the unpolarised Rosenbluth
  // (A, B) only, POLRAD's Eq. (38) supplies no sigma_T at the s-/p-peak, and
  // inventing one from the leading log would be an uncited second definition.
  // SAY WHICH POLRAD: Eq. (18) + Eq. (A.4) DOES carry it and the PolradFull
  // branch a few lines above computes it, so the unqualified sentence is true
  // of Eq. (38) and FALSE of the paper -- which is why `sp_tensor_scale` is
  // refused on that model for the OPPOSITE reason it is refused on `TPeak`.
  // `u_sp`/`qe_sp` therefore go into `tail_ratio_at`'s numerator OUTSIDE the
  // (q_n/6) tensor term.
  if (opt_.tail_model == RcTailModel::TPeakPlusLL) {
    const LlPeaks el = ll_peaks_spin1(*ff_, x_a, y, s_a, m_a);
    out.u_sp = per_nucleon * (el.s + el.p);
    if (opt_.with_qe_tail) {
      const LlPeaks qe =
          ll_peaks_qe(ion_.Z, ion_.N(), x, y, s_, opt_.qe_kf_gev);
      out.qe_sp = (qe.s + qe.p) / a;
    }
  }
  return out;
}

void RcModel::build_tail_tables() {
  // --- the node grid ------------------------------------------------------
  // (ln x, ln y), not (ln x, ln Q^2): the tail carries
  // Y_+ = [1 + (1-y)^2]/(1-y), which diverges as 1/(1-y), and design
  // sec. 1.4.5's refinement rows at y in {0.9, 0.95, 0.97, 0.98, 0.985} are
  // ROWS on this axis and diagonal lines on the other.  x uses the sampler's
  // own cell EDGES (not the centres) so that an event drawn log-uniformly
  // inside the first or last cell is still inside the table.
  const LogGrid& g = dis_->grid();
  node_x_ = g.x_edges;
  if (node_x_.size() < 2) {
    throw std::runtime_error("RcModel: the sampler grid has no x edges");
  }

  const Scenario& sc = dis_->scenario();
  const double y_lo = std::max(1e-6, sc.y_min * (1.0 - 1e-9));
  const double y_want = sc.y_max * (1.0 + 1e-9);
  // FAIL HERE, NOT IN THE EVENT LOOP.  `locate` refuses to extrapolate, so a
  // scenario whose y_max exceeds the table's ceiling would throw out of
  // `Pipeline::event` on the FIRST event above it and kill the run.
  // `Scenario::validate` does not bound y_max and `in_acceptance` happily
  // accepts y = 1.0, so the guard has to be here.  The ceiling is the tail's
  // own: Y_+ = [1 + (1-y)^2]/(1-y) ~ 1/(1-y) diverges, and no table can carry
  // it past y = 0.9995 (where Y_+ is already 2000).
  if (y_want > RC_TAIL_Y_CEILING) {
    throw std::runtime_error(
        "RcModel: the scenario's y_max = " + std::to_string(sc.y_max) +
        " is above the radiative tail's y ceiling " +
        std::to_string(RC_TAIL_Y_CEILING) +
        ".  The t-peak carries Y_+ = [1 + (1-y)^2]/(1-y) ~ 1/(1-y), which "
        "diverges as y -> 1, so the tail tables stop there and "
        "`tail_ratio_at` refuses to extrapolate.  Either lower "
        "Scenario::y_max (the generator scenario's is 0.985) or run with "
        "`--rc off` / RcOptions::with_tail = false");
  }
  const double y_hi = std::min(RC_TAIL_Y_CEILING, y_want);
  if (!(y_hi > y_lo)) {
    throw std::runtime_error("RcModel: the scenario has an empty y window");
  }
  const std::size_t ny = std::max<std::size_t>(24, g.q2_c.size());
  node_y_.clear();
  node_y_.reserve(ny + 5);
  for (std::size_t i = 0; i < ny; ++i) {
    const double f = static_cast<double>(i) / static_cast<double>(ny - 1);
    node_y_.push_back(y_lo * std::pow(y_hi / y_lo, f));
  }
  for (double ye : {0.9, 0.95, 0.97, 0.98, 0.985}) {
    if (ye > y_lo && ye < y_hi) node_y_.push_back(ye);
  }
  std::sort(node_y_.begin(), node_y_.end());
  node_y_.erase(std::unique(node_y_.begin(), node_y_.end(),
                            [](double a, double b) {
                              return std::fabs(a - b) <= 1e-12 * std::fabs(a);
                            }),
                node_y_.end());

  const std::size_t nx = node_x_.size();
  const std::size_t nyv = node_y_.size();
  sigma_u_.assign(nx * nyv, 0.0);
  sigma_t_.assign(nx * nyv, 0.0);
  sigma_qe_.assign(nx * nyv, 0.0);
  sigma_u_sp_.assign(nx * nyv, 0.0);
  sigma_qe_sp_.assign(nx * nyv, 0.0);
  born_pb_.assign(nx * nyv, 0.0);

  std::size_t clipped = 0;
  std::size_t by_n[3] = {0, 0, 0};
  std::size_t by_c[3] = {0, 0, 0};
  for (std::size_t i = 0; i < nx; ++i) {
    for (std::size_t j = 0; j < nyv; ++j) {
      const double x = node_x_[i];
      const double y = node_y_[j];
      const double q2 = x * y * s_;
      const TailTriple tt = tail_sigma_at(x, q2);
      const std::size_t k = i * nyv + j;
      sigma_u_[k] = tt.u;
      sigma_t_[k] = tt.t;
      sigma_qe_[k] = tt.qe;
      sigma_u_sp_[k] = tt.u_sp;
      sigma_qe_sp_[k] = tt.qe_sp;
      born_pb_[k] = born_pb_at(x, q2);
    }
  }
  // The Born is interpolated in ln, and `dsigma_unpol` returns EXACTLY 0 at
  // x = 1 (the last grid edge): F2 -> 0 there and `dsigma_dx_dq2` clamps at
  // zero.  A ln(0) node would poison the whole last x panel, in which real
  // events DO land (the top accepted cell centre is x ~ 0.955 and the draw is
  // log-uniform up to the edge).  Continue the column geometrically instead --
  // the Born is log-linear to a few percent over one edge there, and nothing
  // physical lives at x = 1 anyway.
  for (std::size_t j = 0; j < nyv; ++j) {
    for (std::size_t i = 0; i < nx; ++i) {
      double& b = born_pb_[i * nyv + j];
      if (b > 0.0) continue;
      const double b1 = (i >= 1) ? born_pb_[(i - 1) * nyv + j] : 0.0;
      const double b2 = (i >= 2) ? born_pb_[(i - 2) * nyv + j] : 0.0;
      b = (b1 > 0.0 && b2 > 0.0) ? std::max(b1 * b1 / b2, 1e-300)
                                 : std::max(b1, 1e-300);
    }
  }
  // THE NODE-LEVEL CLIPPING STATISTIC, and what it is NOT.  It counts TABLE
  // NODES, which are not event-weighted, and it is evaluated at q_n = 0 while
  // the per-event clip in `tail_ratio_at` carries the (q_n/6) tensor term as
  // well (a ~1e-3 difference on the ratio, but a real one).  The EVENT-level
  // count is `Event::rc_clipped`'s kRcClipTail bit, which the pipeline sinks
  // reduce into `meta["rc_clipped_tail_event_fraction"]`.  Quote them
  // separately; neither substitutes for the other.
  for (std::size_t i = 0; i < nx; ++i) {
    for (std::size_t j = 0; j < nyv; ++j) {
      const std::size_t k = i * nyv + j;
      // The s+p columns are in this statistic because they are in the
      // event's ratio: under `TPeakPlusLL` the clip bites at low Q^2 in
      // places `TPeak` never reached, and a clipped fraction that ignored
      // them would understate the very thing it exists to report.  They are
      // identically zero under `TPeak`, so this line is bit-for-bit the old
      // one there.
      const double r = (sigma_u_[k] + sigma_u_sp_[k] +
                        opt_.qe_suppression * (sigma_qe_[k] + sigma_qe_sp_[k])) *
                       GEV2_TO_PB / born_pb_[k];
      const std::size_t b = y_band(node_y_[j]);
      ++by_n[b];
      if (!(r < opt_.tail_max)) { ++clipped; ++by_c[b]; }
    }
  }
  const double n_tot = static_cast<double>(nx * nyv);
  clipped_ = (n_tot > 0.0) ? static_cast<double>(clipped) / n_tot : 0.0;
  for (int b = 0; b < 3; ++b) {
    clipped_by_y_[static_cast<std::size_t>(b)] =
        by_n[b] > 0 ? static_cast<double>(by_c[b]) /
                          static_cast<double>(by_n[b])
                    : 0.0;
  }
}

double RcModel::tail_ratio_at(double x, double q2, double q_n,
                              bool* clipped) const {
  if (clipped) *clipped = false;
  if (!tail_applies_) return 0.0;
  const double y = q2 / (x * s_);
  std::size_t ix = 0, iy = 0;
  double tx = 0.0, ty = 0.0;
  locate(node_x_, x, "x", &ix, &tx);
  locate(node_y_, y, "y", &iy, &ty);
  const std::size_t nyv = node_y_.size();

  // sigma_u and sigma_qe are POSITIVE (the Eq. (18) leading minus, T8(0)) and
  // span decades, so they interpolate in ln.  sigma_t CHANGES SIGN with x
  // (check sec. 8: +0.106, +0.064, -0.117 at three deuteron points), so ln is
  // not available for it and its RATIO to sigma_u -- smooth, O(0.1) -- is
  // what is interpolated linearly.  That is an INTERPOLATION device and NOT a
  // physics scale factor: sigma_t is its own independent quadrature at every
  // node, `sigma_tail_t()` returns it, and this reconstruction is EXACT at
  // every node.
  const std::size_t idx[4] = {ix * nyv + iy, ix * nyv + iy + 1,
                              (ix + 1) * nyv + iy, (ix + 1) * nyv + iy + 1};
  double su = 0.0, sqe = 0.0, born = 0.0, ratio_t = 0.0;
  {
    const double w[4] = {(1.0 - tx) * (1.0 - ty), (1.0 - tx) * ty,
                         tx * (1.0 - ty), tx * ty};
    double lu = 0.0, lq = 0.0, lb = 0.0, r = 0.0;
    for (int k = 0; k < 4; ++k) {
      lu += w[k] * safe_log(sigma_u_[idx[k]]);
      lq += w[k] * safe_log(sigma_qe_[idx[k]]);
      lb += w[k] * safe_log(born_pb_[idx[k]]);
      const double u = sigma_u_[idx[k]];
      r += w[k] * ((u > 1e-300) ? sigma_t_[idx[k]] / u : 0.0);
    }
    su = std::exp(lu);
    sqe = std::exp(lq);
    born = std::exp(lb);
    ratio_t = r;
  }
  // THE LEADING-LOG s+p COLUMNS, and the ONE line that makes `--rc off` and
  // `tail_model = TPeak` bit for bit what they were before this branch
  // existed: under `TPeak` the two tables are identically zero and this block
  // is not entered at all, so not one floating-point operation of the shipped
  // default moved.  They are positive and span decades like `sigma_u_`, so
  // they interpolate in ln by the same rule.
  double su_sp = 0.0, sqe_sp = 0.0;
  if (opt_.tail_model == RcTailModel::TPeakPlusLL) {
    const double w[4] = {(1.0 - tx) * (1.0 - ty), (1.0 - tx) * ty,
                         tx * (1.0 - ty), tx * ty};
    double lus = 0.0, lqs = 0.0;
    for (int k = 0; k < 4; ++k) {
      lus += w[k] * safe_log(sigma_u_sp_[idx[k]]);
      lqs += w[k] * safe_log(sigma_qe_sp_[idx[k]]);
    }
    su_sp = std::exp(lus);
    sqe_sp = std::exp(lqs);
  }
  // design sec. 1.4.5's numerator, with Eq. (37)'s Q_N/6 -- dropping that 1/6
  // inflates the tensor tail sixfold and no unpolarised test sees it.
  //
  // THE s+p TERMS ARE OUTSIDE THE TENSOR TERM ON PURPOSE.  The tensor piece is
  // (q_n/6) ratio_t su with ratio_t = sigma_t/sigma_u, i.e. it is
  // (q_n/6) sigma_t -- the t-peak's OWN tensor quadrature and nothing else.
  // Writing `(q_n/6) ratio_t (su + su_sp)` instead would silently assert that
  // the s-/p-peaks carry the t-peak's tensor-to-unpolarised ratio, which is an
  // uncited claim POLRAD does not supply (rc.hpp, RcTailModel::TPeakPlusLL,
  // point (2)).  The consequence -- the tensor FRACTION of the tail falls
  // where the s-/p-peaks matter -- is a physics statement about what is
  // unknown, and is documented rather than absorbed.
  //
  // AND THE QUASI-ELASTIC TENSOR TERM IS A STAND-IN, OFF BY DEFAULT.  The
  // last term below is `RcOptions::qe_tensor_scale` -- zero unless a caller
  // asks -- and it reuses the SAME `ratio_t` the elastic tensor term uses, so
  // `qe_tensor_scale = 1` is literally "the quasi-elastic tensor fraction is
  // the elastic one".  That is a BORROWED magnitude and not a derived bound;
  // read `RcOptions::qe_tensor_scale` before quoting a number from it.  It
  // rides INSIDE `qe_suppression` because that knob is a flat multiplier on
  // the whole quasi-elastic tail, and a tensor stand-in whose unpolarised
  // parent has been switched off would be a tail with no parent.  It is
  // strictly LINEAR in the scale, which is why -- unlike `fq_scale` -- one run
  // rescales to any other value.
  //
  // AND THE s/p TENSOR TERM IS A SECOND STAND-IN, OFF BY DEFAULT AND
  // UNREACHABLE UNDER `TPeak`.  `RcOptions::sp_tensor_scale` -- zero unless a
  // caller asks, and refused unless `tail_model == TPeakPlusLL` -- reuses the
  // SAME `ratio_t`, so `sp_tensor_scale = 1` is literally "the s/p tensor
  // fraction equals the elastic t-peak's".  It sits OUTSIDE `qe_suppression`
  // because `su_sp` does: it is the ELASTIC s+p column.  It is a BOUND WITH
  // NO DERIVATION -- narrower than the quasi-elastic one (it borrows across
  // Q'^2 within ONE coherent vertex, not across coherence) but empty exactly
  // where the coherent form factor is dead; read `RcOptions::sp_tensor_scale`
  // point (3) before quoting a number from it.  Strictly LINEAR, so one run
  // rescales.  The QUASI-ELASTIC s/p column `sqe_sp` is covered by NEITHER
  // scale and stays exactly tensor-blind: point (5).
  const double num = su + (q_n / 6.0) * ratio_t * su + su_sp +
                     (q_n / 6.0) * opt_.sp_tensor_scale * ratio_t * su_sp +
                     opt_.qe_suppression *
                         (sqe + sqe_sp +
                          (q_n / 6.0) * opt_.qe_tensor_scale * ratio_t * sqe);
  const double ratio = num * GEV2_TO_PB / born;
  if (clipped) *clipped = !(ratio < opt_.tail_max);
  return std::min(ratio, opt_.tail_max);
}

double RcModel::tail_ratio(int cell, double q_n) const {
  if (!tail_applies_) return 0.0;
  if (cell < 0 || static_cast<std::size_t>(cell) >= dis_->n_cells()) {
    throw std::runtime_error("RcModel::tail_ratio: cell " +
                             std::to_string(cell) + " is outside the "
                             "sampler's accepted-cell range");
  }
  const std::size_t c = static_cast<std::size_t>(cell);
  return tail_ratio_at(dis_->x_cells()[c], dis_->q2_cells()[c], q_n);
}

// ------------------------------------------------------- the tensor part

const InclusiveSampler::StateTables* RcModel::state_of(
    const SpinLabels& sp) const {
  const std::vector<SpinCategory>& cats = plan_.categories();
  for (std::size_t k = 0; k < cats.size(); ++k) {
    const SpinCategory& c = cats[k];
    if (c.lam_e != sp.lam_e) continue;
    if (std::fabs(c.pe - sp.pe) > 1e-12) continue;
    if (std::fabs(c.j - sp.j) > 1e-12) continue;
    if (std::fabs(c.theta_s - sp.theta_s) > 1e-12) continue;
    if (std::fabs(c.phi_s - sp.phi_s) > 1e-12) continue;
    const std::vector<double>& ms =
        (m_val_.size() == states_[k].size()) ? m_val_ : local_m_values(c.j);
    for (std::size_t im = 0; im < ms.size(); ++im) {
      if (std::fabs(ms[im] - sp.m_ion) <= 1e-12) return states_[k][im];
    }
  }
  throw std::runtime_error(
      "RcModel::tensor_fraction: the event's spin state (j = " +
      std::to_string(sp.j) + ", m = " + std::to_string(sp.m_ion) +
      ", lam_e = " + std::to_string(sp.lam_e) + ", pe = " +
      std::to_string(sp.pe) + ", theta_S = " + std::to_string(sp.theta_s) +
      ", phi_S = " + std::to_string(sp.phi_s) +
      ") is not one this run plan carries.  The StateTables are resolved ONCE "
      "in the constructor (they must never be looked up per event -- that "
      "takes the sampler's cache mutex and serialises Pipeline::for_each), so "
      "an off-plan state has no table to read");
}

double RcModel::q_n_of(double m, double theta_s) const {
  // Q_N == P_zz^eff == 3 Q_NN P_2(cos theta_S): the SAME rank-2 geometry the
  // band uses, and the same `tensor_moments` the sampler's own w_avg is built
  // from -- POLRAD's Eqs. (9)/(10) -1/3 : +1/6 ratio IS P_2(0):P_2(90 deg),
  // and Eq. (43) says the same for the tail (design sec. 1.4.4, 1.4.7).
  return 3.0 * dis_->kernel().tensor_moments(m).first * p2_of(std::cos(theta_s));
}

double RcModel::tagged_tau(const Event& ev,
                           const std::vector<double>* pops) const {
  // design sec. 1.5.1: the struck cluster's DIS kernel has NO b1 when
  // `StruckClusterOptions::inclusive_b1` is false (the default), so the
  // inclusive recipe would return tau = 0 and price nothing.  The tagged
  // tensor structure lives in the M-dependence of the spectator density.
  //
  //   rho_{M,m_S}(k,c; x,Q2,phi) = n_{M,m_S}(k,c) W_DIS(m_S; x,Q2,phi)
  //
  // W_DIS cancels between rho_M and its M-average at lam_e P_e = 0 (which the
  // constructor enforces), leaving
  //
  //   tau_tag = 1 - nbar(k,c)/n_M(k,c),   nbar = (1/(2J+1)) SUM_M n_M(k,c)
  //
  // which is IDENTICALLY (P_zz/2)A_zz^wf/[1 + (P_zz/2)A_zz^wf] with
  // A_zz^wf = `azz_tensor_curve` -- the SAME closed form as the inclusive tau,
  // NOT (1/2)A_zz^wf.  T2 pins it.
  const double k = ev.kin.k;
  const double c = ev.kin.cos_theta_k;
  const std::vector<double>& ms = m_val_;
  double nbar = 0.0;
  for (double m : ms) nbar += tagged_->n_of_kc(m, k, c);
  nbar /= static_cast<double>(ms.size());
  double rho = 0.0;
  if (pops == nullptr) {
    rho = tagged_->n_of_kc(ev.spin.m_ion, k, c);
  } else {
    for (std::size_t i = 0; i < ms.size() && i < pops->size(); ++i) {
      const double p = (*pops)[i];
      if (p <= 0.0) continue;
      rho += p * tagged_->n_of_kc(ms[i], k, c);
    }
  }
  if (!(rho > 0.0)) return 0.0;
  return 1.0 - nbar / rho;
}

double RcModel::tensor_fraction(const Event& ev) const {
  if (!applies_) return 0.0;
  if (is_tagged_channel(ev.channel)) return tagged_tau(ev, nullptr);

  const InclusiveSampler::StateTables* st = state_of(ev.spin);
  const int cell = ev.kin.cell;
  if (cell < 0 || static_cast<std::size_t>(cell) >= st->w_avg.size()) {
    throw std::runtime_error(
        "RcModel::tensor_fraction: Event::kin.cell = " + std::to_string(cell) +
        " is outside the sampler's accepted-cell range");
  }
  const std::size_t c = static_cast<std::size_t>(cell);
  const double phip = ev.kin.phi - ev.spin.phi_s;
  const double cp = std::cos(phip);
  const double c2p = std::cos(2.0 * phip);
  const double w = 1.0 + st->w_avg[c] + st->a1[c] * cp + st->a2[c] * c2p;
  double wt = st->w_tensor[c] + st->a1_tensor[c] * cp + st->a2_tensor[c] * c2p;
  if (opt_.scope == RcScope::TensorAll) {
    // `a2` is entirely rank 2 (the O(gamma^2) b-sector harmonic plus Delta),
    // so `a2 - a2_tensor` is exactly the Delta cos 2phi term the default scope
    // leaves out.  NO LITERATURE SUPPORT -- for PRICING the omission only.
    wt += (st->a2[c] - st->a2_tensor[c]) * c2p;
  }
  if (!(std::fabs(w) > 1e-300)) return 0.0;
  return wt / w;
}

double RcModel::tensor_fraction(const Event& ev, std::size_t k) const {
  if (!applies_) return 0.0;
  const std::vector<SpinCategory>& cats = plan_.categories();
  if (k >= cats.size()) {
    throw std::runtime_error("RcModel::tensor_fraction: spin category " +
                             std::to_string(k) + " is outside the run plan's " +
                             std::to_string(cats.size()) + " categories");
  }
  const SpinCategory& cat = cats[k];
  if (is_tagged_channel(ev.channel)) return tagged_tau(ev, &cat.populations);

  const int cell = ev.kin.cell;
  // BOTH bounds, as the slot-0 overload above already checks and as
  // `tagged_tau` checks: this one tested `cell < 0` alone until 2026-09-16,
  // and `states_[k][im]->w_avg[c]` below then read past the StateTables for
  // any hand-built `Event` whose cell is >= n_cells -- and rc.hpp advertises
  // that a test may build the record by hand.  `weights(ev, k)` calls this
  // first, so the later indexing there is covered by this one check.
  if (cell < 0 || static_cast<std::size_t>(cell) >= dis_->n_cells()) {
    throw std::runtime_error(
        "RcModel::tensor_fraction: Event::kin.cell = " + std::to_string(cell) +
        " is outside the sampler's accepted-cell range (0 .. " +
        std::to_string(dis_->n_cells()) + ")");
  }
  const std::size_t c = static_cast<std::size_t>(cell);
  // The MIXTURE sum_m p_m W_m, mirroring `InclusiveSampler::weights_for` --
  // NOT a pure state.  It equals `tensor_fraction(ev)` only when category k's
  // population vector is pure.
  const double phip = ev.kin.phi - cat.phi_s;
  const double cp = std::cos(phip);
  const double c2p = std::cos(2.0 * phip);
  double w = 0.0, wt = 0.0;
  for (std::size_t im = 0; im < states_[k].size(); ++im) {
    const double p = (im < cat.populations.size()) ? cat.populations[im] : 0.0;
    if (p <= 0.0) continue;
    const InclusiveSampler::StateTables* st = states_[k][im];
    w += p * (1.0 + st->w_avg[c] + st->a1[c] * cp + st->a2[c] * c2p);
    double t = st->w_tensor[c] + st->a1_tensor[c] * cp + st->a2_tensor[c] * c2p;
    if (opt_.scope == RcScope::TensorAll) {
      t += (st->a2[c] - st->a2_tensor[c]) * c2p;
    }
    wt += p * t;
  }
  if (!(std::fabs(w) > 1e-300)) return 0.0;
  return wt / w;
}

// ------------------------------------------------------------ the weights

double RcModel::clamp_tau(double tau, bool* clipped) const {
  // A MONTE-CARLO WEIGHT MUST BE BOUNDED (design line 1327, applied there to
  // the tail and NOT to the band).  On the tagged channels
  // tau_tag = 1 - nbar/n_M diverges at the nodes of the M-dependent spectator
  // density, and the unclamped band published rc_tensor_hi down to -1.79 on
  // 6Li tagged-alpha and -8.35 on 7Li -- i.e. NEGATIVE event weights in the
  // npz.  n_M -> 0 is exactly where the fractional-rescale ansatz breaks: the
  // tensor part of the density cancels the unpolarised part there, so "delta
  // times the tensor fraction" is no longer a small variation of anything.
  // With band_tau_max = 1 both edges stay in [1 - delta, 1 + delta] > 0.
  const double m = opt_.band_tau_max;
  if (!(m > 0.0)) return tau;
  const double out = clip(tau, -m, m);
  if (clipped) *clipped = (out != tau);
  return out;
}

RcWeights RcModel::weights(const Event& ev) const {
  RcWeights w;
  if (!applies_) return w;
  const double d =
      delta(ev.kin.x) * clamp_tau(tensor_fraction(ev), &w.band_clipped);
  w.lo = 1.0 - d;
  w.hi = 1.0 + d;
  if (tail_applies_) {
    const InclusiveSampler::StateTables* st = state_of(ev.spin);
    const std::size_t c = static_cast<std::size_t>(ev.kin.cell);
    const double q_n = q_n_of(ev.spin.m_ion, ev.spin.theta_s);
    // design sec. 1.4.5: the denominator is the PHI-AVERAGED polarised Born,
    // Born_unpol * (1 + w_avg).  The tail itself is the phi-INTEGRATED one, so
    // at theta_S != 0 it dilutes the phi-averaged rate correctly and the
    // phi-differential one only on average (sec. 1.4.7's printed caveat).
    const double dens = 1.0 + st->w_avg[c];
    w.tail = 1.0 + tail_ratio_at(ev.kin.x, ev.kin.q2, q_n, &w.tail_clipped) /
                       std::max(dens, 1e-12);
  }
  return w;
}

RcWeights RcModel::weights(const Event& ev, std::size_t k) const {
  RcWeights w;
  if (!applies_) return w;
  const double d =
      delta(ev.kin.x) * clamp_tau(tensor_fraction(ev, k), &w.band_clipped);
  w.lo = 1.0 - d;
  w.hi = 1.0 + d;
  if (tail_applies_) {
    const SpinCategory& cat = plan_.categories()[k];
    const std::size_t c = static_cast<std::size_t>(ev.kin.cell);
    // Q_N and (1 + w_avg) are both LINEAR in the populations, so the category
    // mixture is just the population average of each.
    double q_n = 0.0, dens = 0.0;
    const std::vector<double>& ms =
        (m_val_.size() == states_[k].size()) ? m_val_ : local_m_values(cat.j);
    for (std::size_t im = 0; im < states_[k].size(); ++im) {
      const double p = (im < cat.populations.size()) ? cat.populations[im] : 0.0;
      if (p <= 0.0) continue;
      q_n += p * q_n_of(ms[im], cat.theta_s);
      dens += p * (1.0 + states_[k][im]->w_avg[c]);
    }
    w.tail = 1.0 + tail_ratio_at(ev.kin.x, ev.kin.q2, q_n, &w.tail_clipped) /
                       std::max(dens, 1e-12);
  }
  return w;
}

void RcModel::fill(Event& ev) const {
  // A no-op when the mode is Off: `Event::rc_weights` stays EMPTY and the
  // HepMC3 weight vector and the npz key set are bit-for-bit today's.
  if (mode_ == RcMode::Off) return;
  const std::size_t n_slot = 1 + ev.spin_weights.size();
  // SLOT-MAJOR: entry `slot * kRcWeightCount + i` is `rc_weight_name(i, slot)`
  // -- see the Event::rc_weights docstring.
  ev.rc_weights.assign(n_slot * kRcWeightCount, 1.0);
  const RcWeights w0 = weights(ev);
  ev.rc_weights[0] = w0.lo;
  ev.rc_weights[1] = w0.hi;
  ev.rc_weights[2] = w0.tail;
  // SLOT 0's ceilings, for the run banner and the npz meta.  The node-level
  // `clipped_cell_fraction()` is a DIFFERENT quantity (nodes are not
  // event-weighted, and it is computed at q_n = 0 while the per-event clip
  // includes the (q_n/6) tensor term -- a 1e-3 difference, but a real one).
  ev.rc_clipped = (w0.tail_clipped ? kRcClipTail : 0u) |
                  (w0.band_clipped ? kRcClipBand : 0u);
  const std::size_t n_cat = plan_.categories().size();
  for (std::size_t k = 0; k < ev.spin_weights.size(); ++k) {
    if (k >= n_cat) continue;   // degrade to 1.0 rather than throw
    const RcWeights wk = weights(ev, k);
    const std::size_t off = (1 + k) * kRcWeightCount;
    ev.rc_weights[off + 0] = wk.lo;
    ev.rc_weights[off + 1] = wk.hi;
    ev.rc_weights[off + 2] = wk.tail;
  }
}

}  // namespace lipolgen
