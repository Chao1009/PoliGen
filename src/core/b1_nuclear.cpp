#include "lipolgen/b1_nuclear.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "lipolgen/asymmetries.hpp"   // gamma_squared
#include "lipolgen/beams.hpp"         // DEUTERON()
#include "lipolgen/constants.hpp"     // HBARC_GEV_FM, M_NUCLEON, LI6_B1_*
#include "lipolgen/numerics.hpp"      // trapezoid, linspace, np_interp
#include "lipolgen/spectator.hpp"     // nuclear_mass, LI6_ALPHA_TAG
#include "lipolgen/tagged.hpp"        // VMC_*, li6_alpha_channel

namespace lipolgen {
namespace {

/// Composite Simpson on a UNIFORM grid with an odd number of points, falling
/// back to `trapezoid` otherwise.  That fallback is UNREACHABLE from
/// `LightConeDensities`, whose constructor bumps an even `n_k` to the next odd
/// value (and reports the bumped value through `options()`) precisely because
/// it would otherwise be a silent 22 % error -- see below.
///
/// DEVIATION from design 3, which says all quadrature is `numerics.hpp`'s
/// trapezoid.  It is a CORRECTNESS item, measured: b1 at small x is a
/// three-decade cancellation (int delta_T f dy = 0 exactly, so at x <~ 0.1 the
/// convolution is dominated by how the density's O(h^2) endpoint bias
/// survives that cancellation), and with the trapezoid at the design's own
/// n_k = 2001 the answer at x = 0.05 was 22 % below its n_k -> infinity limit
/// and still 1.3 % off at n_k = 8001.  Simpson on the SAME grid and the SAME
/// n_k reproduces the Richardson limit to 1e-4.  The y-grid integrals
/// (`norm`, `mean_y`, the reported int delta_T f dy) stay on `trapezoid`,
/// which is what design 5.1 G0e and 5.2 G1c are stated against.
double simpson(const std::vector<double>& y, const std::vector<double>& x) {
  const std::size_t n = y.size();
  if (n < 3 || n % 2 == 0) return trapezoid(y, x);
  const double h = (x.back() - x.front()) / static_cast<double>(n - 1);
  std::vector<double> t(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double w = (i == 0 || i + 1 == n) ? 1.0 : ((i % 2 == 1) ? 4.0 : 2.0);
    t[i] = w * y[i];
  }
  return pairwise_sum(t) * h / 3.0;
}

}  // namespace

/// Natural cubic spline through (x, y).  Written here rather than pulled from
/// `numerics.hpp` because `np_interp` is deliberately LINEAR (it must stay
/// bit-compatible with `np.interp`) and a linear u(k), w(k) biases
/// int k^2 (u^2 + w^2) dk by +2.2 % on the coarse `fdeut.av18` grid -- see
/// `ClusterPartialWave::operator()`.
class CubicSpline {
 public:
  CubicSpline(std::vector<double> x, std::vector<double> y)
      : x_(std::move(x)), y_(std::move(y)), m_(x_.size(), 0.0) {
    const std::size_t n = x_.size();
    if (n < 3) return;
    std::vector<double> a(n, 0.0), b(n, 0.0), c(n, 0.0), d(n, 0.0);
    for (std::size_t i = 1; i + 1 < n; ++i) {
      const double hm = x_[i] - x_[i - 1], hp = x_[i + 1] - x_[i];
      a[i] = hm / 6.0;
      b[i] = (hm + hp) / 3.0;
      c[i] = hp / 6.0;
      d[i] = (y_[i + 1] - y_[i]) / hp - (y_[i] - y_[i - 1]) / hm;
    }
    for (std::size_t i = 2; i + 1 < n; ++i) {   // Thomas, natural boundary
      const double f = a[i] / b[i - 1];
      b[i] -= f * c[i - 1];
      d[i] -= f * d[i - 1];
    }
    for (std::size_t i = n - 2; i >= 1; --i) {
      m_[i] = (d[i] - c[i] * m_[i + 1]) / b[i];
      if (i == 1) break;
    }
  }

  double at(double v, std::size_t hi) const {
    const double h = x_[hi] - x_[hi - 1];
    const double A = (x_[hi] - v) / h, B = (v - x_[hi - 1]) / h;
    return A * y_[hi - 1] + B * y_[hi]
           + ((A * A * A - A) * m_[hi - 1] + (B * B * B - B) * m_[hi]) * h * h
                 / 6.0;
  }

  double operator()(double v) const {
    const std::size_t n = x_.size();
    if (n < 2) return n == 1 ? y_[0] : 0.0;
    std::size_t hi = static_cast<std::size_t>(
        std::upper_bound(x_.begin(), x_.end(), v) - x_.begin());
    if (hi == 0) hi = 1;
    if (hi >= n) hi = n - 1;
    return at(v, hi);
  }

  /// One monotone pass for a sorted query grid.
  void eval_sorted(const std::vector<double>& q, std::vector<double>* out) const {
    const std::size_t n = x_.size(), m = q.size();
    out->assign(m, 0.0);
    if (n < 2) {
      if (n == 1) out->assign(m, y_[0]);
      return;
    }
    std::size_t hi = 1;
    for (std::size_t i = 0; i < m; ++i) {
      const double v = q[i];
      if (v < x_.front() || v > x_.back()) { (*out)[i] = 0.0; continue; }
      while (hi + 1 < n && x_[hi] < v) ++hi;
      (*out)[i] = at(v, hi);
    }
  }

 private:
  std::vector<double> x_, y_, m_;
};

namespace {

/// Spherical Bessel j_0 / j_2, with the series where the closed forms lose
/// their digits to cancellation.  (`cluster.cpp` has the same helper in an
/// anonymous namespace; it is not exported, and duplicating four lines of
/// arithmetic is cheaper than widening that file's API for this one call.)
double sph_bessel(int l, double v) {
  const double v2 = v * v;
  if (l == 0) {
    if (std::fabs(v) < 1e-3) return 1.0 - v2 / 6.0 + v2 * v2 / 120.0;
    return std::sin(v) / v;
  }
  if (l == 2) {
    if (std::fabs(v) < 1e-1) return v2 / 15.0 * (1.0 - v2 / 14.0 + v2 * v2 / 504.0);
    const double v3 = v2 * v;
    return (3.0 / v3 - 1.0 / v) * std::sin(v) - 3.0 / v2 * std::cos(v);
  }
  throw std::runtime_error("b1_nuclear: sph_bessel implemented for l = 0, 2");
}

/// The raw digitized CDKS theory-1 column as a TensorSF (see the header).
///
/// ZERO ABOVE THE TABLE'S OWN x_max (1.59).  `DigitizedTable::interp` is
/// `np_interp`, which CONSTANT-extrapolates, so without this the column's last
/// value 4.04e-6 would be returned as (4.04e-6 x 1.59)/x for ever -- and
/// `LightConeDensities::convolve` now integrates out to x/y = x_max_g on the
/// strength of this promise.  Below the table's first point (x = 0.01) the
/// clamp is KEPT: that is the pre-existing low-x behaviour and nothing here
/// integrates into it.
class RawCdksB1 : public TensorSF {
 public:
  double b1(double x, double, double) const override {
    const DigitizedTable& t = tables::kB1CdksQ2p5();
    if (x > t.x_max()) return 0.0;
    const double xs = std::max(x, 1e-6);
    return t.interp("xb1_theory1_sum", xs) / xs;
  }
};

/// Masses, from the ONE public accessor (`spectator.hpp::nuclear_mass`).  Not
/// `LI6_ALPHA_TAG().m_spec()`, which is the 10-keV-rounded 3.72738.
double m_deuteron() { return nuclear_mass(1, 2); }
double m_alpha() { return nuclear_mass(2, 4); }

}  // namespace

// --------------------------------------------------------------- F1 (Eq. 22)

double f1_cdks(double f2, double x, double q2, const RFunc& r_func) {
  if (!(x > 0.0)) return 0.0;
  return (1.0 + gamma_squared(x, q2)) * f2
         / (2.0 * x * (1.0 + resolve_r(r_func, x, q2)));
}

std::shared_ptr<const TensorSF> cdks_b1_raw_per_nucleon() {
  static const std::shared_ptr<const TensorSF> t = std::make_shared<RawCdksB1>();
  return t;
}

// --------------------------------------------------------- ClusterPartialWave

void ClusterPartialWave::rebuild() {
  spline = (k.size() >= 2) ? std::make_shared<const CubicSpline>(k, phi)
                           : nullptr;
}

double ClusterPartialWave::operator()(double kk) const {
  if (k.size() < 2 || kk < k.front() || kk > k.back()) return 0.0;
  if (!spline) spline = std::make_shared<const CubicSpline>(k, phi);
  return (*spline)(kk);
}

void ClusterPartialWave::eval_sorted(const std::vector<double>& kk,
                                     std::vector<double>* out) const {
  if (k.size() < 2) { out->assign(kk.size(), 0.0); return; }
  if (!spline) spline = std::make_shared<const CubicSpline>(k, phi);
  spline->eval_sorted(kk, out);
}

double ClusterPartialWave::norm2() const {
  std::vector<double> y(k.size());
  for (std::size_t i = 0; i < k.size(); ++i) y[i] = k[i] * k[i] * phi[i] * phi[i];
  return trapezoid(y, k);
}

ClusterPartialWave ClusterPartialWave::scaled(double factor) const {
  ClusterPartialWave out = *this;
  for (double& v : out.phi) v *= factor;
  out.rebuild();
  return out;
}

ClusterPartialWave ClusterPartialWave::from_vmc(const VmcRadial& v, int l) {
  if (l % 2 != 0) {
    throw std::runtime_error("ClusterPartialWave: i^L is real for even L only");
  }
  // phi_L = i^L psi_L.  ONLY that: the global phase belongs to the PAIR and
  // is already fixed by li6_vmc_waves() (see the header's contract).
  const double phase = ((l / 2) % 2 == 0) ? 1.0 : -1.0;
  ClusterPartialWave out;
  out.l = l;
  out.k = v.k();
  out.phi = v.psi();
  for (double& p : out.phi) p *= phase;
  out.provenance = v.provenance();
  out.rebuild();
  return out;
}

std::pair<ClusterPartialWave, ClusterPartialWave>
ClusterPartialWave::from_vmc_pair(const VmcRadial& s, const VmcRadial& d) {
  return {from_vmc(s, s.l()), from_vmc(d, d.l())};
}

ClusterPartialWave ClusterPartialWave::from_uw(std::vector<double> k_gev,
                                               std::vector<double> uw, int l) {
  if (l % 2 != 0) {
    throw std::runtime_error("ClusterPartialWave: i^L is real for even L only");
  }
  const double phase = ((l / 2) % 2 == 0) ? 1.0 : -1.0;
  ClusterPartialWave out;
  out.l = l;
  out.k = std::move(k_gev);
  out.phi = std::move(uw);
  for (double& p : out.phi) p *= phase;
  out.rebuild();
  return out;
}

std::pair<ClusterPartialWave, ClusterPartialWave> li6_alpha_d_partial_waves() {
  const TaggedChannel c =
      li6_alpha_channel(BETA_DEFAULT, P_D_LI6, ClusterWaveSource::VmcAV18);
  if (c.waves.size() != 2 || !c.waves[0].vmc || !c.waves[1].vmc) {
    throw std::runtime_error("li6_alpha_d_partial_waves: no VMC waves");
  }
  auto pair = ClusterPartialWave::from_vmc_pair(*c.waves[0].vmc, *c.waves[1].vmc);
  // Units.  The ANL momentum file tabulates rho_L = A_L^2/(4 pi) [fm^3] and
  // prints 4 pi int rho K^2 dK / (2 pi)^3, so int K^2 rho dK = 2 pi^2 x
  // (printed norm) in fm units; `vmc_from_momentum` returns sqrt(rho) against
  // k in GeV.  One factor puts BOTH waves in the system where
  // int k^2 (|phi_0|^2 + |phi_2|^2) dk = the file's own S + D norm.
  const double u = 1.0 / std::sqrt(2.0 * kPi * kPi * HBARC_GEV_FM * HBARC_GEV_FM
                                   * HBARC_GEV_FM);
  pair.first = pair.first.scaled(u);
  pair.second = pair.second.scaled(u);
  pair.first.provenance = "VMC alpha-d S wave (li6_ad1.momentum x li6.ad sign)";
  pair.second.provenance = "VMC alpha-d D wave (li6_ad1.momentum x li6.ad sign)";
  return pair;
}

// ---------------------------------------------------- ConvolutionKinematics

bool ConvolutionKinematics::k_range(double y, double* lo, double* hi) const {
  const double b = m_struck * (1.0 - y) - separation;
  const double disc = kappa * kappa + 2.0 * b / m_recoil;
  if (!(disc >= 0.0)) return false;
  const double r = std::sqrt(disc);
  const double a = m_recoil * std::fabs(kappa - r);
  const double z = m_recoil * (kappa + r);
  if (!(z > a)) return false;
  if (lo != nullptr) *lo = a;
  if (hi != nullptr) *hi = z;
  return true;
}

double ConvolutionKinematics::cos_star(double k, double y) const {
  if (!(k > 0.0)) return 0.0;
  return (m_struck * (1.0 - y) - separation - k * k / (2.0 * m_recoil))
         / (k * kappa);
}

double ConvolutionKinematics::y_max() const {
  return 1.0 - separation / m_struck
         + (m_recoil / m_struck) * kappa * kappa / 2.0;
}

// ------------------------------------------------------- LightConeDensities

LightConeDensities::LightConeDensities(ClusterPartialWave phi0,
                                       ClusterPartialWave phi2,
                                       ConvolutionKinematics kin, Options opt)
    : phi0_(std::move(phi0)), phi2_(std::move(phi2)), kin_(kin), opt_(opt) {
  if (phi0_.k.size() < 2) {
    throw std::runtime_error("LightConeDensities: empty phi_0 table");
  }
  // The inner k integral is composite SIMPSON, which needs an odd number of
  // points; `simpson()` falls back to the trapezoid otherwise and that is 22 %
  // low at x = 0.05 (the header's `n_k` comment, and `simpson`'s below).  A
  // user refining the grid types 2000 or 4000, so bump rather than throw --
  // and bump `opt_`, which is what `options()` returns, so the object reports
  // the quadrature that actually ran.
  if (opt_.n_k < 3) opt_.n_k = 3;
  if (opt_.n_k % 2 == 0) ++opt_.n_k;
  const double y_top = kin_.y_max();
  const double lo = opt_.y_min, mid_lo = opt_.y_mid_lo, mid_hi = opt_.y_mid_hi;
  const std::vector<double> a = linspace(lo, mid_lo, opt_.n_y_low);
  const std::vector<double> b = linspace(mid_lo, mid_hi, opt_.n_y_mid);
  const std::vector<double> c = linspace(mid_hi, y_top, opt_.n_y_high);
  y_.insert(y_.end(), a.begin(), a.end());
  y_.insert(y_.end(), b.begin() + 1, b.end());
  y_.insert(y_.end(), c.begin() + 1, c.end());

  const std::size_t n = y_.size();
  fs_.assign(n, 0.0); fd_.assign(n, 0.0); sd_.assign(n, 0.0);
  dd_.assign(n, 0.0); p2_.assign(n, 0.0); p4_.assign(n, 0.0);
  // The k table stops where the file does; spending n_k points on the zero
  // region beyond it would be pure waste (0.987 GeV for the alpha-d file,
  // 3.95 GeV for fdeut).
  const double k_table_max = std::min(phi0_.k.back(),
                                      phi2_.k.empty() ? phi0_.k.back()
                                                      : phi2_.k.back());
  std::vector<double> ks(opt_.n_k), w_s(opt_.n_k), w_d(opt_.n_k), w_sd(opt_.n_k),
      w_dd(opt_.n_k), w_p2(opt_.n_k), w_p4(opt_.n_k), sp0(opt_.n_k), sp2(opt_.n_k);
  for (std::size_t i = 0; i < n; ++i) {
    double klo = 0.0, khi = 0.0;
    if (!kin_.k_range(y_[i], &klo, &khi)) continue;
    khi = std::min(khi, k_table_max);
    if (!(khi > klo)) continue;
    ks = linspace(klo, khi, opt_.n_k);
    const double pre = y_[i] * kin_.m_struck / (2.0 * kin_.kappa);
    phi0_.eval_sorted(ks, &sp0);
    if (phi2_.k.empty()) {
      sp2.assign(opt_.n_k, 0.0);
    } else {
      phi2_.eval_sorted(ks, &sp2);
    }
    for (std::size_t j = 0; j < opt_.n_k; ++j) {
      const double k = ks[j];
      const double u = sp0[j];
      const double wv = -sp2[j];                // W(k) = -phi_2(k)
      double c2 = kin_.cos_star(k, y_[i]);
      c2 = clip(c2, -1.0, 1.0);
      const double cc = c2 * c2;
      const double p2v = 1.5 * cc - 0.5;
      const double p4v = (35.0 * cc * cc - 30.0 * cc + 3.0) / 8.0;
      w_s[j] = k * u * u;
      w_d[j] = k * wv * wv;
      // Eq. (21) in the (U, W) form: delta_T f = (3/(2 kappa)) y m
      // int k dk [U W/sqrt2 + W^2/4](3c^2 - 1), i.e. 3 x the same prefactor.
      w_sd[j] = 3.0 * k * (u * wv / std::sqrt(2.0)) * (3.0 * cc - 1.0);
      w_dd[j] = 3.0 * k * (wv * wv / 4.0) * (3.0 * cc - 1.0);
      w_p2[j] = k * wv * wv * p2v;
      w_p4[j] = k * wv * wv * p4v;
    }
    fs_[i] = pre * simpson(w_s, ks);
    fd_[i] = pre * simpson(w_d, ks);
    sd_[i] = pre * simpson(w_sd, ks);
    dd_[i] = pre * simpson(w_dd, ks);
    p2_[i] = pre * simpson(w_p2, ks);
    p4_[i] = pre * simpson(w_p4, ks);
  }

  std::vector<double> tot(n), yf(n);
  for (std::size_t i = 0; i < n; ++i) {
    tot[i] = fs_[i] + fd_[i];
    yf[i] = y_[i] * tot[i];
  }
  raw_norm_ = trapezoid(tot, y_);
  scale_ = 1.0;
  if (opt_.renormalize && raw_norm_ != 0.0) scale_ = opt_.norm_target / raw_norm_;
  if (scale_ != 1.0) {
    for (std::size_t i = 0; i < n; ++i) {
      fs_[i] *= scale_; fd_[i] *= scale_; sd_[i] *= scale_;
      dd_[i] *= scale_; p2_[i] *= scale_; p4_[i] *= scale_;
      tot[i] *= scale_; yf[i] *= scale_;
    }
  }
  norm_ = raw_norm_ * scale_;
  mean_y_ = norm_ != 0.0 ? trapezoid(yf, y_) / norm_ : 0.0;
  p_d_ = norm_ != 0.0 ? trapezoid(fd_, y_) / norm_ : 0.0;
  const double n0 = phi0_.norm2();
  const double n2 = phi2_.k.empty() ? 0.0 : phi2_.norm2();
  p_d_mom_ = (n0 + n2) != 0.0 ? n2 / (n0 + n2) : 0.0;
}

LightConeDensities::LightConeDensities(ClusterPartialWave phi0,
                                       ClusterPartialWave phi2,
                                       ConvolutionKinematics kin)
    : LightConeDensities(std::move(phi0), std::move(phi2), kin, Options()) {}

LightConeDensities LightConeDensities::delta_limit(double p_d, double norm) {
  LightConeDensities d;
  d.delta_ = true;
  d.delta_s_ = (1.0 - p_d) * norm;
  d.delta_d_ = p_d * norm;
  d.y_ = {1.0};
  d.norm_ = norm;
  d.raw_norm_ = norm;
  d.mean_y_ = 1.0;
  d.p_d_ = p_d;
  d.p_d_mom_ = p_d;
  return d;
}

double LightConeDensities::lookup(const std::vector<double>& t, double y) const {
  if (delta_) return 0.0;
  if (y_.size() < 2 || y < y_.front() || y > y_.back()) return 0.0;
  return np_interp(y, y_, t);
}

double LightConeDensities::f_s(double y) const {
  if (delta_) return delta_s_;
  return lookup(fs_, y);
}
double LightConeDensities::f_d(double y) const {
  if (delta_) return delta_d_;
  return lookup(fd_, y);
}
double LightConeDensities::f_unpol(double y) const { return f_s(y) + f_d(y); }
double LightConeDensities::delta_t_f_sd(double y) const { return lookup(sd_, y); }
double LightConeDensities::delta_t_f_dd(double y) const { return lookup(dd_, y); }
double LightConeDensities::delta_t_f(double y) const {
  return delta_t_f_sd(y) + delta_t_f_dd(y);
}
double LightConeDensities::f_d_p2(double y) const { return lookup(p2_, y); }
double LightConeDensities::f_d_p4(double y) const { return lookup(p4_, y); }

double LightConeDensities::norm() const { return norm_; }
double LightConeDensities::mean_y() const { return mean_y_; }
double LightConeDensities::p_d() const { return p_d_; }
double LightConeDensities::p_d_momentum() const { return p_d_mom_; }
const std::vector<double>& LightConeDensities::y_grid() const { return y_; }
const ConvolutionKinematics& LightConeDensities::kinematics() const { return kin_; }
const LightConeDensities::Options& LightConeDensities::options() const { return opt_; }
double LightConeDensities::renormalization() const { return scale_; }

double LightConeDensities::convolve(const std::function<double(double)>& dens_at,
                                    const std::function<double(double)>& g,
                                    double x, double x_max_g,
                                    std::size_t n) const {
  if (delta_) return dens_at(1.0) * g(x);
  if (y_.size() < 2 || n < 2) return 0.0;
  if (!(x_max_g > 0.0)) return 0.0;
  // x / x_max_g, NOT x: the support of g is g's business and x_max_g is how g
  // states it (see the header).  x_max_g = 1 gives back the old lower limit.
  const double lo = std::max(std::max(x / x_max_g, y_.front()), 1e-6);
  const double hi = y_.back();
  if (!(hi > lo)) return 0.0;
  const std::vector<double> yy = linspace(lo, hi, n);
  std::vector<double> v(n);
  for (std::size_t i = 0; i < n; ++i) {
    v[i] = dens_at(yy[i]) * g(x / yy[i]) / yy[i];
  }
  return trapezoid(v, yy);
}

// ------------------------------------------------------- DeuteronConvolutionB1

namespace {

/// The isoscalar nucleon F2 = (F2p + F2n)/2, zero outside 0 < x < 1.  ToyF2's
/// (1 - x)^n factors go NEGATIVE above x = 1, so the cut is not cosmetic --
/// this is the "support of g is g's business" clause of `convolve`.
double f2_isoscalar(const UnpolSF& s, double x, double q2) {
  if (!(x > 0.0) || x >= 1.0) return 0.0;
  return 0.5 * (s.f2p(x, q2) + s.f2n(x, q2));
}

}  // namespace

DeuteronConvolutionB1::DeuteronConvolutionB1()
    : DeuteronConvolutionB1(Options()) {}

DeuteronConvolutionB1::DeuteronConvolutionB1(Options opt) : opt_(std::move(opt)) {
  if (!opt_.unpol) opt_.unpol = std::make_shared<const ToyF2>();
  if (!opt_.r_func) {
    // CDKS use the SLAC/E143 world fit R1998, not the programme's toy R.
    opt_.r_func = [](double x, double q2) { return r1998(x, q2); };
  }
  // The wave function is the ONLY thing `Options::wave` changes: both
  // branches hand the SAME `FdeutTable` contract (GeV abscissa, GeV^-3/2
  // u, w with w = -phi_2, and the wave function's OWN binding energy) to the
  // same two lines below.
  const bool cdbonn = opt_.wave == DeuteronWaveSource::kCdBonn;
  const FdeutTable t =
      cdbonn ? cdbonn_fdeut_table(opt_.cdbonn_k_max_fm, opt_.cdbonn_dk_fm)
             : read_fdeut_k(opt_.fdeut_path);
  eps_ = t.ebind_gev;
  phi0_ = ClusterPartialWave::from_uw(t.k_gev, t.u, 0);
  phi2_ = ClusterPartialWave::from_uw(t.k_gev, t.w, 2);
  const std::string src =
      cdbonn ? "CD-Bonn (Machleidt PRC 63 (2001) 024001, Table XX)"
             : opt_.fdeut_path;
  phi0_.provenance = src + " u(k)";
  phi2_.provenance = src + " w(k)  [phi_2 = -w]";
  ConvolutionKinematics kin;
  kin.m_struck = M_NUCLEON;
  kin.m_recoil = M_NUCLEON;   // the residual system of a deuteron IS a nucleon
  kin.separation = eps_;
  kin.kappa = 1.0;
  LightConeDensities::Options q = opt_.quad;
  q.norm_target = 1.0;
  d_ = std::make_shared<const LightConeDensities>(phi0_, phi2_, kin, q);
}

const LightConeDensities& DeuteronConvolutionB1::dens(double x, double q2) const {
  if (!opt_.finite_q_delta) return *d_;
  if (qd_ && x == q_x_ && q2 == q_q2_) return *qd_;
  ConvolutionKinematics kin = d_->kinematics();
  kin.kappa = std::sqrt(1.0 + gamma_squared(x, q2));
  qd_ = std::make_shared<const LightConeDensities>(phi0_, phi2_, kin,
                                                   d_->options());
  q_x_ = x;
  q_q2_ = q2;
  return *qd_;
}

double DeuteronConvolutionB1::f1_nucleon(double x, double q2) const {
  const double f2 = f2_isoscalar(*opt_.unpol, x, q2);
  if (f2 == 0.0) return 0.0;
  if (!opt_.target_mass) {
    return f2 / (2.0 * x * (1.0 + resolve_r(opt_.r_func, x, q2)));
  }
  return f1_cdks(f2, x, q2, opt_.r_func);
}

double DeuteronConvolutionB1::b1(double x, double q2, double) const {
  const LightConeDensities& d = dens(x, q2);
  return d.convolve([&d](double y) { return d.delta_t_f(y); },
                    [&](double xp) { return f1_nucleon(xp, q2); }, x);
}

double DeuteronConvolutionB1::b1_sd(double x, double q2, double) const {
  const LightConeDensities& d = dens(x, q2);
  return d.convolve([&d](double y) { return d.delta_t_f_sd(y); },
                    [&](double xp) { return f1_nucleon(xp, q2); }, x);
}

double DeuteronConvolutionB1::b1_dd(double x, double q2, double) const {
  const LightConeDensities& d = dens(x, q2);
  return d.convolve([&d](double y) { return d.delta_t_f_dd(y); },
                    [&](double xp) { return f1_nucleon(xp, q2); }, x);
}

const LightConeDensities& DeuteronConvolutionB1::densities() const { return *d_; }
const DeuteronConvolutionB1::Options& DeuteronConvolutionB1::options() const {
  return opt_;
}

// ------------------------------------------------------------ Li6ConvolutionB1

Li6ConvolutionB1::Li6ConvolutionB1()
    : Li6ConvolutionB1(Li6ConvolutionOptions()) {}

Li6ConvolutionB1::Li6ConvolutionB1(Li6ConvolutionOptions opt)
    : opt_(std::move(opt)) {
  auto waves = li6_alpha_d_partial_waves();
  w0_ = waves.first;
  w2_ = waves.second;
  const double target = opt_.norm_target > 0.0
                            ? opt_.norm_target
                            : (opt_.use_spectroscopic_factor ? VMC_N_ALPHA_D_LI6
                                                             : 1.0);
  LightConeDensities::Options q = opt_.quad;
  q.renormalize = true;
  q.norm_target = target;
  ConvolutionKinematics kd;
  kd.m_struck = m_deuteron();
  kd.m_recoil = m_alpha();
  kd.separation = LI6_ALPHA_TAG().separation_energy;
  kd.kappa = 1.0;
  ConvolutionKinematics ka = kd;
  ka.m_struck = m_alpha();
  ka.m_recoil = m_deuteron();
  d_ = std::make_shared<const LightConeDensities>(w0_, w2_, kd, q);
  a_ = std::make_shared<const LightConeDensities>(w0_, w2_, ka, q);

  b1_d_ = opt_.deuteron_b1 ? opt_.deuteron_b1 : cdks_b1_raw_per_nucleon();
  // The support of the injected b1_d, which `convolve` needs and a `TensorSF`
  // cannot be asked for.  The default camp IS the digitized table, whose
  // strength above x = 1 is real (design 1.3); an injected backend is taken as
  // free-nucleon-like unless the caller says otherwise.
  b1_d_x_max_ = opt_.deuteron_b1_x_max > 0.0
                    ? opt_.deuteron_b1_x_max
                    : (opt_.deuteron_b1 ? 1.0 : tables::kB1CdksQ2p5().x_max());
  const std::shared_ptr<const UnpolSF> unpol =
      opt_.unpol ? opt_.unpol : std::make_shared<const ToyF2>();
  const RFunc rf = opt_.r_func;
  // CDKS Eq. (22) on F2_d/2, spelled ONCE.  `NuclearF2::f1a` is the massless
  // form and must never be called from this file.
  const auto f2d = std::make_shared<const NuclearF2>(DEUTERON(), unpol, nullptr, rf);
  f1_d_ = [f2d, rf](double x, double q2) {
    if (!(x > 0.0) || x >= 1.0) return 0.0;
    return f1_cdks(f2d->f2a(x, q2) / 2.0, x, q2, rf);
  };
  if (opt_.alpha_f1) {
    SFFunc3 af = opt_.alpha_f1;
    f1_alpha_ = [af](double x, double q2) { return af(x, q2, 0.0); };
  } else {
    // The isoscalar (F1p + F1n)/2 -- numerically identical to the deuteron's
    // because Z = N in both.  No alpha EMC effect (A9).
    f1_alpha_ = f1_d_;
  }
}

Li6ConvolutionB1::Li6ConvolutionB1(Li6ConvolutionOptions opt,
                                   LightConeDensities d_struck,
                                   LightConeDensities a_struck)
    : Li6ConvolutionB1(std::move(opt)) {
  // `finite_q_delta` rebuilds BOTH density sets per x from the FILE waves
  // (w0_, w2_), which this constructor does not replace -- so the injected
  // densities would be silently discarded at the first b1() call and T5's
  // zero-D-wave pair would grow a D wave back.  Refuse rather than lie.
  if (opt_.finite_q_delta) {
    throw std::runtime_error(
        "Li6ConvolutionB1: prepared densities cannot be rebuilt per x, so "
        "finite_q_delta must be false on the (options, d_struck, a_struck) "
        "constructor -- it would discard the densities you just passed");
  }
  d_ = std::make_shared<const LightConeDensities>(std::move(d_struck));
  a_ = std::make_shared<const LightConeDensities>(std::move(a_struck));
}

const LightConeDensities& Li6ConvolutionB1::dens_d(double x, double q2) const {
  if (!opt_.finite_q_delta) return *d_;
  if (qd_ && x == q_x_ && q2 == q_q2_) return *qd_;
  ConvolutionKinematics kd = d_->kinematics();
  ConvolutionKinematics ka = a_->kinematics();
  const double kap = std::sqrt(1.0 + gamma_squared(x, q2));
  kd.kappa = kap;
  ka.kappa = kap;
  qd_ = std::make_shared<const LightConeDensities>(w0_, w2_, kd, d_->options());
  qa_ = std::make_shared<const LightConeDensities>(w0_, w2_, ka, a_->options());
  q_x_ = x;
  q_q2_ = q2;
  return *qd_;
}

const LightConeDensities& Li6ConvolutionB1::dens_a(double x, double q2) const {
  if (!opt_.finite_q_delta) return *a_;
  dens_d(x, q2);
  return *qa_;
}

Li6ConvolutionB1::Terms Li6ConvolutionB1::terms(double x, double q2) const {
  const LightConeDensities& d = dens_d(x, q2);
  const LightConeDensities& a = dens_a(x, q2);
  const auto b1d = [&](double xp) { return b1_d_->b1(xp, q2, 0.0); };
  const auto f1d = [&](double xp) { return f1_d_(xp, q2); };
  const auto f1a = [&](double xp) { return f1_alpha_(xp, q2); };
  Terms t{};
  t.s = LI6_B1_PER_NUCLEON * opt_.w_embedded_s
        * d.convolve([&d](double y) { return d.f_s(y); }, b1d, x,
                     b1_d_x_max_);
  // w_CG = 1/10 EXACTLY: <P_zz>_{L=2} = sum CG^2 P_zz(m_S) = +1/10 (design
  // 1.4), so the D wave transfers a TENTH of the deuteron's tensor
  // polarization.  `w_cg_dwave` is a knob ON TOP of that exact factor.
  t.cg = LI6_B1_PER_NUCLEON * opt_.w_cg_dwave * 0.1
         * d.convolve([&d](double y) { return d.f_d(y); }, b1d, x,
                      b1_d_x_max_);
  t.dd_d = LI6_B1_PER_NUCLEON * opt_.w_alpha_d_dwave
           * d.convolve([&d](double y) { return d.delta_t_f(y); }, f1d, x);
  // The struck alpha carries TWICE the counting weight of the d terms:
  // 4/6 against 2/6 (design 1.6).
  t.dd_a = 2.0 * LI6_B1_PER_NUCLEON * opt_.w_alpha_d_dwave
           * a.convolve([&a](double y) { return a.delta_t_f(y); }, f1a, x);
  return t;
}

double Li6ConvolutionB1::b1(double x, double q2, double) const {
  const Terms t = terms(x, q2);
  return band_ * (t.s + t.dd_d + t.dd_a + t.cg);
}

double Li6ConvolutionB1::b1_embedded_s(double x, double q2, double) const {
  return band_ * terms(x, q2).s;
}
double Li6ConvolutionB1::b1_alpha_d_dwave_d(double x, double q2, double) const {
  return band_ * terms(x, q2).dd_d;
}
double Li6ConvolutionB1::b1_alpha_d_dwave_alpha(double x, double q2, double) const {
  return band_ * terms(x, q2).dd_a;
}
double Li6ConvolutionB1::b1_cg_dwave(double x, double q2, double) const {
  return band_ * terms(x, q2).cg;
}
double Li6ConvolutionB1::b1_alpha_d_dwave(double x, double q2, double) const {
  const Terms t = terms(x, q2);
  return band_ * (t.dd_d + t.dd_a);
}

double Li6ConvolutionB1::b1_alpha_d_sd(double x, double q2, double) const {
  const LightConeDensities& d = dens_d(x, q2);
  const LightConeDensities& a = dens_a(x, q2);
  const auto f1d = [&](double xp) { return f1_d_(xp, q2); };
  const auto f1a = [&](double xp) { return f1_alpha_(xp, q2); };
  return band_ * LI6_B1_PER_NUCLEON * opt_.w_alpha_d_dwave
         * (d.convolve([&d](double y) { return d.delta_t_f_sd(y); }, f1d, x)
            + 2.0 * a.convolve([&a](double y) { return a.delta_t_f_sd(y); }, f1a, x));
}

double Li6ConvolutionB1::b1_alpha_d_dd(double x, double q2, double) const {
  const LightConeDensities& d = dens_d(x, q2);
  const LightConeDensities& a = dens_a(x, q2);
  const auto f1d = [&](double xp) { return f1_d_(xp, q2); };
  const auto f1a = [&](double xp) { return f1_alpha_(xp, q2); };
  return band_ * LI6_B1_PER_NUCLEON * opt_.w_alpha_d_dwave
         * (d.convolve([&d](double y) { return d.delta_t_f_dd(y); }, f1d, x)
            + 2.0 * a.convolve([&a](double y) { return a.delta_t_f_dd(y); }, f1a, x));
}

double Li6ConvolutionB1::b1_cg_dwave_p2_remainder(double x, double q2,
                                                  double) const {
  // -(1/3) delta_T[sum P_zz f] on the D wave is
  //   (3/16pi)[(72/35) P4 - (4/21) P2] |phi_2|^2 + (1/40pi)|phi_2|^2,
  // i.e. in units of f_D = (1/4pi)|phi_2|^2 -> (3/4)(72/35) P4
  // - (3/4)(4/21) P2 + 1/10.  The 1/10 is term (3); these two are what the
  // angular average drops (design 1.5 truncation item 1).
  const LightConeDensities& d = dens_d(x, q2);
  const auto b1d = [&](double xp) { return b1_d_->b1(xp, q2, 0.0); };
  return band_ * LI6_B1_PER_NUCLEON * opt_.w_cg_dwave * (-1.0 / 7.0)
         * d.convolve([&d](double y) { return d.f_d_p2(y); }, b1d, x,
                      b1_d_x_max_);
}

double Li6ConvolutionB1::b1_cg_dwave_p4_remainder(double x, double q2,
                                                  double) const {
  const LightConeDensities& d = dens_d(x, q2);
  const auto b1d = [&](double xp) { return b1_d_->b1(xp, q2, 0.0); };
  return band_ * LI6_B1_PER_NUCLEON * opt_.w_cg_dwave * (54.0 / 35.0)
         * d.convolve([&d](double y) { return d.f_d_p4(y); }, b1d, x,
                      b1_d_x_max_);
}

std::shared_ptr<const TensorSF> Li6ConvolutionB1::banded(double scale) const {
  auto out = std::make_shared<Li6ConvolutionB1>(*this);
  out->band_ = scale;
  return out;
}

const LightConeDensities& Li6ConvolutionB1::densities() const { return *d_; }
const LightConeDensities& Li6ConvolutionB1::densities_alpha() const { return *a_; }
const Li6ConvolutionOptions& Li6ConvolutionB1::options() const { return opt_; }
double Li6ConvolutionB1::band_scale() const { return band_; }
double Li6ConvolutionB1::f1_deuteron(double x, double q2) const {
  return f1_d_(x, q2);
}
double Li6ConvolutionB1::f1_alpha(double x, double q2) const {
  return f1_alpha_(x, q2);
}
const std::shared_ptr<const TensorSF>& Li6ConvolutionB1::deuteron_b1() const {
  return b1_d_;
}

double Li6ConvolutionB1::close_kumano_integral(double lo, double hi, double q2,
                                               std::size_t n) const {
  const std::vector<double> xs = linspace(lo, hi, n);
  std::vector<double> v(n);
  for (std::size_t i = 0; i < n; ++i) v[i] = b1(xs[i], q2, 0.0);
  return trapezoid(v, xs);
}

// ------------------------------------------------------------------ landmarks

B1Landmarks b1_landmarks(const std::function<double(double)>& xb1,
                         const std::vector<double>& x_grid) {
  B1Landmarks out;
  const std::size_t n = x_grid.size();
  if (n < 2) return out;
  std::vector<double> v(n);
  for (std::size_t i = 0; i < n; ++i) v[i] = xb1(x_grid[i]);
  for (std::size_t i = 1; i < n; ++i) {
    const double a = v[i - 1], b = v[i];
    if (a == 0.0 || b == 0.0 || (a > 0.0) == (b > 0.0)) continue;
    out.zeros.push_back(x_grid[i - 1]
                        - a * (x_grid[i] - x_grid[i - 1]) / (b - a));
    out.zero_slope.push_back(b > a ? 1 : -1);
  }
  const std::size_t imin =
      static_cast<std::size_t>(std::min_element(v.begin(), v.end()) - v.begin());
  const std::size_t imax =
      static_cast<std::size_t>(std::max_element(v.begin(), v.end()) - v.begin());
  out.x_min = x_grid[imin];
  out.xb1_min = v[imin];
  out.x_max = x_grid[imax];
  out.xb1_max = v[imax];
  std::vector<double> b1v(n);
  for (std::size_t i = 0; i < n; ++i) b1v[i] = v[i] / x_grid[i];
  out.integral_b1 = trapezoid(b1v, x_grid);
  return out;
}

B1Landmarks b1_landmarks_of_table() {
  const DigitizedTable& t = tables::kB1CdksQ2p5();
  const std::vector<double> xs(t.x, t.x + t.n);
  const double* v = t.column("xb1_theory1_sum");
  const std::vector<double> ys(v, v + t.n);
  return b1_landmarks([&](double x) { return np_interp(x, xs, ys); }, xs);
}

// ------------------------------------------------------------- quadrupole

double alpha_d_quadrupole_fm2(const ClusterPartialWave& phi0,
                              const ClusterPartialWave& phi2, double r_max_fm,
                              std::size_t n_r, std::size_t n_k) {
  if (phi0.k.size() < 2 || phi2.k.size() < 2) return 0.0;
  const double n2 = phi0.norm2() + phi2.norm2();
  if (!(n2 > 0.0)) return 0.0;
  const double s = 1.0 / std::sqrt(n2);
  // Everything in fm / fm^-1 from here: r^2 in the operator, HBARC_GEV_FM the
  // ONE conversion (design 2.3).
  const double kmax_fm = std::min(phi0.k.back(), phi2.k.back()) / HBARC_GEV_FM;
  const std::vector<double> kf = linspace(0.0, kmax_fm, n_k);
  std::vector<double> p0(n_k), p2(n_k);
  for (std::size_t i = 0; i < n_k; ++i) {
    const double kg = kf[i] * HBARC_GEV_FM;
    p0[i] = s * phi0(kg) * std::pow(HBARC_GEV_FM, 1.5);
    p2[i] = s * phi2(kg) * std::pow(HBARC_GEV_FM, 1.5);
  }
  const std::vector<double> rr = linspace(r_max_fm / static_cast<double>(n_r),
                                          r_max_fm, n_r);
  std::vector<double> ur(n_r), wr(n_r), integ(n_r);
  const double c = std::sqrt(2.0 / kPi);
  std::vector<double> t0(n_k), t2(n_k);
  for (std::size_t i = 0; i < n_r; ++i) {
    const double r = rr[i];
    for (std::size_t j = 0; j < n_k; ++j) {
      const double k2 = kf[j] * kf[j];
      t0[j] = k2 * sph_bessel(0, kf[j] * r) * p0[j];
      t2[j] = k2 * sph_bessel(2, kf[j] * r) * p2[j];
    }
    ur[i] = r * c * trapezoid(t0, kf);           // u(r) = r R_0(r)
    wr[i] = -r * c * trapezoid(t2, kf);          // W(r) = -phi_2 -> w(r)
    integ[i] = r * r * wr[i] * (2.0 * std::sqrt(2.0) * ur[i] - wr[i]);
  }
  return trapezoid(integ, rr) / 20.0;
}

}  // namespace lipolgen
