// SPDX-License-Identifier: GPL-3.0-or-later
#include "lipolgen/cluster_config.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "lipolgen/coherent.hpp"    // gaussian_slope
#include "lipolgen/constants.hpp"   // HBARC_GEV_FM, ALPHA_EM
#include "lipolgen/fsi.hpp"         // cluster_point_a2_fm2
#include "lipolgen/numerics.hpp"    // trapezoid, linspace, clip
#include "lipolgen/rc.hpp"          // LI6_QUADRUPOLE_FM2, LI6_R2_POINT_FM2
#include "lipolgen/spin.hpp"        // clebsch_gordan
// For VMC_LI6_OVERLAP only: the path of `li6.ad` has ONE home and it is
// tagged.hpp's, so it is not retyped here (docs/CONVENTIONS.md).
#include "lipolgen/tagged.hpp"

#ifndef LIPOLGEN_VERSION
#define LIPOLGEN_VERSION "0.0.0"
#endif

namespace lipolgen {
namespace {

/// (1 - 1/4)^-1/2: the radial inflation that makes <s^2> of the RECENTRED
/// core equal <v^2> of the source, exactly and for any shape.
constexpr double kCmInflate = 1.1547005383792515;  // sqrt(4/3)

/// The channel-spin projections, in the library-wide +S ... -S order.
constexpr int kMs[3] = {1, 0, -1};

/// Cores drawn by the deterministic in-constructor Monte Carlo that measures
/// <s^2> when `min_nn_separation_fm > 0` (there is no closed form for the
/// rejection-correlated second moment).  1e5 cores put the statistical error
/// on <s^2> at ~3e-3 fm^2, i.e. ~2e-3 fm^2 on (G6)'s <r^2> -- 40x inside the
/// 5 sigma gate of a 5e4-configuration test.  ONLY paid when the hard core is
/// on: the default path never constructs an Rng in the constructor.
constexpr std::size_t kHardCoreS2Cores = 100000;

double sph_j0(double x) {
  if (std::fabs(x) < 1e-3) {
    const double x2 = x * x;
    return 1.0 - x2 / 6.0 + x2 * x2 / 120.0;
  }
  return std::sin(x) / x;
}

/// A tabulated function that is IDENTICALLY ZERO outside its own abscissa --
/// `VmcRadial::operator()`'s rule, NOT `np_interp`'s clamped extrapolation.
double table_at(const std::vector<double>& x, const std::vector<double>& f,
                double q) {
  if (x.empty()) return 0.0;
  if (q < x.front() || q > x.back()) return 0.0;
  return np_interp(q, x, f);
}

/// Trapezoid of `y` over `x` restricted to x <= x_hi (the nodes' own rule).
double trapz_upto(const std::vector<double>& y, const std::vector<double>& x,
                  double x_hi) {
  std::vector<double> xs, ys;
  xs.reserve(x.size());
  ys.reserve(x.size());
  for (std::size_t i = 0; i < x.size(); ++i) {
    if (x[i] <= x_hi) {
      xs.push_back(x[i]);
      ys.push_back(y[i]);
    }
  }
  if (xs.size() < 2) return 0.0;
  return trapezoid(ys, xs);
}

std::string read_file(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open table: " + path);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

std::vector<std::string> split_lines(const std::string& text) {
  std::vector<std::string> out;
  std::istringstream in(text);
  std::string line;
  while (std::getline(in, line)) out.push_back(line);
  return out;
}

std::vector<double> line_numbers(const std::string& line) {
  std::vector<double> out;
  const char* p = line.c_str();
  while (*p != '\0') {
    if (*p == '(' || *p == ')'
        || std::isspace(static_cast<unsigned char>(*p))) {
      ++p;
      continue;
    }
    char* end = nullptr;
    const double v = std::strtod(p, &end);
    if (end == p) return out;
    out.push_back(v);
    p = end;
  }
  return out;
}

/// The `r u du/dr w dw/dr` block of a `fdeut.*` file, plus the three header
/// numbers the tests anchor on.
///
/// FILE-LOCAL BY DESIGN NOTE.  `cluster.hpp`'s `read_fdeut_k` reads the same
/// file's K-space block and `read_anl_plain` reads the `****`-ruled plain
/// tables; neither matches this block, which is introduced by a bare
/// `r  u  du/dr  w  dw/dr` header and carries five columns.  It lives here
/// rather than beside the other ANL readers because it is the only consumer
/// (design sec. 3.3's escape hatch, taken and stated).
struct FdeutRTable {
  std::vector<double> r, u, w;
  double dstate = 0.0, qm = 0.0, rd = 0.0;   ///< the header's own values
};

FdeutRTable read_fdeut_r(const std::string& path) {
  const std::vector<std::string> lines = split_lines(read_file(path));
  FdeutRTable t;
  std::size_t hdr = lines.size(), ebind_hdr = lines.size();
  for (std::size_t i = 0; i < lines.size(); ++i) {
    std::istringstream probe(lines[i]);
    std::string a, b;
    probe >> a >> b;
    if (hdr == lines.size() && a == "r" && b == "u"
        && lines[i].find("du/dr") != std::string::npos) {
      hdr = i;
    }
    if (ebind_hdr == lines.size() && a == "ebind") ebind_hdr = i;
  }
  if (hdr == lines.size()) {
    throw std::runtime_error("no `r u du/dr w dw/dr` block in " + path);
  }
  if (ebind_hdr != lines.size()) {
    for (std::size_t i = ebind_hdr + 1; i < lines.size(); ++i) {
      const std::vector<double> v = line_numbers(lines[i]);
      if (v.size() < 7) continue;
      t.dstate = v[1];
      t.qm = v[2];
      t.rd = v[6];
      break;
    }
  }
  for (std::size_t i = hdr + 1; i < lines.size(); ++i) {
    const std::vector<double> v = line_numbers(lines[i]);
    if (v.size() != 5) {
      if (!t.r.empty()) break;
      continue;
    }
    t.r.push_back(v[0]);
    t.u.push_back(v[1]);
    t.w.push_back(v[3]);
  }
  if (t.r.size() < 2) throw std::runtime_error("empty fdeut r block: " + path);
  return t;
}

/// Tricomi's U(a, b, z) for a > 0, z > 0, by
/// U = 1/Gamma(a) integral_0^inf e^{-z t} t^{a-1} (1+t)^{b-a-1} dt.
///
/// Two composite-Simpson branches, because ONE substitution cannot serve both
/// regimes.  For a < 2 the integrand's t^{a-1} endpoint behaviour has a
/// singular derivative, and the substitution t = y^{1/a} removes it exactly
/// (t^{a-1} dt = dy/a).  For a >= 2 that same substitution COMPRESSES the
/// peak at t = (a-1)/z into a handful of y cells -- it put U(3.3, 6, 3.69)
/// 20 % high -- so the plain variable is used there instead.  The upper limit
/// is where the log-integrand has fallen 50 below its peak.
double tricomi_u(double a, double b, double z) {
  if (!(a > 0.0) || !(z > 0.0)) {
    throw std::runtime_error("tricomi_u: needs a > 0 and z > 0");
  }
  auto logf = [&](double t) {
    return (a - 1.0) * std::log(t) + (b - a - 1.0) * std::log1p(t) - z * t;
  };
  const double t_peak = std::fmax((a - 1.0) / z, 1e-3);
  double t_max = std::fmax(1.0, 2.0 * t_peak);
  const double lpk = logf(t_peak);
  while (logf(t_max) - lpk > -50.0 && t_max < 1e6) t_max *= 1.5;
  const std::size_t n = 20000;   // even
  double acc = 0.0;
  if (a < 2.0) {
    const double y_max = std::pow(t_max, a);
    const double h = y_max / static_cast<double>(n);
    auto g = [&](double y) {
      const double t = std::pow(y, 1.0 / a);
      return std::exp(-z * t) * std::pow(1.0 + t, b - a - 1.0);
    };
    acc = g(0.0) + g(y_max);
    for (std::size_t i = 1; i < n; ++i) {
      acc += (i % 2 == 1 ? 4.0 : 2.0) * g(static_cast<double>(i) * h);
    }
    return acc * h / 3.0 / (a * std::tgamma(a));
  }
  const double h = t_max / static_cast<double>(n);
  auto g = [&](double t) { return t <= 0.0 ? 0.0 : std::exp(logf(t)); };
  acc = g(0.0) + g(t_max);
  for (std::size_t i = 1; i < n; ++i) {
    acc += (i % 2 == 1 ? 4.0 : 2.0) * g(static_cast<double>(i) * h);
  }
  return acc * h / 3.0 / std::tgamma(a);
}

/// W_{-eta, L+1/2}(z) = e^{-z/2} z^{L+1} U(L+1+eta, 2L+2, z).
double whittaker_w(int l, double eta, double z) {
  return std::exp(-0.5 * z) * std::pow(z, l + 1)
         * tricomi_u(l + 1 + eta, 2 * l + 2, z);
}

}  // namespace

// ------------------------------------------------- the master density (G1)

double cluster_amp2(double f0, double f2, int m, int m_s, double c) {
  const int ml = m - m_s;
  double amp = 0.0;
  if (ml == 0) {
    amp += clebsch_gordan(0.0, 0.0, 1.0, m_s, 1.0, m) * f0 * theta_lm(0, 0, c);
  }
  if (std::abs(ml) <= 2) {
    amp += clebsch_gordan(2.0, ml, 1.0, m_s, 1.0, m) * f2 * theta_lm(2, ml, c);
  }
  return amp * amp;
}

double cluster_density(double f0, double f2, int m, double c) {
  double s = 0.0;
  for (int i = 0; i < 3; ++i) s += cluster_amp2(f0, f2, m, kMs[i], c);
  return s;
}

RadialMoments radial_moments(const std::vector<double>& x,
                             const std::vector<double>& f0,
                             const std::vector<double>& f2) {
  if (x.size() != f0.size() || x.size() != f2.size() || x.size() < 2) {
    throw std::runtime_error("radial_moments: mis-sized radial pair");
  }
  const std::size_t n = x.size();
  std::vector<double> a(n), b(n), r(n), qi(n), qd(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double x2 = x[i] * x[i], x4 = x2 * x2;
    a[i] = f0[i] * f0[i] * x2;
    b[i] = f2[i] * f2[i] * x2;
    r[i] = (f0[i] * f0[i] + f2[i] * f2[i]) * x4;
    qi[i] = 2.0 * std::sqrt(2.0) * f0[i] * f2[i] * x4;
    qd[i] = -f2[i] * f2[i] * x4;
  }
  RadialMoments out;
  out.n0 = trapezoid(a, x);
  out.n2 = trapezoid(b, x);
  out.r2 = trapezoid(r, x);
  out.q_int = 0.2 * trapezoid(qi, x);
  out.q_dd = 0.2 * trapezoid(qd, x);
  return out;
}

// ------------------------------------------------ the quadrupole -> a_2 map

double a2_from_quadrupole(double q_matter_fm2, int a, double t_abs, int m) {
  if (a <= 0) throw std::runtime_error("a2_from_quadrupole: A must be > 0");
  const double delta_fm2 = q_matter_fm2 / (2.0 * a);
  const double delta_gev2 =
      (m == 0 ? -2.0 : 1.0) * delta_fm2 / (HBARC_GEV_FM * HBARC_GEV_FM);
  return -0.25 * delta_gev2 * t_abs;
}

double quadrupole_from_a2_slope(double a2_over_t, int a) {
  if (a <= 0) throw std::runtime_error("quadrupole_from_a2_slope: A must be > 0");
  return -8.0 * static_cast<double>(a) * HBARC_GEV_FM * HBARC_GEV_FM
         * a2_over_t;
}

// ----------------------------------------------------------------- options

void ClusterConfigOptions::validate() const {
  if (exact_coherence) {
    throw std::runtime_error(
        "ClusterConfigOptions::exact_coherence is reserved and not "
        "implemented: the exact 5-D joint density with the m_S coherences "
        "restored needs Metropolis, not inverse-CDF (design sec. 2.6)");
  }
  if (n_r < 8 || n_c < 8) {
    throw std::runtime_error("ClusterConfigOptions: n_r and n_c must be >= 8");
  }
  if (!(r_max_fm > 0.0) || !(rnp_max_fm > 0.0)) {
    throw std::runtime_error("ClusterConfigOptions: r_max_fm and rnp_max_fm "
                             "must be positive");
  }
  if (!(alpha_d_scale > 0.0)) {
    throw std::runtime_error("ClusterConfigOptions: alpha_d_scale must be "
                             "positive");
  }
  if (min_nn_separation_fm < 0.0) {
    throw std::runtime_error("ClusterConfigOptions: min_nn_separation_fm must "
                             "not be negative");
  }
}

// -------------------------------------------------------------- the sampler

namespace {

/// The 9 (M, m_S) slots, M and m_S both in {+1, 0, -1}.
std::size_t slot(int m, int i_ms) {
  return static_cast<std::size_t>((1 - m) * 3 + i_ms);
}

}  // namespace

ClusterConfigSampler::ClusterConfigSampler(ClusterConfigOptions opt)
    : opt_(opt) {
  opt_.validate();
  build_tables();
}

void ClusterConfigSampler::build_tables() {
  std::ostringstream prov;

  // ---------------------------------------------------------- deuteron u, w
  const std::string fdeut = data_path(VMC_DEUTERON_WAVE);
  const FdeutRTable d = read_fdeut_r(fdeut);
  np_x_ = d.r;
  np_f0_.resize(np_x_.size());
  np_f2_.resize(np_x_.size());
  for (std::size_t i = 0; i < np_x_.size(); ++i) {
    np_f0_[i] = d.u[i] / np_x_[i];
    np_f2_[i] = d.w[i] / np_x_[i];
  }
  np_m_ = radial_moments(np_x_, np_f0_, np_f2_);
  r2_np_ = np_m_.r2 / np_m_.norm();
  q_d_fm2_ = 0.25 * np_m_.quadrupole() / np_m_.norm();
  prov << "deuteron u,w: " << fdeut << "  (norm " << np_m_.norm()
       << ", P_D " << np_m_.p_d() << ", Q_d " << q_d_fm2_ << " fm^2, r_d "
       << 0.5 * std::sqrt(r2_np_) << " fm)\n";

  // ------------------------------------------------------ alpha-d R_0, R_2
  const std::string ov_path = data_path(VMC_LI6_OVERLAP);
  const std::string fit_path = data_path(VMC_LI6_AD_FIT);
  std::vector<double> ax, a0, a2;
  if (opt_.alpha_d_source == AlphaDSource::OverlapRaw) {
    const std::vector<AnlTable> t = read_anl_overlap(ov_path);
    ax = t[1].x;
    a0 = t[1].col[0];
    a2 = t[1].col[1];
    prov << "alpha-d R_0,R_2: " << ov_path << " r-block, verbatim\n";
  } else {
    const AnlTable f = read_anl_plain(fit_path, 2);
    ax = f.x;
    a0 = f.col[0];
    a2 = f.col[1];
    if (opt_.alpha_d_source == AlphaDSource::FitRescaled) {
      // Each fitted wave is rescaled to the RAW block's own N_L, so the shape
      // is the smoothed Forest et al. fit and the normalization (hence P_D)
      // is li6.ad's: s_0 = 1.0128, s_2 = 0.8984.
      const std::vector<AnlTable> t = read_anl_overlap(ov_path);
      const RadialMoments rm = radial_moments(t[1].x, t[1].col[0], t[1].col[1]);
      const RadialMoments fm = radial_moments(ax, a0, a2);
      const double s0 = std::sqrt(rm.n0 / fm.n0);
      const double s2 = std::sqrt(rm.n2 / fm.n2);
      for (double& v : a0) v *= s0;
      for (double& v : a2) v *= s2;
      prov << "alpha-d R_0,R_2: " << fit_path << " shapes rescaled to "
           << ov_path << " norms (s_0 " << s0 << ", s_2 " << s2 << ")\n";
    } else {
      prov << "alpha-d R_0,R_2: " << fit_path << ", as published\n";
    }
  }
  {
    const RadialMoments rm = radial_moments(ax, a0, a2);
    s_alpha_d_ = rm.norm();
    const double inv = 1.0 / std::sqrt(s_alpha_d_);
    for (double& v : a0) v *= inv;
    for (double& v : a2) v *= inv;
  }
  // Global phase: unobservable, so fix it ONCE by convention to
  // R_0(R -> infinity) > 0.  The outer region is used rather than a single
  // last point because the raw block's tail is Monte Carlo noise.
  {
    std::vector<double> w(ax.size());
    for (std::size_t i = 0; i < ax.size(); ++i) w[i] = a0[i] * ax[i] * ax[i];
    double tail = 0.0;
    for (std::size_t i = 1; i < ax.size(); ++i) {
      if (ax[i - 1] >= 4.0) tail += 0.5 * (w[i] + w[i - 1]) * (ax[i] - ax[i - 1]);
    }
    if (tail < 0.0) {
      for (double& v : a0) v = -v;
      for (double& v : a2) v = -v;
    }
  }
  ad_x_ = ax;
  ad_f0_ = a0;
  ad_f2_ = a2;
  ad_m_ = radial_moments(ad_x_, ad_f0_, ad_f2_);
  ad_m_base_ = ad_m_;

  // ------------------------------------------- the (G9) quadrupole dial
  const double lam2 = opt_.alpha_d_scale * opt_.alpha_d_scale;
  if (opt_.quadrupole_target_fm2 != 0.0) {
    const double qi = lam2 * ad_m_.q_int, qd = lam2 * ad_m_.q_dd;
    const double pd = ad_m_.p_d();
    auto q_charge = [&](double s) {
      const double n = 1.0 - (1.0 - s * s) * pd;
      const double q = (s * qi + s * s * qd) / n;
      const double pdp = s * s * pd / n;
      return 0.5 * ((4.0 / 3.0) * q + 2.0 * q_d_fm2_ * (1.0 - 0.9 * pdp));
    };
    const double q_hi = q_charge(0.0), q_lo = q_charge(1.0);
    const double t = opt_.quadrupole_target_fm2;
    if (t > std::fmax(q_lo, q_hi) || t < std::fmin(q_lo, q_hi)) {
      std::ostringstream e;
      e << "quadrupole_target_fm2 = " << t << " fm^2 is unreachable: the "
        << "s in [0, 1] range of this source is [" << std::fmin(q_lo, q_hi)
        << ", " << std::fmax(q_lo, q_hi) << "] fm^2 (s = 0 leaves the "
        << "deuteron's own +" << q_d_fm2_ << " fm^2 behind)";
      throw std::runtime_error(e.str());
    }
    double lo = 0.0, hi = 1.0;
    for (int it = 0; it < 200; ++it) {
      const double mid = 0.5 * (lo + hi);
      if ((q_charge(mid) > t) == (q_hi > t)) {
        lo = mid;
      } else {
        hi = mid;
      }
    }
    dial_s_ = 0.5 * (lo + hi);
    // Apply it: R_2 -> s R_2, both waves renormalized by 1/sqrt(n(s)).
    const double n = 1.0 - (1.0 - dial_s_ * dial_s_) * pd;
    const double inv = 1.0 / std::sqrt(n);
    for (double& v : ad_f2_) v *= dial_s_ * inv;
    for (double& v : ad_f0_) v *= inv;
    ad_m_ = radial_moments(ad_x_, ad_f0_, ad_f2_);
    prov << "quadrupole dial: Q_charge target " << t << " fm^2 -> (G9) root s "
         << dial_s_ << ", P_D(alpha-d) now " << ad_m_.p_d()
         << " -- A DEFORMATION DIAL, NOT A WAVE FUNCTION\n";
  }
  if (opt_.alpha_d_scale != 1.0) {
    prov << "alpha-d separation scaled by " << opt_.alpha_d_scale << "\n";
  }

  // ------------------------------------------------------- the alpha core
  if (opt_.alpha_source == AlphaCoreSource::VmcHe4Density) {
    const std::string path = data_path(VMC_HE4_DENSITY);
    const std::vector<AnlTable> t = read_anl_momentum(path);
    core_x_ = t[0].x;
    core_rho_ = t[0].col[0];
    prov << "alpha core: " << path << " (VMC AV18+UX one-body point-nucleon "
         << "density, rho_nucleon = RHORP / Z)\n";
  } else {
    // rho ~ exp(-r^2 / 2 a^2) with the SINGLE copy of a^2 the FSI weight uses.
    const double a2c = cluster_point_a2_fm2(2, 4);
    core_x_ = linspace(0.0, 10.0, 401);
    core_rho_.resize(core_x_.size());
    for (std::size_t i = 0; i < core_x_.size(); ++i) {
      core_rho_[i] = std::exp(-core_x_[i] * core_x_[i] / (2.0 * a2c));
    }
    prov << "alpha core: Gaussian, a^2 = cluster_point_a2_fm2(2, 4) = " << a2c
         << " fm^2 (fsi.hpp)\n";
  }
  // Inflate BEFORE normalizing: rho_src(r) = rho(r / lambda) / lambda^3.
  if (opt_.alpha_cm_inflate) {
    for (double& v : core_x_) v *= kCmInflate;
    const double inv3 = 1.0 / (kCmInflate * kCmInflate * kCmInflate);
    for (double& v : core_rho_) v *= inv3;
  }
  {
    std::vector<double> w(core_x_.size());
    for (std::size_t i = 0; i < core_x_.size(); ++i) {
      w[i] = core_rho_[i] * core_x_[i] * core_x_[i];
    }
    const double n = 4.0 * kPi * trapezoid(w, core_x_);
    if (!(n > 0.0)) throw std::runtime_error("alpha core density has no norm");
    for (double& v : core_rho_) v /= n;
    for (std::size_t i = 0; i < core_x_.size(); ++i) {
      w[i] = core_rho_[i] * core_x_[i] * core_x_[i] * core_x_[i] * core_x_[i];
    }
    // Recentring gives <s^2> = (3/4)<v^2> EXACTLY, for any source shape.
    s2_alpha_ = 0.75 * 4.0 * kPi * trapezoid(w, core_x_);
  }

  // ------------------------------------ norm loss against the grid ceilings
  auto norm_loss = [](const std::vector<double>& x,
                      const std::vector<double>& w, double ceiling) {
    const double all = trapezoid(w, x);
    if (!(std::fabs(all) > 0.0)) return 0.0;
    return 1.0 - trapz_upto(w, x, ceiling) / all;
  };
  {
    std::vector<double> w(ad_x_.size());
    for (std::size_t i = 0; i < ad_x_.size(); ++i) {
      w[i] = (ad_f0_[i] * ad_f0_[i] + ad_f2_[i] * ad_f2_[i]) * ad_x_[i]
             * ad_x_[i];
    }
    const double loss = norm_loss(ad_x_, w, opt_.r_max_fm / opt_.alpha_d_scale);
    if (loss > 1e-3) {
      throw std::runtime_error("r_max_fm drops " + std::to_string(loss)
                               + " of the alpha-d norm (limit 1e-3)");
    }
  }
  {
    std::vector<double> w(np_x_.size());
    for (std::size_t i = 0; i < np_x_.size(); ++i) {
      w[i] = (np_f0_[i] * np_f0_[i] + np_f2_[i] * np_f2_[i]) * np_x_[i]
             * np_x_[i];
    }
    const double loss = norm_loss(np_x_, w, opt_.rnp_max_fm);
    if (loss > 1e-3) {
      throw std::runtime_error("rnp_max_fm drops " + std::to_string(loss)
                               + " of the p-n norm (limit 1e-3)");
    }
  }
  {
    std::vector<double> w(core_x_.size());
    for (std::size_t i = 0; i < core_x_.size(); ++i) {
      w[i] = core_rho_[i] * core_x_[i] * core_x_[i];
    }
    const double loss = norm_loss(core_x_, w, opt_.r_max_fm);
    if (loss > 1e-3) {
      throw std::runtime_error("r_max_fm drops " + std::to_string(loss)
                               + " of the alpha core norm (limit 1e-3)");
    }
  }

  // ------------------------------------------------------------- the grids
  auto make_grid = [](double hi, std::size_t n) {
    Grid g;
    g.dx = hi / static_cast<double>(n);
    g.x.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
      g.x[i] = (static_cast<double>(i) + 0.5) * g.dx;
    }
    return g;
  };
  gr_ad_ = make_grid(std::fmin(opt_.r_max_fm / opt_.alpha_d_scale,
                               ad_x_.back()), opt_.n_r);
  gr_np_ = make_grid(std::fmin(opt_.rnp_max_fm, np_x_.back()), opt_.n_r);
  gr_core_ = make_grid(std::fmin(opt_.r_max_fm, core_x_.back()), opt_.n_r);
  gr_c_.dx = 2.0 / static_cast<double>(opt_.n_c);
  gr_c_.x.resize(opt_.n_c);
  for (std::size_t i = 0; i < opt_.n_c; ++i) {
    gr_c_.x[i] = -1.0 + (static_cast<double>(i) + 0.5) * gr_c_.dx;
  }

  const std::size_t nr = opt_.n_r, nc = opt_.n_c;
  // core radial CDF, weight 4 pi rho v^2
  core_cdf_.assign(nr, 0.0);
  {
    double acc = 0.0;
    for (std::size_t i = 0; i < nr; ++i) {
      const double v = gr_core_.x[i];
      acc += table_at(core_x_, core_rho_, v) * v * v;
      core_cdf_[i] = acc;
    }
    if (!(acc > 0.0)) throw std::runtime_error("empty alpha core CDF");
    for (double& v : core_cdf_) v /= acc;
  }

  // ------------------------------- <s^2> when the hard core is switched on
  //
  // `s2_alpha_` above is the CLOSED FORM (3/4)<v^2>, which is exact for any
  // source shape but ONLY for INDEPENDENT draws.  `min_nn_separation_fm > 0`
  // rejects and redraws, so the four s_i are correlated and the closed form
  // is wrong -- measurably: at 0.9 fm the sampled <r^2> is 6.586 fm^2 against
  // the closed form's 6.445, +2 %, which is 8.7 sigma of a 5e4-configuration
  // set.  Everything downstream of `s2_alpha_` (`r2_analytic_fm2`, (G6),
  // `match_li6_radius`, the sidecar's "moments") would inherit that error and
  // stamp it as exact, so the correlated second moment is MEASURED here
  // instead, once, on a FIXED stream and a fixed count -- the sampler stays
  // immutable, thread-safe and reproducible, and the default path pays
  // nothing.  `rho_alpha_recentred()` is NOT corrected: it is the
  // independent-draw closed form and its header says so.
  if (opt_.min_nn_separation_fm > 0.0) {
    const double s2_free = s2_alpha_;
    Rng rng(0, 0, kConfigStream, 0);
    std::array<std::array<double, 3>, 4> s{};
    double acc = 0.0;
    for (std::size_t i = 0; i < kHardCoreS2Cores; ++i) {
      draw_alpha(rng, s);
      for (int k = 0; k < 4; ++k) {
        acc += s[k][0] * s[k][0] + s[k][1] * s[k][1] + s[k][2] * s[k][2];
      }
    }
    s2_alpha_ = acc / (4.0 * static_cast<double>(kHardCoreS2Cores));
    prov << "alpha <s^2> NUMERICAL (hard core): " << s2_alpha_
         << " fm^2 from " << kHardCoreS2Cores << " cores on the fixed stream "
         << "Rng(0, 0, CONFIG, 0), against the independent-draw closed form "
         << s2_free << " fm^2 -- (G6)'s <r^2>, match_li6_radius() and the "
         << "sidecar's moments are APPROXIMATE at this precision, and "
         << "rho_alpha_recentred() is still the INDEPENDENT-draw density\n";
  }

  // alpha-d: one (R, c) table per (M, m_S), exactly as TaggedModel::build_amp2
  ad_cdf_.assign(9, {});
  ad_weight_.assign(9, 0.0);
  std::vector<double> f0g(nr), f2g(nr);
  for (std::size_t i = 0; i < nr; ++i) {
    f0g[i] = table_at(ad_x_, ad_f0_, gr_ad_.x[i]);
    f2g[i] = table_at(ad_x_, ad_f2_, gr_ad_.x[i]);
  }
  for (int m = 1; m >= -1; --m) {
    for (int i_ms = 0; i_ms < 3; ++i_ms) {
      const int m_s = kMs[i_ms];
      std::vector<double> cdf(nr * nc);
      double acc = 0.0;
      for (std::size_t i = 0; i < nr; ++i) {
        const double r2 = gr_ad_.x[i] * gr_ad_.x[i];
        for (std::size_t j = 0; j < nc; ++j) {
          acc += cluster_amp2(f0g[i], f2g[i], m, m_s, gr_c_.x[j]) * r2;
          cdf[i * nc + j] = acc;
        }
      }
      const std::size_t s = slot(m, i_ms);
      // P(m_S | M) = sum_L N_L |C_L(M, m_S)|^2 -- the CLOSED form, so it does
      // not inherit the grid's quadrature error.
      double p = 0.0;
      if (m - m_s == 0) {
        const double c0 = clebsch_gordan(0.0, 0.0, 1.0, m_s, 1.0, m);
        p += ad_m_.n0 / ad_m_.norm() * c0 * c0;
      }
      if (std::abs(m - m_s) <= 2) {
        const double c2 = clebsch_gordan(2.0, m - m_s, 1.0, m_s, 1.0, m);
        p += ad_m_.n2 / ad_m_.norm() * c2 * c2;
      }
      ad_weight_[s] = p;
      if (acc > 0.0) {
        for (double& v : cdf) v /= acc;
        ad_cdf_[s] = std::move(cdf);
      }
    }
  }

  // p-n: one (r, c) table per m_S
  np_cdf_.assign(3, {});
  {
    std::vector<double> g0(nr), g2(nr);
    for (std::size_t i = 0; i < nr; ++i) {
      g0[i] = table_at(np_x_, np_f0_, gr_np_.x[i]);
      g2[i] = table_at(np_x_, np_f2_, gr_np_.x[i]);
    }
    for (int i_ms = 0; i_ms < 3; ++i_ms) {
      const int m_s = kMs[i_ms];
      std::vector<double> cdf(nr * nc);
      double acc = 0.0;
      for (std::size_t i = 0; i < nr; ++i) {
        const double r2 = gr_np_.x[i] * gr_np_.x[i];
        for (std::size_t j = 0; j < nc; ++j) {
          acc += cluster_density(g0[i], g2[i], m_s, gr_c_.x[j]) * r2;
          cdf[i * nc + j] = acc;
        }
      }
      if (!(acc > 0.0)) throw std::runtime_error("empty p-n CDF");
      for (double& v : cdf) v /= acc;
      np_cdf_[i_ms] = std::move(cdf);
    }
  }

  // ------------------------- the closed-form recentred core (design 4.1(4))
  {
    const double q_max = 25.0;
    const std::size_t nq = 2501;
    std::vector<double> q = linspace(0.0, q_max, nq), ft(nq), fts(nq);
    auto ftilde = [&](double qq) {
      std::vector<double> w(core_x_.size());
      for (std::size_t i = 0; i < core_x_.size(); ++i) {
        w[i] = core_rho_[i] * core_x_[i] * core_x_[i] * sph_j0(qq * core_x_[i]);
      }
      return 4.0 * kPi * trapezoid(w, core_x_);
    };
    for (std::size_t i = 0; i < nq; ++i) ft[i] = ftilde(q[i]);
    for (std::size_t i = 0; i < nq; ++i) {
      const double a = ftilde(0.75 * q[i]);
      const double b = ftilde(0.25 * q[i]);
      fts[i] = a * b * b * b;
    }
    const std::size_t ns = 401;
    core_s_ = linspace(0.0, std::fmin(opt_.r_max_fm, core_x_.back()), ns);
    core_rho_s_.assign(ns, 0.0);
    std::vector<double> w(nq);
    for (std::size_t k = 0; k < ns; ++k) {
      for (std::size_t i = 0; i < nq; ++i) {
        w[i] = q[i] * q[i] * sph_j0(q[i] * core_s_[k]) * fts[i];
      }
      core_rho_s_[k] = trapezoid(w, q) / (2.0 * kPi * kPi);
    }
  }

  prov << "options: alpha_source="
       << (opt_.alpha_source == AlphaCoreSource::VmcHe4Density ? "VmcHe4Density"
                                                              : "Gaussian")
       << " alpha_d_source="
       << (opt_.alpha_d_source == AlphaDSource::FitRescaled   ? "FitRescaled"
           : opt_.alpha_d_source == AlphaDSource::FitRaw      ? "FitRaw"
                                                              : "OverlapRaw")
       << " theta_s=" << opt_.theta_s << " phi_s=" << opt_.phi_s
       << " alpha_cm_inflate=" << (opt_.alpha_cm_inflate ? 1 : 0)
       << " min_nn_separation_fm=" << opt_.min_nn_separation_fm
       << " alpha_d_scale=" << opt_.alpha_d_scale
       << " quadrupole_target_fm2=" << opt_.quadrupole_target_fm2
       << " exact_coherence=0"
       << " n_r=" << opt_.n_r << " n_c=" << opt_.n_c
       << " r_max_fm=" << opt_.r_max_fm << " rnp_max_fm=" << opt_.rnp_max_fm
       << "\n";
  provenance_ = prov.str();
}

// ------------------------------------------------------ analytic layer (G5)

double ClusterConfigSampler::p_d_alpha_d() const { return ad_m_.p_d(); }

double ClusterConfigSampler::tensor_dilution() const {
  return 1.0 - 0.9 * p_d_alpha_d();
}

double ClusterConfigSampler::q_int_fm2() const {
  const double lam2 = opt_.alpha_d_scale * opt_.alpha_d_scale;
  return lam2 * ad_m_.q_int / ad_m_.norm();
}

double ClusterConfigSampler::q_dd_fm2() const {
  const double lam2 = opt_.alpha_d_scale * opt_.alpha_d_scale;
  return lam2 * ad_m_.q_dd / ad_m_.norm();
}

double ClusterConfigSampler::r2_alpha_d_fm2() const {
  const double lam2 = opt_.alpha_d_scale * opt_.alpha_d_scale;
  return lam2 * ad_m_.r2 / ad_m_.norm();
}

double ClusterConfigSampler::q_matter_analytic_fm2(int m) const {
  const double f = 3.0 * static_cast<double>(m) * static_cast<double>(m) - 2.0;
  return f * ((4.0 / 3.0) * (q_int_fm2() + q_dd_fm2())
              + 2.0 * q_d_fm2_ * tensor_dilution());
}

double ClusterConfigSampler::r2_analytic_fm2() const {
  return (2.0 / 3.0) * s2_alpha_ + (1.0 / 12.0) * r2_np_
         + (2.0 / 9.0) * r2_alpha_d_fm2();
}

double ClusterConfigSampler::delta_perp_analytic_fm2(int m) const {
  const double st = std::sin(opt_.theta_s);
  return q_matter_analytic_fm2(m) / 12.0 * st * st * std::cos(2.0 * opt_.phi_s);
}

double ClusterConfigSampler::a2_from_geometry(double t_abs, int m) const {
  return a2_from_quadrupole(q_matter_analytic_fm2(1), 6, t_abs, m);
}

double ClusterConfigSampler::eps_b0_equivalent() const {
  const double delta_gev2 =
      delta_perp_analytic_fm2(1) / (HBARC_GEV_FM * HBARC_GEV_FM);
  return delta_gev2 / gaussian_slope(std::sqrt(LI6_R2_POINT_FM2));
}

double ClusterConfigSampler::quadrupole_for_eta(double eta_target) const {
  const double eta_now = asymptotic_ds_ratio();
  if (!(std::fabs(eta_now) > 0.0)) {
    throw std::runtime_error("quadrupole_for_eta: this source's eta is zero, "
                             "so no dial setting reaches a target");
  }
  // eta(s) = s eta(1) EXACTLY, so the dial that lands on eta_target is the
  // current dial scaled by the ratio of the etas.
  const double s = dial_s_ * (eta_target / eta_now);
  if (!(s >= 0.0 && s <= 1.0)) {
    std::ostringstream e;
    e << "quadrupole_for_eta: eta = " << eta_target << " needs the (G9) dial "
      << "s = " << s << ", outside [0, 1]; this source reaches eta in [0, "
      << eta_now / dial_s_ << "]";
    throw std::runtime_error(e.str());
  }
  // The dial's OWN q_charge(s) (cluster_config.cpp, the (G9) bisection), on
  // the pre-dial moments.  Not a second copy of the map: the same five lines
  // are what the bisection inverts.
  const double lam2 = opt_.alpha_d_scale * opt_.alpha_d_scale;
  const double qi = lam2 * ad_m_base_.q_int, qd = lam2 * ad_m_base_.q_dd;
  const double pd = ad_m_base_.p_d();
  const double n = 1.0 - (1.0 - s * s) * pd;
  const double q = (s * qi + s * s * qd) / n;
  const double pdp = s * s * pd / n;
  return 0.5 * ((4.0 / 3.0) * q + 2.0 * q_d_fm2_ * (1.0 - 0.9 * pdp));
}

double ClusterConfigSampler::match_li6_radius(double target_rms_fm) const {
  const double rest = (2.0 / 3.0) * s2_alpha_ + (1.0 / 12.0) * r2_np_;
  const double want = target_rms_fm * target_rms_fm - rest;
  if (!(want > 0.0)) {
    throw std::runtime_error("match_li6_radius: the alpha and deuteron terms "
                             "alone already exceed the target <r^2>");
  }
  const double base = ad_m_.r2 / ad_m_.norm();   // unscaled <R^2>
  return std::sqrt(want * 9.0 / (2.0 * base));
}

std::array<double, 3> ClusterConfigSampler::quadrupole_band_fm2() const {
  return {{LI6_QUADRUPOLE_FM2, LI6_QUADRUPOLE_GFMC_FM2,
           0.5 * q_matter_analytic_fm2(1)}};
}

double ClusterConfigSampler::asymptotic_ds_ratio() const {
  // kappa and the Sommerfeld parameter come from the channel that already
  // owns them; nothing is retyped.
  const ClusterChannel& ch = LI6_ALPHA_TAG();
  const double kappa_fm = ch.kappa() / HBARC_GEV_FM;
  const double mu = ch.m_spec() * ch.m_partner() / (ch.m_spec() + ch.m_partner());
  const double eta_c = 2.0 * ALPHA_EM * mu / ch.kappa();
  double acc = 0.0;
  int n = 0;
  for (double r = kDsRatioRLoFm; r <= kDsRatioRHiFm + 1e-9; r += 1.0) {
    const double f0 = table_at(ad_x_, ad_f0_, r);
    const double f2 = table_at(ad_x_, ad_f2_, r);
    const double z = 2.0 * kappa_fm * r;
    const double wr = whittaker_w(2, eta_c, z) / whittaker_w(0, eta_c, z);
    acc += (f2 / f0) / wr;
    ++n;
  }
  return acc / static_cast<double>(n);
}

// ----------------------------------------------------------- the densities

double ClusterConfigSampler::rho_alpha(double r_fm) const {
  return table_at(core_x_, core_rho_, r_fm);
}

double ClusterConfigSampler::rho_alpha_recentred(double s_fm) const {
  return table_at(core_s_, core_rho_s_, s_fm);
}

double ClusterConfigSampler::rho_alpha_d(double R_fm, double c, int m,
                                         double m_s) const {
  const double x = R_fm / opt_.alpha_d_scale;
  const double j = opt_.alpha_d_scale * opt_.alpha_d_scale
                   * opt_.alpha_d_scale;
  return cluster_amp2(table_at(ad_x_, ad_f0_, x), table_at(ad_x_, ad_f2_, x), m,
                      static_cast<int>(std::lround(m_s)), c) / j;
}

double ClusterConfigSampler::rho_alpha_d_summed(double R_fm, double c,
                                                int m) const {
  double s = 0.0;
  for (int i = 0; i < 3; ++i) s += rho_alpha_d(R_fm, c, m, kMs[i]);
  return s;
}

double ClusterConfigSampler::rho_np(double r_fm, double c, double m_s) const {
  return cluster_density(table_at(np_x_, np_f0_, r_fm),
                         table_at(np_x_, np_f2_, r_fm),
                         static_cast<int>(std::lround(m_s)), c)
         / np_m_.norm();
}

const std::vector<double>& ClusterConfigSampler::alpha_d_wave(int l) const {
  if (l == 0) return ad_f0_;
  if (l == 2) return ad_f2_;
  throw std::runtime_error("alpha_d_wave: L is 0 or 2");
}

const std::vector<double>& ClusterConfigSampler::np_wave(int l) const {
  if (l == 0) return np_f0_;
  if (l == 2) return np_f2_;
  throw std::runtime_error("np_wave: L is 0 or 2");
}

double ClusterConfigSampler::p_ms(int m, double m_s) const {
  return ad_weight_[slot(m, static_cast<int>(ms_index(m_s)))];
}

std::size_t ClusterConfigSampler::ms_index(double m_s) const {
  for (int i = 0; i < 3; ++i) {
    if (std::fabs(static_cast<double>(kMs[i]) - m_s) < 1e-9) {
      return static_cast<std::size_t>(i);
    }
  }
  throw std::runtime_error("m_S is not a channel-spin projection");
}

const std::vector<double>& ClusterConfigSampler::ad_cdf(int m,
                                                        double m_s) const {
  const std::vector<double>& c = ad_cdf_[slot(m, static_cast<int>(ms_index(m_s)))];
  if (c.empty()) throw std::runtime_error("this (M, m_S) branch is forbidden");
  return c;
}

const std::vector<double>& ClusterConfigSampler::np_cdf(double m_s) const {
  return np_cdf_[ms_index(m_s)];
}

// -------------------------------------------------------- grid quadratures

double ClusterConfigSampler::r2_grid_fm2() const {
  // The alpha-d and p-n RADIAL marginals are m_S-independent (summing (G1)
  // over m_S and integrating over the angles leaves sum_L f_L^2), so these
  // are one-dimensional cell sums with no angular grid in them.
  double wr = 0.0, w = 0.0;
  for (std::size_t i = 0; i < gr_ad_.x.size(); ++i) {
    const double x = gr_ad_.x[i];
    const double f0 = table_at(ad_x_, ad_f0_, x), f2 = table_at(ad_x_, ad_f2_, x);
    const double p = (f0 * f0 + f2 * f2) * x * x;
    w += p;
    wr += p * x * x;
  }
  const double lam2 = opt_.alpha_d_scale * opt_.alpha_d_scale;
  const double r2_ad = lam2 * wr / w;
  double vr = 0.0, v = 0.0;
  for (std::size_t i = 0; i < gr_np_.x.size(); ++i) {
    const double x = gr_np_.x[i];
    const double f0 = table_at(np_x_, np_f0_, x), f2 = table_at(np_x_, np_f2_, x);
    const double p = (f0 * f0 + f2 * f2) * x * x;
    v += p;
    vr += p * x * x;
  }
  const double r2_np = vr / v;
  double ar = 0.0, a = 0.0;
  for (std::size_t i = 0; i < gr_core_.x.size(); ++i) {
    const double x = gr_core_.x[i];
    const double p = table_at(core_x_, core_rho_, x) * x * x;
    a += p;
    ar += p * x * x;
  }
  // 0.75 * <v^2> is the cell sum's own version of the INDEPENDENT-draw closed
  // form; with a hard core the sampler no longer draws from this cell CDF
  // alone, so the constructor's measured value is the honest one here too.
  const double s2 = opt_.min_nn_separation_fm > 0.0 ? s2_alpha_
                                                    : 0.75 * ar / a;
  return (2.0 / 3.0) * s2 + (1.0 / 12.0) * r2_np + (2.0 / 9.0) * r2_ad;
}

double ClusterConfigSampler::q_matter_grid_fm2(int m) const {
  const std::size_t nr = opt_.n_r, nc = opt_.n_c;
  // alpha-d: <R^2 (3 c^2 - 1)> over the m_S-summed joint density of state m.
  double num = 0.0, den = 0.0;
  for (std::size_t i = 0; i < nr; ++i) {
    const double x = gr_ad_.x[i], x2 = x * x;
    const double f0 = table_at(ad_x_, ad_f0_, x), f2 = table_at(ad_x_, ad_f2_, x);
    for (std::size_t j = 0; j < nc; ++j) {
      const double c = gr_c_.x[j];
      double a2 = 0.0;
      for (int k = 0; k < 3; ++k) a2 += cluster_amp2(f0, f2, m, kMs[k], c);
      const double wgt = a2 * x2;
      den += wgt;
      num += wgt * x2 * (3.0 * c * c - 1.0);
    }
  }
  const double lam2 = opt_.alpha_d_scale * opt_.alpha_d_scale;
  double q = (4.0 / 3.0) * lam2 * num / den;
  // deuteron: (1/2) sum_{m_S} P(m_S|M) <r^2 (3 c^2 - 1)>_{m_S}
  for (int k = 0; k < 3; ++k) {
    const int m_s = kMs[k];
    double n2 = 0.0, d2 = 0.0;
    for (std::size_t i = 0; i < nr; ++i) {
      const double x = gr_np_.x[i], x2 = x * x;
      const double f0 = table_at(np_x_, np_f0_, x);
      const double f2 = table_at(np_x_, np_f2_, x);
      for (std::size_t j = 0; j < nc; ++j) {
        const double c = gr_c_.x[j];
        const double wgt = cluster_density(f0, f2, m_s, c) * x2;
        d2 += wgt;
        n2 += wgt * x2 * (3.0 * c * c - 1.0);
      }
    }
    q += 0.5 * ad_weight_[slot(m, k)] * n2 / d2;
  }
  return q;
}

double ClusterConfigSampler::p2_alpha_d_grid(int m, double m_s) const {
  const std::size_t nc = opt_.n_c;
  const int ms = static_cast<int>(std::lround(m_s));
  double num = 0.0, den = 0.0;
  for (std::size_t i = 0; i < gr_ad_.x.size(); ++i) {
    const double x = gr_ad_.x[i], x2 = x * x;
    const double f0 = table_at(ad_x_, ad_f0_, x), f2 = table_at(ad_x_, ad_f2_, x);
    for (std::size_t j = 0; j < nc; ++j) {
      const double c = gr_c_.x[j];
      const double wgt = cluster_amp2(f0, f2, m, ms, c) * x2;
      den += wgt;
      num += wgt * 0.5 * (3.0 * c * c - 1.0);
    }
  }
  if (!(den > 0.0)) throw std::runtime_error("this (M, m_S) branch is forbidden");
  return num / den;
}

// ------------------------------------------------------------- the sampler

void ClusterConfigSampler::draw_alpha(
    Rng& rng, std::array<std::array<double, 3>, 4>& s) const {
  const std::size_t nr = opt_.n_r;
  const double dmin2 = opt_.min_nn_separation_fm * opt_.min_nn_separation_fm;
  for (int i = 0; i < 4; ++i) {
    int tries = 0;
    while (true) {
      const std::size_t cell = static_cast<std::size_t>(
          std::lower_bound(core_cdf_.begin(), core_cdf_.end(), rng.uniform())
          - core_cdf_.begin());
      const std::size_t ir = cell < nr ? cell : nr - 1;
      const double v =
          std::fabs(gr_core_.x[ir] + (rng.uniform() - 0.5) * gr_core_.dx);
      const double c = 2.0 * rng.uniform() - 1.0;
      const double phi = 2.0 * kPi * rng.uniform();
      const double st = std::sqrt(std::fmax(1.0 - c * c, 0.0));
      s[i] = {{v * st * std::cos(phi), v * st * std::sin(phi), v * c}};
      if (dmin2 <= 0.0) break;
      bool ok = true;
      for (int j = 0; j < i; ++j) {
        const double dx = s[i][0] - s[j][0], dy = s[i][1] - s[j][1];
        const double dz = s[i][2] - s[j][2];
        if (dx * dx + dy * dy + dz * dz < dmin2) {
          ok = false;
          break;
        }
      }
      if (ok) break;
      if (++tries > 1000) {
        throw std::runtime_error(
            "min_nn_separation_fm could not be satisfied in 1000 redraws; "
            "the requested separation is too large for this core density");
      }
    }
  }
  // The CM constraint, applied to the core only: sum_i s_i = 0 EXACTLY.
  for (int k = 0; k < 3; ++k) {
    const double mean = 0.25 * (s[0][k] + s[1][k] + s[2][k] + s[3][k]);
    for (int i = 0; i < 4; ++i) s[i][k] -= mean;
  }
}

ClusterConfig ClusterConfigSampler::sample(Rng& rng, int m_ion) const {
  if (m_ion < -1 || m_ion > 1) {
    throw std::runtime_error("m_ion must be +1, 0 or -1");
  }
  const std::size_t nc = opt_.n_c;

  // 1. m_S from P(m_S | M), in the +1, 0, -1 order.
  double tot = 0.0;
  for (int i = 0; i < 3; ++i) tot += ad_weight_[slot(m_ion, i)];
  double u = rng.uniform() * tot, acc = 0.0;
  int i_ms = 2;
  for (int i = 0; i < 3; ++i) {
    acc += ad_weight_[slot(m_ion, i)];
    if (u <= acc) {
      i_ms = i;
      break;
    }
  }
  const int m_s = kMs[i_ms];

  // 2. (R, cos theta_R) from the m_S-CONDITIONED table, uniform phi_R.
  const std::vector<double>& cdf = ad_cdf_[slot(m_ion, i_ms)];
  if (cdf.empty()) throw std::runtime_error("empty alpha-d branch drawn");
  std::size_t cell = static_cast<std::size_t>(
      std::lower_bound(cdf.begin(), cdf.end(), rng.uniform()) - cdf.begin());
  if (cell >= cdf.size()) cell = cdf.size() - 1;
  const double big_r =
      std::fabs(gr_ad_.x[cell / nc] + (rng.uniform() - 0.5) * gr_ad_.dx)
      * opt_.alpha_d_scale;
  const double c_r =
      clip(gr_c_.x[cell % nc] + (rng.uniform() - 0.5) * gr_c_.dx, -1.0, 1.0);
  const double phi_r = 2.0 * kPi * rng.uniform();

  // 3. (r, cos theta_r) of the p-n pair, for THAT m_S.
  const std::vector<double>& ncdf = np_cdf_[i_ms];
  std::size_t ncell = static_cast<std::size_t>(
      std::lower_bound(ncdf.begin(), ncdf.end(), rng.uniform()) - ncdf.begin());
  if (ncell >= ncdf.size()) ncell = ncdf.size() - 1;
  const double r_np =
      std::fabs(gr_np_.x[ncell / nc] + (rng.uniform() - 0.5) * gr_np_.dx);
  const double c_n =
      clip(gr_c_.x[ncell % nc] + (rng.uniform() - 0.5) * gr_c_.dx, -1.0, 1.0);
  const double phi_n = 2.0 * kPi * rng.uniform();

  // 4. the four core nucleons, recentred.
  std::array<std::array<double, 3>, 4> s{};
  draw_alpha(rng, s);

  // 5. assemble, spin frame (axis = +z).
  const double sr = std::sqrt(std::fmax(1.0 - c_r * c_r, 0.0));
  const std::array<double, 3> big{{big_r * sr * std::cos(phi_r),
                                   big_r * sr * std::sin(phi_r), big_r * c_r}};
  const double sn = std::sqrt(std::fmax(1.0 - c_n * c_n, 0.0));
  const std::array<double, 3> rel{{r_np * sn * std::cos(phi_n),
                                   r_np * sn * std::sin(phi_n), r_np * c_n}};
  ClusterConfig out;
  out.m_ion = m_ion;
  out.m_s = static_cast<double>(m_s);
  out.r_ad = big_r;
  out.r_np = r_np;
  std::array<std::array<double, 3>, 6> pos{};
  for (int i = 0; i < 4; ++i) {
    for (int k = 0; k < 3; ++k) pos[i][k] = -big[k] / 3.0 + s[i][k];
  }
  for (int k = 0; k < 3; ++k) {
    pos[4][k] = 2.0 * big[k] / 3.0 + 0.5 * rel[k];
    pos[5][k] = 2.0 * big[k] / 3.0 - 0.5 * rel[k];
  }

  // 6. rotate the spin frame into the lab: R_z(phi_s) R_y(theta_s), the same
  //    composition `boost_spectator` applies to the spectator momentum.
  const double ct = std::cos(opt_.theta_s), stt = std::sin(opt_.theta_s);
  const double cp = std::cos(opt_.phi_s), sp = std::sin(opt_.phi_s);
  for (int i = 0; i < 6; ++i) {
    double x = pos[i][0], y = pos[i][1], z = pos[i][2];
    const double x1 = ct * x + stt * z;
    const double z1 = -stt * x + ct * z;
    x = cp * x1 - sp * y;
    y = sp * x1 + cp * y;
    out.nucleon[i].x = x;
    out.nucleon[i].y = y;
    out.nucleon[i].z = z1;
    out.nucleon[i].isospin = (i < 4) ? ((i % 2 == 0) ? +1 : -1)
                                     : ((i == 4) ? +1 : -1);
    out.nucleon[i].cluster = (i < 4) ? 0 : 1;
  }
  return out;
}

namespace {

/// Accumulate one configuration's three moments about the axis n(theta, phi).
void config_moments(const ClusterConfig& cfg, const std::array<double, 3>& n,
                    double& r2, double& q, double& dperp,
                    std::array<double, 3>& cm) {
  r2 = 0.0;
  q = 0.0;
  dperp = 0.0;
  for (int i = 0; i < 6; ++i) {
    const double x = cfg.nucleon[i].x, y = cfg.nucleon[i].y;
    const double z = cfg.nucleon[i].z;
    const double rr = x * x + y * y + z * z;
    const double za = x * n[0] + y * n[1] + z * n[2];
    r2 += rr;
    q += 3.0 * za * za - rr;
    dperp += x * x - y * y;
    cm[0] += x;
    cm[1] += y;
    cm[2] += z;
  }
  r2 /= 6.0;
  dperp /= 6.0;
}

}  // namespace

ClusterConfigSet ClusterConfigSampler::sample_set(std::size_t n, int m_ion,
                                                  std::uint64_t seed,
                                                  std::uint64_t run) const {
  ClusterConfigSet set;
  set.m_ion = m_ion;
  set.theta_s = opt_.theta_s;
  set.phi_s = opt_.phi_s;
  set.seed = seed;
  set.run = run;
  set.provenance = provenance_;
  set.config.resize(n);
  const std::array<double, 3> axis{{std::sin(opt_.theta_s) * std::cos(opt_.phi_s),
                                    std::sin(opt_.theta_s) * std::sin(opt_.phi_s),
                                    std::cos(opt_.theta_s)}};
  double s1 = 0.0, s2 = 0.0, q1 = 0.0, q2 = 0.0, d1 = 0.0, d2 = 0.0;
  double cm_max = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    Rng rng(seed, run, kConfigStream, static_cast<std::uint64_t>(i));
    set.config[i] = sample(rng, m_ion);
    std::array<double, 3> cm{{0.0, 0.0, 0.0}};
    double r2 = 0.0, q = 0.0, dp = 0.0;
    config_moments(set.config[i], axis, r2, q, dp, cm);
    s1 += r2;
    s2 += r2 * r2;
    q1 += q;
    q2 += q * q;
    d1 += dp;
    d2 += dp * dp;
    for (int k = 0; k < 3; ++k) cm_max = std::fmax(cm_max, std::fabs(cm[k]));
  }
  if (n > 0) {
    const double dn = static_cast<double>(n);
    set.r2_mean_fm2 = s1 / dn;
    set.q_matter_fm2 = q1 / dn;
    set.delta_perp_fm2 = d1 / dn;
    if (n > 1) {
      const double den = dn - 1.0;
      set.r2_sd_fm2 = std::sqrt(std::fmax(s2 - s1 * s1 / dn, 0.0) / den);
      set.q_matter_sd_fm2 = std::sqrt(std::fmax(q2 - q1 * q1 / dn, 0.0) / den);
      set.delta_perp_sd_fm2 = std::sqrt(std::fmax(d2 - d1 * d1 / dn, 0.0) / den);
    }
  }
  set.cm = {{cm_max, cm_max, cm_max}};
  return set;
}

ClusterConfigSet ClusterConfigSampler::sample_set_unpolarized(
    std::size_t n, std::uint64_t seed, std::uint64_t run) const {
  ClusterConfigSet set;
  set.m_ion = -2;
  set.theta_s = opt_.theta_s;
  set.phi_s = opt_.phi_s;
  set.seed = seed;
  set.run = run;
  set.provenance = provenance_;
  set.config.resize(n);
  const std::array<double, 3> axis{{std::sin(opt_.theta_s) * std::cos(opt_.phi_s),
                                    std::sin(opt_.theta_s) * std::sin(opt_.phi_s),
                                    std::cos(opt_.theta_s)}};
  static const int kOrder[3] = {1, 0, -1};
  double s1 = 0.0, s2 = 0.0, q1 = 0.0, q2 = 0.0, d1 = 0.0, d2 = 0.0;
  double cm_max = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    Rng rng(seed, run, kConfigStream, static_cast<std::uint64_t>(i));
    set.config[i] = sample(rng, kOrder[i % 3]);
    std::array<double, 3> cm{{0.0, 0.0, 0.0}};
    double r2 = 0.0, q = 0.0, dp = 0.0;
    config_moments(set.config[i], axis, r2, q, dp, cm);
    s1 += r2;
    s2 += r2 * r2;
    q1 += q;
    q2 += q * q;
    d1 += dp;
    d2 += dp * dp;
    for (int k = 0; k < 3; ++k) cm_max = std::fmax(cm_max, std::fabs(cm[k]));
  }
  if (n > 0) {
    const double dn = static_cast<double>(n);
    set.r2_mean_fm2 = s1 / dn;
    set.q_matter_fm2 = q1 / dn;
    set.delta_perp_fm2 = d1 / dn;
    if (n > 1) {
      const double den = dn - 1.0;
      set.r2_sd_fm2 = std::sqrt(std::fmax(s2 - s1 * s1 / dn, 0.0) / den);
      set.q_matter_sd_fm2 = std::sqrt(std::fmax(q2 - q1 * q1 / dn, 0.0) / den);
      set.delta_perp_sd_fm2 = std::sqrt(std::fmax(d2 - d1 * d1 / dn, 0.0) / den);
    }
  }
  set.cm = {{cm_max, cm_max, cm_max}};
  return set;
}

// -------------------------------------------------------------- the writer

namespace {

std::string json_num(double v) {
  std::ostringstream s;
  s << std::setprecision(17) << v;
  return s.str();
}

/// For numbers that appear inside a human-readable CAVEAT string, where 17
/// significant digits are noise.  Never used for a machine-read field.
std::string num_fixed(double v, int digits) {
  std::ostringstream s;
  s << std::fixed << std::setprecision(digits) << v;
  return s.str();
}

std::size_t file_bytes(const std::string& path) {
  std::ifstream in(path, std::ios::binary | std::ios::ate);
  if (!in) return 0;
  return static_cast<std::size_t>(in.tellg());
}

/// The `ndx  s-wave  d-wave` line `li6.ad` prints under its k-space block.
std::array<double, 3> li6_ad_printed_norms(const std::string& path) {
  const std::vector<std::string> lines = split_lines(read_file(path));
  for (std::size_t i = 0; i < lines.size(); ++i) {
    if (lines[i].find("ndx") != std::string::npos
        && lines[i].find("s-wave") != std::string::npos) {
      for (std::size_t j = i + 1; j < lines.size(); ++j) {
        const std::vector<double> v = line_numbers(lines[j]);
        if (v.size() >= 3) return {{v[0], v[1], v[2]}};
      }
    }
  }
  return {{0.0, 0.0, 0.0}};
}

/// he4.density's own `4*PI*TOTINT(RHORP*R**2:R)` and rms lines.
std::array<double, 2> density_printed_norms(const std::string& path) {
  const std::vector<std::string> lines = split_lines(read_file(path));
  std::array<double, 2> out{{0.0, 0.0}};
  for (const std::string& l : lines) {
    const std::size_t eq = l.find('=');
    if (eq == std::string::npos) continue;
    if (l.find("TOTINT(RHORP*R**2") != std::string::npos) {
      out[0] = std::strtod(l.c_str() + eq + 1, nullptr);
    } else if (l.find("SQRT(4*PI*TOTINT(RHORP*R**4") != std::string::npos) {
      out[1] = std::strtod(l.c_str() + eq + 1, nullptr);
    }
  }
  return out;
}

const char* alpha_source_name(AlphaCoreSource s) {
  return s == AlphaCoreSource::VmcHe4Density ? "VmcHe4Density" : "Gaussian";
}

const char* alpha_d_source_name(AlphaDSource s) {
  switch (s) {
    case AlphaDSource::FitRescaled: return "FitRescaled";
    case AlphaDSource::FitRaw:      return "FitRaw";
    default:                        return "OverlapRaw";
  }
}

std::string sidecar_json(const ClusterConfigSet& set,
                         const ClusterConfigSampler& smp) {
  const ClusterConfigOptions& o = smp.options();
  const std::string fit = data_path(VMC_LI6_AD_FIT);
  const std::string ov = data_path(VMC_LI6_OVERLAP);
  const std::string deut = data_path(VMC_DEUTERON_WAVE);
  const std::string he4 = data_path(VMC_HE4_DENSITY);
  const std::array<double, 3> ad_norms = li6_ad_printed_norms(ov);
  const std::array<double, 2> he4_norms = density_printed_norms(he4);
  const FdeutRTable d = read_fdeut_r(deut);
  const std::array<double, 3> band = smp.quadrupole_band_fm2();
  // `li6.adr.fit` prints no norm of its own, so the entry's "recomputed_norm"
  // is RECOMPUTED -- integrated here from the file, never a typed constant
  // (docs/CONVENTIONS.md: no physics number lives twice).
  const AnlTable fit_t = read_anl_plain(fit, 2);
  const RadialMoments fit_m =
      radial_moments(fit_t.x, fit_t.col[0], fit_t.col[1]);
  // Which tables THIS option set actually reads.  The entries stay so that
  // the sidecar always lists the module's four inputs with their md5s, but
  // "used" says which of them fed the configurations in this file.
  const bool use_fit = o.alpha_d_source != AlphaDSource::OverlapRaw;
  const bool use_ov = o.alpha_d_source == AlphaDSource::OverlapRaw
                      || o.alpha_d_source == AlphaDSource::FitRescaled;
  const bool use_he4 = o.alpha_source == AlphaCoreSource::VmcHe4Density;
  // With a hard core the alpha's <s^2> is a numerical estimate, so (G6)'s
  // <r^2> -- and every moment built on it -- is approximate at that precision.
  const bool hard_core = o.min_nn_separation_fm > 0.0;

  std::ostringstream j;
  j << "{\"generator\":\"LiPolGen cluster_config\",\"version\":\""
    << LIPOLGEN_VERSION << "\",\n \"git\":null,\n"
    << " \"nucleus\":\"6Li\",\"A\":6,\"Z\":3,\"n_config\":"
    << set.config.size() << ",\"m_ion\":"
    << (set.m_ion == -2 ? std::string("\"unpolarized\"")
                        : std::to_string(set.m_ion))
    << ",\n \"axis\":{\"theta_s\":" << json_num(set.theta_s)
    << ",\"phi_s\":" << json_num(set.phi_s)
    << ",\"frame\":\"ion rest frame\"},\n"
    << " \"units\":\"fm\",\"column_order\":[\"x\",\"y\",\"z\"],"
    << "\"cm_convention\":\"sum r_i = 0\",\n"
    << " \"rng\":{\"seed\":" << set.seed << ",\"run\":" << set.run
    << ",\"stream\":\"CONFIG\",\"stream_value\":" << kConfigStream << "},\n"
    << " \"inputs\":[";
  j << "{\"file\":\"" << VMC_LI6_AD_FIT << "\",\"bytes\":" << file_bytes(fit)
    << ",\"used\":" << (use_fit ? "true" : "false")
    << ",\"printed_norm\":null,\"recomputed_norm\":{\"integral\":"
    << json_num(fit_m.norm()) << ",\"p_d\":" << json_num(fit_m.p_d())
    << "},\"md5\":null},\n           ";
  j << "{\"file\":\"" << VMC_LI6_OVERLAP << "\",\"bytes\":" << file_bytes(ov)
    << ",\"used\":" << (use_ov ? "true" : "false")
    << ",\"printed_norm\":{\"ndx\":" << json_num(ad_norms[0]) << ",\"s\":"
    << json_num(ad_norms[1]) << ",\"d\":" << json_num(ad_norms[2])
    << "},\"md5\":null},\n           ";
  j << "{\"file\":\"" << VMC_DEUTERON_WAVE << "\",\"bytes\":"
    << file_bytes(deut) << ",\"used\":true,\"printed_norm\":{\"dstate\":"
    << json_num(d.dstate) << ",\"qm\":" << json_num(d.qm) << ",\"rd\":"
    << json_num(d.rd) << "},\"md5\":null},\n           ";
  j << "{\"file\":\"" << VMC_HE4_DENSITY << "\",\"bytes\":" << file_bytes(he4)
    << ",\"used\":" << (use_he4 ? "true" : "false")
    << ",\"printed_norm\":{\"norm_4pi\":" << json_num(he4_norms[0])
    << ",\"rms\":" << json_num(he4_norms[1]) << "},\"md5\":null}],\n";
  j << " \"options\":{\"alpha_source\":\"" << alpha_source_name(o.alpha_source)
    << "\",\"alpha_d_source\":\"" << alpha_d_source_name(o.alpha_d_source)
    << "\",\n            \"theta_s\":" << json_num(o.theta_s)
    << ",\"phi_s\":" << json_num(o.phi_s)
    << ",\n            \"alpha_cm_inflate\":"
    << (o.alpha_cm_inflate ? "true" : "false")
    << ",\"min_nn_separation_fm\":" << json_num(o.min_nn_separation_fm)
    << ",\n            \"alpha_d_scale\":" << json_num(o.alpha_d_scale)
    << ",\"quadrupole_target_fm2\":" << json_num(o.quadrupole_target_fm2)
    << ",\n            \"quadrupole_dial_s\":"
    << json_num(smp.quadrupole_dial_s())
    << ",\n            \"exact_coherence\":"
    << (o.exact_coherence ? "true" : "false")
    << ",\n            \"n_r\":" << o.n_r << ",\"n_c\":" << o.n_c
    << ",\"r_max_fm\":" << json_num(o.r_max_fm)
    << ",\"rnp_max_fm\":" << json_num(o.rnp_max_fm) << "},\n";
  j << " \"moments\":{\"r2_mean_fm2\":" << json_num(smp.r2_analytic_fm2())
    << ",\"r_rms_fm\":" << json_num(std::sqrt(smp.r2_analytic_fm2()))
    << ",\n             \"q_matter_fm2\":"
    << json_num(smp.q_matter_analytic_fm2(1))
    << ",\"delta_perp_fm2\":" << json_num(smp.delta_perp_analytic_fm2(1))
    << ",\n             \"eps_b0_equivalent\":"
    << json_num(smp.eps_b0_equivalent())
    << ",\"p_d_alpha_d\":" << json_num(smp.p_d_alpha_d())
    << ",\n             \"tensor_dilution\":" << json_num(smp.tensor_dilution())
    << ",\"asymptotic_ds_ratio\":" << json_num(smp.asymptotic_ds_ratio())
    << ",\n             \"alpha_s2_fm2\":" << json_num(smp.r2_alpha_fm2())
    << ",\"alpha_s2_method\":\""
    << (hard_core ? "numerical MC (hard core), see provenance"
                  : "closed form (3/4)<v^2>, exact")
    << "\",\n             \"r2_mean_is_approximate\":"
    << (hard_core ? "true" : "false")
    << ",\n             \"sampled\":{\"r2_mean_fm2\":"
    << json_num(set.r2_mean_fm2) << ",\"r2_sd_fm2\":"
    << json_num(set.r2_sd_fm2)
    << ",\n                        \"q_matter_fm2\":"
    << json_num(set.q_matter_fm2) << ",\"q_matter_sd_fm2\":"
    << json_num(set.q_matter_sd_fm2)
    << ",\n                        \"delta_perp_fm2\":"
    << json_num(set.delta_perp_fm2) << ",\"delta_perp_sd_fm2\":"
    << json_num(set.delta_perp_sd_fm2) << "}},\n";
  j << " \"quadrupole_band_fm2\":{\"measured\":" << json_num(band[0])
    << ",\"gfmc_av18_il7\":" << json_num(band[1])
    << ",\"gfmc_av18_il7_err\":" << json_num(LI6_QUADRUPOLE_GFMC_ERR_FM2)
    << ",\"alpha_d_model\":" << json_num(band[2]) << "},\n";
  // The overshoot factor is COMPUTED from the band this run actually has --
  // never a literal, which would contradict the band above whenever the dial
  // or the alpha-d scale is used (design sec. 2.7: not three numbers for one
  // quantity).  The three sources' own Q_charge live in design sec. 8.
  const double overshoot = std::fabs(band[0]) > 0.0
                               ? std::fabs(band[2] / band[0]) : 0.0;
  const bool tuned = smp.quadrupole_dial_s() != 1.0 || o.alpha_d_scale != 1.0;
  j << " \"caveats\":[\"Q_charge(model) = " << num_fixed(band[2], 4)
    << " fm^2 vs measured " << num_fixed(band[0], 4) << " (TUNL) and GFMC "
       "AV18+IL7\",\n"
    << "            \"" << num_fixed(band[1], 2) << " +- "
    << num_fixed(LI6_QUADRUPOLE_GFMC_ERR_FM2, 2) << " (Pastore et al., PRC "
       "87, 035503 (2013), arXiv:1212.3375):\",\n"
    << "            \"|model/measured| = " << num_fixed(overshoot, 1)
    << ". Do NOT quote a tensor number from this sampler without\",\n"
    << "            \"the band -- see docs/open_items/run_2026-09-02/"
       "design_G_cluster_config.md sec. 2.7\",\n"
    << "            \"(sec. 8 for the three alpha-d sources' own Q_charge) "
       "and docs/OPEN_ITEMS_SOLUTIONS.md sec. 11.\"";
  if (tuned) {
    j << ",\n            \"DIALLED/SCALED: quadrupole_dial_s = "
      << num_fixed(smp.quadrupole_dial_s(), 4) << ", alpha_d_scale = "
      << num_fixed(o.alpha_d_scale, 4) << " -- the factor above is this\",\n"
      << "            \"tuned geometry's, NOT the natural wave functions' "
         "overshoot.\"";
  }
  if (hard_core) {
    j << ",\n            \"min_nn_separation_fm = "
      << num_fixed(o.min_nn_separation_fm, 2) << " fm: the alpha's <s^2> is a "
         "NUMERICAL estimate (the closed\",\n"
      << "            \"form (3/4)<v^2> holds only for independent draws), so "
         "moments.r2_mean_fm2 / r_rms_fm\",\n"
      << "            \"are approximate and rho_alpha_recentred() is still "
         "the independent-draw density.\"";
  }
  j << "]}\n";
  return j.str();
}

}  // namespace

std::size_t write_snd_configs(const ClusterConfigSet& set,
                              const ClusterConfigSampler& sampler,
                              const std::string& path, SndConfigFormat fmt) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write configurations to " + path);
  out << std::setprecision(17);
  const std::string meta = sidecar_json(set, sampler);
  if (fmt == SndConfigFormat::Annotated) {
    std::istringstream in(meta);
    std::string line;
    while (std::getline(in, line)) out << "# " << line << "\n";
  }
  for (const ClusterConfig& c : set.config) {
    for (int i = 0; i < 6; ++i) {
      out << c.nucleon[i].x << " " << c.nucleon[i].y << " " << c.nucleon[i].z
          << " ";
    }
    for (int i = 0; i < 6; ++i) out << c.nucleon[i].isospin << " ";
    out << c.m_ion << "\n";
  }
  out.close();
  std::ofstream side(path + ".meta.json");
  if (!side) {
    throw std::runtime_error("cannot write sidecar " + path + ".meta.json");
  }
  side << meta;
  return set.config.size();
}

}  // namespace lipolgen
