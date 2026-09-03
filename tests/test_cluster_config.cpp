// Phase G(i): the alpha+d configuration sampler for polarized 6Li.
// One TEST_CASE per identity of
// docs/open_items/run_2026-09-02/design_G_cluster_config.md sec. 7.
//
// HOW THE MONTE CARLO TOLERANCES ARE SET.  Every MC gate below is written as
// 5 * set.<moment>_sd / sqrt(n) -- the PER-CONFIGURATION spread the sampler
// itself reports -- never as a hand-typed absolute number.  A single 6Li
// configuration has sum_i (3 z_i^2 - r_i^2) scattered over ~29 fm^2, so an
// absolute gate is meaningless without the variance beside it.  Tolerances
// that must be tight and absolute live on the DETERMINISTIC layer instead.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "check_close.hpp"
#include "doctest.h"
#include "json_min.hpp"
#include "lipolgen/cluster_config.hpp"
#include "lipolgen/coherent.hpp"
#include "lipolgen/constants.hpp"
#include "lipolgen/fsi.hpp"
#include "lipolgen/numerics.hpp"
#include "lipolgen/rc.hpp"
#include "lipolgen/tagged.hpp"

using namespace lipolgen;

namespace {

constexpr double kQdFm2 = 0.26967026805838806;   // (G4)/4 on fdeut.av18

/// The four core nucleons of one configuration, relative to their own c.m.
/// (which is what `s_i` is); |s_i| is rotation invariant.
std::vector<double> core_radii(const ClusterConfigSet& set) {
  std::vector<double> out;
  out.reserve(set.config.size() * 4);
  for (const ClusterConfig& c : set.config) {
    double mx = 0.0, my = 0.0, mz = 0.0;
    for (int i = 0; i < 4; ++i) {
      mx += c.nucleon[i].x;
      my += c.nucleon[i].y;
      mz += c.nucleon[i].z;
    }
    mx *= 0.25;
    my *= 0.25;
    mz *= 0.25;
    for (int i = 0; i < 4; ++i) {
      const double x = c.nucleon[i].x - mx, y = c.nucleon[i].y - my;
      const double z = c.nucleon[i].z - mz;
      out.push_back(std::sqrt(x * x + y * y + z * z));
    }
  }
  return out;
}

/// chi^2 of a radial histogram against a predicted radial pdf, on `nb` equal
/// bins of width `w`.  `pred` is 4 pi rho(r) integrated over each bin by a
/// 21-point trapezoid.  Bins with fewer than 10 expected entries are dropped
/// (the Gaussian approximation to the Poisson does not hold there).
double radial_chi2(const std::vector<double>& r, double w, std::size_t nb,
                   double (*rho)(const void*, double), const void* ctx,
                   std::size_t& ndf) {
  std::vector<double> h(nb, 0.0);
  for (double v : r) {
    const std::size_t b = static_cast<std::size_t>(v / w);
    if (b < nb) h[b] += 1.0;
  }
  const double n = static_cast<double>(r.size());
  double chi2 = 0.0;
  ndf = 0;
  for (std::size_t b = 0; b < nb; ++b) {
    const double lo = static_cast<double>(b) * w;
    double integ = 0.0;
    const std::size_t ns = 21;
    const double dx = w / static_cast<double>(ns - 1);
    for (std::size_t i = 0; i < ns; ++i) {
      const double x = lo + static_cast<double>(i) * dx;
      const double f = 4.0 * kPi * x * x * rho(ctx, x);
      integ += (i == 0 || i + 1 == ns ? 0.5 : 1.0) * f * dx;
    }
    const double exp_n = n * integ;
    if (exp_n < 10.0) continue;
    const double d = h[b] - exp_n;
    chi2 += d * d / exp_n;
    ++ndf;
  }
  return chi2;
}

double rho_recentred_cb(const void* ctx, double r) {
  return static_cast<const ClusterConfigSampler*>(ctx)->rho_alpha_recentred(r);
}

double rho_gauss_cb(const void* ctx, double r) {
  const double a2 = *static_cast<const double*>(ctx);
  const double n = std::pow(2.0 * kPi * a2, 1.5);
  return std::exp(-r * r / (2.0 * a2)) / n;
}

/// The (G2)/(G3) closed forms of arXiv:2408.13213, kept HERE and only here --
/// the library implements (G1) through `theta_lm` / `clebsch_gordan`.
double g2_m1(double u, double v, double c) {
  const double c2 = c * c;
  return (4.0 * u * u - 2.0 * std::sqrt(2.0) * (1.0 - 3.0 * c2) * u * v
          + (5.0 - 3.0 * c2) * v * v) / (16.0 * kPi);
}
double g3_m0(double u, double v, double c) {
  const double c2 = c * c;
  return (2.0 * u * u + 2.0 * std::sqrt(2.0) * (1.0 - 3.0 * c2) * u * v
          + (1.0 + 3.0 * c2) * v * v) / (8.0 * kPi);
}

/// j_l for the no-i^L Fourier-Bessel cross-check of T4(iv).
double sph_j(int l, double x) {
  const double x2 = x * x;
  if (l == 0) {
    if (std::fabs(x) < 1e-3) return 1.0 - x2 / 6.0 + x2 * x2 / 120.0;
    return std::sin(x) / x;
  }
  if (std::fabs(x) < 1e-1) return x2 / 15.0 * (1.0 - x2 / 14.0 + x2 * x2 / 504.0);
  const double x3 = x2 * x;
  return (3.0 / x3 - 1.0 / x) * std::sin(x) - 3.0 / x2 * std::cos(x);
}

double fb_transform(const std::vector<double>& x, const std::vector<double>& f,
                    int l, double k) {
  double acc = 0.0;
  for (std::size_t i = 1; i < x.size(); ++i) {
    const double a = f[i - 1] * x[i - 1] * x[i - 1] * sph_j(l, k * x[i - 1]);
    const double b = f[i] * x[i] * x[i] * sph_j(l, k * x[i]);
    acc += 0.5 * (a + b) * (x[i] - x[i - 1]);
  }
  return 4.0 * kPi * acc;
}

const ClusterConfigSampler& default_sampler() {
  static const ClusterConfigSampler s{};
  return s;
}

}  // namespace

// ------------------------------------------------------------- step 2: reader

TEST_CASE("read_anl_plain reads li6.adr.fit as TWO value columns") {
  const AnlTable t = read_anl_plain(data_path(VMC_LI6_AD_FIT), 2);
  CHECK(t.col.size() == 2);
  CHECK(t.err.empty());
  CHECK(t.x.size() == 100);
  CHECK_CLOSE(t.x.front(), 0.05, 1e-12);
  CHECK_CLOSE(t.x.back(), 9.95, 1e-12);
  // The second VALUE column is R2LI6FIT, not an error bar: read_anl_momentum
  // would have landed it in err[0] because the rows also carry 3 numbers.
  CHECK_CLOSE(t.col[0].front(), 0.56254, 1e-12);
  CHECK_CLOSE(t.col[1].front(), -5.9352e-5, 1e-12);
  CHECK_CLOSE(t.col[0].back(), -0.0039492, 1e-12);

  // he4.density DOES go through the existing momentum reader, unchanged.
  const std::vector<AnlTable> d = read_anl_momentum(data_path(VMC_HE4_DENSITY));
  CHECK(d.size() >= 1);
  CHECK(d[0].col.size() == 1);
  CHECK(d[0].err.size() == 1);
  CHECK_CLOSE(d[0].x.front(), 0.05, 1e-12);
  CHECK_CLOSE(d[0].col[0].front(), 0.1746, 1e-12);
  CHECK_CLOSE(d[0].err[0].front(), 0.00644, 1e-12);
}

// ------------------------------------------------------------------ T1 (G1)

TEST_CASE("T1 (G1) through theta_lm and clebsch_gordan IS Eqs. (7)-(8)") {
  const double pairs[3][2] = {{0.37, -0.11}, {1.0, 0.0}, {-0.2, 0.45}};
  const double cs[7] = {-1.0, -0.6, -0.2, 0.0, 0.35, 0.8, 1.0};
  for (const auto& p : pairs) {
    for (double c : cs) {
      CHECK_CLOSE(cluster_density(p[0], p[1], +1, c), g2_m1(p[0], p[1], c), 1e-14);
      CHECK_CLOSE(cluster_density(p[0], p[1], -1, c), g2_m1(p[0], p[1], c), 1e-14);
      CHECK_CLOSE(cluster_density(p[0], p[1], 0, c), g3_m0(p[0], p[1], c), 1e-14);
    }
  }
  // (G1) has no phi dependence and is normalized: sum over m of the density
  // integrates to 3 with unit radial norms.
  const ClusterConfigSampler& s = default_sampler();
  CHECK_CLOSE(s.rho_np(2.0, 0.3, 1.0), s.rho_np(2.0, 0.3, 1.0), 1e-15);
}

// -------------------------------------------------------------------- T2, T3

TEST_CASE("T2 fdeut.av18 re-read reproduces its own header") {
  const ClusterConfigSampler& s = default_sampler();
  const RadialMoments m = radial_moments(s.np_grid(), s.np_wave(0), s.np_wave(2));
  // 5e-5, not 1e-5: the header numbers come from the Fortran code's own
  // h = 0.005 quadrature, not from re-integrating the printed h = 0.01 table.
  CHECK_CLOSE(m.norm(), 1.0, 5e-5);
  CHECK_CLOSE(m.p_d(), 0.057599, 5e-5);
  CHECK_CLOSE(0.25 * m.quadrupole() / m.norm(), 0.269673, 5e-5);
  CHECK_CLOSE(0.5 * std::sqrt(m.r2 / m.norm()), 1.967364, 5e-5);
  // The sampler's own copies of the same two numbers.
  CHECK_CLOSE(s.r2_np_fm2(), 15.481685, 1e-6);
}

TEST_CASE("T3 li6.ad r-block reproduces the file's printed norm split") {
  const ClusterConfigSampler& s = default_sampler();
  const RadialMoments m = radial_moments(s.alpha_d_grid(), s.alpha_d_wave(0),
                                         s.alpha_d_wave(2));
  const double n0 = s.s_alpha_d() * m.n0 / m.norm();
  const double n2 = s.s_alpha_d() * m.n2 / m.norm();
  CHECK_CLOSE_AT(n0, 0.838, 0.0, 2e-3);
  CHECK_CLOSE_AT(n2, 0.017, 0.0, 2e-3);
  CHECK_CLOSE_AT(n0 + n2, 0.856, 0.0, 2e-3);
  CHECK_CLOSE(s.s_alpha_d(), 0.854232, 1e-5);
}

// ------------------------------------------------------------------------ T4

TEST_CASE("T4 the r-space sign conventions") {
  const ClusterConfigSampler& s = default_sampler();
  const std::vector<double>& x = s.alpha_d_grid();
  std::vector<double> f0 = s.alpha_d_wave(0), f2 = s.alpha_d_wave(2);
  const RadialMoments ref = radial_moments(x, f0, f2);

  SUBCASE("(i) negating BOTH waves is the unobservable global phase") {
    std::vector<double> g0 = f0, g2 = f2;
    for (double& v : g0) v = -v;
    for (double& v : g2) v = -v;
    const RadialMoments m = radial_moments(x, g0, g2);
    CHECK(m.q_int == ref.q_int);
    CHECK(m.q_dd == ref.q_dd);
    CHECK(m.r2 == ref.r2);
    CHECK(m.n0 == ref.n0);
    // and the density itself is untouched
    for (double c : {-0.7, 0.0, 0.4}) {
      CHECK(cluster_amp2(-f0[40], -f2[40], 1, 0, c)
            == cluster_amp2(f0[40], f2[40], 1, 0, c));
    }
  }
  SUBCASE("(ii) negating ONLY the D wave negates the interference term") {
    std::vector<double> g2 = f2;
    for (double& v : g2) v = -v;
    const RadialMoments m = radial_moments(x, f0, g2);
    CHECK(m.q_int == -ref.q_int);
    CHECK(m.q_dd == ref.q_dd);
    CHECK((m.quadrupole() > 0.0) != (ref.quadrupole() > 0.0));
  }
  SUBCASE("(iii) the only OBSERVABLE check: Q has the measured sign") {
    CHECK(s.q_matter_analytic_fm2(+1) < 0.0);
    CHECK(s.q_matter_analytic_fm2(-1) < 0.0);
    CHECK(LI6_QUADRUPOLE_FM2 < 0.0);
    CHECK(s.q_matter_analytic_fm2(0) > 0.0);
    // the global phase is fixed to R_0(R -> infinity) > 0
    CHECK(s.alpha_d_wave(0).back() > 0.0);
  }
  SUBCASE("(iv) the NO-i^L transform reproduces li6.ad's k-block structure") {
    // The physical psi_L carries (-i)^L, so the PHYSICAL k-space S-D
    // interference has the OPPOSITE sign to this one.  Do NOT assert sign
    // identity with the physical amplitude; assert the k-block's structure.
    double node = 0.0, prev_k = 0.05, prev = fb_transform(x, f0, 0, 0.05);
    for (double k = 0.06; k < 1.2; k += 0.005) {
      const double v = fb_transform(x, f0, 0, k);
      if ((v > 0.0) != (prev > 0.0)) {
        node = prev_k - prev * (k - prev_k) / (v - prev);
        break;
      }
      prev = v;
      prev_k = k;
    }
    CHECK_CLOSE_AT(node, 0.678, 0.0, 0.01);
    for (double k : {0.05, 0.2, 0.4, 0.5}) {
      CHECK(fb_transform(x, f2, 2, k) / fb_transform(x, f0, 0, k) < 0.0);
    }
    CHECK(fb_transform(x, f2, 2, 1.0) / fb_transform(x, f0, 0, 1.0) > 0.0);
  }
}

// ------------------------------------------------------------------ T5, T5g

TEST_CASE("T5 the recentred alpha core") {
  const ClusterConfigSampler& s = default_sampler();
  const ClusterConfigSet set = s.sample_set(25000, 1, 5, 0);
  const std::vector<double> r = core_radii(set);
  // (a) <s^2> closes EXACTLY for any source shape: recentring gives
  //     <s^2> = (3/4)<v^2> and lambda^2 = 4/3 puts it back on the table's.
  //     The MC gate is 5 * sd / sqrt(n_config) of the PER-CONFIGURATION mean
  //     of |s_i|^2 over the four core nucleons, like every other gate here.
  const double n5 = static_cast<double>(set.config.size());
  double s2 = 0.0, s2sq = 0.0;
  for (std::size_t k = 0; k < set.config.size(); ++k) {
    double mean = 0.0;
    for (std::size_t i = 0; i < 4; ++i) {
      mean += r[4 * k + i] * r[4 * k + i];
    }
    mean *= 0.25;
    s2 += mean;
    s2sq += mean * mean;
  }
  s2 /= n5;
  const double sd_s2 =
      std::sqrt(std::fmax(s2sq - n5 * s2 * s2, 0.0) / (n5 - 1.0));
  const double gate_s2 = 5.0 * sd_s2 / std::sqrt(n5);
  MESSAGE("T5(a) sampled <s^2> = " << s2 << " vs analytic "
          << s.r2_alpha_fm2() << " (5 sigma " << gate_s2 << ")");
  CHECK_CLOSE_AT(s2, s.r2_alpha_fm2(), 0.0, gate_s2);
  CHECK_CLOSE(s.r2_alpha_fm2(), 2.0774, 1e-4);
  // (b) THE code gate: the sampled histogram against the CLOSED-FORM
  //     recentred prediction ftilde(3q/4) ftilde(q/4)^3, which is what the
  //     sampler is supposed to produce.
  std::size_t ndf = 0;
  const double chi2 = radial_chi2(r, 0.15, 40, rho_recentred_cb, &s, ndf);
  MESSAGE("T5(b) chi2/ndf vs the closed-form recentred density: "
          << chi2 / static_cast<double>(ndf) << " on " << ndf << " bins");
  CHECK(chi2 / static_cast<double>(ndf) < 2.0);
  // (c) RECORDED, not gated: recentring does NOT reproduce he4.density.
  const std::vector<AnlTable> d = read_anl_momentum(data_path(VMC_HE4_DENSITY));
  double norm = 0.0;
  for (std::size_t i = 1; i < d[0].x.size(); ++i) {
    norm += 0.5 * (d[0].col[0][i] * d[0].x[i] * d[0].x[i]
                   + d[0].col[0][i - 1] * d[0].x[i - 1] * d[0].x[i - 1])
            * (d[0].x[i] - d[0].x[i - 1]);
  }
  norm *= 4.0 * kPi;
  double chi2_tab = 0.0;
  std::size_t nb_tab = 0;
  for (std::size_t b = 0; b < 40; ++b) {
    const double lo = static_cast<double>(b) * 0.15;
    double integ = 0.0;
    for (std::size_t i = 0; i < 21; ++i) {
      const double xx = lo + static_cast<double>(i) * 0.15 / 20.0;
      const double rho = (xx < d[0].x.front() || xx > d[0].x.back())
                             ? 0.0
                             : np_interp(xx, d[0].x, d[0].col[0]) / norm;
      integ += (i == 0 || i == 20 ? 0.5 : 1.0) * 4.0 * kPi * xx * xx * rho
               * (0.15 / 20.0);
    }
    const double e = static_cast<double>(r.size()) * integ;
    if (e < 10.0) continue;
    double obs = 0.0;
    for (double v : r) {
      if (v >= lo && v < lo + 0.15) obs += 1.0;
    }
    chi2_tab += (obs - e) * (obs - e) / e;
    ++nb_tab;
  }
  MESSAGE("T5(c) chi2/ndf of the recentred core against he4.density itself: "
          << chi2_tab / static_cast<double>(nb_tab) << " on " << nb_tab
          << " bins -- RECORDED, NOT GATED (the shape is not preserved by "
             "recentring; only <s^2> is)");
}

TEST_CASE("T5g the Gaussian core, where the lambda inflation is exact in SHAPE") {
  ClusterConfigOptions o;
  o.alpha_source = AlphaCoreSource::Gaussian;
  const ClusterConfigSampler s(o);
  const ClusterConfigSet set = s.sample_set(25000, 1, 11, 0);
  const std::vector<double> r = core_radii(set);
  double a2 = cluster_point_a2_fm2(2, 4);
  std::size_t ndf = 0;
  const double chi2 = radial_chi2(r, 0.15, 40, rho_gauss_cb, &a2, ndf);
  MESSAGE("T5g chi2/ndf vs the analytic Gaussian: "
          << chi2 / static_cast<double>(ndf) << " on " << ndf << " bins");
  CHECK(chi2 / static_cast<double>(ndf) < 1.5);
  CHECK_CLOSE(s.r2_alpha_fm2(), 3.0 * a2, 1e-6);
}

// ----------------------------------------------------------------------- T5h

TEST_CASE("T5h the hard core: the analytic layer follows the sampler") {
  // `min_nn_separation_fm > 0` correlates the four s_i, so the closed form
  // <s^2> = (3/4)<v^2> is no longer the sampler's own second moment.  The
  // constructor measures it instead; this is the gate that the ANALYTIC layer
  // (and therefore (G6), match_li6_radius and the sidecar) agrees with what
  // the sampler actually draws.
  ClusterConfigOptions o;
  o.min_nn_separation_fm = 0.9;            // arXiv:2605.00454's value
  const ClusterConfigSampler s(o);
  const ClusterConfigSampler& free = default_sampler();
  // it MOVES: the hard core is not a no-op on the analytic layer
  CHECK(s.r2_alpha_fm2() > free.r2_alpha_fm2() + 0.05);
  MESSAGE("T5h <s^2> = " << s.r2_alpha_fm2() << " fm^2 (hard core 0.9 fm) vs "
          << free.r2_alpha_fm2() << " free; <r^2> " << s.r2_analytic_fm2()
          << " vs " << free.r2_analytic_fm2());
  const ClusterConfigSet set = s.sample_set(50000, 1, 7, 0);
  const double n = static_cast<double>(set.config.size());
  const double gate = 5.0 * set.r2_sd_fm2 / std::sqrt(n);
  CHECK_CLOSE_AT(set.r2_mean_fm2, s.r2_analytic_fm2(), 0.0, gate);
  // and the independent-draw closed form would have FAILED this same set
  CHECK(std::fabs(set.r2_mean_fm2 - free.r2_analytic_fm2()) > gate);
  // the separation is actually enforced
  double dmin = 1e9;
  for (std::size_t k = 0; k < 200; ++k) {
    for (int i = 0; i < 4; ++i) {
      for (int jj = i + 1; jj < 4; ++jj) {
        const double dx = set.config[k].nucleon[i].x - set.config[k].nucleon[jj].x;
        const double dy = set.config[k].nucleon[i].y - set.config[k].nucleon[jj].y;
        const double dz = set.config[k].nucleon[i].z - set.config[k].nucleon[jj].z;
        dmin = std::fmin(dmin, std::sqrt(dx * dx + dy * dy + dz * dz));
      }
    }
  }
  CHECK(dmin >= 0.9);
  // the provenance says the number is numerical, and the run is reproducible
  CHECK(s.provenance().find("alpha <s^2> NUMERICAL (hard core)")
        != std::string::npos);
  const ClusterConfigSampler s2(o);
  CHECK_CLOSE(s2.r2_alpha_fm2(), s.r2_alpha_fm2(), 1e-15);
}

// ------------------------------------------------------------------------ T6

TEST_CASE("T6 the assembled 6Li one-body density vs li6.density") {
  const ClusterConfigSampler& s = default_sampler();
  const ClusterConfigSet set = s.sample_set(100000, 1, 3, 0);
  const double n = static_cast<double>(set.config.size());
  const double rms = std::sqrt(set.r2_mean_fm2);
  const double sig5 = 5.0 * set.r2_sd_fm2 / std::sqrt(n) / (2.0 * rms);
  CHECK_CLOSE_AT(rms, std::sqrt(s.r2_analytic_fm2()), 0.0, sig5);
  // The model is ~4 % LARGER than the ab-initio one-body density.  This gates
  // the DOCUMENTED discrepancy, not agreement.
  const double ratio = rms / LI6_R_POINT_VMC_FM;
  MESSAGE("T6 r_rms(sampled) = " << rms << " fm, ratio to li6.density's "
          << LI6_R_POINT_VMC_FM << " = " << ratio);
  CHECK_CLOSE_AT(ratio, 1.039, 0.0, 0.005);
  CHECK(ratio > 1.0);
}

// ------------------------------------------------------------------ T7, T7a

TEST_CASE("T7 <r^2> closure: the sampled set reproduces (G6)") {
  const ClusterConfigSampler& s = default_sampler();
  const ClusterConfigSet set = s.sample_set(100000, 1, 20260902, 0);
  const double n = static_cast<double>(set.config.size());
  const double gate = 5.0 * set.r2_sd_fm2 / std::sqrt(n);
  MESSAGE("T7 sampled <r^2> = " << set.r2_mean_fm2 << " vs (G6) "
          << s.r2_analytic_fm2() << ", 5 sigma = " << gate);
  CHECK_CLOSE_AT(set.r2_mean_fm2, s.r2_analytic_fm2(), 0.0, gate);
  CHECK_CLOSE(set.r2_sd_fm2, 3.66, 0.05);
}

TEST_CASE("T7a <r^2> on the sampler's own GRID, no Monte Carlo") {
  const ClusterConfigSampler& s = default_sampler();
  // 2e-3, not the design's 1e-6.  The analytic layer is the trapezoid on the
  // TABLES' own 0.1 fm nodes (which is what reproduces every sec. 8 number),
  // while the cell sum converges to the exact integral of the piecewise-linear
  // interpolant.  The two differ by the trapezoid's own O(h^2) error on f^2,
  // -(h^2/6) integral (f')^2 x^2 dx ~ 6e-4 relative, and that gap does NOT
  // shrink with n_r (measured: identical at n_r = 512 and 4096).
  MESSAGE("T7a grid/analytic - 1 = "
          << s.r2_grid_fm2() / s.r2_analytic_fm2() - 1.0);
  CHECK_CLOSE(s.r2_grid_fm2(), s.r2_analytic_fm2(), 2e-3);
}

// ------------------------------------------------------------------ T8, T8a

TEST_CASE("T8 the quadrupole: the sampled set reproduces (G5)") {
  const ClusterConfigSampler& s = default_sampler();
  const double n = 100000.0;
  double q_p1 = 0.0;
  for (int m : {1, 0, -1}) {
    const ClusterConfigSet set = s.sample_set(100000, m, 4242, 0);
    const double gate = 5.0 * set.q_matter_sd_fm2 / std::sqrt(n);
    MESSAGE("T8 m = " << m << ": sampled Q_matter = " << set.q_matter_fm2
            << " vs (G5) " << s.q_matter_analytic_fm2(m) << ", 5 sigma = "
            << gate);
    CHECK_CLOSE_AT(set.q_matter_fm2, s.q_matter_analytic_fm2(m), 0.0, gate);
    if (m == 1) q_p1 = s.q_matter_analytic_fm2(m);
  }
  // (3 m^2 - 2): the m = 0 state is exactly -2 x the m = +-1 one, and the
  // equal-thirds mixture is zero.
  CHECK_CLOSE(s.q_matter_analytic_fm2(0), -2.0 * q_p1, 1e-14);
  CHECK_CLOSE(s.q_matter_analytic_fm2(-1), q_p1, 1e-14);
  const ClusterConfigSet mix = s.sample_set_unpolarized(100000, 4242, 0);
  const double gate = 5.0 * mix.q_matter_sd_fm2 / std::sqrt(n);
  CHECK_CLOSE_AT(mix.q_matter_fm2, 0.0, 0.0, gate);
  CHECK_CLOSE(s.q_matter_analytic_fm2(1), -1.2309, 1e-4);
}

TEST_CASE("T8a the quadrupole on the sampler's own GRID, no Monte Carlo") {
  const ClusterConfigSampler& s = default_sampler();
  // 5e-3 for the same reason as T7a, plus the cos-theta midpoint rule's own
  // O(h_c^2): 3.0e-3 at n_c = 96, falling to 8.6e-4 (the radial floor) by
  // n_c = 3072.  Both are quadrature, not model, differences.
  MESSAGE("T8a grid/analytic - 1 = "
          << s.q_matter_grid_fm2(1) / s.q_matter_analytic_fm2(1) - 1.0);
  CHECK_CLOSE(s.q_matter_grid_fm2(1), s.q_matter_analytic_fm2(1), 5e-3);
  CHECK_CLOSE(s.q_matter_grid_fm2(0), s.q_matter_analytic_fm2(0), 5e-3);
  ClusterConfigOptions o;
  o.n_c = 3072;
  const ClusterConfigSampler fine(o);
  CHECK_CLOSE(fine.q_matter_grid_fm2(1), fine.q_matter_analytic_fm2(1), 1.5e-3);
}

// ------------------------------------------------------------------------ T9

TEST_CASE("T9 the coherent.hpp bridge: a2_from_quadrupole IS the published a2") {
  // Zero free parameters: the deuteron's own quadrupole through (G7) + (G8).
  for (const MantysaariRow& r : mantysaari_a2_deuteron()) {
    const double p1 = a2_from_quadrupole(2.0 * kQdFm2, 2, r.t_abs, 1);
    const double p0 = a2_from_quadrupole(2.0 * kQdFm2, 2, r.t_abs, 0);
    CHECK_CLOSE(p1, r.a2_m1, 0.10);
    CHECK_CLOSE(p0, r.a2_m0, 0.25);
    CHECK_CLOSE(p0, -2.0 * p1, 1e-14);
  }
  // and the 6Li expectation the same map gives
  const ClusterConfigSampler& s = default_sampler();
  CHECK_CLOSE(s.a2_from_geometry(0.3, 1), 0.1976, 1e-3);
  CHECK(s.a2_from_geometry(0.3, 1) > 0.0);   // OPPOSITE sign to the deuteron's
  CHECK(mantysaari_a2_deuteron().back().a2_m1 < 0.0);
  CHECK_CLOSE(s.eps_b0_equivalent(), -0.0506, 1e-3);
  CHECK_CLOSE(gaussian_slope(std::sqrt(LI6_R2_POINT_FM2)), 52.0368, 1e-5);
}

// ----------------------------------------------------------------------- T10

TEST_CASE("T10 D_T: the exact CG identity, and TaggedModel at the SAME P_D") {
  const ClusterConfigSampler& s = default_sampler();
  CHECK_CLOSE(s.tensor_dilution(), 1.0 - 0.9 * s.p_d_alpha_d(), 1e-12);
  CHECK_CLOSE(s.p_d_alpha_d(), 0.020112, 1e-5);
  CHECK_CLOSE(s.tensor_dilution(), 0.98190, 1e-5);
  // P(m_S | M = +1) = N_0 + N_2/10, 3 N_2/10, 6 N_2/10.
  const double pd = s.p_d_alpha_d();
  CHECK_CLOSE(s.p_ms(1, 1.0), 1.0 - pd + 0.1 * pd, 1e-12);
  CHECK_CLOSE(s.p_ms(1, 0.0), 0.3 * pd, 1e-12);
  CHECK_CLOSE(s.p_ms(1, -1.0), 0.6 * pd, 1e-12);
  // The Hulthen path, where P_D is still a knob.  NOT VmcAV18, which pins
  // P_D at VMC_P_D_LI6 = 0.01935 against this sampler's 0.02011.
  const TaggedModel tm(li6_alpha_channel(BETA_DEFAULT, s.p_d_alpha_d()));
  // 1e-4, the precedent constants.hpp sets for the same comparison:
  // population_integrated is a 280 x 96 grid quadrature, not a closed form.
  CHECK_CLOSE(tm.tensor_dilution(), s.tensor_dilution(), 1e-4);
  CHECK(tm.tensor_dilution() != s.tensor_dilution());
}

// ----------------------------------------------------------------------- T11

TEST_CASE("T11 the unpolarized set is isotropic") {
  const ClusterConfigSampler& s = default_sampler();
  const ClusterConfigSet set = s.sample_set_unpolarized(100000, 808, 0);
  const double n = static_cast<double>(set.config.size());
  // equal thirds, counts differing by at most one
  int cnt[3] = {0, 0, 0};
  double sx = 0.0, sy = 0.0, sz = 0.0, p2 = 0.0, p2sq = 0.0;
  for (const ClusterConfig& c : set.config) {
    cnt[1 - c.m_ion] += 1;
    double pc = 0.0;
    for (int i = 0; i < 6; ++i) {
      const double x = c.nucleon[i].x, y = c.nucleon[i].y, z = c.nucleon[i].z;
      sx += x * x;
      sy += y * y;
      sz += z * z;
      const double rr = x * x + y * y + z * z;
      pc += 0.5 * (3.0 * z * z / rr - 1.0);
    }
    pc /= 6.0;
    p2 += pc;
    p2sq += pc * pc;
  }
  // <P2> gets its OWN 5 sigma gate from the per-configuration spread: a
  // hand-typed 0.01 here is ~13 sigma, and a set that was only ~60 %
  // unpolarized (|<P2>| ~ 0.016 when fully m = +1) would still pass it.
  const double sd_p2 = std::sqrt(
      std::fmax(p2sq - p2 * p2 / n, 0.0) / (n - 1.0));
  const double gate_p2 = 5.0 * sd_p2 / std::sqrt(n);
  CHECK(std::abs(cnt[0] - cnt[1]) <= 1);
  CHECK(std::abs(cnt[1] - cnt[2]) <= 1);
  const double k = 1.0 / (6.0 * n);
  sx *= k;
  sy *= k;
  sz *= k;
  p2 *= k;
  const double gate = 5.0 * set.delta_perp_sd_fm2 / std::sqrt(n);
  MESSAGE("T11 <x^2> " << sx << " <y^2> " << sy << " <z^2> " << sz
          << " <P2> " << p2 << " gate " << gate << " P2 gate " << gate_p2);
  CHECK_CLOSE_AT(sx, sy, 0.0, gate);
  CHECK_CLOSE_AT(sy, sz, 0.0, gate);
  CHECK_CLOSE_AT(p2, 0.0, 0.0, gate_p2);
  CHECK_CLOSE_AT(set.delta_perp_fm2, 0.0, 0.0, gate);
}

// ------------------------------------------------------------------ T12, T13

TEST_CASE("T12 every configuration has sum_i r_i = 0") {
  const ClusterConfigSampler& s = default_sampler();
  const ClusterConfigSet set = s.sample_set(20000, 1, 17, 0);
  CHECK(set.cm[0] < 1e-12);
  double worst = 0.0;
  for (const ClusterConfig& c : set.config) {
    double cx = 0.0, cy = 0.0, cz = 0.0;
    for (int i = 0; i < 6; ++i) {
      cx += c.nucleon[i].x;
      cy += c.nucleon[i].y;
      cz += c.nucleon[i].z;
    }
    worst = std::fmax(worst, std::fmax(std::fabs(cx),
                                       std::fmax(std::fabs(cy), std::fabs(cz))));
  }
  MESSAGE("T12 worst |sum_i r_i| component = " << worst << " fm");
  CHECK(worst < 1e-12);
}

TEST_CASE("T13 reproducibility: order, N and per-configuration streams") {
  const ClusterConfigSampler& s = default_sampler();
  const ClusterConfigSet a = s.sample_set(200, 1, 12345, 7);
  const ClusterConfigSet b = s.sample_set(200, 1, 12345, 7);
  const ClusterConfigSet c = s.sample_set(400, 1, 12345, 7);
  for (std::size_t i = 0; i < a.config.size(); ++i) {
    for (int k = 0; k < 6; ++k) {
      CHECK(a.config[i].nucleon[k].x == b.config[i].nucleon[k].x);
      CHECK(a.config[i].nucleon[k].y == c.config[i].nucleon[k].y);
      CHECK(a.config[i].nucleon[k].z == c.config[i].nucleon[k].z);
    }
  }
  // config i is bit-identical when drawn on its own from the header's stream
  for (std::size_t i : {std::size_t(0), std::size_t(3), std::size_t(199)}) {
    Rng rng(12345, 7, kConfigStream, static_cast<std::uint64_t>(i));
    const ClusterConfig one = s.sample(rng, 1);
    for (int k = 0; k < 6; ++k) {
      CHECK(one.nucleon[k].x == a.config[i].nucleon[k].x);
      CHECK(one.nucleon[k].y == a.config[i].nucleon[k].y);
      CHECK(one.nucleon[k].z == a.config[i].nucleon[k].z);
    }
    CHECK(one.m_s == a.config[i].m_s);
  }
  CHECK(kConfigStream == 0x434F4E464947ull);
}

TEST_CASE("T14 a different run is a different stream, same distribution") {
  const ClusterConfigSampler& s = default_sampler();
  const ClusterConfigSet a = s.sample_set(50000, 1, 5150, 0);
  const ClusterConfigSet b = s.sample_set(50000, 1, 5150, 1);
  for (int k = 0; k < 6; ++k) {
    CHECK(a.config[0].nucleon[k].x != b.config[0].nucleon[k].x);
    CHECK(a.config[0].nucleon[k].y != b.config[0].nucleon[k].y);
    CHECK(a.config[0].nucleon[k].z != b.config[0].nucleon[k].z);
  }
  const double sig = std::sqrt(a.r2_sd_fm2 * a.r2_sd_fm2
                               + b.r2_sd_fm2 * b.r2_sd_fm2)
                     / std::sqrt(50000.0);
  CHECK_CLOSE_AT(a.r2_mean_fm2, b.r2_mean_fm2, 0.0, 5.0 * sig);
}

// ------------------------------------------------------------------ T15, T16

TEST_CASE("T15 the writer round-trips through a he3.dat-style reader") {
  const ClusterConfigSampler& s = default_sampler();
  const ClusterConfigSet set = s.sample_set(50, 1, 909, 0);
  const std::string path = "test_cluster_config_rt.dat";
  const std::size_t nrow = write_snd_configs(set, s, path);
  CHECK(nrow == 50);
  std::ifstream in(path);
  REQUIRE(in.good());
  std::string line;
  std::size_t row = 0;
  while (std::getline(in, line)) {
    CHECK(line.find('#') == std::string::npos);
    CHECK(!line.empty());
    // exactly 25 whitespace-separated fields
    std::istringstream fs(line);
    std::vector<double> f;
    double v = 0.0;
    while (fs >> v) f.push_back(v);
    CHECK(f.size() == 25);
    // the consumer's own read: eighteen bare `>>` doubles
    for (int i = 0; i < 6; ++i) {
      CHECK_CLOSE(f[3 * i + 0], set.config[row].nucleon[i].x, 1e-12);
      CHECK_CLOSE(f[3 * i + 1], set.config[row].nucleon[i].y, 1e-12);
      CHECK_CLOSE(f[3 * i + 2], set.config[row].nucleon[i].z, 1e-12);
      CHECK(f[18 + i] == static_cast<double>(set.config[row].nucleon[i].isospin));
    }
    CHECK(f[24] == 1.0);
    ++row;
  }
  CHECK(row == 50);
  in.close();
  std::remove(path.c_str());
  std::remove((path + ".meta.json").c_str());
}

TEST_CASE("T16 the sidecar records every option, and md5/git as null") {
  const ClusterConfigSampler& s = default_sampler();
  const ClusterConfigSet set = s.sample_set(5, 1, 4, 0);
  const std::string path = "test_cluster_config_meta.dat";
  write_snd_configs(set, s, path);
  jsonmin::Value j;
  REQUIRE(jsonmin::load_file(path + ".meta.json", j));
  REQUIRE(j.is_object());
  const jsonmin::Object& o = j.obj();
  REQUIRE(o.count("options") == 1);
  const jsonmin::Object& opt = o.at("options").obj();
  // The hand-maintained list: adding an option without recording it fails.
  const char* fields[13] = {"alpha_source", "alpha_d_source", "theta_s",
                            "phi_s", "alpha_cm_inflate", "min_nn_separation_fm",
                            "alpha_d_scale", "quadrupole_target_fm2",
                            "exact_coherence", "n_r", "n_c", "r_max_fm",
                            "rnp_max_fm"};
  for (const char* f : fields) {
    CHECK_MESSAGE(opt.count(f) == 1, "sidecar options is missing " << f);
  }
  CHECK(opt.count("quadrupole_dial_s") == 1);
  CHECK(o.count("git") == 1);
  CHECK(o.at("git").is_null());
  REQUIRE(o.count("inputs") == 1);
  const jsonmin::Array& in = o.at("inputs").arr();
  CHECK(in.size() == 4);
  for (const jsonmin::Value& v : in) {
    CHECK(v.obj().count("md5") == 1);
    CHECK(v.obj().at("md5").is_null());
    CHECK(v.obj().at("bytes").num() > 0.0);
  }
  const jsonmin::Object& band = o.at("quadrupole_band_fm2").obj();
  CHECK_CLOSE(band.at("measured").num(), LI6_QUADRUPOLE_FM2, 1e-12);
  CHECK_CLOSE(band.at("gfmc_av18_il7").num(), LI6_QUADRUPOLE_GFMC_FM2, 1e-12);
  CHECK_CLOSE(band.at("alpha_d_model").num(),
              0.5 * s.q_matter_analytic_fm2(1), 1e-12);
  CHECK(o.count("caveats") == 1);
  // "recomputed_norm" is RECOMPUTED from li6.adr.fit at write time, not a
  // literal typed into the writer (docs/CONVENTIONS.md: one home per number).
  const AnlTable fit = read_anl_plain(data_path(VMC_LI6_AD_FIT), 2);
  const RadialMoments fm = radial_moments(fit.x, fit.col[0], fit.col[1]);
  const jsonmin::Object& e0 = in.front().obj();
  CHECK(e0.at("file").str() == std::string(VMC_LI6_AD_FIT));
  const jsonmin::Object& rn = e0.at("recomputed_norm").obj();
  CHECK_CLOSE(rn.at("integral").num(), fm.norm(), 1e-12);
  CHECK_CLOSE(rn.at("p_d").num(), fm.p_d(), 1e-12);
  for (const jsonmin::Value& v : in) CHECK(v.obj().count("used") == 1);
  std::remove(path.c_str());
  std::remove((path + ".meta.json").c_str());
}

TEST_CASE("T16b the sidecar's caveats and inputs follow the options") {
  ClusterConfigOptions o;
  o.quadrupole_target_fm2 = LI6_QUADRUPOLE_FM2;
  o.alpha_source = AlphaCoreSource::Gaussian;
  o.alpha_d_source = AlphaDSource::OverlapRaw;
  const ClusterConfigSampler s(o);
  const ClusterConfigSet set = s.sample_set(5, 1, 4, 0);
  const std::string path = "test_cluster_config_meta_dialed.dat";
  write_snd_configs(set, s, path);
  jsonmin::Value j;
  REQUIRE(jsonmin::load_file(path + ".meta.json", j));
  const jsonmin::Object& obj = j.obj();
  // only the tables THIS option set reads are marked used
  for (const jsonmin::Value& v : obj.at("inputs").arr()) {
    const std::string f = v.obj().at("file").str();
    const bool used = v.obj().at("used").boolean();
    if (f == std::string(VMC_LI6_AD_FIT)) CHECK(!used);
    if (f == std::string(VMC_HE4_DENSITY)) CHECK(!used);
    if (f == std::string(VMC_LI6_OVERLAP)) CHECK(used);
    if (f == std::string(VMC_DEUTERON_WAVE)) CHECK(used);
  }
  // the overshoot factor is COMPUTED from this run's own band -- with the
  // dial on the measured value it is 1.0, not the natural geometry's 7.5 --
  // and the fixed "-0.615..-0.730" source range is gone.
  std::string all;
  for (const jsonmin::Value& v : obj.at("caveats").arr()) all += v.str() + " ";
  MESSAGE("T16b caveats: " << all);
  CHECK(all.find("|model/measured| = 1.0") != std::string::npos);
  CHECK(all.find("DIALLED/SCALED") != std::string::npos);
  CHECK(all.find("-0.615..-0.730") == std::string::npos);
  CHECK(all.find("~7.5") == std::string::npos);
  std::remove(path.c_str());
  std::remove((path + ".meta.json").c_str());
}

// ----------------------------------------------------------- T17, T17b, T18

TEST_CASE("T17 the quadrupole dial closes on the analytic layer") {
  ClusterConfigOptions o;
  o.quadrupole_target_fm2 = LI6_QUADRUPOLE_FM2;
  const ClusterConfigSampler s(o);
  CHECK_CLOSE(0.5 * s.q_matter_analytic_fm2(1), LI6_QUADRUPOLE_FM2, 1e-9);
  CHECK_CLOSE(s.quadrupole_dial_s(), 0.4032, 1e-3);
  CHECK_CLOSE(s.p_d_alpha_d(), 3.325e-3, 1e-3);
  // it is NOT the sqrt(target/model) rule, which lands at -0.048
  const ClusterConfigSampler& base = default_sampler();
  const double naive = std::sqrt(std::fabs(LI6_QUADRUPOLE_FM2
                                / (0.5 * base.q_matter_analytic_fm2(1))));
  CHECK_CLOSE(naive, 0.3646, 1e-3);
  CHECK(std::fabs(naive - s.quadrupole_dial_s()) > 0.03);
  // the band's third entry moves with the dial and is never a literal
  CHECK_CLOSE(s.quadrupole_band_fm2()[2], LI6_QUADRUPOLE_FM2, 1e-9);
  CHECK_CLOSE(base.quadrupole_band_fm2()[2], -0.6154, 1e-3);
}

TEST_CASE("T17b the dial's reachable range is enforced") {
  ClusterConfigOptions o;
  o.quadrupole_target_fm2 = 0.35;          // above the s = 0 floor, +Q_d
  CHECK_THROWS_AS(ClusterConfigSampler{o}, std::runtime_error);
  o.quadrupole_target_fm2 = -0.9;          // below the s = 1 value
  CHECK_THROWS_AS(ClusterConfigSampler{o}, std::runtime_error);
  // the floor itself: at s = 0 the deuteron's own +0.270 fm^2 survives
  o.quadrupole_target_fm2 = 0.2696;
  const ClusterConfigSampler at_floor(o);
  CHECK(at_floor.quadrupole_dial_s() < 0.01);
  CHECK_CLOSE(0.5 * at_floor.q_matter_analytic_fm2(1), 0.2696, 1e-9);
}

TEST_CASE("T18 match_li6_radius closes on the sampled radius") {
  const ClusterConfigSampler& base = default_sampler();
  ClusterConfigOptions o;
  o.alpha_d_scale = base.match_li6_radius();
  MESSAGE("T18 alpha_d_scale = " << o.alpha_d_scale);
  const ClusterConfigSampler s(o);
  CHECK_CLOSE(std::sqrt(s.r2_analytic_fm2()), LI6_R_POINT_VMC_FM, 1e-12);
  const ClusterConfigSet set = s.sample_set(100000, 1, 606, 0);
  const double rms = std::sqrt(set.r2_mean_fm2);
  const double gate = 5.0 * set.r2_sd_fm2 / std::sqrt(100000.0) / (2.0 * rms);
  CHECK_CLOSE_AT(rms, LI6_R_POINT_VMC_FM, 0.0, gate);
}

// ----------------------------------------------------------------------- T19

TEST_CASE("T19 the option combinations that must throw") {
  ClusterConfigOptions o;
  o.exact_coherence = true;
  CHECK_THROWS_AS(o.validate(), std::runtime_error);
  CHECK_THROWS_AS(ClusterConfigSampler{o}, std::runtime_error);
  o = ClusterConfigOptions{};
  o.n_r = 4;
  CHECK_THROWS_AS(o.validate(), std::runtime_error);
  o = ClusterConfigOptions{};
  o.n_c = 2;
  CHECK_THROWS_AS(o.validate(), std::runtime_error);
  o = ClusterConfigOptions{};
  o.alpha_d_scale = 0.0;
  CHECK_THROWS_AS(o.validate(), std::runtime_error);
  // a grid that would drop more than 1e-3 of a table's norm
  o = ClusterConfigOptions{};
  o.rnp_max_fm = 2.0;
  CHECK_THROWS_AS(ClusterConfigSampler{o}, std::runtime_error);
  o = ClusterConfigOptions{};
  o.r_max_fm = 3.0;
  CHECK_THROWS_AS(ClusterConfigSampler{o}, std::runtime_error);
  // ... but the DEFAULTS must not throw, even though he4.density runs to
  // 20.05 > r_max_fm = 20 and fdeut.av18 to 100 > rnp_max_fm = 25.
  CHECK_NOTHROW(ClusterConfigSampler{ClusterConfigOptions{}});
  // theta_s = 0 is legal and a_2 == 0 there is the CORRECT value, not an error
  o = ClusterConfigOptions{};
  o.theta_s = 0.0;
  const ClusterConfigSampler lon(o);
  CHECK(lon.delta_perp_analytic_fm2(1) == 0.0);
  CHECK(lon.eps_b0_equivalent() == 0.0);
  CHECK(lon.q_matter_analytic_fm2(1) < 0.0);   // the geometry is unchanged
}

// ----------------------------------------------------------------------- T20

TEST_CASE("T20 throughput") {
  const ClusterConfigSampler& s = default_sampler();
  const auto t0 = std::chrono::steady_clock::now();
  const ClusterConfigSet set = s.sample_set(100000, 1, 1, 0);
  const auto t1 = std::chrono::steady_clock::now();
  const double us = std::chrono::duration<double, std::micro>(t1 - t0).count()
                    / static_cast<double>(set.config.size());
  MESSAGE("T20 " << us << " us per configuration");
  CHECK(us < 5.0);
}

// ----------------------------------------------------------------------- T21

TEST_CASE("T21 R-hat is conditioned on m_S, not on the m_S-summed density") {
  const ClusterConfigSampler& s = default_sampler();
  // (a) the GRID tables: for M = +1 the alpha-d orientation is pure L = 2 in
  //     two of the three branches, so these are exact.
  CHECK_CLOSE_AT(s.p2_alpha_d_grid(1, -1.0), -2.0 / 7.0, 0.0, 1e-6);
  CHECK_CLOSE_AT(s.p2_alpha_d_grid(1, 0.0), 1.0 / 7.0, 0.0, 1e-3);
  // The m_S = 0 branch is |Theta_2^1|^2 ~ (1-c^2)c^2, which has only a SIMPLE
  // zero at c = +-1, so the cos-theta midpoint rule leaves an O(h_c^2)
  // residue; the m_S = -1 branch's (1-c^2)^2 has a double zero and none.
  // Refining n_c drives it away as h_c^2, which is what makes it quadrature.
  ClusterConfigOptions fine_opt;
  fine_opt.n_c = 3072;
  const ClusterConfigSampler fine(fine_opt);
  CHECK_CLOSE_AT(fine.p2_alpha_d_grid(1, 0.0), 1.0 / 7.0, 0.0, 1e-6);
  CHECK_CLOSE_AT(fine.p2_alpha_d_grid(1, -1.0), -2.0 / 7.0, 0.0, 1e-6);
  // (b) the DRAWS.  An m_S-marginal draw gives <P2> = -0.0355 in BOTH
  //     branches; the conditioned one must separate them.
  const ClusterConfigSet set = s.sample_set(200000, 1, 99, 0);
  double n[3] = {0.0, 0.0, 0.0}, p2[3] = {0.0, 0.0, 0.0}, v2[3] = {0, 0, 0};
  for (const ClusterConfig& c : set.config) {
    // R = (r_p + r_n) * 3/4, and the axis is transverse (spin z -> lab x)
    const double rx = (c.nucleon[4].x + c.nucleon[5].x) * 0.75;
    const double ry = (c.nucleon[4].y + c.nucleon[5].y) * 0.75;
    const double rz = (c.nucleon[4].z + c.nucleon[5].z) * 0.75;
    const double r = std::sqrt(rx * rx + ry * ry + rz * rz);
    const double ct = rx / r;
    const double p = 0.5 * (3.0 * ct * ct - 1.0);
    const int k = 1 - static_cast<int>(std::lround(c.m_s));
    n[k] += 1.0;
    p2[k] += p;
    v2[k] += p * p;
  }
  for (int k = 0; k < 3; ++k) {
    const double mean = p2[k] / n[k];
    const double sd = std::sqrt(std::fmax(v2[k] / n[k] - mean * mean, 0.0));
    const double gate = 5.0 * sd / std::sqrt(n[k]);
    MESSAGE("T21 m_S = " << (1 - k) << ": <P2(cos theta_R)> = " << mean
            << " +- " << gate / 5.0 << " (n = " << n[k] << ")");
    if (k == 1) CHECK_CLOSE_AT(mean, 1.0 / 7.0, 0.0, gate);
    if (k == 2) CHECK_CLOSE_AT(mean, -2.0 / 7.0, 0.0, gate);
  }
  // and the m_S marginal itself is the closed form, at the BINOMIAL 5 sigma
  // (a 5 % relative gate on the m_S = 0 count is only 1.7 sigma)
  for (int k = 0; k < 3; ++k) {
    const double p = s.p_ms(1, static_cast<double>(1 - k));
    CHECK_CLOSE_AT(n[k] / 200000.0, p, 0.0,
                   5.0 * std::sqrt(p * (1.0 - p) / 200000.0));
  }
}

// ----------------------------------------------------------------------- T22

TEST_CASE("T22 the asymptotic D/S ratio, Whittaker-divided") {
  const ClusterConfigSampler& s = default_sampler();
  const double eta = s.asymptotic_ds_ratio();
  MESSAGE("T22 eta = " << eta << " (FitRescaled), measured -0.025(12)");
  CHECK_CLOSE_AT(eta, -0.05, 0.0, 0.01);
  // The naive R_2/R_0 is NOT eta: W_2/W_0 is ~3 over this window, so the
  // table's ratio is ~-0.14 while eta is ~-0.048.
  const double naive = np_interp(7.0, s.alpha_d_grid(), s.alpha_d_wave(2))
                       / np_interp(7.0, s.alpha_d_grid(), s.alpha_d_wave(0));
  CHECK_CLOSE_AT(naive, -0.149, 0.0, 0.01);
  CHECK(std::fabs(naive / eta) > 2.5);
  // ~2x the measured value, not 5-15x: a real but MODERATE D-wave excess.
  CHECK(std::fabs(eta / -0.025) < 3.0);
  CHECK(std::fabs(eta / -0.025) > 1.2);
  ClusterConfigOptions o;
  o.alpha_d_source = AlphaDSource::OverlapRaw;
  CHECK_CLOSE_AT(ClusterConfigSampler(o).asymptotic_ds_ratio(), -0.054, 0.0,
                 0.01);
}

// ------------------------------------------------ the sec. 8 numbers, pinned

TEST_CASE("the design's sec. 8 table, all three alpha-d sources") {
  struct Row {
    AlphaDSource src;
    double p_d, r2_ad, q, dt, r2, qm, delta, eps, a2;
  };
  const Row rows[3] = {
      {AlphaDSource::FitRescaled, 0.02011, 16.963, -1.3204, 0.98190, 6.4447,
       -1.2309, -0.1026, -0.0506, 0.197},
      {AlphaDSource::OverlapRaw, 0.02011, 17.548, -1.4009, 0.98190, 6.5747,
       -1.3382, -0.1115, -0.0550, 0.215},
      {AlphaDSource::FitRaw, 0.02542, 16.956, -1.4895, 0.97712, 6.4432,
       -1.4590, -0.1216, -0.0600, 0.234}};
  for (const Row& r : rows) {
    ClusterConfigOptions o;
    o.alpha_d_source = r.src;
    const ClusterConfigSampler s(o);
    CHECK_CLOSE(s.p_d_alpha_d(), r.p_d, 2e-4);
    CHECK_CLOSE(s.r2_alpha_d_fm2(), r.r2_ad, 2e-4);
    CHECK_CLOSE(s.q_int_fm2() + s.q_dd_fm2(), r.q, 2e-4);
    CHECK_CLOSE(s.tensor_dilution(), r.dt, 2e-5);
    CHECK_CLOSE(s.r2_analytic_fm2(), r.r2, 2e-4);
    CHECK_CLOSE(s.q_matter_analytic_fm2(1), r.qm, 2e-4);
    CHECK_CLOSE(s.delta_perp_analytic_fm2(1), r.delta, 1e-3);
    CHECK_CLOSE(s.eps_b0_equivalent(), r.eps, 1e-3);
    CHECK_CLOSE(s.a2_from_geometry(0.3, 1), r.a2, 5e-3);
  }
  // the alpha and deuteron inputs, once
  const ClusterConfigSampler& s = default_sampler();
  CHECK_CLOSE(s.r2_alpha_fm2(), 2.0774, 1e-4);
  CHECK_CLOSE(s.r2_np_fm2(), 15.4817, 1e-5);
  CHECK_CLOSE(s.q_int_fm2(), -1.2570, 3e-4);
  CHECK_CLOSE(s.q_dd_fm2(), -0.0631, 2e-3);
  // the anisotropy budget: alpha-d 76.9 % of |Q|, of which 95.2 % is S-D
  const double ad = std::fabs((4.0 / 3.0) * (s.q_int_fm2() + s.q_dd_fm2()));
  const double deut = std::fabs(2.0 * kQdFm2 * s.tensor_dilution());
  CHECK_CLOSE(ad / (ad + deut), 0.769, 2e-3);
  CHECK_CLOSE(std::fabs(s.q_int_fm2()) / std::fabs(s.q_int_fm2() + s.q_dd_fm2()),
              0.952, 2e-3);
}
