#include "lipolgen/fsi.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "lipolgen/numerics.hpp"

namespace lipolgen {
namespace {

/// J_m(x) for m = 0, 1, 2.  `std::cyl_bessel_j` is a C++17 special function
/// and is exact to the last ulp against SciPy here; it is called only while
/// the Bessel TABLES are built (~1e6 calls per layer), never per event.
double bessel_j(int m, double x) {
  if (x <= 0.0) return m == 0 ? 1.0 : 0.0;
  return std::cyl_bessel_j(static_cast<double>(m), x);
}

/// Midpoint grid x_i = (i + 1/2) dx on (0, x_max].
std::vector<double> midpoints(double x_max, std::size_t n) {
  std::vector<double> v(n);
  const double dx = x_max / static_cast<double>(n);
  for (std::size_t i = 0; i < n; ++i) v[i] = (static_cast<double>(i) + 0.5) * dx;
  return v;
}

Channel channel_of(const ClusterChannel& base) {
  if (base.beam_A == 6) {
    return base.spectator == "alpha" ? Channel::TaggedLi6Alpha
                                     : Channel::TaggedLi6D;
  }
  if (base.beam_A == 7) {
    return base.spectator == "alpha" ? Channel::TaggedLi7Alpha
                                     : Channel::TaggedLi7T;
  }
  if (base.beam_A == 3) return Channel::TaggedHe3P;
  return base.spectator == "p" ? Channel::TaggedDeuteronP
                               : Channel::TaggedDeuteronN;
}

}  // namespace

// ---------------------------------------------------------------- geometry

double cluster_point_a2_fm2(int z, int a) {
  // A single nucleon IS the point scatterer of the Glauber profile: folding
  // it with its own size would double-count what B_XN already carries.
  if (a == 1) return 0.0;
  // Measured charge radii [fm].  The proton's is folded OUT so that what is
  // left is the POINT-NUCLEON distribution the Glauber product needs.
  const double r_p = 0.8409;   // CODATA / PRad-era proton charge radius
  double r_ch = 0.0;
  if (z == 2 && a == 4) r_ch = 1.6755;        // 4He
  else if (z == 1 && a == 2) r_ch = 2.1413;   // deuteron
  else if (z == 1 && a == 3) r_ch = 1.7591;   // triton
  else if (z == 2 && a == 3) r_ch = 1.9661;   // 3He
  else {
    throw std::runtime_error(
        "cluster_point_a2_fm2: no measured charge radius for (Z=" +
        std::to_string(z) + ", A=" + std::to_string(a) + ")");
  }
  const double r_pt2 = r_ch * r_ch - r_p * r_p;
  if (!(r_pt2 > 0.0)) {
    throw std::runtime_error("cluster_point_a2_fm2: non-positive point radius");
  }
  return r_pt2 / 3.0;
}

const char* pipeline_fsi_name(PipelineFsi f) {
  switch (f) {
    case PipelineFsi::Off: return "off";
    case PipelineFsi::GlauberCluster: return "glauber-cluster";
    case PipelineFsi::GlauberNucleon: return "glauber-nucleon";
  }
  return "off";
}

// ------------------------------------------------------------ construction

GlauberFsiWeight::GlauberFsiWeight(TaggedChannel channel, GlauberFsiOptions opt)
    : channel_(std::move(channel)), opt_(opt) {
  build();
}

GlauberFsiWeight::GlauberFsiWeight(const TaggedModel& model,
                                   GlauberFsiOptions opt)
    : channel_(model.channel()), opt_(opt) {
  // Take the normalization grid from the model, so P_L means bit-for-bit
  // what it means in `TaggedModel::radial_table`.
  opt_.k_norm_max = model.k().back();
  opt_.n_norm = model.nk();
  build();
}

void GlauberFsiWeight::build() {
  channel_.validate();
  if (!(opt_.k_max > 0.0) || opt_.n_kz < 2 || opt_.n_kt < 2 || opt_.n_b < 2 ||
      opt_.n_kt_int < 2 || opt_.n_norm < 2) {
    throw std::runtime_error("GlauberFsiWeight: degenerate grid");
  }
  if (!(opt_.sigma_xn_mb >= 0.0) || !(opt_.b_xn > 0.0)) {
    throw std::runtime_error("GlauberFsiWeight: sigma_XN < 0 or B_XN <= 0");
  }
  if (!(opt_.w_max > 0.0)) {
    throw std::runtime_error("GlauberFsiWeight: w_max must be positive");
  }

  const ClusterChannel& base = channel_.base;
  spec_z_ = base.spectator_Z;
  spec_a_ = base.spectator_A;
  a2_fm2_ = cluster_point_a2_fm2(spec_z_, spec_a_);
  event_channel_ = channel_of(base);

  // --- the (L, |m|) components of the UNPOLARIZED density -----------------
  //
  // n(k) = sum_L P_L psihat_L(k)^2 / 4pi = sum_{L,m} (P_L/(2L+1))
  //        |psihat_L(k) Theta_L^m(c)|^2, and |Theta_L^{-m}| = |Theta_L^{+m}|,
  // so only |m| is carried and m != 0 counts twice.  There is NO S-D
  // interference in either the IA or the FSI density: it cancels in the sum
  // over the ion projection M, which is what makes the weight spin
  // independent (see the header).
  const double kappa = base.kappa();
  const std::vector<double> kn = linspace(1e-4, opt_.k_norm_max, opt_.n_norm);
  std::vector<double> integ(kn.size());
  for (const Wave& wv : channel_.waves) {
    if (!(wv.prob > 0.0)) continue;
    const std::vector<double> psi = wv.radial(kn, kappa);
    for (std::size_t i = 0; i < kn.size(); ++i) {
      integ[i] = psi[i] * psi[i] * kn[i] * kn[i];
    }
    const double norm = std::sqrt(trapezoid(integ, kn));
    if (!(norm > 0.0)) {
      throw std::runtime_error("GlauberFsiWeight: a radial wave has zero norm");
    }
    for (int m = 0; m <= wv.l; ++m) {
      Component c;
      c.l = wv.l;
      c.m = m;
      c.wgt = wv.prob / (2.0 * wv.l + 1.0) * (m == 0 ? 1.0 : 2.0);
      c.scale = 1.0 / norm;
      comp_.push_back(c);
    }
  }
  if (comp_.empty()) {
    throw std::runtime_error("GlauberFsiWeight: the channel carries no wave");
  }

  // --- the sigma_XN ladder ------------------------------------------------
  sigma_ladder_.clear();
  if (!opt_.formation_ramp) {
    sigma_ladder_.push_back(opt_.sigma_xn_mb);
  } else {
    const double lo = std::min(opt_.ramp_sigma_lo_mb, opt_.ramp_sigma_hi_mb);
    const double hi = std::max(opt_.ramp_sigma_lo_mb, opt_.ramp_sigma_hi_mb);
    const std::size_t n = std::max<std::size_t>(2, opt_.n_sigma_grid);
    sigma_ladder_ = linspace(lo, hi, n);
  }

  dkz_ = opt_.k_max / static_cast<double>(opt_.n_kz - 1);
  dkt_ = opt_.k_max / static_cast<double>(opt_.n_kt - 1);

  layers_.resize(sigma_ladder_.size());
  std::size_t clipped = 0, cells = 0;
  for (std::size_t i = 0; i < sigma_ladder_.size(); ++i) {
    build_layer(sigma_ladder_[i], layers_[i]);
    for (double v : layers_[i].w) {
      ++cells;
      if (v >= opt_.w_max * (1.0 - 1e-12)) ++clipped;
    }
  }
  clipped_ = cells ? static_cast<double>(clipped) / static_cast<double>(cells)
                   : 0.0;
}

// ------------------------------------------------------------- the profile

std::complex<double> GlauberFsiWeight::gamma_profile(double b_fm,
                                                     double sigma_xn_mb) const {
  // Gamma_N convolved with the cluster's Gaussian point-nucleon density: two
  // Gaussians, so the widths add and the normalisation is preserved,
  //     int d^2b (Gamma_N conv T_a) = sigma_XN (1 - i eps) / 2 .
  const double b_xn_fm2 = opt_.b_xn * HBAR_C_GEV_FM * HBAR_C_GEV_FM;
  const double w2 = b_xn_fm2 + a2_fm2_;
  const double sig_fm2 = sigma_xn_mb * MB_TO_FM2;
  const std::complex<double> amp(1.0, -opt_.eps);
  const std::complex<double> g1 =
      (sig_fm2 * amp / (4.0 * kPi * w2)) * std::exp(-b_fm * b_fm / (2.0 * w2));
  if (opt_.variant == FsiVariant::GlauberNucleon) {
    // (b) single-scattering (optical) limit: no Glauber shadowing at all, so
    // sigma_Xa = A sigma_XN.  The bracket ABOVE (a); see fsi.hpp.
    return static_cast<double>(spec_a_) * g1;
  }
  // (a) coherent cluster: Gamma_a = 1 - (1 - Gamma_N conv T_a)^A.
  std::complex<double> keep(1.0, 0.0);
  const std::complex<double> one_minus = 1.0 - g1;
  for (int i = 0; i < spec_a_; ++i) keep *= one_minus;
  return 1.0 - keep;
}

std::complex<double> GlauberFsiWeight::gtilde(double q_gev,
                                              double sigma_xn_mb) const {
  // Gtil(q) = 2 pi int b db J0(q b) Gamma(b), b in fm, returned in GeV^-2.
  const double db = opt_.b_max / static_cast<double>(opt_.n_b);
  std::complex<double> sum(0.0, 0.0);
  for (std::size_t j = 0; j < opt_.n_b; ++j) {
    const double b = (static_cast<double>(j) + 0.5) * db;      // GeV^-1
    sum += b * bessel_j(0, q_gev * b) * gamma_profile(b * HBAR_C_GEV_FM,
                                                      sigma_xn_mb);
  }
  return 2.0 * kPi * sum * db;
}

double GlauberFsiWeight::sigma_cluster_mb(double sigma_xn_mb) const {
  return 2.0 * gtilde(0.0, sigma_xn_mb).real() * GEV2_TO_MB;
}

double GlauberFsiWeight::sigma_cluster_el_mb(double sigma_xn_mb) const {
  const double db = opt_.b_max / static_cast<double>(opt_.n_b);
  double sum = 0.0;
  for (std::size_t j = 0; j < opt_.n_b; ++j) {
    const double b = (static_cast<double>(j) + 0.5) * db;
    sum += b * std::norm(gamma_profile(b * HBAR_C_GEV_FM, sigma_xn_mb));
  }
  return 2.0 * kPi * sum * db * GEV2_TO_MB;
}

double GlauberFsiWeight::slope_cluster_gev2(double sigma_xn_mb) const {
  const double q = 0.05;
  const double g0 = std::abs(gtilde(0.0, sigma_xn_mb));
  const double gq = std::abs(gtilde(q, sigma_xn_mb));
  return -2.0 * std::log(gq / g0) / (q * q);
}

// -------------------------------------------------------------- the layer

void GlauberFsiWeight::build_layer(double sigma_mb, Layer& out) {
  const std::size_t nb = opt_.n_b, nq = opt_.n_kt_int;
  const std::size_t nkz = opt_.n_kz, nkt = opt_.n_kt;
  const double db = opt_.b_max / static_cast<double>(nb);
  const double dq = opt_.kt_int_max / static_cast<double>(nq);

  const std::vector<double> b = midpoints(opt_.b_max, nb);
  const std::vector<double> q = midpoints(opt_.kt_int_max, nq);

  // Gamma(b) times the b db measure, once.
  std::vector<std::complex<double>> gb(nb);
  for (std::size_t j = 0; j < nb; ++j) {
    gb[j] = b[j] * db * gamma_profile(b[j] * HBAR_C_GEV_FM, sigma_mb);
  }

  int m_max = 0;
  for (const Component& c : comp_) m_max = std::max(m_max, c.m);
  // Jf[m][j * nq + i] = J_m(q_i b_j)   -- inner loop over i is contiguous
  // Jb[m][o * nb + j] = J_m(kt_o b_j)  -- inner loop over j is contiguous
  std::vector<std::vector<double>> jf(m_max + 1), jb(m_max + 1);
  for (int m = 0; m <= m_max; ++m) {
    jf[m].resize(nb * nq);
    for (std::size_t j = 0; j < nb; ++j) {
      for (std::size_t i = 0; i < nq; ++i) {
        jf[m][j * nq + i] = bessel_j(m, q[i] * b[j]);
      }
    }
    jb[m].resize(nkt * nb);
    for (std::size_t o = 0; o < nkt; ++o) {
      const double kt = static_cast<double>(o) * dkt_;
      for (std::size_t j = 0; j < nb; ++j) {
        jb[m][o * nb + j] = bessel_j(m, kt * b[j]);
      }
    }
  }

  const double kappa = channel_.base.kappa();
  // The radial evaluator of each L, hoisted out of the k_z loop.
  std::vector<const Wave*> wave_of(comp_.size(), nullptr);
  for (std::size_t ci = 0; ci < comp_.size(); ++ci) {
    for (const Wave& wv : channel_.waves) {
      if (wv.l == comp_[ci].l && wv.prob > 0.0) { wave_of[ci] = &wv; break; }
    }
  }

  out.sigma_mb = sigma_mb;
  out.w.assign(nkz * nkt, 1.0);

  std::vector<double> g(nq), gfor(nb), num(nkt), den(nkt);
  std::vector<std::complex<double>> hb(nb);
  for (std::size_t p = 0; p < nkz; ++p) {
    const double kz = static_cast<double>(p) * dkz_;
    std::fill(num.begin(), num.end(), 0.0);
    std::fill(den.begin(), den.end(), 0.0);

    for (std::size_t ci = 0; ci < comp_.size(); ++ci) {
      const Component& c = comp_[ci];
      const Wave& wv = *wave_of[ci];

      // g_{Lm}(k_z, k'_T) on the integration grid, times k'_T dk'_T.
      for (std::size_t i = 0; i < nq; ++i) {
        const double kk = std::sqrt(kz * kz + q[i] * q[i]);
        const double cs = kk > 1e-12 ? kz / kk : 1.0;
        g[i] = c.scale * wv.radial(kk, kappa) * theta_lm(c.l, c.m, cs) *
               q[i] * dq;
      }
      // G_m(b) = int k'_T dk'_T J_m(k'_T b) g(k'_T)   -- forward Hankel.
      const double* jfm = jf[c.m].data();
      for (std::size_t j = 0; j < nb; ++j) {
        const double* row = jfm + j * nq;
        double s = 0.0;
        for (std::size_t i = 0; i < nq; ++i) s += row[i] * g[i];
        gfor[j] = s;
      }
      for (std::size_t j = 0; j < nb; ++j) hb[j] = gb[j] * gfor[j];

      // corr(k_T) = int b db J_m(k_T b) Gamma(b) G_m(b)  -- backward Hankel.
      const double* jbm = jb[c.m].data();
      for (std::size_t o = 0; o < nkt; ++o) {
        const double* row = jbm + o * nb;
        double cre = 0.0, cim = 0.0;
        for (std::size_t j = 0; j < nb; ++j) {
          cre += row[j] * hb[j].real();
          cim += row[j] * hb[j].imag();
        }
        const double kt = static_cast<double>(o) * dkt_;
        const double kk = std::sqrt(kz * kz + kt * kt);
        const double cs = kk > 1e-12 ? kz / kk : 1.0;
        const double g0 = c.scale * wv.radial(kk, kappa) *
                          theta_lm(c.l, c.m, cs);
        // |g0 - corr|^2 with the QUADRATIC (gain) term rescaled by
        // `elastic_gain`: the linear term carries sigma_tot, the quadratic
        // one sigma_el (fsi.hpp).
        num[o] += c.wgt * (g0 * g0 - 2.0 * g0 * cre +
                           opt_.elastic_gain * (cre * cre + cim * cim));
        den[o] += c.wgt * g0 * g0;
      }
    }

    for (std::size_t o = 0; o < nkt; ++o) {
      double w = 1.0;
      if (den[o] > 0.0) w = num[o] / den[o];
      if (!(w >= 0.0)) w = 0.0;
      if (w > opt_.w_max) w = opt_.w_max;
      out.w[p * nkt + o] = w;
    }
  }

  // --- the rate-weighted mean = the tagged-cluster survival probability ----
  //
  //   <w> = int w(k, c) n(k) k^2 dk dc / int n(k) k^2 dk dc,
  // with n(k) = sum_L P_L psihat_L(k)^2 / 4pi the ISOTROPIC unpolarized
  // density, on the model's own (k, c) grid.  Below 1 because the
  // rescattering is nearly purely absorptive: LOG IT, it is physical.
  const std::size_t nk_q = 240, nc_q = 96;
  double sw = 0.0, sn = 0.0;
  for (std::size_t ik = 0; ik < nk_q; ++ik) {
    const double kk = (static_cast<double>(ik) + 0.5) * opt_.k_max /
                      static_cast<double>(nk_q);
    double nk = 0.0;
    for (std::size_t ci = 0; ci < comp_.size(); ++ci) {
      if (comp_[ci].m != 0) continue;   // one entry per L
      const double psi = comp_[ci].scale * wave_of[ci]->radial(kk, kappa);
      double prob = 0.0;
      for (const Wave& wv : channel_.waves) {
        if (wv.l == comp_[ci].l) { prob = wv.prob; break; }
      }
      nk += prob * psi * psi;
    }
    const double meas = nk * kk * kk;
    for (std::size_t ic = 0; ic < nc_q; ++ic) {
      const double cs = -1.0 + (static_cast<double>(ic) + 0.5) * 2.0 /
                                    static_cast<double>(nc_q);
      const double kz = std::fabs(kk * cs);
      const double kt = kk * std::sqrt(std::fmax(0.0, 1.0 - cs * cs));
      sw += meas * read(out, kz, kt);
      sn += meas;
    }
  }
  out.mean_w = sn > 0.0 ? sw / sn : 1.0;
}

// ------------------------------------------------------------- the lookup

double GlauberFsiWeight::read(const Layer& lay, double kz, double kt) const {
  const std::size_t nkt = opt_.n_kt;
  double fz = kz / dkz_, ft = kt / dkt_;
  if (fz < 0.0) fz = 0.0;
  if (ft < 0.0) ft = 0.0;
  const double zmax = static_cast<double>(opt_.n_kz - 1);
  const double tmax = static_cast<double>(nkt - 1);
  if (fz > zmax) fz = zmax;
  if (ft > tmax) ft = tmax;
  std::size_t iz = static_cast<std::size_t>(fz);
  std::size_t it = static_cast<std::size_t>(ft);
  if (iz + 1 >= opt_.n_kz) iz = opt_.n_kz - 2;
  if (it + 1 >= nkt) it = nkt - 2;
  const double uz = fz - static_cast<double>(iz);
  const double ut = ft - static_cast<double>(it);
  const double w00 = lay.w[iz * nkt + it];
  const double w01 = lay.w[iz * nkt + it + 1];
  const double w10 = lay.w[(iz + 1) * nkt + it];
  const double w11 = lay.w[(iz + 1) * nkt + it + 1];
  return (1.0 - uz) * ((1.0 - ut) * w00 + ut * w01) +
         uz * ((1.0 - ut) * w10 + ut * w11);
}

void GlauberFsiWeight::ladder_at(double w, std::size_t& i0, std::size_t& i1,
                                 double& f) const {
  i0 = i1 = 0;
  f = 0.0;
  if (sigma_ladder_.size() < 2) return;
  const double s = sigma_eff_mb(w);
  // The ladder is increasing in sigma.
  const double lo = sigma_ladder_.front(), hi = sigma_ladder_.back();
  const double t = (s - lo) / (hi - lo) *
                   static_cast<double>(sigma_ladder_.size() - 1);
  double tc = clip(t, 0.0, static_cast<double>(sigma_ladder_.size() - 1));
  i0 = static_cast<std::size_t>(tc);
  if (i0 + 1 >= sigma_ladder_.size()) i0 = sigma_ladder_.size() - 2;
  i1 = i0 + 1;
  f = tc - static_cast<double>(i0);
}

double GlauberFsiWeight::sigma_eff_mb(double w) const {
  if (!opt_.formation_ramp) return opt_.sigma_xn_mb;
  const double w_lo = opt_.ramp_w_lo, w_hi = opt_.ramp_w_hi;
  if (!(w_hi > w_lo)) return opt_.ramp_sigma_hi_mb;
  const double u = clip((w - w_lo) / (w_hi - w_lo), 0.0, 1.0);
  return opt_.ramp_sigma_lo_mb +
         u * (opt_.ramp_sigma_hi_mb - opt_.ramp_sigma_lo_mb);
}

double GlauberFsiWeight::weight(const FsiKinematics& kin) const {
  // A weight built for one channel must not silently distort another.
  if (kin.spectator_a > 0 &&
      (kin.spectator_a != spec_a_ || kin.spectator_z != spec_z_)) {
    return 1.0;
  }
  if (kin.channel != Channel::Inclusive && kin.channel != event_channel_) {
    return 1.0;
  }
  double k = kin.k;
  if (!(k >= 0.0)) return 1.0;
  if (k > opt_.k_max) k = opt_.k_max;
  const double c = clip(kin.cos_theta_k, -1.0, 1.0);
  const double kz = std::fabs(k * c);
  const double kt = k * std::sqrt(std::fmax(0.0, 1.0 - c * c));
  if (layers_.size() < 2) return read(layers_[0], kz, kt);
  std::size_t i0 = 0, i1 = 0;
  double f = 0.0;
  ladder_at(kin.w, i0, i1, f);
  return (1.0 - f) * read(layers_[i0], kz, kt) + f * read(layers_[i1], kz, kt);
}

double GlauberFsiWeight::weight_normalised(const FsiKinematics& kin) const {
  const double w = weight(kin);
  double mean = layers_[0].mean_w;
  if (layers_.size() >= 2) {
    std::size_t i0 = 0, i1 = 0;
    double f = 0.0;
    ladder_at(kin.w, i0, i1, f);
    mean = (1.0 - f) * layers_[i0].mean_w + f * layers_[i1].mean_w;
  }
  return mean > 0.0 ? w / mean : w;
}

double GlauberFsiWeight::survival(double w) const {
  if (layers_.size() < 2) return layers_[0].mean_w;
  std::size_t i0 = 0, i1 = 0;
  double f = 0.0;
  ladder_at(w, i0, i1, f);
  return (1.0 - f) * layers_[i0].mean_w + f * layers_[i1].mean_w;
}

const std::vector<double>& GlauberFsiWeight::grid(std::size_t i) const {
  if (i >= layers_.size()) {
    throw std::runtime_error("GlauberFsiWeight::grid: index past the ladder");
  }
  return layers_[i].w;
}

}  // namespace lipolgen
