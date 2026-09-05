// SPDX-License-Identifier: GPL-3.0-or-later
// The triton's nucleon spectral function.  Every coefficient, its provenance
// and every physics choice is in triton_sf.hpp.

#include "lipolgen/triton_sf.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "lipolgen/numerics.hpp"

namespace lipolgen {
namespace {

// --- Ciofi degli Atti - Simula, PRC 53 (1996) 1689 -------------------------
// n_0: Table A.1 (A = 2, 3, 4).  A = 3 is their "3He (proton)" column.
const std::vector<CsTerm> kN0A2 = {{157.4, 1.24, 18.3},
                                   {0.234, 1.27, 0.0},
                                   {0.00623, 0.220, 0.0}};
const std::vector<CsTerm> kN0A3 = {{31.7, 1.32, 5.98}, {0.00266, 0.365, 0.0}};
const std::vector<CsTerm> kN0A4 = {{4.33, 1.54, 0.419}, {5.49, 4.90, 0.0}};

// n_1: Eq. (76) for A = 3 -- the one nucleus CS parameterise separately --
// and Eq. (75) with Table A.3 for A = 4.  The deuteron has none.
const std::vector<CsTerm> kN1A2 = {};
const std::vector<CsTerm> kN1A3 = {{7.40, 1.23, 3.21}, {0.0139, 0.234, 0.0}};
const std::vector<CsTerm> kN1A4 = {{0.665, 2.15, 0.0}, {0.0244, 0.22, 0.0}};

/// Composite Simpson of `f` on [0, b] with `n` (even) intervals.  The
/// integrands here are smooth and this is a setup/diagnostic path only --
/// nothing per event calls it.
template <class F>
double simpson(F f, double b, std::size_t n) {
  const double h = b / static_cast<double>(n);
  double s = f(0.0) + f(b);
  for (std::size_t i = 1; i < n; ++i) {
    s += ((i & 1U) ? 4.0 : 2.0) * f(static_cast<double>(i) * h);
  }
  return s * h / 3.0;
}

constexpr std::size_t kQuadN = 100000;

}  // namespace

const std::vector<CsTerm>& cs_n0_terms(int a) {
  switch (a) {
    case 2: return kN0A2;
    case 3: return kN0A3;
    case 4: return kN0A4;
    default: break;
  }
  throw std::runtime_error(
      "cs_n0_terms: Ciofi-Simula gives no n_0 set for A = " + std::to_string(a));
}

const std::vector<CsTerm>& cs_n1_terms(int a) {
  switch (a) {
    case 2: return kN1A2;
    case 3: return kN1A3;
    case 4: return kN1A4;
    default: break;
  }
  throw std::runtime_error(
      "cs_n1_terms: Ciofi-Simula gives no n_1 set for A = " + std::to_string(a));
}

double cs_sum(const std::vector<CsTerm>& terms, double k_fm) {
  const double k2 = k_fm * k_fm;
  double tot = 0.0;
  for (const CsTerm& t : terms) {
    const double d = 1.0 + t.c * k2;
    tot += t.a * std::exp(-t.b * k2) / (d * d);
  }
  return tot;
}

double cs_norm(const std::vector<CsTerm>& terms, double k_max_fm) {
  if (terms.empty()) return 0.0;
  return simpson([&](double k) { return k * k * cs_sum(terms, k); }, k_max_fm,
                 kQuadN);
}

double cs_p_above(const std::vector<CsTerm>& terms, double k_gev,
                  double k_max_fm) {
  if (terms.empty()) return 0.0;
  const double den = cs_norm(terms, k_max_fm);
  if (!(den > 0.0)) return 0.0;
  const double k0 = k_gev / HBARC_GEV_FM;
  if (k0 >= k_max_fm) return 0.0;
  // int_{k0}^{kmax} = int_0^{kmax - k0} shifted, so Simpson stays on [0, b].
  const double num = simpson(
      [&](double u) {
        const double k = k0 + u;
        return k * k * cs_sum(terms, k);
      },
      k_max_fm - k0, kQuadN);
  return num / den;
}

double cs_mean_k(const std::vector<CsTerm>& terms, double k_max_fm) {
  if (terms.empty()) return 0.0;
  const double den = cs_norm(terms, k_max_fm);
  if (!(den > 0.0)) return 0.0;
  const double num =
      simpson([&](double k) { return k * k * k * cs_sum(terms, k); }, k_max_fm,
              kQuadN);
  return HBARC_GEV_FM * num / den;
}

const char* triton_channel_name(TritonChannel c) {
  switch (c) {
    case TritonChannel::NeutronD: return "n+d";
    case TritonChannel::NeutronPnCont: return "n+(pn)";
    case TritonChannel::ProtonNnCont: return "p+(nn)";
  }
  return "?";
}

// ---------------------------------------------------------------------------

CiofiSimulaTriton::CiofiSimulaTriton(CiofiSimulaOptions opt)
    : opt_(std::move(opt)) {
  if (!(opt_.k_max > 0.0)) throw std::runtime_error("CiofiSimulaTriton: k_max <= 0");
  if (!(opt_.q_max > 0.0)) throw std::runtime_error("CiofiSimulaTriton: q_max <= 0");
  if (!(opt_.beta > 0.0)) throw std::runtime_error("CiofiSimulaTriton: beta <= 0");
  if (opt_.n_grid < 16) throw std::runtime_error("CiofiSimulaTriton: n_grid too small");
  if (!(opt_.n1_scale >= 0.0))
    throw std::runtime_error("CiofiSimulaTriton: n1_scale < 0");

  m_p_ = nuclear_mass(1, 1);
  m_n_ = nuclear_mass(0, 1);
  m_d_ = nuclear_mass(1, 2);

  // S_0 and S_1 as this object integrates them -- over the SAMPLING ceiling,
  // not to infinity, so `p_two_body`'s k-integral and the channel fractions a
  // run actually produces are the same number.
  const double k_max_fm = opt_.k_max / HBARC_GEV_FM;
  s0_ = cs_norm(cs_n0_terms(3), k_max_fm);
  s1_ = opt_.n1_scale * cs_norm(cs_n1_terms(3), k_max_fm);

  kgrid_ = linspace(1e-4, opt_.k_max, opt_.n_grid);
  qgrid_ = linspace(1e-4, opt_.q_max, opt_.n_grid);

  // |k| CDFs.  k^2 n(k) on the grid, cumulated left to right and normalised --
  // the same construction `ClusterBreakup::build_cdf` uses, so the two agree
  // wherever they share a density.
  cdf_k_n_.assign(kgrid_.size(), 0.0);
  cdf_k_p_.assign(kgrid_.size(), 0.0);
  double acc_n = 0.0, acc_p = 0.0;
  for (std::size_t i = 0; i < kgrid_.size(); ++i) {
    const double k = kgrid_[i];
    const double n0 = n0_cs(k), n1 = n1_cs(k);
    acc_n += k * k * (n0 + n1);
    acc_p += k * k * (opt_.proton_n1_only ? n1 : (n0 + n1));
    cdf_k_n_[i] = acc_n;
    cdf_k_p_[i] = acc_p;
  }
  for (double& v : cdf_k_n_) v /= acc_n;
  for (double& v : cdf_k_p_) v /= acc_p;

  cdf_q_nn_ = pair_cdf(opt_.kappa_nn, &qnorm_nn_);
  cdf_q_pn_ = pair_cdf(opt_.kappa_pn, &qnorm_pn_);
}

std::vector<double> CiofiSimulaTriton::pair_cdf(double kappa,
                                                double* norm) const {
  std::vector<double> cdf(qgrid_.size(), 0.0);
  double acc = 0.0;
  for (std::size_t i = 0; i < qgrid_.size(); ++i) {
    const double q = qgrid_[i];
    // q^2 |psi_{L=0}(q; kappa, beta)|^2 -- `momentum_density` is |psi|^2
    // without the phase space, exactly as `MomentumSampler` builds its own.
    acc += q * q * momentum_density(q, kappa, opt_.beta, 0);
    cdf[i] = acc;
  }
  // The same sum as an integral, for the continuous `e_rel_density`.
  *norm = acc * (qgrid_[1] - qgrid_[0]);
  if (acc > 0.0) {
    for (double& v : cdf) v /= acc;
  }
  return cdf;
}

double CiofiSimulaTriton::draw_from(const std::vector<double>& cdf,
                                    const std::vector<double>& grid,
                                    double u) const {
  const auto it = std::lower_bound(cdf.begin(), cdf.end(), u);
  if (it == cdf.begin()) return grid.front();
  if (it == cdf.end()) return grid.back();
  const std::size_t i = static_cast<std::size_t>(it - cdf.begin());
  const double c0 = cdf[i - 1], c1 = cdf[i];
  const double t = (c1 > c0) ? (u - c0) / (c1 - c0) : 0.0;
  return grid[i - 1] + t * (grid[i] - grid[i - 1]);
}

double CiofiSimulaTriton::n0_cs(double k_gev) const {
  return cs_sum(cs_n0_terms(3), k_gev / HBARC_GEV_FM);
}

double CiofiSimulaTriton::n1_cs(double k_gev) const {
  return opt_.n1_scale * cs_sum(cs_n1_terms(3), k_gev / HBARC_GEV_FM);
}

double CiofiSimulaTriton::n_of_k(double k, int pdg) const {
  if (!(k >= 0.0)) return 0.0;
  const double n0 = n0_cs(k), n1 = n1_cs(k);
  const bool proton = (pdg == 2212);
  const double num = (proton && opt_.proton_n1_only) ? n1 : (n0 + n1);
  const double den = (proton && opt_.proton_n1_only) ? s1_ : (s0_ + s1_);
  if (!(den > 0.0)) return 0.0;
  // CS normalise int dk k^2 n = S (no 4 pi), so the density over d^3k is
  // n/(4 pi S); the fm -> GeV conversion is one factor of (hbar c)^3 because
  // n has units of volume.
  const double hc3 = HBARC_GEV_FM * HBARC_GEV_FM * HBARC_GEV_FM;
  return num / (4.0 * kPi * den * hc3);
}

double CiofiSimulaTriton::p_two_body(double k, int pdg) const {
  if (pdg == 2212) return 0.0;          // there is no bound nn
  const double n0 = n0_cs(k), n1 = n1_cs(k);
  const double tot = n0 + n1;
  return (tot > 0.0) ? n0 / tot : 0.0;
}

double CiofiSimulaTriton::e_max(int pdg) const {
  const double q = opt_.q_max;
  const double m1 = (pdg == 2212) ? m_n_ : m_p_;   // the pair the SPECIES leaves
  const double m2 = m_n_;
  return std::sqrt(m1 * m1 + q * q) + std::sqrt(m2 * m2 + q * q) - m1 - m2;
}

double CiofiSimulaTriton::e_rel_density(double k, double e_rel, int pdg) const {
  (void)k;   // the model factorises: the pair's own scale does not depend on k
  if (!(e_rel > 0.0) || e_rel > e_max(pdg)) return 0.0;
  const double m1 = (pdg == 2212) ? m_n_ : m_p_;
  const double m2 = m_n_;
  // Invert e(q) = sqrt(m1^2+q^2) + sqrt(m2^2+q^2) - m1 - m2 by bisection: it
  // is monotone in q, and this is a diagnostic path, not the event loop.
  double lo = 0.0, hi = opt_.q_max;
  auto e_of = [&](double q) {
    return std::sqrt(m1 * m1 + q * q) + std::sqrt(m2 * m2 + q * q) - m1 - m2;
  };
  for (int it = 0; it < 200; ++it) {
    const double mid = 0.5 * (lo + hi);
    (e_of(mid) < e_rel ? lo : hi) = mid;
  }
  const double q = 0.5 * (lo + hi);
  if (!(q > 0.0)) return 0.0;
  const double dedq = q / std::sqrt(m1 * m1 + q * q) + q / std::sqrt(m2 * m2 + q * q);
  if (!(dedq > 0.0)) return 0.0;
  const double kappa = (pdg == 2212) ? opt_.kappa_nn : opt_.kappa_pn;
  const double norm = (pdg == 2212) ? qnorm_nn_ : qnorm_pn_;
  if (!(norm > 0.0)) return 0.0;
  // p(q) = q^2 |psi|^2 / int q'^2 |psi|^2 dq', and p(e) = p(q) / (de/dq).
  const double pq = q * q * momentum_density(q, kappa, opt_.beta, 0) / norm;
  return pq / dedq;
}

TritonDraw CiofiSimulaTriton::sample(double p_proton, Rng& rng) const {
  // ---- the fixed six uniforms.  Drawn UP FRONT and unconditionally, so the
  // stream position after this call cannot depend on which channel came out.
  const double u_species = rng.uniform();
  const double u_k = rng.uniform();
  const double u_branch = rng.uniform();
  const double u_cos = rng.uniform();
  const double u_phi = rng.uniform();
  const double u_q = rng.uniform();

  TritonDraw d;
  d.struck_proton = u_species < p_proton;
  const int pdg = d.struck_proton ? 2212 : 2112;
  d.k = draw_from(d.struck_proton ? cdf_k_p_ : cdf_k_n_, kgrid_, u_k);
  d.cos_theta_k = 2.0 * u_cos - 1.0;
  d.phi_k = 2.0 * kPi * u_phi;

  const bool bound = !d.struck_proton && (u_branch < p_two_body(d.k, pdg));
  if (bound) {
    d.channel = TritonChannel::NeutronD;
    d.m_remnant = m_d_;
    d.e_rel = 0.0;
    d.q_pair = 0.0;         // the bound remnant has no internal continuum
    return d;
  }
  d.channel = d.struck_proton ? TritonChannel::ProtonNnCont
                              : TritonChannel::NeutronPnCont;
  const double m1 = d.struck_proton ? m_n_ : m_p_;
  const double m2 = m_n_;
  d.q_pair = draw_from(d.struck_proton ? cdf_q_nn_ : cdf_q_pn_, qgrid_, u_q);
  const double q2 = d.q_pair * d.q_pair;
  d.m_remnant = std::sqrt(m1 * m1 + q2) + std::sqrt(m2 * m2 + q2);
  d.e_rel = d.m_remnant - m1 - m2;
  return d;
}

}  // namespace lipolgen
