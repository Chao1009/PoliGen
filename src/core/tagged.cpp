#include "lipolgen/tagged.hpp"

#include <algorithm>
#include <memory>
#include <cmath>
#include <stdexcept>

#include "lipolgen/constants.hpp"
#include "lipolgen/fsi.hpp"
#include "lipolgen/numerics.hpp"
#include "lipolgen/sampler.hpp"   // w2_from_xq2: the ONE W2(x, Q2) definition
#include "lipolgen/spin.hpp"

namespace lipolgen {
namespace {

/// Map a half-integer projection onto an exact integer map key.
long long mkey(double m) { return std::llround(2.0 * m); }

/// One key for the (M, m_S) pair, exact for half-integers.
long long cdf_key(double m_ion, double m_s) {
  return mkey(m_ion) * 1000 + mkey(m_s) + 500;
}

const double kNaN = std::nan("");

/// `np.clip(np.searchsorted(grid, v) - 1, 0, n - 2)`: the module's
/// nearest-cell-at-or-below lookup, which is part of the model rather than an
/// implementation detail (the reference values are grid CELL values).
std::size_t cell_index(const std::vector<double>& grid, double v) {
  const std::size_t n = grid.size();
  const std::size_t ins = static_cast<std::size_t>(
      std::lower_bound(grid.begin(), grid.end(), v) - grid.begin());
  if (ins == 0) return 0;
  const std::size_t i = ins - 1;
  return i > n - 2 ? n - 2 : i;
}

/// The `Channel` label of a tagged `ClusterChannel` -- used by `fill_event`
/// (the record label) and `apply_fsi` (the FSI channel guard), so the two
/// cannot disagree.
Channel channel_enum_of(const ClusterChannel& base) {
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

const Ion& TRITON() {
  static const Ion i{"t", 3, 1, 0.5, 0.86, -0.028};
  return i;
}

const Ion& NEUTRON_TARGET() {
  static const Ion i{"n", 1, 0, 0.5, 0.0, 1.0};
  return i;
}

void TaggedChannel::validate() const {
  double tot = 0.0;
  for (const Wave& w : waves) tot += w.prob;
  if (std::fabs(tot - 1.0) > 1e-9) {
    throw std::runtime_error("wave probabilities must sum to 1");
  }
}

namespace {

/// The 6Li alpha-d VMC waves, built ONCE and shared.
///
/// MAGNITUDE from `momenta/li6_ad1.momentum` (AV18+UX, 1M samples, an
/// explicit S/D split with MC errors and the file's own normalizations).
/// SIGN from `overlap_old/li6.ad` (AV18+UIX, 2004), which is the only
/// published source of the RELATIVE S-D phase: a momentum density is
/// |psi_L|^2 and carries none.  `vmc_from_momentum` turns the overlap
/// column's zero crossings below 3 fm^-1 into a step function of the sign,
/// which is more robust than copying sign(A_L(k)) point by point (past
/// ~3 fm^-1 the overlap columns are at the MC noise floor and wander, while
/// carrying ~1e-4 of the norm).  The resulting node structure is
///     S wave  one node at 0.678 fm^-1 = 0.134 GeV
///     D wave  one node at 2.25  fm^-1 = 0.444 GeV
/// and both are confirmed independently by minima of the momentum file's own
/// rho_0 / rho_2 in the same bin (docs/open_items/vmc_reconciliation.md).
/// With the global phase fixed by psi_0(k -> 0) > 0, sign(psi_2/psi_0) = -1
/// below the S node -- the OPPOSITE of what the positive-definite Hulthen
/// forms assume.
const std::vector<Wave>& li6_vmc_waves() {
  static const std::vector<Wave> w = [] {
    const VmcRadial s_sign = vmc_from_overlap_k(data_path(VMC_LI6_OVERLAP), 0, 0);
    const VmcRadial d_sign = vmc_from_overlap_k(data_path(VMC_LI6_OVERLAP), 1, 2);
    VmcRadial s_tab = vmc_from_momentum(data_path(VMC_LI6_MOMENTUM), 1, 0, 0,
                                        &s_sign);
    VmcRadial d_tab = vmc_from_momentum(data_path(VMC_LI6_MOMENTUM), 1, 1, 2,
                                        &d_sign);
    // Fix the (unobservable) GLOBAL phase to psi_0(k -> 0) > 0, so that the
    // (observable) relative phase reads off as sign(psi_2).  li6.ad happens
    // to print the alpha-d overlap with A00 < 0 at low k.
    if (s_tab.psi().front() < 0.0) {
      s_tab = s_tab.scaled(-1.0);
      d_tab = d_tab.scaled(-1.0);
    }
    auto s = std::make_shared<const VmcRadial>(std::move(s_tab));
    auto d = std::make_shared<const VmcRadial>(std::move(d_tab));
    Wave ws;
    ws.l = 0;
    ws.prob = 1.0 - VMC_P_D_LI6;
    ws.vmc = s;
    Wave wd;
    wd.l = 2;
    wd.prob = VMC_P_D_LI6;
    wd.vmc = d;
    return std::vector<Wave>{ws, wd};
  }();
  return w;
}

/// The 7Li alpha-t VMC wave.  A single L = 1 channel, so there is no
/// interference term and no observable phase: the magnitude from
/// `momenta/li7_at3.momentum` (the 3/2- GROUND state; `li7_at1` is the 1/2-
/// excited state) is the whole story and no sign reference is needed.
///
/// `overlap_old/li7.at`'s second column `Aat11` is NOT a second partial wave:
/// alpha(0+) x t(1/2+) with L = 1 gives j = 1/2 or 3/2 and only j = 3/2 can
/// build the 3/2- ground state, so it is a selection-rule zero carrying MC
/// leakage at 3.5e-5 of the norm.  It is deliberately not loaded.
const std::vector<Wave>& li7_vmc_waves() {
  static const std::vector<Wave> w = [] {
    auto p = std::make_shared<const VmcRadial>(vmc_from_momentum(
        data_path(VMC_LI7_MOMENTUM), 0, 0, 1, nullptr));
    Wave wp;
    wp.l = 1;
    wp.prob = 1.0;
    wp.vmc = p;
    return std::vector<Wave>{wp};
  }();
  return w;
}

}  // namespace

TaggedChannel li6_alpha_channel(double beta, double p_d,
                                ClusterWaveSource source) {
  TaggedChannel c;
  c.base = LI6_ALPHA_TAG();
  c.j_ion = 1.0;
  c.s_struck = 1.0;
  c.s_spec = 0.0;
  c.s_channel = 1.0;
  if (source == ClusterWaveSource::VmcAV18) {
    c.waves = li6_vmc_waves();
    c.label = "6Li alpha-tag (embedded d, VMC AV18+UX)";
  } else {
    c.waves = {Wave{0, 1.0 - p_d, beta, nullptr}, Wave{2, p_d, beta, nullptr}};
    c.label = "6Li alpha-tag (embedded d)";
  }
  c.dis_target = DEUTERON();
  c.validate();
  return c;
}

TaggedChannel li7_alpha_channel(double beta, ClusterWaveSource source) {
  TaggedChannel c;
  c.base = LI7_ALPHA_TAG();
  c.j_ion = 1.5;
  c.s_struck = 0.5;
  c.s_spec = 0.0;
  c.s_channel = 0.5;
  if (source == ClusterWaveSource::VmcAV18) {
    c.waves = li7_vmc_waves();
    c.label = "7Li alpha-tag (quasi-free t, VMC AV18+UX)";
  } else {
    c.waves = {Wave{1, 1.0, beta, nullptr}};
    c.label = "7Li alpha-tag (quasi-free t)";
  }
  c.dis_target = TRITON();
  c.validate();
  return c;
}

TaggedChannel deuteron_channel(double beta, double p_d) {
  TaggedChannel c;
  c.base = DEUTERON_P_TAG();
  c.j_ion = 1.0;
  c.s_struck = 0.5;
  c.s_spec = 0.5;
  c.s_channel = 1.0;
  c.waves = {Wave{0, 1.0 - p_d, beta, nullptr}, Wave{2, p_d, beta, nullptr}};
  c.dis_target = NEUTRON_TARGET();
  c.label = "d control (n struck, p tagged)";
  c.validate();
  return c;
}

// ------------------------------------------------------------- TaggedModel

TaggedModel::TaggedModel(TaggedChannel channel, double k_max, std::size_t nk,
                         std::size_t nc)
    : channel_(std::move(channel)) {
  channel_.validate();
  k_ = linspace(1e-4, k_max, nk);
  c_ = linspace(-1.0 + 1.0 / static_cast<double>(nc),
                1.0 - 1.0 / static_cast<double>(nc), nc);
  dk_ = k_[1] - k_[0];
  dc_ = 2.0 / static_cast<double>(nc);
  const double kappa = channel_.base.kappa();
  std::vector<double> integrand(nk);
  for (const Wave& w : channel_.waves) {
    std::vector<double> psi = w.radial(k_, kappa);
    for (std::size_t i = 0; i < nk; ++i) integrand[i] = psi[i] * psi[i] * k_[i] * k_[i];
    const double norm = std::sqrt(trapezoid(integrand, k_));
    const double a = std::sqrt(w.prob) / norm;
    for (double& v : psi) v *= a;
    rad_[w.l] = std::move(psi);
  }
  ms_struck_ = m_values(channel_.s_channel);

  // C3: build every grid HERE.  After this the object is immutable, so every
  // accessor is a pure lookup and `TaggedModel` is safe to share between
  // threads with no lock and no warm-up protocol from the caller.
  for (double mi : m_values(channel_.j_ion)) {
    std::vector<std::vector<double>> a2 = build_amp2(mi);
    n_.emplace(mkey(mi), build_n(a2));
    for (std::size_t i = 0; i < ms_struck_.size(); ++i) {
      std::vector<double> cdf = build_cdf(a2[i]);
      if (cdf.empty()) continue;   // a CG-forbidden (M, m_S) pair
      cdf_.emplace(cdf_key(mi, ms_struck_[i]), std::move(cdf));
    }
    amp2_.emplace(mkey(mi), std::move(a2));
  }
}

const std::vector<double>& TaggedModel::radial_table(int l) const {
  const auto it = rad_.find(l);
  if (it == rad_.end()) throw std::runtime_error("no wave with that L");
  return it->second;
}

std::size_t TaggedModel::ms_index(double m_s) const {
  for (std::size_t i = 0; i < ms_struck_.size(); ++i) {
    if (std::fabs(ms_struck_[i] - m_s) < 1e-9) return i;
  }
  throw std::runtime_error("m_S is not a channel-spin projection");
}

const std::vector<std::vector<double>>& TaggedModel::amp2_table(
    double m_ion) const {
  const auto it = amp2_.find(mkey(m_ion));
  if (it == amp2_.end()) {
    throw std::runtime_error("TaggedModel: M is not an ion projection of "
                             "this channel");
  }
  return it->second;
}

std::vector<std::vector<double>> TaggedModel::build_amp2(double m_ion) const {
  const std::size_t nk = k_.size(), nc = c_.size();
  std::vector<std::vector<double>> out(ms_struck_.size(),
                                       std::vector<double>(nk * nc, 0.0));
  for (std::size_t i = 0; i < ms_struck_.size(); ++i) {
    const double m_s = ms_struck_[i];
    const double m_l = m_ion - m_s;
    // per-wave prefactors, so the (k, c) loop is two multiplies deep
    std::vector<const std::vector<double>*> rad;
    std::vector<double> cg_of;
    std::vector<int> l_of;
    for (const Wave& w : channel_.waves) {
      if (std::fabs(m_l) > w.l + 0.5) continue;
      const double cg = clebsch_gordan(w.l, m_l, channel_.s_channel, m_s,
                                       channel_.j_ion, m_ion);
      if (cg == 0.0) continue;
      rad.push_back(&rad_.at(w.l));
      cg_of.push_back(cg);
      l_of.push_back(w.l);
    }
    if (rad.empty()) continue;
    const int ml_int = static_cast<int>(std::lround(m_l));
    std::vector<std::vector<double>> ang(rad.size(), std::vector<double>(nc));
    for (std::size_t w = 0; w < rad.size(); ++w) {
      for (std::size_t ic = 0; ic < nc; ++ic) {
        ang[w][ic] = cg_of[w] * theta_lm(l_of[w], ml_int, c_[ic]);
      }
    }
    std::vector<double>& tab = out[i];
    for (std::size_t ik = 0; ik < nk; ++ik) {
      for (std::size_t ic = 0; ic < nc; ++ic) {
        double amp = 0.0;
        for (std::size_t w = 0; w < rad.size(); ++w) {
          amp += (*rad[w])[ik] * ang[w][ic];
        }
        tab[ik * nc + ic] = amp * amp;
      }
    }
  }
  return out;
}

const std::vector<double>& TaggedModel::n_of_kc(double m_ion) const {
  const auto it = n_.find(mkey(m_ion));
  if (it == n_.end()) {
    throw std::runtime_error("TaggedModel: M is not an ion projection of "
                             "this channel");
  }
  return it->second;
}

std::vector<double> TaggedModel::build_n(
    const std::vector<std::vector<double>>& a2) const {
  std::vector<double> out(k_.size() * c_.size(), 0.0);
  for (const auto& t : a2) {
    for (std::size_t j = 0; j < out.size(); ++j) out[j] += t[j];
  }
  return out;
}

double TaggedModel::n_of_kc(double m_ion, double k, double c) const {
  const std::vector<double>& tab = n_of_kc(m_ion);
  return tab[cell_index(k_, k) * c_.size() + cell_index(c_, c)];
}

std::vector<std::vector<double>> TaggedModel::struck_populations(
    double m_ion) const {
  const auto& a2 = amp2_table(m_ion);
  const std::vector<double>& n = n_of_kc(m_ion);
  std::vector<std::vector<double>> out = a2;
  for (auto& t : out) {
    for (std::size_t j = 0; j < t.size(); ++j) t[j] /= std::fmax(n[j], 1e-300);
  }
  return out;
}

std::vector<double> TaggedModel::population_integrated(double m_ion) const {
  const auto& a2 = amp2_table(m_ion);
  const std::size_t nk = k_.size(), nc = c_.size();
  std::vector<double> raw(a2.size(), 0.0);
  std::vector<double> buf(nk * nc);
  for (std::size_t i = 0; i < a2.size(); ++i) {
    for (std::size_t ik = 0; ik < nk; ++ik) {
      const double k2 = k_[ik] * k_[ik];
      for (std::size_t ic = 0; ic < nc; ++ic) buf[ik * nc + ic] = a2[i][ik * nc + ic] * k2;
    }
    raw[i] = pairwise_sum(buf);
  }
  const double tot = pairwise_sum(raw);
  for (double& v : raw) v /= tot;
  return raw;
}

double TaggedModel::norm(double m_ion) const {
  const std::vector<double>& n = n_of_kc(m_ion);
  const std::size_t nk = k_.size(), nc = c_.size();
  std::vector<double> buf(nk * nc);
  for (std::size_t ik = 0; ik < nk; ++ik) {
    const double k2 = k_[ik] * k_[ik];
    for (std::size_t ic = 0; ic < nc; ++ic) buf[ik * nc + ic] = n[ik * nc + ic] * k2;
  }
  return pairwise_sum(buf) * dk_ * 2.0 * kPi * dc_;
}

std::vector<std::vector<double>> TaggedModel::pair_populations(
    double m_s) const {
  const std::vector<double> m1s = m_values(channel_.s_struck);
  const std::vector<double> m2s = m_values(channel_.s_spec);
  std::vector<std::vector<double>> out(m1s.size(), std::vector<double>(m2s.size(), 0.0));
  for (std::size_t i = 0; i < m1s.size(); ++i) {
    for (std::size_t j = 0; j < m2s.size(); ++j) {
      const double cg = clebsch_gordan(channel_.s_struck, m1s[i],
                                       channel_.s_spec, m2s[j],
                                       channel_.s_channel, m_s);
      out[i][j] = cg * cg;
    }
  }
  return out;
}

double TaggedModel::vector_dilution(double m_ion) const {
  const std::vector<double> p = population_integrated(m_ion);
  double s = 0.0;
  for (std::size_t i = 0; i < p.size(); ++i) s += ms_struck_[i] * p[i];
  return s / channel_.s_channel;
}

double TaggedModel::tensor_dilution(double m_ion) const {
  if (std::fabs(channel_.s_channel - 1.0) > 1e-9) {
    throw std::runtime_error("tensor dilution defined for S_c = 1");
  }
  const std::vector<double> p = population_integrated(m_ion);
  double s = 0.0;
  for (std::size_t i = 0; i < p.size(); ++i) {
    s += (3.0 * ms_struck_[i] * ms_struck_[i] - 2.0) * p[i];
  }
  return s;
}

double TaggedModel::p2_moment(double m_ion) const {
  const std::vector<double>& n = n_of_kc(m_ion);
  const std::size_t nk = k_.size(), nc = c_.size();
  std::vector<double> w(nc, 0.0), wp(nc);
  for (std::size_t ik = 0; ik < nk; ++ik) {
    const double k2 = k_[ik] * k_[ik];
    for (std::size_t ic = 0; ic < nc; ++ic) w[ic] += n[ik * nc + ic] * k2;
  }
  for (std::size_t ic = 0; ic < nc; ++ic) {
    wp[ic] = w[ic] * 0.5 * (3.0 * c_[ic] * c_[ic] - 1.0);
  }
  return pairwise_sum(wp) / pairwise_sum(w);
}

double TaggedModel::p2_moment_mixture(
    const std::vector<double>& populations) const {
  const std::vector<double> ms = m_values(channel_.j_ion);
  if (populations.size() != ms.size()) {
    throw std::runtime_error("population vector has the wrong size");
  }
  double out = 0.0;
  for (std::size_t i = 0; i < ms.size(); ++i) {
    if (populations[i] > 0.0) out += populations[i] * p2_moment(ms[i]);
  }
  return out;
}

const std::vector<double>& TaggedModel::cell_cdf(double m_ion,
                                                 double m_s) const {
  const auto it = cdf_.find(cdf_key(m_ion, m_s));
  // A (M, m_S) pair the CG factors forbid has no density at all and is not
  // built; `rates()` gives it weight zero, so the event loop never asks.
  if (it == cdf_.end()) throw std::runtime_error("empty (M, m_S) density");
  return it->second;
}

std::vector<double> TaggedModel::build_cdf(
    const std::vector<double>& a2) const {
  const std::size_t nk = k_.size(), nc = c_.size();
  std::vector<double> cdf(nk * nc);
  double acc = 0.0;
  for (std::size_t ik = 0; ik < nk; ++ik) {
    const double k2 = k_[ik] * k_[ik];
    for (std::size_t ic = 0; ic < nc; ++ic) {
      acc += a2[ik * nc + ic] * k2;
      cdf[ik * nc + ic] = acc;
    }
  }
  if (!(acc > 0.0)) return std::vector<double>();
  for (double& v : cdf) v /= acc;
  return cdf;
}

void TaggedModel::sample_kc_one(double m_ion, double m_s, Rng& rng,
                                double& k_out, double& c_out,
                                double& phi_out) const {
  const std::vector<double>& cdf = cell_cdf(m_ion, m_s);
  const std::size_t nc = c_.size();
  const std::size_t cell = static_cast<std::size_t>(
      std::lower_bound(cdf.begin(), cdf.end(), rng.uniform()) - cdf.begin());
  const std::size_t cell_c = cell < cdf.size() ? cell : cdf.size() - 1;
  const std::size_t ik = cell_c / nc, ic = cell_c % nc;
  k_out = std::fabs(k_[ik] + (rng.uniform() - 0.5) * dk_);
  c_out = clip(c_[ic] + (rng.uniform() - 0.5) * dc_, -1.0, 1.0);
  phi_out = 2.0 * kPi * rng.uniform();
}

void TaggedModel::sample_kc(double m_ion, double m_s, std::size_t n, Rng& rng,
                            std::vector<double>& k_out,
                            std::vector<double>& c_out,
                            std::vector<double>& phi_out) const {
  k_out.resize(n);
  c_out.resize(n);
  phi_out.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    sample_kc_one(m_ion, m_s, rng, k_out[i], c_out[i], phi_out[i]);
  }
}

// ------------------------------------------------------------------- boost

SpectatorLab boost_spectator(const TaggedChannel& channel, double k, double c,
                             double phi_k, double p_per_nucleon,
                             double theta_s, double phi_s) {
  const ClusterChannel& base = channel.base;
  const double s = std::sqrt(std::fmax(1.0 - c * c, 0.0));
  double kx = k * s * std::cos(phi_k);
  double ky = k * s * std::sin(phi_k);
  double kz = k * c;
  if (theta_s != 0.0 || phi_s != 0.0) {
    const double ct = std::cos(theta_s), st = std::sin(theta_s);
    const double cp = std::cos(phi_s), sp = std::sin(phi_s);
    // R_z(phi_s) R_y(theta_s): spin frame -> lab
    const double kx1 = ct * kx + st * kz;
    const double kz1 = -st * kx + ct * kz;
    kx = kx1;
    kz = kz1;
    const double kx2 = cp * kx - sp * ky;
    const double ky2 = sp * kx + cp * ky;
    kx = kx2;
    ky = ky2;
  }
  const double m = base.m_spec();
  const double e_rest = std::sqrt(m * m + kx * kx + ky * ky + kz * kz);
  const double m_beam = base.m_beam();
  const double p_beam = base.beam_A * p_per_nucleon;
  const double e_beam = std::sqrt(p_beam * p_beam + m_beam * m_beam);
  const double gamma = e_beam / m_beam;
  const double gbeta = p_beam / m_beam;
  SpectatorLab out;
  out.pz_lab = gamma * kz + gbeta * e_rest;
  out.e_lab = gamma * e_rest + gbeta * kz;
  out.pT = std::sqrt(kx * kx + ky * ky);
  out.p_lab = std::sqrt(out.pT * out.pT + out.pz_lab * out.pz_lab);
  out.theta = std::atan2(out.pT, out.pz_lab);
  const double rigidity_beam = p_beam / base.beam_Z;
  out.R = base.spectator_Z > 0 ? (out.p_lab / base.spectator_Z) / rigidity_beam
                               : kNaN;
  out.xL = out.p_lab / (base.spectator_A * p_per_nucleon);
  out.kx = kx;
  out.ky = ky;
  out.kz = kz;
  out.phi_spec = std::atan2(ky, kx);
  return out;
}

StruckCluster struck_cluster(const TaggedChannel& channel,
                             const SpectatorLab& lab, double p_per_nucleon) {
  const ClusterChannel& base = channel.base;
  const double m_beam = base.m_beam();
  const double p_beam = base.beam_A * p_per_nucleon;
  StruckCluster out;
  out.p_ion = Vec4{std::sqrt(p_beam * p_beam + m_beam * m_beam), 0.0, 0.0, p_beam};
  out.p_spectator = Vec4{lab.e_lab, lab.kx, lab.ky, lab.pz_lab};
  out.p = out.p_ion - out.p_spectator;
  out.m2 = out.p.m2();
  out.m_free = base.m_partner();
  out.virtuality = out.m2 - out.m_free * out.m_free;
  const double m_spec = base.m_spec();
  const double k2 = lab.kx * lab.kx + lab.ky * lab.ky + lab.kz * lab.kz;
  const double e_rest = std::sqrt(m_spec * m_spec + k2);
  out.alpha_s = base.beam_A * (e_rest - lab.kz) / m_beam;
  out.alpha_x = base.beam_A - out.alpha_s;
  out.pt_s = std::sqrt(lab.kx * lab.kx + lab.ky * lab.ky);
  out.p_per_nucleon_eff = out.p.pz / base.partner_A();
  return out;
}

// ------------------------------------------------------- wave-function A_zz

namespace {

std::vector<double> azz_from_columns(const std::vector<double>& np1,
                                     const std::vector<double>& n0,
                                     const std::vector<double>& nm1) {
  std::vector<double> out(np1.size());
  for (std::size_t i = 0; i < out.size(); ++i) {
    const double num = np1[i] + nm1[i] - 2.0 * n0[i];
    const double den = np1[i] + nm1[i] + n0[i];
    out[i] = den > 0.0 ? num / den : kNaN;
  }
  return out;
}

}  // namespace

std::vector<double> azz_tensor_curve(const TaggedModel& model, std::size_t ic) {
  const std::size_t nk = model.nk(), nc = model.nc();
  std::vector<double> col[3];
  const double ms[3] = {1.0, 0.0, -1.0};
  for (int j = 0; j < 3; ++j) {
    const std::vector<double>& n = model.n_of_kc(ms[j]);
    col[j].resize(nk);
    for (std::size_t ik = 0; ik < nk; ++ik) col[j][ik] = n[ik * nc + ic];
  }
  return azz_from_columns(col[0], col[1], col[2]);
}

std::vector<double> azz_tensor_curve_weighted(
    const TaggedModel& model, const std::vector<double>& weights) {
  const std::size_t nk = model.nk(), nc = model.nc();
  if (weights.size() != nk * nc) {
    throw std::runtime_error("weights must be an (nk, nc) table");
  }
  std::vector<double> col[3];
  const double ms[3] = {1.0, 0.0, -1.0};
  for (int j = 0; j < 3; ++j) {
    const std::vector<double>& n = model.n_of_kc(ms[j]);
    col[j].assign(nk, 0.0);
    for (std::size_t ik = 0; ik < nk; ++ik) {
      double s = 0.0;
      for (std::size_t ic = 0; ic < nc; ++ic) s += n[ik * nc + ic] * weights[ik * nc + ic];
      col[j][ik] = s;
    }
  }
  return azz_from_columns(col[0], col[1], col[2]);
}

std::vector<double> acceptance_weights(const TaggedModel& model,
                                       double p_per_nucleon,
                                       const Optics& optics,
                                       const std::string& pot_config,
                                       std::size_t n_phi, double theta_s,
                                       double phi_s) {
  const std::size_t nk = model.nk(), nc = model.nc();
  std::vector<double> out(nk * nc, 0.0);
  for (std::size_t ik = 0; ik < nk; ++ik) {
    for (std::size_t ic = 0; ic < nc; ++ic) {
      std::size_t hits = 0;
      for (std::size_t ip = 0; ip < n_phi; ++ip) {
        const double phi = (static_cast<double>(ip) + 0.5) * 2.0 * kPi
                           / static_cast<double>(n_phi);
        const SpectatorLab lab =
            boost_spectator(model.channel(), model.k()[ik], model.c()[ic], phi,
                            p_per_nucleon, theta_s, phi_s);
        const int r = route_charged(lab.R, lab.theta, lab.pT, optics,
                                    lab.phi_spec, kNaN, pot_config);
        if (rp_accepted(r)) ++hits;
      }
      out[ik * nc + ic] = static_cast<double>(hits) / static_cast<double>(n_phi);
    }
  }
  return out;
}

// ----------------------------------------------------------- TaggedSampler

int nuclide_pdg(int z, int a) {
  if (z == 1 && a == 1) return 2212;
  if (z == 0 && a == 1) return 2112;
  return 1000000000 + z * 10000 + a * 10;
}

TaggedSampler::TaggedSampler(const TaggedModel& model, double p_per_nucleon,
                             KinematicsSource* dis)
    : model_(&model), p_u_(p_per_nucleon), dis_(dis),
      ms_ion_(m_values(model.channel().j_ion)) {
  const int beam_z = model.channel().base.beam_Z;
  const int beam_a = model.channel().base.beam_A;
  const std::string ion_name = beam_a == 6 ? "6Li"
                             : beam_a == 7 ? "7Li"
                             : beam_z == 2 ? "3He" : "d";
  optics_ = yr_optics(ion_name, p_per_nucleon, true);
  pot_config_ = yr_config_key(ion_name, p_per_nucleon);
}

TaggedSampler::TaggedSampler(const TaggedModel& model, double p_per_nucleon,
                             KinematicsSource* dis, const Optics& optics,
                             const std::string& pot_config)
    : model_(&model), p_u_(p_per_nucleon), dis_(dis), optics_(optics),
      pot_config_(pot_config), ms_ion_(m_values(model.channel().j_ion)) {}

std::vector<double> TaggedSampler::rates(const IonFill& fill) const {
  const std::vector<double>& ms_ion = ms_ion_;
  const std::vector<double>& ms_c = model_->m_struck_values();
  if (fill.populations.size() != ms_ion.size()) {
    throw std::runtime_error("fill populations do not match the ion spin");
  }
  std::vector<double> out(ms_ion.size() * ms_c.size(), 0.0);
  for (std::size_t a = 0; a < ms_ion.size(); ++a) {
    if (fill.populations[a] <= 0.0) continue;
    const std::vector<double> p_ms = model_->population_integrated(ms_ion[a]);
    for (std::size_t b = 0; b < ms_c.size(); ++b) {
      if (p_ms[b] <= 0.0) continue;
      const double sig = dis_ ? dis_->sigma_tot_pb(ms_c[b], fill.lam_e, fill.pe)
                              : 1.0;
      out[a * ms_c.size() + b] = fill.populations[a] * p_ms[b] * sig;
    }
  }
  return out;
}

double TaggedSampler::sigma_tot_pb(const IonFill& fill) const {
  const std::vector<double> r = rates(fill);
  return pairwise_sum(r);
}

std::vector<double> TaggedSampler::rate_cdf(const IonFill& fill) const {
  std::vector<double> r = rates(fill);
  double acc = 0.0;
  for (double& v : r) {
    acc += v;
    v = acc;
  }
  if (acc <= 0.0) throw std::runtime_error("the fill has no rate at all");
  for (double& v : r) v /= acc;
  return r;
}

void TaggedSampler::set_fsi(std::shared_ptr<const FsiWeight> fsi) {
  fsi_ = std::move(fsi);
}

void TaggedSampler::apply_fsi(TaggedEvent& te) const {
  // A pure function of the DRAWN kinematics: no RNG is consumed, so the
  // event stream, its counter-based determinism and thread safety are
  // exactly what they are without an FSI model.  The weight never moves a
  // four-vector (fsi.hpp: FSI IS A WEIGHT, NEVER A SHIFT).
  if (!fsi_) return;
  const ClusterChannel& base = model_->channel().base;
  FsiKinematics kin;
  kin.k = te.k;
  kin.cos_theta_k = te.cos_theta_k;
  kin.phi_k = te.phi_k;
  kin.w = std::sqrt(std::fmax(w2_from_xq2(te.x, te.q2), 0.0));
  kin.q2 = te.q2;
  kin.x = te.x;
  kin.spectator_z = base.spectator_Z;
  kin.spectator_a = base.spectator_A;
  kin.channel = channel_enum_of(base);
  te.weight = fsi_->weight(kin);
}

TaggedEvent TaggedSampler::sample_one(const IonFill& fill,
                                      const std::vector<double>& cdf,
                                      Rng& rng) const {
  const std::vector<double>& ms_ion = ms_ion_;
  const std::vector<double>& ms_c = model_->m_struck_values();
  std::size_t cell = static_cast<std::size_t>(
      std::lower_bound(cdf.begin(), cdf.end(), rng.uniform()) - cdf.begin());
  if (cell >= cdf.size()) cell = cdf.size() - 1;
  TaggedEvent te;
  te.m_ion = ms_ion[cell / ms_c.size()];
  te.m_struck = ms_c[cell % ms_c.size()];
  model_->sample_kc_one(te.m_ion, te.m_struck, rng, te.k, te.cos_theta_k,
                        te.phi_k);
  te.lab = boost_spectator(model_->channel(), te.k, te.cos_theta_k, te.phi_k,
                           p_u_, fill.theta_s, fill.phi_s);
  te.route = route_charged(te.lab.R, te.lab.theta, te.lab.pT, optics_,
                           te.lab.phi_spec, kNaN, pot_config_);
  if (dis_) {
    // P6.  `KinematicsSource::sample` fills four vectors; allocating them per
    // event cost four malloc/free pairs on the hot path.  THREAD-LOCAL
    // scratch: `sample_one` is const and is called concurrently by
    // `Pipeline::for_each(sink, nthreads)`, so the buffers cannot be members,
    // and their CONTENTS never survive the call, so nothing about the event
    // stream depends on them.
    static thread_local std::vector<double> x, q2, y, phi;
    static thread_local std::vector<int> dis_cell;
    dis_->sample_cells(te.m_struck, fill.lam_e, fill.pe, 1, rng, x, q2, y, phi,
                       dis_cell);
    te.x = x[0];
    te.q2 = q2[0];
    te.y = y[0];
    te.phi = phi[0];
    te.cell = dis_cell[0];
  }
  apply_fsi(te);   // after the DIS draw: the W ramp reads (x, Q2)
  return te;
}

std::vector<TaggedEvent> TaggedSampler::sample_category(const IonFill& fill,
                                                        std::size_t n,
                                                        Rng& rng) const {
  const std::vector<double>& ms_ion = ms_ion_;
  const std::vector<double>& ms_c = model_->m_struck_values();
  const std::vector<double> cdf = rate_cdf(fill);

  // draw the (M, m_S) cell of every event first, so the DIS source is called
  // once per cell rather than once per event (the Python's grouping)
  std::vector<std::size_t> counts(cdf.size(), 0);
  for (std::size_t i = 0; i < n; ++i) {
    std::size_t cell = static_cast<std::size_t>(
        std::lower_bound(cdf.begin(), cdf.end(), rng.uniform()) - cdf.begin());
    if (cell >= cdf.size()) cell = cdf.size() - 1;
    ++counts[cell];
  }

  std::vector<TaggedEvent> out;
  out.reserve(n);
  std::vector<double> k, c, phi_k, x, q2, y, phi;
  std::vector<int> dis_cell;
  for (std::size_t cell = 0; cell < counts.size(); ++cell) {
    const std::size_t cnt = counts[cell];
    if (cnt == 0) continue;
    const double m_ion = ms_ion[cell / ms_c.size()];
    const double m_s = ms_c[cell % ms_c.size()];
    model_->sample_kc(m_ion, m_s, cnt, rng, k, c, phi_k);
    if (dis_)
      dis_->sample_cells(m_s, fill.lam_e, fill.pe, cnt, rng, x, q2, y, phi,
                         dis_cell);
    for (std::size_t i = 0; i < cnt; ++i) {
      TaggedEvent te;
      te.m_ion = m_ion;
      te.m_struck = m_s;
      te.k = k[i];
      te.cos_theta_k = c[i];
      te.phi_k = phi_k[i];
      te.lab = boost_spectator(model_->channel(), te.k, te.cos_theta_k,
                               te.phi_k, p_u_, fill.theta_s, fill.phi_s);
      te.route = route_charged(te.lab.R, te.lab.theta, te.lab.pT, optics_,
                               te.lab.phi_spec, kNaN, pot_config_);
      if (dis_) {
        te.x = x[i];
        te.q2 = q2[i];
        te.y = y[i];
        te.phi = phi[i];
        te.cell = dis_cell[i];
      }
      apply_fsi(te);
      out.push_back(te);
    }
  }
  return out;
}

void TaggedSampler::fill_event(Event& ev, const TaggedEvent& te) const {
  const TaggedChannel& ch = model_->channel();
  const ClusterChannel& base = ch.base;
  const StruckCluster sc = struck_cluster(ch, te.lab, p_u_);

  int mother = -1;
  for (std::size_t i = 0; i < ev.particles.size(); ++i) {
    if (ev.particles[i].role == Role::BeamIon) {
      mother = static_cast<int>(i);
      break;
    }
  }

  Particle spec;
  spec.pdg = nuclide_pdg(base.spectator_Z, base.spectator_A);
  spec.status = Status::Final;
  spec.role = Role::Spectator;
  spec.p = sc.p_spectator;
  spec.mass = base.m_spec();
  spec.charge = base.spectator_Z;
  spec.mother1 = mother;
  ev.particles.push_back(spec);

  Particle x;
  x.pdg = nuclide_pdg(base.partner_Z(), base.partner_A());
  x.status = Status::Intermediate;
  x.role = Role::StruckCluster;
  x.p = sc.p;
  // The struck cluster is off shell by construction (see StruckCluster): the
  // generated mass is its ACTUAL invariant mass, not the free-cluster one.
  x.mass = sc.m2 > 0.0 ? std::sqrt(sc.m2) : -std::sqrt(-sc.m2);
  x.charge = base.partner_Z();
  x.mother1 = mother;
  ev.particles.push_back(x);

  ev.kin.k = te.k;
  ev.kin.cos_theta_k = te.cos_theta_k;
  ev.kin.phi_k = te.phi_k;
  ev.kin.alpha_s = sc.alpha_s;
  ev.kin.pt_s = sc.pt_s;
  // The spectator's LAB block, STORED rather than left to be re-derived by
  // every consumer from the four-vector (`Kinematics`): these are the exact
  // numbers `boost_spectator` computed, so nothing downstream has to invert
  // the beam boost or re-do the rigidity algebra to get them back.
  ev.kin.spec_pt = te.lab.pT;
  ev.kin.spec_theta = te.lab.theta;
  ev.kin.spec_p_lab = te.lab.p_lab;
  ev.kin.spec_r = te.lab.R;
  ev.kin.spec_xl = te.lab.xL;
  ev.kin.spec_kx = te.lab.kx;
  ev.kin.spec_ky = te.lab.ky;
  ev.kin.spec_kz = te.lab.kz;
  ev.kin.phi_spec = te.lab.phi_spec;
  ev.spin.m_ion = te.m_ion;
  ev.spin.m_struck = te.m_struck;
  ev.channel = channel_enum_of(base);
}

}  // namespace lipolgen
