// SPDX-License-Identifier: GPL-3.0-or-later
// Tier T1 -- struck cluster -> struck nucleon + partner spectator(s).
// Every physics choice, formula and separation energy is in breakup.hpp.

#include "lipolgen/breakup.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "lipolgen/numerics.hpp"
#include "lipolgen/spin.hpp"

namespace lipolgen {
namespace {

/// Active boost of a rest-frame vector into the frame where the system `ref`
/// (E > 0) moves.  Uses `ref`'s own invariant mass, so it is exact for the
/// OFF-SHELL P_X as long as P_X^2 > 0 (which the tagged channels guarantee --
/// `Pipeline::make_tagged` rejects the draw otherwise).
Vec4 from_rest_of(const Vec4& v, const Vec4& ref, double m_ref) {
  const double bx = ref.px / ref.e, by = ref.py / ref.e, bz = ref.pz / ref.e;
  const double b2 = bx * bx + by * by + bz * bz;
  if (!(b2 > 0.0)) return v;
  const double gamma = ref.e / m_ref;
  const double bp = bx * v.px + by * v.py + bz * v.pz;
  const double f = (gamma - 1.0) * bp / b2 + gamma * v.e;
  return {gamma * (v.e + bp), v.px + f * bx, v.py + f * by, v.pz + f * bz};
}

/// An on-shell four-vector of mass `m` with three-momentum `p` along
/// (sin t cos f, sin t sin f, cos t).
Vec4 on_shell(double m, double p, double cos_t, double phi) {
  const double st = std::sqrt(std::max(0.0, 1.0 - cos_t * cos_t));
  return {std::sqrt(m * m + p * p), p * st * std::cos(phi),
          p * st * std::sin(phi), p * cos_t};
}

/// One isotropic direction from two uniforms.
void isotropic(Rng& rng, double* cos_t, double* phi) {
  *cos_t = 2.0 * rng.uniform() - 1.0;
  *phi = 2.0 * kPi * rng.uniform();
}

Particle make_fragment(int z, int a, const Vec4& p) {
  Particle f;
  f.pdg = nuclide_pdg(z, a);
  f.status = Status::Final;
  f.role = Role::PartnerSpectator;
  f.p = p;
  f.mass = nuclear_mass(z, a);
  f.charge = static_cast<double>(z);
  f.mother1 = -1;
  f.pol = 9.0;
  return f;
}

/// +-1 helicity label from a vector polarization P in [-1, 1].
double sample_pol(double p_vec, Rng& rng) {
  const double pp = std::max(-1.0, std::min(1.0, p_vec));
  return (rng.uniform() < 0.5 * (1.0 + pp)) ? 1.0 : -1.0;
}

/// Draw the projection m_1 of nucleon 1 of a two-nucleon pair coupled to
/// channel spin S_c with projection m_sc, from the CG factor
/// |<1/2 m_1 1/2 m_2 | S_c m_sc>|^2 (m_2 = m_sc - m_1).  Returns +-1/2; the
/// partner's projection is m_sc - m_1, so the pair stays CORRELATED --
/// drawing the two independently would lose the m_sc = 0 anticorrelation.
double draw_pair_m1(double s_channel, double m_sc, Rng& rng) {
  double w[2] = {0.0, 0.0};
  const double m1s[2] = {0.5, -0.5};
  for (int i = 0; i < 2; ++i) {
    const double m2 = m_sc - m1s[i];
    if (std::fabs(m2) > 0.5 + 1e-9) continue;   // outside spin 1/2
    const double cg = clebsch_gordan(0.5, m1s[i], 0.5, m2, s_channel, m_sc);
    w[i] = cg * cg;
  }
  const double tot = w[0] + w[1];
  if (!(tot > 0.0)) return (rng.uniform() < 0.5) ? 0.5 : -0.5;
  return (rng.uniform() * tot < w[0]) ? 0.5 : -0.5;
}

}  // namespace

ClusterSpecies cluster_species(int z, int a) {
  if (a == 1 && (z == 0 || z == 1)) return ClusterSpecies::Nucleon;
  if (a == 2 && z == 1) return ClusterSpecies::Deuteron;
  if (a == 3 && z == 1) return ClusterSpecies::Triton;
  throw std::runtime_error("ClusterBreakup: no breakup model for (Z, A) = (" +
                           std::to_string(z) + ", " + std::to_string(a) + ")");
}

std::vector<double> ClusterBreakup::build_cdf(const std::vector<double>& grid,
                                              double kappa, double beta) {
  std::vector<double> cdf(grid.size(), 0.0);
  double acc = 0.0;
  for (std::size_t i = 0; i < grid.size(); ++i) {
    const double g = grid[i];
    // k^2 |psi_{L=0}(k)|^2 -- `momentum_density` is |psi|^2 without the
    // phase space, exactly as `MomentumSampler` builds its own CDF.
    acc += g * g * momentum_density(g, kappa, beta, 0);
    cdf[i] = acc;
  }
  const double tot = cdf.back();
  if (tot > 0.0) {
    for (double& v : cdf) v /= tot;
  }
  return cdf;
}

double ClusterBreakup::draw_k(const std::vector<double>& cdf, double u) const {
  // `np.interp(u, cdf, grid)`, the same inverse-CDF read `MomentumSampler`
  // uses, so the two agree wherever they share a (kappa, beta).
  const auto it = std::lower_bound(cdf.begin(), cdf.end(), u);
  if (it == cdf.begin()) return grid_.front();
  if (it == cdf.end()) return grid_.back();
  const std::size_t i = static_cast<std::size_t>(it - cdf.begin());
  const double c0 = cdf[i - 1], c1 = cdf[i];
  const double t = (c1 > c0) ? (u - c0) / (c1 - c0) : 0.0;
  return grid_[i - 1] + t * (grid_[i] - grid_[i - 1]);
}

ClusterBreakup::ClusterBreakup(BreakupOptions opt) : opt_(std::move(opt)) {
  if (!(opt_.beta > 0.0)) throw std::runtime_error("ClusterBreakup: beta <= 0");
  if (!(opt_.k_max > 0.0)) throw std::runtime_error("ClusterBreakup: k_max <= 0");
  if (opt_.n_grid < 16) throw std::runtime_error("ClusterBreakup: n_grid too small");
  f2_ = opt_.f2 ? opt_.f2
                : std::static_pointer_cast<const UnpolSF>(
                      std::make_shared<const ToyF2>());
  tsf_ = opt_.triton_sf;   // null = the sequential triton branch, bit for bit

  // --- the deuteron's own wave function ------------------------------------
  // S + D Hulthen at `p_d` (= P_D_DEUTERON by default), or the exact AV18
  // deuteron on `source` = VmcAV18.  IT FOLLOWS `source` BECAUSE THE RATE
  // DOES: `li6_alpha_channel`'s `dis_target` is `DEUTERON_AV18()` on that
  // setting, and a run whose struck-nucleon spin label came from one deuteron
  // and whose g1 came from another would be inconsistent with itself by
  // 2.0688 % (C5.5b).  `Hulthen` is bit for bit what this line always built.
  dmodel_ = std::unique_ptr<TaggedModel>(new TaggedModel(
      deuteron_channel(opt_.beta, opt_.p_d, opt_.source), opt_.k_max, opt_.nk,
      opt_.nc));
  ms_deuteron_ = m_values(1.0);
  dpop_.resize(ms_deuteron_.size());
  for (std::size_t i = 0; i < ms_deuteron_.size(); ++i) {
    // P(m_sc | M) -- the k- and khat-integrated channel-spin populations.
    // Built ONCE here: `population_integrated` is an O(nk nc) quadrature and
    // must never run on the per-event path.
    dpop_[i] = dmodel_->population_integrated(ms_deuteron_[i]);
  }

  // --- the triton's two sequential channels -------------------------------
  const double m_p = nuclear_mass(1, 1);
  const double m_n = nuclear_mass(0, 1);
  const double m_d = nuclear_mass(1, 2);
  const double m_t = nuclear_mass(1, 3);
  m_nn_ref_ = 2.0 * m_n;
  grid_ = linspace(1e-4, opt_.k_max, opt_.n_grid);

  // t* -> n + d : S_n(3H) = m_n + m_d - m_t = 6.2572 MeV (AME2020).
  {
    const double s = m_n + m_d - m_t;
    const double mu = m_n * m_d / (m_n + m_d);
    cdf_t_nd_ = build_cdf(grid_, std::sqrt(2.0 * mu * s), opt_.beta);
  }
  // t* -> p + (nn) : the three-body threshold m_p + 2 m_n - m_t = 8.4818 MeV,
  // with the nn pair treated as one fragment of mass 2 m_n for the purpose of
  // the reduced mass (the pair's own relative energy is drawn independently
  // below -- the crude step this file flags).
  {
    const double s = m_p + 2.0 * m_n - m_t;
    const double mu = m_p * m_nn_ref_ / (m_p + m_nn_ref_);
    cdf_t_pnn_ = build_cdf(grid_, std::sqrt(2.0 * mu * s), opt_.beta);
  }
  // The unbound nn pair, at its virtual-state pole 1/|a_nn|.
  cdf_nn_ = build_cdf(grid_, opt_.kappa_nn, opt_.beta);
}

double ClusterBreakup::proton_fraction(double x, double q2, int z, int a) const {
  const int n = a - z;
  const double flat = (a > 0) ? static_cast<double>(z) / static_cast<double>(a)
                              : 0.0;
  if (!f2_ || !(x > 0.0) || !(q2 > 0.0)) return flat;
  const double zp = static_cast<double>(z) * f2_->f2p(x, q2);
  const double nn = static_cast<double>(n) * f2_->f2n(x, q2);
  const double tot = zp + nn;
  return (tot > 0.0) ? zp / tot : flat;
}

bool ClusterBreakup::resolve(const BreakupInput& in, Rng& rng,
                             BreakupResult& out) const {
  out.partners.clear();
  out.k = out.cos_theta_k = out.phi_k = out.q_nn = 0.0;
  out.m_remnant = 0.0;
  out.e_rel = 0.0;
  out.channel = TritonChannel::NeutronD;

  if (!(in.p_cluster.e > 0.0)) return false;

  out.struck = Particle();
  out.struck.status = Status::Intermediate;
  out.struck.role = Role::StruckNucleon;
  out.struck.mother1 = -1;

  const double m_p = nuclear_mass(1, 1);
  const double m_n = nuclear_mass(0, 1);
  const double m_d = nuclear_mass(1, 2);

  if (in.species == ClusterSpecies::Nucleon) {
    // Nothing to break up: the "cluster" IS the struck nucleon (the d and 3He
    // control channels).  It keeps P_X verbatim -- off shell, because the
    // tagged spectator was put on shell -- and only the labels are new.
    out.struck.pdg = (in.z >= 1) ? 2212 : 2112;
    out.struck.charge = (in.z >= 1) ? 1.0 : 0.0;
    out.struck.p = in.p_cluster;
    // Spin.  A nucleon that came out of a two-nucleon channel of spin S_c
    // carries the CG projection of that channel spin; a free spin-1/2 struck
    // object simply carries its own m_S.
    if (in.s_cluster > 0.75) {
      out.struck.pol = 2.0 * draw_pair_m1(in.s_cluster, in.m_s, rng);
    } else {
      out.struck.pol = sample_pol(2.0 * in.m_s, rng);
    }
    out.struck.mass = out.struck.p.m2() >= 0.0
                          ? std::sqrt(out.struck.p.m2())
                          : -std::sqrt(-out.struck.p.m2());
    out.virtuality = out.struck.p.m2() - M_NUCLEON * M_NUCLEON;
    out.m_remnant = 0.0;
    return true;
  }

  const double m2_x = in.p_cluster.m2();
  if (!(m2_x > 0.0)) return false;   // no rest frame to break up in
  const double m_x = std::sqrt(m2_x);

  if (in.species == ClusterSpecies::Deuteron) {
    // 1. the channel-spin projection m_sc of the np pair, given the
    //    deuteron's own projection m_S.
    std::size_t i_m = 0;
    for (std::size_t i = 0; i < ms_deuteron_.size(); ++i) {
      if (std::fabs(ms_deuteron_[i] - in.m_s) < 1e-9) { i_m = i; break; }
    }
    const std::vector<double>& pop = dpop_[i_m];
    const double us = rng.uniform();
    double acc = 0.0;
    std::size_t i_sc = pop.size() - 1;
    for (std::size_t i = 0; i < pop.size(); ++i) {
      acc += pop[i];
      if (us < acc) { i_sc = i; break; }
    }
    const double m_sc = dmodel_->m_struck_values()[i_sc];

    // 2. (k, cos theta_k, phi_k) from |A_{m_sc}(m_S; k, khat)|^2 k^2 -- the
    //    D-wave's angular correlation with the spin included.
    dmodel_->sample_kc_one(in.m_s, m_sc, rng, out.k, out.cos_theta_k,
                           out.phi_k);

    // 3. which nucleon is struck: F2p : F2n at the event's own (x, Q2).
    const bool struck_p = rng.uniform() < proton_fraction(in.x, in.q2, 1, 2);
    const int z_struck = struck_p ? 1 : 0;
    const int z_part = 1 - z_struck;
    const double m_part = struck_p ? m_n : m_p;

    // 4. the partner on shell at +k in the P_X rest frame.
    const Vec4 part_rest = on_shell(m_part, out.k, out.cos_theta_k, out.phi_k);
    Particle partner = make_fragment(z_part, 1,
                                     from_rest_of(part_rest, in.p_cluster, m_x));
    // 5. the two nucleon projections, drawn TOGETHER from the pair CG factor
    //    (m_1 + m_2 = m_sc), so the struck nucleon and its partner keep the
    //    m_sc = 0 anticorrelation the deuteron's triplet actually has.
    const double m1 = draw_pair_m1(1.0, m_sc, rng);
    partner.pol = 2.0 * (m_sc - m1);

    out.m_remnant = m_part;
    out.partners.push_back(partner);
    out.struck.pdg = struck_p ? 2212 : 2112;
    out.struck.charge = struck_p ? 1.0 : 0.0;
    out.struck.pol = 2.0 * m1;
  } else if (tsf_) {
    // ---- the triton, the SPECTRAL-FUNCTION branch (see breakup.hpp) -----
    //
    // Exactly nine uniforms, whatever the channel: six inside `sample()`, two
    // for the pair's own direction and one for the polarization label.  The
    // two direction uniforms are drawn even when the remnant is the bound
    // deuteron and there is no pair to point -- that is the whole point.
    const TritonDraw d = tsf_->sample(proton_fraction(in.x, in.q2, 1, 3), rng);
    out.channel = d.channel;
    out.k = d.k;
    out.cos_theta_k = d.cos_theta_k;
    out.phi_k = d.phi_k;
    out.q_nn = d.q_pair;
    out.e_rel = d.e_rel;
    out.m_remnant = d.m_remnant;

    // The remnant cannot be given more energy than the off-shell cluster has:
    // the struck nucleon would come out with E < 0 and the record would be
    // nonsense.  The CS distribution puts 3.3e-5 of its strength above the
    // 1.2 GeV grid ceiling and the cluster carries ~2.8 GeV, so this is a
    // guard, not a physics cut -- but the caller REDRAWS rather than clipping,
    // which is the library's rule for an unusable draw.
    const double e_rem = std::sqrt(d.m_remnant * d.m_remnant + d.k * d.k);
    if (!(e_rem < m_x)) return false;

    const Vec4 rem_rest = on_shell(d.m_remnant, d.k, d.cos_theta_k, d.phi_k);
    const Vec4 rem = from_rest_of(rem_rest, in.p_cluster, m_x);

    double cq = 0.0, pq = 0.0;
    isotropic(rng, &cq, &pq);          // consumed on EVERY channel

    if (d.channel == TritonChannel::NeutronD) {
      out.partners.push_back(make_fragment(1, 2, rem));
    } else {
      // Back to back in the pair's OWN rest frame: the pn pair's two masses
      // differ, so the two on-shell energies do too, but the three-momenta
      // still cancel and M = sum sqrt(m_i^2 + q^2) exactly.
      const bool nn = (d.channel == TritonChannel::ProtonNnCont);
      const int z1 = nn ? 0 : 1;
      const double m1 = nn ? m_n : m_p;
      const Vec4 a_rest = on_shell(m1, d.q_pair, cq, pq);
      const Vec4 b_rest = on_shell(m_n, d.q_pair, -cq, pq + kPi);
      out.partners.push_back(
          make_fragment(z1, 1, from_rest_of(a_rest, rem, d.m_remnant)));
      out.partners.push_back(
          make_fragment(0, 1, from_rest_of(b_rest, rem, d.m_remnant)));
    }

    out.struck.pdg = d.struck_proton ? 2212 : 2112;
    out.struck.charge = d.struck_proton ? 1.0 : 0.0;
    out.struck.pol = sample_pol(
        2.0 * in.m_s * (d.struck_proton ? in.eff_pol_p : in.eff_pol_n), rng);
  } else {
    // ---- the triton, the crude branch (see breakup.hpp) -----------------
    const bool struck_p = rng.uniform() < proton_fraction(in.x, in.q2, 1, 3);
    Vec4 rem_rest;
    double m_rem = 0.0;
    if (struck_p) {
      // t* -> p + (nn): draw the nn relative momentum first, which fixes the
      // remnant's invariant mass exactly.
      out.q_nn = draw_k(cdf_nn_, rng.uniform());
      m_rem = 2.0 * std::sqrt(m_n * m_n + out.q_nn * out.q_nn);
      out.k = draw_k(cdf_t_pnn_, rng.uniform());
    } else {
      m_rem = m_d;
      out.k = draw_k(cdf_t_nd_, rng.uniform());
    }
    isotropic(rng, &out.cos_theta_k, &out.phi_k);
    rem_rest = on_shell(m_rem, out.k, out.cos_theta_k, out.phi_k);
    const Vec4 rem = from_rest_of(rem_rest, in.p_cluster, m_x);
    out.m_remnant = m_rem;

    if (struck_p) {
      // Split the nn pair back to back in ITS own rest frame.
      double cq = 0.0, pq = 0.0;
      isotropic(rng, &cq, &pq);
      const Vec4 n1_rest = on_shell(m_n, out.q_nn, cq, pq);
      const Vec4 n2_rest = on_shell(m_n, out.q_nn, -cq, pq + kPi);
      out.partners.push_back(make_fragment(0, 1, from_rest_of(n1_rest, rem, m_rem)));
      out.partners.push_back(make_fragment(0, 1, from_rest_of(n2_rest, rem, m_rem)));
      out.channel = TritonChannel::ProtonNnCont;
      out.struck.pdg = 2212;
      out.struck.charge = 1.0;
      out.struck.pol = sample_pol(2.0 * in.m_s * in.eff_pol_p, rng);
    } else {
      out.partners.push_back(make_fragment(1, 2, rem));
      out.channel = TritonChannel::NeutronD;
      out.struck.pdg = 2112;
      out.struck.charge = 0.0;
      out.struck.pol = sample_pol(2.0 * in.m_s * in.eff_pol_n, rng);
    }
  }

  // The impulse-approximation rule: the partners are physical, the struck
  // nucleon is whatever is left.  Summed from the EMITTED four-vectors, so
  // the balance is exact to rounding whatever the intermediate algebra did.
  Vec4 p_struck = in.p_cluster;
  for (const Particle& f : out.partners) p_struck = p_struck - f.p;
  out.struck.p = p_struck;
  const double m2 = p_struck.m2();
  out.struck.mass = (m2 >= 0.0) ? std::sqrt(m2) : -std::sqrt(-m2);
  out.virtuality = m2 - M_NUCLEON * M_NUCLEON;
  return true;
}

}  // namespace lipolgen
