#include "lipolgen/spectator.hpp"

#include <cmath>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <utility>

#include "lipolgen/beams.hpp"
#include "lipolgen/constants.hpp"
#include "lipolgen/numerics.hpp"

namespace lipolgen {
namespace {

const std::map<std::string, double>& masses_table() {
  // `spectator.MASSES`, verbatim: the same AME2020-derived evaluation as
  // NUCLEUS_MASS rounded to 10 keV.  Kept rather than re-derived because it
  // is what sets the spectator mass in every published distribution, and
  // 5e-6 relative is far below anything the cluster model resolves.
  static const std::map<std::string, double> t = {
      {"p", 0.93827},   {"n", 0.93957},     {"d", 1.87561},
      {"t", 2.80892},   {"3He", 2.80839},   {"alpha", 3.72738}};
  return t;
}

const std::map<std::pair<int, int>, double>& nucleus_table() {
  // `spectator.NUCLEUS_MASS`, keyed (Z, A).
  static const std::map<std::pair<int, int>, double> t = {
      {{0, 1}, 0.939565420},  {{1, 1}, 0.938272088},  {{1, 2}, 1.875612942},
      {{1, 3}, 2.808921133},  {{2, 3}, 2.808391607},  {{2, 4}, 3.727379407},
      {{3, 6}, 5.601518702},  {{3, 7}, 6.533833028}};
  return t;
}

const double kNaN = std::nan("");

}  // namespace

double cluster_mass(const std::string& name) {
  const auto& t = masses_table();
  const auto it = t.find(name);
  if (it == t.end()) throw std::runtime_error("unknown cluster mass " + name);
  return it->second;
}

bool nuclear_mass_known(int z, int a) {
  return nucleus_table().count(std::make_pair(z, a)) != 0;
}

double nuclear_mass(int z, int a) {
  const auto& t = nucleus_table();
  const auto it = t.find(std::make_pair(z, a));
  if (it != t.end()) return it->second;
  return static_cast<double>(a) * M_U;
}

double ClusterChannel::m_partner() const {
  const int z = beam_Z - spectator_Z;
  const int a = beam_A - spectator_A;
  if (nuclear_mass_known(z, a)) return nuclear_mass(z, a);
  return m_beam() - m_spec() + separation_energy;
}

double ClusterChannel::kappa() const {
  const double ms = m_spec(), mp = m_partner();
  const double mu = ms * mp / (ms + mp);
  return std::sqrt(2.0 * mu * separation_energy);
}

double ClusterChannel::r_at_k_zero() const {
  if (spectator_Z <= 0) return kNaN;
  return (m_spec() / spectator_Z) / (m_beam() / beam_Z);
}

#define LIPOLGEN_CHANNEL(fn, ...)                    \
  const ClusterChannel& fn() {                       \
    static const ClusterChannel c{__VA_ARGS__};      \
    return c;                                        \
  }

LIPOLGEN_CHANNEL(DEUTERON_P_TAG,
                 "d: DIS on neutron, proton spectator (control)", 2, 1, "p", 1,
                 1, 2.2246e-3, 0)
LIPOLGEN_CHANNEL(DEUTERON_N_TAG,
                 "d: DIS on proton, neutron spectator (control)", 2, 1, "n", 1,
                 0, 2.2246e-3, 0)
LIPOLGEN_CHANNEL(HE3_P_TAG,
                 "3He: p spectator (p+d two-body approx., control)", 3, 2, "p",
                 1, 1, 5.49e-3, 0)
LIPOLGEN_CHANNEL(LI6_ALPHA_TAG, "6Li: DIS on d-cluster, alpha spectator", 6, 3,
                 "alpha", 4, 2, 1.4743e-3, 0)
LIPOLGEN_CHANNEL(LI6_D_TAG, "6Li: DIS on alpha-cluster, d spectator", 6, 3, "d",
                 2, 1, 1.4743e-3, 0)
LIPOLGEN_CHANNEL(LI7_ALPHA_TAG, "7Li: DIS on t-cluster, alpha spectator", 7, 3,
                 "alpha", 4, 2, 2.4670e-3, 1)
LIPOLGEN_CHANNEL(LI7_T_TAG, "7Li: DIS on alpha-cluster, triton spectator", 7, 3,
                 "t", 3, 1, 2.4670e-3, 1)

#undef LIPOLGEN_CHANNEL

const ClusterChannel& channel_by_name(const std::string& name) {
  if (name == "DEUTERON_P_TAG") return DEUTERON_P_TAG();
  if (name == "DEUTERON_N_TAG") return DEUTERON_N_TAG();
  if (name == "HE3_P_TAG") return HE3_P_TAG();
  if (name == "LI6_ALPHA_TAG") return LI6_ALPHA_TAG();
  if (name == "LI6_D_TAG") return LI6_D_TAG();
  if (name == "LI7_ALPHA_TAG") return LI7_ALPHA_TAG();
  if (name == "LI7_T_TAG") return LI7_T_TAG();
  throw std::runtime_error("unknown cluster channel " + name);
}

double momentum_density(double k, double kappa, double beta, int l_wave) {
  const double k2 = k * k;
  double psi;
  if (l_wave == 0) {
    psi = 1.0 / (k2 + kappa * kappa) - 1.0 / (k2 + beta * beta);
  } else {
    psi = k / ((k2 + kappa * kappa) * (k2 + beta * beta));
  }
  return psi * psi;
}

MomentumSampler::MomentumSampler(const ClusterChannel& channel, double beta,
                                 double k_max, std::size_t n_grid) {
  grid_ = linspace(1e-4, k_max, n_grid);
  cdf_.resize(n_grid);
  const double kap = channel.kappa();
  double acc = 0.0;
  for (std::size_t i = 0; i < n_grid; ++i) {
    const double g = grid_[i];
    acc += g * g * momentum_density(g, kap, beta, channel.l_wave);
    cdf_[i] = acc;
  }
  const double tot = cdf_.back();
  for (double& v : cdf_) v /= tot;
}

double MomentumSampler::k_of_u(double u) const {
  // `np.interp(u, cdf, grid)`: the cdf is the abscissa, the k grid the
  // ordinate, with clamped extrapolation at both ends.
  return np_interp(u, cdf_, grid_);
}

void MomentumSampler::sample(Rng& rng, double& kx, double& ky,
                             double& kz) const {
  const double k = k_of_u(rng.uniform());
  const double cos_t = 2.0 * rng.uniform() - 1.0;
  const double phi = 2.0 * kPi * rng.uniform();
  const double sin_t = std::sqrt(std::fmax(1.0 - cos_t * cos_t, 0.0));
  kx = k * sin_t * std::cos(phi);
  ky = k * sin_t * std::sin(phi);
  kz = k * cos_t;
}

FragmentLab boost_fragment(const ClusterChannel& channel, double p_per_nucleon,
                           double kx, double ky, double kz, double m,
                           int frag_Z, int frag_A) {
  const double e_rest = std::sqrt(m * m + kx * kx + ky * ky + kz * kz);
  const double m_beam = channel.m_beam();
  const double p_beam = channel.beam_A * p_per_nucleon;
  const double e_beam = std::sqrt(p_beam * p_beam + m_beam * m_beam);
  const double gamma = e_beam / m_beam;
  const double gbeta = p_beam / m_beam;
  FragmentLab out;
  out.pz_lab = gamma * kz + gbeta * e_rest;
  out.e_lab = gamma * e_rest + gbeta * kz;
  out.pT = std::sqrt(kx * kx + ky * ky);
  out.p_lab = std::sqrt(out.pT * out.pT + out.pz_lab * out.pz_lab);
  out.theta = std::atan2(out.pT, out.pz_lab);
  out.phi = std::atan2(ky, kx);
  const double rigidity_beam = p_beam / channel.beam_Z;
  out.R = frag_Z > 0 ? (out.p_lab / frag_Z) / rigidity_beam : kNaN;
  out.xL = out.p_lab / (frag_A * p_per_nucleon);
  out.k = std::sqrt(kx * kx + ky * ky + kz * kz);
  return out;
}

FragmentLab boost_spectator_fragment(const ClusterChannel& channel,
                                     double p_per_nucleon, double kx,
                                     double ky, double kz) {
  return boost_fragment(channel, p_per_nucleon, kx, ky, kz, channel.m_spec(),
                        channel.spectator_Z, channel.spectator_A);
}

BreakupLab boost_breakup(const ClusterChannel& channel, double p_per_nucleon,
                         double kx, double ky, double kz) {
  BreakupLab out;
  out.spectator = boost_spectator_fragment(channel, p_per_nucleon, kx, ky, kz);
  out.partner = boost_fragment(channel, p_per_nucleon, -kx, -ky, -kz,
                               channel.m_partner(), channel.partner_Z(),
                               channel.partner_A());
  out.kx = kx;
  out.ky = ky;
  out.kz = kz;
  out.k = out.spectator.k;
  return out;
}

// --------------------------------------------------------------- farforward

bool Optics::isotropic() const {
  return sigma_theta_v < 0.0 || std::fabs(sigma_theta_v - sigma_theta) < 1e-15;
}

bool Optics::clears(double theta, double phi) const {
  if (std::isnan(phi)) return theta > envelope_x();
  return std::fabs(theta * std::cos(phi)) > envelope_x()
         || std::fabs(theta * std::sin(phi)) > envelope_y();
}

const Optics& HIGH_ACCEPTANCE() {
  static const Optics o{"high-acceptance", 0.20 / (10.0 * 275.0), 10.0, -1.0, 1.0};
  return o;
}
const Optics& HIGH_DIVERGENCE() {
  static const Optics o{"high-divergence", 0.45 / (10.0 * 275.0), 10.0, -1.0, 1.0};
  return o;
}

const PotLevers& pot_levers(const std::string& config) {
  //             R12     R34    D       D2      blind half width
  static const PotLevers l5{19.24, -1.0, 0.311, -0.190, 0.048};
  static const PotLevers l10{21.25, 3.35, 0.287, -0.206, 0.032};
  static const PotLevers l18{29.97, 2.93, 0.292, -0.215, 0.016};
  if (config == "5x41") return l5;
  if (config == "10x100") return l10;
  if (config == "18x275") return l18;
  throw std::runtime_error("no pot levers for configuration " + config);
}

bool over_rigid_route(double r, double theta_x, const std::string& config) {
  const PotLevers& l = pot_levers(config);
  const double delta = r - 1.0;
  double x = std::fabs(delta) <= MEASURED_DELTA_MAX
                 ? l.dispersion * delta + l.dispersion2 * delta * delta
                 : l.dispersion * delta;
  x += l.r12 * theta_x;
  const double ax = std::fabs(x);
  return ax >= l.blind_half_width && ax <= POT_OUTER_HALF_WIDTH;
}

int route_charged(double r, double theta, double pT, const Optics& optics,
                  double phi, double theta_outer,
                  const std::string& pot_config) {
  (void)pT;  // the near-beam decision is angular; pT is kept for call symmetry
  int out = kRouteLost;
  if (theta >= THETA_B0_MIN && theta <= THETA_B0_MAX) out = kRouteB0;
  const double outer = std::isnan(theta_outer) ? THETA_RP_OUTER : theta_outer;
  const bool small = theta < outer;
  const bool near = std::fabs(r - 1.0) < NEAR_BEAM_BAND;
  // Assignment order mirrors the Python's overwrite order exactly:
  // OMD, then RP, then the near-beam tail, then the over-rigid branch.
  if (small && r >= OMD_R_WINDOW_LO && r < OMD_R_WINDOW_HI) out = kRouteOMD;
  if (small && !near && r >= RP_R_WINDOW_LO && r <= RP_R_WINDOW_HI) {
    out = kRouteRomanPots;
  }
  if (small && near && optics.clears(theta, phi)) out = kRouteRPNearBeam;
  const double theta_x = std::isnan(phi) ? 0.0 : theta * std::cos(phi);
  if (small && r > 1.0 + NEAR_BEAM_BAND
      && over_rigid_route(r, theta_x, pot_config)) {
    out = kRouteRPInner;
  }
  return out;
}

int route_neutral(double theta) {
  return theta <= THETA_ZDC_MAX ? kRouteZDC : kRouteLost;
}

std::string yr_config_key(const std::string& ion_name, double p_per_nucleon) {
  const Ion& ion = ion_by_name(ion_name);
  const char* keys[3] = {"5x41", "10x100", "18x275"};
  int best = 0;
  double best_d = 1e300;
  for (int i = 0; i < 3; ++i) {
    const double d = std::fabs(ion.momentum_per_nucleon_at(PROTON_CONFIG_ENERGIES[i])
                               - p_per_nucleon);
    if (d < best_d) {
      best_d = d;
      best = i;
    }
  }
  return keys[best];
}

void sigma_theta_for(const std::string& ion_name, double p_per_nucleon,
                     bool high_acceptance, double& sigma_h, double& sigma_v) {
  // Yellow Report Table 10.1, HADRON beam: (HD_h, HD_v, HA_h, HA_v) [urad].
  const std::string key = yr_config_key(ion_name, p_per_nucleon);
  double hd_h, hd_v, ha_h, ha_v, p_e;
  if (key == "18x275") {
    hd_h = 150.0; hd_v = 150.0; ha_h = 65.0; ha_v = 65.0; p_e = 275.0;
  } else if (key == "10x100") {
    hd_h = 220.0; hd_v = 220.0; ha_h = 180.0; ha_v = 180.0; p_e = 100.0;
  } else {
    hd_h = 220.0; hd_v = 380.0; ha_h = 220.0; ha_v = 380.0; p_e = 41.0;
  }
  const double h = high_acceptance ? ha_h : hd_h;
  const double v = high_acceptance ? ha_v : hd_v;
  const double g_p = std::sqrt(p_e * p_e + PROTON_MASS * PROTON_MASS) / PROTON_MASS;
  const double bg_p = std::sqrt(g_p * g_p - 1.0);
  const double bg_i = p_per_nucleon / ion_by_name(ion_name).mass_per_nucleon();
  const double f = std::sqrt(bg_p / bg_i);
  sigma_h = 1e-6 * h * f;
  sigma_v = 1e-6 * v * f;
}

Optics yr_optics(const std::string& ion_name, double p_per_nucleon,
                 bool high_acceptance, double n_sigma) {
  double h = 0.0, v = 0.0;
  sigma_theta_for(ion_name, p_per_nucleon, high_acceptance, h, v);
  Optics o;
  o.name = yr_config_key(ion_name, p_per_nucleon)
           + (high_acceptance ? " high-acceptance" : " high-divergence");
  o.sigma_theta = h;
  o.n_sigma = n_sigma;
  o.sigma_theta_v = v;
  return o;
}

}  // namespace lipolgen
