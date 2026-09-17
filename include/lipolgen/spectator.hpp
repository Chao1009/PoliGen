// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LIPOLGEN_SPECTATOR_HPP
#define LIPOLGEN_SPECTATOR_HPP

/// \file spectator.hpp
/// Cluster-spectator kinematics and far-forward routing.  C++17 port of
/// `fastsim/polli_fastsim/spectator.py` plus the part of
/// `fastsim/polli_fastsim/farforward.py` that `tagged.rp_accepted` and
/// `boost_spectator` need (Optics envelope, `route_charged`, the over-rigid
/// branch, the per-configuration Yellow Report optics).
///
/// Model: the lithium ground state is a two-cluster system (6Li = alpha + d,
/// S-wave dominant; 7Li = alpha + t, P-wave).  DIS strikes one cluster, the
/// other is a spectator carrying the internal relative momentum k.
///
///   S-wave  psi(k) ~ 1/(k^2+kappa^2) - 1/(k^2+beta^2)          (Hulthen)
///   P-wave  psi(k) ~ k / ((k^2+kappa^2)(k^2+beta^2))
///
/// kappa = sqrt(2 mu S) with mu the reduced mass of the two FREE (separated)
/// clusters -- so kappa is a property of the BEAM, identical for the alpha-tag
/// and the partner-tag of one nucleus.  Masses are the physical AME2020-derived
/// NUCLEAR masses, never A * M_U: the beam mass sets gamma and gamma*beta and
/// with them every lab quantity, and the visible consequence is the rigidity,
///
///     R(k = 0) = (m_spec / Z_spec) / (m_beam / Z_beam),
///
/// a ratio of mass-to-charge ratios -- 0.99813 for the 6Li alpha (not 1) and
/// 0.85571 for the 7Li alpha (not 6/7).

#include <limits>
#include <string>
#include <vector>

#include "lipolgen/rng.hpp"

namespace lipolgen {

/// Atomic mass unit [GeV], as `spectator.M_U`.
inline constexpr double M_U = 0.93149;

/// Cluster masses [GeV] keyed by the `ClusterChannel::spectator` name, the
/// same table rounded to 10 keV that `spectator.MASSES` carries verbatim
/// (it sets the SPECTATOR mass in every published distribution).
/// Throws std::runtime_error for an unknown name.
double cluster_mass(const std::string& name);

/// Ground-state NUCLEAR (not atomic) mass [GeV] of the nuclide (Z, A), from
/// the AME2020-derived table of `spectator.NUCLEUS_MASS`.  Falls back to
/// A * M_U for a nuclide the table does not carry -- low for everything in
/// the table, HIGH at (Z, A) = (6, 12) where u is defined.
double nuclear_mass(int z, int a);
/// True when (z, a) is in the table rather than reached by the A*M_U fallback.
bool nuclear_mass_known(int z, int a);

/// Beam nucleus -> struck cluster + spectator cluster.
struct ClusterChannel {
  std::string name;
  int beam_A = 1;
  int beam_Z = 1;
  std::string spectator;   ///< key into `cluster_mass`
  int spectator_A = 1;
  int spectator_Z = 1;
  double separation_energy = 0.0;  ///< [GeV]
  int l_wave = 0;                  ///< 0 (S) or 1 (P) relative motion

  double m_spec() const { return cluster_mass(spectator); }
  /// Ground-state nuclear mass of the beam [GeV].
  double m_beam() const { return nuclear_mass(beam_Z, beam_A); }
  int partner_A() const { return beam_A - spectator_A; }
  int partner_Z() const { return beam_Z - spectator_Z; }
  /// Mass of the FREE struck cluster [GeV]: the table when it carries the
  /// nuclide, else m_beam - m_spec + S, which keeps the DEFINITION
  /// S = m_spec + m_partner - m_beam exact off the table too.
  double m_partner() const;
  /// Bound-state momentum scale sqrt(2 mu S) [GeV], mu the reduced mass of
  /// the two FREE clusters.
  double kappa() const;
  /// The closed form (m_spec/Z_spec)/(m_beam/Z_beam); NaN for a neutral
  /// spectator.
  double r_at_k_zero() const;
};

/// The seven channels of `spectator.py`.
const ClusterChannel& DEUTERON_P_TAG();
const ClusterChannel& DEUTERON_N_TAG();
const ClusterChannel& HE3_P_TAG();
const ClusterChannel& LI6_ALPHA_TAG();
const ClusterChannel& LI6_D_TAG();
const ClusterChannel& LI7_ALPHA_TAG();
const ClusterChannel& LI7_T_TAG();
/// Lookup by the Python module-level name ("LI6_ALPHA_TAG", ...).
const ClusterChannel& channel_by_name(const std::string& name);

/// Default short-range scale of the momentum densities [GeV]; the model band
/// is 0.20 - 0.40 (`beta_band_lo/hi`), which is the dominant tail systematic.
inline constexpr double BETA_DEFAULT = 0.30;
inline constexpr double BETA_BAND_LO = 0.20;
inline constexpr double BETA_BAND_HI = 0.40;

/// Unnormalized n(k) = |psi(k)|^2, NOT including the k^2 phase space.
double momentum_density(double k, double kappa, double beta, int l_wave);

/// Inverse-CDF sampler of k^2 n(k) on [1e-4, k_max] over a 30000-point grid,
/// exactly the grid `spectator.sample_k` builds (the grid is the model, so it
/// is reproduced even though the random stream cannot be).
class MomentumSampler {
 public:
  MomentumSampler(const ClusterChannel& channel, double beta = BETA_DEFAULT,
                  double k_max = 1.5, std::size_t n_grid = 30000);
  /// |k| from the inverse CDF at a uniform deviate u in (0, 1).
  double k_of_u(double u) const;
  /// One isotropic (kx, ky, kz) draw.
  void sample(Rng& rng, double& kx, double& ky, double& kz) const;
  const std::vector<double>& grid() const { return grid_; }
  const std::vector<double>& cdf() const { return cdf_; }

 private:
  std::vector<double> grid_, cdf_;
};

/// Lab kinematics of one fragment of a two-body breakup.
struct FragmentLab {
  double pT = 0.0;
  double theta = 0.0;   ///< polar angle from the ion axis [rad]
  double phi = 0.0;     ///< lab azimuth, atan2(ky, kx)
  double p_lab = 0.0;
  double R = 0.0;       ///< rigidity ratio vs the beam; NaN for Z = 0
  double xL = 0.0;
  double k = 0.0;       ///< |k| in the beam rest frame
  double e_lab = 0.0;   ///< lab energy [GeV] (not a `spectator.py` key)
  double pz_lab = 0.0;  ///< lab longitudinal momentum [GeV]
};

/// `spectator._boost_fragment`: rest-frame (kx, ky, kz) -> lab, for a fragment
/// of mass `m`, charge `frag_Z` and mass number `frag_A`.
FragmentLab boost_fragment(const ClusterChannel& channel, double p_per_nucleon,
                           double kx, double ky, double kz, double m,
                           int frag_Z, int frag_A);
/// The same for the channel's own spectator (mass, Z, A taken from it).
FragmentLab boost_spectator_fragment(const ClusterChannel& channel,
                                     double p_per_nucleon, double kx,
                                     double ky, double kz);
/// Both fragments of ONE breakup: spectator at +k, partner at -k.
struct BreakupLab {
  FragmentLab spectator, partner;
  double kx = 0.0, ky = 0.0, kz = 0.0, k = 0.0;
};
BreakupLab boost_breakup(const ClusterChannel& channel, double p_per_nucleon,
                         double kx, double ky, double kz);

// --------------------------------------------------------------- farforward
//
// Windows fetch-verified 2026-06-12 (YR detector matrix; arXiv:2108.08314
// Table I).  Only the ANGULAR windows enter the routing.

inline constexpr double THETA_RP_MAX = 5.0e-3;
inline constexpr double THETA_RP_OUTER = THETA_RP_MAX;
inline constexpr double THETA_B0_MIN = 5.5e-3;
inline constexpr double THETA_B0_MAX = 20.0e-3;
inline constexpr double THETA_ZDC_MAX = 4.0e-3;
inline constexpr double RP_R_WINDOW_LO = 0.60;
inline constexpr double RP_R_WINDOW_HI = 0.95;
inline constexpr double OMD_R_WINDOW_LO = 0.45;
inline constexpr double OMD_R_WINDOW_HI = 0.65;
/// |R - 1| below this: inside the beam envelope.
inline constexpr double NEAR_BEAM_BAND = 0.05;
/// The momentum the published near-beam pT cuts are quoted at [GeV].
inline constexpr double PROTON_REFERENCE_MOMENTUM = 275.0;

/// Near-beam envelope of an optics setting, as the ANGLE it subtends.  The
/// Roman Pots are planar, so a fragment at beam rigidity must clear a
/// RECTANGLE of half-widths (n_sigma sigma_h, n_sigma sigma_v): it is tagged
/// when |theta_x| > n sigma_h OR |theta_y| > n sigma_v.  With no azimuth the
/// cut degenerates to the inscribed circle at n sigma_h -- the more generous
/// of the two (1.7x at the tagging optics) and used by nothing published.
struct Optics {
  std::string name;
  double sigma_theta = 0.0;      ///< rad, horizontal
  double n_sigma = 10.0;
  double sigma_theta_v = -1.0;   ///< rad, vertical; < 0 = isotropic
  double lumi_fraction = 1.0;

  double sigma_v() const { return sigma_theta_v < 0.0 ? sigma_theta : sigma_theta_v; }
  bool isotropic() const;
  double envelope_x() const { return n_sigma * sigma_theta; }
  double envelope_y() const { return n_sigma * sigma_v(); }
  /// `phi` NaN = azimuth unknown -> the inscribed circle at n sigma_h.
  bool clears(double theta, double phi) const;
  double pt_cut_for(double momentum) const { return n_sigma * sigma_theta * momentum; }
  double pt_cut_near_beam() const { return pt_cut_for(PROTON_REFERENCE_MOMENTUM); }
};

/// The two proton-derived LEGACY envelopes (0.20 / 0.45 GeV at 275 GeV).
const Optics& HIGH_ACCEPTANCE();
const Optics& HIGH_DIVERGENCE();

/// Far-forward route codes.
enum Route : int {
  kRouteLost = 0,
  kRouteRomanPots = 1,
  kRouteOMD = 2,
  kRouteB0 = 3,
  kRouteRPNearBeam = 4,   ///< R ~ 1, accepted only outside the angular envelope
  kRouteZDC = 5,          ///< neutrals only
  kRouteRPInner = 6       ///< over-rigid: the inner side of the bend
};
/// The Roman-Pot mask `tagged.rp_accepted` applies: main window + near-beam tail.
inline bool rp_accepted(int route) {
  return route == kRouteRomanPots || route == kRouteRPNearBeam;
}

/// Pot transport levers (R12 [m], R34 [m], dispersion D [m]) and the
/// second-order dispersion D2 [m], measured 2026-08-28 per machine
/// configuration ("5x41", "10x100", "18x275"); the 5x41 vertical lever was
/// measured on 2026-08-29.  A negative `r34` would mean "never measured"; no
/// configuration the library carries is in that state any more.
struct PotLevers {
  double r12 = 0.0, r34 = -1.0, dispersion = 0.0, dispersion2 = 0.0;
  double blind_half_width = 0.0;
};
const PotLevers& pot_levers(const std::string& config);
inline constexpr double POT_OUTER_HALF_WIDTH = 0.144;
/// |R - 1| over which the quadratic dispersion was fitted and may be used.
inline constexpr double MEASURED_DELTA_MAX = 0.30;

/// Does an over-rigid fragment (R > 1) land on Roman-pot silicon?  The test is
/// on the pot-plane displacement x = D delta + D2 delta^2 + R12 theta_x with
/// delta = R - 1: silicon between the configuration's blind half-width
/// (48 / 32 / 16 mm) and the last module at 144 mm.  Outside the measured
/// |delta| <= 0.30 the quadratic is dropped (it turns over at delta ~ 0.6,
/// a three-point artefact), which is what keeps the 6Li 3He+t triton at
/// R = 1.5044 "lost".
bool over_rigid_route(double r, double theta_x = 0.0,
                      const std::string& config = "18x275");

/// The "not given" sentinel of `route_charged`'s two optional angles.  It is
/// `quiet_NaN()` and NOT the literal `0.0 / 0.0` these defaults were spelled
/// with until 2026-09-16: that is a division by zero evaluated at every
/// defaulted call, and it is exactly the expression `-ffast-math` /
/// `-ffinite-math-only` turn into an unspecified value -- at which point
/// `std::isnan(phi)` and `std::isnan(theta_outer)` inside `route_charged`
/// stop firing and the near-beam envelope silently changes shape.  Latent
/// today (CMakeLists.txt builds -O2 with no fast-math), and every other NaN
/// in the core is already `std::nan("")` or `quiet_NaN()`.
inline constexpr double kAngleUnknown =
    std::numeric_limits<double>::quiet_NaN();

/// Classify a charged fragment.  `phi` NaN = azimuth unknown, `theta_outer`
/// NaN = THETA_RP_OUTER.  `pT` is not used for the near-beam decision: the
/// envelope is angular.  It is kept in the signature so the call reads like
/// the Python's and so a caller cannot silently pass the wrong triple.
int route_charged(double r, double theta, double pT, const Optics& optics,
                  double phi = kAngleUnknown,
                  double theta_outer = kAngleUnknown,
                  const std::string& pot_config = "18x275");
/// Neutrals: ZDC inside THETA_ZDC_MAX, else lost.
int route_neutral(double theta);

/// Which Yellow Report configuration a per-nucleon ion momentum belongs to:
/// the proton energy whose gamma-matched (or rigidity-capped) per-nucleon
/// momentum reproduces it.  `ion_name` selects the species.
std::string yr_config_key(const std::string& ion_name, double p_per_nucleon);

/// (sigma_theta_h, sigma_theta_v) [rad] of a beam configuration: the Yellow
/// Report PROTON divergence of that machine configuration times the species
/// step sqrt(beta*gamma_p / beta*gamma_ion) -- 1.00 where the ion is
/// gamma-matched, sqrt(2) where the ring rigidity caps it.
void sigma_theta_for(const std::string& ion_name, double p_per_nucleon,
                     bool high_acceptance, double& sigma_h, double& sigma_v);
/// The same as an `Optics` with the 10 sigma envelope.
Optics yr_optics(const std::string& ion_name, double p_per_nucleon,
                 bool high_acceptance = true, double n_sigma = 10.0);

}  // namespace lipolgen

#endif  // LIPOLGEN_SPECTATOR_HPP
