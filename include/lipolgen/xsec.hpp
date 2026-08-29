#ifndef LIPOLGEN_XSEC_HPP
#define LIPOLGEN_XSEC_HPP

/// \file xsec.hpp
/// Inclusive doubly polarized master cross section (the polligen "wheel").
/// Port of `evgen/polligen/xsec.py`.
///
/// Per-nucleon differential cross section for a polarized electron (helicity
/// lam_e = +-1, magnitude P_e) on a spin-J ion in the definite projection
/// state m along the axis n(theta_S, phi_S):
///
///   dsigma/(dx dQ2 dphi) = sigma_unpol(x,Q2)/(2 pi) * W(phi')
///   W = 1 + w_avg + a_1 cos(phi') + a_2 cos(2 phi'),   phi' = phi - phi_S
///
/// with D_phi = F1 + (1-y)/(x y^2) F2 and K = b1 + (1-y)/(x y^2) b2:
///
///   w_avg = t_geo(m, theta_S) * K / D_phi                    [tensor, J>=1]
///         + lam_e P_e (m/J) cos(theta_S) * A_par(x,y)        [vector-L]
///   a_1   = lam_e P_e (m/J) sin(theta_S) * A_perp(x,y)       [vector-T, gT]
///   a_2   = -(1-y)/y^2 * c_eff(m) sin^2(theta_S) * Delta/D_phi [gluonometry]
///
/// RANK-2 GEOMETRY.  Cosyn et al. (arXiv:2410.12764) Eq. (9) writes the
/// alignment tensor of any spin-J state as t_ij = (Q_NN/2)(3 n_i n_j - d_ij)
/// with the single scalar Q_NN(m) = [3 m^2 - J(J+1)]/3, so
///   t_geo = Q_NN P_2(cos theta_S)   and   c_eff = 3 Q_NN
/// is ONE geometry for both spins.  For J = 1, Q_NN = c_m/3 with
/// c_m = 3m^2 - 2, which reproduces the Hoodbhoy-Jaffe-Manohar transcription
/// digit for digit.
///
/// The b1/b2 (tensor RATE) sign is the single constant `TENSOR_LL_SIGN`; the
/// Delta sector does not depend on it.

#include <memory>

#include "lipolgen/asymmetries.hpp"
#include "lipolgen/beams.hpp"
#include "lipolgen/sf.hpp"

namespace lipolgen {

// --- finite-gamma (target-mass) kinematics, E143 PRD 58:112003 -----------
//
// The three factors below are the exact lab-frame ones written in (x, y):
//   eps = 1/[1 + 2(1 + nu^2/Q^2) tan^2(theta/2)]
//   D   = (1 - E' eps/E)/(1 + eps R)
//   eta = eps sqrt(Q^2)/(E - E' eps)
// At gamma -> 0 they collapse to the massless set the fast simulation uses.

/// Target-mass parameter gamma^2 = 4 M^2 x^2/Q^2 (= Q^2/nu^2).  M is the FREE
/// nucleon mass because x is per-nucleon.
double gamma_squared(double x, double q2, double m = M_NUCLEON);
/// Virtual-photon transverse polarization at finite gamma.
double epsilon_gamma(double y, double gamma2);
/// D_gamma = [1 - (1-y) eps]/(1 + eps R) with the finite-gamma eps.
double depolarization_gamma(double y, double gamma2, double r);
/// eta = eps gamma y/[1 - (1-y) eps], the A2 admixture in A_par.
double eta_gamma(double y, double gamma2);

/// Spin configuration of one bunch crossing category / event.
struct EventSpinState {
  int lam_e = 0;        ///< electron helicity sign (+1/-1); 0 = unpolarized
  double pe = 0.0;      ///< electron polarization magnitude in [0, 1]
  double j = 1.0;       ///< ion spin (1 or 3/2; 1/2 supported vector-only)
  double m = 0.0;       ///< projection along the quantization axis
  double theta_s = 0.0; ///< axis polar angle in the lab [rad]
  double phi_s = 0.0;   ///< axis azimuth in the lab [rad]
};

/// Per-nucleon structure functions at one (x, Q2).
struct SFTables {
  double f1 = 0.0;
  double f2 = 0.0;
  double g1 = 0.0;
  double b1 = 0.0;
  double b2 = 0.0;
  double delta = 0.0;
  double g2 = 0.0;
  bool has_g2 = false;
};

/// Modulation amplitudes of W = 1 + w_avg + a_1 cos phi' + a_2 cos 2 phi'.
struct Amplitudes {
  double w_avg = 0.0;
  double a1 = 0.0;
  double a2 = 0.0;
};

enum class G2Mode { kWandzuraWilczek, kZero };

/// Minimum over phi of W/(1 + w_avg) = 1 + A cos phi + B cos 2 phi.
///
/// With c = cos phi in [-1, 1] this is the quadratic
/// f(c) = 2B c^2 + A c + (1 - B), so the minimum is at the vertex
/// c* = -A/(4B) when B > 0 and |c*| <= 1, and at an endpoint otherwise --
/// EXACT, where the accept-reject envelope 1 + |A| + |B| is only a bound.
double density_min(double a1n, double a2n);

/// Doubly polarized inclusive e+A cross section on the lipolgen SF backends.
///
/// All structure functions are per-nucleon (F2A/A).  b1/b2/Delta enter as
/// (x, q2, f1) callables; b2 defaults to 2 x b1.  The spin-3/2 rank-2 slots
/// default to empty -> zero.
///
/// `r_func` is threaded to ALL FOUR places the kernel needs R -- F1 in
/// `NuclearF2::f1a`, F_L in `dsigma_dx_dq2`, D(y) in `depolarization_d`, and
/// the g1/F1 of the default `ToyG1` -- so one argument moves R consistently.
/// An EXPLICIT `g1_model` keeps whatever R it was built with.
///
/// `target_mass` (default false) selects the longitudinal vector kernel:
/// false is the massless A_par = D g1/F1 of the fast simulation, which every
/// published number uses; true is the exact finite-gamma D_gamma (A1 + eta A2)
/// and needs g2, so it is refused with G2Mode::kZero.
class InclusiveKernel {
 public:
  struct Options {
    std::shared_ptr<const UnpolSF> f2_source;   ///< default: ToyF2
    std::shared_ptr<const PolSF> g1_model;      ///< default: ToyG1 on f2_source
    SFFunc3 b1_func, b2_func, delta_func;       ///< spin-1 rank-2 slots
    SFFunc3 b1_32_func, b2_32_func, delta_32_func;  ///< spin-3/2 rank-2 slots
    G2Mode g2_mode = G2Mode::kWandzuraWilczek;
    std::function<double(double)> emc_ratio;
    RFunc r_func;
    bool target_mass = false;
    int g2_npts = 96;
  };

  explicit InclusiveKernel(Ion ion);
  InclusiveKernel(Ion ion, Options options);

  const Ion& ion() const { return ion_; }
  bool target_mass() const { return target_mass_; }

  /// Per-nucleon SF table at one (x, Q2).  g2 is filled when `with_g2` (a_perp
  /// needs it) or whenever `target_mass` is on.
  SFTables tables(double x, double q2, bool with_g2 = false) const;

  /// The phi-averaged bracket D_phi = F1 + (1-y)/(x y^2) F2.
  static double dphi(const SFTables& t, double x, double y);
  /// K = b1 + (1-y)/(x y^2) b2.
  static double tensor_kernel(const SFTables& t, double x, double y);

  /// Longitudinal double-spin asymmetry A_par(x, y).
  double a_parallel(const SFTables& t, double x, double q2, double y) const;
  /// gamma-suppressed transverse-vector amplitude.
  double a_perp(const SFTables& t, double x, double q2, double y) const;

  /// (Q_NN, 3 Q_NN): the rank-2 alignment of a pure state m, in ONE form for
  /// every spin.  Returns (0, 0) below spin 1.
  std::pair<double, double> tensor_moments(double m) const;

  /// (w_avg, a_1, a_2).  a_1 is only computed when `with_perp` (needs g2).
  Amplitudes amplitudes(const SFTables& t, double x, double q2, double s,
                        const EventSpinState& state,
                        bool with_perp = false) const;

  /// Unpolarized per-nucleon d2sigma/dxdQ2 [pb/GeV^2].
  double dsigma_unpol(double x, double q2, double s) const;

  /// Doubly polarized d3sigma/dxdQ2dphi [pb/GeV^2/rad].
  double dsigma(double x, double q2, double phi, double s,
                const EventSpinState& state, bool with_perp = false) const;

  /// The azimuthal density W(phi') itself, unclipped.
  static double density(const Amplitudes& a, double phip);

  /// Exact minimum of W over phi, normalized by (1 + w_avg): negative means
  /// the accept-reject would silently sample max(W, 0).
  static double positivity_margin(const Amplitudes& a);

  const NuclearF2& nuclear_f2() const { return nf2_; }
  const PolSF& g1_model() const { return *g1_model_; }

 private:
  double g1a(double x, double q2) const;

  Ion ion_;
  RFunc r_func_;
  NuclearF2 nf2_;
  std::shared_ptr<const PolSF> g1_model_;
  SFFunc3 b1_func_, b2_func_, delta_func_;
  SFFunc3 b1_32_func_, b2_32_func_, delta_32_func_;
  G2Mode g2_mode_;
  bool target_mass_;
  int g2_npts_;
};

}  // namespace lipolgen

#endif  // LIPOLGEN_XSEC_HPP
