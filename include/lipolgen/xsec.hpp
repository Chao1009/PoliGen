// SPDX-License-Identifier: GPL-3.0-or-later
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
#include <utility>

#include "lipolgen/asymmetries.hpp"
#include "lipolgen/beams.hpp"
#include "lipolgen/sf.hpp"

namespace lipolgen {

// --- finite-gamma (target-mass) kinematics, E143 PRD 58:112003 -----------
//
// `gamma_squared`, `epsilon_gamma`, `depolarization_d_gamma` and `eta_gamma`
// live in `asymmetries.hpp` since 2026-08-29 -- ONE implementation for both
// halves of the library, mirroring the move of the same block into
// `polli_fastsim.asymmetries`.  The alias below is the name this module has
// always used, exactly as `polligen.xsec` re-exports it.
inline double depolarization_gamma(double y, double gamma2, double r) {
  return depolarization_d_gamma(y, gamma2, r);
}

// --- the exact finite-gamma tensor sector (Cosyn et al., plans/08 D2) -----
//
// Cosyn, Roldan Tomei, Sosa and Zec, EPJ A 61 (2025) 83 (arXiv:2410.12764)
// decompose the tensor-polarized inclusive cross section (their Eq. 10) into
// four structure functions whose GEOMETRY is explicit -- F[U T_LL, T],
// F[U T_LL, L], F[U T_LT] and F[U T_TT] -- and give them in terms of b1..b4
// in their Eqs. (17a)-(17e).  The two functions below are those equations and
// Eq. (16) transcribed literally, plus the rest-frame angle theta_q of
// Eq. (24).
//
// CONVENTION MAP.  Cosyn writes a deuteron with its own Bjorken variable x_d
// and mass M_d, and plots against the rescaled x = 2 x_d.  This library is
// per-nucleon throughout: every structure function is the target's divided by
// A and x is the per-nucleon Bjorken variable, i.e. the target is treated as a
// spin-1 object of mass M_NUCLEON with Bjorken variable x.  The map is
// x_d -> x, M_d -> M_NUCLEON, b_i -> the per-nucleon b_i, under which Cosyn's
// gamma = 2 x_d M_d/Q is exactly this program's `gamma_squared`, his tensor
// Callan-Gross b2 = 2 x_d b1 is the default b2 = 2 x b1, and his F2 = 2 x_d F1
// is `NuclearF2`'s F2 = 2 x (1 + R) F1 at R = 0.
//
// THE ANCHOR IS THE PAPER'S OWN TABLE 1, not a re-typing of Eqs. (17).  Both
// rows of that table which retain finite gamma are closed-form functions of
// (gamma^2, eps, theta_q), and both come out of this module exactly (1e-10 in
// tests/test_xsec.cpp): the second row, b1/(F1 A_T) = -9(1 + eps gamma^2)/
// (6 + 5 gamma^2 + 2 eps gamma^4) for a target polarized along q, and the
// third row, the same quantity along the beam, which additionally exercises
// F[U T_LT] and F[U T_TT] and the geometry of Eq. (22).  Two transcription
// errors were caught by exactly that anchor and are recorded here so they are
// not made again: the leading factor 2 of Eq. (17a) multiplies b1 ALONE and
// not the whole bracket, and the incoming beam sits at +sin theta_q in the
// photon frame, which is what makes both rows come out.  Eq. (22b) as printed
// carries the opposite sign for T_LT cos phi_TL; it is the one place where the
// paper is not consistent with itself, since its own Table 1, its Fig. 5 and
// the axis triad of its Fig. 2 all require the frame this module uses.

/// (cos theta_q, sin theta_q): Cosyn Eq. (24).  theta_q is the rest-frame
/// angle between the virtual photon and the incoming electron -- zero for a
/// massless target, O(gamma) otherwise.  It is what makes a spin axis
/// transverse to the BEAM not transverse to q, and so what leaks the b-sector
/// rate term into cos 2phi.
std::pair<double, double> theta_q_cos_sin(double y, double gamma2);

/// The four tensor structure functions of Cosyn Eqs. (17a)-(17e):
/// (F_TLL_T, F_TLL_L, F_TLT, F_TTT).
///
/// At gamma = 0 this collapses to F_TLL_T = -2 b1,
/// F_TLL_L = (2 x b1 - b2)/x and F_TLT = F_TTT = 0 -- b3 and b4 cancel
/// identically -- which is exactly the b-sector of the massless
/// Hoodbhoy-Jaffe-Manohar master formula (`InclusiveKernel::tensor_kernel`
/// over `dphi`), for ANY b2 and not only at the tensor Callan-Gross point.
struct CosynTensorSFs {
  double f_t = 0.0;   ///< F[U T_LL, T]
  double f_l = 0.0;   ///< F[U T_LL, L]
  double f_lt = 0.0;  ///< F[U T_LT]
  double f_tt = 0.0;  ///< F[U T_TT]
};
CosynTensorSFs cosyn_tensor_sfs(double b1, double b2, double b3, double b4,
                                double x, double gamma2);

/// Cosyn Eq. (16): (F_UU_T, F_UU_L) = (2 F1, (1+gamma^2) F2/x - 2 F1).
/// F_UU_T + eps F_UU_L is the exact denominator of every tensor asymmetry; at
/// gamma = 0 it is y^2/(1-y+y^2/2) times the massless D_phi, and
/// F_UU_L/F_UU_T is R there.
std::pair<double, double> cosyn_unpolarized_sfs(double f1, double f2, double x,
                                                double gamma2);

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
  /// The higher-twist tensor slots of Cosyn Eqs. (17); zero unless a
  /// `b3_func` / `b4_func` was given, and read ONLY by the `tensor_gamma`
  /// path (they cancel identically at gamma = 0).
  double b3 = 0.0;
  double b4 = 0.0;
  double delta = 0.0;
  double g2 = 0.0;
  bool has_g2 = false;
};

/// The (constant, cos phi', cos 2phi') harmonics of the exact finite-gamma
/// tensor rate shift (`InclusiveKernel::tensor_harmonics_gamma`).
struct TensorHarmonics {
  double h0 = 0.0;
  double h1 = 0.0;
  double h2 = 0.0;
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
/// `target_mass` (DEFAULT TRUE since 2026-08-29, matching `xsec.py`'s own
/// default) selects the longitudinal vector kernel: true is the exact
/// finite-gamma D_gamma (A1 + eta A2), false the massless A_par = D g1/F1 of
/// the fast simulation that every number published before that date used.
///
/// `g2_mode` chooses the g2 model (Wandzura-Wilczek, the default, or zero) and
/// `g2_scale` MULTIPLIES it.  Together they are the twist-3 handle on the
/// finite-gamma A_par -- the residual systematic that replaced the target-mass
/// bias is measured by re-running an extraction at g2_scale = 0 and 1.5
/// (`evgen/scripts/target_mass_bound.py`).  `target_mass = true` with
/// `G2Mode::kZero` is legitimate and is exactly the g2_scale = 0 variation:
/// g2 is filled with zeros and the finite-gamma kernel is evaluated on them.
class InclusiveKernel {
 public:
  struct Options {
    std::shared_ptr<const UnpolSF> f2_source;   ///< default: ToyF2
    std::shared_ptr<const PolSF> g1_model;      ///< default: ToyG1 on f2_source
    SFFunc3 b1_func, b2_func, delta_func;       ///< spin-1 rank-2 slots
    /// Higher-twist tensor slots b3, b4 of Cosyn Eqs. (17); both default to
    /// empty -> zero, and only the `tensor_gamma` path reads them.
    SFFunc3 b3_func, b4_func;
    /// `tensor_gamma` (DEFAULT FALSE) selects the tensor b-sector kernel:
    /// false is the massless Hoodbhoy-Jaffe-Manohar one, bit for bit what
    /// every published number was made on, and true the exact finite-gamma
    /// Cosyn kernel of the header block above.  The two agree identically at
    /// gamma = 0 -- for any b2, with b3 and b4 cancelling.
    ///
    /// It is off by default because the O(gamma^2) leakage of the rate sector
    /// into cos 2phi is carried as a SYSTEMATIC of the Delta extraction rather
    /// than as a correction the extraction subtracts (at most 0.109 % of the
    /// published amplitude, and model-dependent through the unmeasured b3,
    /// b4), and because switching it on would silently move every published
    /// tensor number by O(gamma^2) with nothing on the analysis side to meet
    /// it.
    bool tensor_gamma = false;
    /// Spin-3/2 rank-2 slots.  They take the HJM/Cosyn-sign b1 (the same
    /// F1^(m) = F1 - Q_NN b1 relation `TENSOR_LL_SIGN` pins for spin 1), which
    /// is MINUS the b1 of arXiv:2209.12161 Eq. (19b).  There is deliberately
    /// no rank-3 (octupole) slot: its leading-twist function `g1_rank3`
    /// (2209.12161 Eq. (19d), renamed from their g2) multiplies the beam
    /// helicity only, so the rank-<=2 truncation is exact for every
    /// unpolarized-beam observable -- docs/theory/SPIN32_FINITE_GAMMA.md
    /// secs. 3-6 give the map, the theorem and the insertion points.
    SFFunc3 b1_32_func, b2_32_func, delta_32_func;
    G2Mode g2_mode = G2Mode::kWandzuraWilczek;
    /// Multiplies g2 (the WW table, or nothing when `g2_mode` is kZero).
    double g2_scale = 1.0;
    std::function<double(double)> emc_ratio;
    RFunc r_func;
    bool target_mass = true;
    int g2_npts = 96;
  };

  explicit InclusiveKernel(Ion ion);
  InclusiveKernel(Ion ion, Options options);

  const Ion& ion() const { return ion_; }
  bool target_mass() const { return target_mass_; }
  bool tensor_gamma() const { return tensor_gamma_; }
  double g2_scale() const { return g2_scale_; }

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

  /// (h0, h1, h2) of the EXACT finite-gamma b-sector: the constant, cos phi'
  /// and cos 2phi' harmonics of the tensor rate shift.
  ///
  /// The massless master formula puts the whole b1/b2 sector in the
  /// phi-independent term, because at gamma = 0 the virtual photon is along
  /// the beam and a spin axis at theta_S to the beam is at theta_S to q.  At
  /// finite gamma it is not: Cosyn Eq. (24) gives the rest-frame angle
  /// theta_q between q and the beam, so the alignment tensor of Eq. (9)
  /// acquires, in the PHOTON frame, components that depend on the
  /// lepton-plane azimuth.  Writing the spin direction in that frame with the
  /// axis at (theta_S, phi_S) to the beam and phi' = phi - phi_S,
  ///
  ///   N = (c s_S cos phi' + s c_S,  s_S sin phi',
  ///        c c_S - s s_S cos phi'),      c, s = cos, sin theta_q,
  ///
  /// (the x axis is the standard one, along the incoming lepton's transverse
  /// projection, so that the beam sits at (s, 0, c) -- the triad of Cosyn
  /// Fig. 2, and the choice under which both finite-gamma rows of Table 1
  /// come out exactly), and the three polarization parameters of Eq. (14) are
  /// T_LL = t_zz, T_LT cos phi_TL = t_xz, T_TT cos 2phi_TT = t_xx - t_yy with
  /// t_ij = (3 Q_NN/2)(N_i N_j - d_ij/3) -- Eq. (9) for a pure state,
  /// normalized on Q_NN so that it is the same geometry `tensor_moments` uses
  /// for every spin.  Each is a quadratic in cos phi', hence exactly a
  /// constant plus cos phi' plus cos 2phi'.
  ///
  /// The three channels are combined with Eq. (17) and the exact denominator
  /// of Eq. (16) as in Eq. (19) and multiplied by -TENSOR_LL_SIGN, so that the
  /// constant harmonic is the same rate shift w_avg the massless path returns
  /// (the program's w is half of Cosyn's A_T, and A_T carries the opposite
  /// sign of b1 to the program's own master formula -- which is what
  /// TENSOR_LL_SIGN is, and why the sign of this leakage was gated on D1).
  ///
  /// Two caveats, both documented rather than hidden.  (i) The cos phi'
  /// harmonic's SIGN depends on which of the two directions in the lepton
  /// plane defines phi' = 0; here it is the incoming lepton's transverse
  /// projection, the convention Eq. (14b) is written in.  The constant and
  /// the cos 2phi' harmonic are even under phi' -> phi' + pi and do not depend
  /// on it -- but they DO depend on the frame, through t_xz, which is why the
  /// Table 1 anchor and not Eq. (22b) fixes it.  (ii) Delta -- the gluon-
  /// transversity term of `a_cos2phi`, which is not part of Cosyn's b1-b4
  /// basis -- keeps its massless kinematic factor in both paths, so that
  /// switching this on does not silently redefine the structure function the
  /// whole programme extracts.
  TensorHarmonics tensor_harmonics_gamma(const SFTables& t, double x, double q2,
                                         double y,
                                         const EventSpinState& state) const;

  /// The b-sector (tensor RATE) contribution to (w_avg, a1, a2) ALONE -- the
  /// piece the tensor RC band rescales (`rc.hpp`).  Exactly the terms
  /// `amplitudes()` adds inside its `j >= 1` branch, on whichever path
  /// `tensor_gamma` selects, with the Delta (gluon-transversity) cos 2phi term
  /// included only when `with_delta` (`RcScope::TensorAll`).
  ///
  /// `amplitudes()` is DEFINED as this plus the vector terms, so the two can
  /// never drift apart and the RC band's tau = W_tensor/W is the rank-2
  /// projection of the very density the sampler drew from.  Returns all zeros
  /// below spin 1 (`tensor_moments` gives (0, 0) and `tables()` leaves
  /// b1 = b2 = Delta = 0 there).
  ///
  /// THROWS on the same spin mismatch `amplitudes()` does -- see below.
  Amplitudes tensor_amplitudes(const SFTables& t, double x, double q2, double s,
                               const EventSpinState& state,
                               bool with_delta = false) const;

  /// (w_avg, a_1, a_2).  The VECTOR a_1 is only computed when `with_perp` (it
  /// needs g2); the tensor sector contributes to a_1 as well, but only on the
  /// exact finite-gamma path (`tensor_gamma`, off by default) and only where
  /// the axis is neither along the beam nor transverse to it, since every
  /// cos phi' coefficient there carries sin theta_S cos theta_S.
  ///
  /// THROWS if `state.j` is not the kernel's own ion spin: the rank-2 branch
  /// is gated on `state.j` but `tensor_moments` reads `ion().spin`, so a
  /// mismatched pair would silently mix the two (P8).  The check lives in
  /// `tensor_amplitudes`, which this calls first and unconditionally.
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
  SFFunc3 b3_func_, b4_func_;
  G2Mode g2_mode_;
  double g2_scale_;
  bool target_mass_;
  bool tensor_gamma_;
  int g2_npts_;
};

}  // namespace lipolgen

#endif  // LIPOLGEN_XSEC_HPP
