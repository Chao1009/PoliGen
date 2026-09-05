// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LIPOLGEN_ASYMMETRIES_HPP
#define LIPOLGEN_ASYMMETRIES_HPP

/// \file asymmetries.hpp
/// Asymmetry formulas for spin-1/2 and spin-1 (tensor) targets.  Port of
/// `fastsim/polli_fastsim/asymmetries.py`.
///
/// Spin-1 master formula (unpolarized e, target spin at angle theta_m,
/// projection lambda_m in {+1, 0, -1}; Hoodbhoy-Jaffe-Manohar NPB 312:571):
///
///   dsigma/(dx dy dphi)(lambda_m) = (2 y alpha^2 / Q^2) *
///       [ F1 + (2/3) a_m b1 + (1-y)/(x y^2) (F2 + (2/3) a_m b2)
///         - (1-y)/y^2 c_m sin^2(theta_m) Delta(x,Q2) cos(2 phi) ]
///
///   a_m = (1/4) c_m (3 cos^2 theta_m - 1),  c_m = 3|lambda_m| - 2 -> (1,-2,1)
///
/// The b1, b2 terms carry the overall sign `TENSOR_LL_SIGN` (constants.hpp),
/// which since 2026-08-29 is the LITERATURE one (Cosyn et al. Eq. 27, HERMES)
/// and not the transcription above: the master formula as written has the
/// opposite sign of b1 to the b1 every published number is quoted in.
///
/// Conventions:
///   Azz  = (s+ + s- - 2 s0)/(s+ + s- + s0)   (longitudinal, theta_m = 0)
///        = -(2/3)[b1 + (1-y)/(x y^2) b2]/[F1 + (1-y)/(x y^2) F2]
///        = -(2/3) (b1/F1)/(1 + eps(y) R)  exactly -> measures b1
///   A2phi(lambda = +-1, theta_m = 90 deg) = -(1-y)/y^2 Delta / D_phi
///
/// The overall sign of the tensor RATE sector is `TENSOR_LL_SIGN`
/// (constants.hpp) -- one constant, defined once, test-guarded.

#include "lipolgen/constants.hpp"
#include "lipolgen/sf.hpp"

namespace lipolgen {

/// Virtual-photon depolarization factor D for A_par ~= D * A1.
/// A null `r_func` means `r_sigma_lt`, bit-for-bit the published value.
double depolarization_d(double y, double x, double q2,
                        const RFunc& r_func = nullptr);

// --- finite-gamma (target-mass) kinematics, E143 PRD 58:112003 -----------
//
// ONE implementation, used by both halves of the library: `xsec.hpp`'s
// `InclusiveKernel` calls these rather than carrying its own copies, so the
// generator kernel and the analytic formulas cannot drift apart in the
// O(gamma^2) term.  They are the exact lab-frame factors written in (x, y):
//   eps = 1/[1 + 2(1 + nu^2/Q^2) tan^2(theta/2)]
//   D   = (1 - E' eps/E)/(1 + eps R)
//   eta = eps sqrt(Q^2)/(E - E' eps).
// At gamma -> 0 they collapse to the massless set: eps -> (1-y)/(1-y+y^2/2)
// and D -> `depolarization_d` above.  (Port of the block that lives in
// `polli_fastsim.asymmetries` since 2026-08-29; `polligen.xsec` imports it
// from there, and `xsec.hpp` re-exports `depolarization_gamma` for the name
// this library has always used.)

/// Target-mass parameter gamma^2 = 4 M^2 x^2/Q^2 (= Q^2/nu^2).  M is the FREE
/// nucleon mass `M_NUCLEON` because x is per-nucleon; the bound-nucleon mass
/// (`Ion::mass_per_nucleon`, 0.9336 GeV for 6Li) would move gamma^2 by 1.0 %,
/// i.e. 1 % of a <= 10 % correction.
double gamma_squared(double x, double q2, double m = M_NUCLEON);
/// Virtual-photon transverse polarization at finite gamma.
double epsilon_gamma(double y, double gamma2);
/// D_gamma = [1 - (1-y) eps]/(1 + eps R) with the finite-gamma eps.
double depolarization_d_gamma(double y, double gamma2, double r);
/// eta = eps gamma y/[1 - (1-y) eps], the A2 admixture in A_par.
double eta_gamma(double y, double gamma2);

/// Longitudinal double-spin asymmetry at finite gamma (E143):
///
///   A_par = D_gamma (A1 + eta A2),
///   A1 = (g1 - gamma^2 g2)/F1,   A2 = gamma (g1 + g2)/F1.
///
/// This is the DEFAULT longitudinal kernel of the programme since 2026-08-29:
/// it is exact given g2, costs one g2^WW table per grid and nothing per call,
/// and removes the O(gamma^2) bias the massless form left on every extracted
/// g1/F1.  Both O(gamma^2) pieces are kept: eta A2 alone is about half the
/// correction and above x ~ 0.5 the half that vanishes, because g2^WW -> -g1
/// there kills A2 while the -gamma^2 g2/F1 inside A1 survives.
double a_parallel_exact(double g1, double g2, double f1, double y, double x,
                        double q2, const RFunc& r_func = nullptr);

/// Longitudinal double-spin asymmetry in the A1 ~= g1/F1 approximation:
/// A_par = D(y) g1/F1.  This is the gamma -> 0 limit of `a_parallel_exact`,
/// bit-for-bit what every figure published before 2026-08-29 was made on and
/// what `InclusiveKernel(target_mass = false)` still computes.
double a_parallel(double g1, double f1, double y, double x, double q2,
                  const RFunc& r_func = nullptr);
/// The same, with a g2 supplied: exactly `a_parallel_exact` (the Python's
/// `a_parallel(..., g2=...)` overload).
double a_parallel(double g1, double f1, double y, double x, double q2,
                  const RFunc& r_func, double g2);

/// The factor D_eff with A_par = D_eff * (g1/F1), and so the divisor that
/// turns delta(A_par) into delta(g1/F1).  Writing rho = g2/g1,
///
///   A_par = D_gamma [1 - gamma^2 rho + eta gamma (1 + rho)] (g1/F1),
///
/// so an extraction that divides a measured A_par by D_eff returns g1/F1 with
/// no O(gamma^2) bias, where dividing by the massless D returned it high by
/// (1 + gamma^2) + O(gamma^2 y).  The step is legitimate BECAUSE rho is a
/// property of the SHAPE of g1 and not of its normalization (g2^WW is linear
/// in g1), so it is the same for the model and for the measurement.  Pass
/// `has_rho = false` to get `depolarization_d` back bit-for-bit.
double depolarization_effective(double y, double x, double q2,
                                double g2_over_g1, bool has_rho = true,
                                const RFunc& r_func = nullptr);

/// gamma-suppressed transverse-vector amplitude d(y)*(A2 - xi*A1).  Both
/// O(gamma) pieces are kept: in massless kinematics gamma - xi = gamma*y/2,
/// so this reduces to d(y) * gamma * ((y/2) g1 + g2)/F1 (E143 conventions).
double a_perp(double g1, double g2, double f1, double y, double x, double q2,
              const RFunc& r_func = nullptr);

/// The phi-averaged bracket D_phi = F1 + (1-y)/(x y^2) F2.
double phi_averaged_density(double f1, double f2, double x, double y);

/// Tensor asymmetry from the master formula.  Pass `b2 = nullptr` for the
/// default b2 = 2x b1.
double azz(double b1, double f1, double f2, double x, double y,
           const double* b2 = nullptr, double theta_m = 0.0);

/// cos(2phi) amplitude for lambda_m = +-1, theta_m = 90 deg.
double a_cos2phi(double delta, double f1, double f2, double x, double y);

// --- statistical uncertainties (per kinematic bin with N events total) ---

/// delta(A_par): two-state +/- flips, equal luminosity halves.
double err_a_parallel(double n, double pe, double pz);
/// delta(Azz) for the (n+ + n- - 2 n0)/(n+ + n- + n0) estimator, equal
/// luminosity thirds: Var(num) ~= 2N -> delta = sqrt(2/N)/Pzz.
double err_azz(double n, double pzz);
/// delta(amplitude) of a cos(2phi) fit: sqrt(2/N), scaled by the tensor
/// polarization of the transverse spin states.
double err_cos2phi_amplitude(double n, double pzz);

}  // namespace lipolgen

#endif  // LIPOLGEN_ASYMMETRIES_HPP
