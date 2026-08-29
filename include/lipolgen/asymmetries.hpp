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
/// Conventions:
///   Azz  = (s+ + s- - 2 s0)/(s+ + s- + s0)   (longitudinal, theta_m = 0)
///        = (2/3)[b1 + (1-y)/(x y^2) b2]/[F1 + (1-y)/(x y^2) F2]
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

/// Longitudinal double-spin asymmetry in the A1 ~= g1/F1 approximation:
/// A_par = D(y) g1/F1.  This is the gamma -> 0 limit of the exact E143 form;
/// the finite-gamma one lives behind `InclusiveKernel(target_mass=true)`.
double a_parallel(double g1, double f1, double y, double x, double q2,
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
