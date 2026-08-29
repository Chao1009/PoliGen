#ifndef LIPOLGEN_SPIN_HPP
#define LIPOLGEN_SPIN_HPP

/// \file spin.hpp
/// Spin-density matrices and multipole moments for spin-1 and spin-3/2 ions.
/// C++17 port of `evgen/polligen/spin.py`.
///
/// The beam is an incoherent mixture of |J m> states along a quantization axis
/// n(theta_S, phi_S): populations p_m with sum(p_m) = 1, ordered m = +J ... -J
/// throughout the library.  The lab-frame density matrix is
///
///     rho_lab = U rho_diag U^dagger,
///     U_{m'm} = D^J_{m'm}(phi_S, theta_S, 0) = exp(-i m' phi_S) d^J_{m'm}(theta_S)
///
/// Normalized moments (dimensionless):
///
///   vector           P    = <J_z>/J                     any J, in [-1, 1]
///   tensor   J=1     Pzz  = <3 J_z^2 - 2> = p+ + p- - 2 p0,  in [-2, +1]
///   tensor   J=3/2   T    = <3 J_z^2 - J(J+1)>/3
///   octupole J=3/2   O    = <J_z^3 - (41/20) J_z>/(3/10)
///
/// The J=1 conventions match Hoodbhoy-Jaffe-Manohar NPB 312:571 (c_m = 3m^2-2).
/// Wigner-d / Clebsch-Gordan follow the standard (Varshalovich) conventions.

#include <array>
#include <complex>
#include <cstddef>
#include <vector>

namespace lipolgen {

/// Dense square matrix, row-major, minimal by design.
template <typename T>
class SquareMatrix {
 public:
  SquareMatrix() : n_(0) {}
  explicit SquareMatrix(std::size_t n) : n_(n), a_(n * n, T{}) {}
  std::size_t size() const { return n_; }
  T& operator()(std::size_t i, std::size_t j) { return a_[i * n_ + j]; }
  const T& operator()(std::size_t i, std::size_t j) const { return a_[i * n_ + j]; }
  const std::vector<T>& data() const { return a_; }

 private:
  std::size_t n_;
  std::vector<T> a_;
};

using RealMatrix = SquareMatrix<double>;
using CplxMatrix = SquareMatrix<std::complex<double>>;

/// Matrix product.
CplxMatrix matmul(const CplxMatrix& a, const CplxMatrix& b);
/// Conjugate transpose.
CplxMatrix adjoint(const CplxMatrix& a);
/// Trace.
std::complex<double> trace(const CplxMatrix& a);

/// Spin projections ordered +J ... -J (library-wide convention).
std::vector<double> m_values(double j);

/// <j1 m1 j2 m2 | j m> via the Racah formula (exact for small spins).
double clebsch_gordan(double j1, double m1, double j2, double m2,
                      double j, double m);

/// Small Wigner matrix d^j_{m'm}(beta), indexed [m', m] with the +J..-J order.
RealMatrix wigner_d(double j, double beta);

/// U_{m'm} = exp(-i m' phi) d^j_{m'm}(theta).
CplxMatrix rotation_matrix(double j, double theta, double phi);

/// (J_x, J_y, J_z) in the +J..-J ordering.
struct AngularMomentumOps {
  CplxMatrix jx, jy, jz;
};
AngularMomentumOps angular_momentum_ops(double j);

/// Irreducible tensor operator T_kq with
/// <j m'|T_kq|j m> = <j m; k q|j m'> sqrt((2k+1)/(2j+1)).
CplxMatrix multipole_operator(double j, int k, int q);

/// Lab-frame density matrix for populations p_m along axis n(theta, phi).
/// Throws std::runtime_error on an unphysical / mis-sized population vector.
CplxMatrix rho_from_populations(double j, const std::vector<double>& populations,
                                double theta = 0.0, double phi = 0.0);

/// <J>/J as a real 3-vector (lab frame).
std::array<double, 3> vector_polarization(const CplxMatrix& rho, double j);

/// Normalized J_z tensor moment: J=1 -> Pzz = <3Jz^2-2>;
/// J=3/2 -> T = <3Jz^2 - J(J+1)>/3.
double tensor_polarization(const CplxMatrix& rho, double j);

/// Normalized rank-3 moment for J=3/2: <Jz^3 - (41/20) Jz>/(3/10).
double octupole_moment(const CplxMatrix& rho, double j = 1.5);

/// (vector, tensor[, octupole]) moments of p_m along its own axis.
struct AxisMoments {
  double vector = 0.0;
  double tensor = 0.0;
  double octupole = 0.0;
  bool has_octupole = false;
};
AxisMoments moments_along_axis(double j, const std::vector<double>& populations);

/// (p+, p0, p-) with <Jz> = pz and <3Jz^2-2> = pzz.
std::array<double, 3> spin1_populations(double pz, double pzz);

/// (p_{3/2}, p_{1/2}, p_{-1/2}, p_{-3/2}) with the normalized moments
/// (vector pz, tensor t, octupole o) of this module.
std::array<double, 4> spin32_populations(double pz, double t = 0.0,
                                         double o = 0.0);

/// Spin-temperature (maximum-entropy) populations p_m ~ exp(beta m) with
/// <J_z>/J = pz.  Any j; beta by bisection, exactly as the Python does it.
std::vector<double> populations_maxent(double j, double pz,
                                       double beta_max = 60.0,
                                       double tol = 1e-13);

/// Populations along an axis + the axis orientation, with lab moments.
struct SpinDensity {
  double j = 1.0;
  std::vector<double> populations;  ///< m = +J ... -J along the axis
  double theta = 0.0;               ///< axis polar angle in the lab [rad]
  double phi = 0.0;                 ///< axis azimuth in the lab [rad]

  CplxMatrix rho_lab() const;

  struct LabMoments {
    std::array<double, 3> vector{{0.0, 0.0, 0.0}};
    double tensor_zz = 0.0;
    double octupole_z = 0.0;
    bool has_octupole = false;
  };
  LabMoments lab_moments() const;
};

}  // namespace lipolgen

#endif  // LIPOLGEN_SPIN_HPP
