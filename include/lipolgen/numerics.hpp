// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LIPOLGEN_NUMERICS_HPP
#define LIPOLGEN_NUMERICS_HPP

/// \file numerics.hpp
/// The handful of NumPy primitives the physics kernel depends on, transcribed
/// so that C++ and Python agree to the last bit rather than merely to rtol.
///
/// Three of them matter numerically:
///   * `pairwise_sum` -- NumPy's summation order (`np.add.reduce`), which is
///     NOT a left-to-right accumulation.  `np.trapz` / `np.trapezoid` go
///     through it, and the Wandzura-Wilczek g2 and the Delta sum rule are
///     both trapezoid integrals compared at rtol 1e-12.
///   * `np_interp` -- `np.interp`, including its clamped extrapolation, its
///     exact-node short circuit and its non-finite fallback.  Every digitized
///     b1 / EMC curve is read through it.
///   * `linspace` / `logspace` -- NumPy's `start + i*step` with the endpoint
///     written back exactly, which is not the same as `start + i*(b-a)/(n-1)`.

#include <cstddef>
#include <cmath>
#include <vector>

namespace lipolgen {

/// NumPy's pairwise summation (numpy/core/src/umath/loops_utils.h.src,
/// `pairwise_sum_@TYPE@`), block size 128.  Reproducing it is what makes a
/// trapezoid integral bit-identical to `np.trapz`.
double pairwise_sum(const double* a, std::size_t n);

inline double pairwise_sum(const std::vector<double>& v) {
  return pairwise_sum(v.data(), v.size());
}

/// `np.trapz(y, x)` = sum(diff(x) * (y[1:] + y[:-1]) / 2), summed pairwise.
double trapezoid(const std::vector<double>& y, const std::vector<double>& x);

/// `np.interp(x, xp, fp)`: linear interpolation, clamped to the end values
/// outside [xp[0], xp[n-1]].  `xp` must be strictly increasing.
double np_interp(double x, const double* xp, const double* fp, std::size_t n);

inline double np_interp(double x, const std::vector<double>& xp,
                        const std::vector<double>& fp) {
  return np_interp(x, xp.data(), fp.data(), xp.size());
}

/// `np.linspace(a, b, n)` with endpoint=True: y_i = a + i*(b-a)/(n-1) and
/// y[n-1] forced to exactly b.
std::vector<double> linspace(double a, double b, std::size_t n);

/// `np.logspace(a, b, n)` = 10 ** linspace(a, b, n).
std::vector<double> logspace(double a, double b, std::size_t n);

/// `np.clip`.
inline double clip(double v, double lo, double hi) {
  return v < lo ? lo : (v > hi ? hi : v);
}

}  // namespace lipolgen

#endif  // LIPOLGEN_NUMERICS_HPP
