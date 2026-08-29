#ifndef LIPOLGEN_CLUSTER_HPP
#define LIPOLGEN_CLUSTER_HPP

/// \file cluster.hpp
/// Cluster relative partial waves: the radial shapes and the azimuth-stripped
/// spherical harmonics the two-cluster amplitude of `tagged.hpp` is built from.
/// C++17 port of `tagged.Wave` / `tagged.theta_lm`
/// (evgen/polligen/tagged.py:126,101).
///
/// The radial forms are UNNORMALIZED; `TaggedModel` normalizes them on its own
/// grid so that integral psihat^2 k^2 dk = 1.  kappa comes from the channel's
/// separation energy (spectator.hpp), beta is the short-range scale scanned as
/// the model band (0.20 - 0.40, default 0.30):
///
///   L = 0  Hulthen   1/(k^2+kappa^2) - 1/(k^2+beta^2)
///   L = 1  P-wave    k / ((k^2+kappa^2)(k^2+beta^2))
///   L = 2  D-wave    k^2 / ((k^2+kappa^2)(k^2+beta^2)^2)
///
/// The L = 0 form is EXACTLY `spectator.momentum_density`'s psi, which is what
/// makes the P_D = 0 limit of the tagged sampler agree with the fast
/// simulation's spectator sampler quantile by quantile.
///
/// A `Backend`-style replacement (VMC two-cluster overlap densities, R. B.
/// Wiringa et al.) lands as a table behind the same `radial()` interface.

#include <vector>

#include "lipolgen/spectator.hpp"

namespace lipolgen {

/// One cluster relative partial wave: orbital L, probability P_L = a_L^2, and
/// the short-range scale beta of its radial shape.
struct Wave {
  int l = 0;
  double prob = 1.0;
  double beta = BETA_DEFAULT;

  /// Unnormalized radial shape psi_L(k).  Throws for L > 2.
  double radial(double k, double kappa) const;
  /// The same on a grid.
  std::vector<double> radial(const std::vector<double>& k, double kappa) const;
};

/// |m|-azimuth-stripped spherical harmonic Theta_l^m(theta) with
/// Condon-Shortley signs: Y_l^m = Theta_l^m(theta) exp(i m phi).  `c` is
/// cos(theta).  Implemented for l <= 2 (throws otherwise), which is every
/// wave the two-cluster expansion carries.
double theta_lm(int l, int m, double c);

}  // namespace lipolgen

#endif  // LIPOLGEN_CLUSTER_HPP
