#include "lipolgen/cluster.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

#include "lipolgen/constants.hpp"

namespace lipolgen {

double Wave::radial(double k, double kappa) const {
  const double k2 = k * k;
  const double b2 = beta * beta;
  const double kap2 = kappa * kappa;
  switch (l) {
    case 0:
      return 1.0 / (k2 + kap2) - 1.0 / (k2 + b2);
    case 1:
      return k / ((k2 + kap2) * (k2 + b2));
    case 2:
      return k2 / ((k2 + kap2) * (k2 + b2) * (k2 + b2));
    default:
      throw std::runtime_error("Wave::l must be 0, 1, or 2");
  }
}

std::vector<double> Wave::radial(const std::vector<double>& k,
                                 double kappa) const {
  std::vector<double> out(k.size());
  for (std::size_t i = 0; i < k.size(); ++i) out[i] = radial(k[i], kappa);
  return out;
}

double theta_lm(int l, int m, double c) {
  const double s = std::sqrt(std::fmax(1.0 - c * c, 0.0));
  if (l == 0 && m == 0) return std::sqrt(1.0 / (4.0 * kPi));
  if (l == 1) {
    if (m == 0) return std::sqrt(3.0 / (4.0 * kPi)) * c;
    if (m == 1 || m == -1) {
      return -static_cast<double>(m > 0 ? 1 : -1) * std::sqrt(3.0 / (8.0 * kPi)) * s;
    }
  }
  if (l == 2) {
    if (m == 0) return std::sqrt(5.0 / (16.0 * kPi)) * (3.0 * c * c - 1.0);
    if (m == 1 || m == -1) {
      return -static_cast<double>(m > 0 ? 1 : -1)
             * std::sqrt(15.0 / (8.0 * kPi)) * s * c;
    }
    if (m == 2 || m == -2) return std::sqrt(15.0 / (32.0 * kPi)) * s * s;
  }
  throw std::runtime_error("theta_lm implemented for l <= 2, got (l="
                           + std::to_string(l) + ", m=" + std::to_string(m) + ")");
}

}  // namespace lipolgen
