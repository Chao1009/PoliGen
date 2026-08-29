#include "lipolgen/asymmetries.hpp"

#include <algorithm>
#include <cmath>

namespace lipolgen {

double depolarization_d(double y, double x, double q2, const RFunc& r_func) {
  const double r = resolve_r(r_func, x, q2);
  return y * (2.0 - y) / (y * y + 2.0 * (1.0 - y) * (1.0 + r));
}

double a_parallel(double g1, double f1, double y, double x, double q2,
                  const RFunc& r_func) {
  const double d = depolarization_d(y, x, q2, r_func);
  return d * g1 / std::max(f1, 1e-30);
}

double a_perp(double g1, double g2, double f1, double y, double x, double q2,
              const RFunc& r_func) {
  const double eps = (1.0 - y) / (1.0 - y + 0.5 * y * y);
  const double d = depolarization_d(y, x, q2, r_func)
                   * std::sqrt(std::max(2.0 * eps / (1.0 + eps), 0.0));
  const double gamma = 2.0 * M_NUCLEON * x / std::sqrt(q2);
  const double amp = gamma * (0.5 * y * g1 + g2) / std::max(f1, 1e-30);
  return d * amp;
}

double phi_averaged_density(double f1, double f2, double x, double y) {
  return f1 + (1.0 - y) / (x * y * y) * f2;
}

double azz(double b1, double f1, double f2, double x, double y,
           const double* b2, double theta_m) {
  const double b2v = (b2 == nullptr) ? 2.0 * x * b1 : *b2;
  const double ct = std::cos(theta_m);
  const double geom = 0.5 * (3.0 * ct * ct - 1.0);  // = 1 at theta_m = 0
  const double num = TENSOR_LL_SIGN * 2.0 / 3.0
                     * (b1 + (1.0 - y) / (x * y * y) * b2v) * geom;
  return num / phi_averaged_density(f1, f2, x, y);
}

double a_cos2phi(double delta, double f1, double f2, double x, double y) {
  return -(1.0 - y) / (y * y) * delta / phi_averaged_density(f1, f2, x, y);
}

double err_a_parallel(double n, double pe, double pz) {
  n = std::max(n, 1e-12);
  return 1.0 / (pe * pz * std::sqrt(n));
}

double err_azz(double n, double pzz) {
  n = std::max(n, 1e-12);
  return std::sqrt(2.0 / n) / pzz;
}

double err_cos2phi_amplitude(double n, double pzz) {
  n = std::max(n, 1e-12);
  return std::sqrt(2.0 / n) / pzz;
}

}  // namespace lipolgen
