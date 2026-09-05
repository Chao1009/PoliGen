// SPDX-License-Identifier: GPL-3.0-or-later
#include "lipolgen/numerics.hpp"

#include <algorithm>

namespace lipolgen {
namespace {
constexpr std::size_t kPwBlocksize = 128;
}  // namespace

double pairwise_sum(const double* a, std::size_t n) {
  if (n < 8) {
    double res = 0.0;
    for (std::size_t i = 0; i < n; ++i) res += a[i];
    return res;
  }
  if (n <= kPwBlocksize) {
    double r[8];
    for (std::size_t k = 0; k < 8; ++k) r[k] = a[k];
    std::size_t i = 8;
    for (; i < n - (n % 8); i += 8) {
      for (std::size_t k = 0; k < 8; ++k) r[k] += a[i + k];
    }
    double res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
    for (; i < n; ++i) res += a[i];
    return res;
  }
  std::size_t n2 = n / 2;
  n2 -= n2 % 8;
  return pairwise_sum(a, n2) + pairwise_sum(a + n2, n - n2);
}

double trapezoid(const std::vector<double>& y, const std::vector<double>& x) {
  const std::size_t n = y.size();
  if (n < 2) return 0.0;
  std::vector<double> terms(n - 1);
  for (std::size_t i = 0; i + 1 < n; ++i) {
    terms[i] = (x[i + 1] - x[i]) * (y[i + 1] + y[i]) / 2.0;
  }
  return pairwise_sum(terms);
}

double np_interp(double x, const double* xp, const double* fp, std::size_t n) {
  if (n == 0) return 0.0;
  if (n == 1) return fp[0];
  if (std::isnan(x)) return x;
  // numpy's binary_search_with_guess: j = -1 below the table, n above it.
  if (x > xp[n - 1]) return fp[n - 1];
  if (x < xp[0]) return fp[0];
  std::size_t j = static_cast<std::size_t>(
      std::upper_bound(xp, xp + n, x) - xp);
  j = (j == 0) ? 0 : j - 1;
  if (j == n - 1) return fp[n - 1];
  if (xp[j] == x) return fp[j];  // exact node: avoid a non-finite slope
  const double slope = (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j]);
  double res = slope * (x - xp[j]) + fp[j];
  if (!std::isfinite(res)) {
    res = slope * (x - xp[j + 1]) + fp[j + 1];
    if (!std::isfinite(res) && fp[j] == fp[j + 1]) res = fp[j];
  }
  return res;
}

std::vector<double> linspace(double a, double b, std::size_t n) {
  std::vector<double> out(n);
  if (n == 0) return out;
  if (n == 1) {
    out[0] = a;
    return out;
  }
  const double step = (b - a) / static_cast<double>(n - 1);
  for (std::size_t i = 0; i < n; ++i) {
    out[i] = static_cast<double>(i) * step + a;
  }
  out[n - 1] = b;
  return out;
}

std::vector<double> logspace(double a, double b, std::size_t n) {
  std::vector<double> out = linspace(a, b, n);
  for (double& v : out) v = std::pow(10.0, v);
  return out;
}

}  // namespace lipolgen
