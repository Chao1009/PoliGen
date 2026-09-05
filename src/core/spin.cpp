#include "lipolgen/spin.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

namespace lipolgen {
namespace {

/// Factorial of a value that must be a non-negative integer (`spin._fact`).
double fact(double n) {
  const long k = std::lround(n);
  if (std::fabs(n - static_cast<double>(k)) > 1e-9 || k < 0) {
    throw std::runtime_error("factorial argument is not a non-negative integer");
  }
  double out = 1.0;
  for (long i = 2; i <= k; ++i) out *= static_cast<double>(i);
  return out;
}

std::size_t dim_of(double j) {
  return static_cast<std::size_t>(std::lround(2.0 * j + 1.0));
}

/// Gaussian elimination with partial pivoting; used only for the 4x4
/// spin-3/2 moment inversion (`np.linalg.solve`).
std::vector<double> solve_linear(std::vector<std::vector<double>> a,
                                 std::vector<double> b) {
  const std::size_t n = b.size();
  for (std::size_t col = 0; col < n; ++col) {
    std::size_t piv = col;
    for (std::size_t r = col + 1; r < n; ++r) {
      if (std::fabs(a[r][col]) > std::fabs(a[piv][col])) piv = r;
    }
    if (std::fabs(a[piv][col]) < 1e-300) {
      throw std::runtime_error("singular moment matrix");
    }
    std::swap(a[col], a[piv]);
    std::swap(b[col], b[piv]);
    for (std::size_t r = col + 1; r < n; ++r) {
      const double f = a[r][col] / a[col][col];
      if (f == 0.0) continue;
      for (std::size_t c = col; c < n; ++c) a[r][c] -= f * a[col][c];
      b[r] -= f * b[col];
    }
  }
  std::vector<double> x(n, 0.0);
  for (std::size_t i = n; i-- > 0;) {
    double s = b[i];
    for (std::size_t c = i + 1; c < n; ++c) s -= a[i][c] * x[c];
    x[i] = s / a[i][i];
  }
  return x;
}

}  // namespace

CplxMatrix matmul(const CplxMatrix& a, const CplxMatrix& b) {
  const std::size_t n = a.size();
  CplxMatrix out(n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t k = 0; k < n; ++k) {
      const std::complex<double> aik = a(i, k);
      if (aik == std::complex<double>(0.0, 0.0)) continue;
      for (std::size_t jj = 0; jj < n; ++jj) out(i, jj) += aik * b(k, jj);
    }
  }
  return out;
}

CplxMatrix adjoint(const CplxMatrix& a) {
  const std::size_t n = a.size();
  CplxMatrix out(n);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t jj = 0; jj < n; ++jj) out(i, jj) = std::conj(a(jj, i));
  }
  return out;
}

std::complex<double> trace(const CplxMatrix& a) {
  std::complex<double> s(0.0, 0.0);
  for (std::size_t i = 0; i < a.size(); ++i) s += a(i, i);
  return s;
}

std::vector<double> m_values(double j) {
  const std::size_t n = dim_of(j);
  std::vector<double> ms(n);
  for (std::size_t i = 0; i < n; ++i) ms[i] = j - static_cast<double>(i);
  return ms;
}

double clebsch_gordan(double j1, double m1, double j2, double m2,
                      double j, double m) {
  if (std::fabs(m1 + m2 - m) > 1e-9) return 0.0;
  if (!(std::fabs(j1 - j2) - 1e-9 <= j && j <= j1 + j2 + 1e-9)) return 0.0;
  double pre = ((2.0 * j + 1.0) * fact(j1 + j2 - j) * fact(j1 - j2 + j)
                * fact(-j1 + j2 + j) / fact(j1 + j2 + j + 1.0));
  pre *= (fact(j1 + m1) * fact(j1 - m1) * fact(j2 + m2) * fact(j2 - m2)
          * fact(j + m) * fact(j - m));
  const long k_min = std::lround(std::max(std::max(0.0, j2 - j - m1), j1 - j + m2));
  const long k_max = std::lround(std::min(std::min(j1 + j2 - j, j1 - m1), j2 + m2));
  double tot = 0.0;
  for (long k = k_min; k <= k_max; ++k) {
    const double kd = static_cast<double>(k);
    tot += (std::pow(-1.0, kd)
            / (fact(kd) * fact(j1 + j2 - j - kd) * fact(j1 - m1 - kd)
               * fact(j2 + m2 - kd) * fact(j - j2 + m1 + kd)
               * fact(j - j1 - m2 + kd)));
  }
  return std::sqrt(pre) * tot;
}

RealMatrix wigner_d(double j, double beta) {
  const std::vector<double> ms = m_values(j);
  const std::size_t n = ms.size();
  RealMatrix d(n);
  const double c = std::cos(beta / 2.0);
  const double s = std::sin(beta / 2.0);
  for (std::size_t a = 0; a < n; ++a) {
    const double mp = ms[a];
    for (std::size_t b = 0; b < n; ++b) {
      const double m = ms[b];
      const double pre = std::sqrt(fact(j + mp) * fact(j - mp)
                                   * fact(j + m) * fact(j - m));
      const long k_min = std::lround(std::max(0.0, m - mp));
      const long k_max = std::lround(std::min(j + m, j - mp));
      double tot = 0.0;
      for (long k = k_min; k <= k_max; ++k) {
        const double kd = static_cast<double>(k);
        const double num = std::pow(-1.0, mp - m + kd);
        const double den = (fact(j + m - kd) * fact(kd)
                            * fact(j - mp - kd) * fact(mp - m + kd));
        tot += (num / den
                * std::pow(c, 2.0 * j + m - mp - 2.0 * kd)
                * std::pow(s, mp - m + 2.0 * kd));
      }
      d(a, b) = pre * tot;
    }
  }
  return d;
}

CplxMatrix rotation_matrix(double j, double theta, double phi) {
  const std::vector<double> ms = m_values(j);
  const RealMatrix d = wigner_d(j, theta);
  const std::size_t n = ms.size();
  CplxMatrix u(n);
  for (std::size_t a = 0; a < n; ++a) {
    const std::complex<double> ph = std::exp(std::complex<double>(0.0, -ms[a] * phi));
    for (std::size_t b = 0; b < n; ++b) u(a, b) = ph * d(a, b);
  }
  return u;
}

AngularMomentumOps angular_momentum_ops(double j) {
  const std::vector<double> ms = m_values(j);
  const std::size_t n = ms.size();
  AngularMomentumOps ops;
  ops.jz = CplxMatrix(n);
  CplxMatrix jp(n);
  for (std::size_t b = 0; b < n; ++b) {
    ops.jz(b, b) = std::complex<double>(ms[b], 0.0);
    if (b > 0) {  // m+1 sits at index b-1 (descending order)
      const double m = ms[b];
      jp(b - 1, b) = std::complex<double>(
          std::sqrt(j * (j + 1.0) - m * (m + 1.0)), 0.0);
    }
  }
  const CplxMatrix jm = adjoint(jp);
  ops.jx = CplxMatrix(n);
  ops.jy = CplxMatrix(n);
  for (std::size_t a = 0; a < n; ++a) {
    for (std::size_t b = 0; b < n; ++b) {
      ops.jx(a, b) = 0.5 * (jp(a, b) + jm(a, b));
      ops.jy(a, b) = std::complex<double>(0.0, -0.5) * (jp(a, b) - jm(a, b));
    }
  }
  return ops;
}

CplxMatrix multipole_operator(double j, int k, int q) {
  const std::vector<double> ms = m_values(j);
  const std::size_t n = ms.size();
  CplxMatrix t(n);
  const double norm = std::sqrt((2.0 * k + 1.0) / (2.0 * j + 1.0));
  for (std::size_t a = 0; a < n; ++a) {
    for (std::size_t b = 0; b < n; ++b) {
      t(a, b) = std::complex<double>(
          clebsch_gordan(j, ms[b], static_cast<double>(k),
                         static_cast<double>(q), j, ms[a]) * norm, 0.0);
    }
  }
  return t;
}

CplxMatrix rho_from_populations(double j, const std::vector<double>& populations,
                                double theta, double phi) {
  const std::size_t n = dim_of(j);
  if (populations.size() != n) {
    throw std::runtime_error("populations must have 2j+1 entries (m=+J..-J)");
  }
  double sum = 0.0;
  for (double p : populations) {
    if (p < -1e-12) throw std::runtime_error("populations must be >= 0");
    sum += p;
  }
  if (std::fabs(sum - 1.0) > 1e-9) {
    throw std::runtime_error("populations must sum to 1");
  }
  const CplxMatrix u = rotation_matrix(j, theta, phi);
  CplxMatrix diag(n);
  for (std::size_t i = 0; i < n; ++i) diag(i, i) = std::complex<double>(populations[i], 0.0);
  return matmul(matmul(u, diag), adjoint(u));
}

std::array<double, 3> vector_polarization(const CplxMatrix& rho, double j) {
  const AngularMomentumOps ops = angular_momentum_ops(j);
  const CplxMatrix* three[3] = {&ops.jx, &ops.jy, &ops.jz};
  std::array<double, 3> out{{0.0, 0.0, 0.0}};
  for (int i = 0; i < 3; ++i) {
    out[static_cast<std::size_t>(i)] =
        trace(matmul(rho, *three[i])).real() / j;
  }
  return out;
}

double tensor_polarization(const CplxMatrix& rho, double j) {
  const AngularMomentumOps ops = angular_momentum_ops(j);
  const std::size_t n = ops.jz.size();
  CplxMatrix op = matmul(ops.jz, ops.jz);
  for (std::size_t a = 0; a < n; ++a) {
    for (std::size_t b = 0; b < n; ++b) {
      op(a, b) = 3.0 * op(a, b) - (a == b ? std::complex<double>(j * (j + 1.0), 0.0)
                                          : std::complex<double>(0.0, 0.0));
    }
  }
  const double q = trace(matmul(rho, op)).real();
  if (std::fabs(j - 1.0) < 1e-9) return q;
  if (std::fabs(j - 1.5) < 1e-9) return q / 3.0;
  throw std::runtime_error("tensor_polarization defined for j = 1, 3/2 only");
}

double octupole_moment(const CplxMatrix& rho, double j) {
  if (std::fabs(j - 1.5) > 1e-9) {
    throw std::runtime_error("octupole moment defined for j = 3/2 only");
  }
  const AngularMomentumOps ops = angular_momentum_ops(j);
  CplxMatrix op = matmul(matmul(ops.jz, ops.jz), ops.jz);
  const std::size_t n = op.size();
  for (std::size_t a = 0; a < n; ++a) {
    for (std::size_t b = 0; b < n; ++b) {
      op(a, b) = op(a, b) - (41.0 / 20.0) * ops.jz(a, b);
    }
  }
  return trace(matmul(rho, op)).real() / 0.3;
}

AxisMoments moments_along_axis(double j, const std::vector<double>& populations) {
  const std::vector<double> ms = m_values(j);
  if (populations.size() != ms.size()) {
    throw std::runtime_error("populations must have 2j+1 entries (m=+J..-J)");
  }
  AxisMoments out;
  double vsum = 0.0;
  for (std::size_t i = 0; i < ms.size(); ++i) vsum += ms[i] * populations[i];
  out.vector = vsum / j;
  if (std::fabs(j - 1.0) < 1e-9) {
    double t = 0.0;
    for (std::size_t i = 0; i < ms.size(); ++i) {
      t += (3.0 * ms[i] * ms[i] - 2.0) * populations[i];
    }
    out.tensor = t;
    return out;
  }
  if (std::fabs(j - 1.5) < 1e-9) {
    double t = 0.0, o = 0.0;
    for (std::size_t i = 0; i < ms.size(); ++i) {
      t += (3.0 * ms[i] * ms[i] - j * (j + 1.0)) * populations[i];
      o += (ms[i] * ms[i] * ms[i] - (41.0 / 20.0) * ms[i]) * populations[i];
    }
    out.tensor = t / 3.0;
    out.octupole = o / 0.3;
    out.has_octupole = true;
    return out;
  }
  throw std::runtime_error("moments defined for j = 1, 3/2 only");
}

std::array<double, 3> spin1_populations(double pz, double pzz) {
  const double p_plus = (2.0 + 3.0 * pz + pzz) / 6.0;
  const double p_zero = (1.0 - pzz) / 3.0;
  const double p_minus = (2.0 - 3.0 * pz + pzz) / 6.0;
  std::array<double, 3> pops{{p_plus, p_zero, p_minus}};
  for (double p : pops) {
    if (p < -1e-12) {
      throw std::runtime_error("unphysical (pz, pzz): negative population");
    }
  }
  for (double& p : pops) p = p < 0.0 ? 0.0 : (p > 1.0 ? 1.0 : p);
  return pops;
}

std::array<double, 4> spin32_populations(double pz, double t, double o) {
  const std::vector<double> ms = m_values(1.5);
  std::vector<std::vector<double>> a(4, std::vector<double>(4, 0.0));
  for (std::size_t i = 0; i < 4; ++i) {
    a[0][i] = 1.0;
    a[1][i] = ms[i] / 1.5;
    a[2][i] = (3.0 * ms[i] * ms[i] - 3.75) / 3.0;
    a[3][i] = (ms[i] * ms[i] * ms[i] - (41.0 / 20.0) * ms[i]) / 0.3;
  }
  const std::vector<double> sol = solve_linear(a, {1.0, pz, t, o});
  std::array<double, 4> pops{};
  for (std::size_t i = 0; i < 4; ++i) {
    if (sol[i] < -1e-12) {
      // SAY WHICH PAIR IS REACHABLE (`phase_D_li7_rank2.md` sec. 1.5, defect
      // F3).  "unphysical (pz, t, o): negative population" was correct and
      // useless: the caller cannot tell whether the vector moment, the
      // alignment or the pair is at fault, and the J = 3/2 domain is SMALLER
      // than the spin-1 one, so ordinary-looking numbers land outside it
      // (`make_plan("helicity-flip", j=1.5, pz=0.7, pzz=0.6)` does).
      //
      // The four populations are fixed UNIQUELY by (1, pz, t, o) -- the 4x4
      // system above is square -- so there is no fill to search for, only a
      // domain.  Inverting the same system by hand:
      //
      //     p(+3/2) + p(-3/2) = (1 + t)/2 ,  p(+3/2) - p(-3/2) = 0.9 pz + 0.1 o
      //     p(+1/2) + p(-1/2) = (1 - t)/2 ,  p(+1/2) - p(-1/2) = 0.3 (pz - o)
      //
      // so positivity is exactly the two inequalities below, and at o = 0 it
      // collapses to 1.8|pz| - 1 <= t <= 1 - 0.6|pz| (hence |pz| <= 5/6).
      // VERIFIED against this function on 1681 (pz, t) points at o = 0 and
      // 4851 (pz, t, o) points: zero disagreements.
      const double lo = 1.8 * std::fabs(pz) - 1.0;
      const double hi = 1.0 - 0.6 * std::fabs(pz);
      char buf[768];
      std::snprintf(
          buf, sizeof(buf),
          "spin32_populations: unphysical (pz = %g, t = %g, o = %g) -- "
          "p(m = %s) = %g < 0.  The four populations are fixed UNIQUELY by "
          "(1, pz, t, o), so this is a DOMAIN, not a solver failure: "
          "p(+3/2) -+ p(-3/2) = (1 + t)/2 and 0.9 pz + 0.1 o, "
          "p(+1/2) -+ p(-1/2) = (1 - t)/2 and 0.3 (pz - o), whence "
          "|0.9 pz + 0.1 o| <= (1 + t)/2 AND |0.3 (pz - o)| <= (1 - t)/2.  "
          "At o = 0 that is %g <= t <= %g at this pz (so |pz| <= 5/6 = "
          "0.833333 at all).  The J = 3/2 domain is SMALLER than the spin-1 "
          "one -- (pz, pzz) = (0.7, 0.6) is inside it for spin 1 and outside "
          "here by 0.02 in t",
          pz, t, o,
          (i == 0 ? "+3/2" : i == 1 ? "+1/2" : i == 2 ? "-1/2" : "-3/2"),
          sol[i], lo, hi);
      throw std::runtime_error(buf);
    }
    pops[i] = sol[i] < 0.0 ? 0.0 : (sol[i] > 1.0 ? 1.0 : sol[i]);
  }
  return pops;
}

std::vector<double> populations_maxent(double j, double pz, double beta_max,
                                       double tol) {
  if (!(-1.0 < pz && pz < 1.0)) {
    throw std::runtime_error("maxent populations need |pz| < 1");
  }
  const std::vector<double> ms = m_values(j);
  const double mmax = ms.front();  // +J, the maximum by construction

  auto mean = [&](double beta) {
    double num = 0.0, den = 0.0;
    for (double m : ms) {
      const double w = std::exp(beta * (m - mmax));
      num += m * w;
      den += w;
    }
    return num / den / j;
  };

  double lo = -beta_max, hi = beta_max;
  if (!(mean(lo) < pz && pz < mean(hi))) {
    throw std::runtime_error("pz out of reach for the bisection bracket");
  }
  while (hi - lo > tol) {
    const double mid = 0.5 * (lo + hi);
    if (mean(mid) < pz) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  const double beta = 0.5 * (lo + hi);
  std::vector<double> w(ms.size());
  double tot = 0.0;
  for (std::size_t i = 0; i < ms.size(); ++i) {
    w[i] = std::exp(beta * (ms[i] - mmax));
    tot += w[i];
  }
  for (double& v : w) v /= tot;
  return w;
}

CplxMatrix SpinDensity::rho_lab() const {
  return rho_from_populations(j, populations, theta, phi);
}

SpinDensity::LabMoments SpinDensity::lab_moments() const {
  const CplxMatrix rho = rho_lab();
  LabMoments out;
  out.vector = vector_polarization(rho, j);
  out.tensor_zz = tensor_polarization(rho, j);
  if (std::fabs(j - 1.5) < 1e-9) {
    out.octupole_z = octupole_moment(rho, j);
    out.has_octupole = true;
  }
  return out;
}

}  // namespace lipolgen
