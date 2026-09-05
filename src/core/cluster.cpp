#include "lipolgen/cluster.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>

#include "lipolgen/constants.hpp"
#include "lipolgen/numerics.hpp"

#ifndef LIPOLGEN_DATA_DIR_DEFAULT
#define LIPOLGEN_DATA_DIR_DEFAULT "data"
#endif

namespace lipolgen {
namespace {

/// GeV fm (CODATA hbar c).  The VMC tables are in fm / fm^-1 throughout.
constexpr double kHbarCGeVfm = 0.1973269804;

/// Spherical Bessel j_l for l = 0, 1, 2, with the small-x series where the
/// closed forms lose all their significant digits to cancellation.
double sph_bessel(int l, double x) {
  const double x2 = x * x;
  if (l == 0) {
    if (std::fabs(x) < 1e-3) return 1.0 - x2 / 6.0 + x2 * x2 / 120.0;
    return std::sin(x) / x;
  }
  if (l == 1) {
    if (std::fabs(x) < 1e-2) return x / 3.0 * (1.0 - x2 / 10.0 + x2 * x2 / 280.0);
    return std::sin(x) / x2 - std::cos(x) / x;
  }
  if (l == 2) {
    if (std::fabs(x) < 1e-1) {
      return x2 / 15.0 * (1.0 - x2 / 14.0 + x2 * x2 / 504.0);
    }
    const double x3 = x2 * x;
    return (3.0 / x3 - 1.0 / x) * std::sin(x) - 3.0 / x2 * std::cos(x);
  }
  throw std::runtime_error("sph_bessel implemented for l <= 2");
}

bool is_number_start(const std::string& tok) {
  if (tok.empty()) return false;
  const char c = tok[0];
  return std::isdigit(static_cast<unsigned char>(c)) || c == '-' || c == '+'
         || c == '.';
}

/// Split a Fortran-output row into numeric fields.  The overlap files print
/// `-0.3212E-01 (.1021E-01)-0.2882E-02` with NO space before a leading minus,
/// so a plain whitespace split loses fields; parentheses are separators too.
std::vector<double> split_numbers(const std::string& line) {
  std::vector<double> out;
  const char* p = line.c_str();
  while (*p) {
    if (*p == '(' || *p == ')' || std::isspace(static_cast<unsigned char>(*p))) {
      ++p;
      continue;
    }
    char* end = nullptr;
    const double v = std::strtod(p, &end);
    if (end == p) return out;   // a non-numeric token: the block has ended
    out.push_back(v);
    p = end;
  }
  return out;
}

std::string read_whole(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open ANL table: " + path);
  std::ostringstream ss;
  ss << in.rdbuf();
  return ss.str();
}

std::vector<std::string> split_lines(const std::string& text) {
  std::vector<std::string> out;
  std::istringstream in(text);
  std::string line;
  while (std::getline(in, line)) out.push_back(line);
  return out;
}

std::string lstrip(const std::string& s) {
  std::size_t i = 0;
  while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
  return s.substr(i);
}

/// Rows of `x v0 (e0) v1 (e1) [ints...]` following `hdr`, stopping at the
/// first line that does not start with a number.  `ncol` value columns.
AnlTable read_pair_block(const std::vector<std::string>& lines, std::size_t hdr,
                         std::size_t ncol) {
  AnlTable t;
  t.col.assign(ncol, {});
  t.err.assign(ncol, {});
  for (std::size_t i = hdr + 1; i < lines.size(); ++i) {
    const std::string s = lstrip(lines[i]);
    std::istringstream probe(s);
    std::string tok;
    probe >> tok;
    if (!is_number_start(tok)) {
      if (!t.x.empty()) break;
      continue;
    }
    const std::vector<double> v = split_numbers(lines[i]);
    if (v.size() < 1 + 2 * ncol) {
      if (!t.x.empty()) break;
      continue;
    }
    t.x.push_back(v[0]);
    for (std::size_t c = 0; c < ncol; ++c) {
      t.col[c].push_back(v[1 + 2 * c]);
      t.err[c].push_back(v[2 + 2 * c]);
    }
  }
  if (t.x.empty()) throw std::runtime_error("ANL block has no rows");
  return t;
}

}  // namespace

// -------------------------------------------------------------- data lookup

const std::string& data_dir() {
  static const std::string dir = [] {
    const char* env = std::getenv("LIPOLGEN_DATA_DIR");
    if (env != nullptr && *env != '\0') return std::string(env);
    return std::string(LIPOLGEN_DATA_DIR_DEFAULT);
  }();
  return dir;
}

std::string data_path(const std::string& relative) {
  const std::string& d = data_dir();
  if (d.empty()) return relative;
  if (d.back() == '/') return d + relative;
  return d + "/" + relative;
}

// --------------------------------------------------------- the ANL readers

std::vector<AnlTable> read_anl_overlap(const std::string& path) {
  const std::vector<std::string> lines = split_lines(read_whole(path));
  std::size_t k_hdr = lines.size(), r_hdr = lines.size();
  for (std::size_t i = 0; i < lines.size(); ++i) {
    const std::string s = lstrip(lines[i]);
    if (k_hdr == lines.size() && s.rfind("k(fm-1)", 0) == 0) k_hdr = i;
    if (r_hdr == lines.size() && s.rfind("rij", 0) == 0) r_hdr = i;
  }
  if (k_hdr == lines.size() || r_hdr == lines.size()) {
    throw std::runtime_error("not an overlap_old raw output (no k(fm-1)/rij "
                             "block): " + path);
  }
  return {read_pair_block(lines, k_hdr, 2), read_pair_block(lines, r_hdr, 2)};
}

std::vector<AnlTable> read_anl_momentum(const std::string& path) {
  const std::vector<std::string> lines = split_lines(read_whole(path));
  std::vector<AnlTable> out;
  for (std::size_t i = 0; i < lines.size(); ++i) {
    // The column rule ` ****  ***********  *********` introduces each block.
    std::string bare;
    for (char c : lines[i]) {
      if (!std::isspace(static_cast<unsigned char>(c))) bare.push_back(c);
    }
    if (bare.size() < 5 || bare.find_first_not_of('*') != std::string::npos) {
      continue;
    }
    AnlTable t;
    std::size_t ncol = 0;
    for (std::size_t j = i + 1; j < lines.size(); ++j) {
      const std::vector<double> v = split_numbers(lines[j]);
      if (v.size() < 3) break;
      if (ncol == 0) {
        ncol = (v.size() - 1) / 2;
        t.col.assign(ncol, {});
        t.err.assign(ncol, {});
      }
      if ((v.size() - 1) / 2 < ncol) break;
      t.x.push_back(v[0]);
      for (std::size_t c = 0; c < ncol; ++c) {
        t.col[c].push_back(v[1 + 2 * c]);
        t.err[c].push_back(v[2 + 2 * c]);
      }
    }
    if (!t.x.empty()) {
      out.push_back(std::move(t));
      i += out.back().x.size();
    }
  }
  if (out.empty()) {
    throw std::runtime_error("no momentum block found in " + path);
  }
  return out;
}

AnlTable read_anl_plain(const std::string& path, std::size_t ncol) {
  if (ncol == 0) throw std::runtime_error("read_anl_plain: ncol must be >= 1");
  const std::vector<std::string> lines = split_lines(read_whole(path));
  AnlTable t;
  t.col.assign(ncol, {});
  for (std::size_t i = 0; i < lines.size(); ++i) {
    // The column rule ` ****  *********  *********` introduces the block --
    // the same marker `read_anl_momentum` uses.
    std::string bare;
    for (char c : lines[i]) {
      if (!std::isspace(static_cast<unsigned char>(c))) bare.push_back(c);
    }
    if (bare.size() < 4 || bare.find_first_not_of('*') != std::string::npos) {
      continue;
    }
    for (std::size_t j = i + 1; j < lines.size(); ++j) {
      const std::vector<double> v = split_numbers(lines[j]);
      if (v.size() < 1 + ncol) break;
      t.x.push_back(v[0]);
      for (std::size_t c = 0; c < ncol; ++c) t.col[c].push_back(v[1 + c]);
    }
    if (!t.x.empty()) return t;
  }
  throw std::runtime_error("no plain `x v1 ... vN` block found in " + path);
}

std::vector<double> read_anl_momentum_norms(const std::string& path) {
  std::vector<double> out;
  for (const std::string& line : split_lines(read_whole(path))) {
    const std::size_t at = line.find("/(2*PI)**3");
    if (at == std::string::npos) continue;
    const std::size_t eq = line.find('=', at);
    if (eq == std::string::npos) continue;
    out.push_back(std::strtod(line.c_str() + eq + 1, nullptr));
  }
  return out;
}

FdeutTable read_fdeut_k(const std::string& path) {
  const std::vector<std::string> lines = split_lines(read_whole(path));
  // The file carries THREE `k ...` blocks (`k ges gms` at line 10018, the wave
  // functions at 10221, and the form-factor blocks below); only one names the
  // wave functions, so key on `u(k)` rather than on the abscissa.
  std::size_t hdr = lines.size(), ebind_hdr = lines.size();
  for (std::size_t i = 0; i < lines.size(); ++i) {
    const std::string s = lstrip(lines[i]);
    if (hdr == lines.size() && s.rfind("k ", 0) == 0
        && s.find("u(k)") != std::string::npos) {
      hdr = i;
    }
    if (ebind_hdr == lines.size() && s.rfind("ebind", 0) == 0) ebind_hdr = i;
  }
  if (hdr == lines.size()) {
    throw std::runtime_error("no `k u(k) w(k)` block in " + path);
  }
  if (ebind_hdr == lines.size()) {
    throw std::runtime_error("no `ebind` header in " + path);
  }
  FdeutTable t;
  // `ebind` is printed in MeV on the first numeric line under its header.
  for (std::size_t i = ebind_hdr + 1; i < lines.size(); ++i) {
    const std::vector<double> v = split_numbers(lines[i]);
    if (v.empty()) continue;
    t.ebind_gev = v[0] * 1e-3;
    break;
  }
  // sqrt(hbar c)^3: the ordinate is in fm^(3/2) and the abscissa in fm^-1, so
  // BOTH move or the norm shifts by hbar c^3 (see the header).
  const double scale = 1.0 / (HBARC_GEV_FM * std::sqrt(HBARC_GEV_FM));
  for (std::size_t i = hdr + 1; i < lines.size(); ++i) {
    const std::string s = lstrip(lines[i]);
    std::istringstream probe(s);
    std::string tok;
    probe >> tok;
    if (!is_number_start(tok)) {
      if (!t.k_gev.empty()) break;
      continue;
    }
    const std::vector<double> v = split_numbers(lines[i]);
    if (v.size() != 3) {
      if (!t.k_gev.empty()) break;
      continue;
    }
    t.k_gev.push_back(v[0] * HBARC_GEV_FM);
    t.u.push_back(v[1] * scale);
    t.w.push_back(v[2] * scale);
  }
  if (t.k_gev.size() < 2) {
    throw std::runtime_error("empty `k u(k) w(k)` block in " + path);
  }
  for (std::size_t i = 1; i < t.k_gev.size(); ++i) {
    if (!(t.k_gev[i] > t.k_gev[i - 1])) {
      throw std::runtime_error("fdeut k column is not increasing: " + path);
    }
  }
  return t;
}

// --------------------------------------------------------- the CD-Bonn deuteron

namespace {

/// Table XX plus the four boundary conditions.  See `CdBonnWave` in
/// cluster.hpp for the parameterisation, the sign and the provenance.
CdBonnWave build_cdbonn_wave() {
  // n = 11 is not a knob: the three-unknown constraint solve below and the
  // `std::array<double, 11>` members are written for it.  If the table ever
  // grows, these fire rather than the wave function quietly changing.
  static_assert(CD_BONN_N == 11, "the CD-Bonn solve below assumes n = 11");
  constexpr std::size_t kNc = sizeof(CD_BONN_C) / sizeof(CD_BONN_C[0]);
  constexpr std::size_t kNd = sizeof(CD_BONN_D) / sizeof(CD_BONN_D[0]);
  static_assert(kNc == 10, "Table XX publishes C_1..C_10; C_11 is Eq. (D23)");
  static_assert(kNd == 8, "Table XX publishes D_1..D_8; D_9..D_11 are Eq. (D24)");
  CdBonnWave w;
  static_assert(std::tuple_size<decltype(w.m)>::value == CD_BONN_N, "size");
  for (int j = 0; j < CD_BONN_N; ++j) {
    w.m[static_cast<std::size_t>(j)] =
        CD_BONN_GAMMA_FM + static_cast<double>(j) * CD_BONN_M0_FM;
  }
  // Eq. (D23), i.e. u_a(0) = 0: sum_j C_j = 0.
  double sum_c = 0.0;
  for (std::size_t j = 0; j < kNc; ++j) {
    w.c[j] = CD_BONN_C[j];
    sum_c += w.c[j];
  }
  w.c[CD_BONN_N - 1] = -sum_c;

  // Eq. (D24) and its two circular permutations, in the equivalent sum-rule
  // form of the header:  sum_j D_j/m_j^2 = sum_j D_j = sum_j D_j m_j^2 = 0.
  // In x_j = m_j^2 and E_j = D_j/x_j the three conditions are
  //     sum_{j=9,10,11} E_j x_j^i = b_i ,   i = 0, 1, 2 ,
  // i.e. a VANDERMONDE system, whose inverse is Lagrange interpolation:
  // with L_k(t) = (t - x_a)(t - x_b)/[(x_k - x_a)(x_k - x_b)] the solution is
  //     E_k = [x_a x_b b_0 - (x_a + x_b) b_1 + b_2] / [(x_k-x_a)(x_k-x_b)] .
  // No matrix solve, no permutation to get wrong, and no pivoting question.
  double b0 = 0.0, b1 = 0.0, b2 = 0.0;
  for (std::size_t j = 0; j < kNd; ++j) {
    const double x = w.m[j] * w.m[j];
    w.d[j] = CD_BONN_D[j];
    b0 -= w.d[j] / x;
    b1 -= w.d[j];
    b2 -= w.d[j] * x;
  }
  const std::size_t last[3] = {kNd, kNd + 1, kNd + 2};
  for (int t = 0; t < 3; ++t) {
    const std::size_t k = last[t];
    const std::size_t a = last[(t + 1) % 3];
    const std::size_t b = last[(t + 2) % 3];
    const double xk = w.m[k] * w.m[k];
    const double xa = w.m[a] * w.m[a];
    const double xb = w.m[b] * w.m[b];
    w.d[k] = xk * (xa * xb * b0 - (xa + xb) * b1 + b2)
             / ((xk - xa) * (xk - xb));
  }
  return w;
}

/// sum_ij A_i B_j/(m_i + m_j) -- the closed form of
/// (2/pi) int_0^inf dp p^2 psi_A psi_B when both are Eq. (D21)/(D22) sums,
/// from (2/pi) int_0^inf dp p^2/[(p^2+a^2)(p^2+b^2)] = 1/(a+b).
double cdbonn_moment(const std::array<double, 11>& a,
                     const std::array<double, 11>& b,
                     const std::array<double, 11>& m) {
  double s = 0.0;
  for (std::size_t i = 0; i < m.size(); ++i)
    for (std::size_t j = 0; j < m.size(); ++j) s += a[i] * b[j] / (m[i] + m[j]);
  return s;
}

}  // namespace

double CdBonnWave::psi_s(double p_fm) const {
  const double p2 = p_fm * p_fm;
  double s = 0.0;
  for (std::size_t j = 0; j < m.size(); ++j) s += c[j] / (p2 + m[j] * m[j]);
  return std::sqrt(2.0 / kPi) * s;
}

double CdBonnWave::psi_d(double p_fm) const {
  const double p2 = p_fm * p_fm;
  double s = 0.0;
  for (std::size_t j = 0; j < m.size(); ++j) s += d[j] / (p2 + m[j] * m[j]);
  // MINUS: w = -psi_2^a, the `fdeut.av18` / CDKS convention.  See the header.
  return -std::sqrt(2.0 / kPi) * s;
}

double CdBonnWave::norm_s() const { return cdbonn_moment(c, c, m); }
double CdBonnWave::norm_d() const { return cdbonn_moment(d, d, m); }
double CdBonnWave::norm() const { return norm_s() + norm_d(); }
double CdBonnWave::a_s() const { return c[0]; }
double CdBonnWave::a_d() const { return d[0]; }
double CdBonnWave::eta() const { return d[0] / c[0]; }

std::array<double, 4> CdBonnWave::constraint_residuals() const {
  std::array<double, 4> r{0.0, 0.0, 0.0, 0.0};
  for (std::size_t j = 0; j < m.size(); ++j) {
    const double x = m[j] * m[j];
    r[0] += c[j];
    r[1] += d[j] / x;
    r[2] += d[j];
    r[3] += d[j] * x;
  }
  return r;
}

const CdBonnWave& cdbonn_wave() {
  static const CdBonnWave w = build_cdbonn_wave();
  return w;
}

FdeutTable cdbonn_fdeut_table(double k_max_fm, double dk_fm) {
  if (!(dk_fm > 0.0)) {
    throw std::runtime_error("cdbonn_fdeut_table: dk_fm must be positive");
  }
  if (!(k_max_fm >= dk_fm)) {
    throw std::runtime_error("cdbonn_fdeut_table: k_max_fm < dk_fm");
  }
  const CdBonnWave& w = cdbonn_wave();
  // The SAME "both or neither" conversion `read_fdeut_k` applies: the
  // abscissa is fm^-1 -> GeV and the ordinate fm^3/2 -> GeV^-3/2, or the norm
  // moves by hbar c^3.
  const double scale = 1.0 / (HBARC_GEV_FM * std::sqrt(HBARC_GEV_FM));
  const std::size_t n =
      static_cast<std::size_t>(std::llround(k_max_fm / dk_fm)) + 1;
  FdeutTable t;
  t.k_gev.reserve(n);
  t.u.reserve(n);
  t.w.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    const double k = static_cast<double>(i) * dk_fm;
    t.k_gev.push_back(k * HBARC_GEV_FM);
    t.u.push_back(w.psi_s(k) * scale);
    t.w.push_back(w.psi_d(k) * scale);
  }
  t.ebind_gev = CD_BONN_BINDING_MEV * 1e-3;
  return t;
}

// ----------------------------------------------------------------- VmcRadial

VmcRadial::VmcRadial(std::vector<double> k_gev, std::vector<double> psi, int l,
                     std::string provenance)
    : k_(std::move(k_gev)), psi_(std::move(psi)), l_(l),
      provenance_(std::move(provenance)) {
  if (k_.size() != psi_.size() || k_.size() < 2) {
    throw std::runtime_error("VmcRadial: need >= 2 matching (k, psi) points");
  }
  for (std::size_t i = 1; i < k_.size(); ++i) {
    if (!(k_[i] > k_[i - 1])) {
      throw std::runtime_error("VmcRadial: k must be strictly increasing");
    }
  }
}

double VmcRadial::operator()(double k) const {
  if (k_.empty() || k < k_.front() || k > k_.back()) return 0.0;
  return np_interp(k, k_, psi_);
}

VmcRadial::VmcRadial(std::vector<double> k_gev, std::vector<double> psi,
                     std::vector<double> dpsi, int l, std::string provenance)
    : VmcRadial(std::move(k_gev), std::move(psi), l, std::move(provenance)) {
  if (!dpsi.empty() && dpsi.size() != psi_.size()) {
    throw std::runtime_error("VmcRadial: dpsi must match psi or be empty");
  }
  dpsi_ = std::move(dpsi);
}

VmcRadial VmcRadial::scaled(double factor) const {
  std::vector<double> p(psi_), d(dpsi_);
  for (double& v : p) v *= factor;
  for (double& v : d) v *= factor;
  return VmcRadial(k_, std::move(p), std::move(d), l_, provenance_);
}

VmcRadial VmcRadial::shifted_by_sigma(double n_sigma) const {
  if (n_sigma == 0.0 || dpsi_.empty()) return *this;
  std::vector<double> p(psi_);
  for (std::size_t i = 0; i < p.size(); ++i) p[i] += n_sigma * dpsi_[i];
  std::ostringstream prov;
  prov << provenance_ << " [+" << n_sigma
       << " sigma_MC, FULLY CORRELATED across the table]";
  return VmcRadial(k_, std::move(p), dpsi_, l_, prov.str());
}

double VmcRadial::norm2() const {
  std::vector<double> y(k_.size());
  for (std::size_t i = 0; i < k_.size(); ++i) y[i] = k_[i] * k_[i] * psi_[i] * psi_[i];
  return trapezoid(y, k_);
}

void VmcRadial::norm2_error(double* correlated, double* quadrature) const {
  double cor = 0.0, quad = 0.0;
  if (!dpsi_.empty()) {
    // d/dn integral k^2 (psi + n dpsi)^2 dk at n = 0, cell by cell, so the
    // two limits use exactly the same cells and the same quadrature.
    for (std::size_t i = 1; i < k_.size(); ++i) {
      const double a = 2.0 * k_[i - 1] * k_[i - 1] * psi_[i - 1] * dpsi_[i - 1];
      const double b = 2.0 * k_[i] * k_[i] * psi_[i] * dpsi_[i];
      const double cell = 0.5 * (a + b) * (k_[i] - k_[i - 1]);
      cor += cell;
      quad += cell * cell;
    }
  }
  if (correlated != nullptr) *correlated = std::fabs(cor);
  if (quadrature != nullptr) *quadrature = std::sqrt(quad);
}

std::vector<double> vmc_sign_steps(const std::vector<double>& x,
                                   const std::vector<double>& amp,
                                   double x_max) {
  std::vector<double> nodes;
  for (std::size_t i = 1; i < x.size(); ++i) {
    if (x[i] > x_max) break;
    const double a = amp[i - 1], b = amp[i];
    if (a == 0.0 || b == 0.0 || (a > 0.0) == (b > 0.0)) continue;
    nodes.push_back(x[i - 1] - a * (x[i] - x[i - 1]) / (b - a));
  }
  return nodes;
}

VmcRadial vmc_from_overlap_k(const std::string& path, int column, int l) {
  const std::vector<AnlTable> b = read_anl_overlap(path);
  const AnlTable& t = b[0];
  if (column < 0 || static_cast<std::size_t>(column) >= t.col.size()) {
    throw std::runtime_error("vmc_from_overlap_k: no such column");
  }
  std::vector<double> k(t.x.size());
  for (std::size_t i = 0; i < t.x.size(); ++i) k[i] = t.x[i] * kHbarCGeVfm;
  // The overlap block's own 1-sigma MC error column travels with the
  // amplitude (C5.2); it is empty only if the file printed none.
  std::vector<double> dpsi;
  const std::size_t c = static_cast<std::size_t>(column);
  if (c < t.err.size() && t.err[c].size() == t.x.size()) dpsi = t.err[c];
  return VmcRadial(std::move(k), t.col[c], std::move(dpsi), l,
                   path + " [k-space block, column " + std::to_string(column)
                       + "]");
}

VmcRadial vmc_from_overlap_r(const std::string& path, int column, int l,
                             double k_max_fm, std::size_t nk) {
  const std::vector<AnlTable> b = read_anl_overlap(path);
  const AnlTable& t = b[1];
  if (column < 0 || static_cast<std::size_t>(column) >= t.col.size()) {
    throw std::runtime_error("vmc_from_overlap_r: no such column");
  }
  const std::vector<double>& r = t.x;
  const std::vector<double>& a = t.col[static_cast<std::size_t>(column)];
  std::vector<double> k_fm = linspace(0.0, k_max_fm, nk);
  std::vector<double> psi(nk), integrand(r.size());
  for (std::size_t i = 0; i < nk; ++i) {
    for (std::size_t j = 0; j < r.size(); ++j) {
      integrand[j] = a[j] * r[j] * r[j] * sph_bessel(l, k_fm[i] * r[j]);
    }
    psi[i] = 4.0 * kPi * trapezoid(integrand, r);
  }
  std::vector<double> k_gev(nk);
  for (std::size_t i = 0; i < nk; ++i) k_gev[i] = k_fm[i] * kHbarCGeVfm;
  return VmcRadial(std::move(k_gev), std::move(psi), l,
                   path + " [r-space block, column " + std::to_string(column)
                       + ", Fourier-Bessel 4pi convention]");
}

VmcRadial vmc_from_momentum(const std::string& path, int block, int column,
                            int l, const VmcRadial* sign_from,
                            double node_search_max_fm) {
  const std::vector<AnlTable> b = read_anl_momentum(path);
  if (block < 0 || static_cast<std::size_t>(block) >= b.size()) {
    throw std::runtime_error("vmc_from_momentum: no such block in " + path);
  }
  const AnlTable& t = b[static_cast<std::size_t>(block)];
  if (column < 0 || static_cast<std::size_t>(column) >= t.col.size()) {
    throw std::runtime_error("vmc_from_momentum: no such column");
  }
  const std::vector<double>& rho = t.col[static_cast<std::size_t>(column)];

  // Nodes of the SIGNED reference amplitude, in fm^-1, plus an ANCHOR: the
  // reference's own sign at the point where it is largest, which is the one
  // place its Monte Carlo sign is beyond doubt.  Counting nodes from a
  // k -> 0 seed instead would throw away exactly the information that
  // matters -- both 6Li columns are node-free at low k, so the S-D RELATIVE
  // phase lives entirely in the anchor.
  std::vector<double> nodes;
  std::size_t anchor_below = 0;
  double anchor_sign = 1.0;
  if (sign_from != nullptr && !sign_from->empty()) {
    std::vector<double> ref_x(sign_from->k().size());
    for (std::size_t i = 0; i < ref_x.size(); ++i) {
      ref_x[i] = sign_from->k()[i] / kHbarCGeVfm;
    }
    nodes = vmc_sign_steps(ref_x, sign_from->psi(), node_search_max_fm);
    std::size_t best = 0;
    for (std::size_t i = 0; i < ref_x.size(); ++i) {
      if (ref_x[i] > node_search_max_fm) break;
      if (std::fabs(sign_from->psi()[i]) > std::fabs(sign_from->psi()[best])) {
        best = i;
      }
    }
    anchor_sign = sign_from->psi()[best] < 0.0 ? -1.0 : 1.0;
    for (double n : nodes) {
      if (ref_x[best] > n) ++anchor_below;
    }
  }

  // psi = s sqrt(rho), so the file's own 1-sigma drho becomes
  // dpsi = s drho / (2 sqrt(rho)) -- SIGNED, so it transforms with psi under
  // the global phase (C5.2).  Empty if the file printed no error column.
  const std::size_t ic = static_cast<std::size_t>(column);
  const bool have_err = ic < t.err.size() && t.err[ic].size() == t.x.size();
  const std::vector<double>* drho = have_err ? &t.err[ic] : nullptr;
  std::vector<double> k(t.x.size()), psi(t.x.size());
  std::vector<double> dpsi(have_err ? t.x.size() : 0);
  for (std::size_t i = 0; i < t.x.size(); ++i) {
    k[i] = t.x[i] * kHbarCGeVfm;
    std::size_t below = 0;
    for (double n : nodes) {
      if (t.x[i] > n) ++below;
    }
    const std::size_t flips = below > anchor_below ? below - anchor_below
                                                   : anchor_below - below;
    const double s = anchor_sign * ((flips % 2 == 0) ? 1.0 : -1.0);
    const double r = std::fmax(rho[i], 0.0);
    psi[i] = s * std::sqrt(r);
    if (have_err) dpsi[i] = r > 0.0 ? s * (*drho)[i] / (2.0 * std::sqrt(r))
                                    : 0.0;
  }
  std::string prov = path + " [block " + std::to_string(block) + ", column "
                     + std::to_string(column) + ", psi = sqrt(rho)";
  if (!nodes.empty()) {
    prov += ", sign from " + sign_from->provenance() + " (nodes at";
    for (double n : nodes) {
      char buf[32];
      std::snprintf(buf, sizeof buf, " %.3f", n);
      prov += buf;
    }
    prov += " fm^-1)";
  }
  prov += have_err ? ", 1-sigma MC errors carried]" : "]";
  return VmcRadial(std::move(k), std::move(psi), std::move(dpsi), l, prov);
}

// ---------------------------------------------------------------- the waves

double Wave::radial(double k, double kappa) const {
  if (vmc) return (*vmc)(k);
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
