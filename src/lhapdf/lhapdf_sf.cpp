// LHAPDF6 backends: `LhapdfSF` / `LhapdfG1` (declared in include/lipolgen/sf.hpp)
// and `Epps21Ratio` (declared in include/lipolgen/lhapdf_sf.hpp).
//
// Ports, formula for formula, of the LHAPDF-grid classes in
// PolarizedLithiumSim/fastsim/polli_fastsim/{structure,polarized}.py:
//   structure.PartonF2   -> LhapdfSF   (CT18NLO default; F2 = x Sum e_q^2 (q+qbar),
//                                       5 flavours d u s c b, neutron by u<->d swap)
//   polarized.PartonG1   -> LhapdfG1   (NNPDFpol11_100 default; g1 = (1/2) Sum e_q^2
//                                       (Dq+Dqbar), 3 flavours d u s -- Delta c =
//                                       Delta b = 0 is NNPDFpol1.1's own convention,
//                                       not an omission here)
// Nucleus combinations (F2A = Z F2p + N F2n, g1A = Z P_p g1p + N P_n g1n) are NOT
// duplicated here: wrap a `LhapdfSF`/`LhapdfG1` in the existing `NuclearF2` /
// `PolSF::g1_nucleus` from sf.hpp (both already isospin/effective-polarization
// generic and compiled into lipolgen_core).
//
// `Epps21Ratio` ports the EPPS21 use of `fastsim/scripts/money_delta_20260729.py`'s
// `NuclearF2FromGrid` (EPPS21_SET = "EPPS21nlo_CT18Anlo_Li6") -- see lhapdf_sf.hpp.
//
// `_safe_xfx` in the Python guards a NumPy>=2 0-d array bug that does not exist in
// C++: `LHAPDF::PDF::xfxQ2` already returns a plain double, so no analogous guard
// is needed here.

#include "lipolgen/lhapdf_sf.hpp"
#include "lipolgen/sf.hpp"

#include <algorithm>
#include <array>
#include <mutex>
#include <utility>

#include "LHAPDF/LHAPDF.h"

namespace lipolgen {

namespace {

void ensure_quiet_impl() { LHAPDF::setVerbosity(0); }

// F2 (PartonF2._E2): five flavours d u s c b, proton e_q^2 weights.
constexpr std::array<std::pair<int, double>, 5> kF2ProtonE2 = {{
    {1, 1.0 / 9.0}, {2, 4.0 / 9.0}, {3, 1.0 / 9.0}, {4, 4.0 / 9.0}, {5, 1.0 / 9.0}}};
// F2 neutron (PartonF2.f2n's e2n dict): isospin u<->d swap of the charges above.
constexpr std::array<std::pair<int, double>, 5> kF2NeutronE2 = {{
    {1, 4.0 / 9.0}, {2, 1.0 / 9.0}, {3, 1.0 / 9.0}, {4, 4.0 / 9.0}, {5, 1.0 / 9.0}}};
// g1 (PartonG1._E2): three flavours d u s only -- Delta c = Delta b = 0.
constexpr std::array<int, 3> kG1Pids = {1, 2, 3};

double e2_dus(int pid) {
  // proton e_q^2 for pid in {1 (d), 2 (u), 3 (s)}, matching PartonG1._E2.
  switch (pid) {
    case 1: return 1.0 / 9.0;
    case 2: return 4.0 / 9.0;
    default: return 1.0 / 9.0;  // pid == 3 (s)
  }
}

double f2_from_weights(const LHAPDF::PDF& pdf, double x, double q2,
                       const std::array<std::pair<int, double>, 5>& weights) {
  // 0 < x < 1 guard: PartonF2._f2p_scalar / f2n's scalar both return 0.0
  // outside the open unit interval before touching the grid.
  if (!(x > 0.0 && x < 1.0)) return 0.0;
  double tot = 0.0;
  for (const auto& pe : weights) {
    tot += pe.second * (pdf.xfxQ2(pe.first, x, q2) + pdf.xfxQ2(-pe.first, x, q2));
  }
  return std::max(tot, 0.0);
}

double g1_scalar(const LHAPDF::PDF& pdf, double x, double q2, bool swap_ud) {
  // PartonG1._g1_scalar: 0 < x < 1 guard, xfxQ2 already returns x*Delta-q,
  // so g1 = 0.5 * sum / x.
  if (!(x > 0.0 && x < 1.0)) return 0.0;
  double tot = 0.0;
  for (int pid : kG1Pids) {
    double e2 = e2_dus(pid);
    if (swap_ud && (pid == 1 || pid == 2)) e2 = e2_dus(3 - pid);  // isospin swap
    tot += e2 * (pdf.xfxQ2(pid, x, q2) + pdf.xfxQ2(-pid, x, q2));
  }
  return 0.5 * tot / x;
}

}  // namespace

void lhapdf_quiet() {
  static std::once_flag flag;
  std::call_once(flag, ensure_quiet_impl);
}

// ------------------------------------------------------------------ LhapdfSF

struct LhapdfSF::Impl {
  std::unique_ptr<LHAPDF::PDF> pdf;
};

LhapdfSF::LhapdfSF(std::string setname, int member) : impl_(new Impl) {
  lhapdf_quiet();
  impl_->pdf.reset(LHAPDF::mkPDF(setname, static_cast<std::size_t>(member)));
}

LhapdfSF::~LhapdfSF() = default;

double LhapdfSF::q2_min() const { return impl_->pdf->q2Min(); }

double LhapdfSF::f2p(double x, double q2) const {
  return f2_from_weights(*impl_->pdf, x, q2, kF2ProtonE2);
}

double LhapdfSF::f2n(double x, double q2) const {
  return f2_from_weights(*impl_->pdf, x, q2, kF2NeutronE2);
}

double LhapdfSF::f2n_over_f2p(double x) const {
  // PartonF2.f2n_over_f2p(x, q2=10.0): the Q2 default the Python signature
  // carries: q2 is optional there and unused by any call site in the
  // Python (there isn't one -- the interface's abstract f2n_over_f2p is
  // single-argument, exactly as ToyF2's is), so it is fixed here, not
  // threaded through.
  constexpr double kQ2 = 10.0;
  const double fp = f2p(x, kQ2);
  if (!(fp > 0.0)) return 1.0;  // PartonF2.f2n_over_f2p: np.where(f2p > 0, ..., 1.0)
  return f2n(x, kQ2) / std::max(fp, 1e-30);
}

// ------------------------------------------------------------------ LhapdfG1

struct LhapdfG1::Impl {
  std::unique_ptr<LHAPDF::PDF> pdf;
};

LhapdfG1::LhapdfG1(std::string setname, int member) : impl_(new Impl) {
  lhapdf_quiet();
  impl_->pdf.reset(LHAPDF::mkPDF(setname, static_cast<std::size_t>(member)));
}

LhapdfG1::~LhapdfG1() = default;

double LhapdfG1::g1p(double x, double q2) const {
  return g1_scalar(*impl_->pdf, x, q2, false);
}

double LhapdfG1::g1n(double x, double q2) const {
  return g1_scalar(*impl_->pdf, x, q2, true);
}

// -------------------------------------------------------------- Epps21Ratio

struct Epps21Ratio::Impl {
  std::unique_ptr<LHAPDF::PDF> nuclear;
  std::unique_ptr<LHAPDF::PDF> proton;
};

Epps21Ratio::Epps21Ratio(std::string nuclear_set, int nuclear_member,
                         std::string proton_set, int proton_member)
    : impl_(new Impl) {
  lhapdf_quiet();
  impl_->nuclear.reset(
      LHAPDF::mkPDF(nuclear_set, static_cast<std::size_t>(nuclear_member)));
  impl_->proton.reset(
      LHAPDF::mkPDF(proton_set, static_cast<std::size_t>(proton_member)));
}

Epps21Ratio::~Epps21Ratio() = default;

double Epps21Ratio::ratio(int pid, double x, double q2) const {
  const double fp = impl_->proton->xfxQ2(pid, x, q2);
  if (fp == 0.0) return 0.0;
  const double fa = impl_->nuclear->xfxQ2(pid, x, q2);
  return fa / fp;
}

double Epps21Ratio::f2_per_nucleon(double x, double q2) const {
  if (!(x > 0.0 && x < 1.0)) return 0.0;
  double tot = 0.0;
  for (const auto& pe : kF2ProtonE2) {
    const int pid = pe.first;
    const double e2 = pe.second;
    const double xfp_q = impl_->proton->xfxQ2(pid, x, q2);
    const double xfp_qbar = impl_->proton->xfxQ2(-pid, x, q2);
    const double r_q = ratio(pid, x, q2);
    const double r_qbar = ratio(-pid, x, q2);
    tot += e2 * (r_q * xfp_q + r_qbar * xfp_qbar);
  }
  return std::max(tot, 0.0);
}

double Epps21Ratio::f2a(double x, double q2, int A) const {
  return f2_per_nucleon(x, q2) * static_cast<double>(A);
}

}  // namespace lipolgen
