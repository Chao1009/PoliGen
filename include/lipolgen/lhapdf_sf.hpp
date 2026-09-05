// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LIPOLGEN_LHAPDF_SF_HPP
#define LIPOLGEN_LHAPDF_SF_HPP

/// \file lhapdf_sf.hpp
/// LHAPDF-only additions that `sf.hpp` does not declare: the quiet-mode
/// switch and the EPPS21 per-flavour nuclear-ratio hook.
///
/// `LhapdfSF` / `LhapdfG1` themselves are declared in `sf.hpp` (P2 leaves
/// them declared, implemented here in `src/lhapdf/lhapdf_sf.cpp`) so the
/// core library headers stay the single point of truth for the backend
/// interface; this header exists only for the LHAPDF-specific extras that
/// have no toy/table counterpart to sit next to.
///
/// Port of the EPPS21 use in `fastsim/scripts/money_delta_20260729.py`
/// (`NuclearF2FromGrid`, `EPPS21_SET = "EPPS21nlo_CT18Anlo_Li6"`): that
/// script evaluates F2^A directly from the bound-nucleon EPPS21 grid and
/// multiplies by A for the whole nucleus (the `parton` package returns
/// nuclear PDFs PER NUCLEON, confirmed there by cross-check against
/// CT18NLO).  `Epps21Ratio` reproduces the identical per-nucleon F2 through
/// an explicit per-flavour ratio R_f(x,Q2) = f_A(x,Q2)/f_p(x,Q2) applied on
/// top of a proton baseline (default CT18NLO): algebraically
/// R_f * xf_p = xf_A, so `f2_per_nucleon` is bit-identical to reading the
/// nuclear grid directly, while `ratio()` exposes the per-flavour EPPS21
/// nuclear-modification factor itself for callers that want R_f(x,Q2) as a
/// standalone hook (e.g. multiplying a *different* proton PDF or a single
/// flavour) rather than only the F2-level combination.

#include <memory>
#include <string>

namespace lipolgen {

/// `LHAPDF::setVerbosity(0)`, called at most once (idempotent) regardless of
/// how many LHAPDF-backed objects are constructed.  All of `LhapdfSF`,
/// `LhapdfG1` and `Epps21Ratio` call this in their constructors, so a caller
/// never needs to call it directly; it is exposed for tests that want
/// silence before the first construction too.
void lhapdf_quiet();

/// EPPS21 nuclear-modification ratio, per LHAPDF/PDG parton id, applied on
/// top of a free-proton baseline -- see the file docstring.  Default grids:
/// nuclear = EPPS21nlo_CT18Anlo_Li6 (6Li, the only nucleus installed),
/// proton = CT18NLO (matching `PartonF2`/`LhapdfSF`, NOT the nominal
/// CT18ANLO EPPS21 was fit on; the two proton baselines differ from each
/// other by a few tenths of a percent at DIS kinematics, see
/// `lhapdf_sf.cpp` for the measured spread -- pass `proton_set` to change
/// it).
class Epps21Ratio {
 public:
  explicit Epps21Ratio(std::string nuclear_set = "EPPS21nlo_CT18Anlo_Li6",
                       int nuclear_member = 0,
                       std::string proton_set = "CT18NLO",
                       int proton_member = 0);
  ~Epps21Ratio();

  /// R_f(pid, x, q2) = f_A(pid,x,q2) / f_p(pid,x,q2).  `pid` is the
  /// LHAPDF/PDG parton id: 1=d, 2=u, 3=s, 4=c, 5=b, 21=g, negative =
  /// antiquark.  Returns 0 where the proton baseline vanishes (guards a
  /// division that `money_delta_20260729.py` never has to take because it
  /// reads the nuclear grid directly rather than a ratio).
  double ratio(int pid, double x, double q2) const;

  /// Per-nucleon F2 of the nucleus the grid was built for, PartonF2's
  /// 5-flavour (d u s c b) e_q^2 convention, built as
  /// Sigma_q e_q^2 [R_f(q) xf_p(q) + R_f(qbar) xf_p(qbar)] -- algebraically
  /// Sigma_q e_q^2 (xf_A(q) + xf_A(qbar)), i.e. `NuclearF2FromGrid.f2a(x,
  /// q2) / ion.A`.
  double f2_per_nucleon(double x, double q2) const;

  /// Whole-nucleus F2A = A * f2_per_nucleon: `NuclearF2FromGrid.f2a`'s own
  /// normalisation (A the nucleon count, e.g. 6 for 6Li).
  double f2a(double x, double q2, int A) const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace lipolgen

#endif  // LIPOLGEN_LHAPDF_SF_HPP
