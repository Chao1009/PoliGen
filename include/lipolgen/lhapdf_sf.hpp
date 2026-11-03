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
/// other by 0.13 % in the valence-window mean of R, i.e. 4.0 % relative on <1 - R> (0.029803 on CT18NLO vs 0.031052 on CT18ANLO, 03_data_nuclear.md sec. 1.2); see
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

// ===================================================== D7 (2026-10-27)
//
// SINCE D7 THE SHIPPED STRUCTURE-FUNCTION DEFAULTS ARE THIS TIER'S GRIDS:
// `UNPOL_SF_DEFAULT` = ct18nlo is `LhapdfSF("CT18NLO", 0)` and `POL_SF_DEFAULT`
// = nnpdfpol is `LhapdfG1("NNPDFpol11_100", 0)` (pipeline.hpp), both evaluated
// at LEADING ORDER (the development item of DECISIONS.md D7).  The core cannot
// build them, so this tier does, and attaches them to the configs that name
// them -- lazily, through the filler a `Pipeline` calls (`install_lhapdf_sf_
// filler`, which the Python module runs at import), or per config
// (`attach_shipped_sf_backends`).  Nothing here is ever replaced by a toy.
#include "lipolgen/pipeline.hpp"

namespace lipolgen {

/// How EVERY LHAPDF grid this tier opens is evaluated, set ON THE OBJECT with
/// LHAPDF's per-PDF setters so that no lhapdf.conf on any machine can move a
/// number (none of the three installed sets' .info files carries these keys;
/// measured 2026-10-27, `Extrapolator: nearest` against `continuation` in the
/// global lhapdf.conf moved the 6Li default cross section +5.03 %).  The
/// author's choice (D7): FREEZE every set at its grid edges -- `nearest`, the
/// knot-array edge, which is what the letter's `--pdf grid` chain does at Q0
/// -- with LHAPDF's standard `logcubic` interpolation and no positivity
/// forcing (ForcePositive 0: Delta q may be negative).
inline constexpr const char* LHAPDF_EXTRAPOLATOR = "nearest";
inline constexpr const char* LHAPDF_INTERPOLATOR = "logcubic";
inline constexpr int LHAPDF_FORCE_POSITIVE = 0;

/// The PRE-FLIGHT: "" when LHAPDF can open member `member` of `set`, else the
/// refusal every constructor of this tier throws -- built from
/// `LHAPDF::findFile` / `paths()` ALONE, before any LHAPDF Config exists, so a
/// missing lhapdf.conf can no longer abort the process at exit (Config's
/// function-static re-throws from its own destructor: std::terminate, exit 134,
/// measured 2026-10-27).  A missing set is reported in LHAPDF's own words,
/// "Info file not found for PDF set '<name>'", which
/// validation/benchmarks/_report.py reads as "blocked (environment)"; the
/// message names the searched paths, LHAPDF_DATA_PATH, lhapdf.conf's state and
/// the escape (`--unpol-sf toy` / `--pol-sf toy`).
std::string lhapdf_set_problem(const std::string& set, int member = 0);

/// The LHAPDF set a grid selector names: "CT18NLO" for `UnpolSfSource::
/// Ct18Nlo`, "NNPDFpol11_100" for `PolSfSource::NnpdfPol`, nullptr otherwise.
const char* lhapdf_set_name(UnpolSfSource source);
const char* lhapdf_set_name(PolSfSource source);

/// The object a grid selector names -- ONE per process, built on first use
/// under a mutex and shared read-only thereafter (a 6Li inclusive `Pipeline`
/// took 0.21 s against 0.06 s on the toy, and every `Pipeline` re-loading the
/// 23-28 ms grid would only repeat that).  Throws the pre-flight message on a
/// missing set and `std::invalid_argument` on a selector that names no grid.
std::shared_ptr<const UnpolSF> make_lhapdf_unpol_sf(UnpolSfSource source);
std::shared_ptr<const PolSF> make_lhapdf_pol_sf(PolSfSource source);

/// Fill every EMPTY slot of `cfg` whose selector names a grid backend --
/// `unpol_sf_obj` on `Ct18Nlo`, `pol_sf_obj` on `NnpdfPol` -- with the shared
/// objects above; a slot already holding an object is left alone.  NOT
/// `b1_unpol_sf`: `b1_unpol` stays toy by default and a C++ program attaches
/// `make_lhapdf_unpol_sf(UnpolSfSource::Ct18Nlo)` there itself, as before D7.
/// What a C++ program calls before building a `PythiaBridge` that must share
/// the run's backend (the T2 species draw).
void attach_shipped_sf_backends(PipelineConfig& cfg);

/// `set_grid_sf_filler(&attach_shipped_sf_backends)`: from then on every
/// `Pipeline` of this process attaches the shipped grids itself.  Loads
/// nothing; the Python module calls it at import.
void install_lhapdf_sf_filler();

}  // namespace lipolgen

#endif  // LIPOLGEN_LHAPDF_SF_HPP
