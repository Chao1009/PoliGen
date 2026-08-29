#ifndef LIPOLGEN_SF_HPP
#define LIPOLGEN_SF_HPP

/// \file sf.hpp
/// Structure-function backends: unpolarized (F2, R, F1, FL), polarized
/// (g1, g2 with a Wandzura-Wilczek default) and tensor (b1, b2, Delta).
///
/// Port of `fastsim/polli_fastsim/structure.py`, `polarized.py` and the toy /
/// moment-ansatz parts of `delta_models.py`.
///
/// Every physics input is an interface with a *toy* implementation always
/// available and a table / LHAPDF implementation optional (docs/CONVENTIONS.md).
/// The digitized theory curves of `polli_fastsim/data/*.csv` are embedded in
/// `src/core/sf_tables.cpp` verbatim and read through `np_interp`, so a C++
/// lookup and a NumPy one return the same bits.
///
/// R = sigma_L/sigma_T is threaded through an `RFunc` hook exactly as the
/// Python threads `r_func`: a null hook means `r_sigma_lt`, the toy default
/// that every published polligen number was made with.

#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "lipolgen/beams.hpp"
#include "lipolgen/constants.hpp"

namespace lipolgen {

// ---------------------------------------------------------------- R hooks

/// `r_func(x, q2) -> R = sigma_L/sigma_T`.  An empty function means "the
/// module default", i.e. `r_sigma_lt`.
using RFunc = std::function<double(double, double)>;

/// A structure-function callable of the (x, q2, f1) signature the tensor
/// slots of `InclusiveKernel` consume.
using SFFunc3 = std::function<double(double, double, double)>;

/// R = sigma_L/sigma_T, simplified R1990-like magnitude (TOY).  The DEFAULT R
/// of the whole program: 0.18/(1 + Q2/50), a placeholder of the right
/// magnitude (0.14-0.18 over the acceptance), not a fit.
double r_sigma_lt(double x, double q2);

inline double resolve_r(const RFunc& f, double x, double q2) {
  return f ? f(x, q2) : r_sigma_lt(x, q2);
}

/// Support of the SLAC/E143 world fit R1998.
inline constexpr double R1998_X_MIN = 0.005;
inline constexpr double R1998_X_MAX = 0.86;
inline constexpr double R1998_Q2_MIN = 0.5;
inline constexpr double R1998_Q2_MAX = 130.0;

/// The three six-parameter forms (R_a, R_b, R_c) of Abe et al. Eq. (2),
/// PLB 452 (1999) 194 (arXiv:hep-ex/9808028), Table II.
struct R1998Forms {
  double a = 0.0, b = 0.0, c = 0.0;
};
R1998Forms r1998_forms(double x, double q2, bool clip = true);

enum class R1998Form { kAverage, kA, kB, kC };

/// R1998 itself is the AVERAGE of the three forms.
double r1998(double x, double q2, R1998Form form = R1998Form::kAverage,
             bool clip = true);
/// max - min of the three forms: the fit's functional-form systematic.
double r1998_spread(double x, double q2, bool clip = true);
/// delta-R of the fit (the unnumbered equation below Abe et al. Eq. 3).
double r1998_fit_error(double x, double q2);

// ------------------------------------------------------- unpolarized SFs

/// Abstract unpolarized backend: F2 on p and n, plus the R-derived F1 and FL.
class UnpolSF {
 public:
  virtual ~UnpolSF() = default;

  virtual double f2p(double x, double q2) const = 0;
  virtual double f2n(double x, double q2) const = 0;
  virtual double f2n_over_f2p(double x) const = 0;

  /// The R this backend converts F2 into F1/FL with.  Null hook = r_sigma_lt.
  void set_r_func(RFunc f) { r_func_ = std::move(f); }
  const RFunc& r_func() const { return r_func_; }
  double r(double x, double q2) const { return resolve_r(r_func_, x, q2); }

  /// F1 via the Callan-Gross relation modified by R (massless):
  /// F1 = F2 / (2 x (1 + R)).
  double f1_from_f2(double f2, double x, double q2) const {
    return f2 / (2.0 * x * (1.0 + r(x, q2)));
  }
  /// FL = F2 R/(1+R).
  double fl_from_f2(double f2, double x, double q2) const {
    const double rr = r(x, q2);
    return f2 * rr / (1.0 + rr);
  }
  double f1p(double x, double q2) const { return f1_from_f2(f2p(x, q2), x, q2); }
  double f1n(double x, double q2) const { return f1_from_f2(f2n(x, q2), x, q2); }
  double flp(double x, double q2) const { return fl_from_f2(f2p(x, q2), x, q2); }
  double fln(double x, double q2) const { return fl_from_f2(f2n(x, q2), x, q2); }

 private:
  RFunc r_func_;
};

/// Crude but smooth F2p / F2n with mild log Q2 evolution (TOY).  Anchored by
/// eye to world ep data: F2p(1e-3, 10) ~ 1.2, F2p(0.1, 10) ~ 0.36,
/// F2p(0.4, 10) ~ 0.12.  Adequate for phase-space maps and factor-1.5 rate
/// estimates ONLY.
class ToyF2 : public UnpolSF {
 public:
  double f2p(double x, double q2) const override;
  double f2n(double x, double q2) const override;
  double f2n_over_f2p(double x) const override;
};

/// Whole-nucleus F2A = Z*F2p + N*F2n with an optional EMC-ratio hook, and the
/// F1A that goes with it.  `emc_ratio(x)` multiplies the isoscalar
/// combination; an empty hook means r = 1 (no medium effect).
class NuclearF2 {
 public:
  NuclearF2(Ion ion, std::shared_ptr<const UnpolSF> base,
            std::function<double(double)> emc_ratio = nullptr,
            RFunc r_func = nullptr);

  double f2a(double x, double q2) const;
  double f1a(double x, double q2) const;

  const Ion& ion() const { return ion_; }
  const std::shared_ptr<const UnpolSF>& base() const { return base_; }
  const RFunc& r_func() const { return r_func_; }

 private:
  Ion ion_;
  std::shared_ptr<const UnpolSF> base_;
  std::function<double(double)> emc_ratio_;
  RFunc r_func_;
};

/// NC DIS double-differential cross section [pb/GeV^2]:
///   d2sigma/dxdQ2 = 4 pi alpha^2/(x Q^4) [(1 - y + y^2/2) F2 - y^2/2 FL],
/// FL = F2 R/(1+R) unless `fl` is given (pass a non-null pointer to override).
double dsigma_dx_dq2(double x, double q2, double s, double f2,
                     const double* fl = nullptr, const RFunc& r_func = nullptr);

// --------------------------------------------------------- polarized SFs

/// Wandzura-Wilczek g2(x) = -g1(x) + int_x^1 du g1(u)/u at fixed Q2.
///
/// Substitution u = x^(1-t) maps the integral to -ln(x) int_0^1 g1(x^(1-t)) dt,
/// evaluated by the trapezoid rule on `npts` points -- the same substitution,
/// the same grid and the same (NumPy pairwise) summation order as
/// `polligen.xsec.g2_ww`, so the two agree bit for bit.
double g2_ww(const std::function<double(double, double)>& g1, double x,
             double q2, int npts = 96);

/// Abstract polarized backend.
class PolSF {
 public:
  virtual ~PolSF() = default;

  virtual double g1p(double x, double q2) const = 0;
  virtual double g1n(double x, double q2) const = 0;

  /// g1A = Z P_p g1p + N P_n g1n, times an optional medium ratio DR(x).
  /// P_p and P_n are PER-NUCLEON effective polarizations, so the nucleon
  /// counts multiply them exactly as `NuclearF2::f2a` multiplies f2p and f2n.
  virtual double g1_nucleus(const Ion& ion, double x, double q2,
                            const std::function<double(double)>& medium_ratio
                                = nullptr) const;

  /// g2 defaults to Wandzura-Wilczek built on this backend's own g1.
  virtual double g2p(double x, double q2, int npts = 96) const;
  virtual double g2n(double x, double q2, int npts = 96) const;
};

/// g1 = A1 * F1 with toy A1(x) shapes (TOY; replace with JAM/DSSV).
class ToyG1 : public PolSF {
 public:
  explicit ToyG1(std::shared_ptr<const UnpolSF> base = nullptr,
                 RFunc r_func = nullptr);

  double a1p(double x) const;
  double a1n(double x) const;
  double g1p(double x, double q2) const override;
  double g1n(double x, double q2) const override;

  const std::shared_ptr<const UnpolSF>& base() const { return base_; }

 private:
  double f1_of(double f2, double x, double q2) const;
  std::shared_ptr<const UnpolSF> base_;
  RFunc r_func_;
};

/// LHAPDF6-backed backends.  DECLARED here, implemented in src/lhapdf/ (P2
/// leaves the implementation out; the core library must not link LHAPDF).
class LhapdfSF : public UnpolSF {
 public:
  explicit LhapdfSF(std::string setname = "CT18NLO", int member = 0);
  ~LhapdfSF() override;
  double f2p(double x, double q2) const override;
  double f2n(double x, double q2) const override;
  double f2n_over_f2p(double x) const override;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

/// LO g1 from a polarized LHAPDF grid (default NNPDFpol11_100).  Declared
/// only; see `LhapdfSF`.
class LhapdfG1 : public PolSF {
 public:
  explicit LhapdfG1(std::string setname = "NNPDFpol11_100", int member = 0);
  ~LhapdfG1() override;
  double g1p(double x, double q2) const override;
  double g1n(double x, double q2) const override;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

// ------------------------------------------------------------ digitized

/// One digitized theory curve set: a strictly increasing x column plus named
/// value columns, embedded verbatim from `polli_fastsim/data/*.csv`.
struct DigitizedTable {
  std::size_t n = 0;
  const double* x = nullptr;
  std::vector<std::string> columns;
  std::vector<const double*> values;

  const double* column(const std::string& name) const;
  /// Table lookup, CONSTANT-EXTRAPOLATED outside the digitized range
  /// (`polarized._interp`).
  double interp(const std::string& name, double xv) const;
  /// Table lookup for a STRUCTURE FUNCTION: constant below the table,
  /// tapered as ((x_end - x)/(x_end - x_max))^power above it
  /// (`polarized._interp_tapered`), so b1/F1 stays frozen instead of
  /// diverging at the generator's top x cell.
  double interp_tapered(const std::string& name, double xv,
                        double x_end = 1.0, double power = 3.0) const;
  double x_min() const { return x[0]; }
  double x_max() const { return x[n - 1]; }
};

namespace tables {
/// Miller PRC 89:045203 Fig. 5, b1(pion) + b1(6q) of the deuteron.
const DigitizedTable& kB1Miller();
/// Cosyn-Dong-Kumano-Sargsian PRD 95:074036 Fig. 4, x*b1 at Q2 = 2.5 GeV2.
const DigitizedTable& kB1CdksQ2p5();
/// Cloet-Bentz-Thomas PLB 642:210 Fig. 6, 7Li at Q2 = 5 GeV2.
const DigitizedTable& kCbtPolemc7liQ5();
/// Tronchin-Matevosyan-Thomas PLB 783:247 Fig. 4, nuclear matter, Q2 = 10.
const DigitizedTable& kTmtPolemcNmQ10();
}  // namespace tables

// -------------------------------------------------------------- EMC hooks

/// Qualitative unpolarized EMC ratio for a light nucleus (SCENARIO): the
/// 12-point canonical shape, linearly interpolated.  This is
/// `polli_fastsim.polarized.unpolarized_emc_ratio(x, mode="table")`, the shape
/// every figure published before 2026-08-29 carried; the data-driven curves
/// (EPPS21, nNNPDF) need nuclear PDFs and live behind the LHAPDF backend, not
/// in the dependency-free core.
double unpolarized_emc_ratio(double x);

/// Digitized 7Li UNPOLARIZED EMC ratio (CBT Fig. 6, blue dashed).
double cbt_unpolarized_emc_ratio(double x);

enum class EmcMode { kDigitized, kConstant };

/// WHICH UNPOLARIZED BASELINE THE POLARIZED-EMC TRANSFER IS REFERENCED TO
/// (`polli_fastsim.polarized.POLEMC_BASELINE`).
///
/// Both published polarized-EMC curves are quoted on a nucleus that is not
/// ours -- CBT computed 7Li, TMT computed nuclear matter -- so each is
/// transferred by the VALENCE SCALE
///     s(table) = <1 - R_unpol,baseline> / <1 - R_unpol,table>
/// averaged over POLEMC_VALENCE_WINDOW, and the transferred curve is
/// 1 - s (1 - R_pol,published).  The scale therefore depends entirely on how
/// deep the baseline's own unpolarized EMC effect is, and that is a CHOICE:
///
///   LegacyTable  the baseline of every number published before 2026-08-29 --
///                CBT's own digitized 7Li R_unpol, <1 - R> = 0.0583.  CBT is
///                then referenced to itself, so its scale is exactly 1, and
///                TMT's is 0.397.
///   Epps21       the default since 2026-08-29 and the default HERE: the
///                EPPS21nlo_CT18Anlo_Li6 / CT18NLO F2 ratio at Q2 = 5,
///                <1 - R> = 0.0298 -- half as deep, so BOTH transferred
///                curves shrink by a factor two (CBT 0.5105, TMT 0.2027).
///
/// The EPPS21 depletion is stored as DATA (`EMC_VALENCE_DEPLETION_EPPS21`)
/// because computing it needs LHAPDF, which the core does not link.
enum class EmcBaseline { LegacyTable, Epps21 };

/// The library default, matching `polli_fastsim.polarized.POLEMC_BASELINE`.
inline constexpr EmcBaseline EMC_BASELINE_DEFAULT = EmcBaseline::Epps21;

/// <1 - R_unpol> of the EPPS21 baseline over POLEMC_VALENCE_WINDOW.
///
/// PROVENANCE: `polli_fastsim.polarized.valence_depletion(mode="epps21")` --
/// the mean of 1 - R over 301 points of x in [0.35, 0.65], with
/// R = F2(EPPS21nlo_CT18Anlo_Li6) / F2(CT18NLO) at Q2 = UNPOL_EMC_Q2 = 5 GeV2
/// through `structure.NuclearF2Ratio`.  Transcribed at full double precision
/// from that call on 2026-08-29; it is the ONLY copy in LiPolGen.
inline constexpr double EMC_VALENCE_DEPLETION_EPPS21 = 0.029788812318099069;

/// <1 - R_unpol> of a baseline over POLEMC_VALENCE_WINDOW.
double emc_valence_depletion(EmcBaseline baseline);

/// Cloet-Bentz-Thomas polarized EMC ratio for 7Li, TRANSFERRED to `baseline`;
/// eq = 23 is the R^{3/2 3/2}_{As} of their Eq. (23), eq = 26 the
/// R^{(3/2 1)}_{As}.
double cbt_polarized_emc_ratio(double x, EmcMode mode = EmcMode::kDigitized,
                               int eq = 23,
                               EmcBaseline baseline = EMC_BASELINE_DEFAULT);
/// Tronchin-Matevosyan-Thomas ratio TRANSFERRED to 7Li valence strength.
double tmt_polarized_emc_ratio(double x, EmcMode mode = EmcMode::kDigitized,
                               EmcBaseline baseline = EMC_BASELINE_DEFAULT);
/// CBT's PUBLISHED 7Li polarized ratio, BEFORE the transfer.
double cbt_published_emc_ratio(double x, int eq = 23);
/// TMT's PUBLISHED nuclear-matter polarized ratio, BEFORE the transfer.
double tmt_published_emc_ratio(double x);
/// Pointwise (1 - R_pol)/(1 - R_unpol) of each camp's OWN figure -- both
/// published curves, so baseline-independent by construction.
double cbt_ratio_of_effects(double x, int eq = 23);
double tmt_ratio_of_effects(double x);

/// <1 - R_unpol,baseline>/<1 - R_unpol,table> over POLEMC_VALENCE_WINDOW, ONE
/// code path for both camps.  On `LegacyTable` CBT is referenced to its own
/// curve, so its scale is exactly 1 and TMT's is 0.397; on `Epps21` they are
/// 0.5105 and 0.2027.
double cbt_valence_scale(EmcBaseline baseline = EMC_BASELINE_DEFAULT);
double tmt_valence_scale(EmcBaseline baseline = EMC_BASELINE_DEFAULT);

// ------------------------------------------------------------- tensor SFs

/// Abstract tensor backend: b1, b2 (default 2x b1) and the double-helicity-
/// flip Delta (default zero).  `f1` is passed through for the scenario models
/// that scale with it; the digitized b1 accessors ignore it.
class TensorSF {
 public:
  virtual ~TensorSF() = default;
  virtual double b1(double x, double q2, double f1) const = 0;
  virtual double b2(double x, double q2, double f1) const {
    return 2.0 * x * b1(x, q2, f1);
  }
  virtual double delta(double x, double q2, double f1) const {
    (void)x; (void)q2; (void)f1;
    return 0.0;
  }
  /// Adaptors to the callable slots of `InclusiveKernel`.  The returned
  /// closures capture `this`; they must not outlive the backend.
  SFFunc3 b1_func() const;
  SFFunc3 b2_func() const;
  SFFunc3 delta_func() const;
};

enum class B1Mode { kDigitized, kToy };

/// The pre-2026-08-28 "HERMES-like" scenario shape times F1.
double toy_b1_shape(double x, double q2, double f1);

/// b1 of the deuteron, per nucleon: Miller's pion + hidden-colour total
/// (PRC 89:045203 Fig. 5), the curve that reproduces HERMES.  Digitized range
/// x = 0.010-0.900; frozen at 0.0647 below, tapered to zero at x = 1 above.
double toy_b1(double x, double q2, double f1, B1Mode mode = B1Mode::kDigitized);

/// b1 of the deuteron, per nucleon: the standard convolution (CDKS SD + DD at
/// Q2 = 2.5 GeV2).  The camp that finds |b1| < 1e-3 at x >~ 0.2, with sign
/// changes at x = 0.06 and 0.42.  Digitized range x = 0.010-1.590.
double b1_convolution(double x, double q2, double f1,
                      B1Mode mode = B1Mode::kDigitized);

/// Integral of b1 over the digitized range -- the Close-Kumano sum rule
/// int b1 dx = 0.  Reported, not enforced.
double close_kumano_integral(bool cdks = true);

/// b1(6Li)/nucleon = transfer * (2/6) * b1(d)/nucleon.  No published 6Li b1
/// exists (plans/04 #9), so the embedded-deuteron scaling is an inference.
double b1_li6_from_deuteron(double b1_d,
                            double transfer = LI6_B1_RANK2_TRANSFER,
                            double per_nucleon = LI6_B1_PER_NUCLEON);

/// TOY double-helicity-flip Delta = scale * F1 * x^0.3 (1-x)^4; `scale` is the
/// peak Delta/F1.  Use 1e-3, 3e-3, 1e-2 as discovery scenarios.
double toy_delta_gluon(double x, double q2, double f1, double scale = 1e-3);

/// Camps as `TensorSF` objects, so a kernel can be handed one backend.
class MillerB1 : public TensorSF {
 public:
  explicit MillerB1(B1Mode mode = B1Mode::kDigitized) : mode_(mode) {}
  double b1(double x, double q2, double f1) const override {
    return toy_b1(x, q2, f1, mode_);
  }

 private:
  B1Mode mode_;
};

class CdksB1 : public TensorSF {
 public:
  explicit CdksB1(B1Mode mode = B1Mode::kDigitized) : mode_(mode) {}
  double b1(double x, double q2, double f1) const override {
    return b1_convolution(x, q2, f1, mode_);
  }

 private:
  B1Mode mode_;
};

/// b1 of 6Li built from a deuteron camp through the rank-2 transfer.
class Li6B1 : public TensorSF {
 public:
  Li6B1(std::shared_ptr<const TensorSF> deuteron,
        double transfer = LI6_B1_RANK2_TRANSFER,
        double per_nucleon = LI6_B1_PER_NUCLEON)
      : d_(std::move(deuteron)), transfer_(transfer), per_nucleon_(per_nucleon) {}
  double b1(double x, double q2, double f1) const override {
    return b1_li6_from_deuteron(d_->b1(x, q2, f1), transfer_, per_nucleon_);
  }

 private:
  std::shared_ptr<const TensorSF> d_;
  double transfer_;
  double per_nucleon_;
};

// ------------------------------------------------------------ Delta models

/// Shape variants of the money_delta convention: {name -> (alpha, beta)}.
struct DeltaVariant {
  double alpha;
  double beta;
};
/// "low_x" (0.3, 4.0), "mid_x" (0.7, 3.0) [production default],
/// "high_x" (1.5, 2.0).
DeltaVariant delta_variant(const std::string& name);

/// Peak value of x^alpha (1-x)^beta at x = alpha/(alpha+beta).
double xab_peak_value(double alpha, double beta);
/// x^a (1-x)^b / peak for the named variant (peak value 1).
double shape_normalized(double x, const std::string& variant);

/// LO analytic alpha_s = 12 pi/[(33 - 2 n_f) ln(Q2/Lambda2)], Lambda = 0.22,
/// n_f = 4 -- the money_delta FALLBACK.  The production convention is the
/// CT18NLO alpha_s table read through `parton`, which has no C++ equivalent
/// in the core library: it makes the moment-constrained Delta ~14 % SMALLER
/// in magnitude than this form.  Pass your own `alphas` to reproduce it.
double alpha_s_lo(double q2, double lambda_qcd = 0.22, double n_f = 4.0);

/// A for Interpretation A: A * int x F1(x, q2_ref) shape(x) dx = c_moment.
/// `f1_func(x, q2)` must return the per-nucleon F1 (F1A/A).  Same grid, same
/// trapezoid and same summation order as `delta_models.solve_A_interp_a`.
double solve_A_interp_a(const std::function<double(double, double)>& f1_func,
                        double q2_ref, const std::string& variant = "mid_x",
                        double c_moment = C_BAG);
/// A for Interpretation B (no F1): the analytic Beta-function integral.
double solve_A_interp_b(const std::string& variant = "mid_x",
                        double c_moment = C_BAG);

/// Callable Delta(x, q2, f1) with metadata for plot annotations.
class DeltaModel : public TensorSF {
 public:
  DeltaModel(std::string name, SFFunc3 func, std::string info)
      : name_(std::move(name)), func_(std::move(func)), info_(std::move(info)) {}
  double b1(double x, double q2, double f1) const override {
    (void)x; (void)q2; (void)f1;
    return 0.0;
  }
  double delta(double x, double q2, double f1) const override {
    return func_(x, q2, f1);
  }
  double operator()(double x, double q2, double f1) const {
    return func_(x, q2, f1);
  }
  const std::string& name() const { return name_; }
  const std::string& info() const { return info_; }

 private:
  std::string name_;
  SFFunc3 func_;
  std::string info_;
};

/// Legacy toy shape; `scale` is the peak Delta/F1.
DeltaModel make_delta_toy(double scale = 1e-3);
/// Interpretation A: Delta = dilution * A * alpha_s(Q2) * F1 * shape.
DeltaModel make_moment_a(const std::function<double(double, double)>& f1_func,
                         double q2_ref, const std::string& variant = "mid_x",
                         double c_moment = C_BAG, double dilution = 1.0,
                         std::function<double(double)> alphas = nullptr);
/// Interpretation B: Delta = dilution * A * alpha_s(Q2) * shape.
DeltaModel make_moment_b(const std::string& variant = "mid_x",
                         double c_moment = C_BAG, double dilution = 1.0,
                         std::function<double(double)> alphas = nullptr);

}  // namespace lipolgen

#endif  // LIPOLGEN_SF_HPP
