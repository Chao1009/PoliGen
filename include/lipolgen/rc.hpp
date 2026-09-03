#ifndef LIPOLGEN_RC_HPP
#define LIPOLGEN_RC_HPP

/// \file rc.hpp
/// Tensor-sector radiative corrections as an OPT-IN, WEIGHT-ONLY family.
///
/// ---------------------------------------------------------------------------
/// RC IS A WEIGHT, AND IT IS NOT EVEN ON THE NOMINAL ONE.
///
/// `fsi.hpp` multiplies `Event::weight` because FSI is a CORRECTION to the
/// model.  These are not: `rc_tensor_lo/hi` is a SYSTEMATIC VARIATION and
/// `rc_tail` is a BACKGROUND, and both must leave the Born sample alone.  They
/// travel in `Event::rc_weights` under their own HepMC3 names and an analysis
/// multiplies them in on purpose.  No four-vector moves, no random number is
/// consumed, and `RcMode::Off` (the default) is today bit for bit.
///
/// ---------------------------------------------------------------------------
/// WHAT IS AND IS NOT CORRECTED (POLRAD 2.0 Eq. (2), sigma = sigma^in +
/// sigma^el + sigma^q + sigma^v; Akushevich, Ilyichev, Shumeiko, Soroko,
/// Tolkachev, arXiv:hep-ph/9706516, CPC 104 (1997) 201).
///
///   sigma^v + sigma^in  the LEPTON-VERTEX correction.  It is spin- and
///     tensor-BLIND: the same radiator multiplies F1, F2 and b1..b4, so the
///     bulk cancels in the ratio A_zz and only an (x, Q^2) shape migration
///     survives.  NOT implemented as a shift -- its residual is exactly what
///     the BAND of `rc_delta` prices.  (A shift would need a photon-energy
///     sampler, i.e. one extra uniform per event, which would break the
///     "RNG untouched" invariant this file rests on.)
///   sigma^el  the ELASTIC radiative tail.  It does NOT cancel -- not because
///     it "carries tensor dependence", but because its own A_zz is not the
///     Born's.  A background with A_zz^tail != A_zz^Born under the DIS bin
///     DILUTES the measured asymmetry; the tail's own (small) tensor part then
///     shifts it back a little.  Implemented as the additive weight `rc_tail`
///     from POLRAD's t-peak closed forms, Eqs. (37)-(39) + (43).
///   sigma^q  the QUASI-ELASTIC tail.  Its POLARISED part is neglected
///     ("there is no net tensor effect by inclusive scattering on weakly-bound
///     spin-1/2 objects", Z.-L. Zhou et al. (NIKHEF), PRL 82 (1999) 687).  Its
///     UNPOLARISED part has A_zz^tail ~ 0, so it does NOT cancel either: it
///     dilutes A_zz exactly as the unpolarised elastic tail does.  HERMES
///     (hep-ex/0506018) subtracted BOTH.  Implemented, unpolarised only, via
///     POLRAD Eq. (44) on the proton elastic tail, folded into `rc_tail`.
///
/// So `rc_tail` is an UNPOLARISED DILUTION with a small tensor correction on
/// top, not "a tensor background": for 6Li the tensor fraction of the elastic
/// tail is ~ 10 % of the deuteron's, because what sets it is Q_A/Z, and
/// Q(6Li)/Q(d) = 0.29 against Z = 3 (design_C_tensor_rc.md sec. 2.1).  It is
/// also QUASI-ELASTIC-DOMINATED at every x >~ 0.03 (sigma^q_U/sigma^el_U =
/// 0.28 / 2.7 / 1000 at x = 0.01 / 0.1 / 0.3), which is why `qe_kf_gev` and
/// `qe_suppression` are the tail's dominant knobs.
///
/// AND IT IS A LOWER BOUND.  `rc_tail` is the t-PEAK ALONE (RcTailModel::
/// TPeak).  T8(c) in tests/test_rc.cpp measures the leading-log s- and
/// p-peaks of the same observable independently: they are negligible for 6Li
/// at EIC Q^2 >= 20 GeV^2 (POLRAD sec. 2.1.3 B confirmed WHERE THIS GENERATOR
/// RUNS) and they are NOT elsewhere -- at the HERMES deuteron point the
/// t-peak is 23 % of the leading-log total, and at the Q^2 ~ 4 GeV^2 corner
/// of the generator window the quasi-elastic s-peak is already 3.3x it.
/// Never quote `rc_tail` as "the" radiative tail, and never port it to
/// fixed-target kinematics without the s-/p-peaks.
///
/// PER-NUCLEON, and the factor is 1/A^2.  Eq. (38) is the WHOLE-NUCLEUS
/// d^2 sigma/(dx_A dy) (Weizsacker-Williams x Compton reproduces it with no
/// 1/A in it, and POLRAD's FORTRAN applies BOTH `ter = m_p/M_A` in `apptai`
/// and `/tara` in the main program), so per nucleon is (1/A) x the Jacobian
/// dx_A/dx = 1/A.  See src/core/rc.cpp's "PER-NUCLEON REDUCTION" block and
/// polrad_transcription_check.md sec. 8b.
///
/// ---------------------------------------------------------------------------
/// THE ONE DEFINITION OF "THE TENSOR PART OF AN EVENT".
///
/// tau  ==  W_tensor(phi') / W(phi')
///
/// with `W` the azimuthal density of `xsec.hpp` and `W_tensor` its b1..b4
/// (T_LL) contribution ALONE -- exactly what `InclusiveKernel::tensor_
/// amplitudes` returns and exactly what `InclusiveSampler::StateTables::
/// w_tensor/a1_tensor/a2_tensor` hold.  At an UNPOLARISED beam and along it
/// this reduces to the supervisor's / HERMES's own coefficient
///
///   tau = (P_zz/2) A_zz / [1 + (P_zz/2) A_zz],   P_zz^eff = 3 Q_NN P_2(cos theta_S)
///
/// (HERMES Eq. (1)), and the SAME closed form covers the tagged channels with
/// A_zz reading the WAVE-FUNCTION asymmetry `azz_tensor_curve` instead of the
/// DIS one -- see `RcModel::tensor_fraction`.  The band is then
///
///   w_+- = 1 +- delta(x) * tau
///        = [1 + (P_zz/2) A_zz (1 +- delta)] / [1 + (P_zz/2) A_zz] .
///
/// NOTE `1 - rhobar/rho` is a CLOSED-FORM SHORTCUT for tau, valid only at
/// lam_e * P_e = 0; it is "everything of rank >= 1" and with a polarised beam
/// it would sweep in the vector A_par (and a1 A_perp) terms this family
/// deliberately excludes.  `RcModel`'s constructor THROWS on a plan with
/// lam_e * P_e != 0 rather than pricing the vector sector as if it were rank 2.
///
/// ---------------------------------------------------------------------------
/// PER-CHANNEL RULE (design_C_tensor_rc.md sec. 1.5; total over `Channel`):
///
///   Inclusive     band AND tail (the tail at EVERY theta_S, through
///                 Q_N -> P_zz^eff -- POLRAD Eq. (43) gives the transverse
///                 case as -1/2 the longitudinal one, i.e. P_2(cos 90 deg))
///   every Tagged* band only; `rc_tail` == 1 exactly, because the elastic
///                 recoil sits at x_L = 1 inside the 10-sigma beam-exclusion
///                 envelope while the tag looks at x_L ~ A_spec/A_beam, so the
///                 tag itself vetoes it.  (The QUASI-elastic half of that
///                 statement is an assumption, not a kinematic fact -- Q3.)
///   CoherentLi6   neither; both weights exactly 1.0 and the run prints why.
///                 Its tensor dependence is entirely AZIMUTHAL (cos 2phi) and
///                 nobody has computed RC for a phi-dependent tensor
///                 observable; and its own M_X >= 1 GeV cut already excludes
///                 the elastic point.
///
/// ---------------------------------------------------------------------------
/// HONEST FLAGS -- repeat them wherever a number from this file is quoted.
///
///   * `RC_DELTA_LOW_X` is the SIZE of a correction this generator does not
///     apply, taken as a 1-sigma band.  It is a CHOICE, and its source has
///     ZERO INSPIRE citations.  BAND IT: run 0.19 and 0.30, never quote one.
///   * No RC calculation exists for a TAGGED tensor asymmetry.  Applying
///     delta(x) to tau_tag is a defensible, conservative EXTRAPOLATION and
///     must never be quoted as a published result.
///   * No RC calculation exists for any PHI-DEPENDENT tensor observable (the
///     Delta cos 2phi sector, the coherent channel) at ANY axis.
///   * Everything cited is DEUTERON.  Whether the deuteron's fractional RC
///     transfers to 6Li is untested, and is the largest unquantified
///     assumption in the band.
///   * On the TAGGED channels tau_tag = 1 - nbar/n_M is UNBOUNDED at the nodes
///     of the M-dependent spectator density, and `RcOptions::band_tau_max`
///     CLAMPS it there rather than publishing negative band edges.  The clamp
///     is a CHOICE; n_M -> 0 is exactly where the fractional-rescale ansatz
///     breaks down (the tensor part of the density cancels the unpolarised
///     part), and the clipped fraction is reported per run and per event
///     (`Event::rc_clipped`) rather than hidden.
///   * `rc_tail` is the t-PEAK ONLY and is therefore a LOWER BOUND on the
///     dilution -- see the RcTailModel::TPeak comment and T8(c).
///
/// Design: docs/open_items/run_2026-09-02/design_C_tensor_rc.md.

#include <array>
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "lipolgen/constants.hpp"
#include "lipolgen/event.hpp"
#include "lipolgen/sampler.hpp"
#include "lipolgen/tagged.hpp"

namespace lipolgen {

// ------------------------------------------------------------- the band

/// PR12-13-011 (E12-13-011 PROPOSAL, unpublished): "no polarized radiative
/// corrections at the lepton vertex, and the unpolarized corrections are known
/// to better than 1.5 %".  `RC_X_HIGH` is that experiment's OWN lower x edge,
/// 0.16 (arXiv:2506.04506 p. 8: 0.16 < x < 0.49, 0.8 < Q^2 < 5.0 GeV^2).
/// Quoting 1.5 % below x ~ 0.1 is unsupported and contradicts HERMES's
/// MEASURED 15 % residual at x = 0.063.
inline constexpr double RC_DELTA_HIGH_X = 0.015;
inline constexpr double RC_X_HIGH       = 0.16;
/// Gakh-Shekhovtsova (hep-ph/0403262, JETP 99 (2004) 898) SIZE of the RC on
/// the spin-dependent cross section at x ~ 1e-3 - 1e-2: "the value of
/// radiative correction changes from 10 % to 30 % as compared with the Born
/// contribution", said of an unpolarised beam on a tensor-polarised deuteron
/// -- exactly this configuration -- in HERMES kinematics.  ZERO INSPIRE
/// CITATIONS: a single unchecked calculation, and its two numbers are the
/// spread of ONE calculation over Q^2, not two independent edges.
/// THIS IS A CORRECTION MAGNITUDE, NOT A MEASURED RESIDUAL: taking it as a
/// 1-sigma band is a CHOICE this generator makes because it does not apply the
/// correction.  Default = the conservative 0.30 edge, which CONTAINS the 0.10
/// one, so nothing is lost by it.
inline constexpr double RC_DELTA_LOW_X            = 0.30;
/// The residual HERMES actually ACHIEVED at the lowest x: 2e-3 against a
/// MEASURED |A_zz| of 1.06e-2 at <x> = 0.012 (hep-ex/0506018 Table II) =>
/// 0.19 fractional (0.19, 0.19, 0.15 in the three low-x bins).  NOT 0.10 --
/// that came from dividing 2e-3 by the RANGE MAXIMUM |A_zz| <= 0.02, which
/// occurs at x = 0.45 where HERMES says RC is negligible.
inline constexpr double RC_DELTA_LOW_X_OPTIMISTIC = 0.19;
inline constexpr double RC_X_LOW                  = 0.01;

/// The ceiling on |tau| the BAND is allowed to see, i.e. `RcOptions::
/// band_tau_max`'s default.  1.0 means "the band never varies more than the
/// whole tensor part of the event", which keeps both edges w_+- = 1 +- delta
/// tau inside [1 - delta, 1 + delta] and therefore STRICTLY POSITIVE for any
/// delta < 1.  It is a CHOICE, and it bites only on the tagged channels: there
/// tau_tag = 1 - nbar/n_M diverges wherever the event's own n_M(k, c) is near
/// a node (the 6Li M = 0 density has them), and without the clamp a 20 k-event
/// tagged-alpha run produces rc_tensor_hi down to -1.8 (7Li: -8.3).  A
/// Monte-Carlo weight must be bounded -- the same rule `RcOptions::tail_max`
/// applies to the tail.
inline constexpr double RC_BAND_TAU_MAX = 1.0;

/// The tail tables' hard y ceiling.  The t-peak carries
/// Y_+ = [1 + (1-y)^2]/(1-y) ~ 1/(1-y), which diverges as y -> 1; the tables
/// stop here and `RcModel`'s constructor REFUSES a scenario whose y_max is
/// above it, rather than letting `tail_ratio_at` throw out of the event loop
/// on the first event past the edge.
inline constexpr double RC_TAIL_Y_CEILING = 0.9995;

/// Fermi momentum of 6Li, E. J. Moniz et al., PRL 26 (1971) 445 (the
/// quasi-elastic e-A scaling fit; POLRAD's own target block carries 0.221 for
/// C and 0.164 for 3He, adgh:696-714).  It is `RcOptions::qe_kf_gev`'s
/// default and drives the de Forest-Walecka Pauli suppression S(q) of the
/// QUASI-ELASTIC tail -- which is NOT a small choice: the t-peak reaches down
/// to t_min ~ (x M_N)^2, so at x <= 0.1 most of the QRT integral sits below
/// 2 k_F.  Set it to 0 to recover POLRAD Eq. (44) at S_E = S_M = 1.
inline constexpr double RC_QE_KF_GEV = 0.169;

/// The band half-width delta(x): log-linear interpolation between the two
/// anchors, clamped outside them.
///
///                { delta_high                                       x >= x_high
///   delta(x)  =  { delta_high + (delta_low - delta_high)
///                {              * ln(x_high/x) / ln(x_high/x_low)
///                { delta_low                                        x <= x_low
///
/// Log-linear in x because both the anchors and the underlying physics (the
/// soft-photon ln(1/x) and the growth of the radiative background towards the
/// kinematic edge) are logarithmic in x.  Monotone NON-INCREASING in x, with
/// both anchor values returned EXACTLY (the clamp branches return the constant
/// unmodified, which is why `tests/test_rc.cpp` may compare them with ==).
///
/// At the defaults it passes through delta(0.063) = 0.1108, i.e. it does not
/// contradict HERMES's measured 15 % residual there.
///
/// THROWS `std::runtime_error` unless 0 < x_low < x_high: the interpolation
/// divides by ln(x_high/x_low) and would otherwise return a silent NaN.
double rc_delta(double x,
                double delta_high = RC_DELTA_HIGH_X,
                double delta_low  = RC_DELTA_LOW_X,
                double x_high     = RC_X_HIGH,
                double x_low      = RC_X_LOW);

// --------------------------------------------- spin-1 elastic form factors

/// POLRAD Eq. (A.4)'s (F_c, F_m, F_q) at the elastic-vertex t [GeV^2].
/// Normalisation is the Rosenbluth one, fixed by Eq. (A.4)'s Q_N = 0 limit
/// (Im^el_2 = A(Q^2), Im^el_1 = B(Q^2)/2 with the standard spin-1
/// A = G_C^2 + (8/9) eta^2 G_Q^2 + (2/3) eta G_M^2, B = (4/3) eta(1+eta) G_M^2):
///
///   F_c(0) = Z,   F_m(0) = (M_A/m_p) mu_A/mu_N,   F_q(0) = M_A^2 Q_A .
///
/// Validated on the DEUTERON with this repository's own constants
/// (M_d = 1.875612942 GeV, HBARC_GEV_FM): G_Q(0) = 25.829, G_M(0) = 1.7139 --
/// the textbook 25.83 and 1.714.  The convention is therefore fixed with no
/// free choice.
class Spin1ElasticFF {
 public:
  virtual ~Spin1ElasticFF() = default;
  virtual double fc(double t_gev2) const = 0;
  virtual double fm(double t_gev2) const = 0;
  virtual double fq(double t_gev2) const = 0;
  virtual std::string provenance() const = 0;   ///< printed by the run
};

/// The v0 6Li model: a two-parameter harmonic-oscillator point-nucleon shape
/// for the charge and quadrupole multipoles, a separate two-parameter shape
/// for the magnetic one, all normalised on the MEASURED moments (NOT on VMC --
/// Wiringa & Schiavilla (nucl-th/9807037) get Q(6Li) = -0.23(9) fm^2 against
/// the measured -0.0818, and `docs/OPEN_ITEMS_SOLUTIONS.md` sec. 5 already
/// forbids deriving a 6Li tensor input from those wave functions):
///
///   F_point(q) = [1 - (alpha/(2+3 alpha)) (q^2 a^2/2)] exp(-q^2 a^2/4)
///   F_mag(q)   = (1 - q^2/q_z^2) exp(-q^2 b^2/4)
///   F_c(t) =  Z       F_point(q) [G_E^p(t) + G_E^n(t)]
///   F_q(t) =  M_A^2 Q F_point(q) [G_E^p(t) + G_E^n(t)] * fq_scale
///   F_m(t) = (M_A/m_p)(mu/mu_N) F_mag(q) * tail_tensor_scale
///
/// with q^2 = t/(hbar c)^2.  `F_m` gets its OWN shape because WS98's F_T has
/// its first peak at q = 0.5 fm^-1, a zero, and a second peak at q = 2 fm^-1 --
/// i.e. structure exactly where the tail lives (q = 2 fm^-1 <=> t ~ 0.155
/// GeV^2) -- and a single monopole has none there.  The nucleon folding uses
/// G_E^p + G_E^n (N = Z = 3, isoscalar), not G_E^p alone, because <r^2>_point
/// is derived by subtracting BOTH <r^2>_p and (N/Z)<r^2>_n.
///
/// NO ION-SPECIFIC DEFAULTS.  Every `0.0` below means "take it from the Ion";
/// `HoSpin1FF::for_ion` fills them.  Retyping M(6Li) or Z here would duplicate
/// the `beams.cpp` mass table and `LI6().Z` (docs/CONVENTIONS.md: no physics
/// number in two places), and an ion-specific DEFAULT would silently give a
/// 7Li or deuteron run 6Li form factors.
struct HoSpin1FFOptions {
  double a_fm     = 0.0;   ///< 0 => from the ion's (LI6_FF_HO_*) block
  double alpha    = 0.0;
  double z        = 0.0;   ///< 0 => ion.Z
  double m_a_gev  = 0.0;   ///< 0 => ion.mass()  (the AME mass, NOT A*M_NUCLEON)
  double mu_n     = 0.0;   ///< 0 => the ion's measured magnetic moment
  double q_fm2    = 0.0;   ///< 0 => the ion's measured quadrupole moment
  /// F_m's OWN shape (WS98 F_T: first peak 0.5, zero, second peak 2 fm^-1).
  /// F_mag(q) = (1 - q^2/q_z^2) exp(-q^2 b^2/4).  A shared monopole shape has
  /// no structure there, which is where the tail lives.
  double fm_qz_fm = 0.0;
  double fm_b_fm  = 0.0;
  /// The +-100 % quadrupole band.  sigma^el_T is QUADRATIC in F_q (Eq. (38)'s
  /// sigma_q^d carries F_q(3F_c + 3 eta F_m + eta F_q) and
  /// F_q(4F_c - 3x F_m + (4/3) eta F_q)), so this must be RUN (0, 1, 2), never
  /// rescaled from one run.  See T12.
  double fq_scale = 1.0;
  /// Flat multiplier on F_m -- the eta F_m^2 tensor sector, which `fq_scale`
  /// does NOT span.  Run 0.5, 1, 2.
  double tail_tensor_scale = 1.0;
  bool   fold_nucleon = true;   ///< multiply by (G_E^p(t) + G_E^n(t)), N = Z
};

/// Starting parameters of the 6Li `F_point` shape, TO BE REFIT by the
/// implementer against the Suelzle-Yearian-Crannell (PR 162 (1967) 992) and
/// Li-Sick-Whitney-Yearian (NPA 162 (1971) 583) elastic data.  They are the
/// (a, alpha) pair satisfying (i) <r^2>_point = 6.078 fm^2 -- from
/// r_ch(6Li) = 2.589(39) fm (Angeli & Marinova, ADNDT 99 (2013) 69) minus
/// <r^2>_p = 0.7071, plus -<r^2>_n = 0.1155, minus the Darwin-Foldy 0.033 fm^2
/// -- and (ii) a first C0 zero at q_0 = 3.1 fm^-1.
///
/// q_0 = 3.1 fm^-1 IS A STARTING GUESS, NOT A SOURCE.  WS98 does not locate
/// the zero; it says only that C2 is "much smaller than C0 below 3 fm^-1",
/// that "for q >= 3 fm^-1 the C2 contribution becomes dominant", and that
/// two-body currents "shift the minimum to lower values of q" with no number.
/// Fit q_0 over [2.9, 3.3] fm^-1 and pin the REFIT result, with the refit's
/// own uncertainty (T11).  Note the shell-model alpha = (Z-2)/3 = 1/3 puts the
/// zero at 2.33 fm^-1, far too low, so alpha must be FREE: this is a
/// phenomenological fit and must be labelled as such, not a shell model.
/// (These two reproduce both targets exactly:
///  <r^2>_point = (3/2) a^2 (2+5 alpha)/(2+3 alpha) = 6.0788 fm^2 and
///  q_0 = (1/a) sqrt(2(2+3 alpha)/alpha) = 3.0998 fm^-1.)
inline constexpr double LI6_FF_HO_A_FM   = 1.9069;
inline constexpr double LI6_FF_HO_ALPHA  = 0.13822;

/// `F_mag`'s zero and width [fm^-1, fm].  UNFITTED STARTING VALUES -- the fit
/// itself is the implementer's (T11, and the "fitted (q_z, b) of F_mag against
/// WS98's three F_T landmarks" row of design sec. 8.3).
///
/// `q_z` is the midpoint of the only range the design states, 1.2-1.4 fm^-1.
/// `b` is a placeholder, and deliberately not quoted to more digits than it
/// deserves: a TWO-parameter shape cannot hit WS98's THREE F_T landmarks, and
/// with F_T ~ q F_m at q_z = 1.30 the two single-landmark solutions of
/// d/dq [q(1 - q^2/q_z^2) exp(-q^2 b^2/4)] = 0 are b = 2.285 fm (first peak at
/// q = 0.5 fm^-1) and b = 1.494 fm (second peak at q = 2 fm^-1); 1.85 is their
/// geometric mean.  Q10: `tail_tensor_scale` bands the NORMALISATION of F_m,
/// not its SHAPE, so this pair is a real model uncertainty until WS98 Fig. 1's
/// F_T is digitised.
inline constexpr double LI6_FF_FM_QZ_FM  = 1.30;
inline constexpr double LI6_FF_FM_B_FM   = 1.85;
/// mu(6Li) = +0.8220473 mu_N -- measured; TUNL A = 6 evaluation.
inline constexpr double LI6_MU_N          = 0.8220473;
/// Q(6Li) = -0.0818(17) fm^2 -- measured; the TUNL A = 6 evaluation prints
/// Q = -0.818(17) mb (1998CE04).  The often-quoted -0.0806(6) fm^2 is the
/// PYYKKO compilation value (Mol. Phys. 106 (2008) 1965, from Cederberg et al.
/// 1998), NOT TUNL's; switching to it is a one-constant change and moves
/// F_q(0) by 1.5 % (-65.914 -> -64.947).
inline constexpr double LI6_QUADRUPOLE_FM2 = -0.0818;
/// Together with `LI6().Z` and `Ion::mass()` these two fix the T11 anchors,
/// which are DERIVED and never retyped:
///   F_c(0) = Z                                          = 3
///   F_m(0) = (M_A/PROTON_MASS) * LI6_MU_N               = 4.90765
///   F_q(0) = (M_A/HBARC_GEV_FM)^2 * LI6_QUADRUPOLE_FM2  = -65.914
/// (M_A = 5.601518702 GeV = 28.38655 fm^-1, the `beams.cpp` AME mass, NOT
/// A*M_NUCLEON -- which would be 5.6296 GeV, 0.5 % high, and would move eta_A
/// by 1 %.)

/// `HoSpin1FFOptions` evaluated on one ion.  BODIES LAND WITH THE MODEL --
/// this file currently defines `rc_delta` and the name helpers only.
class HoSpin1FF : public Spin1ElasticFF {
 public:
  /// The ONLY supported constructor path.  Reads `ion.mass()` and `ion.Z`,
  /// fills every zero field from the ion's measured-moment block, and THROWS
  /// for any ion that has no such block (today: everything but 6Li).
  static std::shared_ptr<HoSpin1FF> for_ion(const Ion& ion,
                                            HoSpin1FFOptions opt = {});

  double fc(double t_gev2) const override;
  double fm(double t_gev2) const override;
  double fq(double t_gev2) const override;
  std::string provenance() const override;

  const HoSpin1FFOptions& options() const { return opt_; }

 private:
  explicit HoSpin1FF(HoSpin1FFOptions opt) : opt_(opt) {}

  HoSpin1FFOptions opt_;
};

/// A digitised (q, F_C0, F_C2, F_M1) table under `data/ff/`, log-linear in q,
/// refusing to extrapolate.  `fq_scale` and `tail_tensor_scale` apply here too.
///
/// The multipole <-> Sachs conversion every source has to be read through is
///   |F_L(q)|^2 = F_C0^2 + F_C2^2,  F_C0 = F_c,  F_C2 = (2 sqrt(2)/3) eta_A F_q
/// which reproduces A(Q^2) = G_C^2 + (8/9) eta^2 G_Q^2 + (2/3) eta G_M^2 term
/// by term.  NOTE no experiment has separated F_C2 for 6Li: the elastic
/// longitudinal data measure F_C0^2 + F_C2^2 and C2 is sub-dominant exactly
/// where those data are good, so the C2 has two checkable ends (the exact
/// low-q limit set by the MEASURED quadrupole moment, and the q >~ 3 fm^-1
/// shoulder WS98 says is entirely C2) and model in between.  Do NOT rescale
/// the whole WS98 C2 curve by ~0.35: that fixes q -> 0 and breaks the high-q
/// end, which is the one place the VMC C2 magnitude is right.
///
/// The path is resolved under `data_dir()` (cluster.hpp -- $LIPOLGEN_DATA_DIR,
/// else the compiled-in prefix), exactly like the VMC tables.
/// docs/CONVENTIONS.md: "Nothing resolves a data path against the working
/// directory."  BODIES LAND WITH THE MODEL.
class TabulatedSpin1FF : public Spin1ElasticFF {
 public:
  /// `relative` is relative to `data_dir()`, e.g. "ff/li6_elastic.csv".
  static std::shared_ptr<TabulatedSpin1FF> from_data_dir(
      const std::string& relative, HoSpin1FFOptions norm);

  double fc(double t_gev2) const override;
  double fm(double t_gev2) const override;
  double fq(double t_gev2) const override;
  std::string provenance() const override;

 private:
  TabulatedSpin1FF() = default;

  std::string path_;
  std::vector<double> q_fm_, c0_, c2_, m1_;
  HoSpin1FFOptions norm_;
};

// ------------------------------------------------ nucleon elastic form factors

/// Proton and neutron magnetic moments [mu_N] and the dipole mass, used ONLY
/// by `nucleon_ff` below.  They are RC-owned (docs/CONVENTIONS.md: a constant
/// lives in the header that owns it) and appear nowhere else in the library.
inline constexpr double MU_PROTON  =  2.7928473;
inline constexpr double MU_NEUTRON = -1.9130427;
/// The dipole mass squared [GeV^2] of G_D(t) = 1/(1 + t/0.71)^2.
inline constexpr double NUCLEON_FF_DIPOLE_GEV2 = 0.71;
/// Galster's G_E^n(t) = -mu_n tau G_D/(1 + a tau), tau = t/4 m_p^2, a = 5.6
/// (S. Galster et al., Nucl. Phys. B32 (1971) 221).
inline constexpr double GALSTER_A = 5.6;

/// The four nucleon Sachs form factors at spacelike t = Q^2 [GeV^2].
struct NucleonFF {
  double ge_p = 0.0, gm_p = 0.0, ge_n = 0.0, gm_n = 0.0;
};

/// Dipole G_E^p, G_M^p/mu_p, G_M^n/mu_n plus Galster's G_E^n.  Two consumers,
/// and they are the reason it is a free function rather than a lambda inside
/// one of them:
///   * the ISOSCALAR NUCLEON FOLDING (G_E^p + G_E^n) of `HoSpin1FF`, which is
///     what turns a POINT-nucleon 6Li density into a charge form factor;
///   * the QUASI-ELASTIC tail of POLRAD Eq. (44), which is the proton elastic
///     tail evaluated once with (G_E^p, G_M^p) and once with (G_E^n, G_M^n).
/// A 5 % dipole is adequate for BOTH: the folding is a small correction on a
/// nuclear form factor that is itself a phenomenological fit, and the QRT is
/// carried with a flat +-100 % `RcOptions::qe_suppression` band anyway.
NucleonFF nucleon_ff(double t_gev2);

// ------------------------------------------- the eta_A integration limits

/// The eta_A = t/(4 M_A^2) window of POLRAD's t-peak quadrature.
struct EtaLimits {
  double lo = 0.0, hi = 0.0;
};

/// POLRAD's OWN EXACT kinematic limits, in the nuclear invariants of design
/// sec. 1.4.1 (`x_a` = x_A, `s_a` = S_A = 2 k1 . p_A, `m_a` = M_A):
///
///   S_x = y S_A,   Q^2 = x_A y S_A,   lambda_Q = S_x^2 + 4 M_A^2 Q^2
///   eta_lo = [ (S_x - Q^2)(S_x - sqrt(lambda_Q)) + 2 M_A^2 Q^2 ]
///            / [ 8 M_A^2 (S_x - Q^2 + M_A^2) ]
///   eta_hi = S_x / (4 M_A^2)
///
/// This is what the POLRAD FORTRAN integrates (deck `apptai`); the paper's
/// eta_min = x_A^2/(4(1-x_A)) and eta_max = infinity of Eq. (39) are the
/// ULTRARELATIVISTIC limit of it and are kept in the shipped code only as a
/// commented-out line.  The two lower limits agree to 0.04 % at x = 0.01 and
/// differ by 8 % at x = 0.3, and since the integrand is proportional to
/// X~, which VANISHES at the u.r. eta_min and does NOT at the exact one, the
/// difference is not a rounding detail at high x.
EtaLimits polrad_eta_limits(double x_a, double y, double s_a, double m_a);

/// The paper's ultrarelativistic lower limit, POLRAD Eq. (39):
/// eta_min = x_A^2/(4(1 - x_A)), i.e. t_min = 4 M_A^2 eta_min
/// = M_A^2 x_A^2/(1 - x_A) ~= (x m_p)^2 -- INDEPENDENT of A, which is what
/// makes the t-peak reach inside the nuclear form factor at every (x, Q^2).
/// Kept as the T8 cross-check on `polrad_eta_limits`.
double polrad_eta_min_ur(double x_a);

// ------------------------------------------ the Eq. (38) t-peak quadratures
//
// PUBLIC, and not private helpers of `RcModel`, for one reason: design sec. 5
// T8 gates them against a LITERAL transcription of Eq. (38) written IN
// `tests/test_rc.cpp` at 1e-10.  That gate is what catches a missing A = 6, a
// swapped (3/4) <-> (4/3), or a dropped Q_N/6 -- none of which any tolerance
// test on `w_tail` would see -- and it can only be written if the test can
// drive the same quadrature with its own integrand.

/// POLRAD's t-peak quadrature, Eq. (38), written ONCE:
///
///   sigma = - (alpha^3 / S_A) Y_+(y) INT_lo^hi (d eta_A / eta_A)
///                                       integrand(eta_A, Xt, X_1)
///
///   X_1 = x_A^2 + 4 x_A eta_A - 4 eta_A      (POLRAD Eq. (39); the last term
///                                             is LINEAR in eta_A -- CONFIRMED
///                                             against polrad2t.tex:936 and
///                                             adgh:8887/8930/8980)
///   Xt  = X_1 / (2 eta_A x_A^2)
///   Y_+ = (1 + (1-y)^2) / (1 - y)
///
/// TWO THINGS THAT ARE NOT IN EQ. (38) AS PRINTED, both from
/// docs/open_items/run_2026-09-02/polrad_transcription_check.md:
///
///   * THE LEADING MINUS (check sec. 8).  X_1 is decreasing in eta_A and
///     vanishes at the ultrarelativistic eta_min, so Xt < 0 over essentially
///     the whole range and Eq. (38) as printed returns a NEGATIVE sigma_u --
///     a tail that REMOVES events from the DIS bin.  Eq. (18) (polrad2t.tex
///     line 580) carries an explicit leading minus that the ultrarelativistic
///     Eq. (38) does not print; with it, sigma^el_U > 0 and w_tail >= 1.
///     T8(0) is the gate.
///   * Y_+ EVERYWHERE (check sec. 1).  The paper prints Y_- on sigma_u^C
///     (polrad2t.tex:928); POLRAD's own `apptai` applies Y_+ to every
///     UNPOLARISED and every TENSOR entry and Y_- only to the vector one
///     (adgh:8614-8648).  Y_- = y(2-y)/(1-y) vanishes as y -> 0, which an
///     unpolarised cross section cannot.  The printed Y_- is a typo.
///
/// The integrand is passed `(eta_A, Xt, X_1)` because Eq. (38)'s sigma_q^d
/// needs all three.  Result in GeV^-2 (it is a d^2 sigma / dx_A dy), PER
/// NUCLEUS in the sense of Eq. (18)'s left-hand side -- i.e. it already
/// carries Eq. (18)'s own 1/A; see `RcModel`'s per-nucleon comment.
///
/// The quadrature is COMPOSITE Gauss-Legendre in ln(eta_A) -- `n_eta` nodes in
/// panels of 8 -- because the integrand lives in the few e-folds of ln(eta_A)
/// where the form factor is alive while the kinematic range spans fifteen.
double polrad_tpeak_quadrature(
    const std::function<double(double eta, double xt, double x1)>& integrand,
    double x_a, double y, double s_a, EtaLimits lim, int n_eta = 128);

/// POLRAD Eq. (38) `sigma_u^d` -- the UNPOLARISED elastic tail of a SPIN-1
/// nucleus, in nuclear invariants, GeV^-2:
///
///   integrand = ( F_c^2 + (8/9) eta_A^2 F_q^2 + (2/3) eta_A F_m^2 ) Xt
///               - (2/3) (1 + eta_A) F_m^2
///
/// = `Xt Im^el_2|_{Q_N=0} - Im^el_1|_{Q_N=0}/eta_A` (check sec. 1's master
/// formula), so a spin-0 nucleus (F_m = F_q = 0, F_c = Z F) reduces it to
/// `Z^2 F^2 Xt` -- the Eq. (38) carbon line, with the CORRECTED Y_+.  That is
/// T8(b), and it is why there is no separate spin-0 entry point.
double polrad_sigma_el_u(const Spin1ElasticFF& ff, double x_a, double y,
                         double s_a, double m_a, int n_eta = 128);

/// POLRAD Eq. (38) `sigma_q^d` -- the QUADRUPOLARISED elastic tail, the
/// partner Eq. (37) carries as `Q_N/6 * sigma_q^A`.  In nuclear invariants,
/// GeV^-2.  The two fractions are `(3/4) x_A^2` inside the Xt bracket and
/// `(4/3) eta_A F_q` in the last line -- CONFIRMED against polrad2t.tex:915-926
/// and adgh:8985-8987; the PDF renders them SWAPPED, so never transcribe this
/// from the PDF.
///
/// It is a SEPARATE TABLE from `polrad_sigma_el_u` and never a scale factor on
/// it: `sigma_q/sigma_u` changes SIGN between x = 0.05 and x = 0.20
/// (check sec. 8: +0.106, +0.064, -0.117 at three deuteron points).
double polrad_sigma_el_t(const Spin1ElasticFF& ff, double x_a, double y,
                         double s_a, double m_a, int n_eta = 128);

/// POLRAD Eq. (44), the UNPOLARISED quasi-elastic tail at S_E = S_M = S_EM = 1:
///
///   sigma_u^q(A) = Z sigma_u^p[G_E^p, G_M^p] + N sigma_u^p[G_E^n, G_M^n]
///
/// with `sigma_u^p` the spin-1/2 line of Eq. (38) in NUCLEON invariants
/// (M -> m_N, x_A -> x, S_A -> s) and Eq. (40) / Eq. (A.5):
///
///   integrand = (G_E^2 + eta G_M^2)/(1 + eta) * Xt - G_M^2
///             = (F_1^2 + eta F_2^2) Xt - (F_1 + F_2)^2 .
///
/// `kf_gev > 0` turns on POLRAD Eq. (44)'s S_E/S_M -- the de Forest-Walecka
/// Fermi-gas Pauli suppression `ffquas` (adgh:5605-5613) codes as
///     S(q) = (3/4)(q/k_F) - (1/16)(q/k_F)^3   (q < 2 k_F),   1 otherwise
/// at the elastic-vertex three-momentum transfer q^2 = t(1 + eta).  Since the
/// integrand is a combination of G_E^2 and G_M^2 alone (no G_E G_M
/// interference), S_E = S_M = S multiplies it as a whole.  `kf_gev = 0` (this
/// function's default, NOT `RcOptions`') is the unsuppressed S = 1 form.
///
/// PER NUCLEUS (the Z and N multiplicities are in it); `RcModel` divides by A.
/// The POLARISED quasi-elastic tail is NOT here and is not implemented at all
/// -- Z.-L. Zhou et al., PRL 82 (1999) 687.
double polrad_sigma_qe_u(int z, int n, double x, double y, double s,
                         double m_n, int n_eta = 128, double kf_gev = 0.0);

/// The de Forest-Walecka Fermi-gas Pauli suppression factor, exposed so a test
/// and a plot can see it: S(q) = (3/4)u - u^3/16 for u = q/k_F < 2, else 1.
/// Continuous at u = 2 (3/2 - 1/2 = 1) and 0 at q = 0.
double pauli_suppression(double q_gev, double kf_gev);

// ------------------------------------------------------------- the model

enum class RcMode : int {
  Off        = 0,   ///< today, bit for bit; every weight identically 1.0
  TensorBand = 1,   ///< the band + the radiative tails
};

/// Which tail formulation.  v0 ships `TPeak` only; `PolradFull` is the
/// documented upgrade path (design_C_tensor_rc.md sec. 1.4.6).
enum class RcTailModel : int {
  TPeak      = 0,   ///< DEFAULT: POLRAD Eqs. (37)-(39), (43) -- one eta_A
                    ///< integral.  IT IS ONE PEAK OF THE ELASTIC TAIL AND ITS
                    ///< ABSOLUTE NORMALISATION IS NOT VALIDATED AGAINST ANY
                    ///< EXTERNAL EXACT TAIL.  POLRAD sec. 2.1.3 B asserts the
                    ///< s- and p-peaks are SUPPRESSED for a tail (the s-peak
                    ///< sits at t ~ (1-y)Q^2, where 6Li's charge form factor
                    ///< is long dead, while the t-peak reaches down to
                    ///< t_min = M_A^2 x_A^2/(1-x_A) ~= (x M_N)^2, inside it),
                    ///< but that is ONE SENTENCE about a `approx` code path
                    ///< that is compiled OUT of the shipped POLRAD build
                    ///< (polrad20.cra:12,17), returns a NEGATIVE tail as
                    ///< printed, and carried an uncaught Y_- typo and a
                    ///< missing Z^2.  Measured against it here (T8(c), the
                    ///< leading-log s-peak): the s-peak is NOT negligible for
                    ///< a deuteron, whose form factor is still alive at
                    ///< t ~ 0.15 GeV^2, and the t-peak alone is LOW against
                    ///< HERMES's quoted "almost 50 % of the statistics in the
                    ///< lowest-x bin" radiative background.  Treat `rc_tail`
                    ///< as a LOWER BOUND on the dilution, band it, and never
                    ///< quote it as "the" radiative tail.
  PolradFull = 1,   ///< Eq. (18) + Appendix B + Eq. (A.4).  NOT IMPLEMENTED.
};
const char* rc_mode_name(RcMode m);          ///< "off", "tensor-band"

/// Which part of the rank-2 sector the band rescales.
enum class RcScope : int {
  TensorRate = 0,   ///< DEFAULT: the b1..b4 (T_LL) sector only -- what POLRAD
                    ///< and Gakh-Shekhovtsova actually compute.  It INCLUDES
                    ///< the O(gamma^2) `tensor_gamma` re-projection of that
                    ///< same sector onto cos phi' / cos 2phi', because
                    ///< excluding it would make the band depend on a purely
                    ///< kinematic option.
  TensorAll  = 1,   ///< also the Delta (gluon-transversity) cos 2phi term.
                    ///< NO LITERATURE SUPPORT: nobody has computed RC for a
                    ///< phi-dependent tensor observable, and Delta is not in
                    ///< POLRAD's b1..b4 basis.  For PRICING the omission
                    ///< only, never for correcting it.
};

struct RcOptions {
  /// NOTE: there is NO `mode` here.  `PipelineConfig::rc` is the single
  /// source of truth and `RcModel` takes the mode as a constructor argument.
  RcScope scope = RcScope::TensorRate;

  // --- the band ---------------------------------------------------------
  double delta_high_x = RC_DELTA_HIGH_X;
  double x_high       = RC_X_HIGH;
  double delta_low_x  = RC_DELTA_LOW_X;   ///< 0.30 default; 0.19 = "as good as HERMES"
  double x_low        = RC_X_LOW;

  // --- the radiative tails ----------------------------------------------
  bool with_tail    = true;   ///< false = band only (a diagnostic run)
  bool with_qe_tail = true;   ///< the UNPOLARISED quasi-elastic tail (Eq. 44).
                              ///< Its A_zz is ~0, so it DILUTES A_zz just as
                              ///< the unpolarised elastic tail does; HERMES
                              ///< subtracted both.  Turning it off prices the
                              ///< elastic tail alone and MUST be labelled so.
  RcTailModel tail_model = RcTailModel::TPeak;
  std::shared_ptr<const Spin1ElasticFF> ff;   ///< null => HoSpin1FF::for_ion
  /// The +-100 % quadrupole band, applied by `RcModel` when it builds the
  /// DEFAULT form factor.  Lives here (not only on `HoSpin1FFOptions`) so that
  /// `make_config(rc_fq_scale=...)` and `--rc-fq-scale` have something to set.
  /// IGNORED, with a printed line, when `ff` is user-supplied.
  /// QUADRATIC in the tail: run it, do not rescale (T12).
  double fq_scale = 1.0;
  /// Flat multiplier on F_m -- the eta F_m^2 tensor sector.  Same rules.
  double tail_tensor_scale = 1.0;
  /// Flat multiplier on the whole quasi-elastic tail, ON TOP of the
  /// `qe_kf_gev` Pauli suppression below -- the band knob, not the physics.
  /// Run 0.0 / 0.5 / 1.0.
  double qe_suppression = 1.0;
  /// POLRAD Eq. (44)'s S_E/S_M, as `ffquas` codes them: the de Forest-Walecka
  /// Fermi-gas factor S(q) = (3/4)(q/k_F) - (q/k_F)^3/16 below q = 2 k_F.
  /// DEFAULT ON, at 6Li's measured k_F (`RC_QE_KF_GEV`, Moniz et al.).  This
  /// is not a small choice: after the per-nucleon fix the QRT is the DOMINANT
  /// piece of `rc_tail` at every x, and the t-peak reaches down to
  /// t_min ~ (x M_N)^2, so at x <~ 0.1 most of its integral sits below 2 k_F.
  /// Set to 0 for the unsuppressed S = 1 edge (which v0 shipped as its
  /// default and which overstates the QRT there).
  double qe_kf_gev = RC_QE_KF_GEV;
  /// eta_A quadrature of Eq. (38): Gauss-Legendre nodes in ln(eta_A).
  int    n_eta       = 128;
  double m_lepton    = M_ELECTRON;     ///< constants.hpp -- NOT a second literal
  /// Ceiling on the returned tail ratio, mirroring `GlauberFsiOptions::w_max`:
  /// a Monte-Carlo weight must be bounded, and the clipped fraction is
  /// reported rather than hidden -- globally AND per y-band, because
  /// Y_+ = [1 + (1-y)^2]/(1-y) ~ 1/(1-y) makes the y -> 1 edge the only place
  /// it bites.
  double tail_max    = 10.0;
  /// Ceiling on |tau| the BAND sees, the same discipline `tail_max` applies to
  /// the tail.  On the TAGGED channels tau_tag = 1 - nbar(k,c)/n_M(k,c) is
  /// UNBOUNDED -- it diverges wherever the event's own n_M is near a node of
  /// the M-dependent spectator density, which the 6Li M = 0 density has -- and
  /// without this the published band edges go NEGATIVE (measured: rc_tensor_hi
  /// down to -1.79 on 6Li tagged-alpha, -8.35 on 7Li, 0.1-0.4 % of events).
  /// n_M -> 0 is exactly where the fractional-rescale ansatz breaks down: the
  /// tensor part of the density cancels the unpolarised part there, so "delta
  /// times the tensor fraction" stops being a small variation.  The clamp is a
  /// CHOICE, the clipped events are counted (`Event::rc_clipped`,
  /// `meta["rc_clipped_band_event_fraction"]`), and the right long-run answer
  /// is probably a tagged band reformulated on the M-averaged density.
  double band_tau_max = RC_BAND_TAU_MAX;
};

/// One event's RC weight triple, plus the two "this event hit a ceiling" flags
/// (`Event::rc_clipped` carries them into the record and the npz meta).
struct RcWeights {
  double lo   = 1.0;   ///< "rc_tensor_lo"  = 1 - delta(x) * tau
  double hi   = 1.0;   ///< "rc_tensor_hi"  = 1 + delta(x) * tau
  double tail = 1.0;   ///< "rc_tail"       = 1 + sigma_tail / sigma_Born
  bool band_clipped = false;   ///< |tau| hit `RcOptions::band_tau_max`
  bool tail_clipped = false;   ///< the tail ratio hit `RcOptions::tail_max`
};

/// `Event::rc_clipped` bits.
inline constexpr unsigned kRcClipTail = 1u;
inline constexpr unsigned kRcClipBand = 2u;

inline constexpr std::size_t kRcWeightCount = 3;
/// "rc_tensor_lo" | "rc_tensor_hi" | "rc_tail" for i = 0, 1, 2.  THROWS for
/// any other i.
const char* rc_weight_name(std::size_t i);
/// The HepMC3 name of entry `i` in SLOT `slot` (0 = the event's own category,
/// 1 + k = spin category k, matching the existing "spin_weight_k"):
/// slot 0 -> "rc_tail", slot 3 -> "rc_tail_3".
std::string rc_weight_name(std::size_t i, std::size_t slot);

/// The tensor-sector RC weight family.  IMMUTABLE after construction and safe
/// to share between threads, exactly like `GlauberFsiWeight`: the constructor
/// does all the work (the eta_A quadratures over the tail-table nodes).
///
/// BODIES LAND WITH THE MODEL (design_C_tensor_rc.md sec. 6, agent 1 steps
/// 4-8); this header is the API those bodies fill.
class RcModel {
 public:
  /// Inclusive / coherent form.  `dis` is the sampler the events were drawn
  /// from; the tail tables are built over ITS accepted (x, Q2) grid and read
  /// by interpolation at the event's own (x, Q2) (design sec. 1.4.5), with
  /// extra rows at y in {0.9, 0.95, 0.97, 0.98, 0.985} so the 1/(1-y) growth
  /// of Y_+ is resolved rather than extrapolated.
  ///
  /// THROWS when any category of `plan` has `lam_e * pe != 0`.  The whole
  /// A_zz programme -- and every tau identity in this header -- assumes an
  /// UNPOLARISED beam (design sec. 1.1); with a polarised one the vector
  /// sector would be priced as if it were rank 2.
  ///
  /// The constructor resolves, ONCE, a `const StateTables*` for every
  /// (category k, projection m) exactly the way `CategoryPlan::states` does,
  /// and stores the pointers.  It must NOT call
  /// `InclusiveSampler::state_tables()` per event: that takes `cache_mutex_`
  /// on EVERY call and would serialise `Pipeline::for_each(sink, 4)`.  Nor may
  /// it call `InclusiveKernel::tables()` in the event loop: `target_mass` is
  /// on by default and that is a 96-point g2^WW quadrature per call.
  RcModel(RcMode mode, RcOptions opt,
          std::shared_ptr<const InclusiveSampler> dis,
          RunPlan plan, Channel channel, Ion ion);
  /// Tagged form: adds the cluster model, whose n_M(k, c) carries the whole
  /// tensor structure when `StruckClusterOptions::inclusive_b1` is false (the
  /// default) -- the struck cluster's DIS kernel then has NO b1 at all and the
  /// inclusive recipe would silently return tau = 0 and price nothing (design
  /// sec. 1.5.1).  On this path the constructor ALSO rebuilds, per
  /// struck-cluster projection m_S, the PURE `SpinCategory(j = s_channel,
  /// e_{m_S}, lam_e, pe, theta_s = 0, phi_s = 0)` key that
  /// `InclusiveKinematicsSource::plan_for` uses, because that is the key the
  /// struck-cluster `StateTables` are cached under.
  RcModel(RcMode mode, RcOptions opt,
          std::shared_ptr<const InclusiveSampler> dis,
          RunPlan plan, Channel channel, Ion ion,
          std::shared_ptr<const TaggedModel> tagged);

  const RcOptions& options() const { return opt_; }
  RcMode mode() const { return mode_; }
  /// The form factor actually in use (the user's `RcOptions::ff` or the
  /// model-built default), so the run banner and Python can print its
  /// provenance -- `RcOptions::ff` is null on the default path and the built
  /// object would otherwise be unreachable.
  const Spin1ElasticFF& ff() const;
  std::string ff_provenance() const;   ///< == ff().provenance()

  /// The per-`Channel` rule, TOTAL over the enum (`event.hpp` has 9 values,
  /// `PipelineChannel` only 5; see the header block and design sec. 1.5):
  ///   Inclusive   -> band AND tail
  ///   every Tagged* (incl. TaggedLi6D/TaggedLi7T/TaggedDeuteronN/TaggedHe3P,
  ///                  which no PipelineChannel builds today) -> band, tail == 1
  ///   CoherentLi6 -> neither; both weights exactly 1.0
  bool applies() const { return applies_; }
  /// False whenever the tail is identically 1: any tagged channel (the elastic
  /// recoil sits at x_L = 1, inside the 10-sigma beam envelope, and is vetoed
  /// by the tag), `with_tail = false`, or `!applies()`.  NOTE it is TRUE at
  /// every theta_S -- the tail is emitted at any axis through
  /// P_zz^eff = 3 Q_NN P_2(cos theta_S) (POLRAD Eq. (43); design sec. 1.4.7).
  /// The run PRINTS the reason whenever this is false.
  bool tail_applies() const { return tail_applies_; }
  const std::string& exclusion_reason() const { return exclusion_reason_; }

  // --- the pieces, exposed so a run and a test can see them --------------
  double delta(double x) const;                        ///< rc_delta at the run's knobs

  /// The rank-2 fraction tau of the header block = W_tensor/W, for the event's
  /// OWN spin state.  READS EXACTLY THESE FIELDS, and nothing else, so a test
  /// can build a record by hand (T1, T2, T14):
  ///   ev.channel
  ///   ev.kin.cell, ev.kin.phi, ev.kin.x, ev.kin.q2,
  ///   ev.kin.k, ev.kin.cos_theta_k            (tagged only)
  ///   ev.spin.j, ev.spin.m_ion, ev.spin.m_struck,
  ///   ev.spin.lam_e, ev.spin.pe, ev.spin.theta_s, ev.spin.phi_s
  ///
  /// UNCLAMPED, on purpose: this is the DIAGNOSTIC that T1/T2 gate the closed
  /// form against, and `RcOptions::band_tau_max` is applied where the published
  /// weight is made (`weights()`), not here.  On a tagged channel it is
  /// unbounded -- see `RcOptions::band_tau_max`.
  double tensor_fraction(const Event& ev) const;
  /// ... for spin category `k` of the run plan (weighted mode).  This is the
  /// category's POPULATION MIXTURE, sum_m p_m W_m, mirroring
  /// `InclusiveSampler::weights_for` -- NOT a pure state.  It equals
  /// `tensor_fraction(ev)` only when category k's population vector is pure.
  /// Also UNCLAMPED.
  double tensor_fraction(const Event& ev, std::size_t k) const;

  /// sigma_tail / sigma_Born at accepted cell `c` and tensor degree `q_n`
  /// (q_n = P_zz^eff = 3 Q_NN P_2(cos theta_S); POLRAD's Q_N, which Eq. (37)
  /// carries as Q_N/6 -- dropping that 1/6 inflates the tensor tail sixfold
  /// and no unpolarised test sees it):
  ///
  ///   [ sigma^el_U + (q_n/6) sigma^el_T + kappa_qe sigma^q_U ]
  ///   ---------------------------------------------------------
  ///                 x * s * dsigma_unpol(x, Q2, s)
  ///
  /// `sigma_Born` here is the UNPOLARISED Born -- the library's own
  /// `InclusiveKernel::dsigma_unpol` and never a second definition of it
  /// (T7).  Design sec. 1.4.5's `(1 + w_avg)` denominator is NOT in this
  /// ratio and IS in `weights()`, which divides by the event's own density.
  /// Two reasons: this signature has no spin state to take `w_avg` from, and
  /// keeping it out makes the ratio EXACTLY LINEAR in `q_n`, which is what
  /// lets T14 assert POLRAD Eq. (43) (`sigma_q,perp = -1/2 sigma_q,par`) at
  /// 1e-12 instead of to a tolerance the unpolarised piece would set.
  ///
  /// Clipped at `RcOptions::tail_max`; `clipped_cell_fraction()` reports how
  /// often, globally and per y-band.
  double tail_ratio(int cell, double q_n) const;
  /// The same at an ARBITRARY (x, Q2), by interpolation of the tables --
  /// public so a test is not hostage to which cells the sampler accepted.
  /// Throws outside the table's (x, y) support.  `clipped`, when non-null, is
  /// set to whether `RcOptions::tail_max` bit on THIS call -- the node-level
  /// `clipped_cell_fraction()` is a different quantity (nodes are not
  /// event-weighted) and the two must never be quoted for each other.
  double tail_ratio_at(double x, double q2, double q_n,
                       bool* clipped = nullptr) const;
  /// The three precomputed tables, for plotting and for the T8 gates.  They
  /// are PER NUCLEON and in GeV^-2, over the (`table_x`, `table_y`) node grid
  /// in row-major (n_x x n_y) order.
  const std::vector<double>& sigma_tail_u() const { return sigma_u_; }
  const std::vector<double>& sigma_tail_t() const { return sigma_t_; }
  const std::vector<double>& sigma_tail_qe() const { return sigma_qe_; }
  /// The node grid the three tables live on: x, and y = Q^2/(x s).  The
  /// SECOND axis is y and not Q^2 because the tail carries
  /// Y_+ = [1 + (1-y)^2]/(1-y), which diverges as 1/(1-y): interpolating in y
  /// puts the design's refinement rows (y in {0.9, 0.95, 0.97, 0.98, 0.985})
  /// on the axis that actually needs them, where in (ln x, ln Q^2) they would
  /// be diagonal lines and not rows at all.
  const std::vector<double>& table_x() const { return node_x_; }
  const std::vector<double>& table_y() const { return node_y_; }

  /// The three tail pieces at an ARBITRARY (x, Q2), PER NUCLEON and in
  /// GeV^-2, straight from the quadrature -- no table, no interpolation.
  /// This is what the tables are built from and what T8's gates compare
  /// against, so that the per-nucleon reduction and the nuclear map are
  /// gated separately from the integrands (which `polrad_sigma_el_*` gate).
  struct TailTriple { double u = 0.0, t = 0.0, qe = 0.0; };
  TailTriple tail_sigma_at(double x, double q2) const;
  /// Fraction of table nodes that hit `RcOptions::tail_max`, globally and per
  /// y-band (y < 0.5, 0.5-0.9, > 0.9).  LOG ALL FOUR: the tail grows like
  /// Y_+ ~ 1/(1-y), so a clipped edge hides inside a small global number.
  double clipped_cell_fraction() const { return clipped_; }
  std::array<double, 3> clipped_fraction_by_y() const { return clipped_by_y_; }

  // --- what the pipeline calls ------------------------------------------
  RcWeights weights(const Event& ev) const;                    ///< slot 0
  RcWeights weights(const Event& ev, std::size_t k) const;     ///< slot 1+k
  /// Fill `ev.rc_weights` with the whole row-major (3 x n_slot) block,
  /// n_slot = 1 + ev.spin_weights.size().  A no-op when `mode == Off`.
  ///
  /// TAKES NO `Rng&`.  That signature is the STRUCTURAL guarantee that
  /// `--rc tensor-band` cannot move the random stream; T6 is the runtime one.
  void fill(Event& ev) const;

 private:
  // POLRAD Eq. (38): three one-dimensional eta_A quadratures at one (x, Q2)
  // -- sigma_u^A, sigma_q^A (the Q_N/6 partner) and the Eq. (44) QRT.
  void build_tail_tables();
  /// The per-nucleon Born d^2 sigma/(dx dy) [pb] the tail is a fraction of:
  /// `x * s * dsigma_unpol(x, Q2, s)`, i.e. the library's OWN unpolarised
  /// Born and never a second definition of it (T7).
  double born_pb_at(double x, double q2) const;
  /// The `StateTables` of one PURE spin state, resolved in the constructor;
  /// throws when the state is not one the run plan carries.
  const InclusiveSampler::StateTables* state_of(const SpinLabels& sp) const;
  /// The rank-2 alignment degree Q_N = P_zz^eff = 3 Q_NN P_2(cos theta_S) of
  /// one pure state, through the library's own `tensor_moments`.
  double q_n_of(double m, double theta_s) const;
  /// tau on a Tagged* channel: 1 - nbar(k, c)/n_M(k, c) (design sec. 1.5.1).
  double tagged_tau(const Event& ev, const std::vector<double>* pops) const;
  /// tau clipped to +-`RcOptions::band_tau_max`, which is what the published
  /// band edges are built from; sets `*clipped` when the clamp bit.
  double clamp_tau(double tau, bool* clipped) const;

  RcMode mode_ = RcMode::Off;
  RcOptions opt_;
  std::shared_ptr<const InclusiveSampler> dis_;
  std::shared_ptr<const TaggedModel> tagged_;
  std::shared_ptr<const Spin1ElasticFF> ff_;
  RunPlan plan_;
  Channel channel_ = Channel::Inclusive;
  Ion ion_;
  bool applies_ = false;
  bool tail_applies_ = false;
  std::string exclusion_reason_;
  // Resolved ONCE in the constructor: [category k][projection index m].
  std::vector<std::vector<const InclusiveSampler::StateTables*>> states_;
  // The tail tables, over the (ln x, ln y) nodes of `node_x_` x `node_y_`.
  std::vector<double> node_x_, node_y_;
  std::vector<double> sigma_u_, sigma_t_, sigma_qe_, born_pb_;
  std::vector<double> m_val_;   ///< m_values(plan J), hoisted out of the loop
  double s_ = 0.0;              ///< the sampler's per-nucleon s
  double clipped_ = 0.0;
  std::array<double, 3> clipped_by_y_{{0.0, 0.0, 0.0}};
};

/// `PipelineConfig::rc`.  Alias of `RcMode` so the CLI, the config struct and
/// the model cannot disagree (the FSI precedent keeps two enums; this one
/// deliberately does not).
using PipelineRc = RcMode;
const char* pipeline_rc_name(PipelineRc r);

}  // namespace lipolgen

#endif  // LIPOLGEN_RC_HPP
