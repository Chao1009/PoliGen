// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LIPOLGEN_B1_NUCLEAR_HPP
#define LIPOLGEN_B1_NUCLEAR_HPP

/// \file b1_nuclear.hpp
/// b1 of a two-cluster nucleus by convolution: Cosyn-Dong-Kumano-Sargsian
/// (CDKS), PRD 95 (2017) 074036, Eqs. (16), (17), (21), applied one level up
/// (alpha + d instead of p + n).  See
/// docs/open_items/run_2026-09-02/design_D_b1_li6.md for the derivation, the
/// symbol table and the FOUR-term truncation argument.
///
/// THE A = 2 VALIDATION GATE PASSES (2026-09-03), AND READ ITS CONDITIONS.
/// `G3a` (shape) and `G3b` (magnitude within a factor 2 of the digitized CDKS
/// Fig. 4 peak) both pass when read where design section 5.4's checklist ends
/// -- MSTW2008 LO, CDKS's own nucleon PDF (checklist item 4, `MstwSF` in
/// mstw_sf.hpp), at CDKS Eq. (21)'s delta-function, which is
/// `DeuteronConvolutionB1`'s default:
///
///   G3a  zeros 0.0279 (falling) / 0.4952 (rising), peak at 0.7716, against
///        the digitized 0.0656 / 0.4572 / 0.7657 -- misses 0.038 / 0.038 /
///        0.0059 against tolerances 0.08 / 0.10 / 0.10.
///   G3b  max|x b1| over [0.10, 0.80] = 9.15096e-4 against 1.08521e-3,
///        ratio 0.843243 (a factor 1.19 low), inside [0.5, 2].
///   G3c  int b1 dx = 2.24896e-4 against 4.59200e-4 -- reported, NOT enforced.
///
/// FOUR CONDITIONS, none of them optional reading.  (1) THE LIFT IS ABOUT A
/// CONFIGURATION, NOT ABOUT A BUILD.  The gate passes for the MSTW2008 LO
/// nucleon input at CDKS Eq. (21)'s delta-function; the SHIPPED DEFAULT
/// unpolarised backend, `ToyF2`, gives 0.440 on the same clause -- OUTSIDE
/// the [0.5, 2] acceptance window.  So: quote numbers made with the
/// unpolarised-backend selector set to `mstw` (`PipelineConfig::b1_unpol`,
/// CLI `--b1-unpol mstw`, which reaches `Li6ConvolutionOptions::unpol`
/// through `default_inclusive_kernel`); the default toy backend is outside
/// the gate's acceptance window and ITS NUMBERS ARE NOT COVERED BY THE LIFT.
/// The two are not a rescaling of each other: on the shipped observable
/// `Li6ConvolutionB1::b1(x, 2.5)` mstw/toy is 1.848 / 1.276 / 0.817 at
/// x = 0.10 / 0.30 / 0.50 -- up to a factor 1.85 and not monotone.  Selecting
/// `mstw` needs the optional PYTHIA tier; a build without it REFUSES the
/// selector at configuration time and is never silently downgraded to the
/// toy, so a build that cannot reproduce the verdict row also cannot emit a
/// number that claims it.  (2) It is comfortable at Eq. (21) and marginal at
/// Eq. (17): MSTW at
/// kappa = 1 gives 0.520, inside by 4 % of its own value.  (3) With CD-Bonn as
/// well (`Options::wave = kCdBonn`, CDKS's own wave function) the ratio is
/// 1.000338 and the zeros 0.0641 / 0.4570 -- below the error of digitizing a
/// published figure, specific to CD-Bonn AND MSTW together, and NOT the
/// default.  (4) The gate is A = 2.  It validates this kernel on the DEUTERON
/// and says nothing about the alpha-d step, for which no measurement exists at
/// any A > 2 -- so the mandatory +-100 % band on every 6Li number stays.  It
/// comes from Q(6Li) vs Q_d, not from the gate.
///
/// Design section 5.4's "Escalation" clause -- "if AFTER the checklist the
/// peak ratio is still outside a factor of 2" -- is therefore not triggered,
/// and the 6Li publication ban it imposed is LIFTED.  Measurements:
/// docs/open_items/run_2026-09-03/phase_A_numbers.md (sections 0-5 the MSTW
/// rerun, section 8 CD-Bonn) and phase_A_cdbonn.md; the superseded failing
/// verdict is docs/open_items/run_2026-09-02/phase_D_gate.md.
///
/// G3a's margins, stated rather than left for the next reader to find: the
/// low-x zero is resolved to a few per cent at best on any of these grids and
/// must not be quoted to more than two significant figures; its counting
/// window is (0, 1.0] and the scan floor 0.001, both stated in the design
/// since 2026-09-03 (the old [0.02, 1.0] window was narrower than the +-0.08
/// tolerance it bracketed, which made that tolerance dead code, and with
/// CT18NLO the zero fell below it and aborted the clause).  The second zero
/// and the peak position are the robust discriminators.
///
/// G3A IS A QUALIFIED PASS, and the qualification is the counting window's
/// UPPER edge.  (0, 1.0] is a scope choice that no position tolerance derives.
/// Over the digitized reference's OWN domain [0.010, 1.590] the reference has
/// TWO sign changes and stays positive to its last point, while EVERY computed
/// configuration crosses zero a THIRD time -- ToyF2 1.220437, MSTW 1.217660,
/// CT18NLO 1.197722 at Eq. (21) with AV18 (T17 pins these).  x > 1 per nucleon
/// is kinematically allowed for a nucleus, so that region is physical and not
/// out of range, and at those crossings the digitized column is still 7.0-8.5 %
/// of its own peak -- not a digitization floor.  Widen the window to (0, 1.59]
/// on the same argument that widened the floor and G3a fails on every
/// configuration.  So: never quote "exactly two sign changes" without the
/// window.  The argument in full, and why the clause is recorded rather than
/// changed: docs/OPEN_ITEMS_SOLUTIONS.md section 10, "G3a's stated limitation
/// -- the counting ceiling", and design_D_b1_li6.md's amendment.
///
/// -------------------------------------------------------------------------
/// THE MODEL (design 1.7), per nucleon, 6Li:
///
///   b1^6Li(x,Q2) = (2/6) int (dz/z) { [ f_S(z) + w_CG f_D(z) ] b1_d(x/z,Q2)
///                                     + w_ad delta_T f_ad(z) F1_d(x/z,Q2) }
///                + (4/6) w_ad int (dz/z) delta_T f_alpha(z) F1_alpha(x/z,Q2)
///
/// with w_CG = 1/10 EXACTLY (the Clebsch-Gordan tensor dilution of the
/// deuteron inside the L = 2 alpha-d component, design 1.4) and w_ad = 1
/// nominal.  The four terms are
///   (1)  embedded deuteron, S wave      -- `b1_embedded_s`
///   (2d) alpha-d D wave, struck d       -- `b1_alpha_d_dwave_d`
///   (2a) alpha-d D wave, struck alpha   -- `b1_alpha_d_dwave_alpha`
///   (3)  CG depolarization, D wave      -- `b1_cg_dwave`
///
/// WHY FOUR AND NOT THREE.  CDKS Eq. (10) sums the spectral function over
/// CONSTITUENTS; one level up that sum runs over {d, alpha}.  The alpha is
/// J = 0 so b1^alpha == 0, but its light-cone density carries the SAME
/// (3cos^2 theta - 1) orbital alignment as the deuteron's -- (3c^2 - 1) is
/// even in k, and the alpha carries -k -- so it contributes to b1 through
/// term (2a).  Its size is fixed by counting and by the 1/M^2 of the
/// P2-weighted density: (2a)/(2d) ~ [(4/6)/(2/6)] (M_d/M_alpha)^2 = 0.506,
/// which the quadrature reproduces (T15).  That is a CONSISTENCY check, not a
/// validation: no published two-cluster b1 convolution exists to check the
/// per-nucleon bookkeeping of the alpha piece against (design 9 Q11).
///
/// WHY NOT SIX.  Two pieces are dropped, both far inside the mandatory 100 %
/// band and both MEASURED rather than asserted:
///   * the anisotropic remainder of the D-wave b1_d weight,
///     -(9/16pi)[(72/35) P4 - (4/21) P2] |phi_2|^2, reported by
///     `b1_cg_dwave_p2_remainder` / `_p4_remainder` and pinned by T14.  It is
///     x-dependent because the P2- and P4-weighted densities integrate to
///     zero, so it enters only through the curvature of b1_d(x/z);
///   * the S-D interference in the b1_d sector: term (2d)'s SD structure with
///     coefficient 1/3 of it and b1_d in place of F1_d, i.e. suppressed by
///     b1_d/F1_d ~ 1e-3, about 3e-4 of term (2d).
///
/// THE 100 % BAND IS MANDATORY.  Q(6Li) = -0.0818(17) fm^2 against
/// Q_d = +0.2859(3) fm^2: the alpha-d relative D wave enters the closest
/// measured observable with the OPPOSITE sign to the deuteron's own D state
/// and nearly cancels it.  (Q(6Li): the repository's single copy is
/// `LI6_QUADRUPOLE_FM2` in rc.hpp, TUNL's A = 6 evaluation, 1998CE04.
/// Pyykko's compilation, Mol. Phys. 106 (2008) 1965, quotes -0.0806(6) fm^2
/// from the molecular-beam measurement of Cederberg et al., Phys. Rev. A 57
/// (1998) 2539 -- a 1.5 % difference that changes nothing here.  Q_d: Bishop
/// and Cheung, Phys. Rev. A 20 (1979) 381.)
/// Never quote a single row: every published number is {0, 1, 2} x b1
/// (`banded`).
///
/// NO RNG anywhere, and every result is a pure function of the constructor
/// arguments: determinism is structural.  The only mutable state is two
/// write-once caches -- `ClusterPartialWave::spline` (built by every factory,
/// so the lazy path is a fallback) and the `finite_q_delta` per-(x, Q2)
/// densities -- neither of which changes any value.  They make the objects
/// SHAREABLE but not concurrently callable from several threads; construct one
/// per thread if you need that.  All quadrature is `numerics.hpp`
/// (`trapezoid`, `linspace`) except the inner k integral, which is Simpson on
/// the same grid -- see `simpson` in b1_nuclear.cpp for the measurement that
/// forced it.

#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "lipolgen/cluster.hpp"   // VmcRadial, read_anl_*, read_fdeut_k, data_path
#include "lipolgen/sf.hpp"        // TensorSF, UnpolSF, RFunc, SFFunc3

namespace lipolgen {

/// Natural cubic spline, defined in b1_nuclear.cpp; `ClusterPartialWave`
/// caches one so the O(n_y n_k) quadrature does not rebuild it per call.
class CubicSpline;

/// F1 of CDKS Eq. (22): (1 + gamma^2) F2 / (2 x (1 + R)).  NOT
/// `UnpolSF::f1_from_f2` and NOT `NuclearF2::f1a`, which are the massless
/// form -- the extra factor is exactly `1 + gamma_squared(x, q2)`
/// (asymmetries.hpp), and gamma^2 = 0.90 at x = 0.8, Q2 = 2.5.  The
/// convolution kernel MUST be fed this one; `f1a` appears nowhere in
/// b1_nuclear.cpp.  A null `r_func` means `r_sigma_lt` (`resolve_r`).
double f1_cdks(double f2, double x, double q2, const RFunc& r_func = nullptr);

/// The RAW digitized CDKS theory-1 column as a `TensorSF`, i.e. b1 = (x b1)/x.
///
/// CDKS Eq. (10) carries an explicit 1/A, the text under their Eq. (16) says
/// "b1 is defined by the one per nucleon", their f(y) is normalised to ONE
/// nucleon, and their Fig. 6 overlays HERMES's per-nucleon b1 on the same
/// axis -- so the digitized curve is ALREADY per nucleon.  Until 2026-09-03
/// `b1_convolution()` (sf.cpp) halved it a second time, which made `CdksB1` a
/// factor 2 low, and this accessor existed to route around that; the constant
/// has since been split per camp (`B1_CDKS_TABLE_TO_PER_NUCLEON` = 1,
/// constants.hpp) and `CdksB1` now agrees with this accessor exactly.  It
/// stays because `Li6ConvolutionOptions::deuteron_b1 == nullptr` uses it and
/// because it names the normalisation at the point of use.
std::shared_ptr<const TensorSF> cdks_b1_raw_per_nucleon();

/// One partial wave of a two-cluster relative wave function in the CDKS
/// convention: phi_L(k) WITH the i^L factor.
///
/// THE i^L PHASE (CDKS below their Eq. 20, and again below their Eq. 40).
/// CDKS define phi_L(p) = 4 pi i^L int dr r^2 j_L(pr) u_L(r) and state
/// explicitly that the i^L makes phi_2(p) < 0 for the deuteron, "although a
/// different convention (phi_2 -> -phi_2), namely without the i^2 factor, is
/// sometimes used".  Below their Eq. (40) they fix phi_0(k) = U(k) and
/// (-i)^2 phi_2(k) = W(k), i.e. phi_2 = -W with U, W >= 0 at low k.  This
/// sign is LOAD-BEARING: Eq. (21)'s SD term is -(3/(4 sqrt2 pi)) phi_0 phi_2,
/// so getting i^L wrong flips terms (2d) and (2a) and moves the whole 6Li
/// answer by ~100 %.
///
/// THE ALPHA-D NODE TABLE (design 2.2), and the sign is k-DEPENDENT -- do not
/// describe it as uniform.  S node at 0.678 fm^-1 = 0.134 GeV, D node at
/// 2.25 fm^-1 = 0.444 GeV (from `li6.ad`, confirmed by the minima of the
/// momentum file's rho_0 / rho_2).  With the unobservable global phase fixed
/// to psi_0(k -> 0) > 0:
///
///   k [fm^-1]     region              sign(psi_2/psi_0)  phi_0 phi_2   weight
///   < 0.678       below the S node          -1            > 0          66.4 %
///   0.678 - 2.25  between the nodes         +1            < 0          33.4 %
///   > 2.25        above the D node          -1            > 0           0.2 %
///
/// ("weight" is the fraction of int k^2 (rho_0 + rho_2) dk.)  So two thirds of
/// the density has phi_0 phi_2 > 0, OPPOSITE to the deuteron, and one third
/// does not; delta_T f_ad(z) integrates ACROSS the node and the net sign of
/// terms (2d)/(2a) is a computed output, not an assumption.  The microscopic
/// statement of the 6Li quadrupole puzzle: `alpha_d_quadrupole_fm2` below is
/// the sign gate.
///
/// CONTRACT of `from_vmc`: it applies ONLY phi_L = i^L psi_L (phi_0 = psi_0,
/// phi_2 = -psi_2).  It does NOT touch the global phase, because the global
/// phase is a property of the PAIR and is already fixed on the pair by
/// `li6_vmc_waves()` (tagged.cpp:107-110, psi_0(k -> 0) > 0).  A `from_vmc`
/// that re-fixed the sign per wave (e.g. forcing psi.front() > 0 on the D
/// wave) would flip the S-D relative sign and hand the sign gate the wrong
/// answer.  `from_vmc_pair` is the recommended entry point because it cannot
/// be called wrongly.
struct ClusterPartialWave {
  int l = 0;
  std::vector<double> k;    ///< GeV, strictly increasing
  std::vector<double> phi;  ///< CDKS-convention amplitude, signed
  std::string provenance;

  /// NATURAL CUBIC SPLINE in k, and 0 OUTSIDE the table -- matching
  /// `VmcRadial::operator()` (cluster.cpp:231), NOT the end-clamped value
  /// `np_interp` returns.  The spline is a CORRECTNESS item, not a style one:
  /// `fdeut.av18` is tabulated on a coarse 0.1 fm^-1 grid and linear
  /// interpolation of u, w biases int k^2 (u^2 + w^2) dk by +2.2 % (1.0218
  /// linear against 1.00035 spline on a fine grid).  The baryon-number
  /// renormalisation then HIDES that error and no grid-refinement test can
  /// see it, because the table spacing is the problem (design 3, 5.2 G1c).
  double operator()(double kk) const;
  /// The same on a SORTED grid, in one monotone pass over the spline's own
  /// nodes -- this is the inner loop of `LightConeDensities` and doing it
  /// point by point costs a binary search each time.
  void eval_sorted(const std::vector<double>& kk, std::vector<double>* out) const;
  /// Rebuild the cached spline.  Call after mutating `k` or `phi`; every
  /// factory below already does.
  void rebuild();

  /// int k^2 phi^2 dk over the table (trapezoid on the table's own nodes).
  double norm2() const;
  /// The same table times `factor`.
  ClusterPartialWave scaled(double factor) const;

  static ClusterPartialWave from_vmc(const VmcRadial& v, int l);
  /// Recommended: takes the pair whose global phase is already fixed.
  static std::pair<ClusterPartialWave, ClusterPartialWave>
      from_vmc_pair(const VmcRadial& s, const VmcRadial& d);
  /// phi_0 = u, phi_2 = -w (the CDKS i^L convention) from a (U, W) column.
  static ClusterPartialWave from_uw(std::vector<double> k_gev,
                                    std::vector<double> uw, int l);

  /// The cached spline; internal, and rebuilt by `rebuild()`.  Held by
  /// shared_ptr so that copying a wave is O(1) and thread-safe to read.
  mutable std::shared_ptr<const CubicSpline> spline;
};

/// The 6Li alpha-d relative partial waves in the CDKS convention, normalised
/// so that int k^2 (|phi_0|^2 + |phi_2|^2) dk = the file's own printed
/// S + D norm (0.8194650 by trapezoid, `VMC_N_ALPHA_D_LI6` = 0.819481 as the
/// file prints it -- 2.0e-5 apart, which is the file's own quadrature spread
/// and is what T3 measures).
///
/// Magnitudes from `momenta/li6_ad1.momentum`, sign structure from
/// `li6_alpha_d/li6.ad`, exactly as `li6_vmc_waves()` builds them -- this
/// calls `li6_alpha_channel(..., ClusterWaveSource::VmcAV18)` rather than
/// re-implementing the pattern.  The unit factor is
/// 1/sqrt(2 pi^2 (hbar c)^3): the ANL momentum file tabulates
/// rho_L = A_L^2/(4 pi) [fm^3] and prints 4 pi int rho K^2 dK/(2 pi)^3, so
/// int K^2 rho dK = 2 pi^2 x (printed norm) in fm units.
std::pair<ClusterPartialWave, ClusterPartialWave> li6_alpha_d_partial_waves();

/// Kinematics of one struck constituent in the CDKS Eq. (12) reduction
/// p^0 = M - eps - k^2/(2 m_recoil).  For 6Li there are TWO of these:
/// {M_d, M_alpha} for the struck deuteron (terms 1, 2d, 3) and
/// {M_alpha, M_d} for the struck alpha (term 2a).
///
/// THE FAIR-SHARE VARIABLE, AND HOW IT DIFFERS FROM THE EXACT PER-NUCLEON x.
/// z here is (E - k_z kappa)/m_struck, i.e. the constituent's light-cone
/// fraction normalised to its own MASS -- CDKS Eq. (18)'s pattern one level
/// up (design 1.6), which is what makes the 6Li convolution literally their
/// Eq. (16) with N -> d (resp. alpha) and puts the 1/3 in exactly one place,
/// the per-nucleon counting factor.  The library's x is Q^2/(2 M_N nu)
/// (design 1.1), so the EXACT map from x to the struck constituent's own
/// per-nucleon Bjorken variable would divide by A_i M_N instead --
/// z_d = p_d.q/(2 M_N nu), z_alpha = p_alpha.q/(4 M_N nu) -- which is the
/// same thing up to 2 M_N/M_d = 1.0005 and 4 M_N/M_alpha = 1.0069, the
/// nuclear binding.  Measured on term (2a) with the denominator alone
/// changed: +0.1 / +0.7 / +2.1 / +1.6 % at x = 0.1 / 0.3 / 0.5 / 0.7, and
/// <= 0.2 % on (2d).  Inside the mandatory 100 % band, and the mass form is
/// what the design pinned, so it stands -- but it IS a choice.
struct ConvolutionKinematics {
  double m_struck = 0.0;    ///< M_d for the struck d, M_alpha for the struck a
  double m_recoil = 0.0;    ///< the OTHER cluster: p^0 = M - eps - k^2/(2 m_r)
  double separation = 0.0;  ///< epsilon [GeV]
  /// |q|/nu = sqrt(1 + gamma^2) of CDKS Eq. (21)'s delta-function.  1.0 is
  /// Eq. (17) as printed (the DEFAULT); `finite_q_delta` makes it
  /// kappa(x, Q2) and then the densities are rebuilt per x.
  double kappa = 1.0;

  /// Allowed |k| interval at fixed y where |cos theta*| <= 1 (design 1.3):
  ///   B = m_struck (1 - y) - eps,  r = sqrt(kappa^2 + 2 B / m_recoil),
  ///   k in [m_recoil |kappa - r|, m_recoil (kappa + r)].
  /// Returns false when B < -m_recoil kappa^2 / 2 (r imaginary).  BOTH
  /// endpoints matter: getting only the lower one right inflates int f dy by
  /// a factor 6 (design 5.1 G0g).
  bool k_range(double y, double* lo, double* hi) const;
  /// [m_struck (1 - y) - eps - k^2/(2 m_recoil)] / (k kappa).
  double cos_star(double k, double y) const;
  /// 1 - eps/m_struck + (m_recoil/m_struck) kappa^2/2, kappa-DEPENDENT.
  ///
  /// DEVIATION from design 3, which writes `1 - eps/m_struck + kappa^2/2`.
  /// That is the A = 2 special case: the mass in the discriminant is the
  /// RECOIL mass (design 1.6's own r = sqrt(kappa^2 + 2B/m_recoil)), and for
  /// the deuteron m_recoil == m_struck so the two coincide (1.49763 at
  /// kappa = 1, design 5.2's G1e).  In 6Li they do not: the struck deuteron
  /// (m_recoil/m_struck = M_alpha/M_d = 1.987) has z_max = 1.99286 where the
  /// design's form would cut it at 1.49921, and the struck alpha
  /// (M_d/M_alpha = 0.503) has 1.25120 where the design's form would run on
  /// to 1.49960 -- too short in one direction and too long in the other.
  double y_max() const;
};

/// The light-cone densities of CDKS Eqs. (17) and (21) for ONE struck
/// constituent, tabulated on a y grid and interpolated.  Deterministic; no
/// RNG; safe to share.  6Li holds two of these.
///
/// The delta-function collapses the angular integral analytically (design
/// 1.3): with c* fixed by the delta and a Jacobian m_struck/(k kappa),
///   f_S(y)       = (1/(2 kappa)) y m_s int k dk |phi_0|^2 ,
///   f_D(y)       = (1/(2 kappa)) y m_s int k dk |phi_2|^2 ,
///   delta_T f(y) = (3/(2 kappa)) y m_s int k dk [U W/sqrt2 + W^2/4](3c*^2 - 1)
/// over the EXACT interval `k_range(y)` (never a mask over a fixed grid --
/// that costs 2 % on int f dy).
class LightConeDensities {
 public:
  struct Options {
    /// y-grid layout, spelled out so two implementers get the same bits.
    /// Three segments, trapezoid on each, breakpoints included in both.
    ///
    /// DEVIATION from design 3, which gives 600 / 2400 / 800.  Those are the
    /// numbers the A = 2 prototype was quoted at, and they are too coarse in
    /// the OUTER segments for 6Li: the struck deuteron's z-width scales with
    /// M_alpha/M_d = 1.987, so a bigger share of its CANCELLING delta_T f
    /// lives outside (0.85, 1.15) than the deuteron's does, and term (2d) came
    /// out 2.3 % high at x = 0.10 and 0.8 % high at x = 0.30 on the design's
    /// grid (measured 2026-09-03 against an independent 2-D Gauss-Legendre
    /// evaluation of Eq. (21) with no y table at all, and against a x4
    /// refinement of this one, which agree).  2400 / 2400 / 3200 lands term
    /// (2d) within 1.2e-3 of the x4-refined value and (2a)/(2d) on the
    /// analytic 2 (M_d/M_alpha)^2 = 0.5064 rather than 0.489-0.497.  It costs
    /// ~0.14 s per density set, once, at construction.
    std::size_t n_y_low = 2400;   ///< [y_min, y_mid_lo]
    std::size_t n_y_mid = 2400;   ///< (y_mid_lo, y_mid_hi]  -- delta_T f's structure
    std::size_t n_y_high = 3200;  ///< (y_mid_hi, y_max]
    double y_mid_lo = 0.85;
    double y_mid_hi = 1.15;
    /// k points per y, between the EXACT endpoints.  MUST BE ODD: the inner k
    /// integral is composite Simpson (see `simpson` in b1_nuclear.cpp for the
    /// measurement that forced it) and an even count would silently fall back
    /// to the trapezoid, which is 22 % low at x = 0.05.  An even value is
    /// bumped to the next odd one by the constructor and `options()` then
    /// reports what actually ran -- a round number like 2000 or 4000 is what a
    /// user refining the grid types, and it must not change the quadrature.
    std::size_t n_k = 2001;
    /// Bottom of the y grid.  DEVIATION from design 3, which hard-codes
    /// 1e-4: it is an `Options` field here so that G0e can be tested on the
    /// FULL support.  int delta_T f dy = 0 exactly (design 5.1 G0e), but only
    /// over the whole support, which reaches NEGATIVE y -- E = m - eps -
    /// k^2/(2 m_r) turns negative above k = sqrt(2 m_s m_r).  On the default
    /// 1e-4 grid the residual sticks at 1.1e-4 x max|delta_T f| and does not
    /// improve with refinement (it is a truncation, not a discretization);
    /// at y_min = -3 it drops to 2.7e-6 x max, which is the real test of the
    /// identity.  The DEFAULT is unchanged from the design.
    double y_min = 1e-4;
    /// Renormalise by the baryon-number condition int f dy = `norm_target`
    /// (CDKS below Eq. 21).  `norm_target` = 1 for the deuteron,
    /// VMC_N_ALPHA_D_LI6 for the alpha-d cluster (the spectroscopic factor is
    /// REAL suppression, design 1.7 A1).  `renormalize = false` exposes the
    /// raw quadrature, which is what G1c / T3 / T4 test.
    bool renormalize = true;
    double norm_target = 1.0;
  };

  LightConeDensities(ClusterPartialWave phi0, ClusterPartialWave phi2,
                     ConvolutionKinematics kin, Options opt);
  /// The same with the default `Options`.  (Two overloads rather than
  /// `opt = {}`: a nested aggregate's default member initializers are not
  /// usable in a default argument of its own enclosing class.)
  LightConeDensities(ClusterPartialWave phi0, ClusterPartialWave phi2,
                     ConvolutionKinematics kin);
  /// Test-only: the no-smearing limit, f_S -> (1 - p_d) norm delta(z - 1),
  /// f_D -> p_d norm delta(z - 1), delta_T f -> 0.  Flagged internally so
  /// that `convolve` short-circuits to `dens_at(1) * g(x)` instead of
  /// quadrature (a delta cannot live on a trapezoid grid).  Used by T6, which
  /// therefore tests the WEIGHT ASSEMBLY only -- the kernel is T2/T5.
  static LightConeDensities delta_limit(double p_d, double norm);

  double f_s(double y) const;        ///< S-wave, Eq. (17) with |phi_0|^2
  double f_d(double y) const;        ///< D-wave, Eq. (17) with |phi_2|^2
  double f_unpol(double y) const;    ///< f_s + f_d
  double delta_t_f(double y) const;  ///< Eq. (21), SD + DD
  double delta_t_f_sd(double y) const;
  double delta_t_f_dd(double y) const;
  /// P2- and P4-weighted D-wave densities (same normalisation as `f_d`) of
  /// design 1.5 truncation item 1, used by T14 to MEASURE the angular-average
  /// error instead of asserting it.
  double f_d_p2(double y) const;
  double f_d_p4(double y) const;

  double norm() const;               ///< int f_unpol dy
  double mean_y() const;             ///< int y f dy / int f dy == <y^2>_phi/<y>_phi
  double p_d() const;                ///< int f_d dy / norm()   (y-WEIGHTED, see T3)
  /// Unweighted momentum-space ratio int k^2 |phi_2|^2 / int k^2 (|phi_0|^2 +
  /// |phi_2|^2): THIS is the quantity that equals the file's D/(S+D), and the
  /// one T3 pins to `VMC_P_D_LI6`.
  double p_d_momentum() const;
  const std::vector<double>& y_grid() const;
  const ConvolutionKinematics& kinematics() const;
  const Options& options() const;
  /// The scale the baryon-number renormalisation applied (1 when off).
  double renormalization() const;

  /// int (dy/y) dens(y) g(x/y), trapezoid on a uniform refinement of the
  /// stored support to `n` points from max(x / x_max_g, y_min) to `y_max`.
  ///
  /// THE SUPPORT OF `g` IS `g`'S BUSINESS, AND `x_max_g` IS HOW `g` SAYS IT.
  /// The lower limit is x / x_max_g, not x: a free-nucleon F1 has x_max_g = 1
  /// (the default) and the integral then starts at x, exactly as before; an
  /// injected b1_d TABLE is defined out to its own `x_max()` (kB1CdksQ2p5 runs
  /// to 1.59, and the deuteron's per-nucleon support genuinely exceeds 1, so
  /// the strength above x/y = 1 is real), and passing that x_max_g is what
  /// keeps it.  Measured on term (1) against an independent quadrature:
  /// cutting at x/y = 1 drops 0 % below x = 0.8, 1.6 % at x = 0.9 and 5.1 % at
  /// x = 0.95 -- the cell the CLI's `--x-max 0.955` sits in.  Do NOT pass an
  /// x_max_g larger than where `g` is actually defined: `RawCdksB1` returns 0
  /// above the table's x_max for exactly that reason, but a `g` that
  /// end-clamps instead would leak forever.
  double convolve(const std::function<double(double)>& dens_at,
                  const std::function<double(double)>& g, double x,
                  double x_max_g = 1.0, std::size_t n = 4001) const;

 private:
  LightConeDensities() = default;
  double lookup(const std::vector<double>& t, double y) const;

  ClusterPartialWave phi0_, phi2_;
  ConvolutionKinematics kin_;
  Options opt_;
  std::vector<double> y_, fs_, fd_, sd_, dd_, p2_, p4_;
  double raw_norm_ = 0.0, scale_ = 1.0, norm_ = 0.0, mean_y_ = 0.0;
  double p_d_ = 0.0, p_d_mom_ = 0.0;
  bool delta_ = false;
  double delta_s_ = 0.0, delta_d_ = 0.0;
};

/// Which deuteron b1 camp and which alpha-d inputs a `Li6ConvolutionB1` is
/// built from, plus the knobs the 100 % band is expressed with.
struct Li6ConvolutionOptions {
  /// b1 of the EMBEDDED deuteron, per nucleon.  Null => the RAW digitized
  /// CDKS theory-1 column (`cdks_b1_raw_per_nucleon()`), because CDKS
  /// Eqs. (10)/(16) and their Fig. 6 say the published curve is already per
  /// nucleon (design 9 Q1).  `MillerB1` is a DIFFERENT CAMP -- mixing them is
  /// legal but must be said out loud in any plot legend.
  std::shared_ptr<const TensorSF> deuteron_b1;
  /// The x above which `deuteron_b1` is ZERO -- its support, which `convolve`
  /// needs and a `TensorSF` cannot be asked for.  0 => automatic: the
  /// digitized table's own `x_max()` = 1.59 when `deuteron_b1` is null (the
  /// raw CDKS column), and 1.0 for any injected backend, which is the
  /// conservative reading for an analytic camp such as `MillerB1`.  Set it
  /// explicitly when injecting a table that runs past x = 1.
  double deuteron_b1_x_max = 0.0;
  /// F1 of the embedded deuteron, per nucleon.  Null =>
  ///   f1_cdks(NuclearF2(DEUTERON(), unpol, nullptr, r_func).f2a(x,q2)/2,
  ///           x, q2, r_func)
  /// -- CDKS Eq. (22), NOT `NuclearF2::f1a` (design 2.4).
  std::shared_ptr<const UnpolSF> unpol;
  /// F1 of the alpha, per nucleon, for term (2a).  Null => the isoscalar
  /// (F1p + F1n)/2 built the same way, i.e. numerically identical to the
  /// deuteron's because Z = N in both; the alpha's own EMC effect (~10 % at
  /// x ~ 0.6) is NOT modelled (design 1.7 A9, 9 Q10).
  SFFunc3 alpha_f1;
  /// R = sigma_L/sigma_T for this backend's own F1 slots.  Null =>
  /// `r_sigma_lt`, DELIBERATELY, and NOT the `r1998` the A = 2 gate defaults
  /// to (decided 2026-09-03; OPEN_ITEMS_SOLUTIONS.md section 10, condition 4).
  /// The shipped observable is a RATIO -- the tensor weight is K/D_phi, and
  /// D_phi's F1 is `InclusiveKernel`'s, built from the SAME `ToyF2` this
  /// object is handed as `unpol`, whose R is `r_sigma_lt` because nothing sets
  /// it -- so a different R here would not cancel:
  /// (1 + r1998)/(1 + r_sigma_lt) is 1.088 at x = 0.1, Q2 = 2.5.
  ///
  /// MEASURED cost of the choice on x*b1 (Q2 = 2.5, defaults): -5.5 % /
  /// +24.6 % / +5.4 % / +2.7 % / -1.9 % at x = 0.05 / 0.10 / 0.20 / 0.30 /
  /// 0.50.  ALL of it is in the two orbital terms: term (1) is BIT-IDENTICAL
  /// under the swap, because it carries b1_d from the injected `TensorSF`,
  /// which has no R in it.  There the swap is x0.54 to x0.94, far more than
  /// the <= 8 % a (1 + R) prefactor allows, because those terms convolve F1
  /// against a density that INTEGRATES TO ZERO and so respond to the SLOPE of
  /// R -- `r_sigma_lt` is x-independent, `r1998` runs 0.30 -> 0.20 over
  /// x = 0.05 -> 0.5.  On the A = 2 gate, where there is no such
  /// cancellation, the same swap is only +2.7 % on G3b.
  ///
  /// THE CLEAN FIX EXISTS SINCE 2026-09-06, AS AN OPT-IN AND NOT AS A NEW
  /// DEFAULT: `PipelineConfig::r_source` (`RSource`, pipeline.hpp; CLI
  /// `--r-source`) threads ONE R object through `default_inclusive_kernel`
  /// into THIS field and `InclusiveKernel::Options::r_func` together, so the
  /// ratio's two halves cannot disagree.  Its default `Unset` does not enter
  /// that branch at all, which is why this field's own null default is still
  /// what every shipped number is made with.  The option is PRICED in
  /// docs/open_items/run_2026-09-06/phase_A_numbers.md sec. A1 and it does
  /// NOT simply undo the mismatch: at y = 0.5, Q2 = 2.5 sharing `r1998`
  /// moves K/D_phi by -3.8357 % / +26.3988 % / +3.3518 % at
  /// x = 0.05 / 0.10 / 0.30 against the numerator-only swap's
  /// -5.5148 % / +24.6099 % / +2.6629 % -- LARGER at x = 0.10 -- and the
  /// shared shift is y-DEPENDENT where the numerator-only one is not
  /// (-5.4705 % at y = 0.1 to +2.3678 % at y = 0.9 for x = 0.05, Q2 = 2.5).
  RFunc r_func;                        ///< null => r_sigma_lt

  /// TERM KNOBS -- each multiplies one term of design 1.7.  All 1 = nominal.
  double w_embedded_s = 1.0;           ///< term (1)
  double w_alpha_d_dwave = 1.0;        ///< terms (2d) AND (2a) <- THE band knob
  double w_cg_dwave = 1.0;             ///< term (3), on top of the exact 1/10

  /// Normalisation.  `use_spectroscopic_factor` = true (default) keeps
  /// `norm_target` = VMC_N_ALPHA_D_LI6 = 0.819481, i.e. the non-alpha-d 18 %
  /// of 6Li is given b1 = 0.  false renormalises the density to 1 (the "all
  /// of 6Li is alpha + d" variant); the two differ by 1/N_ad = 1.22 on all
  /// four terms.  `norm_target` overrides the value, which is how the +-5 %
  /// N_ad systematic of design 2.1 is quoted.
  bool use_spectroscopic_factor = true;
  double norm_target = 0.0;            ///< 0 => VMC_N_ALPHA_D_LI6

  /// CDKS Eq. (21)'s finite-|q| delta-function instead of Eq. (17)'s
  /// kappa = 1 (design 1.3).  FALSE by default -- design A10's stated choice,
  /// and here it costs almost nothing: in 6Li the switch moves x b1 by -1 % to
  /// +7 %, because the term kappa multiplies is the small ORBITAL one and
  /// M_d, M_alpha are large.  (On the A = 2 gate it is worth 1.5-1.7 at
  /// x >= 0.5, which is why `DeuteronConvolutionB1::Options::finite_q_delta`
  /// defaults the other way: that object exists to reproduce CDKS's own
  /// figure, and Eq. (21) is their exact definition.)  true rebuilds both
  /// density sets per x with kappa = sqrt(1 + gamma^2(x, q2)).
  /// THE DEFAULT IS A CHOICE, NOT A DERIVATION.
  bool finite_q_delta = false;

  /// Quadrature.  `norm_target` / `renormalize` are overwritten from the two
  /// fields above.
  LightConeDensities::Options quad;
};

/// b1 of 6Li per nucleon, four-term alpha-d convolution.  OPT-IN; the library
/// default stays `Li6B1(MillerB1)` (pipeline.cpp).
class Li6ConvolutionB1 : public TensorSF {
 public:
  Li6ConvolutionB1();
  explicit Li6ConvolutionB1(Li6ConvolutionOptions opt);
  /// Test-only: build on prepared densities (e.g. `delta_limit`), T6.  It
  /// still runs the full constructor first (so the VMC data files must be
  /// present) and then replaces the two density objects.
  ///
  /// THROWS on `finite_q_delta`: the injected densities would be silently
  /// discarded at the first `b1()` call, because the finite-|q| path rebuilds
  /// both sets per x from the FILE waves (`w0_`, `w2_`), which this
  /// constructor does not replace.  T5's zero-D-wave densities would grow a D
  /// wave back.  Refusing is the honest reading; the alternative (storing the
  /// injected waves alongside) would be a different object.
  Li6ConvolutionB1(Li6ConvolutionOptions opt, LightConeDensities d_struck,
                   LightConeDensities a_struck);

  double b1(double x, double q2, double f1) const override;
  /// b2 = 2 x b1 (base class), delta = 0 (base class).

  /// The four terms separately, same units and conventions as `b1`.
  double b1_embedded_s(double x, double q2, double f1) const;         ///< (1)
  double b1_alpha_d_dwave_d(double x, double q2, double f1) const;    ///< (2d)
  double b1_alpha_d_dwave_alpha(double x, double q2, double f1) const;///< (2a)
  double b1_cg_dwave(double x, double q2, double f1) const;           ///< (3)
  /// (2d) + (2a), the physical "alpha-d orbital" term.
  double b1_alpha_d_dwave(double x, double q2, double f1) const;
  /// SD / DD split of (2d) and (2a), for the Fig.-4-style plot.
  double b1_alpha_d_sd(double x, double q2, double f1) const;
  double b1_alpha_d_dd(double x, double q2, double f1) const;
  /// P2/P4 remainders of design 1.5 truncation item 1 -- REPORTED, not added
  /// to `b1`; T14 records their ratio to term (1).
  double b1_cg_dwave_p2_remainder(double x, double q2, double f1) const;
  double b1_cg_dwave_p4_remainder(double x, double q2, double f1) const;

  /// A scaled copy for the mandatory 100 % band: `scale` in {0, 1, 2}
  /// multiplies the WHOLE b1.  Cheap: shares the densities.
  std::shared_ptr<const TensorSF> banded(double scale) const;

  const LightConeDensities& densities() const;        ///< struck deuteron
  const LightConeDensities& densities_alpha() const;  ///< struck alpha
  const Li6ConvolutionOptions& options() const;
  double band_scale() const;
  /// The F1 slots the convolution actually calls, exposed so that a test (and
  /// a plot legend) can see WHICH F1 went in.  `f1_deuteron` is CDKS Eq. (22)
  /// on F2_d/2 and must exceed `NuclearF2::f1a`/2 by exactly 1 + gamma^2.
  double f1_deuteron(double x, double q2) const;
  double f1_alpha(double x, double q2) const;
  /// The deuteron b1 camp this instance convolves.
  const std::shared_ptr<const TensorSF>& deuteron_b1() const;
  /// int b1 dx over [lo, hi] -- the Close-Kumano number, REPORTED not
  /// enforced, exactly like `close_kumano_integral`.
  double close_kumano_integral(double lo = 0.01, double hi = 1.2,
                               double q2 = 2.5, std::size_t n = 241) const;

 private:
  struct Terms { double s, dd_d, dd_a, cg; };
  Terms terms(double x, double q2) const;
  const LightConeDensities& dens_d(double x, double q2) const;
  const LightConeDensities& dens_a(double x, double q2) const;

  Li6ConvolutionOptions opt_;
  ClusterPartialWave w0_, w2_;
  std::shared_ptr<const LightConeDensities> d_, a_;
  std::shared_ptr<const TensorSF> b1_d_;
  /// The support of `b1_d_`, resolved from `Li6ConvolutionOptions`.
  double b1_d_x_max_ = 1.0;
  std::function<double(double, double)> f1_d_, f1_alpha_;
  double band_ = 1.0;
  /// finite_q_delta only: the last (x, Q2) the densities were rebuilt for.
  mutable std::shared_ptr<const LightConeDensities> qd_, qa_;
  mutable double q_x_ = -1.0, q_q2_ = -1.0;
};

/// Which deuteron wave function the A = 2 gate convolves.
///
/// `kFdeutFile` is the DEFAULT and the pre-existing behaviour, bit for bit:
/// the tabulated AV18 file at `Options::fdeut_path`.  `kCdBonn` swaps in the
/// analytic CD-Bonn parameterisation (`cdbonn_wave`, cluster.hpp), which is
/// what CDKS actually used for Fig. 4 and what gate condition 3 of open
/// item 10 asks for.  The DEFAULT IS A CHOICE and it is the conservative one:
/// every published gate number in docs/open_items/ was measured on AV18, so
/// moving the default would silently move all of them.  The two rows are
/// measured side by side in
/// docs/open_items/run_2026-09-03/phase_A_numbers.md section 8.
///
/// The D-state RESCALING PROXY of open item 5 (AV18's w(k) scaled to
/// P_D = 4.85 %, S renormalised) is NOT an entry here because it never had a
/// committed implementation -- the 2026-09-02 phase-D work built it in
/// scratch by writing an `fdeut`-format file.  It stays exactly as reachable
/// as it was, through `fdeut_path`, and phase_A_cdbonn.md section 7 records
/// the two scale factors that reproduce it.
enum class DeuteronWaveSource {
  kFdeutFile,   ///< the tabulated file at `fdeut_path` (AV18) -- THE DEFAULT
  kCdBonn,      ///< the analytic CD-Bonn parameterisation of cluster.hpp
};

/// A = 2 VALIDATION GATE.  The same kernel fed the AV18 deuteron u(k), w(k)
/// and a nucleon F1: this MUST reproduce CDKS Fig. 4 (`tables::kB1CdksQ2p5`)
/// before any 6Li number is quoted.  It DOES, since 2026-09-03 -- see design
/// section 5, this file's header block for the verdict and its four
/// conditions, and docs/open_items/run_2026-09-03/phase_A_numbers.md for the
/// measurements (docs/open_items/run_2026-09-02/phase_D_gate.md is the
/// superseded failing verdict, kept for its argument).
///
/// The F1 slot is the ISOSCALAR nucleon (F1p + F1n)/2 with R = `r1998`
/// (CDKS's SLAC world fit) and the CDKS Eq. (22) target-mass factor, NOT
/// F1 of the deuteron: F1_d/F1_N differs from 1 by <= 1 % below x = 0.6, so
/// this choice is not where the uncertainty lives, but it is CDKS's.
class DeuteronConvolutionB1 : public TensorSF {
 public:
  struct Options {
    std::string fdeut_path = data_path("vmc/deuteron/fdeut.av18");
    /// WHICH WAVE FUNCTION.  `kFdeutFile` (the default) reads `fdeut_path`;
    /// `kCdBonn` IGNORES `fdeut_path` and builds CD-Bonn analytically on the
    /// grid below.  See `DeuteronWaveSource`.
    DeuteronWaveSource wave = DeuteronWaveSource::kFdeutFile;
    /// `kCdBonn` only: the fm^-1 grid the analytic form is sampled on.  The
    /// defaults are `fdeut.av18`'s OWN grid (0 .. 20 fm^-1 in 0.1, 201 rows),
    /// so an AV18 row and a CD-Bonn row differ in the wave function and in
    /// NOTHING else -- same spline, same node count, same truncation.
    /// Refining to 0.02 moves the gate peak by 5e-4 relative
    /// (phase_A_cdbonn.md section 8.4).  Ignored by `kFdeutFile`.
    double cdbonn_k_max_fm = 20.0;
    double cdbonn_dk_fm = 0.1;
    std::shared_ptr<const UnpolSF> unpol;   ///< null => ToyF2
    RFunc r_func;                           ///< null => r1998 (CDKS's choice)
    bool target_mass = true;                ///< CDKS Eq. (22) factor
    /// CDKS Eq. (21)'s delta-function, y = p.q/(M_N nu) = (E - p_z kappa)/M_N
    /// with kappa = |q|/nu = sqrt(1 + gamma^2).
    ///
    /// TRUE HERE, and false on `Li6ConvolutionOptions` -- the two defaults
    /// differ ON PURPOSE (design 1.3, A10, and the review of 2026-09-03).
    /// Eq. (21) is CDKS's EXACT definition and Eq. (18)'s "~=" is what makes
    /// Eq. (17)'s (E - p_z)/M_N the approximation to it; on the A = 2 gate,
    /// which exists to reproduce THEIR figure, kappa is worth x1.62 at the
    /// peak and x1.97 on int b1 dx.  It is load-bearing for the verdict, not
    /// only for the third digit: with CDKS's own MSTW2008 LO the G3b ratio is
    /// 0.843 at Eq. (21) against 0.520 at Eq. (17), so the gate passes
    /// comfortably at their exact definition and by 4 % at the "~=" of it.
    /// In 6Li the same
    /// switch is a 1 % effect (the term kappa multiplies is the small orbital
    /// one) and the default there stays Eq. (17) as printed, which is what
    /// keeps f(y) a function of y alone -- one cached table reused at every x.
    ///
    /// COST: true rebuilds the densities at EVERY x (kappa depends on x), so a
    /// point costs ~0.24 s against ~0.7 ms.  The cache is one deep.
    bool finite_q_delta = true;
    LightConeDensities::Options quad;
  };
  DeuteronConvolutionB1();
  explicit DeuteronConvolutionB1(Options opt);
  double b1(double x, double q2, double f1) const override;
  double b1_sd(double x, double q2, double f1) const;
  double b1_dd(double x, double q2, double f1) const;
  /// The isoscalar F1 the convolution is fed (CDKS Eq. 22), for the gate's
  /// G1f row and for the checklist.
  double f1_nucleon(double x, double q2) const;
  const LightConeDensities& densities() const;
  const Options& options() const;

 private:
  const LightConeDensities& dens(double x, double q2) const;
  Options opt_;
  std::shared_ptr<const LightConeDensities> d_;
  ClusterPartialWave phi0_, phi2_;
  double eps_ = 0.0;
  mutable std::shared_ptr<const LightConeDensities> qd_;
  mutable double q_x_ = -1.0, q_q2_ = -1.0;
};

/// The A = 2 gate's landmarks, COMPUTED (never typed): the sign changes with
/// their slope, the extremum positions and values of x*b1, and int b1 dx.
/// Both the doctest and `validation/b1_li6_table.py` call this, and it is also
/// applied to the digitized table itself so G3a/G3b compare like with like.
struct B1Landmarks {
  std::vector<double> zeros;   ///< ascending, linearly interpolated crossings
  std::vector<int> zero_slope; ///< +1 rising, -1 falling
  double x_min = 0.0, xb1_min = 0.0;
  double x_max = 0.0, xb1_max = 0.0;
  double integral_b1 = 0.0;    ///< trapezoid of (x*b1)/x over `x_grid`
};
B1Landmarks b1_landmarks(const std::function<double(double)>& xb1,
                         const std::vector<double>& x_grid);
B1Landmarks b1_landmarks_of_table();   ///< tables::kB1CdksQ2p5(), raw column

/// alpha-d contribution to the 6Li quadrupole moment from the SAME phi_0,
/// phi_2 -- the SIGN GATE of design 2.2.  Must come out NEGATIVE.
///
/// The waves are Fourier-transformed back to r (u_L(r) = r sqrt(2/pi) int
/// k^2 j_L(kr) phi_L(k) dk, the transform inverse to CDKS's Eq. (20) pair)
/// and fed the standard quadrupole operator
///   Q = (1/20) int dr r^2 w(r) [2 sqrt2 u(r) - w(r)],   u = phi_0, w = -phi_2,
/// after normalising to int k^2 (|phi_0|^2 + |phi_2|^2) dk = 1, so the return
/// value is the quadrupole moment of the UNIT-NORMALISED relative wave
/// function (multiply by N_ad for the alpha-d contribution to Q(6Li)).
///
/// DEVIATION from design 3, which writes "Q_ad ~ -(sqrt2/5) <r^2>_02" -- the
/// 02 term alone.  The full two-term operator is used because it is the
/// textbook one and because it makes the machinery TESTABLE: on the deuteron
/// pair it must return `fdeut.av18`'s own header value qm = 0.269673 fm^2,
/// which it does to 0.4 %.  The 02 term dominates and carries the sign either
/// way.  fm <-> GeV with `HBARC_GEV_FM` (design 2.3).
double alpha_d_quadrupole_fm2(const ClusterPartialWave& phi0,
                              const ClusterPartialWave& phi2,
                              double r_max_fm = 20.0, std::size_t n_r = 2000,
                              std::size_t n_k = 4001);

/// THE A = 7 GATE, and the ONLY validation 7Li gets offline.
///
/// `Q = -(2/5) Z_eff <r^2>` for the alpha + t cluster state, Eq. (6) of
/// `docs/open_items/run_2026-09-03/phase_D_li7_rank2.md` sec. 4.2, computed
/// from the SAME `momenta/li7_at3.momentum` wave `li7_alpha_channel` already
/// ships to the tagged channel.  Every factor is parameter-free:
///
///   * both clusters are spherical (alpha is 0+, t is 1/2+), so the WHOLE of
///     Q(7Li) is orbital in this picture -- the same statement that makes
///     7Li's b1 purely orbital;
///   * for a pure L = 1, S = 1/2, J = 3/2 state
///     <3cos^2(theta) - 1>_{M=3/2} = 2[L(L+1) - 3 m_L^2]/[(2L-1)(2L+3)]
///     = -2/5 at m_L = +1, which is twice the SAME -1/5 the 7Li polarimeter
///     already uses as <P_2(cos theta_k)> = -T/5;
///   * `Z_eff` = Z_alpha (M_t/M_7)^2 + Z_t (M_alpha/M_7)^2 = 0.695075 with
///     the AME nuclear masses `nuclear_mass(1,3) / (2,4) / (3,7)` (the
///     A-number form 34/49 = 0.693878 is 0.17 % away and is NOT used);
///   * <r^2> comes from u_1(r) = r sqrt(2/pi) int k^2 j_1(kr) phi_1(k) dk,
///     the L = 1 member of the same Fourier-Bessel pair
///     `alpha_d_quadrupole_fm2` uses at L = 0, 2, and is normalised by its own
///     int u^2 dr -- so no norm, no spectroscopic factor and no global phase
///     enters, which is why `VMC_S_ALPHA_T_LI7` = 1.0084 > 1 does not
///     contaminate it.
///
/// WHAT IT VALIDATES, AND WHAT IT DOES NOT.  It validates the alpha-t WAVE
/// FUNCTION's quadrupole -- its <r^2> and its P-wave character -- against a
/// measured rank-2 observable of the same nucleus, `LI7_QUADRUPOLE_FM2`
/// (rc.hpp).  It validates NOTHING about the light-cone convolution, the
/// unpolarised DIS input, or b1(7Li) -- which is NOT IMPLEMENTED (open item
/// 15 of docs/OPEN_ITEMS_SOLUTIONS.md; 7Li's rank-2 sector is exactly zero by
/// construction) and whose sign is not even determined by the inputs in this
/// tree.  There is no A = 3 analogue of the A = 2 b1 gate and there cannot be
/// one: 3H and 3He are J = 1/2 and have no rank-2 structure function at all
/// (`TaggedModel::tensor_dilution` throws for exactly that reason), so this
/// quadrupole is the only offline check the 7Li wave function will get.
///
/// MEASURED (phase_B_numbers.md sec. B3, and pinned in tests):
///   <r^2> = 12.534829 fm^2, r_rms = 3.540456 fm, Z_eff = 0.695075,
///   Q = -3.485059 fm^2 at the defaults below;
///   grid: -3.472004 (r_max 20) / -3.485059 (30) / -3.485797 (40) /
///         -3.486171 (60), i.e. converged to 0.1 % at 30 fm;
///   VMC MC band (`vmc_mc_sigma` = -+1, fully correlated):
///         -3.496427 / -3.485059 / -3.473745, i.e. +-0.33 %;
///   ratio to `LI7_QUADRUPOLE_FM2` = -4.06: 0.858389 (14.2 % low).
/// It is a REPORTED ratio, never a pass/fail on b1.  For contrast, the same
/// construction one cluster level down gives
/// `alpha_d_quadrupole_fm2(li6_alpha_d_partial_waves())` = -0.3333 fm^2
/// against a measured -0.0818 -- a factor 4.08 -- so the 7Li wave-function
/// input is validated by a measured moment an order of magnitude better than
/// 6Li's is, and that asymmetry is the point of the gate.
struct AlphaTQuadrupole {
  double r2_fm2 = 0.0;    ///< <r^2> of the alpha-t relative motion [fm^2]
  double r_rms_fm = 0.0;  ///< sqrt(r2_fm2) [fm]
  double z_eff = 0.0;     ///< Z_alpha (M_t/M_7)^2 + Z_t (M_alpha/M_7)^2
  double q_fm2 = 0.0;     ///< -(2/5) z_eff r2_fm2 [fm^2]
};

/// The gate on an arbitrary L = 1 relative wave.  `phi1` is read for its
/// (k, psi) table only -- the L it carries is not consulted, j_1 is applied
/// here -- and the k grid is the table's own, in GeV.
AlphaTQuadrupole alpha_t_quadrupole(const VmcRadial& phi1,
                                    double r_max_fm = 30.0,
                                    std::size_t n_r = 4000,
                                    std::size_t n_k = 8001);

/// The gate on 7Li's own shipped wave: `li7_alpha_channel(BETA_DEFAULT,
/// ClusterWaveSource::VmcAV18, vmc_mc_sigma)`, i.e. `momenta/li7_at3.momentum`
/// (the 3/2- GROUND state).  `vmc_mc_sigma` is the ANL table's own fully
/// correlated 1-sigma band, 0 by default and then bit for bit today's table.
AlphaTQuadrupole li7_alpha_t_quadrupole(double vmc_mc_sigma = 0.0,
                                        double r_max_fm = 30.0,
                                        std::size_t n_r = 4000,
                                        std::size_t n_k = 8001);

}  // namespace lipolgen

#endif  // LIPOLGEN_B1_NUCLEAR_HPP
