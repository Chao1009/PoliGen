#ifndef LIPOLGEN_CLUSTER_CONFIG_HPP
#define LIPOLGEN_CLUSTER_CONFIG_HPP

/// \file cluster_config.hpp
/// Nucleon-position CONFIGURATIONS of a polarized 6Li in the alpha + d cluster
/// picture -- the input a Good-Walker dipole-model code (subnucleondiffraction,
/// arXiv:2408.13213) averages the coherent amplitude over.
///
/// OPT-IN AND INERT.  Nothing in the generator calls this header.  It adds no
/// default, changes no existing number, and its output is a file.  See
/// docs/open_items/run_2026-09-02/design_G_cluster_config.md for the physics,
/// the validations and the (large) wave-function systematic on the tensor
/// observable; docs/open_items/run_2026-09-02/phase_G_numbers.md carries the
/// reproduced numbers.
///
/// ------------------------------------------------------------------ physics
///
/// 6Li(1+) is a rigid alpha(0+) core plus a deuteron with relative orbital
/// L = 0, 2 coupled to the deuteron spin S = 1.  For ANY (L = 0, 2) x (S = 1)
/// -> J = 1 system with radial functions f_0, f_2 normalized as
/// integral (f_0^2 + f_2^2) x^2 dx = 1, the density of the relative vector in
/// substate m with the quantization axis along +z is
///
///   rho_m(x, cos theta) = sum_{m_S} | sum_L C_L(m, m_S) f_L(x)
///                                     Theta_L^{m-m_S}(cos theta) |^2      (G1)
///   C_L(m, m_S) = <L (m-m_S) 1 m_S | 1 m> ,
///
/// with no phi dependence.  `cluster_density` is (G1) and `cluster_amp2` is
/// one m_S branch of it (G1'); both go through `clebsch_gordan` (spin.hpp) and
/// `theta_lm` (cluster.hpp) -- no hand-expanded closed form lives in the
/// library.  (G1) with (u/r, w/r) IS arXiv:2408.13213 Eqs. (7)-(8), which
/// tests/test_cluster_config.cpp checks to 1e-14.
///
/// THE (-i)^L TRAP, and it is not optional bookkeeping (design sec. 2.1).
/// This module works in r SPACE and applies NO phase: the r-space signs of
/// R_0, R_2 are taken verbatim from the tables.  The physical k-space partial
/// waves are psi_L(k) = (-i)^L 4 pi integral j_L(kr) R_L(r) r^2 dr, so the
/// PHYSICAL k-space S-D interference carries (-i)^2 = -1 and has the OPPOSITE
/// sign to the ANL `li6.ad` k-block, to `VmcRadial` and to
/// docs/open_items/vmc_reconciliation.md -- all three of which are no-phase
/// transforms, exactly like this one.  The only OBSERVABLE check on the
/// r-space signs is sign(Q) == sign(LI6_QUADRUPOLE_FM2) (test T4(iii)); that
/// the no-i^L transform of (R_0, R_2) reproduces the k-block's node at
/// 0.678 fm^-1 and its sign pattern is T4(iv).  DO NOT "fix" either
/// convention to match the other: they differ by exactly this phase.
///
/// Integrating (G1) against (3 cos^2 theta - 1),
///
///   Q[f_0,f_2] = <3 x_z^2 - x^2>_{m=+1}
///              = (1/5) integral x^4 ( 2 sqrt2 f_0 f_2 - f_2^2 ) dx        (G4)
///
/// and <3 x_z^2 - x^2>_m = (3m^2 - 2) Q.  With the alpha at -R/3 + s_i, the
/// proton at +2R/3 + r/2 and the neutron at +2R/3 - r/2 (total c.m. at the
/// origin), the six-nucleon moments are EXACT one-coordinate sums:
///
///   Q_matter(6Li, M) = (3M^2-2) [ (4/3) Q[R_0,R_2] + 2 Q_d D_T ]          (G5)
///   D_T = <3 m_S^2 - 2>_{M=+1} = 1 - 0.9 P_D(alpha-d)   (exact CG identity)
///   <r^2>(6Li) = (2/3)<s^2>_alpha + (1/12)<r_np^2> + (2/9)<R^2>_alphad    (G6)
///
/// For a TRANSVERSE quantization axis the per-nucleon transverse asymmetry and
/// its cos 2Phi coefficient follow with no free parameter:
///
///   delta = <x^2> - <y^2> = Q_matter(m) / (2 A) * sin^2(theta_s) cos(2 phi_s)
///   a_2(m) = -(delta_m / 4) |t|            [arXiv:2408.13213 Eq. (9)]     (G8)
///
/// -- `a2_from_quadrupole`, which reproduces the published polarized-deuteron
/// a_2(+-1) to 8 % at every |t| with zero free parameters.
///
/// ------------------------------------------------------- the approximations
///
/// * DIAGONAL TRUNCATION.  The joint density's m_S coherences are dropped, so
///   only the RELATIVE azimuth of R-hat and r-hat is lost.  This is the same
///   truncation `tagged.hpp` makes and documents.  The conditioning is NOT
///   dropped: (R, cos theta_R) is drawn from the m_S-CONDITIONED table
///   |A_{m_S}(M; R, c)|^2, exactly as `TaggedModel::sample_kc(m_ion, m_s, ...)`
///   draws from `build_amp2`'s per-m_S table.  Drawing R-hat from the m_S-
///   SUMMED density instead is a different model that no moment test can see;
///   `rho_alpha_d_summed` exists for plots and must not be sampled from.
///   `ClusterConfigOptions::exact_coherence` is reserved and THROWS.
/// * The tensor output carries a factor-of-several wave-function systematic:
///   this geometry gives Q_charge(6Li) = -0.615 fm^2 (model range
///   -0.615 .. -0.730) against the MEASURED -0.0818 (`LI6_QUADRUPOLE_FM2`,
///   rc.hpp) and GFMC AV18+IL7's -0.20(6) (`LI6_QUADRUPOLE_GFMC_FM2`, THIS
///   header).  `quadrupole_band_fm2()` returns all three and the writer stamps
///   them; docs/OPEN_ITEMS_SOLUTIONS.md sec. 11's rule -- do not derive a
///   published tensor input from these wave functions -- stands.
///   The 7.52x is TWO factors, not three: 3.317 to the eta-matched dial and
///   2.269 from there to the measurement (T22b).  1/S_alpha-d = 1.171 is NOT
///   the third -- both waves are divided by sqrt(s_alpha_d_) before the
///   moments are taken, so it is already inside the -0.615.
///
/// ------------------------------------------------------- units and hbar c
///
/// Positions are in fm, in the ION REST FRAME, with the total c.m. at the
/// origin.  EVERY fm <-> GeV conversion in this module goes through
/// `HBARC_GEV_FM` (constants.hpp, 0.19733) and never through a literal, so
/// that `eps_b0_equivalent()` is strictly comparable with `gaussian_slope`
/// and `CoherentScenario::slope_b`.  `cluster.cpp`'s file-local
/// `kHbarCGeVfm = 0.1973269804` is a pre-existing second copy and is NOT
/// propagated here.
///
/// EXTENT.  Every tabulated radial function used here is IDENTICALLY ZERO
/// beyond its last abscissa -- `VmcRadial::operator()`'s own rule, and NOT
/// `np_interp`'s, which clamps to the end value (numerics.hpp).  The two
/// defaults collide and the difference is not small: carrying
/// `li6.adr.fit`'s last value R_0(9.95) = -0.0039 flat out to r_max_fm = 20
/// would move <R^2> from 16.97 to 27.39 fm^2.  `r_max_fm` / `rnp_max_fm` are
/// CEILINGS on the sampling grid, never extensions of a table.

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "lipolgen/cluster.hpp"
#include "lipolgen/rng.hpp"

namespace lipolgen {

// --------------------------------------------------------------- constants

/// RNG `bunch` slot owned by this module.  The generator uses the spin
/// CATEGORY INDEX (0, 1, 2, ...) as its bunch slot (generator.hpp), so a
/// configuration stream must sit far away from small integers.  Bound to
/// Python as `CONFIG_STREAM`; the sidecar's "rng":{"stream"} records it.
inline constexpr std::uint64_t kConfigStream = 0x434F4E464947ull;  // 'CONFIG'

/// r_point(6Li) from the ANL VMC one-body density `li6.density` (AV18+UX),
/// header "SQRT(4*PI*TOTINT(RHORP*R**4:R)/3) = 2.4433".  The DEFAULT target of
/// `match_li6_radius()`; a MODEL number, not a measurement -- the measurement
/// is `rc.hpp`'s `LI6_R2_POINT_FM2` = 6.0788 fm^2 -> 2.4655 fm.
inline constexpr double LI6_R_POINT_VMC_FM = 2.4433;

/// Q(6Li) from GFMC with AV18+IL7: -0.20(6) e fm^2.  Pastore, Pieper,
/// Schiavilla & Wiringa, PRC 87, 035503 (2013), arXiv:1212.3375; the repo's
/// own record is docs/open_items/physics_literature.md.  NOT VMC and NOT the
/// AV18+UX Hamiltonian of `li6.density` / `li6.ad`, so it is an
/// independent-Hamiltonian comparison point, not a check on our own tables.
/// The middle entry of `quadrupole_band_fm2()`; nowhere else is it retyped.
inline constexpr double LI6_QUADRUPOLE_GFMC_FM2     = -0.20;
inline constexpr double LI6_QUADRUPOLE_GFMC_ERR_FM2 = 0.06;

/// The MEASURED asymptotic D/S ratio of the alpha-d channel,
/// eta = C_2/C_0 = -0.025 +- 0.006 (stat) +- 0.010 (syst).  George & Knutson,
/// PRC 59, 598 (1999), from d + alpha elastic tensor analysing powers; other
/// determinations lie in -0.01 .. -0.03.  The repo's own record is
/// docs/PHYSICS_CHANNELS.md's reference [GK99] -- NOT physics_literature.md,
/// which predates this citation and does not carry it.  The comparison point
/// for `asymptotic_ds_ratio()`; nowhere else is it retyped.
///
/// READ THE BAND, NOT THE CENTRAL VALUE.  `asymptotic_ds_ratio()` is EXACTLY
/// LINEAR in `quadrupole_dial_s()` (the dial multiplies R_2 by s and both
/// waves by the same 1/sqrt(n(s)), so the wave RATIO carries s alone), which
/// makes eta a one-to-one relabelling of the dial rather than an independent
/// constraint on it.  Mapped through that line, +-1 sigma_comb = +-0.01166
/// spans model Q_charge from -0.4005 fm^2 to +0.0298 fm^2 -- through ZERO, at
/// eta = -0.01497, only 0.86 sigma from the central value.  Anything derived
/// from the central eta alone is a point on a line that the error bar does not
/// even pin the SIGN of.  Test T22b pins both edges.
inline constexpr double LI6_ETA_DS_GK      = -0.025;
inline constexpr double LI6_ETA_DS_GK_STAT = 0.006;
inline constexpr double LI6_ETA_DS_GK_SYST = 0.010;

/// The data tables this module reads, relative to `data_dir()`.  One home
/// each; `tagged.hpp`'s `VMC_LI6_OVERLAP` is the fifth and is NOT retyped.
inline const char* const VMC_HE4_DENSITY   = "vmc/density/he4.density";
inline const char* const VMC_LI6_DENSITY   = "vmc/density/li6.density";
inline const char* const VMC_LI6_AD_FIT    = "vmc/li6_alpha_d/li6.adr.fit";
inline const char* const VMC_DEUTERON_WAVE = "vmc/deuteron/fdeut.av18";

/// The R window [fm] `asymptotic_ds_ratio()` averages over.  A CHOICE: below
/// ~5 fm the Whittaker asymptotics have not set in, above ~9 fm the fit table
/// has ended and the raw block's D wave has run out of signal.  MEASURED on
/// `li6.ad`'s own error column (2026-09-03), median (min) S/N per window:
/// 6-8 fm    R_0 37.7 (24.9), R_2 13.5 (8.9)  -- NOT noise, this is the window
/// 8-9.95 fm R_0 15.5 (9.0),  R_2  6.8 (3.0)
/// 9.95-12   R_0  5.1 (2.8),  R_2  2.8 (1.4)  -- here it is.
inline constexpr double kDsRatioRLoFm = 6.0;
inline constexpr double kDsRatioRHiFm = 8.0;

// ------------------------------------------------------------------ sources

/// Where the alpha core's one-body point-nucleon density comes from.
enum class AlphaCoreSource {
  /// `data/vmc/density/he4.density` (ANL VMC AV18+UX, 500k samples).  DEFAULT.
  VmcHe4Density = 0,
  /// Gaussian rho ~ exp(-r^2/2a^2) with a^2 = `cluster_point_a2_fm2(2, 4)`
  /// = 0.7001 fm^2 -- the SAME single copy the FSI weight uses (fsi.hpp).
  /// rms 1.4492 fm against the table's 1.4413: a 0.5 % size difference and a
  /// visibly different surface (the VMC alpha has a flat-topped interior).
  Gaussian = 1,
};

/// Where the alpha-d relative radial functions R_0(R), R_2(R) come from.
enum class AlphaDSource {
  /// `li6.adr.fit` shapes, each wave rescaled to `li6.ad`'s own N_0, N_2
  /// (s_0 = 1.0128, s_2 = 0.8984).  DEFAULT: the raw block's A_2 column
  /// changes sign four times below 1.2 fm at the Monte Carlo noise level,
  /// which a CDF sampler cannot use, while the fit's own norms are wrong.
  FitRescaled = 0,
  FitRaw      = 1,   ///< `li6.adr.fit` as published (P_D = 0.02542)
  OverlapRaw  = 2,   ///< `li6.ad` r-block verbatim -- noisy below ~1.2 fm
};

// ----------------------------------------------------------------- options

/// Every field here appears in the sidecar's "options" object, and
/// tests/test_cluster_config.cpp pins the list, so adding one without
/// recording it fails.
struct ClusterConfigOptions {
  AlphaCoreSource alpha_source   = AlphaCoreSource::VmcHe4Density;
  AlphaDSource    alpha_d_source = AlphaDSource::FitRescaled;

  /// Quantization axis in the ion rest frame [rad]; the SAME axis as
  /// spin.hpp / tagged.hpp, applied as R_z(phi_s) R_y(theta_s) exactly as
  /// `boost_spectator` does.  A cos 2Phi signal needs theta_s != 0: with the
  /// axis along the beam the projected density is azimuthally symmetric and
  /// a_2 == 0 (arXiv:2408.13213 Fig. 2(c)).
  double theta_s = 1.5707963267948966;   ///< default: transverse (pi/2)
  double phi_s   = 0.0;

  /// Inflate the alpha's SOURCE density by lambda = (1 - 1/4)^-1/2 = 1.154701
  /// in radius before recentring.  What this does and does NOT do:
  ///   * <s^2> closes EXACTLY, for ANY source shape -- recentring four
  ///     independent draws gives <s^2> = (3/4)<v^2>, so lambda^2 = 4/3 puts
  ///     the recentred second moment back on the table's.
  ///   * the SHAPE is NOT restored for a tabulated density; it IS exact in
  ///     shape for `AlphaCoreSource::Gaussian`.  The recentred density has a
  ///     closed form (`rho_alpha_recentred`) which is what the tests gate the
  ///     sampled histogram against -- NOT `he4.density` itself.
  bool alpha_cm_inflate = true;

  /// Minimum nucleon-nucleon separation inside the alpha core [fm].  0 = off
  /// (the default: an uncorrelated product of one-body densities).  0.9 is
  /// what arXiv:2605.00454 imposes on its Woods-Saxon sampling.
  ///
  /// WHAT IT COSTS THE ANALYTIC LAYER.  Rejection correlates the four s_i, so
  /// the closed form <s^2> = (3/4)<v^2> -- exact for INDEPENDENT draws, and
  /// the whole content of `alpha_cm_inflate` -- no longer holds: at 0.9 fm the
  /// sampled <r^2> moves by about +2 %.  The constructor therefore MEASURES
  /// <s^2> once, deterministically (1e5 cores on the fixed stream
  /// Rng(0, 0, kConfigStream, 0)), stamps "alpha <s^2> NUMERICAL (hard core)"
  /// into `provenance()`, and the sidecar sets
  /// `moments.r2_mean_is_approximate`.  `r2_analytic_fm2()`,
  /// `match_li6_radius()` and `r2_grid_fm2()` all use that measured value and
  /// are approximate at its precision (~2e-3 fm^2 on <r^2>).
  /// `rho_alpha_recentred()` is NOT corrected -- it stays the independent-draw
  /// closed form, so it and test T5(b) hold only at 0.
  double min_nn_separation_fm = 0.0;

  /// Multiply the alpha-d separation R by this before assembling -- a pure
  /// coordinate rescaling, so <R^2> and Q[R_0,R_2] both move by its square.
  /// 1.0 = the wave functions as published (DEFAULT, no tuning).
  /// `match_li6_radius()` returns the value that puts <r^2> on a target.
  double alpha_d_scale = 1.0;

  /// If non-zero, scale the alpha-d D-wave amplitude by the ROOT s of (G9)
  /// (design sec. 2.7) so that (G5) returns this Q_charge [fm^2].  NOT
  /// sqrt(target/model): Q is LINEAR in the D amplitude through the
  /// interference term (95 % of it) and the deuteron term 2 Q_d D_T is an
  /// offset, so the sqrt rule misses by 7x (it returns -0.048 for a -0.0818
  /// request).  A DEFORMATION DIAL, not a wave function: the root for the
  /// measured Q is s = 0.4032 with P_D(alpha-d) = 3.3e-3, and the writer
  /// stamps it in the header.  REACHABLE RANGE, s in [0, 1]:
  /// [q_matter_analytic_fm2(1)/2, +Q_d] = [-0.615, +0.270] fm^2 for the
  /// default source -- at s = 0 the alpha-d term vanishes and the deuteron's
  /// own +0.270 fm^2 survives.  `validate()` THROWS outside it.  0 = off.
  double quadrupole_target_fm2 = 0.0;

  /// Reserved.  The exact 5-D joint density with the m_S coherences of design
  /// sec. 2.6 restored (a Metropolis problem, not an inverse-CDF one).
  /// Construction THROWS if set; the flag exists so the approximation is
  /// visible in the API rather than buried in a comment.
  bool exact_coherence = false;

  std::size_t n_r = 512;    ///< radial cells of each (x, cos theta) table
  std::size_t n_c = 96;     ///< cos-theta cells (matches TaggedModel's default)
  double r_max_fm = 20.0;   ///< alpha-d and alpha-core grid CEILING
  double rnp_max_fm = 25.0; ///< p-n grid CEILING

  /// Throws std::runtime_error on a bad combination.  Needs the sampler's own
  /// tables for the reachability and norm-loss rules, so the FREE function
  /// checks only what is knowable without them; `ClusterConfigSampler`'s
  /// constructor calls both.
  void validate() const;
};

// ---------------------------------------------------------------- the data

/// One nucleon of one configuration.  Positions are in the ION REST FRAME,
/// c.m. at the origin, in fm.
struct NucleonPos {
  double x = 0.0, y = 0.0, z = 0.0;
  int isospin = +1;    ///< +1 proton, -1 neutron
  int cluster = 0;     ///< 0 = alpha core, 1 = deuteron
};

/// One 6Li configuration.  Nucleons 0-3 are the alpha core, 4 the proton and
/// 5 the neutron of the deuteron.
struct ClusterConfig {
  std::array<NucleonPos, 6> nucleon{};
  int m_ion = 1;       ///< the substate this configuration was drawn for
  double m_s = 1.0;    ///< the deuteron projection drawn (diagnostic)
  double r_ad = 0.0;   ///< |R| [fm]      (diagnostic)
  double r_np = 0.0;   ///< |r_p - r_n|   (diagnostic)
};

/// N configurations plus everything needed to reproduce and label them.
struct ClusterConfigSet {
  std::vector<ClusterConfig> config;
  int m_ion = 1;       ///< -2 marks the unpolarized (equal thirds) set
  double theta_s = 0.0, phi_s = 0.0;
  std::uint64_t seed = 0, run = 0;
  std::string provenance;    ///< every input file + option, one line each
  /// Moments of the SET, computed on the fly.
  double r2_mean_fm2 = 0.0;      ///< <r^2> per nucleon
  double q_matter_fm2 = 0.0;     ///< sum_i <3 z_i^2 - r_i^2> about the axis
  double delta_perp_fm2 = 0.0;   ///< <x^2> - <y^2> per nucleon, lab axes
  /// PER-CONFIGURATION sample standard deviations of the same three moments.
  /// They are what makes every Monte Carlo gate self-calibrating: the gate is
  /// k * sd / sqrt(n), never a hand-typed absolute number.  A single 6Li
  /// configuration has sum_i (3 z_i^2 - r_i^2) scattered over tens of fm^2,
  /// so an absolute gate is meaningless without the variance beside it.
  double r2_sd_fm2 = 0.0, q_matter_sd_fm2 = 0.0, delta_perp_sd_fm2 = 0.0;
  std::array<double, 3> cm{{0.0, 0.0, 0.0}};   ///< max |sum_i r_i| component
};

// ------------------------------------------------- the master density (G1)

/// (G1'), ONE m_S branch: | sum_L C_L(m, m_s) f_L Theta_L^{m-m_s}(c) |^2,
/// with f_0, f_2 the radial functions ALREADY evaluated at x.  `m` and `m_s`
/// are the J = 1 and S = 1 projections.  Terms with |m - m_s| > L are absent
/// (`theta_lm` throws for them, so they are skipped, not evaluated).
double cluster_amp2(double f0, double f2, int m, int m_s, double c);

/// (G1), the m_S-SUMMED density sum_{m_s} `cluster_amp2`.
double cluster_density(double f0, double f2, int m, double c);

/// The (G4)/(G6) functionals of a tabulated radial pair, by trapezoid on the
/// table's OWN abscissa (the table is zero outside it, so no ceiling enters).
/// This is the analytic layer: `radial_moments` on the committed tables is
/// what reproduces every number of the design's sec. 8.
///
/// FREE, and taking the vectors rather than reading them off a sampler, so
/// that a test can feed it sign-flipped copies: Q is invariant under
/// f_L -> -f_L for BOTH L (the unobservable global phase) and `q_int`
/// negates exactly when only f_2 flips.
struct RadialMoments {
  double n0 = 0.0;     ///< integral f_0^2 x^2 dx
  double n2 = 0.0;     ///< integral f_2^2 x^2 dx
  double r2 = 0.0;     ///< integral (f_0^2 + f_2^2) x^4 dx
  double q_int = 0.0;  ///< (1/5) integral 2 sqrt2 f_0 f_2 x^4 dx
  double q_dd = 0.0;   ///< (1/5) integral (-f_2^2) x^4 dx
  double norm() const { return n0 + n2; }
  double p_d() const { return n2 / (n0 + n2); }
  double quadrupole() const { return q_int + q_dd; }   ///< (G4)
};
RadialMoments radial_moments(const std::vector<double>& x,
                             const std::vector<double>& f0,
                             const std::vector<double>& f2);

// ------------------------------------------------ the quadrupole -> a_2 map

/// a_2(m) at |t| [GeV^2] from a point-matter quadrupole, (G7) + (G8):
///
///   delta_{+-1} = q_matter_fm2 / (2 a)   [fm^2, per nucleon, TRANSVERSE axis]
///               -> GeV^-2 through HBARC_GEV_FM,
///   delta_0     = -2 delta_{+-1},
///   a_2(m)      = -(delta_m / 4) |t|.
///
/// `q_matter_fm2` is the m = +-1 state's sum_i (3 z_i^2 - r_i^2); the m = 0
/// value follows from `m` internally, so ONE quadrupole feeds all three
/// states.  (The design sketch's prose reads "with q_matter_analytic_fm2(m)";
/// that would apply the (3m^2-2) factor twice for m = 0.  Its own gate --
/// a2_from_quadrupole(2 * 0.269670, 2, |t|, m) reproducing every
/// `mantysaari_a2_deuteron()` row -- fixes the convention to the one above.)
///
/// A FREE FUNCTION so the deuteron can be fed through it: with
/// (2 * Q_d = 0.53934 fm^2, a = 2) it reproduces `mantysaari_a2_deuteron()`
/// (coherent.hpp) to 8 % at m = +-1 and 8-21 % at m = 0, with zero free
/// parameters.  `ClusterConfigSampler::a2_from_geometry` is this function at
/// (q_matter_analytic_fm2(1), a = 6).
double a2_from_quadrupole(double q_matter_fm2, int a, double t_abs, int m);

/// The EXACT inverse of `a2_from_quadrupole` at m = +-1: the point-matter
/// quadrupole [fm^2] a coefficient a_2(+-1)/|t| [GeV^-2] implies.
///
///     a_2(+-1) = -(q_matter / (2 a)) / (hbar c)^2 * |t| / 4
///   =>  q_matter = -8 a (hbar c)^2 * (a_2 / |t|).
///
/// It exists so that a `CoherentScenario` can be ASKED what it assumes about
/// the target rather than asserting it in prose (open item O4): with
/// `CoherentScenario::a2_m_state(1.0, 1)` as the argument and a = 6 it returns
/// the point-matter quadrupole of 6Li that `eps_b0` implies, which is TWICE
/// the charge one for N = Z.  At the shipped eps_b0 = -0.08, slope_b = 50 that
/// is -1.8691 fm^2 matter = -0.9345 fm^2 charge, i.e. 11.4x the measured
/// `LI6_QUADRUPOLE_FM2` = -0.0818 fm^2 (T23a; sec. C4 of
/// docs/open_items/run_2026-09-03/phase_C_numbers.md).  No new physics number:
/// it is the same (G7)+(G8) map read backwards, with the same HBARC_GEV_FM.
double quadrupole_from_a2_slope(double a2_over_t, int a);

// -------------------------------------------------------------- the sampler

/// Builds its grids ONCE in the constructor and is immutable afterwards --
/// the same discipline as `TaggedModel`, so it is thread-safe with no lock
/// and every accessor is a pure lookup.
class ClusterConfigSampler {
 public:
  explicit ClusterConfigSampler(ClusterConfigOptions opt = {});

  const ClusterConfigOptions& options() const { return opt_; }
  /// One line per input file and per option -- what the sidecar records.
  const std::string& provenance() const { return provenance_; }

  /// One configuration for ion substate `m_ion` in {+1, 0, -1}.  `rng` is the
  /// caller's; the sampler never owns or seeds one.
  ClusterConfig sample(Rng& rng, int m_ion) const;

  /// N configurations, each drawn from its OWN counter-based stream
  /// Rng(seed, run, kConfigStream, i), so config i is bit-identical whatever
  /// the order, the thread count, or N (docs/CONVENTIONS.md).
  ClusterConfigSet sample_set(std::size_t n, int m_ion, std::uint64_t seed,
                              std::uint64_t run = 0) const;
  /// The unpolarized set: equal thirds of m = +1, 0, -1, INTERLEAVED, so the
  /// counts differ by at most one for any N.  `m_ion` is reported as -2.
  ClusterConfigSet sample_set_unpolarized(std::size_t n, std::uint64_t seed,
                                          std::uint64_t run = 0) const;

  // ---- analytic predictions, for the tests and for the header ------------

  double p_d_alpha_d() const;              ///< P_D of the alpha-d relative wave
  /// S_alpha-d: the SOURCE table's own integral (A_0^2 + A_2^2) R^2 dR, which
  /// is 0.854 for `li6.ad` and the two `FitRescaled` waves rescaled to it.
  double s_alpha_d() const { return s_alpha_d_; }
  /// D_T, eq. (G5).  The CLOSED FORM 1 - 0.9 * p_d_alpha_d(); `TaggedModel`'s
  /// same-named quantity is a (k, cos theta) grid quadrature at a different
  /// default P_D, so the two agree at 1e-4, not exactly.
  double tensor_dilution() const;
  double q_matter_analytic_fm2(int m) const;   ///< (G5)
  /// The (G4) split: the interference (2 sqrt2 f_0 f_2) and pure-D (-f_2^2)
  /// parts of Q[R_0, R_2], AFTER the quadrupole dial and the alpha-d scale.
  double q_int_fm2() const;
  double q_dd_fm2() const;
  double r2_analytic_fm2() const;              ///< (G6), per nucleon
  /// (G7): Q_matter(m)/(2A) * sin^2(theta_s) cos(2 phi_s).  Exactly 0 for a
  /// longitudinal axis, which is the correct value, not an error.
  double delta_perp_analytic_fm2(int m) const;
  /// Asymptotic D/S ratio eta = C_2/C_0 of the Coulomb-bound alpha-d channel:
  /// R_2/R_0 divided by W_{-eta_c,5/2}/W_{-eta_c,1/2}(2 kappa R), averaged
  /// over R = `kDsRatioRLoFm` .. `kDsRatioRHiFm`.  kappa and the Sommerfeld
  /// parameter eta_c are DERIVED from `LI6_ALPHA_TAG()` (spectator.hpp), not
  /// retyped: kappa = 0.3074 fm^-1, eta_c = 0.3002.  The naive R_2/R_0 is NOT
  /// eta -- W_2/W_0 is 3.34 / 2.92 / 2.63 at R = 6 / 7 / 8 fm.  -0.0482 for
  /// the default tables (-0.0538 +- 0.0021 for `OverlapRaw`, propagating
  /// li6.ad's own MC errors) against the MEASURED `LI6_ETA_DS_GK`
  /// = -0.025 +- 0.006 +- 0.010: a real but MODERATE D-wave excess, ~2x, not
  /// the 5-15x a naive ratio suggests.
  ///
  /// IT IS NOT AN INDEPENDENT CHECK ON A DIALLED CONFIGURATION.  The (G9)
  /// dial scales R_2 by s and both waves by the same 1/sqrt(n(s)), so
  /// eta(s) = s * eta(1) EXACTLY; running the dial to a target Q and then
  /// reading eta back recovers the dial, not the wave function.  What that
  /// buys is a budget leg anchored on a MEASUREMENT rather than on GFMC:
  /// eta = -0.025 <-> Q_charge = -0.1842 fm^2, so 3.317x of the 7.52x gap is
  /// "too much D wave" and the remaining 2.269x is everything else (T22b).
  /// What it costs is that the leg inherits GK's error bar, which spans
  /// Q_charge from -0.4005 to +0.0298 fm^2 -- through ZERO.  Quote the band.
  double asymptotic_ds_ratio() const;
  /// The `quadrupole_target_fm2` to request so that `asymptotic_ds_ratio()`
  /// lands on `eta_target` -- a CONVERTER onto the one existing dial, not a
  /// second dial (open item C5.1, decided 2026-09-04).
  ///
  /// WHY THERE IS NO eta DIAL.  eta is EXACTLY LINEAR in `quadrupole_dial_s()`
  /// (the dial multiplies R_2 by s and both waves by the same 1/sqrt(n(s)), so
  /// the wave RATIO carries s alone), and `quadrupole_target_fm2` already
  /// bisects that same s.  An `eta_target` option would therefore be a second
  /// name for the SAME one-parameter family -- two knobs onto one physics
  /// number, which docs/CONVENTIONS.md forbids -- and the two could be set to
  /// contradictory values in one options object.  What was missing was not a
  /// dial but the CONVERSION, which lived only as a number in a document
  /// (sec. C1.4's -0.18557 fm^2 at eta = -0.025); this returns it, so the
  /// George-Knutson band re-runs when an input moves instead of rotting.
  ///
  /// It uses the sampler's own PRE-DIAL moments and its current dial, so it is
  /// correct on a dialled sampler too: s = dial_s * eta_target / eta_now.
  /// THROWS if that s leaves [0, 1], i.e. if the requested eta is outside what
  /// scaling this source's D wave can reach.  READ THE BAND: GK's
  /// +-sqrt(stat^2 + syst^2) = +-0.011662 maps onto Q from -0.4005 fm^2 to
  /// +0.0298 fm^2, THROUGH ZERO at 0.86 sigma (T22b), so a Q from a central
  /// eta is a point on a line whose error bar does not pin its sign.
  double quadrupole_for_eta(double eta_target) const;
  /// a_2(m) at |t| for THIS sampler's 6Li geometry.
  ///
  /// IT IS A COEFFICIENT, NOT A SENSITIVITY.  a_2 is LINEAR in |t| and the
  /// coherent sample is exp(-B|t|) with B ~= 39-55 GeV^-2, so the value at
  /// the |t| = 0.3 GeV^2 that gets quoted is 10.6x the modulation an
  /// experiment would actually weight: at the MEASURED Q(6Li) the
  /// information-weighted a_2 is kappa*sqrt(<t^2>) = 0.25 %, not 2.6 %.
  /// Priced in validation/o5_a2_reach.py (open item O5): coherent J/psi over
  /// the WHOLE Q^2 range with BOTH lepton channels at
  /// Scenario::lumi_fb_per_nucleon = 10 fb^-1/u.  QUOTE IT AS A BAND AND
  /// NEVER AS THE POINT -- S = 2.63 sigma at the band's LOW EDGE and
  /// 2.84-3.29 sigma at its TOP, 3 sigma at 8.3-13.0 fb^-1/u, MARGINAL and
  /// inside the {1, 10, 100} fb^-1/u band at both ends
  /// (docs/OPEN_ITEMS_SOLUTIONS.md sec. 11.3b).  The uncorrected point
  /// estimate, 2.62 sigma / 13.1 fb^-1/u, lands 0.3 % under the band's low
  /// edge and may not be quoted alone: it omits a measured beam-energy
  /// factor (UP, x1.12-1.16) and the decay-lepton acceptance x efficiency
  /// (DOWN, unbounded below in this tree), which happen to cancel to 0.7 %.
  /// WHETHER THE BAND'S TOP CROSSES 3 SIGMA IS NOT ESTABLISHED: the top is a
  /// SPAN because the 7Li -> 6Li efficiency substitution straddles 1 when it
  /// is read off entries that share a beam energy (`o5_a2_reach.py`,
  /// `species_scaling_same_energy`).
  /// (An earlier revision said 0.75 sigma / 160 fb^-1/u; that was ONE lepton
  /// channel in ONE Q^2 window, docs/OPEN_ITEMS_SOLUTIONS.md sec. 11.3a.)
  ///
  /// AND READ IT WITH ITS LIMITATION.  `a2_from_quadrupole`, which this calls,
  /// is a CLOSED FORM AND NOT A GOOD-WALKER DIPOLE-MODEL AMPLITUDE: the
  /// target's quadrupole carried through the deuteron's published |t|
  /// dependence, with no amplitude, no saturation and none of their
  /// uncertainties, and with the MATTER quadrupole standing in for the
  /// transverse GLUON anisotropy.  A dipole-model run could move the number
  /// by x1.15 either way across the 3 sigma line.  FOUR further things are
  /// UNESTABLISHED rather than uncertain: no detection efficiency exists
  /// below Q^2 = 0.1 GeV^2 anywhere in this tree, where 85 % of the coherent
  /// rate sits; no decay-lepton reconstruction efficiency exists in this tree
  /// at all, which is what leaves the band OPEN BELOW (only its geometric
  /// half is bounded, at 0.99); the 7Li -> 6Li efficiency substitution
  /// straddles 1, which is what makes the band's top a span; and the
  /// far-forward working point is unchosen -- at LiPolGen's own de-squeezed
  /// 6Li tagging optics the band is 0.73-0.92 sigma with 3 sigma at
  /// 106-167 fb^-1/u, OUTSIDE the {1, 10, 100} band at both ends, and that
  /// last one is the single correction that on its own restores the NO.
  /// docs/open_items/run_2026-09-03/phase_C_numbers.md
  /// sec. C2 (C2.0 for the limitation, C2.8 for every assumption);
  /// docs/OPEN_ITEMS_SOLUTIONS.md sec. 11.3.
  double a2_from_geometry(double t_abs, int m) const;
  /// eps_b0 equivalent: delta_perp_analytic_fm2(+1) [GeV^-2] divided by
  /// `gaussian_slope(sqrt(LI6_R2_POINT_FM2))` = 52.04 GeV^-2.  The MEASURED
  /// point radius is used, not the model's own 2.539 fm (B = 55.2), so the
  /// number is directly comparable with `CoherentScenario::slope_b`'s 50 and
  /// its `eps_b0` band -(0.04 .. 0.13) -- which this returns -0.0506 against.
  double eps_b0_equivalent() const;
  /// The `alpha_d_scale` that puts (G6)'s <r^2> on `target_rms_fm`.  Throws
  /// if the alpha and deuteron terms alone already exceed the target.
  double match_li6_radius(double target_rms_fm = LI6_R_POINT_VMC_FM) const;
  /// {measured `LI6_QUADRUPOLE_FM2`, GFMC AV18+IL7 `LI6_QUADRUPOLE_GFMC_FM2`,
  /// THIS SOURCE's q_matter_analytic_fm2(1)/2} Q_charge [fm^2] -- the
  /// mandatory band.  The third entry moves with `alpha_d_scale` /
  /// `quadrupole_target_fm2` and is never a literal.
  std::array<double, 3> quadrupole_band_fm2() const;
  /// The (G9) root actually applied to the alpha-d D wave; 1.0 when the dial
  /// is off.  "A deformation dial, not a wave function" -- the writer stamps
  /// it beside the band.
  double quadrupole_dial_s() const { return dial_s_; }

  // ---- the same three moments as the sampler's own GRID quadrature -------
  //
  // These are the cell-CDF sums the sampler actually draws from -- no Monte
  // Carlo noise, so they isolate a grid bug from a sampling bug.  They differ
  // from the analytic layer by the trapezoid's own O(h^2) error at the
  // TABLE's 0.1 fm spacing, which does not shrink with `n_r`: for a
  // piecewise-linear f the node trapezoid is HIGH on integral f^2 x^2 by
  // ~(h^2/6) integral (f')^2 x^2 dx, and it inflates the r^4 moment LESS than
  // the norm, so the analytic <R^2> comes out LOW by ~5e-4 relative to the
  // exact integral of the interpolant -- which is what the cell sum converges
  // to.  (Measured on the FitRescaled table: N_trap/N_exact - 1 = +4.4e-4,
  // M_trap/M_exact - 1 = -3.8e-5, so the <R^2> ratio is +4.8e-4; +6.5e-4 on
  // (G6) once the alpha and p-n pieces are in.)
  double r2_grid_fm2() const;                    ///< (G6) from the grids
  double q_matter_grid_fm2(int m) const;         ///< (G5) from the grids
  /// <P2(cos theta_R) | M, m_s> from the grid table.  For M = +1 the alpha-d
  /// orientation is pure L = 2 in two of the three branches, so this is
  /// EXACTLY -2/7 at m_s = -1 (m_L = 2) and +1/7 at m_s = 0 (m_L = 1); an
  /// m_S-marginal draw gives -0.0355 in both and fails by ~15 sigma.
  double p2_alpha_d_grid(int m, double m_s) const;
  /// P(m_s | M) = sum_L N_L |C_L(M, m_s)|^2.
  double p_ms(int m, double m_s) const;

  // ---- the densities, exposed for the tests and for plots ----------------

  /// The alpha core's SOURCE single-nucleon density [fm^-3], normalized so
  /// that 4 pi integral rho r^2 dr = 1, INCLUDING the lambda inflation when
  /// `alpha_cm_inflate` is set.  This is what step 4 draws v_i from; the
  /// recentred s_i follow `rho_alpha_recentred`.
  double rho_alpha(double r_fm) const;
  /// The CLOSED-FORM recentred alpha single-nucleon density [fm^-3], same
  /// normalization.  With s_1 = (3/4) v_1 - (1/4)(v_2+v_3+v_4) and v
  /// isotropic with 3-D transform ftilde(q) = <j_0(q|v|)>, the recentred
  /// density has transform ftilde(3q/4) ftilde(q/4)^3 and <s^2> = (3/4)<v^2>
  /// EXACTLY, for any shape -- which is the whole content of
  /// `alpha_cm_inflate`.  THIS, not `he4.density`, is what the sampled core
  /// histogram must reproduce: recentring does not preserve the shape of a
  /// tabulated density (chi^2/ndf ~ 200 per 1e6 core entries against
  /// he4.density itself).
  ///
  /// INDEPENDENT DRAWS ONLY.  Both statements assume the four v_i are drawn
  /// independently, so this closed form -- and the histogram gate T5(b) built
  /// on it -- hold only at `min_nn_separation_fm` = 0.  With a hard core the
  /// sampled shape is different and `r2_alpha_fm2()` is measured, not this.
  double rho_alpha_recentred(double s_fm) const;
  /// (G1') -- the alpha-d density CONDITIONED on the deuteron projection m_s,
  /// |A_{m_s}(m; R, c)|^2.  This, not the m_s-summed one, is what the sampler
  /// draws from; the signature mirrors `TaggedModel::sample_kc(m_ion, m_s)`
  /// for the same reason.
  double rho_alpha_d(double R_fm, double c, int m, double m_s) const;
  /// Convenience: the m_s-SUMMED density (G1).  For plots and marginal checks
  /// ONLY -- sampling from it decorrelates R-hat from m_s.
  double rho_alpha_d_summed(double R_fm, double c, int m) const;
  /// (G1)/(G2)/(G3) for the p-n pair in deuteron projection m_s.
  double rho_np(double r_fm, double c, double m_s) const;

  // ---- the tables themselves, for the sign and normalization tests -------

  /// The alpha-d abscissa [fm] and its NORMALIZED radial functions R_0, R_2
  /// (integral (R_0^2 + R_2^2) R^2 dR = 1), after the quadrupole dial.  The
  /// global phase is fixed by convention to R_0(R -> infinity) > 0.
  const std::vector<double>& alpha_d_grid() const { return ad_x_; }
  const std::vector<double>& alpha_d_wave(int l) const;
  /// The deuteron abscissa [fm] and f_0 = u/r, f_2 = w/r, normalized the same
  /// way (so `radial_moments(np_grid(), np_wave(0), np_wave(2))` returns the
  /// file's own header values: norm 1, dstate, qm = quadrupole()/4, and
  /// rd = sqrt(r2)/2).
  const std::vector<double>& np_grid() const { return np_x_; }
  const std::vector<double>& np_wave(int l) const;

  /// The three (G6) pieces [fm^2], for the record and for `match_li6_radius`.
  /// <s^2> of the recentred core: the closed form (3/4)<v^2> at
  /// `min_nn_separation_fm` = 0, the constructor's measured value with a hard
  /// core (see that option).
  double r2_alpha_fm2() const { return s2_alpha_; }
  double r2_np_fm2() const { return r2_np_; }           ///< <r_np^2>
  double r2_alpha_d_fm2() const;                        ///< <R^2>, scaled

 private:
  struct Grid {
    std::vector<double> x;      ///< cell CENTRES
    double dx = 0.0;
  };
  std::size_t ms_index(double m_s) const;
  const std::vector<double>& ad_cdf(int m, double m_s) const;
  const std::vector<double>& np_cdf(double m_s) const;
  void build_tables();
  void draw_alpha(Rng& rng, std::array<std::array<double, 3>, 4>& s) const;

  ClusterConfigOptions opt_;
  std::string provenance_;

  // radial tables (analytic layer): abscissa + normalized waves
  std::vector<double> ad_x_, ad_f0_, ad_f2_;
  std::vector<double> np_x_, np_f0_, np_f2_;
  std::vector<double> core_x_, core_rho_;      ///< SOURCE density (inflated)
  std::vector<double> core_s_, core_rho_s_;    ///< recentred closed form

  RadialMoments ad_m_{}, np_m_{};
  /// The alpha-d moments BEFORE the (G9) quadrupole dial (already normalized
  /// to norm 1, already sign-fixed).  `quadrupole_for_eta` needs them: the
  /// dial's own q_charge(s) closed form is built from the undialled table.
  RadialMoments ad_m_base_{};
  double s_alpha_d_ = 0.0;      ///< the source table's own norm
  double dial_s_ = 1.0;
  double s2_alpha_ = 0.0;       ///< <s^2> after recentring
  double r2_np_ = 0.0;
  double q_d_fm2_ = 0.0;        ///< the deuteron's own quadrupole moment

  // sampling grids
  Grid gr_ad_, gr_np_, gr_c_, gr_core_;
  std::vector<double> core_cdf_;                       ///< 4 pi rho v^2 CDF
  std::vector<std::vector<double>> ad_cdf_;            ///< [(m+1)*3 + i_ms]
  std::vector<std::vector<double>> np_cdf_;            ///< [i_ms]
  std::vector<double> ad_weight_;                      ///< P(m_s | M), 9 slots
};

// ------------------------------------------------------------- the writer

/// Which flavour of the subnucleondiffraction configuration table to write.
enum class SndConfigFormat {
  /// BYTE-COMPATIBLE with `he3.dat`: numeric rows only, NO comments.  The
  /// upstream reader (`src/nucleons.cpp`) does a bare `ss >> x` per field and
  /// would silently read garbage from a '#' line.  Metadata goes to the
  /// sidecar.  DEFAULT.
  He3Compatible = 0,
  /// The same rows with a '#'-prefixed header block.  Requires the patched
  /// reader of the design's sec. 5.3.  THE IN-FILE HEADER IS THE PRE-FILL
  /// COPY: C++ writes it, so its `md5` and `git` fields stay null even after
  /// `python/lipolgen/configs.py` has filled them in the separate
  /// `.meta.json`.  The SIDECAR is authoritative; the in-file block is a
  /// convenience for eyeballing a table, not a second record.
  Annotated = 1,
};

/// Writes `set` as `path` plus the sidecar `path + ".meta.json"`, which
/// always carries m, (theta_s, phi_s), seed/run, EVERY option field, every
/// input file with its size and its own printed normalization,
/// LIPOLGEN_VERSION, and the analytic and sampled moments.  Returns the
/// number of rows written.
///
/// Row layout (fm, ion rest frame, c.m. at the origin, quantization axis
/// already applied), `%.17g`, space separated, LF, no header, no blank line:
///
///     x1 y1 z1 ... x6 y6 z6  t1 ... t6  m
///
/// 18 coordinates, then six isospins (+1 p, -1 n; ordering alpha, alpha,
/// alpha, alpha, p, n), then the ion substate.  Fields 19-25 are past the 3A
/// the consumer reads and are invisible to it, exactly as `he3.dat`'s own
/// trailing fields are.
///
/// NO md5 AND NO GIT SHA FROM C++.  There is no hash function anywhere in
/// src/core or include/lipolgen, no JSON emitter in the core, and CMake
/// compiles in LIPOLGEN_VERSION only.  The writer emits a flat, hand-rolled
/// JSON object with `"md5": null` and `"git": null`; `python/lipolgen/
/// configs.py` fills them from hashlib and `git rev-parse HEAD`.
std::size_t write_snd_configs(const ClusterConfigSet& set,
                              const ClusterConfigSampler& sampler,
                              const std::string& path,
                              SndConfigFormat fmt = SndConfigFormat::He3Compatible);

}  // namespace lipolgen

#endif  // LIPOLGEN_CLUSTER_CONFIG_HPP
