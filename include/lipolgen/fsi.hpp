#ifndef LIPOLGEN_FSI_HPP
#define LIPOLGEN_FSI_HPP

/// \file fsi.hpp
/// Final-state interaction of the DIS debris X with the tagged SPECTATOR
/// cluster, as a per-event WEIGHT.
///
/// ---------------------------------------------------------------------------
/// FSI IS A WEIGHT, NEVER A SHIFT.
///
/// The spectator four-vector is the measurement.  Moving it would break both
/// the "never recoil-correct the light spectator" rule of `tagged.hpp` and the
/// exact whole-nucleus conservation the whole event record rests on, so the
/// rescattering enters as a multiplicative weight on the event and on nothing
/// else.  A null `FsiWeight` is today's plane-wave impulse approximation,
/// bit for bit.
///
/// ---------------------------------------------------------------------------
/// THE MODEL (Cosyn-Weiss / Glauber, transplanted to X-cluster).
///
/// At EIC energies the rescattering is EIKONAL: it transfers TRANSVERSE
/// momentum only, so the spectator light-cone fraction alpha_s -- equivalently
/// k_z in the ion rest frame -- is a spectator of the FSI too.  The distortion
/// is therefore a 2-D convolution at FIXED k_z,
///
///   Psi_FSI(k_z, k_T) = psi(k_z, k_T)
///                       - int d^2k'_T/(2 pi)^2 Gtil(k_T - k'_T) psi(k_z, k'_T)
///
/// which is the momentum-space form of the Glauber statement
/// Psi_FSI(b) = psi(b) [1 - Gamma(b)], and
///
///   Gtil(q) = int d^2b e^{i q.b} Gamma(b),   Gtil(0) = sigma_Xa (1 - i eps)/2 .
///
/// The longitudinal integral over the rescattering propagator has been done by
/// residue (pole part); the principal-value (off-shell) piece is O(Delta/kappa)
/// with Delta = k^2/(2 m_spec) = 0.3 / 1.3 / 5.4 MeV at k = 0.05 / 0.1 / 0.2 GeV
/// against kappa(alpha-d) = 60.7 MeV, and is DROPPED.  That is the one
/// approximation the model stands on, and it is also why the weight comes out
/// symmetric under theta_k -> pi - theta_k (see `tests/test_fsi.cpp`).
/// Independent support: Ciofi degli Atti & Kaptari (nucl-th/0407024) find GEA
/// and plain Glauber differ by < 3-4 % for p_m < 0.6 GeV/c.
///
/// ---------------------------------------------------------------------------
/// HOW THE ANGULAR STRUCTURE IS HANDLED -- read this before quoting a number.
///
/// The distortion is applied AT THE AMPLITUDE LEVEL, wave by wave and
/// azimuthal component by azimuthal component, and the weight is then the
/// ratio of the m-SUMMED (unpolarized) distorted density to the impulse one:
///
///   A_{L m}(k) = psihat_L(k) Y_L^m(khat) = g_{Lm}(k_z, k_T) e^{i m phi}
///   g^FSI_{Lm}(k_z, k_T) = g_{Lm} - int_0^inf b db J_m(k_T b) Gamma(b) G_m(b)
///   G_m(b) = int_0^inf k'_T dk'_T J_m(k'_T b) g_{Lm}(k_z, k'_T)
///
///   w(k, cos theta_k) = sum_{L,m} (P_L/(2L+1)) |g^FSI_{Lm}|^2
///                     / sum_{L,m} (P_L/(2L+1)) |g_{Lm}|^2
///
/// The azimuthal harmonic is preserved by the convolution because Gamma is
/// azimuthally symmetric, so each (L, m) needs one order-|m| Hankel pair and
/// nothing is approximated in the angular dependence of a single partial wave.
/// What IS an approximation, deliberately:
///
///   * the numerator and the denominator are both summed over m and averaged
///     over the ion projection M, i.e. the weight is the UNPOLARIZED shape
///     distortion.  S-D interference cancels in that sum (the unpolarized
///     density Sum_L P_L psihat_L^2 / 4pi is isotropic), so no interference
///     term survives in either the IA or the FSI density.
///   * consequently the weight is SPIN INDEPENDENT by construction, which is
///     exactly what the literature supports: nothing constrains the spin
///     dependence of the rescattering (Cosyn-Weiss arXiv:2603.23700 VI C lists
///     it as an open question).  Quote it as an unpolarized-shape systematic,
///     NEVER as a correction to A_zz.
///
/// psi is the channel's OWN radial: Hulthen or the ANL VMC tables, S+D for the
/// 6Li alpha tag and P for the 7Li one, taken through `Wave::radial` with the
/// same sqrt(P_L)/norm normalization `TaggedModel` applies.
///
/// ---------------------------------------------------------------------------
/// THE PROFILE, AND WHY sigma_Xa IS NOT 4 sigma_XN.
///
/// `FsiVariant::GlauberCluster` (the default, item (a)) builds the X-cluster
/// profile by Glauber from the X-N one over the cluster's own Gaussian
/// point-nucleon density T_a,
///
///   Gamma_a(b) = 1 - [1 - (Gamma_N conv T_a)(b)]^A_spec
///   Gamma_N(b) = (sigma_XN (1 - i eps) / (4 pi B_XN)) exp(-b^2 / 2 B_XN)
///
/// -- both Gaussians, so the convolution is analytic and the widths add.  At
/// sigma_XN = 40 mb, eps = -0.5, B_XN = 6 GeV^-2 and a^2 = 0.700 fm^2 this
/// gives sigma_Xalpha = 131.0 mb (NOT 4 x 40 = 160: the shadowing factor is
/// 0.819) and an effective slope B_alpha = 27.2 GeV^-2.
///
/// The sigma_tot / sigma_el split asked for by the literature review is
/// AUTOMATIC in this form and is not a second knob: expanding
/// |psi - C|^2 = |psi|^2 - 2 Re(psi* C) + |C|^2, the LINEAR (absorptive) term
/// carries 2 Re int d^2b Gamma_a = sigma_tot(Xa) and the QUADRATIC (gain,
/// refraction) term carries int d^2b |Gamma_a|^2 = sigma_el(Xa), because
/// Gamma_a is the coherent ELASTIC profile of the intact cluster.  Measured on
/// this implementation at sigma_XN = 40 mb: sigma_tot = 131.0 mb,
/// sigma_el = 35.2 mb, ratio 0.269 -- against the MEASURED alpha-p
/// 31.4 / 121.5 = 0.258 (Blinov et al., nucl-ex/9910012).  So ~73 % of
/// rescatterings destroy the alpha tag, the distortion is nearly purely
/// absorptive, and int w dGamma is a SURVIVAL PROBABILITY rather than 1.  That
/// deviation is physical and is logged (`survival()`).
/// `GlauberFsiOptions::elastic_gain` rescales the quadratic term alone, for a
/// run that wants to impose a measured sigma_el/sigma_tot instead.
///
/// `FsiVariant::GlauberNucleon` (item (b)) is the PER-NUCLEON variant.  Note
/// what is and is not different: for an UNCORRELATED (independent-particle)
/// cluster density the Ciofi degli Atti-Kaptari per-nucleon product
/// <Prod_i [1 - Gamma_N(b - s_i)]> factorizes into Prod_i <1 - Gamma_N(b - s_i)>
/// = [1 - (Gamma_N conv T_a)(b)]^A -- i.e. it is ALGEBRAICALLY IDENTICAL to
/// the coherent-cluster form above, and a separate code path reproducing it
/// would print the same numbers.  What this variant therefore implements is
/// the genuinely different SINGLE-SCATTERING (optical, unshadowed) limit
///
///   Gamma_a(b) = A_spec (Gamma_N conv T_a)(b)      =>   sigma_Xa = A sigma_XN
///
/// which is the 4 sigma_XN = 160 mb upper bound on the absorption.  Read the
/// bracket PRECISELY: it is a bracket on the POINTWISE low-k_T absorption
/// (measured at sigma_XN = 40 mb: w = 0.309 against the cluster variant's
/// 0.381 at k = 0.1, theta = 90 deg), NOT on the integrated survival.  The
/// unshadowed profile also scales the QUADRATIC gain term by A^2, second
/// order in Gamma_N and formally beyond single scattering, which feeds enough
/// strength back into the tag that the integrated survival comes out ABOVE
/// the shadowed default (0.582 against 0.517 on the pure-S 6Li channel) --
/// pinned in `tests/test_fsi.cpp` so a change is loud.
///
/// THE TWO VARIANTS ARE NOT ONE, AND THE DIFFERENCE IS MEASURED (2026-09-04,
/// D3).  `PLAN.md` carried "the per-nucleon Glauber FSI variant is
/// algebraically identical to the cluster one -- the two 'variants' are one"
/// as an open item; it is RETRACTED, because the identity in the paragraph
/// above is a statement about a per-nucleon product that this file does not
/// implement, not about the single-scattering limit that it does.  On one
/// event stream reweighted three ways (`--channel tagged-6Li-alpha --events
/// 20000 --seed 1234 --fsi {off,glauber-cluster,glauber-nucleon}`, plan
/// `tensor-thirds` at P_z = 0.7 / P_zz = 0.6 / P_e = 0.7, at the DEFAULT
/// sigma_XN = 40 mb -- the whole claim is one stream at one end of the
/// 20-40 mb band; the kinematic columns are bit-identical across the three):
///
///     Sum w             20000.0      10419.07     11632.10
///     Sum w / Sum w_off      1       0.520954     0.581605
///     per-event min..max     1    0.0208..1.404  0.2139..6.952
///
/// and per event w_nucleon/w_cluster has percentiles [1, 5, 25, 50, 75, 95,
/// 99] = 0.821, 0.827, 0.853, 0.897, 0.934, 3.337, 5.935 with a maximum of
/// 68.52.  **99.50 % of the 20 000 events differ by more than 1 %**, and
/// `np.allclose(rtol = 1e-14)` is False.  Different by construction, not by
/// parameter choice.  Pinned in `tests/test_fsi.cpp` ("the two variants
/// differ on essentially every event").
///
/// TODO (open): the two corrections that would make a true per-nucleon
/// product differ from (a) -- the centre-of-mass constraint Sum_i s_i = 0
/// and short-range NN correlations in the cluster density -- are not
/// implemented, and ONE OF THEM CANNOT BE BUILT FROM COMMITTED DATA:
///
///   * the c.m. constraint is a change of T_a alone, i.e. of
///     `cluster_point_a2_fm2` (Gartenhaus-Schwartz on a Gaussian gives
///     a^2_int = a^2 (1 - 1/A) = 0.5250 fm^2 for the alpha).  A better input
///     to the SAME variant, not a different variant;
///   * short-range NN correlations need the TWO-body density
///     rho_2(s_i, s_j) -- that is the whole content of <Prod_i(.)> !=
///     Prod_i<.> -- and the tree has none.  `data/vmc/density/*.density` are
///     ONE-body point-proton densities, which is exactly the T_a the cluster
///     form already convolves, and the ANL page they are fetched from
///     publishes one-body densities only.  So the SRC correction needs an
///     input the repository does not have; it is not a coding task.
///
/// And note what the input is NOT: the rescattering here is the DIS debris X
/// off a SPECTATOR NUCLEON, so the amplitude is X-N and the cross section is
/// `sigma_xn_mb`.  A sigma_NN would enter only for the spectator cluster's
/// own internal absorption, which is not what this weight is; substituting
/// one would change the physics, not the variant.
///
/// ---------------------------------------------------------------------------
/// THE OPEN PHYSICS INPUT.  sigma_XN(W) is the weakest number in the model.
/// At EIC energies the DIS debris has a formation length of many fm, so the
/// effective cross section is well below the free-hadron 40 mb and the 20 mb
/// row is arguably the realistic one.  The default is 40 mb (free hadron) and
/// the documented BAND is 20-40 mb; `GlauberFsiOptions::formation_ramp` turns
/// on a linear-in-W interpolation between two anchors, and grids are then
/// built on a small sigma ladder and interpolated per event so throughput is
/// untouched.  Band it; never quote one row alone.

#include <complex>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "lipolgen/constants.hpp"
#include "lipolgen/event.hpp"
#include "lipolgen/tagged.hpp"

namespace lipolgen {

// ----------------------------------------------------------------- units

/// mb per GeV^-2, derived from `GEV2_TO_PB` so the conversion has ONE
/// definition in the library (docs/CONVENTIONS.md).  0.3894.
inline constexpr double GEV2_TO_MB = GEV2_TO_PB * 1e-9;
/// GeV^-2 per mb.  2.5680.
inline constexpr double MB_TO_GEV2 = 1.0 / GEV2_TO_MB;
/// hbar c [GeV fm].  Alias of the library-wide single definition
/// `HBARC_GEV_FM` (constants.hpp, docs/CONVENTIONS.md).  The prototype
/// (fsi_alpha.py) carried CODATA 0.1973269804; the 3e-5 relative difference
/// moves every profile integral by ~2e-4 absolute, far inside the 2e-3
/// tolerances `tests/test_fsi.cpp` pins the prototype table at.
inline constexpr double HBAR_C_GEV_FM = HBARC_GEV_FM;
/// fm^2 per mb.
inline constexpr double MB_TO_FM2 = 0.1;

/// a^2 [fm^2] of the Gaussian POINT-NUCLEON density of a spectator cluster,
/// rho ~ exp(-r^2/2a^2)  =>  T(s) = (1/2 pi a^2) exp(-s^2/2a^2), from the
/// measured charge radius with the proton's folded out:
///
///     a^2 = (r_ch^2(cluster) - r_ch^2(p)) / 3 .
///
/// alpha: (1.6755^2 - 0.8409^2)/3 = 0.7001 fm^2.  A single nucleon returns 0
/// (it IS the point scatterer), and an unknown cluster throws.
double cluster_point_a2_fm2(int z, int a);

// ------------------------------------------------------------ the interface

/// Which X-cluster profile the weight is built from.
enum class FsiVariant : int {
  /// (a) coherent-cluster amplitude, Glauber-shadowed over the cluster's own
  /// Gaussian profile.  THE DEFAULT.
  GlauberCluster = 0,
  /// (b) per-nucleon, single-scattering (optical, unshadowed) limit --
  /// sigma_Xa = A sigma_XN.  The systematic bracket above (a); see the file
  /// header for why the full Ciofi-Kaptari product coincides with (a).
  GlauberNucleon = 1,
};

/// Everything the weight is allowed to depend on.
///
/// (k, cos_theta_k, phi_k) are the spectator's spherical components in the
/// ION REST FRAME (the spin frame of `TaggedEvent`); the FSI kernel is
/// azimuthally symmetric, so `phi_k` is carried for completeness and is not
/// used by `GlauberFsiWeight`.  (w, q2, x) drive the sigma_XN(W) ramp.
struct FsiKinematics {
  double k = 0.0;              ///< |k| [GeV]
  double cos_theta_k = 0.0;    ///< cos of the polar angle from +z (the photon)
  double phi_k = 0.0;          ///< azimuth [rad]; unused by the eikonal kernel
  double w = 0.0;              ///< hadronic W [GeV] of the DIS subprocess
  double q2 = 0.0;             ///< Q^2 [GeV^2]
  double x = 0.0;              ///< x_Bj
  int spectator_z = 2;         ///< charge of the tagged spectator cluster
  int spectator_a = 4;         ///< mass number of the tagged spectator cluster
  /// The channel the event was drawn from.  `Channel::Inclusive` means
  /// "unspecified" and disables the channel check.
  Channel channel = Channel::Inclusive;
};

/// FSI as a per-event weight.  Implementations must be IMMUTABLE after
/// construction: `Pipeline::for_each(sink, nthreads)` calls them concurrently.
class FsiWeight {
 public:
  virtual ~FsiWeight() = default;

  /// FSI / IA ratio.  >= 0.  1 means no distortion.
  virtual double weight(const FsiKinematics& kin) const = 0;
  /// Norm-preserving variant: `weight()` divided by its RATE-weighted mean,
  /// for a run that wants the SHAPE distortion without changing the total
  /// rate.  Averages to 1 over the model's own spectator distribution.
  virtual double weight_normalised(const FsiKinematics& kin) const = 0;
  /// The X-N cross section [mb] the model used at this W -- what a run
  /// prints so the band it was taken at is on the record.
  virtual double sigma_eff_mb(double w) const = 0;
};

// ------------------------------------------------------ the implementation

struct GlauberFsiOptions {
  /// (a) coherent cluster, or (b) per-nucleon single scattering.
  FsiVariant variant = FsiVariant::GlauberCluster;

  // --- the X-N amplitude -------------------------------------------------
  /// sigma_XN [mb].  40 = free hadron; the documented BAND is 20-40 and the
  /// 20 mb end is arguably the realistic one at EIC formation lengths.
  double sigma_xn_mb = 40.0;
  /// Re/Im of the X-N amplitude (Cosyn-Sargsian Deeps fit).
  double eps = -0.5;
  /// X-N diffractive slope [GeV^-2] (b0 = 0.5 fm => B = 6.4).
  double b_xn = 6.0;
  /// Rescaling of the QUADRATIC (gain / refraction) term alone, i.e. of the
  /// effective sigma_el/sigma_tot.  1 = Glauber's own, which already comes
  /// out at 0.269 against the measured alpha-p 0.258 (see the file header).
  double elastic_gain = 1.0;

  // --- the sigma_XN(W) formation ramp ------------------------------------
  /// OFF by default, so a default run is the reproducible 40 mb row.  On, the
  /// model interpolates sigma_XN linearly in W between the two anchors and
  /// clamps outside them; the grid is built on an `n_sigma_grid`-point ladder
  /// and interpolated per event, so the throughput cost is one extra lerp.
  bool formation_ramp = false;
  double ramp_w_lo = 2.0;            ///< [GeV]
  double ramp_sigma_lo_mb = 40.0;    ///< sigma_XN at (and below) `ramp_w_lo`
  double ramp_w_hi = 10.0;           ///< [GeV]
  double ramp_sigma_hi_mb = 20.0;    ///< sigma_XN at (and above) `ramp_w_hi`
  std::size_t n_sigma_grid = 5;      ///< ladder size when the ramp is on

  // --- numerics ----------------------------------------------------------
  /// Tabulation box: the weight is evaluated on a regular (k_z, k_T) grid on
  /// [0, k_max]^2 and read by BILINEAR interpolation, exactly the shape of
  /// `acceptance_weights`.  (k_z, k_T) rather than (k, cos theta_k) because
  /// the eikonal kernel is a function of k_T alone, so the table is smooth in
  /// these variables and even in k_z; an event at (k, c) reads
  /// (|k c|, k sqrt(1-c^2)).
  double k_max = 1.2;                ///< matches `TaggedModel`'s default
  std::size_t n_kz = 161;
  std::size_t n_kt = 161;
  /// Hankel quadrature: the k'_T integration grid and the impact-parameter
  /// grid (b in GeV^-1; 50 GeV^-1 = 9.87 fm, far past the profile's reach).
  double kt_int_max = 4.0;
  std::size_t n_kt_int = 400;
  double b_max = 50.0;
  std::size_t n_b = 512;
  /// Normalization grid of the radial waves -- `TaggedModel`'s own
  /// linspace(1e-4, k_norm_max, n_norm), so P_L means the same thing here.
  double k_norm_max = 1.2;
  std::size_t n_norm = 280;
  /// Ceiling on the returned weight.  The FSI/IA ratio is UNBOUNDED wherever
  /// the impulse amplitude has a node or a steeply falling tail that the
  /// rescattering feed-in overwhelms (the ratio exceeds 1 past
  /// k_T ~ 0.25 GeV, where it is a sign flip of the amplitude, not a bug).
  /// The rate there is negligible, but a Monte-Carlo weight must not be, so
  /// it is clipped and the clipping is reported by `clipped_grid_fraction()`.
  double w_max = 50.0;
};

/// Cosyn-Weiss / Glauber FSI on the tagged spectator cluster.
///
/// Immutable after construction and safe to share between threads.  The
/// constructor does all the work: one order-|m| Hankel pair per (L, |m|, k_z)
/// per sigma ladder point, ~0.2 s at the default grid.
class GlauberFsiWeight : public FsiWeight {
 public:
  /// From a channel.  The radial waves and kappa are the channel's own; the
  /// spectator (Z, A) fixes the cluster profile.
  explicit GlauberFsiWeight(TaggedChannel channel,
                            GlauberFsiOptions opt = GlauberFsiOptions());
  /// From a built model -- same thing, with the normalization grid taken from
  /// the model so P_L means bit-for-bit what it means there.
  explicit GlauberFsiWeight(const TaggedModel& model,
                            GlauberFsiOptions opt = GlauberFsiOptions());

  double weight(const FsiKinematics& kin) const override;
  double weight_normalised(const FsiKinematics& kin) const override;
  double sigma_eff_mb(double w) const override;

  const GlauberFsiOptions& options() const { return opt_; }
  const TaggedChannel& channel() const { return channel_; }

  // --- the profile, exposed so a run can print what it used ---------------

  /// sigma_tot(X-cluster) [mb] = 2 Re int d^2b Gamma_a(b), at this sigma_XN.
  double sigma_cluster_mb(double sigma_xn_mb) const;
  /// sigma_el(X-cluster) [mb] = int d^2b |Gamma_a(b)|^2 -- the coherent
  /// ELASTIC cross section, i.e. what the quadratic (gain) term carries.
  double sigma_cluster_el_mb(double sigma_xn_mb) const;
  /// Effective slope B_a [GeV^-2] from -2 dln|Gtil_a|/dq^2 at q -> 0.
  double slope_cluster_gev2(double sigma_xn_mb) const;
  /// Gtil_a(q) [GeV^-2], the momentum-space profile itself.
  std::complex<double> gtilde(double q_gev, double sigma_xn_mb) const;
  /// Gamma_a(b), b in fm.
  std::complex<double> gamma_profile(double b_fm, double sigma_xn_mb) const;

  // --- diagnostics --------------------------------------------------------

  /// int w dGamma / int dGamma over the channel's own spectator distribution
  /// -- the tagged-cluster SURVIVAL PROBABILITY, which is what
  /// `weight_normalised` divides by.  Below 1 because the rescattering is
  /// nearly purely absorptive; log it.
  double survival(double w = 0.0) const;
  /// Fraction of grid cells that hit `GlauberFsiOptions::w_max`.
  double clipped_grid_fraction() const { return clipped_; }
  /// The sigma ladder [mb] the grids were built on (one entry when the
  /// formation ramp is off).
  const std::vector<double>& sigma_ladder() const { return sigma_ladder_; }
  /// The raw weight table of ladder point `i`, row-major (n_kz, n_kt).
  const std::vector<double>& grid(std::size_t i = 0) const;

 private:
  struct Layer {
    std::vector<double> w;      ///< (n_kz * n_kt), row-major
    double sigma_mb = 0.0;
    double mean_w = 1.0;        ///< the rate-weighted mean = survival
  };

  void build();
  void build_layer(double sigma_mb, Layer& out);
  /// Bilinear read of one layer at (k_z, k_T), both already clamped.
  double read(const Layer& lay, double kz, double kt) const;
  /// (layer index, blend) for an event at this W.
  void ladder_at(double w, std::size_t& i0, std::size_t& i1, double& f) const;

  TaggedChannel channel_;
  GlauberFsiOptions opt_;
  int spec_z_ = 2, spec_a_ = 4;
  double a2_fm2_ = 0.0;              ///< cluster point-nucleon Gaussian a^2
  Channel event_channel_ = Channel::Inclusive;
  /// (L, |m|, weight P_L/(2L+1) x multiplicity, normalized radial scale).
  struct Component { int l = 0; int m = 0; double wgt = 1.0; double scale = 1.0; };
  std::vector<Component> comp_;
  std::vector<double> sigma_ladder_;
  std::vector<Layer> layers_;
  double dkz_ = 0.0, dkt_ = 0.0;
  double clipped_ = 0.0;
};

/// `PipelineConfig::fsi`: which FSI model a run applies to its tagged
/// spectator, if any.
enum class PipelineFsi : int {
  Off = 0,             ///< today's PWIA, bit for bit
  GlauberCluster = 1,  ///< `FsiVariant::GlauberCluster`
  GlauberNucleon = 2,  ///< `FsiVariant::GlauberNucleon`
};

/// Name of a `PipelineFsi` ("off", "glauber-cluster", "glauber-nucleon").
const char* pipeline_fsi_name(PipelineFsi f);

}  // namespace lipolgen

#endif  // LIPOLGEN_FSI_HPP
