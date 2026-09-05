// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LIPOLGEN_TRITON_SF_HPP
#define LIPOLGEN_TRITON_SF_HPP

/// \file triton_sf.hpp
/// The triton's nucleon SPECTRAL FUNCTION S_N(k, E) -- the joint density of
/// the struck nucleon's momentum k and the residual two-nucleon system's
/// excitation -- replacing the sequential two-body decay that
/// `ClusterBreakup`'s Triton branch does today (`breakup.hpp`, "CRUDE AND
/// FLAGGED").
///
/// ---------------------------------------------------------------------------
/// THE THREE CHANNELS.  A triton is p + n + n, so removing one nucleon leaves
/// a two-nucleon system that either IS bound or is not, and which of the two
/// depends on the species removed:
///
///   struck n -> remnant d, BOUND, E = 0                 "2-body"  (n_0)
///   struck n -> remnant (p n) continuum, E > 0          "3-body"  (n_1)
///   struck p -> remnant (n n) continuum, ALWAYS         "3-body"  (no bound
///                                                        nn state exists)
///
/// The branching between the first two is NOT a constant.  Ciofi degli Atti
/// and Simula (PRC 53 (1996) 1689) split the momentum distribution as
/// n(k) = n_0(k) + n_1(k), where n_0 is the part in which the residual (A-1)
/// system is left in its GROUND state and n_1 the correlated part in which it
/// is not, and their Eq. (28) makes the split a pair of spectroscopic factors,
///
///     S_0 = int_0^inf dk k^2 n_0(k),   S_1 = int_0^inf dk k^2 n_1(k),
///     S_0 + S_1 = 1
///
/// (note: the CS convention has NO 4 pi in that integral).  So the k-dependent
/// branching is a RATIO of two things the parameterization already contains,
///
///     p_2(k) = n_0(k) / (n_0(k) + n_1(k)),
///
/// and its k-integral is S_0 itself, not an assumed constant.
///
/// ---------------------------------------------------------------------------
/// PROVENANCE OF EVERY COEFFICIENT.  All of them are transcribed, no refit.
///
///   n_0(k) = sum_i A_i exp(-B_i k^2) / (1 + C_i k^2)^2       CS Eq. (74)
///            k in fm^-1, n in fm^3; Table A.1 for A = 2, 3, 4.
///
///     A = 2 (2H)          157.4 / 1.24 / 18.3
///                         0.234 / 1.27 / -
///                         0.00623 / 0.220 / -
///     A = 3 (3He, PROTON) 31.7 / 1.32 / 5.98
///                         0.00266 / 0.365 / -
///     A = 4 (4He)         4.33 / 1.54 / 0.419
///                         5.49 / 4.90 / -
///
///   These are the numbers BeAGLE's `DT_KFERMI`
///   (dpmjet3.0-5F-new.f:17216-17403) carries, checked column for column
///   against the paper's own Table A.1 -- including its own in-code comment
///   "These are n0k parametrization not including n1k".
///
///   n_1(k) for A = 3 is CS Eq. (76), the one nucleus for which they do NOT
///   use the two-Gaussian Eq. (75):
///
///     n_1(k) = 7.40 exp(-1.23 k^2) / (1 + 3.21 k^2)^2 + 0.0139 exp(-0.234 k^2)
///
///   which is the SAME functional family as n_0, so one `CsTerm` list carries
///   both.  For A = 4 the two-Gaussian Eq. (75) with Table A.3,
///   0.665 exp(-2.15 k^2) + 0.0244 exp(-0.22 k^2), is carried for the
///   cross-check.  A = 2 has no n_1 at all: the deuteron has no excitable
///   remnant.
///
///   MEASURED HERE from those coefficients alone (`tests/test_triton_sf.cpp`):
///
///     A = 2:  S_0 = 1.00310                          (1, as it must be)
///     A = 3:  S_0 = 0.652548,  S_1 = 0.347113,  S_0 + S_1 = 0.99966
///     A = 4:  S_0 = 0.799576,  S_1 = 0.198250,  S_0 + S_1 = 0.99783
///
///   The sums closing on 1 to 3.4e-4 is the check that the n_0 and n_1 sets
///   belong to each other: CS Eq. (28) says S_0 + S_1 = 1 exactly, and nothing
///   in the transcription was tuned to make it come out.  S_0 = 0.653 for
///   A = 3 is the 3He(e,e'p)d spectroscopic factor, the ~2/3 the literature
///   quotes; the paper's own text says "S_0 is equal to ~0.65 and ~0.8
///   (implying S_1 ~ 0.35 and ~0.2) for 3He and 4He".
///
/// ---------------------------------------------------------------------------
/// 3He OR 3H?  Table A.1's A = 3 column is headed "3He (proton)" -- it is the
/// PROTON momentum distribution in 3He, i.e. the MAJORITY species (3He is
/// p p n).  Under the isospin mirror 3He <-> 3H, p <-> n, that is exactly the
/// NEUTRON momentum distribution in the triton, and the bound remnant of the
/// mirror -- (p n) = d -- is the same deuteron.  So the CS A = 3 set is used
/// here for the struck NEUTRON of a triton, which is the majority species and
/// the one that has a bound remnant.  (BeAGLE branches on the mass number
/// alone and hands 3He and 3H the same numbers with no mirror argument; the
/// argument is the point of writing it down.)
///
/// The MINORITY nucleon -- the proton in 3H, mirror of the neutron in 3He --
/// has no CS parameterization at all.  Its remnant is the unbound nn pair, so
/// its two-body fraction is IDENTICALLY ZERO whatever its n(k) is, and the
/// only thing left to choose is the SHAPE.  The choice made here, and the
/// one thing in this file that is not transcribed: the struck proton is drawn
/// from the SAME total n_0 + n_1 shape as the neutron.  It is the isoscalar
/// approximation, it is right in the k -> 0 and the SRC-dominated k -> large
/// limits (a pn pair carries both), and it is wrong in between by the
/// isovector difference the parameterization does not resolve.
/// `CiofiSimulaOptions::proton_n1_only` is the knob that spans it: true draws
/// the struck proton from the n_1 shape alone -- the hardest defensible
/// alternative, since a struck minority nucleon ALWAYS leaves a continuum.
/// Measured <k>: 126.0 MeV on the default, 171.0 MeV on the knob.
///
/// ---------------------------------------------------------------------------
/// THE CONTINUUM PAIR'S INTERNAL MOMENTUM.  Not from CS -- their n_1 is a
/// one-body distribution and says nothing about how the pair shares its own
/// relative momentum.  It comes from the form `breakup.cpp` already uses for
/// the nn remnant: the L = 0 Hulthen radial shape of `cluster.hpp` evaluated
/// at the pair's VIRTUAL-STATE pole kappa = hbar c / |a|, sampled as
/// q^2 |psi(q; kappa, beta)|^2.  Two poles, one per pair:
///
///     nn:  a_nn = -18.9 fm  -> `KAPPA_NN_VIRTUAL` = 10.44 MeV  (as today)
///     pn:  a_pn = -23.74 fm -> `KAPPA_PN_SINGLET` =  8.31 MeV  (1S0; the
///          3S1 channel of the pn pair is the BOUND deuteron and is the
///          OTHER branch, so the continuum pn pair is the singlet one)
///
/// The pair's invariant mass then follows exactly, M = sum sqrt(m_i^2 + q^2)
/// back to back in its own rest frame, and its excitation energy is
/// E = M - m_1 - m_2, which is what `TritonDraw::e_rel` reports.
///
/// ---------------------------------------------------------------------------
/// THE BeAGLE n_0-ONLY CAVEAT.  `DT_KFERMI` integrates 4 pi k^2 n(k) dk and
/// renormalises the result to unity, so BeAGLE never notices that its A = 3
/// and A = 4 tables are the n_0 piece alone.  It therefore samples the Fermi
/// momentum of 3He/3H and 4He from the GROUND-STATE-REMNANT distribution
/// renormalised to 1, dropping the 34.7 % / 20.0 % correlated continuum
/// entirely: its high-k tail is too soft by construction, and it has no
/// three-body breakup channel to put that strength in.  This file does not
/// renormalise n_0; it adds n_1 and lets the deficit be the branching.
///
/// ---------------------------------------------------------------------------
/// RNG STREAM DISCIPLINE.  `sample()` consumes EXACTLY
/// `TritonSpectralFunction::kUniformsPerSample` uniforms on every branch, so
/// the stream position after a draw does not depend on which channel came
/// out.  `ClusterBreakup` extends the same rule to the whole triton branch
/// (`breakup.hpp`).  This is the rule `CoherentXpomModel::draw` already
/// follows and it is what makes a run reproducible when a channel fraction
/// moves.

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "lipolgen/constants.hpp"
#include "lipolgen/rng.hpp"
#include "lipolgen/spectator.hpp"

namespace lipolgen {

// `HBARC_GEV_FM` (constants.hpp) is the one conversion between the Ciofi
// degli Atti-Simula parameterization (k in fm^-1, n in fm^3) and the
// library's GeV.  It lived here until 2026-09-01; it moved to constants.hpp
// so that it and `coherent.hpp`'s `GEV_PER_FM_INV` are ONE number.

/// nn virtual-state pole momentum [GeV]: hbar c / |a_nn| with a_nn = -18.9 fm.
/// The single definition of the scale the unbound nn remnant is split at.
/// (Lived in `breakup.hpp` until 2026-08-30; moved here so the pn pole can sit
/// beside it and `breakup.hpp` can keep using the name by including this.)
inline constexpr double KAPPA_NN_VIRTUAL = 0.0104399;

/// pn 1S0 virtual-state pole momentum [GeV]: hbar c / |a_pn| with the singlet
/// np scattering length a_pn = -23.74 fm.  The 3S1 pn channel is the BOUND
/// deuteron -- the two-body branch -- so a pn pair that comes out in the
/// CONTINUUM is the singlet one, and this is its scale.
inline constexpr double KAPPA_PN_SINGLET = 0.0083122;

/// One term A e^{-B k^2} / (1 + C k^2)^2 of the Ciofi degli Atti-Simula
/// parameterization, with A in fm^3 and B, C in fm^2.  C = 0 makes it the
/// pure Gaussian of their Eq. (75).
struct CsTerm {
  double a = 0.0, b = 0.0, c = 0.0;
};

/// The n_0 terms of CS Table A.1 for mass number `a` (2, 3 or 4).  Throws for
/// anything else.  A = 3 is the "3He (proton)" column; see the mirror
/// argument at the top of this file.
const std::vector<CsTerm>& cs_n0_terms(int a);

/// The n_1 terms: CS Eq. (76) for A = 3, Table A.3 for A = 4, EMPTY for
/// A = 2 (the deuteron has no excitable remnant).  Throws outside 2, 3, 4.
const std::vector<CsTerm>& cs_n1_terms(int a);

/// sum_i A_i e^{-B_i k^2} / (1 + C_i k^2)^2 at k [fm^-1], in fm^3.
double cs_sum(const std::vector<CsTerm>& terms, double k_fm);

/// int_0^{k_max} dk k^2 n(k) in the CS convention (NO 4 pi): the
/// spectroscopic factor of their Eq. (28).  Composite Simpson; the default
/// ceiling is the 10 fm^-1 = 1.973 GeV beyond which every set is at 1e-16.
double cs_norm(const std::vector<CsTerm>& terms, double k_max_fm = 10.0);

/// P(|k| > k_gev) of k^2 n(k), with k converted by `HBARC_GEV_FM`.
double cs_p_above(const std::vector<CsTerm>& terms, double k_gev,
                  double k_max_fm = 10.0);

/// <|k|> [GeV] of k^2 n(k).
double cs_mean_k(const std::vector<CsTerm>& terms, double k_max_fm = 10.0);

/// Which of the three channels a draw came out in.
enum class TritonChannel : std::uint8_t {
  NeutronD,        ///< struck n, remnant d -- BOUND, E = 0
  NeutronPnCont,   ///< struck n, remnant (p n) continuum
  ProtonNnCont     ///< struck p, remnant (n n) continuum
};

/// Human name of a channel ("n+d", "n+(pn)", "p+(nn)").
const char* triton_channel_name(TritonChannel c);

/// One draw of the spectral function, in the struck cluster's REST FRAME.
struct TritonDraw {
  TritonChannel channel = TritonChannel::NeutronD;
  bool struck_proton = false;
  double k = 0.0;             ///< struck-nucleon momentum [GeV]
  double cos_theta_k = 0.0;   ///< its polar cosine (the remnant recoils at -k)
  double phi_k = 0.0;
  double e_rel = 0.0;         ///< excitation energy of the remnant pair [GeV]
  double q_pair = 0.0;        ///< the pair's internal relative momentum [GeV]
  double m_remnant = 0.0;     ///< invariant mass of the two-nucleon remnant
};

/// The backend interface, in the library's usual shape: an analytic
/// implementation now (`CiofiSimulaTriton`) and a table-driven Faddeev/AV18
/// S(k, E) later, behind one signature.  Immutable after construction, so
/// every method is const and thread safe.
class TritonSpectralFunction {
 public:
  virtual ~TritonSpectralFunction() = default;

  /// Uniforms `sample()` consumes, on EVERY branch.  See the stream-discipline
  /// note at the top of this file.
  static constexpr int kUniformsPerSample = 6;

  /// n(k) of the struck species [GeV^-3], normalised to 1 over d^3k.
  /// `pdg` is 2212 or 2112.
  virtual double n_of_k(double k, int pdg) const = 0;

  /// P(residual BOUND | k, struck species).  Identically 0 for a struck
  /// proton -- there is no bound nn -- and n_0/(n_0 + n_1) for a struck
  /// neutron.
  virtual double p_two_body(double k, int pdg) const = 0;

  /// Density of the continuum pair's excitation energy at this k [GeV^-1],
  /// normalised on [0, e_max()].  Zero for the bound branch.
  virtual double e_rel_density(double k, double e_rel, int pdg) const = 0;

  /// Ceiling of the `e_rel_density` support for the pair left by `pdg`.
  virtual double e_max(int pdg) const = 0;

  /// One draw.  `p_proton` is P(the struck nucleon is the proton), which the
  /// CALLER computes -- `ClusterBreakup::proton_fraction`, i.e. the
  /// Z F2p : N F2n rule of docs/CONVENTIONS.md -- because the species draw is
  /// a structure-function statement and the spectral function has no F2
  /// backend.  (`interfaces_sketch.hpp` §2 sketched this as
  /// `sample(x, q2, m_s, rng)`; passing the already-formed fraction keeps the
  /// one species rule in the one place that owns it.)  Consumes exactly
  /// `kUniformsPerSample` uniforms whatever comes out.
  virtual TritonDraw sample(double p_proton, Rng& rng) const = 0;
};

/// Configuration of `CiofiSimulaTriton`.  Everything here is a documented
/// CHOICE; the transcribed coefficients are not options.
struct CiofiSimulaOptions {
  /// Sampling-grid ceiling on |k| [GeV], matching `BreakupOptions::k_max`.
  /// The CS total n(k) puts 3.33e-5 of its strength above 1.2 GeV and 3.1e-7
  /// above 1.5, so the truncation is invisible; what it buys is that the
  /// remnant can never be given more energy than the off-shell cluster has.
  double k_max = 1.2;
  /// Points of the |k| and |q| inverse-CDF grids.
  std::size_t n_grid = 4096;
  /// Ceiling on the pair's relative momentum |q| [GeV].
  double q_max = 1.2;
  /// Short-range scale of the Hulthen form the continuum pair's |q| is drawn
  /// from -- the run's own `cluster_beta`.
  double beta = BETA_DEFAULT;
  double kappa_nn = KAPPA_NN_VIRTUAL;   ///< nn virtual-state pole [GeV]
  double kappa_pn = KAPPA_PN_SINGLET;   ///< pn 1S0 virtual-state pole [GeV]
  /// Scale on n_1.  1 is the paper.  The knob on the one number the split
  /// rests on: 0 collapses to BeAGLE's n_0-only picture (every struck neutron
  /// leaves a bound d), 2 pushes S_0 to 0.48.
  double n1_scale = 1.0;
  /// Draw the struck PROTON from the n_1 shape alone instead of the total
  /// n_0 + n_1.  False (the isoscalar approximation) by default; see the
  /// minority-nucleon note at the top of this file.
  bool proton_n1_only = false;
};

/// The Ciofi degli Atti-Simula implementation: n_0 and n_1 as transcribed,
/// the 2-body/3-body split taken from their RATIO rather than assumed, and
/// the continuum pair split at its own virtual-state pole.
class CiofiSimulaTriton : public TritonSpectralFunction {
 public:
  explicit CiofiSimulaTriton(CiofiSimulaOptions opt = CiofiSimulaOptions());

  const CiofiSimulaOptions& options() const { return opt_; }

  double n_of_k(double k, int pdg) const override;
  double p_two_body(double k, int pdg) const override;
  double e_rel_density(double k, double e_rel, int pdg) const override;
  double e_max(int pdg) const override;
  TritonDraw sample(double p_proton, Rng& rng) const override;

  /// n_0(k) / n_1(k) of the A = 3 set at k [GeV], in the CS convention
  /// (fm^3, no 4 pi), with `n1_scale` applied to the second.  Public because
  /// the branching is a ratio of them and a test should be able to say so.
  double n0_cs(double k_gev) const;
  double n1_cs(double k_gev) const;

  /// S_0 and S_1 as this object actually integrates them, over [0, k_max].
  double s0() const { return s0_; }
  double s1() const { return s1_; }

 private:
  /// Inverse-CDF read of a cumulative grid, `np.interp(u, cdf, grid)` --
  /// the same read `MomentumSampler` and `ClusterBreakup` use.
  double draw_from(const std::vector<double>& cdf,
                   const std::vector<double>& grid, double u) const;
  std::vector<double> pair_cdf(double kappa, double* norm) const;

  CiofiSimulaOptions opt_;
  std::vector<double> kgrid_;      ///< |k| grid [GeV]
  std::vector<double> qgrid_;      ///< |q| grid [GeV]
  std::vector<double> cdf_k_n_;    ///< struck neutron: n_0 + n_1
  std::vector<double> cdf_k_p_;    ///< struck proton (see `proton_n1_only`)
  std::vector<double> cdf_q_nn_;   ///< nn pair at `kappa_nn`
  std::vector<double> cdf_q_pn_;   ///< pn pair at `kappa_pn`
  double qnorm_nn_ = 0.0, qnorm_pn_ = 0.0;   ///< int q^2 |psi|^2 dq per pair
  double s0_ = 0.0, s1_ = 0.0;
  double m_p_ = 0.0, m_n_ = 0.0, m_d_ = 0.0;
};

}  // namespace lipolgen

#endif  // LIPOLGEN_TRITON_SF_HPP
