#ifndef LIPOLGEN_BREAKUP_HPP
#define LIPOLGEN_BREAKUP_HPP

/// \file breakup.hpp
/// Tier T1: resolve the OFF-SHELL struck cluster of a tagged event into a
/// struck NUCLEON plus its on-shell partner spectator(s), so that
///
///     k + P_ion = k' + p_spec + sum(p_partner) + hadrons
///
/// holds exactly, event by event, with no caller-side hook and no
/// recoil correction anywhere.  `plans/05` 5.1 tier T1 / step 5.D:
/// "cluster-internal nucleon + partner spectators (alpha + p from d*;
/// alpha + d/nn from t*)".
///
/// ---------------------------------------------------------------------------
/// THE ONE RULE.  The impulse approximation says the SPECTATORS are physical
/// and the struck object is not: every partner is put ON SHELL at its
/// AME2020 mass and the struck nucleon takes whatever is left,
///
///     p_N,struck = P_X - sum(p_partner),      p_N,struck^2 != M_NUCLEON^2.
///
/// This is the same rule `tagged.hpp` already applies one level up (the
/// tagged spectator is on shell, the struck cluster P_X = P_ion - p_spec is
/// off shell) and it is BeAGLE's "never recoil-correct the light spectator".
/// Applying it twice makes the whole record close by construction: the only
/// off-shell object in the event is the one the hard process consumes.
///
/// ---------------------------------------------------------------------------
/// THE DEUTERON (6Li alpha tag; the d control's struck object is already a
/// nucleon and needs only relabelling).
///
///   * The internal relative momentum comes from the deuteron's OWN wave
///     function -- the `deuteron_channel()` `TaggedModel`, i.e. exactly the
///     S + D Hulthen forms of `cluster.hpp` at P_D = `P_D_DEUTERON` = 0.045
///     and the deuteron's kappa = sqrt(2 mu S_d) = 45.7 MeV, or, on
///     `BreakupOptions::source` = `VmcAV18`, the exact AV18 deuteron at its
///     own P_D = 0.057600.  The draw is the
///     m_S-dependent joint density |A_{m_sc}(m_S; k, khat)|^2, NOT the
///     spherical marginal: the D wave correlates khat with the spin, and that
///     correlation is what makes the struck nucleon's polarization and its
///     angular distribution consistent with each other.
///   * Which nucleon is struck is drawn `F2p : F2n` at the event's own
///     (x, Q2) -- the P1 species rule of `InclusiveGenerator::proton_fraction`
///     and `NucleonChoice::ByStructureFunctions`, applied here to Z = N = 1.
///     A flat 1 : 1 is its x-independent limit and is wrong wherever
///     F2n/F2p != 1.
///   * The struck nucleon's spin projection is drawn from the CG factor
///     p(m_1, m_2 | m_sc) = |<1/2 m_1 1/2 m_2 | 1 m_sc>|^2, so
///     <2 m_1>_{m_S} = <m_sc>_{m_S} = (1 - (3/2) P_D) m_S: the effective
///     nucleon polarization of the textbook deuteron comes out of the sampling
///     rather than being imposed on it.  `Particle::pol` carries the SAMPLED
///     +-1 helicity label (LHEF SPINUP), not the expectation value.
///     THAT IS ALSO THE CONSISTENCY GATE ON `source`: the dilution this draw
///     implies must equal the one the RATE is computed with,
///     `TaggedChannel::dis_target.eff_pol_*` -- 0.932495 against 0.9325 on
///     `Hulthen` and 0.913595 against 0.913600 on `VmcAV18`, a grid
///     quadrature against a closed form both times (T27).
///
/// ---------------------------------------------------------------------------
/// THE TRITON (7Li alpha tag).  TWO models, selected by
/// `BreakupOptions::triton_sf`:
///
///   null (the DEFAULT)  the sequential two-body Hulthen decay below, kept
///                       bit-for-bit;
///   non-null            the spectral function of `triton_sf.hpp`, which adds
///                       the THIRD channel the sequential model has no room
///                       for (struck n -> a (p n) CONTINUUM remnant) and takes
///                       its branching from Ciofi degli Atti-Simula's own
///                       n_0/(n_0 + n_1) rather than from a separation energy.
///
/// The spectral-function branch, in full:
///
///   * the species is drawn Z F2p : N F2n at the event's own (x, Q2) --
///     UNCHANGED, `proton_fraction` is still the one rule, and the breakup
///     hands the fraction to `TritonSpectralFunction::sample` rather than the
///     other way round, because the species draw is a structure-function
///     statement and the spectral function has no F2 backend;
///   * struck n -> d (bound, E = 0) with weight n_0/(n_0 + n_1);
///   * struck n -> (p n) continuum with weight n_1/(n_0 + n_1);
///   * struck p -> (n n) continuum always;
///   * the continuum pair is split back to back in its OWN rest frame at a
///     relative momentum drawn from the virtual-state Hulthen form -- the same
///     one this file already used for nn, at the nn pole for an nn pair and at
///     the pn 1S0 pole for a pn pair;
///   * the impulse-approximation rule is UNCHANGED (partners on shell, struck
///     nucleon absorbs the difference), so conservation is untouched;
///   * RNG STREAM DISCIPLINE: the branch consumes exactly
///     `TritonSpectralFunction::kUniformsPerSample` + 3 uniforms whatever the
///     channel -- six inside `sample()`, two for the pair's own direction
///     (drawn even when the remnant is the bound d and there is no pair) and
///     one for the polarization label.  The sequential model below does NOT
///     have that property (it consumes 5 or 8), which is one of the reasons
///     the new branch is a separate path rather than a patch.
///
/// The one number the whole split rests on, and where it came from, is at the
/// top of `triton_sf.hpp`.  Measured against the sequential model: <k> goes
/// from 133 MeV (t* -> n + d) / 145 MeV (t* -> p + nn) to 126 MeV over all
/// struck nucleons and 102 MeV on the bound-remnant channel alone.
///
/// THE SEQUENTIAL MODEL (the default) -- CRUDE AND FLAGGED (plans/05 step 5.D:
/// "t* remnant -> d or nn per the triton wave function -- crude, flagged").
///
/// A sequential two-body breakup with the triton's own AME2020 separation
/// energies, isotropic (S-wave) in the relative direction:
///
///   struck n:  t* -> n + d,      S = m_n + m_d   - m_t = 6.2572 MeV
///   struck p:  t* -> p + (nn),   S = m_p + 2 m_n - m_t = 8.4818 MeV
///
/// with kappa = sqrt(2 mu S) (mu the reduced mass of the two fragments) fed
/// into the L = 0 Hulthen radial form of `cluster.hpp` at the same beta as
/// the cluster model.  The unbound nn remnant is then split in its own rest
/// frame into two on-shell neutrons back to back at a relative momentum drawn
/// from the same form at the nn VIRTUAL-STATE pole
/// `KAPPA_NN_VIRTUAL` = 1/|a_nn| = 10.44 MeV (a_nn = -18.9 fm), so the pair
/// comes out nearly collinear -- the crudest defensible statement about a
/// system that has no bound state.  The remnant's invariant mass is then
/// 2 sqrt(m_n^2 + q^2) exactly, which is what the two emitted neutrons sum to.
///
/// What this is NOT: a Faddeev/AV18 three-body triton wave function, a
/// correlated (p, n, n) momentum distribution, or any final-state interaction
/// between the fragments.  The four-momentum balance is exact regardless --
/// the struck nucleon absorbs the whole difference -- so what the crudeness
/// costs is the SHAPE of the partner spectra, not conservation.
///
/// Polarization: the triton's per-nucleon effective polarizations are the
/// `TRITON()` `Ion`'s own, P_p = +0.86 and P_n = -0.028, so a struck nucleon
/// of a triton with projection m_S carries P_N = 2 m_S * P_species, sampled
/// to a +-1 `Particle::pol` label the same way the deuteron's is.
///
/// ---------------------------------------------------------------------------
/// STILL OPEN.  No final-state interaction of any fragment (no FSI, no
/// nuclear transparency, no formation time); no D-wave / tensor structure in
/// either triton branch (both are isotropic in the relative direction); and
/// the spectral-function branch factorises S_N(k, E) into n(k) times a
/// k-INDEPENDENT pair excitation, which a Faddeev/AV18 S(k, E) table would
/// not (`TritonSpectralFunction` is the interface such a table would enter
/// through).

#include <cstdint>
#include <memory>
#include <vector>

#include "lipolgen/cluster.hpp"
#include "lipolgen/constants.hpp"
#include "lipolgen/event.hpp"
#include "lipolgen/rng.hpp"
#include "lipolgen/sf.hpp"
#include "lipolgen/spectator.hpp"
#include "lipolgen/tagged.hpp"
#include "lipolgen/triton_sf.hpp"

namespace lipolgen {

// `KAPPA_NN_VIRTUAL` -- the nn virtual-state pole this file splits the
// unbound nn remnant at -- moved to `triton_sf.hpp` on 2026-08-30 so that it
// and the pn 1S0 pole sit together; the name is still reachable from here.

/// What the struck cluster of a tagged channel is made of.
enum class ClusterSpecies : std::uint8_t {
  Nucleon,    ///< already a nucleon (the d / 3He control channels)
  Deuteron,   ///< p + n            (6Li alpha tag)
  Triton      ///< p + n + n        (7Li alpha tag)
};

/// The species of the (Z, A) a `ClusterChannel::partner_Z/partner_A` names.
/// Throws for anything this module cannot break up.
ClusterSpecies cluster_species(int z, int a);

/// The struck cluster as `ClusterBreakup::resolve` needs it.
struct BreakupInput {
  ClusterSpecies species = ClusterSpecies::Nucleon;
  Vec4 p_cluster;            ///< P_X, OFF shell (`tagged.hpp`)
  int z = 1, a = 1;          ///< cluster charge / mass number
  double m_s = 0.0;          ///< struck-cluster spin projection
  double s_cluster = 0.5;    ///< the spin m_s is a projection of
  double x = 0.0, q2 = 0.0;  ///< event kinematics, for the F2p : F2n draw
  /// Per-nucleon effective polarizations of the cluster (the `Ion` slots);
  /// used for `ClusterSpecies::Triton` only.
  double eff_pol_p = 1.0, eff_pol_n = 1.0;
};

/// One resolved breakup.  `struck` and `partners` are ready to be appended to
/// an `Event`; `struck.p + sum(partners.p) == input.p_cluster` exactly.
struct BreakupResult {
  Particle struck;                  ///< Role::StruckNucleon, Status::Intermediate
  std::vector<Particle> partners;   ///< Role::PartnerSpectator, Status::Final
  double k = 0.0;                   ///< internal relative momentum [GeV]
  double cos_theta_k = 0.0;         ///< its polar cosine in the P_X rest frame
  double phi_k = 0.0;
  double q_nn = 0.0;                ///< the continuum pair's relative momentum
                                    ///< [GeV]; nn on both paths, pn only on
                                    ///< the spectral-function one
  double m_remnant = 0.0;           ///< invariant mass of the partner system
  double virtuality = 0.0;          ///< p_N,struck^2 - M_NUCLEON^2 [GeV^2]
  /// Which triton channel came out.  Meaningful only for
  /// `ClusterSpecies::Triton`; on the sequential path it is `NeutronD` or
  /// `ProtonNnCont` and never `NeutronPnCont`, which is the channel that
  /// model does not have.
  TritonChannel channel = TritonChannel::NeutronD;
  double e_rel = 0.0;               ///< excitation of the remnant pair [GeV]
};

/// Configuration of the breakup model.  Everything here is a documented
/// CHOICE; there is no other copy of any of these numbers.
struct BreakupOptions {
  double beta = BETA_DEFAULT;      ///< short-range scale of the radial forms
  double p_d = P_D_DEUTERON;       ///< deuteron D-state probability (0.045)
  /// WHICH DEUTERON THE T1 BREAKUP RESOLVES INTO (C5.5b, 2026-09-04).
  /// `Hulthen` (the default) is the analytic S + D pair at `p_d` above, bit
  /// for bit what this module always did; `VmcAV18` is the exact AV18
  /// deuteron of `fdeut.av18` at its own P_D = 0.057600 and then IGNORES
  /// `beta` and `p_d` -- the same table `deuteron_channel(.., VmcAV18)` and
  /// `DEUTERON_AV18()` use.
  ///
  /// IT EXISTS BECAUSE THIS MODULE CARRIES A DEUTERON WAVE FUNCTION OF ITS
  /// OWN.  `Pipeline` already forwards `cluster_beta` here for exactly that
  /// reason ("the T1 tier would describe a different nucleus from the T0
  /// one"), and until 2026-09-04 the wave-function FAMILY was not forwarded:
  /// a `--cluster-wave vmc` tagged-alpha run drew its struck-nucleon spin
  /// from the Hulthen deuteron at P_D = 0.045 while its rate was computed on
  /// the AV18 one, 2.0688 % apart, inside one run.
  ///
  /// THE TRITON PATH IS DELIBERATELY NOT AFFECTED: this tree has no AV18
  /// A = 3 wave function to switch to, so the 7Li sequential/spectral
  /// branches stay where they are on either setting.
  ClusterWaveSource source = ClusterWaveSource::Hulthen;
  double kappa_nn = KAPPA_NN_VIRTUAL;
  double k_max = 1.2;              ///< internal-momentum grid ceiling [GeV]
  std::size_t nk = 280, nc = 96;   ///< the deuteron model's own grid
  std::size_t n_grid = 4096;       ///< 1-D inverse-CDF grid of the triton step
  /// Unpolarized backend the `F2p : F2n` species draw is made on.  Null = the
  /// library's `ToyF2`; hand it the struck cluster's own kernel backend so the
  /// two draws stay consistent (docs/CONVENTIONS.md).
  std::shared_ptr<const UnpolSF> f2;
  /// The triton's spectral function.  NULL -- the default -- keeps the
  /// sequential Hulthen breakup above, bit for bit; a `CiofiSimulaTriton`
  /// replaces it with the three-channel model (`triton_sf.hpp`).  Ignored on
  /// every non-triton species.
  std::shared_ptr<const TritonSpectralFunction> triton_sf;
};

/// Struck cluster -> struck nucleon + partner spectator(s).
///
/// Immutable after construction (the deuteron amplitude grid and the two
/// triton inverse CDFs are built once), so `resolve` is const and thread
/// safe: it touches nothing but its arguments and the event's own `Rng`.
class ClusterBreakup {
 public:
  explicit ClusterBreakup(BreakupOptions opt = BreakupOptions());

  const BreakupOptions& options() const { return opt_; }
  /// The deuteron's own `TaggedModel` -- the internal wave function the
  /// `ClusterSpecies::Deuteron` branch draws from.
  const TaggedModel& deuteron_model() const { return *dmodel_; }
  /// The triton spectral function in force, or null on the sequential path.
  const TritonSpectralFunction* triton_sf() const { return tsf_.get(); }

  /// P(proton struck) = Z F2p / (Z F2p + N F2n) at (x, Q2); the flat Z/A when
  /// the backend is not positive there.
  double proton_fraction(double x, double q2, int z, int a) const;

  /// Resolve `in`.  Returns false only if the input is unusable (a
  /// non-positive cluster energy, an unknown species).  `rng` is the event's
  /// own counter-based stream, so the result is a pure function of it.
  bool resolve(const BreakupInput& in, Rng& rng, BreakupResult& out) const;

 private:
  /// One draw of |k| from k^2 |psi_L=0(k; kappa)|^2 by inverse CDF.
  double draw_k(const std::vector<double>& cdf, double u) const;
  static std::vector<double> build_cdf(const std::vector<double>& grid,
                                       double kappa, double beta);

  BreakupOptions opt_;
  std::shared_ptr<const UnpolSF> f2_;
  std::shared_ptr<const TritonSpectralFunction> tsf_;
  std::unique_ptr<TaggedModel> dmodel_;
  /// P(m_sc | m_S) of the deuteron model, keyed by 2*m_S + 2 (m_S = +1, 0, -1).
  std::vector<std::vector<double>> dpop_;
  std::vector<double> ms_deuteron_;
  /// The triton's two channels and the nn split, on one shared |k| grid.
  std::vector<double> grid_;
  std::vector<double> cdf_t_nd_;   ///< t* -> n + d
  std::vector<double> cdf_t_pnn_;  ///< t* -> p + (nn)
  std::vector<double> cdf_nn_;     ///< the nn relative momentum
  double m_nn_ref_ = 0.0;
};

}  // namespace lipolgen

#endif  // LIPOLGEN_BREAKUP_HPP
