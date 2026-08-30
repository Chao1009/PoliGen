#ifndef LIPOLGEN_TAGGED_HPP
#define LIPOLGEN_TAGGED_HPP

/// \file tagged.hpp
/// Tagged mode: spin-correlated (e', spectator) events.  C++17 port of
/// `evgen/polligen/tagged.py`.
///
/// The ion ground state is expanded in cluster relative partial waves,
///
///   |J M> = sum_L a_L sum_{m_L m_S} <L m_L S_c m_S|J M>
///             psihat_L(k) Y_L^{m_L}(khat) |S_c m_S>,
///
/// and everything the tagged observables need follows from the joint amplitude
///
///   A_{m_S}(M; k, khat) = sum_L a_L psihat_L(k) C_L(M, m_S) Y_L^{M-m_S}(khat).
///
/// * spectator density        n_M(k, khat) = sum_{m_S} |A_{m_S}|^2
/// * struck-cluster spin      p(m_S|M,k,khat) = |A_{m_S}|^2 / n_M   (diagonal
///                            truncation: coherences are dropped, documented)
/// * pair decomposition       p(m_struck, m_spec|m_S) from a second CG factor
///
/// Angular conventions: khat is measured in the SPIN frame (quantization axis
/// = z).  All densities are even in cos(theta_k), so the spectator-vs-struck
/// sign convention drops out.  `boost_spectator` rotates the spin frame into
/// the lab before the (longitudinal) beam boost.
///
/// DEPENDENCY NOTE.  This header deliberately does NOT include the inclusive
/// sampler.  The DIS side enters through the abstract `KinematicsSource`
/// below, which the P3 `InclusiveSampler` implements; the tagged model,
/// the boost and the struck-cluster four-vector are usable without it.

#include <cstddef>
#include <map>
#include <string>
#include <vector>

#include "lipolgen/beams.hpp"
#include "lipolgen/cluster.hpp"
#include "lipolgen/event.hpp"
#include "lipolgen/rng.hpp"
#include "lipolgen/spectator.hpp"

namespace lipolgen {

/// alpha-d D-state probability, chosen so the embedded-deuteron VECTOR
/// dilution 1 - (3/2) P_D reproduces the 0.87 of `b1_li6_from_deuteron`
/// (SCENARIO).  It is the DEFAULT because the Hulthen path is the
/// bit-compatibility path; `VMC_P_D_LI6` is the measured replacement and is
/// what `ClusterWaveSource::VmcAV18` uses.
inline constexpr double P_D_LI6 = 0.0867;
/// Deuteron D-state probability (AV18-like).
inline constexpr double P_D_DEUTERON = 0.045;

// ------------------------------------------- the VMC numbers, as DATA
//
// PROVENANCE.  R. B. Wiringa et al. (ANL), `momenta/li6_ad1.momentum` and
// `momenta/li7_at3.momentum`, AV18+UX variational Monte Carlo (1M and 500k
// samples, 22-Mar-14 and 12-Apr-24), fetched via the Internet Archive; see
// `data/vmc/README.md` for URLs and `docs/open_items/vmc_reconciliation.md`
// for the reconciliation against the 2004 `overlap_old` amplitudes.  The
// values below are the files' OWN printed normalizations
// `4*PI*TOTINT(RHO*K**2:K)/(2*PI)**3`, not a re-integration.

/// alpha-d D-state probability, VMC: 0.015861 / (0.80362 + 0.015861).
inline constexpr double VMC_P_D_LI6 = 0.015861 / (0.80362 + 0.015861);
/// alpha-d spectroscopic factor S_ad of 6Li (the file's total block).
inline constexpr double VMC_S_ALPHA_D_LI6 = 0.81971;
/// alpha-t spectroscopic factor S_at of the 7Li 3/2- GROUND state.
/// (`li7_at1.momentum`, the 1/2- excited state, gives 0.98683 -- not this.)
inline constexpr double VMC_S_ALPHA_T_LI7 = 1.0084;

/// The data files `ClusterWaveSource::VmcAV18` opens, relative to
/// `data_dir()` (`$LIPOLGEN_DATA_DIR`, else the compiled-in
/// `${CMAKE_SOURCE_DIR}/data`).
inline const char* const VMC_LI6_MOMENTUM = "vmc/momenta/li6_ad1.momentum";
inline const char* const VMC_LI7_MOMENTUM = "vmc/momenta/li7_at3.momentum";
inline const char* const VMC_LI6_OVERLAP = "vmc/li6_alpha_d/li6.ad";
inline const char* const VMC_LI7_OVERLAP = "vmc/li7_alpha_t/li7.at";

/// Struck-cluster DIS targets that are not beam species.  TRITON mirrors
/// `beams.HE3` under p <-> n and its constants are PER NUCLEON like every
/// other Ion slot, so g1(t) = 0.86 g1p + 2 * (-0.028) g1n.
const Ion& TRITON();
/// The free neutron as a DIS target (deuteron control channel).
const Ion& NEUTRON_TARGET();

/// Spin structure on top of a kinematic `ClusterChannel`.
struct TaggedChannel {
  ClusterChannel base;      ///< masses, separation energy, boost
  double j_ion = 1.0;
  double s_struck = 1.0;    ///< struck-cluster spin
  double s_spec = 0.0;      ///< spectator-cluster spin
  double s_channel = 1.0;   ///< coupled channel spin S_c
  std::vector<Wave> waves;
  Ion dis_target;           ///< structure-function target of the struck cluster
  std::string label;

  /// Throws std::runtime_error unless the wave probabilities sum to 1.
  void validate() const;
};

/// 6Li: DIS on the embedded deuteron, alpha spectator (S+D waves).
///
/// `source` defaults to `Hulthen`, which keeps every existing number
/// BIT-FOR-BIT; `VmcAV18` replaces both radial shapes with the ANL tables and
/// IGNORES `beta` and `p_d` (the D-state probability is then a property of
/// the wave function, `VMC_P_D_LI6`, not a knob).
TaggedChannel li6_alpha_channel(
    double beta = BETA_DEFAULT, double p_d = P_D_LI6,
    ClusterWaveSource source = ClusterWaveSource::Hulthen);
/// 7Li: DIS on the quasi-free triton, alpha spectator (pure P-wave).
/// `source` as above; `VmcAV18` ignores `beta`.
TaggedChannel li7_alpha_channel(
    double beta = BETA_DEFAULT,
    ClusterWaveSource source = ClusterWaveSource::Hulthen);
/// Deuteron control: DIS on the neutron, proton spectator (S+D).  The
/// Cosyn-Weiss tagged limit of the machinery.
TaggedChannel deuteron_channel(double beta = BETA_DEFAULT,
                               double p_d = P_D_DEUTERON);

/// Grid tables of the joint amplitude |A_{m_S}(M; k, cos theta_k)|^2.
///
/// The default grid (k_max = 1.2, nk = 280, nc = 96) is reproduced
/// bit-for-bit from the Python because `n_of_kc(M, k, c)` looks up the
/// NEAREST cell at or below the query point rather than interpolating, so the
/// grid IS part of the model.
class TaggedModel {
 public:
  explicit TaggedModel(TaggedChannel channel, double k_max = 1.2,
                       std::size_t nk = 280, std::size_t nc = 96);

  const TaggedChannel& channel() const { return channel_; }
  const std::vector<double>& k() const { return k_; }
  const std::vector<double>& c() const { return c_; }
  double dk() const { return dk_; }
  double dc() const { return dc_; }
  std::size_t nk() const { return k_.size(); }
  std::size_t nc() const { return c_.size(); }
  /// Channel-spin projections m_S, ordered +S_c ... -S_c.
  const std::vector<double>& m_struck_values() const { return ms_struck_; }
  /// The normalized radial table psihat_L(k) * sqrt(P_L) of wave L.
  const std::vector<double>& radial_table(int l) const;

  /// |A_{m_S}|^2 on the (k, c) grid, row-major [i_ms][ik * nc + ic].
  const std::vector<std::vector<double>>& amp2_table(double m_ion) const;

  /// n_M(k, cos theta_k), row-major [ik * nc + ic], normalized so that
  /// integral n k^2 dk dOmega = 1.
  const std::vector<double>& n_of_kc(double m_ion) const;
  /// The same, at an arbitrary (k, c), through the module's nearest-cell
  /// lookup (`np.clip(np.searchsorted(grid, v) - 1, 0, n - 2)`).
  double n_of_kc(double m_ion, double k, double c) const;

  /// p(m_S | M, k, c), row-major per m_S.
  std::vector<std::vector<double>> struck_populations(double m_ion) const;
  /// P(m_S | M): khat- and k-integrated channel-spin populations.
  std::vector<double> population_integrated(double m_ion) const;
  /// Grid quadrature of integral n_M k^2 dk dOmega (should be ~1; it is a
  /// finite-grid quadrature, so compare it at ~1e-6, not 1e-12).
  double norm(double m_ion) const;
  /// p(m_struck, m_spec | m_S) from the cluster-spin CG factor, indexed
  /// [i_m_struck][i_m_spec] with both orderings +S ... -S.
  std::vector<std::vector<double>> pair_populations(double m_s) const;

  /// khat,k-integrated <m_S>/S_c for the stretched state M = J by default.
  double vector_dilution() const { return vector_dilution(channel_.j_ion); }
  double vector_dilution(double m_ion) const;
  /// khat,k-integrated <3 m_S^2 - 2>, S_c = 1 only (throws otherwise).
  double tensor_dilution() const { return tensor_dilution(channel_.j_ion); }
  double tensor_dilution(double m_ion) const;
  /// <P2(cos theta_k)> of the spectator direction for ion state M.
  double p2_moment(double m_ion) const;
  /// <P2> for a fill with populations over M (linearity in rho).
  double p2_moment_mixture(const std::vector<double>& populations) const;

  /// Sample (k, cos theta_k, phi_k) from |A_{m_S}(M)|^2 k^2: a grid cell by
  /// inverse CDF, then uniform inside the cell, then a uniform azimuth.
  void sample_kc(double m_ion, double m_s, std::size_t n, Rng& rng,
                 std::vector<double>& k_out, std::vector<double>& c_out,
                 std::vector<double>& phi_out) const;
  /// One draw.
  void sample_kc_one(double m_ion, double m_s, Rng& rng, double& k_out,
                     double& c_out, double& phi_out) const;

 private:
  std::size_t ms_index(double m_s) const;
  const std::vector<double>& cell_cdf(double m_ion, double m_s) const;
  /// The three grid builds, called ONCE per ion projection by the
  /// constructor.  See the note on the cache members below.
  std::vector<std::vector<double>> build_amp2(double m_ion) const;
  std::vector<double> build_n(const std::vector<std::vector<double>>& a2) const;
  std::vector<double> build_cdf(const std::vector<double>& a2_ms) const;

  TaggedChannel channel_;
  std::vector<double> k_, c_;
  double dk_ = 0.0, dc_ = 0.0;
  std::vector<double> ms_struck_;
  std::map<int, std::vector<double>> rad_;   ///< L -> normalized radial table
  // C3.  These were `mutable` and filled lazily on first use, which made
  // every accessor a write and left thread safety resting on the caller
  // warming every state first (`Pipeline`'s constructor did; nothing else had
  // to).  They are now built in full by the CONSTRUCTOR and never written
  // again, so the object is immutable after construction, every accessor is a
  // pure lookup, and no lock is needed on any path.
  //
  // A mutex round the lazy caches was the alternative and was MEASURED
  // against this one: `sample_kc + boost + route` runs at 5.18-5.36 Mev/s
  // locked against 5.11-5.40 Mev/s eager, i.e. the lock is free at this
  // granularity (one 26880-entry `lower_bound` dominates the call).  The
  // eager build wins on the guarantee instead: it is the only one that makes
  // the object immutable, so no caller has to know about a warm-up protocol.
  // Its cost is paid once, in the constructor, and `Pipeline` was already
  // paying it in its warm-up loop.
  std::map<long long, std::vector<std::vector<double>>> amp2_;
  std::map<long long, std::vector<double>> n_;
  std::map<long long, std::vector<double>> cdf_;
};

/// Lab kinematics of the spectator cluster.  (k, c, phi_k) are spherical
/// components in the SPIN frame; for a tilted quantization axis they are
/// rotated to the lab (R_z(phi_S) R_y(theta_S)) before the longitudinal beam
/// boost.  `phi_spec` is the spectator's LAB azimuth -- what the planar Roman
/// Pot rectangle needs, and deliberately not called "phi", which is the DIS
/// azimuth on an event record.
struct SpectatorLab {
  double pT = 0.0;
  double theta = 0.0;
  double p_lab = 0.0;
  double R = 0.0;
  double xL = 0.0;
  double kx = 0.0, ky = 0.0, kz = 0.0;   ///< lab-oriented rest-frame components
  double phi_spec = 0.0;
  double e_lab = 0.0, pz_lab = 0.0;      ///< the boosted four-vector's pieces
};

SpectatorLab boost_spectator(const TaggedChannel& channel, double k, double c,
                             double phi_k, double p_per_nucleon,
                             double theta_s = 0.0, double phi_s = 0.0);

/// The Fermi-smeared struck cluster, as the T2 (PYTHIA) tier consumes it.
///
/// CONVENTION (the one physics choice this file makes that the Python does not
/// state).  The spectator is put ON SHELL at its physical mass and the struck
/// cluster takes whatever is left of the ion:
///
///     P_X = P_ion - p_spec,     p_spec^2 = m_spec^2 exactly,
///
/// so the struck cluster is the OFF-SHELL one.  This is the standard impulse
/// approximation of spectator tagging (Cosyn-Weiss; BeAGLE's "never
/// recoil-correct the light spectator" rule) and it is the only choice that
/// conserves four-momentum event by event while keeping the tagged fragment's
/// measured mass exact.  Its virtuality M_X^2 - m_free^2 is negative and is
/// carried explicitly rather than hidden.
///
/// Light-front variables are quoted in the ION REST FRAME, where they are
/// invariant under the longitudinal beam boost:
///
///     alpha_s = A_beam (E_k - k_z) / m_beam,     p_T,s = |k_T|,
///
/// normalized so that alpha_s -> A_spec * (m_spec A_beam)/(m_beam A_spec) ~
/// A_spec at k = 0 and alpha_s + alpha_X = A_beam.  `p_per_nucleon_eff` is
/// P_X,z / A_partner, the per-nucleon momentum the hard process should be
/// injected at.
struct StruckCluster {
  Vec4 p_ion;         ///< head-on frame, ion along +z
  Vec4 p_spectator;   ///< on shell, m^2 = m_spec^2
  Vec4 p;             ///< P_ion - p_spectator, OFF shell
  double m2 = 0.0;          ///< P_X^2 [GeV^2]
  double m_free = 0.0;      ///< free-cluster mass [GeV]
  double virtuality = 0.0;  ///< m2 - m_free^2 [GeV^2], negative in the IA
  double alpha_s = 0.0;     ///< light-front fraction of the SPECTATOR
  double alpha_x = 0.0;     ///< A_beam - alpha_s
  double pt_s = 0.0;        ///< |k_T| of the spectator [GeV]
  double p_per_nucleon_eff = 0.0;
};

StruckCluster struck_cluster(const TaggedChannel& channel,
                             const SpectatorLab& lab, double p_per_nucleon);

/// The wave-function tensor asymmetry of a spin-1 channel vs k,
/// A_zz^wf = (n_+1 + n_-1 - 2 n_0) / (n_+1 + n_-1 + n_0), at one cos(theta_k)
/// cell.  Returns one value per k cell; NaN where the denominator vanishes.
std::vector<double> azz_tensor_curve(const TaggedModel& model, std::size_t ic);
/// The ACCEPTANCE-WEIGHTED version: `weights` is an (nk, nc) row-major table
/// (e.g. `acceptance_weights`) against which both sums are integrated over
/// cos(theta_k).  Weights concentrated in one cell reproduce the cell curve.
std::vector<double> azz_tensor_curve_weighted(const TaggedModel& model,
                                              const std::vector<double>& weights);

/// eps(k, cos theta_k) on the model grid: the fraction of lab azimuths at
/// which a spectator of that (k, c) is Roman-Pot accepted.  theta, R and pT
/// depend on (k, c) alone -- the beam boost is longitudinal -- so the azimuth
/// enters only through the rectangular envelope and a uniform phi average is
/// the exact marginal.  Row-major (nk, nc).
std::vector<double> acceptance_weights(const TaggedModel& model,
                                       double p_per_nucleon,
                                       const Optics& optics,
                                       const std::string& pot_config,
                                       std::size_t n_phi = 64,
                                       double theta_s = 0.0,
                                       double phi_s = 0.0);

/// One ion spin category of a run plan, spelled locally so that this header
/// does not depend on the bookkeeping module.
struct IonFill {
  std::string name;
  double j = 1.0;
  std::vector<double> populations;   ///< m = +J ... -J, sums to 1
  int lam_e = 0;
  double pe = 0.0;
  double theta_s = 0.0, phi_s = 0.0;
};

/// The DIS side of a tagged event, as an interface.
///
/// The P3 inclusive sampler implements this: given the struck-cluster spin
/// projection m_S it draws (x, Q2, y, phi) from its own accepted phase space
/// and reports the accepted cross section of that pure state.  Defining it
/// here rather than including the sampler keeps `tagged` buildable on its own
/// and makes the DIS model swappable (toy source in the tests).
struct KinematicsSource {
  virtual ~KinematicsSource() = default;
  /// Draw `n` DIS kinematics conditioned on the struck-cluster projection
  /// `m_struck`.  The four output vectors are resized to `n`.
  virtual void sample(double m_struck, int lam_e, double pe, std::size_t n,
                      Rng& rng, std::vector<double>& x, std::vector<double>& q2,
                      std::vector<double>& y, std::vector<double>& phi) = 0;
  /// The same, additionally reporting the sampler's accepted-CELL index of
  /// every draw (`Kinematics::cell`, the column Mode-W reweighting needs).
  /// The default forwards to `sample` and reports "unknown" (-1), so a source
  /// that has no cell notion needs no change; `InclusiveKinematicsSource`
  /// overrides it.  It exists as a separate entry point rather than a default
  /// argument because a virtual default argument is bound statically.
  virtual void sample_cells(double m_struck, int lam_e, double pe,
                            std::size_t n, Rng& rng, std::vector<double>& x,
                            std::vector<double>& q2, std::vector<double>& y,
                            std::vector<double>& phi, std::vector<int>& cell) {
    sample(m_struck, lam_e, pe, n, rng, x, q2, y, phi);
    cell.assign(n, -1);
  }
  /// Accepted cross section [pb] of the pure struck-cluster state, which sets
  /// the relative (M, m_S) rates.  Default 1: unweighted rates.
  virtual double sigma_tot_pb(double m_struck, int lam_e, double pe) const {
    (void)m_struck; (void)lam_e; (void)pe;
    return 1.0;
  }
};

/// One tagged event, before it is turned into an `Event` record.
struct TaggedEvent {
  double x = 0.0, q2 = 0.0, y = 0.0, phi = 0.0;
  int cell = -1;              ///< the DIS source's accepted-cell index
  double m_ion = 0.0, m_struck = 0.0;
  double k = 0.0, cos_theta_k = 0.0, phi_k = 0.0;
  SpectatorLab lab;
  int route = kRouteLost;
};

/// Spin-correlated (e', spectator) events for one tagged channel.
///
/// Joint sampling, exactly the Python's: M from the fill populations weighted
/// by the (m_S-summed) DIS rates; then (m_S, k, khat) from |A_{m_S}(M)|^2
/// times the DIS rate of m_S; then (x, Q2, phi) from the `KinematicsSource`
/// conditioned on m_S; then the lab boost and the far-forward route.
class TaggedSampler {
 public:
  /// `optics` defaults to the Yellow Report HIGH-ACCEPTANCE envelope OF THIS
  /// CONFIGURATION (per-configuration and anisotropic since 2026-08-28), not
  /// to the legacy proton-derived 73 microrad: no module-level default can
  /// know the beam.
  TaggedSampler(const TaggedModel& model, double p_per_nucleon,
                KinematicsSource* dis = nullptr);
  TaggedSampler(const TaggedModel& model, double p_per_nucleon,
                KinematicsSource* dis, const Optics& optics,
                const std::string& pot_config);

  const TaggedModel& model() const { return *model_; }
  const Optics& optics() const { return optics_; }
  const std::string& pot_config() const { return pot_config_; }
  double p_per_nucleon() const { return p_u_; }

  /// Accepted tagged cross section [pb] for an ION spin category.  Tagging
  /// probability is 1 in the two-cluster model; purity and efficiency come
  /// from the routing.
  double sigma_tot_pb(const IonFill& fill) const;

  /// Per-(M, m_S) rates, row-major [i_M * n_ms + i_ms].
  std::vector<double> rates(const IonFill& fill) const;

  /// `n` tagged events for one ion spin category.
  std::vector<TaggedEvent> sample_category(const IonFill& fill, std::size_t n,
                                           Rng& rng) const;
  /// One event.  THE THROUGHPUT PATH: it allocates nothing per event -- the
  /// ion-projection ladder is a member and the four DIS draw buffers the
  /// `KinematicsSource` fills are thread-local scratch (P6).
  TaggedEvent sample_one(const IonFill& fill, const std::vector<double>& rate_cdf,
                         Rng& rng) const;
  /// Ion projections M, ordered +J ... -J.
  const std::vector<double>& m_ion_values() const { return ms_ion_; }
  /// Normalized cumulative of `rates`, for `sample_one`.
  std::vector<double> rate_cdf(const IonFill& fill) const;

  /// ADD the spectator and the struck cluster to an event whose beams and e'
  /// are already set, and fill the tagging block of `Event::kin`
  /// (k, cos_theta_k, phi_k, alpha_s, pt_s) and the struck-cluster spin label.
  /// The beam ion is taken from the event's `Role::BeamIon` particle if there
  /// is one, else built from the channel; nothing already in the record is
  /// modified except `kin` and `spin.m_struck` / `spin.m_ion`.
  void fill_event(Event& ev, const TaggedEvent& te) const;

 private:
  const TaggedModel* model_;
  double p_u_;
  KinematicsSource* dis_;
  Optics optics_;
  std::string pot_config_;
  std::vector<double> ms_ion_;   ///< m_values(j_ion), hoisted out of the loop
};

/// PDG code of a nuclide: 10-digit ion code 10LZZZAAAI, with the proton and
/// the neutron given their particle codes instead.
int nuclide_pdg(int z, int a);

}  // namespace lipolgen

#endif  // LIPOLGEN_TAGGED_HPP
