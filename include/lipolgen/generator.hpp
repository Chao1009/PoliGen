#ifndef LIPOLGEN_GENERATOR_HPP
#define LIPOLGEN_GENERATOR_HPP

/// \file generator.hpp
/// T0 event assembly: sampled (x, Q2, phi, m) + a beam configuration ->
/// `lipolgen::Event`.
///
/// FRAME.  Head-on: the ion along +z at A * p_u GeV, the electron along -z at
/// E_e GeV.  The 25 mrad crossing angle is a LAB transform applied downstream
/// (`reco.head_on_to_lab`), never here.  Polar angles theta are measured from
/// +z (the ion direction), so the scattered electron sits at theta ~ pi.
///
/// AZIMUTH -- the one convention this header pins.
///   `Kinematics::phi` is the azimuth of the SCATTERED ELECTRON about the ion
///   (+z) beam axis in the HEAD-ON frame:
///        phi = atan2(p_y(e'), p_x(e')) in [0, 2 pi),
///   and the physics azimuth the master formula modulates is
///        phi' = phi - phi_S,
///   phi_S being the azimuth of the ion alignment axis n(theta_S, phi_S) in
///   the SAME frame (`SpinLabels::phi_s`).  For a massless target this phi'
///   is EXACTLY the covariant Bacchetta et al. (JHEP 02 (2007) 093) azimuth
///   of the spin axis about the virtual photon measured from the lepton
///   plane -- `reco.azimuth_wrt_lepton_plane` is written and tested to that
///   sign convention (`reco.py`: "phi_S = phi_e - phi_s exactly for a
///   massless target").  Only 2 phi enters the tensor observables, so the
///   orientation convention cannot change a result; the point of fixing it
///   is that the generator and the reconstruction agree.
///
/// FOUR-MOMENTUM CONSERVATION.  T0 conserves momentum on the PER-NUCLEON
/// subsystem, which is the system the per-nucleon (x, Q2) describe:
///        k + P_N = k' + X       exactly,
/// with P_N the struck-nucleon four-vector and X a single pseudo-particle
/// (pdg 92) carrying the whole hadronic final state.  The default P_N is the
/// nucleon at rest in the ion rest frame, i.e. moving with the beam at
/// p_u along +z and on shell at M_NUCLEON.  `GeneratorConfig::target` is the
/// hook P4/P6 replace with a Fermi-smeared nucleon drawn from the cluster /
/// spectator model: whatever it returns is used verbatim, so X follows it and
/// the identity above still holds event by event.  The (A-1) remnant is NOT
/// written at T0 -- the tagged tier (`tagged.hpp`, `spectator.hpp`) owns it,
/// and it is what turns the per-nucleon balance into a whole-nucleus one.

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "lipolgen/bookkeeping.hpp"
#include "lipolgen/event.hpp"
#include "lipolgen/sampler.hpp"

namespace lipolgen {

/// 10-digit nuclear PDG code 10LZZZAAAI (1000030060 = 6Li).
int nuclear_pdg(int z, int a);

/// The struck-nucleon four-vector in the head-on frame.  `draw` carries the
/// event's (x, Q2, y, phi, m) and `rng` is the event's own counter-based
/// stream, so a Fermi-smeared implementation stays reproducible.
using TargetSampler = std::function<Vec4(const EventDraw& draw, Rng& rng)>;

/// Callback that consumes finished events without any of them being stored.
using EventSink = std::function<void(const Event&)>;

struct GeneratorConfig {
  /// Null = the default target: one nucleon on shell at M_NUCLEON, at rest in
  /// the ion rest frame, i.e. p_z = p_u (the beam's momentum per nucleon).
  TargetSampler target;
  /// 0 = draw the struck nucleon's species from the STRUCTURE FUNCTIONS:
  /// P(proton) = Z F2p(x, Q2) / (Z F2p(x, Q2) + N F2n(x, Q2)), evaluated at
  /// the event's own (x, Q2) on the sampler kernel's own unpolarized backend
  /// (P1).  A flat Z : A draw -- what this did before 2026-08-29 -- is the
  /// x-independent limit of that and is wrong wherever F2n/F2p is: at x = 0.5
  /// it hands 6Li a proton half the time where the inclusive rate wants
  /// 0.616, and the error grows towards the valence edge.
  /// Set 2212 / 2112 to pin the species.
  int struck_nucleon_pdg = 0;
  /// Emit the virtual photon as a status-3 documentation particle.
  bool with_virtual_photon = true;
  Channel channel = Channel::Inclusive;
  /// Whole-nucleus ion mass [GeV].  <= 0 means the physical AME2020 mass.
  double ion_mass = -1.0;
};

/// Turns sampled kinematics into `lipolgen::Event` records.
class InclusiveGenerator {
 public:
  explicit InclusiveGenerator(std::shared_ptr<const InclusiveSampler> sampler,
                              GeneratorConfig config = GeneratorConfig());

  const InclusiveSampler& sampler() const { return *sampler_; }
  const GeneratorConfig& config() const { return config_; }

  /// Head-on beams.  The electron is massless by convention (`reco.py`
  /// `beam_fourvectors`); the ion carries its physical nuclear mass.
  const Vec4& beam_electron() const { return beam_e_; }
  const Vec4& beam_ion() const { return beam_ion_; }
  int ion_pdg() const { return ion_pdg_; }

  /// The struck nucleon, from `GeneratorConfig::target` or the default.
  Vec4 target_nucleon(const EventDraw& draw, Rng& rng) const;

  /// P(proton) at (x, Q2) -- Z F2p / (Z F2p + N F2n) on the kernel's own
  /// unpolarized backend.  Falls back to Z/A if F2 is not positive there.
  double proton_fraction(double x, double q2) const;

  /// Scattered electron in the head-on frame from (x, y, phi).
  Vec4 scattered_electron_p4(double x, double y, double phi) const;

  /// One T0 event.  `rng` is the event's stream (already used by the sampler
  /// for the kinematics); the generator draws only the nucleon species and
  /// whatever `GeneratorConfig::target` needs from it.
  Event make_event(const SpinCategory& cat, const RunPlan& plan,
                   const EventDraw& draw, Rng& rng, std::uint64_t number,
                   int run, int bunch) const;
  /// The same, reconstructed IN PLACE in `ev`: `Event::reset()` keeps the
  /// record's heap capacity, so a loop over a reused `Event` allocates
  /// nothing per event.  The value-returning overload is this one plus a
  /// move.
  void make_event(const SpinCategory& cat, const RunPlan& plan,
                  const EventDraw& draw, Rng& rng, std::uint64_t number,
                  int run, int bunch, Event& ev) const;

  /// Generate a whole run plan for `total_lumi_pb` and hand every event to
  /// `sink` as it is made -- nothing is stored.  Category k uses bunch index
  /// k of the (seed, run) stream and event counters 0 .. N_k-1, so the whole
  /// run is reproducible from (seed, run) alone.  Returns the total number
  /// of events emitted.
  std::uint64_t run(const RunPlan& plan, double total_lumi_pb,
                    std::uint64_t seed, const EventSink& sink,
                    std::uint64_t run_number = 1, bool poisson = true) const;

  /// Fixed-statistics variant: `n_total` events split across the categories
  /// in proportion to lumi_fraction * sigma (the mixture a luminosity run
  /// would produce on average).
  std::uint64_t run_n(const RunPlan& plan, std::uint64_t n_total,
                      std::uint64_t seed, const EventSink& sink,
                      std::uint64_t run_number = 1) const;

  /// Accepted cross section [pb] per category, in plan order.
  std::vector<double> sigma_per_category(const RunPlan& plan) const;

 private:
  std::shared_ptr<const InclusiveSampler> sampler_;
  GeneratorConfig config_;
  Vec4 beam_e_, beam_ion_;
  int ion_pdg_ = 0;
  double ion_mass_ = 0.0;
  double proton_fraction_ = 0.0;
};

}  // namespace lipolgen

#endif  // LIPOLGEN_GENERATOR_HPP
