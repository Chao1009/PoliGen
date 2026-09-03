#ifndef LIPOLGEN_PIPELINE_HPP
#define LIPOLGEN_PIPELINE_HPP

/// \file pipeline.hpp
/// The integration layer: one configuration object, one `Pipeline`, and
/// `lipolgen::Event` records for every physics channel the library carries.
///
/// It is the only header that knows about all of `sampler`, `generator`,
/// `tagged`, `coherent` and `bookkeeping` at once; every module below it stays
/// independent (docs/CONVENTIONS.md).
///
/// ---------------------------------------------------------------------------
/// WHAT X IS, PER CHANNEL.  Exactly one pseudo-particle (pdg 92,
/// `Role::HadronicX`) carries the whole hadronic final state, and it is always
/// defined as "whatever is left", so the balance below is exact to rounding
/// event by event and needs no recoil correction anywhere:
///
///   Inclusive        k + P_N   = k' + X                (PER-NUCLEON balance,
///                    X = k + P_N - k'.  This is `generator.hpp`'s T0 identity
///                    and is what `InclusiveGenerator` already writes: the
///                    (A-1) remnant is not written at all, so the WHOLE-nucleus
///                    balance and the total charge are deliberately NOT closed.
///
///   Tagged           k + P_ion = k' + p_spec + X       (WHOLE-NUCLEUS balance)
///                    The spectator cluster is put ON SHELL at its physical
///                    (AME2020) mass, the struck cluster takes the remainder
///                    P_X = P_ion - p_spec and is therefore OFF shell
///                    (`tagged.hpp`, the standard impulse approximation), and
///                    X = k + P_X - k' = k + P_ion - p_spec - k'.  Both the
///                    off-shell struck cluster (`Role::StruckCluster`,
///                    status 3) and the spectator (`Role::Spectator`,
///                    status 1) are written, so charge closes too:
///                    q(e') + q(spec) + q(X) = q(e) + Z_ion.
///
///   Coherent         k + P_ion = k' + P_recoil + X     (WHOLE-NUCLEUS balance)
///                    The ground-state nucleus stays intact and on shell at
///                    `nuclear_mass(Z, A)` with p_T = sqrt(|t|) and
///                    p_z = (1 - x_P) A p_u; X = k + P_ion - k' - P_recoil is
///                    the diffractive system and carries charge 0.
///
/// WHAT THE T2 (PYTHIA) TIER MUST CONSUME.  For the tagged channels the hard
/// process belongs to the STRUCK CLUSTER, not to the ion and not to a nucleon
/// at rest in the ion frame: `Role::StruckCluster` carries P_X (off shell,
/// virtuality M_X^2 - m_free^2 < 0), `Kinematics::alpha_s` / `pt_s` carry the
/// light-front variables in the ion rest frame, and
/// `StruckCluster::p_per_nucleon_eff = P_X,z / A_partner` is the per-nucleon
/// momentum a nucleon inside it should be drawn around.  `struck_cluster_of()`
/// below reconstructs the whole `StruckCluster` record from a finished event.
/// `PipelineConfig::hadronizer` is the hook the T2 tier binds itself to; it is
/// called with the finished T0 event and the event's own counter-based stream,
/// which is the signature of `PythiaBridge::hadronize(Event&, Rng&)`.
///
/// DESIGN NOTE -- WHY THE TAGGED PATH DOES NOT GO THROUGH
/// `GeneratorConfig::target`.  That hook returns a struck NUCLEON and lets
/// `InclusiveGenerator` close the per-nucleon balance around it.  Three things
/// a tagged event needs do not fit through it: (a) the spectator cluster has to
/// be WRITTEN, and a `TargetSampler` returning one four-vector cannot add a
/// particle; (b) the DIS kinematics of a tagged event are drawn from the STRUCK
/// CLUSTER's own accepted phase space conditioned on m_S -- a different
/// `InclusiveSampler` from the ion-level one the generator holds -- so the
/// generator's sampler would be the wrong one; (c) the T2 tier wants the
/// cluster, not a nucleon.  The pipeline therefore COMPOSES
/// `InclusiveGenerator` (inclusive channel, used verbatim so that path stays
/// bit-identical to `InclusiveGenerator::run_n`) with `TaggedSampler` /
/// `CoherentSampler`, which own their own `fill_event`, and adds X itself.
///
/// DETERMINISM.  Every event is a PURE FUNCTION of its global index: category
/// k and local counter j are recovered from the index, the stream is
/// `Rng(seed, run, k, j)`, and nothing else is carried between events.  So
/// `event(i)` is thread-safe and `for_each(sink, nthreads)` emits exactly the
/// same events in exactly the same order for any thread count.

#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "lipolgen/beams.hpp"
#include "lipolgen/bookkeeping.hpp"
#include "lipolgen/breakup.hpp"
#include "lipolgen/coherent.hpp"
#include "lipolgen/event.hpp"
#include "lipolgen/fsi.hpp"
#include "lipolgen/generator.hpp"
#include "lipolgen/rc.hpp"
#include "lipolgen/sampler.hpp"
#include "lipolgen/spectator.hpp"
#include "lipolgen/tagged.hpp"
#include "lipolgen/xsec.hpp"

namespace lipolgen {

// ============================================================ DIS adapter

/// `tagged.hpp`'s `KinematicsSource` implemented on top of an
/// `InclusiveSampler` built on the STRUCK CLUSTER's structure functions.
///
/// This is the C++ shape of `polligen.tagged.TaggedSampler.__init__` +
/// `_pure_category`: the DIS side of a tagged event is the Step-5.A inclusive
/// master formula evaluated on the struck cluster (embedded deuteron for the
/// 6Li alpha tag, quasi-free triton for the 7Li alpha tag, free neutron for
/// the deuteron control), and the (x, Q2, phi) draw is CONDITIONED on the
/// struck-cluster projection m_S through a PURE spin category of spin S_c.
///
/// ONE SAMPLER, ONE PURE CATEGORY PER m_S -- not one sampler per m_S.  The
/// (x, Q2) cell grid and the per-cell unpolarized cross sections do not depend
/// on the spin state at all; what depends on m_S is the per-cell
/// (w_avg, a1, a2) triple, and `InclusiveSampler` already caches exactly that
/// per spin state behind a mutex (`state_tables`).  Building N samplers would
/// rebuild the same grid N times and change no number.
///
/// The pure category is built at theta_S = phi_S = 0 whatever the ION fill's
/// axis is, which is what the Python does: the struck-cluster quantization
/// axis is the channel-spin axis of the two-cluster expansion, and the ion
/// axis enters the event only through the spectator rotation in
/// `boost_spectator`.
class InclusiveKinematicsSource : public KinematicsSource {
 public:
  /// `s_channel` is the struck cluster's channel spin S_c (1 for the embedded
  /// deuteron, 1/2 for the triton and for the struck neutron).
  InclusiveKinematicsSource(std::shared_ptr<const InclusiveSampler> sampler,
                            double s_channel);

  const InclusiveSampler& sampler() const { return *sampler_; }
  const std::shared_ptr<const InclusiveSampler>& sampler_ptr() const {
    return sampler_;
  }
  double s_channel() const { return s_channel_; }

  void sample(double m_struck, int lam_e, double pe, std::size_t n, Rng& rng,
              std::vector<double>& x, std::vector<double>& q2,
              std::vector<double>& y, std::vector<double>& phi) override;

  /// The same, reporting the accepted-cell index of each draw
  /// (`Kinematics::cell`).  `sample` forwards to this one, so there is a
  /// single draw path and the two cannot diverge.
  void sample_cells(double m_struck, int lam_e, double pe, std::size_t n,
                    Rng& rng, std::vector<double>& x, std::vector<double>& q2,
                    std::vector<double>& y, std::vector<double>& phi,
                    std::vector<int>& cell) override;

  double sigma_tot_pb(double m_struck, int lam_e, double pe) const override;

  /// The scattered electron of a drawn (x, y, phi) in the head-on frame, at
  /// the STRUCK CLUSTER sampler's own s (which is 4 E_e p_u and therefore the
  /// same as the ion-level one -- `BeamConfig::s_per_nucleon`).
  Vec4 scattered_electron_p4(double x, double y, double phi) const;

  /// Build every pure-state plan of (lam_e, pe) up front.  The per-state cache
  /// is mutex-protected, so this is an optimization for threaded generation,
  /// not a correctness requirement.
  void warm(int lam_e, double pe) const;

 private:
  const InclusiveSampler::CategoryPlan& plan_for(double m_s, int lam_e,
                                                 double pe) const;
  struct Key {
    long long m2;   ///< 2 m_S, exact
    int lam_e;
    double pe;
    bool operator<(const Key& o) const;
  };

  std::shared_ptr<const InclusiveSampler> sampler_;
  double s_channel_;
  std::vector<double> ms_;
  mutable std::mutex mutex_;
  mutable std::map<Key, InclusiveSampler::CategoryPlan> plans_;
};

/// Options of the struck-cluster kernel / DIS source.
struct StruckClusterOptions {
  /// Put an INCLUSIVE b1 in the struck cluster's kernel.  Default FALSE, which
  /// is what `evgen/scripts/money_tagged_azz.py` uses
  /// (`InclusiveKernel(beams.DEUTERON, b1_func=None)`): in the impulse
  /// approximation the embedded deuteron's b1 IS the k-integral of the
  /// m-dependent spectator density `TaggedSampler` already draws from, so an
  /// inclusive b1 on top of it counts the same physics twice.  Set it to
  /// reproduce `evgen/tests/test_tagged.py`'s `li6_sampler` fixture
  /// (`b1_func=toy_b1`), which is what the RATE-asymmetry identity
  ///     Azz(sigma_+1, sigma_0, sigma_-1) = tensor_dilution * <Azz>_sigma
  /// is written against.
  bool inclusive_b1 = false;
  /// Optional double-helicity-flip Delta slot; empty = zero (the default of
  /// both Python scripts).
  SFFunc3 delta_func;
  /// Acceptance window of the struck cluster's sampler.
  Scenario scenario = Scenario();
  /// (x, Q2) grid of the struck cluster's sampler.
  SamplerGrid grid = SamplerGrid();
  /// Compute the transverse-vector amplitude a1 (needs g2).
  bool with_perp = false;
};

/// The kernel the struck cluster of `channel` is generated with:
/// `channel.dis_target` (DEUTERON / TRITON / NEUTRON_TARGET) on the library's
/// default toy F2 and toy g1 backends, with the tensor slots empty unless
/// `opt.inclusive_b1`.
std::shared_ptr<const InclusiveKernel> struck_cluster_kernel(
    const TaggedChannel& channel, const StruckClusterOptions& opt = {});

/// The whole DIS side of a tagged channel: the struck cluster's kernel, its
/// `InclusiveSampler` at the ION beam's per-nucleon momentum (the Python's
/// `beams.BeamConfig(E_e, ch.dis_target, ion_momentum_per_nucleon)`), wrapped
/// as a `KinematicsSource`.
std::shared_ptr<InclusiveKinematicsSource> make_struck_cluster_source(
    const TaggedChannel& channel, const BeamConfig& ion_config,
    const StruckClusterOptions& opt = {});

// ============================================================ configuration

/// Physics channel of a run.
enum class PipelineChannel : std::uint8_t {
  Inclusive,          ///< e + A -> e' + X on the per-nucleon subsystem
  TaggedLi6Alpha,     ///< 6Li(e, e' alpha)X -- DIS on the embedded deuteron
  TaggedLi7Alpha,     ///< 7Li(e, e' alpha)X -- DIS on the quasi-free triton
  TaggedDeuteronP,    ///< d(e, e' p)X control -- the Cosyn-Weiss tagged limit
  CoherentLi6         ///< e + 6Li -> e' + X + 6Li(g.s.)
};
const char* pipeline_channel_name(PipelineChannel c);
/// True for the three spectator-tagged channels.
bool is_tagged(PipelineChannel c);
/// The beam species a channel implies ("6Li", "7Li", "d"); empty for
/// `Inclusive`, which takes whatever `PipelineConfig::isotope` says.
std::string channel_isotope(PipelineChannel c);

/// Which far-forward envelope the ROUTE LABEL on an event is computed with.
enum class OpticsChoice : std::uint8_t {
  YellowReportHighAcceptance,  ///< `yr_optics(..., high_acceptance = true)`
  YellowReportHighDivergence,  ///< `yr_optics(..., high_acceptance = false)`
  Tagging,                     ///< the lithium tagging optics (per-config levers)
  TaggingLegacyLevers,         ///< the same, priced with the 18x275 pot levers
  Custom                       ///< `PipelineConfig::optics` verbatim
};

/// The lithium TAGGING optics of Report 1 Section 6.1 as an `Optics`.
///
/// TABULATED, not computed.  The working point is the maximum of
/// (tagged fraction) x (luminosity) over a de-squeeze scan of beta*_x, which
/// `polli_fastsim.farforward.tagging_optics_point` performs on a 400-point
/// grid; this is that scan's answer for the six (species, configuration)
/// combinations the library carries, transcribed at full double precision.
/// It is the ONLY definition of these numbers in LiPolGen (there is no other
/// C++ copy to disagree with).
///
/// `legacy_levers` reproduces every tagging number published before
/// 2026-08-29, when the dispersive smearing was priced with the 18x275 pot
/// levers at every configuration instead of each configuration's own.
Optics tagging_optics(const std::string& ion_name, double p_per_nucleon,
                      bool legacy_levers = false, double n_sigma = 10.0);

/// The envelope an `OpticsChoice` selects for a beam.
Optics optics_for(OpticsChoice choice, const std::string& ion_name,
                  double p_per_nucleon, double n_sigma = 10.0);

/// Fidelity tier of the final state a run writes (DEVELOPMENT_PLAN.md 2).
///
///   T0  the struck cluster stays a single off-shell pseudo-particle
///       (`Role::StruckCluster`) and nothing inside it is resolved.
///   T1  the struck cluster is broken up into a struck NUCLEON plus its
///       on-shell partner spectator(s) (`breakup.hpp`), so the record names
///       the object the hard process actually consumes and
///       k + P_ion = k' + p_spec + sum(p_partner) + X holds exactly.
///
/// T1 is the DEFAULT on the tagged channels.  T0 is kept because it is what
/// every tagged number published before this tier existed was made with, and
/// because the partner spectators are pure addition: no T0 quantity moves.
/// The switch does nothing on the inclusive and coherent channels, which have
/// no struck cluster to resolve.
enum class Tier : std::uint8_t { T0, T1 };

/// Which spectral function the 7Li alpha tag's T1 triton breakup draws from
/// (`PipelineConfig::triton_sf`; CLI `--triton-sf {hulthen,ciofi-simula}`).
///
///   Hulthen      the DEFAULT: the sequential two-body Hulthen decay of
///                `breakup.hpp`, bit-for-bit what every published 7Li number
///                was made with;
///   CiofiSimula  the Ciofi degli Atti-Simula spectral function of
///                `triton_sf.hpp` -- the k-dependent n_0/(n_0 + n_1)
///                branching and the THIRD channel (struck n -> a (p n)
///                continuum) the sequential model has no room for.  The
///                `Pipeline` constructor builds the `CiofiSimulaTriton`
///                ITSELF, at the run's own `cluster_beta` and the breakup's
///                own `k_max`, so the continuum pair's q-shape and every
///                other radial form share ONE beta (docs/CONVENTIONS.md: no
///                physics number is defined twice).
///
/// Ignored on every non-triton species, and ignored entirely when a C++
/// caller has already put an object in `BreakupOptions::triton_sf`.
enum class TritonSfChoice : std::uint8_t { Hulthen, CiofiSimula };

/// Which b1 backend fills the INCLUSIVE kernel's rank-2 slot when
/// `PipelineConfig::kernel` is null (CLI `--b1-model`).
///
///   Miller          the DEFAULT, and bit-for-bit what every published
///                   inclusive tensor number was made with:
///                   `Li6B1(MillerB1)` through `LI6_B1_RANK2_TRANSFER` and
///                   `LI6_B1_PER_NUCLEON` (2/6).  This IS the run PLAN's
///                   "toy": today's default b1_func is `Li6B1(MillerB1)`,
///                   which is what `toy_b1` reaches.  Pinned at rtol 1e-12
///                   by `validation/reference/b1_default_li6.json`
///                   (tests/test_b1_nuclear.cpp T9).
///   Cdks            the same 6Li rank-2 transfer on the CDKS convolution
///                   camp (`CdksB1`, the digitized PRD 95 (2017) 074036
///                   Fig. 4 column): |b1| two orders of magnitude smaller
///                   than Miller's below x ~ 0.1, COMPARABLE above it (peak
///                   |x b1| 1.7e-4 against Miller's 4.3e-4 at Q2 = 2.5, and
///                   4x LARGER at x = 0.3 with the opposite sign), and a
///                   different sign structure.  Miller (HERMES-like) and
///                   CDKS (convolution) are different CAMPS for b1_d and the
///                   library does not adjudicate between them -- say which
///                   one a plot used.  INCLUSIVE CHANNEL ONLY and 6Li ONLY,
///                   like `Li6Convolution`: it is `Li6B1`'s 6Li transfer,
///                   and elsewhere the flag would not reach the rate but
///                   WOULD reach the metadata (`validate()`).
///   Li6Convolution  `b1_nuclear.hpp`'s four-term alpha-d convolution
///                   (`Li6ConvolutionB1`, design_D_b1_li6.md).  INCLUSIVE
///                   CHANNEL ONLY and 6Li ONLY -- see `validate()`; a spin-1
///                   test is not enough, because the deuteron is spin-1 too
///                   and the alpha-d densities, N_ad and the 2/6 and 4/6
///                   counting factors are all specific to 6Li.  On a TAGGED
///                   channel the alpha-d density is already in the event
///                   weight (`StruckClusterOptions::inclusive_b1`), so
///                   folding this in on top would count the same physics
///                   three times.
///
///   WARNING: `Li6Convolution` has NOT fully passed its A = 2 validation
///   gate (b1_nuclear.hpp's header block, and
///   docs/open_items/run_2026-09-02/phase_D_gate.md).  It is opt-in and
///   band-only: never quote a single row, always {0, 1, 2} x b1.
enum class B1Model : std::uint8_t { Miller, Cdks, Li6Convolution };
const char* b1_model_name(B1Model m);

/// The T2 hand-off.  Called once per finished T0 event with the event's own
/// counter-based stream, AFTER every T0 particle (including the off-shell
/// struck cluster and the spectator) is in place.  Deliberately the signature
/// of `PythiaBridge::hadronize(Event&, Rng&)`.  Must be re-entrant if
/// threaded generation is used.
using HadronizerHook = std::function<void(Event& ev, Rng& rng)>;

struct PipelineConfig {
  // --- beams -------------------------------------------------------------
  std::string isotope = "6Li";   ///< "6Li", "7Li", "d", ... (`ion_by_name`)
  int beam_config = 1;           ///< index into `default_configs(isotope)`

  // --- physics -----------------------------------------------------------
  PipelineChannel channel = PipelineChannel::Inclusive;
  /// Acceptance window.  The default is the GENERATOR window of
  /// docs/DEVELOPMENT_PLAN.md 5 (Q2 >= 0.7, y in [0.004, 0.985], W2 >= 8),
  /// looser than any analysis window so events can migrate in.
  Scenario scenario = generator_scenario(Scenario());
  SamplerGrid grid = SamplerGrid();
  /// Inclusive channel only: the kernel.  Null = `default_inclusive_kernel`.
  std::shared_ptr<const InclusiveKernel> kernel;
  /// Emit the virtual photon as a status-3 documentation particle.
  bool with_virtual_photon = true;

  // --- statistics --------------------------------------------------------
  std::uint64_t seed = 20260713;
  std::uint64_t run = 1;
  /// Integrated luminosity [pb^-1].  Used when `n_events == 0`.
  double lumi_pb = 0.0;
  /// Fixed total event count; splits across categories in proportion to
  /// lumi_fraction * sigma.  Exclusive with `lumi_pb`.
  std::uint64_t n_events = 0;
  /// Luminosity mode only: Poisson-fluctuate the per-category counts.
  bool poisson = true;

  // --- far-forward routing ------------------------------------------------
  OpticsChoice optics_choice = OpticsChoice::YellowReportHighAcceptance;
  Optics optics;                 ///< used only with `OpticsChoice::Custom`
  double n_sigma = 10.0;
  /// Machine configuration the over-rigid branch tests against.  Empty = the
  /// `yr_config_key` of the beam.
  std::string pot_config;

  // --- tagged channels ----------------------------------------------------
  double cluster_beta = BETA_DEFAULT;  ///< short-range scale of the radial waves
  double p_d = P_D_LI6;                ///< D-state probability (6Li alpha tag)
  /// Which family of radial forms the lithium alpha-tag channels use.
  /// `Hulthen` is the default and keeps every published number bit-for-bit;
  /// `VmcAV18` swaps in the ANL VMC tables and then IGNORES `cluster_beta`
  /// and `p_d` for those two channels (`docs/CONVENTIONS.md`).  The deuteron
  /// control channel is always Hulthen -- there is no VMC d -> p+n cluster
  /// table, the deuteron IS the cluster.
  ClusterWaveSource cluster_wave = ClusterWaveSource::Hulthen;
  StruckClusterOptions struck;         ///< struck-cluster DIS options
  /// Final-state interaction of the DIS debris X with the tagged spectator,
  /// as a PER-EVENT WEIGHT on `Event::weight` (never a shift of any
  /// four-vector; spin independent by construction -- see fsi.hpp).  `Off`
  /// (the default) is today's plane-wave impulse approximation bit for bit.
  /// Tagged channels only: `validate()` refuses it elsewhere.
  PipelineFsi fsi = PipelineFsi::Off;
  /// sigma_XN [mb] the FSI weight is built at.  40 = free hadron; the
  /// documented BAND is 20-40 mb and the 20 mb end is arguably the realistic
  /// one at EIC formation lengths.  BAND it (run both ends as a systematic);
  /// never quote one row alone (fsi.hpp).
  double fsi_sigma_mb = 40.0;
  /// Tensor-sector radiative corrections as OPT-IN, WEIGHT-ONLY families
  /// (rc.hpp): the band `rc_tensor_lo`/`rc_tensor_hi` on the tensor part of
  /// the rate, and the 6Li radiative tails `rc_tail` (elastic + unpolarised
  /// quasi-elastic).  They go on `Event::rc_weights` and NOT on
  /// `Event::weight` -- a systematic variation and a background, not a
  /// correction -- so `Off` (the default) is today bit for bit, in every
  /// four-vector, in `Event::weight` and in the RNG stream.  Applies on every
  /// channel except `CoherentLi6`; the tail is additionally identically 1 on
  /// the tagged channels (design_C_tensor_rc.md sec. 1.5).  `validate()` does
  /// NOT refuse the coherent channel -- `RcModel::applies()` decides and the
  /// run prints why, so a channel scan need not special-case `--rc`.
  PipelineRc rc = PipelineRc::Off;
  RcOptions  rc_options;          ///< knobs ONLY -- there is no mode on it

  // --- coherent channel ---------------------------------------------------
  CoherentScenario coherent;
  /// |t| truncation [GeV^2].  0.2 since 2026-08-29 (P4): the azimuthal weight
  /// 1 + c2 cos 2(phi_t - phi_S) goes NEGATIVE beyond |t| = 0.245 at
  /// P_zz = -2, and the deformation input is digitized only to |t| = 0.30.
  /// The `Pipeline` constructor refuses a range in which it would.
  double coherent_t_max = COHERENT_T_MAX_DEFAULT;
  /// How x_P -- and with it the diffractive mass M_X -- is drawn per event.
  /// The default guarantees M_X >= 1 GeV, i.e. a TIMELIKE X (C1).
  CoherentXpomModel coherent_xpom;
  /// Draw phi_t from the modulated density instead of carrying it as an event
  /// weight (`CoherentSampler::set_weighted_azimuth`).
  bool coherent_weighted_azimuth = false;

  // --- tiers --------------------------------------------------------------
  /// Fidelity tier of the TAGGED final state; `Tier::T1` by default (the
  /// struck cluster is resolved into a nucleon + partner spectators).  See
  /// `Tier` above and `breakup.hpp`.
  Tier tier = Tier::T1;
  /// Options of the T1 cluster breakup.  `BreakupOptions::beta` and
  /// `f2` are overwritten by the `Pipeline` constructor with the run's own
  /// `cluster_beta` and the struck cluster's own unpolarized backend, so the
  /// three draws that share them cannot disagree.
  BreakupOptions breakup;
  /// Triton spectral function of the 7Li alpha tag's T1 breakup
  /// (`TritonSfChoice` above; the model itself is `triton_sf.hpp`).
  /// `CiofiSimula` makes the constructor fill `breakup.triton_sf` with a
  /// `CiofiSimulaTriton` built at the run's own `cluster_beta` unless the
  /// caller already set one.
  TritonSfChoice triton_sf = TritonSfChoice::Hulthen;
  HadronizerHook hadronizer;
  // C4 -- CLOSED 2026-08-30.  `hadronize_coherent`, the opt-in that let a
  // hadronizer run on the coherent channel through the broken inclusive
  // fallback, IS GONE together with the fallback itself.  A coherent event
  // now names its own T2 target: `make_coherent` writes the POMERON
  // P_IP = P_ion - P_recoil (`Role::Pomeron`, status 3), which
  // `PythiaBridge` hadronizes on a PYTHIA Pomeron beam (id 990) exactly the
  // way it hadronizes a struck nucleon (docs/PYTHIA_BRIDGE.md sec. 12).  The
  // whole record conserves; the knob that stayed is
  // `PythiaBridgeOptions::coherent_t2` = {Pomeron, Off}.

  /// C6.  Multiply the luminosity of every category by
  /// `Optics::lumi_fraction` -- the share of the machine luminosity the
  /// configured far-forward working point actually delivers.  Default TRUE.
  ///
  /// COUNTS ONLY, NEVER CROSS SECTIONS: `sigma_per_category_pb()` is
  /// share-invariant by the same rule that keeps `Scenario::run_share` out of
  /// it (bookkeeping.hpp).  The tagging optics buy their acceptance by
  /// de-squeezing beta*_x, which costs luminosity -- 0.147 of it at 6Li
  /// 5x41 -- and quoting a tagged yield at the Yellow Report luminosity while
  /// routing it through the de-squeezed envelope overstates the sample by
  /// 7-10x.  Ignored in fixed-`n_events` mode, where the count is given.
  bool apply_optics_lumi_fraction = true;

  // --- inclusive b1 backend ------------------------------------------------
  /// Which b1 backend `default_inclusive_kernel` fills the rank-2 slot with
  /// (`B1Model` above; CLI `--b1-model`).  `Miller` is the DEFAULT and is
  /// today bit for bit.  Ignored when `kernel` is non-null -- and `validate()`
  /// REFUSES that combination for anything but `Miller`, because a
  /// caller-supplied kernel silently wins over the flag and the two would be
  /// saying contradictory things.
  B1Model b1_model = B1Model::Miller;
  /// The MANDATORY 100 % band of design_D_b1_li6.md sec. 4.3: multiplies the
  /// WHOLE b1 of the `Cdks` and `Li6Convolution` backends (all four terms of
  /// the latter).  Run 0 / 1 / 2 and quote the envelope.  Q(6Li) =
  /// -0.0818(17) fm^2 (`LI6_QUADRUPOLE_FM2`, TUNL A = 6, 1998CE04 -- Pyykko's
  /// compilation gives -0.0806(6); the repository constant is TUNL's)
  /// against Q_d = +0.2859(3) fm^2: the alpha-d relative
  /// D wave enters the closest measured observable with the OPPOSITE sign to
  /// the deuteron's own D state and nearly cancels it, so a single row is not
  /// a result.  `Miller`'s numbers are already published and the band is not
  /// applied to them, so `validate()` REFUSES anything but 1.0 there rather
  /// than let the metadata record a variation that did not run.
  double b1_band_scale = 1.0;
  /// `Li6ConvolutionOptions::w_alpha_d_dwave`: the knob on terms (2d) AND
  /// (2a) together -- they are one physical effect (the alpha-d orbital
  /// alignment) and are scaled as one.  The SHAPE variant of the same worry
  /// the band expresses; report 0 / 1 / 2 separately from the band.  Read
  /// ONLY by the `Li6Convolution` branch, so `validate()` refuses anything but
  /// 1.0 on `Miller` and `Cdks` -- same provenance rule as the band.
  double b1_alpha_d_dwave_weight = 1.0;

  /// Throws std::runtime_error on an inconsistent configuration.
  void validate() const;
};

/// The kernel the INCLUSIVE channel uses when `PipelineConfig::kernel` is
/// null: the Miller deuteron b1 through the 6Li rank-2 transfer plus a
/// 1e-2 discovery-scale Delta for a spin-1 ion, empty tensor slots for
/// anything else (there is no published rank-2 input for 7Li).  Same choice as
/// `examples/generate_inclusive.cpp`.
std::shared_ptr<const InclusiveKernel> default_inclusive_kernel(const Ion& ion);
/// The same with the b1 backend chosen (`PipelineConfig::b1_model` and the two
/// knobs).  The one-argument form above is exactly
/// `default_inclusive_kernel(ion, B1Model::Miller, 1.0, 1.0)`, so
/// `examples/generate_inclusive.cpp` and every existing caller are untouched.
///
/// `band_scale` multiplies the whole b1 of the non-Miller backends;
/// `w_alpha_d` is `Li6ConvolutionOptions::w_alpha_d_dwave`.  The Delta slot is
/// the SAME on every `B1Model` -- it is not part of this choice.  The
/// `Li6Convolution` branch shares the kernel's own `ToyF2` with the
/// convolution's `unpol`, so "the kernel's UnpolSF" is literally one object.
std::shared_ptr<const InclusiveKernel> default_inclusive_kernel(
    const Ion& ion, B1Model model, double band_scale = 1.0,
    double w_alpha_d = 1.0);

// ============================================================ event helpers

/// The four-momentum residual of the balance the event's channel guarantees
/// (see the file header): per-nucleon for `Channel::Inclusive`, whole-nucleus
/// for every tagged and the coherent channel.  Should be zero to rounding.
Vec4 momentum_residual(const Event& ev);
/// The four-momentum that ENTERED that balance, for a relative comparison.
Vec4 momentum_scale(const Event& ev);
/// Total final-state charge minus the charge that entered the same balance.
double charge_residual(const Event& ev);

/// Far-forward route of the tagged spectator, or of the coherent intact
/// recoil, of a finished event -- `kRouteLost` when the event carries
/// neither.
///
/// It is RECOMPUTED from the record rather than stored on it, which is what
/// every Python script does ("the sample is drawn ONCE and re-routed per
/// optics", `money_tagged_azz.folded_asymmetry`): the route is a detector
/// statement about a fixed physics event, so one generated sample can be
/// priced at several envelopes without regenerating.  The reconstruction is
/// exact -- the beam boost is longitudinal, so
///   R = (|p_frag| / Z_frag) / (|p_beam| / Z_beam),  theta = atan2(pT, pz),
///   phi = atan2(py, px)
/// are the same numbers `boost_spectator` computed.  Use
/// `route_of(ev, pipeline.optics(), pipeline.pot_config())` for the route at
/// the run's own configured optics.
int route_of(const Event& ev, const Optics& optics,
             const std::string& pot_config = "18x275");
/// `rp_accepted(route_of(...))`: the Roman-Pot mask (main window + near-beam
/// tail) `tagged.rp_accepted` applies.
bool rp_tagged(const Event& ev, const Optics& optics,
               const std::string& pot_config = "18x275");

/// Rebuild the `StruckCluster` record the T2 tier consumes from a finished
/// tagged event.  Returns false (and leaves `out` untouched) if the event
/// carries no `Role::StruckCluster`.
bool struck_cluster_of(const Event& ev, const TaggedChannel& channel,
                       StruckCluster& out);

// ============================================================ the pipeline

/// One configured run: config + run plan -> `lipolgen::Event` records.
///
/// Everything is resolved in the constructor (grids, per-state amplitude
/// tables, spectator amplitude tables, per-category counts), so generation is
/// allocation-light and every public generation entry point is const and
/// thread-safe.
class Pipeline {
 public:
  Pipeline(PipelineConfig config, RunPlan plan);

  const PipelineConfig& config() const { return cfg_; }
  const RunPlan& plan() const { return plan_; }
  const BeamConfig& beam_config() const { return beams_; }
  const Optics& optics() const { return optics_; }
  const std::string& pot_config() const { return pot_config_; }

  /// The DIS sampler the channel draws (x, Q2, phi) from: the ION-level one
  /// for `Inclusive` and `CoherentLi6`, the STRUCK-CLUSTER one for the tagged
  /// channels.
  const InclusiveSampler& dis_sampler() const { return *dis_sampler_; }
  /// Null outside the tagged channels.
  const TaggedModel* tagged_model() const { return model_.get(); }
  const TaggedSampler* tagged_sampler() const { return tsampler_.get(); }
  const TaggedChannel* tagged_channel() const;
  /// The T1 cluster-breakup model; null outside the tagged channels or when
  /// the run is configured at `Tier::T0`.
  const ClusterBreakup* breakup() const { return breakup_.get(); }
  /// The FSI weight model of the run; null when `PipelineConfig::fsi` is
  /// `Off`.  Expose it so a run can print what it used (`sigma_eff_mb`,
  /// `survival`, the profile numbers).
  const GlauberFsiWeight* fsi_weight() const { return fsi_.get(); }
  /// The RC weight model of the run; null when `PipelineConfig::rc` is `Off`.
  /// Expose it so a run can print the band, the form-factor provenance and
  /// the clipped fractions (rc.hpp).
  const RcModel* rc_model() const { return rc_.get(); }
  /// The tier this run actually writes (`Tier::T0` on channels with no
  /// struck cluster, whatever the configuration says).
  Tier tier() const { return tier_; }
  /// Null outside `CoherentLi6`.  One sampler per category (they differ by
  /// P_zz and phi_S); index is the category index.
  const CoherentSampler* coherent_sampler(std::size_t category) const;

  // --- cross-section bookkeeping -----------------------------------------

  /// Accepted cross section [pb] per category, in plan order.  SHARE
  /// INVARIANT: no `lumi_fraction` and no `Scenario::run_share` reaches it
  /// (bookkeeping.hpp, `test_run_share.py`).
  const std::vector<double>& sigma_per_category_pb() const { return sigma_; }
  /// sum_k lumi_fraction_k * sigma_k [pb] -- the mixture one unit of the
  /// plan's luminosity buys.
  double sigma_pb() const;
  /// Integrated luminosity [pb^-1] per category (zero in fixed-count mode).
  /// Already carries `Optics::lumi_fraction` when
  /// `PipelineConfig::apply_optics_lumi_fraction` is on.
  const std::vector<double>& lumi_per_category_pb() const { return lumi_; }
  /// The luminosity factor the optics contributes to the COUNTS: the
  /// configured `Optics::lumi_fraction`, or 1 when the knob is off.
  double optics_lumi_factor() const;

  // --- event counts -------------------------------------------------------

  /// Events per category, in plan order.
  const std::vector<std::uint64_t>& counts() const { return counts_; }
  /// Total number of events this run produces.
  std::uint64_t size() const { return total_; }
  /// Category index of a global event index.
  std::size_t category_of(std::uint64_t index) const;

  // --- generation ---------------------------------------------------------

  /// Event `index` of the run.  A PURE function of the index: safe to call
  /// from any thread, in any order, any number of times.  Throws if
  /// `index >= size()`.
  Event event(std::uint64_t index) const;
  /// The same, reusing `out`'s storage.
  void event(std::uint64_t index, Event& out) const;

  /// Pull API: fills `ev` with the next event and returns true, false at the
  /// end of the run.  The cursor is the only mutable state on the object.
  bool next(Event& ev);
  void rewind() { cursor_ = 0; }
  std::uint64_t cursor() const { return cursor_; }

  /// Streaming API: hands every event to `sink` in index order and stores
  /// nothing.  Returns the number of events emitted.
  std::uint64_t for_each(const EventSink& sink) const;

  /// Streaming over an index RANGE `[first, last)`, one reused `Event` and
  /// no storage -- the pull-free half of `generate_range`.  Prefer this to
  /// `generate_range` whenever the records are consumed and dropped: a
  /// 4096-record buffer costs cache, not allocations (measured 2.19 M ev/s
  /// against 2.87 M streaming, inclusive 6Li at the mid configuration).
  std::uint64_t for_each_range(const EventSink& sink, std::uint64_t first,
                               std::uint64_t last) const;

  /// Threaded streaming.  `nthreads` workers each build a `chunk`-sized block
  /// of the index range; the blocks are handed to `sink` from the CALLING
  /// thread in index order, so `sink` never needs a lock and the output is
  /// bit-identical to `for_each(sink)` for any `nthreads`.
  std::uint64_t for_each(const EventSink& sink, unsigned nthreads,
                         std::size_t chunk = 2048) const;

  /// Events `[first, last)` into `out`, RESIZED but not cleared: the records
  /// already in `out` are reconstructed in place, so a loop that calls this
  /// on the same buffer allocates nothing after the first chunk.
  /// Thread-safe; any decomposition of a range gives the same events.
  void generate_range(std::uint64_t first, std::uint64_t last,
                      std::vector<Event>& out) const;

 private:
  // All three reconstruct IN PLACE (`Event::reset()` keeps the record's heap
  // capacity), so `event(i, ev)` on a reused `ev` -- and therefore
  // `for_each` and `generate_range` -- allocates nothing per event.
  void make_inclusive(std::size_t k, std::uint64_t local, std::uint64_t index,
                      Rng& rng, Event& ev) const;
  void make_tagged(std::size_t k, std::uint64_t local, std::uint64_t index,
                   Rng& rng, Event& ev) const;
  void make_coherent(std::size_t k, std::uint64_t local, std::uint64_t index,
                     Rng& rng, Event& ev) const;
  void label_event(Event& ev, std::size_t k, std::uint64_t index) const;
  void add_beams(Event& ev, const SpinCategory& cat, double m_ion) const;
  void add_hadronic_x(Event& ev, const Vec4& p_x, double charge,
                      int mother) const;
  /// e' from the ION-level sampler's s (inclusive / coherent paths).
  Vec4 gen_scattered_electron(double x, double y, double phi) const;

  PipelineConfig cfg_;
  RunPlan plan_;
  BeamConfig beams_;
  Optics optics_;
  std::string pot_config_;

  // inclusive
  std::shared_ptr<const InclusiveSampler> dis_sampler_;
  std::unique_ptr<InclusiveGenerator> gen_;
  std::vector<InclusiveSampler::CategoryPlan> cplan_;

  // tagged
  std::unique_ptr<TaggedChannel> channel_;
  std::shared_ptr<TaggedModel> model_;
  std::shared_ptr<InclusiveKinematicsSource> dis_source_;
  std::unique_ptr<TaggedSampler> tsampler_;
  std::shared_ptr<GlauberFsiWeight> fsi_;   ///< null when cfg_.fsi == Off
  std::unique_ptr<ClusterBreakup> breakup_;
  Tier tier_ = Tier::T0;
  ClusterSpecies cluster_species_ = ClusterSpecies::Nucleon;
  std::vector<IonFill> fills_;
  std::vector<std::vector<double>> rate_cdf_;

  // coherent
  std::vector<std::unique_ptr<CoherentSampler>> csampler_;
  std::vector<double> coh_cdf_;     ///< over accepted cells, normalized
  double sigma_coh_pb_ = 0.0;

  // rc (rc.hpp) -- built on EVERY channel, so it lives outside the tagged and
  // coherent blocks above.  Null when cfg_.rc == Off.
  std::shared_ptr<RcModel> rc_;

  // beams / bookkeeping
  Vec4 beam_e_, beam_ion_;
  int ion_pdg_ = 0;
  double ion_mass_ = 0.0, ion_charge_ = 0.0;
  std::vector<double> sigma_, lumi_;
  std::vector<std::uint64_t> counts_, offset_;
  std::uint64_t total_ = 0;
  std::uint64_t cursor_ = 0;
};

}  // namespace lipolgen

#endif  // LIPOLGEN_PIPELINE_HPP
