// SPDX-License-Identifier: GPL-3.0-or-later
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

#include <cmath>
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
  /// The UNPOLARISED and POLARISED structure-function backends of the struck
  /// cluster's own kernel -- `InclusiveKernel::Options::f2_source` and
  /// `::g1_model`, forwarded verbatim by `struck_cluster_kernel`.
  ///
  /// BOTH NULL BY DEFAULT, which is the kernel class's own null branch
  /// (`src/core/xsec.cpp`): `ToyF2`, and `ToyG1` built on that same `ToyF2`
  /// with the kernel's own `r_func`.  So a tagged run that sets neither is
  /// bit for bit what it was before these fields existed -- the default is
  /// identical BY CONSTRUCTION, not by inspection.
  ///
  /// These are the tagged channels' ONLY structure-function injection point.
  /// `PipelineConfig::kernel` does not reach them: `Pipeline` reads it on the
  /// non-tagged branch alone, so a kernel attached to a tagged config is
  /// silently ignored.  `Pipeline` fills these two from
  /// `PipelineConfig::unpol_sf_obj` / `pol_sf_obj` (`--unpol-sf` /
  /// `--pol-sf`); a C++ caller driving `make_struck_cluster_source` directly
  /// sets them here.
  ///
  /// R is deliberately NOT here.  It is a third axis with its own measured
  /// justification still to be made (`r_sigma_lt` against `r1998` is up to
  /// 19 % in (1 + R) and 38 % of the shipped window is outside R1998's own
  /// support), and shipping it in the same change as these two would make
  /// two independent 10-40 % movements unattributable.
  std::shared_ptr<const UnpolSF> f2_source;
  std::shared_ptr<const PolSF> g1_model;
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
/// The run-surface name of a triton spectral function -- the two
/// strings `--triton-sf` takes, defined ONCE (the `cluster_wave_name`
/// note in cluster.hpp).
const char* triton_sf_name(TritonSfChoice s);

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
///   `Li6Convolution` PASSED its A = 2 validation gate on 2026-09-03 -- for
///   ONE unpolarised nucleon input, not for the shipped default.  G3b peak
///   ratio 0.843 with CDKS's own MSTW2008 LO at their Eq. (21)
///   delta-function; 0.440 on the DEFAULT `ToyF2`, which is OUTSIDE the
///   gate's [0.5, 2] window.  So the lift covers numbers made with
///   `b1_unpol = B1UnpolSource::Mstw` (CLI `--b1-unpol mstw`) and NOT the
///   default toy backend's.  With the CD-Bonn wave function as well the ratio
///   is 1.000338 -- read that as "the residual is now below the error of
///   digitizing a published figure", never as three-digit agreement with
///   CDKS, and note that it is specific to CD-Bonn AND MSTW together and is
///   not the default.  See b1_nuclear.hpp's header block for the four
///   conditions attached to that and
///   docs/open_items/run_2026-09-03/phase_A_numbers.md for the numbers.
///   That gate is A = 2 and tests NOTHING about the alpha-d step, so this
///   stays opt-in and band-only: never quote a single row, always
///   {0, 1, 2} x b1.
enum class B1Model : std::uint8_t { Miller, Cdks, Li6Convolution };
const char* b1_model_name(B1Model m);

/// Which UNPOLARISED structure-function backend the `Li6Convolution` b1
/// backend folds its own F1 against -- `Li6ConvolutionOptions::unpol`
/// (`PipelineConfig::b1_unpol`; CLI `--b1-unpol`).
///
/// WHY THIS EXISTS.  The A = 2 gate of 2026-09-03 passes at G3b = 0.843
/// with CDKS's OWN MSTW2008 LO and 0.440 on the library `ToyF2`, i.e. the
/// pass is a statement about a configuration, not about a build.  Until this
/// selector existed the generator could not EMIT that configuration:
/// `default_inclusive_kernel` hard-wired `ToyF2` into
/// `Li6ConvolutionOptions::unpol`, and the only other route -- hand-building
/// a kernel and handing it to `PipelineConfig::kernel` -- is refused
/// together with a non-Miller `b1_model` by `validate()` (and stays
/// refused).  So every 6Li b1 number the run surface could produce was made
/// at G3b = 0.440.  The selector closes that gap: the gate-passing
/// backend is now reachable from `PipelineConfig`, the CLI and the Python
/// API, THROUGH the pipeline's own kernel construction.
///
///   Toy      the DEFAULT and bit for bit what every published number was
///            made with: the kernel's own `ToyF2`, shared as ONE object with
///            `InclusiveKernel::Options::f2_source` (see
///            `default_inclusive_kernel`).
///   Mstw     `MstwSF` -- MSTW2008 LO central over PYTHIA 8's own
///            `pdfdata/mstw2008lo.00.dat` (`mstw_sf.hpp`), which is the PDF
///            CDKS computed their b1_d with.  MEASURED on the shipped
///            observable `Li6ConvolutionB1::b1(x, 2.5)`: x1.847766 at
///            x = 0.10, x1.275961 at 0.30, x0.816971 at 0.50 -- up to a
///            factor 1.85, and NOT monotone, so it is not a normalisation.
///   Ct18Nlo  `LhapdfSF("CT18NLO", 0)` -- the stand-in the phase-D numbers
///            were made with, kept selectable so the systematic can be
///            quoted rather than remembered.  x2.221302 / x1.238833 /
///            x1.045870 at the same three points.
///   Custom   whatever object the caller put in `b1_unpol_sf`.  The
///            `OpticsChoice::Custom` arrangement exactly: the enum is the
///            PROVENANCE (it is what `meta["b1_unpol"]` records) and the
///            object is the realisation.
///
/// THE CORE LIBRARY CANNOT BUILD `Mstw` OR `Ct18Nlo` ITSELF, by design:
/// `MstwSF` lives in the optional PYTHIA tier and `LhapdfSF` in the optional
/// LHAPDF tier, and `sf.hpp`'s rule is that the core links neither.  So the
/// enum names the backend and `b1_unpol_sf` carries it, filled by the layer
/// that has the tier -- `_lipolgen.set_b1_unpol(cfg, source)`, which is the
/// `set_pythia_hadronizer` arrangement.  `validate()` REFUSES a named
/// backend with an empty slot and says which tier is missing; it is NEVER
/// silently replaced by `ToyF2`.
///
/// SCOPE, and the price of it.  This reaches `Li6ConvolutionOptions::unpol`
/// and NOTHING else.  `InclusiveKernel::Options::f2_source` -- the F1 of the
/// unpolarised rate, and so the D_phi denominator of the tensor weight -- is
/// set by a SEPARATE flag, `--unpol-sf` (`UnpolSfSource` below).  AT ITS
/// SHIPPED DEFAULT `toy` that F1 is `ToyF2`'s on every `--b1-unpol` setting,
/// which is what keeps the SPIN-BLIND cell cross section
/// (`InclusiveSampler::cell_xsec_pb`) bit for bit under THIS flag -- only
/// the per-state tensor shift moves, and with it the tensor-weighted
/// per-category cross sections, which is the point.  (With `--unpol-sf`
/// set to anything else the rate moves too, but it moves because of THAT
/// flag; the sentence above is about this one.)  The consequence is that
/// with anything but
/// `Toy` the numerator's F1 and the denominator's F1 are no longer the same
/// object, so the partial cancellation `b1_nuclear.hpp` relies on for its
/// `r_func` default is gone: `--b1-unpol mstw` is a statement about b1, not
/// about the rate.  Say which backend a plot used.
///
/// WHY THIS IS A SEPARATE FLAG FROM `--unpol-sf`, and what happens when both
/// are set.  They are different quantities with different conventions: this
/// one is the DEUTERON's per-nucleon F1 inside CDKS Eq. (22), through
/// `f1_cdks`, which is explicitly NOT `NuclearF2::f1a`
/// (`b1_nuclear.hpp`), while `--unpol-sf` is the whole-nucleus F2A/F1A of
/// the rate.  And this one is a CDKS-COMPARABILITY choice that is
/// load-bearing for a published verdict -- the A = 2 gate passes at
/// G3b = 0.843243 with MSTW2008 LO and fails at 0.440 on the toy -- so a
/// user who wants a realistic RATE (`--unpol-sf ct18nlo`) must not be forced
/// to move the b1 gate row off MSTW as a side effect, nor the other way
/// round.  Two flags, therefore.  ONE refusal closes the ambiguity they
/// create: `default_inclusive_kernel` folds the convolution against the
/// KERNEL's own `f2_source` when `b1_unpol` is null, so under
/// `--unpol-sf ct18nlo` a `b1_unpol = Toy` would quietly mean "CT18NLO"
/// while `meta["b1_unpol"]` still wrote "toy".  `validate()` throws on
/// exactly that combination and asks for the b1 backend to be named.
enum class B1UnpolSource : std::uint8_t { Toy, Mstw, Ct18Nlo, Custom };
const char* b1_unpol_name(B1UnpolSource s);

/// WHICH R = sigma_L/sigma_T the 6Li alpha-d convolution and the kernel that
/// carries its denominator are BOTH built with -- ONE hook, threaded through
/// `default_inclusive_kernel` into `Li6ConvolutionOptions::r_func` AND
/// `InclusiveKernel::Options::r_func` (`PipelineConfig::r_source`; CLI
/// `--r-source`).  OFF (`Unset`) by default, and then today bit for bit.
///
/// WHY THIS EXISTS.  The shipped tensor observable is a RATIO: A_zz is
/// -(2/3) K/D_phi with K = b1 + (1-y)/(x y^2) b2 and
/// D_phi = F1 + (1-y)/(x y^2) F2 (`asymmetries.hpp`, `azz`).  Until this
/// selector existed the two halves of that ratio could not be made to agree
/// on R.  The numerator's R is `Li6ConvolutionOptions::r_func`
/// (`b1_nuclear.hpp`), whose null default is `r_sigma_lt`; the denominator's
/// is the kernel's own `Options::r_func`, whose null default is also
/// `r_sigma_lt` -- so they agree TODAY, by coincidence of two independent
/// defaults, and the only way to move one was to move it ALONE.  That is the
/// decision `STATUS.md` row 3 records: the A = 2 gate
/// (`DeuteronConvolutionB1::Options::r_func`) defaults to `r1998`, CDKS's
/// SLAC world fit, and adopting it in the numerator alone does not cancel
/// against a denominator still on `r_sigma_lt`.  This selector is the
/// registry's option (iii), "the better third option": whichever R is named,
/// numerator and denominator are built from the SAME object.
///
///   Unset     the DEFAULT.  NOTHING is threaded: both hooks stay null and
///             each consumer resolves its own `r_sigma_lt` through
///             `resolve_r` (`sf.hpp`), exactly as before this selector
///             existed.  Bit for bit BY CONSTRUCTION -- the branch that
///             sets the hooks is not taken.
///   SigmaLt   ONE `r_sigma_lt` object into both hooks.  MEASURED
///             bit-identical to `Unset` on the tensor observables at every
///             standard point (docs/open_items/run_2026-09-06/
///             phase_A_numbers.md sec. A1) -- which is the wiring's own
///             test, not a physics claim: it says the shared hook reaches
///             both halves and changes nothing when it names the value they
///             already had.
///   R1998     ONE `r1998` object into both -- the A = 2 gate's R, on both
///             halves of the ratio.  This is the value that costs something,
///             and sec. A1 is what it costs.
///
/// SCOPE, and why it is REFUSED elsewhere rather than labelled.  The hook is
/// installed by the `Li6Convolution` branch of `default_inclusive_kernel` and
/// nowhere else, so `validate()` refuses anything but `Unset` off
/// `b1_model = li6-convolution` -- which is already inclusive-only and 6Li-
/// only and refused beside a caller-supplied `kernel`.  The `b1_unpol` rule
/// verbatim, for the `b1_unpol` reason: `meta["r_source"]` is written
/// unconditionally, so a value that did not reach the rate may not be
/// recorded as if it had.  On `Miller` and `Cdks` the numerator has no F1 of
/// its own to give an R to (Miller is a ratio model, Cdks a digitized
/// column), so "one shared R" is not a statement those branches can make.
///
/// WHAT IT MOVES, and the asymmetry with `--b1-unpol`.  `--b1-unpol` moves
/// the NUMERATOR alone and leaves the spin-blind cell cross section bit for
/// bit; this one moves BOTH -- the kernel's `r_func` is threaded to all four
/// places the kernel needs R (F1 in `NuclearF2::f1a`, F_L in
/// `dsigma_dx_dq2`, D(y) in `depolarization_d`, and the default `ToyG1`'s
/// F1), so the unpolarised rate moves too.  That is the POINT: the ratio is
/// what is quoted, and its two halves are now built from one R.  Say which
/// R a plot used.
enum class RSource : std::uint8_t { Unset, SigmaLt, R1998 };
const char* r_source_name(RSource s);

/// Which UNPOLARISED structure-function backend supplies F2 -- and through
/// it F1, F_L and the whole unpolarised rate -- to EVERY kernel the
/// `Pipeline` builds (`PipelineConfig::unpol_sf`; CLI `--unpol-sf`).
///
/// WHY THIS EXISTS.  Until 2026-09-04 the run surface could not select the
/// RATE's unpolarised backend at all.  (`--b1-unpol` / `B1UnpolSource` above
/// landed the same day and is NOT a counter-example: it reaches
/// `Li6ConvolutionOptions::unpol` -- the b1 convolution's own F1 -- and
/// nothing else, which is exactly why the two are separate flags.)
/// `default_inclusive_kernel` hard-wired
/// `ToyF2` into the inclusive (and so the coherent) kernel, and
/// `struck_cluster_kernel` had NO structure-function slot of any kind, so
/// every tagged run was `ToyF2` / `ToyG1` with no way to say otherwise --
/// not even the hand-built `PipelineConfig::kernel` escape hatch, which the
/// `Pipeline` reads on the non-tagged branch alone.  `ToyF2` is labelled
/// TOY in `sf.hpp` and is anchored BY EYE ("adequate for phase-space maps
/// and factor-1.5 rate estimates ONLY"); this selector is what lets a run
/// be made on a real fit instead, and lets the difference be quoted.
///
///   Toy      the DEFAULT and bit for bit what every published number was
///            made with: the kernel class's own `ToyF2` (the pipeline names
///            it explicitly so the `Li6Convolution` b1 branch can share the
///            SAME object).  Leaves both `Options` slots exactly as they
///            were, so the default path is unchanged BY CONSTRUCTION.
///   Mstw     `MstwSF` -- MSTW2008 LO central over PYTHIA 8's own
///            `pdfdata/mstw2008lo.00.dat` (`mstw_sf.hpp`).  Needs the
///            OPTIONAL PYTHIA tier.
///   Ct18Nlo  `LhapdfSF("CT18NLO", 0)`.  Needs the OPTIONAL LHAPDF tier and
///            the CT18NLO set on disk.  READ THE GRID CLAUSE BELOW.
///   Custom   whatever object the caller put in `unpol_sf_obj`.  The
///            `OpticsChoice::Custom` arrangement: the enum is the
///            PROVENANCE (it is what `meta["unpol_sf"]` records) and the
///            object is the realisation.
///
/// MEASURED, 6Li inclusive at `default_configs("6Li")[1]`, the shipped
/// window: the sum of `InclusiveSampler::cell_xsec_pb()` over the 3051
/// accepted cells is 591 846.2 pb on `toy`, 472 571.9 on `ct18nlo`
/// (x0.7985) and 469 556.5 on `mstw` (x0.7934) -- a 20 % rate change.  It is
/// NOT a normalisation: F2p at Q2 = 10 moves by x0.92 / x1.11 / x1.37 /
/// x1.26 at x = 0.01 / 0.10 / 0.30 / 0.50 (ct18nlo) and x0.82 / x1.01 /
/// x1.36 / x1.45 (mstw).  The toy's n/p ratio is the straight line
/// clip(1 - 0.75x, 0.25, 1) and is 0.9625 / 0.8500 / 0.6250 at
/// x = 0.05 / 0.20 / 0.50 against CT18NLO's 0.9218 / 0.7219 / 0.5035, i.e.
/// up to 24 % high -- and that ratio is what the T1 and T2 struck-nucleon
/// SPECIES draws use, not only the rate.
///
/// SCOPE -- it reaches EVERY kernel, which is the whole point.  The
/// inclusive kernel (`default_inclusive_kernel`), the coherent channel
/// (which owns no structure function and rides the inclusive cell cross
/// sections), and the tagged struck-cluster kernel
/// (`StruckClusterOptions::f2_source`).  Three consumers inherit the same
/// ONE object for free and must not grow a knob of their own
/// (docs/CONVENTIONS.md, "no physics number is defined twice"): the T1
/// breakup's species draw (`BreakupOptions::f2`, filled from the sampler's
/// kernel), `InclusiveGenerator::proton_fraction`, and
/// `InclusiveSampler`'s own cell cross sections.  The T2 PYTHIA bridge is
/// the ONE consumer outside the `Pipeline`'s reach -- `python/lipolgen/
/// cli.py` hands it the same backend explicitly when it builds one.
///
/// IT IS NOT ORTHOGONAL TO `PolSfSource`, and saying so is part of using
/// it.  `InclusiveKernel` builds its default `ToyG1` on its OWN base
/// `UnpolSF`, so `--unpol-sf` alone moves g1 as well: at `pol_sf = Toy`,
/// `--unpol-sf ct18nlo` moves 6Li F1 by x0.99 / x1.07 / x1.20 / x1.00 and
/// g1 by x1.05 / x1.13 / x1.31 / x1.06 at x = 0.05 / 0.10 / 0.30 / 0.50, so
/// A1 = g1/F1 moves 5-8 %.  `--pol-sf toy` therefore does NOT mean "g1
/// unchanged".  State BOTH settings next to any plot.
///
/// THE GRID CLAUSE (`ct18nlo`).  CT18NLO's grid starts at Q = 1.295 GeV,
/// i.e. Q2 = 1.677, and the shipped 6Li window's accepted cells start at
/// Q2 = 1.054.  Below its grid LHAPDF does not freeze -- it continues the
/// evolution downward and F2p falls fast (at x = 3e-4 it is a factor 2.34
/// below the toy at Q2 = 0.7).  MEASURED on the shipped window: 42.33 % of
/// the toy run's accepted cell cross section, and 36.18 % of a CT18NLO
/// run's OWN accepted cell cross section, sits below that grid floor.  This
/// is NOT refused -- refusing would make the flag unusable on the shipped
/// scenario -- but it is recorded in `meta["unpol_sf_below_grid_frac"]` and
/// printed at the run banner on every run, because a silent 36 % of the
/// rate on an off-grid extrapolation is exactly the kind of thing a reader
/// must not have to rediscover.  THE FRACTION IS PER CHANNEL, taken against
/// `Pipeline::cell_rate_weights_pb()` and not against the inclusive cells:
/// on the COHERENT channel the same run's own rate is `sigma_cell *
/// f_coh(x)`, and 44.746 % of it is below that floor (measured 2026-09-05,
/// `--channel coherent --unpol-sf ct18nlo` at config 1).  There is NO `--q2-min` flag; raising
/// `Scenario::q2_min` (Python) to 1.7 takes the fraction to 0 and the
/// shipped 6Li ct18nlo cross section from 472571.9 to 301599.8 pb.
///
/// THE CORE LIBRARY CANNOT BUILD `Mstw` OR `Ct18Nlo` ITSELF, by design --
/// the same rule and the same arrangement as `B1UnpolSource`: `MstwSF` is
/// in the optional PYTHIA tier, `LhapdfSF` in the optional LHAPDF tier, and
/// `sf.hpp`'s rule is that the core links neither.  The object is built one
/// layer up by `_lipolgen.set_unpol_sf(cfg, source)`, and `validate()`
/// REFUSES a named backend with an empty slot, naming the missing tier.  It
/// is NEVER silently replaced by `ToyF2`.
enum class UnpolSfSource : std::uint8_t { Toy, Mstw, Ct18Nlo, Custom };
const char* unpol_sf_name(UnpolSfSource s);

/// Which POLARISED structure-function backend supplies g1 -- and through the
/// Wandzura-Wilczek relation g2 -- to every kernel the `Pipeline` builds
/// (`PipelineConfig::pol_sf`; CLI `--pol-sf`).  The `UnpolSfSource` block
/// above carries the shared rules; only what differs is repeated here.
///
///   Toy       the DEFAULT and bit for bit: the kernel's own `ToyG1`, built
///             on the kernel's own `UnpolSF` and `r_func`.
///   NnpdfPol  `LhapdfG1("NNPDFpol11_100", 0)`, which is already
///             `LhapdfG1`'s own declared default.  Needs the OPTIONAL
///             LHAPDF tier.  It is the ONLY polarised set installed in this
///             tree's LHAPDF store, so there is no second row to offer.
///   Custom    whatever object the caller put in `pol_sf_obj`.
///
/// MEASURED, 6Li per-nucleon g1A (and with it A_par, which tracks it to
/// better than 0.1 % at y = 0.5): nnpdfpol/toy = 0.669 / 0.676 / 1.114 /
/// 1.349 / 1.298 / 0.982 / 0.636 at (x, Q2) = (0.01, 2.5), (0.05, 5),
/// (0.10, 10), (0.20, 10), (0.30, 15), (0.50, 25), (0.70, 50).  So x0.64 to
/// x1.35 over the generator window, and NOT monotone.
///
/// THE NEUTRON IS WORSE THAN THE NUCLEUS, AND IT IS A SIGN.  `ToyG1`'s
/// a1n(x) = -0.07(1-x)^2 + 0.8 x^2.2 crosses zero near x ~ 0.25 and is
/// POSITIVE above it; NNPDFpol1.1's g1n stays negative to x ~ 0.6.
/// Measured at Q2 = 10, g1n toy vs NNPDFpol: -0.0801 / -0.1365 at x = 0.10,
/// -0.0114 / -0.0608 at 0.20, +0.00521 / -0.02740 at 0.30, +0.00779 /
/// -0.00037 at 0.50.  **The shipped toy g1n has the WRONG SIGN over roughly
/// 0.25 < x < 0.6.**  On isoscalar 6Li the proton term dominates and this
/// mostly hides; on a NEUTRON-tagged run (`d` + `p` tagging, whose
/// `dis_target` is `NEUTRON_TARGET`) it does not.  That is the strongest
/// single reason this selector exists.
///
/// SCOPE -- IT DOES NOT REACH THE COHERENT CHANNEL, and unlike
/// `UnpolSfSource` it never could.  It reaches the INCLUSIVE kernel
/// (`InclusiveKernel::Options::g1_model`) and the TAGGED struck-cluster
/// kernel (`StruckClusterOptions::g1_model`).  On `CoherentLi6` the rate is
/// spin-independent and no g1 is evaluated anywhere in the run, so the run
/// does not record a backend name there at all: `pol_sf_is_read` is false,
/// `meta["pol_sf"]` carries `pol_sf_unread_label` and
/// `meta["pol_sf_reach"]` the reason.  This is the one asymmetry between the
/// two selectors, and it is stated at every surface rather than smoothed
/// over -- the flag pair reaching "EVERY kernel" is true of `--unpol-sf`
/// alone.
///
/// g2 needs no selector of its own: `PolSF::g2p`/`g2n` is Wandzura-Wilczek
/// built on THIS backend's own g1, so one field moves g1 and g2 together.
/// The polarised-EMC `medium_ratio` argument of `PolSF::g1_nucleus` is
/// still passed by nobody, and this flag does not change that.
enum class PolSfSource : std::uint8_t { Toy, NnpdfPol, Custom };
const char* pol_sf_name(PolSfSource s);

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
  /// The ION-LEVEL kernel, hand-built by the caller.  Null (the default) =
  /// `default_inclusive_kernel`.
  ///
  /// READ ON `Inclusive` AND `CoherentLi6` ONLY -- those are the two channels
  /// whose sampler is the ion-level `InclusiveSampler`, and on the coherent
  /// one it does reach the rate, since the coherent yield rides those very
  /// cell cross sections.  On a TAGGED channel the sampler is the
  /// STRUCK-CLUSTER one (`make_struck_cluster_source`) and this field is not
  /// read at all, so `validate()` REFUSES it there: it used to be accepted
  /// in silence while `meta["unpol_sf"]`, `["pol_sf"]` and `["b1_model"]` all
  /// wrote "caller-supplied kernel" for a run made on `ToyF2`/`ToyG1`
  /// (measured 2026-09-05: a tagged-6Li-alpha config carrying an
  /// `InclusiveKernel` on CT18NLO ran at F2A(0.3, 10) = 0.3691149345, the toy
  /// value, and CT18NLO's 0.4639703442 appeared nowhere).  A tagged run
  /// injects structure functions through `struck.f2_source` /
  /// `struck.g1_model`, which the two selectors fill.
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
  /// Which family of radial forms the tagged channels use.  `Hulthen` is the
  /// default and keeps every published number bit-for-bit; `VmcAV18` swaps in
  /// the ANL VMC tables and then IGNORES `cluster_beta` and `p_d`
  /// (`docs/CONVENTIONS.md`).
  ///
  /// EXACTLY WHERE IT REACHES, AND WHERE IT DOES NOT.  Read this before
  /// quoting anything made with it.
  ///
  /// IT REACHES (all three tagged channels; `validate()` REFUSES it anywhere
  /// else, see below):
  ///   * the CLUSTER RELATIVE wave function -- alpha-d (6Li), alpha-t (7Li),
  ///     p-n (the d control) -- which is what the flag is named for;
  ///   * SINCE 2026-09-04 the deuteron control (open item C5.4).  The comment
  ///     here used to say the control "is always Hulthen -- there is no VMC
  ///     d -> p+n cluster table, the deuteron IS the cluster".  True of the
  ///     `momenta/` files and FALSE of `data/vmc/deuteron/fdeut.av18`, whose
  ///     u(k) and w(k) ARE the p-n relative S and D waves; the flag was being
  ///     SILENTLY IGNORED there, the same defect the b1 knobs are refused for.
  ///     It now selects the exact AV18 deuteron at its own P_D = 0.057600
  ///     instead of the analytic pair at the scenario `P_D_DEUTERON` = 0.045
  ///     (-2.03 % on that channel's vector dilution, -1.18 % on its tensor
  ///     one);
  ///   * SINCE 2026-09-04 the EMBEDDED DEUTERON of the 6Li alpha tag, in BOTH
  ///     places a run reads it (open item C5.5b): `TaggedChannel::dis_target`
  ///     -- the struck cluster's g1 -- is `DEUTERON_AV18()` (tagged.hpp), and
  ///     `BreakupOptions::source` -- the T1 struck-nucleon MOMENTUM draw and
  ///     its spin draw, both from that same deuteron (`sample_kc_one`) -- is
  ///     this field.  It used to reach NEITHER, so a `--cluster-wave vmc` run
  ///     took the alpha-d RELATIVE motion from the ANL VMC AV18+UX overlap and
  ///     the embedded deuteron from the 0.045 scenario: two deuteron
  ///     wave-function families in one run, with every polarized tagged-alpha
  ///     observable 2.069 % HIGH against the wave function this flag says it
  ///     selects.  Fixing it moved that opt-in path by -2.027 % and left the
  ///     Hulthen default bit for bit (T27).
  ///
  /// IT DOES NOT REACH:
  ///   * the 7Li TRITON's internal spin structure (`TRITON()`, a Faddeev-family
  ///     per-nucleon slot).  Not an oversight and not silent: this tree has no
  ///     AV18 A = 3 wave function to switch it to, so on 7Li the flag selects
  ///     the alpha-t relative motion and nothing else;
  ///   * the 6Li INCLUSIVE constants (`LI6_CLUSTER_POLARIZATION` = 0.811228,
  ///     `LI6_B1_RANK2_TRANSFER` = 0.921947), which are not read on a tagged
  ///     channel at all.  A `VmcAV18` tagged row and an inclusive row of one
  ///     PROGRAMME therefore describe 6Li with wave functions whose vector
  ///     dilutions differ by 11.61 % and whose rank-2 transfers differ by
  ///     6.58 % (C5.5) -- band any comparison of the two.  Since 2026-09-04
  ///     `validate()` REFUSES this field on `Inclusive` and `CoherentLi6`
  ///     rather than accepting it unread, because accepting it was the route
  ///     by which that 11.61 % reached someone who thought they had asked for
  ///     a VMC 6Li.
  ///
  /// WHAT ONE `VmcAV18` 6Li RUN THEREFORE SAYS THE WHOLE-NUCLEUS POLARIZATION
  /// IS: 0.887076 = `li6_cluster_polarization(VMC_P_D_LI6,
  /// deuteron_av18_p_d())`, one Hamiltonian end to end, against the shipped
  /// default's 0.811228.  Neither reproduces the ab initio six-body
  /// `LI6_POLARIZATION_VMC_SIX_BODY` = 0.848 (+4.6 % and -4.3 %); the cluster
  /// PRODUCT is what carries that error and C5.5 is the argument.
  ClusterWaveSource cluster_wave = ClusterWaveSource::Hulthen;
  /// Monte Carlo band of the ANL VMC tables, in units of their own printed
  /// 1-sigma error column, fully correlated across k (`li6_alpha_channel`,
  /// open item C5.2).  0 = today bit for bit.  `validate()` REFUSES it unless
  /// the run is a lithium alpha-tag channel on `cluster_wave = VmcAV18`: the
  /// Hulthen forms carry no MC error and `fdeut.av18` prints none, so anywhere
  /// else the knob would be recorded in the metadata without having run.
  /// Measured: +-1 sigma moves the tagged tensor dilution by 0.02 %, i.e. the
  /// ANL statistics are NOT the systematic that matters (the wave-function
  /// choice, worth 6.6 %, is).
  double cluster_vmc_mc_sigma = 0.0;
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
  /// |t| truncation [GeV^2].  0.2 since 2026-08-29 (P4).  The reason it is
  /// 0.2 is the ANCHOR RANGE -- the deformation input is digitized only to
  /// |t| = 0.30 -- and NOT positivity, which binds only at the shipped
  /// eps_b0 and is the contingent second reason
  /// (`COHERENT_T_MAX_DEFAULT`, `CoherentScenario::t_positivity_edge`).
  /// The `Pipeline` constructor still refuses a range in which the truncated
  /// azimuthal weight 1 + c2 cos 2(phi_t - phi_S) would go negative.
  /// Reachable as `--coherent-t-max` and recorded in the npz `meta` since
  /// 2026-09-04: moving it moves the whole |t| spectrum, the tag acceptance
  /// and every c_2 in the file, so a run that moved it must be
  /// distinguishable from one that did not.
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
  /// Which unpolarised backend the `Li6Convolution` b1 folds its F1 against
  /// (`B1UnpolSource` above; CLI `--b1-unpol`).  `Toy` is the DEFAULT and is
  /// today bit for bit.  Read ONLY by the `Li6Convolution` branch of
  /// `default_inclusive_kernel`, so `validate()` refuses anything but `Toy`
  /// on `Miller` and `Cdks` -- same provenance rule as the band and the
  /// alpha-d weight, and the same reason: `meta["b1_unpol"]` is written
  /// unconditionally.
  B1UnpolSource b1_unpol = B1UnpolSource::Toy;
  /// The object `b1_unpol` names.  Empty for `Toy` (the kernel's own `ToyF2`
  /// is used, as one shared object); REQUIRED for every other value, because
  /// the core library links neither the PYTHIA nor the LHAPDF tier and so
  /// cannot construct `MstwSF` / `LhapdfSF` here.  `validate()` refuses a
  /// named backend with an empty slot and names the missing tier -- it is
  /// never silently replaced by `ToyF2`.  Python fills both fields at once
  /// with `_lipolgen.set_b1_unpol(cfg, source)`; a C++ caller sets them
  /// together, the `optics_choice` / `optics` arrangement.
  std::shared_ptr<const UnpolSF> b1_unpol_sf;
  /// Which R = sigma_L/sigma_T is threaded into BOTH the `Li6Convolution`
  /// numerator and the kernel that carries the ratio's denominator
  /// (`RSource` above; CLI `--r-source`).  `Unset` is the DEFAULT and is
  /// today bit for bit -- it installs nothing, and both hooks keep the null
  /// that `resolve_r` turns into `r_sigma_lt`.  Read ONLY by the
  /// `Li6Convolution` branch of `default_inclusive_kernel`, so `validate()`
  /// refuses anything but `Unset` on `Miller` and `Cdks` -- the `b1_unpol`
  /// provenance rule verbatim, and `meta["r_source"]` is written
  /// unconditionally for the same reason.
  RSource r_source = RSource::Unset;

  // --- structure-function backends, every kernel ---------------------------
  /// Which UNPOLARISED backend EVERY kernel this config builds takes its F2
  /// from -- the inclusive kernel, the coherent channel that rides it, and
  /// the tagged struck-cluster kernel (`UnpolSfSource` above; CLI
  /// `--unpol-sf`).  `Toy` is the DEFAULT and is today bit for bit: it
  /// leaves both `Options` slots exactly as they were.  Unlike `b1_unpol`
  /// this one is NOT channel-restricted -- it reaches the rate on every
  /// channel, which is the point of it -- so `validate()` refuses it only
  /// together with a caller-supplied `kernel` (which would silently win) and
  /// in the `b1_unpol` combination below.
  UnpolSfSource unpol_sf = UnpolSfSource::Toy;
  /// The object `unpol_sf` names.  Empty for `Toy`; REQUIRED for every other
  /// value, because the core library links neither the PYTHIA nor the LHAPDF
  /// tier and so cannot construct `MstwSF` / `LhapdfSF` here.  `validate()`
  /// refuses a named backend with an empty slot and names the missing tier
  /// -- it is never silently replaced by `ToyF2`.  Python fills both fields
  /// at once with `_lipolgen.set_unpol_sf(cfg, source)`.
  std::shared_ptr<const UnpolSF> unpol_sf_obj;
  /// Which POLARISED backend the INCLUSIVE and TAGGED kernels take their g1
  /// (and, through Wandzura-Wilczek, their g2) from (`PolSfSource` above;
  /// CLI `--pol-sf`).  NOT every kernel: `PolSfSource`'s own SCOPE note is
  /// the whole rule, and it has TWO axes -- the coherent channel evaluates
  /// no g1, and neither does any unpolarised-beam plan, `tensor-thirds`
  /// included, which is the CLI's default.  The run LABELS the selector
  /// there rather than crediting it (`pol_sf_is_read(config, plan)`), and
  /// "every kernel the pipeline builds" is true of `unpol_sf` alone.
  /// `Toy` is the DEFAULT and is today bit for bit.
  PolSfSource pol_sf = PolSfSource::Toy;
  /// The object `pol_sf` names; the `unpol_sf_obj` rules verbatim, with
  /// `_lipolgen.set_pol_sf(cfg, source)` as the filler.
  std::shared_ptr<const PolSF> pol_sf_obj;

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
/// convolution's `unpol` when `b1_unpol` is null, so "the kernel's UnpolSF"
/// is then literally one object.
///
/// `b1_unpol` is `PipelineConfig::b1_unpol_sf`: the unpolarised backend the
/// `Li6Convolution` branch folds its own F1 against (`B1UnpolSource`).  NULL
/// -- the default, and every existing caller -- means the kernel's own
/// `ToyF2`, so this argument is bit for bit inert unless it is set.  It is
/// read by the `Li6Convolution` branch ONLY; `opt.f2_source`, and with it the
/// unpolarised rate, is `ToyF2` on every path.
///
/// `unpol_sf` and `pol_sf` are `PipelineConfig::unpol_sf_obj` /
/// `pol_sf_obj`: the backends of the kernel's OWN `f2_source` and
/// `g1_model` (`UnpolSfSource`, `PolSfSource`).  NULL -- the default, and
/// every caller that predates them -- leaves `f2_source` at the shared
/// `ToyF2` this function has always named and leaves `g1_model` UNSET, so
/// `InclusiveKernel` builds its own `ToyG1` on that same `ToyF2` exactly as
/// before.  Both are therefore bit for bit inert unless they are set.
/// Non-null `unpol_sf` also becomes what the `Li6Convolution` branch folds
/// against when `b1_unpol` is null -- which is why `PipelineConfig::
/// validate()` refuses `b1_unpol = Toy` there rather than let
/// `meta["b1_unpol"] = "toy"` name an object that is no longer `ToyF2`.
///
/// `r_source` is `PipelineConfig::r_source` (`RSource`): the ONE R hook this
/// function threads into `Li6ConvolutionOptions::r_func` AND
/// `InclusiveKernel::Options::r_func` together, so the tensor weight's
/// numerator and denominator are built from the same R.  `Unset` -- the
/// default, and every existing caller -- takes the branch that installs
/// NOTHING, leaving both hooks null exactly as they were, so it is bit for
/// bit inert by construction and not by an argument about `resolve_r`.  It
/// is read on the `Li6Convolution` branch ONLY (`Miller` is a ratio model
/// with no F1 of its own and `Cdks` a digitized column), which is what
/// `PipelineConfig::validate()` enforces.
std::shared_ptr<const InclusiveKernel> default_inclusive_kernel(
    const Ion& ion, B1Model model, double band_scale = 1.0,
    double w_alpha_d = 1.0,
    std::shared_ptr<const UnpolSF> b1_unpol = nullptr,
    std::shared_ptr<const UnpolSF> unpol_sf = nullptr,
    std::shared_ptr<const PolSF> pol_sf = nullptr,
    RSource r_source = RSource::Unset);

// ------------------------------------------- the 7Li rank-2 zero, out loud
//
// THE DEFECT THESE THREE CLOSE.  `InclusiveKernel::tables` (src/core/xsec.cpp)
// dispatches the rank-2 slots on the ION SPIN: spin 1 reads `b1_func /
// b2_func / delta_func`, spin 3/2 reads `b1_32_func / b2_32_func /
// delta_32_func` (`xsec.hpp`, `Options::b1_32_func`), an unset slot is 0.0 and
// an unset b2 is 2x*b1 and therefore 0 too.  `default_inclusive_kernel` fills
// the SPIN-1 slots only, so every 7Li inclusive run has had a tensor sector
// that is EXACTLY zero -- the tensor term of the phi-averaged rate, the
// cos 2phi (gluon transversity) amplitude and A_zz all of them -- while
// `meta["b1_model"]` recorded `"miller"`, a backend that did not run, and no
// run-surface line said anything at all.
//
// WHAT CHANGED IS THE LOUDNESS, NOT THE NUMBER.  The zero was there before and
// is there after: measured on both sides of this change, the two per-category
// cross sections of a 7Li |m| = 3/2 (T = +1) against |m| = 1/2 (T = -1) plan
// are the SAME DOUBLE, 590952.42641509 pb at the default 7Li config 1, and the
// asymmetry is exactly 0.0
// (docs/open_items/run_2026-09-03/phase_D_li7_rank2.md sec. 1.3).
//
// This is the repository's own rule -- A KNOB THAT DID NOT RUN MAY NOT BE
// RECORDED AS IF IT HAD, which `PipelineConfig::validate()` already enforces
// for `b1_band_scale` -- applied to the model name itself.
//
// WHY A BANNER AND A `meta` KEY RATHER THAN A REFUSAL: a 7Li inclusive run is
// a legitimate unpolarised and VECTOR-sector run (g1, A_parallel and the
// whole spin-1/2-like sector are correct at J = 3/2, and sec. 1.6 of that note
// measures the machinery downstream of the missing input as already right),
// so refusing it would remove working physics to protect one sector.  What the
// user must never get again is the SILENCE.  `python/lipolgen/__init__.py`'s
// `make_plan` does refuse, at the one place a refusal is the honest answer:
// the three spin-1-only tensor plans at J = 3/2, which used to build a spin-1
// plan and throw three frames down inside the sampler.
/// Does this configuration's INCLUSIVE kernel have a rank-2 (tensor) sector
/// that NOTHING fills?  True for exactly one shipped case: `--isotope 7Li
/// --channel inclusive` with no caller-supplied `kernel`.
///
/// FALSE on a caller-supplied `kernel` (its slots are the caller's, and
/// `meta` says "caller-supplied kernel"), false off the inclusive channel
/// (the inclusive rank-2 slots do not reach those rates at all), and false at
/// spin 1/2 or 0, where there is no rank-2 sector to be empty in the first
/// place.  It is the guard the `meta` labels and the CLI run banner key on.
bool inclusive_rank2_is_empty(const PipelineConfig& cfg);

/// What filled -- or did not fill -- this run's rank-2 slots, as one
/// sentence.  ONE DEFINITION, read by `meta["rank2_input"]` and printed by
/// the CLI run banner, so the claim cannot drift between them (the
/// `unpol_sf_grid_report` arrangement; docs/CONVENTIONS.md).
///
/// The `plan` is read for one clause only: when the run's fill carries a
/// rank-2 moment (`RunPlan::pzz_true`, which is T and not P_zz at J = 3/2,
/// bookkeeping.hpp) the sentence says so, because that is the run that asked
/// for the sector that is not there.
std::string rank2_input_report(const PipelineConfig& cfg, const RunPlan& plan);

/// The label `meta["b1_model"]` and `meta["b1_unpol"]` carry instead of a
/// backend name when `inclusive_rank2_is_empty` holds.  A `b1_model` of
/// `"miller"` on such a run names a backend that did not run; the two knobs
/// `b1_band_scale` and `b1_alpha_d_dwave_weight` stay NUMERIC and stay at
/// their defaults, because `validate()`'s Miller branch already refuses any
/// other value, so they cannot record a variation that did not run.
inline const char* rank2_none_label() {
  return "none (spin 3/2: no rank-2 input)";
}

/// Does the POLARISED selector (`PipelineConfig::pol_sf`, CLI `--pol-sf`)
/// reach the rate of this config's channel at all?
///
/// FALSE ON `CoherentLi6` AND NOWHERE ELSE, and that is a physics statement,
/// not a wiring gap: the coherent yield is `f_coh(x)` times the UNPOLARIZED
/// inclusive cell cross sections (`Pipeline`'s coherent branch), and this
/// channel's tensor signal is the recoil azimuth's
/// `1 + c2 cos 2(phi_t - phi_S)`, which `CoherentSampler` owns.  Nothing on
/// that path evaluates g1 -- not the rate, not a column, and not `RcModel`,
/// whose `applies()` is false there so it never asks the sampler for a
/// polarized `StateTables`.  MEASURED 2026-09-05 on `--channel coherent
/// --events 400 --seed 11`: `sigma_pb`, `sigma_per_category_pb` and all 47
/// generated array columns are bit-identical (`array_equal`, `equal_nan`)
/// between `--pol-sf toy` and `--pol-sf nnpdfpol`, under helicity-flip,
/// tensor-thirds and transverse-tensor plans.
///
/// The kernel still CARRIES the selected `g1_model` on that channel (the
/// coherent branch builds the ordinary inclusive kernel), so
/// `dis_sampler().kernel().tables(x, q2).g1` does move -- what never happens
/// is that the run reads it.  This function is about the RUN.
///
/// IT IS TWO AXES, NOT ONE, AND THE SECOND ONE IS THE SHIPPED DEFAULT.  Until
/// 2026-09-05 this predicate keyed on the CHANNEL alone, which made it right
/// on the coherent channel and WRONG on every unpolarised-beam RUN PLAN.  g1
/// enters the rate through exactly one product -- `InclusiveKernel::
/// amplitudes` adds `lam_e * pe * (m/J) * cos(theta_S) * A_par` and nothing
/// else reads `PolSF` -- so a fill with `lam_e * P_e == 0` in every category
/// never evaluates it either.  `tensor_thirds_plan`, `transverse_tensor_plan`
/// and `tensor_flip_plan` build every category at `lam_e = 0, pe = 0`
/// (src/core/bookkeeping.cpp), and `helicity_flip_plan` at `--pe 0` is the
/// same; `tensor-thirds` is the CLI's own default plan.  MEASURED 2026-09-05,
/// 600 events seed 7, sha256 over all 47 columns plus sigma_pb and
/// sigma_per_category_pb: `--pol-sf nnpdfpol` is bit-identical to `toy` on
/// inclusive-6Li, inclusive-d, tagged-6Li-alpha and tagged-d-p under
/// tensor-thirds, and on inclusive-6Li under transverse-tensor, tensor-flip
/// and helicity-flip at `--pe 0` -- while `meta["pol_sf"]` wrote "nnpdfpol"
/// and the banner said it reached the kernel.  Both axes are now here, which
/// is why the plan is an argument.
///
/// WHY IT IS NOT A REFUSAL, on EITHER axis.  `--unpol-sf` reaches the
/// coherent rate (through exactly those cell cross sections: x0.7058 on the
/// shipped 6Li run), so refusing its partner alone would make
/// `--unpol-sf X --pol-sf Y` fail on one channel of a three-channel scan and
/// nowhere else; and on the plan axis a refusal would make `--pol-sf nnpdfpol`
/// fail at the CLI's own DEFAULT plan.  It is labelled instead, the
/// `rank2_input_report` pattern -- and the label is what `meta["pol_sf"]`
/// carries, so nothing records a backend that did not run.  This is the
/// LABEL half of the criterion stated once at `KnobProvenance` below.
bool pol_sf_is_read(const PipelineConfig& cfg, const RunPlan& plan);

/// What the polarised selector reached on this run, as one sentence.  ONE
/// DEFINITION, read by `meta["pol_sf_reach"]`, by the `KnobProvenance` row
/// for `pol_sf` and printed by the CLI run banner (the `rank2_input_report` /
/// `unpol_sf_grid_report` arrangement; docs/CONVENTIONS.md).
std::string pol_sf_reach_report(const PipelineConfig& cfg, const RunPlan& plan);

/// The label `meta["pol_sf"]` carries INSTEAD of a backend name on a run
/// where `pol_sf_is_read` is false -- the `rank2_none_label` rule applied to
/// the polarised selector: a run whose g1 was never evaluated may not record
/// "nnpdfpol", and may not record "toy" either, because `ToyG1` did not run
/// on it any more than `LhapdfG1` did.  `meta["pol_sf_reach"]` carries the
/// reason beside it.  Two shapes, one per axis: "not read on channel
/// coherent-6Li" (the channel) and "not read under this run's fill" (the
/// plan).
std::string pol_sf_unread_label(const PipelineConfig& cfg, const RunPlan& plan);

// ================================================== the knob-provenance table
//
// THE RULE, AND WHY IT NEEDED A MECHANISM.  A KNOB THAT DID NOT RUN MAY NOT BE
// RECORDED IN THE `meta` OR PRINTED IN THE BANNER AS IF IT HAD.  That rule was
// enforced knob by knob -- `validate()` refuses `b1_band_scale` on Miller and
// `rc_qe_tensor_scale` without the quasi-elastic tail; `rank2_input_report`
// labels the 7Li rank-2 zero; `pol_sf_is_read` labels the coherent channel --
// and each of those closed ONE cell of a two-dimensional table it never
// wrote down.  It was then broken five times in one run, each time on an axis
// the previous fix had not looked at (the channel, then the run plan, then the
// T2 tier, then the rc sub-knobs, then the knobs recorded NOWHERE), because
// nothing enumerated the knobs and nothing could be tested exhaustively.
//
// `Pipeline::knob_provenance` is that enumeration: ONE function that returns
// EVERY user-settable knob of a run with what the run did with it.  The npz
// `meta` block, the CLI banner block and `python/tests/test_knob_provenance.py`
// -- which rebuilds the measured (spec x knob) matrix (12 (isotope, channel,
// plan) specs x 76 knob variants = 592 cells, re-measured 2026-09-06; not the
// full channel x plan product -- USAGE.md sec. 7c) and asserts
// the table against the OUTPUT HASH -- all read this one table, so a knob
// cannot be reported without a status and a new knob cannot be added without
// one either.
//
// THE CRITERION, STATED ONCE.  A knob this run does not read is REFUSED when
// the value it would record names a VARIATION OF A PIECE THAT DID NOT RUN --
// a scale, a band edge or a shape on a term the run computes as identically
// 1 -- because such a value is a claim that a systematic was PRICED, and no
// label makes a priced systematic un-priced (the `rc_qe_tensor_scale` and
// `b1_band_scale` precedents).  It is LABELLED when it names a BACKEND, an
// AXIS or a MEMBER OF A FAMILY that a channel-, plan- or set-scan sets
// UNIFORMLY across runs, because refusing one cell of such a scan costs more
// than it buys and the label carries the whole truth anyway (the `--pol-sf`
// on coherent and the `--pzz-mode` under the three tensor plans precedents;
// `--pzz` under helicity-flip was the second of them until 2026-09-06, when
// `--pzz-mode typed` made that flag reachable there and the label moved onto
// the mode).  Either way the
// knob is WRITTEN, with its status and its reason: silence is never an option,
// and that is what closes the fourth class -- knobs accepted and recorded
// nowhere at all.

/// What THIS run did with one user-settable knob.
enum class KnobStatus : std::uint8_t {
  Read,     ///< the run consults it: some quantity it computes is a function
            ///< of this knob, so another value would in general give another
            ///< file
  NotRead,  ///< nothing the run computes is a function of it; `reason` says
            ///< why, and `meta_value()` writes that reason in place of the
            ///< bare value
  Refused   ///< the AXIS is closed on this run: `PipelineConfig::validate()`
            ///< (or the `Pipeline` constructor) throws on any value but the
            ///< one shown, so the bare value cannot mislead and is kept
};
const char* knob_status_name(KnobStatus s);

/// One row of a run's knob-provenance table.
struct KnobProvenance {
  std::string name;    ///< the `meta` key stem / library field ("pol_sf")
  std::string flag;    ///< the CLI switch ("--pol-sf"); empty = API only
  std::string value;   ///< this run's own value, as recorded
  KnobStatus status = KnobStatus::Read;
  /// One sentence about THIS run, never empty.  For `NotRead` it OPENS with
  /// "not read on channel X" / "not read by plan Y" / "not read at Z" and
  /// then gives the reason.
  std::string reason;
  /// The SHORT form of `reason` -- its scope clause alone, everything before
  /// the first ": " or " -- ".  It is what a `meta` key carries in place of a
  /// bare value, so `meta["pol_sf"]` on a coherent run still reads exactly
  /// "not read on channel coherent-6Li" (`pol_sf_unread_label`) while the
  /// whole sentence stays available in the `knob_provenance` block beside it.
  /// Derived from `reason` in ONE place (`Pipeline::knob_provenance`), never
  /// typed twice.  Empty for a `Read` row.
  std::string label;
  bool at_default = true;   ///< `value` is the shipped default

  /// What a `meta` key carries for this row: the bare value where the knob
  /// RAN or where its axis is refused (so the value is necessarily the
  /// default), and the LABEL where it did not -- the `pol_sf_unread_label`
  /// rule, applied to every knob at once.
  const std::string& meta_value() const {
    return status == KnobStatus::NotRead ? label : value;
  }
};

/// The facts a `Pipeline` cannot see, because they are resolved one tier up.
///
/// TWO OF THEM, and neither is an oversight.  The T2 bridge lives in the
/// OPTIONAL PYTHIA tier, which the core library deliberately does not link
/// (`sf.hpp`'s rule), so `PipelineConfig::hadronizer` is a type-erased
/// `std::function` and the core cannot ask it for `PDF:PomSet`.  And a
/// `RunPlan` records its MOMENTS, not which flags produced them -- only the
/// caller knows whether `--pzz` was read or whether the fill came from
/// `helicity_flip_plan`'s max-entropy ladder, which since 2026-09-06 is the
/// `--pzz-mode` switch and is carried here as `pzz_mode`.
/// Default-constructed = "no T2
/// tier bound, plan provenance not supplied", which is what a C++ caller
/// driving `Pipeline` directly has; the rows that need what is missing are
/// then simply not emitted, never guessed.
struct KnobRunContext {
  /// The `--plan` name ("tensor-thirds", "helicity-flip", ...).  Empty = a
  /// caller-supplied `RunPlan`, and the `pz` / `pzz` / `rel_lumi_offset` rows
  /// are omitted (their reach is a property of the plan FACTORY, and `pe`
  /// needs no name: `helicity_flip_plan` is the only factory that produces
  /// `lam_e != 0`).
  std::string plan_name;
  /// What the caller ASKED the plan factory for.  Not the same thing as the
  /// plan's own `pz_true()` / `pzz_true()` / `pe_true()`, and the difference
  /// IS the defect: `transverse_tensor_plan(pzz, phi_s)` never sees `--pz`,
  /// `helicity_flip_plan` sees `--pzz` only at `--pzz-mode typed`, and the
  /// three tensor plans never see `--pe` -- so the fill's moments cannot say
  /// what was typed.
  /// NaN = not supplied; the row then falls back to the plan's own moment.
  double pz = std::nan("");
  double pzz = std::nan("");
  double pe = std::nan("");
  double rel_lumi_offset = std::nan("");
  /// WHICH FILL the caller asked `helicity_flip_plan` for: "ladder" (the
  /// default max-entropy branch, `HelicityFlipOptions::use_explicit_pzz =
  /// false`) or "typed" (`= true`, honouring `--pzz`).  A `RunPlan` cannot
  /// answer this either -- at J = 1 the two branches can in principle reach
  /// the same populations, and at J = 3/2 they differ in R_3, which the
  /// recorded moments do not carry -- so it is the caller's, exactly like
  /// `pz` / `pzz` above.  Empty = not supplied, and the `pzz_mode` row is
  /// then omitted rather than guessed.
  std::string pzz_mode;
  bool t2_bound = false;    ///< a hadronizer is attached to the config
  /// The attached bridge's own options.  The defaults are
  /// `PythiaBridgeOptions`' own, so a context that names no bridge still
  /// records the values a `--hadronize` run would have used rather than a 0
  /// that was never anything.
  bool t2_pomeron = true;   ///< `PythiaBridgeOptions::coherent_t2 == Pomeron`
  int pom_set = 6;          ///< `PDF:PomSet` (H1 2006 Fit B LO)
  double pom_rescale = 1.0; ///< `PDF:PomRescale`
};

/// Do any of this plan's categories carry a POLARISED BEAM -- some category
/// with `lam_e != 0` and `pe != 0`, and a fill with a non-zero population at
/// some m != 0?  That product, and only that product, is what
/// `InclusiveKernel::amplitudes` multiplies g1 by, so it is the exact
/// condition under which a run evaluates a polarised structure function.
/// ONE DEFINITION, read by `pol_sf_is_read` and by the `pe` row of
/// `Pipeline::knob_provenance`.
bool plan_has_beam_helicity(const RunPlan& plan);

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

  /// EVERY user-settable knob of this run, with what the run did with it --
  /// THE mechanism the rule "a knob that did not run may not be recorded as
  /// if it had" is enforced by (`KnobProvenance` above states the criterion
  /// once).  The npz `meta["knob_provenance"]` block, the CLI banner's
  /// provenance block and `python/tests/test_knob_provenance.py` all read
  /// THIS table and nothing else, so a knob cannot be reported without a
  /// status, and adding a knob without a row fails that test.
  ///
  /// IT IS A METHOD ON THE RUN, not on the config, because four of the
  /// rules are properties of the RUN and not of the configuration: the rc
  /// band and tail reaches are `RcModel::applies()` / `tail_applies()` on the
  /// model this run actually built; the tagged breakup rows key on
  /// `tier()`, which is `T0` on channels with no struck cluster whatever the
  /// config says; `x_max`'s reach on the coherent channel is decided against
  /// the upper accepted-x edge of the cells that carry coherent rate, on the
  /// UNCLIPPED cell set (`coherent_rate_x_edge` -- measuring it on the run's
  /// own already-clipped cells is self-referential and said NOT READ of an
  /// `x_max` that had just removed 320 cells); and the three ROUTE knobs are
  /// measured by re-routing this run's own events at the other envelopes
  /// (`RouteReach`).  Reading any of them off the config would be a second
  /// encoding, and in two of the four cases a wrong one.
  ///
  /// `ctx` carries what the core cannot see (`KnobRunContext`).  Rows the
  /// context cannot decide are omitted rather than guessed.
  std::vector<KnobProvenance> knob_provenance(
      const KnobRunContext& ctx = KnobRunContext()) const;
  /// Null outside `CoherentLi6`.  One sampler per category (they differ by
  /// P_zz and phi_S); index is the category index.
  const CoherentSampler* coherent_sampler(std::size_t category) const;

  /// THE PER-CELL WEIGHTS OF THIS RUN'S OWN ACCEPTED RATE, in the cell order
  /// of `dis_sampler().cell_xsec_pb()` -- the denominator any "how much of
  /// this run sits in region R" fraction must be taken against.
  ///
  /// It is NOT always `dis_sampler().cell_xsec_pb()`, which is why it exists.
  /// Per channel:
  ///
  ///   Inclusive  the sampler's own accepted cell cross sections, i.e. the
  ///              same vector (no copy).
  ///   tagged     the STRUCK-CLUSTER sampler's accepted cell cross sections.
  ///              The tagged rate is `sum_(M, m_S) pop_M p(m_S) sigma(m_S)`
  ///              and every one of those `sigma(m_S)` is a sum over THESE
  ///              cells, so the cell-to-cell shape of the rate is this
  ///              vector; what the sum adds on top is a per-(M, m_S) factor,
  ///              cell-independent, which cancels out of any fraction.
  ///   Coherent   `sigma_cell * f_coh(x)` over the cells that admit a
  ///              diffractive mass of at least `CoherentXpomModel::m_x_min`,
  ///              and ZERO on the ones that do not -- i.e. exactly the
  ///              weights `coh_cdf_` is built from, whose sum is
  ///              `sigma_pb()`.  `f_coh` falls by a factor ~26 across the
  ///              window, so this is a LARGE reweighting of the inclusive
  ///              cells: measured on the shipped 6Li ct18nlo run
  ///              (2026-09-05), 44.746 % of the coherent rate is below
  ///              CT18NLO's grid floor against 36.179 % of the inclusive
  ///              cells the coherent channel rides on.  Reporting the
  ///              inclusive number on a coherent run is the defect this
  ///              accessor closes (`unpol_sf_grid_report`).
  ///
  /// The RESIDUAL approximation, stated because it is not zero: the spin
  /// modulation `(1 + w_avg(m))` is per (cell, spin state) and is NOT in
  /// these weights on any channel -- they are the SPIN-BLIND rate, which is
  /// the only per-cell weight a run has that does not depend on which
  /// category is being asked about.
  const std::vector<double>& cell_rate_weights_pb() const;

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
  /// WHAT THE THREE ROUTE KNOBS DO TO **THIS RUN'S OWN** ROUTE LABELS.
  ///
  /// `--optics`, `n_sigma` and `pot_config` reach exactly one quantity, the
  /// per-event `route` column, and they reach it through exactly two branches
  /// of `route_charged`: the envelope is consulted only for a NEAR-BEAM
  /// fragment (|R - 1| < 0.05, theta < THETA_RP_OUTER) and `pot_config` only
  /// for an OVER-RIGID one (R > 1.05).  Whether a different value moves a
  /// LABEL is therefore a property of where THIS run's fragments fall, and it
  /// is not decidable from the channel.  Measured 2026-09-05: `pot_config`
  /// moves nothing on `coherent` at 60, 400, 2000 or 20 000 events (the
  /// intact recoil is never over-rigid) while it moves labels on
  /// `tagged-6Li-alpha` at every one of those sizes, 60 included; and
  /// `--optics yr-high-divergence` moves ONE `tagged-d-p` label at 60 events
  /// with SEED 7 and none at 2000 with SEED 11 (phase D's own pair,
  /// `phase_D_numbers.md` D6).  The seed is half of that statement: at the
  /// SHIPPED DEFAULT seed the same channel is not-read at 60, 400 and 2000
  /// under all four plans -- the one exception being `n_sigma` under
  /// `tensor-flip` at 2000, which is READ -- and at seed 1234 it is not-read
  /// at 60 and 2000 for all three configs and both YR envelopes.  So it is a
  /// property of the SAMPLE (seed and size together), not of the size and not
  /// of the channel, and no tabulated rule can stand in for the probe.
  ///
  /// So it is MEASURED, on this run, the way every script in this tree prices
  /// a sample at a second envelope: the sample is drawn ONCE and re-routed
  /// (`route_of` is a pure function of the finished record and the envelope,
  /// which is what the note above `route_of` is for).  The probe stops early
  /// as soon as all three knobs have been seen to move.
  struct RouteReach {
    bool has_route = false;      ///< the channel writes a far-forward fragment
    std::uint64_t n_probed = 0;  ///< events re-routed
    /// Per knob: did SOME other value move SOME label, which alternative was
    /// the first to do it, and how many labels it moved.
    bool optics_moves = false, n_sigma_moves = false, pot_moves = false;
    std::string optics_alt, n_sigma_alt, pot_alt;
    std::uint64_t optics_n = 0, n_sigma_n = 0, pot_n = 0;
    /// The alternatives actually tried, named, so the reason can say what
    /// "did not move" was measured against.
    std::string optics_tried, n_sigma_tried, pot_tried;
  };
  /// Computed at most once per `Pipeline` (`knob_provenance` is called by the
  /// banner and again by the metadata writer).
  const RouteReach& route_reach() const;
  /// Event `index` WITHOUT the rc weights and WITHOUT the T2 hadronizer.
  /// Neither moves a four-vector, and the probe must not run PYTHIA -- it is
  /// not re-entrant and its `PythiaBridgeStats` are recorded in `meta`.
  void route_probe_event(std::uint64_t index, Event& out) const;
  /// The upper accepted-x edge of the cells that carry coherent rate, on the
  /// UNCLIPPED cell set -- the window at the shipped `Scenario::x_max`, not
  /// this run's own already-clipped one.  Zero off `CoherentLi6`.
  double coherent_rate_x_edge() const;

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
  /// The UNNORMALIZED per-cell coherent rate the line above is the CDF of,
  /// kept because `cell_rate_weights_pb()` needs the weights and not their
  /// cumulative sum, and because differencing a normalized CDF back into
  /// weights loses the low-order bits of exactly the small cells a
  /// below-a-floor fraction is made of.  Empty off `CoherentLi6`.
  std::vector<double> coh_cell_pb_;
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

  // The route-knob probe, memoized.  `mutable` + `once_flag` and not a plain
  // cache, because `knob_provenance` is const and the class promises every
  // const entry point is thread-safe.
  mutable std::once_flag route_reach_once_;
  mutable RouteReach route_reach_;
};

}  // namespace lipolgen

#endif  // LIPOLGEN_PIPELINE_HPP
