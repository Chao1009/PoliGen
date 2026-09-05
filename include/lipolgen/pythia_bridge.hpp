// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef LIPOLGEN_PYTHIA_BRIDGE_HPP
#define LIPOLGEN_PYTHIA_BRIDGE_HPP

/// \file pythia_bridge.hpp
/// Tier T2: PYTHIA 8 showering / hadronization of the gamma*-nucleon system.
///
/// The core library never sees PYTHIA.  This header is PYTHIA-free as well
/// (pimpl), so including it costs nothing but `<memory>`; only
/// `src/pythia/*.cpp` needs the PYTHIA include path.
///
/// The physics and every approximation are written out in
/// docs/PYTHIA_BRIDGE.md.  The two facts that shape the whole design:
///
///   1. `Beams:allowMomentumSpread` runs the hard process at the *initial*
///      sqrt(s) (`ProcessLevel.cc:645` gates `newECM` on `doVarEcm`, which
///      `Pythia.cc:662` forbids with hard processes).  A Fermi-smeared
///      per-nucleon momentum swings sqrt(s) by +-20 %, so the only usable
///      route is an in-memory `Pythia8::LHAup` with `Beams:frameType = 5`
///      (docs/surveys/pythia8_survey.md sections 4(v) and 7).
///
///   2. With LHAup, PYTHIA still forces the total final-state four-momentum
///      to the *initialisation* beam total (`BeamRemnants.cc:662,935`:
///      `wPosRem = eCM - ...`).  It is therefore impossible to hand PYTHIA
///      the physical (k, p_N) of a Fermi-moving nucleon and get the physical
///      total back.  The bridge instead hands PYTHIA a *surrogate* e+N event
///      at the initialisation beams matched in (W^2, Q^2) -- i.e. the same
///      gamma*-nucleon subsystem -- and maps the resulting hadronic system
///      onto the physical one with a pure Lorentz transformation.  Momentum
///      and charge conservation are then exact by construction.
///
/// Because the hard process arrives through LHAup, **all `PhaseSpace:*`
/// settings are irrelevant here** -- including the two silent cuts
/// (`PhaseSpace:pTHatMinDiverge`, `PhaseSpace:mHatMin`) that
/// `PolarizedLithiumSim/tools/pythia8/gen_dis_hfs.py` has to set.  The hard
/// phase space is the core generator's; PYTHIA only showers and hadronizes.

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "lipolgen/beams.hpp"
#include "lipolgen/event.hpp"
#include "lipolgen/rng.hpp"
#include "lipolgen/sf.hpp"

namespace lipolgen {

// ---------------------------------------------------------------- helpers

/// Minkowski product with the (+,-,-,-) metric.
inline double dot4(const Vec4& a, const Vec4& b) {
  return a.e * b.e - a.px * b.px - a.py * b.py - a.pz * b.pz;
}
inline Vec4 scale4(double c, const Vec4& a) {
  return {c * a.e, c * a.px, c * a.py, c * a.pz};
}

/// The head-on-frame scattered electron of an inclusive DIS event.
///
/// Frame (docs/CONVENTIONS.md): ion along +z, electron along -z, so
///   k    = (E_e, 0, 0, -E_e)                       (electron treated massless)
///   p_N  = (E_N, 0, 0, +p_N)  for a nucleon at rest in the ion rest frame.
/// The per-nucleon DIS variables are the invariants
///   Q^2 = -q^2,  q = k - k',   y = (p_N.q)/(p_N.k),   x = Q^2/(2 p_N.q),
/// and `phi` is the azimuth of e' about +z, phi = atan2(k'_y, k'_x) in
/// [0, 2 pi).
///
/// Two linear conditions fix (E', k'_z):
///   k.k'  = Q^2/2            ->  E' + k'_z = Q^2/(2 E_e)
///   p_N.k'= (1-y) p_N.k      ->  E_N E' - p_Nz k'_z = (1-y) p_N.k
/// whose solution for a collinear p_N = (E_N, 0, 0, p) is
///   E'   = E_e (1-y) + p Q^2 / (2 E_e (E_N + p)),
///   k'_z = Q^2/(2 E_e) - E',
///   k'_T = sqrt(E'^2 - k'_z^2).
/// (The massless-target limit is the textbook E' = E_e(1-y) + x y E_p.)
///
/// `y` is derived from (x, Q^2) and the beams: y = Q^2 / (2 x p_N.k).
/// Returns false if the point is outside the physical region
/// (E' < |k'_z|, y outside (0, 1]).
bool dis_scattered_electron(double e_e, const Vec4& p_n, double x, double q2,
                            double phi, Vec4* k_out, double* y_out = nullptr);

/// The light-cone fraction xi of the struck parton that puts the outgoing
/// quark q + xi p_N on the mass shell m_out (zero by default):
///     (q + xi p_N)^2 = m_out^2
///   -> xi^2 p_N^2 + 2 xi (q.p_N) - (Q^2 + m_out^2) = 0,
/// solved in the numerically stable form
///     xi = (Q^2 + m_out^2)
///          / [ (q.p_N) + sqrt((q.p_N)^2 + p_N^2 (Q^2 + m_out^2)) ].
/// This is the root that reduces to xi = x_Bj for a massless target and a
/// massless outgoing quark, and stays the small positive root for p_N^2 < 0.
/// Returns false when no positive root below `xi_max` exists.
bool dis_parton_fraction(const Vec4& q, const Vec4& p_n, double* xi_out,
                         double xi_max = 1.0, double m_out = 0.0);

// -------------------------------------------------------------- summaries

/// Hadronic-final-state sums over `Role::Hadron` particles.
struct HfsSummary {
  double sigma_empz = 0;   ///< Sum (E - p_z)
  double px = 0, py = 0;   ///< Sum p_T, as a vector
  double pt = 0;           ///< |Sum p_T|
  double e = 0, pz = 0;    ///< Sum E, Sum p_z
  double charge = 0;       ///< Sum charge [e]
  int n_total = 0;         ///< multiplicity
  int n_charged = 0;
  int n_neutral = 0;
};

/// Sums over every `Role::Hadron`, `Status::Final` particle of the event.
HfsSummary hfs_summary(const Event& ev);

/// The truth identity the HFS is tested against:
///     Sum (E - p_z)_hadrons = 2 E_e y + m_N^2 / (E_N + p_Nz).
/// Exact up to Q^2 m_N^2 / (2 E_e (E_N + p_Nz)^2) -- 1.1e-5 GeV at
/// 10 x 99.5 GeV/u and Q^2 = 10 GeV^2 -- because the target's own
/// (E - p_z) = m_N^2/(E_N + p_Nz) is not zero.  Valid for a collinear
/// target; for a Fermi-moving one use `hfs_sigma_empz_exact`.
double hfs_sigma_empz_truth(double e_e, double y, const Vec4& p_n);

/// The exact statement, valid for any p_N:
///     Sum (E - p_z)_hadrons = (k + p_N - k')^- .
double hfs_sigma_empz_exact(const Vec4& k, const Vec4& p_n, const Vec4& k_out);

// ---------------------------------------------------------------- options

/// Which nucleon of the target is struck when the event does not name one.
enum class NucleonChoice : std::uint8_t {
  /// P1, AND THE DEFAULT: p with probability Z F2p(x, Q2) / (Z F2p + N F2n)
  /// on `PythiaBridgeOptions::f2_source`.  The inclusive rate off a nucleus
  /// IS Z F2p + N F2n, so this is the mixture the cross section describes.
  ByStructureFunctions,
  /// p with probability Z/A -- the x-independent limit of the above, and what
  /// the bridge did before 2026-08-29.  Wrong wherever F2n/F2p is: at x = 0.5
  /// it draws a 6Li proton half the time against a true 0.616.
  ByZN,
  Proton,
  Neutron
};

/// What the bridge does with a COHERENT event -- one carrying a
/// `Role::Pomeron` particle written by `Pipeline::make_coherent`.
enum class CoherentT2 : std::uint8_t {
  /// THE DEFAULT.  Build a third PYTHIA instance on a POMERON beam
  /// (`Beams:idA = 990`) and run the gamma*-Pomeron system through exactly
  /// the same surrogate + Lorentz map as the gamma*-nucleon one.
  /// docs/PYTHIA_BRIDGE.md sec. 12.
  Pomeron,
  /// Do not build the instance; `hadronize()` returns false on a coherent
  /// event and leaves the record at T0 (the `Role::HadronicX`
  /// pseudo-particle stays `Status::Final` and still carries the whole
  /// hadronic system, so the record still conserves).  Saves the third
  /// `init()` for a run that never touches the coherent channel.
  Off
};

struct PythiaBridgeOptions {
  /// PYTHIA's own `Random:seed` (used only for the fallback stream: the
  /// per-event randomness is driven by the `Rng` handed to `hadronize`).
  std::uint64_t seed = 19780503;

  /// Extra `pythia.readString` lines, applied after the bridge's own
  /// settings and before `init()`, so they win.
  std::vector<std::string> settings;

  /// 0 = `Print:quiet`, no banner, no listings; 1 = banner + init report;
  /// 2 = also `Next:numberShowEvent = 1` for the first event.
  int verbosity = 0;

  /// Re-sample the struck-quark flavour and call `pythia.next()` again this
  /// many times before giving up on an event.
  int max_retries = 4;

  /// The PYTHIA-side nucleon beam energy is `headroom` times the nominal
  /// per-nucleon energy.  The surrogate exists only while the physical
  /// W^2 fits inside the surrogate's s, and Fermi motion can push the
  /// per-nucleon momentum ~20 % above nominal, so the default leaves room.
  /// It is a pure frame choice: the surrogate's (W^2, Q^2, xi, x) do not
  /// depend on it.
  double headroom = 1.5;

  /// Flavours offered to the flavour sampler.
  bool include_strange = true;
  bool include_charm = true;
  bool include_bottom = false;

  /// PDFs are evaluated at max(Q^2, `q2_pdf_min`); PYTHIA's default
  /// NNPDF2.3 LO grid starts at 1 GeV^2 and the generator window reaches
  /// down to Q^2 = 0.7 GeV^2.
  double q2_pdf_min = 1.0;

  /// Build the second (neutron-beam) PYTHIA instance.  Turn it off to halve
  /// the initialisation cost when only protons are hadronized.
  bool with_neutron_instance = true;

  // --- the coherent (gamma*-Pomeron) tier -------------------------------

  /// Whether the third, POMERON-beam instance is built at all.
  CoherentT2 coherent_t2 = CoherentT2::Pomeron;

  /// PYTHIA's `PDF:PomSet`, the Pomeron parton densities the flavour sampler
  /// reads and PYTHIA's backward evolution starts from.  **6 is PYTHIA
  /// 8.317's own default** (H1 2006 Fit B, LO).  The list, read off
  /// `xmldoc/PDFSelection.xml` at PYTHIA 8.317 rather than remembered:
  ///
  ///     1        N x^a (1-x)^b toy, Q^2-INDEPENDENT (not a fit)
  ///     2        pi^0 distributions (not a Pomeron fit)
  ///     3,4,5    H1 2006 Fit A / Fit B / 2007 Jets, NLO
  ///     6        H1 2006 Fit B, LO             <- the default
  ///     7,8,9    ACTW B / D / SG, NLO, eps = 0.14
  ///     10       ACTW D, NLO, eps = 0.19
  ///     11       PomHISASD, the Angantyr rescaled proton -- REFUSED here
  ///     12,13    GKG18-DPDF Fit A / Fit B, LO
  ///     14,15    GKG18-DPDF Fit A / Fit B, NLO
  ///
  /// (6 is NOT "the only LO Q^2-dependent set", as this comment said until
  /// 2026-09-04: 12 and 13 are LO too.  It is the only LO *H1* set.)
  ///
  /// MEASURED SYSTEMATIC (2026-09-04, D4).  Scanned over all fifteen, 20 000
  /// coherent events each at 6Li config 1, seed 4242, `coherent_t2 =
  /// Pomeron`.  The T0 columns t, x_pom, q2, x and the event weight come out
  /// BIT-IDENTICAL across the FOURTEEN sets this constructor still admits --
  /// 1-10 and 12-15, one md5 ffd35a3a62b591c547e9ca2ac4301b5d over the lot,
  /// re-measured 2026-09-05 in 14.3 s.  Set 11 gave the same md5 in the
  /// original scan and is refused now (below), so fourteen is the count that
  /// REPRODUCES; this comment said "every set ... all fifteen" until
  /// 2026-09-05.  The Pomeron
  /// PDF enters only the flavour draw and PYTHIA's backward evolution, while
  /// |t|, x_P, M_X and the rate are fixed upstream by `CoherentSampler` /
  /// `CoherentXpomModel`.  So there is NO PomSet band on |t|, x_P, M_X or
  /// the cross section -- it is identically zero by construction, and
  /// quoting one would be meaningless.  The band is on the HADRONIC FINAL
  /// STATE, over the twelve genuine diffractive-PDF fits (3-10, 12-15),
  /// about the default set 6:
  ///
  ///     <n_charged>      3.964    -2.6 % / +9.9 %   (min set 9, max set 5)
  ///     <n_hadrons>      8.521    -2.6 % / +10.2 %  (min set 10, max set 5)
  ///     <p_T> [GeV]      0.342    -5.9 % / +10.9 %  (min set 5, max set 10)
  ///     kaon fraction    0.0646   -44 % / +55 %     -- a FACTOR 2.8
  ///
  /// Statistical error on <n_charged> is 0.017 at 20 000 events and the
  /// seed-to-seed scatter over 4242 / 777 / 31337 is <= 0.057, against a
  /// 0.50 spread across sets, so 20 000 events per set already resolves the
  /// band and there is no case for more.
  ///
  /// It is an ENVELOPE OVER RE-RUNS, one npz per set: the set changes the
  /// final state event by event and no per-event weight maps one set onto
  /// another.  Sets 1 and 2 are excluded from the band (a toy and a pi^0
  /// PDF, neither a Pomeron fit); set 11 is refused outright by the
  /// constructor.  THE EXCLUSION IS NOT FREE AND IS MEASURED (2026-09-05,
  /// same recipe): including 1 and 2 takes the <n_charged> low edge from
  /// -2.6 % (set 9, 3.8600) to -4.7 % (set 2, 3.7773) and the <n_hadrons>
  /// low edge from -2.6 % to -7.0 % (set 2, 7.9261), while the high edges and
  /// the <p_T> and kaon bands do not move.  So the band above is a
  /// TWELVE-FIT number: never attach it to "all 15 sets".  Set 6 is the only LO H1 set, so the band mixes LO and NLO
  /// DPDFs used in an LO Monte Carlo -- defensible for a systematic
  /// envelope, indefensible for a central value; say which is which.  Set 5
  /// is the only set that puts real weight on a charm initiator: 29.98 % of
  /// its events change when `include_charm` is turned off, against
  /// 10.4-11.1 % on GKG18 (12-15) and EXACTLY 0 on 3, 4, 6, 7, 8, 9, 10.
  int pom_set = 6;
  /// PYTHIA's `PDF:PomRescale`, the overall normalization of the H1/ACTW
  /// Pomeron sets (their momentum sum is arbitrary).  It cancels out of the
  /// bridge's flavour draw, which is normalized per event, and PYTHIA uses it
  /// only inside its own diffractive machinery -- so this is here for
  /// completeness and reproducibility, not because it changes anything the
  /// bridge produces.  1.0 is PYTHIA's default.
  double pom_rescale = 1.0;

  /// Default choice of the struck nucleon when the event carries no
  /// `Role::StruckNucleon`.
  ///
  /// A `Pipeline` INCLUSIVE event always names one (`InclusiveGenerator`
  /// writes `Role::StruckNucleon`, drawn from the same structure functions),
  /// so this fallback only fires for callers driving the bridge directly.
  NucleonChoice nucleon_choice = NucleonChoice::ByStructureFunctions;
  /// Unpolarized structure functions for `ByStructureFunctions`.  Null = the
  /// library's default `ToyF2`, which is also `InclusiveKernel`'s default;
  /// hand it the same backend the kernel was built with to keep the two
  /// draws consistent.
  std::shared_ptr<const UnpolSF> f2_source;
};

/// Per-run counters.
struct PythiaBridgeStats {
  std::uint64_t n_called = 0;      ///< hadronize() entries
  std::uint64_t n_ok = 0;
  std::uint64_t n_failed = 0;      ///< gave up after max_retries
  std::uint64_t n_retries = 0;     ///< extra pythia.next() calls
  std::uint64_t n_no_surrogate = 0;///< kinematics not representable
  /// P2: flavour offers refused because that flavour's own zeta_q had already
  /// run past 1.  Each one used to stay in the pool and burn a
  /// `pythia.next()` retry every time it was picked.
  std::uint64_t n_flavour_dropped = 0;
  std::uint64_t n_proton = 0, n_neutron = 0;
  /// Coherent events hadronized off the Pomeron beam (`Role::Pomeron`).
  std::uint64_t n_pomeron = 0;
  /// Coherent events whose flavour weights e_q^2 x f_q(beta, Q^2) all came
  /// out zero because the Pomeron grid has literally no quarks there (gluon
  /// fraction 1.000 at small beta until Q^2 ~ 1.5-1.75 GeV^2, which is ABOVE
  /// the default `q2_pdf_min` = 1.0 clamp), so the sampler fell back to the
  /// bare charge weights e_q^2 over the LIGHT flavours only.  The three H1
  /// 2006 sets carry no charm or bottom at ANY (beta, Q^2) -- Fit B LO
  /// (`pom_set` 6) and, equally, Fit A and Fit B NLO (3 and 4), because all
  /// three are PYTHIA's one `PomH1FitAB` class and its `xfUpdate` assigns
  /// `xc = xcbar = xb = xbbar = 0.` unconditionally
  /// (`PartonDistributions.cc:2630`) -- so on those the democratic limit
  /// keeps them at zero.  A DEFAULT coherent run sits on this branch for
  /// ~20 % of its events; non-zero is routine whenever the Q^2 window
  /// reaches below the grid's quark-support edge, not a misconfiguration.
  ///
  /// IT IS NOT A SYSTEMATIC.  MEASURED 2026-09-04 (D4.6), 4 000 coherent
  /// events per point at 6Li config 1, seed 4242, diffing the whole final
  /// state (`pid` and `p4` arrays) wholesale: raising `q2_pdf_min` from 1.0
  /// to 1.75 drives the fallback from 20.70 % to 0.00 % (set 6), 31.25 % to
  /// 0.00 % (set 3), 35.52 % to 29.32 % (set 4 -- THE EXCEPTION, and it
  /// travels with the sentence: on Fit B NLO the share stays high), 19.45 %
  /// to 2.65 % (set 12), 20.10 % to 3.23 % (13), 21.65 % to 5.12 % (15) --
  /// and produces a BIT-IDENTICAL final state in every one of those cases,
  /// set 4 included, which is why its residual 29 % costs nothing either.  The reason is structural, not luck: every
  /// Pomeron DPDF PYTHIA ships carries a single light-quark singlet (H1 sets
  /// it explicitly; the GKG18 LHAGrid1 files have columns -3..3 equal row by
  /// row), so e_q^2 x f_q is proportional to e_q^2 exactly over the light
  /// flavours and the normalised "fallback" IS the true draw.  The
  /// light-only restriction has therefore never discarded anything: it fires
  /// only where every weight vanishes, which is below the charm threshold,
  /// where charm is zero anyway.
  ///
  /// The one thing that is NOT free is the CLAMP, and it is a charm effect
  /// rather than a fallback effect: at `q2_pdf_min` = 3.0 the GKG18 sets
  /// (12, 13, 15) stop being bit-identical, because the clamp lifts charm
  /// above its threshold -- with `include_charm = false` the same comparison
  /// is bit-identical again.  So keep the light-only restriction (it is the
  /// correct guard against the `2e404b8` bug) and do not quote the fallback
  /// share as a modelling uncertainty; it is a bookkeeping artefact.
  /// Recorded per run since 2026-09-04 as `meta["n_pom_flavour_fallback"]`
  /// and `meta["pom_flavour_fallback_frac"]`, and printed at the run banner.
  std::uint64_t n_pom_flavour_fallback = 0;
  /// Events that took the DEPRECATED `Role::StruckCluster` branch -- i.e.
  /// arrived with no `Role::StruckNucleon`.  Non-zero on a `Pipeline` run
  /// means the run is at `Tier::T0`, and the whole record does not conserve
  /// (`PythiaBridge::NucleonInCluster`).
  std::uint64_t n_cluster_fallback = 0;
  double max_rescale_dev = 0.0;    ///< max |lambda - 1| of the mass repair
  double sum_w2 = 0.0;             ///< bookkeeping: mean W^2 of accepted events
};

// ---------------------------------------------------------------- bridge

/// One `PythiaBridge` owns one PYTHIA instance per nucleon type.  Not
/// thread-safe: give each thread its own bridge (the per-event randomness
/// comes from the caller's `Rng`, so events stay reproducible).
class PythiaBridge {
 public:
  /// DEPRECATED (2026-08-30, superseded by the T1 tier).  Hook: pick the
  /// nucleon struck inside a `Role::StruckCluster`.  Given the cluster
  /// four-vector and its (A, Z), it must return the nucleon four-vector and
  /// set `pdg_out` to 2212 or 2112.
  ///
  /// Nothing in the library reaches this any more: a `Pipeline` tagged event
  /// is resolved into a struck NUCLEON plus its partner spectators by
  /// `breakup.hpp` before the bridge sees it (`PipelineConfig::tier`), so
  /// `Role::StruckNucleon` -- the bridge's first-priority branch -- is
  /// always there.  The cluster branch survives only for a caller who builds
  /// a record by hand or deliberately runs at `Tier::T0`, and it does NOT
  /// conserve the whole record: it feeds PYTHIA p_cluster/A_c with nothing
  /// carrying the rest (measured 0.186 relative and one charge unit wrong on
  /// half the events, docs/T2_CHAIN.md).  Taking it logs one warning per
  /// bridge and counts every event in
  /// `PythiaBridgeStats::n_cluster_fallback`.
  using NucleonInCluster = std::function<Vec4(const Vec4& p_cluster, int a_c,
                                              int z_c, Rng& rng,
                                              int* pdg_out)>;

  /// Hook: pick the struck nucleon species (2212 / 2112) for an event with
  /// no explicit target.  The default follows `options().nucleon_choice`.
  using NucleonChooser = std::function<int(const Event& ev, Rng& rng)>;

  PythiaBridge(const BeamConfig& beams, PythiaBridgeOptions opt = {});
  ~PythiaBridge();

  PythiaBridge(const PythiaBridge&) = delete;
  PythiaBridge& operator=(const PythiaBridge&) = delete;

  /// Shower and hadronize the event in place.
  ///
  /// Reads: `Role::ScatteredElectron` (mandatory), `Role::BeamElectron`
  /// (mandatory), and a target -- `Role::Pomeron` (the coherent channel),
  /// else `Role::StruckNucleon`, else `Role::StruckCluster`, else the
  /// inclusive fallback P_ion/A.
  /// `Role::Spectator` / `Role::PartnerSpectator` particles are never
  /// touched and never recoil-corrected (the BeAGLE light-nucleus rule).
  ///
  /// Writes: `Role::Hadron`, `Status::Final` particles appended to
  /// `ev.particles`; any `Role::HadronicX` pseudo-particle is demoted to
  /// `Status::Intermediate`; a `Role::StruckNucleon` is appended with
  /// `Status::Intermediate` when the target was implicit.
  ///
  /// Returns false (and counts it) when PYTHIA vetoes or the kinematics are
  /// not representable; the event is then left unmodified.
  bool hadronize(Event& ev, Rng& rng);

  /// The hadrons produced by the last successful `hadronize` call.
  const std::vector<Particle>& last_hadrons() const;

  /// The target four-vector and PDG id used by the last call: the struck
  /// nucleon (2212/2112), or the Pomeron (990) on a coherent event.
  const Vec4& last_struck_nucleon() const;
  int last_struck_nucleon_pdg() const;
  /// The light-cone fraction handed to PYTHIA and the physical xi of the
  /// last call (they differ by the target's off-shellness; see the docs).
  double last_xi_pythia() const;
  double last_xi_physical() const;
  int last_quark_id() const;
  /// The common momentum rescale lambda applied by the mass repair of the
  /// last call, and the hadronic invariant masses it reconciled.
  double last_rescale() const;
  double last_w_pythia() const;
  double last_w_physical() const;

  const PythiaBridgeStats& stats() const;
  const PythiaBridgeOptions& options() const;

  void set_nucleon_in_cluster(NucleonInCluster hook);
  void set_nucleon_chooser(NucleonChooser hook);

  /// PYTHIA's own generated cross section for the surrogate stream, in mb.
  /// Meaningless as a physics number here (the hard process is ours, and
  /// LHAup strategy 3 hands PYTHIA a unit cross section); exposed only so
  /// that the example can print the bookkeeping.
  double pythia_sigma_gen_mb() const;

  /// The exact PYTHIA settings applied, in order -- what the docs quote.
  const std::vector<std::string>& applied_settings() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

}  // namespace lipolgen

#endif  // LIPOLGEN_PYTHIA_BRIDGE_HPP
