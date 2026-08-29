// T2 chain: Pipeline -> PythiaBridge -> (HepMC3) for the three physics
// channels the T2 tier is exercised on: inclusive, tagged 6Li-alpha, and
// coherent 6Li.  See docs/T2_CHAIN.md for the record layout and the two
// interface findings this file exists to pin down:
//
//  * TAGGED.  The v0 default `NucleonInCluster` hook (p_cluster/A_c on
//    shell, flavour drawn Z_c:N_c -- docs/PYTHIA_BRIDGE.md section 6) does
//    not conserve the WHOLE-RECORD balance: it hands PYTHIA only 1/A_c of
//    the struck cluster's momentum, with no compensating particle for the
//    rest, and its stochastic proton/neutron draw does not track the
//    cluster's own integer charge.  Measured on 300 6Li-alpha events with
//    the plain `bridge.hadronize` binding: worst relative 4-momentum
//    residual on (spectator + hadrons + e') vs (beam e + beam ion) is
//    *0.186*, and charge is wrong on 157/300 events.  `set_nucleon_in_cluster`
//    -- a public, documented extension point, not a source patch -- fixes
//    both: hand PYTHIA the cluster's own off-shell four-vector WHOLE
//    (still "no Fermi smearing", the v0 promise) and its exact integer
//    charge.  That is `whole_cluster_hook` below, and with it the tagged
//    checks pass at the same tolerance as inclusive.  A complete fix would
//    also emit the cluster's non-struck nucleon as `Role::PartnerSpectator`
//    (the role exists in event.hpp; nothing currently fills it), which is
//    out of this file's scope (src/ is owned elsewhere).
//
//  * COHERENT.  A coherent event carries neither `Role::StruckNucleon` nor
//    `Role::StruckCluster` (the diffractive system X is not a struck
//    nucleon), so `hadronize()` silently falls back to "a nucleon at rest
//    in the ion frame" -- the ordinary INCLUSIVE fallback.  That fallback
//    nucleon has nothing to do with the true momentum transfer to the
//    recoil (P_ion - P_recoil is small and largely transverse; the
//    fallback carries the full per-nucleon longitudinal momentum p_u), so
//    the call succeeds mechanically (no crash, no excess vetoes) but the
//    whole-record balance is badly broken: measured worst relative
//    4-momentum residual 0.164, charge wrong on 144/300 events.  There is
//    no public hook that lets a caller supply the coherent target's
//    four-vector (unlike the tagged case), so this is reported here rather
//    than patched around: PythiaBridge v0 has no representation of a
//    coherent-diffractive final state.  The coherent test below therefore
//    verifies everything that DOES hold (the call succeeds, the
//    DIS-surrogate HFS identities against the fallback target it actually
//    used, the recoil is untouched, exactly one final electron,
//    determinism) and reports -- without asserting at the tight tolerance
//    -- the whole-record residual, with a loose regression guard.
#ifdef LIPOLGEN_HAVE_PYTHIA8

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "check_close.hpp"
#include "doctest.h"

#include "lipolgen/pipeline.hpp"
#include "lipolgen/pythia_bridge.hpp"

#ifdef LIPOLGEN_HAVE_HEPMC3
#include "lipolgen/hepmc_writer.hpp"
#include <HepMC3/Attribute.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/ReaderAscii.h>
#endif

using namespace lipolgen;

namespace {

constexpr std::uint64_t kSeed = 20260829;
constexpr int kN = 300;

// See the file header: hands PYTHIA the struck cluster's own off-shell
// four-vector whole (not 1/A_c of it) and its exact integer charge, so the
// whole-record balance closes exactly.  Installed through
// `PythiaBridge::set_nucleon_in_cluster` -- the public hook
// docs/PYTHIA_BRIDGE.md section 6 names for exactly this purpose.
Vec4 whole_cluster_hook(const Vec4& p_cluster, int /*a_c*/, int z_c, Rng&,
                        int* pdg_out) {
  *pdg_out = (z_c >= 1) ? 2212 : 2112;
  return p_cluster;
}

// Per-event bookkeeping shared by all three channels.
struct ChainStats {
  int n_attempted = 0;
  int n_ok = 0;               // bridge.hadronize() returned true
  int n_qcheck = 0;           // events the charge check ran on
  double worst_p_rel = 0.0;   // |beams_in - Event::total_final()| / E_in
  double worst_q = 0.0;       // |charge_in - Event::total_charge_final()|
  double worst_empz_exact = 0.0;      // hfs_sigma_empz_exact, absolute [GeV]
  double worst_empz_truth_rel = 0.0;  // hfs_sigma_empz_truth, relative
  double worst_pt = 0.0;              // |Sum p_T,had - (p_N,T - p_T,e')| [GeV]
  int n_untouched_fail = 0;   // "other" role particle changed by hadronize()
  int n_electron_fail = 0;    // not exactly one final e', or a hadron-role e
  double secs = 0.0;
};

// Runs `kN` events of `cfg`/`plan` through `bridge` and fills `st`.  `other`
// names the role that must (a) survive hadronization bit-identically and
// (b) enter the "beams in = Status::Final out" balance alongside e' and the
// hadrons: `Role::Spectator` for tagged, `Role::IntactRecoil` for coherent,
// `Role::Other` (i.e. none -- the per-nucleon balance already excludes it)
// for inclusive.
void run_channel(PipelineConfig cfg, const RunPlan& plan, PythiaBridge& bridge,
                 Role other, ChainStats* st) {
  Particle other_before;
  bool have_other_before = false;

  cfg.hadronizer = [&](Event& ev, Rng& rng) {
    ++st->n_attempted;
    if (other != Role::Other) {
      const Particle* o = ev.find(other);
      REQUIRE(o != nullptr);
      other_before = *o;
      have_other_before = true;
    } else {
      have_other_before = false;
    }
    if (!bridge.hadronize(ev, rng)) return;
    ++st->n_ok;

    // -- untouched spectator / recoil -----------------------------------
    if (have_other_before) {
      const Particle* o = ev.find(other);
      REQUIRE(o != nullptr);
      const bool same =
          o->pdg == other_before.pdg && o->status == other_before.status &&
          o->charge == other_before.charge && o->mass == other_before.mass &&
          o->p.e == other_before.p.e && o->p.px == other_before.p.px &&
          o->p.py == other_before.p.py && o->p.pz == other_before.p.pz;
      if (!same) ++st->n_untouched_fail;
    }

    // -- exactly one Role::ScatteredElectron, and PYTHIA's own copy of it is
    // never carried over as a Role::Hadron duplicate (test_pythia.cpp's
    // check, repeated through the full chain).  A hadronic Dalitz pair or a
    // semileptonic decay electron is real physics and legitimately shares
    // pdg 11 with e' without being a duplicate of it, so the check is on
    // Role + momentum match, not on the pdg code alone.
    const Particle* escat = ev.find(Role::ScatteredElectron);
    REQUIRE(escat != nullptr);
    int n_scat = 0;
    for (const auto& p : ev.particles) if (p.role == Role::ScatteredElectron) ++n_scat;
    bool dup_hadron_e = false;
    for (const auto& p : ev.particles) {
      if (p.role == Role::Hadron && p.pdg == 11 &&
          std::fabs(p.p.e - escat->p.e) < 1e-9 &&
          std::fabs(p.p.pz - escat->p.pz) < 1e-9) {
        dup_hadron_e = true;
      }
    }
    if (n_scat != 1 || dup_hadron_e) ++st->n_electron_fail;

    // -- 4-momentum / charge: "beams in = Status::Final out" -------------
    // The struck nucleon the bridge actually fed to PYTHIA -- appended by
    // hadronize() itself (Role::StruckNucleon, Status::Intermediate) for
    // every implicit target, or already present and unchanged for the
    // explicit-target (inclusive) case.
    const Particle* pn_used = ev.find(Role::StruckNucleon);
    REQUIRE(pn_used != nullptr);
    const Vec4 beam_e = ev.particles[0].p;
    const Vec4 beam_ion = ev.particles[1].p;
    const double q_e = ev.particles[0].charge, q_ion = ev.particles[1].charge;

    // Inclusive's own balance is PER-NUCLEON (pipeline.hpp: the (A-1)
    // remnant is never written), so "beams in" there is e + the struck
    // nucleon, not e + the whole ion; tagged and coherent close the WHOLE
    // nucleus.
    const bool per_nucleon = (other == Role::Other);
    const Vec4 want_p = per_nucleon ? (beam_e + pn_used->p) : (beam_e + beam_ion);
    const double want_q = per_nucleon ? (q_e + pn_used->charge) : (q_e + q_ion);

    const Vec4 got_p = ev.total_final();
    const Vec4 r = want_p - got_p;
    st->worst_p_rel = std::max(
        st->worst_p_rel,
        std::max({std::fabs(r.e), std::fabs(r.px), std::fabs(r.py),
                  std::fabs(r.pz)}) /
            want_p.e);
    ++st->n_qcheck;
    st->worst_q = std::max(st->worst_q, std::fabs(want_q - ev.total_charge_final()));

    // -- HFS truth identities, against the target the bridge actually used
    const HfsSummary h = hfs_summary(ev);
    const double exact = hfs_sigma_empz_exact(beam_e, pn_used->p, escat->p);
    st->worst_empz_exact = std::max(st->worst_empz_exact, std::fabs(h.sigma_empz - exact));
    const double truth = hfs_sigma_empz_truth(beam_e.e, ev.kin.y, pn_used->p);
    st->worst_empz_truth_rel =
        std::max(st->worst_empz_truth_rel, std::fabs(h.sigma_empz - truth) / std::fabs(truth));
    const double want_px = pn_used->p.px - escat->p.px;
    const double want_py = pn_used->p.py - escat->p.py;
    st->worst_pt = std::max(st->worst_pt, std::max(std::fabs(h.px - want_px),
                                                    std::fabs(h.py - want_py)));
  };

  // The real work happens inside `cfg.hadronizer` above, called once per
  // event by `Pipeline::event()`; the sink here has nothing left to do.
  // `cfg.n_events` (set by every caller below) is what caps the run at kN.
  Pipeline p(cfg, plan);
  const auto t0 = std::chrono::steady_clock::now();
  p.for_each([](const Event&) {});
  st->secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

PythiaBridgeOptions bridge_opts(std::uint64_t seed) {
  PythiaBridgeOptions o;
  o.seed = seed;
  return o;
}

}  // namespace

// ---------------------------------------------------------------------------

TEST_CASE("T2 chain: inclusive, 300 events, per-nucleon balance") {
  PipelineConfig cfg;
  cfg.channel = PipelineChannel::Inclusive;
  cfg.isotope = "6Li";
  cfg.beam_config = 1;
  cfg.n_events = kN;
  cfg.seed = kSeed;

  const BeamConfig bc = default_configs(cfg.isotope)[static_cast<std::size_t>(cfg.beam_config)];
  PythiaBridge bridge(bc, bridge_opts(kSeed));

  ChainStats st;
  run_channel(cfg, tensor_thirds_plan(0.7, 0.6), bridge, Role::Other, &st);

  MESSAGE("inclusive: " << st.n_ok << "/" << st.n_attempted
                        << " hadronized, worst 4p rel " << st.worst_p_rel
                        << ", worst charge " << st.worst_q
                        << ", worst Sum(E-pz) exact " << st.worst_empz_exact
                        << " GeV, worst Sum(E-pz) truth rel "
                        << st.worst_empz_truth_rel << ", worst pT "
                        << st.worst_pt << " GeV, " << st.secs << " s -> "
                        << (st.secs > 0 ? st.n_ok / st.secs : 0.0) << " ev/s");
  CHECK(st.n_attempted == kN);
  CHECK(st.n_ok > kN - 30);          // PYTHIA vetoes/no-surrogate: a small tail
  CHECK(st.n_electron_fail == 0);
  CHECK(st.worst_p_rel < 1e-9);      // exact per-nucleon balance (measured ~4e-13)
  CHECK(st.worst_q == 0.0);
  CHECK(st.worst_empz_exact < 1e-6);           // exact identity, any target
  CHECK(st.worst_empz_truth_rel < 1e-3);       // docs/PYTHIA_BRIDGE.md sec.8; measured ~4e-13
  CHECK(st.worst_pt < 1e-8);
}

TEST_CASE("T2 chain: tagged 6Li-alpha, 300 events, whole-record balance "
          "(with the whole-cluster hook -- see file header)") {
  PipelineConfig cfg;
  cfg.channel = PipelineChannel::TaggedLi6Alpha;
  cfg.isotope = channel_isotope(cfg.channel);
  cfg.beam_config = 1;
  cfg.n_events = kN;
  cfg.seed = kSeed;
  cfg.optics_choice = OpticsChoice::Tagging;

  const BeamConfig bc = default_configs(cfg.isotope)[static_cast<std::size_t>(cfg.beam_config)];
  PythiaBridge bridge(bc, bridge_opts(kSeed + 1));
  bridge.set_nucleon_in_cluster(whole_cluster_hook);

  ChainStats st;
  run_channel(cfg, tensor_thirds_plan(0.7, 0.6), bridge, Role::Spectator, &st);

  MESSAGE("tagged 6Li-alpha: " << st.n_ok << "/" << st.n_attempted
                               << " hadronized, worst 4p rel " << st.worst_p_rel
                               << ", worst charge " << st.worst_q
                               << ", worst Sum(E-pz) exact " << st.worst_empz_exact
                               << " GeV, worst Sum(E-pz) truth rel "
                               << st.worst_empz_truth_rel
                               << " (looser than inclusive: P_X carries the "
                                  "spectator's recoil pT, so the collinear-"
                                  "target 'truth' formula is only "
                                  "approximate here), worst pT " << st.worst_pt
                               << " GeV, " << st.secs << " s -> "
                               << (st.secs > 0 ? st.n_ok / st.secs : 0.0)
                               << " ev/s");
  CHECK(st.n_attempted == kN);
  CHECK(st.n_ok > kN - 60);          // the whole-cluster hook pushes some
                                      // events off the surrogate's reach
                                      // (measured ~5%); see docs/T2_CHAIN.md
  CHECK(st.n_untouched_fail == 0);
  CHECK(st.n_electron_fail == 0);
  CHECK(st.worst_p_rel < 1e-9);      // measured ~1e-13 with the whole-cluster hook
  CHECK(st.worst_q == 0.0);
  CHECK(st.worst_empz_exact < 1e-6);
  // Measured ~3.1e-3: P_X is not collinear (the spectator carries recoil
  // p_T), unlike the docs' assumed-collinear target, so the closed-form
  // "truth" formula is only approximate here; the EXACT form above is
  // still exact to the numerical floor.
  CHECK(st.worst_empz_truth_rel < 1e-2);
  CHECK(st.worst_pt < 1e-8);
}

TEST_CASE("T2 chain: coherent 6Li, 300 events -- the call works, but "
          "whole-record conservation does NOT hold (see file header: "
          "PythiaBridge v0 has no coherent-diffractive target)") {
  PipelineConfig cfg;
  cfg.channel = PipelineChannel::CoherentLi6;
  cfg.isotope = "6Li";
  cfg.beam_config = 1;
  cfg.n_events = kN;
  cfg.seed = kSeed;
  cfg.optics_choice = OpticsChoice::Tagging;

  const BeamConfig bc = default_configs(cfg.isotope)[static_cast<std::size_t>(cfg.beam_config)];
  PythiaBridge bridge(bc, bridge_opts(kSeed + 2));
  // No hook applies here: hadronize() never looks for a coherent target
  // (there is no Role for one), so it always takes the plain inclusive
  // fallback (a nucleon at rest in the ion frame) regardless of what a
  // caller might install.

  ChainStats st;
  run_channel(cfg, tensor_thirds_plan(0.7, 0.6), bridge, Role::IntactRecoil, &st);

  MESSAGE("coherent 6Li: " << st.n_ok << "/" << st.n_attempted
                           << " hadronized (the call works mechanically), "
                              "worst 4p rel " << st.worst_p_rel
                           << " and worst charge " << st.worst_q
                           << " -- NOT conserved (measured ~0.16 / 1 with the "
                              "fallback target; see docs/T2_CHAIN.md), worst "
                              "Sum(E-pz) exact " << st.worst_empz_exact
                           << " GeV, worst Sum(E-pz) truth rel "
                           << st.worst_empz_truth_rel << ", worst pT "
                           << st.worst_pt << " GeV, " << st.secs << " s -> "
                           << (st.secs > 0 ? st.n_ok / st.secs : 0.0)
                           << " ev/s");
  CHECK(st.n_attempted == kN);
  CHECK(st.n_ok > kN - 30);
  CHECK(st.n_untouched_fail == 0);   // the recoil itself is still never touched
  CHECK(st.n_electron_fail == 0);
  // These two hold because they are the bridge's OWN internal identities,
  // exact relative to whatever target it used, independent of whether that
  // target has anything to do with the physical coherent process.
  CHECK(st.worst_empz_exact < 1e-6);
  CHECK(st.worst_empz_truth_rel < 1e-3);
  CHECK(st.worst_pt < 1e-8);
  // Whole-record momentum/charge conservation is DELIBERATELY NOT asserted
  // at the tight tolerance the other two channels meet: it does not hold.
  // The bound below is only a regression guard on the size of the known gap.
  CHECK(st.worst_p_rel < 0.30);
  CHECK(st.worst_q <= 1.0);
}

// ---------------------------------------------------------------------------

TEST_CASE("T2 chain: determinism -- same seed gives the identical HepMC3 "
          "content") {
#ifdef LIPOLGEN_HAVE_HEPMC3
  auto run = [&](const std::string& path) {
    PipelineConfig cfg;
    cfg.channel = PipelineChannel::Inclusive;
    cfg.isotope = "6Li";
    cfg.beam_config = 1;
    cfg.n_events = 25;
    cfg.seed = kSeed;
    const BeamConfig bc =
        default_configs(cfg.isotope)[static_cast<std::size_t>(cfg.beam_config)];
    PythiaBridge bridge(bc, bridge_opts(kSeed + 3));
    cfg.hadronizer = [&](Event& ev, Rng& rng) { bridge.hadronize(ev, rng); };
    Pipeline p(cfg, tensor_thirds_plan(0.7, 0.6));
    HepMC3Writer w(path);
    p.for_each([&](const Event& ev) { w.write(ev); });
    w.close();
  };

  const auto a = std::filesystem::temp_directory_path() / "lipolgen_t2_det_a.hepmc";
  const auto b = std::filesystem::temp_directory_path() / "lipolgen_t2_det_b.hepmc";
  run(a.string());
  run(b.string());

  std::ifstream fa(a, std::ios::binary), fb(b, std::ios::binary);
  REQUIRE(fa.good());
  REQUIRE(fb.good());
  std::ostringstream sa, sb;
  sa << fa.rdbuf();
  sb << fb.rdbuf();
  CHECK(sa.str() == sb.str());
  CHECK(sa.str().size() > 0);
  std::filesystem::remove(a);
  std::filesystem::remove(b);
#else
  MESSAGE("built without HepMC3: determinism checked on the Event record "
          "shape only, not HepMC3 content (see test_pythia.cpp)");
#endif
}

#ifdef LIPOLGEN_HAVE_HEPMC3
TEST_CASE("T2 chain: HepMC3 round trip through the full chain (tagged "
          "6Li-alpha) -- counts and attributes") {
  PipelineConfig cfg;
  cfg.channel = PipelineChannel::TaggedLi6Alpha;
  cfg.isotope = channel_isotope(cfg.channel);
  cfg.beam_config = 1;
  cfg.n_events = 40;
  cfg.seed = kSeed;
  cfg.optics_choice = OpticsChoice::Tagging;

  const BeamConfig bc = default_configs(cfg.isotope)[static_cast<std::size_t>(cfg.beam_config)];
  PythiaBridge bridge(bc, bridge_opts(kSeed + 4));
  bridge.set_nucleon_in_cluster(whole_cluster_hook);
  cfg.hadronizer = [&](Event& ev, Rng& rng) { bridge.hadronize(ev, rng); };
  Pipeline p(cfg, tensor_thirds_plan(0.7, 0.6));

  const auto path = std::filesystem::temp_directory_path() / "lipolgen_t2_roundtrip.hepmc";
  std::vector<std::size_t> written_particle_counts;
  {
    HepMC3Writer w(path.string());
    p.for_each([&](const Event& ev) {
      w.write(ev);
      written_particle_counts.push_back(ev.particles.size());
    });
    w.close();
  }
  REQUIRE(written_particle_counts.size() == static_cast<std::size_t>(cfg.n_events));

  HepMC3::ReaderAscii reader(path.string());
  std::size_t n_events_read = 0;
  std::size_t n_with_hadrons = 0;
  while (!reader.failed()) {
    HepMC3::GenEvent gev(HepMC3::Units::GEV, HepMC3::Units::MM);
    if (!reader.read_event(gev)) break;
    if (reader.failed()) break;
    REQUIRE(gev.particles().size() == written_particle_counts[n_events_read]);
    // status 1 (final) beyond e' and the spectator means real T2 hadrons.
    int n_status1 = 0;
    for (const auto& gp : gev.particles()) if (gp->status() == 1) ++n_status1;
    if (n_status1 > 2) ++n_with_hadrons;
    // spot-check one attribute of each kind the writer sets.
    REQUIRE(gev.attribute<HepMC3::StringAttribute>("channel") != nullptr);
    CHECK(gev.attribute<HepMC3::StringAttribute>("channel")->value() ==
          std::string("TaggedLi6Alpha"));
    REQUIRE(gev.attribute<HepMC3::DoubleAttribute>("spin_J") != nullptr);
    REQUIRE(gev.attribute<HepMC3::IntAttribute>("run") != nullptr);
    REQUIRE(gev.run_info() != nullptr);
    ++n_events_read;
  }
  reader.close();
  std::filesystem::remove(path);

  CHECK(n_events_read == static_cast<std::size_t>(cfg.n_events));
  MESSAGE("HepMC3 round trip: " << n_events_read << " events read back, "
                                << n_with_hadrons << " carry T2 hadrons");
  CHECK(n_with_hadrons > 0);
}
#endif  // LIPOLGEN_HAVE_HEPMC3

#endif  // LIPOLGEN_HAVE_PYTHIA8
