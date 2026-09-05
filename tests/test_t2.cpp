// T2 chain: Pipeline -> PythiaBridge -> (HepMC3) for the three physics
// channels the T2 tier is exercised on: inclusive, tagged 6Li-alpha, and
// coherent 6Li.  See docs/T2_CHAIN.md for the record layout and the two
// interface findings this file exists to pin down:
//
//  * TAGGED -- CLOSED by the T1 tier (`breakup.hpp`, 2026-08-30).  A tagged
//    event now names its struck NUCLEON (`Role::StruckNucleon`, the
//    off-shell P_X minus its on-shell partner spectators), so the bridge's
//    first-priority branch takes it verbatim and its `StruckCluster`
//    fallback never fires.  The whole record conserves to the numerical
//    floor with NO caller-side hook of any kind, and the "no surrogate"
//    tail collapses because the surrogate now has to reach one nucleon's
//    worth of momentum instead of a whole cluster's.
//
//    For the record, what this replaced: the v0 default `NucleonInCluster`
//    (p_cluster/A_c on shell, flavour drawn Z_c:N_c) broke the whole-record
//    balance by *0.186* relative and got the charge wrong on 157/300
//    events, and the `set_nucleon_in_cluster` workaround that fixed the
//    balance cost a ~5-15 % no-surrogate tail and always assigned the
//    cluster's own net charge.  Both are gone; the test below runs the
//    plain `bridge.hadronize` binding.
//
//  * COHERENT -- CLOSED by the gamma*-Pomeron tier (2026-08-30).  A
//    coherent event now carries `Role::Pomeron` (P_IP = P_ion - P_recoil,
//    written by `Pipeline::make_coherent`), which the bridge hadronizes on
//    a THIRD PYTHIA instance whose beam A is PYTHIA's Pomeron
//    (`Beams:idA = 990`) through the same surrogate + frame map as every
//    other channel (docs/PYTHIA_BRIDGE.md sec. 12).  The coherent tests
//    below therefore ASSERT at the same tight tolerances as the DIS tiers
//    -- whole-record 4-momentum 1e-9 relative, charge exact, M_had = M_X
//    to 1e-6 -- plus the tier's own guarantees: the intact recoil is never
//    touched, veto rate 0 at M_X >= 1.4 GeV, the light-only e_q^2 flavour
//    fallback below the LO Pomeron grid's quark-support edge, determinism,
//    and CoherentT2::Off degrading to a still-conserving T0 record.
//
//    For the record, what this replaced: a coherent event carried neither
//    `Role::StruckNucleon` nor `Role::StruckCluster`, so `hadronize()`
//    silently fell back to "a nucleon at rest in the ion frame" -- a
//    target unrelated to the true P_ion - P_recoil -- and the whole-record
//    balance was badly broken (worst relative 4-momentum residual 0.164,
//    charge wrong on 144/300 events) with no public hook to supply the
//    coherent target.  The v0 test could only report that residual under a
//    loose regression guard; that guard is gone.
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
#include <stdexcept>

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
  // Coherent bookkeeping: M_X of every ATTEMPTED event, whether it survived,
  // and the charged multiplicity and reconstructed hadronic mass when it did.
  double worst_mhad_rel = 0.0;         // |M_had - M_X| / M_X
  std::vector<double> mx_att;          // M_X [GeV], one per attempt
  std::vector<char> mx_ok;             // 1 = hadronized, 0 = vetoed
  std::vector<int> mx_nch;             // charged multiplicity, -1 when vetoed
};

// Runs `kN` events of `cfg`/`plan` through `bridge` and fills `st`.  `other`
// names the role that must (a) survive hadronization bit-identically and
// (b) enter the "beams in = Status::Final out" balance alongside e' and the
// hadrons: `Role::Spectator` for tagged, `Role::IntactRecoil` for coherent,
// `Role::Other` (i.e. none -- the per-nucleon balance already excludes it)
// for inclusive.
void run_channel(PipelineConfig cfg, const RunPlan& plan, PythiaBridge& bridge,
                 Role other, ChainStats* st,
                 Role target = Role::StruckNucleon, bool check_truth = true) {
  Particle other_before;
  bool have_other_before = false;

  cfg.hadronizer = [&](Event& ev, Rng& rng) {
    ++st->n_attempted;
    // COHERENT bookkeeping: M_X is known BEFORE the bridge runs, so a veto is
    // attributable to the M_X it was vetoed at.
    const double mx_this =
        (ev.kin.m_x2 > 0.0) ? std::sqrt(ev.kin.m_x2) : 0.0;
    st->mx_att.push_back(mx_this);
    st->mx_ok.push_back(0);
    st->mx_nch.push_back(-1);
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
    st->mx_ok.back() = 1;

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
    const Particle* pn_used = ev.find(target);
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
    st->mx_nch.back() = h.n_charged;
    const double exact = hfs_sigma_empz_exact(beam_e, pn_used->p, escat->p);
    st->worst_empz_exact = std::max(st->worst_empz_exact, std::fabs(h.sigma_empz - exact));
    if (check_truth) {
      // The CLOSED-FORM identity assumes a collinear on-shell target and the
      // per-nucleon y; neither holds for the Pomeron, so the coherent case
      // passes check_truth = false and relies on the exact form above.
      const double truth = hfs_sigma_empz_truth(beam_e.e, ev.kin.y, pn_used->p);
      st->worst_empz_truth_rel =
          std::max(st->worst_empz_truth_rel,
                   std::fabs(h.sigma_empz - truth) / std::fabs(truth));
    }
    // M_had == M_X: the hadrons the bridge wrote must reproduce the
    // diffractive mass the core generator drew (coherent), or W (otherwise).
    {
      Vec4 had;
      for (const auto& pp : ev.particles) {
        if (pp.role == Role::Hadron && pp.status == Status::Final) had = had + pp.p;
      }
      const double m2 = had.m2();
      const double m_want = std::sqrt(std::fabs((ev.kin.m_x2 > 0.0)
                                                    ? ev.kin.m_x2
                                                    : (beam_e + pn_used->p - escat->p).m2()));
      if (m2 > 0.0 && m_want > 0.0) {
        st->worst_mhad_rel = std::max(st->worst_mhad_rel,
                                      std::fabs(std::sqrt(m2) - m_want) / m_want);
      }
    }
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

TEST_CASE("T2 chain: tagged at T1, whole-record balance with NO caller-side "
          "hook (6Li-alpha, 7Li-alpha, d control)") {
  struct Row { PipelineChannel ch; const char* name; std::uint64_t seed_off; };
  const Row rows[] = {
      {PipelineChannel::TaggedLi6Alpha, "6Li-alpha", 1},
      {PipelineChannel::TaggedLi7Alpha, "7Li-alpha", 11},
      {PipelineChannel::TaggedDeuteronP, "d control", 21},
  };
  for (const Row& row : rows) {
    PipelineConfig cfg;
    cfg.channel = row.ch;
    cfg.isotope = channel_isotope(cfg.channel);
    cfg.beam_config = 1;
    cfg.n_events = kN;
    cfg.seed = kSeed;
    // The tabulated lithium tagging point exists for 6Li and 7Li only.
    if (row.ch != PipelineChannel::TaggedDeuteronP)
      cfg.optics_choice = OpticsChoice::Tagging;
    REQUIRE(cfg.tier == Tier::T1);        // the default, and the point here

    const BeamConfig bc =
        default_configs(cfg.isotope)[static_cast<std::size_t>(cfg.beam_config)];
    PythiaBridge bridge(bc, bridge_opts(kSeed + row.seed_off));
    // NO set_nucleon_in_cluster: the T1 record names the struck nucleon, so
    // the bridge's cluster branch is never reached.

    const RunPlan plan = (row.ch == PipelineChannel::TaggedLi7Alpha)
                             ? helicity_flip_plan(1.5, 0.7, 0.7)
                             : tensor_thirds_plan(0.7, 0.6);
    ChainStats st;
    run_channel(cfg, plan, bridge, Role::Spectator, &st);

    MESSAGE("tagged " << std::string(row.name) << " (T1): " << st.n_ok << "/"
                      << st.n_attempted << " hadronized, worst 4p rel "
                      << st.worst_p_rel << ", worst charge " << st.worst_q
                      << ", worst Sum(E-pz) exact " << st.worst_empz_exact
                      << " GeV, worst Sum(E-pz) truth rel "
                      << st.worst_empz_truth_rel
                      << " (the struck nucleon carries the Fermi p_T, so the "
                         "collinear-target closed form is only approximate; "
                         "the EXACT form above is not), worst pT "
                      << st.worst_pt << " GeV, " << st.secs << " s -> "
                      << (st.secs > 0 ? st.n_ok / st.secs : 0.0) << " ev/s");
    CHECK(st.n_attempted == kN);
    // The "no surrogate" tail: a few % at most now that the surrogate has to
    // reach ONE NUCLEON's momentum instead of a whole cluster's.
    CHECK(st.n_ok >= kN - kN / 20);
    CHECK(st.n_untouched_fail == 0);
    CHECK(st.n_electron_fail == 0);
    CHECK(st.worst_p_rel < 1e-9);      // measured ~1e-13
    CHECK(st.worst_q == 0.0);
    CHECK(st.worst_empz_exact < 1e-6);
    CHECK(st.worst_empz_truth_rel < 0.2);
    CHECK(st.worst_pt < 1e-8);
  }
}

TEST_CASE("T2 chain: tagged 7Li-alpha at T1 with --triton-sf ciofi-simula -- "
          "the spec's 'conservation T1+T2 1e-9 both options'") {
  // The tagged TEST_CASE above already covers the Hulthen default; this is
  // the SAME chain with the Ciofi degli Atti-Simula spectral function
  // (triton_sf.hpp), whose third channel (struck n -> a (p n) continuum)
  // the sequential model does not have.  The whole-record balance must not
  // care which model resolved the triton.
  PipelineConfig cfg;
  cfg.channel = PipelineChannel::TaggedLi7Alpha;
  cfg.isotope = channel_isotope(cfg.channel);
  cfg.beam_config = 1;
  cfg.n_events = kN;
  cfg.seed = kSeed;
  cfg.optics_choice = OpticsChoice::Tagging;
  cfg.triton_sf = TritonSfChoice::CiofiSimula;   // the one line under test
  REQUIRE(cfg.tier == Tier::T1);

  const BeamConfig bc =
      default_configs(cfg.isotope)[static_cast<std::size_t>(cfg.beam_config)];
  PythiaBridge bridge(bc, bridge_opts(kSeed + 31));

  ChainStats st;
  run_channel(cfg, helicity_flip_plan(1.5, 0.7, 0.7), bridge, Role::Spectator,
              &st);

  MESSAGE("tagged 7Li-alpha (T1, ciofi-simula): " << st.n_ok << "/"
          << st.n_attempted << " hadronized, worst 4p rel " << st.worst_p_rel
          << ", worst charge " << st.worst_q << ", worst Sum(E-pz) exact "
          << st.worst_empz_exact << " GeV, worst Sum(E-pz) truth rel "
          << st.worst_empz_truth_rel << ", worst pT " << st.worst_pt
          << " GeV, " << st.secs << " s -> "
          << (st.secs > 0 ? st.n_ok / st.secs : 0.0) << " ev/s");
  CHECK(st.n_attempted == kN);
  CHECK(st.n_ok >= kN - kN / 20);
  CHECK(st.n_untouched_fail == 0);
  CHECK(st.n_electron_fail == 0);
  CHECK(st.worst_p_rel < 1e-9);
  CHECK(st.worst_q == 0.0);
  CHECK(st.worst_empz_exact < 1e-6);
  CHECK(st.worst_empz_truth_rel < 0.2);
  CHECK(st.worst_pt < 1e-8);
}

TEST_CASE("T2 chain: coherent 6Li -- the gamma*-Pomeron tier, 300 events, "
          "whole-record balance") {
  PipelineConfig cfg;
  cfg.channel = PipelineChannel::CoherentLi6;
  cfg.isotope = "6Li";
  cfg.beam_config = 1;
  cfg.n_events = kN;
  cfg.seed = kSeed;
  cfg.optics_choice = OpticsChoice::Tagging;

  const BeamConfig bc = default_configs(cfg.isotope)[static_cast<std::size_t>(cfg.beam_config)];
  PythiaBridgeOptions bo = bridge_opts(kSeed + 2);
  REQUIRE(bo.coherent_t2 == CoherentT2::Pomeron);   // the default, and the point
  REQUIRE(bo.pom_set == 6);                         // PYTHIA's own default
  PythiaBridge bridge(bc, bo);

  ChainStats st;
  run_channel(cfg, tensor_thirds_plan(0.7, 0.6), bridge, Role::IntactRecoil,
              &st, Role::Pomeron, /*check_truth=*/false);

  MESSAGE("coherent 6Li (gamma*-Pomeron): " << st.n_ok << "/" << st.n_attempted
          << " hadronized, worst 4p rel " << st.worst_p_rel
          << ", worst charge " << st.worst_q
          << ", worst |M_had - M_X|/M_X " << st.worst_mhad_rel
          << ", worst Sum(E-pz) exact " << st.worst_empz_exact
          << " GeV, worst pT " << st.worst_pt << " GeV, " << st.secs
          << " s -> " << (st.secs > 0 ? st.n_ok / st.secs : 0.0) << " ev/s");

  CHECK(st.n_attempted == kN);
  CHECK(st.n_untouched_fail == 0);          // the intact recoil is never touched
  CHECK(st.n_electron_fail == 0);
  // THE POINT OF THE WHOLE TIER: whole-nucleus four-momentum and EXACT charge.
  CHECK(st.worst_p_rel < 1e-9);             // measured ~1e-13
  CHECK(st.worst_q == 0.0);                 // the Pomeron remnant is an antiquark
  CHECK(st.worst_mhad_rel < 1e-6);          // the hadrons ARE the drawn M_X
  CHECK(st.worst_empz_exact < 1e-6);
  CHECK(st.worst_pt < 1e-8);
  CHECK(bridge.stats().n_pomeron == static_cast<std::uint64_t>(st.n_ok));
  CHECK(bridge.stats().n_proton == 0);
  CHECK(bridge.stats().n_neutron == 0);
  // max_rescale_dev is max|lambda - 1| of the bridge's mass-repair rescale.
  // Before the Pomeron-instance surrogate compensation (pythia_bridge.cpp,
  // the `w2_sur` block) this sat at 6.3e-6: PYTHIA closes the meson-like
  // Pomeron-beam record on the electron beam's MASSIVE light-cone minus, a
  // constant Delta(m^2) = -(m_e^2/(2 E_e)) (2 E_A - Q^2/(2 E_e)) =
  // -3.896e-6 GeV^2 deficit on every event, amplified near the 1.2 GeV M_X
  // floor where d(sum E)/d lambda = sum p_i^2/E_i is small (soft, <n_ch> ~ 2
  // rest-frame momenta).  With the compensation the branch sits at the same
  // numerical floor as the nucleon ones (measured ~1e-11), so the tolerance
  // is the SAME 1e-6 the DIS-W tiers get -- do not loosen it: a regression
  // here means PYTHIA moved the electron or the compensation broke.
  MESSAGE("coherent max_rescale_dev = " << bridge.stats().max_rescale_dev);
  CHECK(bridge.stats().max_rescale_dev < 1e-6);
  // The two Pomeron-PDF settings must actually reach the instance: the
  // POMERON instance is built last, so `applied_settings()` carries its list.
  {
    const auto& as = bridge.applied_settings();
    const bool has_pomset =
        std::find(as.begin(), as.end(), std::string("PDF:PomSet = 6")) != as.end();
    const bool has_pomrescale =
        std::find(as.begin(), as.end(), std::string("PDF:PomRescale = 1")) != as.end();
    CHECK(has_pomset);
    CHECK(has_pomrescale);
  }

  // -- <n_ch> and the veto rate against M_X -----------------------------
  // The M_X floor is 1.2 GeV (COHERENT_MX_MIN_DEFAULT); the prototype
  // measured a 6 % PYTHIA veto there, 0 from 1.4 GeV up.  Bin and print.
  const double edges[] = {1.2, 1.4, 1.8, 2.5, 4.0, 7.0, 1e9};
  const int nb = 6;
  int n_bin[nb] = {0}, ok_bin[nb] = {0}, nch_bin[nb] = {0};
  double mx_bin[nb] = {0.0};
  for (std::size_t i = 0; i < st.mx_att.size(); ++i) {
    int b = 0;
    while (b < nb - 1 && st.mx_att[i] >= edges[b + 1]) ++b;
    ++n_bin[b];
    mx_bin[b] += st.mx_att[i];
    if (st.mx_ok[i]) { ++ok_bin[b]; nch_bin[b] += st.mx_nch[i]; }
  }
  std::ostringstream tab;
  tab << "\n  coherent T2: <n_ch> and veto rate vs M_X (" << kN << " events)\n"
      << "    M_X range        N    <M_X>   <n_ch>   veto\n";
  for (int b = 0; b < nb; ++b) {
    if (n_bin[b] == 0) continue;
    char line[160];
    std::snprintf(line, sizeof(line),
                  "    %5.2f - %-6.5g %4d  %7.3f  %7.3f  %5.3f\n", edges[b],
                  (b == nb - 1) ? 99.0 : edges[b + 1], n_bin[b],
                  mx_bin[b] / n_bin[b],
                  ok_bin[b] ? double(nch_bin[b]) / ok_bin[b] : 0.0,
                  1.0 - double(ok_bin[b]) / n_bin[b]);
    tab << line;
  }
  MESSAGE(tab.str());

  // Veto rate is EXACTLY zero from M_X = 1.4 GeV up, and the whole sample is
  // above the 1.2 GeV floor by construction.
  int n_hi = 0, ok_hi = 0, n_below_floor = 0;
  for (std::size_t i = 0; i < st.mx_att.size(); ++i) {
    if (st.mx_att[i] < COHERENT_MX_MIN_DEFAULT * (1.0 - 1e-9)) ++n_below_floor;
    if (st.mx_att[i] >= 1.4) { ++n_hi; ok_hi += st.mx_ok[i]; }
  }
  CHECK(n_below_floor == 0);
  CHECK(n_hi > 0);
  CHECK(ok_hi == n_hi);                    // veto rate 0 at M_X >= 1.4 GeV
  CHECK(st.n_ok > kN - kN / 10);           // and a few % at most overall

  SUBCASE("determinism: the same seed gives the identical coherent record") {
    PythiaBridge b2(bc, bridge_opts(kSeed + 2));
    ChainStats s2;
    run_channel(cfg, tensor_thirds_plan(0.7, 0.6), b2, Role::IntactRecoil, &s2,
                Role::Pomeron, /*check_truth=*/false);
    CHECK(s2.n_ok == st.n_ok);
    CHECK(s2.mx_att.size() == st.mx_att.size());
    for (std::size_t i = 0; i < s2.mx_att.size(); ++i) {
      CHECK(s2.mx_att[i] == st.mx_att[i]);
      CHECK(s2.mx_ok[i] == st.mx_ok[i]);
      CHECK(s2.mx_nch[i] == st.mx_nch[i]);
    }
  }
}

TEST_CASE("T2 chain: coherent 6Li -- CoherentT2::Off leaves the record at T0") {
  PipelineConfig cfg;
  cfg.channel = PipelineChannel::CoherentLi6;
  cfg.isotope = "6Li";
  cfg.beam_config = 1;
  cfg.n_events = 40;
  cfg.seed = kSeed;
  cfg.optics_choice = OpticsChoice::Tagging;

  const BeamConfig bc = default_configs(cfg.isotope)[static_cast<std::size_t>(cfg.beam_config)];
  PythiaBridgeOptions bo = bridge_opts(kSeed + 7);
  bo.coherent_t2 = CoherentT2::Off;
  PythiaBridge bridge(bc, bo);

  std::uint64_t n = 0, n_true = 0;
  double worst_p_rel = 0.0, worst_q = 0.0;
  cfg.hadronizer = [&](Event& ev, Rng& rng) {
    ++n;
    n_true += bridge.hadronize(ev, rng) ? 1 : 0;
  };
  Pipeline p(cfg, tensor_thirds_plan(0.7, 0.6));
  p.for_each([&](const Event& ev) {
    // No hadrons, so the T0 pseudo-particle X is still Status::Final and the
    // WHOLE record still closes -- turning the tier off costs fidelity, never
    // conservation.
    for (const auto& q : ev.particles) CHECK(q.role != Role::Hadron);
    const Vec4 want = ev.particles[0].p + ev.particles[1].p;
    const Vec4 r = want - ev.total_final();
    worst_p_rel = std::max(worst_p_rel,
                           std::max({std::fabs(r.e), std::fabs(r.px),
                                     std::fabs(r.py), std::fabs(r.pz)}) / want.e);
    worst_q = std::max(worst_q,
                       std::fabs(ev.particles[0].charge + ev.particles[1].charge -
                                 ev.total_charge_final()));
  });
  CHECK(n == 40);
  CHECK(n_true == 0);                       // every call refused, none crashed
  CHECK(bridge.stats().n_failed == n);
  CHECK(bridge.stats().n_pomeron == 0);
  CHECK(worst_p_rel < 1e-12);
  CHECK(worst_q == 0.0);
}

TEST_CASE("T2 chain: coherent e_q^2 flavour fallback -- small beta below the "
          "LO Pomeron grid's quark-support edge (Q^2 <~ 1.5), light-only") {
  // The second trap of the coherent prototype: at small beta the H1 2006
  // Fit B LO grid (PDF:PomSet = 6) is pure gluon until Q^2 ~ 1.5-1.75 GeV^2
  // -- every e_q^2 x f_q(beta, Q^2) weight is zero, and the default
  // q2_pdf_min = 1.0 clamp sits BELOW that edge, so a default coherent run
  // takes the e_q^2 fallback on ~20 % of its events (not just at Q^2 = 1
  // exactly).  Build the coherent record by hand in that region and check
  // the bridge (a) hadronizes it, (b) counts the fallback, (c) still
  // reproduces M_X, and (d) NEVER initiates a heavy flavour there: all three
  // H1 2006 grids -- Fit A NLO (3) and Fit B NLO (4) as well as Fit B LO (6),
  // one `PomH1FitAB` class between them -- carry no charm at ANY (beta, Q^2),
  // so e_q^2 x f_q gives charm
  // zero weight everywhere and the democratic fallback must not resurrect
  // it.  (d) is the regression guard for the bug where the fallback reused
  // the DIS offer list -- include_charm defaults to true -- and ~7 % of a
  // default coherent sample came out charm-initiated with zero PDF support.
  const BeamConfig bc = default_configs("6Li")[1];
  PythiaBridge bridge(bc, bridge_opts(kSeed + 9));

  const double e_e = bc.electron_energy;
  const double m_ion = bc.ion.mass();
  const double pz_ion = bc.ion_momentum_per_nucleon * bc.ion.A;
  const Vec4 p_ion{std::sqrt(pz_ion * pz_ion + m_ion * m_ion), 0.0, 0.0, pz_ion};

  // A coherent record at (Q^2, M_X), built exactly the way
  // `Pipeline::make_coherent` guarantees: (q + P_IP)^2 = M_X^2, P_IP
  // spacelike.
  auto make_corner_event = [&](double q2, double mx) {
    // Scattered electron, massless: k.k' = Q^2/2  ->  E' + k'_z = Q^2/(2 E_e).
    const double ep = 8.0;
    const double kz = q2 / (2.0 * e_e) - ep;
    const double kt = std::sqrt(ep * ep - kz * kz);
    const Vec4 k{e_e, 0.0, 0.0, -e_e};
    const Vec4 kp{ep, kt, 0.0, kz};
    const Vec4 q = k - kp;
    REQUIRE(std::fabs(-q.m2() - q2) < 1e-9);

    // The Pomeron: P_IP = eta * P_ion + a transverse kick, with eta the exact
    // positive root of (q + P_IP)^2 = M_X^2.  pt = 0.2 makes P_IP spacelike
    // (P_IP^2 = eta^2 M_A^2 - pt^2 < 0), like a real one.
    const double pt_ip = 0.2;
    // eta^2 M_A^2 + 2 eta (q.e E_A - q.pz p_A) + (q^2 - pt^2 - 2 q.px pt - M_X^2) = 0
    const double qa = m_ion * m_ion;
    const double qb = 2.0 * (q.e * p_ion.e - q.pz * p_ion.pz);
    const double qc = q.m2() - pt_ip * pt_ip - 2.0 * q.px * pt_ip - mx * mx;
    const double eta = (-qb + std::sqrt(qb * qb - 4.0 * qa * qc)) / (2.0 * qa);
    REQUIRE(eta > 0.0);
    REQUIRE(eta < 0.2);
    const Vec4 p_ip{eta * p_ion.e, pt_ip, 0.0, eta * p_ion.pz};
    REQUIRE(std::fabs((q + p_ip).m2() - mx * mx) < 1e-6);
    REQUIRE(p_ip.m2() < 0.0);               // spacelike, like the real thing

    Event ev;
    ev.channel = Channel::CoherentLi6;
    ev.kin.q2 = q2;
    ev.kin.m_x2 = mx * mx;
    Particle pb;
    pb.pdg = 11; pb.status = Status::Beam; pb.role = Role::BeamElectron;
    pb.p = k; pb.charge = -1.0;
    ev.particles.push_back(pb);
    Particle pi;
    pi.pdg = 1000030060; pi.status = Status::Beam; pi.role = Role::BeamIon;
    pi.p = p_ion; pi.mass = m_ion; pi.charge = 3.0;
    ev.particles.push_back(pi);
    Particle pe;
    pe.pdg = 11; pe.status = Status::Final; pe.role = Role::ScatteredElectron;
    pe.p = kp; pe.charge = -1.0;
    ev.particles.push_back(pe);
    Particle pomeron;
    pomeron.pdg = 990; pomeron.status = Status::Intermediate;
    pomeron.role = Role::Pomeron; pomeron.p = p_ip;
    pomeron.mass = -std::sqrt(-p_ip.m2()); pomeron.charge = 0.0;
    ev.particles.push_back(pomeron);
    return ev;
  };

  {
    const double mx = 10.0;
    Event ev = make_corner_event(1.0, mx);
    Rng rng(kSeed + 9, 1, 0, 0);
    const bool ok = bridge.hadronize(ev, rng);
    CHECK(ok);
    CHECK(bridge.stats().n_pom_flavour_fallback > 0);   // (b): the corner fired
    CHECK(bridge.stats().n_pomeron == 1);
    // (c): the hadrons ARE the drawn M_X, fallback or not.
    Vec4 had;
    double q_had = 0.0;
    for (const auto& pp : ev.particles) {
      if (pp.role == Role::Hadron && pp.status == Status::Final) {
        had = had + pp.p;
        q_had += pp.charge;
      }
    }
    CHECK(std::fabs(std::sqrt(had.m2()) - mx) / mx < 1e-6);
    CHECK(q_had == 0.0);                    // gamma* + Pomeron is neutral
  }

  // (d): the fallback pool is LIGHT-ONLY, across the whole quarkless region,
  // not just at the q2_pdf_min floor.  Before the light-only restriction the
  // charm share of one fallback draw was 8/20 (e_q^2 weights with c + cbar
  // in the pool), so 60 events miss it with probability 0.6^60 ~ 5e-14 even
  // if the seed moves.
  const std::uint64_t fb_before = bridge.stats().n_pom_flavour_fallback;
  int n_events = 0;
  const struct { double q2, mx; } corners[] = {{1.0, 10.0}, {1.4, 12.0}};
  for (const auto& c : corners) {
    for (int i = 0; i < 30; ++i) {
      Event ev = make_corner_event(c.q2, c.mx);
      Rng rng(kSeed + 9, 2, static_cast<std::uint64_t>(n_events), 0);
      REQUIRE(bridge.hadronize(ev, rng));   // M_X >> 1.4: veto rate is 0
      ++n_events;
      const int idq = bridge.last_quark_id();
      CHECK(std::abs(idq) <= 3);            // no charm/bottom off a pure-gluon grid
      CHECK(idq != 0);
    }
  }
  // Every one of those events sat below the grid's quark-support edge, so
  // every one took the fallback -- ~20 % of a DEFAULT coherent run does.
  CHECK(bridge.stats().n_pom_flavour_fallback ==
        fb_before + static_cast<std::uint64_t>(n_events));
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

TEST_CASE("T2 chain: HepMC3 round trip through the full chain (coherent "
          "6Li) -- the Role::Pomeron line, counts, attributes, byte "
          "determinism") {
  // The coherent twin of the tagged round trip above: the pdg-990 Pomeron
  // line (Status::Intermediate -> HepMC status 3, mass -sqrt|t|) must
  // survive the writer/reader, and two identical-seed runs must produce
  // byte-identical files -- which extends the full-record determinism
  // guarantee (only proven for the inclusive channel by the byte-compare
  // test above) to the coherent record.
  PipelineConfig cfg;
  cfg.channel = PipelineChannel::CoherentLi6;
  cfg.isotope = "6Li";
  cfg.beam_config = 1;
  cfg.n_events = 40;
  cfg.seed = kSeed;
  cfg.optics_choice = OpticsChoice::Tagging;

  const BeamConfig bc = default_configs(cfg.isotope)[static_cast<std::size_t>(cfg.beam_config)];

  auto run = [&](const std::string& path,
                 std::vector<std::size_t>* counts) {
    PythiaBridge bridge(bc, bridge_opts(kSeed + 8));
    PipelineConfig c = cfg;
    c.hadronizer = [&](Event& ev, Rng& rng) { bridge.hadronize(ev, rng); };
    Pipeline p(c, tensor_thirds_plan(0.7, 0.6));
    HepMC3Writer w(path);
    p.for_each([&](const Event& ev) {
      w.write(ev);
      if (counts) counts->push_back(ev.particles.size());
    });
    w.close();
  };

  const auto a = std::filesystem::temp_directory_path() / "lipolgen_t2_coh_rt_a.hepmc";
  const auto b = std::filesystem::temp_directory_path() / "lipolgen_t2_coh_rt_b.hepmc";
  std::vector<std::size_t> written_particle_counts;
  run(a.string(), &written_particle_counts);
  run(b.string(), nullptr);
  REQUIRE(written_particle_counts.size() == static_cast<std::size_t>(cfg.n_events));

  // -- byte determinism ---------------------------------------------------
  {
    std::ifstream fa(a, std::ios::binary), fb(b, std::ios::binary);
    REQUIRE(fa.good());
    REQUIRE(fb.good());
    std::ostringstream sa, sb;
    sa << fa.rdbuf();
    sb << fb.rdbuf();
    CHECK(sa.str() == sb.str());
    CHECK(sa.str().size() > 0);
  }
  std::filesystem::remove(b);

  // -- read back ----------------------------------------------------------
  HepMC3::ReaderAscii reader(a.string());
  std::size_t n_events_read = 0;
  std::size_t n_with_hadrons = 0;
  while (!reader.failed()) {
    HepMC3::GenEvent gev(HepMC3::Units::GEV, HepMC3::Units::MM);
    if (!reader.read_event(gev)) break;
    if (reader.failed()) break;
    REQUIRE(gev.particles().size() == written_particle_counts[n_events_read]);
    // Exactly one Pomeron line, at HepMC status 3 (documentation, never in a
    // status-1 sum).
    int n_pom = 0, n_status1 = 0;
    for (const auto& gp : gev.particles()) {
      if (gp->pid() == 990) {
        ++n_pom;
        CHECK(gp->status() == 3);
      }
      if (gp->status() == 1) ++n_status1;
    }
    CHECK(n_pom == 1);
    // status 1 beyond e' and the intact recoil means real T2 hadrons.
    if (n_status1 > 2) ++n_with_hadrons;
    REQUIRE(gev.attribute<HepMC3::StringAttribute>("channel") != nullptr);
    CHECK(gev.attribute<HepMC3::StringAttribute>("channel")->value() ==
          std::string("CoherentLi6"));
    REQUIRE(gev.run_info() != nullptr);
    ++n_events_read;
  }
  reader.close();
  std::filesystem::remove(a);

  CHECK(n_events_read == static_cast<std::size_t>(cfg.n_events));
  MESSAGE("coherent HepMC3 round trip: " << n_events_read
          << " events read back, " << n_with_hadrons << " carry T2 hadrons");
  CHECK(n_with_hadrons > 0);
}
#endif  // LIPOLGEN_HAVE_HEPMC3

#endif  // LIPOLGEN_HAVE_PYTHIA8
