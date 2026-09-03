// Round-trip tests for lipolgen::HepMC3Writer against a synthetic event built
// by hand: e + 6Li beams (10 GeV x 99.5 GeV/u), scattered e', an alpha
// spectator, and a status-3 hadronic-system pseudo-particle X. See
// docs/HEPMC3_CONVENTION.md for the attribute/layout convention this checks.
#include "doctest.h"
#include "lipolgen/event.hpp"
#include "lipolgen/hepmc_writer.hpp"
#include "lipolgen/rc.hpp"

#include <HepMC3/Attribute.h>
#include <HepMC3/GenCrossSection.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/ReaderAscii.h>

#include <cmath>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace lipolgen;
using HepMC3::DoubleAttribute;
using HepMC3::IntAttribute;
using HepMC3::StringAttribute;

namespace {

// Ground-state nuclear masses [GeV], AME2020 atomic masses less electrons;
// matches PolarizedLithiumSim tools/fullsim/ion_gun_hepmc.py MASS table.
constexpr double kLi6Mass = 5.601518702;
constexpr double kAlphaMass = 3.727379407;
constexpr double kElectronMass = 0.000511;
constexpr double kTol = 1e-9;

Event make_synthetic_event() {
  Event ev;
  ev.number = 42;
  ev.channel = Channel::TaggedLi6Alpha;
  ev.weight = 1.0;
  ev.spin_weights = {1.05, 0.95, 1.02, 0.98};
  ev.xsec_pb = 123.456;
  ev.xsec_err_pb = 1.234;

  ev.spin.j = 1.0;
  ev.spin.m_ion = 1.0;
  ev.spin.m_struck = 0.5;
  ev.spin.lam_e = -1;
  ev.spin.pe = 0.8;
  ev.spin.theta_s = 0.3;
  ev.spin.phi_s = 1.2;
  ev.spin.pz = 8.0 / 13.0;
  ev.spin.pzz = 4.0 / 13.0;
  ev.spin.category = "pp";
  ev.spin.run = 7;
  ev.spin.bunch = 3;

  ev.kin.x = 0.25;
  ev.kin.q2 = 4.5;
  ev.kin.y = 0.6;
  ev.kin.phi = 1.1;
  ev.kin.k = 0.31;
  ev.kin.cos_theta_k = -0.2;
  ev.kin.phi_k = 2.0;
  ev.kin.alpha_s = 1.0;
  ev.kin.pt_s = 0.05;
  ev.kin.t = -0.02;
  ev.kin.x_pom = 0.0;

  const double p_e = 10.0;
  Particle beam_e;
  beam_e.pdg = 11;
  beam_e.status = Status::Beam;
  beam_e.role = Role::BeamElectron;
  beam_e.p = {std::sqrt(p_e * p_e + kElectronMass * kElectronMass), 0.0, 0.0, -p_e};
  beam_e.mass = kElectronMass;
  beam_e.charge = -1.0;

  const double p_ion = 6.0 * 99.5;  // 6Li, 99.5 GeV/u
  Particle beam_ion;
  beam_ion.pdg = 1000030060;
  beam_ion.status = Status::Beam;
  beam_ion.role = Role::BeamIon;
  beam_ion.p = {std::sqrt(p_ion * p_ion + kLi6Mass * kLi6Mass), 0.0, 0.0, p_ion};
  beam_ion.mass = kLi6Mass;
  beam_ion.charge = 3.0;

  Particle escat;
  escat.pdg = 11;
  escat.status = Status::Final;
  escat.role = Role::ScatteredElectron;
  {
    const double px = 1.5, py = 0.5, pz = -7.8;
    escat.p = {std::sqrt(px * px + py * py + pz * pz + kElectronMass * kElectronMass), px, py, pz};
  }
  escat.mass = kElectronMass;
  escat.charge = -1.0;
  escat.mother1 = 0;  // beam electron
  escat.pol = -1.0;   // helicity label, exercises the "pol" attribute

  Particle gamma;
  gamma.pdg = 22;
  gamma.status = Status::Intermediate;
  gamma.role = Role::VirtualPhoton;
  gamma.p = beam_e.p - escat.p;
  gamma.mass = std::sqrt(std::max(gamma.p.m2(), 0.0));
  gamma.mother1 = 0;  // beam electron

  Particle spectator;
  spectator.pdg = 1000020040;  // alpha
  spectator.status = Status::Final;
  spectator.role = Role::Spectator;
  {
    const double px = 0.3, py = -0.2, pz = 590.0;
    spectator.p = {std::sqrt(px * px + py * py + pz * pz + kAlphaMass * kAlphaMass), px, py, pz};
  }
  spectator.mass = kAlphaMass;
  spectator.charge = 2.0;
  spectator.mother1 = 1;  // beam ion

  Particle x;
  x.pdg = 92;  // generic hadronic-system pseudo-particle
  x.status = Status::Intermediate;
  x.role = Role::HadronicX;
  x.p = (beam_e.p + beam_ion.p) - escat.p - spectator.p;  // exact by construction
  x.mass = std::sqrt(std::max(x.p.m2(), 0.0));
  x.charge = beam_ion.charge - spectator.charge;
  x.mother1 = 3;  // virtual photon (index of `gamma` below)

  ev.particles = {beam_e, beam_ion, escat, gamma, spectator, x};
  return ev;
}

// Minimal event pinning the writer's electron generated-mass rule
// (hepmc_writer.cpp): a |pdg|==11 particle with Particle::mass == 0.0 is
// written with generated_mass = the PDG electron mass; every other particle
// (massive electrons included) keeps its own Particle::mass verbatim.
Event make_electron_mass_rule_event() {
  Event ev;
  ev.number = 7;
  ev.channel = Channel::Inclusive;
  ev.weight = 1.0;
  ev.xsec_pb = 1.0;
  ev.xsec_err_pb = 0.1;

  // Beam electron, built massless (the core's standard DIS convention, per
  // the comment in hepmc_writer.cpp) with E == |p| exactly on-shell.
  Particle beam_e;
  beam_e.pdg = 11;
  beam_e.status = Status::Beam;
  beam_e.role = Role::BeamElectron;
  {
    const double p_e = 10.0;
    beam_e.p = {p_e, 0.0, 0.0, -p_e};
  }
  beam_e.mass = 0.0;
  beam_e.charge = -1.0;

  Particle beam_ion;
  beam_ion.pdg = 1000030060;
  beam_ion.status = Status::Beam;
  beam_ion.role = Role::BeamIon;
  {
    const double p_ion = 6.0 * 99.5;  // 6Li, 99.5 GeV/u
    beam_ion.p = {std::sqrt(p_ion * p_ion + kLi6Mass * kLi6Mass), 0.0, 0.0, p_ion};
  }
  beam_ion.mass = kLi6Mass;
  beam_ion.charge = 3.0;

  // Scattered electron, also built massless with E == |p|: the writer must
  // rewrite generated_mass for this one too.
  Particle escat;
  escat.pdg = 11;
  escat.status = Status::Final;
  escat.role = Role::ScatteredElectron;
  {
    const double px = 1.5, py = 0.5, pz = -7.8;
    escat.p = {std::sqrt(px * px + py * py + pz * pz), px, py, pz};
  }
  escat.mass = 0.0;
  escat.charge = -1.0;

  // A deliberately massive |pdg|==11 particle: Particle::mass != 0.0, so the
  // rule must leave it alone and keep exactly this mass.
  Particle massive_electron;
  massive_electron.pdg = 11;
  massive_electron.status = Status::Final;
  massive_electron.role = Role::Other;
  {
    const double px = 0.2, py = 0.1, pz = 3.0;
    massive_electron.p = {std::sqrt(px * px + py * py + pz * pz + kElectronMass * kElectronMass),
                           px, py, pz};
  }
  massive_electron.mass = kElectronMass;
  massive_electron.charge = -1.0;

  // A non-electron built massless (photon): the rule is keyed on |pdg| == 11
  // so this must keep generated_mass == 0, not get the PDG-mass rewrite.
  Particle gamma;
  gamma.pdg = 22;
  gamma.status = Status::Intermediate;
  gamma.role = Role::VirtualPhoton;
  {
    const double px = 0.05, py = -0.05, pz = 1.0;
    gamma.p = {std::sqrt(px * px + py * py + pz * pz), px, py, pz};
  }
  gamma.mass = 0.0;

  // A massless POSITRON: the rule is keyed on |pdg| == 11, so pdg == -11
  // must get the PDG-mass rewrite exactly like the electrons.
  Particle positron;
  positron.pdg = -11;
  positron.status = Status::Final;
  positron.role = Role::Other;
  {
    const double px = -0.3, py = 0.4, pz = 2.0;
    positron.p = {std::sqrt(px * px + py * py + pz * pz), px, py, pz};
  }
  positron.mass = 0.0;
  positron.charge = 1.0;

  // mother1/mother2 left at -1: every non-beam particle hangs directly off
  // the primary vertex, which is all this test's mass-rule check needs.
  ev.particles = {beam_e, beam_ion, escat, massive_electron, gamma, positron};
  return ev;
}

std::filesystem::path temp_hepmc_path(const std::string& name) {
  return std::filesystem::temp_directory_path() / name;
}

}  // namespace

TEST_CASE("HepMC3Writer round-trips a synthetic tagged event") {
  const Event ev = make_synthetic_event();
  const auto path = temp_hepmc_path("lipolgen_test_hepmc_roundtrip.hepmc");

  {
    HepMC3Writer writer(path.string());
    writer.write(ev);
    writer.close();
  }

  HepMC3::ReaderAscii reader(path.string());
  HepMC3::GenEvent read_ev(HepMC3::Units::GEV, HepMC3::Units::MM);
  REQUIRE(reader.read_event(read_ev));
  REQUIRE_FALSE(reader.failed());
  reader.close();
  std::filesystem::remove(path);

  SUBCASE("particle count, vertex count") {
    CHECK(read_ev.particles().size() == ev.particles.size());
    CHECK(read_ev.vertices().size() == 2);  // primary vertex + photon-absorption vertex
  }

  SUBCASE("statuses, PDG ids, masses in write order") {
    const auto& p = read_ev.particles();
    REQUIRE(p.size() == 6);
    struct Expect { int pdg; int status; double mass; };
    const Expect want[6] = {
        {11, 4, kElectronMass},         // beam e
        {1000030060, 4, kLi6Mass},      // beam ion
        {11, 1, kElectronMass},         // e'
        {22, 3, ev.particles[3].mass},  // virtual photon
        {1000020040, 1, kAlphaMass},    // alpha spectator
        {92, 3, ev.particles[5].mass},  // X (forced status 3)
    };
    for (int i = 0; i < 6; ++i) {
      CAPTURE(i);
      CHECK(p[static_cast<std::size_t>(i)]->pid() == want[i].pdg);
      CHECK(p[static_cast<std::size_t>(i)]->status() == want[i].status);
      CHECK(p[static_cast<std::size_t>(i)]->generated_mass() == doctest::Approx(want[i].mass).epsilon(kTol));
    }
  }

  SUBCASE("4-momentum conservation: status-1 finals + the status-3 X vs the beams") {
    const auto& p = read_ev.particles();
    const auto sum = [](const HepMC3::FourVector& a, const HepMC3::FourVector& b) {
      return HepMC3::FourVector(a.px() + b.px(), a.py() + b.py(), a.pz() + b.pz(), a.e() + b.e());
    };
    HepMC3::FourVector beams = sum(p[0]->momentum(), p[1]->momentum());
    // finals (e', spectator) plus the status-3 X, which is not otherwise
    // represented by any status-1 particle (unlike the virtual photon, whose
    // momentum is already implied by beam_e - e' and must NOT be re-added).
    HepMC3::FourVector out = sum(sum(p[2]->momentum(), p[4]->momentum()), p[5]->momentum());
    CHECK(out.e() == doctest::Approx(beams.e()).epsilon(kTol));
    CHECK(out.px() == doctest::Approx(beams.px()).epsilon(kTol));
    CHECK(out.py() == doctest::Approx(beams.py()).epsilon(kTol));
    CHECK(out.pz() == doctest::Approx(beams.pz()).epsilon(kTol));
  }

  SUBCASE("event attributes round-trip") {
    auto dbl = [&](const char* name) { return read_ev.attribute<DoubleAttribute>(name)->value(); };
    auto itg = [&](const char* name) { return read_ev.attribute<IntAttribute>(name)->value(); };
    auto str = [&](const char* name) { return read_ev.attribute<StringAttribute>(name)->value(); };

    CHECK(dbl("spin_J") == doctest::Approx(ev.spin.j));
    CHECK(dbl("spin_M") == doctest::Approx(ev.spin.m_ion));
    CHECK(dbl("struck_cluster_m") == doctest::Approx(ev.spin.m_struck));
    CHECK(itg("lam_e") == ev.spin.lam_e);
    CHECK(dbl("P_e") == doctest::Approx(ev.spin.pe));
    CHECK(dbl("P_z") == doctest::Approx(ev.spin.pz));
    CHECK(dbl("P_zz") == doctest::Approx(ev.spin.pzz));
    CHECK(dbl("spin_axis_theta") == doctest::Approx(ev.spin.theta_s));
    CHECK(dbl("spin_axis_phi") == doctest::Approx(ev.spin.phi_s));
    CHECK(str("spin_category") == ev.spin.category);
    CHECK(itg("run") == ev.spin.run);
    CHECK(itg("bunch") == ev.spin.bunch);
    CHECK(str("channel") == std::string("TaggedLi6Alpha"));

    CHECK(dbl("dis_x") == doctest::Approx(ev.kin.x));
    CHECK(dbl("dis_Q2") == doctest::Approx(ev.kin.q2));
    CHECK(dbl("dis_y") == doctest::Approx(ev.kin.y));
    CHECK(dbl("dis_phi") == doctest::Approx(ev.kin.phi));
    CHECK(dbl("spectator_k") == doctest::Approx(ev.kin.k));
    CHECK(dbl("spectator_cos_theta") == doctest::Approx(ev.kin.cos_theta_k));
    CHECK(dbl("spectator_phi") == doctest::Approx(ev.kin.phi_k));
    CHECK(dbl("alpha_s") == doctest::Approx(ev.kin.alpha_s));
    CHECK(dbl("pt_s") == doctest::Approx(ev.kin.pt_s));
    CHECK(dbl("t") == doctest::Approx(ev.kin.t));
    CHECK(dbl("x_pom") == doctest::Approx(ev.kin.x_pom));
  }

  SUBCASE("weights and weight names round-trip") {
    REQUIRE(read_ev.weights().size() == 1 + ev.spin_weights.size());
    CHECK(read_ev.weights()[0] == doctest::Approx(ev.weight));
    for (std::size_t i = 0; i < ev.spin_weights.size(); ++i) {
      CHECK(read_ev.weights()[i + 1] == doctest::Approx(ev.spin_weights[i]));
    }
    REQUIRE(read_ev.run_info() != nullptr);
    const auto& names = read_ev.run_info()->weight_names();
    REQUIRE(names.size() == 1 + ev.spin_weights.size());
    CHECK(names[0] == "nominal");
    for (std::size_t i = 0; i < ev.spin_weights.size(); ++i) {
      CHECK(names[i + 1] == "spin_weight_" + std::to_string(i + 1));
    }
  }

  SUBCASE("cross section round-trip") {
    auto cs = read_ev.attribute<HepMC3::GenCrossSection>("GenCrossSection");
    REQUIRE(cs != nullptr);
    // GenCrossSection::from_string pads the vector out to the event's weight
    // count (5 here: nominal + 4 spin weights) on read-back, duplicating
    // index 0 into the extra slots; only index 0 is LiPolGen's cross section.
    REQUIRE_FALSE(cs->xsecs().empty());
    CHECK(cs->xsecs()[0] == doctest::Approx(ev.xsec_pb));
    CHECK(cs->xsec_errs()[0] == doctest::Approx(ev.xsec_err_pb));
  }

  SUBCASE("generator tool info") {
    REQUIRE(read_ev.run_info() != nullptr);
    REQUIRE(read_ev.run_info()->tools().size() == 1);
    CHECK(read_ev.run_info()->tools()[0].name == "LiPolGen");
    CHECK(read_ev.run_info()->tools()[0].version == "0.1.0");
  }

  SUBCASE("Particle::pol -> \"pol\" attribute, only when != 9") {
    const auto& p = read_ev.particles();
    auto escat_pol = p[2]->attribute<DoubleAttribute>("pol");
    REQUIRE(escat_pol != nullptr);
    CHECK(escat_pol->value() == doctest::Approx(-1.0));
    CHECK(p[1]->attribute<DoubleAttribute>("pol") == nullptr);  // beam ion: pol == 9 (unknown)
  }
}

TEST_CASE("HepMC3Writer rejects the unimplemented HepMC2 ascii format") {
  const auto path = temp_hepmc_path("lipolgen_test_hepmc_unsupported.hepmc");
  CHECK_THROWS_AS(HepMC3Writer(path.string(), HepMC3Format::HepMC2Ascii), std::runtime_error);
}

TEST_CASE("HepMC3Writer::close is idempotent and write-after-close throws") {
  const auto path = temp_hepmc_path("lipolgen_test_hepmc_close.hepmc");
  HepMC3Writer writer(path.string());
  writer.write(make_synthetic_event());
  writer.close();
  CHECK_NOTHROW(writer.close());
  CHECK_THROWS_AS(writer.write(make_synthetic_event()), std::runtime_error);
  std::filesystem::remove(path);
}

TEST_CASE("HepMC3Writer writes the PDG mass as the generated mass of a massless electron") {
  constexpr double kPdgElectronMass = 0.51099895e-3;

  const Event ev = make_electron_mass_rule_event();
  const auto path = temp_hepmc_path("lipolgen_test_hepmc_electron_mass_rule.hepmc");

  {
    HepMC3Writer writer(path.string());
    writer.write(ev);
    writer.close();
  }

  HepMC3::ReaderAscii reader(path.string());
  HepMC3::GenEvent read_ev(HepMC3::Units::GEV, HepMC3::Units::MM);
  REQUIRE(reader.read_event(read_ev));
  REQUIRE_FALSE(reader.failed());
  reader.close();
  std::filesystem::remove(path);

  const auto& p = read_ev.particles();
  REQUIRE(p.size() == ev.particles.size());
  REQUIRE(p.size() == 6);

  SUBCASE("massless electrons and positrons (Particle::mass == 0.0) are written with generated_mass == PDG electron mass") {
    // indices 0 (beam e), 2 (scattered e) and 5 (positron, pdg -11), all
    // built with mass == 0.0: the rule is keyed on |pdg| == 11.
    for (std::size_t i : {std::size_t{0}, std::size_t{2}, std::size_t{5}}) {
      CAPTURE(i);
      CHECK(std::abs(p[i]->pid()) == 11);
      CHECK(ev.particles[i].mass == 0.0);
      CHECK(p[i]->generated_mass() ==
            doctest::Approx(kPdgElectronMass).epsilon(1e-12));
    }
  }

  SUBCASE("a massive electron (Particle::mass == 0.000511) keeps its own mass, untouched by the rule") {
    CHECK(std::abs(p[3]->pid()) == 11);
    CHECK(ev.particles[3].mass == doctest::Approx(kElectronMass));
    CHECK(p[3]->generated_mass() == doctest::Approx(kElectronMass).epsilon(kTol));
  }

  SUBCASE("a non-electron with Particle::mass == 0.0 (photon) keeps mass 0, not the PDG rewrite") {
    CHECK(p[4]->pid() == 22);
    CHECK(ev.particles[4].mass == 0.0);
    CHECK(p[4]->generated_mass() == doctest::Approx(0.0).epsilon(kTol));
  }

  SUBCASE("the mass rule does not touch any particle's written four-momentum (E, p bit-identical)") {
    // Exact `==` (not Approx) is deliberate: HepMC3::WriterAscii's default
    // precision (16 digits after the point in %e format) round-trips an IEEE
    // double exactly, so any difference here would be a write-side change.
    for (std::size_t i = 0; i < ev.particles.size(); ++i) {
      CAPTURE(i);
      const HepMC3::FourVector& mom = p[i]->momentum();
      CHECK(mom.px() == ev.particles[i].p.px);
      CHECK(mom.py() == ev.particles[i].p.py);
      CHECK(mom.pz() == ev.particles[i].p.pz);
      CHECK(mom.e() == ev.particles[i].p.e);
    }
  }
}

// --------------------------------------------------------------- rc.hpp

TEST_CASE("T15: the rc.hpp weight block is APPENDED, and round-trips") {
  // A SEPARATE fixture event on purpose: `make_synthetic_event()` is shared
  // with the round-trip case above, whose two `weights().size()` REQUIREs are
  // written against `1 + spin_weights.size()`.  Adding `rc_weights` there
  // would turn them red for no physics reason (design_C_tensor_rc.md T15).
  Event ev = make_synthetic_event();
  REQUIRE(ev.spin_weights.size() == 4);
  // Slot 0 (the event's own pure state) + one slot per spin category, in the
  // row-major (n_slot x kRcWeightCount) order `Event::rc_weights` documents.
  ev.rc_weights = {0.990, 1.010, 1.0300,     // slot 0
                   0.991, 1.009, 1.0301,     // category 1
                   0.992, 1.008, 1.0302,     // category 2
                   0.993, 1.007, 1.0303,     // category 3
                   0.994, 1.006, 1.0304};    // category 4
  REQUIRE(ev.rc_weights.size() == kRcWeightCount * 5);

  const auto path = temp_hepmc_path("lipolgen_test_hepmc_rc.hepmc");
  {
    HepMC3Writer writer(path.string());
    writer.write(ev);
    writer.close();
  }
  HepMC3::ReaderAscii reader(path.string());
  HepMC3::GenEvent read_ev(HepMC3::Units::GEV, HepMC3::Units::MM);
  REQUIRE(reader.read_event(read_ev));
  REQUIRE_FALSE(reader.failed());
  reader.close();
  std::filesystem::remove(path);

  REQUIRE(read_ev.run_info() != nullptr);
  const auto& names = read_ev.run_info()->weight_names();
  const std::vector<std::string> want = {
      "nominal",       "spin_weight_1", "spin_weight_2", "spin_weight_3",
      "spin_weight_4", "rc_tensor_lo",  "rc_tensor_hi",  "rc_tail",
      "rc_tensor_lo_1", "rc_tensor_hi_1", "rc_tail_1",
      "rc_tensor_lo_2", "rc_tensor_hi_2", "rc_tail_2",
      "rc_tensor_lo_3", "rc_tensor_hi_3", "rc_tail_3",
      "rc_tensor_lo_4", "rc_tensor_hi_4", "rc_tail_4"};
  REQUIRE(names.size() == want.size());
  for (std::size_t i = 0; i < want.size(); ++i) CHECK(names[i] == want[i]);

  // "nominal" keeps index 0 and every spin_weight_k keeps its index: the RC
  // block is APPENDED, never inserted (docs/HEPMC3_CONVENTION.md).
  REQUIRE(read_ev.weights().size() ==
          1 + ev.spin_weights.size() + ev.rc_weights.size());
  CHECK(read_ev.weights()[0] == doctest::Approx(ev.weight).epsilon(1e-12));
  for (std::size_t i = 0; i < ev.spin_weights.size(); ++i) {
    CHECK(read_ev.weights()[i + 1] ==
          doctest::Approx(ev.spin_weights[i]).epsilon(1e-12));
  }
  const std::size_t off = 1 + ev.spin_weights.size();
  for (std::size_t i = 0; i < ev.rc_weights.size(); ++i) {
    CHECK(read_ev.weights()[off + i] ==
          doctest::Approx(ev.rc_weights[i]).epsilon(1e-12));
  }
}

TEST_CASE("T15b: --rc off writes exactly today's weight names") {
  // The bit-for-bit half: an event whose `rc_weights` is empty (which is what
  // `Event::reset` leaves and what an `--rc off` run produces) must name
  // nothing new at all.
  const Event ev = make_synthetic_event();
  REQUIRE(ev.rc_weights.empty());
  const auto path = temp_hepmc_path("lipolgen_test_hepmc_rc_off.hepmc");
  {
    HepMC3Writer writer(path.string());
    writer.write(ev);
    writer.close();
  }
  HepMC3::ReaderAscii reader(path.string());
  HepMC3::GenEvent read_ev(HepMC3::Units::GEV, HepMC3::Units::MM);
  REQUIRE(reader.read_event(read_ev));
  reader.close();
  std::filesystem::remove(path);
  REQUIRE(read_ev.run_info() != nullptr);
  CHECK(read_ev.run_info()->weight_names().size() ==
        1 + ev.spin_weights.size());
  CHECK(read_ev.weights().size() == 1 + ev.spin_weights.size());
}
