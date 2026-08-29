// Writes 10 synthetic e+6Li tagged-mode events to a HepMC3 Asciiv3 file, to
// exercise lipolgen::HepMC3Writer end-to-end (including from outside the
// HepMC3 headers themselves -- this file never includes them).
//
// Usage: write_hepmc_example [output.hepmc]
//
// Verify externally with pyhepmc, e.g.:
//   python3 -c "
//   import pyhepmc
//   with pyhepmc.open('write_hepmc_example.hepmc') as f:
//       events = list(f)
//   print(len(events), events[0].particles[0].pid)
//   "
#include "lipolgen/event.hpp"
#include "lipolgen/hepmc_writer.hpp"

#include <cmath>
#include <iostream>
#include <string>

using namespace lipolgen;

namespace {

constexpr double kLi6Mass = 5.601518702;    // AME2020, atomic mass less electrons
constexpr double kAlphaMass = 3.727379407;  // ditto
constexpr double kElectronMass = 0.000511;

// One synthetic e + 6Li -> e' + alpha(spectator) + X tagged event, with the
// scattered-electron and spectator kinematics nudged by `wobble` (in [0,1))
// so the ten output events are not bit-identical copies of each other.
Event make_event(std::uint64_t number, double wobble) {
  Event ev;
  ev.number = number;
  ev.channel = Channel::TaggedLi6Alpha;
  ev.weight = 1.0;
  ev.spin_weights = {1.0 + 0.05 * wobble, 1.0 - 0.05 * wobble};
  ev.xsec_pb = 100.0 + wobble;
  ev.xsec_err_pb = 0.5;

  ev.spin.j = 1.0;
  ev.spin.m_ion = (number % 3 == 0) ? 1.0 : ((number % 3 == 1) ? 0.0 : -1.0);
  ev.spin.m_struck = 0.5;
  ev.spin.lam_e = (number % 2 == 0) ? 1 : -1;
  ev.spin.pe = 0.8;
  ev.spin.theta_s = 0.0;
  ev.spin.phi_s = 0.0;
  ev.spin.pz = 8.0 / 13.0;
  ev.spin.pzz = 4.0 / 13.0;
  ev.spin.category = "pp";
  ev.spin.run = 1;
  ev.spin.bunch = static_cast<int>(number % 4);

  ev.kin.x = 0.2 + 0.01 * wobble;
  ev.kin.q2 = 4.0 + 0.1 * wobble;
  ev.kin.y = 0.55;
  ev.kin.phi = 0.1 * wobble;
  ev.kin.k = 0.31;
  ev.kin.cos_theta_k = -0.2;
  ev.kin.phi_k = 0.0;
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

  const double p_ion = 6.0 * 99.5;
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
    const double px = 1.5 + 0.1 * wobble, py = 0.5, pz = -7.8 - 0.1 * wobble;
    escat.p = {std::sqrt(px * px + py * py + pz * pz + kElectronMass * kElectronMass), px, py, pz};
  }
  escat.mass = kElectronMass;
  escat.charge = -1.0;
  escat.mother1 = 0;
  escat.pol = ev.spin.lam_e;

  Particle gamma;
  gamma.pdg = 22;
  gamma.status = Status::Intermediate;
  gamma.role = Role::VirtualPhoton;
  gamma.p = beam_e.p - escat.p;
  gamma.mass = std::sqrt(std::max(gamma.p.m2(), 0.0));
  gamma.mother1 = 0;

  Particle spectator;
  spectator.pdg = 1000020040;
  spectator.status = Status::Final;
  spectator.role = Role::Spectator;
  {
    const double px = 0.3, py = -0.2 + 0.05 * wobble, pz = 590.0;
    spectator.p = {std::sqrt(px * px + py * py + pz * pz + kAlphaMass * kAlphaMass), px, py, pz};
  }
  spectator.mass = kAlphaMass;
  spectator.charge = 2.0;
  spectator.mother1 = 1;

  Particle x;
  x.pdg = 92;
  x.status = Status::Intermediate;
  x.role = Role::HadronicX;
  x.p = (beam_e.p + beam_ion.p) - escat.p - spectator.p;
  x.mass = std::sqrt(std::max(x.p.m2(), 0.0));
  x.charge = beam_ion.charge - spectator.charge;
  x.mother1 = 3;

  ev.particles = {beam_e, beam_ion, escat, gamma, spectator, x};
  return ev;
}

}  // namespace

int main(int argc, char** argv) {
  const std::string out = (argc > 1) ? argv[1] : "write_hepmc_example.hepmc";
  try {
    HepMC3Writer writer(out);
    for (std::uint64_t i = 0; i < 10; ++i) {
      writer.write(make_event(i, static_cast<double>(i) / 10.0));
    }
    writer.close();
  } catch (const std::exception& e) {
    std::cerr << "write_hepmc_example: " << e.what() << std::endl;
    return 1;
  }
  std::cout << "wrote 10 events to " << out << std::endl;
  return 0;
}
