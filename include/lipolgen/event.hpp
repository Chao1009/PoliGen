// LiPolGen event record (tiers T0/T1; T2 particles are added by the PYTHIA tier).
#pragma once
#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace lipolgen {

struct Vec4 {                       // (E, px, py, pz) in GeV, head-on frame
  double e = 0, px = 0, py = 0, pz = 0;
  double m2() const { return e * e - px * px - py * py - pz * pz; }
  double pt() const;
  double p() const;
  Vec4 operator+(const Vec4& o) const { return {e + o.e, px + o.px, py + o.py, pz + o.pz}; }
  Vec4 operator-(const Vec4& o) const { return {e - o.e, px - o.px, py - o.py, pz - o.pz}; }
};

// HepMC3-compatible status codes.
enum class Status : int {
  Beam = 4,          // incoming beam particle
  Final = 1,         // final-state particle
  Decayed = 2,       // decayed / intermediate physical particle
  Intermediate = 3,  // documentation (virtual photon, struck cluster, pseudo-particle X)
  Spectator = 1      // final-state spectator fragment (same as Final; role carried by Particle::role)
};

enum class Role : std::uint8_t {
  BeamElectron, BeamIon, ScatteredElectron, VirtualPhoton,
  StruckCluster, StruckNucleon, Spectator, PartnerSpectator,
  HadronicX,        // T0 pseudo-particle carrying the whole hadronic final state
  Hadron,           // T2 PYTHIA final-state particle
  IntactRecoil,     // coherent channel: the ground-state nucleus
  Other
};

struct Particle {
  int pdg = 0;              // 10-digit ion codes for nuclei (1000030060 = 6Li)
  Status status = Status::Final;
  Role role = Role::Other;
  Vec4 p;
  double mass = 0;          // generated mass (GeV)
  double charge = 0;        // in units of e
  int mother1 = -1, mother2 = -1;   // indices into Event::particles, -1 = none
  double pol = 9.0;         // helicity label as in PYTHIA/LHEF SPINUP; 9 = unknown
};

// Spin/run labels carried by every event (HepMC3 attributes, plans/04 #17 convention).
struct SpinLabels {
  double j = 0;             // ion spin J (1 or 1.5)
  double m_ion = 0;         // ion projection M on the quantization axis
  double m_struck = 0;      // struck-cluster projection m_S (tagged mode), NaN if inclusive
  int lam_e = 0;            // electron helicity ±1 (0 = unpolarized)
  double pe = 0;            // electron polarization magnitude
  double theta_s = 0, phi_s = 0;   // quantization axis in the head-on frame
  double pz = 0, pzz = 0;   // fill vector / tensor polarization
  std::string category;     // SpinCategory name
  int run = 0, bunch = 0;   // bookkeeping stream identifiers
};

struct Kinematics {
  double x = 0, q2 = 0, y = 0, phi = 0;   // per-nucleon DIS variables, φ = lab azimuth of e' about q... see docs
  double w2 = 0, nu = 0;
  double s = 0;                           // per-nucleon s
  /// Accepted-cell index of the `InclusiveSampler` the (x, Q2) was drawn
  /// from: the ION-level sampler on the inclusive channel, the STRUCK-CLUSTER
  /// sampler on the tagged ones, the f_coh-reweighted one on the coherent
  /// channel, -1 on a record not built by `Pipeline`.  This is the `cell`
  /// column
  /// `polligen.sample.weights_for` (Mode W) reweights on, so it has to travel
  /// with the event rather than staying an `InclusiveSampler` internal.
  int cell = -1;
  // spectator / tagging (tagged mode; zeros otherwise)
  double k = 0, cos_theta_k = 0, phi_k = 0;      // spectator momentum in the spin frame
  double alpha_s = 0, pt_s = 0;                  // light-front α_s and p_T of the spectator
  /// The spectator's LAB kinematics -- `tagged.hpp`'s `SpectatorLab`, stored
  /// on the record instead of being re-derived from the four-vector by every
  /// consumer.  (kx, ky, kz) are the rest-frame components in the LAB
  /// orientation, `spec_r` the rigidity ratio against the beam (NaN for a
  /// neutral fragment), `spec_xl` the longitudinal momentum fraction and
  /// `phi_spec` the lab azimuth the planar Roman-Pot rectangle needs --
  /// deliberately not called `phi`, which is the DIS azimuth.
  /// Zero / NaN on a channel with no spectator.
  double spec_pt = 0, spec_theta = 0, spec_p_lab = 0;
  double spec_r = std::numeric_limits<double>::quiet_NaN();
  double spec_xl = std::numeric_limits<double>::quiet_NaN();
  double spec_kx = 0, spec_ky = 0, spec_kz = 0, phi_spec = 0;
  // coherent channel
  double t = 0;          // |t| the exponential slope was drawn at = pT_recoil^2
  double x_pom = 0;      // PER-NUCLEON pomeron fraction (coherent.hpp)
  double beta_pom = 0;   // x / x_pom = Q^2/(M_X^2 + Q^2), in (0, 1]
  double m_x2 = 0;       // M_X^2 of the diffractive system [GeV^2]
};

enum class Channel : std::uint8_t { Inclusive, TaggedLi6Alpha, TaggedLi6D, TaggedLi7Alpha, TaggedLi7T,
                                    TaggedDeuteronP, TaggedDeuteronN, TaggedHe3P, CoherentLi6 };

struct Event {
  std::uint64_t number = 0;
  Channel channel = Channel::Inclusive;
  double weight = 1.0;                 // 1 for unweighted
  std::vector<double> spin_weights;    // weighted mode: one weight per spin category (may be empty)
  double xsec_pb = 0, xsec_err_pb = 0; // running estimate of the generated cross section
  SpinLabels spin;
  Kinematics kin;
  std::vector<Particle> particles;     // [0]=beam e, [1]=beam ion, then the rest

  const Particle* find(Role r) const;
  Vec4 total_final() const;            // Σ over Status::Final
  double total_charge_final() const;

  /// Back to the default-constructed state, KEEPING the heap capacity of
  /// `particles`, `spin_weights` and `spin.category`.  This is what lets a
  /// generator reconstruct into an existing record instead of building a
  /// fresh one and move-assigning it: a move-assignment frees the target's
  /// buffers and steals the temporary's, so the allocator sees one
  /// free/malloc pair per event either way, while `reset()` sees none.
  void reset();
};

}  // namespace lipolgen
