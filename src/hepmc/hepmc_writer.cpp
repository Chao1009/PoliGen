// SPDX-License-Identifier: GPL-3.0-or-later
#include "lipolgen/hepmc_writer.hpp"

#include "lipolgen/constants.hpp"
#include "lipolgen/rc.hpp"
#include <cstdlib>

#include <HepMC3/Attribute.h>
#include <HepMC3/GenCrossSection.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenRunInfo.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/WriterAscii.h>

#include <stdexcept>
#include <utility>

namespace lipolgen {
namespace {

HepMC3::FourVector to_fourvector(const Vec4& v) {
  // HepMC3::FourVector(x, y, z, t) — momentum-unit fields read as (px,py,pz,e).
  return HepMC3::FourVector(v.px, v.py, v.pz, v.e);
}

// HadronicX is documentation/status-3 by construction (Event::Status
// comment) and is forced here regardless of the stored Status, so it is
// never emitted as a final-state particle even if the event also carries
// real T2 hadrons (Role::Hadron) for the same hadronic system.
int hepmc_status(const Particle& p) {
  if (p.role == Role::HadronicX) return 3;
  return static_cast<int>(p.status);
}

const char* channel_name(Channel c) {
  switch (c) {
    case Channel::Inclusive: return "Inclusive";
    case Channel::TaggedLi6Alpha: return "TaggedLi6Alpha";
    case Channel::TaggedLi6D: return "TaggedLi6D";
    case Channel::TaggedLi7Alpha: return "TaggedLi7Alpha";
    case Channel::TaggedLi7T: return "TaggedLi7T";
    case Channel::TaggedDeuteronP: return "TaggedDeuteronP";
    case Channel::TaggedDeuteronN: return "TaggedDeuteronN";
    case Channel::TaggedHe3P: return "TaggedHe3P";
    case Channel::CoherentLi6: return "CoherentLi6";
  }
  return "Unknown";
}

}  // namespace

struct HepMC3Writer::Impl {
  HepMC3::WriterAscii writer;
  std::shared_ptr<HepMC3::GenRunInfo> run_info;  // built lazily, from event 0
  std::string generator_name;
  std::string generator_version;
  bool closed = false;

  explicit Impl(const std::string& filename, std::string name, std::string version)
      // No GenRunInfo passed here on purpose: WriterAscii only serializes a
      // non-null run_info immediately at construction, before we know the
      // first event's spin_weights size. Instead every GenEvent below carries
      // run_info explicitly, and WriterAscii::write_event() picks it up (and
      // writes it once) from the first event it sees.
      : writer(filename), generator_name(std::move(name)), generator_version(std::move(version)) {}

  void ensure_run_info(const Event& ev) {
    if (run_info) return;
    run_info = std::make_shared<HepMC3::GenRunInfo>();
    HepMC3::GenRunInfo::ToolInfo tool;
    tool.name = generator_name;
    tool.version = generator_version;
    tool.description = "lipolgen::HepMC3Writer (docs/HEPMC3_CONVENTION.md)";
    run_info->tools().push_back(tool);
    std::vector<std::string> names = {"nominal"};
    for (std::size_t i = 0; i < ev.spin_weights.size(); ++i) {
      names.push_back("spin_weight_" + std::to_string(i + 1));
    }
    // rc.hpp: APPENDED after the spin block so that "nominal" stays index 0
    // and every existing spin_weight_k keeps its index -- append, never
    // insert (docs/HEPMC3_CONVENTION.md).  `Event::rc_weights` is EMPTY when
    // the run has `--rc off`, so such a file is byte-identical to today's.
    // The block is row-major (n_slot x kRcWeightCount): slot 0 is the event's
    // own pure spin state, slot 1 + k is spin category k.
    const std::size_t n_slot = ev.rc_weights.size() / kRcWeightCount;
    for (std::size_t s = 0; s < n_slot; ++s) {
      for (std::size_t i = 0; i < kRcWeightCount; ++i) {
        names.push_back(rc_weight_name(i, s));
      }
    }
    run_info->set_weight_names(names);
  }
};

HepMC3Writer::HepMC3Writer(const std::string& filename, HepMC3Format format,
                            std::string generator_name, std::string generator_version)
    : impl_(nullptr) {
  if (format != HepMC3Format::Asciiv3) {
    throw std::runtime_error("HepMC3Writer: only HepMC3Format::Asciiv3 is implemented");
  }
  impl_ = std::make_unique<Impl>(filename, std::move(generator_name), std::move(generator_version));
}

HepMC3Writer::~HepMC3Writer() {
  if (impl_) close();
}

void HepMC3Writer::close() {
  if (!impl_ || impl_->closed) return;
  impl_->writer.close();
  impl_->closed = true;
}

void HepMC3Writer::write(const Event& ev) {
  if (!impl_ || impl_->closed) {
    throw std::runtime_error("HepMC3Writer::write: writer is closed");
  }
  impl_->ensure_run_info(ev);

  HepMC3::GenEvent genevt(HepMC3::Units::GEV, HepMC3::Units::MM);
  genevt.set_run_info(impl_->run_info);
  genevt.set_event_number(static_cast<int>(ev.number));

  // Event weights: index 0 = nominal, 1.. = spin-category weights, matching
  // the names registered in ensure_run_info(). Set explicitly (rather than
  // relying on GenEvent::set_run_info's default fill) so the vector is
  // correctly sized even if this event's spin_weights count differs from
  // the one that established the GenRunInfo names.
  genevt.weights().assign(1 + ev.spin_weights.size() + ev.rc_weights.size(),
                          1.0);
  genevt.weights()[0] = ev.weight;
  for (std::size_t i = 0; i < ev.spin_weights.size(); ++i) {
    genevt.weights()[i + 1] = ev.spin_weights[i];
  }
  // rc.hpp, in the same order `ensure_run_info` named them.
  const std::size_t rc_off = 1 + ev.spin_weights.size();
  for (std::size_t i = 0; i < ev.rc_weights.size(); ++i) {
    genevt.weights()[rc_off + i] = ev.rc_weights[i];
  }

  const std::size_t n = ev.particles.size();
  if (n < 2) {
    throw std::runtime_error("HepMC3Writer::write: event needs at least the two beam particles");
  }

  // Pass 1: build every GenParticle up front so mother indices may point
  // forward or backward in Event::particles without ordering assumptions.
  std::vector<HepMC3::GenParticlePtr> gp(n);
  for (std::size_t i = 0; i < n; ++i) {
    const Particle& p = ev.particles[i];
    gp[i] = std::make_shared<HepMC3::GenParticle>(to_fourvector(p.p), p.pdg, hepmc_status(p));
    // Electrons are built massless in the core (standard DIS kinematics);
    // Geant4/DD4hep assigns the PDG mass anyway and would otherwise nudge E by
    // O(10 ppm) to keep E^2 - p^2 >= m_e^2.  Write the PDG mass as the
    // generated mass so the downstream chain sees a consistent record.
    const double gen_mass = (std::abs(p.pdg) == 11 && p.mass == 0.0) ? M_ELECTRON : p.mass;
    gp[i]->set_generated_mass(gen_mass);
  }

  // Primary vertex: both beams incoming, status 4. HepMC3 requires at least
  // one vertex in the event; since it has incoming particles it is written
  // normally (unlike a bare single-particle gun event with no incoming
  // particles at its vertex, which HepMC3's own writer omits — see
  // docs/HEPMC3_CONVENTION.md).
  auto v0 = std::make_shared<HepMC3::GenVertex>();
  v0->add_particle_in(gp[0]);  // beam electron
  v0->add_particle_in(gp[1]);  // beam ion
  genevt.add_vertex(v0);

  // Pass 2: attach every other particle using Particle::mother1/mother2 when
  // set; particles with no mother hang directly off the primary vertex. A
  // daughter is added as outgoing from its mother's *end vertex* — reused if
  // the mother already has one (e.g. a mother that is a beam particle already
  // ends at v0, so its daughters land back on v0 too), created otherwise.
  // This intentionally has no separate vertex bookkeeping map: the graph
  // state already tracked by GenParticle::end_vertex() is the single source
  // of truth, so two particles sharing a mother land on the same vertex for
  // free, and a particle produced at a beam's own vertex is handled the same
  // way as one produced anywhere else.
  for (std::size_t i = 2; i < n; ++i) {
    const Particle& p = ev.particles[i];
    if (p.mother1 < 0) {
      v0->add_particle_out(gp[i]);
      continue;
    }
    HepMC3::GenParticlePtr mother1 = gp[static_cast<std::size_t>(p.mother1)];
    HepMC3::GenVertexPtr v = mother1->end_vertex();
    if (!v) {
      v = std::make_shared<HepMC3::GenVertex>();
      v->add_particle_in(mother1);
      genevt.add_vertex(v);
    }
    if (p.mother2 >= 0 && p.mother2 != p.mother1) {
      v->add_particle_in(gp[static_cast<std::size_t>(p.mother2)]);
    }
    v->add_particle_out(gp[i]);
  }

  // Particle::pol -> "pol" attribute (PYTHIA/LHEF SPINUP-style helicity
  // label), only when it carries information (9 = unknown, per event.hpp).
  for (std::size_t i = 0; i < n; ++i) {
    if (ev.particles[i].pol != 9.0) {
      gp[i]->add_attribute("pol", std::make_shared<HepMC3::DoubleAttribute>(ev.particles[i].pol));
    }
  }

  // Cross section (pb), attached per event.
  auto cs = std::make_shared<HepMC3::GenCrossSection>();
  cs->set_cross_section(std::vector<double>{ev.xsec_pb}, std::vector<double>{ev.xsec_err_pb});
  genevt.set_cross_section(cs);

  // Spin/run labels (docs/HEPMC3_CONVENTION.md).
  genevt.add_attribute("spin_J", std::make_shared<HepMC3::DoubleAttribute>(ev.spin.j));
  genevt.add_attribute("spin_M", std::make_shared<HepMC3::DoubleAttribute>(ev.spin.m_ion));
  genevt.add_attribute("struck_cluster_m", std::make_shared<HepMC3::DoubleAttribute>(ev.spin.m_struck));
  genevt.add_attribute("lam_e", std::make_shared<HepMC3::IntAttribute>(ev.spin.lam_e));
  genevt.add_attribute("P_e", std::make_shared<HepMC3::DoubleAttribute>(ev.spin.pe));
  genevt.add_attribute("P_z", std::make_shared<HepMC3::DoubleAttribute>(ev.spin.pz));
  genevt.add_attribute("P_zz", std::make_shared<HepMC3::DoubleAttribute>(ev.spin.pzz));
  genevt.add_attribute("spin_axis_theta", std::make_shared<HepMC3::DoubleAttribute>(ev.spin.theta_s));
  genevt.add_attribute("spin_axis_phi", std::make_shared<HepMC3::DoubleAttribute>(ev.spin.phi_s));
  genevt.add_attribute("spin_category", std::make_shared<HepMC3::StringAttribute>(ev.spin.category));
  genevt.add_attribute("run", std::make_shared<HepMC3::IntAttribute>(ev.spin.run));
  genevt.add_attribute("bunch", std::make_shared<HepMC3::IntAttribute>(ev.spin.bunch));
  genevt.add_attribute("channel", std::make_shared<HepMC3::StringAttribute>(channel_name(ev.channel)));

  // Kinematics.
  genevt.add_attribute("dis_x", std::make_shared<HepMC3::DoubleAttribute>(ev.kin.x));
  genevt.add_attribute("dis_Q2", std::make_shared<HepMC3::DoubleAttribute>(ev.kin.q2));
  genevt.add_attribute("dis_y", std::make_shared<HepMC3::DoubleAttribute>(ev.kin.y));
  genevt.add_attribute("dis_phi", std::make_shared<HepMC3::DoubleAttribute>(ev.kin.phi));
  genevt.add_attribute("spectator_k", std::make_shared<HepMC3::DoubleAttribute>(ev.kin.k));
  genevt.add_attribute("spectator_cos_theta", std::make_shared<HepMC3::DoubleAttribute>(ev.kin.cos_theta_k));
  genevt.add_attribute("spectator_phi", std::make_shared<HepMC3::DoubleAttribute>(ev.kin.phi_k));
  genevt.add_attribute("alpha_s", std::make_shared<HepMC3::DoubleAttribute>(ev.kin.alpha_s));
  genevt.add_attribute("pt_s", std::make_shared<HepMC3::DoubleAttribute>(ev.kin.pt_s));
  genevt.add_attribute("t", std::make_shared<HepMC3::DoubleAttribute>(ev.kin.t));
  genevt.add_attribute("x_pom", std::make_shared<HepMC3::DoubleAttribute>(ev.kin.x_pom));

  impl_->writer.write_event(genevt);
}

}  // namespace lipolgen
