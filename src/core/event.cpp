#include "lipolgen/event.hpp"
#include <cmath>
namespace lipolgen {
double Vec4::pt() const { return std::sqrt(px * px + py * py); }
double Vec4::p() const { return std::sqrt(px * px + py * py + pz * pz); }
const Particle* Event::find(Role r) const {
  for (const auto& p : particles) if (p.role == r) return &p;
  return nullptr;
}
Vec4 Event::total_final() const {
  Vec4 s;
  for (const auto& p : particles) if (p.status == Status::Final) s = s + p.p;
  return s;
}
double Event::total_charge_final() const {
  double q = 0;
  for (const auto& p : particles) if (p.status == Status::Final) q += p.charge;
  return q;
}
void Event::reset() {
  number = 0;
  channel = Channel::Inclusive;
  weight = 1.0;
  spin_weights.clear();
  rc_weights.clear();   // rc.hpp; capacity kept, like spin_weights
  rc_clipped = 0;
  xsec_pb = 0.0;
  xsec_err_pb = 0.0;
  kin = Kinematics();
  // Field by field rather than `spin = SpinLabels()`, which would free the
  // category string's buffer -- the one allocation this method exists to keep.
  spin.j = 0.0;
  spin.m_ion = 0.0;
  spin.m_struck = 0.0;
  spin.lam_e = 0;
  spin.pe = 0.0;
  spin.theta_s = 0.0;
  spin.phi_s = 0.0;
  spin.pz = 0.0;
  spin.pzz = 0.0;
  spin.category.clear();
  spin.run = 0;
  spin.bunch = 0;
  particles.clear();
}

}  // namespace lipolgen
