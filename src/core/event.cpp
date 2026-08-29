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
}  // namespace lipolgen
