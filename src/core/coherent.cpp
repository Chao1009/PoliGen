#include "lipolgen/coherent.hpp"

#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "lipolgen/constants.hpp"

namespace lipolgen {
namespace {

const double kNaN = std::nan("");

std::string beam_ion_name(int beam_a, int beam_z) {
  if (beam_z == 3 && beam_a == 6) return "6Li";
  if (beam_z == 3 && beam_a == 7) return "7Li";
  if (beam_z == 2 && beam_a == 3) return "3He";
  if (beam_z == 1 && beam_a == 2) return "d";
  if (beam_z == 1 && beam_a == 1) return "p";
  throw std::runtime_error("no beam species for that (A, Z)");
}

/// Solve for the recoil four-vector that makes (R - P_R)^2 exactly m_x2 with
/// P_R on shell at m_a and transverse momentum (bx, by).
///
/// Writing P_R^+ = w and P_R^- = (m_a^2 + b^2)/w, the mass condition
///     2 R.P_R = R^+ P_R^- + R^- P_R^+ - 2 R_T.b_T = R^2 + m_a^2 - m_x2
/// is a plain QUADRATIC in w,
///     R^- w^2 - D w + R^+ (m_a^2 + b^2) = 0,
///     D = R^2 + m_a^2 - m_x2 + 2 R_T.b_T,
/// whose LARGE root is the forward recoil.  No small-x_P expansion enters, so
/// the balance k + P_ion = k' + P_R + X closes to rounding with M_X^2 = m_x2
/// BY CONSTRUCTION.  False when there is no real forward root -- the
/// diffractive mass does not fit in the available energy.
bool solve_recoil(const Vec4& r, double m_a, double bx, double by, double m_x2,
                  Vec4* out) {
  const double b2 = bx * bx + by * by;
  const double m2 = m_a * m_a;
  const double rp = r.e + r.pz, rm = r.e - r.pz;
  if (!(rm > 0.0) || !(rp > 0.0)) return false;
  const double d = (r.m2() + m2 - m_x2) + 2.0 * (r.px * bx + r.py * by);
  const double disc = d * d - 4.0 * rm * rp * (m2 + b2);
  if (!(disc >= 0.0)) return false;
  const double w = (d + std::sqrt(disc)) / (2.0 * rm);
  if (!(w > 0.0)) return false;
  const double minus = (m2 + b2) / w;
  *out = Vec4{0.5 * (w + minus), bx, by, 0.5 * (w - minus)};
  return true;
}

}  // namespace

double gaussian_slope(double r_rms_fm) {
  const double r = r_rms_fm / GEV_PER_FM_INV;
  return r * r / 3.0;
}

const std::vector<MantysaariRow>& mantysaari_a2_deuteron() {
  static const std::vector<MantysaariRow> t = {{0.05, +0.08, -0.04},
                                               {0.10, +0.15, -0.08},
                                               {0.20, +0.30, -0.17},
                                               {0.30, +0.43, -0.28}};
  return t;
}

double CoherentScenario::coherent_fraction(double x) const {
  const double r = x / x_coh;
  return f0 / (1.0 + r * r);
}

double CoherentScenario::sample_t(Rng& rng, double t_max) const {
  const double u = rng.uniform();
  const double c = 1.0 - std::exp(-slope_b * t_max);
  return -std::log(1.0 - c * u) / slope_b;
}

double CoherentScenario::dsigma_dt(double t_abs) const {
  return slope_b * std::exp(-slope_b * t_abs);
}

double CoherentScenario::tag_acceptance(double pt_cut) const {
  return std::exp(-slope_b * pt_cut * pt_cut);
}

double CoherentScenario::mean_t_tagged(double pt_cut) const {
  return pt_cut * pt_cut + 1.0 / slope_b;
}

double CoherentScenario::tag_acceptance_angular(double sigma_theta,
                                                double p_per_nucleon,
                                                int a_beam,
                                                double n_sigma) const {
  const double cut = n_sigma * sigma_theta * a_beam * p_per_nucleon;
  return std::exp(-slope_b * cut * cut);
}

double CoherentScenario::a2_deformation(double t_abs, double pzz) const {
  return -(pzz / 4.0) * eps_b0 * slope_b * t_abs;
}

double CoherentScenario::cos2phi_coefficient_deformation(double t_abs,
                                                         double pzz) const {
  return 2.0 * a2_deformation(t_abs, pzz);
}

double CoherentScenario::a2_tagged(double pt_cut, double pzz) const {
  return a2_deformation(mean_t_tagged(pt_cut), pzz);
}

double CoherentScenario::a2_m_state(double t_abs, int m) const {
  // The pure-state coefficients that reproduce the population average
  // a_2 = -(P_zz/4) eps_b0 B |t| with P_zz = p+ + p- - 2 p0, i.e.
  // a_2(+-1) = -(1/4) eps_b0 B |t| and a_2(0) = -2 a_2(+-1).
  const double a1 = -0.25 * eps_b0 * slope_b * t_abs;
  return m == 0 ? -2.0 * a1 : a1;
}

double CoherentScenario::cos2phi_coefficient(double t_abs, double pzz) const {
  return cos2phi_coefficient_deformation(t_abs, pzz) + amp * pzz;
}

double CoherentScenario::positivity_margin(double t_max, double pzz) const {
  return 1.0 - std::fabs(cos2phi_coefficient(t_max, pzz));
}

CoherentRecoil recoil_lab(double t_abs, double phi_t, double p_per_nucleon,
                          double x_pom, int a_beam) {
  const double pt = std::sqrt(t_abs);
  const double p_beam = a_beam * p_per_nucleon;
  const double pz = (1.0 - x_pom) * p_beam;
  const double p_lab = std::hypot(pt, pz);
  CoherentRecoil out;
  out.pT = pt;
  out.theta = std::atan2(pt, pz);
  out.R = p_lab / p_beam;   // same Z as the beam: rigidity ratio = momentum ratio
  out.xL = pz / p_beam;
  out.phi_t = phi_t;
  return out;
}

CoherentRecoil recoil_lab_of(const Vec4& p_recoil, double phi_t,
                             double p_per_nucleon, int a_beam) {
  const double p_beam = a_beam * p_per_nucleon;
  const double pt = p_recoil.pt();
  const double p_lab = std::sqrt(pt * pt + p_recoil.pz * p_recoil.pz);
  CoherentRecoil out;
  out.pT = pt;
  out.theta = std::atan2(pt, p_recoil.pz);
  out.R = p_lab / p_beam;   // same Z as the beam: rigidity ratio = |p| ratio
  out.xL = p_recoil.pz / p_beam;
  out.phi_t = phi_t;
  return out;
}

// ------------------------------------------------------------- the pomeron

double CoherentXpomModel::x_pom_of(double m_x2, double q2, double w2) {
  return (m_x2 + q2) / (w2 + q2);
}

double CoherentXpomModel::m_x2_of(double x_pom, double q2, double w2) {
  return x_pom * (w2 + q2) - q2;
}

double CoherentXpomModel::x_pom_min(double q2, double w2) const {
  return x_pom_of(m_x_min * m_x_min, q2, w2);
}

double CoherentXpomModel::draw(double q2, double w2, Rng& rng) const {
  const double lo = x_pom_min(q2, w2);
  const double u = rng.uniform();   // consumed on BOTH branches
  if (!(lo > 0.0)) throw std::runtime_error("CoherentXpomModel: x_P,min <= 0");
  if (!(x_pom_max > lo)) return lo;
  return lo * std::pow(x_pom_max / lo, u);
}

double fragment_rigidity(int a, int z, int beam_a, int beam_z) {
  if (z == 0) return kNaN;
  return (nuclear_mass(z, a) / z) / (nuclear_mass(beam_z, beam_a) / beam_z);
}

std::string fragment_route_label(int a, int z, int beam_a, int beam_z,
                                 const std::string& config) {
  if (z == 0) return "ZDC";
  const double r = fragment_rigidity(a, z, beam_a, beam_z);
  if (std::fabs(r - 1.0) < NEAR_BEAM_BAND) return "RP pT-tail only (beam-blind)";
  if (RP_R_WINDOW_LO <= r && r <= RP_R_WINDOW_HI) return "RomanPots";
  if (OMD_R_WINDOW_LO <= r && r < OMD_R_WINDOW_HI) return "OMD";
  if (r > 1.0 + NEAR_BEAM_BAND) {
    return over_rigid_route(r, 0.0, config) ? "RP-inner (over-rigid)"
                                            : "lost (over-rigid)";
  }
  return "lost";
}

const std::vector<BreakupChannel>& breakup_table(int beam_a, int beam_z) {
  // Lowest-lying particle decompositions with separation energies [MeV].
  // 6Li: TUNL A=6 (Tilley et al. NPA 708:3).  7Li: recomputed from
  // `nuclear_mass`, which reproduces the evaluated values.
  static const std::vector<BreakupChannel> li6 = {
      {"alpha+d", 1.474, {{"alpha", 4, 2}, {"d", 2, 1}}},
      {"alpha+p+n", 3.70, {{"alpha", 4, 2}, {"p", 1, 1}, {"n", 1, 0}}},
      {"3He+t", 15.79, {{"3He", 3, 2}, {"t", 3, 1}}}};
  static const std::vector<BreakupChannel> li7 = {
      {"alpha+t", 2.468, {{"alpha", 4, 2}, {"t", 3, 1}}},
      {"6Li+n", 7.251, {{"6Li", 6, 3}, {"n", 1, 0}}},
      {"alpha+d+n", 8.725, {{"alpha", 4, 2}, {"d", 2, 1}, {"n", 1, 0}}},
      {"6He+p", 9.975, {{"6He", 6, 2}, {"p", 1, 1}}}};
  if (beam_z == 3 && beam_a == 6) return li6;
  if (beam_z == 3 && beam_a == 7) return li7;
  throw std::runtime_error("no breakup table for that (A, Z)");
}

std::vector<VetoRow> veto_table(int beam_a, int beam_z,
                                const std::string& config) {
  std::vector<VetoRow> out;
  for (const BreakupChannel& ch : breakup_table(beam_a, beam_z)) {
    for (const BreakupFragment& f : ch.fragments) {
      VetoRow row;
      row.channel = ch.name;
      row.fragment = f.name;
      row.rigidity = fragment_rigidity(f.a, f.z, beam_a, beam_z);
      row.destination = fragment_route_label(f.a, f.z, beam_a, beam_z, config);
      out.push_back(row);
    }
  }
  return out;
}

CoherentSampler::CoherentSampler(const CoherentScenario& scenario,
                                 double p_per_nucleon, double pzz,
                                 double phi_s, int beam_a, int beam_z)
    : sc_(scenario), p_u_(p_per_nucleon), pzz_(pzz), phi_s_(phi_s),
      beam_a_(beam_a), beam_z_(beam_z),
      m_beam_(nuclear_mass(beam_z, beam_a)) {
  const std::string name = beam_ion_name(beam_a, beam_z);
  optics_ = yr_optics(name, p_per_nucleon, true);
  pot_config_ = yr_config_key(name, p_per_nucleon);
  check_positivity();   // P4: at the default t_max, before any event is drawn
}

void CoherentSampler::set_t_max(double t_max) {
  const double keep = t_max_;
  t_max_ = t_max;
  try {
    check_positivity();
  } catch (...) {
    t_max_ = keep;
    throw;
  }
}

void CoherentSampler::check_positivity() const {
  if (!(t_max_ > 0.0)) {
    throw std::runtime_error("CoherentSampler: t_max must be positive");
  }
  const double margin = sc_.positivity_margin(t_max_, pzz_);
  if (margin < 0.0) {
    char buf[320];
    std::snprintf(buf, sizeof(buf),
                  "CoherentSampler: the cos 2phi weight goes NEGATIVE inside "
                  "the |t| range -- 1 - |c2| = %.6g at |t| = %.6g, P_zz = %.6g "
                  "(c2 = %.6g).  Lower t_max (the default is %g), or lower "
                  "eps_b0 / amp.",
                  margin, t_max_, pzz_, sc_.cos2phi_coefficient(t_max_, pzz_),
                  COHERENT_T_MAX_DEFAULT);
    throw std::runtime_error(buf);
  }
}

void CoherentSampler::set_optics(const Optics& optics,
                                 const std::string& pot_config) {
  optics_ = optics;
  pot_config_ = pot_config;
}

CoherentEvent CoherentSampler::sample(Rng& rng, const CoherentDis& dis) const {
  CoherentEvent ev;
  ev.t = sc_.sample_t(rng, t_max_);
  ev.c2 = sc_.cos2phi_coefficient(ev.t, pzz_);
  if (weighted_) {
    // rejection against the flat envelope 1 + |c2|
    const double env = 1.0 + std::fabs(ev.c2);
    for (;;) {
      const double phi = 2.0 * kPi * rng.uniform();
      const double w = 1.0 + ev.c2 * std::cos(2.0 * (phi - phi_s_));
      if (rng.uniform() * env <= w) {
        ev.phi_t = phi;
        break;
      }
    }
    ev.weight = 1.0;
  } else {
    ev.phi_t = 2.0 * kPi * rng.uniform();
    ev.weight = 1.0 + ev.c2 * std::cos(2.0 * (ev.phi_t - phi_s_));
  }
  // --- the pomeron, the diffractive mass, and the EXACT recoil -----------
  ev.x_pom = xpom_.draw(dis.q2, dis.w2, rng);
  ev.m_x2 = CoherentXpomModel::m_x2_of(ev.x_pom, dis.q2, dis.w2);
  const double pt = std::sqrt(ev.t);
  const double bx = pt * std::cos(ev.phi_t);
  const double by = pt * std::sin(ev.phi_t);
  if (!solve_recoil(dis.residual, m_beam_, bx, by, ev.m_x2, &ev.p_recoil)) {
    // No real forward root at this M_X.  The discriminant is largest at the
    // SMALLEST diffractive mass, so retry there before giving up; it consumes
    // no random number, so the event stream is unchanged.
    ev.x_pom = xpom_.x_pom_min(dis.q2, dis.w2);
    ev.m_x2 = CoherentXpomModel::m_x2_of(ev.x_pom, dis.q2, dis.w2);
    if (!solve_recoil(dis.residual, m_beam_, bx, by, ev.m_x2, &ev.p_recoil)) {
      throw std::runtime_error(
          "CoherentSampler: no on-shell intact recoil exists at this (x, Q2) "
          "even at M_X,min -- the diffractive mass does not fit");
    }
  }
  ev.beta = dis.x / ev.x_pom;
  ev.recoil = recoil_lab_of(ev.p_recoil, ev.phi_t, p_u_, beam_a_);
  // the EXACT Mandelstam t of the constructed recoil, |t| = |t_min| + pT^2
  const double p_beam = beam_a_ * p_u_;
  const Vec4 p_ion{std::sqrt(p_beam * p_beam + m_beam_ * m_beam_), 0.0, 0.0,
                   p_beam};
  ev.t_exact = -(p_ion - ev.p_recoil).m2();
  ev.route = route_charged(ev.recoil.R, ev.recoil.theta, ev.recoil.pT, optics_,
                           ev.phi_t, kNaN, pot_config_);
  return ev;
}

void CoherentSampler::fill_event(Event& ev, const CoherentEvent& ce) const {
  int mother = -1;
  for (std::size_t i = 0; i < ev.particles.size(); ++i) {
    if (ev.particles[i].role == Role::BeamIon) {
      mother = static_cast<int>(i);
      break;
    }
  }
  Particle p;
  p.pdg = 1000000000 + beam_z_ * 10000 + beam_a_ * 10;
  p.status = Status::Final;
  p.role = Role::IntactRecoil;
  p.p = ce.p_recoil;
  p.mass = m_beam_;
  p.charge = beam_z_;
  p.mother1 = mother;
  ev.particles.push_back(p);
  ev.kin.t = ce.t;
  ev.kin.x_pom = ce.x_pom;
  ev.kin.beta_pom = ce.beta;
  ev.kin.m_x2 = ce.m_x2;
  ev.channel = Channel::CoherentLi6;
  ev.weight *= ce.weight;
}

}  // namespace lipolgen
