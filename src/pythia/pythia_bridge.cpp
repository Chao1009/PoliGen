// SPDX-License-Identifier: GPL-3.0-or-later
// Tier T2 -- PYTHIA 8 showering and hadronization of the gamma*-nucleon
// system.  See docs/PYTHIA_BRIDGE.md for the physics and every
// approximation; the header carries the two structural facts that force
// this design.

#include "lipolgen/pythia_bridge.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>

#include "lhaup_dis.hpp"

#ifndef LIPOLGEN_PYTHIA8_XMLDOC
#define LIPOLGEN_PYTHIA8_XMLDOC ""
#endif

namespace lipolgen {

using pythia_detail::DisPayload;
using pythia_detail::LhaupDis;
using pythia_detail::RngEngine;

namespace {

// --------------------------------------------------------------- vectors

inline Vec4 from_p8(const Pythia8::Vec4& v) {
  return {v.e(), v.px(), v.py(), v.pz()};
}

struct Vec3 {
  double x = 0, y = 0, z = 0;
  double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
  double norm() const { return std::sqrt(dot(*this)); }
  Vec3 cross(const Vec3& o) const {
    return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
  }
  Vec3 minus(const Vec3& o, double c) const {
    return {x - c * o.x, y - c * o.y, z - c * o.z};
  }
};
inline Vec3 vec3(const Vec4& v) { return {v.px, v.py, v.pz}; }

/// Active boost by the three-velocity `b`: takes a vector at rest to one
/// moving with velocity `b`.
Vec4 boost_by(const Vec4& v, double bx, double by, double bz) {
  const double b2 = bx * bx + by * by + bz * bz;
  if (b2 <= 0.0) return v;
  const double gamma = 1.0 / std::sqrt(1.0 - b2);
  const double bp = bx * v.px + by * v.py + bz * v.pz;
  const double f = (gamma - 1.0) * bp / b2 + gamma * v.e;
  return {gamma * (v.e + bp), v.px + f * bx, v.py + f * by, v.pz + f * bz};
}
/// Into the rest frame of `ref`.
inline Vec4 to_rest(const Vec4& v, const Vec4& ref) {
  return boost_by(v, -ref.px / ref.e, -ref.py / ref.e, -ref.pz / ref.e);
}
/// Out of the rest frame of `ref`.
inline Vec4 from_rest(const Vec4& v, const Vec4& ref) {
  return boost_by(v, ref.px / ref.e, ref.py / ref.e, ref.pz / ref.e);
}

/// Right-handed triad in the rest frame of the hadronic system:
/// e3 along the virtual photon, e1 in the lepton plane, e2 = e3 x e1.
struct Triad {
  Vec3 e1, e2, e3;
};
Triad make_triad(const Vec4& had, const Vec4& q, const Vec4& k_in) {
  const Vec3 qr = vec3(to_rest(q, had));
  const Vec3 kr = vec3(to_rest(k_in, had));
  Triad t;
  double n = qr.norm();
  if (n <= 0.0) {
    t.e3 = {0, 0, 1};
  } else {
    t.e3 = {qr.x / n, qr.y / n, qr.z / n};
  }
  Vec3 perp = kr.minus(t.e3, kr.dot(t.e3));
  n = perp.norm();
  if (n <= 1e-12 * std::max(1.0, kr.norm())) {
    // Degenerate lepton plane (k parallel to q): any perpendicular will do.
    Vec3 seed = (std::fabs(t.e3.z) < 0.9) ? Vec3{0, 0, 1} : Vec3{1, 0, 0};
    perp = seed.minus(t.e3, seed.dot(t.e3));
    n = perp.norm();
  }
  t.e1 = {perp.x / n, perp.y / n, perp.z / n};
  t.e2 = t.e3.cross(t.e1);
  return t;
}

/// Common factor lambda such that Sum sqrt(m_i^2 + lambda^2 |p_i|^2) = m_want,
/// with the p_i taken in the rest frame of the system.  Monotone in lambda,
/// so a safeguarded Newton converges in a couple of steps from lambda = 1.
bool solve_rescale(const std::vector<double>& m, const std::vector<double>& p,
                   double m_want, double* lambda_out) {
  double m_sum = 0.0;
  for (double mi : m) m_sum += mi;
  if (m_want <= m_sum) return false;
  double lo = 0.0, hi = 1.0;
  auto f = [&](double lam) {
    double s = 0.0;
    for (std::size_t i = 0; i < m.size(); ++i)
      s += std::sqrt(m[i] * m[i] + lam * lam * p[i] * p[i]);
    return s;
  };
  while (f(hi) < m_want) {
    hi *= 2.0;
    if (hi > 1e6) return false;
  }
  double lam = std::min(1.0, hi);
  for (int it = 0; it < 200; ++it) {
    const double val = f(lam) - m_want;
    if (std::fabs(val) <= 1e-14 * m_want) break;
    if (val > 0.0) hi = lam; else lo = lam;
    double der = 0.0;
    for (std::size_t i = 0; i < m.size(); ++i) {
      const double pp = p[i] * p[i];
      der += lam * pp / std::sqrt(m[i] * m[i] + lam * lam * pp);
    }
    double next = (der > 0.0) ? lam - val / der : 0.5 * (lo + hi);
    if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
    lam = next;
  }
  *lambda_out = lam;
  return true;
}

/// A/Z of a 10-digit nuclear code 10LZZZAAAI.
void decode_ion(int pdg, int* a, int* z) {
  const int abspdg = std::abs(pdg);
  *a = (abspdg / 10) % 1000;
  *z = (abspdg / 10000) % 1000;
}

double quark_charge(int id) {
  const int a = std::abs(id);
  const double q = (a == 2 || a == 4 || a == 6) ? 2.0 / 3.0 : -1.0 / 3.0;
  return (id > 0) ? q : -q;
}

}  // namespace

// ------------------------------------------------------ public free funcs

bool dis_scattered_electron(double e_e, const Vec4& p_n, double x, double q2,
                            double phi, Vec4* k_out, double* y_out) {
  if (!(e_e > 0.0) || !(x > 0.0) || !(q2 > 0.0)) return false;
  const Vec4 k{e_e, 0.0, 0.0, -e_e};
  const double pk = dot4(p_n, k);
  if (!(pk > 0.0)) return false;
  const double y = q2 / (2.0 * x * pk);
  if (y_out) *y_out = y;
  if (!(y > 0.0) || y > 1.0) return false;
  // E' + k'_z = Q^2/(2 E_e);  E_N E' - p_Nz k'_z = (1-y) p_N.k
  const double a = q2 / (2.0 * e_e);
  const double den = p_n.e + p_n.pz;
  if (!(den > 0.0)) return false;
  const double ep = ((1.0 - y) * pk + p_n.pz * a) / den;
  const double kz = a - ep;
  const double kt2 = ep * ep - kz * kz;
  if (!(ep > 0.0) || kt2 < 0.0) return false;
  const double kt = std::sqrt(kt2);
  *k_out = {ep, kt * std::cos(phi), kt * std::sin(phi), kz};
  return true;
}

bool dis_parton_fraction(const Vec4& q, const Vec4& p_n, double* xi_out,
                         double xi_max, double m_out) {
  const double q2 = -q.m2() + m_out * m_out;
  const double b = dot4(q, p_n);
  const double disc = b * b + p_n.m2() * q2;
  if (disc < 0.0) return false;
  const double den = b + std::sqrt(disc);
  if (!(den > 0.0)) return false;
  const double xi = q2 / den;
  if (!(xi > 0.0) || xi > xi_max) return false;
  *xi_out = xi;
  return true;
}

HfsSummary hfs_summary(const Event& ev) {
  HfsSummary s;
  for (const auto& p : ev.particles) {
    if (p.role != Role::Hadron || p.status != Status::Final) continue;
    s.e += p.p.e;
    s.pz += p.p.pz;
    s.px += p.p.px;
    s.py += p.p.py;
    s.sigma_empz += p.p.e - p.p.pz;
    s.charge += p.charge;
    ++s.n_total;
    if (std::fabs(p.charge) > 1e-9) ++s.n_charged; else ++s.n_neutral;
  }
  s.pt = std::sqrt(s.px * s.px + s.py * s.py);
  return s;
}

double hfs_sigma_empz_truth(double e_e, double y, const Vec4& p_n) {
  const double den = p_n.e + p_n.pz;
  return 2.0 * e_e * y + ((den > 0.0) ? p_n.m2() / den : 0.0);
}

double hfs_sigma_empz_exact(const Vec4& k, const Vec4& p_n, const Vec4& k_out) {
  const Vec4 h = k + p_n - k_out;
  return h.e - h.pz;
}

// ------------------------------------------------------------------ impl

struct PythiaBridge::Impl {
  BeamConfig beams;
  PythiaBridgeOptions opt;
  std::vector<std::string> applied;

  std::shared_ptr<RngEngine> engine;

  struct Inst {
    int id_beam = 2212;
    double e_a = 0.0, p_a = 0.0, m_a = 0.0;
    std::shared_ptr<LhaupDis> lha;
    std::unique_ptr<Pythia8::Pythia> py;
    Pythia8::PDFPtr pdf;
  };
  Inst proton, neutron;
  /// The COHERENT tier: a POMERON beam (id 990).  PYTHIA counts 990 as a
  /// hadron (`BeamSetup.cc:869-875`), gives it meson-like beam handling
  /// (`BeamParticle.cc:178`, so the remnant is a single ANTIQUARK -- colour-
  /// and charge-neutral against the struck quark) and has a real Pomeron PDF
  /// for it.  Built only when `coherent_t2 == CoherentT2::Pomeron`.
  Inst pomeron;

  double e_e = 0.0;

  PythiaBridgeStats stats;
  std::vector<Particle> hadrons;
  Vec4 struck;
  int struck_pdg = 2212;
  double xi_pythia = 0.0, xi_phys = 0.0;
  double lambda_last = 1.0, w_pythia = 0.0, w_phys = 0.0;
  int quark_id = 0;

  NucleonInCluster cluster_hook;
  bool warned_cluster = false;   ///< the deprecated branch warns once
  NucleonChooser chooser_hook;
  std::shared_ptr<const UnpolSF> f2;   ///< P1: the ByStructureFunctions draw

  /// P1.  P(proton) = Z F2p / (Z F2p + N F2n) at the event's own (x, Q2).
  /// Falls back to Z/A when F2 is not positive there, which is what a flat
  /// draw would have given anyway.
  double proton_fraction_sf(const Event& ev) const {
    const Ion& ion = beams.ion;
    const double flat =
        static_cast<double>(ion.Z) / static_cast<double>(ion.A);
    if (!f2 || !(ev.kin.x > 0.0) || !(ev.kin.q2 > 0.0)) return flat;
    const double zp = static_cast<double>(ion.Z) * f2->f2p(ev.kin.x, ev.kin.q2);
    const double nn = static_cast<double>(ion.N()) * f2->f2n(ev.kin.x, ev.kin.q2);
    const double tot = zp + nn;
    return (tot > 0.0) ? zp / tot : flat;
  }

  Inst* inst_for(int pdg) {
    if (pdg == 990) return pomeron.py ? &pomeron : nullptr;
    if (pdg == 2112) return neutron.py ? &neutron : nullptr;
    return proton.py ? &proton : nullptr;
  }

  void configure(Pythia8::Pythia& py, bool is_pomeron);
  void build(Inst& inst, int id_beam);
  bool run(Event& ev, Rng& rng);
};

void PythiaBridge::Impl::configure(Pythia8::Pythia& py, bool is_pomeron) {
  applied.clear();
  auto set = [&](const std::string& s) {
    applied.push_back(s);
    if (!py.readString(s))
      throw std::runtime_error("PythiaBridge: bad setting '" + s + "'");
  };
  // The hard process arrives through LHAup; frameType 5 is the only mode in
  // which PYTHIA takes per-event four-momenta (survey section 4(v)).
  set("Beams:frameType = 5");
  // BeamSetup::checkBeams accepts a lepton+hadron pair only when a
  // WeakBosonExchange flag is on, or Check:beams is on, or frameType == 4.
  // We want neither an internal DIS process nor frameType 4, so this is the
  // switch that has to be ON -- it is PYTHIA's default, but being explicit
  // keeps a user `settings` override from silently breaking init.
  set("Check:beams = on");
  // The four DIS shower settings of tools/pythia8/gen_dis_hfs.py
  // (= examples/main341.cc).
  set("SpaceShower:dipoleRecoil = on");
  set("SpaceShower:pTmaxMatch = 2");
  set("PDF:lepton = off");
  set("TimeShower:QEDshowerByL = off");
  // A QED ISR photon off the struck quark takes its recoil from the SCATTERED
  // LEPTON (the only other charge in the dipole), which would move the e' the
  // core generator fixed -- measured on 500 events: 0.4 % of them, with the
  // hadronic W shifted by up to 25 %.  Lepton-side QED radiation is out of
  // scope for v0.1 anyway (DEVELOPMENT_PLAN.md section 2), and PYTHIA cannot do
  // it consistently with dipole recoil, so every QED emission that can touch
  // the lepton leg is switched off.  TimeShower:QEDshowerByQ stays ON: photons
  // radiated by the final-state quark belong to the hadronic system and were
  // never observed to move the lepton.
  // The outgoing lepton is handed over exactly massless and must stay that
  // way: with the default (1) PYTHIA reassigns m_e and shuffles the
  // difference onto the struck quark.
  set("SpaceShower:QEDshowerByQ = off");
  set("LesHouches:setLeptonMass = 0");
  // Randomness comes from the RngEngine; this only fixes PYTHIA's own
  // fallback stream.
  set("Random:setSeed = on");
  set("Random:seed = " + std::to_string(1 + (opt.seed % 900000000ULL)));
  set("Next:numberCount = 0");
  set("Next:numberShowInfo = " + std::string(opt.verbosity > 1 ? "1" : "0"));
  set("Next:numberShowProcess = " + std::string(opt.verbosity > 1 ? "1" : "0"));
  set("Next:numberShowEvent = " + std::string(opt.verbosity > 1 ? "1" : "0"));
  if (opt.verbosity < 1) set("Print:quiet = on");
  // The Pomeron parton densities.  6 (H1 2006 Fit B, LO) is PYTHIA's own
  // default and the only LO *H1* set (12 and 13 are LO GKG18; this comment
  // said "the only LO Q^2-dependent set" until 2026-09-04); it is set
  // explicitly so the
  // record of what was run is complete and a PYTHIA default change cannot
  // silently move the coherent flavour draw.  On the nucleon instances these
  // two lines would be inert, so they are simply not applied there.
  if (is_pomeron) {
    set("PDF:PomSet = " + std::to_string(opt.pom_set));
    char buf[64];
    std::snprintf(buf, sizeof(buf), "PDF:PomRescale = %.10g", opt.pom_rescale);
    set(std::string(buf));
  }
  for (const auto& s : opt.settings) set(s);
}

void PythiaBridge::Impl::build(Inst& inst, int id_beam) {
  inst.id_beam = id_beam;
  // Frame choice only: the surrogate's (W^2, Q^2, xi, x) are independent of
  // the PYTHIA-side beam energies, which only have to be large enough for
  // the surrogate to exist.
  const double e_nom = std::sqrt(beams.ion_momentum_per_nucleon *
                                     beams.ion_momentum_per_nucleon +
                                 PROTON_MASS * PROTON_MASS);
  inst.e_a = opt.headroom * e_nom;
  inst.lha = std::make_shared<LhaupDis>(id_beam, inst.e_a, e_e);
  inst.py = std::make_unique<Pythia8::Pythia>(std::string(LIPOLGEN_PYTHIA8_XMLDOC),
                                              opt.verbosity > 0);
  inst.py->setRndmEnginePtr(engine);
  inst.py->setLHAupPtr(inst.lha);
  configure(*inst.py, id_beam == 990);
  if (!inst.py->init())
    throw std::runtime_error("PythiaBridge: pythia.init() failed for beam id " +
                             std::to_string(id_beam));
  inst.m_a = inst.py->particleData.m0(id_beam);
  if (inst.e_a <= inst.m_a)
    throw std::runtime_error("PythiaBridge: beam energy below the nucleon mass");
  // m_0(990) = 0 exactly, which is what makes the coherent zeta come out at
  // beta = Q^2/(M_X^2 + Q^2) with no approximation (docs/PYTHIA_BRIDGE.md
  // sec. 12): with P_A^2 = 0, P_A^+ qtil^- = 2 qtil.P_A = M_X^2 + Q^2.
  inst.p_a = std::sqrt(inst.e_a * inst.e_a - inst.m_a * inst.m_a);
  inst.pdf = inst.py->getPDFPtr(id_beam);
  if (!inst.pdf || !inst.pdf->isSetup())
    throw std::runtime_error("PythiaBridge: no PDF for beam id " +
                             std::to_string(id_beam));
}

// --------------------------------------------------------------- the work

bool PythiaBridge::Impl::run(Event& ev, Rng& rng) {
  ++stats.n_called;
  hadrons.clear();

  const Particle* pbeam = ev.find(Role::BeamElectron);
  const Particle* pscat = ev.find(Role::ScatteredElectron);
  if (!pbeam || !pscat) { ++stats.n_failed; return false; }
  const Vec4 k = pbeam->p;
  const Vec4 kp = pscat->p;

  // ---- resolve the struck nucleon -----------------------------------
  Vec4 p_n;
  int id_n = 2212;
  bool implicit = false;
  if (const Particle* pip = ev.find(Role::Pomeron)) {
    // COHERENT.  `Pipeline::make_coherent` wrote P_IP = P_ion - P_recoil, so
    // q + P_IP IS the diffractive system X and (q + P_IP)^2 = M_X^2 exactly.
    // Everything below then runs unchanged: the surrogate is built at
    // (M_X^2, Q^2) on the POMERON beam instead of (W^2, Q^2) on a nucleon
    // one, and the frame map carries the result onto this X.  Nothing is
    // appended for the target -- the record already names it.
    //
    // P_IP is SPACELIKE (P_IP^2 = t < 0, a few times 1e-2 GeV^2), which
    // `dis_parton_fraction` handles by construction and the frame map never
    // touches; only (q + P_IP)^2 enters.
    p_n = pip->p;
    id_n = 990;
  } else if (const Particle* pn = ev.find(Role::StruckNucleon)) {
    p_n = pn->p;
    id_n = (pn->pdg == 2112) ? 2112 : 2212;
  } else if (const Particle* pc = ev.find(Role::StruckCluster)) {
    // DEPRECATED.  A T1 record names its struck nucleon and never gets here;
    // see `PythiaBridge::NucleonInCluster`.  Warn ONCE per bridge -- this is
    // a configuration mistake, not a per-event condition -- and count every
    // event so a run cannot quietly produce a non-conserving sample.
    ++stats.n_cluster_fallback;
    if (!warned_cluster) {
      warned_cluster = true;
      std::fprintf(stderr,
                   "LiPolGen PythiaBridge: WARNING -- an event carries a "
                   "Role::StruckCluster but no Role::StruckNucleon, so the "
                   "DEPRECATED cluster branch is used and the whole-record "
                   "four-momentum and charge will NOT be conserved.  Run the "
                   "pipeline at PipelineConfig::tier = Tier::T1 (the "
                   "default), which resolves the cluster into a struck "
                   "nucleon plus partner spectators (breakup.hpp).\n");
    }
    int a_c = 1, z_c = 1;
    if (std::abs(pc->pdg) > 1000000000) {
      decode_ion(pc->pdg, &a_c, &z_c);
    } else {
      // A single nucleon written as a "cluster": 2212 / 2112 are NOT
      // 10-digit ion codes and `decode_ion` would read garbage out of them.
      a_c = 1;
      z_c = (pc->pdg == 2212) ? 1 : 0;
    }
    if (a_c < 1) { ++stats.n_failed; return false; }
    if (cluster_hook) {
      p_n = cluster_hook(pc->p, a_c, z_c, rng, &id_n);
    } else {
      // v0: an on-shell nucleon carrying 1/A_c of the cluster three-momentum,
      // with no Fermi smearing at all.
      const double inv = 1.0 / static_cast<double>(a_c);
      const double px = pc->p.px * inv, py = pc->p.py * inv,
                   pz = pc->p.pz * inv;
      id_n = (rng.uniform() * a_c < z_c) ? 2212 : 2112;
      const double m = M_NUCLEON;
      p_n = {std::sqrt(px * px + py * py + pz * pz + m * m), px, py, pz};
    }
    implicit = true;
  } else {
    // Inclusive T0: a nucleon at rest in the ion rest frame, P_ion / A.
    //
    // P3.  ON SHELL AT THE FREE NUCLEON MASS, not at the ion's mass per
    // nucleon (0.9383, not 0.9338).  `InclusiveGenerator::target_nucleon`
    // uses M_NUCLEON, the per-nucleon kinematic labels (W2, nu) are built on
    // M_NUCLEON, and the two implicit targets have to agree or the bridge's
    // own HFS truth identity would be evaluated against a different target
    // from the one the record describes (docs/CONVENTIONS.md).
    const Ion& ion = beams.ion;
    const double pz = beams.ion_momentum_per_nucleon;
    const double m = M_NUCLEON;
    p_n = {std::sqrt(pz * pz + m * m), 0.0, 0.0, pz};
    if (chooser_hook) {
      id_n = chooser_hook(ev, rng);
    } else if (opt.nucleon_choice == NucleonChoice::Proton) {
      id_n = 2212;
    } else if (opt.nucleon_choice == NucleonChoice::Neutron) {
      id_n = 2112;
    } else if (opt.nucleon_choice == NucleonChoice::ByStructureFunctions) {
      // P1: the inclusive rate is Z F2p + N F2n, so draw in that proportion
      // at the event's own (x, Q2) -- the same rule
      // `InclusiveGenerator::proton_fraction` applies, so an event that
      // reaches the bridge with a named struck nucleon and one that does not
      // are drawn from the same mixture.
      id_n = (rng.uniform() < proton_fraction_sf(ev)) ? 2212 : 2112;
    } else {
      id_n = (rng.uniform() * ion.A < ion.Z) ? 2212 : 2112;
    }
    implicit = true;
  }
  Inst* inst = inst_for(id_n);
  if (!inst) {
    // The only way here is a coherent event with `coherent_t2 == Off` (or a
    // neutron with `with_neutron_instance = false`): no instance, no
    // hadronization, the record stays at T0 and still conserves.
    ++stats.n_failed;
    return false;
  }

  // ---- physical invariants -------------------------------------------
  const Vec4 q = k - kp;
  const double q2 = -q.m2();
  const Vec4 h_want = q + p_n;
  const double w2 = h_want.m2();
  if (!(q2 > 0.0) || !(w2 > 0.0)) { ++stats.n_failed; return false; }
  const double w = std::sqrt(w2);
  // The POMERON instance needs a compensated surrogate target.  PYTHIA
  // reconstructs beam B (the electron) at its PHYSICAL mass even though the
  // surrogate hands it over exactly massless, and on a MESON-LIKE beam A
  // (the single-antiquark Pomeron remnant, with no longitudinal freedom) the
  // record closes on the lepton side's MASSIVE light-cone minus: every event
  // comes back short by exactly
  //     Delta(E - pz) = -m_e^2 / (2 E_e)
  // on the hadronic system, i.e. an invariant-mass-squared deficit
  //     Delta(m^2) = -(m_e^2 / (2 E_e)) * (E_h + pz_h)
  //                = -(m_e^2 / (2 E_e)) * (2 E_A - Q^2 / (2 E_e)).
  // Measured on 300 coherent events at config 1: -3.896e-6 GeV^2, CONSTANT
  // to 0.15 % across all M_X and matching this formula to 4 significant
  // digits (the residual spread is its Q^2 term).  A BARYON beam absorbs the
  // same discrepancy in its composite remnant (nucleon-branch residual
  // ~1e-12, no compensation needed).  So on the id 990 branch the surrogate
  // is built at w2 + |Delta(m^2)| while the rescale below still targets the
  // PHYSICAL w: PYTHIA's delivered system then lands on w^2 to ~1e-12
  // instead of carrying a fixed absolute deficit that near-threshold
  // (M_X ~ 1.2 GeV, few heavy hadrons) kinematics amplify into
  // |lambda - 1| ~ 6e-6.
  double w2_sur = w2;
  if (id_n == 990) {
    const double me = inst->py->particleData.m0(11);
    w2_sur += (me * me / (2.0 * e_e)) * (2.0 * inst->e_a - q2 / (2.0 * e_e));
  }
  if (std::sqrt(w2_sur) < inst->m_a + 0.15) { ++stats.n_failed; return false; }
  // A first, massless-quark solve, only to reject kinematics with no root at
  // all; it is redone with the chosen quark's mass once the flavour is known.
  double xi_p = 0.0;
  if (!dis_parton_fraction(q, p_n, &xi_p)) { ++stats.n_failed; return false; }

  // ---- the (W^2, Q^2)-matched surrogate on the PYTHIA beams ----------
  const double e_a = inst->e_a, p_a = inst->p_a, m_a = inst->m_a;
  const Vec4 kb{e_e, 0.0, 0.0, -e_e};
  const double pdotk = e_e * (e_a + p_a);
  const double btil = 0.5 * (w2_sur - m_a * m_a + q2);
  const double acc = q2 / (2.0 * e_e);
  const double etil = (pdotk - btil + p_a * acc) / (e_a + p_a);
  const double kztil = acc - etil;
  const double kt2 = etil * etil - kztil * kztil;
  if (!(etil > 0.0) || kt2 < 0.0) {
    ++stats.n_no_surrogate;
    ++stats.n_failed;
    return false;
  }
  const Vec4 ktil{etil, std::sqrt(kt2), 0.0, kztil};
  const Vec4 qtil = kb - ktil;
  const double lc_plus = e_a + p_a;
  const double qtil_minus = qtil.e - qtil.pz;
  if (!(qtil_minus > 0.0)) { ++stats.n_no_surrogate; ++stats.n_failed; return false; }
  // The struck parton is a massless collinear parton of beam A,
  //     p_in = (zeta P_A^+ / 2) (1, 0, 0, 1),
  // and zeta is fixed by putting the outgoing quark on ITS OWN mass shell:
  //     (qtil + p_in)^2 = m_q^2   ->   zeta = (Q^2 + m_q^2) / (P_A^+ qtil^-).
  // m_q is zero for u/d/s and the PYTHIA table mass for c/b, which is exactly
  // the policy `LesHouches:setQuarkMass = 1` enforces: hand a massless charm
  // over and PYTHIA reassigns m_c and takes the recoil from the *scattered
  // electron*, which would silently move our fixed e'.
  auto zeta_of = [&](double mq) { return (q2 + mq * mq) / (lc_plus * qtil_minus); };

  // ---- flavour: probability ~ e_q^2 x f_q(zeta_q, Q^2) ----------------
  //
  // P2.  Every flavour used to be weighted at the MASSLESS zeta while charm
  // and bottom were then produced at their own, larger zeta_q -- so the pool
  // was priced at one momentum fraction and drawn at another, and a heavy
  // flavour whose zeta_q had run past 1 stayed in the pool and burned a
  // retry every time it was picked.  Each flavour is now weighted at its OWN
  // zeta_q and one with no phase space left is never offered.
  const double zeta_light = zeta_of(0.0);
  if (!(zeta_light > 0.0) || zeta_light >= 1.0) {
    ++stats.n_no_surrogate;
    ++stats.n_failed;
    return false;
  }
  const double q2pdf = std::max(q2, opt.q2_pdf_min);
  int ids[10];
  double cum[10];
  double wgt[10];
  double eq2s[10];
  double zetas[10];
  int nfl = 0;
  auto offer = [&](int id) {
    const int a = std::abs(id);
    const double mq = (a == 4 || a == 5) ? inst->py->particleData.m0(a) : 0.0;
    const double z = zeta_of(mq);
    if (!(z > 0.0) || z >= 1.0) {        // no phase space for this flavour
      ++stats.n_flavour_dropped;
      return;
    }
    const double eq = quark_charge(id);
    const double xf = inst->pdf->xf(id, z, q2pdf);
    ids[nfl] = id;
    eq2s[nfl] = eq * eq;
    wgt[nfl] = eq * eq * std::max(0.0, xf);
    zetas[nfl] = z;
    ++nfl;
  };
  offer(1); offer(2); offer(-1); offer(-2);
  if (opt.include_strange) { offer(3); offer(-3); }
  if (opt.include_charm) { offer(4); offer(-4); }
  if (opt.include_bottom) { offer(5); offer(-5); }
  double tot = 0.0;
  for (int i = 0; i < nfl; ++i) tot += wgt[i];
  if (!(tot > 0.0)) {
    // The second trap of the coherent prototype: at small beta the LO
    // Pomeron grid (PDF:PomSet = 6) has LITERALLY NO QUARKS (gluon momentum
    // fraction 1.000) until Q^2 ~ 1.5-1.75 GeV^2, and the default
    // `q2_pdf_min` = 1.0 clamp sits BELOW that edge -- so this is not a
    // corner case: a default coherent run takes this branch on ~20 % of its
    // events (measured 86/400 at config 1, Q^2 up to 1.59 GeV^2).  Every
    // e_q^2 x f_q weight is zero there and the event would be vetoed for a
    // bookkeeping reason rather than a physical one, so fall back to the
    // bare charge weights e_q^2 -- the flavour-democratic limit of the same
    // formula -- over the LIGHT flavours ONLY.  ALL THREE H1 2006 grids carry
    // no charm or bottom at ANY (beta, Q^2) (xf == 0 identically) -- Fit A
    // NLO (set 3) and Fit B NLO (4) exactly as much as Fit B LO (6), because
    // all three are PYTHIA's one `PomH1FitAB` class and its `xfUpdate`
    // assigns xc = xcbar = xb = xbbar = 0. unconditionally
    // (PartonDistributions.cc:2630; this comment said "the H1 LO grids"
    // until 2026-09-04, which was true but too narrow) -- so
    // e_q^2 x f_q gives the heavy flavours zero weight everywhere and the
    // democratic limit must not resurrect them: before this restriction ~7 %
    // of a default coherent sample came out charm-initiated -- open-charm
    // diffractive final states with zero PDF support, every one on this
    // branch (the prototype's pool, pom_dis.cc, was light-only).  PYTHIA's
    // own backward evolution then still finds the gluon.
    //
    // COUNTED, AND SINCE 2026-09-04 ACTUALLY VISIBLE (D4.7): the counter had
    // been surfaced nowhere -- not in the npz `meta`, not in the run banner.
    // It is now in both.  And MEASURED (D4.6, 4 000 coherent events per
    // point at 6Li config 1, seed 4242), it costs exactly zero: raising
    // `q2_pdf_min` to 1.75 drives the fallback to 0.00 % (sets 6 and 3) or a
    // few per cent (12, 13, 15) -- EXCEPT on `pom_set` 4, H1 2006 Fit B NLO,
    // where it only falls 35.52 % -> 29.32 % and a residual quarter of the
    // sample stays on this branch -- and leaves the final state
    // BIT-IDENTICAL in every one of those cases, set 4 INCLUDED, which is
    // why that residual costs nothing either.  The reason is structural:
    // every Pomeron DPDF PYTHIA ships carries a single light-quark singlet,
    // so e_q^2 x f_q is proportional to e_q^2 exactly over the light
    // flavours and this "fallback" IS the true draw.  See
    // `PythiaBridgeStats::n_pom_flavour_fallback`.
    if (id_n == 990 && nfl > 0) {
      for (int i = 0; i < nfl; ++i) {
        wgt[i] = (std::abs(ids[i]) <= 3) ? eq2s[i] : 0.0;
        tot += wgt[i];
      }
      if (!(tot > 0.0)) { ++stats.n_failed; return false; }
      ++stats.n_pom_flavour_fallback;
    } else {
      ++stats.n_failed;
      return false;
    }
  }
  {
    double acc_w = 0.0;
    for (int i = 0; i < nfl; ++i) { acc_w += wgt[i]; cum[i] = acc_w; }
  }

  // ---- hand it to PYTHIA, retrying with a new flavour on a veto ------
  const double spin_e = (ev.spin.lam_e != 0) ? static_cast<double>(ev.spin.lam_e)
                                             : pbeam->pol;
  double spin_q = 9.0;
  if (const Particle* pn = ev.find(Role::StruckNucleon)) spin_q = pn->pol;

  bool ok = false;
  int id_q = 0;
  double zeta = zeta_light;
  for (int itry = 0; itry <= opt.max_retries; ++itry) {
    if (itry > 0) ++stats.n_retries;
    const double u = rng.uniform() * tot;
    int pick = nfl - 1;
    for (int i = 0; i < nfl; ++i) {
      if (u <= cum[i]) { pick = i; break; }
    }
    id_q = ids[pick];
    // Only c and b carry a mass here, matching LesHouches:setQuarkMass = 1.
    // `zetas[pick]` is the SAME zeta the flavour was weighted at (P2), and it
    // is in (0, 1) by construction, so there is nothing left to reject.
    const int idq_abs = std::abs(id_q);
    const double mq = (idq_abs == 4 || idq_abs == 5)
                          ? inst->py->particleData.m0(idq_abs) : 0.0;
    zeta = zetas[pick];
    const double half = 0.5 * zeta * lc_plus;
    const Vec4 p_in{half, 0.0, 0.0, half};   // exactly massless, along +z
    const Vec4 p_out = qtil + p_in;          // mass m_q by the zeta solve
    DisPayload d;
    d.p_in = Pythia8::Vec4(p_in.px, p_in.py, p_in.pz, p_in.e);
    d.k_in = Pythia8::Vec4(kb.px, kb.py, kb.pz, kb.e);
    d.k_out = Pythia8::Vec4(ktil.px, ktil.py, ktil.pz, ktil.e);
    d.p_out = Pythia8::Vec4(p_out.px, p_out.py, p_out.pz, p_out.e);
    d.m_out = mq;
    d.id_quark = id_q;
    d.x1 = zeta;
    d.q2 = q2;
    d.xf_pdf = std::max(1e-12, inst->pdf->xf(id_q, zeta, q2pdf));
    d.spin_e = spin_e;
    d.spin_q = spin_q;
    d.alpha_em = ALPHA_EM;
    d.alpha_s = inst->py->settings.parm("SigmaProcess:alphaSvalue");
    inst->lha->arm(d);
    engine->set(&rng);
    ok = inst->py->next();
    engine->set(nullptr);
    if (ok) break;
  }
  if (!ok) { ++stats.n_failed; return false; }
  {
    // The physical xi of the target actually struck, with the chosen quark's
    // own mass; it differs from the surrogate's zeta only through
    // p_N^2 != m_A^2 and the light-cone definition (sub-percent).
    const int idq_abs = std::abs(id_q);
    const double mq = (idq_abs == 4 || idq_abs == 5)
                          ? inst->py->particleData.m0(idq_abs) : 0.0;
    double xi_tmp = 0.0;
    if (dis_parton_fraction(q, p_n, &xi_tmp, 1.0, mq)) xi_p = xi_tmp;
  }

  // ---- pull out the final state --------------------------------------
  Pythia8::Event& pev = inst->py->event;
  int i_lep = -1;
  for (int i = 0; i < pev.size(); ++i) {
    if (pev[i].statusAbs() == 23 && pev[i].id() == 11) { i_lep = i; break; }
  }
  int i_lep_final = (i_lep >= 0) ? pev[i_lep].iBotCopyId() : -1;
  if (i_lep_final < 0 || !pev[i_lep_final].isFinal()) {
    // Fallback: the most energetic final electron.  With PDF:lepton = off and
    // TimeShower:QEDshowerByL = off the scattered lepton is inert, so this
    // never fires in practice; keep it so a settings override cannot make us
    // double-count the electron.
    i_lep_final = -1;
    double best = -1.0;
    for (int i = 0; i < pev.size(); ++i) {
      if (pev[i].isFinal() && pev[i].id() == 11 && pev[i].e() > best) {
        best = pev[i].e();
        i_lep_final = i;
      }
    }
    if (i_lep_final < 0) { ++stats.n_failed; return false; }
  }

  std::vector<int> idx;
  Vec4 h_pyth;
  idx.reserve(static_cast<std::size_t>(pev.size()));
  for (int i = 0; i < pev.size(); ++i) {
    if (!pev[i].isFinal() || i == i_lep_final) continue;
    idx.push_back(i);
    h_pyth = h_pyth + from_p8(pev[i].p());
  }
  if (idx.empty()) { ++stats.n_failed; return false; }
  const double m2_pyth = h_pyth.m2();
  if (!(m2_pyth > 0.0) || !(h_pyth.e > 0.0)) { ++stats.n_failed; return false; }

  // ---- map the hadronic system onto the physical one -----------------
  // The surrogate and the physical event share (W^2, Q^2), so this is a pure
  // Lorentz transformation up to PYTHIA-side bookkeeping residuals; a common
  // momentum rescale in the rest frame takes those out, leaving every mass
  // and every charge untouched.  On the NUCLEON instances the residual is
  // the numerical floor (measured max |lambda - 1| ~ 1e-13); on the POMERON
  // instance it would be the constant electron-beam-mass deficit documented
  // at the w2_sur compensation above -- amplified near the M_X floor, where
  // d(Sum E)/d lambda = Sum p_i^2/E_i is small -- which is exactly why that
  // compensation exists.  With it, the Pomeron branch sits at ~1e-11 too.
  const Vec4 q_pyth = kb - from_p8(pev[i_lep_final].p());
  const Triad tp = make_triad(h_pyth, q_pyth, kb);
  const Triad tw = make_triad(h_want, q, k);

  std::vector<Vec4> rest;
  std::vector<double> mrest, prest;
  rest.reserve(idx.size());
  mrest.reserve(idx.size());
  prest.reserve(idx.size());
  for (int i : idx) {
    const Vec4 r = to_rest(from_p8(pev[i].p()), h_pyth);
    rest.push_back(r);
    mrest.push_back(pev[i].m());
    prest.push_back(std::sqrt(r.px * r.px + r.py * r.py + r.pz * r.pz));
  }
  double lambda = 1.0;
  if (!solve_rescale(mrest, prest, w, &lambda)) { ++stats.n_failed; return false; }
  stats.max_rescale_dev = std::max(stats.max_rescale_dev,
                                   std::fabs(lambda - 1.0));
  lambda_last = lambda;
  w_pythia = std::sqrt(m2_pyth);
  w_phys = w;

  hadrons.reserve(idx.size());
  for (std::size_t j = 0; j < idx.size(); ++j) {
    const Vec3 v = vec3(rest[j]);
    const double c1 = lambda * v.dot(tp.e1);
    const double c2 = lambda * v.dot(tp.e2);
    const double c3 = lambda * v.dot(tp.e3);
    Vec4 out;
    out.px = c1 * tw.e1.x + c2 * tw.e2.x + c3 * tw.e3.x;
    out.py = c1 * tw.e1.y + c2 * tw.e2.y + c3 * tw.e3.y;
    out.pz = c1 * tw.e1.z + c2 * tw.e2.z + c3 * tw.e3.z;
    const double mj = mrest[j];
    out.e = std::sqrt(mj * mj + out.px * out.px + out.py * out.py +
                      out.pz * out.pz);
    out = from_rest(out, h_want);

    const Pythia8::Particle& pp = pev[idx[j]];
    Particle h;
    h.pdg = pp.id();
    h.status = Status::Final;
    h.role = Role::Hadron;
    h.p = out;
    h.mass = mj;
    h.charge = pp.charge();
    h.mother1 = -1;
    h.mother2 = -1;
    h.pol = 9.0;
    hadrons.push_back(h);
  }

  // ---- write back -----------------------------------------------------
  for (auto& p : ev.particles) {
    if (p.role == Role::HadronicX) p.status = Status::Intermediate;
  }
  if (implicit) {
    Particle pn;
    pn.pdg = id_n;
    pn.status = Status::Intermediate;
    pn.role = Role::StruckNucleon;
    pn.p = p_n;
    pn.mass = std::sqrt(std::max(0.0, p_n.m2()));
    pn.charge = (id_n == 2212) ? 1.0 : 0.0;
    pn.pol = spin_q;
    ev.particles.push_back(pn);
  }
  ev.particles.insert(ev.particles.end(), hadrons.begin(), hadrons.end());

  struck = p_n;
  struck_pdg = id_n;
  xi_pythia = zeta;
  xi_phys = xi_p;
  quark_id = id_q;
  if (id_n == 990) ++stats.n_pomeron;
  else if (id_n == 2212) ++stats.n_proton;
  else ++stats.n_neutron;
  stats.sum_w2 += w2;
  ++stats.n_ok;
  return true;
}

// ------------------------------------------------------------- interface

PythiaBridge::PythiaBridge(const BeamConfig& beams, PythiaBridgeOptions opt)
    : impl_(new Impl) {
  impl_->beams = beams;
  impl_->opt = std::move(opt);
  impl_->e_e = beams.electron_energy;
  if (!(impl_->e_e > 0.0))
    throw std::runtime_error("PythiaBridge: electron_energy must be positive");
  if (!(impl_->opt.headroom >= 1.0))
    throw std::runtime_error("PythiaBridge: headroom must be >= 1");
  // PDF:PomSet = 11 is REFUSED when the Pomeron instance is built, for the
  // reason `PipelineConfig::validate` refuses any knob that would be recorded
  // without having run.  Set 11 is PYTHIA's `PomHISASD` (a rescaled proton,
  // Angantyr), and it returns its densities only after `setXPom(xPomNow)`,
  // which this bridge never calls -- there is no Angantyr collision system
  // here to supply x_P to the PDF object.  MEASURED, 20 000 coherent events
  // at 6Li config 1: n_pom_flavour_fallback = 20000 / 20000, i.e. 100.00 %
  // -- the Pomeron PDF is not consulted on a single event and the whole run
  // is the charge-democratic e_q^2 fallback, while `meta["pom_set"]` would
  // say 11.  Use any of 1-10 or 12-15 (12-15 are the GKG18 grids), or turn
  // the tier off with `coherent_t2 = Off`.
  if (impl_->opt.coherent_t2 == CoherentT2::Pomeron &&
      impl_->opt.pom_set == 11) {
    throw std::runtime_error(
        "PythiaBridge: PDF:PomSet = 11 (PomHISASD, the Angantyr rescaled "
        "proton) needs PYTHIA's setXPom(x_Pom), which this bridge never "
        "calls -- measured, 100.00 % of coherent events fall back to the "
        "charge-democratic e_q^2 weights and the Pomeron PDF is never "
        "consulted, so recording pom_set = 11 would record a knob that did "
        "not run.  Pick a fitted set (1-10, 12-15) or set coherent_t2 = Off.");
  }
  impl_->f2 = impl_->opt.f2_source
                  ? impl_->opt.f2_source
                  : std::static_pointer_cast<const UnpolSF>(
                        std::make_shared<const ToyF2>());
  impl_->engine = std::make_shared<RngEngine>(impl_->opt.seed);
  impl_->build(impl_->proton, 2212);
  if (impl_->opt.with_neutron_instance) impl_->build(impl_->neutron, 2112);
  // `applied_settings()` returns the LAST-BUILT instance's list.  The
  // Pomeron instance is built last so that, when it exists (the default),
  // the list is the fullest one -- the nucleon settings plus the two
  // `PDF:Pom*` lines.  With `coherent_t2 == Off` the list is the last
  // nucleon instance's and carries no `PDF:Pom*` line, correctly: no
  // instance those settings were applied to exists.
  if (impl_->opt.coherent_t2 == CoherentT2::Pomeron)
    impl_->build(impl_->pomeron, 990);
}

PythiaBridge::~PythiaBridge() = default;

bool PythiaBridge::hadronize(Event& ev, Rng& rng) { return impl_->run(ev, rng); }

const std::vector<Particle>& PythiaBridge::last_hadrons() const {
  return impl_->hadrons;
}
const Vec4& PythiaBridge::last_struck_nucleon() const { return impl_->struck; }
int PythiaBridge::last_struck_nucleon_pdg() const { return impl_->struck_pdg; }
double PythiaBridge::last_xi_pythia() const { return impl_->xi_pythia; }
double PythiaBridge::last_xi_physical() const { return impl_->xi_phys; }
int PythiaBridge::last_quark_id() const { return impl_->quark_id; }
double PythiaBridge::last_rescale() const { return impl_->lambda_last; }
double PythiaBridge::last_w_pythia() const { return impl_->w_pythia; }
double PythiaBridge::last_w_physical() const { return impl_->w_phys; }
const PythiaBridgeStats& PythiaBridge::stats() const { return impl_->stats; }
const PythiaBridgeOptions& PythiaBridge::options() const { return impl_->opt; }
const std::vector<std::string>& PythiaBridge::applied_settings() const {
  return impl_->applied;
}
void PythiaBridge::set_nucleon_in_cluster(NucleonInCluster hook) {
  impl_->cluster_hook = std::move(hook);
}
void PythiaBridge::set_nucleon_chooser(NucleonChooser hook) {
  impl_->chooser_hook = std::move(hook);
}
double PythiaBridge::pythia_sigma_gen_mb() const {
  return impl_->proton.py ? impl_->proton.py->info.sigmaGen() : 0.0;
}

}  // namespace lipolgen
