// Hadronizes N synthetic e + 6Li DIS events through lipolgen::PythiaBridge
// and writes them to a HepMC3 Asciiv3 file, printing the throughput and the
// cross-section bookkeeping.
//
// Usage: hadronize_example [n_events] [output.hepmc]
//
// The (x, Q^2, phi) points are drawn from the generator window by hand -- the
// core sampler (P3) is a separate work package -- with the same head-on
// convention the bridge documents: ion +z at p_u GeV/u, electron -z at E_e,
// phi = azimuth of e' about +z.
//
// Verify externally with pyhepmc:
//   python3 -c "
//   import pyhepmc
//   with pyhepmc.open('hadronize_example.hepmc') as f:
//       evs = list(f)
//   print(len(evs), sum(len(e.particles) for e in evs))
//   "
#include "lipolgen/beams.hpp"
#include "lipolgen/event.hpp"
#include "lipolgen/pythia_bridge.hpp"
#include "lipolgen/rng.hpp"
#ifdef LIPOLGEN_HAVE_HEPMC3
#include "lipolgen/hepmc_writer.hpp"
#endif

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

using namespace lipolgen;

namespace {

constexpr double kEe = 10.0;
constexpr double kPu = 99.5;      // 10 x 99.5 GeV/u, the mid 6Li configuration
constexpr std::uint64_t kSeed = 20260829;

// Generator window of docs/DEVELOPMENT_PLAN.md section 5.
constexpr double kQ2Min = 0.7, kQ2Max = 100.0;
constexpr double kXMin = 1e-3, kXMax = 0.6;
constexpr double kYMin = 0.004, kYMax = 0.985, kW2Min = 8.0;

/// Log-uniform (x, Q^2) and uniform phi, rejected against the window.  This
/// is a stand-in for the real inverse-CDF sampler, not a rate model.
bool sample_point(Rng& rng, const Vec4& p_n, double* x, double* q2,
                  double* phi, Vec4* k_out, double* y) {
  for (int itry = 0; itry < 200; ++itry) {
    *x = kXMin * std::pow(kXMax / kXMin, rng.uniform());
    *q2 = kQ2Min * std::pow(kQ2Max / kQ2Min, rng.uniform());
    *phi = 2.0 * kPi * rng.uniform();
    if (!dis_scattered_electron(kEe, p_n, *x, *q2, *phi, k_out, y)) continue;
    if (*y < kYMin || *y > kYMax) continue;
    const Vec4 k{kEe, 0.0, 0.0, -kEe};
    if ((p_n + k - *k_out).m2() < kW2Min) continue;
    return true;
  }
  return false;
}

Event make_event(const BeamConfig& bc, std::uint64_t number, const Vec4& p_n,
                 int id_n, const Vec4& kp, double x, double q2, double y,
                 double phi) {
  Event ev;
  ev.number = number;
  ev.channel = Channel::Inclusive;
  ev.weight = 1.0;
  ev.spin.lam_e = +1;
  ev.spin.j = 1.0;
  ev.spin.category = "inclusive";
  ev.kin.x = x;
  ev.kin.q2 = q2;
  ev.kin.y = y;
  ev.kin.phi = phi;
  ev.kin.w2 = (p_n + Vec4{kEe, 0.0, 0.0, -kEe} - kp).m2();

  Particle be;
  be.pdg = 11;
  be.status = Status::Beam;
  be.role = Role::BeamElectron;
  be.p = {kEe, 0.0, 0.0, -kEe};
  be.charge = -1.0;
  be.pol = +1.0;
  ev.particles.push_back(be);

  Particle bi;
  bi.pdg = 1000030060;
  bi.status = Status::Beam;
  bi.role = Role::BeamIon;
  bi.p = {std::sqrt(36.0 * kPu * kPu + bc.ion.mass() * bc.ion.mass()), 0.0, 0.0,
          6.0 * kPu};
  bi.mass = bc.ion.mass();
  bi.charge = 3.0;
  ev.particles.push_back(bi);

  Particle sn;
  sn.pdg = id_n;
  sn.status = Status::Intermediate;
  sn.role = Role::StruckNucleon;
  sn.p = p_n;
  sn.mass = std::sqrt(std::max(0.0, p_n.m2()));
  sn.charge = (id_n == 2212) ? 1.0 : 0.0;
  ev.particles.push_back(sn);

  Particle se;
  se.pdg = 11;
  se.status = Status::Final;
  se.role = Role::ScatteredElectron;
  se.p = kp;
  se.charge = -1.0;
  ev.particles.push_back(se);

  Particle hx;
  hx.pdg = 0;
  hx.status = Status::Final;
  hx.role = Role::HadronicX;
  hx.p = p_n + be.p - kp;
  ev.particles.push_back(hx);
  return ev;
}

}  // namespace

int main(int argc, char** argv) {
  const long n_want = (argc > 1) ? std::strtol(argv[1], nullptr, 10) : 2000;
  const std::string out =
      (argc > 2) ? argv[2] : std::string("hadronize_example.hepmc");

  BeamConfig bc;
  bc.electron_energy = kEe;
  bc.ion = LI6();
  bc.ion_momentum_per_nucleon = kPu;

  PythiaBridgeOptions opt;
  opt.seed = kSeed;
  opt.verbosity = 0;

  const auto t_init0 = std::chrono::steady_clock::now();
  PythiaBridge bridge(bc, opt);
  const auto t_init1 = std::chrono::steady_clock::now();

  std::cout << "LiPolGen hadronize_example -- " << bc.label() << "\n"
            << "  PYTHIA settings applied (PhaseSpace:* is irrelevant here:\n"
            << "  the hard process arrives through LHAup, frameType 5):\n";
  for (const auto& s : bridge.applied_settings()) std::cout << "    " << s << "\n";
  std::cout << "  two PYTHIA instances (idA = 2212 and 2112) initialised in "
            << std::fixed << std::setprecision(2)
            << std::chrono::duration<double>(t_init1 - t_init0).count()
            << " s\n\n";

#ifdef LIPOLGEN_HAVE_HEPMC3
  HepMC3Writer writer(out);
#else
  std::cout << "  (built without HepMC3: no file will be written)\n";
#endif

  // The nucleon at rest in the ion rest frame, P_ion / A.
  const double m_per_n = bc.ion.mass_per_nucleon();
  const Vec4 p_n{std::sqrt(kPu * kPu + m_per_n * m_per_n), 0.0, 0.0, kPu};

  long n_written = 0, n_kin_fail = 0;
  long n_charged = 0, n_hadrons = 0;
  double sum_sigma_empz = 0.0, worst_conservation = 0.0;
  const auto t0 = std::chrono::steady_clock::now();
  for (long i = 0; i < n_want; ++i) {
    Rng rng(kSeed, 0, 0, static_cast<std::uint64_t>(i));
    double x = 0, q2 = 0, phi = 0, y = 0;
    Vec4 kp;
    if (!sample_point(rng, p_n, &x, &q2, &phi, &kp, &y)) {
      ++n_kin_fail;
      continue;
    }
    // 6Li has Z = N = 3: alternate p and n so the file holds both.
    const int id_n = (i % 2 == 0) ? 2212 : 2112;
    Event ev = make_event(bc, static_cast<std::uint64_t>(i), p_n, id_n, kp, x,
                          q2, y, phi);
    if (!bridge.hadronize(ev, rng)) continue;

    const HfsSummary h = hfs_summary(ev);
    n_charged += h.n_charged;
    n_hadrons += h.n_total;
    sum_sigma_empz += h.sigma_empz;
    Vec4 sum = ev.find(Role::ScatteredElectron)->p;
    for (const auto& p : ev.particles)
      if (p.role == Role::Hadron) sum = sum + p.p;
    const Vec4 want = ev.find(Role::BeamElectron)->p + p_n;
    worst_conservation =
        std::max(worst_conservation, std::fabs(sum.e - want.e) / want.e);

#ifdef LIPOLGEN_HAVE_HEPMC3
    writer.write(ev);
#endif
    ++n_written;
  }
  const auto t1 = std::chrono::steady_clock::now();
#ifdef LIPOLGEN_HAVE_HEPMC3
  writer.close();
#endif

  const double secs = std::chrono::duration<double>(t1 - t0).count();
  const PythiaBridgeStats& st = bridge.stats();
  std::cout << "  events written        : " << n_written << "\n"
            << "  kinematics rejections : " << n_kin_fail << "\n"
            << "  PYTHIA vetoes (final) : " << st.n_failed
            << "   retries: " << st.n_retries
            << "   no surrogate: " << st.n_no_surrogate << "\n"
            << "  struck p / n          : " << st.n_proton << " / "
            << st.n_neutron << "\n"
            << "  <n_hadrons>           : " << std::setprecision(2)
            << (n_written ? double(n_hadrons) / double(n_written) : 0.0)
            << "   <n_charged>: "
            << (n_written ? double(n_charged) / double(n_written) : 0.0) << "\n"
            << "  <Sum (E - p_z)_had>   : " << std::setprecision(4)
            << (n_written ? sum_sigma_empz / double(n_written) : 0.0)
            << " GeV\n"
            << "  worst |dE|/E          : " << std::scientific
            << worst_conservation << std::fixed << "\n"
            << "  max |lambda - 1|      : " << std::scientific
            << st.max_rescale_dev << std::fixed << "\n"
            << "  time                  : " << std::setprecision(2) << secs
            << " s  ->  " << std::setprecision(0)
            << (secs > 0 ? n_written / secs : 0.0) << " events/s\n";

  // Cross-section bookkeeping.  PYTHIA's own sigmaGen is NOT a physics number
  // here: LHAup strategy 3 hands it a unit cross section per event, so it
  // only reports how many events it was fed.  The physical cross section is
  // the core generator's (Event::xsec_pb, filled by the P3 sampler); this
  // example generates its (x, Q^2) log-uniformly with no rate model, so it
  // has none to report and says so rather than printing a wrong number.
  std::cout << "\n  cross-section bookkeeping\n"
            << "    PYTHIA sigmaGen     : " << std::scientific
            << bridge.pythia_sigma_gen_mb()
            << " mb  (LHAup strategy 3 unit weights -- bookkeeping only,\n"
            << "                            NOT the physical cross section)\n"
            << "    Event::xsec_pb      : "
            << "0 -- this example has no rate model; the P3 sampler fills it\n";

#ifdef LIPOLGEN_HAVE_HEPMC3
  std::cout << "\n  wrote " << out << "\n";
#endif
  return 0;
}
