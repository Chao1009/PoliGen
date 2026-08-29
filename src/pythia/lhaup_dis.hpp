// Internal header of the PYTHIA tier: the in-memory LHAup that carries one
// DIS hard process, and the RndmEngine that routes PYTHIA's randomness
// through lipolgen::Rng.  Not installed, not part of the public interface.
#ifndef LIPOLGEN_SRC_PYTHIA_LHAUP_DIS_HPP
#define LIPOLGEN_SRC_PYTHIA_LHAUP_DIS_HPP

#include "Pythia8/Pythia.h"

#include "lipolgen/rng.hpp"

namespace lipolgen {
namespace pythia_detail {

/// The LHEF-style record of one e + q -> e' + q' scattering.
///
///   1  incoming quark    (status -1, from beam A = the nucleon)
///   2  incoming electron (status -1, from beam B)
///   3  outgoing electron (status  1)
///   4  outgoing quark    (status  1)
///
/// Everything is in the PYTHIA lab frame, i.e. the frame in which beam A is
/// (E_A, 0, 0, +p_A) and beam B is (E_B, 0, 0, -p_B).  `LesHouches:matchInOut`
/// stays on, so PYTHIA recomputes the incoming energies/longitudinal momenta
/// from the outgoing sum; that requires the incoming partons to lie along
/// +-z, which they do by construction.
struct DisPayload {
  Pythia8::Vec4 p_in;    ///< incoming quark, xi * P_A, exactly massless
  Pythia8::Vec4 k_in;    ///< incoming electron, = beam B
  Pythia8::Vec4 k_out;   ///< outgoing electron
  Pythia8::Vec4 p_out;   ///< outgoing quark, massless
  double m_out = 0.0;    ///< outgoing quark mass (0 for u/d/s, m0 for c/b)
  int    id_quark = 2;
  double x1 = 0.1;       ///< light-cone fraction of beam A carried by p_in
  double q2 = 1.0;       ///< hard scale^2 [GeV^2]
  double xf_pdf = 1.0;   ///< x*f(x, Q^2) of the chosen flavour, for setPdf
  double spin_e = 9.0;   ///< LHEF SPINUP of the incoming electron
  double spin_q = 9.0;   ///< LHEF SPINUP of the incoming quark
  double alpha_em = 0.0073;
  double alpha_s = 0.13;
};

class LhaupDis : public Pythia8::LHAup {
 public:
  LhaupDis(int id_beam_a, double e_beam_a, double e_beam_b)
      : LHAup(3), id_a_(id_beam_a), e_a_(e_beam_a), e_b_(e_beam_b) {}

  bool setInit() override {
    setBeamA(id_a_, e_a_, 0, 0);
    setBeamB(11, e_b_, 0, 0);
    // Strategy 3: unweighted events, process choice made here, unit weight.
    setStrategy(3);
    addProcess(kProcessId, 1.0, 0.0, 1.0);
    xSecSumSave = 1.0;
    xErrSumSave = 0.0;
    return true;
  }

  bool setEvent(int = 0) override {
    if (!armed_) return false;
    // The payload stays armed: `Pythia::next()` may call `setEvent` more than
    // once when the parton level has to be retried, and the right answer then
    // is the same hard process again, not a failure.
    const DisPayload& d = pay_;
    setProcess(kProcessId, 1.0, std::sqrt(d.q2), d.alpha_em, d.alpha_s);
    const bool is_quark = (d.id_quark > 0);
    const int col = is_quark ? kColTag : 0;
    const int acol = is_quark ? 0 : kColTag;
    // 1: incoming quark (beam A side; must come first, ProcessContainer.cc
    //    assigns the first status -21 parton to beam A).
    addParticle(d.id_quark, -1, 0, 0, col, acol, d.p_in.px(), d.p_in.py(),
                d.p_in.pz(), d.p_in.e(), 0.0, 0.0, d.spin_q, -1.0);
    // 2: incoming electron.
    addParticle(11, -1, 0, 0, 0, 0, d.k_in.px(), d.k_in.py(), d.k_in.pz(),
                d.k_in.e(), 0.0, 0.0, d.spin_e, -1.0);
    // 3: outgoing electron.
    addParticle(11, 1, 1, 2, 0, 0, d.k_out.px(), d.k_out.py(), d.k_out.pz(),
                d.k_out.e(), 0.0, 0.0, 9.0, -1.0);
    // 4: outgoing quark, colour-connected to the incoming one so that the
    //    string runs from the struck quark to the beam remnant.
    addParticle(d.id_quark, 1, 1, 2, col, acol, d.p_out.px(), d.p_out.py(),
                d.p_out.pz(), d.p_out.e(), d.m_out, 0.0, 9.0, -1.0);
    setIdX(d.id_quark, 11, d.x1, 1.0);
    setPdf(d.id_quark, 11, d.x1, 1.0, std::sqrt(d.q2), d.xf_pdf / d.x1, 1.0,
           true);
    return true;
  }

  /// Load the next event and allow exactly one `setEvent` call.
  void arm(const DisPayload& d) { pay_ = d; armed_ = true; }

  static constexpr int kProcessId = 9911;   ///< our own LHEF process code
  static constexpr int kColTag = 101;

 private:
  int    id_a_;
  double e_a_, e_b_;
  DisPayload pay_;
  bool   armed_ = false;
};

/// PYTHIA's randomness, taken from the per-event counter-based stream so
/// that a given (seed, run, bunch, event) reproduces bit for bit whatever
/// the order or the thread count (docs/CONVENTIONS.md).
class RngEngine : public Pythia8::RndmEngine {
 public:
  explicit RngEngine(std::uint64_t seed) : fallback_(seed, 0, 0, 0) {}
  double flat() override { return (cur_ ? cur_ : &fallback_)->uniform(); }
  void set(Rng* r) { cur_ = r; }

 private:
  Rng  fallback_;
  Rng* cur_ = nullptr;
};

}  // namespace pythia_detail
}  // namespace lipolgen

#endif  // LIPOLGEN_SRC_PYTHIA_LHAUP_DIS_HPP
