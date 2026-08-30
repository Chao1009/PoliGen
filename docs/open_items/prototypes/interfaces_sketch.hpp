// interfaces_sketch.hpp -- PROPOSED public interfaces for the three items.
// Design sketch only: it compiles nowhere, it is here to pin the signatures.

#include <memory>
#include <vector>
#include "lipolgen/event.hpp"
#include "lipolgen/rng.hpp"

namespace lipolgen {

// ===========================================================================
// (1) COHERENT T2 -- gamma*-Pomeron DIS through the EXISTING bridge design.
// ===========================================================================
//
// The pipeline already writes everything needed: `Role::IntactRecoil` and the
// `Role::HadronicX` pseudo-particle, plus `CoherentEvent::{m_x2, beta, x_pom,
// t}` (mirrored into `Event::kin`).  The Pomeron four-vector is
//
//     P_IP = P_ion - P_recoil                     (already conserved exactly)
//
// and the mapping onto the DIS bridge's own variables is EXACT:
//
//     target  p_N  ->  P_IP      (PYTHIA beam id 990, m = 0)
//     W^2          ->  M_X^2     ( = (q + P_IP)^2, the pipeline's m_x2 )
//     xi = x_Bj    ->  beta = Q^2/(M_X^2 + Q^2)   -- exact because P_IP^2 = 0
//
// so PythiaBridge grows a THIRD instance (idA = 990) next to its proton and
// neutron ones, and `hadronize()` grows one more target branch.  No new frame
// map, no new LHAup, no new conservation argument.

/// New role, written by `Pipeline::make_coherent`, read by the bridge.
/// (`Role::Pomeron`, Status::Intermediate, pdg 990, charge 0.)
///   p = P_ion - P_recoil,   and the record then reads
///   k + P_ion = k' + P_recoil + sum(hadrons)  exactly, as the tagged one does.

struct CoherentBridgeOptions {
  /// PYTHIA `PDF:PomSet`.  Default = PYTHIA's own default (H1 2006 Fit B LO);
  /// measured gluon fraction 0.93 at beta = 0.01, 0.22 at beta = 0.93.
  int pom_set = -1;
  /// `PDF:PomRescale`.  Physically it is the Pomeron flux normalisation; the
  /// bridge only samples flavour RATIOS from it, so it is inert here.
  double pom_rescale = 1.0;
  /// Floor on the PDF scale, mirroring PythiaBridgeOptions::q2_pdf_min.  The
  /// LO Pomeron grids have NO quarks at all at (beta < 0.1, Q^2 = 1), so the
  /// flavour sampler must fall back to e_q^2 weights there (measured).
  double q2_pdf_min = 1.0;
  /// Below this M_X the string cannot be built for every flavour and PYTHIA
  /// vetoes: measured 29 % at 1.0 GeV, 6 % at 1.2, 0 % at >= 1.4.  Raise
  /// COHERENT_MX_MIN_DEFAULT to this, or route below it to the exclusive
  /// vector-meson channel.
  double m_x_min = 1.2;
};

// ===========================================================================
// (2) TRITON SPECTRAL FUNCTION
// ===========================================================================
//
// Replaces the sequential two-body decay in `ClusterBreakup`'s Triton branch.
// The object it models is S_N(k, E): the joint density of the struck nucleon's
// momentum k and the residual two-nucleon system's excitation.  Two channels:
//
//   struck n:  remnant = d (BOUND, E = 0)          -- "2-body"
//              remnant = (p n) continuum, E > 0    -- "3-body"
//   struck p:  remnant = (n n) continuum ALWAYS    -- "3-body" (no bound nn)
//
// The branching is NOT a constant: it is the k-dependent ratio
// n_0(k) / [n_0(k) + n_1(k)] of Ciofi degli Atti-Simula.  Its k-integral is
// 0.653 for A = 3 (computed from the CS n_0 coefficients themselves, see
// nk_triton.py), which is the 3He(e,e'p)d spectroscopic factor.

enum class TritonChannel { NeutronD, NeutronPnCont, ProtonNnCont };

struct TritonDraw {
  TritonChannel channel;
  double k = 0.0;            ///< struck-nucleon momentum in the t rest frame
  double cos_theta_k = 0.0, phi_k = 0.0;
  double e_rel = 0.0;        ///< relative energy of the continuum pair [GeV]
  double q_pair = 0.0;       ///< its relative momentum [GeV]
  double m_remnant = 0.0;    ///< invariant mass of the two-nucleon remnant
};

/// Backend, in the library's usual shape: a TOY (analytic CS n_0/n_1 forms)
/// and a TABLE (Faddeev/AV18 S(k,E) grid) implementation behind one interface.
class TritonSpectralFunction {
 public:
  virtual ~TritonSpectralFunction() = default;

  /// n(k) of the struck species [GeV^-3], normalised to 1 over d^3k.
  virtual double n_of_k(double k, int pdg) const = 0;
  /// P(residual bound | k, struck species).  0 for a struck proton (nn is
  /// unbound); the CS n_0/(n_0+n_1) ratio for a struck neutron.
  virtual double p_two_body(double k, int pdg) const = 0;
  /// Continuum relative-energy density at this k, normalised on [0, e_max].
  virtual double e_rel_density(double k, double e_rel, int pdg) const = 0;

  /// One draw.  Must consume a FIXED number of uniforms whatever the branch,
  /// so the RNG stream does not depend on the channel (the rule
  /// CoherentXpomModel::draw already follows).
  virtual TritonDraw sample(double x, double q2, double m_s, Rng& rng) const = 0;
};

/// The CS-parameterised implementation: n_0 and n_1 from Ciofi degli Atti-
/// Simula, the 2-body/3-body split taken from their ratio rather than assumed.
class CiofiSimulaTriton : public TritonSpectralFunction { /* ... */ };

// `BreakupOptions` grows one slot; a null pointer keeps today's behaviour:
//   std::shared_ptr<const TritonSpectralFunction> triton_sf;

// ===========================================================================
// (3) FSI ON THE SPECTATOR -- a WEIGHT, never a momentum shift.
// ===========================================================================
//
// The spectator four-vector is the measurement; moving it would break the
// "never recoil-correct the light spectator" rule and the exact conservation
// the whole record rests on.  So FSI enters as a per-event weight
//
//     w_FSI(k, theta_k; W, Q^2) = |Psi_FSI|^2 / |psi|^2
//
// with Psi_FSI(alpha_s, kT) = psi - int d^2k'T/(2pi)^2 Gtil(kT - k'T) psi,
// Gtil(q) = (sigma_Xa (1 - i eps)/2) exp(-B_a q^2/2)  (see fsi_alpha.py).
// It is evaluated on a (k, cos theta_k) grid at construction and read by
// bilinear interpolation per event -- the same shape as
// `acceptance_weights`.

struct FsiKinematics {
  double k = 0.0, cos_theta_k = 0.0;  ///< spectator, ION REST FRAME
  double w = 0.0, q2 = 0.0, x = 0.0;  ///< to drive sigma_XN(W) / formation time
  int spectator_z = 2, spectator_a = 4;
};

/// The hook.  `TaggedSampler` gains `set_fsi(std::shared_ptr<const FsiWeight>)`
/// and multiplies `TaggedEvent::weight` (a NEW field, defaulting to 1) by it;
/// `Pipeline` multiplies it into `Event::weight` exactly as the coherent
/// channel already multiplies in its azimuthal weight.  Null = today (PWIA).
class FsiWeight {
 public:
  virtual ~FsiWeight() = default;
  /// FSI / IA ratio.  >= 0.  1 means no distortion.
  virtual double weight(const FsiKinematics& kin) const = 0;
  /// Norm-preserving variant: weight() divided by its k-integral, for runs
  /// that want the SHAPE distortion without changing the total rate.
  virtual double weight_normalised(const FsiKinematics& kin) const = 0;
  /// Reported so a run can print what it used.
  virtual double sigma_eff_mb(double w) const = 0;
};

/// Cosyn-Weiss / Glauber implementation, transplanted to X-alpha.
struct GlauberFsiOptions {
  double sigma_xn_mb = 40.0;   ///< sigma_XN at large W; see the W ramp below
  double eps = -0.5;           ///< Re/Im of the XN amplitude
  double b_xn = 6.0;           ///< XN slope [GeV^-2]
  /// Formation-length suppression: sigma_eff(W) = sigma_xn_mb * f(W).  The
  /// piece this model is WEAKEST on and the one that has to be banded.
  bool formation_ramp = true;
};
class GlauberFsiWeight : public FsiWeight { /* ... */ };

}  // namespace lipolgen
