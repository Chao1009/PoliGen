#ifndef LIPOLGEN_COHERENT_HPP
#define LIPOLGEN_COHERENT_HPP

/// \file coherent.hpp
/// Coherent (intact-ground-state) e+6Li channel: scenario rates, recoil
/// tagging and the cos 2phi tensor modulation.  C++17 port of
/// `evgen/polligen/coherent.py`.
///
/// Process: e + 6Li -> e' + X + 6Li(g.s.).  The nucleus stays in its 1+ ground
/// state and recoils with |t| ~= pT^2.  The intact 6Li has A/Z = 2, exactly
/// the beam rigidity, so the only far-forward handle is the Roman-Pot
/// near-beam pT tail; under the exponential coherent t-slope the tagging
/// acceptance is analytic, acc = exp(-B pT_cut^2).
///
/// Everything here is a SCENARIO with explicit bands, to be replaced by a real
/// diffractive model.  The tensor modulation carries TWO mechanisms:
/// a t-linear DEFORMATION term scaled from the polarized-deuteron CGC
/// calculation of Mantysaari et al. (PLB 858:139053) and a FLAT
/// gluon-transversity amplitude bounded by lattice + Drell-Yan estimates.
/// No published calculation exists for any polarized A > 2 nucleus.

#include <string>
#include <vector>

#include "lipolgen/constants.hpp"
#include "lipolgen/event.hpp"
#include "lipolgen/rng.hpp"
#include "lipolgen/spectator.hpp"

namespace lipolgen {

/// 6Li mass quoted in `coherent.py` (6.0151228 u atomic minus 3 m_e) [GeV].
/// It is documentation only: every formula below is mass-free, and the recoil
/// FOUR-VECTOR uses `nuclear_mass(3, 6)`, the AME2020 nuclear mass the rest of
/// the library boosts with.
inline constexpr double M_LI6_DOC = 5.6015;
/// hbar c [GeV fm], as `coherent.GEV_PER_FM_INV` -- an alias of the single
/// definition `HBARC_GEV_FM` (constants.hpp, since 2026-09-01).
inline constexpr double GEV_PER_FM_INV = HBARC_GEV_FM;
/// One-sided rate-weighting model systematic on `a2_tagged`.
inline constexpr double RATE_WEIGHT_SYST = 0.73;

/// |t| truncation of the coherent channel [GeV^2].
///
/// 0.2, NOT the 0.5 carried until 2026-08-29 (P4).  Two reasons, and they
/// agree: the deformation mechanism is scaled from Mantysaari et al.'s
/// polarized-deuteron a_2, digitized over |t| <= 0.30 and LINEAR in |t| only
/// as |t| -> 0, so 0.5 is outside the input; and the linear c_2 crosses -1 at
/// |t| = 0.245 for P_zz = -2, i.e. the azimuthal weight goes NEGATIVE inside
/// the old range.  See `CoherentScenario::positivity_margin`.
inline constexpr double COHERENT_T_MAX_DEFAULT = 0.2;

/// Smallest diffractive mass M_X of the coherent channel [GeV].
///
/// **1.2 since 2026-08-30, raised from 1.0 by the T2 coherent tier.**  The
/// number is now a MEASURED one, not only a scenario choice: the coherent T2
/// path (docs/PYTHIA_BRIDGE.md sec. 12) hands the gamma*-Pomeron system to
/// PYTHIA as an ordinary DIS-like string, and PYTHIA's own hadronization
/// vetoes it when there is not enough mass to make two hadrons out of the
/// struck quark and the Pomeron remnant antiquark.  The veto rate measured on
/// the prototype (Q^2 = 5 GeV^2, 2000 events per point, PomSet 6) is
///
///     M_X   0.8    1.0    1.2    1.4    >= 1.5
///     veto  0.75   0.29   0.06   0.00   0.00
///
/// so 1.2 is where the channel starts costing essentially nothing and 1.4 is
/// where it costs exactly nothing.  1.2 keeps the rate and leaves the
/// residual few-% veto visible in `PythiaBridgeStats::n_failed` rather than
/// hiding it.
///
/// Below it NOTHING IS GENERATED: the coherent channel simply carries no rate
/// there.  That window (M_X < 1.2, i.e. the rho 0.775 / omega 0.783 /
/// phi 1.019 region) is the EXCLUSIVE VECTOR-MESON channel -- a genuinely
/// different process with its own t-slope, spin-density matrix and decay,
/// not a low-mass limit of this one.  It is Phase 2
/// (docs/open_items/code_designs.md sec. 1, option (c)).
///
/// It stays the knob `CoherentXpomModel::m_x_min`.
inline constexpr double COHERENT_MX_MIN_DEFAULT = 1.2;

/// Coherent |F(t)|^2 t-slope B [GeV^-2] for a Gaussian density:
/// F(t) = exp(-R_rms^2 |t|/6) -> |F|^2 = exp(-B|t|), B = R_rms^2/3.
double gaussian_slope(double r_rms_fm);

/// a_2 Fourier coefficients of coherent J/psi photoproduction off transversely
/// polarized deuterons, digitized from arXiv:2408.13213 Fig. 4 (x_P = 1.7e-3):
/// (|t| [GeV^2], a_2(m=0), a_2(m=+-1)).  The wave-function-symmetry relation
/// a_2(0) = -2 a_2(+-1) is the prediction the 6Li scaling inherits.
struct MantysaariRow {
  double t_abs, a2_m0, a2_m1;
};
const std::vector<MantysaariRow>& mantysaari_a2_deuteron();

/// Scenario parameters for coherent diffractive e+6Li (SCENARIO).
struct CoherentScenario {
  /// Coherent fraction of the DIS rate at x -> 0; band {0.02, 0.08}.
  ///
  /// A SCENARIO, not a determination.  eSTARlight (2026-09-02,
  /// docs/open_items/run_2026-09-02/estarlight_li6.md sec. 5) bounds it from
  /// BELOW only: exclusive J/psi inside this channel's M_X >= 1.2 GeV window
  /// gives sigma/sigma_incl = 1.0e-3, 40x under f0 = 0.04, and exclusive
  /// J/psi is a small part of coherent diffraction there.  The all-exclusive-
  /// VM figure 3.0e-2 is NOT an upper bound on f0: ~90 % of it is rho0, with
  /// phi, and both sit BELOW the M_X floor.  Determining f0 needs a coherent
  /// diffractive-DIS calculation (a coherent-A analogue of the H1/ZEUS
  /// diffractive PDFs), which exists in neither eSTARlight nor Sartre.
  double f0 = 0.04;
  /// Coherence falloff scale in x.
  double x_coh = 0.01;
  /// |F(t)|^2 slope [GeV^-2], B = <r^2>/3; band {40, 60}.
  double slope_b = 50.0;
  /// FLAT cos 2phi' modulation of the coherent yield at P_zz = 1 (the
  /// gluon-transversity "exotic glue" scenario); band 3e-3 .. 1e-2.
  double amp = 0.01;
  /// Relative slope modulation Delta B_0/B of the m = 0 state -- the
  /// DEFORMATION mechanism; band -(0.04 .. 0.13).
  double eps_b0 = -0.08;

  /// f_coh(x) = f0 / (1 + (x/x_coh)^2).
  double coherent_fraction(double x) const;
  /// |t| ~ B exp(-B|t|), truncated at t_max [GeV^2].
  double sample_t(Rng& rng, double t_max = COHERENT_T_MAX_DEFAULT) const;
  /// dsigma/d|t| normalized on [0, inf): B exp(-B|t|).
  double dsigma_dt(double t_abs) const;
  /// Fraction of coherent recoils above a near-beam pT cut: exp(-B cut^2).
  double tag_acceptance(double pt_cut) const;
  /// <|t|> of the TAGGED sample: pT_cut^2 + 1/B (exponential tail).
  double mean_t_tagged(double pt_cut) const;
  /// Tag acceptance for an ANGULAR envelope: the pots see the recoil's ANGLE,
  /// so the cut on the nucleus pT is n_sigma sigma_theta A p_u.
  double tag_acceptance_angular(double sigma_theta, double p_per_nucleon,
                                int a_beam = 6, double n_sigma = 10.0) const;

  /// Ensemble cos 2phi' coefficient a_2 of the coherent yield at |t| from the
  /// slope-modulation (deformation) mechanism.  Per pure m state
  /// a_2(m) ~ (Delta B_m/2)|t| with a_2(0) = -2 a_2(+-1); the population
  /// average is a_2 = -(P_zz/4) eps_b0 B |t|.  Exact only as |t| -> 0, where
  /// the m-state phi-averaged rates are equal (see RATE_WEIGHT_SYST).
  double a2_deformation(double t_abs, double pzz) const;
  /// The cos 2phi_t COEFFICIENT of 1 + c_2 cos 2(phi_t - phi_S) from the
  /// deformation mechanism: c_2 = 2 a_2 (the anchor's Eq. (9) normalization).
  double cos2phi_coefficient_deformation(double t_abs, double pzz) const;
  /// <a_2> of the RP-tagged sample in the equal-rate, linear-in-|t|
  /// approximation.
  double a2_tagged(double pt_cut, double pzz) const;
  /// a_2 of a PURE m state: a_2(0) = -2 a_2(+-1), normalized so that the
  /// population average over p_m reproduces `a2_deformation`.  This is the
  /// m-state relation the tests pin.
  double a2_m_state(double t_abs, int m) const;
  /// TOTAL cos 2(phi_t - phi_S) coefficient of the coherent yield: the
  /// deformation coefficient 2 a_2 PLUS the flat gluon-transversity term
  /// amp * P_zz (`money_cos2phi_coherent.py` injects 1 + amp P_zz cos 2phi).
  double cos2phi_coefficient(double t_abs, double pzz) const;

  /// P4.  Minimum over phi_t of the azimuthal weight 1 + c_2 cos 2(phi_t -
  /// phi_S), normalized: 1 - |c_2|.  NEGATIVE means the "weight" is not a
  /// density any more and the sample carries events of negative weight.
  ///
  /// c_2 = -(P_zz/2) eps_B0 B |t| + amp P_zz is LINEAR AND UNBOUNDED in |t|,
  /// so this is the same statement as `InclusiveKernel::positivity_margin`
  /// and is checked the same way -- at setup, over the whole |t| range, and
  /// it throws rather than silently clipping.  |c_2| is monotone in |t|
  /// (2|t| + amp > 0 for every |t| >= 0), so the margin at `t_max` is the
  /// worst one.  At the scenario defaults it is 1 - |P_zz| (2 |t| + 0.01):
  /// zero at |t| = 0.495 for P_zz = +1 and at |t| = 0.245 for P_zz = -2,
  /// which is why `COHERENT_T_MAX_DEFAULT` is 0.2 and not 0.5.
  double positivity_margin(double t_max, double pzz) const;
};

/// Lab kinematics of the intact 6Li recoil.  pT = sqrt(|t|), neglecting
/// t_min ~ (x_P M_A)^2/(1-x_P); the longitudinal momentum keeps (1 - x_pom) of
/// the beam value, so R = p/Z over beam ~ 1.
struct CoherentRecoil {
  double pT = 0.0;
  double theta = 0.0;
  double R = 0.0;
  double xL = 0.0;
  double phi_t = 0.0;
};
/// The geometric form, `reco.recoil_fourvector`'s: `x_pom` is the fraction of
/// the WHOLE-NUCLEUS momentum the pomeron takes.  Kept because it is what the
/// reference tables were dumped from; the generator uses `recoil_lab_of`.
CoherentRecoil recoil_lab(double t_abs, double phi_t, double p_per_nucleon,
                          double x_pom = 0.0, int a_beam = 6);
/// The same observables read off an EXACT recoil four-vector.
CoherentRecoil recoil_lab_of(const Vec4& p_recoil, double phi_t,
                             double p_per_nucleon, int a_beam = 6);

// ------------------------------------------------------------- the pomeron

/// HOW x_P -- AND WITH IT THE DIFFRACTIVE MASS M_X -- IS DRAWN PER EVENT.
///
/// There is no diffractive model to port.  `recopseudo.CoherentResponse`
/// (plans/08 D8) draws x_P log-uniform over a fixed decade and uses it for
/// ONE thing, the t_min kinematic cut; it never forms M_X, and the produced
/// system is not generated at all.  LiPolGen does generate it -- X is the
/// pseudo-particle that closes k + P_ion = k' + P_recoil + X -- so x_P has to
/// be a real per-event variable or X comes out SPACELIKE, which is what
/// happened while x_P was pinned at 0: with P_recoil = P_ion the residual is
/// just the virtual photon and M_X^2 = -Q^2 in 100 % of events.
///
/// CONVENTION.  x_P is the PER-NUCLEON pomeron fraction, the one the
/// diffractive literature quotes and the one whose conventional upper edge is
/// 0.1:
///     x_P = (M_X^2 + Q^2 - t) / (W^2 + Q^2 - M_N^2)  ~  (M_X^2 + Q^2)/(W^2 + Q^2)
/// with W the PER-NUCLEON gamma*N invariant mass the event already carries.
/// The nucleus loses x_P/A of its own light-cone momentum, so the recoil
/// rigidity stays inside the near-beam band: at A = 6 the whole x_P range
/// [x_P,min, 0.1] is a nucleus fraction of at most 1.7e-2, which brackets the
/// [1e-3, 1e-2] decade `recopseudo` draws its (nucleus-fraction) x_pom on.
///
/// beta = x / x_P = Q^2/(M_X^2 + Q^2) <= 1 follows identically.
struct CoherentXpomModel {
  /// Smallest diffractive mass [GeV]; sets the LOWER x_P edge per event.
  double m_x_min = COHERENT_MX_MIN_DEFAULT;
  /// Upper x_P edge: the conventional edge of the diffractive region.
  double x_pom_max = 0.1;

  /// x_P of a diffractive mass, and the inverse.
  static double x_pom_of(double m_x2, double q2, double w2);
  static double m_x2_of(double x_pom, double q2, double w2);
  /// The lower edge x_P(M_X,min) at this (Q^2, W^2).
  double x_pom_min(double q2, double w2) const;
  /// Log-uniform on [x_pom_min, x_pom_max].  When the kinematics cannot fit
  /// M_X,min below `x_pom_max` -- large x, where the coherent weight f_coh(x)
  /// is 1e-4 of its peak anyway -- the draw DEGENERATES to `x_pom_min`, so
  /// M_X = M_X,min exactly and the event is still physical.  One uniform is
  /// consumed either way, so the RNG stream does not depend on the branch.
  double draw(double q2, double w2, Rng& rng) const;
};

/// The DIS side an exact coherent recoil needs.
struct CoherentDis {
  double x = 0.0, q2 = 0.0, w2 = 0.0;
  /// R = k + P_ion - k', the four-momentum the recoil and X share.
  Vec4 residual;
};

/// Rigidity ratio R of a beam-velocity fragment: R = (m/Z)/(m_beam/Z_beam), a
/// ratio of MASS-to-charge ratios, not of mass numbers.  NaN for Z = 0.
double fragment_rigidity(int a, int z, int beam_a = 6, int beam_z = 3);

/// Far-forward destination of a beam-velocity fragment at IP6, from the
/// rigidity windows alone.  Returns the same label strings as
/// `coherent.fragment_route_label`.
std::string fragment_route_label(int a, int z, int beam_a = 6, int beam_z = 3,
                                 const std::string& config = "18x275");

/// One fragment of a breakup channel.
struct BreakupFragment {
  std::string name;
  int a = 0, z = 0;
};
/// One breakup channel: name, threshold [MeV], fragments.
struct BreakupChannel {
  std::string name;
  double threshold_mev = 0.0;
  std::vector<BreakupFragment> fragments;
};
/// Lowest-lying particle decompositions of the beam nucleus (6Li or 7Li).
const std::vector<BreakupChannel>& breakup_table(int beam_a = 6, int beam_z = 3);

/// One row of the veto table: fragment name, rigidity, destination.
struct VetoRow {
  std::string channel, fragment;
  double rigidity = 0.0;
  std::string destination;
};
std::vector<VetoRow> veto_table(int beam_a = 6, int beam_z = 3,
                                const std::string& config = "18x275");

/// One coherent event.
struct CoherentEvent {
  double t = 0.0;       ///< |t| the slope was sampled at; the recoil's pT^2
  double t_exact = 0.0; ///< |(P_ion - P_recoil)^2|, i.e. `t` plus |t_min|
  double phi_t = 0.0;   ///< recoil azimuth about the ion axis [rad]
  double x_pom = 0.0;   ///< PER-NUCLEON pomeron fraction (CoherentXpomModel)
  double m_x2 = 0.0;    ///< M_X^2 of the diffractive system [GeV^2], > 0
  double beta = 0.0;    ///< x / x_P = Q^2/(M_X^2 + Q^2)
  CoherentRecoil recoil;
  Vec4 p_recoil;        ///< head-on frame four-vector of the intact nucleus
  double weight = 1.0;  ///< azimuthal weight 1 + c_2 cos 2(phi_t - phi_S)
  double c2 = 0.0;      ///< the coefficient that weight carries
  int route = kRouteLost;
};

/// Sampler for the coherent channel: |t| from the exponential slope, the
/// recoil azimuth from 1 + c_2 cos 2(phi_t - phi_S), the beam-rigidity recoil
/// four-vector and the far-forward route.
///
/// The azimuth is sampled UNWEIGHTED and the modulation is carried as an
/// event weight by default (`weighted_azimuth = false`), which is what a
/// tensor-modulation fit wants; setting `weighted_azimuth` draws phi_t from
/// the modulated density by rejection instead and leaves the weight at 1.
class CoherentSampler {
 public:
  CoherentSampler(const CoherentScenario& scenario, double p_per_nucleon,
                  double pzz, double phi_s = 0.0, int beam_a = 6,
                  int beam_z = 3);

  const CoherentScenario& scenario() const { return sc_; }
  const CoherentXpomModel& xpom_model() const { return xpom_; }
  double t_max() const { return t_max_; }
  void set_optics(const Optics& optics, const std::string& pot_config);
  void set_weighted_azimuth(bool on) { weighted_ = on; }
  void set_xpom_model(const CoherentXpomModel& m) { xpom_ = m; }
  /// THROWS if the azimuthal weight would go negative anywhere on
  /// [0, t_max] at this sampler's P_zz (`CoherentScenario::positivity_margin`).
  void set_t_max(double t_max);

  /// One coherent event AT this event's DIS kinematics.  `dis` is required:
  /// x_P, M_X and the recoil's longitudinal momentum are all functions of it,
  /// and the recoil is solved so that
  ///     (R - P_recoil)^2 = M_X^2   EXACTLY, with P_recoil^2 = M_A^2,
  /// which is what makes X timelike event by event.
  CoherentEvent sample(Rng& rng, const CoherentDis& dis) const;

  /// ADD the intact recoil (Role::IntactRecoil) to an event whose beams and
  /// e' are already set, and fill `Event::kin.t` / `kin.x_pom`.
  void fill_event(Event& ev, const CoherentEvent& ce) const;

 private:
  void check_positivity() const;

  CoherentScenario sc_;
  CoherentXpomModel xpom_;
  double p_u_;
  double pzz_;
  double phi_s_;
  int beam_a_, beam_z_;
  double m_beam_;
  double t_max_ = COHERENT_T_MAX_DEFAULT;
  bool weighted_ = false;
  Optics optics_;
  std::string pot_config_;
};

}  // namespace lipolgen

#endif  // LIPOLGEN_COHERENT_HPP
