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

#include "lipolgen/event.hpp"
#include "lipolgen/rng.hpp"
#include "lipolgen/spectator.hpp"

namespace lipolgen {

/// 6Li mass quoted in `coherent.py` (6.0151228 u atomic minus 3 m_e) [GeV].
/// It is documentation only: every formula below is mass-free, and the recoil
/// FOUR-VECTOR uses `nuclear_mass(3, 6)`, the AME2020 nuclear mass the rest of
/// the library boosts with.
inline constexpr double M_LI6_DOC = 5.6015;
/// hbar c [GeV fm], as `coherent.GEV_PER_FM_INV`.
inline constexpr double GEV_PER_FM_INV = 0.19733;
/// One-sided rate-weighting model systematic on `a2_tagged`.
inline constexpr double RATE_WEIGHT_SYST = 0.73;

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
  double sample_t(Rng& rng, double t_max = 0.5) const;
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
CoherentRecoil recoil_lab(double t_abs, double phi_t, double p_per_nucleon,
                          double x_pom = 0.0, int a_beam = 6);

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
  double t = 0.0;       ///< |t| [GeV^2]
  double phi_t = 0.0;   ///< recoil azimuth about the ion axis [rad]
  double x_pom = 0.0;
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
  void set_optics(const Optics& optics, const std::string& pot_config);
  void set_weighted_azimuth(bool on) { weighted_ = on; }
  void set_x_pom(double x_pom) { x_pom_ = x_pom; }
  void set_t_max(double t_max) { t_max_ = t_max; }

  CoherentEvent sample(Rng& rng) const;

  /// ADD the intact recoil (Role::IntactRecoil) to an event whose beams and
  /// e' are already set, and fill `Event::kin.t` / `kin.x_pom`.
  void fill_event(Event& ev, const CoherentEvent& ce) const;

 private:
  CoherentScenario sc_;
  double p_u_;
  double pzz_;
  double phi_s_;
  int beam_a_, beam_z_;
  double m_beam_;
  double x_pom_ = 0.0;
  double t_max_ = 0.5;
  bool weighted_ = false;
  Optics optics_;
  std::string pot_config_;
};

}  // namespace lipolgen

#endif  // LIPOLGEN_COHERENT_HPP
