// SPDX-License-Identifier: GPL-3.0-or-later
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
/// 0.2, NOT the 0.5 carried until 2026-08-29 (P4).  It STAYS FIXED, and since
/// 2026-09-04 the reason it is written down is the ANCHOR RANGE, not
/// positivity.  The two are not interchangeable and this header used to
/// present them as "two reasons, and they agree" -- true at the shipped
/// eps_b0 and false at every other 6Li value of it (D5).
///
///   PRIMARY, and INDEPENDENT of every knob: the deformation mechanism is
///   scaled from Mantysaari et al.'s polarized-deuteron a_2, and
///   `mantysaari_a2_deuteron()` carries FOUR digitized rows, |t| = 0.05,
///   0.10, 0.20, 0.30.  The fit is linear in |t| and exact only as
///   |t| -> 0, so 0.2 is inside the input and 0.5 is outside it.  This
///   reason is a property of the INPUT TABLE; it does not move when
///   `CoherentScenario::eps_b0` moves.
///
///   SECONDARY, and CONTINGENT on eps_b0: the linear c_2 crosses -1 at
///   |t| = 0.245 for P_zz = -2 (0.495 at P_zz = +1).  That number is
///   `CoherentScenario::t_positivity_edge`, which DERIVES it rather than
///   repeating it, and it moves with eps_b0: at the measured 6Li quadrupole
///   (eps_b0 = -0.0070 ROUNDED, sec. C4; -0.0070024 as derived) the same
///   formula returns |t| = 2.80 GeV^2 (2.7990 on the derived value,
///   2.8000 on the rounded one -- quote 2.80, never 2.8000, see
///   `t_positivity_edge`), i.e. 9.3x OUTSIDE the anchor's own |t| <= 0.30.  Positivity therefore
///   STOPS BINDING the moment eps_b0 is corrected; the anchor range does
///   not.  Deriving t_max from positivity alone would license extrapolating
///   a linear-in-|t| fit ten times past its data.
///
/// What the ceiling costs, measured on 200 000 generated coherent events at
/// the shipped defaults (seed 99, plan `tensor-thirds` at pz = 0.7,
/// pzz = 0.6, 4 threads -- the plan splits the events across categories, so
/// a single-category run of the same seed gives 0.019974; max |t| and the
/// zero count above 0.20 are the same either way,
/// `phase_D_small_items.md` D5.3): <|t|> = 0.019963 GeV^2, max |t| =
/// 0.197575, and `sample_t` RENORMALIZES on [0, t_max], so the truncation
/// does not lose rate -- it redistributes exp(-B t_max) = 4.54e-5 of it.
/// The ceiling is not a rate question: it is a statement about where the
/// model is DEFINED.
///
/// This depends on the eps_b0 author decision (`CoherentScenario::eps_b0`,
/// STATUS.md decision table row 8) only through the secondary reason.  See
/// also `CoherentScenario::positivity_margin`, which stays exactly as it is:
/// it is the GUARD on the truncated weight the sampler actually uses, not
/// the derivation of this constant.
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

/// One eSTARlight coherent vector-meson row for e 10 GeV x 6Li 99.5 GeV/u
/// (LiPolGen beam configuration 1), the EXTERNAL unpolarized rate baseline of
/// open item 11.1.  Run 2026-09-02 at commit
/// 939b11a24499398392d959db81c7502aeec91046, 2e5 events per channel; the run
/// log, the input files and the adversarial re-check are
/// docs/open_items/run_2026-09-02/estarlight_li6.md (sec. 2a, 2b, 2e, 5).
/// This is the SINGLE code home of those cross sections; every other mention
/// in the tree is prose quoting it.
///
/// `sigma_nb` is over the 0.1 < Q^2 < 100 GeV^2 window (which is
/// arXiv:2511.05638's ACCEPTANCE-STUDY range verbatim -- a kinematic range
/// copied from that paper, NOT a physics window, and NOT eSTARlight's full
/// Q^2 reach) at eSTARlight's DEFAULT light-nucleus
/// density, R_G = 1.2 A^{1/3} = 2.1805 fm.  Most of the coherent rate is
/// BELOW it: see `estarlight_li6_q2_floors()`, where removing the 0.1 GeV^2
/// floor multiplies sigma by 6.75 (J/psi), 21.7 (phi), 35.2 (rho).
/// `sigma_q7_nb` is the same run
/// restricted to Q^2 > 0.7 GeV^2, LiPolGen's own generator window (sec. 5);
/// `sigma_rmeas_nb` is the 0.1 < Q^2 < 100 run with the Angeli-Marinova
/// measured rms charge radius 2.589 fm patched in (sec. 2b) -- the SAME
/// physics with a different density, and the leading rate systematic, +-30 %.
/// `b_fit` is the maximum-likelihood dN/d|t| slope of that same default-
/// density sample over 0 < |t| < 0.10 GeV^2, and it must be used WITH
/// `sigma_nb` (a larger radius raises B by 41 % and lowers sigma by 29 %:
/// sigma and B are not independent knobs).  `branching` is the decay actually
/// reconstructed, as quoted in sec. 2e: J/psi -> e+e- (0.0597), phi -> K+K-
/// (0.49), rho -> pi pi (1.0); the J/psi row deliberately does NOT include
/// mu+mu- -- `branching_all` is the field that does, and it is what a J/psi
/// analysis reconstructing both leptons actually has (1.4137x the reach).
///
/// CAVEATS THAT RIDE WITH EVERY ROW (estarlight_li6.md sec. 6): one
/// spherically symmetric Gaussian density, so NO alpha+d clustering, NO
/// polarization axis and NO diffractive minimum at any A; unpolarized
/// throughout, so it constrains neither of this scenario's two tensor
/// mechanisms (the flat gluon-transversity amplitude, the m = 0 slope
/// modulation); no saturation.  A rate and slope baseline, never an
/// imaging one.
struct EstarlightLi6Row {
  const char* vm;          ///< "jpsi" | "phi" | "rho"
  double sigma_nb;         ///< 0.1 < Q^2 < 100 GeV^2, R_G = 2.1805 fm
  double sigma_q7_nb;      ///< Q^2 > 0.7 GeV^2, same density
  double sigma_rmeas_nb;   ///< 0.1 < Q^2 < 100 GeV^2, R = 2.589 fm
  double b_fit;      ///< fitted B [GeV^-2] at the default density
  double branching;        ///< branching fraction of the quoted decay
  /// B [GeV^-2] of the SAME run at the measured radius R = 2.589 fm, i.e. the
  /// slope that belongs with `sigma_rmeas_nb` (sec. 2b).  Added 2026-09-04:
  /// it had lived only in prose and in a dict inside
  /// validation/o5_a2_reach.py, so a re-fit of sec. 2b could not reach it.
  double b_rmeas;
  /// Branching fraction of EVERY reconstructible decay of the same VM, as
  /// eSTARlight itself defines them (starlightconstants.h): J/psi
  /// JpsiBree + JpsiBrmumu = 0.05971 + 0.05961 = 0.11932, i.e. BOTH lepton
  /// channels; phi -> K+K- and rho -> pi pi have no second channel, so this
  /// equals `branching` for those two.  Verified 2026-09-04 by running
  /// PROD_PID = 443013: eSTARlight reports 106.801 pb generated against a
  /// 1.792 nb total, i.e. 0.05960.  A J/psi analysis that reconstructs both
  /// leptons has sqrt(branching_all / branching) = 1.4137x the reach of one
  /// that does not -- sqrt(2) to 0.05 %.
  double branching_all;
};
const std::vector<EstarlightLi6Row>& estarlight_li6_coherent();

/// GLOBAL far-forward detection efficiency for COHERENT J/psi off 7Li at
/// 18 x 118 GeV/u through the IR-8 secondary focus, arXiv:2511.05638 p. 4
/// (Chang et al., PRD 113 (2026) 032018): "The global detection efficiency as
/// the increasing nucleon number, from deuterium to oxygen, is 47.12 %,
/// 32.23 %, 29.42 %, 17.75 %, ..." -- 7Li is the fourth entry.
///
/// WHAT IT IS.  "tagging efficiency x acceptance" (their Figs. 2 and 5 axis)
/// for the INTACT RECOIL NUCLEUS in the Roman Pots at the secondary focus.
/// It is not an electron-arm efficiency and does not require the scattered
/// electron to be seen, and the paper states it "only accounts for the
/// acceptance effect and does not incorporate the efficiencies of the
/// detector.  Additionally, we did not account for the efficiency and
/// acceptance of the reconstructed distribution" (p. 4, sec. IV).
///
/// SO IT IS A RECOIL-NUCLEUS NUMBER AND NOTHING ELSE.  It carries NO
/// central-detector acceptance and NO reconstruction efficiency for the
/// DECAY LEPTONS of the vector meson.  A chain of the form
/// sigma x BR(l+l-) x COHERENT_JPSI_EFF_IR8_LI7 is therefore missing a
/// factor A_ll x eps_ll <= 1, direction DOWN.  What is bounded: the
/// GEOMETRIC part, which `o5_lepton_pair_acceptance()` in
/// validation/o5_a2_reach.py computes from this tree's own
/// `Scenario::eta_max` = 3.5 and the sample's own <W>, is 0.99 at
/// <W> = 30.2 GeV and never below 0.89 anywhere in the accessible W range --
/// the geometry is NOT where the factor is.  What is NOT bounded anywhere in
/// this tree: the per-lepton tracking/PID/reconstruction efficiency.  See
/// phase_C_numbers.md sec. C2.8 item 3b.
///
/// AND IT IS A TOP-ENERGY NUMBER.  18 x 117.9 GeV/u is 7Li's own top energy
/// (Z/A x 275), NOT the 10 x 99.5 GeV/u of `estarlight_li6_coherent()`; this
/// tree's configuration-identical eSTARlight row for those beams has
/// <W> = 43.2 GeV against 30.2 for the 6Li sample the reach is priced on.
/// The efficiency is NOT flat in beam energy: `chang26_he3_energy_scan()`
/// below is the same paper measuring exactly that dependence, and it rises
/// steeply as the energy falls.  Applying 0.1775 at 10 x 99.5 is therefore
/// CONSERVATIVE, by a factor the scan puts at 1.12-1.16 for the energy step
/// alone (o5_a2_reach.py, `beam_energy_scaling`).
///
/// WHAT IT IS NOT.  There is no 6Li number anywhere; open item O5 uses this
/// as a STAND-IN for 6Li and for all three vector mesons, and every table it
/// enters is also printed at efficiency 1.  THE SIZE AND THE SIGN OF THAT
/// SUBSTITUTION ARE NOT ESTABLISHED.  Until 2026-09-04 this comment said it
/// was "worth a further 1.16-1.22 UP" because the species list is at fixed
/// rigidity and 6Li is lighter than 7Li.  That was a non sequitur and is
/// RETRACTED: fixing the rigidity A/Z x E does NOT hold E/u or the total
/// momentum fixed (see `chang26_species_efficiency()`), and read off the four
/// entries that DO share a beam energy the substitution is x0.99 to x1.33 --
/// it straddles 1, and its direction is undetermined
/// (`o5_a2_reach.py`, `species_scaling_same_energy`).  Stated, not applied.
/// Their sample was generated with 0.1 < Q^2 < 100 GeV^2, so no point of it
/// was ever evaluated in the photoproduction region Q^2 < 0.1 that carries
/// most of the rate (`estarlight_li6_q2_floors()`): applying it there is an
/// EXTRAPOLATION, bounded but unmeasured.  See validation/o5_a2_reach.py and
/// docs/open_items/run_2026-09-03/phase_C_numbers.md sec. C2.
///
/// The single code home of that number; before 2026-09-04 it was a literal
/// inside validation/o5_a2_reach.py.
inline constexpr double COHERENT_JPSI_EFF_IR8_LI7 = 0.1775;
/// COHERENT_JPSI_EFF_IR8_LI7 = 0.1775 is a PURE GEOMETRIC ACCEPTANCE and
/// therefore an UPPER bound on the far-forward intact-recoil tagging efficiency
/// x acceptance: Chang et al. (PRD 113 (2026) 032018, sec. IV, last sentence)
/// state that the simulation 'only accounts for the acceptance effect and does
/// not incorporate the efficiencies of the detector' and that 'the efficiency
/// and acceptance of the reconstructed distribution' were not included. Each
/// omitted factor is <= 1, so restoring any can only LOWER the product;
/// direction DOWN, unquantified, and it COMPOUNDS with the x1.12 beam-energy
/// leg rather than cancelling it (surveyed 2026-09-15,
/// run_2026-09-06/phase_C_survey.md C-S6). It bounds the RECOIL leg only, not
/// the decay leptons.

/// The ION beam energy `COHERENT_JPSI_EFF_IR8_LI7` was measured at [GeV/u]:
/// 7Li's own top energy, Z/A x 275 = 117.857, which arXiv:2511.05638 writes
/// as "18 x 118".  `beams.hpp` derives the same 117.9 from the AME2020
/// masses.  Recorded beside the efficiency because the efficiency is NOT
/// flat in it and open item O5 applies it at 99.5.
inline constexpr double COHERENT_JPSI_EFF_IR8_LI7_E_ION_GEV = 117.9;

/// <W> of THIS tree's configuration-identical eSTARlight run at those beams
/// (7Li J/psi, 18 x 117.9, estarlight_li6.md sec. 2a) [GeV].  The 6Li sample
/// open item O5 prices sits at 30.2 (no Q^2 floor) / 32.2 (Q^2 > 0.1), i.e.
/// LOWER -- and by the same paper's Fig. 2 lower W means HIGHER efficiency,
/// the same direction as the beam-energy step.  Neither is a measurement of
/// the transfer; both say it is conservative.
inline constexpr double COHERENT_JPSI_EFF_IR8_LI7_W_MEAN_GEV = 43.2;

/// One (collision energy, global far-forward efficiency) point of
/// arXiv:2511.05638 sec. V.B -- the SAME quantity as
/// `COHERENT_JPSI_EFF_IR8_LI7`, measured on e+3He at three EIC energy
/// configurations instead of one.  Verbatim (p. 5-6): "At the top collision
/// energy, 32.23 % of the scattered 3He nuclei occur within a safe distance
/// from the beam ... At the energy of 10 x 100 GeV^2, the total detection
/// efficiency is 54.38 % ... For the lowest collision energy, 5 GeV electron
/// beams on 41 GeV 3He beams, the detection efficiency is 99.77 %."
///
/// WHY IT IS HERE.  It is the ONLY measurement of the beam-energy dependence
/// of this efficiency that this tree has seen, and open item O5 applies a
/// TOP-energy 7Li number to a 10 x 99.5 GeV/u 6Li sample.  The dependence is
/// steep and one-directional -- d ln(eff)/d ln(E_ion) = -0.866 between the
/// first two points and -0.681 between the last two -- so the transfer is
/// conservative, and now by a stated amount rather than by assertion.
/// `validation/o5_a2_reach.py` reads this table; `phase_C_numbers.md`
/// sec. C2.8 item 3a is the write-up.
///
/// CAVEATS.  (1) 3He, not Li: this bounds the ENERGY lever only, at fixed
/// species.  (2) The electron energy moves with the ion energy at all three
/// points, so the two cannot be separated here.  (3) The efficiency is
/// bounded above by 1, so the slope must flatten -- which is why the two
/// local slopes differ and why o5_a2_reach.py carries both as a band rather
/// than fitting one power law.  (4) The paper labels the lower two points by
/// the nominal PROTON configuration ("10 x 100", "5 x 41") and its text says
/// "5 GeV electron beams on 41 GeV 3He beams", i.e. those ion energies are
/// NOT Z/A-scaled the way the top one (183 = 2/3 x 275) is; `e_ion_gev` below
/// is what the paper's own text says the ion beam was.
struct Chang26EffEnergyRow {
  double e_electron_gev;   ///< electron beam energy [GeV]
  double e_ion_gev;        ///< 3He beam energy [GeV/nucleon], as the paper states it
  double efficiency;       ///< global far-forward tagging efficiency x acceptance
};
const std::vector<Chang26EffEnergyRow>& chang26_he3_energy_scan();

/// One (nucleus, beam energy, global far-forward efficiency) entry of
/// arXiv:2511.05638's p. 4 species list -- the list
/// `COHERENT_JPSI_EFF_IR8_LI7` is the fourth entry of -- with the beam energy
/// each was measured at, from Fig. 2's own legend.
///
/// THE ONE FACT THIS TABLE EXISTS TO MAKE CHECKABLE: every entry is at the
/// SAME MAGNETIC RIGIDITY.  A/Z x e_ion_gev = 274.0 to 275.3 GeV/e across all
/// seven, because each nucleus is at its own top energy Z/A x 275.  So the
/// fall from 47.12 % (2D) to 1.59 % (16O) is NOT a rigidity effect.
///
/// AND THAT IS ALL IT MAKES CHECKABLE.  Until 2026-09-04 this comment, and
/// four other sites, went on: "so its A-ordering is a species lever on its
/// own".  IT IS NOT, and the sentence is RETRACTED.  Holding A/Z x E fixed
/// eliminates rigidity and NOTHING ELSE: at fixed R the per-nucleon energy is
/// E/u = R Z/A and the total beam momentum is p_z = Z R, so BOTH still vary
/// down the list -- E/u from 118 GeV/u (7Li) to 183 (3He), Z from 1 to 8,
/// p_z from 274 to 2192 GeV.  The list is a JOINT (A, Z, E/u) lever, and
/// nothing in this table decomposes it.  The size of the confound is not
/// small: `chang26_he3_energy_scan()` measures d ln(eff)/d ln(E_ion) =
/// -0.68 to -0.87, so the list's own x1.55 spread in E/u is worth up to
/// x1.46 in the efficiency -- the whole size of the "species gain" that was
/// being read off it.
///
/// WHAT CAN BE READ WITHOUT THAT CONFOUND: the FOUR entries that share a
/// beam energy, 2D, 4He, 12C and 16O at 137 GeV/u.  6Li's own fixed-rigidity
/// energy is Z/A x 275 = 137.5, so those four bracket A = 6 with NO energy
/// step at all, and 4He -> 12C is the adjacent pair.  Read that way the
/// 7Li -> 6Li substitution open item O5 has to make is x0.99 to x1.33 -- it
/// STRADDLES 1 (`o5_a2_reach.py`, `species_scaling_same_energy`), against the
/// x1.16-1.22 the confounded reading gave.
///
/// NOT A FIT.  Every one of the seven entries is carried, because it takes
/// all of them to see both the rigidity statement and the E/u confound that
/// rides with it.  The paper publishes no uncertainty on any of them.
struct Chang26SpeciesEffRow {
  const char* nucleus;   ///< "2D" | "3He" | "4He" | "7Li" | "9Be" | "12C" | "16O"
  int a;                 ///< mass number
  int z;                 ///< proton number
  double e_ion_gev;      ///< beam energy [GeV/nucleon], Fig. 2's legend
  double efficiency;     ///< global tagging efficiency x acceptance, p. 4
};
const std::vector<Chang26SpeciesEffRow>& chang26_species_efficiency();

/// One (vector meson, Q^2 floor) row of the 2026-09-04 PHOTOPRODUCTION scan:
/// the same eSTARlight configuration as `estarlight_li6_coherent()` (same
/// commit 939b11a24499398392d959db81c7502aeec91046, same beams, same seed
/// 5574531, same BREAKUP_MODE / QUANTUM_GLAUBER), re-run with MIN_GAMMA_Q2
/// lowered.  2e5 events per default-density row, 1e5 per measured-radius row.
///
/// WHY IT EXISTS.  `estarlight_li6_coherent()`'s `sigma_nb` is the
/// 0.1 < Q^2 < 100 GeV^2 window, which is arXiv:2511.05638's ACCEPTANCE-STUDY
/// kinematic range copied verbatim -- not a physics window.  Q^2 < 0.1 GeV^2
/// was absent from the whole chain, and it is most of the rate: 6.75x for
/// J/psi, 21.7x for phi, 35.2x for rho.  (The one coherent vector-meson
/// measurement this tree cites, STAR's rho0 arXiv:2204.01625, is
/// ultraperipheral -- Q^2 ~ 0.  Whether photoproduction is where coherent
/// J/psi is "normally" measured is not sourced here; the rate is.)  `q2_floor_gev2 = 0` means NO lower cut, i.e.
/// eSTARlight's own kinematic limit Q^2_min = (m_e E_gamma)^2 /
/// (E_e (E_e - E_gamma)) ~ 1e-9 GeV^2 -- the whole quasi-real region.
///
/// WHAT DOES NOT CHANGE.  The |t| slope: 38.9 -> 38.8 GeV^-2 for J/psi
/// between the 0.1 floor and no floor, 0.4 %, so the recoil pT spectrum the
/// far-forward acceptance cuts on is the same sample.  <W> falls from 32.2 to
/// 30.2 GeV (J/psi), which by arXiv:2511.05638's own Fig. 2 moves the tagging
/// efficiency UP, not down -- but see `COHERENT_JPSI_EFF_IR8_LI7`: no point
/// of that efficiency was ever evaluated below Q^2 = 0.1.
struct EstarlightLi6Q2Row {
  const char* vm;            ///< "jpsi" | "phi" | "rho"
  double q2_floor_gev2;      ///< MIN_GAMMA_Q2; 0 = no floor (kinematic limit)
  double sigma_nb;           ///< default density R_G = 2.1805 fm
  double b_fit;              ///< MLE slope over 0 < |t| < 0.10 GeV^2, same run
  double sigma_rmeas_nb;     ///< same floor at R = 2.589 fm (patched build)
  double b_rmeas;            ///< MLE slope of that run
  /// <W> = sqrt(m_p^2 + 2 m_p E_gamma - Q^2) averaged over generated
  /// events, the W of arXiv:2511.05638's own Eq. (2) [GeV]
  /// (estarlight_li6.md sec. 2f).  In code since 2026-09-04 because the
  /// efficiency transfer argument -- both the Q^2-floor leg and the
  /// beam-energy leg -- runs on it, and it lived only in prose.
  double w_mean_gev;
};
const std::vector<EstarlightLi6Q2Row>& estarlight_li6_q2_floors();

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
  /// The DEFORMATION mechanism's one number.  A SCENARIO, and (open item O4,
  /// closed 2026-09-04) one whose shipped value is NOT a 6Li number.
  ///
  /// DELTA B, DEFINED HERE AND NOWHERE ELSE.  Before 2026-09-04 this header
  /// used the symbol Delta B in three docstrings and never defined it, which
  /// left the label "Delta B_0/B" ambiguous at the level of a sign
  /// (docs/open_items/run_2026-09-02/design_G_cluster_config.md sec. 2.8 lists
  /// the three defensible readings).  The definition, once:
  ///
  ///     |F_m(|t|, Phi)|^2 = exp(-|t| [B + Delta B_m cos 2(Phi - Phi_S)])
  ///
  /// with Phi the azimuth of the momentum transfer about the beam, Phi_S the
  /// spin azimuth, B = `slope_b` the phi-averaged slope and Delta B_m the
  /// cos 2Phi coefficient OF THE SLOPE, in GeV^-2.  For a Gaussian transverse
  /// profile that is exactly Delta B_m = delta_m / 2 with
  /// delta_m = <x^2> - <y^2> per nucleon in state m (design (G7)), and
  /// expanding to first order in Delta B |t| against the anchor's
  /// 1 + 2 a_2 cos 2Phi normalization gives
  ///
  ///     a_2(m) = -(Delta B_m / 2) |t| = -(delta_m / 4) |t|
  ///
  /// which is `a2_m_state` exactly.  `delta_b_m(m)` is the single code home of
  /// Delta B_m and `slope_at_azimuth` is the slope above; nothing else in the
  /// tree recomputes either.
  ///
  /// SO WHAT eps_b0 IS: eps_b0 = delta_{+-1} / B, i.e.
  ///
  ///     eps_b0 = +2 Delta B_{+-1} / B = -1 x Delta B_0 / B.
  ///
  /// The pre-2026-09-04 label "relative slope modulation Delta B_0/B of the
  /// m = 0 state" was therefore off BY A SIGN, not by a factor 2; the code was
  /// right and the label was wrong.  eps_b0 < 0 with P_zz > 0 gives a_2 > 0.
  ///
  /// eps_b0 AND `slope_b` ARE NOT INDEPENDENT.  Every observable here uses the
  /// PRODUCT delta_{+-1} = eps_b0 * slope_b (a2_deformation is
  /// -(P_zz/4) eps_b0 B |t|), and delta is the physics -- a transverse
  /// second-moment difference, fixed by the target's quadrupole and NOT by how
  /// steep its form factor is.  So a `slope_b` scan over the band {40, 60} at
  /// FIXED eps_b0 moves a_2 by +-20 % for no physical reason.  Band delta (or
  /// eps_b0 and slope_b together, anticorrelated), never eps_b0 alone.
  ///
  /// THE DEFAULT IS A DEUTERON-SIZED NUMBER, AND IT IS 11.4x TOO BIG FOR 6Li.
  /// Inverting the map (`quadrupole_from_a2_slope`, cluster_config.hpp) on
  /// eps_b0 = -0.08 at B = 50 gives delta_{+-1} = -4.0 GeV^-2, i.e. an implied
  /// 6Li CHARGE quadrupole of -0.9345 fm^2 and a_2(+-1, |t| = 0.3) = +0.300 --
  /// larger in magnitude than the DEUTERON's own digitized -0.26 at the same
  /// |t|, for a nucleus whose measured quadrupole (`LI6_QUADRUPOLE_FM2` =
  /// -0.0818 fm^2) is 3.5x SMALLER than the deuteron's.  The measured 6Li
  /// band, from `quadrupole_band_fm2()` through `a2_from_quadrupole` at this
  /// B = 50 (docs/open_items/run_2026-09-03/phase_C_numbers.md sec. C4):
  ///
  ///     measured    Q = -0.0818 fm^2  ->  eps_b0 = -0.0070
  ///     GFMC        Q = -0.20(6) fm^2 ->  eps_b0 = -0.0171 (-0.0120..-0.0223)
  ///     alpha+d VMC Q = -0.6154 fm^2  ->  eps_b0 = -0.0527
  ///
  /// so the honest 6Li band is -(0.0070 .. 0.0527) -- the factor 7.5
  /// quadrupole budget of sec. C1, no wider and no narrower -- and the shipped
  /// -0.08 sits 1.52x ABOVE the top of it.  The old band -(0.04 .. 0.13) is a
  /// DEUTERON band: it was set from the deuteron's own eps_b0 = +0.105 at the
  /// deuteron's own B_d = 33.1 GeV^-2, at the opposite SIGN, and never rescaled
  /// to 6Li.  Exactly, because "does not overlap" would be wrong: the old band
  /// contains the alpha+d model row 0.0527 and NOTHING else of 6Li -- the GFMC
  /// 0.0171 and the measured 0.0070 are both below its floor 0.04 (T10b).
  ///
  /// AUTHOR DECISION 2026-09-04, recorded not taken silently: the default
  /// STAYS -0.08, because it is pinned bit-for-bit in
  /// validation/reference/*.json and this run may not move a reference gate.
  /// The cost is stated in one line: every GENERATED coherent tensor number
  /// (a2_deformation, cos2phi_coefficient, the sampled azimuth) is 11.4x the
  /// measured-quadrupole expectation.  What that does NOT any longer imply
  /// is the |t| ceiling: since 2026-09-04 `COHERENT_T_MAX_DEFAULT` = 0.2 is
  /// justified by the ANCHOR RANGE, which does not move with eps_b0, and the
  /// positivity edge is stated as the contingent second reason and DERIVED
  /// by `t_positivity_edge` (|t| = 0.245 GeV^2 here at P_zz = -2; 2.80 at
  /// the measured-quadrupole -0.0070).  Nothing published
  /// from this channel may quote a single eps_b0 row: band it, and say which
  /// quadrupole the row assumes.
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

  /// Delta B_m [GeV^-2], the cos 2(Phi - Phi_S) coefficient of the m state's
  /// |F|^2 slope, in the convention DEFINED on `eps_b0` above and nowhere
  /// else: Delta B_{+-1} = eps_b0 B / 2 and Delta B_0 = -eps_b0 B.  `m` is 0
  /// or +-1 (anything non-zero reads as +-1, as in `a2_m_state`).
  double delta_b_m(int m) const;
  /// The m state's slope at azimuth `phi_rel` = Phi - Phi_S [rad]:
  /// B + Delta B_m cos 2 phi_rel [GeV^-2].  This IS the definition; every
  /// a_2 below is its first-order expansion.
  double slope_at_azimuth(double phi_rel, int m) const;

  /// Ensemble cos 2phi' coefficient a_2 of the coherent yield at |t| from the
  /// slope-modulation (deformation) mechanism.  Per pure m state
  /// a_2(m) = -(Delta B_m/2)|t| with a_2(0) = -2 a_2(+-1) (`delta_b_m`); the
  /// population average is a_2 = -(P_zz/4) eps_b0 B |t|.  Exact only as
  /// |t| -> 0, where the m-state phi-averaged rates are equal (see
  /// RATE_WEIGHT_SYST).
  double a2_deformation(double t_abs, double pzz) const;
  /// The cos 2phi_t COEFFICIENT of 1 + c_2 cos 2(phi_t - phi_S) from the
  /// deformation mechanism: c_2 = 2 a_2 (the anchor's Eq. (9) normalization).
  double cos2phi_coefficient_deformation(double t_abs, double pzz) const;
  /// <a_2> of the RP-tagged sample in the equal-rate, linear-in-|t|
  /// approximation.
  double a2_tagged(double pt_cut, double pzz) const;
  /// a_2 of a PURE m state: a_2(0) = -2 a_2(+-1), normalized so that the
  /// population average over p_m reproduces `a2_deformation`.  This is the
  /// m-state relation the tests pin, and it is identically
  /// -(`delta_b_m(m)`/2) |t| -- the operative definition of Delta B.
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
  /// zero at |t| = 0.495 for P_zz = +1 and at |t| = 0.245 for P_zz = -2 --
  /// `t_positivity_edge` below is that zero in closed form.
  ///
  /// THIS IS A GUARD, NOT THE DERIVATION OF `COHERENT_T_MAX_DEFAULT`
  /// (2026-09-04, D5).  It catches an author who raises `t_max` or `eps_b0`
  /// past the point where the TRUNCATED weight the sampler actually uses
  /// stops being a density.  Both failure modes are real and nothing else
  /// catches them, so it must keep throwing; but the ceiling's own stated
  /// reason is the anchor range (`COHERENT_T_MAX_DEFAULT`), because this one
  /// stops binding as soon as eps_b0 is corrected.
  double positivity_margin(double t_max, double pzz) const;

  /// The |t| [GeV^2] at which `positivity_margin(|t|, pzz)` reaches zero --
  /// i.e. where the linear azimuthal weight 1 + c_2 cos 2(phi - phi_S) first
  /// stops being a density.  DERIVED, so the number moves with `eps_b0` and
  /// is not repeated in four docstrings (2026-09-04, D5).
  ///
  /// With c_2 = A|t| + C, A = -(P_zz/2) eps_b0 B and C = amp P_zz, the edge
  /// |c_2| = 1 is at
  ///
  ///     |t|_pos = (1 - sign(A) C) / |A|
  ///
  /// which for the shipped signs (eps_b0 < 0, amp > 0, so A and C share a
  /// sign) is 2 (1/|P_zz| - amp) / (|eps_b0| B).  Measured against the
  /// sec. C4 quadrupole budget at B = 50, amp = 0.01:
  ///
  ///     eps_b0    P_zz = -2   P_zz = +1
  ///     -0.08     0.2450      0.4950     (shipped -- an EXACT input)
  ///     -0.0527   0.3719      0.7514     (alpha+d model Q, rounded)
  ///     -0.0171   1.1462      2.3158     (GFMC Q, rounded)
  ///     -0.0070   2.80        5.6571     (MEASURED Q, rounded)
  ///
  /// THE LAST THREE eps_b0 ARE ROUNDED, AND THE ROW IS ARITHMETIC ON THE
  /// ROUNDED VALUE.  On the DERIVED band (`ClusterConfigSampler::
  /// quadrupole_band_fm2()` through `a2_from_quadrupole`, at B = 50,
  /// amp = 0.01, 2026-09-05) the eps_b0 are -0.0526846 / -0.0171207 /
  /// -0.0070024 and the P_zz = -2 edges are 0.3720 / 1.1448 / 2.7990.  So
  /// the measured-Q row is "2.80 GeV^2", to the three figures the rounded
  /// input carries -- NEVER "2.8000", which is four figures of arithmetic on
  /// a two-figure input.  `tests/test_coherent.cpp` asserts the DERIVED
  /// 2.7990 on the derived eps_b0, which is where that digit comes from.
  ///
  /// Read the bottom row before using this as a ceiling: it is 9.3x outside
  /// the |t| <= 0.30 the deformation input is digitized over
  /// (`mantysaari_a2_deuteron`), which is exactly why
  /// `COHERENT_T_MAX_DEFAULT` is justified by the ANCHOR RANGE and not by
  /// this function.
  ///
  /// Degenerate cases, stated rather than silently returned: A == 0 (P_zz,
  /// eps_b0 or B zero) means c_2 is constant, so the edge is +infinity when
  /// |C| < 1 and 0 when it is not; and when |C| >= 1 already at |t| = 0 the
  /// weight is not a density anywhere and the edge is 0.
  double t_positivity_edge(double pzz) const;
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
