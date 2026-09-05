#!/usr/bin/env python3
"""Open item O5: does 6Li's tensor a_2 survive EIC statistics, and does the
P_zz flip separate it from the linearly-polarized-photon cos 2phi background?

Numbers and prose: `docs/open_items/run_2026-09-03/phase_C_numbers.md` section
C2, and `docs/OPEN_ITEMS_SOLUTIONS.md` section 11.  Pinned in
`python/tests/test_o5_reach.py` (which imports this module) and in
`tests/test_coherent.cpp` T10a / T10c (the eSTARlight tables only).

WHAT THIS IS AND IS NOT.  It is a closed-form statistical-reach estimate, and
every number it multiplies is read from code, none retyped:

  1. `a2_from_quadrupole` (cluster_config.hpp), the closed-form
     quadrupole -> a_2 map, evaluated at the MEASURED `LI6_QUADRUPOLE_FM2`;
  2. `estarlight_li6_coherent()` and `estarlight_li6_q2_floors()`
     (coherent.hpp), the eSTARlight unpolarized coherent rate baseline at
     e 10 GeV x 6Li 99.5 GeV/u -- the second table is the 2026-09-04
     PHOTOPRODUCTION scan, the Q^2 < 0.1 GeV^2 region the first one omits;
  3. `Scenario` / `tensor_flip_plan` (sampler.hpp / bookkeeping.hpp), the
     repository's own programme luminosity and P_zz run plan;
  4. `COHERENT_JPSI_EFF_IR8_LI7` (coherent.hpp), arXiv:2511.05638's 7Li
     coherent-J/psi far-forward efficiency, used as a stand-in for 6Li;
  5. `tagging_optics` / `yr_optics` (spectator.hpp), whose `lumi_fraction`
     is the luminosity price of a de-squeezed far-forward working point.

(Before 2026-09-04 items 4 and 5 were absent or a literal here, and this
docstring claimed "three in-tree inputs and NOTHING ELSE" while the single
largest multiplier in the chain -- the efficiency -- was typed into this file.
It was not true; it is now.)

It is NOT a dipole-model calculation.  `a2_from_quadrupole` carries only the
target's quadrupole moment through the deuteron's own published |t| dependence
([Mant24] Eq. (9)); it reproduces that published a_2(m = +-1) to 8 % with zero
free parameters, but it does not contain a Good-Walker average, an amplitude,
saturation, or their uncertainties.  Everything below inherits that.  Read it
as "what a MEASURED quadrupole of this size buys, if the map holds", never as
a prediction of the coherent cross section's azimuthal structure.

THE STATISTICS, in one line.  For a normalized azimuthal density
f(phi) = (1 + c_2 cos 2phi)/2pi the per-event Fisher information on c_2 is 1/2,
which is exactly the sqrt(2/N) of `estimators::cos2phi_fit_err`.  The tensor
signal is c_2(|t|) = 2 a_2(m=+-1, |t|) P_zz = 2 kappa |t| P_zz -- the factor 2
is [Mant24]'s Eq. (9) normalization (`cos2phi_coefficient_deformation`), and
the P_zz is exact because a_2(0) = -2 a_2(+-1) makes the population average
sum_m p_m a_2(m) = P_zz a_2(+-1) identically (`a2_m_state`).  So

    delta(kappa) = 1 / sqrt(2 <P_zz^2> N <t^2>),   S = kappa / delta(kappa),

with <t^2> taken over the truncated exponential dN/d|t| ~ B exp(-B|t|) that
eSTARlight's own |t| distributions are (section 2c of estarlight_li6.md: a
Gaussian form factor has no diffractive minimum, so the exponential is exact
in that model, and the fitted <|t|> = 0.0254 GeV^2 reproduces 1/B = 0.0257).

WHICH <P_zz^2>.  TWO of them, and the report prints both.  `PZZ2_OPTIMAL`
= 0.90 is `tensor_flip_plan`'s luminosity-weighted <P_zz^2>, i.e. what an
optimal combination of the two fills would deliver if the photon-polarization
background were absent.  `PZZ2_DIFFERENCE` = 0.81 is the effective value of
the BACKGROUND-IMMUNE two-fill difference -- the estimator section C2.5's
separation argument actually requires, 5.4 % worse in delta.  Every headline
below is quoted on the DIFFERENCE; the tables carry both.

ASSUMPTIONS, all of them, and where each one is stated in the output:
  * transverse (vertical) tensor-polarized 6Li, theta_S = 90 deg -- the signal
    is proportional to sin^2(theta_S) (`delta_perp_analytic_fm2`), so a
    longitudinal axis gives EXACTLY ZERO;
  * uniform acceptance in the azimuth of the momentum transfer, and a |t|
    resolution small enough not to smear it (NO in-tree number exists for the
    coherent-VM p_T resolution, so no smearing dilution is applied -- the
    formula for one is printed);
  * the detection efficiency is flat in |t|, in that azimuth AND IN Q^2 --
    the last of those is an extrapolation: `COHERENT_JPSI_EFF_IR8_LI7` was
    measured on a 0.1 < Q^2 < 100 GeV^2 sample and is UNMEASURED anywhere in
    this tree below Q^2 = 0.1, which is where most of the rate is;
  * ... and it is NOT flat in BEAM ENERGY, which until 2026-09-04 this list
    did not say.  0.1775 is a 7Li number at that nucleus's own TOP energy,
    18 x 117.9 GeV/u (<W> = 43.2 in this tree's own row for those beams);
    it is applied to a 6Li sample at 10 x 99.5 (<W> = 30.2).  Direction UP,
    size x1.1224-1.1582 from arXiv:2511.05638's own 3He energy scan
    (`chang26_he3_energy_scan`).  Section 3b(a);
  * ... and it is a 7Li number applied to 6Li.  THE SIZE AND THE DIRECTION OF
    THAT SUBSTITUTION ARE NOT ESTABLISHED, and until 2026-09-04 this file
    said they were (x1.1603-1.2192 UP, "at fixed rigidity").  Fixing the
    rigidity does not fix E/u or Z, and read off the four
    `chang26_species_efficiency()` entries that DO share a beam energy the
    substitution is x0.9932-1.3327 -- it STRADDLES 1.  Section 3b(b), and it
    is the whole of why the band's TOP is a span rather than an edge;
  * ... and it is a RECOIL-NUCLEUS number, so the chain sigma x BR(l+l-) x
    eps_det has NO central-detector acceptance and NO reconstruction
    efficiency for the DECAY LEPTONS.  Until 2026-09-04 this list did not
    say that either.  Direction DOWN.  The GEOMETRIC part is bounded and is
    NOT the problem (0.9940 at this sample's <W>, never below 0.89 anywhere
    in the accessible W range -- `lepton_pair_acceptance`, on this tree's own
    `Scenario::eta_max`); the per-lepton reconstruction efficiency is
    UNBOUNDED HERE, and the 0.95/track used at the band's low end is the
    sibling ../PolarizedLithiumSim's own stand-in.  Section 3b(c);
  * the two previous items PARTLY CANCEL (x1.1224 against x0.8971 = 1.007),
    which is why neither may be quoted without the other;
  * the programme luminosity carries NO energy dependence: 10 fb^-1/u is a
    scenario number with no Li source at all (`needs_survey.md` sec. 3.8), so
    no penalty is applied for running at 10 x 99.5 rather than at a top
    energy, and none can be -- direction and size both unknown here;
  * the coherent sample is identified as coherent (breakup vetoed);
  * the far-forward working point is the Yellow Report one, `lumi_fraction`
    = 1, NOT LiPolGen's own de-squeezed 6Li tagging point (0.0781 at
    10 x 100).  The price of the other choice is printed.
"""

from __future__ import annotations

import math
import sys

import lipolgen as _l

# --------------------------------------------------------------- the inputs

#: 6Li mass number -- an e+A luminosity quoted PER NUCLEON divides by it.
A_LI6 = 6

#: Detection efficiency used for the "reconstructed" column.  Read from
#: `coherent.hpp`, which carries the provenance and the caveats: it is
#: arXiv:2511.05638's GLOBAL far-forward tagging efficiency x acceptance for
#: coherent J/psi off 7Li at 18 x 118 through the IR-8 secondary focus.  It is
#: a STAND-IN, applied to 6Li and to all three vector mesons, and every table
#: is also printed at efficiency 1 so the reader can see what it costs.
EPS_DET_CHANG = _l.COHERENT_JPSI_EFF_IR8_LI7

#: The Q^2 windows the report prices, in the order it prints them.  "Q2>0.7"
#: is LiPolGen's own generator window; "Q2>0.1" is arXiv:2511.05638's
#: acceptance-study range, which is what the whole chain used to stop at;
#: "Q2>0.01" is a deliberately conservative photoproduction floor; "Q2>0" is
#: no floor at all, i.e. eSTARlight's kinematic m_e limit.
Q2_FLOORS = (0.1, 0.01, 0.0)


def _floor_label(q2: float) -> str:
    return "Q2>0" if q2 == 0.0 else ("Q2>%g" % q2)


# ------------------------------------------- the efficiency chain, corrected
#
# TWO FACTORS WERE MISSING FROM IT UNTIL 2026-09-04, AND THEY PULL OPPOSITE
# WAYS.  Both are stated here, both are priced from `coherent.hpp` tables, and
# NEITHER is folded silently into the headline: section 3b prints the band.

#: The beam configuration the whole rate chain is priced at -- LiPolGen's own
#: e 10 GeV x 6Li 99.5 GeV/u, which is the eSTARlight run of
#: estarlight_li6.md sec. 2a/2f.  Read from `default_configs`, not retyped.
BEAM_LI6 = _l.default_configs("6Li")[1]
#: 6Li's own TOP energy, 18 x 137.5 -- the fixed-RIGIDITY partner of the beam
#: `COHERENT_JPSI_EFF_IR8_LI7` was measured at.
BEAM_LI6_TOP = _l.default_configs("6Li")[2]
#: ... and that beam itself: 7Li's own top energy, 18 x 117.9, which
#: arXiv:2511.05638 writes "18 x 118".
BEAM_LI7_TOP = _l.default_configs("7Li")[2]

# The constant beside the efficiency must be the beam `default_configs` gives,
# or the correction below is being applied against the wrong reference.
assert abs(BEAM_LI7_TOP.ion_momentum_per_nucleon
           - _l.COHERENT_JPSI_EFF_IR8_LI7_E_ION_GEV) < 0.05


def he3_energy_slopes() -> list:
    """d ln(eff) / d ln(E_ion) between consecutive `chang26_he3_energy_scan()`
    points: [-0.8656, -0.6807].

    NEGATIVE -- the far-forward efficiency RISES as the beam energy falls, and
    that is a MEASUREMENT (arXiv:2511.05638 sec. V.B, e+3He at 18 x 183,
    10 x 100 and 5 x 41: 32.23 %, 54.38 %, 99.77 %), not an argument.  The two
    slopes differ because the efficiency is bounded above by 1 and must
    flatten; both are carried, as a band, rather than fitted to one power law.
    """
    rows = list(_l.chang26_he3_energy_scan())
    return [math.log(b.efficiency / a.efficiency)
            / math.log(b.e_ion_gev / a.e_ion_gev)
            for a, b in zip(rows, rows[1:])]


def beam_energy_scaling(e_from: float = None, e_to: float = None) -> tuple:
    """(lo, hi) factor on the far-forward efficiency for moving the ION beam
    from `e_from` to `e_to` GeV/u, on the two local power laws above.

    Defaults are the step the O5 chain actually needs and has never applied:
    117.9 (where 0.1775 was measured) -> 99.5 (where it is used).  Both
    factors exceed 1, so 0.1775 at 10 x 99.5 is CONSERVATIVE.

    THIS IS NOT ARXIV:2511.05638'S 1.687.  That is their 183 -> 100 GeV/u
    ratio, a step 3.56x larger in ln E than the one required here; carrying it
    whole would be scaling by the wrong lever arm.
    """
    if e_from is None:
        e_from = _l.COHERENT_JPSI_EFF_IR8_LI7_E_ION_GEV
    if e_to is None:
        e_to = BEAM_LI6.ion_momentum_per_nucleon
    f = [math.exp(s * math.log(e_to / e_from)) for s in he3_energy_slopes()]
    return (min(f), max(f))


#: The one beam energy four of the seven `chang26_species_efficiency()`
#: entries share [GeV/u]: 2D, 4He, 12C and 16O.  6Li's own fixed-rigidity
#: energy is Z/A x 275 = 137.5, so those four bracket A = 6 with NO energy
#: step -- which is why `species_scaling_same_energy()` exists.
SPECIES_E_SAME_GEV = 137.0


def species_list_is_not_a_species_lever() -> dict:
    """THE RETRACTION, as arithmetic.  What fixing the rigidity does and does
    not eliminate in `chang26_species_efficiency()`.

    Until 2026-09-04 five sites in this tree said "every entry is at the same
    rigidity, so its A-ordering is a species lever on its own".  That is a NON
    SEQUITUR.  Holding R = A E/Z fixed eliminates rigidity and nothing else:
    at fixed R the per-nucleon energy is E/u = R Z/A and the total beam
    momentum is p_z = Z R, so BOTH still vary down the list.  E/u runs 118
    (7Li) to 183 (3He) GeV/u, Z runs 1 to 8, p_z runs 274 to 2192 GeV.

    AND THE CONFOUND IS THE SIZE OF THE EFFECT.  The chain's OTHER leg
    measures d ln(eff)/d ln(E_ion) = -0.68 to -0.87 (`he3_energy_slopes`), so
    the list's own x1.55 spread in E/u is worth x1.35-1.46 in the efficiency
    -- as large as the whole x1.16-1.22 "species gain" that was read off it.
    The list is a JOINT (A, Z, E/u) lever; nothing in the table decomposes it,
    and this chain never tested the decomposition it asserted.
    """
    rows = list(_l.chang26_species_efficiency())
    eu = [r.e_ion_gev for r in rows]
    return dict(
        e_per_u=(min(eu), max(eu)),
        e_per_u_spread=max(eu) / min(eu),
        z=(min(r.z for r in rows), max(r.z for r in rows)),
        p_z=(min(r.a * r.e_ion_gev for r in rows),
             max(r.a * r.e_ion_gev for r in rows)),
        eff_worth_of_e_spread=tuple(
            sorted(math.exp(s * math.log(min(eu) / max(eu)))
                   for s in he3_energy_slopes())),
        n_at_same_energy=sum(1 for r in rows
                             if r.e_ion_gev == SPECIES_E_SAME_GEV))


def _interp_forms(a0: float, a1: float, e0: float, e1: float,
                  target: float = 6.0, anchors=()) -> dict:
    """The closed forms this file interpolates a species list with, between
    (a0, e0) and (a1, e1): log-linear in A, linear in A, log-linear in the
    Gaussian radius R_G = 1.2 A^(1/3), and this repository's OWN acceptance
    model `CoherentScenario::tag_acceptance` = exp(-B pT_cut^2) with
    B = `gaussian_slope` and pT_cut inverted from each named anchor.

    NONE of them is a fit and none has an error bar; they are four (or five)
    readings of two or three points, and their SPREAD is the only uncertainty
    statement available."""
    r_g = lambda a: 1.2 * a ** (1.0 / 3.0)
    b = lambda a: _l.gaussian_slope(r_g(a))
    w_a = (target - a0) / (a1 - a0)
    w_r = (r_g(target) - r_g(a0)) / (r_g(a1) - r_g(a0))
    out = dict(
        log_a=math.exp(math.log(e0) + (math.log(e1) - math.log(e0)) * w_a),
        lin_a=e0 + (e1 - e0) * w_a,
        log_r=math.exp(math.log(e0) + (math.log(e1) - math.log(e0)) * w_r))
    for name, a_anchor, e_anchor in anchors:
        out[name] = math.exp(-b(target) * (math.log(1.0 / e_anchor)
                                           / b(a_anchor)))
    return out


def species_scaling() -> dict:
    """The 7Li -> 6Li substitution as the chain USED to read it: interpolating
    arXiv:2511.05638's 4He (29.42 %, at 137 GeV/u) onto its 7Li (17.75 %, at
    118 GeV/u).  x1.1603-1.2192.

    THIS READING IS CONFOUNDED AND IS KEPT ONLY SO THE CONFOUND IS VISIBLE.
    Its two anchors differ by 19 GeV/u in E/u, and by this chain's own energy
    slopes that difference is worth x1.13-1.17 on its own -- most of the
    x1.16-1.22 it returns.  `species_scaling_same_energy()` is the reading
    that has no energy step in it, and it is the one the band uses.
    See `species_list_is_not_a_species_lever()` for the retracted sentence.

    ALSO RETURNED: the three readings of the ONE place the pT-threshold form
    can be tested against the table, 3He/4He.  See `he3_he4_readings()` --
    the shipped reading applied 7Li's Z = 3 cut to a Z = 2 pair, and the
    "overstates by 12 %, read the low end" that came out of it is withdrawn.
    """
    rows = {r.nucleus: r for r in _l.chang26_species_efficiency()}
    e4, e7 = rows["4He"].efficiency, rows["7Li"].efficiency
    assert abs(e7 - _l.COHERENT_JPSI_EFF_IR8_LI7) < 1e-12
    forms = _interp_forms(4.0, 7.0, e4, e7,
                          anchors=(("pt_threshold", 7.0, e7),))
    r_g = lambda a: 1.2 * a ** (1.0 / 3.0)
    pt2 = math.log(1.0 / e7) / _l.gaussian_slope(r_g(7.0))
    out = dict(eps_li6=forms, pt_cut_gev=math.sqrt(pt2),
               factors={k: v / e7 for k, v in forms.items()},
               lo=min(forms.values()) / e7, hi=max(forms.values()) / e7)
    out.update(he3_he4_readings())
    return out


def he3_he4_readings() -> dict:
    """THREE readings of one table, and the chain picked the one that yielded
    a directional instruction.

    The model is eps = exp(-B(A) pT_cut^2).  The criterion arXiv:2511.05638
    states is "within a safe distance FROM THE BEAM", i.e. a cut on the ANGLE,
    so pT_cut = theta p_z = theta Z R -- IT CARRIES THE CHARGE.

      (1) pT_cut ~ Z.  3He and 4He are Z = 2 and the pT_cut is inverted from
          7Li's Z = 3, so the pair's pT_cut^2 is (2/3)^2 of it.  Prediction
          1.0967 against the measured 1.0955: 0.1 %.
      (2) The SHIPPED reading applied 7Li's Z = 3 cut unchanged to the Z = 2
          pair -- its own comment said "same Z" without noticing it was the
          WRONG Z -- and got 1.2309, "overstating by 12 %".  That, and only
          that, is where the instruction "read the low end" came from.  It is
          withdrawn.
      (3) pT_cut Z-INDEPENDENT.  Over the LIGHT half of the list this fits the
          ABSOLUTE values far better than (1) does: pred/meas 1.003 (2D),
          1.034 (4He), 1.047 (9Be) against (1)'s 1.95, 2.00 and 0.21.  The one
          misfit among the light entries is 3He, by 16 % -- and 3He is the one
          nucleus in the list at E/u = 183 rather than the family's 137, for
          which this chain's own energy scan predicts an 18-22 % shortfall.
          IT IS NOT A DESCRIPTION OF THE WHOLE LIST: the form degrades with A
          and misses 12C by 32 % and 16O by a factor 3.1.  It is a LOCAL
          interpolation, and A = 6 sits inside the range where it works.

    Reading (1) says the form is right, (3) says it is right except where the
    energy confound bites, (2) says it overstates.  Nothing in the table picks
    between them, so nothing licenses a directional instruction either way.
    """
    rows = {r.nucleus: r for r in _l.chang26_species_efficiency()}
    r_g = lambda a: 1.2 * a ** (1.0 / 3.0)
    b = lambda a: _l.gaussian_slope(r_g(a))
    e7 = rows["7Li"].efficiency
    pt2 = math.log(1.0 / e7) / b(7.0)
    meas = rows["3He"].efficiency / rows["4He"].efficiency
    pred_z = math.exp(-(b(3.0) - b(4.0)) * pt2 * (2.0 / 3.0) ** 2)
    pred_wrong = math.exp(-(b(3.0) - b(4.0)) * pt2)
    absolute = {n: math.exp(-b(float(r.a)) * pt2) / r.efficiency
                for n, r in rows.items()}
    absolute_z = {n: math.exp(-b(float(r.a)) * pt2 * (r.z / 3.0) ** 2)
                  / r.efficiency for n, r in rows.items()}
    shortfall = tuple(sorted(
        math.exp(s * math.log(rows["3He"].e_ion_gev / SPECIES_E_SAME_GEV))
        for s in he3_energy_slopes()))
    return dict(he3_he4_measured=meas,
                he3_he4_predicted_z=pred_z,
                he3_he4_ratio_z=pred_z / meas,
                he3_he4_predicted_wrong_z=pred_wrong,
                he3_he4_overstatement_wrong_z=pred_wrong / meas,
                absolute_pred_over_meas=absolute,
                absolute_pred_over_meas_z=absolute_z,
                he3_energy_shortfall=shortfall)


def species_scaling_same_energy() -> dict:
    """The 7Li -> 6Li substitution read off the four entries that share a beam
    energy: 2D, 4He, 12C and 16O, all at 137 GeV/u.  x0.99 to x1.33.

    THIS IS THE READING WITH NO ENERGY STEP IN IT.  6Li's own fixed-rigidity
    energy is Z/A x 275 = 137.5, so 4He and 12C bracket A = 6 at the SAME
    energy, and the interpolation is a pure A interpolation for the first
    time.  Five closed forms (`_interp_forms`, with the pT-threshold anchored
    at each end): 0.1763 to 0.2366 against 7Li's 0.1775, i.e. x0.9932 to
    x1.3327.  IT STRADDLES 1: on this reading the direction of the 7Li -> 6Li
    substitution is not established, let alone its size.

    The band's TOP is quoted as the span this leaves, and the two log forms
    put it below 3 sigma while linear-in-A puts it above.  lin_a is the crude
    end of that -- a straight chord across a 4.6x fall in efficiency between
    A = 4 and A = 12 -- but it is one of the same four forms the confounded
    reading was quoted from, and dropping it now would be choosing the answer.
    """
    rows = {r.nucleus: r for r in _l.chang26_species_efficiency()}
    same = [r for r in _l.chang26_species_efficiency()
            if r.e_ion_gev == SPECIES_E_SAME_GEV]
    assert [r.nucleus for r in same] == ["2D", "4He", "12C", "16O"], same
    e4, e12 = rows["4He"].efficiency, rows["12C"].efficiency
    e7 = rows["7Li"].efficiency
    forms = _interp_forms(4.0, 12.0, e4, e12,
                          anchors=(("pt_threshold@4He", 4.0, e4),
                                   ("pt_threshold@12C", 12.0, e12)))
    return dict(eps_li6=forms, at_energy_gev=SPECIES_E_SAME_GEV,
                nuclei=[r.nucleus for r in same],
                factors={k: v / e7 for k, v in forms.items()},
                lo=min(forms.values()) / e7, hi=max(forms.values()) / e7)


# --- ... and the factor that has never been in the chain AT ALL -------------
#
# arXiv:2511.05638's number is the fraction of scattered NUCLEI inside the
# far-forward acceptance, and its own text (p. 4) says it "only accounts for
# the acceptance effect and does not incorporate the efficiencies of the
# detector.  Additionally, we did not account for the efficiency and
# acceptance of the reconstructed distribution".  So
# sigma x BR(l+l-) x eps_recoil has NO central-detector acceptance and NO
# reconstruction efficiency for the decay leptons.  Direction: DOWN.

#: J/psi mass [GeV] -- PDG, the value eSTARlight itself uses
#: (starlightconstants.h `JpsiMass`); the branchings already in
#: `EstarlightLi6Row` come from the same header.
M_JPSI_GEV = 3.0969
#: proton mass [GeV], for the W <-> E_gamma map.
M_PROTON_GEV = 0.938272


def jpsi_lab_rapidity(w_gev: float, e_ion_per_nucleon: float = None) -> tuple:
    """(y_lab of the J/psi, E_gamma in the lab) for coherent production at
    |t| -> 0 off a nucleon of per-nucleon energy `e_ion_per_nucleon`.

    Exact two-body kinematics: boost the gamma-N centre of mass (which the
    photon and the moving nucleon define) by its own rapidity and subtract the
    VM's rapidity in it.  +z is the ION direction, so a negative y_lab means
    the J/psi leans toward the electron beam.  At <W> = 30.2 GeV and
    99.5 GeV/u this gives E_gamma = 2.29 GeV and y = -0.39: the J/psi is very
    nearly AT REST in the lab, which is the whole reason the next function
    comes out near 1.
    """
    if e_ion_per_nucleon is None:
        e_ion_per_nucleon = BEAM_LI6.ion_momentum_per_nucleon
    e_n = e_ion_per_nucleon
    p_n = math.sqrt(e_n * e_n - M_PROTON_GEV ** 2)
    e_gamma = (w_gev * w_gev - M_PROTON_GEV ** 2) / (2.0 * (e_n + p_n))
    y_cm = 0.5 * math.log((e_n + p_n)
                          / ((e_n - p_n) + 2.0 * e_gamma))
    e_star = (w_gev * w_gev + M_JPSI_GEV ** 2 - M_PROTON_GEV ** 2) / (2.0 * w_gev)
    lam = ((w_gev * w_gev - M_JPSI_GEV ** 2 - M_PROTON_GEV ** 2) ** 2
           - 4.0 * (M_JPSI_GEV * M_PROTON_GEV) ** 2)
    p_star = math.sqrt(lam) / (2.0 * w_gev)
    return (y_cm - math.log((e_star + p_star) / M_JPSI_GEV), e_gamma)


def lepton_pair_acceptance(w_gev: float, eta_max: float = None,
                           pt_min: float = None, schc: bool = True) -> float:
    """GEOMETRIC acceptance for BOTH decay leptons of a coherent J/psi.

    `eta_max` defaults to `Scenario::eta_max` = 3.5, this repository's OWN
    "crude central-detector acceptance" (documented there for the SCATTERED
    ELECTRON -- there is no decay-lepton acceptance anywhere in this tree, and
    that is the point of this function).  `pt_min` is optional and is NOT an
    in-tree number: 0.2 GeV is the sibling ../PolarizedLithiumSim's own
    `HfsModel(pt_min_track=0.2)` stand-in for a 1.7 T solenoid.

    The leptons are back to back at p* = sqrt(M^2/4 - m_l^2) = 1.548 GeV in
    the J/psi frame, so both are inside |eta| < eta_max exactly when the
    rest-frame emission satisfies |eta*| < eta_max - |y_Jpsi|, i.e.
    |cos theta*| < tanh(eta_max - |y_Jpsi|).  `schc` weights the decay by
    s-channel-helicity-conserving (3/8)(1 + cos^2 theta*) -- the transverse
    J/psi of diffractive photoproduction, which pushes leptons TOWARD the beam
    and so is the pessimistic choice; False is isotropic.

    ASSUMPTIONS, all optimistic-to-neutral and none of them a detector
    simulation: p_T(J/psi) ~ sqrt(|t|) ~ 0.16 GeV is neglected (it helps one
    lepton and hurts the other); no material, no magnetic-field sagitta cut
    beyond `pt_min`; and the result is a pure ACCEPTANCE -- the per-lepton
    tracking/PID/reconstruction efficiency is a SEPARATE factor and is
    UNBOUNDED anywhere in this tree.
    """
    if eta_max is None:
        eta_max = _l.Scenario().eta_max
    y_v, _ = jpsi_lab_rapidity(w_gev)
    x = math.tanh(max(eta_max - abs(y_v), 0.0))
    p_star = math.sqrt(0.25 * M_JPSI_GEV ** 2 - 0.000511 ** 2)
    if pt_min is not None and pt_min < p_star:
        x = min(x, math.sqrt(1.0 - (pt_min / p_star) ** 2))
    return (0.75 * x + 0.25 * x ** 3) if schc else x


#: Per-lepton reconstruction efficiency used ONLY to show what a plausible one
#: costs.  It is NOT an in-tree number and NOT a measurement: it is the
#: sibling ../PolarizedLithiumSim's own `HfsModel(eff_track=0.95)`, which that
#: file labels "stand-in" in as many words.  The pair costs its SQUARE.
EPS_TRACK_SIBLING_STANDIN = 0.95


def jpsi_w_mean(q2_floor: float = 0.0) -> float:
    """<W> [GeV] of the 6Li coherent-J/psi sample at that Q^2 floor, from
    `estarlight_li6_q2_floors()` (in code since 2026-09-04)."""
    for r in _l.estarlight_li6_q2_floors():
        if r.vm == "jpsi" and r.q2_floor_gev2 == q2_floor:
            return r.w_mean_gev
    raise KeyError(q2_floor)


def eps_det_chain() -> dict:
    """Every factor on `COHERENT_JPSI_EFF_IR8_LI7`, with its direction, and
    the band they leave.  This is the corrected chain; section 3b prints it.
    """
    beam_lo, beam_hi = beam_energy_scaling()
    sp_conf = species_scaling()              # the CONFOUNDED reading, kept
    sp = species_scaling_same_energy()       # the one with no energy step
    # the two-leg route: species (7Li -> 6Li at 6Li's own top energy), THEN
    # 6Li's own 137.5 -> 99.5 energy step.
    e2_lo, e2_hi = beam_energy_scaling(
        BEAM_LI6_TOP.ion_momentum_per_nucleon,
        BEAM_LI6.ion_momentum_per_nucleon)
    w = jpsi_w_mean(0.0)
    a_geom = lepton_pair_acceptance(w)
    a_geom_pt = lepton_pair_acceptance(w, pt_min=0.2)
    eps_pair_standin = EPS_TRACK_SIBLING_STANDIN ** 2
    return dict(
        beam_lo=beam_lo, beam_hi=beam_hi,
        # THE SPECIES LEG IS A SPAN THAT STRADDLES 1, not a gain.
        species_lo=sp["lo"], species_hi=sp["hi"], species=sp,
        species_confounded=sp_conf,
        species_confounded_lo=sp_conf["lo"], species_confounded_hi=sp_conf["hi"],
        two_leg_lo=sp["lo"] * e2_lo, two_leg_hi=sp["hi"] * e2_hi,
        # the band's TOP is BOTH legs at their most favourable -- but the
        # species leg is a SPAN, so the top is a span too.
        two_leg_top_min=sp["lo"] * e2_hi, two_leg_top_max=sp["hi"] * e2_hi,
        energy_from_li6_top=(e2_lo, e2_hi),
        up_lo=beam_lo, up_hi=sp["hi"] * e2_hi,
        w_mean=w, a_geom=a_geom, a_geom_pt=a_geom_pt,
        eps_pair_standin=eps_pair_standin,
        down_hi=a_geom,                       # geometry alone: an upper bound
        down_standin=a_geom * eps_pair_standin,
        eps_det_shipped=EPS_DET_CHANG)


def kappa_of_quadrupole(q_charge_fm2: float) -> float:
    """a_2(m = +-1) / |t| [GeV^-2] for a 6Li charge quadrupole [fm^2].

    `a2_from_quadrupole` wants the m = +1 point-MATTER quadrupole, which is
    twice the charge one for N = Z (design G eq. (G5)); A = 6.
    """
    return _l.a2_from_quadrupole(2.0 * q_charge_fm2, A_LI6, 1.0, 1)


def t2_mean(slope_b: float, t_max: float) -> float:
    """<|t|^2> of dN/d|t| ~ B exp(-B|t|) truncated at `t_max` [GeV^4]."""
    u = slope_b * t_max
    num = 1.0 - math.exp(-u) * (1.0 + u + 0.5 * u * u)
    return (2.0 / (slope_b * slope_b)) * num / (1.0 - math.exp(-u))


def t_fraction(slope_b: float, t_max: float) -> float:
    """Fraction of the coherent sample with |t| < `t_max`."""
    return 1.0 - math.exp(-slope_b * t_max)


def info_fraction(slope_b: float, t_max: float) -> float:
    """Fraction of the FISHER INFORMATION on kappa kept by |t| < `t_max`."""
    return t2_mean(slope_b, t_max) * t_fraction(slope_b, t_max) \
        / (2.0 / (slope_b * slope_b))


def flip_plan_pzz2(pzz: float = None) -> tuple:
    """(<P_zz^2>, [P_zz per category]) of `tensor_flip_plan`.

    The two fills are m = +-1-rich at P_zz = +pzz and m = 0-rich at -2 pzz,
    equal luminosity shares, both with the axis transverse -- so the
    luminosity-weighted <P_zz^2> is 2.5 pzz^2, which is 2.5x the P_zz^2 of a
    single +pzz fill.  P_zz = 1 - 3 p_0 is read off the populations, not
    retyped.

    THIS IS THE OPTIMAL COMBINATION, not the background-immune one; see
    `flip_plan_pzz2_difference`.
    """
    if pzz is None:
        pzz = _l.Scenario().pol_ion_tensor
    plan = _l.tensor_flip_plan(pzz)
    states = [1.0 - 3.0 * c.populations[1] for c in plan.categories]
    mean = sum(c.lumi_fraction * s * s
               for c, s in zip(plan.categories, states))
    return mean, states


def flip_plan_variances(pzz: float = None) -> dict:
    """Variance on kappa, in units of 1/(2 N <t^2>), for three estimators of
    the SAME two-fill run: the optimal combination, the background-immune
    difference, and a single +pzz fill of the same total luminosity."""
    _mean, states = flip_plan_pzz2(pzz)
    p1, p2 = states
    return dict(
        # each fill carries half the events, so 2/(N/2) = 4/N per amplitude
        difference=(2.0 / 0.5 + 2.0 / 0.5) / ((2.0 * (p1 - p2)) ** 2),
        optimal=1.0 / (2.0 * (0.5 * p1 * p1 + 0.5 * p2 * p2)),
        single=1.0 / (2.0 * p1 * p1),
        states=states)


def flip_plan_pzz2_difference(pzz: float = None) -> float:
    """The EFFECTIVE <P_zz^2> of the background-immune two-fill difference --
    the estimator section C2.5's separation argument requires.

    It is `flip_plan_pzz2()` degraded by the difference's variance penalty,
    1.1111 (5.4 % in delta): 0.90 / (10/9) = 0.81 at the shipped plan.  Quote
    THIS one whenever the photon-polarization cos 2phi has to be removed, and
    the optimal 0.90 only when it is assumed absent.
    """
    mean, _states = flip_plan_pzz2(pzz)
    v = flip_plan_variances(pzz)
    return mean / (v["difference"] / v["optimal"])


#: The two <P_zz^2> the report prints side by side.
PZZ2_OPTIMAL = flip_plan_pzz2()[0]
PZZ2_DIFFERENCE = flip_plan_pzz2_difference()


def delta_kappa(n_events: float, pzz2: float, slope_b: float,
                t_max: float) -> float:
    """1 sigma on kappa = a_2/|t| [GeV^-2] from `n_events` reconstructed."""
    return 1.0 / math.sqrt(2.0 * pzz2 * n_events * t2_mean(slope_b, t_max))


def n_for_significance(n_sigma: float, kappa: float, pzz2: float,
                       slope_b: float, t_max: float) -> float:
    """Reconstructed coherent events needed for an `n_sigma` measurement."""
    return (n_sigma / kappa) ** 2 / (2.0 * pzz2 * t2_mean(slope_b, t_max))


def _windows():
    """(label, sigma_nb, slope_b) for every (Q^2 window, density) combination,
    per vector meson.  sigma and B ALWAYS travel together: a bigger nucleus
    means both a steeper |t| slope and a smaller cross section, so every
    measured-radius row uses the measured-radius slope."""
    out = {}
    floors = {}
    for r in _l.estarlight_li6_q2_floors():
        floors[(r.vm, r.q2_floor_gev2)] = r
    for row in _l.estarlight_li6_coherent():
        combos = [("Q2>0.7", row.sigma_q7_nb, row.b_fit)]
        for q2 in Q2_FLOORS:
            f = floors[(row.vm, q2)]
            lab = _floor_label(q2)
            combos.append((lab, f.sigma_nb, f.b_fit))
            combos.append((lab + ",Rmeas", f.sigma_rmeas_nb, f.b_rmeas))
        out[row.vm] = combos
    return out


def channel_rows(lumi_fb_per_nucleon: float = None, eps_det: float = 1.0,
                 both_leptons: bool = False, pzz2: float = None):
    """One dict per (vector meson, Q^2 window, density), with the event chain
    and the reach.

    `both_leptons` switches the J/psi branching from e+e- only
    (`EstarlightLi6Row::branching`) to e+e- AND mu+mu-
    (`branching_all`); it is a no-op for phi and rho, which have no second
    reconstructible channel.  `pzz2` defaults to the BACKGROUND-IMMUNE
    `PZZ2_DIFFERENCE`, not to the optimal 0.90.
    """
    if lumi_fb_per_nucleon is None:
        lumi_fb_per_nucleon = _l.Scenario().lumi_fb_per_nucleon
    if pzz2 is None:
        pzz2 = PZZ2_DIFFERENCE
    lumi_nucleus_fb = lumi_fb_per_nucleon / A_LI6
    t_max = _l.COHERENT_T_MAX_DEFAULT
    kappa = kappa_of_quadrupole(_l.LI6_QUADRUPOLE_FM2)
    windows = _windows()
    out = []
    for row in _l.estarlight_li6_coherent():
        branching = row.branching_all if both_leptons else row.branching
        for window, sigma_nb, slope_b in windows[row.vm]:
            n_prod = sigma_nb * 1e6 * lumi_nucleus_fb   # 1 nb x 1 fb^-1 = 1e6
            n_br = n_prod * branching * t_fraction(slope_b, t_max)
            n_det = n_br * eps_det
            dk = delta_kappa(n_det, pzz2, slope_b, t_max)
            out.append(dict(vm=row.vm, window=window, sigma_nb=sigma_nb,
                            slope_b=slope_b, branching=branching,
                            pzz2=pzz2, n_produced=n_prod,
                            n_branching=n_br, n_detected=n_det,
                            delta_kappa=dk, delta_a2_at_0p3=0.3 * dk,
                            significance=kappa / dk,
                            lumi_for_3sigma=lumi_fb_per_nucleon
                            * (3.0 * dk / kappa) ** 2))
    return out


def _row(rows, vm, window):
    for r in rows:
        if r["vm"] == vm and r["window"] == window:
            return r
    raise KeyError((vm, window))


# ------------------------------------------------------------------ report

def report(out=sys.stdout) -> dict:
    """Print the whole O5 estimate and return the numbers section C2 quotes."""
    w = out.write
    sc = _l.Scenario()
    t_max = _l.COHERENT_T_MAX_DEFAULT
    b_default = _l.CoherentScenario().slope_b
    pzz2, states = flip_plan_pzz2()
    pzz2_diff = PZZ2_DIFFERENCE
    k_meas = kappa_of_quadrupole(_l.LI6_QUADRUPOLE_FM2)
    k_gfmc = kappa_of_quadrupole(_l.LI6_QUADRUPOLE_GFMC_FM2)
    q_model = 0.5 * _l.ClusterConfigSampler().q_matter_analytic_fm2(1)
    k_model = kappa_of_quadrupole(q_model)
    # The headline row, computed once here because section 3b quotes scaled
    # versions of it before section 5 prints it.
    _head = _row(channel_rows(eps_det=EPS_DET_CHANG, both_leptons=True),
                 "jpsi", "Q2>0")

    def _s_at(f):
        """S of the headline row with eps_det multiplied by `f`."""
        return _head["significance"] * math.sqrt(f)

    def _l3_at(f):
        """3-sigma luminosity [fb^-1/u] of the headline row at that factor."""
        return _head["lumi_for_3sigma"] / f

    w("O5 -- does 6Li's tensor a_2 survive EIC statistics?\n")
    w("=" * 78 + "\n\n")
    w("1. THE SIGNAL.  a_2(m=+-1, |t|) = kappa |t|, kappa = a2_from_quadrupole"
      "(2Q, 6, 1, 1)\n")
    w("   %-34s %12s %12s\n" % ("Q_charge [fm^2]", "kappa [GeV^-2]",
                                "a_2(|t|=0.3)"))
    for name, q, k in (("measured LI6_QUADRUPOLE_FM2", _l.LI6_QUADRUPOLE_FM2,
                        k_meas),
                       ("GFMC AV18+IL7", _l.LI6_QUADRUPOLE_GFMC_FM2, k_gfmc),
                       ("this alpha+d geometry", q_model, k_model)):
        w("   %-34s %12.8f %12.6f\n" % ("%s = %+.6f" % (name, q), k, 0.3 * k))
    w("\n   The |t| SLOPE IS THE ENEMY.  a_2 is linear in |t| but the coherent\n"
      "   sample is exp(-B|t|), so the information-weighted modulation is\n"
      "   kappa sqrt(<t^2>), not kappa x 0.3:\n")
    for b in (38.9, b_default, 55.0):
        w("      B = %5.1f GeV^-2 : sqrt(<t^2>) = %.6f GeV^2 -> a_2,eff = "
          "%.6f  (%.3f %%)\n"
          % (b, math.sqrt(t2_mean(b, t_max)), k_meas * math.sqrt(t2_mean(b, t_max)),
             100.0 * k_meas * math.sqrt(t2_mean(b, t_max))))
    w("   and |t| < COHERENT_T_MAX_DEFAULT = %g GeV^2 costs nothing: it keeps\n"
      "   %.4f of the rate and %.4f of the Fisher information at B = %g.\n"
      % (t_max, t_fraction(b_default, t_max),
         info_fraction(b_default, t_max), b_default))

    w("\n2. THE LUMINOSITY.  Scenario::lumi_fb_per_nucleon = %g fb^-1/u"
      " ('one EIC year',\n"
      "   band {1, 10, 100}); there is NO Li luminosity in any source the\n"
      "   repository has seen (docs/surveys/needs_survey.md section 3.8).\n"
      "   Per NUCLEUS that is %g / %d = %.6f fb^-1 of e+6Li.\n"
      % (sc.lumi_fb_per_nucleon, sc.lumi_fb_per_nucleon, A_LI6,
         sc.lumi_fb_per_nucleon / A_LI6))
    yr = _l.yr_optics("6Li", 99.5)
    tag = _l.tagging_optics("6Li", 99.5)
    w("   THE FAR-FORWARD WORKING POINT IS ASSUMED TO BE THE YELLOW REPORT\n"
      "   ONE: Optics::lumi_fraction = %g ('%s').  LiPolGen's own de-squeezed\n"
      "   6Li tagging point, '%s', delivers %.4f of the machine luminosity,\n"
      "   and estarlight_li6.md section 2e says to multiply by it.  It is NOT\n"
      "   applied here, for a stated reason: the efficiency this chain uses\n"
      "   (%.4f) is arXiv:2511.05638's IR-8 SECONDARY-FOCUS number, and that\n"
      "   paper's own argument for IR-8 is that the secondary focus buys the\n"
      "   low-pT far-forward acceptance WITHOUT the beta*_x de-squeeze that\n"
      "   costs luminosity in IR-6 (its p. 4).  The two are alternatives, not\n"
      "   multipliers.  If the tensor run must instead share the de-squeezed\n"
      "   tagging optics, every significance below multiplies by\n"
      "   sqrt(%.4f) = %.4f and every 3-sigma luminosity divides by it:\n"
      "   the headline falls from the value in section 5 to that x %.4f.\n"
      % (yr.lumi_fraction, yr.name, tag.name, tag.lumi_fraction,
         EPS_DET_CHANG, tag.lumi_fraction, math.sqrt(tag.lumi_fraction),
         math.sqrt(tag.lumi_fraction)))
    w("\n   tensor_flip_plan(%g): P_zz = %s, equal shares, theta_S = 90 deg\n"
      "   -> <P_zz^2> = %.4f (sqrt = %.6f), which is %.2fx the P_zz^2 of a\n"
      "   single +%g fill.  That is the OPTIMAL combination.  The\n"
      "   BACKGROUND-IMMUNE difference of the two fills (section 4, and the\n"
      "   estimator the separation argument requires) has an effective\n"
      "   <P_zz^2> of %.4f -- 5.4 %% worse in delta.  EVERY HEADLINE BELOW IS\n"
      "   QUOTED ON %.4f; the tables print both.\n"
      % (sc.pol_ion_tensor, states, pzz2, math.sqrt(pzz2),
         pzz2 / sc.pol_ion_tensor ** 2, sc.pol_ion_tensor, pzz2_diff,
         pzz2_diff))

    w("\n3. THE Q^2 WINDOW, WHICH IS WHERE MOST OF THE RATE WAS MISSING.\n")
    w("   estarlight_li6_coherent()'s sigma_nb stops at Q^2 > 0.1 GeV^2.  That\n"
      "   is arXiv:2511.05638's ACCEPTANCE-STUDY kinematic range copied\n"
      "   verbatim, not a physics window, and it is not eSTARlight's full Q^2\n"
      "   reach.  estarlight_li6_q2_floors() is the same runs with\n"
      "   MIN_GAMMA_Q2 lowered; 'Q2>0' means no floor at all, i.e.\n"
      "   Q^2_min = (m_e E_gamma)^2/(E_e(E_e-E_gamma)) ~ 1e-9 GeV^2.\n\n")
    w("   %-5s %-9s %12s %8s %12s %8s\n"
      % ("vm", "floor", "sigma[nb]", "B", "sigma_Rmeas", "B_Rmeas"))
    for r in _l.estarlight_li6_q2_floors():
        w("   %-5s %-9s %12.3f %8.1f %12.3f %8.1f\n"
          % (r.vm, _floor_label(r.q2_floor_gev2), r.sigma_nb, r.b_fit,
             r.sigma_rmeas_nb, r.b_rmeas))
    w("\n   Removing the 0.1 GeV^2 floor multiplies the coherent rate by\n")
    floors = {(r.vm, r.q2_floor_gev2): r for r in _l.estarlight_li6_q2_floors()}
    for vm in ("jpsi", "phi", "rho"):
        f0, f1 = floors[(vm, 0.0)], floors[(vm, 0.1)]
        f2 = floors[(vm, 0.01)]
        w("      %-5s x %.3f (no floor)   x %.3f (Q^2 > 0.01 floor)\n"
          % (vm, f0.sigma_nb / f1.sigma_nb, f2.sigma_nb / f1.sigma_nb))
    w("   and it moves the |t| slope by less than 0.5 %% (J/psi 38.9 -> 38.8),\n"
      "   so the recoil p_T spectrum the far-forward acceptance cuts on is\n"
      "   the SAME sample.  What it does move is <Q^2>: 0.906 -> 0.143 GeV^2\n"
      "   for J/psi, and <W> 32.2 -> 30.2 GeV.\n"
      "   THE EFFICIENCY IS NOT MEASURED THERE.  arXiv:2511.05638 generated\n"
      "   0.1 < Q^2 < 100 GeV^2, so no point of its %.4f was ever evaluated\n"
      "   below Q^2 = 0.1.  Applying it to the photoproduction region is an\n"
      "   EXTRAPOLATION.  Two things bound it: the efficiency is a recoil-\n"
      "   nucleus tagging acceptance and does NOT require the scattered\n"
      "   electron (which at Q^2 < 0.1 goes down the beam pipe), and its\n"
      "   published W dependence (their Fig. 2) RISES as W falls, which is\n"
      "   the direction the low-Q^2 sample moves.  Neither is a measurement.\n"
      % EPS_DET_CHANG)

    # ------------------------------------------- 3b: the efficiency chain
    ch = eps_det_chain()
    sp = ch["species"]
    w("\n3b. THE EFFICIENCY CHAIN ITSELF, WHICH HAD TWO DEFECTS PULLING\n"
      "    OPPOSITE WAYS.  Until 2026-09-04 the chain was\n"
      "    sigma x BR(l+l-) x %.4f and nothing else, with the %.4f treated as\n"
      "    flat in |t|, in phi_Delta, in Q^2, across mesons AND -- unstated --\n"
      "    across beam energy, and with NO factor at all for the decay\n"
      "    leptons.  Both omissions are now priced.\n" % (EPS_DET_CHANG,
                                                          EPS_DET_CHANG))

    w("\n    (a) THE BEAM ENERGY.  UP.  %.4f is arXiv:2511.05638's 7Li number\n"
      "    at 18 x %.1f GeV/u -- 7Li's own TOP energy -- and this tree's\n"
      "    configuration-identical eSTARlight row for those beams has\n"
      "    <W> = %.1f GeV.  It is applied to a 6Li sample at %g x %.1f with\n"
      "    <W> = %.1f.  The same paper MEASURES that dependence, on 3He\n"
      "    (its sec. V.B, in code as chang26_he3_energy_scan()):\n"
      % (EPS_DET_CHANG, _l.COHERENT_JPSI_EFF_IR8_LI7_E_ION_GEV,
         _l.COHERENT_JPSI_EFF_IR8_LI7_W_MEAN_GEV,
         BEAM_LI6.electron_energy, BEAM_LI6.ion_momentum_per_nucleon,
         ch["w_mean"]))
    for r in _l.chang26_he3_energy_scan():
        w("        e3He %4.1f x %5.1f GeV/u : %6.2f %%\n"
          % (r.e_electron_gev, r.e_ion_gev, 100.0 * r.efficiency))
    w("      -> d ln(eff)/d ln(E_ion) = %.4f and %.4f.  Efficiency RISES as\n"
      "         the beam energy falls, so applying a top-energy number at\n"
      "         %.1f GeV/u is CONSERVATIVE.  For the step actually required,\n"
      "         %.1f -> %.1f GeV/u, that is a factor %.4f-%.4f UP:\n"
      "         eps_det = %.4f-%.4f, not %.4f.\n"
      % (he3_energy_slopes()[0], he3_energy_slopes()[1],
         BEAM_LI6.ion_momentum_per_nucleon,
         _l.COHERENT_JPSI_EFF_IR8_LI7_E_ION_GEV,
         BEAM_LI6.ion_momentum_per_nucleon, ch["beam_lo"], ch["beam_hi"],
         EPS_DET_CHANG * ch["beam_lo"], EPS_DET_CHANG * ch["beam_hi"],
         EPS_DET_CHANG))
    w("         IT IS NOT THE PAPER'S OWN 1.687.  That is their 183 -> 100\n"
      "         GeV/u ratio (%.4f/%.4f), a step %.2fx larger in ln E than the\n"
      "         one this transfer needs.  Carrying it whole would give\n"
      "         S = %.2f sigma and 3 sigma at %.1f fb^-1/u, which is scaling\n"
      "         by the wrong lever arm, not a more honest number.\n"
      % (_l.chang26_he3_energy_scan()[1].efficiency,
         _l.chang26_he3_energy_scan()[0].efficiency,
         math.log(_l.chang26_he3_energy_scan()[1].e_ion_gev
                  / _l.chang26_he3_energy_scan()[0].e_ion_gev)
         / math.log(BEAM_LI6.ion_momentum_per_nucleon
                    / _l.COHERENT_JPSI_EFF_IR8_LI7_E_ION_GEV),
         _s_at(1.687), _l3_at(1.687)))
    w("         The W axis says the same thing and is NOT independent of it:\n"
      "         <W> %.1f -> %.1f, and by their Fig. 2 lower W means higher\n"
      "         efficiency.  Not applied on top -- it is the same effect.\n"
      % (_l.COHERENT_JPSI_EFF_IR8_LI7_W_MEAN_GEV, ch["w_mean"]))

    conf = ch["species_confounded"]
    nl = species_list_is_not_a_species_lever()
    e4_eff = {r.nucleus: r.efficiency
              for r in _l.chang26_species_efficiency()}["4He"]
    w("\n    (b) THE SPECIES.  DIRECTION UNDETERMINED, AND THIS LEG IS THE\n"
      "    ONE UNQUANTIFIED FACTOR IN THE CHAIN.  Until 2026-09-04 this\n"
      "    section, and four other sites, read: 'Every entry of\n"
      "    chang26_species_efficiency() sits at A/Z x E = 275 GeV/e, so its\n"
      "    A-ordering is a species lever on its own.'  THAT IS A NON SEQUITUR\n"
      "    AND IT IS RETRACTED.  Eliminating rigidity does not leave species\n"
      "    alone: at fixed R = A E/Z both E/u = R Z/A and p_z = Z R still\n"
      "    vary down the list.\n")
    w("        E/u   %.0f (7Li) to %.0f (3He) GeV/u  -- a x%.2f spread\n"
      "        Z     %d to %d\n"
      "        p_z   %.0f to %.0f GeV\n"
      % (nl["e_per_u"][0], nl["e_per_u"][1], nl["e_per_u_spread"],
         nl["z"][0], nl["z"][1], nl["p_z"][0], nl["p_z"][1]))
    w("      -> and the confound is the SIZE of the effect being read off it:\n"
      "         at (a)'s own d ln(eff)/d ln(E) the list's x%.2f spread in E/u\n"
      "         is worth x%.4f-%.4f in the efficiency, as large as the whole\n"
      "         x%.4f-%.4f 'species gain' the confounded reading returns.\n"
      "         The list is a JOINT (A, Z, E/u) lever and nothing in it\n"
      "         decomposes that; this chain asserted the decomposition and\n"
      "         never tested it.\n"
      % (nl["e_per_u_spread"], nl["eff_worth_of_e_spread"][0],
         nl["eff_worth_of_e_spread"][1],
         conf["lo"], conf["hi"]))
    w("\n        THE CONFOUNDED READING, kept only so the confound is visible:\n"
      "        4He (%.4f, at 137 GeV/u) interpolated onto 7Li (%.4f, at 118)\n"
      % (e4_eff, EPS_DET_CHANG))
    for k in ("log_r", "log_a", "pt_threshold", "lin_a"):
        w("        %-18s eps(6Li) = %.4f   (x%.4f)\n"
          % (k, conf["eps_li6"][k], conf["factors"][k]))
    w("        -> x%.4f-%.4f.  Its two anchors are 19 GeV/u apart in E/u,\n"
      "           which on (a)'s own slopes is worth x%.4f-%.4f by itself.\n"
      % (conf["lo"], conf["hi"],
         *sorted(math.exp(sl * math.log(118.0 / 137.0))
                 for sl in he3_energy_slopes())))
    w("\n        THE READING WITH NO ENERGY STEP IN IT, and the one the band\n"
      "        now uses.  FOUR of the seven entries share a beam energy --\n"
      "        %s at %.0f GeV/u -- and 6Li's own\n"
      "        fixed-rigidity energy is Z/A x 275 = 137.5, so 4He -> 12C\n"
      "        brackets A = 6 with NO energy step at all:\n"
      % (", ".join(sp["nuclei"]), sp["at_energy_gev"]))
    for k in sorted(sp["eps_li6"], key=lambda k: sp["eps_li6"][k]):
        w("        %-18s eps(6Li, 137) = %.4f   (x%.4f)\n"
          % (k, sp["eps_li6"][k], sp["factors"][k]))
    w("        -> x%.4f-%.4f.  IT STRADDLES 1.  On this reading the\n"
      "           DIRECTION of the 7Li -> 6Li substitution is not\n"
      "           established, let alone its size, and the band's TOP is the\n"
      "           span it leaves rather than a point (section 7).\n"
      % (sp["lo"], sp["hi"]))
    w("\n        THE ONE PLACE THE pT-THRESHOLD FORM CAN BE TESTED, read three\n"
      "        ways.  eps = exp(-B pT_cut^2) with pT_cut = %.4f GeV inverted\n"
      "        from the 7Li row; the criterion arXiv:2511.05638 states is\n"
      "        'within a safe distance FROM THE BEAM', i.e. a cut on the\n"
      "        ANGLE, so pT_cut = theta p_z = theta Z R -- IT CARRIES Z.\n"
      % conf["pt_cut_gev"])
    w("        (1) pT_cut ~ Z.  3He and 4He are Z = 2, 7Li is Z = 3, so the\n"
      "            pair's pT_cut^2 is (2/3)^2 of the inverted one: predicts\n"
      "            %.4f against the measured %.4f, i.e. %.1f %%.\n"
      % (conf["he3_he4_predicted_z"], conf["he3_he4_measured"],
         100.0 * abs(conf["he3_he4_ratio_z"] - 1.0)))
    w("        (2) THE SHIPPED READING applied 7Li's Z = 3 cut unchanged to\n"
      "            the Z = 2 pair -- its own comment said 'same Z' without\n"
      "            noticing it was the WRONG Z -- and got %.4f, 'overstating\n"
      "            by %.0f %%'.  That, and only that, is where the\n"
      "            instruction 'read the low end' came from.  WITHDRAWN.\n"
      % (conf["he3_he4_predicted_wrong_z"],
         100.0 * (conf["he3_he4_overstatement_wrong_z"] - 1.0)))
    w("        (3) pT_cut Z-INDEPENDENT.  Over the LIGHT half of the list\n"
      "            this fits the ABSOLUTE values far better than (1) does --\n"
      "            pred/meas, every entry:\n")
    for n in ("2D", "3He", "4He", "9Be", "12C", "16O"):
        w("              %-4s  %.4f   (Z-scaled: %.4f)\n"
          % (n, conf["absolute_pred_over_meas"][n],
             conf["absolute_pred_over_meas_z"][n]))
    w("            -- 2D, 4He and 9Be to %.1f %% or better, and 3He off by\n"
      "            %.0f %%.  3He is the one nucleus at E/u = 183 rather than\n"
      "            the family's 137, for which (a)'s own scan predicts an\n"
      "            %.0f-%.0f %% shortfall, so the ONE misfit among the light\n"
      "            entries is explained by the SAME confound as (b).\n"
      "            SAY WHAT THIS DOES NOT DO: the form degrades with A and\n"
      "            misses 12C by %.0f %% and 16O by a factor %.1f, so it is\n"
      "            not a description of the list -- it is a local\n"
      "            interpolation, and 6Li sits inside the range where it\n"
      "            works.  That is the most that can be claimed for it.\n"
      % (100.0 * max(abs(conf["absolute_pred_over_meas"][n] - 1.0)
                     for n in ("2D", "4He", "9Be")),
         100.0 * (conf["absolute_pred_over_meas"]["3He"] - 1.0),
         100.0 * (1.0 - conf["he3_energy_shortfall"][1]),
         100.0 * (1.0 - conf["he3_energy_shortfall"][0]),
         100.0 * (conf["absolute_pred_over_meas"]["12C"] - 1.0),
         conf["absolute_pred_over_meas"]["16O"]))
    w("        Reading (1) says the form is right; (3) says it is right\n"
      "        except where the energy confound bites; (2) says it\n"
      "        overstates.  NOTHING IN THE TABLE PICKS BETWEEN THEM, so\n"
      "        nothing licenses a directional instruction either way.\n")
    w("\n        Chained with 6Li's own %.1f -> %.1f GeV/u energy step\n"
      "        (x%.4f-%.4f) the two legs give x%.4f-%.4f in total.\n"
      % (BEAM_LI6_TOP.ion_momentum_per_nucleon,
         BEAM_LI6.ion_momentum_per_nucleon,
         ch["energy_from_li6_top"][0], ch["energy_from_li6_top"][1],
         ch["two_leg_lo"], ch["two_leg_hi"]))

    w("\n    (c) THE DECAY LEPTONS.  DOWN, AND THIS FACTOR WAS NOT IN THE\n"
      "    CHAIN AT ALL.  arXiv:2511.05638's number is the fraction of\n"
      "    scattered NUCLEI in the far-forward acceptance -- \"32.23 % of the\n"
      "    scattered 3He nuclei occur within a safe distance from the beam\"\n"
      "    -- and its own p. 4 says the simulation \"only accounts for the\n"
      "    acceptance effect and does not incorporate the efficiencies of the\n"
      "    detector.  Additionally, we did not account for the efficiency and\n"
      "    acceptance of the reconstructed distribution\".  So there is no\n"
      "    central-detector acceptance and no reconstruction efficiency for\n"
      "    the e+e-/mu+mu- pair anywhere above.  It is <= 1 by construction.\n")
    w("\n    What CAN be bounded here is the GEOMETRY, and it is not the\n"
      "    problem.  At %g x %.1f the J/psi is nearly at rest in the lab, so\n"
      "    its 1.548 GeV decay leptons are central (Scenario::eta_max = %g,\n"
      "    this tree's own crude central-detector edge -- documented for the\n"
      "    SCATTERED ELECTRON, because no decay-lepton acceptance exists\n"
      "    here at all):\n"
      % (BEAM_LI6.electron_energy, BEAM_LI6.ion_momentum_per_nucleon,
         _l.Scenario().eta_max))
    w("        %8s %8s %10s %12s %12s\n"
      % ("W [GeV]", "E_gam", "y(J/psi)", "A(eta only)", "A(+pT>0.2)"))
    # The largest W eSTARlight can reach at these beams: E_gamma <= E_e.
    _pn = math.sqrt(BEAM_LI6.ion_momentum_per_nucleon ** 2 - M_PROTON_GEV ** 2)
    w_top = math.sqrt(M_PROTON_GEV ** 2
                      + 2.0 * BEAM_LI6.electron_energy
                      * (BEAM_LI6.ion_momentum_per_nucleon + _pn))
    for w_gev in (20.0, jpsi_w_mean(0.0), jpsi_w_mean(0.1), 40.0, 50.0, w_top):
        y_v, e_g = jpsi_lab_rapidity(w_gev)
        w("        %8.1f %8.3f %+10.3f %12.4f %12.4f\n"
          % (w_gev, e_g, y_v, lepton_pair_acceptance(w_gev),
             lepton_pair_acceptance(w_gev, pt_min=0.2)))
    w("      -> A_geom = %.4f at this sample's own <W> = %.1f, and never\n"
      "         below %.2f anywhere in the accessible W range (E_gamma <= E_e\n"
      "         caps W at %.1f GeV).  SCHC transverse decay throughout, the\n"
      "         pessimistic weighting; isotropic gives %.4f at <W>.\n"
      "         p_T(J/psi) ~ 0.16 GeV neglected.\n"
      % (ch["a_geom"], ch["w_mean"], lepton_pair_acceptance(w_top), w_top,
         lepton_pair_acceptance(ch["w_mean"], schc=False)))
    w("         THE RECONSTRUCTION EFFICIENCY IS UNBOUNDED HERE.  No\n"
      "         per-lepton tracking or PID efficiency exists in this tree.\n"
      "         The sibling ../PolarizedLithiumSim assumes eff_track = %.2f,\n"
      "         which IT labels a stand-in; the pair would then cost %.4f,\n"
      "         and A_geom x that is %.4f.  That number is quoted below as a\n"
      "         STAND-IN and is not a measurement of anything.\n"
      % (EPS_TRACK_SIBLING_STANDIN, ch["eps_pair_standin"],
         ch["down_standin"]))

    w("\n    (d) THEY PARTLY CANCEL, WHICH IS WHY BOTH HAD TO BE STATED.\n"
      "    The energy leg alone times the stand-in pair factor is\n"
      "    %.4f x %.4f = %.4f -- unity to %.1f %%.  A reader told only about\n"
      "    the beam energy moves the verdict up; a reader told only about the\n"
      "    leptons moves it down; the pair of them barely move it at all.\n"
      "    With the species leg on top the product is %.4f to %.4f -- a SPAN,\n"
      "    because (b) is a span.\n"
      % (ch["beam_lo"], ch["down_standin"], ch["beam_lo"] * ch["down_standin"],
         100.0 * abs(ch["beam_lo"] * ch["down_standin"] - 1.0),
         ch["two_leg_top_min"] * ch["down_standin"],
         ch["two_leg_top_max"] * ch["down_standin"]))

    w("\n4. THE RATE AND THE REACH, at %g fb^-1/u, <P_zz^2> = %.4f"
      " (background-immune).\n" % (sc.lumi_fb_per_nucleon, pzz2_diff))
    head = ("   %-5s %-14s %10s %10s %10s %10s %10s %8s %10s\n"
            % ("vm", "window", "sigma[nb]", "produced", "x BR", "x eff",
               "d a_2(0.3)", "S", "3sig[fb/u]"))
    for eps, both, label in (
            (EPS_DET_CHANG, False,
             "efficiency %.4f (COHERENT_JPSI_EFF_IR8_LI7, a 7Li stand-in), "
             "J/psi -> e+e- only" % EPS_DET_CHANG),
            (EPS_DET_CHANG, True,
             "efficiency %.4f, J/psi -> e+e- AND mu+mu- (branching_all)"
             % EPS_DET_CHANG),
            (1.0, True,
             "efficiency 1.0 (perfect detection -- an upper bound, not a "
             "projection), both leptons")):
        w("\n   %s\n" % label)
        w(head)
        for r in channel_rows(eps_det=eps, both_leptons=both):
            w("   %-5s %-14s %10.4g %10.4g %10.4g %10.4g %10.5f %8.3f %10.2f\n"
              % (r["vm"], r["window"], r["sigma_nb"], r["n_produced"],
                 r["n_branching"], r["n_detected"], r["delta_a2_at_0p3"],
                 r["significance"], r["lumi_for_3sigma"]))

    w("\n   Events needed, at <P_zz^2> = %.4f:\n" % pzz2_diff)
    for b in (38.9, b_default, 55.0):
        w("      B = %5.1f : 3 sigma at N = %.4g, 5 sigma at N = %.4g\n"
          % (b, n_for_significance(3.0, k_meas, pzz2_diff, b, t_max),
             n_for_significance(5.0, k_meas, pzz2_diff, b, t_max)))

    # ---------------------------------------------------------- the verdict
    diff = channel_rows(eps_det=EPS_DET_CHANG, both_leptons=True)
    opt = channel_rows(eps_det=EPS_DET_CHANG, both_leptons=True,
                       pzz2=pzz2)
    ee = channel_rows(eps_det=EPS_DET_CHANG, both_leptons=False)
    head_row = _row(diff, "jpsi", "Q2>0")
    head_opt = _row(opt, "jpsi", "Q2>0")
    old_row = _row(ee, "jpsi", "Q2>0.1")
    pess = _row(diff, "jpsi", "Q2>0.01,Rmeas")
    cons = _row(diff, "jpsi", "Q2>0.01")
    rmeas = _row(diff, "jpsi", "Q2>0,Rmeas")

    w("\n5. THE DECIDING NUMBER.  Coherent J/psi -- the channel the case is\n"
      "   built on, the only one inside COHERENT_MX_MIN_DEFAULT = %g GeV and\n"
      "   the only one [Mant24] actually computed -- over the WHOLE Q^2 range\n"
      "   with both lepton channels reconstructed gives N = %.4g events at\n"
      "   %g fb^-1/u, hence\n"
      "      delta a_2(|t| = 0.3) = %.5f  against a predicted a_2 = %.6f\n"
      "      S = %.3f sigma;  3 sigma needs %.1f fb^-1/u -- INSIDE the\n"
      "      {1, 10, 100} fb^-1/u band, about ONE EIC year.\n"
      % (_l.COHERENT_MX_MIN_DEFAULT, head_row["n_detected"],
         sc.lumi_fb_per_nucleon, head_row["delta_a2_at_0p3"], 0.3 * k_meas,
         head_row["significance"], head_row["lumi_for_3sigma"]))
    w("   On the OPTIMAL <P_zz^2> = %.2f instead: S = %.3f sigma, 3 sigma at\n"
      "   %.1f fb^-1/u.\n"
      % (pzz2, head_opt["significance"], head_opt["lumi_for_3sigma"]))
    w("\n   THE LADDER, so the reader can see what each correction is worth:\n")
    for label, r in (("Q^2 > 0.1, e+e- only (the number this item shipped "
                      "with, on <P_zz^2> = 0.90)", None),
                     ("Q^2 > 0.1, e+e- only", old_row),
                     ("Q^2 > 0.1, both leptons", _row(diff, "jpsi", "Q2>0.1")),
                     ("Q^2 > 0.01, both leptons", cons),
                     ("Q^2 > 0.01, both leptons, R = 2.589 fm", pess),
                     ("no Q^2 floor, both leptons, R = 2.589 fm", rmeas),
                     ("no Q^2 floor, both leptons", head_row)):
        if r is None:
            r0 = _row(channel_rows(eps_det=EPS_DET_CHANG, pzz2=pzz2),
                      "jpsi", "Q2>0.1")
            w("      %-52s S = %.3f  3 sigma at %7.1f fb^-1/u\n"
              % (label, r0["significance"], r0["lumi_for_3sigma"]))
            continue
        w("      %-52s S = %.3f  3 sigma at %7.1f fb^-1/u\n"
          % (label, r["significance"], r["lumi_for_3sigma"]))
    w("      %-52s S = %.3f  3 sigma at %7.1f fb^-1/u\n"
      % ("the de-squeezed tagging optics, on the last row",
         head_row["significance"] * math.sqrt(tag.lumi_fraction),
         head_row["lumi_for_3sigma"] / tag.lumi_fraction))
    w("   -- so EVERY conservatism on this ladder EXCEPT ONE leaves 3 sigma\n"
      "      inside the {1, 10, 100} band on its own; it takes two of them\n"
      "      stacked to put it back out.  THE ONE EXCEPTION IS THE LAST ROW:\n"
      "      the de-squeezed far-forward optics of section 2 multiply every S\n"
      "      above by %.4f and every 3-sigma luminosity by %.2f, and on their\n"
      "      own that is enough -- S = %.3f at %.0f fb^-1/u, outside the band.\n"
      "      (This ladder used to end 'every single pessimism on its own\n"
      "      leaves 3 sigma inside the band' three lines after saying the\n"
      "      optics row does not.  It was a self-contradiction; the optics\n"
      "      row is now IN the ladder rather than a footnote to it.)\n"
      % (math.sqrt(tag.lumi_fraction), 1.0 / tag.lumi_fraction,
         head_row["significance"] * math.sqrt(tag.lumi_fraction),
         head_row["lumi_for_3sigma"] / tag.lumi_fraction))

    w("\n   AND THE SAME LADDER ON THE CORRECTED eps_det OF SECTION 3b, which\n"
      "   is the honest form: a BAND, not a point.  Every rung is the\n"
      "   headline row (no Q^2 floor, both leptons, default density) with\n"
      "   eps_det = %.4f multiplied by the stated factor:\n" % EPS_DET_CHANG)
    for label, f in (
            ("as shipped -- no beam-energy leg, no lepton factor", 1.0),
            ("+ beam-energy leg only (x%.4f-%.4f)"
             % (ch["beam_lo"], ch["beam_hi"]), ch["beam_lo"]),
            ("+ decay-lepton GEOMETRY only (x%.4f)" % ch["a_geom"],
             ch["a_geom"]),
            ("BAND LOW  = beam lo x geometry x pair stand-in",
             ch["beam_lo"] * ch["down_standin"]),
            ("BAND TOP, species read LOW  (x%.4f)" % ch["species_lo"],
             ch["two_leg_top_min"] * ch["down_standin"]),
            ("BAND TOP, species read HIGH (x%.4f)" % ch["species_hi"],
             ch["two_leg_top_max"] * ch["down_standin"]),
            ("... and the de-squeezed optics on the BAND LOW rung",
             ch["beam_lo"] * ch["down_standin"] * tag.lumi_fraction),
            ("... and the de-squeezed optics on the TOP rung, species lo",
             ch["two_leg_top_min"] * ch["down_standin"] * tag.lumi_fraction),
            ("... and the de-squeezed optics on the TOP rung, species hi",
             ch["two_leg_top_max"] * ch["down_standin"] * tag.lumi_fraction)):
        w("      %-52s S = %.3f  3 sigma at %7.1f fb^-1/u\n"
          % (label, _s_at(f), _l3_at(f)))
    w("      The band is open BELOW its low rung, because the per-lepton\n"
      "      reconstruction efficiency is unbounded in this tree.  What it\n"
      "      takes to leave the band: a pair acceptance x efficiency of\n"
      "      %.3f puts S at 2 sigma, and %.3f puts 3 sigma outside\n"
      "      {1, 10, 100} -- the second is a detector that reconstructs one\n"
      "      J/psi in eight, so 'inside the band' is robust to it in a way\n"
      "      'MARGINAL rather than NO' is not.\n"
      % ((2.0 / _head["significance"]) ** 2 / ch["beam_lo"],
         (_head["lumi_for_3sigma"] / 100.0) / ch["beam_lo"]))

    w("\n   Same J/psi sample (no Q^2 floor, both leptons), same everything,\n"
      "   at the OTHER two quadrupoles:\n")
    for name, k in (("GFMC AV18+IL7 -0.20", k_gfmc),
                    ("this alpha+d geometry %+.4f" % q_model, k_model)):
        w("      %-32s S = %.3f sigma (3 sigma at %.1f fb^-1/u)\n"
          % (name, k / head_row["delta_kappa"],
             sc.lumi_fb_per_nucleon * (3.0 * head_row["delta_kappa"] / k) ** 2))
    w("   -- i.e. the factor 7.5 by which this geometry overshoots Q(6Li) is\n"
      "   still the difference between a comfortable measurement and a\n"
      "   marginal one.  6Li is a near-null test; near-null is expensive.\n")

    w("\n6. THE SEPARATION FROM THE PHOTON-POLARIZATION cos 2phi.\n")
    v = flip_plan_variances()
    p1, p2 = v["states"]
    v_diff, v_opt, v_single = v["difference"], v["optimal"], v["single"]
    w("   The tensor term is ODD in P_zz (c_2 = 2 kappa |t| P_zz); the\n"
      "   photon-polarization term is EVEN in it (it is a property of the\n"
      "   photon and the unpolarized target).  The difference of the two\n"
      "   fills' fitted amplitudes, divided by P_zz(1) - P_zz(2) = %g, is\n"
      "   therefore background-free at first order.  Its cost:\n"
      "      vs the optimal use of the SAME two fills : %.4f in variance,"
      " %.4f in delta\n"
      "      vs putting all the luminosity in one +%g fill: %.4f in variance,"
      " %.4f in delta\n"
      % (p1 - p2, v_diff / v_opt, math.sqrt(v_diff / v_opt),
         p1, v_diff / v_single, math.sqrt(v_diff / v_single)))
    w("   -- the flip is not a cost, it is a %.2fx GAIN, because the m = 0-rich\n"
      "   fill carries |P_zz| = %g against the m = +-1-rich fill's %g.\n"
      % (1.0 / math.sqrt(v_diff / v_single), abs(p2), p1))
    w("   NOTE THAT THIS SURVIVES Q^2 -> 0.  Below Q^2 = 0.1 the scattered\n"
      "   electron is not detected and phi_gamma is unknown event by event,\n"
      "   but neither handle needs it: the spin-axis-referenced histogram\n"
      "   averages the background to zero under a phi_gamma-uniform\n"
      "   acceptance, and the P_zz parity does not reference the lepton plane\n"
      "   at all.  What is lost at low Q^2 is the ability to project the\n"
      "   background out DIRECTLY, which was never the control being used.\n")

    o = _l.ClusterConfigOptions()
    o.quadrupole_target_fm2 = _l.LI6_QUADRUPOLE_FM2
    eps_b0 = abs(_l.ClusterConfigSampler(o).eps_b0_equivalent())
    resid_edge = (eps_b0 * b_default * t_max) ** 2 / 4.0
    resid_mean = eps_b0 * eps_b0 * b_default ** 2 * t2_mean(b_default, t_max) / 4.0
    w("\n   THE RESIDUAL.  What breaks the parity is that the two fills do not\n"
      "   have identical |t| spectra: the deformation that makes a_2 also\n"
      "   modulates the slope, B -> B(1 + eps cos 2phi) with\n"
      "   eps = |eps_b0_equivalent()| = %.6f at the measured Q.  The\n"
      "   phi-AVERAGED spectra differ only at SECOND order, <exp(-eps B t\n"
      "   cos2phi)>_phi = I_0(eps B t) = 1 + (eps B t)^2/4 + ..., i.e. by\n"
      "      %.3e at |t| = t_max, %.3e averaged over the window.\n"
      % (eps_b0, resid_edge, resid_mean))
    w("   A background of amplitude A_gamma leaks at most that fraction of\n"
      "   itself.  Against the per-fill statistical error sqrt(2/N) that\n"
      "   matters only above (no Q^2 floor, both leptons, eff = %.4f):\n"
      % EPS_DET_CHANG)
    for r in diff:
        if r["window"] != "Q2>0":
            continue
        stat = math.sqrt(2.0 / r["n_detected"])
        w("      %-5s N = %10.4g  sqrt(2/N) = %.3e  -> A_gamma > %.3f\n"
          % (r["vm"], r["n_detected"], stat, stat / resid_edge))
    w("   so the P_zz flip is exact enough for J/psi and phi by a wide\n"
      "   margin, and for rho -- the one channel with the statistics -- it\n"
      "   stops being exact at a background amplitude of order 0.02.  NO\n"
      "   MAGNITUDE FOR A_gamma EXISTS IN THIS TREE (physics_literature.md\n"
      "   records the mechanism, STAR arXiv:2204.01625, not a number), which\n"
      "   is why the whole argument above is deliberately magnitude-free: the\n"
      "   parity in P_zz does the work.  Binning in |t| removes the residual\n"
      "   by construction; that is the control, and for rho it is mandatory.\n")
    w("\n   NOT APPLIED, for want of a number: an azimuthal smearing dilution\n"
      "   D = exp(-2 sigma_phi^2) with sigma_phi ~ sigma_pT / p_T.  No in-tree\n"
      "   coherent-VM p_T resolution exists.  At |t| ~ 1/B the transverse\n"
      "   momentum is only %.3f GeV, so this is the assumption most likely to\n"
      "   cost something.  The finite-bin dilution of `cos2phi_fit_err` is\n"
      "   %.4f at 24 bins and %.4f at 8, and is already small.\n"
      % (math.sqrt(1.0 / b_default),
         _l.estimators.cos2phi_fit_err(1.0, 1.0, 24) / math.sqrt(2.0),
         _l.estimators.cos2phi_fit_err(1.0, 1.0, 8) / math.sqrt(2.0)))

    band_lo_f = ch["beam_lo"] * ch["down_standin"]
    band_hi_f_min = ch["two_leg_top_min"] * ch["down_standin"]
    band_hi_f_max = ch["two_leg_top_max"] * ch["down_standin"]
    w("\n7. THE VERDICT.  MARGINAL, A BAND RATHER THAN A POINT -- AND THE\n"
      "   BAND'S TOP IS ITSELF A SPAN, NOT AN EDGE.\n\n"
      "   ON THE REPOSITORY'S OWN DEFAULT DENSITY, one EIC year, the whole\n"
      "   Q^2 range and both lepton channels, coherent-J/psi tensor a_2 is\n\n"
      "      LOW EDGE  S = %.2f sigma, 3 sigma at %.1f fb^-1/u\n"
      "      TOP       S = %.2f to %.2f sigma, 3 sigma at %.1f to %.1f fb^-1/u\n"
      "      -- so 3 sigma sits at %.1f to %.1f fb^-1/u across every form,\n"
      "      INSIDE the {1, 10, 100} fb^-1/u band at BOTH ends and close to\n"
      "      one EIC year.  That much is established.\n\n"
      "      WHETHER THE BAND'S TOP CROSSES 3 SIGMA IS NOT.  Until\n"
      "      2026-09-04 this section printed 'S = 2.63 to 3.15 sigma' and\n"
      "      'the band's TOP just crosses 3 sigma'.  The 3.15 came from the\n"
      "      species leg, and the species leg was read off a list whose\n"
      "      A-ordering is confounded with E/u and Z (section 3b(b)).  Read\n"
      "      off the four entries that share a beam energy the leg straddles\n"
      "      1, so the top straddles 3 sigma: %.2f to %.2f.  The BOTTOM does\n"
      "      not reach 3 sigma on any reading, and the bottom is open\n"
      "      (below).  MARGINAL is the label for a band that straddles the\n"
      "      threshold, not for a point that happens to sit under it.\n\n"
      "   The band is section 3b's, and it is a band because the single\n"
      "   largest multiplier in the chain -- eps_det -- carries three things\n"
      "   that are UNQUANTIFIED rather than uncertain and pull two ways:\n"
      "     UP, and measured: eps_det was taken at 18 x %.1f GeV/u and is\n"
      "        used at %g x %.1f.  arXiv:2511.05638's own 3He energy scan\n"
      "        makes that x%.4f-%.4f.\n"
      "     DIRECTION UNDETERMINED, and the one factor that is presented as\n"
      "        quantified without being quantified anywhere: the 7Li -> 6Li\n"
      "        substitution.  x%.4f-%.4f on the four same-energy entries --\n"
      "        IT STRADDLES 1 -- against the x%.4f-%.4f the confounded\n"
      "        reading gave and this file printed until 2026-09-04.  It is\n"
      "        the whole of the difference between a top that crosses 3\n"
      "        sigma and one that does not.\n"
      "     DOWN, and NOT bounded: there is no decay-lepton acceptance or\n"
      "        reconstruction efficiency in the chain at all.  The GEOMETRY\n"
      "        is %.4f (computed, section 3b(c)); the per-lepton\n"
      "        reconstruction efficiency does not exist anywhere in this\n"
      "        tree, and the %.4f used at the band's low end is the sibling\n"
      "        ../PolarizedLithiumSim's own stand-in, not a measurement.\n"
      "     UNMEASURED EITHER WAY: none of eps_det was evaluated below\n"
      "        Q^2 = 0.1, where 85 %% of the rate sits (section 3).\n"
      "   THE FIRST AND THE LAST NEARLY CANCEL (%.4f x %.4f = %.4f).\n"
      "   That is why both had to be stated: a reader given only one of them\n"
      "   moves the verdict in one direction, and the pair of them barely\n"
      "   move it.\n\n"
      % (_s_at(band_lo_f), _l3_at(band_lo_f),
         _s_at(band_hi_f_min), _s_at(band_hi_f_max),
         _l3_at(band_hi_f_max), _l3_at(band_hi_f_min),
         _l3_at(band_hi_f_max), _l3_at(band_lo_f),
         _s_at(band_hi_f_min), _s_at(band_hi_f_max),
         _l.COHERENT_JPSI_EFF_IR8_LI7_E_ION_GEV,
         BEAM_LI6.electron_energy, BEAM_LI6.ion_momentum_per_nucleon,
         ch["beam_lo"], ch["beam_hi"],
         ch["species_lo"], ch["species_hi"],
         ch["species_confounded_lo"], ch["species_confounded_hi"],
         ch["a_geom"], ch["eps_pair_standin"],
         ch["beam_lo"], ch["down_standin"], band_lo_f))
    w("   THE CLAIM THIS ITEM SHIPPED WITH -- 'the measurement does not exist\n"
      "   at any luminosity the EIC is quoted at' -- IS STILL REFUTED.  It\n"
      "   was computed from the 0.1 < Q^2 < 100 window with one lepton\n"
      "   channel; that row is still %.2f sigma / %.0f fb^-1/u and still\n"
      "   correct FOR THAT ROW, but the row is not the measurement.  Nothing\n"
      "   in section 3b changes that: the band's LOW end is %.2f sigma.\n\n"
      % (_row(channel_rows(eps_det=EPS_DET_CHANG, pzz2=pzz2),
              "jpsi", "Q2>0.1")["significance"],
         _row(channel_rows(eps_det=EPS_DET_CHANG, pzz2=pzz2),
              "jpsi", "Q2>0.1")["lumi_for_3sigma"], _s_at(band_lo_f)))
    w("   WHAT IS STILL NOT ESTABLISHED, and it is not a detail:\n"
      "     * THE ONE CORRECTION THAT ON ITS OWN RESTORES A NO IS\n"
      "       `Optics::lumi_fraction`.  The far-forward working point is\n"
      "       unchosen; at LiPolGen's own de-squeezed 6Li tagging optics\n"
      "       (%.4f) the band becomes %.2f-%.2f sigma with 3 sigma at\n"
      "       %.0f-%.0f fb^-1/u -- OUTSIDE the band at both ends.  Every\n"
      "       other conservatism priced here leaves 3 sigma inside on its\n"
      "       own; this one does not, and stacking the Q^2 > 0.01 floor with\n"
      "       the measured radius (%.2f sigma, %.0f fb^-1/u) is the only\n"
      "       other way out.\n"
      "     * the decay-lepton reconstruction efficiency, which is UNBOUNDED\n"
      "       here and is what makes the band open below.  It would have to\n"
      "       fall to %.2f to cost the 'above 2 sigma' half of MARGINAL and\n"
      "       to %.2f to put 3 sigma outside the band.\n"
      "     * the efficiency below Q^2 = 0.1 -- where 85 %% of the J/psi rate\n"
      "       now sits -- is UNMEASURED anywhere in this tree (section 3);\n"
      "     * a_2 itself is a closed-form map, not an amplitude (top of this\n"
      "       file).  A dipole-model run could move it by a factor.\n\n"
      "   So: MARGINAL at one EIC year, comfortably measurable at a few, and\n"
      "   not established either way until an efficiency exists for\n"
      "   Q^2 < 0.1, a decay-lepton efficiency exists at all, and a working\n"
      "   point is chosen.\n"
      % (tag.lumi_fraction,
         _s_at(band_lo_f * tag.lumi_fraction),
         _s_at(band_hi_f_max * tag.lumi_fraction),
         _l3_at(band_hi_f_max * tag.lumi_fraction),
         _l3_at(band_lo_f * tag.lumi_fraction),
         pess["significance"], pess["lumi_for_3sigma"],
         (2.0 / _head["significance"]) ** 2 / ch["beam_lo"],
         (_head["lumi_for_3sigma"] / 100.0) / ch["beam_lo"]))

    return dict(kappa_measured=k_meas, kappa_gfmc=k_gfmc, kappa_model=k_model,
                pzz2=pzz2, pzz2_difference=pzz2_diff, pzz_states=states,
                a2_eff_at_default_b=k_meas * math.sqrt(t2_mean(b_default, t_max)),
                jpsi=head_row, jpsi_optimal=head_opt,
                jpsi_q2_gt_0p1_ee=old_row,
                jpsi_conservative=cons, jpsi_pessimistic=pess,
                jpsi_rmeas=rmeas,
                lumi_for_3sigma_jpsi=head_row["lumi_for_3sigma"],
                lumi_for_3sigma_jpsi_q2_gt_0p1_ee=old_row["lumi_for_3sigma"],
                tagging_lumi_fraction=tag.lumi_fraction,
                eps_det=EPS_DET_CHANG,
                n_for_3sigma=n_for_significance(3.0, k_meas, pzz2_diff,
                                                b_default, t_max),
                residual_edge=resid_edge, eps_b0=eps_b0,
                eps_chain=ch,
                band_factor_lo=band_lo_f,
                # THE TOP IS A SPAN.  There is deliberately no scalar
                # `band_factor_hi` / `significance_hi` / `lumi_for_3sigma_lo`
                # key any more: quoting the top as one number is the defect
                # this pass removed, and a KeyError is the cheapest way to
                # stop it coming back.
                band_factor_hi_min=band_hi_f_min,
                band_factor_hi_max=band_hi_f_max,
                significance_lo=_s_at(band_lo_f),
                significance_hi_min=_s_at(band_hi_f_min),
                significance_hi_max=_s_at(band_hi_f_max),
                lumi_for_3sigma_top_min=_l3_at(band_hi_f_max),
                lumi_for_3sigma_top_max=_l3_at(band_hi_f_min),
                lumi_for_3sigma_hi=_l3_at(band_lo_f),
                significance_optics_lo=_s_at(band_lo_f * tag.lumi_fraction),
                significance_optics_hi_min=_s_at(band_hi_f_min
                                                 * tag.lumi_fraction),
                significance_optics_hi_max=_s_at(band_hi_f_max
                                                 * tag.lumi_fraction),
                lumi_for_3sigma_optics_lo=_l3_at(band_hi_f_max
                                                 * tag.lumi_fraction),
                eps_pair_for_2sigma=(2.0 / _head["significance"]) ** 2
                / ch["beam_lo"],
                eps_pair_for_band_exit=(_head["lumi_for_3sigma"] / 100.0)
                / ch["beam_lo"])


def _self_check(res: dict) -> None:
    """The assertions section C2 is written against.  Loose enough to survive
    a re-fit of an input, tight enough to catch a sign or a factor.

    RE-KEYED 2026-09-04 (FOURTH pass).  The third pass pinned
    `assert 3.0 < res["significance_hi"] < 3.5` -- a scalar top, taken from a
    species leg that was read off a list confounded with E/u and Z.  The top
    is now a SPAN and the pin is that the span STRADDLES 3 sigma.  It fails
    both ways, which is the point: if a future reading pushes the top's LOW
    reading above 3 sigma, "the top crosses 3 sigma" would have been
    re-asserted as established, and this assertion stops it.

      * the beam-energy factor must exceed 1 (it is a measurement of a
        direction, not a hedge) and stay near the 3He scan's own size;
      * the SPECIES factor must STRADDLE 1 on the same-energy entries -- its
        direction is not established, and every site says so;
      * the lepton factor must NOT exceed 1 (it can only cost);
      * beam energy and leptons must nearly cancel -- if a future input broke
        that, quoting one without the other would start to matter again and
        this file's section 3b(d) would be wrong;
      * 3 sigma must lie inside {1, 10, 100} at the low edge AND at both
        readings of the top;
      * the top must STRADDLE 3 sigma;
      * the de-squeezed optics must still take every rung outside the band.
    """
    assert abs(res["kappa_measured"] - 0.08752978) < 1e-7, res["kappa_measured"]
    assert res["kappa_measured"] > 0.0, "6Li's a_2 is POSITIVE (Q < 0)"
    assert abs(res["pzz2"] - 0.90) < 1e-12, res["pzz2"]
    assert abs(res["pzz2_difference"] - 0.81) < 1e-12, res["pzz2_difference"]
    assert abs(res["a2_eff_at_default_b"] - 0.0024723) < 1e-6
    # eps_det is still the PAPER's number; nothing here silently rescales it.
    assert abs(res["eps_det"] - 0.1775) < 1e-12, res["eps_det"]
    # The restricted branch: it is still under 1 sigma, and that is still the
    # right answer FOR THAT WINDOW AND ONE LEPTON CHANNEL.
    assert res["jpsi_q2_gt_0p1_ee"]["significance"] < 1.0
    assert res["lumi_for_3sigma_jpsi_q2_gt_0p1_ee"] > 100.0
    # The uncorrected headline row, still printed and still the ladder's top.
    assert 2.0 < res["jpsi"]["significance"] < 3.0, res["jpsi"]["significance"]
    assert 1.0 < res["lumi_for_3sigma_jpsi"] < 100.0, \
        "3 sigma must be INSIDE the {1, 10, 100} fb^-1/u band"

    # ---- THE CORRECTED eps_det CHAIN: two omissions, opposite directions.
    ch = res["eps_chain"]
    # (a) beam energy.  UP -- 0.1775 is a TOP-energy number used at 10 x 99.5.
    assert ch["beam_lo"] > 1.0 and ch["beam_hi"] > ch["beam_lo"], ch
    assert 1.10 < ch["beam_lo"] < 1.20, ch["beam_lo"]
    assert 1.10 < ch["beam_hi"] < 1.20, ch["beam_hi"]
    # ... and it is NOT the paper's own 183 -> 100 ratio of 1.687, which would
    # break the band open at the top.  If someone ever pastes that in, this
    # fails.
    assert ch["beam_hi"] < 1.4, "1.687 is the wrong lever arm; see 3b(a)"
    # (b) THE SPECIES LEG.  Read off the four entries that share a beam
    # energy, it STRADDLES 1: neither its size nor its direction is
    # established, and this is the whole of the band's top.
    assert ch["species_lo"] < 1.0 < ch["species_hi"], ch
    assert 0.95 < ch["species_lo"] < 1.05, ch["species_lo"]
    assert 1.25 < ch["species_hi"] < 1.40, ch["species_hi"]
    # ... and the CONFOUNDED reading it replaced is still printed, still
    # clears 1 at both ends, and is still bigger than the honest one at the
    # bottom.  If those ever coincided the retraction would be moot.
    conf = ch["species_confounded"]
    assert 1.10 < conf["lo"] < conf["hi"] < 1.30, conf
    assert conf["lo"] > ch["species_lo"]
    # THE 3He/4He TEST, with the RIGHT Z.  pT_cut = theta Z R, 7Li is Z = 3
    # and the pair is Z = 2, so (2/3)^2.  It then agrees to 0.1 %, and the
    # "overstates by 12 %" that the WRONG Z gave -- and the "read the low
    # end" that came out of it -- is withdrawn.
    assert abs(conf["he3_he4_ratio_z"] - 1.0) < 0.01, conf["he3_he4_ratio_z"]
    assert conf["he3_he4_overstatement_wrong_z"] > 1.10
    # the Z-INDEPENDENT reading fits the list's absolute values better, and
    # under it 3He -- the one entry off the family's E/u -- is the only miss,
    # by less than the shortfall the chain's own energy scan predicts for it.
    for n in ("2D", "4He", "9Be"):
        assert abs(conf["absolute_pred_over_meas"][n] - 1.0) < 0.06, n
    assert conf["absolute_pred_over_meas"]["3He"] > 1.10
    assert (conf["absolute_pred_over_meas"]["3He"]
            < 1.0 / conf["he3_energy_shortfall"][1])
    # (c) decay leptons.  DOWN, and never a gain.
    assert 0.0 < ch["a_geom"] <= 1.0, ch["a_geom"]
    assert 0.98 < ch["a_geom"] < 1.0, "geometry is not where the factor is"
    assert 0.0 < ch["down_standin"] < ch["a_geom"], ch
    # (d) they nearly cancel -- section 3b(d) and the whole point of stating
    # both.  If this ever stopped holding, quoting one alone would move the
    # verdict and the prose would be wrong.
    assert abs(ch["beam_lo"] * ch["down_standin"] - 1.0) < 0.10, \
        "3b(d) claims the two omissions nearly cancel"

    # ---- THE VERDICT, as a BAND whose TOP IS ITSELF A SPAN.
    assert res["significance_lo"] < res["significance_hi_min"]
    assert res["significance_hi_min"] < res["significance_hi_max"]
    assert 2.0 < res["significance_lo"] < 3.0, res["significance_lo"]
    # THE PIN THAT STOPS "the band's TOP just crosses 3 sigma" COMING BACK.
    # The top STRADDLES 3 sigma; it does not clear it.  Re-establishing the
    # crossing means pushing the top's LOW reading above 3, and that fails
    # here.
    assert res["significance_hi_min"] < 3.0 < res["significance_hi_max"], \
        "the band's TOP straddles 3 sigma -- it is NOT established that it " \
        "crosses it (section 3b(b), 7)"
    assert 2.5 < res["significance_hi_min"], res["significance_hi_min"]
    assert res["significance_hi_max"] < 3.6, res["significance_hi_max"]
    # ... and there is deliberately no scalar top to quote.
    for gone in ("significance_hi", "lumi_for_3sigma_lo", "band_factor_hi"):
        assert gone not in res, gone
    for key in ("lumi_for_3sigma_top_min", "lumi_for_3sigma_top_max",
                "lumi_for_3sigma_hi"):
        assert 1.0 < res[key] < 100.0, (key, res[key])
    # THE CANCELLATION, stated as a number: the UNCORRECTED point sits within
    # 1 % of the band's low edge (in fact 0.3 % under it).  If a future input
    # broke that, every site's "they nearly cancel" sentence would be wrong.
    assert abs(res["significance_lo"] / res["jpsi"]["significance"] - 1.0) < 0.01
    # the band is open BELOW, and what it takes to leave the {1,10,100} band
    # on the lepton factor alone is a detector reconstructing 1 J/psi in ~8
    assert 0.4 < res["eps_pair_for_2sigma"] < 0.7
    assert res["eps_pair_for_band_exit"] < 0.15

    # EACH pessimism ON ITS OWN leaves 3 sigma inside the band -- EXCEPT the
    # far-forward optics.  That exception is the correction stated as the one
    # that restores a NO, and it must stay visible at BOTH ends of the band.
    for key in ("jpsi_conservative", "jpsi_rmeas"):
        assert res[key]["lumi_for_3sigma"] < 100.0, key
    # STACKING both of them puts it back outside, and that must stay visible:
    # the verdict is MARGINAL, not YES.
    assert res["jpsi_pessimistic"]["lumi_for_3sigma"] > 100.0
    assert res["tagging_lumi_fraction"] < 0.1
    assert res["lumi_for_3sigma_jpsi"] / res["tagging_lumi_fraction"] > 100.0
    assert res["lumi_for_3sigma_optics_lo"] > 100.0, \
        "the optics must take even the most favourable TOP rung outside " \
        "the {1, 10, 100} band"
    assert res["significance_optics_hi_max"] < 1.0
    assert (res["significance_optics_lo"] < res["significance_optics_hi_min"]
            < res["significance_optics_hi_max"])
    print("\nself-check: OK")


if __name__ == "__main__":
    _self_check(report())
