"""Open item O5 -- the 6Li tensor-a_2 feasibility arithmetic, pinned.

The estimate itself lives in `validation/o5_a2_reach.py` (imported here, never
retyped) and is written up in `docs/open_items/run_2026-09-03/phase_C_numbers.md`
section C2 and `docs/OPEN_ITEMS_SOLUTIONS.md` section 11.  This file exists so
that moving ANY of its inputs -- `LI6_QUADRUPOLE_FM2`, the
`estarlight_li6_coherent()` / `estarlight_li6_q2_floors()` cross sections,
`COHERENT_JPSI_EFF_IR8_LI7`, `Scenario::lumi_fb_per_nucleon` or
`tensor_flip_plan` -- fails a test instead of silently rotting a document.

It also gates the two `estarlight_li6_*` bindings, whose tables are the SINGLE
code home of the eSTARlight numbers `estarlight_li6.md` reports.

2026-09-04 (third pass): the verdict is a BAND, not a point.  Two factors were
missing from the efficiency chain and they pull opposite ways -- the beam
energy `COHERENT_JPSI_EFF_IR8_LI7` was measured at (UP, measured) and the
decay-lepton acceptance x efficiency that was never in the chain at all (DOWN,
geometry bounded here, reconstruction efficiency not).  Both are pinned, in
both directions, together with the fact that they nearly cancel; and the
"every conservatism on its own leaves 3 sigma inside the band" claim, which
was FALSE at three sites (the optics row is a counterexample and was printed
three lines above it), is pinned in its corrected form.

2026-09-04: the verdict tests below used to pin "J/psi never reaches 1 sigma"
and "3 sigma is off the top of the band" from the 0.1 < Q^2 < 100, e+e--only
row -- a row that cannot fail those assertions no matter what the rest of the
chain does.  Both directions are pinned now: that row is still under 1 sigma,
AND the headline (whole Q^2 range, both leptons) is between 2 and 3 sigma with
3 sigma inside the {1, 10, 100} fb^-1/u band.
"""

import io
import math
import os
import sys

import pytest

import lipolgen as lg


@pytest.fixture(scope="module")
def o5(repo_root):
    sys.path.insert(0, os.path.join(repo_root, "validation"))
    try:
        import o5_a2_reach
    finally:
        sys.path.pop(0)
    return o5_a2_reach


# ------------------------------------------------- the new binding's table

def test_estarlight_li6_table_is_the_single_home_of_those_cross_sections():
    rows = {r.vm: r for r in lg.estarlight_li6_coherent()}
    assert set(rows) == {"jpsi", "phi", "rho"}
    # estarlight_li6.md sec. 2a (0.1 < Q^2 < 100, R_G = 1.2 A^(1/3)).
    assert rows["jpsi"].sigma_nb == pytest.approx(1.773)
    assert rows["phi"].sigma_nb == pytest.approx(30.16)
    assert rows["rho"].sigma_nb == pytest.approx(506.4)
    # sec. 5 (the same run restricted to LiPolGen's own Q^2 > 0.7 window).
    assert rows["jpsi"].sigma_q7_nb == pytest.approx(0.605)
    assert rows["phi"].sigma_q7_nb == pytest.approx(1.897)
    assert rows["rho"].sigma_q7_nb == pytest.approx(15.32)
    # sec. 2b (measured Angeli-Marinova charge radius 2.589 fm).
    assert rows["jpsi"].sigma_rmeas_nb == pytest.approx(1.255)
    assert rows["phi"].sigma_rmeas_nb == pytest.approx(22.08)
    assert rows["rho"].sigma_rmeas_nb == pytest.approx(379.5)
    # The slope is VM-INDEPENDENT -- sec. 4's one genuinely new fact.  The
    # spread quoted there as "1.5 %" is (max - min)/mean = 1.54 %.
    b = [r.b_fit for r in lg.estarlight_li6_coherent()]
    assert (max(b) - min(b)) / (sum(b) / 3.0) == pytest.approx(0.0154,
                                                              abs=2e-4)
    # ... and each fitted B undershoots the analytic R_G^2/(3 hbar^2c^2) by
    # 3.7-5.2 % (NOT the "4 %" the prose rounded to: that is the phi row
    # alone; rho is 5.2 % low).  One-sided, as the |t_min| floor and the
    # electroproduction flux require.
    b_analytic = lg.gaussian_slope(1.2 * 6.0 ** (1.0 / 3.0))
    assert b_analytic == pytest.approx(40.7026, abs=1e-3)
    dev = sorted(v / b_analytic - 1.0 for v in b)
    assert dev[0] == pytest.approx(-0.0517, abs=5e-4)     # rho
    assert dev[-1] == pytest.approx(-0.0369, abs=5e-4)    # phi
    assert all(d < 0.0 for d in dev)
    # The measured-radius density LOWERS sigma and RAISES B: the two are not
    # independent knobs (sec. 2b, the +-30 % leading rate systematic).
    for r in lg.estarlight_li6_coherent():
        assert r.sigma_rmeas_nb < r.sigma_nb
    assert rows["jpsi"].branching == pytest.approx(0.0597)   # J/psi -> e+e-
    assert rows["phi"].branching == pytest.approx(0.49)      # phi -> K+K-
    assert rows["rho"].branching == pytest.approx(1.0)       # rho -> pi pi
    # sec. 2b's slope column, which before 2026-09-04 lived only in a dict
    # inside validation/o5_a2_reach.py and could not be reached by a re-fit.
    assert rows["jpsi"].b_rmeas == pytest.approx(55.0)
    assert rows["phi"].b_rmeas == pytest.approx(54.8)
    assert rows["rho"].b_rmeas == pytest.approx(54.1)
    for r in lg.estarlight_li6_coherent():
        assert r.b_rmeas > r.b_fit       # bigger R -> steeper slope
    # ... and the SECOND lepton channel, likewise.  eSTARlight's own
    # JpsiBree + JpsiBrmumu = 0.05971 + 0.05961.
    assert rows["jpsi"].branching_all == pytest.approx(0.11932)
    assert rows["phi"].branching_all == pytest.approx(0.49)   # no 2nd channel
    assert rows["rho"].branching_all == pytest.approx(1.0)
    assert math.sqrt(rows["jpsi"].branching_all
                     / rows["jpsi"].branching) == pytest.approx(1.41374,
                                                                abs=1e-5)
    # The efficiency is a named constant, not a literal in a script.
    assert lg.COHERENT_JPSI_EFF_IR8_LI7 == pytest.approx(0.1775)
    # ... and so are the beam energy and <W> it was MEASURED at, which is the
    # whole subject of the 2026-09-04 third pass: 0.1775 is a TOP-energy 7Li
    # number and the chain applies it to 6Li at 10 x 99.5.
    assert lg.COHERENT_JPSI_EFF_IR8_LI7_E_ION_GEV == pytest.approx(117.9)
    assert lg.COHERENT_JPSI_EFF_IR8_LI7_W_MEAN_GEV == pytest.approx(43.2)
    top7 = lg.default_configs("7Li")[2]
    assert top7.ion_momentum_per_nucleon == pytest.approx(117.9)
    assert top7.electron_energy == pytest.approx(18.0)


def test_the_photoproduction_scan_is_most_of_the_coherent_rate():
    """The 0.1 GeV^2 floor is arXiv:2511.05638's ACCEPTANCE-STUDY range copied
    verbatim, not a physics window.  Below it lives most of the rate."""
    q = list(lg.estarlight_li6_q2_floors())
    assert len(q) == 9                         # 3 mesons x 3 floors
    at = {(r.vm, r.q2_floor_gev2): r for r in q}
    base = {r.vm: r for r in lg.estarlight_li6_coherent()}
    # The Q^2 > 0.1 rows ARE sec. 2a / 2b, bit for bit -- that is what makes
    # the rest of the scan comparable to the numbers already published.
    for vm, b in base.items():
        r = at[(vm, 0.1)]
        assert r.sigma_nb == b.sigma_nb
        assert r.b_fit == b.b_fit
        assert r.sigma_rmeas_nb == b.sigma_rmeas_nb
        assert r.b_rmeas == b.b_rmeas
    # Removing the floor multiplies the rate; the lighter the meson the more,
    # because its Q^2 suppression sets in earlier.
    for vm, want_full, want_001 in (("jpsi", 6.7518, 1.8883),
                                    ("phi", 21.696, 3.4213),
                                    ("rho", 35.195, 4.5119)):
        a, b, c = at[(vm, 0.1)], at[(vm, 0.01)], at[(vm, 0.0)]
        assert c.sigma_nb / a.sigma_nb == pytest.approx(want_full, rel=1e-3)
        assert b.sigma_nb / a.sigma_nb == pytest.approx(want_001, rel=1e-3)
        assert c.sigma_nb > b.sigma_nb > a.sigma_nb
        assert c.sigma_rmeas_nb > b.sigma_rmeas_nb > a.sigma_rmeas_nb
        assert c.sigma_rmeas_nb < c.sigma_nb     # bigger R -> smaller sigma
        assert c.b_rmeas > c.b_fit               # ... and steeper slope
        # The |t| slope barely moves with the floor: the recoil p_T spectrum
        # the far-forward acceptance cuts on is the SAME sample.
        assert abs(c.b_fit / a.b_fit - 1.0) < 0.015
        assert abs(c.b_rmeas / a.b_rmeas - 1.0) < 0.015
        # <W> falls with the floor, on every meson -- the W-axis half of the
        # transfer argument, in code since 2026-09-04 instead of in prose.
        assert c.w_mean_gev < a.w_mean_gev
    assert at[("jpsi", 0.0)].w_mean_gev == pytest.approx(30.20)
    assert at[("jpsi", 0.1)].w_mean_gev == pytest.approx(32.20)
    # 85 % of the coherent J/psi rate is below Q^2 = 0.1 GeV^2.
    j = at[("jpsi", 0.0)], at[("jpsi", 0.1)]
    assert 1.0 - j[1].sigma_nb / j[0].sigma_nb == pytest.approx(0.852,
                                                               abs=2e-3)


# ------------------------------------------- the efficiency chain, corrected

def test_the_efficiency_is_not_flat_in_beam_energy_and_the_paper_says_so(o5):
    """DEFECT 1.  `COHERENT_JPSI_EFF_IR8_LI7` = 0.1775 is a 7Li number at 7Li's
    OWN TOP ENERGY, 18 x 117.9 GeV/u, and the chain applies it to a 6Li sample
    at 10 x 99.5.  arXiv:2511.05638 sec. V.B measures exactly that dependence
    on 3He, and the tree records it as `chang26_he3_energy_scan()`."""
    e = list(lg.chang26_he3_energy_scan())
    assert len(e) == 3
    assert [(r.e_electron_gev, r.e_ion_gev, r.efficiency) for r in e] == [
        (18.0, 183.0, pytest.approx(0.3223)),
        (10.0, 100.0, pytest.approx(0.5438)),
        (5.0, 41.0, pytest.approx(0.9977))]
    # MONOTONE, and that direction is the whole correction: lower beam energy,
    # HIGHER far-forward efficiency.
    for a, b in zip(e, e[1:]):
        assert b.e_ion_gev < a.e_ion_gev
        assert b.efficiency > a.efficiency
        assert b.efficiency <= 1.0
    s1, s2 = o5.he3_energy_slopes()
    assert s1 == pytest.approx(-0.8656, abs=1e-3)
    assert s2 == pytest.approx(-0.6807, abs=1e-3)
    assert s1 < s2 < 0.0          # flattening, as the cap at 1 requires

    lo, hi = o5.beam_energy_scaling()
    assert lo == pytest.approx(1.12243, abs=1e-4)
    assert hi == pytest.approx(1.15821, abs=1e-4)
    assert 1.0 < lo < hi
    # ... and it is NOT the paper's own 183 -> 100 ratio of 1.687.  That is a
    # step 3.56x larger in ln E; using it would be the wrong lever arm, and
    # would push the headline to 3.40 sigma / 7.8 fb^-1/u.
    assert e[1].efficiency / e[0].efficiency == pytest.approx(1.6872, abs=1e-3)
    assert hi < 1.4
    step_ratio = (math.log(e[1].e_ion_gev / e[0].e_ion_gev)
                  / math.log(lg.default_configs("6Li")[1].ion_momentum_per_nucleon
                             / lg.COHERENT_JPSI_EFF_IR8_LI7_E_ION_GEV))
    assert step_ratio == pytest.approx(3.56, abs=0.02)


def test_the_species_list_is_fixed_rigidity_not_fixed_energy(o5):
    """THE RETRACTED SENTENCE.  Until 2026-09-04 five sites read "every entry
    is at the same rigidity, so its A-ordering is a species lever on its own".
    Eliminating rigidity does NOT leave species alone: at fixed R = A E/Z both
    E/u = R Z/A and p_z = Z R still vary down the list, and the chain's own
    energy slope makes the E/u spread worth as much as the whole claimed
    species gain."""
    t = list(lg.chang26_species_efficiency())
    assert len(t) == 7
    for r in t:
        assert r.a / r.z * r.e_ion_gev == pytest.approx(275.0, abs=1.5)
    for a, b in zip(t, t[1:]):
        assert b.a > a.a
        assert b.efficiency < a.efficiency
    by = {r.nucleus: r for r in t}
    assert by["7Li"].efficiency == pytest.approx(lg.COHERENT_JPSI_EFF_IR8_LI7)
    assert by["3He"].efficiency == pytest.approx(
        lg.chang26_he3_energy_scan()[0].efficiency)

    nl = o5.species_list_is_not_a_species_lever()
    assert nl["e_per_u"] == (118.0, 183.0)      # 7Li ... 3He, at fixed R
    assert nl["e_per_u_spread"] == pytest.approx(1.5508, abs=1e-3)
    assert nl["z"] == (1, 8)
    assert nl["p_z"] == (274.0, 2192.0)
    # the confound is as large as the effect that was being read off the list
    assert nl["eff_worth_of_e_spread"][0] == pytest.approx(1.3481, abs=1e-3)
    assert nl["eff_worth_of_e_spread"][1] == pytest.approx(1.4620, abs=1e-3)
    assert nl["eff_worth_of_e_spread"][1] > 1.22
    assert nl["n_at_same_energy"] == 4           # 2D, 4He, 12C, 16O at 137

    # THE CONFOUNDED READING, kept only so the confound stays visible.
    sp = o5.species_scaling()
    assert 1.10 < sp["lo"] < sp["hi"] < 1.30
    assert sp["lo"] == pytest.approx(1.1603, abs=1e-3)
    assert sp["hi"] == pytest.approx(1.2192, abs=1e-3)

    # THE READING WITH NO ENERGY STEP IN IT, and the one the band uses:
    # 4He -> 12C, both at 137 GeV/u, bracketing 6Li's own 137.5.  It
    # STRADDLES 1 -- the direction of the substitution is not established.
    se = o5.species_scaling_same_energy()
    assert se["nuclei"] == ["2D", "4He", "12C", "16O"]
    assert se["at_energy_gev"] == 137.0
    assert se["lo"] == pytest.approx(0.9932, abs=1e-3)
    assert se["hi"] == pytest.approx(1.3327, abs=1e-3)
    assert se["lo"] < 1.0 < se["hi"]
    assert se["eps_li6"]["log_r"] == pytest.approx(0.17823, abs=1e-4)
    assert se["eps_li6"]["lin_a"] == pytest.approx(0.23655, abs=1e-4)
    assert se["eps_li6"]["pt_threshold@12C"] == pytest.approx(0.17629,
                                                              abs=1e-4)


def test_the_one_testable_form_is_read_three_ways_and_none_of_them_wins(o5):
    """A2.  The model is eps = exp(-B(A) pT_cut^2) and the criterion is a safe
    distance FROM THE BEAM, so pT_cut = theta p_z = theta Z R -- it carries
    the CHARGE.  The shipped test applied 7Li's Z = 3 cut to the Z = 2 3He/4He
    pair and called the 12 % that produced an overstatement, then read "the
    low end" off it.  That instruction is withdrawn."""
    sp = o5.species_scaling()
    assert sp["pt_cut_gev"] == pytest.approx(0.19577, abs=1e-4)
    assert sp["he3_he4_measured"] == pytest.approx(1.0955, abs=1e-3)
    # (1) with the RIGHT Z: (2/3)^2 on pT_cut^2, and it agrees to 0.1 %
    assert sp["he3_he4_predicted_z"] == pytest.approx(1.0967, abs=1e-3)
    assert sp["he3_he4_ratio_z"] == pytest.approx(1.0011, abs=1e-3)
    assert abs(sp["he3_he4_ratio_z"] - 1.0) < 0.01
    # (2) the shipped reading, with the WRONG Z
    assert sp["he3_he4_predicted_wrong_z"] == pytest.approx(1.2309, abs=1e-3)
    assert sp["he3_he4_overstatement_wrong_z"] == pytest.approx(1.1236,
                                                                abs=1e-3)
    # (3) Z-INDEPENDENT pT_cut fits the ABSOLUTE values far better than (1)
    a, az = sp["absolute_pred_over_meas"], sp["absolute_pred_over_meas_z"]
    assert a["2D"] == pytest.approx(1.0025, abs=1e-3)
    assert a["4He"] == pytest.approx(1.0336, abs=1e-3)
    assert a["9Be"] == pytest.approx(1.0469, abs=1e-3)
    for n in ("2D", "4He", "9Be"):
        assert abs(a[n] - 1.0) < abs(az[n] - 1.0)
    # ... and the one misfit among the LIGHT entries is 3He, by less than the
    # shortfall the chain's OWN energy scan predicts for its E/u = 183
    assert a["3He"] == pytest.approx(1.1613, abs=1e-3)
    lo, hi = sp["he3_energy_shortfall"]
    assert (1.0 - hi) == pytest.approx(0.179, abs=2e-3)
    assert (1.0 - lo) == pytest.approx(0.222, abs=2e-3)
    assert a["3He"] < 1.0 / hi
    # AND SAY WHAT IT IS NOT.  The form is LOCAL, not a description of the
    # list: it degrades with A and breaks down by 12C / 16O.  Pinned so the
    # "it fits the absolute values" claim can never be widened to the list.
    assert a["12C"] == pytest.approx(1.3217, abs=1e-3)
    assert a["16O"] == pytest.approx(3.1320, abs=1e-3)
    assert a["12C"] > 1.20 and a["16O"] > 2.0


def test_there_is_a_decay_lepton_factor_and_it_can_only_cost(o5):
    """DEFECT 2.  arXiv:2511.05638's number is the fraction of scattered
    NUCLEI inside the far-forward acceptance and its own text says it carries
    no detector efficiency and no reconstructed-distribution acceptance.  So
    sigma x BR(l+l-) x eps_recoil is missing A_ll x eps_ll <= 1 entirely.

    What is BOUNDED here is the geometry, and it is not the problem: at
    10 x 99.5 the J/psi is nearly at rest in the lab.  What is NOT bounded
    anywhere in this tree is the per-lepton reconstruction efficiency."""
    w = o5.jpsi_w_mean(0.0)
    assert w == pytest.approx(30.20)
    y, e_gamma = o5.jpsi_lab_rapidity(w)
    assert e_gamma == pytest.approx(2.289, abs=2e-3)
    assert y == pytest.approx(-0.391, abs=2e-3)
    assert abs(y) < 0.5                        # very nearly at rest
    a = o5.lepton_pair_acceptance(w)
    assert a == pytest.approx(0.9940, abs=1e-3)
    assert a <= 1.0
    # it uses THIS repository's own central-detector edge, not an imported one
    assert o5.lepton_pair_acceptance(w, eta_max=lg.Scenario().eta_max) == a
    assert lg.Scenario().eta_max == 3.5
    # SCHC is the pessimistic weighting; isotropic is larger
    assert o5.lepton_pair_acceptance(w, schc=False) > a
    # a track p_T threshold costs a little more, and it is the SIBLING's
    # stand-in rather than an in-tree number
    assert o5.lepton_pair_acceptance(w, pt_min=0.2) < a
    # ... and it degrades monotonically as the J/psi moves off mid-rapidity,
    # never below 0.89 anywhere the beams can reach (E_gamma <= E_e)
    accs = [o5.lepton_pair_acceptance(x) for x in (20.0, 30.2, 40.0, 50.0, 63.1)]
    assert min(accs) > 0.89
    assert accs[-1] < accs[1]


def test_the_two_omissions_pull_opposite_ways_and_nearly_cancel(o5):
    """The reason both had to be stated: a reader given only one of them moves
    the verdict, and the pair of them barely move it at all."""
    ch = o5.eps_det_chain()
    assert ch["beam_lo"] > 1.0                       # UP
    assert ch["down_standin"] < 1.0                  # DOWN
    assert ch["a_geom"] <= 1.0                       # an acceptance
    assert ch["down_standin"] < ch["a_geom"]         # ... times an efficiency
    prod = ch["beam_lo"] * ch["down_standin"]
    assert prod == pytest.approx(1.0070, abs=2e-3)
    assert abs(prod - 1.0) < 0.10
    # the two-leg ends, which is where the band's top comes from -- and the
    # species leg straddles 1, so the low end is now BELOW the energy leg
    assert ch["two_leg_lo"] == pytest.approx(1.2378, abs=1e-3)
    assert ch["two_leg_hi"] == pytest.approx(1.7633, abs=1e-3)
    assert ch["two_leg_top_min"] == pytest.approx(1.3141, abs=1e-3)
    assert ch["two_leg_top_max"] == pytest.approx(1.7633, abs=1e-3)
    assert ch["species_lo"] < 1.0 < ch["species_hi"]
    assert ch["species_confounded_lo"] == pytest.approx(1.1603, abs=1e-3)
    assert ch["eps_det_shipped"] == pytest.approx(0.1775)


# ------------------------------------------------------------- the signal

def test_the_signal_is_positive_small_and_linear_in_the_measured_quadrupole(o5):
    k = o5.kappa_of_quadrupole(lg.LI6_QUADRUPOLE_FM2)
    assert k == pytest.approx(0.08752978, abs=1e-8)
    assert k > 0.0                       # OPPOSITE sign to the deuteron's
    assert 0.3 * k == pytest.approx(0.026259, abs=1e-6)
    # kappa is exactly linear in Q, so the 7.5x quadrupole gap is a 7.5x
    # signal gap: the model's own Q and GFMC's ride the same line.
    q_model = 0.5 * lg.ClusterConfigSampler().q_matter_analytic_fm2(1)
    assert o5.kappa_of_quadrupole(q_model) / k == pytest.approx(
        q_model / lg.LI6_QUADRUPOLE_FM2, rel=1e-12)
    assert o5.kappa_of_quadrupole(lg.LI6_QUADRUPOLE_GFMC_FM2) == pytest.approx(
        0.21400924, abs=1e-8)


def test_the_coherent_t_slope_is_what_kills_it(o5):
    """a_2 is quoted at |t| = 0.3 but the sample lives at |t| ~ 1/B."""
    b = lg.CoherentScenario().slope_b
    tmax = lg.COHERENT_T_MAX_DEFAULT
    k = o5.kappa_of_quadrupole(lg.LI6_QUADRUPOLE_FM2)
    a2_eff = k * math.sqrt(o5.t2_mean(b, tmax))
    assert a2_eff == pytest.approx(0.0024723, abs=1e-6)
    # a factor 10.6 below the a_2(0.3) = 0.0263 that gets quoted
    assert 0.3 * k / a2_eff == pytest.approx(10.62, abs=0.02)
    # and the |t| window itself costs essentially nothing
    assert o5.t_fraction(b, tmax) > 0.9999
    assert o5.info_fraction(b, tmax) == pytest.approx(0.99723, abs=1e-4)


# ---------------------------------------------------------- the run plan

def test_the_pzz_flip_plan_is_a_gain_not_a_cost(o5):
    pzz2, states = o5.flip_plan_pzz2()
    assert states == pytest.approx([0.6, -1.2])
    # <P_zz^2> = 2.5 pzz^2 exactly: 0.5 (pzz^2 + 4 pzz^2)
    assert pzz2 == pytest.approx(0.90, abs=1e-12)
    assert pzz2 == pytest.approx(2.5 * 0.6 ** 2, rel=1e-12)
    # The background-immune DIFFERENCE estimator against (a) the optimal use
    # of the same two fills and (b) one single +pzz fill of the same total
    # luminosity.  It costs 5.4 % and gains 50 %.
    p1, p2 = states
    v_diff = (2.0 / 0.5 + 2.0 / 0.5) / ((2.0 * (p1 - p2)) ** 2)
    v_opt = 1.0 / (2.0 * (0.5 * p1 * p1 + 0.5 * p2 * p2))
    v_single = 1.0 / (2.0 * p1 * p1)
    assert math.sqrt(v_diff / v_opt) == pytest.approx(1.0541, abs=1e-3)
    assert math.sqrt(v_diff / v_single) == pytest.approx(2.0 / 3.0, abs=1e-12)


def test_delta_kappa_is_the_repositorys_own_cos2phi_error(o5):
    """The Fisher form must agree with `estimators.cos2phi_fit_err` in one
    narrow |t| bin -- otherwise the estimate is using a different statistics
    from the analysis chain it is meant to price."""
    n, pzz, t = 1.0e5, 0.6, 0.05
    # delta(c_2/P_zz) = sqrt(2/N)/P_zz  ->  delta(kappa) = that / (2 t)
    from_repo = lg.estimators.cos2phi_fit_err(n, pzz, 24) / (2.0 * t)
    # ... at 24 bins the repo inflates by 1/dilution; undo it for the compare
    w = math.pi / 24.0
    from_repo *= math.sin(2.0 * w) / (2.0 * w)
    mine = 1.0 / math.sqrt(2.0 * pzz * pzz * n * t * t)
    assert from_repo == pytest.approx(mine, rel=1e-12)


# ---------------------------------------------------------- the verdict

def test_the_restricted_row_is_still_the_number_it_always_was(o5):
    """0.1 < Q^2 < 100, J/psi -> e+e- only, OPTIMAL <P_zz^2> = 0.90: the row
    open item O5 shipped with.  It is still 0.749 sigma and 160 fb^-1/u, and
    it is still arithmetically right FOR THAT ROW.  What was wrong was calling
    it the measurement."""
    rows = o5.channel_rows(eps_det=o5.EPS_DET_CHANG, pzz2=o5.PZZ2_OPTIMAL)
    j = o5._row(rows, "jpsi", "Q2>0.1")
    assert j["n_detected"] == pytest.approx(3.13e4, rel=2e-3)
    assert j["delta_a2_at_0p3"] == pytest.approx(0.03505, abs=1e-4)
    assert j["significance"] == pytest.approx(0.749, abs=2e-3)
    assert j["lumi_for_3sigma"] == pytest.approx(160.3, abs=1.0)
    assert j["lumi_for_3sigma"] > 100.0


def test_the_verdict_is_marginal_and_three_sigma_is_inside_the_band(o5):
    """THE FIX.  Both corrections that were missing point the same way: the
    second lepton channel (sqrt 2) and the photoproduction region (sqrt 6.75).
    Together they take J/psi from 0.75 sigma to 2.6, and 3 sigma from 160
    fb^-1/u to 13 -- inside the {1, 10, 100} band."""
    rows = o5.channel_rows(eps_det=o5.EPS_DET_CHANG, both_leptons=True)
    j = o5._row(rows, "jpsi", "Q2>0")
    assert j["pzz2"] == pytest.approx(0.81)      # background-immune default
    assert j["n_detected"] == pytest.approx(4.224e5, rel=2e-3)
    assert j["significance"] == pytest.approx(2.618, abs=3e-3)
    assert 2.0 < j["significance"] < 3.0         # MARGINAL: not NO, not YES
    assert j["lumi_for_3sigma"] == pytest.approx(13.13, abs=0.05)
    assert 1.0 < j["lumi_for_3sigma"] < 100.0
    # Each correction on its own, so neither can be quietly dropped again.
    both_only = o5._row(rows, "jpsi", "Q2>0.1")
    photo_only = o5._row(o5.channel_rows(eps_det=o5.EPS_DET_CHANG),
                         "jpsi", "Q2>0")
    assert both_only["significance"] == pytest.approx(1.005, abs=3e-3)
    assert photo_only["significance"] == pytest.approx(1.852, abs=3e-3)
    # ... and they multiply: sqrt(2) x sqrt(6.752) on the 0.1-window row
    base = o5._row(o5.channel_rows(eps_det=o5.EPS_DET_CHANG), "jpsi", "Q2>0.1")
    assert (j["significance"] / base["significance"]) == pytest.approx(
        math.sqrt(1.99866 * 6.7518) * math.sqrt(
            o5.t2_mean(38.8, lg.COHERENT_T_MAX_DEFAULT)
            / o5.t2_mean(38.9, lg.COHERENT_T_MAX_DEFAULT)), rel=2e-3)


def test_the_verdict_survives_each_pessimism_but_not_two_of_them(o5):
    """What makes the answer MARGINAL rather than YES: every conservatism on
    its own leaves 3 sigma inside the band EXCEPT the far-forward optics, and
    stacking two of the others puts it back out.

    The exception matters.  Until 2026-09-04 three sites printed "every single
    conservatism on its own leaves 3 sigma inside the band" -- in the script,
    three lines below the bullet saying the optics row does not.  It was
    false, and this test now pins the corrected form."""
    rows = o5.channel_rows(eps_det=o5.EPS_DET_CHANG, both_leptons=True)
    # the conservative Q^2 > 0.01 floor alone: inside
    assert o5._row(rows, "jpsi", "Q2>0.01")["lumi_for_3sigma"] < 100.0
    # the measured-radius density alone: inside
    assert o5._row(rows, "jpsi", "Q2>0,Rmeas")["lumi_for_3sigma"] < 100.0
    # both at once: OUT.  This is why the verdict is not YES.
    assert o5._row(rows, "jpsi", "Q2>0.01,Rmeas")["lumi_for_3sigma"] > 100.0
    # the de-squeezed far-forward tagging optics alone: OUT, and it is the
    # one correction that is a CHOICE rather than a measurement.
    tag = lg.tagging_optics("6Li", 99.5).lumi_fraction
    assert tag < 0.1
    head = o5._row(rows, "jpsi", "Q2>0")
    assert head["lumi_for_3sigma"] / tag > 100.0
    assert lg.yr_optics("6Li", 99.5).lumi_fraction == 1.0
    # THE EXCEPTION, stated as such: the optics alone is the one correction
    # that leaves the band, and it does so at BOTH ends of the corrected
    # eps_det band, not only on the uncorrected point.
    ch = o5.eps_det_chain()
    band_hi = ch["two_leg_top_max"] * ch["down_standin"]
    assert head["lumi_for_3sigma"] / band_hi / tag > 100.0
    assert head["significance"] * math.sqrt(band_hi * tag) < 1.0


def test_the_verdict_is_a_band_whose_top_is_itself_a_span(o5):
    """THE CLAIM, as the FOURTH pass leaves it.  Not a point, and its TOP is
    not an edge: the dominant multiplier eps_det carries a measured UP factor
    (beam energy), an unbounded DOWN factor (the decay-lepton pair) and a
    species factor whose DIRECTION is not established.  The third pass pinned
    `significance_hi == 3.149`; that number came from reading the species leg
    off a list confounded with E/u and Z.  Read off the four entries that
    share a beam energy the leg straddles 1, so the top straddles 3 sigma.

    THE ASSERTION BELOW FAILS IF "the top crosses 3 sigma" IS EVER
    RE-ASSERTED AS ESTABLISHED: that would mean the top's LOW reading had
    been pushed above 3."""
    buf = io.StringIO()
    res = o5.report(out=buf)
    o5._self_check(res)
    assert res["significance_lo"] == pytest.approx(2.627, abs=5e-3)
    assert res["significance_hi_min"] == pytest.approx(2.842, abs=5e-3)
    assert res["significance_hi_max"] == pytest.approx(3.292, abs=5e-3)
    assert 2.0 < res["significance_lo"] < 3.0
    assert res["significance_hi_min"] < 3.0 < res["significance_hi_max"]
    # there is no scalar top left to quote
    for gone in ("significance_hi", "lumi_for_3sigma_lo", "band_factor_hi"):
        assert gone not in res
    assert res["lumi_for_3sigma_top_min"] == pytest.approx(8.31, abs=0.05)
    assert res["lumi_for_3sigma_top_max"] == pytest.approx(11.15, abs=0.05)
    assert res["lumi_for_3sigma_hi"] == pytest.approx(13.04, abs=0.05)
    for k in ("lumi_for_3sigma_top_min", "lumi_for_3sigma_top_max",
              "lumi_for_3sigma_hi"):
        assert 1.0 < res[k] < 100.0        # inside the band on EVERY form
    # THE CANCELLATION: the uncorrected point lands 0.3 % UNDER the band's low
    # edge -- not inside it, and the prose says so at every site.
    ratio = res["significance_lo"] / res["jpsi"]["significance"]
    assert ratio == pytest.approx(1.0035, abs=5e-4)
    assert 1.0 < ratio < 1.01
    # the band is OPEN below: what it takes to leave {1, 10, 100} on the
    # lepton factor alone is a detector reconstructing ~1 J/psi in 8
    assert res["eps_pair_for_2sigma"] == pytest.approx(0.520, abs=5e-3)
    assert res["eps_pair_for_band_exit"] == pytest.approx(0.117, abs=2e-3)
    # and the optics, at both ends, is still a NO
    assert res["lumi_for_3sigma_optics_lo"] > 100.0
    assert res["significance_optics_hi_max"] < 1.0
    assert res["significance_optics_hi_min"] == pytest.approx(0.794, abs=5e-3)
    assert res["significance_optics_hi_max"] == pytest.approx(0.920, abs=5e-3)
    assert res["lumi_for_3sigma_optics_lo"] == pytest.approx(106.4, abs=0.5)


def test_a_perfect_detector_is_an_upper_bound_not_a_projection(o5):
    perfect = o5.channel_rows(eps_det=1.0, both_leptons=True)
    j = o5._row(perfect, "jpsi", "Q2>0")
    assert j["significance"] == pytest.approx(6.213, abs=5e-3)
    # ... exactly 1/sqrt(eff) above the realistic row
    real = o5._row(o5.channel_rows(eps_det=o5.EPS_DET_CHANG,
                                   both_leptons=True), "jpsi", "Q2>0")
    assert j["significance"] / real["significance"] == pytest.approx(
        1.0 / math.sqrt(lg.COHERENT_JPSI_EFF_IR8_LI7), rel=1e-12)


def test_the_verdict_the_statistics_are_in_rho_and_phi_not_jpsi(o5):
    rows = {(r["vm"], r["window"]): r
            for r in o5.channel_rows(eps_det=o5.EPS_DET_CHANG)}
    # rho and phi clear 3 sigma even inside LiPolGen's own Q^2 > 0.7 window
    assert rows[("rho", "Q2>0.7")]["significance"] > 3.0
    assert rows[("rho", "Q2>0.7")]["significance"] == pytest.approx(8.615,
                                                                    abs=0.02)
    assert rows[("phi", "Q2>0.7")]["significance"] == pytest.approx(2.091,
                                                                    abs=0.02)
    # those are on the background-immune <P_zz^2> = 0.81; on the optimal 0.90
    # they are the 9.08 / 2.20 the first write-up quoted, 1.0541x larger
    opt = {(r["vm"], r["window"]): r
           for r in o5.channel_rows(eps_det=o5.EPS_DET_CHANG,
                                    pzz2=o5.PZZ2_OPTIMAL)}
    assert opt[("rho", "Q2>0.7")]["significance"] == pytest.approx(9.08,
                                                                   abs=0.02)
    assert opt[("phi", "Q2>0.7")]["significance"] == pytest.approx(2.20,
                                                                   abs=0.02)
    # ... but BOTH sit below the coherent channel's own M_X floor, which is
    # why this is not a result the shipped generator can consume.
    assert lg.COHERENT_MX_MIN_DEFAULT > 1.019      # phi
    assert lg.COHERENT_MX_MIN_DEFAULT > 0.775      # rho


def test_the_report_runs_and_self_checks(o5):
    buf = io.StringIO()
    res = o5.report(out=buf)
    o5._self_check(res)                     # raises if any headline moved
    text = buf.getvalue()
    assert "THE DECIDING NUMBER" in text
    assert "near-null is expensive" in text
    # the limitation must ride along with the number, in the report itself
    assert "NOT a dipole-model calculation" in o5.__doc__
    # and so must the two things that are NOT established
    assert "UNMEASURED anywhere in this tree" in text
    assert "lumi_fraction" in text
    assert "MARGINAL, A BAND RATHER THAN A POINT" in text
    # the band's TOP is a SPAN, and the report has to say the crossing is not
    # established rather than print a single top
    assert "BAND'S TOP IS ITSELF A SPAN" in text
    assert "WHETHER THE BAND'S TOP CROSSES 3 SIGMA IS NOT" in text
    assert "IT STRADDLES 1" in text
    # the retracted non sequitur, and the withdrawn instruction
    assert "THAT IS A NON SEQUITUR" in text
    assert "Read the low end" not in text
    assert "species lever\n    on its own" not in text
    # ... and the two omissions of the third pass, with their directions
    assert "THE BEAM ENERGY.  UP." in text
    assert "THE DECAY LEPTONS.  DOWN" in text
    assert "THEY PARTLY CANCEL" in text
    assert "UNBOUNDED" in text
    # the self-contradiction that used to sit three lines under the optics
    # row must be GONE, in the corrected form, not merely softened
    assert "Every single pessimism on its own leaves 3 sigma" not in text
    assert "EXCEPT ONE leaves 3 sigma" in text
