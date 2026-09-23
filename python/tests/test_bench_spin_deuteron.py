# SPDX-License-Identifier: GPL-3.0-or-later
"""Benchmark harness group SPIN AND DEUTERON -- BENCHMARK_PLAN.md sec. 4 rows
1, 3, 9 and 10, run as tests.

The harnesses live in `validation/benchmarks/` and are imported here by file
path, never retyped:

  row 1   t5_epios_source_modes.py      EPIOS Table II vs spin1_populations
  row 3   t2_hermes_b1_table2.py        HERMES b1^d Table II, chi2 with no tuning
  row 9   t5_est_identity.py            spin_temperature_pzz == EST closed form
  row 10  t2_bonus_spectator_shape.py   CLAS DB query -- BLOCKED, skips by name

Every row is cheap (the HERMES row, the slowest, is ~5 s), so none sits
behind the `LIPOLGEN_BENCH=1` opt-in.

WHAT IS ASSERTED AND WHAT IS PINNED.  Rows 1 and 9 assert identities (the
spin algebra, the EST closed form).  Row 3 asserts NOTHING about agreement --
a disagreement with HERMES is a measurement, recorded in
docs/open_items/run_2026-09-23/phase_B1_spin_deuteron.md and priced for
registry row 2 -- but it PINS the measured chi2 values (rel. 1e-3), so any
change to a b1 backend, the MSTW reader or the vendored table moves a test
and has to be re-recorded, which is the M2 rule of BENCHMARK_PLAN.md sec. 3.
Row 1's EST distances are pinned the same way.
"""

import importlib.util
import math
import os
from fractions import Fraction

import pytest

import lipolgen._lipolgen as _l

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_BENCH = os.path.join(_ROOT, "validation", "benchmarks")


def _load(name):
    spec = importlib.util.spec_from_file_location(
        "bench_" + name, os.path.join(_BENCH, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


epios = _load("t5_epios_source_modes")
est = _load("t5_est_identity")
hermes = _load("t2_hermes_b1_table2")
bonus = _load("t2_bonus_spectator_shape")


def _report_row(rep):
    for k in ("name", "generator", "reference", "tolerance", "status"):
        assert k in rep, k
    assert rep["status"] in ("pass", "fail", "blocked")


# ------------------------------------------------------------------ row 1

@pytest.fixture(scope="module")
def epios_rep():
    return epios.run(verbose=False)


def test_epios_report_row_passes(epios_rep):
    _report_row(epios_rep)
    assert epios_rep["status"] == "pass"
    assert len(epios_rep["rows"]) == 8
    assert epios_rep["worst_residual"] <= epios.TOL_ROUNDTRIP


def test_epios_populations_are_the_exact_rationals(epios_rep):
    for r in epios_rep["rows"]:
        exact = [float(Fraction(s)) for s in r["pops_exact"]]
        assert max(abs(a - b) for a, b in zip(r["pops"], exact)) <= 1e-15
        assert min(r["pops"]) >= 0.0
        assert abs(sum(r["pops"]) - 1.0) <= 1e-15


def test_epios_every_polarized_mode_is_on_the_domain_edge(epios_rep):
    # measured: min p_m = 0 for modes 1-7; mode 0 (unpolarized) is interior
    for r in epios_rep["rows"]:
        if r["mode"] == 0:
            assert r["min_pop"] == pytest.approx(1.0 / 3.0, abs=1e-15)
        else:
            assert r["min_pop"] == 0.0


def test_epios_ladder_distance_is_recorded_not_asserted(epios_rep):
    """PINS the measured distance of the EST (`--pzz-mode ladder`) fill from
    each mode.  EST reproduces the ideal P_zz only at modes 0, 5, 6 (the last
    two in the |P_z| = 1 limit, where the tree throws); it is >= 0 everywhere
    so it cannot give modes 3, 4, 7 (P_zz < 0) at any P_z."""
    rows = {r["mode"]: r for r in epios_rep["rows"]}
    assert epios_rep["est_modes_matched"] == [0, 5, 6]
    third, half = 2 - math.sqrt(11 / 3), 2 - math.sqrt(13 / 4)  # EST(1/3), EST(1/2)
    want = {1: 2 - math.sqrt(8 / 3),          # 0.367007 vs ideal 0
            2: third - 1.0,                   # -0.914854
            3: third + 1.0,                   # +1.085146
            4: half + 0.5, 7: half + 0.5}     # +0.697224
    for m, d in want.items():
        assert rows[m]["est_minus_ideal"] == pytest.approx(d, abs=1e-9)
    for m in (5, 6):
        assert not rows[m]["est_from_tree"]
        assert "maxent populations need |pz| < 1" in rows[m]["est_note"]
    # the measured/ideal vector ratio of the source: 0.731 .. 0.906
    real = [r["realisation"] for r in epios_rep["rows"] if r["realisation"]]
    assert min(real) == pytest.approx(0.731, abs=5e-4)
    assert max(real) == pytest.approx(0.906, abs=5e-4)


# ------------------------------------------------------------------ row 9

@pytest.fixture(scope="module")
def est_rep():
    return est.run(verbose=False)


def test_est_identity_to_the_bisection_bound(est_rep):
    """spin_temperature_pzz(1, P_z) == 2 - sqrt(4 - 3 P_z^2) over the open
    interval, to the bound the maxent bisection guarantees (2.31e-14).
    MEASURED 2026-09-23: max 1.138e-14 over 19999 points, so the 1e-14 the
    plan asked for does not hold on a dense grid (496 points exceed it)."""
    _report_row(est_rep)
    assert est_rep["status"] == "pass"
    assert est_rep["scan"]["worst"] <= est.TOL
    assert est.TOL == pytest.approx(2.3065e-14, rel=1e-3)


def test_est_five_critic_points_meet_1e14():
    # the five points 06_critic.md sec. 4.1 measured (max 7.7e-15)
    for pz in (0.1, 0.3, 0.5, 0.7, 0.9):
        assert abs(_l.spin_temperature_pzz(1.0, pz) - est.closed_form(pz)) <= 1e-14


def test_est_endpoints_throw(est_rep):
    for pz, msg in est_rep["scan"]["endpoints"].items():
        assert msg.startswith("throws"), (pz, msg)
        assert est.closed_form(pz) == 1.0


def test_est_j32_is_the_brillouin_ladder(est_rep):
    j = est_rep["j32"]
    assert j["brillouin"] <= 1e-14
    assert j["geometric"] <= 1e-14
    assert j["rank2"] <= 1e-13


def test_est_compass_6lid_numbers(est_rep):
    """Koivuniemi SPIN 2004 Fig. 1, from the tree's ladder at P_d = 0.530 and
    a common spin temperature.  Tolerances are the printed digits (inputs
    and outputs both rounded)."""
    c = est_rep["compass"]
    assert c["pops_maxdiff"] <= 0.05              # printed to 0.1 %
    assert c["t_half"] == pytest.approx(c["t_printed"], abs=0.005)
    assert c["pzz_cartesian"] == pytest.approx(2 * c["t_half"], rel=1e-15)
    assert c["t_s_mk"] == pytest.approx(c["t_s_printed"], rel=0.01)
    assert c["p6"] == pytest.approx(c["p6_printed"], abs=1e-3)
    assert c["p7"] == pytest.approx(c["p7_printed"], abs=1e-3)


# ------------------------------------------------------------------ row 3

@pytest.fixture(scope="module")
def hermes_rep():
    return hermes.run(verbose=False)


# chi2 over the 6 bins (b1, A_zz with R1990), measured 2026-09-23
HERMES_CHI2 = {
    "cdks_conv_mstw": (22.264, 22.009),
    "cdks_conv_default": (22.669, 22.392),
    "miller_shipped": (5.262, 5.550),
    "miller_x1": (5.297, 4.514),
    "cdks_digitized": (21.782, 21.559),
    "null_b1_zero": (21.799, 21.586),
}


def test_hermes_report_row(hermes_rep):
    _report_row(hermes_rep)
    assert set(hermes_rep["configs"]) == set(hermes.CONFIGS)
    assert len(hermes_rep["per_bin"]) == 6


def test_hermes_table_vendored_as_published():
    rows = hermes.load_table()
    assert [r["x"] for r in rows] == [0.012, 0.032, 0.063, 0.128, 0.248, 0.452]
    assert rows[0]["b1"] == pytest.approx(0.1120, abs=1e-12)
    assert rows[5]["azz"] == pytest.approx(0.0157, abs=1e-12)
    assert rows[2]["b1_err"] == pytest.approx(math.hypot(1.11, 0.60) * 1e-2, rel=1e-12)


def test_hermes_isoscalar_denominator(hermes_rep):
    """Rule 2: F2^d is (F2p + F2n)/2, recomputed here from the backend."""
    sf = _l.MstwSF() if hermes_rep["f2_source"] == "MstwSF" else _l.ToyF2()
    for p in hermes_rep["per_bin"]:
        assert p["f2d"] == pytest.approx(
            0.5 * (sf.f2p(p["x"], p["q2"]) + sf.f2n(p["x"], p["q2"])), rel=1e-12)


def test_hermes_miller_normalisations_named(hermes_rep):
    """Rule 1: the shipped Miller value keeps the per-deuteron -> per-nucleon
    0.5 (registry row 2); miller_x1 is exactly twice it."""
    assert _l.B1_MILLER_TABLE_TO_PER_NUCLEON == 0.5
    assert _l.B1_CDKS_TABLE_TO_PER_NUCLEON == 1.0
    for p in hermes_rep["per_bin"]:
        assert p["b1"]["miller_x1"] == 2.0 * p["b1"]["miller_shipped"]


def test_hermes_chi2_pinned_no_tuning(hermes_rep):
    for name, (cb, ca) in HERMES_CHI2.items():
        c = hermes_rep["configs"][name]
        if c["status"] == "blocked":   # only cdks_conv_mstw can be
            continue                   # (no PYTHIA 8 tier); the rest run
        assert c["ndf"] == 6
        assert c["chi2_b1"] == pytest.approx(cb, rel=1e-3), name
        assert c["chi2_azz"] == pytest.approx(ca, rel=1e-3), name


def test_hermes_f1_rebuild_priced(hermes_rep):
    """F1^d rebuilt with MSTW + R1990 against the F1^d HERMES's own table
    implies: measured 1.026 .. 1.073 (the ALLM97/NMC -> MSTW difference plus
    3-digit rounding).  R1990 vs R1998 moves F1 by 8.6 % in bin 1 and
    <= 0.7 % elsewhere."""
    if hermes_rep["f2_source"] != "MstwSF":
        pytest.skip("F1 price was measured with MstwSF")
    rat = [p["f1_ratio_to_implied"] for p in hermes_rep["per_bin"]]
    assert min(rat) == pytest.approx(1.0256, abs=5e-4)
    assert max(rat) == pytest.approx(1.0727, abs=5e-4)
    rr = [(1 + p["r1990"]) / (1 + p["r1998"]) for p in hermes_rep["per_bin"]]
    assert rr[0] == pytest.approx(0.9143, abs=5e-4)
    assert max(abs(r - 1) for r in rr[1:]) <= 0.007


def test_r1990_typo_guard():
    r = hermes.r1990_factory()
    assert r(0.3, 3.0) == pytest.approx(0.18698, abs=1e-4)


# ------------------------------------------------------------------ row 10

def test_bonus_row_is_blocked_by_name():
    rep = bonus.run(verbose=False)
    _report_row(rep)
    assert rep["status"] == "blocked"
    f = rep["findings"]
    assert f["bonus"]["quantities"] == ["F2n/F2p"]
    assert f["bonus"]["n_measurements"] == 1
    assert f["deeps_blocks"] == 115
    pytest.skip(bonus.BLOCKED_REASON)
