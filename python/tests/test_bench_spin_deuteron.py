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

VENDORED FILES.  Every harness here re-checks a recorded sha256 of each file
it reads (since 2026-09-26); each file has a guard below -- its count, first
and last values, a byte-identical copy accepted and ONE flipped byte refused
with the sha256 error (so the harness exits 2, never 0 on moved numbers).
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


def _one_flipped_byte_is_refused(tmp_path, loader, src, old, new):
    """A byte-identical copy of `src` loads (the check is on content, not on
    the path); the same copy with `old` -> `new` (ONE byte, at ONE place)
    raises the sha256 RuntimeError."""
    raw = open(src, "rb").read()
    assert raw.count(old) == 1, old
    assert len(old) == len(new) and sum(a != b for a, b in zip(old, new)) == 1
    good = tmp_path / ("good_" + os.path.basename(src))
    good.write_bytes(raw)
    loader(str(good))
    bad = tmp_path / ("bad_" + os.path.basename(src))
    bad.write_bytes(raw.replace(old, new))
    with pytest.raises(RuntimeError, match="sha256"):
        loader(str(bad))


# ------------------------------------------------------------------ row 1

@pytest.fixture(scope="module")
def epios_rep():
    return epios.run(verbose=False)


def test_epios_report_row_passes(epios_rep):
    _report_row(epios_rep)
    assert epios_rep["status"] == "pass"
    assert len(epios_rep["rows"]) == 8
    assert epios_rep["worst_residual"] <= epios.TOL_ROUNDTRIP


def test_epios_table_is_the_vendored_one(tmp_path):
    rows = epios.load_table()
    assert [r["mode"] for r in rows] == list(range(8))
    assert (rows[0]["pz"], rows[0]["pzz"], rows[0]["pz_lep"]) == (0, 0, 0.0)
    assert (rows[-1]["pz"], rows[-1]["pzz"], rows[-1]["pz_lep"]) == (
        Fraction(-1, 2), Fraction(-1, 2), -0.417)
    # the round trip PASSES for any physical table, so this is what stops an
    # altered one from passing: it is refused (the harness exits 2)
    _one_flipped_byte_is_refused(tmp_path, epios.load_table, epios.DATA,
                                 b'"pz_lep": -0.572', b'"pz_lep": -0.573')


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
    each mode.  The tree's ladder reproduces the ideal P_zz only at mode 0 of
    the six modes it computes; at modes 5, 6 (|P_z| = 1) it throws, and only
    the closed-form LIMIT equals the table there (counted as matches, "3 of
    8", until 2026-09-26).  EST is >= 0 everywhere, so it cannot give modes
    3, 4, 7 (P_zz < 0) at any P_z."""
    rows = {r["mode"]: r for r in epios_rep["rows"]}
    assert epios_rep["est_modes_matched"] == [0]
    assert epios_rep["est_modes_computable"] == [0, 1, 2, 3, 4, 7]
    assert epios_rep["est_modes_tree_throws"] == [5, 6]
    assert epios_rep["est_modes_closed_form_limit_matches"] == [5, 6]
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


def test_est_sources_file_is_the_vendored_one(tmp_path):
    d = est.load_sources()
    fig = d["koivuniemi2004"]["fig1_downstream"]
    assert (fig["P_percent"], fig["p_aligned_percent"], fig["p0_percent"],
            fig["p_anti_percent"]) == (53.0, 63.6, 25.9, 10.6)
    assert d["koivuniemi2004"]["T_half_convention_percent"] == 11
    _one_flipped_byte_is_refused(tmp_path, est.load_sources, est.DATA,
                                 b'"P_percent": 53.0', b'"P_percent": 53.1')


def test_est_compass_6lid_numbers(est_rep):
    """Koivuniemi SPIN 2004 Fig. 1, from the tree's ladder at P_d = 0.530 and
    a common spin temperature.  Tolerances are the printed digits (inputs
    and outputs both rounded)."""
    c = est_rep["compass"]
    assert c["pops_maxdiff"] <= 0.05              # printed to 0.1 %
    assert c["t_half"] == pytest.approx(c["t_printed"], abs=0.005)
    # the Cartesian P_zz at P_d = 0.530 IS the EST closed form, to the row's
    # own tolerance (the maxent bisection bound; measured |diff| ~ 3e-15)
    assert c["pzz_cartesian"] == pytest.approx(est.closed_form(0.530), abs=est.TOL)
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


def test_hermes_table_vendored_as_published(tmp_path):
    rows = hermes.load_table()
    assert [r["x"] for r in rows] == [0.012, 0.032, 0.063, 0.128, 0.248, 0.452]
    assert rows[0]["b1"] == pytest.approx(0.1120, abs=1e-12)
    assert rows[5]["azz"] == pytest.approx(0.0157, abs=1e-12)
    assert rows[2]["b1_err"] == pytest.approx(math.hypot(1.11, 0.60) * 1e-2, rel=1e-12)
    _one_flipped_byte_is_refused(tmp_path, hermes.load_table, hermes.DATA,
                                 b"11.20, 5.51", b"11.21, 5.51")


def test_r1990_file_vendored_with_the_typo_corrected(tmp_path):
    _one_flipped_byte_is_refused(tmp_path, hermes.r1990_factory,
                                 hermes.R1990_DATA,
                                 b'"b1": 0.0635', b'"b1": 0.0636')


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
    """chi2_b1 does not depend on the F2 backend and is pinned for every
    configuration that runs.  chi2_azz divides by F1^d, rebuilt from the F2
    backend -- MstwSF, or ToyF2 when MSTW is absent (no PYTHIA 8 tier, or no
    pdfdata/mstw2008lo.00.dat) -- and its pins were measured with MstwSF, so
    they are asserted only then.  Without MSTW the test asserts what it can,
    then SKIPS naming what it could not check."""
    mstw = hermes_rep["f2_source"] == "MstwSF"
    unchecked = []
    for name, (cb, ca) in HERMES_CHI2.items():
        c = hermes_rep["configs"][name]
        if c["status"] == "blocked":   # cdks_conv_mstw, and only when MSTW
            assert name == "cdks_conv_mstw" and not mstw, name   # is absent
            unchecked.append("%s (blocked: %s)" % (name, c["reason"]))
            continue
        assert c["ndf"] == 6
        assert c["chi2_b1"] == pytest.approx(cb, rel=1e-3), name
        if mstw:
            assert c["chi2_azz"] == pytest.approx(ca, rel=1e-3), name
    if not mstw:
        unchecked.append("every chi2_azz pin (measured with MstwSF; this "
                         "run's F1^d uses %s)" % hermes_rep["f2_source"])
        pytest.skip("chi2_b1 pins checked; NOT checked: " + "; ".join(unchecked))


def test_hermes_chi2_b1_does_not_depend_on_the_f2_backend(hermes_rep, monkeypatch):
    """The no-MSTW path, forced: cdks_conv_mstw is blocked by name, F1^d
    falls back to ToyF2, every chi2_b1 is unchanged (only A_zz divides by
    F1^d), and b1 = 0 keeps its A_zz chi2 whatever F1^d is."""
    if hermes_rep["f2_source"] != "MstwSF":
        pytest.skip("needs the MstwSF run to compare with (no MSTW here)")
    monkeypatch.setattr(hermes, "have_mstw", lambda: (False, "forced by the test"))
    toy = hermes.run(verbose=False)
    assert toy["f2_source"] == "ToyF2"
    assert toy["headline"] == "cdks_conv_default"
    assert toy["configs"]["cdks_conv_mstw"] == dict(status="blocked",
                                                    reason="forced by the test")
    for name in hermes.CONFIGS[1:]:
        assert toy["configs"][name]["chi2_b1"] == pytest.approx(
            hermes_rep["configs"][name]["chi2_b1"], rel=1e-12), name
    assert toy["configs"]["null_b1_zero"]["chi2_azz"] == pytest.approx(
        hermes_rep["configs"]["null_b1_zero"]["chi2_azz"], rel=1e-12)


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

def test_bonus_files_are_the_vendored_ones(tmp_path):
    q = bonus.load_query()
    eids = [e["eid"] for e in q["experiments"]]
    assert (len(eids), eids[0], eids[-1]) == (19, 12, 187)
    assert sum(e["n_measurements"] for e in q["experiments"]) == 918
    _one_flipped_byte_is_refused(tmp_path, bonus.load_query, bonus.QUERY,
                                 b'"eid": 135', b'"eid": 136')
    d = bonus.load_deeps()
    blocks = d["blocks"]
    assert [b["mid"] for b in blocks] == list(range(1, 116))
    assert sum(len(b["rows"]) for b in blocks) == 1198
    assert blocks[0]["rows"][0] == [-0.85, 0.0, 0.0, 0.000281081]
    assert blocks[-1]["rows"][-1] == [0.25, 0.00156769, 8.39487e-05, 0.000477365]
    _one_flipped_byte_is_refused(tmp_path, bonus.load_deeps, bonus.DEEPS,
                                 b"0.00156769", b"0.00156768")


def test_bonus_deeps_census_matches_the_docstring():
    """What the vendored Deeps file holds (the harness docstring): 60 blocks
    at the paper's Q2 bins, 55 with unverified labels, 25 duplicated label
    triples, one empty block, 60 value = stat = 0 rows."""
    f = bonus.findings()
    assert f["deeps_blocks"] == 115
    assert f["deeps_blocks_paper_q2"] == 60
    assert f["deeps_blocks_unverified_q2"] == 55
    assert f["deeps_duplicate_label_triples"] == 25
    assert f["deeps_empty_blocks"] == [78]
    assert f["deeps_zero_value_zero_stat_rows"] == 60


def test_bonus_harness_prints_its_row_on_a_non_utf8_stdout():
    """The catalogue's quantity names carry Greek letters; a Latin-1 or ASCII
    stdout gets them escaped, and the REPORT row and exit 0 regardless."""
    import subprocess
    import sys
    for enc in ("latin-1", "ascii"):
        env = dict(os.environ, PYTHONIOENCODING=enc)
        p = subprocess.run([sys.executable, os.path.join(
            _BENCH, "t2_bonus_spectator_shape.py")], env=env,
            capture_output=True)
        assert p.returncode == 0, (enc, p.stderr.decode(enc, "replace"))
        assert b"REPORT | t2_bonus_spectator_shape |" in p.stdout, enc
        assert b"\\u03b7" in p.stdout, enc      # eta, escaped


def test_bonus_row_is_blocked_by_name():
    rep = bonus.run(verbose=False)
    _report_row(rep)
    assert rep["status"] == "blocked"
    f = rep["findings"]
    assert f["bonus"]["quantities"] == ["F2n/F2p"]
    assert f["bonus"]["n_measurements"] == 1
    assert f["deeps_blocks"] == 115
    pytest.skip(bonus.BLOCKED_REASON)
