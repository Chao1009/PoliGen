# SPDX-License-Identifier: GPL-3.0-or-later
"""The T5 MC-closure harnesses under validation/benchmarks/, exercised as
tests (05_chain.md sec. 1.1-1.2, CH-25):

  t5_sum_rules.py            Burkhardt-Cottingham on g2_ww (implementation
                             closure, PASS); Bjorken on ToyG1 (recorded
                             FAIL); Bjorken on LhapdfG1 (BLOCKED unless
                             NNPDFpol11_100 is installed); FSI unitarity
                             (BLOCKED: not definable in the tree as built)
  t5_weighted_unweighted.py  one physics point (6Li inclusive, x >= 0.1,
                             y >= 0.5, helicity-flip apar+) unweighted and in
                             Mode W, the two sigma estimates, |pull| <= 3

Both are cheap (~1.6 s and ~2.7 s), so neither sits behind LIPOLGEN_BENCH.
A harness that DISAGREES with the tree is not a failing test: the tests pin
what the harnesses MEASURE and that each verdict follows from the numbers.
No test here moves, or asserts a new value for, any tree default.
"""
import hashlib
import importlib.util
import json
import math
import os
import subprocess
import sys

import pytest

lg = pytest.importorskip("lipolgen")
np = pytest.importorskip("numpy")
# numpy < 2.0 spells it np.trapz; >= 2.0 np.trapezoid (2.4 removed np.trapz),
# as in validation/vmc_reconcile.py; pyproject.toml declares numpy >= 1.20
_trapz = getattr(np, "trapezoid", None) or np.trapz

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_BENCH = os.path.join(_ROOT, "validation", "benchmarks")


def _load(name):
    path = os.path.join(_BENCH, name + ".py")
    spec = importlib.util.spec_from_file_location("bench_" + name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _needs_scipy():
    pytest.importorskip(
        "scipy.integrate",
        reason="environment: scipy not installed (t5_sum_rules integrates "
               "with scipy.integrate.quad)")


# ============================================================ t5_sum_rules

@pytest.fixture(scope="module")
def sr_mod():
    return _load("t5_sum_rules")


@pytest.fixture(scope="module")
def sr_rep(sr_mod):
    _needs_scipy()
    return sr_mod.run(verbose=False)


def test_sr_reference_constants_are_the_recorded_ones(sr_mod):
    """The in-file stand-in for a vendored table: 7 constants, first and last
    by sorted key, and one moved digit refused with the sha256 error."""
    ref = sr_mod.load_reference()
    keys = sorted(ref)
    assert len(ref) == 7
    assert keys[0] == "alpha_s_mz" and ref[keys[0]] == 0.1180
    assert keys[-1] == "m_z_gev" and ref[keys[-1]] == 91.1879
    # PDG 2026 (pdg-2026.0 wheel): S017AV, Q005M, Q004M
    assert ref["g_A"] == 1.2753 and ref["g_A_err"] == 0.0013
    assert (ref["m_b_gev"], ref["m_c_gev"]) == (4.186, 1.2729)
    assert (hashlib.sha256(json.dumps(sr_mod.REFERENCE, sort_keys=True)
                           .encode("utf-8")).hexdigest()
            == sr_mod.REFERENCE_SHA256)
    bad = dict(ref)
    bad["g_A"] = 1.2754
    with pytest.raises(RuntimeError, match="sha256"):
        sr_mod.load_reference(bad)


def test_sr_harness_exit_codes(tmp_path):
    """0 with its five REPORT rows printed (FAIL recorded and BLOCKED
    included); 2 when a reference constant was altered in the file."""
    _needs_scipy()
    src = os.path.join(_BENCH, "t5_sum_rules.py")
    ok = subprocess.run([sys.executable, src], capture_output=True, text=True,
                        timeout=300)
    assert ok.returncode == 0, ok.stderr
    rows = [l for l in ok.stdout.splitlines() if l.startswith("REPORT | ")]
    assert [l.split(" | ")[1] for l in rows] == [
        "t5_sum_rules[bc_ww]", "t5_sum_rules[bjorken[toyg1]]",
        "t5_sum_rules[bjorken[nnpdfpol11]]", "t5_sum_rules[fsi_unitarity]",
        "t5_sum_rules"]
    assert rows[0].endswith("| PASS") and rows[1].endswith("| FAIL")
    assert "BLOCKED: not definable in the tree as built" in rows[3]
    assert rows[-1].endswith("| FAIL")
    assert "np.float64" not in ok.stdout
    text = open(src, encoding="utf-8").read()
    assert text.count('"g_A": 1.2753,') == 1
    bad = tmp_path / "t5_sum_rules.py"
    bad.write_text(text.replace('"g_A": 1.2753,', '"g_A": 1.2754,'),
                   encoding="utf-8")
    broken = subprocess.run([sys.executable, str(bad)], capture_output=True,
                            text=True, timeout=300)
    assert broken.returncode == 2
    assert "sha256" in broken.stderr


def test_sr_rule3_deuteron_is_per_nucleus(sr_mod):
    """BENCHMARK_PLAN sec. 8 rule 3 as the header states it: the deuteron's
    g1 is PER NUCLEUS, Z P_p g1p + N P_n g1n with P = 1 - 1.5 P_D = 0.9325,
    not divided by A; the WW transform is linear, so g2_d is the same
    combination of the nucleon g2; the binding's g2_ww IS PolSF.g2p."""
    toy, d = lg.ToyG1(), lg.deuteron()
    assert (d.Z, d.N, d.A) == (1, 1, 2)
    assert d.eff_pol_p == d.eff_pol_n == pytest.approx(1.0 - 1.5 * 0.045, abs=1e-15)
    for x, q2 in ((1e-5, 2.0), (1e-3, 5.0), (0.05, 5.0), (0.3, 10.0), (0.8, 5.0)):
        g1d = toy.g1_nucleus(d, x, q2)
        assert g1d == pytest.approx(0.9325 * (toy.g1p(x, q2) + toy.g1n(x, q2)),
                                    rel=1e-14)
        g2d = lg.g2_ww(lambda xx, qq: toy.g1_nucleus(d, xx, qq), x, q2)
        assert g2d == pytest.approx(0.9325 * (toy.g2p(x, q2) + toy.g2n(x, q2)),
                                    rel=1e-12, abs=1e-15)
        assert lg.g2_ww(toy.g1p, x, q2) == toy.g2p(x, q2)
        assert lg.g2_ww(toy.g1n, x, q2) == toy.g2n(x, q2)


def test_sr_bc_closure_verdict_follows_the_numbers(sr_rep):
    bc = sr_rep["subrows"]["bc_ww"]
    cells = [r for v in bc["cells"].values() for r in v["rows"]]
    assert bc["n_cells"] == len(cells) == 45
    for r in cells:
        assert r["tol"] == pytest.approx(r["bound"] + r["quad_err"], rel=1e-15)
        assert r["ok"] == (abs(r["residual"]) <= r["tol"])
        # the residual IS the tree's 96-point trapezoid error: Euler-Maclaurin
        assert r["residual"] / r["em_prediction"] == pytest.approx(1.0, abs=1e-3)
    assert bc["status"] == ("pass" if all(r["ok"] for r in cells) else "fail")
    assert bc["status"] == "pass"
    assert bc["max_ratio"] == pytest.approx(0.667, abs=2e-3)
    p5 = bc["cells"][("p", 5.0)]["rows"][-1]
    assert p5["x_min"] == 1e-6
    assert p5["integral"] == pytest.approx(-2.035413e-04, rel=1e-6)
    assert p5["truncation"] == pytest.approx(-2.111946e-04, rel=1e-6)
    n5 = bc["cells"][("n", 5.0)]["rows"][-1]
    assert n5["truncation"] == pytest.approx(9.946438e-02, rel=1e-6)
    assert n5["residual"] == pytest.approx(-6.517e-04, rel=1e-3)


def test_sr_bc_closure_is_linear_in_g1(sr_rep):
    """Rule 3's 'the normalisation cancels': T_d = 0.9325 (T_p + T_n)."""
    c = sr_rep["subrows"]["bc_ww"]["cells"]
    for q2 in (2.0, 5.0, 10.0):
        for rp, rn, rd in zip(c[("p", q2)]["rows"], c[("n", q2)]["rows"],
                              c[("d", q2)]["rows"]):
            assert rd["truncation"] == pytest.approx(
                0.9325 * (rp["truncation"] + rn["truncation"]), rel=1e-8)


def test_sr_bc_closure_is_not_tuned_to_96_points(sr_mod):
    """The Peano bound follows from npts alone: at npts = 8 the residual is
    ~(95/7)^2 larger and still inside its bound, at ~2/3 of it."""
    _needs_scipy()
    toy = lg.ToyG1()
    g1 = lambda x: toy.g1p(x, 5.0)                                # noqa: E731
    r8 = sr_mod.bc_closure(g1, lambda x: toy.g2p(x, 5.0, 8), npts=8)
    r96 = sr_mod.bc_closure(g1, lambda x: toy.g2p(x, 5.0, 96), npts=96)
    for a, b in zip(r8, r96):
        assert a["ok"] and 0.5 < a["ratio"] < 1.0
        assert 150.0 < a["residual"] / b["residual"] < 220.0


def test_sr_bc_closure_catches_a_wrong_ww(sr_mod):
    """What the closure is FOR: a WW with the Jacobian -ln x dropped, or with
    a g1(u)/u^2 kernel, fails by orders of magnitude."""
    _needs_scipy()
    toy = lg.ToyG1()
    g1 = lambda x: toy.g1p(x, 5.0)                                # noqa: E731
    t = np.linspace(0.0, 1.0, 96)

    def no_jacobian(x):
        return -g1(x) + float(_trapz([g1(x ** (1 - v)) for v in t], t))

    def u2_kernel(x):
        return -g1(x) - math.log(x) * float(_trapz(
            [g1(x ** (1 - v)) / x ** (1 - v) for v in t], t))

    for bad in (no_jacobian, u2_kernel):
        rows = sr_mod.bc_closure(g1, bad, x_mins=(1e-2, 1e-3))
        assert not any(r["ok"] for r in rows)
        assert min(r["ratio"] for r in rows) > 1e3


def test_sr_bc_limit_exists_for_the_proton_only(sr_rep):
    """Recorded, not asserted by the harness: T(x_min) -> 0 for the toy
    proton; for n and d it grows like x_min^-lambda, lambda the toy's own."""
    bc = sr_rep["subrows"]["bc_ww"]
    for q2 in (2.0, 5.0, 10.0):
        lam = 0.045 * math.log(max(q2, 1.1) / 0.04)
        tp = [abs(r["truncation"]) for r in bc["cells"][("p", q2)]["rows"]]
        tn = [r["truncation"] for r in bc["cells"][("n", q2)]["rows"]]
        assert all(b < a for a, b in zip(tp, tp[1:]))
        assert all(b > a > 0 for a, b in zip(tn, tn[1:]))
        assert bc["limit_exponent"]["p"][q2] < -0.4
        assert bc["limit_exponent"]["n"][q2] == pytest.approx(lam, abs=2e-3)
        assert bc["limit_exponent"]["d"][q2] > 0.15


def test_sr_bjorken_reference(sr_mod):
    """C_Bj coefficients (Larin-Vermaseren) and the leading-twist number."""
    _, c2, c3 = sr_mod.c_bjorken_coeffs(3)
    assert (c2, c3) == (pytest.approx(3.5833, abs=1e-4), pytest.approx(20.215, abs=2e-3))
    _, c2, c3 = sr_mod.c_bjorken_coeffs(4)
    assert (c2, c3) == (pytest.approx(3.25, abs=1e-12), pytest.approx(13.850, abs=2e-3))
    ref = sr_mod.load_reference()
    assert sr_mod.alpha_s(ref["m_z_gev"], ref) == pytest.approx(0.1180, abs=1e-12)
    r5 = sr_mod.bjorken_reference(5.0)
    assert r5["nf"] == 4
    assert r5["alpha_s"] == pytest.approx(0.2850, abs=5e-4)
    assert r5["gamma"] == pytest.approx(0.18538, abs=2e-5)
    assert r5["err"] == pytest.approx(0.00233, abs=5e-5)
    gams = [sr_mod.bjorken_reference(q2)["gamma"] for q2 in (2.0, 5.0, 10.0, 100.0)]
    assert all(b > a for a, b in zip(gams, gams[1:]))
    assert gams[-1] < ref["g_A"] / 6.0


def test_sr_bjorken_toy_fails_as_recorded(sr_rep, sr_mod):
    bj = sr_rep["subrows"]["bjorken[toyg1]"]
    for q2, v in bj["per_q2"].items():
        g_ref = v["reference"]["gamma"]
        last, prev = v["rows"][-1], v["rows"][-2]
        step = abs(last["gamma_pn"] - prev["gamma_pn"])
        ok = (step <= sr_mod.BJ_CONV_REL * g_ref
              and abs(last["gamma_pn"] - g_ref) <= sr_mod.BJ_TOL_REL * g_ref)
        assert v["verdict"]["ok"] == ok
        # the proton moment converges, the neutron's does not (g1n ~ x^-1-lambda)
        assert abs(last["gamma_p"] - prev["gamma_p"]) < 1e-3
        assert abs(last["gamma_n"] - prev["gamma_n"]) > 0.1
    assert bj["status"] == ("pass" if all(v["verdict"]["ok"]
                                         for v in bj["per_q2"].values()) else "fail")
    assert bj["status"] == "fail"
    q5 = bj["per_q2"][5.0]
    assert [r["x_min"] for r in q5["rows"]] == [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
    assert q5["rows"][-1]["gamma_pn"] == pytest.approx(0.64302, rel=1e-4)
    assert q5["rows"][0]["gamma_pn"] == pytest.approx(0.13907, rel=1e-4)
    assert q5["rows"][-1]["gamma_p"] == pytest.approx(0.12816, rel=1e-4)
    assert q5["verdict"]["last_decade_step"] == pytest.approx(0.21982, rel=1e-4)
    # the row as a whole: fail if any evaluated sub-row fails
    evaluated = [s["status"] for s in sr_rep["subrows"].values()
                 if s["status"] != "blocked"]
    assert sr_rep["status"] == ("pass" if all(s == "pass" for s in evaluated)
                                else "fail")


def test_sr_bjorken_nnpdfpol_subrow_names_what_is_missing(sr_rep, sr_mod):
    row = sr_rep["subrows"]["bjorken[nnpdfpol11]"]
    lh, why = sr_mod.lhapdf_g1()
    if lh is None:
        assert row["status"] == "blocked"
        assert row["reason"] == why
        assert row["reason"] in (
            "environment: LHAPDF set NNPDFpol11_100 not installed",
            "environment: lipolgen built without the LHAPDF tier "
            "(HAVE_LHAPDF false); needs LHAPDF set NNPDFpol11_100")
    else:
        assert row["status"] in ("pass", "fail")


def test_sr_bjorken_nnpdfpol_config(sr_mod):
    """Runs only where the polarized set is installed; skips naming it."""
    _needs_scipy()
    lh, why = sr_mod.lhapdf_g1()
    if lh is None:
        pytest.skip(why)
    res = sr_mod.measure_bjorken(lh, sr_mod.X_MIN_SCAN_LHAPDF)
    for v in res["per_q2"].values():
        g_ref = v["reference"]["gamma"]
        last, prev = v["rows"][-1], v["rows"][-2]
        assert v["verdict"]["ok"] == (
            abs(last["gamma_pn"] - prev["gamma_pn"]) <= sr_mod.BJ_CONV_REL * g_ref
            and abs(last["gamma_pn"] - g_ref) <= sr_mod.BJ_TOL_REL * g_ref)


def test_sr_fsi_unitarity_is_blocked_not_definable(sr_rep):
    row = sr_rep["subrows"]["fsi_unitarity"]
    assert row["status"] == "blocked"
    assert row["reason"].startswith("not definable in the tree as built")
    ts = row["tree_side"]
    # recorded, compared with nothing: the default is absorptive by design
    assert ts["survival"] == pytest.approx(0.520239, rel=1e-5)
    assert ts["int_w_minus_1"] == pytest.approx(ts["survival"] - 1.0, abs=1e-15)
    assert ts["sigma_el_cluster_mb"] / ts["sigma_tot_cluster_mb"] == pytest.approx(
        0.269, abs=1e-3)
    assert ts["elastic_gain"] == 1.0 and ts["k_max_gev"] == 1.2
    assert ts["w_max"] == 50.0


# =================================================== t5_weighted_unweighted

@pytest.fixture(scope="module")
def wu_mod():
    return _load("t5_weighted_unweighted")


@pytest.fixture(scope="module")
def wu_rep(wu_mod):
    return wu_mod.run(verbose=False, keep_weights=True)


def test_wu_row_passes_and_verdict_follows(wu_rep, wu_mod):
    assert wu_rep["pull"] == pytest.approx(
        (wu_rep["sigma_u"] - wu_rep["sigma_w"])
        / math.hypot(wu_rep["err_u"], wu_rep["err_w"]), rel=1e-12)
    assert wu_rep["status"] == ("pass" if abs(wu_rep["pull"]) <= wu_mod.PULL_MAX
                                else "fail")
    assert wu_rep["status"] == "pass"
    # the estimators are what the header says they are
    assert wu_rep["sigma_u"] == wu_rep["n_unweighted"] / wu_rep["lumi_pb"]
    assert wu_rep["err_u"] == pytest.approx(
        math.sqrt(wu_rep["n_unweighted"]) / wu_rep["lumi_pb"], rel=1e-15)
    w = wu_rep["w_on_weighted"]
    assert w.size == wu_rep["n_weighted"]
    assert wu_rep["sigma_w"] == pytest.approx(float(w.sum()) / wu_rep["lumi_pb"],
                                              rel=1e-15)
    # pinned: fixed seed 20260926, runs 1 / 2, L = 34000 pb^-1, x >= 0.1, y >= 0.5
    assert (wu_rep["seed"], wu_rep["runs"], wu_rep["lumi_pb"]) == (
        20260926, (1, 2), 34000.0)
    assert (wu_rep["x_min"], wu_rep["y_min"]) == (0.1, 0.5)
    assert (wu_rep["n_unweighted"], wu_rep["n_weighted"]) == (3032577, 2987706)
    assert wu_rep["sigma_u"] == pytest.approx(89.193441, rel=1e-8)
    assert wu_rep["sigma_w"] == pytest.approx(89.167357, rel=1e-8)
    assert wu_rep["err_w"] == pytest.approx(0.051588, rel=1e-5)
    assert wu_rep["pull"] == pytest.approx(0.359, abs=2e-3)
    assert wu_rep["sigma_cat_pb"] == pytest.approx(89.167123, rel=1e-8)
    assert wu_rep["sigma_unpol_pb"] == pytest.approx(87.873873, rel=1e-8)


def test_wu_sensitivity_is_as_stated(wu_rep):
    """The window was chosen so that the closure can FAIL: the exact
    polarized shift is +1.47 % and a no-op Mode W would pull ~ +18."""
    assert wu_rep["rate_shift"] == pytest.approx(1.4717e-2, abs=5e-7)
    assert wu_rep["null_pull"] > 10.0
    assert wu_rep["null_pull"] == pytest.approx(18.29, abs=0.02)


def test_wu_compare_catches_mutated_weights(wu_rep, wu_mod):
    """What the closure is FOR: a no-op Mode W, a 1 % normalisation error in
    w, or the polarized part of w with its sign flipped all FAIL."""
    n_u, w, lumi = wu_rep["n_unweighted"], wu_rep["w_on_weighted"], wu_rep["lumi_pb"]
    honest = wu_mod.compare(n_u, w, lumi)
    assert honest["status"] == "pass"
    assert honest["pull"] == pytest.approx(wu_rep["pull"], rel=1e-12)
    for bad in (np.ones_like(w), 1.01 * w, 2.0 - w):
        out = wu_mod.compare(n_u, bad, lumi)
        assert out["status"] == "fail" and abs(out["pull"]) > 5.0


def test_wu_rule3_one_per_nucleon_sampler(wu_mod, wu_rep):
    """Both routes come from ONE per-nucleon sampler: the unpolarized
    category's sigma is the per-nucleon cell sum (vector and tensor parts
    cancel over m), and the category's is sum_m p_m sigma_m -- the two
    numbers the header calls exact."""
    pt = wu_mod.physics_point()
    s, cat = pt["sampler"], pt["category"]
    assert wu_rep["sigma_unpol_pb"] == pytest.approx(float(np.sum(s.cell_xsec_pb)),
                                                     rel=1e-12)
    mix = sum(p * s.sigma_state_pb(cat, m)
              for p, m in zip(cat.populations, (1.0, 0.0, -1.0)))
    assert wu_rep["sigma_cat_pb"] == pytest.approx(mix, rel=1e-12)
    assert pt["beam"] == str(lg.default_configs("6Li")[1])
    assert (cat.lam_e, cat.pe, cat.j) == (1, 0.7, 1.0)
    assert s.scenario.y_min == 0.5 and pt["x_min"] == 0.1
    assert s.n_cells == 550


def test_wu_fixed_seeds_reproduce_and_row_prints(wu_mod, capsys):
    a = wu_mod.run(verbose=False, lumi_pb=1000.0)
    b = wu_mod.run(verbose=True, lumi_pb=1000.0)
    for k in ("n_unweighted", "n_weighted", "sigma_u", "sigma_w", "pull"):
        assert a[k] == b[k]
    assert a["n_unweighted"] > 0
    out = capsys.readouterr().out
    rows = [l for l in out.splitlines() if l.startswith("REPORT | ")]
    assert len(rows) == 1
    assert rows[0].startswith("REPORT | t5_weighted_unweighted | ")
    assert rows[0].endswith("| " + b["status"].upper())
    assert "np.float64" not in out
