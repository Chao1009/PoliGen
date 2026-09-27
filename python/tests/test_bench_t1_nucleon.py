# SPDX-License-Identifier: GPL-3.0-or-later
"""Tier T1 -- the NUCLEON harnesses under validation/benchmarks/, as tests.

  t1_nmc_f2d.py            NMC F2d (survey F-1) vs (F2p + F2n)/2
  t1_nmc_f2d_over_f2p.py   NMC F2d/F2p (survey F-2) vs (1 + F2n/F2p)/2
  t1_g1d_world.py          COMPASS / HERMES / E143 g1d (D-2, D-1, D-6; SMC
                           recorded) vs g1_nucleus(deuteron)/2

Their references are vendored from the nnpdf-data 4.1.5 wheel by
`validation/benchmarks/vendor_nnpdf_commondata.py` (the only file that may
need the network; nothing here does).  All three harnesses are cheap
(< 1 s each), so nothing sits behind LIPOLGEN_BENCH.

A harness that DISAGREES with the tree is not a failing test: the tests pin
what each harness MEASURES, check that its verdict follows from those
numbers, and assert the normalisation it claims (BENCHMARK_PLAN.md sec. 8
rule 2).  The shipped ToyF2 / ToyG1 FAIL all three rows -- a recorded
measurement, never tuned; no test here moves, or asserts a new value for,
any tree default.  A configuration whose LHAPDF set is not installed must
be REPORTED as blocked with its reason (tested either way); a pin that needs
the PYTHIA tier skips naming exactly that.
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

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_BENCH = os.path.join(_ROOT, "validation", "benchmarks")
_DATA = os.path.join(_BENCH, "data")


def _load(name):
    path = os.path.join(_BENCH, name + ".py")
    spec = importlib.util.spec_from_file_location("bench_" + name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _sf(chi2, n):
    """p(chi2, n), computed HERE (not by the harness)."""
    from scipy.stats import chi2 as c
    return float(c.sf(chi2, n))


# ------------------------------------------------------------------ fixtures

@pytest.fixture(scope="module")
def f2d_mod():
    return _load("t1_nmc_f2d")


@pytest.fixture(scope="module")
def dp_mod():
    return _load("t1_nmc_f2d_over_f2p")


@pytest.fixture(scope="module")
def g1_mod():
    return _load("t1_g1d_world")


@pytest.fixture(scope="module")
def f2d(f2d_mod):
    return f2d_mod.run(verbose=False)


@pytest.fixture(scope="module")
def dp(dp_mod):
    return dp_mod.run(verbose=False)


@pytest.fixture(scope="module")
def g1(g1_mod):
    return g1_mod.run(verbose=False)


def _config(rep, name):
    """A sub-row that must have been evaluated -- skipping, with the
    harness's own reason (which names what is missing), if it is blocked."""
    c = rep["configs"][name]
    if c["status"] == "blocked":
        pytest.skip("%s[%s] blocked: %s" % (rep["name"], name, c["reason"]))
    return c


def _mstw():
    """MstwSF, or a skip naming the missing PYTHIA tier / grid file."""
    import lipolgen._lipolgen as _l
    if not getattr(_l, "HAVE_PYTHIA8", False):
        pytest.skip("environment: lipolgen built without the PYTHIA 8 tier "
                    "(MstwSF)")
    d = _l.pythia8_pdfdata_dir()
    if not (d and os.path.isfile(os.path.join(d, "mstw2008lo.00.dat"))):
        pytest.skip("environment: PYTHIA 8 pdfdata/mstw2008lo.00.dat not "
                    "found (MstwSF)")
    return lg.MstwSF()


# ------------------------------------------------- the vendored files

# (file, harness, loader args, npoints, first (x, Q2, value), last (...))
_FILES = [
    ("nnpdf_NMC_NC_NOTFIXED_D.json", "t1_nmc_f2d", (), 158,
     (0.0045, 0.75, 0.2576), (0.5, 65.0, 0.0729)),
    ("nnpdf_NMC_NC_NOTFIXED.json", "t1_nmc_f2d_over_f2p", (), 260,
     (0.0015, 0.16, 0.9815), (0.675, 99.03, 0.7724)),
    ("nnpdf_COMPASS15_NC_NOTFIXED_MUD.json", "t1_g1d_world", ("compass",), 15,
     (0.0046, 1.1, -0.13), (0.567, 60.8, 0.0203)),
    ("nnpdf_HERMES_NC_7GEV_ED.json", "t1_g1d_world", ("hermes",), 15,
     (0.0264, 1.12, 0.114), (0.7248, 12.21, 0.0116)),
    ("nnpdf_E143_NC_NOTFIXED_ED.json", "t1_g1d_world", ("e143",), 28,
     (0.031, 1.27, 0.15), (0.749, 9.52, 0.009)),
    ("nnpdf_SMC_NC_NOTFIXED_MUD.json", "t1_g1d_world", ("smc",), 13,
     (0.002, 0.5, -0.3), (0.479, 54.8, 0.03)),
]
_WHEEL_SHA256 = "6ed209add4d162e3b9e4801c7b538a28867e479b6ffd4781e774526a3490269a"


def _flip_first_value(raw):
    """The file's bytes with ONE byte changed: the last digit of the first
    point's central value (still valid JSON, a different number)."""
    start = raw.index(b'"value": ', raw.index(b'"points": ['))
    end = raw.index(b",", start) - 1
    assert raw[end:end + 1].isdigit()
    new = b"7" if raw[end:end + 1] != b"7" else b"8"
    out = raw[:end] + new + raw[end + 1:]
    assert json.loads(out)["body"]["points"][0]["value"] != \
        json.loads(raw)["body"]["points"][0]["value"]
    return out


@pytest.mark.parametrize("fname,harness,args,n,first,last", _FILES,
                         ids=[f[0] for f in _FILES])
def test_vendored_file_is_the_wheel_and_refuses_a_flipped_byte(
        fname, harness, args, n, first, last, tmp_path):
    mod = _load(harness)
    pts = mod.load_table(*args)
    assert len(pts) == n
    assert (pts[0]["x"], pts[0]["q2"], pts[0]["value"]) == first
    assert (pts[-1]["x"], pts[-1]["q2"], pts[-1]["value"]) == last
    path = os.path.join(_DATA, fname)
    raw = open(path, "rb").read()
    doc = json.loads(raw)
    prov = doc["_provenance"]
    # the provenance the conventions require, on the file itself
    assert prov["source"]["wheel_sha256"] == _WHEEL_SHA256
    assert prov["source"]["fetched"] == "2026-09-26"
    assert prov["source"]["wheel_url"].startswith(
        "https://files.pythonhosted.org/")
    assert all(p.startswith("nnpdf_data/commondata/" + doc["body"]["setname"])
               for p in prov["source"]["files_read"])
    assert prov["record"]["hepdata"].startswith("https://www.hepdata.net/record/ins")
    assert "CC0" in prov["licence"]["data"]
    assert "GPL-3.0-or-later" in prov["licence"]["package"]
    assert prov["publishing_convention"]["normalisation"]
    assert doc["body"]["npoints"] == n == len(doc["body"]["points"])
    # one flipped byte of data is refused, not silently read
    bad = tmp_path / fname
    bad.write_bytes(_flip_first_value(raw))
    with pytest.raises(RuntimeError, match="sha256"):
        mod.load_table(*args, path=str(bad))


def test_the_body_digest_is_checked_on_its_own(f2d_mod, tmp_path, monkeypatch):
    """Even with the whole-file digest re-recorded for an altered file, the
    canonical-body digest (recorded in the file AND the harness) refuses it."""
    raw = _flip_first_value(open(f2d_mod.DATA, "rb").read())
    bad = tmp_path / "f2d.json"
    bad.write_bytes(raw)
    monkeypatch.setattr(f2d_mod, "FILE_SHA256", hashlib.sha256(raw).hexdigest())
    with pytest.raises(RuntimeError, match="body sha256"):
        f2d_mod.load_table(path=str(bad))


def test_g1d_d_state_label_is_unverified_on_every_file():
    """06_critic.md sec. 5.1's D-state trap: the files must not claim a
    convention the metadata does not state."""
    for fname, harness, *_ in _FILES:
        if harness != "t1_g1d_world":
            continue
        conv = json.load(open(os.path.join(_DATA, fname)))[
            "_provenance"]["publishing_convention"]
        assert conv["d_state"].startswith("UNVERIFIED")
        assert "PER NUCLEON" in conv["normalisation"]


def test_vendor_script_refuses_a_wheel_with_another_digest(tmp_path):
    vend = _load("vendor_nnpdf_commondata")
    fake = tmp_path / "nnpdf_data-4.1.5-py3-none-any.whl"
    fake.write_bytes(b"not the wheel")
    with pytest.raises(SystemExit) as exc:
        vend.main(["--wheel", str(fake), "--out-dir", str(tmp_path)])
    assert exc.value.code == 2
    assert not any(p.suffix == ".json" for p in tmp_path.iterdir())


def test_vendor_script_reproduces_the_files_byte_for_byte():
    """Regenerating from the pinned wheel gives the committed files exactly.
    Needs the wheel on disk (the tests never download it)."""
    wheel = os.environ.get("LIPOLGEN_NNPDF_WHEEL", "")
    if not os.path.isfile(wheel):
        pytest.skip("needs nnpdf_data-4.1.5-py3-none-any.whl (sha256 %s, "
                    "https://pypi.org/project/nnpdf-data/4.1.5/) at "
                    "$LIPOLGEN_NNPDF_WHEEL" % _WHEEL_SHA256)
    pytest.importorskip("yaml", reason="PyYAML (vendoring-time only)")
    vend = _load("vendor_nnpdf_commondata")
    assert vend.main(["--check", "--wheel", wheel]) == 0


# ------------------------------------------------- rule 2, asserted

_PTS = ((0.0046, 1.1), (0.0125, 1.25), (0.1, 5.0), (0.35, 23.0), (0.65, 39.0))


def test_f2d_theory_is_the_isoscalar_nucleon_and_the_kernels(f2d_mod, f2d):
    """BENCHMARK_PLAN.md sec. 8 rule 2: the tree side is (F2p + F2n)/2 of
    the free nucleon -- and it IS the generator's per-nucleon deuteron F2."""
    toy = lg.ToyF2()
    kern = lg.default_inclusive_kernel(lg.deuteron())
    for x, q2 in _PTS:
        iso = 0.5 * (toy.f2p(x, q2) + toy.f2n(x, q2))
        assert f2d_mod.isoscalar_nucleon(toy, x, q2) == iso
        assert kern.tables(x, q2).f2 == iso
    # it is NOT the free proton: at x = 0.5 the two differ by 23 %, three
    # times NMC's 8 % error there
    assert toy.f2p(0.5, 20.0) / f2d_mod.isoscalar_nucleon(toy, 0.5, 20.0) \
        == pytest.approx(1.0 / 0.8125, rel=1e-12)
    assert f2d["configs"]["toy"]["kernel_max_rel_diff"] == 0.0


def test_f2d_theory_mstw_is_the_isoscalar_nucleon_and_the_kernels(f2d_mod, f2d):
    m = _mstw()
    c = _config(f2d, "mstw")
    assert c["kernel_max_rel_diff"] == 0.0
    for x, q2 in _PTS[1:]:
        assert f2d_mod.isoscalar_nucleon(m, x, q2) == \
            0.5 * (m.f2p(x, q2) + m.f2n(x, q2))


def test_dp_theory_is_isoscalar_over_proton_and_prices_the_frozen_hook(dp_mod, dp):
    toy = lg.ToyF2()
    for x, q2 in _PTS:
        assert dp_mod.ratio_q2(toy, x, q2) == pytest.approx(
            0.5 * (toy.f2p(x, q2) + toy.f2n(x, q2)) / toy.f2p(x, q2), rel=1e-15)
        # ToyF2's F2n IS F2p times the hook: the two readings coincide
        assert dp_mod.ratio_frozen(toy, x) == pytest.approx(
            dp_mod.ratio_q2(toy, x, q2), abs=1e-15)
    assert dp["configs"]["toy"]["max_abs_frozen_minus_q2"] < 1e-15
    assert dp["configs"]["toy"]["kernel_max_rel_diff"] == 0.0


def test_dp_mstw_hook_is_q2_frozen(dp_mod, dp):
    m = _mstw()
    for x in (0.05, 0.3, 0.6):
        # the hook IS the Q2-dependent ratio frozen at Q2 = 10 (mstw_sf.cpp);
        # away from Q2 = 10 the two differ
        assert dp_mod.ratio_frozen(m, x) == pytest.approx(
            dp_mod.ratio_q2(m, x, 10.0), rel=1e-14)
        assert dp_mod.ratio_q2(m, x, 2.0) != dp_mod.ratio_q2(m, x, 50.0)
    c = _config(dp, "mstw")
    # MstwSF's hook is Q2-frozen: it differs from the Q2-dependent ratio
    assert c["max_abs_frozen_minus_q2"] == pytest.approx(0.0158072, rel=1e-5)
    assert c["kernel_max_rel_diff"] == 0.0


def test_g1d_theory_is_g1_nucleus_over_two_with_the_d_state_factor(g1_mod, g1):
    """The tree's per-nucleon deuteron g1 is g1_nucleus(deuteron)/2 =
    (1 - 1.5 P_D) (g1p + g1n)/2 with P_D = 0.045 -- asserted, together with
    the no-D-state variant it is reported beside."""
    d = lg.deuteron()
    assert (d.A, d.Z) == (2, 1)
    assert lg.P_D_DEUTERON == 0.045
    assert lg.DEUTERON_VECTOR_POLARIZATION == lg.vector_dilution_of(lg.P_D_DEUTERON)
    assert lg.DEUTERON_VECTOR_POLARIZATION == pytest.approx(1 - 1.5 * 0.045, abs=1e-15)
    assert d.eff_pol_p == d.eff_pol_n == lg.DEUTERON_VECTOR_POLARIZATION
    toy = lg.ToyG1()
    kern = lg.default_inclusive_kernel(d)
    for x, q2 in _PTS:
        iso = 0.5 * (toy.g1p(x, q2) + toy.g1n(x, q2))
        tree = g1_mod.g1d_tree(toy, d, x, q2)
        assert tree == pytest.approx(0.9325 * iso, rel=1e-14)
        assert tree == kern.tables(x, q2).g1
        assert g1_mod.g1d_no_dstate(toy, x, q2) == iso
    assert g1["dstate"]["p_d"] == 0.045
    assert g1["dstate"]["factor"] == pytest.approx(0.9325, abs=1e-15)
    assert g1["dstate"]["vector_dilution_of_p_d"] == g1["dstate"]["factor"]
    assert g1["configs"]["toy"]["kernel_max_rel_diff"] == 0.0
    # the no-D-state variant is exactly the tree theory / 0.9325, per point
    e = g1["configs"]["toy"]["experiments"]["compass"]
    for t, t0 in zip(e["theory"], e["theory_no_dstate"]):
        assert t == pytest.approx(0.9325 * t0, rel=1e-14)


def test_g1d_theory_mstw_is_the_kernels(g1_mod, g1):
    m = _mstw()
    c = _config(g1, "mstw")
    assert c["kernel_max_rel_diff"] == 0.0
    d, pol = lg.deuteron(), lg.ToyG1(m)
    for x, q2 in _PTS[1:]:
        assert g1_mod.g1d_tree(pol, d, x, q2) == pytest.approx(
            0.9325 * 0.5 * (pol.g1p(x, q2) + pol.g1n(x, q2)), rel=1e-14)


# ------------------------------------------------- the numbers, independently

def _raw_points(fname, variant=None):
    body = json.load(open(os.path.join(_DATA, fname)))["body"]
    variant = variant or body["primary_uncertainty_variant"]
    return [(p["x"], p["Q2"], p["value"], p["unc"][variant]) for p in body["points"]]


def test_f2d_toy_chi2_recomputed_from_the_file(f2d):
    """The headline chi2 rebuilt here from the raw JSON and ToyF2, with the
    declared DIS cut, without any harness code."""
    toy = lg.ToyF2()
    c, n = 0.0, 0
    for x, q2, v, u in _raw_points("nnpdf_NMC_NC_NOTFIXED_D.json"):
        if q2 >= 1.0 and 0.938272 ** 2 + q2 * (1 - x) / x >= 4.0:
            t = 0.5 * (toy.f2p(x, q2) + toy.f2n(x, q2))
            c += (v - t) ** 2 / (u["stat"] ** 2 + u["sys"] ** 2)
            n += 1
    d = f2d["configs"]["toy"]["dis"]
    assert n == d["n"] == 155
    assert c == pytest.approx(d["chi2"], rel=1e-12)


def _profiled_chi2(points, theory):
    """min over nuisance parameters lambda_k of
    sum_i (r_i - sum_k lambda_k s_ki)^2 / u_i + sum_k lambda_k^2
    -- the Woodbury form, a different computation from the harness's
    Cholesky solve of (diag + S^T S)."""
    import numpy as np
    names = sorted({k for p in points for k in p["corr"]})
    r = np.array([p["value"] - t for p, t in zip(points, theory)])
    dinv = 1.0 / np.array([p["var_uncorr"] for p in points])
    s = np.array([[p["corr"].get(k, 0.0) for p in points] for k in names])
    a = np.eye(len(names)) + (s * dinv) @ s.T
    b = (s * dinv) @ r
    return float(r @ (dinv * r) - b @ np.linalg.solve(a, b))


@pytest.mark.parametrize("exp", ["hermes", "e143"])
def test_g1d_full_covariance_is_the_profiled_nuisance_chi2(g1_mod, g1, exp):
    pts = [p for p in g1_mod.load_table(exp) if g1_mod.dis_cut(p)]
    th = g1["configs"]["toy"]["experiments"][exp]["theory"]
    th = [th[p["i"]] for p in pts]
    assert _profiled_chi2(pts, th) == pytest.approx(
        g1["configs"]["toy"]["experiments"][exp]["dis"]["chi2_cov"], rel=1e-9)


def test_nmc_full_covariances_are_the_profiled_nuisance_chi2(f2d_mod, f2d, dp_mod, dp):
    pts = [p for p in f2d_mod.load_table() if f2d_mod.dis_cut(p)]
    th = [f2d["configs"]["toy"]["theory"][p["i"]] for p in pts]
    assert _profiled_chi2(pts, th) == pytest.approx(
        f2d["configs"]["toy"]["dis"]["chi2_cov"], rel=1e-9)
    for variant, key in (("hepdata", "chi2_cov"), ("legacy", "chi2_legacy_cov")):
        pts = [p for p in dp_mod.load_table(variant=variant) if dp_mod.dis_cut(p)]
        th = [dp["configs"]["toy"]["theory"][p["i"]] for p in pts]
        assert _profiled_chi2(pts, th) == pytest.approx(
            dp["configs"]["toy"]["dis"][key], rel=1e-9)


def test_hermes_statistics_live_in_the_artificial_vectors(g1_mod):
    """HERMES's stat column is 0; its statistical error is carried by the
    15 CORR vectors -- dropped into the diagonal only via quadrature."""
    pts = g1_mod.load_table("hermes")
    assert all(p["unc"]["stat"] == 0.0 for p in pts)
    assert all(len([k for k in p["corr"] if k.startswith("sys_")]) == 15
               for p in pts)
    assert math.sqrt(sum(v * v for k, v in pts[0]["corr"].items()
                         if k.startswith("sys_"))) == pytest.approx(0.0689, abs=1e-6)


# ------------------------------------------------- verdicts follow from numbers

def test_f2d_verdict_follows_from_the_numbers(f2d):
    assert f2d["headline"] == "toy"
    for name, c in f2d["configs"].items():
        if c["status"] == "blocked":
            continue
        d = c["dis"]
        assert d["p"] == pytest.approx(_sf(d["chi2"], d["n"]), rel=1e-12, abs=1e-300)
        assert c["status"] == ("pass" if _sf(d["chi2"], d["n"]) >= 0.01 else "fail")
    assert f2d["status"] == f2d["configs"]["toy"]["status"] == "fail"


def test_dp_verdict_follows_from_the_numbers(dp):
    assert dp["headline"] == "toy"
    for name, c in dp["configs"].items():
        if c["status"] == "blocked":
            continue
        d = c["dis"]
        assert c["status"] == ("pass" if _sf(d["chi2"], d["n"]) >= 0.01 else "fail")
    assert dp["status"] == dp["configs"]["toy"]["status"] == "fail"


def test_g1d_verdict_follows_from_the_numbers(g1):
    assert g1["headline"] == "toy"
    for name, c in g1["configs"].items():
        if c["status"] == "blocked":
            continue
        for e, s in c["experiments"].items():
            d = s["dis"]
            assert s["status"] == ("pass" if _sf(d["chi2"], d["n"]) >= 0.01
                                   else "fail")
            p_cov = d["p"] if d["p_cov"] is None else d["p_cov"]
            assert s["verdict_differs_under_full_cov"] == (
                (d["p"] >= 0.01) != (p_cov >= 0.01))
        # the verdict is COMPASS, HERMES and E143; SMC is recorded only
        assert c["status"] == ("pass" if all(
            c["experiments"][e]["status"] == "pass"
            for e in ("compass", "hermes", "e143")) else "fail")
        assert c["experiments"]["smc"]["in_verdict"] is False
    assert g1["status"] == g1["configs"]["toy"]["status"] == "fail"


def test_chi2_sf_fallback_is_scipy_where_the_verdicts_and_headlines_live(
        f2d_mod, dp_mod, g1_mod):
    """Without scipy the harnesses use `chi2_sf_fallback`: the continued
    fraction for Q itself above chi2/2 = ndf/2 + 1, so a p of 1e-16 or 1e-26
    is the number scipy prints, not 1 - P noise, and chi2 = 0 is p = 1
    (until 2026-09-27: 'compass 107.8/15 (p = 0)' without scipy, and a math
    domain error at chi2 = 0).  One function, three harnesses."""
    pytest.importorskip("scipy", reason="scipy (the reference for the "
                                        "no-scipy chi2_sf fallback)")
    f = f2d_mod.chi2_sf_fallback
    for m in (dp_mod, g1_mod):
        assert m.chi2_sf_fallback.__code__.co_code == f.__code__.co_code
    assert f(0.0, 5) == 1.0 and f(float("inf"), 5) == 0.0
    with pytest.raises(ValueError):
        f(-1.0, 5)
    cases = [(1e-8, 1), (0.5, 1), (1.0, 1), (3.2, 4), (10.0, 10),
             (30.578, 15), (247.18, 211), (4.6269, 4), (26.941, 15),
             (201.0, 200), (107.84918016285502, 15), (416.46, 155),
             (69.4, 15), (54.9, 28), (53.2, 15), (500.0, 3), (5000.0, 4000)]
    for chi2, n in cases:
        assert f(chi2, n) == pytest.approx(_sf(chi2, n), rel=1e-10), (chi2, n)
    # the threshold side of the verdict is the same either way
    assert (f(30.578, 15) >= 0.01) == (_sf(30.578, 15) >= 0.01)
    assert f(8307.9, 155) == 0.0 == _sf(8307.9, 155)


# ------------------------------------------------- pinned measurements

def test_f2d_pinned_toy(f2d):
    c = f2d["configs"]["toy"]
    assert (c["dis"]["n"], c["all"]["n"]) == (155, 158)
    assert c["dis"]["chi2"] == pytest.approx(8307.8965, rel=1e-6)
    assert c["dis"]["chi2_cov"] == pytest.approx(21100.704, rel=1e-6)
    assert c["all"]["chi2"] == pytest.approx(9582.5732, rel=1e-6)
    assert c["dis"]["mean_ratio"] == pytest.approx(1.0601823, rel=1e-6)
    assert c["dis"]["worst_pull"] == pytest.approx(21.461662, rel=1e-6)
    assert (c["dis"]["worst_x"], c["dis"]["worst_q2"]) == (0.0125, 1.25)


def test_f2d_pinned_mstw(f2d):
    c = _config(f2d, "mstw")
    assert c["dis"]["n"] == 155
    assert c["dis"]["chi2"] == pytest.approx(416.45940, rel=1e-6)
    assert c["dis"]["chi2_cov"] == pytest.approx(858.53643, rel=1e-6)
    assert c["all"]["chi2"] == pytest.approx(485.45076, rel=1e-6)
    assert c["status"] == "fail"


def test_dp_pinned_toy(dp):
    c = dp["configs"]["toy"]
    assert (c["dis"]["n"], c["all"]["n"]) == (211, 260)
    assert c["dis"]["chi2"] == pytest.approx(2161.4387, rel=1e-6)
    assert c["dis"]["chi2_cov"] == pytest.approx(1154.2873, rel=1e-6)
    assert c["dis"]["chi2_legacy_cov"] == pytest.approx(950.70699, rel=1e-6)
    assert c["all"]["chi2"] == pytest.approx(2251.9861, rel=1e-6)


def test_dp_pinned_mstw(dp):
    c = _config(dp, "mstw")
    assert c["dis"]["chi2"] == pytest.approx(247.17894, rel=1e-6)
    assert c["dis"]["p"] == pytest.approx(0.0444218, rel=1e-5)
    assert c["status"] == "pass"             # a sub-row; the row is toy's
    assert c["dis"]["chi2_cov"] == pytest.approx(251.63526, rel=1e-6)
    assert c["dis"]["chi2_legacy_cov"] == pytest.approx(224.67659, rel=1e-6)
    assert c["dis_frozen"]["chi2"] == pytest.approx(260.70160, rel=1e-6)
    assert c["all"]["chi2"] == pytest.approx(297.89628, rel=1e-6)
    assert c["all_frozen"]["chi2"] == pytest.approx(310.17827, rel=1e-6)


_G1_TOY = dict(compass=(15, 107.84918, None, 126.36996),
               hermes=(15, 69.383903, 104.45971, 71.167262),
               e143=(28, 54.927009, 51.526179, 55.509263),
               smc=(12, 16.652517, None, 16.503757))
_G1_MSTW = dict(compass=(15, 53.237005, None, 66.298574),
                hermes=(15, 47.134153, 68.827103, 64.040126),
                e143=(28, 50.243680, 51.339501, 61.496627),
                smc=(12, 14.401418, None, 14.756620))


def _check_g1(c, table):
    for e, (n, chi2, cov, nod) in table.items():
        s = c["experiments"][e]
        assert s["dis"]["n"] == n
        assert s["dis"]["chi2"] == pytest.approx(chi2, rel=1e-6)
        if cov is None:
            assert s["dis"]["chi2_cov"] is None
        else:
            assert s["dis"]["chi2_cov"] == pytest.approx(cov, rel=1e-6)
        assert s["dis_no_dstate"]["chi2"] == pytest.approx(nod, rel=1e-6)


def test_g1d_pinned_toy(g1):
    c = g1["configs"]["toy"]
    _check_g1(c, _G1_TOY)
    assert (c["world_chi2"], c["world_ndf"]) == (pytest.approx(232.16009, rel=1e-6), 58)
    # the scout's diagnosis: COMPASS's x < 0.01 points drive the toy's chi2
    e = c["experiments"]["compass"]
    assert e["theory"][0] == pytest.approx(-0.96277, abs=5e-5)
    assert e["dis"]["worst_x"] < 0.01


def test_g1d_pinned_mstw(g1):
    c = _config(g1, "mstw")
    _check_g1(c, _G1_MSTW)
    assert c["world_chi2"] == pytest.approx(150.61484, rel=1e-6)
    assert c["status"] == "fail"


# ------------------------------------------------- blocked, never dropped

def _lhapdf_missing(setname):
    """True iff LHAPDF reports `setname` not installed; any other error
    propagates."""
    if not getattr(lg, "HAVE_LHAPDF", False):
        return True
    try:
        lg.LhapdfSF(setname)
    except RuntimeError as exc:
        if "Info file not found for PDF set" not in str(exc):
            raise
        return True
    return False


@pytest.mark.parametrize("rep_name,config,setname", [
    ("f2d", "ct18nlo", "CT18NLO"), ("dp", "ct18nlo", "CT18NLO"),
    ("g1", "ct18nlo", "CT18NLO"), ("g1", "nnpdfpol", "NNPDFpol11_100")])
def test_lhapdf_configs_are_reported_not_dropped(request, rep_name, config, setname):
    rep = request.getfixturevalue(rep_name)
    c = rep["configs"][config]
    if not getattr(lg, "HAVE_LHAPDF", False):
        assert c["status"] == "blocked"
        assert c["reason"] == ("environment: lipolgen built without the LHAPDF "
                               "tier (HAVE_LHAPDF false)")
    elif _lhapdf_missing(setname):
        assert c["status"] == "blocked"
        assert c["reason"] == "environment: LHAPDF set %s not installed" % setname
    else:                                   # the set is there: a real sub-row
        assert c["status"] in ("pass", "fail")
        assert c["q2_floor"] >= 1.0


@pytest.mark.parametrize("harness", ["t1_nmc_f2d", "t1_nmc_f2d_over_f2p",
                                     "t1_g1d_world"])
def test_harness_prints_one_report_row_and_every_subrow(harness, capsys):
    mod = _load(harness)
    rep = mod.run(verbose=True)
    out = capsys.readouterr().out
    rows = [l for l in out.splitlines() if l.startswith("REPORT |")]
    assert len(rows) == 1
    fields = [f.strip() for f in rows[0].split("|")]
    assert fields[1] == rep["name"] == harness and len(fields) == 6
    assert fields[5] == rep["status"].upper() == "FAIL"
    assert "np." not in out                        # numpy scalars never leak
    sub = [l for l in out.splitlines() if l.startswith("SUBROW |")]
    per = len(mod.EXPERIMENTS) if hasattr(mod, "EXPERIMENTS") else 1
    assert len(sub) == len(mod.CONFIGS) * per      # blocked ones included
    for name, c in rep["configs"].items():
        if c["status"] == "blocked":
            assert any(("[%s" % name) in l and l.endswith("BLOCKED (%s)" % c["reason"])
                       for l in sub)
    assert rep["runtime_s"] < 30.0                 # no LIPOLGEN_BENCH needed


# ------------------------------------------------- the exit convention

@pytest.mark.parametrize("harness", ["t1_nmc_f2d", "t1_nmc_f2d_over_f2p",
                                     "t1_g1d_world"])
def test_harness_exits_0_on_a_recorded_fail(harness):
    env = dict(os.environ)
    here = os.path.dirname(os.path.dirname(os.path.abspath(lg.__file__)))
    env["PYTHONPATH"] = os.pathsep.join(
        [here] + [p for p in env.get("PYTHONPATH", "").split(os.pathsep) if p])
    res = subprocess.run([sys.executable, os.path.join(_BENCH, harness + ".py")],
                         capture_output=True, text=True, env=env, timeout=120)
    assert res.returncode == 0, res.stderr
    assert [l for l in res.stdout.splitlines()
            if l.startswith("REPORT | %s |" % harness)][0].endswith("| FAIL")


def test_harness_exits_2_on_an_altered_file(f2d_mod, dp_mod, g1_mod, tmp_path,
                                            monkeypatch, capsys):
    for mod, fname in ((f2d_mod, "nnpdf_NMC_NC_NOTFIXED_D.json"),
                       (dp_mod, "nnpdf_NMC_NC_NOTFIXED.json")):
        bad = tmp_path / fname
        bad.write_bytes(_flip_first_value(open(mod.DATA, "rb").read()))
        monkeypatch.setattr(mod, "DATA", str(bad))
        assert mod.main() == 2
    for fname, *_ in _FILES[2:]:
        (tmp_path / fname).write_bytes(open(os.path.join(_DATA, fname), "rb").read())
    hermes = tmp_path / "nnpdf_HERMES_NC_7GEV_ED.json"
    hermes.write_bytes(_flip_first_value(hermes.read_bytes()))
    monkeypatch.setattr(g1_mod, "DATA_DIR", str(tmp_path))
    assert g1_mod.main() == 2
    err = capsys.readouterr().err.splitlines()
    assert len([l for l in err if l.startswith("RuntimeError: ")
                and "vendored file sha256" in l]) == 3
