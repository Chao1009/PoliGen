# SPDX-License-Identifier: GPL-3.0-or-later
"""BENCHMARK_PLAN.md sec. 4 rows 4, 5, 6 -- the NUCLEAR-INPUTS harnesses under
validation/benchmarks/ (tier T3), exercised as tests.

All three are cheap (< 2 s together), so none sits behind LIPOLGEN_BENCH.
A harness that DISAGREES with the tree is not a failing test: the test pins
what the harness MEASURES and that its verdict follows from the numbers.  No
test here moves, or asserts a new value for, any tree default.

  row 4  t3_nmc_li6_over_d.py      NMC F2(6Li)/F2(D) vs EPPS21 (isoscalar)
  row 5  t3_li6_charge_ff_fb.py    UVa FB charge density vs HoSpin1FF C0
  row 6  t3_li_magnetization_rfy.py  BLOCKED until RFY66 / ADNDT 14 lands
"""
import hashlib
import importlib.util
import math
import os

import pytest

lg = pytest.importorskip("lipolgen")

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_BENCH = os.path.join(_ROOT, "validation", "benchmarks")


def _load(name):
    path = os.path.join(_BENCH, name + ".py")
    spec = importlib.util.spec_from_file_location("bench_" + name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ------------------------------------------------------------------ row 4

@pytest.fixture(scope="module")
def nmc_mod():
    """The row-4 harness module alone: its vendored-table check needs neither
    LHAPDF nor the sets, so it runs everywhere."""
    return _load("t3_nmc_li6_over_d")


def _missing_lhapdf_sets(mod):
    """The row's sets LHAPDF reports as not installed ("Info file not found
    for PDF set ..."), probed one by one.  Any OTHER error propagates."""
    missing = []
    for name in (mod.NUCLEAR_SET, mod.PROTON_SET):
        try:
            lg.LhapdfSF(name)
        except RuntimeError as exc:
            if "Info file not found for PDF set" not in str(exc):
                raise
            missing.append(name)
    return missing


@pytest.fixture(scope="module")
def nmc(nmc_mod):
    """Row 4 needs LHAPDF with EPPS21nlo_CT18Anlo_Li6 and CT18NLO.  It skips,
    naming them, ONLY when the LHAPDF tier is not built or LHAPDF reports a
    set not installed; any other error in `Tree()` or `run()` (a renamed or
    re-signatured binding, a harness bug) FAILS instead of hiding as a skip."""
    mod = nmc_mod
    if not getattr(lg, "HAVE_LHAPDF", False):
        pytest.skip("environment: lipolgen built without the LHAPDF tier "
                    "(HAVE_LHAPDF false); row 4 needs LHAPDF sets %s and %s"
                    % (mod.NUCLEAR_SET, mod.PROTON_SET))
    missing = _missing_lhapdf_sets(mod)
    if missing:
        pytest.skip("; ".join("environment: LHAPDF set %s not installed" % n
                              for n in missing)
                    + " (row 4 needs %s and %s)" % (mod.NUCLEAR_SET,
                                                    mod.PROTON_SET))
    tree = mod.Tree()
    return mod, tree, mod.run(verbose=False)


def test_nmc_vendored_table_is_the_download(nmc_mod, tmp_path):
    mod = nmc_mod
    pts = mod.load_table()
    assert len(pts) == 24
    assert pts[0]["x"] == 0.00014 and pts[-1]["x"] == 0.65
    # one flipped byte in the body is refused, not silently read
    raw = open(mod.DATA, "rb").read()
    assert raw.count(b"0.801,") == 1
    raw = raw.replace(b"0.801,", b"0.802,")
    bad = tmp_path / "t1.csv"
    bad.write_bytes(raw)
    with pytest.raises(RuntimeError, match="sha256"):
        mod.load_table(str(bad))


def test_nmc_rule2_denominator_is_isoscalar(nmc):
    """BENCHMARK_PLAN.md sec. 8 rule 2, ASSERTED: the denominator is
    (F2p + F2n)/2 of the free nucleon -- and it is NOT the free proton."""
    mod, tree, _ = nmc
    free = lg.LhapdfSF(mod.PROTON_SET)
    for x, q2 in ((0.0125, 1.8), (0.1, 5.0), (0.35, 23.0), (0.65, 39.0)):
        iso = 0.5 * (free.f2p(x, q2) + free.f2n(x, q2))
        assert tree.isoscalar_denominator(x, q2) == pytest.approx(iso, rel=1e-15)
        assert tree.ratio(x, q2) == pytest.approx(
            tree.epps.f2_per_nucleon(x, q2) / iso, rel=1e-15)
    # at large x the two denominators differ by far more than the EMC effect
    assert tree.ratio(0.65, 5.0) - tree.ratio_WRONG_free_proton(0.65, 5.0) > 0.2


def test_nmc_chi2_reproduces_the_survey_and_passes(nmc):
    _, _, rep = nmc
    # 03_data_nuclear.md sec. 1.2: chi2 = 4.63 for the 4 points x >= 0.30
    assert rep["ndf_valence"] == 4
    assert rep["chi2_valence"] == pytest.approx(4.626864, abs=5e-6)
    assert rep["ndf_ongrid"] == 15 and rep["n_below_grid"] == 9
    assert rep["chi2_ongrid"] == pytest.approx(26.94125, abs=5e-5)
    assert rep["status"] == ("pass" if min(rep["p_valence"], rep["p_ongrid"]) >= 0.01
                             else "fail")
    assert rep["status"] == "pass"
    assert rep["nmc_valence_mean"] == pytest.approx(0.961837, abs=5e-6)
    assert rep["nmc_valence_mean_err"] == pytest.approx(0.021863, abs=5e-6)


def test_nmc_wrong_denominator_trap_is_nine_times(nmc):
    """03 sec. 1.3: the free-proton denominator gives <1 - R> ~ 0.268, ~9x the
    isoscalar value -- recorded, and never in the verdict."""
    _, _, rep = nmc
    assert rep["depletion_WRONG"] == pytest.approx(0.26763, abs=5e-5)
    assert rep["depletion_iso_ct18nlo"] == pytest.approx(0.029803, abs=5e-6)
    assert 8.0 < rep["depletion_WRONG"] / rep["depletion_iso_ct18nlo"] < 10.0
    assert rep["chi2_valence_WRONG"] > 50.0
    assert rep["r_u_065"] == pytest.approx(0.578, abs=5e-4)
    assert rep["r_d_065"] == pytest.approx(2.829, abs=5e-4)
    # the library constant is untouched by this harness
    assert lg.EMC_VALENCE_DEPLETION_EPPS21 == 0.031052077003862335


# ------------------------------------------------------------------ row 5

@pytest.fixture(scope="module")
def fb():
    mod = _load("t3_li6_charge_ff_fb")
    return mod, mod.run(verbose=False)


def test_fb_row_is_the_vendored_row(fb, tmp_path):
    mod, rep = fb
    row = rep["row"]
    assert (row["A"], row["Z"], row["R"], len(row["a"])) == (6, 3, 6.0, 7)
    raw = open(mod.DATA, "rb").read().replace(b"1.6353e-02", b"1.6354e-02")
    bad = tmp_path / "fb.dat"
    bad.write_bytes(raw)
    with pytest.raises(RuntimeError, match="sha256"):
        mod.load_row(str(bad))
    # since 2026-09-27 every byte is pinned, the provenance header included:
    # one digit of the A/Z note, or the licence tag, is refused too
    orig = open(mod.DATA, "rb").read()
    assert hashlib.sha256(orig).hexdigest() == mod.FILE_SHA256
    for old, new in ((b"The harness reads A = 6, Z = 3.",
                      b"The harness reads A = 7, Z = 3."),
                     (b"SPDX-License-Identifier: GPL-3.0-or-later",
                      b"SPDX-License-Identifier: LicenseRef-NoneStated")):
        assert orig.count(old) == 1, old
        bad.write_bytes(orig.replace(old, new))
        with pytest.raises(RuntimeError, match="vendored file sha256"):
            mod.load_row(str(bad))


def test_fb_density_numbers(fb):
    mod, rep = fb
    assert rep["fb_charge"] == pytest.approx(2.991170, abs=5e-7)
    assert rep["fb_rms"] == pytest.approx(2.52061, abs=5e-5)
    assert rep["q_fb"] == pytest.approx(2.694413, abs=2e-6)
    # closed form vs direct quadrature of rho(r) at one q
    q = 1.3
    n = 20000
    h = rep["row"]["R"] / n
    f = mod.FourierBessel(rep["row"]["a"], rep["row"]["R"])
    quad = 4 * math.pi * h * sum(
        f.rho((i + 0.5) * h) * math.sin(q * (i + 0.5) * h) / q * (i + 0.5) * h
        for i in range(n))
    assert f.ff(q) == pytest.approx(quad, rel=1e-6)


def test_fb_tree_side_and_verdict(fb):
    mod, rep = fb
    t = rep["tree"]
    # the shipped zero, inside the (unmoved) T11 window
    assert mod.T11_WINDOW[0] <= t["q0"] <= mod.T11_WINDOW[1]
    assert t["q0"] == pytest.approx(3.0998, abs=1e-4)
    # the tree's own |F_L|^2 minimum IS its C0 zero: C2 shares F_point
    assert abs(t["q_min_fl2"] - t["q0"]) < 1e-6
    assert t["rms_folded"] == pytest.approx(2.5710, abs=1e-4)
    # the verdict follows from the numbers, and carries its distance
    lo, hi = rep["band"]
    assert hi == pytest.approx(math.sqrt(8.0))
    assert rep["inside"] == (lo <= t["q0"] <= hi)
    assert rep["status"] == ("pass" if rep["inside"] else "fail")
    assert rep["distance_above_band"] == pytest.approx(t["q0"] - hi, abs=1e-15)
    assert rep["reason"] is None and t["no_zero"] is None


def test_fb_tree_c0_without_a_zero_is_a_recorded_fail(monkeypatch, capsys):
    """A tree C0 with no zero below the scan limit (rc.hpp's VmcFt shape has
    none) contradicts the band: a recorded FAIL carrying that reason, with
    its REPORT row printed -- not a crash, which would read as a broken
    harness (exit 2).  The reference side does not move."""
    mod = _load("t3_li6_charge_ff_fb")
    orig, opt = lg.HoSpin1FF, lg.HoSpin1FFOptions()
    opt.c0_shape = lg.C0Shape.VmcFt

    class _VmcFtEdge:
        @staticmethod
        def for_ion(ion):
            return orig.for_ion(ion, opt)

    monkeypatch.setattr(lg, "HoSpin1FF", _VmcFtEdge)
    rep = mod.run(verbose=True)
    out = capsys.readouterr().out
    assert rep["status"] == "fail" and rep["inside"] is False
    assert rep["tree"]["q0"] is None and rep["distance_above_band"] is None
    assert rep["reason"] == "the tree's C0 has no zero below 5 fm^-1"
    assert rep["q_fb"] == pytest.approx(2.694413, abs=2e-6)
    assert ("REPORT | t3_li6_charge_ff_fb | HoSpin1FF C0: NO zero below 5 "
            "fm^-1" in out)
    assert out.rstrip().endswith("| FAIL (recorded, the tree's C0 has no zero "
                                 "below 5 fm^-1; nothing moved)")


# ------------------------------------------------------------------ row 6

def test_rfy_tree_side_closed_form():
    """Measured even while the row is blocked: the tree's F_m slope radius is
    the closed form 6/q_z^2 + (3/2) b^2 of rc.hpp's F_mag."""
    mod = _load("t3_li_magnetization_rfy")
    t = mod.tree_side()
    assert t["r_mag_li6"] == pytest.approx(t["r_mag_li6_closed"], rel=1e-5)
    assert t["r_mag_li6_closed"] == pytest.approx(2.94688, abs=5e-5)
    assert t["li7"].startswith("not implemented in the tree")


def test_rfy_magnetization_benchmark():
    mod = _load("t3_li_magnetization_rfy")
    rep = mod.run(verbose=False)
    if rep["status"] == "blocked":
        assert "Phys. Rev. 144 (1966) 859" in rep["reason"]
        assert "10.1016/S0092-640X(74)80002-1" in rep["reason"]
        pytest.skip(rep["reason"])
    assert rep["status"] in ("pass", "fail")      # pragma: no cover
