# SPDX-License-Identifier: GPL-3.0-or-later
"""The proposed new T2 row t2_deeps_spectator_tail.py -- CLAS Deeps
F2N x P(p_s) spectator tail against the tree's deuteron n(k) -- exercised as
tests.  Cheap (~2 s), so not behind LIPOLGEN_BENCH.

A harness that DISAGREES with the tree is not a failing test: these pin what
the harness MEASURES, that the verdict follows from the numbers, that the
vendored file is the one recorded, and that the shape normalisation really
cancels a common factor (the F2N argument of sec. 8 rule 3).
"""
import importlib.util
import os

import pytest

lg = pytest.importorskip("lipolgen")

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_BENCH = os.path.join(_ROOT, "validation", "benchmarks")


def _load(name):
    spec = importlib.util.spec_from_file_location(
        "bench_" + name, os.path.join(_BENCH, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


deeps = _load("t2_deeps_spectator_tail")


@pytest.fixture(scope="module")
def rep():
    return deeps.run(verbose=False)


# ------------------------------------------------------------ the data file

def test_digest_is_the_one_the_bonus_row_records():
    bonus = _load("t2_bonus_spectator_shape")
    assert deeps.DEEPS_SHA256 == bonus.DEEPS_SHA256
    assert os.path.samefile(deeps.DEEPS, bonus.DEEPS)


def test_vendored_file_guard_refuses_one_flipped_byte(tmp_path):
    raw = open(deeps.DEEPS, "rb").read()
    old, new = b"0.00156769", b"0.00156768"
    assert raw.count(old) == 1
    good = tmp_path / "good.json"
    good.write_bytes(raw)
    deeps.load_deeps(str(good))           # content, not path, is checked
    bad = tmp_path / "bad.json"
    bad.write_bytes(raw.replace(old, new))
    with pytest.raises(RuntimeError, match="sha256"):
        deeps.load_deeps(str(bad))


def test_block_census_pinned():
    c = deeps.census(deeps.load_deeps())
    assert c["blocks_total"] == 115
    assert c["blocks_usable"] == 60
    assert c["blocks_not_read"] == 55
    assert c["q2"] == [1.8, 2.8]
    assert c["w_star"] == [0.94, 1.25, 1.5, 1.73, 2.02, 2.4]
    assert c["ps_labels_MeV"] == [300.0, 340.0, 390.0, 460.0, 560.0]
    assert c["rows_usable_blocks"] == 643
    assert c["rows_missing_dropped"] == 37
    assert c["rows_negative_stat0_dropped"] == 2
    assert c["rows_negative_kept"] == 0
    assert c["rows_kept"] == 604


def test_points_use_paper_ps_and_drop_every_stat_zero_row():
    pts = deeps.points(deeps.load_deeps())
    assert len(pts) == 604
    assert all(p["stat"] > 0 for p in pts)
    assert sorted({p["ps"] for p in pts}) == [0.30, 0.34, 0.39, 0.46, 0.53]
    assert all(p["q2"] in (1.8, 2.8) for p in pts)


def test_duplicate_paper_block_is_refused():
    d = deeps.load_deeps()
    d["blocks"].append(dict(d["blocks"][0]))
    with pytest.raises(RuntimeError, match="duplicate"):
        deeps.usable_blocks(d)


# ------------------------------------------------------------ the method

def test_scaled_chi2_cancels_a_common_factor():
    """A common factor on the data (F2N at fixed x*) or on the model (the
    n(k) normalisation) leaves chi2 unchanged -- the shape argument."""
    d = [3.0, 2.1, 1.2, 0.55, 0.2]
    m = [2.9, 2.2, 1.1, 0.6, 0.18]
    e = [0.2, 0.15, 0.1, 0.06, 0.03]
    c0, ndf, s0 = deeps.scaled_chi2(d, m, e)
    assert ndf == 4
    c1, _, s1 = deeps.scaled_chi2([7.3 * x for x in d], m,
                                  [7.3 * x for x in e])
    c2, _, s2 = deeps.scaled_chi2(d, [0.01 * x for x in m], e)
    assert c1 == pytest.approx(c0, rel=1e-12)
    assert c2 == pytest.approx(c0, rel=1e-12)
    assert s1 == pytest.approx(7.3 * s0, rel=1e-12)
    assert s2 == pytest.approx(100 * s0, rel=1e-12)
    # a perfect shape at any scale is chi2 = 0
    assert deeps.scaled_chi2([4 * x for x in m], m, e)[0] == pytest.approx(0, abs=1e-20)


def test_shape_chi2_is_blind_to_a_per_q2_w_factor():
    """Multiplying the data of every (Q2, W*) by its own F2N-like factor (and
    its errors with it) does not move the combined shape chi2."""
    pts = deeps.points(deeps.load_deeps())
    tree = deeps.Tree()
    model = lambda p: tree.p_model("hulthen", p["ps"], p["c"])   # noqa: E731
    base = deeps.shape_chi2(pts, model)
    scaled = []
    for p in pts:
        f = 0.3 + p["q2"] * p["w"]          # any function of (Q2, W*) alone
        scaled.append(dict(p, val=f * p["val"], err=f * p["err"]))
    assert deeps.shape_chi2(scaled, model)["chi2"] == pytest.approx(
        base["chi2"], rel=1e-10)


def test_light_cone_map():
    t = deeps.Tree()
    # backward spectator: alpha > 1 and k_LC > p_s
    a, pt, k, flux = t.lc(0.30, -1.0)
    assert a > 1 and pt == 0.0 and k == pytest.approx(0.373, abs=2e-3)
    assert t.lc(0.53, -1.0)[2] == pytest.approx(0.95, abs=0.01)
    # alpha_s = 1, p_T = 0 is k = 0
    assert t.lc(1e-9, 0.0)[2] == pytest.approx(0.0, abs=1e-3)


def test_tree_density_isotropic_and_waves_differ():
    t = deeps.Tree()
    assert t.isotropy("hulthen", 0.4) == pytest.approx(1.0, abs=1e-12)
    assert t.isotropy("av18", 0.4) == pytest.approx(1.0, abs=1e-12)
    # CD-Bonn read directly: psi_s^2 + psi_d^2 falls through the window
    r = [t.n_of_k("cdbonn", k) for k in (0.3, 0.4, 0.5)]
    assert r[0] > r[1] > r[2] > 0


def test_fsi_weight_is_live_for_the_deuteron_channel():
    """The Glauber weight must not be the channel-guard's silent 1."""
    t = deeps.Tree()
    assert t.fsi_weight(0.46, 0.0, 2.0, 2.8) > 5.0
    assert t.fsi_weight(0.46, -0.8, 2.0, 2.8) < 1.0


# ------------------------------------------------------------ measured / verdict

def test_report_row_shape(rep):
    for k in ("name", "generator", "reference", "tolerance", "status",
              "subrows"):
        assert k in rep, k
    assert rep["name"] == "t2_deeps_spectator_tail"
    assert rep["status"] in ("pass", "fail", "blocked")
    assert [s for s, v in rep["subrows"].items() if v["in_verdict"]] == [
        "hulthen/lc"]


def test_verdict_follows_from_the_numbers(rep):
    from scipy.stats import chi2 as chi2_dist
    head = rep["measured"]["waves"]["hulthen"]["lc"]
    assert head["p"] == pytest.approx(chi2_dist.sf(head["chi2"], head["ndf"]))
    assert rep["status"] == ("pass" if head["p"] >= deeps.P_MIN else "fail")
    for name, s in rep["subrows"].items():
        if "p" in s:
            assert s["status"] == ("pass" if s["p"] >= deeps.P_MIN else "fail"), name


MEASURED = {  # (wave, prescription): chi2 over 187 dof, 2026-09-27
    ("hulthen", "lc"): 627.848820063069,
    ("hulthen", "instant"): 322.5414096505045,
    ("av18", "lc"): 579.1162236521793,
    ("av18", "instant"): 805.15668556147,
    ("cdbonn", "lc"): 278.8119735316571,
    ("cdbonn", "instant"): 359.8814090175042,
}


def test_measured_values_pinned(rep):
    m = rep["measured"]
    assert m["n_points"] == 604 and m["n_points_pwia"] == 245
    for (wv, form), chi2 in MEASURED.items():
        r = m["waves"][wv][form]
        assert r["ndf"] == 187 and r["groups_used"] == 57
        assert r["chi2"] == pytest.approx(chi2, rel=1e-6), (wv, form)
    assert m["lc_times_glauber"]["chi2"] == pytest.approx(422.92334145, rel=1e-6)
    assert rep["status"] == "fail"
    by_ps = m["waves"]["hulthen"]["lc"]["data_over_model_by_ps"]
    assert by_ps == pytest.approx([0.790, 0.954, 1.271, 1.715, 1.409], abs=2e-3)


def test_transverse_window_recorded_not_asserted(rep):
    tr = rep["measured"]["transverse"]
    assert rep["subrows"]["transverse window (recorded)"]["in_verdict"] is False
    assert tr["n_points"] == 359
    d = [v["data_over_pwia"] for v in tr["per_ps"].values()]
    g = [v["tree_glauber"] for v in tr["per_ps"].values()]
    assert d == pytest.approx([0.760, 1.198, 1.760, 2.642, 2.957], abs=2e-3)
    assert g == pytest.approx([1.681, 3.115, 6.048, 12.404, 18.936], abs=2e-3)
    assert [v["same_side_of_1"] for v in tr["per_ps"].values()] == [
        False, True, True, True, True]
    assert sum(v["n"] for v in tr["per_ps"].values()) == 359


def test_harness_exit_zero_and_report_line():
    import subprocess
    import sys
    p = subprocess.run([sys.executable, os.path.join(
        _BENCH, "t2_deeps_spectator_tail.py")], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    lines = [l for l in p.stdout.splitlines() if l.startswith("REPORT | ")]
    assert len(lines) == 1
    assert lines[0].startswith("REPORT | t2_deeps_spectator_tail | ")
    assert lines[0].endswith("| FAIL")
