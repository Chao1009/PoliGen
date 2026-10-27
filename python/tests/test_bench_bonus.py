# SPDX-License-Identifier: GPL-3.0-or-later
"""BENCHMARK_PLAN.md sec. 4 row 10 -- t2_bonus_spectator_shape.py, CLAS
BONuS R_D/S against cos theta_pq in the DIS blocks, held against the tree's
GlauberFsiWeight on the deuteron control -- exercised as tests (run
2026-10-26, phase B1; the row was BLOCKED until then).

The harness costs ~2.5 s (a library evaluation of the weight, no
generation), so nothing here sits behind LIPOLGEN_BENCH.

A harness that DISAGREES with the tree is not a failing test: these pin what
the harness MEASURES (the row's recorded FAIL is pinned as a FAIL), check
that every verdict follows from those numbers by the rule the header
declares, guard every vendored file it reads (count, first and last value,
one flipped byte refused with exit 2), and check the evaluation itself (the
weight live and symmetric, converged in its table; the density the tree's
own; the light-cone map the Deeps row's).  Nothing here moves, or asserts a
new value for, any tree default; nothing needs the network.
"""
import hashlib
import importlib.util
import json
import math
import os
import shutil
import subprocess
import sys

import pytest

lg = pytest.importorskip("lipolgen")

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_BENCH = os.path.join(_ROOT, "validation", "benchmarks")
_DATA = os.path.join(_BENCH, "data")
NAME = "t2_bonus_spectator_shape"


def _load(name):
    spec = importlib.util.spec_from_file_location(
        "bench_" + name, os.path.join(_BENCH, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


bonus = _load(NAME)


@pytest.fixture(scope="module")
def rep():
    return bonus.run(verbose=False)


@pytest.fixture(scope="module")
def tree():
    return bonus.Tree()


def _env():
    env = dict(os.environ)
    here = os.path.dirname(os.path.dirname(os.path.abspath(lg.__file__)))
    env["PYTHONPATH"] = os.pathsep.join(
        [here] + [p for p in env.get("PYTHONPATH", "").split(os.pathsep) if p])
    return env


def _run_script(script, *args, **env_extra):
    env = _env()
    env.update(env_extra)
    return subprocess.run([sys.executable, script] + list(args),
                          capture_output=True, env=env, timeout=300)


def _one_flipped_byte_is_refused(tmp_path, loader, src, old, new):
    """A byte-identical copy loads (content, not path, is checked); the same
    copy with ONE byte changed at ONE place raises the sha256 error."""
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


# ------------------------------------------------------------ the files

def test_bonus_file_is_the_record_and_refuses_a_flipped_byte(tmp_path):
    raw = open(bonus.BONUS, "rb").read()
    assert hashlib.sha256(raw).hexdigest() == bonus.BONUS_SHA256
    d = bonus.load_bonus()
    prov = d["_provenance"]
    assert prov["doi"] == "10.1103/PhysRevC.89.045206"
    assert prov["url"].endswith("10.1103/PhysRevC.89.045206/BONuSprcSupMat2.txt")
    assert prov["fetched"] == "2026-10-26"
    assert prov["original_file"] == dict(
        name="BONuSprcSupMat2.txt", bytes=88845, lines=5452,
        format="ASCII, tab-separated, LF line ends, blank line between rows",
        sha256=bonus.ORIGINAL_TXT_SHA256)
    assert "facts quoted with attribution" in prov["licence"]
    secs = d["ratio_blocks"]["sections"]
    assert [s["section"] for s in secs] == [1, 2]
    assert [len(s["blocks"]) for s in secs] == [60, 40]
    assert [s["beam_label"] for s in secs] == ["4 GeV", "5 GeV"]
    assert all(b["beam"] == s["beam_label"] for s in secs for b in s["blocks"])
    assert all(len(b["rows"]) == 10 for s in secs for b in s["blocks"])
    first, last = secs[0]["blocks"][0], secs[1]["blocks"][-1]
    assert (first["q2"], first["w_star"], first["ps"], first["src_line"]) == (
        0.93, 1.17, 0.078, 13)
    assert first["rows"][0] == [-0.9, 1.262, 0.068, 0.177]
    assert (last["q2"], last["w_star"], last["ps"], last["src_line"]) == (
        3.38, 2.44, 0.135, 2793)
    assert last["rows"][-1] == [0.9, 1.259, 0.063, 0.154]
    tabs = d["dis_f2n_ratios"]["tables"]
    assert [t["cut_verbatim"] for t in tabs] == ["W > 1.4 GeV", "W > 1.6 GeV",
                                                 "W > 1.8 GeV"]
    assert [len(t["rows"]) for t in tabs] == [13, 12, 10]
    assert tabs[0]["rows"][0] == [2.4338, 1.0867, 0.177, 0.4255, 0.0035,
                                  0.0214, 0.7314, 0.0062, 0.0383]
    assert tabs[-1]["rows"][-1] == [1.8377, 4.0755, 0.6155, 0.2697, 0.0424,
                                    0.0176, 0.4188, 0.0684, 0.0284]
    _one_flipped_byte_is_refused(tmp_path, bonus.load_bonus, bonus.BONUS,
                                 b"1.8377", b"1.8378")


def test_the_recorded_original_digest_is_checked_on_its_own(tmp_path,
                                                            monkeypatch):
    """A file whose recorded ORIGINAL digest moved is refused even when its
    own sha256 is (re)accepted."""
    raw = open(bonus.BONUS, "rb").read()
    old = bonus.ORIGINAL_TXT_SHA256.encode()
    assert raw.count(old) == 1
    bad = raw.replace(old, old[:-1] + (b"1" if old[-1:] != b"1" else b"2"))
    p = tmp_path / "bonus_bad_original.json"
    p.write_bytes(bad)
    monkeypatch.setattr(bonus, "BONUS_SHA256", hashlib.sha256(bad).hexdigest())
    with pytest.raises(RuntimeError, match="records original sha256"):
        bonus.load_bonus(str(p))


def test_history_files_are_the_vendored_ones(tmp_path):
    """The 2026-09-23 query and the Deeps table the harness still reads (its
    HISTORY census): the same guards as test_bench_spin_deuteron.py's."""
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
    f = bonus.findings()
    assert f["bonus"]["quantities"] == ["F2n/F2p"]
    assert (f["deeps_blocks"], f["deeps_blocks_paper_q2"],
            f["deeps_blocks_unverified_q2"], f["deeps_duplicate_label_triples"],
            f["deeps_empty_blocks"], f["deeps_zero_value_zero_stat_rows"]) == (
        115, 60, 55, 25, [78], 60)


# ------------------------------------------------------------ the blocks

def test_census_is_the_header_s():
    c = bonus.census(bonus.load_bonus())
    assert c["blocks"] == 100 and c["blocks_by_beam"] == {"4GeV": 60, "5GeV": 40}
    assert c["rows"] == 1000 and c["missing_rows"] == 22
    four = [("4GeV", 3.38, 2.44, ps) for ps in (0.078, 0.093, 0.11, 0.135)]
    assert c["blocks_with_missing"] == four
    assert c["dis_window"] == 32
    assert c["verdict_blocks"] == 28
    assert c["verdict_blocks_by_beam"] == {"4GeV": 12, "5GeV": 16}
    assert c["dis_no_backward"] == four
    assert (c["out_of_model"], c["out_of_model_q2_below"],
            c["out_of_model_resonance"]) == (68, 20, 48)
    assert c["f2n_ratio_tables"] == {"W > 1.4 GeV": 13, "W > 1.6 GeV": 12,
                                     "W > 1.8 GeV": 10}


def test_window_is_on_the_labels_and_missing_rows_are_dropped():
    bl = {(b["beam"], b["q2"], b["w_star"], b["ps"]): b
          for b in bonus.blocks(bonus.load_bonus())}
    assert bonus.in_dis_window(bl[("5GeV", 1.66, 2.03, 0.078)])
    assert not bonus.in_dis_window(bl[("5GeV", 3.38, 1.73, 0.078)])
    assert not bonus.in_dis_window(bl[("4GeV", 0.93, 2.44, 0.078)])
    b = bl[("4GeV", 3.38, 2.44, 0.078)]
    assert b["n_missing"] == 5 and [r[0] for r in b["rows"]] == [
        0.1, 0.3, 0.5, 0.7, 0.9]
    assert bonus.backward(b) == []
    full = bl[("5GeV", 1.66, 2.44, 0.093)]
    assert [r[0] for r in bonus.backward(full)] == [-0.9, -0.7, -0.5, -0.3]
    assert [r[0] for r in bonus.forward(full)] == [0.1, 0.3, 0.5, 0.7, 0.9]


def test_a_duplicate_block_is_refused():
    d = json.loads(open(bonus.BONUS).read())
    sec = d["ratio_blocks"]["sections"][1]
    sec["blocks"].append(dict(sec["blocks"][0]))
    with pytest.raises(RuntimeError, match="duplicate"):
        bonus.blocks(d)


# ------------------------------------------------------------ statistics

def test_flat_band_by_hand():
    # chi2 below ndf: S = 1, delta = sigma / Rbar
    b = bonus.flat_band([1.0, 1.02, 0.98, 1.0], [0.02] * 4)
    assert b["rbar"] == pytest.approx(1.0)
    assert b["sigma"] == pytest.approx(0.01)
    assert b["chi2"] == pytest.approx(2.0) and b["ndf"] == 3
    assert b["scale_factor"] == 1.0 and b["delta"] == pytest.approx(0.01)
    # chi2 above ndf: the band widens by sqrt(chi2 / ndf)
    b = bonus.flat_band([1.1, 1.0, 1.0, 0.9], [0.02] * 4)
    assert b["chi2"] == pytest.approx(50.0)
    assert b["scale_factor"] == pytest.approx(math.sqrt(50.0 / 3.0))
    assert b["delta"] == pytest.approx(math.sqrt(50.0 / 3.0) * 0.01)


def test_scale_fit_absorbs_a_common_factor():
    v, e, m = [1.1, 1.0, 0.95, 0.9], [0.02, 0.03, 0.02, 0.04], [1.05, 1.0, 0.97, 0.9]
    s1, c1, n1 = bonus.scale_fit(v, e, m)
    s2, c2, n2 = bonus.scale_fit(v, e, [7.0 * x for x in m])
    assert c2 == pytest.approx(c1, rel=1e-12) and s2 == pytest.approx(s1 / 7.0)
    assert n1 == n2 == 3


def test_chi2_sf_fallback_is_t1s():
    f = _load("t1_nmc_f2d").chi2_sf_fallback
    assert bonus.chi2_sf_fallback.__code__.co_code == f.__code__.co_code
    for chi2, n in ((40.5, 10), (8.3, 9), (1.7, 4), (739.4, 84)):
        assert bonus.chi2_sf_fallback(chi2, n) == pytest.approx(f(chi2, n),
                                                               rel=1e-12)


# ------------------------------------------------------------ the tree

def test_weight_is_live_symmetric_and_reads_only_k_and_abs_cos(tree):
    w = lambda ps, c, ww=2.44, q2=1.66: tree.weight(       # noqa: E731
        bonus.DEFAULT_CONFIG, ps, c, ww, q2)
    # a live absorption, not the channel guard's silent 1
    assert w(0.135, -0.3) < 0.7 and w(0.078, -0.9) < 0.9
    for ps in (0.078, 0.11, 0.135):
        for c in (0.1, 0.5, 0.9):
            assert w(ps, c) == w(ps, -c)                 # exactly symmetric
        # the ramp is off: W*, Q2 (and x) do not reach the weight
        assert w(ps, -0.5, 2.03, 3.38) == w(ps, -0.5, 2.44, 1.66)
    # more absorption toward 90 deg, deeper at larger p_s
    assert w(0.135, -0.9) > w(0.135, -0.3) > w(0.135, -0.1)
    assert w(0.135, -0.9) / w(0.135, -0.3) > w(0.078, -0.9) / w(0.078, -0.3)
    assert tree.fsi["hulthen-20mb"].sigma_eff_mb(2.44) == 20.0


def test_weight_variation_is_converged_in_its_table(tree):
    opt = lg.GlauberFsiOptions()
    opt.n_kz = opt.n_kt = 321
    opt.k_max = 0.6
    fine = lg.GlauberFsiWeight(tree.models["hulthen"], opt)
    cs = (-0.9, -0.7, -0.5, -0.3)

    def v(f, ps):
        ws = [f.weight(lg.FsiKinematics(
            k=ps, cos_theta_k=c, phi_k=0.0, w=2.44, q2=1.66, x=0.3,
            spectator_z=1, spectator_a=1,
            channel=lg.Channel.TaggedDeuteronP)) for c in cs]
        mu = sum(ws) / len(ws)
        return max(abs(x / mu - 1.0) for x in ws)
    for ps in (0.078, 0.093, 0.11, 0.135):
        assert abs(v(fine, ps) - v(tree.fsi[bonus.DEFAULT_CONFIG], ps)) < 5e-4


def test_bin_averaging_moves_the_variation_by_under_two_percent(tree):
    """The weight is read at the bin labels; averaged uniformly over each
    p_s bin (70-85, 85-100, 100-120, 120-150 MeV/c) and each cos bin
    (width 0.2), its backward variation v moves by < 2 % of itself."""
    bins = {0.078: (0.070, 0.085), 0.093: (0.085, 0.100),
            0.11: (0.100, 0.120), 0.135: (0.120, 0.150)}
    cs = (-0.9, -0.7, -0.5, -0.3)

    def v(ws):
        mu = sum(ws) / len(ws)
        return max(abs(x / mu - 1.0) for x in ws)
    w = lambda ps, c: tree.weight(bonus.DEFAULT_CONFIG, ps, c, 2.44, 1.66)  # noqa: E731
    for lab, (lo, hi) in bins.items():
        at_label = v([w(lab, c) for c in cs])
        avg = v([sum(w(lo + (hi - lo) * (i + 0.5) / 10, c - 0.1 + 0.02 * (j + 0.5))
                     for i in range(10) for j in range(10)) / 100.0 for c in cs])
        assert abs(avg / at_label - 1.0) < 0.02, (lab, at_label, avg)


def test_density_is_the_trees_own_n_of_k(tree):
    """sum_L (a_L psi_L)^2 / 4 pi equals the M-averaged n_of_kc table on the
    grid (n_of_kc reads the cell below, so query mid-cell)."""
    for wave in ("hulthen", "av18"):
        mod = tree.models[wave]
        ks = list(mod.k)
        for i in (10, 25, 60, 120, 200):
            mid = 0.5 * (ks[i] + ks[i + 1])
            table = sum(mod.n_of_kc(m, mid, 0.3) for m in (1.0, 0.0, -1.0)) / 3
            assert tree.density(wave, ks[i]) == pytest.approx(table, rel=1e-10)


def test_light_cone_map_is_the_deeps_rows(tree):
    deeps = _load("t2_deeps_spectator_tail").Tree()
    for ps in (0.078, 0.135, 0.30):
        for c in (-0.9, -0.3, 0.0, 0.7):
            assert tree.lc(ps, c) == pytest.approx(deeps.lc(ps, c), rel=1e-14)
    # backward: alpha_s > 1 and k_LC > p_s; forward k_LC < p_s
    a, pt, k, _ = tree.lc(0.135, -0.9)
    assert a > 1.0 and k > 0.135
    assert tree.lc(0.135, 0.9)[2] < 0.135


def test_toyf2_hook_is_the_shipped_formula(tree):
    for x in (0.177, 0.3235, 0.6155, 0.95):
        assert tree.toy.f2n_over_f2p(x) == pytest.approx(
            min(1.0, max(0.25, 1.0 - 0.75 * x)), rel=1e-15)


# ------------------------------------------------------------ the row

def test_report_row_shape(rep):
    for k in ("name", "generator", "reference", "tolerance", "status",
              "subrows", "config", "history", "findings"):
        assert k in rep, k
    assert rep["name"] == NAME and rep["status"] in ("pass", "fail", "blocked")
    assert [k for k, s in rep["subrows"].items() if s["in_verdict"]] == [
        "4GeV/hulthen-40mb", "5GeV/hulthen-40mb"]
    assert rep["history"]["blocked_2026_09_23"] == bonus.BLOCKED_REASON


def test_verdict_follows_from_the_numbers(rep):
    m = rep["measured"]
    for name, judged in list(m["judged"].items()) + [
            ("stat+syst", m["judged_stat_syst"])]:
        for j in judged:
            s = max(1.0, math.sqrt(j["chi2_flat"] / j["ndf"]))
            assert j["scale_factor"] == pytest.approx(s, rel=1e-12)
            assert j["delta"] == pytest.approx(s * j["sigma"] / j["rbar"],
                                               rel=1e-12)
            assert j["passed"] == (j["v"] <= j["delta"]), (name, j)
    for key, sub in rep["subrows"].items():
        if "n_blocks" not in sub:
            assert sub["status"] == "recorded" and not sub["in_verdict"]
            continue
        assert sub["status"] == ("pass" if sub["n_pass"] == sub["n_blocks"]
                                 else "fail"), key
    heads = [rep["subrows"]["%s/hulthen-40mb" % b] for b in ("4GeV", "5GeV")]
    assert rep["status"] == ("pass" if all(h["status"] == "pass"
                                           for h in heads) else "fail")


#: what the harness measured on 2026-10-26 (this tree): per beam sub-row,
#: (n_pass, n_blocks, sum chi2_flat, sum chi2 against the weight's shape)
MEASURED = {
    "4GeV/hulthen-40mb": (4, 12, 205.71, 109.74),
    "5GeV/hulthen-40mb": (0, 16, 533.66, 116.59),
    "4GeV/av18-40mb": (3, 12, 205.71, 112.43),
    "5GeV/av18-40mb": (0, 16, 533.66, 109.19),
    "4GeV/hulthen-20mb": (8, 12, 205.71, 129.18),
    "5GeV/hulthen-20mb": (13, 16, 533.66, 274.12),
    "4GeV/hulthen-40mb stat+syst": (8, 12, 23.49, 19.02),
    "5GeV/hulthen-40mb stat+syst": (8, 16, 33.89, 9.13),
}


def test_measured_values_pinned_a_recorded_fail(rep):
    assert rep["status"] == "fail"
    sub = rep["subrows"]
    for key, (n_pass, n, chi2_flat, chi2_tree) in MEASURED.items():
        s = sub[key]
        assert (s["n_pass"], s["n_blocks"]) == (n_pass, n), key
        assert s["chi2_flat_sum"] == pytest.approx(chi2_flat, abs=0.01), key
        assert s["chi2_tree_sum"] == pytest.approx(chi2_tree, abs=0.01), key
        assert s["ndf_sum"] == 3 * n
    for b in ("4GeV", "5GeV"):
        s = sub["%s/hulthen-40mb" % b]
        assert s["ps"] == [0.078, 0.093, 0.11, 0.135]
        assert s["v_by_ps"] == pytest.approx([0.0344, 0.0548, 0.0848, 0.1406],
                                             abs=1e-4)
    assert sub["4GeV/hulthen-40mb"]["min_delta_by_ps"] == pytest.approx(
        [0.0138, 0.0255, 0.0506, 0.0536], abs=1e-4)
    assert sub["5GeV/hulthen-40mb"]["min_delta_by_ps"] == pytest.approx(
        [0.0100, 0.0257, 0.0505, 0.0672], abs=1e-4)
    assert sub["4GeV/hulthen-40mb"]["worst_v_over_delta"] == pytest.approx(
        2.622, abs=1e-3)
    assert sub["4GeV/hulthen-40mb"]["worst_block"] == "Q2=1.66,W*=2.44,ps=0.135"
    assert sub["5GeV/hulthen-40mb"]["worst_v_over_delta"] == pytest.approx(
        3.441, abs=1e-3)
    assert sub["5GeV/hulthen-40mb"]["worst_block"] == "Q2=1.66,W*=2.03,ps=0.078"
    # the four passing 4 GeV blocks are the Q2 = 3.38, W* = 2.03 ones
    passed = [(j["beam"], j["q2"], j["w_star"])
              for j in rep["measured"]["judged"]["hulthen-40mb"] if j["passed"]]
    assert passed == [("4GeV", 3.38, 2.03)] * 4


def test_forward_rise_recorded_not_asserted(rep):
    sub = rep["subrows"]["forward rise (recorded)"]
    assert sub["status"] == "recorded" and sub["in_verdict"] is False
    fwd = rep["measured"]["forward"]
    assert list(fwd) == ["%s/ps=%g" % (b, p) for b in ("4GeV", "5GeV")
                         for p in (0.078, 0.093, 0.11, 0.135)]
    assert [v["rho_data"] for v in fwd.values()] == pytest.approx(
        [1.036, 0.973, 1.020, 1.026, 1.031, 0.994, 0.974, 1.030], abs=1e-3)
    for v, rt in zip(fwd.values(), [0.9902, 0.9838, 0.9741, 0.9552] * 2):
        assert v["rho_tree_min"] == pytest.approx(rt, abs=1e-4)
        assert v["rho_tree_max"] == pytest.approx(rt, abs=1e-4)
    assert sum(v["same_sign_blocks"] for v in fwd.values()) == 12
    assert sum(v["n"] for v in fwd.values()) == 28


def test_prescription_and_f2n_f2p_recorded(rep):
    pr = rep["subrows"]["PWIA prescription instant/LC (recorded, not FSI)"]
    assert pr["status"] == "recorded" and pr["in_verdict"] is False
    assert pr["max_dev_by_ps"] == pytest.approx([0.0227, 0.0309, 0.0409, 0.0571],
                                                abs=1e-4)
    assert pr["ends_ratio_by_ps"] == pytest.approx(
        [1.0457, 1.0626, 1.0835, 1.1177], abs=1e-4)
    assert pr["weight_ends_ratio_by_ps"] == pytest.approx(
        [1.0710, 1.1140, 1.1790, 1.3103], abs=1e-4)
    f2 = rep["measured"]["f2n_f2p"]
    assert f2["all"]["n"] == 10 and f2["w_ge_2"]["n"] == 5
    assert f2["all"]["chi2_abs"] == pytest.approx(40.5, abs=0.05)
    assert f2["all"]["chi2_scaled"] == pytest.approx(8.30, abs=0.01)
    assert f2["all"]["scale"] == pytest.approx(0.9020, abs=1e-4)
    assert f2["w_ge_2"]["chi2_abs"] == pytest.approx(35.5, abs=0.05)
    assert f2["w_ge_2"]["chi2_scaled"] == pytest.approx(1.70, abs=0.01)
    assert f2["all"]["tree_over_data"][0] == pytest.approx(1.190, abs=1e-3)
    assert f2["all"]["tree_over_data"][-1] == pytest.approx(1.286, abs=1e-3)
    sub = rep["subrows"]["F2n/F2p ToyF2 vs W*>1.8 table (recorded)"]
    assert sub["status"] == "recorded" and sub["in_verdict"] is False


def test_header_declares_the_rule_the_code_applies():
    doc = bonus.__doc__
    assert (bonus.W_STAR_MIN, bonus.Q2_MIN, bonus.BACKWARD_COS_MAX,
            bonus.FORWARD_COS_GT) == (2.0, 1.0, -0.3, 0.0)
    for text in ("W* >= 2 GeV AND Q2 >= 1 GeV^2", "v <= delta",
                 "BEFORE any fit, band, chi2", "previewed",
                 "28 VERDICT BLOCKS: 12 at 4 GeV, 16 at 5 GeV",
                 "CORRECTED 2026-10-26", "THE FRAME CAVEAT",
                 "WHAT THIS ROW DOES NOT VALIDATE"):
        assert text in doc, text
    assert doc.index("THE RULE AND ITS TOLERANCE") < doc.index("MEASURED (2026-10-26")


def test_report_normalises_for_make_report(rep):
    sys.path.insert(0, _BENCH)
    import _report
    row = _report.normalise_row(rep, NAME)
    assert row["status"] == "fail" and not row["notes"]
    names = [s["name"] for s in row["subrows"]]
    assert names[:2] == ["4GeV/hulthen-40mb", "5GeV/hulthen-40mb"]
    assert len(names) == 12
    verdict = [s["name"] for s in row["subrows"] if s.get("in_verdict")]
    assert verdict == names[:2]
    line = _report.report_line(rep["name"], rep["generator"], rep["reference"],
                               rep["tolerance"], rep["status"])
    assert line.startswith("REPORT | %s | " % NAME) and line.endswith("| FAIL")


# ------------------------------------------------------------ exit codes

def test_harness_exits_0_and_prints_one_report_row():
    p = _run_script(os.path.join(_BENCH, NAME + ".py"))
    assert p.returncode == 0, p.stderr.decode("utf-8", "replace")
    out = p.stdout.decode("utf-8")
    rows = [l for l in out.splitlines() if l.startswith("REPORT | ")]
    assert len(rows) == 1
    assert rows[0].startswith("REPORT | %s | " % NAME) and rows[0].endswith("| FAIL")


def test_harness_prints_its_row_on_a_non_utf8_stdout():
    """The history's catalogue names carry Greek letters: an ASCII stdout
    gets them escaped, and the REPORT row and exit 0 regardless."""
    p = _run_script(os.path.join(_BENCH, NAME + ".py"), PYTHONIOENCODING="ascii")
    assert p.returncode == 0, p.stderr.decode("ascii", "replace")
    assert b"REPORT | t2_bonus_spectator_shape |" in p.stdout
    assert b"\\u03b7" in p.stdout


def _copy_harness(tmp_path, altered=None, drop=None):
    d = tmp_path / "bench"
    (d / "data").mkdir(parents=True)
    shutil.copy(os.path.join(_BENCH, NAME + ".py"), d / (NAME + ".py"))
    for path in (bonus.BONUS, bonus.QUERY, bonus.DEEPS):
        f = os.path.basename(path)
        if f == drop:
            continue
        raw = open(path, "rb").read()
        if f == altered:
            raw = raw.replace(b"1.8377", b"1.8378", 1)
        (d / "data" / f).write_bytes(raw)
    return str(d / (NAME + ".py"))


def test_harness_exits_2_on_an_altered_or_missing_file(tmp_path):
    f = os.path.basename(bonus.BONUS)
    p = _run_script(_copy_harness(tmp_path / "a", altered=f))
    assert p.returncode == 2
    assert b"sha256" in p.stderr and b"REPORT |" not in p.stdout
    p = _run_script(_copy_harness(tmp_path / "m", drop=f))
    assert p.returncode == 2 and b"REPORT |" not in p.stdout
    # the copy itself runs: the refusal above is the file, not the copy
    p = _run_script(_copy_harness(tmp_path / "ok"))
    assert p.returncode == 0 and b"REPORT | t2_bonus_spectator_shape |" in p.stdout


# ------------------------------------------------------------ the figure

def test_plot_writes_a_png(rep, tmp_path):
    pytest.importorskip("matplotlib")
    path = bonus.plot(rep, out_dir=str(tmp_path))
    assert path == str(tmp_path / (NAME + ".png"))
    raw = open(path, "rb").read()
    assert len(raw) > 10000 and raw[:8] == b"\x89PNG\r\n\x1a\n"
    fs = rep["measured"]["figure"]
    assert list(fs["blocks"]) == ["4GeV/W*=2.03", "4GeV/W*=2.44",
                                  "5GeV/W*=2.03", "5GeV/W*=2.44"]
    assert list(fs["lines"]) == list(bonus.CONFIGS)
