# SPDX-License-Identifier: GPL-3.0-or-later
"""validation/benchmarks/t3_li6_ff_wiringa_schiavilla.py exercised as tests --
docs/benchmarking/04_theory.md T-20, Wiringa & Schiavilla, PRL 81 (1998) 4317,
Fig. 1 read back from the PDF's vector drawing (run 2026-10-26 phase C4;
< 1 s per harness run, so not behind LIPOLGEN_BENCH).

The harness holds the shipped HoSpin1FF C0 zero (and the vmc-ft edge) to the
one bound the figure gives -- WS98's IA+MEC C0 is drawn non-zero up to
q_last, so the zero lies above it -- under the plan's tolerance (delta +
0.1 fm^-1), declared in a hashed design block; it records the rest (zero
extrapolations, the slope radius, the C0 shape, F_M1 against fm, F_C2 against
fq, and closures of the reading against the paper's printed mu and Q, the
figure's own line/symbol sum and RFY Table I).  These tests pin what the
harness MEASURES and that its verdict follows from the numbers: the closures
that do not hold are pinned as recorded FAILs; nothing here moves or asserts
a new value for any tree default.  Nothing needs the network.  The one test
that re-reads the PDF (`--redigitize`) skips unless PyMuPDF and the PDF are
present ($LIPOLGEN_WS98_PDF, else the sibling corpus's refs/).
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
_NAME = "t3_li6_ff_wiringa_schiavilla"
_DATA = "ws1998_li6_formfactors.csv"
_RFY = "rfy1966_li_magnetic_ff.csv"


def _load(name=_NAME):
    path = os.path.join(_BENCH, name + ".py")
    spec = importlib.util.spec_from_file_location("bench_" + name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _flip(path, old, new, tmp_path):
    """A copy of `path` with the one occurrence of `old` replaced."""
    raw = open(path, "rb").read()
    assert raw.count(old) == 1, old
    bad = tmp_path / os.path.basename(path)
    bad.write_bytes(raw.replace(old, new))
    return str(bad)


def _env(extra=None):
    env = dict(os.environ)
    pkg_parent = os.path.dirname(os.path.dirname(os.path.abspath(lg.__file__)))
    env["PYTHONPATH"] = os.pathsep.join(
        [pkg_parent] + [p for p in env.get("PYTHONPATH", "").split(os.pathsep) if p])
    env.update(extra or {})
    return env


def _run_harness(script, args=()):
    """The harness as a script, with the imported lipolgen on its path."""
    return subprocess.run([sys.executable, script] + list(args), capture_output=True,
                          text=True, env=_env(), timeout=300)


def _copy_bench(tmp_path, data_name=None, old=None, new=None, harness_edit=None):
    """The harness and both files it reads in a fresh directory, one of them
    (or the harness's own text) altered by one replacement."""
    d = tmp_path / "bench"
    (d / "data").mkdir(parents=True)
    src = open(os.path.join(_BENCH, _NAME + ".py"), encoding="utf-8").read()
    if harness_edit is not None:
        assert src.count(harness_edit[0]) == 1
        src = src.replace(*harness_edit)
    (d / (_NAME + ".py")).write_text(src, encoding="utf-8")
    for name in (_DATA, _RFY):
        raw = open(os.path.join(_BENCH, "data", name), "rb").read()
        if name == data_name:
            assert raw.count(old) == 1, old
            raw = raw.replace(old, new)
        (d / "data" / name).write_bytes(raw)
    return str(d / (_NAME + ".py"))


@pytest.fixture(scope="module")
def mod():
    return _load()


@pytest.fixture(scope="module")
def rep(mod):
    return mod.run(verbose=False)


# ------------------------------------------------------------ the vendored file

def test_vendored_file_guards(mod, tmp_path):
    """Point counts, the first and last values, the provenance, and one
    flipped byte refused -- in the body and in the provenance block."""
    ref = mod.load_reference()
    rec, slot = ref["record"], ref["slot"]
    assert (len(rec), len(slot)) == (503, 68)
    counts = {}
    for r in rec:
        k = (r["curve"], r["calc"], r["kind"])
        counts[k] = counts.get(k, 0) + 1
    counts[("slot", "IA+MEC", "slot")] = len(slot)
    assert counts == mod.EXPECTED_COUNTS
    first, last = rec[0], rec[-1]
    assert (first["curve"], first["calc"], first["kind"], first["q"], first["F2"],
            first["x"], first["y"]) == ("C0", "IA", "vertex", 0.0, 1.0, 104.40, 214.56)
    assert (last["curve"], last["calc"], last["q"], last["F2"], last["F2_lo"],
            last["F2_hi"]) == ("FT2", "IA+MEC", 4.0, 1.35154e-08, 1.35356e-10, 4.35112e-08)
    assert (slot[0]["q"], slot[0]["F_C0"], slot[0]["F_C2"], slot[0]["F_M1"]) == (
        0.119048, 2.93801, -0.00161562, 4.9768)
    assert (slot[-1]["q"], slot[-1]["F_C0"], slot[-1]["F_C2"], slot[-1]["F_M1"]) == (
        2.919048, 0.00109798, 0.00259789, -0.0502486)
    raw = open(mod.DATA, encoding="utf-8").read()
    assert raw.startswith("# SPDX-License-Identifier: GPL-3.0-or-later\n")
    for needle in ("fc51cb930e81b16238bb02561e81d4f9160a7d1708e813f8380ba669991d20e8",
                   "Phys. Rev. Lett. 81 (1998)", "(104.40, 214.56, 306.00, 585.12)",
                   "--redigitize", "tools/digitize_figure.py",
                   "PUBLISHING CONVENTION (per nucleus"):
        assert needle in raw
    assert hashlib.sha256(raw.encode("utf-8")).hexdigest() == mod.DATA_SHA256
    bad = _flip(mod.DATA, b"2.919048,0.00109798,", b"2.919048,0.00109799,", tmp_path)
    with pytest.raises(RuntimeError, match="sha256"):
        mod.load_reference(bad)
    bad = _flip(mod.DATA, b"# Licence:   facts", b"# Licence:   Facts", tmp_path)
    with pytest.raises(RuntimeError, match="sha256"):
        mod.load_reference(bad)


def test_every_row_is_rederived_from_its_coordinates(mod, tmp_path, monkeypatch):
    """Even with its sha256 re-recorded, an edited number is refused: q and
    F2 come back from (x_pt, y_pt) through the calibration, and every slot
    row from the record rows by the stated rule."""
    def rehashed(old, new):
        bad = _flip(mod.DATA, old, new, tmp_path)
        monkeypatch.setattr(mod, "DATA_SHA256",
                            hashlib.sha256(open(bad, "rb").read()).hexdigest())
        return bad
    bad = rehashed(b"C0,IA+MEC,vertex,90,0,251.52,578.40,,,1.3395e-07,",
                   b"C0,IA+MEC,vertex,90,0,251.52,578.40,,,1.4395e-07,")
    with pytest.raises(RuntimeError, match="calibrated"):
        mod.load_reference(bad)
    bad = rehashed(b"2.919048,0.00109798,", b"2.919048,0.00119798,")
    with pytest.raises(RuntimeError, match="rule"):
        mod.load_reference(bad)


def test_calibration_and_device_unit(mod):
    assert mod.to_q(104.40, "left") == 0.0 and mod.to_q(306.00, "left") == 4.0
    assert mod.to_q(306.00, "right") == 0.0 and mod.to_q(507.60, "right") == 4.0
    assert mod.to_log10f2(214.56, "left") == 0.0
    assert mod.to_log10f2(585.12, "left") == -7.0
    assert mod.to_log10f2(214.56, "right") == -3.0
    assert mod.to_log10f2(585.12, "right") == -10.0
    # the major tick marks of the left panel (decades) read back within 0.0019
    for y, dec in ((532.08, -6), (479.28, -5), (426.24, -4), (373.44, -3),
                   (320.40, -2), (267.60, -1)):
        assert abs(mod.to_log10f2(y, "left") - dec) < 0.002
    assert mod.DQ_FM == pytest.approx(0.24 / 50.4, rel=1e-12)
    assert mod.DLOG10F2 == pytest.approx(0.24 / (370.56 / 7.0), rel=1e-12)


def test_promotion_is_one_copy(mod, tmp_path):
    """data/ff/li6_elastic.csv's reader, TabulatedSpin1FF, run on a COPY of
    the vendored file in a scratch LIPOLGEN_DATA_DIR (data/ff/ in the
    repository stays empty: tests/test_rc.cpp asserts no 6Li table ships).
    It must skip every record row and read exactly the 68 slot rows."""
    ff = tmp_path / "ff"
    ff.mkdir()
    shutil.copy(mod.DATA, ff / "li6_elastic.csv")
    code = r'''
import json, math, lipolgen as lg
norm = lg.HoSpin1FFOptions(); norm.m_a_gev = lg.li6().mass()
t = lg.TabulatedSpin1FF.from_data_dir("ff/li6_elastic.csv", norm)
hc = lg.HBARC_GEV_FM
out = dict(prov=t.provenance(), nodes=[], refuse=[])
for q in json.loads(input()):
    tt = (q * hc) ** 2
    out["nodes"].append([t.fc(tt), t.fm(tt), t.fq(tt), tt / (4 * norm.m_a_gev ** 2)])
for q in (0.119048 * 0.999, 2.919048 * 1.001):
    try:
        t.fc((q * hc) ** 2); out["refuse"].append(False)
    except Exception as e:
        out["refuse"].append("refuses to extrapolate" in str(e))
print(json.dumps(out))
'''
    slot = mod.load_reference()["slot"]
    inner = slot[1:-1]
    res = subprocess.run([sys.executable, "-c", code], input=json.dumps([s["q"] for s in inner]),
                         capture_output=True, text=True, timeout=120,
                         env=_env({"LIPOLGEN_DATA_DIR": str(tmp_path)}))
    assert res.returncode == 0, res.stderr
    out = json.loads(res.stdout)
    assert "68 rows" in out["prov"] and "0.119" in out["prov"] and "2.919" in out["prov"]
    for s, (fc, fm, fq, eta) in zip(inner, out["nodes"]):
        assert fc == pytest.approx(s["F_C0"], rel=1e-9, abs=1e-15)
        assert fm == pytest.approx(s["F_M1"], rel=1e-9, abs=1e-15)
        assert fq == pytest.approx(3.0 * s["F_C2"] / (2.0 * math.sqrt(2.0) * eta), rel=1e-9)
    assert out["refuse"] == [True, True]
    assert not os.path.exists(os.path.join(_ROOT, "data", "ff", "li6_elastic.csv"))


def test_slot_rule(mod):
    """The slot columns are the tree's normalisation: F_C0 = Z |F_C0|; F_C2
    negative before the C2 floor gap and positive after; F_M1 the Sachs F_m,
    + below the F_T zero and - above, near (M_A/m_p) mu_WS at q -> 0."""
    slot = mod.load_reference()["slot"]
    assert all(s["F_C0"] > 0 for s in slot)
    assert all(s["F_C2"] < 0 for s in slot if s["q"] < 2.53)
    assert all(s["F_C2"] > 0 for s in slot if s["q"] > 2.67)
    sign = [s["F_M1"] > 0 for s in slot]
    flip = sign.index(False)
    assert all(sign[:flip]) and not any(sign[flip:])
    assert 1.27 < slot[flip - 1]["q"] < slot[flip]["q"] < 1.33
    ratio = lg.li6().mass() / lg.PROTON_MASS
    assert slot[0]["F_M1"] / ratio == pytest.approx(0.8336, abs=5e-4)
    rec = mod.load_reference()["record"]
    c0 = [r for r in mod.pick(rec, "C0", "IA+MEC") if r["x"] == slot[0]["x"]][0]
    assert slot[0]["F_C0"] == pytest.approx(3.0 * math.sqrt(c0["F2"]), rel=2e-6)


# ------------------------------------------------------------ the design block

def test_design_block_is_hashed_and_states_the_rules(mod):
    block = mod.design_block()
    assert hashlib.sha256(block.encode("utf-8")).hexdigest() == mod.DESIGN_BLOCK_SHA256
    assert mod.design_block_sha256() == mod.DESIGN_BLOCK_SHA256
    rules = mod.rules_from_design_block(block)
    assert rules["ZERO_MARGIN_FM"] == 0.1 and rules["DELTA_FM"] == 0.0047619
    assert rules["MU_TOL"] == 0.02 and rules["Q_TOL"] == 0.05
    assert rules["SLOPE_QMAX_VARIANTS"] == (0.6, 1.0, 1.4)
    assert "THE UPPER SIDE IS NOT\nTESTED" in block
    assert "WHAT WAS SEEN BEFORE THIS BLOCK WAS WRITTEN" in block
    src = open(os.path.join(_BENCH, _NAME + ".py"), encoding="utf-8").read()
    with pytest.raises(RuntimeError, match="design block"):
        mod.design_block_sha256(src.replace("(delta + 0.1 fm^-1)", "(delta + 0.2 fm^-1)", 1))


def test_rule_constants_must_match_the_block(mod, monkeypatch):
    monkeypatch.setattr(mod, "ZERO_MARGIN_FM", 0.2)
    with pytest.raises(RuntimeError, match="ZERO_MARGIN_FM"):
        mod.design_block_sha256()


# ------------------------------------------------------------ the verdict

def test_verdict_is_pinned(rep):
    assert rep["status"] == "pass"
    v = rep["verdict"]
    assert rep["q_last"] == {"IA+MEC": 2.919048, "IA": 3.2}
    assert rep["q0"]["ho"] == pytest.approx(3.0997711, abs=2e-7)
    assert rep["q0"]["vmc-ft"] is None
    assert v["ho"]["bound"] == pytest.approx(2.919048 - 0.1047619, abs=1e-9)
    assert v["ho"]["distance_above_bound"] == pytest.approx(0.285485, abs=2e-6)
    assert (v["ho"]["status"], v["vmc-ft"]["status"]) == ("pass", "pass")
    assert "no C0 zero below 5" in v["vmc-ft"]["reason"]
    # against the IA line (recorded): the shipped zero meets the same rule
    # by 0.0046 fm^-1 only
    w = rep["verdict_ia"]["ho"]
    assert w["bound"] == pytest.approx(3.2 - 0.1047619, abs=1e-9)
    assert w["distance_above_bound"] == pytest.approx(0.004533, abs=2e-6)


def test_verdict_follows_the_numbers(mod):
    q_last = 2.919048
    bound = q_last - mod.DELTA_FM - mod.ZERO_MARGIN_FM
    assert mod.zero_verdict(bound + 1e-9, q_last)["status"] == "pass"
    assert mod.zero_verdict(bound - 1e-6, q_last)["status"] == "fail"
    # row 5's data band [2.6944, 2.8284] fm^-1 straddles the bound: a tree
    # zero refitted to the FB zero would contradict WS98's drawn C0, one at
    # Li71's minimum would not (by 0.014 fm^-1)
    assert mod.zero_verdict(2.6944, q_last)["status"] == "fail"
    assert mod.zero_verdict(2.8284, q_last)["status"] == "pass"
    assert mod.zero_verdict(2.8284, q_last)["distance_above_bound"] == pytest.approx(
        0.014114, abs=2e-6)
    assert mod.zero_verdict(None, q_last)["status"] == "pass"
    # the tree's zero is row 5's first_zero of the same fc
    e = mod.Edge(lg.HoSpin1FF.for_ion(lg.li6()), mod.tree_side())
    assert mod.first_zero(e.c0) == pytest.approx(3.0997711, abs=2e-7)


# ------------------------------------------------------------ recorded items

def test_recorded_zero_items_are_pinned(rep):
    f = rep["extrapolation"]["IA+MEC"]["fits"]
    assert f["e1_line_last2"]["root"] == pytest.approx(3.025151, abs=2e-6)
    assert f["e2_line_last5"]["root"] == pytest.approx(2.973489, abs=2e-6)
    assert f["e3_quad_last5"]["root"] is None and f["e4_quad_last10"]["root"] is None
    g = rep["extrapolation"]["IA"]["fits"]
    assert g["e1_line_last2"]["root"] == pytest.approx(3.413816, abs=2e-6)
    c = rep["ceiling"]
    assert c["ho"]["max_f2"] == pytest.approx(6.30969e-10, rel=1e-4)
    assert c["vmc-ft"]["max_f2"] == pytest.approx(7.00349e-07, rel=1e-4)
    assert (c["ho"]["above_floor"], c["vmc-ft"]["above_floor"]) == (False, True)


def test_recorded_slope_and_shape_are_pinned(rep):
    s = rep["slope"]["IA+MEC"]
    assert s["ws"]["r"] == pytest.approx(2.569905, abs=2e-6)
    assert s["ws"]["r_err_dig"] == pytest.approx(0.009131, abs=2e-6)
    assert (s["ws"]["r_min_variants"], s["ws"]["r_max_variants"]) == pytest.approx(
        (2.524345, 2.694871), abs=2e-6)
    assert s["ho"]["r_same_fit"] == pytest.approx(2.570813, abs=2e-6)
    assert s["ho"]["r_richardson"] == pytest.approx(2.5710008, abs=2e-7)
    assert s["vmc-ft"]["r_richardson"] == pytest.approx(2.5709924, abs=2e-6)
    assert rep["slope"]["IA"]["ws"]["r"] == pytest.approx(2.557652, abs=2e-6)
    sh = rep["shape"]["ho vs IA+MEC"]
    assert (sh["chi2"], sh["ndf"], sh["s"]) == (pytest.approx(1.893953e6, rel=1e-5), 72,
                                                pytest.approx(4.797828, abs=2e-6))
    assert sh["ratio_q25"]["ratio"] == pytest.approx(10.99717, rel=1e-5)
    assert sh["ratio_qlast"]["ratio"] == pytest.approx(43.20022, rel=1e-5)
    sv = rep["shape"]["vmc-ft vs IA+MEC"]
    assert sv["s"] == pytest.approx(1.058195, abs=2e-6)
    assert (sv["ratio_q1"]["ratio"], sv["ratio_q2"]["ratio"], sv["ratio_q25"]["ratio"],
            sv["ratio_qlast"]["ratio"]) == pytest.approx((1.018712, 1.253943, 1.203200,
                                                          0.179386), abs=2e-6)


def test_recorded_m1_and_c2_are_pinned(rep):
    m = rep["m1"]
    assert (m["chi2"], m["ndf"]) == (pytest.approx(88.90743, abs=1e-4), 19)
    assert m["p"] == pytest.approx(5.1772e-11, rel=1e-3)
    assert (m["s_hat"], m["s_err"]) == pytest.approx((1.024359, 0.015620), abs=2e-6)
    assert m["chi2_s1"] == pytest.approx(91.33945, abs=1e-4)
    assert m["ws_zero_bracket"] == (1.0, 1.247619, 1.5)
    n = rep["m1_nominal_q"]
    assert (n["chi2"], n["s_hat"]) == pytest.approx((77.70252, 1.009766), abs=1e-5)
    c = rep["c2"]
    assert c["ratio_low_q_mean"] == pytest.approx(5.872197, abs=2e-6)
    assert (c["ratio_low_q_min"], c["ratio_low_q_max"]) == pytest.approx(
        (5.722994, 6.063628), abs=2e-6)
    assert c["ws_zero_bracket"] == (2.519048, 2.67619)
    assert c["tree_c2_zero_equals_c0_zero"] is True
    assert (c["ws_second_lobe_q"], c["ws_second_lobe_f2"]) == (3.238095, 1.08063e-06)


def test_closures_are_pinned(rep):
    cl = rep["closures"]
    # declared: the mu closure does NOT hold at the digitized q (a recorded
    # FAIL); the post-hoc nominal-q diagnostic does
    assert (cl["mu[IA+MEC]"]["mu"], cl["mu[IA+MEC]"]["holds"]) == (
        pytest.approx(0.895187, abs=2e-6), False)
    assert (cl["mu[IA]"]["mu"], cl["mu[IA]"]["holds"]) == (
        pytest.approx(0.879008, abs=2e-6), False)
    assert (cl["mu_nominal_q[IA+MEC]"]["mu"], cl["mu_nominal_q[IA+MEC]"]["holds"]) == (
        pytest.approx(0.841427, abs=2e-6), True)
    assert (cl["mu_nominal_q[IA]"]["mu"], cl["mu_nominal_q[IA]"]["holds"]) == (
        pytest.approx(0.826842, abs=2e-6), True)
    q = cl["Q[IA+MEC]"]
    assert (q["abs_q_fm2"], q["holds"]) == (pytest.approx(0.490144, abs=2e-6), False)
    assert q["ratio_sqrt2_over_3"] == pytest.approx(0.231056, abs=2e-6)
    ls = cl["lines_vs_symbols[IA+MEC]"]
    assert (ls["n"], ls["holds"]) == (11, True)
    assert ls["max_abs_pull"] == pytest.approx(0.542001, abs=2e-6)
    r = cl["rfy_table1[Stanford F_T^2]"]
    assert (r["q_holds"], r["constant"]) == (True, True)
    assert r["max_abs_dq"] == pytest.approx(0.004969, abs=2e-6)
    assert (r["ratio_mean"], r["ratio_min"], r["ratio_max"]) == pytest.approx(
        (0.711609, 0.709189, 0.714383), abs=2e-6)
    # the digitized bars reproduce RFY's printed relative errors
    assert [round(p["rel_err_dig"], 2) for p in r["points"]] == [
        0.30, 0.21, 0.41, 0.30, 0.60, 0.50, 0.93]


def test_rule3_per_nucleus_reading(mod, rep):
    """Both sides per nucleus: the tree's fc/Z is 1 at q = 0 like WS98's
    F_C0 (the line starts at the frame corner); F_T^2 -> (4/3) eta (1 + eta)
    [(M_A/m_p) mu]^2/Z^2 at q -> 0; the zero does not depend on the scale."""
    tree = mod.tree_side()
    e = mod.Edge(tree["edges"]["ho"], tree)
    assert e.c0(0.0) == 1.0
    q = 1e-3
    eta = q * q / (4.0 * tree["m_a_fm"] ** 2)
    fm0 = tree["m_a_gev"] / lg.PROTON_MASS * tree["options"]["mu_n"]
    assert e.ft2(q) == pytest.approx(4.0 / 3.0 * eta * (1 + eta) * fm0 ** 2 / 9.0, rel=1e-5)
    assert e.c2(q) == pytest.approx(math.sqrt(2.0) / 6.0 * q * q * (-0.0818) / 3.0, rel=1e-5)
    assert mod.first_zero(lambda x: 7.0 * e.c0(x)) == mod.first_zero(e.c0)
    rec = mod.load_reference()["record"]
    assert mod.pick(rec, "C0", "IA+MEC")[0]["F2"] == 1.0


# ------------------------------------------------------------ the row

def test_report_line(mod, capsys):
    r = mod.run(verbose=True)
    out = capsys.readouterr().out.splitlines()
    rows = [l for l in out if l.startswith("REPORT |")]
    assert len(rows) == 1
    assert rows[0].startswith(
        "REPORT | t3_li6_ff_wiringa_schiavilla | HoSpin1FF C0 zero q0 = 3.0998 fm^-1 "
        "(shipped ho); vmc-ft: no zero below 5 fm^-1 | Wiringa-Schiavilla PRL 81 "
        "(1998) 4317 Fig. 1")
    assert rows[0].endswith("| PASS")
    assert "q_last = 2.9190 fm^-1" in rows[0] and "= 2.8143 fm^-1" in rows[0]
    subs = [l for l in out if l.startswith("SUBROW |")]
    assert len(subs) == len(r["subrows"]) == 24
    assert "np.float64" not in "\n".join(out)
    assert r["runtime_s"] < 30.0


def test_exit_codes(tmp_path):
    ok = _run_harness(os.path.join(_BENCH, _NAME + ".py"))
    assert ok.returncode == 0, ok.stderr
    assert "REPORT | %s |" % _NAME in ok.stdout
    bad = _run_harness(_copy_bench(tmp_path / "a", _DATA, b"2.919048,0.00109798,",
                                   b"2.919048,0.00109799,"))
    assert bad.returncode == 2 and "sha256" in bad.stderr
    assert "REPORT |" not in bad.stdout
    bad = _run_harness(_copy_bench(tmp_path / "b", _RFY, b"6Li,I,1,865,85,0.72,0.10,0.03",
                                   b"6Li,I,1,865,85,0.72,0.11,0.03"))
    assert bad.returncode == 2 and "sha256" in bad.stderr
    bad = _run_harness(_copy_bench(tmp_path / "c", harness_edit=(
        "must not exceed delta + 0.1 fm^-1, i.e.", "must not exceed delta + 0.2 fm^-1, i.e.")))
    assert bad.returncode == 2 and "design block" in bad.stderr


def test_redigitize_reproduces_the_body(mod):
    """DEV TIME: re-read the PDF's drawing operators and compare with the
    vendored body, line for line.  Needs PyMuPDF and the arXiv PDF (not in
    this repository): $LIPOLGEN_WS98_PDF, else the sibling corpus."""
    pytest.importorskip("pymupdf", reason="PyMuPDF (dev-time re-digitization only)")
    pdf = os.environ.get("LIPOLGEN_WS98_PDF") or os.path.join(
        os.path.dirname(_ROOT), "PolarizedLithiumSim", "refs", "nucl-th_9807037.pdf")
    if not os.path.isfile(pdf):
        pytest.skip("the WS98 PDF (arXiv:nucl-th/9807037v2, sha256 %s...) is not at %s; "
                    "set LIPOLGEN_WS98_PDF" % (mod.WS_PDF_SHA256[:12], pdf))
    res = _run_harness(os.path.join(_BENCH, _NAME + ".py"), ["--redigitize", pdf])
    assert res.returncode == 0, res.stdout + res.stderr
    assert "identical (572 lines)" in res.stdout


def test_plot_writes_png(mod, rep, tmp_path):
    """plot() draws the row's figure -- a view of the report, never the
    record -- as a non-empty PNG where it is told to (here a tmp dir)."""
    pytest.importorskip("matplotlib", reason="matplotlib (the figure only)")
    path = mod.plot(rep, out_dir=str(tmp_path))
    assert path == os.path.join(str(tmp_path), _NAME + ".png")
    with open(path, "rb") as f:
        assert f.read(8) == b"\x89PNG\r\n\x1a\n"
    assert os.path.getsize(path) > 0
