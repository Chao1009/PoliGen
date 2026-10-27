# SPDX-License-Identifier: GPL-3.0-or-later
"""validation/benchmarks/t3_li_magnetization_rfy.py exercised as tests --
BENCHMARK_PLAN.md sec. 4 row 6, Rand-Frosch-Yearian, Phys. Rev. 144 (1966)
859 (run 2026-10-26 phase B2; < 1 s per harness run, so not behind
LIPOLGEN_BENCH).

The harness compares RFY Table I's seven measured 6Li F_M^2 points with the
tree's unfitted F_mag shape (one free scale, p >= 0.01, declared in a hashed
design block), records the 24 7Li points of Table II (the tree has no 7Li
F_m), re-evaluates the paper's own fits on the vendored points as a
transcription closure, and vendors a0 verbatim (the radius rule of
2026-09-23 is retired).  These tests pin what the harness MEASURES and that
its verdict follows from the numbers: the 7Li printed-parameter closure that
does not hold is pinned as a recorded FAIL; nothing here moves or asserts a
new value for any tree default.  Nothing needs the network.

Supersedes python/tests/test_bench_nuclear.py's two test_rfy_* tests (the
BLOCKED-era skip and the tree-side closed form), which the integration stage
retires or re-points; until then they still pass against this harness.
"""
import importlib.util
import math
import os
import shutil
import subprocess
import sys

import pytest

lg = pytest.importorskip("lipolgen")

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_BENCH = os.path.join(_ROOT, "validation", "benchmarks")
_NAME = "t3_li_magnetization_rfy"
_FF = "rfy1966_li_magnetic_ff.csv"
_A0 = "rfy1966_li_magnetization.csv"
_PDF_SHA256 = "7cbb2d21615d5e500a9c3235a026382e8529fc6c7828fa12a8fa6cfd31d5bbe2"


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


def _run_harness(script, args=()):
    """The harness as a script, with the imported lipolgen on its path."""
    env = dict(os.environ)
    pkg_parent = os.path.dirname(os.path.dirname(os.path.abspath(lg.__file__)))
    env["PYTHONPATH"] = os.pathsep.join(
        [pkg_parent] + [p for p in env.get("PYTHONPATH", "").split(os.pathsep) if p])
    return subprocess.run([sys.executable, script] + list(args), capture_output=True,
                          text=True, env=env, timeout=300)


def _copy_bench(tmp_path, data_name=None, old=None, new=None, harness_edit=None):
    """The harness and both vendored files in a fresh directory, one of them
    (or the harness's own text) altered by one replacement."""
    d = tmp_path / "bench"
    (d / "data").mkdir(parents=True)
    src = open(os.path.join(_BENCH, _NAME + ".py"), encoding="utf-8").read()
    if harness_edit is not None:
        assert src.count(harness_edit[0]) == 1
        src = src.replace(*harness_edit)
    (d / (_NAME + ".py")).write_text(src, encoding="utf-8")
    for name in (_FF, _A0):
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


def _sub(rep, name):
    hits = [s for s in rep["subrows"] if s["name"] == name]
    assert len(hits) == 1, name
    return hits[0]


# ------------------------------------------------------------ the vendored files

def test_ff_file_is_tables_i_and_ii(mod, tmp_path):
    li6, li7 = mod.load_ff()
    assert (len(li6), len(li7)) == (7, 24)
    first, last = li6[0], li6[-1]
    assert (first["e0_mev"], first["q2_fm2"], first["fm2"], first["fm2_err"]) == (
        85.0, 0.72, 0.10, 0.03)
    assert (last["e0_mev"], last["q2_fm2"], last["fm2"], last["fm2_err"]) == (
        140.0, 1.92, 0.001, 0.001)
    assert (li7[0]["e0_mev"], li7[0]["q2_fm2"], li7[0]["fm2"], li7[0]["fm2_err"]) == (
        70.0, 0.49, 0.43, 0.03)
    assert (li7[-1]["e0_mev"], li7[-1]["q2_fm2"], li7[-1]["fm2"], li7[-1]["fm2_err"]) == (
        200.0, 3.87, 0.016, 0.001)
    # the printed text is kept: trailing zeros, the paper's two 100-MeV 7Li rows
    assert [p["as_printed"] for p in li6[:2]] == ["0.10 +- 0.03", "0.050 +- 0.010"]
    assert [p["as_printed"] for p in li7 if p["e0_mev"] == 100.0] == [
        "0.130 +- 0.010", "0.153 +- 0.013"]
    assert all(p["page"] == 865 for p in li6 + li7)
    assert [p["row"] for p in li6] == list(range(1, 8))
    head = open(mod.DATA_FF, encoding="utf-8").read().split("nucleus,table,")[0]
    assert head.startswith("# SPDX-License-Identifier: GPL-3.0-or-later")
    assert "READ FROM THE PAGE IMAGE" in head and _PDF_SHA256 in head
    assert "doi:10.1103/PhysRev.144.859" in head and "NO digit" in head
    assert "PER NUCLEUS" in head and "FOUR-momentum transfer" in head
    # one flipped byte in the body is refused ...
    bad = _flip(mod.DATA_FF, b"6Li,I,4,865,120,1.42,0.010,0.003",
                b"6Li,I,4,865,120,1.42,0.011,0.003", tmp_path)
    with pytest.raises(RuntimeError, match="sha256"):
        mod.load_ff(bad)
    # ... and one in the provenance (every byte is digested)
    bad = _flip(mod.DATA_FF, b"READ FROM THE PAGE IMAGE", b"READ FROM THE PAGE IMAGe",
                tmp_path)
    with pytest.raises(RuntimeError, match="sha256"):
        mod.load_ff(bad)


def test_a0_file_is_tables_vii_viii_xiv(mod, tmp_path):
    rows = mod.load_a0()
    assert len(rows) == 24
    assert (rows[0]["nucleus"], rows[0]["fit"], rows[0]["params"]) == (
        "6Li", "theory:L-S", {"R1": (0.0, None)})
    assert (rows[-1]["nucleus"], rows[-1]["fit"], rows[-1]["a0_fm"],
            rows[-1]["a0_err_fm"], rows[-1]["table"], rows[-1]["page"]) == (
        "7Li", "XIV:mag", 1.70, 0.08, "XIV", 873)
    a = [r for r in rows if r["nucleus"] == "6Li" and r["fit"] == "A"][0]
    assert (a["a0_fm"], a["a0_err_fm"], a["params"]["R1"], a["params"]["chi2"],
            a["params"]["dof"], a["page"]) == (2.16, 0.04, (0.53, 0.03), (1.7, None),
                                              (6.0, None), 869)
    b = [r["a0_fm"] for r in rows if r["nucleus"] == "6Li" and r["fit"] == "B"]
    assert b == [1.78, 2.10]
    rng = [r for r in rows if r["nucleus"] == "6Li" and r["fit"] == "XIV:mag"][0]
    assert rng["a0_fm"] is None and rng["a0_range_fm"] == (1.6, 2.2)
    ad = [r for r in rows if r["fit"] == "alpha+d(1s)"][0]
    assert ad["params"]["q0"] == (0.39, 0.01) and ad["params"]["chi2"] == (2.8, None)
    a7 = [r for r in rows if r["nucleus"] == "7Li" and r["fit"] == "A"][0]
    assert (a7["a0_fm"], a7["a0_err_fm"]) == (1.70, 0.01)
    assert [a7["params"][k] for k in ("A", "B", "C", "R1", "R2", "Omega0/mu0")] == [
        (1.24, 0.05), (0.23, 0.07), (1.18, 0.02), (0.09, 0.04), (0.96, 0.05),
        (2.83, 0.15)]
    ip = [r for r in rows if r["fit"] == "theory:I.P."][0]
    assert [ip["params"][k][0] for k in ("A", "B", "C", "R1", "R2")] == [
        1.29, 0.28, 1.35, 0.116, 1.08]
    head = open(mod.DATA_A0, encoding="utf-8").read().split("nucleus,fit,")[0]
    assert "VERBATIM" in head and _PDF_SHA256 in head
    assert "prints NO magnetization radius" in head
    bad = _flip(mod.DATA_A0, b'"R1=0.53+-0.03; chi2=1.7; dof=6',
                b'"R1=0.54+-0.03; chi2=1.7; dof=6', tmp_path)
    with pytest.raises(RuntimeError, match="sha256"):
        mod.load_a0(bad)


def test_note_parser_refuses_a_malformed_item(mod):
    assert mod.parse_note("R1=-0.02+-0.03; chi2=8.4 | text") == (
        {"R1": (-0.02, 0.03), "chi2": (8.4, None)}, "text")
    assert mod.parse_note(" | only words") == ({}, "only words")
    with pytest.raises(RuntimeError):
        mod.parse_note("R1 0.53 | no equals sign")


# ------------------------------------------------------------ the design block

def test_design_block_is_hashed_and_states_the_tolerance(mod):
    assert mod.design_block_sha256() == mod.DESIGN_BLOCK_SHA256 == (
        "d00b1307a0e397dcad2830513abdcbe21da9fb7eef3d1931e5976a7a0b3a3589")
    block = mod.design_block()
    assert mod.tolerances_from_design_block(block) == (6, 0.01, 1.0, 0.006)
    assert "PASS iff p(chi2_min, ndf = 6) >= 0.01" in block
    assert "WHAT WAS SEEN BEFORE THIS BLOCK WAS WRITTEN" in block
    assert "IS RETIRED" in block
    src = open(os.path.join(_BENCH, _NAME + ".py"), encoding="utf-8").read()
    edited = src.replace("ndf = 6) >= 0.01", "ndf = 6) >= 0.02")
    assert edited != src
    with pytest.raises(RuntimeError, match="design block sha256"):
        mod.design_block_sha256(edited)


def test_tolerance_constants_must_match_the_block(mod, monkeypatch):
    monkeypatch.setattr(mod, "P_MIN", 0.02)
    with pytest.raises(RuntimeError, match="tolerance constants"):
        mod.design_block_sha256()


# ------------------------------------------------------------ the measurement

def test_verdict_is_pinned(rep):
    assert rep["status"] == "pass"
    v = _sub(rep, "li6_shape[free scale]")
    assert v["status"] == "pass" and v["in_verdict"] is True
    assert (v["ndf"], v["n_points"]) == (6, 7)
    assert v["chi2"] == pytest.approx(13.350910518221088, rel=1e-9)
    assert v["p"] == pytest.approx(0.03778982540199176, rel=1e-7)
    assert v["s_hat"] == pytest.approx(1.4188382001858673, rel=1e-9)
    assert v["s_err"] == pytest.approx(0.20576480378243867, rel=1e-9)
    assert v["pulls"] == pytest.approx(
        [-1.2111806717059193, 0.5269786027400948, 0.24267356835826326,
         2.2705591320436564, 1.5576803985964967, 1.991298557777409,
         0.016674762261493666], rel=1e-8)
    # the three points just below the tree's zero carry most of the chi2
    assert sum(p * p for p in v["pulls"][3:6]) == pytest.approx(11.5, abs=0.05)
    model = [p["tree_fmag2"] for p in rep["points"]["li6"]]
    assert model == pytest.approx(
        [0.09608947668121545, 0.031525944231512386, 0.01078437627312804,
         0.0022471361452287945, 0.0002304412188561588, 1.8398381622621635e-05,
         0.0006930495933995089], rel=1e-10)
    assert rep["shape_check_abs"] <= 1e-12


def test_recorded_subrows_are_pinned(rep):
    a = _sub(rep, "li6_absolute[s = 1]")
    assert (a["status"], a["in_verdict"], a["ndf"]) == ("recorded", False, 7)
    assert a["chi2"] == pytest.approx(17.494248761688574, rel=1e-9)
    assert a["p"] == pytest.approx(0.014472632435389845, rel=1e-7)
    e = _sub(rep, "li6_shape[q2 from E0]")
    assert (e["status"], e["in_verdict"]) == ("recorded", False)
    assert e["chi2"] == pytest.approx(13.248113308223, rel=1e-8)
    assert e["s_hat"] == pytest.approx(1.4226764407859973, rel=1e-8)
    q = _sub(rep, "q2_column_check")
    assert (q["status"], q["n_rows"], q["limit_fm2"]) == ("pass", 31, 0.006)
    assert q["max_abs_dq2_fm2"] == pytest.approx(0.005808978315370972, rel=1e-6)
    li7 = _sub(rep, "li7[Table II]")
    assert li7["status"] == "not compared" and li7["n_points"] == 24
    assert li7["reason"].startswith("not implemented in the tree: "
                                    "HoSpin1FF.for_ion(li7()) throws")


def test_transcription_closures_are_pinned(rep):
    c6 = _sub(rep, "closure[6Li Table VII]")
    assert c6["status"] == "pass" and c6["in_verdict"] is False
    assert c6["chi2_paper"] == [1.7, 8.4, 6.0, 3.4, 2.3, 1.8, 1.7, 2.8]
    assert c6["chi2_table"] == pytest.approx(
        [1.436000371530528, 8.430229799355143, 6.0441219296979245, 3.2651587457469473,
         2.0168427562010156, 1.526092999375521, 1.4651782955823165,
         2.4414275667180307], rel=1e-8)
    c7 = _sub(rep, "closure[7Li Table VIII]")
    # RECORDED FAIL: the (C) a0 = 1.66 row's printed parameters give 26.91 > 25.6 + 1
    assert c7["status"] == "fail" and c7["in_verdict"] is False
    assert c7["chi2_paper"] == [23.3, 30.8, 25.6, 23.5, 30.0]
    assert c7["chi2_table"] == pytest.approx(
        [21.814495545018524, 29.374313077461693, 26.90876068115511,
         21.475294971946692, 25.683577437852303], rel=1e-8)
    holds = [t <= p + 1.0 for t, p in zip(c7["chi2_table"], c7["chi2_paper"])]
    assert holds == [True, True, False, True, True]
    r = _sub(rep, "closure_refit[7Li Table VIII, post-hoc]")
    assert r["status"] == "pass" and r["in_verdict"] is False
    assert r["chi2_min"] == pytest.approx(
        [21.3578018272403, 28.816478048174915, 24.531326409315085,
         20.32774246919651, 20.19522569922619], rel=1e-8)
    for cmin, ctab in zip(r["chi2_min"], c7["chi2_table"]):
        assert cmin <= ctab            # a minimum over (A, B, C) never exceeds a point
    d = [x for x in rep["closure_refit_detail"] if x["name"] == "7Li VIII (C) a0=1.66"][0]
    assert d["abc_shift_in_printed_errors"] == pytest.approx([0.499, 0.020, 1.406],
                                                             abs=2e-3)


def test_rfy_model_functions(mod):
    """Eqs. (17), (18) and the alpha+d form at q = 0 and in closed form."""
    assert mod.rfy_eq17_fm2(0.0, 2.16, 0.53, 6) == 1.0
    assert mod.rfy_eq18_fm2(0.0, 1.70, 1.24, 0.23, 1.18, 7) == pytest.approx(1.24,
                                                                           rel=1e-15)
    assert mod.rfy_alpha_d_fm2(0.0, 0.39, 6) == 1.0
    q2, a0, r1 = 1.2, 2.16, 0.53
    x = q2 * a0 * a0 / 4.0
    f = (1.0 - (2.0 / 3.0) * (1.0 - r1) * x) * math.exp(-x * (1.0 - 1.0 / 6.0)) / (
        1.0 + q2 / 7.1)
    assert mod.rfy_eq17_fm2(q2, a0, r1, 6) == pytest.approx(f * f, rel=1e-13)
    g = math.exp(-q2 * 1.8 ** 2 / 6.0 - q2 / (4.0 * 0.39 ** 2) * (2.0 / 3.0))
    assert mod.rfy_alpha_d_fm2(q2, 0.39, 6) == pytest.approx(g * g, rel=1e-13)


def test_kinematics_at_180_degrees(mod):
    """Record (r3): q^2 = 4 E0 E'/(hbar c)^2, E' = E0/(1 + 2 E0/M_A)."""
    m6 = lg.li6().mass()
    q2 = mod.q2_at_pi(85.0, m6, lg.HBARC_GEV_FM)
    e1 = 0.085 / (1.0 + 0.17 / m6)
    assert q2 == pytest.approx(4.0 * 0.085 * e1 / lg.HBARC_GEV_FM ** 2, rel=1e-15)
    assert round(q2, 2) == 0.72
    # the 100-MeV rows: 0.99 for 6Li and 1.00 for 7Li, the recoil telling them apart
    assert round(mod.q2_at_pi(100.0, m6, lg.HBARC_GEV_FM), 2) == 0.99
    assert round(mod.q2_at_pi(100.0, lg.li7().mass(), lg.HBARC_GEV_FM), 2) == 1.00


def test_tree_side_and_rule3(rep):
    t = rep["tree"]
    assert (t["fm_qz_fm"], t["fm_b_fm"]) == (1.30, 1.85)    # rc.hpp, unfitted
    assert t["r_mag_li6_closed"] == pytest.approx(2.94687, abs=5e-5)
    assert t["r_mag_li6"] == pytest.approx(t["r_mag_li6_closed"], rel=1e-5)
    assert t["li7"].startswith("not implemented in the tree")
    assert t["li7_implemented"] is False
    assert rep["config"]["tree_r_mag_li6_fm_uncompared"] == t["r_mag_li6"]
    assert rep["config"]["radius_rule"].startswith("retired 2026-10-26")
    # rule 3: F_m(t)/F_m(0) does not see tail_tensor_scale or M_A
    hc2 = lg.HBARC_GEV_FM ** 2
    ref = lg.HoSpin1FF.for_ion(lg.li6())
    o = lg.HoSpin1FFOptions()
    o.tail_tensor_scale = 2.0
    o.m_a_gev = 7.0
    alt = lg.HoSpin1FF.for_ion(lg.li6(), o)
    assert alt.fm(0.0) == pytest.approx(2.0 * 7.0 / ref.options.m_a_gev * ref.fm(0.0),
                                        rel=1e-14)
    for q2 in (0.72, 1.42, 1.92):
        assert alt.fm(q2 * hc2) / alt.fm(0.0) == pytest.approx(
            ref.fm(q2 * hc2) / ref.fm(0.0), rel=1e-14)


def test_verdict_follows_the_numbers(mod, rep):
    li6, _ = mod.load_ff()
    d = [p["fm2"] for p in li6]
    sig = [p["fm2_err"] for p in li6]
    m = [p["tree_fmag2"] for p in rep["points"]["li6"]]
    v = _sub(rep, "li6_shape[free scale]")
    # s_hat minimises chi2(s): a nudge either way raises it
    for ds in (-1e-4, 1e-4):
        assert mod.chi2_fixed(d, sig, m, v["s_hat"] + ds) > v["chi2"]
    assert v["chi2"] == pytest.approx(sum(p * p for p in v["pulls"]), rel=1e-14)
    assert v["p"] == pytest.approx(mod.chi2_sf(v["chi2"], 6), rel=1e-14)
    assert v["p"] >= mod.P_MIN
    # the scipy-free fallback agrees with scipy where both exist
    pytest.importorskip("scipy")
    for chi2, ndf in ((13.350910518221088, 6), (17.494248761688574, 7), (0.5, 3),
                      (80.0, 6)):
        assert mod.chi2_sf_fallback(chi2, ndf) == pytest.approx(
            mod.chi2_sf(chi2, ndf), rel=1e-9)


def test_a0_record_is_the_1p_orbital_rms(rep):
    """Record (r6): sqrt(5/2) a0 reproduces, at the printed digit, the four
    lithium 'magnetization radii' docs/benchmarking/06_critic.md sec. 4.6
    transcribed from ADNDT 14 Table V -- they are RFY's 1p-orbital rms."""
    rec = {(r["nucleus"], r["fit"], r["a0_as_printed"]): r for r in rep["a0_record"]}
    a6 = rec[("6Li", "A", "2.16")]
    assert (round(a6["r_1p_orbital_rms_fm"], 2), round(a6["r_1p_orbital_rms_err_fm"], 2)) \
        == (3.42, 0.06)
    assert round(rec[("6Li", "B", "1.78")]["r_1p_orbital_rms_fm"], 2) == 2.81
    x7 = rec[("7Li", "XIV:mag", "1.70")]
    assert (round(x7["r_1p_orbital_rms_fm"], 2), round(x7["r_1p_orbital_rms_err_fm"], 2)) \
        == (2.69, 0.13)
    assert round(rec[("7Li", "B", "1.72")]["r_1p_orbital_rms_fm"], 2) == 2.72
    assert rec[("6Li", "XIV:mag", "1.6-2.2")]["r_1p_orbital_rms_range_fm"] == \
        pytest.approx([math.sqrt(2.5) * 1.6, math.sqrt(2.5) * 2.2], rel=1e-15)
    # never a sub-row: no radius is compared
    assert not [s for s in rep["subrows"] if "radius" in s["name"] or "r_mag" in s["name"]]


def test_report_line(mod, capsys):
    r = mod.run(verbose=True)
    out = capsys.readouterr().out.splitlines()
    rows = [l for l in out if l.startswith("REPORT |")]
    subs = [l for l in out if l.startswith("SUBROW |")]
    assert len(rows) == 1 and len(subs) == 8
    assert rows[0].startswith(
        "REPORT | t3_li_magnetization_rfy | HoSpin1FF F_mag^2 (q_z = 1.30 fm^-1, "
        "b = 1.85 fm; unfitted starting values) at RFY Table I's 7 points: chi2 = "
        "13.35 / 6 with one free scale s = 1.419 +- 0.206 (p = 0.0378);")
    assert rows[0].endswith("| PASS")
    assert subs[0].startswith("SUBROW | t3_li_magnetization_rfy[li6_shape[free scale]] "
                              "| chi2 = 13.3509 / 6")
    assert "np.float64" not in "\n".join(out)
    assert r["runtime_s"] < 30.0 and mod.EXPENSIVE is False


def test_exit_codes(tmp_path):
    ok = _run_harness(os.path.join(_BENCH, _NAME + ".py"))
    assert ok.returncode == 0, ok.stderr
    assert "REPORT | t3_li_magnetization_rfy |" in ok.stdout
    cases = [
        dict(data_name=_FF, old=b"7Li,II,24,865,200,3.87,0.016,0.001",
             new=b"7Li,II,24,865,200,3.87,0.017,0.001"),
        dict(data_name=_A0, old=b'"q0=0.39+-0.01;', new=b'"q0=0.38+-0.01;'),
        dict(harness_edit=("about 13 for 6 dof", "about 14 for 6 dof")),
    ]
    for i, case in enumerate(cases):
        script = _copy_bench(tmp_path / str(i), **case)
        bad = _run_harness(script)
        assert bad.returncode == 2, (case, bad.stdout[-500:])
        assert "sha256" in bad.stderr
        assert "REPORT |" not in bad.stdout


def test_plot_writes_a_png(mod, rep, tmp_path):
    pytest.importorskip("matplotlib")
    path = mod.plot(rep, out_dir=str(tmp_path))
    assert path == str(tmp_path / (_NAME + ".png"))
    raw = open(path, "rb").read()
    assert len(raw) > 10000 and raw[:8] == b"\x89PNG\r\n\x1a\n"
