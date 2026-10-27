# SPDX-License-Identifier: GPL-3.0-or-later
"""validation/benchmarks/t2_boeglin_nk.py exercised as tests --
docs/benchmarking/02_data_nucleon_deuteron.md sec. 5.4, Boeglin et al., PRL
107 (2011) 262501, Fig. 2 read back from the PDF's vector drawing (run
2026-10-26 phase C5; < 1 s per harness run, so not behind LIPOLGEN_BENCH).

The harness holds the three deuteron waves' n(k = p_m) to the reduced cross
section at theta_nq = 35 and 45 deg (22 points, p_m <= 0.54 GeV/c) with one
free scale, the bars and the paper's own FSI / off-shell allowance (20 %),
p >= 0.01, declared in a hashed design block; it records the rest (the
absolute closure against the paper's CD-Bonn PWIA, theory against theory,
no allowance, per panel, the MS FSI factor, and a post-hoc smooth reading of
the TaggedModel table).  These tests pin what the harness MEASURES and that
its verdict follows from the numbers; nothing here moves or asserts a new
value for any tree default.  Nothing needs the network.  The one test that
re-reads the PDF (`--redigitize`) skips unless PyMuPDF and the PDF are
present ($LIPOLGEN_BOEGLIN_PDF, else the sibling corpus's refs/).
"""
import hashlib
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
_NAME = "t2_boeglin_nk"
_DATA = "boeglin2011_deuteron_pmiss.csv"


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


def _copy_bench(tmp_path, old=None, new=None, harness_edit=None):
    """The harness and its vendored file in a fresh directory, one of them
    altered by one replacement."""
    d = tmp_path / "bench"
    (d / "data").mkdir(parents=True)
    src = open(os.path.join(_BENCH, _NAME + ".py"), encoding="utf-8").read()
    if harness_edit is not None:
        assert src.count(harness_edit[0]) == 1
        src = src.replace(*harness_edit)
    (d / (_NAME + ".py")).write_text(src, encoding="utf-8")
    raw = open(os.path.join(_BENCH, "data", _DATA), "rb").read()
    if old is not None:
        assert raw.count(old) == 1, old
        raw = raw.replace(old, new)
    (d / "data" / _DATA).write_bytes(raw)
    return str(d / (_NAME + ".py"))


@pytest.fixture(scope="module")
def mod():
    return _load()


@pytest.fixture(scope="module")
def rep(mod):
    return mod.run(verbose=False)


# ------------------------------------------------------------ the vendored file

def test_vendored_file_guards(mod, tmp_path):
    rows = mod.load_reference()
    assert len(rows) == 35
    assert {p: sum(1 for r in rows if r["panel"] == p) for p in "abc"} == {
        "a": 10, "b": 12, "c": 13}
    first, last = rows[0], rows[-1]
    assert (first["panel"], first["theta_nq"], first["point"], first["pm"], first["d"],
            first["lo"], first["hi"], first["cdbonn_pwia"]) == (
        "a", 35, 1, 0.17833, 0.0148438, 0.0141679, 0.0155189, 0.0217689)
    assert (last["panel"], last["theta_nq"], last["point"], last["pm"], last["d"],
            last["paris_fsi"]) == ("c", 75, 13, 0.53785, 0.000287569, 0.000226531)
    raw = open(mod.DATA, encoding="utf-8").read()
    assert raw.startswith("# SPDX-License-Identifier: GPL-3.0-or-later\n")
    for needle in ("b9abd4132e4736ac195d608d43354d3774c2ca233a587495adb1e31395c70f91",
                   "Phys. Rev. Lett.", "(114.973, 63.567, 251.737, 258.945)",
                   "--redigitize", "tools/digitize_figure.py",
                   "PUBLISHING CONVENTION (per deuteron)"):
        assert needle in raw
    assert hashlib.sha256(raw.encode("utf-8")).hexdigest() == mod.DATA_SHA256
    bad = _flip(mod.DATA, b"a,35,1,0.17833,0.0148438,", b"a,35,1,0.17833,0.0148439,", tmp_path)
    with pytest.raises(RuntimeError, match="sha256"):
        mod.load_reference(bad)
    bad = _flip(mod.DATA, b"# Licence:   facts", b"# Licence:   Facts", tmp_path)
    with pytest.raises(RuntimeError, match="sha256"):
        mod.load_reference(bad)


def test_every_number_is_rederived_from_its_coordinates(mod, tmp_path, monkeypatch):
    bad = _flip(mod.DATA, b"a,35,1,0.17833,0.0148438,", b"a,35,1,0.17833,0.0158438,", tmp_path)
    monkeypatch.setattr(mod, "DATA_SHA256", hashlib.sha256(open(bad, "rb").read()).hexdigest())
    with pytest.raises(RuntimeError, match="calibrated"):
        mod.load_reference(bad)


def test_calibration_and_bars(mod):
    for p, (x0, y0, x1, y1) in mod.FRAMES.items():
        assert mod.to_pm(x0, p) == 0.0 and mod.to_pm(x1, p) == pytest.approx(0.7, abs=1e-12)
        assert mod.to_log10(y1, p) == -5.0 and mod.to_log10(y0, p) == pytest.approx(1.0, abs=1e-12)
    # the six major x ticks of panel (a) read back 0.1 .. 0.6
    for i, x in enumerate((134.511, 154.048, 173.586, 193.124, 212.662, 232.200), 1):
        assert mod.to_pm(x, "a") == pytest.approx(0.1 * i, abs=2e-5)
    # the bars are symmetric in sigma_red (linear)
    for r in mod.load_reference():
        assert abs((r["hi"] - r["d"]) - (r["d"] - r["lo"])) <= 2e-3 * (r["hi"] - r["lo"]) + 1e-12
    assert mod.DPM_GEV == pytest.approx(0.01 * 0.7 / 136.764, rel=1e-6)
    assert mod.DLOG10 == pytest.approx(0.01 * 6.0 / 195.378, rel=1e-6)


# ------------------------------------------------------------ the design block

def test_design_block_is_hashed_and_states_the_rules(mod):
    block = mod.design_block()
    assert hashlib.sha256(block.encode("utf-8")).hexdigest() == mod.DESIGN_BLOCK_SHA256
    assert mod.design_block_sha256() == mod.DESIGN_BLOCK_SHA256
    assert mod.rules_from_design_block(block) == dict(
        ALLOWANCE=0.20, NDF_VERDICT=21, P_MIN=0.01, DIG_LOG10=0.000307)
    assert "WHAT WAS SEEN BEFORE THIS BLOCK WAS WRITTEN" in block
    src = open(os.path.join(_BENCH, _NAME + ".py"), encoding="utf-8").read()
    with pytest.raises(RuntimeError, match="design block"):
        mod.design_block_sha256(src.replace("A = 0.20 -- the upper", "A = 0.30 -- the upper", 1))


def test_rule_constants_must_match_the_block(mod, monkeypatch):
    monkeypatch.setattr(mod, "ALLOWANCE", 0.3)
    with pytest.raises(RuntimeError, match="ALLOWANCE"):
        mod.design_block_sha256()


# ------------------------------------------------------------ the verdict

def test_verdict_is_pinned(rep):
    v = rep["verdict"]
    assert rep["status"] == "pass"
    assert (v["hulthen"]["chi2"], v["hulthen"]["ndf"]) == (pytest.approx(21.98372, abs=1e-4), 21)
    assert v["hulthen"]["p"] == pytest.approx(0.400449, abs=1e-5)
    assert (v["hulthen"]["s"], v["hulthen"]["s_err"]) == pytest.approx((0.481148, 0.024203), abs=2e-6)
    assert (v["av18"]["chi2"], v["av18"]["p"]) == pytest.approx((11.53753, 0.951270), abs=1e-5)
    assert (v["cdbonn"]["chi2"], v["cdbonn"]["p"]) == pytest.approx((8.88738, 0.990073), abs=1e-5)
    assert [v[w]["status"] for w in ("hulthen", "av18", "cdbonn")] == ["pass"] * 3
    assert [s["in_verdict"] for s in rep["subrows"][:3]] == [True, False, False]


def test_verdict_follows_the_numbers(mod, rep):
    v = rep["verdict"]
    for w in mod.WAVES:
        assert v[w]["status"] == ("pass" if v[w]["p"] >= 0.01 else "fail")
        assert v[w]["p"] == pytest.approx(mod.chi2_sf(v[w]["chi2"], 21), rel=1e-12)
    # the allowance decides it: without it every wave fails (record (r3))
    na = rep["no_allowance"]
    assert na["hulthen"]["chi2"] == pytest.approx(185.9066, abs=1e-3)
    assert (na["av18"]["chi2"], na["cdbonn"]["chi2"]) == pytest.approx((66.4996, 60.3476), abs=1e-3)
    assert all(na[w]["p"] < 0.01 for w in mod.WAVES)
    # the fallback p equals scipy's where both exist
    assert mod._gammq(10.5, 0.5 * 21.98372) == pytest.approx(0.400449, abs=2e-6)
    # a model equal to the data passes with chi2 = 0 and s = 1
    win = [r for r in mod.load_reference() if r["panel"] in mod.WINDOW_PANELS]
    f = mod.scaled_chi2(win, lambda r: r["d"])
    assert f["chi2"] == pytest.approx(0.0, abs=1e-20) and f["s"] == pytest.approx(1.0, rel=1e-12)


def test_recorded_items_are_pinned(rep):
    c = rep["closure"]
    assert (c["mean"], c["min"], c["max"]) == pytest.approx((0.901930, 0.875820, 0.934274), abs=2e-6)
    pc = rep["paris_over_cdbonn"]
    assert (pc["mean"], pc["min"], pc["max"]) == pytest.approx((1.160545, 0.939071, 1.714466), abs=2e-6)
    t = rep["theory"]
    assert (t["cdbonn"]["min"], t["cdbonn"]["max"]) == pytest.approx((0.971230, 1.036052), abs=2e-6)
    assert (t["hulthen"]["min"], t["hulthen"]["max"]) == pytest.approx((0.823027, 1.299910), abs=2e-6)
    pp = rep["per_panel"]
    assert (pp["a"]["hulthen"]["chi2"], pp["b"]["hulthen"]["chi2"]) == pytest.approx(
        (4.243175, 17.471651), abs=1e-5)
    assert pp["b"]["hulthen"]["p"] == pytest.approx(0.094678, abs=2e-6)
    f = rep["ms_fsi"]
    assert (f["min"], f["max"], f["mean"]) == pytest.approx((0.746105, 1.108599, 0.855088), abs=2e-6)
    r = rep["ratio_route"]
    assert (r["a"]["mean"], r["b"]["mean"]) == pytest.approx((0.734752, 0.836996), abs=2e-6)
    fc = rep["fsi_corrected"]
    assert (fc["hulthen"]["chi2"], fc["av18"]["chi2"], fc["cdbonn"]["chi2"]) == pytest.approx(
        (95.63767, 76.69770, 24.39897), abs=1e-4)
    assert fc["cdbonn"]["p"] == pytest.approx(0.274131, abs=2e-6)
    sm = rep["smooth_table"]
    assert (sm["hulthen"]["chi2"], sm["av18"]["chi2"]) == pytest.approx((23.42066, 14.60900), abs=1e-4)
    cb = rep["cell_bias"]
    assert cb["hulthen"]["at_lowest_two"] == pytest.approx([1.160480, 1.080289], abs=2e-6)
    assert cb["av18"]["max"] == pytest.approx(1.180757, abs=2e-6)


def test_rule3_and_the_tree_side(mod):
    """Per deuteron: the CD-Bonn n is normalised to INT n d^3p = 1 (fm^3);
    the TaggedModel waves' grid quadrature is 1 to 1e-4; the lookup is the
    cell AT OR BELOW k, so a query exactly at a node returns the cell below
    and the node value is the density there."""
    tree = mod.Tree()
    hc = lg.HBARC_GEV_FM
    ks = [0.0005 * i for i in range(1, 6001)]               # GeV, to 3 GeV
    integ = sum(4 * math.pi * (k / hc) ** 2 * tree.n("cdbonn", k) for k in ks) * 0.0005 / hc
    assert integ == pytest.approx(1.0, abs=2e-4)
    t = tree.models["hulthen"]
    assert all(abs(t.norm(m) - 1.0) < 1e-4 for m in (1.0, 0.0, -1.0))
    nodes = tree.nodes()
    j = 12
    assert tree.n("hulthen", nodes[j]) == tree.n_node("hulthen", j - 1)
    assert tree.n("hulthen", 0.5 * (nodes[j] + nodes[j + 1])) == tree.n_node("hulthen", j)
    r0, r2 = list(t.radial_table(0)), list(t.radial_table(2))
    assert tree.n_node("hulthen", j) == pytest.approx(
        (r0[j] ** 2 + r2[j] ** 2) / (4 * math.pi) * hc ** 3, rel=1e-12)


# ------------------------------------------------------------ the row

def test_report_line(mod, capsys):
    r = mod.run(verbose=True)
    out = capsys.readouterr().out.splitlines()
    rows = [l for l in out if l.startswith("REPORT |")]
    assert len(rows) == 1
    assert rows[0].startswith(
        "REPORT | t2_boeglin_nk | n(k = p_m) x s, one free scale: hulthen (default) "
        "chi2 = 21.98/21 (p = 0.4, s = 0.481); av18 11.54 (p = 0.951); cdbonn 8.89 "
        "(p = 0.99) | Boeglin et al. PRL 107 (2011) 262501 Fig. 2")
    assert rows[0].endswith("| PASS")
    subs = [l for l in out if l.startswith("SUBROW |")]
    assert len(subs) == len(r["subrows"]) == 29
    assert "np.float64" not in "\n".join(out)
    assert r["runtime_s"] < 30.0


def test_exit_codes(tmp_path):
    ok = _run_harness(os.path.join(_BENCH, _NAME + ".py"))
    assert ok.returncode == 0, ok.stderr
    assert "REPORT | %s |" % _NAME in ok.stdout
    bad = _run_harness(_copy_bench(tmp_path / "a", b"a,35,1,0.17833,0.0148438,",
                                   b"a,35,1,0.17833,0.0148439,"))
    assert bad.returncode == 2 and "sha256" in bad.stderr
    assert "REPORT |" not in bad.stdout
    bad = _run_harness(_copy_bench(tmp_path / "b", harness_edit=(
        "TOLERANCE: PASS iff p(chi2_min, ndf = 21) >= 0.01.",
        "TOLERANCE: PASS iff p(chi2_min, ndf = 21) >= 0.001.")))
    assert bad.returncode == 2 and "design block" in bad.stderr


def test_redigitize_reproduces_the_body(mod):
    """DEV TIME: re-read the PDF's drawing operators and compare with the
    vendored body, line for line.  Needs PyMuPDF and the arXiv PDF (not in
    this repository): $LIPOLGEN_BOEGLIN_PDF, else the sibling corpus."""
    pytest.importorskip("pymupdf", reason="PyMuPDF (dev-time re-digitization only)")
    pdf = os.environ.get("LIPOLGEN_BOEGLIN_PDF") or os.path.join(
        os.path.dirname(_ROOT), "PolarizedLithiumSim", "refs", "1106.0275.pdf")
    if not os.path.isfile(pdf):
        pytest.skip("the Boeglin PDF (arXiv:1106.0275v3, sha256 %s...) is not at %s; "
                    "set LIPOLGEN_BOEGLIN_PDF" % (mod.BG_PDF_SHA256[:12], pdf))
    res = _run_harness(os.path.join(_BENCH, _NAME + ".py"), ["--redigitize", pdf])
    assert res.returncode == 0, res.stdout + res.stderr
    assert "identical (36 lines)" in res.stdout


def test_plot_writes_png(mod, rep, tmp_path):
    """plot() draws the row's figure -- a view of the report, never the
    record -- as a non-empty PNG where it is told to (here a tmp dir)."""
    pytest.importorskip("matplotlib", reason="matplotlib (the figure only)")
    path = mod.plot(rep, out_dir=str(tmp_path))
    assert path == os.path.join(str(tmp_path), _NAME + ".png")
    with open(path, "rb") as f:
        assert f.read(8) == b"\x89PNG\r\n\x1a\n"
    assert os.path.getsize(path) > 0
