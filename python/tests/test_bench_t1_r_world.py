# SPDX-License-Identifier: GPL-3.0-or-later
"""Run 2026-10-26 phase C1 -- the two harnesses it added, as tests:

  t1_r_world.py     R = sigma_L/sigma_T: the shipped r_sigma_lt and r1998
                    against NMC's R (PRIMARY), H1's F_L (EXTRAPOLATION),
                    Whitlow's global SLAC R and E140's R (CLOSURE for r1998),
                    and the (1 + R) bias table at the letter's sweet spots
  t1_hermes_f2d.py  HERMES F2d / F2p against ToyF2 / MstwSF / CT18NLO, and
                    the HERMES b1 inversion against HERMES's own F2d

Their references are vendored from HEPData (CC0) under
validation/benchmarks/data/; nothing here needs the network.  Both harnesses
cost < 3 s, so nothing sits behind LIPOLGEN_BENCH.

A harness that DISAGREES with the tree is not a failing test: the tests pin
what each harness MEASURES (a recorded FAIL is pinned as a FAIL), check that
each verdict follows from those numbers by the rule its header declares, and
guard every vendored file.  No test here moves, or asserts a new value for,
any tree default.  A cell that needs the PYTHIA tier or an LHAPDF set skips
naming exactly what is missing, or checks that it is REPORTED as blocked.
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


def _load(name):
    path = os.path.join(_BENCH, name + ".py")
    spec = importlib.util.spec_from_file_location("bench_" + name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _sf(chi2, n):
    """p(chi2, n) computed HERE -- scipy, or without it (CI) the fallback,
    itself checked against scipy where present."""
    try:
        from scipy.stats import chi2 as c
    except ImportError:
        return _load("t1_nmc_f2d").chi2_sf_fallback(chi2, n)
    return float(c.sf(chi2, n))


def _env():
    env = dict(os.environ)
    here = os.path.dirname(os.path.dirname(os.path.abspath(lg.__file__)))
    env["PYTHONPATH"] = os.pathsep.join(
        [here] + [p for p in env.get("PYTHONPATH", "").split(os.pathsep) if p])
    return env


def _run_script(script, *args):
    return subprocess.run([sys.executable, script] + list(args),
                          capture_output=True, text=True, env=_env(),
                          timeout=300)


def _have_mstw():
    import lipolgen._lipolgen as _l
    if not getattr(_l, "HAVE_PYTHIA8", False):
        return False, "environment: lipolgen built without the PYTHIA 8 tier (MstwSF)"
    d = _l.pythia8_pdfdata_dir()
    if not (d and os.path.isfile(os.path.join(d, "mstw2008lo.00.dat"))):
        return False, "environment: PYTHIA 8 pdfdata/mstw2008lo.00.dat not found (MstwSF)"
    return True, ""


def _lhapdf_missing(setname):
    """True iff the LHAPDF tier is off or reports `setname` not installed."""
    if not getattr(lg, "HAVE_LHAPDF", False):
        return True
    try:
        lg.LhapdfSF(setname)
    except RuntimeError as exc:
        if "Info file not found for PDF set" not in str(exc):
            raise
        return True
    return False


# ------------------------------------------------------------------ fixtures

@pytest.fixture(scope="module")
def rmod():
    return _load("t1_r_world")


@pytest.fixture(scope="module")
def rrep(rmod):
    return rmod.run(verbose=False)


@pytest.fixture(scope="module")
def hmod():
    return _load("t1_hermes_f2d")


@pytest.fixture(scope="module")
def hrep(hmod):
    return hmod.run(verbose=False)


def _flip_first_value(raw, marker):
    """The bytes with ONE byte changed: the last digit of the first
    `"value": ` after `marker` (still valid JSON, a different number)."""
    start = raw.index(b'"value": ', raw.index(marker))
    end = raw.index(b",", start) - 1
    if not raw[end:end + 1].isdigit():          # value is the last key
        end = raw.index(b"\n", start) - 1
    assert raw[end:end + 1].isdigit()
    new = b"7" if raw[end:end + 1] != b"7" else b"8"
    return raw[:end] + new + raw[end + 1:]


# ----------------------------------------------------- the vendored files

# (sub-row key, file, list path in body, n, first (x, q2, value), last (...))
_R_FILES = [
    ("nmc", "hepdata_ins424154_nmc_r.json", "points", 10,
     (0.0045, 1.38, 0.537), (0.11, 20.6, 0.043)),
    ("h1", "hepdata_ins786161_h1_fl.json", "averaged", 8,
     (0.00028, 12.0, 0.22), (0.0036, 90.0, 0.48)),
    ("h1", "hepdata_ins786161_h1_fl.json", "per_bin", 25,
     (0.000236, 12.0, 0.34), (0.003617, 90.0, 0.48)),
    ("r1990", "hepdata_ins296980_whitlow_r.json", "points", 89,
     (0.1, 0.63, 0.256), (0.86, 20.0, 0.082)),
    ("e140", "hepdata_ins360765_e140_r.json", "points", 11,
     (0.2, 1.0, 0.376), (0.5, 10.0, 0.046)),
]


@pytest.mark.parametrize("key,fname,lst,n,first,last", _R_FILES)
def test_r_files_are_the_records_and_refuse_a_flipped_byte(
        rmod, key, fname, lst, n, first, last, tmp_path):
    path = os.path.join(_DATA, fname)
    assert rmod.FILES[key][0] == fname
    raw = open(path, "rb").read()
    assert hashlib.sha256(raw).hexdigest() == rmod.FILES[key][1]
    doc = rmod.read_vendored(key)
    pts = doc["body"][lst]
    assert len(pts) == n
    assert (pts[0]["x"], pts[0]["q2"], pts[0]["value"]) == first
    assert (pts[-1]["x"], pts[-1]["q2"], pts[-1]["value"]) == last
    prov = doc["_provenance"]
    assert prov["licence"].startswith("CC0-1.0")
    assert prov["record"]["inspire"] == "ins" + fname.split("_")[1][3:]
    for t in prov["tables_used"]:
        assert len(t["yaml_sha256"]) == 64 and t["table_license"].startswith("CC0")
    # one flipped byte in the body is refused ...
    bad = tmp_path / fname
    bad.write_bytes(_flip_first_value(raw, b'"%s": [' % lst.encode()))
    with pytest.raises(RuntimeError, match="sha256"):
        rmod.read_vendored(key, str(bad))
    # ... and so is one in the provenance (every byte is digested)
    i = raw.index(b"CC0-1.0")
    bad.write_bytes(raw[:i] + b"CC1" + raw[i + 3:])
    with pytest.raises(RuntimeError, match="sha256"):
        rmod.read_vendored(key, str(bad))


def test_r_files_carry_their_r1998_fit_status_and_the_paper_facts(rmod):
    """The CLOSURE labels rest on the papers: NMC, Whitlow and E140 are R1998
    fit inputs (Abe et al. [12], [3], via [3]); H1 is not."""
    yes = {k: rmod.read_vendored(k)["_provenance"]["in_r1998_fit"] for k in rmod.FILES}
    assert yes["nmc"].startswith("YES") and "[12]" in yes["nmc"]
    assert yes["r1990"].startswith("YES") and "[3]" in yes["r1990"]
    assert yes["e140"].startswith("YES")
    assert yes["h1"].startswith("NO")
    e140 = rmod.read_vendored("e140")
    assert e140["body"]["paper_normalisation_error_r"] == 0.03
    q = e140["_provenance"]["paper_quotes"]["ra_minus_rd_average"]
    assert (q["value"], q["stat"], q["sys"]) == (0.001, 0.018, 0.016)
    w = rmod.read_vendored("r1990")["body"]["rd_minus_rp"]
    assert (w["value"], w["stat"], w["sys"]) == (-0.001, 0.009, 0.009)
    nmc = rmod.read_vendored("nmc")["_provenance"]
    assert len(nmc["f2d_tables_checked"]) == 16
    h1 = rmod.read_vendored("h1")["body"]["averaged"]
    assert [p["total_as_printed"] for p in h1] == [0.11, 0.11, 0.10, 0.10, 0.13,
                                                   0.18, 0.27, 0.39]


def test_hermes_file_is_the_record_and_refuses_a_flipped_byte(hmod, tmp_path):
    raw = open(hmod.DATA, "rb").read()
    assert hashlib.sha256(raw).hexdigest() == hmod.FILE_SHA256
    doc = hmod.read_vendored()
    for target, norm, first, last in (
            ("deuteron", 7.5, (0.0069, 0.291, 0.164), (0.2759, 7.501, 0.249)),
            ("proton", 7.6, (0.0069, 0.291, 0.162), (0.2759, 7.501, 0.296))):
        t = doc["body"][target]
        assert len(t["points"]) == 81 and t["norm_pct"] == norm
        p0, p1 = t["points"][0], t["points"][-1]
        assert (p0["x"], p0["q2"], p0["value"]) == first
        assert (p1["x"], p1["q2"], p1["value"]) == last
    prov = doc["_provenance"]
    assert prov["licence"].startswith("CC0-1.0")
    assert prov["record_defect"].startswith("RECORD DEFECT")
    assert "per nucleon" in prov["publishing_convention"]
    bad = tmp_path / "hermes.json"
    bad.write_bytes(_flip_first_value(raw, b'"deuteron": {'))
    with pytest.raises(RuntimeError, match="sha256"):
        hmod.load_table("deuteron", str(bad))
    i = raw.index(b"RECORD DEFECT")
    bad.write_bytes(raw[:i] + b"RECORD DEFECt" + raw[i + 13:])
    with pytest.raises(RuntimeError, match="sha256"):
        hmod.read_vendored(str(bad))


def test_the_body_digest_is_checked_on_its_own(rmod, hmod, tmp_path, monkeypatch):
    """A file whose every-byte digest is updated but whose body moved is still
    refused by the canonical-body digest."""
    raw = open(os.path.join(_DATA, rmod.FILES["nmc"][0]), "rb").read()
    bad = _flip_first_value(raw, b'"points": [')
    p = tmp_path / "nmc.json"
    p.write_bytes(bad)
    monkeypatch.setitem(rmod.FILES, "nmc", (rmod.FILES["nmc"][0],
                                            hashlib.sha256(bad).hexdigest(),
                                            rmod.FILES["nmc"][2]))
    with pytest.raises(RuntimeError, match="body sha256"):
        rmod.read_vendored("nmc", str(p))
    raw = open(hmod.DATA, "rb").read()
    bad = _flip_first_value(raw, b'"proton": {')
    p = tmp_path / "hermes.json"
    p.write_bytes(bad)
    monkeypatch.setattr(hmod, "FILE_SHA256", hashlib.sha256(bad).hexdigest())
    with pytest.raises(RuntimeError, match="body sha256"):
        hmod.read_vendored(str(p))


def test_b1_inversion_inputs_are_guarded(hmod, tmp_path):
    """The two older files the inversion reads are sha256-checked with the
    digests t2_hermes_b1_table2.py records for them."""
    t2 = _load("t2_hermes_b1_table2")
    assert hmod.B1_SHA256 == t2.DATA_SHA256 and hmod.R1990_SHA256 == t2.R1990_SHA256
    pts = hmod.load_table("deuteron")
    for src, kw in ((hmod.B1_DATA, "b1_path"), (hmod.R1990_DATA, "r1990_path")):
        raw = open(src, "rb").read()
        bad = tmp_path / os.path.basename(src)
        bad.write_bytes(raw.replace(b"0.", b"0,", 1) if kw == "b1_path"
                        else raw.replace(b'"b2": 0.5747', b'"b2": 0.5748'))
        assert bad.read_bytes() != raw
        with pytest.raises(RuntimeError, match="sha256"):
            hmod.b1_inversion(pts, **{kw: str(bad)})


# ------------------------------------------------------ t1_r_world: design

def test_r1998_support_is_the_shipped_clip(rmod, monkeypatch):
    sup = rmod.check_r1998_support()
    assert sup == dict(x=(0.005, 0.86), q2=(0.5, 130.0))
    monkeypatch.setattr(rmod, "R1998_X_MIN", 0.004)
    with pytest.raises(RuntimeError, match="no longer match"):
        rmod.check_r1998_support()


def test_dis_cut_and_support_counts(rrep):
    d = rrep["detail"]
    assert [(d[k]["n_points"], d[k]["n_dis"]) for k in ("nmc", "h1", "r1990", "e140")] \
        == [(10, 10), (8, 8), (89, 68), (11, 9)]
    # outside R1998's support: NMC's lowest bin and every H1 point
    assert [d[k]["n_extrapolated"] for k in ("nmc", "h1", "r1990", "e140")] == [1, 8, 0, 0]
    assert [p["x"] for p in d["nmc"]["points"] if not p["in_r1998_support"]] == [0.0045]
    for p in d["e140"]["points"]:
        assert p["dis"] == (p["q2"] >= 1.0 and p["w2"] >= 4.0)


def test_r_theory_is_the_tree_backends_at_each_point(rmod, rrep):
    for key in ("nmc", "r1990", "e140"):
        m = rrep["detail"][key]
        for p, a, b, c in zip(m["points"], m["backends"]["r_sigma_lt"]["theory_all"],
                              m["backends"]["r1998"]["theory_all"],
                              m["backends"]["r1998_unclipped"]["theory_all"]):
            assert a == lg.r_sigma_lt(p["x"], p["q2"]) == 0.18 / (1.0 + p["q2"] / 50.0)
            assert b == lg.r1998(p["x"], p["q2"])
            assert c == lg.r1998(p["x"], p["q2"], clip=False)
            if p["in_r1998_support"]:
                assert b == c


def test_h1_fl_is_the_trees_flp_with_the_backends_r(rmod, rrep):
    """F_L compared directly: the tree side is UnpolSF.flp = F2p R/(1+R);
    the null hook IS r_sigma_lt (the shipped F_L), the others are installed."""
    pts = [p for p in rrep["detail"]["h1"]["points"] if p["dis"]]
    fl, f2 = rmod.tree_fl(lambda: lg.ToyF2(), "r_sigma_lt", pts)
    for p, v, f in zip(pts, fl, f2):
        assert v == lg.ToyF2().flp(p["x"], p["q2"])
        r = lg.r_sigma_lt(p["x"], p["q2"])
        assert v == pytest.approx(f * r / (1.0 + r), rel=1e-14)
    fl98, _ = rmod.tree_fl(lambda: lg.ToyF2(), "r1998", pts)
    for p, v, f in zip(pts, fl98, f2):
        r = lg.r1998(p["x"], p["q2"])
        assert v == pytest.approx(f * r / (1.0 + r), rel=1e-14)
        # below the support r1998 is frozen at x = 0.005
        assert lg.r1998(p["x"], p["q2"]) == lg.r1998(0.005, p["q2"])


def test_h1_cells_block_by_name_without_mstw(rmod, monkeypatch):
    real = rmod.f2_backend

    def fake(name):
        if name == "mstw":
            return None, "MstwSF", "environment: lipolgen built without the PYTHIA 8 tier (MstwSF)"
        return real(name)

    monkeypatch.setattr(rmod, "f2_backend", fake)
    m = rmod.measure_h1()
    for b in rmod.BACKENDS:
        assert m["backends"][b] == dict(status="blocked", reason=(
            "environment: lipolgen built without the PYTHIA 8 tier (MstwSF)"))
    assert m["f2_backends"]["toy"]["status"] == "ran"
    assert "r_display" not in m


# --------------------------------------------- t1_r_world: the measurement

def test_r_chi2_recomputed_from_the_file(rmod, rrep):
    """The verdict cell's chi2, p and pulls recomputed here, independently."""
    pts = rmod.read_vendored("nmc")["body"]["points"]
    chi2 = sum(((0.18 / (1.0 + p["q2"] / 50.0) - p["value"]) ** 2)
               / (p["stat"] ** 2 + p["sys"] ** 2) for p in pts)
    c = rrep["detail"]["nmc"]["backends"]["r_sigma_lt"]
    assert c["chi2"] == pytest.approx(chi2, rel=1e-13)
    assert c["p"] == pytest.approx(_sf(chi2, 10), rel=1e-9)
    # the recorded fully-correlated reading is the profiled-nuisance chi2
    np = pytest.importorskip("numpy")
    r = np.array([0.18 / (1.0 + p["q2"] / 50.0) - p["value"] for p in pts])
    cov = np.diag([p["stat"] ** 2 for p in pts]) + np.outer(
        [p["sys"] for p in pts], [p["sys"] for p in pts])
    assert c["correlated"]["chi2"] == pytest.approx(float(r @ np.linalg.solve(cov, r)),
                                                    rel=1e-10)


_PINS = {  # (sub-row, backend): (chi2, ndf, status, recorded correlated chi2 or None)
    ("nmc", "r_sigma_lt"): (12.543626, 10, "pass", 60.215962),
    ("nmc", "r1998"): (6.2021095, 10, "pass", 21.281325),
    ("nmc", "r1998_unclipped"): (6.1910616, 10, "pass", 21.220963),
    ("r1990", "r_sigma_lt"): (209.75779, 68, "fail", None),
    ("r1990", "r1998"): (63.051159, 68, "pass", None),
    ("r1990", "r1998_unclipped"): (63.051159, 68, "pass", None),
    ("e140", "r_sigma_lt"): (29.657126, 9, "fail", 29.459187),
    ("e140", "r1998"): (12.714682, 9, "pass", 12.06818),
    ("e140", "r1998_unclipped"): (12.714682, 9, "pass", 12.06818),
}


@pytest.mark.parametrize("cell", sorted(_PINS))
def test_r_cells_pinned(rrep, cell):
    chi2, n, status, cor = _PINS[cell]
    c = rrep["detail"][cell[0]]["backends"][cell[1]]
    assert c["n"] == n and c["status"] == status
    assert c["chi2"] == pytest.approx(chi2, rel=1e-6)
    if cor is None:
        assert c["correlated"] is None
    else:
        assert c["correlated"]["chi2"] == pytest.approx(cor, rel=1e-6)


def test_r_worst_points_pinned(rrep):
    d = rrep["detail"]
    c = d["nmc"]["backends"]["r_sigma_lt"]
    assert c["worst_pull"] == pytest.approx(-2.8093118, rel=1e-6)
    assert (c["worst_x"], c["worst_q2"]) == (0.0045, 1.38)
    c = d["r1990"]["backends"]["r_sigma_lt"]
    assert c["worst_pull"] == pytest.approx(-4.0413414, rel=1e-6)
    assert (c["worst_x"], c["worst_q2"]) == (0.175, 1.25)
    c = d["e140"]["backends"]["r_sigma_lt"]
    assert (c["worst_x"], c["worst_q2"]) == (0.2, 1.0)


def test_h1_pinned_with_mstw_f2(rrep):
    ok, why = _have_mstw()
    if not ok:
        pytest.skip(why)
    m = rrep["detail"]["h1"]
    assert m["f2_verdict"] == "mstw"
    for b, chi2, cor, per_bin in (("r_sigma_lt", 8.8180530, 13.370260, 38.703918),
                                  ("r1998", 8.1812396, 13.429308, 37.650653),
                                  ("r1998_unclipped", 8.1052913, 13.527882, 37.477525)):
        c = m["backends"][b]
        assert (c["n"], c["status"]) == (8, "pass")
        assert c["chi2"] == pytest.approx(chi2, rel=1e-6)
        assert c["correlated"]["chi2"] == pytest.approx(cor, rel=1e-6)
        assert (c["per_bin"]["chi2"], c["per_bin"]["n"]) == (pytest.approx(per_bin, rel=1e-6), 25)
    assert [round(d["r"], 6) for d in m["r_display"]][:2] == [0.193227, 0.061883]


def test_h1_recorded_toyf2_variant_pinned(rrep):
    """The shipped generator's own proton F_L (ToyF2 F2p; no PYTHIA needed)."""
    cells = rrep["detail"]["h1"]["f2_backends"]["toy"]["cells"]
    assert cells["r_sigma_lt"]["chi2"] == pytest.approx(6.9606813, rel=1e-6)
    assert cells["r1998"]["chi2"] == pytest.approx(6.4793472, rel=1e-6)


def test_h1_ct18nlo_variant_reported_not_dropped(rrep):
    fb = rrep["detail"]["h1"]["f2_backends"]["ct18nlo"]
    if _lhapdf_missing("CT18NLO"):
        assert fb["status"] == "blocked" and fb["reason"].startswith("environment: ")
    else:
        assert fb["status"] == "ran"
        assert fb["cells"]["r_sigma_lt"]["chi2"] == pytest.approx(8.4836604, rel=1e-6)


def test_r_verdict_follows_from_the_numbers(rmod, rrep):
    for key in rmod.SUBROWS:
        for b in rmod.BACKENDS:
            c = rrep["detail"][key]["backends"][b]
            if c.get("status") == "blocked":
                continue
            assert c["p"] == pytest.approx(_sf(c["chi2"], c["n"]), rel=1e-9)
            assert c["status"] == ("pass" if c["p"] >= 0.01 else "fail")
    head = rrep["detail"]["nmc"]["backends"]["r_sigma_lt"]
    assert rrep["status"] == head["status"] == "pass"
    assert head["p"] == pytest.approx(0.25031862, rel=1e-6)
    # the verdict is the PRIMARY cell alone: the closure FAILs of r_sigma_lt
    # do not fail the row, and a synthetic NMC that r_sigma_lt misses would
    assert rrep["detail"]["r1990"]["backends"]["r_sigma_lt"]["status"] == "fail"
    pts = rmod.load_points("nmc")
    for p in pts:
        p["value"] += 0.4
    s = rmod.stats(pts, [lg.r_sigma_lt(p["x"], p["q2"]) for p in pts])
    assert s["p"] < 0.01


# ------------------------------------------------ t1_r_world: the bias table

def test_bias_cells_pinned_and_the_kernel_is_the_one_over_one_plus_r_law(rrep):
    cells = rrep["bias"]["cells"]
    assert [(c["label"], c["x"], c["q2"]) for c in cells] == [
        ("spot 2", 0.0112, 1.14), ("spot 1", 0.0282, 1.14),
        ("spot 3", 0.0708, 3.13), ("spot 4", 0.141, 14.3)]
    for c, ratio, kern in zip(cells, (1.1847575, 1.1779585, 1.0762386, 0.9563193),
                              (0.84409212, 0.84893202, 0.92916548, 1.0456619)):
        assert c["ratio"] == pytest.approx(ratio, rel=1e-6)
        assert c["kernel_ratio"] == pytest.approx(kern, rel=1e-7)
        # the kernel IS the closed form a2 ~ 1/[1 + 2(1-y)(1+R)/y^2] ...
        assert c["kernel_ratio"] == pytest.approx(c["closed_form"], rel=1e-12)
        # ... and (1 + R_lt)/(1 + R_98) up to O(y^2) at these y (0.010-0.026)
        assert abs(c["kernel_ratio"] - 1.0 / c["ratio"]) < 5e-5
        assert not c["extrapolated"] and c["y"] < 0.03
    assert rrep["bias"]["kernel"]["s"] == pytest.approx(3980.0, rel=1e-12)


def test_bias_grid_marks_and_freezes_the_extrapolated_region(rmod, rrep):
    g = rrep["bias"]["grid"]
    assert len(g) == len(rmod.BIAS_GRID_X) * len(rmod.BIAS_GRID_Q2)
    for q2 in rmod.BIAS_GRID_Q2:
        row = {c["x"]: c for c in g if c["q2"] == q2}
        for x in (0.001, 0.002, 0.003):
            assert row[x]["extrapolated"]
            assert row[x]["ratio"] == row[0.005]["ratio"]      # frozen at the edge
        assert not row[0.005]["extrapolated"]
    s = rrep["bias"]["summary"]
    assert s["in_support_min"] == pytest.approx(0.92066746, rel=1e-7)
    assert s["in_support_max"] == pytest.approx(1.1918173, rel=1e-7)
    assert s["extrapolated_min"] == pytest.approx(1.0036624, rel=1e-7)


# -------------------------------------------------- t1_hermes_f2d: design

def test_hermes_theory_is_the_isoscalar_nucleon_and_the_kernel(hmod, hrep):
    toy = lg.ToyF2()
    pts = hrep["points"]["deuteron"]
    th = hrep["configs_detail"]["toy"]["experiments"]["deuteron"]["theory"]
    for p, t in zip(pts, th):
        assert t == 0.5 * (toy.f2p(p["x"], p["q2"]) + toy.f2n(p["x"], p["q2"]))
    thp = hrep["configs_detail"]["toy"]["experiments"]["proton"]["theory"]
    for p, t in zip(hrep["points"]["proton"], thp):
        assert t == toy.f2p(p["x"], p["q2"])
    for c in hrep["configs_detail"].values():
        if c["status"] != "blocked":
            assert c["kernel_max_rel_diff"] < 1e-14


def test_hermes_covariance_chi2_is_the_profiled_normalisation(hmod, hrep):
    np = pytest.importorskip("numpy")
    e = hrep["configs_detail"]["toy"]["experiments"]["deuteron"]
    pts = [hrep["points"]["deuteron"][j] for j in e["selected"]]
    th = [e["theory"][j] for j in e["selected"]]
    r = np.array([p["value"] - t for p, t in zip(pts, th)])
    n = np.array([0.075 * p["value"] for p in pts])
    cov = np.diag([p["err_pp"] ** 2 for p in pts]) + np.outer(n, n)
    assert e["chi2"] == pytest.approx(float(r @ np.linalg.solve(cov, r)), rel=1e-10)
    # the per-bin error is the five components in quadrature, from percent
    p0 = hrep["points"]["deuteron"][0]
    assert p0["err_pp"] == pytest.approx(
        0.164 * math.sqrt(1.41 ** 2 + 0.69 ** 2 + 3.14 ** 2 + 5.36 ** 2 + 0.36 ** 2) / 100,
        rel=1e-12)


def test_hermes_dis_cut_counts(hrep):
    c = hrep["configs_detail"]
    assert c["toy"]["experiments"]["deuteron"]["n"] == 68
    assert c["toy"]["experiments"]["proton"]["n"] == 68
    if c["ct18nlo"]["status"] != "blocked":
        assert c["ct18nlo"]["experiments"]["deuteron"]["n"] == 45   # Q2 >= 1.677


# ------------------------------------------- t1_hermes_f2d: the measurement

def test_hermes_toy_pinned_as_a_recorded_fail(hrep):
    d = hrep["configs_detail"]["toy"]["experiments"]
    assert d["deuteron"]["chi2"] == pytest.approx(5303.4587, rel=1e-6)
    assert d["deuteron"]["norm_pull"] == pytest.approx(0.69046757, rel=1e-6)
    assert d["deuteron"]["chi2_pp"] == pytest.approx(5718.4490, rel=1e-6)
    assert d["deuteron"]["chi2_quad"] == pytest.approx(455.16476, rel=1e-6)
    assert d["deuteron"]["mean_ratio"] == pytest.approx(0.96543688, rel=1e-7)
    assert (d["deuteron"]["worst_x"], d["deuteron"]["worst_q2"]) == (0.49, 6.026)
    assert d["proton"]["chi2"] == pytest.approx(5678.4736, rel=1e-6)
    assert d["deuteron"]["status"] == d["proton"]["status"] == "fail"
    assert hrep["status"] == "fail"


def test_hermes_mstw_pinned(hrep):
    ok, why = _have_mstw()
    if not ok:
        pytest.skip(why)
    d = hrep["configs_detail"]["mstw"]["experiments"]
    assert d["deuteron"]["chi2"] == pytest.approx(318.86840, rel=1e-6)
    assert d["deuteron"]["chi2_quad"] == pytest.approx(31.430255, rel=1e-6)
    assert d["deuteron"]["mean_ratio"] == pytest.approx(1.0126662, rel=1e-7)
    assert d["proton"]["chi2"] == pytest.approx(247.79117, rel=1e-6)
    assert d["deuteron"]["status"] == "fail"


def test_hermes_ct18nlo_reported_not_dropped(hrep):
    c = hrep["configs_detail"]["ct18nlo"]
    if not getattr(lg, "HAVE_LHAPDF", False):
        assert c["status"] == "blocked"
        assert c["reason"] == ("environment: lipolgen built without the LHAPDF "
                               "tier (HAVE_LHAPDF false)")
    elif _lhapdf_missing("CT18NLO"):
        assert c["status"] == "blocked"
        assert c["reason"] == "environment: LHAPDF set CT18NLO not installed"
    else:
        d = c["experiments"]["deuteron"]
        assert d["chi2"] == pytest.approx(1788.2315, rel=1e-6)
        assert c["q2_floor"] == pytest.approx(1.677025, rel=1e-6)


def test_hermes_verdict_follows_from_the_numbers(hmod, hrep):
    for c in hrep["configs_detail"].values():
        if c["status"] == "blocked":
            continue
        for t in hmod.TARGETS:
            e = c["experiments"][t]
            assert e["p"] == pytest.approx(_sf(e["chi2"], e["n"]), rel=1e-9) or e["p"] == 0.0
            assert e["status"] == ("pass" if e["p"] >= 0.01 else "fail")
        assert c["status"] == c["experiments"]["deuteron"]["status"]
    assert hrep["status"] == hrep["configs_detail"]["toy"]["status"]


def test_hermes_b1_inversion_reads_per_nucleon(hmod, hrep):
    inv = hrep["b1_inversion"]
    ratios = [b["ratio"] for b in inv["bins"]]
    assert ratios == pytest.approx([0.92856163, 1.0385010, 0.98582960, 0.95780597,
                                    1.0056957, 1.0057174], rel=1e-7)
    assert inv["reading"] == "per nucleon"
    assert inv["geometric_mean"] == pytest.approx(0.98637018, rel=1e-7)
    b = inv["bins"][0]
    assert b["f1_implied"] == pytest.approx(-(2.0 / 3.0) * 11.2 / -1.06, rel=1e-14)
    assert (b["hermes_x"], b["hermes_q2"], b["hermes_f2d"]) == (0.0117, 0.477, 0.226)
    # the declared rule, both ways and in between
    assert hmod.classify([2.0 * v for v in ratios]) == "per deuteron"
    assert hmod.classify(ratios + [1.4]) == "undetermined"


# --------------------------------------------- lines, exit codes, figures

def test_r_world_prints_one_report_row_and_every_subrow(rmod, capsys):
    rep = rmod.run(verbose=True)
    out = capsys.readouterr().out
    rows = [l for l in out.splitlines() if l.startswith("REPORT |")]
    assert len(rows) == 1
    f = [x.strip() for x in rows[0].split("|")]
    assert f[1] == "t1_r_world" and len(f) == 6 and f[5] == "PASS"
    assert f[2].startswith("r_sigma_lt (shipped R) on NMC R: chi2/ndf = 12.54/10, p = 0.25")
    sub = [l for l in out.splitlines() if l.startswith("SUBROW |")]
    assert len(sub) == 8
    assert [l.split("|")[1].strip() for l in sub[:4]] == [
        "t1_r_world[nmc]", "t1_r_world[h1]", "t1_r_world[r1990]", "t1_r_world[e140]"]
    assert sub[0].split("|")[2].strip().startswith("PRIMARY:")
    assert sub[1].split("|")[2].strip().startswith("EXTRAPOLATION:")
    assert all(l.endswith("| RECORDED") for l in sub[4:])
    assert "np." not in out
    assert rep["runtime_s"] < 30.0


def test_hermes_prints_one_report_row_and_every_subrow(hmod, capsys):
    rep = hmod.run(verbose=True)
    out = capsys.readouterr().out
    rows = [l for l in out.splitlines() if l.startswith("REPORT |")]
    assert len(rows) == 1 and rows[0].endswith("| FAIL")
    assert "-> PER NUCLEON" in rows[0]
    sub = [l for l in out.splitlines() if l.startswith("SUBROW |")]
    blocked = [n for n, c in rep["configs_detail"].items() if c["status"] == "blocked"]
    assert len(sub) == 2 * (3 - len(blocked)) + len(blocked) + 1
    assert sub[-1].startswith("SUBROW | t1_hermes_f2d[b1_inversion]")
    assert "np." not in out
    assert rep["runtime_s"] < 30.0


def test_reports_normalise_for_make_report(rrep, hrep):
    sys.path.insert(0, _BENCH)
    import _report
    row = _report.normalise_row(rrep, "t1_r_world")
    assert row["status"] == "pass" and not row["notes"]
    names = [s["name"] for s in row["subrows"]]
    assert names[:4] == ["nmc", "h1", "r1990", "e140"] and len(names) == 8
    row = _report.normalise_row(hrep, "t1_hermes_f2d")
    assert row["status"] == "fail" and row["subrows"][-1]["name"] == "b1_inversion"
    assert ("toy/deuteron", "fail") in [(s["name"], s["status"]) for s in row["subrows"]]


@pytest.mark.parametrize("harness,last", [("t1_r_world", "| PASS"),
                                          ("t1_hermes_f2d", "| FAIL")])
def test_harness_exits_0_on_a_run(harness, last):
    res = _run_script(os.path.join(_BENCH, harness + ".py"))
    assert res.returncode == 0, res.stderr
    rows = [l for l in res.stdout.splitlines() if l.startswith("REPORT | %s |" % harness)]
    assert len(rows) == 1 and rows[0].endswith(last)


def _altered(raw, how):
    """One byte (or one number) of a vendored file changed."""
    if how == "rows":                       # hermes2005_b1d_table2.json: no body
        out = raw.replace(b"0.012", b"0.013", 1)
    else:
        out = _flip_first_value(raw, how)
    assert out != raw
    return out


@pytest.mark.parametrize("harness,altered,how", [
    ("t1_r_world", "hepdata_ins360765_e140_r.json", b'"points": ['),
    ("t1_r_world", "hepdata_ins424154_nmc_r.json", b'"points": ['),
    ("t1_hermes_f2d", "hepdata_ins894309_hermes_f2.json", b'"deuteron": {'),
    ("t1_hermes_f2d", "hermes2005_b1d_table2.json", "rows")])
def test_harness_exits_2_on_an_altered_file(tmp_path, harness, altered, how):
    """The harness copied next to a data/ holding every file it reads, one of
    them altered: it refuses it and exits 2, printing no REPORT line."""
    files = ([v[0] for v in _load("t1_r_world").FILES.values()]
             if harness == "t1_r_world" else
             ["hepdata_ins894309_hermes_f2.json", "hermes2005_b1d_table2.json",
              "whitlow1990_r1990.json"])
    d = tmp_path / harness
    (d / "data").mkdir(parents=True)
    shutil.copy(os.path.join(_BENCH, harness + ".py"), d / (harness + ".py"))
    for f in files:
        raw = open(os.path.join(_DATA, f), "rb").read()
        (d / "data" / f).write_bytes(_altered(raw, how) if f == altered else raw)
    script = str(d / (harness + ".py"))
    res = _run_script(script)
    assert res.returncode == 2
    assert "sha256" in res.stderr and "REPORT |" not in res.stdout


@pytest.mark.parametrize("harness", ["t1_r_world", "t1_hermes_f2d"])
def test_plot_writes_a_png(harness, rrep, hrep, tmp_path):
    pytest.importorskip("matplotlib")
    mod = _load(harness)
    rep = rrep if harness == "t1_r_world" else hrep
    path = mod.plot(rep, out_dir=str(tmp_path))
    assert path == str(tmp_path / (harness + ".png"))
    raw = open(path, "rb").read()
    assert len(raw) > 10000 and raw[:8] == b"\x89PNG\r\n\x1a\n"


def test_chi2_sf_fallback_is_t1s_and_scipys(rmod, hmod):
    f = _load("t1_nmc_f2d").chi2_sf_fallback
    for m in (rmod, hmod):
        assert m.chi2_sf_fallback.__code__.co_code == f.__code__.co_code
    pytest.importorskip("scipy", reason="scipy (the reference for the fallback)")
    for chi2, n in ((12.543626, 10), (6.2021095, 10), (8.818053, 8), (209.75779, 68),
                    (29.657126, 9), (318.8684, 68), (60.215962, 10)):
        assert rmod.chi2_sf_fallback(chi2, n) == pytest.approx(_sf(chi2, n), rel=1e-10)
