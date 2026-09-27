# SPDX-License-Identifier: GPL-3.0-or-later
"""T4 row "PYTHIA 8 ep closure, broadened" and T6 row "HepMC3 convention
conformance" -- the harnesses under validation/benchmarks/, exercised as
tests.

  t4_pythia_ep_closure.py   LiPolGen's PythiaBridge vs stock PYTHIA 8 e p at
                            two (x, Q2) windows (01_generators.md D-1).
                            EXPENSIVE by kind: the full-size row runs only
                            with LIPOLGEN_BENCH=1; the config checks and a
                            200-event smoke run are always on (~1.5 s).
  t6_hepmc3_convention.py   eight lipolgen-run HepMC3 files read back with
                            pyhepmc and checked clause by clause against
                            docs/HEPMC3_CONVENTION.md (~4.5 s, always on).

A harness that DISAGREES with the tree is not a failing test: the tests pin
what the harness MEASURES and that its verdict follows from the numbers.
T6 is a recorded FAIL today (the inclusive record does not carry the (A-1)
remnant the document's conservation section implies), and the test pins
exactly that.  No test here moves a tree default.
"""
import hashlib
import importlib.util
import math
import os
import re
import subprocess
import sys

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


def _have_pythia8_module():
    try:
        import pythia8                                     # noqa: F401
    except ImportError:
        return False
    return True


def _skip_without_pythia_tier(what):
    """Skip naming EXACTLY what is missing."""
    if not getattr(lg, "HAVE_PYTHIA8", False):
        pytest.skip("environment: lipolgen built without PYTHIA 8 "
                    "(HAVE_PYTHIA8 false); %s needs the PythiaBridge" % what)
    if not _have_pythia8_module():
        pytest.skip("environment: pythia8 Python module not importable; %s "
                    "needs it" % what)


# =================================================================== T4

@pytest.fixture(scope="module")
def t4():
    return _load("t4_pythia_ep_closure")


def test_t4_config_is_well_formed(t4):
    """The row's declared constants: D-1's 20 %, two windows, the observable
    list, the eta bins and the compare threshold, all fixed in the module."""
    assert t4.NAME == "t4_pythia_ep_closure"
    assert t4.EXPENSIVE is True
    assert t4.TOL == 0.20
    assert t4.N_EVENTS >= 10000 and t4.N_SUB == 4 and t4.MIN_STOCK_COUNTS == 100
    assert [w["label"] for w in t4.WINDOWS] == ["W1", "W2"]
    w1, w2 = t4.WINDOWS
    assert (w1["x"], w1["q2"]) == ((0.005, 0.10), (4.0, 30.0))   # D-1's own
    assert (w2["x"], w2["q2"]) == ((0.10, 0.50), (30.0, 300.0))
    for w in t4.WINDOWS:
        assert 0 < w["x"][0] < w["x"][1] <= 1 and 0 < w["q2"][0] < w["q2"][1]
    edges = t4.ETA_EDGES
    assert edges[0] == -math.inf and edges[-1] == math.inf
    assert all(a < b for a, b in zip(edges[:-1], edges[1:]))
    assert t4.OBSERVABLES[:3] == ["n_ch", "sum_empz", "sum_pt_ch"]
    assert len(t4.OBSERVABLES) == 3 + len(edges) - 1
    assert t4.PYTHIA_PINNED == "8.317"


def test_t4_stock_settings_are_d1_verbatim(t4):
    """The stock run IS tests/test_pythia.cpp's D-1 configuration, plus
    exactly one setting -- the bridge's own SpaceShower:QEDshowerByQ = off
    (docs/PYTHIA_BRIDGE.md sec. 7) -- and the per-window Q2Min / seed."""
    src = open(os.path.join(_ROOT, "tests", "test_pythia.cpp")).read()
    block = src[src.index("Pythia8::Pythia stock("):]
    block = block[block.index("{"):block.index("})")]
    d1 = re.findall(r'"([^"]+)"', block)
    assert len(d1) == 18
    for s in t4.STOCK_SETTINGS_D1:
        assert s in d1, s
    rest = [s for s in d1 if s not in t4.STOCK_SETTINGS_D1]
    assert sorted(re.sub(r"=.*", "", s).strip() for s in rest) == \
        ["PhaseSpace:Q2Min", "Random:seed"]
    assert t4.STOCK_SETTINGS_ADDED == ("SpaceShower:QEDshowerByQ = off",)
    bridge_doc = open(os.path.join(_ROOT, "docs", "PYTHIA_BRIDGE.md")).read()
    assert "SpaceShower:QEDshowerByQ = off" in bridge_doc
    st = t4.stock_settings(t4.WINDOWS[1], 7, 100.0044, 10.0)
    assert "PhaseSpace:Q2Min = 30" in st and "Random:seed = 8" in st
    assert "SpaceShower:QEDshowerByQ = off" in st


def test_t4_normalisation_rules(t4):
    """Sec. 8, as they apply to this row: rule 1 -- an UNPOLARISED spin-1/2
    target (P_e = 0 in every category, flat populations); rules 2 and 3 --
    A = 1, a FREE proton, so per nucleus == per nucleon and there is no
    average nucleon to divide by.  The grid is the window."""
    w = t4.WINDOWS[0]
    cfg = t4.make_config(w, 10, 1)
    assert cfg.isotope == "p"
    beams = lg.default_configs(cfg.isotope)[cfg.beam_config]
    assert (beams.ion.A, beams.ion.Z) == (1, 1)
    assert beams.electron_energy == 10.0 and beams.ion_momentum_per_nucleon == 100.0
    assert (cfg.grid.x_min, cfg.grid.x_max) == w["x"]
    assert (cfg.grid.q2_min, cfg.grid.q2_max) == w["q2"]
    plan = t4.make_plan()
    for c in plan.categories:
        assert c.j == 0.5 and c.pe == 0.0
        assert c.populations[0] == pytest.approx(0.5, abs=1e-12)
    cuts = t4.scenario_cuts(cfg)
    for w in t4.WINDOWS:                   # both windows inside the Scenario
        assert w["x"][1] <= cuts["x_max"]
        assert w["q2"][0] >= cfg.scenario.q2_min


def test_t4_kinematics_and_observables(t4):
    """The stock side's (x, Q2, y, W2) from a scattered electron, against the
    tree's own dis_scattered_electron; and the per-event observables."""
    if not getattr(lg, "HAVE_PYTHIA8", False):
        pytest.skip("environment: lipolgen built without PYTHIA 8 "
                    "(HAVE_PYTHIA8 false); lg.dis_scattered_electron lives there")
    m = 0.938272088
    p = (math.sqrt(100.0 ** 2 + m * m), 0.0, 0.0, 100.0)
    for x, q2, phi in ((0.01, 5.0, 0.3), (0.3, 120.0, 2.0)):
        kp, y = lg.dis_scattered_electron(10.0, lg.Vec4(*p), x, q2, phi)
        got = t4.dis_kinematics((10.0, 0.0, 0.0, -10.0), p,
                                (kp.e, kp.px, kp.py, kp.pz))
        assert got[0] == pytest.approx(x, rel=1e-10)
        assert got[1] == pytest.approx(q2, rel=1e-10)
        assert got[2] == pytest.approx(y, rel=1e-10)
        assert got[3] == pytest.approx(m * m + q2 * (1 - x) / x, rel=1e-9)
    # one charged particle at eta = 1, one neutral, one charged with pT = 0
    sh = math.sinh(1.0)
    row = t4.hfs_row([(math.sqrt(1 + sh * sh), 1.0, 0.0, sh, +1.0),
                      (2.0, 0.0, 0.0, -2.0, 0.0),
                      (50.0, 0.0, 0.0, 49.9, -1.0)])
    assert row[0] == 2.0                                   # n_ch
    assert row[1] == pytest.approx((math.sqrt(1 + sh * sh) - sh) + 4.0 + 0.1)
    assert row[2] == pytest.approx(1.0)                    # scalar sum pT_ch
    assert row[3:] == [0.0, 0.0, 1.0, 0.0, 1.0]            # eta 1 -> [0,2); +inf -> last


def _t4_verdict_follows(rep, t4):
    comp = [s for s in rep["subrows"] if s["compared"]]
    for s in rep["subrows"]:
        if s["observable"].startswith("n_ch_eta"):
            assert s["compared"] == (s["n_stock_particles"] >= t4.MIN_STOCK_COUNTS)
        else:
            assert s["compared"]
        if s["compared"]:
            assert s["ratio"] == pytest.approx(s["bridge"] / s["stock"], rel=1e-12)
            assert s["status"] == ("pass" if abs(s["ratio"] - 1) <= t4.TOL else "fail")
    assert rep["status"] == ("pass" if all(s["status"] == "pass" for s in comp)
                             else "fail")
    worst = max(comp, key=lambda s: s["distance"])
    assert rep["worst"]["distance"] == worst["distance"]


def test_t4_smoke_small_n(t4):
    """200 events per window per generator: the whole row end to end, fast.
    The verdict at this N is NOT the row's; what is pinned is that the row
    is well formed, follows from its numbers, and ran the configuration its
    header claims (one PYTHIA library, PYTHIA's own PDF, protons only)."""
    _skip_without_pythia_tier("the T4 smoke run")
    import pythia8
    rep = t4.run(verbose=False, n_events=200)
    for k in ("name", "generator", "reference", "tolerance", "status"):
        assert k in rep
    assert rep["status"] in ("pass", "fail")
    assert len(rep["subrows"]) == 2 * len(t4.OBSERVABLES)
    _t4_verdict_follows(rep, t4)
    cfg = rep["config"]
    ver = float(pythia8.Pythia("", False).settings.parm("Pythia:versionNumber"))
    assert cfg["pythia_version_used"] == ver
    assert cfg["pythia_version_pinned"] == "8.317"
    if cfg["pythia_libraries"] is not None:
        assert len(cfg["pythia_libraries"]) == 1     # ONE library, both sides
    assert "SpaceShower:QEDshowerByQ = off" in cfg["bridge_applied_settings"]
    for w in rep["windows"]:
        b, s = w["bridge"], w["stock"]
        assert b["pdf_setting_applied"] == []        # the bridge keeps PYTHIA's PDF
        assert b["A"] == 1 and b["struck_all_proton"]  # rules 2 / 3: a free proton
        assert b["n_generated"] == 200 and s["n_accepted"] == 200
        assert s["pdf_pset"] == cfg["pdf_pset"]
        assert b["max_rescale_dev"] < 1e-6
        assert w["match"]["n_stock_matched"] + w["match"]["n_stock_dropped"] == 200
    assert rep["runtime_s"] < 30.0


def test_t4_blocked_row_names_the_environment(t4, monkeypatch, capsys):
    reason = "environment: pythia8 Python module not importable (stock side)"
    monkeypatch.setattr(t4, "environment", lambda: reason)
    rep = t4.run(verbose=True)
    assert rep["status"] == "blocked" and rep["reason"] == reason
    out = capsys.readouterr().out
    assert "REPORT | t4_pythia_ep_closure |" in out and "BLOCKED" in out


@pytest.mark.parametrize("harness", ["t4_pythia_ep_closure",
                                     "t6_hepmc3_convention"])
def test_t4_t6_exit_2_when_lipolgen_is_not_importable(harness, tmp_path):
    """A missing or broken lipolgen is the harness broken (exit 2, README
    convention 6), not a BLOCKED row: only a build tier that is off or the
    pythia8 / pyhepmc module is an environment reason (until 2026-09-27
    both harnesses printed a BLOCKED row and exited 0 here)."""
    fake = tmp_path / "lipolgen"
    fake.mkdir()
    (fake / "__init__.py").write_text(
        'raise ImportError("simulated: lipolgen not importable")\n')
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(tmp_path)] + [p for p in [env.get("PYTHONPATH")] if p])
    p = subprocess.run([sys.executable, os.path.join(_BENCH, harness + ".py")],
                       env=env, capture_output=True, text=True, timeout=300)
    assert p.returncode == 2, (p.stdout, p.stderr)
    assert "simulated: lipolgen not importable" in p.stderr
    assert "REPORT |" not in p.stdout


# The full-size row: ~10 s measured, but a generator-vs-generator sample is an
# EXPENSIVE row by kind (BENCHMARK_PLAN.md sec. 3).
T4_PINNED = {  # measured 2026-09-26, PYTHIA 8.312, N = 20 000 per side/window
    ("W1", "n_ch"): (8.0267, 8.0539),
    ("W1", "sum_empz"): (2.7229, 2.7405),
    ("W1", "sum_pt_ch"): (3.8271, 3.8492),
    ("W2", "n_ch"): (7.2606, 7.2428),
    ("W2", "sum_empz"): (2.0012, 2.0060),
    ("W2", "sum_pt_ch"): (6.1278, 6.1313),
}


@pytest.mark.bench_expensive
def test_t4_full_row_passes_and_is_pinned(t4):
    if os.environ.get("LIPOLGEN_BENCH") != "1":
        pytest.skip("expensive benchmark row (a generator-vs-generator "
                    "sample): set LIPOLGEN_BENCH=1 to run it")
    _skip_without_pythia_tier("the T4 row")
    rep = t4.run(verbose=False)
    _t4_verdict_follows(rep, t4)
    assert rep["status"] == "pass"
    by = {(s["window"], s["observable"]): s for s in rep["subrows"]}
    # rel 3 %: several times each mean's statistical error (<= 0.7 %), loose
    # enough for the pinned PYTHIA 8.317 against the 8.312 that measured it.
    for key, (stock, bridge) in T4_PINNED.items():
        assert by[key]["stock"] == pytest.approx(stock, rel=0.03), key
        assert by[key]["bridge"] == pytest.approx(bridge, rel=0.03), key
        assert by[key]["ratio"] == pytest.approx(bridge / stock, abs=0.03), key
    assert sum(1 for s in rep["subrows"] if s["compared"]) == 15
    assert not by[("W2", "n_ch_eta[-inf,-2)")]["compared"]
    assert rep["worst"]["distance"] == pytest.approx(0.023, abs=0.03)
    for w in rep["windows"]:
        assert w["match"]["n_sub_kept"] == 16 and w["match"]["n_stock_dropped"] == 0
        assert w["stock"]["n_accepted"] == t4.N_EVENTS
        assert w["bridge"]["n_vetoed"] == 0 and w["bridge"]["struck_all_proton"]


# =================================================================== T6

@pytest.fixture(scope="module")
def t6():
    return _load("t6_hepmc3_convention")


def _skip_without_hepmc3():
    if not getattr(lg, "HAVE_HEPMC3", False):
        pytest.skip("environment: lipolgen built without HepMC3 "
                    "(HAVE_HEPMC3 false); T6 needs the HepMC3 writer")
    try:
        import pyhepmc                                     # noqa: F401
    except ImportError:
        pytest.skip("environment: pyhepmc not importable; T6 reads the files "
                    "with it")


@pytest.fixture(scope="module")
def t6_rep(t6):
    _skip_without_hepmc3()
    return t6.run(verbose=False)


def test_t6_clause_table_is_complete(t6):
    ids = [c[0] for c in t6.CLAUSES]
    assert ids == ["C%02d" % i for i in range(1, len(ids) + 1)]
    assert set(t6.CHECKS) == set(ids)
    n_lines = len(open(t6.DOC, encoding="utf-8").read().splitlines())
    for cid, lines, (aline, phrase), what in t6.CLAUSES:
        lo, hi = (int(v) for v in re.match(r"(\d+)-(\d+)", lines).groups())
        assert 1 <= lo <= aline <= hi <= n_lines, cid
        assert phrase.isascii() and what
    assert len(t6.ATTRIBUTES) == 24                 # doc lines 147-170
    assert len(t6.RUNS) == 8
    assert {r["kind"] for r in t6.RUNS} == {"inclusive", "alpha-tag", "d-tag",
                                           "coherent"}


def test_t6_doc_anchors_stand_and_sha_is_the_authored_one(t6):
    """The reference is the DOCUMENT.  Every clause's anchor phrase stands on
    its cited line, and the document is the one the clauses were written
    against -- a doc edit fails HERE, visibly, until someone re-reads the
    clauses against it and updates DOC_SHA256_AT_AUTHORING."""
    st = t6.doc_state()
    moved = {c: a for c, a in st["anchors"].items() if not a["found"]}
    assert not moved, moved
    raw = open(t6.DOC, "rb").read()
    assert st["sha256"] == hashlib.sha256(raw).hexdigest()
    assert st["sha256"] == t6.DOC_SHA256_AT_AUTHORING, (
        "docs/HEPMC3_CONVENTION.md changed: re-read the T6 clauses against it, "
        "then update DOC_SHA256_AT_AUTHORING")


def test_t6_nuclear_code_parser(t6):
    assert t6.nuclear_code(1000030060) == (3, 6)
    assert t6.nuclear_code(1000020040) == (2, 4)
    assert t6.nuclear_code(1000010020) == (1, 2)
    assert t6.nuclear_code(1000030061) is None           # I != 0
    assert t6.nuclear_code(1010030060) is None           # L != 0
    assert t6.nuclear_code(1000040030) is None           # Z > A
    assert t6.nuclear_code(2212) is None


def test_t6_row_and_verdict_follow_from_the_clauses(t6, t6_rep):
    rep = t6_rep
    for k in ("name", "generator", "reference", "tolerance", "status"):
        assert k in rep
    assert rep["name"] == "t6_hepmc3_convention"
    subs = rep["subrows"]
    assert [s["clause"] for s in subs] == [c[0] for c in t6.CLAUSES]
    fails = [s for s in subs if s["status"] == "fail"]
    if fails:
        assert rep["status"] == "fail"
        assert rep["reason"].startswith(fails[0]["clause"])
        assert fails[0]["first_violation"] in rep["reason"]
    else:
        assert rep["status"] in ("pass", "blocked")
    for s in subs:
        assert s["status"] in ("pass", "fail", "blocked", "n/a")
        assert (s["status"] == "fail") == (s["n_violations"] > 0)
        assert s["anchor_found"]
    assert rep["doc_sha256"] == hashlib.sha256(open(t6.DOC, "rb").read()).hexdigest()
    assert rep["doc_sha256"][:16] in rep["reference"]


def test_t6_measured_row_is_pinned(t6, t6_rep):
    """Measured 2026-09-27, 100 events per file: all 23 clauses hold.  Until
    decision D1 (docs/open_items/run_2026-09-27/DECISIONS.md) C13 and C14
    failed on the four INCLUSIVE files, whose whole deficit was the (A-1)
    remnant USAGE.md sec. 2 says is not written; the document now states the
    inclusive balance per nucleon (docs/HEPMC3_CONVENTION.md lines 140-147)."""
    if not getattr(lg, "HAVE_PYTHIA8", False):
        pytest.skip("environment: lipolgen built without PYTHIA 8 "
                    "(HAVE_PYTHIA8 false); six of T6's eight files are T2")
    if not _have_pythia8_module():
        pytest.skip("environment: pythia8 Python module not importable; T6 "
                    "clause C14 reads PYTHIA's PDG charge table")
    rep = t6_rep
    assert rep["status"] == "pass"
    by = {s["clause"]: s for s in rep["subrows"]}
    assert all(s["status"] == "pass" for s in by.values()), \
        {c: s["status"] for c, s in by.items()}
    assert by["C12"]["n_checked"] == 200 and by["C22"]["files"] == 1
    assert by["C23"]["files"] == 7
    # every record balances to far inside the doc's 1e-9: 700 on status 1
    # alone, the 100 T0 inclusive events with X; no remnant-shaped deficit left
    c13 = rep["conservation"]
    assert c13["n_final"] == 700 and c13["n_x"] == 100 and c13["worst"] < 1e-12
    assert rep["charge"]["n_checked"] == 700 and rep["charge"]["n_unknown"] == 0
    assert not rep["remnant_diagnostic"] and not rep["remnant_charge_diagnostic"]
    assert rep["n_spin_weights_seen"] == [0]
    assert rep["runtime_s"] < 30.0


def test_t6_inclusive_balance_is_per_nucleon(t6):
    """D1: an inclusive record balances against e + the struck nucleon, a
    tagged / coherent one against the full beams."""
    e_ = {"parts": [dict(id=1, pid=11, status=4, p4=(10.0, 0.0, 0.0, -10.0)),
                    dict(id=2, pid=1000030060, status=4, p4=(600.0, 0, 0, 599.0)),
                    dict(id=3, pid=2212, status=3, p4=(100.0, 0, 0, 99.9)),
                    dict(id=4, pid=92, status=3, p4=(1.0, 0, 0, 0), prod=-2)],
          "verts": {-2: dict(ins=[3], outs=[4])}}
    inc = dict(kind="inclusive")
    got = t6._initial_state(e_, inc)
    assert got == [(10.0, 0.0, 0.0, -10.0), (100.0, 0, 0, 99.9)]
    assert t6._initial_state(e_, dict(kind="alpha-tag")) == [
        (10.0, 0.0, 0.0, -10.0), (600.0, 0, 0, 599.0)]


_TAGGED_T0 = dict(label="tagged-6Li-alpha-T0", ion=1000030060,
                  enum="TaggedLi6Alpha", kind="alpha-tag", cluster=1000010020,
                  t2=False, rc=False, argv=["--channel", "tagged-6Li-alpha"])


def _tamper(src, dst, fn):
    lines = open(src).read().splitlines()
    first = next(i for i, s in enumerate(lines) if s.startswith("E "))
    nxt = next((i for i, s in enumerate(lines) if s.startswith("E ") and i > first),
               len(lines))
    out = fn(lines, first, nxt)
    open(dst, "w").write("\n".join(out) + "\n")


def _violations(t6, run, path, cid):
    raw, events = t6.read_file(path)
    ctx = dict(n_events=len(events), version=lg.__version__, charge=None,
               charge_reason="not needed")
    return t6.CHECKS[cid](run, raw, events, ctx)[1]


def test_t6_checks_can_fail(t6, tmp_path):
    """Not tautologies: a T0 tagged file (no PYTHIA needed) passes every
    clause it is subject to, and one tampered line in its first event fails
    the clause that reads it."""
    _skip_without_hepmc3()
    path = str(tmp_path / "tag.hepmc")
    t6.generate(_TAGGED_T0, path, n_events=5, seed=3)
    for cid in t6.CHECKS:
        if cid == "C14":
            continue                       # needs the charge table; n/a at T0
        assert _violations(t6, _TAGGED_T0, path, cid) == [], cid

    def units(lines, a, b):
        return [s.replace("U GEV MM", "U MEV MM") if a <= i < b else s
                for i, s in enumerate(lines)]

    def drop_spin_j(lines, a, b):
        return [s for i, s in enumerate(lines)
                if not (a <= i < b and s.startswith("A 0 spin_J "))]

    def x_final(lines, a, b):
        out = []
        for i, s in enumerate(lines):
            f = s.split()
            if a <= i < b and f[0] == "P" and f[3] == "92":
                s = " ".join(f[:-1] + ["1"])
            out.append(s)
        return out

    def beam_e_mass(lines, a, b):
        out = []
        for i, s in enumerate(lines):
            f = s.split()
            if a <= i < b and f[0] == "P" and f[3] == "11" and f[-1] == "4":
                s = " ".join(f[:-2] + ["5.1100000000000000e-04", "4"])
            out.append(s)
        return out

    for fn, cids in ((units, ("C02",)), (drop_spin_j, ("C15",)),
                     (x_final, ("C08", "C13")), (beam_e_mass, ("C10",))):
        bad = str(tmp_path / (fn.__name__ + ".hepmc"))
        _tamper(path, bad, fn)
        for cid in cids:
            v = _violations(t6, _TAGGED_T0, bad, cid)
            assert v and v[0].startswith("event 0"), (fn.__name__, cid, v)


def test_t6_blocked_paths_name_the_environment(t6, monkeypatch):
    reason = "environment: pyhepmc not importable"
    monkeypatch.setattr(t6, "environment", lambda: dict(all=reason))
    rep = t6.run(verbose=False)
    assert rep["status"] == "blocked" and rep["reason"] == reason


def test_t6_missing_charge_table_blocks_only_c14(t6, monkeypatch):
    _skip_without_hepmc3()
    real = t6.environment()
    env = dict(real, charge="environment: pythia8 Python module (the PDG "
                            "charge table) not importable")
    monkeypatch.setattr(t6, "environment", lambda: env)
    rep = t6.run(verbose=False, n_events=5, runs=[_TAGGED_T0])
    by = {s["clause"]: s for s in rep["subrows"]}
    assert by["C14"]["status"] == "blocked"
    assert by["C14"]["reason"] == env["charge"]
    assert rep["status"] == "blocked" and rep["reason"] == env["charge"]
    assert all(s["status"] in ("pass", "n/a") for c, s in by.items() if c != "C14")
