# SPDX-License-Identifier: GPL-3.0-or-later
"""The REPORT.md generator (`validation/benchmarks/make_report.py`), its
shared helpers (`_report.py`), the sec. 7 change -> tier map (`tiers.py`) and
the LIPOLGEN_BENCH opt-in (`conftest.py`) -- BENCHMARK_PLAN.md sec. 3:
"writes one row of docs/benchmarking/REPORT.md ... The report is
regenerated, never hand-edited."

Everything except the last test runs always and needs neither lipolgen nor
LHAPDF: the harnesses here are FAKE modules written into `tmp_path`, one per
path through the generator (pass, fail, blocked, error, a missing LHAPDF
set, an import error, a `sys.exit`, a crashed child process, an opt-in row,
numpy scalars, sub-rows), plus a fake BENCHMARK_PLAN.md.  The last test
(`bench_expensive`, i.e. only with LIPOLGEN_BENCH=1) regenerates the real
report and requires the committed docs/benchmarking/REPORT.md / REPORT.json
to be what it gives (`make_report.py --check`), skipping -- and naming the
difference -- when the environment is not the one the committed report was
made in (an LHAPDF set installed since, another PYTHIA).

The generator judges nothing: no test here asserts a physics number; the
numbers below are the fake harnesses' own.
"""
import importlib.util
import io
import json
import os
import subprocess
import sys
import textwrap

import pytest

np = pytest.importorskip("numpy")

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_BENCH = os.path.join(_ROOT, "validation", "benchmarks")
_REPORT_DIR = os.path.join(_ROOT, "docs", "benchmarking")
_PLAN = os.path.join(_REPORT_DIR, "BENCHMARK_PLAN.md")


def _load(name, path=None):
    path = path or os.path.join(_BENCH, name + ".py")
    spec = importlib.util.spec_from_file_location("test_bench_report_" + name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def mr():
    return _load("make_report")


@pytest.fixture(scope="module")
def rp():
    return _load("_report")


@pytest.fixture(scope="module")
def tiers():
    return _load("tiers")


FAKE_ENV = dict(python="3.x", numpy="2.x", scipy="1.x", pyhepmc="2.x",
                lipolgen=dict(version="0.1.0", have_lhapdf=True,
                              have_pythia8=True),
                pythia="8.312", pythia_source="fake", pythia_pinned="8.317",
                lhapdf="6.5.6", hepmc3="3.03.01", lhapdf_search_path=["/x"],
                lhapdf_sets_present=[])
FAKE_NOW = "2026-09-27T00:00:00Z"
FAKE_GIT = ("abc1234", False)

FAKE_PLAN = """\
# plan

## 4. The first ten

| # | benchmark | tier | effort | what it buys | status |
|---|---|---|---|---|---|
| 1 | **the pass row** | T5 | hours | text | **wired** (`validation/benchmarks/t5_fake_pass.py`) |
| 2 | **resolved in code** | T2 | hours | text | **resolved** (not a harness) |
| 3 | **the fail row** | T1 | low | text | wired as `t1_fake_fail` |

## 5. Corrections

nothing

## 6. What cannot be benchmarked

| # | observable | closest proxy | distance |
|---|---|---|---|
| U-1 | b1(6Li) | b1d | A = 2 -> 6 |
| U-2 | a2 | a2(|t|) of the deuteron | 3.5x |

## 7. The map
"""

FAKES = {
    "t5_fake_pass": """
        NAME = "t5_fake_pass"
        def run(verbose=True):
            return dict(name=NAME, generator="residual 0.000e+00 over 8 modes",
                        reference="the table", tolerance=1e-12, status="pass")
    """,
    "t1_fake_fail": """
        import numpy as np
        def run(verbose=True):
            return dict(name="t1_fake_fail", generator="chi2/ndf = 12.5/6",
                        reference="data", tolerance="p >= 0.01", status="fail",
                        configs={"toy": dict(status="fail", chi2=np.float64(12.5),
                                             ndf=np.int64(6), theory=list(range(20))),
                                 "ct18": dict(status="blocked",
                                              reason="environment: LHAPDF set "
                                                     "CT18NLO not installed")})
    """,
    "t2_fake_blocked": """
        def run(verbose=True):
            return dict(name="t2_fake_blocked", generator="not evaluated",
                        reference="none on disk", tolerance=None,
                        status="blocked", reason="needs Doe et al., PRX 1 (2030) 1")
    """,
    "t3_fake_error": """
        def run(verbose=True):
            raise ValueError("boom\\nsecond line never shown")
    """,
    "t3_fake_env": """
        def run(verbose=True):
            raise RuntimeError("Info file not found for PDF set 'CT18NLO'")
    """,
    "t3_fake_sysexit": """
        import sys
        def run(verbose=True):
            sys.exit(2)
    """,
    "t4_fake_import": """
        raise ImportError("no module named frobnicate")
    """,
    "t4_fake_expensive": """
        import os
        EXPENSIVE = True
        SENTINEL = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "ran.flag")
        def run(verbose=True):
            open(SENTINEL, "w").close()
            return dict(name="t4_fake_expensive", generator="ratio 1.003",
                        reference="stock", tolerance="20 %", status="pass")
    """,
    "t6_fake_numpy": """
        import os
        import numpy as np
        HERE = os.path.dirname(os.path.abspath(__file__))
        DATA = os.path.join(HERE, "data", "table.json")
        MISSING = os.path.join(HERE, "data", "not_vendored.csv")
        TABLES = {"a": ("other.json", "not a file name in data/ either")}
        def run(verbose=True):
            print("noise that must not reach the report")
            return dict(name="t6_fake_numpy",
                        generator_value={"t-peak": [np.float64(1.0), np.float32(2.5)]},
                        reference_value=[np.float64(1.25)],
                        tolerance=np.float64(0.2), status=np.str_("fail"),
                        reason=np.str_("recorded"),
                        config=dict(seed=np.int64(3), build_s=16.2,
                                    nested=dict(n=np.int32(4), seconds=1.0)),
                        subrows=[dict(window="W1", observable="n_ch",
                                      status="pass", ratio=np.float64(1.003)),
                                 dict(clause="C13", status="fail",
                                      first_violation="residual 0.82")])
    """,
}


def _write_fakes(bench, names=None):
    bench.mkdir(exist_ok=True)
    (bench / "data").mkdir(exist_ok=True)
    (bench / "data" / "table.json").write_text('{"x": 1}\n')
    for name, body in FAKES.items():
        if names is None or name in names:
            (bench / (name + ".py")).write_text(textwrap.dedent(body))
    return bench


@pytest.fixture()
def fake(tmp_path):
    bench = _write_fakes(tmp_path / "bench")
    plan = tmp_path / "PLAN.md"
    plan.write_text(FAKE_PLAN)
    return bench, plan


def _gen(mr, bench, plan, opt_in=False, **kw):
    return mr.generate(str(bench), opt_in=opt_in, isolate=False,
                       plan_path=str(plan), now=FAKE_NOW, git=FAKE_GIT,
                       env=FAKE_ENV, echo=False, **kw)


def _rows(rep):
    return {r["name"]: r for r in rep["rows"]}


# ------------------------------------------------------------- the paths

def test_every_path_gets_its_status_and_the_generator_survives(mr, fake):
    bench, plan = fake
    rep = _gen(mr, bench, plan)
    rows = _rows(rep)
    assert [r["name"] for r in rep["rows"]] == sorted(FAKES), "sorted discovery"
    assert {n: r["status"] for n, r in rows.items()} == {
        "t1_fake_fail": "fail",
        "t2_fake_blocked": "blocked",
        "t3_fake_env": "blocked (environment)",
        "t3_fake_error": "error",
        "t3_fake_sysexit": "error",
        "t4_fake_expensive": "skipped (opt-in)",
        "t4_fake_import": "error",
        "t5_fake_pass": "pass",
        "t6_fake_numpy": "fail",
    }
    assert rep["summary"] == {"pass": 1, "fail": 2, "blocked": 1,
                              "blocked (environment)": 1, "error": 3,
                              "skipped (opt-in)": 1}
    assert {n: r["tier"] for n, r in rows.items()}["t6_fake_numpy"] == "T6"
    assert rows["t5_fake_pass"]["tolerance"] == "1e-12"
    assert rows["t2_fake_blocked"]["reason"] == "needs Doe et al., PRX 1 (2030) 1"
    assert rows["t2_fake_blocked"]["tolerance"] is None


def test_an_exception_is_error_with_its_first_line_only(mr, fake):
    rows = _rows(_gen(mr, *fake))
    assert rows["t3_fake_error"]["reason"] == "ValueError: boom"
    assert rows["t4_fake_import"]["reason"] == "ImportError: no module named frobnicate"
    assert rows["t3_fake_sysexit"]["reason"] == "SystemExit: 2"
    assert rows["t3_fake_error"]["runtime_s"] is not None
    assert rows["t4_fake_import"]["runtime_s"] is None


def test_a_missing_lhapdf_set_is_blocked_by_the_environment(mr, rp, fake):
    rows = _rows(_gen(mr, *fake))
    assert rows["t3_fake_env"]["reason"] == \
        "environment: LHAPDF set CT18NLO not installed"
    assert rp.classify_exception(RuntimeError(
        "Info file not found for PDF set 'EPPS21nlo_CT18Anlo_Li6'")) == (
        "blocked (environment)",
        "environment: LHAPDF set EPPS21nlo_CT18Anlo_Li6 not installed")
    assert rp.classify_exception(KeyError("x"))[0] == "error"


def test_the_opt_in_row_runs_only_when_asked(mr, fake, monkeypatch):
    bench, plan = fake
    flag = bench / "ran.flag"
    r = _rows(_gen(mr, bench, plan))["t4_fake_expensive"]
    assert r["status"] == "skipped (opt-in)" and not flag.exists()
    assert "LIPOLGEN_BENCH=1" in r["reason"] and r["expensive"] is True
    assert r["runtime_s"] is None
    r = _rows(_gen(mr, bench, plan, opt_in=True))["t4_fake_expensive"]
    assert r["status"] == "pass" and flag.exists()
    # the generator's own list: the concurrent PYTHIA closure row
    assert "t4_pythia_ep_closure" in mr.EXPENSIVE
    # the environment variable is the same switch as --all (in main())
    flag.unlink()
    out = bench.parent / "out"
    monkeypatch.setenv("LIPOLGEN_BENCH", "1")
    assert mr.main(["--bench-dir", str(bench), "--plan", str(plan),
                    "--in-process", "--out-dir", str(out)]) == 0
    assert flag.exists()
    assert json.loads((out / "REPORT.json").read_text())["generated"]["opt_in"]


def test_numpy_scalars_are_coerced_everywhere(mr, rp, fake):
    rep = _gen(mr, *fake)
    r = _rows(rep)["t6_fake_numpy"]
    assert r["generator"] == "t-peak: 1 / 2.5"          # DJANGOH-style dict
    assert r["reference"] == "1.25"
    assert r["tolerance"] == "0.2" and r["status"] == "fail"
    assert r["config"] == {"seed": 3, "nested.n": 4}     # timings dropped
    assert type(r["config"]["seed"]) is int
    text = mr.dump_json(rep) + mr.render_md(rep)
    for token in ("np.float", "np.int", "np.str_", "numpy.", "array("):
        assert token not in text, token
    assert "noise that must not reach the report" not in text
    assert rp.report_line("x", np.float64(0.5), "r", None, np.str_("pass")) == \
        "REPORT | x | 0.5 | r | n/a | PASS"
    buf = io.StringIO()
    line = rp.emit(dict(name="t4_djangoh_rad_noRad", status="blocked",
                        generator_value={"t-peak": [np.float64(1.01616918)]},
                        reference_value=[np.float64(1.07245359)],
                        tolerance="none definable"), stream=buf)
    assert line == ("REPORT | t4_djangoh_rad_noRad | t-peak: 1.01617 | 1.07245 "
                    "| none definable | BLOCKED")
    assert buf.getvalue() == line + "\n"
    assert rp.to_builtin({(1, 2): {np.int64(1)}}) == {"(1, 2)": [1]}
    assert rp.round_sig(float("nan")) == "nan"


def test_subrows_are_extracted_and_rendered(mr, fake):
    rep = _gen(mr, *fake)
    rows = _rows(rep)
    subs = {s["name"]: s for s in rows["t1_fake_fail"]["subrows"]}
    assert list(subs) == ["toy", "ct18"]
    assert subs["toy"]["values"] == {"chi2": 12.5, "ndf": 6}   # long list left out
    assert subs["ct18"]["status"] == "blocked"
    names = [s["name"] for s in rows["t6_fake_numpy"]["subrows"]]
    assert names == ["W1/n_ch", "C13"]
    md = mr.render_md(rep)
    assert "| `toy` | **FAIL** | chi2 = 12.5; ndf = 6 |" in md
    assert ("| `ct18` | **BLOCKED** | **reason:** environment: LHAPDF set "
            "CT18NLO not installed |") in md
    blocked = [b["row"] for b in rep["planned_not_wired"]["blocked"]]
    assert blocked == ["t1_fake_fail[ct18]", "t2_fake_blocked", "t3_fake_env"]
    # an error row's reason stands in its value column
    assert "| **ERROR** | — (*ValueError: boom*) |" in md


def test_vendored_files_named_by_constants_are_hashed(mr, fake):
    bench, plan = fake
    r = _rows(_gen(mr, bench, plan))["t6_fake_numpy"]
    import hashlib
    sha = hashlib.sha256((bench / "data" / "table.json").read_bytes()).hexdigest()
    assert r["data_sha256"] == {"data/not_vendored.csv": None,
                                "data/table.json": sha}


def test_output_is_deterministic(mr, fake):
    a, b = _gen(mr, *fake), _gen(mr, *fake)
    for x in (a, b):                          # runtimes are the one volatile
        for r in x["rows"]:
            r["runtime_s"] = 0.0 if r["runtime_s"] is not None else None
    assert mr.dump_json(a) == mr.dump_json(b)
    assert mr.render_md(a) == mr.render_md(b)
    md = mr.render_md(a)
    assert md.splitlines()[0] == "<!-- SPDX-License-Identifier: GPL-3.0-or-later -->"
    assert "DO NOT EDIT" in md.splitlines()[1]
    assert "- **generated:** %s" % FAKE_NOW in md
    assert "`abc1234` (clean)" in md
    assert a["generated"]["spdx_license_identifier"] == "GPL-3.0-or-later"


# ------------------------------------------------------------- the plan

def test_plan_rows_and_u_rows_are_read_not_invented(mr, fake):
    rep = _gen(mr, *fake)
    rows = _rows(rep)
    assert rows["t5_fake_pass"]["plan_row"] == 1
    assert rows["t1_fake_fail"]["plan_row"] == 3
    assert rows["t6_fake_numpy"]["plan_row"] is None
    pnw = rep["planned_not_wired"]
    assert [p["row"] for p in pnw["plan_rows_without_harness"]] == [2]
    assert pnw["plan_rows_without_harness"][0]["status"] == \
        "**resolved** (not a harness)"
    u = {x["id"]: x for x in pnw["u_rows"]}
    assert list(u) == ["U-1", "U-2"]
    assert u["U-1"] == dict(id="U-1", observable="b1(6Li)", closest_proxy="b1d",
                            distance="A = 2 -> 6")
    # the unescaped pipe in U-2's text: re-joined and said so
    assert u["U-2"]["closest_proxy"] == "a2(\\|t\\|) of the deuteron"
    assert u["U-2"]["distance"] == "3.5x" and "unescaped" in u["U-2"]["parse_note"]


def test_the_real_plan_maps_the_first_ten(mr):
    plan = mr.read_plan(_PLAN)
    assert plan is not None
    assert [r["row"] for r in plan["rows4"]] == list(range(1, 11))
    assert mr.plan_row_map(plan) == {
        "t5_epios_source_modes": 1, "t2_hermes_b1_table2": 3,
        "t3_nmc_li6_over_d": 4, "t3_li6_charge_ff_fb": 5,
        "t3_li_magnetization_rfy": 6, "t4_djangoh_rad_noRad": 7,
        "t5_est_identity": 9, "t2_bonus_spectator_shape": 10}
    assert [r["row"] for r in plan["rows4"] if not r["harnesses"]] == [2, 8]
    assert [u["id"] for u in plan["u_rows"]] == ["U-%d" % i for i in range(1, 13)]


# ------------------------------------------------------------- --compare

def test_compare_prints_every_change_and_counts_it(mr, fake):
    old = _gen(mr, *fake)
    new = json.loads(json.dumps(old))
    buf = io.StringIO()
    assert mr.compare(old, new, out=buf) == 0
    rows = _rows(new)
    rows["t5_fake_pass"]["runtime_s"] = 99.0                 # not a change
    new["generated"]["date"] = "2027-01-01T00:00:00Z"        # not a change
    assert mr.compare(old, new, out=io.StringIO()) == 0
    rows["t5_fake_pass"]["status"] = "fail"
    rows["t1_fake_fail"]["generator"] = "chi2/ndf = 12.6/6"
    rows["t1_fake_fail"]["subrows"][0]["values"]["chi2"] = 12.6
    rows["t1_fake_fail"]["subrows"][1]["status"] = "pass"
    buf = io.StringIO()
    assert mr.compare(old, new, out=buf) == 4
    text = buf.getvalue()
    assert "STATUS   t5_fake_pass: PASS -> FAIL" in text
    assert "MOVED    t1_fake_fail:" in text and "12.6/6" in text
    assert "SUBMOVED t1_fake_fail[toy]: chi2 12.5 -> 12.6" in text
    assert "SUBROW   t1_fake_fail[ct18]: BLOCKED -> PASS" in text
    # a row in the baseline but gone is a change; a new one is not
    newer = json.loads(json.dumps(old))
    newer["rows"] = [r for r in newer["rows"] if r["name"] != "t2_fake_blocked"]
    assert mr.compare(old, newer, out=io.StringIO()) == 1
    assert mr.compare(newer, old, out=io.StringIO()) == 0
    # a row that did not run is ONE change, with its reason, not a cascade
    skipped = json.loads(json.dumps(old))
    r = _rows(skipped)["t1_fake_fail"]
    r.update(status="skipped (opt-in)", generator=None, subrows=[],
             reason=mr.OPT_IN_REASON)
    buf = io.StringIO()
    assert mr.compare(old, skipped, out=buf) == 1
    assert "not run here; rerun with --all" in buf.getvalue()


def test_compare_sees_a_move_of_a_tiny_number(mr, fake):
    """The M2 rule is relative: a 1e-14 headline that doubles, or a 4e-16
    sub-row p-value that becomes 9e-10, is a move (until 2026-09-27 the
    scale was max(1, |x|, |y|) and both passed as 'no change')."""
    old = _gen(mr, *fake)
    rows = _rows(old)
    rows["t1_fake_fail"]["generator"] = "max residual 1.138e-14"
    rows["t1_fake_fail"]["subrows"][0]["values"]["p"] = 4.167883349e-16
    rows["t1_fake_fail"]["subrows"][0]["values"]["zero"] = 0.0
    new = json.loads(json.dumps(old))
    assert mr.compare(old, new, out=io.StringIO()) == 0
    _rows(new)["t1_fake_fail"]["generator"] = "max residual 2.290e-14"
    buf = io.StringIO()
    assert mr.compare(old, new, out=buf) == 1
    assert "MOVED    t1_fake_fail:" in buf.getvalue()
    new = json.loads(json.dumps(old))
    _rows(new)["t1_fake_fail"]["subrows"][0]["values"]["p"] = 9e-10
    buf = io.StringIO()
    assert mr.compare(old, new, out=buf) == 1
    assert "SUBMOVED t1_fake_fail[" in buf.getvalue()
    # the rule itself: relative, zero is not a move, inf / nan changes are
    assert not mr._moved([0.0, 1.0e-300], [0.0, 1.0e-300])
    assert not mr._moved([12.5], [12.5 * (1 + 5e-10)])
    assert mr._moved([12.5], [12.5 * (1 + 2e-9)])
    assert mr._moved([1e-14], [2e-14]) and mr._moved([0.0], [1e-300])
    assert mr._moved([1.0], [float("inf")]) and mr._moved([1.0], [float("nan")])
    assert not mr._moved([float("nan")], [float("nan")])


def test_compare_exit_codes_through_main(mr, fake, tmp_path):
    bench, plan = fake
    common = ["--bench-dir", str(bench), "--plan", str(plan), "--in-process"]
    base = tmp_path / "base"
    assert mr.main(common + ["--out-dir", str(base)]) == 0
    js = base / "REPORT.json"
    assert mr.main(common + ["--compare", str(js)]) == 0      # nothing moved
    old = json.loads(js.read_text())
    _rows(old)["t2_fake_blocked"]["status"] = "pass"
    moved = tmp_path / "moved.json"
    moved.write_text(json.dumps(old))
    assert mr.main(common + ["--compare", str(moved)]) == 1
    assert mr.main(common + ["--compare", str(tmp_path / "absent.json")]) == 2
    # a partial report is never written over the default location
    assert mr.main(common + ["--only", "t5_fake_pass"]) == 2
    assert mr.main(common + ["--only", "t5_fake_pass", "--out-dir",
                             str(tmp_path / "part")]) == 0
    part = json.loads((tmp_path / "part" / "REPORT.json").read_text())
    assert [r["name"] for r in part["rows"]] == ["t5_fake_pass"]
    assert part["generated"]["only"] == ["t5_fake_pass"]
    assert mr.main(common + ["--only", "t9_nope", "--out-dir",
                             str(tmp_path / "x")]) == 2


# ------------------------------------------------------------- --check

def test_check_passes_on_a_fresh_report_and_fails_on_a_hand_edit(mr, fake, tmp_path,
                                                                 monkeypatch):
    bench, plan = fake
    import tempfile
    monkeypatch.setattr(tempfile, "tempdir", str(tmp_path))  # --check's scratch
    common = ["--bench-dir", str(bench), "--plan", str(plan), "--in-process"]
    out = tmp_path / "committed"
    assert mr.main(common + ["--all", "--out-dir", str(out)]) == 0
    assert mr.main(common + ["--check", "--out-dir", str(out)]) == 0
    js, md = out / "REPORT.json", out / "REPORT.md"
    # the volatile header fields and runtimes are ignored ...
    rep = json.loads(js.read_text())
    rep["generated"]["date"] = "1999-01-01T00:00:00Z"
    rep["generated"]["git_rev"] = "0000000"
    for r in rep["rows"]:
        r["runtime_s"] = 123.0 if r["runtime_s"] is not None else None
    js.write_text(mr.dump_json(rep))
    md.write_text(mr.render_md(rep))
    assert mr.main(common + ["--check", "--out-dir", str(out)]) == 0
    # ... a hand edit of the markdown is not
    md.write_text(md.read_text().replace("**PASS**", "**FAIL**", 1))
    assert mr.main(common + ["--check", "--out-dir", str(out)]) == 1
    md.write_text(mr.render_md(rep))
    # ... nor a moved number in the json
    _rows(rep)["t1_fake_fail"]["generator"] = "chi2/ndf = 1.0/6"
    js.write_text(mr.dump_json(rep))
    md.write_text(mr.render_md(rep))
    assert mr.main(common + ["--check", "--out-dir", str(out)]) == 1
    assert mr.main(common + ["--check", "--out-dir", str(tmp_path / "none")]) == 2


# ------------------------------------------------------------- isolation

def test_a_crashing_harness_costs_only_its_own_row(mr, tmp_path):
    bench = _write_fakes(tmp_path / "bench", names={"t5_fake_pass"})
    (bench / "t2_fake_exit.py").write_text(
        "import os\ndef run(verbose=True):\n    os._exit(3)\n")
    (bench / "t3_fake_killed.py").write_text(
        "import os, signal\ndef run(verbose=True):\n"
        "    os.kill(os.getpid(), signal.SIGKILL)\n")
    (bench / "t1_fake_loud.py").write_text(
        "import os\ndef run(verbose=True):\n"
        "    os.write(1, b'C++ writes to fd 1\\n')\n"
        "    return dict(name='t1_fake_loud', generator='g 1', reference='r',"
        " tolerance='t', status='pass')\n")
    rep = mr.generate(str(bench), isolate=True, plan_path=_PLAN, now=FAKE_NOW,
                      git=FAKE_GIT, env=FAKE_ENV, echo=False, timeout=120)
    rows = _rows(rep)
    assert rows["t5_fake_pass"]["status"] == "pass"
    assert rows["t1_fake_loud"]["status"] == "pass"     # fd-1 noise tolerated
    assert rows["t2_fake_exit"]["status"] == "error"
    assert rows["t2_fake_exit"]["reason"] == \
        "harness process exited with status 3 before returning its row"
    assert rows["t3_fake_killed"]["reason"] == \
        "harness process killed by signal 9 (SIGKILL)"


# ------------------------------------------------------------- tiers.py

def test_section7_as_data_matches_the_plan(tiers):
    """The plan is the spec: every sec. 7 bullet, its backticked patterns,
    its non-file items and its tiers (parentheticals aside) as transcribed."""
    import re
    text = open(_PLAN, encoding="utf-8").read()
    sec = text[text.index("## 7."):]
    sec = sec[:sec.index("\n## ", 1)] if "\n## " in sec[1:] else sec
    bullets = [l[2:] for l in sec.splitlines() if l.startswith("- ")]
    assert len(bullets) == len(tiers.SECTION7)
    for line, (patterns, items, tiers_, _note) in zip(bullets, tiers.SECTION7):
        lhs, rhs = line.split("→")
        assert tuple(re.findall(r"`([^`]+)`", lhs)) == patterns, line
        rest = [x.strip() for x in re.sub(r"`[^`]+`", "", lhs).split(",")]
        assert tuple(x for x in rest if x) == items, line
        rhs = re.sub(r"\([^)]*\)", "", rhs)
        assert tuple(re.findall(r"T\d(?:-[A-Za-z]+)?", rhs)) == tiers_, line


def test_a_diff_maps_to_tiers_and_harnesses(tiers):
    t2 = [h for h in tiers.harness_files() if h.startswith("t2_")]
    t4 = [h for h in tiers.harness_files() if h.startswith("t4_")]
    p = tiers.plan_for(["src/core/rc.cpp"])
    assert p["tiers"] == ["T2", "T4-RC"]
    assert p["harnesses"] == sorted(t2 + t4)
    assert p["matches"] == [("src/core/rc.cpp", "rc.*", ("T4-RC", "T2"))]
    p = tiers.plan_for(["include/lipolgen/spin.hpp"])
    assert p["tiers"] == ["T5"] and any("T1 #1" in n for n in p["notes"])
    assert tiers.plan_for(["data/vmc/deuteron/fdeut.av18"])["tiers"] == ["T2", "T3"]
    assert tiers.plan_for(["python/lipolgen/cli.py"])["tiers"] == ["T6"]
    p = tiers.plan_for(["src/core/sf.cpp"])              # sec. 7 names sf.hpp only
    assert p["unmapped"] == ["src/core/sf.cpp"] and p["harnesses"] == []
    p = tiers.plan_for(["docs/USAGE.md", "python/tests/test_cli.py"])
    assert p["tiers"] == [] and p["harnesses"] == [] and p["unmapped"] == []
    p = tiers.plan_for(["validation/benchmarks/t5_est_identity.py"])
    assert p["harness_changes"] == ["t5_est_identity"]
    assert p["harnesses"] == ["t5_est_identity"]
    p = tiers.plan_for(["validation/benchmarks/data/est_identity_sources.json"])
    assert p["data_changes"] == {
        "validation/benchmarks/data/est_identity_sources.json": ["t5_est_identity"]}
    assert any("any PDF-set version" in n for n in p["notes"])
    p = tiers.plan_for(["validation/benchmarks/_report.py",
                        "validation/benchmarks/make_report.py"])
    assert list(p["helper_changes"]) == ["validation/benchmarks/_report.py"]
    assert p["data_changes"] == {} and p["tiers"] == []
    assert any(n.startswith("validation/benchmarks/make_report.py changed")
               for n in p["notes"])
    # _report.py alone maps to no harness, and still asks for the report to
    # be regenerated (make_report.py imports it: it shapes every row)
    p = tiers.plan_for(["validation/benchmarks/_report.py"])
    assert p["helper_changes"] == {"validation/benchmarks/_report.py": []}
    assert p["harnesses"] == [] and p["tiers"] == []
    assert any(n.startswith("validation/benchmarks/_report.py changed")
               and "regenerate the report" in n for n in p["notes"])
    p = tiers.plan_for(["include/lipolgen/coherent.hpp"])
    assert p["tiers"] == ["T4-inclusive", "T4-VM", "T6"]


def test_tiers_cli(tiers, capsys):
    assert tiers.main(["--files", "src/core/rc.cpp", "--names"]) == 0
    names = capsys.readouterr().out.strip().split(",")
    assert names and all(n[:3] in ("t2_", "t4_") for n in names)
    assert tiers.main(["--map"]) == 0
    assert "rc.* -> T4-RC, T2" in capsys.readouterr().out
    assert tiers.main([]) == 2


def test_tiers_since_reads_git_without_writing(tiers):
    if not os.path.isdir(os.path.join(_ROOT, ".git")):
        pytest.skip("not a git checkout: %s has no .git" % _ROOT)
    files = tiers.changed_since("HEAD")
    assert isinstance(files, list) and all(isinstance(f, str) for f in files)
    assert tiers.main(["--since", "no-such-revision-xyz"]) == 2


# ------------------------------------------------------------- the opt-in

def test_the_bench_expensive_marker_is_registered(request):
    lines = request.config.getini("markers")
    assert any(l.startswith("bench_expensive:") and "LIPOLGEN_BENCH=1" in l
               for l in lines), lines


class _Item:
    def __init__(self, marked):
        self.marked, self.added = marked, []

    def get_closest_marker(self, name):
        return object() if (self.marked and name == "bench_expensive") else None

    def add_marker(self, m):
        self.added.append(m)


def test_the_opt_in_hook_skips_unless_the_variable_is_exactly_1(monkeypatch):
    conftest = _load("conftest", os.path.join(os.path.dirname(__file__),
                                              "conftest.py"))
    for value, skipped in ((None, True), ("0", True), ("yes", True), ("1", False)):
        if value is None:
            monkeypatch.delenv("LIPOLGEN_BENCH", raising=False)
        else:
            monkeypatch.setenv("LIPOLGEN_BENCH", value)
        marked, plain = _Item(True), _Item(False)
        conftest.pytest_collection_modifyitems(None, [marked, plain])
        assert plain.added == []
        assert bool(marked.added) is skipped, value
        if skipped:
            mark = marked.added[0].mark
            assert mark.name == "skip"
            assert "LIPOLGEN_BENCH=1" in mark.kwargs["reason"]


# ------------------------------------------------------------- the real report

_ENV_KEYS = ("python", "numpy", "scipy", "pyhepmc", "pythia", "lhapdf",
             "hepmc3", "lhapdf_sets_present")


@pytest.mark.bench_expensive
def test_the_committed_report_is_what_regeneration_gives(mr):
    """`make_report.py --check` on the committed docs/benchmarking/REPORT.*:
    every harness re-run (in the committed report's opt-in mode) and the
    result diffed, ignoring only the date, the revision, the runtimes and
    the environment header.  ~45 s."""
    pytest.importorskip("lipolgen")
    js = os.path.join(_REPORT_DIR, "REPORT.json")
    assert os.path.isfile(js), "docs/benchmarking/REPORT.json is not committed"
    with open(js, encoding="utf-8") as f:
        committed = json.load(f)["generated"]["environment"]
    here = mr.environment()
    differs = ["%s: committed %r, here %r" % (k, committed.get(k), here.get(k))
               for k in _ENV_KEYS if committed.get(k) != here.get(k)]
    if differs:
        pytest.skip("the committed report was made in another environment "
                    "(regenerate it here with make_report.py --all): "
                    + "; ".join(differs))
    p = subprocess.run([sys.executable, os.path.join(_BENCH, "make_report.py"),
                        "--check"], capture_output=True, text=True, timeout=1800)
    assert p.returncode == 0, p.stdout[-6000:] + p.stderr[-3000:]
    assert "check: OK" in p.stdout
