# SPDX-License-Identifier: GPL-3.0-or-later
"""Benchmark harness group CHAIN AND RADIATIVE (run 2026-09-23, phase B3;
re-wired by run 2026-10-26, phase B3).

T4 row 7 of docs/benchmarking/BENCHMARK_PLAN.md: DJANGOH and LiPolGen's
radiative tail (`validation/benchmarks/t4_djangoh_rad_noRad.py`).  Since
2026-10-26 the row's PRIMARY comparison is DJANGOH's isolated elastic
radiative tail (`data/djangoh_eltail_ep18x275.json`), pinned in
test_bench_djangoh_eltail.py.  THIS file keeps the guards of the ePIC
Rad/noRad four-bin table (`data/djangoh_rad_norad_ep18x275.json`), which is
now the row's RECORDED HISTORY: the ePIC samples ran with the elastic
radiative tail switched off (IEL2 = IEL31 = IEL32 = IEL33 = 0) and with a
hadronic W_h >= 3 GeV cut that excludes the tail whatever IEL is, so the
table shares no O(alpha) term with LiPolGen's `rc_tail` and is reported as
recorded sub-rows, never as the verdict.  The 2026-09-23 numbers
(docs/open_items/run_2026-09-23/phase_B3_chain_rc.md sec. 1) are pinned on
the history rows so that a move in the proton tail (or in the vendored
table) fails a test instead of silently rotting the record.  The vendored
table's sha256 is re-checked on every read (one flipped byte is refused).

The full row (both tail models, 22-26 s measured 2026-10-26) runs always: it
is under the 30 s threshold of BENCHMARK_PLAN.md sec. 3.
"""

import importlib.util
import json
import math
import os
import subprocess
import sys

import pytest

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(
    __file__))))
_BENCH = os.path.join(_ROOT, "validation", "benchmarks")


def _load(name):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(_BENCH, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def dj():
    return _load("t4_djangoh_rad_noRad")


@pytest.fixture(scope="module")
def dj_row(dj):
    return dj.run(verbose=False)


def test_djangoh_reference_is_vendored_with_provenance(dj, tmp_path):
    # a byte-identical copy loads; ONE flipped byte is refused (sha256), so
    # the harness exits 2 instead of printing moved numbers
    raw = open(dj.DATA, "rb").read()
    assert raw.count(b"0.6733778912") == 1
    good = tmp_path / "good.json"
    good.write_bytes(raw)
    dj.load_reference(good)
    bad = tmp_path / "bad.json"
    bad.write_bytes(raw.replace(b"0.6733778912", b"0.6733778913"))
    with pytest.raises(RuntimeError, match="sha256"):
        dj.load_reference(bad)
    doc, rows = dj.load_reference()
    prov = doc["provenance"]
    for key in ("source_table", "date_fetched", "licence", "papers",
                "readme_sha256", "repository_commit"):
        assert prov[key], key
    assert "9869d9a0d21b9191a1bbec7294b4851feb5ee20a" in prov["source_table"]
    assert prov["date_fetched"] == "2026-09-23"
    # The four bins and their ratios, computed from the two columns.
    assert [(r["q2_lo"], r["q2_hi"]) for r in rows] == [
        (1.0, 10.0), (10.0, 100.0), (100.0, 1000.0), (1000.0, 10000.0)]
    for r, want in zip(rows, (1.0724536, 1.0894069, 1.1606256, 1.2279651)):
        assert r["ratio"] == pytest.approx(want, rel=1e-7)
    # The seven logs that carry the HERACLES cross section agree with the
    # README column to every digit they print (5 significant figures); the
    # Rad Q2 1-10 log is truncated and carries none.
    by_file = {lg["file"]: lg for lg in doc["log_cross_check"]["logs"]}
    tags = {"Rad": 2, "norad": 3}
    n = 0
    for f, lg in by_file.items():
        kind = "Rad" if ".Rad." in f else "norad"
        lo = f.split("_Q2_")[1].split("_")[0]
        row = [r for r in doc["rows"] if int(r[0]) == int(lo)][0]
        if lg["sigma_nb"] is None:
            assert kind == "Rad" and lo == "1" and "TRUNCATED" in lg["note"]
            continue
        readme_nb = row[tags[kind]] * 1e3
        assert float("%.5g" % readme_nb) == pytest.approx(lg["sigma_nb"],
                                                          rel=1e-12), f
        n += 1
    assert n == 7
    # The fact the row turns on is IN the vendored file, not only in prose --
    # and so is the hadronic cut that excludes the elastic tail whatever IEL is.
    assert "IEL2 = IEL31 = IEL32 = IEL33 = 0" in \
        doc["run_settings"]["IMPORTANT_elastic_tail_off"]
    assert "W_h >= 3.0 GeV" in doc["run_settings"]["kinem_cuts"]
    assert "hadronic final-state mass" in doc["run_settings"]["kinem_cuts"]


def test_djangoh_epic_table_is_recorded_history_beside_the_tail(dj, dj_row,
                                                                capsys):
    """The ePIC Rad/noRad table is no longer the verdict: it is four
    RECORDED history sub-rows, with the reason it cannot be a comparison."""
    row = dj_row
    assert row["name"] == "t4_djangoh_rad_noRad"
    hist = [s for s in row["subrows"] if s["name"].startswith("history_epic[")]
    assert len(hist) == 4
    for s in hist:
        assert s["status"] == "recorded" and s["in_verdict"] is False
        assert "elastic tail OFF and W_h >= 3 GeV" in s["reason"]
        assert "disjoint pieces, no tolerance" in s["reason"]
    h = row["history"]
    assert set(h["lipolgen_ratio"]) == {"t-peak", "polrad-full"}
    for tm, vals in h["lipolgen_ratio"].items():
        assert len(vals) == 4
        # the one thing that CAN fail on the history: a tail is >= 0, finite
        assert all(math.isfinite(v) and v >= 1.0 for v in vals), tm
        assert all(type(v) is float for v in vals), tm
        # measured 2026-09-23 and again 2026-10-26: on this window the proton
        # elastic tail is 9-37 % of DJANGOH's inelastic excess -- never more
        for v, e in zip(vals, h["epic_ratio"]):
            assert 0.0 < (v - 1.0) / (e - 1.0) < 0.4
    dj.print_row(row)
    out = capsys.readouterr().out
    assert "REPORT | t4_djangoh_rad_noRad |" in out
    assert "np.float64" not in out
    assert out.count("[history_epic[") == 4


def test_djangoh_harness_exits_2_without_numpy(tmp_path):
    """A missing dependency is a broken harness (exit 2), numpy included --
    its import is the first thing the module does."""
    fake = tmp_path / "numpy"
    fake.mkdir()
    (fake / "__init__.py").write_text(
        'raise ImportError("simulated: numpy not installed")\n')
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(tmp_path)] + [p for p in [env.get("PYTHONPATH")] if p])
    p = subprocess.run([sys.executable, os.path.join(
        _BENCH, "t4_djangoh_rad_noRad.py")], env=env, capture_output=True,
        text=True)
    assert p.returncode == 2, p.stderr
    assert "simulated: numpy not installed" in p.stderr


def test_djangoh_row_numbers_are_the_recorded_ones(dj_row):
    """The numbers of phase_B3_chain_rc.md sec. 1, now the row's history
    (the W_l^2 >= 9 window, the harness's own quadrature, converged to
    <= 1.5e-5 absolute on 2026-09-23)."""
    h = dj_row["history"]
    want = {
        "t-peak": (1.01616918, 1.01317215, 1.01497026, 1.02068610),
        "polrad-full": (1.02709574, 1.01507182, 1.01513763, 1.02075698),
    }
    for tm, vals in want.items():
        for g, w in zip(h["lipolgen_ratio"][tm], vals):
            assert g == pytest.approx(w, rel=2e-7), tm
    assert h["epic_ratio"] == pytest.approx(
        [1.07245359, 1.08940691, 1.16062561, 1.2279651], rel=1e-8)
    # The Q^2 TREND of (ratio - 1): DJANGOH's inelastic excess rises, the
    # proton tail is flat to within +-0.04 in d ln / d ln Q^2.
    tr = h["trend_dlnexcess_dlnQ2"]
    assert tr["djangoh"] == pytest.approx(0.1748, abs=5e-4)
    assert abs(tr["t-peak"]) < 0.05 and abs(tr["polrad-full"]) < 0.05
    # The window the numbers belong to.
    cfg = dj_row["config"]
    assert cfg["y_max"] == pytest.approx(0.985)
    assert cfg["x_min"] == 1e-5 and cfg["y_min"] == 1e-4
    assert cfg["W2_min"] == 9.0
