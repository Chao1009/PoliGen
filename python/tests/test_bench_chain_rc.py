# SPDX-License-Identifier: GPL-3.0-or-later
"""Benchmark harness group CHAIN AND RADIATIVE (run 2026-09-23, phase B3).

T4 row 7 of docs/benchmarking/BENCHMARK_PLAN.md: DJANGOH's Rad/noRad
four-bin table (`validation/benchmarks/t4_djangoh_rad_noRad.py`, reference
vendored in `validation/benchmarks/data/djangoh_rad_norad_ep18x275.json`).

The row is BLOCKED, and these tests pin WHY as well as WHAT: the ePIC DJANGOH
samples ran with the elastic radiative tail switched off (IEL2 = IEL31 =
IEL32 = IEL33 = 0), and the elastic tail is the only O(alpha) term LiPolGen's
`rc_tail` computes -- so the two share no term and no tolerance exists.  The
samples' hadronic W_h >= 3 GeV cut excludes that tail as well (W_h = M_p in
ep -> e p gamma) whatever the IEL flags are, so the row unblocks only with a
DJANGOH run at IEL31..33 != 0 AND WMIN <= M_p.  The
measured numbers are pinned so that a move in the proton tail (or in the
vendored table) fails a test instead of silently rotting the record
docs/open_items/run_2026-09-23/phase_B3_chain_rc.md.  The vendored table's
sha256 is re-checked on every read (one flipped byte is refused).

The full row (both tail models, 17.3-19.6 s measured 2026-09-23) runs always: it is under the 30 s
threshold of BENCHMARK_PLAN.md sec. 3 / the run's ground rules.
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


def test_djangoh_row_is_blocked_and_prints_both_sides(dj, dj_row, capsys):
    row = dj_row
    assert row["name"] == "t4_djangoh_rad_noRad"
    assert row["status"] == "blocked"
    assert "IEL2=IEL31=IEL32=IEL33=0" in row["reason"]
    # the unblock recipe names the W_h cut, not the IEL flags alone
    assert "W_h >= 3 GeV cut also excludes the elastic tail" in row["reason"]
    assert "IEL31..33 != 0 AND WMIN <= M_p" in row["reason"]
    assert "none definable" in row["tolerance"]
    assert set(row["generator_value"]) == {"t-peak", "polrad-full"}
    for tm, vals in row["generator_value"].items():
        assert len(vals) == 4
        # the one thing that CAN fail: a tail is non-negative and finite
        assert all(math.isfinite(v) and v >= 1.0 for v in vals), tm
        # plain floats, so the ROW line prints numbers, not np.float64(...)
        assert all(type(v) is float for v in vals), tm
    # A tail can never exceed... nothing: the two share no term.  What IS a
    # measured fact (2026-09-23) is that on this window the proton elastic
    # tail is 9-37 % of DJANGOH's inelastic excess -- never more.
    for d in row["detail"]:
        for tm in ("t-peak", "polrad-full"):
            assert 0.0 < d[tm + ":excess_over_djangoh_excess"] < 0.4
    dj.print_row(row, dj.TAIL_MODELS)
    out = capsys.readouterr().out
    assert "ROW | t4_djangoh_rad_noRad" in out and "BLOCKED" in out
    assert "np.float64" not in out
    assert "'t-peak': [1.01616918, 1.01317215, 1.01497026, 1.0206861]" in out


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
    """The numbers of phase_B3_chain_rc.md sec. 1, at the harness's own
    quadrature (converged to <= 1.5e-5 absolute; see `convergence`)."""
    got = dj_row["generator_value"]
    want = {
        "t-peak": (1.01616918, 1.01317215, 1.01497026, 1.02068610),
        "polrad-full": (1.02709574, 1.01507182, 1.01513763, 1.02075698),
    }
    for tm, vals in want.items():
        for g, w in zip(got[tm], vals):
            assert g == pytest.approx(w, rel=2e-7), tm
    assert dj_row["reference_value"] == pytest.approx(
        [1.07245359, 1.08940691, 1.16062561, 1.2279651], rel=1e-8)
    # The Q^2 TREND of (ratio - 1): DJANGOH rises, the proton tail is flat to
    # within +-0.04 in d ln / d ln Q^2 -- different pieces, different trends.
    tr = dj_row["trend_dlnexcess_dlnQ2"]
    assert tr["djangoh"] == pytest.approx(0.1748, abs=5e-4)
    assert abs(tr["t-peak"]) < 0.05 and abs(tr["polrad-full"]) < 0.05
    # The window the numbers belong to.
    cfg = dj_row["config"]
    assert cfg["y_max"] == pytest.approx(0.985)
    assert cfg["x_min"] == 1e-5 and cfg["y_min"] == 1e-4
    assert cfg["W2_min"] == 9.0
