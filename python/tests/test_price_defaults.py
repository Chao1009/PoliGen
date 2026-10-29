# SPDX-License-Identifier: GPL-3.0-or-later
"""validation/price_defaults.py -- the default-input price table of the
ranked observables (run 2026-10-26, PLAN.md phase A3) -- runs, returns its
table and draws its figure.

ONE fast test: N = 2000 and the inclusive-6Li family only (about 5 s on the
2026-10-26 machine).  The full table -- N = 200 000, six families, about
3 min -- is the record's: docs/open_items/run_2026-10-26/phase_A_baseline.md,
section "A3 -- the price table".  The script is a RECIPE, not a harness: it
has no tolerance and no verdict.  This test checks that it runs (exit 0) and
that the table it returns carries what the recipe promises -- the shipped
default as numbers, every other READ setting as a ratio, a NOT-READ knob as
its label with the run bit-identical to the default, a REFUSED knob as its
refusal text, a backend whose tier or LHAPDF set is missing as `environment`
-- and that `--plot` writes a non-empty PNG.

The default row's numbers are N-INDEPENDENT (the run's own cross sections and
amplitude tables), so they are pinned at the values the script measured on
2026-10-26 at N = 200 000; the two backend ratios are README.md's selector-
table costs (x0.7985 / x0.7934 on sigma, x1.2524 / x1.2604 on A_zz), which
the script reproduces.  Nothing here asserts that a default should change.
Not network-dependent.
"""
import importlib.util
import json
import os
import subprocess
import sys

import pytest

lg = pytest.importorskip("lipolgen")

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRIPT = os.path.join(_ROOT, "validation", "price_defaults.py")

# measured by the script, 2026-10-26: 6Li config 1, the shipped window
SIGMA_TOY_PB = 591846.1688784091          # tensor-thirds, meta["sigma_pb"]
AZZ_TOY = -5.193230890419342e-04          # thirds estimator on meta's per-category sigma
# cos 2phi' amplitude per unit P_zz, super-bins S1..S4 (sibling fig1.json).
# Re-taken at D4 (2026-10-27, the Jaffe-Manohar coefficient): exactly 1/2 of
# the 2026-10-26 -1.302913804135417e-03 / -1.0625622769273495e-03 /
# -1.4290788342523734e-03 / -1.3290877454563259e-03 (a sum of halved terms).
A2_TOY = (-6.514569020677085e-04,
          -5.312811384636748e-04,
          -7.145394171261867e-04,
          -6.645438727281629e-04)
#: backend -> (sigma ratio, A_zz ratio); README.md quotes x0.7934 / x1.2604
#: (mstw) and x0.7985 / x1.2524 (ct18nlo)
RATIOS = {"mstw": (0.7933758399382593, 1.2604366677843946),
          "ct18nlo": (0.7984708573984506, 1.2523938610005776)}


def test_price_defaults_li6_inclusive_runs_returns_table_and_png(tmp_path):
    have_mpl = importlib.util.find_spec("matplotlib") is not None
    table_json = tmp_path / "price_defaults.json"
    cmd = [sys.executable, SCRIPT, "--events", "2000", "--families", "li6-incl",
           "--json", str(table_json)]
    if have_mpl:
        cmd += ["--plot", "--out-dir", str(tmp_path)]
    env = dict(os.environ)
    pkg_parent = os.path.dirname(os.path.dirname(os.path.abspath(lg.__file__)))
    env["PYTHONPATH"] = os.pathsep.join(
        [pkg_parent] + [p for p in env.get("PYTHONPATH", "").split(os.pathsep) if p])
    r = subprocess.run(cmd, capture_output=True, text=True, env=env,
                       cwd=str(tmp_path), timeout=600)
    assert r.returncode == 0, r.stdout[-2000:] + r.stderr[-4000:]
    assert "| setting | provenance |" in r.stdout

    table = json.loads(table_json.read_text())
    assert table["events"] == 2000 and table["seed"] == 20261026
    fam = table["families"]["li6-incl"]
    rows = {row["setting"]: row for row in fam["rows"]}

    # the shipped default: numbers, no ratio
    d = fam["rows"][0]
    assert d["provenance"] == "default" and not d["ratios"]
    assert d["values"]["sigma"] == pytest.approx(SIGMA_TOY_PB, rel=1e-9)
    assert d["values"]["azz"] == pytest.approx(AZZ_TOY, rel=1e-9)
    for i, a2 in enumerate(A2_TOY, 1):
        assert d["values"]["A2_S%d" % i] == pytest.approx(a2, rel=1e-9)

    # the READ backends: ratios -- or `environment`, never a number
    for backend, (r_sigma, r_azz) in RATIOS.items():
        row = rows["unpol_sf = %s" % backend]
        if row["provenance"] == "environment":
            assert row["label"] and not row["values"] and not row["ratios"]
            continue
        assert row["provenance"] == "read"
        assert row["ratios"]["sigma"] == pytest.approx(r_sigma, rel=1e-6)
        assert row["ratios"]["azz"] == pytest.approx(r_azz, rel=1e-6)
        # F2 cancels in the cos 2phi' amplitude cell by cell
        assert row["values"]["cell_max_rel_a2"] < 1e-12
        for i in range(1, 5):
            assert abs(row["ratios"]["A2_S%d" % i] - 1.0) < 5e-3

    # NOT READ under a tensor fill: the label, and the run bit-identical
    for plan in ("transverse-tensor", "tensor-thirds"):
        row = rows["pol_sf = nnpdfpol (%s)" % plan]
        if row["provenance"] == "environment":
            continue
        assert row["provenance"] == "not-read"
        assert row["label"] == "not read under this run's fill"
        assert row["identical"] is True
        assert not row["values"] and not row["ratios"]

    # REFUSED on the shipped miller b1: the refusal text, never a number
    row = rows["r_source = r1998"]
    assert row["provenance"] == "refused"
    assert "read ONLY by b1_model = li6-convolution" in row["label"]
    assert not row["values"] and not row["ratios"]

    if not have_mpl:
        pytest.skip("matplotlib not installed: the table was checked, the "
                    "price_defaults.png figure was not drawn")
    png = tmp_path / "price_defaults.png"
    assert png.is_file() and png.stat().st_size > 1000
    assert png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
