# SPDX-License-Identifier: GPL-3.0-or-later
"""D7 (2026-10-27): a pin of the kernel the DEFAULT run builds.

Since D7 the shipped structure-function defaults are the LHAPDF grid
backends -- CT18NLO (`--unpol-sf ct18nlo`) and NNPDFpol1.1 (`--pol-sf
nnpdfpol`), frozen at their grid edges (LHAPDF `nearest`, set on the object)
and evaluated at leading order.  Every reference JSON of validation/reference/
pins the CLASS-LEVEL TOY kernels (`default_inclusive_kernel(ion)` with one
argument, `ToyF2` / `ToyG1`), which D7 leaves alone -- so no reference pins the
kernel a default run actually builds (the adversarial check of
plan_D7_backend_defaults.md, underrated risk 3).  This file does: the 6Li
inclusive kernel of `lipolgen.shipped_inclusive_kernel(li6())`, asserted to be
the default Pipeline's own kernel bit for bit, at six (x, Q2) points -- one of
them BELOW CT18NLO's grid floor (Q2 = 1.677), which pins the freeze -- and at
the four sweet-spot cells of the letter's money plot.

Measured 2026-10-29 on the I7 build (docs/open_items/run_2026-10-27/
phase_I7_backends.md), rtol 1e-12 like the toy pins.  A move here is a change
of the shipped physics (a grid, its evaluation, the freeze, the Miller b1 or
the 1e-2 Delta on CT18NLO's F1) and must be re-taken with a dated comment.
The sweet-spot a2 is backend-INDEPENDENT to 1e-15: F2 cancels in
-(1/2)(1-y)/y^2 Delta/D_phi at fixed Delta/F1 (plan_D7 sec. 3.2), so the four
values equal the toy kernel's (test_delta_convention.py) -- asserted too.
"""
import json
import math
import os
import subprocess
import sys

import pytest

lg = pytest.importorskip("lipolgen")
_l = lg._lipolgen

pytestmark = pytest.mark.skipif(not _l.HAVE_LHAPDF,
                                reason="no LHAPDF tier: the shipped default cannot be built")

FIELDS = ("f2", "f1", "g1", "b1", "delta")
#: D7 (2026-10-29): (x, Q2) -> f2, f1, g1, b1, delta of the shipped 6Li kernel
#: (per nucleon).  (0.001, 1.2) lies below CT18NLO's Q2 = 1.677: frozen there.
POINTS = {
    (0.001, 1.2): (0.4317738907688184, 183.6114884332517, -0.2506218987136328,
                   0.01987702366216667, 0.2302299420058581),
    (0.01, 2.5): (0.4147549238048744, 17.7029540648422, -0.05673406297858672,
                  0.01986783507271072, 0.04271560096809542),
    (0.05, 5.0): (0.430302869035778, 3.697915280776217, 0.02185486496094642,
                  0.007300617234311208, 0.01226146547757892),
    (0.1, 10.0): (0.397170686049529, 1.726829069780561, 0.03828123801819102,
                  0.004122185209144863, 0.00567831379454835),
    (0.3, 20.0): (0.2166356151253902, 0.3199260138982556, 0.02351256874732484,
                  2.238900853839869e-05, 0.0005352763943447789),
    (0.5, 50.0): (0.06725151184864621, 0.06169863472352863, 0.007063863573879983,
                  -0.0003730035712867323, 3.132178994130887e-05),
}
#: D7 (2026-10-29): the sweet-spot cells at default_configs('6Li')[1]
#: (s = 3980 GeV2), m = +1, theta_S = 90 deg: a2, and the m = 0 rate's w_avg
#: at theta_S = 0 (which carries b1 and so F1's backend)
SWEET = {
    (0.0112, 1.14): (-5.279741070760268e-04, 7.882966883962785e-04),
    (0.0282, 1.14): (-6.499805707128155e-04, 9.705568999934926e-04),
    (0.0708, 3.13): (-7.201153771106359e-04, 1.389263315361509e-03),
    (0.141, 14.3): (-6.632213403199799e-04, 1.180851442079627e-03),
}


def _state(m, theta_s=math.pi / 2):
    return lg.EventSpinState(lam_e=0, pe=0.0, j=1.0, m=m, theta_s=theta_s,
                             phi_s=0.0)


@pytest.fixture(scope="module")
def kernels():
    shipped = lg.shipped_inclusive_kernel(_l.li6())
    p = _l.Pipeline(lg.make_config(events=10), lg.make_plan("tensor-thirds", j=1.0))
    return shipped, p.dis_sampler.kernel


def test_it_is_the_default_pipelines_own_kernel(kernels):
    shipped, pipe = kernels
    for (x, q2) in POINTS:
        a, b = shipped.tables(x, q2), pipe.tables(x, q2)
        for f in FIELDS + ("g2", "b2"):
            assert getattr(a, f) == getattr(b, f), (x, q2, f)


def test_the_shipped_kernel_at_six_points(kernels):
    shipped, _ = kernels
    for (x, q2), want in POINTS.items():
        t = shipped.tables(x, q2)
        for f, w in zip(FIELDS, want):
            assert getattr(t, f) == pytest.approx(w, rel=1e-12), (x, q2, f)


def test_below_the_grid_the_shipped_kernel_is_frozen(kernels):
    shipped, _ = kernels
    floor = float(_l.LhapdfSF("CT18NLO", 0).q2_min)
    assert floor == pytest.approx(1.677025, rel=1e-12)
    for x in (0.001, 0.01, 0.1):
        assert shipped.tables(x, 0.8).f2 == shipped.tables(x, floor).f2


def test_the_sweet_spot_cells(kernels):
    shipped, _ = kernels
    toy = _l.default_inclusive_kernel(_l.li6())
    s = lg.default_configs("6Li")[1].s_per_nucleon()
    for (x, q2), (a2, w0) in SWEET.items():
        t = shipped.tables(x, q2)
        got = shipped.amplitudes(t, x, q2, s, _state(1.0)).a2
        assert got == pytest.approx(a2, rel=1e-12), (x, q2)
        assert shipped.amplitudes(t, x, q2, s, _state(0.0, 0.0)).w_avg == pytest.approx(
            w0, rel=1e-12), (x, q2)
        tt = toy.tables(x, q2)                   # F2 cancels in a2 cell by cell
        assert got == pytest.approx(toy.amplitudes(tt, x, q2, s, _state(1.0)).a2,
                                    rel=1e-14), (x, q2)


# ------------------------------------------------- the default vs lhapdf.conf
#
# D7 amendment 1 (2026-10-29).  CT18NLO.info and NNPDFpol11_100.info pin no
# Extrapolator, Interpolator or ForcePositive, so before D7 a `ct18nlo` number
# was whatever the FIRST lhapdf.conf LHAPDF found said (measured by the
# adversarial check: `nearest` moved the 6Li rate +5.03 %, `ForcePositive: 1`
# the 7Li A_par x1.488).  The tier now sets all three on every grid it opens
# and clamps at the knot edges, so the default numbers below are computed in
# two FRESH processes -- one on this environment's store, one on the same set
# files under an lhapdf.conf that says the opposite of the tier on every key --
# and must agree to the bit.  `sigma_pb` and the per-category sigma are cell
# integrals, independent of the event count.

_ENV_PROBE = r"""
import json, lipolgen as lg
_l = lg._lipolgen
out = {}
for iso, ch, plan in (("6Li", "inclusive", "tensor-thirds"),
                      ("6Li", "coherent", "tensor-thirds"),
                      ("d", "tagged-d-p", "tensor-thirds"),
                      ("7Li", "inclusive", "helicity-flip")):
    cfg = lg.make_config(isotope=iso, channel=ch, events=10, seed=20261026)
    p = _l.Pipeline(cfg, lg.make_plan(plan, j=lg.ion_spin(cfg.isotope)))
    out[iso + "/" + ch] = [repr(v) for v in [p.sigma_pb()] + list(p.sigma_per_category_pb())]
    m = p.generate(0, False, 1)["meta"]
    out[iso + "/" + ch].append(m["unpol_sf_extrapolation"])
print("PROBE " + json.dumps(out, sort_keys=True))
"""

_OPPOSITE_CONF = """Verbosity: 0
Interpolator: linear
Extrapolator: continuation
ForcePositive: 1
AlphaS_Type: analytic
"""


def _store():
    for d in os.environ.get("LHAPDF_DATA_PATH", "").split(":"):
        if d and os.path.isfile(os.path.join(d, "CT18NLO", "CT18NLO.info")) and \
                os.path.isfile(os.path.join(d, "NNPDFpol11_100", "NNPDFpol11_100.info")):
            return d
    return None


def _probe(data_path):
    env = dict(os.environ)
    if data_path is not None:
        env["LHAPDF_DATA_PATH"] = data_path
    r = subprocess.run([sys.executable, "-c", _ENV_PROBE], capture_output=True,
                       text=True, env=env, timeout=300)
    assert r.returncode == 0, r.stderr[-2000:]
    line = [l for l in r.stdout.splitlines() if l.startswith("PROBE ")]
    assert len(line) == 1, r.stdout[-2000:]
    return json.loads(line[0][len("PROBE "):])


def test_the_default_numbers_do_not_depend_on_lhapdf_conf(tmp_path):
    store = _store()
    if store is None:
        pytest.skip("LHAPDF_DATA_PATH names no directory holding both shipped "
                    "sets, so the opposite environment cannot be assembled")
    opp = tmp_path / "opposite"
    opp.mkdir()
    (opp / "lhapdf.conf").write_text(_OPPOSITE_CONF)
    for s in ("CT18NLO", "NNPDFpol11_100"):
        os.symlink(os.path.join(store, s), str(opp / s))
    here = _probe(None)
    there = _probe(str(opp) + "::")          # `::` keeps the compiled-in prefix out
    assert here == there
    # ... and they are the frozen numbers, not merely two equal ones
    assert float(here["6Li/inclusive"][0]) == pytest.approx(496334.155954048, rel=1e-12)
    assert all(v[-1] == _l.LHAPDF_EXTRAPOLATOR == "nearest" for v in here.values())


_PREFLIGHT_PROBE = """
import lipolgen._lipolgen as _l
print("PROBE " + repr(_l.lhapdf_set_problem("CT18NLO")))
"""


def test_the_preflight_advice_follows_what_is_missing(tmp_path):
    """F7 (2026-10-30), the verifier's V7-10: with the sets present and ONLY
    LHAPDF's lhapdf.conf missing, the pre-flight said "Install the set"; it now
    says to put lhapdf.conf on the data path, and still names the escape.  A
    missing SET keeps the install advice and LHAPDF's own "Info file not
    found" (make_report's `blocked (environment)` contract)."""
    import ast
    store = _store()
    if store is None:
        pytest.skip("LHAPDF_DATA_PATH names no directory holding both shipped "
                    "sets, so the environment cannot be assembled")
    noconf = tmp_path / "noconf"
    noconf.mkdir()
    for s in ("CT18NLO", "NNPDFpol11_100"):
        os.symlink(os.path.join(store, s), str(noconf / s))
    empty = tmp_path / "empty"
    empty.mkdir()
    got = {}
    for name, d in (("noconf", noconf), ("empty", empty)):
        env = dict(os.environ, LHAPDF_DATA_PATH=str(d) + "::")
        r = subprocess.run([sys.executable, "-c", _PREFLIGHT_PROBE], capture_output=True,
                           text=True, env=env, timeout=120)
        assert r.returncode == 0, r.stderr[-2000:]          # never the abort, 134
        line = [l for l in r.stdout.splitlines() if l.startswith("PROBE ")]
        got[name] = ast.literal_eval(line[0][len("PROBE "):])   # a repr'd str
    assert "lhapdf.conf was not found" in got["noconf"]
    assert "Put LHAPDF's own lhapdf.conf" in got["noconf"]
    assert "Install the set" not in got["noconf"]
    assert "--unpol-sf toy" in got["noconf"]
    assert "Info file not found for PDF set 'CT18NLO'" in got["empty"]
    assert "Install the set" in got["empty"] and "--unpol-sf toy" in got["empty"]
