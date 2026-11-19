# SPDX-License-Identifier: GPL-3.0-or-later
"""K1 (2026-11-17): the digitized deuteron b1 tables read as the ratio b1/F1 at
a reference scale on the run surface -- `--b1-ref-q2` / `PipelineConfig.
b1_reference_q2` (default B1_REFERENCE_Q2 = 2.5 GeV2; 0 = the frozen table).

docs/open_items/run_2026-11-17/DECISIONS.md K1 is the AUTHOR's decision; its
details (a)-(g) are the supervisor's choices under it.  Record:
docs/open_items/run_2026-11-17/phase_P1_b1_ratio.md.  Tolerances declared
before the values were seen: the b1/F1 identity at 1e-12 relative; the frozen
reading against a caller-supplied frozen kernel bit for bit.
"""
import subprocess
import sys

import numpy as np
import pytest

lg = pytest.importorskip("lipolgen")
_l = lg._lipolgen

pytestmark = pytest.mark.skipif(not _l.HAVE_LHAPDF,
                                reason="needs the LHAPDF tier: the shipped run surface")
N = 200
TT = lg.tensor_thirds_plan(0.7, 0.6)


def _cfg(x_max=None, **kw):
    cfg = lg.make_config(seed=3, **dict(dict(events=N), **kw))
    if x_max is not None:                    # the Scenario comes back by value
        sc = cfg.scenario
        sc.x_max = x_max
        cfg.scenario = sc
    return cfg


def _sig(p):
    return np.asarray(p.sigma_per_category_pb(), dtype=float).tobytes()


def _run(plan=TT, x_max=None, **kw):
    cfg = _cfg(x_max, **kw)
    p = _l.Pipeline(cfg, plan)
    return p, p.generate(0, False, 1)


def _rows(p):
    return {r.name: r for r in p.knob_provenance()}


def test_the_constant_and_the_default():
    assert _l.B1_REFERENCE_Q2 == 2.5
    assert lg.make_config().b1_reference_q2 == 2.5
    assert _l.InclusiveKernel.Options().b1_reference_q2 == 0.0          # class level: frozen
    assert lg.shipped_inclusive_kernel(_l.li6()).b1_reference_q2 == 0.0  # the factory: frozen unless asked


@pytest.mark.parametrize("isotope", ["6Li", "d"])
def test_the_default_run_reads_it_and_says_so(isotope):
    p, cols = _run(isotope=isotope)
    k = p.dis_sampler.kernel
    assert k.b1_reference_q2 == 2.5
    m = cols["meta"]
    assert m["b1_reference_q2"] == 2.5
    assert "K1 (2026-11-17): the digitized table is read as the ratio b1/F1 at Q0^2 = 2.5 GeV2" in m["rank2_input"]
    r = _rows(p)["b1_reference_q2"]
    assert r.status == _l.KnobStatus.Read and r.flag == "--b1-ref-q2" and r.value == "2.5"
    # the identity, on the run's own kernel: b1/F1 the same at every Q2 at fixed x
    for x in (0.1, 0.3, 0.5, 0.7):
        ref = k.tables(x, 2.5)
        for q2 in (10.0, 100.0, 1897.0):
            t = k.tables(x, q2)
            if isotope == "d" and k.unconfined_b1(x, q2, t.f1) != t.b1:
                continue                                   # the confinement binds (not at these x)
            assert t.b1 / t.f1 == pytest.approx(ref.b1 / ref.f1, rel=1e-12), (x, q2)


@pytest.mark.parametrize("isotope", ["6Li", "d"])
def test_zero_is_the_frozen_table_bit_for_bit(isotope):
    """`--b1-ref-q2 0` reproduces the frozen reading: the same arrays as a run on
    the frozen shipped kernel handed over as a caller-supplied kernel."""
    p0, c0 = _run(isotope=isotope, b1_reference_q2=0.0)
    assert p0.dis_sampler.kernel.b1_reference_q2 == 0.0
    assert c0["meta"]["b1_reference_q2"] == 0.0
    assert "--b1-ref-q2 0, the frozen table" in c0["meta"]["rank2_input"]
    cfg = lg.make_config(events=N, seed=3, isotope=isotope)
    cfg.kernel = lg.shipped_inclusive_kernel(lg.ion_by_name(cfg.isotope))
    pk = _l.Pipeline(cfg, TT)
    ck = pk.generate(0, False, 1)
    arr = [k for k, v in c0.items() if isinstance(v, np.ndarray)]
    assert arr
    for k in arr:
        assert c0[k].tobytes() == ck[k].tobytes(), k
    assert _sig(p0) == _sig(pk)
    # ... and the default (2.5) moves the tensor weights
    p1, _ = _run(isotope=isotope)
    assert _sig(p1) != _sig(p0)


@pytest.mark.parametrize("kw,match", [
    (dict(b1_reference_q2=-1.0), "must be finite and >= 0"),
    (dict(b1_reference_q2=float("nan")), "must be finite and >= 0"),
    (dict(b1_reference_q2=float("inf")), "must be finite and >= 0"),
    (dict(b1_reference_q2=1.0), "outside the unpolarised grid's own Q2 range"),
    (dict(b1_reference_q2=1.6), "outside the unpolarised grid's own Q2 range"),
    (dict(isotope="7Li", b1_reference_q2=0.0), r"not read on 7Li \(no rank-2 input\)"),
    (dict(channel="coherent", b1_reference_q2=0.0), "not read on channel coherent-6Li"),
    (dict(channel="tagged-6Li-alpha", b1_reference_q2=10.0),
     "not read on channel tagged-6Li-alpha without --inclusive-b1"),
    (dict(channel="tagged-d-p", b1_reference_q2=10.0), "not read on channel tagged-d-p"),
    (dict(b1_model="li6-convolution", b1_unpol="ct18nlo", b1_reference_q2=0.0, x_max=0.95),
     "not read at b1_model = li6-convolution"),
])
def test_every_refusal(kw, match):
    plan = lg.make_plan("helicity-flip", j=1.5) if kw.get("isotope") == "7Li" else TT
    with pytest.raises((RuntimeError, ValueError), match=match):
        _l.Pipeline(_cfg(**dict(kw, events=10)), plan)


def test_a_value_inside_the_grid_is_accepted_and_read():
    floor = _l.make_lhapdf_unpol_sf(_l.UNPOL_SF_DEFAULT).q2_min
    p, c = _run(b1_reference_q2=floor)
    assert p.dis_sampler.kernel.b1_reference_q2 == floor
    assert c["meta"]["b1_reference_q2"] == floor
    cfg = lg.make_config(events=10, b1_reference_q2=1.0, unpol_sf="toy", pol_sf="toy")
    _l.Pipeline(cfg, TT)                                  # ToyF2 reports no grid: accepted


@pytest.mark.parametrize("kw,label", [
    (dict(isotope="7Li"), "not read on 7Li (no rank-2 input)"),
    (dict(channel="coherent"), "not read on channel coherent-6Li"),
    (dict(channel="tagged-6Li-alpha"), "not read on channel tagged-6Li-alpha without --inclusive-b1"),
    (dict(b1_model="li6-convolution", b1_unpol="ct18nlo", x_max=0.95), "not read at b1_model = li6-convolution"),
])
def test_the_not_read_label_where_it_is_not_read(kw, label):
    plan = lg.make_plan("helicity-flip", j=1.5) if kw.get("isotope") == "7Li" else TT
    for explicit in ({}, dict(b1_reference_q2=2.5)):           # A17: the explicit default too
        p = _l.Pipeline(_cfg(**dict(kw, events=10, **explicit)), plan)
        r = _rows(p)["b1_reference_q2"]
        assert r.status == _l.KnobStatus.NotRead and r.label == label, r.reason[:200]
        assert p.generate(0, False, 1)["meta"]["b1_reference_q2"] == label
        assert p.dis_sampler.kernel.b1_reference_q2 == 0.0


def test_the_li6_convolution_path_is_untouched():
    """Not a digitized table: the knob is not read there and the kernel is built
    with the option off; the factory refuses a positive value with it."""
    p, c = _run(b1_model="li6-convolution", b1_unpol="ct18nlo", x_max=0.95)   # refused at 1.0 (the rank-2 check, as before K1)
    assert p.dis_sampler.kernel.b1_reference_q2 == 0.0
    assert "K1" not in c["meta"]["rank2_input"]
    with pytest.raises(RuntimeError, match="li6-convolution is refused"):
        _l.default_inclusive_kernel(_l.li6(), model=_l.B1Model.Li6Convolution,
                                    b1_reference_q2=2.5)


def test_tagged_6li_alpha_with_inclusive_b1_reads_it():
    """The struck deuteron's MillerB1 slot.  Measured 2026-11-19: the frozen
    table is refused at the default x_max 1.0 (the sampler's rank-2 check at
    x = 0.955, Q2 = 167.3, the struck deuteron being unconfined), as it was
    before K1; under K1 the run builds.  Compared at x_max 0.85."""
    p, c = _run(channel="tagged-6Li-alpha", inclusive_b1=True)
    assert p.dis_sampler.kernel.b1_reference_q2 == 2.5
    assert c["meta"]["b1_reference_q2"] == 2.5
    with pytest.raises(RuntimeError, match="negative phi-averaged density"):
        _run(channel="tagged-6Li-alpha", inclusive_b1=True, b1_reference_q2=0.0)
    p0, _ = _run(channel="tagged-6Li-alpha", inclusive_b1=True, b1_reference_q2=0.0, x_max=0.85)
    p1, _ = _run(channel="tagged-6Li-alpha", inclusive_b1=True, x_max=0.85)
    assert p0.dis_sampler.kernel.b1_reference_q2 == 0.0
    assert _sig(p1) != _sig(p0)
    with pytest.raises(RuntimeError, match="struck.b1_reference_q2"):
        cfg = lg.make_config(events=10)
        s = cfg.struck
        s.b1_reference_q2 = 2.5
        cfg.struck = s
        cfg.validate()


def test_the_cli_flag_and_banner(tmp_path):
    base = [sys.executable, "-m", "lipolgen.cli", "--events", "50", "--seed", "1"]
    r = subprocess.run(base + ["--npz", str(tmp_path / "a.npz")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "b1 reference scale (K1 2026-11-17, --b1-ref-q2): b1/F1 at Q0^2 = 2.5 GeV2" in r.stdout
    r = subprocess.run(base + ["--b1-ref-q2", "0", "--quiet"], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    r = subprocess.run(base + ["--b1-ref-q2", "1.0"], capture_output=True, text=True)
    assert r.returncode != 0
    assert "outside the unpolarised grid's own Q2 range" in (r.stderr + r.stdout)
    assert "Traceback" not in r.stderr
