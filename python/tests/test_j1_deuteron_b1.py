# SPDX-License-Identifier: GPL-3.0-or-later
"""J1 (2026-11-10; decision J1, THE AUTHOR'S, docs/open_items/run_2026-11-10/
DECISIONS.md; record phase_P2_deuteron_b1.md): on an inclusive DEUTERON run the
kernel's rank-2 slot is the deuteron camp itself, MillerB1 per nucleon,
unmapped (until then Li6B1(MillerB1) = 0.32753 x MillerB1, the 6Li transfer);
6Li keeps Li6B1(MillerB1).  A moved b1_transfer on the deuteron RAISES, the
default is accepted and `meta["li6_b1_rank2_transfer"]` carries a `not read`
label (A17); the pre-J1 reading stays constructible from Python through
`InclusiveKernel.Options.b1_func`."""
import pytest

import lipolgen as lg
import lipolgen._lipolgen as _l

PTS = [(x, q2) for x in (0.012, 0.1, 0.5, 0.9, 0.955) for q2 in (2.5, 100.0, 1897.0)]
needs_lhapdf = pytest.mark.skipif(not _l.HAVE_LHAPDF, reason="no LHAPDF tier: the shipped default cannot be built")
FACTOR = 1.0 / (_l.LI6_B1_PER_NUCLEON * _l.LI6_B1_RANK2_TRANSFER)   # 3.0531989383...


def test_the_class_toy_deuteron_reads_millerb1_itself_and_6li_the_transfer():
    m, li6 = _l.MillerB1(), _l.Li6B1(_l.MillerB1())
    kd, k6 = _l.default_inclusive_kernel(lg.deuteron()), _l.default_inclusive_kernel(lg.li6())
    for x, q2 in PTS:
        td, t6 = kd.tables(x, q2), k6.tables(x, q2)
        assert td.b1 == m.b1(x, q2, td.f1)
        assert td.b2 == 2.0 * x * td.b1
        assert t6.b1 == li6.b1(x, q2, t6.f1)
        assert td.b1 / li6.b1(x, q2, td.f1) == pytest.approx(FACTOR, rel=1e-14)
    assert FACTOR == pytest.approx(3.0531989383, rel=1e-10)


@needs_lhapdf
def test_the_shipped_deuteron_kernel_is_the_pipelines_and_reads_millerb1():
    """The kernel a default deuteron run builds (D7 grids) carries MillerB1."""
    m = _l.MillerB1()
    k = lg.shipped_inclusive_kernel(lg.deuteron())
    cfg = lg.make_config(isotope="d", events=10)
    sc = cfg.scenario
    sc.x_max = 0.95   # J1 (2026-11-10): the default x_max 1.0 is refused on the deuteron (phase P2 sec. 1, M3)
    cfg.scenario = sc
    kp = _l.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6)).dis_sampler.kernel
    for x, q2 in PTS:
        assert k.tables(x, q2).b1 == m.b1(x, q2, 0.0) == kp.tables(x, q2).b1


@needs_lhapdf
def test_a_moved_transfer_on_the_deuteron_raises_and_the_default_is_accepted():
    lg.shipped_inclusive_kernel(lg.deuteron())                                  # None = the default
    lg.shipped_inclusive_kernel(lg.deuteron(), b1_transfer=_l.LI6_B1_RANK2_TRANSFER)
    with pytest.raises(RuntimeError, match="b1_transfer = 0.92194699999999996 on the deuteron is refused"):
        lg.shipped_inclusive_kernel(lg.deuteron(), b1_transfer=_l.LI6_B1_RANK2_TRANSFER_HULTHEN)
    with pytest.raises(RuntimeError, match="b1_model = cdks on the deuteron is refused"):
        _l.default_inclusive_kernel(lg.deuteron(), _l.B1Model.Cdks)
    # 6Li: the polligen gates' explicit 0.921947 still builds, and is Li6B1 at it
    k = _l.default_inclusive_kernel(lg.li6(), b1_transfer=_l.LI6_B1_RANK2_TRANSFER_HULTHEN)
    want = _l.Li6B1(_l.MillerB1(), _l.LI6_B1_RANK2_TRANSFER_HULTHEN).b1(0.3, 5.0, 0.0)
    assert k.tables(0.3, 5.0).b1 == want


def test_the_pre_j1_reading_is_constructible_from_python():
    old = _l.Li6B1(_l.MillerB1())
    o = _l.InclusiveKernelOptions() if hasattr(_l, "InclusiveKernelOptions") else _l.InclusiveKernel.Options()
    o.b1_func = lambda x, q2, f1: old.b1(x, q2, f1)
    k = _l.InclusiveKernel(lg.deuteron(), o)
    assert k.tables(0.3, 5.0).b1 == old.b1(0.3, 5.0, 0.0)
    assert k.tables(0.3, 5.0).b1 * FACTOR == pytest.approx(
        _l.default_inclusive_kernel(lg.deuteron()).tables(0.3, 5.0).b1, rel=1e-14)


@needs_lhapdf
def test_the_deuteron_run_labels_the_transfer_not_read_and_6li_keeps_the_number():
    def meta(iso):
        cfg = lg.make_config(isotope=iso, events=20)
        sc = cfg.scenario
        sc.x_max = 0.95   # J1 (2026-11-10): the deuteron is refused at 1.0
        cfg.scenario = sc
        p = _l.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6))
        return p, p.generate(0, False, 1)["meta"]
    p, m = meta("d")
    assert m["li6_b1_rank2_transfer"].startswith("not read on isotope d")
    assert "J1, 2026-11-10" in m["li6_b1_rank2_transfer"]
    f = _l.inclusive_polarization_factors(p)
    assert f["li6_b1_rank2_transfer"] == m["li6_b1_rank2_transfer"]
    assert "the 6Li rank-2 transfer is not read (J1, 2026-11-10" in f["report"]
    assert "on the deuteron b1_func is MillerB1 itself" in m["rank2_input"]
    p6, m6 = meta("6Li")
    assert m6["li6_b1_rank2_transfer"] == _l.LI6_B1_RANK2_TRANSFER
    assert "deuteron" not in m6["rank2_input"]


@needs_lhapdf
def test_every_default_deuteron_plan_is_refused_at_x_max_1_by_the_rank2_term():
    """M3 of phase P2, pinned: the exposure J1 opens (the treatment is the
    supervisor's, by a dated note under J1).  Default tensor-thirds, config 1."""
    cfg = lg.make_config(isotope="d", events=10)
    with pytest.raises(RuntimeError) as e:
        _l.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6))
    msg = str(e.value)
    assert "x = 0.955, Q2 = 167.3" in msg and "-1.13 from the RANK-2 term" in msg
    assert "24 cell(s) of this spin state are refused" in msg
