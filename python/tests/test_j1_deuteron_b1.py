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
        u = kd.unconfined_b1(x, q2, td.f1)   # J1 note (2026-11-12): the slot before the spin-1 confinement
        assert u == m.b1(x, q2, td.f1)
        assert td.b1 == min(max(u, -1.5 * td.f1), 3.0 * (td.f1 - abs(td.g1)))   # J1 note (2026-11-12): binds at x 0.955, Q2 2.5 on the toy
        assert td.b2 == 2.0 * x * td.b1
        assert t6.b1 == li6.b1(x, q2, t6.f1)
        assert u / li6.b1(x, q2, td.f1) == pytest.approx(FACTOR, rel=1e-14)
    assert FACTOR == pytest.approx(3.0531989383, rel=1e-10)


@needs_lhapdf
def test_the_shipped_deuteron_kernel_is_the_pipelines_and_reads_millerb1():
    """The kernel a default deuteron run builds (D7 grids) carries MillerB1."""
    m = _l.MillerB1()
    k = lg.shipped_inclusive_kernel(lg.deuteron())
    cfg = lg.make_config(isotope="d", events=10)
    sc = cfg.scenario
    sc.x_max = 1.0   # J1 note (2026-11-12): builds at the default again (J1 (2026-11-10): 0.95, refused at 1.0)
    cfg.scenario = sc
    kp = _l.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6)).dis_sampler.kernel
    for x, q2 in PTS:
        t = k.tables(x, q2)
        assert k.unconfined_b1(x, q2, t.f1) == m.b1(x, q2, 0.0) == kp.unconfined_b1(x, q2, t.f1)   # J1 note (2026-11-12)
        assert t.b1 == kp.tables(x, q2).b1 == min(max(m.b1(x, q2, 0.0), -1.5 * t.f1), 3.0 * (t.f1 - abs(t.g1)))


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
        sc.x_max = 1.0   # J1 note (2026-11-12): the default again (J1 (2026-11-10): 0.95, the deuteron was refused at 1.0)
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
    # J1 note (2026-11-12): the confinement's key, on the deuteron only
    assert m["b1_clipped_rate_fraction"] == pytest.approx(3.871134487837413e-07, rel=1e-9)   # config 1, 59 cells (phase P2b N4)
    assert "confined to the spin-1 leading-order positivity domain" in m["rank2_input"]
    assert "b1_clipped_rate_fraction" not in m6


@needs_lhapdf
def test_every_default_deuteron_plan_is_refused_at_x_max_1_by_the_rank2_term():
    """M3 of phase P2, pinned: the exposure J1 opens.  J1 note (2026-11-12):
    the supervisor's treatment confines the default deuteron slot, so the
    DEFAULT run builds now (its key has its value) and the exposure is pinned on
    the UNCONFINED kernel built through the public API -- the same CT18NLO /
    NNPDFpol1.1 objects, MillerB1, the 1e-2 Delta, the option off -- which
    reproduces P2's refusal.  Default tensor-thirds, config 1."""
    p = _l.Pipeline(lg.make_config(isotope="d", events=10), lg.tensor_thirds_plan(0.7, 0.6))
    assert _l.b1_rank2_confinement(p)["n_cells"] == 59
    o = _l.InclusiveKernel.Options()
    o.f2_source = _l.make_lhapdf_unpol_sf(_l.UNPOL_SF_DEFAULT)
    o.g1_model = _l.make_lhapdf_pol_sf(_l.POL_SF_DEFAULT)
    mb = _l.MillerB1()
    o.b1_func = lambda x, q2, f1: mb.b1(x, q2, f1)
    o.delta_func = lambda x, q2, f1: _l.toy_delta_gluon(x, q2, f1, 1e-2)
    assert o.confine_spin1_rank2 is False
    cfg = lg.make_config(isotope="d", events=10)
    cfg.kernel = _l.InclusiveKernel(lg.deuteron(), o)
    with pytest.raises(RuntimeError) as e:
        _l.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6))
    msg = str(e.value)
    assert "x = 0.955, Q2 = 167.3" in msg and "-1.13 from the RANK-2 term" in msg
    assert "24 cell(s) of this spin state are refused" in msg


def test_the_confinement_option_from_python():
    """J1 note (2026-11-12): `InclusiveKernel.Options.confine_spin1_rank2` is bound,
    off by default, and on it confines both halves; the factory sets it on the
    deuteron only."""
    o = _l.InclusiveKernel.Options()
    assert o.confine_spin1_rank2 is False
    o.b1_func = lambda x, q2, f1: 10.0 * f1
    off = _l.InclusiveKernel(lg.deuteron(), o)
    o.confine_spin1_rank2 = True
    on = _l.InclusiveKernel(lg.deuteron(), o)
    t0, t1 = off.tables(0.5, 10.0), on.tables(0.5, 10.0)
    assert t0.b1 == 10.0 * t0.f1 and not off.confine_spin1_rank2
    assert t1.b1 == 3.0 * (t1.f1 - abs(t1.g1)) and on.confine_spin1_rank2
    o.b1_func = lambda x, q2, f1: -10.0 * f1
    assert _l.InclusiveKernel(lg.deuteron(), o).tables(0.5, 10.0).b1 == -1.5 * t1.f1
    assert _l.default_inclusive_kernel(lg.deuteron()).confine_spin1_rank2
    assert not _l.default_inclusive_kernel(lg.li6()).confine_spin1_rank2
    with pytest.raises(RuntimeError, match="SPIN-1 positivity domain"):
        _l.InclusiveKernel(lg.li7(), o)
