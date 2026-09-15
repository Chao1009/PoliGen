"""The Python half of the tensor-sector radiative corrections (rc.hpp).

T16  the binding and the CLI surface: `lipolgen.RC`, the classes, the
     `make_config` round trip, the three npz columns on BOTH paths, and the
     `--rc-delta-low-x 0.19 <= 0.30` band ordering.
T6   the pipeline half of the bit-for-bit gate: an `--rc tensor-band` run is
     an `--rc off` run in every column but the three new ones.
T13  `--rc off` emits NO rc key at all, and `export.rc_columns` then returns
     exact ones.
T17  the weighted-mode block shape.

The rule these all guard: the RC weights are a SYSTEMATIC BAND and a
BACKGROUND, so they live on `Event.rc_weights` and never on `Event.weight`.
"""

import numpy as np
import pytest

import lipolgen as lg
from lipolgen import export
from lipolgen import _lipolgen as _l

N = 4000


def _cfg(rc="off", **kw):
    return lg.make_config(isotope="6Li", channel="inclusive", config=1,
                          events=N, seed=99, rc=rc, **kw)


@pytest.fixture(scope="module")
def rc_run():
    cfg = _cfg("tensor-band")
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
    return p, p.generate(0, True)


@pytest.fixture(scope="module")
def off_run():
    cfg = _cfg("off")
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
    return p, p.generate(0, True)


# ------------------------------------------------------------------- T16

def test_rc_name_table():
    assert set(lg.RC) == {"off", "tensor-band"}
    assert lg.RC["off"] == _l.RcMode.Off
    assert lg.RC["tensor-band"] == _l.RcMode.TensorBand
    assert "RC" in lg.__all__
    assert "RC_C0_SHAPES" in lg.__all__ and set(lg.RC_C0_SHAPES) == {"ho",
                                                                     "vmc-ft"}


def test_rc_symbols_are_importable():
    for name in ("RcMode", "RcScope", "RcTailModel", "RcOptions", "RcWeights",
                 "RcModel", "HoSpin1FF", "HoSpin1FFOptions",
                 "TabulatedSpin1FF", "Spin1ElasticFF", "rc_delta",
                 "rc_mode_name", "rc_weight_name", "nucleon_ff",
                 "C0Shape", "c0_shape_name"):
        assert hasattr(_l, name), name


def test_rc_delta_anchors_and_monotonicity():
    # The clamp branches return the constant unmodified, so `==` is legitimate.
    assert _l.rc_delta(0.01) == 0.30
    assert _l.rc_delta(1e-5) == 0.30
    assert _l.rc_delta(0.16) == 0.015
    assert _l.rc_delta(0.9) == 0.015
    assert _l.rc_delta(0.063) == pytest.approx(0.1108061822113555, rel=1e-12)
    xs = np.geomspace(1e-5, 0.9, 2000)
    d = np.array([_l.rc_delta(x) for x in xs])
    assert np.all(np.diff(d) <= 1e-15)


def test_rc_weight_names():
    assert _l.rc_weight_name(0, 0) == "rc_tensor_lo"
    assert _l.rc_weight_name(1, 0) == "rc_tensor_hi"
    assert _l.rc_weight_name(2, 0) == "rc_tail"
    assert _l.rc_weight_name(2, 3) == "rc_tail_3"
    assert export.RC_KEYS == ("rc_tensor_lo", "rc_tensor_hi", "rc_tail")


def test_make_config_round_trips_every_knob():
    cfg = _cfg("tensor-band", rc_delta_low_x=0.19, rc_delta_high_x=0.02,
               rc_a_transfer_frac=0.75,
               rc_fq_scale=2.0, rc_tail_tensor_scale=0.5,
               rc_qe_suppression=0.0, rc_qe_tensor_scale=2.5,
               rc_sp_tensor_scale=3.5, rc_tail_model="t-peak+ll",
               rc_c0_shape="vmc-ft")
    assert cfg.rc == _l.RcMode.TensorBand
    o = cfg.rc_options
    assert o.delta_low_x == 0.19
    assert o.delta_high_x == 0.02
    assert o.a_transfer_frac == 0.75
    assert o.fq_scale == 2.0
    assert o.tail_tensor_scale == 0.5
    assert o.qe_suppression == 0.0
    assert o.qe_tensor_scale == 2.5
    # The DEFAULT is 0: the shipped quasi-elastic tail is tensor-blind, and a
    # borrowed magnitude may not become one by accident.
    assert _cfg("tensor-band").rc_options.qe_tensor_scale == 0.0
    # ... and the s/p stand-in the same way.  It needs `t-peak+ll` above,
    # because the shipped `t-peak` computes no s-/p-peaks and validate()
    # refuses a price on a piece that did not run.
    assert o.sp_tensor_scale == 3.5
    assert _cfg("tensor-band").rc_options.sp_tensor_scale == 0.0
    assert o.c0_shape == _l.C0Shape.VmcFt
    # The DEFAULT edge is `Ho`: every published number was made with it.
    assert _cfg("tensor-band").rc_options.c0_shape == _l.C0Shape.Ho
    with pytest.raises(ValueError):
        _cfg("tensor-band", rc_c0_shape="vmc")
    # There is NO mode on RcOptions -- PipelineConfig.rc is the one source.
    assert not hasattr(o, "mode")
    with pytest.raises(ValueError):
        _cfg("nope")


def test_form_factor_anchors_and_refusals():
    ff = _l.HoSpin1FF.for_ion(_l.ion_by_name("6Li"))
    assert ff.fc(0.0) == pytest.approx(3.0, rel=1e-12)
    assert ff.fm(0.0) == pytest.approx(4.90765, rel=1e-5)
    assert ff.fq(0.0) == pytest.approx(-65.914, rel=1e-4)
    assert "MEASURED" in ff.provenance()
    # No measured-moment block for anything but 6Li -- a 7Li run must NOT
    # silently get 6Li form factors.
    with pytest.raises(RuntimeError):
        _l.HoSpin1FF.for_ion(_l.ion_by_name("7Li"))


def test_c0_shape_band_shares_its_anchors_and_flips_the_tail_sign():
    """B1: the C0 (monopole) shape is a BAND, and the band is not decorative.

    The two edges are built to agree exactly where the MEASURED input is --
    F_c(0) = Z, F_q(0) = M_A^2 Q and <r^2>_point = 6.0788 fm^2 -- and to
    disagree only in the SHAPE, which is the unfitted part.  What that
    disagreement costs is the sign of the tensor fraction of the elastic tail
    at x = 0.1, so this test gates both halves: the shared anchors, and the
    fact that `sigma^el_T` changes sign between the edges (which is why
    `phase_B_numbers.md` sec. B1 publishes that sign as INDETERMINATE).
    """
    ion = _l.ion_by_name("6Li")
    ho_o, vm_o = _l.HoSpin1FFOptions(), _l.HoSpin1FFOptions()
    vm_o.c0_shape = _l.C0Shape.VmcFt
    assert ho_o.c0_shape == _l.C0Shape.Ho          # the default
    ho, vm = (_l.HoSpin1FF.for_ion(ion, o) for o in (ho_o, vm_o))

    # (i) the MEASURED anchors are identical on both edges.
    for f in ("fc", "fm", "fq"):
        assert getattr(vm, f)(0.0) == pytest.approx(getattr(ho, f)(0.0),
                                                    rel=1e-12)
    # (ii) <r^2>_point, on the UNFOLDED shape, is the same 6.0788 fm^2.
    pt_o = [_l.HoSpin1FFOptions(), _l.HoSpin1FFOptions()]
    pt_o[1].c0_shape = _l.C0Shape.VmcFt
    for o in pt_o:
        o.fold_nucleon = False
        f = _l.HoSpin1FF.for_ion(ion, o)
        h = 1e-6
        t = h * _l.HBARC_GEV_FM ** 2
        slope = (f.fc(t) - f.fc(0.0)) / h
        assert -6.0 * slope / f.fc(0.0) == pytest.approx(6.0788, rel=1e-4)
    # (iii) and they disagree by a factor 2.5 in F_point where the tail lives.
    t2 = (2.0 * _l.HBARC_GEV_FM) ** 2
    assert vm.fc(t2) / ho.fc(t2) == pytest.approx(2.508, rel=1e-3)
    # (iv) the Ho edge has a C0 zero at q_0 = 3.0999 fm^-1 and the VmcFt edge
    # has NONE at any q -- the difference that moves sigma^el_T.
    t0 = (3.0998 * _l.HBARC_GEV_FM) ** 2
    assert abs(ho.fc(t0)) < 1e-6 * ho.fc(0.0)
    assert ho.fc((3.5 * _l.HBARC_GEV_FM) ** 2) < 0.0
    for q in (3.0, 3.0998, 3.5, 10.0, 88.0):
        assert vm.fc((q * _l.HBARC_GEV_FM) ** 2) > 0.0
    # ... and at the very top of the t-peak's reach it UNDERFLOWS to +0, which
    # is the continuation dying, not a sign change.
    assert vm.fc((278.0 * _l.HBARC_GEV_FM) ** 2) == 0.0
    assert "vmc-ft" in vm.provenance() and "NO C0 ZERO" in vm.provenance()
    assert "C0 shape = ho" in ho.provenance()
    assert _l.c0_shape_name(_l.C0Shape.VmcFt) == "vmc-ft"
    # (v) VmcFt is 6Li's OWN density and is REFUSED for any other ion rather
    # than silently handed over.
    d = _l.HoSpin1FFOptions()
    d.a_fm, d.alpha, d.mu_n, d.q_fm2 = 1.0, 0.1, 0.857, 0.2857
    d.fm_qz_fm, d.fm_b_fm = 1.3, 1.85
    _l.HoSpin1FF.for_ion(_l.ion_by_name("d"), d)          # legal on Ho
    d.c0_shape = _l.C0Shape.VmcFt
    with pytest.raises(RuntimeError):
        _l.HoSpin1FF.for_ion(_l.ion_by_name("d"), d)


def test_c0_shape_reaches_the_tail_and_the_npz_meta():
    """The knob is plumbed end to end: it reaches `sigma^el_T` through the
    Pipeline (not only `HoSpin1FF`), it FLIPS ITS SIGN at x = 0.1, and the npz
    `meta` records which edge ran -- without that key the two files would be
    indistinguishable."""
    out = {}
    for name in ("ho", "vmc-ft"):
        cfg = _cfg("tensor-band", rc_c0_shape=name)
        p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
        out[name] = p.rc_model.tail_sigma_at(0.10, 5.0)
        assert p.rc_model.options.c0_shape == lg.RC_C0_SHAPES[name]
        assert ("C0 shape = " + name) in p.rc_model.ff_provenance
    # sigma^el_U rises by 10 %; sigma^el_T CHANGES SIGN.
    assert out["vmc-ft"][0] / out["ho"][0] == pytest.approx(1.1043, rel=1e-3)
    assert out["ho"][1] > 0.0 > out["vmc-ft"][1]
    # the quasi-elastic tail is untouched: this knob is the ELASTIC shape only
    assert out["vmc-ft"][2] == pytest.approx(out["ho"][2], rel=1e-12)


def test_npz_columns_on_the_columnar_path(rc_run):
    p, cols = rc_run
    assert p.rc_model is not None
    for k in export.RC_KEYS:
        assert k in cols
    lo, hi, tail = (cols[k] for k in export.RC_KEYS)
    assert lo.shape == cols["x"].shape
    # w_lo + w_hi == 2 exactly: they are the two SIGNS of one delta.
    assert np.all(lo + hi == 2.0)
    assert np.all((lo > 0.5) & (lo < 2.0))
    assert np.all((hi > 0.5) & (hi < 2.0))
    assert np.all(tail >= 1.0)
    # The band is not trivially 1 -- otherwise this whole file gates nothing.
    assert np.any(hi != 1.0)
    # THE rule: the RC weights are NOT on the nominal weight.
    assert np.all(cols["weight"] == 1.0)
    assert cols["meta"]["rc"] == "tensor-band"
    assert cols["meta"]["rc_weight_names"] == list(export.RC_KEYS)


def test_meta_records_every_rc_knob_that_moves_a_column(rc_run):
    """PROVENANCE: two npz files that differ in any RC knob must differ in
    `meta`.  The seven numeric knobs were recorded from the start; `scope`,
    `with_qe_tail`, `tail_model`, `n_eta` and `tail_max` were not, so a
    `RcScope.TensorAll` / `with_qe_tail = False` run -- both API-reachable,
    and the QRT is 73-99.9 % of the tail at x >= 0.1 -- was byte-
    indistinguishable in `meta` from the default one.  Same rule the Miller
    branch enforces for `b1_band_scale`."""
    _, cols = rc_run
    meta = cols["meta"]
    defaults = _l.RcOptions()
    assert set(k for k in meta if k.startswith("rc")) == {
        "rc", "rc_weight_names", "rc_applies", "rc_tail_applies",
        "rc_exclusion_reason", "rc_ff_provenance",
        "rc_delta_low_x", "rc_delta_high_x", "rc_a_transfer_frac",
        "rc_fq_scale",
        "rc_c0_shape",
        "rc_tail_tensor_scale", "rc_qe_suppression", "rc_qe_tensor_scale",
        "rc_sp_tensor_scale", "rc_qe_kf_gev",
        "rc_band_tau_max",
        "rc_scope", "rc_with_qe_tail", "rc_tail_model", "rc_n_eta",
        "rc_tail_max",
        "rc_clipped_cell_fraction", "rc_clipped_fraction_by_y",
        "rc_clipped_band_events", "rc_clipped_tail_events",
        "rc_clipped_band_event_fraction", "rc_clipped_tail_event_fraction",
    }
    assert meta["rc_c0_shape"] == "ho"
    # `rc_scope` carries the LABEL under this fixture's `tensor_thirds_plan`,
    # whose categories sit at theta_S = 0 where the two scopes are one run
    # (test_meta_distinguishes_a_tensor_all_no_qe_run measures it).  The name
    # goes in wherever the knob is read.
    assert meta["rc_scope"] == meta["knob_provenance"]["rc_scope"]["label"]
    assert meta["knob_provenance"]["rc_scope"]["value"] == "tensor-rate"
    assert meta["rc_with_qe_tail"] is True
    assert meta["rc_tail_model"] == "t-peak"
    assert meta["rc_n_eta"] == defaults.n_eta
    assert meta["rc_tail_max"] == defaults.tail_max


def test_meta_distinguishes_a_tensor_all_no_qe_run(rc_run):
    """The two knobs the numeric block could not see, each written down --
    and `rc_scope` only where it RAN.

    2026-09-05.  This test used to set `scope = TensorAll` under
    `tensor_thirds_plan` and assert `meta["rc_scope"] == "tensor-all"` and
    "it really is a different file, not only a different label" -- while
    checking only the label.  It is NOT a different file there: `TensorAll`
    differs from `TensorRate` only in the cos 2phi amplitude, which
    `tensor_amplitudes` builds with a sin^2(theta_S) factor, and every
    `tensor_thirds_plan` category sits at theta_S = 0.  Measured (400 events,
    seed 7, inclusive 6Li): bit-identical in all 47 columns and every sigma.
    So the knob-provenance table calls it NOT READ under that plan and
    `meta["rc_scope"]` carries the label; under `transverse_tensor_plan`
    (theta_S = pi/2) it runs, the file moves, and the name goes in.
    """
    _, base = rc_run
    # (a) theta_S = 0: it did NOT run, and the file says so.
    cfg = _cfg("tensor-band")
    cfg.rc_options.scope = _l.RcScope.TensorAll
    cfg.rc_options.with_qe_tail = False
    flat = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
    meta = flat.generate(0)["meta"]
    assert meta["rc_with_qe_tail"] is False
    row = meta["knob_provenance"]["rc_scope"]
    assert row["status"] == "not-read"
    assert meta["rc_scope"] == row["label"] != "tensor-all"
    # ... and the run really is bit for bit the tensor-rate one, which is why.
    ref = _cfg("tensor-band")
    ref.rc_options.with_qe_tail = False
    flat0 = lg.Pipeline(ref, lg.tensor_thirds_plan(0.0, 0.6)).generate(0)
    for k in export.RC_KEYS:
        assert np.array_equal(flat.generate(0)[k], flat0[k]), k

    # (b) theta_S = pi/2: it DID run, and now the name goes in.
    tilt = _cfg("tensor-band")
    tilt.rc_options.scope = _l.RcScope.TensorAll
    tmeta = lg.Pipeline(tilt, lg.transverse_tensor_plan(0.6)).generate(0)["meta"]
    assert tmeta["rc_scope"] == "tensor-all"
    assert tmeta["knob_provenance"]["rc_scope"]["status"] == "read"
    assert tmeta["rc_scope"] != base["meta"]["rc_scope"]


def test_m_lepton_is_refused_on_the_t_peak_tails_and_READ_on_polrad_full():
    """`RcOptions::m_lepton` was RESERVED for `PolradFull` and refused on
    every model until 2026-09-06, because nothing read it.  `PolradFull`
    reads it now -- it is the m^2 of Eq. (B.13)'s C_1,2(tau), of
    F_IR = m^2 F_2+ - Q_m^2 F_d and of lambda_s, i.e. what sets the WIDTH of
    the s- and p-peaks -- so the rule is unchanged and its answer moved: a
    knob is refused where it does not run and read where it does."""
    for model in ("t-peak", "t-peak+ll"):
        cfg = _cfg("tensor-band", rc_tail_model=model)
        assert cfg.rc_options.m_lepton == _l.RcOptions().m_lepton
        cfg.rc_options.m_lepton = 0.1056583755      # a muon beam, say
        with pytest.raises(RuntimeError, match="m_lepton"):
            lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
    plan = lg.tensor_thirds_plan(0.0, 0.6)
    cfg = _cfg("tensor-band", rc_tail_model="polrad-full")
    base = lg.Pipeline(cfg, plan).rc_model
    cfg2 = _cfg("tensor-band", rc_tail_model="polrad-full")
    cfg2.rc_options.m_lepton = 0.1056583755
    muon = lg.Pipeline(cfg2, plan).rc_model
    # A heavier lepton NARROWS the s-/p-peaks, so the tail falls.  Asserted as
    # an inequality: what is gated is that the field reaches the quadrature.
    assert muon.tail_ratio_at(1e-4, 0.2, 0.0) < base.tail_ratio_at(1e-4, 0.2, 0.0)
    # ... and a non-positive lepton mass is refused everywhere: it is a mass
    # squared in C_1,2(tau) and in lambda_s = S^2 - 4 m^2 M^2.
    cfg3 = _cfg("tensor-band", rc_tail_model="polrad-full")
    cfg3.rc_options.m_lepton = 0.0
    with pytest.raises(RuntimeError, match="m_lepton"):
        lg.Pipeline(cfg3, plan)


def test_qe_tensor_scale_prices_the_polarised_quasi_elastic_tail(rc_run):
    """B3.  The polarised quasi-elastic tail is not computed anywhere; the
    quasi-elastic piece is 73 % of `rc_tail` at x = 0.1 and 99.9 % at x = 0.3
    and was treated as exactly tensor-blind.  `qe_tensor_scale` lends it the
    ELASTIC tail's own tensor fraction -- a BORROWED magnitude, not a derived
    bound (see rc.hpp).  Gated here as a chain, because the default is 0 and
    a knob whose npz `meta` did not record it would be invisible."""
    _, base = rc_run
    assert base["meta"]["rc_qe_tensor_scale"] == 0.0
    assert _l.RcOptions().qe_tensor_scale == 0.0

    cfg = _cfg("tensor-band", rc_qe_tensor_scale=1.0)
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
    meta = p.generate(0)["meta"]
    assert meta["rc_qe_tensor_scale"] == 1.0

    b = lg.Pipeline(_cfg("tensor-band"),
                    lg.tensor_thirds_plan(0.0, 0.6)).rc_model
    m = p.rc_model
    for x in (0.01, 0.10, 0.30):
        # It carries q_n and NOTHING else: the unpolarised tail is untouched.
        assert m.tail_ratio_at(x, 5.0, 0.0) == b.tail_ratio_at(x, 5.0, 0.0)
        # ... and the tensor part moves, in the direction the elastic tensor
        # fraction points (whose SIGN is a C0-shape band edge -- magnitude
        # only).
        d0 = b.tail_ratio_at(x, 5.0, 1.0) - b.tail_ratio_at(x, 5.0, 0.0)
        d1 = m.tail_ratio_at(x, 5.0, 1.0) - m.tail_ratio_at(x, 5.0, 0.0)
        assert d1 != d0
        assert abs(d1) > abs(d0)
    # EXACTLY LINEAR in the scale -- one run rescales, unlike fq_scale (T12).
    m2 = lg.Pipeline(_cfg("tensor-band", rc_qe_tensor_scale=2.0),
                     lg.tensor_thirds_plan(0.0, 0.6)).rc_model
    for x in (0.01, 0.30):
        b1 = b.tail_ratio_at(x, 5.0, 1.0)
        a1 = m.tail_ratio_at(x, 5.0, 1.0) - b1
        a2 = m2.tail_ratio_at(x, 5.0, 1.0) - b1
        assert a2 / a1 == pytest.approx(2.0, rel=1e-10)


def test_qe_tensor_scale_is_refused_when_the_qe_tail_did_not_run():
    """The `m_lepton` rule: a knob that did not run may not be recorded as if
    it had.  `qe_tensor_scale` is a fraction OF Eq. (44)'s sigma^q_U, so with
    `with_qe_tail = False` there is nothing for it to be a fraction of."""
    cfg = _cfg("tensor-band", rc_qe_tensor_scale=1.0)
    cfg.rc_options.with_qe_tail = False
    with pytest.raises(RuntimeError, match="qe_tensor_scale"):
        lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
    # ... zero with the tail off is fine.
    ok = _cfg("tensor-band")
    ok.rc_options.with_qe_tail = False
    lg.Pipeline(ok, lg.tensor_thirds_plan(0.0, 0.6))
    # ... and a negative scale is refused like every other >= 0 knob, by
    # `PipelineConfig::validate()` itself -- `make_config` already calls it.
    with pytest.raises(RuntimeError, match="qe_tensor_scale"):
        _cfg("tensor-band", rc_qe_tensor_scale=-1.0)


def test_sp_tensor_scale_bounds_the_tensor_fraction_of_the_sp_peaks(rc_run):
    """B1 (2026-09-06).  `--rc-tail-model t-peak+ll` puts the leading-log s-
    and p-peaks in the UNPOLARISED numerator only -- POLRAD's *Eq. (38)*
    supplies no tensor s/p peak, and Eq. (18) + Eq. (A.4) does supply one
    which `polrad-full` computes, so the unqualified sentence is true of
    Eq. (38) and false of the paper, and this scale is refused on
    `polrad-full` for the OPPOSITE reason (the term ran) -- so the upper edge
    of the tail band LOWERS the tensor fraction
    of the tail purely by growing its denominator, and the tensor fraction of
    that piece was left at exactly zero.  `sp_tensor_scale` prices it by
    lending the ELASTIC s-/p-peaks the elastic t-peak's own tensor fraction.
    IT IS A BOUND WITH NO DERIVATION (see rc.hpp), and it is EMPTY wherever
    6Li's coherent form factor is dead at the s/p vertex -- which includes all
    three standard points.  Both halves are gated here."""
    _, base = rc_run
    assert base["meta"]["rc_sp_tensor_scale"] == 0.0
    assert _l.RcOptions().sp_tensor_scale == 0.0

    cfg = _cfg("tensor-band", rc_tail_model="t-peak+ll",
               rc_sp_tensor_scale=1.0)
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
    meta = p.generate(0)["meta"]
    assert meta["rc_sp_tensor_scale"] == 1.0
    assert meta["rc_tail_model"] == "t-peak+ll"

    b = lg.Pipeline(_cfg("tensor-band", rc_tail_model="t-peak+ll"),
                    lg.tensor_thirds_plan(0.0, 0.6)).rc_model
    m = p.rc_model
    s = p.dis_sampler.s

    # (a) THE BOUND IS EMPTY AT THE THREE STANDARD POINTS, bit for bit.  The
    # s-/p-peak's own elastic vertex sits at Q'^2 = z_s Q^2 = 4.37 / 4.94 /
    # 4.98 GeV^2 there, four decades above the t-peak's t ~ t_min ~ 8.7e-05
    # GeV^2, and 6Li's coherent form factor is dead at that scale
    # (F_c = -3.75e-45 against +2.99).  This is not "small"; it is nothing.
    for x in (0.01, 0.10, 0.30):
        for q_n in (0.0, 1.0, -2.0):
            assert m.tail_ratio_at(x, 5.0, q_n) == b.tail_ratio_at(x, 5.0, q_n)

    # (b) ... AND IT IS NOT EMPTY WHERE THE COHERENT s/p PEAK IS ALIVE: the
    # low-x, y -> 1 corner.  It carries q_n and nothing else, so the
    # UNPOLARISED tail never moves.
    xs = np.asarray(p.dis_sampler.x_cells)
    q2s = np.asarray(p.dis_sampler.q2_cells)
    moved = 0
    for x, q2 in zip(xs, q2s):
        assert m.tail_ratio_at(x, q2, 0.0) == b.tail_ratio_at(x, q2, 0.0)
        if m.tail_ratio_at(x, q2, 1.0) != b.tail_ratio_at(x, q2, 1.0):
            moved += 1
            assert q2 / (x * s) > 0.3        # every one of them is at high y
            assert x < 0.01                  # ... and at low x
    assert moved > 0

    # (c) EXACTLY LINEAR in the scale -- one run rescales, unlike fq_scale.
    m2 = lg.Pipeline(_cfg("tensor-band", rc_tail_model="t-peak+ll",
                          rc_sp_tensor_scale=2.0),
                     lg.tensor_thirds_plan(0.0, 0.6)).rc_model
    checked = 0
    for x, q2 in zip(xs, q2s):
        b1 = b.tail_ratio_at(x, q2, 1.0)
        a1 = m.tail_ratio_at(x, q2, 1.0) - b1
        # `a1` is a DIFFERENCE of two tail ratios whose common size is the
        # whole numerator, so isolating it costs one digit per decade by which
        # the numerator exceeds it.  The algebra is linear everywhere; it is
        # gated where a double can still see that.
        if abs(a1) < 1e-6 * abs(b1):
            continue
        a2 = m2.tail_ratio_at(x, q2, 1.0) - b1
        # rel 1e-10: MEASURED worst 8.2435e-11 on 2026-09-06 over the 6 cells
        # of this fixture's grid that clear the guard.  The C++ case gates the
        # same identity AT A NODE, over 270 of them, where no interpolation
        # enters at all.
        assert a2 / a1 == pytest.approx(2.0, rel=1e-10)
        checked += 1
    assert checked > 0


def test_sp_tensor_scale_is_refused_when_the_sp_peaks_did_not_run():
    """The `m_lepton` rule again, one level below `qe_tensor_scale`'s.
    `sp_tensor_scale` is a fraction OF the leading-log s-/p-peak column, and
    only `tail_model = t-peak+ll` computes that column: under the shipped
    `t-peak` the table is identically zero, so a price there would be recorded
    in `meta` without a single floating-point operation behind it."""
    with pytest.raises(RuntimeError, match="sp_tensor_scale"):
        _cfg("tensor-band", rc_sp_tensor_scale=1.0)
    # ... zero under the shipped tail is fine: nothing is recorded that did
    # not run.
    assert _cfg("tensor-band").rc_options.sp_tensor_scale == 0.0
    # ... and a negative scale is refused like every other >= 0 knob.
    with pytest.raises(RuntimeError, match="sp_tensor_scale"):
        _cfg("tensor-band", rc_tail_model="t-peak+ll",
             rc_sp_tensor_scale=-1.0)
    # The C++ API path is refused by `RcModel`'s own constructor too, so a
    # caller who reaches past `make_config` gets the same answer.
    cfg = _cfg("tensor-band", rc_tail_model="t-peak+ll",
               rc_sp_tensor_scale=1.0)
    cfg.rc_options.tail_model = _l.RcTailModel.TPeak
    with pytest.raises(RuntimeError, match="sp_tensor_scale"):
        lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))


def test_the_quasi_elastic_sp_column_is_bounded_by_NEITHER_scale():
    """What B1 does NOT cover, gated so it cannot be quietly forgotten.
    `qe_tensor_scale` multiplies sigma^q_U alone and `sp_tensor_scale` the
    ELASTIC s/p column alone, so the QUASI-ELASTIC s/p column keeps a tensor
    part of exactly zero -- and it is the one that survives at high Q^2, where
    the coherent form factor is dead.  The test asserts the gap by showing
    that at the three standard points, where the coherent s/p peak is dead but
    `t-peak+ll` still grows the tail, BOTH scales at once leave the tensor
    term exactly where `qe_tensor_scale` alone puts it."""
    plan = lg.tensor_thirds_plan(0.0, 0.6)
    ll = lg.Pipeline(_cfg("tensor-band", rc_tail_model="t-peak+ll"),
                     plan).rc_model
    qe = lg.Pipeline(_cfg("tensor-band", rc_tail_model="t-peak+ll",
                          rc_qe_tensor_scale=1.0), plan).rc_model
    both = lg.Pipeline(_cfg("tensor-band", rc_tail_model="t-peak+ll",
                            rc_qe_tensor_scale=1.0,
                            rc_sp_tensor_scale=1.0), plan).rc_model
    tpk = lg.Pipeline(_cfg("tensor-band"), plan).rc_model
    for x in (0.01, 0.10, 0.30):
        # the s/p scale adds nothing here ...
        assert both.tail_ratio_at(x, 5.0, 1.0) == qe.tail_ratio_at(x, 5.0, 1.0)
        # ... while `t-peak+ll` HAS grown the unpolarised tail, so the tensor
        # FRACTION of the tail really did fall and stays fallen.
        assert ll.tail_ratio_at(x, 5.0, 0.0) > tpk.tail_ratio_at(x, 5.0, 0.0)
        f_tpk = ((tpk.tail_ratio_at(x, 5.0, 1.0) -
                  tpk.tail_ratio_at(x, 5.0, 0.0)) /
                 tpk.tail_ratio_at(x, 5.0, 0.0))
        f_ll = ((ll.tail_ratio_at(x, 5.0, 1.0) -
                 ll.tail_ratio_at(x, 5.0, 0.0)) /
                ll.tail_ratio_at(x, 5.0, 0.0))
        assert abs(f_ll) < abs(f_tpk)


def test_a_transfer_frac_prices_the_A2_to_A6_transfer_of_the_band(rc_run):
    """B6, design Q8.  `delta(x)` is the DOMINANT RC systematic and every
    anchor it interpolates between is a DEUTERON number -- HERMES's measured
    low-x residual, Gakh-Shekhovtsova's 10-30 %, E12-13-011's 1.5 %.  No
    A > 2 tensor RC calculation exists at all, so the shipped band silently
    ASSUMES the deuteron fractional RC transfers to 6Li.  `a_transfer_frac`
    adds that doubt in quadrature: delta_eff = delta * sqrt(1 + f^2)."""
    _, base = rc_run
    # The default is 0 and it is RECORDED, so a widened band is never
    # mistaken in `meta` for the shipped one.
    assert base["meta"]["rc_a_transfer_frac"] == 0.0
    assert _l.RcOptions().a_transfer_frac == 0.0

    b = lg.Pipeline(_cfg("tensor-band"),
                    lg.tensor_thirds_plan(0.0, 0.6)).rc_model
    xs = (1e-4, 0.005, 0.01, 0.02, 0.05, 0.063, 0.1, 0.16, 0.3, 0.9)
    # At the default the model's delta IS the published rc_delta, bit for bit.
    for x in xs:
        assert b.delta(x) == _l.rc_delta(x)

    for f in (0.5, 1.0, 2.0):
        cfg = _cfg("tensor-band", rc_a_transfer_frac=f)
        p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
        assert p.generate(0)["meta"]["rc_a_transfer_frac"] == f
        m = p.rc_model
        k = (1.0 + f * f) ** 0.5
        for x in xs:
            assert m.delta(x) == pytest.approx(k * _l.rc_delta(x), rel=1e-14)
        # It moves the BAND only: the tail is a background and does not see it.
        for x in (0.01, 0.10, 0.30):
            for q_n in (0.0, 1.0, -2.0):
                assert m.tail_ratio_at(x, 5.0, q_n) == \
                    b.tail_ratio_at(x, 5.0, q_n)
        # ... and the free function -- the PUBLISHED shape, what T3 gates --
        # is untouched by any of this.
        assert _l.rc_delta(0.01) == _l.RC_DELTA_LOW_X

    # A negative fraction is refused rather than squared away by the hypot.
    with pytest.raises(RuntimeError, match="a_transfer_frac"):
        _cfg("tensor-band", rc_a_transfer_frac=-0.5)


def test_tagged_tail_exclusion_names_the_quasi_elastic_omission():
    """B4.  `rc_tail == 1` on a tagged channel is HALF a kinematic fact.  The
    ELASTIC recoil is vetoed (x_L = 1, inside the beam envelope); the
    QUASI-ELASTIC one is not -- the A-1 remnant is unbound and its alpha lands
    at x_L ~ 2/3, the tag window itself.  The run has to say both, because an
    analysis reading `rc_tail == 1` out of the npz would otherwise take it for
    a veto on the whole tail."""
    cfg = lg.make_config(isotope="6Li", channel="tagged-alpha", config=1,
                         events=200, seed=11, rc="tensor-band")
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
    why = p.rc_model.exclusion_reason
    assert not p.rc_model.tail_applies
    for token in ("x_L", "FACT", "OMISSION", "QUASI-elastic", "2/3",
                  "SPECTATOR"):
        assert token in why, token
    # It travels into the npz, where the analysis actually reads it.
    assert p.generate(0)["meta"]["rc_exclusion_reason"] == why


def test_tail_model_knob_is_a_seven_file_chain(rc_run):
    """`RcTailModel::TPeakPlusLL` -- the OPT-IN leading-log s-/p-peak tail.

    Gated end to end because the knob is a chain and any missing link makes an
    npz that is indistinguishable from a default run: the enum, `make_config`,
    the CLI's own choices tuple, the `meta` string, and the numbers actually
    changing.  The DEFAULT is asserted unchanged in the same test, because
    that is the whole point of making it opt-in."""
    assert set(lg.RC_TAIL_MODELS) == {"t-peak", "t-peak+ll", "polrad-full"}
    assert "RC_TAIL_MODELS" in lg.__all__
    assert lg.RC_TAIL_MODELS["t-peak"] == _l.RcTailModel.TPeak
    assert lg.RC_TAIL_MODELS["t-peak+ll"] == _l.RcTailModel.TPeakPlusLL
    assert lg.RC_TAIL_MODELS["polrad-full"] == _l.RcTailModel.PolradFull
    # ONE spelling, and it is the C++ one.
    for name, e in lg.RC_TAIL_MODELS.items():
        assert _l.rc_tail_model_name(e) == name
    assert _l.rc_tail_model_name(_l.RcTailModel.PolradFull) == "polrad-full"

    _, base = rc_run
    assert base["meta"]["rc_tail_model"] == "t-peak"

    cfg = _cfg("tensor-band", rc_tail_model="t-peak+ll")
    assert cfg.rc_options.tail_model == _l.RcTailModel.TPeakPlusLL
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
    cols = p.generate(0, True)
    assert cols["meta"]["rc_tail_model"] == "t-peak+ll"

    # The kinematics are IDENTICAL -- `RcModel::fill` takes no Rng& and the
    # tail model may not move the random stream (T6's rule).
    assert np.array_equal(cols["x"], base["x"])
    assert np.array_equal(cols["q2"], base["q2"])
    assert np.array_equal(cols["weight"], base["weight"])
    # ... and the tail is strictly larger, event by event.  Never smaller:
    # the s- and p-peaks are positive wherever the radiator is defined.
    assert np.all(cols["rc_tail"] >= base["rc_tail"] - 1e-15)
    assert cols["rc_tail"].mean() > base["rc_tail"].mean()

    # The three t-peak tables are untouched and the two s+p tables are
    # non-zero -- and identically zero on the default.
    m0, m1 = lg.Pipeline(_cfg("tensor-band"),
                         lg.tensor_thirds_plan(0.0, 0.6)).rc_model, p.rc_model
    for k in ("sigma_tail_u", "sigma_tail_t", "sigma_tail_qe"):
        assert np.array_equal(getattr(m1, k), getattr(m0, k))
    for k in ("sigma_tail_u_sp", "sigma_tail_qe_sp"):
        assert not np.any(getattr(m0, k))
        assert np.any(getattr(m1, k))
    # tail_sigma_at is five wide, and the last two are the s+p peaks.
    assert len(m1.tail_sigma_at(0.10, 5.0)) == 5
    assert m0.tail_sigma_at(0.10, 5.0)[3:] == (0.0, 0.0)

    with pytest.raises(ValueError, match="unknown rc_tail_model"):
        _cfg("tensor-band", rc_tail_model="eq-18")

    # ... and the last link: the CLI flag, its default, and its refusal of a
    # name that is not a model.  Without this an npz written by `lipolgen-run`
    # could not reach the knob at all.
    from lipolgen import cli
    assert cli.DEFAULTS["rc_tail_model"] == "t-peak"
    assert cli.resolve(["--rc", "tensor-band"])["rc_tail_model"] == "t-peak"
    assert cli.resolve(["--rc", "tensor-band", "--rc-tail-model",
                        "t-peak+ll"])["rc_tail_model"] == "t-peak+ll"
    assert cli.resolve(["--rc", "tensor-band", "--rc-tail-model",
                        "polrad-full"])["rc_tail_model"] == "polrad-full"
    with pytest.raises(SystemExit):
        cli.resolve(["--rc-tail-model", "eq-18"])


def test_polrad_full_is_the_exact_tail_and_does_not_bracket_the_t_peak_band(
        rc_run):
    """B2.  `RcTailModel::PolradFull` -- POLRAD Eq. (18) + Appendix B +
    Eq. (A.4), the exact tau_A quadrature, opt-in since 2026-09-06.

    THREE THINGS THIS ASSERTS THAT THE C++ SUITE CANNOT.  (i) the chain --
    enum, `RC_TAIL_MODELS`, `make_config`, the CLI choices tuple, the `meta`
    string -- because a missing link makes an npz indistinguishable from a
    default run.  (ii) that the DEFAULT is untouched, event by event.
    (iii) that the model does NOT lie inside the [t-peak, t-peak+ll] band on
    the production grid, which is the one thing a reader would assume from
    "the two are a band"."""
    plan = lg.tensor_thirds_plan(0.0, 0.6)
    _, base = rc_run
    cfg = _cfg("tensor-band", rc_tail_model="polrad-full")
    assert cfg.rc_options.tail_model == _l.RcTailModel.PolradFull
    p = lg.Pipeline(cfg, plan)
    cols = p.generate(0, True)
    assert cols["meta"]["rc_tail_model"] == "polrad-full"
    # The kinematics are IDENTICAL: `RcModel::fill` takes no Rng&, so no tail
    # model may move the random stream (T6's rule).
    assert np.array_equal(cols["x"], base["x"])
    assert np.array_equal(cols["q2"], base["q2"])
    assert np.array_equal(cols["weight"], base["weight"])
    # It ADDS the s- and p-peaks, so the tail never falls below the t-peak.
    assert np.all(cols["rc_tail"] >= base["rc_tail"] - 1e-15)
    assert cols["rc_tail"].mean() > base["rc_tail"].mean()

    # Eq. (18) has NO s/p decomposition: all three peaks are in one tau_A
    # integral, so the two leading-log columns stay zero for a DIFFERENT
    # reason than under t-peak -- there is nothing to read there, not
    # something empty.
    mf = p.rc_model
    for k in ("sigma_tail_u_sp", "sigma_tail_qe_sp"):
        assert not np.any(getattr(mf, k))
    assert mf.tail_sigma_at(0.10, 5.0)[3:] == (0.0, 0.0)
    # THE PIN.  `RcModel::tail_sigma_at` under this model is where the
    # historical factor-A bug lived (the shipped `per_nucleon = m_p/M_A` was
    # 6.00x too large; polrad_transcription_check.md sec. 4), and until
    # T19(i)/this block NOTHING pinned a VALUE there -- replacing the branch's
    # `1/A^2` by `1/A` left the whole C++ suite and both polrad-full pytests
    # green.  These are the C++ T19(i) numbers read back THROUGH the bindings,
    # so a reduction that is right in C++ and wrong in the binding, or a
    # binding that hands back the wrong column, is caught here and only here.
    # 6Li config 1, RcOptions() defaults: sigma per nucleon [GeV^-2], the
    # elastic pair = Eq. (18)/A^2 and the quasi-elastic = Eq. (A.5)/A.
    for x, q2, u, t, qe in (
            (0.01,  5.0,  2.2301252856028388e-05, -6.8964416598016638e-08,
                          1.4082400494855565e-05),
            (0.10, 20.0,  1.3474277075174561e-08,  1.0398154385161343e-11,
                          4.2214155843320782e-08),
            (0.30,  5.0,  2.5848200248542323e-12,  1.6976986824052815e-13,
                          9.5240888233400158e-06)):
        got = mf.tail_sigma_at(x, q2)
        assert got[0] == pytest.approx(u, rel=1e-9), (x, q2, "sigma^el_U")
        assert got[1] == pytest.approx(t, rel=1e-9), (x, q2, "sigma^el_T")
        assert got[2] == pytest.approx(qe, rel=1e-9), (x, q2, "sigma^q_U")
    # 2.230125e-05 is 8.0284510282e-04/36 and NOT /6: the x6 the pin exists
    # to catch is bigger than every tolerance above by seven decades.
    assert mf.tail_sigma_at(0.01, 5.0)[0] * 36.0 == pytest.approx(
        8.0284510281702199e-04, rel=1e-9)

    # ... and the s/p tensor STAND-IN is refused on it, for the OPPOSITE
    # reason it is refused on t-peak: the term it stands in for RAN.
    with pytest.raises(RuntimeError, match="double-count"):
        lg.Pipeline(_cfg("tensor-band", rc_tail_model="polrad-full",
                         rc_sp_tensor_scale=1.0), plan)

    # A resolution that cannot resolve the peaks is refused, not degraded.
    bad = _cfg("tensor-band", rc_tail_model="polrad-full")
    bad.rc_options.n_eta = 32
    with pytest.raises(RuntimeError, match="n_eta"):
        lg.Pipeline(bad, plan)

    # THE HEADLINE.  Over the sampler's own accepted cells the exact tail is
    # OUTSIDE the [t-peak, t-peak+ll] interval on a large minority of them --
    # measured 1326 of 3051 (43.5 %) on the production grid -- so the two
    # t-peak models are a PRICE RANGE and not a confidence interval.
    mt = lg.Pipeline(_cfg("tensor-band"), plan).rc_model
    ml = lg.Pipeline(_cfg("tensor-band", rc_tail_model="t-peak+ll"),
                     plan).rc_model
    ds = p.dis_sampler
    outside = 0
    for i in range(ds.n_cells):
        x, q2 = ds.x_cells[i], ds.q2_cells[i]
        a = mt.tail_ratio_at(x, q2, 0.0)
        b = ml.tail_ratio_at(x, q2, 0.0)
        f = mf.tail_ratio_at(x, q2, 0.0)
        if f < min(a, b) or f > max(a, b):
            outside += 1
    assert outside > 0.2 * ds.n_cells
    # ... and the QUASI-ELASTIC tensor gap is NOT closed by it.  A nucleon has
    # no tensor structure function (Eq. (A.5) -> Im_5..8 = 0), so the
    # quasi-elastic tail is tensor-blind at all three peaks here too and
    # `qe_tensor_scale` is still the only stand-in for it -- which is exactly
    # what makes it still ACCEPTED on this model.
    both = lg.Pipeline(_cfg("tensor-band", rc_tail_model="polrad-full",
                            rc_qe_tensor_scale=1.0), plan).rc_model
    assert both.tail_ratio_at(0.10, 5.0, 1.0) != mf.tail_ratio_at(0.10, 5.0, 1.0)


def test_npz_columns_on_the_records_path(rc_run):
    # `export.columns_from_events` builds its own dict and does NOT go through
    # the C++ `columns_to_dict`; without its own rule a records-path sample
    # from an RC-on run would carry no rc_* columns while `rc_columns()`
    # silently returned ones.
    p, cols = rc_run
    rec = export.columns_from_events(cols["events"], p.optics, p.pot_config)
    for k in export.RC_KEYS:
        assert k in rec
        assert np.allclose(rec[k], cols[k], rtol=0, atol=0)


def test_rc_columns_helper(rc_run, off_run):
    _, cols = rc_run
    lo, hi, tail = export.rc_columns(cols)
    assert np.all(lo + hi == 2.0)
    _, off_cols = off_run
    lo0, hi0, tail0 = export.rc_columns(off_cols)
    # Exactly ones, so an analysis can multiply them in unconditionally.
    assert np.all(lo0 == 1.0) and np.all(hi0 == 1.0) and np.all(tail0 == 1.0)
    assert len(lo0) == len(off_cols["x"])


def test_optimistic_delta_never_widens_the_band():
    """`--rc-delta-low-x 0.19` is inside the 0.30 band, everywhere at x < 0.01.

    Strictly inside only where tau != 0: at a b1 zero crossing both edges are
    exactly 0 and the strict form would fail on those events.
    """
    plan = lg.tensor_thirds_plan(0.0, 0.6)
    wide = lg.Pipeline(_cfg("tensor-band", rc_delta_low_x=0.30), plan)
    tight = lg.Pipeline(_cfg("tensor-band", rc_delta_low_x=0.19), plan)
    a = wide.generate(0)
    b = tight.generate(0)
    assert np.all(a["x"] == b["x"])
    lowx = a["x"] < 0.01
    assert lowx.sum() > 100
    da = np.abs(a["rc_tensor_hi"] - 1.0)[lowx]
    db = np.abs(b["rc_tensor_hi"] - 1.0)[lowx]
    assert np.all(db <= da + 1e-15)
    nonzero = da > 0.0
    assert nonzero.sum() > 0
    assert np.all(db[nonzero] < da[nonzero])


# -------------------------------------------------------------------- T6

def test_rc_on_moves_nothing_the_run_already_had(rc_run, off_run):
    p_on, on = rc_run
    p_off, off = off_run
    keys = [k for k in off
            if isinstance(off[k], np.ndarray) and k not in export.RC_KEYS]
    assert "weight" in keys and "x" in keys and "kp" in keys
    for k in keys:
        a, b = off[k], on[k]
        assert a.shape == b.shape, k
        if a.dtype.kind == "f":
            assert np.array_equal(a, b, equal_nan=True), k
        else:
            assert np.array_equal(a, b), k
    # ... and the records too, four-vector by four-vector.
    for ea, eb in zip(off["events"], on["events"]):
        assert ea.weight == eb.weight
        assert len(ea.particles) == len(eb.particles)
        for pa, pb in zip(ea.particles, eb.particles):
            assert (pa.p.e, pa.p.px, pa.p.py, pa.p.pz) == \
                   (pb.p.e, pb.p.px, pb.p.py, pb.p.pz)
        assert ea.rc_weights == []
        assert len(eb.rc_weights) == 3


# ------------------------------------------------------------------- T13

def test_off_emits_no_rc_key_at_all(off_run):
    p, cols = off_run
    assert p.rc_model is None
    for k in export.RC_KEYS:
        assert k not in cols
    # The whole rc meta block is absent, so an --rc off npz is byte-identical
    # to today's and not merely key-compatible with it.
    assert "rc" not in cols["meta"]
    assert not [k for k in cols["meta"] if k.startswith("rc")]
    rec = export.columns_from_events(cols["events"], p.optics, p.pot_config)
    for k in export.RC_KEYS:
        assert k not in rec


def test_coherent_channel_is_allowed_and_prices_nothing():
    cfg = lg.make_config(isotope="6Li", channel="coherent", config=1,
                         events=500, seed=3, rc="tensor-band")
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
    assert p.rc_model is not None
    assert not p.rc_model.applies
    assert not p.rc_model.tail_applies
    assert p.rc_model.exclusion_reason
    cols = p.generate(0)
    for k in export.RC_KEYS:
        assert np.all(cols[k] == 1.0)


def test_tagged_band_applies_but_the_tail_is_exactly_one():
    cfg = lg.make_config(isotope="6Li", channel="tagged-alpha", config=1,
                         events=20000, seed=4, rc="tensor-band")
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
    assert p.rc_model.applies
    assert not p.rc_model.tail_applies
    assert p.rc_model.exclusion_reason
    cols = p.generate(0)
    assert np.all(cols["rc_tail"] == 1.0)
    assert np.any(cols["rc_tensor_hi"] != 1.0)
    assert np.all(cols["rc_tensor_lo"] + cols["rc_tensor_hi"] == 2.0)
    # THE PUBLISHED BAND EDGES ARE POSITIVE AND BOUNDED.  tau_tag =
    # 1 - nbar/n_M is unbounded at the nodes of the M-dependent spectator
    # density, and without RcOptions.band_tau_max this run shipped
    # rc_tensor_hi down to -1.79 -- negative weights in the npz.  An analysis
    # multiplying weight * rc_tensor_hi must never get one.
    for k in ("rc_tensor_lo", "rc_tensor_hi"):
        assert cols[k].min() > 0.0, (k, cols[k].min())
        assert cols[k].max() <= 2.0, (k, cols[k].max())
        # ... and, tighter, inside [1 - delta, 1 + delta] with delta <= 0.30.
        assert np.all(np.abs(cols[k] - 1.0) <= _l.RC_DELTA_LOW_X + 1e-12)
    # The clamp actually bites here, and the run says how often.
    meta = cols["meta"]
    assert meta["rc_band_tau_max"] == _l.RC_BAND_TAU_MAX
    assert meta["rc_clipped_band_events"] > 0
    assert 0.0 < meta["rc_clipped_band_event_fraction"] < 0.05
    assert meta["rc_clipped_tail_events"] == 0        # no tail on this channel


def test_7li_tagged_band_is_reachable_from_the_API_and_bounded():
    """The `--rc tensor-band` band on 7Li.

    From the CLI there is exactly ONE route -- `--plan helicity-flip --pe 0`,
    where lam_e*pe = 0 satisfies RcModel's unpolarised-beam check; every
    other CLI plan is refused (helicity plans at pe != 0 by RcModel, spin-1
    tensor plans by Pipeline against a J = 3/2 channel ion).  An explicit
    (P_z, T) J = 3/2 fill needs the API, and the band is then bounded exactly
    as on 6Li -- the unclamped version reached -8.35 here."""
    cats = [_l.SpinCategory("m32", 1.5, [0.5, 0.0, 0.0, 0.5]),
            _l.SpinCategory("m12", 1.5, [0.0, 0.5, 0.5, 0.0])]
    plan = _l.RunPlan(cats, 0.0, 0.0, 0.6)
    cfg = lg.make_config(isotope="7Li", channel="tagged-alpha", config=1,
                         events=20000, seed=11, rc="tensor-band")
    p = lg.Pipeline(cfg, plan)
    assert p.rc_model.applies
    assert not p.rc_model.tail_applies
    cols = p.generate(0)
    assert np.all(cols["rc_tail"] == 1.0)
    assert np.any(cols["rc_tensor_hi"] != 1.0)
    for k in ("rc_tensor_lo", "rc_tensor_hi"):
        assert cols[k].min() > 0.0, (k, cols[k].min())
        assert np.all(np.abs(cols[k] - 1.0) <= _l.RC_DELTA_LOW_X + 1e-12)
    assert cols["meta"]["rc_clipped_band_events"] > 0
    # THE CLI ROUTE, exactly as USAGE sec. 7b and the --rc help now state it:
    # helicity-flip at pe = 0 DOES build a working band.
    zero = lg.Pipeline(cfg, lg.make_plan("helicity-flip", j=1.5, pz=0.7,
                                         pzz=0.6, pe=0.0))
    assert zero.rc_model.applies
    assert np.any(zero.generate(0)["rc_tensor_hi"] != 1.0)
    # ... and it is the ONLY one: every CLI plan at pe = 0.7 is refused.
    for name in sorted(set(lg.PLANS)):
        try:
            q = lg.make_plan(name, j=1.5, pz=0.7, pzz=0.6, pe=0.7)
        except Exception:
            continue                       # not a J = 3/2 plan at all
        try:
            lg.Pipeline(cfg, q)
        except RuntimeError:
            continue                       # refused, as documented
        raise AssertionError(
            "USAGE sec. 7b says every CLI plan at pe != 0 is refused on the "
            "7Li band, but --plan %s builds one" % name)


def test_a_python_subclass_form_factor_survives_the_pipeline():
    """T16: the documented use case of the Spin1ElasticFF trampoline.

    `RcOptions.ff` must keep the PYTHON object alive, not just the C++
    trampoline: with a plain def_readwrite the natural `opt.ff = MyFF()`
    dropped the Python half and the first `fc()` call threw
    `Tried to call pure virtual function "Spin1ElasticFF::fc"`."""
    import gc

    class ConstantFF(_l.Spin1ElasticFF):
        def fc(self, t):
            return 3.0

        def fm(self, t):
            return 4.9

        def fq(self, t):
            return -65.9

        def provenance(self):
            return "test: a constant Python-subclass form factor"

    opt = _l.RcOptions()
    opt.ff = ConstantFF()
    cfg = _cfg("tensor-band")
    cfg.rc_options = opt
    del opt                       # the ONLY strong reference left is the C++ one
    gc.collect()
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
    assert "Python-subclass" in p.rc_model.ff_provenance
    assert p.rc_model.ff_provenance == p.rc_model.options.ff.provenance()
    cols = p.generate(0)          # would raise "pure virtual" before the fix
    assert np.all(np.isfinite(cols["rc_tail"]))
    assert np.all(cols["rc_tail"] >= 1.0)
    # None puts the default back.
    opt2 = _l.RcOptions()
    opt2.ff = None
    assert opt2.ff is None


# ------------------------------------------------------------------- T17

def test_weight_block_shape_is_three_times_one_slot(rc_run):
    """`Pipeline` never fills `Event.spin_weights` (a latent capability), so
    n_slot is 1 and the block is exactly the three slot-0 weights.  Slot 0 is
    the event's OWN PURE state; the per-category mixture slots are gated in
    C++ (`tests/test_rc.cpp`, T17), which is where a plan with populated
    `spin_weights` can be built by hand."""
    p, cols = rc_run
    for ev in cols["events"][:50]:
        assert ev.spin_weights == []
        assert len(ev.rc_weights) == 3 * (1 + len(ev.spin_weights))
        w = p.rc_model.weights(ev)
        assert ev.rc_weights[0] == w.lo
        assert ev.rc_weights[1] == w.hi
        assert ev.rc_weights[2] == w.tail
        # The per-category mixture is reachable and agrees on a pure category.
        for k in range(len(p.plan.categories)):
            wk = p.rc_model.weights(ev, k)
            assert wk.lo + wk.hi == 2.0


# ------------------------------------------------------------------- CLI

def test_cli_exposes_the_rc_switches():
    from lipolgen import cli
    # CONVENTIONS.md "no physics number defined twice": assert against the
    # MODULE CONSTANTS, not against re-typed literals, so a move of the band
    # anchor in rc.hpp cannot leave the CLI silently overriding it.
    for k, v in (("rc", "off"),
                 ("rc_delta_low_x", _l.RC_DELTA_LOW_X),
                 ("rc_delta_high_x", _l.RC_DELTA_HIGH_X),
                 ("rc_fq_scale", 1.0),
                 ("rc_tail_tensor_scale", 1.0), ("rc_qe_suppression", 1.0),
                 ("rc_qe_tensor_scale", 0.0),
                 ("rc_sp_tensor_scale", 0.0)):
        assert cli.DEFAULTS[k] == v
    opts = cli.resolve(["--rc", "tensor-band", "--rc-delta-low-x", "0.19",
                        "--rc-fq-scale", "2", "--rc-tail-tensor-scale", "0.5",
                        "--rc-qe-suppression", "0",
                        "--rc-qe-tensor-scale", "1",
                        "--rc-sp-tensor-scale", "0.5"])
    assert opts["rc"] == "tensor-band"
    assert opts["rc_delta_low_x"] == 0.19
    assert opts["rc_fq_scale"] == 2.0
    assert opts["rc_tail_tensor_scale"] == 0.5
    assert opts["rc_qe_suppression"] == 0.0
    assert opts["rc_qe_tensor_scale"] == 1.0
    assert opts["rc_sp_tensor_scale"] == 0.5
