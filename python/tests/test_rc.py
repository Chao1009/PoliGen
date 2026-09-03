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


def test_rc_symbols_are_importable():
    for name in ("RcMode", "RcScope", "RcTailModel", "RcOptions", "RcWeights",
                 "RcModel", "HoSpin1FF", "HoSpin1FFOptions",
                 "TabulatedSpin1FF", "Spin1ElasticFF", "rc_delta",
                 "rc_mode_name", "rc_weight_name", "nucleon_ff"):
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
               rc_fq_scale=2.0, rc_tail_tensor_scale=0.5,
               rc_qe_suppression=0.0)
    assert cfg.rc == _l.RcMode.TensorBand
    o = cfg.rc_options
    assert o.delta_low_x == 0.19
    assert o.delta_high_x == 0.02
    assert o.fq_scale == 2.0
    assert o.tail_tensor_scale == 0.5
    assert o.qe_suppression == 0.0
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
        "rc_delta_low_x", "rc_delta_high_x", "rc_fq_scale",
        "rc_tail_tensor_scale", "rc_qe_suppression", "rc_qe_kf_gev",
        "rc_band_tau_max",
        "rc_scope", "rc_with_qe_tail", "rc_tail_model", "rc_n_eta",
        "rc_tail_max",
        "rc_clipped_cell_fraction", "rc_clipped_fraction_by_y",
        "rc_clipped_band_events", "rc_clipped_tail_events",
        "rc_clipped_band_event_fraction", "rc_clipped_tail_event_fraction",
    }
    assert meta["rc_scope"] == "tensor-rate"
    assert meta["rc_with_qe_tail"] is True
    assert meta["rc_tail_model"] == "t-peak"
    assert meta["rc_n_eta"] == defaults.n_eta
    assert meta["rc_tail_max"] == defaults.tail_max


def test_meta_distinguishes_a_tensor_all_no_qe_run(rc_run):
    """The two knobs the numeric block could not see, each written down."""
    _, base = rc_run
    cfg = _cfg("tensor-band")
    cfg.rc_options.scope = _l.RcScope.TensorAll
    cfg.rc_options.with_qe_tail = False
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
    meta = p.generate(0)["meta"]
    assert meta["rc_scope"] == "tensor-all"
    assert meta["rc_with_qe_tail"] is False
    # ... and it really is a different file, not only a different label.
    assert meta["rc_scope"] != base["meta"]["rc_scope"]


def test_m_lepton_is_reserved_and_refused_on_the_shipped_tail():
    """`RcOptions::m_lepton` is RESERVED for the unimplemented PolradFull
    tail: nothing in src/core/rc.cpp reads it, so setting it would be a
    silent no-op that `meta` does not record.  Refused, exactly as
    `b1_band_scale != 1` is on the Miller branch."""
    cfg = _cfg("tensor-band")
    assert cfg.rc_options.m_lepton == _l.RcOptions().m_lepton
    cfg.rc_options.m_lepton = 0.1056583755          # a muon beam, say
    with pytest.raises(RuntimeError, match="m_lepton"):
        lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))


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

def test_cli_exposes_the_six_switches():
    from lipolgen import cli
    # CONVENTIONS.md "no physics number defined twice": assert against the
    # MODULE CONSTANTS, not against re-typed literals, so a move of the band
    # anchor in rc.hpp cannot leave the CLI silently overriding it.
    for k, v in (("rc", "off"),
                 ("rc_delta_low_x", _l.RC_DELTA_LOW_X),
                 ("rc_delta_high_x", _l.RC_DELTA_HIGH_X),
                 ("rc_fq_scale", 1.0),
                 ("rc_tail_tensor_scale", 1.0), ("rc_qe_suppression", 1.0)):
        assert cli.DEFAULTS[k] == v
    opts = cli.resolve(["--rc", "tensor-band", "--rc-delta-low-x", "0.19",
                        "--rc-fq-scale", "2", "--rc-tail-tensor-scale", "0.5",
                        "--rc-qe-suppression", "0"])
    assert opts["rc"] == "tensor-band"
    assert opts["rc_delta_low_x"] == 0.19
    assert opts["rc_fq_scale"] == 2.0
    assert opts["rc_tail_tensor_scale"] == 0.5
    assert opts["rc_qe_suppression"] == 0.0
