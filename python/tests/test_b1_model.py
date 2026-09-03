"""The Python half of design D: the `--b1-model` surface of `b1_nuclear.hpp`.

P1  the binding and the CLI surface: `lg.B1Model`, `lg.B1_MODELS`, the
    `make_config` round trip, the three `cli.DEFAULTS` keys and the three
    `meta` keys a generated file must carry.
P2  the four terms sum to `b1`, and the two term knobs kill exactly their own.
P3  the band: `--b1-band-scale 0` vs `1` on the same seed leaves the
    SPIN-BLIND cell cross sections bit-identical, and at band 0 the
    tensor-thirds categories collapse onto each other while at band 1 they do
    not.
P4  no RNG in the kernel -- the grep the doctest binary cannot do, because
    only `LIPOLGEN_REFERENCE_DIR` is compiled in, not the source directory.

THE RULE ALL OF THESE GUARD.  `Li6ConvolutionB1` has NOT passed its A = 2
magnitude gate (docs/open_items/run_2026-09-02/phase_D_gate.md), and even if
it had, Q(6Li) = -0.0806(6) fm^2 against Q_d = +0.2859(3) fm^2 says the alpha-d
D wave enters the closest measured observable with the opposite sign to the
deuteron's and nearly cancels it.  So every number is {0, 1, 2} x b1 and the
band is not a formality -- which is what P1's meta keys and P3 exist to make
mechanically checkable.
"""

import json
import os

import numpy as np
import pytest

import lipolgen as lg
from lipolgen import cli
from lipolgen import _lipolgen as _l

N = 400
SEED = 7

# The CDKS camp's b1_d is a Q2 = 2.5 DIGITIZATION with no Q2 evolution, and in
# the topmost cell of the default window (x = 0.955) its b1/F1 reaches 3.3 for
# `cdks` and 5.6 for `li6-convolution` -- past where 1 + w_avg stays positive,
# so `InclusiveSampler` refuses the run.  Both opt-in backends therefore need
# a window that stops below it; Miller's b1 is a ratio model and does not.
X_MAX = 0.95


def _cfg(model="miller", band=1.0, w_ad=1.0, x_max=X_MAX, isotope="6Li",
         channel="inclusive", **kw):
    cfg = lg.make_config(isotope=isotope, channel=channel, config=1,
                         events=N, seed=SEED, b1_model=model,
                         b1_band_scale=band, b1_alpha_d_dwave_weight=w_ad,
                         **kw)
    if x_max is not None:
        sc = cfg.scenario           # Scenario comes back BY VALUE
        sc.x_max = x_max
        cfg.scenario = sc
    return cfg


def _pipeline(**kw):
    return _l.Pipeline(_cfg(**kw), lg.tensor_thirds_plan(0.0, 0.6))


@pytest.fixture(scope="module")
def conv():
    return lg.Li6ConvolutionB1()


# ------------------------------------------------------------------- P1

def test_b1_model_name_table():
    assert set(lg.B1_MODELS) == {"miller", "cdks", "li6-convolution"}
    assert lg.B1_MODELS["miller"] == _l.B1Model.Miller
    assert lg.B1_MODELS["cdks"] == _l.B1Model.Cdks
    assert lg.B1_MODELS["li6-convolution"] == _l.B1Model.Li6Convolution
    assert "B1_MODELS" in lg.__all__
    for key, val in lg.B1_MODELS.items():
        assert _l.b1_model_name(val) == key


def test_b1_symbols_are_importable():
    for name in ("B1Model", "Li6ConvolutionB1", "Li6ConvolutionOptions",
                 "DeuteronConvolutionB1", "LightConeDensities",
                 "ClusterPartialWave", "ConvolutionKinematics", "B1Landmarks",
                 "b1_landmarks", "b1_landmarks_of_table", "f1_cdks",
                 "cdks_b1_raw_per_nucleon", "alpha_d_quadrupole_fm2",
                 "li6_alpha_d_partial_waves", "b1_model_name",
                 "VMC_N_ALPHA_D_LI6"):
        assert hasattr(_l, name), name


def test_n_alpha_d_has_one_home():
    # design 2.1: VMC_P_D_LI6 is now EXPRESSED THROUGH VMC_N_ALPHA_D_LI6, so
    # the spectroscopic factor is defined once.  Both are the momentum file's
    # OWN printed norms; `VMC_S_ALPHA_D_LI6` is the file's TOTAL block and is a
    # third, different value that must not be "unified" with them.
    assert _l.VMC_N_ALPHA_D_LI6 == 0.80362 + 0.015861
    assert _l.VMC_P_D_LI6 == 0.015861 / _l.VMC_N_ALPHA_D_LI6
    assert _l.VMC_S_ALPHA_D_LI6 != _l.VMC_N_ALPHA_D_LI6


def test_make_config_round_trip():
    cfg = _cfg()
    assert cfg.b1_model == _l.B1Model.Miller          # the DEFAULT
    assert cfg.b1_band_scale == 1.0
    assert cfg.b1_alpha_d_dwave_weight == 1.0

    cfg = _cfg("li6-convolution", band=2.0, w_ad=0.5)
    assert cfg.b1_model == _l.B1Model.Li6Convolution
    assert cfg.b1_band_scale == 2.0
    assert cfg.b1_alpha_d_dwave_weight == 0.5
    cfg.validate()

    # the enum straight through, not only the string key
    assert (lg.make_config(events=N, b1_model=_l.B1Model.Cdks).b1_model
            == _l.B1Model.Cdks)
    with pytest.raises(ValueError):
        lg.make_config(events=N, b1_model="nope")


@pytest.mark.parametrize("model", ["cdks", "li6-convolution"])
def test_validate_refuses_the_illegal_combinations(model):
    # 6Li ONLY -- and the deuteron is the trap, because it is spin 1 too.
    # BOTH opt-in models: "cdks" is Li6B1's 6Li rank-2 transfer as much as the
    # convolution is, and on any other isotope or channel the flag would not
    # reach the rate but WOULD reach meta["b1_model"].
    with pytest.raises(RuntimeError):
        _cfg(model, isotope="d").validate()
    with pytest.raises(RuntimeError):
        _cfg(model, isotope="7Li").validate()
    # INCLUSIVE ONLY -- on a tagged channel the alpha-d density is already in
    # the event weight
    with pytest.raises(RuntimeError):
        _cfg(model, channel="tagged-alpha").validate()
    with pytest.raises(RuntimeError):
        _cfg(model, channel="coherent").validate()
    # a caller-supplied kernel and a non-default model contradict each other
    cfg = _cfg(model)
    cfg.kernel = _l.default_inclusive_kernel(_l.li6())
    with pytest.raises(RuntimeError):
        cfg.validate()
    # ... but the DEFAULT model with a caller-supplied kernel is legal, as it
    # has always been
    cfg.b1_model = _l.B1Model.Miller
    cfg.validate()


def test_validate_refuses_a_knob_that_did_not_run():
    """A knob the chosen backend never reads may not be recorded as if it had.

    `--b1-band-scale` is deliberately NOT applied to Miller (those numbers are
    the published ones) and the Cdks branch of `default_inclusive_kernel`
    never reads the alpha-d D-wave weight -- but the npz/HFS `meta` records
    all three keys unconditionally, so a run with either set would claim a
    variation it never made.  That is exactly the false provenance the three
    keys were added to prevent (design 4.2).
    """
    for band in (0.0, 2.0):
        with pytest.raises(RuntimeError):
            _cfg("miller", band=band).validate()
    with pytest.raises(RuntimeError):
        _cfg("miller", w_ad=2.0).validate()
    with pytest.raises(RuntimeError):
        _cfg("cdks", w_ad=2.0).validate()
    # ... and the same knobs ARE legal where they run
    _cfg("cdks", band=0.0).validate()
    _cfg("li6-convolution", band=0.0, w_ad=2.0).validate()
    # the >= 0 range check applies to every model, not only the opt-in ones
    # (it used to sit inside the non-Miller branch, so miller took -1)
    for model in ("miller", "cdks", "li6-convolution"):
        with pytest.raises(RuntimeError):
            _cfg(model, band=-1.0).validate()
        with pytest.raises(RuntimeError):
            _cfg(model, w_ad=-5.0).validate()


def test_cli_defaults_and_flags():
    for key in ("b1_model", "b1_band_scale", "b1_alpha_d_dwave_weight"):
        assert key in cli.DEFAULTS, key
    assert cli.DEFAULTS["b1_model"] == "miller"
    assert cli.DEFAULTS["b1_band_scale"] == 1.0
    assert cli.DEFAULTS["b1_alpha_d_dwave_weight"] == 1.0

    opts = cli.resolve(["--b1-model", "li6-convolution",
                        "--b1-band-scale", "2", "--x-max", "0.95"])
    assert opts["b1_model"] == "li6-convolution"
    assert opts["b1_band_scale"] == 2.0
    assert opts["b1_alpha_d_dwave_weight"] == 1.0
    assert opts["x_max"] == 0.95
    # and the default resolve is untouched
    assert cli.resolve([])["b1_model"] == "miller"
    assert cli.resolve([])["x_max"] is None

    help_text = cli.build_parser().format_help()
    for flag in ("--b1-model", "--b1-band-scale", "--b1-alpha-d-dwave-weight",
                 "--x-max"):
        assert flag in help_text, flag


def test_cli_says_which_flag_fixes_the_top_x_cell():
    """`--b1-model cdks` with no `--x-max` always failed, opaquely.

    What the user saw was the sampler's "negative phi-averaged density for
    m=1 at x = 0.955", which names neither `--x-max` nor the CDKS table.  A
    shipped flag whose plain invocation fails with an unrelated-looking
    message is a defect; `cli.main` now says which flag fixes it, before the
    pipeline is built.
    """
    for model in ("cdks", "li6-convolution"):
        with pytest.raises(SystemExit) as e:
            cli.main(["--b1-model", model, "--events", "10", "--quiet"])
        msg = str(e.value)
        assert "--x-max 0.95" in msg, msg
        assert model in msg
    # ... and the DEFAULT is untouched: no x_max rule for miller
    opts = cli.resolve([])
    assert opts["b1_model"] == "miller" and opts["x_max"] is None


def test_meta_records_which_b1_made_the_file():
    p = _pipeline(model="li6-convolution", band=2.0, w_ad=0.5)
    meta = p.generate(0, True)["meta"]
    assert meta["b1_model"] == "li6-convolution"
    assert meta["b1_band_scale"] == 2.0
    assert meta["b1_alpha_d_dwave_weight"] == 0.5
    # the default path still says what it is, so no npz is anonymous
    meta = _pipeline()  .generate(0, True)["meta"]
    assert meta["b1_model"] == "miller"
    assert meta["b1_band_scale"] == 1.0
    # ... and a hand-built kernel is NOT labelled with a flag that did not run
    cfg = _cfg()
    cfg.kernel = _l.default_inclusive_kernel(_l.li6())
    p = _l.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
    assert p.generate(0, True)["meta"]["b1_model"] == "caller-supplied kernel"


# ------------------------------------------------------------------- P2

def test_four_terms_sum_to_b1(conv):
    q2 = 2.5
    for x in (0.05, 0.10, 0.20, 0.30, 0.50):
        parts = (conv.b1_embedded_s(x, q2), conv.b1_alpha_d_dwave_d(x, q2),
                 conv.b1_alpha_d_dwave_alpha(x, q2), conv.b1_cg_dwave(x, q2))
        assert conv.b1(x, q2, 0.0) == pytest.approx(sum(parts), rel=1e-12)
        # (2d) + (2a) is the physical "alpha-d orbital" term
        assert (conv.b1_alpha_d_dwave(x, q2)
                == pytest.approx(parts[1] + parts[2], rel=1e-12))
        # ... and the struck-alpha piece is HALF of it, not a rounding error
        assert 0.45 < parts[2] / parts[1] < 0.60


def test_term_knobs_kill_exactly_their_own_term(conv):
    q2, x = 2.5, 0.30
    o = _l.Li6ConvolutionOptions()
    o.w_alpha_d_dwave = 0.0
    off = _l.Li6ConvolutionB1(o)
    assert off.b1_alpha_d_dwave_d(x, q2) == 0.0
    assert off.b1_alpha_d_dwave_alpha(x, q2) == 0.0
    assert off.b1_embedded_s(x, q2) == conv.b1_embedded_s(x, q2)
    assert off.b1_cg_dwave(x, q2) == conv.b1_cg_dwave(x, q2)


def test_banded_scales_the_whole_b1(conv):
    q2, x = 2.5, 0.30
    assert conv.banded(0.0).b1(x, q2, 0.0) == 0.0
    assert (conv.banded(2.0).b1(x, q2, 0.0)
            == 2.0 * conv.banded(1.0).b1(x, q2, 0.0))


# ------------------------------------------------------------------- P3

@pytest.mark.parametrize("model", ["cdks", "li6-convolution"])
def test_band_leaves_the_spin_blind_rate_bit_identical(model):
    p0 = _pipeline(model=model, band=0.0)
    p1 = _pipeline(model=model, band=1.0)
    c0 = np.asarray(p0.dis_sampler.cell_xsec_pb)
    c1 = np.asarray(p1.dis_sampler.cell_xsec_pb)
    # `cell_xsec_pb` is the SPIN-BLIND cell cross section: b1 enters only
    # through the per-state tensor shift, so the band must not move it at all.
    assert np.array_equal(c0, c1)


@pytest.mark.parametrize("model", ["cdks", "li6-convolution"])
def test_band_zero_collapses_the_tensor_categories(model):
    names = [c.name for c in _pipeline().plan.categories]
    s0 = np.asarray(_pipeline(model=model, band=0.0).sigma_per_category_pb())
    s1 = np.asarray(_pipeline(model=model, band=1.0).sigma_per_category_pb())
    s2 = np.asarray(_pipeline(model=model, band=2.0).sigma_per_category_pb())
    # At band 0 the tensor shift is identically zero, so the three
    # tensor-thirds categories are the same cross section (to summation
    # order); at band 1 the m = 0 category has separated from the m = +-1 pair.
    assert s0 == pytest.approx(s0[0], rel=1e-13)
    assert names.index("azz0") == 2
    assert s1[2] != s1[0]
    assert s1[0] == s1[1]                     # +-1 are P_zz-degenerate
    # ... and the band is LINEAR in the scale, which is what makes
    # "b1 +- 100 %" the envelope of the three runs rather than three
    # unrelated numbers.
    assert (s2 - s1) == pytest.approx(s1 - s0, rel=1e-6, abs=1e-9)


def test_default_kernel_is_the_miller_kernel(repo_root):
    """The one-argument overload IS `(Miller, 1, 1)`, through the binding.

    The SAME 600 rows T9 pins on the C++ side, at the same rtol 1e-12, read
    from `validation/reference/b1_default_li6.json` -- not a smoke test: any
    finite kernel passes `isfinite`, which is what this used to check.
    """
    path = os.path.join(repo_root, "validation", "reference",
                        "b1_default_li6.json")
    with open(path) as fh:
        doc = json.load(fh)
    assert doc["kernel"]["b1_model"] == "miller"
    k = _l.default_inclusive_kernel(_l.li6())
    xs = doc["x"]
    assert len(xs) == 200 and len(doc["tables"]) == 3
    n = 0
    for blk in doc["tables"]:
        q2 = blk["q2"]
        for i, x in enumerate(xs):
            t = k.tables(x, q2)
            for name in ("f1", "b1", "b2", "delta"):
                assert getattr(t, name) == pytest.approx(blk[name][i],
                                                         rel=1e-12, abs=0.0), (
                    "%s at x=%r q2=%r" % (name, x, q2))
            n += 1
    assert n == 600


# ------------------------------------------------------------------- P4

def test_no_rng_in_the_kernel(repo_root):
    """The grep T8 cannot do: the doctest binary knows only the reference dir.

    A structure function must be a pure function of (x, Q2) -- docs/
    CONVENTIONS.md, "no RNG in a structure function".  `b1_nuclear.cpp` is a
    quadrature over tabulated wave functions and has no business drawing a
    random number; if one ever appears, two runs with the same seed stop
    agreeing and nothing else in the suite would notice.
    """
    src = os.path.join(repo_root, "src", "core", "b1_nuclear.cpp")
    assert os.path.isfile(src), src
    text = open(src).read().lower()
    for token in ("rng", "random", "rand(", "drand", "mt19937", "uniform("):
        assert token not in text, "%s appears in b1_nuclear.cpp" % token
