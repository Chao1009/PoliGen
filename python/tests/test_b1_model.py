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

P6  the b1 normalisation split of 2026-09-03: the two camps' constants, and
    that `CdksB1` is now the raw column while `MillerB1` is still halved.
P7  the r-default decision: the null hook IS `r_sigma_lt`, term (1) is
    bit-identical under the swap to `r1998`, and the measured cost of the
    choice.

THE RULE ALL OF THESE GUARD.  `Li6ConvolutionB1` PASSED its A = 2 magnitude
gate on 2026-09-03 (G3b ratio 0.843 with CDKS's own MSTW2008 LO at their
Eq. (21); docs/open_items/run_2026-09-03/phase_A_numbers.md) -- but that is an
A = 2 gate and it tests NOTHING about the alpha-d step, and
Q(6Li) = -0.0818(17) fm^2 (LI6_QUADRUPOLE_FM2, TUNL; Pyykko's compilation
gives -0.0806(6)) against Q_d = +0.2859(3) fm^2 says the alpha-d D wave enters
the closest measured observable with the opposite sign to the deuteron's and
nearly cancels it.  So every number is {0, 1, 2} x b1 and the band is not a
formality -- which is what P1's meta keys and P3 exist to make mechanically
checkable.
"""

import json
import os
import sys

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


@pytest.fixture(scope="module")
def band_pipelines():
    """The 2 models x 3 band rows P3 needs, built ONCE.

    A `li6-convolution` pipeline costs ~4.9 s to construct, and the two P3
    tests below between them asked for five of them with two duplicated -- 29
    of this file's 43 s.  The assertions are unchanged; only the construction
    is shared."""
    return {(m, b): _pipeline(model=m, band=b)
            for m in ("cdks", "li6-convolution")
            for b in (0.0, 1.0, 2.0)}


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
                 "VMC_N_ALPHA_D_LI6", "CdBonnWave", "cdbonn_wave",
                 "cdbonn_fdeut_table", "DeuteronWaveSource"):
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
def test_band_leaves_the_spin_blind_rate_bit_identical(model, band_pipelines):
    p0 = band_pipelines[(model, 0.0)]
    p1 = band_pipelines[(model, 1.0)]
    c0 = np.asarray(p0.dis_sampler.cell_xsec_pb)
    c1 = np.asarray(p1.dis_sampler.cell_xsec_pb)
    # `cell_xsec_pb` is the SPIN-BLIND cell cross section: b1 enters only
    # through the per-state tensor shift, so the band must not move it at all.
    assert np.array_equal(c0, c1)


@pytest.mark.parametrize("model", ["cdks", "li6-convolution"])
def test_band_zero_collapses_the_tensor_categories(model, band_pipelines):
    names = [c.name for c in band_pipelines[(model, 1.0)].plan.categories]
    s0 = np.asarray(band_pipelines[(model, 0.0)].sigma_per_category_pb())
    s1 = np.asarray(band_pipelines[(model, 1.0)].sigma_per_category_pb())
    s2 = np.asarray(band_pipelines[(model, 2.0)].sigma_per_category_pb())
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


# ------------------------------------------------------------------- P5
# The CD-Bonn deuteron wave function (A5): `lg.CdBonnWave`,
# `lg.cdbonn_fdeut_table` and `DeuteronConvolutionB1.Options.wave`.
#
# WHAT THESE GUARD.  `cluster.hpp` types twenty published coefficients
# (Machleidt, PRC 63 (2001) 024001, Table XX) and computes four more from the
# r -> 0 boundary conditions.  A mistyped digit would not fail anything else
# in the library -- it would just be a slightly different deuteron, and the
# A = 2 gate would move by an amount a reader would take for physics.  The
# doctest in tests/test_cluster.cpp is the primary gate; these repeat its
# sharpest clauses across the binding, and add the one thing C++ cannot see:
# that the OPTION is reachable from Python and that the DEFAULT did not move.

def test_cdbonn_reproduces_machleidts_published_deuteron_properties():
    """The coefficients against CD-Bonn's own Table XV, and against a pin.

    The momentum moments are CLOSED FORM -- with
    (2/pi) int dp p^2/[(p^2+a^2)(p^2+b^2)] = 1/(a+b) the norms are double sums
    over the coefficients -- so this tests the twenty numbers with no
    quadrature anywhere to hide an error in.
    """
    w = _l.cdbonn_wave()
    # Table XV: P_D = 4.85 %, eta = 0.0256(4), A_S = 0.8846(9) fm^-1/2.  The
    # parameterisation is a FIT to Machleidt's numerical wave function (his
    # quoted L2 quality is 2.2e-4 in u, 1.1e-4 in w), so it reproduces those
    # to the fit residual and no better: 4.8562 % does not round to 4.85, and
    # A_S misses by 1.3e-4, a seventh of the experimental error it is
    # compared against.  Saying so is the point of the two-sided test below.
    assert abs(w.norm() - 1.0) < 5e-7
    assert abs(w.norm_d() - 0.0485) < 1.5e-4
    assert abs(w.eta() - 0.0256) < 5e-5
    assert abs(w.a_s() - 0.8846) < 2e-4
    # The PINS -- what actually catches a mistyped digit.  Perturbing any one
    # published coefficient by a unit in its last figure moves the norm by
    # 4.6e-9 to 2.4e-7 relative (and P_D by 2e-8 to 1.4e-6 for the D_j), so
    # 1e-9 catches all of them but the last digit of the smallest coefficient
    # C_3, which is worth 2.1e-10 and is not resolvable in double precision.
    assert w.norm() == pytest.approx(0.99999982615384, rel=1e-9)
    assert w.norm_d() == pytest.approx(0.048562073644234, rel=1e-9)
    assert w.norm_s() == pytest.approx(0.95143775250961, rel=1e-9)
    assert w.eta() == pytest.approx(0.0255713786530431, rel=1e-13)
    assert w.a_s() == pytest.approx(0.88472985, rel=1e-15)
    assert w.a_d() == pytest.approx(w.eta() * w.a_s(), rel=1e-15)
    # C_11 and D_9..D_11 are COMPUTED, never typed: the four boundary sums
    # must vanish to Machleidt's "about 15 decimal digits".
    r = list(w.constraint_residuals())
    assert r[0] == 0.0
    assert abs(r[1]) < 1e-12
    assert abs(r[2]) < 1e-10
    assert abs(r[3]) < 1e-8
    assert len(w.c) == 11 and len(w.d) == 11 and len(w.m) == 11
    # m_j = gamma + (j-1) m_0, and m_1 = gamma is what fixes the asymptotics.
    assert w.m[0] == pytest.approx(0.2315380, rel=1e-15)
    assert w.m[10] == pytest.approx(0.2315380 + 10 * 0.9, rel=1e-15)


def test_cdbonn_fdeut_table_is_the_av18_grid_in_the_av18_units():
    """The adapter: same grid, same hbar c convention, CD-Bonn's own e_bind."""
    w = _l.cdbonn_wave()
    t = _l.cdbonn_fdeut_table()
    assert len(t["k_gev"]) == 201            # fdeut.av18's own 0 .. 20 in 0.1
    assert t["k_gev"][0] == 0.0
    # 2.224575 MeV, CD-Bonn's Table XV -- NOT the AV18 file's 2.224574.
    assert t["ebind_gev"] == pytest.approx(2.224575e-3, rel=1e-15)
    # "Both or neither": abscissa fm^-1 -> GeV and ordinate fm^3/2 ->
    # GeV^-3/2, so the GeV table carries the same normalisation as the fm one.
    hbarc = 0.19733
    scale = 1.0 / (hbarc * np.sqrt(hbarc))
    for p in (0.1, 0.7, 2.0, 5.0):
        i = int(round(p / 0.1))
        assert t["k_gev"][i] == pytest.approx(p * hbarc, rel=1e-14)
        assert t["u"][i] == pytest.approx(w.psi_s(p) * scale, rel=1e-12)
        assert t["w"][i] == pytest.approx(w.psi_d(p) * scale, rel=1e-12)
    k = np.asarray(t["k_gev"])
    n = np.trapz(k ** 2 * (np.asarray(t["u"]) ** 2 + np.asarray(t["w"]) ** 2), k)
    assert n == pytest.approx(0.999977985, rel=1e-7)
    # The D wave vanishes at p = 0 -- that IS the sum_j D_j/m_j^2 = 0
    # constraint, seen from momentum space.
    assert abs(t["w"][0]) < 1e-12 * t["u"][0]
    # SIGN: w = -psi_2^a > 0 at small p, the fdeut.av18 / CDKS convention.
    assert t["w"][1] > 0.0 and t["u"][1] > 0.0
    # A finer grid changes the TABLE, not the wave function.
    f = _l.cdbonn_fdeut_table(20.0, 0.02)
    assert len(f["k_gev"]) == 1001
    assert f["u"][50] == pytest.approx(t["u"][10], rel=1e-15)
    with pytest.raises(RuntimeError):
        _l.cdbonn_fdeut_table(20.0, 0.0)


def test_deuteron_wave_option_is_opt_in_and_the_default_did_not_move():
    """`wave` defaults to the AV18 file, and `kCdBonn` ignores `fdeut_path`.

    THE DEFAULT IS A CHOICE and the conservative one: every published A = 2
    gate number in docs/open_items/ was measured on AV18, so moving the
    default would silently move all of them.  This asserts it stayed put.
    """
    o = _l.DeuteronConvolutionB1.Options()
    assert o.wave == _l.DeuteronWaveSource.kFdeutFile
    assert o.fdeut_path.endswith("vmc/deuteron/fdeut.av18")
    assert o.cdbonn_k_max_fm == 20.0
    assert o.cdbonn_dk_fm == 0.1
    # Setting `wave` explicitly to the default must be a NO-OP, bit for bit.
    a, b = _l.DeuteronConvolutionB1.Options(), _l.DeuteronConvolutionB1.Options()
    a.finite_q_delta = b.finite_q_delta = False
    b.wave = _l.DeuteronWaveSource.kFdeutFile
    da, db = _l.DeuteronConvolutionB1(a), _l.DeuteronConvolutionB1(b)
    for x in (0.1, 0.3, 0.5, 0.7):
        assert da.b1(x, 2.5, 0.0) == db.b1(x, 2.5, 0.0)
    # `kCdBonn` is analytic and reads no file at all: a nonsense `fdeut_path`
    # must not stop it (and must still stop `kFdeutFile`).
    c = _l.DeuteronConvolutionB1.Options()
    c.wave = _l.DeuteronWaveSource.kCdBonn
    c.fdeut_path = "/nonexistent/not-a-file"
    c.finite_q_delta = False
    assert _l.DeuteronConvolutionB1(c).options().wave == _l.DeuteronWaveSource.kCdBonn
    bad = _l.DeuteronConvolutionB1.Options()
    bad.fdeut_path = "/nonexistent/not-a-file"
    with pytest.raises(RuntimeError):
        _l.DeuteronConvolutionB1(bad)


def test_cdbonn_moves_the_gate_and_in_the_measured_direction():
    """The swap is worth a factor ~2 at the peak, with a SIGN CHANGE at low x.

    Pinned at CDKS Eq. (17)'s kappa = 1 with the toy F2, which is the cheap
    configuration (the default finite-|q| delta-function rebuilds the
    densities at every x).  The verdict rows -- CDKS's own MSTW2008 LO at
    CDKS's own Eq. (21) delta-function, where CD-Bonn takes the G3b peak
    ratio from 0.8432 to 1.0003 -- are measured in the doctest, checklist
    item 5 of tests/test_b1_nuclear.cpp, and recorded in
    docs/open_items/run_2026-09-03/phase_A_numbers.md section 8.
    """
    def kernel(wave):
        o = _l.DeuteronConvolutionB1.Options()
        o.wave = wave
        o.finite_q_delta = False
        return _l.DeuteronConvolutionB1(o)

    av = kernel(_l.DeuteronWaveSource.kFdeutFile)
    cd = kernel(_l.DeuteronWaveSource.kCdBonn)
    got = {x: (av.b1(x, 2.5, 0.0), cd.b1(x, 2.5, 0.0))
           for x in (0.1, 0.3, 0.5, 0.7)}
    want = {0.1: (-1.9927419369e-04, +1.0386245242e-04),
            0.3: (-9.2802539910e-05, +9.7049151139e-05),
            0.5: (+1.6899772527e-04, +3.2457276551e-04),
            0.7: (+4.0824484064e-04, +5.0442973677e-04)}
    for x, (a, c) in got.items():
        assert a == pytest.approx(want[x][0], rel=1e-9), x
        assert c == pytest.approx(want[x][1], rel=1e-9), x
    # The two structural statements, so the pins above are not the only
    # content: CD-Bonn is LARGER at the large-x peak (it is what closes G3b's
    # deficit) and it has flipped the sign at low x (its S node sits 13 %
    # higher in k, which removes the negative lobe of the S-D interference).
    assert got[0.7][1] > got[0.7][0] > 0.0
    assert got[0.1][0] < 0.0 < got[0.1][1]


# ------------------------------------------------------------------- P6
#
# The b1 NORMALISATION split of 2026-09-03.  The two digitized deuteron camps
# stopped sharing one applied factor: the Miller table is per DEUTERON and is
# halved, the CDKS table is already per NUCLEON and is not.  Arbiter:
# docs/open_items/run_2026-09-03/phase_A_miller_normalisation.md.

def test_the_two_b1_camps_have_their_own_normalisation_constants():
    # 1/A at A = 2 still has exactly ONE home, and Miller's factor IS it --
    # not a second 0.5 typed somewhere else (docs/CONVENTIONS.md).
    assert lg.B1_PER_DEUTERON_TO_PER_NUCLEON == 0.5
    assert lg.B1_MILLER_TABLE_TO_PER_NUCLEON == lg.B1_PER_DEUTERON_TO_PER_NUCLEON
    # CDKS is NOT converted: Eq. (10)'s explicit 1/A, the text under Eq. (16)
    # ("b1 is defined by the one per nucleon"), F1^N = (F1p + F1n)/2.
    assert lg.B1_CDKS_TABLE_TO_PER_NUCLEON == 1.0


def test_cdks_camp_is_the_raw_column_and_miller_is_half_of_his():
    cdks, miller = _l.CdksB1(), _l.MillerB1()
    raw = _l.cdks_b1_raw_per_nucleon()
    for x in (0.05, 0.1, 0.3, 0.5, 0.766, 1.2):
        assert cdks.b1(x, 2.5, 0.0) == pytest.approx(raw.b1(x, 2.5, 0.0),
                                                     rel=1e-15), x
    # the two values tests/test_sf.cpp pins, which DOUBLED on 2026-09-03
    assert cdks.b1(0.3, 2.5, 0.0) == pytest.approx(-5.632397950173448e-04,
                                                   rel=1e-12)
    assert cdks.b1(0.05, 2.5, 0.0) == pytest.approx(+1.1703432536002871e-04,
                                                    rel=1e-12)
    # Miller's is still halved, and the DEFAULT b1 model is his, so this is
    # the assertion that says the default did not move.
    assert (miller.b1(0.012, 2.5, 0.0) * 2.0
            == pytest.approx(0.11429317074113018, rel=1e-12))


# ------------------------------------------------------------------- P7
#
# The r-default decision of 2026-09-03 (OPEN_ITEMS_SOLUTIONS.md section 10,
# condition 4): `Li6ConvolutionOptions` keeps `r_sigma_lt` and does NOT adopt
# the A = 2 gate's `r1998`, because the shipped observable is a ratio whose
# denominator carries `InclusiveKernel`'s R.

def test_li6_convolution_r_default_is_r_sigma_lt_and_the_swap_is_measured():
    def kernel(r_func=None):
        o = _l.Li6ConvolutionOptions()
        if r_func is not None:
            o.r_func = r_func
        return _l.Li6ConvolutionB1(o)

    lt = kernel()                                          # the default
    lt_explicit = kernel(lambda x, q2: _l.r_sigma_lt(x, q2))
    r98 = kernel(lambda x, q2: _l.r1998(x, q2))

    for x in (0.05, 0.10, 0.20, 0.30, 0.50):
        # the null hook IS r_sigma_lt, bit for bit
        assert lt.b1(x, 2.5, 0.0) == lt_explicit.b1(x, 2.5, 0.0), x
        # term (1) does not see R at all -- it carries b1_d from the injected
        # TensorSF (the raw digitized CDKS column), which has no R in it.
        # This is the assertion that explains why the swap is NOT a (1+R)
        # prefactor, and it must stay bit-for-bit.
        assert (lt.b1_embedded_s(x, 2.5, 0.0)
                == r98.b1_embedded_s(x, 2.5, 0.0)), x

    # the measured cost of the decision, pinned so it cannot drift silently
    want = {0.05: (+3.650716e-06, +3.449387e-06),
            0.10: (-2.796612e-06, -3.484857e-06),
            0.20: (-2.120522e-05, -2.235704e-05),
            0.30: (-4.173248e-05, -4.284379e-05),
            0.50: (+5.173934e-05, +5.075777e-05)}
    for x, (a, b) in want.items():
        assert x * lt.b1(x, 2.5, 0.0) == pytest.approx(a, rel=1e-5), x
        assert x * r98.b1(x, 2.5, 0.0) == pytest.approx(b, rel=1e-5), x
    # the orbital term is where the whole effect lives, and it is far larger
    # than any (1 + R) prefactor: x0.54 to x0.94.
    for x in (0.10, 0.20, 0.30):
        ratio = (r98.b1_alpha_d_dwave_d(x, 2.5, 0.0)
                 / lt.b1_alpha_d_dwave_d(x, 2.5, 0.0))
        assert 0.50 < ratio < 0.95, (x, ratio)


# ------------------------------------------------------------------- P8
# `--b1-unpol`: THE GATE-PASSING CONFIGURATION HAS TO BE EMITTABLE.
#
# The A = 2 gate passes at G3b = 0.843 with CDKS's own MSTW2008 LO and at
# 0.440 on the library ToyF2, so the pass is a statement about a
# CONFIGURATION, not about a build.  Before `B1UnpolSource` the run surface
# could not produce it: `default_inclusive_kernel` hard-wired `ToyF2` into
# `Li6ConvolutionOptions::unpol`, and the one other route -- a hand-built
# `PipelineConfig.kernel` -- is refused together with a non-Miller `b1_model`
# by `validate()`.  Every 6Li b1 number the generator could emit was made at
# 0.440, on a curve up to a factor 1.85 away from the one the gate cleared.
#
# THE KERNEL/FLAG GUARD IS NOT LOOSENED to fix that, and P8 pins it: the
# selector works THROUGH the pipeline's own kernel construction, so the
# passing configuration no longer needs a hand-built kernel to be expressed,
# and a hand-built kernel plus the flag still throws.

#: mstw / toy on `Li6ConvolutionB1.b1(x, 2.5)`, measured 2026-09-04.  Up to a
#: factor 1.85 and NOT monotone -- a shape change, not a normalisation.
MSTW_OVER_TOY = {0.10: 1.847766, 0.30: 1.275961, 0.50: 0.816971}
CT18_OVER_TOY = {0.10: 2.221302, 0.30: 1.238833, 0.50: 1.045870}

needs_pythia = pytest.mark.skipif(not _l.HAVE_PYTHIA8,
                                  reason="MstwSF needs the PYTHIA 8 tier")
needs_lhapdf = pytest.mark.skipif(not _l.HAVE_LHAPDF,
                                  reason="LhapdfSF needs the LHAPDF tier")


def _unpol_cfg(unpol=None, model="li6-convolution", coarse=True, **kw):
    """A cheap li6-convolution config; the coarse grid is a SPEED knob only."""
    cfg = _cfg(model, **kw) if unpol is None else \
        lg.make_config(isotope="6Li", channel="inclusive", config=1, events=N,
                       seed=SEED, b1_model=model, b1_unpol=unpol, **kw)
    if coarse:
        g = cfg.grid            # SamplerGrid comes back BY VALUE
        g.nx, g.nq2 = 12, 6
        cfg.grid = g
    sc = cfg.scenario
    sc.x_max = X_MAX
    cfg.scenario = sc
    return cfg


def test_b1_unpol_name_table():
    assert set(lg.B1_UNPOL) == {"toy", "mstw", "ct18nlo"}
    assert lg.B1_UNPOL["toy"] == _l.B1UnpolSource.Toy
    assert lg.B1_UNPOL["mstw"] == _l.B1UnpolSource.Mstw
    assert lg.B1_UNPOL["ct18nlo"] == _l.B1UnpolSource.Ct18Nlo
    assert "B1_UNPOL" in lg.__all__
    for key, val in lg.B1_UNPOL.items():
        assert _l.b1_unpol_name(val) == key
    # Custom is deliberately NOT a CLI name: it names an OBJECT, so it is set
    # by assigning `config.b1_unpol_sf`, never by a flag.
    assert _l.b1_unpol_name(_l.B1UnpolSource.Custom) == "custom"
    assert "custom" not in lg.B1_UNPOL


def test_the_default_did_not_move():
    """THE DEFAULT IS BIT FOR BIT, which is the whole licence for this flag.

    Three clauses, each on a different surface:
      * the config default is Toy with an EMPTY slot, so the pipeline's call
        into `default_inclusive_kernel` is the pre-flag four-argument call;
      * the li6-convolution kernel built that way is bit-identical to one
        built with `b1_unpol` spelled out as None, and its b1 IS
        `Li6ConvolutionB1()`'s -- the default-options object this change
        never touched;
      * an end-to-end default run and an explicit `--b1-unpol toy` run agree
        in every event weight, not merely in the cross section.
    `test_default_kernel_is_the_miller_kernel` above pins the SHIPPED default
    (miller) against validation/reference/b1_default_li6.json at rtol 1e-12
    and is the fourth clause; it is unchanged by this flag.
    """
    cfg = lg.make_config(events=N)
    assert cfg.b1_unpol == _l.B1UnpolSource.Toy
    assert cfg.b1_unpol_sf is None

    a = _l.default_inclusive_kernel(_l.li6(), _l.B1Model.Li6Convolution)
    b = _l.default_inclusive_kernel(_l.li6(), _l.B1Model.Li6Convolution,
                                    1.0, 1.0, None)
    conv = lg.Li6ConvolutionB1()          # default options: unpol = null
    for x in (0.10, 0.30, 0.50):
        ta, tb = a.tables(x, 2.5), b.tables(x, 2.5)
        assert ta.b1 == tb.b1, x          # bit for bit, not merely close
        assert ta.f1 == tb.f1, x
        assert ta.delta == tb.delta, x
        assert ta.b1 == conv.b1(x, 2.5, ta.f1), x

    plan = lg.tensor_thirds_plan(0.0, 0.6)
    p_implicit = _l.Pipeline(_unpol_cfg(), plan)
    p_explicit = _l.Pipeline(_unpol_cfg("toy"), plan)
    assert np.array_equal(np.asarray(p_implicit.dis_sampler.cell_xsec_pb),
                          np.asarray(p_explicit.dis_sampler.cell_xsec_pb))
    assert np.array_equal(np.asarray(p_implicit.sigma_per_category_pb()),
                          np.asarray(p_explicit.sigma_per_category_pb()))
    d_implicit = p_implicit.generate(N, False)
    d_explicit = p_explicit.generate(N, False)
    for col in ("x", "q2", "y", "weight"):
        assert np.array_equal(np.asarray(d_implicit[col]),
                              np.asarray(d_explicit[col])), col


def test_make_config_and_binding_keep_the_two_fields_in_step():
    """The enum is the PROVENANCE and the object is the realisation.

    They are one choice (the `optics_choice` / `optics` arrangement), so
    neither surface lets them drift: `set_b1_unpol` builds the object for a
    named backend, and assigning an object names it `Custom`.
    """
    cfg = lg.make_config(events=N, b1_model="li6-convolution", b1_unpol="toy")
    assert cfg.b1_unpol == _l.B1UnpolSource.Toy and cfg.b1_unpol_sf is None

    with pytest.raises(ValueError):
        lg.make_config(events=N, b1_unpol="nope")

    cfg.b1_unpol_sf = lg.ToyF2()             # an object IS Custom
    assert cfg.b1_unpol == _l.B1UnpolSource.Custom
    cfg.b1_unpol_sf = None                   # ... and clearing it is Toy
    assert cfg.b1_unpol == _l.B1UnpolSource.Toy

    # Custom cannot be asked for by name -- it names an object this function
    # cannot build.
    with pytest.raises(RuntimeError):
        _l.set_b1_unpol(cfg, _l.B1UnpolSource.Custom)


def test_validate_refuses_a_named_backend_with_an_empty_slot():
    """A NAMED BACKEND IS NEVER SILENTLY REPLACED BY ToyF2.

    That fallback would reproduce exactly the defect this flag exists to fix
    -- a run labelled `mstw` whose numbers are the `toy` ones.  The core
    library links neither the PYTHIA nor the LHAPDF tier, so it cannot build
    the backend itself; it throws and says which tier is missing.
    """
    for src in (_l.B1UnpolSource.Mstw, _l.B1UnpolSource.Ct18Nlo,
                _l.B1UnpolSource.Custom):
        cfg = _cfg("li6-convolution")
        cfg.b1_unpol = src                   # the enum ALONE, no object
        with pytest.raises(RuntimeError) as e:
            cfg.validate()
        assert _l.b1_unpol_name(src) in str(e.value)
        cfg.b1_unpol_sf = lg.ToyF2()         # ... legal once attached
        cfg.b1_unpol = src                   # (the setter above said Custom)
        cfg.validate()

    # `toy` with an object attached is the same contradiction the other way
    cfg = _cfg("li6-convolution")
    cfg.b1_unpol_sf = lg.ToyF2()
    cfg.b1_unpol = _l.B1UnpolSource.Toy
    with pytest.raises(RuntimeError):
        cfg.validate()


@pytest.mark.parametrize("model", ["miller", "cdks"])
def test_validate_refuses_the_unpol_flag_where_it_does_not_run(model):
    """A KNOB THAT DID NOT RUN MAY NOT BE RECORDED AS IF IT HAD.

    Only the li6-convolution branch reads `Li6ConvolutionOptions::unpol` --
    miller is a ratio model with no F1 of its own and cdks carries the
    digitized column -- but `meta["b1_unpol"]` is written unconditionally, so
    on either of the other two the flag would reach the file and not the rate.
    """
    cfg = _cfg(model)
    cfg.b1_unpol_sf = lg.ToyF2()             # -> Custom
    with pytest.raises(RuntimeError) as e:
        cfg.validate()
    assert "li6-convolution" in str(e.value)


def test_the_caller_supplied_kernel_guard_is_unchanged():
    """The guard is NOT loosened to let the passing configuration through.

    A caller-supplied kernel wins over every b1 flag and would make the
    metadata say something the numbers do not, so the two stay refused
    together.  What changed is that the configuration no longer NEEDS a
    hand-built kernel: `--b1-model li6-convolution --b1-unpol mstw` builds it
    through `default_inclusive_kernel`.
    """
    cfg = _cfg("li6-convolution")
    cfg.b1_unpol_sf = lg.ToyF2()
    cfg.kernel = _l.default_inclusive_kernel(_l.li6())
    with pytest.raises(RuntimeError):
        cfg.validate()
    # and the flagless default with a hand-built kernel is still legal
    cfg.b1_unpol_sf = None
    cfg.b1_model = _l.B1Model.Miller
    cfg.validate()


def test_cli_b1_unpol_flag():
    assert "b1_unpol" in cli.DEFAULTS
    assert cli.DEFAULTS["b1_unpol"] == "toy"
    assert cli.resolve([])["b1_unpol"] == "toy"
    assert cli.resolve(["--b1-unpol", "mstw"])["b1_unpol"] == "mstw"
    help_text = cli.build_parser().format_help()
    assert "--b1-unpol" in help_text
    for name in sorted(lg.B1_UNPOL):
        assert name in help_text, name


@needs_pythia
def test_mstw_moves_b1_through_the_pipelines_own_kernel():
    """The measured gap on the SHIPPED observable, at the run surface.

    `default_inclusive_kernel` is what the `Pipeline` calls when no kernel is
    supplied, so these are the numbers a `--b1-model li6-convolution
    --b1-unpol mstw` run actually generates with.
    """
    toy = _l.default_inclusive_kernel(_l.li6(), _l.B1Model.Li6Convolution)
    mstw = _l.default_inclusive_kernel(_l.li6(), _l.B1Model.Li6Convolution,
                                       1.0, 1.0, _l.MstwSF())
    for x, want in MSTW_OVER_TOY.items():
        tt, tm = toy.tables(x, 2.5), mstw.tables(x, 2.5)
        assert tt.b1 != 0.0
        assert tm.b1 / tt.b1 == pytest.approx(want, rel=1e-5), x
        # F1 DOES NOT MOVE: the kernel's own f2_source is ToyF2 on every
        # setting of this flag, so F1 -- and with it the spin-blind cell
        # cross section and the D_phi denominator of the tensor weight -- is
        # bit-identical.  The price is that numerator and denominator no
        # longer share one object; see B1UnpolSource in pipeline.hpp.
        assert tm.f1 == tt.f1, x
        assert tm.delta == tt.delta, x


@needs_lhapdf
def test_ct18nlo_is_selectable_too():
    toy = _l.default_inclusive_kernel(_l.li6(), _l.B1Model.Li6Convolution)
    ct18 = _l.default_inclusive_kernel(_l.li6(), _l.B1Model.Li6Convolution,
                                       1.0, 1.0, _l.LhapdfSF("CT18NLO", 0))
    for x, want in CT18_OVER_TOY.items():
        assert (ct18.tables(x, 2.5).b1 / toy.tables(x, 2.5).b1
                == pytest.approx(want, rel=1e-5)), x


@needs_pythia
def test_end_to_end_the_tensor_observable_moves_and_the_rate_does_not():
    """A 6Li inclusive run under each selector, through the CLI's own config.

    The SPIN-BLIND cell cross sections are bit-identical -- b1 enters the rate
    only through the per-state tensor shift, and this flag moves b1 alone --
    while the tensor-thirds categories separate by a different amount.  That
    is the designed scope of `--b1-unpol`, and the same invariant P3 uses for
    the band.

    (`sigma_per_category_pb` and `sigma_pb` DO move: they carry the tensor
    shift.  The last-ULP agreement of `sigma_pb` on this coarse grid is an
    accident of summation and is not asserted.)
    """
    plan = lg.tensor_thirds_plan(0.0, 0.6)
    out = {}
    for u in ("toy", "mstw"):
        p = _l.Pipeline(_unpol_cfg(u), plan)
        out[u] = (np.asarray(p.dis_sampler.cell_xsec_pb),
                  np.asarray(p.sigma_per_category_pb()),
                  p.sigma_pb(), p.generate(1, True)["meta"]["b1_unpol"])

    assert out["toy"][3] == "toy" and out["mstw"][3] == "mstw"
    assert np.array_equal(out["toy"][0], out["mstw"][0])   # the rate: bit for bit
    # the tensor signal: it MOVED, and by 16 %
    def azz(row):
        return (row[1][2] - row[1][0]) / row[2]
    a_toy, a_mstw = azz(out["toy"]), azz(out["mstw"])
    assert a_toy == pytest.approx(2.4281032e-05, rel=1e-5)
    assert a_mstw == pytest.approx(2.0423864e-05, rel=1e-5)
    assert a_mstw / a_toy == pytest.approx(0.841145, rel=1e-4)


def test_a_missing_tier_or_grid_fails_loudly_and_names_what_is_missing():
    """Point of honesty: NEVER a quiet downgrade to the toy.

    The toy's A = 2 gate ratio is 0.440 against MSTW's 0.843, so a silent
    fallback would put the toy's numbers under the label the gate cleared --
    which is the exact defect this flag exists to remove.  Three loud
    failures, one per layer:
      * the CLI refuses the flag when the tier is absent (both branches are
        exercised here by passing the availability flags in, so the test does
        not need a build that lacks a tier);
      * `set_b1_unpol` refuses the same thing at the binding;
      * with the tier present, a missing grid FILE throws out of `MstwSF`
        naming the file and the directory it was looked for in.
    """
    with pytest.raises(SystemExit) as e:
        cli.require_b1_unpol_tier("mstw", have_pythia=False)
    assert "PYTHIA" in str(e.value) and "mstw" in str(e.value)
    with pytest.raises(SystemExit) as e:
        cli.require_b1_unpol_tier("ct18nlo", have_lhapdf=False)
    assert "LHAPDF" in str(e.value)
    # ... and it is a no-op where the tier is there, and for the toy always
    cli.require_b1_unpol_tier("toy", have_pythia=False, have_lhapdf=False)
    cli.require_b1_unpol_tier("mstw", have_pythia=True)
    cli.require_b1_unpol_tier("ct18nlo", have_lhapdf=True)

    cfg = _cfg("li6-convolution")
    if _l.HAVE_PYTHIA8:
        with pytest.raises(RuntimeError) as e:
            _l.MstwSF(3, "/nonexistent-pdfdata-dir")
        assert "mstw2008lo.00.dat" in str(e.value)
        assert "/nonexistent-pdfdata-dir" in str(e.value)
    else:
        with pytest.raises(RuntimeError) as e:
            _l.set_b1_unpol(cfg, _l.B1UnpolSource.Mstw)
        assert "PYTHIA" in str(e.value)
    if not _l.HAVE_LHAPDF:
        with pytest.raises(RuntimeError) as e:
            _l.set_b1_unpol(cfg, _l.B1UnpolSource.Ct18Nlo)
        assert "LHAPDF" in str(e.value)


def test_meta_records_which_unpol_backend_made_the_file():
    """Without this key a toy npz and an mstw npz are indistinguishable.

    They differ by up to a factor 1.85 in b1 and by NOTHING else in the file
    -- the cross section, the kinematics and the seed are identical -- so the
    key is the only thing that tells them apart.
    """
    plan = lg.tensor_thirds_plan(0.0, 0.6)
    meta = _l.Pipeline(_unpol_cfg("toy"), plan).generate(1, True)["meta"]
    assert meta["b1_unpol"] == "toy"
    # the shipped default says what it is too, so no npz is anonymous
    meta = _pipeline().generate(1, True)["meta"]
    assert meta["b1_unpol"] == "toy"
    # ... and a hand-built kernel is not labelled with a flag that did not run
    cfg = _cfg()
    cfg.kernel = _l.default_inclusive_kernel(_l.li6())
    meta = _l.Pipeline(cfg, plan).generate(1, True)["meta"]
    assert meta["b1_unpol"] == "caller-supplied kernel"


# ------------------------------------------------------------------- P9
# The gate's CLAIM, as opposed to the gate's wiring: three places where the
# tree used to say more than it had measured.  Each one is a pure function or
# a source grep on purpose -- none of them needs a pipeline, so all three run
# in milliseconds and there is no excuse to drop them.


@pytest.mark.parametrize("name, covered", [
    ("mstw", True), ("toy", False), ("ct18nlo", False),
])
def test_the_run_banner_says_which_side_of_the_gate_this_run_is_on(name,
                                                                   covered):
    """A "PASSED" banner over a run the pass does not cover is an overclaim.

    The A = 2 verdict is a statement about ONE unpolarised nucleon input --
    MSTW2008 LO at CDKS Eq. (21) -- and the shipped default is not it.  The
    banner therefore has to say, per run, whether THIS run is in that
    configuration, and it must not call ct18nlo covered either: 0.719 is
    inside the window but CT18NLO is the retired stand-in, not the row the
    lift was read off.
    """
    lines = cli.b1_gate_lines(name)
    assert len(lines) == 2
    head, verdict = lines
    # the standing claim always carries its configuration and the default's
    # own number, so the line can never be quoted as an unconditional pass
    assert "MSTW2008 LO" in head and "Eq. (21)" in head
    assert "0.440" in head and "OUTSIDE" in head
    # ... and never as three-digit agreement with CDKS
    assert "1.000338" in head and "NOT three-digit agreement" in head
    assert "1.000 " not in head
    if covered:
        assert verdict.startswith("^ THIS RUN IS in that configuration")
        assert "NOT covered" not in verdict
    else:
        assert verdict.startswith("^ THIS RUN IS NOT in that configuration")
        assert name in verdict
        assert "NOT covered by the 2026-09-03 lift" in verdict
        assert "--b1-unpol mstw" in verdict


def test_the_validation_report_footer_is_conditional_on_the_verdict(repo_root):
    """`b1_li6_table.py` printed the ban lift even when it measured nothing.

    The `verdict is None` branch prints "MSTW2008 LO is not available in this
    build, so G3b's VERDICT row was not measured here" -- and then the footer
    announced the lift anyway, into a file that outlives the run.
    """
    sys.path.insert(0, os.path.join(repo_root, "validation"))
    try:
        import b1_li6_table
    finally:
        sys.path.pop(0)

    missing = b1_li6_table.gate_footer(None, 0.440001)
    assert "DID NOT MEASURE THE VERDICT ROW" in missing
    assert "WAS LIFTED" not in missing
    assert "0.4400" in missing               # says what it DID measure

    measured = b1_li6_table.gate_footer(0.843243, 0.440001)
    assert "THE 6Li PUBLICATION BAN WAS LIFTED ON 2026-09-03" in measured
    assert "NOT for the ToyF2 default row" in measured

    # THE THIRD BRANCH, and the one the first fix still got wrong: a verdict
    # row that was MEASURED and FAILED.  Keying the footer on "was it
    # measured?" printed "OUTSIDE [0.5, 2]: G3b FAILS" in the body and the
    # ban lift in the footer of the same report.  It is keyed on the ratio.
    failing = b1_li6_table.gate_footer(0.44, 0.440001)
    assert "G3b FAILS at 0.440000" in failing
    assert "WAS LIFTED" not in failing
    assert "contradicts it" in failing

    # both edges of the acceptance window are inside, neither outside
    assert "WAS LIFTED" in b1_li6_table.gate_footer(0.5, 0.440001)
    assert "WAS LIFTED" in b1_li6_table.gate_footer(2.0, 0.440001)
    assert "WAS LIFTED" not in b1_li6_table.gate_footer(2.0001, 0.440001)

    # the band does not depend on the branch: it never came from this gate
    for text in (missing, measured, failing):
        assert "MANDATORY +-100 % BAND" in text


def test_no_shipped_python_source_still_says_the_gate_fails(repo_root):
    """A1: the `B1_MODELS` docstring said "NOT passed" a day after the lift.

    It ships in the installed package and describes the same enum as
    `python/bindings.cpp`, so the two shipped descriptions of one object said
    opposite things.  Nothing caught it because every ban-lift pass grepped
    the BAN wording ("may be published"), and this site carries none.  This
    test greps the FAIL wording instead, over everything the wheel installs.
    """
    pkg = os.path.join(repo_root, "python", "lipolgen")
    stale = ("magnitude gate is\n#: NOT passed", "gate is NOT passed",
             "has NOT passed", "fails its magnitude", "fails on magnitude",
             "gate is not passed")
    hits = []
    for fname in sorted(os.listdir(pkg)):
        if not fname.endswith(".py"):
            continue
        with open(os.path.join(pkg, fname)) as fh:
            text = fh.read()
        flat = " ".join(text.split())
        for needle in stale:
            if " ".join(needle.split()) in flat:
                hits.append((fname, needle))
    assert not hits, (
        "shipped Python source still asserts the A = 2 gate fails: %r.  The "
        "gate passes for the MSTW2008 LO input at CDKS Eq. (21) (G3b 0.843243) "
        "and is 0.440 on the toy default -- say the configuration, do not say "
        "it failed." % (hits,))
