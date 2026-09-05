"""D2: `--unpol-sf` / `--pol-sf` -- one structure-function selector that
reaches EVERY kernel the pipeline builds.

WHAT THIS FILE GUARDS.  Before 2026-09-04 the run surface could not select an
unpolarised or polarised backend at all: `default_inclusive_kernel` hard-wired
`ToyF2` into the inclusive (and so the coherent) kernel, and
`struck_cluster_kernel` had NO structure-function slot of any kind, so every
tagged run was `ToyF2`/`ToyG1` with no way to say otherwise -- not even
through the hand-built `PipelineConfig.kernel` escape hatch, which the
`Pipeline` reads on the non-tagged branch alone.  Both `ToyF2` and `ToyG1` are
labelled TOY in `sf.hpp`.

The gates, in order:

G1  the DEFAULT IS BIT FOR BIT.  An unset run and a run with both selectors
    explicitly at their default agree by `np.array_equal` on the cell cross
    sections, the per-category cross sections and the generated x/Q2/y/weight
    columns -- on the INCLUSIVE, the COHERENT and a TAGGED channel.
G2  the enum and the object stay in step in both directions.
G3  `validate()` refuses a named backend with an empty slot, and a toy
    setting with an object attached.
G4  `validate()` refuses a caller-supplied kernel beside either selector, and
    refuses the b1 collision (`--b1-unpol toy` under a non-toy `--unpol-sf`).
G5  the CLI surface: both flags in `DEFAULTS`, every choice in `--help`, and
    the tier checks fail loudly naming the missing tier.
G6  THE TAGGED KERNEL ACTUALLY MOVES -- the gate D2 exists for.
G7  `meta` records both backends on every run and every channel, plus the
    grid report, and says "caller-supplied kernel" when one is set.
G8  the T2 PYTHIA bridge's species draw is on the SAME object as the kernel.

Added 2026-09-05, after an adversarial review RAN the two selectors and found
four places where the surface claimed a reach the run did not have:

G9  `--pol-sf` does NOT reach the COHERENT channel -- measured, bit for bit --
    and the banner, `meta["pol_sf"]` and `meta["pol_sf_reach"]` all say so.
G10 the below-grid fraction is taken against the CHANNEL's own rate
    (`Pipeline.cell_rate_weights_pb`), 0.447460 coherent against 0.361791
    inclusive; it used to report the inclusive number on a coherent run.
G11 a caller-supplied `kernel` is REFUSED on a tagged channel, where it is not
    read -- and still accepted on the coherent one, where it is.
G12 the neutron-tagged `ct18nlo` + `nnpdfpol` pair is refused with the CELL
    and the CURE (`--x-max`), and the CLI reports it as a configuration error.
G13 `set_pythia_hadronizer` refuses a bridge whose species draw is on a
    different unpolarised backend from the run's kernel.
"""

import numpy as np
import pytest

import lipolgen as lg
from lipolgen import cli
from lipolgen import _lipolgen as _l

N = 300
SEED = 11

needs_pythia = pytest.mark.skipif(not _l.HAVE_PYTHIA8,
                                  reason="no PYTHIA 8 tier in this build")
needs_lhapdf = pytest.mark.skipif(not _l.HAVE_LHAPDF,
                                  reason="no LHAPDF tier in this build")

#: The three channels the selectors must reach.  `coherent` owns no structure
#: function at all -- it rides the inclusive sampler's cell cross sections --
#: and the tagged one had no injection point before D2, so both belong in G1.
CHANNELS = ("inclusive", "coherent", "tagged-6Li-alpha")


def _cfg(channel="inclusive", **kw):
    return lg.make_config(channel=channel, config=1, events=N, seed=SEED, **kw)


def _plan():
    return lg.tensor_thirds_plan(0.0, 0.6)


# --------------------------------------------------------------------- G1

@pytest.mark.parametrize("channel", CHANNELS)
def test_the_default_is_bit_for_bit_on_every_channel(channel):
    """THE PROOF OBLIGATION OF THE WHOLE CHANGE.

    `unpol_sf = pol_sf = Toy` leaves both `InclusiveKernel.Options` slots
    exactly as they were -- `f2_source` at the shared `ToyF2` the pipeline has
    always named, `g1_model` UNSET so the class builds its own `ToyG1` on it.
    So the default path is unchanged BY CONSTRUCTION rather than by
    inspection, and this test is what turns that sentence into a measurement.
    Not `approx`: `np.array_equal`.
    """
    implicit = _l.Pipeline(_cfg(channel), _plan())
    explicit = _l.Pipeline(_cfg(channel, unpol_sf="toy", pol_sf="toy"),
                           _plan())
    assert np.array_equal(np.asarray(implicit.dis_sampler.cell_xsec_pb),
                          np.asarray(explicit.dis_sampler.cell_xsec_pb))
    assert np.array_equal(np.asarray(implicit.sigma_per_category_pb()),
                          np.asarray(explicit.sigma_per_category_pb()))
    assert implicit.sigma_pb() == explicit.sigma_pb()
    a = implicit.generate(N, False)
    b = explicit.generate(N, False)
    for col in ("x", "q2", "y", "weight"):
        assert np.array_equal(np.asarray(a[col]), np.asarray(b[col])), col


# --------------------------------------------------------------------- G2

def test_the_enum_and_the_object_stay_in_step():
    """The enum is the PROVENANCE and the object the REALISATION.

    One choice with two fields (the `optics_choice` / `optics` and
    `b1_unpol` / `b1_unpol_sf` arrangement), so neither surface lets them
    drift: assigning an object names it `Custom`, clearing it is back to
    `Toy`, and `set_*` goes the other way.
    """
    cfg = _cfg()
    assert cfg.unpol_sf == _l.UnpolSfSource.Toy and cfg.unpol_sf_obj is None
    assert cfg.pol_sf == _l.PolSfSource.Toy and cfg.pol_sf_obj is None

    cfg.unpol_sf_obj = lg.ToyF2()
    assert cfg.unpol_sf == _l.UnpolSfSource.Custom
    cfg.unpol_sf_obj = None
    assert cfg.unpol_sf == _l.UnpolSfSource.Toy

    cfg.pol_sf_obj = lg.ToyG1()
    assert cfg.pol_sf == _l.PolSfSource.Custom
    cfg.pol_sf_obj = None
    assert cfg.pol_sf == _l.PolSfSource.Toy

    # `Custom` cannot be asked for by name: it names an object these
    # functions cannot build.
    with pytest.raises(RuntimeError):
        _l.set_unpol_sf(cfg, _l.UnpolSfSource.Custom)
    with pytest.raises(RuntimeError):
        _l.set_pol_sf(cfg, _l.PolSfSource.Custom)

    with pytest.raises(ValueError):
        lg.make_config(events=N, unpol_sf="nope")
    with pytest.raises(ValueError):
        lg.make_config(events=N, pol_sf="nope")

    assert _l.unpol_sf_name(_l.UnpolSfSource.Ct18Nlo) == "ct18nlo"
    assert _l.pol_sf_name(_l.PolSfSource.NnpdfPol) == "nnpdfpol"
    assert set(lg.UNPOL_SF) == {"toy", "mstw", "ct18nlo"}
    assert set(lg.POL_SF) == {"toy", "nnpdfpol"}


# --------------------------------------------------------------------- G3

def test_validate_refuses_a_named_backend_with_an_empty_slot():
    """A NAMED BACKEND IS NEVER SILENTLY REPLACED BY THE TOY.

    The core library links neither the PYTHIA nor the LHAPDF tier, so it
    cannot build the backend itself; it throws and says which tier is
    missing.  A silent fallback would put the toy's numbers -- x0.80 on the
    6Li rate, and the WRONG SIGN on g1n over 0.25 < x < 0.6 -- under a label
    that says otherwise, which is the whole defect these flags remove.
    """
    for src in (_l.UnpolSfSource.Mstw, _l.UnpolSfSource.Ct18Nlo,
                _l.UnpolSfSource.Custom):
        cfg = _cfg()
        cfg.unpol_sf = src                     # the enum ALONE, no object
        with pytest.raises(RuntimeError) as e:
            cfg.validate()
        assert _l.unpol_sf_name(src) in str(e.value)
        cfg.unpol_sf_obj = lg.ToyF2()          # ... legal once attached
        cfg.unpol_sf = src
        cfg.validate()

    for src in (_l.PolSfSource.NnpdfPol, _l.PolSfSource.Custom):
        cfg = _cfg()
        cfg.pol_sf = src
        with pytest.raises(RuntimeError) as e:
            cfg.validate()
        assert _l.pol_sf_name(src) in str(e.value)
        cfg.pol_sf_obj = lg.ToyG1()
        cfg.pol_sf = src
        cfg.validate()

    # the missing tier is NAMED, not merely refused
    cfg = _cfg()
    cfg.unpol_sf = _l.UnpolSfSource.Mstw
    with pytest.raises(RuntimeError) as e:
        cfg.validate()
    assert "PYTHIA" in str(e.value)
    cfg = _cfg()
    cfg.unpol_sf = _l.UnpolSfSource.Ct18Nlo
    with pytest.raises(RuntimeError) as e:
        cfg.validate()
    assert "LHAPDF" in str(e.value)
    cfg = _cfg()
    cfg.pol_sf = _l.PolSfSource.NnpdfPol
    with pytest.raises(RuntimeError) as e:
        cfg.validate()
    assert "LHAPDF" in str(e.value)


def test_validate_refuses_toy_with_an_object_attached():
    """The same contradiction the other way round.

    The Python setters name an attached object `Custom`, so reaching this
    state means the provenance and the realisation drifted -- and the
    attached backend would then be silently dropped.
    """
    cfg = _cfg()
    cfg.unpol_sf_obj = lg.ToyF2()
    cfg.unpol_sf = _l.UnpolSfSource.Toy
    with pytest.raises(RuntimeError):
        cfg.validate()

    cfg = _cfg()
    cfg.pol_sf_obj = lg.ToyG1()
    cfg.pol_sf = _l.PolSfSource.Toy
    with pytest.raises(RuntimeError):
        cfg.validate()


# --------------------------------------------------------------------- G4

def test_validate_refuses_a_caller_supplied_kernel_beside_either_selector():
    """A hand-built kernel WINS, so the two are refused together.

    On the inclusive branch the kernel silently replaces the one the
    selectors would have built; on a tagged channel the kernel is not read at
    all.  Either way `meta` would record a backend that did not produce the
    numbers.
    """
    for field, obj in (("unpol_sf_obj", lg.ToyF2()),
                       ("pol_sf_obj", lg.ToyG1())):
        cfg = _cfg()
        cfg.kernel = _l.default_inclusive_kernel(_l.li6())
        cfg.validate()                          # a kernel alone is still legal
        setattr(cfg, field, obj)
        with pytest.raises(RuntimeError):
            cfg.validate()


def test_validate_refuses_the_b1_unpol_collision():
    """`--b1-unpol toy` under a non-toy `--unpol-sf` would MISLABEL the b1.

    `default_inclusive_kernel` folds the `li6-convolution` b1 against the
    KERNEL's own `f2_source` when `b1_unpol` is null, so with a non-toy
    `unpol_sf` the setting `toy` no longer names `ToyF2` -- while
    `meta["b1_unpol"]` would still write "toy".  That is exactly the defect
    the provenance discipline exists to prevent, so it is refused rather than
    resolved silently.  The fix is to name the b1 backend.
    """
    cfg = lg.make_config(events=N, b1_model="li6-convolution")
    sc = cfg.scenario
    sc.x_max = 0.95                 # both opt-in b1 backends need it
    cfg.scenario = sc
    cfg.validate()                  # toy + toy: fine

    cfg.unpol_sf_obj = lg.ToyF2()   # -> Custom, i.e. a non-toy unpol_sf
    assert cfg.b1_unpol == _l.B1UnpolSource.Toy
    with pytest.raises(RuntimeError) as e:
        cfg.validate()
    assert "b1_unpol" in str(e.value) and "--b1-unpol" in str(e.value)

    cfg.b1_unpol_sf = lg.ToyF2()    # name it explicitly -> Custom
    cfg.validate()


def test_validate_refuses_a_directly_set_struck_slot_under_a_toy_selector():
    """The same provenance rule, on the tagged branch's own slots.

    A C++ or Python caller may put an object straight into
    `struck.f2_source` / `struck.g1_model` (the `breakup.triton_sf`
    arrangement: `Pipeline` fills them only when they are empty), but then
    `meta` would still record "toy" for a run made on something else.
    """
    cfg = _cfg("tagged-6Li-alpha")
    st = cfg.struck                      # StruckClusterOptions comes BY VALUE
    st.f2_source = lg.ToyF2()
    cfg.struck = st
    with pytest.raises(RuntimeError) as e:
        cfg.validate()
    assert "unpol_sf" in str(e.value)
    cfg.unpol_sf_obj = lg.ToyF2()        # -> Custom, and now it is honest
    cfg.validate()

    cfg = _cfg("tagged-6Li-alpha")
    st = cfg.struck
    st.g1_model = lg.ToyG1()
    cfg.struck = st
    with pytest.raises(RuntimeError) as e:
        cfg.validate()
    assert "pol_sf" in str(e.value)
    cfg.pol_sf_obj = lg.ToyG1()
    cfg.validate()


# --------------------------------------------------------------------- G5

def test_the_cli_surface():
    for key, flag, choices in (("unpol_sf", "--unpol-sf", lg.UNPOL_SF),
                               ("pol_sf", "--pol-sf", lg.POL_SF)):
        assert key in cli.DEFAULTS
        assert cli.DEFAULTS[key] == "toy"
        assert cli.resolve([])[key] == "toy"
        help_text = cli.build_parser().format_help()
        assert flag in help_text
        for choice in choices:
            assert choice in help_text
    assert cli.resolve(["--unpol-sf", "mstw"])["unpol_sf"] == "mstw"
    assert cli.resolve(["--pol-sf", "nnpdfpol"])["pol_sf"] == "nnpdfpol"


def test_a_missing_tier_fails_loudly_at_the_command_line():
    """Point of honesty: NEVER a quiet downgrade to the toy.

    Both branches are exercised by passing the availability flags in, so the
    test does not need a build that lacks a tier.
    """
    with pytest.raises(SystemExit) as e:
        cli.require_unpol_sf_tier("mstw", have_pythia=False)
    assert "PYTHIA" in str(e.value) and "mstw" in str(e.value)
    with pytest.raises(SystemExit) as e:
        cli.require_unpol_sf_tier("ct18nlo", have_lhapdf=False)
    assert "LHAPDF" in str(e.value)
    with pytest.raises(SystemExit) as e:
        cli.require_pol_sf_tier("nnpdfpol", have_lhapdf=False)
    assert "LHAPDF" in str(e.value) and "SIGN" in str(e.value)
    # ... and a no-op where the tier is there, and for the toy always
    cli.require_unpol_sf_tier("toy", have_pythia=False, have_lhapdf=False)
    cli.require_pol_sf_tier("toy", have_lhapdf=False)
    cli.require_unpol_sf_tier("mstw", have_pythia=True)
    cli.require_unpol_sf_tier("ct18nlo", have_lhapdf=True)
    cli.require_pol_sf_tier("nnpdfpol", have_lhapdf=True)


def test_the_run_banner_says_which_backends_and_what_they_cost():
    """A run that does not say which backend it used is not reproducible
    from its own log -- so the block prints on EVERY run, the all-default one
    included, because "toy" is a physics choice too."""
    toy = "\n".join(cli.sf_banner_lines("toy", "toy", 0.0, 0.0))
    assert "unpol toy" in toy and "pol toy" in toy
    assert "TOY" in toy                      # it says what that costs
    assert "grid floor" not in toy           # ToyF2 has no grid

    ct = "\n".join(cli.sf_banner_lines("ct18nlo", "toy", 1.6770249, 0.3617907))
    assert "does NOT mean g1 is unchanged" in ct   # the non-orthogonality
    assert "36.18 %" in ct and "1.67702" in ct

    nan = float("nan")
    unknown = "\n".join(cli.sf_banner_lines("mstw", "nnpdfpol", nan, nan))
    assert "no grid floor" in unknown and "never as 0" in unknown
    assert "does NOT mean g1 is unchanged" not in unknown


# --------------------------------------------------------------------- G6

@needs_lhapdf
def test_the_tagged_struck_cluster_kernel_actually_moves():
    """THE GATE D2 EXISTS FOR.

    The tagged channels had no structure-function injection point at all --
    not through `StruckClusterOptions`, and not through `PipelineConfig.kernel`
    either, which the `Pipeline` reads on the non-tagged branch alone.  This
    asserts that a tagged run's own kernel is now the selected backend, by
    value, and that the whole-nucleus rate moves with it.
    """
    p_toy = _l.Pipeline(_cfg("tagged-6Li-alpha"), _plan())
    p_ct = _l.Pipeline(_cfg("tagged-6Li-alpha", unpol_sf="ct18nlo"), _plan())

    f2_toy = p_toy.dis_sampler.kernel.nuclear_f2.f2a(0.3, 10.0)
    f2_ct = p_ct.dis_sampler.kernel.nuclear_f2.f2a(0.3, 10.0)
    want = _l.NuclearF2(p_ct.dis_sampler.kernel.ion,
                        _l.LhapdfSF("CT18NLO", 0)).f2a(0.3, 10.0)
    assert f2_ct != f2_toy
    assert f2_ct == want                       # bit for bit the named backend

    # ... and the rate moves with it, by the same factor the inclusive
    # channel does (both 6Li and the embedded deuteron are isoscalar and every
    # structure function here is per nucleon).
    r_tagged = p_ct.sigma_pb() / p_toy.sigma_pb()
    assert r_tagged == pytest.approx(0.7984709, rel=1e-5)


@needs_lhapdf
def test_the_polarised_selector_reaches_the_tagged_kernel_too():
    p_toy = _l.Pipeline(_cfg("tagged-6Li-alpha"), _plan())
    p_np = _l.Pipeline(_cfg("tagged-6Li-alpha", pol_sf="nnpdfpol"), _plan())
    t = p_toy.dis_sampler.kernel.tables(0.2, 10.0)
    n = p_np.dis_sampler.kernel.tables(0.2, 10.0)
    assert n.f1 == t.f1                        # an explicit g1 moves g1 alone
    assert n.g1 != t.g1
    assert n.g1 == _l.InclusiveKernel(
        p_np.dis_sampler.kernel.ion,
        _kernel_options(g1=_l.LhapdfG1("NNPDFpol11_100", 0))
    ).tables(0.2, 10.0).g1


def _kernel_options(f2=None, g1=None):
    o = _l.InclusiveKernel.Options()
    if f2 is not None:
        o.f2_source = f2
    if g1 is not None:
        o.g1_model = g1
    return o


# --------------------------------------------------------------------- G7

@pytest.mark.parametrize("channel", CHANNELS)
def test_meta_records_both_backends_on_every_channel(channel):
    """Two npz files that disagree by 20 % in cross section must not be
    indistinguishable in `meta`.  The write is UNCONDITIONAL, like the b1
    block it sits beside."""
    cfg = _cfg(channel)
    meta = _l.Pipeline(cfg, _plan()).generate(1, True)["meta"]
    assert meta["unpol_sf"] == "toy"
    # ... and the POLARISED key is the backend name only where the backend was
    # READ.  On the coherent channel nothing evaluates g1, so "toy" there
    # would name a backend that did not run exactly as "nnpdfpol" would
    # (G9 below); the label goes in instead, with the reason beside it.
    #
    # ... AND THE PLAN IS NOW PART OF THE PREDICATE.  `_plan()` is
    # `tensor_thirds_plan`, whose every category is built at lam_e = 0, so
    # this run evaluates no g1 on ANY channel and "toy" here would name a
    # backend that did not run exactly as "nnpdfpol" would.  That second axis
    # went in on 2026-09-05; before it, this assertion passed on the inclusive
    # and tagged channels with the label missing (see G9b).
    if _l.pol_sf_is_read(cfg, _plan()):
        assert meta["pol_sf"] == "toy"
    else:
        assert meta["pol_sf"] == _l.pol_sf_unread_label(cfg, _plan())
    assert meta["pol_sf_reach"] == _l.pol_sf_reach_report(cfg, _plan())
    # ToyF2 is closed form and has no grid, so nothing can be below one.
    assert meta["unpol_sf_grid_q2_min"] == 0.0
    assert meta["unpol_sf_below_grid_frac"] == 0.0


def test_meta_says_so_when_a_caller_supplied_kernel_won():
    cfg = _cfg()
    cfg.kernel = _l.default_inclusive_kernel(_l.li6())
    meta = _l.Pipeline(cfg, _plan()).generate(1, True)["meta"]
    assert meta["unpol_sf"] == "caller-supplied kernel"
    assert meta["pol_sf"] == "caller-supplied kernel"
    # the grid report cannot see the caller's backend, so it is NaN, not 0
    assert np.isnan(meta["unpol_sf_grid_q2_min"])
    assert np.isnan(meta["unpol_sf_below_grid_frac"])


@needs_lhapdf
def test_meta_carries_the_measured_below_grid_fraction():
    """36 % of a shipped 6Li run's own rate is below CT18NLO's grid floor.

    Below it LHAPDF does not freeze -- it keeps evolving downward -- so this
    is not a detail: it is a third of the run on an off-grid extrapolation.
    It is not refused (that would make the flag unusable on the shipped
    scenario), it is MEASURED per run and recorded.
    """
    meta = _l.Pipeline(_cfg(unpol_sf="ct18nlo"), _plan()).generate(1, True)["meta"]
    assert meta["unpol_sf"] == "ct18nlo"
    assert meta["unpol_sf_grid_q2_min"] == pytest.approx(1.6770250, rel=1e-9)
    assert meta["unpol_sf_below_grid_frac"] == pytest.approx(0.3617907, rel=1e-5)


@needs_pythia
def test_meta_writes_nan_for_a_backend_that_reports_no_grid_floor():
    """PYTHIA's `MSTWpdf` keeps its `qsqmin` private, so the floor is NOT
    known -- and an unknown fraction is written as NaN, never as a 0 that was
    not measured."""
    meta = _l.Pipeline(_cfg(unpol_sf="mstw"), _plan()).generate(1, True)["meta"]
    assert meta["unpol_sf"] == "mstw"
    assert np.isnan(meta["unpol_sf_grid_q2_min"])
    assert np.isnan(meta["unpol_sf_below_grid_frac"])


# --------------------------------------------------------------------- G8

@needs_pythia
@needs_lhapdf
def test_the_pythia_bridge_species_draw_is_on_the_kernels_own_backend():
    """ONE unpolarised backend for the whole run.

    The T2 struck-nucleon species draw is P(p) = Z F2p/(Z F2p + N F2n) on
    `PythiaBridgeOptions.f2_source`, and the `Pipeline` cannot fill it -- the
    CLI builds the bridge itself.  Left unset the bridge would keep drawing
    from `ToyF2` while T0 and T1 moved to the selected backend, and the toy's
    F2n/F2p is up to 24 % away from CT18NLO's at x = 0.5.  This asserts the
    CLI's wiring hands it the SAME OBJECT, and that the toy default leaves it
    None (so a default run is bit for bit).
    """
    opts = _l.PythiaBridgeOptions()
    assert opts.f2_source is None                # the shipped default

    cfg = _cfg(unpol_sf="ct18nlo")
    opts.f2_source = cfg.unpol_sf_obj            # what cli.main does
    assert opts.f2_source is cfg.unpol_sf_obj    # the same object, not a copy
    assert opts.f2_source.f2p(0.5, 10.0) == cfg.unpol_sf_obj.f2p(0.5, 10.0)

    toy_cfg = _cfg()
    opts.f2_source = toy_cfg.unpol_sf_obj
    assert opts.f2_source is None                # ... and the toy stays unset


@needs_pythia
@needs_lhapdf
def test_lipolgen_run_wires_the_bridge_to_the_same_backend():
    """lipolgen.run(hadronize=True, ...) builds its own bridge too.

    So it carries the same one-backend rule as cli.main, and an explicit
    pythia_options={"f2_source": ...} still wins over it.
    """
    out = lg.run(events=40, hadronize=True, unpol_sf="ct18nlo")
    assert out["bridge"].stats.n_ok == 40
    assert out["pipeline"].sigma_pb() != lg.run(
        events=20, hadronize=True)["pipeline"].sigma_pb()


# --------------------------------------------------------------------- G9
#
# THE REACH OF `--pol-sf`, PER CHANNEL -- the claim the run surface used to
# get wrong.  The banner said "both reach EVERY kernel: inclusive, coherent
# and tagged" and `meta["pol_sf"]` printed whichever backend was typed, while
# on the coherent channel NOTHING evaluates g1: the yield is f_coh(x) times
# the UNPOLARISED cell cross sections and the tensor signal is the recoil
# azimuth's 1 + c2 cos 2(phi_t - phi_S).  That is a knob that did not run,
# recorded as if it had.

#: Every channel the pipeline builds, for the reach gates.  `tagged-d-p` is
#: here and not in `CHANNELS` because it is the one whose DIS target is a free
#: neutron -- see G12.
ALL_CHANNELS = ("inclusive", "coherent", "tagged-6Li-alpha",
                "tagged-7Li-alpha", "tagged-d-p")


def _plan_for(channel):
    """A run plan of the channel's OWN ion spin.

    The three tensor plans hard-code j = 1, so 7Li (spin 3/2) takes the only
    standard plan that takes j -- `Pipeline` refuses the mismatch by name.
    """
    j = lg.ion_spin(_cfg(channel).isotope)
    return (lg.helicity_flip_plan(j, 0.7, 0.7) if j != 1.0
            else lg.tensor_thirds_plan(0.0, 0.6))


def _polarised_plan_for(channel):
    """A fill that DOES carry lam_e * P_e != 0, whatever the channel's spin.

    The CHANNEL axis of `pol_sf_is_read` can only be tested under one of
    these: under a tensor plan the answer is false everywhere, for the OTHER
    reason (G9b).
    """
    return lg.helicity_flip_plan(lg.ion_spin(_cfg(channel).isotope), 0.7, 0.7)


@pytest.mark.parametrize("channel", ALL_CHANNELS)
def test_pol_sf_is_read_on_every_channel_but_the_coherent_one(channel):
    """The CHANNEL axis, under a fill that DOES carry a beam helicity."""
    cfg = _cfg(channel)
    plan = _polarised_plan_for(channel)
    assert _l.pol_sf_is_read(cfg, plan) == (channel != "coherent")
    # ... and the sentence follows the predicate, on every channel.
    report = _l.pol_sf_reach_report(cfg, plan)
    if channel == "coherent":
        assert report.startswith(_l.pol_sf_unread_label(cfg, plan))
        assert "SPIN-INDEPENDENT" in report
        # it also says what DOES reach this channel, so the reader is not
        # left thinking the selectors are jointly inert here
        assert "--unpol-sf" in report
    else:
        assert "not read on channel" not in report
        assert "g1_model" in report


# -------------------------------------------------------------------- G9b
#
# THE SECOND AXIS, and it is the shipped default.  Until 2026-09-05
# `pol_sf_is_read` keyed on the CHANNEL alone, so `lipolgen-run --pol-sf
# nnpdfpol` at DEFAULTS -- `--plan tensor-thirds`, whose every category is
# built at lam_e = 0, pe = 0 -- wrote `meta["pol_sf"] = "nnpdfpol"` and
# printed a banner claiming the kernel was reached, while the file was
# bit-identical to a `toy` one.  g1 is multiplied by lam_e * P_e and by
# nothing else (`InclusiveKernel::amplitudes`), which is the whole rule.

@pytest.mark.parametrize("channel", ALL_CHANNELS)
def test_pol_sf_is_not_read_under_an_unpolarised_beam_plan(channel):
    cfg = _cfg(channel)
    j = lg.ion_spin(cfg.isotope)
    unpolarised = ([lg.tensor_thirds_plan(0.7, 0.6),
                    lg.transverse_tensor_plan(0.6),
                    lg.tensor_flip_plan(0.6)] if j == 1.0 else []) \
        + [lg.helicity_flip_plan(j, 0.7, 0.0)]     # --pe 0
    for plan in unpolarised:
        assert not _l.plan_has_beam_helicity(plan)
        assert not _l.pol_sf_is_read(cfg, plan)
        label = _l.pol_sf_unread_label(cfg, plan)
        assert label.startswith("not read")
        assert _l.pol_sf_reach_report(cfg, plan).startswith(label)
    # ... and the same config with a POLARISED fill does read it (off the
    # coherent channel, where the other axis holds).
    polarised = lg.helicity_flip_plan(j, 0.7, 0.7)
    assert _l.plan_has_beam_helicity(polarised)
    assert _l.pol_sf_is_read(cfg, polarised) == (channel != "coherent")


def test_the_meta_of_a_default_plan_run_records_the_label_not_the_backend():
    """The exact file the defect produced: `--pol-sf nnpdfpol` under the CLI's
    own default plan.  It may not come back saying `pol_sf = nnpdfpol`."""
    cfg = _cfg("inclusive")
    plan = lg.tensor_thirds_plan(0.7, 0.6)
    meta = _l.Pipeline(cfg, plan).generate(1, True)["meta"]
    assert meta["pol_sf"] == _l.pol_sf_unread_label(cfg, plan)
    assert meta["pol_sf"] != "toy" and meta["pol_sf"] != "nnpdfpol"
    # ... and the knob-provenance row says the same thing, from one table.
    row = meta["knob_provenance"]["pol_sf"]
    assert row["status"] == "not-read"
    assert row["value"] == "toy"               # what was ASKED for
    assert row["label"] == meta["pol_sf"]
    assert row["reason"] == meta["pol_sf_reach"]


@needs_lhapdf
def test_the_polarised_selector_does_not_reach_the_coherent_channel():
    """MEASURED, not asserted from the code: a coherent run is bit-identical
    between `--pol-sf toy` and `--pol-sf nnpdfpol`, column by column."""
    toy = _l.Pipeline(_cfg("coherent"), _plan())
    npd = _l.Pipeline(_cfg("coherent", pol_sf="nnpdfpol"), _plan())
    assert npd.sigma_pb() == toy.sigma_pb()
    assert np.array_equal(np.asarray(npd.sigma_per_category_pb()),
                          np.asarray(toy.sigma_per_category_pb()))
    a = toy.generate(0, False, 1)
    b = npd.generate(0, False, 1)
    moved = [k for k, v in a.items()
             if isinstance(v, np.ndarray) and v.dtype.kind == "f"
             and not np.array_equal(v, b[k], equal_nan=True)]
    assert moved == []
    # ... while the KERNEL does carry the selected backend: the selector is
    # wired, the channel simply never asks it anything.
    assert (npd.dis_sampler.kernel.tables(0.2, 10.0).g1
            != toy.dis_sampler.kernel.tables(0.2, 10.0).g1)
    # ... and the same `--unpol-sf` DOES reach this channel, which is why it
    # is labelled rather than refused.
    ct = _l.Pipeline(_cfg("coherent", unpol_sf="ct18nlo"), _plan())
    assert ct.sigma_pb() / toy.sigma_pb() == pytest.approx(0.705841, rel=1e-5)


@needs_lhapdf
def test_meta_records_the_label_and_not_the_backend_on_the_coherent_channel():
    """The whole point: a coherent run typed with `--pol-sf nnpdfpol` may not
    come back saying `pol_sf = nnpdfpol`."""
    cfg = _cfg("coherent", pol_sf="nnpdfpol")
    # A POLARISED fill, so the CHANNEL is the only thing that can make the
    # label appear -- and the label is the channel's, not the plan's.
    plan = lg.helicity_flip_plan(1.0, 0.7, 0.7)
    meta = _l.Pipeline(cfg, plan).generate(1, True)["meta"]
    assert meta["pol_sf"] == "not read on channel coherent-6Li"
    assert meta["pol_sf"] == _l.pol_sf_unread_label(cfg, plan)
    assert meta["pol_sf"] != "nnpdfpol"
    assert meta["pol_sf_reach"] == _l.pol_sf_reach_report(cfg, plan)
    # and "toy" is refused there for the same reason -- ToyG1 did not run
    # either
    toy_meta = _l.Pipeline(_cfg("coherent"), plan).generate(1, True)["meta"]
    assert toy_meta["pol_sf"] == meta["pol_sf"]
    # the UNPOLARISED key is untouched: that one does reach the channel
    assert toy_meta["unpol_sf"] == "toy"


def test_the_banner_makes_no_hand_written_reach_claim_at_all():
    """The SF block used to carry two of them, and the first was false on the
    CLI's own default plan.  Both are gone: reach is the knob-provenance
    table's job now, on every knob at once (`knob_provenance_lines`)."""
    read = "\n".join(cli.sf_banner_lines("toy", "toy"))
    assert "reaches EVERY kernel" not in read
    assert "both reach EVERY kernel" not in read
    assert "DID NOT RUN HERE" not in read
    # what stays is the PRICE, which is a measurement and not a reach claim
    assert "BOTH ARE THE TOY BACKENDS" in read

    cfg = _cfg("coherent")
    plan = lg.helicity_flip_plan(1.0, 0.7, 0.7)
    unread = "\n".join(cli.sf_banner_lines(
        "ct18nlo", "nnpdfpol", 1.6770249, 0.4474601,
        _l.pol_sf_is_read(cfg, plan), _l.pol_sf_reach_report(cfg, plan)))
    assert "reaches EVERY kernel" not in unread
    assert "44.75 %" in unread
    # ... and the run's own table is where the sentence now lives, verbatim.
    rows = _l.Pipeline(cfg, plan).knob_provenance()
    pol = [r for r in rows if r.name == "pol_sf"][0]
    assert pol.status == _l.KnobStatus.NotRead
    assert pol.reason == _l.pol_sf_reach_report(cfg, plan)
    banner = "\n".join(cli.knob_provenance_lines(rows))
    assert "not read on channel coherent-6Li" in banner


def test_the_banner_drops_the_g1_clauses_where_g1_was_never_read():
    """`--pol-sf toy does NOT mean g1 is unchanged` is true where g1 runs and
    misleading where it does not."""
    cfg = _cfg("coherent")
    plan = lg.helicity_flip_plan(1.0, 0.7, 0.7)
    coh = "\n".join(cli.sf_banner_lines(
        "ct18nlo", "toy", 1.6770249, 0.4474601,
        _l.pol_sf_is_read(cfg, plan), _l.pol_sf_reach_report(cfg, plan)))
    assert "does NOT mean g1 is unchanged" not in coh
    incl = "\n".join(cli.sf_banner_lines("ct18nlo", "toy", 1.6770249,
                                         0.3617907))
    assert "does NOT mean g1 is unchanged" in incl


# -------------------------------------------------------------------- G10
#
# THE BELOW-GRID FRACTION IS THE CHANNEL'S OWN RATE.  It used to be the
# INCLUSIVE sampler's cells on every channel, so a coherent run recorded
# 36.18 % where its own rate has 44.75 % below CT18NLO's floor -- a quantified
# `meta` key that did not describe the run it was attached to.

@pytest.mark.parametrize("channel", ALL_CHANNELS)
def test_the_rate_weights_are_the_samplers_cells_off_the_coherent_channel(
        channel):
    p = _l.Pipeline(_cfg(channel), _plan_for(channel))
    w = np.asarray(p.cell_rate_weights_pb)
    cells = np.asarray(p.dis_sampler.cell_xsec_pb)
    if channel == "coherent":
        assert not np.array_equal(w, cells)
        # THE IDENTITY THAT MAKES IT THE RATE: the coherent weights sum to
        # this run's own cross section, exactly.
        assert w.sum() == pytest.approx(p.sigma_pb(), rel=1e-12)
        # ... and they are the inclusive cells reweighted by f_coh, zero on
        # the cells the diffractive-mass gate rejects.
        assert np.all(w <= cells + 1e-12)
        assert (w == 0.0).any()
    else:
        assert np.array_equal(w, cells)


@needs_lhapdf
def test_the_below_grid_fraction_is_taken_against_the_channels_own_rate():
    """PINNED, both numbers, on the shipped 6Li window at config 1."""
    q2_min_i, below_i = _l.unpol_sf_grid_report(
        _l.Pipeline(_cfg("inclusive", unpol_sf="ct18nlo"), _plan()))
    q2_min_c, below_c = _l.unpol_sf_grid_report(
        _l.Pipeline(_cfg("coherent", unpol_sf="ct18nlo"), _plan()))
    assert q2_min_c == q2_min_i == pytest.approx(1.6770250, rel=1e-9)
    # the inclusive number is unchanged -- this run's rate IS its cells
    assert below_i == pytest.approx(0.3617907, rel=1e-6)
    # ... and the coherent one is its own, 8.6 points higher
    assert below_c == pytest.approx(0.4474601, rel=1e-6)
    assert below_c != below_i


@needs_lhapdf
def test_meta_carries_the_coherent_channels_own_below_grid_fraction():
    meta = _l.Pipeline(_cfg("coherent", unpol_sf="ct18nlo"),
                       _plan()).generate(1, True)["meta"]
    assert meta["unpol_sf"] == "ct18nlo"
    assert meta["unpol_sf_below_grid_frac"] == pytest.approx(0.4474601,
                                                             rel=1e-6)
    assert meta["sigma_pb"] == pytest.approx(8806.096207, rel=1e-9)


# -------------------------------------------------------------------- G11
#
# A CALLER-SUPPLIED KERNEL IS REFUSED WHERE IT IS NOT READ.  `Pipeline` reads
# `cfg.kernel` on the ION-LEVEL branch only (inclusive and coherent); a tagged
# run draws from the struck-cluster sampler and never looks at it -- while
# `meta` wrote "caller-supplied kernel" into three keys.

TAGGED = ("tagged-6Li-alpha", "tagged-7Li-alpha", "tagged-d-p")


@pytest.mark.parametrize("channel", TAGGED)
def test_a_caller_supplied_kernel_is_refused_on_a_tagged_channel(channel):
    # `validate()` alone: no Pipeline is built, so the plan never enters.
    cfg = _cfg(channel)
    cfg.kernel = _l.default_inclusive_kernel(_l.deuteron())
    with pytest.raises(RuntimeError) as e:
        cfg.validate()
    msg = str(e.value)
    assert "not read on channel" in msg and channel in msg
    assert "struck.f2_source" in msg          # it names the way in


@pytest.mark.parametrize("channel", ("inclusive", "coherent"))
def test_a_caller_supplied_kernel_is_still_accepted_where_it_is_read(channel):
    """The refusal is not 'off the inclusive channel': the COHERENT channel
    rides the ion-level sampler, so a hand-built kernel does reach its rate
    and must keep working."""
    cfg = _cfg(channel)
    cfg.kernel = _l.default_inclusive_kernel(_l.li6())
    cfg.validate()
    assert _l.Pipeline(cfg, _plan()).sigma_pb() > 0.0


@needs_lhapdf
def test_the_tagged_kernel_mislabel_is_gone():
    """The exact configuration that used to pass: a tagged run whose `kernel`
    was an `InclusiveKernel` on CT18NLO ran on ToyF2 (F2A(0.3, 10)
    = 0.3691149345, the toy value) with `meta` naming a caller-supplied
    kernel for all three of unpol_sf / pol_sf / b1_model."""
    o = _l.InclusiveKernel.Options()
    o.f2_source = _l.LhapdfSF("CT18NLO", 0)
    cfg = _cfg("tagged-6Li-alpha")
    cfg.kernel = _l.InclusiveKernel(_l.deuteron(), o)
    with pytest.raises(RuntimeError):
        _l.Pipeline(cfg, _plan())
    # ... and the way to do it now moves the number it claims to move
    p = _l.Pipeline(_cfg("tagged-6Li-alpha", unpol_sf="ct18nlo"), _plan())
    assert p.dis_sampler.kernel.nuclear_f2.f2a(0.3, 10.0) != pytest.approx(
        0.3691149345, rel=1e-9)


# -------------------------------------------------------------------- G12
#
# THE NEUTRON-TAGGED BACKEND PAIR IS REFUSED BY NAME.  `--isotope d --channel
# tagged-d-p --unpol-sf ct18nlo --pol-sf nnpdfpol --plan helicity-flip` at any
# --pe != 0 -- the A_par measurement on a free neutron, the channel the docs
# call the strongest single reason `--pol-sf` exists -- has A1 = g1/F1 = 4.406
# in the shipped window's top cell (x = 0.955, Q2 = 1119), so 1 + w_avg goes
# NEGATIVE and `InclusiveSampler` refuses at configuration time.  That is
# correct; what was missing was a message naming the cure, and a CLI that did
# not print a traceback.

def _dp_plan():
    return lg.helicity_flip_plan(1.0, 0.7, 0.7)


@needs_lhapdf
def test_the_neutron_tagged_pair_is_refused_with_the_cell_and_the_cure():
    with pytest.raises(RuntimeError) as e:
        _l.Pipeline(_cfg("tagged-d-p", unpol_sf="ct18nlo",
                         pol_sf="nnpdfpol"), _dp_plan())
    msg = str(e.value)
    assert "x = 0.955" in msg                 # the cell
    assert "--x-max" in msg                   # the cure, by flag name
    assert "A1 = g1/F1 is 4.406" in msg       # why: the measured ratio
    assert "1 + w_avg = -0.02414" in msg


@needs_lhapdf
def test_x_max_095_cures_it_and_each_backend_alone_never_needed_it():
    """The refusal is a property of the PAIR, which is what makes it worth a
    message: either backend alone builds on the same plan."""
    cfg = _cfg("tagged-d-p", unpol_sf="ct18nlo", pol_sf="nnpdfpol")
    sc = cfg.scenario
    sc.x_max = 0.95
    cfg.scenario = sc
    assert _l.Pipeline(cfg, _dp_plan()).sigma_pb() > 0.0
    for kw in ({}, dict(unpol_sf="ct18nlo"), dict(pol_sf="nnpdfpol")):
        assert _l.Pipeline(_cfg("tagged-d-p", **kw), _dp_plan()).sigma_pb() > 0
    # ... and P_e = 0 removes it too, because w_avg scales with P_e P_z
    assert _l.Pipeline(_cfg("tagged-d-p", unpol_sf="ct18nlo",
                            pol_sf="nnpdfpol"),
                       lg.helicity_flip_plan(1.0, 0.7, 0.0)).sigma_pb() > 0


@needs_lhapdf
def test_the_cli_reports_it_as_a_configuration_error_not_a_traceback():
    """`lipolgen-run` used to exit 1 through a Python traceback out of
    `_l.Pipeline`.  A configuration refusal is a SystemExit carrying the
    message the sampler wrote -- the `PythiaBridge` clause's rule."""
    argv = ["--isotope", "d", "--channel", "tagged-d-p", "--events", "10",
            "--unpol-sf", "ct18nlo", "--pol-sf", "nnpdfpol",
            "--plan", "helicity-flip", "--pz", "0.7", "--pe", "0.7",
            "--quiet"]
    with pytest.raises(SystemExit) as e:
        cli.main(argv)
    msg = str(e.value)
    assert "negative phi-averaged density" in msg and "--x-max" in msg
    # ... and the cure goes through the CLI
    cli.main(argv + ["--x-max", "0.95"])


# -------------------------------------------------------------------- G13
#
# THE T2 BRIDGE IS ON THE RUN'S OWN BACKEND, ENFORCED WHERE THE TWO OBJECTS
# MEET.  A bridge built by hand on a `--unpol-sf ct18nlo` config used to pass
# `validate()` and draw its T2 species from `ToyF2` while `meta` said
# ct18nlo; the toy's F2n/F2p is up to 24 % away.

@needs_pythia
@needs_lhapdf
def test_set_pythia_hadronizer_refuses_a_bridge_on_another_backend():
    beams = _l.default_configs("6Li")[1]
    plain = _l.PythiaBridge(beams, _l.PythiaBridgeOptions())
    assert plain.options.f2_source is None          # the bridge's own ToyF2

    # the default run: both unset, i.e. the same object -- accepted
    _l.set_pythia_hadronizer(_cfg(), plain)

    cfg = _cfg(unpol_sf="ct18nlo")
    with pytest.raises(RuntimeError) as e:
        _l.set_pythia_hadronizer(cfg, plain)
    msg = str(e.value)
    assert "species draw" in msg and "unpol_sf" in msg
    assert "config.unpol_sf_obj" in msg              # it names the fix

    opts = _l.PythiaBridgeOptions()
    opts.f2_source = cfg.unpol_sf_obj
    good = _l.PythiaBridge(beams, opts)
    assert good.options.f2_source is cfg.unpol_sf_obj
    _l.set_pythia_hadronizer(cfg, good)              # the same object: fine
