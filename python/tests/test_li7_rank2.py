"""D1: the ⁷Li rank-2 (tensor) sector is identically zero -- OUT LOUD.

WHAT THIS FILE GUARDS.  Before 2026-09-04 a ⁷Li inclusive run's whole rank-2
sector -- the tensor term of the φ-averaged rate, the cos 2φ (gluon
transversity) amplitude and therefore A_zz -- was exactly, bit-for-bit zero,
and NOTHING on the run surface said so: `default_inclusive_kernel` fills
`b1_func`/`delta_func` inside `if (kernel_fills_rank2(ion))`, which is spin 1,
so `b1_32_func / b2_32_func / delta_32_func` stayed unset and
`InclusiveKernel::tables` returned 0.0 for all three (an unset `b2` is
2x*b1 and therefore 0 too).  Meanwhile the npz `meta` recorded
`b1_model = "miller"` -- a backend that did not run -- and the run banner
printed the b1 block only when `b1_model != Miller`, which a ⁷Li run can never
reach, so it printed nothing at all.  Measured then and re-measured here: the
two per-category cross sections of a T = +1 against T = -1 fill are the SAME
DOUBLE and the asymmetry is exactly 0.0.

The physics -- why ⁷Li's rank-2 input is purely orbital, the measured leading
α-t estimate, and why its SIGN is not predicted by anything in this tree (the
`ToyF2 -> MSTW2008 LO` swap flips it at four of six x points), which is why no
b1 ships -- is `docs/open_items/run_2026-09-03/phase_D_li7_rank2.md`, and the
deferred design is open item 15 of `docs/OPEN_ITEMS_SOLUTIONS.md`.

The gates, in order:

G1  the zero is REAL and is now RECORDED: the two categories are the same
    double, `meta["rank2_input"]` says so in words, and `meta["b1_model"]` /
    `meta["b1_unpol"]` no longer name a backend that did not run.
G2  the BANNER prints the same sentence the `meta` records, on every ⁷Li
    inclusive run, tensor plan or not.
G3  nothing else moved: ⁶Li, the deuteron, a caller-supplied kernel and every
    non-inclusive channel keep their old labels.
G4  the two b1 SCALES cannot record a variation on ⁷Li -- `validate()`
    already refuses them, which is why they stay numeric.
G5  `make_plan` REFUSES the three spin-1-only tensor plans at J = 3/2 instead
    of building a spin-1 plan that throws three frames down (defect F2), and
    its message names what does work.
G6  the `--unpol-sf` / `--pol-sf` selectors REACH ⁷Li -- inclusive and tagged
    -- which is the nucleon input a ⁷Li b1 would fold against; and the tensor
    sector stays exactly zero under every one of them.
G7  the two adjacent run-surface defects: `ClusterPartialWave`'s half-assigned
    segfault (F4) and `spin32_populations`' uninformative refusal (F3).
"""

import subprocess
import sys

import numpy as np
import pytest

import lipolgen as lg
from lipolgen import cli
from lipolgen import _lipolgen as _l

N = 400
SEED = 11

#: What `meta["b1_model"]` and `meta["b1_unpol"]` say when nothing filled a
#: rank-2 slot.  The C++ `rank2_none_label()` is the definition; this is the
#: test's copy of it, and G1 asserts they agree.
NONE_LABEL = "none (spin 3/2: no rank-2 input)"


def li7_tensor_plan(name="T"):
    """The honest ⁷Li A_zz plan: pure |m| = 3/2 (T = +1) against |m| = 1/2.

    Built by hand because there is no spin-3/2 tensor plan in the tree --
    `tensor_thirds_plan`, `transverse_tensor_plan` and `tensor_flip_plan` all
    hard-code j = 1 (`src/core/bookkeeping.cpp`), which is defect F1 of the
    note and is what G5 refuses rather than papers over.  For J = 3/2
    Q_NN(±3/2) = +1 and Q_NN(±1/2) = -1, so these two fills are T = ±1 and the
    contrast is the whole rank-2 observable.
    """
    return _l.RunPlan(
        [_l.SpinCategory(name + "+", 1.5, [0.5, 0.0, 0.0, 0.5],
                         0, 0.0, 0.0, 0.0, 0.5),
         _l.SpinCategory(name + "-", 1.5, [0.0, 0.5, 0.5, 0.0],
                         0, 0.0, 0.0, 0.0, 0.5)],
        0.0, 0.0, 1.0)


def _cfg(isotope="7Li", **kw):
    kw.setdefault("events", N)
    kw.setdefault("seed", SEED)
    return lg.make_config(isotope=isotope, **kw)


# ------------------------------------------------------------------- G1

def test_the_7li_tensor_sector_is_exactly_zero_and_the_run_records_it():
    """The zero, the same double, and a `meta` that no longer names Miller."""
    cfg = _cfg()
    assert _l.inclusive_rank2_is_empty(cfg)
    plan = li7_tensor_plan()
    p = _l.Pipeline(cfg, plan)
    sig = p.sigma_per_category_pb()
    # NOT approx: the two are the same double.  There is no "small" here --
    # the structure functions are absent, not suppressed.
    assert sig[0] == sig[1]
    assert (sig[0] - sig[1]) == 0.0
    meta = p.generate(0, False)["meta"]
    assert meta["b1_model"] == NONE_LABEL == _l.rank2_none_label()
    assert meta["b1_unpol"] == NONE_LABEL
    r = meta["rank2_input"]
    assert r == _l.rank2_input_report(cfg, plan)
    assert r.startswith("EMPTY -- 7Li is spin 3/2")
    assert "IDENTICALLY ZERO, not small" in r
    assert "b1_32 = b2_32 = delta_32 = 0" in r
    # ... and the run that ASKED for the sector is told that it asked.
    assert "THIS RUN'S PLAN CARRIES A RANK-2 FILL (T = 1)" in r


def test_the_kernel_itself_returns_zero_for_all_three_rank2_slots():
    """One level below the run: the tables, not the cross sections."""
    k = _l.default_inclusive_kernel(_l.ion_by_name("7Li"))
    t = k.tables(0.2, 5.0)
    assert (t.b1, t.b2, t.delta) == (0.0, 0.0, 0.0)
    # the SAME call on ⁶Li is not zero -- so this is the input, not the code
    t6 = _l.default_inclusive_kernel(_l.ion_by_name("6Li")).tables(0.2, 5.0)
    assert t6.b1 != 0.0 and t6.delta != 0.0


def test_the_zero_is_recorded_even_without_a_tensor_plan():
    """A vector (helicity-flip) ⁷Li run gets the sentence too.

    The user who has to be told is the one who did not know to ask, so the
    key and the banner are unconditional on a ⁷Li inclusive run -- and a
    J = 3/2 vector fill carries a rank-2 moment anyway (T = 0.4 at
    P_z = 0.7 on the max-entropy ladder), so the "asked for it" clause fires
    even here.
    """
    cfg = _cfg()
    plan = lg.make_plan("helicity-flip", j=1.5, pz=0.7, pe=0.7)
    meta = _l.Pipeline(cfg, plan).generate(0, False)["meta"]
    assert meta["rank2_input"].startswith("EMPTY -- 7Li")
    assert meta["b1_model"] == NONE_LABEL
    assert plan.pzz_true == pytest.approx(0.4, abs=1e-12)


# ------------------------------------------------------------------- G2

def test_the_banner_prints_the_sentence_the_meta_records():
    """One definition, two surfaces (`rank2_zero_lines` wraps, never rewrites)."""
    cfg = _cfg()
    plan = li7_tensor_plan()
    report = _l.rank2_input_report(cfg, plan)
    lines = cli.rank2_zero_lines(report, "7Li")
    assert lines[0].startswith("  7Li rank-2: ")
    # the report survives the wrap word for word
    joined = " ".join(l.strip() for l in lines)
    assert report.split(" (docs/OPEN_ITEMS")[0] in joined
    # ... and the advice names both escape routes, plus the measured
    # hand-supplied-slot check that shows the machinery is fine
    assert "--isotope 6Li" in joined
    assert "tagged-7Li-alpha" in joined
    assert "A_T = -0.043258" in joined


def test_the_cli_prints_the_block_on_a_7li_inclusive_run():
    out = subprocess.run(
        [sys.executable, "-m", "lipolgen.cli", "--isotope", "7Li",
         "--channel", "inclusive", "--plan", "helicity-flip",
         "--events", "100"],
        capture_output=True, text=True, check=True).stdout
    assert "7Li rank-2: EMPTY" in out
    assert "IDENTICALLY ZERO" in out
    assert "tagged-7Li-alpha" in out
    # ... and a ⁶Li run does not print it
    out6 = subprocess.run(
        [sys.executable, "-m", "lipolgen.cli", "--events", "100"],
        capture_output=True, text=True, check=True).stdout
    assert "rank-2: EMPTY" not in out6


def test_the_run_banner_says_which_moments_the_fill_actually_has():
    """`--pzz` is read by the three spin-1 tensor plans always, and by
    `helicity-flip` only at `--pzz-mode typed`.

    At the DEFAULT `--pzz-mode ladder` `helicity_flip_plan` leaves
    `use_explicit_pzz` false and takes the max-entropy ladder at `--pz`, so
    `--plan helicity-flip --pzz 0.6` produced T = 0.4 at J = 3/2 with nothing
    saying so until 2026-09-05.  The banner prints the plan's OWN recorded
    moments -- and, since 2026-09-06, the mode and the fill the OTHER mode
    would have built -- which closes that without moving any fill.  The
    DEFAULT fill is unchanged, and that is what this test checks first.
    """
    out = subprocess.run(
        [sys.executable, "-m", "lipolgen.cli", "--isotope", "7Li",
         "--plan", "helicity-flip", "--pzz", "0.6", "--events", "100"],
        capture_output=True, text=True, check=True).stdout
    assert "fill helicity-flip: J = 1.5, P_z = 0.7, T = 0.4" in out
    assert "--pzz-mode ladder (the default)" in out
    assert "--pzz = 0.6 is NOT read" in out
    out6 = subprocess.run(
        [sys.executable, "-m", "lipolgen.cli", "--events", "100"],
        capture_output=True, text=True, check=True).stdout
    assert "fill tensor-thirds: J = 1, P_z = 0.7, P_zz = 0.6" in out6
    assert "NOT read" not in out6
    assert "--pzz-mode" not in out6


# ------------------------------------------------------------------- G3

def test_nothing_else_is_relabelled():
    """⁶Li, the deuteron, a caller kernel and the other channels are untouched."""
    six = _cfg("6Li")
    assert not _l.inclusive_rank2_is_empty(six)
    meta = _l.Pipeline(six, lg.tensor_thirds_plan(0.0, 0.6)) \
             .generate(0, False)["meta"]
    assert meta["b1_model"] == "miller" and meta["b1_unpol"] == "toy"
    assert meta["rank2_input"].startswith("b1_model = miller (spin 1")

    d = _cfg("d")
    assert not _l.inclusive_rank2_is_empty(d)
    assert _l.rank2_input_report(d, lg.tensor_thirds_plan(0.0, 0.6)) \
             .startswith("b1_model = miller (spin 1")

    # a caller-supplied kernel is the caller's business, on either isotope
    ck = _cfg()
    ck.kernel = _l.default_inclusive_kernel(_l.ion_by_name("7Li"))
    assert not _l.inclusive_rank2_is_empty(ck)
    meta = _l.Pipeline(ck, li7_tensor_plan()).generate(0, False)["meta"]
    assert meta["b1_model"] == "caller-supplied kernel"
    assert meta["rank2_input"] == "caller-supplied kernel"


def test_a_tagged_7li_run_is_not_labelled_with_the_inclusive_zero():
    """The tagged α-t alignment IS carried -- in the event weight.

    Saying "identically zero" there would be a false claim: the ⟨P₂⟩ = -T/5
    polarimeter of that channel is a passing test today.
    """
    cfg = _cfg(channel="tagged-7Li-alpha")
    assert not _l.inclusive_rank2_is_empty(cfg)
    plan = lg.make_plan("helicity-flip", j=1.5, pz=0.7, pe=0.7)
    meta = _l.Pipeline(cfg, plan).generate(0, False)["meta"]
    r = meta["rank2_input"]
    assert r.startswith("not read on channel tagged-7Li-alpha")
    assert "in the event weight" in r
    assert "EMPTY" not in r
    # the b1 labels are the flag's own defaults there, not the ⁷Li label
    assert meta["b1_model"] == "miller"


def test_the_coherent_channel_says_where_its_tensor_signal_is():
    cfg = _cfg("6Li", channel="coherent")
    r = _l.rank2_input_report(cfg, lg.tensor_thirds_plan(0.0, 0.6))
    assert r.startswith("not read on channel coherent-6Li")
    assert "recoil azimuth" in r


# ------------------------------------------------------------------- G4

@pytest.mark.parametrize("setter,value", [
    ("b1_band_scale", 2.0),
    ("b1_alpha_d_dwave_weight", 0.5),
])
def test_the_two_b1_scales_cannot_record_a_variation_on_7li(setter, value):
    """Why they stay NUMERIC in `meta` while the two names became a label.

    `validate()`'s Miller branch already refuses any value but 1.0, and ⁷Li
    can only ever be on Miller (the next test), so neither scale can record a
    variation that did not run.  Changing their TYPE per isotope would break
    every consumer that reads them as floats.
    """
    cfg = _cfg()
    setattr(cfg, setter, value)
    with pytest.raises(RuntimeError, match="do not apply to b1_model = miller"):
        cfg.validate()


def test_7li_cannot_reach_a_b1_backend_or_a_b1_unpol_at_all():
    cfg = _cfg()
    cfg.b1_model = _l.B1Model.Li6Convolution
    with pytest.raises(RuntimeError, match="is 6Li ONLY, got isotope 7Li"):
        cfg.validate()
    cfg = _cfg()
    _l.set_b1_unpol(cfg, _l.B1UnpolSource.Toy)
    cfg.b1_unpol = _l.B1UnpolSource.Mstw
    with pytest.raises(RuntimeError, match="read ONLY by b1_model"):
        cfg.validate()


# ------------------------------------------------------------------- G5

@pytest.mark.parametrize("name", ["tensor-thirds", "azz", "transverse-tensor",
                                  "cos2phi", "tensor-flip", "flip"])
def test_make_plan_refuses_the_spin1_tensor_plans_at_j_32(name):
    """Defect F2: two of them USED to build a spin-1 plan and throw later.

    `make_plan("transverse-tensor", j=1.5)` returned categories at j = 1.0 and
    the run then died inside `InclusiveKernel::amplitudes` with "spin state
    J = 1.000000 is not the kernel's ion spin 1.500000" -- and
    `tensor-thirds`' own refusal ADVISED that plan.
    """
    with pytest.raises(ValueError) as e:
        lg.make_plan(name, j=1.5, pz=0.0, pzz=0.6)
    msg = str(e.value)
    assert "SPIN-1 pattern" in msg
    assert "helicity-flip" in msg           # what does take j
    assert "tagged-7Li-alpha" in msg        # what does carry the alignment
    assert "590952.42641509" in msg         # the measured zero, with its config


def test_make_plan_still_builds_all_four_at_spin_1_and_helicity_flip_at_32():
    for name in ("tensor-thirds", "transverse-tensor", "tensor-flip",
                 "helicity-flip"):
        plan = lg.make_plan(name, j=1.0, pz=0.7, pzz=0.6, pe=0.7)
        assert all(c.j == 1.0 for c in plan.categories)
    plan = lg.make_plan("helicity-flip", j=1.5, pz=0.7, pe=0.7)
    assert [c.j for c in plan.categories] == [1.5, 1.5]


# ------------------------------------------------------------------- G6

def _sigma(**kw):
    cfg = _cfg(**kw)
    p = _l.Pipeline(cfg, li7_tensor_plan())
    return p, np.asarray(p.sigma_per_category_pb())


@pytest.mark.skipif(not _l.HAVE_LHAPDF, reason="no LHAPDF tier in this build")
def test_the_sf_selectors_reach_the_7li_inclusive_kernel():
    """D2's selectors ARE the nucleon input a ⁷Li b1 would fold against.

    That is the whole reason the b1 is deferred rather than shipped: the
    unpolarised backend decides the SIGN of the α-t estimate at four of six
    x points (`phase_D_li7_rank2.md` §5.4), and this is the flag that
    chooses it.  Here it moves the ⁷Li rate and is recorded -- while the
    tensor sector stays exactly zero under every setting, because nothing
    fills the slot the backend would feed.
    """
    p0, s0 = _sigma()
    p1, s1 = _sigma(unpol_sf="ct18nlo")
    assert p1.sigma_pb() < 0.85 * p0.sigma_pb()      # measured ratio 0.795936
    assert s0[0] == s0[1] and s1[0] == s1[1]          # still exactly zero
    meta = p1.generate(0, False)["meta"]
    assert meta["unpol_sf"] == "ct18nlo"
    assert meta["b1_model"] == NONE_LABEL             # and still no backend


@pytest.mark.skipif(not _l.HAVE_LHAPDF, reason="no LHAPDF tier in this build")
def test_the_sf_selectors_reach_the_tagged_7li_kernel():
    def run(**kw):
        cfg = _cfg(channel="tagged-7Li-alpha", **kw)
        return _l.Pipeline(cfg, lg.make_plan("helicity-flip", j=1.5, pz=0.7,
                                             pe=0.7))
    p0, p1 = run(), run(unpol_sf="ct18nlo")
    assert p1.sigma_pb() < 0.85 * p0.sigma_pb()      # measured ratio 0.792543
    assert p1.generate(0, False)["meta"]["unpol_sf"] == "ct18nlo"


# ------------------------------------------------------------------- G7

def test_a_half_assigned_cluster_wave_raises_instead_of_crashing():
    """Defect F4: three lines of documented Python used to SEGFAULT.

        w = _lipolgen.ClusterPartialWave(); w.k = [0.1, 0.2, 0.3]

    `rebuild()` built `CubicSpline(k, phi)` with `phi` still empty and the
    constructor indexed `y_[i + 1]` off the end.  `rebuild()` now DEFERS, so
    either assignment order works, and every path that would read the
    mismatched pair throws.
    """
    w = _l.ClusterPartialWave()
    w.k = [0.1, 0.2, 0.3]
    for call in (lambda: w(0.15), lambda: w.norm2()):
        with pytest.raises(RuntimeError, match="half-assigned wave"):
            call()
    w.phi = [1.0, 2.0, 1.0]
    assert w(0.15) == pytest.approx(1.6875)
    # ... and the other order, which is the one a reader would write first
    v = _l.ClusterPartialWave()
    v.phi = [1.0, 2.0, 1.0]
    v.k = [0.1, 0.2, 0.3]
    assert v(0.15) == pytest.approx(1.6875)


def test_spin32_refusal_names_the_reachable_domain():
    """Defect F3, with the domain DERIVED and checked, not asserted.

    The four populations are fixed uniquely by (1, pz, t, o), so this is a
    domain and not a solver failure:

        p(+3/2) ± p(-3/2) = (1 + t)/2 , 0.9 pz + 0.1 o
        p(+1/2) ± p(-1/2) = (1 - t)/2 , 0.3 (pz - o)

    whence |0.9 pz + 0.1 o| <= (1 + t)/2 and |0.3 (pz - o)| <= (1 - t)/2, and
    at o = 0, 1.8|pz| - 1 <= t <= 1 - 0.6|pz| (so |pz| <= 5/6).
    """
    with pytest.raises(RuntimeError) as e:
        _l.spin32_populations(0.7, 0.6, 0.0)
    msg = str(e.value)
    assert "unphysical (pz = 0.7, t = 0.6, o = 0)" in msg
    assert "0.26 <= t <= 0.58" in msg
    assert "|pz| <= 5/6" in msg
    # the edge is where the message says it is
    _l.spin32_populations(0.7, 0.58, 0.0)
    # ... and the stated domain IS the function's, on a grid
    for pz in np.linspace(-1.0, 1.0, 21):
        for t in np.linspace(-1.0, 1.0, 21):
            inside = (1.8 * abs(pz) - 1.0 <= t + 1e-12
                      and t <= 1.0 - 0.6 * abs(pz) + 1e-12)
            try:
                _l.spin32_populations(float(pz), float(t), 0.0)
                ok = True
            except RuntimeError:
                ok = False
            assert ok == inside, (pz, t)
