# SPDX-License-Identifier: GPL-3.0-or-later
"""The 6Li b1 tables as a BAND over --b1-unpol x --unpol-sf.

`run_2026-09-02/phase_D_numbers.md` published the four-term 6Li convolution
tables on the DEFAULT `ToyF2` unpolarised input, and
`docs/OPEN_ITEMS_SOLUTIONS.md` sec. 10 has said since 2026-09-03 that they
"must be REGENERATED under --b1-unpol mstw before they are quoted as
physics".  They were regenerated on 2026-09-06 as the band over BOTH
selectors that reach them:
`docs/open_items/run_2026-09-06/phase_A_li6_tables.md`.

WHAT THIS FILE GATES:

  B1  ONE ROW of the band -- x = 0.30, Q^2 = 2.5, all three --b1-unpol
      backends -- pinned at rtol 1e-12.  This is table A4.3's x = 0.30 row
      and the `toy` entry is ALSO `phase_D_numbers.md`'s own published value,
      so a drift in either document fails here.
  B2  the band's two STRUCTURAL facts, which are why it is three wide and not
      seven:
        (a) terms (1) and (3) are BIT-IDENTICAL across the --b1-unpol axis
            (they convolve b1_d, not F1), so the axis moves the orbital
            sector (2d)+(2a) and nothing else;
        (b) the --unpol-sf axis is EXACTLY FLAT on x*b1 through the pipeline
            -- `Li6ConvolutionB1::b1` ignores its `f1` argument -- and the two
            cells where the kernel's own object WOULD leak in are refused by
            `PipelineConfig::validate()` BY NAME.
  B3  the mandatory +-100 % band is still exactly linear on every backend
      (scale 2 == 2 x scale 1 as doubles), and 0 is exactly 0.

NOTHING HERE IS A DEFAULT.  `Li6ConvolutionOptions.unpol = None` is the
shipped default and resolves to `ToyF2`; every non-`toy` row is opt-in and
skips loudly in a build without the matching optional tier, exactly as the
A = 2 gate's own MSTW rows do.

THE +-100 % A > 2 BAND RIDES WITH EVERY NUMBER BELOW.  It comes from
Q(6Li) = -0.0818(17) fm^2 against Q_d = +0.2859(3) fm^2, NOT from the A = 2
gate, and the gate says nothing about the alpha-d step.  These pins are the
CENTRES of {0, 1, 2} x b1 bands and are never quotable alone.

Run: `python -m pytest python/tests/test_li6_unpol_band.py -q`
"""

import pytest

import lipolgen as lg
from lipolgen import _lipolgen as _l

Q2 = 2.5
XS = (0.05, 0.10, 0.20, 0.30, 0.50)

needs_pythia = pytest.mark.skipif(not _l.HAVE_PYTHIA8,
                                  reason="MstwSF needs the PYTHIA 8 tier")
needs_lhapdf = pytest.mark.skipif(not _l.HAVE_LHAPDF,
                                  reason="LhapdfSF needs the LHAPDF tier")

# phase_A_li6_tables.md sec. A4.3, the x = 0.30 row: x*b1 per nucleon, 6Li,
# Q^2 = 2.5 GeV^2, DEFAULT options (raw digitized CDKS theory-1 b1_d,
# r_sigma_lt, kappa = 1, norm_target = VMC_N_ALPHA_D_LI6, w_CG = 1/10,
# w_alphad = 1, the 2400/2400/3200 y grid).  The document prints seven
# digits; pinned here to the full double, because rtol 1e-12 is what stops a
# silent drift.  The `toy` entry is phase_D_numbers.md's own -4.1732e-05.
BAND_X030 = {
    "toy":     -4.1732477875044357e-05,
    "ct18nlo": -5.1699559085229882e-05,
    "mstw":    -5.3249012164690352e-05,
}


def unpol_object(name):
    """The object `--b1-unpol NAME` installs, or None for the default."""
    if name == "toy":
        return None
    if name == "ct18nlo":
        lg.lhapdf_quiet()
        return lg.LhapdfSF("CT18NLO")
    if name == "mstw":
        return lg.MstwSF()
    raise AssertionError(name)


def model(name, **kw):
    o = lg.Li6ConvolutionOptions()
    u = unpol_object(name)
    if u is not None:
        o.unpol = u
    for k, v in kw.items():
        setattr(o, k, v)
    return lg.Li6ConvolutionB1(o)


# ------------------------------------------------------------------ B1

@pytest.mark.parametrize("name", [
    "toy",
    pytest.param("ct18nlo", marks=needs_lhapdf),
    pytest.param("mstw", marks=needs_pythia),
])
def test_band_row_x030(name):
    """ONE row of the band, pinned: x*b1 at x = 0.30, Q^2 = 2.5."""
    got = 0.30 * model(name).b1(0.30, Q2, 0.0)
    assert got == pytest.approx(BAND_X030[name], rel=1e-12)


def test_band_row_x030_is_a_band_and_not_a_row():
    """The three pins are a BAND: the spread is a factor, and it is stated.

    Guards the reading, not just the numbers: at x = 0.30 the unpolarised
    backend moves |x*b1| by x1.275961 (max/min), which is INSIDE the
    mandatory +-100 % band there -- unlike x = 0.05 and x = 0.10, where it is
    NOT (x2.318 and x2.221, both > 2).  phase_A_li6_tables.md sec. A4.3.
    """
    if not (_l.HAVE_LHAPDF and _l.HAVE_PYTHIA8):
        pytest.skip("the band needs both optional tiers")
    v = [abs(BAND_X030[n]) for n in ("toy", "ct18nlo", "mstw")]
    assert max(v) / min(v) == pytest.approx(1.2759610, rel=1e-6)
    assert max(v) / min(v) < 2.0          # contained at x = 0.30 ...
    lo = [abs(0.05 * model(n).b1(0.05, Q2, 0.0))
          for n in ("toy", "ct18nlo", "mstw")]
    assert max(lo) / min(lo) > 2.0        # ... and NOT contained at x = 0.05


# ------------------------------------------------------------------ B2 (a)

@pytest.mark.skipif(not (_l.HAVE_LHAPDF and _l.HAVE_PYTHIA8),
                    reason="needs both optional tiers")
def test_axis_moves_only_the_orbital_terms():
    """Terms (1) and (3) are BIT-IDENTICAL across the --b1-unpol axis.

    They convolve b1_d; only (2d) and (2a) fold against F1.  Bit equality,
    not a tolerance: there is no arithmetic in between.
    """
    ref = model("toy")
    r1 = [ref.b1_embedded_s(x, Q2) for x in XS]
    r3 = [ref.b1_cg_dwave(x, Q2) for x in XS]
    for name in ("ct18nlo", "mstw"):
        b = model(name)
        assert [b.b1_embedded_s(x, Q2) for x in XS] == r1
        assert [b.b1_cg_dwave(x, Q2) for x in XS] == r3
        # ... and the orbital sector is what DOES move.
        assert [b.b1_alpha_d_dwave(x, Q2) for x in XS] != \
               [ref.b1_alpha_d_dwave(x, Q2) for x in XS]
    # Same fact from the other side: w_alpha_d_dwave = 0 deletes the axis.
    r0 = [model("toy", w_alpha_d_dwave=0.0).b1(x, Q2, 0.0) for x in XS]
    for name in ("ct18nlo", "mstw"):
        assert [model(name, w_alpha_d_dwave=0.0).b1(x, Q2, 0.0)
                for x in XS] == r0


# ------------------------------------------------------------------ B2 (b)

def _pipeline(b1_unpol, unpol_sf):
    """The recipe of phase_A_li6_tables.md sec. A4.1, coarse grid.

    `scenario.x_max = 0.95` is a REFUSAL and not a taste: at the shipped 1.0
    every cell of this matrix -- the toy/toy default included -- is refused by
    `InclusiveSampler` ("negative phi-averaged density for m = 1 at
    x = 0.955"), which is the pre-existing refused-base condition and has
    nothing to do with either selector.
    """
    cfg = lg.make_config(isotope="6Li", channel="inclusive", config=1,
                         events=2000, seed=99, b1_model="li6-convolution",
                         b1_unpol=b1_unpol, unpol_sf=unpol_sf)
    g = cfg.grid
    g.nx, g.nq2 = 12, 6                  # SPEED only
    cfg.grid = g
    sc = cfg.scenario
    sc.x_max = 0.95
    cfg.scenario = sc
    return lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))


@pytest.mark.skipif(not (_l.HAVE_LHAPDF and _l.HAVE_PYTHIA8),
                    reason="needs both optional tiers")
@pytest.mark.parametrize("b1_unpol", ["ct18nlo", "mstw"])
def test_unpol_sf_axis_is_flat_on_b1(b1_unpol):
    """--unpol-sf moves the kernel's F1 and NOT b1: bit-identical, and the
    pipeline's b1 is bit-identical to the direct Li6ConvolutionB1's."""
    ref = _pipeline(b1_unpol, "toy").dis_sampler.kernel
    ref_b1 = [ref.tables(x, Q2).b1 for x in XS]
    for sf in ("ct18nlo", "mstw"):
        k = _pipeline(b1_unpol, sf).dis_sampler.kernel
        assert [k.tables(x, Q2).b1 for x in XS] == ref_b1
        # ... while the kernel's own F1 DOES move, which is the other half.
        assert k.tables(0.10, Q2).f1 != ref.tables(0.10, Q2).f1
    direct = model(b1_unpol)
    assert ref_b1 == [direct.b1(x, Q2, 0.0) for x in XS]


@pytest.mark.skipif(not (_l.HAVE_LHAPDF and _l.HAVE_PYTHIA8),
                    reason="needs both optional tiers")
@pytest.mark.parametrize("sf", ["ct18nlo", "mstw"])
def test_toy_b1_under_a_non_toy_unpol_sf_is_refused_by_name(sf):
    """The two cells that make the 3 x 3 matrix seven-wide, not nine.

    They are a PROVENANCE refusal: `toy` there means "the kernel's own
    UnpolSF, shared as one object", and off `--unpol-sf toy` that object is
    not ToyF2 -- so meta["b1_unpol"] would record "toy" for something else.
    """
    with pytest.raises(Exception) as e:
        _pipeline("toy", sf)
    msg = str(e.value)
    assert "b1_unpol = toy" in msg and ("unpol_sf = %s" % sf) in msg
    assert 'meta["b1_unpol"]' in msg


# ------------------------------------------------------------------ B3

@pytest.mark.parametrize("name", [
    "toy",
    pytest.param("ct18nlo", marks=needs_lhapdf),
    pytest.param("mstw", marks=needs_pythia),
])
def test_mandatory_band_is_exactly_linear_on_every_backend(name):
    """`{0, 1, 2} x b1` -- 0 is exactly 0 and 2 is bit-for-bit twice 1."""
    b = model(name)
    at = [0.30 * b.banded(s).b1(0.30, Q2, 0.0) for s in (0.0, 1.0, 2.0)]
    assert at[0] == 0.0
    assert at[1] == pytest.approx(BAND_X030[name], rel=1e-12)
    assert at[2] == 2.0 * at[1]
