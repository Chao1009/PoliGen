"""MSTW2008 LO through the binding -- the same numbers the doctest pins.

Companion to `tests/test_mstw_sf.cpp`.  MSTW2008 LO is the PDF CDKS computed
their b1 predictions with; every document in the tree records it as "not
installed", which is true only of LHAPDF's set store -- PYTHIA ships the
central grid as `<pdfdata>/mstw2008lo.00.dat` and `Pythia8::MSTWpdf` to read
it.  `lg.MstwSF` is that reader behind the ordinary `UnpolSF` surface, built
from the SAME `f2_from_weights` body as `lg.LhapdfSF` (copied verbatim
between `src/lhapdf/lhapdf_sf.cpp` and `src/pythia/mstw_sf.cpp`), so an
MSTW/CT18NLO ratio taken through the two is a PDF comparison and not a
convention comparison.

Skips wholesale when the PYTHIA 8 tier is not built, and per-test when the
tier is built without its `pdfdata` tree.
"""

import os

import pytest

import lipolgen as lg

pytestmark = pytest.mark.skipif(not lg._lipolgen.HAVE_PYTHIA8,
                                reason="no PYTHIA 8 tier")

# F2p / F2n at Q2 = 2.5 GeV^2, from the independent probe described in
# tests/test_mstw_sf.cpp: a standalone translation unit that linked
# libpythia8 and called Pythia8::MSTWpdf::xf directly, no LiPolGen in the
# path.  Full precision so the binding is pinned, not merely sampled.
Q2P5 = [
    (0.01, 0.39012889888888885, 0.3779788988888888),
    (0.10, 0.3683144888888889, 0.32580782222222227),
    (0.30, 0.32790599999999992, 0.21192599999999995),
    (0.50, 0.18034921000000001, 0.089512543333333333),
    (0.70, 0.049053195137777769, 0.01764119513777777),
    (0.80, 0.015238453801333329, 0.0046093538013333331),
    (0.95, 0.00020443430577674784, 5.1962239110081217e-05),
]

# ... and at Q2 = 10, the scale the single-argument f2n_over_f2p is nailed to.
Q10 = [
    (0.01, 0.54330254088888885, 0.52937587422222221),
    (0.10, 0.39035346666666659, 0.34016679999999999),
    (0.30, 0.28382388888888888, 0.17948722222222221),
    (0.50, 0.13249083633333331, 0.064338169666666653),
    (0.70, 0.030178569177777777, 0.010666569177777778),
    (0.80, 0.0083307632742444446, 0.002497113274244444),
    (0.95, 7.872514577744692e-05, 2.0044952444113585e-05),
]


def _grid_present():
    d = lg._lipolgen.pythia8_pdfdata_dir()
    return bool(d) and os.path.isfile(os.path.join(d, "mstw2008lo.00.dat"))


needs_grid = pytest.mark.skipif(
    not _grid_present(),
    reason="PYTHIA 8 built without its pdfdata/mstw2008lo.00.dat")


def test_mstw_is_bound_beside_lhapdf():
    assert hasattr(lg, "MstwSF")
    assert issubclass(lg.MstwSF, lg.UnpolSF)


@needs_grid
def test_construction_and_provenance():
    sf = lg.MstwSF()
    assert sf.i_fit == 3                       # MSTW 2008 LO, central member
    assert sf.pdfdata_path                     # non-empty
    assert os.path.isfile(os.path.join(sf.pdfdata_path, "mstw2008lo.00.dat"))
    # A directory with no grid in it is a setup error, not a silent zero.
    # (PYTHIA prints "Error in MSTWpdf::init: did not find data file" first.)
    with pytest.raises(RuntimeError):
        lg.MstwSF(3, "/nonexistent/pdfdata")


@needs_grid
@pytest.mark.parametrize("x,f2p,f2n", Q2P5)
def test_f2_at_q2_2p5(x, f2p, f2n):
    sf = lg.MstwSF()
    assert sf.f2p(x, 2.5) == pytest.approx(f2p, rel=1e-14)
    assert sf.f2n(x, 2.5) == pytest.approx(f2n, rel=1e-14)


@needs_grid
@pytest.mark.parametrize("x,f2p,f2n", Q10)
def test_f2_and_ratio_at_q2_10(x, f2p, f2n):
    sf = lg.MstwSF()
    assert sf.f2p(x, 10.0) == pytest.approx(f2p, rel=1e-14)
    assert sf.f2n(x, 10.0) == pytest.approx(f2n, rel=1e-14)
    # f2n_over_f2p carries LhapdfSF's fixed Q2 = 10.
    assert sf.f2n_over_f2p(x) == pytest.approx(f2n / f2p, rel=1e-14)


@needs_grid
def test_the_quoted_five_digit_numbers():
    """The digits docs/open_items/run_2026-09-03/PLAN.md quotes, by eye."""
    sf = lg.MstwSF()
    want = [3.9013e-1, 3.6831e-1, 3.2791e-1, 1.8035e-1,
            4.9053e-2, 1.5238e-2, 2.0443e-4]
    got = [sf.f2p(x, 2.5) for x, _, _ in Q2P5]
    assert got == pytest.approx(want, rel=1e-4)


@needs_grid
def test_open_unit_interval_guard():
    # Verbatim from f2_from_weights: outside 0 < x < 1 the grid is never
    # touched and the answer is exactly zero; f2n_over_f2p then falls back
    # to 1.0, as PartonF2's np.where(f2p > 0, ..., 1.0) does.
    sf = lg.MstwSF()
    for x in (0.0, 1.0, -0.1, 1.5):
        assert sf.f2p(x, 2.5) == 0.0
        assert sf.f2n(x, 2.5) == 0.0
    assert sf.f2n_over_f2p(0.0) == 1.0
    assert sf.f2n_over_f2p(1.0) == 1.0


@needs_grid
def test_usable_as_an_unpolsf_backend():
    # The point of the class: it plugs into everything that takes a UnpolSF.
    sf = lg.MstwSF()
    f2 = sf.f2p(0.3, 2.5)
    assert sf.f1p(0.3, 2.5) == pytest.approx(
        f2 / (2.0 * 0.3 * (1.0 + sf.r(0.3, 2.5))), rel=1e-15)
    li6 = lg.li6()
    nuc = lg.NuclearF2(li6, sf)
    # F2A = Z F2p + N F2n with no EMC ratio: 3 * (F2p + F2n) for 6Li.
    assert nuc.f2a(0.3, 2.5) == pytest.approx(
        3.0 * (sf.f2p(0.3, 2.5) + sf.f2n(0.3, 2.5)), rel=1e-14)


@needs_grid
@pytest.mark.skipif(not lg._lipolgen.HAVE_LHAPDF, reason="no LHAPDF backend")
@pytest.mark.parametrize("x,want", [(0.01, 0.93223106571500802),
                                    (0.10, 0.88493472787478378),
                                    (0.30, 0.98386258869146770),
                                    (0.50, 1.16832690499277360),
                                    (0.70, 1.41758597117899640),
                                    (0.80, 1.55631406269775230)])
def test_mstw_over_ct18nlo_on_f2p(x, want):
    """The ratio the 6Li b1 gate turns on, through both wrappers.

    PLAN.md:54 quotes 0.93 / 0.89 / 0.98 / 1.17 / 1.42 / 1.56 to two
    decimals.  Five of six round to exactly that; x = 0.10 measures 0.884935,
    which is on the rounding boundary (0.88 half-even, 0.89 half-up) -- a
    display artifact of 6.5e-5 in the ratio, not a physics difference.  F2p
    itself reproduces all seven of the probe's 5-digit values exactly, and
    the CT18NLO half is the one that is not pinned bit-for-bit anywhere in
    this repository, hence rel=1e-9 rather than 1e-14.
    """
    mstw = lg.MstwSF()
    ct18 = lg.LhapdfSF("CT18NLO", 0)
    assert mstw.f2p(x, 2.5) / ct18.f2p(x, 2.5) == pytest.approx(want, rel=1e-9)
