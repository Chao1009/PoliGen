// A1 -- the MSTW2008 LO unpolarized backend (`include/lipolgen/mstw_sf.hpp`,
// `src/pythia/mstw_sf.cpp`).
//
// WHY THIS FILE EXISTS.  Every document in the tree records MSTW2008 LO as
// "not installed" and treats CT18NLO as a stand-in for the PDF CDKS actually
// used.  It is installed -- as a PYTHIA grid rather than an LHAPDF set --
// and `MstwSF` reads it.  What is pinned here:
//
//   1. the backend constructs, over a HEAP-allocated Pythia8::MSTWpdf
//      (sizeof == 7 988 336 bytes: a stack-local instance segfaults during
//      construction, measured under gdb);
//   2. F2p and F2n at Q2 = 2.5 GeV^2 at seven x, against an INDEPENDENT
//      probe -- a standalone translation unit that linked libpythia8 and
//      called Pythia8::MSTWpdf::xf directly, with no LiPolGen code in the
//      path (docs/open_items/run_2026-09-03/PLAN.md:54);
//   3. the 0 < x < 1 guard, verbatim from LhapdfSF;
//   4. F2n/F2p at the fixed Q2 = 10 the single-argument interface carries;
//   5. (only when LHAPDF is also built) the MSTW/CT18NLO ratio on F2p,
//      which is the number the 6Li b1 gate turns on -- it must be computed
//      through BOTH wrappers, since the point of copying `f2_from_weights`
//      verbatim is that the ratio is a PDF comparison, not a convention
//      comparison.
//
// The whole file compiles out when PYTHIA 8 is not configured (the
// `#ifdef LIPOLGEN_HAVE_PYTHIA8` pattern of tests/test_t2.cpp), and the
// grid-file presence is re-checked at run time so a PYTHIA build without
// its pdfdata tree reports a skip instead of a failure.

#ifdef LIPOLGEN_HAVE_PYTHIA8

#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>

#include "check_close.hpp"
#include "doctest.h"

#include "lipolgen/mstw_sf.hpp"
#include "lipolgen/sf.hpp"

#ifdef LIPOLGEN_HAVE_LHAPDF
#include "lipolgen/lhapdf_sf.hpp"
#endif

using namespace lipolgen;

namespace {

/// The MSTW2008 LO central grid, i.e. what `MstwSF{}` (iFit = 3) reads.
bool mstw_grid_present() {
  const std::string& d = pythia8_pdfdata_dir();
  if (d.empty()) return false;
  return std::filesystem::exists(std::filesystem::path(d) /
                                 "mstw2008lo.00.dat");
}

// The independent probe's numbers, at %.17g so they round-trip exactly.  The
// probe is the one described in the file header: it builds the SAME 5-flavour
// e_q^2 combination out of raw Pythia8::MSTWpdf::xf calls, so agreement here
// is "LiPolGen's wrapper adds nothing", which is exactly what must be true
// for MSTW and CT18NLO to be comparable.
struct Row { double x, f2p, f2n; };

// Q2 = 2.5 GeV^2 -- the gate's own working scale.
constexpr Row kQ2p5[] = {
    {0.01, 0.39012889888888885, 0.3779788988888888},
    {0.10, 0.3683144888888889, 0.32580782222222227},
    {0.30, 0.32790599999999992, 0.21192599999999995},
    {0.50, 0.18034921000000001, 0.089512543333333333},
    {0.70, 0.049053195137777769, 0.01764119513777777},
    {0.80, 0.015238453801333329, 0.0046093538013333331},
    {0.95, 0.00020443430577674784, 5.1962239110081217e-05},
};

// Q2 = 10 GeV^2 -- the scale f2n_over_f2p is nailed to.
constexpr Row kQ10[] = {
    {0.01, 0.54330254088888885, 0.52937587422222221},
    {0.10, 0.39035346666666659, 0.34016679999999999},
    {0.30, 0.28382388888888888, 0.17948722222222221},
    {0.50, 0.13249083633333331, 0.064338169666666653},
    {0.70, 0.030178569177777777, 0.010666569177777778},
    {0.80, 0.0083307632742444446, 0.002497113274244444},
    {0.95, 7.872514577744692e-05, 2.0044952444113585e-05},
};

// Same interpolator, same grid, same arithmetic order: the wrapper must be
// bit-for-bit the probe.  1e-14 is a floating-point-reassociation allowance,
// not a physics tolerance.
constexpr double kExact = 1e-14;

}  // namespace

TEST_CASE("MstwSF constructs off PYTHIA's own MSTW2008 LO grid") {
  if (!mstw_grid_present()) {
    MESSAGE("SKIP: no mstw2008lo.00.dat under '" << pythia8_pdfdata_dir()
            << "' -- PYTHIA 8 built without its pdfdata tree");
    return;
  }
  // Not a stack local.  sizeof(Pythia8::MSTWpdf) is 7 988 336 bytes, so
  // `MstwSF` must (and does) keep it behind a unique_ptr; the class itself
  // is small enough to live anywhere.
  CHECK(sizeof(MstwSF) <= 64);

  std::unique_ptr<MstwSF> sf;
  CHECK_NOTHROW(sf = std::make_unique<MstwSF>());
  REQUIRE(sf != nullptr);
  CHECK(sf->i_fit() == 3);  // MSTW 2008 LO, central member
  CHECK(!sf->pdfdata_path().empty());

  // It IS a UnpolSF: the F1/FL machinery of the base class works on it.
  const UnpolSF& base = *sf;
  const double f2 = base.f2p(0.3, 2.5);
  CHECK(f2 > 0.0);
  CHECK_CLOSE(base.f1p(0.3, 2.5), f2 / (2.0 * 0.3 * (1.0 + base.r(0.3, 2.5))),
              1e-15);

  // A grid that does not exist is a setup error, not a silent zero.  PYTHIA
  // reports the miss by printing "Error in MSTWpdf::init: did not find data
  // file" on stdout and clearing isSetup(); the line below is that expected
  // printout, not a test failure.
  MESSAGE("expect one 'Error in MSTWpdf::init' line from PYTHIA next:");
  CHECK_THROWS_AS(MstwSF(3, "/nonexistent/pdfdata"), std::runtime_error);
}

TEST_CASE("MstwSF::f2p / f2n reproduce a raw Pythia8::MSTWpdf probe") {
  if (!mstw_grid_present()) {
    MESSAGE("SKIP: no mstw2008lo.00.dat under '" << pythia8_pdfdata_dir()
            << "'");
    return;
  }
  MstwSF sf;

  for (const Row& r : kQ2p5) {
    CAPTURE(r.x);
    CHECK_CLOSE(sf.f2p(r.x, 2.5), r.f2p, kExact);
    CHECK_CLOSE(sf.f2n(r.x, 2.5), r.f2n, kExact);
  }
  for (const Row& r : kQ10) {
    CAPTURE(r.x);
    CHECK_CLOSE(sf.f2p(r.x, 10.0), r.f2p, kExact);
    CHECK_CLOSE(sf.f2n(r.x, 10.0), r.f2n, kExact);
  }

  // The same numbers as the run-2026-09-03 inventory wrote them down, to the
  // 5 significant figures it quoted: F2p(Q2 = 2.5) = 3.9013e-1, 3.6831e-1,
  // 3.2791e-1, 1.8035e-1, 4.9053e-2, 1.5238e-2, 2.0443e-4.  Redundant with
  // the table above by construction, and kept because THESE are the digits a
  // reader of the plan can check by eye.
  CHECK_CLOSE(sf.f2p(0.01, 2.5), 3.9013e-1, 1e-4);
  CHECK_CLOSE(sf.f2p(0.10, 2.5), 3.6831e-1, 1e-4);
  CHECK_CLOSE(sf.f2p(0.30, 2.5), 3.2791e-1, 1e-4);
  CHECK_CLOSE(sf.f2p(0.50, 2.5), 1.8035e-1, 1e-4);
  CHECK_CLOSE(sf.f2p(0.70, 2.5), 4.9053e-2, 1e-4);
  CHECK_CLOSE(sf.f2p(0.80, 2.5), 1.5238e-2, 1e-4);
  CHECK_CLOSE(sf.f2p(0.95, 2.5), 2.0443e-4, 1e-4);
}

TEST_CASE("MstwSF carries LhapdfSF's 0 < x < 1 guard and Q2 = 10 ratio") {
  if (!mstw_grid_present()) {
    MESSAGE("SKIP: no mstw2008lo.00.dat under '" << pythia8_pdfdata_dir()
            << "'");
    return;
  }
  MstwSF sf;

  // Verbatim from f2_from_weights: outside the OPEN unit interval the grid
  // is never touched and the answer is exactly zero.
  for (double x : {0.0, 1.0, -0.1, 1.5}) {
    CAPTURE(x);
    CHECK(sf.f2p(x, 2.5) == 0.0);
    CHECK(sf.f2n(x, 2.5) == 0.0);
  }
  // ... and there f2n_over_f2p falls back to 1.0, as PartonF2's
  // np.where(f2p > 0, ..., 1.0) does.
  CHECK(sf.f2n_over_f2p(0.0) == 1.0);
  CHECK(sf.f2n_over_f2p(1.0) == 1.0);

  // f2n_over_f2p is evaluated at the fixed Q2 = 10 the single-argument
  // interface carries (identical to LhapdfSF::f2n_over_f2p).
  for (const Row& r : kQ10) {
    CAPTURE(r.x);
    CHECK_CLOSE(sf.f2n_over_f2p(r.x), r.f2n / r.f2p, kExact);
  }
  // Sanity on the physics: the valence-driven fall of F2n/F2p with x, from
  // ~0.97 near x = 0.01 down to ~0.25 at x = 0.95.
  CHECK(sf.f2n_over_f2p(0.01) > sf.f2n_over_f2p(0.5));
  CHECK(sf.f2n_over_f2p(0.5) > sf.f2n_over_f2p(0.95));
}

#ifdef LIPOLGEN_HAVE_LHAPDF
TEST_CASE("MSTW2008 LO / CT18NLO on F2p -- the ratio the b1 gate turns on") {
  if (!mstw_grid_present()) {
    MESSAGE("SKIP: no mstw2008lo.00.dat under '" << pythia8_pdfdata_dir()
            << "'");
    return;
  }
  MstwSF mstw;
  LhapdfSF ct18("CT18NLO", 0);

  // The inventory's measured ratios at Q2 = 2.5 GeV^2 (PLAN.md:54, phase A):
  // 0.93 / 0.89 / 0.98 / 1.17 / 1.42 / 1.56 at x = 0.01 ... 0.80, quoted to
  // two decimals.  What this backend measures, in full: 0.932231 / 0.884935
  // / 0.983863 / 1.168327 / 1.417586 / 1.556314.  Five of the six round to
  // the quoted digits exactly; x = 0.10 lands on the rounding BOUNDARY --
  // 0.884935 is 0.88 under round-half-even and 0.89 under the round-half-up
  // of a printf("%.2f") applied to a value the probe evaluated as 0.885.
  // That is a display artifact of 6.5e-5 in the ratio, not a physics
  // difference: F2p itself reproduces all seven of the probe's 5-digit
  // values to the last quoted digit (previous test case), and CT18NLO's own
  // two independent interpolators bracket it -- LHAPDF's C++
  // LogBicubicInterpolator gives F2p(0.10, 2.5) = 0.416205 and the pure
  // Python `parton` package 0.416313, a 2.6e-4 spread that straddles
  // nothing here.  So the bound below is ONE unit in the last quoted digit,
  // and the full-precision pins that follow carry the real tolerance.
  struct RatioRow { double x, want; };
  const RatioRow rows[] = {{0.01, 0.93}, {0.10, 0.89}, {0.30, 0.98},
                           {0.50, 1.17}, {0.70, 1.42}, {0.80, 1.56}};
  for (const RatioRow& r : rows) {
    CAPTURE(r.x);
    const double got = mstw.f2p(r.x, 2.5) / ct18.f2p(r.x, 2.5);
    CHECK_CLOSE_AT(got, r.want, 0.0, 0.01);
  }

  // Full precision of the same six, so a future PDF or interpolator change
  // is visible rather than absorbed by the two-decimal bound above.  Both
  // sides are built by the SAME `f2_from_weights` body (copied verbatim
  // between src/lhapdf/lhapdf_sf.cpp and src/pythia/mstw_sf.cpp), so this
  // is a PDF ratio and nothing else.  1e-9 rather than 1e-14: the CT18NLO
  // half comes from LHAPDF's own interpolator, which this repository does
  // not pin bit-for-bit anywhere.
  const RatioRow exact[] = {{0.01, 0.93223106571500802},
                            {0.10, 0.88493472787478378},
                            {0.30, 0.98386258869146770},
                            {0.50, 1.16832690499277360},
                            {0.70, 1.41758597117899640},
                            {0.80, 1.55631406269775230}};
  for (const RatioRow& r : exact) {
    CAPTURE(r.x);
    CHECK_CLOSE(mstw.f2p(r.x, 2.5) / ct18.f2p(r.x, 2.5), r.want, 1e-9);
  }
  // The G3b peak sits at x = 0.736, where MSTW is ~1.45x CT18NLO on F2p --
  // the number that moves the gate's measured 0.7193.
  MESSAGE("MSTW/CT18NLO F2p at x = 0.736, Q2 = 2.5: "
          << mstw.f2p(0.736, 2.5) / ct18.f2p(0.736, 2.5));
}
#endif  // LIPOLGEN_HAVE_LHAPDF

#endif  // LIPOLGEN_HAVE_PYTHIA8
