// Structure-function backends against the Python they are ported from:
// PolarizedLithiumSim/fastsim/polli_fastsim/{structure,polarized,delta_models}.py
// and their tests (test_structure_r.py, test_digitized_curves.py,
// test_delta_models.py, test_polarized_normalisation.py).
//
// Every literal in this file is a number produced by that Python and printed
// at %.17g, so the comparisons are genuine cross-language pins rather than
// self-consistency checks.

#include <cmath>
#include <memory>
#include <vector>

#include "check_close.hpp"
#include "lipolgen/beams.hpp"
#include "lipolgen/numerics.hpp"
#include "lipolgen/sf.hpp"

using namespace lipolgen;

namespace {
constexpr double kRtol = 1e-12;
}

TEST_CASE("np_interp reproduces numpy's clamped linear interpolation") {
  const double xp[4] = {0.0, 1.0, 2.0, 3.0};
  const double fp[4] = {0.0, 10.0, 20.0, 25.0};
  CHECK(np_interp(-1.0, xp, fp, 4) == 0.0);   // clamped below
  CHECK(np_interp(4.0, xp, fp, 4) == 25.0);   // clamped above
  CHECK(np_interp(2.0, xp, fp, 4) == 20.0);   // exact node
  CHECK(np_interp(3.0, xp, fp, 4) == 25.0);   // last node
  CHECK_CLOSE(np_interp(2.5, xp, fp, 4), 22.5, 1e-15);
}

TEST_CASE("the NumPy pairwise summation order is reproduced") {
  // 1/i over 300 terms: a left-to-right accumulation differs from NumPy's
  // eight-accumulator pairwise sum in the last bits, which is exactly the
  // difference a trapezoid integral compared at rtol 1e-12 must not carry.
  std::vector<double> v(300);
  for (std::size_t i = 0; i < v.size(); ++i) v[i] = 1.0 / static_cast<double>(i + 1);
  double naive = 0.0;
  for (double d : v) naive += d;
  const double pw = pairwise_sum(v);
  CHECK_CLOSE(pw, naive, 1e-14);        // same value to 14 digits
  CHECK_CLOSE(pw, 6.2826638802995038, 1e-16);  // ... and this is NumPy's
}

TEST_CASE("ToyF2 reproduces the Python F2p / F2n / R") {
  const ToyF2 f2;
  struct Row { double x, q2, f2p, f2n, r; };
  const Row rows[] = {
      {0.0224, 1.14, 0.4240130750444826, 0.41688965538373529, 0.17598748533437622},
      {0.056, 1.14, 0.38728987162222994, 0.37102369701409627, 0.17598748533437622},
      {0.141, 3.13, 0.32795748003196507, 0.29327597651858478, 0.16939582156973462},
      {0.30, 30.0, 0.20931444996438417, 0.16221869872239775, 0.11249999999999999},
  };
  for (const Row& row : rows) {
    CHECK_CLOSE(f2.f2p(row.x, row.q2), row.f2p, kRtol);
    CHECK_CLOSE(f2.f2n(row.x, row.q2), row.f2n, kRtol);
    CHECK_CLOSE(r_sigma_lt(row.x, row.q2), row.r, kRtol);
    // F1 = F2/(2x(1+R)) and FL = F2 R/(1+R) are the two R-derived slots
    CHECK_CLOSE(f2.f1p(row.x, row.q2),
                row.f2p / (2.0 * row.x * (1.0 + row.r)), kRtol);
    CHECK_CLOSE(f2.flp(row.x, row.q2), row.f2p * row.r / (1.0 + row.r), kRtol);
  }
  CHECK_CLOSE(f2.f2n_over_f2p(0.5), 1.0 - 0.75 * 0.5, kRtol);
  CHECK(f2.f2n_over_f2p(1.0) == 0.25);  // clipped
  CHECK(f2.f2n_over_f2p(0.0) == 1.0);
}

TEST_CASE("R1998 is the average of its three published forms") {
  const R1998Forms f = r1998_forms(0.1, 10.0);
  CHECK_CLOSE(f.a, 0.12367334583889994, kRtol);
  CHECK_CLOSE(f.b, 0.12624333788877903, kRtol);
  CHECK_CLOSE(f.c, 0.12972133778614703, kRtol);
  CHECK_CLOSE(r1998(0.1, 10.0), 0.12654600717127532, kRtol);
  CHECK_CLOSE(r1998_spread(0.1, 10.0), f.c - f.a, kRtol);
  CHECK_CLOSE(r1998_fit_error(0.1, 10.0), 0.009747863247863248, kRtol);
  // the fit's support is frozen at the boundary with clip=true
  CHECK(r1998(1e-4, 10.0, R1998Form::kAverage, true)
        == r1998(R1998_X_MIN, 10.0, R1998Form::kAverage, true));
  // over the whole support the forms stay in 0.012 <= R <= 0.45
  for (double lx = -2.3; lx <= -0.07; lx += 0.1) {
    for (double lq = -0.3; lq <= 2.11; lq += 0.1) {
      const double v = r1998(std::pow(10.0, lx), std::pow(10.0, lq));
      CHECK(v > 0.01);
      CHECK(v < 0.46);
    }
  }
}

TEST_CASE("NuclearF2 builds F2A = Z F2p + N F2n with the EMC hook") {
  const auto base = std::make_shared<const ToyF2>();
  const NuclearF2 nf(LI6(), base);
  const double x = 0.056, q2 = 1.14;
  CHECK_CLOSE(nf.f2a(x, q2),
              3 * base->f2p(x, q2) + 3 * base->f2n(x, q2), kRtol);
  CHECK_CLOSE(nf.f1a(x, q2),
              nf.f2a(x, q2) / (2.0 * x * (1.0 + r_sigma_lt(x, q2))), kRtol);
  // the EMC hook multiplies the isoscalar combination and nothing else
  const NuclearF2 nf_emc(LI6(), base, [](double xx) { return 1.0 - 0.1 * xx; });
  CHECK_CLOSE(nf_emc.f2a(x, q2), nf.f2a(x, q2) * (1.0 - 0.1 * x), kRtol);
  // ... and an explicit r_func moves F1 only
  const NuclearF2 nf_r(LI6(), base, nullptr,
                       [](double xx, double qq) { return r1998(xx, qq); });
  CHECK(nf_r.f2a(x, q2) == nf.f2a(x, q2));
  CHECK_CLOSE(nf_r.f1a(x, q2),
              nf.f2a(x, q2) / (2.0 * x * (1.0 + r1998(x, q2))), kRtol);
}

TEST_CASE("dsigma_dx_dq2 is the NC DIS double-differential cross section") {
  const ToyF2 f2;
  const double x = 0.056, q2 = 1.14, s = 3980.0000000000005;
  CHECK_CLOSE(dsigma_dx_dq2(x, q2, s, f2.f2p(x, q2)), 1379598.3747138162, kRtol);
  // an explicit FL overrides R entirely
  const double fl = 0.0;
  const double y = q2 / (s * x);
  CHECK_CLOSE(dsigma_dx_dq2(x, q2, s, f2.f2p(x, q2), &fl),
              4.0 * kPi * ALPHA_EM * ALPHA_EM / (x * q2 * q2)
                  * (1.0 - y + 0.5 * y * y) * f2.f2p(x, q2) * GEV2_TO_PB,
              kRtol);
}

TEST_CASE("ToyG1 and the nuclear g1A = Z P_p g1p + N P_n g1n") {
  const auto base = std::make_shared<const ToyF2>();
  const ToyG1 g1(base);
  const double x = 0.141, q2 = 3.13;
  const double r = r_sigma_lt(x, q2);
  CHECK_CLOSE(g1.a1p(x), std::pow(x, 0.7), kRtol);
  CHECK_CLOSE(g1.a1n(x), -0.07 * std::pow(1.0 - x, 2.0) + 0.8 * std::pow(x, 2.2),
              kRtol);
  CHECK_CLOSE(g1.g1p(x, q2),
              g1.a1p(x) * base->f2p(x, q2) / (2.0 * x * (1.0 + r)), kRtol);
  CHECK_CLOSE(g1.g1_nucleus(LI7(), x, q2),
              LI7().Z * LI7().eff_pol_p * g1.g1p(x, q2)
                  + LI7().N() * LI7().eff_pol_n * g1.g1n(x, q2), kRtol);
  // the per-nucleon storage convention: Z*P_p and N*P_n are the whole-nucleus
  // VMC sums 0.866 and -0.037 (test_polarized_normalisation.py)
  CHECK_CLOSE(g1.g1_nucleus(LI7(), x, q2),
              0.866 * g1.g1p(x, q2) - 0.037 * g1.g1n(x, q2), 1e-13);
  // the medium ratio DR(x) multiplies the whole thing
  CHECK_CLOSE(g1.g1_nucleus(LI6(), x, q2, [](double) { return 0.9; }),
              0.9 * g1.g1_nucleus(LI6(), x, q2), kRtol);
}

TEST_CASE("g2_ww on a power law matches its analytic Wandzura-Wilczek form") {
  // g1 = x^a -> g2_WW = -x^a + (1 - x^a)/a
  const double a = 0.7;
  const auto g1 = [a](double x, double) { return std::pow(x, a); };
  for (double x : {0.01, 0.1, 0.4, 0.8}) {
    const double got = g2_ww(g1, x, 10.0, 400);
    const double want = -std::pow(x, a) + (1.0 - std::pow(x, a)) / a;
    CHECK_CLOSE(got, want, 2e-3);
  }
  // ... and bit-for-bit the Python's default 96-point trapezoid
  CHECK_CLOSE(g2_ww(g1, 0.1, 10.0), 0.94403515478179245, kRtol);
}

TEST_CASE("the digitized Miller b1 curve reproduces the published figure") {
  // The tables hold the published PER-DEUTERON b1 and the accessor halves it,
  // because every consumer pairs b1 with a per-nucleon F1.  At x = 0.012 the
  // digitized total is 0.114 per deuteron, against HERMES's 0.112 +- 0.055.
  CHECK_CLOSE(2.0 * toy_b1(0.012, 2.5, 0.0), 0.11429317074113018, kRtol);
  CHECK_CLOSE_AT(2.0 * toy_b1(0.012, 2.5, 0.0), 0.114, 0.0, 0.01);
  // constant extrapolation below the table's x = 0.01
  CHECK_CLOSE(toy_b1(0.001, 4.0, 1.0), 0.064679500000000001, kRtol);
  CHECK(toy_b1(1e-5, 4.0, 1.0) == toy_b1(0.001, 4.0, 1.0));
  // ... and the (1-x)^3 taper above x = 0.900, which is what keeps b1/F1 from
  // diverging at the generator's x = 0.955 cell
  CHECK_CLOSE(toy_b1(0.95, 4.0, 1.0), 1.9532280777319135e-05, kRtol);
  CHECK(toy_b1(0.95, 4.0, 1.0) < toy_b1(0.9, 4.0, 1.0));
}

TEST_CASE("the digitized CDKS convolution is the other b1 camp") {
  CHECK_CLOSE(b1_convolution(0.3, 2.5, 0.0), -0.0002816198975086724, kRtol);
  CHECK_CLOSE(b1_convolution(0.05, 2.5, 0.0), 5.8517162680014353e-05, kRtol);
  CHECK(std::fabs(b1_convolution(0.3, 2.5, 0.0)) < 1e-3);
  CHECK(std::fabs(b1_convolution(0.25, 2.5, 0.0)) < 1e-3);
  // two sign changes in 0.01 < x < 0.5 (at 0.06 and 0.42), which the old
  // `0.1 * toy_b1` never had
  const DigitizedTable& t = tables::kB1CdksQ2p5();
  const double* v = t.column("xb1_theory1_sum");
  int flips = 0;
  double prev = 0.0;
  bool have_prev = false;
  for (std::size_t i = 0; i < t.n; ++i) {
    if (t.x[i] <= 0.02 || t.x[i] >= 0.5) continue;
    const double sgn = (v[i] > 0) - (v[i] < 0);
    if (have_prev && sgn != prev) ++flips;
    prev = sgn;
    have_prev = true;
  }
  CHECK(flips == 2);
  // the Close-Kumano integral is respected by the convolution and not by
  // Miller's total (his Sec. V) -- reported, not enforced
  CHECK_CLOSE(close_kumano_integral(true), 0.00045920030669150986, 1e-14);
  CHECK_CLOSE(close_kumano_integral(false), 0.0059147403185172498, 1e-14);
  CHECK(std::fabs(close_kumano_integral(true)) < 1e-3);
  CHECK(std::fabs(close_kumano_integral(false))
        > 5.0 * std::fabs(close_kumano_integral(true)));
}

TEST_CASE("the legacy analytic b1 shapes are unchanged") {
  const double x = 0.05, q2 = 4.0, f1 = 6.0;
  const double shape = 0.01 * std::pow(x, -0.2) * (1.0 - x / 0.20)
                       * std::exp(-3.0 * x) * f1;
  CHECK_CLOSE(toy_b1(x, q2, f1, B1Mode::kToy), shape, kRtol);
  CHECK_CLOSE(b1_convolution(x, q2, f1, B1Mode::kToy), 0.1 * shape, kRtol);
}

TEST_CASE("the 6Li b1 transfer constants") {
  CHECK_CLOSE(LI6_B1_RANK2_TRANSFER, 0.921947, 1e-6);
  CHECK_CLOSE(LI6_B1_PER_NUCLEON, 1.0 / 3.0, 1e-15);
  CHECK_CLOSE(b1_li6_from_deuteron(1.0), 0.30731566666666665, kRtol);
  CHECK_CLOSE(b1_li6_from_deuteron(3.0, LI6_B1_LEGACY_TRANSFER, 1.0), 2.61,
              1e-14);
}

TEST_CASE("the digitized polarized-EMC curves and the valence transfer") {
  const double xs[4] = {0.10, 0.30, 0.50, 0.70};
  const double unpol[4] = {1.0239332568856121, 0.99446484083658282,
                           0.93914054732364405, 0.91028333333333333};
  const double pol23[4] = {0.9270243359682393, 0.93043464976958523,
                           0.91428491634172282, 0.90516809929078013};
  const double tmt[4] = {0.92581573098751413, 0.92510085569517631,
                         0.84892190839075798, 0.88122116416700425};
  const EmcBaseline kLegacy = EmcBaseline::LegacyTable;
  for (int i = 0; i < 4; ++i) {
    CHECK_CLOSE(cbt_unpolarized_emc_ratio(xs[i]), unpol[i], kRtol);
    CHECK_CLOSE(cbt_published_emc_ratio(xs[i]), pol23[i], kRtol);
    // on the LEGACY baseline CBT is referenced to itself, so the transferred
    // curve IS the published one
    CHECK(cbt_polarized_emc_ratio(xs[i], EmcMode::kDigitized, 23, kLegacy)
          == pol23[i]);
    CHECK_CLOSE(tmt_published_emc_ratio(xs[i]), tmt[i], kRtol);
  }
  // CBT computed 7Li itself, so on its own baseline the valence scale is
  // exactly 1; TMT's nuclear matter is scaled down to 7Li strength by 0.397
  CHECK(cbt_valence_scale(kLegacy) == 1.0);
  CHECK_CLOSE(tmt_valence_scale(kLegacy), 0.39700861081338656, kRtol);
  CHECK_CLOSE(tmt_polarized_emc_ratio(0.45, EmcMode::kDigitized, kLegacy),
              0.95204010409497919, kRtol);
  CHECK_CLOSE(tmt_polarized_emc_ratio(0.45, EmcMode::kDigitized, kLegacy),
              1.0 - tmt_valence_scale(kLegacy)
                        * (1.0 - tmt_published_emc_ratio(0.45)),
              kRtol);
  // "about twice" vs "about equal" over the valence region
  const double rx[4] = {0.40, 0.45, 0.50, 0.60};
  const double cbt_r[4] = {2.25, 1.69, 1.41, 1.14};
  const double tmt_r[4] = {1.01, 0.98, 1.00, 1.08};
  for (int i = 0; i < 4; ++i) {
    CHECK_CLOSE_AT(cbt_ratio_of_effects(rx[i]), cbt_r[i], 0.0, 0.05);
    CHECK_CLOSE_AT(tmt_ratio_of_effects(rx[i]), tmt_r[i], 0.0, 0.05);
  }
  // the legacy constant modes reproduce the pre-digitization curves
  CHECK_CLOSE(cbt_polarized_emc_ratio(0.3, EmcMode::kConstant),
              1.0 - 2.0 * (1.0 - unpolarized_emc_ratio(0.3)), kRtol);
  CHECK_CLOSE(tmt_polarized_emc_ratio(0.3, EmcMode::kConstant),
              unpolarized_emc_ratio(0.3), kRtol);
  CHECK_CLOSE(unpolarized_emc_ratio(0.5), 0.93666666666666665, kRtol);
  CHECK_THROWS(cbt_polarized_emc_ratio(0.3, EmcMode::kDigitized, 24));
}

// C5.  The polarized-EMC transfer is referenced to an explicit unpolarized
// BASELINE, and the library default is the Python's default (epps21) -- not
// the pre-2026-08-29 CBT-on-CBT one that was hard-coded as `1.0` and
// `0.397009`.  Every number below is
// `polli_fastsim.polarized.valence_scale / *_polarized_emc_ratio` run at that
// baseline on 2026-08-29.
TEST_CASE("the polarized-EMC valence transfer names its unpolarized baseline") {
  // --- the two baselines' own valence depletion <1 - R_unpol>
  CHECK_CLOSE(emc_valence_depletion(EmcBaseline::LegacyTable),
              0.05834952032138685, kRtol);
  CHECK(emc_valence_depletion(EmcBaseline::Epps21)
        == 0.031052077003862335);
  // EPPS21 is just over HALF as deep as the digitized CBT curve, which is the
  // whole content of the change: both transferred curves shrink by that factor
  CHECK_CLOSE_AT(emc_valence_depletion(EmcBaseline::Epps21)
                     / emc_valence_depletion(EmcBaseline::LegacyTable),
                 0.5322, 0.0, 5e-4);

  // --- both scales, from the SAME code path, against the Python
  CHECK(cbt_valence_scale(EmcBaseline::LegacyTable) == 1.0);
  CHECK_CLOSE(tmt_valence_scale(EmcBaseline::LegacyTable),
              0.39700861081338656, kRtol);
  CHECK_CLOSE(cbt_valence_scale(EmcBaseline::Epps21),
              0.5321736465497698, kRtol);
  CHECK_CLOSE(tmt_valence_scale(EmcBaseline::Epps21),
              0.2112775201282183, kRtol);

  // --- the DEFAULT is the Python's default
  CHECK(EMC_BASELINE_DEFAULT == EmcBaseline::Epps21);
  CHECK(cbt_valence_scale() == cbt_valence_scale(EmcBaseline::Epps21));
  CHECK(tmt_valence_scale() == tmt_valence_scale(EmcBaseline::Epps21));

  // --- the transferred curves themselves, both baselines, both camps
  struct Row { double x, cbt23_legacy, cbt23_epps, tmt_legacy, tmt_epps; };
  const Row rows[5] = {
      {0.10, 0.9270243359682393,  0.9611642747628271,
             0.9705482064151464,  0.9843265316105173},
      {0.30, 0.9304346497695852,  0.9629791538943683,
             0.9702643947684306,  0.9841754945315513},
      {0.45, 0.9183139996455157,  0.9565288633192883,
             0.9520401040949792,  0.9744770073080777},
      {0.50, 0.9142849163417228,  0.954384691365256,
             0.9400206967258773,  0.9680805954590955},
      {0.70, 0.9051680992907801,  0.9495329615903287,
             0.952843779391911,   0.9749047021214879}};
  for (const Row& r : rows) {
    CHECK_CLOSE(cbt_polarized_emc_ratio(r.x, EmcMode::kDigitized, 23,
                                        EmcBaseline::LegacyTable),
                r.cbt23_legacy, kRtol);
    CHECK_CLOSE(cbt_polarized_emc_ratio(r.x, EmcMode::kDigitized, 23,
                                        EmcBaseline::Epps21),
                r.cbt23_epps, kRtol);
    CHECK_CLOSE(tmt_polarized_emc_ratio(r.x, EmcMode::kDigitized,
                                        EmcBaseline::LegacyTable),
                r.tmt_legacy, kRtol);
    CHECK_CLOSE(tmt_polarized_emc_ratio(r.x, EmcMode::kDigitized,
                                        EmcBaseline::Epps21),
                r.tmt_epps, kRtol);
    // the default really is the epps21 column
    CHECK(cbt_polarized_emc_ratio(r.x)
          == cbt_polarized_emc_ratio(r.x, EmcMode::kDigitized, 23,
                                     EmcBaseline::Epps21));
    CHECK(tmt_polarized_emc_ratio(r.x)
          == tmt_polarized_emc_ratio(r.x, EmcMode::kDigitized,
                                     EmcBaseline::Epps21));
  }

  // --- the ratio-of-effects statements are on each camp's PUBLISHED figure,
  // so they do not move with the baseline at all
  for (double x : {0.40, 0.45, 0.50, 0.60}) {
    const double c = cbt_ratio_of_effects(x);
    const double t = tmt_ratio_of_effects(x);
    CHECK(c == (1.0 - cbt_published_emc_ratio(x))
                   / (1.0 - cbt_unpolarized_emc_ratio(x)));
    CHECK(std::isfinite(t));
  }
  CHECK_CLOSE(cbt_ratio_of_effects(0.40), 2.2487868592171982, kRtol);
  CHECK_CLOSE(cbt_ratio_of_effects(0.60), 1.1366929792889406, kRtol);
  CHECK_CLOSE(tmt_ratio_of_effects(0.45), 0.98408316213139446, kRtol);
}

TEST_CASE("toy_delta_gluon is the scenario cos-2phi source") {
  CHECK_CLOSE(toy_delta_gluon(0.1, 4.0, 2.0), 0.00065765788796570657, kRtol);
  CHECK_CLOSE(toy_delta_gluon(0.1, 4.0, 2.0, 3e-3),
              3.0 * toy_delta_gluon(0.1, 4.0, 2.0), kRtol);
  CHECK(toy_delta_gluon(1.0, 4.0, 2.0) == 0.0);  // vanishes at the elastic edge
}

TEST_CASE("the Delta moment ansatz A and B") {
  // shapes are peak-normalized, so A * alpha_s is the peak Delta/F1
  CHECK_CLOSE(shape_normalized(0.3, "mid_x"), 0.88857273914267154, kRtol);
  CHECK_CLOSE(shape_normalized(0.3, "low_x"), 0.4966684278915946, kRtol);
  const DeltaVariant mid = delta_variant("mid_x");
  CHECK_CLOSE(shape_normalized(mid.alpha / (mid.alpha + mid.beta), "mid_x"), 1.0,
              1e-14);
  CHECK_THROWS(shape_normalized(0.3, "nonsense"));

  // Interpretation B is the analytic Beta-function integral
  CHECK_CLOSE(solve_A_interp_b("low_x"), -0.18355179570297245, kRtol);
  CHECK_CLOSE(solve_A_interp_b("mid_x"), -0.088951691418711606, kRtol);
  CHECK_CLOSE(solve_A_interp_b("high_x"), -0.047616034334759805, kRtol);

  CHECK_CLOSE(alpha_s_lo(4.47), 0.33320440488583836, kRtol);

  // Interpretation A on the per-nucleon toy F1 of 6Li at <Q2> = 4.47
  const auto base = std::make_shared<const ToyF2>();
  const NuclearF2 nf(LI6(), base);
  const auto f1_per_nucleon = [&nf](double x, double q2) {
    return nf.f1a(x, q2) / 6.0;
  };
  CHECK_CLOSE(solve_A_interp_a(f1_per_nucleon, 4.47), -0.29241049840786215,
              1e-14);

  // and the models themselves close the sum rule they were solved from
  const DeltaModel ma = make_moment_a(f1_per_nucleon, 4.47, "mid_x");
  const DeltaModel mb = make_moment_b("mid_x");
  const double q2 = 4.47;
  CHECK_CLOSE(ma(0.3, q2, 2.0),
              solve_A_interp_a(f1_per_nucleon, 4.47) * alpha_s_lo(q2) * 2.0
                  * shape_normalized(0.3, "mid_x"),
              1e-14);
  CHECK_CLOSE(mb(0.3, q2, 2.0),
              solve_A_interp_b("mid_x") * alpha_s_lo(q2)
                  * shape_normalized(0.3, "mid_x"),
              kRtol);
  CHECK(mb(0.3, q2, 2.0) == mb(0.3, q2, 99.0));  // B does not carry F1
  // the dilution factor is the 6Li 2-of-6 convention of plans/04 #6
  const DeltaModel mb3 = make_moment_b("mid_x", C_BAG, 1.0 / 3.0);
  CHECK_CLOSE(mb3(0.3, q2, 1.0), mb(0.3, q2, 1.0) / 3.0, kRtol);
  CHECK(ma.name() == "moment_A");
  CHECK(mb.info().find("moment_B") == 0);
}

TEST_CASE("TensorSF backends adapt to the kernel's callable slots") {
  const MillerB1 miller;
  const CdksB1 cdks;
  CHECK_CLOSE(miller.b1_func()(0.1, 4.0, 1.0), toy_b1(0.1, 4.0, 1.0), kRtol);
  CHECK_CLOSE(miller.b2_func()(0.1, 4.0, 1.0), 2.0 * 0.1 * toy_b1(0.1, 4.0, 1.0),
              kRtol);
  CHECK(miller.delta_func()(0.1, 4.0, 1.0) == 0.0);
  CHECK_CLOSE(cdks.b1_func()(0.1, 4.0, 1.0), b1_convolution(0.1, 4.0, 1.0),
              kRtol);
  const auto d = std::make_shared<const MillerB1>();
  const Li6B1 li6(d);
  CHECK_CLOSE(li6.b1(0.1, 4.0, 1.0),
              b1_li6_from_deuteron(toy_b1(0.1, 4.0, 1.0)), kRtol);
}
