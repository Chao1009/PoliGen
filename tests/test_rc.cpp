// Tensor-sector radiative corrections (rc.hpp), the OPT-IN weight-only family
// of docs/open_items/run_2026-09-02/design_C_tensor_rc.md.
//
// This file grows with the model.  Landed so far:
//
//   T3  the band delta(x): monotone non-increasing in x, both anchors returned
//       EXACTLY (the clamp branches return the constant unmodified, which is
//       what makes `==` legitimate there), and the design's own delta table
//       reproduced at 1e-12 against the full-precision doubles it prints.
//
// Still to come with the model bodies: T1/T2 (the band identity, inclusive and
// tagged), T4 (signs), T5 (the tail's y and t_min behaviour), T7 (the Born
// normalisation against POLRAD Eq. (9)), T8/T8' (the tail normalisation gates),
// T9 (SKIP in v0), T10-T12' (the form factors and their bands), T13 (Off is
// identically 1), T14 (the per-channel table) and T18 (threading).

#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "check_close.hpp"
#include "doctest.h"
#include "lipolgen/beams.hpp"
#include "lipolgen/bookkeeping.hpp"
#include "lipolgen/rc.hpp"
#include "lipolgen/sampler.hpp"
#include "lipolgen/sf.hpp"
#include "lipolgen/xsec.hpp"

using namespace lipolgen;

namespace {

/// A spin-1 kernel with BOTH rank-2 slots alive -- the b-sector (b1, b2) that
/// the RC band prices AND the Delta (cos 2phi) term it deliberately leaves out
/// -- so the `with_delta` split of `tensor_amplitudes` is actually exercised.
std::shared_ptr<const InclusiveKernel> band_kernel(bool tensor_gamma) {
  InclusiveKernel::Options o;
  o.b1_func = [](double x, double q2, double f1) { return toy_b1(x, q2, f1); };
  o.delta_func = [](double x, double q2, double f1) {
    return toy_delta_gluon(x, q2, f1, 1.0);
  };
  o.tensor_gamma = tensor_gamma;
  return std::make_shared<const InclusiveKernel>(LI6(), o);
}

/// A small (x, Q2) grid: this file is checking an algebraic identity cell by
/// cell, not a rate, so the 100x72 production grid buys nothing.
const InclusiveSampler& small_sampler() {
  static const InclusiveSampler s = [] {
    InclusiveSampler::GridSpec g;
    g.nx = 20; g.nq2 = 12;
    return InclusiveSampler(band_kernel(false), default_configs("6Li")[1],
                            Scenario(), g, /*with_perp=*/true);
  }();
  return s;
}

}  // namespace

TEST_CASE("T3: the RC band delta(x) is monotone, clamped and anchored") {
  SUBCASE("both anchors and both clamps are EXACT") {
    // Legitimate `==`: `x >= x_high` and `x <= x_low` return the anchor
    // constant unmodified, with no arithmetic in between.
    CHECK(rc_delta(RC_X_HIGH) == RC_DELTA_HIGH_X);
    CHECK(rc_delta(0.9) == RC_DELTA_HIGH_X);
    CHECK(rc_delta(RC_X_LOW) == RC_DELTA_LOW_X);
    CHECK(rc_delta(1e-5) == RC_DELTA_LOW_X);
    // ... and with non-default anchors too, since the CLI can move them.
    CHECK(rc_delta(0.2, 0.02, 0.25, 0.2, 0.005) == 0.02);
    CHECK(rc_delta(0.005, 0.02, 0.25, 0.2, 0.005) == 0.25);
  }

  SUBCASE("monotone non-increasing on a 2000-point log grid in [1e-5, 0.9]") {
    const int n = 2000;
    const double lo = std::log(1e-5), hi = std::log(0.9);
    double prev = rc_delta(1e-5);
    for (int i = 1; i < n; ++i) {
      const double x = std::exp(lo + (hi - lo) * i / (n - 1));
      const double d = rc_delta(x);
      CHECK_MESSAGE(d <= prev, "delta rose at x = " << x << ": " << prev
                                                    << " -> " << d);
      CHECK(d >= RC_DELTA_HIGH_X);
      CHECK(d <= RC_DELTA_LOW_X);
      prev = d;
    }
  }

  SUBCASE("the design's delta table, at the printed double precision") {
    // design_C_tensor_rc.md section 1.3.  A 5-significant-figure table cannot
    // meet 1e-12, which is why these are the full doubles.
    struct Row { double x, delta; };
    const Row rows[] = {
        {0.001, 0.3},
        {0.005, 0.3},
        {0.01,  0.3},
        {0.02,  0.22874999999999995},
        {0.03,  0.18707142182361763},
        {0.05,  0.13456262323927543},
        {0.063, 0.1108061822113555},
        {0.08,  0.08625000000000000},
        {0.1,   0.06331262323927542},
        {0.16,  0.015},
        {0.3,   0.015},
    };
    for (const Row& r : rows) {
      CHECK_CLOSE(rc_delta(r.x), r.delta, 1e-12);
    }
    // The interpolation must NOT contradict HERMES's measured 15 % residual at
    // x = 0.063 -- which the old x_high = 0.05 version did (it gave 1.5 %).
    CHECK(rc_delta(0.063) > 0.10);
  }

  SUBCASE("the optimistic low-x edge is contained by the conservative one") {
    for (double x : {1e-4, 0.003, 0.01, 0.02, 0.05, 0.1, 0.2}) {
      const double conservative =
          rc_delta(x, RC_DELTA_HIGH_X, RC_DELTA_LOW_X);
      const double optimistic =
          rc_delta(x, RC_DELTA_HIGH_X, RC_DELTA_LOW_X_OPTIMISTIC);
      CHECK(optimistic <= conservative);
    }
  }

  SUBCASE("bad anchors throw instead of returning a silent NaN") {
    CHECK_THROWS_AS(rc_delta(0.05, 0.015, 0.30, 0.16, 0.0), std::runtime_error);
    CHECK_THROWS_AS(rc_delta(0.05, 0.015, 0.30, 0.01, 0.01),
                    std::runtime_error);
    CHECK_THROWS_AS(rc_delta(0.05, 0.015, 0.30, 0.005, 0.16),
                    std::runtime_error);
  }
}

TEST_CASE("rc.hpp names: the weight block and the mode strings") {
  CHECK(std::string(rc_mode_name(RcMode::Off)) == "off");
  CHECK(std::string(rc_mode_name(RcMode::TensorBand)) == "tensor-band");
  CHECK(std::string(pipeline_rc_name(PipelineRc::TensorBand)) == "tensor-band");

  REQUIRE(kRcWeightCount == 3);
  CHECK(std::string(rc_weight_name(0)) == "rc_tensor_lo");
  CHECK(std::string(rc_weight_name(1)) == "rc_tensor_hi");
  CHECK(std::string(rc_weight_name(2)) == "rc_tail");
  CHECK_THROWS_AS(rc_weight_name(3), std::runtime_error);

  // Slot 0 is the event's own state and keeps the bare name, so an unweighted
  // run's HepMC3 names are exactly {nominal, rc_tensor_lo, rc_tensor_hi,
  // rc_tail}; slot 1 + k carries the "_k" suffix of "spin_weight_k".
  CHECK(rc_weight_name(2, 0) == "rc_tail");
  CHECK(rc_weight_name(0, 1) == "rc_tensor_lo_1");
  CHECK(rc_weight_name(2, 3) == "rc_tail_3");
}

// ---------------------------------------------------------------------------
// The `tensor_amplitudes` extraction (xsec.hpp) and the tensor StateTables
// (sampler.hpp) the RC band reads.  These are the guards on "the refactor
// moves no number": `amplitudes()` is DEFINED as the tensor part plus the
// vector terms, so the two cannot drift apart, and the sampler's tensor
// vectors are the rank-2 projection of the very density it drew from.

TEST_CASE("xsec: amplitudes() == tensor_amplitudes() + the vector terms") {
  const double s = 3980.0;
  const double xs[] = {1e-3, 0.01, 0.1, 0.4};
  const double q2s[] = {1.5, 8.0, 60.0};
  const double ths[] = {0.0, 0.7, kPi / 2.0, 2.2};

  for (int tg = 0; tg < 2; ++tg) {
    const std::shared_ptr<const InclusiveKernel> kp = band_kernel(tg == 1);
    const InclusiveKernel& k = *kp;
    for (double x : xs) {
      for (double q2 : q2s) {
        const double y = q2 / (s * x);
        if (!(y > 0.0 && y < 0.99)) continue;
        const SFTables t = k.tables(x, q2, true);
        for (double m : {1.0, 0.0, -1.0}) {
          for (double th : ths) {
            const EventSpinState unpol{0, 0.0, 1.0, m, th, 0.3};

            // (1) At an UNPOLARISED beam -- the configuration the whole A_zz
            //     programme, and every tau identity of rc.hpp, runs at -- the
            //     density IS its tensor part.  BIT for bit, not to a tolerance.
            const Amplitudes a = k.amplitudes(t, x, q2, s, unpol, true);
            const Amplitudes ta =
                k.tensor_amplitudes(t, x, q2, s, unpol, /*with_delta=*/true);
            CHECK(a.w_avg == ta.w_avg);
            CHECK(a.a1 == ta.a1);
            CHECK(a.a2 == ta.a2);

            // (2) `with_delta = false` (RcScope::TensorRate, the default)
            //     removes the gluon-transversity cos 2phi term and NOTHING
            //     else: Delta is not in POLRAD's b1..b4 basis.
            const Amplitudes tr =
                k.tensor_amplitudes(t, x, q2, s, unpol, /*with_delta=*/false);
            CHECK(tr.w_avg == ta.w_avg);
            CHECK(tr.a1 == ta.a1);
            const double delta_term =
                -(1.0 - y) / (y * y) * k.tensor_moments(m).second *
                std::sin(th) * std::sin(th) * t.delta /
                std::max(InclusiveKernel::dphi(t, x, y), 1e-30);
            CHECK_CLOSE(ta.a2 - tr.a2, delta_term, 1e-12);
            // ... and it really is a live term wherever the axis is off-beam.
            if (std::fabs(std::sin(th)) > 1e-9 && m != 0.0) {
              CHECK(tr.a2 != ta.a2);
            }

            // (3) With a POLARISED beam the difference is exactly the vector
            //     sector -- which is what tau must NOT contain.
            const EventSpinState pol{+1, 0.7, 1.0, m, th, 0.3};
            const Amplitudes ap = k.amplitudes(t, x, q2, s, pol, true);
            const Amplitudes tp =
                k.tensor_amplitudes(t, x, q2, s, pol, /*with_delta=*/true);
            CHECK(tp.a2 == ap.a2);
            if (m != 0.0) {
              // atol is scaled to w_avg itself: at theta_S = pi/2 the vector
              // term is O(1e-20) (cos(pi/2) is 6e-17 in double, not 0) and
              // adding it to an O(1e-2) w_avg is a no-op, so a pure RELATIVE
              // comparison there would be comparing two roundoffs.
              CHECK_CLOSE_AT(ap.w_avg - tp.w_avg,
                             0.7 * m * std::cos(th) * k.a_parallel(t, x, q2, y),
                             1e-12, 1e-15 * std::fabs(ap.w_avg));
            } else {
              CHECK(ap.w_avg == tp.w_avg);
            }
          }
        }
      }
    }
  }
}

TEST_CASE("xsec: tensor_amplitudes refuses the same spin mismatch amplitudes does") {
  const std::shared_ptr<const InclusiveKernel> kp = band_kernel(false);
  const SFTables t = kp->tables(0.2, 5.0, true);
  const double s = 3980.0;
  CHECK_THROWS_AS(
      kp->tensor_amplitudes(t, 0.2, 5.0, s, EventSpinState{0, 0.0, 1.5, 1.5}),
      std::runtime_error);
  CHECK_NOTHROW(
      kp->tensor_amplitudes(t, 0.2, 5.0, s, EventSpinState{0, 0.0, 1.0, 1.0}));
}

TEST_CASE("sampler: the tensor StateTables are tensor_amplitudes, cell by cell") {
  const InclusiveSampler& s = small_sampler();
  const SpinCategory cat("t+", 1.0, {1.0, 0.0, 0.0}, 0, 0.0, 0.6, 0.2);
  const InclusiveSampler::StateTables& st = s.state_tables(cat, 1.0);
  REQUIRE(st.w_tensor.size() == st.w_avg.size());
  REQUIRE(st.w_tensor.size() == s.n_cells());

  EventSpinState state;
  state.lam_e = cat.lam_e; state.pe = cat.pe; state.j = cat.j; state.m = 1.0;
  state.theta_s = cat.theta_s; state.phi_s = cat.phi_s;
  bool any_nonzero = false;
  for (std::size_t i = 0; i < s.n_cells(); ++i) {
    const Amplitudes ta = s.kernel().tensor_amplitudes(
        s.tables()[i], s.x_cells()[i], s.q2_cells()[i], s.s(), state);
    CHECK(st.w_tensor[i] == ta.w_avg);
    CHECK(st.a1_tensor[i] == ta.a1);
    CHECK(st.a2_tensor[i] == ta.a2);
    // At an unpolarised beam the RATE shift is PURELY rank 2, so `w_tensor`
    // must be the whole of `w_avg` -- the one place the two coincide, and the
    // check that the extraction did not lose a term.
    CHECK(st.w_tensor[i] == st.w_avg[i]);
    if (st.w_tensor[i] != 0.0) any_nonzero = true;
  }
  CHECK(any_nonzero);
}

TEST_CASE("sampler: tensor_weights_for is the rank-2 part of weights_for") {
  const InclusiveSampler& s = small_sampler();
  const SpinCategory unpol("unpol", 1.0, {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0});
  // A tensor-thirds-like pair: one m0-enriched mixture and one pure state, so
  // both the MIXTURE branch and the pure one are exercised.
  const std::vector<SpinCategory> cats = {
      SpinCategory("t0", 1.0, {0.15, 0.7, 0.15}),
      SpinCategory("m+", 1.0, {1.0, 0.0, 0.0}, 0, 0.0, 0.5, 0.25)};
  const EventBatch ev = s.sample_n(unpol, 500, 11, 0, 0);
  const std::vector<double> w = s.weights_for(ev, cats);
  const std::vector<double> tw = s.tensor_weights_for(ev, cats);
  REQUIRE(tw.size() == w.size());

  const std::vector<double> ms = m_values(1.0);
  bool any_nonzero = false;
  for (std::size_t k = 0; k < cats.size(); ++k) {
    for (std::size_t i = 0; i < ev.size(); ++i) {
      const std::size_t c = static_cast<std::size_t>(ev.cell[i]);
      const double phip = ev.phi[i] - cats[k].phi_s;
      double want = 0.0;
      for (std::size_t im = 0; im < ms.size(); ++im) {
        const double p_m = cats[k].populations[im];
        if (p_m <= 0.0) continue;
        const InclusiveSampler::StateTables& st =
            s.state_tables(cats[k], ms[im]);
        want += p_m * (st.w_tensor[c] + st.a1_tensor[c] * std::cos(phip) +
                       st.a2_tensor[c] * std::cos(2.0 * phip));
      }
      CHECK(tw[i * cats.size() + k] == want);
      // tau = t/w is the rank-2 fraction of the density: bounded, and it is
      // what the band multiplies delta(x) by.
      const double tau = tw[i * cats.size() + k] / w[i * cats.size() + k];
      CHECK(std::fabs(tau) < 1.0);
      if (tau != 0.0) any_nonzero = true;
    }
  }
  CHECK(any_nonzero);
}
