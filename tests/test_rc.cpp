// SPDX-License-Identifier: GPL-3.0-or-later
// Tensor-sector radiative corrections (rc.hpp), the OPT-IN weight-only family
// of docs/open_items/run_2026-09-02/design_C_tensor_rc.md.
//
// Design sec. 5's numbering is kept, so a reader can go row by row:
//
//   T1   the band identity against the PUBLISHED A_zz          T8'  the QRT, Eq. (44)
//   T2   the tagged identity (the SAME closed form)            T9   Eq. (A.4)'s Q_N split
//   T3   the band delta(x)                                     T10  the deuteron FF normalisation
//   T4   the signs, incl. the TENSOR_LL_SIGN-free form         T11  the 6Li FF anchors
//   T5   the tail's y and t_min behaviour                      T12  fq_scale is QUADRATIC
//   T7   the Born IS dsigma_unpol                              T12' tail_tensor_scale moves it
//   T8(0) the tail ADDS events (the sign gate)                 T13  Off is identically 1
//   T8(a) the deuteron gate, Eq. (38) written HERE             T14  the per-Channel table
//   T8(a') the per-nucleon reduction is m_p/M_A                T17  weighted mode
//   T8(b) the carbon Z^2                                       T18  threading
//   T8(c) the A = 1 spin-1/2 limit  (+ T8(c'), POLRAD's own numbers)
//
// T6, T15, T16 and the pipeline halves of T13/T17 are the wiring agent's
// (design sec. 6, agent 2) and live in test_pipeline / test_hepmc / pytest.

#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <thread>
#include <vector>

#include "check_close.hpp"
#include "doctest.h"
#include "lipolgen/beams.hpp"
#include "lipolgen/bookkeeping.hpp"
#include "lipolgen/cluster_config.hpp"
#include "lipolgen/rc.hpp"
#include "lipolgen/sampler.hpp"
#include "lipolgen/asymmetries.hpp"
#include "lipolgen/sf.hpp"
#include "lipolgen/spin.hpp"
#include "lipolgen/tagged.hpp"
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

// ===========================================================================
// The tensor RC model itself (design_C_tensor_rc.md secs. 1.4, 1.5, 2.1;
// POLRAD 2.0 Eqs. (37)-(40), (43), (44); every transcription verdict of
// docs/open_items/run_2026-09-02/polrad_transcription_check.md).
//
// WHERE THE CHECK OVERRULES THE DESIGN, the check wins and the test says so.

namespace {

/// A 6Li sampler for the RC model.  Coarser than production (the tail tables
/// are 3 one-dimensional quadratures per node) but on the GENERATOR scenario,
/// so the y window reaches 0.985 and T5's y = 0.97 point exists.
const InclusiveSampler& rc_sampler() {
  static const InclusiveSampler s = [] {
    InclusiveSampler::GridSpec g;
    g.nx = 40; g.nq2 = 24;
    return InclusiveSampler(band_kernel(false), default_configs("6Li")[1],
                            generator_scenario(Scenario()), g,
                            /*with_perp=*/true);
  }();
  return s;
}

std::shared_ptr<const InclusiveSampler> rc_sampler_ptr() {
  static const std::shared_ptr<const InclusiveSampler> p(
      &rc_sampler(), [](const InclusiveSampler*) {});   // non-owning
  return p;
}

/// The default 6Li inclusive model at the run's defaults.
const RcModel& rc_model() {
  static const RcModel m(RcMode::TensorBand, RcOptions(), rc_sampler_ptr(),
                         tensor_thirds_plan(0.0, 0.6), Channel::Inclusive,
                         LI6());
  return m;
}

/// One inclusive Event built BY HAND, reading only the fields
/// `RcModel::tensor_fraction` documents.
Event inclusive_event(const InclusiveSampler& s, int cell, double m,
                      double phi = 0.0, double theta_s = 0.0) {
  Event ev;
  ev.channel = Channel::Inclusive;
  ev.kin.cell = cell;
  ev.kin.x = s.x_cells()[static_cast<std::size_t>(cell)];
  ev.kin.q2 = s.q2_cells()[static_cast<std::size_t>(cell)];
  ev.kin.s = s.s();
  ev.kin.y = ev.kin.q2 / (ev.kin.x * ev.kin.s);
  ev.kin.phi = phi;
  ev.spin.j = 1.0;
  ev.spin.m_ion = m;
  ev.spin.lam_e = 0;
  ev.spin.pe = 0.0;
  ev.spin.theta_s = theta_s;
  ev.spin.phi_s = 0.0;
  return ev;
}

/// The deuteron form factor of T8/T10: the SAME HO shape, normalised on the
/// deuteron's MEASURED moments.  `HoSpin1FF::for_ion` fills M_A and Z from the
/// Ion and throws on any field left at its 0 sentinel, so a deuteron object
/// exists only when the caller supplies the moments explicitly -- which is
/// exactly the "no ion-specific default" rule (T11 asserts the throw).
std::shared_ptr<HoSpin1FF> deuteron_ff(double fq_scale = 1.0,
                                       double fm_scale = 1.0) {
  HoSpin1FFOptions o;
  o.mu_n = 0.857406;      ///< POLRAD's own `data dmu/0.857406d0/` (adgh:5489)
  o.q_fm2 = 0.2859;       ///< design sec. 2.1's Q_d
  // The SHAPE is a stand-in for POLRAD's `ffdeu` -- every T8(a)/T10 gate is
  // shape-blind -- but T8(c') compares MAGNITUDES against POLRAD's own
  // compiled numbers, so it is set to reproduce the deuteron's own
  // <r^2>_point = r_d^2 - <r^2>_p = 2.1413^2 - 0.7071 = 3.88 fm^2 rather than
  // left arbitrary.  alpha is near zero: the deuteron's C0 has no zero in the
  // range that matters, unlike 6Li's.
  o.a_fm = 1.593;
  o.alpha = 0.02;
  o.fm_qz_fm = 1.3;
  o.fm_b_fm = 1.85;
  o.fq_scale = fq_scale;
  o.tail_tensor_scale = fm_scale;
  return HoSpin1FF::for_ion(DEUTERON(), o);
}

/// A GAUSSIAN spin-1 form factor whose WIDTH the caller chooses, so that a
/// test can put the s-/p-peak vertex outside it and the t-peak inside.  That
/// is the switch T19(a) needs: Eq. (38) is Eq. (18)'s ULTRARELATIVISTIC
/// t-PEAK reduction, so the two agree only where the other two peaks are
/// absent, and a form factor 20+ decades down at Q'^2 ~ z_s Q^2 kills them
/// without touching either quadrature.  NOT a physical deuteron -- the
/// magnitudes here are meaningless and only the RATIO of the two models is
/// read.
class GaussianSpin1FF : public Spin1ElasticFF {
 public:
  GaussianSpin1FF(double b_gev2, double mu, double q, double fc0)
      : b_(b_gev2), mu_(mu), q_(q), fc0_(fc0) {}
  double fc(double t) const override { return fc0_ * std::exp(-t / b_); }
  double fm(double t) const override { return mu_ * std::exp(-t / b_); }
  double fq(double t) const override { return q_ * std::exp(-t / b_); }
  std::string provenance() const override { return "test Gaussian"; }

 private:
  double b_, mu_, q_, fc0_;
};

/// A SPIN-0 nucleus: F_m = F_q = 0 and F_c = Z F(t), which is how POLRAD
/// Eq. (38)'s carbon line sits inside the same master formula (T8(b)).
class Spin0FF : public Spin1ElasticFF {
 public:
  Spin0FF(double z, double b_fm) : z_(z), b_(b_fm) {}
  double shape(double t) const {
    const double q2 = std::max(t, 0.0) / (HBARC_GEV_FM * HBARC_GEV_FM);
    return std::exp(-0.25 * q2 * b_ * b_);
  }
  double fc(double t) const override { return z_ * shape(t); }
  double fm(double) const override { return 0.0; }
  double fq(double) const override { return 0.0; }
  std::string provenance() const override { return "test spin-0 gaussian"; }

 private:
  double z_, b_;
};

}  // namespace

// ---------------------------------------------------------------- T10, T11

TEST_CASE("T10: the DEUTERON fixes the (F_c, F_m, F_q) normalisation") {
  // POLRAD Eq. (A.4) at Q_N = 0 is the spin-1 Rosenbluth pair
  //   Im_1|_0 = (2/3) eta (1+eta) F_m^2 == B/2,
  //   Im_2|_0 = F_c^2 + (2/3) eta F_m^2 + (8/9) eta^2 F_q^2 == A,
  // so F_c == G_C, F_m == G_M, F_q == G_Q in the textbook convention and the
  // three t -> 0 values are FIXED, not chosen.  POLRAD's own `ffdeu`, compiled
  // standalone, gives G_C(0) = 1, G_M(0) = 1.7139610634, G_Q(0) = 25.84
  // (polrad_transcription_check.md sec. 5b).
  const std::shared_ptr<HoSpin1FF> ff = deuteron_ff();
  CHECK_CLOSE(ff->fc(0.0), 1.0, 1e-12);          // = Z, and the deuteron has Z = 1
  CHECK_CLOSE(ff->fm(0.0), 1.714, 1e-3);         // textbook G_M(0)
  CHECK_CLOSE(ff->fq(0.0), 25.83, 1e-3);         // textbook G_Q(0)
  // ... and against POLRAD's own number, which is tighter on G_M.  The 2e-8
  // residual is POLRAD's own deuteron mass against this repository's AME one
  // (`beams.cpp`: M_d = 1.875612942 GeV), not a convention difference.
  CHECK_CLOSE(ff->fm(0.0), 1.7139610634, 1e-7);
  // The DERIVED expressions, so that a change to beams.cpp's mass table or to
  // HBARC_GEV_FM moves the test and the code together.
  const double m_d = DEUTERON().mass();
  CHECK_CLOSE(ff->fm(0.0), (m_d / PROTON_MASS) * 0.857406, 1e-14);
  CHECK_CLOSE(ff->fq(0.0), std::pow(m_d / HBARC_GEV_FM, 2) * 0.2859, 1e-14);
}

TEST_CASE("T11: the 6Li form-factor anchors, and for_ion's refusals") {
  const std::shared_ptr<HoSpin1FF> ff = HoSpin1FF::for_ion(LI6());
  const double m_a = LI6().mass();

  SUBCASE("F_c(0) = Z, F_m(0), F_q(0) -- DERIVED, never retyped") {
    // ON BOTH EDGES OF THE C0 BAND.  `C0Shape` moves the SHAPE and never the
    // q -> 0 normalisation -- that is exactly what makes it a band on the
    // unfitted part alone, and the VmcFt edge would otherwise be smuggling
    // WS98's Q(6Li) = -0.23(9) in through the back door.
    for (C0Shape sh : {C0Shape::Ho, C0Shape::VmcFt}) {
      HoSpin1FFOptions o;
      o.c0_shape = sh;
      const std::shared_ptr<HoSpin1FF> f = HoSpin1FF::for_ion(LI6(), o);
      CAPTURE(c0_shape_name(sh));
      CHECK_CLOSE(f->fc(0.0), static_cast<double>(LI6().Z), 1e-14);
      CHECK_CLOSE(f->fm(0.0), (m_a / PROTON_MASS) * LI6_MU_N, 1e-14);
      CHECK_CLOSE(f->fq(0.0),
                  std::pow(m_a / HBARC_GEV_FM, 2) * LI6_QUADRUPOLE_FM2, 1e-14);
      // The design's printed values (sec. 2.1), at its own precision.
      CHECK_CLOSE(f->fc(0.0), 3.0, 1e-14);
      CHECK_CLOSE(f->fm(0.0), 4.90765, 1e-4);
      CHECK_CLOSE(f->fq(0.0), -65.914, 1e-4);
      // TUNL's -0.0818(17) fm^2, NOT Pyykko's -0.0806 (which gives -64.947).
      CHECK(f->fq(0.0) < -65.0);
    }
  }

  SUBCASE("<r^2>_point = 6.078 fm^2 from the POINT shape") {
    // -6 dF/dq^2|_0 / F(0) on the UNFOLDED shape: with `fold_nucleon` the
    // slope would also carry <r^2>_p and <r^2>_n, which is what <r^2>_point
    // was derived by SUBTRACTING (design sec. 2.1).  ON BOTH EDGES: the
    // VmcFt shape is r-rescaled precisely so that this number does not move,
    // and the rescale is computed from the transform's OWN <r^2> (5.9713
    // fm^2) rather than from the file's printed rms (2.4433 fm), because the
    // two differ by 1.4e-4 relative and this gate is 1e-4.
    for (C0Shape sh : {C0Shape::Ho, C0Shape::VmcFt}) {
      HoSpin1FFOptions o;
      o.fold_nucleon = false;
      o.c0_shape = sh;
      CAPTURE(c0_shape_name(sh));
      const std::shared_ptr<HoSpin1FF> pt = HoSpin1FF::for_ion(LI6(), o);
      const double h = 1e-6;                      // in q^2 [fm^-2]
      const double t = h * HBARC_GEV_FM * HBARC_GEV_FM;
      const double slope = (pt->fc(t) - pt->fc(0.0)) / h;
      // 6.0788 is the design's own full value; its table rounds it to 6.078.
      CHECK_CLOSE(-6.0 * slope / pt->fc(0.0), 6.0788, 1e-4);
    }
    // ... and the closed form the design derives it from.
    const double a = LI6_FF_HO_A_FM, al = LI6_FF_HO_ALPHA;
    CHECK_CLOSE(1.5 * a * a * (2.0 + 5.0 * al) / (2.0 + 3.0 * al), 6.0788, 1e-4);
  }

  SUBCASE("the first C0 zero lies in [2.9, 3.3] fm^-1") {
    // NOT pinned to 3.10 +- 0.02: WS98 does not locate the zero (it says only
    // that C2 dominates above 3 fm^-1 and that two-body currents shift the
    // minimum lower, with no number), so the design asks for a REFIT and this
    // test pins only the range the refit must land in.
    //
    // WHERE THE WINDOW COMES FROM, since it is now the repository's ONE
    // statement about this feature (B1.9, 2026-09-04).  [2.9, 3.3] is
    // design_C_tensor_rc.md sec. 2.1's own bracket around its own admitted
    // STARTING GUESS q_0 = 3.1; it is not a measurement, and no elastic data
    // are in this repository to make it one.  `docs/PHYSICS_CHANNELS.md`
    // used to state, uncited, that "the real 6Li form factor dips near
    // |t| ~ 0.31 GeV^2" -- that is q = 2.82 fm^-1, 9.6 % below the shipped
    // q_0 and OUTSIDE this window, i.e. the two sites contradicted each
    // other.  The 0.31 was withdrawn rather than the window widened: it had
    // no source, and the one piece of ab-initio evidence this repository
    // owns disfavours it MORE than it disfavours 3.1 -- the j0 transform of
    // `data/vmc/density/li6.density` is 2.454e-3 +- 1.76e-4 at 2.82 fm^-1
    // (14 sigma from zero) against 1.076e-3 +- 1.60e-4 at 3.0999 (6.7
    // sigma), and has no zero below q ~ 4.3 fm^-1 at all.  Which is why the
    // VmcFt edge below has NO zero, and why a dip location is a MODEL
    // statement on either edge.
    const double a = LI6_FF_HO_A_FM, al = LI6_FF_HO_ALPHA;
    const double q0 = std::sqrt(2.0 * (2.0 + 3.0 * al) / al) / a;
    CHECK(q0 > 2.9);
    CHECK(q0 < 3.3);
    const double t0 = std::pow(q0 * HBARC_GEV_FM, 2);
    CHECK(std::fabs(ff->fc(t0)) < 1e-10 * ff->fc(0.0));
    // The shell-model alpha = (Z-2)/3 = 1/3 would put it at 2.33 fm^-1 -- far
    // too low, which is why alpha is FREE and this is a phenomenological fit.
    CHECK(std::sqrt(2.0 * (2.0 + 1.0) / (1.0 / 3.0)) / a < 2.5);
    // |t| of that zero, the number PHYSICS_CHANNELS.md's coherent row quotes.
    // The whole gated window sits above COHERENT_T_MAX_DEFAULT = 0.2 GeV^2,
    // which is the only thing that row's argument needs.
    CHECK_CLOSE(std::pow(q0 * HBARC_GEV_FM, 2), 0.3742, 2e-4);
    CHECK(std::pow(2.9 * HBARC_GEV_FM, 2) > 0.2);
  }

  SUBCASE("C0Shape::VmcFt -- the other band edge") {
    // The BAND, not a refit: a two-parameter HO cannot reach these values
    // while holding <r^2> (rc.hpp's LI6_VMC_C0_Q_CUT_FM block).
    HoSpin1FFOptions o;
    o.fold_nucleon = false;
    o.c0_shape = C0Shape::VmcFt;
    const std::shared_ptr<HoSpin1FF> v = HoSpin1FF::for_ion(LI6(), o);
    HoSpin1FFOptions oh;
    oh.fold_nucleon = false;
    const std::shared_ptr<HoSpin1FF> h = HoSpin1FF::for_ion(LI6(), oh);
    auto t_of = [](double q) { return std::pow(q * HBARC_GEV_FM, 2); };
    // The two edges where the tail lives: x2.51 in F, x6.3 in F^2 at
    // q = 2 fm^-1; x6.15 and x37.8 at 2.5.
    CHECK_CLOSE(v->fc(t_of(2.0)) / h->fc(t_of(2.0)), 2.508, 1e-3);
    CHECK_CLOSE(v->fc(t_of(2.5)) / h->fc(t_of(2.5)), 6.150, 1e-3);
    // ... and agreeing to 0.6 % where the data are good and both are fitted
    // to the same <r^2>.
    CHECK_CLOSE(v->fc(t_of(0.5)) / h->fc(t_of(0.5)), 1.0, 7e-3);
    // NO C0 ZERO ANYWHERE: positive and strictly falling out to the top of
    // the t-peak's reach, where the Ho edge has already changed sign.  This
    // difference -- not the size of the shape -- is what flips sigma^el_T.
    CHECK(h->fc(t_of(3.5)) < 0.0);
    double prev = v->fc(0.0);
    for (double q = 0.1; q < 90.0; q += 0.1) {
      const double f = v->fc(t_of(q));
      CHECK(f > 0.0);
      CHECK(f < prev);
      prev = f;
    }
    // The continuation is continuous in VALUE at the cut (rc.hpp says it is
    // deliberately not continuous in slope).
    const double eps = 1e-9;
    CHECK_CLOSE(v->fc(t_of(LI6_VMC_C0_Q_CUT_FM + eps)),
                v->fc(t_of(LI6_VMC_C0_Q_CUT_FM - eps)), 1e-7);
    // The provenance names the edge, the file and the continuation.
    const std::string p = v->provenance();
    CHECK(p.find("vmc-ft") != std::string::npos);
    CHECK(p.find("NO C0 ZERO") != std::string::npos);
    CHECK(p.find(VMC_LI6_DENSITY) != std::string::npos);
    CHECK(std::string(c0_shape_name(C0Shape::Ho)) == "ho");
    CHECK(std::string(c0_shape_name(C0Shape::VmcFt)) == "vmc-ft");
    // The DEFAULT is Ho: every published number was made with it.
    CHECK(HoSpin1FFOptions().c0_shape == C0Shape::Ho);
    CHECK(RcOptions().c0_shape == C0Shape::Ho);
    // ... and VmcFt is 6Li's OWN density, so it is REFUSED for any other ion
    // rather than silently handed over (the same rule as the measured-moment
    // block).  A deuteron WITH its moments supplied is otherwise legal --
    // T10 builds exactly that -- so this is the only thing stopping it.
    HoSpin1FFOptions d;
    d.a_fm = 1.0; d.alpha = 0.1; d.mu_n = 0.857; d.q_fm2 = 0.2857;
    d.fm_qz_fm = 1.3; d.fm_b_fm = 1.85;
    CHECK_NOTHROW(HoSpin1FF::for_ion(DEUTERON(), d));
    d.c0_shape = C0Shape::VmcFt;
    CHECK_THROWS_AS(HoSpin1FF::for_ion(DEUTERON(), d), std::runtime_error);
  }

  SUBCASE("for_ion THROWS for any ion with no measured-moment block") {
    CHECK_THROWS_AS(HoSpin1FF::for_ion(LI7()), std::runtime_error);
    CHECK_THROWS_AS(HoSpin1FF::for_ion(DEUTERON()), std::runtime_error);
    CHECK_THROWS_AS(HoSpin1FF::for_ion(PROTON()), std::runtime_error);
    // ... but a deuteron WITH its moments supplied is fine (T10 needs it).
    CHECK_NOTHROW(deuteron_ff());
    // The provenance string names the model and the flags the run must repeat.
    const std::string p = ff->provenance();
    CHECK(p.find("MEASURED") != std::string::npos);
    CHECK(p.find("NOT VMC") != std::string::npos);
  }
}

// ---------------------------------------------------------------- T7, T8

TEST_CASE("T7: the Born denominator IS the library's own dsigma_unpol") {
  // POLRAD Eq. (9) at Q_N = 0, P_L = 0:
  //   d^2 sigma/(dx dy) = (4 pi alpha^2 S/Q^4) [ x y^2 F1 + (1-y) F2 ] .
  // `x s dsigma_unpol(x, Q2, s)` is IDENTICALLY that -- with F_L = F2 - 2x F1
  // the two brackets are the same polynomial -- which is why rc.cpp divides by
  // the library's own Born and does not transcribe Eq. (9) a second time
  // (docs/CONVENTIONS.md: no physics number defined twice).
  const InclusiveKernel& k = rc_sampler().kernel();
  const double s = rc_sampler().s();
  for (double x : {1e-3, 0.01, 0.1, 0.3, 0.6}) {
    for (double q2 : {1.5, 5.0, 20.0, 100.0}) {
      const double y = q2 / (x * s);
      if (!(y > 1e-4 && y < 0.985)) continue;
      const SFTables t = k.tables(x, q2, true);
      const double want = 4.0 * kPi * ALPHA_EM * ALPHA_EM * s / (q2 * q2) *
                          (x * y * y * t.f1 + (1.0 - y) * t.f2) * GEV2_TO_PB;
      CHECK_CLOSE(x * s * k.dsigma_unpol(x, q2, s), want, 1e-12);
    }
  }
}

namespace {

/// THE LITERAL EQ. (38) TRANSCRIPTIONS, written HERE and not in rc.cpp.
/// They go through the SHARED quadrature `polrad_tpeak_quadrature`, so what
/// they gate is the INTEGRAND, the fractions, the sign and the prefactor --
/// which is exactly where a missing A = 6, a swapped (3/4) <-> (4/3) or a
/// dropped Q_N/6 would hide, and where no 5 %-tolerance test could see it.

/// POLRAD Eq. (38), sigma_u^d  (polrad2t.tex:905-914; adgh:8893-8894).
double eq38_sigma_u_d(const Spin1ElasticFF& ff, double x_a, double y,
                      double s_a, double m_a, int n_eta) {
  return polrad_tpeak_quadrature(
      [&](double eta, double xt, double) {
        const double t = 4.0 * m_a * m_a * eta;
        const double fc = ff.fc(t), fm = ff.fm(t), fq = ff.fq(t);
        return (fc * fc + (8.0 / 9.0) * fq * fq * eta * eta +
                (2.0 / 3.0) * fm * fm * eta) * xt
               - (2.0 / 3.0) * (1.0 + eta) * fm * fm;
      },
      x_a, y, s_a, polrad_eta_limits(x_a, y, s_a, m_a), n_eta);
}

/// POLRAD Eq. (38), sigma_q^d  (polrad2t.tex:915-926; adgh:8985-8987).
/// (3/4) x_A^2 inside the Xt bracket, (4/3) eta_A F_q in the last line -- the
/// PDF renders these two SWAPPED, the .tex and the FORTRAN do not.
double eq38_sigma_q_d(const Spin1ElasticFF& ff, double x_a, double y,
                      double s_a, double m_a, int n_eta) {
  return polrad_tpeak_quadrature(
      [&](double eta, double xt, double x1) {
        const double t = 4.0 * m_a * m_a * eta;
        const double fc = ff.fc(t), fm = ff.fm(t), fq = ff.fq(t);
        return (1.0 + eta + ((3.0 / 4.0) * x_a * x_a - eta) * xt) * fm * fm
               - (xt * x1 / (1.0 + eta)) * fq *
                     (3.0 * fc + 3.0 * eta * fm + eta * fq)
               - 2.0 * eta * xt * fq *
                     (4.0 * fc - 3.0 * x_a * fm + (4.0 / 3.0) * eta * fq);
      },
      x_a, y, s_a, polrad_eta_limits(x_a, y, s_a, m_a), n_eta);
}

/// POLRAD Eq. (38), sigma_u^C -- the SPIN-0 line, with the explicit Z^2 that
/// pins the charge normalisation, and with Y_+ (CORRECTED: the paper prints
/// Y_-, POLRAD's own `apptai` applies Y_+ to every unpolarised entry).
double eq38_sigma_u_c(double z, const Spin0FF& ff, double x_a, double y,
                      double s_a, double m_a, int n_eta) {
  return polrad_tpeak_quadrature(
      [&](double eta, double xt, double) {
        const double f = ff.shape(4.0 * m_a * m_a * eta);
        return z * z * xt * f * f;
      },
      x_a, y, s_a, polrad_eta_limits(x_a, y, s_a, m_a), n_eta);
}

/// The SPIN-1/2 line, written through Eq. (40)'s (F_1, F_2) rather than the
/// (G_E, G_M) the code uses -- an independent algebraic path to the same
/// number (F_1 + F_2 = G_M, F_1^2 + eta F_2^2 = (G_E^2 + eta G_M^2)/(1+eta)).
double eq38_sigma_u_nucleon(bool proton, double x, double y, double s,
                            double m_n, int n_eta) {
  return polrad_tpeak_quadrature(
      [&](double eta, double xt, double) {
        const NucleonFF f = nucleon_ff(4.0 * m_n * m_n * eta);
        const double ge = proton ? f.ge_p : f.ge_n;
        const double gm = proton ? f.gm_p : f.gm_n;
        const double f1 = (ge + eta * gm) / (1.0 + eta);
        const double f2 = (gm - ge) / (1.0 + eta);
        return (f1 * f1 + eta * f2 * f2) * xt - (f1 + f2) * (f1 + f2);
      },
      x, y, s, polrad_eta_limits(x, y, s, m_n), n_eta);
}


// --------------------------------------------------------------------------
// AN INDEPENDENT QUADRATURE, so that T8(a)'s 1e-10 gate covers more than the
// integrand.  `eq38_sigma_*_d` above go through the SAME
// `polrad_tpeak_quadrature` as rc.cpp, so the leading minus, the
// -(alpha^3/S_A) Y_+ prefactor and the eta limits were gated only loosely.
// This one shares NOTHING with rc.cpp: a plain composite trapezoid in
// ln(eta_A), its own prefactor, and eta_1/eta_2 rebuilt from POLRAD Eq. (14)
// (tau_min/max = (S_x -+ sqrt(lambda_Q))/(2M^2)) and Eq. (17)
// (R_el = (S_xA - Q^2)/(1 + tau_A), eta_A = (Q^2 + R_el tau_A)/(4 M_A^2))
// rather than from the rationalised closed form `polrad_eta_limits` uses.
double trapz_eq38(const Spin1ElasticFF& ff, double x_a, double y, double s_a,
                  double m_a, int n, bool tensor) {
  const double m2 = m_a * m_a;
  const double q2 = x_a * y * s_a;
  const double sx = y * s_a;                       // S_x = S_A - X_A
  const double sq = std::sqrt(sx * sx + 4.0 * m2 * q2);      // sqrt(lambda_Q)
  const double tau1 = (sx - sq) / (2.0 * m2);      // Eq. (14)
  const double tau2 = (sx + sq) / (2.0 * m2);
  auto eta_of = [&](double tau) {                  // Eq. (17)
    return (q2 + (sx - q2) / (1.0 + tau) * tau) / (4.0 * m2);
  };
  const double a = std::log(eta_of(tau1)), b = std::log(eta_of(tau2));
  double acc = 0.0;
  for (int i = 0; i <= n; ++i) {
    const double eta = std::exp(a + (b - a) * i / n);
    const double x1 = x_a * x_a + 4.0 * x_a * eta - 4.0 * eta;   // Eq. (39)
    const double xt = x1 / (2.0 * eta * x_a * x_a);
    const double t = 4.0 * m2 * eta;
    const double fc = ff.fc(t), fm = ff.fm(t), fq = ff.fq(t);
    const double f =
        tensor ? ((1.0 + eta + (0.75 * x_a * x_a - eta) * xt) * fm * fm
                  - xt * x1 / (1.0 + eta) * fq *
                        (3.0 * fc + 3.0 * eta * fm + eta * fq)
                  - 2.0 * eta * xt * fq *
                        (4.0 * fc - 3.0 * x_a * fm + (4.0 / 3.0) * eta * fq))
               : ((fc * fc + (8.0 / 9.0) * eta * eta * fq * fq +
                   (2.0 / 3.0) * eta * fm * fm) * xt
                  - (2.0 / 3.0) * (1.0 + eta) * fm * fm);
    acc += ((i == 0 || i == n) ? 0.5 : 1.0) * f;
  }
  acc *= (b - a) / n;
  const double yp = (1.0 + (1.0 - y) * (1.0 - y)) / (1.0 - y);   // Y_+
  return -(ALPHA_EM * ALPHA_EM * ALPHA_EM / s_a) * yp * acc;     // Eq. (18)'s minus
}

// --------------------------------------------------------------------------
// THE LEADING-LOG s- AND p-PEAKS now SHIP: `ll_radiator`, `rosenbluth_spin1`,
// `rosenbluth_nucleon`, `dsigma_el_dq2`, `ll_peaks_spin1` and `ll_peaks_qe`
// moved from this anonymous namespace into src/core/rc.cpp (declared in
// rc.hpp beside `polrad_sigma_qe_u`) so that T8(c) below gates the SHIPPED
// code path and not a construction private to this file, and so that
// `RcTailModel::TPeakPlusLL` can tabulate them.  Their derivation is in the
// rc.hpp block that declares them.

}  // namespace

TEST_CASE("T8(0): the tail ADDS events -- sigma^el_U > 0 and w_tail >= 1") {
  // polrad_transcription_check.md sec. 8, the BLOCKING flag.  X_1 is decreasing
  // in eta_A and vanishes at the ultrarelativistic eta_min, so Xt < 0 over
  // essentially the whole range and Eq. (38) AS PRINTED returns a NEGATIVE
  // sigma_u -- a radiative tail that REMOVES events from the DIS bin and
  // therefore ANTI-dilutes A_zz.  Eq. (18) (polrad2t.tex:580) carries an
  // explicit leading minus that Eq. (38) does not print; rc.cpp applies it.
  // This is the single cheapest guard against shipping that sign error, and
  // it runs BEFORE the normalisation gates on purpose.
  const RcModel& m = rc_model();
  const double s = rc_sampler().s();
  for (double x : {0.005, 0.01, 0.1, 0.3}) {
    for (double q2 : {2.0, 5.0, 20.0}) {
      const double y = q2 / (x * s);
      if (!(y > 0.0041 && y < 0.984)) continue;
      const RcModel::TailTriple tt = m.tail_sigma_at(x, q2);
      CHECK_MESSAGE(tt.u > 0.0, "sigma^el_U <= 0 at x = " << x << ", Q2 = "
                                                          << q2);
      CHECK_MESSAGE(tt.qe > 0.0, "sigma^q_U <= 0 at x = " << x << ", Q2 = "
                                                          << q2);
      CHECK(m.tail_ratio_at(x, q2, 0.0) > 0.0);
      // ... and at both tensor extremes of a spin-1 fill (P_zz = +1, -2).
      for (double q_n : {0.0, 1.0, -2.0}) {
        CHECK_MESSAGE(m.tail_ratio_at(x, q2, q_n) > 0.0,
                      "w_tail < 1 at x = " << x << ", Q2 = " << q2
                                           << ", Q_N = " << q_n);
      }
    }
  }
  // Every entry of the built table obeys it too.
  const std::vector<double>& su = m.sigma_tail_u();
  const std::vector<double>& sq = m.sigma_tail_qe();
  const std::size_t ny = m.table_y().size();
  REQUIRE(su.size() == m.table_x().size() * ny);
  for (std::size_t i = 0; i < m.table_x().size(); ++i) {
    // The QRT is in NUCLEON invariants, so at the grid's last edge x = 1 it
    // vanishes identically (eta_min -> infinity: that IS the elastic point).
    // The ELASTIC tail does not: its x_A = x/A is 1/6 there.
    const bool elastic_point = m.table_x()[i] >= 1.0;
    for (std::size_t j = 0; j < ny; ++j) {
      CHECK(su[i * ny + j] > 0.0);
      // (numerically 0 up to the eta-limit rounding: 2e-33 at worst)
      if (elastic_point) CHECK(sq[i * ny + j] < 1e-30);
      else CHECK(sq[i * ny + j] > 0.0);
    }
  }
}

TEST_CASE("T8(a): the DEUTERON gate -- Eq. (38) transcribed IN the test") {
  const std::shared_ptr<HoSpin1FF> ff = deuteron_ff();
  const double m_d = DEUTERON().mass();
  const int n = 128;
  for (double s_n : {2.0 * PROTON_MASS * 27.6, 3980.0}) {
    for (double x : {0.005, 0.05, 0.2}) {
      for (double y : {0.1, 0.5, 0.9}) {
        // THIS library's map (rc.cpp, "THE NUCLEAR MAP"): S_A = A s and
        // x_A = x/A, both EXACT for its collider variables.  POLRAD's
        // fixed-target S_A = S m_A/m_p, x_A = x m_p/M_A differs by
        // A m_p/M_A = 1.001 for the deuteron; the gate below is an identity
        // between two transcriptions and holds under either.
        const double s_a = 2.0 * s_n;
        const double x_a = x / 2.0;
        CHECK_CLOSE(polrad_sigma_el_u(*ff, x_a, y, s_a, m_d, n),
                    eq38_sigma_u_d(*ff, x_a, y, s_a, m_d, n), 1e-10);
        CHECK_CLOSE(polrad_sigma_el_t(*ff, x_a, y, s_a, m_d, n),
                    eq38_sigma_q_d(*ff, x_a, y, s_a, m_d, n), 1e-10);
      }
    }
  }

  SUBCASE("an INDEPENDENT quadrature: prefactor, leading minus and limits") {
    // `eq38_sigma_*_d` share `polrad_tpeak_quadrature` with rc.cpp, so the
    // 1e-10 gate above covers the INTEGRAND only.  `trapz_eq38` shares
    // nothing: its own composite trapezoid in ln(eta_A), its own
    // -(alpha^3/S_A) Y_+ prefactor and its own eta_1/eta_2 from Eqs. (14) and
    // (17).  A 10 % error in the prefactor, a dropped leading minus or a
    // wrong eta_hi would pass every other gate in this file and fails here.
    for (double x : {0.005, 0.05, 0.2}) {
      for (double y : {0.1, 0.5, 0.9}) {
        const double s_a = 2.0 * 2.0 * PROTON_MASS * 27.6;
        const double x_a = x / 2.0;
        CHECK_CLOSE(polrad_sigma_el_u(*ff, x_a, y, s_a, m_d, n),
                    trapz_eq38(*ff, x_a, y, s_a, m_d, 20000, false), 1e-6);
        CHECK_CLOSE(polrad_sigma_el_t(*ff, x_a, y, s_a, m_d, n),
                    trapz_eq38(*ff, x_a, y, s_a, m_d, 20000, true), 1e-6);
      }
    }
    // ... and the 6Li production point, at the map the model actually uses.
    {
      const std::shared_ptr<HoSpin1FF> f6 = HoSpin1FF::for_ion(LI6());
      const double m6 = LI6().mass(), s6 = 6.0 * 3980.0;
      CHECK_CLOSE(polrad_sigma_el_u(*f6, 0.01 / 6.0, 0.5, s6, m6, n),
                  trapz_eq38(*f6, 0.01 / 6.0, 0.5, s6, m6, 20000, false),
                  1e-6);
    }
  }

  SUBCASE("X_1 vanishes at the ULTRARELATIVISTIC eta_min, and Xt < 0 above it") {
    // The identity that fixes X_1's LAST term as LINEAR in eta_A: a
    // transcription with 4 eta_A^2 fails it.
    for (double x_a : {0.001, 0.02, 0.2, 0.5}) {
      const double em = polrad_eta_min_ur(x_a);
      CHECK_CLOSE(x_a * x_a + 4.0 * x_a * em - 4.0 * em, 0.0, 1e-12);
      CHECK(x_a * x_a + 4.0 * x_a * (1.5 * em) - 4.0 * (1.5 * em) < 0.0);
    }
  }

  SUBCASE("t_min = M_A^2 x_A^2/(1-x_A) ~= (x M_A/A)^2, INDEPENDENT of A") {
    // If it scaled with A the nuclear map would be wrong somewhere.  With THIS
    // library's map (x_A = x/A, exact for its collider variables -- see
    // rc.cpp's "THE NUCLEAR MAP") the cancellation is exact up to
    // 1/(1 - x_A) and up to M_A/A vs m_p, which is 0.5 % for 6Li.
    for (const Ion* ion : {&DEUTERON(), &LI6(), &LI7()}) {
      const double m_a = ion->mass();
      const double a = static_cast<double>(ion->A);
      for (double x : {0.01, 0.1, 0.3}) {
        const double x_a = x / a;
        const double t_min = 4.0 * m_a * m_a * polrad_eta_min_ur(x_a);
        // EXACT for every A: M_A^2 x_A^2 = (x M_A/A)^2, and M_A/A is the
        // nucleus's own mass per nucleon (0.9336 GeV for 6Li).
        const double want = std::pow(x * m_a / a, 2);
        CHECK_CLOSE(t_min, want / (1.0 - x_a), 1e-12);
        // Nothing scales with A: t_min sits within 1.5 % of (x m_p)^2 for
        // EVERY ion here, the whole spread being the binding energy per
        // nucleon plus the 1/(1 - x_A).
        CHECK_CLOSE(t_min, std::pow(x * PROTON_MASS, 2),
                    0.02 + x_a / (1.0 - x_a));
      }
    }
    // The transcription check's own number: with the EXACT lower limit at
    // 6Li, x = 0.1 and an 11 GeV fixed target, t_min = 0.00880 GeV^2 against
    // (x m_p)^2 = 0.00880 -- the design's claim confirmed to 3 figures.  The
    // map change moves it by A m_p/M_A = 1.005 in x_A, i.e. 1 % in t_min,
    // which is inside the 3-figure statement.
    {
      const double m_a = LI6().mass();
      const double x_a = 0.1 / 6.0;
      const double s_a = 6.0 * 2.0 * PROTON_MASS * 11.0;
      const EtaLimits lim = polrad_eta_limits(x_a, 0.5, s_a, m_a);
      CHECK_CLOSE(4.0 * m_a * m_a * lim.lo, 0.00880, 1.5e-2);
      CHECK_CLOSE(std::pow(0.1 * PROTON_MASS, 2), 0.00880, 1e-3);
    }
  }

  SUBCASE("the EXACT eta limits, against the ultrarelativistic eta_min") {
    // POLRAD integrates eta1..eta2 (adgh:8611-8613) and keeps the paper's
    // eta_min = x^2/(4(1-x)) only as a COMMENTED-OUT line.  The ratio measured
    // for 6Li in the transcription check: 0.9996 at x = 0.01 y = 0.8, 0.983 at
    // x = 0.1 y = 0.5, 0.921 at x = 0.3 y = 0.3.
    const double m_a = LI6().mass();
    // The check did not print the beam energy those three ratios were measured
    // at; an 11 GeV fixed target (JLab) reproduces the last two to three
    // figures and the first to 7e-4, and no other energy reproduces any of
    // them, so that is the configuration pinned here.
    const double s_n = 2.0 * PROTON_MASS * 11.0;
    // This library's map: S_A = A s, x_A = x/A (rc.cpp, "THE NUCLEAR MAP").
    // It moves these three ratios by 1e-5 .. 8e-4 against POLRAD's
    // fixed-target one, i.e. inside the check's own quoted precision.
    const double s_a = 6.0 * s_n;
    struct Row { double x, y, ratio, rtol; };
    for (const Row& r : {Row{0.01, 0.8, 0.9996, 1e-3},
                         Row{0.1, 0.5, 0.983, 2e-3},
                         Row{0.3, 0.3, 0.921, 2e-3}}) {
      const double x_a = r.x / 6.0;
      const EtaLimits lim = polrad_eta_limits(x_a, r.y, s_a, m_a);
      CHECK(lim.hi > lim.lo);
      CHECK_CLOSE(lim.lo / polrad_eta_min_ur(x_a), r.ratio, r.rtol);
      // The EXACT limit is always BELOW the ultrarelativistic one, and the
      // integrand is proportional to Xt, which VANISHES at the u.r. eta_min
      // and does NOT at the exact one -- so this is not a rounding detail.
      CHECK(lim.lo < polrad_eta_min_ur(x_a));
    }
    // The upper limit is NOT infinity: eta_2 = S_x/(4 M_A^2).
    const double x_a = 0.1 / 6.0;
    CHECK_CLOSE(polrad_eta_limits(x_a, 0.5, s_a, m_a).hi,
                0.5 * s_a / (4.0 * m_a * m_a), 1e-12);
  }
}

TEST_CASE("T8(a'): the per-nucleon reduction is 1/A^2 -- ONE Jacobian and ONE 1/A") {
  // The single most likely error in this file, and no tolerance test can see
  // it.  Eq. (38) is the WHOLE-NUCLEUS d^2 sigma^el/(dx_A dy) -- the paper's
  // "sigma_1^el = (1/A) d^2 sigma/dx_A dy" (polrad2t.tex:579) is loose
  // notation, and POLRAD's own FORTRAN applies TWO factors to INT elu:
  // `ter = m_p/M_A` in `apptai` (adgh:8607-8620) AND 1/`tara` in the main
  // program (adgh:489, 497), against a PER-NUCLEON Born (adgh:4330-4339).
  // The independent derivation is Weizsacker-Williams x Compton, which
  // reproduces -Eq. (38) as (1 - x_A) 2 alpha^3 Z^2/(x_A^2 S_A) Y_+
  // INT (d eta/eta)(1 - eta_min/eta) F^2 -- with no 1/A in it.  So
  //
  //   per nucleon = (1/A) * (dx_A/dx) = (1/A)(1/A) = 1/A^2
  //
  // with this library's own map x_A = x/A.  The earlier m_p/M_A applied ONE
  // of the two and was 6.00x too large for 6Li.
  const RcModel& m = rc_model();
  const double s = rc_sampler().s();
  const double m_a = LI6().mass();
  const double a = static_cast<double>(LI6().A);
  const std::shared_ptr<HoSpin1FF> ff = HoSpin1FF::for_ion(LI6());
  for (double x : {0.01, 0.1}) {
    for (double q2 : {2.0, 8.0}) {
      const double y = q2 / (x * s);
      const double s_a = a * s;
      const double x_a = x / a;
      const double raw = polrad_sigma_el_u(*ff, x_a, y, s_a, m_a, 128);
      const RcModel::TailTriple tt = m.tail_sigma_at(x, q2);
      CHECK_CLOSE(tt.u / raw, 1.0 / (a * a), 1e-12);
      CHECK_CLOSE(tt.t / polrad_sigma_el_t(*ff, x_a, y, s_a, m_a, 128),
                  1.0 / (a * a), 1e-12);
    }
  }
  // ... and that number is 1/36 = 0.0277778, which is NOT m_p/M_A = 0.16750
  // (the factor v0 shipped, 6.00x too large), NOT 1/A = 0.16667 and NOT 1.
  CHECK_CLOSE(1.0 / (a * a), 0.027777777777777776, 1e-15);
  CHECK(std::fabs(PROTON_MASS / m_a - 1.0 / (a * a)) > 0.1);
  CHECK_CLOSE((PROTON_MASS / m_a) / (1.0 / (a * a)), 6.0301, 1e-4);

  SUBCASE("the ELASTIC and QUASI-ELASTIC per-nucleon factors are CONSISTENT") {
    // The bug the blocker was: the ERT carried m_p/M_A and the QRT carried
    // 1/A, so the two pieces of the SAME weight were normalised differently
    // and their ratio was wrong by a factor A.  Both now carry exactly one
    // 1/A on top of their own Jacobian (the QRT's is 1, it is already in
    // nucleon invariants), so dividing the model's own tail_sigma_at by the
    // raw quadratures must give the SAME 1/A on both.
    const double x = 0.05, q2 = 6.0;
    const double y = q2 / (x * s);
    const RcModel::TailTriple tt = m.tail_sigma_at(x, q2);
    const double raw_el = polrad_sigma_el_u(*ff, x / a, y, a * s, m_a, 128);
    const double raw_qe = polrad_sigma_qe_u(LI6().Z, LI6().N(), x, y, s,
                                            PROTON_MASS, 128,
                                            m.options().qe_kf_gev);
    // ERT: (1/A) x Jacobian(1/A).  QRT: (1/A) x Jacobian(1).
    CHECK_CLOSE((tt.u / raw_el) * a, 1.0 / a, 1e-12);
    CHECK_CLOSE(tt.qe / raw_qe, 1.0 / a, 1e-12);
  }
}

TEST_CASE("T8(b): the CARBON line -- the explicit Z^2 pins the charge") {
  // A spin-0 nucleus has Im_1 = 0 and Im_2 = Z^2 F^2, so the master formula
  // collapses to Eq. (38)'s carbon line.  NOTE it cannot be diffed against
  // POLRAD's `elu`: the shipped `approx` carbon branch (adgh:8903-8905) is
  // `elu = xxt*ff**2/eta` with `ffco` normalised to ff(0) = 1 and NO `tarz`
  // anywhere -- POLRAD's own closed-form carbon entry is low by Z^2 = 36.  The
  // paper (tex:928) and POLRAD's EXACT path (adgh:4251,
  // `f2=4.*amp2*tau*(tarz*ff)**2`) both carry Z^2, so the gate stands as
  // written (polrad_transcription_check.md sec. 7.2).
  const double m_c = 12.0 * M_NUCLEON;   // a stand-in mass; the gate is in Z
  const double s_a = 3980.0 * m_c / PROTON_MASS;
  const int n = 128;
  for (double z : {1.0, 6.0, 8.0}) {
    const Spin0FF ff(z, 2.0);
    for (double x : {0.01, 0.1}) {
      for (double y : {0.2, 0.7}) {
        const double x_a = x * PROTON_MASS / m_c;
        CHECK_CLOSE(polrad_sigma_el_u(ff, x_a, y, s_a, m_c, n),
                    eq38_sigma_u_c(z, ff, x_a, y, s_a, m_c, n), 1e-12);
        // THE SECOND CLAUSE OF THE DESIGN'S T9: sigma_q's integrand is
        // IDENTICALLY ZERO when F_m = F_q = 0.  Every term of Eq. (38)'s
        // sigma_q carries F_m or F_q (see `eq38_sigma_q_d` above), so the
        // quadrature sums exact zeros and `== 0.0` is the right assertion --
        // anything weaker would let a sign or term error leak F_c into the
        // tensor tail and still pass every gate in this file.
        CHECK(polrad_sigma_el_t(ff, x_a, y, s_a, m_c, n) == 0.0);
      }
    }
  }
  // The Z^2, explicitly: doubling Z quadruples the tail.
  const Spin0FF f1(1.0, 2.0), f3(3.0, 2.0);
  const double x_a = 0.05 * PROTON_MASS / m_c;
  CHECK_CLOSE(polrad_sigma_el_u(f3, x_a, 0.5, s_a, m_c, n) /
                  polrad_sigma_el_u(f1, x_a, 0.5, s_a, m_c, n),
              9.0, 1e-12);
  // ... and it is Y_+ and NOT the printed Y_-: Y_- vanishes as y -> 0, which
  // an unpolarised cross section cannot, so the ratio of the tail at two y at
  // fixed (x_A, S_A) must follow Y_+.
  const double ya = 0.15, yb = 0.75;
  auto yp = [](double y) { return (1.0 + (1.0 - y) * (1.0 - y)) / (1.0 - y); };
  auto ym = [](double y) { return y * (2.0 - y) / (1.0 - y); };
  const double ra = polrad_sigma_el_u(f3, x_a, ya, s_a, m_c, n);
  const double rb = polrad_sigma_el_u(f3, x_a, yb, s_a, m_c, n);
  // The eta integral moves a little with y through the exact limits, so this
  // is a 5 % statement -- but Y_+/Y_- differ by 6x here, so it is decisive.
  CHECK_CLOSE(ra / rb, yp(ya) / yp(yb), 5e-2);
  // Y_- would be off by a factor 5 here, so the two are not confusable.
  CHECK(std::fabs(ra / rb - ym(ya) / ym(yb)) / (ym(ya) / ym(yb)) > 1.0);
}

TEST_CASE("T8(c): the A = 1 spin-1/2 limit, by an INDEPENDENT algebraic path") {
  // POLRAD Eq. (38)'s proton entry is written in Eq. (40)'s Dirac/Pauli
  // (F_1, F_2); rc.cpp evaluates it in (G_E, G_M) through Eq. (A.5).  The two
  // are the same number only because F_1 + F_2 = G_M and
  // F_1^2 + eta F_2^2 = (G_E^2 + eta G_M^2)/(1+eta) -- so running both is a
  // real check on the spin-1/2 line, not a restatement of it.
  const double s = 2.0 * PROTON_MASS * 10.0;   // E = 10 GeV fixed target
  for (double x : {0.05, 0.3}) {
    for (double y : {0.2, 0.6}) {
      CHECK_CLOSE(polrad_sigma_qe_u(1, 0, x, y, s, PROTON_MASS, 128),
                  eq38_sigma_u_nucleon(true, x, y, s, PROTON_MASS, 128), 1e-12);
      CHECK_CLOSE(polrad_sigma_qe_u(0, 1, x, y, s, PROTON_MASS, 128),
                  eq38_sigma_u_nucleon(false, x, y, s, PROTON_MASS, 128),
                  1e-12);
    }
  }
  // -----------------------------------------------------------------------
  // THE ABSOLUTE GATE the design asked for, and the answer is NOT "the t-peak
  // is the dominant piece".
  //
  // The design's T8(c) asked for a 10 % comparison against the MO-TSAI exact
  // elastic tail (RMP 41 (1969) 205, App. B).  Mo-Tsai is not on arXiv and
  // was not obtained for this run, and a gate written against a formula
  // transcribed from memory would certify nothing.  What lands instead is the
  // review's own fallback: the LEADING-LOG s- and p-peaks of the SAME
  // observable, built above from the elastic cross section and the
  // Weizsacker-Williams lepton radiator and sharing no line with rc.cpp.  It
  // MEASURES what POLRAD sec. 2.1.3 B only asserts.
  SUBCASE("(c) the leading-log s-/p-peaks: WHERE the t-peak is the whole tail") {
    // First, the radiator machinery's own normalisation, against the
    // LABORATORY Rosenbluth by a completely different algebraic route:
    //   dsigma/dOmega = sigma_Mott [A + B tan^2(theta/2)],
    //   sigma_Mott = alpha^2 cos^2(theta/2)/(4 E^2 sin^4(theta/2)) (E'/E),
    //   dQ^2/du = 4 E'^2 with u = sin^2(theta/2)  =>  dOmega/dQ^2 = pi/E'^2 .
    {
      const double e = 10.0, q2 = 1.0, m = PROTON_MASS;
      const double u = q2 / (4.0 * e * e - 2.0 * e * q2 / m);   // Q^2 = 4E^2u/(1+2Eu/M)
      const double ep = e / (1.0 + 2.0 * e * u / m);
      const RosenbluthAB ab = rosenbluth_nucleon(true, q2);
      const double mott = ALPHA_EM * ALPHA_EM * (1.0 - u) /
                          (4.0 * e * e * u * u) * (ep / e);
      const double lab =
          mott * (ab.a + ab.b * u / (1.0 - u)) * kPi / (ep * ep);
      CHECK_CLOSE(dsigma_el_dq2(ab, q2, 2.0 * m * e, m), lab, 1e-10);
    }

    // ... and that `rosenbluth_spin1` is NOT a second definition of Eq. (38)'s
    // unpolarised integrand.  The two ARE the same object:
    //     A Xt - B/(2 eta)
    //   = [F_C^2 + (8/9) eta^2 F_Q^2 + (2/3) eta F_M^2] Xt
    //     - (2/3)(1 + eta) F_M^2 ,
    // so running the SHIPPED (A, B) through the SHIPPED quadrature has to
    // reproduce the SHIPPED `polrad_sigma_el_u` to round-off.  CONVENTIONS.md
    // forbids a physics number written twice; this is the gate that says it
    // was not.
    {
      const std::shared_ptr<HoSpin1FF> f6 = HoSpin1FF::for_ion(LI6());
      const double m6 = LI6().mass(), s6 = 6.0 * 3980.0;
      for (const double xa : {0.01 / 6.0, 0.1 / 6.0, 0.3 / 6.0}) {
        for (const double y : {0.3, 0.7, 0.95}) {
          const double via_ab = polrad_tpeak_quadrature(
              [&](double eta, double xt, double) {
                const RosenbluthAB ab =
                    rosenbluth_spin1(*f6, 4.0 * m6 * m6 * eta, m6);
                return ab.a * xt - ab.b / (2.0 * eta);
              },
              xa, y, s6, polrad_eta_limits(xa, y, s6, m6), 128);
          CHECK_CLOSE(via_ab, polrad_sigma_el_u(*f6, xa, y, s6, m6, 128),
                      1e-14);
        }
      }
    }

    // (i) 6Li at EIC config 1, Q^2 >= 20 GeV^2, AT FOUR POINTS WITH y >= 0.5.
    //     Read what the four rows below are: the ELASTIC half of POLRAD's
    //     claim, plus four SAMPLED points of the quasi-elastic half.  They are
    //     NOT a window statement, and the comment that used to stand here --
    //     "POLRAD's claim CONFIRMED -- the t-peak is > 99 % of t + s + p, for
    //     BOTH tails" -- read them as one.  NARROWED 2026-09-04; the
    //     assertions are unchanged and still hold AT THESE POINTS.  The
    //     mechanism is real: the nuclear charge form factor is dead by
    //     t ~ 0.25 GeV^2 and the nucleon dipole by t ~ 1, while the s-/p-peaks
    //     sit at t ~ z Q^2 ~ Q^2.
    //
    //     WHAT THE WINDOW ACTUALLY DOES, measured cell by cell over the
    //     101 x 77 CLI grid's 1356 accepted cells at Q^2 >= 20 and y <= 0.9,
    //     in this subcase's own unsuppressed kf_gev = 0 convention:
    //       * ELASTIC t/(t+s+p): minimum 1.0000.  The claim HOLDS, and more
    //         strongly than "> 0.999" -- 6Li's coherent form factor is dead at
    //         Q'^2 ~ Q^2 >= 20, so the elastic s+p is nothing there.
    //       * QUASI-ELASTIC: minimum 1.6458e-04, and 275 of the 1356 cells are
    //         BELOW 0.96 (282 at the shipped k_F = 0.169 GeV).
    //       * TOTAL: 346 of the 1356 are below 0.99.
    //     MOST of the failure is at y -> 0 -- but not all of it, and the
    //     earlier wording "every failure is at y -> 0, which no row here
    //     reaches" was WITHDRAWN 2026-09-04 as overstated: 21 of the 346
    //     total-tail failures and 9 of the 275 quasi-elastic ones sit at
    //     y >= 0.1, and FIVE total-tail failures sit at y = 0.845-0.892
    //     (x = 0.0060-0.0087) -- the y -> 1 end of the window, right beside
    //     the (x = 0.01, y = 0.9) row this subcase does assert.  One
    //     quasi-elastic failure (x = 0.0060256, y = 0.84484, f_qe = 0.9576)
    //     is a direct counterexample.  The hedged form -- 318 of the 331 at
    //     y < 0.1 -- is the one every doc site carries and the one to quote.
    //     NOTE ON GATING, corrected 2026-09-04: this census is NOT gated by
    //     T8(d)(i).  That subcase runs the test file's own 40 x 24 sampler,
    //     compares only the combined tail_ratio, and never forms t/(t+s+p) at
    //     all -- so nothing in the suite asserts the elastic minimum 1.0000,
    //     the 275/282 quasi-elastic counts, or 346 of 1356.  These four
    //     bullets are a MEASUREMENT RECORD on the 101 x 77 CLI grid, not a
    //     gated claim, and they are labelled as such rather than re-asserted
    //     here: re-running them in the suite would multiply its cost for a
    //     statement no shipped number depends on.
    {
      const std::shared_ptr<HoSpin1FF> f6 = HoSpin1FF::for_ion(LI6());
      const double m6 = LI6().mass(), s6n = 3980.0, s6 = 6.0 * s6n;
      struct Row { double x, y; };
      for (const Row& r : {Row{0.01, 0.5}, Row{0.01, 0.9}, Row{0.1, 0.5},
                           Row{0.3, 0.5}}) {
        const double x_a = r.x / 6.0;
        const double tp = polrad_sigma_el_u(*f6, x_a, r.y, s6, m6, 256);
        const LlPeaks sp = ll_peaks_spin1(*f6, x_a, r.y, s6, m6);
        const double f_el = tp / (tp + sp.s + sp.p);
        const double qtp = polrad_sigma_qe_u(3, 3, r.x, r.y, s6n, PROTON_MASS,
                                             256, 0.0);
        const LlPeaks qsp = ll_peaks_qe(3, 3, r.x, r.y, s6n);
        const double f_qe = qtp / (qtp + qsp.s + qsp.p);
        MESSAGE("6Li EIC x = " << r.x << ", y = " << r.y << ", Q2 = "
                               << r.x * r.y * s6n
                               << ": t-peak fraction elastic " << f_el
                               << ", quasi-elastic " << f_qe);
        CHECK(f_el > 0.99);
        if (r.x * r.y * s6n > 20.0) CHECK(f_qe > 0.96);
      }
    }

    // (ii) ... and where it is NOT.  These two rows are the honest half of the
    //      gate: they FAIL a "t-peak is the dominant piece" claim and are
    //      asserted so that a later change cannot quietly pretend otherwise.
    //      (a) the PROTON at E = 10 GeV, x = 0.3, Q^2 = 2 GeV^2 -- the design's
    //          own Mo-Tsai point.  The dipole is still alive at Q'^2 = 1.4
    //          GeV^2, so the s-peak is ~26x the t-peak.
    {
      const double s = 2.0 * PROTON_MASS * 10.0, x = 0.3, q2 = 2.0;
      const double y = q2 / (x * s);
      const double tp = polrad_sigma_qe_u(1, 0, x, y, s, PROTON_MASS, 256, 0.0);
      const LlPeaks sp = ll_peaks_qe(1, 0, x, y, s);
      MESSAGE("proton E = 10 GeV, x = 0.3, Q2 = 2: (s+p)/t = "
              << (sp.s + sp.p) / tp);
      CHECK((sp.s + sp.p) / tp > 10.0);
    }
    //      (b) the HERMES deuteron point, x = 0.012, y = 0.85 (their lowest-x
    //          bin sits at y ~ 0.8-0.85 because of the Q^2 > 0.5 cut).  The
    //          t-peak is 23 % of the leading-log total, i.e. LOW BY A FACTOR
    //          ~4.4 -- which is what makes rc_tail a LOWER BOUND there and is
    //          consistent with HERMES's "almost 50 % of the statistics in the
    //          lowest-x bin" radiative background against the ~1.4 % the
    //          t-peak alone gives.
    {
      const std::shared_ptr<HoSpin1FF> fd = deuteron_ff();
      const double m = DEUTERON().mass(), sn = 2.0 * PROTON_MASS * 27.6;
      const double x = 0.012, y = 0.85, x_a = x / 2.0, s_a = 2.0 * sn;
      const double te = polrad_sigma_el_u(*fd, x_a, y, s_a, m, 256) / 4.0;
      const LlPeaks se = ll_peaks_spin1(*fd, x_a, y, s_a, m);
      const double tq = polrad_sigma_qe_u(1, 1, x, y, sn, PROTON_MASS, 256,
                                          0.0) / 2.0;
      const LlPeaks sq = ll_peaks_qe(1, 1, x, y, sn);
      const double t_tot = te + tq;
      const double all =
          t_tot + (se.s + se.p) / 4.0 + (sq.s + sq.p) / 2.0;
      MESSAGE("HERMES deuteron x = 0.012, y = 0.85: t-peak fraction "
              << t_tot / all << " (elastic "
              << te / (te + (se.s + se.p) / 4.0) << ", quasi-elastic "
              << tq / (tq + (sq.s + sq.p) / 2.0) << ")");
      CHECK(t_tot / all < 0.30);
      CHECK(t_tot / all > 0.10);
    }
  }
}

TEST_CASE("T8(c'): POLRAD's OWN compiled deuteron numbers") {
  // polrad_transcription_check.md sec. 8 evaluated POLRAD's `elu`/`elq` with
  // POLRAD's own `ffdeu`, its own limits and its own prefactors, at three
  // points.  Reproduced here with THIS repository's deuteron form factor,
  // which is a DIFFERENT shape (the HO stand-in of T10, not `ffdeu`) -- so the
  // magnitudes are gated loosely and the SIGN STRUCTURE tightly.  What this
  // catches is a wrong power of A, a wrong prefactor, a missing Y_+ or a
  // dropped leading minus: any of those moves the magnitude by 6x or more.
  const std::shared_ptr<HoSpin1FF> ff = deuteron_ff();
  const double m_d = DEUTERON().mass();
  struct Row { double e, x, y, sigma_u_polrad, ratio_polrad; };
  // `sigma_u_polrad` is the check's own number with its `ter = m_p/M_d`
  // divided back out, i.e. Eq. (38) itself, and with the leading minus applied
  // so that it is POSITIVE here.
  const Row rows[] = {
      {27.6, 0.050, 0.60, 2.352e-05 / (PROTON_MASS / 1.8756280),  0.1060},
      {27.6, 0.012, 0.50, 1.001e-03 / (PROTON_MASS / 1.8756280),  0.0642},
      {11.0, 0.200, 0.50, 1.570e-07 / (PROTON_MASS / 1.8756280), -0.1165},
  };
  for (const Row& r : rows) {
    const double s_n = 2.0 * PROTON_MASS * r.e;
    const double s_a = 2.0 * s_n;      // S_A = A s (rc.cpp's map)
    const double x_a = r.x / 2.0;      // x_A = x/A
    const double su = polrad_sigma_el_u(*ff, x_a, r.y, s_a, m_d, 256);
    const double st = polrad_sigma_el_t(*ff, x_a, r.y, s_a, m_d, 256);
    CHECK(su > 0.0);                                    // the sign gate again
    // Measured with this repository's HO stand-in: 0.883, 0.983, 0.390 of
    // POLRAD's own numbers at the three points (the map change moved them by
    // 0.1-0.3 %).  The x = 0.2 point is the
    // loose one because its t-peak starts at q ~ 1 fm^-1, exactly where a
    // two-parameter HO and `ffdeu` differ most; the two low-x points, where
    // the t-peak lives at q < 0.3 fm^-1 and the shape hardly matters, agree to
    // 2-12 %, which is the real normalisation statement.
    CHECK_MESSAGE(su > 0.25 * r.sigma_u_polrad,
                  "sigma_u is " << su / r.sigma_u_polrad
                                << "x POLRAD's at x = " << r.x);
    CHECK_MESSAGE(su < 4.0 * r.sigma_u_polrad,
                  "sigma_u is " << su / r.sigma_u_polrad
                                << "x POLRAD's at x = " << r.x);
    if (r.x < 0.1) CHECK_CLOSE(su, r.sigma_u_polrad, 0.2);
    // sigma_q/sigma_u CHANGES SIGN between x = 0.05 and x = 0.20.  That is
    // the second reason sigma^el_T is a SEPARATE table and never a scale
    // factor on sigma^el_U -- a scale factor cannot change sign with x.
    CHECK_MESSAGE((st / su) * r.ratio_polrad > 0.0,
                  "sigma_q/sigma_u sign disagrees with POLRAD at x = " << r.x
                      << ": " << st / su << " vs " << r.ratio_polrad);
  }
}

TEST_CASE("T8': the quasi-elastic tail, Eq. (44)") {
  const double s = 3980.0;
  const int n = 128;
  SUBCASE("sigma^q_U = Z sigma_u^p[G_E^p,G_M^p] + N sigma_u^p[G_E^n,G_M^n]") {
    for (double x : {0.01, 0.1, 0.3}) {
      for (double y : {0.1, 0.5, 0.9}) {
        const double want =
            3.0 * eq38_sigma_u_nucleon(true, x, y, s, PROTON_MASS, n) +
            3.0 * eq38_sigma_u_nucleon(false, x, y, s, PROTON_MASS, n);
        CHECK_CLOSE(polrad_sigma_qe_u(LI6().Z, LI6().N(), x, y, s,
                                      PROTON_MASS, n),
                    want, 1e-12);
      }
    }
    // (Z, N) = (3, 3) for 6Li, and the neutron term is NOT negligible -- it is
    // the G_M^n one, mu_n = -1.913.
    CHECK(LI6().Z == 3);
    CHECK(LI6().N() == 3);
    // measured 0.1156 -- the G_M^n term, mu_n = -1.913, carries it.
    CHECK(eq38_sigma_u_nucleon(false, 0.1, 0.5, s, PROTON_MASS, n) >
          0.05 * eq38_sigma_u_nucleon(true, 0.1, 0.5, s, PROTON_MASS, n));
  }

  SUBCASE("qe_suppression = 0 reproduces the elastic-only tail BIT FOR BIT") {
    RcOptions off_qe;
    off_qe.qe_suppression = 0.0;
    RcOptions no_qe;
    no_qe.with_qe_tail = false;
    const RunPlan plan = tensor_thirds_plan(0.0, 0.6);
    const RcModel a(RcMode::TensorBand, off_qe, rc_sampler_ptr(), plan,
                    Channel::Inclusive, LI6());
    const RcModel b(RcMode::TensorBand, no_qe, rc_sampler_ptr(), plan,
                    Channel::Inclusive, LI6());
    for (double x : {0.01, 0.1}) {
      for (double q2 : {3.0, 9.0}) {
        for (double q_n : {0.0, 1.0, -2.0}) {
          CHECK(a.tail_ratio_at(x, q2, q_n) == b.tail_ratio_at(x, q2, q_n));
        }
      }
    }
    // ... and the QRT really is a large part of the default tail.
    const RcModel full(RcMode::TensorBand, RcOptions(), rc_sampler_ptr(), plan,
                       Channel::Inclusive, LI6());
    CHECK(full.tail_ratio_at(0.1, 5.0, 0.0) >
          1.2 * a.tail_ratio_at(0.1, 5.0, 0.0));
  }
}

TEST_CASE("B3: the POLARISED quasi-elastic stand-in, RcOptions::qe_tensor_scale") {
  // The polarised quasi-elastic tail is not computed anywhere -- POLRAD has no
  // tensor partner to Eq. (44) and no such calculation exists for an A = 6
  // spin-1 nucleus.  `qe_tensor_scale` lends it the ELASTIC tail's own
  // sigma^el_T/sigma^el_U, which is a BORROWED magnitude and not a derived
  // bound (rc.hpp, RcOptions::qe_tensor_scale).  What is testable is the
  // ALGEBRA and the plumbing, and that is what this case gates; the physics
  // caveats live in the header and in phase_B_numbers.md sec. B3.
  const RunPlan plan = tensor_thirds_plan(0.0, 0.6);
  const RcModel& base = rc_model();
  // In-support (x, Q^2): the tables refuse to extrapolate, and
  // y = Q^2/(x s) has to land inside [0.004, 0.985] at s = 3980.
  struct XQ { double x, q2; };
  const XQ pts[] = {{0.01, 3.0}, {0.01, 5.0}, {0.01, 9.0},
                    {0.1, 3.0},  {0.1, 5.0},  {0.1, 9.0},
                    {0.3, 5.0},  {0.3, 9.0}};

  SUBCASE("the default is 0 and the shipped tail does not move") {
    CHECK(RcOptions().qe_tensor_scale == 0.0);
    RcOptions zero;
    zero.qe_tensor_scale = 0.0;
    const RcModel m(RcMode::TensorBand, zero, rc_sampler_ptr(), plan,
                    Channel::Inclusive, LI6());
    for (const XQ& p : pts) {
      for (double q_n : {0.0, 1.0, -2.0}) {
        CHECK(m.tail_ratio_at(p.x, p.q2, q_n) ==
              base.tail_ratio_at(p.x, p.q2, q_n));
      }
    }
  }

  SUBCASE("it is EXACTLY LINEAR in the scale -- no T12 needed") {
    // Unlike `fq_scale`, which enters sigma^el_T quadratically through F_q and
    // therefore has to be RUN at each edge (T12), this term is one product.
    // One run rescales, and the tests say so rather than leaving a reader to
    // assume it.
    RcOptions o1, o2, o4;
    o1.qe_tensor_scale = 1.0;
    o2.qe_tensor_scale = 2.0;
    o4.qe_tensor_scale = 4.0;
    const RcModel m1(RcMode::TensorBand, o1, rc_sampler_ptr(), plan,
                     Channel::Inclusive, LI6());
    const RcModel m2(RcMode::TensorBand, o2, rc_sampler_ptr(), plan,
                     Channel::Inclusive, LI6());
    const RcModel m4(RcMode::TensorBand, o4, rc_sampler_ptr(), plan,
                     Channel::Inclusive, LI6());
    for (const XQ& p : pts) {
      const double base1 = base.tail_ratio_at(p.x, p.q2, 1.0);
      const double d1 = m1.tail_ratio_at(p.x, p.q2, 1.0) - base1;
      const double d2 = m2.tail_ratio_at(p.x, p.q2, 1.0) - base1;
      const double d4 = m4.tail_ratio_at(p.x, p.q2, 1.0) - base1;
      REQUIRE(std::fabs(d1) > 0.0);
      CHECK_CLOSE(d2 / d1, 2.0, 1e-10);
      CHECK_CLOSE(d4 / d1, 4.0, 1e-10);
      // ... and it is LINEAR IN q_n too, like every other tensor term here,
      // which is what keeps T14's Eq. (43) identity exact.
      const double a =
          m1.tail_ratio_at(p.x, p.q2, 1.0) - m1.tail_ratio_at(p.x, p.q2, 0.0);
      const double b =
          m1.tail_ratio_at(p.x, p.q2, -2.0) - m1.tail_ratio_at(p.x, p.q2, 0.0);
      CHECK_CLOSE(b / a, -2.0, 1e-9);
    }
  }

  SUBCASE("what it adds IS (q_n/6) scale (sigma^el_T/sigma^el_U) sigma^q_U") {
    // Gated AT A TABLE NODE, where the interpolation is exact, so the identity
    // is the algebra and not an interpolation tolerance.  The added tensor
    // term divided by the ELASTIC tensor term must be kappa_qe sigma^q_U /
    // sigma^el_U -- both read from `tail_sigma_at`, i.e. from the quadrature
    // and not from the ratio under test.
    RcOptions o;
    o.qe_tensor_scale = 1.0;
    const RcModel m(RcMode::TensorBand, o, rc_sampler_ptr(), plan,
                    Channel::Inclusive, LI6());
    const std::vector<double>& nx = m.table_x();
    const std::vector<double>& ny = m.table_y();
    REQUIRE(nx.size() > 8);
    REQUIRE(ny.size() > 8);
    const double s_pn = rc_sampler().s();
    int checked = 0;
    for (std::size_t i = 4; i < nx.size(); i += 7) {
      for (std::size_t j = 3; j < ny.size(); j += 9) {
        const double x = nx[i];
        const double q2 = x * s_pn * ny[j];
        const RcModel::TailTriple t = m.tail_sigma_at(x, q2);
        if (!(t.u > 0.0) || !(t.qe > 0.0)) continue;
        // Skip the top corner, where sigma^q_U/sigma^el_U reaches 1e7: both
        // sides of the identity are then DIFFERENCES of two nearly equal tail
        // ratios and the cancellation, not the algebra, sets the agreement
        // (measured 1.4e-8 there against 1e-12 in the bulk).  The identity is
        // gated where it can be gated.
        if (t.qe > 1e4 * t.u) continue;
        const double el_t = base.tail_ratio_at(x, q2, 1.0) -
                            base.tail_ratio_at(x, q2, 0.0);
        const double added = m.tail_ratio_at(x, q2, 1.0) -
                             base.tail_ratio_at(x, q2, 1.0);
        if (!(std::fabs(el_t) > 1e-300)) continue;
        CHECK_CLOSE(added / el_t, t.qe / t.u, 1e-8);
        ++checked;
      }
    }
    CHECK(checked >= 8);
  }

  SUBCASE("it rides INSIDE qe_suppression -- no parent, no stand-in") {
    // `qe_suppression` is documented as a flat multiplier on the WHOLE
    // quasi-elastic tail.  A tensor stand-in that survived
    // `qe_suppression = 0` would be a tail with no unpolarised parent, which
    // is not a band edge but a contradiction.
    RcOptions o;
    o.qe_suppression = 0.0;
    o.qe_tensor_scale = 5.0;
    RcOptions ref;
    ref.qe_suppression = 0.0;
    const RcModel m(RcMode::TensorBand, o, rc_sampler_ptr(), plan,
                    Channel::Inclusive, LI6());
    const RcModel r(RcMode::TensorBand, ref, rc_sampler_ptr(), plan,
                    Channel::Inclusive, LI6());
    for (double x : {0.01, 0.1, 0.3}) {
      for (double q_n : {0.0, 1.0, -2.0}) {
        CHECK(m.tail_ratio_at(x, 5.0, q_n) == r.tail_ratio_at(x, 5.0, q_n));
      }
    }
    // ... and HALF the suppression is half the stand-in, exactly.
    RcOptions h, f;
    h.qe_suppression = 0.5; h.qe_tensor_scale = 1.0;
    f.qe_suppression = 1.0; f.qe_tensor_scale = 1.0;
    const RcModel mh(RcMode::TensorBand, h, rc_sampler_ptr(), plan,
                     Channel::Inclusive, LI6());
    const RcModel mf(RcMode::TensorBand, f, rc_sampler_ptr(), plan,
                     Channel::Inclusive, LI6());
    RcOptions hb;
    hb.qe_suppression = 0.5;
    const RcModel bh(RcMode::TensorBand, hb, rc_sampler_ptr(), plan,
                     Channel::Inclusive, LI6());
    for (double x : {0.01, 0.1}) {
      const double half = mh.tail_ratio_at(x, 5.0, 1.0) -
                          bh.tail_ratio_at(x, 5.0, 1.0);
      const double full = mf.tail_ratio_at(x, 5.0, 1.0) -
                          base.tail_ratio_at(x, 5.0, 1.0);
      REQUIRE(std::fabs(full) > 0.0);
      CHECK_CLOSE(half / full, 0.5, 1e-9);
    }
  }

  SUBCASE("a knob that did not run is REFUSED, and a negative one too") {
    RcOptions bad;
    bad.with_qe_tail = false;
    bad.qe_tensor_scale = 1.0;
    CHECK_THROWS_AS(RcModel(RcMode::TensorBand, bad, rc_sampler_ptr(), plan,
                            Channel::Inclusive, LI6()),
                    std::runtime_error);
    // ... but zero with the tail off is fine: nothing is recorded that did
    // not run.
    RcOptions ok;
    ok.with_qe_tail = false;
    CHECK_NOTHROW(RcModel(RcMode::TensorBand, ok, rc_sampler_ptr(), plan,
                          Channel::Inclusive, LI6()));
    RcOptions neg;
    neg.qe_tensor_scale = -1e-12;
    CHECK_THROWS_AS(RcModel(RcMode::TensorBand, neg, rc_sampler_ptr(), plan,
                            Channel::Inclusive, LI6()),
                    std::runtime_error);
  }

  SUBCASE("the SIZE of what was unpriced, printed") {
    // The point of the whole exercise: at x = 0.30 the quasi-elastic tail is
    // 99.9 % of `rc_tail`, so lending it the elastic tensor fraction changes
    // the tensor part of the tail by ~1e3, while at x = 0.01 -- where the
    // elastic tail still dominates -- it is a tens-of-percent effect.  These
    // are the numbers phase_B_numbers.md sec. B3 publishes.
    RcOptions o;
    o.qe_tensor_scale = 1.0;
    const RcModel m(RcMode::TensorBand, o, rc_sampler_ptr(), plan,
                    Channel::Inclusive, LI6());
    struct Row { double x, lo, hi; };
    const Row rows[] = {{0.01, 1.1, 1.6}, {0.1, 3.0, 5.0}, {0.3, 3e2, 3e3}};
    for (const Row& r : rows) {
      const double el = base.tail_ratio_at(r.x, 5.0, 1.0) -
                        base.tail_ratio_at(r.x, 5.0, 0.0);
      const double both = m.tail_ratio_at(r.x, 5.0, 1.0) -
                          m.tail_ratio_at(r.x, 5.0, 0.0);
      const double f = both / el;
      MESSAGE("B3 x = " << r.x << " Q2 = 5: tensor tail x" << f
                        << " (r_T " << el << " -> " << both << ")");
      CHECK(f > r.lo);
      CHECK(f < r.hi);
    }
  }
}

TEST_CASE("B1: the s-/p-peak tensor stand-in, RcOptions::sp_tensor_scale") {
  // `RcTailModel::TPeakPlusLL` puts the leading-log s- and p-peaks in the
  // UNPOLARISED numerator only, because POLRAD's Eq. (38) supplies no tensor
  // s/p peak and inventing one from the leading log would be an uncited second
  // definition.  (Eq. (18) + Eq. (A.4) DOES supply one and
  // `RcTailModel::PolradFull` computes it -- T19 -- which is why this scale is
  // refused there for the OPPOSITE reason: the term RAN.  Unqualified, the
  // sentence is true of Eq. (38) and false of the paper.)  That left the
  // tensor fraction of that piece at EXACTLY ZERO -- a choice, not a
  // measurement.
  // `sp_tensor_scale` prices it by lending the ELASTIC s-/p-peaks the elastic
  // t-peak's own sigma^el_T/sigma^el_U.  IT IS A BOUND WITH NO DERIVATION
  // (rc.hpp, RcOptions::sp_tensor_scale); what is testable is the ALGEBRA and
  // the plumbing, and that is what this case gates.
  const RunPlan plan = tensor_thirds_plan(0.0, 0.6);
  const RcModel& base = rc_model();          // the shipped TPeak default
  RcOptions ll_opt;
  ll_opt.tail_model = RcTailModel::TPeakPlusLL;
  const RcModel ll(RcMode::TensorBand, ll_opt, rc_sampler_ptr(), plan,
                   Channel::Inclusive, LI6());

  SUBCASE("the default is 0 and the SHIPPED tail cannot see it") {
    CHECK(RcOptions().sp_tensor_scale == 0.0);
    // Under `TPeak` the u_sp table is identically zero, which is why a
    // non-zero scale is REFUSED there rather than silently ignored.
    for (double x : {0.01, 0.1, 0.3}) {
      for (double q2 : {3.0, 5.0, 9.0}) {
        CHECK(base.tail_sigma_at(x, q2).u_sp == 0.0);
        CHECK(base.tail_sigma_at(x, q2).qe_sp == 0.0);
      }
    }
  }

  SUBCASE("a knob that did not run is REFUSED, and a negative one too") {
    RcOptions bad;                            // tail_model = TPeak
    bad.sp_tensor_scale = 1.0;
    CHECK_THROWS_AS(RcModel(RcMode::TensorBand, bad, rc_sampler_ptr(), plan,
                            Channel::Inclusive, LI6()),
                    std::runtime_error);
    // ... zero under TPeak is fine: nothing is recorded that did not run.
    CHECK_NOTHROW(RcModel(RcMode::TensorBand, RcOptions(), rc_sampler_ptr(),
                          plan, Channel::Inclusive, LI6()));
    // ... and under TPeakPlusLL it is accepted, because the peaks ran.
    RcOptions ok = ll_opt;
    ok.sp_tensor_scale = 1.0;
    CHECK_NOTHROW(RcModel(RcMode::TensorBand, ok, rc_sampler_ptr(), plan,
                          Channel::Inclusive, LI6()));
    RcOptions neg = ll_opt;
    neg.sp_tensor_scale = -1e-12;
    CHECK_THROWS_AS(RcModel(RcMode::TensorBand, neg, rc_sampler_ptr(), plan,
                            Channel::Inclusive, LI6()),
                    std::runtime_error);
  }

  SUBCASE("what it adds IS (q_n/6) scale (sigma^el_T/sigma^el_U) u_sp") {
    // Gated AT A TABLE NODE, where the interpolation is exact, so the identity
    // is the algebra and not an interpolation tolerance.  The added tensor
    // term divided by the elastic t-peak's own tensor term must be
    // u_sp/sigma^el_U -- both read from `tail_sigma_at`, i.e. from the
    // quadrature and not from the ratio under test.
    RcOptions o = ll_opt;
    o.sp_tensor_scale = 1.0;
    const RcModel m(RcMode::TensorBand, o, rc_sampler_ptr(), plan,
                    Channel::Inclusive, LI6());
    const std::vector<double>& nx = m.table_x();
    const std::vector<double>& ny = m.table_y();
    REQUIRE(nx.size() > 8);
    REQUIRE(ny.size() > 8);
    const double s_pn = rc_sampler().s();
    int checked = 0, alive = 0, dead = 0;
    double worst_rel = 0.0, wrx = 0.0, wry = 0.0;
    for (std::size_t i = 0; i < nx.size(); ++i) {
      for (std::size_t j = 0; j < ny.size(); ++j) {
        const double x = nx[i];
        const double q2 = x * s_pn * ny[j];
        const RcModel::TailTriple t = m.tail_sigma_at(x, q2);
        if (!(t.u > 0.0)) continue;
        (t.u_sp > 0.0) ? ++alive : ++dead;
        if (!(t.u_sp > 0.0)) continue;
        const double el_t = ll.tail_ratio_at(x, q2, 1.0) -
                            ll.tail_ratio_at(x, q2, 0.0);
        const double added =
            m.tail_ratio_at(x, q2, 1.0) - ll.tail_ratio_at(x, q2, 1.0);
        if (!(std::fabs(el_t) > 1e-300)) continue;
        // TWO CANCELLATION GUARDS, both the B3 subcase's guard in a new
        // dress.  `added` and `el_t` are each DIFFERENCES of two tail ratios
        // whose common size is the WHOLE numerator, so the identity loses one
        // digit for every decade by which the whole numerator exceeds the
        // piece being isolated: once where sigma^q_U swamps sigma^el_U, and
        // again where u_sp is a tiny fraction of sigma^el_U.  Gated where it
        // can be gated; the ALGEBRA is the same everywhere.
        if (t.qe > 1e2 * t.u) continue;
        if (t.u_sp > 1e4 * t.u || t.u_sp < 1e-3 * t.u) continue;
        const double got = added / el_t, want = t.u_sp / t.u;
        const double rel = std::fabs(got - want) / std::fabs(want);
        if (rel > worst_rel) { worst_rel = rel; wrx = x; wry = ny[j]; }
        // 1e-5 and not B3's 1e-8 because of those two guards' own subject:
        // even inside them the isolation costs digits.  MEASURED worst on
        // this grid 2026-09-06: 3.60359e-06 at x = 0.0251189, y = 0.004 (the
        // y-grid floor, where sigma^q_U is closest to its 1e2 ceiling), over
        // 270 gated nodes.
        CHECK_CLOSE(got, want, 1e-5);
        ++checked;
      }
    }
    MESSAGE("B1: u_sp alive at " << alive << " of " << (alive + dead)
                                 << " nodes; identity gated at " << checked
                                 << "; worst rel " << worst_rel << " at x = "
                                 << wrx << ", y = " << wry);
    CHECK(checked >= 8);
    // AND THE COHERENT s/p PEAK IS DEAD OVER MOST OF THE GRID.  That is the
    // whole of rc.hpp point (3), measured: the s-/p-peak's own elastic vertex
    // sits at Q'^2 ~ Q^2, where 6Li's coherent form factor is numerically
    // zero, while the t-peak's sits at t ~ t_min, where it is alive.
    CHECK(dead > 0);
    CHECK(alive > 0);
  }

  SUBCASE("it is EXACTLY LINEAR in the scale and in q_n") {
    RcOptions o1 = ll_opt, o2 = ll_opt, o4 = ll_opt;
    o1.sp_tensor_scale = 1.0;
    o2.sp_tensor_scale = 2.0;
    o4.sp_tensor_scale = 4.0;
    const RcModel m1(RcMode::TensorBand, o1, rc_sampler_ptr(), plan,
                     Channel::Inclusive, LI6());
    const RcModel m2(RcMode::TensorBand, o2, rc_sampler_ptr(), plan,
                     Channel::Inclusive, LI6());
    const RcModel m4(RcMode::TensorBand, o4, rc_sampler_ptr(), plan,
                     Channel::Inclusive, LI6());
    const std::vector<double>& nx = m1.table_x();
    const std::vector<double>& ny = m1.table_y();
    const double s_pn = rc_sampler().s();
    int checked = 0;
    for (std::size_t i = 0; i < nx.size() && checked < 12; ++i) {
      for (std::size_t j = 0; j < ny.size() && checked < 12; ++j) {
        const double x = nx[i];
        const double q2 = x * s_pn * ny[j];
        if (!(m1.tail_sigma_at(x, q2).u_sp > 0.0)) continue;
        const double b1 = ll.tail_ratio_at(x, q2, 1.0);
        const double d1 = m1.tail_ratio_at(x, q2, 1.0) - b1;
        const double d2 = m2.tail_ratio_at(x, q2, 1.0) - b1;
        const double d4 = m4.tail_ratio_at(x, q2, 1.0) - b1;
        if (!(std::fabs(d1) > 0.0)) continue;
        CHECK_CLOSE(d2 / d1, 2.0, 1e-10);
        CHECK_CLOSE(d4 / d1, 4.0, 1e-10);
        const double a = m1.tail_ratio_at(x, q2, 1.0) -
                         m1.tail_ratio_at(x, q2, 0.0);
        const double b = m1.tail_ratio_at(x, q2, -2.0) -
                         m1.tail_ratio_at(x, q2, 0.0);
        CHECK_CLOSE(b / a, -2.0, 1e-9);
        ++checked;
      }
    }
    CHECK(checked >= 8);
  }

  SUBCASE("it is OUTSIDE qe_suppression -- u_sp is the ELASTIC column") {
    // The mirror image of B3's "it rides INSIDE qe_suppression": that knob is
    // a flat multiplier on the QUASI-ELASTIC tail, and `u_sp` is not part of
    // it, so switching the quasi-elastic tail off must leave this stand-in
    // standing.
    RcOptions a = ll_opt, b = ll_opt;
    a.qe_suppression = 0.0;
    b.qe_suppression = 0.0;
    a.sp_tensor_scale = 1.0;
    const RcModel ma(RcMode::TensorBand, a, rc_sampler_ptr(), plan,
                     Channel::Inclusive, LI6());
    const RcModel mb(RcMode::TensorBand, b, rc_sampler_ptr(), plan,
                     Channel::Inclusive, LI6());
    const std::vector<double>& nx = ma.table_x();
    const std::vector<double>& ny = ma.table_y();
    const double s_pn = rc_sampler().s();
    int moved = 0;
    for (std::size_t i = 0; i < nx.size(); ++i) {
      for (std::size_t j = 0; j < ny.size(); ++j) {
        const double x = nx[i];
        const double q2 = x * s_pn * ny[j];
        if (!(ma.tail_sigma_at(x, q2).u_sp > 0.0)) continue;
        if (ma.tail_ratio_at(x, q2, 1.0) != mb.tail_ratio_at(x, q2, 1.0)) {
          ++moved;
        }
        // ... and the UNPOLARISED tail is untouched at every node, q_n = 0.
        CHECK(ma.tail_ratio_at(x, q2, 0.0) == mb.tail_ratio_at(x, q2, 0.0));
      }
    }
    CHECK(moved > 0);
  }

  SUBCASE("THE PRICE, and where the bound is EMPTY -- printed") {
    // The measurement phase_B_numbers.md sec. B1 publishes.  At the three
    // standard points the term is EXACTLY ZERO, because the coherent s/p
    // vertex Q'^2 = z_s Q^2 = 4.37 / 4.94 / 4.98 GeV^2 is four decades above
    // the t-peak's own t ~ t_min ~ 8.7e-05 GeV^2 and 6Li's coherent form
    // factor is dead there (F_c = -3.8e-45 against +2.99).  The bound is not
    // conservative at those points; it is EMPTY there.
    RcOptions o = ll_opt;
    o.sp_tensor_scale = 1.0;
    const RcModel m(RcMode::TensorBand, o, rc_sampler_ptr(), plan,
                    Channel::Inclusive, LI6());
    for (double x : {0.01, 0.1, 0.3}) {
      // The priced tail ratio is EQUAL BIT FOR BIT to the unpriced one: the
      // added term is so far below the numerator's ulp that it is not a small
      // correction, it is no correction at all.
      CHECK(m.tail_ratio_at(x, 5.0, 1.0) == ll.tail_ratio_at(x, 5.0, 1.0));
      const RcModel::TailTriple t = m.tail_sigma_at(x, 5.0);
      MESSAGE("B1 x = " << x << " Q2 = 5: u_sp = " << t.u_sp
                        << ", sigma^el_U = " << t.u
                        << ", u_sp/sigma^el_U = " << (t.u_sp / t.u));
      // NOT `== 0.0`: the coherent s/p peak is not switched off there, it is
      // KILLED BY THE FORM FACTOR -- 70+ decades below the t-peak it is being
      // lent a tensor fraction from.  That is rc.hpp point (3), measured.
      CHECK(t.u_sp < 1e-60 * t.u);
    }
    // ... and where it is ALIVE -- low x, y -> 1 -- it is not small.  The
    // largest tensor-tail multiplier on this grid is reported, never asserted
    // as a physics claim: it is a PRICE TAG on an omission.
    const std::vector<double>& nx = m.table_x();
    const std::vector<double>& ny = m.table_y();
    const double s_pn = rc_sampler().s();
    double worst = 1.0, wx = 0.0, wy = 0.0;
    for (std::size_t i = 0; i < nx.size(); ++i) {
      for (std::size_t j = 0; j < ny.size(); ++j) {
        const double x = nx[i];
        const double q2 = x * s_pn * ny[j];
        if (!(m.tail_sigma_at(x, q2).u_sp > 0.0)) continue;
        const double a = ll.tail_ratio_at(x, q2, 1.0) -
                         ll.tail_ratio_at(x, q2, 0.0);
        const double b = m.tail_ratio_at(x, q2, 1.0) -
                         m.tail_ratio_at(x, q2, 0.0);
        if (!(std::fabs(a) > 1e-300)) continue;
        if (std::fabs(b / a) > std::fabs(worst)) {
          worst = b / a;
          wx = x;
          wy = ny[j];
        }
      }
    }
    MESSAGE("B1 sp_tensor_scale = 1: tensor tail x" << worst << " at x = "
            << wx << ", y = " << wy << "; EXACTLY 1 at x = 0.01/0.1/0.3, "
            "Q2 = 5 (the coherent s/p vertex is dead there)");
    CHECK(worst > 1.0);
  }
}

TEST_CASE("Q9: the QUASI-ELASTIC Pauli suppression S(q), POLRAD Eq. (44)") {
  // Design Q9, closed.  v0 shipped S_E = S_M = 1, which after the per-nucleon
  // fix is the S = 1 EDGE of a band on what is now the DOMINANT piece of
  // `rc_tail`.  `ffquas` (adgh:5605-5613) codes de Forest-Walecka's Fermi-gas
  // factor, and the t-peak's reach down to t_min ~ (x M_N)^2 makes it big.
  SUBCASE("the factor itself: shape, continuity and the two limits") {
    const double kf = RC_QE_KF_GEV;
    CHECK(pauli_suppression(0.0, kf) == 1.0);            // guarded, not 0/0
    CHECK_CLOSE(pauli_suppression(2.0 * kf, kf), 1.0, 1e-15);   // 3/2 - 1/2
    // S is STATIONARY at u = 2 (dS/du = 0.75 - 3u^2/16 = 0 there), so it
    // joins the q >= 2 k_F plateau smoothly; step back far enough to see it.
    CHECK(pauli_suppression(2.0 * kf * 0.9, kf) < 1.0);
    CHECK(pauli_suppression(10.0, kf) == 1.0);           // q >> 2 k_F
    CHECK_CLOSE(pauli_suppression(kf, kf), 0.75 - 1.0 / 16.0, 1e-15);
    // monotone rising on (0, 2 k_F)
    double prev = 0.0;
    for (int i = 1; i <= 200; ++i) {
      const double v = pauli_suppression(2.0 * kf * i / 200.0, kf);
      CHECK(v >= prev);
      prev = v;
    }
    // kf <= 0 disables it -- that is how `polrad_sigma_qe_u`'s own default
    // reproduces POLRAD Eq. (44) at S_E = S_M = 1.
    CHECK(pauli_suppression(0.05, 0.0) == 1.0);
  }

  SUBCASE("QRT(S)/QRT(1) at 6Li: 0.47 at x = 0.01, 0.87 at 0.1, 1.00 at 0.3") {
    // The measured answer to Q9.  These are the numbers that make the S = 1
    // default 15-60 % high on the now-dominant piece of the tail at x <= 0.1.
    const double s = 3980.0;
    struct Row { double x, r169, r221, rtol; };
    for (const Row& r : {Row{0.01, 0.4721, 0.4053, 2e-3},
                         Row{0.1,  0.8683, 0.7831, 2e-3},
                         Row{0.3,  1.0000, 0.9903, 2e-3}}) {
      const double y = 5.0 / (r.x * s);
      const double one = polrad_sigma_qe_u(3, 3, r.x, y, s, PROTON_MASS, 256,
                                           0.0);
      CHECK_CLOSE(polrad_sigma_qe_u(3, 3, r.x, y, s, PROTON_MASS, 256, 0.169) /
                      one, r.r169, r.rtol);
      CHECK_CLOSE(polrad_sigma_qe_u(3, 3, r.x, y, s, PROTON_MASS, 256, 0.221) /
                      one, r.r221, r.rtol);
    }
    // A LARGER k_F suppresses MORE (more of the integral sits below 2 k_F).
    const double y = 5.0 / (0.01 * s);
    CHECK(polrad_sigma_qe_u(3, 3, 0.01, y, s, PROTON_MASS, 256, 0.221) <
          polrad_sigma_qe_u(3, 3, 0.01, y, s, PROTON_MASS, 256, 0.169));
  }

  SUBCASE("it is ON by default in RcModel, and qe_suppression is the BAND") {
    // `qe_kf_gev` is the physics; `qe_suppression` stays a flat multiplier on
    // top of it, so the two knobs compose rather than duplicate.
    RcOptions on;                       // qe_kf_gev = RC_QE_KF_GEV
    CHECK(on.qe_kf_gev == RC_QE_KF_GEV);
    RcOptions off;
    off.qe_kf_gev = 0.0;
    const RcModel m_on(RcMode::TensorBand, on, rc_sampler_ptr(),
                       tensor_thirds_plan(0.0, 0.6), Channel::Inclusive, LI6());
    const RcModel m_off(RcMode::TensorBand, off, rc_sampler_ptr(),
                        tensor_thirds_plan(0.0, 0.6), Channel::Inclusive,
                        LI6());
    const double x = 0.01, q2 = 5.0;
    const double qe_on = m_on.tail_sigma_at(x, q2).qe;
    const double qe_off = m_off.tail_sigma_at(x, q2).qe;
    CHECK(qe_on < qe_off);
    CHECK_CLOSE(qe_on / qe_off, 0.4721, 2e-3);
    // ... and the ELASTIC tail is untouched by it.
    CHECK_CLOSE(m_on.tail_sigma_at(x, q2).u, m_off.tail_sigma_at(x, q2).u,
                1e-14);
  }
}

TEST_CASE("the tail tables REFUSE a scenario whose y_max is above the ceiling") {
  // The latent run abort: `locate` throws rather than extrapolate, and
  // `Scenario::validate` does not bound y_max, so before this guard the FIRST
  // event with y > 0.9995 killed the whole run from inside `Pipeline::event`.
  // Fail at CONSTRUCTION instead.
  InclusiveSampler::GridSpec g;
  g.nx = 8;
  g.nq2 = 6;
  Scenario sc = generator_scenario(Scenario());
  sc.y_max = 1.0;
  std::shared_ptr<const InclusiveSampler> wide =
      std::make_shared<const InclusiveSampler>(band_kernel(false),
                                               default_configs("6Li")[1], sc, g,
                                               /*with_perp=*/false);
  CHECK_THROWS_WITH_AS(
      RcModel(RcMode::TensorBand, RcOptions(), wide,
              tensor_thirds_plan(0.0, 0.6), Channel::Inclusive, LI6()),
      doctest::Contains("above the radiative tail's y ceiling"),
      std::runtime_error);
  // ... and the same run with the tail off is fine: only the tail tables have
  // the ceiling, the band does not.
  RcOptions band_only;
  band_only.with_tail = false;
  const RcModel ok(RcMode::TensorBand, band_only, wide,
                   tensor_thirds_plan(0.0, 0.6), Channel::Inclusive, LI6());
  CHECK(ok.applies());
  CHECK(!ok.tail_applies());
  // The generator scenario itself is comfortably under it.
  CHECK(generator_scenario(Scenario()).y_max < RC_TAIL_Y_CEILING);
}

TEST_CASE("T9: the Q_N = 0 Rosenbluth limit of Eq. (A.4), and the Q_N part") {
  // THIS TEST CARRIED `doctest::skip(true)` FROM 2026-09-02 TO 2026-09-06 --
  // the ONE unconditional skip in the suite -- because Eq. (A.4) at Q_N != 0
  // was never evaluated: `RcTailModel::PolradFull` threw.  It is implemented
  // now (`polrad_full_sigma_el`), Eq. (A.4) is `polrad_im_el_spin1`, and the
  // design's T9 is exactly this: Im^el_2|_0 == A(Q^2), Im^el_1|_0 == B(Q^2)/2,
  // Im^el_{5,6,7,8}|_0 == 0.
  //
  // WHY IT IS A REAL GATE AND NOT A RESTATEMENT.  `rosenbluth_spin1` is
  // written from the (A, B) side -- it is what the collinear peaks need, a
  // cross section at a shifted Q'^2 -- while `polrad_im_el_spin1` is
  // polrad2t.tex:2698-2709 transcribed as printed.  Getting the Q_N split
  // wrong silently doubles the UNPOLARISED tail into the tensor one, which no
  // unpolarised test can see (design sec. 1.4.6's own warning).
  const std::shared_ptr<HoSpin1FF> ff = HoSpin1FF::for_ion(LI6());
  const double m_a = LI6().mass();

  SUBCASE("(a) Im_1|_0 = B/2, Im_2|_0 = A, Im_{5..8}|_0 = 0") {
    for (double t : {1e-4, 1e-3, 1e-2, 0.05, 0.2, 1.0, 5.0}) {
      const PolradIm im = polrad_im_el_spin1(*ff, t, m_a);
      const RosenbluthAB ab = rosenbluth_spin1(*ff, t, m_a);
      CHECK_CLOSE(im.u[2], ab.a, 1e-15);
      CHECK_CLOSE(im.u[1], 0.5 * ab.b, 1e-15);
      // EXACT zeros, not "small": Eq. (A.4)'s Im_5..8 carry an explicit
      // overall Q_N, so their Q_N = 0 halves are structurally absent and
      // anything weaker would let a tensor term leak into the unpolarised
      // tail and still pass.
      for (int i : {3, 4, 5, 6, 7, 8}) CHECK(im.u[i] == 0.0);
      // ... and the vector (P_N) entries are deliberately left at zero:
      // Im_3, Im_4 reach every tail formula through m M P_L alone.
      CHECK(im.t[3] == 0.0);
      CHECK(im.t[4] == 0.0);
    }
  }

  SUBCASE("(b) the Q_N part, against a LITERAL Eq. (A.4) written here") {
    // polrad2t.tex:2698-2709, transcribed in this file rather than read from
    // rc.cpp, with the CORRECTED Im_6 (4/(1+eta), not 4 eta/(1+eta) --
    // polrad_transcription_check.md sec. 7.1, which the design carried wrong
    // from 2026-09-02 to 2026-09-04).
    for (double t : {1e-4, 1e-2, 0.2, 2.0}) {
      const double eta = t / (4.0 * m_a * m_a);
      const double fc = ff->fc(t), fm = ff->fm(t), fq = ff->fq(t);
      const double want1 = eta * eta * fm * fm;
      const double want2 = eta * fm * fm + (4.0 * eta * eta / (1.0 + eta)) *
                                               ((eta / 3.0) * fq + fc - fm) * fq;
      const double want5 = fm * fm / 4.0;
      const double want6 = (fm * fm + (4.0 / (1.0 + eta)) *
                                          ((eta / 3.0) * fq + fc + eta * fm) *
                                          fq) / 4.0;
      const double want7 = eta * (1.0 + eta) * fm * fm;
      const double want8 = -eta * fm * (fm + 2.0 * fq);
      const PolradIm im = polrad_im_el_spin1(*ff, t, m_a);
      CHECK_CLOSE(im.t[1], want1, 1e-14);
      CHECK_CLOSE(im.t[2], want2, 1e-14);
      CHECK_CLOSE(im.t[5], want5, 1e-14);
      CHECK_CLOSE(im.t[6], want6, 1e-14);
      CHECK_CLOSE(im.t[7], want7, 1e-14);
      CHECK_CLOSE(im.t[8], want8, 1e-14);
      // AND Im_1, Im_2 ARE NOT PURELY UNPOLARISED.  design sec. 1.4.6 warns
      // that splitting them wrong doubles the unpolarised tail into the
      // tensor one; these two `!=` are what would catch it.
      CHECK(im.t[1] != 0.0);
      CHECK(im.t[2] != 0.0);
    }
  }

  SUBCASE("(c) the SHIPPED tail reads this function, so the gate is not "
          "cosmetic") {
    // Scaling F_m alone must move BOTH columns of `polrad_full_sigma_el`,
    // because Im_1|_0, Im_2|_0 and Im_{1,2,5,6,7,8}^T all carry F_m.  If the
    // quadrature had its own copy of Eq. (A.4) this subcase would still pass
    // -- what it pins is that the copy does not exist: the Q_N = 0 column is
    // EXACTLY the Rosenbluth pair, so setting F_m = F_q = 0 must leave
    // sigma_u alive (F_c^2 survives) and sigma_q identically zero.
    HoSpin1FFOptions o;
    o.fq_scale = 0.0;
    o.tail_tensor_scale = 0.0;
    const std::shared_ptr<HoSpin1FF> flat = HoSpin1FF::for_ion(LI6(), o);
    for (double t : {1e-3, 0.1, 1.0}) {
      const PolradIm im = polrad_im_el_spin1(*flat, t, m_a);
      CHECK(im.u[2] > 0.0);
      for (int i = 1; i <= 8; ++i) CHECK(im.t[i] == 0.0);
      CHECK(im.u[1] == 0.0);
    }
    const double s_a = 6.0 * 3980.0;
    const PolradFullPair p =
        polrad_full_sigma_el(*flat, 0.01 / 6.0, 0.5, s_a, m_a, 128);
    CHECK(p.u > 0.0);
    CHECK(p.t == 0.0);
  }
}

// ------------------------------------------------------------ T1, T4, T17

TEST_CASE("B6: the A = 2 -> A = 6 transfer price, RcOptions::a_transfer_frac") {
  // design_C_tensor_rc.md Q8.  Every anchor `rc_delta` interpolates between is
  // a DEUTERON number; this knob widens the band by sqrt(1 + f^2) to price the
  // transfer, and does not correct it.
  const RunPlan plan = tensor_thirds_plan(0.0, 0.6);
  const double xs[] = {1e-4, 0.005, 0.01, 0.02, 0.05, 0.063,
                       0.1,  0.16,  0.3,  0.9};

  SUBCASE("the default is 0 and RcModel::delta is rc_delta BIT FOR BIT") {
    CHECK(RcOptions().a_transfer_frac == 0.0);
    const RcModel& m = rc_model();
    // Legitimate `==`: hypot(d, 0) returns |d| exactly and d >= 0 here, so the
    // shipped band is unchanged to the last bit -- which is what keeps the
    // reference JSONs and T3's own `==` anchors valid.
    for (double x : xs) CHECK(m.delta(x) == rc_delta(x));
    CHECK(m.delta(RC_X_HIGH) == RC_DELTA_HIGH_X);
    CHECK(m.delta(RC_X_LOW) == RC_DELTA_LOW_X);
  }

  SUBCASE("delta widens by EXACTLY sqrt(1 + f^2), and rc_delta itself does not") {
    for (double f : {0.25, 0.5, 1.0, 2.0}) {
      RcOptions o;
      o.a_transfer_frac = f;
      const RcModel m(RcMode::TensorBand, o, rc_sampler_ptr(), plan,
                      Channel::Inclusive, LI6());
      const double k = std::sqrt(1.0 + f * f);
      for (double x : xs) {
        CHECK_CLOSE(m.delta(x), k * rc_delta(x), 1e-14);
        // The PUBLISHED shape is untouched: the free function is what T3 gates
        // and what the Python `rc_delta` binding exposes.
        CHECK(rc_delta(x) == rc_delta(x, RC_DELTA_HIGH_X, RC_DELTA_LOW_X,
                                      RC_X_HIGH, RC_X_LOW));
      }
    }
  }

  SUBCASE("it moves the BAND and NOTHING else -- not the tail, not tau") {
    RcOptions o;
    o.a_transfer_frac = 1.0;
    const RcModel m(RcMode::TensorBand, o, rc_sampler_ptr(), plan,
                    Channel::Inclusive, LI6());
    const RcModel& base = rc_model();
    const InclusiveSampler& s = rc_sampler();
    for (std::size_t c = 0; c < s.n_cells(); c += 11) {
      const double x = s.x_cells()[c];
      const double q2 = s.q2_cells()[c];
      for (double q_n : {0.0, 1.0, -2.0})
        CHECK(m.tail_ratio_at(x, q2, q_n) == base.tail_ratio_at(x, q2, q_n));
      for (double mm : {1.0, 0.0, -1.0}) {
        const Event ev = inclusive_event(s, static_cast<int>(c), mm);
        CHECK(m.tensor_fraction(ev) == base.tensor_fraction(ev));
        // ... and the band it DOES move stays the same closed form, with
        // delta -> sqrt(2) delta.  w_lo + w_hi == 2 still, by T4's argument.
        const RcWeights w = m.weights(ev);
        const RcWeights b = base.weights(ev);
        const double tau = m.tensor_fraction(ev);
        // 1e-9, not 1e-12: `w.hi - 1` is a CANCELLATION (w.hi is within 1e-4
        // of 1 on most cells), so the ratio carries ~1e-16/|w.hi - 1| of
        // relative noise -- 3e-12 at the smallest tau here.  The identity is
        // exact; the subtraction is not.
        CHECK_CLOSE(w.hi - 1.0, std::sqrt(2.0) * (b.hi - 1.0), 1e-9);
        CHECK_CLOSE(w.lo - 1.0, std::sqrt(2.0) * (b.lo - 1.0), 1e-9);
        if (std::fabs(tau) > 1e-12) {
          CHECK_CLOSE(w.hi + w.lo, 2.0, 1e-12);
        }
      }
    }
  }

  SUBCASE("a negative fraction is REFUSED, not silently squared away") {
    RcOptions bad;
    bad.a_transfer_frac = -0.5;
    CHECK_THROWS_AS(RcModel(RcMode::TensorBand, bad, rc_sampler_ptr(), plan,
                            Channel::Inclusive, LI6()),
                    std::runtime_error);
  }
}

TEST_CASE("T1: the band IS [1 + (P_zz/2)A_zz(1+delta)]/[1 + (P_zz/2)A_zz]") {
  // The anchor that ties the code's tau to the PUBLISHED A_zz.  A_zz is
  // `asymmetries::azz` at its DEFAULT theta_m = 0 -- the LONGITUDINAL
  // asymmetry -- because all the axis geometry is in
  // P_zz^eff = 3 Q_NN P_2(cos theta_S); passing theta_m = theta_S as well
  // would apply P_2(cos theta_S) TWICE.
  const InclusiveSampler& s = rc_sampler();
  const RcModel& m = rc_model();
  const InclusiveKernel& k = s.kernel();
  std::size_t checked = 0;
  for (std::size_t c = 0; c < s.n_cells(); c += 7) {
    const double x = s.x_cells()[c];
    const double q2 = s.q2_cells()[c];
    const double y = q2 / (x * s.s());
    const SFTables t = k.tables(x, q2, true);
    const double azz_l = azz(t.b1, t.f1, t.f2, x, y, &t.b2);
    for (double mm : {1.0, 0.0, -1.0}) {
      const Event ev = inclusive_event(s, static_cast<int>(c), mm);
      const double p_zz = k.tensor_moments(mm).second;   // 3 Q_NN at theta_S = 0
      const double d = m.delta(x);
      const double den = 1.0 + 0.5 * p_zz * azz_l;
      REQUIRE(den > 0.0);
      const RcWeights w = m.weights(ev);
      CHECK_CLOSE(w.hi, (1.0 + 0.5 * p_zz * azz_l * (1.0 + d)) / den, 1e-13);
      CHECK_CLOSE(w.lo, (1.0 + 0.5 * p_zz * azz_l * (1.0 - d)) / den, 1e-13);
      // ... and tau itself is the supervisor's / HERMES's own coefficient.
      CHECK_CLOSE(m.tensor_fraction(ev), 0.5 * p_zz * azz_l / den, 1e-13);
      ++checked;
    }
  }
  CHECK(checked > 100);
}

TEST_CASE("T4: the band's signs") {
  const InclusiveSampler& s = rc_sampler();
  const RcModel& m = rc_model();
  const InclusiveKernel& k = s.kernel();

  SUBCASE("(a) w_lo + w_hi == 2.0, EXACTLY, on every event") {
    // Legitimate `==`: lo = 1 - d and hi = 1 + d with |d| <= 0.5 (delta <= 0.3,
    // |tau| < 1), and for such d the two rounding errors sum to less than half
    // an ulp of 2.
    for (std::size_t c = 0; c < s.n_cells(); c += 7) {
      for (double mm : {1.0, 0.0, -1.0}) {
        for (double phi : {0.0, 1.3, 4.0}) {
          const Event ev = inclusive_event(s, static_cast<int>(c), mm, phi);
          const RcWeights w = m.weights(ev);
          CHECK(w.lo + w.hi == 2.0);
        }
      }
    }
  }

  SUBCASE("(b) flipping m = +-1 -> 0 flips tau and scales it by P_zz W") {
    for (std::size_t c = 0; c < s.n_cells(); c += 31) {
      const Event e1 = inclusive_event(s, static_cast<int>(c), 1.0);
      const Event e0 = inclusive_event(s, static_cast<int>(c), 0.0);
      const double t1 = m.tensor_fraction(e1);
      const double t0 = m.tensor_fraction(e0);
      if (std::fabs(t1) < 1e-14) continue;          // a b1 zero crossing
      CHECK(t1 * t0 < 0.0);
      // [w_hi(0) - 1]/[w_hi(+1) - 1] = P_zz(0)/P_zz(+1) * W(+1)/W(0)
      //                              = -2 W(+1)/W(0)
      const double num1 = t1 / (1.0 - t1);   // = w_avg(m=+1)
      const double num0 = t0 / (1.0 - t0);   // = w_avg(m=0)
      CHECK_CLOSE(num0 / num1, -2.0, 1e-11);
      const double hi1 = m.weights(e1).hi - 1.0;
      const double hi0 = m.weights(e0).hi - 1.0;
      CHECK_CLOSE(hi0 / hi1, -2.0 * (1.0 + num1) / (1.0 + num0), 1e-11);
    }
  }

  SUBCASE("(c) the LITERATURE form: sign(w_hi - 1) == -sign(b1 P_zz)") {
    // Written WITHOUT naming TENSOR_LL_SIGN, so that recompiling with +1
    // FAILS it: Cosyn et al. Eq. (27) / HERMES have A_zz = -(2/3) b1/F1, i.e.
    // b1 > 0 means the m = 0 state has the LARGER rate.  An assertion phrased
    // as sign(TENSOR_LL_SIGN * b1 * P_zz) holds for EITHER value of the
    // constant and guards nothing.
    int seen = 0;
    for (std::size_t c = 0; c < s.n_cells(); c += 13) {
      const double x = s.x_cells()[c];
      const double q2 = s.q2_cells()[c];
      const SFTables t = k.tables(x, q2, true);
      if (std::fabs(t.b1) < 1e-12) continue;
      for (double mm : {1.0, 0.0}) {
        const Event ev = inclusive_event(s, static_cast<int>(c), mm);
        const double p_zz = k.tensor_moments(mm).second;
        const double hi = m.weights(ev).hi - 1.0;
        if (std::fabs(hi) < 1e-14) continue;
        CHECK_MESSAGE(hi * (t.b1 * p_zz) < 0.0,
                      "sign(w_hi - 1) is not -sign(b1 P_zz) at x = " << x
                          << ", Q2 = " << q2 << ", m = " << mm);
        ++seen;
      }
    }
    CHECK(seen > 20);
  }
}

TEST_CASE("T17: weighted mode -- slot 0 is PURE, slot 1+k is the MIXTURE") {
  const InclusiveSampler& s = rc_sampler();
  // Three PURE categories, so that slot 1+k and slot 0 can be compared at all.
  const std::vector<SpinCategory> pure = {
      SpinCategory("m+", 1.0, {1.0, 0.0, 0.0}),
      SpinCategory("m0", 1.0, {0.0, 1.0, 0.0}),
      SpinCategory("m-", 1.0, {0.0, 0.0, 1.0})};
  const RcModel mp(RcMode::TensorBand, RcOptions(), rc_sampler_ptr(),
                   RunPlan(pure), Channel::Inclusive, LI6());
  const std::vector<double> ms = m_values(1.0);

  SUBCASE("(a) slot 1+k == slot 0 for a PURE category with the event's m") {
    for (std::size_t c = 0; c < s.n_cells(); c += 41) {
      for (std::size_t k = 0; k < pure.size(); ++k) {
        Event ev = inclusive_event(s, static_cast<int>(c), ms[k], 0.9);
        ev.spin_weights.assign(pure.size(), 1.0);
        mp.fill(ev);
        REQUIRE(ev.rc_weights.size() == kRcWeightCount * (1 + pure.size()));
        const std::size_t off = (1 + k) * kRcWeightCount;
        // SLOT-MAJOR: [0..2] is slot 0, [3*(1+k) ..] is category k.
        CHECK_CLOSE(ev.rc_weights[off + 0], ev.rc_weights[0], 1e-14);
        CHECK_CLOSE(ev.rc_weights[off + 1], ev.rc_weights[1], 1e-14);
        CHECK_CLOSE(ev.rc_weights[off + 2], ev.rc_weights[2], 1e-14);
      }
    }
  }

  SUBCASE("(b) the 1 : -2 : 1 ratio is on the NUMERATOR, not on tau") {
    // W_tensor is linear in P_zz; tau = W_tensor/W is not, because its
    // denominator moves with P_zz too.  So the design's "1 : 1 : -2 across
    // the tensor-thirds fills" is a statement about tau_k * W_k.
    for (std::size_t c = 0; c < s.n_cells(); c += 37) {
      const Event ev = inclusive_event(s, static_cast<int>(c), 1.0);
      double num[3], tau[3];
      for (std::size_t k = 0; k < 3; ++k) {
        tau[k] = mp.tensor_fraction(ev, k);
        num[k] = tau[k] / (1.0 - tau[k]);        // = W_tensor, since W = 1 + it
      }
      if (std::fabs(num[0]) < 1e-14) continue;
      CHECK_CLOSE(num[1] / num[0], -2.0, 1e-11);   // P_zz = -2 vs +1
      CHECK_CLOSE(num[2] / num[0], 1.0, 1e-11);    // P_zz = +1 vs +1
      // ... and tau's own ratio is NOT -2.
      CHECK(std::fabs(tau[1] / tau[0] + 2.0) > 1e-6);
    }
  }
}

// ------------------------------------------------------------------- T2

namespace {

std::shared_ptr<const TaggedModel> tagged_model() {
  static const std::shared_ptr<const TaggedModel> t =
      std::make_shared<const TaggedModel>(li6_alpha_channel());
  return t;
}

}  // namespace

TEST_CASE("T2: the TAGGED tau is the SAME closed form, not (1/2) A_zz^wf") {
  // StruckClusterOptions::inclusive_b1 is false by default, so the struck
  // cluster's DIS kernel carries NO b1 and the inclusive recipe would return
  // tau = 0 and price nothing.  The tagged tensor structure lives in the
  // M-dependence of n_M(k, c).
  const std::shared_ptr<const TaggedModel> tm = tagged_model();
  const RcModel m(RcMode::TensorBand, RcOptions(), rc_sampler_ptr(),
                  tensor_thirds_plan(0.0, 0.6), Channel::TaggedLi6Alpha, LI6(),
                  tm);
  CHECK(m.applies());
  CHECK_FALSE(m.tail_applies());          // the tag vetoes the elastic recoil

  // The (k, c) must be ON the model's grid nodes: `n_of_kc` is a NEAREST-cell
  // lookup, so an off-node point would compare two different cells.
  std::size_t checked = 0;
  for (std::size_t ik = 4; ik < tm->nk(); ik += 37) {
    for (std::size_t ic = 3; ic < tm->nc(); ic += 29) {
      const double k = tm->k()[ik];
      const double c = tm->c()[ic];
      const double n1 = tm->n_of_kc(1.0, k, c);
      const double n0 = tm->n_of_kc(0.0, k, c);
      const double nm = tm->n_of_kc(-1.0, k, c);
      if (!(n1 > 0.0) || !(n0 > 0.0)) continue;
      CHECK_CLOSE(nm, n1, 1e-12);         // n is EVEN in M (parity eigenstate)
      const double azz_wf = (n1 + nm - 2.0 * n0) / (n1 + nm + n0);
      for (double mm : {1.0, 0.0, -1.0}) {
        Event ev;
        ev.channel = Channel::TaggedLi6Alpha;
        ev.kin.k = k;
        ev.kin.cos_theta_k = c;
        ev.kin.x = 0.05;
        ev.spin.j = 1.0;
        ev.spin.m_ion = mm;
        const double p_zz = (mm == 0.0) ? -2.0 : 1.0;
        const double want = 0.5 * p_zz * azz_wf / (1.0 + 0.5 * p_zz * azz_wf);
        CHECK_CLOSE(m.tensor_fraction(ev), want, 1e-12);
        // The design's explicit forms, and NOT the first draft's (1/2)A_zz^wf.
        if (mm != 0.0) CHECK_CLOSE(m.tensor_fraction(ev),
                                   (n1 - n0) / (3.0 * n1), 1e-12);
        else CHECK_CLOSE(m.tensor_fraction(ev), 2.0 * (n0 - n1) / (3.0 * n0),
                         1e-12);
        if (std::fabs(azz_wf) > 0.05) {
          CHECK(std::fabs(m.tensor_fraction(ev) - 0.5 * azz_wf) > 1e-6);
        }
        ++checked;
      }
      // The tail is EXACTLY 1 on every tagged channel, and the run says why.
      // The published BAND EDGES are strictly positive on every one of them:
      // tau_tag itself is unbounded (n_M has nodes), and `band_tau_max` is
      // what keeps w_+- inside [1 - delta, 1 + delta].
      for (double mm : {1.0, 0.0, -1.0}) {
        Event ev;
        ev.channel = Channel::TaggedLi6Alpha;
        ev.kin.k = k; ev.kin.cos_theta_k = c; ev.kin.x = 0.05;
        ev.spin.j = 1.0; ev.spin.m_ion = mm;
        const RcWeights w = m.weights(ev);
        CHECK(w.tail == 1.0);
        REQUIRE(w.lo > 0.0);
        REQUIRE(w.hi > 0.0);
        CHECK(w.lo <= 2.0);
        CHECK(w.hi <= 2.0);
        CHECK_CLOSE(w.lo + w.hi, 2.0, 1e-15);
      }
    }
  }
  CHECK(checked > 20);
  // B4.  `rc_tail == 1` on a tagged channel is HALF a kinematic fact: the
  // ELASTIC recoil really is vetoed (x_L = 1, inside the beam envelope), the
  // QUASI-ELASTIC one is not (the A-1 remnant breaks up and its alpha lands at
  // x_L ~ 2/3, the tag window).  The run must say BOTH, so that a reader of
  // the npz cannot take the weight for a veto on the whole tail.
  const std::string why = m.exclusion_reason();
  CHECK(why.find("x_L") != std::string::npos);
  CHECK(why.find("FACT") != std::string::npos);
  CHECK(why.find("OMISSION") != std::string::npos);
  CHECK(why.find("QUASI-elastic") != std::string::npos);
  CHECK(why.find("2/3") != std::string::npos);
  CHECK(why.find("SPECTATOR") != std::string::npos);
  // The OLD text called the quasi-elastic half an unquantified "ASSUMPTION".
  // It is not: it is a background whose acceptance is now argued, and the
  // word must not come back without the argument.
  CHECK(why.find("ASSUMPTION") == std::string::npos);

  // Without a TaggedModel the tagged constructor REFUSES, rather than
  // silently returning tau = 0.
  CHECK_THROWS_AS(RcModel(RcMode::TensorBand, RcOptions(), rc_sampler_ptr(),
                          tensor_thirds_plan(0.0, 0.6),
                          Channel::TaggedLi6Alpha, LI6()),
                  std::runtime_error);
}

// ------------------------------------------------------------------- T5

TEST_CASE("the BAND is bounded: tau_tag diverges at the n_M nodes, the clamp holds") {
  // The published band edges must be non-negative.  On a TAGGED channel
  // tau_tag = 1 - nbar(k,c)/n_M(k,c) is unbounded -- wherever the event's own
  // n_M is near a node of the M-dependent spectator density (the 6Li M = 0
  // density has them) the ratio blows up -- and before `band_tau_max` the
  // shipped build produced rc_tensor_hi down to -1.79 on 6Li tagged-alpha and
  // -8.35 on 7Li, i.e. NEGATIVE weights in the npz.
  const std::shared_ptr<const TaggedModel> tm = tagged_model();
  const RcModel m(RcMode::TensorBand, RcOptions(), rc_sampler_ptr(),
                  tensor_thirds_plan(0.0, 0.6), Channel::TaggedLi6Alpha, LI6(),
                  tm);
  CHECK(m.options().band_tau_max == RC_BAND_TAU_MAX);

  bool saw_unbounded = false, saw_clip = false;
  double worst_lo = 1.0, worst_hi = 1.0, worst_tau = 0.0;
  for (std::size_t ik = 0; ik < tm->nk(); ik += 3) {
    for (std::size_t ic = 0; ic < tm->nc(); ic += 5) {
      for (double mm : {1.0, 0.0, -1.0}) {
        Event ev;
        ev.channel = Channel::TaggedLi6Alpha;
        ev.kin.k = tm->k()[ik];
        ev.kin.cos_theta_k = tm->c()[ic];
        ev.kin.x = 0.01;                       // delta(0.01) = 0.30, the max
        ev.spin.j = 1.0;
        ev.spin.m_ion = mm;
        // `tensor_fraction` is the UNCLAMPED diagnostic -- T1/T2 gate the
        // closed form against it, so the clamp lives in `weights()` only.
        const double tau = m.tensor_fraction(ev);
        if (std::fabs(tau) > RC_BAND_TAU_MAX) saw_unbounded = true;
        worst_tau = std::max(worst_tau, std::fabs(tau));
        const RcWeights w = m.weights(ev);
        if (w.band_clipped) saw_clip = true;
        REQUIRE(w.lo > 0.0);
        REQUIRE(w.hi > 0.0);
        CHECK(w.lo <= 2.0);
        CHECK(w.hi <= 2.0);
        // The clamp is EXACTLY band_tau_max * delta wide at its edge.
        CHECK(std::fabs(w.hi - 1.0) <= m.delta(ev.kin.x) * RC_BAND_TAU_MAX +
                                           1e-15);
        worst_lo = std::min(worst_lo, std::min(w.lo, w.hi));
        worst_hi = std::max(worst_hi, std::max(w.lo, w.hi));
        CHECK((w.band_clipped == (std::fabs(tau) > RC_BAND_TAU_MAX)));
      }
    }
  }
  MESSAGE("tagged-6Li-alpha band: max |tau| = " << worst_tau
          << " (UNCLAMPED); with the clamp every published edge lies in ["
          << worst_lo << ", " << worst_hi << "] = [1 - delta, 1 + delta]");
  // The clamp is SATURATED, i.e. it is doing work here and not decorative.
  CHECK_CLOSE(worst_lo, 1.0 - RC_DELTA_LOW_X, 1e-12);
  CHECK_CLOSE(worst_hi, 1.0 + RC_DELTA_LOW_X, 1e-12);
  CHECK(saw_unbounded);      // the divergence is real, not hypothetical
  CHECK(saw_clip);

  SUBCASE("band_tau_max = 0 disables the band entirely rather than NaN-ing") {
    RcOptions o;
    o.band_tau_max = 0.0;                 // "not > 0" => no clamp at all
    const RcModel unbounded(RcMode::TensorBand, o, rc_sampler_ptr(),
                            tensor_thirds_plan(0.0, 0.6),
                            Channel::TaggedLi6Alpha, LI6(), tm);
    bool saw_negative = false;
    for (std::size_t ik = 0; ik < tm->nk(); ik += 3) {
      for (std::size_t ic = 0; ic < tm->nc(); ic += 5) {
        Event ev;
        ev.channel = Channel::TaggedLi6Alpha;
        ev.kin.k = tm->k()[ik];
        ev.kin.cos_theta_k = tm->c()[ic];
        ev.kin.x = 0.01;
        ev.spin.j = 1.0;
        ev.spin.m_ion = 0.0;
        const RcWeights w = unbounded.weights(ev);
        if (w.lo < 0.0 || w.hi < 0.0) saw_negative = true;
      }
    }
    // This is the behaviour the clamp exists to prevent, pinned so nobody
    // "simplifies" band_tau_max away without seeing it.
    CHECK(saw_negative);
  }
}

TEST_CASE("T5: the tail's y and t_min behaviour") {
  const RcModel& m = rc_model();
  const double s = rc_sampler().s();
  const double y_max = rc_sampler().scenario().y_max;

  SUBCASE("(i) the y-scan at Q2 = 5 GeV^2 tracks Y_+") {
    const double q2 = 5.0;
    const double ys[] = {0.1, 0.3, 0.5, 0.7, 0.9, 0.97};
    double scaled[6];
    for (int i = 0; i < 6; ++i) {
      const double y = ys[i];
      const double x = q2 / (y * s);
      // (iii) every grid point is inside the sampler's own y window.
      CHECK_MESSAGE(y <= y_max, "y = " << y << " exceeds y_max = " << y_max
                                       << " at s = " << s);
      const double yp = (1.0 + (1.0 - y) * (1.0 - y)) / (1.0 - y);
      scaled[i] = m.tail_ratio_at(x, q2, 0.0) / yp;
      CHECK(scaled[i] > 0.0);
    }
    for (int i = 1; i < 6; ++i) {
      CHECK_MESSAGE(scaled[i] > scaled[i - 1],
                    "(w_tail-1)/Y_+ fell between y = " << ys[i - 1] << " and "
                                                       << ys[i]);
    }
    // Over the LAST THREE points it varies by less than 3x.
    CHECK(scaled[5] / scaled[3] < 3.0);
    // DESIGN vs CODE: sec. 5 T5(i) asks for "< 3x across the scan".  A correct
    // implementation does NOT meet that, and cannot: at FIXED Q^2 the scan
    // also sweeps x by a factor 10 (x = Q^2/(y s) runs 0.0126 -> 0.0013), and
    // the design's own sec. 8.1 says "the x-dependence is the strong one".
    // Measured spread over the full scan: ~290x.  The claim is kept where it
    // is true -- the last three points -- and the monotonicity, which IS the
    // Y_+ statement, is asserted everywhere.
    CHECK(scaled[5] / scaled[0] > 10.0);
  }

  SUBCASE("(ii) at fixed y the tail FALLS with x as t_min climbs") {
    for (double y : {0.05, 0.2, 0.6}) {
      double prev = 1e300;
      for (double x : {0.005, 0.01, 0.03, 0.1, 0.3}) {
        const double q2 = x * y * s;
        const double r = m.tail_ratio_at(x, q2, 0.0);
        CHECK_MESSAGE(r < prev, "the tail rose with x at y = " << y
                                                              << ", x = " << x);
        prev = r;
        // t_min = (x m_p)^2/(1 - x_A) climbs into the 6Li form factor, which
        // is dead above t ~ 0.25 GeV^2 <=> x ~ 0.53.
        const double x_a = x * PROTON_MASS / LI6().mass();
        const double t_min =
            4.0 * LI6().mass() * LI6().mass() * polrad_eta_min_ur(x_a);
        CHECK((t_min < 0.26 || x > 0.5));
      }
    }
  }

  SUBCASE("the table refuses to extrapolate") {
    CHECK_THROWS_AS(m.tail_ratio_at(1e-9, 5.0, 0.0), std::runtime_error);
    CHECK_THROWS_AS(m.tail_ratio_at(0.1, 1e-6, 0.0), std::runtime_error);
  }
}

// -------------------------------------------------------------- T12, T12'

namespace {

RcModel band_model(double fq_scale, double fm_scale = 1.0) {
  RcOptions o;
  o.fq_scale = fq_scale;
  o.tail_tensor_scale = fm_scale;
  return RcModel(RcMode::TensorBand, o, rc_sampler_ptr(),
                 tensor_thirds_plan(0.0, 0.6), Channel::Inclusive, LI6());
}

/// The Q_N-DEPENDENT part of the tail ratio at a TABLE NODE, where the
/// interpolation is exact and the whole thing is (q_n/6) sigma^el_T / Born.
///
/// `q_n` is deliberately UNPHYSICALLY LARGE at the call sites below (a real
/// spin-1 fill has |Q_N| <= 2).  The function is EXACTLY linear in Q_N -- that
/// is the whole point of keeping (1 + w_avg) out of `tail_ratio_at` -- and the
/// tensor piece is O(1e-3) of the unpolarised one, so at |Q_N| ~ 1 the
/// difference below is a cancellation that costs 13 of the 16 digits and no
/// algebraic identity could be asserted at 1e-12 through it.  At Q_N ~ 600 the
/// two terms are comparable and the identity is asserted where it lives.
double tensor_part(const RcModel& m, double x, double q2, double q_n) {
  return m.tail_ratio_at(x, q2, q_n) - m.tail_ratio_at(x, q2, 0.0);
}

}  // namespace

TEST_CASE("T12: the quadrupole band is QUADRATIC in fq_scale, not linear") {
  // Eq. (38)'s sigma_q^d carries F_q(3F_c + 3 eta F_m + eta F_q) and
  // F_q(4F_c - 3x F_m + (4/3) eta F_q), so sigma^el_T is a QUADRATIC in
  // fq_scale.  "Exactly linear to 1e-12" -- the first draft's assertion --
  // fails against a correct implementation, which is why the band must be RUN
  // at 0/1/2 and never rescaled from one run.
  const RcModel m0 = band_model(0.0);
  const RcModel m1 = band_model(1.0);
  const RcModel m2 = band_model(2.0);
  const RcModel mh = band_model(0.5);
  const RcModel m3 = band_model(3.0);
  const double q_n = 600.0;   // see `tensor_part`: linear in Q_N by design
  // Evaluate AT A NODE, where the bilinear reconstruction is exact and the
  // Q_N-dependent part is (q_n/6) sigma^el_T / Born to the last bit.
  const std::vector<double>& xs = m1.table_x();
  const std::vector<double>& ys = m1.table_y();
  const double s = rc_sampler().s();
  int quad_seen = 0, points = 0;
  for (std::size_t i = 8; i < xs.size(); i += 9) {
    for (std::size_t j = 4; j < ys.size(); j += 7) {
      const double x = xs[i], y = ys[j], q2 = x * y * s;
      if (!(q2 > 0.2)) continue;
      const double t0 = tensor_part(m0, x, q2, q_n);
      const double t1 = tensor_part(m1, x, q2, q_n);
      const double t2 = tensor_part(m2, x, q2, q_n);
      // T(s) = T0 + s T1 + s^2 T2, exactly determined by three points.
      const double c2 = 0.5 * (t2 - 2.0 * t1 + t0);
      const double c1 = t1 - t0 - c2;
      CHECK_CLOSE(tensor_part(mh, x, q2, q_n),
                  t0 + 0.5 * c1 + 0.25 * c2, 1e-11);
      CHECK_CLOSE(tensor_part(m3, x, q2, q_n), t0 + 3.0 * c1 + 9.0 * c2, 1e-11);
      if (std::fabs(c2) > 1e-9 * std::fabs(c1)) ++quad_seen;
      ++points;
    }
  }
  CHECK(points > 10);
  CHECK_MESSAGE(quad_seen > 0,
                "the quadratic coefficient vanished everywhere -- a LINEAR "
                "implementation would pass a linearity test and fail the "
                "physics");
}

TEST_CASE("T12': tail_tensor_scale moves what fq_scale does NOT span") {
  // F_m carries the eta F_m^2 tensor term and two interferences (3 eta F_m,
  // -3 x F_m), so an fq_scale-only band leaves it unpriced.
  const double s = rc_sampler().s();
  const double x = 0.1, q2 = 5.0;
  REQUIRE(q2 / (x * s) < 0.985);
  const RcModel a = band_model(1.0, 0.5);
  const RcModel b = band_model(1.0, 1.0);
  const RcModel c = band_model(1.0, 2.0);
  const double ta = tensor_part(a, x, q2, 600.0);
  const double tb = tensor_part(b, x, q2, 600.0);
  const double tc = tensor_part(c, x, q2, 600.0);
  CHECK(std::fabs(ta - tb) > 1e-6 * std::fabs(tb));
  CHECK(std::fabs(tc - tb) > 1e-6 * std::fabs(tb));
  // ... and it is NOT a rescaling of the fq band either.
  CHECK(std::fabs((tc - tb) / (tb - ta)) > 1e-6);
}

// ---------------------------------------------------------------- T13, T14

TEST_CASE("T13: RcMode::Off is identically 1, and fills NOTHING") {
  // The mode is a CONSTRUCTOR ARGUMENT, not an RcOptions field: there is one
  // source of truth (`PipelineConfig::rc`) and the two cannot disagree.
  const RcModel m(RcMode::Off, RcOptions(), rc_sampler_ptr(),
                  tensor_thirds_plan(0.0, 0.6), Channel::Inclusive, LI6());
  CHECK(m.mode() == RcMode::Off);
  CHECK_FALSE(m.applies());
  CHECK_FALSE(m.tail_applies());
  CHECK(m.exclusion_reason().find("off") != std::string::npos);
  Event ev = inclusive_event(rc_sampler(), 3, 1.0, 1.1);
  ev.spin_weights.assign(3, 1.0);
  const RcWeights w = m.weights(ev);
  CHECK(w.lo == 1.0);
  CHECK(w.hi == 1.0);
  CHECK(w.tail == 1.0);
  m.fill(ev);
  CHECK(ev.rc_weights.empty());          // and so the npz / HepMC are today's
  // No form factor is built at all, and asking for one says why.
  CHECK_THROWS_AS(m.ff(), std::runtime_error);
  CHECK(m.ff_provenance().find("does not apply") != std::string::npos);
}

TEST_CASE("T14: the per-channel table, TOTAL over Channel") {
  const RunPlan plan = tensor_thirds_plan(0.0, 0.6);
  const Channel all[] = {Channel::Inclusive,      Channel::TaggedLi6Alpha,
                         Channel::TaggedLi6D,     Channel::TaggedLi7Alpha,
                         Channel::TaggedLi7T,     Channel::TaggedDeuteronP,
                         Channel::TaggedDeuteronN, Channel::TaggedHe3P,
                         Channel::CoherentLi6};
  int n_band = 0, n_tail = 0, n_none = 0;
  for (Channel ch : all) {
    const bool tagged = (ch != Channel::Inclusive && ch != Channel::CoherentLi6);
    const RcModel m(RcMode::TensorBand, RcOptions(), rc_sampler_ptr(), plan, ch,
                    LI6(), tagged ? tagged_model() : nullptr);
    if (ch == Channel::CoherentLi6) {
      CHECK_FALSE(m.applies());
      CHECK_FALSE(m.tail_applies());
      CHECK(m.exclusion_reason().find("AZIMUTHAL") != std::string::npos);
      Event ev = inclusive_event(rc_sampler(), 5, 1.0);
      ev.channel = ch;
      const RcWeights w = m.weights(ev);
      CHECK(w.lo == 1.0);
      CHECK(w.hi == 1.0);
      CHECK(w.tail == 1.0);
      ++n_none;
    } else if (tagged) {
      CHECK(m.applies());
      CHECK_FALSE(m.tail_applies());
      CHECK_FALSE(m.exclusion_reason().empty());
      ++n_band;
    } else {
      CHECK(m.applies());
      CHECK(m.tail_applies());
      CHECK(m.exclusion_reason().empty());
      ++n_tail;
    }
  }
  CHECK(n_band == 7);
  CHECK(n_tail == 1);
  CHECK(n_none == 1);

  SUBCASE("the tail is emitted at EVERY theta_S -- POLRAD Eq. (43)") {
    // sigma_q,perp = -(1/2) sigma_q,par is the theta_S = 90 deg instance of
    // Q_N -> P_zz^eff = 3 Q_NN P_2(cos theta_S).
    const RcModel& m = rc_model();
    CHECK(m.tail_applies());
    const double s = rc_sampler().s();
    for (double x : {0.01, 0.1}) {
      for (double q2 : {3.0, 8.0}) {
        // Q_N = 600 rather than 1: the ratio is EXACTLY linear in Q_N (that
        // is why `tail_ratio_at` keeps (1 + w_avg) out), and at |Q_N| ~ 1 the
        // tensor piece is 1e-3 of the unpolarised one, so the difference below
        // would be a 13-digit cancellation.  Both scalings are asserted.
        const double q_par = 600.0;               // theta_S = 0,  P_2 = 1
        const double q_perp = -0.5 * q_par;       // theta_S = 90, P_2 = -1/2
        const double base = m.tail_ratio_at(x, q2, 0.0);
        const double par = m.tail_ratio_at(x, q2, q_par) - base;
        const double perp = m.tail_ratio_at(x, q2, q_perp) - base;
        REQUIRE(std::fabs(par) > 0.0);
        CHECK_CLOSE(perp / par, -0.5, 1e-12);
        // ... and at the PHYSICAL Q_N = +1 / -1/2, to what the cancellation
        // leaves.
        CHECK_CLOSE((m.tail_ratio_at(x, q2, -0.5) - base) /
                        (m.tail_ratio_at(x, q2, 1.0) - base),
                    -0.5, 1e-8);
        // NOTE the design's T14 says "rc_tail - 1 at theta_S = 90 deg is -1/2
        // its theta_S = 0 value".  That is true of the TENSOR part only: the
        // unpolarised sigma^el_U and sigma^q_U carry no Q_N at all and do not
        // halve.  Asserting it on the whole of w_tail - 1 would be asserting
        // that the unpolarised tail vanishes.
        const double whole = m.tail_ratio_at(x, q2, -0.5) /
                             m.tail_ratio_at(x, q2, 1.0);
        CHECK(std::fabs(whole + 0.5) > 0.1);
      }
    }
    // ... and it really is emitted at a transverse axis, through a real event.
    const std::vector<SpinCategory> perp_cats = {
        SpinCategory("t+perp", 1.0, {1.0, 0.0, 0.0}, 0, 0.0, kPi / 2.0, 0.0)};
    const RcModel mp(RcMode::TensorBand, RcOptions(), rc_sampler_ptr(),
                     RunPlan(perp_cats), Channel::Inclusive, LI6());
    CHECK(mp.tail_applies());
    Event ev = inclusive_event(rc_sampler(), 9, 1.0, 0.0, kPi / 2.0);
    CHECK(mp.weights(ev).tail != 1.0);
  }

  SUBCASE("a polarised beam is REFUSED, not silently priced as rank 2") {
    const std::vector<SpinCategory> hel = {
        SpinCategory("h+", 1.0, {1.0, 0.0, 0.0}, +1, 0.7)};
    CHECK_THROWS_AS(RcModel(RcMode::TensorBand, RcOptions(), rc_sampler_ptr(),
                            RunPlan(hel), Channel::Inclusive, LI6()),
                    std::runtime_error);
  }

  SUBCASE("PolradFull BUILDS now, and its two preconditions are refused") {
    // This subcase asserted the OPPOSITE until 2026-09-06 ("PolradFull is
    // refused, with the design section that owns it").  Eq. (18) + Appendix B
    // + Eq. (A.4) is implemented; what stays refused is a resolution that
    // cannot resolve the s-/p-peaks it exists to carry.
    RcOptions o;
    o.tail_model = RcTailModel::PolradFull;
    CHECK_NOTHROW(RcModel(RcMode::TensorBand, o, rc_sampler_ptr(), plan,
                          Channel::Inclusive, LI6()));
    RcOptions coarse = o;
    coarse.n_eta = 32;
    CHECK_THROWS_AS(RcModel(RcMode::TensorBand, coarse, rc_sampler_ptr(), plan,
                            Channel::Inclusive, LI6()),
                    std::runtime_error);
    // ... and the s/p tensor STAND-IN is refused on it, for the OPPOSITE
    // reason it is refused on TPeak: PolradFull computes the s-/p-peaks'
    // tensor content, so the stand-in would double-count a term that ran.
    RcOptions sp = o;
    sp.sp_tensor_scale = 1.0;
    CHECK_THROWS_AS(RcModel(RcMode::TensorBand, sp, rc_sampler_ptr(), plan,
                            Channel::Inclusive, LI6()),
                    std::runtime_error);
    // The `long double` precondition is a PLATFORM fact, so it is asserted
    // rather than provoked: on a platform where it fails the constructor
    // refuses instead of shipping a two-figure tensor tail.
    CHECK(std::numeric_limits<long double>::digits >
          std::numeric_limits<double>::digits);
  }
}

// ------------------------------------------------------------------- T18

TEST_CASE("T18: RcModel is immutable and thread-safe after construction") {
  const RcModel& m = rc_model();
  const InclusiveSampler& s = rc_sampler();
  std::vector<Event> evs;
  for (std::size_t c = 0; c < s.n_cells(); c += 5) {
    for (double mm : {1.0, 0.0, -1.0}) {
      Event ev = inclusive_event(s, static_cast<int>(c), mm, 0.3 * c);
      ev.spin_weights.assign(3, 1.0);
      evs.push_back(ev);
    }
  }
  REQUIRE(evs.size() > 100);
  std::vector<std::vector<double>> serial(evs.size());
  for (std::size_t i = 0; i < evs.size(); ++i) {
    Event ev = evs[i];
    m.fill(ev);
    serial[i] = ev.rc_weights;
  }
  std::vector<std::vector<double>> par(evs.size());
  const unsigned nthreads = 4;
  std::vector<std::thread> pool;
  for (unsigned t = 0; t < nthreads; ++t) {
    pool.emplace_back([&, t] {
      for (std::size_t i = t; i < evs.size(); i += nthreads) {
        Event ev = evs[i];
        m.fill(ev);
        par[i] = ev.rc_weights;
      }
    });
  }
  for (std::thread& th : pool) th.join();
  for (std::size_t i = 0; i < evs.size(); ++i) {
    REQUIRE(par[i].size() == serial[i].size());
    for (std::size_t j = 0; j < par[i].size(); ++j) {
      CHECK(par[i][j] == serial[i][j]);   // BIT for bit
    }
  }
}

TEST_CASE("T8(d): RcTailModel::TPeakPlusLL -- what it adds, and what it does NOT") {
  // The promotion of the s-/p-peaks from a test-local construction to a
  // shipped, OPT-IN tail model.  Four things are gated: that the default did
  // not move, that the new columns are what `ll_peaks_*` returns, that they
  // enter the UNPOLARISED numerator only, and that the tensor FRACTION falls
  // as a consequence -- which is a physics change and is asserted, not hidden.
  RcOptions plus;
  plus.tail_model = RcTailModel::TPeakPlusLL;
  const RcModel mp(RcMode::TensorBand, plus, rc_sampler_ptr(),
                   tensor_thirds_plan(0.0, 0.6), Channel::Inclusive, LI6());
  const RcModel& mt = rc_model();          // the shipped default, TPeak
  const double s = rc_sampler().s();

  SUBCASE("(a) the DEFAULT is untouched: the s+p columns are identically 0") {
    const std::vector<double>& usp = mt.sigma_tail_u_sp();
    const std::vector<double>& qsp = mt.sigma_tail_qe_sp();
    REQUIRE(usp.size() == mt.sigma_tail_u().size());
    REQUIRE(qsp.size() == mt.sigma_tail_qe().size());
    for (std::size_t k = 0; k < usp.size(); ++k) {
      REQUIRE(usp[k] == 0.0);
      REQUIRE(qsp[k] == 0.0);
    }
    // ... and the three OLD tables are bit for bit the same on both models,
    // because `tail_sigma_at` computes them from the same quadratures.  This
    // is the structural half of the "--rc off / TPeak stays byte-identical"
    // promise; the whole-file half is the npz diff in phase_B_numbers.md B2.
    REQUIRE(mp.sigma_tail_u().size() == mt.sigma_tail_u().size());
    for (std::size_t k = 0; k < mt.sigma_tail_u().size(); ++k) {
      REQUIRE(mp.sigma_tail_u()[k] == mt.sigma_tail_u()[k]);
      REQUIRE(mp.sigma_tail_t()[k] == mt.sigma_tail_t()[k]);
      REQUIRE(mp.sigma_tail_qe()[k] == mt.sigma_tail_qe()[k]);
    }
  }

  SUBCASE("(b) the new columns ARE `ll_peaks_*`, with the right 1/A powers") {
    // The one way this block can be wrong by a factor A and still look
    // plausible: `ll_peaks_spin1` is WHOLE-NUCLEUS d^2 sigma/(dx_A dy) and
    // takes 1/A^2 (per nucleon x the Jacobian dx_A/dx), while `ll_peaks_qe`
    // is already in nucleon invariants and takes 1/A -- exactly as
    // `polrad_sigma_el_u` and `polrad_sigma_qe_u` do.
    const std::shared_ptr<HoSpin1FF> ff = HoSpin1FF::for_ion(LI6());
    const double m6 = LI6().mass(), s6 = 6.0 * s;
    for (double x : {0.01, 0.1}) {
      for (double q2 : {3.0, 8.0}) {
        const double y = q2 / (x * s);
        const RcModel::TailTriple tt = mp.tail_sigma_at(x, q2);
        const LlPeaks el = ll_peaks_spin1(*ff, x / 6.0, y, s6, m6);
        const LlPeaks qe = ll_peaks_qe(3, 3, x, y, s, RcOptions().qe_kf_gev);
        CHECK_CLOSE(tt.u_sp, (el.s + el.p) / 36.0, 1e-14);
        CHECK_CLOSE(tt.qe_sp, (qe.s + qe.p) / 6.0, 1e-14);
        // and TPeak fills neither
        const RcModel::TailTriple t0 = mt.tail_sigma_at(x, q2);
        CHECK(t0.u_sp == 0.0);
        CHECK(t0.qe_sp == 0.0);
      }
    }
  }

  SUBCASE("(c) the s+p enter the UNPOLARISED numerator ONLY") {
    // THE PHYSICS TRAP.  The tensor term is (q_n/6) (sigma_t/sigma_u) sigma_u
    // = (q_n/6) sigma_t, the t-peak's OWN tensor quadrature.  If the s+p had
    // been folded into `su` before that product, the tensor piece would have
    // grown with them -- silently asserting that the s-/p-peaks carry the
    // t-peak's tensor-to-unpolarised ratio, which POLRAD does not supply.
    // Gate: the Q_N-DERIVATIVE of the ratio is IDENTICAL on the two models.
    //
    // THE TOLERANCE IS 1e-9 AND IT IS NOT SLACK.  The quantity is a
    // DIFFERENCE of two ratios, and adding the s+p makes the ratio itself up
    // to ~4000x larger at x = 0.1, Q^2 = 3 while the tensor term does not
    // move -- so the subtraction loses log10(base/|d|) ~ 5 digits there and
    // the double-precision floor is eps * base/|d| ~ 2e-16 * 1.6e5 ~ 3e-11.
    // Measured worst case over these four points: 1.5e-11 (the rest are
    // 1e-15 or exactly 0).  A 1e-12 gate here would be gating round-off, not
    // the tensor numerator.
    for (double x : {0.01, 0.1}) {
      for (double q2 : {3.0, 8.0}) {
        const double qn = 600.0;                 // linear in q_n; see T14
        const double dt = mt.tail_ratio_at(x, q2, qn) -
                          mt.tail_ratio_at(x, q2, 0.0);
        const double dp = mp.tail_ratio_at(x, q2, qn) -
                          mp.tail_ratio_at(x, q2, 0.0);
        REQUIRE(std::fabs(dt) > 0.0);
        CHECK_CLOSE(dp, dt, 1e-9);
      }
    }
  }

  SUBCASE("(d) ... so the TENSOR FRACTION of the tail FALLS.  Stated, not hidden") {
    // The documented consequence of (c): the denominator grows and the tensor
    // numerator does not, so (w_tail(q_n) - w_tail(0))/(w_tail(0) - 1) drops.
    // The tensor part of the s-/p-peaks is BOUNDED (RcOptions::sp_tensor_scale,
    // 2026-09-06), not computed, and this is what makes TPeakPlusLL a
    // systematic to run BESIDE TPeak and never instead.  That bound is EMPTY
    // at these very (x, Q^2): see the B1 case's own "THE PRICE" subcase.
    bool saw_a_drop = false;
    for (double x : {0.01, 0.1}) {
      for (double q2 : {3.0, 8.0}) {
        const double qn = 600.0;
        const double ft = (mt.tail_ratio_at(x, q2, qn) -
                           mt.tail_ratio_at(x, q2, 0.0)) /
                          mt.tail_ratio_at(x, q2, 0.0);
        const double fp = (mp.tail_ratio_at(x, q2, qn) -
                           mp.tail_ratio_at(x, q2, 0.0)) /
                          mp.tail_ratio_at(x, q2, 0.0);
        CHECK(std::fabs(fp) <= std::fabs(ft) * (1.0 + 1e-12));
        if (std::fabs(fp) < std::fabs(ft) * 0.999) saw_a_drop = true;
        MESSAGE("x = " << x << ", Q2 = " << q2 << ": tensor fraction of the "
                "tail " << ft << " (t-peak) -> " << fp << " (t-peak+ll)");
      }
    }
    CHECK(saw_a_drop);
  }

  SUBCASE("(e) the tail GROWS, and by how much -- the acceptance numbers") {
    // Never smaller: the s+p are positive everywhere the radiator is defined.
    struct Row { double x, q2; };
    for (const Row& r : {Row{0.01, 3.0}, Row{0.01, 8.0}, Row{0.03, 8.0},
                         Row{0.1, 8.0}, Row{0.1, 40.0}, Row{0.3, 40.0}}) {
      const double a = mt.tail_ratio_at(r.x, r.q2, 0.0);
      const double b = mp.tail_ratio_at(r.x, r.q2, 0.0);
      CHECK(b >= a * (1.0 - 1e-12));
      MESSAGE("x = " << r.x << ", Q2 = " << r.q2 << ": tail ratio " << a
                     << " -> " << b << "  (x" << (a > 0.0 ? b / a : 0.0)
                     << ")");
    }
    // The SURVIVING per-cell agreement window, and nothing wider.  These three
    // rows sit at y = 0.5, i.e. inside Q^2 >= 20 GeV^2 AND 0.15 <= y <= 0.7 --
    // the ONLY window in which the two models agree cell by cell (660 accepted
    // cells, worst 0.55 %).  Q^2 >= 20 alone does NOT buy 1 %: over the CLI
    // grid's own cells 346 of 1399 disagree by more than that, 318 of them at
    // y < 0.1.  NARROWED 2026-09-04 -- this comment used to read "Q^2 >= 20
    // GeV^2 AND y < 0.9 ... agree to better than 1 %", which is three y = 0.5
    // points read as a window.  The census is gated by (i) and the y -> 1 edge
    // by (f); the three CHECKs below are unchanged.  The 346/1399 census is a
    // measurement record on the CLI grid, NOT something (i) gates -- (i) runs
    // the 40 x 24 test sampler and reports 45 of 188 there.
    for (double x : {0.03, 0.1, 0.3}) {
      const double y = 0.5, q2 = x * y * s;
      REQUIRE(q2 > 20.0);
      const double a = mt.tail_ratio_at(x, q2, 0.0);
      const double b = mp.tail_ratio_at(x, q2, 0.0);
      CHECK(b / a < 1.01);
    }
  }

  SUBCASE("(f) the y -> 1 edge, where the s-peak BEATS Y_+ -- the honest half") {
    // At y -> 1 the t-peak grows like Y_+ = [1+(1-y)^2]/(1-y) ~ 1/(1-y).  The
    // s-peak grows FASTER: z_s = (1-y)/(1-x_A y) -> 0 puts the elastic vertex
    // at Q'^2 = z_s Q^2 -> 0, where the form factor is 1 and the Born-like
    // 1/Q'^4 is enormous.  So "the t-peak is > 99 % at Q^2 >= 20 GeV^2" is a
    // statement about the BULK and NOT about the y -> 1 rows, and this
    // subcase exists so a later reader cannot mistake one for the other.
    const double x = 0.01, y = 0.985, q2 = x * y * s;
    REQUIRE(q2 > 20.0);
    const double a = mt.tail_ratio_at(x, q2, 0.0);
    const double b = mp.tail_ratio_at(x, q2, 0.0);
    MESSAGE("x = 0.01, y = 0.985, Q2 = " << q2 << ": tail ratio " << a
                                         << " -> " << b);
    CHECK(b / a > 1.3);
  }

  SUBCASE("(g) the node clip is REPORTED on the new model, not hidden") {
    // `build_tail_tables` counts a node clipped on the SAME sum the event
    // sees, s+p included, so switching model moves the statistic.  It is zero
    // on the default and NOT zero on TPeakPlusLL, concentrated in y > 0.9 --
    // which is exactly (f)'s edge.
    CHECK(mt.clipped_cell_fraction() == 0.0);
    CHECK(mp.clipped_cell_fraction() > 0.0);
    MESSAGE("clipped_cell_fraction: " << mt.clipped_cell_fraction() << " -> "
                                      << mp.clipped_cell_fraction());
    const std::array<double, 3> bt = mt.clipped_fraction_by_y();
    const std::array<double, 3> bp = mp.clipped_fraction_by_y();
    for (int i = 0; i < 3; ++i) CHECK(bp[i] >= bt[i]);
    CHECK(bp[2] > bp[0]);            // the y > 1 edge, not the bulk
    MESSAGE("clipped_fraction_by_y (t-peak+ll): " << bp[0] << " / " << bp[1]
                                                  << " / " << bp[2]);
  }

  SUBCASE("(h) PolradFull BUILDS since 2026-09-06; the name table is total") {
    RcOptions full;
    full.tail_model = RcTailModel::PolradFull;
    CHECK_NOTHROW(RcModel(RcMode::TensorBand, full, rc_sampler_ptr(),
                          tensor_thirds_plan(0.0, 0.6), Channel::Inclusive,
                          LI6()));
    CHECK(std::string(rc_tail_model_name(RcTailModel::TPeak)) == "t-peak");
    CHECK(std::string(rc_tail_model_name(RcTailModel::PolradFull)) ==
          "polrad-full");
    CHECK(std::string(rc_tail_model_name(RcTailModel::TPeakPlusLL)) ==
          "t-peak+ll");
  }

  SUBCASE("(i) the LOW-y corner: where TPeakPlusLL BREAKS DOWN, and the gate "
          "that was missing") {
    // WHY THIS SUBCASE EXISTS.  Until 2026-09-04 four shipped documents said
    // the two models "agree to 0.16 % at Q^2 >= 20 GeV^2 and y <= 0.7", and
    // (e) above gated it at y = 0.5 only, (c)(i) at y in {0.5, 0.9} and (f) at
    // y = 0.985.  NOTHING probed y -> 0, and the claim is FALSE there: over
    // the sampler's own accepted cells with Q^2 >= 20 and y <= 0.9, a quarter
    // of them disagree by more than 1 % and the worst by a factor of
    // thousands.  The claim survives only as an EVENT-WEIGHTED statement.
    //
    // THE MECHANISM, and it is NOT (f)'s.  At y -> 1 the s-peak wins because
    // z_s -> 0 drags the elastic vertex to Q'^2 -> 0.  At y -> 0 the opposite
    // happens: z_s -> 1, Q'^2 -> Q^2, 6Li's coherent form factor is DEAD
    // there, and what wins is the UNCANCELLED SOFT RADIATOR --
    //   1 - z_s = y(1 - x_A)/(1 - x_A y),   1 - z_p = y(1 - x_A)
    // so D(z) = (alpha/pi) L (1+z^2)/(1-z) ~ (2 alpha/pi) L/[y(1 - x_A)]
    // diverges like 1/y.  The t-peak gets no matching growth (Y_+ -> 2) and is
    // simultaneously CRUSHED, because Q^2 >= 20 forces x y >= 5.0e-3 -- low y
    // means HIGH x -- and its own elastic vertex starts at
    // t_min = M_A^2 x_A^2/(1 - x_A), which grows like x^2.
    const double x_bad = 0.794328, y_bad = 0.00879665;
    const double q2_bad = x_bad * y_bad * s;
    REQUIRE(q2_bad > 20.0);
    REQUIRE(y_bad < 0.9);

    // (1) THE RAW QUADRATURE at that point -- no table, no interpolation, so
    //     this number does not move with the sampler grid.
    {
      const RcModel::TailTriple a = mt.tail_sigma_at(x_bad, q2_bad);
      const RcModel::TailTriple b = mp.tail_sigma_at(x_bad, q2_bad);
      const double t_only = a.u + a.qe;
      const double t_plus = b.u + b.qe + b.u_sp + b.qe_sp;
      REQUIRE(t_only > 0.0);
      MESSAGE("x = " << x_bad << ", y = " << y_bad << ", Q2 = " << q2_bad
                     << ": t-peak " << t_only << " -> t+s+p " << t_plus
                     << "  (x" << t_plus / t_only << ")");
      CHECK(t_plus / t_only > 3000.0);
      // ... and the excess is ENTIRELY QUASI-ELASTIC.  6Li's coherent form
      // factor is dead at Q'^2 ~ Q^2 ~ 28 GeV^2, so the ELASTIC s+p is not
      // small there -- it is EXACTLY ZERO.  This is the half of the mechanism
      // that (f)'s "the form factor is 1 at Q'^2 -> 0" reading gets backwards.
      CHECK(b.u_sp == 0.0);
      CHECK(b.qe_sp > 0.0);
    }

    // (2) THE CENSUS over the sampler's own accepted cells.  This is the
    //     statement the documents now carry, and it is measured here rather
    //     than sampled at four hand-picked (x, y).
    {
      std::size_t n_win = 0, n_bad = 0, n_bad_low_y = 0;
      double worst = 1.0, worst_y = -1.0, worst_x = -1.0;
      for (std::size_t c = 0; c < rc_sampler().n_cells(); ++c) {
        const double x = rc_sampler().x_cells()[c];
        const double q2 = rc_sampler().q2_cells()[c];
        const double y = q2 / (x * s);
        if (!(q2 >= 20.0) || !(y <= 0.9)) continue;
        ++n_win;
        const double a = mt.tail_ratio(static_cast<int>(c), 0.0);
        const double b = mp.tail_ratio(static_cast<int>(c), 0.0);
        if (!(a > 0.0)) continue;
        const double r = b / a;
        if (r > 1.01) { ++n_bad; if (y < 0.1) ++n_bad_low_y; }
        if (r > worst) { worst = r; worst_x = x; worst_y = y; }
      }
      REQUIRE(n_win > 0);
      MESSAGE("Q^2 >= 20, y <= 0.9: " << n_bad << " of " << n_win
                                      << " accepted cells differ by > 1 % ("
                                      << (100.0 * static_cast<double>(n_bad) /
                                          static_cast<double>(n_win))
                                      << " %), " << n_bad_low_y
                                      << " of them at y < 0.1; worst x"
                                      << worst << " at x = " << worst_x
                                      << ", y = " << worst_y);
      // The failure EXISTS, is not a rare outlier, is ORDERS of magnitude and
      // sits at y -> 0.  A future change that quietly restores the old
      // "agrees to 0.16 %" reading has to break one of these four.
      CHECK(n_bad > 0);
      CHECK(static_cast<double>(n_bad) > 0.05 * static_cast<double>(n_win));
      CHECK(worst > 100.0);
      CHECK(worst_y < 0.1);
      CHECK(n_bad_low_y > n_bad / 2);

      // ... and the window that DOES hold per cell, which is the only
      // per-cell agreement statement any document may quote.
      for (std::size_t c = 0; c < rc_sampler().n_cells(); ++c) {
        const double x = rc_sampler().x_cells()[c];
        const double q2 = rc_sampler().q2_cells()[c];
        const double y = q2 / (x * s);
        if (!(q2 >= 20.0) || !(y >= 0.15) || !(y > 0.0) || !(y <= 0.7)) {
          continue;
        }
        const double a = mt.tail_ratio(static_cast<int>(c), 0.0);
        const double b = mp.tail_ratio(static_cast<int>(c), 0.0);
        REQUIRE(a > 0.0);
        CHECK_MESSAGE(b / a < 1.01,
                      "0.15 <= y <= 0.7 is the surviving per-cell window and "
                      "it just failed at x = " << x << ", y = " << y
                          << ": ratio " << b / a);
      }
    }

    // (3) THE MECHANISM, gated: along a FIXED Q^2 the disagreement leaves 1 %
    //     exactly where the leading-log radiator D(z_s) passes 1 -- i.e. where
    //     one emission stops being the right expansion.  `TPeakPlusLL` is
    //     BREAKING DOWN there; the factor 6444 is not evidence that `TPeak` is
    //     low by 6444.
    {
      const double xy = 0.006;             // Q^2 = 23.88 GeV^2 on this beam
      const double q2 = xy * s;
      REQUIRE(q2 > 20.0);
      double last_d = 0.0, last_r = 0.0;
      for (const double x : {0.01, 0.03, 0.05, 0.10, 0.20, 0.30, 0.50, 0.72}) {
        const double y = xy / x, x_a = x / 6.0;
        const double z_s = (1.0 - y) / (1.0 - x_a * y);
        // the soft factor, written from the algebra and checked against z_s
        CHECK_CLOSE(1.0 - z_s, y * (1.0 - x_a) / (1.0 - x_a * y), 1e-12);
        const double d = ll_radiator(z_s, q2);
        const RcModel::TailTriple a = mt.tail_sigma_at(x, q2);
        const RcModel::TailTriple b = mp.tail_sigma_at(x, q2);
        const double r = (b.u + b.qe + b.u_sp + b.qe_sp) / (a.u + a.qe);
        MESSAGE("Q^2 = " << q2 << ", x = " << x << ", y = " << y
                         << ": D(z_s) = " << d << ", ratio = " << r);
        CHECK((d > 1.0) == (r > 1.01));    // the mechanism, as an iff
        CHECK(d > last_d);                 // both monotone down in y
        CHECK(r > last_r);
        last_d = d;
        last_r = r;
      }
    }

    // (4) THE TENSOR CONSEQUENCE at the bad cell.  The tensor numerator is the
    //     t-peak's own and does not move, so the tensor FRACTION of the tail
    //     collapses by the same factor the denominator grew.  That is a
    //     four-orders-of-magnitude change INSIDE the window the documents used
    //     to declare safe to 0.16 %.
    {
      const double qn = 600.0;
      const double dt = mt.tail_ratio_at(x_bad, q2_bad, qn) -
                        mt.tail_ratio_at(x_bad, q2_bad, 0.0);
      const double dp = mp.tail_ratio_at(x_bad, q2_bad, qn) -
                        mp.tail_ratio_at(x_bad, q2_bad, 0.0);
      const double ft = dt / mt.tail_ratio_at(x_bad, q2_bad, 0.0);
      const double fp = dp / mp.tail_ratio_at(x_bad, q2_bad, 0.0);
      MESSAGE("tensor fraction of the tail at the bad cell: " << ft
              << " (t-peak) -> " << fp << " (t-peak+ll), factor " << fp / ft);
      REQUIRE(std::fabs(ft) > 0.0);
      CHECK(std::fabs(fp) < 1e-3 * std::fabs(ft));
      // The ABSOLUTE tensor term is untouched -- (c)'s gate, at the point
      // where it loses the most digits (base/|d| ~ 4.5e11, so the
      // double-precision floor on this subtraction is ~1e-4 relative).
      CHECK_CLOSE(dp, dt, 1e-3);
    }
  }
}

// -------------------------------------------------------- the odds and ends

TEST_CASE("T19: RcTailModel::PolradFull -- POLRAD Eq. (18) + App. B + Eq. (A.4)") {
  // THE THREE GATES design sec. 5 / the B2 brief ask for, in order:
  //   (a) the limit in which Eq. (18) MUST reduce to Eq. (38) -- stated, and
  //       it is not "switch the s-/p-peaks off" by hand: it is x_A -> 0 with a
  //       form factor DEAD at the s-/p-peak vertex, which is what makes
  //       Eq. (38)'s ultrarelativistic t-peak extraction exact.
  //   (b) T9 (its own TEST_CASE above), un-skipped: the Q_N = 0 Rosenbluth
  //       limit of Eq. (A.4).
  //   (c) the transcription check's own tabulated numbers -- T8(c') for the
  //       t-peak, extended here.
  // and then what the two models actually differ by.
  const double m_d = DEUTERON().mass();
  const double s_a_d = 2.0 * (2.0 * PROTON_MASS * 27.6);   // S_A = A s

  SUBCASE("(a) THE LIMIT: x_A -> 0 with the s-/p-peaks dead, Eq. (18) -> "
          "Eq. (38), UNPOLARISED AND TENSOR") {
    // WHY THE FORM FACTOR IS THE SWITCH.  Eq. (38) is POLRAD's own
    // ultrarelativistic t-peak reduction of Eq. (18) (sec. 2.1.3 B), so the
    // two agree only where the OTHER two peaks are absent.  A Gaussian with
    // b = 20 t_min is alive across the whole t-peak and 20+ decades down at
    // the s-peak's Q'^2 ~ z_s Q^2, which isolates the t-peak WITHOUT touching
    // the quadrature.  What is left is the ultrarelativistic error, and it is
    // O(x_A): MEASURED 1 + 1.64 x_A on the spin-0 sector.
    struct Row { double x_a, want_u, want_t; };
    // spin-0 (F_c alone): no tensor column at all.
    const Row spin0[] = {{0.003, 1.00492, 0.0}, {0.006, 1.00990, 0.0},
                         {0.012, 1.01997, 0.0}};
    for (const Row& r : spin0) {
      const double b = 20.0 * (r.x_a * 0.938) * (r.x_a * 0.938);
      const GaussianSpin1FF ff(b, 0.0, 0.0, 1.0);
      const double e38 = polrad_sigma_el_u(ff, r.x_a, 0.5, s_a_d, m_d, 1024);
      const PolradFullPair p =
          polrad_full_sigma_el(ff, r.x_a, 0.5, s_a_d, m_d, 128);
      CHECK_CLOSE(p.u / e38, r.want_u, 2e-4);
      // ... and sigma_q is IDENTICALLY zero when F_m = F_q = 0: every Q_N
      // term of Eq. (A.4) carries one of them, so `== 0.0` is the right
      // assertion and anything weaker would let F_c leak into the tensor tail.
      CHECK(p.t == 0.0);
      CHECK(polrad_sigma_el_t(ff, r.x_a, 0.5, s_a_d, m_d, 1024) == 0.0);
    }
    // F_m alone, and F_q alone -- the two sectors of sigma_q^d that do NOT
    // cancel against each other.  (The full deuteron does: sigma_q^d is a
    // ~7x cancellation between them there, so its RATIO is not a gate.)
    const Row fm[] = {{0.003, 1.00230, 1.00313}, {0.006, 1.00462, 1.00627},
                      {0.012, 1.00931, 1.01257}};
    const Row fq[] = {{0.003, 1.00670, 1.01872}, {0.006, 1.01347, 1.03850},
                      {0.012, 1.02722, 1.08168}};
    for (int sector = 0; sector < 2; ++sector) {
      const Row* rows = sector == 0 ? fm : fq;
      for (int i = 0; i < 3; ++i) {
        const Row& r = rows[i];
        const double b = 20.0 * (r.x_a * 0.938) * (r.x_a * 0.938);
        const GaussianSpin1FF ff(b, sector == 0 ? 1.7139610634 : 0.0,
                                 sector == 0 ? 0.0 : 25.84, 0.0);
        const double e38u = polrad_sigma_el_u(ff, r.x_a, 0.5, s_a_d, m_d, 1024);
        const double e38t = polrad_sigma_el_t(ff, r.x_a, 0.5, s_a_d, m_d, 1024);
        const PolradFullPair p =
            polrad_full_sigma_el(ff, r.x_a, 0.5, s_a_d, m_d, 128);
        CHECK_CLOSE(p.u / e38u, r.want_u, 2e-4);
        CHECK_CLOSE(p.t / e38t, r.want_t, 2e-4);
      }
    }
    // THE LIMIT ITSELF, as a monotone statement and not three numbers: the
    // deviation from Eq. (38) must SHRINK as x_A does, in both columns.
    double prev_u = 0.0, prev_t = 0.0;
    for (double xa : {0.012, 0.006, 0.003, 0.0015}) {
      const double b = 20.0 * (xa * 0.938) * (xa * 0.938);
      const GaussianSpin1FF ff(b, 1.7139610634, 0.0, 0.0);
      const double e38u = polrad_sigma_el_u(ff, xa, 0.5, s_a_d, m_d, 1024);
      const double e38t = polrad_sigma_el_t(ff, xa, 0.5, s_a_d, m_d, 1024);
      const PolradFullPair p = polrad_full_sigma_el(ff, xa, 0.5, s_a_d, m_d, 128);
      const double du = std::fabs(p.u / e38u - 1.0);
      const double dt = std::fabs(p.t / e38t - 1.0);
      if (prev_u > 0.0) { CHECK(du < prev_u); CHECK(dt < prev_t); }
      prev_u = du; prev_t = dt;
    }
    CHECK(prev_u < 2e-3);            // measured 1.15e-03 at x_A = 0.0015
    CHECK(prev_t < 3e-3);            // measured 1.57e-03
  }

  SUBCASE("(b) the tau_A range and where the peaks sit inside it") {
    // tau_s = -Q^2/S and tau_p = Q^2/X are PANEL EDGES of the quadrature, so
    // if either fell outside [tau_min, tau_max] a whole peak would be
    // integrated with no clustering and silently under-resolved.
    const double s_a = 6.0 * 3980.0, m_a = LI6().mass();
    for (double x : {1e-4, 0.01, 0.3, 0.9}) {
      for (double y : {0.01, 0.5, 0.985}) {
        const double x_a = x / 6.0, q2 = x_a * y * s_a;
        const TauLimits t = polrad_tau_limits(x_a, y, s_a, m_a);
        const double tas = -q2 / s_a, tap = q2 / ((1.0 - y) * s_a);
        CHECK(t.hi > t.lo);
        CHECK(t.lo > -1.0);                     // R_el = (S_x-Q^2)/(1+tau)
        CHECK(t.lo < tas);
        CHECK(tas < 0.0);
        CHECK(tap > 0.0);
        CHECK(tap < t.hi);
        // tau_min tau_max = -Q^2/M^2 exactly (Eq. (14)); the stable form.
        CHECK_CLOSE(t.lo * t.hi, -q2 / (m_a * m_a), 1e-12);
      }
    }
  }

  SUBCASE("(c) POLRAD's own three deuteron points -- the SIGN and the scale") {
    // T8(c') gates the t-peak against polrad_transcription_check.md sec. 8.
    // The exact tail must land on the SAME SIDE of zero in the tensor column
    // at every one of them (that column CHANGES SIGN with x, which is the
    // whole reason sigma^el_T is a separate table) and within a factor of a
    // few in the unpolarised one -- it adds two peaks, it does not replace
    // the object.  MEASURED full/t-peak: 1.065 / 1.176 / 1.547 unpolarised,
    // 1.197 / 2.571 / 1.076 tensor.
    const std::shared_ptr<HoSpin1FF> ff = deuteron_ff();
    struct Row { double e, x, y, ru, rt; };
    const Row rows[] = {{27.6, 0.050, 0.60, 1.0647, 1.1967},
                        {27.6, 0.012, 0.50, 1.1760, 2.5712},
                        {11.0, 0.200, 0.50, 1.5474, 1.0758}};
    for (const Row& r : rows) {
      const double s_n = 2.0 * PROTON_MASS * r.e, s_a = 2.0 * s_n;
      const double x_a = r.x / 2.0;
      const double eu = polrad_sigma_el_u(*ff, x_a, r.y, s_a, m_d, 256);
      const double et = polrad_sigma_el_t(*ff, x_a, r.y, s_a, m_d, 256);
      const PolradFullPair p = polrad_full_sigma_el(*ff, x_a, r.y, s_a, m_d, 128);
      CHECK(p.u > 0.0);                                    // T8(0)'s sign gate
      CHECK(p.t * et > 0.0);                               // same side of zero
      CHECK_CLOSE(p.u / eu, r.ru, 2e-3);
      CHECK_CLOSE(p.t / et, r.rt, 2e-3);
      CHECK(p.u > eu);                        // it ADDS the s- and p-peaks
    }
  }

  SUBCASE("(d) the LEADING-LOG fallback, at the HERMES deuteron point") {
    // x = 0.012, y = 0.85, Q^2 = 0.53 -- the point the header block quotes as
    // "the t-peak is 23 % of the leading-log total (low by 4.36x)".  The exact
    // tail settles it: it is low by 2.36x, i.e. the leading-log edge
    // OVERSHOOTS by 1.85x.  That is the whole reason this model was built, so
    // it is gated rather than only reported.
    const std::shared_ptr<HoSpin1FF> ff = deuteron_ff();
    const double s_n = 2.0 * PROTON_MASS * 27.6, s_a = 2.0 * s_n;
    const double x = 0.012, y = 0.85, x_a = x / 2.0;
    const double te = polrad_sigma_el_u(*ff, x_a, y, s_a, m_d, 256) / 4.0;
    const LlPeaks se = ll_peaks_spin1(*ff, x_a, y, s_a, m_d);
    const double tq = polrad_sigma_qe_u(1, 1, x, y, s_n, PROTON_MASS, 256, 0.0)
                      / 2.0;
    const LlPeaks sq = ll_peaks_qe(1, 1, x, y, s_n);
    const PolradFullPair fe = polrad_full_sigma_el(*ff, x_a, y, s_a, m_d, 128);
    const double fq = polrad_full_sigma_qe_u(1, 1, x, y, s_n, PROTON_MASS, 128,
                                             0.0) / 2.0;
    const double tot_t = te + tq;
    const double tot_ll = tot_t + (se.s + se.p) / 4.0 + (sq.s + sq.p) / 2.0;
    const double tot_f = fe.u / 4.0 + fq;
    MESSAGE("HERMES deuteron point: t-peak " << tot_t << ", t-peak+ll "
            << tot_ll << ", polrad-full " << tot_f << "  (full/t = "
            << tot_f / tot_t << ", full/(t+ll) = " << tot_f / tot_ll << ")");
    CHECK_CLOSE(tot_f / tot_t, 2.3562, 2e-3);
    CHECK_CLOSE(tot_f / tot_ll, 0.5401, 2e-3);
    // ... and it is BETWEEN the two edges here, which is the case the band was
    // built for.  (T19(g) is the case it is NOT.)
    CHECK(tot_f > tot_t);
    CHECK(tot_f < tot_ll);
  }

  SUBCASE("(e) the quadrature CONVERGES, and the tensor column is the one "
          "that needs the resolution") {
    const std::shared_ptr<HoSpin1FF> ff = HoSpin1FF::for_ion(LI6());
    const double s_a = 6.0 * 3980.0, m_a = LI6().mass();
    for (double x : {1e-3, 1e-4}) {
      for (double y : {0.5, 0.9}) {
        const PolradFullPair a =
            polrad_full_sigma_el(*ff, x / 6.0, y, s_a, m_a, 128);
        const PolradFullPair b =
            polrad_full_sigma_el(*ff, x / 6.0, y, s_a, m_a, 512);
        // 1e-8: MEASURED worst spread over the four points is 2.13e-09
        // (x = 1e-4, y = 0.5, 35.95009232 at n_tau = 128 against 35.95009240
        // at 512).  The unpolarised column is converged; the tensor one is
        // not, and the two tolerances below say by how much.
        CHECK_CLOSE(a.u, b.u, 1e-8);
        // 1e-3, not 1e-9: Eq. (B.3)'s a_ik cancel to ~5 decimal digits in the
        // collinear region and `long double` leaves ~1e-4 of that on the
        // tensor column.  rc.cpp says so where the type is chosen; this is
        // the number.
        CHECK_CLOSE(a.t, b.t, 1e-3);   // measured worst 6.0e-04
        // ... and the UNPOLARISED column is genuinely converged, which is what
        // makes the tensor spread arithmetic and not quadrature.
        CHECK(std::fabs(a.u / b.u - 1.0) < 1e-8);
      }
    }
  }

  SUBCASE("(f) THE QUASI-ELASTIC TAIL IS STILL TENSOR-BLIND -- the gap that "
          "the exact tail does NOT close") {
    // THIS SUBCASE ASSERTS AN OMISSION, exactly like B1's
    // `..._bounded_by_NEITHER_scale` pytest, and for the same reason: a
    // reader could reasonably expect Eq. (18) to close the hole
    // `RcOptions::sp_tensor_scale` only bounds.  It closes the ELASTIC half
    // and not the quasi-elastic one, because a nucleon has no tensor
    // structure function to put at any peak (Eq. (A.5) -> Im_{5..8} = 0).
    //
    // The assertion: `polrad_full_sigma_qe_u` has no tensor entry point AT
    // ALL, and the elastic tensor column is untouched by every quasi-elastic
    // knob -- so the quasi-elastic tensor tail is exactly zero on this model
    // as on the other two, and only `qe_tensor_scale` prices it.
    const double s_n = 3980.0;
    const double a = polrad_full_sigma_qe_u(3, 3, 0.1, 0.5, s_n, PROTON_MASS,
                                            128, 0.169);
    const double b = polrad_full_sigma_qe_u(3, 3, 0.1, 0.5, s_n, PROTON_MASS,
                                            128, 0.0);
    CHECK(a > 0.0);
    CHECK(b > a);                       // Pauli suppression bites, as Eq. (44)
    // ... and it is BIG: the quasi-elastic s+p is what the t-peak misses at
    // x >= 0.1, and none of it is tensor.  MEASURED at 6Li config 1,
    // x = 0.1, Q^2 = 5, k_F = 0.169: full/t-peak = 202.2.
    const double y = 5.0 / (0.1 * s_n);
    const double tqe = polrad_sigma_qe_u(3, 3, 0.1, y, s_n, PROTON_MASS, 128,
                                         0.169);
    const double fqe = polrad_full_sigma_qe_u(3, 3, 0.1, y, s_n, PROTON_MASS,
                                              128, 0.169);
    MESSAGE("quasi-elastic at x = 0.1, Q^2 = 5: t-peak " << tqe
            << " -> polrad-full " << fqe << " (x" << fqe / tqe << "), and "
            "ALL of that growth is tensor-blind");
    CHECK(fqe / tqe > 100.0);
  }

  SUBCASE("(g) it does NOT sit inside the t-peak band, and the model moves "
          "the tail through RcModel") {
    RcOptions o;
    o.tail_model = RcTailModel::PolradFull;
    const RunPlan plan = tensor_thirds_plan(0.0, 0.6);
    const RcModel mf(RcMode::TensorBand, o, rc_sampler_ptr(), plan,
                     Channel::Inclusive, LI6());
    RcOptions ol;
    ol.tail_model = RcTailModel::TPeakPlusLL;
    const RcModel ml(RcMode::TensorBand, ol, rc_sampler_ptr(), plan,
                     Channel::Inclusive, LI6());
    const RcModel& mt = rc_model();
    // The tables are NOT the s+p decomposition: under PolradFull `u` and `qe`
    // already contain the s- and p-peaks, so u_sp/qe_sp stay zero for a
    // DIFFERENT reason than under TPeak.
    const RcModel::TailTriple tt = mf.tail_sigma_at(0.01, 5.0);
    CHECK(tt.u_sp == 0.0);
    CHECK(tt.qe_sp == 0.0);
    CHECK(tt.u > mt.tail_sigma_at(0.01, 5.0).u);
    // ... and there is at least one accepted cell where the exact tail is
    // OUTSIDE the [TPeak, TPeakPlusLL] interval, so the pair is a price range
    // and not a confidence interval.  MEASURED on the PRODUCTION grid (this
    // loop runs on the small test sampler, so it gates the sign of the
    // statement and not its size): 1326 of 3051 cells are outside, and the
    // worst is x4518.3 the t-peak at x = 0.954993, Q^2 = 206.68, y = 0.0544
    // with `tail_max` raised so the census measures the MODEL.  At the
    // shipped `tail_max` = 10 that cell and five others clip to x11 = 1 +
    // `tail_max` -- the CEILING, which is what "the worst x11 the t-peak at
    // x = 0.955" meant here until 2026-09-06; that cell is x3449.9 unclipped.
    int outside = 0, n = 0;
    for (std::size_t c = 0; c < rc_sampler().n_cells(); ++c) {
      const double x = rc_sampler().x_cells()[c];
      const double q2 = rc_sampler().q2_cells()[c];
      const double a = mt.tail_ratio_at(x, q2, 0.0);
      const double b = ml.tail_ratio_at(x, q2, 0.0);
      const double f = mf.tail_ratio_at(x, q2, 0.0);
      ++n;
      if (f < std::min(a, b) || f > std::max(a, b)) ++outside;
    }
    MESSAGE("polrad-full outside the [t-peak, t-peak+ll] interval on "
            << outside << " of " << n << " test-sampler cells");
    CHECK(outside > 0);
  }

  SUBCASE("(h) m_lepton is READ here and by nothing else") {
    // The knob `PipelineConfig::validate()` reserved for this model since the
    // model was a comment.  It is the m^2 of C_{1,2}(tau), so it sets the
    // WIDTH of the s-/p-peaks: a heavier lepton makes them narrower and the
    // tail SMALLER.  Gated as an inequality, not a number -- what matters is
    // that the field reaches the quadrature at all.
    const std::shared_ptr<HoSpin1FF> ff = HoSpin1FF::for_ion(LI6());
    const double s_a = 6.0 * 3980.0, m_a = LI6().mass();
    const double x_a = 1e-4 / 6.0, y = 0.9;
    const PolradFullPair e =
        polrad_full_sigma_el(*ff, x_a, y, s_a, m_a, 128, M_ELECTRON);
    const PolradFullPair mu =
        polrad_full_sigma_el(*ff, x_a, y, s_a, m_a, 128, 0.1056583755);
    CHECK(e.u != mu.u);
    CHECK(mu.u < e.u);
    MESSAGE("m_lepton e -> mu at x = 1e-4, y = 0.9: sigma^el_U "
            << e.u << " -> " << mu.u << " (x" << mu.u / e.u << ")");
    // ... and the two t-peak entry points take no lepton mass at all, which
    // is why validate() refuses the knob on them.
    CHECK(polrad_sigma_el_u(*ff, x_a, y, s_a, m_a, 128) ==
          polrad_sigma_el_u(*ff, x_a, y, s_a, m_a, 128));
  }

  SUBCASE("(i) THE PIN -- tail_sigma_at under PolradFull is pinned to a "
          "VALUE, not to a ratio or an inequality") {
    // T8(a') does this for the t-peak, and the reason it exists is that the
    // per-nucleon reduction at THIS call site is where the historical factor-A
    // bug lived (polrad_transcription_check.md sec. 4: the shipped
    // `per_nucleon = m_p/M_A` was 6.00x too large).  Until this subcase,
    // NOTHING pinned a value on `RcTailModel::PolradFull`: T19(a)-(h) and T9
    // test `polrad_full_sigma_el` / `..._qe_u` directly, T14 and
    // test_rc_pipeline check only refusals, and the one assertion that touched
    // the model's own PolradFull branch was T19(g)'s inequality
    // `tt.u > mt.tail_sigma_at(0.01, 5.0).u`.  MEASURED: replacing the
    // branch's `1.0/(a*a)` by `1.0/a` -- exactly the historical bug, at the
    // exact line it lived on -- left T9 + T19 at 316/316 and both polrad-full
    // pytests passing.  It cannot now: every number below moves by x6.
    //
    // THE SECOND SLIP THE SAME CALL CAN MAKE is passing `x` where `x_A = x/A`
    // belongs (and `s` where `S_A = A s` does).  That was caught by exactly
    // one inequality; it is now caught by the values, and by how far the wrong
    // argument lands from them -- printed below, x0.0052 / x3.6e-10 / x9.4e-62
    // at the three points.
    RcOptions o;
    o.tail_model = RcTailModel::PolradFull;
    const RcModel mf(RcMode::TensorBand, o, rc_sampler_ptr(),
                     tensor_thirds_plan(0.0, 0.6), Channel::Inclusive, LI6());
    const double s = rc_sampler().s();
    const double a = static_cast<double>(LI6().A);
    const double m_a = LI6().mass();
    const std::shared_ptr<HoSpin1FF> ff = HoSpin1FF::for_ion(LI6());
    // 6Li config 1 (`default_configs("6Li")[1]`, S = 3980 GeV^2), RcOptions()
    // defaults but for the model: n_eta = 128, m_lepton = m_e,
    // k_F = RC_QE_KF_GEV, with_qe_tail = true.  All three columns are
    // sigma per nucleon [GeV^-2], i.e. Eq. (18)/A^2 and Eq. (A.5)/A.
    //   (0.01,  5) the reviewer's point: u = 2.230125e-05 is the raw
    //              8.0284510282e-04 divided by A^2 = 36, and NOT by A = 6.
    //   (0.1,  20) x_A = 1/60 is far from x = 0.1 on a form factor that has
    //              fallen 7.8 decades between them: the point where the
    //              x-vs-x_A slip is unmissable.
    //   (0.3,   5) the QUASI-ELASTIC point: qe is 3.7e+06 times u here, so
    //              the third column is pinned where it is the whole tail --
    //              and it carries ONE 1/A, not two, which is the other half
    //              of the same historical bug (the ERT and the QRT were
    //              normalised differently, so their ratio was wrong by A).
    struct Pin { double x, q2, u, t, qe; };
    const Pin pins[] = {
        {0.01,  5.0,  2.2301252856028388e-05, -6.8964416598016638e-08,
                      1.4082400494855565e-05},
        {0.10, 20.0,  1.3474277075174561e-08,  1.0398154385161343e-11,
                      4.2214155843320782e-08},
        {0.30,  5.0,  2.5848200248542323e-12,  1.6976986824052815e-13,
                      9.5240888233400158e-06}};
    for (const Pin& p : pins) {
      const RcModel::TailTriple tt = mf.tail_sigma_at(p.x, p.q2);
      CHECK_CLOSE(tt.u,  p.u,  1e-9);
      CHECK_CLOSE(tt.t,  p.t,  1e-9);
      CHECK_CLOSE(tt.qe, p.qe, 1e-9);
      // ... and the pinned values ARE the raw quadratures over A^2 (elastic)
      // and A (quasi-elastic), so a reader can rebuild them from the header's
      // own reduction rather than trusting three literals.
      const double y = p.q2 / (p.x * s);
      const PolradFullPair raw =
          polrad_full_sigma_el(*ff, p.x / a, y, a * s, m_a, o.n_eta,
                               o.m_lepton);
      const double raw_qe =
          polrad_full_sigma_qe_u(LI6().Z, LI6().N(), p.x, y, s, PROTON_MASS,
                                 o.n_eta, o.qe_kf_gev, o.m_lepton);
      CHECK_CLOSE(tt.u  * (a * a), raw.u,  1e-12);
      CHECK_CLOSE(tt.t  * (a * a), raw.t,  1e-12);
      CHECK_CLOSE(tt.qe * a,       raw_qe, 1e-12);
      // The two wrong reductions and the wrong argument, so the failure a
      // future edit produces is legible rather than just red.
      const PolradFullPair wrong_x =
          polrad_full_sigma_el(*ff, p.x, y, a * s, m_a, o.n_eta, o.m_lepton);
      MESSAGE("PIN x = " << p.x << ", Q^2 = " << p.q2 << ": u = " << tt.u
              << " (1/A would give " << raw.u / a << ", x" << a << "), "
              "x-for-x_A would give " << wrong_x.u / (a * a) << " (x"
              << wrong_x.u / raw.u << ")");
      CHECK(std::fabs(raw.u / a / p.u - a) < 1e-9);      // 1/A is exactly x6
      CHECK(wrong_x.u / raw.u < 0.01);                   // and x_A is not x
    }
    // The s/p columns are zero HERE for the reason T19(g) states, and the
    // pin says so as a value too: nothing leaks into them from Eq. (18).
    CHECK(mf.tail_sigma_at(0.01, 5.0).u_sp == 0.0);
    CHECK(mf.tail_sigma_at(0.01, 5.0).qe_sp == 0.0);
  }
}

TEST_CASE("TabulatedSpin1FF: the API exists, and NO 6Li table is shipped") {
  // design sec. 2.1's fallback (F) is a digitised `data/ff/li6_elastic.csv`.
  // NOTHING IS SHIPPED, on purpose: the recommended source, Wiringa &
  // Schiavilla (nucl-th/9807037), is FIGURES ONLY -- "there is no table
  // anywhere in the paper" -- so there are no digitised numbers WITH
  // PROVENANCE to put in a file, and inventing some would be worse than
  // having none.  The reader is implemented and its refusals are gated; it is
  // NOT gated against a file, and this test is where that is recorded.
  HoSpin1FFOptions norm;
  norm.m_a_gev = LI6().mass();
  CHECK_THROWS_AS(TabulatedSpin1FF::from_data_dir("ff/li6_elastic.csv", norm),
                  std::runtime_error);
  try {
    TabulatedSpin1FF::from_data_dir("ff/li6_elastic.csv", norm);
  } catch (const std::runtime_error& e) {
    const std::string msg = e.what();
    CHECK(msg.find("FIGURES ONLY") != std::string::npos);
    CHECK(msg.find("HoSpin1FF::for_ion") != std::string::npos);
  }
  // M_A is not optional: F_C2 -> F_q is eta_A-dependent and eta_A needs it.
  CHECK_THROWS_AS(
      TabulatedSpin1FF::from_data_dir("ff/li6_elastic.csv",
                                      HoSpin1FFOptions()),
      std::runtime_error);
  // The path goes through data_dir(), never the working directory
  // (docs/CONVENTIONS.md), so a bare relative name is resolved under it.
  CHECK_THROWS_AS(TabulatedSpin1FF::from_data_dir("no/such/table.csv", norm),
                  std::runtime_error);
}

TEST_CASE("RcScope::TensorAll buys back the Delta cos 2phi term, and only it") {
  // NO LITERATURE SUPPORT: nobody has computed RC for a phi-dependent tensor
  // observable and Delta is not in POLRAD's b1..b4 basis.  The switch exists
  // to PRICE the omission, never to correct it -- so what it must do is add
  // exactly the Delta term and nothing else.
  RcOptions all;
  all.scope = RcScope::TensorAll;
  const InclusiveSampler& s = rc_sampler();
  // theta_S = 0 kills Delta (it carries sin^2 theta_S), so this needs an axis
  // off the beam -- which is also where the tail's phi-integrated caveat bites.
  const std::vector<SpinCategory> cats = {
      SpinCategory("tilt", 1.0, {1.0, 0.0, 0.0}, 0, 0.0, 0.6, 0.25)};
  const RcModel base(RcMode::TensorBand, RcOptions(), rc_sampler_ptr(),
                     RunPlan(cats), Channel::Inclusive, LI6());
  const RcModel wide(RcMode::TensorBand, all, rc_sampler_ptr(), RunPlan(cats),
                     Channel::Inclusive, LI6());
  int moved = 0;
  for (std::size_t c = 0; c < s.n_cells(); c += 23) {
    for (double phi : {0.0, 1.1, 2.7}) {
      Event ev = inclusive_event(s, static_cast<int>(c), 1.0, phi, 0.6);
      ev.spin.phi_s = 0.25;
      const double tb = base.tensor_fraction(ev);
      const double tw = wide.tensor_fraction(ev);
      const InclusiveSampler::StateTables& st =
          s.state_tables(cats[0], 1.0);
      const double phip = phi - 0.25;
      const double den = 1.0 + st.w_avg[c] + st.a1[c] * std::cos(phip) +
                         st.a2[c] * std::cos(2.0 * phip);
      const double delta_term =
          (st.a2[c] - st.a2_tensor[c]) * std::cos(2.0 * phip) / den;
      CHECK_CLOSE(tw - tb, delta_term, 1e-12);
      if (std::fabs(delta_term) > 1e-14) ++moved;
    }
  }
  CHECK(moved > 5);
}

TEST_CASE("tail_ratio(cell, Q_N) is tail_ratio_at at that cell's centre") {
  const RcModel& m = rc_model();
  const InclusiveSampler& s = rc_sampler();
  for (std::size_t c = 0; c < s.n_cells(); c += 53) {
    for (double q_n : {0.0, 1.0, -2.0}) {
      CHECK(m.tail_ratio(static_cast<int>(c), q_n) ==
            m.tail_ratio_at(s.x_cells()[c], s.q2_cells()[c], q_n));
    }
  }
  CHECK_THROWS_AS(m.tail_ratio(-1, 0.0), std::runtime_error);
  CHECK_THROWS_AS(m.tail_ratio(static_cast<int>(s.n_cells()), 0.0),
                  std::runtime_error);
}
