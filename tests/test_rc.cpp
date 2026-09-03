// Tensor-sector radiative corrections (rc.hpp), the OPT-IN weight-only family
// of docs/open_items/run_2026-09-02/design_C_tensor_rc.md.
//
// Design sec. 5's numbering is kept, so a reader can go row by row:
//
//   T1   the band identity against the PUBLISHED A_zz          T8'  the QRT, Eq. (44)
//   T2   the tagged identity (the SAME closed form)            T9   SKIP-ped, with the reason
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
#include <memory>
#include <string>
#include <thread>
#include <vector>

#include "check_close.hpp"
#include "doctest.h"
#include "lipolgen/beams.hpp"
#include "lipolgen/bookkeeping.hpp"
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
    CHECK_CLOSE(ff->fc(0.0), static_cast<double>(LI6().Z), 1e-14);
    CHECK_CLOSE(ff->fm(0.0), (m_a / PROTON_MASS) * LI6_MU_N, 1e-14);
    CHECK_CLOSE(ff->fq(0.0),
                std::pow(m_a / HBARC_GEV_FM, 2) * LI6_QUADRUPOLE_FM2, 1e-14);
    // The design's printed values (sec. 2.1), at its own precision.
    CHECK_CLOSE(ff->fc(0.0), 3.0, 1e-14);
    CHECK_CLOSE(ff->fm(0.0), 4.90765, 1e-4);
    CHECK_CLOSE(ff->fq(0.0), -65.914, 1e-4);
    // TUNL's -0.0818(17) fm^2, NOT Pyykko's -0.0806 (which gives -64.947).
    CHECK(ff->fq(0.0) < -65.0);
  }

  SUBCASE("<r^2>_point = 6.078 fm^2 from the POINT shape") {
    // -6 dF/dq^2|_0 / F(0) on the UNFOLDED shape: with `fold_nucleon` the
    // slope would also carry <r^2>_p and <r^2>_n, which is what <r^2>_point
    // was derived by SUBTRACTING (design sec. 2.1).
    HoSpin1FFOptions o;
    o.fold_nucleon = false;
    const std::shared_ptr<HoSpin1FF> pt = HoSpin1FF::for_ion(LI6(), o);
    const double h = 1e-6;                      // in q^2 [fm^-2]
    const double t = h * HBARC_GEV_FM * HBARC_GEV_FM;
    const double slope = (pt->fc(t) - pt->fc(0.0)) / h;
    // 6.0788 is the design's own full value; its table rounds it to 6.078.
    CHECK_CLOSE(-6.0 * slope / pt->fc(0.0), 6.0788, 1e-4);
    // ... and the closed form the design derives it from.
    const double a = LI6_FF_HO_A_FM, al = LI6_FF_HO_ALPHA;
    CHECK_CLOSE(1.5 * a * a * (2.0 + 5.0 * al) / (2.0 + 3.0 * al), 6.0788, 1e-4);
  }

  SUBCASE("the first C0 zero lies in [2.9, 3.3] fm^-1") {
    // NOT pinned to 3.10 +- 0.02: WS98 does not locate the zero (it says only
    // that C2 dominates above 3 fm^-1 and that two-body currents shift the
    // minimum lower, with no number), so the design asks for a REFIT and this
    // test pins only the range the refit must land in.
    const double a = LI6_FF_HO_A_FM, al = LI6_FF_HO_ALPHA;
    const double q0 = std::sqrt(2.0 * (2.0 + 3.0 * al) / al) / a;
    CHECK(q0 > 2.9);
    CHECK(q0 < 3.3);
    const double t0 = std::pow(q0 * HBARC_GEV_FM, 2);
    CHECK(std::fabs(ff->fc(t0)) < 1e-10 * ff->fc(0.0));
    // The shell-model alpha = (Z-2)/3 = 1/3 would put it at 2.33 fm^-1 -- far
    // too low, which is why alpha is FREE and this is a phenomenological fit.
    CHECK(std::sqrt(2.0 * (2.0 + 1.0) / (1.0 / 3.0)) / a < 2.5);
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
// THE LEADING-LOG s- AND p-PEAKS -- T8(c)'s ABSOLUTE gate on the t-peak.
//
// POLRAD sec. 2.1.3 B asserts that for a TAIL the s- and p-peaks are
// suppressed, and rc.cpp's whole tail rests on that one sentence.  Nothing in
// this file measured it, so this block does: the standard collinear
// (equivalent-radiator) estimate of the SAME observable, built from the
// elastic cross section and the Weizsacker-Williams lepton radiator, sharing
// no line with rc.cpp.
//
//   ISR (s-peak): the incoming lepton radiates, k1 -> z k1, and elasticity
//   fixes z = (1-y)/(1 - x_A y).  With Q'^2 = z Q^2, S' = z S_A and the
//   Jacobian |d(Q^2,y)/d(Q'^2,z)| = (1 - x_A y)/z,
//       d^2 sigma/(dx_A dy) = y S_A D(z) [d sigma_el/dQ'^2] z/(1 - x_A y) .
//   FSR (p-peak): k2 -> k2/z' with z' = 1 - y(1 - x_A), Q'^2 = Q^2/z',
//   S' = S_A and |J| = z', so the same product with 1/z'.
//   D(z) = (alpha/pi) ln(Q^2/m_e^2) (1 + z^2)/(1 - z).
double ll_radiator(double z, double q2) {
  if (!(z > 0.0) || !(z < 1.0)) return 0.0;
  return ALPHA_EM / kPi * std::log(q2 / (M_ELECTRON * M_ELECTRON)) *
         (1.0 + z * z) / (1.0 - z);
}

/// Rosenbluth (A, B) of a spin-1 target -- POLRAD Eq. (A.4) at Q_N = 0.
void rosenbluth_spin1(const Spin1ElasticFF& ff, double q2, double m,
                      double* a, double* b) {
  const double eta = q2 / (4.0 * m * m);
  const double fc = ff.fc(q2), fm = ff.fm(q2), fq = ff.fq(q2);
  *a = fc * fc + (8.0 / 9.0) * eta * eta * fq * fq +
       (2.0 / 3.0) * eta * fm * fm;
  *b = (4.0 / 3.0) * eta * (1.0 + eta) * fm * fm;
}
void rosenbluth_nucleon(bool proton, double q2, double* a, double* b) {
  const NucleonFF f = nucleon_ff(q2);
  const double ge = proton ? f.ge_p : f.ge_n;
  const double gm = proton ? f.gm_p : f.gm_n;
  const double tau = q2 / (4.0 * PROTON_MASS * PROTON_MASS);
  *a = (ge * ge + tau * gm * gm) / (1.0 + tau);
  *b = 2.0 * tau * gm * gm;
}

/// d sigma_el/dQ'^2 = (4 pi alpha^2/Q'^4)[ (y'^2/2) B/(2 eta')
///                                       + (1 - y' - m^2 y'^2/Q'^2) A ],
/// y' = Q'^2/S', eta' = Q'^2/(4m^2).  (Elastic F_1 = B/(4 eta), F_2 = A.)
double dsigma_el_dq2(double a, double b, double q2p, double sp, double m) {
  if (!(q2p > 0.0) || !(sp > 0.0)) return 0.0;
  const double yp = q2p / sp;
  const double etap = q2p / (4.0 * m * m);
  const double kin = 1.0 - yp - m * m * yp * yp / q2p;
  if (!(kin > 0.0)) return 0.0;
  return 4.0 * kPi * ALPHA_EM * ALPHA_EM / (q2p * q2p) *
         (0.5 * yp * yp * b / (2.0 * etap) + kin * a);
}

/// (s-peak, p-peak) of a spin-1 nucleus, WHOLE-NUCLEUS d^2 sigma/(dx_A dy).
std::pair<double, double> ll_peaks_spin1(const Spin1ElasticFF& ff, double x_a,
                                         double y, double s_a, double m) {
  const double q2 = x_a * y * s_a;
  const double zs = (1.0 - y) / (1.0 - x_a * y);
  const double zp = 1.0 - y * (1.0 - x_a);
  double as, bs, ap, bp;
  rosenbluth_spin1(ff, zs * q2, m, &as, &bs);
  rosenbluth_spin1(ff, q2 / zp, m, &ap, &bp);
  const double s = y * s_a * ll_radiator(zs, q2) *
                   dsigma_el_dq2(as, bs, zs * q2, zs * s_a, m) * zs /
                   (1.0 - x_a * y);
  const double p = y * s_a * ll_radiator(zp, q2) *
                   dsigma_el_dq2(ap, bp, q2 / zp, s_a, m) / zp;
  return {s, p};
}
/// ... and of the QUASI-ELASTIC tail, Z protons + N neutrons, per nucleus.
std::pair<double, double> ll_peaks_qe(int z, int n, double x, double y,
                                      double s) {
  const double q2 = x * y * s;
  const double zs = (1.0 - y) / (1.0 - x * y);
  const double zp = 1.0 - y * (1.0 - x);
  double sp = 0.0, pp = 0.0;
  for (int k = 0; k < 2; ++k) {
    const bool proton = (k == 0);
    const double mult = proton ? z : n;
    double as, bs, ap, bp;
    rosenbluth_nucleon(proton, zs * q2, &as, &bs);
    rosenbluth_nucleon(proton, q2 / zp, &ap, &bp);
    sp += mult * y * s * ll_radiator(zs, q2) *
          dsigma_el_dq2(as, bs, zs * q2, zs * s, PROTON_MASS) * zs /
          (1.0 - x * y);
    pp += mult * y * s * ll_radiator(zp, q2) *
          dsigma_el_dq2(ap, bp, q2 / zp, s, PROTON_MASS) / zp;
  }
  return {sp, pp};
}

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
      double a, b;
      rosenbluth_nucleon(true, q2, &a, &b);
      const double mott = ALPHA_EM * ALPHA_EM * (1.0 - u) /
                          (4.0 * e * e * u * u) * (ep / e);
      const double lab = mott * (a + b * u / (1.0 - u)) * kPi / (ep * ep);
      CHECK_CLOSE(dsigma_el_dq2(a, b, q2, 2.0 * m * e, m), lab, 1e-10);
    }

    // (i) 6Li at EIC config 1, Q^2 >= 20 GeV^2: POLRAD's claim CONFIRMED --
    //     the t-peak is > 99 % of t + s + p, for BOTH tails.  The nuclear
    //     charge form factor is dead by t ~ 0.25 GeV^2 and the nucleon dipole
    //     by t ~ 1, while the s-/p-peaks sit at t ~ z Q^2 ~ Q^2.
    {
      const std::shared_ptr<HoSpin1FF> f6 = HoSpin1FF::for_ion(LI6());
      const double m6 = LI6().mass(), s6n = 3980.0, s6 = 6.0 * s6n;
      struct Row { double x, y; };
      for (const Row& r : {Row{0.01, 0.5}, Row{0.01, 0.9}, Row{0.1, 0.5},
                           Row{0.3, 0.5}}) {
        const double x_a = r.x / 6.0;
        const double tp = polrad_sigma_el_u(*f6, x_a, r.y, s6, m6, 256);
        const std::pair<double, double> sp =
            ll_peaks_spin1(*f6, x_a, r.y, s6, m6);
        const double f_el = tp / (tp + sp.first + sp.second);
        const double qtp = polrad_sigma_qe_u(3, 3, r.x, r.y, s6n, PROTON_MASS,
                                             256, 0.0);
        const std::pair<double, double> qsp =
            ll_peaks_qe(3, 3, r.x, r.y, s6n);
        const double f_qe = qtp / (qtp + qsp.first + qsp.second);
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
      const std::pair<double, double> sp = ll_peaks_qe(1, 0, x, y, s);
      MESSAGE("proton E = 10 GeV, x = 0.3, Q2 = 2: (s+p)/t = "
              << (sp.first + sp.second) / tp);
      CHECK((sp.first + sp.second) / tp > 10.0);
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
      const std::pair<double, double> se = ll_peaks_spin1(*fd, x_a, y, s_a, m);
      const double tq = polrad_sigma_qe_u(1, 1, x, y, sn, PROTON_MASS, 256,
                                          0.0) / 2.0;
      const std::pair<double, double> sq = ll_peaks_qe(1, 1, x, y, sn);
      const double t_tot = te + tq;
      const double all = t_tot + (se.first + se.second) / 4.0 +
                         (sq.first + sq.second) / 2.0;
      MESSAGE("HERMES deuteron x = 0.012, y = 0.85: t-peak fraction "
              << t_tot / all << " (elastic "
              << te / (te + (se.first + se.second) / 4.0)
              << ", quasi-elastic "
              << tq / (tq + (sq.first + sq.second) / 2.0) << ")");
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

TEST_CASE("T9: the Q_N = 0 Rosenbluth limit of Eq. (38)'s integrand" *
          doctest::skip(true)) {
  // The `PolradFull` half of design sec. 5's T9 -- Im^el_2 == A(Q^2),
  // Im^el_1 == B(Q^2)/2, Im^el_{5,6,7,8} == 0 at Q_N = 0 -- is SKIPPED in v0
  // because Eq. (A.4) is never evaluated: `RcTailModel::PolradFull` is not
  // implemented (rc.cpp throws on it) and the t-peak closed forms of
  // Eqs. (37)-(39) have those contractions already folded in.  Kept in the
  // file, with the reason, so that whoever builds the upgrade path finds it.
  //
  // NOTE for that implementer: design sec. 1.4.6's Im^el_6 has a SPURIOUS
  // eta_A.  polrad2t.tex:2706-2707 and POLRAD's own b2, b3, b4 (adgh:4235-4242
  // through Eq. (A.3), Im_6 = (Q_N/6) eps^3 (b2/3 + b3 + b4)) both give
  //   Im^el_6 = (Q_N/24) ( F_m^2 + (4/(1 + eta_A))
  //                        ((eta_A/3) F_q + F_c + eta_A F_m) F_q )
  // with 4/(1+eta_A), NOT 4 eta_A/(1+eta_A).  Im_2's 4 eta_A^2/(1+eta_A),
  // Im_5, Im_7 and Im_8 all check out (polrad_transcription_check.md sec. 7.1).
  //
  // The v0 half of T9 -- that Eq. (38)'s sigma_u^d integrand IS
  // A(Q^2) Xt - (2/3)(1+eta) F_m^2 -- is covered by T8(a), which compares the
  // code's integrand against that expression written literally in this file.
  CHECK(false);
}

// ------------------------------------------------------------ T1, T4, T17

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
  CHECK(m.exclusion_reason().find("x_L") != std::string::npos);
  CHECK(m.exclusion_reason().find("ASSUMPTION") != std::string::npos);

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

  SUBCASE("PolradFull is refused, with the design section that owns it") {
    RcOptions o;
    o.tail_model = RcTailModel::PolradFull;
    CHECK_THROWS_AS(RcModel(RcMode::TensorBand, o, rc_sampler_ptr(), plan,
                            Channel::Inclusive, LI6()),
                    std::runtime_error);
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

// -------------------------------------------------------- the odds and ends

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
