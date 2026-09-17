// SPDX-License-Identifier: GPL-3.0-or-later
// b1 of a two-cluster nucleus by convolution (b1_nuclear.hpp): the A = 2
// validation gate of design D section 5, the analytic limits of section 6,
// and the number placeholders of section 7.
//
// Design: docs/open_items/run_2026-09-02/design_D_b1_li6.md
// Measured results:  docs/open_items/run_2026-09-03/phase_A_numbers.md
//                      (the CURRENT gate: G3b passes at 0.843243 with CDKS's
//                       own MSTW2008 LO at their Eq. (21); 1.000338 with the
//                       CD-Bonn wave function as well) and phase_A_cdbonn.md
//                    docs/open_items/run_2026-09-02/phase_D_numbers.md
//                    docs/open_items/run_2026-09-02/phase_D_gate.md
//                      (SUPERSEDED on G3b -- it recorded the ToyF2 row, 0.440,
//                       as the verdict; kept for the argument, not the number)
//
// TWO DELIBERATE DEPARTURES FROM THIS REPO'S TEST CONVENTIONS, stated here so
// that a reviewer does not "fix" them.
//  (i) Every other data-dependent test SKIPS when data/vmc is absent
//      (test_cluster.cpp:180).  T1 must NOT: a blocking gate that silently
//      skips is not a blocking gate, so it uses REQUIRE(have(...)) and the
//      suite goes RED when `data/vmc/deuteron/fdeut.av18` is deleted.
// (ii) The "no RNG in the kernel" grep is not expressible in the doctest
//      binary (only LIPOLGEN_REFERENCE_DIR is compiled in, not the source
//      directory); it lives in python/tests/test_b1_model.py as P4.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <stdexcept>
#include <memory>
#include <string>
#include <thread>
#include <vector>

#include "check_close.hpp"
#include "doctest.h"
#include "json_min.hpp"
#include "lipolgen/asymmetries.hpp"
#include "lipolgen/b1_nuclear.hpp"
#include "lipolgen/beams.hpp"
#include "lipolgen/cluster.hpp"
#include "lipolgen/constants.hpp"
#include "lipolgen/numerics.hpp"
#include "lipolgen/pipeline.hpp"
#include "lipolgen/sf.hpp"
#include "lipolgen/spectator.hpp"
#include "lipolgen/tagged.hpp"

// Checklist item 4 (the nucleon PDF) needs a real global fit, so it needs an
// OPTIONAL tier.  MSTW2008 LO rides the PYTHIA tier (`MstwSF` wraps
// `Pythia8::MSTWpdf` over the grid PYTHIA ships), CT18NLO the LHAPDF one.
#ifdef LIPOLGEN_HAVE_PYTHIA8
#include <filesystem>

#include "lipolgen/mstw_sf.hpp"
#endif
#ifdef LIPOLGEN_HAVE_LHAPDF
#include "lipolgen/lhapdf_sf.hpp"
#endif

#ifndef LIPOLGEN_REFERENCE_DIR
#define LIPOLGEN_REFERENCE_DIR "validation/reference"
#endif

using namespace lipolgen;

namespace {

const std::string kFdeut = data_path("vmc/deuteron/fdeut.av18");
const std::string kLi6Momentum = data_path("vmc/momenta/li6_ad1.momentum");
const std::string kLi6Overlap = data_path("vmc/li6_alpha_d/li6.ad");
const std::string kLi7Momentum = data_path("vmc/momenta/li7_at3.momentum");

bool have(const std::string& path) {
  std::ifstream f(path);
  return static_cast<bool>(f);
}

#ifdef LIPOLGEN_HAVE_PYTHIA8
/// Is the MSTW2008 LO grid the gate's VERDICT row needs on disk?
///
/// Called at doctest REGISTRATION time by the `doctest::skip()` decorator of
/// T1v, so that a build without the grid reports that case as SKIPPED in the
/// tally.  It used to be an `if` inside T1 whose else-branch printed a
/// MESSAGE: with the grid absent the case still reported "1 passed" and the
/// clause the 6Li ban lift rests on was never evaluated, while every doc said
/// PASSED unconditionally.  A verdict that can silently not be measured is
/// not a verdict.
bool mstw_grid_present() {
  return std::filesystem::exists(
      std::filesystem::path(pythia8_pdfdata_dir()) / "mstw2008lo.00.dat");
}
#endif

constexpr double kQ2 = 2.5;

/// The Clebsch-Gordan table of design 1.4, <2 m_L : 1 m_S | 1 H>, as EXACT
/// rationals under a square root (sign, numerator, denominator).
struct Cg { int h, ml, ms, sign, num, den; };
const Cg kCg[9] = {
    {+1, 0, +1, +1, 1, 10}, {+1, 1, 0, -1, 3, 10}, {+1, 2, -1, +1, 3, 5},
    {0, -1, +1, +1, 3, 10}, {0, 0, 0, -1, 2, 5},   {0, 1, -1, +1, 3, 10},
    {-1, -2, +1, +1, 3, 5}, {-1, -1, 0, -1, 3, 10}, {-1, 0, -1, +1, 1, 10},
};
double cg(const Cg& c) { return c.sign * std::sqrt(double(c.num) / c.den); }

/// A deuteron kernel built once: the densities cost ~0.27 s each.
///
/// THIS IS THE DEFAULT, i.e. `finite_q_delta = true` (CDKS Eq. (21)'s exact
/// delta-function) since the review of 2026-09-03.  That path REBUILDS the
/// densities at every x -- 0.24 s a point against 0.7 ms -- so every scan
/// below evaluates x*b1 exactly once per grid point and shares the result;
/// doing it the obvious way (once for the landmarks, once for the peak, once
/// for G3b) would triple a 73 s test.
const DeuteronConvolutionB1& gate() {
  static const DeuteronConvolutionB1 d;
  return d;
}

/// The same at CDKS Eq. (17)'s kappa = 1, for the checklist-item-0 column and
/// for everything that is about the QUADRATURE rather than the kinematics
/// (where kappa is irrelevant and 300x cheaper).
const DeuteronConvolutionB1& gate_k1() {
  static const DeuteronConvolutionB1 d([] {
    DeuteronConvolutionB1::Options o;
    o.finite_q_delta = false;
    return o;
  }());
  return d;
}

/// How many threads a landmark scan uses.
unsigned gate_threads() {
  return std::min(8u, std::max(1u, std::thread::hardware_concurrency()));
}

/// x*b1 on `xs`, one evaluation per grid point, spread over threads, with
/// ONE `DeuteronConvolutionB1` PER THREAD: the model is deterministic and
/// every result is a pure function of the constructor arguments (T8 pins two
/// independently built objects bit-identical), but a single object is NOT
/// shareable across threads -- its one-deep density cache
/// (`qd_`/`q_x_`/`q_q2_`, src/core/b1_nuclear.cpp) is mutable per-object
/// state that every new x invalidates.  Same values as a sequential scan over
/// the same Options, bit for bit; measured 72 s -> 12 s on 8 threads, which
/// is half of the whole suite's wall clock.
///
/// `opts` carries one entry per thread, and when it carries a BACKEND
/// structure function (checklist item 4) each entry needs its OWN `unpol`:
/// neither `MstwSF` (Pythia8::PDF::xf is non-const and memoizes the last
/// (x, Q2)) nor `LhapdfSF` (LHAPDF's own interpolator caches) is safe to
/// share across threads.  The Options are built in the CALLING thread, so
/// grid loading is serialized too.
std::vector<double> xb1_threaded(
    const std::vector<double>& xs,
    const std::vector<DeuteronConvolutionB1::Options>& opts) {
  const unsigned nt = static_cast<unsigned>(opts.size());
  std::vector<double> v(xs.size());
  std::vector<std::thread> pool;
  pool.reserve(nt);
  for (unsigned t = 0; t < nt; ++t) {
    pool.emplace_back([&, t] {
      const DeuteronConvolutionB1 d(opts[t]);
      for (std::size_t i = t; i < xs.size(); i += nt)
        v[i] = xs[i] * d.b1(xs[i], kQ2, 0.0);
    });
  }
  for (std::thread& th : pool) th.join();
  return v;
}

/// The same at the DEFAULT (finite-|q|) `DeuteronConvolutionB1`, i.e. the
/// library `ToyF2` -- the gate's layer-3 scan.
std::vector<double> xb1_on_gate_threaded(const std::vector<double>& xs) {
  return xb1_threaded(
      xs, std::vector<DeuteronConvolutionB1::Options>(gate_threads()));
}

/// Everything the G3 clauses ask of ONE x*b1 column, on the layer-3 grid:
/// the sign changes inside G3a's (0, 1.0] counting window with their slopes,
/// the peak position inside [0.5, 1.0], G3b's max|x*b1| over [0.10, 0.80] and
/// its ratio to the digitized peak, and G3c's integral on the x >= 0.01
/// sub-grid.  Hoisted out of checklist item 4 so that item 5 (the deuteron
/// wave function) measures its rows with the SAME code and the two cases
/// cannot drift apart.
struct GateRow {
  std::size_t nz = 0;
  double z0 = 0.0, z1 = 0.0;
  int sl0 = 0, sl1 = 0;
  double xpk = 0.0, win = 0.0, ratio = 0.0;
  double x_min = 0.0, xb1_min = 0.0, x_max = 0.0, xb1_max = 0.0;
  double integral_g3c = 0.0;
};

GateRow gate_row(const std::vector<double>& xs, const std::vector<double>& xb1,
                 const B1Landmarks& ref) {
  const auto at = [&](double x) {
    const std::size_t i = static_cast<std::size_t>(
        std::lower_bound(xs.begin(), xs.end(), x - 1e-12) - xs.begin());
    REQUIRE(i < xs.size());
    return xb1[i];
  };
  const B1Landmarks lm = b1_landmarks(at, xs);
  GateRow r;
  for (std::size_t i = 0; i < lm.zeros.size(); ++i) {
    if (lm.zeros[i] <= 0.0 || lm.zeros[i] > 1.0) continue;
    if (r.nz == 0) { r.z0 = lm.zeros[i]; r.sl0 = lm.zero_slope[i]; }
    if (r.nz == 1) { r.z1 = lm.zeros[i]; r.sl1 = lm.zero_slope[i]; }
    ++r.nz;
  }
  double vpk = -1e30;
  for (std::size_t i = 0; i < xs.size(); ++i) {
    if (xs[i] >= 0.5 && xs[i] <= 1.0 && xb1[i] > vpk) {
      vpk = xb1[i];
      r.xpk = xs[i];
    }
    if (xs[i] >= 0.10 && xs[i] <= 0.80)
      r.win = std::max(r.win, std::fabs(xb1[i]));
  }
  r.ratio = r.win / ref.xb1_max;
  std::vector<double> sub;
  for (std::size_t i = 0; i < xs.size(); ++i)
    if (xs[i] >= 0.01) sub.push_back(xs[i]);
  r.integral_g3c = b1_landmarks(at, sub).integral_b1;
  r.x_min = lm.x_min;
  r.xb1_min = lm.xb1_min;
  r.x_max = lm.x_max;
  r.xb1_max = lm.xb1_max;
  return r;
}

/// G3a on one row.  `who` is a std::string, not a const char*: doctest's
/// MESSAGE stringifies a raw char pointer as a pointer.
void check_g3a(const std::string& who, const GateRow& r,
               const B1Landmarks& ref) {
  MESSAGE(who << " G3a: zeros in (0,1] = " << r.z0 << " (" << r.sl0 << "), "
              << r.z1 << " (" << r.sl1 << ") ; peak in [0.5,1] at x = "
              << r.xpk);
  REQUIRE(r.nz == 2);
  CHECK(r.sl0 == -1);
  CHECK(r.sl1 == +1);
  CHECK(std::fabs(r.z0 - ref.zeros[0]) < 0.08);
  CHECK(std::fabs(r.z1 - ref.zeros[1]) < 0.10);
  CHECK(std::fabs(r.xpk - ref.x_max) < 0.10);
}

const Li6ConvolutionB1& li6() {
  static const Li6ConvolutionB1 m;
  return m;
}

}  // namespace

// ===================================================================== T1
// The A = 2 gate, layers 0-3 (design section 5).  REQUIRE, never a skip.

TEST_CASE("b1_nuclear T1 gate layer 0: the analytic identities (G0a-G0g)") {
  // ---- G0a  <P_zz>_{L=2} = sum_{m_L} CG^2 P_zz(m_S) = +1/10 for H = +-1
  //      (design 1.4's table row), and -1/5 for H = 0.  P_zz = +1 for
  //      m_S = +-1 and -2 for m_S = 0 (CDKS Eq. 46 on a pure state).  Exact
  //      rational arithmetic: every CG^2 is n/den with den | 10, so
  //      "units of 1/10" are integers.
  //
  //      What the D-wave WEIGHT 1/10 actually is, is the delta_T RATIO, and
  //      that is checked here too: delta_T<P_zz> = <P_zz>_{H=0} -
  //      (<P_zz>_{H=+1} + <P_zz>_{H=-1})/2 is -3/10 for the D wave against
  //      -3 for the S wave (m_S = H there), and -0.3/-3 = 1/10 EXACTLY.
  {
    long tenths[3] = {0, 0, 0};                  // H = -1, 0, +1
    for (const Cg& c : kCg) {
      const int pzz = (c.ms == 0) ? -2 : 1;
      tenths[c.h + 1] += pzz * (10 / c.den) * c.num;
    }
    CHECK(tenths[2] == 1);                       // H = +1 -> +1/10
    CHECK(tenths[0] == 1);                       // H = -1 -> +1/10
    CHECK(tenths[1] == -2);                      // H =  0 -> -1/5
    const long dt_d = 2 * tenths[1] - (tenths[2] + tenths[0]);   // 2 x, in 1/10
    CHECK(dt_d == -6);                           // delta_T = -3/10
    const long dt_s = 2 * (-2) * 10 - (1 * 10 + 1 * 10);         // 2 x (-3), 1/10
    CHECK(dt_s == -60);
    CHECK(dt_d * 10 == dt_s);                    // the weight is 1/10 EXACTLY
  }
  // ---- G0b  <S_z>_{L=2} = -1/2 (H = +1; -1/2 x H by symmetry).
  {
    long tenths = 0;                             // in units of 1/10
    for (const Cg& c : kCg) {
      if (c.h != 1) continue;
      tenths += c.ms * (10 / c.den) * c.num;
    }
    CHECK(tenths == -5);                         // -1/2 EXACTLY
    // ... which is `ALPHA_D_VECTOR_POLARIZATION`'s 1 - (3/2) P_D.
    CHECK_CLOSE(1.0 + (-0.5) * P_D_LI6 - 0.0, 1.0 - 0.5 * P_D_LI6, 1e-15);
    CHECK_CLOSE(ALPHA_D_VECTOR_POLARIZATION, 1.0 - 1.5 * P_D_LI6, 1e-15);
  }
  // ---- G0c  Eq. (21)'s SD coefficient -3/(4 sqrt2 pi) from the Y-algebra.
  //      delta_T of <2 0 : 1 H|1 H> is C^0_0 - (C^{+1}_0 + C^{-1}_0)/2, and
  //      the interference is 2 phi_0 phi_2 Y_00 Y_20 x that.
  {
    const double c00 = cg(kCg[4]), cp = cg(kCg[0]), cm = cg(kCg[8]);
    const double dt = c00 - 0.5 * (cp + cm);
    CHECK_CLOSE(dt, -0.9486832980505138, 1e-12);
    const double y00 = 1.0 / std::sqrt(4.0 * kPi);
    const double y20 = std::sqrt(5.0 / (16.0 * kPi));   // Y_20 / (3c^2 - 1)
    CHECK_CLOSE(2.0 * y00 * y20 * dt, -3.0 / (4.0 * std::sqrt(2.0) * kPi), 1e-12);
  }
  // ---- G0d  Eq. (21)'s DD coefficient +3/(16 pi).
  //      delta_T sum_{m_L} |C^H_{m_L}|^2 |Y_{2 m_L}|^2
  //        = (3/10)|Y20|^2 + (3/10)|Y21|^2 - (3/5)|Y22|^2 = (3/16pi)(3c^2-1).
  {
    // |Y_{2m}|^2 as a polynomial in c^2, azimuth-averaged:
    //   |Y20|^2 = (5/16pi)(3c^2-1)^2, |Y21|^2 = (15/8pi) c^2(1-c^2),
    //   |Y22|^2 = (15/32pi)(1-c^2)^2.
    for (double c : {0.0, 0.3, 0.5, 0.77, 1.0}) {
      const double c2 = c * c;
      const double y20 = 5.0 / (16.0 * kPi) * (3 * c2 - 1) * (3 * c2 - 1);
      const double y21 = 15.0 / (8.0 * kPi) * c2 * (1 - c2);
      const double y22 = 15.0 / (32.0 * kPi) * (1 - c2) * (1 - c2);
      const double dt = 0.3 * y20 + 0.3 * y21 - 0.6 * y22;
      CHECK_CLOSE(dt, 3.0 / (16.0 * kPi) * (3 * c2 - 1), 1e-12);
    }
  }
  // ---- G0f  the P_zz-weighted S-D CG sum: P_zz(0) C^0_0 - (C^{+1}_0 +
  //      C^{-1}_0)/2 = +0.9486833, so the DROPPED b1_d-sector SD coefficient
  //      is exactly 1/3 of Eq. (21)'s -0.9486833.
  {
    const double c00 = cg(kCg[4]), cp = cg(kCg[0]), cm = cg(kCg[8]);
    const double s = (-2.0) * c00 - 0.5 * (cp + cm);
    CHECK_CLOSE(s, 0.9486832980505138, 1e-12);
    const double dt = c00 - 0.5 * (cp + cm);       // Eq. (21)'s -0.9486833
    CHECK_CLOSE(s, -dt, 1e-12);
    // The dropped b1_d-sector SD term carries the -(1/3) of
    // sum_m p(m) A_m = A_avg - (P_zz/3) b1_d, so its coefficient is
    // -(1/3) s = +(1/3) dt: exactly 1/3 of Eq. (21)'s, with b1_d in place of
    // F1_d and therefore suppressed by b1_d/F1_d ~ 1e-3 on top of that.
    CHECK_CLOSE((-s / 3.0) / dt, 1.0 / 3.0, 1e-12);
  }
  // ---- G0g  k_range against a brute-force scan of |c*| <= 1, at kappa = 1
  //      and kappa = 1.38, for BOTH mass assignments.
  const double md = nuclear_mass(1, 2), ma = nuclear_mass(2, 4);
  const double eps = LI6_ALPHA_TAG().separation_energy;
  for (double kap : {1.0, 1.38}) {
    for (int swap = 0; swap < 3; ++swap) {
      ConvolutionKinematics kin;
      kin.kappa = kap;
      if (swap == 0) { kin.m_struck = md; kin.m_recoil = ma; kin.separation = eps; }
      else if (swap == 1) { kin.m_struck = ma; kin.m_recoil = md; kin.separation = eps; }
      else { kin.m_struck = M_NUCLEON; kin.m_recoil = M_NUCLEON;
             kin.separation = 2.224574e-3; }
      CAPTURE(kap);
      CAPTURE(swap);
      for (double y : {0.05, 0.4, 0.8, 0.95, 1.0, 1.05, 1.2, 1.4}) {
        double lo = 0.0, hi = 0.0;
        CAPTURE(y);
        if (y > kin.y_max()) {           // struck alpha: y_max = 1.2512
          CHECK(!kin.k_range(y, &lo, &hi));
          continue;
        }
        REQUIRE(kin.k_range(y, &lo, &hi));
        // both endpoints sit exactly on |c*| = 1 ...
        CHECK_CLOSE_AT(std::fabs(kin.cos_star(lo, y)), 1.0, 0.0, 1e-9);
        CHECK_CLOSE_AT(std::fabs(kin.cos_star(hi, y)), 1.0, 0.0, 1e-9);
        // ... and a brute-force scan finds no allowed k outside them.
        const std::vector<double> ks = linspace(1e-6, 3.0 * hi, 200001);
        double blo = 1e30, bhi = -1e30;
        for (double k : ks) {
          if (std::fabs(kin.cos_star(k, y)) <= 1.0) {
            blo = std::min(blo, k);
            bhi = std::max(bhi, k);
          }
        }
        const double dk = ks[1] - ks[0];
        CHECK_CLOSE_AT(blo, lo, 0.0, 2.0 * dk);
        CHECK_CLOSE_AT(bhi, hi, 0.0, 2.0 * dk);
      }
      // no real interval above y_max
      double lo = 0.0, hi = 0.0;
      CHECK(!kin.k_range(kin.y_max() + 1e-6, &lo, &hi));
      CHECK(kin.k_range(kin.y_max() - 1e-6, &lo, &hi));
    }
  }
}

TEST_CASE("b1_nuclear T1 gate layer 0: int delta_T f dy = 0 (G0e)") {
  REQUIRE(have(kFdeut));
  REQUIRE(have(kLi6Momentum));
  REQUIRE(have(kLi6Overlap));
  const FdeutTable t = read_fdeut_k(kFdeut);
  const ClusterPartialWave d0 = ClusterPartialWave::from_uw(t.k_gev, t.u, 0);
  const ClusterPartialWave d2 = ClusterPartialWave::from_uw(t.k_gev, t.w, 2);
  const auto ad = li6_alpha_d_partial_waves();
  const double md = nuclear_mass(1, 2), ma = nuclear_mass(2, 4);
  const double eps = LI6_ALPHA_TAG().separation_energy;

  struct Case { const char* tag; ClusterPartialWave p0, p2; ConvolutionKinematics kin; };
  const std::vector<Case> cases = {
      {"deuteron", d0, d2, {M_NUCLEON, M_NUCLEON, t.ebind_gev, 1.0}},
      {"struck d", ad.first, ad.second, {md, ma, eps, 1.0}},
      {"struck a", ad.first, ad.second, {ma, md, eps, 1.0}},
  };
  for (const Case& c : cases) {
    CAPTURE(c.tag);
    for (double y_min : {1e-4, -3.0}) {
      LightConeDensities::Options o;
      o.renormalize = false;
      o.y_min = y_min;
      const LightConeDensities L(c.p0, c.p2, c.kin, o);
      const std::vector<double>& yy = L.y_grid();
      std::vector<double> a(yy.size()), b(yy.size());
      double mx = 0.0;
      for (std::size_t i = 0; i < yy.size(); ++i) {
        const double v = L.delta_t_f(yy[i]);
        a[i] = v;
        b[i] = v / yy[i];
        mx = std::max(mx, std::fabs(v));
      }
      const double i1 = trapezoid(a, yy), i2 = trapezoid(b, yy);
      CAPTURE(y_min);
      MESSAGE("G0e " << std::string(c.tag) << " y_min=" << y_min << ": int dT f dy = " << i1
                     << " (" << i1 / mx << " x max), int dT f dy/y = " << i2
                     << " (" << i2 / mx << " x max)");
      // On the FULL support the identity holds to a few 1e-5 x max for every
      // wave function and any k cutoff, because int dOmega (3c^2-1) = 0 and
      // int dOmega c (3c^2-1) = 0.  On the design's default grid the deuteron
      // residual is 1.0e-4 / 4.6e-4 x max instead, and does NOT improve with
      // refinement: the y = 1e-4 floor truncates the p > sqrt(2) M_N tail,
      // where E = M - eps - p^2/2M is NEGATIVE and y with it.  That is a
      // support truncation, not a quadrature error -- see the y_min comment
      // in the header.
      const double tol = (y_min < 0.0) ? 1e-4 : 1e-3;
      CHECK(std::fabs(i1) < tol * mx);
      CHECK(std::fabs(i2) < tol * mx);
    }
  }
}

TEST_CASE("b1_nuclear T1 gate layer 1: the unpolarized convolution (G1a-G1f)") {
  REQUIRE(have(kFdeut));
  const FdeutTable t = read_fdeut_k(kFdeut);
  // G1a -- int k^2 (u^2 + w^2) dk on the file's own 0.1 fm^-1 grid.
  std::vector<double> tot(t.k_gev.size()), dw(t.k_gev.size());
  for (std::size_t i = 0; i < t.k_gev.size(); ++i) {
    const double k2 = t.k_gev[i] * t.k_gev[i];
    tot[i] = k2 * (t.u[i] * t.u[i] + t.w[i] * t.w[i]);
    dw[i] = k2 * t.w[i] * t.w[i];
  }
  const double n_table = trapezoid(tot, t.k_gev);
  CHECK_CLOSE(n_table, 0.999976, 1e-4);
  // The SPLINE-refined value is a DIFFERENT number, 1.00046, and the design's
  // G1a row conflates the two: 0.999976 is the table's own trapezoid.  Both
  // are recorded; the linear-interpolation value on the same fine grid is
  // 1.0218, which is the 2.2 % bias the spline exists to avoid.
  const ClusterPartialWave p0 = ClusterPartialWave::from_uw(t.k_gev, t.u, 0);
  const ClusterPartialWave p2 = ClusterPartialWave::from_uw(t.k_gev, t.w, 2);
  {
    const std::vector<double> kk = linspace(0.0, t.k_gev.back(), 100001);
    std::vector<double> a, b, z(kk.size()), zl(kk.size());
    p0.eval_sorted(kk, &a);
    p2.eval_sorted(kk, &b);
    for (std::size_t i = 0; i < kk.size(); ++i) {
      z[i] = kk[i] * kk[i] * (a[i] * a[i] + b[i] * b[i]);
      const double ul = np_interp(kk[i], t.k_gev, t.u);
      const double wl = np_interp(kk[i], t.k_gev, t.w);
      zl[i] = kk[i] * kk[i] * (ul * ul + wl * wl);
    }
    const double spline = trapezoid(z, kk), linear = trapezoid(zl, kk);
    MESSAGE("G1a: table trapezoid " << n_table << ", spline-refined " << spline
                                    << ", LINEAR-refined " << linear);
    CHECK_CLOSE(spline, 1.00046, 1e-4);
    CHECK(linear > 1.02);            // the bias the spline removes
  }
  // G1b -- the file's own D-state probability.
  CHECK_CLOSE(trapezoid(dw, t.k_gev) / n_table, 0.0575999, 1e-5);
  // The header's ebind, and its 2.6e-5 distance from the rounded constant.
  CHECK_CLOSE(t.ebind_gev, 2.224574e-3, 1e-12);
  CHECK_CLOSE(t.ebind_gev, DEUTERON_P_TAG().separation_energy, 3e-5);

  const DeuteronConvolutionB1& d = gate();
  const LightConeDensities& L = d.densities();
  // G1c -- int f dy = 1 - eps/M_N - <p^2>/(2 M_N^2), BOTH sides from the same
  //        table.  Necessarily < 1: the p_z term averages to zero.
  {
    LightConeDensities::Options o;
    o.renormalize = false;
    const LightConeDensities raw(p0, p2, L.kinematics(), o);
    const std::vector<double> kk = linspace(0.0, t.k_gev.back(), 100001);
    std::vector<double> a, b, m2(kk.size()), m0(kk.size());
    p0.eval_sorted(kk, &a);
    p2.eval_sorted(kk, &b);
    for (std::size_t i = 0; i < kk.size(); ++i) {
      const double k2 = kk[i] * kk[i], r = a[i] * a[i] + b[i] * b[i];
      m0[i] = k2 * r;
      m2[i] = k2 * k2 * r;
    }
    const double p2mom = trapezoid(m2, kk) / trapezoid(m0, kk);
    const double ident = 1.0 - t.ebind_gev / M_NUCLEON
                         - p2mom / (2.0 * M_NUCLEON * M_NUCLEON);
    MESSAGE("G1c: raw int f dy = " << raw.norm() << ", moment identity "
                                   << ident << " (design 0.98707)");
    CHECK_CLOSE(raw.norm(), 0.98707, 1e-3);
    CHECK_CLOSE(raw.norm(), 0.9876619, 1e-5);
    CHECK_CLOSE(raw.norm(), ident, 1e-3);
    CHECK(raw.norm() < 1.0);
    // G1d -- mean_y = <y^2>_phi/<y>_phi = [<E^2> + kappa^2 <p^2>/3]/(M <E>).
    // (Design 5.2 writes <k^2>/(3 M^2) inside the bracket; that is
    // dimensionally inconsistent -- the bracket is in GeV^2 -- and the form
    // used here is the one that reproduces its own quoted 0.99546.)
    std::vector<double> e1(kk.size()), e2(kk.size());
    for (std::size_t i = 0; i < kk.size(); ++i) {
      const double r = kk[i] * kk[i] * (a[i] * a[i] + b[i] * b[i]);
      const double e = M_NUCLEON - t.ebind_gev
                       - kk[i] * kk[i] / (2.0 * M_NUCLEON);
      e1[i] = r * e;
      e2[i] = r * e * e;
    }
    const double nrm = trapezoid(m0, kk);
    const double me = trapezoid(e1, kk) / nrm, me2 = trapezoid(e2, kk) / nrm;
    const double my = (me2 + p2mom / 3.0) / (M_NUCLEON * me);
    MESSAGE("G1d: mean_y = " << raw.mean_y() << ", identity " << my
                             << " (design 0.99546)");
    CHECK_CLOSE(raw.mean_y(), 0.99546, 1e-3);
    CHECK_CLOSE(raw.mean_y(), 0.9952583, 1e-5);
    CHECK_CLOSE(raw.mean_y(), my, 2e-3);
  }
  // G1e -- the y support's upper end, kappa = 1 and kappa(x = 0.8, Q2 = 2.5).
  CHECK_CLOSE(L.kinematics().y_max(), 1.4976, 1e-3);
  {
    ConvolutionKinematics kin = L.kinematics();
    kin.kappa = std::sqrt(1.0 + gamma_squared(0.8, kQ2));
    MESSAGE("G1e: y_max(kappa) = " << kin.y_max() << " (design 1.9502)");
    CHECK_CLOSE(kin.y_max(), 1.9502, 1e-3);
  }
  // G1f -- F1D/F1N per nucleon from the SAME convolution: the classic
  //        Fermi-motion/EMC shape (dip ~0.984 near x = 0.5, rise above 0.65).
  //        The design quotes 0.997/0.990/0.986/1.108 at 1 %; the first three
  //        agree, x = 0.8 comes out 1.071 (3.4 % from 1.108) and is pinned to
  //        the measured value -- the large-x point is dominated by the k tail
  //        and by the R choice, which the design's row does not specify.
  {
    const double want[4] = {0.99754, 0.99162, 0.98429, 1.07082};
    const double xs[4] = {0.1, 0.3, 0.5, 0.8};
    for (int i = 0; i < 4; ++i) {
      const double f1d = L.convolve([&L](double y) { return L.f_unpol(y); },
                                    [&](double xp) { return d.f1_nucleon(xp, kQ2); },
                                    xs[i]);
      const double r = f1d / d.f1_nucleon(xs[i], kQ2);
      MESSAGE("G1f x=" << xs[i] << " F1D/F1N = " << r);
      CHECK_CLOSE(r, want[i], 1e-4);
    }
    CHECK(want[2] < want[1]);            // the dip
    CHECK(want[3] > 1.05);               // the rise
  }
}

TEST_CASE("b1_nuclear T1 gate layer 3: CDKS Fig. 4 (G3a hard, G3b/G3c recorded)") {
  REQUIRE(have(kFdeut));
  const B1Landmarks ref = b1_landmarks_of_table();
  // The reference landmarks are COMPUTED from the digitized column, never
  // typed -- a re-digitization moves the target automatically.
  REQUIRE(ref.zeros.size() == 2);
  CHECK(ref.zero_slope[0] == -1);
  CHECK(ref.zero_slope[1] == +1);
  CHECK_CLOSE(ref.zeros[0], 0.06564, 1e-3);
  CHECK_CLOSE(ref.zeros[1], 0.45718, 1e-3);
  CHECK_CLOSE(ref.x_min, 0.3324, 1e-3);
  CHECK_CLOSE(ref.xb1_min, -1.76814e-4, 1e-4);
  CHECK_CLOSE(ref.x_max, 0.7657, 1e-3);
  CHECK_CLOSE(ref.xb1_max, 1.08521e-3, 1e-4);
  CHECK_CLOSE(ref.integral_b1, 4.5920e-4, 1e-4);

  // The scan grid.  Its FLOOR is part of G3a -- a zero below it is invisible
  // to `b1_landmarks` no matter what the counting window says -- so design
  // 5.4 states it (amended 2026-09-03, phase A: 0.01 -> 0.001, a decade below
  // the lowest crossing any PDF tried produces).  Still 300 points, still a
  // plain linspace, so the step moves 0.0052843 -> 0.0053144 and every
  // landmark below shifts by less than that.
  const std::vector<double> xs = linspace(0.001, 1.59, 300);
  // ONE evaluation per grid point, shared by the landmarks, the peak scan and
  // G3b: the default is the finite-|q| delta-function, which rebuilds the
  // light-cone densities at every x (see `gate()`).  300 of those is the
  // single most expensive case in the whole suite, so the scan is threaded --
  // one object per thread, identical values (see `xb1_on_gate_threaded`).
  // `gate()` itself stays the object of the checklist tests below.
  const std::vector<double> xb1 = xb1_on_gate_threaded(xs);
  const auto at = [&](double x) {
    const std::size_t i = static_cast<std::size_t>(
        std::lower_bound(xs.begin(), xs.end(), x - 1e-12) - xs.begin());
    REQUIRE(i < xs.size());
    return xb1[i];
  };
  const B1Landmarks got = b1_landmarks(at, xs);

  // ---- G3a, HARD.  Exactly two sign changes in (0, 1.0], the first
  //      falling within 0.08 of the digitized first zero and the second
  //      rising within 0.10 of the second, and the maximum in [0.5, 1.0]
  //      within 0.10 of the digitized peak position.
  //
  //      The counting window's FLOOR was 0.02 until 2026-09-03; design 5.4
  //      now says (0, 1.0].  0.02 was narrower than the +-0.08 position
  //      tolerance this same clause grants (which reaches down to x = 0), so
  //      it was a second, tighter, unstated position cut and the tolerance
  //      was dead code below it.  (0, 1.0] introduces no new number: it is
  //      0.0656 - 0.08 clipped at the physical floor.  The CEILING stays 1.0
  //      and is load-bearing -- there is a real third crossing at x = 1.2204
  //      here (1.1977 with CT18NLO, 1.2177 with MSTW2008 LO) and "exactly
  //      two" is the only clause that excludes it.
  //
  //      THIS FLOOR IS REGISTRY ROW 18 (`AUTHOR_DECISIONS.md` §B17 ->
  //      `run_2026-09-03/STATUS.md` row 18), APPLIED -- CONFIRM OR REVERT.
  //      The 0.02 -> 0 change was taken on design 5.4's authority and the
  //      author has not ruled on it; until 2026-09-16 no site in the tree
  //      outside the run records pointed at row 18 at all, so the one place a
  //      reader meets the change said only "design 5.4 now says (0, 1.0]".
  std::vector<double> z;
  std::vector<int> sl;
  for (std::size_t i = 0; i < got.zeros.size(); ++i) {
    if (got.zeros[i] > 0.0 && got.zeros[i] <= 1.0) {
      z.push_back(got.zeros[i]);
      sl.push_back(got.zero_slope[i]);
    }
  }
  MESSAGE("G3a: zeros in (0,1] = " << (z.empty() ? 0.0 : z[0]) << ", "
                                   << (z.size() < 2 ? 0.0 : z[1])
                                   << " ; peak at x = " << got.x_max);
  REQUIRE(z.size() == 2);
  CHECK(sl[0] == -1);
  CHECK(sl[1] == +1);
  CHECK(std::fabs(z[0] - ref.zeros[0]) < 0.08);
  CHECK(std::fabs(z[1] - ref.zeros[1]) < 0.10);
  double xpk = 0.0, vpk = -1e30;
  for (std::size_t i = 0; i < xs.size(); ++i) {
    if (xs[i] < 0.5 || xs[i] > 1.0) continue;
    if (xb1[i] > vpk) { vpk = xb1[i]; xpk = xs[i]; }
  }
  CHECK(std::fabs(xpk - ref.x_max) < 0.10);

  // ---- G3b, SOFT but RECORDED.  max|x b1| over [0.10, 0.80] against the
  //      RAW digitized column's maximum.
  //
  //      THIS ROW IS NOT THE GATE'S VERDICT.  It is the LIBRARY ToyF2, and
  //      design 5.4's Escalation clause reads the ratio "after the
  //      checklist" -- checklist item 4 is the nucleon PDF, the dominant
  //      term of the residual.  ToyF2 gives 0.440 (a factor 2.27 low; it was
  //      0.272, a factor 3.68, at the kappa = 1 form this object no longer
  //      defaults to), and that number is pinned here as a REGRESSION guard
  //      on the kernel, not as a verdict.  The verdict is measured with
  //      CDKS's own MSTW2008 LO in the checklist-item-4 case below, where
  //      the ratio is 0.843 -- inside the factor of 2.  See
  //      docs/open_items/run_2026-09-03/phase_A_numbers.md.
  double win = 0.0;
  for (std::size_t i = 0; i < xs.size(); ++i) {
    if (xs[i] < 0.10 || xs[i] > 0.80) continue;
    win = std::max(win, std::fabs(xb1[i]));
  }
  const double ratio = win / ref.xb1_max;
  MESSAGE("G3b (ToyF2, NOT the verdict): max|x b1| over [0.10,0.80] = "
          << win << " ; ratio to the digitized peak = " << ratio
          << " (factor " << 1.0 / ratio << " low)");
  CHECK_CLOSE(win, 4.77494e-4, 2e-3);
  CHECK_CLOSE(ratio, 0.440001, 2e-3);
  // ToyF2 alone is outside G3b's [0.5, 2] window, and the nucleon PDF is why
  // -- pinned as the statement it is, about ToyF2 and not about the gate.
  CHECK(ratio < 0.5);

  // ---- G3c, REPORTED.  Close-Kumano, next to the two library numbers.
  //      Taken on the x >= 0.01 SUB-GRID, per design 5.4's G3c as amended
  //      2026-09-03: `integral_b1` is a trapezoid of b1 = (x*b1)/x, so
  //      running it from the scan floor 0.001 would move it +23 % on 1/x
  //      weighting alone, while the reference `close_kumano_integral(true)`
  //      is over the digitized table's own [0.0100, 1.590].  The
  //      evaluations are REUSED (`at` indexes the full scan), not repeated.
  std::vector<double> xs_g3c;
  for (std::size_t i = 0; i < xs.size(); ++i)
    if (xs[i] >= 0.01) xs_g3c.push_back(xs[i]);
  const B1Landmarks g3c = b1_landmarks(at, xs_g3c);
  MESSAGE("G3c: int b1 dx (this kernel, x >= 0.01 sub-grid, "
          << xs_g3c.front() << "-" << xs_g3c.back() << ") = "
          << g3c.integral_b1 << " ; digitized CDKS "
          << close_kumano_integral(true) << " ; Miller "
          << close_kumano_integral(false));
  CHECK_CLOSE(g3c.integral_b1, 2.14607e-4, 3e-3);
  CHECK_CLOSE(close_kumano_integral(true), ref.integral_b1, 1e-9);
  // The whole-scan integral is NOT G3c and is 24 % larger for the 1/x reason
  // above; pinned so the two are never confused.
  CHECK_CLOSE(got.integral_b1, 2.65316e-4, 3e-3);
  // Landmarks of the computed curve, pinned.
  CHECK_CLOSE(got.x_min, 0.23483, 1e-3);
  CHECK_CLOSE(got.xb1_min, -3.30902e-5, 3e-3);
  CHECK_CLOSE(got.x_max, 0.75564, 1e-3);
  CHECK_CLOSE(got.xb1_max, 4.77494e-4, 2e-3);
  CHECK_CLOSE(z[0], 0.022040, 5e-3);
  CHECK_CLOSE(z[1], 0.377374, 2e-3);
  // G3a's ONE remaining fragile clause, flagged rather than hidden.  The
  // low-x zero's POSITION is the weak landmark: it is 0.0220 here, 0.0098
  // with CT18NLO and 0.0279 with MSTW2008 LO -- a factor 2.8 spread across
  // three nucleon PDFs -- and on a uniform 300-point grid it is bracketed by
  // two points and interpolated, so it is resolved to a few per cent at
  // best.  Its +-0.08 tolerance swallows all of that (the widest miss is
  // 0.056, CT18NLO), which is why the clause is weak but no longer BRITTLE:
  // until 2026-09-03 the counting window's floor of 0.02 sat 0.0021 below
  // this zero and one PDF change turned "exactly two sign changes" into one.
  // The second zero and the peak position are the discriminating landmarks.
  //
  // The second zero clears its +-0.10 window by 0.020 here.  That margin is
  // a property of ToyF2, not of the kernel: MSTW2008 LO puts the second zero
  // at 0.4952, only 0.038 from the digitized 0.4572.
  CHECK(std::fabs(z[1] - ref.zeros[1]) < 0.10);
  CHECK(std::fabs(z[1] - ref.zeros[1]) > 0.05);     // ... by 0.020
}

TEST_CASE("b1_nuclear T1 gate checklist item 0/1/2: kappa, target mass, R") {
  REQUIRE(have(kFdeut));
  const DeuteronConvolutionB1& d = gate();        // Eq. (21), the DEFAULT
  const DeuteronConvolutionB1& d1 = gate_k1();    // Eq. (17), kappa = 1
  // item 0 -- CDKS Eq. (21)'s finite-|q| delta-function against Eq. (17)'s
  // kappa = 1.  1.5-1.7 at large x, the single largest identified effect, and
  // since the review of 2026-09-03 it is the DEFAULT of this object: Eq. (21)
  // is CDKS's exact definition and Eq. (18)'s "~=" is what makes Eq. (17) the
  // approximation to it, so the gate -- which exists to reproduce THEIR
  // figure -- is quoted there.  (`Li6ConvolutionOptions` still defaults to
  // kappa = 1; in 6Li the same switch is worth 1 %.)
  CHECK(d.options().finite_q_delta);
  CHECK(!d1.options().finite_q_delta);
  const double want[3] = {1.5692, 1.5791, 1.6876};
  const double xs[3] = {0.5, 0.7, 0.8};
  for (int i = 0; i < 3; ++i) {
    const double r = d.b1(xs[i], kQ2, 0.0) / d1.b1(xs[i], kQ2, 0.0);
    MESSAGE("checklist 0: kappa/kappa=1 at x=" << xs[i] << " -> " << r);
    CHECK_CLOSE(r, want[i], 3e-3);
    CHECK(r > 1.4);
    CHECK(r < 1.8);
  }
  // item 1 -- CDKS Eq. (22)'s target-mass factor, worth ~1.5 at x = 0.8.
  // Measured at BOTH delta-functions, because the factor is not the same
  // number at the two (1.487 at kappa, 1.520 at kappa = 1).
  DeuteronConvolutionB1::Options om;
  om.target_mass = false;
  const DeuteronConvolutionB1 dm(om);
  const double rm = d.b1(0.8, kQ2, 0.0) / dm.b1(0.8, kQ2, 0.0);
  MESSAGE("checklist 1: target-mass factor at x=0.8 -> " << rm);
  CHECK_CLOSE(rm, 1.48729, 3e-3);
  om.finite_q_delta = false;
  const DeuteronConvolutionB1 dm1(om);
  CHECK_CLOSE(d1.b1(0.8, kQ2, 0.0) / dm1.b1(0.8, kQ2, 0.0), 1.52049, 3e-3);
  // item 2 -- R.  r1998 is the DEFAULT here (CDKS's choice); r_sigma_lt moves
  // the second zero AWAY from the digitized 0.4572.
  //
  // MEASURED AT kappa = 1, deliberately: R is a property of F1 and is
  // independent of the delta-function, while a 300-point landmark scan at the
  // default kappa costs 73 s (the densities are rebuilt at every x).  The
  // r1998 column it is compared against is the kappa = 1 one, so the two are
  // like for like.
  DeuteronConvolutionB1::Options orr;
  orr.finite_q_delta = false;
  orr.r_func = [](double x, double q2) { return r_sigma_lt(x, q2); };
  const DeuteronConvolutionB1 dr(orr);
  // Same grid as the layer-3 case above, floor included -- changing one and
  // not the other would split the two cases and make phase_D_gate.md's
  // kappa column no longer like for like.
  const std::vector<double> xg = linspace(0.001, 1.59, 300);
  const B1Landmarks lr =
      b1_landmarks([&](double x) { return x * dr.b1(x, kQ2, 0.0); }, xg);
  const B1Landmarks l1 =
      b1_landmarks([&](double x) { return x * d1.b1(x, kQ2, 0.0); }, xg);
  REQUIRE(lr.zeros.size() >= 2);
  REQUIRE(l1.zeros.size() >= 2);
  MESSAGE("checklist 2 (kappa = 1): second zero with r_sigma_lt = "
          << lr.zeros[1] << " against " << l1.zeros[1]
          << " with r1998 (digitized 0.45718); peak " << lr.xb1_max << " at "
          << lr.x_max << " against " << l1.xb1_max << " at " << l1.x_max);
  CHECK_CLOSE(lr.zeros[1], 0.365202, 3e-3);
  CHECK_CLOSE(l1.zeros[1], 0.392320, 3e-3);
  CHECK_CLOSE(lr.xb1_max, 2.99431e-4, 3e-3);
  CHECK_CLOSE(l1.xb1_max, 2.94839e-4, 3e-3);
  // ... and the kappa = 1 column of the whole gate, pinned next to the
  // default's, so phase_D_gate.md's two rows are both regression-guarded.
  MESSAGE("checklist 0 (kappa = 1 column): zeros " << l1.zeros[0] << ", "
          << l1.zeros[1] << " ; peak " << l1.xb1_max << " at " << l1.x_max
          << " ; int b1 dx " << l1.integral_b1);
  CHECK_CLOSE(l1.zeros[0], 0.022060, 5e-3);
  CHECK_CLOSE(l1.x_max, 0.73970, 1e-3);
  // G3c on the x >= 0.01 sub-grid, exactly as in the layer-3 case: the
  // whole-scan value is 47 % larger and is the 1/x weight of the extension
  // to 0.001, not physics.  Both are pinned so the two never get confused.
  std::vector<double> xg_g3c;
  for (std::size_t i = 0; i < xg.size(); ++i)
    if (xg[i] >= 0.01) xg_g3c.push_back(xg[i]);
  const B1Landmarks l1c =
      b1_landmarks([&](double x) { return x * d1.b1(x, kQ2, 0.0); }, xg_g3c);
  CHECK_CLOSE(l1c.integral_b1, 1.07152e-4, 3e-3);
  CHECK_CLOSE(l1.integral_b1, 1.57861e-4, 3e-3);
  CHECK_CLOSE(l1.xb1_max / b1_landmarks_of_table().xb1_max, 0.271689, 2e-3);
}

// Checklist item 4, THE VERDICT ROW: MSTW2008 LO, CDKS's OWN nucleon PDF, at
// CDKS's OWN Eq. (21) delta-function.  This is the row design 5.4's
// Escalation clause is read off and the row the 6Li publication ban was
// lifted on, so it is the one clause in this file that must never be able to
// not-run quietly.
//
// It lives in its own case, decorated with `doctest::skip(...)`, for exactly
// that reason.  Until 2026-09-04 it was an `if (!exists(grid)) { MESSAGE(); }
// else { ... }` inside the item-4 case: with `mstw2008lo.00.dat` absent that
// case reported "1 passed", the verdict was never measured, and every doc in
// the tree said "PASSED" in that build too.  As a decorated case, doctest's
// own tally reports it SKIPPED by name -- the build states what it could not
// evaluate, in the one place a reader always looks.
//
// So: layers 0-3 REQUIRE their data (a blocking gate that silently skips is
// not a blocking gate); this row SKIPS, visibly, because it additionally
// needs an OPTIONAL tier.  The two are not in tension.
#ifdef LIPOLGEN_HAVE_PYTHIA8
TEST_CASE("b1_nuclear T1v gate checklist item 4, THE VERDICT ROW: MSTW2008 LO"
          " at CDKS Eq. (21) -- G3a/G3b/G3c" *
          doctest::skip(!mstw_grid_present())) {
  REQUIRE(have(kFdeut));
  REQUIRE(mstw_grid_present());   // the decorator already guaranteed it
  const B1Landmarks ref = b1_landmarks_of_table();
  const std::vector<double> xs = linspace(0.001, 1.59, 300);   // as layer 3
  const auto measure = [&](const std::vector<double>& xb1) {
    return gate_row(xs, xb1, ref);
  };
  const auto g3a = [&](const std::string& who, const GateRow& r) {
    check_g3a("item 4 " + who, r, ref);
  };

  // One MstwSF PER THREAD: Pythia8::PDF::xf is non-const and memoizes the
  // last (x, Q2), so a shared instance is a data race.  8 MB of grid each
  // (sizeof(Pythia8::MSTWpdf) == 7 988 336), built here in the calling
  // thread so the file reads are serialized.
  const unsigned nt = gate_threads();
  std::vector<DeuteronConvolutionB1::Options> opts(nt);
  for (unsigned t = 0; t < nt; ++t) opts[t].unpol = std::make_shared<MstwSF>();
  const GateRow m = measure(xb1_threaded(xs, opts));
  g3a("MSTW2008 LO, kappa", m);
  MESSAGE("item 4 MSTW2008 LO, kappa  G3b: max|x b1| over [0.10,0.80] = "
          << m.win << " ; ratio to the digitized peak = " << m.ratio
          << " -- G3b asks for [0.5, 2.0], so this PASSES (factor "
          << 1.0 / m.ratio << " low)");
  // G3b, THE GATE'S MAGNITUDE VERDICT, at CDKS's own PDF and CDKS's own
  // Eq. (21) delta-function.  0.843 is inside the factor of 2 with room:
  // the ToyF2 row of the layer-3 case is 0.440 and the CT18NLO stand-in
  // 0.719, so the nucleon PDF is worth x1.92 and x1.17 respectively.
  CHECK_CLOSE(m.win, 9.15096e-4, 2e-3);
  CHECK_CLOSE(m.ratio, 0.843243, 2e-3);
  CHECK(m.ratio > 0.5);
  CHECK(m.ratio < 2.0);
  // G3a's landmarks, pinned.  Note that CDKS's own PDF reproduces CDKS's
  // own figure BETTER than either stand-in on every one of them: the
  // second zero misses by 0.038 (ToyF2 0.080, CT18NLO 0.018) and the peak
  // position by 0.0059 (ToyF2 0.0101, CT18NLO 0.0313).
  CHECK_CLOSE(m.z0, 0.027886, 5e-3);
  CHECK_CLOSE(m.z1, 0.495239, 2e-3);
  CHECK_CLOSE(m.x_min, 0.35706, 1e-3);
  CHECK_CLOSE(m.xb1_min, -2.34606e-4, 3e-3);
  CHECK_CLOSE(m.x_max, 0.77159, 1e-3);
  CHECK_CLOSE(m.xb1_max, 9.15096e-4, 2e-3);
  // G3c, x >= 0.01 sub-grid, against the digitized +4.592e-4.
  MESSAGE("item 4 MSTW2008 LO, kappa  G3c: int b1 dx (x >= 0.01) = "
          << m.integral_g3c << " ; digitized CDKS "
          << close_kumano_integral(true));
  CHECK_CLOSE(m.integral_g3c, 2.24896e-4, 3e-3);

  // --- the same at Eq. (17)'s kappa = 1, the checklist-item-0 column.
  //     300x cheaper, so no threads.  It lands at 0.520 -- INSIDE G3b's
  //     window but by only 4 %, which is the honest measure of how much of
  //     the gate's pass the finite-|q| delta-function is carrying.
  DeuteronConvolutionB1::Options o1;
  o1.unpol = std::make_shared<MstwSF>();
  o1.finite_q_delta = false;
  const DeuteronConvolutionB1 d1m(o1);
  std::vector<double> v1(xs.size());
  for (std::size_t i = 0; i < xs.size(); ++i)
    v1[i] = xs[i] * d1m.b1(xs[i], kQ2, 0.0);
  const GateRow m1 = measure(v1);
  g3a("MSTW2008 LO, kappa = 1", m1);
  MESSAGE("item 4 MSTW2008 LO, kappa = 1  G3b ratio = " << m1.ratio);
  CHECK_CLOSE(m1.ratio, 0.520332, 2e-3);
  CHECK(m1.ratio > 0.5);
  CHECK_CLOSE(m1.z0, 0.027946, 5e-3);
  CHECK_CLOSE(m1.z1, 0.507853, 2e-3);
  CHECK_CLOSE(m1.xb1_max, 5.64669e-4, 2e-3);
  CHECK_CLOSE(m1.integral_g3c, 5.79583e-5, 3e-3);
}
#endif  // LIPOLGEN_HAVE_PYTHIA8

// The one case that runs in EVERY build and says whether the tree's
// "the A = 2 magnitude gate PASSED" claim was actually evaluated here.
//
// T1v's `doctest::skip` shows up only as the skipped COUNT, which is easy to
// read past.  This prints the sentence, so `build/lipolgen_tests` itself
// answers the question the docs now promise it answers (docs/USAGE.md sec.
// 2a, condition 1: "the run prints SKIPPED (MSTW2008 LO unavailable)").
TEST_CASE("b1_nuclear T1r: was the G3b VERDICT row measurable in this build?") {
  // Not vacuous: it pins the selector name the docs tell a user to set in
  // order to GENERATE in the configuration the verdict was read off.  A
  // verdict row that a build cannot measure and a run surface that cannot
  // emit it are the same defect seen from two sides.
  CHECK(std::string(b1_unpol_name(B1UnpolSource::Mstw)) == "mstw");
#ifdef LIPOLGEN_HAVE_PYTHIA8
  if (mstw_grid_present()) {
    MESSAGE("VERDICT ROW MEASURED: MSTW2008 LO is on disk under "
            << pythia8_pdfdata_dir()
            << ", so T1v ran and G3b's 0.843243 was checked in THIS build.");
  } else {
    MESSAGE("SKIPPED (MSTW2008 LO unavailable): no mstw2008lo.00.dat under "
            << pythia8_pdfdata_dir()
            << " -- T1v did NOT run, so G3b's verdict row (0.843243, the "
               "clause the 6Li publication ban was lifted on) was NOT "
               "evaluated in this build.  Every 'the gate passes' statement "
               "in the tree is a statement about the MSTW2008 LO "
               "configuration, which this build cannot reproduce and, by "
               "PipelineConfig::validate(), also cannot emit.");
  }
#else
  MESSAGE("SKIPPED (MSTW2008 LO unavailable): built without the PYTHIA tier "
          "-- T1v is not even compiled in, so G3b's verdict row (0.843243) "
          "was NOT evaluated in this build.  This build's only measured G3b "
          "is the ToyF2 row, 0.440001, which is OUTSIDE the gate's [0.5, 2] "
          "window; it also cannot emit --b1-unpol mstw.");
#endif
}

// Checklist item 4 of design 5.4 -- the NUCLEON PDF, CT18NLO half.
//
// CDKS built the Fig. 4 curve this gate compares against on MSTW2008 LO, so
// the MSTW row is the LIKE-FOR-LIKE measurement and the one G3b's verdict is
// read off -- it is T1v above, in its own skip-decorated case.  The CT18NLO
// row is here beside it -- not as the answer but so the PDF DEPENDENCE
// itself is pinned; `MstwSF` and `LhapdfSF` share the charge weights and
// `f2_from_weights` VERBATIM (mstw_sf.hpp, "COMPARABILITY"), so the two rows
// differ in the grid and in nothing else.
//
// Compiled away, not skipped, when LHAPDF is absent: unlike the verdict row
// there is no claim in the tree that rests on CT18NLO, so a build without
// the tier has nothing to report about it.
#ifdef LIPOLGEN_HAVE_LHAPDF
TEST_CASE("b1_nuclear T1 gate checklist item 4: the nucleon PDF (G3a/G3b/G3c)") {
  REQUIRE(have(kFdeut));
  const B1Landmarks ref = b1_landmarks_of_table();
  const std::vector<double> xs = linspace(0.001, 1.59, 300);   // as layer 3

  // Everything G3 asks of one x*b1 column, on the layer-3 grid and with the
  // layer-3 clauses: the sign changes inside G3a's (0, 1.0] counting window,
  // the peak position inside [0.5, 1.0], G3b's max|x*b1| over [0.10, 0.80],
  // and G3c's integral on the x >= 0.01 sub-grid.
  const auto measure = [&](const std::vector<double>& xb1) {
    return gate_row(xs, xb1, ref);
  };
  const auto g3a = [&](const std::string& who, const GateRow& r) {
    check_g3a("item 4 " + who, r, ref);
  };

  // --- CT18NLO: the modern cross-check, kept so the PDF DEPENDENCE is
  //     pinned and not only MSTW's answer.  It is also the row that made the
  //     G3a counting window a live question: its low-x zero is at 0.0098,
  //     below BOTH the old counting floor of 0.02 and the old scan floor of
  //     0.01, so under the pre-2026-09-03 clause this REQUIRE(nz == 2) would
  //     abort.  Under the amended clause it passes, by 0.024 of its +-0.08.
  lhapdf_quiet();
  {
    const unsigned nt = gate_threads();
    std::vector<DeuteronConvolutionB1::Options> opts(nt);
    for (unsigned t = 0; t < nt; ++t)
      opts[t].unpol = std::make_shared<LhapdfSF>("CT18NLO", 0);
    const GateRow c = measure(xb1_threaded(xs, opts));
    g3a("CT18NLO, kappa", c);
    MESSAGE("item 4 CT18NLO, kappa  G3b ratio = " << c.ratio
            << " ; MSTW2008 LO is the like-for-like row");
    CHECK_CLOSE(c.win, 7.80734e-4, 2e-3);
    CHECK_CLOSE(c.ratio, 0.719432, 2e-3);
    CHECK_CLOSE(c.z0, 0.009795, 5e-3);
    CHECK_CLOSE(c.z1, 0.438761, 2e-3);
    CHECK_CLOSE(c.x_min, 0.29861, 1e-3);
    CHECK_CLOSE(c.xb1_min, -1.86774e-4, 3e-3);
    CHECK_CLOSE(c.x_max, 0.73438, 1e-3);
    CHECK_CLOSE(c.xb1_max, 7.80734e-4, 2e-3);
    CHECK_CLOSE(c.integral_g3c, 2.09620e-4, 3e-3);
    // The old floor, stated as the measurement that motivated moving it.
    CHECK(c.z0 < 0.01);
    CHECK(std::fabs(c.z0 - ref.zeros[0]) < 0.08);
  }
}
#endif  // LIPOLGEN_HAVE_LHAPDF

// Checklist item 5 of design 5.4 -- the DEUTERON WAVE FUNCTION, i.e. gate
// condition 3 of open item 10: "a real CD-Bonn u, w instead of the D-state
// rescaling proxy of item 5".  CDKS built the Fig. 4 curve this gate is
// compared against on CD-Bonn; every gate number before this case was
// measured on AV18.
//
// THE OPTION IS OPT-IN AND THE DEFAULT DOES NOT MOVE.
// `DeuteronConvolutionB1::Options::wave` defaults to `kFdeutFile`, so layer 3
// and item 4 above are untouched and every published gate number keeps its
// meaning.  What this case does is MEASURE what the swap is worth, on the
// same grid and with the same clauses.
//
// THE RESULT, and it is the largest single move the gate has seen: with
// CD-Bonn AND MSTW2008 LO -- CDKS's own wave function and CDKS's own PDF --
// the kernel reproduces the digitized CDKS Fig. 4 on EVERY G3 landmark:
//
//   landmark              digitized CDKS   AV18 (item 4)   CD-Bonn (here)
//   G3b peak ratio            1              0.8432          1.0003
//   peak position x           0.765663       0.771585        0.766271
//   first zero                0.065645       0.027886        0.064129
//   second zero               0.457177       0.495239        0.457018
//   dip                      -1.76814e-4    -2.34606e-4     -1.76907e-4
//   dip position              0.33236        0.35706         0.33049
//   G3c int b1 dx (x>=0.01)   4.59200e-4     2.24896e-4      4.48580e-4
//
// Nothing was tuned: the coefficients are Machleidt's Table XX (typed once,
// in cluster.hpp, and gated against his own published deuteron properties in
// tests/test_cluster.cpp), the PDF is PYTHIA's shipped MSTW grid, and the
// kernel is unchanged.  The mechanism is in
// docs/open_items/run_2026-09-03/phase_A_cdbonn.md section 7: CD-Bonn's S
// node sits 13 % higher in k than AV18's and its u is half as big beyond it,
// so the NEGATIVE lobe of the S-D interference that cancels part of the
// positive one in AV18 is almost entirely gone.
//
// WHAT THIS DOES NOT SAY.  It does not make `Li6ConvolutionB1` validated --
// the A = 2 gate is a necessary condition and 6Li adds the alpha-d wave
// function, the four-term truncation and N_ad on top.  It does not remove the
// digitization uncertainty of the reference column.  And a 0.03 % agreement
// on a curve read off a published figure is BETTER THAN THE REFERENCE
// DESERVES; read it as "the residual is now below the digitization error",
// not as three-digit agreement with CDKS.
TEST_CASE("b1_nuclear T1 gate checklist item 5: the deuteron wave function") {
  REQUIRE(have(kFdeut));
  const B1Landmarks ref = b1_landmarks_of_table();
  const std::vector<double> xs = linspace(0.001, 1.59, 300);   // as layer 3

  // ---- THE SIGN GATE, and it is the cheapest one available.  CD-Bonn is
  //      published with a BARE j_L Fourier kernel for both L (his Eq. D13),
  //      which is inconsistent at L = 2 with his own Eqs. (D20)/(D22); the
  //      drop-in convention is w = -psi_2^a (cluster.hpp, `CdBonnWave`).  Get
  //      it wrong and the S-D interference flips.  `alpha_d_quadrupole_fm2`
  //      settles it against CD-Bonn's own published Q_d = 0.270 fm^2, from
  //      the same code path that returns +0.2694 on the AV18 file (whose own
  //      header prints qm = 0.269673).
  const FdeutTable cdb = cdbonn_fdeut_table();
  const ClusterPartialWave c0 = ClusterPartialWave::from_uw(cdb.k_gev, cdb.u, 0);
  const ClusterPartialWave c2 = ClusterPartialWave::from_uw(cdb.k_gev, cdb.w, 2);
  std::vector<double> w_flipped = cdb.w;
  for (double& v : w_flipped) v = -v;
  const ClusterPartialWave c2_wrong =
      ClusterPartialWave::from_uw(cdb.k_gev, w_flipped, 2);
  const FdeutTable av = read_fdeut_k(kFdeut);
  const ClusterPartialWave a0 = ClusterPartialWave::from_uw(av.k_gev, av.u, 0);
  const ClusterPartialWave a2 = ClusterPartialWave::from_uw(av.k_gev, av.w, 2);
  const double qd_cdb = alpha_d_quadrupole_fm2(c0, c2);
  const double qd_bad = alpha_d_quadrupole_fm2(c0, c2_wrong);
  const double qd_av = alpha_d_quadrupole_fm2(a0, a2);
  MESSAGE("item 5 Q_d: CD-Bonn " << qd_cdb << " fm^2 (published 0.270), "
          << "wrong sign " << qd_bad << ", AV18 " << qd_av
          << " fm^2 (fdeut.av18 header qm = 0.269673)");
  CHECK(qd_cdb > 0.0);
  CHECK(std::fabs(qd_cdb - 0.270) < 5e-4);
  CHECK_CLOSE(qd_cdb, 0.270178, 2e-5);
  CHECK(qd_bad < 0.0);                       // the wrong sign is unmistakable
  CHECK_CLOSE(qd_av, 0.269362, 2e-5);

  // ---- The two wave functions ON THE SAME GRID, so that the gate rows below
  //      have a stated cause.  CD-Bonn is a SOFTER wave function overall
  //      (P_D 4.86 % against AV18's 5.76 %) but its D wave is LARGER at small
  //      k and smaller only in the tail -- which is exactly the structure the
  //      item-5 rescaling proxy (one factor multiplying w everywhere) cannot
  //      have, and why the proxy moved the gate the other way.
  double nd_c = 0.0, nd_a = 0.0;
  for (std::size_t i = 1; i < av.k_gev.size(); ++i) {
    const double dk = av.k_gev[i] - av.k_gev[i - 1];
    nd_c += 0.5 * dk * (cdb.k_gev[i - 1] * cdb.k_gev[i - 1] * cdb.w[i - 1]
                            * cdb.w[i - 1]
                        + cdb.k_gev[i] * cdb.k_gev[i] * cdb.w[i] * cdb.w[i]);
    nd_a += 0.5 * dk * (av.k_gev[i - 1] * av.k_gev[i - 1] * av.w[i - 1]
                            * av.w[i - 1]
                        + av.k_gev[i] * av.k_gev[i] * av.w[i] * av.w[i]);
  }
  MESSAGE("item 5 P_D on the shared 0.1 fm^-1 grid: CD-Bonn " << nd_c
          << " (published 4.85 %), AV18 " << nd_a
          << " (fdeut.av18 header dstate = 0.057599)");
  CHECK_CLOSE(nd_c, 0.0485621, 1e-5);
  CHECK_CLOSE(nd_a, 0.0575985, 1e-5);
  // w(CD-Bonn)/w(AV18) runs from 1.02 at 0.1 fm^-1 through 1 near 0.8 to 0.36
  // at 5 fm^-1: a rescaling proxy is a CONSTANT here and cannot be either.
  const auto ratio_at = [&](double p_fm) {
    const std::size_t i = static_cast<std::size_t>(std::llround(p_fm / 0.1));
    return cdb.w[i] / av.w[i];
  };
  MESSAGE("item 5 w(CD-Bonn)/w(AV18) at p = 0.1, 1.0, 5.0 fm^-1 = "
          << ratio_at(0.1) << ", " << ratio_at(1.0) << ", " << ratio_at(5.0));
  CHECK_CLOSE(ratio_at(0.1), 1.02026, 1e-4);
  CHECK_CLOSE(ratio_at(1.0), 0.988783, 1e-4);
  CHECK_CLOSE(ratio_at(5.0), 0.3597447, 1e-4);
  // The S node, the feature that drives the b1 result.
  CHECK(av.u[20] > 0.0);            // AV18's node is above p = 2.0 fm^-1 ...
  CHECK(av.u[21] < 0.0);            // ... and below 2.1
  CHECK(cdb.u[23] > 0.0);           // CD-Bonn's is above 2.3 ...
  CHECK(cdb.u[24] < 0.0);           // ... and below 2.4, i.e. 13 % higher

  const auto cdbonn_opts = [&](unsigned nt) {
    std::vector<DeuteronConvolutionB1::Options> o(nt);
    for (unsigned t = 0; t < nt; ++t) o[t].wave = DeuteronWaveSource::kCdBonn;
    return o;
  };

  // ---- ToyF2, the row that runs with NO optional tier, so this case is
  //      never vacuous.  Against layer 3's AV18 value 0.440001, CD-Bonn is
  //      0.528662 -- and its curve has NO sign change in (0, 1].  That is a
  //      property of ToyF2, not of the wave function: both real PDFs below
  //      put two zeros back, close to the digitized ones.  So this row is
  //      RECORDED and G3a's counting clause is deliberately NOT applied to
  //      it; applying it would gate the wave function on the toy.
  {
    const GateRow t = gate_row(xs, xb1_threaded(xs, cdbonn_opts(gate_threads())),
                               ref);
    MESSAGE("item 5 ToyF2, kappa, CD-Bonn: G3b ratio = " << t.ratio
            << " (layer 3's AV18 row is 0.440001), zeros in (0,1] = " << t.nz
            << ", peak " << t.xb1_max << " at x = " << t.x_max);
    CHECK(t.nz == 0);
    CHECK_CLOSE(t.ratio, 0.528662, 2e-3);
    CHECK_CLOSE(t.xb1_max, 5.737089e-4, 2e-3);
    CHECK_CLOSE(t.integral_g3c, 4.499038e-4, 3e-3);
  }

#ifdef LIPOLGEN_HAVE_LHAPDF
  // ---- CT18NLO: the same swap on a DIFFERENT modern PDF, so that "CD-Bonn
  //      lifts the peak" is not an MSTW artefact.  0.719432 -> 0.864233, the
  //      same direction and 80 % of the same size.  This is also the row the
  //      research stage predicted (phase_A_cdbonn.md section 8.2, 0.8642 on
  //      its own grid) BEFORE any of this was written, which is why it is
  //      kept rather than dropped as redundant.
  lhapdf_quiet();
  {
    const unsigned nt = gate_threads();
    std::vector<DeuteronConvolutionB1::Options> o = cdbonn_opts(nt);
    for (unsigned t = 0; t < nt; ++t)
      o[t].unpol = std::make_shared<LhapdfSF>("CT18NLO", 0);
    const GateRow c = gate_row(xs, xb1_threaded(xs, o), ref);
    check_g3a("item 5 CT18NLO, kappa, CD-Bonn", c, ref);
    MESSAGE("item 5 CT18NLO, kappa, CD-Bonn  G3b ratio = " << c.ratio
            << " (AV18 on the same row is 0.719432)");
    CHECK_CLOSE(c.ratio, 0.864233, 2e-3);
    CHECK_CLOSE(c.z0, 0.043538, 5e-3);
    CHECK_CLOSE(c.z1, 0.391923, 2e-3);
    CHECK_CLOSE(c.xb1_max, 9.378742e-4, 2e-3);
    CHECK_CLOSE(c.integral_g3c, 4.389734e-4, 3e-3);
    // CT18NLO's low-x zero was BELOW the old 0.02 counting floor on AV18
    // (item 4 pins it at 0.0098); CD-Bonn lifts it to 0.0435, inside even the
    // pre-2026-09-03 window.  Recorded because the amended window was argued
    // on that AV18 row.
    CHECK(c.z0 > 0.02);
  }
#endif  // LIPOLGEN_HAVE_LHAPDF
}

// ===================================================== item 5v, GATE COND. 3
//
// CD-Bonn AND MSTW2008 LO TOGETHER -- the 1.000338 peak ratio that is the
// THIRD condition of the b1 ban lift (docs/OPEN_ITEMS_SOLUTIONS.md sec. 10
// item 3, STATUS.md decision row 4).  It used to be an `if (!mstw_grid_
// present()) { MESSAGE(...) } else { ... }` INSIDE item 5, and with the grid
// absent doctest tallied that case as PASSED while the condition the lift
// rests on was never evaluated -- the very defect T1v was split out for
// (docs/USAGE.md sec. 4, "the doctest's MSTW rows are skipped loudly").  As
// its own decorated case the build now states in its tally what it could not
// measure.  Split out 2026-09-05 (phase F); every pin below is unchanged.
#ifdef LIPOLGEN_HAVE_PYTHIA8
TEST_CASE("b1_nuclear T1 gate checklist item 5v: CD-Bonn + MSTW2008 LO --"
          " gate condition 3" * doctest::skip(!mstw_grid_present())) {
  REQUIRE(have(kFdeut));
  REQUIRE(mstw_grid_present());   // the decorator already guaranteed it
  const B1Landmarks ref = b1_landmarks_of_table();
  const std::vector<double> xs = linspace(0.001, 1.59, 300);   // as item 5
  const auto cdbonn_opts = [](unsigned nt) {
    std::vector<DeuteronConvolutionB1::Options> o(nt);
    for (unsigned t = 0; t < nt; ++t) o[t].wave = DeuteronWaveSource::kCdBonn;
    return o;
  };
  const unsigned nt = gate_threads();
  std::vector<DeuteronConvolutionB1::Options> o = cdbonn_opts(nt);
  for (unsigned t = 0; t < nt; ++t) o[t].unpol = std::make_shared<MstwSF>();
  const GateRow m = gate_row(xs, xb1_threaded(xs, o), ref);
  check_g3a("item 5 MSTW2008 LO, kappa, CD-Bonn", m, ref);
  MESSAGE("item 5 MSTW2008 LO, kappa, CD-Bonn  G3b: max|x b1| over "
          "[0.10,0.80] = " << m.win << " ; ratio to the digitized peak = "
          << m.ratio << " (AV18 on the same row is 0.843243)");
  CHECK_CLOSE(m.win, 1.085577e-3, 2e-3);
  CHECK_CLOSE(m.ratio, 1.000338, 2e-3);
  CHECK(m.ratio > 0.5);
  CHECK(m.ratio < 2.0);
  // G3a's landmarks.  Each is CLOSER to the digitized column than the AV18
  // row of item 4, by an order of magnitude on the two zeros and the dip.
  CHECK_CLOSE(m.z0, 0.064129, 5e-3);
  CHECK_CLOSE(m.z1, 0.457018, 2e-3);
  CHECK_CLOSE(m.x_min, 0.33049, 1e-3);
  CHECK_CLOSE(m.xb1_min, -1.769065e-4, 3e-3);
  CHECK_CLOSE(m.x_max, 0.766271, 1e-3);
  CHECK(std::fabs(m.z0 - ref.zeros[0]) < 0.005);
  CHECK(std::fabs(m.z1 - ref.zeros[1]) < 0.005);
  CHECK(std::fabs(m.x_max - ref.x_max) < 0.005);
  CHECK(std::fabs(m.xb1_min / ref.xb1_min - 1.0) < 0.01);
  MESSAGE("item 5 MSTW2008 LO, kappa, CD-Bonn  G3c: int b1 dx (x >= 0.01) = "
          << m.integral_g3c << " ; digitized CDKS "
          << close_kumano_integral(true) << " (AV18 2.24896e-4)");
  CHECK_CLOSE(m.integral_g3c, 4.485801e-4, 3e-3);

  // THE GRID IS NOT DOING IT.  The analytic form is sampled on
  // `fdeut.av18`'s own 0.1 fm^-1 spacing so that the AV18 and CD-Bonn rows
  // differ in the wave function alone; refining to 0.02 fm^-1 (5x, 1001
  // rows) moves the peak by 5e-4 relative and does not move either zero.
  // Measured here on ONE x -- the peak's own grid point -- rather than on
  // the whole scan, because each finite-|q| point costs ~0.24 s.
  std::vector<DeuteronConvolutionB1::Options> of = cdbonn_opts(1);
  of[0].unpol = std::make_shared<MstwSF>();
  of[0].cdbonn_dk_fm = 0.02;
  const DeuteronConvolutionB1 dfine(of[0]);
  const double fine = m.x_max * dfine.b1(m.x_max, kQ2, 0.0);
  MESSAGE("item 5 dk = 0.02 fm^-1 at the peak: " << fine << " against "
          << m.xb1_max << " (relative " << fine / m.xb1_max - 1.0 << ")");
  CHECK(std::fabs(fine / m.xb1_max - 1.0) < 2e-3);
}
#endif  // LIPOLGEN_HAVE_PYTHIA8

// ===================================================================== T2
TEST_CASE("b1_nuclear T2: quadrature convergence") {
  REQUIRE(have(kFdeut));
  // A BLOCKING GATE THAT SILENTLY SKIPS IS NOT A BLOCKING GATE (the header's
  // departure (i)).  The 6Li clause at the bottom is the one number the
  // default y grid can get wrong without any other test seeing it, so it must
  // not pass vacuously when the VMC tables are absent.
  REQUIRE(have(kLi6Momentum));
  REQUIRE(have(kLi6Overlap));
  const FdeutTable t = read_fdeut_k(kFdeut);
  const ClusterPartialWave p0 = ClusterPartialWave::from_uw(t.k_gev, t.u, 0);
  const ClusterPartialWave p2 = ClusterPartialWave::from_uw(t.k_gev, t.w, 2);
  const ConvolutionKinematics kin{M_NUCLEON, M_NUCLEON, t.ebind_gev, 1.0};
  // At kappa = 1 THROUGHOUT: this test is about the QUADRATURE, on which the
  // delta-function's kappa has no bearing, and the default (Eq. 21) rebuilds
  // the densities at every x, which would make a 3 x 5-point convergence scan
  // cost minutes rather than seconds.
  const DeuteronConvolutionB1& d = gate_k1();

  double prev_zero = 0.0;
  for (int f : {1, 2, 4}) {
    DeuteronConvolutionB1::Options o;
    o.finite_q_delta = false;
    o.quad.n_k = 1 + (2001 - 1) * f;
    o.quad.n_y_low *= f;
    o.quad.n_y_mid *= f;
    o.quad.n_y_high *= f;
    const DeuteronConvolutionB1 dd(o);
    for (double x : {0.05, 0.1, 0.2, 0.35, 0.5}) {
      CAPTURE(f);
      CAPTURE(x);
      CHECK_CLOSE(dd.b1(x, kQ2, 0.0), d.b1(x, kQ2, 0.0), 1e-3);
    }
    // int delta_T f dy -> 0 faster than 1/n ... except that on the default
    // y_min it does not: it is a SUPPORT truncation (see G0e), so what is
    // asserted here is that refinement does not make it worse and that the
    // raw norm converges onto the G1c moment identity.
    //
    // TWO LEVELS ONLY (f = 1, 2).  The un-renormalised densities have to be
    // built a SECOND time here -- `dd`'s own are renormalised, and `norm()`
    // reports the post-scale value -- and that second build costs 0.24 / 0.96
    // / 3.84 s at f = 1 / 2 / 4.  The "refinement does not make it worse"
    // clause is a comparison between consecutive levels and the spline pin is
    // level-independent, so both survive on {1, 2} and the f = 4 rebuild buys
    // nothing but 3.8 s.  f = 4 still runs the b1 convergence checks above.
    if (f <= 2) {
      LightConeDensities::Options q = o.quad;
      q.renormalize = false;
      const LightConeDensities raw(p0, p2, kin, q);
      const std::vector<double>& yy = raw.y_grid();
      std::vector<double> a(yy.size());
      double mx = 0.0;
      for (std::size_t i = 0; i < yy.size(); ++i) {
        a[i] = raw.delta_t_f(yy[i]);
        mx = std::max(mx, std::fabs(a[i]));
      }
      const double z = std::fabs(trapezoid(a, yy)) / mx;
      MESSAGE("T2 x" << f << ": raw norm = " << raw.norm()
                     << ", |int dT f|/max = " << z);
      CHECK(z < 2e-4);
      if (f > 1) CHECK(z < prev_zero * 1.5);
      prev_zero = z;
      // The only clause of T2 that can see an interpolation bias: refining
      // the y/k grids cannot fix the TABLE spacing, so this pins the spline.
      CHECK_CLOSE(raw.norm(), 0.98707, 1e-3);
    }
  }

  // ---- n_k PARITY.  The inner k integral is composite Simpson, which needs
  //      an ODD count; `simpson()` falls back to the trapezoid otherwise and
  //      that is 22 % low at x = 0.05 (the very bias Simpson was introduced
  //      for).  `n_k` is a plain read-write field, also exposed in Python, and
  //      a user refining the grid types 2000 or 4000 -- so the constructor
  //      bumps an even value to the next odd one and `options()` reports what
  //      actually ran.  Without that, n_k = 2000 gave -5.4649e-6 at x = 0.05
  //      against 2001's -6.9675e-6, silently.
  {
    DeuteronConvolutionB1::Options oe;
    oe.finite_q_delta = false;
    oe.quad.n_k = 2000;
    const DeuteronConvolutionB1 de(oe);
    CHECK(de.densities().options().n_k == 2001);
    MESSAGE("T2 parity: n_k = 2000 -> " << de.densities().options().n_k
            << ", b1(0.05) = " << de.b1(0.05, kQ2, 0.0) << " against "
            << d.b1(0.05, kQ2, 0.0));
    CHECK(de.b1(0.05, kQ2, 0.0) == d.b1(0.05, kQ2, 0.0));   // bit-for-bit
  }

  // ---- THE 6Li CLAUSE (added by the review of 2026-09-03).  The struck
  //      deuteron's z-width scales with M_alpha/M_d = 1.987, so more of its
  //      CANCELLING delta_T f lives outside the dense (0.85, 1.15) segment
  //      than the A = 2 deuteron's does, and term (2d) is the one number in
  //      this design that the default y grid can get wrong without any test
  //      seeing it (T2 above exercises only the deuteron; T14/T15 pin at
  //      rtol 3e-2).  On the design's 600/2400/800 it was 2.3 % high at
  //      x = 0.10; on the 2400/2400/3200 default it is within 1.2e-3 of the
  //      x4-refined value, which is what this pins.
  {
    const Li6ConvolutionB1& m = li6();
    Li6ConvolutionOptions o4;
    o4.quad.n_y_low *= 4;
    o4.quad.n_y_mid *= 4;
    o4.quad.n_y_high *= 4;
    const Li6ConvolutionB1 m4(o4);
    for (double x : {0.05, 0.10, 0.30}) {
      CAPTURE(x);
      const double a = m.b1_alpha_d_dwave_d(x, kQ2, 0.0);
      const double b = m4.b1_alpha_d_dwave_d(x, kQ2, 0.0);
      MESSAGE("T2 6Li x=" << x << ": (2d) = " << a << " against " << b
                          << " on a x4 y grid, " << (a / b - 1.0) << " apart");
      CHECK_CLOSE(a, b, 3e-3);
      // ... and (2a)/(2d) then sits on the analytic 2 (M_d/M_alpha)^2 =
      // 0.5064 rather than the 0.489-0.497 the coarse grid gave (T15).  The
      // ratio is a ratio of two cancellations and converges more slowly than
      // (2d) itself: 0.5024 here against 0.5053 on the x4 grid at x = 0.05.
      CHECK_CLOSE(m.b1_alpha_d_dwave_alpha(x, kQ2, 0.0) / a,
                  m4.b1_alpha_d_dwave_alpha(x, kQ2, 0.0) / b, 1e-2);
    }
  }
}

// ===================================================================== T3, T4
TEST_CASE("b1_nuclear T3/T4: the alpha-d z-density normalisation and <z>") {
  REQUIRE(have(kLi6Momentum));
  REQUIRE(have(kLi6Overlap));
  const auto w = li6_alpha_d_partial_waves();
  const double md = nuclear_mass(1, 2), ma = nuclear_mass(2, 4);
  const double eps = LI6_ALPHA_TAG().separation_energy;
  const ConvolutionKinematics kin{md, ma, eps, 1.0};
  LightConeDensities::Options o;
  o.renormalize = false;
  const LightConeDensities raw(w.first, w.second, kin, o);

  // The trapezoid S + D of the table against the file's own printed norms.
  // NOT 1e-9 (that would be tautological once `renormalize` sets it): 3e-4 is
  // the file's own quadrature spread, design 2.1.
  const double tab = w.first.norm2() + w.second.norm2();
  MESSAGE("T3: table trapezoid S+D = " << tab << " against VMC_N_ALPHA_D_LI6 = "
                                       << VMC_N_ALPHA_D_LI6);
  CHECK_CLOSE(tab, VMC_N_ALPHA_D_LI6, 3e-4);
  const std::vector<double> printed = read_anl_momentum_norms(kLi6Momentum);
  REQUIRE(printed.size() == 3);
  CHECK_CLOSE(printed[1] + printed[2], VMC_N_ALPHA_D_LI6, 1e-12);

  // raw norm() = N_ad x <E>/M_d, both sides from the same table.
  const std::vector<double> kk = linspace(0.0, w.first.k.back(), 100001);
  std::vector<double> a, b, m0(kk.size()), m2(kk.size());
  w.first.eval_sorted(kk, &a);
  w.second.eval_sorted(kk, &b);
  for (std::size_t i = 0; i < kk.size(); ++i) {
    const double r = kk[i] * kk[i] * (a[i] * a[i] + b[i] * b[i]);
    m0[i] = r;
    m2[i] = kk[i] * kk[i] * r;
  }
  const double k2 = trapezoid(m2, kk) / trapezoid(m0, kk);
  const double mean_e = 1.0 - eps / md - k2 / (2.0 * ma * md);
  MESSAGE("T3: raw norm = " << raw.norm() << " ; N_ad x <E>/M_d = "
                            << VMC_N_ALPHA_D_LI6 * mean_e << " (design 0.8174)");
  CHECK_CLOSE(raw.norm(), 0.8174, 1e-3);
  CHECK_CLOSE(raw.norm(), 0.8176782, 1e-5);
  CHECK_CLOSE(raw.norm(), VMC_N_ALPHA_D_LI6 * mean_e, 1e-3);

  // p_d_momentum() is the file's own D/(S+D).  The design asks rtol 1e-4; the
  // table's trapezoid is 0.0193516 against VMC_P_D_LI6 = 0.0193549, i.e.
  // 1.7e-4 -- which design 2.1 itself states.  Pinned at 3e-4.
  MESSAGE("T3: p_d_momentum = " << raw.p_d_momentum()
                                << " against VMC_P_D_LI6 = " << VMC_P_D_LI6);
  CHECK_CLOSE(raw.p_d_momentum(), VMC_P_D_LI6, 3e-4);
  // p_d() is the y-WEIGHTED one, P_D <y>_D/<y>; ABSOLUTE tolerance, because
  // the two differ by ~1e-3 relative BY CONSTRUCTION.
  CHECK_CLOSE_AT(raw.p_d(), VMC_P_D_LI6, 0.0, 1e-4);
  CHECK_CLOSE(raw.p_d(), 0.0193299, 1e-4);

  // T4 -- <z> is <z^2>/<z>, NOT <z>.
  std::vector<double> e1(kk.size()), e2(kk.size());
  for (std::size_t i = 0; i < kk.size(); ++i) {
    const double e = md - eps - kk[i] * kk[i] / (2.0 * ma);
    e1[i] = m0[i] * e;
    e2[i] = m0[i] * e * e;
  }
  const double nn = trapezoid(m0, kk);
  const double me = trapezoid(e1, kk) / nn, me2 = trapezoid(e2, kk) / nn;
  const double ident = (me2 + k2 / 3.0) / (md * me);
  MESSAGE("T4: mean_y = " << raw.mean_y() << " ; identity " << ident
                          << " ; closed form 1 - eps/M_d - <k2>/(2 M_a M_d) = "
                          << mean_e);
  CHECK_CLOSE(raw.mean_y(), 0.99979, 2e-3);
  CHECK_CLOSE(raw.mean_y(), 0.9997871, 1e-5);
  CHECK_CLOSE(raw.mean_y(), ident, 1e-5);
  // The closed form is a DIFFERENT number: it differs by Var/<z> ~ <k2>/(3M_d^2).
  CHECK(std::fabs(raw.mean_y() - mean_e) > 1e-3);
  CHECK_CLOSE(raw.mean_y() - mean_e, k2 / (3.0 * md * md), 5e-2);
}

// ===================================================================== T5
TEST_CASE("b1_nuclear T5: the D-wave terms vanish as P_D(alpha-d) -> 0") {
  REQUIRE(have(kLi6Momentum));
  const auto w = li6_alpha_d_partial_waves();
  ClusterPartialWave zero2 = w.second;
  for (double& v : zero2.phi) v = 0.0;
  zero2.rebuild();
  const double md = nuclear_mass(1, 2), ma = nuclear_mass(2, 4);
  const double eps = LI6_ALPHA_TAG().separation_energy;
  LightConeDensities::Options q;
  q.norm_target = VMC_N_ALPHA_D_LI6;
  const LightConeDensities dd(w.first, zero2, {md, ma, eps, 1.0}, q);
  const LightConeDensities da(w.first, zero2, {ma, md, eps, 1.0}, q);
  const Li6ConvolutionB1 m(Li6ConvolutionOptions(), dd, da);
  for (double x : {0.05, 0.1, 0.3, 0.5, 0.7}) {
    CAPTURE(x);
    CHECK(m.b1_alpha_d_dwave_d(x, kQ2, 0.0) == 0.0);
    CHECK(m.b1_alpha_d_dwave_alpha(x, kQ2, 0.0) == 0.0);
    CHECK(m.b1_cg_dwave(x, kQ2, 0.0) == 0.0);
    CHECK(m.b1(x, kQ2, 0.0) == m.b1_embedded_s(x, kQ2, 0.0));   // bit-for-bit
  }
  // The test-only constructor REFUSES `finite_q_delta`: that path rebuilds
  // both density sets per x from the FILE waves, so the zero-D-wave pair just
  // injected would be silently discarded and the D wave would grow back.
  Li6ConvolutionOptions oq;
  oq.finite_q_delta = true;
  CHECK_THROWS_AS(Li6ConvolutionB1(oq, dd, da), std::runtime_error);
}

// ============================================================== T5b (review)
// The support of b1_d is b1_d's business, and `convolve` now asks for it.
//
// The header promised "this routine does not cut at x/y = 1" while `lo =
// max(x, ...)` did exactly that; the digitized CDKS column runs to x = 1.59
// and the deuteron's per-nucleon support genuinely exceeds 1, so the cut
// dropped real strength from terms (1) and (3) near the top of the window --
// 5 % at x = 0.95, which is the cell `--x-max 0.955` sits in.
TEST_CASE("b1_nuclear T5b: convolve keeps b1_d's support above x/y = 1") {
  REQUIRE(have(kLi6Momentum));
  const Li6ConvolutionB1& m = li6();
  CHECK(m.deuteron_b1()->b1(1.60, kQ2, 0.0) == 0.0);      // zero past the table
  CHECK(m.deuteron_b1()->b1(1.50, kQ2, 0.0) != 0.0);
  Li6ConvolutionOptions oc;
  oc.deuteron_b1_x_max = 1.0;                             // the old behaviour
  const Li6ConvolutionB1 cut(oc);
  const double xs[5] = {0.30, 0.50, 0.80, 0.90, 0.95};
  const double want[5] = {0.0, 0.0, 3.00e-4, 1.5504e-2, 5.0788e-2};
  for (int i = 0; i < 5; ++i) {
    CAPTURE(xs[i]);
    const double full = m.b1_embedded_s(xs[i], kQ2, 0.0);
    const double none = cut.b1_embedded_s(xs[i], kQ2, 0.0);
    const double dropped = 1.0 - none / full;
    MESSAGE("T5b x=" << xs[i] << ": cutting at x/y = 1 drops "
                     << 100.0 * dropped << " % of term (1)");
    CHECK_CLOSE_AT(dropped, want[i], 2e-2, 1e-5);
  }
  // Nothing below x = 0.8 (the density's own z support ends there), and the
  // whole effect is inside terms (1) and (3): the F1 slots have x_max_g = 1
  // and are untouched.
  CHECK(m.b1_alpha_d_dwave_d(0.95, kQ2, 0.0)
        == cut.b1_alpha_d_dwave_d(0.95, kQ2, 0.0));
  CHECK(m.b1_alpha_d_dwave_alpha(0.95, kQ2, 0.0)
        == cut.b1_alpha_d_dwave_alpha(0.95, kQ2, 0.0));
}

// ===================================================================== T6
TEST_CASE("b1_nuclear T6: the no-smearing limit IS b1_li6_from_deuteron") {
  REQUIRE(have(kLi6Momentum));
  Li6ConvolutionOptions o;
  const auto b1d = cdks_b1_raw_per_nucleon();
  o.deuteron_b1 = b1d;
  for (double pd : {0.0, VMC_P_D_LI6, P_D_LI6, 0.2}) {
    const Li6ConvolutionB1 m(o, LightConeDensities::delta_limit(pd, 1.0),
                             LightConeDensities::delta_limit(pd, 1.0));
    for (double x : {0.05, 0.2, 0.5, 0.9}) {
      CAPTURE(pd);
      CAPTURE(x);
      CHECK_CLOSE(m.b1(x, kQ2, 0.0),
                  b1_li6_from_deuteron(b1d->b1(x, kQ2, 0.0), 1.0 - 0.9 * pd,
                                       LI6_B1_PER_NUCLEON),
                  1e-12);
    }
  }
  // ... and at P_D = P_D_LI6 the weight equals LI6_B1_RANK2_TRANSFER, which is
  // the analytic bridge from this backend to the constant the default uses.
  CHECK_CLOSE(1.0 - 0.9 * P_D_LI6, LI6_B1_RANK2_TRANSFER, 1e-4);
  // The VMC P_D gives 0.98258 instead -- design 1.7's table.
  CHECK_CLOSE(1.0 - 0.9 * VMC_P_D_LI6, 0.98258, 1e-4);
}

// ===================================================================== T7
TEST_CASE("b1_nuclear T7: the band and the term knobs") {
  REQUIRE(have(kLi6Momentum));
  const Li6ConvolutionB1& m = li6();
  const auto b0 = m.banded(0.0), b1 = m.banded(1.0), b2 = m.banded(2.0);
  for (double x : {0.05, 0.1, 0.3, 0.5, 0.7}) {
    CAPTURE(x);
    CHECK(b0->b1(x, kQ2, 0.0) == 0.0);
    CHECK(b2->b1(x, kQ2, 0.0) == 2.0 * b1->b1(x, kQ2, 0.0));   // bit-for-bit
    CHECK(b1->b1(x, kQ2, 0.0) == m.b1(x, kQ2, 0.0));
    // the four terms sum to b1, and (2d)+(2a) is the orbital term
    const double s = m.b1_embedded_s(x, kQ2, 0.0);
    const double dd = m.b1_alpha_d_dwave_d(x, kQ2, 0.0);
    const double da = m.b1_alpha_d_dwave_alpha(x, kQ2, 0.0);
    const double cgv = m.b1_cg_dwave(x, kQ2, 0.0);
    CHECK_CLOSE(s + dd + da + cgv, m.b1(x, kQ2, 0.0), 1e-12);
    CHECK(m.b1_alpha_d_dwave(x, kQ2, 0.0) == dd + da);
    // SD + DD == (2d) + (2a)
    CHECK_CLOSE(m.b1_alpha_d_sd(x, kQ2, 0.0) + m.b1_alpha_d_dd(x, kQ2, 0.0),
                dd + da, 1e-12);
  }
  struct Knob { const char* tag; Li6ConvolutionOptions o; };
  Li6ConvolutionOptions ks, ka, kc;
  ks.w_embedded_s = 0.0;
  ka.w_alpha_d_dwave = 0.0;
  kc.w_cg_dwave = 0.0;
  const Li6ConvolutionB1 ms(ks), ma(ka), mc(kc);
  for (double x : {0.1, 0.3, 0.5}) {
    CAPTURE(x);
    CHECK(ms.b1_embedded_s(x, kQ2, 0.0) == 0.0);
    CHECK_CLOSE(ms.b1_cg_dwave(x, kQ2, 0.0), m.b1_cg_dwave(x, kQ2, 0.0), 1e-12);
    CHECK(mc.b1_cg_dwave(x, kQ2, 0.0) == 0.0);
    CHECK_CLOSE(mc.b1_embedded_s(x, kQ2, 0.0), m.b1_embedded_s(x, kQ2, 0.0), 1e-12);
    // w_alpha_d_dwave kills BOTH (2d) and (2a) and nothing else
    CHECK(ma.b1_alpha_d_dwave_d(x, kQ2, 0.0) == 0.0);
    CHECK(ma.b1_alpha_d_dwave_alpha(x, kQ2, 0.0) == 0.0);
    CHECK_CLOSE(ma.b1_embedded_s(x, kQ2, 0.0), m.b1_embedded_s(x, kQ2, 0.0), 1e-12);
    CHECK_CLOSE(ma.b1_cg_dwave(x, kQ2, 0.0), m.b1_cg_dwave(x, kQ2, 0.0), 1e-12);
  }
  // The spectroscopic-factor knob: 1/N_ad = 1.22 on the whole answer.
  Li6ConvolutionOptions kn;
  kn.use_spectroscopic_factor = false;
  const Li6ConvolutionB1 mn(kn);
  for (double x : {0.1, 0.5}) {
    CAPTURE(x);
    CHECK_CLOSE(mn.b1(x, kQ2, 0.0) / m.b1(x, kQ2, 0.0), 1.0 / VMC_N_ALPHA_D_LI6,
                1e-9);
  }
  // ... and the +-5 % N_ad systematic, quoted through norm_target.
  Li6ConvolutionOptions kp;
  kp.norm_target = 1.05 * VMC_N_ALPHA_D_LI6;
  const Li6ConvolutionB1 mp(kp);
  CHECK_CLOSE(mp.b1(0.3, kQ2, 0.0) / m.b1(0.3, kQ2, 0.0), 1.05, 1e-9);
}

// ===================================================================== T8
TEST_CASE("b1_nuclear T8: determinism") {
  REQUIRE(have(kLi6Momentum));
  const Li6ConvolutionB1 a{}, b{};
  const std::vector<double> xs = linspace(0.02, 0.95, 200);
  for (double x : xs) {
    CAPTURE(x);
    CHECK(a.b1(x, kQ2, 0.0) == b.b1(x, kQ2, 0.0));   // bit-identical
  }
}

// ===================================================================== T10
TEST_CASE("b1_nuclear T10: the alpha-d sign gate") {
  REQUIRE(have(kLi6Momentum));
  REQUIRE(have(kLi6Overlap));
  REQUIRE(have(kFdeut));
  const auto w = li6_alpha_d_partial_waves();
  const double q = alpha_d_quadrupole_fm2(w.first, w.second);
  MESSAGE("T10: Q(alpha-d relative motion) = " << q << " fm^2 -- MUST be < 0");
  CHECK(q < 0.0);
  // The node crossing, pinned at both signs: phi_0 phi_2 (CDKS convention) is
  // POSITIVE below the S node (0.678 fm^-1) and NEGATIVE between the nodes.
  const double k02 = 0.2 * HBARC_GEV_FM, k10 = 1.0 * HBARC_GEV_FM;
  CHECK(w.first(k02) * w.second(k02) > 0.0);
  CHECK(w.first(k10) * w.second(k10) < 0.0);
  // ... and above the D node (2.25 fm^-1) it is positive again.
  const double k30 = 3.0 * HBARC_GEV_FM;
  CHECK(w.first(k30) * w.second(k30) > 0.0);
  // The DEUTERON's phi_0 phi_2 is negative at both, from the SAME code path.
  const FdeutTable t = read_fdeut_k(kFdeut);
  const ClusterPartialWave p0 = ClusterPartialWave::from_uw(t.k_gev, t.u, 0);
  const ClusterPartialWave p2 = ClusterPartialWave::from_uw(t.k_gev, t.w, 2);
  CHECK(p0(k02) * p2(k02) < 0.0);
  CHECK(p0(k10) * p2(k10) < 0.0);
  // and the same quadrupole code on the deuteron reproduces the file's own
  // header value qm = 0.269673 fm^2 -- the machinery's validation.
  const double qd = alpha_d_quadrupole_fm2(p0, p2);
  MESSAGE("T10: Q(deuteron) from the same code = " << qd
          << " fm^2 against fdeut.av18's own qm = 0.269673");
  CHECK_CLOSE(qd, 0.269673, 3e-3);
  CHECK(qd > 0.0);
}

// ===================================================================== T13
// THE A = 7 GATE.  Q(7Li) from the alpha-t overlap against the MEASURED
// moment, and it validates the alpha-t WAVE FUNCTION's quadrupole ONLY --
// its <r^2> and its P-wave character.  It is NOT a pass/fail on b1(7Li):
// b1(7Li) IS NOT IMPLEMENTED (open item 15; 7Li's whole rank-2 sector is
// exactly zero by construction because `default_inclusive_kernel` fills no
// spin-3/2 slot), it waits on the unpolarised-backend decision, and nothing
// here licenses one.  Nor does the gate say anything about the light-cone
// convolution or the DIS input, which is where the real 7Li uncertainty
// lives; and there is NO A = 3 analogue of the A = 2 b1 gate, because 3H and
// 3He are J = 1/2 and carry no rank-2 structure function at all.  So this is
// the only offline validation the 7Li wave function gets.
//
// Reference: docs/open_items/run_2026-09-03/phase_D_li7_rank2.md sec. 4.2
// (the construction) and docs/open_items/run_2026-09-06/phase_B_numbers.md
// sec. B3 (this measurement, and the sourcing of `LI7_QUADRUPOLE_FM2`).
TEST_CASE("b1_nuclear T13: the A = 7 alpha-t quadrupole gate") {
  REQUIRE(have(kLi7Momentum));

  // (0) the measured moment has ONE home and this is its value: TUNL's
  //     A = 5, 6, 7 evaluation (Tilley et al., NPA 708 (2002) 3) prints
  //     Q(7Li) = -40.6 +- 0.8 mb, and 1 mb = 0.1 fm^2.  The SAME paper's
  //     A = 6 half is where LI6_QUADRUPOLE_FM2 = -0.0818 comes from, so the
  //     two lithium quadrupoles share a document, a sign convention and a
  //     unit rule.  See rc.hpp for the eight-entry Stone compilation spread.
  CHECK_CLOSE(LI7_QUADRUPOLE_FM2, -40.6 * 0.1, 1e-12);
  CHECK(LI7_QUADRUPOLE_FM2 < 0.0);

  const AlphaTQuadrupole g = li7_alpha_t_quadrupole();

  // (1) Z_eff is the MASS-weighted charge, not the A-number one.
  const double m_t = nuclear_mass(1, 3), m_a = nuclear_mass(2, 4),
               m_7 = nuclear_mass(3, 7);
  CHECK_CLOSE(g.z_eff,
              2.0 * (m_t / m_7) * (m_t / m_7) + (m_a / m_7) * (m_a / m_7),
              1e-14);
  CHECK_CLOSE(g.z_eff, 0.695075101362, 1e-9);
  // the A-number form 34/49 is 0.172 % away, MEASURED -- and is not used.
  CHECK_CLOSE_AT(34.0 / 49.0 / g.z_eff, 1.0, 0.0, 2.0e-3);
  CHECK(std::fabs(34.0 / 49.0 / g.z_eff - 1.0) > 1.0e-3);

  // (2) the pin, at the shipped defaults (r_max = 30 fm, n_r = 4000,
  //     n_k = 8001).  These reproduce the offline note's own numpy recipe
  //     (phase_D_li7_rank2.md sec. 9 `q_at`) to the LAST BIT -- measured
  //     2026-09-06, relative difference 0.000e+00 on all three.
  CHECK_CLOSE(g.r2_fm2, 12.534828961030, 1e-9);
  CHECK_CLOSE(g.r_rms_fm, 3.540456038568, 1e-9);
  CHECK_CLOSE(g.q_fm2, -3.485059004257, 1e-9);
  CHECK_CLOSE(g.q_fm2, -0.4 * g.z_eff * g.r2_fm2, 1e-14);
  CHECK(g.q_fm2 < 0.0);                        // same sign as the measurement

  // (3) grid/truncation, PINNED and converged: 30 fm is within 0.032 % of
  //     60 fm, so the 14 % discrepancy below is physics, not the r cut.
  const double q20 = li7_alpha_t_quadrupole(0.0, 20.0).q_fm2;
  const double q40 = li7_alpha_t_quadrupole(0.0, 40.0).q_fm2;
  const double q60 = li7_alpha_t_quadrupole(0.0, 60.0).q_fm2;
  CHECK_CLOSE(q20, -3.472004271671, 1e-9);
  CHECK_CLOSE(q40, -3.485797071382, 1e-9);
  CHECK_CLOSE(q60, -3.486171274781, 1e-9);
  CHECK(std::fabs(q60 - g.q_fm2) / std::fabs(g.q_fm2) < 1e-3);

  // (4) the ANL Monte Carlo band, fully correlated, +-0.33 % -- an order of
  //     magnitude smaller than the discrepancy, so it does not explain it.
  const double qm = li7_alpha_t_quadrupole(-1.0).q_fm2;
  const double qp = li7_alpha_t_quadrupole(+1.0).q_fm2;
  CHECK_CLOSE(qm, -3.496426569222, 1e-9);
  CHECK_CLOSE(qp, -3.473744615578, 1e-9);
  CHECK(std::fabs(qm - g.q_fm2) / std::fabs(g.q_fm2) < 4e-3);
  CHECK(std::fabs(qp - g.q_fm2) / std::fabs(g.q_fm2) < 4e-3);

  // (5) THE GATE ITSELF, REPORTED AS A RATIO AND NEVER AS A b1 VERDICT.
  //     0.858389 against TUNL's -4.06 fm^2 (14.2 % low); 0.871265 against
  //     the Voelk CER -4.00(3) fm^2 that Stone also lists (12.9 % low).  The
  //     choice of compilation moves it by 1.5 % and does not change the
  //     verdict, which is "13-14 % low, converged to 0.03 %, MC band 0.33 %".
  const double ratio = g.q_fm2 / LI7_QUADRUPOLE_FM2;
  MESSAGE("T13: Q(7Li)_alpha-t = " << g.q_fm2 << " fm^2 against the measured "
          << LI7_QUADRUPOLE_FM2 << " fm^2 -- ratio " << ratio
          << ".  This validates the alpha-t WAVE FUNCTION's quadrupole only; "
             "b1(7Li) remains UNIMPLEMENTED (open item 15).");
  CHECK_CLOSE(ratio, 0.858388917305, 1e-9);
  CHECK_CLOSE(g.q_fm2 / -4.00, 0.871264751064, 1e-9);
  CHECK(ratio > 0.80);
  CHECK(ratio < 1.00);            // the cluster overlap UNDERSHOOTS, as noted

  // (6) the contrast that makes the gate worth committing: the SAME
  //     construction one cluster level down (alpha-d, 6Li) misses its own
  //     measured moment by a factor 4.07, so the 7Li wave-function input is
  //     validated by a measured moment an ORDER OF MAGNITUDE better than
  //     6Li's is.  That asymmetry, not the 14 %, is the result.
  if (have(kLi6Momentum) && have(kLi6Overlap)) {
    const auto w6 = li6_alpha_d_partial_waves();
    const double f6 =
        alpha_d_quadrupole_fm2(w6.first, w6.second) / LI6_QUADRUPOLE_FM2;
    MESSAGE("T13: the 6Li contrast from the same code family -- 6Li misses "
            "by a factor " << f6 << " (" << 100.0 * (f6 - 1.0)
            << " %) against 7Li's " << 100.0 * (1.0 - ratio) << " %");
    CHECK(f6 > 4.0);
    // MEASURED 2026-09-06: 307.5 % against 14.2 %, a factor 21.7 -- so
    // "an order of magnitude better" is a measurement, not a manner of
    // speaking.  The threshold is 10, well inside it.
    CHECK((f6 - 1.0) > 10.0 * std::fabs(1.0 - ratio));
  }
}

// ===================================================================== T12
TEST_CASE("b1_nuclear T12: the VMC spread, documented not explained") {
  REQUIRE(have(kLi6Momentum));
  CHECK_CLOSE(VMC_N_ALPHA_D_LI6, 0.819481, 1e-12);
  CHECK_CLOSE(VMC_P_D_LI6, 0.015861 / VMC_N_ALPHA_D_LI6, 0.0);
  // 1. Wiringa, Schiavilla, Pieper, Carlson, PRC 89 (2014) 024305 sec. III:
  //    "The integrated N_ad = 0.86 is a sum of S- and D-wave parts of 0.846
  //    and 0.017" -> 0.863, P_D = 0.0197.
  const double n_pub = 0.846 + 0.017, pd_pub = 0.017 / n_pub;
  // 2. the 2004 overlap file li6.ad prints ndx s-wave d-wave = 0.856 0.838
  //    0.017.
  const double n_2004 = 0.856;
  MESSAGE("T12: three tabulations of N_ad -- 0.819481 (2014 momenta, used), "
          << n_2004 << " (2004 overlap), " << n_pub << " (published); 5 % "
          "spread, UNEXPLAINED, carried as a systematic (not a version "
          "difference: 1 and 3 are the same year and Hamiltonian family)");
  CHECK(std::fabs(VMC_N_ALPHA_D_LI6 / n_pub - 1.0) < 0.06);
  CHECK(std::fabs(VMC_N_ALPHA_D_LI6 / n_2004 - 1.0) < 0.05);
  CHECK(std::fabs(VMC_P_D_LI6 / pd_pub - 1.0) < 0.02);
  // VMC_S_ALPHA_D_LI6 is the file's TOTAL block: a THIRD value, deliberately
  // not unified with N_ad (P_D must be S_2/N with the same N).
  CHECK(VMC_S_ALPHA_D_LI6 != VMC_N_ALPHA_D_LI6);
  CHECK_CLOSE(VMC_S_ALPHA_D_LI6, VMC_N_ALPHA_D_LI6, 4e-4);
}

// ===================================================================== T13
TEST_CASE("b1_nuclear T13: CDKS Eq. (22) is NOT f1_from_f2") {
  // REQUIRE, not skip (the header's departure (i)): the second half of this
  // case is the check that `f1a` did not sneak back into the kernel, and it
  // must not pass vacuously when the VMC tables are absent.
  REQUIRE(have(kLi6Momentum));
  const ToyF2 toy;
  for (double x : {0.1, 0.3, 0.5, 0.8}) {
    for (double q2 : {1.5, 2.5, 10.0}) {
      CAPTURE(x);
      CAPTURE(q2);
      const double f2 = toy.f2p(x, q2);
      CHECK_CLOSE(f1_cdks(f2, x, q2) / toy.f1_from_f2(f2, x, q2),
                  1.0 + gamma_squared(x, q2), 1e-14);
    }
  }
  CHECK_CLOSE(1.0 + gamma_squared(0.8, kQ2), 1.90, 1e-2);
  // ... and the kernel's own F1_d is that same factor above NuclearF2::f1a/2,
  // which is the check that f1a did not sneak back in.
  {
    const Li6ConvolutionB1& m = li6();
    const auto unp = std::make_shared<const ToyF2>();
    const NuclearF2 nd(DEUTERON(), unp, nullptr, nullptr);
    for (double x : {0.3, 0.5, 0.8}) {
      CAPTURE(x);
      CHECK_CLOSE(m.f1_deuteron(x, kQ2) / (nd.f1a(x, kQ2) / 2.0),
                  1.0 + gamma_squared(x, kQ2), 1e-12);
    }
    CHECK_CLOSE(1.0 + gamma_squared(0.5, kQ2), 1.35, 1e-2);
    // F1_alpha is the SAME isoscalar (no alpha EMC effect, A9).
    CHECK(m.f1_alpha(0.5, kQ2) == m.f1_deuteron(0.5, kQ2));
  }
}

// ===================================================================== T14
TEST_CASE("b1_nuclear T14: the angular-average truncation, MEASURED") {
  REQUIRE(have(kLi6Momentum));
  const Li6ConvolutionB1& m = li6();
  // Pinned to what they come out at, so that the number in the header comment
  // and in OPEN_ITEMS_SOLUTIONS.md is a MEASUREMENT and not the first
  // revision's asserted "0.2 %".
  const double xs[5] = {0.05, 0.10, 0.20, 0.30, 0.50};
  const double p2w[5] = {-2.3131e-5, -8.7351e-6, -1.7072e-6, 1.4106e-5, -5.1059e-5};
  const double p4w[5] = {5.5027e-7, -9.6501e-6, -1.0349e-6, -4.6105e-7, -1.3614e-5};
  for (int i = 0; i < 5; ++i) {
    const double s = m.b1_embedded_s(xs[i], kQ2, 0.0);
    const double r2 = m.b1_cg_dwave_p2_remainder(xs[i], kQ2, 0.0) / s;
    const double r4 = m.b1_cg_dwave_p4_remainder(xs[i], kQ2, 0.0) / s;
    MESSAGE("T14 x=" << xs[i] << ": P2 remainder / term(1) = " << r2
                     << " ; P4 remainder / term(1) = " << r4);
    CAPTURE(xs[i]);
    // rtol 3e-2, NOT the design's 1e-6: both remainders are themselves
    // cancellations (the P2- and P4-weighted densities integrate to zero), so
    // the P4 ratio at x = 0.05 moves by 2 % between the default grid and a
    // 4x-refined one.  What is physics here is the ORDER OF MAGNITUDE, and
    // that is what the |.| < 1e-3 clause pins.
    CHECK_CLOSE(r2, p2w[i], 3e-2);
    CHECK_CLOSE(r4, p4w[i], 3e-2);
    // both are far inside the mandatory 100 % band
    CHECK(std::fabs(r2) < 1e-3);
    CHECK(std::fabs(r4) < 1e-3);
  }
}

// ===================================================================== T15
TEST_CASE("b1_nuclear T15: the struck-alpha term is real") {
  REQUIRE(have(kLi6Momentum));
  const Li6ConvolutionB1& m = li6();
  const double xs[6] = {0.05, 0.10, 0.20, 0.30, 0.50, 0.70};
  // Re-pinned 2026-09-03 with the 2400/2400/3200 y grid: the design's
  // 600/2400/800 gave 0.489-0.497 here, i.e. 2 % below the analytic 0.5064,
  // because the struck deuteron's wider z-distribution put a bigger share of
  // its cancelling delta_T f in the coarse outer segments (T2's 6Li clause).
  const double want[6] = {0.50245, 0.50349, 0.50326, 0.50113, 0.49361,
                          0.50294};
  for (int i = 0; i < 6; ++i) {
    const double dd = m.b1_alpha_d_dwave_d(xs[i], kQ2, 0.0);
    const double da = m.b1_alpha_d_dwave_alpha(xs[i], kQ2, 0.0);
    const double r = da / dd;
    MESSAGE("T15 x=" << xs[i] << ": (2a)/(2d) = " << r);
    CAPTURE(xs[i]);
    CHECK(dd * da > 0.0);                       // same sign
    CHECK_CLOSE(r, want[i], 3e-2);
    // The design's window, unchanged: the analytic scaling is
    // 2 (M_d/M_alpha)^2 = 0.5064 and the quadrature now gives 0.494-0.503 on
    // it (0.489-0.503 on the design's coarser y grid).  A regression that
    // drops term (2a) then fails loudly rather than quietly shifting b1 by
    // 30 %.
    CHECK(r > 0.45);
    CHECK(r < 0.60);
  }
  const double md = nuclear_mass(1, 2), ma = nuclear_mass(2, 4);
  CHECK_CLOSE(2.0 * (md / ma) * (md / ma), 0.506, 3e-3);
}

// ===================================================================== T9
// The library DEFAULT is unchanged, bit for bit.
//
// This is the test the whole `B1Model` wiring is built behind: it pins what
// `default_inclusive_kernel(LI6())` -- the function `Pipeline` calls when
// `PipelineConfig::kernel` is null -- puts in the tensor slots, so that
// "`--b1-model miller` is today's default bit for bit" is a MEASURED
// statement and not a claim.
//
// Why a new reference file was needed (design section 8, Agent B step 0).
// `validation/reference/*.json` are dumped from the PYTHON generator by
// `validation/dump_polligen_reference.py`, whose kernels use
// `toy_b1(mode = kToy)`; `tests/test_reference.cpp::build_kernel` builds its
// own `InclusiveKernel::Options` and never calls `default_inclusive_kernel`.
// So the existing rtol-1e-12 gates do not exercise this path at all, and
// without `b1_default_li6.json` this test would be tautological (it would
// only compare the function with itself).
//
// Regenerate with `python3 validation/dump_b1_default_li6.py` -- and only
// when a default is deliberately changed.

TEST_CASE("b1_nuclear T9: default_inclusive_kernel(LI6()) is unchanged") {
  jsonmin::Value doc;
  const std::string path =
      std::string(LIPOLGEN_REFERENCE_DIR) + "/b1_default_li6.json";
  // NOT a skip: this file is in the tree and it is the only guard on the
  // default b1 path (see the block comment above).
  REQUIRE_MESSAGE(jsonmin::load_file(path, doc),
                  "validation/reference/b1_default_li6.json is missing -- "
                  "regenerate with validation/dump_b1_default_li6.py");

  CHECK(doc["kernel"]["ion"].str() == "6Li");
  CHECK(doc["kernel"]["b1_model"].str() == "miller");

  const auto kernel = default_inclusive_kernel(LI6());
  REQUIRE(kernel);
  CHECK(kernel->ion().A == 6);
  CHECK(kernel->target_mass() == doc["kernel"]["target_mass"].boolean());
  CHECK(kernel->tensor_gamma() == doc["kernel"]["tensor_gamma"].boolean());

  const jsonmin::Value& xs = doc["x"];
  REQUIRE(xs.size() == 200);
  REQUIRE(doc["tables"].size() == 3);
  std::size_t compared = 0;
  for (std::size_t b = 0; b < doc["tables"].size(); ++b) {
    const jsonmin::Value& blk = doc["tables"][b];
    const double q2 = blk["q2"].num();
    for (std::size_t i = 0; i < xs.size(); ++i) {
      const double x = xs[i].num();
      const SFTables t = kernel->tables(x, q2);
      CAPTURE(x);
      CAPTURE(q2);
      CHECK_CLOSE_AT(t.f1, blk["f1"][i].num(), 1e-12, 0.0);
      CHECK_CLOSE_AT(t.b1, blk["b1"][i].num(), 1e-12, 0.0);
      CHECK_CLOSE_AT(t.b2, blk["b2"][i].num(), 1e-12, 0.0);
      CHECK_CLOSE_AT(t.delta, blk["delta"][i].num(), 1e-12, 0.0);
      ++compared;
    }
  }
  CHECK(compared == 600);
  MESSAGE("T9: " << compared << " (x, Q2) rows of the DEFAULT kernel's "
                    "f1/b1/b2/delta at rtol 1e-12");
}

// ==================================================================== T11
// Channel legality: `PipelineConfig::validate()` is the ONLY guard.
//
// design_D_b1_li6.md sec. 4.1 decided this deliberately, and it is worth
// restating because the first revision of the design specified something that
// cannot be written: `struck_cluster_kernel(const TaggedChannel&, const
// StruckClusterOptions&)` takes NO `B1Model`, so there is nothing there to
// reject, and a "struck_cluster_kernel must refuse B1Model::Li6Convolution"
// clause would not compile.  The guard is therefore
//  * `validate()`, which is the only route from the CLI, and
//  * `StruckClusterOptions::inclusive_b1` installing `toy_b1` and nothing
//    else (the comment in `struck_cluster_kernel` says why).
// A caller who hand-builds an `InclusiveKernel` with a `Li6ConvolutionB1`
// b1_func and feeds it to a tagged sampler is out of the library's reach,
// exactly as they are today for any other hand-built kernel.
//
// THE RULE KEYS ON THE ISOTOPE, NOT ON THE SPIN.  The inclusive channel
// accepts isotope "d", which is spin 1, so a spin test would let a deuteron
// beam run the 6Li alpha-d convolution -- N_ad, the alpha-d densities, the
// 2/6 and 4/6 counting factors -- on a deuteron.

TEST_CASE("b1_nuclear T11: b1_model channel and isotope legality") {
  auto base = [] {
    PipelineConfig c;
    c.channel = PipelineChannel::Inclusive;
    c.isotope = "6Li";
    c.n_events = 100;
    return c;
  };

  // ---- the legal combination
  {
    PipelineConfig c = base();
    c.b1_model = B1Model::Li6Convolution;
    CHECK_NOTHROW(c.validate());
    c.b1_band_scale = 0.0;
    c.b1_alpha_d_dwave_weight = 2.0;
    CHECK_NOTHROW(c.validate());
  }
  // ... and the default is legal everywhere it is today
  {
    PipelineConfig c = base();
    CHECK(c.b1_model == B1Model::Miller);
    CHECK(c.b1_band_scale == 1.0);
    CHECK(c.b1_alpha_d_dwave_weight == 1.0);
    CHECK_NOTHROW(c.validate());
    c.channel = PipelineChannel::TaggedLi6Alpha;
    CHECK_NOTHROW(c.validate());
  }

  // ---- 7Li: spin 3/2, no rank-2 input here
  {
    PipelineConfig c = base();
    c.isotope = "7Li";
    c.b1_model = B1Model::Li6Convolution;
    CHECK_THROWS_AS(c.validate(), std::runtime_error);
  }
  // ---- the deuteron: SPIN 1, and that is exactly the trap.  BOTH opt-in
  //      models are refused: `Cdks` is `Li6B1(CdksB1)`, i.e. the 6Li rank-2
  //      transfer (2/6 x LI6_B1_RANK2_TRANSFER), which on a deuteron beam is
  //      the same category error.
  {
    PipelineConfig c = base();
    c.isotope = "d";
    c.b1_model = B1Model::Li6Convolution;
    CHECK_THROWS_AS(c.validate(), std::runtime_error);
    c.b1_model = B1Model::Cdks;
    CHECK_THROWS_AS(c.validate(), std::runtime_error);
    // ... and the DEFAULT is of course still legal there
    c.b1_model = B1Model::Miller;
    CHECK_NOTHROW(c.validate());
  }
  // ---- 7Li with the CDKS camp: spin 3/2 has no b1 slot at all, so the flag
  //      would reach the metadata and nothing else
  {
    PipelineConfig c = base();
    c.isotope = "7Li";
    c.b1_model = B1Model::Cdks;
    CHECK_THROWS_AS(c.validate(), std::runtime_error);
  }
  // ---- tagged and coherent channels, for BOTH opt-in models
  for (B1Model model : {B1Model::Li6Convolution, B1Model::Cdks}) {
    PipelineConfig c = base();
    c.b1_model = model;
    c.channel = PipelineChannel::TaggedLi6Alpha;
    CHECK_THROWS_AS(c.validate(), std::runtime_error);
    c.channel = PipelineChannel::CoherentLi6;
    CHECK_THROWS_AS(c.validate(), std::runtime_error);
    c.isotope = "d";
    c.channel = PipelineChannel::TaggedDeuteronP;
    CHECK_THROWS_AS(c.validate(), std::runtime_error);
  }
  // ---- a caller-supplied kernel and a non-default b1_model contradict
  {
    PipelineConfig c = base();
    c.kernel = default_inclusive_kernel(LI6());
    CHECK_NOTHROW(c.validate());              // ... with the DEFAULT model
    c.b1_model = B1Model::Cdks;
    CHECK_THROWS_AS(c.validate(), std::runtime_error);
    c.b1_model = B1Model::Li6Convolution;
    CHECK_THROWS_AS(c.validate(), std::runtime_error);
  }
  // ---- the band scale is a scale, on EVERY model (the >= 0 check used to
  //      sit inside the non-Miller branch, so `miller` accepted -1)
  for (B1Model model : {B1Model::Miller, B1Model::Cdks,
                        B1Model::Li6Convolution}) {
    PipelineConfig c = base();
    c.b1_model = model;
    c.b1_band_scale = -1.0;
    CHECK_THROWS_AS(c.validate(), std::runtime_error);
    c.b1_band_scale = 1.0;
    c.b1_alpha_d_dwave_weight = -5.0;
    CHECK_THROWS_AS(c.validate(), std::runtime_error);
  }
  // ---- A KNOB THAT DID NOT RUN MAY NOT BE RECORDED AS IF IT HAD.  Neither
  //      scale reaches the Miller branch of `default_inclusive_kernel`, and
  //      the Cdks branch never reads the alpha-d weight -- but the npz/HFS
  //      `meta` records all three unconditionally, so a run with either set
  //      would claim a variation it never made.  That is the false provenance
  //      the three meta keys exist to prevent (design 4.2).
  {
    PipelineConfig c = base();
    c.b1_band_scale = 2.0;                       // miller + band
    CHECK_THROWS_AS(c.validate(), std::runtime_error);
    c.b1_band_scale = 0.0;
    CHECK_THROWS_AS(c.validate(), std::runtime_error);
    c.b1_band_scale = 1.0;
    c.b1_alpha_d_dwave_weight = 2.0;             // miller + alpha-d weight
    CHECK_THROWS_AS(c.validate(), std::runtime_error);
    c.b1_model = B1Model::Cdks;                  // cdks + alpha-d weight
    CHECK_THROWS_AS(c.validate(), std::runtime_error);
    c.b1_model = B1Model::Li6Convolution;        // ... and it is legal HERE
    CHECK_NOTHROW(c.validate());
  }
  // ---- and the names the CLI prints
  CHECK(std::string(b1_model_name(B1Model::Miller)) == "miller");
  CHECK(std::string(b1_model_name(B1Model::Cdks)) == "cdks");
  CHECK(std::string(b1_model_name(B1Model::Li6Convolution)) ==
        "li6-convolution");
}

// ==================================================================== T16
// `--b1-unpol`: THE GATE-PASSING CONFIGURATION HAS TO BE EMITTABLE.
//
// The A = 2 gate passes at G3b = 0.843 with CDKS's own MSTW2008 LO
// (checklist item 4 above) and at 0.440 on the library `ToyF2` (layer 3), so
// the pass is a statement about a CONFIGURATION.  Until `B1UnpolSource`
// existed the run surface could not produce that configuration:
// `default_inclusive_kernel` hard-wired `ToyF2` into
// `Li6ConvolutionOptions::unpol`, and the only other route -- hand-building a
// kernel and putting it in `PipelineConfig::kernel` -- is refused together
// with a non-Miller `b1_model` by `validate()`.  Every 6Li b1 number the
// generator could emit was therefore made at 0.440.
//
// THE KERNEL/FLAG GUARD IS NOT TOUCHED HERE, and this test pins that too: a
// caller-supplied kernel plus `b1_model = li6-convolution` still throws.  It
// does not need to be loosened, because the selector works THROUGH
// `default_inclusive_kernel` -- the passing configuration no longer needs a
// hand-built kernel to be expressed.
TEST_CASE("b1_nuclear T16: --b1-unpol reaches the kernel and the default "
          "does not move") {
  const double q2 = 2.5;
  const std::array<double, 3> xs = {0.10, 0.30, 0.50};

  // ---- the names the CLI and `meta["b1_unpol"]` print
  CHECK(std::string(b1_unpol_name(B1UnpolSource::Toy)) == "toy");
  CHECK(std::string(b1_unpol_name(B1UnpolSource::Mstw)) == "mstw");
  CHECK(std::string(b1_unpol_name(B1UnpolSource::Ct18Nlo)) == "ct18nlo");
  CHECK(std::string(b1_unpol_name(B1UnpolSource::Custom)) == "custom");

  // ---- THE DEFAULT DID NOT MOVE, BIT FOR BIT.  Two clauses: the new
  //      argument defaulted away is the old four-argument call, and an
  //      EXPLICIT null is the same object again.  (T9 above pins the Miller
  //      default against validation/reference/b1_default_li6.json at rtol
  //      1e-12; this pins the li6-convolution branch, which has no reference
  //      file, against itself.)
  {
    const auto a = default_inclusive_kernel(LI6(), B1Model::Li6Convolution,
                                            1.0, 1.0);
    const auto b = default_inclusive_kernel(LI6(), B1Model::Li6Convolution,
                                            1.0, 1.0, nullptr);
    REQUIRE(a);
    REQUIRE(b);
    for (double x : xs) {
      const SFTables ta = a->tables(x, q2);
      const SFTables tb = b->tables(x, q2);
      CAPTURE(x);
      CHECK(ta.b1 == tb.b1);          // bit for bit, not merely close
      CHECK(ta.f1 == tb.f1);
      CHECK(ta.delta == tb.delta);
    }
    // ... and the config default is `Toy` with an EMPTY slot, which is what
    // makes `default_inclusive_kernel(..., cfg.b1_unpol_sf)` the old call.
    PipelineConfig c;
    CHECK(c.b1_unpol == B1UnpolSource::Toy);
    CHECK(!c.b1_unpol_sf);
  }

  // ---- validate(): the provenance rules, and the guard that STAYS
  {
    auto base = [] {
      PipelineConfig c;
      c.channel = PipelineChannel::Inclusive;
      c.isotope = "6Li";
      c.n_events = 100;
      c.b1_model = B1Model::Li6Convolution;
      return c;
    };
    const auto toy = std::make_shared<const ToyF2>();   // stands for any UnpolSF

    // a named backend with an EMPTY slot is an error, never a fallback
    for (B1UnpolSource s : {B1UnpolSource::Mstw, B1UnpolSource::Ct18Nlo,
                            B1UnpolSource::Custom}) {
      PipelineConfig c = base();
      c.b1_unpol = s;
      CHECK_THROWS_AS(c.validate(), std::runtime_error);
      c.b1_unpol_sf = toy;                    // ... and legal once attached
      CHECK_NOTHROW(c.validate());
    }
    // `Toy` with an object attached is the same contradiction the other way
    {
      PipelineConfig c = base();
      c.b1_unpol_sf = toy;
      CHECK_THROWS_AS(c.validate(), std::runtime_error);
    }
    // A KNOB THAT DID NOT RUN MAY NOT BE RECORDED AS IF IT HAD: only the
    // li6-convolution branch reads `unpol`, but meta["b1_unpol"] is written
    // unconditionally.
    for (B1Model m : {B1Model::Miller, B1Model::Cdks}) {
      PipelineConfig c = base();
      c.b1_model = m;
      c.b1_unpol = B1UnpolSource::Custom;
      c.b1_unpol_sf = toy;
      CHECK_THROWS_AS(c.validate(), std::runtime_error);
      c.b1_unpol = B1UnpolSource::Toy;        // ... and the default is legal
      c.b1_unpol_sf.reset();
      CHECK_NOTHROW(c.validate());
    }
    // THE KERNEL/FLAG GUARD IS UNCHANGED.  It is not loosened to let the
    // gate-passing configuration through -- the selector routes around the
    // need for it.
    {
      PipelineConfig c = base();
      c.kernel = default_inclusive_kernel(LI6());
      c.b1_unpol = B1UnpolSource::Custom;
      c.b1_unpol_sf = toy;
      CHECK_THROWS_AS(c.validate(), std::runtime_error);
    }
  }

#ifdef LIPOLGEN_HAVE_PYTHIA8
  // ---- MSTW2008 LO through the PIPELINE'S OWN kernel construction: the
  //      configuration the gate passes on, now emittable.
  if (!mstw_grid_present()) {
    MESSAGE("SKIPPED (MSTW2008 LO unavailable): no mstw2008lo.00.dat under "
            << pythia8_pdfdata_dir()
            << " -- T16's MSTW row did NOT run, so this build did not check "
               "that --b1-unpol mstw reaches the kernel.  It also cannot "
               "SELECT it: PipelineConfig::validate() refuses the flag here "
               "rather than downgrading to the toy.");
  } else {
    const auto toyk = default_inclusive_kernel(LI6(), B1Model::Li6Convolution,
                                               1.0, 1.0);
    const auto mstwk = default_inclusive_kernel(
        LI6(), B1Model::Li6Convolution, 1.0, 1.0, std::make_shared<MstwSF>());
    REQUIRE(mstwk);
    // Measured 2026-09-04.  Up to a factor 1.85, and NOT monotone in x -- so
    // the choice of nucleon PDF is a shape change in the shipped 6Li
    // observable and not a normalisation that a band would absorb.
    const std::array<double, 3> want = {1.847766, 1.275961, 0.816971};
    for (std::size_t i = 0; i < xs.size(); ++i) {
      const SFTables tt = toyk->tables(xs[i], q2);
      const SFTables tm = mstwk->tables(xs[i], q2);
      CAPTURE(xs[i]);
      CHECK(tt.b1 != 0.0);
      CHECK_CLOSE(tm.b1 / tt.b1, want[i], 1e-5);
      // F1 DOES NOT MOVE.  `opt.f2_source` is set by a SEPARATE selector
      // (`--unpol-sf`, left at its `toy` default here), so under THIS flag
      // F1 -- and with it the spin-blind cell cross section and the D_phi
      // denominator of the tensor weight -- is bit-identical; what moves is
      // the tensor shift, which is the point.
      // The price is that the numerator's F1 and the denominator's F1 are no
      // longer the same object; that is documented at `B1UnpolSource`.
      CHECK(tm.f1 == tt.f1);
      CHECK(tm.delta == tt.delta);
      MESSAGE("T16 x = " << xs[i] << ": b1 toy " << tt.b1 << " -> mstw "
                         << tm.b1 << "  (x" << tm.b1 / tt.b1 << ")");
    }
  }
#endif  // LIPOLGEN_HAVE_PYTHIA8
}

// ===================================================================== T17
// G3a's CEILING, pinned.  The clause counts sign changes in (0, 1.0]; this
// case measures what that ceiling excludes.
//
// The 2026-09-03 amendment widened G3a's counting FLOOR on the argument that
// "bookkeeping that brackets a physics tolerance has to be at least as wide
// as the tolerance, or it is a second, tighter, unstated cut".  The ceiling
// was left at 1.0.  Over the digitized reference's OWN domain [0.010, 1.590]
// the reference has TWO sign changes and stays positive from 0.4572 to its
// last point, while EVERY computed configuration crosses zero a THIRD time.
// "Exactly two" is therefore true of a windowed curve, and the window's upper
// edge is a scope choice that no position tolerance derives.  x > 1 per
// nucleon is kinematically allowed for a nucleus (for A = 2, x runs to 2), so
// the crossing is in a physically meaningful region and cannot be dismissed
// as out of range.
//
// This case exists so that the exclusion is a MEASUREMENT in the suite and
// not a footnote: if a change moves the third crossing, or removes it, or
// the reference acquires one, the suite says so.  The clause itself is NOT
// changed here -- see design_D_b1_li6.md's amendment and
// docs/OPEN_ITEMS_SOLUTIONS.md sec. 10 ("G3a's stated limitation") for the
// argument and for why the status is recorded as a QUALIFIED pass.
//
// Grid: the x >= 0.9 sub-grid of the gate's own `linspace(0.001, 1.59, 300)`.
// `b1_landmarks` locates a zero by linear interpolation between CONSECUTIVE
// grid points, so dropping the points below 0.9 -- all of them far from this
// crossing -- leaves every pair above it intact and the position identical to
// the full-scan one.  Cheap on purpose: three columns, not thirty.
TEST_CASE("b1_nuclear T17: G3a's counting ceiling -- the THIRD sign change") {
  REQUIRE(have(kFdeut));
  const std::vector<double> full = linspace(0.001, 1.59, 300);
  std::vector<double> xs;
  for (double x : full)
    if (x >= 0.9) xs.push_back(x);
  REQUIRE(xs.size() == 130);

  // The reference, over its own domain: TWO zeros, and positive at the top.
  const B1Landmarks ref = b1_landmarks_of_table();
  REQUIRE(ref.zeros.size() == 2);
  CHECK(ref.zero_slope[0] == -1);
  CHECK(ref.zero_slope[1] == +1);
  const std::shared_ptr<const TensorSF> raw = cdks_b1_raw_per_nucleon();
  CHECK(1.59 * raw->b1(1.59, kQ2, 0.0) > 0.0);
  // ... and NOT at its noise floor where the AV18 rows cross: the digitized
  // x*b1 is still 7.0 % of its own peak at x = 1.2204.  That is the number
  // that kills the "the reference is unreadable up there" defence of the
  // ceiling for those rows; it only becomes an honest defence near x = 1.5,
  // where the reference is 0.8 % of its peak (the CD-Bonn rows).
  const double at_toy_zero = 1.220437 * raw->b1(1.220437, kQ2, 0.0);
  CHECK_CLOSE(at_toy_zero, 7.606130e-5, 1e-3);
  CHECK_CLOSE(at_toy_zero / ref.xb1_max, 0.070089, 2e-3);

  const auto third = [&](const std::string& who,
                         const std::vector<double>& xb1) {
    std::map<double, double> tab;
    for (std::size_t i = 0; i < xs.size(); ++i) tab[xs[i]] = xb1[i];
    const B1Landmarks lm =
        b1_landmarks([&](double x) { return tab.at(x); }, xs);
    MESSAGE("T17 " << who << ": zeros above x = 0.9 = " << lm.zeros.size()
                   << " ; x*b1(0.9) = " << xb1.front()
                   << " ; x*b1(1.59) = " << xb1.back());
    REQUIRE(lm.zeros.size() == 1);       // the THIRD of the full scan
    CHECK(lm.zero_slope[0] == -1);       // falling: positive -> negative
    CHECK(xb1.front() > 0.0);            // positive where G3a stops counting
    CHECK(xb1.back() < 0.0);             // and negative at the top, unlike
                                         // the reference, which is positive
    return lm.zeros[0];
  };

  // ToyF2 -- the SHIPPED default, so this case is never vacuous.
  CHECK_CLOSE(third("ToyF2, kappa [Eq. 21], AV18",
                    xb1_on_gate_threaded(xs)), 1.220437, 2e-4);

#ifdef LIPOLGEN_HAVE_PYTHIA8
  if (!mstw_grid_present()) {
    MESSAGE("SKIPPED (MSTW2008 LO unavailable): no mstw2008lo.00.dat under "
            << pythia8_pdfdata_dir()
            << " -- the third crossing of the VERDICT row was not measured "
               "in this build.");
  } else {
    const unsigned nt = gate_threads();
    std::vector<DeuteronConvolutionB1::Options> opts(nt);
    for (unsigned t = 0; t < nt; ++t) opts[t].unpol = std::make_shared<MstwSF>();
    // The verdict row -- the configuration the 6Li ban was lifted on -- has
    // it too, and closer to the reference's readable region than the toy.
    CHECK_CLOSE(third("MSTW2008 LO, kappa [Eq. 21], AV18",
                      xb1_threaded(xs, opts)), 1.217660, 2e-4);
  }
#endif  // LIPOLGEN_HAVE_PYTHIA8

#ifdef LIPOLGEN_HAVE_LHAPDF
  lhapdf_quiet();
  {
    const unsigned nt = gate_threads();
    std::vector<DeuteronConvolutionB1::Options> opts(nt);
    for (unsigned t = 0; t < nt; ++t)
      opts[t].unpol = std::make_shared<LhapdfSF>("CT18NLO", 0);
    CHECK_CLOSE(third("CT18NLO, kappa [Eq. 21], AV18",
                      xb1_threaded(xs, opts)), 1.197722, 2e-4);
  }
#endif  // LIPOLGEN_HAVE_LHAPDF
}

TEST_CASE("b1_nuclear D1/F4: a half-assigned ClusterPartialWave throws "
          "instead of walking off the end of phi") {
  // `k` and `phi` are two public vectors that Python assigns ONE AT A TIME,
  // so between the assignments they have different lengths; `rebuild()` used
  // to build the spline anyway and the constructor indexed `y_[i + 1]` over
  // an empty `phi` -- a segmentation fault from three lines of documented
  // API (docs/open_items/run_2026-09-03/phase_D_li7_rank2.md sec. 1.5, F4).
  ClusterPartialWave w;
  w.k = {0.1, 0.2, 0.3};
  w.rebuild();                       // DEFERS, so either order works
  CHECK(!w.spline);
  CHECK_THROWS_AS(w(0.15), std::runtime_error);
  CHECK_THROWS_AS(w.norm2(), std::runtime_error);
  std::vector<double> out;
  CHECK_THROWS_AS(w.eval_sorted({0.15}, &out), std::runtime_error);
  w.phi = {1.0, 2.0, 1.0};
  w.rebuild();
  CHECK(!!w.spline);
  CHECK_CLOSE_AT(w(0.15), 1.6875, 0.0, 1e-12);
  // the other assignment order, which is the one a reader writes first
  ClusterPartialWave v;
  v.phi = {1.0, 2.0, 1.0};
  v.rebuild();
  CHECK_THROWS_AS(v(0.15), std::runtime_error);
  v.k = {0.1, 0.2, 0.3};
  v.rebuild();
  CHECK_CLOSE_AT(v(0.15), 1.6875, 0.0, 1e-12);
}

// -----------------------------------------------------------------------
// T18: `--r-source`, the registry's option (iii) -- ONE R hook threaded into
// BOTH halves of the 6Li tensor weight.
//
// WHAT THIS PINS.  The shipped observable is a RATIO, K/D_phi (`azz`,
// asymmetries.hpp): its numerator's R is `Li6ConvolutionOptions::r_func` and
// its denominator's is `InclusiveKernel::Options::r_func`.  The decision of
// 2026-09-03 (STATUS.md row 3) kept the numerator on `r_sigma_lt` BECAUSE
// moving it alone does not cancel against a denominator still on
// `r_sigma_lt`, and named the fix it did not make: one hook into both.
// `RSource` is that hook.  Three clauses, one per status of the axis:
//
//   Unset     does not enter the wiring branch -- bit for bit, and the two
//             `default_inclusive_kernel` overloads agree;
//   SigmaLt   installs the SHARED object and is bit-identical anyway,
//             because `resolve_r`'s null branch IS `r_sigma_lt`.  This is
//             the wiring's own test: it says the hook reaches both halves;
//   R1998     moves b1 exactly as the numerator-only option (ii) does, AND
//             moves F1 -- which is what (ii) cannot do.
//
// The measured cost is in docs/open_items/run_2026-09-06/phase_A_numbers.md
// sec. A1; python/tests/test_b1_model.py P10 pins the numbers, this case
// pins the WIRING.
TEST_CASE("b1_nuclear T18: --r-source puts ONE R in numerator and denominator") {
  const double q2 = 2.5;
  const std::array<double, 3> xs = {0.05, 0.10, 0.30};

  CHECK(std::string(r_source_name(RSource::Unset)) == "unset");
  CHECK(std::string(r_source_name(RSource::SigmaLt)) == "sigma-lt");
  CHECK(std::string(r_source_name(RSource::R1998)) == "r1998");
  CHECK(PipelineConfig().r_source == RSource::Unset);

  const auto base = default_inclusive_kernel(LI6(), B1Model::Li6Convolution);
  const auto unset = default_inclusive_kernel(
      LI6(), B1Model::Li6Convolution, 1.0, 1.0, nullptr, nullptr, nullptr,
      RSource::Unset);
  const auto shared_lt = default_inclusive_kernel(
      LI6(), B1Model::Li6Convolution, 1.0, 1.0, nullptr, nullptr, nullptr,
      RSource::SigmaLt);
  const auto shared_98 = default_inclusive_kernel(
      LI6(), B1Model::Li6Convolution, 1.0, 1.0, nullptr, nullptr, nullptr,
      RSource::R1998);
  REQUIRE(base);
  REQUIRE(unset);
  REQUIRE(shared_lt);
  REQUIRE(shared_98);

  // The NUMERATOR-ONLY option (ii), hand-built: no flag reaches it and none
  // is being added.  Everything but `o.r_func` is the Li6Convolution branch
  // of `default_inclusive_kernel` verbatim, including the ONE shared ToyF2.
  const auto f2 = std::make_shared<const ToyF2>();
  Li6ConvolutionOptions o;
  o.unpol = f2;
  o.r_func = [](double x, double qq) { return r1998(x, qq); };
  const auto num_only_b1 = std::make_shared<const Li6ConvolutionB1>(std::move(o));
  InclusiveKernel::Options nopt;
  nopt.f2_source = f2;
  nopt.b1_func = [num_only_b1](double x, double qq, double f1) {
    return num_only_b1->b1(x, qq, f1);
  };
  nopt.delta_func = [](double x, double qq, double f1) {
    return toy_delta_gluon(x, qq, f1, 1e-2);
  };
  const InclusiveKernel num_only(LI6(), nopt);

  for (double x : xs) {
    CAPTURE(x);
    const SFTables tb = base->tables(x, q2);
    const SFTables tu = unset->tables(x, q2);
    const SFTables tl = shared_lt->tables(x, q2);
    const SFTables ts = shared_98->tables(x, q2);
    const SFTables tn = num_only.tables(x, q2);

    // (1) the default did not move, on either overload
    CHECK(tb.b1 == tu.b1);
    CHECK(tb.f1 == tu.f1);
    CHECK(tb.f2 == tu.f2);
    CHECK(tb.delta == tu.delta);
    // (2) the shared r_sigma_lt hook is installed and changes nothing
    CHECK(tb.b1 == tl.b1);
    CHECK(tb.f1 == tl.f1);
    CHECK(tb.b2 == tl.b2);
    CHECK(tb.delta == tl.delta);
    // (3) r1998 moves the NUMERATOR exactly as (ii) does -- same object,
    //     same output, bit for bit -- and the DENOMINATOR as well, which is
    //     the whole difference between the two options.
    CHECK(ts.b1 == tn.b1);
    CHECK(ts.b2 == tn.b2);
    CHECK(ts.b1 != tb.b1);
    CHECK(tn.f1 == tb.f1);            // (ii) leaves F1 alone
    CHECK(ts.f1 != tb.f1);            // (iii) does not
    CHECK(ts.f2 == tb.f2);            // F2 has no R in it, on either option
    // ... and F1 moves by exactly the (1 + R) ratio the registry quoted.
    CHECK_CLOSE_AT(tb.f1 / ts.f1,
                   (1.0 + r1998(x, q2)) / (1.0 + r_sigma_lt(x, q2)),
                   0.0, 1e-12);
  }

  // (4) the provenance rule: refused where the hook is not installed.
  {
    PipelineConfig c;
    c.channel = PipelineChannel::Inclusive;
    c.isotope = "6Li";
    c.n_events = 100;
    c.b1_model = B1Model::Li6Convolution;
    c.r_source = RSource::R1998;
    CHECK_NOTHROW(c.validate());
    for (B1Model m : {B1Model::Miller, B1Model::Cdks}) {
      c.b1_model = m;
      CHECK_THROWS_AS(c.validate(), std::runtime_error);
      c.r_source = RSource::Unset;         // ... and the default is legal
      CHECK_NOTHROW(c.validate());
      c.r_source = RSource::R1998;
    }
  }
}
