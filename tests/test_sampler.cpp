// SPDX-License-Identifier: GPL-3.0-or-later
// Sampler-level checks (evgen/tests/test_sampler.py) and the Gate-3
// pseudo-experiment estimator closure (evgen/tests/test_pseudoexp.py):
// distributions, rates, reproducibility, Mode-W weights, the exact positivity
// guard, and pulls / spreads / relative-luminosity biases of the analysis
// estimators.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

#include "check_close.hpp"
#include "lipolgen/asymmetries.hpp"
#include "lipolgen/beams.hpp"
#include "lipolgen/bookkeeping.hpp"
#include "lipolgen/sampler.hpp"
#include "lipolgen/sf.hpp"
#include "lipolgen/xsec.hpp"

using namespace lipolgen;

namespace {

const SpinCategory& unpol_category() {
  static const SpinCategory c("unpol", 1.0,
                              {1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0});
  return c;
}

BeamConfig mid_6li() { return default_configs("6Li")[1]; }

std::shared_ptr<const InclusiveKernel> make_kernel(bool b1, double delta_scale) {
  InclusiveKernel::Options opt;
  if (b1) {
    opt.b1_func = [](double x, double q2, double f1) {
      return toy_b1(x, q2, f1);
    };
  }
  if (delta_scale != 0.0) {
    opt.delta_func = [delta_scale](double x, double q2, double f1) {
      return toy_delta_gluon(x, q2, f1, delta_scale);
    };
  }
  return std::make_shared<InclusiveKernel>(LI6(), opt);
}

InclusiveSampler::GridSpec spec_of(int nx, int nq2, double xlo, double xhi,
                                   double qlo, double qhi) {
  InclusiveSampler::GridSpec g;
  g.nx = nx; g.nq2 = nq2;
  g.x_min = xlo; g.x_max = xhi;
  g.q2_min = qlo; g.q2_max = qhi;
  return g;
}

/// The module-scope fixture of test_sampler.py.
const InclusiveSampler& shared_sampler() {
  static const InclusiveSampler s(make_kernel(true, 1e-2), mid_6li(),
                                  Scenario(),
                                  spec_of(24, 18, 1e-3, 0.5, 1.5, 100.0));
  return s;
}

double sigma_weighted(const InclusiveSampler& s,
                      const std::vector<double>& values) {
  double num = 0.0, den = 0.0;
  for (std::size_t i = 0; i < values.size(); ++i) {
    num += s.cell_xsec_pb()[i] * values[i];
    den += s.cell_xsec_pb()[i];
  }
  return num / den;
}

double mean_of(const std::vector<double>& v) {
  return std::accumulate(v.begin(), v.end(), 0.0) /
         static_cast<double>(v.size());
}

double stddev_of(const std::vector<double>& v) {
  const double mu = mean_of(v);
  double s = 0.0;
  for (double a : v) s += (a - mu) * (a - mu);
  return std::sqrt(s / static_cast<double>(v.size() - 1));
}

}  // namespace

TEST_CASE("sampler: unpolarized cell distribution follows the cell rates") {
  const InclusiveSampler& s = shared_sampler();
  const std::size_t n = 150000;
  const EventBatch ev = s.sample_n(unpol_category(), n, 1, 0, 0);
  REQUIRE(ev.size() == n);
  std::vector<double> counts(s.n_cells(), 0.0);
  for (int c : ev.cell) counts[static_cast<std::size_t>(c)] += 1.0;
  double xsum = 0.0;
  for (double v : s.cell_xsec_pb()) xsum += v;

  // population-averaged rates = unpolarized cell cross sections
  double zmax = 0.0, z2 = 0.0;
  std::size_t used = 0;
  for (std::size_t i = 0; i < s.n_cells(); ++i) {
    const double expected =
        static_cast<double>(n) * s.cell_xsec_pb()[i] / xsum;
    if (expected <= 50.0) continue;
    const double z = (counts[i] - expected) / std::sqrt(expected);
    zmax = std::max(zmax, std::fabs(z));
    z2 += z * z;
    ++used;
  }
  REQUIRE(used > 20);
  CHECK(zmax < 5.0);
  CHECK_CLOSE_AT(std::sqrt(z2 / static_cast<double>(used)), 1.0, 0.0, 0.25);
}

TEST_CASE("sampler: event kinematics are internally consistent") {
  const InclusiveSampler& s = shared_sampler();
  const EventBatch ev = s.sample_n(unpol_category(), 20000, 2, 0, 0);
  for (std::size_t i = 0; i < ev.size(); ++i) {
    CHECK_CLOSE(ev.y[i], ev.q2[i] / (s.s() * ev.x[i]), 1e-12);
    CHECK(ev.phi[i] >= 0.0);
    CHECK(ev.phi[i] < 2.0 * kPi);
    const std::size_t c = static_cast<std::size_t>(ev.cell[i]);
    // events stay inside their cells
    CHECK(ev.x[i] >= std::exp(s.logx_lo()[c]));
    CHECK(ev.x[i] <= std::exp(s.logx_hi()[c]));
    CHECK(ev.q2[i] >= std::exp(s.logq2_lo()[c]));
    CHECK(ev.q2[i] <= std::exp(s.logq2_hi()[c]));
  }
}

TEST_CASE("sampler: sampled events respect the acceptance window") {
  // 2026-08-11 audit: boundary cells must not emit y > 1 or out-of-window
  // events -- log-uniform in-cell placement is redrawn into acceptance.
  const InclusiveSampler& s = shared_sampler();
  const EventBatch ev = s.sample_n(unpol_category(), 200000, 11, 0, 0);
  double ymin = 1e9, ymax = -1e9;
  for (double v : ev.y) { ymin = std::min(ymin, v); ymax = std::max(ymax, v); }
  CHECK(ymax <= s.scenario().y_max + 1e-12);
  CHECK(ymin >= s.scenario().y_min - 1e-12);
}

TEST_CASE("sampler: pure-state cross sections reproduce the analytic Azz") {
  const InclusiveSampler& s = shared_sampler();
  std::vector<double> azz_cells(s.n_cells());
  for (std::size_t i = 0; i < s.n_cells(); ++i) {
    const SFTables& t = s.tables()[i];
    const double y = y_from_xq2(s.x_cells()[i], s.q2_cells()[i], s.s());
    azz_cells[i] = azz(t.b1, t.f1, t.f2, s.x_cells()[i], y);
  }
  const double truth = sigma_weighted(s, azz_cells);
  const double sp = s.sigma_state_pb(unpol_category(), 1.0);
  const double s0 = s.sigma_state_pb(unpol_category(), 0.0);
  const double sm = s.sigma_state_pb(unpol_category(), -1.0);
  CHECK_CLOSE((sp + sm - 2.0 * s0) / (sp + sm + s0), truth, 1e-9);
}

TEST_CASE("sampler: Poisson counts normalize to the expected rate") {
  const InclusiveSampler& s = shared_sampler();
  const double lumi_pb = 1e-4 * s.scenario().lumi_fb_per_nucleon * 1e3;
  const double mu = s.expected_events(unpol_category(), lumi_pb);
  double sum = 0.0;
  const int trials = 60;
  for (int t = 0; t < trials; ++t) {
    sum += static_cast<double>(
        s.sample_lumi(unpol_category(), lumi_pb, 3, static_cast<std::uint64_t>(t),
                      0).size());
  }
  CHECK_CLOSE_AT(sum / trials, mu, 0.0, 4.0 * std::sqrt(mu / trials));
}

TEST_CASE("sampler: (seed, run, bunch) fixes the events") {
  const InclusiveSampler& s = shared_sampler();
  const EventBatch a = s.sample_n(unpol_category(), 5000, 7, 1, 42);
  const EventBatch b = s.sample_n(unpol_category(), 5000, 7, 1, 42);
  const EventBatch c = s.sample_n(unpol_category(), 5000, 7, 1, 43);
  CHECK(a.x == b.x);
  CHECK(a.phi == b.phi);
  CHECK(a.m == b.m);
  CHECK(a.cell == b.cell);
  CHECK(a.x != c.x);
}

TEST_CASE("sampler: vector rate sign and the m = +J..-J labels") {
  // Deterministic guard on the population <-> m pairing: category rates must
  // equal the kernel evaluated at EXPLICITLY constructed spin states,
  // including the SIGN of the vector-sector rate shift.  A mutation that
  // reverses the m-ordering convention passes every statistical closure test
  // -- the A_par truths are smaller than the MC tolerances -- but fails here
  // exactly.
  const InclusiveSampler& s = shared_sampler();
  const SpinCategory plus("v+", 1.0, {1.0, 0.0, 0.0}, +1, 1.0);
  const SpinCategory minus("v-", 1.0, {0.0, 0.0, 1.0}, +1, 1.0);
  double expect_plus = 0.0, expect_minus = 0.0, apar_w = 0.0;
  for (std::size_t i = 0; i < s.n_cells(); ++i) {
    const double x = s.x_cells()[i], q2 = s.q2_cells()[i];
    const double y = y_from_xq2(x, q2, s.s());
    EventSpinState st;
    st.lam_e = +1; st.pe = 1.0; st.j = 1.0;
    st.m = 1.0;
    expect_plus += s.cell_xsec_pb()[i] *
                   (1.0 + s.kernel().amplitudes(s.tables()[i], x, q2, s.s(),
                                                st).w_avg);
    st.m = -1.0;
    expect_minus += s.cell_xsec_pb()[i] *
                    (1.0 + s.kernel().amplitudes(s.tables()[i], x, q2, s.s(),
                                                 st).w_avg);
    apar_w += s.cell_xsec_pb()[i] *
              a_parallel(s.tables()[i].g1, s.tables()[i].f1, y, x, q2);
  }
  CHECK_CLOSE(s.sigma_tot_pb(plus), expect_plus, 1e-12);
  CHECK_CLOSE(s.sigma_tot_pb(minus), expect_minus, 1e-12);
  CHECK(expect_plus != expect_minus);
  CHECK(((expect_plus - expect_minus) > 0.0) == (apar_w > 0.0));

  // event spin labels carry the m = +J ... -J convention
  const EventBatch ev = s.sample_n(plus, 2000, 9, 0, 0);
  for (double m : ev.m) CHECK(m == 1.0);
}

TEST_CASE("sampler: Mode-W weights reproduce the polarized rate ratio") {
  const InclusiveSampler& s = shared_sampler();
  const EventBatch ev = s.sample_n(unpol_category(), 200000, 4, 0, 0);
  const SpinCategory tensor_cat("t0", 1.0, {0.15, 0.7, 0.15});  // m0-enriched
  const std::vector<double> w = s.weights_for(ev, {tensor_cat});
  const double expected =
      s.sigma_tot_pb(tensor_cat) / s.sigma_tot_pb(unpol_category());
  const double mu = mean_of(w);
  const double se = stddev_of(w) / std::sqrt(static_cast<double>(w.size()));
  CHECK_CLOSE_AT(mu, expected, 0.0, 4.0 * se);
}

TEST_CASE("sampler: the positivity guard is exact and refuses a negative W") {
  const BeamConfig cfg = mid_6li();
  // An oversized Delta makes 1 + a2 cos 2phi' dip below zero; the
  // accept-reject would silently sample max(W, 0), diluting the modulation
  // and skewing the (x, Q2) mixture.
  InclusiveKernel::Options big;
  big.delta_func = [](double, double, double f1) { return 3.0 * f1; };
  const InclusiveSampler bad(std::make_shared<InclusiveKernel>(LI6(), big), cfg,
                             Scenario(), spec_of(12, 9, 1e-4, 1.0, 1.0, 2e3));
  const RunPlan plan = transverse_tensor_plan(0.6);
  CHECK_THROWS_AS(bad.sigma_tot_pb(plan.categories()[0]), std::runtime_error);

  // Nothing in the repository's own scenarios comes near the bound.
  InclusiveKernel::Options prod;
  prod.b1_func = [](double x, double q2, double f1) {
    return toy_b1(x, q2, f1);
  };
  prod.delta_func = [](double, double, double f1) { return -1e-2 * f1; };
  const InclusiveSampler ok(std::make_shared<InclusiveKernel>(LI6(), prod), cfg,
                            Scenario(), spec_of(20, 15, 1e-4, 1.0, 1.0, 2e3));
  const SpinCategory& cat = transverse_tensor_plan(0.6).categories()[0];
  CHECK(ok.state_tables(cat, 1.0).margin > 0.9);
}

TEST_CASE("sampler: phi_histogram_pseudo refuses negative bin means") {
  // This path bypasses the sampler's own guard and feeds the estimator
  // directly, so it needs its own.
  const PhiHistogram h = phi_histogram_pseudo(1e6, 0.5, 36, nullptr, 0.0, false);
  for (double c : h.counts) CHECK(c > 0.0);
  CHECK_THROWS_AS(phi_histogram_pseudo(1e6, 1.5, 36, nullptr, 0.0, false),
                  std::runtime_error);
  CHECK_THROWS_AS(phi_histogram_pseudo(1e6, 0.2, 36, nullptr, 1.5, false),
                  std::runtime_error);
}

TEST_CASE("sampler: a run-plan share moves counts, never cross sections") {
  // When run_share was introduced it divided by the PROGRAMME luminosity
  // while n_events already carried the share, so every cross section came out
  // scaled by it and the share was then applied a second time by the caller's
  // lumi_pb.  Dividing by the EFFECTIVE luminosity is the fix.
  const BeamConfig cfg = mid_6li();
  const auto kern = make_kernel(true, 0.0);
  const auto spec = spec_of(24, 18, 1e-4, 1.0, 0.7, 2e3);
  Scenario quarter_scenario;
  quarter_scenario.run_share = 0.25;
  const InclusiveSampler full(kern, cfg, Scenario(), spec);
  const InclusiveSampler quarter(kern, cfg, quarter_scenario, spec);
  REQUIRE(full.n_cells() == quarter.n_cells());
  for (std::size_t i = 0; i < full.n_cells(); ++i) {
    CHECK_CLOSE(quarter.cell_xsec_pb()[i], full.cell_xsec_pb()[i], 1e-12);
  }
  // ... and the counts the caller asks for do carry it, exactly once
  const SpinCategory& cat = tensor_flip_plan(0.6).categories()[0];
  CHECK_CLOSE(quarter.expected_events(cat, 0.25 * 1.0e4),
              0.25 * full.expected_events(cat, 1.0e4), 1e-12);

  Scenario bad;
  bad.run_share = 0.0;
  CHECK_THROWS_AS(InclusiveSampler(kern, cfg, bad, spec), std::runtime_error);
}

TEST_CASE("sampler: the generator window is looser than the analysis one") {
  Scenario analysis;
  const Scenario gen = generator_scenario(analysis);
  CHECK(gen.q2_min == 0.7);
  CHECK(gen.y_min == 0.004);
  CHECK(gen.y_max == 0.985);
  CHECK(gen.w2_min == 8.0);
  CHECK(gen.eta_min == analysis.eta_min - 0.3);
  CHECK(gen.eta_max == analysis.eta_max + 0.3);
  CHECK(gen.e_prime_min == 0.3);
  CHECK(gen.lumi_fb_per_nucleon == analysis.lumi_fb_per_nucleon);
  // every analysis-accepted point is generator-accepted
  const auto kern = make_kernel(false, 0.0);
  const InclusiveSampler a(kern, mid_6li(), analysis,
                           spec_of(30, 20, 1e-4, 1.0, 0.7, 2e3));
  const InclusiveSampler g(kern, mid_6li(), gen,
                           spec_of(30, 20, 1e-4, 1.0, 0.7, 2e3));
  CHECK(g.n_cells() > a.n_cells());
  for (std::size_t i = 0; i < a.n_cells(); ++i) {
    CHECK(g.in_acceptance(a.x_cells()[i], a.q2_cells()[i]));
  }
}

// ------------------------------------------------ Gate 3: estimator closure

namespace {

/// The nx=10, nq2=8 windows of test_pseudoexp.py.
struct ClosureSetup {
  std::shared_ptr<const InclusiveSampler> sampler;
  double truth = 0.0;
};

const ClosureSetup& vector_setup() {
  static ClosureSetup c = [] {
    ClosureSetup out;
    out.sampler = std::make_shared<InclusiveSampler>(
        make_kernel(false, 0.0), mid_6li(), Scenario(),
        spec_of(10, 8, 0.03, 0.3, 4.0, 60.0));
    std::vector<double> v(out.sampler->n_cells());
    for (std::size_t i = 0; i < v.size(); ++i) {
      const double x = out.sampler->x_cells()[i], q2 = out.sampler->q2_cells()[i];
      const double y = y_from_xq2(x, q2, out.sampler->s());
      v[i] = a_parallel(out.sampler->tables()[i].g1,
                        out.sampler->tables()[i].f1, y, x, q2);
    }
    out.truth = sigma_weighted(*out.sampler, v);
    return out;
  }();
  return c;
}

const ClosureSetup& tensor_setup() {
  static ClosureSetup c = [] {
    ClosureSetup out;
    out.sampler = std::make_shared<InclusiveSampler>(
        make_kernel(true, 0.0), mid_6li(), Scenario(),
        spec_of(10, 8, 0.03, 0.3, 4.0, 60.0));
    std::vector<double> v(out.sampler->n_cells());
    for (std::size_t i = 0; i < v.size(); ++i) {
      const SFTables& t = out.sampler->tables()[i];
      const double x = out.sampler->x_cells()[i], q2 = out.sampler->q2_cells()[i];
      const double y = y_from_xq2(x, q2, out.sampler->s());
      v[i] = azz(t.b1, t.f1, t.f2, x, y);
    }
    out.truth = sigma_weighted(*out.sampler, v);
    return out;
  }();
  return c;
}

const ClosureSetup& cos2phi_setup() {
  static ClosureSetup c = [] {
    ClosureSetup out;
    out.sampler = std::make_shared<InclusiveSampler>(
        make_kernel(false, 0.05), mid_6li(), Scenario(),
        spec_of(10, 8, 0.02, 0.12, 8.0, 80.0));
    std::vector<double> v(out.sampler->n_cells());
    for (std::size_t i = 0; i < v.size(); ++i) {
      const SFTables& t = out.sampler->tables()[i];
      const double x = out.sampler->x_cells()[i], q2 = out.sampler->q2_cells()[i];
      const double y = y_from_xq2(x, q2, out.sampler->s());
      v[i] = a_cos2phi(t.delta, t.f1, t.f2, x, y);
    }
    out.truth = sigma_weighted(*out.sampler, v);
    return out;
  }();
  return c;
}

double lumi_for(const InclusiveSampler& s, const RunPlan& plan,
                double n_target) {
  double sigma = 0.0;
  for (const SpinCategory& c : plan.categories()) {
    sigma += s.sigma_tot_pb(c) * c.lumi_fraction;
  }
  return n_target / sigma;
}

}  // namespace

TEST_CASE("closure: A_par pulls unbiased, spread within 15% of err_a_parallel") {
  const ClosureSetup& setup = vector_setup();
  const InclusiveSampler& s = *setup.sampler;
  const double pe = 0.7, pz = 0.6;
  const RunPlan plan = helicity_flip_plan(1.0, pz, pe);
  const double lumi = lumi_for(s, plan, 60000);
  std::vector<double> vals;
  double nsum = 0.0;
  const int trials = 240;
  for (int t = 0; t < trials; ++t) {
    const PseudoExperiment px =
        run_pseudo_experiment(s, plan, lumi, 11, static_cast<std::uint64_t>(t));
    const double np = static_cast<double>(px.count("apar+"));
    const double nm = static_cast<double>(px.count("apar-"));
    vals.push_back(estimators::apar_flip(np, nm, pe, pz));
    nsum += np + nm;
  }
  const double sd = stddev_of(vals);
  const double expected_err = err_a_parallel(nsum / trials, pe, pz);
  CHECK_CLOSE_AT(mean_of(vals), setup.truth, 0.0,
                 4.0 * sd / std::sqrt(static_cast<double>(trials)));
  CHECK(sd / expected_err > 0.85);
  CHECK(sd / expected_err < 1.15);
}

TEST_CASE("closure: Azz pulls unbiased, spread within 15% of err_azz") {
  const ClosureSetup& setup = tensor_setup();
  const InclusiveSampler& s = *setup.sampler;
  const double pz = 0.7, pzz = 0.6;
  const RunPlan plan = tensor_thirds_plan(pz, pzz);
  const double lumi = lumi_for(s, plan, 50000);
  std::vector<double> vals;
  double nsum = 0.0;
  const int trials = 280;
  for (int t = 0; t < trials; ++t) {
    const PseudoExperiment px =
        run_pseudo_experiment(s, plan, lumi, 12, static_cast<std::uint64_t>(t));
    const double np = static_cast<double>(px.count("azz+"));
    const double nm = static_cast<double>(px.count("azz-"));
    const double n0 = static_cast<double>(px.count("azz0"));
    vals.push_back(estimators::azz_thirds(np, nm, n0, pzz));
    nsum += np + nm + n0;
  }
  const double sd = stddev_of(vals);
  const double expected_err = err_azz(nsum / trials, pzz);
  CHECK_CLOSE_AT(mean_of(vals), setup.truth, 0.0,
                 4.0 * sd / std::sqrt(static_cast<double>(trials)));
  CHECK(sd / expected_err > 0.85);
  CHECK(sd / expected_err < 1.15);
}

TEST_CASE("closure: the Azz relative-luminosity bias is -(2/3) delta / Pzz") {
  const ClosureSetup& setup = tensor_setup();
  const InclusiveSampler& s = *setup.sampler;
  const double pz = 0.7, pzz = 0.6, delta = 0.02;
  const RunPlan plan = tensor_thirds_plan(pz, pzz, delta);
  const double lumi = lumi_for(s, plan, 90000);
  std::vector<double> naive, corrected;
  const int trials = 140;
  for (int t = 0; t < trials; ++t) {
    const PseudoExperiment px =
        run_pseudo_experiment(s, plan, lumi, 13, static_cast<std::uint64_t>(t));
    const double np = static_cast<double>(px.count("azz+"));
    const double nm = static_cast<double>(px.count("azz-"));
    const double n0 = static_cast<double>(px.count("azz0"));
    naive.push_back(estimators::azz_thirds(np, nm, n0, pzz));
    corrected.push_back(estimators::azz_thirds(
        np, nm, n0, pzz,
        {px.lumi_of("azz+"), px.lumi_of("azz-"), px.lumi_of("azz0")}));
  }
  const double se = stddev_of(naive) / std::sqrt(static_cast<double>(trials));
  const double bias = azz_rel_lumi_bias(delta, pzz);
  CHECK(std::fabs(bias) > 10.0 * se);   // the test can resolve the bias
  CHECK_CLOSE_AT(mean_of(naive) - setup.truth, bias, 0.0, 5.0 * se);
  CHECK_CLOSE_AT(mean_of(corrected), setup.truth, 0.0, 5.0 * se);
}

TEST_CASE("closure: the A_par relative-luminosity bias is delta/(2 Pe Pz)") {
  const ClosureSetup& setup = vector_setup();
  const InclusiveSampler& s = *setup.sampler;
  const double pe = 0.7, pz = 0.6, delta = 0.02;
  HelicityFlipOptions opt;
  opt.rel_lumi_offset = delta;
  const RunPlan plan = helicity_flip_plan(1.0, pz, pe, opt);
  const double lumi = lumi_for(s, plan, 80000);
  std::vector<double> naive, corrected;
  const int trials = 140;
  for (int t = 0; t < trials; ++t) {
    const PseudoExperiment px =
        run_pseudo_experiment(s, plan, lumi, 14, static_cast<std::uint64_t>(t));
    const double np = static_cast<double>(px.count("apar+"));
    const double nm = static_cast<double>(px.count("apar-"));
    naive.push_back(estimators::apar_flip(np, nm, pe, pz));
    corrected.push_back(estimators::apar_flip(np, nm, pe, pz,
                                              px.lumi_of("apar+"),
                                              px.lumi_of("apar-")));
  }
  const double se = stddev_of(naive) / std::sqrt(static_cast<double>(trials));
  const double bias = apar_rel_lumi_bias(delta, pe, pz);
  CHECK(std::fabs(bias) > 10.0 * se);
  CHECK_CLOSE_AT(mean_of(naive) - setup.truth, bias, 0.0, 5.0 * se);
  CHECK_CLOSE_AT(mean_of(corrected), setup.truth, 0.0, 5.0 * se);
}

TEST_CASE("closure: cos2phi moment unbiased, spread within 15%") {
  const ClosureSetup& setup = cos2phi_setup();
  const InclusiveSampler& s = *setup.sampler;
  const double pzz = 0.8;
  const RunPlan plan = transverse_tensor_plan(pzz, 0.6);
  const SpinCategory& cat = plan.categories()[0];
  const std::size_t n = 40000;
  std::vector<double> vals;
  const int trials = 260;
  for (int t = 0; t < trials; ++t) {
    const EventBatch ev =
        s.sample_n(cat, n, 15, static_cast<std::uint64_t>(t), 0);
    std::vector<double> phip(ev.size());
    for (std::size_t i = 0; i < ev.size(); ++i) phip[i] = ev.phi[i] - cat.phi_s;
    vals.push_back(estimators::cos2phi_moment(phip, pzz));
  }
  const double sd = stddev_of(vals);
  const double expected_err = err_cos2phi_amplitude(static_cast<double>(n), pzz);
  CHECK_CLOSE_AT(mean_of(vals), setup.truth, 0.0,
                 4.0 * sd / std::sqrt(static_cast<double>(trials)));
  CHECK(sd / expected_err > 0.85);
  CHECK(sd / expected_err < 1.15);
}

TEST_CASE("closure: the binned cos2phi fit survives two dead phi sectors") {
  // Gate: phi-modulation recovery with uniform AND holey acceptance.  The
  // holes are asymmetric in 2 phi, so the naive moment breaks and the fit
  // must not.
  const ClosureSetup& setup = cos2phi_setup();
  const InclusiveSampler& s = *setup.sampler;
  const double pzz = 0.8;
  const RunPlan plan = transverse_tensor_plan(pzz);
  const SpinCategory& cat = plan.categories()[0];
  auto acceptance = [](double p) {
    return !((p >= 0.4 && p <= 1.0) || (p >= 3.3 && p <= 4.1));
  };
  const std::size_t n = 60000;
  std::vector<double> fit_vals, moment_vals;
  const int trials = 120;
  for (int t = 0; t < trials; ++t) {
    const EventBatch ev =
        s.sample_n(cat, n, 16, static_cast<std::uint64_t>(t), 0);
    std::vector<double> kept;
    kept.reserve(ev.size());
    for (std::size_t i = 0; i < ev.size(); ++i) {
      const double p = ev.phi[i] - cat.phi_s;
      double w = std::fmod(p, 2.0 * kPi);
      if (w < 0.0) w += 2.0 * kPi;
      if (acceptance(w)) kept.push_back(p);
    }
    fit_vals.push_back(estimators::cos2phi_fit(kept, pzz, 36, acceptance));
    moment_vals.push_back(estimators::cos2phi_moment(kept, pzz));
  }
  const double se = stddev_of(fit_vals) / std::sqrt(static_cast<double>(trials));
  CHECK_CLOSE_AT(mean_of(fit_vals), setup.truth, 0.0, 5.0 * se);
  // the naive moment IS biased under these holes (sanity of the gate)
  CHECK(std::fabs(mean_of(moment_vals) - setup.truth) > 10.0 * se);
}

TEST_CASE("estimators: nothing downstream divides a flip plan by pzz") {
  // The estimator a helicity-flip plan feeds divides by pe * pz alone, so the
  // spin-1/2 fill's exactly-zero pzz_true can never reach a denominator.
  const RunPlan plan = helicity_flip_plan(0.5, 0.6, 0.7);
  const double a_par = 0.037;
  const double scale = plan.pe_true() * plan.pz_true();
  const double n_plus = 1.0e6 * (1.0 + scale * a_par);
  const double n_minus = 1.0e6 * (1.0 - scale * a_par);
  CHECK_CLOSE(estimators::apar_flip(n_plus, n_minus, plan.pe_true(),
                                    plan.pz_true()),
              a_par, 1e-12);
  // the luminosity-corrected form is the same number at equal shares
  CHECK_CLOSE(estimators::apar_flip(n_plus, n_minus, plan.pe_true(),
                                    plan.pz_true(), 500.0, 500.0),
              a_par, 1e-12);
}

TEST_CASE("closure: the binned fit reproduces an injected amplitude") {
  // Full-statistics path: exact bin means, no sampling, so the only thing
  // under test is the finite-bin dilution correction.
  const double pzz = 0.6, amp = 0.05;
  const PhiHistogram h =
      phi_histogram_pseudo(1e8, pzz * amp, 24, nullptr, 0.0, false);
  CHECK_CLOSE(estimators::cos2phi_fit_binned(h.counts, h.edges, pzz), amp,
              1e-10);
  // ... and the binned error inflates the analytic one by the same dilution
  const double w = kPi / 24.0;
  CHECK_CLOSE(estimators::cos2phi_fit_err(1e6, pzz, 24),
              err_cos2phi_amplitude(1e6, pzz) / (std::sin(2.0 * w) / (2.0 * w)),
              1e-12);
  // too few live bins is refused, not silently fitted
  auto dead = [](double p) { return p < 0.2; };
  CHECK_THROWS_AS(
      estimators::cos2phi_fit_binned(h.counts, h.edges, pzz, dead),
      std::runtime_error);
}

TEST_CASE("sampler: single-core throughput") {
  const InclusiveSampler& s = shared_sampler();
  const std::size_t n = 400000;
  s.sample_n(unpol_category(), 1000, 1, 0, 0);  // warm the caches
  const auto t0 = std::chrono::steady_clock::now();
  const EventBatch ev = s.sample_n(unpol_category(), n, 99, 0, 0);
  const auto t1 = std::chrono::steady_clock::now();
  const double secs = std::chrono::duration<double>(t1 - t0).count();
  const double rate = static_cast<double>(ev.size()) / secs;
  MESSAGE("T0 sampling throughput: " << rate / 1e6 << " Mevents/s single core");
  CHECK(rate > 1e5);   // the P8 performance gate
}
