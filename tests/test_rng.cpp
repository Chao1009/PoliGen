// SPDX-License-Identifier: GPL-3.0-or-later
// The reproducibility discipline of docs/CONVENTIONS.md: a counter-based
// stream keyed by (seed, run, bunch, event), never a global RNG.  Same
// (seed, run, bunch) -> identical events on any thread count.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <thread>
#include <vector>

#include "check_close.hpp"
#include "lipolgen/beams.hpp"
#include "lipolgen/bookkeeping.hpp"
#include "lipolgen/generator.hpp"
#include "lipolgen/rng.hpp"
#include "lipolgen/sampler.hpp"
#include "lipolgen/sf.hpp"
#include "lipolgen/xsec.hpp"

using namespace lipolgen;

namespace {

std::shared_ptr<const InclusiveSampler> rng_test_sampler() {
  static std::shared_ptr<const InclusiveSampler> s = [] {
    InclusiveKernel::Options opt;
    opt.b1_func = [](double x, double q2, double f1) {
      return toy_b1(x, q2, f1);
    };
    opt.delta_func = [](double, double, double f1) { return -1e-2 * f1; };
    InclusiveSampler::GridSpec g;
    g.nx = 24; g.nq2 = 18;
    g.x_min = 1e-3; g.x_max = 0.5;
    g.q2_min = 1.5; g.q2_max = 100.0;
    return std::make_shared<InclusiveSampler>(
        std::make_shared<InclusiveKernel>(LI6(), opt), default_configs("6Li")[1],
        Scenario(), g);
  }();
  return s;
}

}  // namespace

TEST_CASE("rng: the counter alone fixes the stream") {
  Rng a(12345, 7, 3, 99);
  Rng b(12345, 7, 3, 99);
  for (int i = 0; i < 64; ++i) CHECK(a.next_u64() == b.next_u64());

  // every counter slot is a different stream
  const std::uint64_t keys[5][4] = {{12345, 7, 3, 99}, {12346, 7, 3, 99},
                                    {12345, 8, 3, 99}, {12345, 7, 4, 99},
                                    {12345, 7, 3, 100}};
  std::vector<std::uint64_t> first;
  for (const auto& k : keys) {
    Rng r(k[0], k[1], k[2], k[3]);
    first.push_back(r.next_u64());
  }
  std::sort(first.begin(), first.end());
  CHECK(std::adjacent_find(first.begin(), first.end()) == first.end());
}

TEST_CASE("rng: uniform() is uniform on (0, 1) and normal() is standard") {
  Rng r(2026, 0, 0, 0);
  const int n = 1000000;
  const int nbins = 20;
  std::vector<double> counts(nbins, 0.0);
  double sum = 0.0, sum2 = 0.0, lo = 2.0, hi = -1.0;
  for (int i = 0; i < n; ++i) {
    const double u = r.uniform();
    CHECK(u > 0.0);
    CHECK(u < 1.0);
    lo = std::min(lo, u);
    hi = std::max(hi, u);
    sum += u;
    sum2 += u * u;
    counts[std::min(nbins - 1, static_cast<int>(u * nbins))] += 1.0;
  }
  const double mean = sum / n;
  const double var = sum2 / n - mean * mean;
  CHECK_CLOSE_AT(mean, 0.5, 0.0, 5.0 / std::sqrt(12.0 * n));
  CHECK_CLOSE_AT(var, 1.0 / 12.0, 0.0, 1e-3);
  CHECK(lo < 1e-4);
  CHECK(hi > 1.0 - 1e-4);
  // Pearson chi2 on 20 equal bins: mean nbins-1 = 19, sd sqrt(2*19) = 6.2
  double chi2 = 0.0;
  const double expected = static_cast<double>(n) / nbins;
  for (double c : counts) chi2 += (c - expected) * (c - expected) / expected;
  CHECK(chi2 < 19.0 + 5.0 * std::sqrt(38.0));

  double gsum = 0.0, gsum2 = 0.0, gsum4 = 0.0;
  const int ng = 400000;
  for (int i = 0; i < ng; ++i) {
    const double g = r.normal();
    gsum += g;
    gsum2 += g * g;
    gsum4 += g * g * g * g;
  }
  CHECK_CLOSE_AT(gsum / ng, 0.0, 0.0, 5.0 / std::sqrt(static_cast<double>(ng)));
  CHECK_CLOSE_AT(gsum2 / ng, 1.0, 0.0, 0.01);
  CHECK_CLOSE_AT(gsum4 / ng, 3.0, 0.0, 0.1);   // kurtosis of a Gaussian
}

TEST_CASE("rng: rng_poisson has the right mean and variance in both branches") {
  for (double mu : {0.7, 12.0, 250.0}) {
    Rng r(31415, 0, 0, static_cast<std::uint64_t>(mu));
    const int n = 200000;
    double sum = 0.0, sum2 = 0.0;
    for (int i = 0; i < n; ++i) {
      const double k = static_cast<double>(rng_poisson(r, mu));
      sum += k;
      sum2 += k * k;
    }
    const double mean = sum / n;
    const double var = sum2 / n - mean * mean;
    CAPTURE(mu);
    CHECK_CLOSE_AT(mean, mu, 0.0, 6.0 * std::sqrt(mu / n));
    CHECK_CLOSE_AT(var, mu, 0.0, 0.05 * mu);
  }
  Rng r(1, 0, 0, 0);
  CHECK(rng_poisson(r, 0.0) == 0u);
}

TEST_CASE("rng: the sampled batch is identical at 1 and 8 threads") {
  const InclusiveSampler& s = *rng_test_sampler();
  const SpinCategory cat("cos2phi", 1.0, {0.4, 0.2, 0.4}, 0, 0.0, kPi / 2.0,
                         0.3);
  const std::size_t n = 120000;
  const EventBatch one = s.sample_n(cat, n, 424242, 3, 5, 0, 1);
  const EventBatch eight = s.sample_n(cat, n, 424242, 3, 5, 0, 8);
  REQUIRE(one.size() == eight.size());
  CHECK(one.x == eight.x);
  CHECK(one.q2 == eight.q2);
  CHECK(one.y == eight.y);
  CHECK(one.phi == eight.phi);
  CHECK(one.m == eight.m);
  CHECK(one.cell == eight.cell);
  // ... and a Poisson-count run is identical too (the count itself comes off
  // its own reserved counter, so it cannot move with the thread count)
  const EventBatch pa = s.sample_lumi(cat, 5.0, 7, 1, 2, true, 1);
  const EventBatch pb = s.sample_lumi(cat, 5.0, 7, 1, 2, true, 8);
  REQUIRE(pa.size() == pb.size());
  CHECK(pa.phi == pb.phi);
}

TEST_CASE("rng: whole events are identical at 1 and 8 threads") {
  // The generator's own per-event stream, exercised through the full T0
  // assembly (kinematics + nucleon species + four-vectors).
  const auto sampler = rng_test_sampler();
  const InclusiveGenerator gen(sampler);
  const RunPlan plan = tensor_flip_plan(0.6);
  const SpinCategory& cat = plan.categories()[0];
  const InclusiveSampler::CategoryPlan cplan = sampler->make_plan(cat);
  const std::size_t n = 20000;
  const std::uint64_t seed = 5150, run = 2, bunch = 1;

  auto build = [&](std::size_t lo, std::size_t hi, std::vector<Event>& out) {
    for (std::size_t i = lo; i < hi; ++i) {
      Rng rng(seed, run, bunch, i);
      const EventDraw d = sampler->draw_event(cplan, rng);
      out[i] = gen.make_event(cat, plan, d, rng, i, static_cast<int>(run),
                              static_cast<int>(bunch));
    }
  };

  std::vector<Event> serial(n), parallel(n);
  build(0, n, serial);
  {
    std::vector<std::thread> pool;
    const std::size_t chunk = (n + 7) / 8;
    for (unsigned t = 0; t < 8; ++t) {
      const std::size_t lo = std::min(n, t * chunk);
      const std::size_t hi = std::min(n, lo + chunk);
      if (lo >= hi) break;
      pool.emplace_back(build, lo, hi, std::ref(parallel));
    }
    for (std::thread& th : pool) th.join();
  }
  for (std::size_t i = 0; i < n; ++i) {
    const Event& a = serial[i];
    const Event& b = parallel[i];
    REQUIRE(a.particles.size() == b.particles.size());
    CHECK(a.kin.x == b.kin.x);
    CHECK(a.kin.q2 == b.kin.q2);
    CHECK(a.kin.phi == b.kin.phi);
    CHECK(a.spin.m_ion == b.spin.m_ion);
    for (std::size_t k = 0; k < a.particles.size(); ++k) {
      CHECK(a.particles[k].pdg == b.particles[k].pdg);
      CHECK(a.particles[k].p.e == b.particles[k].p.e);
      CHECK(a.particles[k].p.px == b.particles[k].p.px);
      CHECK(a.particles[k].p.py == b.particles[k].p.py);
      CHECK(a.particles[k].p.pz == b.particles[k].p.pz);
    }
  }
}
