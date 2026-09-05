#include "lipolgen/sampler.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <thread>
#include <utility>

#include "lipolgen/numerics.hpp"
#include "lipolgen/sf.hpp"

namespace lipolgen {

// ------------------------------------------------------------- kinematics

ScatteredElectron scattered_electron(double x, double y, double s,
                                     double electron_energy) {
  // From k = (Ee, 0, 0, -Ee), P = (EN, 0, 0, EN):
  //   Q2 = 2 Ee E' (1 + cos theta),  y = 1 - (E'/2Ee)(1 - cos theta)
  // => E'(1 + cos) = Q2/(2 Ee),      E'(1 - cos) = 2 Ee (1 - y)
  const double q2 = s * x * y;
  const double a = q2 / (2.0 * electron_energy);
  const double b = 2.0 * electron_energy * (1.0 - y);
  ScatteredElectron out;
  out.e_prime = 0.5 * (a + b);
  double cos_theta = out.e_prime > 0.0 ? (a - b) / (2.0 * out.e_prime) : -1.0;
  cos_theta = clip(cos_theta, -1.0, 1.0);
  out.theta = std::acos(cos_theta);
  const double tan_half = std::tan(std::min(out.theta, kPi - 1e-9) / 2.0);
  out.eta = -std::log(std::max(tan_half, 1e-12));
  return out;
}

LogGrid log_grid(double x_min, double x_max, double q2_min, double q2_max,
                 int nx, int nq2) {
  if (nx < 1 || nq2 < 1) throw std::runtime_error("log_grid: nx, nq2 >= 1");
  LogGrid g;
  g.x_edges = logspace(std::log10(x_min), std::log10(x_max),
                       static_cast<std::size_t>(nx) + 1);
  g.q2_edges = logspace(std::log10(q2_min), std::log10(q2_max),
                        static_cast<std::size_t>(nq2) + 1);
  g.x_c.resize(static_cast<std::size_t>(nx));
  g.q2_c.resize(static_cast<std::size_t>(nq2));
  for (int i = 0; i < nx; ++i) {
    g.x_c[static_cast<std::size_t>(i)] =
        std::sqrt(g.x_edges[static_cast<std::size_t>(i)] *
                  g.x_edges[static_cast<std::size_t>(i) + 1]);
  }
  for (int j = 0; j < nq2; ++j) {
    g.q2_c[static_cast<std::size_t>(j)] =
        std::sqrt(g.q2_edges[static_cast<std::size_t>(j)] *
                  g.q2_edges[static_cast<std::size_t>(j) + 1]);
  }
  return g;
}

void Scenario::validate() const {
  if (!(run_share > 0.0)) {
    throw std::runtime_error(
        "run_share must be positive (it multiplies the luminosity)");
  }
}

Scenario generator_scenario(const Scenario& analysis, double q2_min,
                            double y_min, double y_max, double w2_min,
                            double eta_pad, double e_prime_min) {
  Scenario s = analysis;
  s.q2_min = q2_min;
  s.y_min = y_min;
  s.y_max = y_max;
  s.w2_min = w2_min;
  s.eta_min = analysis.eta_min - eta_pad;
  s.eta_max = analysis.eta_max + eta_pad;
  s.e_prime_min = e_prime_min;
  return s;
}

// --------------------------------------------------------------- utilities

std::uint64_t rng_poisson(Rng& rng, double mean) {
  if (!(mean > 0.0)) return 0;
  if (mean < 30.0) {
    // Knuth: multiply uniforms until the product drops below exp(-mean).
    const double limit = std::exp(-mean);
    double p = 1.0;
    std::uint64_t k = 0;
    while (true) {
      p *= rng.uniform();
      if (p <= limit) return k;
      ++k;
      if (k > 1000000) return k;  // unreachable for mean < 30
    }
  }
  // Hoermann's PTRS transformed rejection (the algorithm NumPy uses).
  const double smu = std::sqrt(mean);
  const double b = 0.931 + 2.53 * smu;
  const double a = -0.059 + 0.02483 * b;
  const double inv_alpha = 1.1239 + 1.1328 / (b - 3.4);
  const double v_r = 0.9277 - 3.6224 / (b - 2.0);
  const double log_mean = std::log(mean);
  while (true) {
    const double u = rng.uniform() - 0.5;
    const double v = rng.uniform();
    const double us = 0.5 - std::fabs(u);
    const double kf = std::floor((2.0 * a / us + b) * u + mean + 0.43);
    if (us >= 0.07 && v <= v_r) return static_cast<std::uint64_t>(kf);
    if (kf < 0.0 || (us < 0.013 && v > us)) continue;
    if (std::log(v * inv_alpha / (a / (us * us) + b)) <=
        -mean + kf * log_mean - std::lgamma(kf + 1.0)) {
      return static_cast<std::uint64_t>(kf);
    }
  }
}

void EventBatch::reserve(std::size_t n) {
  x.reserve(n); q2.reserve(n); y.reserve(n); phi.reserve(n);
  m.reserve(n); weight.reserve(n); cell.reserve(n);
}

void EventBatch::resize(std::size_t n) {
  x.resize(n); q2.resize(n); y.resize(n); phi.resize(n);
  m.resize(n); weight.resize(n); cell.resize(n);
}

// ----------------------------------------------------------------- sampler

bool InclusiveSampler::StateKey::operator<(const StateKey& o) const {
  if (lam_e != o.lam_e) return lam_e < o.lam_e;
  if (pe != o.pe) return pe < o.pe;
  if (j != o.j) return j < o.j;
  if (m != o.m) return m < o.m;
  if (theta_s != o.theta_s) return theta_s < o.theta_s;
  return phi_s < o.phi_s;
}

InclusiveSampler::InclusiveSampler(std::shared_ptr<const InclusiveKernel> kernel,
                                   BeamConfig config, Scenario scenario,
                                   GridSpec grid, bool with_perp)
    : kernel_(std::move(kernel)), config_(std::move(config)),
      scenario_(scenario), spec_(grid), with_perp_(with_perp) {
  if (!kernel_) throw std::runtime_error("InclusiveSampler: null kernel");
  scenario_.validate();
  s_ = config_.s_per_nucleon();
  grid_ = log_grid(spec_.x_min, spec_.x_max, spec_.q2_min, spec_.q2_max,
                   spec_.nx, spec_.nq2);

  const NuclearF2& nf2 = kernel_->nuclear_f2();
  const double a = static_cast<double>(kernel_->ion().A);
  const RFunc& r_func = nf2.r_func();
  // The ONE luminosity line: programme luminosity x this observable's share.
  // Dividing `n_events` by the EFFECTIVE luminosity is what makes the per-cell
  // pb number share-invariant (test_run_share.py).
  const double lumi_pb = scenario_.lumi_effective_pb_per_nucleon();

  for (int i = 0; i < spec_.nx; ++i) {
    const std::size_t ii = static_cast<std::size_t>(i);
    const double xc = grid_.x_c[ii];
    const double dx = grid_.x_edges[ii + 1] - grid_.x_edges[ii];
    for (int j = 0; j < spec_.nq2; ++j) {
      const std::size_t jj = static_cast<std::size_t>(j);
      const double q2c = grid_.q2_c[jj];
      if (!in_acceptance(xc, q2c)) continue;
      const double dq2 = grid_.q2_edges[jj + 1] - grid_.q2_edges[jj];
      const double f2 = nf2.f2a(xc, q2c) / a;
      const double xsec = dsigma_dx_dq2(xc, q2c, s_, f2, nullptr, r_func);
      const double n_events = xsec * dx * dq2 * lumi_pb;
      const double cell_pb = n_events / lumi_pb;
      if (!(cell_pb > 0.0)) continue;
      x_cells_.push_back(xc);
      q2_cells_.push_back(q2c);
      xsec_flat_.push_back(cell_pb);
      logx_lo_.push_back(std::log(grid_.x_edges[ii]));
      logx_hi_.push_back(std::log(grid_.x_edges[ii + 1]));
      logq2_lo_.push_back(std::log(grid_.q2_edges[jj]));
      logq2_hi_.push_back(std::log(grid_.q2_edges[jj + 1]));
      tables_.push_back(kernel_->tables(xc, q2c, with_perp_));
    }
  }
  if (x_cells_.empty()) {
    throw std::runtime_error(
        "InclusiveSampler: no accepted (x, Q2) cell -- the grid and the "
        "acceptance window do not overlap");
  }
}

bool InclusiveSampler::in_acceptance(double x, double q2) const {
  const Scenario& sc = scenario_;
  if (!(q2 >= sc.q2_min)) return false;
  if (!(x <= sc.x_max)) return false;
  const double y = y_from_xq2(x, q2, s_);
  if (!(y >= sc.y_min) || !(y <= sc.y_max)) return false;
  if (!(w2_from_xq2(x, q2) >= sc.w2_min)) return false;
  const double yc = clip(y, 1e-9, 1.0 - 1e-12);
  const ScatteredElectron e = scattered_electron(x, yc, s_,
                                                 config_.electron_energy);
  return e.eta >= sc.eta_min && e.eta <= sc.eta_max &&
         e.e_prime >= sc.e_prime_min;
}

const InclusiveSampler::StateTables& InclusiveSampler::build_state(
    const StateKey& key) const {
  const std::size_t n = x_cells_.size();
  StateTables st;
  st.w_avg.resize(n); st.a1.resize(n); st.a2.resize(n);
  st.w_tensor.resize(n); st.a1_tensor.resize(n); st.a2_tensor.resize(n);
  st.a1n.resize(n); st.a2n.resize(n); st.bound.resize(n); st.cdf.resize(n);

  EventSpinState state;
  state.lam_e = key.lam_e;
  state.pe = key.pe;
  state.j = key.j;
  state.m = key.m;
  state.theta_s = key.theta_s;
  state.phi_s = key.phi_s;

  double worst = 1e300;
  std::size_t worst_cell = 0;
  std::size_t bad = 0;
  double acc = 0.0;
  for (std::size_t i = 0; i < n; ++i) {
    const Amplitudes amp = kernel_->amplitudes(tables_[i], x_cells_[i],
                                               q2_cells_[i], s_, state,
                                               with_perp_);
    const double den = 1.0 + amp.w_avg;
    if (!(den > 0.0)) {
      // NAME THE CELL AND THE CURE.  This is a CONFIGURATION-TIME refusal --
      // it fires while the `Pipeline` is being built, before a single event
      // is drawn -- and the thing the caller has to change is the acceptance
      // WINDOW, `Scenario::x_max` (the CLI's `--x-max`), which is what drops
      // the offending cell.  Until 2026-09-05 the message stopped after the
      // parenthesis and a user met it as a bare traceback naming no way out;
      // the combination that first showed it is documented under "the top x
      // cell" in docs/USAGE.md.  The A1 = g1/F1 of the cell is printed
      // because it is the quantity that is too large: the phi-averaged
      // modulation is P_e P_z D(y) A1 plus the tensor term, so a backend
      // pair whose A1 runs away at high x is refused here and a milder one
      // is not.  NO CLAMP: max(g1, ...) here would be a silent physics
      // change, and the sampler's own accept-reject would then draw
      // max(W, 0), diluting the modulation AND skewing the (x, Q2) mixture.
      const double f1 = tables_[i].f1;
      char msg[1024];
      std::snprintf(msg, sizeof(msg),
                    "negative phi-averaged density for m=%g at x = %.4g, "
                    "Q2 = %.4g (1 + w_avg = %.4g).  That cell's A1 = g1/F1 "
                    "is %.4g (F1 = %.4g, g1 = %.4g) at lam_e = %d, "
                    "P_e = %.4g, J = %g: the phi-averaged rate would go "
                    "NEGATIVE there, and the sampler refuses rather than "
                    "draw max(W, 0), which would dilute the modulation and "
                    "skew the (x, Q2) mixture.  Cure: lower the acceptance "
                    "window's x_max below %.4g (Scenario::x_max; the CLI's "
                    "--x-max, e.g. --x-max 0.95 on the shipped grid, whose "
                    "top cell is x = 0.955), or reduce P_e / P_z, which "
                    "scale w_avg linearly.  It is a property of the (x, Q2) "
                    "CELL and the structure-function backends, not of the "
                    "event count.",
                    key.m, x_cells_[i], q2_cells_[i], den,
                    f1 != 0.0 ? tables_[i].g1 / f1 : 0.0, f1, tables_[i].g1,
                    key.lam_e, key.pe, key.j, x_cells_[i]);
      throw std::runtime_error(msg);
    }
    st.w_avg[i] = amp.w_avg;
    st.a1[i] = amp.a1;
    st.a2[i] = amp.a2;
    // The rank-2 (b1..b4) part of the SAME amplitudes, for the RC band of
    // rc.hpp.  Computed here, once per (spin state, cell), so that no consumer
    // ever has to re-enter the kernel in the event loop -- `tables()` alone is
    // a 96-point g2^WW quadrature per call at the default `target_mass = true`.
    // Nothing in the sampler reads these; they change no drawn number.
    const Amplitudes tam = kernel_->tensor_amplitudes(tables_[i], x_cells_[i],
                                                      q2_cells_[i], s_, state);
    st.w_tensor[i] = tam.w_avg;
    st.a1_tensor[i] = tam.a1;
    st.a2_tensor[i] = tam.a2;
    st.a1n[i] = amp.a1 / den;
    st.a2n[i] = amp.a2 / den;
    st.bound[i] = 1.0 + std::fabs(st.a1n[i]) + std::fabs(st.a2n[i]);
    const double margin = density_min(st.a1n[i], st.a2n[i]);
    if (margin < worst) { worst = margin; worst_cell = i; }
    if (margin < 0.0) ++bad;
    acc += xsec_flat_[i] * den;
    st.cdf[i] = acc;
  }
  st.sigma_pb = acc;
  st.margin = worst;
  if (bad > 0) {
    char msg[512];
    std::snprintf(msg, sizeof(msg),
                  "negative phi density for m=%g: %zu of %zu cells, worst "
                  "min(W)/(1+w_avg) = %.4f at x = %.4g, Q2 = %.4g "
                  "(a1/den = %.3f, a2/den = %.3f).  The sampler would "
                  "silently draw max(W, 0), diluting the modulation and "
                  "skewing the (x, Q2) mixture; reduce the scenario "
                  "amplitude.",
                  key.m, bad, n, worst, x_cells_[worst_cell],
                  q2_cells_[worst_cell], st.a1n[worst_cell],
                  st.a2n[worst_cell]);
    throw std::runtime_error(msg);
  }

  std::lock_guard<std::mutex> lock(cache_mutex_);
  auto it = cache_.find(key);
  if (it == cache_.end()) {
    it = cache_.emplace(key, std::move(st)).first;
  }
  return it->second;  // std::map nodes are stable, so the reference outlives
                      // any later insertion.
}

const InclusiveSampler::StateTables& InclusiveSampler::state_tables(
    const SpinCategory& cat, double m) const {
  const StateKey key{cat.lam_e, cat.pe, cat.j, m, cat.theta_s, cat.phi_s};
  {
    std::lock_guard<std::mutex> lock(cache_mutex_);
    const auto it = cache_.find(key);
    if (it != cache_.end()) return it->second;
  }
  return build_state(key);
}

double InclusiveSampler::sigma_tot_pb(const SpinCategory& cat) const {
  const std::vector<double> ms = m_values(cat.j);
  if (cat.populations.size() != ms.size()) {
    throw std::runtime_error("populations must have 2j+1 entries (m=+J..-J)");
  }
  double tot = 0.0;
  for (std::size_t i = 0; i < ms.size(); ++i) {
    if (cat.populations[i] <= 0.0) continue;
    tot += cat.populations[i] * state_tables(cat, ms[i]).sigma_pb;
  }
  return tot;
}

InclusiveSampler::EffectiveModulation InclusiveSampler::effective_modulation(
    const SpinCategory& cat, const std::vector<char>* mask) const {
  const std::vector<double> ms = m_values(cat.j);
  if (cat.populations.size() != ms.size()) {
    throw std::runtime_error("populations must have 2j+1 entries (m=+J..-J)");
  }
  if (mask && mask->size() != x_cells_.size()) {
    throw std::runtime_error("effective_modulation: mask size mismatch");
  }
  double den = 0.0, num1 = 0.0, num2 = 0.0;
  for (std::size_t k = 0; k < ms.size(); ++k) {
    const double p_m = cat.populations[k];
    if (p_m <= 0.0) continue;
    const StateTables& st = state_tables(cat, ms[k]);
    for (std::size_t i = 0; i < x_cells_.size(); ++i) {
      if (mask && !(*mask)[i]) continue;
      den += p_m * xsec_flat_[i] * (1.0 + st.w_avg[i]);
      num1 += p_m * xsec_flat_[i] * st.a1[i];
      num2 += p_m * xsec_flat_[i] * st.a2[i];
    }
  }
  EffectiveModulation out;
  if (!(den > 0.0)) return out;
  out.sigma_pb = den;
  out.a1 = num1 / den;
  out.a2 = num2 / den;
  return out;
}

InclusiveSampler::CategoryPlan InclusiveSampler::make_plan(
    const SpinCategory& cat) const {
  const std::vector<double> ms = m_values(cat.j);
  if (cat.populations.size() != ms.size()) {
    throw std::runtime_error("populations must have 2j+1 entries (m=+J..-J)");
  }
  CategoryPlan plan;
  plan.phi_s = cat.phi_s;
  double acc = 0.0;
  for (std::size_t i = 0; i < ms.size(); ++i) {
    const double p_m = cat.populations[i];
    if (p_m <= 0.0) continue;
    const StateTables& st = state_tables(cat, ms[i]);
    acc += p_m * st.sigma_pb;
    plan.m_val.push_back(ms[i]);
    plan.m_cdf.push_back(acc);
    plan.states.push_back(&st);
  }
  if (plan.m_val.empty()) {
    throw std::runtime_error("make_plan: category has no populated m state");
  }
  plan.sigma_pb = acc;
  return plan;
}

EventDraw InclusiveSampler::draw_event(const CategoryPlan& plan,
                                       Rng& rng) const {
  EventDraw ev;
  // 1. the spin projection, from the population-weighted RATES
  const double um = rng.uniform() * plan.sigma_pb;
  std::size_t k = static_cast<std::size_t>(
      std::lower_bound(plan.m_cdf.begin(), plan.m_cdf.end(), um) -
      plan.m_cdf.begin());
  if (k >= plan.m_val.size()) k = plan.m_val.size() - 1;
  ev.m = plan.m_val[k];
  const StateTables& st = *plan.states[k];

  // 2. the (x, Q2) cell by inverse CDF, then log-uniform inside it
  const double uc = rng.uniform() * st.sigma_pb;
  std::size_t c = static_cast<std::size_t>(
      std::lower_bound(st.cdf.begin(), st.cdf.end(), uc) - st.cdf.begin());
  if (c >= st.cdf.size()) c = st.cdf.size() - 1;
  ev.cell = static_cast<int>(c);

  double x = 0.0, q2 = 0.0;
  bool ok = false;
  for (int tries = 0; tries < 21; ++tries) {
    const double u1 = rng.uniform();
    const double u2 = rng.uniform();
    x = std::exp(logx_lo_[c] + u1 * (logx_hi_[c] - logx_lo_[c]));
    q2 = std::exp(logq2_lo_[c] + u2 * (logq2_hi_[c] - logq2_lo_[c]));
    if (in_acceptance(x, q2)) { ok = true; break; }
  }
  if (!ok) {  // guaranteed termination: fall back to the accepted cell center
    x = x_cells_[c];
    q2 = q2_cells_[c];
  }
  ev.x = x;
  ev.q2 = q2;
  ev.y = y_from_xq2(x, q2, s_);

  // 3. phi by accept-reject on the cell-center amplitudes
  const double a1n = st.a1n[c], a2n = st.a2n[c], bound = st.bound[c];
  double phip = 0.0;
  while (true) {
    const double cand = rng.uniform() * 2.0 * kPi;
    const double u = rng.uniform() * bound;
    if (u < 1.0 + a1n * std::cos(cand) + a2n * std::cos(2.0 * cand)) {
      phip = cand;
      break;
    }
  }
  double phi = std::fmod(phip + plan.phi_s, 2.0 * kPi);
  if (phi < 0.0) phi += 2.0 * kPi;
  ev.phi = phi;
  return ev;
}

EventBatch InclusiveSampler::sample_n(const SpinCategory& cat, std::size_t n,
                                      std::uint64_t seed, std::uint64_t run,
                                      std::uint64_t bunch,
                                      std::uint64_t event0,
                                      unsigned nthreads) const {
  const CategoryPlan plan = make_plan(cat);  // fills the cache single-threaded
  EventBatch out;
  out.category = cat.name;
  out.lam_e = cat.lam_e;
  out.resize(n);
  if (n == 0) return out;

  auto worker = [&](std::size_t lo, std::size_t hi) {
    for (std::size_t i = lo; i < hi; ++i) {
      Rng rng(seed, run, bunch, event0 + i);
      const EventDraw d = draw_event(plan, rng);
      out.x[i] = d.x;
      out.q2[i] = d.q2;
      out.y[i] = d.y;
      out.phi[i] = d.phi;
      out.m[i] = d.m;
      out.cell[i] = d.cell;
      out.weight[i] = 1.0;
    }
  };

  if (nthreads <= 1 || n < 2 * nthreads) {
    worker(0, n);
    return out;
  }
  std::vector<std::thread> pool;
  pool.reserve(nthreads);
  const std::size_t chunk = (n + nthreads - 1) / nthreads;
  for (unsigned t = 0; t < nthreads; ++t) {
    const std::size_t lo = std::min(n, static_cast<std::size_t>(t) * chunk);
    const std::size_t hi = std::min(n, lo + chunk);
    if (lo >= hi) break;
    pool.emplace_back(worker, lo, hi);
  }
  for (std::thread& th : pool) th.join();
  return out;
}

EventBatch InclusiveSampler::sample_lumi(const SpinCategory& cat,
                                         double lumi_pb, std::uint64_t seed,
                                         std::uint64_t run,
                                         std::uint64_t bunch, bool poisson,
                                         unsigned nthreads) const {
  const double mu = lumi_pb * sigma_tot_pb(cat);
  std::size_t n = 0;
  if (poisson) {
    Rng counter(seed, run, bunch, kCountStreamEvent);
    n = static_cast<std::size_t>(rng_poisson(counter, mu));
  } else {
    n = static_cast<std::size_t>(std::llround(mu));
  }
  return sample_n(cat, n, seed, run, bunch, 0, nthreads);
}

std::vector<double> InclusiveSampler::weights_for(
    const EventBatch& batch, const std::vector<SpinCategory>& cats) const {
  const std::size_t n = batch.size();
  const std::size_t nk = cats.size();
  std::vector<double> out(n * nk, 0.0);
  for (std::size_t k = 0; k < nk; ++k) {
    const SpinCategory& cat = cats[k];
    const std::vector<double> ms = m_values(cat.j);
    if (cat.populations.size() != ms.size()) {
      throw std::runtime_error("populations must have 2j+1 entries");
    }
    for (std::size_t im = 0; im < ms.size(); ++im) {
      const double p_m = cat.populations[im];
      if (p_m <= 0.0) continue;
      const StateTables& st = state_tables(cat, ms[im]);
      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t c = static_cast<std::size_t>(batch.cell[i]);
        const double phip = batch.phi[i] - cat.phi_s;
        out[i * nk + k] += p_m * (1.0 + st.w_avg[c] +
                                  st.a1[c] * std::cos(phip) +
                                  st.a2[c] * std::cos(2.0 * phip));
      }
    }
  }
  return out;
}

std::vector<double> InclusiveSampler::tensor_weights_for(
    const EventBatch& batch, const std::vector<SpinCategory>& cats) const {
  // The exact mirror of `weights_for` above, on the tensor vectors of the same
  // StateTables and WITHOUT the unpolarized 1.0: t[i, k] is the rank-2 part of
  // the density w[i, k], so their ratio is the tau the RC band scales.
  const std::size_t n = batch.size();
  const std::size_t nk = cats.size();
  std::vector<double> out(n * nk, 0.0);
  for (std::size_t k = 0; k < nk; ++k) {
    const SpinCategory& cat = cats[k];
    const std::vector<double> ms = m_values(cat.j);
    if (cat.populations.size() != ms.size()) {
      throw std::runtime_error("populations must have 2j+1 entries");
    }
    for (std::size_t im = 0; im < ms.size(); ++im) {
      const double p_m = cat.populations[im];
      if (p_m <= 0.0) continue;
      const StateTables& st = state_tables(cat, ms[im]);
      for (std::size_t i = 0; i < n; ++i) {
        const std::size_t c = static_cast<std::size_t>(batch.cell[i]);
        const double phip = batch.phi[i] - cat.phi_s;
        out[i * nk + k] += p_m * (st.w_tensor[c] +
                                  st.a1_tensor[c] * std::cos(phip) +
                                  st.a2_tensor[c] * std::cos(2.0 * phip));
      }
    }
  }
  return out;
}

// ------------------------------------------------------- pseudo-experiments

const EventBatch& PseudoExperiment::operator[](const std::string& name) const {
  for (std::size_t i = 0; i < names.size(); ++i) {
    if (names[i] == name) return batches[i];
  }
  throw std::runtime_error("PseudoExperiment: no category " + name);
}

double PseudoExperiment::lumi_of(const std::string& name) const {
  for (std::size_t i = 0; i < names.size(); ++i) {
    if (names[i] == name) return lumi_pb[i];
  }
  throw std::runtime_error("PseudoExperiment: no category " + name);
}

std::size_t PseudoExperiment::count(const std::string& name) const {
  return (*this)[name].size();
}

PseudoExperiment run_pseudo_experiment(const InclusiveSampler& sampler,
                                       const RunPlan& plan,
                                       double total_lumi_pb,
                                       std::uint64_t seed, std::uint64_t run,
                                       bool poisson, unsigned nthreads) {
  PseudoExperiment out;
  const std::vector<double> shares = plan.lumi_share_vector(total_lumi_pb);
  for (std::size_t i = 0; i < plan.categories().size(); ++i) {
    const SpinCategory& cat = plan.categories()[i];
    out.names.push_back(cat.name);
    out.lumi_pb.push_back(shares[i]);
    out.batches.push_back(sampler.sample_lumi(cat, shares[i], seed, run, i,
                                              poisson, nthreads));
  }
  return out;
}

PhiHistogram phi_histogram_pseudo(double n_expected, double a2, int nbins,
                                  Rng* rng, double a1, bool poisson) {
  if (nbins < 1) throw std::runtime_error("phi_histogram_pseudo: nbins >= 1");
  const double margin = density_min(a1, a2);
  if (margin < 0.0) {
    char msg[256];
    std::snprintf(msg, sizeof(msg),
                  "negative phi' density: min(1 + %.4g cos phi' + %.4g "
                  "cos 2phi') = %.4f", a1, a2, margin);
    throw std::runtime_error(msg);
  }
  PhiHistogram out;
  out.edges = linspace(0.0, 2.0 * kPi, static_cast<std::size_t>(nbins) + 1);
  out.counts.resize(static_cast<std::size_t>(nbins));
  Rng fallback(20260713, 0, 0, 0);
  Rng& r = rng ? *rng : fallback;
  for (int i = 0; i < nbins; ++i) {
    const double lo = out.edges[static_cast<std::size_t>(i)];
    const double hi = out.edges[static_cast<std::size_t>(i) + 1];
    const double mu = (n_expected / (2.0 * kPi)) *
                      ((hi - lo) + a1 * (std::sin(hi) - std::sin(lo)) +
                       0.5 * a2 * (std::sin(2.0 * hi) - std::sin(2.0 * lo)));
    out.counts[static_cast<std::size_t>(i)] =
        poisson ? static_cast<double>(rng_poisson(r, mu)) : mu;
  }
  return out;
}

// --------------------------------------------------------------- estimators

namespace estimators {

std::vector<double> yields(const std::vector<double>& counts,
                           const std::vector<double>& lumis) {
  if (lumis.empty()) return counts;
  if (lumis.size() != counts.size()) {
    throw std::runtime_error("yields: counts and lumis must align");
  }
  double mean = 0.0;
  for (double l : lumis) mean += l;
  mean /= static_cast<double>(lumis.size());
  std::vector<double> out(counts.size());
  for (std::size_t i = 0; i < counts.size(); ++i) {
    out[i] = counts[i] / lumis[i] * mean;
  }
  return out;
}

double apar_flip(double n_plus, double n_minus, double pe, double pz,
                 double l_plus, double l_minus) {
  std::vector<double> lumis;
  if (l_plus > 0.0 || l_minus > 0.0) lumis = {l_plus, l_minus};
  const std::vector<double> y = yields({n_plus, n_minus}, lumis);
  return (y[0] - y[1]) / (y[0] + y[1]) / (pe * pz);
}

double azz_thirds(double n_plus, double n_minus, double n_zero, double pzz,
                  const std::vector<double>& lumis) {
  const std::vector<double> y = yields({n_plus, n_minus, n_zero}, lumis);
  return (y[0] + y[1] - 2.0 * y[2]) / (y[0] + y[1] + y[2]) / pzz;
}

double cos2phi_moment(const std::vector<double>& phi_prime, double pzz) {
  if (phi_prime.empty()) throw std::runtime_error("cos2phi_moment: no events");
  double s = 0.0;
  for (double p : phi_prime) s += std::cos(2.0 * p);
  return 2.0 * (s / static_cast<double>(phi_prime.size())) / pzz;
}

double cos2phi_fit_binned(const std::vector<double>& counts,
                          const std::vector<double>& edges, double pzz,
                          const std::function<bool(double)>& acceptance,
                          int nsub) {
  const std::size_t nbins = counts.size();
  if (edges.size() != nbins + 1 || nbins == 0) {
    throw std::runtime_error("cos2phi_fit_binned: edges must be counts+1");
  }
  const double width = edges[1] - edges[0];
  std::vector<char> use(nbins, 1);
  if (acceptance) {
    const std::vector<double> frac = linspace(0.0, 1.0,
                                              static_cast<std::size_t>(nsub));
    for (std::size_t i = 0; i < nbins; ++i) {
      for (double f : frac) {
        if (!acceptance(edges[i] + f * width)) { use[i] = 0; break; }
      }
    }
  }
  std::size_t nlive = 0;
  for (char u : use) nlive += (u != 0);
  if (nlive < 3) {
    char msg[160];
    std::snprintf(msg, sizeof(msg),
                  "cos2phi_fit_binned: only %zu live phi bins after the "
                  "acceptance probe -- amplitude unidentifiable", nlive);
    throw std::runtime_error(msg);
  }
  // Linear LSQ for N = C + B cos 2phi through the 2x2 normal equations.
  double s11 = 0.0, s12 = 0.0, s22 = 0.0, t1 = 0.0, t2 = 0.0;
  for (std::size_t i = 0; i < nbins; ++i) {
    if (!use[i]) continue;
    const double c2 = std::cos(2.0 * 0.5 * (edges[i] + edges[i + 1]));
    s11 += 1.0;
    s12 += c2;
    s22 += c2 * c2;
    t1 += counts[i];
    t2 += c2 * counts[i];
  }
  const double det = s11 * s22 - s12 * s12;
  if (std::fabs(det) < 1e-300) {
    throw std::runtime_error("cos2phi_fit_binned: singular design matrix");
  }
  const double c_coef = (s22 * t1 - s12 * t2) / det;
  const double b_coef = (s11 * t2 - s12 * t1) / det;
  const double amp = b_coef / c_coef;
  const double w = 0.5 * width;              // half-width; dilution in 2phi
  const double dilution = std::sin(2.0 * w) / (2.0 * w);
  return amp / dilution / pzz;
}

double cos2phi_fit(const std::vector<double>& phi_prime, double pzz, int nbins,
                   const std::function<bool(double)>& acceptance, int nsub) {
  if (nbins < 1) throw std::runtime_error("cos2phi_fit: nbins >= 1");
  const std::vector<double> edges =
      linspace(0.0, 2.0 * kPi, static_cast<std::size_t>(nbins) + 1);
  std::vector<double> counts(static_cast<std::size_t>(nbins), 0.0);
  const double two_pi = 2.0 * kPi;
  for (double p : phi_prime) {
    double v = std::fmod(p, two_pi);
    if (v < 0.0) v += two_pi;
    std::size_t b = static_cast<std::size_t>(v / two_pi *
                                             static_cast<double>(nbins));
    if (b >= counts.size()) b = counts.size() - 1;
    counts[b] += 1.0;
  }
  return cos2phi_fit_binned(counts, edges, pzz, acceptance, nsub);
}

double cos2phi_fit_err(double n, double pzz, int nbins) {
  const double nn = std::max(n, 1e-12);
  const double w = kPi / static_cast<double>(nbins);
  const double dilution = std::sin(2.0 * w) / (2.0 * w);
  return std::sqrt(2.0 / nn) / pzz / dilution;
}

}  // namespace estimators

}  // namespace lipolgen
