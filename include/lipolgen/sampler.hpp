#ifndef LIPOLGEN_SAMPLER_HPP
#define LIPOLGEN_SAMPLER_HPP

/// \file sampler.hpp
/// Inclusive spin-labeled event sampler (Mode G, inclusive tier) plus the
/// acceptance window, the (x, Q2) cell grid and the analysis-side estimators.
/// C++17 port of `evgen/polligen/sample.py`, `evgen/polligen/estimators.py`,
/// the `Scenario` / `project_rates` half of `fastsim/polli_fastsim/fom.py`
/// and `fastsim/polli_fastsim/kinematics.py`.
///
/// SAMPLING SCHEME (plans/05 5.2.4)
///
/// 1. per spin category (`SpinCategory`) and per spin projection m (drawn
///    from the fill populations), the event count is Poisson with mean
///    L_cat * p_m * sum_cells sigma_cell * (1 + w_avg(m)) -- the phi-averaged
///    polarized modulation shifts the RATE per spin state, which is exactly
///    what counting-asymmetry estimators measure;
/// 2. (x, Q2) by inverse-CDF over a log-log cell grid of the unpolarized
///    cross section times (1 + w_avg), log-uniform inside the cell, redrawn
///    until the event lands inside the acceptance window;
/// 3. phi by accept-reject on 1 + a1 cos(phi') + a2 cos(2 phi'),
///    phi' = phi - phi_S, refused outright if the density can go negative.
///
/// AZIMUTH.  `phi` is the azimuth of the SCATTERED ELECTRON about the ion
/// (+z) beam axis in the HEAD-ON frame, phi = atan2(p_y, p_x) of e', in
/// [0, 2 pi).  `phi' = phi - phi_S` is then the covariant Bacchetta et al.
/// phi_S of the alignment axis exactly (massless target; `reco.py`
/// `azimuth_wrt_lepton_plane`).  generator.hpp builds e' from this angle.
///
/// LUMINOSITY.  A run-plan share moves COUNTS, never CROSS SECTIONS: the
/// per-cell pb number divides `n_events` by the EFFECTIVE luminosity
/// (programme x `run_share`), so `cell_xsec_pb()` is share-invariant and the
/// caller's `lumi_pb` applies the share exactly once (test_run_share.py).
///
/// DETERMINISM.  Every event i of a (seed, run, bunch) stream is drawn from
/// its own counter-based `Rng(seed, run, bunch, i)`, so the batch is
/// identical at any thread count and independent of generation order.

#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "lipolgen/beams.hpp"
#include "lipolgen/bookkeeping.hpp"
#include "lipolgen/constants.hpp"
#include "lipolgen/rng.hpp"
#include "lipolgen/xsec.hpp"

namespace lipolgen {

// ------------------------------------------------------------- kinematics

/// y = Q2/(s x) with the per-nucleon s.
inline double y_from_xq2(double x, double q2, double s) { return q2 / (s * x); }

/// Invariant mass^2 of the hadronic system on a per-nucleon target.
inline double w2_from_xq2(double x, double q2, double m_n = M_NUCLEON) {
  return m_n * m_n + q2 * (1.0 - x) / (x > 1e-12 ? x : 1e-12);
}

/// Head-on-frame scattered electron.  theta is measured from +z (the ION
/// direction), so the scattered electron sits at theta ~ pi (negative eta).
struct ScatteredElectron {
  double e_prime = 0.0;
  double theta = 0.0;
  double eta = 0.0;
};
ScatteredElectron scattered_electron(double x, double y, double s,
                                     double electron_energy);

/// Bin edges and (geometric) centers of a log-log (x, Q2) grid.
struct LogGrid {
  std::vector<double> x_edges, q2_edges, x_c, q2_c;
};
LogGrid log_grid(double x_min, double x_max, double q2_min, double q2_max,
                 int nx, int nq2);

// --------------------------------------------------------------- scenario

/// Acceptance window + luminosity bookkeeping.  Port of `fom.Scenario`.
///
/// `lumi_fb_per_nucleon` is the PROGRAMME luminosity ("one EIC year" =
/// 10 fb^-1/u) and `run_share` is the fraction of it this observable is
/// given.  The two are kept apart so that a projection can say "10 fb^-1/u of
/// programme, 1/3 of it here".  Every published number is at run_share = 1.
struct Scenario {
  double lumi_fb_per_nucleon = 10.0;
  double run_share = 1.0;
  double pol_electron = 0.70;
  double pol_ion_vector = 0.70;   ///< P_z in the ring (placeholder)
  double pol_ion_tensor = 0.60;   ///< P_zz in the ring (placeholder)
  double q2_min = 1.0;
  double y_min = 0.01;
  double y_max = 0.95;
  double w2_min = 10.0;
  double x_max = 1.0;
  /// crude central-detector acceptance for the scattered electron
  double eta_min = -3.5;
  double eta_max = 3.5;
  double e_prime_min = 0.5;  ///< GeV

  /// Luminosity this observable actually receives [fb^-1/nucleon].
  double lumi_effective_fb_per_nucleon() const {
    return lumi_fb_per_nucleon * run_share;
  }
  /// ... in pb^-1, the unit every cross section here is quoted in.
  double lumi_effective_pb_per_nucleon() const {
    return lumi_effective_fb_per_nucleon() * 1e3;
  }
  /// Throws if `run_share` is not positive (it multiplies the luminosity).
  void validate() const;
};

/// Looser copy of an analysis Scenario for the GENERATOR: events outside the
/// analysis cuts must exist so that they can migrate in
/// (`recopseudo.generator_scenario`).  Defaults are the conventions of
/// docs/DEVELOPMENT_PLAN.md 5: Q2 >= 0.7, y in [0.004, 0.985], W2 >= 8.
Scenario generator_scenario(const Scenario& analysis, double q2_min = 0.7,
                            double y_min = 0.004, double y_max = 0.985,
                            double w2_min = 8.0, double eta_pad = 0.3,
                            double e_prime_min = 0.3);

// --------------------------------------------------------------- utilities

/// Poisson deviate: Knuth's product method below mean 30, Hoermann's
/// transformed rejection (PTRS, the algorithm NumPy uses) above it.
std::uint64_t rng_poisson(Rng& rng, double mean);

// ----------------------------------------------------------------- sampler

/// Column-oriented batch of sampled events.  One row per event; `cell` is the
/// index into the sampler's ACCEPTED-cell arrays.
struct EventBatch {
  std::vector<double> x, q2, y, phi, m, weight;
  std::vector<int> cell;
  std::string category;
  int lam_e = 0;

  std::size_t size() const { return x.size(); }
  void reserve(std::size_t n);
  void resize(std::size_t n);
};

/// Log-log (x, Q2) cell grid the inverse-CDF is built on.  The default
/// 100x72 is ~4x finer than the 40x30 FOM analysis binning.
struct SamplerGrid {
  int nx = 100;
  int nq2 = 72;
  double x_min = 1e-4, x_max = 1.0;
  double q2_min = 1.0, q2_max = 2e3;
};

/// One sampled event, before it becomes a `lipolgen::Event`.
struct EventDraw {
  double x = 0.0, q2 = 0.0, y = 0.0, phi = 0.0, m = 0.0;
  int cell = -1;
};

/// Grid-backed sampler for one beam configuration + kernel.
///
/// Fidelity note: modulation amplitudes are evaluated at CELL CENTERS (the
/// default 100x72 grid is ~4x finer than the 40x30 FOM analysis binning);
/// in-cell (x, Q2) placement is log-uniform.  Adequate for binned Phase-1
/// estimators by construction; not for unbinned in-cell shapes.
class InclusiveSampler {
 public:
  using GridSpec = SamplerGrid;

  InclusiveSampler(std::shared_ptr<const InclusiveKernel> kernel,
                   BeamConfig config, Scenario scenario = Scenario(),
                   GridSpec grid = GridSpec(), bool with_perp = false);

  // --- the grid -----------------------------------------------------------
  const InclusiveKernel& kernel() const { return *kernel_; }
  const std::shared_ptr<const InclusiveKernel>& kernel_ptr() const {
    return kernel_;
  }
  const BeamConfig& config() const { return config_; }
  const Scenario& scenario() const { return scenario_; }
  double s() const { return s_; }
  bool with_perp() const { return with_perp_; }

  std::size_t n_cells() const { return x_cells_.size(); }
  const std::vector<double>& x_cells() const { return x_cells_; }
  const std::vector<double>& q2_cells() const { return q2_cells_; }
  /// Per-cell UNPOLARIZED accepted cross section [pb].  Share-invariant.
  const std::vector<double>& cell_xsec_pb() const { return xsec_flat_; }
  const std::vector<double>& logx_lo() const { return logx_lo_; }
  const std::vector<double>& logx_hi() const { return logx_hi_; }
  const std::vector<double>& logq2_lo() const { return logq2_lo_; }
  const std::vector<double>& logq2_hi() const { return logq2_hi_; }
  const std::vector<SFTables>& tables() const { return tables_; }
  const LogGrid& grid() const { return grid_; }

  /// The kinematic + scattered-electron window `project_rates` applies at
  /// cell centers, evaluated event by event.
  bool in_acceptance(double x, double q2) const;

  // --- per-(spin state) amplitude tables ----------------------------------

  /// Cell-wise (w_avg, a1, a2) plus the sampling tables derived from them.
  struct StateTables {
    std::vector<double> w_avg, a1, a2;   ///< per accepted cell
    std::vector<double> a1n, a2n;        ///< divided by (1 + w_avg)
    std::vector<double> bound;           ///< 1 + |a1n| + |a2n|
    std::vector<double> cdf;             ///< cumulative sigma*(1 + w_avg)
    double sigma_pb = 0.0;               ///< the cdf's total
    double margin = 0.0;                 ///< exact min over phi of W/(1+w_avg)
  };

  /// Cached amplitudes for one category's spin configuration at projection m.
  /// Throws std::runtime_error if the phi density can go negative anywhere on
  /// the grid: the accept-reject draws u in [0, bound) and accepts on
  /// u < density, so a negative density would be silently sampled as
  /// max(W, 0) -- the modulation comes out diluted AND the (x, Q2) mixture is
  /// skewed, because the cell weights assume the per-cell phi integral is
  /// 2 pi (1 + w_avg).
  const StateTables& state_tables(const SpinCategory& cat, double m) const;

  // --- expected rates -----------------------------------------------------

  /// Accepted cross section [pb] of one pure spin state m.
  double sigma_state_pb(const SpinCategory& cat, double m) const {
    return state_tables(cat, m).sigma_pb;
  }
  /// Accepted cross section [pb] of the category's population mixture.
  double sigma_tot_pb(const SpinCategory& cat) const;
  double expected_events(const SpinCategory& cat, double lumi_pb) const {
    return lumi_pb * sigma_tot_pb(cat);
  }

  /// Fill-averaged accepted cross section and phi amplitudes:
  ///   sigma_pb = sum_m p_m sum_c sigma_c (1 + w_avg)
  ///   a_k      = [sum_m p_m sum_c sigma_c a_k] / sigma_pb
  /// so the category's phi' distribution is proportional to
  /// 1 + a1 cos phi' + a2 cos 2phi'.  `mask` (optional, one entry per
  /// accepted cell) restricts the sum to a super-bin in (x, Q2).
  struct EffectiveModulation {
    double sigma_pb = 0.0, a1 = 0.0, a2 = 0.0;
  };
  EffectiveModulation effective_modulation(
      const SpinCategory& cat, const std::vector<char>* mask = nullptr) const;

  // --- sampling -----------------------------------------------------------

  /// Everything the per-event draw needs for one category, resolved once.
  struct CategoryPlan {
    std::vector<double> m_val;              ///< m = +J ... -J
    std::vector<double> m_cdf;              ///< cumulative p_m * sigma_m
    std::vector<const StateTables*> states; ///< aligned with m_val
    double phi_s = 0.0;
    double sigma_pb = 0.0;                  ///< sum_m p_m sigma_m
  };
  CategoryPlan make_plan(const SpinCategory& cat) const;

  /// One event from one counter-based stream.  Draw order:
  /// m, cell, (x, Q2) [+ redraws], phi.
  EventDraw draw_event(const CategoryPlan& plan, Rng& rng) const;

  /// `n` events of one category on the (seed, run, bunch) stream, event
  /// counters `event0 .. event0+n-1`.  `nthreads > 1` splits the counter
  /// range over std::thread; the result is bit-identical to nthreads == 1.
  EventBatch sample_n(const SpinCategory& cat, std::size_t n,
                      std::uint64_t seed, std::uint64_t run,
                      std::uint64_t bunch, std::uint64_t event0 = 0,
                      unsigned nthreads = 1) const;

  /// Events for an integrated luminosity.  The count is Poisson with mean
  /// `lumi_pb * sigma_tot_pb(cat)` (independent Poissons per spin state are
  /// the same thing as one Poisson total plus a multinomial over m, which is
  /// what the per-event m draw does), from the reserved counter
  /// `kCountStreamEvent` so that adding events never moves the count.
  /// `poisson = false` rounds the mean instead.
  EventBatch sample_lumi(const SpinCategory& cat, double lumi_pb,
                         std::uint64_t seed, std::uint64_t run,
                         std::uint64_t bunch, bool poisson = true,
                         unsigned nthreads = 1) const;

  // --- Mode-W style weights ------------------------------------------------

  /// Per-event weight matrix w[i, k] = W(event_i | category_k), the
  /// polarized/unpolarized density ratio at the event's cell (mean over the
  /// fill populations).  Row-major, `batch.size()` rows by `cats.size()`
  /// columns.  This is the reweighting kernel Mode W (step 5.C) applies to
  /// external unpolarized samples.
  std::vector<double> weights_for(const EventBatch& batch,
                                  const std::vector<SpinCategory>& cats) const;

 private:
  struct StateKey {
    int lam_e;
    double pe, j, m, theta_s, phi_s;
    bool operator<(const StateKey& o) const;
  };
  const StateTables& build_state(const StateKey& key) const;

  std::shared_ptr<const InclusiveKernel> kernel_;
  BeamConfig config_;
  Scenario scenario_;
  GridSpec spec_;
  bool with_perp_;
  double s_ = 0.0;
  LogGrid grid_;
  std::vector<double> x_cells_, q2_cells_, xsec_flat_;
  std::vector<double> logx_lo_, logx_hi_, logq2_lo_, logq2_hi_;
  std::vector<SFTables> tables_;

  mutable std::mutex cache_mutex_;
  mutable std::map<StateKey, StateTables> cache_;
};

/// One pseudo-experiment: per-category batches plus the luminosities the
/// bookkeeper hands the analysis.  Category k uses bunch index k of the
/// (seed, run) stream, so two categories never share a counter.
struct PseudoExperiment {
  std::vector<std::string> names;  ///< aligned with plan.categories()
  std::vector<EventBatch> batches;
  std::vector<double> lumi_pb;

  /// Throws std::runtime_error if the name is not in the plan.
  const EventBatch& operator[](const std::string& name) const;
  double lumi_of(const std::string& name) const;
  /// Event count of a named category.
  std::size_t count(const std::string& name) const;
};

PseudoExperiment run_pseudo_experiment(const InclusiveSampler& sampler,
                                       const RunPlan& plan,
                                       double total_lumi_pb,
                                       std::uint64_t seed,
                                       std::uint64_t run = 1,
                                       bool poisson = true,
                                       unsigned nthreads = 1);

/// Binned phi' pseudo-experiment at full projected statistics.
///
/// Expected counts per bin are the EXACT integral of
/// n_expected/(2 pi) * (1 + a1 cos phi' + a2 cos 2phi') over each of `nbins`
/// uniform bins, Poisson-fluctuated.  For binned estimators this carries
/// statistics identical to event-level sampling, so a 1e8-event projection
/// costs `nbins` Poisson draws instead of 1e8 rows.
///
/// Throws if 1 + a1 cos phi' + a2 cos 2phi' goes negative anywhere: this path
/// bypasses the sampler's own guard and would otherwise hand negative bin
/// means to the estimator.
struct PhiHistogram {
  std::vector<double> counts;
  std::vector<double> edges;
};
PhiHistogram phi_histogram_pseudo(double n_expected, double a2, int nbins = 36,
                                  Rng* rng = nullptr, double a1 = 0.0,
                                  bool poisson = true);

// --------------------------------------------------------------- estimators

/// Analysis-side counting / fit estimators for spin-labeled pseudo-experiments
/// (`evgen/polligen/estimators.py`).  These are what an experiment would
/// apply; the Step-5.A closure tests check that their spreads reproduce the
/// analytic error formulas of `asymmetries.hpp` and that their pulls are
/// unbiased -- including with luminosity-corrected yields when the run plan
/// carries relative-luminosity offsets.
namespace estimators {

/// Luminosity-normalized yields counts/lumis * mean(lumis).  An empty `lumis`
/// means the naive equal-share assumption (yields = counts).
std::vector<double> yields(const std::vector<double>& counts,
                           const std::vector<double>& lumis = {});

/// Helicity-flip estimator (y+ - y-)/(y+ + y-)/(pe pz).  Pass
/// (l_plus, l_minus) > 0 for the luminosity-corrected form.
double apar_flip(double n_plus, double n_minus, double pe, double pz,
                 double l_plus = 0.0, double l_minus = 0.0);

/// Tensor thirds estimator (y+ + y- - 2 y0)/(y+ + y- + y0)/pzz for the
/// (+pzz, +pzz, -2 pzz) fill pattern of `tensor_thirds_plan`.  `lumis` empty
/// = naive.
double azz_thirds(double n_plus, double n_minus, double n_zero, double pzz,
                  const std::vector<double>& lumis = {});

/// Moment estimator 2 <cos 2phi'>/pzz -- uniform acceptance ONLY.
double cos2phi_moment(const std::vector<double>& phi_prime, double pzz);

/// Binned least-squares amplitude of N(phi') = C (1 + A cos 2phi') from a
/// precomputed phi' histogram, robust to phi acceptance holes.
///
/// `acceptance`: optional predicate marking live phi values.  Each bin is
/// probed at `nsub` points across its FULL width and dropped unless ALL are
/// live -- a bin merely straddling a hole edge would enter the fit at full
/// weight with depleted counts and bias the amplitude (by more than the
/// gluonometry signal for unlucky hole placements).
/// Returns A/pzz with the finite-bin-width dilution sin(2w)/(2w)
/// (w = bin half-width) corrected.
///
/// Preconditions: uniform bins; the counts follow C(1 + A cos 2phi') -- a
/// residual cos phi' modulation is orthogonal only at full coverage.  The
/// gluonometry fills satisfy this by design (m-symmetric tensor populations
/// give a1_eff = 0 identically).
double cos2phi_fit_binned(const std::vector<double>& counts,
                          const std::vector<double>& edges, double pzz,
                          const std::function<bool(double)>& acceptance = nullptr,
                          int nsub = 9);

/// Event-level wrapper: histogram phi' on `nbins` uniform bins over [0, 2pi)
/// and fit the cos 2phi' amplitude.
double cos2phi_fit(const std::vector<double>& phi_prime, double pzz,
                   int nbins = 36,
                   const std::function<bool(double)>& acceptance = nullptr,
                   int nsub = 9);

/// Statistical error of the binned amplitude on n events over `nbins` uniform
/// full-coverage bins: sqrt(2/n)/pzz inflated by the same 1/dilution the fit
/// divides out (+1.1 % at 24 bins, +11 % at 8).  The bare sqrt(2/n)/pzz of
/// `err_cos2phi_amplitude` is the unbinned/analytic-FOM convention.
double cos2phi_fit_err(double n, double pzz, int nbins = 24);

/// (estimate - truth)/expected_err.
inline double pull(double estimate, double truth, double expected_err) {
  return (estimate - truth) / expected_err;
}

}  // namespace estimators

}  // namespace lipolgen

#endif  // LIPOLGEN_SAMPLER_HPP
