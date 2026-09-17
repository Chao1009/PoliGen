// SPDX-License-Identifier: GPL-3.0-or-later
#include "lipolgen/pipeline.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>

#include "lipolgen/b1_nuclear.hpp"
#include "lipolgen/constants.hpp"
#include "lipolgen/sf.hpp"
#include "lipolgen/spin.hpp"
#include "lipolgen/triton_sf.hpp"

namespace lipolgen {
namespace {

const double kNaN = std::nan("");

/// Exact integer key for a half-integer projection.
long long mkey2(double m) { return std::llround(2.0 * m); }

/// How many times a tagged draw may be repeated to land a TIMELIKE X.
/// 0.02 % of d-control draws need one repeat, none needs two; the cap is a
/// runaway guard, not a tuning knob.
constexpr int kTaggedMaxRedraw = 64;

/// Log-uniform (x, Q2) inside accepted cell `c`, redrawn until the event is
/// inside the window -- the same 21-try loop, in the same draw order, as
/// `InclusiveSampler::draw_event` (sampler.cpp).
void draw_in_cell(const InclusiveSampler& s, std::size_t c, Rng& rng,
                  double& x, double& q2) {
  bool ok = false;
  for (int tries = 0; tries < 21; ++tries) {
    const double u1 = rng.uniform();
    const double u2 = rng.uniform();
    x = std::exp(s.logx_lo()[c] + u1 * (s.logx_hi()[c] - s.logx_lo()[c]));
    q2 = std::exp(s.logq2_lo()[c] + u2 * (s.logq2_hi()[c] - s.logq2_lo()[c]));
    if (s.in_acceptance(x, q2)) { ok = true; break; }
  }
  if (!ok) {
    x = s.x_cells()[c];
    q2 = s.q2_cells()[c];
  }
}

std::size_t pick(const std::vector<double>& cdf, double u) {
  std::size_t i = static_cast<std::size_t>(
      std::lower_bound(cdf.begin(), cdf.end(), u) - cdf.begin());
  if (i >= cdf.size()) i = cdf.size() - 1;
  return i;
}

/// The `Kinematics` spectator-lab block of a far-forward fragment that did
/// NOT come out of `boost_spectator` (the coherent intact recoil).  The beam
/// boost is longitudinal, so (kx, ky) = (px, py) and kz is the exact inverse
/// boost of (E, pz) -- the same algebra `boost_spectator` runs forwards.
void fill_fragment_lab(Event& ev, const Particle& frag, int frag_a,
                       double p_per_nucleon) {
  const Particle* bi = ev.find(Role::BeamIon);
  const double pt = frag.p.pt(), plab = frag.p.p();
  ev.kin.spec_pt = pt;
  ev.kin.spec_theta = std::atan2(pt, frag.p.pz);
  ev.kin.spec_p_lab = plab;
  ev.kin.phi_spec = std::atan2(frag.p.py, frag.p.px);
  ev.kin.spec_kx = frag.p.px;
  ev.kin.spec_ky = frag.p.py;
  if (bi && bi->mass > 0.0) {
    const double gamma = bi->p.e / bi->mass;
    const double gbeta = bi->p.pz / bi->mass;
    ev.kin.spec_kz = gamma * frag.p.pz - gbeta * frag.p.e;
    ev.kin.spec_r = (frag.charge != 0.0 && bi->charge > 0.0)
                        ? (plab / frag.charge) / (bi->p.pz / bi->charge)
                        : kNaN;
  }
  ev.kin.spec_xl = (frag_a > 0 && p_per_nucleon > 0.0)
                       ? plab / (frag_a * p_per_nucleon)
                       : kNaN;
}

Vec4 sum_p(const std::vector<Particle>& v) {
  Vec4 s;
  for (const Particle& p : v) s = s + p.p;
  return s;
}

int find_role(const Event& ev, Role r) {
  for (std::size_t i = 0; i < ev.particles.size(); ++i) {
    if (ev.particles[i].role == r) return static_cast<int>(i);
  }
  return -1;
}

/// One row of the tabulated lithium tagging optics.
struct TaggingRow {
  const char* ion;
  const char* key;
  double sigma_h, sigma_v, lumi_fraction;
  const char* name;
};

// `polli_fastsim.farforward.tagging_optics_point(config)` -- the de-squeeze
// working point that maximises (tagged fraction) x (luminosity), with the
// dispersive angular smearing priced on EACH CONFIGURATION's own pot levers
// (the default since 2026-08-29).
const TaggingRow kTaggingPerConfig[] = {
  {"6Li", "5x41",   3.6304912778772464e-05, 3.7997787958413641e-04, 0.14665506513094445, "5x41 tagging (beta*_x x 46)"},
  {"6Li", "10x100", 1.9210026303918496e-05, 1.8000055169177715e-04, 0.078054153005126567, "10x100 tagging (beta*_x x 164)"},
  {"6Li", "18x275", 1.1747140156192393e-05, 9.1694064905087008e-05, 0.10579263907143123, "18x275 tagging (beta*_x x 89)"},
  {"7Li", "5x41",   3.2460516257702269e-05, 3.7994089451492735e-04, 0.12668141170555919, "5x41 tagging (beta*_x x 62)"},
  {"7Li", "10x100", 1.7857468391993752e-05, 1.7998303137488586e-04, 0.067423585290711935, "10x100 tagging (beta*_x x 220)"},
  {"7Li", "18x275", 1.1821474197437985e-05, 9.9013268802557072e-05, 0.098880102921711496, "18x275 tagging (beta*_x x 102)"},
};

// The same, priced with the 18x275 pot levers everywhere: the single-lever
// behaviour every tagging number published before 2026-08-29 was made with.
const TaggingRow kTaggingLegacy[] = {
  {"6Li", "5x41",   3.2765085114524507e-05, 3.7997787958413641e-04, 0.14178288087065763, "5x41 tagging (beta*_x x 50)"},
  {"6Li", "10x100", 1.6547385211948192e-05, 1.8000055169177715e-04, 0.075461032778545642, "10x100 tagging (beta*_x x 176)"},
  {"6Li", "18x275", 1.1747140156192393e-05, 9.1694064905087008e-05, 0.10579263907143123, "18x275 tagging (beta*_x x 89)"},
  {"7Li", "5x41",   2.8748264721914243e-05, 3.7994089451492735e-04, 0.12247279347862211, "5x41 tagging (beta*_x x 67)"},
  {"7Li", "10x100", 1.5065056501145624e-05, 1.7998303137488586e-04, 0.065183634486883898, "10x100 tagging (beta*_x x 235)"},
  {"7Li", "18x275", 1.1821474197437985e-05, 9.9013268802557072e-05, 0.098880102921711496, "18x275 tagging (beta*_x x 102)"},
};

// THE ONE SPIN TEST `default_inclusive_kernel` MAKES before it fills the
// rank-2 (tensor) slots.  Written once and read twice -- by the kernel
// factory and by `inclusive_rank2_is_empty` below -- so the loud zero and the
// thing it reports on cannot drift apart (docs/CONVENTIONS.md).
//
// It keys on the SPIN and not on the isotope on purpose: that is what makes an
// `--isotope d --channel inclusive` run pick up the 6Li rank-2 transfer, which
// is documented behaviour and a phase-F decision (OPEN_ITEMS_SOLUTIONS.md
// sec. 14, "Adjacent defects"), not something this function may quietly change.
bool kernel_fills_rank2(const Ion& ion) {
  return std::fabs(ion.spin - 1.0) < 1e-9;
}

// The spin at which `InclusiveKernel::tables` (src/core/xsec.cpp) opens a
// rank-2 sector that `kernel_fills_rank2` then leaves EMPTY: its dispatch is
// `spin 1 -> b1_func/b2_func/delta_func`, `spin 3/2 ->
// b1_32_func/b2_32_func/delta_32_func`, and any other spin gets no rank-2
// sector at all (`rank2 = false`), which is physics rather than a gap.
bool kernel_has_rank2_sector_unfilled(const Ion& ion) {
  return std::fabs(ion.spin - 1.5) < 1e-9;
}

// %g of a double for the report sentences below.
std::string fmt_g(double v) {
  char buf[32];
  std::snprintf(buf, sizeof(buf), "%g", v);
  return std::string(buf);
}

}  // namespace

// ============================================================ DIS adapter

bool InclusiveKinematicsSource::Key::operator<(const Key& o) const {
  if (m2 != o.m2) return m2 < o.m2;
  if (lam_e != o.lam_e) return lam_e < o.lam_e;
  return pe < o.pe;
}

InclusiveKinematicsSource::InclusiveKinematicsSource(
    std::shared_ptr<const InclusiveSampler> sampler, double s_channel)
    : sampler_(std::move(sampler)), s_channel_(s_channel) {
  if (!sampler_) throw std::runtime_error("InclusiveKinematicsSource: null sampler");
  ms_ = m_values(s_channel_);
}

const InclusiveSampler::CategoryPlan& InclusiveKinematicsSource::plan_for(
    double m_s, int lam_e, double pe) const {
  const Key key{mkey2(m_s), lam_e, pe};
  std::lock_guard<std::mutex> lock(mutex_);
  const auto it = plans_.find(key);
  if (it != plans_.end()) return it->second;

  // The PURE struck-cluster category: spin S_c, all the population on m_S,
  // quantization axis along z -- the BEAM axis, NOT the ion fill's (the
  // Python's `_pure_category` does not pass the fill's axis either).  The
  // fill's axis reaches the event only through the spectator rotation in
  // `boost_spectator`, so this category is the right one exactly as long as
  // the struck-cluster DIS kernel is SPIN-BLIND; `Pipeline`'s constructor
  // refuses a tilted fill on a tagged channel when it is not (P_e != 0 on any
  // of the three, or `StruckClusterOptions::inclusive_b1` on the one channel
  // that reads it).  See the class comment in pipeline.hpp for the
  // measurement and for the alternative.
  SpinCategory cat;
  cat.j = s_channel_;
  cat.lam_e = lam_e;
  cat.pe = pe;
  cat.populations.assign(ms_.size(), 0.0);
  bool found = false;
  for (std::size_t i = 0; i < ms_.size(); ++i) {
    if (std::fabs(ms_[i] - m_s) < 1e-9) { cat.populations[i] = 1.0; found = true; }
  }
  if (!found) throw std::runtime_error("InclusiveKinematicsSource: unknown m_S");
  char buf[64];
  std::snprintf(buf, sizeof(buf), "struck m=%g lam=%d", m_s, lam_e);
  cat.name = buf;
  return plans_.emplace(key, sampler_->make_plan(cat)).first->second;
}

void InclusiveKinematicsSource::warm(int lam_e, double pe) const {
  for (double m : ms_) (void)plan_for(m, lam_e, pe);
}

double InclusiveKinematicsSource::sigma_tot_pb(double m_struck, int lam_e,
                                               double pe) const {
  return plan_for(m_struck, lam_e, pe).sigma_pb;
}

void InclusiveKinematicsSource::sample(double m_struck, int lam_e, double pe,
                                       std::size_t n, Rng& rng,
                                       std::vector<double>& x,
                                       std::vector<double>& q2,
                                       std::vector<double>& y,
                                       std::vector<double>& phi) {
  static thread_local std::vector<int> cell;
  sample_cells(m_struck, lam_e, pe, n, rng, x, q2, y, phi, cell);
}

void InclusiveKinematicsSource::sample_cells(double m_struck, int lam_e,
                                             double pe, std::size_t n, Rng& rng,
                                             std::vector<double>& x,
                                             std::vector<double>& q2,
                                             std::vector<double>& y,
                                             std::vector<double>& phi,
                                             std::vector<int>& cell) {
  const InclusiveSampler::CategoryPlan& plan = plan_for(m_struck, lam_e, pe);
  x.resize(n); q2.resize(n); y.resize(n); phi.resize(n); cell.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    const EventDraw d = sampler_->draw_event(plan, rng);
    x[i] = d.x; q2[i] = d.q2; y[i] = d.y; phi[i] = d.phi; cell[i] = d.cell;
  }
}

Vec4 InclusiveKinematicsSource::scattered_electron_p4(double x, double y,
                                                      double phi) const {
  const double e_e = sampler_->config().electron_energy;
  const ScatteredElectron e = scattered_electron(x, y, sampler_->s(), e_e);
  const double st = std::sin(e.theta), ct = std::cos(e.theta);
  return {e.e_prime, e.e_prime * st * std::cos(phi),
          e.e_prime * st * std::sin(phi), e.e_prime * ct};
}

std::shared_ptr<const InclusiveKernel> struck_cluster_kernel(
    const TaggedChannel& channel, const StruckClusterOptions& opt) {
  InclusiveKernel::Options o;
  if (opt.inclusive_b1) {
    // DELIBERATELY `toy_b1` AND NOTHING ELSE -- there is no `B1Model` here,
    // and design_D_b1_li6.md sec. 4.1 decided there must not be one.
    // `Li6ConvolutionB1` models the b1 of the WHOLE 6Li as seen INCLUSIVELY:
    // the z-smearing of its embedded-deuteron term and the orbital alignment
    // of its alpha-d D-wave terms are both properties of the alpha-d relative
    // motion.  On a tagged channel `TaggedSampler` already draws the alpha-d
    // momentum AND its m-dependent angular correlation from n_M(k, k_hat)
    // (tagged.hpp), so that same wave function is in the event weight
    // already -- which is the double counting `inclusive_b1` warns about
    // just above.  Putting the convolution backend here would count it
    // THREE times (the sampler's k_hat correlation, the orbital term, and
    // the z-smearing of the same density).  The guard is
    // `PipelineConfig::validate()`, which is the only route from the CLI; a
    // caller who hand-builds an `InclusiveKernel` and feeds it to a tagged
    // sampler is out of the library's reach, exactly as they are today for
    // any other hand-built kernel.
    o.b1_func = [](double x, double q2, double f1) {
      return toy_b1(x, q2, f1);
    };
  }
  o.delta_func = opt.delta_func;
  // The tagged channels' ONLY structure-function injection point (D2).
  // BOTH SLOTS ARE NULL BY DEFAULT and are handed over unchanged, so the
  // null branch of `InclusiveKernel`'s constructor -- a fresh `ToyF2` and a
  // `ToyG1` built on it with this kernel's own `r_func` -- is taken exactly
  // as it was before these fields existed.  A default tagged run is
  // therefore bit for bit BY CONSTRUCTION.  `Pipeline` fills them from
  // `PipelineConfig::unpol_sf_obj` / `pol_sf_obj`; note that
  // `PipelineConfig::kernel` does NOT reach here (the `Pipeline` reads it on
  // the non-tagged branch alone), which is why these two fields, and not
  // that escape hatch, are what a tagged run has.
  o.f2_source = opt.f2_source;
  o.g1_model = opt.g1_model;
  return std::make_shared<InclusiveKernel>(channel.dis_target, o);
}

std::shared_ptr<InclusiveKinematicsSource> make_struck_cluster_source(
    const TaggedChannel& channel, const BeamConfig& ion_config,
    const StruckClusterOptions& opt) {
  // The DIS beam configuration: the same electron energy and the same
  // PER-NUCLEON ion momentum, with the struck cluster as the target species
  // (the Python's `beams.BeamConfig(E_e, ch.dis_target, p_u)`).  s = 4 E_e p_u
  // is species-independent, so the DIS phase space is the ion's; what the
  // target species changes is Z/N in F2A and the effective polarizations in
  // g1A.
  BeamConfig dis = ion_config;
  dis.ion = channel.dis_target;
  const auto sampler = std::make_shared<InclusiveSampler>(
      struck_cluster_kernel(channel, opt), dis, opt.scenario, opt.grid,
      opt.with_perp);
  return std::make_shared<InclusiveKinematicsSource>(sampler,
                                                     channel.s_channel);
}

// ============================================================ configuration

const char* pipeline_channel_name(PipelineChannel c) {
  switch (c) {
    case PipelineChannel::Inclusive:       return "inclusive";
    case PipelineChannel::TaggedLi6Alpha:  return "tagged-6Li-alpha";
    case PipelineChannel::TaggedLi7Alpha:  return "tagged-7Li-alpha";
    case PipelineChannel::TaggedDeuteronP: return "tagged-d-p";
    case PipelineChannel::CoherentLi6:     return "coherent-6Li";
  }
  return "unknown";
}

const char* b1_model_name(B1Model m) {
  switch (m) {
    case B1Model::Miller:         return "miller";
    case B1Model::Cdks:           return "cdks";
    case B1Model::Li6Convolution: return "li6-convolution";
  }
  return "unknown";
}

const char* triton_sf_name(TritonSfChoice s) {
  switch (s) {
    case TritonSfChoice::Hulthen: return "hulthen";
    case TritonSfChoice::CiofiSimula: return "ciofi-simula";
  }
  return "?";
}

const char* knob_status_name(KnobStatus s) {
  switch (s) {
    case KnobStatus::Read: return "read";
    case KnobStatus::NotRead: return "not-read";
    case KnobStatus::Refused: return "refused";
  }
  return "?";
}

const char* b1_unpol_name(B1UnpolSource s) {
  switch (s) {
    case B1UnpolSource::Toy:     return "toy";
    case B1UnpolSource::Mstw:    return "mstw";
    case B1UnpolSource::Ct18Nlo: return "ct18nlo";
    case B1UnpolSource::Custom:  return "custom";
  }
  return "unknown";
}

const char* r_source_name(RSource s) {
  switch (s) {
    case RSource::Unset:   return "unset";
    case RSource::SigmaLt: return "sigma-lt";
    case RSource::R1998:   return "r1998";
  }
  return "unknown";
}

const char* unpol_sf_name(UnpolSfSource s) {
  switch (s) {
    case UnpolSfSource::Toy:     return "toy";
    case UnpolSfSource::Mstw:    return "mstw";
    case UnpolSfSource::Ct18Nlo: return "ct18nlo";
    case UnpolSfSource::Custom:  return "custom";
  }
  return "unknown";
}

const char* pol_sf_name(PolSfSource s) {
  switch (s) {
    case PolSfSource::Toy:      return "toy";
    case PolSfSource::NnpdfPol: return "nnpdfpol";
    case PolSfSource::Custom:   return "custom";
  }
  return "unknown";
}

bool is_tagged(PipelineChannel c) {
  return c == PipelineChannel::TaggedLi6Alpha ||
         c == PipelineChannel::TaggedLi7Alpha ||
         c == PipelineChannel::TaggedDeuteronP;
}

std::string channel_isotope(PipelineChannel c) {
  switch (c) {
    case PipelineChannel::TaggedLi6Alpha:  return "6Li";
    case PipelineChannel::TaggedLi7Alpha:  return "7Li";
    case PipelineChannel::TaggedDeuteronP: return "d";
    case PipelineChannel::CoherentLi6:     return "6Li";
    case PipelineChannel::Inclusive:       break;
  }
  return std::string();
}

Optics tagging_optics(const std::string& ion_name, double p_per_nucleon,
                      bool legacy_levers, double n_sigma) {
  const std::string key = yr_config_key(ion_name, p_per_nucleon);
  const TaggingRow* table = legacy_levers ? kTaggingLegacy : kTaggingPerConfig;
  const std::size_t n = 6;
  for (std::size_t i = 0; i < n; ++i) {
    if (ion_name == table[i].ion && key == table[i].key) {
      Optics o;
      o.name = table[i].name;
      o.sigma_theta = table[i].sigma_h;
      o.sigma_theta_v = table[i].sigma_v;
      o.n_sigma = n_sigma;
      o.lumi_fraction = table[i].lumi_fraction;
      return o;
    }
  }
  throw std::runtime_error("tagging_optics: no tabulated working point for " +
                           ion_name + " at " + key +
                           " (the scan is only published for 6Li and 7Li)");
}

Optics optics_for(OpticsChoice choice, const std::string& ion_name,
                  double p_per_nucleon, double n_sigma) {
  switch (choice) {
    case OpticsChoice::YellowReportHighAcceptance:
      return yr_optics(ion_name, p_per_nucleon, true, n_sigma);
    case OpticsChoice::YellowReportHighDivergence:
      return yr_optics(ion_name, p_per_nucleon, false, n_sigma);
    case OpticsChoice::Tagging:
      return tagging_optics(ion_name, p_per_nucleon, false, n_sigma);
    case OpticsChoice::TaggingLegacyLevers:
      return tagging_optics(ion_name, p_per_nucleon, true, n_sigma);
    case OpticsChoice::Custom:
      break;
  }
  throw std::runtime_error("optics_for: OpticsChoice::Custom has no table");
}

// rc.hpp is keyed on `Channel` (event.hpp) and CANNOT include pipeline.hpp --
// pipeline.hpp includes rc.hpp for `PipelineRc` and `RcOptions`.  So the
// 5 -> 9 widening happens here, once (design_C_tensor_rc.md sec. 3.3).
static Channel event_channel_of(PipelineChannel c) {
  switch (c) {
    case PipelineChannel::Inclusive:       return Channel::Inclusive;
    case PipelineChannel::TaggedLi6Alpha:  return Channel::TaggedLi6Alpha;
    case PipelineChannel::TaggedLi7Alpha:  return Channel::TaggedLi7Alpha;
    case PipelineChannel::TaggedDeuteronP: return Channel::TaggedDeuteronP;
    case PipelineChannel::CoherentLi6:     return Channel::CoherentLi6;
  }
  throw std::runtime_error("event_channel_of: unhandled PipelineChannel");
}

void PipelineConfig::validate() const {
  if (beam_config < 0 || beam_config > 2) {
    throw std::runtime_error("PipelineConfig: beam_config must be 0, 1 or 2");
  }
  if (lumi_pb > 0.0 && n_events > 0) {
    throw std::runtime_error("PipelineConfig: lumi_pb and n_events are exclusive");
  }
  if (!(lumi_pb > 0.0) && n_events == 0) {
    throw std::runtime_error("PipelineConfig: set either lumi_pb or n_events");
  }
  const std::string want = channel_isotope(channel);
  if (!want.empty() && want != isotope) {
    throw std::runtime_error(std::string("PipelineConfig: channel ") +
                             pipeline_channel_name(channel) + " needs isotope " +
                             want + ", got " + isotope);
  }
  if (fsi != PipelineFsi::Off && !is_tagged(channel)) {
    throw std::runtime_error(
        std::string("PipelineConfig: fsi = ") + pipeline_fsi_name(fsi) +
        " needs a tagged channel (the weight is a distortion of the tagged "
        "spectator's spectrum; " + pipeline_channel_name(channel) +
        " has no tagged spectator)");
  }
  if (!(fsi_sigma_mb >= 0.0)) {
    throw std::runtime_error("PipelineConfig: fsi_sigma_mb must be >= 0");
  }
  // C5.5b.  THE SAME RULE, ONE LEVEL UP: `cluster_wave` itself is read ONLY
  // where a channel has cluster relative waves, i.e. on the three tagged
  // channels (`Pipeline`'s `is_tagged` branch is the only reader).  On
  // `Inclusive` and `CoherentLi6` it was ACCEPTED and never read -- exactly the
  // "silently ignored flag" C5.4 named as a defect on the deuteron control,
  // and the mechanism by which C5.5's 11.61 % inclusive-vs-tagged drift used to
  // reach a user: they set the flag wanting a "VMC 6Li" inclusive row, nothing
  // refused it, and they got `LI6_CLUSTER_POLARIZATION` = 0.811228 with no
  // warning.  The asymmetry with `cluster_vmc_mc_sigma` (refused for a 0.02 %
  // effect two lines below) had no defence: this one is worth 11.61 %.
  //
  // WHAT THE REFUSAL DOES NOT CLAIM.  It does not make an inclusive VMC 6Li
  // available -- there is none, and C5.5 is the argument that substituting
  // `VMC_P_D_LI6` into `LI6_CLUSTER_POLARIZATION` would move it AWAY from the
  // one ab initio anchor.  It converts a silent 11.61 % into an error message.
  if (cluster_wave != ClusterWaveSource::Hulthen && !is_tagged(channel)) {
    throw std::runtime_error(
        std::string("PipelineConfig: cluster_wave selects the CLUSTER RELATIVE "
                    "wave function of a tagged channel and is never read on ")
        + pipeline_channel_name(channel)
        + " -- the inclusive 6Li effective polarization "
          "(LI6_CLUSTER_POLARIZATION = 0.811228) and the coherent form factors "
          "do not follow it, so accepting it would promise a VMC row and "
          "deliver the shipped one, 11.61 % away in the vector sector (open "
          "item C5.5).  Leave it at Hulthen, or run a tagged channel");
  }
  // C5.2.  A KNOB THAT DID NOT RUN MAY NOT BE RECORDED AS IF IT HAD -- the
  // same rule the b1 band scales are refused under.  The MC band exists only
  // where a table carries printed 1-sigma errors, which is the two `momenta/`
  // files and nothing else: the Hulthen forms are analytic and `fdeut.av18`
  // (the deuteron control's VmcAV18 wave) prints none.
  if (cluster_vmc_mc_sigma != 0.0) {
    if (cluster_wave != ClusterWaveSource::VmcAV18) {
      throw std::runtime_error(
          "PipelineConfig: cluster_vmc_mc_sigma is the ANL VMC tables' own "
          "Monte Carlo band and needs cluster_wave = VmcAV18 (the analytic "
          "Hulthen forms carry no MC error)");
    }
    if (channel != PipelineChannel::TaggedLi6Alpha
        && channel != PipelineChannel::TaggedLi7Alpha) {
      throw std::runtime_error(
          std::string("PipelineConfig: cluster_vmc_mc_sigma needs a lithium "
                      "alpha-tag channel, got ")
          + pipeline_channel_name(channel)
          + " (only li6_ad1.momentum and li7_at3.momentum print MC errors; "
            "the deuteron control's fdeut.av18 does not, so the band would be "
            "recorded without having run)");
    }
  }
  // rc.hpp.  NOTE what is deliberately NOT here: a channel refusal.  Unlike
  // FSI, `--rc` is legal on every channel -- `RcModel::applies()` returns
  // false on `CoherentLi6`, the weights are exactly 1.0 and the run PRINTS
  // the reason, so a scan over channels does not have to special-case it
  // (design_C_tensor_rc.md sec. 1.5.3, sec. 3.2).
  if (rc != PipelineRc::Off) {
    if (!(rc_options.delta_low_x >= 0.0 && rc_options.delta_high_x >= 0.0)) {
      throw std::runtime_error("PipelineConfig: rc delta must be >= 0");
    }
    if (!(rc_options.x_low > 0.0 && rc_options.x_high > rc_options.x_low)) {
      throw std::runtime_error(
          "PipelineConfig: rc needs 0 < x_low < x_high (the band anchors)");
    }
    if (rc_options.n_eta < 8) {
      throw std::runtime_error("PipelineConfig: rc n_eta must be >= 8");
    }
    if (!(rc_options.fq_scale >= 0.0 && rc_options.tail_tensor_scale >= 0.0 &&
          rc_options.qe_suppression >= 0.0 &&
          rc_options.qe_tensor_scale >= 0.0 &&
          rc_options.sp_tensor_scale >= 0.0)) {
      throw std::runtime_error(
          "PipelineConfig: rc fq_scale / tail_tensor_scale / qe_suppression / "
          "qe_tensor_scale / sp_tensor_scale must be >= 0");
    }
    // Its sign is erased by the quadrature in `RcModel::delta`, so a negative
    // value would be recorded in `meta` as a price it did not charge.
    if (!(rc_options.a_transfer_frac >= 0.0)) {
      throw std::runtime_error(
          "PipelineConfig: rc a_transfer_frac must be >= 0 (it is added in "
          "QUADRATURE to delta(x), which erases its sign)");
    }
    // THE SAME "a knob that did not run may not be recorded" RULE as
    // `m_lepton` below.  `qe_tensor_scale` is a fraction OF Eq. (44)'s
    // sigma^q_U; with `with_qe_tail = false` there is no sigma^q_U, so the
    // npz `meta` would carry `rc_qe_tensor_scale` as if the polarised
    // quasi-elastic omission had been priced when nothing was computed.
    // `RcModel`'s constructor makes the identical check for the C++ API path.
    if (!rc_options.with_qe_tail && rc_options.qe_tensor_scale != 0.0) {
      throw std::runtime_error(
          "PipelineConfig: rc qe_tensor_scale != 0 needs with_qe_tail = true "
          "-- the polarised stand-in is a fraction of the UNPOLARISED "
          "quasi-elastic tail this run switches off, so it would record a "
          "price that was never paid (rc.hpp, RcOptions::qe_tensor_scale)");
    }
    // ... AND THE OTHER TWO KNOBS ON THE SAME TERM.  `qe_tensor_scale` was
    // the only member of this clause until 2026-09-05, and the two sub-knobs
    // BESIDE it were in exactly the situation it exists to refuse:
    // `qe_suppression` is a flat multiplier ON Eq. (44)'s sigma^q_U and
    // `qe_kf_gev` is the Fermi momentum of the Pauli suppression S(q) INSIDE
    // it, so with `with_qe_tail = false` neither is a factor of anything --
    // measured (inclusive 6Li, --rc tensor-band, 300 events seed 7)
    // `qe_suppression = 0.5` and `qe_kf_gev = 0.25` are BIT-IDENTICAL to
    // `with_qe_tail = false` alone, and `meta` recorded them as 0.5 and 0.25
    // with row status `read`, at_default false.  The criterion is the one
    // stated at `KnobProvenance`: a scale on a term the run computes as
    // identically 1 is REFUSED, not labelled.
    {
      const RcOptions od;
      const char* qe_knob =
          !rc_options.with_qe_tail
              ? (rc_options.qe_suppression != od.qe_suppression
                     ? "--rc-qe-suppression"
                     : (rc_options.qe_kf_gev != od.qe_kf_gev
                            ? "rc_options.qe_kf_gev"
                            : nullptr))
              : nullptr;
      if (qe_knob != nullptr) {
        throw std::runtime_error(
            std::string("PipelineConfig: ") + qe_knob +
            " is a knob on the UNPOLARISED QUASI-ELASTIC TAIL, and "
            "rc_options.with_qe_tail = false switches that tail off -- "
            "POLRAD Eq. (44)'s sigma^q_U is not computed at all, so the "
            "multiplier has nothing to multiply and the Fermi momentum "
            "nothing to suppress, and the setting is bit-identical to "
            "with_qe_tail = false alone.  Recording it would claim a "
            "systematic was PRICED that was never computed (the "
            "rc_qe_tensor_scale rule, one level up).  Drop the knob, or "
            "leave the quasi-elastic tail on");
      }
    }
    // `RcTailModel::PolradFull` -- POLRAD Eq. (18) + Appendix B + Eq. (A.4),
    // implemented 2026-09-06 (phase_B_numbers.md sec. B2).  Its two
    // preconditions are `RcModel`'s own and are repeated here for the same
    // reason every other rc check is: `validate()` has to refuse BEFORE any
    // model exists.
    if (rc_options.tail_model == RcTailModel::PolradFull) {
      if (std::numeric_limits<long double>::digits <=
          std::numeric_limits<double>::digits) {
        throw std::runtime_error(
            "PipelineConfig: rc tail_model polrad-full needs a `long double` "
            "wider than `double` (Eq. (B.3)'s a_ik cancel to ~5 decimal "
            "digits and the TENSOR column loses 14 % to it); this platform's "
            "is not.  Run --rc-tail-model t-peak or t-peak+ll");
      }
      if (rc_options.n_eta < 64) {
        throw std::runtime_error(
            "PipelineConfig: rc tail_model polrad-full needs n_eta >= 64 -- "
            "under it n_eta is the tanh-sinh node count PER PANEL of the "
            "tau_A quadrature, and sigma^el_U is 0.51 % low at 32 (measured, "
            "6Li config 1, x = 1e-3, y = 0.5)");
      }
    }
    // THE SAME RULE AGAIN, ONE LEVEL DOWN: `sp_tensor_scale` is a fraction OF
    // the leading-log s-/p-peak column `TailTriple::u_sp`, and only
    // `tail_model = TPeakPlusLL` computes that column -- under `TPeak` the
    // table is identically zero, so the scale multiplies nothing and the run
    // is bit-identical to the default while `meta` would carry
    // `rc_sp_tensor_scale` as if the s/p tensor omission had been priced.
    // `RcModel`'s constructor makes the identical check for the C++ API path.
    if (rc_options.tail_model != RcTailModel::TPeakPlusLL &&
        rc_options.sp_tensor_scale != 0.0) {
      throw std::runtime_error(
          std::string("PipelineConfig: rc sp_tensor_scale != 0 needs "
                      "tail_model = TPeakPlusLL (--rc-tail-model t-peak+ll) "
                      "-- the s/p tensor stand-in is a fraction of the "
                      "UNPOLARISED leading-log s-/p-peaks, and ") +
          (rc_options.tail_model == RcTailModel::PolradFull
               ? "polrad-full does not need one: Eq. (18) carries the s- and "
                 "p-peaks in its own tau_A integral WITH their Eq. (A.4) "
                 "tensor content, so the stand-in would double-count a term "
                 "that ran"
               : "the t-peak-only tail does not compute them, so it would "
                 "record a price that was never paid") +
          " (rc.hpp, RcOptions::sp_tensor_scale)");
    }
    // A KNOB THAT DID NOT RUN MAY NOT BE RECORDED AS IF IT HAD -- the same
    // rule the Miller branch enforces for `b1_band_scale` below.  `m_lepton`
    // is RESERVED for `RcTailModel::PolradFull` (F_IR, l_m).
    //
    // IT IS REFUSED ON `TPeakPlusLL` TOO, and the reason is not the same as on
    // `TPeak`.  `TPeak` carries no lepton mass at all.  `TPeakPlusLL` DOES --
    // `ll_radiator` is (alpha/pi) ln(Q^2/m_e^2)(1+z^2)/(1-z) -- but it reads
    // `M_ELECTRON` from constants.hpp directly and NOT `rc_options.m_lepton`,
    // deliberately: a muon radiator would need the muon's own elastic
    // kinematics (the `dsigma_el_dq2` kinematic factor carries the LEPTON mass
    // nowhere and the TARGET mass everywhere), so honouring `m_lepton` in the
    // log alone would be a half-change dressed as a whole one.  Refusing it
    // keeps the rule intact: what `meta` records is what ran.
    //
    // SINCE 2026-09-06 THE RESERVATION IS DISCHARGED ON ONE MODEL.
    // `PolradFull` reads `m_lepton` and nothing else does: it is the m^2 of
    // C_{1,2}(tau) (Eq. (B.13)), of F_IR = m^2 F_2+ - Q_m^2 F_d, of Q_m^2
    // and of lambda_s, i.e. it is what REGULATES the s- and p-peaks.  So the
    // row is three-way now -- read under polrad-full, refused on the two
    // t-peak models -- and the refusal below applies only to them.
    if (rc_options.tail_model != RcTailModel::PolradFull &&
        rc_options.m_lepton != M_ELECTRON) {
      throw std::runtime_error(
          "PipelineConfig: rc_options.m_lepton is read ONLY by tail_model = "
          "PolradFull (POLRAD's F_IR, C_{1,2} and lambda_s); the shipped "
          "TPeak tail has no lepton-mass dependence and TPeakPlusLL's "
          "leading-log radiator reads constants.hpp's M_ELECTRON directly, "
          "so setting it here would record a variation that did not run -- "
          "leave it at M_ELECTRON or run --rc-tail-model polrad-full");
    }
    if (!(rc_options.m_lepton > 0.0)) {
      throw std::runtime_error(
          "PipelineConfig: rc_options.m_lepton must be > 0 -- it is a mass "
          "squared in C_{1,2}(tau) and in lambda_s = S^2 - 4 m^2 M^2");
    }
    // A KNOB THAT DID NOT RUN MAY NOT BE RECORDED AS IF IT HAD, applied to
    // the rc sub-knobs AS A CLASS rather than one at a time.  The
    // `qe_tensor_scale` clause above is the precedent and was, until
    // 2026-09-05, the only member of it: the identical situation one level up
    // was accepted in silence.  MEASURED (600 events, seed 7, sha256 over all
    // 47 columns plus sigma_pb and sigma_per_category_pb, against the same
    // plan's `--rc tensor-band` baseline): on tagged-6Li-alpha,
    // tagged-7Li-alpha and tagged-d-p, where `rc_tail_applies` is false BY
    // CONSTRUCTION, `--rc-fq-scale 2`, `--rc-tail-tensor-scale 2`,
    // `--rc-qe-suppression 0.5`, `--rc-qe-tensor-scale 1`, `--rc-c0-shape
    // vmc-ft` and `--rc-tail-model t-peak+ll` are each BIT-IDENTICAL to that
    // baseline, and every one of them was written into the npz `meta` as a
    // number or a name; on the coherent channel, where `rc_applies` is false
    // too, so were `--rc-delta-low-x 0.19` and `--rc-a-transfer-frac 0.5`.
    //
    // WHY REFUSED AND NOT LABELLED -- the criterion is stated once at
    // `KnobProvenance` (pipeline.hpp).  Each of these names a VARIATION OF A
    // PIECE THAT DID NOT RUN: a scale, a band edge or a shape on a term this
    // run computes as identically 1.0.  Recording such a value claims a
    // systematic was PRICED, and no label makes a priced systematic
    // un-priced.  `--rc` ITSELF stays accepted on every channel (a channel
    // scan needs no special case, design sec. 1.5.3) and every sub-knob stays
    // available wherever its piece runs; what is refused is the combination.
    // `rc_band_applies` / `rc_tail_applies` are `RcModel`'s OWN two
    // predicates (rc.hpp), read here rather than re-derived.
    {
      const Channel ec = event_channel_of(channel);
      const RcOptions od;
      struct SubKnob { const char* flag; bool moved; };
      const SubKnob band_knobs[] = {
          {"--rc-delta-low-x", rc_options.delta_low_x != od.delta_low_x},
          {"--rc-delta-high-x", rc_options.delta_high_x != od.delta_high_x},
          {"rc_options.x_low", rc_options.x_low != od.x_low},
          {"rc_options.x_high", rc_options.x_high != od.x_high},
          {"--rc-a-transfer-frac",
           rc_options.a_transfer_frac != od.a_transfer_frac},
          {"rc_options.band_tau_max",
           rc_options.band_tau_max != od.band_tau_max},
          {"rc_options.scope", rc_options.scope != od.scope},
      };
      const SubKnob tail_knobs[] = {
          {"--rc-fq-scale", rc_options.fq_scale != od.fq_scale},
          {"--rc-tail-tensor-scale",
           rc_options.tail_tensor_scale != od.tail_tensor_scale},
          {"--rc-c0-shape", rc_options.c0_shape != od.c0_shape},
          {"--rc-qe-suppression",
           rc_options.qe_suppression != od.qe_suppression},
          {"--rc-qe-tensor-scale",
           rc_options.qe_tensor_scale != od.qe_tensor_scale},
          {"--rc-sp-tensor-scale",
           rc_options.sp_tensor_scale != od.sp_tensor_scale},
          {"rc_options.qe_kf_gev", rc_options.qe_kf_gev != od.qe_kf_gev},
          {"--rc-tail-model", rc_options.tail_model != od.tail_model},
          {"rc_options.with_qe_tail",
           rc_options.with_qe_tail != od.with_qe_tail},
          {"rc_options.n_eta", rc_options.n_eta != od.n_eta},
          {"rc_options.tail_max", rc_options.tail_max != od.tail_max},
      };
      if (!rc_band_applies(ec)) {
        for (const SubKnob& k : band_knobs) {
          if (!k.moved) continue;
          throw std::runtime_error(
              std::string("PipelineConfig: ") + k.flag +
              " is a knob on the tensor-RC BAND, and the band does not apply "
              "on channel " + pipeline_channel_name(channel) +
              " -- RcModel::applies() is false there (its tensor dependence "
              "is entirely AZIMUTHAL and nobody has computed RC for a "
              "phi-dependent tensor observable), every rc_tensor_* weight is "
              "exactly 1.0, and the setting is bit-identical to the default. "
              " Recording it would claim a systematic was PRICED that was "
              "never computed (the rc_qe_tensor_scale rule).  Drop the knob "
              "-- `--rc tensor-band` itself stays accepted here and prints "
              "why it prices nothing -- or run it on a channel where the "
              "band applies");
        }
      }
      if (!rc_tail_applies(ec, rc_options.with_tail)) {
        for (const SubKnob& k : tail_knobs) {
          if (!k.moved) continue;
          throw std::runtime_error(
              std::string("PipelineConfig: ") + k.flag +
              " is a knob on the RADIATIVE TAIL, and rc_tail == 1 exactly on "
              "this run (channel " + pipeline_channel_name(channel) +
              (rc_options.with_tail
                   ? std::string("")
                   : std::string(", rc_options.with_tail = false")) +
              "), so there is nothing for it to scale and the setting is "
              "bit-identical to the default -- measured on all three tagged "
              "channels.  Recording it would claim a systematic was PRICED "
              "that was never computed (the rc_qe_tensor_scale rule, which is "
              "this rule's own first member).  Drop the knob, or run the "
              "inclusive channel, where the tail applies");
        }
      }
    }
    // DESIGN vs CODE: sec. 3.2's snippet also loops over `plan.categories()`
    // here to refuse a polarised beam.  `PipelineConfig` has no run plan --
    // `RunPlan` is the Pipeline's SECOND constructor argument, not a config
    // field -- so that loop cannot live in this function.  `RcModel`'s own
    // constructor makes exactly that check (rc.cpp, "--rc assumes an
    // UNPOLARISED beam"), which is the single source of truth the header
    // already promises; the run therefore still throws, one frame later.
  }
  // design_D_b1_li6.md sec. 4.1/4.2.  The rule KEYS ON THE ISOTOPE, not on
  // the spin: the inclusive channel accepts isotope == "d", which is spin 1,
  // so a spin test would let a deuteron beam run with
  // `--b1-model li6-convolution` and silently apply the 6Li alpha-d
  // convolution -- N_ad, the alpha-d densities, the 2/6 and 4/6 counting
  // factors -- to a deuteron.
  //
  // The two knobs are checked for EVERY model, not only inside the non-Miller
  // branch: `--b1-band-scale -1 --b1-model miller` used to pass because the
  // range check sat inside the branch.
  if (!(b1_band_scale >= 0.0)) {
    throw std::runtime_error("PipelineConfig: b1_band_scale must be >= 0");
  }
  if (!(b1_alpha_d_dwave_weight >= 0.0)) {
    throw std::runtime_error(
        "PipelineConfig: b1_alpha_d_dwave_weight must be >= 0");
  }
  if (b1_model == B1Model::Miller) {
    // A KNOB THAT DID NOT RUN MAY NOT BE RECORDED AS IF IT HAD.  Neither
    // scale reaches the Miller branch of `default_inclusive_kernel` -- the
    // band is deliberately not applied to the published numbers and the
    // alpha-d weight belongs to `Li6ConvolutionB1` -- but the npz/HFS `meta`
    // records both unconditionally, so a run with either set would claim a
    // variation it never made.  Refusing is the same rule as the
    // kernel/flag conflict below (design 4.2's provenance promise).
    if (b1_band_scale != 1.0 || b1_alpha_d_dwave_weight != 1.0) {
      throw std::runtime_error(
          "PipelineConfig: --b1-band-scale / --b1-alpha-d-dwave-weight do not "
          "apply to b1_model = miller (the band is deliberately not applied "
          "to the published numbers, and the alpha-d D-wave knob belongs to "
          "li6-convolution), so setting either would record a variation that "
          "did not run -- drop them, or choose another b1_model");
    }
  } else {
    // Everything below is the rule for BOTH opt-in models.  `Cdks` goes
    // through `Li6B1`, i.e. the 6Li rank-2 transfer (2/6 x
    // LI6_B1_RANK2_TRANSFER), exactly as `Li6Convolution` does, so
    // "inclusive only, 6Li only" is as true of it as of the convolution --
    // and on any other channel or isotope `default_inclusive_kernel` never
    // fills the b1 slot from the flag at all, while the metadata would still
    // say it did.
    if (kernel) {
      throw std::runtime_error(
          std::string("PipelineConfig: b1_model = ") + b1_model_name(b1_model) +
          " with a caller-supplied kernel: the kernel wins and the flag would "
          "be silently ignored, so the two are refused together (drop one)");
    }
    if (channel != PipelineChannel::Inclusive) {
      throw std::runtime_error(
          std::string("PipelineConfig: b1_model = ") + b1_model_name(b1_model) +
          " needs the inclusive channel, got " +
          pipeline_channel_name(channel) +
          " (on a tagged channel the alpha-d density is already in the event "
          "weight; on the coherent channel the tensor signal is in the recoil "
          "azimuth and the inclusive b1 is not folded in at all -- either way "
          "the flag would not reach the rate but WOULD reach the metadata)");
    }
    if (isotope == "d") {
      throw std::runtime_error(
          std::string("PipelineConfig: b1_model = ") + b1_model_name(b1_model) +
          " is 6Li ONLY, got isotope d -- the deuteron is spin 1 too, but "
          "both opt-in models carry the 6Li rank-2 transfer (li6-convolution "
          "the alpha-d densities, N_ad and the 2/6 and 4/6 counting factors; "
          "cdks the Li6B1 2/6 transfer).  For the A = 2 kernel use "
          "DeuteronConvolutionB1 directly (b1_nuclear.hpp); that is the "
          "validation gate, not a beam species");
    }
    if (isotope != "6Li") {
      throw std::runtime_error(
          std::string("PipelineConfig: b1_model = ") + b1_model_name(b1_model) +
          " is 6Li ONLY, got isotope " + isotope +
          (isotope == "7Li"
               ? " -- spin 3/2 has no rank-2 input here (the 7Li rank-2 slots "
                 "are empty by design; there is no published b1 for it).  "
                 "AND NOTE WHAT DROPPING THE FLAG DOES NOT BUY: a 7Li "
                 "inclusive run's tensor and cos 2phi sector is identically "
                 "zero with or without it -- the run now says so in its "
                 "banner and in meta[\"rank2_input\"] "
                 "(inclusive_rank2_is_empty, docs/USAGE.md sec. 2c)"
               : ""));
    }
    if (b1_model == B1Model::Cdks && b1_alpha_d_dwave_weight != 1.0) {
      throw std::runtime_error(
          "PipelineConfig: --b1-alpha-d-dwave-weight is a knob on terms (2d) "
          "and (2a) of li6-convolution; the cdks branch never reads it, so "
          "setting it with b1_model = cdks would record a variation that did "
          "not run");
    }
  }
  // --b1-unpol: the UNPOLARISED backend the `Li6Convolution` b1 folds its own
  // F1 against.  Two rules, both of them the ones already enforced above.
  //
  // (1) A KNOB THAT DID NOT RUN MAY NOT BE RECORDED AS IF IT HAD.  Only the
  //     `Li6Convolution` branch of `default_inclusive_kernel` reads it
  //     (`Miller` is a ratio model with no F1 in it and `Cdks` carries the
  //     digitized column), but `meta["b1_unpol"]` is written unconditionally,
  //     so any other model would claim a PDF that never touched the rate.
  //     The caller-supplied-kernel case needs no separate clause: a non-Toy
  //     backend requires b1_model = li6-convolution, and THAT is already
  //     refused together with `kernel` a few lines above.  The guard is not
  //     loosened anywhere -- the point of this flag is that the passing
  //     configuration no longer needs a hand-built kernel to be expressed.
  //
  // (2) A NAMED BACKEND WITH AN EMPTY SLOT IS AN ERROR, NEVER A FALLBACK.
  //     `MstwSF` lives in the optional PYTHIA tier and `LhapdfSF` in the
  //     optional LHAPDF tier; `sf.hpp`'s rule is that the core library links
  //     neither, so this function cannot build either one.  Falling back to
  //     `ToyF2` here would reproduce exactly the defect this flag exists to
  //     fix -- a run labelled `mstw` whose numbers are the `toy` ones -- so
  //     it throws instead, and says which tier is missing.
  if (b1_unpol == B1UnpolSource::Toy) {
    if (b1_unpol_sf) {
      throw std::runtime_error(
          "PipelineConfig: b1_unpol = toy with an object in b1_unpol_sf -- "
          "the toy setting IS the kernel's own ToyF2, shared as one object, "
          "so the attached backend would be silently dropped.  Use "
          "B1UnpolSource::Custom to fold against your own UnpolSF (Python: "
          "config.b1_unpol_sf = obj, which sets Custom for you), or clear "
          "the slot");
    }
  } else {
    if (b1_model != B1Model::Li6Convolution) {
      throw std::runtime_error(
          std::string("PipelineConfig: b1_unpol = ") + b1_unpol_name(b1_unpol) +
          " is read ONLY by b1_model = li6-convolution (miller is a ratio "
          "model with no F1 of its own, and cdks carries the digitized "
          "column), got b1_model = " + b1_model_name(b1_model) +
          " -- the flag would not reach the rate but WOULD reach "
          "meta[\"b1_unpol\"], so it is refused rather than recorded as a "
          "variation that did not run");
    }
    if (!b1_unpol_sf) {
      throw std::runtime_error(
          std::string("PipelineConfig: b1_unpol = ") + b1_unpol_name(b1_unpol) +
          " but b1_unpol_sf is empty, and it is NEVER silently replaced by "
          "ToyF2" +
          (b1_unpol == B1UnpolSource::Mstw
               ? " -- MstwSF reads PYTHIA 8's own "
                 "pdfdata/mstw2008lo.00.dat and lives in the OPTIONAL PYTHIA "
                 "tier (include/lipolgen/mstw_sf.hpp, src/pythia/"
                 "mstw_sf.cpp), which the core library deliberately does not "
                 "link.  Configure with -DLIPOLGEN_WITH_PYTHIA=ON (and a "
                 "pdfdata grid on disk; $LIPOLGEN_PYTHIA8_PDFDATA overrides "
                 "the compiled-in path), then attach it: Python "
                 "lipolgen.make_config(b1_unpol='mstw') / "
                 "_lipolgen.set_b1_unpol(cfg, B1UnpolSource.Mstw), C++ "
                 "cfg.b1_unpol_sf = std::make_shared<const MstwSF>()"
               : b1_unpol == B1UnpolSource::Ct18Nlo
                     ? " -- LhapdfSF lives in the OPTIONAL LHAPDF tier "
                       "(src/lhapdf/lhapdf_sf.cpp), which the core library "
                       "deliberately does not link.  Configure with "
                       "-DLIPOLGEN_WITH_LHAPDF=ON (and the CT18NLO set in "
                       "LHAPDF's store), then attach it: Python "
                       "lipolgen.make_config(b1_unpol='ct18nlo'), C++ "
                       "cfg.b1_unpol_sf = std::make_shared<const LhapdfSF>("
                       "\"CT18NLO\", 0)"
                     : " -- Custom means the object YOU attach; put it in "
                       "b1_unpol_sf or leave b1_unpol at Toy"));
    }
  }
  // --r-source: the ONE R hook `default_inclusive_kernel` threads into the
  // `Li6Convolution` numerator (`Li6ConvolutionOptions::r_func`) AND the
  // kernel that carries the tensor weight's denominator
  // (`InclusiveKernel::Options::r_func`).  The `b1_unpol` provenance rule
  // verbatim, and only its clause (1) applies: `r_sigma_lt` and `r1998` both
  // live in the CORE (`sf.hpp`), so there is no optional tier to be missing
  // and no clause (2).  Only the `Li6Convolution` branch installs the hook --
  // `Miller` is a ratio model with no F1 of its own and `Cdks` carries the
  // digitized column, so neither has a numerator R to share -- while
  // `meta["r_source"]` is written unconditionally, so any other model would
  // claim an R that never touched the rate.  Nothing more is needed for the
  // channel, the isotope or the caller-supplied kernel: `b1_model =
  // li6-convolution` is ALREADY refused off the inclusive channel, off 6Li
  // and beside a `kernel` a few clauses above, so this one condition is the
  // whole reach rule.
  if (r_source != RSource::Unset && b1_model != B1Model::Li6Convolution) {
    throw std::runtime_error(
        std::string("PipelineConfig: r_source = ") + r_source_name(r_source) +
        " is read ONLY by b1_model = li6-convolution -- it is the ONE R the "
        "alpha-d convolution's own F1 and the kernel's F1/F_L/D(y) are both "
        "built from, and miller (a ratio model with no F1 of its own) and "
        "cdks (a digitized column) have no numerator R to share -- got "
        "b1_model = " + b1_model_name(b1_model) +
        ".  The flag would not reach the rate but WOULD reach "
        "meta[\"r_source\"], so it is refused rather than recorded as a "
        "variation that did not run");
  }
  // --unpol-sf / --pol-sf: the structure-function backend of EVERY kernel
  // this config builds (`UnpolSfSource`, `PolSfSource`).  Five clauses, all
  // of them the rules already enforced above, and the last two new because
  // these two selectors are the first ones that reach the tagged channels.
  //
  // (1) A NAMED BACKEND WITH AN EMPTY SLOT IS AN ERROR, NEVER A FALLBACK --
  //     the `b1_unpol` rule verbatim.  The core library links neither the
  //     PYTHIA nor the LHAPDF tier, so it cannot build `MstwSF` /
  //     `LhapdfSF` / `LhapdfG1`; falling back to the toy would put the toy's
  //     numbers under a label that says otherwise, and the toy is 20 % away
  //     on the 6Li rate and has the WRONG SIGN on g1n over 0.25 < x < 0.6.
  // (2) `Toy` WITH AN OBJECT is the same contradiction the other way: the
  //     Python setters name an attached object `Custom`, so reaching this
  //     state means the provenance and the realisation drifted apart.
  // (3) A CALLER-SUPPLIED KERNEL WINS over both selectors on the inclusive
  //     branch and is not even read on the tagged one, while `meta` would
  //     print both -- the `b1_model` argument verbatim.
  // (4) THE b1 COLLISION.  See `B1UnpolSource`: with a non-toy `unpol_sf`
  //     the `Li6Convolution` branch's null `b1_unpol` no longer means
  //     `ToyF2`, so `meta["b1_unpol"] = "toy"` would name something else.
  // (5) A DIRECTLY SET `struck.f2_source` / `struck.g1_model` with the
  //     matching selector still at `Toy` -- same rule, tagged branch.
  //
  // (6) A CALLER-SUPPLIED KERNEL ON A TAGGED CHANNEL, selectors or no
  //     selectors -- the `cluster_wave` rule (a knob that is never read on
  //     this channel is refused, not accepted in silence) applied to the
  //     escape hatch itself.  See `PipelineConfig::kernel`.
  //
  // What is deliberately NOT refused: either selector on a TAGGED or the
  // COHERENT channel.  Unlike `b1_model`, `unpol_sf` reaches the rate on
  // every channel -- that is the whole point of it -- so the "inclusive
  // only" guard above has no analogue here.  `pol_sf` does NOT reach the
  // coherent rate: nothing on that channel evaluates g1, which is a physics
  // property of it and not a wiring gap (`pol_sf_is_read`).  It is not
  // refused there either, because refusing one half of a flag PAIR on one
  // channel of a scan costs more than it buys; it is LABELLED instead, so
  // `meta["pol_sf"]` records `pol_sf_unread_label` rather than a backend
  // name and `meta["pol_sf_reach"]` carries the reason.
  if (unpol_sf == UnpolSfSource::Toy) {
    if (unpol_sf_obj) {
      throw std::runtime_error(
          "PipelineConfig: unpol_sf = toy with an object in unpol_sf_obj -- "
          "the toy setting IS the kernel's own ToyF2, so the attached "
          "backend would be silently dropped.  Use UnpolSfSource::Custom to "
          "run on your own UnpolSF (Python: config.unpol_sf_obj = obj, which "
          "sets Custom for you), or clear the slot");
    }
  } else if (!unpol_sf_obj) {
    throw std::runtime_error(
        std::string("PipelineConfig: unpol_sf = ") + unpol_sf_name(unpol_sf) +
        " but unpol_sf_obj is empty, and it is NEVER silently replaced by "
        "ToyF2" +
        (unpol_sf == UnpolSfSource::Mstw
             ? " -- MstwSF reads PYTHIA 8's own pdfdata/mstw2008lo.00.dat and "
               "lives in the OPTIONAL PYTHIA tier (include/lipolgen/"
               "mstw_sf.hpp, src/pythia/mstw_sf.cpp), which the core library "
               "deliberately does not link.  Configure with "
               "-DLIPOLGEN_WITH_PYTHIA=ON (and a pdfdata grid on disk; "
               "$LIPOLGEN_PYTHIA8_PDFDATA overrides the compiled-in path), "
               "then attach it: Python lipolgen.make_config(unpol_sf='mstw') "
               "/ _lipolgen.set_unpol_sf(cfg, UnpolSfSource.Mstw), C++ "
               "cfg.unpol_sf_obj = std::make_shared<const MstwSF>()"
             : unpol_sf == UnpolSfSource::Ct18Nlo
                   ? " -- LhapdfSF lives in the OPTIONAL LHAPDF tier "
                     "(src/lhapdf/lhapdf_sf.cpp), which the core library "
                     "deliberately does not link.  Configure with "
                     "-DLIPOLGEN_WITH_LHAPDF=ON (and the CT18NLO set in "
                     "LHAPDF's store), then attach it: Python "
                     "lipolgen.make_config(unpol_sf='ct18nlo'), C++ "
                     "cfg.unpol_sf_obj = std::make_shared<const LhapdfSF>("
                     "\"CT18NLO\", 0)"
                   : " -- Custom means the object YOU attach; put it in "
                     "unpol_sf_obj or leave unpol_sf at Toy"));
  }
  if (pol_sf == PolSfSource::Toy) {
    if (pol_sf_obj) {
      throw std::runtime_error(
          "PipelineConfig: pol_sf = toy with an object in pol_sf_obj -- the "
          "toy setting IS the kernel's own ToyG1, built on the kernel's own "
          "UnpolSF and r_func, so the attached backend would be silently "
          "dropped.  Use PolSfSource::Custom to run on your own PolSF "
          "(Python: config.pol_sf_obj = obj, which sets Custom for you), or "
          "clear the slot");
    }
  } else if (!pol_sf_obj) {
    throw std::runtime_error(
        std::string("PipelineConfig: pol_sf = ") + pol_sf_name(pol_sf) +
        " but pol_sf_obj is empty, and it is NEVER silently replaced by "
        "ToyG1" +
        (pol_sf == PolSfSource::NnpdfPol
             ? " -- LhapdfG1 lives in the OPTIONAL LHAPDF tier "
               "(src/lhapdf/lhapdf_sf.cpp), which the core library "
               "deliberately does not link.  Configure with "
               "-DLIPOLGEN_WITH_LHAPDF=ON (and the NNPDFpol11_100 set in "
               "LHAPDF's store -- it is the only POLARISED set installed "
               "here), then attach it: Python "
               "lipolgen.make_config(pol_sf='nnpdfpol'), C++ cfg.pol_sf_obj "
               "= std::make_shared<const LhapdfG1>(\"NNPDFpol11_100\", 0)"
             : " -- Custom means the object YOU attach; put it in pol_sf_obj "
               "or leave pol_sf at Toy"));
  }
  // (5) ... and the same rule for the tagged branch's own slots.  A C++
  //     caller may put an object straight into `struck.f2_source` /
  //     `struck.g1_model` (the `breakup.triton_sf` arrangement: `Pipeline`
  //     fills them only when they are empty), but then `meta["unpol_sf"]`
  //     would still write "toy" for a run made on something else.  Naming
  //     the selector `Custom` alongside costs one line and keeps the file
  //     honest.
  if (struck.f2_source && unpol_sf == UnpolSfSource::Toy) {
    throw std::runtime_error(
        "PipelineConfig: struck.f2_source is set but unpol_sf = toy -- the "
        "tagged kernel would run on your object while meta[\"unpol_sf\"] "
        "recorded \"toy\".  Set unpol_sf = Custom as well (Python: "
        "config.unpol_sf_obj = obj, which sets Custom for you, and Pipeline "
        "then fills struck.f2_source from it), or clear the slot");
  }
  if (struck.g1_model && pol_sf == PolSfSource::Toy) {
    throw std::runtime_error(
        "PipelineConfig: struck.g1_model is set but pol_sf = toy -- the "
        "tagged kernel would run on your object while meta[\"pol_sf\"] "
        "recorded \"toy\".  Set pol_sf = Custom as well (Python: "
        "config.pol_sf_obj = obj), or clear the slot");
  }
  if (kernel && (unpol_sf != UnpolSfSource::Toy ||
                 pol_sf != PolSfSource::Toy)) {
    throw std::runtime_error(
        std::string("PipelineConfig: unpol_sf = ") + unpol_sf_name(unpol_sf) +
        " / pol_sf = " + pol_sf_name(pol_sf) +
        " with a caller-supplied kernel: on the inclusive channel the kernel "
        "wins and the selectors would be silently ignored, and on a tagged "
        "channel the kernel is not read at all -- either way meta would "
        "record a backend that did not run, so the two are refused together "
        "(drop one; a hand-built kernel already carries its own f2_source "
        "and g1_model)");
  }
  if (kernel && is_tagged(channel)) {
    throw std::runtime_error(
        std::string("PipelineConfig: a caller-supplied kernel is not read on "
                    "channel ") + pipeline_channel_name(channel) +
        " -- the tagged channels draw from the STRUCK-CLUSTER sampler "
        "(make_struck_cluster_source), and `kernel` reaches the ion-level "
        "sampler of the inclusive and coherent channels alone.  Accepting it "
        "would run ToyF2/ToyG1 while meta[\"unpol_sf\"], meta[\"pol_sf\"] "
        "and meta[\"b1_model\"] all recorded \"caller-supplied kernel\": "
        "measured, a tagged-6Li-alpha config carrying an InclusiveKernel "
        "built on CT18NLO ran at F2A(0.3, 10) = 0.3691149345, the TOY value, "
        "and CT18NLO's 0.4639703442 appeared nowhere.  Inject structure "
        "functions here through unpol_sf / pol_sf, which fill "
        "struck.f2_source / struck.g1_model, or put your own objects "
        "straight into those two slots and name the selectors Custom");
  }
  if (b1_model == B1Model::Li6Convolution &&
      unpol_sf != UnpolSfSource::Toy && b1_unpol == B1UnpolSource::Toy) {
    throw std::runtime_error(
        std::string("PipelineConfig: unpol_sf = ") + unpol_sf_name(unpol_sf) +
        " with b1_unpol = toy on b1_model = li6-convolution: `toy` there "
        "means \"the kernel's own UnpolSF, shared as one object\", and that "
        "object is no longer ToyF2 -- so meta[\"b1_unpol\"] would record "
        "\"toy\" for a b1 folded against " + unpol_sf_name(unpol_sf) +
        ".  Name the b1 backend explicitly instead: " +
        (unpol_sf == UnpolSfSource::Custom
             ? std::string("attach the same object (Python: "
                           "config.b1_unpol_sf = config.unpol_sf_obj)")
             : std::string("--b1-unpol ") + unpol_sf_name(unpol_sf) +
                   " to fold against the same one") +
        ", or --b1-unpol mstw for the configuration the A = 2 gate passes "
        "on (G3b 0.843243)");
  }
  scenario.validate();
}

std::shared_ptr<const InclusiveKernel> default_inclusive_kernel(const Ion& ion) {
  return default_inclusive_kernel(ion, B1Model::Miller, 1.0, 1.0, nullptr,
                                  nullptr, nullptr, RSource::Unset);
}

std::shared_ptr<const InclusiveKernel> default_inclusive_kernel(
    const Ion& ion, B1Model model, double band_scale, double w_alpha_d,
    std::shared_ptr<const UnpolSF> b1_unpol,
    std::shared_ptr<const UnpolSF> unpol_sf,
    std::shared_ptr<const PolSF> pol_sf, RSource r_source) {
  InclusiveKernel::Options opt;
  // ONE UnpolSF, shared: the core must not link LHAPDF or PYTHIA (sf.hpp),
  // so an `LhapdfSF` / `MstwSF` only ever enters through the `unpol_sf`
  // argument (filled one layer up by `_lipolgen.set_unpol_sf`), a
  // caller-supplied `cfg.kernel`, or a validation script.  Naming the object
  // here rather than letting `InclusiveKernel` default it is what lets the
  // `Li6Convolution` branch hand the SAME object to
  // `Li6ConvolutionOptions::unpol`, so "the kernel's UnpolSF" is literally
  // one object and not a second one that merely looks like it.  ToyF2 is
  // stateless closed form, so the default is bit for bit the old null
  // default (validation/reference/b1_default_li6.json, T9).
  // `--unpol-sf` / `--pol-sf` (`UnpolSfSource`, `PolSfSource`).  NULL -- the
  // default, and every caller that predates them -- keeps the shared `ToyF2`
  // named on the line below and leaves `g1_model` UNSET, so
  // `InclusiveKernel` takes its own null branch and builds `ToyG1` on that
  // same object with this kernel's own `r_func`, exactly as before: bit for
  // bit BY CONSTRUCTION, which is what `validation/reference/*.json` at
  // rtol 1e-12 relies on.  Non-null replaces the object and NOTHING else --
  // in particular the `Li6Convolution` branch below still shares whatever
  // `f2` ends up being, which is the collision `PipelineConfig::validate()`
  // refuses rather than mislabel.
  const std::shared_ptr<const UnpolSF> f2 =
      unpol_sf ? std::move(unpol_sf)
               : std::static_pointer_cast<const UnpolSF>(
                     std::make_shared<const ToyF2>());
  opt.f2_source = f2;
  if (pol_sf) opt.g1_model = std::move(pol_sf);
  // ONE definition of the spin test (`kernel_fills_rank2`), read here and by
  // `inclusive_rank2_is_empty` below: at spin 3/2 NOTHING is filled and the
  // whole tensor and cos 2phi sector of the run is identically zero, which is
  // what the run surface now says out loud instead of returning silent zeros.
  if (kernel_fills_rank2(ion)) {
    if (model == B1Model::Miller) {
      // UNCHANGED DEFAULT PATH, including the function-local `static`:
      // `TensorSF::b1_func()` returns a closure capturing raw `this` and
      // "must not outlive the backend" (sf.hpp), and this is how that
      // lifetime has always been met here.  `band_scale` is deliberately NOT
      // applied -- Miller's numbers are the published ones and stay bit for
      // bit; the band belongs to the two opt-in backends.
      static const auto li6_b1 =
          std::make_shared<Li6B1>(std::make_shared<MillerB1>());
      opt.b1_func = li6_b1->b1_func();
    } else if (model == B1Model::Cdks) {
      // NOT a `static`: that would freeze the first band scale for the whole
      // process, so a `--b1-band-scale 0` run after a `1` run would silently
      // get the wrong object.  Capturing the shared_ptr BY VALUE in the
      // closure keeps the backend alive exactly as long as the closure --
      // which is what `b1_func()`'s raw-`this` capture cannot do.
      const auto b = std::make_shared<const Li6B1>(std::make_shared<CdksB1>());
      opt.b1_func = [b, band_scale](double x, double q2, double f1) {
        return band_scale * b->b1(x, q2, f1);
      };
    } else {
      Li6ConvolutionOptions o;
      o.w_alpha_d_dwave = w_alpha_d;
      // `--b1-unpol`.  Null (the default, and every caller that predates the
      // flag) is the kernel's own `f2` object, handed over as the SAME
      // object, so this line is bit for bit what it was.  Non-null is the
      // selected backend -- MSTW2008 LO is what CDKS computed their b1_d
      // with, and it is worth up to a factor 1.85 on `Li6ConvolutionB1::b1`
      // (pipeline.hpp, `B1UnpolSource`).  AT `--unpol-sf toy` `opt.f2_source`
      // above is ToyF2 either way, so this flag then moves b1 and only b1:
      // the SPIN-BLIND cell cross section
      // (`InclusiveSampler::cell_xsec_pb`) is bit-identical across settings,
      // and the tensor-weighted per-category cross sections move, which is
      // the point.  That asymmetry is deliberate and documented at
      // `B1UnpolSource`.  Off that default the shared object is no longer
      // `ToyF2`, and `b1_unpol = Toy` would then record "toy" for something
      // else -- which `PipelineConfig::validate()` refuses outright rather
      // than resolve silently.
      o.unpol = b1_unpol ? std::move(b1_unpol) : f2;
      // `--r-source`: ONE R OBJECT INTO BOTH HALVES OF THE RATIO.  The
      // shipped tensor observable is K/D_phi (`azz`, asymmetries.hpp); its
      // numerator's R is `Li6ConvolutionOptions::r_func` and its
      // denominator's is this kernel's own `Options::r_func`, and until this
      // branch existed the only way to move one was to move it ALONE.  Here
      // one `RFunc` is built and COPIED into both, so whichever R is named
      // the two halves are the same choice -- registry option (iii) of
      // STATUS.md row 3, priced in
      // docs/open_items/run_2026-09-06/phase_A_numbers.md sec. A1.
      //
      // `Unset` DOES NOT ENTER HERE AT ALL.  That is deliberate: leaving both
      // hooks null is bit for bit the pre-flag tree by construction, not by
      // an argument that an explicit `r_sigma_lt` closure resolves to the
      // same double as `resolve_r`'s null branch.  (It does -- `SigmaLt` is
      // measured bit-identical to `Unset` on every tensor observable at the
      // standard points, sec. A1 -- but the DEFAULT does not rely on it.)
      if (r_source != RSource::Unset) {
        const RFunc r =
            (r_source == RSource::R1998)
                ? RFunc([](double x, double q2) { return r1998(x, q2); })
                : RFunc([](double x, double q2) { return r_sigma_lt(x, q2); });
        o.r_func = r;
        opt.r_func = r;
      }
      const auto b = std::make_shared<const Li6ConvolutionB1>(std::move(o));
      opt.b1_func = [b, band_scale](double x, double q2, double f1) {
        return band_scale * b->b1(x, q2, f1);
      };
    }
    // The Delta slot is the SAME for every `B1Model` -- it is not part of
    // this choice and must not be dropped on the new branches.
    opt.delta_func = [](double x, double q2, double f1) {
      return toy_delta_gluon(x, q2, f1, 1e-2);
    };
  }
  return std::make_shared<InclusiveKernel>(ion, opt);
}

bool inclusive_rank2_is_empty(const PipelineConfig& cfg) {
  if (cfg.channel != PipelineChannel::Inclusive) return false;
  if (cfg.kernel) return false;
  const Ion& ion = ion_by_name(cfg.isotope);
  return !kernel_fills_rank2(ion) && kernel_has_rank2_sector_unfilled(ion);
}

std::string rank2_input_report(const PipelineConfig& cfg, const RunPlan& plan) {
  if (cfg.channel != PipelineChannel::Inclusive) {
    return std::string("not read on channel ") +
           pipeline_channel_name(cfg.channel) +
           " -- the inclusive rank-2 slots reach the rate on the inclusive "
           "channel only (a tagged channel's alignment is in the event "
           "weight through TaggedModel, and the coherent channel's tensor "
           "signal is the recoil azimuth's 1 + c2 cos 2(phi_t - phi_S))";
  }
  if (cfg.kernel) return "caller-supplied kernel";
  const Ion& ion = ion_by_name(cfg.isotope);
  if (kernel_fills_rank2(ion)) {
    return std::string("b1_model = ") + b1_model_name(cfg.b1_model) +
           " (spin 1: b1_func filled, b2 = 2x*b1 by default, delta_func = "
           "toy_delta_gluon at the 1e-2 discovery scale)";
  }
  if (!kernel_has_rank2_sector_unfilled(ion)) {
    return std::string("no rank-2 sector at spin ") + fmt_g(ion.spin) +
           " (InclusiveKernel::tables carries rank-2 slots for spin 1 and "
           "spin 3/2 only, and a spin-1/2 or spin-0 target has no rank-2 "
           "structure function to carry)";
  }
  std::string s =
      "EMPTY -- " + cfg.isotope +
      " is spin 3/2 and default_inclusive_kernel fills a rank-2 slot for "
      "spin 1 ONLY, so b1_32 = b2_32 = delta_32 = 0 and the tensor term of "
      "the phi-averaged rate, the cos 2phi (gluon transversity) amplitude "
      "and therefore A_zz of this run are IDENTICALLY ZERO, not small: two "
      "categories that differ only in their alignment come back as the same "
      "double. There is no published b1 for 7Li and no backend in this build "
      "fills it";
  if (plan.pzz_true() != 0.0) {
    s += "; THIS RUN'S PLAN CARRIES A RANK-2 FILL (T = " +
         fmt_g(plan.pzz_true()) +
         ") AND NONE OF IT REACHES THE RATE";
  }
  s += " (docs/OPEN_ITEMS_SOLUTIONS.md open item 15; the physics, the "
       "measured leading estimate and why it is not shipped are in "
       "docs/open_items/run_2026-09-03/phase_D_li7_rank2.md)";
  return s;
}

bool plan_has_beam_helicity(const RunPlan& plan) {
  // The EXACT condition, and it is one product: `InclusiveKernel::amplitudes`
  // adds `lam_e * pe * (m / J) * cos(theta_S) * A_par(g1, g2, F1)` and there
  // is no other reader of a `PolSF` anywhere in a run.  So a fill needs a
  // beam helicity AND a non-zero population at some m != 0 before any
  // polarised structure function is evaluated at all.
  for (const SpinCategory& c : plan.categories()) {
    if (c.lam_e == 0 || c.pe == 0.0) continue;
    const std::size_t n = c.populations.size();
    for (std::size_t i = 0; i < n; ++i) {
      // m runs +J ... -J, so the m = 0 entry (odd n only) is the middle one.
      const bool m_is_zero = (2 * i + 1 == n);
      if (!m_is_zero && c.populations[i] != 0.0) return true;
    }
  }
  return false;
}

bool pol_sf_is_read(const PipelineConfig& cfg, const RunPlan& plan) {
  // TWO AXES, ONE DEFINITION, read by `meta["pol_sf"]`, by
  // `pol_sf_reach_report`, by the `pol_sf` row of `Pipeline::knob_provenance`
  // and through it by the CLI banner.
  //
  // THE CHANNEL.  Every channel but `CoherentLi6` evaluates g1: the inclusive
  // one through `InclusiveKernel::Options::g1_model`, the three tagged ones
  // through `StruckClusterOptions::g1_model`.
  //
  // THE RUN PLAN.  `plan_has_beam_helicity` is the whole of it -- see there.
  // This clause is the one that was missing until 2026-09-05, and it is the
  // one that bites at the CLI's own default plan.
  if (cfg.channel == PipelineChannel::CoherentLi6) return false;
  return plan_has_beam_helicity(plan);
}

std::string pol_sf_unread_label(const PipelineConfig& cfg,
                                const RunPlan& plan) {
  // The CHANNEL label is checked first and its text is UNCHANGED, because it
  // is quoted verbatim in docs/PHYSICS_CHANNELS.md and pinned by
  // python/tests/test_sf_backend.py: on the coherent channel no fill of any
  // kind would make g1 run, so that is the stronger statement of the two.
  if (cfg.channel == PipelineChannel::CoherentLi6) {
    return std::string("not read on channel ") +
           pipeline_channel_name(cfg.channel);
  }
  (void)plan;
  // SHORT on both axes, deliberately: this is the string a `meta` key carries
  // IN PLACE OF a backend name, and `KnobProvenance::label` derives the same
  // string from the report's own opening clause.  Two spellings of one label
  // is how the file and the table start disagreeing.  The explanation is
  // `pol_sf_reach_report`'s, beside it in `meta["pol_sf_reach"]`.
  return "not read under this run's fill";
}

std::string pol_sf_reach_report(const PipelineConfig& cfg,
                                const RunPlan& plan) {
  if (cfg.channel != PipelineChannel::CoherentLi6 &&
      !plan_has_beam_helicity(plan)) {
    return pol_sf_unread_label(cfg, plan) +
           ": no category carries lam_e * P_e != 0, and lam_e * P_e is the "
           "only thing g1 is multiplied by (InclusiveKernel::amplitudes adds "
           "helicity * (m/J) * cos(theta_S) * A_par and nothing else reads a "
           "PolSF).  The three TENSOR plans "
           "build every category at lam_e = 0, pe = 0 "
           "(bookkeeping.cpp: tensor_thirds_plan, transverse_tensor_plan, "
           "tensor_flip_plan) and helicity_flip_plan at --pe 0 is the same; "
           "tensor-thirds is the CLI's DEFAULT plan.  MEASURED 2026-09-05, "
           "600 events seed 7, sha256 over all 47 columns plus sigma_pb and "
           "sigma_per_category_pb: --pol-sf nnpdfpol is bit-identical to toy "
           "on inclusive-6Li, inclusive-d, tagged-6Li-alpha and tagged-d-p "
           "under tensor-thirds, and on inclusive-6Li under "
           "transverse-tensor, tensor-flip and helicity-flip at --pe 0.  The "
           "kernel still CARRIES the selected g1_model, so "
           "dis_sampler().kernel().tables(x, q2).g1 does move -- the RUN "
           "never asks for it.  --unpol-sf, by contrast, reaches this run "
           "(it is the unpolarised rate).  Run --plan helicity-flip at "
           "--pe != 0 to make this selector matter";
  }
  if (!pol_sf_is_read(cfg, plan)) {
    return pol_sf_unread_label(cfg, plan) +
           " -- the coherent rate is SPIN-INDEPENDENT (sigma_coh is the sum "
           "over accepted cells of the UNPOLARIZED sigma_cell(x, Q2) times "
           "f_coh(x)), and this channel's tensor signal is the recoil "
           "azimuth's 1 + c2 cos 2(phi_t - phi_S), which CoherentSampler "
           "owns.  Nothing in the run reads g1: not the rate, not a column, "
           "and not RcModel, whose applies() is false here.  Measured "
           "2026-09-05, --events 400 --seed 11: sigma_pb, "
           "sigma_per_category_pb and all 47 generated array columns are "
           "bit-identical between --pol-sf toy and --pol-sf nnpdfpol.  The "
           "kernel still carries the selected g1_model, so "
           "dis_sampler().kernel().tables(x, q2).g1 does move -- the RUN "
           "never asks for it.  --unpol-sf, by contrast, DOES reach this "
           "channel (x0.705841 on the shipped 6Li ct18nlo run), through the "
           "very cell cross sections above";
  }
  if (cfg.kernel) return "caller-supplied kernel";
  return std::string("pol_sf = ") + pol_sf_name(cfg.pol_sf) +
         (is_tagged(cfg.channel)
              ? " -> StruckClusterOptions::g1_model (the tagged "
                "struck-cluster kernel)"
              : " -> InclusiveKernel::Options::g1_model") +
         "; g2 follows it through Wandzura-Wilczek, and at pol_sf = toy the "
         "kernel's own ToyG1 is built on the kernel's own UnpolSF, so "
         "--unpol-sf moves g1 too";
}

// ============================================================ event helpers

Vec4 momentum_residual(const Event& ev) {
  const Vec4 in = momentum_scale(ev);
  Vec4 out;
  for (const Particle& p : ev.particles) {
    if (p.status == Status::Final) out = out + p.p;
  }
  return in - out;
}

Vec4 momentum_scale(const Event& ev) {
  const Particle* be = ev.find(Role::BeamElectron);
  if (!be) throw std::runtime_error("momentum_scale: no beam electron");
  if (ev.channel == Channel::Inclusive) {
    const Particle* n = ev.find(Role::StruckNucleon);
    if (!n) throw std::runtime_error("momentum_scale: no struck nucleon");
    return be->p + n->p;
  }
  const Particle* bi = ev.find(Role::BeamIon);
  if (!bi) throw std::runtime_error("momentum_scale: no beam ion");
  return be->p + bi->p;
}

double charge_residual(const Event& ev) {
  const Particle* be = ev.find(Role::BeamElectron);
  if (!be) throw std::runtime_error("charge_residual: no beam electron");
  double in = be->charge;
  if (ev.channel == Channel::Inclusive) {
    const Particle* n = ev.find(Role::StruckNucleon);
    if (!n) throw std::runtime_error("charge_residual: no struck nucleon");
    in += n->charge;
  } else {
    const Particle* bi = ev.find(Role::BeamIon);
    if (!bi) throw std::runtime_error("charge_residual: no beam ion");
    in += bi->charge;
  }
  return ev.total_charge_final() - in;
}

int route_of(const Event& ev, const Optics& optics,
             const std::string& pot_config) {
  const Particle* f = ev.find(Role::Spectator);
  if (!f) f = ev.find(Role::IntactRecoil);
  if (!f) return kRouteLost;
  const double pt = f->p.pt();
  const double theta = std::atan2(pt, f->p.pz);
  if (f->charge == 0.0) return route_neutral(theta);
  const Particle* bi = ev.find(Role::BeamIon);
  if (!bi || bi->charge <= 0.0) throw std::runtime_error("route_of: no beam ion");
  const double r = (f->p.p() / f->charge) / (bi->p.pz / bi->charge);
  return route_charged(r, theta, pt, optics, std::atan2(f->p.py, f->p.px),
                       kNaN, pot_config);
}

bool rp_tagged(const Event& ev, const Optics& optics,
               const std::string& pot_config) {
  return rp_accepted(route_of(ev, optics, pot_config));
}

bool struck_cluster_of(const Event& ev, const TaggedChannel& channel,
                       StruckCluster& out) {
  const Particle* x = ev.find(Role::StruckCluster);
  const Particle* s = ev.find(Role::Spectator);
  const Particle* i = ev.find(Role::BeamIon);
  if (!x || !s || !i) return false;
  const ClusterChannel& base = channel.base;
  out.p_ion = i->p;
  out.p_spectator = s->p;
  out.p = x->p;
  out.m2 = x->p.m2();
  out.m_free = base.m_partner();
  out.virtuality = out.m2 - out.m_free * out.m_free;
  out.alpha_s = ev.kin.alpha_s;
  out.alpha_x = base.beam_A - out.alpha_s;
  out.pt_s = ev.kin.pt_s;
  out.p_per_nucleon_eff = x->p.pz / base.partner_A();
  return true;
}

// ============================================================ the pipeline

Pipeline::Pipeline(PipelineConfig config, RunPlan plan)
    : cfg_(std::move(config)), plan_(std::move(plan)) {
  cfg_.validate();
  if (plan_.categories().empty()) {
    throw std::runtime_error("Pipeline: the run plan has no category");
  }
  beams_ = default_configs(cfg_.isotope)[static_cast<std::size_t>(cfg_.beam_config)];
  const double p_u = beams_.ion_momentum_per_nucleon;
  optics_ = cfg_.optics_choice == OpticsChoice::Custom
                ? cfg_.optics
                : optics_for(cfg_.optics_choice, cfg_.isotope, p_u, cfg_.n_sigma);
  pot_config_ = cfg_.pot_config.empty() ? yr_config_key(cfg_.isotope, p_u)
                                        : cfg_.pot_config;

  const std::size_t nc = plan_.categories().size();
  sigma_.assign(nc, 0.0);

  if (is_tagged(cfg_.channel)) {
    channel_.reset(new TaggedChannel(
        cfg_.channel == PipelineChannel::TaggedLi6Alpha
            ? li6_alpha_channel(cfg_.cluster_beta, cfg_.p_d, cfg_.cluster_wave,
                                cfg_.cluster_vmc_mc_sigma)
            : cfg_.channel == PipelineChannel::TaggedLi7Alpha
                  ? li7_alpha_channel(cfg_.cluster_beta, cfg_.cluster_wave,
                                      cfg_.cluster_vmc_mc_sigma)
                  : deuteron_channel(cfg_.cluster_beta, P_D_DEUTERON,
                                     cfg_.cluster_wave)));
    model_ = std::make_shared<TaggedModel>(*channel_);
    StruckClusterOptions sopt = cfg_.struck;
    sopt.scenario = cfg_.scenario;
    sopt.grid = cfg_.grid;
    // --unpol-sf / --pol-sf on a TAGGED run.  This is D2's injection point:
    // before it existed the struck-cluster kernel had no structure-function
    // slot at all, and `cfg_.kernel` -- read on the non-tagged branch alone,
    // a few dozen lines below -- never reached here.  Both are NULL at the
    // shipped default, so a default tagged run is bit for bit.  A C++ caller
    // who already put an object in `cfg_.struck` keeps it, the
    // `breakup.triton_sf` rule.
    if (!sopt.f2_source) sopt.f2_source = cfg_.unpol_sf_obj;
    if (!sopt.g1_model) sopt.g1_model = cfg_.pol_sf_obj;
    dis_source_ = make_struck_cluster_source(*channel_, beams_, sopt);
    dis_sampler_ = dis_source_->sampler_ptr();
    tsampler_.reset(new TaggedSampler(*model_, p_u, dis_source_.get(), optics_,
                                      pot_config_));

    // FSI as a per-event weight (fsi.hpp).  Built from the MODEL so the
    // normalization grid -- and with it P_L -- is bit-for-bit the sampler's
    // own; every draw then carries `TaggedEvent::weight`, which `make_tagged`
    // multiplies into `Event::weight`.  Off = today's PWIA, bit for bit.
    if (cfg_.fsi != PipelineFsi::Off) {
      GlauberFsiOptions fopt;
      fopt.variant = cfg_.fsi == PipelineFsi::GlauberNucleon
                         ? FsiVariant::GlauberNucleon
                         : FsiVariant::GlauberCluster;
      fopt.sigma_xn_mb = cfg_.fsi_sigma_mb;
      fsi_ = std::make_shared<GlauberFsiWeight>(*model_, fopt);
      tsampler_->set_fsi(fsi_);
    }

    // NO WARM-UP LOOP HERE, and the comment that asked for one is gone with
    // it (2026-09-16).  It read "`TaggedModel`'s amplitude, density and
    // cell-CDF maps are `mutable` and NOT mutex-protected, so building them
    // all here is what makes threaded generation safe" and then walked every
    // (M, m_S) through `n_of_kc`, `population_integrated` and
    // `sample_kc_one`.  Since C3 those maps are built in FULL by
    // `TaggedModel`'s constructor and never written again (tagged.hpp:
    // "the object is immutable after construction ... no lock is needed on
    // any path"), so every one of those calls was a pure lookup on an
    // already-built table: the loop cost a constructor pass and changed no
    // number.  The sampler's OWN per-state cache is still warmed below, per
    // category -- that one is mutex-guarded and warming it is a real
    // anti-contention measure.
    fills_.reserve(nc);
    rate_cdf_.reserve(nc);
    for (std::size_t k = 0; k < nc; ++k) {
      const SpinCategory& c = plan_.categories()[k];
      if (std::fabs(c.j - channel_->j_ion) > 1e-9) {
        throw std::runtime_error("Pipeline: run-plan spin " +
                                 std::to_string(c.j) + " != channel ion spin");
      }
      // A TILTED FILL AND A SPIN-DEPENDENT STRUCK-CLUSTER KERNEL DO NOT
      // COMPOSE ON A TAGGED CHANNEL.  `InclusiveKinematicsSource::plan_for`
      // builds the struck cluster's pure category at theta_S = phi_S = 0
      // whatever this fill's axis is (the Python's `_pure_category` does the
      // same), so the ion axis reaches the event through the SPECTATOR
      // rotation alone.  That is exact while the DIS side is spin-blind and
      // wrong the moment it is not: MEASURED 2026-09-16 (6Li config 1, seed
      // 1), tagged-6Li-alpha `sigma_per_category_pb` is bit-identical at
      // theta_S = 0 and theta_S = pi/2 under helicity-flip at P_e = 0.7,
      // where the inclusive channel's asymmetry correctly collapses from
      // -5.741406e-04 to 0; and `--plan transverse-tensor --inclusive-b1`
      // used to run and record `inclusive_b1` as READ while applying the
      // struck deuteron's b1 with P2(cos 0) = +1 on a fill at
      // P2(cos 90 deg) = -1/2.  Refused here, at configuration time, rather
      // than silently computed for the wrong orientation.
      //
      // ONLY THE TERMS THIS CHANNEL ACTUALLY READS COUNT, which is why the
      // two halves of the test are not the same shape.  P_e is read on all
      // three tagged channels -- MEASURED 2026-09-16, `sigma_per_category_pb`
      // moves between P_e = 0 and 0.7 on every one (asymmetry -1.72e-03
      // tagged-6Li-alpha, +6.41e-04 tagged-7Li-alpha, -5.73e-03 tagged-d-p).
      // `inclusive_b1` is read on tagged-6Li-alpha ALONE: the other two
      // struck clusters are spin 1/2 (dis_target = triton resp. free
      // neutron), `InclusiveKernel::tables` opens no rank-2 sector there and
      // the slot is never filled -- which is what `knob_provenance`'s own
      // `ib1_read` says, and MEASURED the same day, the knob moves
      // tagged-6Li-alpha (5.9184617e+05 -> 5.9156954e+05 / 5.9239941e+05 per
      // category) and is bit-for-bit inert on tagged-7Li-alpha and tagged-d-p.
      // Refusing a tilted fill for an INERT knob would be a false refusal
      // with a false reason attached, so it is not refused.
      //
      // NOTHING SPIN-BLIND MOVES: P_e = 0 without a read `inclusive_b1` is
      // every reference gate, every CLI default, and
      // `transverse_tensor_plan` / `tensor_flip_plan`, whose P_e is 0 by
      // construction -- the only two shipped plans that tilt at all.
      // Threading the axis into `plan_for`'s cache key and into
      // `KinematicsSource`'s signature -- so that a tilted tagged fill
      // COMPUTES instead of being refused -- is the alternative, and it is
      // the maintainer's call, not a registry row.
      const bool ib1_read =
          cfg_.struck.inclusive_b1 &&
          cfg_.channel == PipelineChannel::TaggedLi6Alpha;
      if (std::fabs(std::sin(c.theta_s)) > 1e-12 &&
          (c.pe != 0.0 || ib1_read)) {
        const std::string terms =
            (c.pe != 0.0 ? "P_e = " + fmt_g(c.pe) : std::string("")) +
            (c.pe != 0.0 && ib1_read ? ", " : "") +
            (ib1_read ? std::string("inclusive_b1 = true") : std::string(""));
        throw std::runtime_error(
            "Pipeline: run-plan category \"" + c.name + "\" sits at "
            "theta_S = " + fmt_g(c.theta_s) + " rad (sin theta_S = " +
            fmt_g(std::sin(c.theta_s)) + "), and channel " +
            pipeline_channel_name(cfg_.channel) + " reads a SPIN-DEPENDENT "
            "struck-cluster DIS term here (" + terms +
            ").  The struck cluster's |S_c m_S> is evaluated with its "
            "quantization axis along the BEAM whatever the fill's axis is, so "
            "the tagged rate would be computed for the wrong spin "
            "orientation -- tagged sigma_per_category_pb is bit-identical at "
            "theta_S = 0 and pi/2, where the inclusive channel's asymmetry "
            "correctly goes to zero.  Either run the tilted fill SPIN-BLIND "
            "(P_e = 0 and --inclusive-b1 off, which is what "
            "transverse-tensor and tensor-flip already are), or run the "
            "spin-dependent term on an UNTILTED fill (theta_S = 0), or use "
            "the inclusive channel, whose kernel does carry the axis.");
      }
      dis_source_->warm(c.lam_e, c.pe);
      IonFill f;
      f.name = c.name;
      f.j = c.j;
      f.populations = c.populations;
      f.lam_e = c.lam_e;
      f.pe = c.pe;
      f.theta_s = c.theta_s;
      f.phi_s = c.phi_s;
      fills_.push_back(f);
      rate_cdf_.push_back(tsampler_->rate_cdf(f));
      sigma_[k] = tsampler_->sigma_tot_pb(f);
    }

    // --- T1: the struck cluster is resolved into a nucleon + partners -----
    //
    // The species is read off the CHANNEL's own partner (Z, A), so nothing
    // here has to know which channel it is: (1, 2) -> deuteron, (1, 3) ->
    // triton, (Z, 1) -> a nucleon that only needs relabelling.
    cluster_species_ = cluster_species(channel_->base.partner_Z(),
                                       channel_->base.partner_A());
    tier_ = cfg_.tier;
    if (tier_ == Tier::T1) {
      BreakupOptions bo = cfg_.breakup;
      // One beta, one F2 backend.  The breakup's internal wave function and
      // its species draw MUST be the ones the rest of the run is made with,
      // or the T1 tier would describe a different nucleus from the T0 one
      // (docs/CONVENTIONS.md: no physics number is defined twice).
      bo.beta = cfg_.cluster_beta;
      // ... and ONE WAVE-FUNCTION FAMILY (C5.5b).  Same rule, same sentence
      // as the comment above: the breakup carries a deuteron of its own and
      // it must be the run's deuteron.  Before 2026-09-04 only `beta` was
      // forwarded, so `--cluster-wave vmc` left the T1 struck-nucleon spin
      // draw on the Hulthen deuteron while the rate moved to the AV18 one.
      bo.source = cfg_.cluster_wave;
      if (!bo.f2) bo.f2 = dis_sampler_->kernel().nuclear_f2().base();
      // --triton-sf ciofi-simula: the spectral function is built HERE, at
      // the run's own `cluster_beta` and the breakup's own `k_max`, so the
      // continuum pair's q-shape shares the one beta everything else in the
      // run is made with -- the same rule as the two lines above.  A C++
      // caller who already put an object in `breakup.triton_sf` keeps it.
      if (cfg_.triton_sf == TritonSfChoice::CiofiSimula && !bo.triton_sf) {
        CiofiSimulaOptions tso;
        tso.beta = cfg_.cluster_beta;
        tso.k_max = bo.k_max;
        bo.triton_sf = std::make_shared<const CiofiSimulaTriton>(tso);
      }
      breakup_.reset(new ClusterBreakup(bo));
    }
  } else {
    const Ion& ion = ion_by_name(cfg_.isotope);
    const auto kernel = cfg_.kernel
                            ? cfg_.kernel
                            : default_inclusive_kernel(
                                  ion, cfg_.b1_model, cfg_.b1_band_scale,
                                  cfg_.b1_alpha_d_dwave_weight,
                                  cfg_.b1_unpol_sf, cfg_.unpol_sf_obj,
                                  cfg_.pol_sf_obj, cfg_.r_source);
    dis_sampler_ = std::make_shared<InclusiveSampler>(kernel, beams_,
                                                      cfg_.scenario, cfg_.grid);
  }

  if (cfg_.channel == PipelineChannel::Inclusive) {
    GeneratorConfig gc;
    gc.with_virtual_photon = cfg_.with_virtual_photon;
    gc.channel = Channel::Inclusive;
    gen_.reset(new InclusiveGenerator(dis_sampler_, gc));
    cplan_.reserve(nc);
    // RUN-PLAN SPIN vs ION SPIN, on the INCLUSIVE branch too.  The tagged
    // branch has checked this since day one ("Pipeline: run-plan spin ... !=
    // channel ion spin"); here the mismatch used to surface two frames down
    // as `InclusiveKernel::amplitudes: spin state J = 1.000000 is not the
    // kernel's ion spin 1.500000`, which names neither the plan that built
    // the category nor the way out.  It is reachable because three of the
    // four standard plans hard-code j = 1 (`tensor_thirds_plan`,
    // `transverse_tensor_plan`, `tensor_flip_plan`), so
    // `transverse_tensor_plan(...)` on a 7Li config is a plausible thing to
    // write -- `make_plan` now refuses it, and this catches the C++ and the
    // hand-built route (phase_D_li7_rank2.md sec. 1.5, defect F2).
    const double j_ion = ion_by_name(cfg_.isotope).spin;
    for (const SpinCategory& c : plan_.categories()) {
      if (std::fabs(c.j - j_ion) > 1e-9) {
        throw std::runtime_error(
            "Pipeline: run-plan category \"" + c.name + "\" is J = " +
            fmt_g(c.j) + " but the inclusive kernel's ion " + cfg_.isotope +
            " is spin " + fmt_g(j_ion) +
            " -- tensor_thirds_plan, transverse_tensor_plan and "
            "tensor_flip_plan all hard-code j = 1, so they cannot run on a "
            "spin-3/2 beam; helicity_flip_plan is the only standard plan that "
            "takes j.  There is no spin-3/2 tensor plan in this tree, and on "
            "the inclusive channel there would be nothing for one to measure "
            "(inclusive_rank2_is_empty: the 7Li rank-2 slots are unset and "
            "the whole tensor sector is exactly 0)");
      }
    }
    for (std::size_t k = 0; k < nc; ++k) {
      cplan_.push_back(dis_sampler_->make_plan(plan_.categories()[k]));
      sigma_[k] = dis_sampler_->sigma_tot_pb(plan_.categories()[k]);
    }
  } else if (cfg_.channel == PipelineChannel::CoherentLi6) {
    // The coherent yield is f_coh(x) times the UNPOLARIZED inclusive rate
    // (`coherent.project_coherent`: n_coh = proj.n_events * f_coh(x)), so the
    // (x, Q2) cell weights are the sampler's own accepted cross sections
    // reweighted by f_coh at the cell center, and sigma_coh is spin
    // independent.  The tensor signal of this channel lives in the RECOIL
    // azimuth, 1 + c2 cos 2(phi_t - phi_S), which `CoherentSampler` owns;
    // folding the inclusive w_avg in on top would count the polarization twice.
    const std::vector<double>& sc = dis_sampler_->cell_xsec_pb();
    const std::vector<double>& xc = dis_sampler_->x_cells();
    const std::vector<double>& q2c = dis_sampler_->q2_cells();
    coh_cdf_.resize(sc.size());
    // The same weights, UNNORMALIZED and un-accumulated: `cell_rate_weights_pb`
    // hands them out as THIS channel's own per-cell rate, which is what any
    // "how much of this run sits below X" fraction must be taken against.
    // Before 2026-09-05 `unpol_sf_grid_report` took the inclusive cells
    // instead and a coherent run recorded 36.179 % where its own rate has
    // 44.746 % below CT18NLO's floor.
    coh_cell_pb_.assign(sc.size(), 0.0);
    double acc = 0.0;
    for (std::size_t c = 0; c < sc.size(); ++c) {
      // C1.  A cell only carries coherent rate if a diffractive system of at
      // least M_X,min FITS in it at a pomeron fraction inside the diffractive
      // region: x_P(M_X,min) <= x_P,max.  That is a kinematic statement, and
      // it cuts exactly where the coherent scenario has almost no rate left
      // anyway -- x >~ 0.05, where f_coh(x) is already 1/26 of its peak.
      // Without it the largest-x cells would be forced onto the degenerate
      // x_P = x_P,min branch, at which the nucleus loses ~9 % of its momentum
      // and the "intact recoil" leaves the near-beam band altogether.
      const double w2c = w2_from_xq2(xc[c], q2c[c]);
      const bool fits = cfg_.coherent_xpom.x_pom_min(q2c[c], w2c)
                        <= cfg_.coherent_xpom.x_pom_max;
      if (fits) {
        coh_cell_pb_[c] = sc[c] * cfg_.coherent.coherent_fraction(xc[c]);
        acc += coh_cell_pb_[c];
      }
      coh_cdf_[c] = acc;
    }
    if (!(acc > 0.0)) {
      throw std::runtime_error(
          "Pipeline: no coherent rate -- no accepted (x, Q2) cell admits a "
          "diffractive mass of at least CoherentXpomModel::m_x_min");
    }
    sigma_coh_pb_ = acc;
    for (double& v : coh_cdf_) v /= acc;
    coh_ms_.reserve(nc);
    for (std::size_t k = 0; k < nc; ++k) {
      const SpinCategory& c = plan_.categories()[k];
      coh_ms_.push_back(m_values(c.j));
      const double pzz = c.j >= 1.0 ? c.moments().tensor : 0.0;
      csampler_.emplace_back(new CoherentSampler(cfg_.coherent, p_u, pzz,
                                                 c.phi_s, beams_.ion.A,
                                                 beams_.ion.Z));
      csampler_.back()->set_optics(optics_, pot_config_);
      csampler_.back()->set_xpom_model(cfg_.coherent_xpom);
      // P4: THROWS if 1 + c2 cos 2phi would go negative anywhere on [0, t_max]
      // at this category's P_zz.
      csampler_.back()->set_t_max(cfg_.coherent_t_max);
      csampler_.back()->set_weighted_azimuth(cfg_.coherent_weighted_azimuth);
      sigma_[k] = sigma_coh_pb_;
    }
  }

  // --- rc (rc.hpp) --------------------------------------------------------
  // ANCHOR: this sits AFTER the channel if/else chain (tagged / Inclusive /
  // CoherentLi6) and BEFORE `// --- beams ---`, because that is the first
  // point at which `dis_sampler_`, `model_` and `plan_` are resolved on EVERY
  // channel.  It cannot go beside the FSI block: that block is INSIDE
  // `if (is_tagged(cfg_.channel))`, and on the inclusive and coherent channels
  // `dis_sampler_` is not assigned until the non-tagged else-branch -- building
  // `RcModel` there would dereference a null sampler on the channel this
  // design calls the home case.
  if (cfg_.rc != PipelineRc::Off) {
    const Channel rc_channel = event_channel_of(cfg_.channel);
    if (model_) {
      rc_ = std::make_shared<RcModel>(cfg_.rc, cfg_.rc_options, dis_sampler_,
                                      plan_, rc_channel, beams_.ion, model_);
    } else {
      rc_ = std::make_shared<RcModel>(cfg_.rc, cfg_.rc_options, dis_sampler_,
                                      plan_, rc_channel, beams_.ion);
    }
  }

  // --- beams --------------------------------------------------------------
  const Ion& ion = beams_.ion;
  ion_pdg_ = nuclear_pdg(ion.Z, ion.A);
  ion_charge_ = ion.Z;
  // For the tagged channels the beam ion mass MUST be the one the struck
  // cluster is built against (`struck_cluster` uses `ClusterChannel::m_beam`),
  // or P_ion - p_spec would not be the record's own P_X.  The two tables agree
  // to the last bit for every species the library carries; taking the
  // channel's copy makes the identity structural rather than lucky.
  ion_mass_ = channel_ ? channel_->base.m_beam() : ion.mass();
  beam_e_ = {beams_.electron_energy, 0.0, 0.0, -beams_.electron_energy};
  const double p_a = static_cast<double>(ion.A) * p_u;
  beam_ion_ = {std::sqrt(p_a * p_a + ion_mass_ * ion_mass_), 0.0, 0.0, p_a};

  // --- counts -------------------------------------------------------------
  counts_.assign(nc, 0);
  lumi_.assign(nc, 0.0);
  if (cfg_.n_events > 0) {
    // The largest-remainder split of `InclusiveGenerator::run_n`, reproduced
    // digit for digit so the inclusive path stays bit-identical to it.
    std::vector<double> rate(nc, 0.0);
    double tot = 0.0;
    for (std::size_t k = 0; k < nc; ++k) {
      rate[k] = plan_.categories()[k].lumi_fraction * sigma_[k];
      tot += rate[k];
    }
    if (!(tot > 0.0)) throw std::runtime_error("Pipeline: the plan has no rate");
    std::uint64_t assigned = 0;
    for (std::size_t k = 0; k < nc; ++k) {
      counts_[k] = static_cast<std::uint64_t>(
          std::floor(static_cast<double>(cfg_.n_events) * rate[k] / tot));
      assigned += counts_[k];
    }
    for (std::size_t k = 0; assigned < cfg_.n_events; k = (k + 1) % nc) {
      ++counts_[k];
      ++assigned;
    }
  } else {
    // C6.  `Optics::lumi_fraction` is the share of the machine luminosity the
    // configured far-forward working point delivers -- 1 for the Yellow
    // Report envelopes, 0.147 for the 6Li 5x41 tagging point, which buys its
    // acceptance by de-squeezing beta*_x.  It multiplies the COUNTS and never
    // the cross sections (`sigma_per_category_pb` stays share-invariant).
    lumi_ = plan_.lumi_share_vector(cfg_.lumi_pb * optics_lumi_factor());
    for (std::size_t k = 0; k < nc; ++k) {
      const double mu = lumi_[k] * sigma_[k];
      if (cfg_.poisson) {
        Rng counter(cfg_.seed, cfg_.run, k, kCountStreamEvent);
        counts_[k] = rng_poisson(counter, mu);
      } else {
        counts_[k] = static_cast<std::uint64_t>(std::llround(mu));
      }
    }
  }
  offset_.assign(nc + 1, 0);
  for (std::size_t k = 0; k < nc; ++k) offset_[k + 1] = offset_[k] + counts_[k];
  total_ = offset_[nc];
}

const TaggedChannel* Pipeline::tagged_channel() const { return channel_.get(); }

const CoherentSampler* Pipeline::coherent_sampler(std::size_t category) const {
  return category < csampler_.size() ? csampler_[category].get() : nullptr;
}

const std::vector<double>& Pipeline::cell_rate_weights_pb() const {
  // The coherent branch is the ONLY one whose rate is not the sampler's own
  // accepted cell cross sections -- see the header for the tagged case, whose
  // extra factors are per (M, m_S) and cell-independent, so they cancel out
  // of every fraction taken against this vector.
  if (cfg_.channel == PipelineChannel::CoherentLi6) return coh_cell_pb_;
  return dis_sampler_->cell_xsec_pb();
}

// ------------------------------------------------ the route-knob measurement
//
// `--optics`, `n_sigma` and `pot_config` reach ONE quantity of a run, the
// per-event `route` column (`fill_row` calls `route_of(ev, optics,
// pot_config)`; the coherent sampler stores the same label from the same
// classifier), and they reach it through exactly two branches of
// `route_charged`: `Optics::clears` is consulted only for a NEAR-BEAM fragment
// (|R - 1| < NEAR_BEAM_BAND, theta < THETA_RP_OUTER) and `over_rigid_route`
// only for an OVER-RIGID one (R > 1 + NEAR_BEAM_BAND).  So "would another
// value give another file" is a question about where THIS run's fragments
// fall, and MEASURED 2026-09-05 it is not answerable from the channel:
//
//   channel            60 ev seed 7                  2000 ev seed 11
//   coherent-6Li       optics tagging MOVES,         optics tagging MOVES,
//                      n_sigma 1 MOVES,              n_sigma 1/3 MOVE,
//                      pot_config NEITHER            pot_config NEITHER
//   tagged-6Li-alpha   all three MOVE                all three MOVE
//   tagged-7Li-alpha   only n_sigma 30 MOVES         all three MOVE
//   tagged-d-p         optics yr-high-divergence     NONE of the eight
//                      and n_sigma 30 MOVE           alternatives moves
//
// -- which is why the row is measured on the run rather than tabulated by
// channel, and why the old hand-written caveat ("all four envelope changes
// move the route column on tagged-6Li-alpha and NONE of them does on
// tagged-d-p at config 1") was a claim about one (channel, config, seed,
// event count) printed on every run.
//
// The measurement is the idiom the header states above `route_of`: the sample
// is drawn ONCE and re-routed per optics, because `route_of` is a pure
// function of the finished record and the envelope.  It costs one extra
// generation pass over the run, on the channels that write a far-forward
// fragment only, memoized per `Pipeline`, and it stops as soon as all three
// knobs have been seen to move.
void Pipeline::route_probe_event(std::uint64_t index, Event& out) const {
  const std::size_t k = category_of(index);
  const std::uint64_t local = index - offset_[k];
  Rng rng(cfg_.seed, cfg_.run, k, local);
  switch (cfg_.channel) {
    case PipelineChannel::Inclusive:
      make_inclusive(k, local, index, rng, out);
      break;
    case PipelineChannel::CoherentLi6:
      make_coherent(k, local, index, rng, out);
      break;
    default:
      make_tagged(k, local, index, rng, out);
      break;
  }
  // Deliberately NOT `rc_->fill(out)` and NOT `cfg_.hadronizer(out, rng)`:
  // neither moves a four-vector, the hadronizer is not re-entrant and its
  // `PythiaBridgeStats` are recorded in `meta`, so running it twice would
  // make the file describe a run that did not happen.
}

const Pipeline::RouteReach& Pipeline::route_reach() const {
  std::call_once(route_reach_once_, [this] {
    RouteReach& r = route_reach_;
    r.has_route = is_tagged(cfg_.channel) ||
                  cfg_.channel == PipelineChannel::CoherentLi6;
    if (!r.has_route || total_ == 0) return;
    const double p_u = beams_.ion_momentum_per_nucleon;
    // THE ALTERNATIVES, named once here and repeated in
    // python/tests/test_knob_provenance.py's route cells, which vary exactly
    // these values against the output hash.
    struct Alt { std::string name; Optics optics; };
    std::vector<Alt> opt_alts, ns_alts;
    const OpticsChoice choices[4] = {OpticsChoice::YellowReportHighAcceptance,
                                     OpticsChoice::YellowReportHighDivergence,
                                     OpticsChoice::Tagging,
                                     OpticsChoice::TaggingLegacyLevers};
    const char* choice_names[4] = {"yr-high-acceptance", "yr-high-divergence",
                                   "tagging", "tagging-legacy"};
    for (int i = 0; i < 4; ++i) {
      if (cfg_.optics_choice == choices[i]) continue;
      try {
        opt_alts.push_back({choice_names[i],
                            optics_for(choices[i], cfg_.isotope, p_u,
                                       cfg_.n_sigma)});
      } catch (const std::exception&) {
        // Not tabulated for this species (the tagging scan is published for
        // 6Li and 7Li only): not an alternative this run HAS.
      }
    }
    // `n_sigma` is not read at all under `OpticsChoice::Custom` -- the
    // `Optics` object is taken verbatim and no envelope is rebuilt -- so
    // there is nothing to probe there.
    if (cfg_.optics_choice != OpticsChoice::Custom) {
      for (double ns : {1.0, 3.0, 30.0}) {
        if (ns == cfg_.n_sigma) continue;
        try {
          ns_alts.push_back({fmt_g(ns), optics_for(cfg_.optics_choice,
                                                   cfg_.isotope, p_u, ns)});
        } catch (const std::exception&) {
        }
      }
    }
    std::vector<std::string> pot_alts;
    for (const char* pc : {"5x41", "10x100", "18x275"})
      if (pot_config_ != pc) pot_alts.emplace_back(pc);
    auto join = [](const std::vector<std::string>& v) {
      std::string s;
      for (std::size_t i = 0; i < v.size(); ++i)
        s += (i ? ", " : "") + v[i];
      return s.empty() ? std::string("none") : s;
    };
    std::vector<std::string> on, nn;
    for (const Alt& a : opt_alts) on.push_back(a.name);
    for (const Alt& a : ns_alts) nn.push_back(a.name);
    r.optics_tried = join(on);
    r.n_sigma_tried = join(nn);
    r.pot_tried = join(pot_alts);

    Event ev;
    for (std::uint64_t i = 0; i < total_; ++i) {
      if (r.optics_moves && r.n_sigma_moves && r.pot_moves) break;
      route_probe_event(i, ev);
      ++r.n_probed;
      const int base = route_of(ev, optics_, pot_config_);
      for (const Alt& a : opt_alts) {
        if (route_of(ev, a.optics, pot_config_) == base) continue;
        if (!r.optics_moves) { r.optics_moves = true; r.optics_alt = a.name; }
        if (a.name == r.optics_alt) ++r.optics_n;
      }
      for (const Alt& a : ns_alts) {
        if (route_of(ev, a.optics, pot_config_) == base) continue;
        if (!r.n_sigma_moves) { r.n_sigma_moves = true; r.n_sigma_alt = a.name; }
        if (a.name == r.n_sigma_alt) ++r.n_sigma_n;
      }
      for (const std::string& pc : pot_alts) {
        if (route_of(ev, optics_, pc) == base) continue;
        if (!r.pot_moves) { r.pot_moves = true; r.pot_alt = pc; }
        if (pc == r.pot_alt) ++r.pot_n;
      }
    }
  });
  return route_reach_;
}

double Pipeline::coherent_rate_x_edge() const {
  if (cfg_.channel != PipelineChannel::CoherentLi6) return 0.0;
  // THE EDGE IS A CELL CENTRE, NOT A CELL BOUNDARY.  `InclusiveSampler`
  // admits a cell by `in_acceptance(x_c, q2_c)` on its CENTRE
  // (src/core/sampler.cpp) and the coherent branch never re-tests x per
  // event, so the quantity an `x_max` is actually compared against is the
  // largest rate-carrying centre (0.0954992586 at config 1), not the top of
  // that cell (exp(logx_hi) = 0.1 + 1 ulp).  Measured 2026-09-05: every
  // x_max in (0.09549926, 0.1] is bit-identical to the default (1822 cells,
  // 1726 carrying, sigma 12476.025113185518), and 0.09549 is the first value
  // that moves (1770 / 1702, 12475.913921); the boundary rule called 0.097
  // "read" on a run identical to the default.
  auto edge_of = [](const std::vector<double>& w,
                    const std::vector<double>& centre) {
    double e = 0.0;
    for (std::size_t i = 0; i < w.size() && i < centre.size(); ++i)
      if (w[i] > 0.0) e = std::max(e, centre[i]);
    return e;
  };
  const PipelineConfig d;
  // THE UNCLIPPED CELL SET, and the whole point of this function.  Measuring
  // the edge on `cell_rate_weights_pb()` -- this run's own accepted cells --
  // is self-referential: at `--x-max 0.05` the surviving cells all lie below
  // 0.05 by construction, so the edge comes back 0.047863 and the rule
  // concluded "0.05 is above the edge, it clips nothing" about a value that
  // had just removed 320 of 1726 rate-carrying cells and moved sigma from
  // 12476.025113 to 12446.109030 pb.  The edge belongs to the WINDOW AT THE
  // SHIPPED x_max, which is a property of the channel and of the rest of the
  // scenario, and not of the knob being judged.
  const double ref_x_max = std::max(d.scenario.x_max, cfg_.scenario.x_max);
  if (cfg_.scenario.x_max >= ref_x_max) {
    return edge_of(cell_rate_weights_pb(), dis_sampler_->x_cells());
  }
  Scenario s = cfg_.scenario;
  s.x_max = ref_x_max;
  const InclusiveSampler ref(dis_sampler_->kernel_ptr(), beams_, s, cfg_.grid);
  const std::vector<double>& sc = ref.cell_xsec_pb();
  const std::vector<double>& xc = ref.x_cells();
  const std::vector<double>& q2c = ref.q2_cells();
  std::vector<double> w(sc.size(), 0.0);
  for (std::size_t c = 0; c < sc.size(); ++c) {
    // The SAME gate the constructor's coherent branch applies, read from the
    // same two objects: a cell carries coherent rate only if a diffractive
    // system of at least M_X,min fits in it at x_P <= x_P,max.
    const double w2c = w2_from_xq2(xc[c], q2c[c]);
    if (cfg_.coherent_xpom.x_pom_min(q2c[c], w2c) <=
        cfg_.coherent_xpom.x_pom_max) {
      w[c] = sc[c] * cfg_.coherent.coherent_fraction(xc[c]);
    }
  }
  return edge_of(w, xc);
}

// --------------------------------------------------- the knob-provenance table
//
// One row per user-settable knob, in the order a reader meets them: beams and
// statistics, the spin fill, the acceptance window, the structure-function
// selectors, the tagged cluster, FSI, RC, the coherent channel and the T2
// tier.  Every `NotRead` reason opens with "not read ..." because it is the
// string `meta` writes IN PLACE OF the value (`KnobProvenance::meta_value`),
// so it has to read as a sentence about this run on its own.
//
// THE REACH RULES ARE MEASURED, not inferred.  Each one is the (channel x
// plan x knob) cell of the 2026-09-05 matrix -- 12 channel x plan specs, up
// to 38 knob cells each, 600 events at seed 7, sha256 over all 47 ndarray
// columns plus sigma_pb and sigma_per_category_pb against a same-plan
// baseline -- rebuilt as an executable assertion in
// python/tests/test_knob_provenance.py.  The numbers behind each sentence are
// in docs/open_items/run_2026-09-03/phase_D_numbers.md sec. D6.
std::vector<KnobProvenance> Pipeline::knob_provenance(
    const KnobRunContext& ctx) const {
  const PipelineConfig d;          // the shipped defaults, read not retyped
  const PipelineConfig& c = cfg_;
  std::vector<KnobProvenance> rows;
  auto add = [&rows](const char* name, const char* flag, std::string value,
                     KnobStatus st, std::string reason, bool at_default) {
    KnobProvenance r;
    r.name = name;
    r.flag = flag;
    r.value = std::move(value);
    r.status = st;
    r.reason = std::move(reason);
    // THE LABEL IS DERIVED, ONCE, HERE.  Every `NotRead` reason below opens
    // with its scope clause -- "not read on channel X", "not read by plan Y",
    // "not read at Z" -- and then explains; the clause alone is what a `meta`
    // key carries in place of the value, which is how `meta["pol_sf"]` on a
    // coherent run still reads exactly "not read on channel coherent-6Li".
    // Two separators, and only these two: ": " and " -- " (a bare ':' would
    // split "Tier::T0" in half).
    if (st == KnobStatus::NotRead) {
      const std::size_t a = r.reason.find(": ");
      const std::size_t b = r.reason.find(" -- ");
      const std::size_t cut = std::min(a, b);
      r.label = (cut == std::string::npos) ? r.reason : r.reason.substr(0, cut);
    }
    r.at_default = at_default;
    rows.push_back(std::move(r));
  };
  auto yn = [](bool b) { return std::string(b ? "true" : "false"); };
  const std::string chn = pipeline_channel_name(c.channel);
  const std::string off_ch = "not read on channel " + chn + ": ";
  // Which of `--pz` / `--pzz` / `--rel-lumi-offset` the NAMED plan factory
  // takes.  bookkeeping.hpp's four factories are the whole of it:
  // tensor_thirds_plan(pz, pzz, rel_lumi_offset, ...),
  // helicity_flip_plan(j, pz, pe, opt{rel_lumi_offset, ...}),
  // transverse_tensor_plan(pzz, phi_s),
  // tensor_flip_plan(pzz, phi_s, share_plus, rel_lumi_offset).
  // `--pe` needs no entry: `helicity_flip_plan` is the only factory that
  // produces lam_e != 0, so the FILL answers for it exactly.
  const std::string pn = ctx.plan_name;
  const bool is_thirds = (pn == "tensor-thirds" || pn == "azz");
  const bool is_flip = (pn == "helicity-flip" || pn == "apar");
  const bool is_perp = (pn == "transverse-tensor" || pn == "cos2phi");
  const bool is_tflip = (pn == "tensor-flip" || pn == "flip");
  const bool named = is_thirds || is_flip || is_perp || is_tflip;
  const bool plan_reads_pz = is_thirds || is_flip;
  const bool plan_reads_rel = is_thirds || is_flip || is_tflip;
  // DID THIS FILL COME FROM THE LADDER?  Measured from the plan itself, not
  // taken from the context, because a C++ caller supplies no `pzz_mode` and
  // a guess is exactly what `KnobRunContext`'s header forbids: the ladder
  // branch IS `populations_maxent(J, P_z)`, so comparing this run's own
  // populations against it answers the question with no second source of
  // truth.  It also gets the coincidence right in the only way the matrix
  // accepts: a typed P_zz that reproduces the ladder's populations double
  // for double did not move the file, and `--pzz` is then NOT read whatever
  // was typed.
  //
  // Since 2026-09-06 `--pzz-mode typed` makes `helicity_flip_plan` read
  // `--pzz` (`HelicityFlipOptions::use_explicit_pzz`), so this row is no
  // longer "the three tensor plans and nobody else".
  const double flip_j =
      plan_.categories().empty() ? 0.0 : plan_.categories().front().j;
  bool fill_is_ladder = false;
  if (is_flip && !plan_.categories().empty()) {
    try {
      fill_is_ladder = plan_.categories().front().populations ==
                       populations_maxent(flip_j, plan_.pz_true());
    } catch (const std::exception&) {
      fill_is_ladder = false;
    }
  }
  // A fill below spin 1 has NO rank-2 moment at all, so neither branch of
  // `helicity_flip_plan` reads `--pzz` there (the j = 1/2 branch ignores it).
  const bool flip_reads_pzz =
      is_flip && flip_j >= 1.0 - 1e-9 && !fill_is_ladder;
  const bool plan_reads_pzz = is_thirds || is_perp || is_tflip || flip_reads_pzz;
  const std::string flip_fill_clause =
      ctx.pzz_mode.empty() ? std::string("a max-entropy fill")
                           : ("--pzz-mode " + ctx.pzz_mode);
  const std::string by_plan =
      named ? ("not read by plan " + pn + ": ")
            : std::string("not read by this run's plan: ");

  // ------------------------------------------------ beams and statistics
  add("isotope", "--isotope", c.isotope, KnobStatus::Read,
      "the beam species: every kernel, every mass, the optics row and the "
      "channel's own ion are built from it",
      c.isotope == d.isotope);
  add("beam_config", "--config", std::to_string(c.beam_config),
      KnobStatus::Read,
      "index into default_configs(" + c.isotope + "): it sets E_e = " +
          fmt_g(beams_.electron_energy) + " GeV, p/u = " +
          fmt_g(beams_.ion_momentum_per_nucleon) + " GeV and s/u = " +
          fmt_g(beams_.s_per_nucleon()) + " GeV^2",
      c.beam_config == d.beam_config);
  add("channel", "--channel", chn, KnobStatus::Read,
      "the physics channel -- it is what every reach rule below is scoped by",
      c.channel == d.channel);
  // THE THREE ROUTE KNOBS, MEASURED ON THIS RUN (`RouteReach`, above).
  //
  // They are consulted by `route_of`, and `route_of` returns `kRouteLost`
  // BEFORE it looks at the envelope when the event carries neither a tagged
  // spectator nor an intact recoil -- which is every INCLUSIVE event.
  // Measured 2026-09-05: on `--channel inclusive` (6Li, 7Li and d, 60 and
  // 2000 events) the `route` column is 0 = Route::Lost for every event, and
  // all eight alternative envelopes / machine configurations are BIT-IDENTICAL
  // to the default in all 47 columns and every sigma.
  //
  // ON A CHANNEL THAT DOES WRITE A FRAGMENT the answer is no longer a
  // property of the channel, so it is re-routed rather than asserted: the old
  // row said READ on all three of them everywhere a route exists and hung a
  // hand-written caveat about tagged-d-p on the read side, which the matrix
  // then flagged as 22 "did not move but read" cells -- `--pot-config` moves
  // NOTHING on the coherent channel at any statistics (the intact recoil is
  // never over-rigid, so `over_rigid_route` is never reached) and nothing on
  // tagged-d-p, whose proton spectator is outside the beam band by rigidity.
  const RouteReach& rr = route_reach();
  const bool has_route = rr.has_route;
  const std::string no_route =
      off_ch +
      "no event of this channel carries a far-forward fragment (no tagged "
      "spectator, no intact recoil), so route_of returns Route::Lost before "
      "it looks at an envelope -- measured, the whole `route` column is 0";
  // THE ONE NON-ROUTE READER OF `--optics`, and it is a different question:
  // in LUMINOSITY mode `apply_optics_lumi_fraction` multiplies every
  // category's count by `Optics::lumi_fraction`, which the tagging working
  // point carries and the Yellow Report envelopes do not.  That reaches every
  // column of the file on every channel, the inclusive one included, so it is
  // tested before the route probe and not after it.
  const bool optics_sets_counts = (c.n_events == 0) &&
                                  c.apply_optics_lumi_fraction;
  const std::string probe_tail =
      " -- MEASURED on this run by re-routing its own " +
      std::to_string(rr.n_probed) +
      " events (the sample is drawn ONCE and re-routed per envelope: "
      "route_of is a pure function of the finished record)";
  add("optics", "--optics", optics_.name,
      (optics_sets_counts || rr.optics_moves) ? KnobStatus::Read
                                              : KnobStatus::NotRead,
      optics_sets_counts
          ? std::string("in this run's LUMINOSITY mode it reaches every "
                        "column: apply_optics_lumi_fraction multiplies every "
                        "category's count by Optics::lumi_fraction = " +
                        fmt_g(optics_lumi_factor()) +
                        (rr.optics_moves
                             ? std::string(", and it moves the route label "
                                           "of ") + std::to_string(rr.optics_n)
                                   + " of this run's events at " +
                                   rr.optics_alt
                             : std::string("")))
      : rr.optics_moves
          ? std::string("the far-forward envelope the route label of every "
                        "event is priced at (route_of); it moves no "
                        "four-vector") +
                probe_tail + ": " + rr.optics_alt + " moves " +
                std::to_string(rr.optics_n) + " of them"
      : !has_route
          ? no_route
          // THE NAME STAYS IN THE LABEL, and therefore in `meta["optics"]`
          // and in the banner header (2026-09-05, phase F).  This branch is
          // the one where the classifier IS consulted -- the channel writes a
          // far-forward fragment and every event's route label is priced at
          // THIS envelope -- and only the SENSITIVITY to the envelope is
          // absent on this sample.  A bare "not read on this run" said the
          // stronger thing and dropped the envelope name that the a94fd6e
          // files carried and that the tag fraction printed three lines below
          // it in the banner is quoted at.  The status stays NotRead, because
          // the output hash is what the matrix measures; the label carries
          // the name AND the scope, which is what a reader needs.
          : std::string("not read on this run at ") + optics_.name +
                " (the classifier IS consulted here; this run's route labels "
                "are insensitive to " + rr.optics_tried + ")" + probe_tail +
                ", every one of them keeps its route label under "
                "every other tabulated envelope (" +
                rr.optics_tried +
                ").  The classifier IS consulted here; what it is "
                "not is sensitive to the envelope on this sample",
      c.optics_choice == d.optics_choice);
  add("n_sigma", "", fmt_g(c.n_sigma),
      rr.n_sigma_moves ? KnobStatus::Read : KnobStatus::NotRead,
      rr.n_sigma_moves
          ? std::string("beam-exclusion half-width the tabulated envelope is "
                        "built at (optics_for)") +
                probe_tail + ": n_sigma = " + rr.n_sigma_alt + " moves " +
                std::to_string(rr.n_sigma_n) + " route labels"
      : !has_route ? no_route
      : c.optics_choice == OpticsChoice::Custom
          ? std::string("not read at optics_choice = Custom: the Optics "
                        "object in PipelineConfig::optics is taken verbatim "
                        "and no envelope is rebuilt")
          // The same rule as `optics` above: the envelope every route label
          // IS priced at was built at THIS n_sigma, so the label names it.
          : std::string("not read on this run at n_sigma = ") +
                fmt_g(c.n_sigma) +
                " (the envelope the classifier consults was built at it; "
                "this run's route labels are insensitive to " +
                rr.n_sigma_tried + ")" + probe_tail +
                ", every one of them keeps its route label at n_sigma = " +
                rr.n_sigma_tried,
      c.n_sigma == d.n_sigma);
  add("pot_config", "", pot_config_,
      rr.pot_moves ? KnobStatus::Read : KnobStatus::NotRead,
      rr.pot_moves
          ? std::string("the machine configuration the OVER-RIGID branch of "
                        "the route classification tests against") +
                probe_tail + ": pot_config = " + rr.pot_alt + " moves " +
                std::to_string(rr.pot_n) + " route labels"
      : !has_route
          ? no_route
          // The same rule as `optics` above, and the same reason: this run's
          // route classification consults the machine configuration, so the
          // label names the one it consulted (`meta["pot_config"]` carried a
          // bare "10x100" before 2026-09-05 and must not lose it).
          : std::string("not read on this run at ") + pot_config_ +
                " (the branch it enters is never reached on this run)" +
                probe_tail +
                ", none of them is over-rigid enough to reach the "
                "only branch it has: pot_config enters route_charged "
                "through over_rigid_route alone, which is tested at "
                "R > 1 + NEAR_BEAM_BAND, and every machine "
                "configuration (" + rr.pot_tried +
                ") leaves this run's route column bit-identical",
      c.pot_config == d.pot_config);
  add("seed", "--seed", std::to_string(c.seed), KnobStatus::Read,
      "every event is Rng(seed, run, category, local index): the whole "
      "sample is a pure function of it", c.seed == d.seed);
  add("run", "--run", std::to_string(c.run), KnobStatus::Read,
      "the second word of every event's counter-based stream", c.run == d.run);
  const bool fixed_count = c.n_events > 0;
  add("events", "--events", std::to_string(c.n_events),
      fixed_count ? KnobStatus::Read : KnobStatus::NotRead,
      fixed_count ? std::string("fixed-count mode: the total splits across "
                                "categories in proportion to lumi_fraction x "
                                "sigma")
                  : std::string("not read in LUMINOSITY mode (lumi_pb = " +
                                fmt_g(c.lumi_pb) +
                                "): the counts come from the luminosity, and "
                                "validate() refuses the two together"),
      c.n_events == d.n_events);
  add("lumi_pb", "--lumi", fmt_g(c.lumi_pb),
      fixed_count ? KnobStatus::NotRead : KnobStatus::Read,
      fixed_count ? std::string("not read in FIXED-COUNT mode (events = " +
                                std::to_string(c.n_events) +
                                "): the count is given, and validate() "
                                "refuses the two together")
                  : std::string("luminosity mode: the per-category counts are "
                                "drawn from lumi_fraction x lumi x sigma"),
      c.lumi_pb == d.lumi_pb);
  add("poisson", "", yn(c.poisson),
      fixed_count ? KnobStatus::NotRead : KnobStatus::Read,
      fixed_count ? std::string("not read in fixed-count mode: there is "
                                "nothing to fluctuate when the count is given")
                  : std::string("Poisson-fluctuates the per-category counts of "
                                "the luminosity mode"),
      c.poisson == d.poisson);
  add("apply_optics_lumi_fraction", "", yn(c.apply_optics_lumi_fraction),
      fixed_count ? KnobStatus::NotRead : KnobStatus::Read,
      fixed_count
          ? std::string("not read in fixed-count mode: it multiplies the "
                        "COUNTS and never a cross section (bookkeeping.hpp), "
                        "and the count is given here")
          : std::string("multiplies every category's luminosity by "
                        "Optics::lumi_fraction = " +
                        fmt_g(optics_lumi_factor()) + " -- COUNTS only"),
      c.apply_optics_lumi_fraction == d.apply_optics_lumi_fraction);
  add("with_virtual_photon", "", yn(c.with_virtual_photon), KnobStatus::Read,
      "writes the exchanged photon as a status-3 documentation particle in "
      "every event record", c.with_virtual_photon == d.with_virtual_photon);

  // ------------------------------------------------------- the spin fill
  //
  // A fill moment has no "shipped default": every value of --pz / --pzz /
  // --pe is a deliberate choice, so `at_default` is false on all three and
  // the banner always names them.
  if (named) {
    const double pz_v = std::isnan(ctx.pz) ? plan_.pz_true() : ctx.pz;
    add("pz", "--pz", fmt_g(pz_v),
        plan_reads_pz ? KnobStatus::Read : KnobStatus::NotRead,
        plan_reads_pz
            ? ("the fill's vector moment; this run's plan carries P_z = " +
               fmt_g(plan_.pz_true()))
            : (by_plan + "its factory takes P_zz and the azimuth only "
                         "(bookkeeping.hpp), so the fill's own vector moment "
                         "is " + fmt_g(plan_.pz_true())),
        false);
    const double pzz_v = std::isnan(ctx.pzz) ? plan_.pzz_true() : ctx.pzz;
    add("pzz", "--pzz", fmt_g(pzz_v),
        plan_reads_pzz ? KnobStatus::Read : KnobStatus::NotRead,
        plan_reads_pzz
            ? (std::string(is_flip
                               ? "the fill's rank-2 moment, read here because "
                                 "this fill is NOT the max-entropy ladder: "
                                 "--pzz-mode typed built it with "
                                 "spin1_populations / spin32_populations at "
                                 "the typed value, and this run's plan "
                                 "carries "
                               : "the fill's rank-2 moment; this run's "
                                 "plan carries ") +
               std::string(plan_.categories().empty() ||
                                   std::fabs(plan_.categories().front().j -
                                             1.5) > 1e-9
                               ? "P_zz = "
                               : "T = ") +
               fmt_g(plan_.pzz_true()))
            : (is_flip
                   ? ("not read by plan " + pn + " at " + flip_fill_clause +
                      ": helicity_flip_plan leaves HelicityFlipOptions::"
                      "use_explicit_pzz false and builds its fill from the "
                      "MAX-ENTROPY ladder at --pz, so the alignment this run "
                      "carries is " + fmt_g(plan_.pzz_true()) + " and not "
                      "what was typed" +
                      (fill_is_ladder
                           ? std::string("; MEASURED, this run's populations "
                                         "ARE populations_maxent(J, --pz) "
                                         "double for double")
                           : std::string("; this fill carries no rank-2 "
                                         "moment at all below spin 1")) +
                      ".  --pzz-mode typed honours the typed value, and "
                      "refuses it outside the plan's domain")
                   : (by_plan + "its factory builds no rank-2 moment "
                                "(bookkeeping.hpp) -- unreachable today: the "
                                "other three factories all take P_zz")),
        false);
    // `--pzz-mode`, THE FILL SELECTOR (2026-09-06).  Which of
    // `helicity_flip_plan`'s two branches built this run's populations:
    // "ladder" (`use_explicit_pzz = false`, the max-entropy fill at --pz --
    // the shipped default and bit for bit the tree before the flag existed)
    // or "typed" (`= true`, the typed --pzz through spin1_populations /
    // spin32_populations, which REFUSE a value outside the plan's domain
    // rather than clamp to its edge).
    //
    // MEASURED, NOT TABULATED -- the `RouteReach` rule applied to a fill.
    // The row rebuilds the OTHER mode's plan at this run's own (J, P_z, P_e)
    // and compares the populations it would have used, doubles to doubles,
    // because "another value would give another file" is a question about
    // THIS fill and not about plan names: at J = 1 the populations are fixed
    // uniquely by (P_z, P_zz), so a typed P_zz that reproduces the ladder's
    // own alignment reproduces the whole fill and the row must say not-read,
    // while at J = 3/2 the two differ in R_3 even at equal T.  A counter-mode
    // the domain REFUSES is read: it does not give this file either.
    if (!ctx.pzz_mode.empty()) {
      // A CONTEXT THAT NAMES A MODE THAT DOES NOT EXIST IS A CALLER ERROR,
      // and it is refused here rather than reported: every value this table
      // writes is a claim about what the run did, so a third string silently
      // reported as the ladder would be exactly the defect the table exists
      // to prevent.  Empty stays legal -- it means "not supplied", and the
      // row is then omitted.
      if (ctx.pzz_mode != "ladder" && ctx.pzz_mode != "typed") {
        throw std::runtime_error(
            "KnobRunContext::pzz_mode is \"" + ctx.pzz_mode +
            "\"; the only fills helicity_flip_plan has are \"ladder\" "
            "(HelicityFlipOptions::use_explicit_pzz false, the default) and "
            "\"typed\" (true).  Leave it empty to omit the row instead");
      }
      const bool typed_now = (ctx.pzz_mode == "typed");
      KnobStatus st = KnobStatus::NotRead;
      std::string why;
      if (!is_flip) {
        why = by_plan +
              "only helicity_flip_plan has two fills to choose between "
              "(HelicityFlipOptions::use_explicit_pzz, bookkeeping.hpp); the "
              "three tensor factories build their categories from the typed "
              "P_zz directly and have no max-entropy branch, so both modes "
              "leave this run's fill exactly as it is";
      } else {
        HelicityFlipOptions other_opt;
        other_opt.use_explicit_pzz = !typed_now;   // the OTHER mode
        other_opt.pzz = std::isnan(ctx.pzz) ? plan_.pzz_true() : ctx.pzz;
        const std::string rank2 =
            (std::fabs(flip_j - 1.5) < 1e-9) ? "T" : "P_zz";
        std::string other_side;
        std::string counted = "the other mode's fill cannot be built";
        bool moves = true;
        try {
          const RunPlan other = helicity_flip_plan(
              flip_j, plan_.pz_true(), plan_.pe_true(), other_opt);
          const std::vector<double>& a =
              plan_.categories().front().populations;
          const std::vector<double>& b = other.categories().front().populations;
          std::size_t ndiff = (a.size() != b.size()) ? a.size() : 0;
          if (a.size() == b.size()) {
            for (std::size_t i = 0; i < a.size(); ++i)
              if (a[i] != b[i]) ++ndiff;
          }
          moves = ndiff != 0;
          // THE COUNT, not the rounded moment: the two fills can print the
          // same P_zz at %g and still be different doubles (a typed P_zz set
          // to the ladder's own moment is built by the closed form of
          // `spin1_populations` and the ladder by `populations_maxent`'s
          // bisection), and it is the POPULATIONS that the weights are built
          // from.
          counted = std::to_string(ndiff) + " of " +
                    std::to_string(a.size()) +
                    " populations are different doubles";
          other_side = std::string(typed_now ? "the max-entropy ladder at "
                                               "this P_z would give "
                                             : "the typed P_zz would give ") +
                       rank2 + " = " + fmt_g(other.pzz_true());
        } catch (const std::exception& e) {
          // The counter-value is outside the plan's domain.  That is not
          // "did not move": it is a fill this run cannot have at all.
          other_side = std::string(typed_now ? "the max-entropy ladder"
                                             : "the typed P_zz") +
                       " is REFUSED at this fill -- " + std::string(e.what());
          moves = true;
        }
        st = moves ? KnobStatus::Read : KnobStatus::NotRead;
        why = moves
                  ? ("which of helicity_flip_plan's two fills this run "
                     "carries: " +
                     std::string(typed_now
                                     ? "typed reads --pzz through "
                                       "spin1_populations / "
                                       "spin32_populations"
                                     : "ladder builds populations_maxent at "
                                       "--pz and does NOT read --pzz") +
                     ", and this run's alignment is " + rank2 + " = " +
                     fmt_g(plan_.pzz_true()) + " while " + other_side +
                     " -- MEASURED here by rebuilding the other mode's fill: " +
                     counted)
                  : ("not read at this fill: both modes build the SAME "
                     "populations here, double for double (" + counted +
                     "; " + other_side +
                     "), so the fill -- and every weight built from it -- is "
                     "the same whichever mode is named");
      }
      add("pzz_mode", "--pzz-mode", ctx.pzz_mode, st, why,
          ctx.pzz_mode == "ladder");
    }
    if (!std::isnan(ctx.rel_lumi_offset)) {
      add("rel_lumi_offset", "--rel-lumi-offset", fmt_g(ctx.rel_lumi_offset),
          plan_reads_rel ? KnobStatus::Read : KnobStatus::NotRead,
          plan_reads_rel
              ? std::string("the relative-luminosity offset on this plan's "
                            "own category (bookkeeping.hpp sec. 5.0)")
              : (by_plan + "transverse_tensor_plan takes P_zz and phi_S only "
                           "and has no offset argument"),
          ctx.rel_lumi_offset == 0.0);
    }
  }
  {
    const bool pe_read = plan_has_beam_helicity(plan_);
    // NOTE which predicate this is.  `--pe` is read wherever the fill carries
    // a beam HELICITY at all, `--pol-sf` only where lam_e * P_e is non-zero:
    // at --pe 0 under helicity-flip the plan DID consult --pe (it is what
    // made the product zero) while g1 was never evaluated.  Two rules, one
    // for each question.
    bool any_lam = false;
    for (const SpinCategory& sc : plan_.categories())
      if (sc.lam_e != 0) any_lam = true;
    const double pe_v = std::isnan(ctx.pe) ? plan_.pe_true() : ctx.pe;
    add("pe", "--pe", fmt_g(pe_v),
        any_lam ? KnobStatus::Read : KnobStatus::NotRead,
        any_lam
            ? ("the beam polarization of this fill's helicity categories; "
               "lam_e * P_e is what g1 is multiplied by, and it is " +
               std::string(pe_read ? "non-zero here" : "ZERO here"))
            : (by_plan +
               "every category is built at lam_e = 0 (bookkeeping.cpp: the "
               "three tensor plans hard-code an unpolarised beam), so P_e "
               "multiplies nothing -- the per-event `pe` column records 0"),
        false);
  }

  // -------------------------------------------------- the acceptance window
  {
    KnobStatus st = KnobStatus::Read;
    std::string why =
        "the upper x edge of the accepted (x, Q2) window, applied before any "
        "rate is computed (Scenario::x_max; every other Scenario field is "
        "API-only and is read on every channel)";
    if (c.channel == PipelineChannel::CoherentLi6) {
      // THE EDGE, measured, and measured on the UNCLIPPED CELL SET
      // (`coherent_rate_x_edge`): the coherent channel carries rate only in
      // cells whose diffractive-mass gate x_P(M_X,min) <= x_P,max passes, and
      // that gate is (x, Q2)-dependent.
      //
      // TWO THINGS WERE WRONG WITH THE 2026-09-05 FORM OF THIS RULE, one at
      // each edge, and both are why the comparison is written the way it is.
      // (1) It read the edge off `cell_rate_weights_pb()`, i.e. off THIS
      // RUN's already-clipped cells, so it was self-referential: at
      // `--x-max 0.05` the surviving cells lie below 0.05 by construction,
      // the edge came back 0.047863, and the row said NOT READ of a value
      // that had just removed 320 of the 1726 rate-carrying cells and moved
      // sigma from 12476.025113 to 12446.109030 pb.  (2) It compared
      // exactly, and the edge is `exp(logx_hi)` of a grid boundary: at
      // `--x-max 0.1` the edge comes back 0.10000000000000002, one ulp above
      // the value that produced it, so a run BIT-IDENTICAL to the default
      // said READ.  Measured on the unclipped set with a relative tolerance,
      // 0.05 and 0.09 are READ (both move sigma) and 0.10, 0.95 and 1.0 are
      // NOT READ (all three bit-identical to the default).
      // (3) 2026-09-05, later: the edge was the top of the last carrying
      // CELL (exp(logx_hi)), but the sampler admits a cell on its CENTRE and
      // the coherent branch never re-tests x per event, so every x_max in
      // (largest carrying centre, cell top] -- 0.0955 ... 0.1 at config 1, a
      // 4.5 %-wide window -- was bit-identical to the default and said READ.
      // `coherent_rate_x_edge` now returns the largest rate-carrying centre.
      const double edge = coherent_rate_x_edge();
      const double tol = 1e-12 * std::max(1.0, edge);
      if (c.scenario.x_max >= edge - tol) {
        st = KnobStatus::NotRead;
        why = off_ch +
              "every cell that carries coherent rate has its CENTRE at or "
              "below x = " +
              fmt_g(edge) +
              " (the sampler admits a cell on its centre; the diffractive-mass "
              "gate x_P(M_X,min) <= x_P,max, Pipeline's coherent branch, "
              "measured on the UNCLIPPED cell set -- the window at the shipped "
              "x_max and not this run's own already-clipped one), and x_max = " +
              fmt_g(c.scenario.x_max) +
              " is at or above that centre, so it clips nothing.  A value "
              "BELOW it would move this channel";
      } else {
        why = "the upper x edge of the accepted (x, Q2) window: x_max = " +
              fmt_g(c.scenario.x_max) +
              " is BELOW this channel's largest rate-carrying cell centre x = " +
              fmt_g(edge) +
              " (measured on the unclipped cell set), so it removes cells "
              "that carry coherent rate and moves sigma_pb with them";
      }
    }
    add("x_max", "--x-max", fmt_g(c.scenario.x_max), st, why,
        c.scenario.x_max == d.scenario.x_max);
  }

  // ------------------------------------------- structure-function selectors
  if (c.kernel) {
    add("kernel", "", "caller-supplied InclusiveKernel",
        is_tagged(c.channel) ? KnobStatus::Refused : KnobStatus::Read,
        is_tagged(c.channel)
            ? (off_ch +
               "the tagged channels draw from the STRUCK-CLUSTER sampler and "
               "`kernel` reaches the ion-level one alone; validate() refuses "
               "it here")
            : std::string("the ion-level kernel of this run: it WINS over "
                          "--unpol-sf / --pol-sf / --b1-model, which "
                          "validate() therefore refuses beside it"),
        false);
  } else {
    add("kernel", "", "null (default_inclusive_kernel)",
        is_tagged(c.channel) ? KnobStatus::Refused : KnobStatus::Read,
        is_tagged(c.channel)
            ? (off_ch +
               "a caller-supplied kernel is not read on a tagged channel and "
               "validate() refuses one; the struck-cluster kernel takes its "
               "backends from --unpol-sf / --pol-sf")
            : std::string("null means the pipeline builds the kernel itself "
                          "from --b1-model / --unpol-sf / --pol-sf"),
        true);
  }
  add("unpol_sf", "--unpol-sf",
      c.kernel ? std::string("caller-supplied kernel")
               : std::string(unpol_sf_name(c.unpol_sf)),
      c.kernel ? KnobStatus::Refused : KnobStatus::Read,
      c.kernel
          ? std::string("a caller-supplied kernel carries its own f2_source, "
                        "so validate() refuses any non-toy selector beside "
                        "one")
          : std::string("-> InclusiveKernel::Options::f2_source on EVERY "
                        "kernel this run builds: the unpolarised rate, the "
                        "D_phi denominator of the tensor weight, the T1 "
                        "species draw and (through the CLI) the T2 one"),
      c.unpol_sf == d.unpol_sf);
  {
    const bool read = pol_sf_is_read(c, plan_);
    add("pol_sf", "--pol-sf",
        c.kernel ? std::string("caller-supplied kernel")
                 : std::string(pol_sf_name(c.pol_sf)),
        c.kernel ? KnobStatus::Refused
                 : (read ? KnobStatus::Read : KnobStatus::NotRead),
        // ONE sentence, from `pol_sf_reach_report` on EVERY branch -- the
        // caller-supplied-kernel one included -- so `meta["pol_sf_reach"]`,
        // this row and the banner cannot say three different things.
        pol_sf_reach_report(c, plan_),
        c.pol_sf == d.pol_sf);
  }
  {
    const bool no_rank2 = inclusive_rank2_is_empty(c);
    const bool incl = (c.channel == PipelineChannel::Inclusive) && !c.kernel;
    const bool b1_read = incl && !no_rank2;
    std::string val = c.kernel ? std::string("caller-supplied kernel")
                      : no_rank2 ? std::string(rank2_none_label())
                                 : std::string(b1_model_name(c.b1_model));
    // `no_rank2` (7Li inclusive) is REFUSED and not merely NotRead: the value
    // shown is `rank2_none_label()`, which names no backend and so cannot
    // mislead, AND `validate()` throws on cdks and li6-convolution there.
    // `Refused` is the status that says both.
    add("b1_model", "--b1-model", val,
        b1_read ? KnobStatus::Read : KnobStatus::Refused,
        b1_read ? rank2_input_report(c, plan_)
        : no_rank2 ? rank2_input_report(c, plan_)
                   : (c.kernel
                          ? std::string("the caller's kernel carries the "
                                        "rank-2 slots; validate() refuses any "
                                        "value but miller beside one")
                          : off_ch +
                                "the inclusive rank-2 slot reaches the "
                                "inclusive rate alone (a tagged channel's "
                                "alignment is in the event weight, the "
                                "coherent one's in the recoil azimuth), and "
                                "validate() refuses any value but miller "
                                "here"),
        c.b1_model == d.b1_model);
    const bool band_read = b1_read && c.b1_model != B1Model::Miller;
    add("b1_band_scale", "--b1-band-scale", fmt_g(c.b1_band_scale),
        band_read ? KnobStatus::Read : KnobStatus::Refused,
        band_read
            ? std::string("multiplies the WHOLE b1 of the opt-in backend -- "
                          "the MANDATORY 0/1/2 band of design_D_b1_li6.md "
                          "sec. 4.3")
            : std::string("the band is deliberately not applied to the "
                          "published Miller numbers and does not exist off "
                          "the inclusive 6Li rank-2 path, so validate() "
                          "refuses any value but 1.0 here rather than let "
                          "meta record a variation that did not run"),
        c.b1_band_scale == d.b1_band_scale);
    const bool conv = b1_read && c.b1_model == B1Model::Li6Convolution;
    add("b1_alpha_d_dwave_weight", "--b1-alpha-d-dwave-weight",
        fmt_g(c.b1_alpha_d_dwave_weight),
        conv ? KnobStatus::Read : KnobStatus::Refused,
        conv ? std::string("Li6ConvolutionOptions::w_alpha_d_dwave -- terms "
                           "(2d) and (2a) together, one physical effect")
             : std::string("only the li6-convolution branch reads it (miller "
                           "is a ratio model, cdks a digitized column), so "
                           "validate() refuses any value but 1.0 here"),
        c.b1_alpha_d_dwave_weight == d.b1_alpha_d_dwave_weight);
    add("b1_unpol", "--b1-unpol",
        c.kernel ? std::string("caller-supplied kernel")
        : no_rank2 ? std::string(rank2_none_label())
                   : std::string(b1_unpol_name(c.b1_unpol)),
        conv ? KnobStatus::Read : KnobStatus::Refused,
        conv ? std::string("the UNPOLARISED F1 the alpha-d convolution folds "
                           "against (Li6ConvolutionOptions::unpol); mstw is "
                           "the configuration the A = 2 gate passes on")
             : std::string("read ONLY by b1_model = li6-convolution, so "
                           "validate() refuses any value but toy here"),
        c.b1_unpol == d.b1_unpol);
    // `--r-source`, the ONE R hook into BOTH halves of the tensor weight.
    // Same reach as `b1_unpol` and for the same reason -- only the
    // `Li6Convolution` branch of `default_inclusive_kernel` installs it --
    // but a WIDER effect where it is read: it sets the kernel's own
    // `Options::r_func` as well, so it moves the unpolarised rate too, which
    // `b1_unpol` deliberately does not.
    // THREE STATUSES ON ONE AXIS, and the middle one is measured, not
    // reasoned: `sigma-lt` NAMES the value both null hooks already resolve
    // to (`resolve_r`, sf.hpp), so it installs the shared object and changes
    // nothing.  It is LABELLED rather than refused -- it names a member of
    // the R family, not a variation of a term that did not run -- which is
    // the `--pol-sf` precedent, and the criterion stated once at
    // `KnobProvenance`.
    const bool r_inert = conv && c.r_source == RSource::SigmaLt;
    add("r_source", "--r-source",
        c.kernel ? std::string("caller-supplied kernel")
        : no_rank2 ? std::string(rank2_none_label())
                   : std::string(r_source_name(c.r_source)),
        !conv ? KnobStatus::Refused
        : r_inert ? KnobStatus::NotRead
                  : KnobStatus::Read,
        r_inert
            ? std::string("not read at r_source = sigma-lt: the hook IS "
                          "installed in both halves, but resolve_r's null "
                          "branch (sf.hpp) IS r_sigma_lt, so naming it "
                          "reproduces the unset run exactly -- MEASURED "
                          "2026-09-06 bit-identical in sigma_pb, in all three "
                          "per-category cross sections and in every generated "
                          "column (6Li inclusive, li6-convolution, 2000 "
                          "events, seed 7, x_max 0.95), and in b1, K/D_phi, "
                          "A_zz and the cos 2phi amplitude at all six "
                          "standard points x = 0.05/0.10/0.30 x Q2 = 2.5/5 "
                          "for y = 0.1, 0.5 and 0.9")
        : conv ? std::string("ONE R = sigma_L/sigma_T into BOTH the alpha-d "
                             "convolution's F1 (Li6ConvolutionOptions::"
                             "r_func, the tensor weight's numerator) and this "
                             "kernel's F1/F_L/D(y)/ToyG1 (InclusiveKernel::"
                             "Options::r_func, its denominator), so the "
                             "ratio's two halves are the same choice; unset "
                             "installs neither and is today bit for bit")
               : std::string("read ONLY by b1_model = li6-convolution (miller "
                             "is a ratio model with no F1 of its own, cdks a "
                             "digitized column), so validate() refuses any "
                             "value but unset here"),
        c.r_source == d.r_source);
  }

  // ----------------------------------------------------- the tagged cluster
  {
    const bool tag = is_tagged(c.channel);
    const bool vmc = c.cluster_wave == ClusterWaveSource::VmcAV18;
    const bool li_tag = (c.channel == PipelineChannel::TaggedLi6Alpha ||
                         c.channel == PipelineChannel::TaggedLi7Alpha);
    add("cluster_wave", "--cluster-wave",
        cluster_wave_name(c.cluster_wave),
        tag ? KnobStatus::Read : KnobStatus::Refused,
        tag ? std::string("the family of radial forms of this channel's "
                          "cluster relative wave -- and, on the 6Li alpha "
                          "tag and the d control, of the embedded deuteron "
                          "too (open items C5.4, C5.5b)")
            : (off_ch +
               "the inclusive 6Li constants (LI6_CLUSTER_POLARIZATION = "
               "0.811228) and the coherent form factors do not follow it, so "
               "validate() refuses it here rather than promise a VMC row and "
               "deliver the shipped one, 11.61 % away (C5.5)"),
        c.cluster_wave == d.cluster_wave);
    const bool mc_read = vmc && li_tag;
    add("cluster_vmc_mc_sigma", "", fmt_g(c.cluster_vmc_mc_sigma),
        mc_read ? KnobStatus::Read : KnobStatus::Refused,
        mc_read ? std::string("the ANL VMC tables' own printed 1-sigma band, "
                              "fully correlated across k (open item C5.2); "
                              "+-1 sigma is 0.02 % on the tagged tensor "
                              "dilution")
                : std::string("only li6_ad1.momentum and li7_at3.momentum "
                              "print MC errors -- the analytic Hulthen forms "
                              "and fdeut.av18 print none -- so validate() "
                              "refuses any value but 0 here"),
        c.cluster_vmc_mc_sigma == d.cluster_vmc_mc_sigma);
    const bool beta_read =
        tag && (!vmc || (c.channel == PipelineChannel::TaggedLi7Alpha &&
                         tier_ == Tier::T1));
    add("cluster_beta", "--cluster-beta", fmt_g(c.cluster_beta),
        beta_read ? KnobStatus::Read : KnobStatus::NotRead,
        beta_read
            ? (vmc ? std::string("read on THIS run through the T1 breakup "
                                 "alone: cluster_wave = vmc takes the alpha-t "
                                 "RELATIVE motion from the ANL table, but the "
                                 "triton's own sequential two-body decay is "
                                 "still analytic at this beta (breakup.hpp)")
                   : std::string("the short-range scale of every analytic "
                                 "radial form this run builds -- the cluster "
                                 "relative wave and, at T1, the breakup "
                                 "(BreakupOptions::beta is overwritten with "
                                 "it so the draws cannot disagree)"))
            : (tag ? std::string("not read at cluster_wave = vmc on channel ") +
                         chn +
                         ": the ANL tables carry their own radial scale and "
                         "IGNORE this knob (cluster.hpp, docs/CONVENTIONS.md)"
                         + (c.channel == PipelineChannel::TaggedLi7Alpha
                                ? "; the 7Li T1 triton breakup would read it, "
                                  "but this run is at T0"
                                : "")
                   : off_ch +
                         "the cluster relative waves live on the tagged "
                         "channels only"),
        c.cluster_beta == d.cluster_beta);
    const bool pd_read =
        (c.channel == PipelineChannel::TaggedLi6Alpha) && !vmc;
    add("p_d", "--p-d", fmt_g(c.p_d),
        pd_read ? KnobStatus::Read : KnobStatus::NotRead,
        pd_read ? std::string("the alpha-d relative D-state probability of "
                              "the 6Li alpha tag's analytic wave")
        : (c.channel == PipelineChannel::TaggedLi6Alpha
               ? std::string("not read at cluster_wave = vmc: the ANL alpha-d "
                             "table carries its own D-state weight and "
                             "IGNORES this knob (docs/CONVENTIONS.md)")
               : off_ch +
                     "it is the 6Li alpha tag's alpha-d D state; the 7Li "
                     "alpha-t and d control waves and the deuteron's own "
                     "P_D (the scenario constant) are not this field"),
        c.p_d == d.p_d);
    const bool triton_read =
        (c.channel == PipelineChannel::TaggedLi7Alpha) && tier_ == Tier::T1;
    add("triton_sf", "--triton-sf", triton_sf_name(c.triton_sf),
        triton_read ? KnobStatus::Read : KnobStatus::NotRead,
        triton_read
            ? std::string("the spectral function the 7Li alpha tag's T1 "
                          "triton breakup draws from (triton_sf.hpp)")
        : (c.channel == PipelineChannel::TaggedLi7Alpha
               ? std::string("not read at Tier::T0: nothing inside the triton "
                             "is resolved, so no spectral function is "
                             "consulted")
               : off_ch +
                     "there is no triton to break up -- it is read on "
                     "tagged-7Li-alpha at T1 and nowhere else"),
        c.triton_sf == d.triton_sf);
    add("tier", "", tier_ == Tier::T1 ? "T1" : "T0",
        tag ? KnobStatus::Read : KnobStatus::NotRead,
        tag ? std::string("T1 resolves the struck cluster into a nucleon plus "
                          "its partner spectator(s); T0 leaves it one "
                          "off-shell pseudo-particle")
            : (off_ch +
               "there is no struck cluster to resolve, so Pipeline::tier() "
               "reports T0 whatever the configuration says"),
        c.tier == d.tier);
    const bool ib1_read = (c.channel == PipelineChannel::TaggedLi6Alpha);
    add("inclusive_b1", "--inclusive-b1", yn(c.struck.inclusive_b1),
        ib1_read ? KnobStatus::Read : KnobStatus::NotRead,
        ib1_read
            ? std::string("puts an inclusive b1 in the EMBEDDED DEUTERON's "
                          "kernel (StruckClusterOptions::inclusive_b1); it is "
                          "off by default because the alpha-d density is "
                          "already in the event weight")
        : (tag ? std::string("not read on channel ") + chn +
                     ": its struck cluster is spin 1/2 (dis_target = triton "
                     "resp. free neutron) and InclusiveKernel::tables opens "
                     "no rank-2 sector there, so the slot is never filled"
               : off_ch +
                     "StruckClusterOptions is the TAGGED struck-cluster "
                     "kernel's option block and no such kernel is built here"),
        c.struck.inclusive_b1 == d.struck.inclusive_b1);
  }

  // ------------------------------------------------------------------ FSI
  {
    const bool tag = is_tagged(c.channel);
    add("fsi", "--fsi", pipeline_fsi_name(c.fsi),
        tag ? KnobStatus::Read : KnobStatus::Refused,
        tag ? std::string("a per-event WEIGHT on Event::weight (never a shift "
                          "of any four-vector), spin independent by "
                          "construction")
            : (off_ch +
               "the weight is a distortion of the TAGGED spectator's "
               "spectrum and there is no tagged spectator here, so validate() "
               "refuses it"),
        c.fsi == d.fsi);
    const bool on = c.fsi != PipelineFsi::Off;
    add("fsi_sigma_mb", "--fsi-sigma-mb", fmt_g(c.fsi_sigma_mb),
        on ? KnobStatus::Read : KnobStatus::NotRead,
        on ? std::string("sigma_XN the Glauber weight is built at; the "
                         "documented band is 20-40 mb and a single row is "
                         "never a result (fsi.hpp)")
           : std::string("not read at fsi = off: no GlauberFsiWeight is "
                         "built, so there is nothing for sigma_XN to scale"),
        c.fsi_sigma_mb == d.fsi_sigma_mb);
  }

  // ------------------------------------------------------------------- RC
  {
    const RcModel* rc = rc_.get();
    const bool on = (c.rc != PipelineRc::Off);
    const bool band = on && rc && rc->applies();
    const bool tail = on && rc && rc->tail_applies();
    add("rc", "--rc", rc_mode_name(c.rc), KnobStatus::Read,
        on ? ("rc_tensor_lo / rc_tensor_hi / rc_tail are written as extra "
              "columns and go on Event::rc_weights, never on Event::weight; "
              "on this run the band " +
              std::string(band ? "APPLIES" : "does NOT apply") +
              " and the tail " +
              std::string(tail ? "APPLIES" : "does NOT apply") +
              (rc && !(band && tail)
                   ? std::string(" -- ") + rc->exclusion_reason().substr(
                                               0, 160) + " ..."
                   : std::string("")))
           : std::string("off: no rc_* column is written, Event::rc_weights "
                         "is empty and the run is today bit for bit"),
        c.rc == d.rc);
    // The two families, and which piece of the model each one scales.  A
    // sub-knob whose piece did not run is REFUSED, not labelled: it would
    // record a systematic as PRICED that was never computed -- the
    // `rc_qe_tensor_scale` precedent, generalised (see `KnobProvenance`).
    const RcOptions& o = c.rc_options;
    const RcOptions od;
    auto band_row = [&](const char* name, const char* flag, std::string value,
                        std::string what, bool at_def) {
      add(name, flag, std::move(value),
          band ? KnobStatus::Read
               : (on ? KnobStatus::Refused : KnobStatus::NotRead),
          band ? ("BAND knob: " + what)
          : on  ? ("the band does not apply on channel " + chn +
                  " (RcModel::applies() is false, every rc_tensor_* weight is "
                  "exactly 1.0), so this knob would record a systematic as "
                  "PRICED that was never computed -- validate() refuses any "
                  "value but the default here")
                : std::string("not read at rc = off: no RcModel is built"),
          at_def);
    };
    auto tail_row = [&](const char* name, const char* flag, std::string value,
                        std::string what, bool at_def) {
      add(name, flag, std::move(value),
          tail ? KnobStatus::Read
               : (on ? KnobStatus::Refused : KnobStatus::NotRead),
          tail ? ("TAIL knob: " + what)
          : on  ? ("the radiative tail does not apply on this run (rc_tail == "
                  "1 exactly" +
                  std::string(is_tagged(c.channel)
                                  ? " by construction on channel " + chn
                                  : "") +
                  "), so this knob would record a systematic as PRICED that "
                  "was never computed -- validate() refuses any value but the "
                  "default here")
                : std::string("not read at rc = off: no RcModel is built"),
          at_def);
    };
    band_row("rc_delta_low_x", "--rc-delta-low-x", fmt_g(o.delta_low_x),
             "the low-x band edge; 0.30 is the conservative end of "
             "Gakh-Shekhovtsova's uncited 10-30 %, 0.19 the residual HERMES "
             "achieved -- and 0.113 is what that source's panel READS at "
             "x = 0.00966, the x nearest this anchor's 0.01, so the default "
             "errs WIDE by x2.65 (priced 2026-09-06: half-width on A_zz at "
             "x = 0.01, Q2 = 5 is 1.358e-04 / 1.205e-04 / 5.117e-05 at "
             "0.30 / 0.266 / 0.113, unmoved at x >= 0.16)",
             o.delta_low_x == od.delta_low_x);
    band_row("rc_delta_high_x", "--rc-delta-high-x", fmt_g(o.delta_high_x),
             "the high-x band edge (E12-13-011); the record names NO "
             "alternative value to band it against -- the 1.5 % is in the "
             "unpublished proposal only and the published companion has no "
             "RC discussion, so the only alternative stated anywhere is "
             "'cite it by page or drop the anchor' (checked 2026-09-06)",
             o.delta_high_x == od.delta_high_x);
    band_row("rc_x_low", "", fmt_g(o.x_low), "the low-x band anchor",
             o.x_low == od.x_low);
    band_row("rc_x_high", "", fmt_g(o.x_high), "the high-x band anchor",
             o.x_high == od.x_high);
    band_row("rc_a_transfer_frac", "--rc-a-transfer-frac",
             fmt_g(o.a_transfer_frac),
             "prices the A = 2 -> A = 6 transfer of delta(x) in QUADRATURE, "
             "i.e. widens the band by sqrt(1 + f^2)",
             o.a_transfer_frac == od.a_transfer_frac);
    band_row("rc_band_tau_max", "", fmt_g(o.band_tau_max),
             "the ceiling on |tau| the band sees", o.band_tau_max == od.band_tau_max);
    // rc_scope is the ONE band knob with a second condition, and it is the
    // fill's own AXIS.  `TensorAll` adds the band to the cos 2phi (delta)
    // amplitude as well as to the rate, and `tensor_amplitudes` builds that
    // amplitude with a sin^2(theta_S) factor -- so at theta_S = 0, which is
    // every `tensor_thirds_plan` category, `tensor-all` and `tensor-rate` are
    // the same run.  MEASURED 2026-09-05 (400 events, seed 7, inclusive 6Li
    // --rc tensor-band): bit-identical under tensor-thirds, and it MOVES
    // under transverse-tensor and tensor-flip, whose categories sit at
    // theta_S = pi/2.
    bool tilted = false;
    for (const SpinCategory& sc : plan_.categories())
      if (std::fabs(std::sin(sc.theta_s)) > 1e-12) tilted = true;
    const bool scope_read = band && tilted;
    add("rc_scope", "",
        o.scope == RcScope::TensorRate ? "tensor-rate" : "tensor-all",
        scope_read ? KnobStatus::Read
                   : (band ? KnobStatus::NotRead
                           : (on ? KnobStatus::Refused
                                 : KnobStatus::NotRead)),
        scope_read
            ? std::string("BAND knob: tensor-all adds the band to the "
                          "cos 2phi (delta) amplitude as well as to the rate")
        : band ? (by_plan +
                  "tensor-all differs from tensor-rate only in the cos 2phi "
                  "amplitude, which tensor_amplitudes builds with a "
                  "sin^2(theta_S) factor, and every category of this fill "
                  "sits at theta_S = 0 -- measured bit-identical under "
                  "tensor-thirds, and it does move under transverse-tensor "
                  "and tensor-flip")
        : on   ? std::string("the band does not apply on channel " + chn +
                            " (RcModel::applies() is false, every rc_tensor_* "
                            "weight is exactly 1.0), so this knob would "
                            "record a systematic as PRICED that was never "
                            "computed -- validate() refuses any value but the "
                            "default here")
               : std::string("not read at rc = off: no RcModel is built"),
        o.scope == od.scope);
    tail_row("rc_fq_scale", "--rc-fq-scale", fmt_g(o.fq_scale),
             "+-100 % systematic on the 6Li quadrupole form factor; the "
             "tensor tail is QUADRATIC in it", o.fq_scale == od.fq_scale);
    tail_row("rc_tail_tensor_scale", "--rc-tail-tensor-scale",
             fmt_g(o.tail_tensor_scale),
             "multiplier on the 6Li MAGNETIC form factor -- the eta F_m^2 "
             "tensor sector --rc-fq-scale does not span",
             o.tail_tensor_scale == od.tail_tensor_scale);
    tail_row("rc_c0_shape", "--rc-c0-shape", c0_shape_name(o.c0_shape),
             "which 6Li monopole shape F_c AND F_q share; a SHAPE, and its "
             "sign at x = 0.1 is a band edge", o.c0_shape == od.c0_shape);
    // THE QUASI-ELASTIC SUB-FAMILY, one level below the tail.  A `tail_row`
    // asks whether the RADIATIVE TAIL ran; these three additionally need
    // POLRAD Eq. (44)'s QUASI-ELASTIC piece to have run, and
    // `with_qe_tail = false` switches exactly that piece off while leaving
    // the elastic t-peak in place.  Until 2026-09-05 they were `tail_row`s,
    // so on an inclusive `--rc tensor-band` run with `with_qe_tail = false`
    // `rc_qe_suppression = 0.5` and `rc_qe_kf_gev = 0.25` were accepted,
    // bit-identical, and recorded read / non-default -- the criterion at
    // `KnobProvenance` (a scale on a term computed as identically 1 is
    // REFUSED) applied to `qe_tensor_scale` alone and not to the two knobs
    // beside it.  `validate()` and `RcModel`'s constructor now refuse all
    // three there, and this row says so.
    auto qe_row = [&](const char* name, const char* flag, std::string value,
                      std::string what, bool at_def) {
      add(name, flag, std::move(value),
          (tail && o.with_qe_tail)
              ? KnobStatus::Read
              : (on ? KnobStatus::Refused : KnobStatus::NotRead),
          (tail && o.with_qe_tail) ? ("QUASI-ELASTIC TAIL knob: " + what)
          : !on ? std::string("not read at rc = off: no RcModel is built")
          : !tail
              ? ("the radiative tail does not apply on this run (rc_tail == 1 "
                 "exactly" +
                 std::string(is_tagged(c.channel)
                                 ? " by construction on channel " + chn
                                 : "") +
                 "), so this knob would record a systematic as PRICED that "
                 "was never computed -- validate() refuses any value but the "
                 "default here")
              : std::string("rc_options.with_qe_tail = false switches POLRAD "
                            "Eq. (44)'s unpolarised quasi-elastic tail off, "
                            "so sigma^q_U is not computed at all and this "
                            "knob scales nothing -- validate() refuses any "
                            "value but the default here rather than record a "
                            "systematic as PRICED that was never computed"),
          at_def);
    };
    qe_row("rc_qe_suppression", "--rc-qe-suppression",
           fmt_g(o.qe_suppression),
           "flat multiplier on the unpolarised quasi-elastic tail, on top "
           "of the de Forest-Walecka Pauli suppression",
           o.qe_suppression == od.qe_suppression);
    qe_row("rc_qe_tensor_scale", "--rc-qe-tensor-scale",
           fmt_g(o.qe_tensor_scale),
           "prices the POLARISED quasi-elastic tail with a BORROWED "
           "magnitude; 0 is the shipped tensor-blind tail",
           o.qe_tensor_scale == od.qe_tensor_scale);
    qe_row("rc_qe_kf_gev", "", fmt_g(o.qe_kf_gev),
           "6Li's measured Fermi momentum in the Pauli suppression S(q)",
           o.qe_kf_gev == od.qe_kf_gev);
    tail_row("rc_tail_model", "--rc-tail-model", rc_tail_model_name(o.tail_model),
             "t-peak alone (a LOWER bound), t-peak plus the leading-log s- "
             "and p-peaks, or POLRAD Eq. (18)'s exact tau_A quadrature",
             o.tail_model == od.tail_model);
    // THE s-/p-PEAK SUB-FAMILY, the same shape as `qe_row` one branch over.
    // `sp_tensor_scale` additionally needs the LEADING-LOG s-/p-peaks to have
    // run, and only `tail_model = TPeakPlusLL` computes them; under the
    // shipped `TPeak` the `u_sp` table is identically zero, so the scale is a
    // multiplier on a term this run computes as exactly 0 -- REFUSED by
    // `validate()` and by `RcModel`'s constructor, not labelled, under the
    // criterion stated once at `KnobProvenance`.
    {
      const bool sp = tail && o.tail_model == RcTailModel::TPeakPlusLL;
      const bool pf = o.tail_model == RcTailModel::PolradFull;
      add("rc_sp_tensor_scale", "--rc-sp-tensor-scale",
          fmt_g(o.sp_tensor_scale),
          sp ? KnobStatus::Read
             : (on ? KnobStatus::Refused : KnobStatus::NotRead),
          sp ? std::string("LEADING-LOG s-/p-PEAK knob: prices the tensor "
                           "fraction of the s-/p-peaks with a BORROWED "
                           "magnitude -- a BOUND WITH NO DERIVATION; 0 is the "
                           "shipped tensor-blind s+p")
          : !on ? std::string("not read at rc = off: no RcModel is built")
          : !tail
              ? ("the radiative tail does not apply on this run (rc_tail == 1 "
                 "exactly" +
                 std::string(is_tagged(c.channel)
                                 ? " by construction on channel " + chn
                                 : "") +
                 "), so this knob would record a systematic as PRICED that "
                 "was never computed -- validate() refuses any value but the "
                 "default here")
          : pf ? std::string("rc_options.tail_model = polrad-full needs no "
                             "s/p tensor stand-in: Eq. (18) carries the s- "
                             "and p-peaks inside its own tau_A integral WITH "
                             "their Eq. (A.4) tensor content, so this knob "
                             "would double-count a term that RAN -- "
                             "validate() refuses any value but the default "
                             "here.  The OPPOSITE reason to the t-peak one")
              : std::string("rc_options.tail_model = t-peak computes no "
                            "leading-log s-/p-peaks at all (the u_sp table is "
                            "identically zero), so this knob scales nothing "
                            "-- validate() refuses any value but the default "
                            "here rather than record a systematic as PRICED "
                            "that was never computed"),
          o.sp_tensor_scale == od.sp_tensor_scale);
    }
    tail_row("rc_with_qe_tail", "", yn(o.with_qe_tail),
             "whether POLRAD Eq. (44)'s unpolarised quasi-elastic tail is "
             "computed at all", o.with_qe_tail == od.with_qe_tail);
    // `rc_with_tail` is NOT a `tail_row`, and could not be: it is the knob
    // that DECIDES whether the tail runs, so refusing it "because the tail
    // did not run" would be circular, and a `Refused` row would be a false
    // claim -- `with_tail = false` is exactly the value that makes the tail
    // not apply, and it is accepted.  Its own rule: it is read wherever the
    // tail WOULD apply if it were true, i.e. `rc_tail_applies(channel,
    // true)`, and labelled where the channel already forces rc_tail == 1.
    {
      const bool would = on && rc_tail_applies(event_channel_of(c.channel),
                                               /*with_tail=*/true);
      add("rc_with_tail", "", yn(o.with_tail),
          would ? KnobStatus::Read : KnobStatus::NotRead,
          would ? std::string("whether the radiative tail is computed at all "
                              "(false = a band-only diagnostic run)")
          : on  ? (off_ch +
                   "rc_tail == 1 exactly here whatever this knob says, so "
                   "switching the tail off changes nothing -- it is the "
                   "CHANNEL that already did")
                : std::string("not read at rc = off: no RcModel is built"),
          o.with_tail == od.with_tail);
    }
    tail_row("rc_n_eta", "", std::to_string(o.n_eta),
             o.tail_model == RcTailModel::PolradFull
                 ? "tanh-sinh nodes PER PANEL of Eq. (18)'s tau_A integral "
                   "(four panels, split at tau_s, 0 and tau_p)"
                 : "quadrature nodes of the eta_A integral of the tail",
             o.n_eta == od.n_eta);
    tail_row("rc_tail_max", "", fmt_g(o.tail_max),
             "ceiling on the returned tail ratio; the clipped fraction is "
             "reported, not hidden", o.tail_max == od.tail_max);
    // `rc_m_lepton` WAS `Refused` UNCONDITIONALLY UNTIL 2026-09-06, with the
    // reason "RESERVED for tail_model = PolradFull".  PolradFull exists now
    // and reads it -- it is the m^2 of C_{1,2}(tau), of F_IR and of
    // lambda_s, i.e. what regulates the s- and p-peaks -- so the row is
    // three-way, exactly like `rc_sp_tensor_scale` one branch up.
    {
      const bool ml = tail && o.tail_model == RcTailModel::PolradFull;
      add("rc_m_lepton", "", fmt_g(o.m_lepton),
          ml ? KnobStatus::Read : (on ? KnobStatus::Refused
                                      : KnobStatus::NotRead),
          ml ? std::string("the LEPTON MASS of POLRAD Eq. (18): the m^2 of "
                           "C_{1,2}(tau) (Eq. (B.13)), of F_IR = m^2 F_2+ - "
                           "Q_m^2 F_d and of lambda_s -- it is what regulates "
                           "the s- and p-peaks")
          : !on ? std::string("not read at rc = off: no RcModel is built")
          : !tail
              ? (off_ch +
                 "the radiative tail does not apply on this run, so no "
                 "quadrature reads a lepton mass at all")
              : std::string("read ONLY by tail_model = polrad-full.  The "
                            "t-peak forms of Eqs. (37)-(39) carry no lepton "
                            "mass and TPeakPlusLL's leading-log radiator "
                            "takes M_ELECTRON from constants.hpp directly, "
                            "so validate() refuses any value but M_ELECTRON "
                            "here"),
          o.m_lepton == od.m_lepton);
    }
  }

  // ------------------------------------------------------ coherent channel
  {
    const bool coh = (c.channel == PipelineChannel::CoherentLi6);
    const CoherentScenario& s = c.coherent;
    const CoherentScenario sd;
    auto coh_row = [&](const char* name, const char* flag, std::string value,
                       std::string what, bool at_def) {
      add(name, flag, std::move(value),
          coh ? KnobStatus::Read : KnobStatus::NotRead,
          coh ? what
              : (off_ch +
                 "the coherent scenario is built on CoherentLi6 alone; no "
                 "other channel has an intact recoil"),
          at_def);
    };
    coh_row("coherent_t_max", "--coherent-t-max", fmt_g(c.coherent_t_max),
            "the |t| ceiling: it moves the whole |t| spectrum, the tag "
            "acceptance and every c_2 in the file",
            c.coherent_t_max == d.coherent_t_max);
    coh_row("coherent_f0", "--coherent-f0", fmt_g(s.f0),
            "the coherent fraction's normalisation", s.f0 == sd.f0);
    coh_row("coherent_slope_b", "--coherent-slope-b", fmt_g(s.slope_b),
            "the |F(t)|^2 t-slope B [GeV^-2]", s.slope_b == sd.slope_b);
    coh_row("coherent_amp", "--coherent-amp", fmt_g(s.amp),
            "the tensor amplitude of the recoil azimuth", s.amp == sd.amp);
    coh_row("coherent_eps_b0", "", fmt_g(s.eps_b0),
            "the deformation input the positivity edge is derived from",
            s.eps_b0 == sd.eps_b0);
    coh_row("coherent_m_x_min", "", fmt_g(c.coherent_xpom.m_x_min),
            "the smallest diffractive mass; it is the gate that decides "
            "which cells carry coherent rate at all",
            c.coherent_xpom.m_x_min == d.coherent_xpom.m_x_min);
    coh_row("coherent_x_pom_max", "", fmt_g(c.coherent_xpom.x_pom_max),
            "the upper x_P edge of the diffractive region",
            c.coherent_xpom.x_pom_max == d.coherent_xpom.x_pom_max);
    coh_row("coherent_weighted_azimuth", "",
            yn(c.coherent_weighted_azimuth),
            "draw phi_t from the modulated density instead of carrying it as "
            "an event weight",
            c.coherent_weighted_azimuth == d.coherent_weighted_azimuth);
  }

  // -------------------------------------------------------------- T2 tier
  {
    const bool coh = (c.channel == PipelineChannel::CoherentLi6);
    // WHETHER A T2 TIER IS BOUND IS SOMETHING THE CORE CAN SEE.  It is
    // `PipelineConfig::hadronizer`, a `std::function` on the config, and this
    // row used to key on `ctx.t2_bound` alone -- the CONTEXT flag, which only
    // a caller who went through `set_pythia_hadronizer` (or the CLI) sets.  A
    // hadronizer bound as a plain callable -- `cfg.hadronizer = f`, the route
    // python/README.md:246 documents, and the C++ lambda of USAGE.md:1741 --
    // was then called on every event while this row said "off" and "no T2
    // tier is bound: every record stops at T0".  What the context adds is not
    // WHETHER a hook is bound but WHICH bridge it is: the core deliberately
    // does not link the PYTHIA tier (sf.hpp's rule), so it cannot ask a
    // type-erased hook for `PDF:PomSet`.  Two facts, two sources, and the
    // three Pomeron rows below are careful to keep them apart.
    const bool bound = static_cast<bool>(c.hadronizer);
    const bool named_bridge = bound && ctx.t2_bound;
    add("hadronize", "--hadronize", bound ? "on" : "off",
        KnobStatus::Read,
        named_bridge
            ? std::string("a T2 hadronizer is bound to "
                          "PipelineConfig::hadronizer and runs on every event")
        : bound
            ? std::string("a hadronizer is bound to "
                          "PipelineConfig::hadronizer and runs on every event "
                          "-- the core sees the std::function, and the "
                          "KnobRunContext does NOT name a PythiaBridge behind "
                          "it, so the three Pomeron rows below describe "
                          "nothing this run built")
            : std::string("no T2 tier is bound: every record stops at T0 "
                          "(and T1 on the tagged channels)"),
        !bound);
    // THE MEASUREMENT BEHIND THESE THREE.  The Pomeron PYTHIA instance is
    // built whenever coherent_t2 == Pomeron REGARDLESS of channel
    // (pythia_bridge.cpp) and then hadronizes ZERO events off the coherent
    // channel.  Measured 2026-09-05 with a deterministic hadron hash, 150
    // events seed 4242: pom_set 6 vs 5 vs coherent_t2 off moves NOTHING --
    // not a T0 column, not one of 3278 particles -- on inclusive 6Li/7Li/d,
    // tagged-6Li-alpha, tagged-7Li-alpha and tagged-d-p, while meta recorded
    // pom_set = 5 and coherent_t2 = "pomeron" as if they had run.
    // ... AND THE "not read without --hadronize" SENTENCE IS NOW CONSISTENT
    // WITH THE ROW ABOVE.  It keys on `bound`, the same fact the `hadronize`
    // row states, so a run with a callable hadronizer can no longer be told
    // "no PYTHIA instance of any kind is built" by one row while the row
    // above it says a hadronizer runs on every event.  The middle case --
    // bound, but not through `set_pythia_hadronizer` -- gets its own
    // sentence, because there the values shown are `PythiaBridgeOptions`'
    // defaults and describe no object this run holds.
    const std::string t2_off =
        named_bridge
            ? (off_ch +
               "the Pomeron PYTHIA instance is built whenever coherent_t2 = "
               "pomeron regardless of channel, and it hadronizes ZERO events "
               "here (measured: pom_set 6 vs 5 vs coherent_t2 off moves "
               "neither a T0 column nor one hadron off the coherent channel)")
        : bound
            ? std::string("not read through this run's hadronizer: it is a "
                          "callable bound to PipelineConfig::hadronizer and "
                          "not a PythiaBridge attached by "
                          "set_pythia_hadronizer, so the core cannot ask it "
                          "for PDF:PomSet and the value shown is "
                          "PythiaBridgeOptions' own default, not something "
                          "this run configured")
            : std::string("not read without --hadronize: no PYTHIA instance "
                          "of any kind is built, so there is no Pomeron beam "
                          "for it to configure");
    const bool t2_read = named_bridge && coh;
    add("coherent_t2", "--coherent-t2", ctx.t2_pomeron ? "pomeron" : "off",
        t2_read ? KnobStatus::Read : KnobStatus::NotRead,
        t2_read ? std::string("what --hadronize does with a coherent event: "
                              "hadronize the gamma*-Pomeron system on a "
                              "PYTHIA Pomeron beam (id 990), or leave the "
                              "record at T0")
                : t2_off,
        ctx.t2_pomeron);
    const bool pom_read = t2_read && ctx.t2_pomeron;
    add("pom_set", "--pom-set", std::to_string(ctx.pom_set),
        pom_read ? KnobStatus::Read : KnobStatus::NotRead,
        pom_read ? std::string("PDF:PomSet, the coherent T2 tier's LARGEST "
                               "model systematic: it moves nothing at T0 and "
                               "the kaon fraction by a factor 2.8 over the 12 "
                               "DPDF fits about set 6")
                 : (t2_read ? std::string("not read at coherent_t2 = off: no "
                                          "Pomeron instance is built")
                            : t2_off),
        ctx.pom_set == 6);
    add("pom_rescale", "--pom-rescale", fmt_g(ctx.pom_rescale),
        pom_read ? KnobStatus::Read : KnobStatus::NotRead,
        pom_read ? std::string("PDF:PomRescale, the overall Pomeron-PDF "
                               "normalisation; it cancels out of the bridge's "
                               "per-event flavour draw and is recorded for "
                               "reproducibility")
                 : (t2_read ? std::string("not read at coherent_t2 = off: no "
                                          "Pomeron instance is built")
                            : t2_off),
        ctx.pom_rescale == 1.0);
  }
  return rows;
}

double Pipeline::optics_lumi_factor() const {
  if (!cfg_.apply_optics_lumi_fraction) return 1.0;
  const double f = optics_.lumi_fraction;
  return (f > 0.0 && f <= 1.0) ? f : 1.0;
}

double Pipeline::sigma_pb() const {
  double s = 0.0;
  for (std::size_t k = 0; k < plan_.categories().size(); ++k) {
    s += plan_.categories()[k].lumi_fraction * sigma_[k];
  }
  return s;
}

std::size_t Pipeline::category_of(std::uint64_t index) const {
  if (index >= total_) throw std::runtime_error("Pipeline: index out of range");
  const std::size_t k = static_cast<std::size_t>(
      std::upper_bound(offset_.begin(), offset_.end(), index) - offset_.begin());
  return k - 1;
}

void Pipeline::label_event(Event& ev, std::size_t k, std::uint64_t index) const {
  const SpinCategory& cat = plan_.categories()[k];
  ev.number = index;
  ev.spin.j = cat.j;
  ev.spin.lam_e = cat.lam_e;
  ev.spin.pe = cat.pe;
  ev.spin.theta_s = cat.theta_s;
  ev.spin.phi_s = cat.phi_s;
  ev.spin.pz = plan_.pz_true();
  ev.spin.pzz = plan_.pzz_true();
  ev.spin.category = cat.name;
  ev.spin.run = static_cast<int>(cfg_.run);
  ev.spin.bunch = static_cast<int>(k);
  ev.xsec_pb = sigma_[k];
}

void Pipeline::add_beams(Event& ev, const SpinCategory& cat,
                         double m_ion) const {
  Particle be;
  be.pdg = 11;
  be.status = Status::Beam;
  be.role = Role::BeamElectron;
  be.p = beam_e_;
  be.mass = 0.0;
  be.charge = -1.0;
  be.pol = cat.lam_e != 0 ? static_cast<double>(cat.lam_e) : 9.0;

  Particle bi;
  bi.pdg = ion_pdg_;
  bi.status = Status::Beam;
  bi.role = Role::BeamIon;
  bi.p = beam_ion_;
  bi.mass = ion_mass_;
  bi.charge = ion_charge_;
  bi.pol = m_ion;

  ev.particles.push_back(be);
  ev.particles.push_back(bi);
}

void Pipeline::add_hadronic_x(Event& ev, const Vec4& p_x, double charge,
                              int mother) const {
  // C1.  X is the whole hadronic final state and MUST be timelike.  This used
  // to be sqrt(max(m2, 0)), which silently turned a spacelike residual into a
  // massless X -- and the coherent channel produced one in 100 % of events
  // while x_P was pinned at zero.  A hard check, not a clip: it is an
  // invariant of the channel's own balance, not a condition on the event.
  const double m2 = p_x.m2();
  if (!(m2 >= 0.0)) {
    char buf[224];
    std::snprintf(buf, sizeof(buf),
                  "Pipeline: the hadronic system X came out SPACELIKE on the "
                  "%s channel, M_X^2 = %.6g GeV^2 (x = %.6g, Q2 = %.6g)",
                  pipeline_channel_name(cfg_.channel), m2, ev.kin.x, ev.kin.q2);
    throw std::runtime_error(buf);
  }
  Particle x;
  x.pdg = 92;
  x.status = Status::Final;
  x.role = Role::HadronicX;
  x.p = p_x;
  x.mass = std::sqrt(m2);
  x.charge = charge;
  x.mother1 = mother;
  x.mother2 = find_role(ev, Role::VirtualPhoton);
  ev.particles.push_back(x);
}

void Pipeline::make_inclusive(std::size_t k, std::uint64_t local,
                              std::uint64_t index, Rng& rng, Event& ev) const {
  (void)local;
  const SpinCategory& cat = plan_.categories()[k];
  const EventDraw draw = dis_sampler_->draw_event(cplan_[k], rng);
  gen_->make_event(cat, plan_, draw, rng, index, static_cast<int>(cfg_.run),
                   static_cast<int>(k), ev);
  ev.xsec_pb = sigma_[k];
}

void Pipeline::make_tagged(std::size_t k, std::uint64_t local,
                           std::uint64_t index, Rng& rng, Event& ev) const {
  (void)local;
  const SpinCategory& cat = plan_.categories()[k];

  // C1, and a defect the hard M_X^2 check below exposed.  The tagged struck
  // cluster is DELIBERATELY off shell (P_X = P_ion - p_spec, the impulse
  // approximation), and at a spectator momentum approaching the model's own
  // k_max = 1.2 GeV it goes so far off shell that X = k + P_X - k' comes out
  // SPACELIKE.  Measured on 200k events at the mid configuration: 0.020 % of
  // the d(e,e'p) control (spectator proton, m_spec = 0.94 against the
  // deuteron's 1.88 -- at k = 1.2 the struck neutron's P_X^2 is -1.3 GeV^2),
  // and 0 % of the 6Li and 7Li alpha tags, where the spectator is an alpha and
  // P_X^2 stays positive over the whole grid.  Before the check it was
  // silently clipped to a MASSLESS X, which is not a hadronic system.
  //
  // The draw is REJECTED and repeated on the event's own stream: the tagged
  // phase space is restricted to the region where the impulse approximation
  // closes, which is a statement about the model, not about the event.  The
  // event stays a pure function of its index, so determinism is untouched.
  // The test is on the FINISHED record, so the fast path -- 99.98 % of draws --
  // pays nothing for it; the rebuild on the rare rejection costs one more pass
  // through the same code.
  TaggedEvent te = tsampler_->sample_one(fills_[k], rate_cdf_[k], rng);
  for (int tries = 0;; ++tries) {
    if (tries > 0) {
      if (tries >= kTaggedMaxRedraw) {
        throw std::runtime_error(
            "Pipeline: no timelike hadronic system on the tagged channel "
            "after " + std::to_string(kTaggedMaxRedraw) + " draws");
      }
      te = tsampler_->sample_one(fills_[k], rate_cdf_[k], rng);
    }
  ev.reset();
  label_event(ev, k, index);
  ev.kin.x = te.x;
  ev.kin.q2 = te.q2;
  ev.kin.y = te.y;
  ev.kin.phi = te.phi;
  ev.kin.cell = te.cell;
  ev.kin.s = dis_sampler_->s();
  ev.kin.w2 = w2_from_xq2(te.x, te.q2);
  ev.kin.nu = te.q2 / (2.0 * M_NUCLEON * te.x);

  add_beams(ev, cat, te.m_ion);

  Particle esc;
  esc.pdg = 11;
  esc.status = Status::Final;
  esc.role = Role::ScatteredElectron;
  esc.p = dis_source_->scattered_electron_p4(te.x, te.y, te.phi);
  esc.mass = 0.0;
  esc.charge = -1.0;
  esc.mother1 = 0;
  esc.pol = ev.particles[0].pol;
  ev.particles.push_back(esc);

  if (cfg_.with_virtual_photon) {
    Particle g;
    g.pdg = 22;
    g.status = Status::Intermediate;
    g.role = Role::VirtualPhoton;
    g.p = beam_e_ - esc.p;
    g.mass = -std::sqrt(std::max(-g.p.m2(), 0.0));
    g.charge = 0.0;
    g.mother1 = 0;
    ev.particles.push_back(g);
  }

  // Adds the ON-SHELL spectator and the OFF-SHELL struck cluster
  // P_X = P_ion - p_spec, the tagging block of `kin` and the struck-cluster
  // spin label; also sets `Event::channel`.
  tsampler_->fill_event(ev, te);

  // FSI as a per-event WEIGHT and nothing else (fsi.hpp) -- the exact mirror
  // of the coherent channel's `ev.weight *= ce.weight`.  1.0 when
  // `PipelineConfig::fsi` is Off.  Inside the redraw loop is right:
  // `ev.reset()` restores weight = 1, so a redrawn event carries its OWN
  // spectator's weight, never a stale one.
  ev.weight *= te.weight;

  const int i_spec = find_role(ev, Role::Spectator);
  const int i_clus = find_role(ev, Role::StruckCluster);
  const Particle& spec = ev.particles[static_cast<std::size_t>(i_spec)];

  if (tier_ == Tier::T1) {
    // ---- T1: resolve the struck cluster ---------------------------------
    //
    // The partner spectators go on shell and the struck NUCLEON takes the
    // remainder of P_X, so the whole-nucleus balance becomes
    //     k + P_ion = k' + p_spec + sum(p_partner) + X,
    // and X is the PER-NUCLEON remainder X = k + p_N,struck - k'.  Naming
    // the struck nucleon is also what lets `PythiaBridge` find a target it
    // can use verbatim: its cluster branch never fires on a T1 event.
    BreakupInput bi;
    bi.species = cluster_species_;
    bi.p_cluster = ev.particles[static_cast<std::size_t>(i_clus)].p;
    bi.z = channel_->base.partner_Z();
    bi.a = channel_->base.partner_A();
    bi.m_s = te.m_struck;
    bi.s_cluster = channel_->s_channel;
    bi.x = te.x;
    bi.q2 = te.q2;
    bi.eff_pol_p = channel_->dis_target.eff_pol_p;
    bi.eff_pol_n = channel_->dis_target.eff_pol_n;
    BreakupResult br;
    if (!breakup_->resolve(bi, rng, br)) continue;

    const Vec4 p_x = (beam_e_ + beam_ion_) - esc.p - spec.p - sum_p(br.partners);
    // The T1 X is SMALLER than the T0 one by every partner spectator, so the
    // timelike test has to be redone on it -- a draw that was fine at T0 can
    // leave a spacelike per-nucleon remainder here.
    if (!(p_x.m2() >= 0.0)) continue;

    for (Particle& f : br.partners) {
      f.mother1 = i_clus;
      ev.particles.push_back(f);
    }
    br.struck.mother1 = i_clus;
    ev.particles.push_back(br.struck);
    const int i_n = static_cast<int>(ev.particles.size()) - 1;
    add_hadronic_x(ev, p_x, br.struck.charge, i_n);
    return;
  }

  const Particle& clus = ev.particles[static_cast<std::size_t>(i_clus)];
  // WHOLE-NUCLEUS balance: X = k + P_ion - k' - p_spec = k + P_X - k'.
  const Vec4 p_x = (beam_e_ + beam_ion_) - esc.p - spec.p;
  if (!(p_x.m2() >= 0.0)) continue;   // the rejection described above
  add_hadronic_x(ev, p_x, clus.charge, i_clus);
  return;
  }
}

void Pipeline::make_coherent(std::size_t k, std::uint64_t local,
                             std::uint64_t index, Rng& rng, Event& ev) const {
  (void)local;
  const SpinCategory& cat = plan_.categories()[k];

  // 1. the ion spin projection (label only: the coherent rate is spin
  //    independent and the tensor modulation is an ENSEMBLE coefficient
  //    carried by the recoil azimuth -- `CoherentScenario::a2_m_state` is the
  //    per-m relation, kept analytic).
  const std::vector<double>& ms = coh_ms_[k];
  double acc = 0.0;
  const double um = rng.uniform();
  double m_ion = ms.back();
  for (std::size_t i = 0; i < ms.size(); ++i) {
    acc += cat.populations[i];
    if (um < acc) { m_ion = ms[i]; break; }
  }
  // 2. (x, Q2) from the f_coh-reweighted cell CDF, then phi uniform.
  const std::size_t c = pick(coh_cdf_, rng.uniform());
  double x = 0.0, q2 = 0.0;
  draw_in_cell(*dis_sampler_, c, rng, x, q2);
  const double phi = 2.0 * kPi * rng.uniform();
  const double y = y_from_xq2(x, q2, dis_sampler_->s());
  ev.reset();
  label_event(ev, k, index);
  ev.spin.m_ion = m_ion;
  ev.spin.m_struck = kNaN;
  ev.kin.x = x;
  ev.kin.q2 = q2;
  ev.kin.y = y;
  ev.kin.phi = phi;
  ev.kin.cell = static_cast<int>(c);
  ev.kin.s = dis_sampler_->s();
  ev.kin.w2 = w2_from_xq2(x, q2);
  ev.kin.nu = q2 / (2.0 * M_NUCLEON * x);

  add_beams(ev, cat, m_ion);

  Particle esc;
  esc.pdg = 11;
  esc.status = Status::Final;
  esc.role = Role::ScatteredElectron;
  esc.p = gen_scattered_electron(x, y, phi);
  esc.mass = 0.0;
  esc.charge = -1.0;
  esc.mother1 = 0;
  esc.pol = ev.particles[0].pol;
  ev.particles.push_back(esc);

  if (cfg_.with_virtual_photon) {
    Particle g;
    g.pdg = 22;
    g.status = Status::Intermediate;
    g.role = Role::VirtualPhoton;
    g.p = beam_e_ - esc.p;
    g.mass = -std::sqrt(std::max(-g.p.m2(), 0.0));
    g.charge = 0.0;
    g.mother1 = 0;
    ev.particles.push_back(g);
  }

  // 3. the recoil.  It is solved AGAINST this event's own residual
  // R = k + P_ion - k', so that X = R - P_recoil comes out at exactly the
  // diffractive mass the drawn x_P implies (C1).
  CoherentDis dis;
  dis.x = x;
  dis.q2 = q2;
  dis.w2 = ev.kin.w2;
  dis.residual = (beam_e_ + beam_ion_) - esc.p;
  const CoherentEvent ce = csampler_[k]->sample(rng, dis);

  // Adds the intact ground-state recoil and sets kin.t / kin.x_pom / channel,
  // and multiplies in the azimuthal weight 1 + c2 cos 2(phi_t - phi_S).
  csampler_[k]->fill_event(ev, ce);
  const int i_rec = find_role(ev, Role::IntactRecoil);
  const Particle& rec = ev.particles[static_cast<std::size_t>(i_rec)];
  fill_fragment_lab(ev, rec, beams_.ion.A, beams_.ion_momentum_per_nucleon);
  // COPIED, not aliased: the push_back below may reallocate `ev.particles`
  // and leave `rec` dangling.
  const Vec4 p_rec = rec.p;

  // 4. the POMERON, P_IP = P_ion - P_recoil.  This is the T2 target of the
  // coherent channel and the exact analogue of `Role::StruckNucleon` on the
  // other channels: q + P_IP is the diffractive system X, so
  // (q + P_IP)^2 = M_X^2 identically and `PythiaBridge` can run its ordinary
  // (W^2, Q^2)-matched surrogate on a PYTHIA Pomeron beam (id 990) with
  // W^2 -> M_X^2 (docs/PYTHIA_BRIDGE.md sec. 12).  Documentation-only,
  // `Status::Intermediate`: it is not final state and never enters the
  // whole-record sum.
  //
  // It is SPACELIKE -- P_IP^2 = t < 0 -- so its `mass` follows the virtual
  // photon's convention above and is written NEGATIVE, -sqrt(|t|).
  {
    Particle ip;
    ip.pdg = 990;
    ip.status = Status::Intermediate;
    ip.role = Role::Pomeron;
    ip.p = beam_ion_ - p_rec;
    ip.mass = -std::sqrt(std::max(-ip.p.m2(), 0.0));
    ip.charge = 0.0;
    ip.mother1 = 1;                       // the beam ion emitted it
    ev.particles.push_back(ip);
  }
  const int i_pom = find_role(ev, Role::Pomeron);

  // WHOLE-NUCLEUS balance: X = k + P_ion - k' - P_recoil; the diffractive
  // system is neutral.
  add_hadronic_x(ev, (beam_e_ + beam_ion_) - esc.p - p_rec, 0.0, i_pom);
}

Vec4 Pipeline::gen_scattered_electron(double x, double y, double phi) const {
  const ScatteredElectron e = scattered_electron(x, y, dis_sampler_->s(),
                                                 beams_.electron_energy);
  const double st = std::sin(e.theta), ct = std::cos(e.theta);
  return {e.e_prime, e.e_prime * st * std::cos(phi),
          e.e_prime * st * std::sin(phi), e.e_prime * ct};
}

void Pipeline::event(std::uint64_t index, Event& out) const {
  const std::size_t k = category_of(index);
  const std::uint64_t local = index - offset_[k];
  Rng rng(cfg_.seed, cfg_.run, k, local);
  switch (cfg_.channel) {
    case PipelineChannel::Inclusive:
      make_inclusive(k, local, index, rng, out);
      break;
    case PipelineChannel::CoherentLi6:
      make_coherent(k, local, index, rng, out);
      break;
    default:
      make_tagged(k, local, index, rng, out);
      break;
  }
  // rc.hpp: the tensor-sector RC weight family.  A PURE FUNCTION of the
  // FINISHED record and of tables built in the constructor -- it consumes no
  // random number, moves no four-vector and does not touch `Event::weight`,
  // which is exactly what makes `--rc tensor-band` bit-for-bit identical to
  // `--rc off` in every one of those (T6).  Null when `PipelineConfig::rc` is
  // Off, and then `Event::rc_weights` stays empty (`Event::reset` cleared it).
  // It sits AFTER the switch because RC applies on every channel, and BEFORE
  // the hadronizer so that a T2 event carries the same weights a T0 one does.
  if (rc_) rc_->fill(out);

  if (cfg_.hadronizer) cfg_.hadronizer(out, rng);
}

Event Pipeline::event(std::uint64_t index) const {
  Event ev;
  event(index, ev);
  return ev;
}

bool Pipeline::next(Event& ev) {
  if (cursor_ >= total_) return false;
  event(cursor_, ev);
  ++cursor_;
  return true;
}

std::uint64_t Pipeline::for_each(const EventSink& sink) const {
  if (!sink) throw std::runtime_error("Pipeline::for_each: null sink");
  Event ev;
  for (std::uint64_t i = 0; i < total_; ++i) {
    event(i, ev);
    sink(ev);
  }
  return total_;
}

std::uint64_t Pipeline::for_each_range(const EventSink& sink,
                                       std::uint64_t first,
                                       std::uint64_t last) const {
  if (!sink) throw std::runtime_error("Pipeline::for_each_range: null sink");
  if (last > total_) last = total_;
  if (first >= last) return 0;
  Event ev;
  for (std::uint64_t i = first; i < last; ++i) {
    event(i, ev);
    sink(ev);
  }
  return last - first;
}

void Pipeline::generate_range(std::uint64_t first, std::uint64_t last,
                              std::vector<Event>& out) const {
  if (last > total_) last = total_;
  if (first >= last) { out.clear(); return; }
  // RESIZE, never clear-then-resize.  `clear()` destroys every record and
  // frees its particle vector, so a chunked loop over one buffer paid a
  // free/malloc pair per event; a plain resize keeps the records that are
  // already there and `Pipeline::event` reconstructs into them in place.
  out.resize(static_cast<std::size_t>(last - first));
  for (std::uint64_t i = first; i < last; ++i) {
    event(i, out[static_cast<std::size_t>(i - first)]);
  }
}

std::uint64_t Pipeline::for_each(const EventSink& sink, unsigned nthreads,
                                 std::size_t chunk) const {
  if (!sink) throw std::runtime_error("Pipeline::for_each: null sink");
  if (nthreads <= 1) return for_each(sink);
  if (chunk == 0) chunk = 1;
  std::vector<std::vector<Event>> buf(nthreads);
  std::vector<std::thread> workers;
  workers.reserve(nthreads);
  std::uint64_t pos = 0;
  while (pos < total_) {
    workers.clear();
    for (unsigned t = 0; t < nthreads; ++t) {
      const std::uint64_t lo = pos + static_cast<std::uint64_t>(t) * chunk;
      const std::uint64_t hi = lo + chunk;
      if (lo >= total_) { buf[t].clear(); continue; }
      workers.emplace_back([this, t, lo, hi, &buf] {
        generate_range(lo, hi, buf[t]);
      });
    }
    for (std::thread& th : workers) th.join();
    for (unsigned t = 0; t < nthreads; ++t) {
      for (const Event& ev : buf[t]) sink(ev);
    }
    pos += static_cast<std::uint64_t>(nthreads) * chunk;
  }
  return total_;
}

}  // namespace lipolgen
