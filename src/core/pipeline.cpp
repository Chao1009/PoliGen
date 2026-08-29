#include "lipolgen/pipeline.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>

#include "lipolgen/constants.hpp"
#include "lipolgen/sf.hpp"
#include "lipolgen/spin.hpp"

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
  // quantization axis along z (the Python's `_pure_category`, which does not
  // pass the ion fill's axis either -- the axis enters the event through the
  // spectator rotation in `boost_spectator`).
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
  const InclusiveSampler::CategoryPlan& plan = plan_for(m_struck, lam_e, pe);
  x.resize(n); q2.resize(n); y.resize(n); phi.resize(n);
  for (std::size_t i = 0; i < n; ++i) {
    const EventDraw d = sampler_->draw_event(plan, rng);
    x[i] = d.x; q2[i] = d.q2; y[i] = d.y; phi[i] = d.phi;
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
    o.b1_func = [](double x, double q2, double f1) {
      return toy_b1(x, q2, f1);
    };
  }
  o.delta_func = opt.delta_func;
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
  if (hadronizer && channel == PipelineChannel::CoherentLi6 &&
      !hadronize_coherent) {
    throw std::runtime_error(
        "PipelineConfig: a hadronizer on the COHERENT channel is refused "
        "(C4).  A coherent event carries no struck nucleon and no struck "
        "cluster, so PythiaBridge v0 falls through to its inclusive branch "
        "and invents a nucleon at rest in the ion frame; that nucleon is not "
        "in the record's balance, so the hadronized event loses "
        "P_ion (1 - 1/A) of four-momentum and Z - 1 of charge.  Set "
        "hadronize_coherent to reproduce that known-broken behaviour "
        "deliberately.");
  }
  const std::string want = channel_isotope(channel);
  if (!want.empty() && want != isotope) {
    throw std::runtime_error(std::string("PipelineConfig: channel ") +
                             pipeline_channel_name(channel) + " needs isotope " +
                             want + ", got " + isotope);
  }
  scenario.validate();
}

std::shared_ptr<const InclusiveKernel> default_inclusive_kernel(const Ion& ion) {
  InclusiveKernel::Options opt;
  if (std::fabs(ion.spin - 1.0) < 1e-9) {
    static const auto li6_b1 =
        std::make_shared<Li6B1>(std::make_shared<MillerB1>());
    opt.b1_func = li6_b1->b1_func();
    opt.delta_func = [](double x, double q2, double f1) {
      return toy_delta_gluon(x, q2, f1, 1e-2);
    };
  }
  return std::make_shared<InclusiveKernel>(ion, opt);
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
            ? li6_alpha_channel(cfg_.cluster_beta, cfg_.p_d)
            : cfg_.channel == PipelineChannel::TaggedLi7Alpha
                  ? li7_alpha_channel(cfg_.cluster_beta)
                  : deuteron_channel(cfg_.cluster_beta, P_D_DEUTERON)));
    model_ = std::make_shared<TaggedModel>(*channel_);
    StruckClusterOptions sopt = cfg_.struck;
    sopt.scenario = cfg_.scenario;
    sopt.grid = cfg_.grid;
    dis_source_ = make_struck_cluster_source(*channel_, beams_, sopt);
    dis_sampler_ = dis_source_->sampler_ptr();
    tsampler_.reset(new TaggedSampler(*model_, p_u, dis_source_.get(), optics_,
                                      pot_config_));

    // Warm every cache the event loop reads.  `TaggedModel`'s amplitude,
    // density and cell-CDF maps are `mutable` and NOT mutex-protected, so
    // building them all here is what makes threaded generation safe; the
    // sampler's own per-state cache is guarded, so warming it is only an
    // anti-contention measure.
    const std::vector<double> ms_ion = m_values(channel_->j_ion);
    const std::vector<double>& ms_c = model_->m_struck_values();
    Rng warm_rng(1, 0, 0, 0);
    for (double mi : ms_ion) {
      (void)model_->n_of_kc(mi);
      const std::vector<double> p_ms = model_->population_integrated(mi);
      for (std::size_t b = 0; b < ms_c.size(); ++b) {
        // A (M, m_S) pair the CG factors forbid (|M - m_S| > L_max) has no
        // density at all and `sample_kc_one` would throw on its empty CDF;
        // `rates()` gives it weight zero, so the event loop never asks.
        if (!(p_ms[b] > 0.0)) continue;
        double kk = 0.0, cc = 0.0, pp = 0.0;
        model_->sample_kc_one(mi, ms_c[b], warm_rng, kk, cc, pp);
      }
    }
    fills_.reserve(nc);
    rate_cdf_.reserve(nc);
    for (std::size_t k = 0; k < nc; ++k) {
      const SpinCategory& c = plan_.categories()[k];
      if (std::fabs(c.j - channel_->j_ion) > 1e-9) {
        throw std::runtime_error("Pipeline: run-plan spin " +
                                 std::to_string(c.j) + " != channel ion spin");
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
  } else {
    const Ion& ion = ion_by_name(cfg_.isotope);
    const auto kernel = cfg_.kernel ? cfg_.kernel : default_inclusive_kernel(ion);
    dis_sampler_ = std::make_shared<InclusiveSampler>(kernel, beams_,
                                                      cfg_.scenario, cfg_.grid);
  }

  if (cfg_.channel == PipelineChannel::Inclusive) {
    GeneratorConfig gc;
    gc.with_virtual_photon = cfg_.with_virtual_photon;
    gc.channel = Channel::Inclusive;
    gen_.reset(new InclusiveGenerator(dis_sampler_, gc));
    cplan_.reserve(nc);
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
      if (fits) acc += sc[c] * cfg_.coherent.coherent_fraction(xc[c]);
      coh_cdf_[c] = acc;
    }
    if (!(acc > 0.0)) {
      throw std::runtime_error(
          "Pipeline: no coherent rate -- no accepted (x, Q2) cell admits a "
          "diffractive mass of at least CoherentXpomModel::m_x_min");
    }
    sigma_coh_pb_ = acc;
    for (double& v : coh_cdf_) v /= acc;
    for (std::size_t k = 0; k < nc; ++k) {
      const SpinCategory& c = plan_.categories()[k];
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

Event Pipeline::make_inclusive(std::size_t k, std::uint64_t local,
                               std::uint64_t index, Rng& rng) const {
  (void)local;
  const SpinCategory& cat = plan_.categories()[k];
  const EventDraw draw = dis_sampler_->draw_event(cplan_[k], rng);
  Event ev = gen_->make_event(cat, plan_, draw, rng, index,
                              static_cast<int>(cfg_.run), static_cast<int>(k));
  ev.xsec_pb = sigma_[k];
  return ev;
}

Event Pipeline::make_tagged(std::size_t k, std::uint64_t local,
                            std::uint64_t index, Rng& rng) const {
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
  TaggedEvent te = tsampler_->sample_one(fills_[k], rate_cdf_[k], rng);
  {
    int tries = 0;
    for (; tries < kTaggedMaxRedraw; ++tries) {
      const Vec4 p_e = dis_source_->scattered_electron_p4(te.x, te.y, te.phi);
      const Vec4 p_s = tsampler_->spectator_p4(te);
      if (((beam_e_ + beam_ion_) - p_e - p_s).m2() >= 0.0) break;
      te = tsampler_->sample_one(fills_[k], rate_cdf_[k], rng);
    }
    if (tries == kTaggedMaxRedraw) {
      throw std::runtime_error(
          "Pipeline: no timelike hadronic system on the tagged channel after "
          + std::to_string(kTaggedMaxRedraw) + " draws");
    }
  }

  Event ev;
  label_event(ev, k, index);
  ev.kin.x = te.x;
  ev.kin.q2 = te.q2;
  ev.kin.y = te.y;
  ev.kin.phi = te.phi;
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

  const int i_spec = find_role(ev, Role::Spectator);
  const int i_clus = find_role(ev, Role::StruckCluster);
  const Particle& spec = ev.particles[static_cast<std::size_t>(i_spec)];
  const Particle& clus = ev.particles[static_cast<std::size_t>(i_clus)];
  // WHOLE-NUCLEUS balance: X = k + P_ion - k' - p_spec = k + P_X - k'.
  add_hadronic_x(ev, (beam_e_ + beam_ion_) - esc.p - spec.p, clus.charge,
                 i_clus);
  return ev;
}

Event Pipeline::make_coherent(std::size_t k, std::uint64_t local,
                              std::uint64_t index, Rng& rng) const {
  (void)local;
  const SpinCategory& cat = plan_.categories()[k];

  // 1. the ion spin projection (label only: the coherent rate is spin
  //    independent and the tensor modulation is an ENSEMBLE coefficient
  //    carried by the recoil azimuth -- `CoherentScenario::a2_m_state` is the
  //    per-m relation, kept analytic).
  const std::vector<double> ms = m_values(cat.j);
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
  Event ev;
  label_event(ev, k, index);
  ev.spin.m_ion = m_ion;
  ev.spin.m_struck = kNaN;
  ev.kin.x = x;
  ev.kin.q2 = q2;
  ev.kin.y = y;
  ev.kin.phi = phi;
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
  // WHOLE-NUCLEUS balance: X = k + P_ion - k' - P_recoil; the diffractive
  // system is neutral.
  add_hadronic_x(ev, (beam_e_ + beam_ion_) - esc.p - rec.p, 0.0, 1);
  return ev;
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
      out = make_inclusive(k, local, index, rng);
      break;
    case PipelineChannel::CoherentLi6:
      out = make_coherent(k, local, index, rng);
      break;
    default:
      out = make_tagged(k, local, index, rng);
      break;
  }
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

void Pipeline::generate_range(std::uint64_t first, std::uint64_t last,
                              std::vector<Event>& out) const {
  if (last > total_) last = total_;
  out.clear();
  if (first >= last) return;
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
