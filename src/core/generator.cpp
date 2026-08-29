#include "lipolgen/generator.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace lipolgen {

int nuclear_pdg(int z, int a) {
  return 1000000000 + z * 10000 + a * 10;
}

InclusiveGenerator::InclusiveGenerator(
    std::shared_ptr<const InclusiveSampler> sampler, GeneratorConfig config)
    : sampler_(std::move(sampler)), config_(std::move(config)) {
  if (!sampler_) throw std::runtime_error("InclusiveGenerator: null sampler");
  const BeamConfig& bc = sampler_->config();
  const Ion& ion = bc.ion;
  ion_pdg_ = nuclear_pdg(ion.Z, ion.A);
  ion_mass_ = config_.ion_mass > 0.0 ? config_.ion_mass : ion.mass();
  proton_fraction_ = static_cast<double>(ion.Z) / static_cast<double>(ion.A);

  // Massless electron along -z, whole-nucleus ion along +z (reco.py
  // `beam_fourvectors`).
  beam_e_ = {bc.electron_energy, 0.0, 0.0, -bc.electron_energy};
  const double p_a = static_cast<double>(ion.A) * bc.ion_momentum_per_nucleon;
  beam_ion_ = {std::sqrt(p_a * p_a + ion_mass_ * ion_mass_), 0.0, 0.0, p_a};
}

Vec4 InclusiveGenerator::target_nucleon(const EventDraw& draw, Rng& rng) const {
  if (config_.target) return config_.target(draw, rng);
  const double p_u = sampler_->config().ion_momentum_per_nucleon;
  return {std::sqrt(p_u * p_u + M_NUCLEON * M_NUCLEON), 0.0, 0.0, p_u};
}

Vec4 InclusiveGenerator::scattered_electron_p4(double x, double y,
                                               double phi) const {
  const BeamConfig& bc = sampler_->config();
  const ScatteredElectron e =
      scattered_electron(x, y, sampler_->s(), bc.electron_energy);
  const double st = std::sin(e.theta), ct = std::cos(e.theta);
  return {e.e_prime, e.e_prime * st * std::cos(phi),
          e.e_prime * st * std::sin(phi), e.e_prime * ct};
}

Event InclusiveGenerator::make_event(const SpinCategory& cat,
                                     const RunPlan& plan,
                                     const EventDraw& draw, Rng& rng,
                                     std::uint64_t number, int run,
                                     int bunch) const {
  Event ev;
  ev.number = number;
  ev.channel = config_.channel;
  ev.weight = 1.0;

  ev.spin.j = cat.j;
  ev.spin.m_ion = draw.m;
  ev.spin.m_struck = std::nan("");  // inclusive: no struck-cluster label
  ev.spin.lam_e = cat.lam_e;
  ev.spin.pe = cat.pe;
  ev.spin.theta_s = cat.theta_s;
  ev.spin.phi_s = cat.phi_s;
  ev.spin.pz = plan.pz_true();
  ev.spin.pzz = plan.pzz_true();
  ev.spin.category = cat.name;
  ev.spin.run = run;
  ev.spin.bunch = bunch;

  ev.kin.x = draw.x;
  ev.kin.q2 = draw.q2;
  ev.kin.y = draw.y;
  ev.kin.phi = draw.phi;
  ev.kin.s = sampler_->s();
  ev.kin.w2 = w2_from_xq2(draw.x, draw.q2);
  ev.kin.nu = draw.q2 / (2.0 * M_NUCLEON * draw.x);

  // --- particles ---------------------------------------------------------
  Particle beam_e;
  beam_e.pdg = 11;
  beam_e.status = Status::Beam;
  beam_e.role = Role::BeamElectron;
  beam_e.p = beam_e_;
  beam_e.mass = 0.0;
  beam_e.charge = -1.0;
  beam_e.pol = cat.lam_e != 0 ? static_cast<double>(cat.lam_e) : 9.0;

  Particle beam_ion;
  beam_ion.pdg = ion_pdg_;
  beam_ion.status = Status::Beam;
  beam_ion.role = Role::BeamIon;
  beam_ion.p = beam_ion_;
  beam_ion.mass = ion_mass_;
  beam_ion.charge = static_cast<double>(sampler_->config().ion.Z);
  beam_ion.pol = draw.m;

  Particle escat;
  escat.pdg = 11;
  escat.status = Status::Final;
  escat.role = Role::ScatteredElectron;
  escat.p = scattered_electron_p4(draw.x, draw.y, draw.phi);
  escat.mass = 0.0;
  escat.charge = -1.0;
  escat.mother1 = 0;
  escat.pol = beam_e.pol;

  const Vec4 pn = target_nucleon(draw, rng);
  int nucleon_pdg = config_.struck_nucleon_pdg;
  if (nucleon_pdg == 0) {
    nucleon_pdg = (rng.uniform() < proton_fraction_) ? 2212 : 2112;
  }

  Particle nucleon;
  nucleon.pdg = nucleon_pdg;
  nucleon.status = Status::Intermediate;
  nucleon.role = Role::StruckNucleon;
  nucleon.p = pn;
  nucleon.mass = std::sqrt(std::max(pn.m2(), 0.0));
  nucleon.charge = (nucleon_pdg == 2212) ? 1.0 : 0.0;
  nucleon.mother1 = 1;

  ev.particles = {beam_e, beam_ion, escat, nucleon};
  const int i_nucleon = 3;
  int i_gamma = -1;
  if (config_.with_virtual_photon) {
    Particle gamma;
    gamma.pdg = 22;
    gamma.status = Status::Intermediate;
    gamma.role = Role::VirtualPhoton;
    gamma.p = beam_e_ - escat.p;
    gamma.mass = -std::sqrt(std::max(-gamma.p.m2(), 0.0));  // spacelike
    gamma.charge = 0.0;
    gamma.mother1 = 0;
    ev.particles.push_back(gamma);
    i_gamma = static_cast<int>(ev.particles.size()) - 1;
  }

  // X carries the remainder EXACTLY: k + P_N - k'.  On the per-nucleon
  // balance with an ON-SHELL target its mass IS W, so a spacelike X means the
  // `GeneratorConfig::target` hook handed back something the acceptance window
  // cannot support -- a hard check, not a clip (C1).
  Particle x;
  x.pdg = 92;
  x.status = Status::Final;
  x.role = Role::HadronicX;
  x.p = (beam_e_ + pn) - escat.p;
  const double m2_x = x.p.m2();
  if (!(m2_x >= 0.0)) {
    throw std::runtime_error(
        "InclusiveGenerator: the hadronic system X came out SPACELIKE "
        "(M_X^2 = " + std::to_string(m2_x) + " GeV^2) -- check the target hook");
  }
  x.mass = std::sqrt(m2_x);
  x.charge = nucleon.charge;
  x.mother1 = i_nucleon;
  x.mother2 = i_gamma;
  ev.particles.push_back(x);
  return ev;
}

std::vector<double> InclusiveGenerator::sigma_per_category(
    const RunPlan& plan) const {
  std::vector<double> out;
  out.reserve(plan.categories().size());
  for (const SpinCategory& c : plan.categories()) {
    out.push_back(sampler_->sigma_tot_pb(c));
  }
  return out;
}

std::uint64_t InclusiveGenerator::run(const RunPlan& plan,
                                      double total_lumi_pb, std::uint64_t seed,
                                      const EventSink& sink,
                                      std::uint64_t run_number,
                                      bool poisson) const {
  if (!sink) throw std::runtime_error("InclusiveGenerator::run: null sink");
  const std::vector<double> shares = plan.lumi_share_vector(total_lumi_pb);
  std::uint64_t total = 0;
  for (std::size_t k = 0; k < plan.categories().size(); ++k) {
    const SpinCategory& cat = plan.categories()[k];
    const double sigma = sampler_->sigma_tot_pb(cat);
    const double mu = shares[k] * sigma;
    std::uint64_t n = 0;
    if (poisson) {
      Rng counter(seed, run_number, k, kCountStreamEvent);
      n = rng_poisson(counter, mu);
    } else {
      n = static_cast<std::uint64_t>(std::llround(mu));
    }
    const InclusiveSampler::CategoryPlan cplan = sampler_->make_plan(cat);
    for (std::uint64_t i = 0; i < n; ++i) {
      Rng rng(seed, run_number, k, i);
      const EventDraw draw = sampler_->draw_event(cplan, rng);
      Event ev = make_event(cat, plan, draw, rng, total + i,
                            static_cast<int>(run_number), static_cast<int>(k));
      ev.xsec_pb = sigma;
      sink(ev);
    }
    total += n;
  }
  return total;
}

std::uint64_t InclusiveGenerator::run_n(const RunPlan& plan,
                                        std::uint64_t n_total,
                                        std::uint64_t seed,
                                        const EventSink& sink,
                                        std::uint64_t run_number) const {
  if (!sink) throw std::runtime_error("InclusiveGenerator::run_n: null sink");
  const std::size_t nc = plan.categories().size();
  std::vector<double> rate(nc, 0.0);
  double tot = 0.0;
  for (std::size_t k = 0; k < nc; ++k) {
    rate[k] = plan.categories()[k].lumi_fraction *
              sampler_->sigma_tot_pb(plan.categories()[k]);
    tot += rate[k];
  }
  if (!(tot > 0.0)) throw std::runtime_error("run_n: plan has no rate");

  // Largest-remainder split, so the counts sum to n_total exactly.
  std::vector<std::uint64_t> n(nc, 0);
  std::uint64_t assigned = 0;
  for (std::size_t k = 0; k < nc; ++k) {
    n[k] = static_cast<std::uint64_t>(
        std::floor(static_cast<double>(n_total) * rate[k] / tot));
    assigned += n[k];
  }
  for (std::size_t k = 0; assigned < n_total; k = (k + 1) % nc) {
    ++n[k];
    ++assigned;
  }

  std::uint64_t total = 0;
  for (std::size_t k = 0; k < nc; ++k) {
    const SpinCategory& cat = plan.categories()[k];
    const double sigma = sampler_->sigma_tot_pb(cat);
    const InclusiveSampler::CategoryPlan cplan = sampler_->make_plan(cat);
    for (std::uint64_t i = 0; i < n[k]; ++i) {
      Rng rng(seed, run_number, k, i);
      const EventDraw draw = sampler_->draw_event(cplan, rng);
      Event ev = make_event(cat, plan, draw, rng, total + i,
                            static_cast<int>(run_number), static_cast<int>(k));
      ev.xsec_pb = sigma;
      sink(ev);
    }
    total += n[k];
  }
  return total;
}

}  // namespace lipolgen
