#include "lipolgen/bookkeeping.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace lipolgen {

double SpinCategory::vector_moment() const {
  const std::vector<double> ms = m_values(j);
  if (populations.size() != ms.size()) {
    throw std::runtime_error("populations must have 2j+1 entries (m=+J..-J)");
  }
  double v = 0.0;
  for (std::size_t i = 0; i < ms.size(); ++i) v += ms[i] * populations[i];
  return v / j;
}

RunPlan::RunPlan(std::vector<SpinCategory> categories, double pe_true,
                 double pz_true, double pzz_true, double delta_p_over_p,
                 std::uint64_t polarimetry_seed)
    : categories_(std::move(categories)), pe_true_(pe_true), pz_true_(pz_true),
      pzz_true_(pzz_true), delta_p_over_p_(delta_p_over_p),
      polarimetry_seed_(polarimetry_seed) {
  double tot = 0.0;
  for (const SpinCategory& c : categories_) tot += c.lumi_fraction;
  if (!(tot > 0.0)) throw std::runtime_error("run plan has no luminosity");

  measured_pe_ = pe_true_;
  measured_pz_ = pz_true_;
  measured_pzz_ = pzz_true_;
  if (delta_p_over_p_ != 0.0) {
    // One fixed stream, three draws in the order (pe, pz, pzz) -- the same
    // order as the Python's dict literal, so the three scales stay distinct
    // and reproducible for a given seed.
    Rng rng(polarimetry_seed_, 0, 0, 0);
    measured_pe_ = pe_true_ * (1.0 + delta_p_over_p_ * rng.normal());
    measured_pz_ = pz_true_ * (1.0 + delta_p_over_p_ * rng.normal());
    measured_pzz_ = pzz_true_ * (1.0 + delta_p_over_p_ * rng.normal());
  }
}

std::map<std::string, double> RunPlan::lumi_shares(double total_lumi_pb) const {
  std::map<std::string, double> out;
  for (const SpinCategory& c : categories_) {
    out[c.name] = total_lumi_pb * c.lumi_fraction;
  }
  return out;
}

std::vector<double> RunPlan::lumi_share_vector(double total_lumi_pb) const {
  std::vector<double> out;
  out.reserve(categories_.size());
  for (const SpinCategory& c : categories_) {
    out.push_back(total_lumi_pb * c.lumi_fraction);
  }
  return out;
}

int RunPlan::index_of(const std::string& name) const {
  for (std::size_t i = 0; i < categories_.size(); ++i) {
    if (categories_[i].name == name) return static_cast<int>(i);
  }
  return -1;
}

// --------------------------------------------------------- standard plans

namespace {

std::vector<double> to_vec(const std::array<double, 3>& a) {
  return {a[0], a[1], a[2]};
}
std::vector<double> to_vec(const std::array<double, 4>& a) {
  return {a[0], a[1], a[2], a[3]};
}

}  // namespace

RunPlan helicity_flip_plan(double j, double pz, double pe,
                           const HelicityFlipOptions& opt) {
  std::vector<double> pops;
  if (!opt.use_explicit_pzz) {
    pops = populations_maxent(j, pz);
  } else if (std::fabs(j - 1.0) < 1e-9) {
    pops = to_vec(spin1_populations(pz, opt.pzz));
  } else if (std::fabs(j - 1.5) < 1e-9) {
    // Fixes (P_z, T) and leaves the octupole moment R_3 at its default 0,
    // whereas the max-entropy ladder above carries R_3 != 0 for every P_z
    // (docs/theory/SPIN32_FINITE_GAMMA.md sec. 2.4).  R_3 reaches only
    // beam-helicity-odd observables (sec. 4), so no unpolarized-beam number
    // depends on which of the two a 7Li run uses.
    pops = to_vec(spin32_populations(pz, opt.pzz));
  } else if (std::fabs(j - 0.5) < 1e-9) {
    pops = {(1.0 + pz) / 2.0, (1.0 - pz) / 2.0};
  } else {
    throw std::runtime_error("helicity_flip_plan: unsupported spin");
  }
  const double lp = 0.5 * (1.0 + opt.rel_lumi_offset);
  const double lm = 0.5;
  std::vector<SpinCategory> cats{
      SpinCategory(opt.name + "+", j, pops, +1, pe, opt.theta_s, opt.phi_s, lp),
      SpinCategory(opt.name + "-", j, pops, -1, pe, opt.theta_s, opt.phi_s, lm)};

  // A spin-1/2 fill has no rank-2 moment at all: exactly 0.0, never a divisor.
  const double pzz_true =
      (j >= 1.0 - 1e-9) ? moments_along_axis(j, pops).tensor : 0.0;
  return RunPlan(std::move(cats), pe, pz, pzz_true);
}

RunPlan tensor_thirds_plan(double pz, double pzz, double rel_lumi_offset,
                           double theta_s, double phi_s,
                           const std::string& name) {
  const std::vector<double> pops_p = to_vec(spin1_populations(pz, pzz));
  const std::vector<double> pops_m = to_vec(spin1_populations(-pz, pzz));
  const std::vector<double> pops_0 = to_vec(spin1_populations(0.0, -2.0 * pzz));
  const double f0 = (1.0 + rel_lumi_offset) / 3.0;
  const double fpm = 1.0 / 3.0;
  std::vector<SpinCategory> cats{
      SpinCategory(name + "+", 1.0, pops_p, 0, 0.0, theta_s, phi_s, fpm),
      SpinCategory(name + "-", 1.0, pops_m, 0, 0.0, theta_s, phi_s, fpm),
      SpinCategory(name + "0", 1.0, pops_0, 0, 0.0, theta_s, phi_s, f0)};
  return RunPlan(std::move(cats), 0.0, pz, pzz);
}

RunPlan transverse_tensor_plan(double pzz, double phi_s,
                               const std::string& name) {
  const std::vector<double> pops = to_vec(spin1_populations(0.0, pzz));
  std::vector<SpinCategory> cats{
      SpinCategory(name, 1.0, pops, 0, 0.0, kPi / 2.0, phi_s, 1.0)};
  return RunPlan(std::move(cats), 0.0, 0.0, pzz);
}

RunPlan tensor_flip_plan(double pzz, double phi_s, double share_plus,
                         double rel_lumi_offset, const std::string& name) {
  const std::vector<double> pops_p = to_vec(spin1_populations(0.0, pzz));
  const std::vector<double> pops_0 = to_vec(spin1_populations(0.0, -2.0 * pzz));
  std::vector<SpinCategory> cats{
      SpinCategory(name + "+", 1.0, pops_p, 0, 0.0, kPi / 2.0, phi_s,
                   share_plus * (1.0 + rel_lumi_offset)),
      SpinCategory(name + "0", 1.0, pops_0, 0, 0.0, kPi / 2.0, phi_s,
                   1.0 - share_plus)};
  return RunPlan(std::move(cats), 0.0, 0.0, pzz);
}

RunPlan with_offset(const RunPlan& plan, const std::string& category_name,
                    double offset) {
  if (plan.index_of(category_name) < 0) {
    throw std::runtime_error("with_offset: no category " + category_name);
  }
  std::vector<SpinCategory> cats = plan.categories();
  for (SpinCategory& c : cats) {
    if (c.name == category_name) c.lumi_fraction *= (1.0 + offset);
  }
  return RunPlan(std::move(cats), plan.pe_true(), plan.pz_true(),
                 plan.pzz_true(), plan.delta_p_over_p(),
                 plan.polarimetry_seed());
}

double azz_rel_lumi_bias(double offset, double pzz) {
  return -(2.0 / 3.0) * offset / pzz;
}

double apar_rel_lumi_bias(double offset, double pe, double pz) {
  return offset / (2.0 * pe * pz);
}

// --------------------------------------------------- spin-temperature ladder

SpinTemperatureLadder spin_temperature_ladder(double j, double pz,
                                              int iterations) {
  const std::vector<double> ms = m_values(j);
  auto vector_of = [&](double t) {
    double num = 0.0, den = 0.0;
    for (double m : ms) {
      const double w = std::pow(t, m);
      num += m * w;
      den += w;
    }
    return num / den / j;
  };
  double lo = 1e-12, hi = 1e12;
  for (int i = 0; i < iterations; ++i) {
    const double mid = std::sqrt(lo * hi);  // bisection in log t
    if (vector_of(mid) < pz) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  SpinTemperatureLadder out;
  out.t = std::sqrt(lo * hi);
  double sum = 0.0;
  out.populations.reserve(ms.size());
  for (double m : ms) {
    const double w = std::pow(out.t, m);
    out.populations.push_back(w);
    sum += w;
  }
  for (double& p : out.populations) p /= sum;
  return out;
}

double spin_temperature_pzz(double j, double pz) {
  if (j < 1.0 - 1e-9) return 0.0;
  return moments_along_axis(j, populations_maxent(j, pz)).tensor;
}

}  // namespace lipolgen
