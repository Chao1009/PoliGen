// SPDX-License-Identifier: GPL-3.0-or-later
#include "lipolgen/beams.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <stdexcept>
#include <tuple>

namespace lipolgen {
namespace {

/// Python's `round(v, 1)`: correct decimal rounding of the double, ties to
/// even.  `std::round(v*10)/10` is not the same object; going through the
/// correctly-rounded decimal conversion is.
double round1(double v) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.1f", v);
  return std::strtod(buf, nullptr);
}

const std::map<std::tuple<std::string, int, int>, double>& mass_table() {
  static const std::map<std::tuple<std::string, int, int>, double> t = {
      {{"p", 1, 1}, 0.938272088},
      {{"d", 2, 1}, 1.875612942},
      {{"3He", 3, 2}, 2.808391607},
      {{"6Li", 6, 3}, 5.601518702},
      {{"7Li", 7, 3}, 6.533833028},
  };
  return t;
}

}  // namespace

double nucleus_mass(const std::string& name, int a, int z) {
  const auto& t = mass_table();
  const auto it = t.find(std::make_tuple(name, a, z));
  if (it == t.end()) throw std::runtime_error("unknown nucleus " + name);
  return it->second;
}

double Ion::momentum_per_nucleon_at(double proton_energy) const {
  const double gamma = std::sqrt(proton_energy * proton_energy
                                 + PROTON_MASS * PROTON_MASS) / PROTON_MASS;
  const double p_u = mass_per_nucleon() * std::sqrt(gamma * gamma - 1.0);
  return std::min(p_u, momentum_per_nucleon_max());
}

double gamma_of(double proton_energy) {
  return std::sqrt(proton_energy * proton_energy + PROTON_MASS * PROTON_MASS)
         / PROTON_MASS;
}

std::string epios_window_of(double proton_energy, double shift_gamma,
                            double tol) {
  const double g = gamma_of(proton_energy);
  if (std::fabs(g - EPIOS_GAMMA_BYPASS) <= shift_gamma) return "bypass";
  if (EPIOS_GAMMA_SHIFT_LO * (1.0 - tol) <= g
      && g <= EPIOS_GAMMA_SHIFT_HI * (1.0 + tol)) {
    return "radial-shift";
  }
  return "";
}

const Ion& PROTON() {
  static const Ion i{"p", 1, 1, 0.5, 1.0, 0.0};
  return i;
}
const Ion& DEUTERON() {
  // The deuteron's slot IS the expression 6Li's is built from, so the two
  // agree bit for bit and the ratio of their per-nucleon g1 is
  // (1 - 1.5 P_D_LI6)/3 exactly.
  static const Ion i{"d", 2, 1, 1.0, DEUTERON_VECTOR_POLARIZATION,
                     DEUTERON_VECTOR_POLARIZATION};
  return i;
}
const Ion& HE3() {
  // Bissey PRC 65:064317, per-nucleon already.
  static const Ion i{"3He", 3, 2, 0.5, -0.028, 0.86};
  return i;
}
const Ion& LI6() {
  // The CLUSTER PICTURE, since the author decision of 2026-08-29 that closed
  // plans/04 #6: the 6Li spin is carried by the deuteron cluster, so a nucleon
  // of that deuteron is polarized along the 6Li spin by the product of the
  // alpha-d and deuteron vector dilutions, LI6_CLUSTER_POLARIZATION =
  // 0.86995 x 0.9325 = 0.81123 whole-nucleus (Schellingerhout PRC 48:2714),
  // and the slots hold a THIRD of it each so that Z*P_p = N*P_n = 0.81123.
  // Built from the same two D-state probabilities the tagged sector uses
  // (`P_D_LI6`, `P_D_DEUTERON` in this header), so the inclusive and tagged
  // 6Li share one wave function.  Per-nucleon g1(6Li)/g1(d) is therefore
  // (1 - 1.5 P_D_LI6)/3 = 0.290 -- the deuteron's own dilution cancels
  // between the two isoscalar ions -- against the 0.358 the retired
  // `LI6_NAIVE_ONE_THIRD` gave.  NOT adopted, and the upper end of the band:
  // the six-body VMC of Wiringa PRC 89:024305 Table I reads the same
  // whole-nucleus quantity ab initio as 0.848; 0.81-0.85 is the band.
  static const Ion i{"6Li", 6, 3, 1.0, LI6_CLUSTER_POLARIZATION / 3.0,
                     LI6_CLUSTER_POLARIZATION / 3.0};
  return i;
}
const Ion& LI7() {
  // The verified whole-nucleus VMC sums P_p = +0.866, P_n = -0.037 stored
  // DIVIDED BY Z = 3 and N = 4, so Z*P_p and N*P_n return them exactly.
  static const Ion i{"7Li", 7, 3, 1.5, 0.866 / 3.0, -0.037 / 4.0};
  return i;
}

const Ion& ion_by_name(const std::string& name) {
  if (name == "p") return PROTON();
  if (name == "d") return DEUTERON();
  if (name == "3He") return HE3();
  if (name == "6Li") return LI6();
  if (name == "7Li") return LI7();
  throw std::runtime_error("unknown ion " + name);
}

const double PROTON_CONFIG_ENERGIES[3] = {41.0, 100.0, 275.0};
const double ELECTRON_ENERGIES[3] = {5.0, 10.0, 18.0};

double BeamConfig::sqrt_s_per_nucleon() const {
  return std::sqrt(4.0 * electron_energy * ion_momentum_per_nucleon);
}

double BeamConfig::s_per_nucleon() const {
  const double r = sqrt_s_per_nucleon();
  return r * r;
}

std::string BeamConfig::label() const {
  char buf[128];
  std::snprintf(buf, sizeof(buf), "e(%g) x %s(%g/u)", electron_energy,
                ion.name.c_str(), ion_momentum_per_nucleon);
  return std::string(buf);
}

std::vector<BeamConfig> default_configs(const std::string& ion_name) {
  const Ion& ion = ion_by_name(ion_name);
  std::vector<BeamConfig> out;
  out.reserve(3);
  for (int i = 0; i < 3; ++i) {
    BeamConfig c;
    c.electron_energy = ELECTRON_ENERGIES[i];
    c.ion = ion;
    c.ion_momentum_per_nucleon =
        round1(ion.momentum_per_nucleon_at(PROTON_CONFIG_ENERGIES[i]));
    out.push_back(c);
  }
  return out;
}

}  // namespace lipolgen
