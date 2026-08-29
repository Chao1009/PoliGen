#ifndef LIPOLGEN_BEAMS_HPP
#define LIPOLGEN_BEAMS_HPP

/// \file beams.hpp
/// Beam species and EIC energy configurations.  Port of
/// `fastsim/polli_fastsim/beams.py`.
///
/// TWO DIFFERENT CONSTRAINTS SET AN ION ENERGY.
///
/// TOP ENERGY -- magnetic RIGIDITY.  The ring caps at the 275 GeV proton
/// rigidity, so p/nucleon <= 275 GeV * (Z/A).  This reproduces d (137.5),
/// Au (110 GeV/u), 3He (183 GeV/u) and the Li values of the EPIOS white paper
/// (Atoian et al., arXiv:2510.10794): ~138 GeV/u 6Li, ~117 GeV/u 7Li.
///
/// LOWER ENERGIES -- revolution period, i.e. SPEED.  The HSR and the ESR must
/// have equal revolution period and the electrons are ultrarelativistic, so
/// the hadron GAMMA is fixed by the ring circumference.  Ions are therefore
/// GAMMA-MATCHED to the proton configuration, NOT rigidity-scaled.  Yellow
/// Report Table 10.2's gold at 41 GeV/u (gamma 44.02 against the 41 GeV
/// proton's 43.70) is what settles it; rigidity scaling would have put gold
/// at 16.4 GeV/u.
///
/// Nothing here is a rounded literal: 6Li's 40.8 / 99.5 / 137.5 GeV/u and
/// 7Li's 40.8 / 99.5 / 117.9 are DERIVED from the physical (AME2020) nuclear
/// masses and the two constraints, then rounded to 0.1 GeV/u exactly as
/// `default_configs` does.

#include <string>
#include <vector>

#include "lipolgen/constants.hpp"

namespace lipolgen {

/// Physical nuclear masses [GeV] -- the same AME2020-derived values as
/// `polli_fastsim.spectator.NUCLEUS_MASS`.  Mass per NUCLEON is what converts
/// a gamma into a per-nucleon momentum, and it is not A * M_U: a 6Li is bound
/// by ~5.3 MeV/nucleon relative to a free proton.
/// Throws std::runtime_error for an unknown (name, A, Z).
double nucleus_mass(const std::string& name, int a, int z);

struct Ion {
  std::string name;
  int A = 1;
  int Z = 1;
  double spin = 0.5;
  /// PER-NUCLEON effective polarizations P_p, P_n, from which
  ///   g1A = Z*P_p*g1p + N*P_n*g1n
  /// exactly as NuclearF2 builds F2A from Z*f2p + N*f2n (plans/08 D7).
  double eff_pol_p = 0.0;
  double eff_pol_n = 0.0;

  int N() const { return A - Z; }
  double mass() const { return nucleus_mass(name, A, Z); }
  double mass_per_nucleon() const { return mass() / static_cast<double>(A); }
  /// Rigidity-limited top momentum per nucleon [GeV].
  double momentum_per_nucleon_max() const {
    return PROTON_TOP_MOMENTUM * static_cast<double>(Z) / static_cast<double>(A);
  }
  /// Gamma-matched to `proton_energy`, then capped by the ring rigidity.
  double momentum_per_nucleon_at(double proton_energy) const;
};

/// Lorentz factor of the ring at the configuration whose proton energy is
/// `proton_energy` [GeV].  Species-independent.
double gamma_of(double proton_energy);

/// Which EPIOS synchronisation window a configuration falls in:
/// "bypass", "radial-shift", or "" if neither.
std::string epios_window_of(double proton_energy, double shift_gamma = 2.7,
                            double tol = 0.005);

/// The five species of `polli_fastsim.beams`.
const Ion& PROTON();
const Ion& DEUTERON();
const Ion& HE3();
const Ion& LI6();
const Ion& LI7();
/// Lookup by name ("p", "d", "3He", "6Li", "7Li").
const Ion& ion_by_name(const std::string& name);

/// Proton beam energies of the three reference EIC configurations [GeV].
extern const double PROTON_CONFIG_ENERGIES[3];
/// Electron energies of the same three configurations [GeV].
extern const double ELECTRON_ENERGIES[3];

struct BeamConfig {
  double electron_energy = 0.0;      ///< GeV
  Ion ion;
  double ion_momentum_per_nucleon = 0.0;  ///< GeV

  /// sqrt(s) of the electron-nucleon system [GeV] (massless approximation).
  double sqrt_s_per_nucleon() const;
  /// s per nucleon [GeV^2], the argument every kernel call takes.
  double s_per_nucleon() const;
  std::string label() const;
};

/// Reference energy scan low/mid/top, mirroring ep 5x41, 10x100, 18x275.
/// Each point is gamma-matched then rigidity-capped, and rounded to 0.1 GeV/u
/// with Python's round-half-even decimal rounding.
std::vector<BeamConfig> default_configs(const std::string& ion_name = "7Li");

}  // namespace lipolgen

#endif  // LIPOLGEN_BEAMS_HPP
