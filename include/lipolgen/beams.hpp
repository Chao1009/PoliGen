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

// --- the 6Li cluster wave function ---------------------------------------
//
// ONE SOURCE OF TRUTH for the two D-state probabilities the 6Li cluster
// picture is built from.  They live HERE, in the module every spin consumer
// already includes, and `tagged.hpp` uses these names rather than keeping its
// own copies (`polli_fastsim.beams`, which `polligen.tagged` re-exports from).
//
// "THE INCLUSIVE AND TAGGED 6Li CANNOT DRIFT APART" IS TRUE OF THE DEFAULT AND
// FALSE OF `ClusterWaveSource::VmcAV18` (open item C5.5, measured 2026-09-04).
// On the Hulthen default the two ARE one wave function: `LI6_CLUSTER_POLA-
// RIZATION`'s alpha-d factor 1 - 1.5 P_D_LI6 = 0.869950 and the tagged
// `TaggedModel(li6_alpha_channel()).vector_dilution()` = 0.869939 agree to
// 1.22e-5 (a grid quadrature against a closed form).  Under `--cluster-wave
// vmc` the tagged channel's alpha-d wave is the ANL VMC overlap with
// P_D = `VMC_P_D_LI6` = 0.019355, whose vector dilution is 0.970966 -- and the
// inclusive constant below does NOT move, so the two differ by +11.61 %.  The
// same drift in the rank-2 sector is +6.58 % (`LI6_B1_RANK2_TRANSFER` = 0.9219
// against the VMC channel's tensor_dilution 0.9826).
//
// IT IS NOT CLOSED BY SUBSTITUTION, AND THAT IS THE POINT.  Feeding
// `VMC_P_D_LI6` into the product below gives 0.905427, which is 6.8 % ABOVE
// the one ab initio number for this very observable (Wiringa PRC 89:024305
// Table I, 0.848 -- see `LI6()` in beams.cpp), while the shipped 0.811228 is
// 4.3 % BELOW it.  Adding the AV18 deuteron's own P_D too (`deuteron_av18_p_d`
// = 0.057600, open item C5.4) gives 0.887076, still 4.6 % above.  So the
// product 1 - 1.5 P_D x 1 - 1.5 P_D spans 0.811 .. 0.905 depending on which
// wave functions it is fed, the ab initio answer sits in the middle of that,
// and the FORMULA is what carries the error, not the choice of P_D.
// AUTHOR DECISION 2026-09-04: the default stays 0.811228 (it is the closer of
// the two to 0.848 and it is pinned bit for bit in validation/reference/), the
// drift under `--cluster-wave vmc` is DOCUMENTED at 11.61 % rather than
// "fixed", and nothing quotes an inclusive 6Li polarization without the band
// 0.81 .. 0.91.  Numbers: docs/open_items/run_2026-09-03/phase_C_numbers.md
// sec. C5.5.
//
// WHAT *IS* CLOSED, AND WAS NOT (sec. C5.5b, 2026-09-04).  The 11.61 % above is
// between two RUNS and is left standing.  Inside ONE `--cluster-wave vmc`
// tagged-alpha run there was a SECOND, smaller and strictly wrong split: the
// alpha-d relative motion came from the ANL VMC AV18+UX overlap while the
// EMBEDDED deuteron stayed on `P_D_DEUTERON` = 0.045 in both places a run reads
// it, so every polarized tagged-alpha observable was 2.069 % high against the
// AV18 deuteron of that same Hamiltonian.  That one is fixed: `DEUTERON_AV18()`
// and `BreakupOptions::source` (tagged.hpp, breakup.hpp) follow the flag, and
// such a run's whole-nucleus reading is now
// `li6_cluster_polarization(VMC_P_D_LI6, deuteron_av18_p_d())` = 0.887076, one
// Hamiltonian end to end.  `LI6_CLUSTER_POLARIZATION` below is untouched.

/// alpha-d relative D-state probability.  Chosen so that the embedded
/// deuteron's vector dilution 1 - (3/2) P_D reproduces the 0.87 of
/// `b1_li6_from_deuteron` (SCENARIO -- VMC overlaps are the scheduled
/// replacement, plans/04 #15).  It is also the DEFAULT of the Hulthen tagged
/// channel, which is the bit-compatibility path; `VMC_P_D_LI6` (tagged.hpp)
/// is the measured replacement `ClusterWaveSource::VmcAV18` uses.
inline constexpr double P_D_LI6 = 0.0867;
/// The deuteron's own D-state probability (AV18-like).
inline constexpr double P_D_DEUTERON = 0.045;

/// Vector depolarization 1 - (3/2) P_D of a spin-1 system with D-state
/// probability P_D.  THE ONE HOME OF THE EXPRESSION (C5.5b, 2026-09-04): it
/// used to be retyped at four sites and is now written once, so a wave
/// function's P_D and the dilution it implies cannot drift apart.  Every value
/// below is bit for bit what the literal gave -- it is the same expression,
/// constexpr-evaluated.
constexpr double vector_dilution_of(double p_d) { return 1.0 - 1.5 * p_d; }

/// The two dilutions the shipped cluster picture is built from.  The deuteron
/// slot below carries the second of these verbatim, so the two ions are built
/// from one expression and their ratio is exact rather than rounded.
inline constexpr double ALPHA_D_VECTOR_POLARIZATION =
    vector_dilution_of(P_D_LI6);
inline constexpr double DEUTERON_VECTOR_POLARIZATION =
    vector_dilution_of(P_D_DEUTERON);

/// WHOLE-NUCLEUS vector polarization of the two polarized nucleons of 6Li in
/// the cluster picture (author decision 2026-08-29, plans/04 #6).  The 6Li
/// spin is carried by the alpha-d relative motion and by the deuteron inside
/// it, so a nucleon of that deuteron is polarized along the 6Li spin by the
/// PRODUCT of the two dilutions.  The alpha contributes nothing (J = 0).
///
/// THE ONE HOME OF THE FORMULA (C5.5).  It exists as a function so that the
/// alternative wave-function readings above can be written down without
/// retyping (1 - 1.5 P) anywhere: `LI6_CLUSTER_POLARIZATION` is this at
/// (P_D_LI6, P_D_DEUTERON) and `LI6_CLUSTER_POLARIZATION_VMC` (tagged.hpp) is
/// the same function at the VMC alpha-d P_D.  Multiplication order is the
/// pre-2026-09-04 one, so the default is bit for bit unchanged.
constexpr double li6_cluster_polarization(double p_d_alpha_d,
                                          double p_d_deuteron) {
  return vector_dilution_of(p_d_alpha_d) * vector_dilution_of(p_d_deuteron);
}

/// 0.86995 x 0.9325 = 0.81123.
inline constexpr double LI6_CLUSTER_POLARIZATION =
    li6_cluster_polarization(P_D_LI6, P_D_DEUTERON);

/// The SIX-BODY VMC answer for the same whole-nucleus quantity, ab initio:
/// R. B. Wiringa et al., PRC 89 (2014) 024305, Table I -> 0.848.  It is NOT
/// what the library uses (the cluster product is), and it is the reference
/// point that decides C5.5: the shipped 0.811228 is 4.3 % below it and the
/// "consistent" VMC alpha-d reading 0.905427 is 6.8 % ABOVE it, so making the
/// cluster formula consistent with the tagged VMC wave function would move the
/// inclusive number AWAY from the only ab initio anchor there is.  Quoted in
/// `LI6()` (beams.cpp) since 2026-08-29; given a name here so the band can be
/// checked instead of read.
inline constexpr double LI6_POLARIZATION_VMC_SIX_BODY = 0.848;

/// The retired alternative, kept reachable and pinned: Cloet's slides use
/// P_p = P_n = 1/3, i.e. a whole-nucleus Z*P_p = N*P_n = 1 -- one fully
/// polarized proton and neutron out of three each.  It was the default before
/// 2026-08-29 and is 1/0.81123 = 1.233 times the cluster value.
inline constexpr double LI6_NAIVE_ONE_THIRD = 1.0 / 3.0;

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
