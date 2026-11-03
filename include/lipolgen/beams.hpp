// SPDX-License-Identifier: GPL-3.0-or-later
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
// "THE INCLUSIVE AND TAGGED 6Li CANNOT DRIFT APART" IS TRUE OF THE HULTHEN WAVE
// AND FALSE OF `ClusterWaveSource::VmcAV18` -- THE RUN DEFAULT SINCE D6 (C5.5).
// On the Hulthen option the two ARE one wave function: `LI6_CLUSTER_POLA-
// RIZATION`'s alpha-d factor 1 - 1.5 P_D_LI6 = 0.869950 and the tagged
// `TaggedModel(li6_alpha_channel()).vector_dilution()` = 0.8699431789 agree to
// 7.84e-6 (a grid quadrature against a closed form; 6.82e-6 absolute).
// Re-measured 2026-09-06 on the fixed tagged S-D interference sign, which
// moves the quadrature residual by +4.3e-6: it read 0.869939 / 1.22e-5
// before.  Under `--cluster-wave
// vmc` (THE DEFAULT SINCE D6) the tagged alpha-d wave is the ANL VMC overlap,
// P_D = `VMC_P_D_LI6` = 0.019355, vector dilution 0.970966 -- and until A1 the
// inclusive constant did NOT move, so the two DEFAULT runs differed by +11.61 %
// (rank-2: +6.58 %, `LI6_B1_RANK2_TRANSFER` 0.9219 against tensor_dilution 0.9826),
// registry row 10; A1 (2026-11-03): inclusive 0.848 (tag 0.887076) and 0.982576.
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
// AUTHOR DECISION 2026-09-04: the default stayed 0.811228 (the closer cluster value
// to 0.848).  A1 (2026-11-03; docs/open_items/run_2026-11-03/DECISIONS.md, the
// supervisor's choice under the author's delegation) MOVES IT TO THE ANCHOR:
// `LI6_CLUSTER_POLARIZATION` = `LI6_POLARIZATION_VMC_SIX_BODY` = 0.848; the band is
// the two cluster readings 0.811 ... 0.887 around 0.848 (the ab initio value is the
// default).  Numbers: run_2026-09-03/phase_C_numbers.md sec. C5.5.  Row 10: closed.
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
// Hamiltonian end to end.  (A1: the inclusive deuteron is AV18's too, below.)

/// alpha-d relative D-state probability.  Chosen so that the embedded
/// deuteron's vector dilution 1 - (3/2) P_D reproduces the 0.87 of
/// `b1_li6_from_deuteron` (SCENARIO -- VMC overlaps are the scheduled
/// replacement, plans/04 #15).  It is also the `p_d` of the Hulthen tagged
/// channel (the factories' default; the run default until D6); `VMC_P_D_LI6`
/// (tagged.hpp) is what `ClusterWaveSource::VmcAV18`, the run default, uses.
inline constexpr double P_D_LI6 = 0.0867;
inline constexpr double P_D_DEUTERON_SCENARIO = 0.045;  ///< the pre-A1 "AV18-like" scenario (notes below)
inline constexpr double P_D_DEUTERON = 0.0576;  ///< A1 (2026-11-03): AV18's, the inclusive deuteron

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
    vector_dilution_of(P_D_DEUTERON);  // A1 (2026-11-03): 0.9325 -> 0.9136

/// WHOLE-NUCLEUS vector polarization of the two polarized nucleons of 6Li in
/// the cluster picture (author decision 2026-08-29, plans/04 #6).  The 6Li
/// spin is carried by the alpha-d relative motion and by the deuteron inside
/// it, so a nucleon of that deuteron is polarized along the 6Li spin by the
/// PRODUCT of the two dilutions.  The alpha contributes nothing (J = 0).
///
/// THE ONE HOME OF THE FORMULA (C5.5).  It exists as a function so that the
/// alternative wave-function readings above can be written down without
/// retyping (1 - 1.5 P) anywhere: `LI6_CLUSTER_POLARIZATION_HULTHEN` is this at
/// (P_D_LI6, P_D_DEUTERON_SCENARIO) and `LI6_CLUSTER_POLARIZATION_VMC` (tagged.hpp)
/// the same function at the VMC alpha-d P_D.  Multiplication order is the
/// pre-2026-09-04 one, so 0.811228375 is bit for bit unchanged.
constexpr double li6_cluster_polarization(double p_d_alpha_d,
                                          double p_d_deuteron) {
  return vector_dilution_of(p_d_alpha_d) * vector_dilution_of(p_d_deuteron);
}

/// 0.86995 x 0.9325 = 0.81123: the pre-A1 LI6_CLUSTER_POLARIZATION (the Hulthen reading).
inline constexpr double LI6_CLUSTER_POLARIZATION_HULTHEN =
    li6_cluster_polarization(P_D_LI6, P_D_DEUTERON_SCENARIO);

/// The SIX-BODY VMC answer for the same whole-nucleus quantity, ab initio:
/// R. B. Wiringa et al., PRC 89 (2014) 024305, Table I -> 0.848 -- C5.5's
/// reference point (0.811228 4.3 % below, 0.905427 6.8 % and 0.887076 4.6 % above).
/// Quoted in `LI6()` (beams.cpp) since 2026-08-29.
inline constexpr double LI6_POLARIZATION_VMC_SIX_BODY = 0.848;
/// WHOLE-NUCLEUS vector polarization the default INCLUSIVE 6Li run reads.  A1
/// (2026-11-03): 0.811228 -> 0.848, the six-body value, one home.  Band 0.811 ...
/// 0.887 around 0.848 (the two cluster readings); the pre-A1 value is reached
/// through the API, `li6_ion(LI6_CLUSTER_POLARIZATION_HULTHEN)`.
inline constexpr double LI6_CLUSTER_POLARIZATION = LI6_POLARIZATION_VMC_SIX_BODY;

/// The retired alternative, kept reachable and pinned: Cloet's slides use
/// P_p = P_n = 1/3, i.e. a whole-nucleus Z*P_p = N*P_n = 1 -- one fully
/// polarized proton and neutron out of three each.  It was the default before
/// 2026-08-29 and is 1/0.81123 = 1.233 times the Hulthen cluster value.
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

// --- A1 (2026-11-03): the inclusive factors' named pre-A1 values and their
// parameters (docs/open_items/run_2026-11-03/DECISIONS.md A1, the supervisor's
// choice under the author's delegation).
//
// `P_D_DEUTERON_SCENARIO` = 0.045 (above) is the sibling's "AV18-like"
// scenario, chosen with P_D_LI6 to reproduce the 0.87 of `b1_li6_from_deuteron`;
// it stays the Hulthen deuteron of the TAGGED channels (`deuteron_channel`'s
// and `BreakupOptions::p_d`'s default, `DEUTERON_SCENARIO()` on the Hulthen
// alpha tag) and the polligen port gates' value -- it was `P_D_DEUTERON` until
// A1.  `P_D_DEUTERON` = 0.0576 (above) is AV18's, what `deuteron_av18_p_d()`
// derives from fdeut.av18 (0.0575998920; pinned to it at 5e-7, half a unit of
// the last printed digit, tests/test_tagged.cpp); CD-Bonn's 0.0485 and the
// measured omega_D = 0.05 +- 0.01 bound its band.  `DEUTERON()` is built from it.
/// 0.9325, the scenario deuteron's (the Hulthen tag's embedded deuteron).
inline constexpr double DEUTERON_VECTOR_POLARIZATION_SCENARIO =
    vector_dilution_of(P_D_DEUTERON_SCENARIO);

/// The 6Li Ion at a chosen WHOLE-NUCLEUS vector polarization, a third of it in
/// each per-nucleon slot (Z*P_p = N*P_n = P).  The parameter through which a
/// caller (the polligen port gates) reaches the pre-A1 0.811228; `LI6()` is
/// this at the default.
Ion li6_ion(double whole_nucleus_polarization = LI6_CLUSTER_POLARIZATION);
/// The deuteron Ion at a chosen D-state probability, eff_pol = 1 - 1.5 P_D in
/// both slots.  `DEUTERON()` is this at `P_D_DEUTERON` (AV18's since A1),
/// `DEUTERON_SCENARIO()` at `P_D_DEUTERON_SCENARIO`.
Ion deuteron_ion(double p_d = P_D_DEUTERON);
/// The scenario deuteron (eff_pol 0.9325): the Hulthen alpha tag's embedded
/// deuteron and polligen's `DEUTERON` -- the pre-A1 `DEUTERON()` bit for bit.
const Ion& DEUTERON_SCENARIO();

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
