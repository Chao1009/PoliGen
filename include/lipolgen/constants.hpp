#ifndef LIPOLGEN_CONSTANTS_HPP
#define LIPOLGEN_CONSTANTS_HPP

/// \file constants.hpp
/// Single definition of every physics convention constant of LiPolGen.
///
/// Mirrors `polli_fastsim.asymmetries`, `polli_fastsim.structure`,
/// `polli_fastsim.beams`, `polli_fastsim.polarized`, `polli_fastsim.delta_models`
/// and `polligen.xsec`.  Nothing physical is defined twice anywhere else in
/// the library (docs/CONVENTIONS.md).

namespace lipolgen {

/// Sign of the tensor RATE (b1, b2) sector -- the single place it is set.
///
///   +1  the program's transcription of Hoodbhoy-Jaffe-Manohar
///       (docs/Discussions.pptx p.5), giving Azz = +(2/3) b1/F1
///   -1  the HJM/HERMES convention as written by Cosyn, Roldan Tomei, Sosa
///       and Zec, EPJ A 61 (2025) 83 (arXiv:2410.12764) Eq. (27),
///       Azz = -(2/3) b1/F1
///
/// The two differ by the sign of b1 itself, so |Azz| and the whole Delta
/// (cos 2phi) sector are unaffected.  plans/08 D1 is the open decision;
/// tests/test_xsec.cpp pins the identity it controls.
inline constexpr double TENSOR_LL_SIGN = +1.0;

/// Fine-structure constant, exactly as `polli_fastsim.structure.ALPHA_EM`.
inline constexpr double ALPHA_EM = 1.0 / 137.036;

/// (hbar c)^2 in GeV^2 * pb.
inline constexpr double GEV2_TO_PB = 0.3894e9;

/// Per-nucleon target mass used for gamma = 2 M x / Q (`polligen.xsec`).
/// FREE-nucleon mass, because x is per-nucleon; the bound-nucleon mass
/// would move gamma^2 by 1.0 %.
inline constexpr double M_NUCLEON = 0.9383;

/// Proton mass used by the gamma-matching of `polli_fastsim.beams`.
inline constexpr double PROTON_MASS = 0.938272088;

/// Ring rigidity cap, expressed as the top proton momentum [GeV].
inline constexpr double PROTON_TOP_MOMENTUM = 275.0;

/// EPIOS synchronisation windows (arXiv:2510.10794 pp. 12-13).
inline constexpr double EPIOS_GAMMA_BYPASS = 43.5;
inline constexpr double EPIOS_GAMMA_SHIFT_LO = 118.0;
inline constexpr double EPIOS_GAMMA_SHIFT_HI = 293.0;

/// The published b1 curves are per DEUTERON; every consumer here pairs b1
/// with a per-NUCLEON F1, so the tables are halved on the way out.
inline constexpr double B1_PER_DEUTERON_TO_PER_NUCLEON = 0.5;

/// Rank-2 transfer of the embedded deuteron's tensor polarization to 6Li,
/// 1 - (9/10) P_D at P_D = 0.0867 (plans/08 D9).
inline constexpr double LI6_B1_RANK2_TRANSFER = 0.921947;
/// The pre-2026-08-28 VECTOR dilution 1 - (3/2) P_D -- the wrong rank for b1.
inline constexpr double LI6_B1_LEGACY_TRANSFER = 0.87;
/// 6Li carries two polarized nucleons out of six.
inline constexpr double LI6_B1_PER_NUCLEON = 2.0 / 6.0;

/// Sather-Schmidt bag-model sum-rule coefficient (PRD 42:1424):
/// int_0^1 x Delta dx = C_BAG * alpha_s(Q2).
inline constexpr double C_BAG = -0.012;

/// The valence window over which the two polarized-EMC camps are put on a
/// common 7Li baseline (`polli_fastsim.polarized.POLEMC_VALENCE_WINDOW`).
inline constexpr double POLEMC_VALENCE_WINDOW_LO = 0.35;
inline constexpr double POLEMC_VALENCE_WINDOW_HI = 0.65;

/// pi, spelled once (M_PI is not standard C++).
inline constexpr double kPi = 3.14159265358979323846;

}  // namespace lipolgen

#endif  // LIPOLGEN_CONSTANTS_HPP
