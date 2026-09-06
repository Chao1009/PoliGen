// SPDX-License-Identifier: GPL-3.0-or-later
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
///   -1  THE CONVENTION OF THIS PROGRAM SINCE 2026-08-29 (author decision,
///       plans/08 D1): the LITERATURE one, as written by Cosyn, Roldan Tomei,
///       Sosa and Zec, EPJ A 61 (2025) 83 (arXiv:2410.12764) Eq. (27) and
///       used by HERMES,
///
///           Azz = -(2/3) b1/F1     (axis along q, Bjorken limit)
///
///       with P_zz = n+ + n- - 2 n0 the tensor polarization, i.e. b1 > 0
///       means the m = 0 state is the one with the LARGER cross section.
///   +1  the repository's own private convention until that date -- the
///       transcription of Hoodbhoy-Jaffe-Manohar in docs/Discussions.pptx
///       p.5, giving Azz = +(2/3) b1/F1.  Kept reachable by setting this
///       constant back, which is the whole of the change: nothing else in
///       the program knows the sign.
///
/// The decision was taken on the literature and not on a new derivation.
/// The two differ by the sign of b1 itself, so |Azz| and the whole Delta
/// (cos 2phi) sector are unaffected; what flips is the sign of Azz at fixed
/// b1, of the by-product kappa of the spin-state ratio, and of any b-sector
/// subtraction built on kappa -- including the O(gamma^2) tensor leakage into
/// cos 2phi (`InclusiveKernel::Options::tensor_gamma`), which is why that
/// switch was gated on this decision.  `tests/test_xsec.cpp` pins the
/// identity it controls, in Cosyn's own form.
inline constexpr double TENSOR_LL_SIGN = -1.0;

/// Fine-structure constant, exactly as `polli_fastsim.structure.ALPHA_EM`.
inline constexpr double ALPHA_EM = 1.0 / 137.036;

/// (hbar c)^2 in GeV^2 * pb.
inline constexpr double GEV2_TO_PB = 0.3894e9;

/// hbar c [GeV fm]: THE one conversion between fm^-1 and GeV (0.19733, the
/// value BeAGLE and the Ciofi degli Atti-Simula transcription carry).  Single
/// definition; `coherent.hpp`'s `GEV_PER_FM_INV` is an alias of it and
/// `triton_sf.hpp` converts the CS coefficients with it.
inline constexpr double HBARC_GEV_FM = 0.19733;

/// Per-nucleon target mass used for gamma = 2 M x / Q (`polligen.xsec`).
/// FREE-nucleon mass, because x is per-nucleon; the bound-nucleon mass
/// would move gamma^2 by 1.0 %.
inline constexpr double M_NUCLEON = 0.9383;

/// Proton mass used by the gamma-matching of `polli_fastsim.beams`.
inline constexpr double PROTON_MASS = 0.938272088;

/// Electron mass [GeV], the PDG value.  ONE definition, per
/// docs/CONVENTIONS.md ("no physics number is hard-coded in two places"):
///   * `src/hepmc/hepmc_writer.cpp` stamps it as the GENERATED mass of an
///     electron the core builds massless, so that Geant4/DD4hep does not
///     nudge E by O(10 ppm) to satisfy E^2 - p^2 >= m_e^2;
///   * `rc.hpp`'s `RcOptions::m_lepton` RESERVES it for the unimplemented
///     `RcTailModel::PolradFull`, where the LEPTON mass is what keeps
///     POLRAD's infrared factor F_IR finite and what sets l_m = ln(Q^2/m^2).
///     The shipped `TPeak` tail carries no lepton mass, so `hepmc_writer` is
///     the only consumer that reads this constant today, and
///     `PipelineConfig::validate()` refuses a non-default `m_lepton`.
/// The literal is unchanged from the one it replaces (tests/test_hepmc.cpp
/// and python/tests/test_hepmc.py pin the written value).
inline constexpr double M_ELECTRON = 0.51099895e-3;

/// Ring rigidity cap, expressed as the top proton momentum [GeV].
inline constexpr double PROTON_TOP_MOMENTUM = 275.0;

/// EPIOS synchronisation windows (arXiv:2510.10794 pp. 12-13).
inline constexpr double EPIOS_GAMMA_BYPASS = 43.5;
inline constexpr double EPIOS_GAMMA_SHIFT_LO = 118.0;
inline constexpr double EPIOS_GAMMA_SHIFT_HI = 293.0;

/// b1 NORMALISATION.  Every consumer here pairs b1 with a per-NUCLEON F1, but
/// the two digitized deuteron curves are NOT normalised to the same thing, so
/// since 2026-09-03 they no longer share one applied factor.  The arbiter is
/// docs/open_items/run_2026-09-03/phase_A_miller_normalisation.md, whose
/// pivot is that HERMES's PUBLISHED b1_d -- the data both camps plot against
/// -- is per NUCLEON, because their Eq. (5) divides by an F1_d built from
/// F2_d = (F2_p + F2_n)/2 (confirmed by inverting their own Table II in all
/// six bins).
///
/// 1/A at A = 2: the conversion of a per-DEUTERON b1 into a per-NUCLEON one.
inline constexpr double B1_PER_DEUTERON_TO_PER_NUCLEON = 0.5;

/// Applied to `tables::kB1Miller` by `toy_b1()` (sf.cpp).  Miller (PRC
/// 89:045203) is PER DEUTERON: his Eq. (1) number densities are "in a target
/// hadron", his Eq. (5) is a light-cone correlator in the normalised deuteron
/// state, and Eq. (6)'s 1/2 is fully consumed by the quark-spin average of a
/// SPINLESS pion (q_up = q_down = Delta q / 2), so no 1/A is left anywhere in
/// Eqs. (1)/(5)/(6)/(20).  LIKELY, NOT CERTAIN -- his Table I transcribes
/// HERMES's per-nucleon numbers unrescaled, his Fig. 5 overlays them on this
/// curve and he tunes P_6q to one of them, so the paper is self-inconsistent
/// by exactly this factor and only the author (or a numerical reproduction of
/// his Eq. (20)) can close it.  Keeping 0.5 is the status quo; it is recorded
/// as an AUTHOR DECISION in docs/OPEN_ITEMS_SOLUTIONS.md section 10.
inline constexpr double B1_MILLER_TABLE_TO_PER_NUCLEON =
    B1_PER_DEUTERON_TO_PER_NUCLEON;

/// Applied to `tables::kB1CdksQ2p5` by `b1_convolution()` (sf.cpp).  CDKS
/// (PRD 95:074036) is ALREADY per nucleon -- their Eq. (10) spectral function
/// carries an explicit 1/A, the text under their Eq. (16) says in words "the
/// structure function b1 is defined by the one per nucleon", their f(y) is
/// normalised to ONE nucleon and their F1^N = (F1_p + F1_n)/2 -- so the
/// column is NOT converted.  CERTAIN.  This was 0.5 until 2026-09-03, which
/// made `CdksB1` a factor 2 low; `cdks_b1_raw_per_nucleon()` (b1_nuclear.hpp)
/// existed to route around it and is now merely a second name for the same
/// normalisation.
inline constexpr double B1_CDKS_TABLE_TO_PER_NUCLEON = 1.0;

/// Rank-2 transfer of the embedded deuteron's tensor polarization to 6Li:
/// `TaggedModel(li6_alpha_channel()).tensor_dilution()` evaluated at
/// `P_D_LI6` (beams.hpp), pinned in tests/test_tagged.cpp (plans/08 D9).
/// It is a QUADRATURE over the channel and not a closed form: the closed form
/// 1 - (9/10) P_D gives 0.921970 against the 0.9219490 measured, so the two
/// are pinned to each other at 1e-4 and not asserted equal.  (0.9219467 until
/// the 2026-09-06 tagged S-D sign fix, which moves this quadrature residual by
/// +2.46e-6 -- against a 1e-4 pin, so the constant itself does not move.)
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
