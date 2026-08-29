# Coding and physics conventions

- Namespace `lipolgen`; headers in `include/lipolgen/<module>.hpp`; sources in
  `src/core/`, `src/lhapdf/`, `src/hepmc/`, `src/pythia/`; tests `tests/test_<module>.cpp`
  (doctest, one `TEST_CASE` per identity, reference JSON in `validation/reference/`).
- Units: GeV, GeV², fm only where stated; angles rad; azimuth φ ∈ [0, 2π).
- Frame: head-on. Ion +z, electron −z. Lab (25 mrad crossing) only in an explicit transform.
- Spin: populations ordered m = +J … −J. Quantization axis n̂(θ_S, φ_S) in the head-on frame.
  `TENSOR_LL_SIGN = +1.0` — one constant, defined once (`constants.hpp`), test-guarded.
- Structure-function inputs are `Backend` interfaces: toy implementation always available,
  table/LHAPDF implementations optional. No physics number is hard-coded in two places.
- Randomness: counter-based stream keyed by (seed, run, bunch, event); never a global RNG.
- Every double compared to `polligen` is compared at rtol 1e-12 unless the reference itself is MC.
- No exceptions for control flow in the event loop; errors in setup throw `std::runtime_error`.

## Physics defaults that are a CHOICE, and where the single copy lives

- **Struck nucleon mass.** The implicit struck nucleon of the per-nucleon
  subsystem is on shell at the FREE nucleon mass `M_NUCLEON` (0.9383), never at
  `Ion::mass_per_nucleon()` (0.9338 for 6Li).  Every per-nucleon label on the
  record — `w2_from_xq2`, `Kinematics::nu` — is built on `M_NUCLEON`, so a
  target at the ion's mass per nucleon would describe a different object from
  its own kinematics.  The binding energy is not lost: it is carried by the
  (A−1) remnant, which the per-nucleon balance deliberately does not write.
  Both places that build that nucleon —
  `InclusiveGenerator::target_nucleon` and `PythiaBridge`'s implicit-target
  branch — use `M_NUCLEON` (P3).
- **Struck nucleon species.** Drawn `Z F2p(x, Q²) : N F2n(x, Q²)` at the
  event's own kinematics, because that is what the inclusive rate
  `Z F2p + N F2n` is made of — never flat `Z : N`, which is its x-independent
  limit and is wrong wherever `F2n/F2p ≠ 1` (0.5 against a true 0.615 for 6Li
  at x = 0.5).  `InclusiveGenerator::proton_fraction` is the one rule;
  `NucleonChoice::ByStructureFunctions` applies it in the bridge (P1).
- **Target mass.** `InclusiveKernel::Options::target_mass` defaults to TRUE,
  matching `xsec.py`.  Identities written against the massless
  `A_par = D(y) g1/F1` must construct the massless kernel explicitly.
  `target_mass = true` with `G2Mode::kZero` is legitimate — it is the
  `g2_scale = 0` twist-3 variation, and `g2_scale` is the knob that spans it.
- **Polarized-EMC baseline.** The CBT and TMT curves are quoted on nuclei that
  are not ours, so each is transferred by a valence scale
  `⟨1 − R_unpol,baseline⟩ / ⟨1 − R_unpol,table⟩` over
  `POLEMC_VALENCE_WINDOW`.  That baseline is a CHOICE and is named:
  `EmcBaseline::Epps21` (the default, matching
  `polli_fastsim.polarized.POLEMC_BASELINE`) or `EmcBaseline::LegacyTable`
  (CBT on itself — the pre-2026-08-29 constants 1 and 0.397009).  Both scales
  come out of one code path; the EPPS21 depletion is the single stored number
  `EMC_VALENCE_DEPLETION_EPPS21`, because computing it needs LHAPDF.
- **Coherent |t| range.** `COHERENT_T_MAX_DEFAULT = 0.2 GeV²`.  The cos 2φ
  coefficient is linear and unbounded in |t| and crosses −1 at |t| = 0.245 for
  P_zz = −2, so a larger range makes the azimuthal weight negative;
  `CoherentSampler` checks `CoherentScenario::positivity_margin` at SETUP and
  throws, mirroring `InclusiveKernel::positivity_margin`.
- **Coherent x_P.** Per-nucleon pomeron fraction,
  `x_P = (M_X² + Q²)/(W² + Q²)` with the per-nucleon W², drawn log-uniform on
  `[x_P(M_X,min), 0.1]`; the nucleus loses `x_P/A` of its own light-cone
  momentum.  `β = x/x_P` and `M_X²` are on the record.  X must be TIMELIKE on
  every channel — a hard check in `Pipeline::add_hadronic_x` and in
  `InclusiveGenerator`, never a clip.
- **Luminosity shares.** `Optics::lumi_fraction` multiplies COUNTS and never
  cross sections (`PipelineConfig::apply_optics_lumi_fraction`, default true) —
  the same share rule that keeps `Scenario::run_share` out of
  `sigma_per_category_pb()`.
