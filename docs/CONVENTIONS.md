# Coding and physics conventions

- Namespace `lipolgen`; headers in `include/lipolgen/<module>.hpp`; sources in
  `src/core/`, `src/lhapdf/`, `src/hepmc/`, `src/pythia/`; tests `tests/test_<module>.cpp`
  (doctest, one `TEST_CASE` per identity, reference JSON in `validation/reference/`).
- Units: GeV, GeV², fm only where stated; angles rad; azimuth φ ∈ [0, 2π).
- Frame: head-on. Ion +z, electron −z. Lab (25 mrad crossing) only in an explicit transform.
- Spin: populations ordered m = +J … −J. Quantization axis n̂(θ_S, φ_S) in the head-on frame.
  `TENSOR_LL_SIGN = -1.0` — one constant, defined once (`constants.hpp`), test-guarded.
- Structure-function inputs are `Backend` interfaces: toy implementation always available,
  table/LHAPDF implementations optional. No physics number is hard-coded in two places.
- Randomness: counter-based stream keyed by (seed, run, bunch, event); never a global RNG.
- Every double compared to `polligen` is compared at rtol 1e-12 unless the reference itself is MC.
- No exceptions for control flow in the event loop; errors in setup throw `std::runtime_error`.

## Physics defaults that are a CHOICE, and where the single copy lives

- **Tensor-sector sign.** `TENSOR_LL_SIGN = -1.0` since 2026-08-29 (author
  decision, plans/08 D1): the LITERATURE convention, Cosyn et al. EPJ A 61
  (2025) 83 Eq. (27) and HERMES,

      A_zz(θ_S = 0) (1 + ε(y) R) = −(2/3) b₁/F₁   exactly, at every y,

  with P_zz = n₊ + n₋ − 2n₀, i.e. b₁ > 0 means the m = 0 state has the LARGER
  cross section.  It was `+1` — the repository's own transcription of
  Hoodbhoy–Jaffe–Manohar — until that date, and setting the constant back is
  the whole of the change: nothing else in the library knows the sign.  The
  guard test (`tests/test_xsec.cpp`, "the program sign IS the literature
  sign") is written against the LITERATURE relation itself, with no reference
  to the constant, so flipping it back fails it.  What flips with it: A_zz at
  fixed b₁, the by-product κ of the spin-state ratio, and any b-sector
  subtraction built on κ — including the O(γ²) tensor leakage into cos 2φ.
  What does NOT: |A_zz|, the whole Δ (cos 2φ) sector, and `A_zz^tag(k)`
  (`azz_tensor_curve*`), which is a ratio of cluster-wave populations and
  carries no b₁ at all.
- **⁶Li effective polarization.** The CLUSTER PICTURE since 2026-08-29
  (plans/04 #6, closed).  `LI6_CLUSTER_POLARIZATION = (1 − 1.5 P_D_LI6)
  (1 − 1.5 P_D_DEUTERON) = 0.81123` whole-nucleus, and `LI6().eff_pol_p =
  eff_pol_n = LI6_CLUSTER_POLARIZATION/3`, so `Z·P_p = N·P_n = 0.81123`.  The
  two D-state probabilities live in **`beams.hpp`** — one source of truth, as
  in `polli_fastsim.beams` — and `tagged.hpp` uses those names rather than
  keeping copies, so the INCLUSIVE effective polarization and the TAGGED S/D
  interference of `li6_alpha_channel` are the same wave function seen in two
  experiments.  The deuteron slot carries `DEUTERON_VECTOR_POLARIZATION =
  1 − 1.5 P_D_DEUTERON = 0.9325` verbatim, so per-nucleon
  g₁(⁶Li)/g₁(d) = (1 − 1.5 P_D_LI6)/3 = 0.290 exactly (the deuteron's own D
  state cancels between the two isoscalar ions).  The retired Cloet
  convention stays reachable as `LI6_NAIVE_ONE_THIRD` — a whole-nucleus 1.0,
  1.233× this one and above the 0.81–0.85 band whose top is the Wiringa VMC
  0.848.
- **Exact finite-γ tensor sector.** `InclusiveKernel::Options::tensor_gamma`
  defaults to FALSE, matching `xsec.py`.  True replaces the massless
  Hoodbhoy–Jaffe–Manohar b-sector with the Cosyn Eqs. (9)/(10)/(14)/(16)/(17)/
  (24) kernel (`theta_q_cos_sin`, `cosyn_tensor_sfs`,
  `cosyn_unpolarized_sfs`, `InclusiveKernel::tensor_harmonics_gamma`), which
  carries the unmeasured `b3_func`/`b4_func` slots and leaks the rate sector
  into cos 2φ at O(γ²).  The two paths agree IDENTICALLY at γ = 0 for any b₂,
  with b₃ and b₄ cancelling, which is what makes the switch reversible; it is
  off by default because that leakage is carried as a SYSTEMATIC of the Δ
  extraction rather than as a correction the extraction subtracts, and
  because switching it on would silently move every published tensor number.
  Its size is Δ_fake/(γ² b₁) = 0.14–0.16 — the twist-4 Eq. (17e) term almost
  alone, the leading-twist T_LL and twist-3 T_LT channels standing 3 : −3 : 1
  and cancelling — which RETIRES the old bound "γ² b₁ × 1.15".
  `tests/test_tensor_gamma.cpp` is the port of
  `evgen/tests/test_tensor_gamma.py`, anchored on the two finite-γ rows of the
  paper's own Table 1 at 1e-10.

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
  `A_par = D(y) g1/F1` must construct the massless kernel explicitly.  The
  finite-γ kinematics themselves (`gamma_squared`, `epsilon_gamma`,
  `depolarization_d_gamma`, `eta_gamma`, `a_parallel_exact`,
  `depolarization_effective`) live in **`asymmetries.hpp`** — ONE
  implementation for both halves of the library, mirroring their move into
  `polli_fastsim.asymmetries`; `xsec.hpp` re-exports only the alias
  `depolarization_gamma`, exactly as `polligen.xsec` does.  `M_NUCLEON`
  (0.9383, the FREE nucleon mass) stays in `constants.hpp`, which is this
  library's single-definition rule for a convention constant.
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
  `EMC_VALENCE_DEPLETION_EPPS21` = 0.031052077003862335, because computing it
  needs LHAPDF.  Its free-nucleon denominator is **CT18ANLO**, EPPS21's own
  proton baseline, so the fit cancels and the ratio is the nuclear
  modification alone; against CT18NLO it was 0.02979, 4.2 % shallower.  On
  that baseline the two transferred camps are CBT 0.5322 and TMT 0.2113.
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
- **Cluster radial forms.** `ClusterWaveSource::Hulthen` is the DEFAULT
  everywhere and is bit-compatible with every published number; the analytic
  two-parameter forms and their `beta` band (0.20–0.40, default 0.30) are in
  `cluster.hpp`.  `ClusterWaveSource::VmcAV18` replaces them, on the two
  LITHIUM alpha-tag channels only, with the ANL VMC tables:

  | wave | magnitude |ψ_L(k)| | sign |
  |---|---|---|
  | ⁶Li α+d, L = 0 and L = 2 | `momenta/li6_ad1.momentum` S/D blocks, ψ_L = √ρ_L | `overlap_old/li6.ad`, k-space amplitudes |
  | ⁷Li α+t, L = 1 | `momenta/li7_at3.momentum` (3/2⁻ ground state) | none needed — a single wave has no interference |

  Two files per wave because they carry different things.  A momentum density
  is |ψ_L|² and has no phase, but the **S–D relative sign is observable**: the
  interference term of `n_M(k, k̂)` goes as ψ₀ψ₂ and, at P_D ≈ 0.019, it
  DOMINATES the tensor asymmetry (negating the D table flips A_zz^tag from
  +0.45 to −0.52 at k = 0.20 GeV).  The signed `overlap_old` amplitudes supply
  it.  `vmc_from_momentum` does not copy sign(A_L(k)) point by point — it
  reduces the reference to its zero CROSSINGS below 3 fm⁻¹ (⁶Li: one S node at
  0.678 fm⁻¹ = 0.134 GeV, one D node at 2.25 fm⁻¹ = 0.444 GeV, both confirmed
  bin-for-bin by minima of the momentum file's own ρ_L) and anchors the phase
  at the reference's largest |A|, where its Monte Carlo sign is beyond doubt.
  Above ~3 fm⁻¹ both columns are at the noise floor and carry ~1e-4 of the
  norm, so a noise-driven flip there would be all cost and no signal.  The
  global phase is then fixed to ψ₀(k → 0) > 0, which is unobservable but makes
  sign(ψ₂) read as the relative sign.
  On the VMC path `beta` and `p_d` are IGNORED: the shape is the table's and
  P_D is a property of the wave function (`VMC_P_D_LI6` = 0.01935, the file's
  own 0.015861/(0.80362+0.015861), against the 0.0867 SCENARIO placeholder).
  Tables are ZERO outside 5 fm⁻¹ = 0.9866 GeV, so the sampler cannot draw a
  spectator past it.  The deuteron control channel is always Hulthen — the
  deuteron IS the cluster and no d → p+n two-cluster table exists — which is
  what keeps the Cosyn–Weiss tensor gate on the analytic path.
  Reconciliation, provenance and every number:
  `docs/open_items/vmc_reconciliation.md`, `data/vmc/README.md`.
- **Triton spectral function.** `TritonSfChoice::Hulthen` is the DEFAULT
  everywhere and keeps the sequential two-body triton decay of `breakup.hpp`
  bit-compatible with every published ⁷Li number.  `TritonSfChoice::CiofiSimula`
  (CLI `--triton-sf ciofi-simula`) replaces it, on the ⁷Li α tag only, with
  the Ciofi degli Atti–Simula spectral function (`triton_sf.hpp`).
  Provenance, all transcribed and none refit: n₀(k) is CS PRC 53 (1996) 1689
  Eq. (74) with Table A.1 (A = 2, 3, 4 carried as data — the same numbers
  BeAGLE's `DT_KFERMI` holds, checked column for column against the paper);
  n₁(k) is Eq. (76) for A = 3 and Eq. (75)/Table A.3 for A = 4; k in fm⁻¹ is
  converted at ħc = `HBARC_GEV_FM` = 0.19733 GeV·fm — a `constants.hpp`
  single definition since 2026-09-01, with `coherent.hpp`'s `GEV_PER_FM_INV`
  now an alias of it.  S₀ = ∫dk k²n₀ = 0.6525 is the ³He(e,e′p)d
  spectroscopic factor (the literature's ~2/3), S₀ + S₁ closes on 1 to
  3.4 × 10⁻⁴ untuned, and the 2-body/3-body branching is the RATIO
  n₀/(n₀ + n₁) at the event's own k, never an assumed constant.  Table A.1's
  A = 3 column is the ³He **proton** distribution; under the isospin mirror
  ³He ↔ ³H it is the triton's struck-**neutron** distribution and the bound
  remnant of the mirror is the same deuteron (BeAGLE branches on A alone with
  no mirror argument).  The continuum pair is split at its own virtual-state
  pole: nn at `KAPPA_NN_VIRTUAL` = ħc/|a_nn| = 10.44 MeV (a_nn = −18.9 fm),
  pn at `KAPPA_PN_SINGLET` = ħc/|a_pn| = 8.31 MeV (a_pn = −23.74 fm, the
  ¹S₀ channel — the ³S₁ pn channel IS the bound deuteron and is the other
  branch); both constants live in `triton_sf.hpp`, one definition each.
  **The BeAGLE n₀-only caveat**: `DT_KFERMI` renormalizes its A = 3 and A = 4
  tables — the n₀ piece ALONE — to unity, silently dropping the 34.7 %
  (A = 3) / 20.0 % (A = 4) correlated continuum: its high-k tails there are
  too soft by construction and it has no three-body breakup channel to put
  that strength in.  LiPolGen does not renormalize n₀; it adds n₁ and lets
  the deficit BE the branching.  `sample()` consumes a fixed 6 uniforms on
  every branch (the breakup branch a fixed 9), so the stream position never
  depends on which channel came out.  The `Pipeline` constructs the
  `CiofiSimulaTriton` itself at the run's own `cluster_beta`, so the
  continuum pair's q-shape and every other radial form share one β.
  Full derivation and every measured number: the `triton_sf.hpp` header and
  `tests/test_triton_sf.cpp`.
- **Spectator FSI.** `PipelineConfig::fsi = PipelineFsi::Off` is the DEFAULT:
  every tagged number is the plane-wave impulse approximation bit for bit
  unless a run turns the Glauber weight on.  Three rules when it is on
  (`fsi.hpp`, `tests/test_fsi.cpp`):
  (1) **FSI is a WEIGHT, never a shift** — the rescattering multiplies
  `Event::weight` (the FSI/IA density ratio at the drawn (k, cos θ_k)) and
  moves NO four-vector; the spectator is the measurement, and the tagged
  record of an FSI-on run is bit-identical to the FSI-off run.
  (2) **The weight is SPIN INDEPENDENT by construction** — numerator and
  denominator are m-summed, so it is the unpolarized-shape distortion, which
  is what the literature supports (Cosyn–Weiss VI C leaves the spin
  dependence open).  Quote its effect as an **unpolarized-shape systematic of
  the tagged spectrum, never as a correction to A_zz**.
  (3) **σ_XN is a BAND, 20–40 mb, never one number** — 40 mb (the default,
  `--fsi-sigma-mb`) is the free hadron, 20 mb the formation-length end; run
  both.  The survival probability (∫w dΓ/∫dΓ ≈ 0.52 at 40 mb on the ⁶Li α
  tag) is logged by the run summary and exposed as
  `Pipeline::fsi_weight()->survival()`.  Variant (a) `GlauberCluster`
  shadows σ_Xα over the α's own profile (131.0 mb at σ_XN = 40, not
  4 × 40); variant (b) `GlauberNucleon` is the unshadowed A σ_XN = 160 mb
  single-scattering bracket — a POINTWISE low-k bracket, not an integrated
  one (its survival lands above (a)'s; `fsi.hpp` header).  The weight table
  is built on a (k_z, k_T) grid — a deliberate deviation from the design
  note's literal (k, cos θ_k): the eikonal kernel transfers k_T only, so the
  table is smooth and even in k_z there — and read per event by bilinear
  interpolation, consuming no randomness, so stream discipline and
  event-index determinism are untouched.
- **Where data files live at run time.** `data_dir()` (`cluster.hpp`) is
  `$LIPOLGEN_DATA_DIR` when set and non-empty, else the compiled-in
  `LIPOLGEN_DATA_DIR_DEFAULT`, which CMake sets to `${CMAKE_SOURCE_DIR}/data`.
  Nothing resolves a data path against the working directory.
- **Luminosity shares.** `Optics::lumi_fraction` multiplies COUNTS and never
  cross sections (`PipelineConfig::apply_optics_lumi_fraction`, default true) —
  the same share rule that keeps `Scenario::run_share` out of
  `sigma_per_category_pb()`.
