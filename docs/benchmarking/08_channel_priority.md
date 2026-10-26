<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# 08 — Channel priority: physics weight → benchmarking and development order

*Written 2026-10-26 from the benchmark report as it stood on 2026-09-27
(`REPORT.md`: 18 harnesses — 6 PASS, 8 FAIL recorded, 3 BLOCKED, 1 blocked
(environment)), the plan (`BENCHMARK_PLAN.md`), the surveys `00`–`07`,
`../PHYSICS_CHANNELS.md`, the decision registry
(`../open_items/run_2026-09-03/STATUS.md`, 27 rows) and the sibling
programme's plans (`PolarizedLithiumSim/plans/01`, `05`, `06`, `07`). It is a
LIVING ranking: §6 names what moves a row and §7 dates every revision. A row
here is a priority, not a decision — no default moves because of this file.*

## 0. The verdict in one paragraph

The generator exists to produce three observables for one programme: the
double-helicity-flip gluon distribution Δ(x, Q²) through the cos 2φ modulation
on a transversely tensor-polarized ⁶Li beam (the PLB letter's claim 1,
`plans/07` §7.0), the same modulation in **coherent** diffraction with the ⁶Li
detected intact (claim 2, the small-quadrupole null test), and the polarized
EMC effect in g₁(⁷Li) (paper 3) — with the spectator-tagged channels as the
programme's novelty (paper 2) and as the validation vehicle for the A = 2
machinery everything tensor is anchored on. Eighteen harnesses later the
picture is sharp: **the bookkeeping, the Monte-Carlo closures and the event
record are sound** (EPIOS, EST, weighted/unweighted, PYTHIA closure, HepMC3
convention, r_ch, NMC ⁶Li/D — every PASS is one of these), **every data
comparison of a shipped default physics input fails** (ToyF2 against NMC
F₂ᵈ at χ²/ndf = 8307.9/155; ToyG1 against world g₁ᵈ at 107.8/15, 69.4/15,
54.9/28 and a Bjorken integral that does not exist; the Hulthén deuteron's
η at +17.2σ; the ⁶Li C0 zero +0.27 fm⁻¹ above the data band; the Deeps
spectator tail at p = 5e−49 for every wave; HERMES b₁ᵈ matched in every bin
by no configuration), and **nothing validates the tensor observable itself**
(`BENCHMARK_PLAN.md` §6, U-1 … U-12 — no polarized e + Li data exist and no
other Monte Carlo does a tensor-polarized nucleus). The priority that follows
is therefore not "more rows". It is, in order: (a) price and then replace the
defaults the two top channels ride on — the unpolarized nucleon inputs F₂ and
R first, the coherent scenario's `eps_b0` second; (b) close the three rows
that were BLOCKED on a document, now that the documents are reachable
(`../open_items/run_2026-10-26/PLAN.md` phase B); (c) add the few gates that
bear on the letter's own observable — positivity of the b₃/b₄/Δ scenario
shapes, the gluonic Soffer bound on Δ, and the fill model behind the two-fill
estimator; and (d) keep the tagged sector's A = 2 anchors honest (η_d, the
Deeps tail, BONuS, d(e,e′p)n), because that is where b₁(⁶Li) is built from.

## 1. The criteria, and their weights

| id | criterion | weight | what scores 3 | what scores 0 |
|---|---|---|---|---|
| **C1** | leverage on a deliverable | 3 | a headline number of the letter (`plans/07` claims 1–3, window 22 Mar – 2 Apr 2027) | no paper names it |
| **C2** | statistical reach at the EIC | 2 | ≥ 5σ per bin in year 1 (10 fb⁻¹/u), or ε × D² ≳ 0.5 on a tagged channel | < 1σ per bin in year 1, or ε × D² < 0.01 |
| **C3** | known-wrong default or unpriced systematic | 3 | the shipped default moves the channel's observable by a measured factor ≥ 20 %, or a systematic the channel needs has no calculation at all | every input is anchored and every systematic priced |
| **C4** | benchmarkability | 1 | an external anchor exists, is machine-readable and costs ≤ 1 day to wire | nothing external touches the channel (`06_critic.md` §2.1) |
| **C5** | schedule | 1 | the isotope the source commissions first (⁷Li), or a number M2 (Jan 2027) freezes | needs Phase 2 full simulation (M3) |

Score = Σ weight × points, maximum 30. C3 carries the same weight as C1 on
purpose: a channel that feeds a headline number through an input known to be
wrong by ×1.25 is the one where benchmarking buys the most, and that is the
reading the author asked for — significant *and* uncertain goes first.

## 2. The ranking

| rank | channel / sector | C1 | C2 | C3 | C4 | C5 | score | the one-line reason |
|---|---|---|---|---|---|---|---|---|
| **1** | **Inclusive ⁶Li, tensor sector** — cos 2φ (Δ), A_zz / b₁ | 3 | 3 | 3 | 2 | 3 | **29** | the letter's claim 1 at 21–43σ per bin in year 1 (32–64σ on the two-fill estimator), riding on F₁ from a toy F₂ (σ ×0.80, A_zz ×1.25 against CT18NLO), a placeholder R (K/D_φ +26 % at x = 0.10), an ungated Δ model, a b₁(⁶Li) with a mandatory ±100 % band, and no radiative-correction calculation for any φ-dependent tensor observable |
| **2** | **Coherent ⁶Li** — a₂, cos 2(φ_t − φ_S) with the recoil intact | 3 | 2 | 3 | 1 | 2 | **25** | claim 2 and the far-forward requirement of claim 3; N_tag ≈ 1.7 × 10⁷ per year at the 0.20 GeV envelope with best-bin δA 1.8 × 10⁻³; the shipped `eps_b0` = −0.08 implies a ⁶Li quadrupole **11.4×** the measured one (registry row 8), `f0` is anchored by nothing, no Good–Walker amplitude exists for A = 6, and eSTARlight is an exclusive-VM rate scale, not a check |
| **3** | **Tagged ⁷Li α** — quasi-free polarized triton, A_∥; polarimetry ⟨P₂⟩ = −T/5 | 2 | 3 | 2 | 1 | 3 | **22** | the best tagged figure of merit by far (ε = 0.97–0.99 at IP6, ε × D² = 0.72–0.74 — `plans/05` §5.E), the isotope commissioned first, the ECRP's purity deliverable; rides on a sequential triton model by default (Ciofi–Simula opt-in, S₀ untuned), no A = 3 tensor structure, and an incoherent-background purity that is FLUKA-gated |
| **4** | **Inclusive ⁷Li** — g₁, the polarized EMC effect | 2 | 1 | 3 | 2 | 2 | **21** | paper 3's flagship, but 0.45σ per bin at 10 fb⁻¹/u in the discriminating window (1.43σ at 100); the toy g₁ⁿ has the wrong sign over 0.25 < x < 0.6 and is not integrable at small x; the polarized-EMC transfer curves are **not reachable from event generation at all**; the rank-2 sector is identically zero |
| **5** | **Tagged ⁶Li α** — embedded polarized deuteron, tagged A_zz | 2 | 1 | 2 | 2 | 1 | **17** | the programme's novelty and the object b₁(⁶Li) is built from, but beam-blind at IP6: ε = 0.017–0.026 at the YR high-acceptance optics (ε × D² = 0.0016–0.0024), 0.22–0.32 only at the tagging optics; the default deuteron's η is +17σ off, `--cluster-wave vmc` moves every polarized observable −2.03 %, the FSI weight is spin-blind |
| **6** | **Tagged d–p control** — the Cosyn–Weiss limit | 1 | 2 | 2 | 3 | 1 | **17** | not a lithium deliverable (and the EIC's polarized-deuteron option is what ⁶Li replaces), but the channel every A = 2 anchor runs on — CW Eq. (6.12) to 8.9e−16, HERMES, Deeps, BONuS, d(e,e′p)n, t₂₀ — and the one where the toy g₁ⁿ's sign bites hardest (a free neutron is the target); ranked below row 5 because it is a control |

The cross-cutting ingredients are ranked in §4 by their leverage on rows 1–3.

## 3. Channel dossiers

Each dossier: what the channel buys, what it rides on, what validates it
today (with the harness row and its status), what is missing, and the
development order that follows. Numbers are quoted with the record they come
from; none was re-measured for this file.

### 3.1 Inclusive ⁶Li, tensor sector (rank 1)

- **Buys.** The letter: xΔ(x, Q²) extracted from the cos 2φ′ modulation at
  four sweet-spot super-bins (x ≈ 0.011–0.14, Q² = 1.1–14 GeV²), amplitudes
  (0.44–0.95) × 10⁻² against δA = (1.4–4.5) × 10⁻⁴ in year 1 (`plans/07` §7.1,
  re-derived 2026-09-15 on the two-fill estimator, factor 0.667 on every
  error). Also A_zz and b₁(⁶Li), paper 2's inclusive half.
- **Rides on.** F₁ = F₂/(2x(1 + R)) on the shipped `ToyF2` and `r_sigma_lt`;
  the toy Δ (`toy_delta_gluon`, scale 10⁻³, the `moment_A` registry on the
  sibling side); `Li6B1(MillerB1)` with the ×0.5 Miller normalisation
  (registry row 2) and the rank-2 transfer 0.921949 × 2/6; the fill (the
  `ladder` EST fill, which cannot produce P_zz < 0, so the two-fill estimator's
  m = 0-rich fill at P_zz = −1.2 exists only under `--pzz-mode typed`); the
  exact finite-γ tensor kernel off by default (leakage Δ_fake/(γ²b₁) =
  0.14–0.16 carried as a systematic); the RC band on A_zz (deuteron-calibrated,
  low-x anchor row 17).
- **Validated today.** The kernel against `polligen` at rtol 1e−12 (internal);
  the relation A_zz = −(2/3) b₁/F₁ (a convention gate, U-9); `t5_est_identity`
  PASS; `t5_epios_source_modes` PASS (an atomic-beam deuteron source, mode 0
  of 6 for the EST fill); `t5_weighted_unweighted` PASS; `t5_sum_rules`
  BC closure PASS. **Against data: `t1_nmc_f2d` FAIL (8307.9/155),
  `t1_nmc_f2d_over_f2p` FAIL (2161.4/211; MSTW passes at p = 0.044),
  `t1_g1d_world` FAIL, `t2_hermes_b1_table2` FAIL (Miller ×0.5 at χ²/6 = 5.26
  is the closest; the convolution camp is indistinguishable from b₁ = 0).**
  The CT18NLO and NNPDFpol11 sub-rows have **never been measured** (blocked
  (environment) on 2026-09-26; this run's phase A).
- **Missing.** Any gate on Δ (the gluonic Soffer bound is the only one that
  exists, `04_theory.md` T-18); positivity of the b₃/b₄/Δ scenario shapes
  (T-13); a measured R inside the EIC window (H1 F_L, `02` §4 R-5) and at the
  fixed-target anchor (R1990, E140); any RC statement for a φ-dependent tensor
  observable (U-6 — not even a placeholder exists); b₁(⁶Li) beyond the α–d
  picture (U-1).
- **Development order.** (1) Measure CT18NLO and NNPDFpol11 against the T1
  data and price the default switch for registry decision — the toy costs
  ×0.80 on σ(⁶Li) and ×1.25 on A_zz, so the letter's amplitudes move with it;
  (2) wire R data and price `r_sigma_lt` against `r1998` across the window;
  (3) the positivity and Soffer gates; (4) state the RC on cos 2φ as an
  unpriced systematic in the letter's own words until a calculation exists;
  (5) resolve registry row 2 by the HERMES price (data do not decide it:
  Δχ² = +0.035).

### 3.2 Coherent ⁶Li (rank 2)

- **Buys.** Claim 2: the a₂(|t|) decomposition into a nuclear-shape term and
  a gluon-transversity term, with ⁶Li's small quadrupole (Q = −0.0818 fm², 3.5×
  smaller than the deuteron's and of the opposite sign) as a null test; the
  O5 reach, S = 2.63σ at one EIC year at the band's low edge, 3σ at
  8.3–13.0 fb⁻¹/u. Claim 3's far-forward requirement (no published IP6
  optics delivers the 0.20 GeV envelope; `plans/10`).
- **Rides on.** A SCENARIO with bands: `f0` = 0.04 (bounds the repository
  does not cite), `slope_b` = 50 ± 10 (the one anchored number), `amp`,
  `eps_b0` = −0.08 (a deuteron-scaled number implying Q_charge(⁶Li) =
  −0.9345 fm², 11.42× measured — kept as a reference gate, registry row 8;
  the honest ⁶Li band is −(0.0070 … 0.0527)); M_X ≥ 1.2 GeV with no
  exclusive-VM channel below it; x_P sampling; the γ*–Pomeron T2 tier
  (`PDF:PomSet` moves nothing at T0).
- **Validated today.** `tests/test_coherent.cpp` closed forms (⟨|t|⟩ = 1/B,
  a₀ = −2a₁); the port gate; `t6_hepmc3_convention` PASS on the coherent
  record; `t4_pythia_ep_closure` PASS for the hadronizer. **Against anything
  external: nothing** (U-3). eSTARlight was retracted as a benchmark
  (`BENCHMARK_PLAN.md` §5.3): 85–97 % of its rate is at Q² < 0.1 GeV².
- **Missing.** A dipole-model (Good–Walker) amplitude for polarized ⁶Li —
  `cluster_config.hpp` writes the configurations such a run averages over,
  and the Mäntysaari-group letter that would start it is unsent (row 21);
  an anchor for `f0`; the radiative correction on the recoil-tagged
  cross section (the RC family is refused on this channel by design).
- **Development order.** (1) Decide `eps_b0` at M2 (row 8): every generated
  coherent tensor number is 11.4× the measured-quadrupole expectation, and the
  letter quotes a₂; the cost is one re-pin of `coherent.json`; (2) send the
  dipole-model request (row 21); (3) price `f0`'s band against the only
  diffractive-continuum estimate available (the eSTARlight lower bound is 40×
  below the default and is a lower bound, not a value); (4) the exclusive-VM
  window below M_X = 1.2 GeV as a Phase 2 item.

### 3.3 Tagged ⁷Li α (rank 3)

- **Buys.** Paper 3's tagged half — a quasi-free polarized triton (P_p =
  +0.866, P_n = −0.037) as an in-medium polarized proton, A_∥ with the
  tagging purity the ECRP proposal explicitly defers to simulation; the
  polarimetry identity ⟨P₂⟩ = −T/5; the first isotope the source delivers.
- **Rides on.** The pure L = 1 α–t wave (Hulthén P by default, the ANL VMC
  α+t table under `--cluster-wave vmc`); the triton breakup — a sequential
  two-body decay at AME2020 separation energies by default, the three-channel
  Ciofi–Simula spectral function opt-in (S₀ = 0.6525, untuned, isoscalar
  approximation for the struck proton); g₁ᵖ from `ToyG1` (`--pol-sf nnpdfpol`
  opt-in); the RC band reached only at `--pe 0`.
- **Validated today.** `tests/test_tagged.cpp`: P_p = 0.866 and the forward
  limit, ⟨P₂⟩ = −T/5 (internal identities); `t6_hepmc3_convention` PASS on the
  ⁷Li-α record. **Against data: nothing at A = 3 or A = 7.** The Ciofi–Simula
  parameterization is a fit to a many-body calculation and a ³He column stands
  in for ³H (`02` §6).
- **Missing.** Any A = 3 tagged datum; the A = 3 tensor/D-wave structure (both
  triton branches isotropic); the FLUKA-gated incoherent-background purity
  (`plans/05` §5.E — left empty rather than guessed); an A = 3 wave function
  with MC errors (the AV18 A = 3 overlap does not exist in-tree).
- **Development order.** (1) g₁ backend: on this channel the struck proton
  carries the signal, so the toy g₁ᵖ's shape is the systematic — measure
  NNPDFpol11 against E143/HERMES/COMPASS g₁ᵖ (the T1 harness compares g₁ᵈ;
  g₁ᵖ is one vendoring away); (2) the triton SF default (Ciofi–Simula vs
  sequential) priced on the tagged A_∥ and the spectator spectrum; (3) the
  polarized-EMC transfer as a reachable kernel option (today `medium_ratio` is
  empty in every default kernel), since paper 3's observable is exactly that
  ratio.

### 3.4 Inclusive ⁷Li (rank 4)

- **Buys.** Paper 3's inclusive observable ΔR_A(x) = g₁ᴬ/(P_p g₁ᵖ + P_n g₁ⁿ):
  CBT against TMT, a separation of 0.021 at x = 0.36 and 0.006 at 0.65 on the
  transferred curves, 0.45σ per bin at 10 fb⁻¹/u (`plans/01` §2.3).
- **Rides on.** `ToyG1` (`a1n = −0.07(1 − x)² + 0.8x^2.2`, whose g₁ⁿ ~ x^−1.22
  below x ≈ 10⁻³ is not integrable), the effective polarizations (VMC), a
  rank-2 sector that is exactly zero and says so (`meta["rank2_input"]`), and
  **no** polarized-medium modification: `cbt_`/`tmt_polarized_emc_ratio` exist
  but no kernel passes a `medium_ratio` (`PHYSICS_CHANNELS.md` §13).
- **Validated today.** `t1_g1d_world` FAIL on both backends measured
  (ToyG1 and ToyG1-on-MSTW-F₁); NNPDFpol11 never measured; the Bjorken
  integral FAIL on ToyG1; the rank-2 zero is pinned (`tests/test_pipeline.cpp`
  D1).
- **Missing.** The reachable polarized-EMC ratio; a ⁷Li rank-2 structure
  function (U-4, no calculation exists anywhere; deferred on row 20); a g₁ᵖ
  gate; any ⁷Li DIS datum (none exists, `03` §8).
- **Development order.** (1) the NNPDFpol11 measurement (phase A here);
  (2) the `medium_ratio` hook made reachable and priced — without it the
  generator cannot produce paper 3's observable at all; (3) the rank-2 sector
  stays zero until row 20 is decided.

### 3.5 Tagged ⁶Li α (rank 5)

- **Buys.** Paper 2: the tagged tensor asymmetry on the embedded deuteron,
  O(1) at k ≈ 0.3 GeV from D-wave dominance (A_zz^tag(0.20 GeV) = −1.2069 on
  Hulthén, −0.5191 on AV18, acceptance-weighted at the YR high-acceptance
  optics); the α–d D-wave itself, which `LI6_B1_RANK2_TRANSFER` = 0.921949 and
  the inclusive b₁(⁶Li) are built from.
- **Rides on.** The α–d S + D wave (Hulthén by default, VMC AV18 opt-in, the
  two differing by 6.6 % on the tagged tensor observable and by +11.61 % /
  +6.58 % between the tagged and inclusive dilutions under `vmc`); the
  embedded deuteron (same flag since 2026-09-04); the T1 breakup; the
  spin-blind Glauber FSI weight; the i^L phase fixed on 2026-09-06 (row 26).
- **Validated today.** The Cosyn–Weiss gate as an identity on the deuteron
  control (8.88e−16 over 26 880 cells); the port gate; **`t2_deuteron_static`
  FAIL for the shipped Hulthén (η = 0.0325 against 0.0256(4), +17.2σ;
  decision D2 keeps it and quotes physics on AV18); `t2_deeps_spectator_tail`
  FAIL for every wave and both prescriptions (CD-Bonn/LC 278.8/187 is the
  best); `t2_bonus_spectator_shape` BLOCKED (unblocked this run)**.
- **Missing.** Anything with an α spectator (U-2: no calculation exists); the
  (e,e′d) / (e,e′α) α–d momentum distributions (Ent, Mitchell — paywalled,
  digitizable); the D wave from data (U-8); an FSI with spin dependence.
- **Development order.** (1) Understand the Deeps failure — the harness codes
  its own light-cone map and leaves the x* off-shell drift uncorrected, so the
  verdict is a joint test of n(k), the map and FSI; separate them before
  reading it as a wave-function verdict; (2) BONuS as the low-p_s complement;
  (3) the D2 revisit at M2 (switch the default deuteron to AV18 if any quoted
  number still rides on the Hulthén D wave); (4) IP6 optics: the physics case
  for this channel is the IR-8 secondary focus, outside this generator.

### 3.6 Tagged d–p control (rank 6)

- **Buys.** Validation only: the Cosyn–Weiss tagged limit, the row every A = 2
  anchor runs on. Its own physics (tagged spin-1 DIS on a polarized deuteron
  beam) is what the EIC may not be able to run, which is the premise of the
  lithium programme.
- **Rides on, validated, missing.** As 3.5 at A = 2, plus: the toy g₁ⁿ's
  wrong sign over 0.25 < x < 0.6 bites hardest here (a free neutron is the
  DIS target; the pair `--unpol-sf ct18nlo|mstw --pol-sf nnpdfpol` needs
  `--x-max 0.95`); the RC elastic-tail deuteron form factor is a shape-blind HO
  stand-in while Abbott's t₂₀ / A(Q²) world data exist on HEPData (`02` §5.3).
- **Development order.** (1) t₂₀ / form-factor data into the RC control
  (phase C here); (2) d(e,e′p)n momentum distributions (Boeglin) to turn
  T25's model-vs-model into model-vs-data; (3) BONuS and Deeps as above.

## 4. Cross-cutting ingredients, ranked by leverage on ranks 1–3

| rank | ingredient | moves | status | the development item |
|---|---|---|---|---|
| **A** | **unpolarized nucleon inputs F₂, R** (`--unpol-sf`, `--r-source`) | every F₁, hence every A_zz and cos 2φ normalisation on every channel: toy → CT18NLO is σ ×0.7985 and A_zz ×1.2524 at ⁶Li config 1; `r1998` against `r_sigma_lt` is K/D_φ −3.8 / +26.4 / +3.4 % at x = 0.05 / 0.10 / 0.30 | T1: toy FAIL, MSTW FAIL/PASS, **CT18NLO unmeasured**; R: no data comparison anywhere | measure CT18NLO (phase A), wire R data (phase C), then a registry decision on the default backend at M2 |
| **B** | **the fill model** (`--pzz-mode`, EST ladder) | the letter's 1.5× two-fill gain needs an m = 0-rich fill at P_zz = −1.2, which the default EST ladder cannot represent (P_zz ≥ 0 always); EPIOS matches the EST fill at 1 of 6 modes | `t5_epios_source_modes` PASS on the round-trip; the fill question recorded, not asserted | state in the letter which fill the generator assumed; a Li source fill table does not exist (N-2) |
| **C** | **the Δ sector** (`toy_delta_gluon`, b₃/b₄ scenario shapes, `tensor_gamma`) | claim 1's extracted xΔ directly | **no gate at all**; the finite-γ leakage is a carried systematic | the gluonic Soffer bound and the spin-1 positivity bounds as T5 closures (phase C) |
| **D** | **b₁(⁶Li) model** (`--b1-model`, row 2, `LI6_B1_RANK2_TRANSFER`) | A_zz(⁶Li) by ×2 (row 2) and by the ±100 % band | HERMES does not decide row 2; the α–d step is validated by nothing (U-1) | keep the band mandatory; a theory input (Mäntysaari group / Cosyn) is the only way forward |
| **E** | **coherent scenario** (`eps_b0`, `f0`, `amp`, `slope_b`) | claim 2's a₂ by 11.4× (`eps_b0`) and N_tag by `f0`'s factor-2 band | `slope_b` anchored; the rest not | row 8 at M2; row 21 (the dipole-model request) |
| **F** | **radiative corrections** (`--rc`, `--rc-tail-model`) | A_zz by the band's half-width (row 17: the low-x anchor errs wide by ×2.65); cos 2φ by an amount **no calculation gives** (U-6); nothing on the coherent channel by design | POLRAD rows reproduce (2 of 3); the band is deuteron-calibrated; the polarized QE tail is unpriced | row 7's DJANGOH run (phase B) is the only generator-vs-generator RC check possible; the cos 2φ RC stays a stated unknown |
| **G** | **deuteron wave function default** (Hulthén) | the tagged S–D interference and η; −2.03 % on the α-tag path under `vmc` | `t2_deuteron_static` FAIL +17.2σ; D2 keeps the default until M2 | D2 revisit at M2 |
| **H** | **polarized nucleon input g₁** (`--pol-sf`) | A_∥ on ranks 3 and 4; nothing on claim 1 (unpolarized electrons) | `t1_g1d_world` FAIL; **NNPDFpol11 unmeasured** | phase A measures it |
| **I** | **spectator FSI** (`--fsi`) | the tagged spectrum's shape, never A_zz (spin-blind by construction) | unitarity not definable; the two variants differ on 99.5 % of events | the Deeps transverse window and BONuS as direction checks |
| **J** | **T2 hadronization** | reco-level FOMs on the sibling side | `t4_pythia_ep_closure` PASS (worst 0.023 against 0.20) | none for the letter; the PomSet band is hadronic only |
| **K** | **output and chain** | nothing physical | `t6_hepmc3_convention` 23/23; the ePIC legs (abconv, npsim, EICrecon, Rivet) unrun since 2026-09-02 | M3; a light-ion `abconv` preset does not exist upstream |
| **L** | **⁷Li rank-2 sector** | nothing: identically zero | pinned zero; deferred on row 20 | waits on theory input |

## 5. What the ranking does to the benchmark plan

1. **The M2 rule gains a clause.** "No benchmark moved beyond its stated
   tolerance since M1" is necessary and not sufficient: at M2 the report must
   have been regenerated on a machine with the three LHAPDF sets installed, so
   that the CT18NLO and NNPDFpol11 sub-rows carry numbers (they never have),
   and each of ranks 1–3 must carry a **default-input price table** — the
   observable's movement under every opt-in backend — so that the registry can
   decide the defaults from numbers.
2. **New rows, in this order** (wired in `../open_items/run_2026-10-26/`):
   R data (H1 F_L, R1990, E140) for ingredient A; the positivity and Soffer
   closures for C; BONuS, RFY and the DJANGOH elastic-tail run for the three
   BLOCKED rows; Abbott's t₂₀ / A(Q²) for the RC control's form factor;
   HERMES F₂ᵈ and the final COMPASS g₁ᵈ as T1 extensions; the
   Wiringa–Schiavilla ⁶Li form factors (a theory anchor for the C0 zero);
   Boeglin's d(e,e′p)n momentum distribution.
3. **Rows that stay FAIL by decision are read against this ranking.** The
   Hulthén η (rank 5/6, D2) and the C0 zero (ingredient F, ≤ 2.2e−7 on
   ΔA_zz) can stay recorded. The toy F₂ / g₁ failures (ingredient A/H) and
   `eps_b0` (ingredient E) cannot: they sit under ranks 1 and 2 and are the
   decisions this file asks the author to take at M2.
4. **What is deliberately NOT raised.** The ePIC chain legs (K): they
   validate the record, not the physics, and need Phase 2's container stack;
   ⁷Li rank-2 (L): no input exists; FSI unitarity: not definable in the tree.

## 6. Triggers that re-rank

- A change of the letter's scope (`plans/07` D1): if the coherent channel is
  dropped from the letter, rank 2 falls to rank 4 and ingredient E to the
  bottom of §4.
- An IP6 optics decision for ⁶Li tags (`plans/10`): the tagging optics lift
  rank 5's ε × D² by 10×, to 0.02–0.03, and with the IR-8 secondary focus the
  channel would join rank 3.
- A theory input: a dipole-model a₂ for polarized ⁶Li (row 21) moves
  ingredient E from "scenario" to "calculation"; a b₁(⁶Li) calculation beyond
  the α–d picture does the same for D.
- An EIC lithium luminosity or polarization number (gap #4 of `plans/01`):
  every C2 entry is at 10 fb⁻¹/u and P_zz = 0.6 by assumption.
- Data: JLab E12-13-011 (b₁ᵈ, 0.16 < x < 0.49) changes the HERMES-only state
  of row 3; E12-14-001 (polarized EMC, ⁷LiD) gives rank 4 its first datum.
- Registry decisions on rows 2, 3, 8, 9 (D2), 17 and 20: each one removes a
  C3 point from the channel it sits under.

## 7. Revision log

| date | what changed | by |
|---|---|---|
| 2026-10-26 | first version, from the 2026-09-27 report and the 2026-09-27 decisions | run `../open_items/run_2026-10-26/` |
