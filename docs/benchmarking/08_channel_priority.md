<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# 08 — Channel priority: physics weight → benchmarking and development order

*Written 2026-10-26 from the benchmark report as it stood on 2026-09-27
(`REPORT.md`: 18 harnesses — 6 PASS, 8 FAIL recorded, 3 BLOCKED, 1 blocked
(environment)) and re-measured the same day on a machine with the LHAPDF
sets (`../open_items/run_2026-10-26/phase_0_environment.md`); from the plan
(`BENCHMARK_PLAN.md`), the surveys `00`–`07`, `../PHYSICS_CHANNELS.md`, the
decision registry (`../open_items/run_2026-09-03/STATUS.md`, 27 rows) and
the sibling programme's plans (`PolarizedLithiumSim/plans/01`, `05`, `06`,
`07`, `10`). Reviewed the same day by two adversarial readers (49 findings;
§7). It is a LIVING ranking: §6 names what moves a row and §7 dates every
revision. A row here is a priority, not a decision — no default moves
because of this file, and §5.1's change to the milestone rule is a
proposal for the registry, not an edit of the plan.*

## 0. The verdict in one paragraph

The generator exists to produce three observables for one programme: the
double-helicity-flip gluon distribution Δ(x, Q²) through the cos 2φ
modulation on a transversely tensor-polarized ⁶Li beam (the PLB letter's
claim 1, `plans/07` §7.0), the same modulation in **coherent** diffraction
with the ⁶Li detected intact (claim 2, the small-quadrupole null test), and
the polarized EMC effect in g₁(⁷Li) (paper 3) — with b₁(⁶Li)/A_zz and the
spectator-tagged channels as paper 2, and the tagged d–p and inclusive d,
³He, p controls as the validation vehicles of the A = 2 machinery everything
tensor is anchored on. Eighteen harnesses later the picture is sharp. **The
bookkeeping, the Monte-Carlo closures and the event record are sound**: of
the six PASS rows of the 2026-09-27 report, four are those (EPIOS source
modes, the EST identity, weighted/unweighted, the PYTHIA closure), one is a
document conformance that passed by re-scoping the document (HepMC3, D1),
and one is the single data comparison of a shipped default that passes —
the ⁶Li charge radius of `HoSpin1FF` against de Vries (a consistency check
of an input, previewed); NMC F₂(⁶Li)/F₂(D) passed on 2026-09-23 and again
in phase 0 here (blocked by the environment on 2026-09-27). **Every other
data comparison of a shipped default physics input fails**: ToyF2 against
NMC F₂ᵈ at χ²/ndf = 8307.9/155 (and no opt-in backend passes it either —
MSTW 416.5/155, CT18NLO evaluated at leading order 3401.8/146, while the
ratio F₂ᵈ/F₂ᵖ passes for both); ToyG1 against world g₁ᵈ at 107.8/15, 69.4/15,
54.9/28 with a Bjorken integral that does not exist; the Hulthén deuteron's η
at +17.2σ; the ⁶Li C0 zero +0.27 fm⁻¹ above the data band; the Deeps
spectator tail for every wave (p = 4.9e−49 for the Hulthén default, 1.5e−5
for CD-Bonn/LC); HERMES b₁ᵈ matched in every bin by no configuration. And
**nothing validates the tensor observable itself** (`BENCHMARK_PLAN.md`
§6, U-1 … U-12 — no polarized e + Li data exist and no other Monte Carlo
does a tensor-polarized nucleus; no harness in the tree exercises the
cos 2φ amplitude at all). The priority that follows is therefore not "more
rows". It is, in order: (a) price the inputs the two top channels ride on
for the observable each one quotes — for the cos 2φ amplitude that is R and
the Δ model, not F₂; for A_zz it is F₂, R, the b₁ camp and the α–d transfer;
for the coherent a₂ it is `eps_b0`; (b) close the two rows that were
BLOCKED on a document and the one blocked on a run, now that the documents
are reachable and `gfortran` exists (`../open_items/run_2026-10-26/PLAN.md`
phase B); (c) put the one bound that exists on the Δ sector into the tree,
labelled as a bound under stated assumptions; and (d) keep the A = 2
anchors honest for what they validate — the d–p control and the opt-in
convolution path (η_d, the Deeps tail, BONuS, d(e,e′p)n), not the shipped
b₁(⁶Li), which is Miller's b₁ᵈ times an α–d transfer that nothing external
touches.

## 1. The criteria, their levels, and the three override rules

Score = Σ weight × points, maximum 30. Every level is defined so that a
second reader can reproduce the table from the dossiers' numbers.

| id | criterion | weight | 3 | 2 | 1 | 0 |
|---|---|---|---|---|---|---|
| **C1** | leverage on a deliverable | 3 | a headline number of the letter (`plans/07` claims 1–3, window 22 Mar – 2 Apr 2027) | the main observable of paper 2 or 3 (`plans/07` §7.0) | a control, or an input to those papers | no paper names it |
| **C2** | statistical reach, scored at 10 fb⁻¹/u, P_zz = 0.6 / P_z = 0.7, the published default optics (YR high-acceptance for tags) and the toy backend unless a grid number exists | 2 | ≥ 5σ per bin in year 1, or ε × D² ≥ 0.5 on a tagged channel | 1–5σ per bin, or ε × D² 0.1–0.5 | 0.5–1σ per bin, or ε × D² 0.01–0.1 | < 0.5σ per bin, or ε × D² < 0.01 |
| **C3** | known-wrong default or unpriced systematic, for THIS observable | 3 | the shipped default moves the observable by ≥ 20 % against an opt-in backend or a measured number, or the observable's central model is a scenario with no external calculation | a 5–20 % shift, or a band-mandatory input that has an anchor | shifts < 5 %, every input anchored but banded | anchored and priced |
| **C4** | benchmarkability | 1 | a machine-readable external anchor exists for the observable itself, ≤ 1 day to wire | anchors exist for its inputs | anchors only at another A or another process | nothing external touches it (`06_critic.md` §2.1) |
| **C5** | schedule | 1 | the isotope the source commissions first (⁷Li), or a number M2 (Jan 2027) freezes for the letter | a paper-2/3 number with its own date | waits on an external input (optics, theory, data) | needs Phase 2 full simulation (M3) |

**Override rules, applied after the scores and stated as rules.** R1 — a
*control* species (the tagged d–p channel; inclusive d, ³He, p) ranks below
every lithium channel whatever its score: the author's criterion is the
channel's contribution to the physics, and a control contributes none of its
own. R2 — ties among lithium channels are broken by C1, then C5. R3 — C2 is
scored at the published default optics; an optics that trades luminosity for
acceptance is a §6 trigger, not a score.

C3 carries the same weight as C1 — this file's choice, not the author's,
who asked that *significant* channels go first: a channel that feeds a
headline number through an input known to be wrong is where benchmarking
buys the most. The order of the top four is unchanged with C3 weighted 1.

## 2. The ranking

| rank | channel / sector | C1 | C2 | C3 | C4 | C5 | score | the one-line reason |
|---|---|---|---|---|---|---|---|---|
| **1** | **Inclusive ⁶Li, cos 2φ — Δ(x, Q²)** (claim 1) | 3 | 3 | 3 | 1 | 3 | **28** | 18–28σ per sweet-spot bin in year 1 on the grid backend (21–43σ on the toy; `plans/07` §7.1); Δ is a scenario (`toy_delta_gluon` at the shipped scale 10⁻², no external calculation, U-5); R moves the amplitude −15.6 / −15.1 / −7.1 / +4.6 % at the four sweet-spot cells x = 0.011 / 0.028 / 0.071 / 0.141 generator-side (`r1998` vs `r_sigma_lt`; −15.2 / −14.7 / −7.8 / +4.2 % over the generated super-bins; the sibling's whole chain 0.842 / 0.849 / 0.930 / 1.039 — *2026-10-26:* the "−3 … −8 %" written here until then were function-level points, not the sweet spots; `../open_items/run_2026-10-26/phase_C1_r_world.md` §4, `phase_A_baseline.md` §4) and F₂ cancels in it exactly; no anchor beyond a moment-level bound under assumptions |
| **2** | **Coherent ⁶Li — a₂, cos 2(φ_t − φ_S)** with the recoil intact (claim 2) | 3 | 2 | 3 | 0 | 3 | **25** | N_tag ≈ 1.7 × 10⁷ per year at the 0.20 GeV envelope with best-bin δA 1.8 × 10⁻³ (1.2 × 10⁻³ two-fill): at the *measured* quadrupole the deformation a₂ (≈ 0.006 at the tag ⟨\|t\|⟩) is a 3–4σ-class signal and the gluon term 3 × 10⁻³–10⁻² a 2–6σ one; the shipped `eps_b0` implies a quadrupole **11.4×** the measured one (row 8), `f0` is anchored by nothing, no Good–Walker amplitude exists for A = 6 (U-3). *The coherent J/ψ a₂ reach (O5, 2.63σ, a MARGINAL band whose top is not established to cross 3σ) is another process, not generated here.* |
| **3** | **Tagged ⁷Li α** — A_∥ on the quasi-free polarized triton; the polarimetry identity ⟨P₂⟩ = −T/5 (paper 2) | 2 | 3 | 3 | 1 | 3 | **25** | the best tagged figure of merit (ε = 0.97–0.99 at IP6, ε × D² = 0.72–0.74, `plans/05` §5.E); the isotope commissioned first; the triton breakup is a sequential model by default (Ciofi–Simula opt-in, S₀ untuned), the P_n gate is OPEN, and the incoherent-background purity is FLUKA-gated — no calculation |
| **4** | **Inclusive ⁶Li, A_zz / b₁(⁶Li)** (paper 2's inclusive half) | 2 | 2 | 3 | 2 | 2 | **23** | the sibling's A_zz money plot exists but its per-bin significance was not read for this file (scored 2, UNSCORED from a tree number); toy → CT18NLO is ×1.2524 on the window-level A_zz (0.86–1.44 per x bin), Miller's normalisation is a factor 2 (row 2), the α–d transfer is validated by nothing (U-1) and the ±100 % band is mandatory; HERMES b₁ᵈ and the A = 2 calculations anchor the input |
| **5** | **Inclusive ⁷Li — g₁, the polarized EMC effect** (paper 3) | 2 | 0 | 3 | 2 | 3 | **20** | 0.45σ per bin at 10 fb⁻¹/u in the discriminating window, nothing reaches 1σ (`plans/01` §2.3); ToyG1's g₁ⁿ has the wrong sign over 0.25 < x < 0.6 and is not integrable; the polarized-EMC transfer curves are **not reachable from event generation**; the rank-3 (octupole) contamination of A_∥ by the default ladder fill is unpriced and has no structure function |
| **6** | **Tagged ⁶Li α** — embedded polarized deuteron, tagged A_zz (paper 2) | 2 | 0 | 3 | 2 | 1 | **18** | the programme's novelty, beam-blind at IP6: ε = 0.017–0.026 at the YR optics, ε × D² = 0.0016–0.0024; the wave-function default is the largest measured lever on the observable — A_zz^tag(0.20 GeV) = −1.2069 (Hulthén) against −0.5191 (AV18), ×4.8 at 0.30 GeV, opposite signs by 0.50 GeV; waits on the IR-8 secondary focus |
| **7** | **Tagged d–p control** — the Cosyn–Weiss limit (R1: a control) | 1 | 3 | 2 | 3 | 1 | 19 → last | ε × D² = 0.88 *if* the EIC runs a polarized deuteron, which is the premise ⁶Li replaces; the channel the spectator anchors run on — CW Eq. (6.12) to 8.9e−16, Deeps, BONuS, d(e,e′p)n; the default Hulthén's η is +17σ off and `--cluster-wave vmc` moves its vector/tensor dilutions −2.03 % / −1.18 %; the toy g₁ⁿ's sign bites hardest here |
| — | **Inclusive control species d, ³He, p** (R1: controls, not scored) | | | | | | | the vehicles of ingredients A and H: HERMES b₁ᵈ, the T1 F₂ᵈ / g₁ᵈ rows and t₂₀ are *inclusive* or elastic deuteron objects; ³He (P_n = 0.86) is the direct probe of the toy g₁ⁿ's sign (HERMES ³He g₁ⁿ, E06-014 at 0.25 ≤ x ≤ 0.90, `02` §6) |

## 3. Channel dossiers

Each dossier: what the channel buys, what it rides on, what validates it
today (harness row and status), what is missing, and the development order
that follows. Numbers carry the record they come from; none was re-measured
for this file.

### 3.1 Inclusive ⁶Li, cos 2φ — Δ (rank 1)

- **Buys.** The letter: xΔ(x, Q²) extracted from the cos 2φ′ modulation at
  four sweet-spot super-bins (x = 0.011 / 0.028 / 0.071 / 0.141, Q² ≈
  1.1–14 GeV²), amplitudes (0.44–0.95) × 10⁻² on the toy backend,
  (0.30–0.82) × 10⁻² on the grid (`plans/07` WP1), against δA =
  (1.4–4.5) × 10⁻⁴ single-fill, (0.95–3.0) × 10⁻⁴ two-fill (which footing the
  letter quotes is the author's call, `plans/07` §7.1).
- **Rides on.** The shipped `toy_delta_gluon` at scale 10⁻² (Δ/F₁ fixed by
  the scale; the letter quotes the sibling's `moment_A` model instead);
  R through D_φ (the amplitude is Δ/D_φ, F₂ cancels exactly); the exact
  finite-γ tensor kernel off by default (leakage Δ_fake/(γ²b₁) = 0.14–0.16
  carried as a systematic); the two-fill pair, which the generator builds by
  default (`tensor-flip` is (+P_zz, −2P_zz); `tensor-thirds` carries the
  −2P_zz category) — whether a lithium source delivers P_zz = −1.2 at
  matched purity is an external input (`plans/07` §7.1; no EPIOS mode
  reaches below −1). *2026-10-26:* at the four super-bins (x = 0.011 /
  0.028 / 0.071 / 0.141) the shipped Δ gives A₂ = −1.06 / −1.30 / −1.43 /
  −1.33 × 10⁻³ per unit P_zz — the opposite sign of the letter's `moment_A`
  amplitudes (+4.35 / +7.42 / +9.49 / +9.52 × 10⁻³) and 4–7× smaller
  (`../open_items/run_2026-10-26/phase_A_baseline.md` §4.3).
- **Validated today.** Only the EPIOS round trip of `spin1_populations`
  touches the tensor plans' fill, and U-9 is a convention gate on a
  relation. **No harness exercises the cos 2φ amplitude** *(2026-10-26: two now
  gate its inputs — `t1_r_world`, row 12, R against NMC's R, PASS; and
  `t5_delta_moment_bound`, row 13, Δ's moment against the gluonic Soffer
  bound, PASS for the shipped toy — and the first prices R's lever on it;
  the amplitude itself still has no external anchor)*; `t5_est_identity`
  concerns the helicity-flip fill, `t5_weighted_unweighted` is a
  φ-integrated closure at θ_S = 0 where the Δ term vanishes, and g₁ does not
  enter this observable (unpolarized electrons).
- **Missing.** Any gate on Δ beyond the moment-level gluonic Soffer-type
  bound, wired 2026-10-26 under stated assumptions (`t5_delta_moment_bound`:
  PASS for the shipped toy; the letter's scenario at 0.67–0.74 or 1.35–1.48
  of the bound depending on the factor-2 normalisation of Δ in the kernel's
  cos 2φ term, open);
  a measured R at the sweet spots (NMC's R at 0.002 < x < 0.12, `02` §4
  R-3, is the one data set that covers them; H1 F_L is below them, the SLAC
  sets above — *2026-10-26:* compared, row 12; NMC fixes (1 + R) at the
  sweet-spot x to ±7.5–9.3 %, and at higher Q² than the cells); a full O(α) polarized RC — the letter carries a leading-log
  ISR bound of ≤ 2.9 % on the amplitude (`plans/07` WP4) and
  `RcScope::TensorAll` prices the Δ cos 2φ omission generator-side, but no
  hadronic-side calculation exists (U-6).
- **Development order.** (1) R against NMC R, R1990 and E140 (the latter two
  are `r1998`'s own fit input — a closure) and H1 at low x, with the R
  price evaluated at the four sweet-spot cells; (2) the moment bound on Δ,
  for the shipped toy at 10⁻² and the letter's `moment_A`/`moment_B`;
  (3) the finite-γ kernel's leakage and the ISR bound stated in the
  letter's own words; (4) F₂ enters only through the rate (N_DIS ×0.766 on
  the grid, δA ×0.95–1.16) — a rate statement, not an amplitude one.

### 3.2 Coherent ⁶Li (rank 2)

- **Buys.** Claim 2: the a₂(|t|) decomposition into a nuclear-shape term and
  a gluon-transversity term, with ⁶Li's small quadrupole (Q = −0.0818 fm²,
  3.5× smaller than the deuteron's and of the opposite sign) as a null test;
  claim 3's far-forward requirement (no published IP6 optics delivers the
  0.20 GeV envelope, `plans/10`). The coherent J/ψ reach (O5: S = 2.63σ at
  one year at the band's low edge, 3σ at 8.3–13.0 fb⁻¹/u, answered as a
  band — MARGINAL, the top not established to cross 3σ,
  `../OPEN_ITEMS_SOLUTIONS.md` §11.3c) is exclusive vector-meson
  production, which this generator does not produce.
- **Rides on.** A SCENARIO with bands: `f0` = 0.04 (bounds the repository
  does not cite), `slope_b` = 50 ± 10 (the one anchored number), `amp`,
  `eps_b0` = −0.08 (a deuteron-scaled number implying Q_charge(⁶Li) =
  −0.9345 fm², 11.42× measured — kept as a reference gate, registry row 8;
  the honest ⁶Li band is −(0.0070 … 0.0527)); M_X ≥ 1.2 GeV with no
  exclusive-VM channel; the γ*–Pomeron T2 tier (`PDF:PomSet` moves nothing
  at T0).
- **Validated today.** `tests/test_coherent.cpp` closed forms (⟨|t|⟩ = 1/B,
  a₀ = −2a₁); the port gate; `t6_hepmc3_convention` on the coherent record;
  `t4_pythia_ep_closure` for the hadronizer. **Against anything external:
  nothing** (U-3). eSTARlight was retracted as a benchmark
  (`BENCHMARK_PLAN.md` §5.3).
- **Missing.** A dipole-model (Good–Walker) amplitude for polarized ⁶Li —
  `cluster_config.hpp` writes the configurations such a run averages over,
  and the Mäntysaari-group letter that would start it is unsent (row 21);
  an anchor for `f0`; any RC on the recoil-tagged cross section (refused on
  this channel by design).
- **Development order.** (1) Decide `eps_b0` at M2 (row 8): the deformation
  a₂ is 11.4× too large, and the generated modulation c₂ = 2a₂ + amp·P_zz
  moves by about 6× at the tag ⟨|t|⟩ under the default `amp`; the letter
  quotes a₂; the cost is one re-pin of `coherent.json`; (2) send the
  dipole-model request (row 21); (3) price `f0`'s band (the eSTARlight
  lower bound is 40× below the default and is a bound, not a value);
  (4) the exclusive-VM window as a Phase 2 item, with J/ψ above M_X = 1.2 GeV
  named as a separate process.

### 3.3 Tagged ⁷Li α (rank 3)

- **Buys.** Paper 2's tagged half on the first isotope the source delivers:
  A_∥ on a quasi-free polarized triton (the generator's triton: P_p = +0.86,
  P_n = −0.028, the Bissey ³He mirror; the ⁷Li whole-nucleus VMC sums are
  +0.866 / −0.037), the ECRP's purity deliverable, the polarimetry identity
  ⟨P₂⟩ = −T/5.
- **Rides on.** The pure L = 1 α–t wave (Hulthén P by default, the ANL VMC
  α+t table under `--cluster-wave vmc`); the triton breakup — a sequential
  two-body decay by default, the three-channel Ciofi–Simula spectral
  function opt-in (S₀ = 0.6525, untuned, isoscalar approximation for the
  struck proton); g₁ᵖ from `ToyG1` (`--pol-sf nnpdfpol` opt-in); the RC band
  reached only at `--pe 0`.
- **Validated today.** `tests/test_tagged.cpp`: P_p within ±0.02 of 0.866
  (met) and the forward limit; **the P_n gate (−0.037) is OPEN** — the model
  does not reach it; ⟨P₂⟩ = −T/5 (internal identity);
  `t6_hepmc3_convention` on the ⁷Li-α record. **Against data: nothing at
  A = 3 or A = 7.**
- **Missing.** Any A = 3 tagged datum; the A = 3 tensor/D-wave structure;
  the FLUKA-gated incoherent-background purity (`plans/05` §5.E — left empty
  rather than guessed); an A = 3 wave function with MC errors; the neutron
  polarization gate.
- **Development order.** (1) g₁ backend — the struck proton carries the
  signal, so g₁ᵖ's shape is the systematic: measure NNPDFpol11 against
  E143/HERMES/COMPASS g₁ᵖ (one vendoring away from the g₁ᵈ harness); (2) the
  triton SF default priced on the tagged A_∥ and the spectator spectrum;
  (3) the polarized-EMC transfer as a reachable kernel option (`medium_ratio`
  is empty in every default kernel).

### 3.4 Inclusive ⁶Li, A_zz / b₁(⁶Li) (rank 4)

- **Buys.** Paper 2's inclusive observable, A_zz = −(2/3) b₁/F₁ (the
  relation gated at 1e−12, U-9), and b₁(⁶Li) itself — the first A > 2 b₁ of
  any kind.
- **Rides on.** `Li6B1(MillerB1)` — Miller's digitized b₁ᵈ × 0.5 (row 2,
  likely not certain) × `LI6_B1_RANK2_TRANSFER` = 0.921949 (the α–d tensor
  dilution at P_D_LI6 = 0.0867, a fitted effective depolarisation, row 10) ×
  2/6 (an inert α, a modelling assumption); the ±100 % band from Q(⁶Li) vs
  Q_d; F₁ from `ToyF2` and `r_sigma_lt`; the RC band (row 17).
- **Validated today.** `t2_hermes_b1_table2` FAIL — Miller ×0.5 at χ²/6 =
  5.26 is the closest, the convolution camp indistinguishable from b₁ = 0,
  and the data do not decide row 2 (Δχ² = +0.035); the T1 rows FAIL for
  every backend on absolute F₂ᵈ (phase 0: toy 1.060, MSTW 1.025, CT18NLO-LO
  1.086 model/data) while the ratio passes (MSTW p = 0.044, CT18NLO 0.090).
- **Missing.** b₁(⁶Li) beyond the α–d picture (U-1: no calculation for any
  A > 2); any anchor on the transfer; a measured R.
- **Development order.** (1) the default-input price table for A_zz
  (phase A here) — done: ×1.2524 (CT18NLO) / ×1.2604 (MSTW), ×0.9555 under
  `--r-source r1998` on the `li6-convolution` path
  (`../open_items/run_2026-10-26/phase_A_baseline.md` §4); (2) row 2 by a theory input, since HERMES cannot decide
  it; (3) the opt-in `li6-convolution` with `--b1-unpol mstw` as the quotable
  configuration, as the lift already requires.

### 3.5 Inclusive ⁷Li (rank 5)

- **Buys.** Paper 3's inclusive observable ΔR_A(x) = g₁ᴬ/(P_p g₁ᵖ + P_n g₁ⁿ):
  CBT against TMT, a separation of 0.021 at x = 0.36 and 0.006 at 0.65 on the
  transferred curves, 0.45σ per bin at 10 fb⁻¹/u (`plans/01` §2.3).
- **Rides on.** `ToyG1` (whose g₁ⁿ ~ x^−1.22 below x ≈ 10⁻³ is not
  integrable); the effective polarizations (VMC); a rank-2 sector that is
  exactly zero and says so; a rank-3 (octupole) moment R₃ ≠ 0 in the
  default spin-temperature ladder fill that A_∥ does pick up and **no
  rank-3 structure function prices** (`PHYSICS_CHANNELS.md` §13); **no**
  polarized-medium modification reachable from generation.
- **Validated today.** `t1_g1d_world` FAIL on ToyG1, ToyG1-on-MSTW,
  ToyG1-on-CT18NLO and NNPDFpol11 (phase 0: COMPASS 45.0/15) — the last
  pointing at the harness's D-state / F₁ construction before the PDF
  (phase D); the Bjorken integral FAIL on ToyG1, within 10–19 % on NNPDFpol;
  the rank-2 zero pinned.
- **Missing.** The reachable polarized-EMC ratio; a ⁷Li rank-2 structure
  function (U-4); the rank-3 contamination priced (ladder vs `--pzz-mode
  typed`, where R₃ = 0, is the generator-side price); a g₁ᵖ and a g₁ⁿ(³He)
  gate; any ⁷Li DIS datum (none exists, `03` §8).
- **Development order.** (1) the D-state question and the NNPDFpol verdict
  (phase D); (2) the `medium_ratio` hook made reachable; (3) the rank-3
  price; (4) the rank-2 sector stays zero until row 20 is decided.

### 3.6 Tagged ⁶Li α (rank 6)

- **Buys.** Paper 2: the tagged tensor asymmetry on the embedded deuteron,
  O(1) at k ≈ 0.3 GeV from D-wave dominance; the α–d D-wave itself.
- **Rides on.** The α–d S + D wave (Hulthén by default, VMC AV18 opt-in);
  the embedded deuteron (same flag since 2026-09-04); the T1 breakup; the
  spin-blind Glauber FSI weight; the i^L phase fixed on 2026-09-06 (row 26).
- **Validated today.** The Cosyn–Weiss gate as an identity on the deuteron
  control (8.88e−16 over 26 880 cells); the port gate. **Against data:
  `t2_deuteron_static` FAIL for the shipped Hulthén (η +17.2σ; D2 keeps it
  and quotes physics on AV18); `t2_deeps_spectator_tail` FAIL for every
  wave and both prescriptions; `t2_bonus_spectator_shape` FAIL recorded
  2026-10-26: the 40 mb Glauber weight's backward-window variation exceeds
  the data's 1σ flat band in 24 of 28 DIS blocks, while the data fall
  toward 90° as the weight does (χ² 739.4 flat against 226.3 for the
  weight's shape, 84 dof) — the rule's premise fails, not the weight's
  direction; `t2_boeglin_nk` (row 16, 2026-10-26) PASS for all three waves
  only inside the paper's 20 % FSI/off-shell allowance, Hulthén the worst
  shape (21.98/21 against AV18 11.54, CD-Bonn 8.89).**
- **Missing.** Anything with an α spectator (U-2); the (e,e′d) / (e,e′α)
  α–d momentum distributions (Ent, Mitchell — paywalled, digitizable); the
  D wave from data (U-8); an FSI with spin dependence.
- **Development order.** (1) the wave-function default — A_zz^tag(k) is
  −1.2069 / −0.9533 / −0.5851 on Hulthén against −0.5191 / −0.1993 / +0.1170
  on AV18 at k = 0.20 / 0.30 / 0.50 GeV (most of it P_D: Hulthén at the VMC
  P_D gives −0.6001 at 0.20), so D2's revisit at M2 is this channel's first
  item; (2) read the Deeps failure for what it tests — n(k), the
  harness-coded light-cone map and FSI together; (3) BONuS as the low-p_s
  complement, with the resonance-region blocks declared out of model;
  (4) IP6 optics: the physics case is the IR-8 secondary focus, outside
  this generator.

### 3.7 Tagged d–p control (rank 7 by rule R1)

- **Buys.** Validation: the Cosyn–Weiss tagged limit, the row the spectator
  anchors run on (CW, Deeps, BONuS, d(e,e′p)n). Its own physics is what the
  EIC may not be able to run, which is the premise of the lithium programme.
- **Rides on, validated, missing.** As 3.6 at A = 2, plus: the toy g₁ⁿ's
  wrong sign over 0.25 < x < 0.6 bites hardest here (a free neutron is the
  DIS target); `--cluster-wave vmc` moves its vector/tensor dilutions
  −2.03 % / −1.18 % (row 9). The RC elastic tail does **not** reach this
  channel (`rc_tail` ≡ 1 on every tagged channel); the deuteron elastic form
  factor exists only as the T10 normalisation check and an inclusive
  deuteron RC run throws today — so Abbott's t₂₀ / A(Q²) data are a
  development item (a deuteron `TabulatedSpin1FF` backend), not a benchmark
  of any generation path (§5.4).
- **Development order.** (1) d(e,e′p)n momentum distributions (Boeglin),
  with the paper's own 10–20 % FSI allowance declared; (2) BONuS and Deeps
  as above. *(2026-10-26: (1) wired as row 16, PASS under the 20 % allowance,
  all three waves failing without it; BONuS wired as row 10, FAIL recorded.)*

## 4. Cross-cutting ingredients, ranked by leverage on ranks 1–3

| rank | ingredient | moves | status | the development item |
|---|---|---|---|---|
| **A** | **unpolarized nucleon inputs F₂, R** (`--unpol-sf`, `--r-source`) | on the cos 2φ amplitude (rank 1) **F₂ cancels exactly** and R is the lever: −15.6 / −15.1 / −7.1 / +4.6 % at the four sweet-spot cells x = 0.011 / 0.028 / 0.071 / 0.141 generator-side (`r1998` vs `r_sigma_lt`, `t1_r_world`'s kernel factor; −15.2 / −14.7 / −7.8 / +4.2 % over the generated super-bins on the `li6-convolution` path, `price_defaults.py`), the sibling's whole chain (0.842 / 0.849 / 0.930 / 1.039) to ≤ 0.9 %; the −8.3 / −6.7 / −3.1 % quoted until 2026-10-26 are `docs/USAGE.md`'s function-level points x = 0.05 / 0.10 / 0.30 at y = 0.5, Q² = 2.5, not the sweet spots (`../open_items/run_2026-10-26/phase_A_baseline.md` §4, `phase_C1_r_world.md` §4); F₂ enters only the rate (N_DIS ×0.766, δA ×0.95–1.16). On A_zz (rank 4): toy → CT18NLO ×1.2524 window-level (0.86–1.44 per x bin); `r1998` vs `r_sigma_lt` is +26.4 % on K/D_φ at x = 0.10, y = 0.5, Q² = 2.5 **on the opt-in `li6-convolution` path, the only one that accepts `--r-source`** (y-dependent) | phase 0: **no backend passes absolute F₂ᵈ** (toy 8307.9/155, MSTW 416.5/155, CT18NLO at LO 3401.8/146 — the 8.6 % excess is of the size of the missing NLO coefficient functions); the ratio passes for MSTW and CT18NLO; *2026-10-26:* NMC's 2.5 % normalisation is outside 'sys' (phase D) and cannot absorb CT18NLO's 8.6 % or the toy's 6.0 % excess; HERMES F₂ᵈ (row 14) fails for every backend too (ToyF2 5303.5/68, MSTW 318.9/68, CT18NLO-LO 1788.2/45); R: `t1_r_world` (row 12) PASS on NMC R (12.54/10) but `r_sigma_lt` fails the SLAC sets (209.76/68, 29.66/9) where `r1998` passes them (its own fit input); NMC fixes (1 + R) at the sweet-spot x to ±7.5–9.3 % | R data wired 2026-10-26 (row 12, phase C1); next, a registry decision at M2 on the default backend, framed as "toy vs LO-evaluated grids: both fail absolute F₂ᵈ, the ratio passes, claim 1 sees only R and the rate" |
| **C** | **the Δ sector** (`toy_delta_gluon` at 10⁻², the finite-γ kernel) | claim 1's extracted xΔ directly | **one moment-level bound under stated assumptions** (`t5_delta_moment_bound`, row 13, PASS 2026-10-26): the shipped toy uses 6–9 % of the gluonic Soffer room (13–19 % under DS16's own normalisation of Δ, V6); the letter's `moment_A` / `moment_B` at dilution 1/3 use 0.67–0.74 under the first-pass normalisation and 1.35–1.48 (beyond the bound) under V6; the finite-γ leakage 0.14–0.16 γ²b₁ is a carried systematic | the factor-2 cross-section normalisation of Δ in the kernel's cos 2φ term (DS16 Eqs. (3)–(5) and LOI Eq. (6) against LOI Eq. (9) and the kernel; a registry candidate, `../open_items/run_2026-10-26/phase_C2_delta_bound.md` §9); f₁LL^g for ⁶Li (nothing calculates it). A bound says the scenario is *allowed*, never that it is right |
| **E** | **coherent scenario** (`eps_b0`, `f0`, `amp`, `slope_b`) | the deformation a₂ by 11.4× and the generated c₂ by ≈ 6× at the tag ⟨\|t\|⟩; N_tag by `f0`'s factor-2 band | `slope_b` anchored; the rest not | row 8 at M2; row 21 |
| **G** | **the deuteron wave function default** (Hulthén) | A_zz^tag ×2.3–4.8 at k = 0.2–0.3 GeV with a sign change by 0.5 GeV (rank 6); −2.03 % / −1.18 % on the d–p control's dilutions | `t2_deuteron_static` FAIL +17.2σ against Rodning–Knutson's measured η_d (read, phase D); D2 keeps the default until M2; *2026-10-26:* Boeglin's d(e,e′p)n (row 16) passes all three waves only inside the paper's 20 % FSI allowance, Hulthén the worst shape (21.98/21; AV18 11.54, CD-Bonn 8.89; only CD-Bonn passes with the MS FSI factor) | D2 at M2 — the largest measured lever on paper 2's tagged observable |
| **D** | **b₁(⁶Li) model** (`--b1-model`, row 2, the transfer) | A_zz(⁶Li) by ×2 (row 2) and by the ±100 % band | HERMES does not decide row 2; the α–d step is validated by nothing (U-1) | a theory input is the only way forward; keep the band mandatory |
| **F** | **radiative corrections** (`--rc`, `--rc-tail-model`) | cos 2φ: the leading-log ISR bound ≤ 2.9 % (`plans/07` WP4), no full O(α) calculation (U-6), `RcScope::TensorAll` prices the omission; A_zz: the band's half-width (row 17, the low-x anchor errs wide by ×2.65); nothing on the coherent channel by design | LiPolGen's elastic tail is within 0.25×–4× of POLRAD's (0.88 / 0.95 / 0.39; HO form factor); POLRAD's own rows re-drive 2 of 3; the polarized QE tail is unpriced; row 7 (2026-10-26): against DJANGOH 4.6.21's exact O(α) proton elastic tail at 18 × 275, Q² 1–10⁴ GeV², `polrad-full` R = 1.020 / 1.026 / 1.022 / 1.010 (PASS at 3 % + 2δ) and 1.000 / 1.001 / 1.000 / 0.987 straight from the quadrature; the shipped `t-peak` R = 0.61 / 0.90 / 1.01 / 1.01 (0.60 / 0.88 / 0.99 / 0.98 direct); the (x, y) tail table adds +1.9 … +2.5 %; *2026-10-26:* F_m's shape against RFY's ⁶Li F_M² (row 6): PASS at p = 0.038 with a free scale over q = 0.85–1.39 fm⁻¹ (0.0145 at s = 1); the C0 zero against Wiringa–Schiavilla (row 15): PASS one-sided (q₀ = 3.0998 ≥ 2.8143 fm⁻¹), while row 5's data band sits below it | row 7 run and PASS 2026-10-26 (phase B3); the tail table's interpolation (+2 %) and the `TPeak` default (−39 % / −10 % of the exact proton tail at Q² 1–10 / 10–100 GeV²) are candidates for `DECISIONS_PRICED.md` (phase E) |
| **H** | **polarized nucleon input g₁** (`--pol-sf`) | A_∥ on ranks 3 and 5; nothing on claim 1 | `t1_g1d_world` FAIL on every backend including NNPDFpol11 (phase 0); phase D read the four papers: all publish g₁ᵈ per nucleon, NOT D-state corrected, so the verdict column stands and the D-state is cleared (the other reading is worse for NNPDFpol11 on all three sets); `LhapdfG1`'s leading-order evaluation of an NLO set is the open reading (the Bjorken integral sits at 98 % of g_A/6, without C_Bj) | phase D; a g₁ᵖ and a g₁ⁿ(³He) gate |
| **B** | **beam and fill assumptions** | the two-fill pair (+0.6, −1.2) is built by default; whether a Li source delivers P_zz = −1.2 at matched purity is external (no source table, N-2; EPIOS modes reach −1); on ⁷Li the default ladder fill's R₃ ≠ 0 contaminates A_∥ with no structure function to price it | recorded, not asserted | state the assumed fill in the letter; price R₃ by ladder vs `--pzz-mode typed` |
| **I** | **spectator FSI** (`--fsi`) | the tagged spectrum's shape, never A_zz (spin-blind) | unitarity not definable; the two variants differ on 99.5 % of events; BONuS (row 10, 2026-10-26): the weight's direction at p_s = 70–150 MeV/c is the data's, its 40 mb backward variation exceeds the flat band (4/12, 0/16 blocks); the forward excess is what its cos → −cos symmetry cannot carry | the Deeps transverse window and BONuS (DIS blocks only) as direction checks |
| **J** | **T2 hadronization** | reco-level FOMs on the sibling side | `t4_pythia_ep_closure` PASS (worst 0.076 on the pinned 8.317; 0.023 on 8.312) | none for the letter |
| **K** | **output and chain** | nothing physical | `t6_hepmc3_convention` 23/23; the ePIC legs unrun since 2026-09-02 | M3 |
| **L** | **⁷Li rank-2 sector** | nothing: identically zero | pinned zero; deferred on row 20 | waits on theory input |

## 5. What the ranking does to the benchmark plan

1. **A proposed clause for the M2 rule** (a registry candidate, not an
   edit of `BENCHMARK_PLAN.md` §3): "no benchmark moved beyond its stated
   tolerance since M1" is necessary and not sufficient — at M2 the report
   must have been regenerated on a machine with the three LHAPDF sets
   (done 2026-10-26), and each of ranks 1–4 must carry a **default-input
   price table of its own observable** (the cos 2φ amplitude at the sweet
   spots, a₂, A_∥, A_zz) under every knob the knob-provenance table reports
   as read, so that the registry can decide the defaults from numbers.
2. **Rows to wire, in rank order** (this run, `../open_items/run_2026-10-26/`,
   future tense until its integration stage records them): R data with
   NMC R primary (ingredient A); the Δ moment bound (C); BONuS (row 10) and
   RFY (row 6), blocked on a document until today; the DJANGOH elastic-only
   run (row 7, run 2026-10-26: PASS); HERMES F₂ᵈ as a T1 extension; the
   Wiringa–Schiavilla ⁶Li form factors as a theory anchor for the C0 zero;
   Boeglin's d(e,e′p)n with its FSI allowance. *Recorded 2026-10-26 by the
   integration stage: all wired — R row 12 PASS, the Δ bound row 13 PASS,
   BONuS row 10 FAIL recorded, RFY row 6 PASS, DJANGOH row 7 PASS, HERMES
   F₂ᵈ row 14 FAIL recorded, Wiringa–Schiavilla row 15 PASS (one-sided),
   Boeglin row 16 PASS (under its allowance).*
3. **Rows that stay FAIL by decision are read against this ranking.** The
   Hulthén η (D2) and the C0 zero (≤ 2.2e−7 on ΔA_zz) can stay recorded
   until M2. The toy F₂ / g₁ failures (A/H) and `eps_b0` (E) sit under ranks
   1, 2 and 4 and are the decisions this file asks the author to take at M2.
4. **What is deliberately NOT raised.** The ePIC chain legs (K): they
   validate the record, not the physics, and need Phase 2's container
   stack; ⁷Li rank-2 (L): no input exists; FSI unitarity: not definable;
   Abbott's deuteron t₂₀ / A(Q²): no generation path reads a deuteron elastic
   form factor (an inclusive deuteron RC run throws), so it is a development
   item — a deuteron `TabulatedSpin1FF` from the Abbott EPJ A 7
   parametrization, which would also unblock that run.

## 6. Triggers that re-rank

- A change of the letter's scope (`plans/07` D1): if the coherent channel is
  dropped from the letter, rank 2 falls below rank 4 and ingredient E to the
  bottom of §4.
- An IP6 optics decision for ⁶Li tags (`plans/10`): the tagging optics reach
  ε = 0.25–0.37 only at 1/7–1/13 of the luminosity, so tagged events per
  year gain ×1–3, not ×10; the IR-8 secondary focus (≈ 20 % acceptance)
  gives ε × D² ≈ 0.02, still C2 = 1. Rank 6 moves only with a rescore.
- A theory input: a dipole-model a₂ for polarized ⁶Li (row 21) moves
  ingredient E from "scenario" to "calculation"; a b₁(⁶Li) calculation
  beyond the α–d picture does the same for D.
- An EIC lithium luminosity or polarization number (gap #4 of `plans/01`):
  every C2 entry is at 10 fb⁻¹/u and P_zz = 0.6 by assumption.
- Data: JLab E12-13-011 (b₁ᵈ, 0.16 < x < 0.49) changes the HERMES-only state
  of row 3; E12-14-001 (polarized EMC, ⁷LiD) gives rank 5 its first datum.
- Registry decisions on rows 2, 3, 8, 9 (D2), 17 and 20: each one removes a
  C3 point from the channel it sits under.

## 7. Revision log

| date | what changed | by |
|---|---|---|
| 2026-10-26 | first version, from the 2026-09-27 report and the 2026-09-27 decisions | run `../open_items/run_2026-10-26/` |
| 2026-10-26 | revised after two adversarial readers (49 findings, physics and honesty/feasibility lenses; `../open_items/run_2026-10-26/STATUS.md` phase 0r): the rubric's levels defined and applied literally with three stated override rules; rank 1 split into cos 2φ (claim 1) and A_zz/b₁ (paper 2); F₂ removed as a lever on the cos 2φ amplitude (it cancels; R is the lever); the fill model demoted (the two-fill pair is built by default); the deuteron wave function priced on A_zz^tag (×2.3–4.8); the Abbott t₂₀ item re-filed as development; the PASS list, the Deeps p-values, the triton polarizations, the Δ scale (10⁻²), the δA footings and the O5 qualifiers corrected; phase 0's CT18NLO / NNPDFpol11 numbers folded in | the supervising session |
| 2026-10-26 | the run's measurements folded in (phases A, B1, B2, B3, C1, C2, C4, C5, D, P; `../open_items/run_2026-10-26/STATUS.md`): the R price at the sweet spots corrected from the function-level −8.3 / −6.7 / −3.1 % (USAGE.md's x = 0.05 / 0.10 / 0.30, y = 0.5, Q² = 2.5) to −15.6 / −15.1 / −7.1 / +4.6 % at the cells (−15.2 / −14.7 / −7.8 / +4.2 % over the generated super-bins) in §2 rank 1 and §4 A; the shipped Δ's amplitude recorded as the opposite sign of `moment_A` and 4–7× smaller (§3.1); row C wired (`t5_delta_moment_bound` PASS for the shipped toy at 6–9 % of the gluonic Soffer room; DS16's own definition of Δ makes the kernel's cos 2φ coefficient twice DS16's per unit Δ, so the letter's scenario uses 0.67–0.74 or 1.35–1.48 of the bound — a registry candidate); R (row 12) and HERMES F₂ᵈ (row 14) in §4 A; BONuS (row 10) FAIL, Boeglin (row 16) PASS under its allowance, in §3.6, §3.7, §4 G and I; RFY (row 6) and Wiringa–Schiavilla (row 15) in §4 F; the g₁ᵈ D-state convention cleared and η_d read in §4 G and H; §5 item 2 recorded. **No rank moved**: no worker's measurement changed a criterion score — rank 1 keeps C3 = 3 (Δ is still a scenario; a bound is not an external calculation) and C4 = 1 (R's new anchor is at A = 1, 2, another A; the Δ bound is a bound, not an anchor); ranks 3 and 5 keep C3 = 3 (`--pol-sf nnpdfpol` moves A_∥ ×0.51–0.53); rank 6 keeps C3 = 3 (A_zz^tag ×0.43 / ×0.21 under `--cluster-wave vmc`); the row-7 RC check (ingredient F) scores no channel | the integration stage (phase E) |
