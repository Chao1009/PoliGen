<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-26 — decisions priced, none taken

The decision registry is still `../run_2026-09-03/STATUS.md` (27 rows); a
row is cited here by its number there. This file lists what run 2026-10-26
**priced** — every registry row a worker measured a cost for, the asks of
`../../benchmarking/08_channel_priority.md` §5, and the candidates this run
files — as options with their measured cost and the window the cost was
measured in. **No option was taken, no default, constant, tolerance, registry
label or reference number moved.** Each number is copied from the phase
record named with it; where two phases measured the same quantity in
different windows, both are given. Written by the integration stage (phase E).

Shorthand for the price-table window (phase A3, `validation/price_defaults.py`,
`phase_A_baseline.md` §4.1): config 1 (e 10 GeV × ion 99.5 GeV/u), `lg.run`
fixed-count mode, **N = 200 000, seed 20261026**, P_z = 0.7, P_zz = 0.6,
P_e = 0.7 (helicity-flip only), YR high-acceptance optics, x_max 1.0 (0.95 on
the `li6-convolution` rows), one knob at a time. The four sweet-spot cells /
super-bins of claim 1 are listed in **x order, x = 0.011 / 0.028 / 0.071 /
0.141** (Q² = 1.14 / 1.14 / 3.13 / 14.3 GeV²) everywhere below; the phase A
table calls them S2 / S1 / S3 / S4.

## 1. Registry rows priced by this run

### Row 2 — Miller's ×0.5 stays (`B1_MILLER_TABLE_TO_PER_NUCLEON` = 0.5)

- **What was measured** (`t1_hermes_f2d`, phase C1, `phase_C1_r_world.md` §6):
  the premise of the row's arbiter, that HERMES's published b₁ᵈ is per nucleon,
  checked against HERMES's **own** measured F₂ᵈ (JHEP 05 (2011) 126, HEPData
  ins894309) — the F₂ᵈ implied by HERMES PRL 95 Table II (with R1990, as
  HERMES used) over HERMES's F₂ᵈ in the nearest bin = 0.929 / 1.039 / 0.986 /
  0.958 / 1.006 / 1.006 at x = 0.012 / 0.032 / 0.063 / 0.128 / 0.248 / 0.452
  (geometric mean 0.986) → **per nucleon** by the declared rule (all six in
  [0.8, 1.25]). Row 3 of the plan still does not decide it: χ²(b₁ᵈ)/6 = 5.262
  (×0.5, shipped) vs 5.297 (×1).
- **Options:** keep ×0.5 (shipped) / ×1. **Cost:** the premise of ×0.5 now
  stands on HERMES's own F₂ᵈ, not only on ALLM97; Miller's own normalisation
  is untouched by it. **None taken.**

### Row 3 — `Li6ConvolutionOptions` keeps `r_sigma_lt`; more broadly, the default R (`r_sigma_lt` vs `r1998`)

- **Against data** (`t1_r_world`, phase C1; Q² ≥ 1, W² ≥ 4 GeV², stat ⊕ sys
  diagonal unless stated): NMC R **12.54/10 (p = 0.25) vs 6.20/10 (p = 0.80)**;
  NMC's sys taken 100 % correlated 60.22/10 (p = 3.3e−9) vs 21.28/10
  (p = 0.019); Whitlow's SLAC R 209.76/68 (p = 2.3e−16) vs 63.05/68 (p = 0.65);
  E140 29.66/9 (p = 5.0e−4) vs 12.71/9 (p = 0.18); H1 F_L 8.82/8 vs 8.18/8.
  Whitlow and E140 are R1998 fit inputs, so for `r1998` they are a closure.
- **On claim 1's observable:** the cos 2φ amplitude at the four sweet spots
  moves **×0.844 / 0.849 / 0.929 / 1.046** at the cell centres (the kernel's
  factor, `t1_r_world`; (1 + r1998)/(1 + r_sigma_lt) = 1.1848 / 1.1780 /
  1.0762 / 0.9563) and **×0.8485 / 0.8534 / 0.9221 / 1.0420** over the
  generated super-bins (phase A, on the `li6-convolution` path — `--r-source`
  is refused on the shipped `miller` b₁ and accepted only there); the bias
  spans 0.921–1.192 over the window (C1).
- **On paper 2's observable** (phase A, `li6-convolution`, x_max 0.95):
  run-level A_zz ×0.95554 (bands 1 and 2; zero at band 0 either way),
  ×0.962683 with `--b1-unpol mstw`; σ ×0.993153. The README's function-level
  points reproduce: K/D_φ −3.8357 / +26.3988 / +3.3518 % and Δ/D_φ −8.3268 /
  −6.7268 / −3.1441 % at x = 0.05 / 0.10 / 0.30, y = 0.5, Q² = 2.5 GeV².
- **Options:** `unset` / `sigma-lt` (bit-identical) vs `r1998`. **None taken.**
  Note (phase A): the R choice on claim 1 cannot be exercised on the shipped
  configuration at all — registry row 3 is framed on `Li6ConvolutionOptions`;
  the claim-1 R choice is a separate question.

### Row 5 — the `DeuteronConvolutionB1` wave default stays AV18 (`kFdeutFile`)

- η_d against Rodning–Knutson's 0.0256 ± 0.0004 (now read as a measurement,
  phase D): AV18 −1.39σ, CD-Bonn −0.07σ — both inside 2σ; unchanged.
  **None taken.**

### Row 8 — `CoherentScenario::eps_b0` stays −0.08 (against the measured-quadrupole −0.0070)

- Phase A, coherent ⁶Li, `transverse-tensor`, P_zz 0.6, B = 50, |t| ≤ 0.2:
  c₂ = 0.029989 → **0.008099 (×0.27007)** at the generated spectrum's exact
  mean |t| = 0.019991 GeV², and 0.0780 → **0.0123 (×0.15769)** at the tag
  ⟨|t|⟩ = 0.06 GeV²; the generated c₂ 0.0336 ± 0.0032 → 0.0117 ± 0.0032;
  σ_coh = 12476.03 pb unchanged; 0 of 200 000 recoils RP-accepted at the YR
  envelope (either setting). **None taken.**

### Row 9 (the Hulthén deuteron control at P_D = 0.045), with decision D2 of `../run_2026-09-27/DECISIONS.md` (Hulthén stays, physics on AV18, revisit at M2) and row 10

- **On paper 2's tagged observable** (phase A, `--cluster-wave vmc`): tagged
  ⁶Li α A_zz^tag −1.2069 → −0.5191 (**×0.430**) at k = 0.1979 GeV and
  −0.9533 → −0.1993 (**×0.209**) at 0.3012 GeV (YR high acceptance,
  `n_phi` 32); integrated dilutions ×1.11613 / ×1.06576 (row 10's +11.61 % /
  +6.58 %, reproduced); the ⁶Li α RP tag fraction 0.02466 → 0.03422
  (sampled). Tagged ⁷Li α: ⟨k⟩ 0.28852 → 0.18608 GeV (×0.645), A_∥
  bit-identical, tag fraction 0.97386 → 0.99839.
- **Against data:** η_d (phase D) Hulthén 0.032478 **+17.20σ** (FAIL at
  \|pull\| ≤ 2), AV18 −1.39σ, CD-Bonn −0.07σ; the impulse-approximation Q_d
  against the corrected 0.2860(15) fm²: +1.52 / −5.71 / −5.42 % (recorded).
  Boeglin's d(e,e′p)n at 35–45° (phase C5, 22 points, p_m 0.055–0.538 GeV/c,
  one free scale): under the paper's 20 % allowance Hulthén 21.98/21
  (p = 0.40), AV18 11.54 (0.95), CD-Bonn 8.89 (0.99); with the MS FSI factor
  instead 95.6 / 76.7 / 24.4 (only CD-Bonn passes, p = 0.27); with no
  allowance 185.9 / 66.5 / 60.3 (all fail). BONuS (phase B1, σ_XN = 40 mb):
  the weight's backward variation 0.0344 / 0.0548 / 0.0848 / 0.1406
  (Hulthén) vs 0.0355 / 0.0555 / 0.0871 / 0.1505 (AV18) at p_s = 0.078 /
  0.093 / 0.11 / 0.135 GeV/c; blocks in the band 4/12 + 0/16 vs 3/12 + 0/16;
  χ² against the weight's shape 226.3 vs 221.6 on 84 dof — the wave barely
  moves that row. The Deeps row (unchanged) reads 627.85 / 579.12 / 278.81
  on 187 dof.
- **Options:** Hulthén (shipped) / AV18 (`--cluster-wave vmc`) / CD-Bonn.
  The switch's own release-wide cost (every reference pin, the port gate, the
  knob tables) was not re-measured here. **None taken.**

### Row 23 — `LI6_CLUSTER_POLARIZATION` = 0.811228 against the naive 1/3

- On `t5_delta_moment_bound` (phase C2) it enters only through
  g₁^g = P_eff × Δg: P_eff = 1/3 instead of 0.2704 widens the bound B by
  +0.63 / +0.58 / +0.55 % at Q² = 2 / 5 / 10 GeV² and lowers every
  saturation by the same factor (toy 0.06404 / 0.07980 / 0.09355). No
  verdict moves. **None taken.**

## 2. The asks of `08_channel_priority.md` §5 and PLAN.md E3

### The default unpolarized backend (`--unpol-sf` toy / mstw / ct18nlo)

- **Against data** (Q² ≥ 1, W² ≥ 4 GeV²; CT18NLO on Q² ≥ 1.677): NMC F₂ᵈ
  8307.9/155 / 416.5/155 / 3401.8/146 (all FAIL; mean model/data 1.060 /
  1.025 / 1.086); NMC F₂ᵈ/F₂ᵖ 2161.4/211 (FAIL) / 247.2/211 (p = 0.044) /
  215.5/189 (p = 0.090); HERMES F₂ᵈ (per-bin diagonal + correlated 7.5 %
  normalisation) 5303.5/68 / 318.9/68 / 1788.2/45 (all FAIL). No backend
  passes absolute F₂ᵈ; the ratio passes for both grids. NMC's 2.5 %
  normalisation is outside 'sys' (phase D) and cannot absorb CT18NLO's 8.6 %
  or the toy's 6.0 % excess (MSTW's 2.5 % is 1.0 normalisation σ: open).
- **On the observables** (phase A, ratio to toy, mstw / ct18nlo): inclusive
  ⁶Li σ ×0.79338 / ×0.79847, run-level A_zz ×1.26044 / ×1.25239; the cos 2φ
  amplitude identical per cell (8.9e−16 over 3051 cells) and ≤ 0.23 % at the
  super-bins; inclusive ⁷Li A_∥ ×1.0262 / ×1.0602 (σ ×0.7913 / ×0.7959);
  tagged ⁷Li α A_∥ ×1.0248 / ×1.0552; σ_coh ×0.7439 / ×0.7058; A_zz^tag and
  c₂ unchanged. **None taken.**

### The default polarized backend (`--pol-sf` toy / nnpdfpol)

- **Against data:** world g₁ᵈ FAILs for every backend — ToyG1 107.8/15,
  69.4/15, 54.9/28; NNPDFpol11_100 45.0/15, 52.6/15, 51.0/28 (COMPASS,
  HERMES, E143); the papers publish g₁ᵈ per nucleon, uncorrected, so the
  verdict column stands (phase D). Bjorken: NNPDFpol11 at 98.0–98.2 % of
  g_A/6 at Q² = 2 / 5 / 10 — about g_A/6 without C_Bj, the leading-order g₁
  formula's signature (phase A).
- **On the observables** (phase A): inclusive ⁷Li A_∥ ×0.5316 (×0.6717 /
  ×0.6678 together with MSTW / CT18NLO F₁); tagged ⁷Li α A_∥ ×0.5142; not
  read on any tensor plan. **None taken.**

### `r_sigma_lt` vs `r1998` — registry row 3, §1 above. D2 at M2 — row 9, §1 above. `eps_b0` — row 8, §1 above. Row 2 on HERMES — §1 above.

### Not registry rows, priced on the way

- **The coherent `f0` band** (phase A): f0 = 0.02 / 0.08 gives σ_coh ×0.5 / ×2
  exactly; c₂ unchanged (the same events). None taken.
- **`--triton-sf ciofi-simula`** (S₀ = 0.6525, untuned; phase A, tagged ⁷Li α):
  A_∥ bit-identical, ⟨k⟩ ×0.99985 (sampling); only the T1 final state moves —
  `n_partner` on 45 969 of 200 000 events (mean 1.338 → 1.568),
  `struck_virtuality` on all (−0.12356 → −0.12518 GeV²), `struck_pol` on
  74 843 (0.19006 → 0.18779), while the per-category cross sections stay
  bit-identical (an open model question for the author: should the tagged
  rate read the struck polarization the breakup assigns?). None taken.
- **The σ_XN band 20–40 mb of `GlauberFsiOptions`** (phase B1, BONuS,
  W* ≈ 2–2.7 GeV, p_s 70–150 MeV/c): χ² against the weight's shape 226.3/84
  (40 mb) vs 403.3/84 (20 mb) vs flat 739.4/84; blocks in the band 4/12 + 0/16
  vs 8/12 + 13/16. Far below EIC W. None taken.
- **The `rc_tail` default model, `TPeak` (shipped) vs `PolradFull`** (phase B3,
  proton at e 18 × p 275, row 7's window; exact = DJANGOH 4.6.21 IEL31–33):
  straight from the quadrature TPeak/exact = 0.599 / 0.878 / 0.989 / 0.985 and
  PolradFull/exact = 1.000 / 1.001 / 1.000 / 0.987 at Q² = 1–10 / 10–100 /
  100–10³ / 10³–10⁴ GeV²; through the shipped table 0.612 / 0.897 / 1.011 /
  1.007 and 1.020 / 1.026 / 1.022 / 1.010. Switching costs the tail-table build
  ~2 s → 10–21 s and moves every published number. Row 17 untouched.
  None taken.
- **The RC tail table's interpolation** (phase B3): the log-bilinear (ln x,
  ln y) interpolation of the tail column on the 101 × 77 grid raises the
  window-integrated proton tail by +2.0 … +3.7 % (+1.9 … +2.5 % net of the
  `tail_max` clip), both models; options keep / refine the y rows /
  interpolate a smoother variable; the cost of a fix not measured; likely
  present on ⁶Li's quasi-elastic column, not measured. None taken.
- **`TaggedModel`'s nearest-cell-below lookup** (`n_of_kc`; phases C5, B1):
  +16 % (Hulthén) / +18 % (AV18) on n(k) at k = 0.055 GeV/c against the
  smoothly read table; Boeglin's verdict χ² 21.98 → 23.42 with the smooth
  reading; at the Deeps momenta k = 0.30 / 0.373 / 0.46 / 0.60 / 0.95 GeV the
  lookup reads n(k) +5.7 / +4.8 / +5.6 / +2.4 / +3.0 % high (Hulthén), its
  effect on the Deeps row's 627.85/187 unmeasured. Listed, not priced as a
  decision. None taken.

## 3. Candidates this run files (priced options only)

1. **The proposed M2 clause** (`08_channel_priority.md` §5.1, a proposal, not
   an edit of `BENCHMARK_PLAN.md` §3's adopted rule): "at M2 the report must
   have been regenerated on a machine with the three LHAPDF sets, and each of
   ranks 1–4 must carry a default-input price table of its own observable
   under every knob the knob-provenance table reports as read". Cost, as
   measured this run: the regeneration with the sets is done
   (`phase_E_report_compare.md`: 23 harnesses, 13 / 10 / 0, 77 s with
   `--plots`); the price table exists for ranks 1–4 and more
   (`validation/price_defaults.py`, 183–190 s on 4 shared cores, two runs
   identical over all 289 common values). Options: adopt the clause into the
   M2 rule / keep the adopted rule. **None taken.**
2. **A deuteron `TabulatedSpin1FF` from Abbott *et al.*, EPJ A 7 (2000) 421**
   (C3, deferred with design): no generation path reads a deuteron elastic
   form factor today (`rc_tail` ≡ 1 on the tagged channels; the inclusive
   deuteron `HoSpin1FF::for_ion(DEUTERON())` throws), so a t₂₀ / A(Q²)
   comparison would move nothing and `HoSpin1FF`'s t₂₀ fails by construction
   (one monopole for F_c and F_q). Options: build the backend (it would also
   unblock the inclusive-deuteron RC run) / leave it. Cost: not measured.
   **None taken.**
3. **Shipping the Wiringa–Schiavilla table into `data/ff/li6_elastic.csv`**
   (phase C4; `tests/test_rc.cpp`'s "NO 6Li table is shipped" case would flip).
   (a) keep `HoSpin1FF` (status quo); (b) promote the 68 slot rows of
   `ws1998_li6_formfactors.csv` by one copy (tested: `TabulatedSpin1FF` reads
   exactly those rows from a scratch `LIPOLGEN_DATA_DIR`). Cost of (b): support
   only q = 0.1190 … 2.9190 fm⁻¹ and the reader refuses to extrapolate, so
   every elastic-tail evaluation beyond 2.919 fm⁻¹ throws (the t-peak reaches
   q ≈ 88–278 fm⁻¹) — a continuation law is needed; the F_C2 column carries
   the unresolved 3/√2 (the Q closure 0.4901 vs 0.23 fm²); F_M1 is interpolated
   between 20 Monte Carlo symbols through its zero. Gain: a C0 shape inside the
   calculation (`ho` is ×11 too steep at 2.5 fm⁻¹ after a free scale) and a C2
   second lobe (1.1e−6 vs 4e−14). Beside it, the constraint row 15 adds on any
   HO refit (Q1): the C0 zero ≥ 2.8143 fm⁻¹ (≥ 3.0952 against the IA line) —
   row 5's FB zero (2.6944) would fail it, Li71's minimum (2.8284) passes by
   0.0141, the shipped 3.0998 passes. **None taken.**
4. **The cross-section normalisation of Δ in the kernel's cos 2φ term**
   (phase C2 §4.2, §9): `src/core/xsec.cpp`'s a₂ = LOI18 Eq. (9) carries
   −(1−y)/y² Δ, twice what DS16 Eqs. (3)–(5) and LOI18 Eq. (6) give per unit
   Δ. (i) Read the kernel's Δ as Δ_DS16/2 and carry C_BAG/2 for the bag
   scenario: every letter cos 2φ amplitude at the bag normalisation halves;
   the scenario uses 0.67–0.74 of the gluonic Soffer room. (ii) Halve the
   kernel's coefficient: the same physics as (i); the toy's 10⁻² gives half
   its amplitude. (iii) Keep both: the letter's amplitudes encode twice the
   bag value, 1.35–1.48 of the room under the central assumptions, allowed
   only for f₁LL^g ≳ +0.7 … +1.0 f₁^g. A copy of Jaffe–Manohar PLB 223
   (`REFERENCES.md` §1 #4) would settle which printed formula is theirs.
   **None taken.** (Not a registry row, priced alongside: carrying the ⁶Li
   rank-2 transfer 0.921947 on Δ as on b₁ would scale the letter's `moment_A`
   saturations to 0.659 / 0.644 / 0.647 (A7) and 1.318 / 1.288 / 1.294 (V6) at
   Q² = 2 / 5 / 10 GeV².)
5. **The shipped Δ against the letter's `moment_A`** (phase A §4.3, under
   ingredient C): at the four super-bins the shipped `toy_delta_gluon` at
   10⁻² gives A₂ = −1.06 / −1.30 / −1.43 / −1.33 × 10⁻³ per unit P_zz — the
   opposite sign of `moment_A`'s +4.35 / +7.42 / +9.49 / +9.52 × 10⁻³ and 4–7×
   smaller. Options: keep the toy (a scenario with no CLI flag) / make the
   letter's model the generator default. Cost of a switch not measured.
   **None taken.**
6. **`F_mag`'s (q_z, b) = (1.30 fm⁻¹, 1.85 fm)**, `rc.hpp`'s "UNFITTED
   STARTING VALUES", assigned to a fit against Wiringa–Schiavilla's F_T
   landmarks: whether that fit should also answer to RFY Table I (row 6:
   p = 0.038 with a free scale, 0.0145 at s = 1; the free scale 1.419 ± 0.206
   matches the 1/0.7116 at which WS98 drew the same points). Price on any
   observable: not measured. **None taken.**
7. **Filing the fetched documents into `PolarizedLithiumSim/refs/`** (an
   author action; the session scratch directory is not permanent). Fetched
   2026-10-26, outside the repository: the BONuS Supplemental Material
   (`BONuSprcSupMat2.txt`, sha256 `18c7eb38…b3f65f0e`); Rand–Frosch–Yearian
   PR 144 (1966) 859 (`7cbb2d21…31d5bbe2`); Suelzle–Yearian–Crannell PR 162
   (1967) 992 (`ce366617…cc6e22d7`); Rodning–Knutson PRC 41 (1990) 898
   (`4dbfb75f…ee0ff161`); Bishop–Cheung PRA 20 (1979) 381 (`a28b358d…e4091950`);
   Close–Kumano PRD 42 (1990) 2377 (`105c9554…ec762ffb`); Sather–Schmidt PRD 42
   (1990) 1424 (`63378c24…439dc37a`); Moniz *et al.* PRL 26 (1971) 445
   (`2bdde956…5bc199e4`); SMC PRD 58 (1998) 112001 (`b6b4ca7b…c61eb195`); COMPASS
   PLB 769 (2017) 34, the journal's SCOAP3 copy (`de96c817…dd8d0751`); PDG 2026
   Table 1.1 (`f54f2e23…a9ac23c1`) and the QCD review (`c04c628d…2baa386f`);
   the DJANGOH 4.6.21 source and manual (`phase_B3_djangoh.md` §2, §8) and the
   `cteq61` set. Options: file them / leave them. **None taken** (no file was
   written into the sibling repository).
8. **The licence question.** The vendored numbers of three documents come
   from APS-copyrighted pages read through the institution's licence — the
   BONuS Supplemental Material's tables (`bonus2014_tkachenko_supplemental.json`)
   and Rand–Frosch–Yearian's Tables I, II, VII, VIII, XIV
   (`rfy1966_li_magnetic_ff.csv`, `rfy1966_li_magnetization.csv`) — vendored as
   facts with attribution; and the DJANGOH input cards vendored verbatim in
   `djangoh_eltail_ep18x275.json` come from repositories that state no licence
   (`spiesber/DJANGOH`, `eic/InclusiveDjangohSamples`). Options: keep them
   vendored / replace them by a pointer and a hash. **None taken.**
