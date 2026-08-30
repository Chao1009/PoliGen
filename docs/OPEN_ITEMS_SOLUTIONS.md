# Open items — solutions explored (2026-08-29)

Synthesis of four investigations (full reports in `docs/open_items/`:
`physics_literature.md`, `vmc_overlaps.md`, `code_designs.md`, `engineering.md`,
prototypes in `docs/open_items/prototypes/`). Each item below: what was found,
the recommended solution, effort, and status. Ordered by leverage.

| # | item | verdict | effort | status |
|---|---|---|---|---|
| 1 | VMC α+d / α+t cluster wave functions (plans/04 #15) | **closed & implemented** — ANL AV18 VMC tables (overlaps 2004; momentum distributions 2024) in `data/vmc/`, `VmcRadial` backend, `--cluster-wave vmc` | done | published tag fractions must be re-run |
| 2 | ePIC chain gate (HepMC3 → abconv → npsim) | **passed** — 10/10 events through `npsim` directly and via `abconv -p ip6_hiacc_100x10` | done | one cosmetic writer fix (electron generated_mass) |
| 3 | Tensor sign `TENSOR_LL_SIGN` (plans/08 D1) | **decided by literature: −1** (Cosyn Eq. 27, HERMES Eq. 6, HJM derivation, POLRAD Eqs. 9/10 all give A_zz = −(2/3) b₁/F₁) | 1 line + test | author to confirm; unblocks D2 |
| 4 | ⁶Li/⁷Li effective polarizations (plans/04 #6) | **closed** — 1/3 vs 0.81 was a convention mismatch; VMC (Wiringa 2014 Table I, Piarulli 2023) gives whole-nucleus P_p = P_n = 0.85 ± 0.03 (⁶Li), 0.87 / −0.03 (⁷Li) | 0.5 d | author to confirm |
| 5 | Coherent T2 final state | **solved in design** — PYTHIA accepts a Pomeron beam (id 990) through `LHAup`; γ*–IP → X reuses the whole bridge with ζ = β exactly; prototype: Q = 0, M_had = M_X, 40k ev/s | 2–3 d | prototype `pom_dis.cc` |
| 6 | Triton remnant (t* → N + …) | **solved in design** — Ciofi–Simula A=3 n₀(k) from BeAGLE; its norm 0.653 *is* the two-body/three-body split; current Hulthén is 2.7–5.5× too hard at 0.2–0.3 GeV | 2 d | prototype `nk_triton.py` |
| 7 | Spectator FSI (plans/04 #16) | **model chosen** — Glauber rescattering weight, per-nucleon product (Ciofi–Kaptari) with measured σ_tot(αp) = 121.5 mb, σ_el = 31.4 mb, B = 31 GeV⁻²; enters as a weight, never a momentum shift | 5–8 d | prototype `fsi_alpha.py` |
| 8 | Spin-3/2 SF basis (plans/04 #14) | **exists** — Jaffe–Manohar NPB 321 (1989); explicit J=3/2 functions arXiv:2209.12161 Eqs. 19a–d; rank-≤2 truncation is *exact* for unpolarized-beam inclusive observables | 5–10 d note | theory note |
| 9 | Tensor-sector RC (plans/04 #10) | **formulas exist** — POLRAD 2.0 tensor sector + elastic tail; Gakh–Shekhovtsova (uncited); precedent E12-13-011: 1.5 % | 5–8 d | band adoptable |
| 10 | b₁ for A > 2 (plans/04 #9) | **first-mover** — nothing exists; three-term α–d convolution on the Cosyn–Dong–Kumano–Sargsian kernel; 100 % band mandatory (⁶Li quadrupole puzzle) | 10–15 d | design only |
| 11 | Coherent ⁶Li amplitude (plans/04 #18) | **route changed** — Sartre ruled out (hard-coded nuclei, no polarization axis); use eSTARlight for unpolarized rates (1 d) and `hejajama/subnucleondiffraction` (code of arXiv:2408.13213) with an α+d configuration sampler for the tensor cos 2φ | 1 d / 10–15 d / collab | Mäntysaari-group ask |
| 12 | Packaging | scikit-build-core `pip install -e .` works (96 s); portable wheel needs `auditwheel` + GPL-3 terms | 1 d | prototype in scratchpad |
| 13 | License | **GPL-3.0-or-later** (forced by HepMC3/LHAPDF; matches MCnet norms) | 0 | author to confirm |

## 1. Cluster wave functions — the biggest physics correction

The genuine two-cluster VMC overlaps live on ANL's *older* page (`overlap_old/`:
`li6.ad`, `li7.at`, AV18+UIX, 2004, r- and k-space with MC errors); the
`momenta/` page (2024, AV18+UX) has α–d and α–t relative-momentum distributions
with an S/D block split. Both were fetched via the Wayback mirror (the ANL site is
behind a Cloudflare challenge). Findings:

- ⁶Li α+d: P_D = 1.9–2.0 % (the scenario value 0.0867 is 4× too high — it was
  tuned to reproduce a vector dilution of 0.87 that is itself a convention error,
  see item 4). S_αd = 0.82.
- The k-space density has a Pauli node (α–d 2S relative state) at k ≈ 0.15–0.2 GeV
  that no nodeless Hulthén form carries; VMC has *more* strength than Hulthén
  β = 0.30 in the 0.2–0.35 GeV Roman-Pot window and much less above 0.45 GeV.
- Reconciled (`validation/vmc_reconcile.py`, `docs/open_items/vmc_reconciliation.md`):
  the "β band biased low" claim in `physics_literature.md` is withdrawn — it evaluated
  the S-wave Hulthén for the P-wave ⁷Li channel. Against the right form VMC α+t is
  softer than every β in the far tail; ⁶Li VMC is harder at 0.2–0.3 GeV (P(k>0.2)
  0.148 → 0.241) and softer above 0.45. Both ANL file families agree bin-for-bin;
  production input = `momenta/` magnitudes + `overlap_old/` S–D sign.
- **Implemented**: `VmcRadial` backend, `ClusterWaveSource::VmcAV18`
  (`--cluster-wave vmc`); default stays Hulthén. Measured at 10×99.5: ⁶Li α tag
  0.0249 → **0.0348** (YR high-acceptance, ×1.40), 0.253 → 0.249 (tagging optics);
  ⁷Li 0.973 → 0.998. Tagged A_zz^tag(k) roughly halves (P_D 0.087 → 0.019 plus shape);
  the S–D interference sign flips below the S node at 0.134 GeV, currently outside
  every Roman-Pot acceptance.
- Consequence for the physics case: the ⁶Li α-tag acceptance is entirely a
  p_T-tail measurement, so the published tag fractions and the tagged A_zz
  curves must be re-run on VMC; the β-band of plans/04 #15 does not bracket
  the truth.

## 2. Chain gate — passed

`lipolgen-run … --hadronize` → `npsim --compactFile epic_craterlake_10x100.xml`
accepts the file directly (10-digit ion codes for beam ⁶Li status 4 and the α
spectator status 1); `abconv -p 1` cannot decode a ⁶Li ion (per-nucleon energy
0), `abconv -p ip6_hiacc_100x10` works and its output also runs through npsim.
EDM4hep `MCParticles` carry the full role chain. Cosmetic: write m_e = 0.511 MeV
as the electron `generated_mass` to silence DD4hep's ppm energy fix-ups.

## 3–4. Conventions now decided by sources

- `TENSOR_LL_SIGN`: four independent conventions (Cosyn 2410.12764 Eq. 27, HERMES
  hep-ex/0506018 Eq. 6, a two-line HJM parton-model derivation, POLRAD 2.0 Eqs.
  9/10) all give A_zz = −(2/3) b₁/F₁; the program's +1 is opposite, exactly as its
  own test `test_the_program_sign_is_opposite_to_the_literature` says. Flipping
  the constant also flips κ and the O(γ²) subtraction (D2).
- ⁶Li effective polarization: whole-nucleus convention (Cloët–Bentz–Thomas Eq. 24,
  no ÷Z): VMC 0.848 (Wiringa 2014 Table I), 0.86 chiral (Piarulli 2023);
  adopt 0.85 ± 0.03 (per-nucleon 0.283 if the slot divides by Z/N — state which).
  ⁷Li 0.866 / −0.037 confirmed. Do **not** derive tensor inputs from these wave
  functions (Q(⁶Li) is off by 2.5×).

## 5. Coherent T2 — Pomeron beam through the existing bridge

`Diffraction:doHard` cannot be driven externally and never fires for virtual
photons, but `Beams:idA = 990` is a legal LHAup beam (`BeamSetup.cc:870`). With
P_IP = P_ion − P_recoil (massless), ζ = β = Q²/(M_X²+Q²) exactly; the Pomeron
remnant is a single antiquark, so charge and colour close automatically. Traps
found: colour-tag orientation for an incoming antiquark (else 50 % silent veto);
LO Pomeron grids have no quarks at Q² = 1, β < 0.1 (need the `q2_pdf_min` floor).
Raise `COHERENT_MX_MIN_DEFAULT` to 1.2 GeV; below it, an exclusive-VM channel
(ρ, φ, J/ψ) is a separate Phase-2 item. Fallback: hand-filled qq̄ string with
`ProcessLevel:all = off` (130k ev/s, 1e-15 conservation).

## 6. Triton remnant — Ciofi–Simula A = 3

BeAGLE's `DT_KFERMI` carries n₀(k) = Σᵢ Aᵢ e^{−Bᵢk²}/(1+Cᵢk²)² with the CS
coefficients (A=3: 31.7/1.32/5.98 + 0.00266/0.365/0); its norm ∫n₀k²dk = 0.653 is
the ground-state-remnant (two-body) fraction — the split the breakup needs.
Design: `TritonSpectralFunction` backend with k-dependent branching
p₂(k) = n₀/(n₀+n₁), a third channel (struck n → (pn) continuum, absent today),
fixed uniform consumption per draw. Also recorded: BeAGLE renormalizes n₀ to 1 and
therefore drops the 35 % continuum for A = 3 (its tail is too soft by construction).

## 7. FSI — weight, not a shift

Eikonal rescattering transfers transverse momentum only, so the distortion is a
function of k_T at fixed k_z; pole dominance = small k_T, which is exactly the
Roman-Pot α-tag selection. Prototype (σ_XN = 40 mb): ≤ 20 % for k ≲ 0.02 GeV,
~27 % at 0.05, 40–60 % at 0.10; σ_XN(W) formation-length ramp is the weakest
input and must be banded (the 20 mb row is arguably realistic at EIC). Use σ_tot
in the absorptive term and σ_el in the gain term (74 % of α rescatterings destroy
the tag). Hook: `FsiWeight` interface on `TaggedSampler`, multiplied into
`Event::weight`; null = today's PWIA bit-for-bit.

## 8–10. Theory notes

- Spin-3/2: the complete basis is Jaffe–Manohar (1989); 2209.12161 writes the four
  leading-twist J = 3/2 functions (their "g₂" is the rank-3 partner of g₁ — rename
  `g1_rank3`); Cosyn–Weiss 2603.23699 App. D is the on-ramp. Time reversal makes
  odd-l multipoles unobservable with an unpolarized beam in inclusive DIS, so the
  generator's rank-≤2 truncation is a theorem, not an approximation. Gap: the
  finite-γ inclusive decomposition for J = 3/2 (5–10 d, ideal Cosyn/Weiss co-authorship).
- Tensor RC: adopt ISR shift (spin-blind) + POLRAD Eq. (A.4) tensor elastic tail
  with VMC ⁶Li form factors (Wiringa–Schiavilla 1998); quote 1.5 % (x ≳ 0.05) and a
  10–30 % band on the tensor part at x ≲ 0.01 (Gakh–Shekhovtsova, uncited).
- b₁(⁶Li): three-term α–d convolution (embedded deuteron b₁ ⊗ f_{d/Li}(z), α–d
  D-wave term with F₁ᵈ, CG depolarization); validate on A = 2 (Cosyn Figs. 4/5)
  first; 100 % band.

## 11. Coherent ⁶Li amplitude

Sartre: nuclei hard-coded (`Nucleus.cpp` switch, no A = 6), spherical sampling
(no polarization axis), tables CPU-years (the 2026 speed-up code is unreleased) —
ruled out. Outdated project claims: 2605.00454 publishes coherent J/ψ down to
A = 3, 4 with α-clustering; 2511.05638 has e+⁷Li coherent J/ψ |t| distributions
with tagging efficiency. Path: eSTARlight ⁶Li (1 d) for rates; then an α+d
configuration sampler (α core from VMC density, p–n pair from AV18 u/w oriented
by the polarization axis, α–d separation from the VMC momentum distribution)
grafted into `subnucleondiffraction` (author asks to be contacted). Warn in
plans/06: photon-polarization cos 2φ (STAR 2204.01625) is a distinct mechanism
and a background at Q² > 0.

## 12–13. Engineering

- `pyproject.toml` + scikit-build-core builds and installs in 96 s with
  `--config-settings=cmake.define.LIPOLGEN_DEPS_PREFIX=…`; RPATH is absolute
  (non-relocatable); a shippable wheel needs `auditwheel` vendoring and, with
  it, GPL-3 distribution terms; data (xmldoc, LHAPDF grids) not vendored.
- License: GPL-3.0-or-later for LiPolGen's own code (the linked combination is
  GPL-3 regardless; permissive headers would mislead; MCnet-consistent).
