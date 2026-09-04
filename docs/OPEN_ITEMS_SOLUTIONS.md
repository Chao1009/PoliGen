# Open items — solutions explored (2026-08-29)

Synthesis of four investigations (full reports in `docs/open_items/`:
`physics_literature.md`, `vmc_overlaps.md`, `code_designs.md`, `engineering.md`,
prototypes in `docs/open_items/prototypes/`). Each item below: what was found,
the recommended solution, effort, and status. Ordered by leverage.

| # | item | verdict | effort | status |
|---|---|---|---|---|
| 1 | VMC α+d / α+t cluster wave functions (plans/04 #15) | **closed & implemented** — ANL AV18 VMC tables (overlaps 2004; momentum distributions 2024) in `data/vmc/`, `VmcRadial` backend, `--cluster-wave vmc` | done | re-run on VMC at 5×41, 10×100, 18×275 (USAGE cluster-wave section; vmc_reconciliation.md "Impact on the tagged pipeline"); β band retired |
| 2 | ePIC chain gate (HepMC3 → abconv → npsim) | **passed** — 10/10 events through `npsim` directly and via `abconv -p ip6_hiacc_100x10` | done | writer fix in (m_e written as generated_mass; pinned by tests) |
| 3 | Tensor sign `TENSOR_LL_SIGN` (plans/08 D1) | **decided by literature: −1** (Cosyn Eq. 27, HERMES Eq. 6, HJM derivation, POLRAD Eqs. 9/10 all give A_zz = −(2/3) b₁/F₁) | 1 line + test | author to confirm; unblocks D2 |
| 4 | ⁶Li/⁷Li effective polarizations (plans/04 #6) | **closed** — 1/3 vs 0.81 was a convention mismatch; VMC (Wiringa 2014 Table I, Piarulli 2023) gives whole-nucleus P_p = P_n = 0.85 ± 0.03 (⁶Li), 0.87 / −0.03 (⁷Li) | 0.5 d | author to confirm |
| 5 | Coherent T2 final state | **implemented 2026-09-01** — γ*–Pomeron tier is the default coherent T2 (`CoherentT2::{Pomeron, Off}`, `docs/PYTHIA_BRIDGE.md` §12); ζ = β exact, 300-event chain conserves to 2.8e-14, M_had = M_X to 2.3e-11, veto 0 at M_X ≥ 1.4 GeV | done | `--coherent-t2`, `--pom-set` |
| 6 | Triton remnant (t* → N + …) | **implemented 2026-09-01** — `triton_sf.hpp` CS n₀+n₁ model, S₀ = 0.6525 untuned, third channel n → (pn); opt-in `--triton-sf ciofi-simula`, sequential Hulthén stays the default bit for bit | done | numbers in §6 below |
| 7 | Spectator FSI (plans/04 #16) | **implemented 2026-09-01** — `GlauberFsiWeight` per-event weight on `Event::weight` (`--fsi`, tagged channels; never a momentum shift); σ_Xα = 131.0 mb at σ_XN = 40, band 20–40 mb mandatory | done | numbers in §7 below |
| 8 | Spin-3/2 SF basis (plans/04 #14) | **implemented as a theory note** — Jaffe–Manohar NPB 321 (1989); explicit J=3/2 functions arXiv:2209.12161 Eqs. 19a–d; rank-≤2 truncation is *exact* for unpolarized-beam inclusive observables; the finite-γ J = 3/2 decomposition and the full list of code changes needed if the rank-3 sector is ever switched on are written up, with no default behaviour change | done | `docs/theory/SPIN32_FINITE_GAMMA.md` |
| 9 | Tensor-sector RC (plans/04 #10) | **implemented 2026-09-03** — `rc.hpp`/`rc.cpp`, opt-in `--rc tensor-band`: the two-sided band `rc_tensor_lo`/`rc_tensor_hi` on the tensor part of the rate (δ log-linear, 0.30 at x = 0.01 → 0.015 at x = 0.16) plus `rc_tail`, POLRAD's t-peak elastic tail with its tensor part and the unpolarised quasi-elastic tail. Weight-only, on `Event::rc_weights` and never on `Event::weight`; `--rc off` is byte-identical | done | numbers in §9 below |
| 10 | b₁ for A > 2 (plans/04 #9) | **implemented 2026-09-03 as an OPT-IN backend; its A = 2 gate CLOSED the same day** — `b1_nuclear.hpp`/`b1_nuclear.cpp`, `--b1-model li6-convolution`: the **four**-term α–d convolution on the Cosyn–Dong–Kumano–Sargsian kernel (the struck-α orbital term is not optional, it is ≈ 0.5 × the struck-d one). The default stays `Li6B1(MillerB1)`, bit for bit. **The A = 2 validation gate PASSES**: read where design §5.4's checklist ends — MSTW2008 LO (CDKS's own PDF, now reachable as `MstwSF` from PYTHIA's own grid) at CDKS Eq. (21)'s δ-function — G3a passes on all three landmarks and G3b's peak ratio is **0.843243**, inside the factor-2 window; with the real CD-Bonn wave function as well, **1.000338**. Four conditions travel with that: **the lift is about a configuration, not about a build** — the shipped default `ToyF2` gives **0.440**, outside the window, so quote numbers made with `--b1-unpol mstw` and not the default backend's; it is marginal at Eq. (17)'s κ = 1 (0.520); the CD-Bonn agreement is below the reference's own digitization error and is not three-digit agreement with CDKS; and **the gate is A = 2** — so the **mandatory ±100 % band stays**. G3a's own pass is **qualified** by its counting ceiling (§10, "G3a's stated limitation"). The ⁶Li publication ban is lifted. Of the seven close conditions five are done and two are author decisions (Miller's normalisation; the ⁶Li R default) | 10–15 d, ~9 spent | **gate closed**, band mandatory — §10 below |
| 11 | Coherent ⁶Li amplitude (plans/04 #18) | **both halves in-tree** — eSTARlight ⁶Li unpolarized rates/slope (2026-09-02, `estarlight_li6.md`) settle `slope_b = 50 ± 10` GeV⁻² as citable and put a **lower bound** of 1.0e-3 on `f0` (the all-VM 3.0e-2 is not an upper bound — ρ/φ sit below the M_X floor); both recorded in `coherent.hpp`; the α+d configuration sampler (2026-09-03, `cluster_config.hpp`, `phase_G_numbers.md`) predicts ⁶Li's tensor a₂ ~10× smaller than the deuteron's and of **opposite sign** | done / done / collab | §11 below; Mäntysaari-group ask drafted, not yet sent |
| 12 | Packaging | **implemented 2026-09-02** — `pyproject.toml` (scikit-build-core) in-tree, `pip install -e .` works (66 s); one copy of each `.so` in `lipolgen/`, `$ORIGIN`+deps-prefix RPATH, data/vmc vendored; portable wheel still needs `auditwheel` + GPL-3 terms | done | see §12–13 below |
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
  curves were re-run on VMC (`validation/vmc_tag_fractions.py`,
  `docs/open_items/vmc_reconciliation.md` "Impact on the tagged pipeline",
  `docs/USAGE.md`): ⁶Li tag fraction 0.0249 → 0.0348 (YR high-acceptance),
  0.2530 → 0.2485 (tagging optics); ⁷Li 0.9730 → 0.9981; ⁶Li A_zz^tag at
  k = 0.20 GeV +0.845 → +0.452; P_D 0.0867 → 0.01935. Because no single β
  reproduces the VMC shape (Pauli node + window-dependent tail), the Hulthén
  β band is retired as the cluster-wave systematic rather than widened — the
  `--cluster-beta` knob itself stays (Hulthén is still the default radial
  form), but the quoted systematic is now the VMC-vs-Hulthén difference
  shown above, not a scan over β. Confirmed at the other two reference
  configurations (`validation/vmc_tag_fractions.py --configs 0,1,2`,
  `docs/USAGE.md`): the ⁶Li YR high-acceptance tag fraction moves
  0.0286 → 0.0365 at 5×41 and 0.0266 → 0.0349 at 18×275, so the
  Hulthén→VMC shift is not a single-energy artifact.

## 2. Chain gate — passed

`lipolgen-run … --hadronize` → `npsim --compactFile epic_craterlake_10x100.xml`
accepts the file directly (10-digit ion codes for beam ⁶Li status 4 and the α
spectator status 1); `abconv -p 1` cannot decode a ⁶Li ion (per-nucleon energy
0), `abconv -p ip6_hiacc_100x10` works and its output also runs through npsim.
EDM4hep `MCParticles` carry the full role chain. Cosmetic: write m_e = 0.511 MeV
as the electron `generated_mass` to silence DD4hep's ppm energy fix-ups.

**Done (2026-09-02)**: `src/hepmc/hepmc_writer.cpp` writes
`generated_mass = 0.51099895e-3` GeV for any massless electron (`|pdg| == 11`,
`p.mass == 0.0`); every other particle's generated mass is unchanged. Pinned
by a test in the doctest suite.

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

**IMPLEMENTED (2026-09-01), measured.** The design above is the shipped tier
(`docs/PYTHIA_BRIDGE.md` §12): `make_coherent` writes `Role::Pomeron`
(P_IP = P_ion − P_recoil, pdg 990, status 3), the bridge's third instance
(`Beams:idA = 990`, behind `PythiaBridgeOptions::coherent_t2`) hadronizes it
with W² → M_X², ζ = β verified exact to 10 digits, flavour drawn
e_q² × f_q(β, Q²) from `PDF:PomSet` (default 6) with the e_q² fallback. Both
predicted traps were hit and closed (antiquark colour-tag orientation;
quark-free LO grids → `q2_pdf_min` floor). 300-event chain: conservation
2.8e-14 relative, charge exact, recoil untouched, M_had = M_X to 2.3e-11,
veto 0 at M_X ≥ 1.4 GeV (the measured veto table sits in `coherent.hpp` and
set `COHERENT_MX_MIN_DEFAULT = 1.2`). One extra find: PYTHIA closes the
meson-like Pomeron-beam record on the electron beam's massive light-cone
minus — a constant Δ(m²) = −3.9e-6 GeV² the baryon-beam remnant path absorbs —
so the surrogate is built at a compensated w2_sur, keeping the mass-repair
rescale at ~1e-11. `PDF:PomSet` remains the tier's largest systematic.

## 6. Triton remnant — Ciofi–Simula A = 3

BeAGLE's `DT_KFERMI` carries n₀(k) = Σᵢ Aᵢ e^{−Bᵢk²}/(1+Cᵢk²)² with the CS
coefficients (A=3: 31.7/1.32/5.98 + 0.00266/0.365/0); its norm ∫n₀k²dk = 0.653 is
the ground-state-remnant (two-body) fraction — the split the breakup needs.
Design: `TritonSpectralFunction` backend with k-dependent branching
p₂(k) = n₀/(n₀+n₁), a third channel (struck n → (pn) continuum, absent today),
fixed uniform consumption per draw. Also recorded: BeAGLE renormalizes n₀ to 1 and
therefore drops the 35 % continuum for A = 3 (its tail is too soft by construction).

**IMPLEMENTED (2026-09-01), measured** (`triton_sf.hpp` / `triton_sf.cpp`,
opt-in via `BreakupOptions::triton_sf` / `PipelineConfig::triton_sf` /
`--triton-sf ciofi-simula`; the sequential Hulthén model stays the default
bit for bit). S₀ = 0.652548 untuned against the CS 0.6525; n₀ tail
P(k > 0.1/0.2/0.3/0.45 GeV) = 0.4267 / 0.06276 / 0.009621 / 0.002365 (all at
the transcription's own digits); sampled bound fraction 0.6524 at 10⁵ draws;
⟨k⟩ = 102 MeV on the bound-d channel, 126 MeV over all struck nucleons
(sequential model: 133/145). The (pn) continuum pair splits at the pn ¹S₀
pole `KAPPA_PN_SINGLET` beside the nn one; sample() consumes exactly 6
uniforms whatever the branch and the breakup exactly 9, verified empirically.
The Pipeline builds the model itself at the run's own `cluster_beta`, so no
physics number is defined twice. T1+T2 chain with the hadronizer on
conserves < 1e-9 on both options (`tests/test_t2.cpp`).

## 7. FSI — weight, not a shift

Eikonal rescattering transfers transverse momentum only, so the distortion is a
function of k_T at fixed k_z; pole dominance = small k_T, which is exactly the
Roman-Pot α-tag selection. Prototype (σ_XN = 40 mb): ≤ 20 % for k ≲ 0.02 GeV,
~27 % at 0.05, 40–60 % at 0.10; σ_XN(W) formation-length ramp is the weakest
input and must be banded (the 20 mb row is arguably realistic at EIC). Use σ_tot
in the absorptive term and σ_el in the gain term (74 % of α rescatterings destroy
the tag). Hook: `FsiWeight` interface on `TaggedSampler`, multiplied into
`Event::weight`; null = today's PWIA bit-for-bit.

**IMPLEMENTED (2026-09-01), measured** (`fsi.hpp` / `fsi.cpp`,
`TaggedSampler::set_fsi` → `TaggedEvent::weight` → `Event::weight`;
`PipelineConfig::fsi` / `--fsi {off,glauber-cluster,glauber-nucleon}` +
`--fsi-sigma-mb`, tagged channels only, `validate()` refuses it elsewhere).
Profile numbers at σ_XN = 40 mb: σ_Xα = 131.0 mb (Glauber-shadowed, not
4 × 40), σ_el = 35.2 mb, B_α = 27.2 GeV⁻²; at 20 mb: 72.5 mb. The pinned
FSI/IA table reproduces the `fsi_alpha.py` prototype to ≤ 3.3e-3 (0.787 at
k = 0.02; 0.607/0.381 at k = 0.10, θ = 0/90°; 0.398 at 0.20; the 20 mb row
0.877/0.844/0.764), production S+D channel 0.617 at k = 0.10, θ = 0;
survival 0.517 (cluster) / 0.582 (nucleon variant) on the ⁶Li α tag at
40 mb. Exact θ → π−θ symmetry, ratio → 1 as σ → 0, the wrong-spectator
guard, the ⁷Li P-wave channel and the σ_XN(W) formation ramp are all in
`tests/test_fsi.cpp`; the weight reaches HepMC3 `weights()[0]` and the npz
`weight` column, and every four-vector stays bit-identical to the PWIA run.
Quote it as an unpolarized-shape systematic banded over 20–40 mb, never as
a correction to A_zz.

## 8–10. Theory notes

- Spin-3/2: the complete basis is Jaffe–Manohar (1989); 2209.12161 writes the four
  leading-twist J = 3/2 functions (their "g₂" is the rank-3 partner of g₁ — rename
  `g1_rank3`); Cosyn–Weiss 2603.23699 App. D is the on-ramp. Time reversal makes
  odd-l multipoles unobservable with an unpolarized beam in inclusive DIS, so the
  generator's rank-≤2 truncation is a theorem, not an approximation. Gap: the
  finite-γ inclusive decomposition for J = 3/2 (5–10 d, ideal Cosyn/Weiss co-authorship).
  → **Written up 2026-09-03 as `docs/theory/SPIN32_FINITE_GAMMA.md`**: the
  rank-≤2 truncation is proved exact for unpolarised-beam observables (§4), the
  finite-γ J = 3/2 decomposition is derived there (§5 — two rank-3 structure
  functions, master formula (46)–(47); the repository's own derivation, still
  unpublished in the literature, with the Appendix listing every unsourced
  claim), and §6 lists the code changes if the rank-3 sector is ever switched
  on. **No code behaviour changed.**
- Tensor RC: adopt ISR shift (spin-blind) + POLRAD Eq. (A.4) tensor elastic tail
  with VMC ⁶Li form factors (Wiringa–Schiavilla 1998); quote 1.5 % (x ≳ 0.05) and a
  10–30 % band on the tensor part at x ≲ 0.01 (Gakh–Shekhovtsova, uncited).
  → **Implemented as §9, with three deliberate departures from this note**: the
  ISR *shift* is **not** applied (it would cost an extra uniform per event and
  break the "RC moves nothing" invariant — its residual is what the band
  prices); the tail is built on POLRAD's **t-peak closed forms Eqs. (37)–(39),
  (43)**, not the Eq. (A.4)/Appendix-B machinery, because POLRAD §2.1.3 B says
  the s- and p-peaks are *suppressed* for a tail; and the ⁶Li form-factor
  **normalisations are the measured moments, not VMC** (whose Q(⁶Li) is 3× the
  measured one). The 1.5 % anchor also moved to x = 0.16, E12-13-011's own
  lower kinematic edge.
- b₁(⁶Li): three-term α–d convolution (embedded deuteron b₁ ⊗ f_{d/Li}(z), α–d
  D-wave term with F₁ᵈ, CG depolarization); validate on A = 2 (Cosyn Figs. 4/5)
  first; 100 % band.
  → **Implemented as §10, with one structural correction to this note: it is
  FOUR terms, not three.** CDKS Eq. (10) sums the spectral function over
  *constituents*, and one level up that sum runs over {d, α}. The α is J = 0 so
  b₁^α ≡ 0, but its light-cone density carries the *same* (3cos²θ − 1) orbital
  alignment as the deuteron's — (3c² − 1) is even in k⃗ and the α carries −k⃗ —
  so the **struck-α orbital term (2α) exists and is ≈ 0.5 × the struck-d one
  (2d)**, fixed by counting (4/6 against 2/6) and by the 1/M² of the
  P₂-weighted density: 2(M_d/M_α)² = 0.5064 against a measured **0.494–0.504**
  (0.489–0.497 on the design's coarser y grid — a quadrature artefact, since
  fixed by the 2400/2400/3200 default; `phase_D_numbers.md`).
  A regression that silently drops it moves b₁ by ~30 %. The A = 2 validation
  was done and it **passes on magnitude in one configuration** — MSTW2008 LO
  as the unpolarised nucleon input (`--b1-unpol mstw`) at CDKS Eq. (21)'s
  δ-function, G3b = 0.843243. On the shipped default `ToyF2` the same clause
  is **0.440, outside** the [0.5, 2] window. See §10, which states the
  configuration and the limitation on G3a's counting window with it.

## 9. Tensor-sector RC — a band and a background, never a shift

**IMPLEMENTED (2026-09-03), measured** (`include/lipolgen/rc.hpp`,
`src/core/rc.cpp`; `PipelineConfig::rc` / `--rc {off,tensor-band}` plus
`--rc-delta-low-x`, `--rc-delta-high-x`, `--rc-fq-scale`,
`--rc-tail-tensor-scale`, `--rc-qe-suppression`; `docs/USAGE.md` §7b; the raw
run in `docs/open_items/run_2026-09-02/phase_C_numbers.md`).

> **Revised 2026-09-03 after review.** Every tail number below was
> regenerated: the elastic tail's per-nucleon reduction was **6× too large**
> (Eq. (38) is the whole-nucleus `d²σ/dx_A dy`, so the factor is `1/A²` and not
> `m_p/M_A`), the quasi-elastic **Pauli suppression `S(q)` is now on by
> default** (design Q9, closed), the nuclear map moved to this library's exact
> `x_A = x/A`, and the band is now **clamped** (`RcOptions::band_tau_max`)
> because `tau_tag` diverges at the nodes of the tagged spectator density and
> was publishing negative weights.

Three weights per event — `rc_tensor_lo = 1 − δ(x)τ`,
`rc_tensor_hi = 1 + δ(x)τ`, `rc_tail = 1 + σ_tail/σ_Born` — on
`Event::rc_weights` and **never** on `Event::weight`: the band is a systematic
variation and the tail a background, so neither may touch the Born sample.
No four-vector moves, no RNG is consumed (`RcModel::fill` takes no `Rng&`),
and an `--rc off` npz and HepMC3 file are **byte-identical** to what the same
seed wrote before `rc.hpp` existed. Applies on every channel except coherent
⁶Li (a φ-dependent tensor observable — nothing to cite); the tail is
identically 1 on the tagged channels, where the tag itself vetoes the elastic
recoil at x_L = 1.

**Configuration for every number below.** ⁶Li, `--config 1`
(s per nucleon = 3980 GeV²), `--channel inclusive`,
`--plan tensor-thirds --pzz 0.6`, θ_S = 0, defaults, `n_eta = 128`, tail node
grid 101 × 77 in (ln x, ln y).

**The band, at P_zz = +1, Q² = 5 GeV²** (`w_lo + w_hi == 2.0` bit-for-bit):

| x | δ(x) | A_zz | τ | `w_hi − 1` |
|---|---|---|---|---|
| 0.010 | 0.30000 | −1.4735e−03 | −7.3728e−04 | −2.2119e−04 |
| 0.063 | 0.11081 | −4.2848e−03 | −2.1470e−03 | **−2.3790e−04** |
| 0.100 | 0.06331 | −4.9751e−03 | −2.4937e−03 | −1.5789e−04 |
| 0.160 | 0.01500 | −4.1666e−03 | −2.0876e−03 | −3.1315e−05 |
| 0.300 | 0.01500 | −1.5852e−04 | −7.9268e−05 | −1.1890e−06 |

The band **peaks near x = 0.063, not at the lowest x**: δ(x) is still rising
there while |A_zz| has not yet fallen. The largest RC systematic on `A_zz`
sits in the *middle* of the low-x range.

**The tails, `w_tail − 1`** (elastic + unpolarised quasi-elastic, as a
fraction of the Born, at m = ±1):

| x | Q² = 2 | Q² = 5 | Q² = 10 |
|---|---|---|---|
| 0.01 | 8.9654e−05 | 5.18398e−04 | 2.14577e−03 |
| 0.10 | 2.18453e−07 | 1.32298e−06 | 5.19642e−06 |
| 0.30 | — (y < 0.004) | 6.8624e−08 | 2.72302e−07 |

**and the headline that is not what a fixed-target intuition expects:** at EIC
collider kinematics these cells sit at small y (`y = Q²/(x s)`), where
`Y₊ = [1+(1−y)²]/(1−y) ≈ 2`. The radiative-tail dilution is therefore **five
to eight orders of magnitude smaller than at HERMES**. The same model run at
HERMES-like y gives 1.8 % at y = 0.5 and 33 % at y = 0.9. **`rc_tail` is a
small effect at the EIC and a large one at a fixed target, and the difference
is entirely `Y₊` and the Born's 1/Q⁴.**

**And `rc_tail` is a LOWER BOUND: it is the t-peak ALONE.** T8(c)
(`tests/test_rc.cpp`) now measures the leading-log s-/p-peaks of the same
observable with an independent Weizsäcker–Williams × elastic construction.
POLRAD §2.1.3 B's "the s- and p-peaks are suppressed" holds **where this
generator runs** — > 99 % of the total for ⁶Li at Q² ≥ 20 GeV² — and fails
elsewhere: at the HERMES deuteron point (x = 0.012, y = 0.85, Q² = 0.53) the
t-peak is only **23 %** of the total, low by a factor **4.4**, and at the
low-Q² corner of the generator window (Q² ≈ 4 GeV²) the quasi-elastic s-peak
is already 3.3× the t-peak. Never quote `rc_tail` as "the" radiative tail.

**The derived quantities that settle three design estimates:**

| x | Q² | (1/6)σ^el_T/σ^el_U | σ^q_U/σ^el_U (S(q) on / S = 1) | ΔA_zz (tail) |
|---|---|---|---|---|
| 0.01 | 5 | −5.0943e−04 | 0.28237 / 0.5981 | +3.5167e−07 |
| 0.10 | 5 | +1.5595e−04 | 2.7475 / 3.164 | +6.6751e−09 |
| 0.30 | 5 | +1.3287e−02 | **1000.2** / 1000.2 | +1.2665e−11 |

1. **The tensor fraction of the elastic tail changes SIGN with x** (−5.1e−04
   at x = 0.01, +1.3e−02 at x = 0.30), which is why `σ^el_T` is a separate
   table and never a scale factor on `σ^el_U`.
2. **`σ^q_U/σ^el_U` runs 0.28 → 2.75 → 1000** over x = 0.01 → 0.30 (0.60 →
   3.16 → 1000 with the Pauli suppression off). The design's "30–70 % of the
   ERT" estimate is right only at the bottom of the range: **`rc_tail` is
   quasi-elastic-dominated at every x ≳ 0.03** (22 % of the tail at x = 0.01,
   73 % at 0.1, **99.9 %** at 0.30), so `--rc-qe-suppression` and
   `RcOptions::qe_kf_gev` are the dominant tail knobs almost everywhere.
3. **The tail is not the leading RC systematic at the EIC.** ΔA_zz from the
   whole tail is 3.5e−07 at x = 0.01, Q² = 5, against a **band half-width of
   4.4e−04** — a factor 1300.

**Run-level:**

| quantity | value |
|---|---|
| `RcModel` construction (101 × 77 nodes, `n_eta = 128`) | **0.106 s** — the `Pipeline` setup difference, median of 7. ONE number, shared with `phase_C_numbers.md` §8.3 and `USAGE.md` §7b; the earlier 0.16 s / 0.134 s pair is withdrawn. |
| inclusive throughput, `--rc off` → `--rc tensor-band` (1 core, end to end) | 574 852 → 495 017 ev/s (**−13.9 %**), same build and machine as the row above |
| clipped tail nodes at `tail_max = 10`, global / y < 0.5 / 0.5–0.9 / y > 0.9 | **0 / 0 / 0 / 0** (was 1.71 % / 0 / 0 / 21.95 % before the per-nucleon fix) |
| clipped EVENTS — tail / band, default 2000-event inclusive run | **0 / 0** (`meta["rc_clipped_tail_event_fraction"]`, `..._band_...`; a DIFFERENT quantity from the node fractions — nodes are not event-weighted) |
| clipped EVENTS on the BAND, 20 k tagged-alpha | **0.62 %** (⁶Li), **2.6 %** (⁷Li). `tau_tag` reaches 30.7 at the M = 0 density nodes; `band_tau_max = 1` holds every edge in [0.7, 1.3] instead of the −1.79 / −8.35 v0 published. |
| δ(A_zz) at x = 0.01, Q² = 5 from the band, δ_low = 0.19 vs 0.30 | 2.7996e−04 vs **4.4204e−04** |
| … from `fq_scale` 0 vs 2 | ΔA_zz +7.54154e−07 → −4.9925e−08, spread 8.041e−07 |
| … from `tail_tensor_scale` 0.5 vs 2 | +3.59193e−07 → +3.23118e−07, spread 3.61e−08 |
| … from `qe_suppression` 0 vs 1 | +1.83644e−07 → +3.51674e−07, spread 1.680e−07 |
| … from `qe_kf_gev` 0 (S = 1) vs 0.169 | +5.39491e−07 → +3.51674e−07, spread 1.878e−07 |

**The `fq_scale` band changes the SIGN of ΔA_zz** and is not symmetric about
the nominal — σ^el_T is *quadratic* in F_q. That is why it must be **RUN**
(0, 1, 2) and never rescaled from one row, and why T12 fits a quadratic
through three runs instead of asserting linearity. **Nothing clips on the tail
any more** — dividing the elastic tail by the missing A = 6 removed the whole
21.95 % that used to sit in the `y > 0.9` band — but the four node fractions,
the two per-event fractions and `Event.rc_clipped` are all still reported,
because `Y₊ ~ 1/(1−y)` is still where a wider `y` window would bite. What does
clip is the **band on the tagged channels**, 0.6–2.6 % of events, where
`tau_tag = 1 − n̄/n_M` diverges at the nodes of `n_M`.

**Honest flags, mandatory wherever any of this is quoted.** `δ_low = 0.30` is
the **size of a correction this generator does not apply**, taken as a 1σ
band; its source (Gakh–Shekhovtsova, hep-ph/0403262) has **zero INSPIRE
citations**. `RC_DELTA_LOW_X_OPTIMISTIC = 0.19` is HERMES's *measured*
fractional residual at its lowest-x bin (2×10⁻³ on A_zz = −1.06×10⁻²).
Everything cited is **deuteron** — whether the deuteron's fractional RC
transfers to ⁶Li is the largest unquantified assumption here. The ⁶Li
form-factor shape parameters `(a, α, q₀) = (1.9069 fm, 0.13822,
3.0999 fm⁻¹)` and `(q_z, b) = (1.30 fm⁻¹, 1.85 fm)` are **unfitted starting
values** (the Suelzle–Yearian–Crannell and Li–Sick–Whitney–Yearian elastic
data are not in this repository in machine-readable form — **Q1 and Q10 stay
open**); their *normalisations* are the **measured** moments (μ = +0.822047
μ_N, Q = −0.0818 fm², TUNL A = 6) and never VMC, whose Q(⁶Li) = −0.23(9) fm²
is 3× the measured one. **No RC calculation exists for a tagged tensor
asymmetry** (the tagged band is an uncited extrapolation) and none for any
φ-dependent tensor observable at any axis. The **polarised** quasi-elastic
tail is not priced at all (Zhou *et al.*, PRL **82** (1999) 687), and at
x = 0.30 the tail is **99.9 %** quasi-elastic. **`rc_tail` is the t-peak only
and is a LOWER BOUND** (T8(c): fine for ⁶Li at Q² ≥ 20 GeV², low by 4.4× at
the HERMES deuteron point). The **tagged band is clamped** at
`|τ| ≤ band_tau_max = 1`, which is a *choice*: `n_M → 0` is exactly where the
fractional-rescale ansatz breaks down, and the clipped fraction is reported
rather than hidden.

## 10. b₁ of ⁶Li — a four-term α–d convolution, opt-in, and the A = 2 gate now PASSES

**IMPLEMENTED (2026-09-03) AS AN OPT-IN BACKEND; ITS A = 2 GATE CLOSED THE SAME
DAY; THE ITEM ITSELF IS ALL BUT CLOSED — two author decisions and one
unquantified systematic remain, listed at the end**
(`include/lipolgen/b1_nuclear.hpp`, `src/core/b1_nuclear.cpp`;
`PipelineConfig::b1_model` / `--b1-model {miller,cdks,li6-convolution}` plus
`--b1-band-scale` and `--b1-alpha-d-dwave-weight`; `docs/USAGE.md` §2a; design
`docs/open_items/run_2026-09-02/design_D_b1_li6.md`; the gate as it now stands
in `docs/open_items/run_2026-09-03/phase_A_numbers.md` and `phase_A_cdbonn.md`;
the superseded failing verdict in `phase_D_gate.md`; the ⁶Li numbers in
`phase_D_numbers.md`).

> ### ⚠ READ THIS BEFORE QUOTING ANY NUMBER BELOW
>
> Design §5 makes the A = 2 validation gate **blocking**: the same kernel, fed
> the AV18 deuteron u(k), w(k) instead of the α–d waves, must reproduce the
> digitized CDKS Fig. 4 before any ⁶Li number is quoted. **Since 2026-09-03 it
> does.** Read at the end of §5.4's checklist — MSTW2008 LO, CDKS's own nucleon
> PDF (item 4), at CDKS Eq. (21)'s δ-function — G3a passes on all three
> landmarks and **G3b's peak ratio is 0.843243**, inside the factor-2 window
> with the lower edge cleared by a factor 1.69. With the real CD-Bonn wave
> function (item 5) as well it is **1.000338**. Design §5.4's *Escalation*
> clause is not triggered, so **the ⁶Li publication ban is LIFTED** — for that
> configuration. The tables below are **not** in it; see condition 1 and
> "The ⁶Li numbers" heading further down.
>
> **Four conditions travel with that, and none of them is optional.**
>
> 1. **The lift is about a CONFIGURATION, not about a build.** The gate passes
>    for the **MSTW2008 LO** nucleon input at CDKS Eq. (21)'s δ-function. On
>    the **shipped default** unpolarised backend, the library `ToyF2`, the same
>    gate gives **0.440** — *outside* the [0.5, 2] acceptance window. That row
>    is still pinned in the doctest as a regression guard and is *not* the
>    verdict: §5.4's Escalation clause reads the ratio *after* the checklist,
>    and item 4 is the nucleon PDF. Since 2026-09-04 the passing configuration
>    is emittable from the run surface — `PipelineConfig::b1_unpol` /
>    `--b1-unpol mstw` threads `MstwSF` into `Li6ConvolutionOptions::unpol`
>    through `default_inclusive_kernel` — so the rule is **quote numbers made
>    with `--b1-unpol mstw`; the default `toy` backend is outside the window
>    and its numbers are not covered by the lift**, and they cannot be rescaled
>    into covered ones (mstw/toy on `Li6ConvolutionB1::b1(x, 2.5)` is
>    **1.848 / 1.276 / 0.817** at x = 0.10 / 0.30 / 0.50 — a shape change).
>    Selecting `mstw` needs the optional PYTHIA tier and is **refused at
>    configuration time, never downgraded**; a build that cannot reproduce the
>    verdict row therefore cannot emit a number claiming it, and its doctest
>    reports the verdict case as **skipped by name** (`b1_nuclear T1v`, plus
>    `T1r`, which prints in every build whether the row was measured).
> 2. **Comfortable at Eq. (21), marginal at Eq. (17).** MSTW at κ = 1 gives
>    **0.520** — inside by 4 % of its own value. "The gate passes" without "at
>    CDKS Eq. (21)'s δ-function" oversells it.
> 3. **The CD-Bonn agreement is better than the reference deserves.** The
>    target is a *digitization of a published figure* and this gate cannot
>    resolve the low-x zero better than a few per cent. Read 1.000338 as "the
>    residual is now below the digitization error", never as three-digit
>    agreement with CDKS — and it is specific to CD-Bonn **and** MSTW together.
> 4. **The gate is A = 2.** It validates the kernel on the **deuteron** and
>    tests nothing about the α–d step, for which no measurement exists at any
>    A > 2. **The mandatory ±100 % band on every ⁶Li number therefore stays** —
>    it comes from Q(⁶Li) against Q_d, not from the gate — and so does the rule
>    that a number is quoted with the configuration that made it.

### What was built

Per nucleon, ⁶Li (design §1.7):

> b₁^{⁶Li}(x,Q²) = (2/6) ∫ (dz/z) { [f_S(z) + w_CG f_D(z)] b₁ᵈ(x/z,Q²)
>   + w_αd δ_T f_αd(z) F₁ᵈ(x/z,Q²) } + (4/6) w_αd ∫ (dz/z) δ_T f_α(z) F₁^α(x/z,Q²)

with **w_CG = 1/10 exactly** — the Clebsch–Gordan tensor dilution of the
deuteron inside the L = 2 α–d component, which is the closed form quoted next
to `LI6_B1_RANK2_TRANSFER` — and w_αd = 1 nominal. The four terms are
`b1_embedded_s` (1), `b1_alpha_d_dwave_d` (2d), `b1_alpha_d_dwave_alpha` (2α)
and `b1_cg_dwave` (3); they sum to `b1` to 1e-12 (T7, P2). Inputs: α–d
magnitudes from `data/vmc/momenta/li6_ad1.momentum`, the S–D **sign structure**
from `data/vmc/li6_alpha_d/li6.ad`, the deuteron b₁ᵈ from whichever `TensorSF`
camp the caller injects (default: the **raw** digitized CDKS theory-1 column,
`cdks_b1_raw_per_nucleon()`), F₁ from **CDKS Eq. (22)**
`(1 + γ²) F₂/(2x(1+R))` — *not* `NuclearF2::f1a`, which is the massless form
and differs by 1 + γ² = 1.35 at x = 0.5, Q² = 2.5.

### The A = 2 gate, measured

| clause | result |
|---|---|
| G0a–G0g — CG algebra, Eq. (21)'s two coefficients, ∫δ_T f dy = 0, `k_range` against a brute-force scan | **PASS** (exact / 1e-12 / endpoints to the scan step) |
| G1a–G1f — ∫k²(u²+w²)dk = 0.99998, P_D = 0.05760, ∫f dy = 0.98766 against the moment identity's 0.98707, ⟨y⟩ = 0.99526, y_max = 1.4976, the F₁ᴰ/F₁ᴺ EMC shape | **PASS** |
| layer 2 — the α–d sign gate | **PASS**: Q(α–d) = **−0.3333 fm²** < 0, and the same code on the deuteron pair returns `fdeut.av18`'s own header qm = 0.269673 as **0.269362** (0.12 %) |
| **G3a** — exactly two sign changes in (0, 1.0], positions, peak position | **PASS, QUALIFIED.** It passes on every configuration the clause is *applied* to, and must be read with **"G3a's stated limitation — the counting ceiling"** below, which measures what the (0, 1.0] ceiling excludes. One measured configuration is **exempt by construction and is not a pass**: **ToyF2 + CD-Bonn** has **zero** sign changes in (0, 1.0] — `tests/test_b1_nuclear.cpp` checklist item 5 pins it as `CHECK(t.nz == 0)` and deliberately does *not* call `check_g3a` on it, because applying a shape clause there would gate the *wave function* on the *toy* PDF; both real PDFs put the two zeros back. Verdict row (MSTW2008 LO, Eq. (21)): zeros **0.0279** (falling) / **0.4952** (rising) and peak at **0.7716** against the digitized 0.0656 / 0.4572 / 0.7657 — misses of 0.038 / 0.038 / 0.0059 against tolerances 0.08 / 0.10 / 0.10. CDKS's own PDF reproduces CDKS's own figure better than either stand-in on all three landmarks, with nothing tuned. **Read the margin**: the low-x zero is resolved to a few per cent at best on any of these grids (0.0220 ToyF2, 0.0098 CT18NLO, 0.0279 MSTW) and must not be quoted to more than two significant figures; the second zero and the peak position are the discriminators. The counting window is **(0, 1.0]** and the scan floor **0.001**, both stated in the design since 2026-09-03 — the old [0.02, 1.0] was narrower than the ±0.08 tolerance it bracketed, and with CT18NLO the zero fell below it and aborted the clause |
| **G3b** — peak magnitude within a factor 2 | **PASS**: max\|x·b₁\| over [0.10, 0.80] = **9.15096e−4** against the digitized 1.08521e−3, **ratio 0.843243**, a factor **1.19 low**, inside [0.5, 2] with the lower edge cleared by a factor 1.69. With CD-Bonn as well: **1.000338**. Same clause on the other rows: ToyF2 at Eq. (21) 0.440 (outside), ToyF2 at κ = 1 0.272, CT18NLO at Eq. (21) 0.719, MSTW at κ = 1 0.520 |
| G3c — Close–Kumano, **reported, not enforced** | this kernel **+2.24896e−4** (MSTW, Eq. (21)) against the digitized CDKS **+4.59200e−4** — 0.49 of it, essentially where ToyF2 was (0.467). **The PDF that fixes the peak does not fix the integral.** With CD-Bonn as well it is +4.48580e−4 = 0.98, which is a *new* fact and not a restatement of the peak one. Miller's is +5.9147e−3. None is zero; none is enforced |

**G3a's stated limitation — the counting ceiling.** G3a counts sign changes in
**(0, 1.0]**. That upper edge is a **scope choice, not a derived tolerance**,
and it is load-bearing. Over the digitized reference's **own** domain
[0.010, 1.590] the reference has **two** sign changes and stays **positive**
from 0.4572 to its last point (+4.040e−6 at x = 1.590), while **every**
computed configuration crosses zero a **third** time and is **negative** at
the top. So over the reference's own domain the computed curve has **three**
sign changes and the reference has **two**; "exactly two" holds only because
of the ceiling. Measured on the gate's own grid `linspace(0.001, 1.59, 300)`
and pinned by `tests/test_b1_nuclear.cpp` **T17**:

| configuration | third zero (falling) | digitized x·b₁ there | % of the digitized peak |
|---|---|---|---|
| ToyF2, Eq. (21), AV18 — **the shipped default** | 1.220437 | +7.6061e−5 | 7.01 % |
| **MSTW2008 LO, Eq. (21), AV18 — the verdict row** | **1.217660** | +7.7860e−5 | 7.17 % |
| CT18NLO, Eq. (21), AV18 | 1.197722 | +9.1937e−5 | 8.47 % |
| ToyF2, Eq. (17) κ = 1, AV18 | 1.142281 | +1.4599e−4 | 13.45 % |
| MSTW2008 LO, Eq. (17) κ = 1, AV18 | 1.137076 | +1.5245e−4 | 14.05 % |
| CT18NLO, Eq. (17) κ = 1, AV18 | 1.133695 | +1.5680e−4 | 14.45 % |
| ToyF2, Eq. (21), CD-Bonn | 1.500803 | +8.3317e−6 | 0.77 % |
| MSTW2008 LO, Eq. (21), CD-Bonn | 1.493052 | +8.8438e−6 | 0.81 % |
| CT18NLO, Eq. (21), CD-Bonn | 1.473742 | +1.0256e−5 | 0.95 % |

**Widen the window to (0, 1.59] on exactly the argument used for the floor and
G3a fails on every configuration, MSTW + CD-Bonn included.** The floor was
widened in phase A on the argument that *"bookkeeping that brackets a physics
tolerance has to be at least as wide as the tolerance, or it is a second,
tighter, unstated cut."* Applied to the ceiling that argument does **not**
bite the same way — ±0.10 about 0.4572 reaches 0.5572 and the peak clause's
±0.10 about 0.766 reaches 0.866, so (0, 1.0] already brackets every tolerance
this clause writes down, with room. What the floor argument does not reach,
and what the ceiling actually is, is the **count**: a counting clause is a
statement about a domain, and the domain decides the answer. G3a's domain is
chosen, and it is narrower than the reference it compares against. **The phase
that re-opened the window argument moved the edge that cost it nothing and
stopped at the edge that would have cost it the clause.** That is recorded
here rather than left for the next reader to find.

**Two defences of the ceiling, and what survives of them.**

1. *"x > 1 is out of range."* **False, and it must not be used.** x per nucleon
   for a nucleus runs to A — to 2 for the deuteron this gate is about — so
   x ≈ 1.13–1.22 is kinematically allowed and physically meaningful (the Fermi
   motion / short-range-correlation region), and it is **inside** the digitized
   reference's own domain, which CDKS themselves plotted to 1.59.
2. *"the reference is unreadable that high."* **Partly true, and only for the
   CD-Bonn rows.** At the AV18 crossings the digitized x·b₁ is still **7.0–8.5 %**
   of its own peak, and at the κ = 1 crossings **13.5–14.5 %** — far above
   anything one can call a digitization floor, and the digitized column there
   is a smooth monotone decay, not jitter. So for those rows the third crossing
   is a **real sign disagreement in the several-per-cent-of-peak region**, not
   noise. The defence becomes honest only near x ≈ 1.47–1.50 (the CD-Bonn rows),
   where the reference is 0.8–1.0 % of its peak.

**Status: G3a is a QUALIFIED pass, and the clause is not changed here.**
Recorded, not resolved, because changing an acceptance criterion is the
design's and the author's call, not a correction pass's — and because moving
it in either direction after seeing the answer is exactly the failure mode
this note is about. What a reader must take from it:

* The ban lift is triggered by design §5.4's **Escalation** clause, which reads
  the **peak ratio** — that is G3b, and the ceiling does not touch it.
* But G3a is listed as *"shape, **hard**"* in an acceptance criterion whose own
  words are *"three parts, all must hold"*. **If the ceiling were widened to the
  reference's own domain, G3a would fail on every configuration and the
  criterion as a whole would not hold**, so item 10's close would have to be
  re-argued. Anyone quoting "the A = 2 gate passes" is quoting a pass that
  includes this scope.
* Nobody may quote "exactly two sign changes" without the window. The honest
  sentence is **"exactly two sign changes in (0, 1.0], and a third at
  x ≈ 1.22 that the window excludes."**

**The residual is attributed, and after this phase there is essentially none
left.** The §5.4 checklist was worked in order with no fudge factor tuned.
Starting point: the κ = 1 `ToyF2` column, ratio 0.271689. All effects measured
on the doctest grid `linspace(0.001, 1.59, 300)`.

| item | measured effect on the peak | direction |
|---|---|---|
| 0. finite-\|q⃗\| δ-function, κ = 1 → √(1+γ²) | **×1.6195** (ToyF2) / ×1.6206 (MSTW) | closes — **IN**: the gate's default since 2026-09-03 |
| 4. nucleon PDF, `ToyF2` → **MSTW2008 LO**, CDKS's own | **×1.9152** (at κ = 1) / ×1.9165 (at κ) | closes — **the verdict row**, opt-in through `MstwSF` |
| *(the retired CT18NLO stand-in, for reference)* | ×1.6680 at κ = 1; CT18NLO → MSTW is ×1.1721 | superseded |
| 1. target mass (CDKS Eq. 22, already in) | ×1.49 at x = 0.8 (1.52 at κ = 1) | in |
| 2. R (`r1998`, CDKS's own; in **on the gate** — the ⁶Li backend uses `r_sigma_lt`, see below) | ×0.98; swapping it back is +2.7 % on G3b | in |
| 5. wave function, AV18 → **the real CD-Bonn** | **×1.18631** (MSTW, Eq. (21)); ×1.198–1.202 on the other three rows | closes — opt-in through `Options::wave = kCdBonn` |
| *(the retired D-state rescaling proxy for item 5)* | ×0.88 — it had the **sign of the effect backwards** | superseded |
| **the gate's default wave function + MSTW + Eq. (21)** | **ratio 0.843243 (factor 1.19)** | **INSIDE G3b** |
| **CD-Bonn + MSTW + Eq. (21)** | **ratio 1.000338** | — |

Items 0 and 4 are independent to **0.07 %**: 0.271689 × 1.91519 × 1.61951 =
0.842682 predicted against 0.843243 measured, so there is no third unexplained
factor at the peak. The 1.186 that was left over after them is supplied, to
1.3 % across two PDFs and both δ-functions, by item 5 — **the wave function,
not a fudge**. `Li6ConvolutionB1`'s own `finite_q_delta` stays `false` (A10):
in ⁶Li the switch is worth −1 % to +7 %.

**Two things this budget does not say.** The convolution is not a pointwise
multiplication: MSTW/CT18NLO on the isoscalar F₁ is 1.40–1.44 at the peak x
but the peak x·b₁ moved only ×1.1721, and the peak itself moves
0.7344 → 0.7716. And **x = 0.50 got worse**: MSTW gives +1.65e−5 there against
the digitized +1.48e−4, a factor 9 low, because MSTW's second zero sits at
0.4952 and x = 0.50 is 0.005 past it while the digitized curve crossed at
0.4572. No clause of G3 sees that — G3a gates the zero's *position* (cleared by
0.062) and G3b the peak — but a reader comparing curves at x = 0.5 will.

### The ⁶Li numbers — a ToyF2 measurement, NOT covered by the lift until rerun

x·b₁ per nucleon, Q² = 2.5 GeV², default options — which include **ToyF2 F₁
with `r_sigma_lt`** (not the gate's `r1998`; see below) and **κ = 1** (not the
gate's Eq. (21)), and the 2400 / 2400 / 3200 y grid the review of 2026-09-03
set. Every row is the centre of a mandatory {0, 1, 2} × b₁ band:

> **These values must be REGENERATED under `--b1-unpol mstw` before they are
> quoted as physics.** The 2026-09-03 lift is a statement about a
> configuration: the MSTW2008 LO nucleon input at CDKS Eq. (21). The table
> below is the **`ToyF2`** configuration, whose own G3b is **0.440**, outside
> the [0.5, 2] window the lift was read off — and the difference is a shape
> change, not a normalisation (mstw/toy on `Li6ConvolutionB1::b1(x, 2.5)` is
> **1.848 / 1.276 / 0.817** at x = 0.10 / 0.30 / 0.50; ×1.92 on the gate's own
> peak), so no rescaling converts one into the other. What survives untouched
> is only what is algebraic: the four terms summing to `b1` to 1e−12 and the
> band's exact linearity. Every **value** below, and every ratio between terms
> and derived systematic read off it, is a ToyF2 measurement until it is rerun.
> `docs/USAGE.md` §2a carries the same decision and the same reasons.

*(The `Li6B1(CdksB1)` column DOUBLED on 2026-09-03 — `b1_convolution()` stopped
halving a column that is already per nucleon; see Q1b below. Every other column
is unchanged, because `Li6ConvolutionB1` already reached the raw column through
`cdks_b1_raw_per_nucleon()`.)*

| x | (1) embedded d, S | (3) CG D-wave | (2d) struck d | (2α) struck α | **total** | `Li6B1(MillerB1)` | `Li6B1(CdksB1)` |
|---|---|---|---|---|---|---|---|
| 0.05 | +1.5686e−6 | +3.0979e−9 | +1.3838e−6 | +6.9528e−7 | **+3.6507e−6** | +3.6503e−4 | +1.7983e−6 |
| 0.10 | −4.6926e−6 | −9.2965e−9 | +1.2673e−6 | +6.3805e−7 | **−2.7966e−6** | +4.1222e−4 | −5.3132e−6 |
| 0.20 | −2.3664e−5 | −4.6736e−8 | +1.6669e−6 | +8.3886e−7 | **−2.1205e−5** | +2.6700e−4 | −2.6982e−5 |
| 0.30 | −4.4967e−5 | −8.8342e−8 | +2.2133e−6 | +1.1092e−6 | **−4.1732e−5** | +6.7167e−6 | −5.1928e−5 |
| 0.50 | +4.2448e−5 | +8.5642e−8 | +6.1632e−6 | +3.0422e−6 | **+5.1739e−5** | −1.8650e−4 | +4.5585e−5 |

Reading it:

* **Term (3) is 0.197 % of term (1)** at every x — exactly the CG algebra's
  0.1 × P_D^{αd}. It is kept because it is free and because it is the analytic
  bridge to `LI6_B1_RANK2_TRANSFER`, not because it matters numerically.
* **Terms (2d)+(2α) are 7–133 % of term (1), with the OPPOSITE sign** over the
  whole window, because the α–d S–D relative sign is opposite to the
  deuteron's below the S node (66 % of the density). They nearly cancel term (1)
  at x ≈ 0.1 and dominate below. **This is the ⁶Li quadrupole puzzle showing up
  in b₁ and it is the quantitative reason the 100 % band is mandatory.**
* **The struck-α term is half the orbital sector**: (2α)/(2d) = 0.494–0.503
  against the analytic 2(M_d/M_α)² = 0.5064 (T15). Dropping it is a 30 % error.
  *(On the design's 600/2400/800 y grid this column read 0.489–0.497 and term
  (2d) itself was 2.3 % high at x = 0.10 — a quadrature artefact of the struck
  deuteron's wider z-distribution, fixed by the 2400/2400/3200 default and now
  guarded by a ⁶Li clause in T2.)*
* Against the production default `Li6B1(MillerB1)` the new model differs by two
  orders of magnitude **and by sign**, because Miller (HERMES-like) and CDKS
  (convolution) are different *camps* for b₁ᵈ — a pre-existing disagreement this
  work does not resolve. **The default does not change.**

Densities (struck d): raw `norm()` **0.8176766** against the identity
N_αd⟨E⟩/M_d = 0.817421 (3.1e−4 apart, T3's own tolerance); `mean_y()` =
⟨z²⟩/⟨z⟩ **0.9997871** (⟨z⟩ ≈ 1, *not* 1/3 — the fair-share normalisation of
design §1.6 is what makes this literally CDKS Eq. (16) one level up); `p_d()` 0.0193299, `p_d_momentum()` 0.0193516
against `VMC_P_D_LI6` = 0.0193549 (1.7e−4, the file's own quadrature spread).
∫b₁ dx over [0.01, 1.2]: **+1.4588e−4**, against `Li6B1(MillerB1)`'s +9.1243e−4
and `Li6B1(CdksB1)`'s +1.38596e−4 (that one doubled on 2026-09-03; it was
+6.9298e−5). `finite_q_delta` moves ⁶Li's x·b₁ by only
−1 % to +7 % (against 1.57–1.69 on the A = 2 gate), because in ⁶Li the term κ
multiplies is the small orbital one — which is why the two objects default
differently.

**The ⁶Li backend and the gate use different R — decided 2026-09-03, and both
defaults stay.** `Li6ConvolutionOptions::r_func` null ⇒ **`r_sigma_lt`** (the
kernel's own toy R); `DeuteronConvolutionB1::Options::r_func` null ⇒ **`r1998`**
(CDKS's SLAC world fit, because the gate reproduces *their* figure). So the R
the gate was validated with is not the R the shipped ⁶Li numbers above carry.

**What decides it is that the shipped observable is a RATIO.** The tensor
weight is K/D_φ with K ∋ b₁ and D_φ ∋ F₁, and that F₁ is `InclusiveKernel`'s,
built from the `ToyF2` the pipeline *shares* with `Li6ConvolutionOptions::unpol`
— whose R is `r_sigma_lt`, because nothing sets it. Give the b₁ backend a
different R and the (1 + R) in the numerator stops cancelling the one in the
denominator: (1 + r1998)/(1 + r_sigma_lt) at Q² = 2.5 is 1.088 at x = 0.1 and
1.039 at x = 0.3, the same size as the effect being chased. Fidelity to CDKS is
worth having on the *gate*, which has no denominator; inside the pipeline
consistency wins. **This is an author decision, not a derivation** — the third
option, and the better one, is to thread ONE R hook through
`default_inclusive_kernel` into both objects so they cannot disagree by
accident. That is a wiring change and it is not made here.

**Measured cost of the choice** (Q² = 2.5, default options), x·b₁:

| x | `r_sigma_lt` (the default) | `r1998` | change |
|---|---|---|---|
| 0.05 | +3.650716e−6 | +3.449387e−6 | −5.5 % |
| 0.10 | −2.796612e−6 | −3.484857e−6 | +24.6 % |
| 0.20 | −2.120522e−5 | −2.235704e−5 | +5.4 % |
| 0.30 | −4.173248e−5 | −4.284379e−5 | +2.7 % |
| 0.50 | +5.173934e−5 | +5.075777e−5 | −1.9 % |

**And it is not the (1 + R) prefactor that does it.** Term (1), the embedded
deuteron, is **bit-identical** under the swap — it carries b₁ᵈ from the injected
`TensorSF`, the raw digitized column, which has no R in it at all. The whole
effect sits in the two orbital terms, whose F₁ᵈ slot *is* `f1_cdks`, and there
it is **×0.54 to ×0.94**, far more than the ≤ 8 % a prefactor allows: those
terms convolve F₁ᵈ(x/z) against a density that **integrates to zero**, so they
respond to the *slope* of R — and `r_sigma_lt` is x-independent by construction
while `r1998` runs from 0.30 at x = 0.05 to 0.20 at x = 0.5. The +24.6 % at
x = 0.10 is that ×0.64 seen through the near-cancellation between term (1) and
terms (2d)+(2α). On the **gate**, where there is no such cancellation, the same
swap is worth **+2.7 %** on G3b (0.8432 → 0.8662, MSTW at Eq. (21), max over a
0.01 grid on [0.10, 0.80]) and, at κ = 1, moves the second zero from 0.392 to
0.365 — *away* from the digitized 0.457. **The gate's verdict does not rest on
this choice either way.**

**The truncation is MEASURED, not asserted.** The dropped P₂ and P₄ remainders
of the D-wave b₁ᵈ weight are **below 6e−5 of term (1) everywhere** (T14) — four
orders of magnitude inside the band. The reason is *not* "P_D × 1/10" (the P₄
coefficient is 15× the isotropic one term (3) keeps): they are small because the
P₂- and P₄-weighted densities integrate to zero, so they enter only through the
*curvature* of b₁ᵈ(x/z), and because they multiply b₁ᵈ, not F₁ᵈ. The other
dropped piece — the S–D interference in the b₁ᵈ sector — is term (2d)'s SD
structure with coefficient exactly 1/3 of it (G0f, exact) and b₁ᵈ in place of
F₁ᵈ: ~3e−4 of term (2d).

### Systematics that are stated, not hidden

* **The mandatory 100 % band.** Q(⁶Li) = −0.0818(17) fm² against
  Q_d = +0.2859(3) fm²: the α–d D wave enters the closest measured observable
  with the opposite sign to the deuteron's own D state and nearly cancels it.
  *(Q(⁶Li) is `LI6_QUADRUPOLE_FM2` in `include/lipolgen/rc.hpp`, TUNL's A = 6
  evaluation, 1998CE04; Pyykkö's compilation, Mol. Phys. **106** (2008) 1965,
  gives −0.0806(6) fm² from Cederberg et al., Phys. Rev. A **57** (1998) 2539
  — the repository constant is TUNL's. Q_d: Bishop & Cheung, Phys. Rev. A
  **20** (1979) 381.)* Every number is `--b1-band-scale 0/1/2`,
  and `--b1-alpha-d-dwave-weight 0/1/2` is the shape variant reported next to
  it.
* **±5 % on N_αd, unexplained.** Three tabulations of the α–d spectroscopic
  factor span 5 %: **0.819481** (2014 `li6_ad1.momentum`, the default, and now
  the single home `VMC_N_ALPHA_D_LI6` that `VMC_P_D_LI6` is expressed through),
  0.856 (2004 `li6.ad`) and 0.863 (Wiringa et al., PRC **89** (2014) 024305
  §III: 0.846 + 0.017). Entries 1 and 3 are the same year and the same
  Hamiltonian family, so this is **not** a version difference anyone can name.
  b₁ is exactly linear in `norm_target`, so it is a flat ±5 % (T12).
* **Q4 — suppression or renormalisation?** The default gives the 18 % of ⁶Li
  that is not α+d **b₁ = 0**; `use_spectroscopic_factor = false` renormalises
  instead and is ×1/N_αd = **1.2203** on the whole answer. Both are defensible,
  nobody has published either, and the difference is well inside the band.
* **Q1b — RESOLVED on the CDKS half, an author decision on the Miller half
  (2026-09-03).** The two camps no longer share one applied constant.
  * **The arbiter is HERMES, not either theory paper.** Their published b₁ᵈ —
    the data both camps plot against — is **per nucleon**, because their Eq. (5)
    divides by an F₁ᵈ built from F₂ᵈ = (F₂ᵖ + F₂ⁿ)/2. Inverting their own
    Table II reproduces that F₁ᵈ in all six bins (mean ratio **0.946**, against
    **0.473** for the per-deuteron one) across five decades of b₁ᵈ and a factor
    9 in Q². *(Their own QPM definition table, read literally, is the
    per-deuteron object — and that inconsistency is exactly the trap that
    caught the two theory camps in opposite directions.)*
  * **CDKS: certain, and `CdksB1` stopped halving.**
    `B1_CDKS_TABLE_TO_PER_NUCLEON` = **1**. Eq. (10)'s spectral function carries
    an explicit 1/A, the text under Eq. (16) says in words *"the structure
    function b₁ is defined by the one per nucleon"*, f(y) is normalised to one
    nucleon and F₁ᴺ = (F₁ᵖ + F₁ⁿ)/2. **`b1_convolution` / `CdksB1` /
    `--b1-model cdks` therefore DOUBLED.** Nothing on the default path moved:
    the only b₁ reference JSON records `"b1_model": "miller"`, so the rtol-1e−12
    gate is untouched, and `Li6ConvolutionB1` already reached the raw column
    through `cdks_b1_raw_per_nucleon()`.
  * **Miller: likely, not certain — the 0.5 stays, as an author decision.**
    `B1_MILLER_TABLE_TO_PER_NUCLEON` = **0.5**, and it *is*
    `B1_PER_DEUTERON_TO_PER_NUCLEON` rather than a second copy of the number.
    His Eq. (1) densities are "in a target hadron", his Eq. (5) is a light-cone
    correlator in the normalised deuteron state, and Eq. (6)'s ½ is fully
    consumed by the quark-spin average of a *spinless* pion — every ½ in
    Eqs. (1)/(5)/(6)/(20) is spoken for and none of them is a 1/A. **Against
    that**, his Table I transcribes HERMES's per-nucleon numbers unrescaled,
    his Fig. 5 overlays them on the curve, and he tunes P₆q to one of them at
    face value. The paper is self-inconsistent by exactly this factor, and no
    third party adjudicates. Keeping the 0.5 is the status quo, moves no
    default and touches no gate; **the numerical consequence of the other
    reading is that every `MillerB1` number would double** — including the
    default inclusive tensor rate and its rtol-1e−12 reference. What would
    settle it: ask Miller, or evaluate his Eq. (20) at x = 0.012 and see
    whether it lands on 10.5 × 10⁻² (per nucleon, the 0.5 is wrong) or
    5.25 × 10⁻² (per deuteron, the 0.5 is right) — which needs the pion PDF set
    of his ref [29], model 1, that this tree does not have.
  * **A retraction while here.** This document previously said the answer to
    Miller's normalisation "is on the axis of Miller's Fig. 5". It is not: that
    ordinate is `100 b₁(x)` and nothing else, his caption is bare, and the
    strings "per nucleon" and "1/A" do not occur anywhere in the paper.
  * Full argument, with the six-bin inversion and its reproduction script:
    `docs/open_items/run_2026-09-03/phase_A_miller_normalisation.md`. One
    consequence to carry: `close_kumano_integral(false)` and
    `close_kumano_integral(true)` integrate the RAW columns and are therefore on
    **different scales** — never compare them to each other.
* **Q6.** `LI6_B1_RANK2_TRANSFER` = 0.921947 is tied to the Hulthén-scenario
  P_D = 0.0867; the VMC value 0.0193549 would give 0.982581. Changing it would
  move every published inclusive tensor number, so it stays — but this backend
  makes the inconsistency visible (the tagged channel already uses VMC at
  `--cluster-wave vmc`; the inclusive b₁ never does).
* **Q8 — there is nothing to validate against for A = 6.** No published b₁
  exists for any A > 2. The A = 2 gate plus the analytic limits (T5, T6, T7) are
  the *only* checks there are, which is why the gate is blocking and why the
  band is not a formality. **The gate closing on 2026-09-03 does not change
  this**: it says the *kernel* reproduces CDKS on the deuteron, and says nothing
  about the α–d step it is applied to one level up.

### The seven close conditions, and where each one landed

All seven were worked on 2026-09-03. Five are **done**, two are recorded as
**author decisions** — a decision, with its evidence and its measured cost, is
a close and not a deferral.

| # | condition | outcome |
|---|---|---|
| 1 | Decide item 0's default (the δ-function) | **DONE.** The A = 2 gate is quoted at CDKS Eq. (21), their exact definition; the ⁶Li backend keeps Eq. (17)'s κ = 1, where the switch is worth 1 %. It turned out to be load-bearing for the verdict, not only for the third digit: MSTW gives 0.843 at Eq. (21) against 0.520 at Eq. (17) |
| 2 | Get MSTW2008 LO (or Kumano's tabulated curves) and rerun checklist item 4 with CDKS's own PDF | **DONE.** MSTW2008 LO was on disk all along — the CENTRAL member ships with PYTHIA 8 as `pdfdata/mstw2008lo.00.dat`; only LHAPDF's set store lacked it. `MstwSF` (`include/lipolgen/mstw_sf.hpp`) reads it with the same charge weights as `LhapdfSF`, so an MSTW/CT18NLO ratio is a PDF comparison and not a convention one. **G3b 0.719 (CT18NLO stand-in) → 0.843243.** Kumano's curves were not needed |
| 3 | Get a real CD-Bonn u, w instead of the D-state rescaling proxy | **DONE.** Machleidt's Appendix-D parameterisation is `cdbonn_wave()` in `cluster.hpp`, with C₁₁ and D₉..D₁₁ **computed** from the r → 0 boundary conditions rather than typed; reached through `DeuteronConvolutionB1::Options::wave = kCdBonn`. It reproduces his Table XV (P_D 4.8562 % against 4.85, η 0.0255714 against 0.0256(4), Q_d +0.270178 fm² against 0.270) and is **×1.18631** on the G3b peak — it **closes**, which is the opposite sign to the retired proxy's ×0.88. The default stays AV18 |
| 4 | Decide whether `Li6ConvolutionOptions` should default to `r1998` like the gate | **DECIDED — it stays `r_sigma_lt`**, because the shipped observable is a ratio whose denominator carries `InclusiveKernel`'s R. Measured cost of the choice and the term-by-term reason are above ("The ⁶Li backend and the gate use different R"). The better third option — one R hook threaded through both — is named there and not taken here |
| 5 | Decide whether `CdksB1` should stop halving the per-nucleon column (Q5(a)) | **DECIDED — yes, and done.** `B1_CDKS_TABLE_TO_PER_NUCLEON` = 1; `--b1-model cdks` doubled. Certain, from CDKS's own words plus the HERMES inversion (Q1b above). Two pinned values in `tests/test_sf.cpp` doubled with it; no reference JSON moved |
| 6 | Check **Miller's** own normalisation independently (Q5(b)) | **DECIDED AS AN AUTHOR DECISION — the 0.5 stays, and it is *likely*, not certain.** His derivation is per deuteron and his presentation is per nucleon; the paper cannot be both. Evidence for each reading and the numerical consequence of the other (every `MillerB1` number, and the default inclusive tensor rate, would double) are in Q1b above. **The framing this document gave the question was wrong and is retracted**: the answer is not on the axis of his Fig. 5 |
| 7 | Only then re-open G3b | **RE-OPENED, AND IT PASSES.** Ratio **0.843243** with MSTW at Eq. (21); **1.000338** with CD-Bonn as well. Design §5.4's *Escalation* clause is not triggered, so **the ⁶Li publication ban is lifted** — read the four conditions in the warning block at the top of this section before quoting anything |

**What is still open, stated as what it is.**

* **The Miller normalisation (condition 6)** is a decision, not a measurement.
  Closing it properly needs the author or a numerical evaluation of his
  Eq. (20), which needs a pion PDF set this tree does not have.
* **The R wiring (condition 4)** is a decision too; the clean fix is one R hook
  threaded through `default_inclusive_kernel` into both objects.
* **G3c is not enforced and did not improve with the PDF** (0.467 → 0.490 of
  the digitized integral; 0.98 only with CD-Bonn as well). The peak agreeing
  does not make the integral agree.
* **x = 0.50 is a factor 9 low** on the verdict row, for the zero-position
  reason recorded above. No clause of G3 sees it.
* **The ±5 % N_αd spread is still unexplained** (0.819481 / 0.856 / 0.863 from
  three tabulations, two of them the same year and Hamiltonian family).
* **Q6, Q4 and Q8 are unchanged** by this phase: `LI6_B1_RANK2_TRANSFER` is
  still tied to the Hulthén-scenario P_D, `use_spectroscopic_factor` is still a
  1.22 fork, and there is still nothing to validate an A = 6 b₁ against.
* **The default wave function is still AV18**, deliberately: making `kCdBonn`
  the default of `DeuteronConvolutionB1` is a strong argument — that object
  exists to reproduce CDKS's figure and CDKS used CD-Bonn — but it would move
  every pinned gate number, so it belongs to a task that says so.

## 11. Coherent ⁶Li amplitude

Sartre: nuclei hard-coded (`Nucleus.cpp` switch, no A = 6), spherical sampling
(no polarization axis), tables CPU-years (the 2026 speed-up code is unreleased) —
ruled out. Outdated project claims: 2605.00454 publishes coherent J/ψ down to
A = 3, 4 with α-clustering; 2511.05638 has e+⁷Li coherent J/ψ |t| distributions
with tagging efficiency. Path, both halves now **in-tree**: (11.1) eSTARlight
⁶Li for the unpolarized rate and slope baseline; (11.2) an α+d configuration
sampler (α core from the VMC density, p–n pair from AV18 u/w oriented by the
polarization axis, α–d separation from the VMC α–d overlap) whose output grafts
into `hejajama/subnucleondiffraction` (author asks to be contacted) for the
tensor cos 2φ. Warn in plans/06: photon-polarization cos 2φ (STAR 2204.01625)
is a distinct mechanism and a background at Q² > 0.

### 11.1 eSTARlight — the unpolarized rate and slope baseline

**In-tree since 2026-09-02.** eSTARlight (`github.com/eic/estarlight`, commit
`939b11a24499398392d959db81c7502aeec91046`, same code arXiv:2511.05638 cites)
runs e+⁶Li today with no code change — `nucleus::init()` has no Z = 3 case, so
it falls to the generic light-nucleus branch (`_Radius = 1.2·A^{1/3}`, Gaussian
form factor). Full run log, input files and an adversarial re-check are
`docs/open_items/run_2026-09-02/estarlight_li6.md`; nothing there was
committed. 2×10⁵ events per channel, `0.1 < Q² < 100 GeV²` (arXiv:2511.05638's
own window), e 10 GeV × ⁶Li 99.5 GeV/u:

| channel | σ_coh (default R = 2.1805 fm) | B fitted [GeV⁻²] | σ_coh (measured R = 2.589 fm) | B fitted |
|---|---|---|---|---|
| coherent ρ⁰ | **506.4 nb** | 38.6 | 379.5 nb | 54.1 |
| coherent φ | **30.16 nb** | 39.2 | 22.08 nb | 54.8 |
| coherent J/ψ | **1.773 nb** | 38.9 | 1.255 nb | 55.0 |

The fitted slope is VM-independent to 1.5 % (there is no VM-dependent slope in
a coherent Gaussian-form-factor calculation) and reproduces the analytic
B = R_G²/(3ħc²) to 4 %. The ⁷Li cross-check against arXiv:2511.05638 (identical
configuration, same commit-era code, same Q² window, same beam) agrees with
every *shape* statement the paper makes — no diffractive minimum, the |t| slope
and the Fourier-transform width read off its Figs. 6–7 — but the paper quotes
no absolute σ, so no number-to-number check is possible.

**Conclusion for `CoherentScenario::slope_b`.** `gaussian_slope(r_rms)` in
`coherent.cpp` (B = R_rms²/(3ħ²c²)) is algebraically the same object as
eSTARlight's light-nucleus form factor; only the radius differs. The scenario
band `slope_b ∈ {40, 60}` GeV⁻² corresponds to R_rms ∈ {2.162, 2.647} fm, which
brackets almost exactly the two defensible ⁶Li densities — eSTARlight's own
`1.2·A^{1/3}` (fitted **38.6–39.2**) at the bottom and the Angeli–Marinova
measured charge radius 2.589 fm (fitted **54.1–55.0**) at the top, i.e. a
citable **B = 39–55 GeV⁻²** band from simulation rather than a hand-tuned one.
**`slope_b = 50 ± 10` survives contact with eSTARlight unchanged, and the band
could be tightened to 45–58 GeV⁻² but there is no reason to** — the ±30 % rate
swing between the two radii shows the density assumption, not the dynamics, is
the leading systematic.

**Conclusion for `CoherentScenario::f0`.** eSTARlight generates *exclusive*
vector-meson production only — M_X is the meson mass exactly, with no
diffractive continuum — so **it cannot measure f0**, LiPolGen's coherent
fraction of the DIS rate at x → 0. What it *can* bound: running ⁶Li at
`MIN_GAMMA_Q2 = 0.7` against `generate_inclusive`'s σ_incl = 592 nb in the same
window gives σ(ρ+φ+J/ψ)/σ_incl = 3.0 × 10⁻² (all three channels; f0 = 0.04 is
not refuted at the order-of-magnitude level, but 90 % of this is ρ⁰, which sits
below LiPolGen's own M_X ≥ 1.2 GeV coherent floor) and
σ(J/ψ only)/σ_incl = **1.0 × 10⁻³** — the one channel actually inside
LiPolGen's window, and therefore a genuine **lower bound**, 40× below f0 = 0.04,
on the part of coherent diffraction eSTARlight's exclusive-VM model can see.
Determining f0 itself needs a coherent diffractive-DIS calculation (a
coherent-A analogue of the H1/ZEUS diffractive PDFs) that exists in neither
eSTARlight nor Sartre; **recommendation: leave `f0 = 0.04` and its {0.02, 0.08}
band as a scenario**, and do not treat eSTARlight's numbers as a
determination. **Done 2026-09-03**: the comment on `CoherentScenario::f0`
(`include/lipolgen/coherent.hpp`) now records a **lower bound of 1.0 × 10⁻³**
(exclusive J/ψ inside the M_X ≥ 1.2 GeV window). The all-VM figure 3.0 × 10⁻²
is *not* an upper bound on f0 — ~90 % of it is ρ⁰, which with φ sits below that
floor.

Caveats carried forward unchanged from `estarlight_li6.md` §6: one spherically
symmetric Gaussian density (no α+d clustering, no polarization axis, no
diffractive minimum — not an imaging baseline), and no saturation.

### 11.2 The α+d configuration sampler — the tensor cos 2Φ half

**In-tree since 2026-09-03**: `include/lipolgen/cluster_config.hpp` +
`src/core/cluster_config.cpp`, `lipolgen-configs`,
`tests/test_cluster_config.cpp`, `python/tests/test_cluster_config.py`, design
`docs/open_items/run_2026-09-02/design_G_cluster_config.md`, numbers
`.../phase_G_numbers.md`. **Opt-in and inert**: nothing in the generator calls
it, `CoherentScenario` is untouched, and its output is a file. It writes
`he3.dat`-compatible position tables (fm, ion rest frame, c.m. at the origin,
the polarization axis applied) plus a sidecar carrying the axis, the substate,
every option and every input table's md5. `docs/USAGE.md` §9 has the CLI,
the Python API and the output format.

**What it is.** ⁶Li(1⁺) as a rigid α(0⁺) core plus a deuteron with relative
L = 0, 2 coupled to S = 1: the α's four nucleons come from the ANL VMC ⁴He
one-body density (recentred so the configuration's own c.m. sits at the
origin, exactly), the α–d separation R and its orientation from the ANL VMC
α–d overlap (S+D wave, `li6.ad`/`li6.adr.fit`), and the p–n pair from AV18
u(r)/w(r) — all three sharing one drawn deuteron projection m_S per
configuration, so the R̂-to-r̂ correlation the design calls out is kept and
only the *relative azimuth* of the two is dropped (`tagged.hpp`'s own
truncation, in r space).

**The three m-state densities.** One (R, cos θ_R) table per (m, m_S) pair —
nine tables total, the same discipline `TaggedModel::build_amp2` uses in
k-space — because sampling R̂ from the m_S-*summed* density instead would
decorrelate it from the deuteron spin a moment test cannot see: for m = +1
the exact ⟨P₂(cos θ_R)⟩ is −2/7 at m_S = −1 and +1/7 at m_S = 0 —
reproduced by the sampler's grid to 1e-6 for m_S = −1 and to 2.3e-4 at the
default 96 cos θ cells for m_S = 0 (1e-6 at n_c = 3072: a midpoint-rule
residue at the simple zero of |Θ₂¹|², `phase_G_numbers.md` §3) — while an
m_S-marginal draw would give the m_S-**summed** ⟨P₂⟩ = **−0.0370** in
*every* branch (the design's own earlier estimate of it is −0.0355), which
the two conditioned branches sit 46σ and 17σ away from in a
2×10⁵-configuration test.
The three substates m = +1, 0, −1 (plus the interleaved unpolarized mix) all
read the *same* radial tables; only the Clebsch–Gordan recoupling changes.

**The quadrupole puzzle, and how the sampler handles it rather than hides
it.** This geometry's own point-matter quadrupole overshoots the measured
one by a factor ≈ 7.5: Q_charge(⁶Li) = **−0.615 fm²** (model range
−0.615…−0.730 across the three α–d sources) against the measured
**−0.0818 fm²** (`LI6_QUADRUPOLE_FM2`) and GFMC AV18+IL7's **−0.20(6) fm²**
(Pastore *et al.*, PRC 87, 035503 (2013)) — an independent-Hamiltonian
comparison point, not a check on the same tables. The asymptotic α–d D/S
ratio, correctly divided by the Coulomb Whittaker ratio (which is 3.3 at
R = 6 fm and 2.6 at R = 8 fm — the naive R₂/R₀ is *not* the asymptotic ratio),
gives **η = −0.048** against the measured **η = −0.025 ± 0.006 ± 0.010**
(George & Knutson, PRC 59, 598 (1999)): a real but *moderate* ≈ 2× D-wave
excess, not the 5–15× a naive ratio suggests — which demotes an ANL
normalization/phase-convention error as the leading explanation (design O1).
The sampler does not paper over the gap: `quadrupole_band_fm2()` returns all
three numbers together and the writer stamps them on every output, and
**the rule stands — do not derive a published tensor input from this
geometry alone.** A `quadrupole_target_fm2` **deformation dial** (not a
wave-function fit) can rescale the α–d D-wave to land the geometry's own
Q_charge on any target in the reachable range [−0.615, +0.270] fm²; the root
for the measured Q(⁶Li) is s = 0.4032, which drags P_D(α–d) down to
3.3 × 10⁻³ from 0.02011 — the dial trades away the natural D-state
probability to match the tensor moment, and the sidecar always records which
choice was made.

**The measured numbers** (default `FitRescaled` source; full table, all
three sources and the grid-vs-analytic quadrature study in
`phase_G_numbers.md`): ⟨r²⟩(⁶Li) = **6.4447 fm²** (r_rms 2.5386 fm — **3 %**
above the measured point radius 2.4655 fm, `LI6_R2_POINT_FM2`, and **4 %**
above the VMC `li6.density` value 2.4433 fm, `LI6_R_POINT_VMC_FM`); the
point-matter quadrupole per
substate, Q_matter(m) = (3m² − 2)·[(4/3) Q[R₀,R₂] + 2 Q_d D_T] (eq. (G5)),
is **Q_matter(±1) = −1.2309 fm²** and **Q_matter(0) = +2.4618 fm²** — the two
non-zero substates carry opposite sign by construction, and Q_matter(+1)/2 is
exactly the Q_charge(+1) quoted above. A 2×10⁵-configuration Monte Carlo
closure reproduces both to within 5σ of the sampler's own per-configuration
variance (sampled ⟨r²⟩ = 6.4483 vs 6.4447, 5σ = 0.041; sampled
Q_matter(+1) = −1.2216 vs −1.2309, 5σ = 0.322), and Σᵢ r⃗ᵢ = 0 to 2.2 × 10⁻¹⁵ fm
worst component over 20 000 configurations. Feeding Q_matter(+1) through the
closed-form (G7)+(G8) quadrupole → a₂ map — which reproduces the only
published polarized coherent calculation, Mäntysaari *et al.*'s digitized
deuteron a₂, to 8 % at m = ±1 and 8–21 % at m = 0 with **zero free
parameters** — gives a₂(±1) = **+0.1976** at |t| = 0.3 GeV² for this
geometry's own (overshot) Q, and **+0.026** for the *measured* Q(⁶Li) — about
10× smaller than the deuteron's and of the **opposite sign**, because
Q(⁶Li) < 0 while Q_d > 0. **Timing**: 0.54 µs/config in-process, 1.31 µs/config
through the `lipolgen-configs` CLI, both well inside the design's < 5 µs/config
target; a 10⁵-configuration table with its sidecar writes in ≈ 1 s.

**How the graft would consume the table.** The writer emits one line per
configuration in `he3.dat`'s own layout, in fm, ion rest frame, c.m. at the
origin, with the polarization axis already applied — so the consumer needs no
knowledge of our conventions beyond "these are the six nucleon positions."
Their side then needs the `-configfile/-configid` generalization described
below. Their Good–Walker loop is unchanged; the amplitude uses only x and y.
One set per m ∈ {+1, 0, −1} plus one unpolarized set; a₂ comes out of the Φ
dependence of ⟨A⟩ exactly as in their Fig. 2.

**The exact ask to send the Mäntysaari group** (`design_G_cluster_config.md`
§10, verbatim):

> We have built an α+d configuration sampler for polarized ⁶Li that emits
> nucleon-position tables in exactly the format `Nucleons::InitializeTarget`
> reads for ³He (positions in fm, c.m. at the origin, one configuration per
> line), with the polarization axis and the substate m = +1, 0, −1 applied and
> recorded. The α core comes from your ANL VMC/GFMC ⁴He one-body density, the
> p–n pair from AV18 u(r), w(r) with the full m_S angular correlation (your
> Eqs. (7)–(8), which we reproduce exactly), and the α–d separation from the
> ANL VMC α–d overlap with the L = 2 orientation correlated to m through the
> CG recoupling. Three concrete requests, in increasing size. **(a)** Would you
> share the ⁴He GFMC nucleon configurations used in arXiv:2605.00454? Our α
> core is currently an uncorrelated product of one-body densities, which is the
> weakest part of the sampler and the one part you already have solved.
> **(b)** Would you accept a ~30-line generalization of the `A == 3` branch to
> an `-configfile/-configid` option for any A? We can send the patch; it makes
> the code A-agnostic and removes the 13698-configuration bound. **(c)** The
> polarized-deuteron machinery of arXiv:2408.13213 is not in any public branch
> of the repository — would you run your existing polarized setup on our ⁶Li
> configuration tables, or release that branch? Two honest caveats we would
> want your view on before you spend time: our α+d model reproduces the ⁶Li
> point radius to 3 % but overshoots Q(⁶Li) by a factor ≈ 7.5 (measured
> −0.0818 fm², our α+d geometry −0.615…−0.730 fm², **GFMC AV18+IL7 −0.20(6) fm²**,
> Pastore *et al.*, PRC 87, 035503 (2013)), so the tensor amplitude carries a
> factor-of-several systematic that we would carry explicitly as a band. The
> α–d asymptotic D/S ratio of the overlap we use is η ≈ −0.05 against the
> measured −0.025(12), so the excess is a real but moderate D-wave effect rather
> than a convention error. And a closed-form quadrupole → a₂ map that reproduces
> your published deuteron a₂(m = ±1) to 8 % at every |t| predicts that ⁶Li's
> a₂ is ~10× smaller than the deuteron's **and of the opposite sign**, because
> Q(⁶Li) < 0. That sign flip is the interesting measurement and it is also the
> reason the effect is small — we would rather establish that together, before
> either side commits person-months.

**What remains open.** The tensor cos 2Φ prediction above (a₂(±1) ≈ +0.026 at
the measured Q, opposite sign and ~10× smaller than the deuteron's) is a
closed-form estimate from the geometry's own point-matter quadrupole — **it is
not yet a coherent-diffraction amplitude from a dipole-model run**, and it
cannot become one without the actual `subnucleondiffraction` graft (ask (b)
and (c) above): the Good–Walker amplitude, its Φ-averaging and its own
statistical and saturation-model uncertainties are not reproduced by
`a2_from_quadrupole`, which only carries the target's quadrupole moment
through the deuteron's own published |t| dependence. Also still open, in the
design's own numbering: **O1** how a moderate (≈2×) D-wave excess and a
missing ≈15 % non-α+d component together produce the observed 7.5× gap in Q;
**O2** whether the smoothed `li6.adr.fit` R₂ node near 1.1 fm is a real
short-range effect or a fit artefact (decides the honest default between
`FitRescaled` and `OverlapRaw`); **O3** how much an uncorrelated α core (vs.
GFMC ⁴He configurations with correlations) moves the incoherent/coherent
split; **O4** the `eps_b0` convention (`coherent.hpp` uses ΔB without
defining it — a separate, reviewed decision, not taken here); **O5** whether
the ⁶Li a₂ survives EIC statistics at all and separates from the
linearly-polarized-photon cos 2φ background (`physics_literature.md`'s
two-mechanism warning) — answerable today from the numbers above plus a rate
estimate, and the question that decides whether the collaboration is worth
proposing.

**The rule of this section is unchanged and is now quantified rather than
repealed: do not derive a tensor input for a published observable from these
wave functions.** The α+d truncation reproduces the ⁶Li point radius to 3 %
(4 % against the VMC `li6.density` value) but
overshoots Q(⁶Li) by **7.5×** (model −0.615…−0.730 fm² against the measured
−0.0818 and GFMC AV18+IL7's −0.20(6), Pastore *et al.*, PRC 87, 035503
(2013)). `quadrupole_band_fm2()` returns all three and the writer stamps them;
the `quadrupole_target_fm2` dial can put the geometry on the measured Q, but
it is a **deformation dial, not a wave function**, and the sidecar labels it
as one.

## 12–13. Engineering

- **In-tree now** (2026-09-02; the prototype in `engineering.md` §B was
  redone for real, prototype files themselves are gone): `pyproject.toml`
  (`[build-system] requires = ["scikit-build-core>=0.9", "pybind11>=2.10"]` —
  CMakeLists.txt's own pybind11 discovery is `find_package(pybind11 CONFIG)`
  via `python -m pybind11 --cmakedir`, no vendored copy, so this is the only
  pybind11 declaration and it is what the isolated build env installs).
  `[project].version` is `dynamic`, read out of CMakeLists.txt's
  `project(LiPolGen VERSION 0.1.0 ...)` by
  `scikit_build_core.metadata.regex` — CMake's `project()` version is the
  one source of truth, not retyped in `pyproject.toml`.
- `wheel.packages = ["python/lipolgen"]`; `cmake.args =
  ["-DLIPOLGEN_BUILD_TESTS=OFF"]` (the option already existed, no new CMake
  option added); `build-dir = ".skbuild/{wheel_tag}"` (repo-local, gitignored,
  never the in-tree `build/`); `editable.mode = "redirect"` (`rebuild` left
  at its default/off — scikit-build-core's ninja/cmake are only guaranteed
  present in pip's *build*-isolation env, not in the venv doing the
  importing afterwards, so on-import rebuild is not reliable there).
  `LIPOLGEN_DEPS_PREFIX` honoured both via
  `--config-settings=cmake.define.LIPOLGEN_DEPS_PREFIX=…` (scikit-build-core,
  no extra wiring needed) and, as a convenience, from the environment
  (`CMakeLists.txt`: `if(NOT DEFINED CACHE{LIPOLGEN_DEPS_PREFIX} AND DEFINED
  ENV{LIPOLGEN_DEPS_PREFIX}) …`).
- `CMakeLists.txt`: `install()` is now branched on `SKBUILD`. The `SKBUILD`
  branch installs `_lipolgen` and all four `lipolgen_*` bridge libraries into
  **one** destination, `lipolgen/` (the package dir `wheel.packages` also
  populates) — the prototype's `lib/` + `lipolgen/` duplication is gone, one
  copy of each `.so`, confirmed via `python -m zipfile -l` on a built wheel.
  The non-`SKBUILD` branch keeps the previous `include/` +
  `${CMAKE_INSTALL_LIBDIR}` install for conventional `cmake --install
  --prefix …` use, untouched. `CMAKE_INSTALL_RPATH = "$ORIGIN:
  ${LIPOLGEN_DEPS_PREFIX}/lib"`, `CMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE` and
  `CMAKE_BUILD_WITH_INSTALL_RPATH=TRUE` all live inside `if(SKBUILD)` —
  **not** unconditional. That guard is load-bearing, not cosmetic: CMake
  pads a target's BUILD-tree RUNPATH with empty (`:`) entries at configure
  time for any target with an `install()` rule once `CMAKE_INSTALL_RPATH` is
  set at all, regardless of whether an install step ever runs, so setting it
  outside `if(SKBUILD)` — as an earlier revision briefly did — put
  CWD-lookup RUNPATH entries into `build/libLiPolGenCore.so` and the other
  in-tree `.so`s even though the in-tree flow never installs anything. With
  the fix, the non-`SKBUILD` configure sets no `CMAKE_INSTALL_RPATH` at all,
  so `build/lipolgen_tests` and the in-tree `build/python` module keep
  CMake's default build-tree RPATH exactly as before this packaging support
  existed (`readelf -d build/libLiPolGenCore.so` shows no RUNPATH entry, as
  at HEAD pre-packaging). `readelf -d` on the built wheel's `.so`s confirms
  `RUNPATH: $ORIGIN:<deps>/lib`; `ldd` (with `LD_LIBRARY_PATH` unset)
  resolves `libHepMC3.so.4`/`libpythia8.so`/`libLHAPDF.so` from the deps
  prefix alone.
- Data: `data/vmc` (1.1 MB) is installed to `lipolgen/data` inside the wheel
  under the `SKBUILD` branch; `python/lipolgen/__init__.py` points
  `$LIPOLGEN_DATA_DIR` at it when that env var is unset and the directory
  exists next to the installed module (no-op for the in-tree build, whose
  `build/python/lipolgen/` never has a `data/` sibling, so the compiled-in
  `$CMAKE_SOURCE_DIR/data` default keeps resolving exactly as before).
  PYTHIA8's `xmldoc` and LHAPDF's grids are **not** vendored (much larger,
  own licensing questions for the LHAPDF grids) — documented in
  `docs/USAGE.md` as coming from `$LIPOLGEN_DEPS_PREFIX` via
  `PYTHIA8DATA`/`LHAPDF_DATA_PATH`, same as `env.sh`.
- **Gate, measured 2026-09-02** (fresh venv, `numpy`+`pytest`+`pyhepmc`, all
  installed cleanly from PyPI):
  - `pip install -e .` with `LIPOLGEN_DEPS_PREFIX=…`: **66 s**, wheel
    `lipolgen-0.1.0-cp311-cp311-linux_x86_64.whl`.
  - From `/tmp` (no `env.sh`, no `PYTHONPATH`, only `PYTHIA8DATA`/
    `LHAPDF_DATA_PATH` exported): `import lipolgen; lipolgen.ion_spin('6Li')`
    → `1.0`; `python -m pytest python/tests -q` with the new
    `LIPOLGEN_TESTS_USE_INSTALLED=1` opt-out (`conftest.py`, default
    behaviour unchanged) → 144 passed, 2 skipped (missing optional `yaml`/
    `scipy` in the minimal venv — unrelated to packaging).
  - `pip wheel . --no-deps` + `python -m zipfile -l`: one `_lipolgen*.so`
    and one each of the four `libLiPolGen*.so`, all under `lipolgen/`, no
    duplicates.
  - In-tree flow re-verified bit for bit after all of the above:
    `build/lipolgen_tests` 276/276, `pytest python/tests` 146 passed.
- Portable wheel still needs `auditwheel repair` to vendor
  HepMC3/PYTHIA8/LHAPDF and rewrite RPATHs `$ORIGIN`-relative; documented in
  `docs/USAGE.md` together with the GPL-3 consequence of doing so (see §C in
  `engineering.md` and §13 below).
- License: GPL-3.0-or-later for LiPolGen's own code (the linked combination is
  GPL-3 regardless; permissive headers would mislead; MCnet-consistent).
  `pyproject.toml`'s `license = "GPL-3.0-or-later"` (SPDX string) matches.
