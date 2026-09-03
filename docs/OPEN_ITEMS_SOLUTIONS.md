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
| 8 | Spin-3/2 SF basis (plans/04 #14) | **exists** — Jaffe–Manohar NPB 321 (1989); explicit J=3/2 functions arXiv:2209.12161 Eqs. 19a–d; rank-≤2 truncation is *exact* for unpolarized-beam inclusive observables | 5–10 d note | theory note |
| 9 | Tensor-sector RC (plans/04 #10) | **implemented 2026-09-03** — `rc.hpp`/`rc.cpp`, opt-in `--rc tensor-band`: the two-sided band `rc_tensor_lo`/`rc_tensor_hi` on the tensor part of the rate (δ log-linear, 0.30 at x = 0.01 → 0.015 at x = 0.16) plus `rc_tail`, POLRAD's t-peak elastic tail with its tensor part and the unpolarised quasi-elastic tail. Weight-only, on `Event::rc_weights` and never on `Event::weight`; `--rc off` is byte-identical | done | numbers in §9 below |
| 10 | b₁ for A > 2 (plans/04 #9) | **implemented 2026-09-03 as an OPT-IN backend, and STILL OPEN** — `b1_nuclear.hpp`/`b1_nuclear.cpp`, `--b1-model li6-convolution`: the **four**-term α–d convolution on the Cosyn–Dong–Kumano–Sargsian kernel (the struck-α orbital term is not optional, it is ≈ 0.5 × the struck-d one). The default stays `Li6B1(MillerB1)`, bit for bit. **The A = 2 validation gate FAILS its magnitude clause** (a factor 2.27 below the digitized CDKS Fig. 4 peak at CDKS Eq. (21)'s δ-function, which the gate now defaults to — 3.68 at Eq. (17)'s κ = 1; the nucleon PDF, the one remaining identified item, closes it to 1.39), so **no ⁶Li number from it may be published** and the item does not close | 10–15 d, ~8 spent | code done, **gate open** — §10 below |
| 11 | Coherent ⁶Li amplitude (plans/04 #18) | **both halves in-tree** — eSTARlight ⁶Li unpolarized rates/slope (2026-09-02, `estarlight_li6.md`) settle `slope_b = 50 ± 10` GeV⁻² as citable and bracket `f0` one-sided, [1.0e-3, 3.0e-2]; the α+d configuration sampler (2026-09-03, `cluster_config.hpp`, `phase_G_numbers.md`) predicts ⁶Li's tensor a₂ ~10× smaller than the deuteron's and of **opposite sign** | done / done / collab | §11 below; Mäntysaari-group ask drafted, not yet sent |
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
  P₂-weighted density: 2(M_d/M_α)² = 0.5064 against a measured 0.489–0.503.
  A regression that silently drops it moves b₁ by ~30 %. The A = 2 validation
  was done and it **fails on magnitude** — see §10.

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

## 10. b₁ of ⁶Li — a four-term α–d convolution, opt-in, and the gate is NOT passed

**IMPLEMENTED (2026-09-03) AS AN OPT-IN BACKEND, AND THE ITEM STAYS OPEN**
(`include/lipolgen/b1_nuclear.hpp`, `src/core/b1_nuclear.cpp`;
`PipelineConfig::b1_model` / `--b1-model {miller,cdks,li6-convolution}` plus
`--b1-band-scale` and `--b1-alpha-d-dwave-weight`; `docs/USAGE.md` §2a; design
`docs/open_items/run_2026-09-02/design_D_b1_li6.md`; the measured gate in
`phase_D_gate.md` and the measured numbers in `phase_D_numbers.md`).

> ### ⚠ READ THIS BEFORE QUOTING ANY NUMBER BELOW
>
> Design §5 makes the A = 2 validation gate **blocking**: the same kernel, fed
> the AV18 deuteron u(k), w(k) instead of the α–d waves, must reproduce the
> digitized CDKS Fig. 4 before any ⁶Li number is quoted. **It fails clause
> G3b** — a factor **2.27** low at CDKS Eq. (21)'s δ-function, which is what
> the gate now defaults to (3.68 at Eq. (17)'s κ = 1 form), with the nucleon
> PDF the dominant remaining item.
>
> Under design §5.4's *Escalation* clause the backend stays in the tree
> behind its flag with the warning in its header and in `--help`, and **no ⁶Li
> number from it may be published**. The ⁶Li tables below are recorded so the
> phase is **reproducible and its regressions catchable** — they are not
> results, and every one of them additionally carries the mandatory ±100 %
> band. **Item 10 does not close.**

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
| **G3a** — exactly two sign changes, positions, peak position | **PASS**: zeros 0.0221 / 0.3774 against 0.0656 / 0.4572 (windows ±0.08 / ±0.10), peak 0.755 against 0.766 (±0.10). **Read the margin**: the low-x zero clears its counting window's edge (x = 0.02) by 0.0021 and, with CT18NLO, drops below the scan floor so the counting clause *fails* — it is not a robust discriminator; the second zero and the peak position are |
| **G3b** — peak magnitude within a factor 2 | **FAIL**: 4.7748e−4 against the digitized 1.08521e−3, **ratio 0.440**, a factor **2.27 low** (2.9484e−4, ratio 0.272, factor 3.68, at Eq. (17)'s κ = 1) |
| G3c — Close–Kumano, reported | this kernel +2.1537e−4 (κ = 1: +1.0792e−4); digitized CDKS +4.5920e−4; Miller +5.9147e−3. None is zero; none is enforced |

**The residual is attributed, not dangling.** The §5.4 checklist was worked in
order with no fudge factor tuned. Starting point: the κ = 1 column, ratio
0.2717.

| item | measured effect on the peak | direction |
|---|---|---|
| 0. finite-\|q⃗\| δ-function, κ = 1 → √(1+γ²) | ×1.620 | closes — **now IN**: the gate's default since 2026-09-03 |
| 4. nucleon PDF, `ToyF2` → `LhapdfSF("CT18NLO")` | ×1.669 | closes — **still open** |
| 1. target mass (CDKS Eq. 22, already in) | ×1.49 at x = 0.8 (1.52 at κ = 1) | in |
| 2. R (`r1998`, CDKS's own; already in **on the gate** — the ⁶Li backend defaults to `r_sigma_lt` instead, see below) | ×0.98 | in |
| 5. wave function, AV18 → rescaled to P_D = 4.85 % | ×0.88 | **opens** |
| **the gate's default today** | **ratio 0.440 (factor 2.27)** | **outside G3b** |
| **with item 4 as well, measured** | **ratio 0.719 (factor 1.39)** | **inside G3b** |

The two are very nearly independent: separately they predict
0.2717 × 1.669 × 1.620 = 0.735 against the measured 0.7193, a 2 % overlap, so
there is no third unexplained factor at the peak. With CT18NLO the second zero
also lands at **0.449** against the digitized 0.457 and the mid-x dip at
−1.70e−4 against −1.77e−4 — a 4 % agreement where `ToyF2` was a factor 5.7 low.
**Item 0's default was changed by the review of 2026-09-03 and item 4 was
not**: Eq. (21) is CDKS's own exact δ-function and Eq. (17) is the "≃" of it
(design §5.4 item 0: a gate-driven change there "is a finding, not a fudge"),
while item 4 needs LHAPDF, which the core must not link. `Li6ConvolutionB1`'s
own `finite_q_delta` stays `false` (A10) — in ⁶Li the switch is worth −1 % to
+7 %. What remains unexplained is the **low-x tail**, where even CT18 + κ gives
x·b₁(0.10) = −3.8e−5 against the digitized −1.73e−5.

### The ⁶Li numbers — RECORDED, NOT PUBLISHED (see the warning above)

x·b₁ per nucleon, Q² = 2.5 GeV², default options — which include **ToyF2 F₁
with `r_sigma_lt`**, not the gate's `r1998` (below), and the
2400 / 2400 / 3200 y grid the review of 2026-09-03 set:

| x | (1) embedded d, S | (3) CG D-wave | (2d) struck d | (2α) struck α | **total** | `Li6B1(MillerB1)` | `Li6B1(CdksB1)` |
|---|---|---|---|---|---|---|---|
| 0.05 | +1.5686e−6 | +3.0979e−9 | +1.3838e−6 | +6.9528e−7 | **+3.6507e−6** | +3.6503e−4 | +8.9916e−7 |
| 0.10 | −4.6926e−6 | −9.2965e−9 | +1.2673e−6 | +6.3805e−7 | **−2.7966e−6** | +4.1222e−4 | −2.6566e−6 |
| 0.20 | −2.3664e−5 | −4.6736e−8 | +1.6669e−6 | +8.3886e−7 | **−2.1205e−5** | +2.6700e−4 | −1.3491e−5 |
| 0.30 | −4.4967e−5 | −8.8342e−8 | +2.2133e−6 | +1.1092e−6 | **−4.1732e−5** | +6.7167e−6 | −2.5964e−5 |
| 0.50 | +4.2448e−5 | +8.5642e−8 | +6.1632e−6 | +3.0422e−6 | **+5.1739e−5** | −1.8650e−4 | +2.2792e−5 |

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
and `Li6B1(CdksB1)`'s +6.9298e−5. `finite_q_delta` moves ⁶Li's x·b₁ by only
−1 % to +7 % (against 1.57–1.69 on the A = 2 gate), because in ⁶Li the term κ
multiplies is the small orbital one — which is why the two objects default
differently.

**The ⁶Li backend and the gate use different R, and that is a choice, not an
oversight.** `Li6ConvolutionOptions::r_func` null ⇒ **`r_sigma_lt`** (the
kernel's own toy R, for consistency with `InclusiveKernel`'s F₁);
`DeuteronConvolutionB1::Options::r_func` null ⇒ **`r1998`** (CDKS's SLAC world
fit, because the gate reproduces *their* figure). So the R the gate was
validated with is not the R the shipped ⁶Li numbers above carry. Measured on
the gate at κ = 1, like for like: r1998 → r_sigma_lt moves the second zero from
0.392 to **0.365**, *away* from the digitized 0.457, and the peak up by
**1.6 %**. Inside the band, not cosmetic; whether the ⁶Li default should become
`r1998` is on the close-out list.

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

* **The mandatory 100 % band.** Q(⁶Li) = −0.0806(6) fm² against
  Q_d = +0.2859(3) fm²: the α–d D wave enters the closest measured observable
  with the opposite sign to the deuteron's own D state and nearly cancels it.
  *(Q(⁶Li): Pyykkö, Mol. Phys. **106** (2008) 1965, whose ⁶Li entry is
  Cederberg et al., Phys. Rev. A **57** (1998) 2539; Q_d: Bishop & Cheung,
  Phys. Rev. A **20** (1979) 381.)* Every number is `--b1-band-scale 0/1/2`,
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
* **Q1b, a separate item that this work routed around and did NOT change.**
  `b1_convolution()` multiplies the raw digitized column by
  `B1_PER_DEUTERON_TO_PER_NUCLEON = 0.5`, so `CdksB1` is a factor 2 low against
  a curve that CDKS Eq. (10)'s explicit 1/A, their text under Eq. (16) and their
  Fig. 6 all say is **already per nucleon**. The new backend reaches the raw
  column through `cdks_b1_raw_per_nucleon()`; the constant is untouched, because
  every published number carries the current convention. Two close-out items:
  (a) decide whether `CdksB1` should stop halving; (b) check **Miller's** own
  normalisation independently — the answer is on the axis of Miller's Fig. 5 and
  is not assumable from CDKS.
* **Q6.** `LI6_B1_RANK2_TRANSFER` = 0.921947 is tied to the Hulthén-scenario
  P_D = 0.0867; the VMC value 0.0193549 would give 0.982581. Changing it would
  move every published inclusive tensor number, so it stays — but this backend
  makes the inconsistency visible (the tagged channel already uses VMC at
  `--cluster-wave vmc`; the inclusive b₁ never does).
* **Q8 — there is nothing to validate against for A = 6.** No published b₁
  exists for any A > 2. The A = 2 gate plus the analytic limits (T5, T6, T7) are
  the *only* checks there are, which is why the gate is blocking and why the
  band is not a formality.

### What has to happen before item 10 can close

1. ~~Decide item 0's default~~ — **done, 2026-09-03**: the A = 2 gate is
   quoted at CDKS Eq. (21)'s δ-function, which is their exact definition; the
   ⁶Li backend keeps Eq. (17)'s κ = 1, where the switch is worth 1 %.
2. Get MSTW2008 LO, or the tabulated curves from S. Kumano's group (design Q2),
   and rerun checklist item 4 with CDKS's own PDF rather than CT18NLO. **This
   is the only identified item still open at the peak.**
3. Get a real CD-Bonn u, w instead of the D-state rescaling proxy of item 5.
4. Decide whether `Li6ConvolutionOptions` should default to `r1998` like the
   gate rather than to `r_sigma_lt` (checklist item 2).
5. Only then re-open G3b. **Until it passes, no ⁶Li number ships.**

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
band as a scenario**, record eSTARlight's one-sided bracket
[1.0 × 10⁻³, 3.0 × 10⁻²] in `coherent.hpp`, and do not treat it as a
determination.

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
`phase_G_numbers.md`): ⟨r²⟩(⁶Li) = **6.4447 fm²** (r_rms 2.5386 fm, 4 % above
the measured point radius 2.4655 fm); the point-matter quadrupole per
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
> point radius to 4 % but overshoots Q(⁶Li) by a factor ≈ 7.5 (measured
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
wave functions.** The α+d truncation reproduces the ⁶Li radius to 3–4 % but
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
