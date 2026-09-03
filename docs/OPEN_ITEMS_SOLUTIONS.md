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
| 10 | b₁ for A > 2 (plans/04 #9) | **first-mover** — nothing exists; three-term α–d convolution on the Cosyn–Dong–Kumano–Sargsian kernel; 100 % band mandatory (⁶Li quadrupole puzzle) | 10–15 d | design only |
| 11 | Coherent ⁶Li amplitude (plans/04 #18) | **route changed** — Sartre ruled out (hard-coded nuclei, no polarization axis); use eSTARlight for unpolarized rates (1 d) and `hejajama/subnucleondiffraction` (code of arXiv:2408.13213) with an α+d configuration sampler for the tensor cos 2φ | 1 d / 10–15 d / collab | Mäntysaari-group ask |
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
