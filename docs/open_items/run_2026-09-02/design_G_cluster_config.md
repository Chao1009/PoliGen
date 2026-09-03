# Design G(i) — the α+d configuration sampler for polarized ⁶Li

`include/lipolgen/cluster_config.hpp` + `src/core/cluster_config.cpp`
· Phase G(i) of `docs/open_items/run_2026-09-02/PLAN.md`
· open item 11, `docs/OPEN_ITEMS_SOLUTIONS.md` §11
· authored 2026-09-02, design only — no code written, no existing file touched.

---

## 0. What this module is, and what it is not

**Is.** A generator of *nucleon-position configurations* of a ⁶Li nucleus in a
chosen magnetic substate `m ∈ {+1, 0, −1}` (or unpolarized) along a chosen
quantization axis: N configurations × 6 nucleons × (x, y, z in fm, isospin).
The configurations are what a Good–Walker dipole-model code
(`hejajama/subnucleondiffraction`, the code of arXiv:2408.13213) averages the
scattering amplitude over to obtain the coherent |t| distribution and its
cos 2Φ tensor modulation.

**Is not.** It is not a diffractive amplitude, not a cross section, and not a
new physics default anywhere in LiPolGen. Nothing in the existing generator
calls it. `CoherentScenario` keeps its scenario numbers unchanged
(`coherent.hpp`); this module *supplies a number that can be compared to
`eps_b0`*, and §2.8 does exactly that, but the comparison is a document, not a
code path.

**The honest framing, up front.** The α+d two-body truncation of the ANL VMC
overlap reproduces the ⁶Li *size* to 3–4 % and overshoots the ⁶Li *quadrupole*
by a factor 7.5–8.9 against measurement (§2.7, §8). The sampler is therefore a **geometry generator whose
tensor output is a model input carrying an explicit factor-of-several band**,
not a prediction of Q(⁶Li). `docs/OPEN_ITEMS_SOLUTIONS.md` §11's rule — do not
derive tensor inputs from these wave functions — is *quantified* here rather
than repealed.

---

## 1. Conventions inherited (nothing is redefined)

| quantity | single copy already in the repo | value |
|---|---|---|
| ħc | `constants.hpp:51` `HBARC_GEV_FM` (aliased `GEV_PER_FM_INV`, `coherent.hpp:34`) | **0.19733** GeV fm |
| data root | `cluster.hpp:41` `data_dir()` / `data_path()` | `$LIPOLGEN_DATA_DIR`, else compiled-in `${CMAKE_SOURCE_DIR}/data` |
| Condon–Shortley Θ_l^m | `cluster.hpp:173` `theta_lm(l, m, cosθ)` | Y_l^m = Θ_l^m(θ) e^{imφ}, l ≤ 2 |
| Clebsch–Gordan | `spin.hpp:62` `clebsch_gordan` (Racah, Varshalovich) | — |
| ANL overlap reader | `cluster.hpp:59` `read_anl_overlap` | `[0]` k-block, `[1]` r-block, two signed columns + MC errors |
| α+d D-state probability | `tagged.hpp:72` `VMC_P_D_LI6` | 0.019353 (`momenta/li6_ad1`) |
| α+d spectroscopic factor | `tagged.hpp:74` `VMC_S_ALPHA_D_LI6` | 0.81971 |
| α Gaussian point-nucleon size | `fsi.hpp:177` / `fsi.cpp:47` `cluster_point_a2_fm2(2,4)` | a² = 0.7001 fm², ⟨r²⟩ = 3a² = 2.1003 fm², rms 1.4492 fm |
| Q(⁶Li), measured | `rc.hpp:290` `LI6_QUADRUPOLE_FM2` | **−0.0818 fm²** (TUNL A=6, 1998CE04) |
| ⟨r²⟩_point(⁶Li) | `rc.hpp:248–252` carries it only in prose and through `LI6_FF_HO_A_FM`/`LI6_FF_HO_ALPHA`; **this design adds the named constant `LI6_R2_POINT_FM2` to `rc.hpp`** next to `LI6_QUADRUPOLE_FM2` (r_ch = 2.589(39) fm, Angeli–Marinova ADNDT 99 (2013) 69, minus the nucleon terms) | 6.0788 fm² → **r_point = 2.4655 fm** |
| r_point(⁶Li), VMC | **new** `cluster_config.hpp` `LI6_R_POINT_VMC_FM` (`li6.density` header, AV18+UX) | 2.4433 fm |
| Gaussian \|F\|² slope | `coherent.hpp:83` `gaussian_slope(r_rms)` = r²/3 | B(⁶Li) = 52.0 GeV⁻² at r = 2.4655 fm |
| digitized deuteron a₂ | `coherent.cpp:58` `mantysaari_a2_deuteron()` | (\|t\|, a₂(0), a₂(±1)) at \|t\| = 0.05…0.30 |
| RNG | `rng.hpp` `Rng(seed, run, bunch, event)` — counter-based, never global | — |
| interp / quadrature | `numerics.hpp` `np_interp`, `trapezoid`, `linspace` | — |

**Three notes the implementer must not slip on.**

0. **ħc is 0.19733, not 0.1973269804.** `HBARC_GEV_FM` (`constants.hpp:51`, "THE
   one conversion between fm^-1 and GeV") is 0.19733; 0.1973269804 is
   `cluster.cpp:23`'s file-local `kHbarCGeVfm`, a **pre-existing second copy**
   that this module must not propagate. **Rule for `cluster_config.cpp`: every
   fm ↔ GeV conversion goes through `HBARC_GEV_FM` / `GEV_PER_FM_INV`, never
   through a literal**, so that `eps_b0_equivalent` is strictly comparable with
   `gaussian_slope` and `CoherentScenario`. `HBARC_GEV_FM` is **not** exposed to
   Python, so the CLI must take every GeV⁻² quantity from the C++ predictors and
   convert nothing itself. (The §2.8 and §8 numbers below were computed at
   0.19732698; at 0.19733 they move by 2.8 × 10⁻⁵ relative — immaterial, and not
   re-tabulated.)
1. The task brief quotes Q(⁶Li) = −0.0806 fm². That is the Pyykkö compilation
   value; the repo's single copy is TUNL's **−0.0818** (`rc.hpp:290`, which
   documents the difference explicitly). Use `LI6_QUADRUPOLE_FM2`. The two
   differ by 1.5 % and nothing here resolves at that level.
2. The task brief quotes r_ch(⁴He) = 1.678 fm (muonic helium, Krauth et al.,
   Nature 589 (2021) 527 → 1.67824(83)). The repo's single copy is
   `fsi.cpp:53–55`'s **1.6755** with the proton folded out as
   a² = (r_ch² − r_p²)/3. Do **not** add a second ⁴He radius. The 0.16 %
   difference moves the α core rms by 0.002 fm, far below the ≈ 5 % ⁶Li radius
   error this module already carries.

Frame and axis. Positions are in the **ion rest frame**; the quantization axis
n̂(θ_S, φ_S) is the same axis `spin.hpp` and `tagged.hpp` use. The sampler
generates in the *spin frame* (axis = +z) and rotates by R_z(φ_S) R_y(θ_S) at
the end, exactly as `boost_spectator` does for the tagged channel. For the
coherent cos 2Φ the axis must have a **transverse** component: with the axis
along the beam the projected density is azimuthally symmetric and a₂ ≡ 0
(arXiv:2408.13213, Fig. 2(c)).

---

## 2. Physics

### 2.1 The state

⁶Li(1⁺) is written as a rigid α(0⁺) core plus a deuteron, with relative
orbital L = 0, 2 coupled to the deuteron spin S = 1:

    |⁶Li; 1 M⟩ = Σ_{L=0,2} R_L(R) [ Y_L(R̂) ⊗ φ_d ]^{1M}
               = Σ_{L,m_S} ⟨L (M−m_S) 1 m_S | 1 M⟩ R_L(R) Y_L^{M−m_S}(R̂) φ_d^{(m_S)}(r⃗)

with R⃗ = (deuteron c.m.) − (α c.m.), r⃗ = r⃗_p − r⃗_n, and the deuteron's own
internal state carrying the *same* (L = 0, 2) ⊗ (S = 1) → J = 1 structure:

    φ_d^{(m_S)}(r⃗) = Σ_{λ=0,2} ⟨λ (m_S−σ) 1 σ | 1 m_S⟩ (χ_λ(r)/r) Y_λ^{m_S−σ}(r̂) |1σ⟩,
    χ_0 = u,  χ_2 = w.

This is **the same coupling `tagged.hpp` already implements** — its joint
amplitude `A_{m_S}(M; k, k̂) = Σ_L a_L ψ̂_L(k) C_L(M,m_S) Y_L^{M−m_S}(k̂)`
(`tagged.hpp:12–17`) with C_L = ⟨L (M−m_S) S_c m_S | J M⟩ — read in **r space
instead of k space**, and *up to the (−i)^L phase of the partial-wave Fourier
transform*, which is not optional bookkeeping:

    ψ_L^phys(k) = (−i)^L 4π ∫ j_L(kr) R_L(r) r² dr,

so the physical k-space S–D interference term carries the **opposite sign** to
the plain (no-i^L) Fourier–Bessel transform. The ANL k-block, `VmcRadial` and
`vmc_reconciliation.md`'s "sign(ψ₂/ψ₀) = −1 below the S node" are all the
**no-phase real transform**; the sibling `design_D_b1_li6.md` §1.2 applies
φ₂ = −W to the deuteron for exactly this reason (its convention table, line 60).
**This module works entirely in r space and therefore needs no phase at all** —
but an implementer cross-checking this sampler's oblate m = ±1 (𝒬 < 0) against
the tagged channel's positive A_zz must not "fix" one of the two to match the
other. §3.1's sign trap and §7 T4 make this a gate. (Out of scope here, but
worth flagging to `tagged.hpp`'s owner: `li6_vmc_waves()` carries no i^L.)

The recoupling is not re-derived and not re-coded: the module calls
`clebsch_gordan` and `theta_lm`.

Positions, with the total c.m. at the origin (A_α = 4, A_d = 2, A = 6):

    r⃗_α,i = −(1/3) R⃗ + s⃗_i     (i = 1…4,  Σ_i s⃗_i = 0)
    r⃗_p   = +(2/3) R⃗ + r⃗/2
    r⃗_n   = +(2/3) R⃗ − r⃗/2

### 2.2 The master angular density (one formula, three uses)

For any (L = 0, 2) ⊗ (S = 1) → J = 1 system with radial functions f₀, f₂
normalized as ∫(f₀² + f₂²) x² dx = 1, the probability density of the relative
vector in substate m, with the quantization axis along +z, is

    ρ_m(x, cosθ) = Σ_{m_S=+1,0,−1} | Σ_{L=0,2} C_L(m,m_S) f_L(x) Θ_L^{m−m_S}(cosθ) |²,
    C_L(m,m_S) = ⟨L (m−m_S) 1 m_S | 1 m⟩,                                     (G1)

with ∫ ρ_m x² dx dΩ = 1 and **no φ dependence** (|e^{i m_L φ}|² = 1).

Written out for the deuteron (f₀ = u/r, f₂ = w/r) this is *exactly*
arXiv:2408.13213 Eqs. (7)–(8):

    ρ_d(r,θ; m_S=±1) = (1/16π) [ 4U² − 2√2 (1−3cos²θ) UV + (5−3cos²θ) V² ]     (G2)
    ρ_d(r,θ; m_S= 0) = (1/8π)  [ 2U² + 2√2 (1−3cos²θ) UV + (1+3cos²θ) V² ]     (G3)
    U = u(r)/r,  V = w(r)/r,  ∫(u²+w²) dr = 1.

**Verified**: (G1) evaluated through the repo's own `theta_lm` and
`clebsch_gordan` reproduces (G2)/(G3) to machine precision at every cosθ
tested (−1, −0.6, 0, 0.35, +1) — this is a doctest (§7 T1). Implement (G1);
keep (G2)/(G3) only as the test's independent closed form.

The same (G1) with f₀ = R₀, f₂ = R₂ is the **α–d relative** density
ρ_{αd}(R, cosθ_R) — the ⁶Li m-state anisotropy proper.

### 2.3 The quadrupole master formula (and its validation)

Integrating (G1) against (3cos²θ − 1):

    𝒬[f₀,f₂] ≡ ⟨3x_z² − x²⟩_{m=+1} = (1/5) ∫ x⁴ ( 2√2 f₀ f₂ − f₂² ) dx        (G4)

and ⟨3x_z²−x²⟩_m = (3m²−2) 𝒬, i.e. m = ±1 → 𝒬, m = 0 → −2𝒬, sum over m = 0.
𝒬 is invariant under f_L → −f_L for **both** L (the unobservable global
phase); only the *relative* sign of f₀ and f₂ is physical.

Two independent validations, both done during this design pass:

* **Against the AV18 deuteron's own printed quadrupole.**
  𝒬[u/r, w/r]/4 = **0.269670 fm²** vs the `fdeut.av18` header's
  `qm = 0.269673`. (The /4 is because Q_d is the *proton-only* moment and the
  proton sits at r⃗/2.) The same read of the file also reproduces
  `dstate = 0.057599` and `rd = 1.967364` exactly. This pins the sign
  convention: **fdeut.av18's u > 0 and w > 0 fed straight into (G1)/(G2) give
  the measured positive Q_d — no sign flip anywhere.**
* **Against 2-D quadrature of (G1)** on a toy (f₀, f₂) pair: 0.8206111 vs
  (G4)'s 0.8205472, 8 × 10⁻⁵ relative.

### 2.4 ⁶Li point-matter quadrupole from the cluster geometry

Summing 3z_i² − r_i² over the six nucleons with §2.1's positions, using
Σ_i s⃗_i = 0 and the α's 0⁺ (spherical, zero internal quadrupole):

    Q_matter(⁶Li, M) = (3M²−2) [ (4/3) 𝒬[R₀,R₂] + 2 Q_d D_T ]                  (G5)
    D_T ≡ ⟨3 m_S² − 2⟩_{M=+1} = P(+1) + P(−1) − 2P(0),
    P(m_S) = Σ_L N_L |⟨L (1−m_S) 1 m_S | 1 1⟩|²,   N_0 = 1−P_D^{αd}, N_2 = P_D^{αd}.

D_T is the *same quantity* as `TaggedModel::tensor_dilution()`
(`tagged.hpp:176`), but **it is not the same number by default and the two must
not be gated as equal**:

* `TaggedModel::tensor_dilution()` is a **grid quadrature** — `population_integrated`
  (`tagged.cpp:340`) is a pairwise sum of |A_{m_S}|² k² over the 280 × 96 (k, cosθ)
  table, not a closed form. `constants.hpp:85–91` already records the size of the
  gap for the sibling constant: the closed form gives 0.921970 against the
  0.9219467 quadrature, and the repo's own rule is that the two are **"pinned to
  each other at 1e-4 and not asserted equal."** Follow that precedent.
* `li6_alpha_channel(..., VmcAV18)` fixes P_D at `VMC_P_D_LI6` = **0.01935**
  (the `momenta/li6_ad1` file) and *ignores* its `p_d` argument
  (`tagged.hpp:113–115`), whereas this sampler's default P_D is the `li6.ad`
  r-block's **0.02011**. D_T differs by 7 × 10⁻⁴ between the two (0.98259 vs
  0.98190) for that reason alone.

So §7 T10 gates two different things: the **exact CG identity**
D_T = 1 − 0.9 P_D^{αd} at 1e-12 (it follows from P(+1) = N₀ + N₂/10,
P(0) = 3N₂/10, P(−1) = 6N₂/10), and agreement with a `TaggedModel` built at the
*sampler's own* P_D through the Hulthén path where P_D is still a knob,
`li6_alpha_channel(BETA_DEFAULT, sampler.p_d_alpha_d())`, at 1e-4.

(G5) is **exact under the diagonal truncation of §2.6**, because both operators
are one-coordinate operators and ∫d³R A_{m_S}A*_{m_S'} = δ_{m_S m_S'} P(m_S) by
orthogonality of Y_L^{m_L}. So the sampled quadrupole is a clean MC-statistics
gate, not an approximation test.

Analogously (and also exact, no cross term because Σ_i s⃗_i = 0 and
⟨R⃗·r⃗⟩ = 0):

    ⟨r²⟩(⁶Li) = (2/3)⟨s²⟩_α + (1/12)⟨r_np²⟩ + (2/9)⟨R²⟩_αd                     (G6)

### 2.5 Numbers this geometry gives (computed here, from the committed tables)

| input | source | value |
|---|---|---|
| ⟨r_np²⟩ | `fdeut.av18` | 15.4817 fm² (r_d = 1.96734 fm) |
| P_D(d) | `fdeut.av18` | 0.057599 |
| Q_d | (G4)/4 | +0.269670 fm² |
| S_αd = ∫(R₀²+R₂²)R²dR | `li6.ad` r-block | 0.85423 (file prints `ndx 0.856`, `s 0.838`, `d 0.017`; recomputed 0.83705 / 0.01718) |
| P_D(αd) | `li6.ad` r-block | 0.02011 (`li6.adr.fit`: 0.02542; `momenta`: `VMC_P_D_LI6` = 0.01935) |
| ⟨R²⟩_αd | `li6.ad` r-block | 17.548 fm² → rms **4.189 fm** (`li6.adr.fit`: 16.956, 4.118 fm) |
| 𝒬[R₀,R₂] | (G4) | **−1.4009 fm²** raw (interference −1.3325, pure-D −0.0684); −1.4895 from the smoothed fit |
| D_T | (G5) | 0.98190 at P_D = 0.02011 |
| ⟨s²⟩_α | `he4.density` | 2.0774 fm² → rms **1.4413 fm** (file prints 1.4404) |

⇒ **Q_matter(⁶Li, M=+1) = (4/3)(−1.4009) + 2(0.26967)(0.98190) = −1.868 + 0.530
= −1.338 fm²**, i.e. Q_charge = −0.669 fm² (N = Z, point nucleons).
⇒ **⟨r²⟩ = 6.575 fm² → r_point = 2.564 fm** (G6).

Those two lines are the worked example for the **raw** `li6.ad` r-block; the
module's *default* source combination (`FitRescaled`, §3.1) gives
Q_matter = −1.231 fm² and r_point = 2.539 fm. §8 tabulates all three variants,
and every gate in §7 is written against the default.

Anisotropy budget: the α–d term is **77.9 %** of |Q|, the deuteron's intrinsic
quadrupole 22.1 %; within the α–d term the S–D *interference* is **95.1 %** and
the pure D² term 4.9 %. Everything is D-wave-driven — with f₂ ≡ 0 in both
places Q ≡ 0 — and it is **linear in the D amplitude, not in P_D**.

### 2.6 The one approximation, named and sized

The joint density ρ_M(R⃗, r⃗) contains coherences between different deuteron
projections m_S, of the form A_{m_S} A*_{m_S'} with m_S ≠ m_S'. Because they
carry e^{i(m_L−m_L')(φ_R−φ_r)}, they are precisely a **correlation between the
azimuth of R̂ and the azimuth of r̂ about the quantization axis**.

The sampler drops them: it draws m_S from P(m_S|M), then **(R, cosθ_R) from
the m_S-CONDITIONED density** |A_{m_S}(M; R, cosθ_R)|², then (r, cosθ_r) from
(G1) with (u/r, w/r) **for that same m_S**, with independent uniform azimuths.
This is *the identical truncation `tagged.hpp` already makes and documents*
("diagonal truncation: coherences are dropped, documented", `tagged.hpp:23`),
so the two channels stay one model.

**The conditioning is not optional.** The diagonal joint density is

    ρ_M(R⃗, r⃗) = Σ_{m_S} |A_{m_S}(M; R⃗)|² Σ_σ |φ^{m_S}_σ(r⃗)|² ,               (G1′)

i.e. R̂ and r̂ are *both* correlated to the drawn m_S; only their relative
azimuth is thrown away. Drawing R̂ from the **m_S-summed** density and then r̂
from the drawn m_S's deuteron density decorrelates R̂ from m_S and is a
different model. For M = +1 the α–d orientation is pure L = 2 in two of the
three branches: m_S = −1 (P = 0.6 P_D) forces m_L = 2, so R̂ follows
|Y₂²|² ∝ sin⁴θ_R, and m_S = 0 (P = 0.3 P_D) forces m_L = 1, |Y₂¹|². This is
exactly why `TaggedModel::sample_kc(m_ion, m_s, ...)` (`tagged.hpp:185`) takes
m_S as an argument and draws from the per-m_S table `build_amp2` builds; mirror
that structure. (G5)/(G6) are **blind** to the error — only the m_S marginal and
the R marginal enter them — so no moment test would catch it; what it moves is
the p/n one-body density and the R̂–r̂ correlation, i.e. precisely the
transverse geometry the dipole model integrates. §7 T21 gates it with two exact
numbers: ⟨P₂(cosθ_R) | M=+1, m_S=−1⟩ = **−2/7** and
⟨P₂(cosθ_R) | M=+1, m_S=0⟩ = **+1/7** (from ⟨3cos²θ−1⟩ = 2[L(L+1)−3m_L²] /
[(2L+3)(2L−1)] at L = 2), which an m_S-marginal draw fails by construction
(a 2 × 10⁵-configuration MC of the corrected algorithm reproduces −0.2870 and
+0.162 here, against an m_S-summed ⟨P₂⟩ of −0.0355).

Size: the dropped terms are O(√(P_D^{αd}) × √(P_D^{d})) = 0.143 × 0.240 ≈
**3.4 %** of the joint angular correlation. They do **not** affect ⟨r²⟩ (G6) or
Q (G5) at all (§2.4). They can only move the p/n one-body density at the
per-cent level, through the mixed coordinate (2/3)R⃗ ± r⃗/2.

`ClusterConfigOptions::exact_coherence` is reserved as an **off-by-default,
not-implemented-in-v1** flag; if it is ever wanted, the exact route is a 5-D
(R, cosθ_R, r, cosθ_r, φ_r−φ_R) density and Metropolis rather than inverse-CDF.
Do not implement it in this phase; do declare the flag and throw on `true`, so
the approximation is visible in the API rather than buried in a comment.

### 2.7 Where the model fails, quantitatively

| observable | α+d sampler (default `FitRescaled`) | full α+d range (§8) | ab-initio | measurement |
|---|---|---|---|---|
| r_point(⁶Li) | 2.539 fm | 2.538–2.564 fm | 2.4433 fm (`li6.density`, VMC AV18+UX) | 2.4655 fm (from r_ch = 2.589(39)) |
| Q_charge(⁶Li) | **−0.615 fm²** | −0.615…−0.730 fm² | −0.20(6) fm² (**GFMC AV18+IL7**) | **−0.0818 fm²** (`LI6_QUADRUPOLE_FM2`) |

**Source of the −0.20(6).** Pastore, Pieper, Schiavilla & Wiringa, PRC 87,
035503 (2013), [arXiv:1212.3375](https://arxiv.org/abs/1212.3375) — the repo's
own record of it is `docs/open_items/physics_literature.md:258`. It is **GFMC
with AV18+IL7**, *not* VMC and *not* the AV18+UX Hamiltonian of `li6.density`
and `li6.ad` that this sampler reads, so it is an independent-Hamiltonian
comparison point, not a consistency check on our own tables. Label it
"GFMC AV18+IL7" everywhere it appears — the band, the sidecar caveat, §8, §10
and the collaboration letter — because "independent VMC" is a wrong
attribution to send to the group that published it.

The radius overshoots by **+3.0 %** against measurement and **+3.9 %** against
`li6.density` (⟨r²⟩ by +6.0 % and +7.9 %); the quadrupole overshoots by a
factor **7.5** against measurement and **3.1** against the GFMC AV18+IL7 value
(**8.9** and **3.6** for the noisier `FitRaw` shapes). Both have the same
root: S_αd = 0.854, so 15 % of the ⁶Li wave function is not α+d at all, and the
model gives the two clusters their *free* internal wave functions inside a
nucleus that compresses and polarizes them.

**The asymptotic D/S ratio, done properly.** R₂/R₀ at finite R is *not* the
asymptotic D/S ratio η ≡ C₂/C₀, which is defined through the Whittaker
functions of the Coulomb-bound α–d channel,

    R_L(R) → C_L W_{−η_c, L+1/2}(2κR) / R ,
    κ = √(2μE_sep)/ħc = 0.3074 fm⁻¹   (E_sep = 1.4743 MeV, μ = 1247.7 MeV),
    η_c = Z_α Z_d α μ /(ħc κ) = 0.3002   (Sommerfeld parameter),

and W₂/W₀ is **3.34 / 2.92 / 2.63** at R = 6 / 7 / 8 fm — an order-unity factor
the naive ratio misses entirely. Dividing it out (recomputed in this pass from
the committed tables):

| R [fm] | R₂/R₀ (`FitRescaled`) | η | R₂/R₀ (`OverlapRaw`) | η |
|---|---|---|---|---|
| 6 | −0.139 | −0.042 | −0.145 | −0.043 |
| 7 | −0.148 | −0.051 | −0.139 | −0.048 |
| 8 | −0.139 | −0.053 | −0.189 | −0.072 |

so **η ≈ −0.05(1)**, not −0.15. And the literature value is **not** unsettled:
it is measured — **η = −0.025 ± 0.006 ± 0.010** (George & Knutson, PRC 59, 598
(1999), from d + α elastic tensor analysing powers; other determinations lie in
−0.01…−0.03). The table's D-wave tail is therefore **≈ 2× the measured η**, not
5–15×. That is a genuine but *moderate* D-wave excess, entirely compatible with
the missing 15 % non-α+d component and with cluster distortion — it does **not**
indicate a convention or normalization error in the ANL L = 2 column. §10 O1 is
re-worded accordingly. Expose the number as `asymptotic_ds_ratio()` on the
analytic layer, with a doctest that reproduces −0.05 from the `FitRescaled`
tables (κ, η_c and the 6–8 fm window fixed in the header, the Whittaker ratio
from the standard `W_{−η,L+1/2}(z) = e^{−z/2} z^{L+1} U(L+1+η, 2L+2, z)`).

**Which one is the sampler expected to reproduce? Neither, and that is the
point.** The sampler must reproduce, to MC precision, the quadrupole *of the
wave functions it was built from* — (G5) = −1.231 fm² for the default source combination — because that is a
statement about the code, not about ⁶Li. The physics number is carried
separately: the module exposes `quadrupole_band_fm2()` returning
{measured −0.0818, GFMC AV18+IL7 −0.20(6), **α+d model = the configured
source's own** `q_matter_analytic_fm2(1)/2`} and the writer records which one
the run used, because §11's rule stands: **do not derive a tensor input for a
published observable from these wave functions.** The third entry must **not**
be a hard-coded −0.669 (that is the `OverlapRaw` value, while the default source
gives −0.615 and the §5.2 sidecar caveat prints −0.62 — three numbers for one
quantity); it is computed from the *configured* source, after any
`alpha_d_scale` / `quadrupole_target_fm2` dial, and the dial is stamped beside
it. Prose quotes the full model range **−0.615…−0.730 fm²**, not one number.

**The `quadrupole_target` dial, and why √(target/model) does not close.** 𝒬 is
**linear** in the D amplitude through the interference term (95 % of it), and
the deuteron term 2 Q_d D_T = +0.530 fm² is an offset that barely moves with the
dial. Scaling the α–d D wave by s (and renormalizing) gives

    n(s) = 1 − (1 − s²) P_D^{αd},        P_D′(s) = s² P_D^{αd} / n(s),
    2 Q_target = (4/3)·[ s·𝒬_int + s²·𝒬_DD ] / n(s)
               + 2 Q_d · [ 1 − 0.9 P_D^{αd} s² / n(s) ] ,                  (G9)

a **quadratic in s** (bracketed root-find, or the closed form; `𝒬_int` and
`𝒬_DD` are the split §2.5 already quotes and the header exposes). The naive
√ rule fails by 7×: for the default `FitRescaled` source
(𝒬_int = −1.2570, 𝒬_DD = −0.0631, P_D = 0.020111, Q_d = 0.269670),
s = √(0.0818/0.6154) = 0.3646 gives Q_charge = **−0.048 fm²**, not −0.0818 —
outside even the loosened T17 gate. The root of (G9) at
Q_target = `LI6_QUADRUPOLE_FM2` is **s = 0.4032**, and the D-state probability it
implies is **P_D′ = 3.3 × 10⁻³** (the √ rule would give 2.7 × 10⁻³; the "3 × 10⁻⁴"
of the first draft was wrong by an order of magnitude). Physically absurd as a
wave function, honest as a *deformation dial*, and it must be labelled as one in
the output header.

**The dial has a floor and `validate()` must enforce it.** At s = 0 the α–d
term vanishes and the deuteron's own moment survives:
Q_charge(s=0) = +Q_d = **+0.270 fm²**. The reachable range with s ∈ [0, 1] is
therefore Q_charge ∈ [−0.615, +0.270] fm² for the default source; **a
`quadrupole_target_fm2` above +0.270 (or below the s = 1 value) is unreachable
and `validate()` throws** rather than silently returning the nearest root.

**Systematics list for the model quadrupole** (each already sized in this
document, none of them removable inside the α+d truncation): the missing 15 %
non-α+d component (§2.7); free vs in-medium cluster wave functions (§2.7);
`FitRescaled` vs `FitRaw` vs `OverlapRaw` shape choice (§8, ±9 %); the
`li6.adr.fit` truncation at 9.95 fm (§3.1, +0.6 % of ⟨R²⟩ missing relative to the
raw block's tail); the diagonal truncation of §2.6 (no effect on 𝒬, ≈ 3 % on the
angular correlation).

### 2.8 The bridge to `coherent.hpp`, and a free validation

For a transverse spin axis (x̂), a nearly-Gaussian transverse profile gives
|F|² ∝ exp(−|t|[B̄ + (Δ/2)cos2Φ]) with B̄ = (σ_x²+σ_y²)/2 = ⟨b²⟩/2 = ⟨r²⟩/3 =
`gaussian_slope`, and Δ = σ_x² − σ_y². With the axis along x̂, (G5) gives per
nucleon

    δ ≡ ⟨x²⟩ − ⟨y²⟩ = Q_matter(m) / (2A)                                       (G7)
    a₂(m) = −(δ_m/4) |t|   [arXiv:2408.13213 Eq. (9) normalization]            (G8)
    ⇒ eps_b0 ≡ δ_{±1}/B, matching `CoherentScenario::a2_m_state`
      (`coherent.cpp:110–116`, a₂(±1) = −¼ eps_b0 B |t|).

**(G8) is validated with zero free parameters against the digitized
`mantysaari_a2_deuteron()` table.** For the deuteron, δ = Q_d/2 = 3.4632 GeV⁻²:

| \|t\| [GeV²] | a₂(±1) pred | digitized | a₂(0) pred | digitized |
|---|---|---|---|---|
| 0.05 | −0.043 | −0.04 | +0.087 | +0.08 |
| 0.10 | −0.087 | −0.08 | +0.173 | +0.15 |
| 0.20 | −0.173 | −0.17 | +0.346 | +0.30 |
| 0.30 | −0.260 | −0.28 | +0.519 | +0.43 |

m = ±1 agrees to ≤ 8 % at every point; m = 0 to 8–21 %, degrading with |t|
exactly as the linear-in-|t| limit should. **A closed-form quadrupole → a₂ map
therefore reproduces the only published polarized coherent calculation in
existence.** That makes §7 T9 a real gate and gives the ⁶Li expectation before
any dipole-model run:

| Q_charge(⁶Li) assumed | δ [GeV⁻²] | a₂(±1) | a₂(±1) at \|t\| = 0.3 | eps_b0 |
|---|---|---|---|---|
| measured −0.0818 | −0.350 | +0.0875 \|t\| | **+0.026** | **−0.0067** |
| GFMC AV18+IL7 −0.20(6) | −0.856 | +0.214 \|t\| | +0.064 | −0.0165 (−0.012…−0.021) |
| α+d model −0.669 (`OverlapRaw`) | −2.864 | +0.716 \|t\| | +0.215 | −0.055 |
| α+d model −0.615 (`FitRescaled`, default) | −2.635 | +0.659 \|t\| | +0.198 | −0.0506 |
| deuteron +0.2697 (ref) | +3.463 | −0.866 \|t\| | −0.260 | +0.105 ‡ |

The `eps_b0` column is δ/B at **B(⁶Li) = 52.04 GeV⁻²** =
`gaussian_slope(√LI6_R2_POINT_FM2)` — the *measured* ⁶Li point radius, so the
number is directly comparable with `CoherentScenario::slope_b`.
**‡ The deuteron row does not use B = 52.0**: +0.105 is δ_d/B_d with
B_d = `gaussian_slope(1.967 fm)` = **33.1 GeV⁻²**, the deuteron's own point
radius, which is the only slope that makes "the deuteron's eps_b0" mean
anything. At B = 52.0 the same δ_d would read +0.067. Conclusion 3 below
compares the ⁶Li numbers *at B = 52.0* with a band whose provenance is the
deuteron *at B_d = 33.1* — that mismatch is part of why the band is
deuteron-sized, and it is the first thing O4 must settle.

Three conclusions the design hands to §11 and to `plans/06`:

1. **The ⁶Li tensor cos 2Φ has the opposite sign to the deuteron's** — Q(⁶Li)
   is negative — which is exactly the "near-null test against the deuteron's
   +0.286" that `physics_literature.md` §7 anticipated, now with a number.
2. **Its magnitude is ~10× smaller than the deuteron's** at the measured
   quadrupole (a₂(±1) = +0.026 vs −0.28 at |t| = 0.3).
3. **`CoherentScenario::eps_b0`'s band −(0.04…0.13) is a deuteron-sized band.**
   The geometric estimate for ⁶Li is −0.0067 (measured Q) to −0.055 (the α+d
   VMC overlap). The band's *magnitude* is defensible only at the top of the
   wave-function systematic; at the measured quadrupole it is **6× too large**.
   This design does **not** change `eps_b0` — that is a separate, reviewed
   decision — but it is now a documented finding with a derivation.

Also flag for that decision: `coherent.hpp:105–107`'s docstring calls `eps_b0`
"ΔB₀/B of the m = 0 state" while `a2_m_state` (`coherent.cpp:110–116`) ties it to
the m = ±1 coefficient. **How big the mismatch is depends on a ΔB convention
that is nowhere written down**, so the discrepancy cannot be stated as a bare
factor:

| reading of ΔB_m | relation to the code | mismatch |
|---|---|---|
| (A) ΔB_m is the cos 2Φ coefficient of the **slope**, B(Φ) = B̄ + ΔB_m cos 2Φ, so a₂(m) = −ΔB_m\|t\| | ΔB_{±1} = ¼ eps_b0 B, ΔB₀ = −½ eps_b0 B | eps_b0 = **−2** ΔB₀/B |
| (B) the docstring's own literal a₂(m) ≈ (ΔB_m/2)\|t\| (`coherent.hpp:124–128`) | ΔB₀ = eps_b0 B | eps_b0 = **+1** × ΔB₀/B — *no* discrepancy |
| (C) ΔB_m ≡ δ_m/2 with \|F\|² ∝ exp(−\|t\|[B̄ + (δ_m/2)cos 2Φ]) | ΔB₀ = δ₀/2 = −δ_{±1} | eps_b0 = **−1** × ΔB₀/B |

**Reviewer note / response.** The first draft asserted reading (A)'s factor −2
as if it were the only one; the review proposed (C)'s −1. Under the docstring's
*own* stated relation (B) the label and the code are already consistent, so
"the label and the code differ by a factor −2" is not established. The correct
finding is therefore weaker and more useful: **`coherent.hpp` uses the symbol
ΔB without defining it**, and O4 must fix the *convention* before touching any
factor. Use the code (`a2_m_state`) as the operative definition, and phrase the
docstring repair in the one quantity that carries no convention at all:
**eps_b0 = δ_{±1}/B with δ_{±1} = ⟨x²⟩ − ⟨y²⟩ per nucleon in the m = ±1 state**
(= Q_matter(±1)/(2A), (G7)), which is exactly what
`ClusterConfigSampler::eps_b0_equivalent()` returns.

---

## 3. Data inputs

### 3.1 Already committed (`data/vmc/`, see `data/vmc/README.md`)

| file | block/columns | units | normalization | used for |
|---|---|---|---|---|
| `li6_alpha_d/li6.ad` | r-block after the `rij` header: `r A_ad00 (err) A_ad22 (err) nsamples`, r = 0.05…19.95 fm step 0.1 | fm, fm^−3/2 | ∫(A₀²+A₂²)r²dr = S_αd = 0.85423 (**verified**, matches the file's printed `ndx/s-wave/d-wave` 0.856/0.838/0.017) | R₀, R₂ **and the relative S–D sign** |
| `li6_alpha_d/li6.adr.fit` | 3 plain columns `R R0LI6FIT R2LI6FIT`, R = 0.05…9.95 step 0.1, no errors | fm | ∫ = 0.83730, P_D = 0.02542 | smoothed R₀, R₂ (the **default**, see below) |
| `deuteron/fdeut.av18` | `r u du/dr w dw/dr` to r = 100 fm | fm | ∫(u²+w²)dr = 1 (**verified** 0.999998) | u, w; header `dstate`, `qm`, `rd` are the test anchors |

**Shape source, default = `li6.adr.fit`, normalization = `li6.ad`.** The raw
block's A₂ column changes sign four times below 1.2 fm at the MC noise level
(errors ±0.026 on a 0.014 value at r = 0.05), which a rejection/CDF sampler
cannot use; the Forest et al. smoothed fit (PRC 54 (1996) 646, Fig. 21) is
smooth and node-structured. But the fit's own norms differ (P_D 0.02542 vs
0.02011), so **rescale each fitted wave to the raw block's N₀ and N₂**:
s₀ = 1.0128, s₂ = 0.8984. That combination gives ⟨R²⟩ = 16.963 fm²,
𝒬 = −1.3204 fm², Q_matter = −1.231 fm², r_point = 2.539 fm — the "hybrid" row
of §8. Both pure variants stay reachable
(`AlphaDSource::{FitRescaled, FitRaw, OverlapRaw}`).

**Extent and extrapolation — state it or every default number is wrong.**
Every tabulated radial function used by this module is **identically zero
beyond its last abscissa**. That is `VmcRadial::operator()`'s own rule
(`cluster.cpp:232–235`: `if (k < k_.front() || k > k_.back()) return 0.0;`) and
it must be repeated in the header comment, because the module's *other* named
interpolation primitive does the opposite: **`np_interp` CLAMPS to the end
value** (`numerics.hpp:13–15`, "including its clamped extrapolation"), and the
two defaults collide. `AlphaDSource::FitRescaled` reads `li6.adr.fit`, which
ends at **R = 9.95 fm with R₀ = −0.0039**, while `r_max_fm` defaults to
**20 fm**. A straightforward `np_interp` implementation therefore carries
R₀ = −0.0039 flat from 9.95 to 20 fm, and — recomputed in this pass — that
alone moves ⟨R²⟩ from **16.97 to 27.39 fm²** and 𝒬 from **−1.320 to −1.104
fm²**, i.e. *every* default number in §8 and §2.8. Rules:

* tables are zero outside their own abscissa range; `np_interp` may be called
  **only inside** the table;
* the working grid runs over [0, min(`r_max_fm`, last abscissa)] — `r_max_fm`
  and `rnp_max_fm` are ceilings, never extensions;
* the §8 `Fit*` rows are for the fit **truncated at 9.95 fm**, and that
  truncation is a listed systematic: the raw `li6.ad` block's R > 10 fm tail
  carries **0.087 % of the norm and 0.106 fm² of ⟨R²⟩ (0.60 %)**, which the fit
  does not have.

**Sign trap.** `vmc_from_momentum`'s k-space node/anchor machinery
(`cluster.cpp:311`) must **not** be reused here. The k-space and r-space sign
structures are different functions (they are Fourier–Bessel transforms of one
another): in k the ratio ψ₂/ψ₀ is −1 below 0.678 fm⁻¹, +1 to 2.25 fm⁻¹, −1
above (`vmc_reconciliation.md`); in r, R₀ has its single node at 1.85 fm and R₂
a small-r sign change at ≈ 1.1 fm, so R₀R₂ < 0 through the whole dominant
2–8 fm region. **Take the r-space signs from the r-space columns, verbatim.**
The global phase is unobservable and must be fixed once, by convention, as
R₀(r → ∞) > 0 (i.e. negate both raw columns); the tests must show 𝒬 is
invariant under that flip.

*And the second half of the trap, which the first draft left open:* "Fourier–
Bessel transforms of one another" is true only for the **no-i^L** transform. The
physical relation carries (−i)^L (§2.1), so the *physical* k-space S–D
interference has the **opposite** sign to the ANL k-block's, to `VmcRadial`'s,
and to `vmc_reconciliation.md`'s. Two consequences, both testable:

* This module works in r space and applies **no phase**. Its own sign
  convention is validated against an **observable**, not against the k-block:
  `sign(q_matter_analytic_fm2(±1)) == sign(LI6_QUADRUPOLE_FM2)` (both negative
  — oblate m = ±1). §7 T4(i).
* The r-tables are tied to the existing k-space validation by the **no-phase**
  Fourier–Bessel transform of (R₀, R₂), which must reproduce `li6.ad`'s k-block
  sign structure without asserting sign identity with the physical amplitude.
  Verified here on the `FitRescaled` tables: ψ̃₀ crosses zero at **0.678 fm⁻¹**
  (the k-block's own S node) and ψ̃₂/ψ̃₀ = **−0.0006, −0.0115, −0.069, −0.160**
  at k = 0.05, 0.2, 0.4, 0.5 fm⁻¹ — negative below the node, positive above
  (+0.21 at 1.0 fm⁻¹), exactly the k-block's structure. §7 T4(ii). *Physically*
  those ratios are +0.0006…+0.160, because of the (−i)².

### 3.2 To be fetched — the ⁴He one-body density

`data/vmc/bonus_other_clusters/he4.dd` and `he4.tp` are **two-cluster overlaps**
(⁴He → d + d and ⁴He → t + p). They are *not* a one-body density and cannot be
used for the α core: `he4.dd` gives the d–d relative motion of a *different*
decomposition and would double-count. Nothing currently in `data/` is a ⁴He
nucleon density.

**Fetch (verified reachable and plain text during this design pass, 2026-09-02;
the ANL site itself 403s behind Cloudflare, so use the Wayback `id_` raw
form exactly as `data/vmc/README.md` documents):**

    curl -o data/vmc/density/he4.density \
      "http://web.archive.org/web/20150905165459id_/http://www.phy.anl.gov/theory/research/density/he4.density"
    curl -o data/vmc/density/li6.density \
      "http://web.archive.org/web/20150905180021id_/http://www.phy.anl.gov/theory/research/density/li6.density"

| file | bytes | md5 | content |
|---|---|---|---|
| `he4.density` | 4751 | `5d717cf0ed642dc7779aef3873e06317` | `⁴He(0⁺) — AV18+UX — 18-Dec-12`, VMC 500k samples; `R RHORP DRHORP`, R = 0.05…20.05 fm step 0.1; header prints `4π∫ρR²dR = 1.9974` (= Z) and rms 1.4404 fm |
| `li6.density` | 5189 | `96c6b15f045785db0593c00e9a32f9b1` | `⁶Li(1⁺) — AV18+UX — 20-Mar-14`, VMC 200k; same columns; `4π∫ρR²dR = 2.9991`, rms **2.4433 fm** |

`li6.density` is **not an input** — it is the independent validation target for
the assembled ⁶Li one-body density (§7 T6). Commit both (≈ 10 kB) with a
`data/vmc/density/README.md` row in the same table format as
`data/vmc/README.md`, recording the live URL, the Wayback timestamp, the fetch
date and the md5. Optional third file for a Hamiltonian systematic:
`he4_v18.density` (AV18 without UX, snapshot `20170427210207`).

These are **point-nucleon** densities (ANL ρ₁, no nucleon form factor folded
in) — exactly what a dipole model wants — and for N = Z = 2 the neutron density
equals the proton density by isospin, so ρ_nucleon = ρ_RHORP/Z.

### 3.3 Parsers

* `he4.density` parses correctly with the **existing** `read_anl_momentum`
  (`cluster.cpp:159`): its `****` block rule matches, rows have 3 numbers →
  ncol = 1, `col[0]` = RHORP, `err[0]` = DRHORP. Reuse it; add nothing.
* `li6.adr.fit` **must not** go through `read_anl_momentum`: its rows also have
  3 numbers, so the reader would silently land `R2LI6FIT` in `err[0]`. Add one
  small sibling to the ANL readers in `cluster.hpp`/`cluster.cpp`, next to
  `read_anl_overlap`, purely additive:

      /// A plain `x v1 v2 ... vN` table with no MC-error columns, introduced
      /// by the `****` column rule (`li6.adr.fit`, `*.rho1` fits).
      AnlTable read_anl_plain(const std::string& path, std::size_t ncol);

  This is the only edit outside the two new files, and it changes no existing
  behaviour. If even that is unwanted, make it a file-local static in
  `cluster_config.cpp` — but then say so in the header comment, because it
  belongs with the other ANL readers by convention.

---

## 4. API sketch — `include/lipolgen/cluster_config.hpp`

```cpp
#ifndef LIPOLGEN_CLUSTER_CONFIG_HPP
#define LIPOLGEN_CLUSTER_CONFIG_HPP

/// \file cluster_config.hpp
/// Nucleon-position CONFIGURATIONS of a polarized 6Li in the alpha + d cluster
/// picture -- the input a Good-Walker dipole-model code (subnucleondiffraction,
/// arXiv:2408.13213) averages the coherent amplitude over.
///
/// OPT-IN AND INERT.  Nothing in the generator calls this header.  It adds no
/// default, changes no existing number, and its output is a file.  See
/// docs/open_items/run_2026-09-02/design_G_cluster_config.md for the physics,
/// the validations and the (large) wave-function systematic on the tensor
/// observable.

#include <array>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "lipolgen/cluster.hpp"
#include "lipolgen/rng.hpp"

namespace lipolgen {

// --------------------------------------------------------------- constants

/// RNG `bunch` slot owned by this module.  The generator uses the spin
/// CATEGORY INDEX (0, 1, 2, ...) as its bunch slot (generator.hpp:125), so a
/// configuration stream must sit far away from small integers.  Bound to
/// Python as `CONFIG_STREAM`; the sidecar's "rng":{"stream"} records it.
inline constexpr std::uint64_t kConfigStream = 0x434F4E464947ull;  // 'CONFIG'

/// r_point(6Li) from the ANL VMC one-body density li6.density (AV18+UX),
/// header "SQRT(4*PI*TOTINT(RHORP*R**4:R)/3) = 2.4433".  The DEFAULT target of
/// match_li6_radius(); a MODEL number, not a measurement (the measurement is
/// rc.hpp's LI6_R2_POINT_FM2 = 6.0788 fm^2 -> 2.4655 fm).
inline constexpr double LI6_R_POINT_VMC_FM = 2.4433;

/// Q(6Li) from GFMC with AV18+IL7: -0.20(6) e fm^2.  Pastore, Pieper,
/// Schiavilla & Wiringa, PRC 87, 035503 (2013), arXiv:1212.3375; the repo's
/// own record is docs/open_items/physics_literature.md:258.  NOT VMC and NOT
/// the AV18+UX Hamiltonian of li6.density / li6.ad, so it is an
/// independent-Hamiltonian comparison point, not a check on our own tables.
/// The middle entry of quadrupole_band_fm2(); nowhere else is it retyped.
inline constexpr double LI6_QUADRUPOLE_GFMC_FM2     = -0.20;
inline constexpr double LI6_QUADRUPOLE_GFMC_ERR_FM2 = 0.06;

// ------------------------------------------------------------------ sources

/// Where the alpha core's one-body point-nucleon density comes from.
enum class AlphaCoreSource {
  /// data/vmc/density/he4.density (ANL VMC AV18+UX, 500k samples).  DEFAULT.
  VmcHe4Density = 0,
  /// Gaussian rho ~ exp(-r^2/2a^2) with a^2 = cluster_point_a2_fm2(2, 4)
  /// = 0.7001 fm^2 -- the SAME single copy the FSI weight uses (fsi.hpp).
  /// rms 1.4492 fm vs the table's 1.4413: a 0.5 % size difference and a
  /// visibly different surface (the VMC alpha has a flat-topped interior).
  Gaussian = 1,
};

/// Where the alpha-d relative radial functions R_0(R), R_2(R) come from.
enum class AlphaDSource {
  /// li6.adr.fit shapes, each wave rescaled to li6.ad's own N_0, N_2.  DEFAULT.
  FitRescaled = 0,
  FitRaw      = 1,   ///< li6.adr.fit as published (P_D = 0.02542)
  OverlapRaw  = 2,   ///< li6.ad r-block verbatim -- noisy below ~1.2 fm
};

// ----------------------------------------------------------------- options

struct ClusterConfigOptions {
  AlphaCoreSource alpha_source = AlphaCoreSource::VmcHe4Density;
  AlphaDSource    alpha_d_source = AlphaDSource::FitRescaled;

  /// Quantization axis in the ion rest frame [rad]; the SAME axis as
  /// spin.hpp / tagged.hpp.  A cos 2Phi signal needs theta_s != 0.
  double theta_s = 1.5707963267948966;   ///< default: transverse (pi/2)
  double phi_s   = 0.0;

  /// Inflate the alpha's SOURCE density by lambda = (1 - 1/4)^-1/2 = 1.154701
  /// in radius before recentring.  What this does and does NOT do:
  ///   * <s^2> closes EXACTLY, for ANY source shape -- recentring four
  ///     independent draws gives <s^2> = (3/4)<v^2>, so lambda^2 = 4/3 puts
  ///     the recentred second moment back on the table's.  Measured: 2.0784
  ///     against he4.density's 2.0771 fm^2.
  ///   * the SHAPE is NOT restored for a tabulated density.  Measured against
  ///     he4.density: +73 % at r = 0.05 fm, +25 % at 0.25, -7 % at 0.85-1.15,
  ///     +5..10 % at 1.5-2.5, -15 % at 3.05 fm; chi^2/ndf ~ 170 on 40 bins
  ///     with 10^6 entries.  It IS exact in shape for AlphaCoreSource::Gaussian.
  /// The recentred density has a closed form the tests use instead
  /// (sec. 4.1 step 4): its 3-D transform is ftilde(3q/4) * ftilde(q/4)^3.
  bool alpha_cm_inflate = true;

  /// Minimum nucleon-nucleon separation inside the alpha core [fm].  0 = off
  /// (the default: an uncorrelated product of one-body densities).  0.9 is
  /// what arXiv:2605.00454 imposes on its Woods-Saxon sampling.
  double min_nn_separation_fm = 0.0;

  /// Multiply the alpha-d separation R by this before assembling.  1.0 = the
  /// wave functions as published (DEFAULT, no tuning).  `match_li6_radius()`
  /// returns the value that puts <r^2> on LI6_R_POINT_VMC_FM (li6.density).
  double alpha_d_scale = 1.0;

  /// If non-zero, scale the alpha-d D-wave amplitude by the ROOT s of (G9)
  /// (sec. 2.7) so that (G5) returns this Q_charge [fm^2].  NOT
  /// sqrt(target/model): Q is LINEAR in the D amplitude through the
  /// interference term (95 % of it) and the deuteron term is an offset, so the
  /// sqrt rule misses by 7x (it gives -0.048 for a -0.0818 request).  A
  /// DEFORMATION DIAL, not a wave function: the root for the measured Q is
  /// s = 0.4032, P_D(alpha-d) = 3.3e-3, and the writer stamps it in the header.
  /// REACHABLE RANGE, s in [0, 1]: [q_matter_analytic_fm2(1)/2, +Q_d], i.e.
  /// [-0.615, +0.270] fm^2 for the default source -- s = 0 leaves the
  /// deuteron's own +0.270 fm^2 behind.  validate() THROWS outside it.
  /// 0 = off (DEFAULT).
  double quadrupole_target_fm2 = 0.0;

  /// Reserved.  The exact 5-D joint density with the m_S coherences of sec.
  /// 2.6 restored.  Construction THROWS if set; the flag exists so the
  /// approximation is visible in the API.
  bool exact_coherence = false;

  std::size_t n_r = 512;    ///< radial grid cells of each (x, cos theta) table
  std::size_t n_c = 96;     ///< cos-theta cells (matches TaggedModel's default)
  double r_max_fm = 20.0;   ///< alpha-d and alpha-core grid extent
  double rnp_max_fm = 25.0; ///< p-n grid extent

  void validate() const;    ///< throws std::runtime_error on a bad combination
};

// ---------------------------------------------------------------- the data

/// One nucleon of one configuration.  Positions are in the ION REST FRAME,
/// c.m. at the origin, in fm.
struct NucleonPos {
  double x = 0.0, y = 0.0, z = 0.0;
  int isospin = +1;    ///< +1 proton, -1 neutron
  int cluster = 0;     ///< 0 = alpha core, 1 = deuteron
};

/// One 6Li configuration.
struct ClusterConfig {
  std::array<NucleonPos, 6> nucleon{};
  int m_ion = 1;       ///< the substate this configuration was drawn for
  double m_s = 1.0;    ///< the deuteron projection drawn (diagnostic)
  double r_ad = 0.0;   ///< |R| [fm]      (diagnostic)
  double r_np = 0.0;   ///< |r_p - r_n|   (diagnostic)
};

/// N configurations plus everything needed to reproduce and label them.
struct ClusterConfigSet {
  std::vector<ClusterConfig> config;
  int m_ion = 1;
  double theta_s = 0.0, phi_s = 0.0;
  std::uint64_t seed = 0, run = 0;
  std::string provenance;    ///< every input file + option, one line each
  /// Moments of the SET, computed on the fly (sec. 7 uses these directly).
  double r2_mean_fm2 = 0.0;      ///< <r^2> per nucleon
  double q_matter_fm2 = 0.0;     ///< sum_i <3 z_i^2 - r_i^2> about the axis
  double delta_perp_fm2 = 0.0;   ///< <x^2> - <y^2> per nucleon, axis = x
  /// PER-CONFIGURATION sample standard deviations of the same three moments.
  /// They are what makes every MC gate in sec. 7 self-calibrating: the gate is
  /// k * <sd> / sqrt(n), never a hand-typed absolute number.  Measured at the
  /// default source: 3.66 fm^2, 28.9 fm^2, 2.87 fm^2 -- so at n = 1e5 one
  /// sigma is 0.012, 0.091 and 0.009 fm^2 respectively.
  double r2_sd_fm2 = 0.0, q_matter_sd_fm2 = 0.0, delta_perp_sd_fm2 = 0.0;
  std::array<double, 3> cm{{0.0, 0.0, 0.0}};   ///< must be 0 to 1e-12
};

// -------------------------------------------------------------- the sampler

/// Builds its three (value, cos theta) grids ONCE in the constructor and is
/// immutable afterwards -- the same discipline as TaggedModel (tagged.hpp:205),
/// so it is thread-safe with no lock and every accessor is a pure lookup.
class ClusterConfigSampler {
 public:
  explicit ClusterConfigSampler(ClusterConfigOptions opt = {});

  const ClusterConfigOptions& options() const;

  /// One configuration for ion substate `m_ion` in {+1, 0, -1}.  `rng` is the
  /// caller's; the sampler never owns or seeds one.
  ClusterConfig sample(Rng& rng, int m_ion) const;

  /// N configurations, each drawn from its OWN counter-based stream
  /// Rng(seed, run, CONFIG_STREAM, i), so config i is bit-identical whatever
  /// the order, the thread count, or N (docs/CONVENTIONS.md).
  ClusterConfigSet sample_set(std::size_t n, int m_ion, std::uint64_t seed,
                              std::uint64_t run = 0) const;
  /// The unpolarized set: equal thirds of m = +1, 0, -1, interleaved.
  ClusterConfigSet sample_set_unpolarized(std::size_t n, std::uint64_t seed,
                                          std::uint64_t run = 0) const;

  // ---- analytic predictions, for the tests and for the header ------------
  double p_d_alpha_d() const;              ///< P_D of the alpha-d relative wave
  /// D_T, eq. (G5).  The CLOSED FORM 1 - 0.9 * p_d_alpha_d(); TaggedModel's
  /// same-named quantity is a (k, cos theta) grid quadrature at a different
  /// default P_D, so the two agree at 1e-4, not exactly (sec. 2.4, T10).
  double tensor_dilution() const;
  double q_matter_analytic_fm2(int m) const;   ///< (G5)
  /// The (G4) split sec. 2.5 quotes: interference 2 sqrt2 f0 f2 and pure-D
  /// -f2^2 parts of Q[R_0, R_2].  Needed by the (G9) root-find and by T4.
  double q_int_fm2() const;
  double q_dd_fm2() const;
  double r2_analytic_fm2() const;              ///< (G6)
  double delta_perp_analytic_fm2(int m) const; ///< (G7); 0 when theta_s == 0
  /// Asymptotic D/S ratio eta = C_2/C_0 of the alpha-d channel: R_2/R_0
  /// divided by the Whittaker ratio W_{-eta_c,5/2}/W_{-eta_c,1/2}(2 kappa R),
  /// averaged over R = 6..8 fm (kappa = 0.3074 fm^-1, eta_c = 0.3002).
  /// -0.05(1) for the committed tables, against the MEASURED -0.025(12)
  /// (George & Knutson, PRC 59, 598 (1999)).  Sec. 2.7, O1.
  double asymptotic_ds_ratio() const;
  /// a_2(m) at |t| for THIS sampler's 6Li geometry -- a2_from_quadrupole()
  /// with q_matter_analytic_fm2(m) and A = 6.
  double a2_from_geometry(double t_abs, int m) const;
  /// eps_b0 equivalent: delta_perp(+1) [GeV^-2] / gaussian_slope(r_point)
  /// with r_point = sqrt(LI6_R2_POINT_FM2) (rc.hpp), B = 52.04 GeV^-2.  The
  /// MEASURED point radius is used, not the model's own 2.539 fm (B = 55.2),
  /// so the number is directly comparable with CoherentScenario::slope_b's 50
  /// and its eps_b0 band.  Every fm <-> GeV step goes through HBARC_GEV_FM.
  double eps_b0_equivalent() const;
  /// alpha_d_scale that puts <r^2> on `target_rms_fm` (default: li6.density).
  double match_li6_radius(double target_rms_fm = LI6_R_POINT_VMC_FM) const;
  /// {LI6_QUADRUPOLE_FM2 (measured), LI6_QUADRUPOLE_GFMC_FM2 (GFMC AV18+IL7),
  /// THIS SOURCE's q_matter_analytic_fm2(1)/2} Q_charge [fm^2] -- the
  /// mandatory band of sec. 2.7.  Every entry is a named constant or a
  /// computed number; the third one moves with alpha_d_scale /
  /// quadrupole_target_fm2 and is never a literal.
  std::array<double, 3> quadrupole_band_fm2() const;

  /// The three densities, exposed for the tests and for plots.
  double rho_alpha(double r_fm) const;                     ///< per nucleon
  /// (G1') -- the alpha-d density CONDITIONED on the deuteron projection m_s,
  /// |A_{m_s}(m; R, c)|^2.  This, not the m_s-summed one, is what sec. 4.1
  /// step 2 draws from (sec. 2.6); the signature mirrors
  /// TaggedModel::sample_kc(m_ion, m_s, ...) for the same reason.
  double rho_alpha_d(double R_fm, double c, int m, double m_s) const;
  /// Convenience: the m_s-SUMMED density (G1).  For plots and for the
  /// marginal checks only -- sampling from it decorrelates R-hat from m_s.
  double rho_alpha_d_summed(double R_fm, double c, int m) const;
  double rho_np(double r_fm, double c, double m_s) const;  ///< (G1)/(G2)/(G3)

 private:
  /* grids, CDFs, tables, options */
};

// ------------------------------------------------ the quadrupole -> a_2 map

/// a_2(m) at |t| [GeV^2] from a point-matter quadrupole, (G7) + (G8):
///   delta = q_matter_fm2 / (2 a)   [fm^2, per nucleon, axis transverse]
///         -> GeV^-2 through HBARC_GEV_FM,
///   a_2(m) = -(delta_m / 4) |t|,  with delta_0 = -2 delta_{+-1}.
/// A FREE FUNCTION so the deuteron can be fed through it: sec. 7 T9 calls it
/// with (2 * 0.269670, a = 2) and reproduces mantysaari_a2_deuteron() to 8 %
/// at m = +-1 and 8-21 % at m = 0.  The member a2_from_geometry() is this
/// function called with q_matter_analytic_fm2(m) and a = 6.
double a2_from_quadrupole(double q_matter_fm2, int a, double t_abs, int m);

// ------------------------------------------------------------- the writer

/// Which flavour of the subnucleondiffraction configuration table to write.
enum class SndConfigFormat {
  /// BYTE-COMPATIBLE with he3.dat: numeric rows only, NO comments.  The
  /// upstream reader (src/nucleons.cpp:118-152) does a bare `ss >> x` per
  /// field and would silently read garbage from a '#' line.  Metadata goes to
  /// the sidecar.  DEFAULT.
  He3Compatible = 0,
  /// The same rows with a '#'-prefixed header block.  Requires the patched
  /// reader of sec. 5.3.
  Annotated = 1,
};

/// Writes `set` as `path` plus the sidecar `path + ".meta.json"`, which always
/// carries m, (theta_s, phi_s), seed/run, EVERY option field, every input file
/// with its size and its own printed normalization, LIPOLGEN_VERSION, and the
/// analytic and sampled moments.  Returns the number of rows written.
///
/// NO md5 AND NO GIT SHA FROM C++.  There is no hash function anywhere in
/// src/core or include/lipolgen, no JSON emitter in the core (tests/json_min.hpp
/// is a READER), and CMake compiles in LIPOLGEN_VERSION only -- no git SHA, and
/// none is definable in a wheel install.  Adding an MD5 implementation to the
/// core library for a metadata field is out of scope for an opt-in, inert
/// module.  The C++ writer emits a flat, hand-rolled JSON object (no library);
/// python/lipolgen/configs.py adds "md5" via hashlib and "git" via
/// `git rev-parse HEAD` when it is inside a checkout, else null.  The reference
/// md5s of the input tables live in data/vmc/density/README.md (sec. 3.2).
std::size_t write_snd_configs(const ClusterConfigSet& set,
                              const std::string& path,
                              SndConfigFormat fmt = SndConfigFormat::He3Compatible);

}  // namespace lipolgen
#endif
```

### 4.1 Sampling algorithm (per configuration)

1. `m_S` from the tabulated P(m_S | M) (3-way inverse CDF; § 2.4).
   P(+1) = N₀ + N₂/10, P(0) = 3N₂/10, P(−1) = 6N₂/10.
2. (R, cosθ_R) from the **m_S-conditioned** table
   |A_{m_S}(M; R, cosθ_R)|² R² on the (n_r × n_c) grid by cell inverse-CDF +
   uniform inside the cell — **the same primitive, and the same per-m_S
   structure, as `TaggedModel::build_amp2` / `sample_kc(m_ion, m_s, …)`**
   (`tagged.hpp:185, 197`), which is measured at 5.2 M draws/s. Uniform φ_R.
   There are 3 (m_S) × 3 (M) such tables, built once in the constructor.
   **Not** the m_S-summed density: see §2.6 and T21. (Equivalently: draw the
   triple (m_S, R, cosθ_R) from one stacked CDF; steps 1 and 2 are then one
   inverse-CDF lookup and the m_S marginal is automatically right.)
3. (r, cosθ_r) from ρ_{np}(r, cosθ_r; m_S) r² the same way. Uniform φ_r.
4. Four α nucleons: r from the 1-D inverse CDF of 4πρ_α(r)r² (inflated by
   λ = 1.154701 if `alpha_cm_inflate`), isotropic direction; if
   `min_nn_separation_fm > 0`, redraw any nucleon violating it (cap the
   retries and throw on exhaustion). Then subtract their mean so Σ s⃗_i = 0
   **exactly** (this is the CM constraint, applied to the core only — the α's
   c.m. is then placed at −R⃗/3 and the whole 6-body c.m. is 0 by construction,
   which the test checks at 1e-12).
   **The recentred one-body density has a closed form and the tests use it.**
   With s⃗₁ = (3/4)v⃗₁ − (1/4)(v⃗₂+v⃗₃+v⃗₄) and v⃗ isotropic with 3-D transform
   f̃(q) = ⟨j₀(q|v⃗|)⟩, the recentred single-nucleon density has transform
   **f̃_s(q) = f̃(3q/4) · f̃(q/4)³**, and ⟨s²⟩ = (3/4)⟨v²⟩ *exactly, for any
   shape* — which is the whole content of `alpha_cm_inflate`. Numerically
   inverting f̃_s is a one-dimensional deterministic transform, and it is the
   reference T5 gates the sampled histogram against. Recentring does **not**
   reproduce `he4.density` itself (§7 T5).
5. Assemble §2.1's six positions; assign isospin (p, n in the deuteron;
   2p + 2n in the α).
6. Rotate by R_z(φ_S) R_y(θ_S).

RNG discipline: one `Rng` per configuration, `Rng(seed, run, kConfigStream, i)`
with `kConfigStream = 0x434F4E464947` ('CONFIG'), declared in the header (§4) and
bound as `CONFIG_STREAM`. It must be far from small integers because the
generator uses the **spin category index** (0, 1, 2, …) as its `bunch` slot
(`generator.hpp:125`). Every draw inside a configuration comes from that one
object in a fixed order, so a configuration is reproducible independently of N
and of thread count. **No global RNG, no `std::rand`, no GSL.**

---

## 5. Output format

### 5.1 What the consumer actually reads (determined from the code, not guessed)

Cloned `github.com/hejajama/subnucleondiffraction` at master
(`1db86bb866b2f7cdf7556f5849fb62032c01717e`) during this pass:

* `src/nucleons.cpp:112–155`, the A == 3 branch, opens the **hard-coded
  relative path `"he3.dat"`** in the process CWD, walks lines with a 0-based
  index, and on `index == he3_id` reads **nine whitespace-separated doubles**
  `x1 y1 z1 x2 y2 z2 x3 y3 z3`, then builds `Vec n1(x1*FMGEV, y1*FMGEV)` —
  a **two-argument** `Vec`, so **z is read and then discarded** (the eikonal
  amplitude only needs transverse positions).
* `FMGEV = 5.068` (`src/subnucleon_config.hpp:20`), applied as
  `position * FMGEV`, so **the file is in fm** and the code converts to GeV⁻¹.
* `Nucleons::SetHeId` (`src/nucleons.cpp:407`) hard-bounds the id to
  **0…13698**; the shipped `he3.dat` has 13700 lines.
* The shipped `he3.dat` rows carry **13** fields: the 9 coordinates, then three
  numbers that sum to 1.000 (0.4522 0.3592 0.1886 …) — evidently per-nucleon
  weights — and a trailing integer in {1,2,3}, plausibly the sampled neutron
  index. **The reader ignores fields 10–13 entirely.** That is the slot our
  isospin column goes in.
* There is **no header, no comment syntax, and no tolerance for one**: a `#`
  line would leave the nine doubles at their previous/indeterminate values and
  produce a silently wrong configuration.
* All nucleons are the same `Ipsat_Proton` dipole object — the code has **no
  proton/neutron distinction at all**.
* `-He3 <id>` is the CLI hook (`src/main.cpp:333`), applied at
  `src/main.cpp:443–446` under `if (A == 3)`.
* **A = 6 currently falls to the Woods–Saxon `else` branch**
  (`src/nucleons.cpp:156–170`). Reading a ⁶Li table needs a code change
  upstream — see §5.3.
* **No branch of the public repository contains any polarization code.** All
  20 remote heads were fetched and grepped: `polariz` appears nowhere in
  `nucleons.{cpp,hpp}` or `main.cpp` on any of them. The arXiv:2408.13213
  polarized-deuteron implementation is **unreleased**. This is the single most
  important fact for the collaboration ask (§10).

### 5.2 The format we write

**`SndConfigFormat::He3Compatible` (default)** — one line per configuration,
space-separated, no header, no blank lines, LF endings, `%.17g`:

    x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4 x5 y5 z5 x6 y6 z6  t1 t2 t3 t4 t5 t6  m

18 coordinates in **fm**, ion rest frame, c.m. at the origin, quantization axis
already rotated to (θ_S, φ_S); then six isospins (+1 p, −1 n; ordering is
α, α, α, α, p, n — nucleons 1–4 are the core); then the ion substate m. Fields
19–25 are past the 3A the reader consumes and are therefore invisible to it,
exactly as `he3.dat`'s own trailing fields are.

**The sidecar `<path>.meta.json` is mandatory** and is where "the polarization
axis and m-state recorded in the header" actually lives. **Written in two
passes** (§4, `write_snd_configs`): the C++ writer emits everything it can
compute from the core — including each input file's size and its *own printed
normalization*, which the existing readers already parse — and the Python CLI
adds the two fields the core has no business implementing, `md5` (hashlib) and
`git` (`git rev-parse HEAD`, `null` outside a checkout). **All 13
`ClusterConfigOptions` fields appear under `options`** (T16 checks the list),
with `axis` kept as a duplicate, human-facing view of `theta_s`/`phi_s`:

```json
{"generator":"LiPolGen cluster_config","version":"0.1.0",
 "git":null,
 "nucleus":"6Li","A":6,"Z":3,"n_config":100000,"m_ion":1,
 "axis":{"theta_s":1.5707963267948966,"phi_s":0.0,"frame":"ion rest frame"},
 "units":"fm","column_order":["x","y","z"],"cm_convention":"sum m_i r_i = 0",
 "rng":{"seed":20260902,"run":0,"stream":"CONFIG","stream_value":74007894706503},
 "inputs":[{"file":"vmc/li6_alpha_d/li6.adr.fit","bytes":3400,
            "printed_norm":{"integral":0.83730,"p_d":0.02542},"md5":null},
           {"file":"vmc/li6_alpha_d/li6.ad","bytes":36000,
            "printed_norm":{"ndx":0.856,"s":0.838,"d":0.017},"md5":null},
           {"file":"vmc/deuteron/fdeut.av18","bytes":1100000,
            "printed_norm":{"dstate":0.057599,"qm":0.269673,"rd":1.967364},
            "md5":null},
           {"file":"vmc/density/he4.density","bytes":4751,
            "printed_norm":{"norm_4pi":1.9974,"rms":1.4404},"md5":null}],
 "options":{"alpha_source":"VmcHe4Density","alpha_d_source":"FitRescaled",
            "theta_s":1.5707963267948966,"phi_s":0.0,
            "alpha_cm_inflate":true,"min_nn_separation_fm":0.0,
            "alpha_d_scale":1.0,"quadrupole_target_fm2":0.0,
            "quadrupole_dial_s":1.0,
            "exact_coherence":false,
            "n_r":512,"n_c":96,"r_max_fm":20.0,"rnp_max_fm":25.0},
 "moments":{"r2_mean_fm2":6.4447,"r_rms_fm":2.5386,"q_matter_fm2":-1.2309,
            "delta_perp_fm2":-0.1026,"eps_b0_equivalent":-0.0506,
            "p_d_alpha_d":0.02011,"tensor_dilution":0.98190,
            "asymptotic_ds_ratio":-0.05,
            "sampled":{"r2_mean_fm2":6.4436,"r2_sd_fm2":3.66,
                       "q_matter_fm2":-1.204,"q_matter_sd_fm2":28.9,
                       "delta_perp_fm2":-0.1027,"delta_perp_sd_fm2":2.87}},
 "quadrupole_band_fm2":{"measured":-0.0818,"gfmc_av18_il7":-0.20,
                        "alpha_d_model":-0.6154},
 "caveats":["Q_charge(model) = -0.615 fm^2 (range -0.615..-0.730) vs measured",
            "-0.0818 (TUNL) and GFMC AV18+IL7 -0.20(6) (Pastore et al.,",
            "PRC 87, 035503 (2013), arXiv:1212.3375);",
            "see docs/open_items/run_2026-09-02/design_G_cluster_config.md sec. 2.7"]}
```

`"md5": null` and `"git": null` are what the **C++** writer emits; the CLI fills
them in. `quadrupole_dial_s` is the root of (G9) actually applied (1.0 when the
dial is off) — the "labelled as a deformation dial" requirement of §2.7.

**`SndConfigFormat::Annotated`** writes the same rows with the JSON as `#`
comment lines on top. Use it only with the patched reader.

### 5.3 The upstream patch the graft needs (≈ 30 lines)

    // nucleons.hpp
    void SetConfigFile(std::string f);   // path, default "he3.dat"
    void SetConfigId(int i);             // replaces SetHeId, no hard bound

    // nucleons.cpp: replace the `A==3` special case by a general branch taken
    // whenever a config file was set:
    //   read line `config_id`, parse 3*A doubles, build Vec(x*FMGEV, y*FMGEV),
    //   error out if fewer than 3*A numeric fields were consumed.
    // main.cpp: -configfile <path> and -configid <n>, applied for any A.

That single generalization makes the code A-agnostic and removes the 13698
bound, the hard-coded filename and the `A == 3` gate all at once. Nothing else
in the code cares about A.

---

## 6. Wiring

* **CMake.** `src/core/cluster_config.cpp` joins `LIPOLGEN_CORE_SOURCES`
  (`CMakeLists.txt:59` — confirmed a `file(GLOB … CONFIGURE_DEPENDS
  src/core/*.cpp)`, so nothing to edit for the source itself). No new option, no
  new target, no new dependency. **One line is needed, though**, and
  `CMakeLists.txt` therefore joins §9's edited-files list: in-tree runs get
  `build/lipolgen-run` from a `file(GENERATE …)` shim
  (`CMakeLists.txt:217–226`), not from `[project.scripts]`, so
  `lipolgen-configs` does not exist in a build tree unless the sibling shim is
  added next to it —

      file(GENERATE OUTPUT ${CMAKE_BINARY_DIR}/lipolgen-configs CONTENT
      "#!/bin/sh
      # generated by CMake -- LiPolGen configuration-table writer
      PYTHONPATH=\"${CMAKE_BINARY_DIR}/python:\${PYTHONPATH}\" \\
        exec \"${Python3_EXECUTABLE}\" -m lipolgen.configs \"$@\"
      " FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
        GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)

  (copy the exact quoting and the flush-left `CONTENT` body from the existing
  `lipolgen-run` shim at `CMakeLists.txt:217–226`; the indentation above is the
  document's, not the file's.)

* **Data install.** `data/vmc/density/` is inside the already-installed
  `data/vmc` tree (`OPEN_ITEMS_SOLUTIONS.md` §12), so packaging needs nothing.
* **Bindings** (`python/bindings.cpp`). The file is organised as **one
  `static void bind_<module>(py::module_&)` per header**, forward-declared at
  lines 530–542 and called at 601–613 — *not* as blocks appended inside another
  binder. So: add `static void bind_cluster_config(py::module_& m);` to the
  forward list, `bind_cluster_config(m);` to the module body, and the definition
  next to `bind_tagged`'s. It binds `AlphaCoreSource`, `AlphaDSource`,
  `SndConfigFormat`, `ClusterConfigOptions` (all 13 fields),
  `NucleonPos`, `ClusterConfig`, `ClusterConfigSet` (with `positions()`
  returning an (N, 6, 3) numpy array and `isospin()` an (N, 6) int array — the
  tests want arrays, not object lists), `ClusterConfigSampler` with every
  analytic predictor of §4, the free `a2_from_quadrupole`, `write_snd_configs`,
  and the constants `CONFIG_STREAM`, `LI6_R_POINT_VMC_FM`,
  `LI6_QUADRUPOLE_GFMC_FM2`/`_ERR_FM2` (plus `LI6_R2_POINT_FM2` from
  `rc.hpp`). **`positions()` needs a new 3-D copy
  helper**: the existing `move_array2` (`bindings.cpp:102`) is 2-D only, so add
  a `py::array_t<double>({n, 6, 3})` variant beside it. Per `PLAN.md`'s ground
  rules every new option and constant gets a binding and a pytest.
* **CLI.** A **new console entry point, `lipolgen-configs`**, not a flag on
  `lipolgen-run`: `python/lipolgen/configs.py` with `main()`, registered as
  `lipolgen-configs = "lipolgen.configs:main"` under `[project.scripts]`
  alongside the existing `lipolgen-run` (`pyproject.toml:30–33`), plus
  `python -m lipolgen.configs`. `lipolgen-run`'s flat argparse is untouched —
  adding subparsers there would change its parse behaviour, which the
  "no default changes" rule forbids.

      lipolgen-configs --m {+1,0,-1,unpolarized} --n 100000 \
          --theta-s 1.5707963 --phi-s 0 --seed 20260902 --run 0 \
          --alpha-source {vmc,gaussian} --alpha-d-source {fit-rescaled,fit-raw,overlap-raw} \
          --alpha-d-scale 1.0 --match-li6-radius --quadrupole-target -0.0818 \
          --min-nn-separation 0.0 --format {he3,annotated} \
          --out li6_m+1.dat [--moments-json li6_m+1.moments.json]

  The flag is spelled **`--moments-json <path>`** everywhere (the first draft
  also called it `--moments`; one spelling only). It writes the analytic and
  sampled moments side by side — the fast path for filling §8's table — and it
  is also where the sidecar's `md5` and `git` fields are filled in (§5.2). The
  CLI is additionally the **only** place a human-facing warning can live (there
  is no `std::cerr` anywhere in `src/core`; setup errors throw, per
  `CONVENTIONS.md`), so `--theta-s 0` prints "a₂ ≡ 0 for a longitudinal axis"
  from Python while C++ simply returns 0.
* **Docs.** A `docs/USAGE.md` section ("⁶Li configuration tables for a
  dipole-model code"); the §11 row of `docs/OPEN_ITEMS_SOLUTIONS.md` gains the
  measured numbers and the §2.8 conclusions; `docs/PHYSICS_CHANNELS.md` gains a
  row under §8 (coherent) pointing at this module as *input generation, not a
  channel*; `data/vmc/density/README.md` as in §3.2.

---

## 7. Tests

C++ doctest `tests/test_cluster_config.cpp`, one `TEST_CASE` per identity, plus
`python/tests/test_cluster_config.py` for the bindings and the CLI.

**How the MC tolerances are set (this replaces the first draft's absolute
numbers, which were 0.1–0.3 σ and would have failed for most seeds).** The
per-configuration spread of every moment this module gates is large — a single
⁶Li configuration has Σ_i(3z_i²−r_i²) scattered over tens of fm² — so an
absolute gate is meaningless without the variance beside it. Measured here, on
a 2 × 10⁵-configuration MC of the design's own densities (default
`FitRescaled` + `VmcHe4Density`, and reproducing the analytic ⟨r²⟩ = 6.444 fm²
and 𝒬-based Q_matter = −1.231 fm² to well inside their errors):

| moment | per-configuration sd | 1 σ at N = 10⁵ | N for 5 σ at the *first draft's* gate |
|---|---|---|---|
| Σ_i(3z_i²−r_i²) [fm²] | **28.9** | 0.091 | T8's 0.02 → **5.2 × 10⁷** |
| Q_charge = ½Σ_i(…) [fm²] | 14.5 | 0.046 | T17's 0.005 → **2.1 × 10⁸** |
| ⟨r²⟩ per nucleon [fm²] | **3.66** | 0.0116 | T7's 0.003 → **3.7 × 10⁷** |
| r_rms [fm] | — | 0.0023 | T6/T18's 0.005 → 2.2 σ, i.e. *tighter* than 5 σ needs |
| ⟨x²⟩−⟨y²⟩ per nucleon [fm²] | **2.87** | 0.0091 | T11's 0.01 → 1.1 σ |

**Rule.** `ClusterConfigSet` carries `r2_sd_fm2`, `q_matter_sd_fm2`,
`delta_perp_sd_fm2` (§4), and **every MC gate below is written as
`5 · set.<moment>_sd / sqrt(n)`, never as a hand-typed absolute number**; the
value quoted in the table is what that evaluates to at N = 10⁵, for
orientation. A tolerance that must be tight and absolute belongs on the
**deterministic** layer instead: the sampler's own grid quadrature of (G1)
against the closed forms (G4)/(G5)/(G6) is a cell-CDF sum with no MC noise and
is gated at 1e-6 (T7a/T8a), which is where the first draft's 10⁻³ intent
actually lives.

**Reviewer note / response — the α core (T5).** The review proposed two ways
out of "the recentred VMC α does not reproduce `he4.density`": (a) gate ⟨s²⟩ and
merely record the χ², or (b) make the source *self-consistent* by iterating
ρ_src ← ρ_src · (ρ_target/ρ_recentred), "a 3–4 pass fixed point at
construction". **(b) was tried in this pass and is not adopted**: a five-pass
fixed point on the 1-D radial histogram, driven by a 2.5 × 10⁵-configuration MC
of the recentring, moved χ²/ndf from 162 to only ≈ 75 (against a target of ~1)
and left ⟨s²⟩ oscillating between 2.058 and 2.153 fm² — the ratio estimate is
itself MC noise, so the iteration does not converge at any statistics an
`O(1 s)` constructor can afford. The exact deconvolution is a genuine numerical
problem (solve f̃(3q/4)·f̃(q/4)³ = f̃_target(q) for f̃), not a three-line fixed
point, and it is out of scope for an opt-in module. **Adopted instead**: (a) for
the physics comparison — ⟨s²⟩ gated (it closes *exactly*, for any shape),
χ² against `he4.density` recorded — **plus a new, strictly better code gate**:
the sampled histogram is compared to the **closed-form recentred prediction**
f̃_s(q) = f̃(3q/4)·f̃(q/4)³, which is what the sampler is actually supposed to
produce and which a bug in the recentring would break. `alpha_cm_inflate` is
therefore kept as specified (it is exact in the second moment and exact in shape
for `AlphaCoreSource::Gaussian`), not replaced by an iterated source.

| # | what | gate |
|---|---|---|
| **T1** | (G1) through `theta_lm` + `clebsch_gordan` equals the closed forms (G2)/(G3) at cosθ ∈ {−1, −0.6, −0.2, 0, 0.35, 0.8, 1} for three (U,V) pairs | rtol **1e-14** (verified exact in this pass) |
| **T2** | `fdeut.av18` re-read: ∫(u²+w²)dr, ∫w²/∫, 𝒬/4, ½√⟨r_np²⟩ vs the file's own header `1`, `dstate`, `qm`, `rd` | **5e-5 relative.** Not 1e-5: the header numbers come from the Fortran code's own quadrature on h = 0.005, not from re-integrating the printed h = 0.01 table, and re-integration reproduces 0.999998 / 0.057599 / **0.2696703 vs 0.269673** (1.2e-5) / **1.9673386 vs 1.967364** (1.3e-5) |
| **T3** | `li6.ad` r-block: N₀, N₂, N₀+N₂ vs the file's printed `s-wave 0.838`, `d-wave 0.017`, `ndx 0.856` | 2e-3 absolute (measured 0.83709 / 0.01718 / 0.85428) |
| **T4** | **Sign conventions**, four separate assertions: (i) 𝒬, ⟨R²⟩, ρ_{αd} **unchanged** when *both* R₀, R₂ are negated; (ii) `q_int_fm2()` **exactly negates** and `q_dd_fm2()` is **unchanged** when only R₂ is negated, and `sign(𝒬)` flips; (iii) `sign(q_matter_analytic_fm2(±1)) == sign(LI6_QUADRUPOLE_FM2)` — the only *observable* check of the r-space S–D convention; (iv) the **no-i^L** Fourier–Bessel transform of (R₀, R₂) reproduces `li6.ad`'s k-block structure: ψ̃₀ node at 0.678 fm⁻¹, ψ̃₂/ψ̃₀ < 0 below it, > 0 above | (i),(ii) exact; (iii) exact; (iv) node within 0.01 fm⁻¹, signs exact. **Do not assert sign identity with the *physical* k-space amplitude** — it differs by (−i)² (§2.1) |
| **T5** | **α core**: 10⁵ recentred cores. Two gates and one record. (a) ⟨s²⟩ = 2.0774 ± 0.015 fm² — *exact* for any source shape, since recentring gives ⟨s²⟩ = (3/4)⟨v²⟩ and λ² = 4/3; (b) the radial histogram against the **closed-form recentred prediction** f̃_s(q) = f̃(3q/4)f̃(q/4)³ (§4.1 step 4), χ²/ndf < 2 on 40 grid-aligned bins — a code-correctness gate with no physics in it; (c) χ² against `he4.density` **recorded, not gated** | (a) 0.015 fm²; (b) χ²/ndf < 2; (c) recorded (measured ≈ 170 on 40 bins / 10⁶ entries, with +25 % at 0.25 fm, −7 % at 0.85–1.15, +5–10 % at 1.5–2.5, −15 % at 3.05) |
| **T5g** | The same with `AlphaCoreSource::Gaussian`, where the λ inflation is exact in **shape** as well | χ²/ndf < 1.5 against the analytic Gaussian |
| **T6** | **⁶Li one-body density**: 10⁵ configs (unpolarized), radial histogram of all 6 nucleons vs `li6.density`. A **shape** test at fixed normalization, **expected to fail a strict χ²** — the model is 4 % too large. Gate the *documented* discrepancy, not agreement | r_rms(sampled) = 2.5386 ± **0.012 fm** (5σ), ratio to `li6.density`'s 2.4433 = 1.039 ± 0.005; χ²/ndf recorded, not gated |
| **T7** | **⟨r²⟩ closure**: sampled ⟨r²⟩ equals (G6) | 5σ = 5·`r2_sd_fm2`/√n ≈ **0.058 fm²** at N = 10⁵ (analytic 6.4447 fm²) |
| **T7a** | The same closure on the **grid quadrature**, no MC: ∫ρ_αd R⁴ + the α and deuteron terms vs (G6) | **1e-6 relative** |
| **T8** | **Quadrupole**: sampled Σ_i(3z_i²−r_i²) for m = +1, 0, −1 equals (G5); the m = 0 value is −2× the m = ±1 value; the equal-thirds mixture is 0 | 5σ ≈ **0.46 fm²** at N = 10⁵ (m=+1 → −1.2309, m=0 → +2.4618, mixture → 0). Raise to N = 10⁶ for a 0.15 fm² gate if the run budget allows |
| **T8a** | The same on the **grid quadrature** — (G4) from the cell CDFs vs the analytic (G5) | **1e-6 relative** |
| **T9** | **The `coherent.hpp` bridge.** The **free function** `a2_from_quadrupole(2·0.269670, 2, |t|, m)` reproduces every row of `mantysaari_a2_deuteron()`. (The member `a2_from_geometry` takes no quadrupole and no A, so it cannot be fed the deuteron — §4 adds the free function for exactly this test) | m = ±1 within **10 %** at all four \|t\|; m = 0 within **25 %** (measured 8 % / 8–21 %) |
| **T10** | **D_T**, two separate assertions: (i) the exact CG identity `tensor_dilution() == 1 − 0.9 · p_d_alpha_d()`; (ii) agreement with `TaggedModel(li6_alpha_channel(BETA_DEFAULT, sampler.p_d_alpha_d()))` — the **Hulthén** path, where P_D is still a knob | (i) **1e-12**; (ii) **1e-4 relative**, the precedent `constants.hpp:85–91` sets for the same comparison (`population_integrated` is a 280 × 96 grid quadrature, not a closed form). **Not** `VmcAV18`, which pins P_D at `VMC_P_D_LI6` = 0.01935 against the sampler's 0.02011 (D_T 0.98259 vs 0.98190, 7e-4 apart) |
| **T11** | **Isotropy of the unpolarized set**: equal thirds of m = ±1, 0 give ⟨P₂(cosθ_i)⟩ = 0 over all nucleons, and ⟨x²⟩ = ⟨y²⟩ = ⟨z²⟩ | each \|⟨P₂⟩\| < 5σ; the three second moments equal within 5σ ≈ **0.045 fm²** at N = 10⁵ |
| **T12** | **CM constraint**: Σ_i r⃗_i = 0 for **every** configuration | max component < **1e-12 fm** |
| **T13** | **Reproducibility**: `sample_set(N, m, seed, run)` is bit-identical when re-run, when N is changed to 2N (first N rows), and when the configs are drawn individually via `Rng(seed, run, kConfigStream, i)` with the header's `kConfigStream` (§4) | exact double equality |
| **T14** | **Stream independence**: two different `run` values give sets whose ⟨r²⟩ differ by less than 5 MC σ but whose row 0 differs in every coordinate | statistical + inequality |
| **T15** | **Writer round-trip**: read back `He3Compatible` output with a Python parser that mimics `nucleons.cpp` (nine — here eighteen — bare `>>` reads); positions match to 1e-12; **assert the file contains no `#`, no blank line, and exactly 25 fields per row** | exact |
| **T16** | **Sidecar completeness**: all **13** `ClusterConfigOptions` fields appear under `options` in the `.meta.json` (reflection-free: the hand-maintained list `{alpha_source, alpha_d_source, theta_s, phi_s, alpha_cm_inflate, min_nn_separation_fm, alpha_d_scale, quadrupole_target_fm2, exact_coherence, n_r, n_c, r_max_fm, rnp_max_fm}`, so adding an option without recording it fails), and `md5`/`git` are present as keys with `null` from the C++ writer | exact |
| **T17** | **`quadrupole_target` closes.** Gate it on the **analytic** layer, where it is exact: with `quadrupole_target_fm2 = LI6_QUADRUPOLE_FM2`, `q_matter_analytic_fm2(1)/2 == −0.0818` to 1e-9, the applied root is s = 0.4032 and `p_d_alpha_d()` = 3.3e-3. The MC check is 5σ only, and note that at N = 10⁵ it **cannot resolve −0.0818 from 0** | analytic **1e-9**; MC 5σ ≈ **0.23 fm²** on Q_charge. Also: the √(target/model) rule returns −0.048, so a test that the implementation is *not* the √ rule is worth one line |
| **T17b** | **`quadrupole_target` bounds**: `validate()` throws for a target above +0.270 fm² (= +Q_d, the s = 0 floor) or below `q_matter_analytic_fm2(1)/2` (s = 1) | `std::runtime_error` |
| **T18** | **`match_li6_radius` closes**: with that scale the sampled r_rms is `LI6_R_POINT_VMC_FM` | 5σ ≈ **0.012 fm** at N = 10⁵ |
| **T19** | **Options validation throws**: `exact_coherence = true`; `n_r` or `n_c` < 8; a grid that would drop more than **1e-3 of any table's norm** | `std::runtime_error`. **Two rules removed from the first draft**: (a) "θ_S = 0 with a request for a₂ (warn)" — there is no warning channel in the core (no `std::cerr` in `src/core`), a₂ ≡ 0 is the *correct* value there, and the human-facing warning belongs in `configs.py` (§6); (b) "r_max below the table extent throws" — that fires on the **defaults** (`he4.density` runs to 20.05 > `r_max_fm` = 20, `fdeut.av18` to 100 > `rnp_max_fm` = 25). The norm-loss rule is the right form of the same check |
| **T20** | **Timing**: `sample_set(10⁵)` wall time recorded and asserted under a generous ceiling | < 5 µs/config (expect 1–2; `TaggedModel::sample_kc` measures 5.2 M draws/s for the same primitive, `tagged.hpp:214`) |
| **T21** | **The m_S conditioning of R̂** (§2.6) — the one error no moment test can see. For M = +1: ⟨P₂(cosθ_R) \| m_S = −1⟩ = **−2/7** and ⟨P₂(cosθ_R) \| m_S = 0⟩ = **+1/7** exactly (pure L = 2, m_L = 2 and 1). Gate the **grid tables** at 1e-6 (deterministic) and the MC draws at 5σ | grid 1e-6; MC 5σ. An m_S-marginal draw gives ⟨P₂⟩ = −0.0355 in *both* branches and fails by ~15 σ even at N = 10⁴ |
| **T22** | **`asymptotic_ds_ratio()`**: the Whittaker-divided η of §2.7 from the `FitRescaled` tables | −0.05 ± 0.01, and **> the measured −0.025(12) in magnitude by ≈ 2×, not 5–15×** — the assertion is the *number*, the factor is documented in O1 |

pytest mirrors: bindings return the right array shapes and dtypes ((N, 6, 3)
float64 and (N, 6) int); the CLI writes a file with the right row count;
`--m unpolarized` splits N into thirds whose **counts differ by at most 1**
(N = 10⁵ gives 33334/33333/33333, not equal thirds); the moments JSON matches
the C++ analytic predictors; `CONFIG_STREAM` and the two radius constants are
importable.

---

## 8. Numbers to fill (all values below were computed during this design pass; the implementer must reproduce them, not re-derive them)

| quantity | `FitRescaled` (default) | `OverlapRaw` | `FitRaw` |
|---|---|---|---|
| P_D(α–d) | 0.02011 | 0.02011 | 0.02542 |
| ⟨R²⟩_αd [fm²] (rms) | 16.963 (4.119) | 17.548 (4.189) | 16.956 (4.118) |
| 𝒬[R₀,R₂] [fm²] | −1.3204 | −1.4009 | −1.4895 |
| D_T | 0.98190 | 0.98190 | 0.97712 |
| **⟨r²⟩ [fm²] → r_rms [fm]** | **6.4447 → 2.5386** | 6.5747 → 2.5641 | 6.4432 → 2.5383 |
| **Q_matter(m = ±1) [fm²]** | **−1.2309** | −1.3382 | −1.4590 |
| Q_matter(m = 0) [fm²] | +2.4618 | +2.6764 | +2.9180 |
| Q_charge(m = +1) [fm²] | −0.6154 | −0.6691 | −0.7295 |
| δ⊥ = ⟨x²⟩−⟨y²⟩ per nucleon [fm²] | −0.1026 | −0.1115 | −0.1216 |
| eps_b0 equivalent (B = 52.0 GeV⁻²) | −0.0506 | −0.0550 | −0.0600 |
| a₂(±1) at \|t\| = 0.3 | +0.197 | +0.215 | +0.234 |
| α–d share of \|Q\| / of which S–D interference | 76.9 % / 95.2 % | 77.9 % / 95.1 % | 79.0 % / 94.6 % |
| 𝒬 split: interference / pure-D [fm²] | −1.2570 / −0.0631 | −1.3324 / −0.0684 | −1.4094 / −0.0798 |
| asymptotic D/S ratio η (Whittaker-divided, R = 6–8 fm) | −0.042…−0.053 | −0.043…−0.072 | as `FitRescaled` |
| `quadrupole_target` root s for Q_charge = −0.0818 (→ P_D′) | **0.4032** (3.3 × 10⁻³) | 0.3806 (3.0 × 10⁻³) | 0.3577 (3.3 × 10⁻³) |
| per-config sd of Σ_i(3z_i²−r_i²) [fm²] | 28.9 | ≈ 29 | ≈ 29 |
| per-config sd of ⟨r²⟩ per nucleon [fm²] | 3.66 | ≈ 3.7 | ≈ 3.7 |

Every `Fit*` column is for the fit **truncated at its last abscissa, 9.95 fm**
(§3.1): tables are zero beyond their own range, never clamped. The raw block's
R > 10 fm tail — 0.087 % of the norm, 0.106 fm² (0.60 %) of ⟨R²⟩ — is the size
of that truncation and is listed among §2.7's systematics. The `η` and
`quadrupole_target` rows are new in this revision; the `OverlapRaw` / `FitRaw`
roots are indicative (recomputed from the same (G9) with each source's own
𝒬_int, 𝒬_DD and P_D) and the implementer reproduces them, as with every other
number here.

Reference points: r_point(⁶Li) measured **2.4655 fm** (`LI6_R2_POINT_FM2`),
VMC `li6.density` **2.4433 fm** (`LI6_R_POINT_VMC_FM`); Q_charge measured
**−0.0818 fm²** (`LI6_QUADRUPOLE_FM2`), **GFMC AV18+IL7 −0.20(6) fm²** (Pastore
et al., PRC 87, 035503 (2013), arXiv:1212.3375 — *not* VMC, and a different
Hamiltonian from the AV18+UX tables this module reads); deuteron Q_d
**+0.26967 fm²**, ⟨r_np²⟩ **15.4817 fm²**, P_D(d) **0.057599**, asymptotic
η_d **0.025045** (`fdeut.av18` header, for scale against the α–d η above);
α ⟨s²⟩ **2.0774 fm²** (`he4.density`) or **2.1003 fm²** (Gaussian,
`cluster_point_a2_fm2(2,4)`).

**Timing** — to be measured and written back: target < 5 µs per configuration
single-threaded, 10⁵ configurations in well under a second; report
µs/config and the wall time of the reference `lipolgen-configs --n 100000` run.

---

## 9. Implementation steps (one Opus agent, in this order)

1. **Data.** Fetch `he4.density` and `li6.density` per §3.2 into
   `data/vmc/density/`; write `data/vmc/density/README.md` in the same table
   shape as `data/vmc/README.md` (live URL, Wayback timestamp, fetch date,
   md5, printed normalization and rms). Verify the two md5s in §3.2.
2. **Reader.** Add `read_anl_plain` to `cluster.hpp` + `cluster.cpp` (§3.3);
   one doctest that it reads `li6.adr.fit` into 2 value columns and that
   `read_anl_momentum` on `he4.density` gives `col[0]` = ρ, `err[0]` = δρ.
   Existing tests must stay green bit-for-bit.
3. **Header.** Write `include/lipolgen/cluster_config.hpp` as §4, with the
   physics of §2 in the file comment and a pointer to this document.
4. **Densities and analytics.** `cluster_config.cpp`: table loading (**zero
   outside every table's own abscissa range**, §3.1) + the ρ builders from
   (G1)/(G1′) (via `theta_lm` / `clebsch_gordan`, no hand-expanded formulas),
   then `p_d_alpha_d`, `tensor_dilution`, `q_matter_analytic_fm2`,
   `q_int_fm2`, `q_dd_fm2`, `r2_analytic_fm2`, `delta_perp_analytic_fm2`,
   `asymptotic_ds_ratio`, the free `a2_from_quadrupole` and the member
   `a2_from_geometry`, `eps_b0_equivalent`, `quadrupole_band_fm2`. Add
   `LI6_R2_POINT_FM2` to `rc.hpp` (§1). **Land T1–T4, T7a, T8a, T9, T10, T22
   here** — the whole analytic layer is testable before a single configuration
   is drawn.
5. **Sampler.** Grids + cell CDFs (mirror `TaggedModel`'s eager,
   immutable-after-construction pattern), **one (R, cosθ_R) table per
   (M, m_S)** as `build_amp2` does (§2.6), the six-step algorithm of §4.1, the
   CM constraint, the rotation, and the per-configuration moment **variances**
   `ClusterConfigSet` reports. Land T5, T5g, T6–T8, T11–T14, T19, T20, T21.
6. **Writer.** `write_snd_configs` + the sidecar: a **hand-rolled flat JSON
   object in `cluster_config.cpp`** (no library, no hash, no git — §4/§5.2;
   `tests/json_min.hpp` is a reader and there is no emitter in the core), with
   `md5` and `git` emitted as `null` for `configs.py` to fill. Land T15, T16.
7. **Scaling knobs.** `alpha_d_scale`, `match_li6_radius`, and
   `quadrupole_target_fm2` as the **(G9) root-find with its reachability
   check** (§2.7) — not √(target/model). Land T17, T17b, T18.
8. **Bindings + CLI** (`python/bindings.cpp`, `python/lipolgen/configs.py`,
   `pyproject.toml` entry point) + pytest.
9. **Docs.** `docs/USAGE.md` section; `OPEN_ITEMS_SOLUTIONS.md` §11 rewritten
   with §2.8's three conclusions and §8's table; a `PHYSICS_CHANNELS.md` row;
   fill §8's timing row from the reference run.
10. **Gate.** `build/lipolgen_tests` and `pytest python/tests` green with the
    pre-existing counts unchanged plus the new cases; every rtol 1e-12
    reference gate untouched; `git diff` shows no edit to any existing physics
    default. One commit.

Files created: `include/lipolgen/cluster_config.hpp`,
`src/core/cluster_config.cpp`, `tests/test_cluster_config.cpp`,
`python/lipolgen/configs.py`, `python/tests/test_cluster_config.py`,
`data/vmc/density/{he4.density,li6.density,README.md}`.
Files edited (additively only): `include/lipolgen/cluster.hpp`,
`src/core/cluster.cpp` (one reader), `include/lipolgen/rc.hpp` (one constant,
`LI6_R2_POINT_FM2`), `python/bindings.cpp` (one `bind_cluster_config` + a 3-D
array helper), `CMakeLists.txt` (the one-line `lipolgen-configs` shim, §6),
`pyproject.toml` (one script entry), `docs/USAGE.md`,
`docs/OPEN_ITEMS_SOLUTIONS.md`, `docs/PHYSICS_CHANNELS.md`.

---

## 10. Open questions, and the exact ask to send the Mäntysaari group

**Open, and not settled in this pass.**

* **O1. The quadrupole gap is a D-wave excess, and it is moderate.** The α+d
  truncation gives Q_charge = −0.615…−0.730 fm² against a measured −0.0818.
  The decisive check the first draft deferred **was done in this revision**
  (§2.7): dividing the table's R₂/R₀ by the Whittaker ratio W₂/W₀ — which is
  3.3 at 6 fm and 2.6 at 8 fm, so the naive ratio is *not* the asymptotic
  D/S ratio — gives **η = −0.05(1)** against the measured **η = −0.025 ± 0.006
  ± 0.010** (George & Knutson, PRC 59, 598 (1999), from d + α elastic tensor
  analysing powers; other determinations −0.01…−0.03). The tail D-wave is
  therefore **≈ 2× the measured value, not 5–15×**, which changes the reading:
  hypothesis **(c)** — an ANL L = 2 normalization or phase convention different
  from the standard [Y_L ⊗ χ_1]^{JM} one — is **no longer the leading
  candidate** and is demoted; what remains is (a) the missing 15 % non-α+d
  component plus (b) core/cluster polarization (a free deuteron's own quadrupole
  and a free α–d D wave inside a nucleus that compresses and polarizes both),
  amplified into Q by the r⁴ weight of (G4) and by the fact that 𝒬 is *linear*
  in the D amplitude. Still open: how a moderate D-wave excess and a 15 %
  missing component together make a factor 7.5 in Q. **The rule stands
  regardless — do not quote a tensor number from the sampler without the band.**
* **O2.** Does the smoothed `li6.adr.fit` R₂ node near 1.1 fm reflect
  antisymmetrization at short α–d distance, or is it a fit artefact? It
  contributes little to 𝒬 (the r⁴ weight kills it) but it decides whether
  `FitRescaled` or `OverlapRaw` is the honest default.
* **O3.** The α core is sampled as an uncorrelated product of one-body
  densities. GFMC ⁴He *configurations* (with correlations) exist and are what
  arXiv:2605.00454 actually uses; how much does the incoherent/coherent split
  move when correlations are restored? Unanswerable from our side.
* **O4.** `CoherentScenario::eps_b0`'s band, and the fact that **`coherent.hpp`
  uses the symbol ΔB without defining it** (§2.8): under three defensible
  readings of ΔB_m the docstring's "ΔB₀/B" is off from `a2_m_state` by −2, by
  −1, or not at all. Fix the *convention* first, and phrase the repaired
  docstring in the unambiguous quantity **eps_b0 = δ_{±1}/B with δ_{±1} =
  ⟨x²⟩−⟨y²⟩ per nucleon** — which is exactly `eps_b0_equivalent()`. Note also
  that the ⁶Li estimates are quoted at B = 52.0 GeV⁻² while the deuteron
  reference point that motivated the band uses B_d = 33.1 GeV⁻². A separate,
  reviewed decision — flagged, not taken here.
* **O5.** Whether the ⁶Li a₂ at the measured quadrupole (+0.026 at |t| = 0.3,
  ~10× below the deuteron's and opposite in sign) survives the EIC's
  statistics at all, and how it separates from the linearly-polarized-photon
  cos 2φ background (`physics_literature.md`'s two-mechanism warning). This is
  the question that decides whether the whole collaboration is worth
  proposing, and it is answerable *today* from §2.8 plus a rate estimate — do
  that before writing the letter.

**How the graft would consume the table.** Our writer emits one line per
configuration in `he3.dat`'s own layout, in fm, ion rest frame, c.m. at the
origin, with the polarization axis already applied — so the consumer needs no
knowledge of our conventions beyond "these are the six nucleon positions."
Their side then needs the §5.3 patch: replace the `A == 3` special case and the
hard-coded `"he3.dat"` with `-configfile <path> -configid <n>` for any A. Their
Good–Walker loop is unchanged; the amplitude uses only x and y. One set per
m ∈ {+1, 0, −1} plus one unpolarized set; a₂ comes out of the Φ dependence of
⟨A⟩ exactly as in their Fig. 2.

**The ask, in one paragraph.**

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

