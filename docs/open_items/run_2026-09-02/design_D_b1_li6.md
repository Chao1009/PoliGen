# Design D — b₁(⁶Li) from a three-term α–d convolution

**Open item 10** (`docs/OPEN_ITEMS_SOLUTIONS.md` row 10 and §"8–10"; recommendation in
`docs/open_items/physics_literature.md` §2). Design only — no code changed by this
document. Author: supervising session (Fable), 2026-09-02.

**Kernel paper (local).** W. Cosyn, Y.-B. Dong, S. Kumano, M. Sargsian,
*Tensor-polarized structure function b₁ in standard convolution description of
deuteron*, PRD **95** (2017) 074036, [arXiv:1702.05337](https://arxiv.org/abs/1702.05337)
— `/home/cpeng/Projects/polli/PolarizedLithiumSim/refs/1702.05337.pdf`. Cited below
as **CDKS**, by their equation numbers.

**What the supervisor asked for, restated as the deliverable.** One opt-in
`TensorSF` backend, working name `Li6ConvolutionB1`, that computes b₁ of ⁶Li per
nucleon from the α–d convolution — (1) the embedded deuteron's own b₁ convolved with
the α–d light-cone density, (2) the α–d relative D-wave (CDKS Eq. 21 with α–d
overlaps), which has **two** pieces, (2d) struck deuteron × F₁ᵈ and (2α) struck α ×
F₁^α, and (3) the Clebsch–Gordan depolarization of that deuteron inside the L = 2
component — normalised to the VMC α–d spectroscopic factor, with each term and the
α–d D-wave weight as knobs and a **mandatory 100 % band**; blocked behind an
**A = 2 validation gate** against the digitized CDKS Fig. 4.

**Revision 2 (2026-09-02, after adversarial review).** The recommendation quoted below
says "three terms". It is **four**: CDKS Eq. (10) sums the spectral function over
constituents, and one level up that sum runs over {d, α}, so the struck-α orbital term
(2α) exists and is ≈ 0.5 × term (2d) (§1.5, computed). The other structural changes of
this revision are the finite-|q⃗| δ-function (§1.3), the light-cone density
normalisation identities (§5.2), the resolution of the per-nucleon question (§9 Q1),
and one home for N_{αd} (§2.1).

> `physics_literature.md` §2, verbatim: *"Recouple |⁶Li; 1H⟩ in the α–d basis exactly
> as Cosyn Eq. (19) does for the deuteron, with L = α–d relative orbital angular
> momentum and S = deuteron spin. Three terms result: (1) embedded deuteron — Eq. (16)
> with F₁ᴺ → b₁ᵈ, convolved with the α–d light-cone density f_{d/⁶Li}(z) from Eq. (17)
> (peaked near z = 1/3); (2) α–d D-wave — Eq. (21) unchanged in form with φ₀,φ₂ →
> φ₀^{αd}, φ₂^{αd} and F₁ᴺ → F₁ᵈ; (3) CG depolarization of the deuteron inside the
> L=2 component, an analytic P₂-weighted Clebsch–Gordan sum."*

---

## 1. Physics

### 1.1 Symbol table

| symbol | meaning | where it comes from |
|---|---|---|
| `x` | per-nucleon Bjorken variable of the **target**, x = Q²/(2 M_N ν) | the library's convention everywhere (`xsec.hpp`: "All structure functions are per-nucleon (F2A/A)") |
| `Q²`, `ν` | photon virtuality, energy transfer, target rest frame | kernel arguments |
| `M_N` | free nucleon mass 0.9383 GeV | `constants.hpp::M_NUCLEON` |
| `M_d`, `M_α`, `M_{6Li}` | 1.875612942, 3.727379407, 5.601518702 GeV | `nuclear_mass(1,2)`, `nuclear_mass(2,4)`, `nuclear_mass(3,6)` (`spectator.hpp:49`) — the public accessor. **Not** `nucleus_table()` (anonymous namespace) and **not** `LI6_ALPHA_TAG().m_spec()`, which is the 10-keV-rounded 3.72738 |
| `ε_d` | deuteron separation energy, 2.224574 MeV | the `fdeut.av18` header's own `ebind`, returned by `read_fdeut_k` (§2.3). `DEUTERON_P_TAG().separation_energy` = 2.2246e-3 is the same number rounded; the 2.6e-5 relative difference is recorded in the gate report, not silently chosen |
| `ε_{αd}` | α–d separation energy of ⁶Li, 1.4743 × 10⁻³ GeV | `LI6_ALPHA_TAG().separation_energy` (`spectator.cpp:91`) — **do not re-type** |
| `p⃗` (A = 2) | nucleon momentum in the deuteron rest frame | integration variable |
| `k⃗` (A = 6) | **α–d relative momentum** = deuteron momentum in the ⁶Li rest frame (the α carries −k⃗) | VMC tables |
| `θ` | angle between p⃗ (resp. k⃗) and q⃗ | integration variable |
| `κ` | \|q⃗\|/ν = √(1 + γ²), γ² = 4M_N²x²/Q² | `asymmetries.hpp::gamma_squared`; **§1.3** — it sits inside CDKS's Eq. (21) δ-function and is 1.38 at x = 0.8, Q² = 2.5 |
| `y` (A = 2) | nucleon light-cone fraction, CDKS Eq. (18) | see §1.3 |
| `z_d`, `z_α` (A = 6) | deuteron / α light-cone fraction of ⁶Li, **each normalised to its own fair share** (1/3, 2/3) | see §1.6 |
| `φ₀, φ₂` | S- and D-state momentum-space radial wave functions, CDKS convention (with the iᴸ factor, so φ₂ < 0 for the deuteron) | §1.2, §2 |
| `U(k), W(k)` | the same without the iᴸ factor, φ₀ = U, φ₂ = −W; both > 0 at low k for the deuteron | `fdeut.av18` |
| `f(y)`, `δ_T f(y)` | unpolarized and tensor light-cone densities | CDKS Eqs. (17), (16) |
| `f_α(z)`, `δ_T f_α(z)` | the same for the **struck α** (m_struck = M_α, m_recoil = M_d) | §1.6, term (2α) |
| `P_D^{αd}` | α–d relative D-state probability | VMC, `VMC_P_D_LI6` = 0.0193549 (§2.1) |
| `N_{αd}` | α–d spectroscopic factor | VMC, `VMC_N_ALPHA_D_LI6` = 0.819481 (§2.1) |
| `H` | spin projection of the spin-1 system (⁶Li, or the deuteron in the A = 2 case) | ±1, 0 |
| `m_S` | spin projection of the **deuteron** inside ⁶Li | CG sum |
| `m_L` | projection of the α–d relative orbital angular momentum | CG sum |
| `P_zz` | tensor polarization p⁺ + p⁻ − 2p⁰ | CDKS Eq. (46) |

### 1.2 CDKS, transcribed

**Eq. (9)–(12), the impulse-approximation setup.** The nuclear hadron tensor is the
nucleonic one convolved with the spectral function,
W^A_{μν}(P_A, q) = ∫d⁴p S(p) W_{μν}(p, q) (Eq. 9), with, in a simple shell model,
S(p) = (1/A) Σ_i |φ_i(p⃗)|² δ(p⁰ − M_A + √(M²_{A−i} + p⃗²)) (Eq. 10), separation energy
ε_i = (M_{A−i} + M_N) − M_A (Eq. 11), and — this is the approximation CDKS actually
evaluate — the **non-relativistic energy**

> **p⁰ = M_N − ε − p⃗²/(2M_N).**  (CDKS Eq. 12)

*(The p⃗²/(2M_N) is really p⃗²/(2M_{A−i}); for the deuteron the residual system is a
nucleon, so the two coincide. For α–d in ⁶Li the residual system is the α, and the
recoil term is k⃗²/(2M_α) — see §1.6.)*

**Eq. (13)–(15), what b₁ is.** With the helicity amplitudes
A_{hH,hH}(x,Q²) = ε*^μ_h ε^ν_h W^D_{μν}(p_D,q) (Eq. 13) and the photon polarization
vectors ε^μ_{h=±1} = (1/√2)(0, ∓1, −i, 0), ε^μ_{h=0} = (1/√Q²)(|q⃗|, 0, 0, q⁰) (Eq. 14),

> **b₁ = A_{+0,+0} − (A_{++,++} + A_{+−,+−})/2 |_LT ,
> F₁ᴺ = (A_{+↑,+↑} + A_{+↓,+↓})/2.**  (CDKS Eq. 15)

Note the asymmetry that fixes every factor of 2 downstream: **b₁ is a difference of
un-averaged amplitudes, F₁ᴺ an average over the nucleon spin.**

**Eq. (16), the convolution — the master equation of this design.**

> **b₁(x, Q²) = ∫ (dy/y) δ_T f(y) F₁ᴺ(x/y, Q²),   δ_T f(y) ≡ f⁰(y) − [f⁺(y) + f⁻(y)]/2.**  (CDKS Eq. 16)

CDKS: *"Here, the structure function b₁ is defined by the one per nucleon."*

**Eq. (17), the light-cone density.**

> **f^H(y) = ∫ d³p y |φ^H(p⃗)|² δ( y − (E − p_z)/M_N ),**  (CDKS Eq. 17)

E ≡ p⁰ of Eq. (12); f^H sums over the struck nucleon's spin,
f^H(y) ≡ f^H_↑(y) + f^H_↓(y), because b₁ is an unpolarized-quark observable.

**Eq. (18), the light-cone variable.**

> **y = M p·q /(M_N P·q) ≃ 2p⁻/P⁻,  p⁻ ≡ (p⁰ − p³)/√2,**  (CDKS Eq. 18)

collinear frame with q⃗ along +z, so the **minus** components survive; M is the
deuteron mass, so y is the nucleon fraction **normalised to its fair share** (⟨y⟩ ≈ 1),
not to the whole nucleus.

**A discrepancy inside the paper that this design must decide.** Eq. (17)'s
δ-argument is (E − p_z)/M_N; Eq. (21)'s (transcribed below, verified in the PDF) is
p·q/(M_N ν) = (E − p_z κ)/M_N with κ ≡ |q⃗|/ν = √(1 + γ²). The first is the κ → 1
Bjorken limit of the second, which is also the "≃" in Eq. (18). κ = 1.06 at x = 0.3
and 1.38 at x = 0.8 (Q² = 2.5), so this is **not** a rounding choice. §1.3 implements
both and makes the choice an option; the κ = 1 form is the default and the κ form is
item 0 of the §5.4 checklist.

**Eq. (19)–(20), the wave function.**

> **φ^H(p⃗) = φ₀(p) Y₀₀(p̂) χ_H + Σ_{m_L} ⟨2 m_L : 1 m_S | 1 H⟩ φ₂(p) Y_{2m_L}(p̂) χ_{m_S},**  (CDKS Eq. 19)
> **ψ^H(r⃗) = [u₀(r)/r] Y₀₀(r̂) χ_H + Σ_{m_L} ⟨2 m_L : 1 m_S | 1 H⟩ [u₂(r)/r] Y_{2m_L}(r̂) χ_{m_S},**  (CDKS Eq. 20)

with m_S = H − m_L, D-state probability ∫dp p²|φ₂(p)|², and
**φ_L(p) = 4π iᴸ ∫dr r² j_L(pr) u_L(r).** CDKS state explicitly that the iᴸ makes
**φ₂(p) < 0** for the deuteron, "although a different convention (φ₂ → −φ₂), namely
without the i² factor, is sometimes used". Later (below their Eq. 40) they fix
**φ₀(k) = U(k) and (−i)²φ₂(k) = W(k), i.e. φ₂ = −W** with U, W ≥ 0 at low k. **This
sign is load-bearing** (§1.5, §2.3).

**Eq. (21), the tensor density — the kernel of terms (2d) AND (2α).**

> **δ_T f(y) = ∫ d³p y [ −(3/(4√2 π)) φ₀(p)φ₂(p) + (3/(16π)) |φ₂(p)|² ] (3cos²θ − 1) δ( y − p·q/(M_N ν) ).**  (CDKS Eq. 21)

The first bracket term is the **SD interference**, the second the **DD** term. In the
sign-unambiguous (U, W) form, using φ₀φ₂ = −UW,

> **δ_T f(y) = ∫ d³p y (3/(4π)) [ U(p)W(p)/√2 + W(p)²/4 ] (3cos²θ − 1) δ(…),**  (Eq. 21′)

which matches the bracket of their VNA Eq. (42)/(44) — a useful internal consistency
check that the two theories share one convention.

**Verification of Eq. (21) done for this design (reproduce it in the unit tests).**
Both coefficients follow from Eq. (19) and the CG table of §1.4:
* SD: δ_T of ⟨2 0 : 1 H|1 H⟩ is C⁰₀ − (C⁺¹₀ + C⁻¹₀)/2 = −√(2/5) − √(1/10) = −0.9486833;
  the interference is 2 φ₀φ₂ Y₀₀Y₂₀ × that = 2·(1/√4π)·√(5/16π)·(−0.9486833)(3c²−1)
  = **−0.1688094 (3c²−1) φ₀φ₂**, and −3/(4√2π) = **−0.1688094**. ✓
* DD: δ_T of Σ_{m_L}|C^H_{m_L}|²|Y_{2m_L}|² = (3/10)|Y₂₀|² + (3/10)|Y₂₁|² − (3/5)|Y₂₂|²
  = **(3/16π)(3c² − 1)**, exactly CDKS's coefficient. ✓

**Eq. (22), the nucleon F₁ the convolution must be fed.**

> **F₁ᴺ(x,Q²) = [1 + 4 M_N² x²/Q²] / (2x[1 + R(x,Q²)]) · F₂ᴺ(x,Q²),  F₂ᴺ|_LO = x Σ_i e_i²[q_i + q̄_i]|_LO.**  (CDKS Eq. 22)

**This is not the library's `UnpolSF::f1_from_f2`,** which is the massless
F₁ = F₂/(2x(1+R)). The extra factor is exactly `1 + gamma_squared(x, q2)`
(`asymmetries.hpp:60`), and it is **large** at the kinematics where b₁ is biggest:
γ² = 0.90 at x = 0.8, Q² = 2.5. Feeding the massless F₁ into the A = 2 gate changes
the peak of x·b₁ by ~50 % (measured in the prototype of §5.4). **The convolution's
nucleon input is CDKS Eq. (22), spelled once, in the new header, as
`(1.0 + gamma_squared(x, q2)) * f2 / (2*x*(1+R))` with a comment pointing at
`f1_from_f2` and saying why they differ.** CDKS use the SLAC-R1998 R — the library
already has it (`sf.hpp::r1998`, `R1998Form::kAverage`) — and MSTW2008 LO PDFs.
F₁ᴺ = (F₁ᵖ + F₁ⁿ)/2 (CDKS, below Eq. 18).

**Eqs. (34)–(46), the light-front virtual-nucleon variant ("theory 2").** Transcribed
because the design must state which theory it implements and because Eqs. (45)–(48)
are the tensor conventions the rest of LiPolGen already carries.

> **√(2/3) F_{UT_LL,T}/F_{UU,T} = [A_{++,++} − 2A_{+0,+0} + A_{+−,+−}]/[A_{++,++} + A_{+0,+0} + A_{+−,+−}] = −(2/3) b₁^EPW/F₁.**  (Eq. 34)
> **W^{λ'λ}_{μν}(P,q) = 4(2π)³ ∫ dΓ_N (α_N/α_i) W^N_{μν}(p_i,q) ρ_D(λ',λ),**  (Eq. 35)
> **α_i = 2p⁻_i/P⁻, α_N = 2p⁻_N/P⁻ = 2 − α_i,**  (Eq. 36)
> **k = √(E_k² − M_N²),  E_k² = (m_N² + k⃗_⊥²)/(α_i(2−α_i)),  k³ = (1−α_i)E_k,  k⃗_⊥ = p⃗_i^⊥ + (α_i/2)P⃗^⊥,**  (Eq. 37)
> **dΓ_N = d³p_N/(2E_{p_N}(2π)³) = dα_i d p⃗_i^⊥/(2α_i(2π)³) = α_i d³k⃗/(2E_k(2π)³),**  (Eq. 38)
> **ρ_D(λ',λ) = Σ_{λ_N,λ'_N} [Ψ^D_λ'(k⃗,λ'_N,λ_N)]† Ψ^D_λ(k⃗,λ'_N,λ_N)/(α_N α_i),**  (Eq. 39)
> **Ψ^D_λ(k⃗,λ₁,λ₂) = √E_k Σ_{λ'₁λ'₂} D^{1/2}_{λ₁λ'₁}[R(k₁/m_N)] D^{1/2}_{λ₂λ'₂}[R(k₂/m_N)] Σ_{l=0,2; λ_S} ⟨l λ_l : S λ_S | j λ⟩ ⟨s₁λ'₁ : s₂λ'₂|1λ_S⟩ Y_{lλ_l}(Ω_k)(−i)ˡ φ_l(k),**  (Eq. 40)

with D^{1/2} the **Melosh rotation** matrices, and the wave functions
"approximated by the nonrelativistic ones φ₀(k) = U(k) and (−i)²φ₂(k) = W(k)". The
light-front wave function satisfies the baryon-number and momentum sum rules

> **Σ_{λ₁λ₂} ∫ dα_i dk⃗_⊥/(α_i(2−α_i)) |Ψ^D_λ(α_i,k⃗_⊥,λ₁,λ₂)|² = 1,
> Σ_{λ₁λ₂} ∫ dα_i dk⃗_⊥/(α_i(2−α_i)) α_i |Ψ^D_λ(α_i,k⃗_⊥,λ₁,λ₂)|² = 1.**  (Eq. 41)
> **F_{UT_LL,T} = −∫ (k²dk/α_i) d(cosθ_k)[F₁ᴺ(x_i,Q²) − (T²/(2p_i·q))F₂ᴺ(x_i,Q²)] √(3/2)[U(k)W(k)/√2 + W(k)²/4][3cos(2θ_k) + 1],
> F^{cos2φ_{T⊥}}_{UT_TT} = −∫ (k²dk/α_i) d(cosθ_k)(−T²/(2p_i·q))F₂ᴺ(x_i,Q²) √(3/2)[U(k)W(k)/√2 + W(k)²/4] sin²θ_k,**  (Eq. 42)
> **T^μ = p^μ_N + ((p_N·q)/Q²)q^μ − ((p_N·L)/L²)L^μ,   L^μ = P^μ + ((P·q)/Q²)q^μ,**  (Eq. 43)
> **b₁(x,Q²) = (3/(4(1+γ²))) ∫ (k²dk/α_i) d(cosθ_k) [F₁ᴺ(x_i,Q²)(6cos²θ_k − 2) − (T²/(2p_i·q))F₂ᴺ(x_i,Q²)(5cos²θ_k − 1)] [U(k)W(k)/√2 + W(k)²/4],**  (Eq. 44)

with x_i = Q²/(2p_i·q) ≃ x/α_i.

> **dσ/dxdQ² = (dσ^U/dxdQ²)(1 + ½ P_zz A_zz),**  (Eq. 45)
> **P_zz = √6 T̃_zz = p⁺ + p⁻ − 2p⁰,**  (Eq. 46)

and, in the scaling limit with the Callan–Gross relations (their Eqs. 47–49),
**A_zz = −(2/3) b₁/F₁** (Eq. 48) — which is exactly this repo's
`TENSOR_LL_SIGN = −1` (`docs/CONVENTIONS.md`, `OPEN_ITEMS_SOLUTIONS.md` §3). **This
design implements theory 1 (Eqs. 16/17/21) only.** Theory 2 is transcribed here so
that a later implementer has the finite-γ variant in one place, and because Eq. (42)'s
bracket is the sign cross-check used in §1.2.

### 1.3 Non-relativistic reduction and the y support (A = 2)

CDKS write the δ-function two ways (§1.2). Keeping **κ ≡ |q⃗|/ν = √(1 + γ²)** as
Eq. (21) writes it, with c ≡ cos θ,

  y = (E − p_z κ)/M_N = 1 − ε/M_N − p²/(2M_N²) − (p κ/M_N) c ,
  c*(p, y) = [ M_N(1 − y) − ε − p²/(2M_N) ] / (p κ) .

The δ fixes c and contributes a Jacobian M_N/(p κ), so **every angular integral
collapses analytically**:

  δ_T f(y) = (3/(2κ)) y M_N ∫ p dp [U W/√2 + W²/4] (3c*² − 1) ,
  f(y)     = (1/(2κ)) y M_N ∫ p dp [U² + W²] ,

both over the **interval** in p where |c*| ≤ 1. With B ≡ M_N(1−y) − ε and
r ≡ √(κ² + 2B/M_N),

> **p ∈ [ M_N |κ − r| , M_N (κ + r) ] , real only for B ≥ −M_N κ²/2 ,
> y ≤ y_max = 1 − ε/M_N + κ²/2 .**

*Both* endpoints matter: getting only the lower one wrong inflates ∫f dy by a factor
6 (seen in the prototype). At κ = 1 this is the familiar
p ∈ [M_N|1 − r|, M_N(1 + r)] with y_max = 3/2 − ε/M_N = **1.4976**; at x = 0.8,
Q² = 2.5 it runs instead to **1.9502** — i.e. with κ the support is **x- and
Q²-dependent**, and `y_max()` must take κ as an argument.

**Which one to implement: both, behind `finite_q_delta` (§3), default `false`.** The
κ = 1 form is Eq. (17) as printed, it is what makes f(y) a function of y alone — one
cached table reused at every x, which is the whole reason `LightConeDensities` exists
— and it keeps this design's default reproducible by anyone reading Eq. (17).
`finite_q_delta = true` rebuilds the densities per x with κ(x, Q²) and is **item 0 of
the §5.4 checklist**, because it is worth a factor 1.5–1.7 at the large-x peak.

**Measured for this revision** (AV18 u, w; CDKS Eq. (22) F₁ᴺ on `ToyF2` with
R = `r_sigma_lt`; the κ = 1 column reproduces §5.4's prototype table to 1–3 %, which is
how this recomputation was validated):

| x | x·b₁, κ = 1 | x·b₁, κ = √(1+γ²) | ratio |
|---|---|---|---|
| 0.10 | −1.638e−5 | −1.672e−5 | 1.02 |
| 0.20 | −2.030e−5 | −2.130e−5 | 1.05 |
| 0.30 | −1.596e−5 | −1.507e−5 | 0.94 |
| 0.50 | +9.364e−5 | +1.445e−4 | 1.54 |
| 0.70 | +2.843e−4 | +4.523e−4 | 1.59 |
| 0.80 | +2.876e−4 | +4.802e−4 | 1.67 |
| ∫b₁ dx, x ∈ [0.02, 1.2] | +1.240e−4 | +2.409e−4 | 1.94 |

so κ alone closes about **half** the ∫b₁ dx gap to the digitized +4.59e−4.

> **Reviewer note / response.** The review also states that κ moves the second zero
> from 0.28 to 0.34, "about a third of the unexplained log-gap at the peak". This
> session's independent re-run does not reproduce that part: the κ = 1 zeros come out
> at x = 0.0217 and 0.363 (§5.4's prototype paragraph says 0.03 and 0.28, which this
> re-run also does not reproduce), and κ moves the second one **down** to 0.348,
> slightly *away* from the digitized 0.457. The magnitude effect in the table above is
> confirmed twice over and is why the option is mandatory; the zero shift is small, of
> uncertain sign, and must be **measured by the implementer, not assumed**. G3a's
> windows (§5.4) are set wide enough to survive either answer.

**On "y_max = 2".** The sentence on CDKS p. 10 — *"the upper limit of their convolution
integral is y_max = 2"* — refers to **Ref. [4], Khan–Hoodbhoy**, not to CDKS's own
calculation; the earlier revision misattributed it. CDKS's own support is the y_max
above. Nor does the figure settle it: the digitized theory-1 table does **not** die at
x ≈ 1.4 — it is still +4.04e−6 at its last point, x = 1.590 — so "the curve dies at
1.4, consistent with 1.5" was not an argument and has been removed (§9 Q7).

### 1.4 The Clebsch–Gordan algebra (term 3), computed exactly

Coupling L = 2 with S = 1 to J = 1 (the α–d D-wave of ⁶Li, and identically the
deuteron's own D-wave), ⟨2 m_L : 1 m_S | 1 H⟩:

| H | (m_L, m_S) | CG | CG² |
|---|---|---|---|
| +1 | (0, +1) | +√(1/10) = +0.3162278 | 1/10 |
| +1 | (1, 0) | −√(3/10) = −0.5477226 | 3/10 |
| +1 | (2, −1) | +√(3/5) = +0.7745967 | 3/5 |
| 0 | (−1, +1) | +√(3/10) | 3/10 |
| 0 | (0, 0) | −√(2/5) = −0.6324555 | 2/5 |
| 0 | (1, −1) | +√(3/10) | 3/10 |
| −1 | (−2, +1) | +√(3/5) | 3/5 |
| −1 | (−1, 0) | −√(3/10) | 3/10 |
| −1 | (0, −1) | +√(1/10) | 1/10 |

With P_zz(m_S) = +1 for m_S = ±1 and −2 for m_S = 0 (Eq. 46 on a pure state):

* ⟨S_z⟩_{L=2} = Σ CG² m_S = **−1/2** ⟹ vector dilution 1 − (3/2)P_D. *(This is
  `ALPHA_D_VECTOR_POLARIZATION` in `beams.hpp`.)*
* ⟨P_zz⟩_{L=2} = Σ CG² P_zz(m_S) = **+1/10** ⟹ tensor dilution **1 − (9/10)P_D**.
  *(This is the closed form quoted next to `LI6_B1_RANK2_TRANSFER` in
  `constants.hpp:79`, which the tagged quadrature measures as 0.921947 at
  P_D = `P_D_LI6` = 0.0867.)*

So **term (3) is term (1) with the D-wave density and the weight 1/10** — the
"P₂-weighted CG sum" the recommendation asks for, and the analytic bridge between the
new backend and the constant the default uses today.

Writing it once more the way it must appear in the code comment: for ⁶Li in |1,H⟩,
b₁ picks up the deuteron's own tensor polarization weighted by
δ_T[⟨P_zz^d⟩] ≡ ⟨P_zz^d⟩_{H=0} − (⟨P_zz^d⟩_{H=+1} + ⟨P_zz^d⟩_{H=−1})/2, which is
**−3** for the S-wave (m_S = H, so ⟨P_zz⟩_H = P_zz(H)) and **−0.3** for the D-wave;
the ratio is 1/10 and the −3 is absorbed by the −1/3 in
Σ_m p(m)A_m = A_avg − (P_zz/3)b₁.

### 1.5 Recoupling ⁶Li, and why exactly four terms

⁶Li(1⁺) in the α–d basis, Eq. (19) verbatim with L = α–d relative orbital angular
momentum, S = deuteron spin, k⃗ = α–d relative momentum:

> **Φ^H(k⃗) = φ₀^{αd}(k) Y₀₀(k̂) χ^d_H + Σ_{m_L} ⟨2 m_L : 1 m_S | 1 H⟩ φ₂^{αd}(k) Y_{2m_L}(k̂) χ^d_{m_S},   m_S = H − m_L.**

The α is J = 0, so it has **no intrinsic b₁**. It is nevertheless a *constituent* of
the ⁶Li spectral function: CDKS Eq. (10) sums S(p) over constituents i, and one level
up that sum runs over {d, α}. The α carries −k⃗, and the H-dependence of its momentum
distribution is the **same**
δ_T ρ(k⃗) = [−(3/(4√2π))φ₀φ₂ + (3/(16π))|φ₂|²](3cos²θ − 1) as the deuteron's, because
(3cos²θ − 1) is even in k⃗. **So the α contributes to b₁ through the orbital alignment
of its own light-cone density.** This is the term the first revision of this design
missed.

Writing the deuteron sub-amplitude for a spin population p(m_S) as
**Σ_m p(m) A_m = A_avg − (P_zz^d/3) b₁ᵈ** (elementary, from
b₁ᵈ = A₀ − (A₊+A₋)/2 and P_zz = 1 − 3p₀), and taking δ_T over H, the ⁶Li b₁ separates
**exactly** into

* **δ_T of the spin-summed density × F₁ of the struck constituent** — this is CDKS
  Eq. (21) with the α–d overlaps, and it has **two** pieces:
  * **term (2d)** — struck deuteron, F₁ᵈ, z_d = (E_d − k_z κ)/M_d (§1.6);
  * **term (2α)** — struck α, F₁^α, z_α = (E_α + k_z κ)/M_α. The sign of k_z is
    immaterial: the integrand depends on cos θ only through cos²θ, and the δ maps
    c → −c, so f_α(z) is the same function either way and the *same*
    `ConvolutionKinematics` code serves both with the masses swapped;
* **−(1/3) δ_T[Σ_{m_S} P_zz(m_S) f^H_{m_S}] × b₁ᵈ** — struck deuteron only (b₁^α ≡ 0),
  whose pieces are
  * S-wave: δ_T = −3 |φ₀^{αd}Y₀₀|², isotropic ⟹ weight **+1** ⟹ **term (1)**;
  * D-wave: δ_T = −(9/16π)(9c⁴ − 8c² + 1)|φ₂^{αd}|², whose **angular average is
    −3/(40π)** ⟹ weight **+1/10** ⟹ **term (3)**;
  * S–D interference: the P_zz-weighted CG sum is
    P_zz(0)C⁰₀ − (C⁺¹₀ + C⁻¹₀)/2 = (−2)(−0.6324555) − (+0.3162278) = **+0.9486833**,
    so δ_T = **+0.9486833 × 2 φ₀φ₂ Y₀₀Y₂₀** ⟹ a term ∝ (3c²−1) φ₀φ₂ **b₁ᵈ** with
    coefficient exactly **1/3** of term (2d)'s SD structure (−0.9486833). *(The first
    revision wrote +1.8974 here, double-counting the explicit 2; only with 0.9486833
    does the stated 1/3 follow.)*

**Size of term (2α), computed for this revision.** With the kinematics of §1.6
(m_struck = M_α, m_recoil = M_d) the struck-α term is
(4/6) ∫(dz_α/z_α) δ_T f_α(z_α) F₁^{α, per nucleon}(x/z_α). Its P₂-weighted density
scales as 1/M², and it carries **twice** the counting factor of the d terms, so

  term (2α)/term (2d) ≈ [(4/6)/(2/6)] × (M_d/M_α)² = 2 (M_d/M_α)² = **0.506**,

and a direct quadrature over the VMC α–d waves with a common per-nucleon F₁ gives
**0.502 / 0.502 / 0.503 / 0.504 / 0.506 / 0.511 at x = 0.05 / 0.10 / 0.20 / 0.30 /
0.50 / 0.70, with the same sign as term (2d)**. It is *half of the term the whole
100 % band is built around*, not a truncation error, and it is cheap: in the default
F₁ it is literally the same isoscalar F₁ as term (2d) (§2.4, A9).

**The truncation is controlled, and this is the justification for "four terms".**
The two dropped pieces are
1. the anisotropic remainder of the D-wave b₁ᵈ weight,
   −(9/16π)[(72/35)P₄ − (4/21)P₂]|φ₂|². Its P₄ coefficient is ~15× the isotropic
   2/15 that term (3) keeps, so the first revision's *reasoning* ("≲ P_D × 1/10 ≈
   0.2 %") was wrong even where the number was small. It is small for a different,
   **x-dependent** reason: the P₂- and P₄-weighted densities integrate to zero
   (§5.1 G0e), so they enter only through the curvature of b₁ᵈ(x/z), and they
   multiply **b₁ᵈ, not F₁ᵈ** — crudely ~ (b₁ᵈ/F₁ᵈ) √P_D^{αd} × term (1).
   **Do not assert a number: compute it.** T14 (§6) convolves the P₂- and P₄-weighted
   D-wave densities with b₁ᵈ and records their ratio to term (1) at the §7 x points;
   that ratio, not 0.2 %, is what goes in the header comment and in
   `OPEN_ITEMS_SOLUTIONS.md` §10;
2. the S–D interference in the b₁ᵈ sector, which is term (2d)'s SD structure with
   coefficient 1/3 of it and **b₁ᵈ in place of F₁ᵈ**, i.e. suppressed by
   b₁ᵈ/F₁ᵈ ~ 10⁻³ — about **3 × 10⁻⁴ of term (2d)**.

Both are far inside the mandatory 100 % band. **They must nevertheless be written
down in the header comment**, because "why four and not six" is the first question a
reviewer will ask.

### 1.6 The α–d light-cone variables and the per-nucleon convention

**Definition (pin this; it is the single choice the whole model hangs on).** Mirroring
CDKS Eq. (18) one level up — y there is the nucleon fraction normalised to *its* fair
share 1/2 — define each constituent's fraction normalised to **its** fair share
(1/3 for the deuteron, 2/3 for the α):

> **z_d ≡ (M_{6Li}/M_d)·(p_d·q)/(P_{6Li}·q) ≃ 3 p_d⁻/P⁻_{6Li} = (p_d⁰ − k_z κ)/M_d ,
> p_d⁰ = M_d − ε_{αd} − k⃗²/(2M_α) ;**
> **z_α ≡ (M_{6Li}/M_α)·(p_α·q)/(P_{6Li}·q) ≃ (3/2) p_α⁻/P⁻_{6Li} = (p_α⁰ + k_z κ)/M_α ,
> p_α⁰ = M_α − ε_{αd} − k⃗²/(2M_d) .**

⟨z⟩ ≈ 1 for both (§5.2), **not** 1/3. The recommendation's "peaked near z = 1/3"
describes the *unnormalised* fraction p_d⁻/P⁻; this design uses the fair-share
normalisation because it makes the ⁶Li convolution **literally CDKS Eq. (16) with
N → d (resp. N → α)**, and because the 1/3 then appears exactly once, in the
per-nucleon factor. State this explicitly in the header so no one "fixes" it back.

With those z and the same δ-collapse as §1.3 — Jacobian M_struck/(k κ),
B ≡ M_struck(1−z) − ε_{αd}, r ≡ √(κ² + 2B/M_recoil),
k ∈ [M_recoil|κ−r|, M_recoil(κ+r)] — the **struck deuteron**
(M_struck = M_d, M_recoil = M_α) gives

  f_S(z) = (1/(2κ)) z M_d ∫ k dk |φ₀^{αd}|² ,
  f_D(z) = (1/(2κ)) z M_d ∫ k dk |φ₂^{αd}|² ,
  δ_T f_{αd}(z) = (3/(2κ)) z M_d ∫ k dk [ −φ₀^{αd}φ₂^{αd}/√2 + |φ₂^{αd}|²/4 ] (3c*² − 1) ,

and the **struck α** (M_struck = M_α, M_recoil = M_d) gives the same three integrals
with the masses swapped; only δ_T f_α is used (term 2α), because b₁^α ≡ 0:

  δ_T f_α(z) = (3/(2κ)) z M_α ∫ k dk [ −φ₀^{αd}φ₂^{αd}/√2 + |φ₂^{αd}|²/4 ] (3c*² − 1) ,
  c*(k, z) = [ M_α(1 − z) − ε_{αd} − k²/(2M_d) ] / (k κ) .

*(f_α uses the α's own recoil mass M_d in p_α⁰, per the parenthesis under CDKS
Eq. (12): the p⃗²/(2M) is really p⃗²/(2M_{A−i}).)*

**Per-nucleon factor.** Per nucleus,
b₁^{6Li} = Σ_i ∫dz′ δ_T f_i(z′) [structure function of i at x_A/z′] with z′ the raw
fraction; dividing by A = 6, converting to the fair-share z, and using
b₁^{d,per-deuteron} = 2 b₁^{d,per-nucleon} gives **2/6 = `LI6_B1_PER_NUCLEON`
multiplying the three deuteron terms** — the same constant the default `Li6B1` already
carries (`constants.hpp:83`; reuse it, do not define a new one) — and
**4/6 = 2 × `LI6_B1_PER_NUCLEON` multiplying the struck-α term**. The α term
therefore carries twice the counting weight of the d terms; that factor 2 is half of
why it is not negligible (§1.5).

### 1.7 The model, assembled

> **b₁^{6Li}(x,Q²)|_per nucleon =
> (2/6) ∫ (dz/z) { [ f_S(z) + w_{CG} f_D(z) ] b₁ᵈ(x/z, Q²)|_per nucleon
>                  + w_{αd} δ_T f_{αd}(z) F₁ᵈ(x/z, Q²)|_per nucleon }
> + (4/6) w_{αd} ∫ (dz/z) δ_T f_α(z) F₁^α(x/z, Q²)|_per nucleon**

with **w_CG = 1/10 exactly** (§1.4) and **w_{αd} = 1** nominal (the α–d D-wave weight
knob, which scales terms (2d) **and** (2α) together because they are one physical
effect, §4). b₂ = 2x b₁ (the `TensorSF` default); Δ = 0.

**Relation to what the default does today.** In the no-smearing limit
f_S → (1−P_D)δ(z−1), f_D → P_D δ(z−1), δ_Tf_{αd} → 0, δ_Tf_α → 0, the model collapses
to

  b₁^{6Li} = (2/6) [1 − (9/10)P_D^{αd}] b₁ᵈ(x) = `LI6_B1_PER_NUCLEON` ×
  `LI6_B1_RANK2_TRANSFER` × b₁ᵈ(x)  ⟸ **exactly `b1_li6_from_deuteron`** at
  P_D^{αd} = `P_D_LI6` = 0.0867.

So `LI6_B1_RANK2_TRANSFER = 0.921947` **encodes terms (1) + (3) with no smearing, no
α–d orbital term, no spectroscopic factor, at the scenario P_D**, and
`LI6_B1_PER_NUCLEON = 2/6` encodes the counting dilution. The new model changes four
things and nothing else:

| | default `Li6B1` | `Li6ConvolutionB1` | effect on the terms-(1)+(3) weight |
|---|---|---|---|
| P_D^{αd} | 0.0867 (Hulthén scenario) | 0.0193549 (VMC) | 0.921947 → 0.98258 |
| smearing in z | none | VMC α–d density | 0.98258 → 0.98304 (tiny; ⟨z⟩ ≈ 1) |
| N_{αd} | 1 (implicit) | 0.819481 | 0.98304 → **0.80555** |
| α–d orbital terms | absent | terms (2d) **and** (2α) | separate; (2α) ≈ 0.5 × (2d), same sign |

**Assumptions, every one of them.**
A1. Impulse approximation; ⁶Li = α ⊗ d in a single relative-motion wave function,
    plus a non-α–d remainder (1 − N_{αd} ≈ 18 %) whose b₁ is set to **zero**.
A2. Non-relativistic energy (CDKS Eq. 12) for each constituent inside ⁶Li, with the
    *other* cluster's mass as the recoil mass (M_α for the struck d, M_d for the
    struck α).
A3. Leading-twist Eq. (15) relation between b₁ and the helicity amplitudes; no b₃, b₄.
A4. No medium modification of b₁ᵈ or F₁ᵈ inside ⁶Li (CDKS make the same choice for
    F₁ᴺ inside the deuteron, their §III A).
A5. The α is J = 0, so **b₁^α ≡ 0** — but it is a constituent of the ⁶Li spectral
    function and its light-cone density carries the same (3cos²θ − 1) alignment as the
    deuteron's, so it **does** contribute to b₁, through term (2α) ≈ 0.5 × term (2d)
    (§1.5). *The first revision's "the α contributes to F₁ but not to b₁" was wrong;
    it is the reason this design has four terms, not three.*
A6. Angular average taken for everything multiplying b₁ᵈ (§1.5). The remainder is
    **measured** (T14), not asserted: the first revision's "≲ 0.2 %" rested on wrong
    reasoning.
A7. The deuteron's internal b₁ᵈ comes from a **separate** `TensorSF` (Miller or CDKS
    camp) and is *not* recomputed here: this design does not re-derive b₁ᵈ, it
    convolves whatever camp the caller injects.
A8. Q² is a spectator: b₁ᵈ, F₁ᵈ and F₁^α are evaluated at the *same* Q² as the ⁶Li
    kernel (no Q² rescaling with z). Standard in convolution models; note it.
A9. F₁^{α, per nucleon} is the **isoscalar nucleon F₁**, which in the default is
    numerically *identical* to F₁^{d, per nucleon} (both are (F₁ᵖ + F₁ⁿ)/2, because
    Z = N in both) — so the α's own EMC effect (≈ −10 % at x ≈ 0.6 for A = 4) is
    **not** modelled. It enters, if anyone wants it, through the `emc_ratio` hook of
    `NuclearF2` on a caller-supplied `alpha_f1`; it is a ≤ 10 % effect on a term that
    is itself inside the 100 % band.
A10. The δ-function is CDKS Eq. (17)'s κ = 1 form by default; Eq. (21)'s finite-|q⃗|
    form is `finite_q_delta` (§1.3) and is worth a factor 1.5–1.7 at large x. **The
    default is a choice, not a derivation**, and it is checklist item 0.

## 2. Data inputs

### 2.1 α–d magnitudes — `data/vmc/momenta/li6_ad1.momentum`

Second block (`RHOKA0, DRHOKA0, RHOKA2, DRHOKA2`), read with
`read_anl_momentum(path)[1]`, columns 0 and 1. K in fm⁻¹ (51 rows: 0.001, then
0.1 → 5.0 in steps of 0.1), ρ_L in fm³, **ρ_L = A_L²/(4π)** — i.e. the file gives
|φ_L|² up to normalisation and **carries no sign**. Header: *"6Li(1+) -- AV18+UX --
22-Mar-14, VMC 1M samples"* — this **is** the 2014 tabulation. *(The first revision
called it "2024"; the 12-Apr-24 date in `tagged.hpp`'s provenance note belongs to the
⁷Li file.)*

**One home for N_{αd}, and it is not a re-typed literal.** The file's own printed
normalisations `4*PI*TOTINT(RHOKA0*K**2:K)/(2*PI)**3 = 0.80362` and
`… RHOKA2 … = 0.015861` (lines 69–70) are the S- and D-wave spectroscopic factors, and
they are already the two numbers `VMC_P_D_LI6` is built from (`tagged.hpp:72`).
Therefore **add one constant next to it and rewrite the other in terms of it**:

```cpp
// tagged.hpp, replacing the current VMC_P_D_LI6 line
/// alpha-d spectroscopic factor from the momentum file's OWN printed norms.
inline constexpr double VMC_N_ALPHA_D_LI6 = 0.80362 + 0.015861;   // 0.819481
/// alpha-d D-state probability, VMC.
inline constexpr double VMC_P_D_LI6 = 0.015861 / VMC_N_ALPHA_D_LI6;
```

* **N_{αd} = `VMC_N_ALPHA_D_LI6` = 0.819481**, **P_D^{αd} = `VMC_P_D_LI6` =
  0.01935493**. These are the numbers the code uses (`norm_target`, §3). They are
  reachable at run time as `read_anl_momentum_norms(data_path(VMC_LI6_MOMENTUM))[1]
  + [2]`, which is the check T3 makes.
* The **trapezoid** recomputation of the same table gives 0.8036070 and 0.0158580,
  hence 0.8194650 and P_D = 0.0193516. These are *different numbers* — by 2.0e−5 and
  1.7e−4 relative — and they are the file's own quadrature spread, not an error.
  The first revision quoted "N_{αd} = 0.8194650 (file)" and "P_D = 0.0193516, which is
  `VMC_P_D_LI6` to all printed digits": the first is the trapezoid value presented as
  the file's, the second is simply false. **Neither appears in this design any more
  except as the tolerance T3 measures against.**
* `VMC_S_ALPHA_D_LI6 = 0.81971` (`tagged.hpp:74`) is the file's *total* block and is a
  **third** value. It is not used here, because P_D must be S₂/N with the same N;
  say so in the header comment so the next reader does not "unify" them.
* **Published cross-checks — three values, a 5 % spread, unexplained:**
  1. Wiringa, Schiavilla, Pieper, Carlson, PRC 89 (2014) 024305,
     [arXiv:1309.3794](https://arxiv.org/abs/1309.3794) §III, verified in the PDF:
     *"The integrated N_αd = 0.86 is a sum of S- and D-wave parts of 0.846 and 0.017"*
     ⟹ **0.863**, P_D = 0.0197;
  2. the 2004 overlap file `li6.ad` prints its own `ndx s-wave d-wave` =
     **0.856 0.838 0.017** (22-Apr-04, AV18+UIX);
  3. this momentum file (22-Mar-14, AV18+UX): **0.819481**, P_D = 0.019355.
  They span **5 %** in N_{αd} and 2 % in P_D. The spread is *not* a version
  difference this design can name (1 and 3 are the same year and the same
  Hamiltonian family), so it is carried as a **stated ±5 % systematic on terms (1)
  and (3)**, quoted next to the 100 % band in `OPEN_ITEMS_SOLUTIONS.md` §10 and
  reachable as a knob (`Li6ConvolutionOptions::norm_target`, §3). Default: value 3,
  the file the rest of the library already reads. T12 documents the spread; it does
  not explain it.

### 2.2 α–d S–D relative sign — `data/vmc/li6_alpha_d/li6.ad`

k-space block (`read_anl_overlap(path)[0]`), columns `Aad00(k)`, `Aad22(k)`, AV18+UIX
2004, signed amplitudes with MC errors. The existing machinery already does exactly
what is needed:

```
psi0 = vmc_from_momentum(momfile, 1, 0, 0, &vmc_from_overlap_k(ovfile, 0, 0));
psi2 = vmc_from_momentum(momfile, 1, 1, 2, &vmc_from_overlap_k(ovfile, 1, 2));
```
(the pattern in `src/core/tagged.cpp::li6_vmc_waves()` — **call that, do not
re-implement it**; and note that `li6_vmc_waves()` fixes the global phase **on the
pair**, `tagged.cpp:107-110`, which is why `ClusterPartialWave::from_vmc` must not
re-fix it per wave, §3). Magnitudes from `momenta/`, sign structure from
`overlap_old/`, per `docs/CONVENTIONS.md` and `docs/open_items/vmc_reconciliation.md`.

**Nodes, and the sign is k-dependent — do not describe it as uniform.**
S node at 0.678 fm⁻¹ = 0.134 GeV, D node at 2.25 fm⁻¹ = 0.444 GeV (from the file,
confirmed by minima of the momentum file's ρ₀/ρ₂, `vmc_reconciliation.md`). Raw file
phase: `Aad00 < 0` below the S node, `Aad22 > 0` below the D node. With the
(unobservable) global phase fixed to ψ₀(k→0) > 0, sign(ψ₂/ψ₀) is

| k [fm⁻¹] | region | sign(ψ₂/ψ₀) | φ₀φ₂ in the CDKS convention (φ₂ = −ψ₂) | fraction of ∫k²(ρ₀+ρ₂)dk |
|---|---|---|---|---|
| < 0.678 | below the S node | **−1** | **> 0** (opposite to the deuteron) | **66.4 %** |
| 0.678 – 2.25 | between the nodes | +1 | < 0 (deuteron-like) | 33.4 % |
| > 2.25 | above the D node | −1 | > 0 | 0.2 % |

So two thirds of the density has the "opposite to the deuteron" sign and one third
does not, and **δ_T f_{αd}(z) integrates across the node**: the header comment must
state the node structure, not "term (2) has the opposite sign to the deuteron's".
What survives is the *net* sign of the term, which the code measures (§7) and which
comes out opposite over most of the gate window because the low-k region dominates.

**Consequence, and the sign gate.** The α–d SD interference in Eq. (21),
−(3/(4√2π))φ₀φ₂, is therefore **negative at low k for α–d and positive for the
deuteron**. This is the microscopic statement of the ⁶Li quadrupole puzzle
(§4.3 for the primary sources of Q(⁶Li) and Q_d) and it is *checkable*: the α–d
contribution to the ⁶Li quadrupole moment, Q_{αd} ∝ −(√2/5)⟨r²⟩_{02} built from the
same φ₀, φ₂, **must come out negative**. Make that a test (§6, T10), together with
**both** signs of φ₀φ₂ — positive at k = 0.2 fm⁻¹ and negative at k = 1.0 fm⁻¹ — so
the node crossing itself is pinned. If the code accidentally uses the ANL amplitudes
as φ₀, φ₂ directly (i.e. with the iᴸ convention mismatched), terms (2d) and (2α) flip
and the whole ⁶Li answer moves by ~100 %.

### 2.3 Deuteron u(k), w(k) for the A = 2 gate — `data/vmc/deuteron/fdeut.av18`

Header block: `ebind dstate qm mu as eta rd r4` → `2.224574 0.057599 …`, so the file
carries **its own ε_d = 2.224574 MeV** and D-state probability **0.057599**; the next
block's `ti = 19.810` MeV is the kinetic energy. The k-space block begins at the line
`   k          u(k)                w(k)` (line 10221; 201 rows, k = 0 → 20 fm⁻¹, step
0.1). Verified for this design: **∫k²[u(k)² + w(k)²]dk = 0.999976** and
∫k²w²dk/∫k²(u²+w²)dk = **0.0575999** — i.e. the file is already in CDKS's
normalisation ∫dp p²[|φ₀|²+|φ₂|²] = 1 with φ₀ = u(k), φ₂ = −w(k). No rescaling
needed beyond the baryon-number renormalisation of §3.

There is **no parser for this file in C++ yet** — only `validation/vmc_overlap_prototype.py`
and `validation/vmc_reconcile.py` have `parse_fdeut`. The kernel agent must add
`read_fdeut_k(path) -> {k_gev, u, w, ebind_gev}` next to the other ANL readers in
`src/core/cluster.cpp` (same `split_numbers` machinery; the block is plain 3-column).
It returns the header's `ebind` so the gate's ε_d has **one** source and no literal is
re-typed; `DEUTERON_P_TAG().separation_energy` = 2.2246e-3 is the same number rounded
and the 2.6e-5 relative difference is recorded in the gate report.

**Units.** `read_fdeut_k` converts k from fm⁻¹ to GeV with **`HBARC_GEV_FM`**
(`constants.hpp:51`, 0.19733), which `docs/CONVENTIONS.md` names as *the* conversion,
and so does everything in `b1_nuclear.cpp` (`alpha_d_quadrupole_fm2`, the k = 0.2 and
k = 1.0 fm⁻¹ sign checks). `cluster.cpp` has a **file-local** `kHbarCGeVfm =
0.1973269804` (line 23) that the existing ANL readers use and that disagrees with
`HBARC_GEV_FM` in the 5th digit; `read_fdeut_k` lives in that file, so it will pick up
the file-local one unless told otherwise. **Decide once and say so in the reader's
comment**: use `HBARC_GEV_FM` for the new reader and note the 1.5e−5 relative
inconsistency with its neighbours as a pre-existing item (it is 20 keV at k = 5 fm⁻¹,
far inside everything here), so the §8 review grep has something to find.

**Wave-function caveat for the gate.** CDKS used **CD-Bonn** (P_D = 4.85 %); the repo
has **AV18** (P_D = 5.76 %). That is a 19 % difference in the D-state probability and
a larger one in the high-k tail, which CDKS explicitly say "play an important role"
(§IV). The gate's tolerance must absorb it (§5).

### 2.4 Nucleon / deuteron / α structure functions

* **b₁ᵈ per nucleon** — injected as a `std::shared_ptr<const TensorSF>`, exactly as
  `Li6B1` already takes one. The default is **the raw digitized CDKS theory-1 column**,
  i.e. `CdksB1` built with `per_nucleon_table = true` (§9 Q1: CDKS Eqs. (10)/(16) and
  Fig. 6 establish that the digitized curve is **already per nucleon**, so
  `b1_convolution()`'s extra `B1_PER_DEUTERON_TO_PER_NUCLEON = 0.5` would halve it a
  second time). `MillerB1` and the default `CdksB1` stay reachable, and the choice is
  printed in the plot legend, because Miller and CDKS are different *camps*.
* **F₁ᵈ per nucleon** — CDKS Eq. (22), spelled **once and unambiguously**:

  ```cpp
  f1_cdks(NuclearF2(DEUTERON(), unpol, nullptr, r_func).f2a(x, q2) / 2.0, x, q2, r_func)
  ```

  `NuclearF2::f1a` is the **massless** F₂/(2x(1+R)) and must never be called from
  `b1_nuclear.cpp`; the two differ by 1 + γ² = 1.35 at x = 0.5, Q² = 2.5, which is the
  very factor §1.2 calls load-bearing. `unpol` null ⇒ one `ToyF2` shared with the
  kernel (§4.2); `r_func` null ⇒ `r_sigma_lt`.
* **F₁^α per nucleon** — the **same expression with `Ion{"4He", 4, 2, …}`**, which
  because Z = N is numerically identical to the deuteron's:
  (F₁ᵖ + F₁ⁿ)/2. Implement it as the isoscalar directly (no new `Ion`, no
  `nuclear_mass` lookup — `NuclearF2::f2a` uses only Z and N) and say in the comment
  that the α's own EMC effect is **not** modelled (A9); `Li6ConvolutionOptions::alpha_f1`
  overrides it for anyone who wants an EMC-corrected α.
* For the **A = 2 gate** the F₁ᴺ slot is instead **(F₁ᵖ + F₁ⁿ)/2 with R = `r1998`**,
  matching CDKS. F₁ᵈ/F₁ᴺ differs from 1 by ≤ 1 % below x = 0.6, so this choice is not
  where the uncertainty lives — but spell it in the header anyway.

## 3. API

New pair of files, **`include/lipolgen/b1_nuclear.hpp` + `src/core/b1_nuclear.cpp`**
(the name the run PLAN Phase D already uses). `sf.hpp` is not touched except for one
`#include`-free forward mention in a comment; `cluster.hpp`/`cluster.cpp` gain only
`read_fdeut_k`; `tagged.hpp` gains `VMC_N_ALPHA_D_LI6` and rewrites `VMC_P_D_LI6` in
terms of it (§2.1). All numerics use `numerics.hpp` (`trapezoid`, `linspace`,
`pairwise_sum`) — no new quadrature library, and no RNG anywhere. **`np_interp` is
NOT used for the wave functions**: see "Where the numerics live".

```cpp
#ifndef LIPOLGEN_B1_NUCLEAR_HPP
#define LIPOLGEN_B1_NUCLEAR_HPP

/// \file b1_nuclear.hpp
/// b1 of a two-cluster nucleus by convolution: Cosyn-Dong-Kumano-Sargsian
/// PRD 95 (2017) 074036 Eqs. (16), (17), (21), one level up.  See
/// docs/open_items/run_2026-09-02/design_D_b1_li6.md for the derivation,
/// the symbol table and the FOUR-term truncation argument (the struck-alpha
/// orbital term is NOT optional: it is ~0.5 x the struck-deuteron one).

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "lipolgen/cluster.hpp"   // VmcRadial, read_anl_*, read_fdeut_k, data_path
#include "lipolgen/sf.hpp"        // TensorSF, UnpolSF, RFunc, SFFunc3

namespace lipolgen {

/// F1 of CDKS Eq. (22): (1 + gamma^2) F2 / (2x (1+R)).  NOT
/// `UnpolSF::f1_from_f2` and NOT `NuclearF2::f1a`, which are the massless
/// form -- the target-mass factor is `gamma_squared(x, q2)` (asymmetries.hpp)
/// and is 0.90 at x = 0.8, Q2 = 2.5.  The convolution kernel MUST be fed this
/// one; `f1a` must not appear anywhere in b1_nuclear.cpp.
double f1_cdks(double f2, double x, double q2, const RFunc& r_func = nullptr);

/// One partial wave of a two-cluster relative wave function in the CDKS
/// convention: phi_L(k) with the i^L factor, so phi_2 < 0 wherever the
/// "no-i^L" amplitude W(k) > 0.
///
/// CONTRACT of `from_vmc`: it applies ONLY phi_L = i^L psi_L (i.e. phi_2 =
/// -psi_2, phi_0 = psi_0).  It does NOT touch the global phase, because the
/// global phase is a property of the PAIR and is already fixed on the pair by
/// `li6_vmc_waves()` (tagged.cpp:107-110, psi_0(k->0) > 0).  A `from_vmc`
/// that re-fixed the sign per wave (e.g. forcing psi.front() > 0 on the D
/// wave) would flip the S-D relative sign and hand the sign gate the wrong
/// answer.  Callers therefore pass `li6_vmc_waves()[0].vmc` and `[1].vmc`;
/// the pair overload below is the recommended entry point because it cannot
/// be called wrongly.
struct ClusterPartialWave {
  int l = 0;
  std::vector<double> k;    ///< GeV, strictly increasing
  std::vector<double> phi;  ///< CDKS-convention amplitude, signed
  std::string provenance;

  /// Cubic-spline interpolation in k (see "Where the numerics live"),
  /// and 0 OUTSIDE the table -- matching `VmcRadial::operator()`
  /// (cluster.cpp:231), NOT the end-clamped value `np_interp` returns.
  double operator()(double kk) const;

  static ClusterPartialWave from_vmc(const VmcRadial& v, int l);
  /// Recommended: takes the pair whose global phase is already fixed.
  static std::pair<ClusterPartialWave, ClusterPartialWave>
      from_vmc_pair(const VmcRadial& s, const VmcRadial& d);
  static ClusterPartialWave from_uw(std::vector<double> k_gev,
                                    std::vector<double> uw, int l);
};

/// Kinematics of one struck constituent in the CDKS Eq. (12) reduction.
/// For 6Li there are TWO of these: {M_d, M_alpha} for the struck deuteron
/// (terms 1, 2d, 3) and {M_alpha, M_d} for the struck alpha (term 2a).
struct ConvolutionKinematics {
  double m_struck = 0.0;    ///< M_d for the struck d, M_alpha for the struck a
  double m_recoil = 0.0;    ///< the OTHER cluster: enters p^0 = M - eps - k^2/(2 m_recoil)
  double separation = 0.0;  ///< epsilon [GeV]
  /// |q|/nu = sqrt(1 + gamma^2) of CDKS Eq. (21)'s delta-function.  1.0 is
  /// Eq. (17) as printed (the DEFAULT); `Li6ConvolutionOptions::finite_q_delta`
  /// makes it kappa(x, Q2) and then the densities are rebuilt per x.
  double kappa = 1.0;

  /// Allowed |k| interval at fixed y where |cos theta*| <= 1 (design 1.3):
  /// [m_recoil |kappa - r|, m_recoil (kappa + r)],
  /// r = sqrt(kappa^2 + 2B/m_recoil), B = m_struck (1 - y) - epsilon.
  /// Returns false when B < -m_recoil kappa^2 / 2.
  bool k_range(double y, double* lo, double* hi) const;
  double cos_star(double k, double y) const;   ///< [.. - k^2/(2 m_recoil)]/(k kappa)
  double y_max() const;     ///< 1 - eps/m_struck + kappa^2/2   (kappa-DEPENDENT)
};

/// The light-cone densities of CDKS Eqs. (17) and (21) for ONE struck
/// constituent, tabulated on a y grid and interpolated.  Deterministic; no
/// RNG; safe to share.  6Li holds two of these.
class LightConeDensities {
 public:
  struct Options {
    /// y-grid layout, spelled out so two implementers get the same bits.
    /// Three segments, trapezoid on each, breakpoints included in both.
    std::size_t n_y_low = 600;    ///< [1e-4, y_mid_lo]
    std::size_t n_y_mid = 2400;   ///< (y_mid_lo, y_mid_hi]  -- delta_T f's structure
    std::size_t n_y_high = 800;   ///< (y_mid_hi, y_max]
    double y_mid_lo = 0.85;
    double y_mid_hi = 1.15;
    std::size_t n_k = 2001;       ///< k points per y, between the EXACT endpoints
    /// Renormalise by the baryon-number condition int f dy = `norm_target`
    /// (CDKS below Eq. 21).  `norm_target` = 1 for the deuteron,
    /// VMC_N_ALPHA_D_LI6 for the alpha-d cluster (the spectroscopic factor is
    /// REAL suppression, design 1.7 A1).  `renormalize = false` exposes the
    /// raw quadrature, which is what G1c / T3 / T4 test.
    bool renormalize = true;
    double norm_target = 1.0;
  };

  LightConeDensities(ClusterPartialWave phi0, ClusterPartialWave phi2,
                     ConvolutionKinematics kin, Options opt = {});
  /// Test-only: the no-smearing limit, f_S -> (1-p_d) delta(z-1),
  /// f_D -> p_d delta(z-1), delta_T f -> 0.  Flagged internally so that
  /// `convolve` short-circuits to `g(x)` instead of quadrature (a delta
  /// cannot live on the trapezoid grid).  Used by T6, which therefore tests
  /// the WEIGHT ASSEMBLY only -- the kernel is covered by T2/T5.
  static LightConeDensities delta_limit(double p_d, double norm);

  double f_s(double y) const;        ///< S-wave, Eq. (17) with |phi_0|^2
  double f_d(double y) const;        ///< D-wave, Eq. (17) with |phi_2|^2
  double f_unpol(double y) const;    ///< f_s + f_d
  double delta_t_f(double y) const;  ///< Eq. (21), SD + DD
  double delta_t_f_sd(double y) const;
  double delta_t_f_dd(double y) const;
  /// P2- and P4-weighted D-wave densities of design 1.5 truncation item 1,
  /// used by T14 to MEASURE the angular-average error instead of asserting it.
  double f_d_p2(double y) const;
  double f_d_p4(double y) const;

  double norm() const;               ///< int f_unpol dy
  double mean_y() const;             ///< int y f dy / int f dy  ==  <y^2>_phi/<y>_phi
  double p_d() const;                ///< int f_d dy / norm()   (y-WEIGHTED, see T3)
  /// Unweighted momentum-space ratio int k^2|phi_2|^2 / int k^2(|phi_0|^2+|phi_2|^2):
  /// THIS is the quantity that equals the file's D/(S+D), and the one T3 pins
  /// to VMC_P_D_LI6.
  double p_d_momentum() const;
  const std::vector<double>& y_grid() const;

  /// int (dy/y) dens(y) g(x/y), trapezoid on the stored grid refined to `n`
  /// points.  The SUPPORT OF `g` IS `g`'S BUSINESS: this routine does not cut
  /// at x/y = 1.  A free-nucleon F1 returns 0 above 1 on its own; an injected
  /// b1_d table is defined out to its own x_max() (kB1CdksQ2p5 runs to 1.59,
  /// and the deuteron's per-nucleon support genuinely exceeds 1, design 1.3),
  /// and cutting it at 1 would silently drop that strength from terms (1)/(3).
  double convolve(const std::function<double(double)>& dens_at,
                  const std::function<double(double)>& g, double x,
                  std::size_t n = 4001) const;
 private:
  /* grids + cached tables */
};

/// Which deuteron b1 camp and which alpha-d inputs a `Li6ConvolutionB1` is
/// built from, plus the knobs the 100 % band is expressed with.
struct Li6ConvolutionOptions {
  /// b1 of the EMBEDDED deuteron, per nucleon.  Null => the RAW digitized
  /// CDKS theory-1 column (CdksB1 with per_nucleon_table = true), because
  /// CDKS Eqs. (10)/(16) and their Fig. 6 say the published curve is already
  /// per nucleon (design 9 Q1).  MillerB1 is a different camp -- mixing them
  /// is legal but must be said out loud in any plot legend.
  std::shared_ptr<const TensorSF> deuteron_b1;
  /// F1 of the embedded deuteron, per nucleon.  Null =>
  ///   f1_cdks(NuclearF2(DEUTERON(), unpol, nullptr, r_func).f2a(x,q2)/2,
  ///           x, q2, r_func)
  /// -- CDKS Eq. (22), NOT NuclearF2::f1a (design 2.4).
  std::shared_ptr<const UnpolSF> unpol;
  /// F1 of the alpha, per nucleon, for term (2a).  Null => the isoscalar
  /// (F1p + F1n)/2 built the same way, i.e. numerically identical to the
  /// deuteron's; no alpha EMC effect (design 1.7 A9).
  SFFunc3 alpha_f1;
  RFunc r_func;                        ///< null => r_sigma_lt

  /// TERM KNOBS -- each multiplies one term of design 1.7.  All 1 = nominal.
  double w_embedded_s = 1.0;           ///< term (1)
  double w_alpha_d_dwave = 1.0;        ///< terms (2d) AND (2a) <- THE band knob
  double w_cg_dwave = 1.0;             ///< term (3), on top of the exact 1/10

  /// Normalisation.  `use_spectroscopic_factor` = true (default) keeps
  /// `norm_target` = VMC_N_ALPHA_D_LI6 = 0.819481, i.e. the non-alpha-d 18 %
  /// of 6Li is given b1 = 0.  false renormalises the density to 1 (the "all
  /// of 6Li is alpha+d" variant); the two differ by 1/N_ad = 1.22 on terms
  /// (1), (2d), (2a) and (3).  `norm_target` overrides the value, which is
  /// how the +-5 % N_ad systematic of design 2.1 is quoted.
  bool use_spectroscopic_factor = true;
  double norm_target = 0.0;            ///< 0 => VMC_N_ALPHA_D_LI6

  /// CDKS Eq. (21)'s finite-|q| delta-function instead of Eq. (17)'s
  /// kappa = 1 (design 1.3).  FALSE by default.  true rebuilds both density
  /// sets per x with kappa = sqrt(1 + gamma^2(x, q2)); worth a factor
  /// 1.5-1.7 at x >= 0.5 and it is checklist item 0 of design 5.4.
  bool finite_q_delta = false;

  /// Quadrature.
  LightConeDensities::Options quad;
};

/// b1 of 6Li per nucleon, four-term alpha-d convolution.  OPT-IN; the
/// library default stays `Li6B1(MillerB1)` (pipeline.cpp:344).
class Li6ConvolutionB1 : public TensorSF {
 public:
  explicit Li6ConvolutionB1(Li6ConvolutionOptions opt = {});
  /// Test-only: build on prepared densities (e.g. `delta_limit`), T6.
  Li6ConvolutionB1(Li6ConvolutionOptions opt, LightConeDensities d_struck,
                   LightConeDensities a_struck);

  double b1(double x, double q2, double f1) const override;
  /// b2 = 2 x b1 (base class), delta = 0 (base class).

  /// The four terms separately, same units and conventions as `b1`.
  double b1_embedded_s(double x, double q2, double f1) const;        ///< (1)
  double b1_alpha_d_dwave_d(double x, double q2, double f1) const;   ///< (2d)
  double b1_alpha_d_dwave_alpha(double x, double q2, double f1) const;///< (2a)
  double b1_cg_dwave(double x, double q2, double f1) const;          ///< (3)
  /// (2d) + (2a), the physical "alpha-d orbital" term.
  double b1_alpha_d_dwave(double x, double q2, double f1) const;
  /// SD / DD split of (2d) and (2a), for the Fig.-4-style plot.
  double b1_alpha_d_sd(double x, double q2, double f1) const;
  double b1_alpha_d_dd(double x, double q2, double f1) const;
  /// P2/P4 remainders of design 1.5 truncation item 1 -- REPORTED, not added
  /// to `b1`; T14 records their ratio to term (1).
  double b1_cg_dwave_p2_remainder(double x, double q2, double f1) const;
  double b1_cg_dwave_p4_remainder(double x, double q2, double f1) const;

  /// A scaled copy for the mandatory 100 % band: `scale` in {0, 1, 2}
  /// multiplies the WHOLE b1 (design 4.3).  Cheap: shares the densities.
  std::shared_ptr<const TensorSF> banded(double scale) const;

  const LightConeDensities& densities() const;        ///< struck deuteron
  const LightConeDensities& densities_alpha() const;  ///< struck alpha
  const Li6ConvolutionOptions& options() const;
  /// int b1 dx over [lo, hi] -- the Close-Kumano number, REPORTED not
  /// enforced, exactly like `close_kumano_integral`.
  double close_kumano_integral(double lo = 0.01, double hi = 1.2,
                               double q2 = 2.5, std::size_t n = 241) const;
 private:
  /* opt_, densities_d_, densities_a_, cached b1_d / F1_d / F1_alpha closures */
};

/// A = 2 VALIDATION GATE.  The same kernel fed the AV18 deuteron u(k), w(k)
/// and a nucleon F1: this MUST reproduce CDKS Fig. 4 (tables::kB1CdksQ2p5)
/// before any 6Li number is quoted.  See design section 5.
class DeuteronConvolutionB1 : public TensorSF {
 public:
  struct Options {
    std::string fdeut_path = data_path("vmc/deuteron/fdeut.av18");
    std::shared_ptr<const UnpolSF> unpol;   ///< null => ToyF2
    RFunc r_func;                           ///< null => r1998 (CDKS's choice)
    bool target_mass = true;                ///< CDKS Eq. (22) factor
    bool finite_q_delta = false;            ///< design 1.3, checklist item 0
    LightConeDensities::Options quad;
  };
  explicit DeuteronConvolutionB1(Options opt = {});
  double b1(double x, double q2, double f1) const override;
  double b1_sd(double x, double q2, double f1) const;
  double b1_dd(double x, double q2, double f1) const;
  const LightConeDensities& densities() const;
};

/// The A = 2 gate's landmarks, computed (never typed): the two sign changes,
/// the extremum positions and values, and int b1 dx.  Both the doctest and
/// `validation/b1_li6_table.py` call this, and it is also applied to the
/// digitized table itself so G3a/G3b compare like with like.
struct B1Landmarks {
  std::vector<double> zeros;   ///< ascending, with the sign of the crossing
  std::vector<int> zero_slope; ///< +1 rising, -1 falling
  double x_min = 0.0, xb1_min = 0.0;
  double x_max = 0.0, xb1_max = 0.0;
  double integral_b1 = 0.0;
};
B1Landmarks b1_landmarks(const std::function<double(double)>& xb1,
                         const std::vector<double>& x_grid);
B1Landmarks b1_landmarks_of_table();   ///< tables::kB1CdksQ2p5(), raw column

/// alpha-d contribution to the 6Li quadrupole moment from the SAME phi_0,
/// phi_2 -- the SIGN GATE of design 2.2.  Must come out NEGATIVE.  fm <-> GeV
/// with `HBARC_GEV_FM` (design 2.3).
double alpha_d_quadrupole_fm2(const ClusterPartialWave& phi0,
                              const ClusterPartialWave& phi2);

}  // namespace lipolgen
#endif
```

**Where the numerics live.** `LightConeDensities` owns *all* quadrature:
* the k integral at fixed y, on a **uniform grid between the two exact endpoints**
  `k_range(y)` (never a mask over a fixed grid — that costs 2 % on ∫f dy, measured),
  with the upper endpoint clipped to `min(hi, phi.k.back())` so the `n_k` points are
  not spent on the zero region beyond the table edge (0.987 GeV for the α–d file,
  3.95 GeV for `fdeut`);
* the y grid, **non-uniform and fully specified by `Options`**:
  `n_y_low` points on [10⁻⁴, `y_mid_lo`], `n_y_mid` on (`y_mid_lo`, `y_mid_hi`],
  `n_y_high` on (`y_mid_hi`, `y_max()`] — dense where δ_T f has a factor-10 structure.
  The prototype's 600 / 2400 / 800 with breakpoints 0.85 / 1.15 are the defaults, and
  the §7 table and the T2 pins are quoted **for those values**;
* the convolution `∫(dy/y) dens(y) g(x/y)` on a refined uniform y grid from
  max(x, y_min) to `y_max()`, with the support of `g` left to `g` (see the header).

**Interpolation of the wave functions is a cubic spline, not `np_interp`.** This is a
correctness item, not a style one: `fdeut.av18` is tabulated on a coarse 0.1 fm⁻¹ grid,
and linear interpolation of u, w biases ∫k²(u²+w²)dk by **+2.2 %**
(1.0218 linear vs 1.0003 cubic on a fine grid, measured for this revision) — which is
exactly the discrepancy between the first revision's raw ∫f dy = 1.0088 and the exact
0.98707 (1.0088/0.98707 = 1.0220). The baryon-number renormalisation then *hides* that
2.2 % error, and no grid-refinement test can see it because the **table spacing** is the
problem. Spline u(k), w(k) and √ρ_L(k) (or interpolate in k²); `np_interp` stays fine
for the y-grid lookups of the already-smooth densities.

Everything is `const` after construction; no mutable state, no RNG, no thread-local
anything — determinism is structural.

## 4. Wiring

### 4.1 On which channel the model is legal

**Inclusive only.** `Li6ConvolutionB1` models the b₁ of the **whole ⁶Li as seen
inclusively**: the z-smearing of terms (1)/(3) and the orbital alignment of terms
(2d)/(2α) are both properties of the α–d relative motion.

In the **tagged** channels (`TaggedLi6Alpha`, `TaggedLi7Alpha`, `TaggedDeuteronP`) the
sampler already draws the α–d momentum **and** its m-dependent angular correlation
from n_M(k, k̂) (`tagged.hpp`), so the same α–d wave function is in the event weight
already. That is exactly the double counting `StruckClusterOptions::inclusive_b1`
guards against (`pipeline.hpp:176-186`: *"in the impulse approximation the embedded
deuteron's b₁ IS the k-integral of the m-dependent spectator density"*). Putting
`Li6ConvolutionB1` in a struck-cluster kernel would count it **three** times (the
sampler's k̂ correlation, term (2)'s orbital alignment, and the z-smearing of the same
density).

**How that is enforced — decided, because the first revision specified something that
cannot be written.** `struck_cluster_kernel(const TaggedChannel&, const
StruckClusterOptions&)` (`pipeline.hpp:202`) takes **no `B1Model`**, so there is
nothing there to reject; the first revision's "`struck_cluster_kernel` must reject a
`B1Model::Li6Convolution`" and T11's first half could not compile. **Resolution: do
not add a `B1Model` to `StruckClusterOptions`.** The guard is

* `PipelineConfig::validate()` (below), which is the only route from the CLI, and
* the fact that `StruckClusterOptions::inclusive_b1` installs `toy_b1` and nothing
  else — a comment in `struck_cluster_kernel` says why, pointing at this paragraph.

A caller who hand-builds an `InclusiveKernel` with a `Li6ConvolutionB1` b1_func and
feeds it to a tagged sampler is out of the library's reach, exactly as they are today
for any other hand-built kernel. T11 therefore tests `validate()` alone.

On the **coherent** channel the tensor signal lives in the recoil azimuth
(`pipeline.cpp:571-578`), so the inclusive b₁ is not folded in at all — nothing to do.

### 4.2 `PipelineConfig` field and CLI flag

```cpp
// pipeline.hpp, next to TritonSfChoice
/// Which b1 backend the INCLUSIVE kernel's rank-2 slot is filled with when
/// `PipelineConfig::kernel` is null (CLI `--b1-model`).
///   Miller        the DEFAULT, bit-for-bit what every published number used:
///                 Li6B1(MillerB1) through LI6_B1_RANK2_TRANSFER (2/6).
///                 This IS the run PLAN's "toy": today's default b1_func is
///                 `Li6B1(MillerB1)`, which is what `toy_b1` reaches.
///   Cdks          the same rank-2 transfer on the CDKS convolution camp --
///                 |b1| ~ 100x smaller, sign structure different.
///   Li6Convolution  b1_nuclear.hpp's four-term alpha-d convolution.
///                 6Li ONLY (see validate(); a spin-1 test is not enough --
///                 the deuteron is spin-1 too).
enum class B1Model : std::uint8_t { Miller, Cdks, Li6Convolution };

// in PipelineConfig
B1Model b1_model = B1Model::Miller;          ///< default unchanged
double  b1_band_scale = 1.0;                 ///< 0 / 1 / 2, design 4.3
double  b1_alpha_d_dwave_weight = 1.0;       ///< terms (2d)+(2a) knob
```

**`validate()` keys on the ISOTOPE, not the spin.** The inclusive channel accepts
`isotope == "d"`, which is spin-1, so a spin test would let a deuteron beam run with
`--b1-model li6-convolution` and silently apply the ⁶Li α–d convolution (N_{αd}, α–d
densities, the 2/6 and 4/6 factors) to a deuteron. The rule is

* `b1_model == Li6Convolution` requires `channel == Inclusive` **and**
  `isotope == "6Li"`; `"d"` and `"7Li"` throw, each with its own message
  (`"d"`: use `DeuteronConvolutionB1` directly, that is the A = 2 gate;
  `"7Li"`: spin 3/2 has no rank-2 input here);
* `cfg.kernel != nullptr` **and** `b1_model != Miller` also throws, because a
  caller-supplied kernel silently wins over the flag today and there is no rule in
  `validate()` at all. Throwing is the safe reading: the two say contradictory things.

*(`default_inclusive_kernel` keys on `spin == 1` today, which has the same latent
oddity for Miller. Out of scope here — but note it in the commit message.)*

`default_inclusive_kernel(const Ion&)` gains an overload
`default_inclusive_kernel(const Ion&, B1Model, double band_scale, double w_alpha_d)`;
the one-argument form keeps calling it with `(Miller, 1.0, 1.0)` so
`examples/generate_inclusive.cpp` and every existing caller are untouched.
`pipeline.cpp:556` becomes
`cfg_.kernel ? cfg_.kernel : default_inclusive_kernel(ion, cfg_.b1_model, cfg_.b1_band_scale, cfg_.b1_alpha_d_dwave_weight)`.

**Spell the wiring, because copying the existing pattern gives either a stale kernel
or a use-after-free.** `TensorSF::b1_func()` returns a closure that captures raw
`this` and "must not outlive the backend" (`sf.hpp:370-371`); today's
`default_inclusive_kernel` keeps the `Li6B1` alive with a function-local `static`. A
`static` is wrong for the new branches — it would freeze the first band scale for the
whole process, so `banded(0)` after `banded(1)` returns the wrong object — but
without it the closure dangles as soon as the local `shared_ptr` dies. So:

```cpp
std::shared_ptr<const InclusiveKernel> default_inclusive_kernel(
    const Ion& ion, B1Model model, double band_scale, double w_alpha_d) {
  InclusiveKernel::Options opt;
  if (std::fabs(ion.spin - 1.0) < 1e-9) {
    if (model == B1Model::Miller) {
      static const auto li6_b1 =            // UNCHANGED default path
          std::make_shared<Li6B1>(std::make_shared<MillerB1>());
      opt.b1_func = li6_b1->b1_func();
    } else if (model == B1Model::Cdks) {
      auto b = std::make_shared<const Li6B1>(std::make_shared<CdksB1>());
      opt.b1_func = [b, band_scale](double x, double q2, double f1) {
        return band_scale * b->b1(x, q2, f1);   // capture the shared_ptr
      };
    } else {
      Li6ConvolutionOptions o;
      o.w_alpha_d_dwave = w_alpha_d;
      auto b = std::make_shared<const Li6ConvolutionB1>(std::move(o));
      opt.b1_func = [b, band_scale](double x, double q2, double f1) {
        return band_scale * b->b1(x, q2, f1);
      };
    }
    // The Delta slot is the SAME for every B1Model -- it is not part of this
    // change and must not be dropped on the new branches.
    opt.delta_func = [](double x, double q2, double f1) {
      return toy_delta_gluon(x, q2, f1, 1e-2);
    };
  }
  // ONE ToyF2, shared: the core must not link LHAPDF (sf.hpp).  LhapdfSF only
  // ever enters through a caller-supplied cfg.kernel or the validation script.
  auto f2 = std::make_shared<ToyF2>();
  opt.f2_source = f2;                      // (as today: null would also mean ToyF2)
  return std::make_shared<InclusiveKernel>(ion, opt);
}
```
and the same `f2` is handed to `Li6ConvolutionOptions::unpol`, so "the kernel's
`UnpolSF`" is literally one object and not a second `ToyF2` that merely looks like it.

CLI (`python/lipolgen/cli.py`, mirroring `--triton-sf`/`--cluster-wave`):

```
--b1-model {miller,cdks,li6-convolution}   default miller
--b1-band-scale FLOAT                      default 1.0; run 0 and 2 as the
                                           mandatory 100 % band
--b1-alpha-d-dwave FLOAT                   default 1.0; terms (2d)+(2a) knob
```
plus `B1_MODELS = {"miller": …, "cdks": …, "li6-convolution": …}` in
`python/lipolgen/__init__.py`, three new `make_config` keyword arguments and three
`DEFAULTS` entries — **all three are mandatory**, because `cli.resolve()` rejects any
config-file key that is not in `DEFAULTS` (`cli.py:168`). Help text must carry the
sentence *"inclusive channel only, ⁶Li only; on a tagged channel the α–d density is
already in the event weight"*.

**Provenance in the outputs.** The band is produced by re-running with
`--b1-band-scale 0/1/2`, but neither the npz/HFS `meta` dict
(`bindings.cpp:387-406`: channel, isotope, seed, optics …) nor the CLI summary records
which b₁ backend or which band scale produced a file, so three otherwise identical
runs are indistinguishable — against this design's own "never quote a single row".
Add to `generate`'s metadata

```
meta["b1_model"]                 = "miller" | "cdks" | "li6-convolution"
meta["b1_band_scale"]            = cfg.b1_band_scale
meta["b1_alpha_d_dwave_weight"]  = cfg.b1_alpha_d_dwave_weight
```
and a summary line in `cli.py` mirroring the FSI one (`cli.py:236-242`):
*"b1 li6-convolution: band scale 1.0, alpha-d D-wave 1.0 (band 0/1/2; never quote one
row alone)"*. Both go in Agent B's file list and in P1.

### 4.3 The mandatory 100 % band

`--b1-band-scale {0,1,2}` multiplies the **whole** b₁ (all four terms). Rationale, to
be reproduced verbatim in `docs/USAGE.md`: **Q(⁶Li) = −0.0806(6) fm² against
Q_d = +0.2859(3) fm²** — the α–d relative D-wave enters the closest measured
observable with the *opposite* sign to the deuteron's own D-state and nearly cancels
it. *(Primary sources, to be cited in the header and in `USAGE.md`, not via
`physics_literature.md`: Q(⁶Li) from the Pyykkö 2008 nuclear-quadrupole-moment
compilation, Mol. Phys. **106** (2008) 1965, whose ⁶Li entry is the molecular-beam
measurement of Cederberg et al., Phys. Rev. A **57** (1998) 2539; Q_d = 0.2859(3) fm²
from Bishop & Cheung, Phys. Rev. A **20** (1979) 381, as used by Ericson &
Rosa-Clot.)*

The two are different operators (charge quadrupole vs light-cone momentum alignment)
so the cancellation need not carry over, and this design **confirms the sign flip is
real in the VMC overlaps** (§2.2) — with the caveat that the α–d S–D sign is
**k-dependent**: opposite to the deuteron's below the S node (66 % of the density),
deuteron-like between the nodes (33 %). Terms (2d)/(2α) integrate across that node,
so the net sign is a computed output, not an assumption. Therefore: **never quote a
single row.** Every published number from this backend is `{0, 1, 2} × b₁`, or
equivalently "b₁ ± 100 %". The `w_alpha_d_dwave` knob is the *shape* variant of the
same worry and is reported separately (0, 1, 2 on terms (2d)+(2α) together).

### 4.4 pybind11

`python/bindings.cpp`, next to the `Li6B1` block (line 1086):
* `py::enum_<B1Model>` with the three values and the docstring above;
* `py::class_<Li6ConvolutionOptions>` with the plain fields `def_readwrite`, but the
  two `std::shared_ptr<const T>` fields (`deuteron_b1`, `unpol`) as **`def_property`
  with `std::const_pointer_cast`** — `def_readwrite` on a `shared_ptr<const T>` does
  not round-trip with the registered `shared_ptr<TensorSF>` / `shared_ptr<UnpolSF>`
  holders. This is exactly what `PipelineConfig::kernel` already does
  (`bindings.cpp:2331`);
* `py::class_<Li6ConvolutionB1, TensorSF, shared_ptr<…>>` with a **lambda `py::init`**
  taking `std::shared_ptr<TensorSF>` etc. (the pattern `Li6B1` uses at
  `bindings.cpp:1087`), exposing `b1`, `b1_embedded_s`, `b1_alpha_d_dwave`,
  `b1_alpha_d_dwave_d`, `b1_alpha_d_dwave_alpha`, `b1_cg_dwave`, `b1_alpha_d_sd`,
  `b1_alpha_d_dd`, the two P₂/P₄ remainders, `banded`, `close_kumano_integral`, and
  `densities()` / `densities_alpha()` / `options()` with
  **`py::return_value_policy::reference_internal`** (they return references into an
  object the `shared_ptr` owns);
* the `LightConeDensities` accessors (`f_s`, `f_d`, `delta_t_f`, `norm`, `mean_y`,
  `p_d`, `p_d_momentum`, `y_grid`);
* `py::class_<DeuteronConvolutionB1, TensorSF, …>` (the gate is scriptable) and
  `b1_landmarks` / `b1_landmarks_of_table` / `B1Landmarks`;
* `alpha_d_quadrupole_fm2`, `f1_cdks` as free functions;
* `PipelineConfig.b1_model`, `.b1_band_scale`, `.b1_alpha_d_dwave_weight` via
  `def_readwrite`;
* `m.attr("VMC_N_ALPHA_D_LI6") = VMC_N_ALPHA_D_LI6;` next to the existing
  `VMC_P_D_LI6` / `VMC_S_ALPHA_D_LI6` exports (`bindings.cpp:1868`).

### 4.5 Docs

* `docs/USAGE.md` — new §2a "b₁ backends for the inclusive channel" under §2
  Inclusive: the three models, the flags, the band rule, the "inclusive only, ⁶Li
  only" restriction, the primary sources for the two quadrupole moments, and the
  measured table of §7.
* `docs/OPEN_ITEMS_SOLUTIONS.md` — row 10 status `design only` → `implemented
  (opt-in)`, plus a §10 with the measured numbers, the G3b ratio, the ±5 % N_{αd}
  systematic and the resolution of Q1.
* `docs/CONVENTIONS.md` — one line under "Physics defaults that are a CHOICE":
  N_{αd} and where it lives (`VMC_N_ALPHA_D_LI6`, `tagged.hpp`), and one line that
  the κ = 1 δ-function is a default, not a derivation.
* `docs/PHYSICS_CHANNELS.md` (Phase A) row 3 gains the new backend.

## 5. Validation gate (A = 2) — **blocking**

**Nothing about ⁶Li may be quoted, plotted or merged until this passes.** The gate is
layered, cheapest first, so a failure localises itself.

### 5.1 Layer 0 — analytic identities (no data, no figure)

| id | statement | tolerance |
|---|---|---|
| G0a | ⟨P_zz⟩_{L=2} from the CG table = 1/10 | exact (rational arithmetic in the test) |
| G0b | ⟨S_z⟩_{L=2} = −1/2 | exact |
| G0c | Eq. (21) SD coefficient −3/(4√2π) reproduced by the Y-algebra of §1.2 | 1e-12 |
| G0d | Eq. (21) DD coefficient +3/(16π) reproduced likewise | 1e-12 |
| G0e | **∫ δ_T f(y) dy = 0** and **∫ (dy/y) δ_T f(y) = 0** — exact, because ∫dΩ(3c²−1) = 0 and ∫dΩ c(3c²−1) = 0, for *any* wave function and any k cutoff. Assert it for **both** density sets: struck d **and** struck α | ≤ 1e-4 × max\|δ_T f\| (quadrature-limited) |
| G0f | the P_zz-weighted S–D CG sum of §1.5, P_zz(0)C⁰₀ − (C⁺¹₀+C⁻¹₀)/2 = **+0.9486833**, and hence that the dropped b₁ᵈ-sector SD coefficient is exactly **1/3** of Eq. (21)'s −0.9486833 | 1e-12 |
| G0g | `ConvolutionKinematics::k_range` against a brute-force scan of \|c*\| ≤ 1 on a fine k grid, at κ = 1 **and** κ = 1.38, for both mass assignments (M_d/M_α and M_α/M_d) | endpoints to 1e-9 |

G0e is the single best quadrature diagnostic in the whole design: it is an exact zero
that the code must approach as the grids refine, and it caught a factor-6 endpoint bug
in the prototype for this document. G0g is the second: getting only one endpoint of
`k_range` right inflates ∫f dy by a factor 6.

### 5.2 Layer 1 — the unpolarized convolution (data, still no figure)

Fed the AV18 u(k), w(k) with the file's own ε_d = 2.224574 MeV, **before** any
renormalisation (`renormalize = false`):

| id | quantity | expected | tolerance |
|---|---|---|---|
| G1a | ∫k²[u²+w²]dk from `fdeut.av18` (spline) | 0.999976 | 1e-4 |
| G1b | P_D of the file | 0.0575999 (header says 0.057599) | 1e-5 |
| G1c | **∫ f(y) dy = 1 − ε_d/M_N − ⟨p²⟩/(2M_N²)**, both sides from the same table | **0.98707** | 1e-3 |
| G1d | `mean_y()` = ⟨y²⟩_φ/⟨y⟩_φ = [⟨E²⟩ + ⟨p²⟩/3]/(M_N⟨E⟩), E = 1 − ε/M_N − p²/(2M_N²) | **0.99546** | 1e-3 |
| G1e | y support upper end = `y_max()` = 1 − ε/M_N + κ²/2 | **1.4976** at κ = 1; **1.9502** at x = 0.8, Q² = 2.5 | 1e-3 |
| G1f | F₁ᴰ/F₁ᴺ per nucleon at x = 0.1 / 0.3 / 0.5 / 0.8, Q² = 2.5, ToyF2 | **0.997 / 0.990 / 0.986 / 1.108** | 1 % |

**G1c is the load-bearing one and it replaces a wrong pin.** With CDKS Eq. (12),
∫f dy = ⟨E⟩/M_N = 1 − ε/M_N − ⟨p²⟩/(2M_N²) **necessarily < 1** (the p_z term averages
to zero); the file's own `ti = 19.810` MeV gives 0.98707, and re-integrating ⟨p²⟩ from
the table gives 0.987065. The first revision pinned **1.0088**, which is reproducible
— but only as the **linear-interpolation artefact** of the coarse 0.1 fm⁻¹ table
(1.0088/0.98707 = 1.0220, and ∫k²(u²+w²)dk on a fine grid is 1.0218 linear against
1.0003 cubic). The baryon-number renormalisation then hides that 2.2 % quadrature
error, and T2's n_k/n_y doubling **cannot** detect it because the table spacing is the
problem. Spline the wave functions (§3) and pin the identity: **G1c is the only check
in this design sensitive to the interpolation bias.** Likewise G1d's 0.99546 is the
exact ⟨y²⟩_φ/⟨y⟩_φ, not ⟨y⟩ — the first revision's 0.9953 agreed only by cancellation.

G1f is the classic Fermi-motion/EMC shape (dip ~0.986 near x = 0.5, rise above
x ≈ 0.65). If it is not there, the kernel is wrong and the tensor numbers are
meaningless.

### 5.3 Layer 2 — the α–d sign gate

`alpha_d_quadrupole_fm2(phi0, phi2)` < 0, and — pinning the **node crossing**, not
just one point — φ₀φ₂ in the CDKS convention **> 0 at k = 0.2 fm⁻¹** (below the S
node) and **< 0 at k = 1.0 fm⁻¹** (between the S and D nodes), i.e. the SD bracket of
Eq. (21) is negative at low k, opposite to the deuteron, and changes sign at
0.678 fm⁻¹. The deuteron's own φ₀φ₂ must be < 0 at both, from the same code path.
See §2.2.

### 5.4 Layer 3 — the figure: `tables::kB1CdksQ2p5`

Digitized column `xb1_theory1_sum`, CDKS Fig. 4 theory 1 (SD+DD) at Q² = 2.5 GeV²,
x ∈ [0.0100, 1.590], 300 points. **Its landmarks are COMPUTED by
`b1_landmarks_of_table()` (§3), never typed** — the first revision typed them and got
four of five wrong. For reference, what that function returns:

| landmark | digitized value |
|---|---|
| sign change, **falling** | x = **0.06564** |
| minimum | x·b₁ = **−1.7681 × 10⁻⁴** at x = **0.3324** |
| sign change, **rising** | x = **0.45718** |
| maximum | x·b₁ = **+1.0852 × 10⁻³** at x = **0.7657** |
| ∫b₁ dx over the digitized range (Close–Kumano) | **+4.592 × 10⁻⁴** |
| last tabulated point | +4.04 × 10⁻⁶ at x = 1.590 (the curve does **not** die at 1.4) |

*(The first revision said "sign change (rising) x ≈ 0.062" — it is falling, at 0.0656;
"minimum −1.69e−4 at x ≈ 0.30" — that is just the table's value at x = 0.30, the
minimum is −1.768e−4 at 0.332; "sign change 0.42" — it is 0.457; "maximum 1.066e−3 at
0.80" — that is the value at x = 0.80, the maximum is 1.0852e−3 at 0.766.)*

**x range.** CDKS's convolution and Umnikov's warning
([hep-ph/9605291](https://arxiv.org/abs/hep-ph/9605291), PLB 391 (1997) 177) put the
validity floor at **x ≳ 0.1**; above x ≈ 0.8 CDKS themselves say the answer is
dominated by the high-momentum tail of the wave function. **Gate window: x ∈ [0.10,
0.80].** Report x < 0.1 and x > 0.8, gate nothing there.

**Acceptance criterion (three parts, all must hold).**

* **G3a — shape, hard.** The computed x·b₁ has exactly **two** sign changes in
  [0.02, 1.0], the first **falling** within **Δx = ±0.08** of 0.0656 and the second
  **rising** within **Δx = ±0.10** of 0.4572, and its maximum in [0.5, 1.0] sits
  within **Δx = ±0.10** of x = 0.766. Both reference values come from
  `b1_landmarks_of_table()`, so a re-digitization moves the target automatically.
* **G3b — magnitude, soft but recorded.** max |x·b₁| over [0.10, 0.80] agrees with
  **1.0852 × 10⁻³** — the raw-column maximum, compared against the **raw column**,
  because the per-nucleon question is *resolved* (§9 Q1) and the gate no longer has to
  absorb it — **within a factor of 2**, and the ratio is written into
  `docs/OPEN_ITEMS_SOLUTIONS.md` §10 to three digits. The factor 2 is budgeted, and
  each piece is reported separately so the residual is attributable: AV18 vs CD-Bonn
  is ~20 % on P_D and more on the tail; MSTW2008 LO vs ToyF2/CT18 is a factor ~1.5 on
  F₁ at x ≳ 0.6 (measured); the κ = 1 vs κ choice is a factor 1.5–1.7 (§1.3); and the
  curve is a *figure digitization*. **A factor 2 is not permission to be a factor 2
  out — it is the width at which the gate stops being able to distinguish those four.**
* **G3c — Close–Kumano, reported.** ∫b₁ dx from the new code over the same range,
  next to `close_kumano_integral(true)` = **+4.592 × 10⁻⁴** and
  `close_kumano_integral(false)` (Miller). Neither is zero; neither is enforced
  (`sf.hpp`: "Reported, not enforced"). Record both.

**Honest warning to the implementer — read this before starting.** A prototype written
for this design (AV18 u,w; CDKS Eqs. 16/17/21 exactly as transcribed at κ = 1; the
library ToyF2 with R = `r_sigma_lt`; target-mass factor on/off) reproduces the
**shape** — positive at very low x, a negative dip in the middle, a positive bump at
large x — but lands a factor **3.7 (at the large-x peak) to ~10 (at the mid-x dip)
BELOW** the digitized curve, with sign changes at x ≈ 0.02 and ≈ 0.36 (the independent
re-run for revision 2; the first revision reported 0.03 and 0.28) rather than 0.066
and 0.457:

| x | prototype x·b₁, massless F₁ | prototype x·b₁, CDKS Eq. (22) F₁, κ = 1 | the same with κ = √(1+γ²) | digitized |
|---|---|---|---|---|
| 0.10 | −2.17e−5 | −1.64e−5 | −1.67e−5 | −1.73e−5 |
| 0.20 | −2.07e−5 | −2.03e−5 | −2.13e−5 | −8.78e−5 |
| 0.30 | +1.19e−5 | −1.60e−5 | −1.51e−5 | −1.69e−4 |
| 0.50 | +1.56e−4 | +9.36e−5 | +1.45e−4 | +1.48e−4 |
| 0.70 | +2.33e−4 | +2.84e−4 | +4.52e−4 | +9.97e−4 |
| 0.80 | +1.89e−4 | +2.88e−4 | +4.80e−4 | +1.07e−3 |
| ∫b₁dx | +1.28e−4 | +1.24e−4 | +2.41e−4 | +4.59e−4 |

**Therefore G3 as stated will probably FAIL on the first pass, and that is the point of
having it.** Work the following checklist, in order, recording each answer in the
commit message; do **not** tune a fudge factor:

0. **The δ-function.** Eq. (17)'s κ = 1 or Eq. (21)'s κ = √(1+γ²)? (§1.3.) Worth
   **1.54 / 1.59 / 1.67 at x = 0.5 / 0.7 / 0.8** and a factor 1.94 on ∫b₁ dx — the
   single largest identified effect, and the reason `finite_q_delta` exists. Run the
   gate both ways and report both columns; the default stays κ = 1 unless the gate
   says otherwise, and if it does, **that is a finding, not a fudge**.
1. **Target mass.** Is CDKS Eq. (22) in? (Worth ×1.5 at x = 0.8; the table above.)
2. **R.** `r1998` (CDKS's SLAC-R1998), not `r_sigma_lt`. Worth a few %.
3. **Per nucleon vs per deuteron — RESOLVED, do not re-litigate.** CDKS Eq. (10)
   carries 1/A, Eq. (16)'s text says *"b₁ is defined by the one per nucleon"*, f is
   normalised to one nucleon, and their Fig. 6 overlays HERMES's per-nucleon b₁ on the
   same axis. The digitized column **is per nucleon**, so the gate compares against the
   **raw table column** and `Li6ConvolutionB1`'s default b₁ᵈ is the raw column
   (§2.4, §9 Q1). `b1_convolution()`'s `B1_PER_DEUTERON_TO_PER_NUCLEON = 0.5` makes
   `CdksB1` a factor 2 low; that is a separate open item against `constants.hpp:69-71`
   (§9 Q1) and **must not** be absorbed into this gate's tolerance.
4. **Nucleon PDF.** ToyF2 vs MSTW2008 LO. Rerun the gate with
   `LhapdfSF("CT18NLO")` (installed: `../deps/install/share/LHAPDF/CT18NLO`) — from
   the validation script, not from the core, which must not link LHAPDF. The
   large-x valence exponent moves x·b₁(0.8) by a factor 2.6 between (1−x)³ and (1−x)⁵
   (measured).
5. **Wave function.** AV18 (5.76 %) vs CD-Bonn (4.85 %). If everything else is
   settled and a residual factor remains, quantify it by scaling W(k) to P_D = 4.85 %
   and re-running; report the residual.
6. **Interpolation.** Confirm G1c on the spline (0.98707) *before* believing any of
   the above; a linear-interpolated u, w is 2.2 % high on the norm and the
   renormalisation hides it (§3, §5.2).

**Escalation.** If after the checklist the peak ratio is still outside a factor of 2,
**stop and hand back**: publish the layer-0/1/2 results, the checklist answers and
the prototype table, leave `Li6ConvolutionB1` in the tree behind the flag with a
`WARNING: A=2 gate not passed` in its header and in `--help`, and put **no ⁶Li number
in `OPEN_ITEMS_SOLUTIONS.md`**. Item 10 stays open. That is a legitimate and expected
outcome for a first-mover calculation and is much better than a tuned number.

## 6. Tests

C++ (doctest) in **`tests/test_b1_nuclear.cpp`** — picked up automatically by the
`file(GLOB tests/test_*.cpp)` at `CMakeLists.txt:144`. Python in
**`python/tests/test_b1_model.py`**.

**Two deliberate departures from the repo's test conventions, stated here so a
reviewer does not "fix" them.** (i) Every existing data-dependent test *skips* when
`data/vmc` is absent (`if (!have(kLi6Momentum)) { MESSAGE(...); return; }`,
`test_cluster.cpp:180`). **T1 must not skip**: §8's review step (iii) requires the
suite to go RED when the gate's data file is removed, so T1 uses
`REQUIRE(have(...))` and lets `read_fdeut_k` throw. A blocking gate that silently
skips is not a blocking gate. (ii) The "no RNG" grep of T8 is **not** expressible in
the doctest binary (only `LIPOLGEN_REFERENCE_DIR` is compiled in, not the source
directory); it moves to `python/tests/test_b1_model.py`, which knows `repo_root`.

| id | test | assertion |
|---|---|---|
| T1 | **A = 2 gate** (§5, layers 0–3) | G0a–g, G1a–f exact/tight; G3a hard against `b1_landmarks_of_table()`; G3b/G3c `MESSAGE`d, and `CHECK`ed against whatever the implementer records, so a later regression is caught. **`REQUIRE(have(...))`, never a skip** |
| T2 | **quadrature convergence** | `n_k` and the three `n_y_*` ×2 and ×4: b₁(x) at x = 0.05…0.5 stable to **1e-3 relative**; ∫δ_T f dy → 0 faster than 1/n; **and** the raw `norm()` (renormalize = false) matches G1c's moment identity to 1e-3 — the only clause of T2 that can see an interpolation bias, because refining the y/k grids cannot fix the table spacing (§3) |
| T3 | **z-density normalisation** | raw `norm()` (renormalize = false) == `VMC_N_ALPHA_D_LI6 × mean_E/M_d` = **0.8174** to 1e-3; the *trapezoid* S+D sum of the same table == `VMC_N_ALPHA_D_LI6` to **3e-4** (the file's own quadrature spread, §2.1 — not 1e-9, which would be tautological once `renormalize` sets it); `p_d_momentum()` == `VMC_P_D_LI6` at rtol **1e-4**; `p_d()` (the y-weighted one, = P_D⟨y⟩_D/⟨y⟩ = **0.019337**) compared with `CHECK_CLOSE_AT(..., VMC_P_D_LI6, 0.0, 1e-4)` — **absolute**, because `CHECK_CLOSE` is relative and the two differ by 9e-4 relative by construction |
| T4 | **⟨z⟩** | `mean_y()` = **0.99979 ± 0.002** (a pin, and it is ⟨z²⟩/⟨z⟩, **not** ⟨z⟩); and it equals the identity it actually satisfies, [⟨E²⟩ + ⟨k²⟩/(3M_d²)]/(M_d⟨E⟩) with E = M_d − ε_{αd} − k²/(2M_α), to 1e-6. The closed form 1 − ε/M_d − ⟨k²⟩/(2M_αM_d) = **0.997486** is a *different* number: it differs by Var/⟨z⟩ ≈ ⟨k²⟩/(3M_d²) = 2.3e−3 and the first revision's "equals … to 1e-3" would fail a correct implementation |
| T5 | **D-wave terms vanish as P_D(αd) → 0** | with φ₂ ≡ 0: `b1_alpha_d_dwave_d`, `b1_alpha_d_dwave_alpha` and `b1_cg_dwave` == 0 exactly, and `b1` == `b1_embedded_s` bit-for-bit |
| T6 | **CG term reduces to the known deuteron limit** | with `LightConeDensities::delta_limit` (§3) `b1` == `b1_li6_from_deuteron(b1_d, 1 − 0.9*P_D, 2/6)` to 1e-12; and at P_D = `P_D_LI6` the weight equals `LI6_B1_RANK2_TRANSFER` to **1e-4** (the same 1e-4 pinning the tagged quadrature already uses, `constants.hpp:74-78`). **This exercises the weight assembly only** — the δ short-circuits `convolve`; the kernel is covered by T2/T5 |
| T7 | **band and term knobs** | `banded(0)` ≡ 0; `banded(2)` == 2×`banded(1)` bit-for-bit; `w_embedded_s = 0` and `w_cg_dwave = 0` kill exactly their own term; `w_alpha_d_dwave = 0` kills **both** (2d) and (2α) and nothing else; `b1` == `b1_embedded_s + b1_alpha_d_dwave_d + b1_alpha_d_dwave_alpha + b1_cg_dwave` to 1e-12, and `b1_alpha_d_dwave` == (2d)+(2α) exactly |
| T8 | **determinism** | two constructions from the same options give bit-identical b₁ on a 200-point x grid (the grep clause moves to P4) |
| T9 | **default unchanged** | `default_inclusive_kernel(LI6())` with `B1Model::Miller` gives bit-identical b₁ to `validation/reference/b1_default_li6.json` at rtol **1e-12**. **That file does not exist yet and is Agent B's step 0** (§8): the existing `validation/reference/*.json` are dumped from the Python with `toy_b1(mode=kToy)` and `test_reference.cpp::build_kernel` never calls `default_inclusive_kernel`, so **the existing rtol-1e-12 gates do not cover this path at all** and cannot catch a regression in it. Without step 0 this test is tautological |
| T10 | **α–d sign gate** | `alpha_d_quadrupole_fm2` < 0; φ₀φ₂ (CDKS convention) **> 0 at k = 0.2 fm⁻¹ and < 0 at k = 1.0 fm⁻¹**, pinning the S node at 0.678 fm⁻¹; the *deuteron's* φ₀φ₂ < 0 at both, from the same code path |
| T11 | **channel legality** | `PipelineConfig::validate()` throws for `isotope == "7Li"`, for `isotope == "d"` (the spin-1 trap of §4.2), for a tagged/coherent channel, and for a non-null `cfg.kernel` with `b1_model != Miller`; and passes for `("6Li", Inclusive)`. **No `struck_cluster_kernel` clause** — it takes no `B1Model` (§4.1) |
| T12 | **VMC spread, documented not explained** | `VMC_N_ALPHA_D_LI6` = 0.819481 against the published 0.846 + 0.017 = 0.863 and against `li6.ad`'s own printed 0.856 — agreement to **6 %**; P_D 0.019355 vs 0.0197 to **2 %**. The test's `MESSAGE` says *"three tabulations, 5 % spread, unexplained; carried as a systematic"*, not "a version difference" |
| T13 | **CDKS Eq. (22) is not `f1_from_f2`** | `f1_cdks(f2,x,q2)` / `f1_from_f2(f2,x,q2)` == 1 + `gamma_squared(x,q2)` to 1e-14, and == 1.90 at x = 0.8, Q² = 2.5 to 1e-2; **and** `Li6ConvolutionB1`'s internal F₁ᵈ over `NuclearF2(DEUTERON(), unpol).f1a/2` == 1 + γ² at x = 0.5, Q² = 2.5 (1.35) to 1e-12 — the check that `f1a` did not sneak back in |
| T14 | **the angular-average truncation, measured** | `b1_cg_dwave_p2_remainder` and `_p4_remainder` over `b1_embedded_s` at the §7 x points, `MESSAGE`d and pinned to whatever they come out at (rtol 1e-6) so the number in the header comment and in `OPEN_ITEMS_SOLUTIONS.md` §10 is a *measurement*, replacing the first revision's asserted "0.2 %" |
| T15 | **struck-α term is real** | `b1_alpha_d_dwave_alpha` / `b1_alpha_d_dwave_d` ∈ [0.45, 0.60] at x = 0.05…0.7 and has the **same sign**, against the analytic 2(M_d/M_α)² = 0.506 (§1.5). A regression that drops term (2α) then fails loudly rather than quietly shifting b₁ by 30 % |
| P1 | pytest: enum + config round-trip | `lg.B1Model.Li6Convolution`, `make_config(b1_model="li6-convolution")`, `cfg.b1_band_scale`, `cli.resolve(["--b1-model","li6-convolution"])`, all three `DEFAULTS` keys, and `meta["b1_model"] / ["b1_band_scale"] / ["b1_alpha_d_dwave_weight"]` present in a `generate` dict (§4.2) |
| P2 | pytest: the four terms sum | `b1 == b1_embedded_s + b1_alpha_d_dwave_d + b1_alpha_d_dwave_alpha + b1_cg_dwave` at 5 x values |
| P3 | pytest: band | with `--b1-band-scale 0` vs `1` on the same seed: **`cell_xsec_pb()` (spin-blind) is bit-identical**, and so is the P_zz = 0 category's `sigma_per_category_pb`. At band 0 the P_zz ≠ 0 categories' `sigma_per_category_pb` are all **equal to each other** (the tensor shift is zero); at band 1 they are **not**. *(The first revision said "identical unpolarized rate"; with the tensor-thirds plan `sigma_pb()` includes Σ_m p_m(1 + w_avg(m)), which depends on b₁ through the P_zz-weighted categories, so the totals are NOT identical and that phrasing would fail.)* |
| P4 | pytest: no RNG in the kernel | `repo_root/src/core/b1_nuclear.cpp` contains no `rng`, `random`, `rand` (the grep T8 cannot do) |

---

## 7. Number placeholders the implementer must fill

Same table into `docs/OPEN_ITEMS_SOLUTIONS.md` §10 and `docs/USAGE.md` §2a.
**x·b₁ per nucleon at Q² = 2.5 GeV², ⁶Li.**

| x | (1) embedded d, S | (3) CG D-wave | (2d) α–d D-wave, struck d | (2α) α–d D-wave, struck α | total | today's default `Li6B1(MillerB1)` | `Li6B1(CdksB1)` |
|---|---|---|---|---|---|---|---|
| 0.05 | ☐ | ☐ | ☐ | ☐ | ☐ | ☐ | ☐ |
| 0.10 | ☐ | ☐ | ☐ | ☐ | ☐ | ☐ | ☐ |
| 0.20 | ☐ | ☐ | ☐ | ☐ | ☐ | ☐ | ☐ |
| 0.30 | ☐ | ☐ | ☐ | ☐ | ☐ | ☐ | ☐ |
| 0.50 | ☐ | ☐ | ☐ | ☐ | ☐ | ☐ | ☐ |

plus: ∫b₁ dx (Close–Kumano) for the new model and for both `Li6B1` variants; the band
rows (`--b1-band-scale 0/1/2`); the term-(2) knob rows (0/1/2); the `finite_q_delta`
on/off pair (§1.3); `norm()` (raw and renormalised), `p_d()`, `p_d_momentum()`,
`mean_y()`; T14's P₂/P₄ remainders as a fraction of term (1); the ±5 % N_{αd}
systematic; and the A = 2 gate's G3b ratio.

**Provisional values from the design prototype** (α–d VMC densities + the digitized
CDKS b₁ᵈ **per nucleon, raw column** + ToyF2 F₁ᵈ, κ = 1; **illustrative only — they
inherit the unresolved A = 2 factor of §5.4 and are not a substitute for the measured
table**). *The (2α) column is a revision-2 addition; the (2d) column is the first
revision's "term (2)". Note that the first revision's numbers were computed with the
half-sized `CdksB1` default for b₁ᵈ, so terms (1) and (3) below are a factor 2 low
relative to what §2.4 now specifies — one more reason to treat this table as a shape,
not a value.*

| x | (1) | (3) | (2d) | (2α) ≈ 0.5 × (2d) | total | `Li6B1(MillerB1)` | `Li6B1(CdksB1)` |
|---|---|---|---|---|---|---|---|
| 0.05 | +7.85e−7 | +1.51e−9 | +1.28e−6 | +6.4e−7 | +2.71e−6 | +3.65e−4 | +8.99e−7 |
| 0.10 | −2.35e−6 | −4.54e−9 | +1.10e−6 | +5.5e−7 | −7.0e−7 | +4.12e−4 | −2.66e−6 |
| 0.20 | −1.18e−5 | −2.28e−8 | +1.76e−6 | +8.8e−7 | −9.2e−6 | +2.67e−4 | −1.35e−5 |
| 0.30 | −2.25e−5 | −4.31e−8 | +3.37e−6 | +1.70e−6 | −1.75e−5 | +6.72e−6 | −2.60e−5 |
| 0.50 | +2.12e−5 | +4.18e−8 | +9.35e−6 | +4.73e−6 | +3.53e−5 | −1.87e−4 | +2.28e−5 |

Reading of the prototype table, for the discussion section:
* **Term (3) is 0.2 % of term (1)** — as the CG algebra says (0.1 × P_D = 0.0019).
  It is kept because it is free and because it is the analytic bridge to
  `LI6_B1_RANK2_TRANSFER`, not because it matters numerically. *(This 0.2 % is the
  measured size of term (3); it is not the truncation bound, which T14 measures
  separately — the first revision conflated the two.)*
* **Terms (2d) + (2α) together are 25–200 % of term (1)** and have the **opposite
  sign** over most of the window, because the α–d S–D relative sign is opposite to the
  deuteron's below the S node, which carries 66 % of the density (§2.2). They cancel
  term (1) near x ≈ 0.1 and dominate below. **This is the quadrupole puzzle showing up
  in b₁, and it is the quantitative reason the 100 % band is mandatory** — and the
  struck-α piece is *half* of it, which is why revision 2 exists.
* Terms (1)+(3) come to **0.8056 × b₁ᵈ × (2/6)** against the default's
  0.921947 × (2/6): the ratio **0.874** decomposes as
  N_{αd}(1 − 0.9 P_D^VMC)/(1 − 0.9 P_D^Hulthén) = 0.8195 × 0.98258/0.92195 = 0.874.
  *(The first revision wrote "entirely N_{αd}/[1 − 0.9 P_D^{αd}] = 0.8195/0.9830";
  that quotient is 0.834, not 0.874.)*
* Against the **production** default `Li6B1(MillerB1)` the new model differs by two
  orders of magnitude and by sign, because Miller (HERMES-like) and CDKS
  (convolution) are different *camps* for b₁ᵈ, a pre-existing disagreement this design
  does not resolve. **Do not change the default.**

---

## 8. Implementation steps — two Opus agents

Agent A and Agent B work on **disjoint files within design D**; agent A must land
first because agent B's tests import its header.

**Coordination with design C (tensor RC), which runs in the same worktree in this
run.** Design C edits `include/lipolgen/pipeline.hpp`, `src/core/pipeline.cpp`
(`PipelineConfig` + `validate()`), `python/bindings.cpp` (`PipelineConfig` block,
`Columns`, `columns_to_dict`), `python/lipolgen/__init__.py` and
`python/lipolgen/cli.py` (`DEFAULTS`, `make_config`, the summary banner) — **the same
five files as design D's Agent B**. "Disjoint files" therefore holds only *inside*
design D. The rule for this run:

* **Design C's Agent B lands first; design D's Agent B rebases onto it.** C's change
  is larger in those five files and its `DEFAULTS` keys are already specified.
* Both add keys to `cli.py::DEFAULTS`, and `cli.resolve()` rejects any config-file key
  not in `DEFAULTS` (`cli.py:168`), so **every key from both designs is mandatory** —
  a partial merge makes existing config files fail with "unknown option(s)".
* Both add a clause to `PipelineConfig::validate()` and a block to the
  `py::class_<PipelineConfig>` chain; those are append-only and merge cleanly if C is
  in first.
* If the two are handed to one wiring agent instead, that agent does C's block then
  D's, in that order, and runs both pytest files before either commit.

**Names.** `PLAN.md` Phase D calls the flag `--b1-model {toy,li6-convolution}` and the
class `B1Li6Convolution`. This design uses `--b1-model {miller,cdks,li6-convolution}`
and `Li6ConvolutionB1`. **`miller` IS the PLAN's `toy`**: today's default b1_func is
`Li6B1(MillerB1)`, which is what `toy_b1` reaches, and `cdks` is the second camp that
already exists in `sf.hpp` and had no flag. Update `PLAN.md` and `STATUS.md` to the
chosen names as part of Agent B's docs step, rather than leaving two vocabularies.

### Agent A — kernel and convolution (physics-bearing, Opus)

**Files created**
* `include/lipolgen/b1_nuclear.hpp` (§3)
* `src/core/b1_nuclear.cpp`

**Files edited**
* `include/lipolgen/cluster.hpp`, `src/core/cluster.cpp` — add `read_fdeut_k` only
  (returning `ebind` as well as k, u, w; converting with `HBARC_GEV_FM`, §2.3).
* `include/lipolgen/tagged.hpp` — add `VMC_N_ALPHA_D_LI6 = 0.80362 + 0.015861` and
  rewrite `VMC_P_D_LI6` as `0.015861 / VMC_N_ALPHA_D_LI6` (§2.1). This is a **value-
  preserving** edit (`VMC_P_D_LI6` is bit-identical) and it is the only way N_{αd} has
  one home; `python/tests/test_module.py:252` pins `VMC_P_D_LI6` and must still pass.

**Steps**
1. `read_fdeut_k`: locate the `k  u(k)  w(k)` header (line 10221), read 3-column
   rows with the existing `split_numbers`, convert k to GeV with `HBARC_GEV_FM`, and
   return the header's `ebind` alongside. Assert ∫k²(u²+w²)dk = 1 to 1e-4 in a test,
   not in the reader.
2. `ClusterPartialWave::from_uw` (deuteron: φ₀ = u, φ₂ = −w) and `::from_vmc` /
   `::from_vmc_pair` (α–d: take the `VmcRadial`s from `li6_vmc_waves()`, apply
   **only** φ_L = iᴸψ_L, and **do not** re-fix the global phase — §2.2, §3). Zero
   outside the table, cubic-spline inside. **Write the iᴸ paragraph of §1.2 and the
   node table of §2.2 into the comment.**
3. `ConvolutionKinematics::k_range` with **both** endpoints and with κ (§1.3).
   Unit-test it against a brute-force scan (G0g) before anything else.
4. `LightConeDensities`: the densities (including `f_d_p2`, `f_d_p4`), the
   `Options`-specified y grid, the baryon-number renormalisation,
   `norm/mean_y/p_d/p_d_momentum`, `convolve`, and `delta_limit`.
5. `f1_cdks` (reusing `gamma_squared`), `DeuteronConvolutionB1`, `b1_landmarks` /
   `b1_landmarks_of_table`, and the **A = 2 gate driver** as a standalone function
   returning the landmarks of §5.4 so both the doctest and a validation script can
   call it.
6. `alpha_d_quadrupole_fm2`.
7. `Li6ConvolutionOptions` / `Li6ConvolutionB1` with **two** `LightConeDensities`
   (struck d and struck α), the four term accessors plus the two P₂/P₄ remainders,
   `banded`, `close_kumano_integral`.
8. Run the gate. Record the checklist of §5.4 — **including item 0, the κ column** —
   in the commit message.

**Do not touch**: `sf.hpp`, `constants.hpp` (`LI6_B1_PER_NUCLEON`,
`LI6_B1_RANK2_TRANSFER`, `M_NUCLEON`, `HBARC_GEV_FM`,
`B1_PER_DEUTERON_TO_PER_NUCLEON` — the last is §9 Q1's separate item and **must not**
change in this phase), `spectator.*`, `pipeline.*`, bindings.

### Agent B — wiring, bindings, tests, docs (Opus)

**Step 0, before touching `pipeline.*` — the reference T9 compares against.**
Nothing in the tree stores today's `default_inclusive_kernel(LI6())` values:
`validation/reference/*.json` are dumped from the Python with `toy_b1(mode=kToy)` and
`test_reference.cpp::build_kernel` builds its own `InclusiveKernel::Options`, so the
existing rtol-1e-12 gates never exercise that function. Dump
`default_inclusive_kernel(LI6())`'s b₁ (and Δ) on 200 x points × 3 Q² through
`InclusiveKernel::tables` at repr precision into
`validation/reference/b1_default_li6.json`, add it to `_manifest.json`'s `files` list,
and only then edit `pipeline.cpp`. T9 reads that file at rtol 1e-12.

**Files edited**
* `include/lipolgen/pipeline.hpp`, `src/core/pipeline.cpp` — `B1Model`,
  `PipelineConfig::{b1_model, b1_band_scale, b1_alpha_d_dwave_weight}`,
  `default_inclusive_kernel` overload **with the shared_ptr-capturing closure of
  §4.2** (no `static` on the non-Miller branches; the Miller branch's `static` stays,
  untouched), `validate()`'s isotope + channel + kernel-conflict rules, and the
  `struck_cluster_kernel` **comment** of §4.1 (no code change there).
* `python/bindings.cpp` (§4.4, including the three `meta` keys of §4.2)
* `python/lipolgen/__init__.py` (`B1_MODELS`, `make_config`)
* `python/lipolgen/cli.py` (three flags, three `DEFAULTS` keys, the summary line)
* `docs/USAGE.md` (new §2a), `docs/OPEN_ITEMS_SOLUTIONS.md` (row 10 + §10),
  `docs/CONVENTIONS.md` (two lines), `PLAN.md` / `STATUS.md` (the naming, above)

**Files created**
* `validation/reference/b1_default_li6.json` (step 0)
* `tests/test_b1_nuclear.cpp` (T1–T15)
* `python/tests/test_b1_model.py` (P1–P4)
* `validation/b1_li6_table.py` — regenerates the §7 table and the A = 2 gate report,
  runs the LHAPDF variant of checklist step 4; writes nothing outside `docs/` when run
  with `--write`.

**Order**: step 0, then T9 (default unchanged) and T11 (channel legality) — they are
the two that protect everything already published.

### Review (Fable, after both)

Adversarial verify with ≥3 lenses: (i) every equation number in the header comments
resolves to the transcription in this document and to the PDF; (ii) no physics number
appears twice — `grep` for `0.8195`, `0.819481`, `0.80362`, `0.015861`, `0.01935`,
`0.9219`, `0.3333`, `1.4743`, `2.2245`, `0.19733`, `0.1973269804` across
`include/ src/ python/`, and confirm `VMC_P_D_LI6` is now expressed through
`VMC_N_ALPHA_D_LI6`; (iii) the A = 2 gate is genuinely blocking — delete
`data/vmc/deuteron/fdeut.av18`, confirm the suite goes **red, not skipped** (§6);
(iv) the four term accessors sum to `b1` and term (2α) is present and ≈ 0.5 × (2d)
(a regression that silently drops it is this design's specific failure mode).
Completeness critic: is every new CLI flag in `--help`, `docs/USAGE.md` **and** a
pytest, and does every generated file record `b1_model` and `b1_band_scale`?
---

## 9. Open questions

**Q1 (RESOLVED by the paper; the follow-up is a separate item). Is the digitized CDKS
curve per nucleon?** **Yes.** CDKS Eq. (10) carries the explicit 1/A; the text under
Eq. (16) says *"the structure function b₁ is defined by the one per nucleon"*; their
f(y) is normalised to **one** nucleon (∫f dy = 1); and their Fig. 6 overlays HERMES's
per-nucleon b₁ᵈ on the same axis as the curve digitized into
`tables::kB1CdksQ2p5`. The first revision left this open and widened G3b's tolerance
to "a factor 2.5 … the per-nucleon vs per-deuteron convention is a factor 2" — **a
gate whose tolerance is designed to contain a known factor-2 ambiguity is not a
physics gate**, so that is gone: G3b compares against the raw column at a factor 2
(§5.4), and `Li6ConvolutionB1`'s default b₁ᵈ is the raw column (§2.4).

**Q1b (NEW, and it is not this phase's to fix). Is
`B1_PER_DEUTERON_TO_PER_NUCLEON = 0.5` right where it is applied?**
`b1_convolution()` multiplies the raw column by 0.5 (`sf.cpp:356-361`), so `CdksB1` —
which the first revision made this model's default deuteron b₁ᵈ — is **a factor 2
low**, and terms (1) and (3) inherited that. This design routes around it (§2.4) but
does **not** change the constant, because every published number carries the current
convention. Two things must happen in the close-out phase, separately: (a) decide
whether `CdksB1` should stop halving, which moves published CDKS-camp numbers; (b)
check **Miller's** normalisation independently — it may differ from CDKS's, and the
answer is on the axis of Miller's own Fig. 5, not assumable from CDKS. Open a distinct
item for both.

**Q2 (RESOLVED). Which of the Fig. 4 curves is `xb1_theory1_sum`?** The sum. In
Fig. 4 both SD and DD are positive at the large-x peak, so the SD+DD curve must be the
**topmost** of the three; the first revision's reading ("solid ≈ 9e−4, dashed ≈
1.15e−3, and the digitized 1.066e−3 sits between them") is impossible — the sum cannot
be below a component of the same sign. The digitized maximum, 1.0852e−3 at x = 0.766,
**is** the sum. The first revision's checklist item 6 ("re-digitize if the residual is
~1.2") is therefore deleted; there is no digitization ambiguity of that size to chase.
*(Asking S. Kumano's group for the tabulated curves is still worth doing, but it is
not a gate dependency.)*

**Q3. The residual A = 2 factor.** The prototype is 3.7–10× below the digitized curve
at κ = 1 with CDKS Eq. (22) F₁. Checklist item 0 (the finite-|q⃗| δ-function) recovers
1.5–1.7 at large x and a factor 1.94 on ∫b₁ dx on its own — the largest single
identified piece — and items 1–6 are the rest of the plan. If it does not close, item
10 stays open and no ⁶Li number ships. **This is the single largest risk in the
phase.**

**Q4. Normalisation to N_{αd}: suppression or renormalisation?** Default here: the
18 % of ⁶Li that is not α+d gets b₁ = 0 (assumption A1), so b₁ carries a factor 0.82.
The alternative — renormalise the density to 1, i.e. give the non-cluster
configurations the *same* b₁ per unit weight — is 22 % higher. Both are defensible;
the knob exists (`use_spectroscopic_factor`), the default is the conservative one, and
the difference is well inside the 100 % band. Nobody has published either. **On top of
that sits the ±5 % spread between the three tabulations of N_{αd} itself (§2.1),
carried as a stated systematic through `norm_target`.**

**Q5. Q² evolution (A8).** b₁ᵈ, F₁ᵈ and F₁^α are evaluated at the ⁶Li kernel's Q²,
with no rescaling by z. Standard, but at Q² = 2.5 and z-smearing of ±10 % it is a
percent-level effect that has never been checked for b₁.

**Q6. Should `LI6_B1_RANK2_TRANSFER` be re-evaluated at the VMC P_D?** The constant is
0.921947, tied to the **Hulthén scenario** P_D = 0.0867; the VMC value 0.0193549 gives
0.98258. Changing it would move every published inclusive tensor number, so it is out
of scope here — but the new backend makes the inconsistency visible (the tagged channel
already uses VMC when `--cluster-wave vmc`, the inclusive b₁ never does). Flag for the
close-out phase.

**Q7. The δ-function and the upper limit — now a live choice, not a curiosity.**
Eq. (17) (κ = 1) and Eq. (21) (κ = √(1+γ²)) are different equations in the same paper
and give answers differing by 1.5–1.7 at large x (§1.3); this design defaults to
Eq. (17) and exposes Eq. (21) as `finite_q_delta`, and **the default is a choice, not
a derivation**. The "y_max = 2" on CDKS p. 10 belongs to **Ref. [4],
Khan–Hoodbhoy**, not to CDKS; CDKS's own support is 1 − ε/M_N + κ²/2 (= 1.4976 at
κ = 1, 1.9502 at x = 0.8, Q² = 2.5), and the digitized curve's non-zero tail out to
x = 1.59 says nothing about it either way. If a later implementer wants the light-front
variant (their Eqs. 40–44, transcribed in §1.2), the α–d version needs Melosh rotations
for a spin-1 constituent, which is *not* in the paper. Out of scope; noted so it is not
rediscovered.

**Q8. Nothing to validate against for A = 6.** There is no published b₁ for any
A > 2 (`physics_literature.md` §2, confirmed by title/fulltext/citation search). The
A = 2 gate plus the analytic limits of §6 are the *only* checks that exist. Every ⁶Li
number from this backend is a prediction, and the 100 % band is not a formality.

**Q9. Second-order α–d effects not modelled.** ⁶Li also has α + p + n and
³He + ³H configurations, and the deuteron inside ⁶Li is itself distorted (its own P_D
need not be the free 5.76 %). Both fold into the (1 − N_{αd}) hole of Q4.

**Q10 (NEW). The α's own EMC effect (A9).** F₁^{α, per nucleon} is taken as the
isoscalar nucleon F₁, i.e. identical to F₁^{d, per nucleon}. The real α is EMC-
suppressed by ≈ 10 % near x = 0.6 and enhanced below x ≈ 0.3. That is a ≤ 10 % error
on term (2α), which is itself ≈ 0.5 × term (2d) — so ≤ 5 % of the orbital sector,
inside the band, but it is the first thing to add if term (2α) turns out to drive the
answer. `Li6ConvolutionOptions::alpha_f1` is the hook.

**Q11 (NEW). Is the struck-α term derived correctly, or only plausibly?** §1.5's
argument is that CDKS Eq. (10) sums over constituents and that the α's density carries
the same (3cos²θ − 1), which is solid. What is *not* checked against any paper is the
per-nucleon bookkeeping of the α piece (the 4/6 and the fair-share z_α), because no
published two-cluster b₁ convolution exists to check it against (Q8). The independent
handle is the analytic scaling 2(M_d/M_α)² = 0.506 against the quadrature's 0.502–0.511
(T15) — a consistency check, not a validation. **State it that way in the header.**
