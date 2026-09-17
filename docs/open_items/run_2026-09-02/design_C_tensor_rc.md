# Design C — the tensor-sector radiative-correction band (`rc.hpp`)

Open item 9 (`docs/OPEN_ITEMS_SOLUTIONS.md` row 9 and §"8–10 Theory notes";
recommendation verbatim in `docs/open_items/physics_literature.md` §3).
Design only — no code was written for this document.

**One sentence.** Add an opt-in, weight-only radiative-correction family
`--rc tensor-band` that (i) prices the *lepton-vertex* RC uncertainty on the
**tensor part alone** of every event as a two-sided band `rc_tensor_lo` /
`rc_tensor_hi`, and (ii) carries the **⁶Li radiative tails** — the elastic
tail with its tensor part, plus the *unpolarised* quasi-elastic tail — as an
additive per-event weight `rc_tail`, because a tail whose own `A_zz` differs
from the Born's does not cancel in `A_zz`: **it dilutes it.** The tail is
overwhelmingly an unpolarised dilution; its tensor piece is a small correction
on top (≈ 10 % of the deuteron's, §2.1). The nominal weight, every four-vector
and the RNG stream are untouched, so `--rc off` (the default) is bit-for-bit
today.

> **Revision note (review of 2026-09-02).** This document was revised against a
> physics + integration review. The substantive changes: the tagged τ algebra
> of §1.5.1 (it is the *same* `(P_zz/2)A_zz/[1+(P_zz/2)A_zz]` as the inclusive
> one, not `½ A_zz`); the elastic tail is now built on POLRAD's **t-peak closed
> forms Eqs. (37)–(39)** rather than the full Eq. (18)/Appendix-B machinery
> (POLRAD §2.1.3 B says the s- and p-peaks are *suppressed* for a tail, so the
> old leading-log normalisation gate was measuring the wrong thing); the
> unpolarised **quasi-elastic** tail is now priced, not dropped; the tail is
> emitted at **every** `θ_S` via `P₂(cos θ_S)` (POLRAD Eq. (43)); and the band
> anchors moved to sourced values. Points where this revision *disagrees* with a
> review finding carry a **Reviewer note / response** line in place.

---

## Contents

1. [Physics](#1-physics)
2. [Data inputs](#2-data-inputs)
3. [API](#3-api)
4. [Output & bindings](#4-output--bindings)
5. [Tests](#5-tests)
6. [Implementation steps](#6-implementation-steps-two-opus-agents)
7. [Open questions](#7-open-questions)
8. [Placeholders the implementer must fill](#8-placeholders-the-implementer-must-fill)

---

## 1. Physics

### 1.0 What is and is not being corrected

Radiative corrections to inclusive `e + A` split (POLRAD 2.0 Eq. (2)) into

> **σ = σ^in + σ^el + σ^q + σ^v**   — Akushevich, Ilyichev, Shumeiko, Soroko,
> Tolkachev, *POLRAD 2.0*, [arXiv:hep-ph/9706516](https://arxiv.org/abs/hep-ph/9706516),
> CPC **104** (1997) 201, Eq. (2)

with `σ^in` the inelastic radiative tail (IRT), `σ^el` the elastic radiative
tail (ERT), `σ^q` the quasi-elastic tail (QRT) and `σ^v` the virtual (loop)
correction. This design treats them as three different objects:

| piece | what it does to `A_zz` | what LiPolGen does |
|---|---|---|
| `σ^v` + `σ^in` (lepton-vertex ISR/FSR, vacuum polarisation, soft/hard bremsstrahlung off the *inelastic* target) | **spin- and tensor-blind at the lepton vertex**: it multiplies `F₁`, `F₂` and `b₁…b₄` by the *same* radiator, so the bulk cancels in the ratio `A_zz` and only an `(x, Q²)` **shape migration** survives | **not implemented as a shift.** Its residual is exactly what the **band** of §1.3 prices. §1.6 says what implementing the shift would require. |
| `σ^el` (elastic radiative tail) | **does not cancel** — not because it "carries tensor dependence", but because **its own `A_zz` is not the Born's**. A background with `A_zz^tail ≠ A_zz^Born` under the DIS bin *dilutes* the measured asymmetry by `σ^tail_U/(σ^Born + σ^tail)`; the tail's own (small) tensor part then shifts it back a little. | **implemented**, POLRAD Eqs. (37)–(39) + (43) (the t-peak closed forms), as the additive weight `rc_tail` (§1.4) |
| `σ^q` (quasi-elastic tail) | **polarised part** neglected: *"there is no net tensor effect by inclusive scattering on weakly-bound spin-½ objects"* (**Z.-L. Zhou et al. (NIKHEF), PRL 82 (1999) 687** — a quasi-elastic `e-d` T₂₀ measurement, cited as HERMES ref. [13] in [hep-ex/0506018](https://arxiv.org/abs/hep-ex/0506018)). Its **unpolarised** part is `A_zz^tail ≈ 0`, so it **does NOT cancel** — it dilutes `A_zz` exactly as the unpolarised elastic tail does. HERMES subtracted **both** coherent and quasi-elastic tails, and *"the radiative background … reaches almost 50 % of the statistics in the lowest-x bin"*. | **unpolarised part implemented** via POLRAD Eq. (44) on the proton ERT (`σ^q_U`, folded into `rc_tail`); **polarised part not implemented**, citing Zhou et al. |

**Why the QRT cannot be waved away for ⁶Li.** The ERT scales as `Z² F_c²`
with a *nuclear* charge form factor that is dead by `t ≈ 0.25 GeV²`; the QRT
scales as `Z F_e² + N …` with *nucleon* form factors still alive at
`t ∼ 1 GeV²`. Against the `∫dt/t` (really `∫dη/η²`) weighting from
`t_min = (x M_N)²` the two are of the same order — a back-of-envelope 30–70 %
of the ERT over `x = 0.01–0.3`. Pricing the ERT alone and calling it "the
radiative-tail dilution" would be misleading, which is why `rc_tail` carries
both. See §1.4.7.

Consequence: the whole design is **multiplicative weights on a finished
event**, never a momentum shift. That is the FSI precedent
(`include/lipolgen/fsi.hpp` file header; `src/core/pipeline.cpp:469-480` builds
the model, `src/core/pipeline.cpp:3551` applies `ev.weight *= te.weight`), with
one deliberate difference: **the RC weights do not multiply `Event::weight`.**
FSI is a *correction to the model*, so it belongs on the nominal weight; the RC
band is a *systematic variation family* and the tail is a *background*, and
both must leave the Born sample alone. They travel as separate named weights
(§4) that an analysis multiplies in on purpose.

### 1.1 The rank-2 part of an event's rate — one definition for every channel

`include/lipolgen/xsec.hpp` writes the polarised inclusive rate as

```
  dsigma/(dx dQ2 dphi) = sigma_unpol(x,Q2)/(2 pi) * W(phi'),
  W(phi') = 1 + w_avg + a1 cos(phi') + a2 cos(2 phi'),   phi' = phi - phi_S
```

and (`src/core/xsec.cpp:255-267`, massless default path) the **b-sector** part
of `w_avg` is

```
  w_tensor = TENSOR_LL_SIGN * Q_NN * P2(cos theta_S) * K / D_phi
  Q_NN     = [3 m^2 - J(J+1)]/3                (InclusiveKernel::tensor_moments)
  K        = b1 + (1-y)/(x y^2) b2             (InclusiveKernel::tensor_kernel)
  D_phi    = F1 + (1-y)/(x y^2) F2             (InclusiveKernel::dphi)
```

With `TENSOR_LL_SIGN = -1` (`include/lipolgen/constants.hpp:40`) and
`azz()` (`src/core/asymmetries.cpp:79-87`)
`A_zz = TENSOR_LL_SIGN·(2/3)·(K/D_phi)·P2(cos θ_m)`, this is **identically**

```
  w_tensor = (P_zz^eff / 2) * A_zz ,     P_zz^eff = 3 Q_NN P2(cos theta_S)
```

— i.e. for a pure spin-1 state along the beam (`θ_S = 0`), `P_zz^eff = 3m²−2`
= +1 for `m = ±1` and −2 for `m = 0`, and the rate is `σ_U[1 + (P_zz/2)A_zz]`,
exactly the supervisor's form and exactly **HERMES Eq. (1)**
([hep-ex/0506018](https://arxiv.org/abs/hep-ex/0506018);
their Eq. (2) is the *extraction* formula `A_zz = (2σ¹ − 2σ⁰)/(3σ_U P_zz^eff)`).
**This identity is the design's anchor test** (§5, T1).

Two things must be said precisely, or the geometry is double counted and the
band leaks into the vector sector:

* **`A_zz` here means the LONGITUDINAL asymmetry**, i.e.
  `A_zz ≡ azz(b1, f1, f2, x, y)` at its default `θ_m = 0`
  (`src/core/asymmetries.cpp:80-88` (anchor read 2026-09-15)). All of the axis geometry sits in
  `P_zz^eff = 3 Q_NN P₂(cos θ_S)`. Substituting `θ_m = θ_S` into `azz()` *as
  well* would apply `P₂(cos θ_S)` twice.
* **The whole `A_zz` programme runs at `λ_e P_e = 0`** (an unpolarised beam).
  This is stated once here and assumed everywhere below: in the τ definition,
  in the tagged factorisation of §1.5.1 and in the tail denominator `1 + w_avg`
  of §1.4.7. `RcModel`'s constructor **throws** when any category of the run
  plan has `λ_e · P_e ≠ 0` (§3.1) rather than silently pricing the vector
  sector as if it were rank 2.

Define the **tensor fraction of an event** as the **rank-2 (b-sector)
projection of its own density**:

> **τ ≡ W_tensor(φ′) / W(φ′)**
> with `W_tensor` the `b₁…b₄` (`T_LL`) contribution to `W` alone — exactly what
> `InclusiveKernel::tensor_amplitudes` (§3.2) returns, and exactly what the new
> `StateTables::w_tensor/a1_tensor/a2_tensor` hold.

*Corollary, valid only at `P_e = 0`:* `τ = 1 − ρ̄/ρ`, with ρ the density at the
event's own spin state and ρ̄ the same density averaged over the `2J+1`
projections at everything else fixed (`ρ̄ = 1`, since both `Σ_M Q_NN = 0` and
`Σ_M m = 0`). **`1 − ρ̄/ρ` is not the definition** — it is "everything of rank
≥ 1", and with `λ_e P_e ≠ 0` it would sweep in the vector `A_par` (and `a₁
A_perp`) terms that §1.2 explicitly excludes. Use it as the closed-form
shortcut it is, not as the specification.

At `P_e = 0` and along the beam, τ reduces to
`(P_zz/2)A_zz / [1 + (P_zz/2)A_zz]` — the supervisor's `δ`-coefficient exactly.
§1.5 shows the *same* expression works on the tagged channels, where the tensor
structure lives in the cluster wave function and not in the DIS kernel.

The RC band is then, for a stated `δ(x)`,

> **w_± = 1 ± δ(x) · τ**   ⟺   **w_± = [1 + (P_zz/2) A_zz (1 ± δ)] / [1 + (P_zz/2) A_zz]**

with `w_lo = w_−` and `w_hi = w_+`. The two forms are algebraically identical;
the code computes the first, because it is the one that generalises.

### 1.2 What "the tensor part" contains — and what it deliberately does not

| term of `W(φ')` | rank | in the band by default? | why |
|---|---|---|---|
| `w_avg` b-sector (`b₁…b₄`, `T_LL`) | 2 | **yes** | POLRAD's ℑ₁,ℑ₂,ℑ₅–ℑ₈ and Gakh–Shekhovtsova compute exactly this |
| `a1`, `a2` from `tensor_harmonics_gamma` (`Options::tensor_gamma = true`) | 2 | **yes** | same b-sector, only re-projected onto φ' by the O(γ²) frame rotation; excluding it would make the band depend on a kinematic option |
| `a2` from `Δ` (gluon transversity, `a_cos2phi`) | 2 | **no** (opt-in `RcScope::TensorAll`) | **nobody has computed RC for a φ-dependent tensor observable** (`physics_literature.md` §3 verdict, verbatim). `Δ` is not in POLRAD's `b₁…b₄` basis. Keeping it out preserves the honest "no band" for the cos 2φ programme. |
| `w_avg` vector (`A_par`), `a1` vector (`A_perp`) | 1 | **no** | a vector RC band is a different systematic (`g₁`), out of scope for item 9 |

### 1.3 The band δ(x) — two anchors, a log-linear interpolation, and their provenance

The two anchors are **not the same kind of object**, and the doc says so
rather than pretending otherwise: the high-x one is a quoted *uncertainty*, the
low-x one is the *size of a correction that this design does not apply*, taken
as a 1σ band. That choice is stated in `CONVENTIONS.md` and repeated in the
banner.

* **High x — 1.5 %, at `x_high = 0.16`.** JLab **E12-13-011** (proposal
  **PR12-13-011**, Slifer et al.): *"no polarized radiative corrections at the
  lepton vertex, and the unpolarized corrections are known to better than
  1.5 %"*. **Provenance flag:** this quote and the accompanying systematics
  budget live in the **unpublished proposal only**. The published companion
  paper, Poudel, Bacchetta, Chen, Santiesteban,
  [arXiv:2506.04506](https://arxiv.org/abs/2506.04506), EPJ A **61** (2025),
  contains no radiative-correction discussion at all — the string "radiative"
  does not occur in it — so the earlier draft's *"RC is 1.5 % of a 9.2 % total
  (polarimetry 8 %, dilution 4 %)"* and its label *"published, adopted"* are
  **withdrawn** (they were neither in that paper nor internally consistent:
  √(8² + 4² + 1.5²) = 9.07). What 2506.04506 *does* fix is the experiment's
  **kinematic reach: `0.16 < x < 0.49`, `0.8 < Q² < 5.0 GeV²`, `W ≥ 1.85 GeV`**
  (p. 8). The anchor therefore sits at **`x_high = 0.16`**, the lower edge of
  the region the number was quoted for — the previous `x_high = 0.05` was
  outside E12-13-011's phase space and had no source. An implementer who finds
  the proposal should cite it **by page** for the 1.5 %, or drop the anchor.
* **Low x — 10–30 %, a CORRECTION SIZE, not a residual.** **G. I. Gakh,
  O. Shekhovtsova, [arXiv:hep-ph/0403262](https://arxiv.org/abs/hep-ph/0403262),
  JETP 99 (2004) 898**, §4 *Numerical estimations*, text, verbatim
  (**CORRECTED 2026-09-04**: this said "§5". The arXiv source has four numbered
  sections plus three appendices, and what the PDF renders as "§5" is
  Appendix A; `../run_2026-09-03/phase_B_numbers.md` §B6.3(d)): *"In the range of low x
  (x ∼ 10⁻³–10⁻²) the value of radiative correction changes from 10 % to 30 %
  as compared with the Born contribution"*, said of **the spin-dependent part of
  the cross section** for an unpolarised beam on a tensor-polarised deuteron —
  exactly our configuration — in HERMES kinematics, with the elastic tail
  included. Their Fig. 2 also shows the RC **shifting the `b₁` zero-crossing to
  smaller x**; its caption lists **four** panels,
  `a — 0.1, b — 1, c — 4, d — 10 GeV²` (the earlier draft said two). **Two
  honest flags, to be repeated wherever the number is quoted:** (i) this is a
  single calculation with **zero citations on INSPIRE**, never independently
  checked or used, and its two numbers are — **CORRECTED 2026-09-04**, measured
  off the paper's own arXiv figure arrays — **the two ends of ONE panel's `x`
  window at a single `Q² = 0.1`**: `|δ| = 0.1133` at `x = 0.00966` and `0.2662`
  at `x = 0.00226`, which the sentence quoted just above scopes to **x**
  ("in the range of low x"), not to `Q²`. They are **not** "the spread of *one*
  calculation over `Q²`, not two independent edges" — the wording this bullet
  carried until then, now **withdrawn** — and the consequence is worth its own
  line: **`δ_low = 0.30` is anchored at `x_low = 0.01`, which lies ABOVE that
  panel's top `x = 0.00966`, where the paper itself reads `|δ| = 0.113`**; the
  `0.266` that rounds to "30 %" sits at `x = 0.00226`, a factor **4.3 lower in
  x**. So the anchor's real basis is *the panel's lowest-`x` value carried
  upward in `x`*, not *the conservative end of a `Q²` spread*
  (`../run_2026-09-03/phase_B_numbers.md` §B6.3); (ii) **10–30 % is how big the
  correction is, not how well it is known.** Using it as a 1σ band is the *choice* this design
  makes — the choice a generator that does not apply the correction has to make
  — and it is deliberately conservative.
* **What a real experiment's residual actually was.** HERMES
  ([hep-ex/0506018](https://arxiv.org/abs/hep-ex/0506018), PRL **95** (2005)
  242001) quotes a residual RC systematic on `A_zz` of **≈ 2×10⁻³ "for the three
  bins at low x" and negligible at high x**. The measured `A_zz` in those three
  bins (their Table II, ×10⁻²) are **−1.06 at ⟨x⟩ = 0.012, −1.07 at 0.032,
  −1.32 at 0.063** — so the *fractional* residual is **18.9 %, 18.7 %, 15.2 %**,
  not the 10 % the earlier draft inferred by dividing 2×10⁻³ by the range
  maximum `|A_zz| ≤ 0.02` (which occurs at x = 0.45, where HERMES says RC is
  negligible). `RC_DELTA_LOW_X_OPTIMISTIC` is therefore **0.19**, the value at
  the lowest-x bin where the anchor sits, with 0.15–0.19 quoted as its range.
  Note this **contradicts a δ of 1.5 % anywhere below x ≈ 0.1**: at x = 0.063
  the only measured residual in the literature is 15 %.

**The interpolation.**

```
                 { delta_high                                              x >= x_high
  delta(x)   =   { delta_high + (delta_low - delta_high) * ln(x_high/x)/ln(x_high/x_low)
                 { delta_low                                               x <= x_low

  defaults:  delta_high = 0.015   at   x_high = 0.16     (PR12-13-011, at ITS OWN lower x edge)
             delta_low  = 0.30    at   x_low  = 0.01     (Gakh-Shekhovtsova, conservative edge)
```

This interpolation now passes through **δ(0.063) = 0.1108**, i.e. it is no
longer in conflict with HERMES's measured 15 % residual at that x — which the
old `x_high = 0.05` version was (it gave 1.5 % there).

Log-linear in `x` because both the anchors and the underlying physics (the
soft-photon `ln(1/x)` and the growth of the radiative background towards the
kinematic edge) are logarithmic in `x`. Clamped outside, so `δ` is **monotone
non-increasing in x with both anchor values exact** (§5, T3).

**Why `lo`/`hi` are the two SIGNS of one δ and not the 10 % and 30 % edges.**
The supervisor's own formula is `w_± = [1 + (P_zz/2)A_zz(1 ± δ)]/[...]`, i.e.
`±` is the direction of the variation. A one-sided pair (10 %, 30 %) would not
be a *band on `A_zz`* — it would be two upward variations, and an analysis
could not symmetrise it. Defaulting `δ_low = 0.30` makes the pair
`(w_lo, w_hi)` **contain** the `δ_low = 0.10` result, so nothing is lost; the
optimistic edge is a one-switch run (`--rc-delta-low-x 0.19`,
`RcOptions::delta_low_x`) and both must be quoted in any writeup, as with the
FSI `σ_XN` band (`fsi.hpp`: *"Band it; never quote one row alone."*).

**δ table at the defaults** (implementer: pin these in the doctest, at the
**full double precision printed here** — T3 asserts 1e-12, which a 5-digit
table cannot meet):

| x | δ(x) |
|---|---|
| 0.001 | 0.3 (clamp, exact `==`) |
| 0.005 | 0.3 (clamp, exact `==`) |
| **0.01** | **0.3** (anchor, exact `==`) |
| 0.02 | 0.22874999999999995 |
| 0.03 | 0.18707142182361763 |
| 0.05 | 0.13456262323927543 |
| 0.063 | 0.1108061822113555 |
| 0.08 | 0.08625000000000000 |
| 0.1 | 0.06331262323927542 |
| **0.16** | **0.015** (anchor, exact `==`) |
| 0.3 | 0.015 (clamp, exact `==`) |

The four `==` rows are legitimate `==` because they take the **clamp branch**
(`x ≤ x_low`, `x ≥ x_high`) and return the constant unmodified; the
interpolated rows are pinned at 1e-12 against the printed doubles.

### 1.4 The radiative tails — POLRAD 2.0, the t-peak closed forms

This is the piece that **does not cancel in `A_zz`**: a background under the DIS
bin from `e + ⁶Li → e' + γ + ⁶Li(g.s.)` (elastic) and
`e + ⁶Li → e' + γ + N + ⁵Li/⁵He` (quasi-elastic), whose own `A_zz` is not the
Born's — mostly because it is ≈ 0. It **dilutes**; the tail's small tensor part
then shifts the dilution slightly.

**Design change from the first draft, and why.** The first draft built the tail
from POLRAD **Eq. (18)** (the general `τ_A` quadrature) with the **Appendix-B**
`θ_ij` machinery and the **Eq. (A.4)** generalised structure functions, and
proposed to fix its normalisation on a leading-log equivalent-radiator limit.
That is the wrong gate, and POLRAD says so on p. 10, §2.1.3 B, verbatim:

> *"In the case of ultrarelativistic approximation in calculation of RC from
> elastic radiative tail, due to strong dependence of formfactors `F_i` on the
> square of transfer momenta `Q²`, the leading contribution to the total cross
> section gives only **t−peak (Compton peak)** and contributions of s− and p−
> peaks are **suppressed** [14]."*

and again in §2.1.4: the t-peak *"is extremely important in the cases of elastic
and quasielastic radiative tails"*. The s-peak sits at `t ≈ (1−y)Q²`, where the
⁶Li charge form factor is dead for any `y` not → 1; the **t-peak** reaches down
to `t_min = 4M_A² η_min = M_A² x_A²/(1−x_A) ≈ (x M_N)²`, i.e. **inside** the
form factor at every `(x, Q²)`. So the collinear/equivalent-radiator limit is a
*sub-dominant* piece of Eq. (18), it cannot agree with it to 5 %, and a gate
built on it fixes nothing.

**Therefore v0 ships POLRAD's own t-peak closed forms — Eqs. (37)–(39)
plus (43) — which are ONE `η_A` integral and need no Appendix B at all.** The
full Eq. (18) + Appendix-B transcription is retained below as the documented
**upgrade path** (§1.4.6), gated behind `RcTailModel::PolradFull`, and is where
the `+5 person-days` of `physics_literature.md` §3 actually go.

#### 1.4.1 Invariants (POLRAD §2.1.1, Eq. (4); §2.1.2, Eqs. (14), (17))

| symbol | definition | note |
|---|---|---|
| `k₁, k₂` | incoming / outgoing lepton four-momenta | |
| `p` (`p_A`) | target four-momentum; `p² = M²` (`p_A² = M_A²`) | POLRAD writes *nucleon* invariants without the index and *nucleus* invariants with `A` |
| `m` | lepton mass | kept: it is what makes `F_IR` finite |
| `S = 2 k₁ p`, `X = 2 k₂ p = (1−y)S` | | POLRAD Eq. (4) |
| `Q² = −(k₁−k₂)² = x y S` | | |
| `S_x = S − X`, `S_p = S + X`, `Q_m² = Q² + 2m²`, `λ_s = S² − 4m²M²` | | Eq. (4) |
| `k` | **real (radiated) photon** four-momentum | |
| `R = 2 p k` | | integration variable of Eq. (13) |
| `τ = k(k₁−k₂)/(p k)` | | second integration variable |
| `λ_Q = S_x² + 4M²Q²`, `W² = S_x − Q² + M²` | | Eq. (14) |
| `τ_max,min = (S_x ± √λ_Q)/(2M²)` | | Eq. (14) |
| `R_el = (S_xA − Q²)/(1 + τ_A)` | **the elastic constraint** — the nucleus stays in the ground state, so `R` is no longer free | Eq. (17) |
| `η_A = t/(4M_A²) = (Q² + R_el τ_A)/(4M_A²)` | the elastic-vertex momentum transfer, in units of `4M_A²` | stated under Eq. (A.4) |
| `P_N`, `Q_N` | target **vector** and **tensor** polarisation degrees | `Q_N` is the object we map to `P_zz` — see §1.4.4 |
| `P_L` | initial-lepton polarisation degree | **0** in the whole `A_zz` programme |
| **the nuclear map** | `S_A = 2k₁p_A = A·s_per-nucleon`, `X_A = (1−y)S_A`, `x_A = Q²/(y S_A) = x/A`, `y` unchanged, `M → M_A` | **Every** POLRAD elastic-tail formula (Eqs. (18), (37)–(39) and every Appendix-A/B function inside them) is evaluated in these *nuclear* invariants. The source is POLRAD **Eq. (21)**, whose elastic exponentiation factor reads `(y²(1−x/A)²/(1−xy/A))^{t_r}` against the quasi-elastic `(y²(1−x)²/(1−xy))^{t_r}` — so the elastic formulas' `x` is `x_A = x/A`. |
| `M_A` | **the AME/`beams.cpp` nuclear mass**, `M(⁶Li) = 5.601518702 GeV` | **not** `A·M_NUCLEON` (which would be 5.6296 GeV, 0.5 % high and would move `η_A` by 1 %). Passed in from `Ion::mass()`, never retyped. |

#### 1.4.2 The v0 tail — POLRAD's t-peak closed forms, Eqs. (37)–(39), (43)

POLRAD Eq. (37) splits the whole elastic tail by polarisation degree:

```
  sigma_1^el  =  sigma_u^A  +  P_L P_N sigma_p^A  +  (Q_N / 6) sigma_q^A       (37)
```

**Note the explicit `Q_N/6`** — the tensor tail's normalisation is fixed here,
not chosen. With `P_L = 0` the whole `A_zz` programme needs `σ_u^A` and
`σ_q^A` only; `σ_p^A` (and Eq. (41)/(42), the `M/√S`-suppressed transverse
vector piece) is never evaluated.

Eq. (38), transcribed for the **deuteron** (spin 1 — the shape ⁶Li inherits)
and for **carbon** (spin 0 — the `Z²` gate):

```
                 alpha^3        INF  d eta_A [                                                ]
  sigma_u^d  =  --------- Y_+   INT  ------- [ (F_c^2 + (8/9) F_q^2 eta_A^2 + (2/3) F_m^2 eta_A) Xt
                    S          eta_m  eta_A  [                          - (2/3)(1 + eta_A) F_m^2 ]

                 alpha^3        INF  d eta_A     {                                                    }
  sigma_q^d  =  --------- Y_+   INT  -------  {  ( 1 + eta_A + ((3/4) x^2 - eta_A) Xt ) F_m^2
                    S          eta_m  eta_A     {
                                                 - ( Xt X_1 / (1 + eta_A) ) F_q ( 3 F_c + 3 eta_A F_m + eta_A F_q )
                                                 - 2 eta_A Xt F_q ( 4 F_c - 3 x F_m + (4/3) eta_A F_q )  }

                 alpha^3            INF  d eta_A
  sigma_u^C  =  --------- Z^2 Y_-   INT  ------- Xt F^2                                          (38)
                    S              eta_m  eta_A
```

with (Eq. (39))

```
  X_1  = x^2 + 4 x eta_A - 4 eta_A          (vanishes identically at eta_A = eta_min)
  Xt   = X_1 / (2 eta_A x^2)
  Y_+- = ( 1 +- (1 - y)^2 ) / (1 - y)
  eta_A   = t / (4 M_A^2)
  eta_min = x^2 / (4 (1 - x))               =>  t_min = M_A^2 x^2/(1-x) ~= (x M_N)^2
```

and, for a spin-½ nucleus, `F₁ = (G_E + η_A G_M)/(1 + η_A)`,
`F₂ = (G_M − G_E)/(1 + η_A)` (Eq. (40)).

**Every `x` in Eqs. (38)/(39) is the NUCLEAR `x_A = x/A`** and every `M` is
`M_A`, per the §1.4.1 map. That is what makes `t_min = 4M_A² η_min ≈ (x M_N)²`
independent of `A`: `x_A M_A = x M_N`.

**Three transcription flags the implementer must clear before trusting a
number** (they are cheap, and Q2 makes the FORTRAN check a *prerequisite*, not
a test):

1. `σ_u^C` is printed with **`Y₋`**, not `Y₊`, unlike every other unpolarised
   entry in Eq. (38). Transcribe it as printed, then check it against the
   FORTRAN; if it is a typo the carbon gate below must use `Y₊`.
2. `X₁ = x² + 4xη_A − 4η_A` — the last term is **linear** in `η_A`. The check
   that fixes it: `X₁(η_min) = x² + 4η_min(x−1) = 0` identically, which is what
   makes `X̃` vanish at the lower limit and the `∫dη_A/η_A` converge there.
   A transcription with `4η_A²` fails that identity.
3. The two fractions inside `σ_q^d` are `(3/4)x²` (inside the `X̃` bracket) and
   `(4/3)η_A F_q` (in the last line). They are *not* the same fraction and
   `pdftotext` renders both as bare digit pairs.

**The transverse and general axis, Eq. (43):**

```
  sigma_q,perp^d  =  -(1/2) sigma_q,parallel^d                                  (43)
```

i.e. `P₂(cos 90°) : P₂(cos 0°) = −½ : 1`. §1.4.7 turns that into the general
`θ_S` rule.

#### 1.4.3 The quasi-elastic tail, POLRAD Eq. (44) — unpolarised only

POLRAD builds the QRT in the same t-peak form by **reusing the proton elastic
tail** with three substitutions:

```
  F_e^2   ->  S_e(Q^2)  F_e^2 ,
  F_m^2   ->  S_m(Q^2)  F_m^2 ,
  F_e F_m ->  S_em(Q^2) F_e F_m                                                 (44)
```

`S_{E,M,EM}` are the Pauli/sum-rule **suppression factors** (POLRAD §2.1.2,
Eq. (20) and its ref. [10]). So:

```
  sigma_u^q(A)  =  Z * sigma_u^p[F_e -> G_E^p, F_m -> G_M^p, S_*]
                +  N * sigma_u^p[F_e -> G_E^n, F_m -> G_M^n, S_*]
```

with `σ_u^p` the first line of Eq. (38) (in nucleon invariants: `M → M_N`,
`x_A → x`, `η_min = x²/(4(1−x))`, i.e. the **same** `t_min = (x M_N)²`), and
`(Z, N) = (3, 3)` for ⁶Li.

~~**v0 takes `S_e = S_m = S_em = 1`** (no Pauli suppression) and exposes
`RcOptions::qe_suppression` as a flat multiplier on the whole QRT, to be run at
`0.0 / 0.5 / 1.0` as a band.~~

> **Q9 CLOSED 2026-09-03.** `S = 1` stopped being a small choice the moment the
> per-nucleon factor was fixed: the QRT is then **0.6× / 3.2× / 1000×** the ERT
> at `x = 0.01 / 0.1 / 0.3`, i.e. `rc_tail` is quasi-elastic-**dominated**
> everywhere. `src/core/rc.cpp` now implements `S_E = S_M = S(q)` exactly as
> POLRAD's `ffquas` codes it (adgh:5605-5613), the de Forest–Walecka Fermi-gas
> factor `S(q) = (3/4)(q/k_F) − (q/k_F)³/16` below `q = 2k_F`, at the
> elastic-vertex three-momentum transfer `q² = t(1 + η)`, with
> `RcOptions::qe_kf_gev = 0.169 GeV` (⁶Li, Moniz et al., PRL **26** (1971) 445)
> as the **default**. Measured `QRT(S)/QRT(1)`: **0.472** (`k_F = 0.169`) /
> 0.405 (0.221) at `x = 0.01`, **0.868** / 0.783 at `x = 0.1`, **1.000** /
> 0.990 at `x = 0.3` — so the old default was 15–60 % high on the tail's
> dominant piece at `x ≤ 0.1`. `qe_suppression` survives as the flat band knob
> *around* it (run `0.0 / 0.5 / 1.0`), and `qe_kf_gev = 0` recovers the old
> `S = 1` edge. `S_EM` never enters: the Eq. (38) spin-½ integrand is a
> combination of `G_E²` and `G_M²` alone.

~~**The polarised QRT stays at zero**, citing Z.-L. Zhou et al. (§1.0). That is
the one place the *tensor* part of the tail is knowingly incomplete, and it is
also the smallest part of it.~~

> **CORRECTED 2026-09-04 (task B3). "Also the smallest part of it" is
> REFUTED by this document's own §8.1 measurement.** After the per-nucleon fix
> and Q9's Pauli factor the quasi-elastic term is **22 % / 73 % / 99.9 %** of
> `rc_tail` at `x = 0.01 / 0.10 / 0.30` (`σ^q_U/σ^el_U` = 0.28 / 2.75 / 1000).
> Setting its tensor part to zero is therefore not a small omission but the
> **largest unpriced piece of the tail**, and the sentence above was written
> before the piece it dismisses grew by a factor `A = 6`.
>
> The polarised QRT is still **not computed** — POLRAD supplies no tensor
> partner to Eq. (44), and no polarised quasi-elastic radiative-tail
> calculation exists for an A = 6 spin-1 nucleus — but it is now **priced**,
> opt-in, by `RcOptions::qe_tensor_scale` (default **0.0**, the shipped
> tensor-blind tail bit for bit), which lends the quasi-elastic tail the
> **elastic** tail's own `σ^el_T/σ^el_U`. That is a **borrowed magnitude, not
> a derived bound**: it puts a *coherent* nuclear quadrupole fraction on an
> *incoherent* nucleon process, and ⁶Li's elastic tensor fraction is
> anomalously small for a reason the quasi-elastic piece has no reason to
> share, so scale = 1 may be ~10² too small at `x ≤ 0.1`. Measured sizes,
> caveats and the full argument: `rc.hpp`'s `RcOptions::qe_tensor_scale` and
> `../run_2026-09-03/phase_B_numbers.md` §B3.

**Why this is in the design at all.** The QRT is not a rounding error: its
nucleon form factors are alive at `t ∼ 1 GeV²` where ⁶Li's `Z² F_c²` is long
dead, and against the `∫dη_A/η_A` weighting from `t_min = (x M_N)²` the two are
comparable — a first estimate puts the QRT at **30–70 % of the ERT** over
`x = 0.01–0.3`. HERMES subtracted both. A `rc_tail` carrying only the elastic
piece would be read as "the radiative-tail dilution" and would be wrong by that
factor. §8.1 measures the ratio and replaces this estimate.

#### 1.4.4 The Born normalisation and the `Q_N ↔ P_zz` map — POLRAD Eqs. (9)/(10)

```
  d sigma / dx dy  =  (4 pi alpha^2 S / Q^4) [ (F1 - (Q_N/3) b1) x y^2
                                             + (F2 - (Q_N/3) b2) (1 - y)
                                             - P_L P_N x y (2 - y) g1 ]              (9)   longitudinal target

  d sigma / dx dy  =  (4 pi alpha^2 S / Q^4) [ (F1 + (Q_N/6) b1) x y^2
                                             + (F2 + (Q_N/6) b2) (1 - y)
                                             - 2 P_L P_N (x sqrt(x y (1-y)) M / sqrt(S)) (y g1 + 2 g2) ]  (10)  transverse
```

Two things follow, and both are already in the repository:

1. **The sign and the geometry.** Eq. (9) gives `σ/σ_U = 1 − (P_zz/3)b₁/F₁ =
   1 + (P_zz/2)A_zz ⟹ A_zz = −(2/3)b₁/F₁` with `Q_N = P_zz`, which is the
   fourth independent confirmation of `TENSOR_LL_SIGN = −1`
   (`docs/CONVENTIONS.md`; `OPEN_ITEMS_SOLUTIONS.md` §6(d)). And the
   **−1/3 : +1/6 ratio between Eqs. (9) and (10) is exactly
   `P₂(cos 0°) : P₂(cos 90°) = 1 : −½`**, i.e. POLRAD's `Q_N` is our
   `P_zz^eff = 3 Q_NN P₂(cos θ_S)` — the *same* rank-2 geometry
   `InclusiveKernel::tensor_moments` implements. **Eq. (43) says the same thing
   for the TAIL** (`σ_q⊥ = −½σ_q∥`), which is what lets §1.4.7 emit `rc_tail` at
   every axis. So **one geometry drives both the band and the tail**, and the
   design defines `Q_N` nowhere new: `Q_N ≡ P_zz^eff` of §1.1. Note the tail
   carries it as **`Q_N/6`** (Eq. (37)), not `Q_N` — that factor is fixed by
   POLRAD, not chosen here, and dropping it inflates the tensor tail sixfold
   with no unpolarised test seeing it.
2. **The Born is not re-implemented.** `docs/CONVENTIONS.md` forbids defining a
   physics number twice. The unpolarised part of Eq. (9),
   `(4πα²S/Q⁴)[x y² F₁ + (1−y)F₂]`, is **identically** the library's own
   `InclusiveKernel::dsigma_unpol` after `d²σ/dx dy = x s · d²σ/dx dQ²`, so
   `rc.cpp` uses `dsigma_unpol(x, q2, s) * x * s` and multiplies by the
   library's own `1 + w_avg` for the polarised denominator. A doctest pins the
   equality (§5, T7).

#### 1.4.5 The event weight

```
                            sigma^el_U(x,Q2) + (Q_N/6) sigma^el_T(x,Q2) + kappa_qe * sigma^q_U(x,Q2)
  w_tail(event) = 1 + ----------------------------------------------------------------------------------
                                  [ x * s * dsigma_unpol(x, Q2, s) ] * ( 1 + w_avg )
```

with `Q_N = P_zz^eff` of the event's spin state (§1.4.7 for `θ_S ≠ 0`),
`σ^el_U ≡ σ_u^A`, `σ^el_T ≡ σ_q^A` of Eq. (38), `σ^q_U` the Eq. (44) QRT of
§1.4.3, and `κ_qe = RcOptions::qe_suppression`.

**Units and the per-nucleon convention — stated, not left to a test.** The
library is per nucleon throughout (`CONVENTIONS.md`; `xsec.hpp:313`), so the
**denominator** is the per-nucleon Born `d²σ/dx dy` in per-nucleon `x` and `s`.
The **numerator** must be the same object. POLRAD's tail formulas are per
nucleus in `x_A = x/A` with `S_A = A·s`, and

```
  d^2 sigma / (dx dy)  =  (1/A) d^2 sigma / (dx_A dy)      [ since x = A x_A ]
```

so the per-nucleon tail is `(1/A)·[Eq. (38) evaluated at (x_A, y, S_A, M_A)]`.

> **CORRECTED 2026-09-03 (physics review).** The line above is **half of the
> answer and shipping it alone was a factor-`A` bug.** `(1/A) d²σ/(dx_A dy)` is
> the **Jacobian** `dx_A/dx`; it converts the *whole-nucleus* rate in `x_A` to
> the *whole-nucleus* rate in `x`. Going per **nucleon** needs a **second**
> `1/A`. Eq. (38) is the whole-nucleus `d²σ/(dx_A dy)` — the paper's
> `σ₁^el = (1/A) d²σ/dx_A dy` (polrad2t.tex:579) is loose notation, POLRAD's
> FORTRAN applies **both** `ter = m_p/M_A` (`apptai`, adgh:8607-8620) **and**
> `/tara` (main, adgh:489, 497), and a Weizsäcker–Williams × Compton
> construction reproduces `−`Eq. (38) with **no** `1/A` in it
> (`polrad_transcription_check.md` §8b). So
>
> > **per-nucleon tail = `(1/A) · (dx_A/dx) · [Eq. (38)] = [Eq. (38)] / A²`**
>
> with this document's own map `x_A = x/A`. `src/core/rc.cpp` applies `1/A²`
> and `tests/test_rc.cpp` T8(a') pins it at 1e-12, together with a gate that
> the elastic and quasi-elastic per-nucleon factors agree (v0 shipped
> `m_p/M_A` on the elastic tail and `1/A` on the quasi-elastic one, which was
> 6.00× too large on the former and internally inconsistent between them).
> `polrad_transcription_check.md` §4's "net difference 0.5 %" sentence, which
> the implementer cited to justify the single factor, is corrected in place.
Both numerator and denominator are then per nucleon, in **pb** after the
library's own `GEV2_TO_PB`. **A missing factor `A = 6` here is exactly the kind
of error a loosely specified 5 % test would absorb into a compensating
prefactor**, which is why §5's tail gates are written against POLRAD's own
Eq. (38) rather than against a radiator the implementer also writes.

The numerator is the tail *into the same reconstructed `(x, Q²)` bin*: the
design's stated approximation is that the reconstructed `(x, Q²)` of a radiative
elastic event is the one computed from the scattered electron alone, which is
what POLRAD's `(x, y)` already means (POLRAD §2: *"on a level of RC radiated
real photon momentum is indefinite, hence ν and Q² are arbitrary"* — the tail is
by construction quoted at the **observed** `(x, y)`).

**Precomputation, and where it is NOT good enough.** The `η_A` quadrature is too
slow per event, so `RcModel` builds three tables at construction —
`σ^el_U`, `σ^el_T`, `σ^q_U` — over the sampler's accepted `(x, Q²)` cell grid
(`InclusiveSampler::x_cells()/q2_cells()`, 100 × 72). **But the tail is not as
smooth across a cell as the modulation amplitudes are**, and the first draft's
claim that "the weight cannot be more precise than the sample and cannot be
less" is wrong for it: the tail carries `Y₊ = [1 + (1−y)²]/(1−y)`, which
**diverges as `1/(1−y)`**, and a lower limit `t_min = (x M_N)²` that moves
steeply with `x`. Inside an edge cell at `y → y_max = 0.985` a cell-centre value
can be off by `O(1)`.

Therefore:

* tabulate the three `η_A` integrals **per cell corner**, and evaluate the tail
  at the **event's own `(x, Q²)`** by bilinear interpolation in
  `(ln x, ln Q²)` of `ln σ_tail` (the tail is smooth in the logs away from
  `y → 1`);
* refine the `y` direction: add extra table rows at `y ∈ {0.9, 0.95, 0.97,
  0.98, 0.985}` so the `1/(1−y)` growth is resolved rather than extrapolated;
* report `clipped_cell_fraction()` **per `y`-band** (`y < 0.5`, `0.5–0.9`,
  `> 0.9`) as well as globally, so a run cannot hide a clipped edge inside a
  small global number.

Construction cost: `≈ 3 × n_nodes` one-dimensional Gauss–Legendre quadratures
(≈ 128 nodes each, no Appendix B) — well under a second, once, immutable
afterwards and therefore thread-safe like `GlauberFsiWeight`.

#### 1.4.6 The upgrade path — Eq. (18) + Appendix B + Eq. (A.4)

> **IMPLEMENTED 2026-09-06 (task B2).** This section was written as a
> specification for something that threw; `RcTailModel::PolradFull` now runs,
> opt-in as `--rc-tail-model polrad-full`, and is measured in
> `../run_2026-09-06/phase_B_numbers.md` §B2. **What the section below got
> right:** the ingredient list, `x_A = x/A` in every Appendix-A/B function, the
> `1/A²` per-nucleon reduction (now *measured* — Eq. (18)/Eq. (38) → 1 as
> `x_A → 0`, so both are the whole-nucleus rate), and `ℑ^el_6`'s corrected
> `4/(1+η_A)`. **What it could not have got right, because Appendix B was
> never transcribed here:** `polrad2t.tex`'s Appendix B is wrong in **five**
> places against POLRAD's own FORTRAN, and the code follows the FORTRAN — the
> `a_ik` M-powers (Eq. (B.3)), the level `q_ik` acts at (Eq. (B.7)), the
> `T_821` pairing (Eq. (B.4)), and two terms of Eq. (B.8)'s second lift. Each
> was found by the `x_A → 0` gate, not by inspection; §B2.2 tabulates what
> each wrong reading gives.

`RcTailModel::PolradFull` (specified here first, implemented 2026-09-06).
Eq. (18) as printed on POLRAD p. 7:

```
                1  d^2 sigma^el      alpha^3 y   tau_A,max      8   k_i                      2 M_A^2 R_el^{j-2}
  sigma_1^el = ---  -------------  = - --------- INT   dtau_A  SUM SUM  theta_ij(tau_A) ------------------------------ Im^el_i(R_el, tau_A)
                A    dx_A dy             A^2     tau_A,min     i=1 j=1                  (1+tau_A)(Q^2 + R_el tau_A)^2
```

**The prefactor is `1/A²` on the RIGHT and `1/A` on the LEFT** — not "`1/A` on
both sides" as the first draft printed and as Q2 described it.

> **CORRECTED 2026-09-03.** What follows below — *"`σ₁^el` is already
> `(1/A)d²σ/dx_A dy`, i.e. already the per-nucleon-normalised object"* — is
> the reading that made §1.4.5's single `1/A` look sufficient, and it is
> **wrong**: the left-hand `1/A` is the `dx_A/dx` Jacobian, not a per-nucleon
> division, and Eq. (38) (which is what v0 actually ships) carries neither.
> See §1.4.5's correction box and `polrad_transcription_check.md` §8b. Compare Eq. (19)
(the QRT) and Eq. (20), which carry `α³y/A`, with no square. The two are
consistent with §1.4.5's per-nucleon reduction: `σ₁^el` is *already*
`(1/A)d²σ/dx_A dy`, i.e. already the per-nucleon-normalised object, and the
extra `1/A` on the right is what turns the whole-nucleus `dτ_A` integral into
it. This is settled by the printed equation, by the `A = 1` limit and by the
carbon `Z²` gate of §5 — **not** by tuning against a radiator.

Ingredients, for `i ∈ {1, 2, 5, 6, 7, 8}` only (`ℑ₃`, `ℑ₄` enter Eq. (3) and
every tail formula through `m M P_L` alone and vanish at `P_L = 0`, so the
longest block of Appendix B — `T_{3jk}`, `T_{4jk}` — need not be transcribed):

* `k_i = (3, 3, 4, 5, 5, 5, 3, 4)`, `l_i = (1, 1, 1, 2, 3, 3, 1, 2)` — Eq. (B.2);
* `θ_ij(τ) = Σ_k a_ik T_ijk(τ)` — Eq. (B.1), summed over
  `k = max(1, j + l_i − k_i) … min(j, l_i)`;
* `a_ik` — Eq. (B.3): `1` for `k = 1, i = 1,2,3,7`; `{η·q, −1}` for
  `k = {1,2}, i = 4,8`; `{Q² − 3(η·q)², 6 η·q, −3}` for `k = {1,2,3}, i = 5,6`;
* `T_{ij1}` — Eq. (B.4); `T_{ij2}`, `T_{ij3}` — the substitution rule Eq. (B.6);
* **Eq. (B.7)**, the `q_ik = δ_{k2}(δ_{i5} + δ_{i6}) τ/M` term that Eq. (B.6)
  uses for `i = 5, 6`, `k = 2` — i.e. *inside* the tensor sector we keep;
* **Eq. (B.8)**, the upper-index `ξ/η` `F`-functions that Eq. (B.6) substitutes;
* `F, F_IR, F_d, F_{1+}, F_{2±}, F_i, F_ii` — Eq. (B.12) on `B_{1,2}(τ)`,
  `C_{1,2}(τ)` of Eq. (B.13), with the numerically-stable Eq. (B.14) for `F_d`
  at `τ = 0`;
* `a_η, b_η, c_η` — Eq. (B.10); `s_η, r_η` — Eq. (B.9); the scalar products —
  Eq. (B.11).
* **Eq. (B.5) saves half of it:** `T_{5j1} = T_{1j1}` and `T_{6j1} = T_{2j1}`,
  the two purely-unpolarised kernels, so half the tensor sector re-uses the
  unpolarised transcription verbatim.

**Every one of these functions is evaluated in the NUCLEAR invariants** of the
§1.4.1 map (`S_A`, `X_A`, `M_A`, `x_A = x/A`), not the per-nucleon ones.

The Eq. (A.4) spin-1 generalised structure functions, transcribed for reference
(`η_A = (Q² + R_el τ_A)/(4M_A²)`):

```
  Im^el_1 = (1/6) eta_A F_m^2 ( 4(1 + eta_A) + eta_A Q_N )

  Im^el_2 = ( F_c^2 + (2/3) eta_A F_m^2 + (8/9) eta_A^2 F_q^2 )
            + (Q_N/6) [ eta_A F_m^2 + (4 eta_A^2/(1 + eta_A)) ( (eta_A/3) F_q + F_c - F_m ) F_q ]

  Im^el_3 = -(P_N/2) (1 + eta_A) F_m ( (eta_A/3) F_q + F_c )          [ P_L = 0 => unused ]
  Im^el_4 =  (P_N/4) F_m ( (1/2) F_m - F_c - (eta_A/3) F_q )          [ P_L = 0 => unused ]

  Im^el_5 = (Q_N/24) F_m^2
  Im^el_6 = (Q_N/24) ( F_m^2 + (4/(1 + eta_A)) ( (eta_A/3) F_q + F_c + eta_A F_m ) F_q )
  Im^el_7 = (Q_N/6) eta_A (1 + eta_A) F_m^2
  Im^el_8 = -(Q_N/6) eta_A F_m ( F_m + 2 F_q )
```

> **CORRECTED 2026-09-04 (task B5). `ℑ^el_6` above carried a spurious `η_A` in
> the `F_q` numerator and has been fixed in place.** What this document printed
> from 2026-09-02 until now was
>
> ```
>   Im^el_6 = (Q_N/24) ( F_m^2 + (4 eta_A/(1 + eta_A)) ( ... ) F_q )     [ WRONG ]
> ```
>
> The source has `4/(1+η_A)`, not `4η_A/(1+η_A)` — `polrad2t.tex` lines
> 2705–2706, verbatim:
>
> ```tex
> \Im ^{el}_{6}= {Q_N\over 24} \pmatrix{F_m^{2}
>     +{4\over 1+\eta_A}({\eta_A\over 3}F_q+F_c+\eta_A F_m)F_q} ,
> ```
>
> and POLRAD's FORTRAN says the same independently. `ℑ₆ = (Q_N/6)ε³(b₂/3+b₃+b₄)`
> (Eq. (A.3)) with `strf`'s elastic `b₂, b₃, b₄` (adgh:4235-4242) and `ε = 1/(2τ)`:
> the `f_m²` pieces collapse as `(4/3)τ²[1 − (3τ+2) + (6τ+1)] = 4τ³` and the `f_q`
> pieces as `(16/9)(τ³/τ₁)[3τF_q + 9F_c + 9τF_m]F_q = (16/3)(τ³/τ₁)(τF_q+3F_c+3τF_m)F_q`,
> so `ε³(b₂/3+b₃+b₄) = (1/2)[F_m² + (4/(1+η))((η/3)F_q + F_c + ηF_m)F_q]·2M_A²`
> — the `4/(1+η_A)` of the source, with no second `η_A` anywhere.
> (`ℑ₂`'s `4η_A²/(1+η_A)` **is** correct as printed; only `ℑ₆` was wrong.)
>
> **The CODE never had the error, because the code never had `ℑ^el_6`.**
> Eq. (A.4) at `Q_N ≠ 0` was evaluated only on the `RcTailModel::PolradFull`
> path, which refused with `NOT IMPLEMENTED` until 2026-09-06. What v0 shipped
> is Eq. (38) (`polrad_sigma_el_u` / the `σ_q` integrand) and the `Q_N = 0`
> Rosenbluth pair (`rosenbluth_spin1`) — both independently checked, neither
> containing `ℑ₆`. So this was a **document-only** defect, found by
> `polrad_transcription_check.md` §7.1 and corrected here; nothing in the
> shipped tail moved, and no test number changed.
> **AND IT IS NOW EXERCISED.** Since 2026-09-06 Eq. (A.4) at `Q_N ≠ 0` is
> `polrad_im_el_spin1`, `PolradFull` calls it, and T9 — un-skipped for this —
> gates the corrected `4/(1+η_A)` against a literal transcription in the test
> file at 1e-14. The correction is no longer dead code.
— POLRAD Eq. (A.4). Every `Q_N` appears linearly, so
`ℑ^el_i(Q_N) = U_i + Q_N T_i` with

```
  U_1 = (2/3) eta_A (1+eta_A) F_m^2                      = B(Q^2)/2
  U_2 = F_c^2 + (2/3) eta_A F_m^2 + (8/9) eta_A^2 F_q^2  = A(Q^2)
  U_5 = U_6 = U_7 = U_8 = 0
```

`U_1` and `U_2` are the **standard Rosenbluth `B/2` and `A`** of a spin-1 target
(`A = G_C² + (8/9)η²G_Q² + (2/3)ηG_M²`, `B = (4/3)η(1+η)G_M²`), which **fixes
the form-factor normalisation convention with no ambiguity** — see §2.1. The
spin-½ line, Eq. (A.5) `ℑ^el_1 = Z²η_A G_m²`,
`ℑ^el_2 = Z²(G_e² + η_A G_m²)/(1 + η_A)`, is the familiar `W₁^el, W₂^el` pair
with the same convention **and shows the `Z²` explicitly** — which is the second
independent pin on the `A`/`Z` normalisation (§5, T8b/T8c).

**Note `ℑ^el_1` and `ℑ^el_2` are NOT purely tensor** — only their `Q_N` pieces
are. Splitting them wrong silently doubles the unpolarised tail into the tensor
one; T9 catches it.

#### 1.4.7 Axis scope: the tail is emitted at EVERY `θ_S`


**The first draft zeroed `rc_tail` for `θ_S ≠ 0`. That was wrong and is
reversed.** Two independent arguments give the general axis with no new
machinery:

1. **POLRAD states the answer for the transverse case.** Eq. (43):
   `σ_q⊥^d = −½ σ_q∥^d`. The transverse quadrupolarised tail is exactly `−½`
   the longitudinal one — which is `P₂(cos 90°)/P₂(cos 0°)`. (POLRAD's basis
   also *does* carry `η_⊥`: Eq. (6) is `η = (η·p)p/M² − (η·η_L)η_L −
   (η·η_t)η_t − (η·η_⊥)η_⊥`, and the text under Eq. (5) says the fourth basis
   vector is added precisely for the normal-to-plane case. The earlier claim
   that "nothing exists" for it was a misreading.)
2. **Rotational invariance about the beam settles the general case.** The tail
   is *inclusive* — integrated over the scattered-lepton azimuth about `k₁` and
   over the real photon. Such a quantity can depend on a **rank-2** target
   polarisation only through the one rotational invariant available,
   `P₂(cos θ_S)`. So for any axis

```
     sigma^el(theta_S)  =  sigma^el_U  +  P_zz^eff * (1/6) sigma_q^A ,
     P_zz^eff = 3 Q_NN P2(cos theta_S)
```

   — the **same** `P_zz^eff` the band already uses (§1.1) and the same
   `Q_N ≡ P_zz^eff` identification of §1.4.4. Eq. (43) is the `θ_S = 90°`
   instance of it.

Therefore:

* **`rc_tail` is emitted at every `θ_S`**, with `Q_N → P_zz^eff`. `θ_S = 0`
  (`--plan tensor-thirds`, the default) is the special case `P_zz^eff = P_zz`.
* **The one caveat, printed by the run at `θ_S ≠ 0`:** the tail is the
  *φ-integrated* one, applied per event to a density that at `θ_S ≠ 0` carries
  `cos φ′` and `cos 2φ′` modulations. The tail therefore dilutes the
  φ-*averaged* rate correctly and the φ-*differential* rate only on average.
  Analyses binning in φ at a non-longitudinal axis must treat `rc_tail` as a
  bin-integrated correction, not a per-φ one.
* The **band** (`rc_tensor_lo/hi`) is likewise emitted for every `θ_S`: it only
  rescales whatever tensor part the library already computes.
* What still does **not** exist at any axis: RC for the *φ-dependent* tensor
  observables themselves (the `Δ` cos 2φ sector, the coherent channel). That
  exclusion is unchanged — see §1.5.3 and Q6.

### 1.5 Per-channel applicability

Written in `Channel` (`event.hpp`), which is what `rc.hpp` takes — see §3.3 for
the `PipelineChannel → Channel` mapping and why the API is keyed on `Channel`.

| `Channel` | `PipelineChannel` | `rc_tensor_lo/hi` | `rc_tail` | reason |
|---|---|---|---|---|
| `Inclusive` | `Inclusive` | **full** | **full** (every `θ_S`) | the home case; τ from `W_tensor/W` |
| `TaggedLi6Alpha` | `TaggedLi6Alpha` | **full** (τ from the cluster density, §1.5.1) | **≡ 1** | see §1.5.2 |
| `TaggedLi7Alpha` | `TaggedLi7Alpha` | **full** | **≡ 1** | see §1.5.2 |
| `TaggedDeuteronP` | `TaggedDeuteronP` | **full** | **≡ 1** | see §1.5.2 |
| `TaggedLi6D`, `TaggedLi7T`, `TaggedDeuteronN`, `TaggedHe3P` | *(unreachable from `PipelineChannel`)* | **full** | **≡ 1** | the rule is "every `Tagged*`: band only" — stated so the API is total over `Channel`, not just over the five the pipeline can build |
| `CoherentLi6` | `CoherentLi6` | **≡ 1** | **≡ 1** | see §1.5.3 |

#### 1.5.1 Tagged channels — where the tensor part actually lives

**This is the trap in the tagged channels and it must be handled explicitly.**
`StruckClusterOptions::inclusive_b1` defaults to **false**
(`include/lipolgen/pipeline.hpp:176-186`), so the struck cluster's DIS kernel
carries **no `b₁` at all** and `W_tensor ≡ 0`. Applying the inclusive recipe
would silently return `w_± = 1` and price nothing. The tagged tensor asymmetry
lives instead in the **spin ⊗ cluster correlation**, the `M`-dependence of the
spectator density `n_M(k, cos θ_k)` (`TaggedModel::n_of_kc`,
`azz_tensor_curve`).

The **general definition of §1.1 covers it with no new physics**. For a tagged
event drawn at ion projection `M`, struck-cluster projection `m_S` and
spectator `(k, c)`, the density factorises as

```
  rho_{M,m_S}(k,c; x,Q2,phi)  =  n_{M,m_S}(k,c) * W_DIS(m_S; x,Q2,phi)
```

and `τ_tag` is the rank-2 part of `ρ` in `M` at fixed `m_S`, divided by `ρ`.
**Two assumptions, both stated rather than hidden:** (i) `λ_e · P_e = 0` (§1.1)
— the sampler draws `(k, c)` from `|A_{m_S}(M)|²`
(`sample_kc_one`, `pipeline.cpp:500`), so on a *helicity* plan the vector term
`λ_e P_e (m_S/J) A_par` inside `W_DIS` does **not** cancel between `ρ_M` and its
`M`-average; (ii) `StruckClusterOptions::inclusive_b1 = false`, so `W_DIS` is
`m_S`-independent in its rank-2 sector. `tensor_thirds_plan` satisfies both
(`bookkeeping.cpp:109-117` sets `lam_e = 0`). With them, `W_DIS(m_S)` cancels
between `ρ_M` and `ρ̄` and

```
  tau_tag = 1 - nbar(k,c) / n_M(k,c) ,   nbar = (1/(2J+1)) SUM_M n_M(k,c)
```

— one `TaggedModel::n_of_kc` lookup per projection, 3 (or 4) nearest-cell
lookups per event, free.

**The algebra, done correctly this time.** For `J = 1`, `n_{+1} = n_{−1} ≡ n₁`
(the density is even in `M`: it is a parity eigenstate, so only *even*
rank-`P_L(cos θ_k)` moments survive), so `n̄ = (2n₁ + n₀)/3`, and with

```
  A_zz^wf = (n_{+1} + n_{-1} - 2 n_0)/(n_{+1} + n_{-1} + n_0) = 2(n_1 - n_0)/(2 n_1 + n_0)
```

the population itself is **exactly the §1.1 rate form**:

```
  n_M  =  nbar * [ 1 + (P_zz(M)/2) A_zz^wf ] ,   P_zz(+-1) = +1,  P_zz(0) = -2
     check M = +-1:  nbar (1 + (n_1-n_0)/(2n_1+n_0)) = ((2n_1+n_0)/3)(3n_1/(2n_1+n_0)) = n_1   OK
     check M = 0  :  nbar (1 - A_zz^wf)              = ((2n_1+n_0)/3)(3n_0/(2n_1+n_0)) = n_0   OK
```

Therefore

```
                    (P_zz/2) A_zz^wf
  tau_tag  =  ----------------------------
               1 + (P_zz/2) A_zz^wf

     tau_tag(M = +-1) = (n_1 - n_0)/(3 n_1)         [ P_zz = +1 ]
     tau_tag(M = 0)   = 2(n_0 - n_1)/(3 n_0)        [ P_zz = -2 ]
```

**This is the SAME expression as the inclusive τ of §1.1 — not `½ A_zz^wf`.**
The first draft's `τ_tag(M = ±1) = (n₁−n₀)/(2n₁+n₀) = ½A_zz^wf` dropped the
`1 + (P_zz/2)A_zz` denominator; with `|A_zz^wf|` up to **0.45–0.52** on the VMC
path (`CONVENTIONS.md`) the two differ by **20–35 %**, and shipping the wrong
one would have forced the implementer to make the tagged code inconsistent with
the inclusive code (or to "fix" the inclusive one). `A_zz^wf` here is exactly
`azz_tensor_curve` (`include/lipolgen/tagged.hpp:287-291`). **So `τ` has ONE
closed form on every channel**, with `A_zz` reading either the DIS asymmetry or
the wave-function asymmetry. (§5, T2 pins this.)

**Honest flag to carry in the header and the docs:** Gakh–Shekhovtsova and
POLRAD compute RC for *inclusive* tensor DIS. **No RC calculation exists for a
tagged tensor asymmetry.** Applying `δ(x)` to `τ_tag` is an extrapolation. It
is defensible — the lepton-vertex radiator is blind to the spectator, and
`A_zz^tag` is a ratio of counts in fills that share the same radiator, so the
cancellation is at least as good as inclusive — but it is **not** a published
result and must never be quoted as one. If anything, the tagged band is
*conservative*.

#### 1.5.2 Why the tagged `rc_tail` is exactly 1

`e + A → e' + γ + A(g.s.)` leaves the nucleus **intact**. The tag requires a
*fragment* in the Roman Pots. The two are cleanly separated — **not by
rigidity**: `⁶Li` and `α` have the *same* `Z/A = 0.5`, so `R = 1` for both, and
`⁷Li` (3/7) vs `α` (1/2) differ but the deuteron control's `d` (1/2) vs `p`
(1) does — **but by longitudinal fraction**: the elastic recoil sits at
`x_L = 1` and inside the 10σ beam-exclusion envelope (`PipelineConfig::n_sigma
= 10`), while the tagged fragment sits at `x_L ≈ A_spec/A_beam` = 2/3 (⁶Li α),
4/7 (⁷Li α), 1/2 (d control). The elastic tail is therefore **vetoed by the
tag itself**, `rc_tail ≡ 1.0` exactly, and the reason is a property of the
route classification, not an approximation.

**The quasi-elastic tail is a different matter and `rc_tail ≡ 1` is an
assumption there, not a fact.** A quasi-elastic event knocks a nucleon out and
**can** leave an intact α at `x_L ≈ 2/3` — exactly where the tag looks. Now that
the *unpolarised* QRT is priced inclusively (§1.4.3), its absence on the tagged
channels is the one place the tagged `rc_tail ≡ 1` line is unsupported.
~~It is probably small (a radiative tail of a quasi-elastic peak into a DIS
bin, with the α-tag's own acceptance on top) but nobody has computed it.~~
**Q3.**

> **WORKED OUT 2026-09-04 (task B4), and "probably small" is withdrawn.** The
> kinematics are sharper than "can leave an intact α", and they point the
> other way:
>
> 1. **The remnant has to break up.** Quasi-elastic knockout removes one
>    nucleon; the A−1 remnant of either lithium channel is **unbound** — ⁵Li
>    and ⁵He are resonances above the α + N threshold with no
>    particle-stable state — so it decays to α + N and the α emerges at
>    `x_L ≈ (4/5)(5/6) = 2/3`, the tag window itself, with a `p_T` spread of
>    the order of the Fermi motion the tag already accepts.
> 2. **In the α + d picture the generator actually uses, the α is a TRUE
>    SPECTATOR** whenever the struck nucleon is one of the embedded deuteron's
>    two: nothing touches it, and its momentum distribution is the same
>    `n_M(k, c)` the tagged Born is built on. On the deuteron control the
>    statement is bare — the quasi-elastic tail there **is** elastic `e`–`n`
>    scattering with a spectator proton, the classic spectator-tagging
>    background.
> 3. **So the tag does not suppress it in the RATIO.** `rc_tail` is a ratio,
>    the same spectator density and the same Roman-Pot acceptance multiply its
>    numerator and its tagged Born, and both select the **2 of 6** nucleons
>    inside the deuteron (a nucleon knocked out of the **α** destroys the α and
>    *is* vetoed — the other 4 of 6). To leading order the two factors cancel
>    and the omitted dilution is of the **same order as the inclusive
>    quasi-elastic one**, which is 22 % / 73 % / 99.9 % of the inclusive tail
>    at `x = 0.01 / 0.10 / 0.30`.
>
> That is an **order-of-magnitude argument, not a computed number**, and it is
> why the omission stays an omission: pricing it needs a **tagged** Born
> denominator (`RcModel::born_pb_at` is the inclusive one) and the tag
> acceptance folded into Eq. (44) — plus a fourth piece Eq. (44)'s free-nucleon
> sum does not contain at all, the **cluster-elastic**
> `e + A → e' + γ + d + α`, whose tag acceptance would be perfect. What
> changed is that `RcModel::exclusion_reason()` now says all of this in the run
> banner and in the npz `meta`, so `rc_tail ≡ 1` cannot be read as a veto on
> the whole tail. `../run_2026-09-03/phase_B_numbers.md` §B4.

#### 1.5.3 Why the coherent channel is excluded

Two independent reasons, both structural:

1. **The coherent channel's tensor dependence is entirely azimuthal.** It is
   the `1 + c₂ cos 2(φ_t − φ_S)` coefficient of `CoherentSampler::fill_event`
   (`include/lipolgen/coherent.hpp`, `a2_deformation` / `a2_m_state` /
   `cos2phi_coefficient`), a **φ-dependent tensor observable**. The literature
   verdict is unambiguous: *"What does not exist is any RC treatment for the
   φ-dependent (cos φ_TL, cos 2φ_TT) observables"* (`physics_literature.md`
   §3). There is nothing to interpolate and nothing to cite.
2. **The elastic peak is already outside the channel by construction.** The
   coherent channel writes the intact ground-state recoil *in the final state*
   and requires `M_X ≥ 1 GeV` (`CoherentXpomModel`, C1). The elastic point is
   `M_X = 0`. It is not a background to this measurement; it is the `x_P → 0`
   boundary the channel already excludes.

`RcModel::applies()` returns `false`, the weights are exactly 1.0, and the
`Pipeline` prints one line naming the exclusion (rather than throwing, so that
a scan over channels does not have to special-case `--rc`).

### 1.6 The spin-blind ISR/FSR shift — what it is, why it is NOT here, and what it would need

The lepton-vertex correction multiplies `F₁`, `F₂` and `b₁…b₄` by the *same*
`(x, Q²)` radiator; because `A_zz` is a ratio taken between spin fills at the
*same* beam energy and the *same* acceptance, the radiator cancels and only the
**`(x, Q²)` shape migration** survives. That is what both existing experiments
did (HERMES: RADGEN with the *unpolarised* photon spectrum, joint
detector+radiative unfolding; E12-13-011: no polarised RC at the lepton vertex
at all), and it is why the recommended treatment prices the residual as a band
instead of shifting momenta.

**If it were ever implemented as a shift**, it would need all of:

1. **A photon-energy sampler.** Either POLRAD's own infrared-free IRT
   `σ_Fin` (Eq. (13), the `R`, `τ` double integral over `ℑ_i(R, τ)`) or a
   leading-log equivalent-radiator / structure-function form (RADGEN's route).
   `physics_literature.md` §3 prices the full ℑ-based transcription at **+5
   person-days** over the leading-log radiator.
2. **One extra uniform per event.** This *breaks the "RNG untouched" invariant
   of this design* (§5, T6). It would therefore have to be its own option
   (`--isr`) drawing from its **own counter offset** in the
   `Rng(seed, run, bunch, event)` stream, so that turning it on does not move
   any other draw — the same discipline as `kCountStreamEvent`
   (`bookkeeping.hpp`).
3. **True vs observed kinematics.** The kernel would be evaluated at
   `(x_true, Q²_true)` before radiation and the record written at
   `(x_obs, Q²_obs)` from the radiated electron. The generator window is
   already looser than any analysis window (`generator_scenario`: Q² ≥ 0.7,
   y ∈ [0.004, 0.985], W² ≥ 8), which is exactly the migration headroom this
   needs — but the *sampler grid* `SamplerGrid{x_min = 1e-4, q2_min = 1}` would
   have to be widened as well.
4. **An explicit ISR photon in the record** (`Role::IsrPhoton`, pdg 22,
   status 1) so that `k + P_ion = k' + γ + p_spec + X` still closes exactly;
   the HepMC3 convention doc and the whole-nucleus balance tests would move.
5. **The multiplicative virtual factor** `exp(δ_vert + δ_vac^l + δ_vac^h)`
   (POLRAD Eqs. (12), (15), (16)) as an overall normalisation.
6. **A mutual exclusion with `--rc tensor-band`.** Applying both would double
   count: the band exists precisely because the shift is not applied.

---

## 2. Data inputs

### 2.1 ⁶Li elastic form factors — the only genuinely new physics input

POLRAD carries form factors for **d, ³He, C, O only**
(`physics_literature.md` §3, "Limits"). ⁶Li must be added.

**Normalisation convention, fixed by the Rosenbluth identity of §1.4.6** (not
by guesswork; the same identity is visible in Eq. (38)'s `σ_u^d`, whose
integrand is `A(Q²)X̃ − (2/3)(1+η_A)F_m²`):

```
  F_c(0) = Z = 3                                   [ POLRAD's deuteron has Z = 1, so its F_c(0) = 1 ]
  F_m(0) = (M_A / m_p) * mu_A / mu_N               [ = G_M(0) ]
  F_q(0) = M_A^2 * Q_A                             [ = G_Q(0), M_A in fm^-1, Q_A in fm^2 ]
```

**Validated against the deuteron** with the repository's own constants
(`beams.cpp` mass table, `HBARC_GEV_FM = 0.19733`):
`M_d = 1.875612942 GeV = 9.5049 fm⁻¹`, `Q_d = 0.2859 fm²`, `μ_d = 0.8574 μ_N`
⟹ `G_Q(0) = 25.829`, `G_M(0) = 1.7139` — the **textbook deuteron values**
(25.83, 1.714). The convention is therefore right, with no free choice.

**⁶Li numbers** (`M(⁶Li) = 5.601518702 GeV` from `beams.cpp`, i.e.
28.3866 fm⁻¹):

| quantity | value | provenance |
|---|---|---|
| `M_A` | 5.601518702 GeV | `src/core/beams.cpp` mass table — **reuse, do not retype** |
| `Z` | 3 | `LI6().Z` |
| `μ(⁶Li)` | **+0.8220473 μ_N** | measured; TUNL A = 6 evaluation (`refs/TUNL_A6_2002.pdf`) |
| `Q(⁶Li)` | **−0.0818(17) fm²** | measured; **TUNL A = 6 evaluation prints `Q = −0.818(17) mb (1998CE04)`** (`refs/TUNL_A6_2002.pdf`, the ⁶Li g.s. moment block). The often-quoted **−0.0806(6) fm²** is the *Pyykkö* compilation value (Mol. Phys. **106** (2008) 1965, from Cederberg et al. 1998), **not** TUNL's — the first draft printed −0.0806 with a TUNL citation, which was mis-sourced. The design uses the in-repo TUNL number; switching to Pyykkö's is a one-constant change and moves `F_q(0)` by 1.5 %. |
| `F_c(0)` | **3** | `= Z` |
| `F_m(0)` | **4.90765** | `(5.601518702/0.938272088)·0.8220473` |
| `F_q(0)` | **−65.914** | `(28.38655)²·(−0.0818)` (Pyykkö's −0.0806 would give −64.947) |

**How big is the ⁶Li tensor tail really? — the first draft's argument was
wrong.** It said *"`|F_q(0)|` for ⁶Li is 2.5× the deuteron's … so the ⁶Li tensor
elastic tail is not obviously smaller"*. `F_q` never appears alone. In
Eq. (A.4) it appears as `η_A F_q` or `η_A² F_q`, and

```
  eta_A F_q(0)  =  (t / 4 M_A^2) * M_A^2 Q_A  =  t Q_A / 4          <- M_A CANCELS
```

In the t-peak form Eq. (38) the same cancellation happens through `x_A`: the
leading `F_c F_q` interference carries `X̃X₁ ∝ x_A^{-2}` against `σ_u`'s
`F_c² X̃`, giving

```
  (Q_N/6) sigma_q / sigma_u   ~   (Q_N/6) * (Delta_eta) * F_q / F_c
                              =   (Q_N/6) * [t_cut/(4 M_A^2)] * M_A^2 Q_A / Z
                              =   (Q_N/24) * t_cut * Q_A / Z            <- M_A CANCELS AGAIN
```

**What actually matters is `Q_A / Z`, and `t_cut`, and both go the wrong way for
⁶Li:** `Q(⁶Li)/Q(d) = 0.0818/0.2859 = 0.29` (and opposite in sign), `Z = 3`
against 1, and the ⁶Li charge form factor dies at `t ≈ 0.25 GeV²` against the
deuteron's `t ∼ 1 GeV²`. So

> **the tensor fraction of the ⁶Li elastic tail is ≈ 10 % of the deuteron's**
> (`0.29/3 = 0.096`), of the opposite sign, and further suppressed by the
> smaller `t` reach.

The consequence for the whole design: **`rc_tail` is an unpolarised dilution
with a small tensor correction on top**, not "a tensor background". §8.1
measures `σ^el_T/σ^el_U` and settles the **magnitude**: `O(10⁻²)` per unit
`Q_N` is met at `x = 0.30` and missed by 20–60× at `x ≤ 0.1`, on **both** edges of
the C0 shape band added 2026-09-04, against a tail/Born ratio that is `O(10⁻¹)`.
It does **not** settle the **sign**. §8.1c item 1 withdraws the `+1.5595e−04`
published at `x = 0.10` as a *result*: the two edges read `+1.5595e−04` (`ho`) and
`−4.0702e−05` (`vmc-ft`) there, so at that `x` both the sign and the
"opposite in sign" half of the box above are **band edges**, not measurements.
They do hold on both edges at `x = 0.01` and `x = 0.30`, where the deuteron
reads `+0.106` and `−0.117` against ⁶Li's negative and positive.

> **Reviewer note / response.** The review argued from Eq. (A.4) that the ⁶Li
> tensor tail is *`F_m`-driven* (`ℑ₂`, `ℑ₇` carry `η F_m²`; `ℑ₈` carries
> `η F_m(F_m + 2F_q)`), and that the `F_c F_q` interference is `O(η²)` and
> negligible. That is true term-by-term in `ℑ`, but the `ℑ_i` are weighted by
> `θ_ij` and the Eq. (18) kernel, which supply `1/η` enhancements. In POLRAD's
> own *reduced* t-peak form — Eq. (38)'s `σ_q^d`, where those weights are
> already folded in — the `F_q(3F_c + …)` and `−2η X̃ F_q(4F_c − …)` terms carry
> `X̃X₁ ∝ 1/x_A²` and **dominate** the `F_m²` term by an order of magnitude at
> `x ≲ 0.1`. So the design keeps `F_c F_q` as the leading tensor term. **The
> review's substantive conclusions are adopted in full**: the `2.5×` argument is
> wrong, `M_A` cancels, the scale is `Q_A` (3.5× smaller) over `Z` (3× larger),
> and the ⁶Li tensor tail is ≈ 10 % of the deuteron's. The two `F_m`-specific
> fixes are also adopted below (its own shape; its own scale knob), because
> `F_m²` is the *sub*-leading tensor term and an `F_q`-only band would leave it
> unpriced either way.

**Shapes.** The recommended source in `physics_literature.md` §3 is
**R. B. Wiringa & R. Schiavilla,
[arXiv:nucl-th/9807037](https://arxiv.org/abs/nucl-th/9807037), PRL 81 (1998)
4317** (VMC AV18+UIX, IA and IA+MEC). It was read for this design. Two findings
change the recommendation:

* **It is figures only.** Fig. 1 plots `F_L(q)` and `F_T(q)` with the C0 and C2
  contributions shown separately (caption: *"the monopole (C0) and quadrupole
  (C2) contributions to the longitudinal form factor are also shown by the
  dashed…"*). **There is no table anywhere in the paper.** The same is true of
  the later ab-initio symmetry-adapted work
  ([arXiv:1502.03066](https://arxiv.org/abs/1502.03066)), which is C0 only.
* **Its quadrupole normalisation is wrong by 3×**, and the paper says so:
  *"Our prediction for the latter is −0.23(9) fm², larger (though with a 50 %
  statistical error) in absolute value than the measured value of −0.08 fm²"*.
  This is the same ⁶Li quadrupole failure already recorded in
  `OPEN_ITEMS_SOLUTIONS.md` §5 (`Q(⁶Li) = −0.33(18)` vs exp. −0.083, Pudliner
  et al.) with the standing instruction **"do not derive a ⁶Li tensor
  polarization or b₁ input from these wave functions"**. Taking the VMC
  normalisation would overstate the tensor tail by a factor ≈ 3.
  **`σ^el_T` is QUADRATIC in `F_q`, not linear** — Eq. (38)'s `σ_q^d` contains
  `F_q(3F_c + 3η F_m + η F_q)` and `F_q(4F_c − 3x F_m + (4/3)η F_q)`, and
  Eq. (A.4)'s `ℑ^el_2` and `ℑ^el_6` likewise contain `F_q · F_q` pieces (only
  `ℑ^el_8` is linear). So the ±100 % band **must be run, not rescaled**
  (§5, T12), and its two edges are not symmetric about the nominal.
* **Its `F_T` (i.e. `F_m`) shape is qualitatively different from a monopole.**
  WS98's text: *"The experimental `F_T(q)` is well reproduced by our
  calculations in the **first peak at q = 0.5 fm⁻¹**, but the zero comes a
  little too early and the **second peak at q = 2 fm⁻¹** is somewhat
  overpredicted."* A single harmonic-oscillator monopole shape with its first
  zero near 3 fm⁻¹ (the v0 `F_point` below) has **no** structure in
  `0.5–2 fm⁻¹` — which is precisely the `t` range
  (`q = 2 fm⁻¹ ⟺ t ≈ 0.155 GeV²`) where the tail lives. Giving `F_m` the shared
  `F_point` shape is therefore the largest *shape* error in v0, and it is a
  tensor-sector error (`F_m²` is the sub-leading tensor term and `F_m` also
  enters `σ_q^d`'s `3η F_m` and `−3x F_m` interferences).

**Therefore, the design's data plan:**

**(P) Primary, v0 — analytic, no data file, everything in `rc.hpp`.**
Two-parameter harmonic-oscillator (HO) point-nucleon density
`ρ(r) ∝ (1 + α r²/a²) e^{−r²/a²}`, whose form factor is exact:

```
  F_point(q) = [ 1 - (alpha/(2+3 alpha)) * (q^2 a^2 / 2) ] * exp(- q^2 a^2 / 4)
  <r^2>_point = (3/2) a^2 (2 + 5 alpha) / (2 + 3 alpha)
  first zero at  q_0 = (1/a) sqrt( 2(2+3 alpha)/alpha )
```

with the *normalisations* taken from the **measured** moments above:

```
  F_c(t) = 3        * F_point(q) * [ G_E^p(t) + G_E^n(t) ]     (nucleon folding, N = Z)
  F_q(t) = -65.914  * F_point(q) * [ G_E^p(t) + G_E^n(t) ]     * fq_scale
  F_m(t) = 4.90765  * F_mag(q)                                 * tail_tensor_scale
                                                  q^2 = t / (hbar c)^2
```

**Two departures from the first draft, both forced by the paragraphs above:**

* **`F_m` gets its own shape `F_mag`, not the shared `F_point`.** v0 uses a
  two-parameter form with a zero between the WS98 peaks,
  `F_mag(q) = (1 − q²/q_z²) exp(−q² b²/4)`, fitted to reproduce a first peak at
  `q ≈ 0.5 fm⁻¹`, a zero at `q_z ≈ 1.2–1.4 fm⁻¹` and a second peak at
  `q ≈ 2 fm⁻¹` — i.e. WS98 Fig. 1's `F_T`, digitised if the fallback table
  (F) lands, and fitted to those three landmarks otherwise. `q_z` and `b` are
  `HoSpin1FFOptions::fm_qz_fm`, `fm_b_fm`.
* **The nucleon folding uses `G_E^p + G_E^n`, not `G_E^p` alone.** `⟨r²⟩_point`
  below is derived by subtracting *both* `⟨r²⟩_p` and `(N/Z)⟨r²⟩_n`, so folding
  with the proton form factor alone is internally inconsistent. With `N = Z = 3`
  the isoscalar combination is the right one, and `G_E^n(0) = 0` keeps
  `F_c(0) = Z` exact.

**Starting parameters, derived here, to be refit by the implementer:**
`a = 1.9069 fm`, `α = 0.13822`. They are the HO pair satisfying
(i) `⟨r²⟩_point = 6.078 fm²` — from `r_ch(⁶Li) = 2.589(39) fm`
(Angeli & Marinova, ADNDT 99 (2013) 69) minus `⟨r²⟩_p = 0.7071` plus
`−⟨r²⟩_n = 0.1155` minus the Darwin–Foldy `0.033` fm² — and (ii) a first C0
zero at `q₀ = 3.1 fm⁻¹`.

**`q₀ = 3.1 fm⁻¹` is a STARTING GUESS, not a source.** WS98 does not give the
zero's position: it says only that C2 *"is much smaller than C0 below
3 fm⁻¹"*, that *"for q ≥ 3 fm⁻¹ the C2 contribution becomes dominant, and the
shoulder seen in the data is entirely due to this component"*, and — separately
— that two-body currents *"shift the minimum to lower values of q"*, with no
number. So the implementer must **fit `q₀` to the Li/Sick/Whitney/Yearian and
Suelzle data over a stated range, `q₀ ∈ [2.9, 3.3] fm⁻¹`**, and T11 pins the
**refit result** with its own quoted uncertainty — not `3.10 ± 0.02`, which the
first draft attributed to WS98 and which WS98 does not say. **Note the
shell-model value `α = (Z−2)/3 = 1/3` puts the zero at 2.33 fm⁻¹, far too low —
so `α` must be free; this is a phenomenological fit, not a shell model, and must
be labelled as such.**

**(F) Fallback / mandatory upgrade — a digitised table under `data/`.**
`data/ff/li6_elastic.csv`, columns `q_fm^-1, F_C0, F_C2, F_M1, source`, read
by `TabulatedSpin1FF` with log-linear interpolation and a hard extrapolation
refusal. Sources, in order of preference:

1. **WS98 Fig. 1** (IA+MEC curves) for the **shapes** of C0, C2 and F_T —
   digitised. **Do NOT rescale the whole C2 curve by ≈ 0.35.** The first draft
   proposed exactly that, and it destroys the one high-`q` constraint the next
   paragraph relies on: WS98 says `F_L` is *"in excellent agreement with
   experiment"* and that the `q ≥ 3 fm⁻¹` shoulder *"is entirely due to"* C2 —
   i.e. the VMC C2 **magnitude is right at high `q`** while its `q → 0` limit
   (`Q = −0.23(9)`) is 3× too big. A global 0.35 makes the high-`q` C2 wrong by
   3× instead. **Constrain the two ends separately** (next paragraph) and let
   the *interpolation* between them carry the model uncertainty.
2. **The classical elastic data** for the C0 normalisation and low-`q` slope:
   **L. R. Suelzle, M. R. Yearian, H. Crannell, Phys. Rev. 162 (1967) 992**
   (`|F(q)|` for ⁶Li and ⁷Li) and **G. C. Li, I. Sick, R. R. Whitney,
   M. R. Yearian, Nucl. Phys. A162 (1971) 583** (high-`q` ⁶Li). These are the
   "Li/Sick/Suelzle-type data" of the brief. They measure `|F_L|² = F_C0² +
   F_C2²`, so they constrain the C0 below `q ≈ 2.5 fm⁻¹` (where C2 is
   negligible) and only the *sum* above it.

**What the quadrupole one needs, stated precisely.** No experiment has
separated `F_C2(q)` for ⁶Li: the elastic longitudinal data measure
`F_C0² + F_C2²` and the C2 is sub-dominant exactly where the data are good. So
the C2 has **two** independently checkable constraints and nothing in between:

* **low `q`** — `F_C2 → (2√2/3) η_A G_Q(0)` with `G_Q(0) = M_A² Q(⁶Li)` fixed
  by the **measured** quadrupole moment. Exact, no model.
* **high `q` (≳ 3 fm⁻¹)** — the shoulder in `|F_L|`, which WS98 states is
  *entirely* C2. This is a genuine data constraint on the C2 magnitude there,
  extractable from the Li et al. data by subtracting a C0 model.
* **in between** — model. Build `F_C2(q)` as an interpolation *pinned at both
  ends*: the exact low-`q` limit above, the WS98 shoulder magnitude at
  `q ≳ 3 fm⁻¹` (or Li et al. minus a C0 model), and a smooth shape between them
  (WS98's, renormalised **locally** so both ends are honoured — e.g. a
  `q`-dependent factor running from `0.35` at `q → 0` to `1` above `3 fm⁻¹`).
  Carry a **±100 % band on `F_q`** as the model uncertainty **of that
  interpolation only**, because (a) the only ab-initio calculation gets `Q`
  wrong by 3× and (b) `physics_literature.md` §5 already forbids deriving a ⁶Li
  tensor input from these wave functions. Expose it as `RcOptions::fq_scale`
  (default 1.0; run 0.0 and 2.0). **`σ^el_T` is QUADRATIC in `F_q`, so the band
  must be RUN, not rescaled** (§5, T12) — the first draft's "costs nothing to
  evaluate" was wrong.
* **the magnetic one** — `F_m` carries the sub-leading tensor term (`F_m²`) and
  two interferences, and its WS98 shape has a zero and a second peak inside the
  `t` range that matters. Expose **`RcOptions::tail_tensor_scale`** (default
  1.0; run 0.5 and 2.0) as a flat multiplier on `F_m`, so the `η F_m²` sector
  is banded too. `fq_scale` alone does not span the tensor tail.

The conversion between the multipole and Sachs conventions, needed to read any
of these sources into POLRAD's `F_q`:

```
  |F_L(q)|^2 = F_C0^2 + F_C2^2 ,   F_C0 = F_c ,   F_C2 = (2 sqrt(2)/3) eta_A F_q
```

which reproduces `A(Q²) = G_C² + (8/9)η²G_Q² + (2/3)ηG_M²` term by term.

### 2.2 Inputs reused, not redefined

| input | where it already lives |
|---|---|
| `M(⁶Li)`, `Z`, `A`, `spin` | `src/core/beams.cpp` mass table, `beams.hpp` `LI6()` |
| `α_EM`, `M_NUCLEON`, `ħc`, `GeV⁻²→pb`, `TENSOR_LL_SIGN` | `include/lipolgen/constants.hpp` |
| `F₁, F₂, R, b₁, b₂, Δ` and the Born `σ_U` | `sf.hpp`, `xsec.hpp` (`dsigma_unpol`) |
| the rank-2 geometry `Q_NN`, `P₂(cos θ_S)` | `InclusiveKernel::tensor_moments` |
| the accepted `(x, Q²)` cell grid and `Event::kin.cell` | `sampler.hpp`, `event.hpp` |
| the spectator density `n_M(k, c)` | `TaggedModel::n_of_kc` |

**New numbers, and they exist in exactly one place — `rc.hpp`:**
`RC_DELTA_HIGH_X`, `RC_X_HIGH`, `RC_DELTA_LOW_X`, `RC_DELTA_LOW_X_OPTIMISTIC`,
`RC_X_LOW`, and the ⁶Li elastic form-factor block
(`LI6_MU_N`, `LI6_QUADRUPOLE_FM2`, `LI6_FF_HO_A_FM`, `LI6_FF_HO_ALPHA`,
`LI6_FF_FM_QZ_FM`, `LI6_FF_FM_B_FM`). These are RC-owned, not library-wide
(`docs/CONVENTIONS.md`: constants live in the *owning* header).

**One number DOES move into `constants.hpp`: `M_ELECTRON = 0.51099895e-3` GeV.**
It is currently a bare literal in `src/hepmc/hepmc_writer.cpp:136`, and
`RcOptions` needs it (the lepton mass is what makes `F_IR` finite and what sets
`l_m = ln(Q²/m²)`). Defining it a second time in `rc.hpp` would violate
`CONVENTIONS.md` "No physics number is hard-coded in two places", so the design
adds it to `constants.hpp` and **changes `hepmc_writer.cpp` to use it** in the
same commit.

**Three numbers the first draft retyped and this one does not:**
`HoSpin1FFOptions::m_a_gev = 5.601518702` (duplicates the `beams.cpp` mass
table — its *own comment* said "pass it, never retype"), `z = 3.0` (duplicates
`LI6().Z`) and `RcOptions::m_lepton` (above). §3.1 replaces the first two with
`0` sentinels filled by an `HoSpin1FF::for_ion(const Ion&, …)` factory. That
also fixes a correctness bug: `RcModel` takes an arbitrary `Ion` (⁷Li and d are
reachable) while `HoSpin1FFOptions` was ⁶Li-specific **by default**, so a ⁷Li
run would have silently used ⁶Li form factors.

---

## 3. API

### 3.1 `include/lipolgen/rc.hpp` — header sketch

> **DATED NOTE, 2026-09-06 (task B2).** The sketch below still reads
> `PolradFull = 1, ... NOT IMPLEMENTED` and *"v0 ships `TPeak` only"*. **That
> is no longer true.** `RcTailModel::PolradFull` is **implemented and shipped**
> (opt-in, `--rc-tail-model polrad-full`); the enum gained a third value
> `TPeakPlusLL` on 2026-09-03; and the refusal reproduced in §3.2 below
> (*"tail_model PolradFull is not implemented"*) was replaced by the `n_eta >=
> 64` and `m_lepton` rules. Read the sketch as the v0 specification it was, and
> `include/lipolgen/rc.hpp` for what ships. §1.4.6 above and
> `../run_2026-09-06/phase_B_numbers.md` §B2 carry the measurements.

```cpp
#ifndef LIPOLGEN_RC_HPP
#define LIPOLGEN_RC_HPP

/// \file rc.hpp
/// Tensor-sector radiative corrections as an OPT-IN, WEIGHT-ONLY family.
///
/// ---------------------------------------------------------------------------
/// RC IS A WEIGHT, AND IT IS NOT EVEN ON THE NOMINAL ONE.
///
/// `fsi.hpp` multiplies `Event::weight` because FSI is a CORRECTION to the
/// model.  These are not: `rc_tensor_lo/hi` is a SYSTEMATIC VARIATION and
/// `rc_tail` is a BACKGROUND, and both must leave the Born sample alone.  They
/// travel in `Event::rc_weights` under their own HepMC3 names and an analysis
/// multiplies them in on purpose.  No four-vector moves, no random number is
/// consumed, and `RcMode::Off` (the default) is today bit for bit.
/// ...
#include <memory>
#include <string>
#include <vector>

#include "lipolgen/constants.hpp"
#include "lipolgen/event.hpp"
#include "lipolgen/sampler.hpp"
#include "lipolgen/tagged.hpp"

namespace lipolgen {

// ------------------------------------------------------------- the band

/// PR12-13-011 (E12-13-011 PROPOSAL, unpublished): "the unpolarized
/// corrections are known to better than 1.5 %".  x_high is that experiment's
/// OWN lower x edge, 0.16 (arXiv:2506.04506 p. 8: 0.16 < x < 0.49,
/// 0.8 < Q^2 < 5.0 GeV^2).  Quoting 1.5 % below x ~ 0.1 is unsupported and
/// contradicts HERMES's measured 15 % residual at x = 0.063.
inline constexpr double RC_DELTA_HIGH_X = 0.015;
inline constexpr double RC_X_HIGH       = 0.16;
/// Gakh-Shekhovtsova (hep-ph/0403262) SIZE of the RC on the spin-dependent
/// cross section at x ~ 1e-3 - 1e-2: "changes from 10% to 30% as compared
/// with the Born contribution".  ZERO INSPIRE CITATIONS -- a single unchecked
/// calculation.  THIS IS A CORRECTION MAGNITUDE, NOT A MEASURED RESIDUAL:
/// taking it as a 1-sigma band is a CHOICE this generator makes because it
/// does not apply the correction.  Default = the conservative 0.30 edge.
inline constexpr double RC_DELTA_LOW_X            = 0.30;
/// The residual HERMES actually ACHIEVED at the lowest x: 2e-3 against a
/// MEASURED |A_zz| of 1.06e-2 at <x> = 0.012 (hep-ex/0506018 Table II) =>
/// 0.19 fractional (0.19, 0.19, 0.15 in the three low-x bins).  NOT 0.10 --
/// that came from dividing 2e-3 by the RANGE MAXIMUM |A_zz| <= 0.02, which
/// occurs at x = 0.45 where HERMES says RC is negligible.
inline constexpr double RC_DELTA_LOW_X_OPTIMISTIC = 0.19;
inline constexpr double RC_X_LOW                  = 0.01;

/// Log-linear interpolation between the two anchors, clamped outside them;
/// monotone non-increasing in x, exact at both anchors.
double rc_delta(double x,
                double delta_high = RC_DELTA_HIGH_X,
                double delta_low  = RC_DELTA_LOW_X,
                double x_high     = RC_X_HIGH,
                double x_low      = RC_X_LOW);

// --------------------------------------------- spin-1 elastic form factors

/// POLRAD Eq. (A.4)'s (F_c, F_m, F_q) at the elastic-vertex t [GeV^2].
/// Normalisation is the Rosenbluth one, fixed by Eq. (A.4)'s Q_N = 0 limit:
///   F_c(0) = Z,  F_m(0) = (M_A/m_p) mu_A/mu_N,  F_q(0) = M_A^2 Q_A .
class Spin1ElasticFF {
 public:
  virtual ~Spin1ElasticFF() = default;
  virtual double fc(double t_gev2) const = 0;
  virtual double fm(double t_gev2) const = 0;
  virtual double fq(double t_gev2) const = 0;
  virtual std::string provenance() const = 0;   ///< printed by the run
};

/// The v0 6Li model: one two-parameter harmonic-oscillator shape for all
/// three multipoles, normalised on the MEASURED moments (NOT on VMC -- see
/// the header block and docs/open_items/run_2026-09-02/design_C_tensor_rc.md
/// section 2.1: WS98's VMC Q(6Li) = -0.23(9) fm^2 is 3x the measured -0.0806).
/// NO ION-SPECIFIC DEFAULTS.  `m_a_gev = 0` and `z = 0` mean "take them from
/// the Ion"; `HoSpin1FF::for_ion` fills them.  Retyping M(6Li) or Z here
/// would duplicate the beams.cpp mass table and `LI6().Z`
/// (CONVENTIONS.md: no physics number in two places), and an ion-specific
/// DEFAULT would silently give a 7Li or deuteron run 6Li form factors.
struct HoSpin1FFOptions {
  double a_fm     = 0.0;   ///< 0 => from the ion's (LI6_FF_HO_*) block
  double alpha    = 0.0;
  double z        = 0.0;   ///< 0 => ion.Z
  double m_a_gev  = 0.0;   ///< 0 => ion.mass()  (the AME mass, NOT A*M_NUCLEON)
  double mu_n     = 0.0;   ///< 0 => the ion's measured moment
  double q_fm2    = 0.0;   ///< 0 => the ion's measured quadrupole moment
  /// F_m's OWN shape (WS98 F_T: first peak 0.5, zero, second peak 2 fm^-1).
  /// F_mag(q) = (1 - q^2/q_z^2) exp(-q^2 b^2/4).  A shared monopole shape has
  /// no structure there, which is where the tail lives.
  double fm_qz_fm = 0.0;
  double fm_b_fm  = 0.0;
  /// The +-100 % quadrupole band.  sigma^el_T is QUADRATIC in F_q, so this
  /// must be RUN (0, 1, 2), never rescaled from one run.  See T12.
  double fq_scale = 1.0;
  /// Flat multiplier on F_m -- the eta F_m^2 tensor sector, which fq_scale
  /// does NOT span.  Run 0.5, 1, 2.
  double tail_tensor_scale = 1.0;
  bool   fold_nucleon = true;   ///< multiply by (G_E^p(t) + G_E^n(t)), N = Z
};

class HoSpin1FF : public Spin1ElasticFF {
 public:
  /// The ONLY supported constructor path.  Reads `ion.mass()` and `ion.Z`,
  /// fills every zero field from the ion's measured-moment block, and THROWS
  /// for any ion that has no such block (today: everything but 6Li).
  static std::shared_ptr<HoSpin1FF> for_ion(const Ion& ion,
                                            HoSpin1FFOptions opt = {});
  ...
};

/// A digitised (q, F_C0, F_C2, F_M1) table under data/ff/, log-linear in q,
/// refusing to extrapolate.  `fq_scale` and `tail_tensor_scale` apply here too.
///
/// The path is resolved under `data_dir()` (cluster.hpp -- $LIPOLGEN_DATA_DIR,
/// else the compiled-in prefix), exactly like the VMC tables.  CONVENTIONS.md:
/// "Nothing resolves a data path against the working directory."
class TabulatedSpin1FF : public Spin1ElasticFF {
 public:
  /// `relative` is relative to `data_dir()`, e.g. "ff/li6_elastic.csv".
  static std::shared_ptr<TabulatedSpin1FF> from_data_dir(
      const std::string& relative, HoSpin1FFOptions norm);
  ...
};

// ------------------------------------------------------------- the model

enum class RcMode : int {
  Off        = 0,   ///< today, bit for bit; every weight identically 1.0
  TensorBand = 1,   ///< the band + the radiative tails
};

/// Which tail formulation.  v0 ships `TPeak` only; `PolradFull` is the
/// documented upgrade path (design section 1.4.6).
enum class RcTailModel : int {
  TPeak      = 0,   ///< DEFAULT: POLRAD Eqs. (37)-(39), (43) -- one eta_A
                    ///< integral.  POLRAD sec. 2.1.3 B: for a TAIL the t-peak
                    ///< is the leading contribution and the s-/p-peaks are
                    ///< SUPPRESSED, so this is not an approximation to the
                    ///< collinear limit -- it is the dominant piece.
  PolradFull = 1,   ///< Eq. (18) + Appendix B + Eq. (A.4).  NOT IMPLEMENTED.
                    ///< [2026-09-06: IMPLEMENTED -- see the note above this
                    ///< block.  The shipped enum keeps this value (TPeak = 0,
                    ///< PolradFull = 1) and adds TPeakPlusLL = 2.]
};
const char* rc_mode_name(RcMode m);          ///< "off", "tensor-band"

/// Which part of the rank-2 sector the band rescales.
enum class RcScope : int {
  TensorRate = 0,   ///< DEFAULT: the b1..b4 (T_LL) sector only -- what POLRAD
                    ///< and Gakh-Shekhovtsova actually compute
  TensorAll  = 1,   ///< also the Delta (gluon-transversity) cos 2phi term.
                    ///< NO LITERATURE SUPPORT: nobody has computed RC for a
                    ///< phi-dependent tensor observable.  For pricing only.
};

struct RcOptions {
  /// NOTE: there is NO `mode` here.  `PipelineConfig::rc` is the single
  /// source of truth and `RcModel` takes the mode as a constructor argument.
  /// The first draft had both and let `rc` silently overwrite
  /// `rc_options.mode`, so `cfg.rc_options.mode = TensorBand` alone did
  /// nothing.  (The FSI precedent likewise keeps one config enum.)
  RcScope scope = RcScope::TensorRate;

  // --- the band ---------------------------------------------------------
  double delta_high_x = RC_DELTA_HIGH_X;
  double x_high       = RC_X_HIGH;
  double delta_low_x  = RC_DELTA_LOW_X;   ///< 0.30 default; 0.10 = "as good as HERMES"
  double x_low        = RC_X_LOW;

  // --- the radiative tails ----------------------------------------------
  bool with_tail    = true;   ///< false = band only (a diagnostic run)
  bool with_qe_tail = true;   ///< the UNPOLARISED quasi-elastic tail (Eq. 44).
                              ///< Its A_zz is ~0, so it DILUTES A_zz just as
                              ///< the unpolarised elastic tail does; HERMES
                              ///< subtracted both.  Turning it off prices the
                              ///< elastic tail alone and MUST be labelled so.
  RcTailModel tail_model = RcTailModel::TPeak;
  std::shared_ptr<const Spin1ElasticFF> ff;   ///< null => HoSpin1FF::for_ion
  /// The +-100 % quadrupole band, applied by RcModel when it builds the
  /// DEFAULT form factor.  Lives here (not only on HoSpin1FFOptions) so that
  /// `make_config(rc_fq_scale=...)` and `--rc-fq-scale` have something to
  /// set.  IGNORED, with a printed line, when `ff` is user-supplied.
  /// QUADRATIC in the tail: run it, do not rescale (T12).
  double fq_scale = 1.0;
  /// Flat multiplier on F_m -- the eta F_m^2 tensor sector.  Same rules.
  double tail_tensor_scale = 1.0;
  /// Flat multiplier on the whole quasi-elastic tail, standing in for
  /// POLRAD Eq. (44)'s S_E/S_M/S_EM suppression factors, which v0 sets to 1
  /// (the conservative direction for a dilution).  Run 0.0 / 0.5 / 1.0.
  double qe_suppression = 1.0;
  /// eta_A quadrature of Eq. (38): Gauss-Legendre nodes in ln(eta_A).
  int    n_eta       = 128;
  double m_lepton    = M_ELECTRON;     ///< constants.hpp -- NOT a second literal
  /// Ceiling on the returned tail ratio, mirroring GlauberFsiOptions::w_max:
  /// a Monte-Carlo weight must be bounded, and the clipped fraction is
  /// reported rather than hidden -- globally AND per y-band, because
  /// Y_+ ~ 1/(1-y) makes the y -> 1 edge the only place it bites.
  double tail_max    = 10.0;
};

/// One event's RC weight triple.
struct RcWeights {
  double lo   = 1.0;   ///< "rc_tensor_lo"  = 1 - delta(x) * tau
  double hi   = 1.0;   ///< "rc_tensor_hi"  = 1 + delta(x) * tau
  double tail = 1.0;   ///< "rc_tail"       = 1 + sigma_tail / sigma_Born
};

inline constexpr std::size_t kRcWeightCount = 3;
/// "rc_tensor_lo" | "rc_tensor_hi" | "rc_tail" for i = 0, 1, 2.
const char* rc_weight_name(std::size_t i);
/// The HepMC3 name of entry `i` in SLOT `slot` (0 = the event's own category,
/// 1 + k = spin category k, matching the existing "spin_weight_k"):
/// slot 0 -> "rc_tail", slot 3 -> "rc_tail_3".
std::string rc_weight_name(std::size_t i, std::size_t slot);

/// The tensor-sector RC weight family.  IMMUTABLE after construction and safe
/// to share between threads, exactly like `GlauberFsiWeight`: the constructor
/// does all the work (one tau_A quadrature per accepted (x, Q2) cell).
class RcModel {
 public:
  /// Inclusive / coherent form.  `dis` is the sampler the events were drawn
  /// from; the tail tables are built over ITS accepted (x, Q2) grid and read
  /// by interpolation at the event's own (x, Q2) (design section 1.4.5).
  ///
  /// THROWS when any category of `plan` has `lam_e * pe != 0`.  The whole
  /// A_zz programme -- and every tau identity in this header -- assumes an
  /// UNPOLARISED beam (design section 1.1); with a polarised one the vector
  /// sector would be priced as if it were rank 2.
  ///
  /// The constructor resolves, ONCE, a `const StateTables*` for every
  /// (category k, projection m) exactly the way `CategoryPlan::states` does,
  /// and stores the pointers.  It must NOT call
  /// `InclusiveSampler::state_tables()` per event: that takes `cache_mutex_`
  /// on EVERY call (sampler.cpp:271-279) and would serialise
  /// `Pipeline::for_each(sink, 4)`.
  RcModel(RcMode mode, RcOptions opt,
          std::shared_ptr<const InclusiveSampler> dis,
          RunPlan plan, Channel channel, Ion ion);
  /// Tagged form: adds the cluster model, whose n_M(k, c) carries the whole
  /// tensor structure when `StruckClusterOptions::inclusive_b1` is false
  /// (the default) -- see the design doc section 1.5.1.  On this path the
  /// constructor ALSO rebuilds, per struck-cluster projection m_S, the PURE
  /// `SpinCategory(j = s_channel, e_{m_S}, lam_e, pe, theta_s = 0,
  /// phi_s = 0)` key that `InclusiveKinematicsSource::plan_for` uses
  /// (pipeline.hpp:158, private), because that is the key the struck-cluster
  /// StateTables are cached under.
  RcModel(RcMode mode, RcOptions opt,
          std::shared_ptr<const InclusiveSampler> dis,
          RunPlan plan, Channel channel, Ion ion,
          std::shared_ptr<const TaggedModel> tagged);

  const RcOptions& options() const;
  RcMode mode() const;
  /// The form factor actually in use (the user's `RcOptions::ff` or the
  /// model-built default), so the run banner and Python can print its
  /// provenance -- `RcOptions::ff` is null on the default path and the built
  /// object was otherwise unreachable.
  const Spin1ElasticFF& ff() const;
  std::string ff_provenance() const;   ///< == ff().provenance()

  /// The per-`Channel` rule, TOTAL over the enum (event.hpp has 9 values,
  /// PipelineChannel only 5; see design section 1.5):
  ///   Inclusive   -> band AND tail
  ///   every Tagged* (incl. TaggedLi6D/TaggedLi7T/TaggedDeuteronN/TaggedHe3P,
  ///                  which no PipelineChannel builds today) -> band, tail == 1
  ///   CoherentLi6 -> neither; both weights exactly 1.0
  bool applies() const;
  /// False whenever the tail is identically 1: any tagged channel (the
  /// elastic recoil sits at x_L = 1, inside the 10-sigma beam envelope, and
  /// is vetoed by the tag), `with_tail = false`, or `!applies()`.  NOTE it is
  /// TRUE at every theta_S -- the tail is emitted at any axis through
  /// P_zz^eff = 3 Q_NN P2(cos theta_S) (POLRAD Eq. (43); design 1.4.7).
  /// The run PRINTS the reason whenever this is false.
  bool tail_applies() const;
  const std::string& exclusion_reason() const;

  // --- the pieces, exposed so a run and a test can see them --------------
  double delta(double x) const;                        ///< rc_delta at the run's knobs

  /// The rank-2 fraction tau of section 1.1 = W_tensor/W, for the event's OWN
  /// spin state.  READS EXACTLY THESE FIELDS, and nothing else, so a test can
  /// build a record by hand (T1, T2, T14):
  ///   ev.channel
  ///   ev.kin.cell, ev.kin.phi, ev.kin.x, ev.kin.q2,
  ///   ev.kin.k, ev.kin.cos_theta_k            (tagged only)
  ///   ev.spin.j, ev.spin.m_ion, ev.spin.m_struck,
  ///   ev.spin.lam_e, ev.spin.pe, ev.spin.theta_s, ev.spin.phi_s
  double tensor_fraction(const Event& ev) const;
  /// ... for spin category `k` of the run plan (weighted mode).  This is the
  /// category's POPULATION MIXTURE, sum_m p_m W_m, mirroring
  /// `InclusiveSampler::weights_for` (sampler.cpp:461-483) -- NOT a pure
  /// state.  It equals `tensor_fraction(ev)` only when category k's
  /// population vector is pure.
  double tensor_fraction(const Event& ev, std::size_t k) const;

  /// sigma_tail / sigma_Born at accepted cell `c` and tensor degree `q_n`.
  double tail_ratio(int cell, double q_n) const;
  /// The same at an ARBITRARY (x, Q2), by interpolation of the tables --
  /// public so a test is not hostage to which cells the sampler accepted.
  /// Throws outside the table's (x, Q2) support.
  double tail_ratio_at(double x, double q2, double q_n) const;
  /// The three precomputed tables, for plotting and for the T8 gates.
  const std::vector<double>& sigma_tail_u() const;   ///< elastic, unpolarised
  const std::vector<double>& sigma_tail_t() const;   ///< elastic, tensor (sigma_q^A)
  const std::vector<double>& sigma_tail_qe() const;  ///< quasi-elastic, unpolarised
  /// Fraction of table nodes that hit `RcOptions::tail_max`, globally and per
  /// y-band (y < 0.5, 0.5-0.9, > 0.9).  LOG ALL FOUR: the tail grows like
  /// Y_+ ~ 1/(1-y), so a clipped edge hides inside a small global number.
  double clipped_cell_fraction() const;
  std::array<double, 3> clipped_fraction_by_y() const;

  // --- what the pipeline calls ------------------------------------------
  RcWeights weights(const Event& ev) const;                    ///< slot 0
  RcWeights weights(const Event& ev, std::size_t k) const;     ///< slot 1+k
  /// Fill `ev.rc_weights` with the whole row-major (3 x n_slot) block,
  /// n_slot = 1 + ev.spin_weights.size().  A no-op when `mode == Off`.
  ///
  /// TAKES NO `Rng&`.  That signature is the STRUCTURAL guarantee that
  /// `--rc tensor-band` cannot move the random stream; T6 is the runtime one.
  void fill(Event& ev) const;

 private:
  // POLRAD Eq. (38): three one-dimensional eta_A quadratures at one (x, Q2)
  // -- sigma_u^A, sigma_q^A (the Q_N/6 partner) and the Eq. (44) QRT.
  void build_tail_tables();
  struct TailTriple { double u, t, qe; };
  TailTriple ert_at(double x, double q2) const;
  // Resolved ONCE in the constructor: [category k][projection index m].
  std::vector<std::vector<const InclusiveSampler::StateTables*>> states_;
  ...
};

/// `PipelineConfig::rc`.  Alias of `RcMode` so the CLI, the config struct and
/// the model cannot disagree (the FSI precedent keeps two enums; this one
/// deliberately does not).
using PipelineRc = RcMode;
const char* pipeline_rc_name(PipelineRc r);

}  // namespace lipolgen
#endif
```

### 3.2 Changes to existing headers

**`include/lipolgen/event.hpp`** — one field on `Event`, plus one line in
`reset()`:

```cpp
  /// Opt-in radiative-correction weights (rc.hpp), row-major
  /// (kRcWeightCount x n_slot) with n_slot = 1 + spin_weights.size():
  /// slot 0 is the event's OWN, PURE spin state m; slot 1+k is spin category
  /// k's population MIXTURE, sum_m p_m W_m, mirroring
  /// `InclusiveSampler::weights_for`.  The two agree only when category k's
  /// population vector is pure -- on the tensor-thirds plan (P_zz = 0.6
  /// mixtures) they do not.  P_zz differs by category, so the band
  /// coefficient does too.
  /// EMPTY when `PipelineConfig::rc` is Off -- and then the HepMC3 weight
  /// vector and the npz are bit-for-bit what they are today.
  std::vector<double> rc_weights;
```
`Event::reset()` gains `rc_weights.clear();` beside `spin_weights.clear();`
(`src/core/event.cpp:26`) — capacity kept, so the allocator still sees nothing
per event.

**`include/lipolgen/xsec.hpp`** — one method, and `amplitudes()` refactored to
call it so the b-sector geometry exists **once**:

```cpp
  /// The b-sector (tensor RATE) contribution to (w_avg, a1, a2) ALONE -- the
  /// piece the tensor RC band rescales (rc.hpp).  Exactly the terms
  /// `amplitudes()` adds inside its `j >= 1` branch, on whichever path
  /// `tensor_gamma` selects, with the Delta (gluon-transversity) cos 2phi
  /// term included only when `with_delta` (RcScope::TensorAll).
  Amplitudes tensor_amplitudes(const SFTables& t, double x, double q2,
                               double s, const EventSpinState& state,
                               bool with_delta = false) const;
```
`amplitudes()` then reads `out += tensor_amplitudes(...)` and adds the vector
terms. **No number moves** — it is a pure extraction, and `test_xsec.cpp` /
`test_tensor_gamma.cpp` / the `validation/reference/*.json` rtol-1e-12 gates
prove it.

**`include/lipolgen/sampler.hpp`** — three vectors on `StateTables` and one
method, the exact mirror of `weights_for`:

```cpp
  struct StateTables {
    std::vector<double> w_avg, a1, a2;
    /// The b-sector part of the same three, for the RC band (rc.hpp).
    /// Filled by the SAME `build_state` call; nothing else reads them.
    std::vector<double> w_tensor, a1_tensor, a2_tensor;
    ...
  };

  /// The TENSOR-only counterpart of `weights_for`: t[i, k] is the rank-2 part
  /// of the same density W_k(event_i).  The RC band is
  /// w_+- = 1 +- delta(x) * t[i,k] / w[i,k]  (rc.hpp).
  std::vector<double> tensor_weights_for(
      const EventBatch& batch, const std::vector<SpinCategory>& cats) const;
```

**`include/lipolgen/pipeline.hpp`** — the config field, mirroring `fsi`:

```cpp
  /// Tensor-sector radiative corrections as OPT-IN, WEIGHT-ONLY families
  /// (rc.hpp): the band `rc_tensor_lo`/`rc_tensor_hi` on the tensor part of
  /// the rate, and the 6Li tensor ELASTIC radiative tail `rc_tail`.  They go
  /// on `Event::rc_weights` and NOT on `Event::weight` -- a systematic
  /// variation and a background, not a correction -- so `Off` (the default)
  /// is today bit for bit, in every four-vector, in `Event::weight` and in
  /// the RNG stream.  Applies on every channel except `CoherentLi6`; the
  /// elastic tail additionally applies only on `Inclusive` at theta_S = 0
  /// (design_C_tensor_rc.md sections 1.5, 1.4.6).
  PipelineRc rc = PipelineRc::Off;
  RcOptions  rc_options;          ///< knobs ONLY -- there is no mode on it
```

plus, on `Pipeline`, the accessor and the member (mirroring `fsi_weight()` /
`fsi_`):

```cpp
  /// The RC weight model of the run; null when `PipelineConfig::rc` is Off.
  const RcModel* rc_model() const { return rc_.get(); }
  ...
  std::shared_ptr<RcModel> rc_;   ///< null when cfg_.rc == Off
```

**`PipelineConfig::validate()`** (`src/core/pipeline.cpp:456`) gains, after the
FSI block:

```cpp
  if (rc != PipelineRc::Off) {
    if (!(rc_options.delta_low_x >= 0.0 && rc_options.delta_high_x >= 0.0)) {
      throw std::runtime_error("PipelineConfig: rc delta must be >= 0");
    }
    if (!(rc_options.x_low > 0.0 && rc_options.x_high > rc_options.x_low)) {
      throw std::runtime_error(
          "PipelineConfig: rc needs 0 < x_low < x_high (the band anchors)");
    }
    if (rc_options.n_eta < 8) {
      throw std::runtime_error("PipelineConfig: rc n_eta must be >= 8");
    }
    if (!(rc_options.fq_scale >= 0.0 && rc_options.tail_tensor_scale >= 0.0 &&
          rc_options.qe_suppression >= 0.0)) {
      throw std::runtime_error(
          "PipelineConfig: rc fq_scale / tail_tensor_scale / qe_suppression "
          "must be >= 0");
    }
    if (rc_options.tail_model != RcTailModel::TPeak) {          // v0 only
      throw std::runtime_error(
          "PipelineConfig: rc tail_model PolradFull is not implemented "
          "(design_C_tensor_rc.md section 1.4.6)");
    }
    // [2026-09-06: GONE.  Every tail model runs; what validate() refuses now
    //  is n_eta < 64 under PolradFull and a non-default m_lepton under the
    //  two t-peak models.  See the dated note in section 3.1.]
    for (const SpinCategory& c : plan.categories()) {
      if (c.lam_e != 0 && c.pe != 0.0) {
        throw std::runtime_error(
            "PipelineConfig: --rc assumes an UNPOLARISED beam (lam_e*pe == 0); "
            "the tensor fraction tau is the rank-2 projection and a polarised "
            "beam would put the vector sector in it "
            "(design_C_tensor_rc.md section 1.1)");
      }
    }
  }
```
It does **not** refuse the coherent channel: `RcModel::applies()` returns
false, the weights are 1.0 and the run prints the reason, so a channel scan
does not have to special-case `--rc`.

### 3.3 The pipeline hook — exact location

`src/core/pipeline.cpp`, the **dispatcher** `Pipeline::event(index, out)` at
lines 1038–1053 (not inside `make_tagged`, because RC applies on every
channel), immediately after the channel switch and **before** the hadronizer:

```cpp
void Pipeline::event(std::uint64_t index, Event& out) const {
  const std::size_t k = category_of(index);
  const std::uint64_t local = index - offset_[k];
  Rng rng(cfg_.seed, cfg_.run, k, local);
  switch (cfg_.channel) {
    ...
  }
  // rc.hpp: the tensor-sector RC weight family.  A PURE FUNCTION of the
  // FINISHED record and of tables built in the constructor -- it consumes no
  // random number, moves no four-vector and does not touch `Event::weight`,
  // which is exactly what makes `--rc tensor-band` bit-for-bit identical to
  // `--rc off` in every one of those.  Null when `PipelineConfig::rc` is Off,
  // and then `Event::rc_weights` stays empty (`Event::reset` cleared it).
  // It sits AFTER the switch because RC applies on every channel, and BEFORE
  // the hadronizer so that a T2 event carries the same weights a T0 one does.
  if (rc_) rc_->fill(out);

  if (cfg_.hadronizer) cfg_.hadronizer(out, rng);
}
```

and, in the constructor, **after the whole channel `if`/`else` chain** — i.e.
immediately before the `// --- beams ---` block, where `dis_sampler_`,
`model_`, `plan_` and `cplan_` are resolved on **every** channel:

```cpp
  // --- rc (rc.hpp) --------------------------------------------------------
  // ANCHOR: this must sit AFTER the channel if/else chain (tagged /
  // Inclusive / CoherentLi6) and BEFORE `// --- beams ---`.  It cannot go
  // beside the FSI block: that block is INSIDE `if (is_tagged(cfg_.channel))`,
  // and on the inclusive and coherent channels `dis_sampler_` is not assigned
  // until the non-tagged else-branch and `cplan_` not until the Inclusive
  // branch -- building RcModel there would dereference a null `dis_sampler_`
  // on the channel this design calls the home case.
  if (cfg_.rc != PipelineRc::Off) {
    rc_ = std::make_shared<RcModel>(cfg_.rc, cfg_.rc_options, dis_sampler_,
                                    plan_, event_channel_of(cfg_.channel),
                                    beams_.ion,
                                    model_ /* null off the tagged channels */);
  }
```

**The location is given by that anchor comment, not by a line number** — the
first draft said "beside the FSI block (`pipeline.cpp:469-480`)", which is
inside the tagged branch (`pipeline.cpp:452-481`) and would have crashed the
inclusive channel.

**`event_channel_of` does not exist yet and this design adds it**, file-local in
`pipeline.cpp` (the first draft called a `channel_of` that exists nowhere; the
`PipelineChannel → Channel` map is only implicit inside the generators, e.g.
`gc.channel = Channel::Inclusive` at `pipeline.cpp:564`):

```cpp
// pipeline.cpp, file-local.  `rc.hpp` is keyed on `Channel` (event.hpp) and
// CANNOT include pipeline.hpp -- pipeline.hpp includes rc.hpp for PipelineRc
// and RcOptions.  So the 5 -> 9 widening happens here, once.
static Channel event_channel_of(PipelineChannel c) {
  switch (c) {
    case PipelineChannel::Inclusive:      return Channel::Inclusive;
    case PipelineChannel::TaggedLi6Alpha: return Channel::TaggedLi6Alpha;
    case PipelineChannel::TaggedLi7Alpha: return Channel::TaggedLi7Alpha;
    case PipelineChannel::TaggedDeuteronP:return Channel::TaggedDeuteronP;
    case PipelineChannel::CoherentLi6:    return Channel::CoherentLi6;
  }
  throw std::runtime_error("event_channel_of: unhandled PipelineChannel");
}
```

`rc.hpp`'s `applies()`/`tail_applies()` are documented over **all nine**
`Channel` values (§1.5's table), so the four that no `PipelineChannel` reaches
still have a defined answer and T14 is written over that mapping.

`dis_sampler_` is the ION-level sampler on the inclusive and coherent channels
and the **struck-cluster** sampler on the tagged ones — which is exactly the
sampler `Event::kin.cell` indexes on each channel (`event.hpp`
`Kinematics::cell` docstring), so the cell lookup is correct on all five with
no branch.

---

## 4. Output & bindings

### 4.1 HepMC3 weight names

`src/hepmc/hepmc_writer.cpp::Impl::ensure_run_info` (line 64; the name list at line 72) currently builds
`{"nominal", "spin_weight_1", …}`. **Append** the RC block, so no existing
index moves:

```cpp
    std::vector<std::string> names = {"nominal"};
    for (std::size_t i = 0; i < ev.spin_weights.size(); ++i) {
      names.push_back("spin_weight_" + std::to_string(i + 1));
    }
    // rc.hpp: appended AFTER the spin block so that "nominal" stays index 0
    // and every existing spin_weight_k keeps its index.  Emitted only when
    // the run has RC on -- an --rc off file is byte-identical to today's.
    const std::size_t n_slot = ev.rc_weights.size() / kRcWeightCount;
    for (std::size_t s = 0; s < n_slot; ++s) {
      for (std::size_t i = 0; i < kRcWeightCount; ++i) {
        names.push_back(rc_weight_name(i, s));
      }
    }
```

and in `write()` (line 114 of the same file), extend the per-event vector the same way:

```cpp
  genevt.weights().assign(1 + ev.spin_weights.size() + ev.rc_weights.size(), 1.0);
  genevt.weights()[0] = ev.weight;
  for (...) genevt.weights()[i + 1] = ev.spin_weights[i];
  const std::size_t off = 1 + ev.spin_weights.size();
  for (std::size_t i = 0; i < ev.rc_weights.size(); ++i)
    genevt.weights()[off + i] = ev.rc_weights[i];
```

Resulting names, unweighted mode: `nominal`, `rc_tensor_lo`, `rc_tensor_hi`,
`rc_tail`. Weighted mode with N categories: `nominal`, `spin_weight_1..N`,
`rc_tensor_lo`, `rc_tensor_hi`, `rc_tail`, `rc_tensor_lo_1`, `rc_tensor_hi_1`,
`rc_tail_1`, … `rc_tail_N`.

### 4.2 npz columns

`python/bindings.cpp` `struct Columns` (line 173) gains
`std::vector<double> rc_tensor_lo, rc_tensor_hi, rc_tail;` filled in
`fill_row` (line 201) from `ev.rc_weights[0..2]` (slot 0), and `columns_to_dict`
(line 328) emits the three keys **only when the run has RC on**
(`p.rc_model() != nullptr`), so an `--rc off` npz has exactly today's key set.

**`export.columns_from_events` must be changed too** (`export.py:70-108`, the
`Event`-record path). It builds its own dict and does not go through
`columns_to_dict`, so without a change a records-path dict from an RC-on run
would carry no `rc_*` columns while `rc_columns()` silently returned ones. Add
the three keys there whenever `ev.rc_weights` is non-empty (and leave them out
entirely when it is empty, matching the C++ path).

`python/lipolgen/export.py` gains

```python
#: rc.hpp weight columns; present only when the run had --rc on.
RC_KEYS = ("rc_tensor_lo", "rc_tensor_hi", "rc_tail")

def rc_columns(columns):
    """(lo, hi, tail) arrays, exactly 1.0 where the run had --rc off."""
```

and `INCLUSIVE_KEYS` / `TAGGED_KEYS` are **not** changed (they are polligen
schemas). The per-category weighted-mode block, when it exists, is exported as
a **`(n, 3 * n_slot)`** array under the key `rc_weights` via the existing
`move_array2` (`bindings.cpp:102`), with `rc_weight_names` (a list of str) in
the npz `meta` giving the column order. *Not* `(n, 3, n_slot)`: only
`move_array` and `move_array2` exist, and adding a `move_array3` for one
consumer is not worth it — reshape in Python from `rc_weight_names`.

### 4.3 CLI

`python/lipolgen/cli.py`:

```python
    p.add_argument("--rc", choices=sorted(RC), default=None,
                   help="tensor-sector radiative corrections as OPT-IN "
                        "WEIGHTS (rc.hpp): 'off' (default) is today bit for "
                        "bit; 'tensor-band' adds rc_tensor_lo/rc_tensor_hi "
                        "(the band on the tensor part of the rate) and "
                        "rc_tail (the 6Li tensor elastic radiative tail). "
                        "Never a momentum shift; never on Event.weight.")
    p.add_argument("--rc-delta-low-x", type=float, default=None,
                   help="the low-x band edge (default 0.30, the conservative "
                        "end of Gakh-Shekhovtsova's UNCITED 10-30 %%; 0.10 is "
                        "the residual HERMES actually achieved). BAND IT: "
                        "run both, never quote one row alone.")
    p.add_argument("--rc-delta-high-x", type=float, default=None,
                   help="the high-x band edge (default 0.015, E12-13-011)")
    p.add_argument("--rc-fq-scale", type=float, default=None,
                   help="+-100 %% systematic on the 6Li quadrupole form "
                        "factor.  The tensor tail is QUADRATIC in it, so the "
                        "band must be RUN (0.0, 1.0, 2.0), never rescaled "
                        "from one run")
    p.add_argument("--rc-tail-tensor-scale", type=float, default=None,
                   help="multiplier on the 6Li MAGNETIC form factor -- the "
                        "eta*F_m^2 tensor sector, which --rc-fq-scale does "
                        "NOT span (default 1.0; run 0.5 and 2.0)")
    p.add_argument("--rc-qe-suppression", type=float, default=None,
                   help="multiplier on the UNPOLARISED quasi-elastic "
                        "radiative tail, standing in for POLRAD Eq. (44)'s "
                        "S_E/S_M/S_EM factors (default 1.0 = no suppression, "
                        "the conservative direction; run 0.0 and 0.5)")
```

`DEFAULTS` gains `rc="off", rc_delta_low_x=0.30, rc_delta_high_x=0.015,
rc_fq_scale=1.0, rc_tail_tensor_scale=1.0, rc_qe_suppression=1.0`;
`make_config(...)` gains the same keyword set (each writing the matching
`RcOptions` field); and the run banner gains, beside the FSI line
(`cli.py:236`):

```python
    if p.rc_model is not None:
        r = p.rc_model
        say("  RC %s: delta(x) = %.4g at x=0.01, %.4g at x=0.1; tail %s%s"
            % (_l.rc_mode_name(cfg.rc), r.delta(0.01), r.delta(0.1),
               "on" if r.tail_applies else "OFF",
               "" if r.tail_applies else " (%s)" % r.exclusion_reason))
        say("     6Li FF: %s;  clipped %.3g%% (by y-band %.3g/%.3g/%.3g%%)"
            % ((r.ff_provenance, 100.0 * r.clipped_cell_fraction)
               + tuple(100.0 * f for f in r.clipped_fraction_by_y)))
        say("     tails: elastic + UNPOLARISED quasi-elastic (Eq. 44, "
            "S=%.2g); the POLARISED QE tail is NOT priced (Zhou et al.)"
            % cfg.rc_options.qe_suppression)
        say("     Gakh-Shekhovtsova hep-ph/0403262 has ZERO INSPIRE "
            "citations: band it, never quote one row alone.")
```

`python/lipolgen/__init__.py` gains, next to `FSI`:

```python
#: RC weight families accepted on the command line.  "off" (the default) is
#: today bit for bit; "tensor-band" emits rc_tensor_lo / rc_tensor_hi /
#: rc_tail on Event.rc_weights and on NOTHING else -- never Event.weight,
#: never a four-vector, never a random number.
RC = {"off": _lipolgen.RcMode.Off, "tensor-band": _lipolgen.RcMode.TensorBand}
```
and `__all__` gains `"RC"`.

### 4.4 pybind11

`python/bindings.cpp`, beside the FSI block (lines 2009–2126):

```cpp
  py::enum_<RcMode>(m, "RcMode", "PipelineConfig.rc: the RC weight family.");
  py::enum_<RcScope>(m, "RcScope", "Which rank-2 terms the band rescales.");
  py::class_<RcOptions>(m, "RcOptions")            // every field def_readwrite
  py::class_<RcWeights>(m, "RcWeights")            // lo, hi, tail
  py::class_<Spin1ElasticFF, PySpin1ElasticFF, std::shared_ptr<...>>(...)
  py::class_<HoSpin1FFOptions>(m, "HoSpin1FFOptions")
  py::class_<HoSpin1FF, Spin1ElasticFF, std::shared_ptr<HoSpin1FF>>(...)
  py::class_<TabulatedSpin1FF, ...>(...).def_static("from_data_dir", ...)
  py::class_<RcModel, std::shared_ptr<RcModel>>(m, "RcModel")
      .def("delta", &RcModel::delta, py::arg("x"))
      .def("tensor_fraction", ...)
      .def("tail_ratio", &RcModel::tail_ratio, py::arg("cell"), py::arg("q_n"))
      .def("tail_ratio_at", &RcModel::tail_ratio_at,
           py::arg("x"), py::arg("q2"), py::arg("q_n"))
      .def("weights", ...)
      .def_property_readonly("ff_provenance", &RcModel::ff_provenance)
      .def_property_readonly("clipped_fraction_by_y", ...)
      .def_property_readonly("applies", &RcModel::applies)
      .def_property_readonly("tail_applies", &RcModel::tail_applies)
      .def_property_readonly("exclusion_reason", &RcModel::exclusion_reason)
      .def_property_readonly("clipped_cell_fraction", ...)
      .def_property_readonly("sigma_tail_u", ...)   // numpy
      .def_property_readonly("sigma_tail_t", ...);
  m.def("rc_delta", &rc_delta, ...);
  m.def("rc_mode_name", &rc_mode_name);
  m.def("rc_weight_name", (std::string(*)(std::size_t, std::size_t)) &rc_weight_name);
```
plus `Event`: `.def_readonly("rc_weights", &Event::rc_weights)`;
`PipelineConfig`: `.def_readwrite("rc", ...)`, `.def_readwrite("rc_options", ...)`;
`Pipeline`: `.def_property_readonly("rc_model", ...)` returning `None` when null,
mirroring `fsi_weight`.

### 4.5 Docs

* `docs/USAGE.md` — a "Radiative corrections" section: what the three weights
  are, that they are **not** on `Event::weight`, how to apply them
  (`w = weight * rc_tensor_hi`), the per-channel table of §1.5, the band
  provenance with the uncited-calculation flag, and the mandatory
  `--rc-delta-low-x 0.19 / 0.30`, `--rc-fq-scale 0 / 1 / 2`,
  `--rc-tail-tensor-scale 0.5 / 1 / 2` and `--rc-qe-suppression 0 / 1` band
  runs.
* `docs/OPEN_ITEMS_SOLUTIONS.md` — row 9 → `implemented 2026-09-0x`, plus a
  `§9` block with the measured numbers of §8.
* `docs/PHYSICS_CHANNELS.md` **§11 "Not handled (explicitly)"** — gains the
  RC bullets: the *polarised* quasi-elastic tail is not priced; there is no RC
  for the φ-dependent tensor observables (`Δ`, the coherent channel) at any
  axis; the tagged band is an uncited extrapolation; `RcTailModel::PolradFull`
  (implemented 2026-09-06, and STILL not checked against Mo-Tsai)
  is unimplemented.
  > **Reviewer note / response.** The review said `docs/PHYSICS_CHANNELS.md`
  > "does not exist". It does — 144 KB, with §11 *"Not handled (explicitly)"*
  > at line 397, which is the right home. What is **not** in it is any row
  > saying *"RC (until Phase C)"*; that string appears only in this design, so
  > the first draft's instruction to *remove* such a row is dropped and
  > replaced by the additions above.
* **`docs/HEPMC3_CONVENTION.md`** — it owns the weight-name rule that §4.1
  changes (lines 194-198). It must gain the `rc_*` block, the append-never-
  insert rule and the "empty when `--rc off`" statement, or the convention doc
  and the writer disagree.
* `docs/CONVENTIONS.md` — a bullet under "Physics defaults that are a CHOICE":
  the band anchors; that `δ_low = 0.30` is the **size of an unapplied
  correction taken as a 1σ band**, not a measured residual; that
  `RC_DELTA_LOW_X_OPTIMISTIC = 0.19` is HERMES's *measured* fractional residual
  at its lowest-x bin; and that the ⁶Li quadrupole form-factor normalisation is
  the **measured** moment and not VMC.
* `docs/USAGE.md` — additionally: that `rc_tail` carries the elastic **and
  unpolarised quasi-elastic** tails, and that the polarised QE tail is not
  priced.

---

## 5. Tests

C++ in a new `tests/test_rc.cpp` (picked up by the `GLOB` at
`CMakeLists.txt:144`); Python in `python/tests/test_rc.py`.

| # | test | expected value / limit |
|---|---|---|
| **T1** | **The band identity.** For the ⁶Li inclusive kernel, `tensor_thirds_plan(0, 0.6)`, `θ_S = 0`, a pure state `m`, at a grid of `(x, Q²)`: `RcModel::weights(ev).hi` equals `[1 + (P_zz/2)A_zz(1+δ)] / [1 + (P_zz/2)A_zz]` built from **`asymmetries::azz(b1, f1, f2, x, y)` at its default `θ_m = 0`** (the *longitudinal* asymmetry — all the axis geometry is in `P_zz^eff`, so passing `θ_m = θ_S` would apply `P₂(cos θ_S)` twice) and `P_zz = 3m²−2`. | equality to **1e-13** relative; this is the anchor that ties the code's `τ` to the published `A_zz` |
| **T2** | **The tagged identity.** Same, `TaggedLi6Alpha`, `inclusive_b1 = false`, **`λ_e = 0`** (as `tensor_thirds_plan` already sets, `bookkeeping.cpp:115-117` — with a helicity plan the vector term `λ_e P_e (m_S/J) A_par` does not cancel between `ρ_M` and `ρ̄` and this identity is false): `τ_tag = (P_zz/2)A_zz^wf / [1 + (P_zz/2)A_zz^wf]` against `azz_tensor_curve` at the event's `(k, c)`, i.e. `(n₁−n₀)/(3n₁)` at `M = ±1` and `2(n₀−n₁)/(3n₀)` at `M = 0`. **Not** `½A_zz^wf` / `−A_zz^wf`. The `(k, c)` must be chosen **on the `TaggedModel` grid nodes** — `n_of_kc` is a nearest-cell lookup (`tagged.hpp:131-133`), so an off-node point compares two different cells. | **1e-12** relative |
| **T3** | **Band monotone in x, anchors exact.** `rc_delta` on a 2000-point log grid `x ∈ [1e-5, 0.9]`. | non-increasing everywhere; the four **clamp-branch** identities exactly (bit-for-bit `==`, legitimate because the branch returns the constant unmodified): `rc_delta(0.16) == 0.015`, `rc_delta(0.9) == 0.015`, `rc_delta(0.01) == 0.30`, `rc_delta(1e-5) == 0.30`; and the §1.3 δ-table reproduced to **1e-12 against the full-precision doubles printed there** (a 5-significant-figure table cannot meet 1e-12 — that was a self-contradiction in the first draft) |
| **T4** | **Signs.** (a) `w_lo + w_hi == 2.0` exactly, every event. (b) Flipping the state `m = ±1 → 0` flips the sign of `τ` and swaps the roles: `w_hi(m=0) − 1` and `w_lo(m=±1) − 1` have the same sign, and `[w_hi(m=0)−1]/[w_hi(m=±1)−1] = P_zz(0)/P_zz(±1) · [W(±1)/W(0)] = −2·W(±1)/W(0)`. (c) **The literature form, so that recompiling with `TENSOR_LL_SIGN = +1` FAILS it:** at `θ_S = 0`, `sign(w_hi − 1) == −sign(b1 · P_zz)` (i.e. `b1 > 0` ⟹ the `m = 0` state has the larger rate). The first draft asserted `sign(w_hi−1) == sign(TENSOR_LL_SIGN · b1 · P_zz)`, which holds for *either* value of the constant and therefore guards nothing — the opposite of the `test_xsec.cpp` precedent it invoked, which is written against the literature relation without naming the constant. | (a) exact; (b), (c) 1e-12 |
| **T5** | **The tail's `y` and `t_min` behaviour.** Replaces the first draft's "strictly decreasing in `Q²`, faster than `1/Q²`" claim, which is **not implied by the physics**: the t-peak's `t_min = (x M_N)²` is independent of `Q²`, and `Y₊ = [1 + (1−y)²]/(1−y)` **grows** toward `y → 1`, so at fixed `x` the tail need not fall at all. (The first draft's grid also left the phase space — `Q² = 50` at `x = 0.01` needs `y = 1.26` at config 1's `s ≈ 3980 GeV²` — and evaluated through the private `ert_at`.) **(i) `y`-scan:** at fixed `Q² = 5 GeV²` and **stated `s`**, `w_tail − 1` over `y ∈ {0.1, 0.3, 0.5, 0.7, 0.9, 0.97}` tracks `Y₊(y)` to within a slowly varying factor: `[w_tail−1](y)·(1−y)/[1+(1−y)²]` varies by less than 3× across the scan and is **monotone increasing in `y`** over the last three points. **(ii) `t_min` scan:** at fixed `y`, `w_tail − 1` falls with `x` as `t_min = (x M_N)²` climbs into the form factor. **(iii)** every grid point is inside the sampler's scenario: `Q² ≤ y_max · s · x` with `y_max = 0.985`, computed per `x`, and `s` printed in the failure message. Evaluated through the **public** `RcModel::tail_ratio_at(x, q2, q_n)`. | (i), (ii) as stated; no point outside `y ≤ 0.985` |
| **T6** | **RNG untouched, four-vectors untouched, `Event::weight` untouched.** Generate 5000 events at `--rc off` and at `--rc tensor-band` with the same `(seed, run)`; compare `Event::number`, every particle's `(pdg, status, role, mass, charge, mother1, mother2)` and every four-vector component with `==` (not `Approx`), and `Event::weight`, `xsec_pb`, `kin.*`, `spin.*` likewise. **The stream check uses a `HadronizerHook`** (`pipeline.hpp:300`), not a "counting wrapper": no such wrapper exists — `Rng` (`rng.hpp`) exposes only `uniform()/normal()/next_u64()` and `Pipeline::event` constructs its own `Rng` internally, so the stream position is not observable from outside. Install `[&](Event&, Rng& r){ log.push_back(r.uniform()); }` on **both** pipelines and compare the two `log`s with `==`. Comparing T0 four-vectors (the `test_fsi.cpp:332` precedent) does **not** prove this on its own, because RC runs *after* every T0 draw. The structural guarantee is `RcModel::fill(Event&) const` taking no `Rng&`; this is the runtime one. | **bit-for-bit identical**, zero extra uniforms consumed |
| **T7** | **Born normalisation against POLRAD Eq. (9).** `x s · dsigma_unpol(x, q2, s)` equals `(4πα²S/Q⁴)[x y² F₁ + (1−y)F₂]` built directly from `SFTables`, at `Q_N = 0`, `P_L = 0`. | 1e-12 relative; this is the gate that keeps the tail's denominator from being a second definition of the Born |
| **T8** | **Tail normalisation — POLRAD-INTERNAL gates, not a radiator.** The first draft asked the leading-log equivalent-radiator limit to reproduce Eq. (18) to 5 %. It cannot: POLRAD §2.1.3 B says the **s- and p-peaks are suppressed** for a tail and only the **t-peak** leads, and the t-peak reaches `t_min = (x M_N)²` — inside the form factor — while the s-peak sits at `t ≈ (1−y)Q²` where `F_c(⁶Li)` is dead. Worse, that test would have *determined* the prefactor it was meant to check (the implementer tunes until 5 % appears), which is exactly how a missing `A = 6` gets absorbed. **Three internal gates instead, each with an unambiguous reference:** **(a)** `RcModel`'s `σ_u` quadrature at `Q_N = 0`, run with **deuteron** inputs (`Ion` = d, `A = 1` nuclear map degenerate to `M_A = M_d`), must reproduce the design's own literal transcription of Eq. (38) `σ_u^d` (a 15-line function written *in the test file*, not in `rc.cpp`) to **1e-10**; likewise `σ_q` at `Q_N ≠ 0` against Eq. (38) `σ_q^d`. **(b)** the **carbon** line: with a spin-0 ion (`P_N = Q_N = 0`) and a single form factor `F`, the model must reproduce `σ_u^C = (α³/S) Z² Y₋ ∫dη_A/η_A X̃ F²` — the **explicit `Z²`** is what pins the charge normalisation, and the `A`-normalisation follows from the §1.4.5 per-nucleon reduction. **(c)** the `A = 1`, spin-½ limit (POLRAD Eq. (A.5) / Eqs. (38) `σ_u^p`, (40)) against the **Mo–Tsai exact proton elastic tail** at one fixed-target point (`E = 10 GeV`, `x = 0.3`, `Q² = 2 GeV²`), to **10 %** — the residual there is genuinely the ultrarelativistic approximation, and it is a *cross-check*, not the definition. | (a), (b) 1e-10; (c) 10 % |
| **T8'** | **The QRT.** `σ^q_U` at `S_e = S_m = S_em = 1` must equal `Z·σ_u^p[G_E^p, G_M^p] + N·σ_u^p[G_E^n, G_M^n]` built in the test from the same Eq. (38) transcription, in **nucleon** invariants. And `qe_suppression = 0` must reproduce the elastic-only tail bit for bit. | 1e-12 |
| **T9** | **Rosenbluth limit.** *v0 form:* Eq. (38)'s `σ_u^d` integrand is `A(Q²)X̃ − (2/3)(1+η_A)F_m²` with `A = F_c² + (2/3)η F_m² + (8/9)η² F_q²` — assert the code's `σ_u` integrand against that, and that `σ_q`'s integrand is identically 0 when `F_m = F_q = 0`. *`PolradFull` form (`SKIP`-ped from 2026-09-02 to 2026-09-06, **UN-SKIPPED** when the model landed):* with `Q_N = 0`, `ℑ^el_2 == A(Q²)`, `ℑ^el_1 == B(Q²)/2 = (2/3)η(1+η)F_m²`, `ℑ^el_{5,6,7,8} == 0` identically — asserted through the SHIPPED `polrad_im_el_spin1`, at 1e-15 against `rosenbluth_spin1` and with `== 0.0` (not "small") on the four tensor entries, plus the six `Q_N` parts against a literal Eq. (A.4) at 1e-14. It was the suite's ONE unconditional skip. | exact; 1e-15 / 1e-14 |
| **T19** | **`RcTailModel::PolradFull` itself** (2026-09-06, and it is where the design's own "upgrade path" is discharged). **(a)** THE LIMIT, stated: `x_A → 0` with a form factor **dead at the s-/p-peak vertex**, so only the t-peak survives and Eq. (38)'s ultrarelativistic extraction is exact — measured `Eq. (18)/Eq. (38)` = 1.00230 (unpol.) / 1.00313 (tensor) at `x_A = 0.003`, tending to 1, on three form-factor sectors run separately (the full deuteron `σ_q^d` is a ~7× cancellation between two of them, so *its* ratio is not a gate). **(b)** the τ_A range and the peak positions inside it. **(c)** the transcription check's own three deuteron points. **(d)** the leading-log fallback at the HERMES point. **(e)** convergence in `n_eta`. **(f)** THE GAP IT DOES NOT CLOSE: the quasi-elastic tail is tensor-blind on this model too. **(g)** it is **outside** the `[TPeak, TPeakPlusLL]` interval on a large minority of cells. **(h)** `m_lepton` is read. | (a) 2e-4; (c) 2e-3; (d) 2e-3; (e) 1e-8 / 1e-3 |
| **T10** | **Deuteron form-factor normalisation.** `HoSpin1FF` built with the *deuteron's* moments returns `F_q(0) = 25.83` and `F_m(0) = 1.714`. | 1e-3 relative against the textbook values — proves the convention |
| **T11** | **⁶Li form-factor anchors.** `F_c(0) = 3`, `F_m(0) = 4.90765`, `F_q(0) = −65.914` (TUNL `Q = −0.0818(17) fm²`; pin whichever moment §2.1 finally adopts — Pyykkö's −0.0806 gives −64.947); `−6 dF_c/dq²|₀ / F_c(0)` reproduces `⟨r²⟩_point = 6.078 fm²`; `HoSpin1FF::for_ion` **throws** for ⁷Li and for the deuteron-by-name (no measured-moment block). The first `F_c` zero is pinned to **the refit result with the refit's own uncertainty**, and asserted only to lie in `[2.9, 3.3] fm⁻¹` — *not* to `3.10 ± 0.02`, which the first draft attributed to WS98 and which WS98 does not state. | 1e-4 relative on the values, the zero in `[2.9, 3.3]` |
| **T12** | **Quadrupole band is QUADRATIC, not linear.** The first draft asserted "exactly linear in `fq_scale` to 1e-12", which **fails against a correct implementation**: Eq. (38)'s `σ_q^d` contains `F_q(3F_c + 3ηF_m + ηF_q)` and `F_q(4F_c − 3xF_m + (4/3)ηF_q)`, and Eq. (A.4)'s `ℑ^el_2`/`ℑ^el_6` likewise carry `F_q·F_q` (only `ℑ^el_8` is linear). Instead: take the `Q_N`-dependent part of `w_tail − 1` at `s = fq_scale ∈ {0, 1, 2}`, fit `T(s) = T₀ + s·T₁ + s²·T₂` exactly (three points, three unknowns), and require the **prediction at `s = 0.5` and `s = 3`** to match the model. Also assert `T₂ ≠ 0` at some `(x, Q²)`, so a linear implementation fails. | 1e-12 |
| **T12'** | **The `F_m` band moves something.** `tail_tensor_scale ∈ {0.5, 1, 2}` changes the `Q_N`-dependent part of `w_tail − 1` by more than 1e-6 relative at `x = 0.1`, `Q² = 5` — i.e. `fq_scale` alone does **not** span the tensor tail. | as stated |
| **T13** | **Off is identically 1.** `RcModel` constructed with `RcMode::Off` (the mode is a **constructor argument**, not an `RcOptions` field — there is only one source of truth, §3.2): `weights()` returns `{1.0, 1.0, 1.0}` with `==`; `fill()` leaves `rc_weights` empty; `applies()` and `tail_applies()` are false; a full `--rc off` pipeline emits **no** `rc_*` npz key and a HepMC file whose `weight_names` are exactly `{"nominal"}` (+ the spin block). | exact |
| **T14** | **Per-channel table, over `Channel`.** For each of the **nine** `event.hpp` `Channel` values (not the five `PipelineChannel`s — the API is keyed on `Channel`, §3.3), assert the §1.5 row: `applies()`, `tail_applies()` and that the emitted weights are exactly 1.0 where the table says `≡ 1`, with `exclusion_reason()` non-empty and naming the reason. Separately assert `event_channel_of` maps all five `PipelineChannel`s and throws on nothing. **Also assert `tail_applies()` is TRUE at `θ_S = 90°`** on the inclusive channel, and that `rc_tail − 1` there is `−½` of its `θ_S = 0` value at the same `(x, Q², Q_NN)` (POLRAD Eq. (43); §1.4.7). | exact; the Eq. (43) ratio to 1e-12 |
| **T15** | **HepMC3 weight names round trip.** In `tests/test_hepmc.cpp`, use a **separate fixture event** for this (or update the two existing `REQUIRE`s at lines 324 and 331 to `1 + spin_weights.size() + rc_weights.size()`). Adding `rc_weights` to the shared fixture — whose `spin_weights` are set at line 40 — breaks both of those subcases, and the first draft's "extend" did not say so; an implementer following it gets a red suite. Then: `spin_weights = {1.05, 0.95, 1.02, 0.98}` **and** `rc_weights` of 3 × 5 entries, read back with `ReaderAscii`, `run_info->weight_names()` exactly `{"nominal", "spin_weight_1..4", "rc_tensor_lo", "rc_tensor_hi", "rc_tail", "rc_tensor_lo_1", …, "rc_tail_4"}` in that order, every value round-tripping. | exact names; values to 1e-12 |
| **T16** | **pytest gating the binding.** `python/tests/test_rc.py`: `lipolgen.RC` maps both strings; `RcMode`, `RcOptions`, `RcWeights`, `RcModel`, `HoSpin1FF`, `rc_delta` are importable; `make_config(rc="tensor-band", rc_fq_scale=…, rc_tail_tensor_scale=…, rc_qe_suppression=…)` round-trips through `PipelineConfig.rc` / `.rc_options`; a 2000-event run emits the three npz columns with `0.5 < w < 2` and `lo + hi == 2`; the **records path** (`export.columns_from_events`) emits them too; the same run at `rc="off"` emits none and `export.rc_columns` returns exact ones. Final clause: `--rc-delta-low-x 0.19` gives `|hi − 1| **<=** |hi − 1|(0.30)` **everywhere** at `x < 0.01`, and **strictly `<`** on the subset with `τ ≠ 0` — the strict form alone fails on any event at a `b₁` zero crossing of `Li6B1(MillerB1)`, where both are exactly 0. | as stated |
| **T17** | **Weighted mode.** With `spin_weights` filled, `rc_weights` has `3 × (1 + N)` entries. **(a)** slot `1+k` equals slot 0 only for a category with a **pure** population vector — build a plan of three pure `SpinCategory`s for this assertion. Slot 0 uses the event's own **pure** state `m`; a per-category slot mirrors `weights_for` (`sampler.cpp:461-483`), i.e. the population **mixture** `Σ_m p_m W_m`, and on the `tensor-thirds` plan (`P_zz = 0.6` mixtures, `bookkeeping.cpp:109-117`) those differ. `Event::rc_weights`'s docstring states this. **(b)** the `1 : 1 : −2` ratio across the three `tensor_thirds` fills holds for the **numerator** `W_tensor` (linear in `P_zz`), not for `τ = W_tensor/W` whose denominator also depends on `P_zz`. So assert the ratio on `tensor_weights_for(...)` (or on `τ_k · w_k`), and assert `τ` itself against `τ_k = (P_zz,k/2)A_zz / [1 + (P_zz,k/2)A_zz]` per category — whose ratio is `−2·W₁/W₃`, not `−2`. | 1e-12 |
| **T18** | **Threading.** `Pipeline::for_each(sink, 4)` with RC on returns the same weights as `nthreads = 1`. | bit-for-bit |

---

## 6. Implementation steps (two Opus agents)

**Sequencing.** Agent 1 lands `include/lipolgen/rc.hpp` (API + constants +
doc-comments, no bodies) **first and alone**, on its own commit; agent 2 starts
from that header. After that the two work on disjoint file sets and no worktree
isolation is needed.

### Agent 1 — the core model (physics)

Files it owns:

* `include/lipolgen/rc.hpp` **(new)** — §3.1 in full.
* `src/core/rc.cpp` **(new)** — `rc_delta`; `HoSpin1FF` (+ `for_ion`),
  `TabulatedSpin1FF`; **the three Eq. (38) `η_A` quadratures** (`σ_u^A`,
  `σ_q^A`, and the Eq. (44) QRT built on `σ_u^p` with `(Z, N)` nucleon form
  factors); `RcModel` with the table build + interpolation, the two
  `tensor_fraction` forms, `tail_ratio`, `tail_ratio_at`, `weights`, `fill`.
  **No Appendix B in v0** — that is `RcTailModel::PolradFull`, §1.4.6,
  **implemented 2026-09-06** (`../run_2026-09-06/phase_B_numbers.md` §B2).
* `include/lipolgen/constants.hpp` — `M_ELECTRON`, and the one-line change in
  `src/hepmc/hepmc_writer.cpp:135` (as of 0145885) to use it.
* `include/lipolgen/xsec.hpp` + `src/core/xsec.cpp` —
  `InclusiveKernel::tensor_amplitudes` and the `amplitudes()` refactor.
* `include/lipolgen/sampler.hpp` + `src/core/sampler.cpp` — the three
  `StateTables` tensor vectors and `tensor_weights_for`.
* `tests/test_rc.cpp` **(new)** — T1–T5, T7–T14, T18.
* `data/ff/li6_elastic.csv` **(new, optional in v0)** + its README row.

Order of work:
0. **Prerequisite, not a test: obtain the POLRAD FORTRAN** (CPC Program
   Library, `ADXQ`) and check the three transcription flags of §1.4.2 —
   `σ_u^C`'s `Y₋`, `X₁`'s linear `−4η_A`, and the `(3/4)x²` / `(4/3)η_A F_q`
   fractions — plus Eq. (18)'s `1/A²`. An hour of work that removes the whole
   Q2 class of error. If it cannot be obtained, say so in the commit message
   and lean on T8(a)–(c).
1. `rc.hpp` (API only) → commit → hand to agent 2.
2. `rc_delta` + T3. Cheap, and it unblocks the CLI.
3. `tensor_amplitudes` refactor + `StateTables` vectors + `tensor_weights_for`,
   then **re-run the whole existing suite and the `validation/reference/*.json`
   rtol-1e-12 gates before writing anything else** — this refactor must move
   no number.
4. `tensor_fraction` (inclusive, then tagged) + T1, T2, T4, T17.
5. `HoSpin1FF::for_ion` (+ `F_mag`) + T10, T11.
6. The Eq. (38) `η_A` quadratures + **T8(a) first** (the deuteron gate is the
   one that fixes the normalisation), then T8(b) carbon `Z²`, T8(c) `A = 1`
   Mo–Tsai, then T7, T5, T12, T12'.
7. The Eq. (44) QRT + T8'.
8. `RcModel::fill`, the per-channel table + T13, T14, T18. T9 belongs with
   §1.4.6 and is skipped in v0 (Eq. (A.4) is not evaluated) — keep it in the
   file, `SKIP`-ped, with the reason.

Watch-outs, in the order they will bite:
* **`inclusive_b1` defaults false on the tagged channels** — the inclusive
  recipe returns `τ = 0` there. §1.5.1.
* **`Amplitudes` must not be recomputed per event** — read the sampler's
  `StateTables` at `Event::kin.cell`, exactly like the sampler drew it.
* **`target_mass = true` is the default** and needs a 96-point `g2^WW`
  quadrature per `tables()` call — another reason never to call `tables()` in
  the event loop.
* **`σ_u^A` and `σ_q^A` are separate integrals with the SAME kernel weights** —
  and `Q_N` enters Eq. (37) as `Q_N/6`, not `Q_N`. Dropping the `1/6` inflates
  the tensor tail sixfold and no unpolarised test sees it; T8(a) does.
* **Per-nucleon vs per-nucleus.** POLRAD's tail is per nucleus in
  `x_A = x/A`, `S_A = A·s`, `M_A`; the library's `dsigma_unpol` is per nucleon
  in per-nucleon `x`, `s`. `d²σ/dx dy = (1/A) d²σ/dx_A dy`. **A missing `A = 6`
  is the single most likely error in this file** and no 5 %-tolerance test can
  see it — which is why T8(a) compares against a literal Eq. (38) transcription
  written *in the test*, at 1e-10.
* **Eq. (18)'s prefactor is `α³y/A²` on the right and `1/A` on the left**
  (§1.4.6). Relevant since 2026-09-06, when `PolradFull` was implemented; the
  `ℑ₆` correction is now exercised by T9.
* **`t_min = 4M_A²η_min = M_A²x_A²/(1−x_A)`** must come out `≈ (x M_N)²`. If it
  scales with `A`, the nuclear map is wrong somewhere. Assert it.
* `RcModel` must be **immutable after construction** — `Pipeline::for_each`
  calls it concurrently, and `GlauberFsiWeight` is the precedent for how.

### Agent 2 — wiring, output, bindings, tests, docs

Files it owns:

* `include/lipolgen/event.hpp`, `src/core/event.cpp` — `rc_weights` + `reset`.
* `include/lipolgen/pipeline.hpp`, `src/core/pipeline.cpp` — `PipelineConfig::rc`
  and `rc_options`, `validate()`, `event_channel_of`, the constructor build
  **after the channel if/else chain, immediately before `// --- beams ---`**
  (§3.3 — NOT beside the FSI block, which is inside the tagged branch), the
  `Pipeline::event` hook after the channel switch and before the hadronizer,
  the `rc_model()` accessor and the `rc_` member.
* `src/hepmc/hepmc_writer.cpp` — `ensure_run_info` names and the per-event
  weight vector (§4.1).
* `python/bindings.cpp` — the `Columns` block, `fill_row` (line 201), `columns_to_dict`,
  and the whole binding block of §4.4.
* `python/lipolgen/__init__.py` (`RC`, `make_config`), `cli.py` (six switches:
  `--rc`, `--rc-delta-low-x`, `--rc-delta-high-x`, `--rc-fq-scale`,
  `--rc-tail-tensor-scale`, `--rc-qe-suppression`; `DEFAULTS`; the banner),
  `export.py` (`RC_KEYS`, `rc_columns`, **and the three keys in
  `columns_from_events`**).
* `src/hepmc/hepmc_writer.cpp:135` (as of 0145885) — replace the bare `0.51099895e-3` with
  `M_ELECTRON` from `constants.hpp` (agent 1 adds the constant; this is the
  one line of the change that lives in agent 2's file set, and it must land
  after agent 1's header commit).
* `tests/test_hepmc.cpp` — T15.
* `python/tests/test_rc.py` **(new)** — T16, and the pipeline half of T6, T13,
  T17.
* `docs/USAGE.md`, `docs/OPEN_ITEMS_SOLUTIONS.md` row 9 + §9,
  `docs/CONVENTIONS.md`, `docs/HEPMC3_CONVENTION.md` (the weight-name rule,
  lines 194-198), `docs/PHYSICS_CHANNELS.md` §11.

Order of work:
1. `Event::rc_weights` + `reset` (T13 half).
2. `PipelineConfig::rc` + `validate()` + the constructor build + the hook.
3. **T6 immediately** — the bit-for-bit gate. Everything after this is
   additive, and T6 is what proves it.
4. HepMC3 names + T15.
5. Bindings + npz columns + `export.rc_columns`.
6. CLI + `__init__.RC` + the banner + T16.
7. Docs, then the §8 numbers once agent 1's model is in.

Watch-outs:
* **Append, never insert**, in the HepMC weight vector — `nominal` stays index
  0 and every `spin_weight_k` keeps its index, or every downstream reader
  breaks.
* **Emit the npz columns only when RC is on**, so an `--rc off` run's key set
  is unchanged and the existing reference gates do not have to move — **on
  BOTH paths**, `columns_to_dict` (C++) and `export.columns_from_events`
  (Python). Missing the second one makes `rc_columns()` return ones for a run
  that had RC on.
* **`tests/test_hepmc.cpp` lines 324 and 331 will go red** if `rc_weights` is
  added to the shared fixture event. Use a separate event or move both
  `REQUIRE`s (T15).
* **`RcOptions` has no `mode`.** `PipelineConfig::rc` is the one source of
  truth and `RcModel` takes it as a constructor argument.
* `PipelineConfig::validate()` must **not** refuse the coherent channel (unlike
  FSI): `RcModel::applies()` handles it and the run prints the reason.
* `Pipeline` currently never fills `Event::spin_weights` — weighted mode is a
  latent capability. Slot 0 must be correct on its own, and the `1 + k` slots
  must degrade gracefully to nothing when `spin_weights` is empty.

**Effort.** `physics_literature.md` §3 prices item 9 at **5–8 person-days**
(+5 if POLRAD's full ℑ-based *inelastic*-tail integrals are transcribed, which
this design deliberately does not do). Agent 1 ≈ 4–6 d, agent 2 ≈ 2 d.

---

## 7. Open questions

* **Q1 — the ⁶Li C2 form factor has no measurement.** Elastic longitudinal
  data give `F_C0² + F_C2²` only, and C2 is sub-dominant exactly where the data
  are good. The only ab-initio C2 (WS98) has a quadrupole moment 3× the
  measured one **at `q → 0` while its high-`q` magnitude is right** (WS98: `F_L`
  is *"in excellent agreement with experiment"* and the `q ≥ 3 fm⁻¹` shoulder
  *"is entirely due to"* C2). The design therefore pins the two ends separately
  and carries the ±100 % `fq_scale` band on the **interpolation between them**.
  **Unresolved:** the shape in between. Would need either a digitisation of WS98
  Fig. 1 with the C0 subtracted from the Li et al. high-`q` data, or a new
  calculation.
* **Q2 — POLRAD transcription checks (RESOLVED for the prefactor, OPEN for
  three fractions).** Eq. (18) reads `−(α³y/A²)` on the right and
  `(1/A)d²σ^el/dx_A dy` on the left — **`1/A²` and `1/A`, not "`1/A` on both
  sides"** as the first draft said; Eqs. (19)/(20) carry a bare `1/A`, which is
  the contrast that settles it, and Eq. (21)'s `(1 − x/A)` fixes `x_A = x/A`.
  What remains open is the v0 formulas' fine print: `σ_u^C`'s `Y₋` (every other
  unpolarised entry has `Y₊`), `X₁`'s linear `−4η_A`, and `σ_q^d`'s `(3/4)x²`
  vs `(4/3)η_A F_q`. A copy of the POLRAD FORTRAN (CPC Program Library,
  `ADXQ`) settles all three in an hour and is listed as a **prerequisite**
  (§6 step 0), not a test.
* **Q3 — the quasi-elastic tail on a TAGGED channel. STILL OPEN, and its
  "probably small" is WITHDRAWN (2026-09-04, task B4).** The *unpolarised* QRT
  is priced inclusively (§1.4.3); its *polarised* part is neglected under Zhou
  et al. everywhere. On the α-tagged channel **`rc_tail ≡ 1` is half a
  kinematic fact**: the elastic half is one, the quasi-elastic half is an
  omission. The A−1 remnant is unbound, so the α comes out at `x_L ≈ 2/3` —
  the tag window — and when the struck nucleon belongs to the embedded
  deuteron the α is a **true spectator** carrying the same `n_M(k, c)` the
  tagged Born does. The tag therefore neither vetoes it nor suppresses it in
  the ratio, and the omitted dilution is of the **same order as the inclusive
  quasi-elastic one** (99.9 % of the inclusive tail at x = 0.30), not a
  negligible one. Nobody has computed it; pricing it needs a tagged Born
  denominator and the tag acceptance inside Eq. (44), plus the cluster-elastic
  `e + A → e' + γ + d + α` that Eq. (44) does not contain. Full kinematics in
  §1.5.2 and `../run_2026-09-03/phase_B_numbers.md` §B4.
* **Q4 — the tagged band is an extrapolation.** No RC calculation exists for a
  tagged tensor asymmetry (§1.5.1). Applying `δ(x)` to `τ_tag` is defensible
  and conservative but uncited.
* **~~Q5 — the transverse / general spin axis.~~ CLOSED.** The first draft
  zeroed the tail for `θ_S ≠ 0`, claiming POLRAD "does not carry" a
  non-longitudinal tensor axis. It does: Eq. (43) gives
  `σ_q⊥^d = −½ σ_q∥^d` and Eq. (6)'s basis includes `η_⊥`. More generally an
  azimuthally-integrated inclusive cross section can depend on a rank-2 target
  polarisation only through `P₂(cos θ_S)`, so the general axis is
  `Q_N → P_zz^eff = 3 Q_NN P₂(cos θ_S)` with no new machinery. §1.4.7. The one
  residual caveat — the φ-integrated tail applied per event to a density that
  carries `cos φ′`/`cos 2φ′` modulations at `θ_S ≠ 0` — is printed by the run
  and stated in USAGE, not left as an open question.
* **Q6 — φ-dependent tensor observables have no RC at all.** The `Δ` (cos 2φ)
  sector and the whole coherent channel are excluded for this reason, and the
  exclusion is permanent until somebody computes it. `RcScope::TensorAll`
  exists only to *price* the omission, never to correct it.
* **Q7 — is the band the right shape?** Log-linear between two anchors is a
  choice. Gakh–Shekhovtsova's Fig. 2 has the actual `x` shape at four `Q²`
  values, and a digitisation would replace the interpolation with the
  calculation's own — at the cost of inheriting an uncited calculation's shape
  as well as its magnitude. The design deliberately keeps the interpolation and
  the anchors visible.

  > **DECIDED 2026-09-04 (task B6): the Gakh–Shekhovtsova SHAPE IS REJECTED,
  > with numbers. Not overlooked — investigated and refused.** The digitisation
  > is *available*: the arXiv source of hep-ph/0403262 ships the four panels as
  > EPS files (`rrc01.eps`, `rrc1.eps`, `rrc4.eps`, `rrc10.eps` at
  > `Q² = 0.1, 1, 4, 10 GeV²`), each carrying a 30-point `/x` array and **two**
  > `/y` arrays — the dashed one is the RC-included `Δσ`, the solid one the
  > Born, per the Fig. 2 caption. Every number below is read off those arrays,
  > not off a picture. Four reasons it may not be adopted:
  >
  > **(a) `δ` is SINGULAR, not a bounded fraction.** The Born `Δσ` **crosses
  > zero** — the `b₁` zero-crossing — inside the `Q² = 4` and `Q² = 10` panels,
  > and the RC **moves the crossing**, which is the paper's own stated result
  > (*"the inclusion of the radiative correction shifts the zero value of `b₁`
  > and `b₂` to the smaller `x`-value region"*). Measured on the arrays:
  > `x₀ = 0.20118 → 0.18467` at `Q² = 4` and `0.20161 → 0.18331` at `Q² = 10`.
  > A ratio across a zero is unbounded. `δ ≡ (Δσ_RC − Δσ_Born)/Δσ_Born`
  > throughout — the RC as a fraction **of the Born**, which is the base the
  > paper's own quoted number names and the base the `w = 1 ∓ δτ` ansatz
  > applies; it is this record's construction from the two curves,
  > **inferred, not a symbol the paper defines** (`phase_B_numbers.md` §B6.3
  > sets it out, corrected 2026-09-04 from `Δσ_Born/Δσ_RC − 1`, and every
  > number here has been reconverted). It runs
  >
  > | `Q²` | `x` window | `δ` range on the panel |
  > |---|---|---|
  > | 0.1 | 0.00226 – 0.00966 | −0.266 … −0.113 |
  > | 1 | 0.01348 – 0.09662 | **−0.171 … +0.240** (changes SIGN) |
  > | 4 | 0.05394 – 0.38647 | **−1.691 … +7.810** |
  > | 10 | 0.13097 – 0.85000 | **−0.887 … +4.091** |
  >
  > and at `Q² = 4` specifically −0.170 (`x = 0.111`) → **−0.537**
  > (`x = 0.169`) → **+1.026** (`x = 0.214`) → +0.030 (`x = 0.283`). The band
  > ansatz is `w = 1 ± δτ` with `band_tau_max = 1`, which assumes a **small,
  > bounded, one-signed fractional rescale** and needs `|δ| < 1` for
  > `w_lo > 0` and `w_hi > 0` at all. A paper-shaped `δ` is none of those
  > things: it is
  > unbounded, it changes sign, and it would emit negative weights.
  >
  > **(b) The shape is NOT monotone non-increasing in `x`**, which is this
  > design's stated justification for the log-linear form (§1.3) and is
  > *asserted* by `tests/test_rc.cpp:70-83` — T3 checks monotone
  > non-increasing on a 2000-point log grid and compares **both anchors with
  > `==`**. Beyond the sign changes of (a), the panels give **−0.128 at
  > `x = 0.386`** (`Q² = 4`) and **−0.400 at `x = 0.85`** (`Q² = 10`) against
  > `δ_high = 0.015` for **all** `x ≥ 0.16` — **8.5×** and **26.7×** past the
  > E12-13-011 high-`x` anchor in magnitude, and negative, which a half-width
  > cannot be. The two anchors and the shape cannot both
  > stand; adopting the shape would silently delete the only *measured*
  > statement in the band.
  >
  > **(c) The `Q²` support is patchy and not rectangular.** Four `Q²` values
  > with barely-overlapping `x` windows (above), and panels **a and b do not
  > overlap at all** — there is a gap at `0.00966 < x < 0.01348`. This
  > generator's own quoted point, **`x = 0.01, Q² = 5`**, falls *between* the
  > `Q² = 4` and `Q² = 10` panels **and** outside the `x` window of both
  > (`0.01 < 0.05394` and `0.01 < 0.13097`) — and its `x` lands exactly in the
  > a–b gap. There is no panel to interpolate from at the point the design
  > publishes. **Related, and worth its own flag (corrected 2026-09-04):** the
  > *magnitude* actually taken, the paper's "10 % to 30 %", is quoted there for
  > `x ∼ 10⁻³–10⁻²`, and **only panel (a) lies in that window** — all 30 of its
  > points do, while panel (b) starts at `x = 0.01348`, above `10⁻²`, and
  > contributes none. So `δ_low = 0.30` is a **`Q² = 0.1 GeV²`** number quoted
  > at `Q² = 5`, not the `0.1–1 GeV²` this block previously said, and
  > `|δ| = 0.113…0.266` on that panel. The `0.171…0.240` published here as a
  > `Q² = 1` `|δ|` range was **not a magnitude range at all**: `−0.171` and
  > `+0.240` are that panel's two *signed* extremes under this `δ`, at
  > opposite ends of its `x` window, while its `|δ|` spans **0.006 … 0.240**. Taking the magnitude is already an
  > extrapolation; taking the shape would compound it.
  >
  > **(d) Their Born is not our Born.** The ratio is against **their** tensor
  > model, not this generator's `b₁`: their Eq. (60) HERMES `A_zz`
  > parametrisation `A_zz = −1.56·10⁻²(1 − 1.74x − 1.45√x)`, `F₂^d` from
  > ALLM97 [21] with `F₂ⁿ/F₂ᵖ` [22], and `b₂` from the Callan–Gross-type
  > relation Eq. (59). A *fractional* RC computed against a different `b₁`
  > shape does not transfer to a different `b₁` shape — that is a fresh,
  > unpriced assumption on top of the A = 2 → A = 6 one of Q8.
  >
  > **What is kept.** The log-linear interpolation, the two visible anchors,
  > and the *magnitude* `δ_low = 0.30` with its flags. What Q7 is really
  > asking — "is the shape right?" — has an answer that is **not** "digitise
  > Fig. 2": it is *"the source cannot supply a shape this band can use"*.
  > Q7 stays open; the digitisation route is now **closed**, and reopening it
  > requires a `δ` definition that survives a zero crossing.
* **Q8 — no `A > 2` tensor RC exists at all.** Everything cited is deuteron.
  ⁶Li enters only through its form factors (which we supply) and its `b₁`
  (Phase D). Whether the deuteron's *fractional* RC transfers to ⁶Li is
  untested and is the largest unquantified assumption in the band.

  > **PRICED 2026-09-04 (task B6), still not answered.** The assumption now has
  > a dial: `RcOptions::a_transfer_frac` / `--rc-a-transfer-frac`, a fraction of
  > `δ(x)` added **in quadrature** by `RcModel::delta` — the model's single call
  > site of `rc_delta`, so it is applied exactly once and on every path:
  >
  > ```
  >   delta_eff(x) = hypot( delta(x), a_transfer_frac * delta(x) )
  >                = delta(x) * sqrt(1 + a_transfer_frac^2)
  > ```
  >
  > **Default 0.0**, which is the assumption v0 shipped *unstated* — "the
  > deuteron fraction transfers exactly" — and the shipped output is bit for
  > bit unchanged (`hypot(d, 0) == |d|`, and `rc_delta ≥ 0`; measured, §B6).
  > `0.5` = "known to 50 % of itself"; `1.0` = "as uncertain as it is large",
  > a `√2`-wider band. **There is no measurement to prefer any of them** — this
  > is a price tag, not a correction, and `meta["rc_a_transfer_frac"]` records
  > which one ran. Measured half-widths on `A_zz` at `x = 0.01, Q² = 5`:
  > **1.3585e−04 / 1.5188e−04 / 1.9212e−04** at `f = 0 / 0.5 / 1`
  > (`run_2026-09-03/phase_B_numbers.md` §B6).
  >
  > It is deliberately **multiplicative and `δ`-proportional**: it must vanish
  > where `δ` does, because the E12-13-011 high-`x` anchor is an `A = 2`
  > measurement too — the doubt about the transfer cannot be larger than the
  > correction it doubts. Q8 stays **open**: no `A > 2` tensor RC calculation
  > exists, and none of `0`, `0.5`, `1` is derived from one.
* **Q9 — the quasi-elastic suppression factors.** v0 sets POLRAD Eq. (44)'s
  `S_E = S_M = S_EM = 1`, i.e. no Pauli blocking, banded by
  `qe_suppression ∈ {0, 0.5, 1}`. POLRAD's own peak prescription (its §2.1.2
  Eq. (20) and ref. [10]'s sum rules) would replace the knob with a `Q²`-
  dependent factor. Until it does, the QRT is overestimated at low `Q²` — the
  conservative direction for a dilution, but a real ≈ ±50 % on a piece that is
  30–70 % of the elastic tail.
* **Q10 — the ⁶Li magnetic form-factor shape.** `F_m` drives the sub-leading
  tensor terms and enters two interferences, and WS98's `F_T` has a zero and a
  second peak inside the `t` range that matters. v0 uses a three-landmark fit
  (peak 0.5, zero ≈ 1.2–1.4, peak 2 fm⁻¹) and bands it with
  `tail_tensor_scale`, which is a *normalisation* band, not a *shape* band. A
  digitisation of WS98 Fig. 1's `F_T` closes this.

---

## 8. Placeholders the implementer must fill

Fill these from a real run, then copy them into `docs/OPEN_ITEMS_SOLUTIONS.md`
§9 and `docs/USAGE.md`. **Configuration for every table: ⁶Li, `--config 1`
(10 GeV e × 99.5 GeV/u, so `s_per_nucleon = 4 E_e p_u ≈ 3980 GeV²`
— `beams.cpp:121-126`), `--channel inclusive`,
`--plan tensor-thirds --pzz 0.6`, `--rc tensor-band` at the defaults,
`θ_S = 0`.** **Every `(x, Q²)` below must satisfy `y = Q²/(x s) ≤ 0.985`**;
`s` is printed with each table. (The first draft's `x = 0.01, Q² = 50` needed
`y = 1.26` and does not exist.)

### 8.1 The radiative tails — `w_tail − 1` (elastic + unpolarised quasi-elastic, as a fraction of the Born)

Evaluate at the `m = +1` state (`P_zz = +1`) and at `m = 0` (`P_zz = −2`); the
difference between the two rows is the *tensor* part of the tail, which is the
whole point.

| x | Q² = 2 GeV² (y) | Q² = 5 GeV² (y) | Q² = 10 GeV² (y) |
|---|---|---|---|
| 0.01, `m = ±1` | **TBD** (0.050) | **TBD** (0.126) | **TBD** (0.251) |
| 0.01, `m = 0` | **TBD** | **TBD** | **TBD** |
| 0.10, `m = ±1` | **TBD** (0.005) | **TBD** (0.013) | **TBD** (0.025) |
| 0.10, `m = 0` | **TBD** | **TBD** | **TBD** |
| 0.30, `m = ±1` | **TBD** (0.002) | **TBD** (0.004) | **TBD** (0.008) |
| 0.30, `m = 0` | **TBD** | **TBD** | **TBD** |

and the `y`-scan that T5 is written against, at fixed `Q² = 5 GeV²`:

| y | 0.1 | 0.3 | 0.5 | 0.7 | 0.9 | 0.97 |
|---|---|---|---|---|---|---|
| x (= Q²/(y s)) | 0.0126 | 0.0042 | 0.0025 | 0.0018 | 0.0014 | 0.0013 |
| `w_tail − 1` | **TBD** | **TBD** | **TBD** | **TBD** | **TBD** | **TBD** |
| `Y₊(y)` | 10.9 | 4.13 | 3.50 | 3.63 | 4.51 | 5.83 |

and the derived quantity that actually matters:

| x | Q² = 2 | Q² = 5 | Q² = 10 | quantity |
|---|---|---|---|---|
| 0.01 | **TBD** | **TBD** | **TBD** | `(1/6)σ^el_T / σ^el_U` — the tensor fraction of the elastic tail (expect `O(10⁻²)` per unit `Q_N`, §2.1) |
| 0.10 | **TBD** | **TBD** | **TBD** | |
| 0.30 | **TBD** | **TBD** | **TBD** | |
| 0.01 | **TBD** | **TBD** | **TBD** | `σ^q_U / σ^el_U` — the QRT against the ERT (the §1.4.3 estimate is 0.3–0.7; **this table replaces it**) |
| 0.10 | **TBD** | **TBD** | **TBD** | |
| 0.30 | **TBD** | **TBD** | **TBD** | |
| 0.01 | **TBD** | **TBD** | **TBD** | `ΔA_zz` induced by the tail (Born → Born+tail) — **dominated by the unpolarised dilution** |
| 0.10 | **TBD** | **TBD** | **TBD** | |
| 0.30 | **TBD** | **TBD** | **TBD** | |
| 0.01 | **TBD** | **TBD** | **TBD** | the same `ΔA_zz` with `qe_suppression = 0` — i.e. how much of it the QRT is |
| 0.10 | **TBD** | **TBD** | **TBD** | |
| 0.30 | **TBD** | **TBD** | **TBD** | |

*Expected qualitative behaviour, to sanity-check against — rewritten; the first
draft's version was not implied by the physics.* The ⁶Li elastic form factors
are dead above `q ≈ 2.5 fm⁻¹` ⟺ `t ≈ 0.25 GeV²` (`q = 5.07 fm⁻¹` at
`t = 1 GeV²`). The tail therefore lives entirely in the **t-peak**, whose
lower limit is

```
  t_min = 4 M_A^2 eta_min = M_A^2 x_A^2/(1 - x_A) ~= (x M_N)^2       INDEPENDENT OF Q^2
```

so:

* the tail's `t` window is `[(x M_N)², ~0.25 GeV²]` and it **closes with
  increasing `x`**, not with increasing `Q²`. It shuts entirely at
  `x M_N ≈ 0.5 GeV`, i.e. `x ≈ 0.53`. **The `x`-dependence is the strong one.**
* at fixed `x` the `Q²`-dependence enters only through `Y₊` and the Born
  denominator, so the ratio is **not** "strictly decreasing in `Q²`" and
  certainly not "faster than `1/Q²`" — the first draft's T5, built on the
  s-peak picture POLRAD says is *suppressed*, asserted both.
* `Y₊ = [1 + (1−y)²]/(1−y)` **diverges as `y → 1`**, so the tail *grows* toward
  the `y` edge. "If any entry exceeds ~1 at `Q² ≥ 5` something is wrong" is
  false there; the honest statement is that `clipped_fraction_by_y()`'s
  `y > 0.9` band is where clipping is expected, and it must be reported.
* the **QRT** window is wider (nucleon form factors reach `t ∼ 1 GeV²`), so
  `σ^q_U/σ^el_U` should **grow with `x`** as the ERT window closes faster.

### 8.2 The band — `w_hi − 1` (per event, at the sampler's own `(x, Q²)`)

| x | `δ(x)` | `τ = (P_zz/2)A_zz / [1 + (P_zz/2)A_zz]` at `P_zz = +1` | `w_hi − 1` | `w_lo − 1` |
|---|---|---|---|---|
| 0.01 | 0.30000 | **TBD** | **TBD** | **TBD** |
| 0.063 | 0.11081 | **TBD** | **TBD** | **TBD** |
| 0.10 | 0.06331 | **TBD** | **TBD** | **TBD** |
| 0.16 | 0.01500 | **TBD** | **TBD** | **TBD** |
| 0.30 | 0.01500 | **TBD** | **TBD** | **TBD** |

**Note the τ column header.** The first draft's `τ = (P_zz/2)A_zz` omitted the
denominator; with `|A_zz|` of a few per cent that is a small error inclusively
and a **20–35 %** one on the tagged channels (§1.5.1).

### 8.3 Run-level summary numbers

| quantity | value |
|---|---|
| `RcModel` construction time (100 × 72 grid + y-refinement, `n_eta = 128`) | **TBD** s |
| throughput with `--rc tensor-band` vs `--rc off` | **TBD** vs **TBD** ev/s |
| `clipped_cell_fraction()` at `tail_max = 10`, global | **TBD** % |
| `clipped_fraction_by_y()` — `y < 0.5` / `0.5–0.9` / `> 0.9` | **TBD** / **TBD** / **TBD** % |
| fitted `(a, α, q₀)` and the fit's χ²/ndf against the Li/Suelzle data | **TBD** |
| fitted `(q_z, b)` of `F_mag` against WS98's three `F_T` landmarks | **TBD** |
| δ(A_zz) systematic at `x = 0.01` from the band, `δ_low = 0.19` vs `0.30` | **TBD** |
| δ(A_zz) systematic at `x = 0.01` from `fq_scale = 0` vs `2` | **TBD** |
| δ(A_zz) systematic at `x = 0.01` from `tail_tensor_scale = 0.5` vs `2` | **TBD** |
| δ(A_zz) systematic at `x = 0.01` from `qe_suppression = 0` vs `1` | **TBD** |

---

## References

* I. Akushevich, A. Ilyichev, N. Shumeiko, A. Soroko, A. Tolkachev, **POLRAD 2.0**,
  [arXiv:hep-ph/9706516](https://arxiv.org/abs/hep-ph/9706516), CPC **104** (1997) 201.
  **v0 uses:** §2.1.3 B (p. 10, the t-peak / suppressed s- and p-peaks statement),
  §2.1.4 (p. 11, the same for higher orders), Eqs. (2), (4), (5), (6), (9), (10),
  (21), **(37)**, **(38)**, **(39)**, (40), **(43)**, **(44)**.
  **The upgrade path (§1.4.6) additionally uses:** Eqs. (3), (13), (14), (17),
  **(18)**, (19), (20), (A.1)–(A.3), **(A.4)**, (A.5), (B.1)–(B.14) —
  including (B.5), (B.7) and (B.8).
* G. I. Gakh, O. Shekhovtsova, [arXiv:hep-ph/0403262](https://arxiv.org/abs/hep-ph/0403262),
  JETP **99** (2004) 898 — the only dedicated tensor-DIS RC calculation, Eqs. (57)–(60),
  Fig. 2, **§4** *Numerical estimations* ("10 % to 30 %"; **corrected 2026-09-04**
  from "§5" — in the arXiv source's own numbering §5 is Appendix A).
  **Zero INSPIRE citations.**
* HERMES, [arXiv:hep-ex/0506018](https://arxiv.org/abs/hep-ex/0506018), PRL **95** (2005)
  242001 — **Eq. (1)** (the rate form `σ_U[1 − P_z P_B D A₁ + ½ P_zz A_zz]`), **Eq. (2)**
  (the extraction, `A_zz = (2σ¹ − 2σ⁰)/(3σ_U P_zz^eff)`), Eq. (5), **Table II** (the
  measured `A_zz` per x bin); RADGEN unpolarised photon spectrum; *"the radiative
  background … reaches almost 50 % of the statistics in the lowest-x bin"*; **both**
  coherent and quasi-elastic tails subtracted with deuteron form-factor
  parameterisations; residual RC systematic ≈ 2×10⁻³ "for the three bins at low x"
  (= 19 %, 19 %, 15 % of the measured `A_zz` there).
* JLab **E12-13-011**, proposal **PR12-13-011** (Slifer et al.) — the 1.5 % RC line.
  **Unpublished; cite it by page or drop the anchor.** Poudel, Bacchetta, Chen,
  Santiesteban, [arXiv:2506.04506](https://arxiv.org/abs/2506.04506), EPJ A **61**
  (2025) — the experiment's **kinematics, `0.16 < x < 0.49`, `0.8 < Q² < 5.0 GeV²`,
  `W ≥ 1.85 GeV`** (p. 8), which is where `RC_X_HIGH = 0.16` comes from. It contains
  **no** radiative-correction discussion; the "9.2 % total" of the first draft is
  withdrawn.
* **Z.-L. Zhou et al.** (NIKHEF), PRL **82** (1999) 687 — a **quasi-elastic `e-d`
  tensor-analysing-power (`T₂₀`) measurement**, HERMES's ref. [13] for neglecting
  the polarised QRT. (Not "Zhou, Beane, Ji" — the first draft's author list was a
  guess, left with a question mark.)
* R. B. Wiringa, R. Schiavilla, [arXiv:nucl-th/9807037](https://arxiv.org/abs/nucl-th/9807037),
  PRL **81** (1998) 4317 — ⁶Li `F_L` (C0 + C2) and `F_T`, **figures only, no table**;
  `Q(⁶Li)_VMC = −0.23(9) fm²` vs measured −0.08; `F_L` *"in excellent agreement with
  experiment"* and the `q ≥ 3 fm⁻¹` shoulder *"entirely due to"* C2 (so the VMC C2 is
  right at high q and 3× too big at q → 0); `F_T` has its **first peak at
  `q = 0.5 fm⁻¹`, a zero, and a second peak at `q = 2 fm⁻¹`**; two-body currents
  *"shift the minimum to lower values of q"* with **no number given** — WS98 does not
  locate the first C0 zero.
* **TUNL A = 6 evaluation** (`refs/TUNL_A6_2002.pdf`) — `μ(⁶Li) = +0.8220473 μ_N`,
  `Q(⁶Li) = −0.818(17) mb = −0.0818(17) fm²` (1998CE04). The commonly quoted
  **−0.0806(6) fm²** is **P. Pyykkö**, Mol. Phys. **106** (2008) 1965 (from
  Cederberg et al. 1998), **not** TUNL.
* B. S. Pudliner et al., [arXiv:nucl-th/9705009](https://arxiv.org/abs/nucl-th/9705009),
  PRC **56** (1997) 1720, Table XIV — `Q(⁶Li) = −0.33(18)` vs exp. −0.083; the standing
  "do not derive a ⁶Li tensor input from these wave functions".
* L. R. Suelzle, M. R. Yearian, H. Crannell, Phys. Rev. **162** (1967) 992;
  G. C. Li, I. Sick, R. R. Whitney, M. R. Yearian, Nucl. Phys. **A162** (1971) 583 —
  the ⁶Li elastic electron-scattering data.
* I. Angeli, K. P. Marinova, ADNDT **99** (2013) 69 — `r_ch(⁶Li) = 2.589(39) fm`.
* W. Cosyn, A. Roldan Tomei, D. Sosa, A. Zec, [arXiv:2410.12764](https://arxiv.org/abs/2410.12764),
  EPJ A **61** (2025) 83, Eq. (27) — `A_T = −(2/3) b₁/F₁`, the `TENSOR_LL_SIGN` anchor.
* In-repo: `docs/open_items/physics_literature.md` §3 (the recommendation this design
  implements), `docs/CONVENTIONS.md`, `include/lipolgen/fsi.hpp` (the weight-only
  precedent), `docs/OPEN_ITEMS_SOLUTIONS.md` row 9 and §5.
