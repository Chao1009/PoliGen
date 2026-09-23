<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# 07 — The Cosyn–Weiss tensor gate: the mapping, the S–D sign, and what moves

**Investigation date 2026-09-06.** Seventh file in the benchmarking series, and the
resolution of `BENCHMARK_PLAN.md` row 2 / `00_in_tree_checks.md` row `E-1` —
"the single strongest external physics check in the tree", returned **refuted at
blocker severity** by `06`.

Everything below was re-derived from the sources on the investigation date. The
repository was **not modified**: every number is either a `pdftotext` extraction
of the papers, a read of the shipped `data/`, or a scratch probe compiled
against the shipped `build/libLiPolGenCore.so`. The probes reproduce
`TaggedModel::n_of_kc` **bit for bit** (max relative difference `0.000e+00` on
every channel and every `M`) before changing anything, so the "fixed" column is
a like-for-like measurement of the same code with one phase applied.

---

## 0. Verdict, in four lines

| question | answer |
|---|---|
| Is `A_T∥ = −2 · A_zz^wf` (tests/test_tagged.cpp:492-548) right? | **No. The mapping is `+1`, exactly, for every θ_k, every f₂/f₀, every α_p.** |
| Is the S–D interference sign in `n_M(k)` the CW/CDKS sign? | **No — it is the opposite.** It is not a convention that cancels: it flips a shipped observable. |
| Bug or registry decision? | **Bug.** It contradicts a published formula, the repository's *own* quadrupole sign gate, and the repository's *own* b₁ sector. |
| One-line fix | `src/core/tagged.cpp:413` — multiply the per-wave angular factor by `((l/2)%2==0) ? 1.0 : -1.0`, i.e. the relative part of `φ_L = i^L ψ_L`, exactly as `ClusterPartialWave::from_vmc` (`src/core/b1_nuclear.cpp:262`) already does. |

The headline consequence, measured on the shipped defaults:

> **⁶Li `A_zz^tag` at k = 0.20 GeV moves from +0.845 to −1.207 (Hulthén β = 0.30)
> and from +0.452 to −0.519 (VMC AV18).** The corrected VMC column is,
> digit for digit, the column the tree currently ships as the **`vmc-flipD`
> control** — `docs/open_items/vmc_reconciliation.md:180-188`. The control is
> the physics; the physics is the control.

---

## 1. What Cosyn–Weiss II actually says

`arXiv:2603.23700`, extracted with `pdftotext -layout` and `pdftotext` (raw) over
pp. 5-7, 16-22 and 33-37. Verbatim, with the raw extraction reassembled:

**Eq. (2.28)–(2.30), p. 6 — a pure spin state along `N`:**

```
S_D = Λ N ,        T_D = W(Λ) ( −1/6 + 1/2 N⊗N ) ,      W(Λ) ≡ (1, −2, 1) for Λ = (+1, 0, −1)
N = e_z:   S_D = Λ e_z ,  T_D = W(Λ) [ (1/√6) e_LL ] ,  {T_LL, T_LT, T_TT} = W(Λ) × {1/3, 0, 0}
```

**Eq. (3.22b), p. 9 — the spin-1 LF wave function:**

```
Ψ~_D(k, λ'_n, λ'_p | Λ) = (1/√2) ε^i_D(Λ) [ δ^ij f₀(k) + (1/√2)(3 k^i k^j/|k|² − δ^ij) f₂(k) ]
                          χ⁺(λ'_n) σ^j (iσ²) χ*(λ'_p)
```

There is **no i^L anywhere in Eq. (3.22b)**; the D-form enters with coefficient
`+1/√2`; Eq. (3.24) fixes `f_L(k) = √E(k) f_{L,nr}(k)`, so `f₂/f₀` is the
non-relativistic ratio exactly.

**Eq. (4.17)–(4.20a), p. 18:**

```
P_[U](α_p, p_pT | T_D) = (f₀² + f₂²)/(2−α_p) − [3/(2−α_p)] (k T_D k)/|k|² (2f₀ + f₂/√2)(f₂/√2)
P_[U,U]      = (f₀² + f₂²)/(2−α_p)                                                    (4.18)
P_[T_LL,U]   = −[3/(2−α_p)] (2f₀ + f₂/√2)(f₂/√2) (k e_LL k)/|k|²                      (4.19)
(k e_LL k)/|k|² = (1/√6)(3cos²θ_k − 1)                                                (4.20a)
P_[U](T_D)   = P_[U,U] + √(3/2) T_LL P_[T_LL,U] + [T_LT, T_TT terms]                  (4.22)
```

**Eq. (6.11)–(6.14) and TABLE II, p. 35:**

```
A_T∥ = (2/3) [F_(U T_LL,T)D + ε F_(U T_LL,L)D] / [F_(U U,T)D + ε F_(U U,L)D]          (6.11)
     = √(2/3) P_[T_LL,U] / P_[U,U]
     = [ (2f₀ + f₂/√2)(f₂/√2) / (f₀² + f₂²) ] (1 − 3cos²θ_k)                          (6.12)
   bracket max = +1  at  f₂/f₀ = +√2                                                  (6.13)
   bracket min = −1/2 at  f₂/f₀ = −1/√2                                               (6.14)
   "One can easily verify that the expression takes values in [−2, 1]."
   angular factor: +1 at θ_k = π/2, −2 at θ_k = 0, π
TABLE II (AV18 radial wave functions):
   A_T∥   f₂/f₀     k        θ_k    |1−α_p|   p_pT
   −2     +√2      0.3 GeV   0       0.3       0
   +1     +√2      0.3 GeV   π/2     0        0.3 GeV
   +1     −1/√2    1 GeV     0       0.7       0
   "the condition Eq. (6.13) is satisfied for k = 0.30 GeV"
   "The condition Eq. (6.14) is satisfied only at k ≈ 1 GeV"
```

**Internal consistency of the extraction.** Substituting (4.19)+(4.20a) into
`√(2/3) P_[T_LL,U]/P_[U,U]` needs the coefficient `√(2/3) × 3/√6`, which is
**exactly 1** — so the two lines of Eq. (6.12) agree and nothing was lost in the
PDF layout. `Q(+√2) = +1.0000000` and `Q(−1/√2) = −0.5000000` reproduce
Eqs. (6.13)/(6.14) to machine precision, and the numeric extremum search over
`r ∈ [−50, 50]` finds `+1.0000000 at r = +1.41422` and `−0.5000000 at
r = −0.70710`. The extraction is sound.

---

## 2. `A_T∥ = +1 · A_zz^wf`, derived twice

### 2.1 From CW's own definitions (analytic)

Put the polarization axis along `z`. `n_Λ` is the unpolarized-neutron density in
deuteron spin state `Λ`, i.e. `P_[U]` at `T_D = W(Λ)(1/√6) e_LL`. From (4.22) and
(2.30b), `√(3/2) · T_LL = √(3/2) · W(Λ)/3 = W(Λ)/√6`, so

```
n_Λ = P_[U,U] + (W(Λ)/√6) P_[T_LL,U]
```

Now use `W(+1) + W(−1) − 2W(0) = 1 + 1 + 4 = 6` and `W(+1) + W(0) + W(−1) = 0`:

```
n_+1 + n_−1 − 2 n_0 = (6/√6) P_[T_LL,U] = √6 P_[T_LL,U]
n_+1 + n_−1 +   n_0 = 3 P_[U,U]                       (the tensor part cancels identically)

A_zz^wf ≡ (n_+1 + n_−1 − 2n_0)/(n_+1 + n_−1 + n_0) = (√6/3) P_[T_LL,U]/P_[U,U]
        = √(2/3) P_[T_LL,U]/P_[U,U]
        = A_T∥                                        [Eq. (6.12), second line, verbatim]
```

`√6/3 = 0.816496580927726 = √(2/3)`. **The coefficient is `+1`.** The `−2` the
test applies is the *angular factor at θ_k = 0*, which `A_zz^wf` already carries
through its own `P₂(cos θ_k)`; applying it again double-counts.

### 2.2 From CW's Eq. (3.22b) directly (numeric, convention-free)

Build `n_Λ = Σ_{m_s} |ε^i_Λ M^{ij} ε^{j*}_{m_s}|²` with
`M^{ij} = (f₀ − f₂/√2) δ^ij + 3(f₂/√2) k̂^i k̂^j`. Because `M` is real symmetric
and the spherical basis is unitary, `n_Λ` is *independent of the spin-basis
convention*. Then form the `(+1, +1, −2)` ratio:

| f₂/f₀ | θ_k | A_zz from Eq. (3.22b) | CW Eq. (6.12) | diff |
|---|---|---|---|---|
| −0.70711 | 0 | +1.000000000 | +1.000000000 | +1.1e−16 |
| −0.70711 | π/2 | −0.500000000 | −0.500000000 | −3.9e−16 |
| −0.30000 | 0 | +0.695897374 | +0.695897374 | −1.1e−16 |
| +0.30000 | 0 | −0.861034988 | −0.861034988 | −3.3e−16 |
| +0.70711 | 0 | −1.666666667 | −1.666666667 | 0 |
| +1.00000 | π/3 | +0.239276695 | +0.239276695 | −1.9e−16 |
| +1.41421 | 0 | −2.000000000 | −2.000000000 | −2.2e−16 |
| +1.41421 | π/2 | +1.000000000 | +1.000000000 | +1.1e−16 |
| +3.00000 | 0 | −1.748528137 | −1.748528137 | −4.4e−16 |

`A_T∥ = A_zz^wf` to machine precision, everywhere. **Coefficient `+1`, not `−2`.**

### 2.3 What each Table II row therefore demands of `A_zz^wf`

| row | k, θ_k | CW A_T∥ | required A_zz^wf (mapping +1) | required under the tree's −2 |
|---|---|---|---|---|
| 1 | 0.30 GeV, 0 | −2 | **−2.000** | +1.000 |
| 2 | 0.30 GeV, π/2 | +1 | **+1.000** | −0.500 |
| 3 | 1 GeV, 0 | +1 | **+1.000** | −0.500 |

At `f₂/f₀ = +√2` the `Λ = ±1` densities vanish identically
(`n_±1 = (f₀ − f₂/√2)² = 0`, `n_0 = (f₀ + √2 f₂)² = 9 f₀²`, giving
`A_zz = −18/9 = −2.000000`) — CW's own "this is the location where the polarized
neutron distributions have a node", under Eq. (6.13). The `−2` is a **node**, not
a conversion factor.

---

## 3. Which sign is CW's `f₂`? — pinned on the tree's own AV18 table

This is the only place the argument could still hide, and CW pin it themselves
by quoting two `k` landmarks for the AV18 wave function. Measured on the shipped
`data/vmc/deuteron/fdeut.av18` k-block (`ħc = 0.1973269804 GeV·fm`):

| condition | with f₂ = +W(k) | with f₂ = −W(k) | CW's stated k |
|---|---|---|---|
| f₂/f₀ = +√2 (Eq. 6.13) | **k = 0.2984 GeV** | k = 0.6455 GeV | "k = 0.30 GeV" |
| f₂/f₀ = −1/√2 (Eq. 6.14) | **k = 1.0257 GeV** | k = 0.2373 GeV | "k ≈ 1 GeV" |

`f₂ = +W` reproduces both landmarks; `f₂ = −W` reproduces neither.
**CW's `f₂` is `+W(k)`, the plain (i^L-stripped) spherical-Bessel transform,
positive at low k.** (U's node is at k = 0.4130 GeV = 2.093 fm⁻¹.)

Supporting checks on the same file, so the convention is nailed on both sides:

```
r-block:   u(r) > 0 and w(r) > 0 at every one of the 10000 rows
           Q_d = (1/20)∫dr r² w (2√2 u − w) = +0.269670 fm²   [header qm = 0.269673]
           P_D = ∫w² dr = 0.057599                            [header dstate = 0.057599]
           with w → −w:  Q_d = −0.269670 fm²
k-block:   file/(bare Bessel transform of the r-block) = 0.79788 for BOTH L=0 and L=2,
           at k = 0.3, 1.0, 2.0, 3.0, 5.0 fm⁻¹;  √(2/π) = 0.79788
           ⇒ the k-block carries NO i^L and W(k) > 0 at low k.
```

So `w(r) > 0` is pinned by the **measured** positive quadrupole moment, and
`W(k) > 0` is the plain transform of it.

---

## 4. What `build_amp2` actually computes

`src/core/tagged.cpp:388-428`:

```cpp
ang[w][ic] = cg_of[w] * theta_lm(l_of[w], ml_int, c_[ic]);        // line 413
...
amp += (*rad[w])[ik] * ang[w][ic];                                 // line 421
tab[ik * nc + ic] = amp * amp;
```

i.e. `A(M, m_s) = Σ_L a_L(k) ⟨L, M−m_s; S, m_s | J, M⟩ Θ_L^{M−m_s}(c)`,
`n_M = Σ_{m_s} A²` — CDKS Eq. (19) with **no i^L**, with
`a_0 = ψ_0 = U`, `a_2 = ψ_2 = +W` (`src/core/tagged.cpp:201-220` for AV18; the
analytic Hulthén D-form `k²/((k²+κ²)(k²+β²)²)`, `cluster.hpp:19`, is
positive-definite).

Measured decomposition, using the library's own `theta_lm` and
`clebsch_gordan` (10-digit fit):

```
n_M = (1/4π) [ a₀² + a₂² + β_M P₂(c) ]
β_+1 = β_−1 = +1.4142135624 a₀a₂ − 0.5000000000 a₂²  =  √2 a₀a₂ − a₂²/2
β_0            = −2.8284271247 a₀a₂ + 1.0000000000 a₂²  =  −2 β_+1
coefficient of a₀²:  0.0000000000  (exactly)

⇒ A_zz^wf = 6 β_+1 P₂ / (3(a₀²+a₂²)) = P₂(c) · g(x),   x ≡ a₂/a₀,
   g(x) = (2√2 x − x²)/(1 + x²)
```

CW Eq. (6.12) in the same algebra is `A_T∥ = P₂(c) · (−2) Q(r)` with
`Q(r) = (√2 r + r²/2)/(1+r²)`, and

```
g(x) = −2 Q(−x)      ⇒   A_zz^wf(a₂ = +f₂) = A_T∥ evaluated at MINUS f₂/f₀
                     ⇒   A_zz^wf(a₂ = −f₂) = A_T∥ exactly
```

Verified against the shipped library: with the L=2 amplitude negated,

```
max | A_zz^wf − CW Eq.(6.12) | over the FULL (k, c) grid  (280 × 96 cells)
    deuteron AV18 control :  shipped 2.740e+00   i^L-applied 8.882e-16
    deuteron Hulthén ctrl :  shipped 2.740e+00   i^L-applied 8.882e-16
```

**With the i^L phase the identity `A_T∥ = A_zz^wf` holds to machine precision on
every cell of the shipped grid. Without it, the curve is CW's formula at the
wrong sign of `f₂/f₀`.** That settles both questions at once: the mapping is
`+1`, and the tree's D amplitude carries the wrong sign.

Why `i^L` is not negotiable: `ψ̃(k) = ∫d³r e^{∓ik·r} ψ(r)` gives `(∓i)^L`, and
`i² = (−i)² = −1`, so for `L = 2` the factor is `−1` under **either** Fourier
sign. It is not a global phase — it is `+1` on `L = 0` and `−1` on `L = 2`,
a *relative* phase inside one sum. The only escape would be redefining
`Y_{211}` so that `w(r) < 0`, which flips `Q_d` negative and contradicts §3.
CDKS state it directly (`1702.05337`, below Eq. 21): *"the D-state wave function
has the negative sign (φ₂(p) < 0) due to the i^L factor"*.

### 4.1 The tell that needs no paper at all

Two facts internal to the tree prove the sign independently.

**(i) The gate's own landmark is at 1/√2, not √2.** `tests/test_tagged.cpp:521-527`
asserts *"the k envelope reaches its maximum 1 at f₂/f₀ = √2"* and finds the peak
at k = 0.3098 GeV. Measured on the very model it runs on
(`deuteron_channel()`, Hulthén β = 0.30):

```
at k = 0.3098 GeV the model's f₂/f₀ = 0.705343      (1/√2 = 0.707107;  √2 = 1.414214)
Hulthén f₂/f₀ over the whole grid: 0.0000 … 1.2866  — it NEVER reaches √2
```

`g` peaks at `x = 1/√2` precisely because `g(x) = −2Q(−x)` and `Q` *minimises*
at `r = −1/√2`. The gate is matching CW's Eq. (6.14) **minimum** with a flipped
`f₂` and calling it Eq. (6.13)'s **maximum**. The agreement with "k = 0.30 GeV
for AV18" is a numerical coincidence of β = 0.30 on a channel that is not AV18.

**(ii) The tree's own quadrupole sign gate refuses it.** `alpha_d_quadrupole_fm2`
(`include/lipolgen/b1_nuclear.hpp:740-758`) is documented as *"on the deuteron
pair it must return `fdeut.av18`'s own header value qm = 0.269673 fm²"*. Fed the
two candidate amplitudes:

```
φ₂ = −W   (CDKS / b1_nuclear sign)   →  Q_d = +0.269362 fm²    [header +0.269673]
φ₂ = +W   (what build_amp2 sums)     →  Q_d = −0.303639 fm²
```

The tagged sector's D amplitude is the one that gives the deuteron a **negative**
quadrupole moment. (−0.3036 is the same number `cluster.hpp` records for CD-Bonn
"with the other sign".)

---

## 5. The two sectors disagree with each other

`src/core/b1_nuclear.cpp:262` and `:277` apply the phase explicitly
(`phase = ((l/2)%2==0) ? 1.0 : -1.0`, guard message *"i^L is real for even L
only"*), then `:379` recovers `const double wv = -sp2[j];  // W(k) = -phi_2(k)`
and forms CDKS Eq. (21):

```
δ_T f  ∝  3 [ U W/√2 + W²/4 ] (3c*² − 1)          SD and DD terms with the SAME sign
```

The tagged sector's tensor combination, from §4's exact decomposition:

```
n_+1 + n_−1 − 2n_0 = (3/4π) ( √2 a₀a₂ − a₂²/2 ) (3c² − 1)
   a₂ = +W  (shipped) :  (6/4π) [ +U W/√2 − W²/4 ] (3c² − 1)   SD and DD OPPOSITE
   a₂ = −W  (fixed)   :  −(6/4π) [ U W/√2 + W²/4 ] (3c² − 1)   SD and DD SAME  ✓
```

The `W²` term is phase-blind, so it is the fixed reference: **as shipped the two
sectors carry opposite S–D interference signs, and after the fix they agree.**
`b1_nuclear` is the correct one. The reviewer's `b1_nuclear.hpp:335` observation
is confirmed, with the polarity resolved: b₁ is right, tagged is wrong.

The false step is recorded in the tree in two places, and it is the same
sentence both times — `include/lipolgen/tagged.hpp:228-230` and
`src/core/tagged.cpp:197-200`:

> *"CDKS fix φ_L = i^L ψ_L with φ₂ = −W and U, W ≥ 0 at low k, so ψ₂ = +W"*

— true, and then used to feed `ψ₂` into a partial-wave sum that needs `φ₂`.
`ψ₂ = +W` is right; `build_amp2` is simply not the place it belongs.

---

## 6. What moves — measured on the shipped defaults

Probes rebuild `build_amp2` from `TaggedModel::radial_table(l)` and reproduce
`n_of_kc` **bit for bit** (`max rel diff 0.000e+00`, every channel, every `M`)
before applying the phase.

### 6.1 Does not move (angle-integrated; L-orthogonality kills the cross term)

| quantity | shipped | i^L-applied | relative |
|---|---|---|---|
| deuteron `tensor_dilution` | 0.9594880739 | 0.9594888782 | +8.4e−07 |
| deuteron `vector_dilution` | 0.9324947691 | 0.9324961093 | +1.4e−06 |
| ⁶Li Hulthén `tensor_dilution` | 0.9219467089 | 0.9219489770 | +2.5e−06 |
| ⁶Li Hulthén `vector_dilution` | 0.8699393995 | 0.8699431789 | +4.3e−06 |
| ⁶Li VMC `tensor_dilution` | 0.9825758170 | 0.9825758723 | +5.6e−08 |
| deuteron AV18 `tensor_dilution` | 0.9481456181 | 0.9481463394 | +7.6e−07 |
| `population_integrated`, all channels, all M | — | — | ≤ 5e−07 abs |
| `p2_moment_mixture_uniform` | −5.425347e−05 | −5.425347e−05 | 0 |

> **CORRECTION, 2026-09-06 (verification pass) — the last two rows are wrong at
> the magnitudes they quote, and the paragraph below oversold three claims.**
> Re-measured as the re-pinned `validation/reference/tagged.json` against
> `git show HEAD:validation/reference/tagged.json` (the pre-fix dump), and
> independently on the library with the stored L = 2 table negated:
>
> | quantity | channel | shipped | i^L-applied | move |
> |---|---|---|---|---|
> | `population_integrated` | ⁶Li Hulthén, M = 0, m_S = 0 | 0.947966373524 | 0.947963349333 | **3.02e−06 abs** |
> | `population_integrated` | deuteron Hulthén, M = 0, m_S = 0 | 0.972992754631 | 0.972991682225 | 1.07e−06 abs |
> | `norm` | ⁶Li Hulthén, M = 0 | 1.000026689626 | 0.999968571520 | **5.81e−05 rel** |
> | `norm` | deuteron Hulthén, M = 0 | 1.000018994706 | 0.999979287452 | 3.97e−05 rel |
> | `p2_moment_mixture_uniform` | ⁶Li Hulthén | −5.3160743e−05 | −5.2153460e−05 | **1.89e−02 rel** |
> | `p2_moment_mixture_uniform` | deuteron Hulthén | −5.3694953e−05 | −5.3337762e−05 | 6.65e−03 rel |
> | `p2_moment_mixture_uniform` | ⁶Li VMC | −5.4245146e−05 | −5.4220571e−05 | 4.53e−04 rel |
> | `p2_moment_mixture_uniform` | ⁷Li | −4.3408664e−05 | −4.3408664e−05 | **0** |
>
> So `population_integrated` is ≤ **3.0e−06** abs, not ≤ 5e−07; `norm` is up to
> **5.8e−05** rel, not 2e−5; and `p2_moment_mixture_uniform` moves on every
> spin-1 channel — the row's −5.425347e−05 pair is the ⁶Li VMC value to seven
> figures and is not a "does not move" entry at all. All four quantities had to
> be **re-pinned** in `validation/reference/tagged.json`, which is what the
> re-pin record now says (`validation/README.md`, generated from
> `dump_polligen_reference.py`). ⁷Li is the one exact zero.

The residuals are the 96-cell midpoint quadrature of `∫Θ_0 Θ_2 dc = 0`, nothing
more. **The tagged tensor dilution does not move** beyond that residual
(+2.5e−6 on the ⁶Li Hulthén channel, against the 1e−4 at which
`LI6_B1_RANK2_TRANSFER` is pinned), so `LI6_B1_RANK2_TRANSFER`, the effective
polarizations and the inclusive sector are untouched. `norm` moves by up to
**5.8e−05** relative (see §7 — it trips a 1e−9 gate and was re-pinned).

> **The tag fractions are NOT in that list — corrected 2026-09-06.** A tag
> fraction is spin-blind only for a spin-blind or category-averaged fill. The
> uniform-M mix accepted fraction is unmoved to all 15 digits
> (**0.024675932148828** Hulthén β = 0.30, **0.033810227625842** VMC AV18,
> acceptance-weighted Σ_M n_M k² at `n_phi = 32`, YR high-acceptance,
> `default_configs("6Li")[1]`), and the equal-thirds `tensor-thirds` average
> reproduces it to +3.1e−6 (a per-M norm-residual effect — Σ n_M k² differs by 5e−5 between M = 0 and ±1 — not summation order; the pooled category average equals the uniform mix to 2.2e−16 on both builds). For a **tensor-polarised** fill the tag fraction is
> the same acceptance-weighted integral as A_zz^tag and moves with it: measured
> on the ⁶Li Hulthén channel at those optics, `tensor-thirds` categories
> **0.028127 → 0.018403** (azz±) and **0.017775 → 0.037222** (azz0), and the
> `--pz 0.7` max-entropy ladder of `helicity-flip` **0.027030 → 0.020396**
> (−24.5 %); ⁶Li VMC 0.035296 → 0.032155, deuteron Hulthén 0.025802 → 0.021492,
> deuteron AV18 0.022175 → 0.019407. ⁷Li does not move at all.

**⁷Li moves by exactly zero**: one L = 1 wave, no interference. `n_of_kc`
identical to the bit; `p2_moment(M) = −0.1999349117 / +0.1998480944` unchanged.

### 6.2 Moves, and by how much

| quantity | channel | shipped | i^L-applied | change |
|---|---|---|---|---|
| `p2_moment(M=+1)` | deuteron Hulthén | +0.0320226969 | **−0.0411257151** | sign flip, ×1.284 |
| `p2_moment(M=0)` | deuteron Hulthén | −0.0642064787 | **+0.0820914169** | sign flip, ×1.279 |
| `p2_moment(M=+1)` | ⁶Li Hulthén | +0.0448133305 | **−0.0622510957** | sign flip, ×1.389 |
| `p2_moment(M=0)` | ⁶Li Hulthén | −0.0897861433 | **+0.1243457311** | sign flip, ×1.385 |
| `p2_moment(M=+1)` | ⁶Li VMC | +0.0038621121 | **−0.0078394851** | sign flip, ×2.030 |
| `p2_moment(M=+1)` | deuteron AV18 | +0.0198147059 | **−0.0314365982** | sign flip, ×1.587 |
| `n_of_kc` at the reference points | deuteron Hulthén | — | — | up to **+719 %** |
| `n_of_kc` at the reference points | ⁶Li Hulthén | — | — | up to **+725 %** |
| `n_of_kc` at the reference points | ⁶Li VMC | — | — | up to **+1453 %** |
| `n_of_kc` at the reference points | deuteron AV18 | — | — | up to **+19215 %** |
| `struck_populations` p(m_S\|M,k,c) | deuteron Hulthén | — | — | up to **0.860** absolute |
| `struck_populations` | ⁶Li Hulthén | — | — | up to **0.839** absolute |
| `norm(M=+1)` | deuteron Hulthén | 0.9999928764 | 1.0000127301 | +2.0e−05 |

### 6.3 `A_zz^wf(k)` at the θ_k = 0 and θ_k = π/2 cells

Deuteron control, Hulthén β = 0.30 (`deut_model()`), `c = +0.989583` / `−0.010417`:

| k [GeV] | θ≈0 shipped | θ≈0 fixed | θ≈90° shipped | θ≈90° fixed |
|---|---|---|---|---|
| 0.0990 | +0.344326 | −0.378651 | −0.177629 | +0.195336 |
| 0.1979 | +0.827412 | −1.111630 | −0.426841 | +0.573462 |
| 0.2968 | +0.967648 | −1.575777 | −0.499185 | +0.812903 |
| 0.4001 | +0.937986 | −1.778248 | −0.483883 | +0.917353 |
| 0.5979 | +0.837952 | −1.891958 | −0.432278 | +0.976013 |
| 0.9979 | +0.743817 | −1.928114 | −0.383716 | +0.994664 |

Deuteron control, **AV18** (`ClusterWaveSource::VmcAV18`) — CW's own wave function:

| k [GeV] | f₂/f₀ | θ≈0 shipped | θ≈0 fixed | CW Eq. (6.12) here |
|---|---|---|---|---|
| 0.1979 | +0.7676 | +0.875862 | −1.220872 | −1.220872 |
| 0.2968 | +1.3942 | +0.658141 | **−1.937694** | −1.937694 |
| 0.3012 | +1.4619 | **+0.617018** | **−1.937119** | −1.937119 |
| 0.4001 | (past the U node) | −0.799965 | −1.130776 | −1.130776 |
| 0.5979 | — | −1.918141 | +0.480791 | +0.480791 |
| 0.9979 | −0.7426 | −1.656113 | **+0.967340** | +0.967340 |

The reviewer's headline number is confirmed exactly: **AV18, k = 0.30 GeV,
θ_k ≈ 0: LiPolGen ships +0.617 where Eq. (6.12) gives −1.937** (−2 at the exact
θ_k = 0, `f₂/f₀ = √2` point; the outermost cell centre is `c = 0.98958`,
`P₂ = 0.968913`).

### 6.4 The headline: `A_zz^tag(k)`, ⁶Li, YR high-acceptance, 10 × 99.5 GeV/u

Acceptance-weighted `azz_tensor_curve_weighted` at the shipped optics
(`default_configs("6Li")[1]`, `yr_optics(..., high_acceptance = true)`,
`n_phi = 32`):

| k [GeV] | Hulthén β=0.30 shipped | Hulthén β=0.30 **fixed** | Hulthén P_D=VMC shipped | Hulthén P_D=VMC **fixed** | VMC AV18 shipped | VMC AV18 **fixed** |
|---|---|---|---|---|---|---|
| 0.1979 | +0.8450 | **−1.2069** | +0.5108 | **−0.6001** | +0.4518 | **−0.5191** |
| 0.2495 | +0.6551 | **−1.0639** | +0.4639 | **−0.5767** | +0.2604 | **−0.2899** |
| 0.3012 | +0.5226 | **−0.9533** | +0.4296 | **−0.5615** | +0.1819 | **−0.1993** |
| 0.4001 | +0.3325 | **−0.7267** | +0.3435 | **−0.4839** | +0.0779 | **−0.0822** |
| 0.4990 | +0.2349 | **−0.5851** | +0.2840 | **−0.4214** | −0.1302 | **+0.1170** |

The "shipped" columns are `docs/open_items/vmc_reconciliation.md:180-188` exactly.
The **fixed VMC column is that document's `vmc-flipD` column, digit for digit**
(−0.5191, −0.2899, −0.1993, −0.0822, +0.1170). The tree's designated *control on
the S–D relative sign* is the physical answer.

Documentation rows that move with it:
`docs/USAGE.md:1353` "⁶Li A_zz^tag at k = 0.20 GeV | +0.845 | +0.452"
→ **−1.207 | −0.519**; `docs/CONVENTIONS.md:328` and
`docs/open_items/vmc_reconciliation.md:133` "negating the D table flips A_zz^tag
from +0.45 to −0.52" → the sentence survives but the two numbers swap roles;
`docs/PHYSICS_CHANNELS.md:268` reconciliation values 0.845 / 0.511 / 0.452 /
−0.519 → −1.207 / −0.600 / −0.519 / +0.452.

### 6.5 What does **not** move outside `tagged`

`b1_nuclear` builds its own amplitudes through `ClusterPartialWave::from_vmc`,
which already applies `φ_L = i^L ψ_L`, and only *reads* the tagged channel's
`Wave::vmc` tables (which the fix leaves untouched — see §7). **No b₁ number
moves**: not `b1_default_li6.json`, not the CDKS shape gate, not
`alpha_d_quadrupole_fm2`, not `LI6_B1_RANK2_TRANSFER`. The inclusive sector,
the coherent channel, `boost_spectator` and every dilution beyond its
quadrature residual are unaffected.

> **CORRECTION, 2026-09-06 (verification pass) — three items were wrongly in
> that list.** (1) **The RC sector is spin-dependent and does move.**
> `RcModel::tagged_tau` (`src/core/rc.cpp:1565-1580`) reads
> `TaggedModel::n_of_kc(M, k, c)` per ion projection — only the M-average is
> isotropic — so under `--rc tensor-band` on a spin-1 tagged channel
> `rc_tensor_lo/hi` move in **every** event (max 0.488 absolute) and the
> band-clip count moves 96 → 248 (⁶Li Hulthén), 32 → 20 (⁶Li VMC) and 57 → 129
> (deuteron), measured over 20 000-event runs on a pre-fix build against this
> tree. RC on ⁷Li is unreachable (refused under any ⁷Li plan). (2) **Tag
> fractions** are invariant only for a spin-blind or category-averaged fill —
> see the box in §6.1. (3) **`population_integrated` itself moves** by up to
> 3.0e−06 absolute, and the polarised per-category cross sections it normalises
> move with it: σ_per_category_pb 590826.8808411484 → 590826.8764129529 (⁶Li
> Hulthén, 7.5e−9 relative) and 581694.0988441816 → 581694.0940292849 (deuteron
> Hulthén, 8.3e−9). Physically null, but it is a bound, not a zero.

What *does* change at event level: `TaggedModel::sample_kc` draws from
`|A_{m_S}(M; k, c)|² k²`, so the generated correlation between the ion
projection `M`, the struck-cluster projection `m_S` and the spectator direction
`cos θ_k` flips. The SPIN-BLIND rate does not move — unmoved to ≤ 4.0e−16
relative, 1–2 ulp, the VMC channel bit-identical — while a polarised category's
cross section moves at the 1e−8 level through the `population_integrated`
quadrature residual above; the angular correlation is what changes. That is the
whole point of the observable.

---

## 7. The fix

**One line.** `src/core/tagged.cpp:410-414`, inside `build_amp2`:

```cpp
    for (std::size_t w = 0; w < rad.size(); ++w) {
      for (std::size_t ic = 0; ic < nc; ++ic) {
        ang[w][ic] = cg_of[w] * theta_lm(l_of[w], ml_int, c_[ic]);   // line 413
      }
    }
```

becomes

```cpp
    for (std::size_t w = 0; w < rad.size(); ++w) {
      // phi_L = i^L psi_L.  Parity fixes L mod 2 inside one channel, so the
      // common i^(L mod 2) is an unobservable GLOBAL phase and the observable
      // RELATIVE phase is real: (-1)^floor(L/2).  Identical rule to
      // ClusterPartialWave::from_vmc / from_uw (src/core/b1_nuclear.cpp:262,
      // :269), which is why the b1 and tagged sectors then agree.
      const double phase = ((l_of[w] / 2) % 2 == 0) ? 1.0 : -1.0;
      for (std::size_t ic = 0; ic < nc; ++ic) {
        ang[w][ic] = phase * cg_of[w] * theta_lm(l_of[w], ml_int, c_[ic]);
      }
    }
```

`L = 0 → +1`, `L = 1 → +1` (⁷Li: a single wave, so the common `i` is global and
unobservable — measured: ⁷Li `n_of_kc` and `p2_moment` identical to the bit),
`L = 2 → −1`.

**Do the fix here, not in the radial tables.** Negating `rad_[2]` or the stored
`Wave::vmc->psi()` would be numerically equivalent but would break two tests that
*correctly* assert the stored `ψ` convention — `"tagged: the VMC S-D interference
flips sign below the S node"` (`test_tagged.cpp:1146`, checks
`radial_table(2)` signs about the S node) and `T25` (`test_tagged.cpp:1268-1269`,
checks `waves[1].vmc->psi()[1] > 0`). `ψ₂ = +W` is right; only the amplitude sum
needs `φ`.

Alongside it, two doc comments must be corrected because they carry the wrong
inference verbatim: `include/lipolgen/tagged.hpp:228-230` and
`src/core/tagged.cpp:197-200` ("…so ψ₂ = +W: the physical deuteron carries
sign(ψ₂/ψ₀) = +1, which is exactly what the positive-definite Hulthén forms
assume"). The Hulthén forms' positive-definiteness is fine; what was missing is
that `build_amp2` consumes `φ`, not `ψ`.

### 7.1 Tests that would fail (and are *supposed* to)

| test | file:line | why |
|---|---|---|
| `tagged: the model grid and its tables against polligen` | test_tagged.cpp:132 | `n_of_kc` (rtol), `struck_populations` (rtol), `p2_moment` (rtol), `norm` (1e−9 vs a 2e−5 move) all read `validation/reference/tagged.json` |
| `tagged: the Cosyn-Weiss deuteron tensor gate (CW TABLE II)` | test_tagged.cpp:492 | (a) mean ratio 0.99940 → −1.639538; (b) envelope peak → no interior peak; (c) both extremes; (d) the ×(−2) range check |
| `tagged: the deuteron S/D interference shape` | test_tagged.cpp:444 | `CHECK(m.k()[jbest] < 0.45)` — the corrected \|A_zz\| is monotone in k on this channel, so `jbest` runs to k = 1.2 |
| `tagged: the acceptance-weighted curve reduces to the cell curve` | test_tagged.cpp:551 | the last two lines, `CHECK(ref[jk] < -0.4)` and `CHECK(...weighted...[jk] > 0.4)`, both flip sign (→ +0.891 and −0.953) |

Tests that keep passing, checked by construction: the pure-S and pure-D limits
(no interference), the ⁷Li block (one wave), `every channel is normalized`,
`sum_M n_M is isotropic`, `boost_spectator`, the dilution tests, `T24`–`T27`,
`the VMC S-D interference flips sign below the S node`, `sample_kc reproduces its
own density` (it compares against the model's own `amp2_table`), and
`pipeline: A_zz^tag(k) closes on the acceptance-weighted truth`
(test_pipeline.cpp:350 — a generator-vs-own-model closure, sign-blind).

### 7.2 Reference JSON

`validation/reference/tagged.json` must be regenerated with
`validation/dump_polligen_reference.py`. Entries that move:
`channels.{deuteron,li6_alpha}.model.{n_of_kc, struck_populations, p2_moment,
norm}`. Entries that do **not**: everything under `li7_alpha.model`,
`population_integrated`, `vector_dilution`, `tensor_dilution`,
`p2_moment_mixture_uniform`, `boost_spectator`, `waves`, `base`,
`beam_configs`. `validation/reference/b1_default_li6.json` does not move.
`_manifest.json` hashes update accordingly.

> **Caution — the Python side is the reference generator.** `tagged.json` is a
> port gate against `polligen`'s `tagged.py`. If that Python carries the same
> missing `i^L`, regenerating from it will re-bake the wrong sign. The fix must
> land in *both* implementations, or the JSON must be regenerated from the
> corrected C++ and the port note updated to say so.
>
> **Resolved.** The fix landed in both implementations (LiPolGen 2026-09-06;
> `polligen` in PolarizedLithiumSim `1066555`, 2026-09-15); the JSON was
> regenerated from the corrected C++ from 2026-09-06 to 2026-09-23 and from the
> corrected `polligen` since, and `validation/dump_polligen_reference.py` now
> refuses a `polligen` without the phase
> (`../open_items/run_2026-09-23/phase_A_port_gate.md`).

### 7.3 Bug fix, not a registry decision

A phase convention that no observable sees is a registry decision. This is not
one:

1. It disagrees with a **published formula** — CW Eq. (6.12) — by up to 2.74 in
   an asymmetry whose whole range is `[−2, 1]`.
2. It disagrees with the repository's **own quadrupole sign gate**: the shipped
   tagged amplitude gives the deuteron `Q_d = −0.3036 fm²` against the measured
   `+0.2697`.
3. It disagrees with the repository's **own b₁ sector**, which applies `i^L`
   explicitly and documents the sign as *"LOAD-BEARING"*.
4. It **flips a shipped number**: `A_zz^tag(0.20 GeV)` for ⁶Li, +0.845 → −1.207
   and +0.452 → −0.519.

Any one of these is sufficient. Together they are decisive.

---

## 8. The gate as it should be written — like-for-like, on the AV18 control

CW's TABLE II is *"based on the AV18 radial wave functions"*. The current gate
runs on `deuteron_channel()`, which is the **Hulthén** analytic pair at the
scenario `P_D = 0.045` — a channel whose `f₂/f₀` never reaches `√2` anywhere on
the grid (§4.1). That is the like-for-likeness failure `06` flagged; the fix is
to run the CW rows on `ClusterWaveSource::VmcAV18` and keep only the
wave-function-independent parts on Hulthén.

All targets below are **measured** on the shipped `build/` with the phase applied.

```cpp
TEST_CASE("tagged: the Cosyn-Weiss deuteron tensor gate (CW TABLE II)"
          * doctest::skip(!vmc_data_present())) {
  // Cosyn-Weiss II (arXiv:2603.23700) p. 35, Eqs. (6.11)-(6.14) and TABLE II.
  //
  // THE MAPPING IS +1, NOT -2.  A_T|| = sqrt(2/3) P_[T_LL,U]/P_[U,U], and with
  // CW Eq. (2.30b) {T_LL,T_LT,T_TT} = W(Lambda) x {1/3,0,0}, W = (1,-2,1), the
  // (+1,+1,-2) combination of P_[U](T_D) gives sqrt(6)/3 = sqrt(2/3) times the
  // same ratio.  A_T|| IS our A_zz^wf.  The (1 - 3cos^2 theta_k) = -2 P2 factor
  // is INSIDE A_zz^wf already; CW's "-2" is the value at theta_k = 0, where the
  // Lambda = +-1 densities have a node (Eq. 6.13), not a conversion factor.
  //
  // TABLE II is quoted for the AV18 radial wave functions, so the CW rows run
  // on the AV18 control.  The Hulthen pair carries only the two statements that
  // do not depend on which wave function it is.
  const TaggedModel v(deuteron_channel(BETA_DEFAULT, P_D_DEUTERON,
                                       ClusterWaveSource::VmcAV18));
  const TaggedModel& h = deut_model();

  // ---- (a) the mapping, as an IDENTITY on every cell of the grid.
  // A_zz^wf(k,c) == [(2 f0 + f2/sqrt2)(f2/sqrt2)/(f0^2+f2^2)] (1 - 3 cos^2),
  // with f_L the model's own normalized radial tables.  Machine precision:
  // this is the whole of Eq. (6.12), not a shape comparison.
  for (const TaggedModel* m : {&v, &h}) {
    double worst = 0.0;
    for (std::size_t ik = 0; ik < m->nk(); ++ik) {
      const double f0 = m->radial_table(0)[ik], f2 = m->radial_table(2)[ik];
      const double r = f2 / f0;
      if (!std::isfinite(r)) continue;
      const double q = (2.0 + r / std::sqrt(2.0)) * (r / std::sqrt(2.0)) / (1.0 + r * r);
      const std::vector<double>& n1 = m->n_of_kc(1.0);
      const std::vector<double>& n0 = m->n_of_kc(0.0);
      const std::vector<double>& nm = m->n_of_kc(-1.0);
      for (std::size_t ic = 0; ic < m->nc(); ++ic) {
        const std::size_t j = ik * m->nc() + ic;
        const double a = (n1[j] + nm[j] - 2.0 * n0[j]) / (n1[j] + nm[j] + n0[j]);
        const double c = m->c()[ic];
        worst = std::fmax(worst, std::fabs(a - q * (1.0 - 3.0 * c * c)));
      }
    }
    CHECK(worst < 1e-12);          // MEASURED 8.882e-16 on both channels
  }

  // ---- (b) the P2 factorization is exact: A_zz^wf/P2 is a function of k alone.
  // (Unchanged from the old gate (a); it is a real property and it survives.)
  std::vector<double> p2(v.nc());
  for (std::size_t i = 0; i < v.nc(); ++i) p2[i] = 0.5 * (3.0 * v.c()[i] * v.c()[i] - 1.0);
  const std::size_t ik30 = argmin_abs(v.k(), 0.30);
  CHECK_CLOSE_AT(v.k()[ik30], 0.3012, 0.0, 5e-4);
  double rmin = 1e300, rmax = -1e300;
  for (std::size_t ic = 0; ic < v.nc(); ++ic) {
    if (std::fabs(p2[ic]) <= 1e-3) continue;
    const double r = azz_tensor_curve(v, ic)[ik30] / p2[ic];
    rmin = std::fmin(rmin, r); rmax = std::fmax(rmax, r);
  }
  CHECK(rmax - rmin < 1e-5);
  CHECK_CLOSE_AT(0.5 * (rmin + rmax), -1.99928, 0.0, 1e-4);   // MEASURED -1.999276

  // ---- (c) CW's OWN k landmarks, on CW's own wave function.
  // "Eq. (6.13) [f2/f0 = +sqrt2] is satisfied for k = 0.30 GeV"
  // "Eq. (6.14) [f2/f0 = -1/sqrt2] is satisfied only at k ~ 1 GeV"
  //   MEASURED on the model grid: 0.29812 GeV and 1.03487 GeV
  //   (raw fdeut.av18 k-block: 0.2984 GeV and 1.0257 GeV)
  CHECK_CLOSE_AT(k_where_ratio(v, +std::sqrt(2.0)), 0.30, 0.0, 0.01);
  CHECK_CLOSE_AT(k_where_ratio(v, -1.0 / std::sqrt(2.0)), 1.00, 0.0, 0.05);

  // ---- (d) TABLE II's two extremes, through the P2 factorization, with the
  // +1 mapping and NO extra factor anywhere.
  const std::size_t ic0 = argmax_abs(v.c());          // c = +0.989583
  const std::vector<double> env = [&] {
    std::vector<double> e(v.nk());
    const std::vector<double> curve = azz_tensor_curve(v, ic0);
    for (std::size_t i = 0; i < v.nk(); ++i) e[i] = curve[i] / p2[ic0];
    return e;
  }();
  std::size_t jlo = 0, jhi = 0;
  for (std::size_t i = 1; i < v.nk(); ++i) {
    if (env[i] < env[jlo]) jlo = i;
    if (env[i] > env[jhi]) jhi = i;
  }
  CHECK_CLOSE_AT(env[jlo], -2.0, 0.0, 1e-3);          // MEASURED -1.999864, CW: -2
  CHECK_CLOSE_AT(v.k()[jlo], 0.2968, 0.0, 5e-4);      //  at f2/f0 = +1.3942 (-> sqrt2)
  CHECK_CLOSE_AT(env[jhi], +1.0, 0.0, 1e-3);          // MEASURED +0.999997, CW: +1
  CHECK_CLOSE_AT(v.k()[jhi], 1.0366, 0.0, 5e-4);      //  at f2/f0 = -0.7056 (-> -1/sqrt2)
  CHECK_CLOSE_AT(env[jlo] * (-0.5), +1.0, 0.0, 1e-3); // the theta = pi/2 row of TABLE II

  // ---- (e) TABLE II row by row, at the actual cell centres.
  const std::size_t ic90 = argmin_abs(v.c(), 0.0);    // c = -0.010417
  CHECK_CLOSE_AT(azz_tensor_curve(v, ic0)[ik30],  -1.93712, 0.0, 2e-3);  // CW -2
  CHECK_CLOSE_AT(azz_tensor_curve(v, ic90)[ik30], +0.99931, 0.0, 2e-3);  // CW +1
  const std::size_t ik100 = argmin_abs(v.k(), 1.00);
  CHECK_CLOSE_AT(azz_tensor_curve(v, ic0)[ik100], +0.96734, 0.0, 2e-3);  // CW +1

  // ---- (f) the whole curve stays inside CW's stated range [-2, 1].
  for (std::size_t ic = 0; ic < v.nc(); ++ic)
    for (double x : azz_tensor_curve(v, ic)) {
      if (std::isnan(x)) continue;
      CHECK(x > -2.001);
      CHECK(x <  1.001);
    }

  // ---- (g) REGRESSION GUARD.  The pre-2026-09-06 amplitude summed psi_2 = +W
  // with no i^L, which is CW Eq. (6.12) evaluated at MINUS f2/f0.  Pin the sign
  // that distinguishes them, so the old behaviour cannot come back silently:
  // at k = 0.30 GeV, theta_k ~ 0, on AV18, CW is NEGATIVE (-1.937); the old
  // code gave +0.617.
  CHECK(azz_tensor_curve(v, ic0)[ik30] < -1.5);
}
```

Helper used above (`k_where_ratio`): linear interpolation of the first crossing
of `radial_table(2)/radial_table(0)` through the target, restricted to
`|ratio − target| < 1` so the `U`-node blow-up at k = 0.4130 GeV is not picked up.

Notes on what this gate now buys, and what it still does not:

- It is a **true external check**: (a) is CW Eq. (6.12) as an identity to 1e−12,
  on the wave function CW used, and (c)/(d)/(e) are CW's own quoted `k` values
  and TABLE II entries. Nothing in it is a self-consistency restatement.
- It still does **not** touch CW's Fig. 13, which is in light-front `(α_p, p_pT)`
  variables the generator does not carry — the caveat `00 §5` already records.
- The Hulthén channel keeps only (a) and (b), which is honest: it is not AV18
  and cannot carry AV18's landmarks. The old gate's "k = 0.3098 GeV vs CW's
  0.30" was a coincidence at `f₂/f₀ = 0.7053`, not a wave-function statement.
- (g) is the regression guard `06 §1` asks every restored gate to carry: it fails
  loudly if the `i^L` is ever dropped again.

---

## 9. Reproduction

Scratch probes (compiled read-only against the shipped library; nothing in the
repository was written except this file):

```
g++ -std=c++17 -O2 -I/home/cpeng/Projects/polli/LiPolGen/include -o pN pN.cpp \
    -L/home/cpeng/Projects/polli/LiPolGen/build -lLiPolGenCore \
    -Wl,-rpath,/home/cpeng/Projects/polli/LiPolGen/build
LIPOLGEN_DATA_DIR=/home/cpeng/Projects/polli/LiPolGen/data ./pN
```

| probe | what it shows |
|---|---|
| `p1.cpp` | the exact `β_M` quadratic form and the closed form `A_zz^wf = P₂ g(a₂/a₀)` |
| `cwwf.py` | `n_Λ` built straight from CW Eq. (3.22b) reproduces Eq. (6.12) to 1e−16 |
| `p2.cpp` | shipped vs i^L-applied on all five channels; the bit-for-bit rebuild check |
| `p3.cpp` | `max\|A_zz − CW Eq.(6.12)\|` over the full grid; TABLE II row by row |
| `p4.cpp` | every `tagged.json` model quantity, shipped vs fixed |
| `p5.cpp` | the acceptance-weighted `A_zz^tag(k)` table of §6.4 |
| `p6.cpp` | `alpha_d_quadrupole_fm2` on both candidate D-amplitude signs |
| `p7.cpp` | the `f₂/f₀` landmarks, and the 1/√2-vs-√2 tell |

Paper text: `pdftotext -layout` and `pdftotext` (raw, for the equation bodies) on
`/home/cpeng/Projects/polli/PolarizedLithiumSim/refs/2603.23700.pdf` pp. 5-7,
16-22, 33-37, and `.../1702.05337.pdf` (CDKS Eqs. 19-21).
Data: `data/vmc/deuteron/fdeut.av18` r-block and k-block, read directly.
