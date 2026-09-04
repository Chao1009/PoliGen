# Phase A3 — the b₁ normalisation arbiter: what each camp is per

**Item.** `docs/OPEN_ITEMS_SOLUTIONS.md`:543-545 (Q5) and :568-571 ("What has to
happen" items 5 and 6). Is `B1_PER_DEUTERON_TO_PER_NUCLEON = 0.5`
(`include/lipolgen/constants.hpp:86`) right at **both** sites it is applied —
`src/core/sf.cpp:350-354` (Miller table, reached by `MillerB1`/`toy_b1`) and
`src/core/sf.cpp:356-361` (CDKS convolution, reached by `CdksB1`/`b1_convolution`)?

**Read-only on code.** This document is the whole deliverable; nothing under
`src/`, `include/`, `tests/`, `python/` or `validation/` was touched. Baseline
re-measured before and after with no edits in between: `build/lipolgen_tests` →
**363 cases / 17203863 assertions / 1 skipped, SUCCESS**.

**Verdict in one line.** The arbiter is HERMES, and HERMES's published b₁ᵈ is
**per nucleon**. **CDKS is per nucleon too and must stop being halved.** **Miller,
by his own Eqs. (1)/(5)/(6)/(20), is per deuteron and the 0.5 stays** — but his
paper contradicts itself on exactly this factor, so the Miller half is *likely*,
not *certain*, and the residue needs the author. The two camps therefore **do not
share** `constants.hpp:86` any more.

| camp | code site | today | verdict | confidence |
|---|---|---|---|---|
| Miller | `sf.cpp:350-354` | ×0.5 | **keep ×0.5** | likely — his paper is self-inconsistent by exactly 2 |
| CDKS | `sf.cpp:356-361` | ×0.5 | **drop to ×1** | certain — CDKS say "per nucleon" in words |
| HERMES data | — | — | **per nucleon** | certain — their Eq. (5) recipe, confirmed numerically in all six bins |

---

## 1. Retraction

`docs/OPEN_ITEMS_SOLUTIONS.md`:545 says the answer to Miller's normalisation

> "is on the axis of Miller's Fig. 5"

and `:571` repeats it ("Check **Miller's** own normalisation independently, off
the axis of his Fig. 5"). **That is false and I am retracting it.** The framing
sent this item looking in a place where no answer exists:

* Miller's Fig. 5 caption is bare — *"FIG. 5: Computed values of b₁ = b₁^π + b₁^{6q}
  from Eq. (20) and Eq. (26). The pion structure function is that of [29], model 1"*
  (arXiv:1311.4561 p. 10, verified against the rendered page, not only against
  `pdftotext`).
* The ordinate is **`100 b₁(x)`** and nothing else. A scale factor of 100 carries
  no A-normalisation information whatsoever.
* The string "per nucleon" **never appears anywhere in the paper**
  (`grep -n "per nucleon\|per-nucleon" ` over the full `pdftotext -layout` dump
  returns only "non-nucleonic", "nucleon variable", "nucleon mass",
  "pion-nucleon", "nucleon-nucleon"). Neither does "1/A".

So the axis is silent. What Fig. 5 *does* carry, and what turned out to matter, is
the **six HERMES data points overlaid on it** — evidence of a different kind,
handled in §5. The right arbiter was never Miller's axis; it was always the
normalisation of the HERMES b₁ᵈ that *both* camps plot against.

A second, smaller correction while I am here: `include/lipolgen/constants.hpp:86`'s
comment —

> `/// The published b1 curves are per DEUTERON; every consumer here pairs b1`
> `/// with a per-NUCLEON F1, so the tables are halved on the way out.`

— states as a blanket fact something that is true for at most one of the two
curves. §4 shows it is flatly false for CDKS.

---

## 2. The arbiter: what HERMES normalises to

**Paper.** A. Airapetian *et al.* (HERMES), *First Measurement of the Tensor
Structure Function b₁ of the Deuteron*, PRL **95** (2005) 242001,
[arXiv:hep-ex/0506018v2](https://arxiv.org/abs/hep-ex/0506018). It was **not** in
`refs/`; I fetched it (5 pages) and read it. It is worth adding to `refs/` — it is
the only measurement either camp is calibrated against.

### 2.1 The QPM definition table (p. 2) — silent on A

HERMES write the three leading-twist structure functions in the QPM:

| | Nucleon | Deuteron |
|---|---|---|
| F₁ | ½ Σ_q e_q² [q↑^{1/2} + q↑^{−1/2}] | ⅓ Σ_q e_q² [q↑^1 + q↑^{−1} + q↑^0] |
| g₁ | ½ Σ_q e_q² [q↑^{1/2} − q↓^{1/2}] | ½ Σ_q e_q² [q↑^1 − q↓^1] |
| b₁ | — | ½ Σ_q e_q² [2q↑^0 − (q↑^1 + q↑^{−1})] |

with the gloss *"where q↑^m (q↓^m) is the number density of quarks with spin
up(down) along the z axis in a hadron(nucleus) with helicity m"*, and in the text
*"Because b₁ depends only on the spin averaged quark distributions
b₁ = ½(q⁰ − q¹), its measurement does not require a polarized beam"*.

Read literally ("in a … nucleus"), those densities are **per deuteron**. But note
what this table actually fixes: it fixes the *helicity combination*, not the
A-scaling — the same formulae hold verbatim if you feed them per-nucleon densities.
**The table is not where the published numbers get their normalisation.** Eq. (5)
is.

### 2.2 Eq. (5) — where the normalisation is actually set

HERMES p. 4, verbatim (the leading minus sign is a dropped glyph in the arXiv PDF's
broken embedded fonts; it is restored below — see §2.4):

> "The tensor structure function b₁ᵈ is extracted from the tensor asymmetry using
> the relations [18,27]
>
>   b₁ᵈ = −(3/2) A_zz^d F₁ᵈ ;  F₁ᵈ = (1 + Q²/ν²) F₂ᵈ / [2x(1 + R)]   (5)
>
> … **The structure function F₂ᵈ is calculated as F₂ᵈ = F₂^p (1 + F₂ⁿ/F₂^p)/2**
> using the parameterizations of the precisely measured structure function F₂^p [16]
> and F₂ⁿ/F₂^p ratio [17]."

**F₂ᵈ = F₂^p(1 + F₂ⁿ/F₂^p)/2 = (F₂^p + F₂ⁿ)/2.** That is the isoscalar *average*
of the free proton and neutron — the per-nucleon deuteron structure function, the
universal DIS convention. The per-*deuteron* object would be F₂^p + F₂ⁿ, with no
½. There is no reading of an arithmetic mean under which it is per deuteron.
(Note also that no nuclear correction — shadowing, EMC, Fermi smearing — is
applied at all: it is the bare free-nucleon average, which is only sensible as a
per-nucleon quantity.)

A_zz is a ratio of cross sections and carries no normalisation. Therefore **the
entire A-normalisation of the published b₁ᵈ is inherited from F₁ᵈ, and F₁ᵈ is per
nucleon. HERMES's published b₁ᵈ is PER NUCLEON.**

Refs [16] and [17] confirm the reading: [16] = H. Abramowicz *et al.*,
hep-ph/9712415 — the **ALLM97 parameterisation of σ_tot(γ\*p)**, a *proton*
parameterisation; [17] = NMC, P. Amaudruz *et al.*, NPB **371** (1992) 3 — the
*F₂ⁿ/F₂^p ratio*. Both are per-nucleon inputs. Nothing in the chain ever
multiplies by A.

### 2.3 Numerical confirmation, all six bins

Eq. (5) is invertible: HERMES publish A_zz^d and b₁ᵈ in the same Table II, so their
own numbers hand you the F₁ᵈ they used, F₁ᵈ = −b₁ᵈ/(3/2 · A_zz^d). I compared that
against F₁ᵈ rebuilt from their own stated recipe with ALLM97 (parameters taken
verbatim from Table 2 of hep-ph/9712415 — i.e. their ref [16], the same
parameterisation they used), Whitlow R₁₉₉₀ (their ref [18]) and NMC F₂ⁿ/F₂^p:

| ⟨x⟩ | ⟨Q²⟩ | A_zz [10⁻²] | b₁ᵈ [10⁻²] | F₁ᵈ implied | F₂^p (ALLM97) | R₁₉₉₀ | F₁ᵈ **per nucleon** | implied / per-N | implied / per-D |
|---|---|---|---|---|---|---|---|---|---|
| 0.012 | 0.51 | −1.06 | +11.20 | 7.044 | 0.2257 | 0.243 | 7.576 | **0.930** | 0.465 |
| 0.032 | 1.06 | −1.07 | +5.50 | 3.427 | 0.2995 | 0.383 | 3.380 | **1.014** | 0.507 |
| 0.063 | 1.65 | −1.32 | +3.82 | 1.929 | 0.3392 | 0.341 | 1.993 | **0.968** | 0.484 |
| 0.128 | 2.33 | −0.19 | +0.29 | 1.018 | 0.3583 | 0.262 | 1.091 | **0.933** | 0.466 |
| 0.248 | 3.11 | −0.39 | +0.29 | 0.496 | 0.3226 | 0.190 | 0.535 | **0.927** | 0.463 |
| 0.452 | 4.69 | +1.57 | −0.38 | 0.161 | 0.1881 | 0.129 | 0.179 | **0.904** | 0.452 |

**Mean ratio to per-nucleon F₁ᵈ: 0.946. Mean ratio to per-deuteron F₁ᵈ: 0.473.**

Across five decades of b₁ᵈ magnitude and a factor 9 in Q², the F₁ᵈ that HERMES
actually divided by is the per-nucleon one, to 5-10 %. It is not the per-deuteron
one, which it misses by a clean factor 2 in every bin. The residual ~5 % is the
size of the inputs *I* supplied and HERMES did not tabulate: R₁₉₉₀ enters as
(1 + R) so a ΔR = 0.1 is 7-8 %; F₂ⁿ/F₂^p was read off the NMC ratio rather than
re-fitted; and HERMES round A_zz to two decimals in units of 10⁻², which at
x = 0.128 (A_zz = −0.19) is already ±2.6 %. None of that is a factor 2.

Reproduce with (Python 3, no dependencies):

```python
import math
# ALLM97 (Abramowicz & Levy, hep-ph/9712415, Table 2) -- HERMES ref [16]
m02,mP2,mR2,Q02,L2 = 0.31985,49.457,0.15052,0.52544,0.06527
cP=(0.28067,0.22291,2.1979); aP=(-0.0808,-0.44812,1.1709); bP=(0.36292,1.8917,1.8439)
cR=(0.80107,0.97307,3.4942); aR=(0.58400,0.37888,2.6063); bR=(0.01147,3.7582,0.49338)
M=0.93827
def g(p,t): return p[0]+(p[0]-p[1])*(1.0/(1.0+t**p[2])-1.0)
def f(p,t): return p[0]+p[1]*t**p[2]
def F2p(x,Q2):
    t=math.log(math.log((Q2+Q02)/L2)/math.log(Q02/L2))
    W2=M*M+Q2*(1.0-x)/x
    xP=1.0/(1.0+(W2-M*M)/(Q2+mP2)); xR=1.0/(1.0+(W2-M*M)/(Q2+mR2))
    return Q2/(Q2+m02)*(g(cP,t)*xP**g(aP,t)*(1.0-x)**f(bP,t)
                       +f(cR,t)*xR**f(aR,t)*(1.0-x)**f(bR,t))
def R1990(x,Q2):  # Whitlow et al., PLB 250 (1990) 193 -- HERMES ref [18]
    th=1.0+12.0*(Q2/(Q2+1.0))*(0.125**2/(0.125**2+x*x))
    return 0.0635/math.log(Q2/0.04)*th + 0.5747/Q2 - 0.3534/(Q2*Q2+0.09)
rnp={0.012:1.00,0.032:0.99,0.063:0.97,0.128:0.92,0.248:0.83,0.452:0.68}  # NMC, ref [17]
rows=[(0.012,0.51,-1.06,11.20),(0.032,1.06,-1.07,5.50),(0.063,1.65,-1.32,3.82),
      (0.128,2.33,-0.19,0.29),(0.248,3.11,-0.39,0.29),(0.452,4.69,1.57,-0.38)]  # HERMES Table II
for x,Q2,azz,b1 in rows:
    azz*=1e-2; b1*=1e-2
    F1imp = -b1/(1.5*azz)                                   # invert HERMES Eq. (5)
    F1perN = (1+4*M*M*x*x/Q2)*F2p(x,Q2)*(1+rnp[x])/2/(2*x*(1+R1990(x,Q2)))
    print(f"{x:6.3f} {F1imp:8.3f} {F1perN:8.3f}  N:{F1imp/F1perN:5.3f}  D:{F1imp/(2*F1perN):5.3f}")
```

### 2.4 Two wrinkles worth recording, neither of which changes the verdict

1. **The sign.** Eq. (5) as rendered by `pdftotext` and by `pdftoppm` shows
   `b₁ᵈ = [gap] (3/2) A_zz^d F₁ᵈ`; the arXiv PDF drops the minus glyph (the same
   font damage eats most minus signs, `σ`, `ν` and `γ` throughout the paper). The
   minus is real: with it the back-out of §2.3 gives a positive F₁ᵈ in every bin,
   and CDKS quote HERMES's relation explicitly as *"In the HERMES analysis, b₁ was
   then extracted from A_zz using A_zz = −2b₁/(3F₁)"* (CDKS Eq. 48). Irrelevant to
   normalisation, but do not copy the unsigned form out of the PDF.
2. **HERMES's own table is off by 2 from HERMES's own extraction.** §2.1's QPM
   table, read with its gloss "in a … nucleus", is a per-deuteron definition;
   §2.2's Eq. (5) delivers a per-nucleon number. The published values follow
   Eq. (5). This is not a criticism of the measurement — the QPM table is there to
   define *which helicity combination* b₁ is, and the A-scaling is set at
   extraction — but it is exactly the trap that has caught the two theory camps in
   opposite directions, and it is why this item could not be closed by reading one
   paper.

---

## 3. Camp CDKS — definitional chain

**Paper.** W. Cosyn, Y.-B. Dong, S. Kumano, M. Sargsian, PRD **95** (2017) 074036,
[arXiv:1702.05337](https://arxiv.org/abs/1702.05337), on disk at
`/home/cpeng/Projects/polli/PolarizedLithiumSim/refs/1702.05337.pdf`.

| step | what it says | normalisation |
|---|---|---|
| Eq. (9) | W^A_{μν} = ∫d⁴p S(p) W_{μν}(p,q) | — |
| **Eq. (10)** | **S(p) = (1/A) Σ_i \|φ_i(p)\|² δ(…)** — explicit 1/A | per nucleon |
| Eq. (16) | b₁(x,Q²) = ∫(dy/y) δ_T f(y) F₁^N(x/y,Q²), δ_T f = f⁰ − (f⁺+f⁻)/2 | per nucleon |
| text under Eq. (16) | *"Here, the structure function b₁ is defined by the one per nucleon"* | **stated in words** |
| Eq. (17) | f^H(y) = ∫d³p y \|φ^H(p)\|² δ(y − (E−p_z)/M_N), with φ^H normalised to 1 ⇒ ∫f dy ≈ 1 | one nucleon |
| text after Eq. (18) | *"F₁^N is defined for the nucleon by the average of the proton and neutron structure functions: F₁^N = (F₁^p + F₁^n)/2"* | per nucleon |

Four independent statements, all agreeing, one of them literal. CDKS's b₁ is
**per nucleon**, on exactly the same footing as HERMES's published b₁ᵈ (§2), which
is why their Fig. 6 overlay is coherent — *"In Fig. 6, our xb₁ curves are shown at
Q² = 2.5 GeV² for comparison with the HERMES data … the magnitude of xb₁ is much
smaller than the HERMES data at x < 0.5"* (CDKS §IV). Their curve is genuinely
~400× below the data; that is their physics conclusion, not a normalisation
artefact. Their own Table I (theory-2 b₁ = 2.81×10⁻⁴ at x = 0.012 against HERMES's
1.12×10⁻¹) says the same.

CDKS Eqs. (45)-(49) — the block the task flagged as "where the per-nucleon-ness of
F₁ enters the extraction" — turn out **not** to carry the answer. Eq. (47) is a
ratio of cross sections; Eq. (49) (F_{UU,T} = 2F₁, F^{T_LL}_{U,T} = −2√(2/3) b₁)
and Eq. (48) (A_zz = −2b₁/3F₁) are all homogeneous in the b₁/F₁ ratio, so they are
invariant under a common 1/A on both. They tell you the two must share a
normalisation; they cannot tell you which one. Only HERMES's choice of F₁ᵈ (§2.2)
breaks the degeneracy — which is why this item needed the HERMES paper and not
CDKS's.

**Verdict: CDKS is PER NUCLEON. `b1_convolution()` must stop halving it.**
Certain — the paper says so in words, and every other element of the chain agrees.
This confirms, from the primary source rather than by inheritance, what
`docs/open_items/run_2026-09-02/design_D_b1_li6.md`:1590-1600 (design 9 Q1) already
concluded, and closes `OPEN_ITEMS_SOLUTIONS.md`:568 item 5 (Q5(a)) in the
affirmative.

---

## 4. Camp Miller — definitional chain

**Paper.** G. A. Miller, PRC **89** (2014) 045203,
[arXiv:1311.4561](https://arxiv.org/abs/1311.4561), on disk at
`/home/cpeng/Projects/polli/PolarizedLithiumSim/refs/1311.4561.pdf`.

### 4.1 The chain says per deuteron

| step | what it says | normalisation |
|---|---|---|
| **Eq. (1)** | b₁ = Σ_q e_q² [q↑⁰ − ½(q↑¹ + q↑^{−1})], with *"q↑^m (q↓^m) is the number density of quarks with spin up(down) along the z axis in a **target hadron** with helicity m"* | per deuteron; **no 1/A** |
| Eq. (3) | Δ^π q^{(m)}(x) = ∫_x (dy/y) q^π(x/y) f_π^{(m)}(y) | inherits f_π |
| Eq. (4) | q^π = (5/9)u_v^π + (10/9)ū^π + (2/9)s^π — already **charge-weighted and spin-summed** | — |
| **Eq. (5)** | f_π^{(m)}(y_A) = ∫(dξ⁻/2π) e^{−iy_A P_D⁺ξ⁻} ⟨D,m\|φ_π(ξ⁻)φ_π(0)\|D,m⟩_c — a light-cone correlator in the **normalised deuteron state**; ∫dy f_π = ⟨n_π⟩ **in the deuteron** | per deuteron |
| Eqs. (7)-(8) | the explicit evaluation, F_m(q) = ∫d³r ⟨D,m\|e^{−iq·r}σ₁·q σ₂·q\|D,m⟩ with the Eq. (10) deuteron wave function normalised to 1; the leading 3 is the isospin τ₁·τ₂ | per deuteron; **no ½** |
| **Eq. (6)** | b₁^π = ½[Δ^π q^{(0)} − Δ^π q^{(1)}] | the ½ is the **quark-spin average** — see below |
| Eq. (20) | b₁^π(x) = ½ ∫_x (dy/y) q^π(x/y) δf_π(y) | per deuteron |
| Eq. (24) | b₁^{6q} = ½ · (2) · (F₀(x) − F₁(x)) P_{6q}; his own text attributes the **(2)** to *"either state can be in the d-wave"* and is silent on the ½ | same ½ as Eq. (6) |

**Eq. (6)'s ½ is fully consumed by the quark-spin average and has no room left to
be a 1/A.** Derivation from Miller's own Eq. (1), which is the crux of the whole
item:

- The pion is spinless, so its quark distributions are spin-summed: for every
  target helicity m, q↑^{(m)} = q↓^{(m)} = ½ Δ^π q^{(m)}.
- Reflection symmetry gives q↑^{(1)} + q↑^{(−1)} = ½(Δ^π q^{(1)} + Δ^π q^{(−1)}) = Δ^π q^{(1)}.
- Substituting into Eq. (1):
  b₁^π = ½Δ^π q^{(0)} − ½ · Δ^π q^{(1)} = **½[Δ^π q^{(0)} − Δ^π q^{(1)}]** — exactly Eq. (6).

The ½ is accounted for to the last factor by Eq. (1) alone. If it were *also*
doing duty as 1/A = ½ for A = 2, Eq. (6) would be a factor 2 too small relative to
Eq. (1), which is the definition Eq. (6) is derived from. **So Miller's Eqs. (6),
(20), (24) and (26) are per deuteron, in the same sense as his Eq. (1).**

And Miller's Eq. (1) is, symbol for symbol, HERMES's own QPM definition (§2.1):
Miller's `q↑⁰ − ½(q↑¹ + q↑^{−1})` is algebraically identical to HERMES's
`½[2q↑⁰ − (q↑¹ + q↑^{−1})]`, and Miller's gloss *"in a target hadron with helicity
m"* is HERMES's *"in a hadron(nucleus) with helicity m"*. Miller evaluated the
shared **formal** definition — which is per nucleus (§2.1). HERMES's published
**numbers** are not that object; they are Eq. (5)'s, which is half of it (§2.2).

### 4.2 But his presentation says otherwise

Against that, three facts about how Miller *used* his result:

1. **Table I is HERMES's Table II verbatim.** Miller's columns
   `b₁ ± δb₁^stat ± δb₁^sys [10⁻²]` are 11.20/5.51/2.77, 5.50/2.53/1.84,
   3.82/1.11/0.60, 0.29/0.53/0.44, 0.29/0.28/0.24, −0.38/0.16/0.03 — an exact,
   unrescaled transcription of HERMES Table II. His own b₁^π and b₁^{6q} sit in
   adjacent columns **in the same 10⁻² units**, and he writes *"The differences
   obtained by using different structure functions are generally not larger than
   the experimental error bars."*
2. **Fig. 5 overlays the six HERMES points on his curve** (confirmed by rendering
   p. 10 at 200 dpi — the points and error bars are visibly HERMES's), and the text
   says *"This is also shown in Fig. 5, where very good agreement between data and
   our model can be observed."*
3. **He tuned to a HERMES number at face value.** *"We choose P_{6q} = 0.0015 to
   reproduce the central value"* of *"the measured value, b₁ = −3.8 ± 0.16 × 10⁻³"*
   at x = 0.452. Incidentally he mis-transcribed the uncertainty by a decade
   (HERMES Table II gives 0.16 × 10⁻², i.e. ±1.6 × 10⁻³, not ±0.16 × 10⁻³) — a
   small slip, but it is the slip of someone reading numbers straight off a
   published table without interrogating what they are normalised to.

The digitized table in-tree agrees with (2): `tables::kB1Miller` interpolated at
x = 0.012 is **0.11429**, against HERMES's 0.1120 — Miller's plotted curve passes
through the HERMES point.

### 4.3 The two readings, and which one the repository should hold

The paper cannot be both. Exactly one of these is true:

* **(A) The chain is right and the comparison is off by 2.** Miller faithfully
  evaluated Eq. (1)'s per-deuteron b₁ and compared it, unaware, to HERMES's
  per-nucleon numbers — an easy mistake to make, because HERMES's *own definition
  table* is the per-deuteron one (§2.4 wrinkle 2). Then the digitized column is per
  deuteron and **the 0.5 is correct**.
* **(B) The comparison is right and there is an unlocated 1/A in the chain.** Then
  the digitized column is per nucleon and the 0.5 is wrong.

**(A) is much the stronger.** (B) requires an extra factor ½ that §4.1 shows is not
in Eqs. (1), (5), (6), (7), (8) or (20) — every ½ in the chain is spoken for, and
the deuteron matrix element in Eq. (5) is manifestly a per-deuteron pion number
density. (A) requires only that Miller trusted HERMES's definition table over
HERMES's extraction recipe, which is what the definition table invites. Miller's
Fig. 5 caption also states the curve *is* Eq. (20) + Eq. (26) — so the figure is
the chain, and the chain is per deuteron.

Under (A), Miller's physics survives the factor; only his adjectives move. His
b₁^π = 10.5 × 10⁻² at x = 0.012 becomes 5.25 × 10⁻² per nucleon against HERMES's
(11.20 ± 5.51 ± 2.77) × 10⁻², which is 0.97 σ low — still consistent, but "very
good agreement" overstates it. And P_{6q} = 0.0015, tuned to the x = 0.452 point,
would become 0.00075 — a model parameter he already calls *"very, very small … can
not be ruled out by any observations"*, so halving it costs him nothing.

**Verdict: Miller is PER DEUTERON; keep the ×0.5 at `sf.cpp:350-354`.** Likely, not
certain (§7).

---

## 5. Why "both camps sit on the same HERMES points" is not the contradiction it looked like

`OPEN_ITEMS_SOLUTIONS.md` framed this as *"Both cannot be true if both sit on the
same HERMES points."* They can, if one of them is wrong about the points — and that
is what happened. CDKS overlay a per-nucleon curve on per-nucleon data and say
correctly that they undershoot by ~400×. Miller overlays a per-deuteron curve on
per-nucleon data and reports agreement that is a factor 2 too flattering. The
overlay is evidence about *intent*, and intent is not normalisation. The chain is.

---

## 6. Code sites: exactly what must stop sharing `constants.hpp:86`

The two camps now need different factors, so they must stop sharing one applied
constant. The minimal fix that respects "no physics number is defined twice"
(`docs/CONVENTIONS.md`) is **not** to add a second constant — 0.5 is a single,
correct, well-named quantity (the deuteron→nucleon conversion, 1/A at A = 2) — but
to **stop applying it at the site that does not need it**, and to fix the comment
that asserts both need it.

### 6.1 `src/core/sf.cpp:356-361` — `b1_convolution()` / `CdksB1` — **DROP the factor**

```cpp
  const double xb1 = tables::kB1CdksQ2p5().interp("xb1_theory1_sum", xs);
  return B1_PER_DEUTERON_TO_PER_NUCLEON * xb1 / xs;   // <-- the 0.5 is wrong
```

should return `xb1 / xs`. The digitized CDKS column is already per nucleon (§3);
halving it makes `CdksB1` a factor 2 low. Blast radius, measured:

* **`tests/test_sf.cpp:175-176`** — the only rtol-pinned CDKS values. Both double:
  `b1_convolution(0.3, 2.5, 0.0)`: −2.816198975086724e-04 → −5.632397950173448e-04;
  `b1_convolution(0.05, 2.5, 0.0)`: +5.8517162680014353e-05 → +1.1703432536002871e-04.
* The `|b1_convolution| < 1e-3` bounds at `tests/test_sf.cpp:177-178` **survive**
  doubling — I evaluated them: x = 0.25 → 5.328e-04, x = 0.3 → 5.632e-04, both
  still under 1e-3. No assertion has to be relaxed.
* `tests/test_sf.cpp:405` (`cdks.b1_func()(…) == b1_convolution(…)`) is relative
  and unaffected.
* **`python/tests/`: no numeric pin on `CdksB1` exists** (`grep -n -i cdks
  python/tests/test_b1_model.py` shows only structural, CLI and validation tests —
  registry membership, `--x-max` behaviour, band handling). Nothing there moves.
* **`validation/reference/*.json` is untouched**: the only b₁ reference,
  `b1_default_li6.json`, records `"b1_model": "miller"`. **The rtol 1e-12 gate does
  not move.**
* `docs`-level: `validation/b1_li6_table.py:61` builds `lg.Li6B1(lg.CdksB1())` for
  a comparison table; that column doubles and the table must be regenerated.
* `src/core/b1_nuclear.cpp:573` already routes around the bug via
  `cdks_b1_raw_per_nucleon()`, so the ⁶Li convolution backend is **unaffected** —
  and once `b1_convolution()` is fixed, that accessor's reason for existing
  becomes documentation rather than a workaround. Its doc comment
  (`include/lipolgen/b1_nuclear.hpp:115-125`) and the binding docstring
  (`python/bindings.cpp:1247-1251`), which both say "`CdksB1` halves it a second
  time", then need updating too.

This *is* a behaviour change for the opt-in `--b1-model cdks` path. Under the
"new physics is OPT-IN / defaults stay bit-for-bit" rule it is a supervisor
decision whether to flip it directly (the default `miller` path does not move, so
no reference JSON regenerates) or to stage it behind a flag. The physics says it
must be flipped; the mechanism is not mine to choose.

### 6.2 `src/core/sf.cpp:350-354` — `toy_b1()` / `MillerB1` — **KEEP the factor**

No change. This is the default b₁ path (`Li6B1(MillerB1)`, `cli.py:164`,
`validation/reference/b1_default_li6.json`), so keeping it means **the rtol 1e-12
gate is untouched by this item**. What must change is the *justification*: the
comment must stop pointing at Fig. 5's axis (§1) and point at Eqs. (1), (5) and (6)
instead, and must record that the factor rests on reading (A) of §4.3.

### 6.3 `include/lipolgen/constants.hpp:86` — **fix the comment, keep the number**

`B1_PER_DEUTERON_TO_PER_NUCLEON = 0.5` stays, one home, unchanged value. The two
lines above it must stop claiming it applies to both curves. Something like: *the
Miller table (`kB1Miller`) is per deuteron — Miller Eq. (1)'s densities are "in a
target hadron", and Eq. (6)'s ½ is the quark-spin average, not 1/A — so it is
halved on the way out; the CDKS column (`kB1CdksQ2p5`) is already per nucleon by
CDKS Eq. (10)'s explicit 1/A and their text under Eq. (16), and is NOT halved.*

### 6.4 `src/core/sf.cpp:363-371` — `close_kumano_integral()` — **flag, do not fix here**

It integrates the **raw** columns for both camps (no factor either way), so it is
normalisation-blind by construction. But that means the two numbers it returns are
on *different* scales — a per-deuteron integral for `cdks == false` and a
per-nucleon one for `cdks == true` — and they are not comparable to each other.
Whether that matters depends on what consumes it; not in this item's scope, but it
should not be left unrecorded.

---

## 7. `tests/test_sf.cpp:159-164` — the premise holds; one comment is wrong

```cpp
  // The tables hold the published PER-DEUTERON b1 and the accessor halves it,
  // because every consumer pairs b1 with a per-nucleon F1.  At x = 0.012 the
  // digitized total is 0.114 per deuteron, against HERMES's 0.112 +- 0.055.
  CHECK_CLOSE(2.0 * toy_b1(0.012, 2.5, 0.0), 0.11429317074113018, kRtol);
  CHECK_CLOSE_AT(2.0 * toy_b1(0.012, 2.5, 0.0), 0.114, 0.0, 0.01);
```

**No assertion has to move.** Both are digitization checks on `kB1Miller` — they
assert that reading the table at x = 0.012 gives 0.114, which is a fact about the
digitization and is true whatever the table is per. Sentence one of the comment
("the tables hold the published PER-DEUTERON b1") is right for the Miller table
under §4.3(A) — though "the tables", plural, is wrong: it is false for
`kB1CdksQ2p5` (§3).

**Sentence two is wrong and must be corrected.** *"the digitized total is 0.114 per
deuteron, against HERMES's 0.112 ± 0.055"* compares a per-deuteron number to a
per-nucleon datum — it silently reproduces Miller's own factor-2 slip (§4.2) inside
our test suite, and reads as confirmation that the halving is *unnecessary*, which
is the opposite of what it is there to defend. The honest comparison is
**0.114/2 = 0.057 per nucleon against HERMES's 0.112 ± 0.055 ± 0.028 per nucleon,
i.e. 0.97 σ low** — consistent, but not the "sits on the data point" agreement the
comment implies. Whoever fixes it should also note that the *coincidence* the
comment is pointing at (the digitized curve passing through the HERMES point) is
evidence for reading (B), not (A), and is precisely the thing §4.3 decides against.

---

## 8. Is this decidable offline?

**Two thirds yes, one third no.**

* **HERMES's normalisation — DECIDED, offline, certain.** The paper states the
  recipe (§2.2) and its own Table II inverts to confirm it in all six bins (§2.3).
  Nothing further is needed. Recommend adding hep-ex/0506018 to
  `PolarizedLithiumSim/refs/`.
* **CDKS — DECIDED, offline, certain.** *"the structure function b₁ is defined by
  the one per nucleon"*, plus Eq. (10)'s 1/A, plus F₁^N = (F₁^p+F₁^n)/2, plus a
  coherent Fig. 6. `OPEN_ITEMS_SOLUTIONS.md`:568 item 5 can close: **yes, `CdksB1`
  must stop halving.**
* **Miller — NOT fully decidable offline. It needs the author, or a day of
  numerics.** The paper is self-inconsistent by exactly the factor in question:
  the derivation (§4.1) says per deuteron, the presentation (§4.2) says per
  nucleon, and no third party adjudicates — CDKS cite Miller (their ref [7]) once,
  in the introduction, with no comment on his normalisation. My verdict, reading
  (A), is a derivation from his equations plus a judgement about which of two
  mistakes is likelier. It is not a measurement of his code.

  **What would settle it, in order of cost.** (i) **Ask Miller.** One sentence from
  him closes it. (ii) **Reproduce Eq. (20) numerically**: δf_π(y) from Eqs. (17)-(19)
  with an AV18 u, w, convolved with the GRS/SMRS pion PDF he used (his ref [29],
  "model 1"), evaluated at x = 0.012 and compared against his tabulated
  10.5 × 10⁻². If it lands on 10.5, the chain is per nucleon and reading (B) wins;
  if it lands on 5.25, reading (A) wins and the current 0.5 is confirmed. This is a
  real piece of work — it needs a pion PDF set the repository does not have (the
  same gap as `OPEN_ITEMS_SOLUTIONS.md`:562 item 2, which wants MSTW2008 LO or
  Kumano's tabulated curves) — and it is the *only* offline route. I did not
  attempt a shortcut estimate: at x = 0.012 the convolution is dominated by the
  small-argument sea tail of q^π, where a factor-few uncertainty is unavoidable
  without the actual parameterisation, and a factor-few answer cannot discriminate
  a factor 2. **A number I did not run is not a result, so I am not reporting one.**

* **The safe posture in the meantime.** Reading (A) *is* the status quo, so keeping
  the 0.5 for Miller changes nothing, moves no default, and touches no rtol 1e-12
  gate. The item stays open with a much sharper question than it had — no longer
  "check Miller's normalisation" but "does Eq. (20) evaluate to 10.5 × 10⁻² or
  5.25 × 10⁻² at x = 0.012?" — and with the CDKS half closed.

---

## 9. What I measured vs. what I assumed

**Measured.**
* `build/lipolgen_tests` before and after: 363 / 17203863 / 1 skipped, SUCCESS. No
  code touched.
* `kB1Miller` interpolated at x = 0.012 = 0.11429317074113018 (from the table
  arrays in `src/core/sf_tables.cpp`, reproducing the pinned test value); table
  neighbours x = 0.00999 → 0.129359, x = 0.01297 → 0.107025.
* `lg.b1_convolution` at x = 0.05, 0.1, 0.25, 0.3 (Q² = 2.5), and that all four
  stay under the 1e-3 bound when doubled.
* `lg.B1_PER_DEUTERON_TO_PER_NUCLEON == 0.5` through the Python bindings.
* The §2.3 six-bin inversion of HERMES Eq. (5), with ALLM97 parameters read out of
  hep-ph/9712415 Table 2.
* The absence of "per nucleon" / "1/A" anywhere in Miller, by grep over the full
  text dump.
* Miller p. 10 and HERMES p. 4 rendered at 200 dpi and read as images, because the
  arXiv PDFs' embedded fonts are damaged and `pdftotext` alone is not trustworthy
  for signs and Greek letters in either.

**Assumed.**
* R₁₉₉₀'s functional form and coefficients (0.0635, 0.5747, −0.3534 with
  Θ = 1 + 12·(Q²/(Q²+1))·(0.125²/(0.125²+x²)), Λ² = 0.04) are from memory, not from
  Whitlow's paper, which I did not fetch. F₂ⁿ/F₂^p was read off the NMC ratio
  approximately, not refitted. Both enter §2.3 at the 5-10 % level and neither can
  manufacture or hide a factor 2 — but the §2.3 residual of 0.946 rather than 1.000
  is plausibly them and should not be read as a physics statement.
* That Miller's deuteron wave function (his Eq. 10) is normalised
  ∫dr(u² + w²) = 1, the universal convention, which he does not state.
* Reading (A) over reading (B) in §4.3 — a judgement, not a measurement. Everything
  in §6.2 and §7 that keeps the Miller 0.5 rests on it.
