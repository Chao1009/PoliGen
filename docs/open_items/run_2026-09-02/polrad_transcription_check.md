# POLRAD 2.0 transcription check — design C §6 step 0 / open question Q2

**Status: step 0 is DISCHARGED.** Both primary sources were obtained: the
paper's **LaTeX source** (so no `pdftotext` fraction ambiguity survives) and the
**POLRAD 2.0 FORTRAN** (CPC catalogue **`ADGH_v1_0`**, *not* `ADXQ` as the
design says at lines 1964 and 2096; `ADXQ` appears nowhere else under `docs/`).
All five questions are settled: four confirmed, one corrected. Two further errors were found while checking and are
recorded in §7 — one in the design (`ℑ^el_6`), one in POLRAD itself (the
approx-path carbon `Z²`) — plus one **new blocking flag** (§8, the overall sign
of Eq. (38)) that the implementer must resolve before `rc_tail` is trusted.

---

## 0. Sources obtained, and how to get them again

| source | how | check |
|---|---|---|
| **Paper LaTeX** | `curl -sL -o src.tar.gz https://arxiv.org/e-print/hep-ph/9706516` → `tar xzf` → **`polrad2t.tex`** (4201 lines) | the fraction macro is `\def\ot#1#2{\textstyle \frac{#1}{#2}}` (line 9), so `\ot 3 4` **is** 3/4 and `\ot 4 3` **is** 4/3 — no ordering ambiguity |
| **FORTRAN** | Mendeley Data, dataset `37vgvzgr2w`: `curl -sL "https://data.mendeley.com/public-files/datasets/37vgvzgr2w/files/50bc36c4-d146-4487-b697-218850c0a057/file_downloaded" -o adgh_v1_0.gz` (97 335 B, `sha256 = 83703668bacfbffe073a63c75374be261b4b36dab8ea8de951762120fba1ea7b`) → `gunzip` → one PATCHY archive `adgh`, 11 593 lines | `+deck,elu` @ 8878, `+deck,elp` @ 8915, `+deck,elq` @ 8972, `+deck,apptai` @ 8559, `+deck,strf` @ 4147, `+deck,ffdeu` @ 5479, `+deck,ffco` @ 5576. Line numbers below are lines of `adgh`. |

`pdftotext -layout` on the PDF renders `\frac{8}{9}` as `98` but `\frac{2}{3}` as
`23` — the reading order of a stacked fraction flips with the glyph widths. **The
PDF is not usable evidence for any fraction in this paper.** The `.tex` is.

**One caveat on the FORTRAN, stated up front.** The Eq. (37)–(39) closed forms
live in patch `polrad_add` (decks `apptai`, `elu`, `elp`, `elq`) under the
PATCHY switch `approx`, and the shipped steering file `polrad20.cra` has **both
commented out**:

> `c+use,polrad_add,if=polrad.` (line 12) and `c+use,approx,if=polrad.` (line 17)

so the **default build of POLRAD 2.0 does not compile the t-peak closed forms at
all** — it uses the exact Eq. (18) + Appendix-B path (`+use,exact,if=polrad.`,
line 16). The `approx` decks are therefore an *authoritative transcription by the
authors* but not a *routinely exercised* one. Where the two disagree (§7.2) the
exact path wins. Design C v0 is built entirely on the closed forms, i.e. on
POLRAD's own less-travelled road; §8 is the price of that.

---

## 1. Does `σ_u^C` carry `Y₋`?  → **CORRECTED: it is `Y₊`.**

**Design (line 387–389, §1.4.2 flag 1):**

> *"1. `σ_u^C` is printed with **`Y₋`**, not `Y₊`, unlike every other unpolarised
> entry in Eq. (38). Transcribe it as printed, then check it against the FORTRAN;
> if it is a typo the carbon gate below must use `Y₊`."*

and in the design's Eq. (38) block:

> ```
>   sigma_u^C  =  (alpha^3/S) Z^2 Y_-   INT dEta_A/Eta_A  Xt F^2      (38)
> ```

**Source, `polrad2t.tex` line 928 — the design's reading of the paper is right:**

```tex
\sigma^C_u&=&
\od{\alpha^3}{S}Z^2Y_-
\int\limits_{\eta_{min}}^{\infty}\od{d\eta_A}{\eta_A}
{\tilde X}F^2
```

**But the FORTRAN says `Y₊`.** `Y±` is not inside `elu` at all; `elu` returns the
`η_A` integrand only, and `apptai` applies the `Y` factor once, by *polarisation
degree*, for **every** target:

```fortran
c  adgh:8614-8648   (+deck,apptai)
      yy1=(1.+(1.-ys)**2)/(1.-ys)                      ! = Y_+
      ...
        siau=ter*alfa**3/s*barn*yy1*relu               ! unpolarised  -> Y_+
      ...
        yy2=ys*(2.-ys)/(1.-ys)                         ! = Y_-
        siap=ter*alfa**3/s*barn*yy2*relp               ! vector       -> Y_-
      ...
        siaq=ter*alfa**3/s*barn*yy1*relq               ! tensor       -> Y_+
      ...
      apptai=un*siau + pl*pn*siap + qn/6.*siaq
```

and the carbon branch of `elu` (adgh:8903-8905) supplies only the integrand:

```fortran
+self,if=targ_c,targ_o.
         call ffco(t,ff)
         elu=xxt*ff**2/eta
```

The second-order deck `elual2` (adgh:9201 ff.) corroborates: it defines both
`ypl` and `ymi` and then uses **`ypl` only** in the unpolarised kernels
`rr1u`, `rr2u`, with `ff1u=ff**2, ff2u=0d0` for carbon (adgh:9267-9270).

**Independent internal check (no code needed).** Eq. (38)'s own proton and
deuteron entries obey one master formula,

> `σ_u^A = (α³/S) Y₊ ∫dη_A/η_A [ X̃ ℑ₂^el|_{Q_N=0} − ℑ₁^el|_{Q_N=0}/η_A ]`

— for the proton, Eq. (40) gives `F₁+F₂ = G_M` and
`F₁²+η F₂² = (G_E²+ηG_M²)/(1+η)`, so `σ_u^p` is exactly that with Eq. (A.5); for
the deuteron `ℑ₁|₀/η = (2/3)(1+η)F_m²` and `ℑ₂|₀ = F_c²+(2/3)ηF_m²+(8/9)η²F_q²`
reproduce `σ_u^d` term for term. A spin-0 nucleus has `ℑ₁ = 0`,
`ℑ₂ = Z²F²`, which gives `σ_u^C = (α³/S) Z² **Y₊** ∫dη_A/η_A X̃F²`. Physically the
same: `Y₋ = y(2−y)/(1−y)` is POLRAD's *vector-polarised* combination (it is what
Eq. (38) puts on `σ_p^p` and `σ_p^d`) and it **vanishes as `y → 0`**, which an
unpolarised cross section cannot.

> **Verdict: CORRECTED.** The `Y₋` in the paper's `σ_u^C` is a typo.
> **Replacement:**
> ```
>   sigma_u^C  =  (alpha^3 / S) * Z^2 * Y_+ * INT_{eta_min}^{eta_max} (d eta_A / eta_A) * Xt * F^2
>   Y_+ = (1 + (1-y)^2) / (1 - y)
> ```
> Design line 387–389 should be rewritten from "transcribe as printed, then
> check" to "the paper prints `Y₋`; POLRAD's own `apptai` applies `Y₊`; use `Y₊`",
> and the carbon gate **T8(b)** must be written against `Y₊`.

---

## 2. Is the last term of `X₁` linear in `η_A`?  → **CONFIRMED.**

**Design (line 390–393, flag 2):**

> *"2. `X₁ = x² + 4xη_A − 4η_A` — the last term is **linear** in `η_A`. The check
> that fixes it: `X₁(η_min) = x² + 4η_min(x−1) = 0` identically, which is what
> makes `X̃` vanish at the lower limit and the `∫dη_A/η_A` converge there. A
> transcription with `4η_A²` fails that identity."*

**Source, `polrad2t.tex` line 936 (Eq. (39)):**

```tex
X_1=x^2 + 4 x \eta_A  - 4 \eta_A, \quad {\tilde X}=\od{X_1}{2\eta_Ax^2}, \quad
Y_{\pm}=\od{1\pm (1-y)^2}{1-y},
\eta_A=\od{t}{4M^2_A}, \quad \eta_{min}=\od{x^2}{4(1-x)},
```

**FORTRAN, three independent copies, all identical** (adgh:8887 in `elu`,
adgh:8930 in `elp`, adgh:8980 in `elq`; also 9213 in `elual2`, 9292 in `elpal2`):

```fortran
      xx1=xa**2 + 4.*xa*eta - 4.*eta
      xxt=xx1/(2.*eta*xa**2)
```

The design's own identity check is also correct as written:
`X₁(η_min) = x² + 4·[x²/(4(1−x))]·(x−1) = x² − x² = 0`.

> **Verdict: CONFIRMED.** No replacement. `X₁ = x_A² + 4 x_A η_A − 4 η_A`,
> `X̃ = X₁/(2 η_A x_A²)`.

---

## 3. `(3/4)x²` vs `(4/3)η_A F_q` in `σ_q^d`  → **CONFIRMED, exactly as the design has them.**

**Design (line 394–396, flag 3):**

> *"3. The two fractions inside `σ_q^d` are `(3/4)x²` (inside the `X̃` bracket) and
> `(4/3)η_A F_q` (in the last line). They are *not* the same fraction and
> `pdftotext` renders both as bare digit pairs."*

**Source, `polrad2t.tex` lines 915–926 (Eq. (38), `σ^d_q`):**

```tex
\sigma ^d_q &=&
\od{\alpha^3}{S}  Y_+
\int\limits_{\eta_{min}}^{\infty}\od{d\eta_A}{\eta_A}
\biggl\{
\Bigl(1+\eta_A+\Bigl(\ot 3 4 x^2 - \eta_A\Bigr) {\tilde X}\Bigr) F_m^2
-\od {{\tilde X} X_1}{(1+\eta_A)}  F_q  \bigl(3 F_c+3 \eta_A F_m+\eta_A F_q\bigr)
-2 \eta_A {\tilde X}  F_q  \Bigl(4 F_c-3 x F_m+\ot 4 3 \eta_A F_q\Bigr)
\biggr\}
```

with `\ot 3 4` = ¾ and `\ot 4 3` = 4/3 by the macro on line 9.
(For the record, `pdftotext` prints these two as `43 x2` and `34 ηA Fq` —
i.e. **backwards**. Anyone who transcribed from the PDF would have swapped them.)

**FORTRAN, `+deck,elq` (adgh:8985-8987), decimal and integer, unambiguous:**

```fortran
         elq=((1.+eta+(.75*xa**2-eta)*xxt)*fm**2
     .        -xxt*xx1/(1.+eta)*fq*(3.*fc+3.*eta*fm+eta*fq)
     .        -2.*eta*xxt*fq*(4.*fc-3.*xa*fm+4./3.*eta*fq))/eta
```

`.75*xa**2` ⇒ **(3/4)x_A²**; `4./3.*eta*fq` ⇒ **(4/3)η_A F_q**. The remaining
coefficients (`3F_c + 3η_A F_m + η_A F_q`, `4F_c − 3x_A F_m`) match the design's
block character for character. The companion `σ_u^d` is equally confirmed
(adgh:8893-8894: `8./9.*fq**2*eta**2`, `2./3.*fm**2*eta`, `-2./3.*(1.+eta)*fm**2`).

> **Verdict: CONFIRMED.** No replacement. Design §1.4.2's `σ_q^d` block is
> correct as printed. Keep flag 3's warning about `pdftotext` — it is a real trap
> and it points the wrong way.

---

## 4. Eq. (18)'s `1/A²` right / `1/A` left, and `x_A = x/A` (Eq. (21))  → **CONFIRMED.**

**Design (line 558–563, §1.4.6):**

> *"**The prefactor is `1/A²` on the RIGHT and `1/A` on the LEFT** — not "`1/A` on
> both sides" as the first draft printed and as Q2 described it. Compare Eq. (19)
> (the QRT) and Eq. (20), which carry `α³y/A`, with no square."*

**Design (line 331, §1.4.1 nuclear-map row):**

> *"The source is POLRAD **Eq. (21)**, whose elastic exponentiation factor reads
> `(y²(1−x/A)²/(1−xy/A))^{t_r}` against the quasi-elastic `(y²(1−x)²/(1−xy))^{t_r}`
> — so the elastic formulas' `x` is `x_A = x/A`."*

**Source, `polrad2t.tex` lines 578–586 (Eq. (18)):**

```tex
\sigma_1^{el}= {1\over A}{d^2\sigma ^{el}\over dx_Ady}=
- {\alpha ^{3}y\over A^2}\int\limits^{\tau_{Amax}}_{\tau_{Amin}}d\tau_A
 \sum^{8}_{i=1}\sum^{k_{i}}_{j=1} \theta _{ij}(\tau_A){
 2M^{2}_A R^{j-2}_{el}\over (1+\tau_A)(Q^2+R_{el}\ta_A)^{2}}\Im^{el}_{i}(R_{el},\tau_A).
```

`{1\over A}` on the left, `{\alpha^3 y \over A^2}` on the right. **Verbatim.**
Eqs. (19) and (20) (lines 596 and 620) both read `- {\alpha^{3}y\over A}` — the
contrast the design relies on is real.

**Source, `polrad2t.tex` lines 631–636 (Eq. (21)):**

```tex
 \si^{el}_1\rightarrow \left({y^2(1-x/A)^2\over 1-xy/A}\right)^{t_r}\si^{el}_1,
\qquad
 \si^{q}_1\rightarrow \left({y^2(1-x)^2\over 1-xy}\right)^{t_r}\si^{q}_1,
```

**FORTRAN (adgh:479-480), which is Eq. (21) with `tara` = A:**

```fortran
      extai2=((sx-y/tara)**2/s/(s-y/tara))**tr     ! elastic
      extai3=((sx-y)**2/s/(s-y))**tr               ! quasi-elastic
```

With POLRAD's `y ≡ Q²`, `sx ≡ S_x`, `s ≡ S`: `y/(A·sx) = x/A` and
`y/(A·s) = x y_Bj / A`, and `sx²/s² = y_Bj²`, so `extai2 =
(y²(1−x/A)²/(1−xy/A))^{t_r}` exactly. Confirmed both ways.

> **Verdict: CONFIRMED.** No replacement. Design §1.4.6's Eq. (18) block and
> §1.4.1's `x_A = x/A` row are right.

**Three refinements the FORTRAN adds, which the design should absorb:**

* **POLRAD's own `x_A` is `x·m_p/M_A`, not `x/A`.** `conkin` (adgh:747) sets
  `s = snuc*amp/amh`, i.e. `S_A = S·M_A/m_p`, and the tail code uses
  `xa = y/sx = Q²/S_xA` (adgh:8886). For ⁶Li `M_A/(A m_p) = 1.00502`, so
  `x_A(POLRAD) = x/A × 0.99500` — a **0.5 % offset** from the design's `x/A`.
  > **REVISED 2026-09-03.** POLRAD's map is exact for POLRAD's variables: a
  > **target at rest**, with the per-nucleon `x` defined against the FREE
  > NUCLEON mass. LiPolGen's are not those. `beams.cpp` carries the
  > per-nucleon momentum `p_u = p_A/A`, so `s = 4E_e p_u` and
  > `S_A = 2k₁·p_A = A·s` and `x_A = Q²/(y S_A) = x/A` — **both exact**, for
  > *this* library. `src/core/rc.cpp` therefore uses the design's `x_A = x/A`,
  > `S_A = A·s` and keeps POLRAD's form only as a documented fixed-target
  > alternative. The difference, `A m_p/M_A = 1.005`, is not negligible on the
  > tail because it is steep in `x_A` through `t_min`: measured **0.26 % at
  > `x = 0.01`, 1.4 % at 0.1, 6.7 % at 0.3** on `σ^el_U`.
* ~~**The per-nucleon reduction in the code is `(m_p/M_A) × (1/A)`, not
  `1/A`** … Net difference for ⁶Li: **0.5 %**~~
  > **THIS PARAGRAPH WAS WRONG AND IT CAUSED A FACTOR-6 BUG. REWRITTEN
  > 2026-09-03.** The two factors are real — `apptai` carries
  > `ter = amh/amp = m_p/M_A` for `ita=2` (adgh:8607-8608) **and** the main
  > program divides again by `tara` (adgh:489, 497:
  > `(tai(2)*extai2+tai(3)*extai3)/tara`), against a Born that is itself per
  > nucleon (`f2sfun`, adgh:4330-4339) — but the design's §1.4.5 single `1/A`
  > is **not the same object**: it is the **Jacobian `dx_A/dx`**, not a
  > per-nucleon division. So POLRAD applies **two** factors where the design
  > applies **one**, and the difference is
  >
  > > **`A · m_p/M_A = 5.97`, NOT 0.5 %.**
  >
  > The old "net difference 0.5 %" sentence compared `(m_p/M_A)(1/A)` against
  > the design's `1/A` as if both were per-nucleon divisions. `rc.cpp`'s first
  > implementation cited exactly that sentence to justify applying **one**
  > factor (`per_nucleon = m_p/M_A`), which made `σ^el_U` and `σ^el_T`
  > **6.00× too large** for ⁶Li (2× for the deuteron) while the QRT, which had
  > always carried its `/A`, was right — so the ERT : QRT ratio was internally
  > inconsistent by a factor `A`. §8's Weizsäcker–Williams × Compton identity
  > settles it independently: **Eq. (38) is the WHOLE-NUCLEUS
  > `d²σ/(dx_A dy)`**, and the paper's `σ₁^el = (1/A) d²σ/dx_A dy` notation
  > (polrad2t.tex:579) is loose. The correct per-nucleon factor is
  > `(1/A)·(dx_A/dx)`, i.e. **`1/A²`** with this library's map — which is what
  > `src/core/rc.cpp` now applies, and what `tests/test_rc.cpp` T8(a') pins at
  > 1e-12 together with a gate that the ERT and QRT factors agree.
* **The upper limit is NOT `∞`.** The paper writes `∫_{η_min}^{∞}`; the code
  integrates `eta1 → eta2` with (adgh:8610-8613, comment at 8610)
  ```fortran
  c     eta1=xs**2/(4.*(1.-xs))                                  ! the paper's eta_min, commented OUT
        um=sx-y
        eta1=(um*(sx-sqly)+2.*amp2*y)/(8.*amp2*(um+amp2))        ! exact kinematic lower limit
        eta2=sx/4./amp2                                          ! exact kinematic UPPER limit
  ```
  i.e. POLRAD uses the **exact** limits and keeps the ultrarelativistic
  `η_min = x²/(4(1−x))` only as a commented-out reference. Measured spread
  (⁶Li, `M_A = 5.601518702`): `η_min(exact)/η_min(u.r.) =` 0.9996 at
  `x = 0.01, y = 0.8`; 0.983 at `x = 0.1, y = 0.5`; **0.921 at `x = 0.3, y = 0.3`.**
  Since the integrand is `∝ X̃`, which **vanishes at the u.r. `η_min` and is
  non-zero at the exact one**, the choice of lower limit is not a rounding
  detail at high `x`. Design §1.4.2's `INT eta_m ... INF` should become the
  exact pair, with `η_min(u.r.)` kept as the T8 cross-check. (`t_min = 4M_A²η_min
  = 0.00880 GeV²` vs `(x m_N)² = 0.00880 GeV²` at `x = 0.1` — the design's
  `t_min ≈ (x M_N)²` claim is confirmed to 3 figures.)

---

## 5. `Q_N/6` in Eq. (37), and the `Q_N = 0` Rosenbluth limit of Eq. (A.4)  → **CONFIRMED.**

**Design (line 342, §1.4.2):**

> *"**Note the explicit `Q_N/6`** — the tensor tail's normalisation is fixed here,
> not chosen."*

**Design (§1.4.4):** *"the tail carries it as **`Q_N/6`** (Eq. (37)), not `Q_N` —
that factor is fixed by POLRAD, not chosen here, and dropping it inflates the
tensor tail sixfold with no unpolarised test seeing it."*

**Design (line 618–620, §1.4.6):**

> *"`U_1` and `U_2` are the **standard Rosenbluth `B/2` and `A`** of a spin-1
> target (`A = G_C² + (8/9)η²G_Q² + (2/3)ηG_M²`, `B = (4/3)η(1+η)G_M²`), which
> **fixes the form-factor normalisation convention with no ambiguity**"*

**Design (§2.1, line 866–868):** `F_c(0) = Z = 3`, `F_m(0) = (M_A/m_p)·μ_A/μ_N`,
`F_q(0) = M_A²·Q_A`.

### 5a. `Q_N/6` — confirmed four times over

`polrad2t.tex` lines 877–882 (Eq. (37)):

```tex
\sigma ^{el}_1= \sigma ^A_u + P_LP_N\sigma ^A_p + \ot {Q_N} 6\sigma ^A_q ,
```
> *"where index `A` corresponds to the considered nuclei and `u`, `p`, `q` define
> the unpolarized, polarized and quadrupolarized contributions. For spin 1/2
> nuclei `Q_N = 0` and for spin 0 nuclei `P_N = Q_N = 0`."*

and in the code the same `qn/6` appears in **every** place a tensor structure
function is contracted, never `qn`:

```fortran
adgh:8648   apptai = un*siau + pl*pn*siap + qn/6.*siaq        ! t-peak tail   (Eq. 37)
adgh:4271   sfm(1)=un*f1+qn/6.*b1                             ! Im_1          (Eq. A.3)
adgh:4272   sfm(2)=epsi*(un*f2+qn/6.*b2)                      ! Im_2
adgh:1335   if(isf.ge.5)ppol=qn/6                             ! exact tail    (Eq. 18)
adgh:0830   if(isf.ge.5)ppol=qn/6                             ! Born          (Eqs. 9,10)
```

> **Verdict: CONFIRMED.** No replacement. `Q_N/6` is POLRAD's, not a choice, and
> the design's `Q_N ≡ P_zz^eff` identification is used consistently in the Born
> (Eqs. (9)/(10)) and in the tail (Eq. (37)) by the same code path.

### 5b. Eq. (A.4) at `Q_N = 0` — confirmed, and it does fix the normalisation

`polrad2t.tex` lines 2698–2709 (Eq. (A.4), spin-1):

```tex
\Im ^{el}_{1}= {1\over 6} \eta_A F^{2}_{m}(4(1+\eta_A)+\eta_A Q_N),
\Im ^{el}_{2}=\pmatrix{F^{2}_{c}+{2\over 3}\eta_A F^{2}_{m}+{8\over 9}\eta_A^2F^{2}_{q}}
   +{Q_N\over 6}\pmatrix{\eta_A F_m^2+ {4\eta_A^{2}\over 1+\eta_A}({\eta_A\over 3}F_q+F_c-F_m)F_q} ,
```

At `Q_N = 0`:

```
  Im^el_1|_0 = (2/3) eta_A (1 + eta_A) F_m^2
  Im^el_2|_0 = F_c^2 + (2/3) eta_A F_m^2 + (8/9) eta_A^2 F_q^2
```

which are exactly `W₁^el = B/2` and `W₂^el = A` of the spin-1 Rosenbluth pair
`dσ/dΩ = σ_Mott[A + B tan²(θ/2)]`, `A = G_C² + (8/9)η²G_Q² + (2/3)ηG_M²`,
`B = (4/3)η(1+η)G_M²`, `η = Q²/4M_A²`. So `F_c ≡ G_C`, `F_m ≡ G_M`, `F_q ≡ G_Q`
in the textbook convention — **the design's claim, verbatim, is right.**

The FORTRAN says the same in a different variable set. `strf` for `ita=2`,
`targ_d` (adgh:4231-4232) sets `fc = fcdeu*amt` etc. and

```fortran
           f1=4./3.*tau*tau1*fm**2
           f2=4./9.*tau*(8.*tau**2*fq**2+6.*tau*fm**2+9.*fc**2)
```

with `sfm(1)=f1`, `sfm(2)=epsi*f2`, `epsi = 2M_A²/t = 1/(2τ)`. Then
`sfm(1) = 2M_A²·[(2/3)τ(1+τ)F_m²]` and
`sfm(2) = 2M_A²·[F_c² + (2/3)τF_m² + (8/9)τ²F_q²]` — i.e. `ℑ_code = 2M_A²·ℑ_paper`
for both, so `ℑ₁|₀ = B/2` and `ℑ₂|₀ = A` are reproduced exactly.

**And POLRAD's own `ffdeu` pins the three normalisations numerically.** It
carries `data dmu/0.857406d0/ dqu/25.84d0/` (adgh:5489) and imposes the `t → 0`
constraints `bbb=(2−μ_d M_d/m_p)/(2√2 M_d)` (adgh:5514),
`ccc=(1−μ_d M_d/m_p−Q_d)/(4M_d²)` (adgh:5519). Compiled standalone with
`gfortran` and evaluated at `t → 0`:

| `t` (GeV²) | `G_C` | `G_M` | `G_Q` |
|---|---|---|---|
| 1e-4 | 0.9982025906 | 1.7108112026 | 25.7829474869 |
| 1e-8 | 0.9999998200 | 1.7139607479 | 25.8399942818 |
| 1e-12 | **1.0000000000** | **1.7139610634** | **25.8399999994** |

against `μ_d·M_d/m_p = 1.7139610634` and `Q_d M_d² = 25.84`. So

```
  F_c(0) = G_C(0) = 1      (= Z; the deuteron has Z = 1)
  F_m(0) = G_M(0) = (M_A / m_p) * mu_A / mu_N          [ = 1.7139610634 for d ]
  F_q(0) = G_Q(0) = M_A^2 * Q_A                        [ = 25.84         for d ]
```

The design's §2.1 validation numbers (`G_Q(0) = 25.829`, `G_M(0) = 1.7139`) match
POLRAD's own to 6 digits on `G_M` and 0.04 % on `G_Q` (POLRAD hard-codes
`Q_d M_d² = 25.84`; the design derives 25.829 from `Q_d = 0.2859 fm²`).

> **Verdict: CONFIRMED.** No replacement. **T9** and **T10** as written in
> design §5 are the right gates and will pass.
>
> **One `Z` caveat, worth a comment in `HoSpin1FF`.** POLRAD factors the charge
> out **explicitly** for spin-½ (Eq. (A.5), `ℑ₁ = Z²η G_m²`) and for spin-0
> (`strf` adgh:4251, `f2=4.*amp2*tau*(tarz*ff)**2` with `ffco` returning
> `ff(0) = 1`), but **not** for spin-1 — Eq. (A.4) has no `Z²` because POLRAD's
> only spin-1 target is the deuteron, `Z = 1`. The design's choice —
> `F_c(0) = Z = 3` folded into the form factor, no separate `Z²` — is the correct
> way to extend it, since `A(0) = G_C(0)² = Z²` either way. **Do not do both.**

---

## 6. Compact verdict table

| # | item | verdict | replacement |
|---|---|---|---|
| 1 | `σ_u^C` carries `Y₋` | **corrected** | `σ_u^C = (α³/S) Z² **Y₊** ∫_{η_min}^{η_max} (dη_A/η_A) X̃ F²`, `Y₊ = (1+(1−y)²)/(1−y)` |
| 2 | `X₁`'s linear `−4η_A` | **confirmed** | — (`X₁ = x_A² + 4x_Aη_A − 4η_A`, `X̃ = X₁/(2η_A x_A²)`) |
| 3 | `(3/4)x²` / `(4/3)η_A F_q` in `σ_q^d` | **confirmed** | — (design's block is correct; the *PDF* has them swapped) |
| 4 | Eq. (18) `1/A²` right, `1/A` left; `x_A = x/A` | **confirmed** | — (refinements, as REVISED 2026-09-03: LiPolGen keeps `x_A = x/A`, `S_A = A s` — both exact for its collider variables; the per-nucleon factor is **`1/A²`**, i.e. `(1/A)×Jacobian`, and the old "`(m_p/M_A)(1/A)` differs by 0.5 %" line was wrong by `A m_p/M_A = 5.97`; exact `η` limits — §4, §8b) |
| 5 | `Q_N/6` in Eq. (37); `Q_N=0` Rosenbluth limit of (A.4) | **confirmed** | — (`ℑ₁\|₀ = B/2`, `ℑ₂\|₀ = A`; `F_c(0)=Z`, `F_m(0)=(M_A/m_p)μ_A/μ_N`, `F_q(0)=M_A²Q_A`) |

---

## 7. Two further errors found while checking

### 7.1 DESIGN ERROR — `ℑ^el_6` in §1.4.6 has a spurious `η_A`

The design (§1.4.6, Eq. (A.4) transcription) prints

> `Im^el_6 = (Q_N/24) ( F_m^2 + (4 eta_A/(1 + eta_A)) ( (eta_A/3) F_q + F_c + eta_A F_m ) F_q )`

The source, `polrad2t.tex` lines 2706–2707, has **no `η_A` in that numerator**:

```tex
\Im ^{el}_{6}= {Q_N\over 24}\pmatrix{F_m^{2}
    +{4\over 1+\eta_A}({\eta_A\over 3}F_q+F_c+\eta_A F_m)F_q} ,
```

Confirmed independently from the FORTRAN: `ℑ₆ = (Q_N/6)ε³(b₂/3 + b₃ + b₄)`
(Eq. (A.3)) with `strf`'s elastic `b₂, b₃, b₄` (adgh:4235-4242) and `ε = 1/(2τ)`
gives, after the `fm²` terms collapse to `4τ³F_m²` and the `F_q` terms to
`(16/3)(τ³/τ₁)(τF_q + 3F_c + 3τF_m)F_q`,

```
  eps^3 (b2/3 + b3 + b4) = (1/2)[ F_m^2 + (4/(1+eta)) ((eta/3) F_q + F_c + eta F_m) F_q ] * 2 M_A^2
```

i.e. exactly `2M_A²·(Q_N/24)(…)` with **`4/(1+η_A)`**. (`ℑ₂`'s `4η_A²/(1+η_A)`,
`ℑ₅`, `ℑ₇`, `ℑ₈` all check out against `b₁…b₄` the same way — only `ℑ₆` is wrong
in the design.) **Replacement:**

```
  Im^el_6 = (Q_N/24) ( F_m^2 + (4/(1 + eta_A)) ( (eta_A/3) F_q + F_c + eta_A F_m ) F_q )
```

This is dead code in v0 (Eq. (A.4) is only evaluated on the `PolradFull` path,
and T9's `PolradFull` half is `SKIP`-ped) — but it is exactly the kind of error
that would be invisible until the upgrade path is built, so fix the doc now.

> **APPLIED 2026-09-04 (task B5).** `design_C_tensor_rc.md` §1.4.6 now prints
> `4/(1 + eta_A)`, with a dated correction box beside it recording what it said
> before. Re-verified before changing anything: the `\Im ^{el}_{6}` line of
> `polrad2t.tex` (lines 2705–2706) has no `η_A` in that numerator, the FORTRAN
> collapse of `ε³(b₂/3 + b₃ + b₄)` reproduces `4/(1+η)` exactly, and the CODE
> was confirmed **not** to contain the error — it does not contain `ℑ^el_6` at
> all. `RcTailModel::PolradFull` throws `NOT IMPLEMENTED` (`src/core/rc.cpp`
> constructor), and what v0 evaluates is Eq. (38) (`polrad_sigma_el_u` and the
> `σ_q` integrand) plus the `Q_N = 0` Rosenbluth pair (`rosenbluth_spin1`).
> Document-only defect; no test number moved.

### 7.2 POLRAD BUG — the `approx` carbon tail is missing `Z²`

`elu`'s carbon branch (adgh:8903-8905) is `elu = xxt*ff**2/eta` with `ffco`
normalised to `ff(0) = 1` (adgh:5576-5590, no `tarz` anywhere in it), and
`apptai` adds no charge factor. The **paper** has `Z²` explicitly
(`\od{\alpha^3}{S}Z^2Y_-`, tex:928) and the **exact** path has it
(`strf`, adgh:4251: `f2=4.*amp2*tau*(tarz*ff)**2`). So the closed-form carbon
entry in POLRAD 2.0's `approx` decks is low by `Z² = 36` for carbon.
`elual2` (adgh:9269-9270) repeats the omission.

The physics is not in doubt — Eq. (A.5)'s explicit `Z²`, `strf`'s `(tarz*ff)²`
and Eq. (38)'s printed `Z²` agree — so **design C's T8(b) carbon `Z²` gate stands
as written**; it just cannot be validated by diffing against `elu`. Note the two
POLRAD errors point the same way: the `approx` carbon entry was evidently never
exercised (see the `c+use,approx` in §0), which is also the most likely origin of
the `Y₋` typo in §1.

---

## 8. NEW BLOCKING FLAG — the overall sign of Eq. (38)

Not one of the five questions, found while checking them, and it gates
`rc_tail`'s sign.

**`X̃` is negative over essentially the whole integration range.**
`X₁(η_A) = x_A² − 4η_A(1 − x_A)` is *decreasing* in `η_A` and vanishes at the
ultrarelativistic `η_min = x_A²/(4(1−x_A))`, so `X̃ = X₁/(2η_Ax_A²) < 0` for every
`η_A > η_min`. Consequently **Eq. (38) as printed returns a negative number** for
every unpolarised entry.

Measured, using POLRAD's own `elu`/`elq`/`ffdeu` compiled standalone with POLRAD's
own limits and prefactors (deuteron, `M_d = 1.8756280`):

| beam / `x` / `y` | `∫dη elu` | `σ_u^d` | `σ_q^d` | `σ_q/σ_u` |
|---|---|---|---|---|
| 27.6 GeV, 0.05, 0.60 | −4.320e+03 | **−2.352e−05** | −2.492e−06 | +0.1060 |
| 27.6 GeV, 0.012, 0.50 | −2.133e+05 | **−1.001e−03** | −6.429e−05 | +0.0642 |
| 11.0 GeV, 0.20, 0.50 | −1.333e+01 | **−1.570e−07** | +1.829e−08 | −0.1165 |

(`α³/S · Y₊ · ter`, `ter = m_p/M_d`, no `barn`; the sign is the point, not the units.)

This is *consistent* with Eq. (18), which carries an explicit leading minus
(`σ₁^el = −(α³y/A²)∫…`, tex:580) that the ultrarelativistic Eq. (38) does not
print — i.e. POLRAD's `σ₁^el`, and therefore its `σ_u^A`, `σ_q^A`, is a
**negative-signed object in this convention**. The main program nevertheless adds
it (`sig = sib*extai1*(…) + sia + tai(1) + (tai(2)*extai2 + tai(3)*extai3)/tara`,
adgh:488-490) with no sign flip, and `bornin` returns a positive `sibor`
(adgh:833, `sibor=ssum*an/y/y*2.` with `an > 0` and `tm(1) = y − 2m² > 0`, adgh:817).
**Whether POLRAD's shipped total therefore subtracts the elastic tail was not
resolved here** — resolving it needs a PATCHY build and a run, which is outside
this note's read-only scope, and the `approx` path it would exercise is not the
default build anyway (§0).

**What design C must do about it.** §1.4.5 writes

```
  w_tail = 1 + [ sigma^el_U + (Q_N/6) sigma^el_T + kappa_qe * sigma^q_U ] / [ x s dsigma_unpol (1 + w_avg) ]
```

which needs `σ^el_U > 0` — a radiative tail *adds* events to the DIS bin, it
cannot remove them. A literal transcription of Eq. (38) gives `σ^el_U < 0` and
`w_tail < 1`, i.e. a tail that *reduces* the yield and therefore *anti*-dilutes
`A_zz`. So:

* **take `σ_u^A = −(α³/S) Y₊ ∫ dη_A/η_A [ … ]`** (equivalently, integrate `|X̃|`),
  matching Eq. (18)'s leading minus;
* **add T8(0), a sign gate, before T8(a):** assert `σ^el_U > 0` and
  `w_tail ≥ 1` at three `(x, Q²)` points spanning the grid. It is one line and it
  is the single cheapest guard against shipping a tail that anti-dilutes;
* keep the **ratio** `σ_q^d/σ_u^d` as the tensor-tail gate rather than the
  absolute normalisation — it is sign-convention-free, and the three values above
  (+0.106, +0.064, −0.117) are usable reference numbers for a **deuteron**
  regression before ⁶Li form factors exist. Note it **changes sign** between
  `x = 0.05` and `x = 0.20`: the tensor tail is not a fixed fraction of the
  unpolarised one, which is a second reason §1.4.5 must carry `σ^el_T`
  separately rather than as a scale factor on `σ^el_U`.

### 8b. The derivation that fixes BOTH the sign and the whole-nucleus normalisation (added 2026-09-03)

The sign flag above was closed by inspection ("a tail cannot remove events").
There is a stronger statement available, and it also settles §4's per-nucleon
question, so it is recorded here rather than left to the next reader.

**The t-peak IS Weizsäcker–Williams × Compton.** At small photon virtuality the
nucleus radiates an equivalent photon flux

```
  dn = (Z^2 alpha / pi) (dz/z) (dt/t) (1 - t_min/t) F^2(t)          [ z << 1 ]
```

and the lepton Compton-scatters off it with

```
  d sigma / d t^  =  (2 pi alpha^2 / s^^2) [ (1 + (1-y)^2)/(1-y) ]  =  (2 pi alpha^2/s^^2) Y_+ .
```

Folding the two with `dz dQ^2 = x_A S_A dx_A dy` gives the **whole-nucleus**

```
  d^2 sigma / (dx_A dy)  =  (2 alpha^3 Z^2 / (x_A^2 S_A)) Y_+
                            INT (d eta / eta) (1 - eta_min/eta) F^2(4 M_A^2 eta) .
```

Evaluated numerically against `polrad_sigma_el_u` for a **spin-0 Gaussian** at
eight `(Z, x, y)` points, the ratio to `−`Eq. (38) is **exactly `(1 − x_A)`**
(0.99917 at `x_A = 8.3e−4`, 0.99167 at `8.3e−3`) — the `O(x_A)` difference being
the exact-vs-ultrarelativistic `η` limits and the `X̃` bracket's own `x_A`
dependence. Two conclusions, both used by `src/core/rc.cpp`:

1. **The leading minus is real.** Eq. (38) as printed integrates a manifestly
   negative `X̃`; the WW × Compton form is manifestly positive. `rc.cpp` applies
   Eq. (18)'s minus and `tests/test_rc.cpp` T8(0) gates it.
2. **Eq. (38) is the WHOLE-NUCLEUS `d²σ/(dx_A dy)` — there is no `1/A` in it.**
   The paper's `σ₁^el = (1/A) d²σ/dx_A dy` (polrad2t.tex:579) is loose notation,
   and POLRAD's own FORTRAN confirms it by applying **two** reductions to the
   same `∫elu`: `ter = m_p/M_A` in `apptai` (adgh:8607-8620) *and* `/tara` in
   the main program (adgh:489, 497), against a per-nucleon Born. §4's third
   refinement, which said the two conventions "differ by 0.5 %", is corrected
   in place: they differ by **`A m_p/M_A = 5.97`**, and design §1.4.5's single
   `1/A` is a **Jacobian**, not a per-nucleon division. The per-nucleon tail is
   `(1/A)·(dx_A/dx)·[Eq. (38)]`, i.e. **`[Eq. (38)]/A²`** with `x_A = x/A`.

**Design §1.4.5 and §1.4.6 are marked CORRECTED** for this: §1.4.5's "the
per-nucleon tail is `(1/A)·[Eq. (38)]`" needs the second `1/A`, and §1.4.6's
reading of Eq. (18)'s left-hand `1/A` as "already the per-nucleon-normalised
object" is what made the single factor look sufficient.

### 8c. A second flag, of the same kind: the t-peak is not always the tail (added 2026-09-03)

§0 already warns that the `approx` decks are "not a routinely exercised" path.
`tests/test_rc.cpp` T8(c) now measures POLRAD §2.1.3 B's *"contributions of s−
and p−peaks are suppressed"* with an independent leading-log construction:

| target, point | `t / (t + s + p)`, elastic | quasi-elastic |
|---|---|---|
| ⁶Li, EIC config 1, `Q² ≥ 20 GeV²`, `y ≤ 0.9` — **per cell**, 1356 cells (CORRECTED 2026-09-04) | **1.0000** (minimum) | **1.6458e−04** (minimum); **275 of the 1356 cells are below 0.96** (282 at the shipped `k_F`) |
| ⁶Li, EIC config 1, `x = 0.01, y = 0.985` (`Q² = 39`) | > 0.999 | **0.44** |
| ⁶Li, EIC config 1, `x = 0.01, y = 0.1` (`Q² = 4`) | > 0.999 | **0.23** |
| deuteron, `E = 27.6`, `x = 0.012, y = 0.85` (`Q² = 0.53`) | 0.65 | **0.17** |
| proton, `E = 10`, `x = 0.3, Q² = 2` | — | **0.037** |

**Row 1 was rewritten on 2026-09-04 and the old version is withdrawn.** It
read `> 0.999` / `> 0.96` as if both held over the whole window; the second
half is false. The **elastic** column survives verbatim and is in fact
stronger than it looked — ⁶Li's coherent form factor is dead at `Q′² ≈ Q² ≥ 20`,
so `ll_peaks_spin1` returns essentially nothing there and the minimum of
`t/(t+s+p)` over all 1356 accepted cells is **1.0000**. The **quasi-elastic**
column is not a window statement at all: its minimum is **1.6458e−04** (at
`x = 0.7244`, `y = 0.0070`) and **275** of the 1356 cells sit below 0.96 at
this section's own `kf_gev = 0` convention, **282** at the shipped
`k_F = 0.169 GeV`. On the *total* tail at this section's own `kf_gev = 0`,
`t/(t+s+p) < 0.99` in **346** of the 1356; 21 of those 346 sit at `y ≥ 0.1`,
five of them at `y = 0.845–0.892`, i.e. at the `y → 1` end as well, so they
are not confined to the low-y corner.

**Two censuses, on two conventions, and they are not nested.** The 346 above is
`t/(t+s+p) < 0.99` at `kf_gev = 0` (unsuppressed; this section's convention, and
a measurement record). The figure every doc site carries is the *other* one:
`phase_B_numbers.md` §B2.3(i)'s **331 of 1356 cells with
`RcModel::tail_ratio > 1.01` at the shipped `k_F = 0.169 GeV`, 318 of them at
`y < 0.1`**. 346 > 331 under a *stricter* cut is only possible because the two
conventions differ — the 318/331 are **not** a subset of the 346, and an
earlier form of this paragraph read them as one (corrected 2026-09-05). (An earlier revision of this paragraph read "the failures
are at `y → 0`, not at `y → 1`"; that contradicted this document's own
`y → 1` rows below and is withdrawn, 2026-09-04.) Measured on the 101 × 77 CLI
grid and set out in `../run_2026-09-03/phase_B_numbers.md` §B2.3(i). It is a
measurement record, **not** a gated claim: `tests/test_rc.cpp` T8(d)(i) gates
the combined `tail_ratio` census on the test file's own 40 × 24 sampler
(45 of 188 there) and never forms `t/(t+s+p)`, so no assertion in the suite
pins the numbers in this row.

So POLRAD's sentence is right **where the nuclear form factor is dead at the
s-peak's `t ≈ z Q²` and the quasi-elastic nucleon form factor is dead with
it**, and badly wrong where either one is alive. For ⁶Li at EIC `Q²` the
t-peak is the whole tail **event-weighted** and **not cell by cell**; at
fixed-target `Q²`, for the *nucleon* form factors of the quasi-elastic tail,
at the `y → 1` edge of this generator's own window and — the wider failure —
all along its `y → 0` edge, it is not.

Every row is at T8(c)'s own convention, the **unsuppressed** quasi-elastic tail
(`kf_gev = 0`), on both the t-peak and the peaks it is compared with. At the
shipped Pauli `k_F = 0.169 GeV` the t-peak's own quasi-elastic piece is roughly
halved and every quasi-elastic entry falls with it — 0.44 → **0.27** in row 2,
0.23 → **0.12** in row 3 — so these are the *conservative* readings.

**Row 2 is the one added on 2026-09-04 and it narrows the claim.** At
`y = 0.985` the t-peak carries only **68 %** of the *total* tail (63 % at the
shipped `k_F`) even though `Q² = 39 GeV² ≥ 20`, because `z_s = (1−y)/(1−x_A y) → 0` drives the elastic
vertex to `Q′² → 0`, where the form factor is 1 and the Born-like `1/Q′⁴` is
enormous — the s-peak grows **faster** than the t-peak's own `Y₊ ~ 1/(1−y)`.
So "the t-peak is > 99 % at `Q² ≥ 20 GeV²`" is an **event-weighted** statement
and never a per-cell one. It fails at the `y → 1` rows for the reason just
given, and it fails far more widely at `y → 0`: of the 331 cells in the
`Q² ≥ 20`, `y ≤ 0.9` window whose two tail models differ by more than 1 %,
**318 sit at `y < 0.1`** — a different mechanism (the uncancelled soft
radiator, `D(z_s) ∝ 1/y`, against a t-peak crushed by `t_min ∝ x²` because
`Q² ≥ 20` forces low `y` to mean high `x`). The `y → 1` edge is the smaller
exception, not the only one; this paragraph said otherwise until 2026-09-04.

**What changed in the code (B2, 2026-09-04).** The leading-log construction is
no longer test-local: `ll_radiator`, `rosenbluth_spin1`, `rosenbluth_nucleon`,
`dsigma_el_dq2`, `ll_peaks_spin1` and `ll_peaks_qe` moved into
`src/core/rc.cpp`, so T8(c) now measures the **shipped** code path, and
`RcTailModel::TPeakPlusLL` tabulates them beside the t-peak. The default stays
`TPeak`, so `rc_tail` is still a **lower bound** unless the run asks otherwise.
Two things this section must not be read as saying:

* the sum is a **stated model**, not a controlled O(α) expansion — an η_A
  quadrature plus a single-z leading log is accurate to the worse of the two
  (~5–10 %), and the radiator's soft `1/(1−z)` is **uncancelled**, so it
  overshoots as `y → 0`;
* there is **no tensor s/p peak** in POLRAD and none is invented here, so the
  s+p enter the unpolarised numerator only and `TPeakPlusLL` **lowers the
  tensor fraction** of the tail. That fraction is unknown, not zero.

`rc.hpp`'s `RcTailModel` comments, the run banner and `USAGE.md` §7b all say so.

---

## 9. What to change in `design_C_tensor_rc.md`

1. **§1.4.2 Eq. (38) block, `σ_u^C` line:** `Z^2 Y_-` → `Z^2 Y_+`.
2. **§1.4.2 flag 1 (line 387–389):** replace "transcribe as printed, then check
   it against the FORTRAN" with the resolved statement + the `apptai` citation.
3. **§1.4.2 flags 2 and 3:** mark **confirmed**, keep the `pdftotext` warning of
   flag 3 but note the PDF renders the two fractions *swapped*.
4. **§1.4.2 integration limits:** `INT eta_m … INF` → the exact `(eta1, eta2)`
   of adgh:8611-8613, with `η_min(u.r.)` retained as a cross-check (§4).
5. ~~**§1.4.1 nuclear-map row:** `x_A = x/A` → `x_A = x·m_p/M_A`~~ **WITHDRAWN
   2026-09-03.** `x_A = x/A` and `S_A = A·s` are EXACT for this library's
   collider variables (`beams.cpp`'s per-nucleon `p_u = p_A/A`), and POLRAD's
   `x·m_p/M_A` is exact only for its own fixed-target ones. Keep the design's
   map. **§1.4.5's `1/A` is a Jacobian and needs a SECOND `1/A`** — the
   per-nucleon tail is `[Eq. (38)]/A²`; see §4's revised third bullet and §8b.
6. **§1.4.6 `ℑ^el_6`:** drop the spurious `η_A` (§7.1).
7. **§5:** add **T8(0)**, the tail-sign gate (§8); note that **T8(b)** cannot be
   diffed against `elu` because POLRAD's `approx` carbon branch omits `Z²` (§7.2).
8. **§7 Q2:** mark **CLOSED**. **§6 step 0:** mark **DONE**, and correct the
   catalogue identifier `ADXQ` → **`ADGH_v1_0`**, dataset `37vgvzgr2w` on
   Mendeley Data (only lines 1964 and 2096 of this design carry `ADXQ`; no other
   file under `docs/` mentions it).
