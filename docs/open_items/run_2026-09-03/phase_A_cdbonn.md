# A5 — the CD-Bonn deuteron wave function in momentum space

**Task.** Gate condition 3 of open item 10 (`docs/OPEN_ITEMS_SOLUTIONS.md` §10,
"The seven close conditions, and where each one landed", row 3 — a line
citation into that file was stale by three sections and is replaced by the
section name, 2026-09-05)
wants "a real CD-Bonn u, w instead of the D-state rescaling proxy of item 5"
(`docs/open_items/run_2026-09-02/phase_D_gate.md:322-334`). This document
supplies it: the coefficients, read from Machleidt's own source; every
verification number computed rather than quoted; and the difference against
both AV18 and the proxy **measured through the repository's own A = 2 gate**,
not asserted.

**Status.** The coefficients are obtained and verified. **No repository file was
edited** — this file is the only addition. Nothing is wired in; A5's code step
is still to be done, and §9 says exactly what it should contain.

**Headline.** The proxy has the **wrong sign of the effect**. Against
`fdeut.av18` + CT18NLO at CDKS's own κ (Eq. 21), the real CD-Bonn wave function
moves the G3b peak ratio from **0.7193 to 0.8642** (a factor 1.39 low becomes a
factor 1.16 low) and **restores the low-x zero that CT18NLO had pushed below
the scan floor** (back to x = 0.043, giving two again) — so G3a's counting clause, which `phase_D_gate.md` recorded as
regressing, passes again. The proxy moved the same number the other way, to
0.6326, and left the counting clause failing. The sentence at
`phase_D_gate.md:330-331` — *"The residual factor is therefore **not** explained
by the wave function; if anything CDKS's softer D state should make their curve
smaller, not larger"* — is **contradicted by the real wave function**, and the
reason is §7: CD-Bonn's D wave is not uniformly smaller, it is *larger* below
k ≈ 1 fm⁻¹ and only smaller above it.

---

## 1. Provenance

| what | value |
|---|---|
| paper | R. Machleidt, *The high-precision, charge-dependent Bonn nucleon-nucleon potential (CD-Bonn)*, **Phys. Rev. C 63, 024001 (2001)** |
| e-print | `arXiv:nucl-th/0006014v1`, 9 Jun 2000 |
| how fetched | `curl -sL https://arxiv.org/e-print/nucl-th/0006014` → gzipped tar of the LaTeX source (HTTP 200, 488053 bytes, md5 `e0c0896bfa68a3173b01a8ec68248c93`); single source file `bonn.tex`, 4237 lines, md5 `d701fabc639fae71e26a72c0f9484c46` |
| independent cross-check | `curl -sL https://arxiv.org/pdf/nucl-th/0006014v1` (69 pages, md5 `ac0ce5f5c4a92c26291f8b544669d5c1`) → `pdftotext -layout` |
| coefficients | **Table XX**, "Coefficients for the parametrized deuteron wave functions (n = 11)". LaTeX label `tab_dwpar`, `bonn.tex:4005-4029`; PDF text line 2749 |
| parameterisation | **Appendix D**, Eqs. **(D19)–(D25)**; `bonn.tex:2817-2857`; PDF text lines 1905-1960 |
| γ | Eq. **(D6)**, `bonn.tex:2708`; PDF line 1809 |
| deuteron properties (for checking) | **Table XV**, label `tab_deu`, `bonn.tex:3820-3843`; PDF line 2610 |
| numerical r-space u(r), w(r) (for checking) | **Table XIX**, label `tab_dwaves`, `bonn.tex:3953-4001`; PDF line 2707 |

**Correction to the plan.** `docs/open_items/run_2026-09-03/PLAN.md:70-71` and
the A5 task statement say *"Tables XVII/XVIII"*. That is **wrong**. In both the
e-print and its PDF, Table XVII is "Parameters of the scalar isoscalar bosons
… for the T = 0 np potential" and Table XVIII is "Coupling constants of the
scalar isoscalar bosons …". The deuteron material is **Table XIX** (numerical
u, w) and **Table XX** (the parameterisation coefficients). I counted the
twenty `\begin{table}` environments in `bonn.tex` in order and confirmed the
numbering against the PDF's own printed `TABLE XX` heading. *I have not seen
the published PRC 63 024001 page proofs, so I cannot rule out a renumbering
between e-print and journal; the LaTeX labels `tab_dwpar` / `tab_dwaves` are
the unambiguous handles and are what a header comment should cite.*

**Extraction method.** The LaTeX source is machine-readable, so no digit was
OCR'd or transcribed by eye. The table body was read out of `bonn.tex` and
then **the same twenty numbers were read independently out of `pdftotext`
output** and compared — they agree character for character. Both listings are
reproduced verbatim in §3.

---

## 2. The parameterisation, verbatim (Appendix D)

r-space ansatz, Eqs. (D19)/(D20):

```
u_a(r) = sum_{j=1}^{n} C_j exp(-m_j r)
w_a(r) = sum_{j=1}^{n} D_j exp(-m_j r) [ 1 + 3/(m_j r) + 3/(m_j r)^2 ]
```

momentum space, Eqs. (D21)/(D22) — **this is the object the gate needs**:

```
psi_0^a(q) = (2/pi)^{1/2} sum_{j=1}^{n} C_j / (q^2 + m_j^2)
psi_2^a(q) = (2/pi)^{1/2} sum_{j=1}^{n} D_j / (q^2 + m_j^2)
```

masses, Eq. (D25), with m_0 = 0.9 fm⁻¹ and γ from Eq. (D6):

```
m_j = gamma + (j-1) m_0 ,   gamma = 0.2315380 fm^-1 ,   n = 11
```

The boundary conditions u_a(r) → r and w_a(r) → r³ as r → 0 give one
constraint on the C_j (Eq. D23) and three on the D_j (Eq. D24 plus two
circular permutations of n−2, n−1, n):

```
C_n = - sum_{j=1}^{n-1} C_j                                          (D23)

                m_{n-2}^2
D_{n-2} = --------------------------- x
          (m_n^2 - m_{n-2}^2)(m_{n-1}^2 - m_{n-2}^2)

          [ -m_{n-1}^2 m_n^2 sum_{j=1}^{n-3} D_j/m_j^2
            + (m_{n-1}^2 + m_n^2) sum_{j=1}^{n-3} D_j
            - sum_{j=1}^{n-3} D_j m_j^2 ]                            (D24)
```

**Equivalent and simpler form of (D24), derived and checked here.** Expanding
`f(x) = e^{-x}(1 + 3/x + 3/x²) = 3/x² − 1/2 + x²/8 − x³/15 + …` (the 1/x term
cancels identically), the three conditions are just

```
sum_j D_j / m_j^2 = 0      (kills the 3/r^2 singularity)
sum_j D_j         = 0      (kills the constant)
sum_j D_j m_j^2   = 0      (kills the r^2 term; the r^1 term vanishes by itself)
```

i.e. a 3×3 linear system for D_9, D_10, D_11. **Both routes were implemented
and agree to 1.0e−14 relative** (§4), so a C++ header may use either; the
linear system is the one that cannot be mis-permuted. (Analogously,
Eq. (D23) is `sum_j C_j = 0`, from u_a(0) = 0.)

Machleidt's warning, quoted because it is a real trap: *"The constraints
Eqs. (D23) and (D24) must be enforced by double precision (i.e., to about 15
decimal digits), otherwise the wave function is not reproduced correctly for
r ≤ 0.5 fm. This applies, particularly, to the D wave."*

---

## 3. The coefficient table, verbatim

**As it appears in `bonn.tex` (`tab_dwpar`)**, Fortran `D` exponents:

```
      1  $   0.88472985D+00 $  $   0.22623762D-01 $
      2  $  -0.26408759D+00 $  $  -0.50471056D+00 $
      3  $  -0.44114404D-01 $  $   0.56278897D+00 $
      4  $  -0.14397512D+02 $  $  -0.16079764D+02 $
      5  $   0.85591256D+02 $  $   0.11126803D+03 $
      6  $  -0.31876761D+03 $  $  -0.44667490D+03 $
      7  $   0.70336701D+03 $  $   0.10985907D+04 $
      8  $  -0.90049586D+03 $  $  -0.16114995D+04 $
      9  $   0.66145441D+03 $     Eq. (D24)
     10  $  -0.25958894D+03 $     Eq. (D24)
     11     Eq. (D23)              Eq. (D24)
```

**The same twenty numbers as `pdftotext` reads them from the PDF** (independent
path, identical digits):

```
 1   0.88472985D + 00     0.22623762D − 01
 2  −0.26408759D + 00    −0.50471056D + 00
 3  −0.44114404D − 01     0.56278897D + 00
 4  −0.14397512D + 02    −0.16079764D + 02
 5   0.85591256D + 02     0.11126803D + 03
 6  −0.31876761D + 03    −0.44667490D + 03
 7   0.70336701D + 03     0.10985907D + 04
 8  −0.90049586D + 03    −0.16114995D + 04
 9   0.66145441D + 03        Eq. (D24)
10  −0.25958894D + 03        Eq. (D24)
11      Eq. (D23)             Eq. (D24)
```

### 3.1 In a form a C++ header can be typed from

Units: C_j, D_j in fm^(−1/2); m_j and q in fm⁻¹; ψ in fm^(3/2).

```cpp
// arXiv:nucl-th/0006014 (Phys. Rev. C 63, 024001 (2001)), Table XX
// (LaTeX label tab_dwpar).  ONLY j = 1..10 of C and j = 1..8 of D are
// published; the rest come from Eqs. (D23)/(D24), computed at run time.
constexpr double kCdBonnGammaFm  = 0.2315380;  // Eq. (D6)
constexpr double kCdBonnM0Fm     = 0.9;        // Eq. (D25)
// m_j = kCdBonnGammaFm + (j-1) * kCdBonnM0Fm,  j = 1..11
constexpr double kCdBonnC[10] = {   //  fm^-1/2, j = 1..10
   0.88472985e+00, -0.26408759e+00, -0.44114404e-01, -0.14397512e+02,
   0.85591256e+02, -0.31876761e+03,  0.70336701e+03, -0.90049586e+03,
   0.66145441e+03, -0.25958894e+03 };
constexpr double kCdBonnD[8]  = {   //  fm^-1/2, j = 1..8
   0.22623762e-01, -0.50471056e+00,  0.56278897e+00, -0.16079764e+02,
   0.11126803e+03, -0.44667490e+03,  0.10985907e+04, -0.16114995e+04 };
```

**Derived, NOT published** — compute them; these values are for a pin test
only. Doubles, 17 significant digits, from the constraints of §2:

```
m_j (fm^-1):  0.2315380  1.1315380  2.0315380  2.9315380  3.8315380  4.7315380
              5.6315380  6.5315380  7.4315380  8.3315380  9.2315380

C_11 =    42.260718143999895      (exact decimal: 42.260718144)
D_9  =  1374.0642726388592
D_10 =  -630.43718361517256
D_11 =   120.6876428043135
```

The 60-decimal-digit values (Python `decimal`, prec 60, from the published
8-digit inputs) are

```
D_9  =  1374.06427263886393575448289251797204699873969261256402572772
D_10 =  -630.437183615181024105121698354854794514989447431622218746441
D_11 =   120.687642804317088350638805836882747516249754819058193018765
```

so the double-precision solve is right to **3.4e−15 / 1.3e−14 / 3.0e−14**
relative — comfortably inside Machleidt's "about 15 decimal digits" demand.
`C_11` summed in double lands 1.05e−13 (2.5e−15 relative) below the exact
decimal 42.260718144; irrelevant, but recorded so nobody "fixes" it.

Residuals of the four constraints at these doubles:

```
sum_j C_j       =  0.000e+00
sum_j D_j/m_j^2 =  2.665e-15
sum_j D_j       =  1.705e-13
sum_j D_j m_j^2 = -1.819e-12
```

---

## 4. Sign and normalisation — **measured, and it is not what Eq. (D13) says**

This is the `(−i)^L` trap that `docs/PHYSICS_CHANNELS.md:246` and
`b1_nuclear.hpp:159` already warn about, and CD-Bonn walks straight into it.

Machleidt's Eq. (D13) is printed with a **bare** j_L kernel,
`u_L(r)/r = sqrt(2/pi) ∫ dk k² j_L(kr) ψ_L(k)`, with u_0 ≡ u and u_2 ≡ w.
Taken literally together with (D20)/(D22) it is **inconsistent for L = 2**.
Measured, not argued (`scratchpad/ftcheck.py`, `scipy.special.spherical_jn`):

| r [fm] | u_a(r), Eq. (D19) | bare-j₀ FT of ψ₀ᵃ | w_a(r), Eq. (D20) | bare-j₂ FT of ψ₂ᵃ |
|---|---|---|---|---|
| 1.00 | 0.431118 | 0.431118 | 0.119132 | **−0.119132** |
| 2.00 | 0.510438 | 0.510438 | 0.141346 | **−0.141346** |
| 5.00 | 0.277063 | 0.277063 | 0.038602 | **−0.038602** |

The S wave round-trips; the D wave comes back **negated**. The r-space forms
are the trustworthy pair: §5 shows `u_a`, `w_a` built from the published C_j,
D_j reproduce Machleidt's own Table XIX (where w(r) > 0 beyond r ≈ 0.15 fm),
so the D_j are certainly right *for Eqs. (D19)/(D20)*. It is the momentum-space
sign that carries the `i^L`.

**And `fdeut.av18` uses the other convention.** Measured on the repository's
own file (`scratchpad/conv.py`) — its r block and its k block are related by
Eq. (D13) with a **bare** kernel for *both* L:

| r [fm] | fdeut w(r) | bare-j₂ FT of fdeut w(k) | ratio |
|---|---|---|---|
| 1.00 | 0.156324 | 0.156328 | 1.000 |
| 2.00 | 0.146061 | 0.146058 | 1.000 |
| 5.00 | 0.037839 | 0.037840 | 1.000 |

(the S wave likewise, ratio 1.000). Therefore, **to drop CD-Bonn in where
`read_fdeut_k` reads `fdeut.av18`:**

```
u(p) = + sqrt(2/pi) * sum_j C_j / (p^2 + m_j^2)      // = + psi_0^a
w(p) = - sqrt(2/pi) * sum_j D_j / (p^2 + m_j^2)      // = - psi_2^a  <-- the sign
```

with `int_0^inf dp p^2 [u^2 + w^2] = 1`, which is exactly the convention
`design_D_b1_li6.md:601` states for the gate (φ₀ = u, φ₂ = −w). Get this wrong
and the S–D interference flips; §8 shows what that does.

**Third, independent confirmation — the repository's own sign gate.**
`alpha_d_quadrupole_fm2` on `ClusterPartialWave::from_uw` columns:

| (u, w) fed in | Q returned [fm²] | reference |
|---|---|---|
| `fdeut.av18` (control) | **+0.269362** | its own header `qm` = 0.269673 (0.12 %) |
| CD-Bonn, **w = −ψ₂ᵃ** | **+0.270178** | CD-Bonn Table XV, Q_d = 0.270 |
| CD-Bonn, w = +ψ₂ᵃ | −0.303719 | wrong sign **and** wrong magnitude |

---

## 5. Verification — every number below was computed here

All momentum integrals are done in **closed form**, using
`(2/pi) ∫_0^inf dp p²/[(p²+a²)(p²+b²)] = 1/(a+b)`, so

```
(2/pi) int dp p^2 u^2 = sum_{i,j} C_i C_j / (m_i + m_j)
(2/pi) int dp p^2 w^2 = sum_{i,j} D_i D_j / (m_i + m_j)
```

— no quadrature error at all. r-space integrals are 500-point Gauss–Legendre
on five panels out to 150 fm, with the wave functions evaluated in 60-digit
`decimal` arithmetic (the D wave near r → 0 is a catastrophic cancellation;
Machleidt says so).

| quantity | **computed here** | CD-Bonn's published value | source of the published value |
|---|---|---|---|
| normalisation (2/π)∫dp p²(u²+w²) | **0.999999826159561** | 1 | Eq. (D14) |
| same, in r space ∫dr(u_a²+w_a²) | **0.999999826** | 1 | Parseval — the two agree to all 9 digits |
| S-state part | 0.951437753 | — | |
| **P_D** = (2/π)∫dp p²w² | **0.048562073645887 → 4.856207 %** | **4.85 %** | Table XV |
| **η** = A_D/A_S = D_1/C_1 | **0.0255713787** | **0.0256** | Table XV |
| A_S = C_1 | **0.88472985** | **0.8846** | Table XV |
| A_D = D_1 | 0.022623762 | — | |
| **Q_d** = (1/20)∫dr r²w(√8 u − w) | **0.270493 fm²** | **0.270 fm²** | Table XV, Eq. (D16) |
| **r_d** = ½[∫dr r²(u²+w²)]^½ | **1.965994 fm** | **1.966 fm** | Table XV, Eq. (D17) |
| B_d | (input) | 2.224575 MeV | Table XV |

A_S and η are read off the parameterisation's own asymptotics: m_1 = γ is the
slowest-decaying mass, so u_a → C_1 e^{−γr} and w_a → D_1 e^{−γr}[1+3/(γr)+3/(γr)²],
which is Eq. (D15) with A_S = C_1, A_D = D_1.

### 5.1 The analytic form against Machleidt's own numerical wave function

Table XIX (70 tabulated r from 0.01 to 14 fm) versus `u_a`, `w_a` built from
the coefficients of §3:

| r [fm] | u tabulated | u_a | Δ | w tabulated | w_a | Δ |
|---|---|---|---|---|---|---|
| 0.01 | 3.040610e−03 | 3.082652e−03 | +4.20e−05 | −1.372760e−06 | −1.168489e−06 | +2.04e−07 |
| 0.30 | 9.938760e−02 | 9.938706e−02 | −5.36e−07 | 3.350710e−03 | 3.322099e−03 | −2.86e−05 |
| 1.00 | 4.310720e−01 | 4.311183e−01 | +4.63e−05 | 1.191650e−01 | 1.191323e−01 | −3.27e−05 |
| 1.70 | 5.201580e−01 | 5.201685e−01 | +1.05e−05 | 1.522870e−01 | 1.522621e−01 | −2.49e−05 |
| 3.00 | 4.312750e−01 | 4.312573e−01 | −1.77e−05 | 9.374530e−02 | 9.376886e−02 | +2.36e−05 |
| 5.00 | 2.771260e−01 | 2.770627e−01 | −6.33e−05 | 3.854870e−02 | 3.860173e−02 | +5.30e−05 |
| 14.00 | 3.459290e−02 | 3.459639e−02 | +3.49e−06 | 1.955820e−03 | 1.955942e−03 | +1.22e−07 |

**max |u_a − u_tab| = 1.355e−04, max |w_a − w_tab| = 7.212e−05** over all 70
points. Machleidt quotes the fit quality as the L² norms
{∫dr[u−u_a]²}^½ = 2.2e−4 and {∫dr[w−w_a]²}^½ = 1.1e−4; the pointwise maxima
above are the right size for those. **This is the check that says the C_j, D_j
were transcribed correctly and the constraints implemented correctly** — a
single wrong digit anywhere would not reproduce a 70-point table to 1e−4.

### 5.2 m_1 and the binding energy

The task asked for `m_1 = sqrt(M_N |E_d|)`. That form is the **non-relativistic
approximation and it does not reproduce γ**; Machleidt's γ carries a
relativistic correction. With his own constants — ħc = 197.327053 MeV fm
(`bonn.tex:2577`), B_d = 2.224575 MeV, M_p = 938.27231, M_n = 939.56563 MeV
(from M̄ = 938.91897 and δM = 1.29332 MeV, `bonn.tex:2719-2722`),
M̂ = 2M_pM_n/(M_p+M_n) = 938.918525 MeV:

| formula | γ [fm⁻¹] | relative deviation from Table's 0.2315380 |
|---|---|---|
| √(M̂ B_d)/ħc — the naive `sqrt(M_N|E_d|)` | 0.2316066 | **+2.961e−04** |
| Eq. (D11), γ² = M̂B_d(1 − B_d/4M̄) | **0.2315380** | −2.005e−07 |
| Eq. (D6)/(D5) exact, γ² = [4M_p²M_n² − (M_d²−M_p²−M_n²)²]/4M_d² | **0.2315380** | −2.010e−07 |

Inverting: γ = 0.2315380 fm⁻¹ implies B_d = **2.223258 MeV** non-relativistically
(0.06 % low) and **2.224576 MeV** through Eq. (D11) — against the published
2.224575 MeV. **m_1 is consistent with the deuteron binding energy to 5e−7,
but only with the relativistic factor.** A header that writes
`m_1 = sqrt(M_N*Ed)/hbarc` would be 3e−4 wrong; write the literal 0.2315380
and cite Eq. (D6).

---

## 6. CD-Bonn against AV18, point by point

Repository reference: `data/vmc/deuteron/fdeut.av18`, k block at line 10221,
201 rows, 0 → 20 fm⁻¹ in steps of 0.1. Its own header:
`ebind 2.224574, dstate 0.057599, qm 0.269673, as 0.885056, eta 0.025045,
rd 1.967364`. Measured on that block: ∫k²(u²+w²)dk = 0.9999764 and
P_D = 0.0575999 — i.e. AV18 is already in the same normalisation, so the two
wave functions are directly comparable with no rescaling.

| deuteron property | AV18 (`fdeut.av18` header) | CD-Bonn (computed, §5) |
|---|---|---|
| B_d [MeV] | 2.224574 | 2.224575 |
| P_D | 0.057599 | **0.048562** |
| Q_d [fm²] | 0.269673 | 0.270493 |
| A_S [fm^−1/2] | 0.885056 | 0.884730 |
| η | 0.025045 | 0.025571 |
| r_d [fm] | 1.967364 | 1.965994 |

**u and w on the shared grid** (CD-Bonn in `fdeut.av18`'s sign convention of
§4, i.e. w = −ψ₂ᵃ > 0 at low p; "proxy" is item 5's, reconstructed in §7):

| p [fm⁻¹] | p [GeV] | u AV18 | u CD-Bonn | CDB/AV18 | w AV18 | w CD-Bonn | CDB/AV18 | w proxy | CDB/proxy |
|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 0.0197 | 1.0693e+01 | 1.0702e+01 | 1.0008 | 4.9364e−02 | 5.0365e−02 | 1.0203 | 4.530e−02 | **1.1119** |
| 0.2 | 0.0395 | 7.1434e+00 | 7.1510e+00 | 1.0011 | 1.3129e−01 | 1.3383e−01 | 1.0193 | 1.205e−01 | 1.1109 |
| 0.3 | 0.0592 | 4.5262e+00 | 4.5337e+00 | 1.0017 | 1.8586e−01 | 1.8919e−01 | 1.0179 | 1.705e−01 | 1.1093 |
| 0.5 | 0.0987 | 1.9580e+00 | 1.9670e+00 | 1.0046 | 2.1980e−01 | 2.2267e−01 | 1.0131 | 2.017e−01 | 1.1040 |
| 0.7 | 0.1381 | 9.5928e−01 | 9.6980e−01 | 1.0110 | 2.0923e−01 | 2.1036e−01 | 1.0054 | 1.920e−01 | 1.0957 |
| 1.0 | 0.1973 | 3.7581e−01 | 3.8790e−01 | 1.0322 | 1.7197e−01 | 1.7004e−01 | 0.9888 | 1.578e−01 | 1.0776 |
| 1.5 | 0.2960 | 8.2175e−02 | 9.4898e−02 | 1.1548 | 1.1226e−01 | 1.0653e−01 | 0.9490 | 1.030e−01 | 1.0342 |
| 2.0 | 0.3947 | 6.0157e−03 | 1.7928e−02 | *2.98* | 7.0965e−02 | 6.3458e−02 | 0.8942 | 6.512e−02 | 0.9745 |
| 3.0 | 0.5920 | −1.6126e−02 | −7.7239e−03 | 0.4790 | 2.7633e−02 | 2.0772e−02 | 0.7517 | 2.536e−02 | 0.8192 |
| 4.0 | 0.7893 | −1.0131e−02 | −5.4036e−03 | 0.5334 | 1.0203e−02 | 5.9334e−03 | 0.5815 | 9.363e−03 | 0.6337 |
| 5.0 | 0.9867 | −4.5165e−03 | −2.3786e−03 | 0.5266 | 3.3634e−03 | 1.2100e−03 | **0.3597** | 3.086e−03 | 0.3920 |

*The 2.98 at p = 2 fm⁻¹ is not a magnitude ratio — it straddles the S node and
should be read as the node moving, not as u being three times bigger.*

**Nodes** (Brent on the analytic form; on AV18, linear interpolation of the
file grid):

| | AV18 | CD-Bonn |
|---|---|---|
| first u node [fm⁻¹] | 2.0929 | **2.3632** (+12.9 %) |
| second u node | 7.6249 | 8.3947 |
| first w node | between 7.3 and 7.4 | 5.9139 |

**The tail, which is the whole point** — fraction of the norm above p:

| p [fm⁻¹] | AV18 | proxy | CD-Bonn | CDB/AV18 |
|---|---|---|---|---|
| 1.0 | 0.08168 | 0.07498 | 0.07435 | 0.910 |
| 1.5 | 0.03666 | 0.03204 | 0.02657 | 0.725 |
| 2.0 | 0.02208 | 0.01936 | 0.01227 | **0.556** |
| 3.0 | 0.00819 | 0.00748 | 0.00298 | **0.364** |

and the kinetic energy ⟨T⟩ = (ħc)²/M̂ ∫dk k⁴(u²+w²), a pure tail diagnostic:

| | ⟨T⟩ [MeV] |
|---|---|
| AV18, file grid | **19.811** — against `fdeut.av18`'s own header `ti = 19.810`, so the estimator is validated |
| proxy | 18.570 |
| CD-Bonn (identical on the 0.1 grid and on a 0.02 grid to 200 fm⁻¹) | **15.605** |

*Caveat: Machleidt's Table XV does not quote ⟨T⟩, so 15.605 MeV is unchecked
against any published CD-Bonn number, and the parameterisation is fitted in
r space (§5.1) with no stated momentum-space fit quality. Only 0.9 % of the
integral comes from above 5 fm⁻¹ (against 2.6 % for AV18), so it is not
tail-dominated, but treat the third digit as unwarranted.*

---

## 7. What the proxy actually did, and by how much it was wrong

`phase_D_gate.md:324-325` describes item 5's proxy as *"AV18's w(k) rescaled so
that P_D = 4.850 %, S renormalised"*. **There is no committed implementation**
(the 2026-09-02 phase-D work was done in scratch), so it was reconstructed from
that sentence:

```
s_w = sqrt(0.04850 / 0.0575999)   = 0.9176140     w_proxy = s_w * w_AV18
s_u = sqrt(0.95150 / 0.9424001)   = 1.0048164     u_proxy = s_u * u_AV18
```

which gives ∫k²(u²+w²)dk = 0.9999764 (unchanged) and P_D = 0.0485000 exactly,
as advertised. §8 confirms the reconstruction is right by reproducing the
proxy's published landmarks to five digits.

**The proxy is a single number multiplying w everywhere.** CD-Bonn's D wave is
not a rescaled AV18 D wave — the column "CDB/proxy" in §6 runs from **1.11 at
p = 0.1 fm⁻¹ down to 0.39 at p = 5 fm⁻¹**, so the proxy is ~11 % too *small*
where the D wave is largest and 2.6× too *large* in the tail. It got both
directions wrong at once:

| ∫dk k² w² over | AV18 | proxy | CD-Bonn | CDB/AV18 | **CDB/proxy** |
|---|---|---|---|---|---|
| 0.0 – 0.5 | 1.704237e−03 | 1.434994e−03 | 1.757048e−03 | 1.0310 | **1.2244** |
| 0.5 – 1.0 | 1.121565e−02 | 9.443752e−03 | 1.122931e−02 | 1.0012 | 1.1891 |
| 1.0 – 1.5 | 1.501409e−02 | 1.264209e−02 | 1.413619e−02 | 0.9415 | 1.1182 |
| 1.5 – 2.0 | 1.219912e−02 | 1.027185e−02 | 1.042805e−02 | 0.8548 | 1.0152 |
| 2.0 – 3.0 | 1.275711e−02 | 1.074168e−02 | 8.985620e−03 | 0.7044 | 0.8365 |
| 3.0 – 5.0 | 4.582433e−03 | 3.858479e−03 | 2.009601e−03 | 0.4385 | 0.5208 |
| 5.0 – 20.0 | 1.258962e−04 | 1.060066e−04 | 1.623857e−05 | **0.1290** | 0.1532 |

and the S–D interference, which is what actually drives b₁:

| ∫dk k² u w over | AV18 | proxy | CD-Bonn | CDB/AV18 |
|---|---|---|---|---|
| 0.0 – 0.5 | 2.714650e−02 | 2.502998e−02 | 2.766283e−02 | 1.0190 |
| 0.5 – 1.0 | 4.544295e−02 | 4.189993e−02 | 4.622737e−02 | 1.0173 |
| 1.0 – 1.5 | 2.022717e−02 | 1.865014e−02 | 2.101947e−02 | 1.0392 |
| 1.5 – 2.0 | 4.838986e−03 | 4.461708e−03 | 6.064740e−03 | **1.2533** |
| 2.0 – 3.0 | −2.911346e−03 | −2.684359e−03 | **+5.418111e−06** | ≈ 0 |
| 3.0 – 5.0 | −3.738589e−03 | −3.447105e−03 | −1.236020e−03 | 0.3306 |
| 5.0 – 20.0 | −1.883304e−04 | −1.736469e−04 | −1.783146e−05 | 0.0947 |

**This is the mechanism of §8's result.** CD-Bonn's u node sits 13 % higher
(2.363 vs 2.093 fm⁻¹) and its u is half as big beyond it, so the *negative*
lobe of the S–D interference — which in AV18 partially cancels the positive
one — is almost entirely gone (−2.9e−3 → +5.4e−6 in the 2–3 fm⁻¹ bin). The
proxy, which only rescales, keeps AV18's node and keeps the full negative lobe;
it can never see this. That is precisely the failure `phase_D_gate.md:331-333`
predicted it would have and could not quantify.

---

## 8. The measured effect on the A = 2 gate

**How.** `DeuteronConvolutionB1::Options::fdeut_path` is bound to Python
(`python/bindings.cpp:1449`), so CD-Bonn and the proxy were fed to the gate by
writing them as `fdeut`-format files **in the scratchpad**, with no repository
file touched. The curve, the grid and the landmarks are
`validation/b1_li6_table.py`'s own — `x·b1(x, 2.5, 0)` on
`sorted(x_grid(0.01,1.0,100) | x_grid(0.10,0.80,100))`, `lg.b1_landmarks`,
G3b = peak/|table peak| — the script was imported, not modified.

**Four controls, all passed, before any new number was believed:**

1. Rewriting `fdeut.av18`'s own columns through my writer and reading them back
   reproduces the repository file's landmarks to **every printed digit**.
2. AV18 + ToyF2 + κ = 1 gives zeros (0.02231, 0.39230), min −3.12114e−05 @ 0.2343,
   max 2.94836e−04 @ 0.7400, G3b 0.2717 — `phase_D_gate.md`'s own row is
   *0.39230 / −3.1211e−5 @ 0.2343 / 2.9484e−4 @ 0.7400 / 0.2717*. Exact.
3. AV18 + CT18NLO + default κ gives zeros (0.43877), min −1.86794e−04 @ 0.3000,
   max 7.80616e−04 @ 0.7364, G3b 0.7193 — the doc's row is
   *0.43877 / −1.8679e−4 @ 0.3000 / 7.8062e−4 @ 0.7364 / 0.7193*. Exact.
   And the proxy reproduces item 5's *2.5855e−4* and *0.40952* as
   **2.58546e−04 @ 0.7400** and **0.40954** — so §7's reconstruction of the
   proxy is the proxy.

Digitized CDKS reference (`lg.b1_landmarks_of_table()`): zeros 0.06564 and
0.45718, min −1.76814e−04 @ 0.3324, max **1.08521e−03 @ 0.7657**.

### 8.1 ToyF2 (the repository's stand-in F₂)

| wave function | κ | zeros in [0.01, 1] | min | max | ∫b₁dx | **G3b** |
|---|---|---|---|---|---|---|
| AV18 | 1 | 0.02231, 0.39230 | −3.12114e−05 @ 0.2343 | 2.94836e−04 @ 0.7400 | +1.06203e−04 | 0.2717 |
| proxy | 1 | 0.01333, 0.40954 | −3.31084e−05 @ 0.2556 | 2.58561e−04 @ 0.7400 | +7.56451e−05 | 0.2382 |
| **CD-Bonn** | 1 | **none** | +7.86463e−06 @ 0.1354 | 3.57271e−04 @ 0.7293 | +3.00741e−04 | **0.3292** |
| AV18 | Eq. 21 | 0.02228, 0.37737 | −3.30889e−05 @ 0.2343 | 4.77515e−04 @ 0.7576 | +2.07311e−04 | 0.4400 |
| proxy | Eq. 21 | 0.01332, 0.39372 | −3.52817e−05 @ 0.2556 | 4.19941e−04 @ 0.7576 | +1.63073e−04 | 0.3870 |
| **CD-Bonn** | Eq. 21 | **none** | +8.04692e−06 @ 0.1354 | 5.73711e−04 @ 0.7505 | +4.33992e−04 | **0.5287** |

The AV18 / Eq. 21 row is control 4: `phase_D_gate.md`'s own entry is
*0.37737 / −3.3089e−5 @ 0.2343 / 4.7751e−4 @ 0.7576 / 0.4400*. Exact.

### 8.2 CT18NLO — the row that matters, because a real PDF is the known dominant systematic

| wave function | κ | zeros in [0.01, 1] | min | max | ∫b₁dx | **G3b peak ratio** |
|---|---|---|---|---|---|---|
| AV18 | 1 | 0.44920 | −1.70253e−04 @ 0.2980 | 4.91781e−04 @ 0.7152 | +5.35877e−05 | 0.4532 |
| proxy | 1 | 0.45656 | −1.58894e−04 @ 0.3051 | 4.31128e−04 @ 0.7200 | +3.20904e−05 | 0.3972 |
| **CD-Bonn** | 1 | **0.04340, 0.39824** | −1.21512e−04 @ 0.2697 | 5.96736e−04 @ 0.7010 | +2.16730e−04 | **0.5499** |
| AV18 | **Eq. 21 (default)** | 0.43877 | −1.86794e−04 @ 0.3000 | 7.80616e−04 @ 0.7364 | +2.01428e−04 | **0.7193** |
| proxy | Eq. 21 | 0.44557 | −1.74823e−04 @ 0.3051 | 6.86489e−04 @ 0.7364 | +1.59986e−04 | **0.6326** |
| **CD-Bonn** | **Eq. 21** | **0.04336, 0.39191** | −1.31611e−04 @ 0.2700 | **9.37865e−04 @ 0.7293** | +4.11024e−04 | **0.8642** |
| *digitized CDKS* | — | *0.06564, 0.45718* | *−1.76814e−04 @ 0.3324* | *1.08521e−03 @ 0.7657* | — | *1* |

**Three findings, in order of consequence.**

1. **G3b: the deficit halves.** The peak ratio goes 0.7193 → **0.8642**, i.e. a
   factor 1.39 low becomes a factor **1.157** low. `phase_D_gate.md` closed
   with CT18 as "the only identified item still open at the peak"; on this
   measurement the wave function was a second one of comparable size, and in
   the opposite direction to what the proxy said.

2. **G3a: the counting clause un-regresses.** `phase_D_gate.md:317-321` records
   that with CT18 "the first sign change moves below x = 0.01 (the scan floor),
   so 'exactly two sign changes in [0.02, 1.0]' becomes one and G3a's counting
   clause fails." With CD-Bonn the low-x zero comes back **inside the window,
   at x = 0.04336** (against CDKS's own 0.06564, and well inside the design's
   stated tolerance of ±0.08 about 0.0656). The `PLAN.md` A2 item —
   "G3a will abort before G3b is reached … widening the counting window" —
   **may not be needed at all** once A5 lands, and A2 should be re-measured
   after A5 rather than before it.

3. **The mid-x dip gets worse, and this is not hidden.** −1.868e−04 @ 0.300
   (AV18) → −1.316e−04 @ 0.270 (CD-Bonn), against CDKS's −1.768e−04 @ 0.3324.
   AV18 was **closer** on the dip. CD-Bonn improves the peak, the zero count
   and ∫b₁dx and degrades the dip; it is not a uniform improvement and must
   not be reported as one.

*Caveat carried from the existing docs, undiminished: CDKS used **MSTW2008 LO**,
not CT18NLO. Every row above is a stand-in. A5's number and A1/A2's number are
not independent, and the honest statement of the residual after A5 needs MSTW
in hand.*

### 8.3 Sign diagnostics (ToyF2, κ = 1) — why §4's sign is not a matter of taste

| wave function | zeros | min | max | G3b |
|---|---|---|---|---|
| AV18 (correct) | 0.02231, 0.39230 | −3.12e−05 @ 0.2343 | 2.94836e−04 @ 0.7400 | 0.2717 |
| AV18, w → −w *(wrong, diagnostic)* | 0.81759, 0.91005 | −5.85e−06 @ 0.8600 | 2.15647e−04 @ 0.4500 | 0.1987 |
| CD-Bonn (correct) | none | +7.86e−06 @ 0.1354 | 3.57271e−04 @ 0.7293 | 0.3292 |
| CD-Bonn, w → −w *(wrong, diagnostic)* | 0.55510 | −1.17555e−04 @ 0.7788 | 5.38743e−05 @ 0.3828 | 0.1083 |

Either wrong-sign curve is unrecognisable as CDKS Fig. 4. Combined with the
quadrupole gate of §4, the sign is settled three independent ways.

### 8.4 Robustness

The CD-Bonn table was fed to the gate on **`fdeut.av18`'s own 0.1 fm⁻¹ grid**
and on a **0.02 fm⁻¹ grid** (5× finer, 1001 rows). ToyF2, κ = 1:
peak 3.57271e−04 vs 3.57443e−04, min +7.86463e−06 vs +7.86853e−06, G3b 0.3292
vs 0.3294 — **5e−4 relative**, so the coarse-grid spline that
`b1_nuclear.hpp:176` warns about is not what is driving any of this. The gate's
own `LightConeDensities` reports `p_d_momentum` = 0.0485631 for the CD-Bonn
file against the closed-form 0.0485621 (2e−6), confirming the file is read as
intended.

---

## 9. What the A5 code step should do

1. A new header (`include/lipolgen/cdbonn.hpp` or a section of `cluster.hpp`)
   carrying **only** the twenty published numbers of §3.1, γ, m₀ and n = 11.
   C₁₁ and D₉₋₁₁ are **computed** from Eqs. (D23)/(D24), never typed — they are
   not published numbers, and typing them would be the second definition of a
   physics quantity that `docs/CONVENTIONS.md` forbids.
2. Return the pair in **`fdeut.av18`'s convention** (§4): `u = +ψ₀ᵃ`,
   `w = −ψ₂ᵃ`, normalised as ∫dp p²(u²+w²) = 1. Say so in the header, and
   point at §4's measurement rather than at Eq. (D13), which is misleading as
   printed.
3. Feed the gate through a **new option**, not by changing
   `DeuteronConvolutionB1::Options::fdeut_path`'s default. Item 5's rescaling
   proxy stays reachable, as `PLAN.md` A5 requires.
4. Units: the header should be in fm⁻¹ / fm^(3/2) — the paper's — and convert
   with `HBARC_GEV_FM` at the one boundary, as `read_fdeut_k` already does.
   Do not introduce a third ħc (`cluster.cpp:23`'s file-local copy is already
   one too many).
5. Tests to pin, all computed above and all cheap and closed-form:
   norm = 0.999999826159561, P_D = 0.048562073645887, η = 0.0255713787,
   A_S = 0.88472985, `alpha_d_quadrupole_fm2` = +0.270178 fm², and the four
   constraint residuals of §3.1. A pytest on the new option, per the ground rules.
6. **Re-run A2 after A5, not before.** §8.2 finding 2.

---

## 10. What is verified and what is not

**Verified here, by computation:** the coefficients reproduce Machleidt's own
Table XIX to 1.4e−4 over 70 points (§5.1); normalisation, P_D, η, A_S, Q_d and
r_d all reproduce Table XV (§5); the four boundary constraints hold to 1.8e−12
(§3.1); the closed-form and linear-system routes to D₉₋₁₁ agree to 1e−14 (§2);
m₁ = γ is consistent with B_d to 5e−7 through Eq. (D11) (§5.2); the sign
convention is fixed three independent ways (§4, §8.3); the gate reproduction
reproduces three published rows exactly (§8).

**Not verified, and not assumable:**

* **The published-journal table numbering.** §1 — I read the e-print, not
  PRC 63 024001. The plan's "Tables XVII/XVIII" is wrong for the e-print;
  whether it is right for the journal I cannot say.
* **The parameterisation's momentum-space fidelity above ≈ 3 fm⁻¹.** Machleidt
  quotes fit quality only as an r-space L² norm. The analytic form has
  u ∼ p⁻⁴, w ∼ p⁻⁶ tails (measured: p⁴u → 0.349 at p = 100 against the
  predicted −√(2/π)ΣC_jm_j² = 0.362; p⁶w → −667 against 710, still converging),
  which is a property of the ansatz, not of CD-Bonn. §6's ⟨T⟩ = 15.605 MeV has
  nothing published to check it against. **If the b₁ result of §8 turns out to
  hinge on p > 3 fm⁻¹, ask Machleidt for the momentum-space data file** — the
  paper says (Appendix D, last line) they are available from the author.
* **MSTW2008 LO.** Every §8 row uses ToyF2 or CT18NLO. A5 and A1/A2 interact.
* **Whether CD-Bonn is the right thing to hold fixed.** CDKS used CD-Bonn *and*
  MSTW *and* their own kernel. Swapping in one of the three does not make this
  a reproduction of their calculation.

---

## 11. Scratch files (not in the repository)

Under `/tmp/claude-1000/-home-cpeng-Projects-polli/3ac6e55b-3bf7-4951-8170-ad408d0b635b/scratchpad/`:

| file | what |
|---|---|
| `cdbonn/eprint.tar.gz`, `cdbonn/bonn.tex` | the fetched e-print source |
| `cdbonn.pdf`, `cdbonn.txt` | the PDF and its `pdftotext -layout` rendering |
| `cdbonn_coef.py` | the coefficients + constraints + ψ₀, ψ₂ (the reference implementation) |
| `cdbonn_rwaves.json` | Table XIX parsed out of `bonn.tex` |
| `rcheck.py` | 60-digit r-space check against Table XIX; Q_d, r_d, P_D |
| `ftcheck.py`, `conv.py` | the §4 Fourier-transform sign measurements |
| `compare.py`, `diag.py`, `diag2.py` | §5.2, §6, §7 |
| `gate2.py`, `gate3.py`, `gate4.py` (+ `.out`) | §8 |
| `cdbonn_uw_fixture.txt` | u(p), w(p) at 17 p points, 15 digits, in the repo's sign convention — a ready-made regression fixture |

**Baseline.** No repository file was created, edited or deleted other than this
one, so nothing can have regressed. The suite was re-run anyway, after all the
measurements above:

```
[doctest] test cases:      363 |      363 passed | 0 failed | 1 skipped
[doctest] assertions: 17203863 | 17203863 passed | 0 failed |
211 passed in 34.28s
```

— the required baseline, unmoved.
