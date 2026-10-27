<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-26 — phase D: the "least verified inputs" read against their papers

`PLAN.md` phase D (D1, D2), the two items `../run_2026-09-26/SUMMARY.md`
"What this run did NOT establish" left open (3: the T1 references are a
second-hand vendoring; 4: η_d, Q_d and α_s(M_Z) are the least verified
inputs). Every paper behind them was read on 2026-10-26, page by page; the
copies stay outside the repository (§1.1, §2), only numbers are vendored.

**Outcome in one paragraph.** The g₁ᵈ D-state convention **held**: COMPASS,
HERMES, E143 and SMC all publish g₁ᵈ per nucleon and **not** D-state
corrected, and apply (1 − 1.5 ω_D) only to extract the neutron or the nucleon,
so `t1_g1d_world`'s verdict column (`g1d_tree`, with the factor) stands and
nothing was swapped. The E143 arXiv identifier is **hep-ph/9802357** (the
NNPDF metadata's hep-ex/9705012 is E154's g₁ⁿ letter); the SMC set is
**PRD 58 (1998) 112001**, Table IX (still recorded only). Every vendored T1
point equals its paper's table (COMPASS 11 of 15, HERMES 15/15, E143 28/28,
SMC 13/13, NMC F₂ᵈ 158/158, NMC F₂ᵈ/F₂ᵖ 260/260) except **four high-x COMPASS
points, where HEPData differs from the printed table and HEPData is the
self-consistent one** (§1.3). NMC's 2.5 % normalisation is **not** in its F₂ᵈ
'sys' and nowhere in the file; NMC extracted F₂ᵈ with its own measured R below
x = 0.12 and R1990 above; the ratio assumes R_d = R_p. η_d = 0.0256(4) is
Rodning–Knutson's measured conclusion, transcribed right. **Q_d was
transcribed wrong**: Bishop–Cheung conclude 0.2860 ± 0.0015 fm² (a calculation
confirming Reid–Vaida's semi-empirical value); 0.2859(3) does not occur in the
paper — it is what the AV18 and CD-Bonn papers attribute to it. The vendored
JSON is corrected with a dated note; **no verdict moves** (Q_d is recorded,
not asserted). α_s(M_Z) = 0.1180(9) is PDG 2026's value (Table 1.1 and the QCD
review, Eq. (9.25)); nothing moved. Suelzle's ⁷Li HO fit is a₀ = 1.755 ± 0.02 F
(the survey's (1.77(2), 0.327) is de Vries's re-expression of it); ⁶Li has no
HO fit in that paper; Moniz's ⁶Li k_F = 169 ± 5 MeV/c is the tree's
`RC_QE_KF_GEV` (recorded only).

## 1. D1 — the T1 references

### 1.1 The copies read

| paper | identifier | copy (never in the repository) | sha256 | read |
|---|---|---|---|---|
| COMPASS final g₁ᵈ | arXiv:1612.00620**v1** (the only version) | sibling `refs/1612.00620.pdf` | `851679b1…2398c4de` | pp. 5–6, 8, 10, 14–16 of the PDF |
| the same, as published | PLB 769 (2017) 34–41, SCOAP3 CC BY 4.0 | INSPIRE's file store (`inspirehep.net/files/28462197…`), fetched with curl's own identity, scratch `refs/compass2017_PLB769_34_scoap3.pdf` | `de96c817…d8dd0751` | pp. 36, 38–41 |
| HERMES g₁ | hep-ex/0609039, PRD 75 (2007) 012007 | sibling `refs/hep-ex_0609039.pdf` | `c83001c2…d0ee1f90` | PDF pp. 6, 14–15, 22, 45 |
| E143 g₁, g₂ | **hep-ph/9802357**, SLAC-PUB-7753, PRD 58 (1998) 112003 | sibling `refs/hep-ph_9802357.pdf` | `ae0b7a94…2a010633ec` | PDF pp. 1–2, 7, 14, 37–40, 46, 49–51, 75–76 |
| SMC A₁, g₁ | PRD 58 (1998) 112001, CERN-EP-98-085 (no arXiv eprint) | APS, curl's own identity (the institution's licence), scratch `refs/smc1998_PhysRevD.58.112001.pdf` | `b6b4ca7b…61eb195` | pp. 112001-3, -6, -12, -13, -15, -16 |
| NMC F₂ᵖ, F₂ᵈ, R | hep-ph/9610231, NPB 483 (1997) 3 | sibling `refs/hep-ph_9610231.pdf` | `ebbccbce…d9f41b01` | PDF pp. 6–10, 14, 18–19, 25, 33–36 |
| NMC F₂ᵈ/F₂ᵖ, R_d − R_p | hep-ex/9611022, NPB 487 (1997) 3 | sibling `refs/hep-ex_9611022.pdf` | `54f64f81…544db76c6` | PDF pp. 4, 6, 8, 10–13, 15–20 |
| NNPDFpol1.1 | 1406.5539, NPB 887 (2014) 276 | sibling `refs/1406.5539.pdf` | `9a8de21b…38cf69a3` | searched for its deuteron convention: none stated there |

Identifiers were resolved through INSPIRE's API (standard `requests` client):
record 467140 = "Measurements of the proton and deuteron spin structure
functions g1 and g2", E143, arXiv hep-ph/9802357, DOI
10.1103/PhysRevD.58.112003; record 471981 = SMC, "Spin asymmetries A1 and
structure functions g1 of the proton and the deuteron from polarized
high-energy muon scattering", PRD 58 (1998) 112001, no arXiv eprint; arXiv
hep-ex/9705012 = INSPIRE 443170, E154, "Precision determination of the
neutron spin structure function g1n", PRL 79 (1997) 26. HEPData
`ins1501480` (COMPASS) was read through its JSON API: version 1, last updated
17 July 2017, table 1 "Values of A1d and g1d … in x bins averaged over Q²".

### 1.2 The rule, and the D-state convention the papers state

**The rule was declared before any paper was read**: `PLAN.md` phase D1
(committed at `87787dd`) — "declare, before reading, that the verdict follows
the convention the papers state". This worker wrote it into the harness header
after the reading (the order of the edit, not of the decision; the header says
so). The verdict column would have been swapped to `dis_no_dstate` if the
papers divided the factor out.

| experiment | what the paper states | where |
|---|---|---|
| COMPASS | g₁ᵈ = F₂ᵈ A₁ᵈ / (2x(1 + R)), Eq. (7), F₂ᵈ from SMC [2], R from E143 R1998 [20]; the factor enters only the **nucleon** moment Γ₁ᴺ = ∫ g₁ᵈ/(1 − 1.5 ω_D) dx with **ω_D = 0.05 ± 0.02** [Machleidt, Phys. Rep. 149] | PLB pp. 38, 39 (arXiv v1 PDF p. 10 prints "ω_D = −0.05 ± 0.01" — a misprint the journal corrects) |
| HERMES | g₁ᵈ = ½(g₁ᵖ + g₁ⁿ)(1 − 1.5 ω_D), Eq. (23), **ω_D = 0.05 ± 0.01**, "In this paper the neutron structure function g₁ⁿ is evaluated according to Eq. (23)"; ω_D enters the systematics of g₁ⁿ and g₁ᴺˢ only | arXiv PDF pp. 6, 15 |
| E143 | g₁ⁿ = 2 g₁ᵈ/(1 − 1.5 ω_D) − g₁ᵖ, Eq. (56), **ω_D = 0.05 ± 0.01** ("given by N-N potential calculations"); Γ₁ᵈ carries (1 − 1.5 ω_D), Eq. (11) | arXiv PDF pp. 37, 7 |
| SMC | g₁ᴺˢ = 2[g₁ᵖ − g₁ᵈ/(1 − 1.5 ω_D)], Eq. (13), **ω_D = 0.05 ± 0.01** "as in our previous publications" | p. 112001-12 |

**All four publish the deuteron's g₁ per nucleon, uncorrected.** The tree's
`g1_nucleus(deuteron)/2 = (1 − 1.5 P_D)(g₁ᵖ + g₁ⁿ)/2` is therefore the
like-for-like comparison; **the verdict column stands; no column was
swapped; no test pin moved.** One number is not like for like and is stated,
not corrected: the experiments' ω_D = 0.05 against the tree's
`P_D_DEUTERON = 0.045` (factors 0.925 and 0.9325, 0.8 % of g₁ᵈ) — a property of
the tree's deuteron model on the theory side; ω_D does not enter the published
g₁ᵈ.

**What phase 0 had previewed, and what it means now.** NNPDFpol1.1 fails on
this column — COMPASS 45.0/15 (p = 7.7e−5), HERMES 52.6/15 (p = 4.5e−6), E143
51.0/28 (p = 0.0050), DIS cut Q² ≥ 1, W² ≥ 4 GeV², diagonal stat ⊕ sys — and
the no-D-state reading is **worse** on all three (56.4, 67.5, 61.9;
`phase_A_baseline.md` §3). With the convention now confirmed, the D-state is
**cleared** as the cause. The F₁ division phase 0 also named concerns only the
ToyG1 configurations (toy, mstw, ct18nlo), not `LhapdfG1`. What remains for
NNPDFpol1.1 is the reading phase A §3 already offered: `LhapdfG1` is the
leading-order g₁ of an NLO set (`include/lipolgen/sf.hpp`, "LO g1 from a
polarized LHAPDF grid"; no ΔC_q ⊗ Δq, no ΔC_g ⊗ Δg) — **a reading, not tested
here.** The NNPDFpol1.1 paper on disk (1406.5539) states no deuteron convention
of its own (it refers its DIS treatment to NNPDFpol1.0, not on disk).

### 1.3 Every vendored point against its paper's table

Each comparison parsed the paper's table from `pdftotext -layout` and compared
every point of the vendored file in x, Q², value and each error component
(the scripts are in the session scratch directory; the outcomes, which anyone
can repeat from the copies in §1.1, are):

| set | the paper's table | points equal to the file | the file's components, identified |
|---|---|---|---|
| HERMES | Table XXI (PDF p. 45: 15 x-bins after a Q² > 1 GeV² cut, each the average of its Q² bins evolved to its ⟨Q²⟩, PDF p. 14) | **15/15** | `exp` = the paper's syst, **which includes the 5 % deuteron normalisation** (caption); `param` = par; `evol` = evol; the 15 artificial vectors give the paper's stat in quadrature to 5e−11 |
| E143 | Table XIV (PDF p. 75) | **28/28** | `stat`, `sys` = the table's; `sys_beam` = 4.9 % × value (signed) = **Table XV's total deuteron normalisation** (beam polarisation 2.4, target polarisation 4.0, dilution 1.4 %), which the table's syst excludes ("Not included in the systematic errors listed are normalization uncertainties shown in Table XV"); NNPDF's label "beam Normalization" names one of its three parts |
| SMC | Table IX, g₁ᵈ at the measured Q² (p. 112001-13) | **13/13** | `stat`, `sys` = the table's |
| COMPASS | Table A.1 (PLB p. 40; arXiv v1 PDF p. 14) | **11/15**: equal for x ≤ 0.171, **different at the four x ≥ 0.222** | `stat`, `sys` = the table's where they agree |
| NMC F₂ᵈ | Table 6 (arXiv PDF pp. 35–36: F₂ᵈ averaged over the four energies and both triggers, the Table 2 normalisation shifts applied) | **158/158** | `stat`, `sys` = the table's ΔF₂^stat, ΔF₂^syst |
| NMC F₂ᵈ/F₂ᵖ | Table 3 (arXiv PDF pp. 15–20) | **260/260** | `hepdata`: stat, sys = the table's; `legacy` `sys_corr_1..5` = the table's columns 6–10 **VX, SM, RC, E, E′** in % × value, signs kept, at all 260 points (max deviation 1.7e−18) |

**The COMPASS high-x discrepancy.** At x = 0.222 / 0.290 / 0.405 / 0.567 the
paper prints g₁ᵈ = 0.0744 / 0.076 / 0.0576 / 0.0254 (both arXiv v1 and the
journal), the HEPData record — and so nnpdf-data and the vendored file —
carries 0.0930 / 0.0675 / 0.0500 / 0.0203, the errors scaled alike; the A₁ᵈ
column agrees (one stat digit apart at x = 0.222). The paper's own Eq. (7)
decides which is self-consistent: on the paper's A₁ᵈ, with SMC's F₂ᵈ (PRD 58
112001 Eq. (A1), Table XIII, p. 112001-16; signs read from the page image) and
the tree's `lg.r1998` for E143's R1998:

```python
import math
import lipolgen as lg
# SMC F2d, PRD 58 (1998) 112001 Eq. (A1), Table XIII central values (p. 112001-16)
a = [-0.28151, 1.0115, 0.08415, -0.72973, 2.8647, -2.5328, 0.47477]
b = [0.20040, -2.5154, 0.02599, 0.01858]
c = [-1.3569, 7.8938, -29.117, 37.657]
def f2d(x, q2, q02=20.0, lam2=0.25 ** 2):
    A = x ** a[0] * (1 - x) ** a[1] * sum(a[2 + k] * (1 - x) ** k for k in range(5))
    B = b[0] + b[1] * x + b[2] / (x + b[3])
    C = sum(c[k] * x ** (k + 1) for k in range(4))
    return A * (math.log(q2 / lam2) / math.log(q02 / lam2)) ** B * (1 + C / q2)
# COMPASS PLB 769 (2017) 34, Table A.1 (p. 40): <x>, <Q2>, A1d, printed g1d; HEPData ins1501480 t1 g1d
rows = [(0.0766, 9.00, 0.0527, 0.111, 0.111), (0.121, 13.5, 0.095, 0.123, 0.123),
        (0.171, 18.6, 0.121, 0.101, 0.101), (0.222, 23.8, 0.160, 0.0744, 0.0930),
        (0.290, 31.1, 0.190, 0.076, 0.0675), (0.405, 43.9, 0.317, 0.0576, 0.0500),
        (0.567, 60.8, 0.494, 0.0254, 0.0203)]
for x, q2, a1, printed, hepdata in rows:
    g = f2d(x, q2) * a1 / (2 * x * (1 + lg.r1998(x, q2)))      # Eq. (7), PLB p. 38
    print("x = %.4f  Eq.(7) %.4f  printed/Eq.(7) %.3f  HEPData/Eq.(7) %.3f"
          % (x, g, printed / g, hepdata / g))
```

| x | Eq. (7) | printed / Eq. (7) | HEPData / Eq. (7) |
|---|---|---|---|
| 0.0766 | 0.1108 | 1.002 | 1.002 |
| 0.121 | 0.1229 | 1.001 | 1.001 |
| 0.171 | 0.1025 | 0.986 | 0.986 |
| **0.222** | 0.0927 | **0.803** | 1.003 |
| **0.290** | 0.0682 | **1.114** | 0.990 |
| **0.405** | 0.0500 | **1.152** | 1.000 |
| **0.567** | 0.0204 | **1.244** | 0.994 |

HEPData's values follow from the paper's A₁ᵈ by the paper's Eq. (7) to ≤ 1.4 %
at all seven bins; the printed ones miss by −20 / +11 / +15 / +24 % at the four.
**The vendored (HEPData) values are the self-consistent ones and are kept** —
recorded, not re-vendored; the effect of using the printed values instead on
any χ² is not computed (no second vendored table exists, and none is made
here). The cause is not known; the author may wish to tell the COMPASS
contacts or HEPData.

### 1.4 Identifiers and beam energies

* **E143**: the paper is PRD 58 (1998) 112003 = **hep-ph/9802357** =
  SLAC-PUB-7753 = INSPIRE 467140 = HEPData ins467140 (the survey's identifier
  was right). nnpdf-data's metadata field `arXiv: hep-ex/9705012` is **wrong**:
  that identifier is E154's g₁ⁿ letter. The vendored file keeps NNPDF's field
  verbatim (it is regenerated byte for byte from the pinned wheel); the
  harness label now names the paper read. **Beam energies**: 29.13 GeV (122 M
  events), 16.18 GeV (56 M) and 9.71 GeV (58 M) (PDF p. 14); Table XIV is g₁
  at the averaged measured Q² from g₁/F₁ averaged over all three energies and
  both spectrometers (4.5°, 7°) for Q² ≥ 1 (PDF pp. 38, 40), with F₂ from NMC
  PLB 364 (1995) 107 and R from R1990 (PDF p. 37).
* **SMC**: B. Adeva *et al.*, PRD 58 (1998) 112001 (COMPASS's reference [2]).
  Its deuteron data are on **deuterated butanol C₄D₉OD** (Table I,
  p. 112001-3), not ⁶LiD — the `t1_g1d_world` header said ⁶LiD for both COMPASS
  and SMC and is corrected. R for x < 0.12 from NMC, above from R1990
  (p. 112001-6); F₂ᵈ from SMC's own fit (Appendix). SMC **stays recorded
  only**: the survey has no SMC row, and adding it to the verdict now would
  change a declared rule after its χ² was seen (toy 16.7/12).
* **HERMES** and **COMPASS**: identifiers as the harness had them; COMPASS's
  only arXiv version is v1.

### 1.5 NMC: which R, the normalisation, the nuclear effects

* **F₂ᵈ (NPB 483)**: for x < 0.12, NMC's **own** R, one Q²-averaged value per x
  bin determined with F₂ (Table 7, PDF p. 14: R = 0.537 at x = 0.0045 … 0.043
  at x = 0.11); for x > 0.12, **R1990** (Sect. 4, PDF p. 7). Table 6's R column
  equals its R1990 column at every x > 0.12 row. Table 6 also prints F̃₂ᵈ,
  extracted with R1990 throughout: on the DIS-cut points it is **6.3–9.1 %
  below** F₂ᵈ at x = 0.0045, 0.9–3.8 % below at x = 0.008, and within
  −1.4 … +2.8 % for x ≥ 0.0125 — the size of the R choice in the data. The tree's
  F₂ is R-independent, so `t1_nmc_f2d` inherits NMC's R through the data alone.
* **The normalisation (NPB 483)**: "The total systematic uncertainties on F2,
  obtained by adding in quadrature all the above contributions apart from the
  normalisation uncertainty, range from 1% to 5%"; the per-energy 2 % "is
  included in the 2.5% total normalisation uncertainty of the combined data"
  (Sect. 4.1, PDF p. 8; Figs. 5–6 captions, PDF pp. 18–19). **The 2.5 % is not
  in 'sys' and nowhere in the vendored file** (its definitions: `stat` UNCORR
  ADD, `sys` CORR MULT). Neither χ² carries it; it was not added (a new error
  model after the numbers would be a re-tuned tolerance). For scale: a common
  2.5 % shift cannot absorb CT18NLO's mean excess of 8.6 % (3.4 normalisation
  σ with every other residual zero) nor the toy's 6.0 %; MSTW's 2.5 % is one σ,
  and what it would do to MSTW's 416.5/155 (worst pull +5.0) is not computed.
* **F₂ᵈ/F₂ᵖ (NPB 487)**: "structure functions are always given per nucleon"
  (PDF p. 4). **No R enters the ratio**: "As ΔR is compatible with zero we have
  taken the structure function ratio, F2d/F2p, to be equal to the cross
  section ratio, σd/σp" (Sect. 6, PDF p. 11) — R_d = R_p; an absolute R (R1990,
  ref. [11]) enters only the radiative corrections, beside the NMC F₂ᵖ fit PLB
  364 (1995) 107 (Sect. 3.1, PDF p. 6). The ratio's own normalisation error is
  < 0.1 % (target densities and lengths, Sect. 3.3, PDF p. 8), with no
  additional normalisation at the 0.15 % level (Sect. 5, PDF p. 10) — not in
  'sys', negligible against the 0.4 % average systematic (PDF p. 11). The
  paper calls the radiative, vertex, smearing and functional-form errors fully
  correlated between its eleven data sets (PDF p. 8); NNPDF's CORR across all
  points is NNPDF's extension of that.
* **Nuclear effects**: neither paper corrects F₂ᵈ for them. NPB 483's analysis
  (Sects. 3–4, PDF pp. 6–8) applies only the higher-order electroweak
  correction (whose deuteron quasi-elastic suppression concerns the radiative
  tail); NPB 487 forms F₂ⁿ/F₂ᵖ = 2F₂ᵈ/F₂ᵖ − 1 "Neglecting nuclear effects in the
  deuteron" and discusses shadowing and the EMC effect as possible deviations
  (Sect. 6, PDF pp. 11–12).

### 1.6 The T1 numbers (unchanged) and what they mean

Re-run 2026-10-26 after the header edits, bit for bit the phase 0 / phase A
numbers (DIS cut Q² ≥ 1 GeV², W² ≥ 4 GeV², and ≥ the backend's grid floor —
CT18NLO 1.677 GeV²; diagonal stat ⊕ sys; tolerance p(χ², ndf) ≥ 0.01, unchanged):

| row | toy (verdict) | mstw | ct18nlo | nnpdfpol |
|---|---|---|---|---|
| `t1_g1d_world` COMPASS / HERMES / E143 | 107.8/15, 69.4/15, 54.9/28 — **FAIL** | 53.2, 47.1, 50.2 — FAIL | 20.2/11 (p = 0.043), 64.3/11, 59.0/24 — FAIL | 45.0/15, 52.6/15, 51.0/28 — FAIL (no-D-state 56.4, 67.5, 61.9) |
| `t1_nmc_f2d` | 8307.9/155, ⟨th/data⟩ 1.060 — **FAIL** | 416.5/155, 1.025 — FAIL | 3401.8/146, 1.086 — FAIL | — |
| `t1_nmc_f2d_over_f2p` | 2161.4/211 — **FAIL** | 247.2/211 (p = 0.044) — PASS | 215.5/189 (p = 0.090) — PASS | — |

The LHAPDF sub-rows were previewed by phase 0 before this phase; each harness
header now says so. No tolerance was touched.

**Physically, per `../../benchmarking/08_channel_priority.md` §4.**
`t1_g1d_world` sits under ingredient **H** (polarized g₁: A_∥ on ranks 3 and
5, nothing on claim 1). Its FAILs are now FAILs **against the like-for-like
data**: the D-state construction is cleared, so the A_∥ price phase A measured
toy → NNPDFpol1.1 (×0.514 tagged ⁷Li α, ×0.532 inclusive ⁷Li) remains a price
between two backends that both fail it, and nothing here decides which side
is right. The NMC rows sit under ingredient **A**: claim 1's cos 2φ amplitude
does not see F₂ (R is the lever, phase C1), A_zz (rank 4) moves ×1.25 from
toy to CT18NLO (window level, `08` row A). The normalisation reading makes "no shipped unpolarised
backend passes absolute F₂ᵈ" robust for toy and CT18NLO (6.0 % and 8.6 %
excesses are 2.4 and 3.4 normalisation σ); for MSTW it is open. The F₂ᵈ data at
x < 0.12 carry NMC's own R — the same R data phase C1 compares the generator's
R with — so the absolute F₂ᵈ comparison and the R comparison share one input.

## 2. D2 — η_d and Q_d

| paper | copy (outside the repository) | sha256 | read |
|---|---|---|---|
| N. L. Rodning, L. D. Knutson, PRC 41 (1990) 898 | APS, curl's own identity, scratch `refs/rodning_knutson1990_PhysRevC.41.898.pdf` (706 610 bytes) | `4dbfb75f5d67e437f9115a8427e6db0508b0f3984e59a2d1eebee79dee0ff161` | pp. 898–900, 906–908 |
| D. M. Bishop, L. M. Cheung, PRA 20 (1979) 381 | APS, the same way, scratch `refs/bishop_cheung1979_PhysRevA.20.381.pdf` (227 994 bytes) | `a28b358d7223447c83f13c6855a6eb551a0a1a5e79350f3c29e75669e4091950` | pp. 381–384 |
| R. B. Wiringa *et al.*, AV18, PRC 51 (1995) 38 | sibling `refs/nucl-th_9408016.pdf` | `44f551a2…10996f5` | Table X |
| R. Machleidt, CD-Bonn, PRC 63 (2001) 024001 | sibling `refs/nucl-th_0006014.pdf` | `208c7a6a…9298de` | Table XV, refs. |

### 2.1 η_d — Rodning & Knutson: a measurement, transcribed right

**Concluded**: η = 0.0256 ± 0.0004, Eq. (10), p. 906 (abstract p. 898; Sec. VII
p. 908), "where the uncertainty includes the contributions from statistics,
calculations, and normalization" (no confidence level stated; read as 1σ, as
before). **What it is**: a MEASUREMENT — an empirical determination through a
reaction model: T₂₀ of ²⁰⁸Pb(d,p)²⁰⁹Pb at 6.0 and 7.0 MeV, T₂₁ of the same at
7.0 MeV, and T₂₀ of ¹³⁶Xe(d,p)¹³⁷Xe at 4.5 and 5.5 MeV, two final states each,
ten angular distributions (Sec. III.A, p. 900), compared with DWBA; ten values,
weighted. **The transcription 0.0256(4) is right; nothing moved.** The file's
2026-09-26 method note ("d + ⁴He sub-Coulomb tensor analysing powers", from
the survey) was wrong about the targets and is corrected in its provenance.

### 2.2 Q_d — Bishop & Cheung: a calculation that confirms someone else's value

**Concluded**: Q = **0.2860 ± 0.0015 fm²** — the abstract (p. 381: "…also
confirms the earlier value of the electric field gradient and necessarily,
therefore, the predicted quadrupole moment of the deuteron (i.e., Q = 0.2860 ±
0.0015 fm²)") and the Conclusions (pp. 382–383: "Our principal conclusion is
that Reid and Vaida's calculation of the electric field gradient … is correct,
as is their value (0.2860 ± 0.0015 fm²) derived therefrom … If anything, their
error limits on Q are on the conservative side."). **What it is**: a
THEORETICAL quantum-chemistry calculation of the electric field gradient q at
the nucleus in D₂ (Kołos–Wolniewicz 54-term Born–Oppenheimer wave function at
R = 1.4 a.u., and a 540-term nonadiabatic one) that confirms Reid & Vaida's q
(PRA 7 (1973) 1841; PRL 34 (1975) 1064E) and hence their **semi-empirical**
Q = (eqQ/h, measured by Ramsey and co-workers in molecular-beam magnetic
resonance on HD and D₂) / q. "There is no purely experimental value" (p. 381).

**The number 0.2859(3) does not occur in the paper.** It is what the NN
potential papers attribute to it: AV18, Table X, "Q_d 0.2859(3)" with footnote
f = Bishop & Cheung; CD-Bonn, Table XV, "0.2859(3) [76, 73]" = Bishop & Cheung
and Ericson & Rosa-Clot, NPA 405 (1983) 497 (not on disk). Where the ±0.0003
comes from is **not established here**. The tree's text
(`include/lipolgen/b1_nuclear.hpp`, the "THE 100 % BAND IS MANDATORY" block)
carried the literature's attribution, and `deuteron_static_measured.json`
transcribed that block.

### 2.3 The correction, and what moved

`validation/benchmarks/data/deuteron_static_measured.json`, revised 2026-10-26
as a **dated correction**: Q_d value 0.2859 → **0.2860**, σ 0.0003 → **0.0015**,
as printed "+0.2860 +- 0.0015 fm^2", reference and `transcribed_from` naming the
pages; η_d's `as_printed` now the paper's "0.0256 +- 0.0004" (value and σ
unchanged); `_provenance.transcription` starts "FIRST-HAND since 2026-10-26";
the 2026-09-26 statement is kept verbatim (`transcription_2026_09_26`), the old
`transcribed_from` likewise (`transcribed_from_2026_09_26`), and
`_provenance.corrections` holds the old value (0.2859, 0.0003, "+0.2859(3)
fm^2") with its source (b1_nuclear.hpp lines 121–128 at `f0a8f1e`) and the
reason. File sha256 `6433c31e…f692d5` → `3e62013f573e76f8e622bff87cca8d954b6616ae3578d99b17368d9e06062e0f`;
body sha256 `18f805e5…6ede849` → `51adce98ee6ebc80fa40278244cffe9d92d94f22ed8c3213b9e62eba75ea4acb`;
both re-recorded in `t2_deuteron_static.py` (`FILE_SHA256`, `BODY_SHA256`,
the old digests in a comment).

| model | η (pull) — **unchanged** | Q_d (IA) | vs 0.2859(3), 2026-09-26 | vs 0.2860(15), now | in σ, now |
|---|---|---|---|---|---|
| CD-Bonn | 0.0255714 (−0.07σ) PASS | 0.270493 fm² | −5.39 % | **−5.42 %** | −10.3 |
| AV18 | 0.0250447 (−1.39σ) PASS | 0.269670 fm² | −5.68 % | **−5.71 %** | −10.9 |
| Hulthén (tagged default) | 0.0324784 (+17.20σ) **FAIL** | 0.290349 fm² | +1.56 % | **+1.52 %** | +2.9 |

**No verdict moves**: Q_d is recorded, not asserted (the status is η's alone,
FAIL, Hulthén +17.2σ outside); no tolerance was touched. In units of the
corrected error the impulse-approximation shortfall is ≈ −10σ for the realistic
waves (it read ≈ −50σ on the old 0.0003) — the size of the two-body and
relativistic terms an IA Q_d omits. The tree's own Q_d input (the RC tail's
deuteron F_q(0), `tests/test_rc.cpp`'s fixture, the band prose) is 0.2859, the
literature's attribution, 0.035 % below the paper's 0.2860 — a constant this
phase does not move.

**Physically** (`08_channel_priority.md` ingredient **G**, ranks 6–7): the
Hulthén FAIL is a FAIL against the right number, measured; η_d constrains the
asymptotic S–D mixture, which CD-Bonn and AV18 have and the tagged channel's
default does not. A_zz^tag's ×2.3–4.8 lever and D2 at M2 are unchanged.

## 3. α_s(M_Z) — PDG 2026, verified

`t5_sum_rules.py` typed α_s(M_Z) = 0.1180 ± 0.0009 from memory on 2026-09-26.
`pip install --user pdg` installed the PDG API **2026.0** (wheel sha256
`681e11f8…b45c8e8a`, the same file the harness cites; `pdg.sqlite` data release
2026-05-30; it also installed SQLAlchemy 2.0.55 and greenlet 3.5.6 into
`~/.local`): its database has **no α_s(M_Z) row** (searched by description and
by value), while S017AV (−1.2753 ± 0.0013), S044M (91.1879 ± 0.0020 GeV), Q005M
(4.186 ± 0.006 GeV) and Q004M (1.2729 ± 0.0045 GeV) are the harness's. PDG's
web site answers the standard Python client from this machine (no identity
change), so the source was read directly:

* Table 1.1 "Physical Constants" (revised May 2026 by D. Robinson and
  P. A. Zyla), `https://pdg.lbl.gov/2026/reviews/rpp2026-rev-phys-constants.pdf`
  (252 551 bytes, sha256
  `f54f2e23ca28f945158a01649e199653ef443773337c3fb99e1482baa9ac23c1`):
  "strong coupling constant α_s(m_Z) 0.1180(9)".
* Sec. 9 "Quantum Chromodynamics" (revised August 2025, pages dated 1 June
  2026), `https://pdg.lbl.gov/2026/reviews/rpp2026-rev-qcd.pdf` (sha256
  `c04c628d76b18610c5fa2a919c6081918a25b55fb971b6af5829f4ca2baa386f`): the world
  average α_s(m_Z²) = 0.1180 ± 0.0009, Eq. (9.25), p. 42 (0.1178 ± 0.0010
  without lattice, Eq. (9.24)).

**The value typed is the value printed**: no constant moved,
`REFERENCE_SHA256` is unchanged, Γ_ref = 0.18538 ± 0.00233 at Q² = 5 GeV²
unchanged, no number of the row moved. The header's three lines that called
α_s "NOT re-read … the least verified input" now say "VERIFIED 2026-10-26 in
PDG 2026 itself: Table 1.1 "0.1180(9)" and the QCD review, Eq. (9.25), p. 42"
and point here for the sources and digests — written line for line in the
same three lines (and the "403" remark re-dated to the 2026-09-26 machine in
place), so that the line citations `phase_A_baseline.md` §3 makes into this
file stay valid (the records gate reads 0 broken); a first, longer version
shifted them and was withdrawn. `test_bench_t5_closures.py`'s reference test
gains a comment and an assertion on the error (0.0009), no re-pin.

## 4. Suelzle 1967 and Moniz 1971 — recorded only, nothing implemented

`../../benchmarking/03_data_nuclear.md` §7 rank 6 ("Fit a ⁷Li `HoSpin1FF` on
Suelzle's published HO parameters (a = 1.77(2), α = 0.327)") and the
`RC_QE_KF_GEV` row. Copies: Phys. Rev. 162 (1967) 992 (scratch
`refs/suelzle1967_PhysRev.162.992.pdf`, sha256 `ce366617…cc6e22d7`, pages read
from the page images where the scan's text is noisy) and PRL 26 (1971) 445
(scratch `refs/moniz1971_PhysRevLett.26.445.pdf`, sha256 `2bdde956…5bc199e4`).

### 4.1 ⁷Li (Suelzle, Yearian, Crannell)

* **Simple harmonic-well shell model + undeformed p-shell quadrupole** (the
  paper's model of record for ⁷Li): F_C0 = [1 − (Z − 2)/(6Z) q²a₀²]
  exp[−¼q²a₀²(1 − 1/A)] F_N(q), Eq. (41); best fit **a₀ = 1.755 F, |Q| = 4.20 F²,
  χ² = 23.1 for 27 d.o.f.** (p. 1003); χ² ≤ χ²_min + 1: |Q| = 4.35 F² at
  a₀ = 1.767 F and 4.00 F² at a₀ = 1.745 F; all F_ch² −5 %: a₀ = 1.792 F,
  |Q| = 4.70 F² (χ² 25.0); +5 %: a₀ = 1.728 F, |Q| = 3.90 F² (χ² 22.8) (p. 1003).
  Conclusion (p. 1004): **a₀ = 1.755 ± 0.02 F, |Q| = 4.20 ± 0.25 F², rms charge
  radius 2.39 ± 0.03 F** (Eq. (42) with r_p² = 0.726 F², Z = 3, A = 7; abstract
  p. 992). The deformed-oscillator quadrupole model does not fit.
* **α is not a fitted parameter** in this model: the monopole's q² coefficient
  is fixed by the p shell, (Z − 2)/(6Z).
* The **modified harmonic well** (C0 only) also fits: a₁ = 1.55 ± 0.015 F,
  a₂ = 2.02 ± 0.06 F, rms 2.43 ± 0.03 F (χ² 19.8/27, p. 1002) — "a physically
  unlikely model" (p. 1004).

**The survey's (a, α) = (1.77(2), 0.327) is not printed by Suelzle.** It is de
Vries's re-expression (ADNDT 36 (1987) 495, Table I, as `03_data_nuclear.md` §2
reads it) of this same fit in the de Vries HO form
ρ ∝ [1 + α(r/a)²] e^{−(r/a)²}, F = [1 − αq²a²/(2(2 + 3α))] e^{−q²a²/4}, with the
c.m. factor and the proton folded in. Arithmetic check (a Gaussian proton of
⟨r_p²⟩ = 0.726 F² assumed): a² = a₀²(1 − 1/A) + (2/3)⟨r_p²⟩ = 3.0800 × 6/7 +
0.4840 = 3.124 F², **a = 1.7675 F**; matching the q² coefficient,
α/(2(2 + 3α)) = (1/18)a₀²/a² = 0.05478 ⇒ **α = 0.3264**; and ⟨r²⟩ =
(3/2)a²(2 + 5α)/(2 + 3α) = 5.713 F², **2.390 F** — Suelzle's 2.39 F. So the
two pairs are one fit in two conventions. **For the rank-6 action this
matters**: `HoSpin1FF` takes a **point** shape and folds the nucleon form
factor separately (`include/lipolgen/rc.hpp`, the `fold_nucleon` option, and
the `LI6_FF_HO_A_FM` block's ⟨r²⟩_point convention), so a ⁷Li fit in the
tree's convention must start from Suelzle's a₀ and Eq. (41), not from de
Vries's folded (a, α). Nothing is implemented.

### 4.2 ⁶Li (the same paper)

* **No harmonic-well fit**: the simple HO (a₁ = a₂) fits nowhere, not even for
  q² ≤ 3 F⁻²; the modified HO fits only q² ≤ 3 F⁻²: a₁ = 1.632 F, a₂ = 1.98 F,
  rms 2.426 F (χ² = 19 for 17 d.o.f.); with F_ch² ×1.05: a₁ = 1.625, a₂ = 1.94 F,
  2.398 F; ×0.95: a₁ = 1.645, a₂ = 2.02 F, 2.458 F (Eq. (27), r_p² = 0.726 F²;
  p. 1000).
* **Phenomenological fit** over the full range, Eq. (28):
  F(q²) = exp(−a²q²) − c²q² exp(−b²q²), ⟨r²⟩ = 6(a² + c²): **a² = 0.87 F²,
  b² = 1.70 F², c² = 0.205 F², rms 2.540 F, χ² = 26 for 32 d.o.f.** (p. 1001);
  within χ²_min + 1 the rms runs 2.510–2.575 F; F_ch² ×1.05: (0.865, 1.67, 0.160)
  → 2.480 F; ×0.95: (0.875, 1.73, 0.250) → 2.598 F; the paper suggests **2.54 ±
  0.05 F** (pp. 1001, 1004; abstract 2.54 ± 0.05 F). A low-q² expansion
  (q² ≤ 1.2 F⁻²) gives 2.5 ± 0.1 F (p. 1001).
* **The ratio** ⟨r²⟩^½(⁶Li)/⟨r²⟩^½(⁷Li) = 1.055 ± 0.008 from the 100 MeV data
  (0.1 < q² < 0.6 F⁻², p. 1004).
* The paper quotes Bernheim's simple-HO a₀ = 1.78 F (rms 2.41 F) for ⁶Li over
  q² ≤ 1.5 F⁻² (p. 1004) — not its own fit.

Against the tree: `LI6_FF_HO_A_FM = 1.9069`, `LI6_FF_HO_ALPHA = 0.13822`
(`rc.hpp`) are a **phenomenological point-shape** pair built on Angeli &
Marinova's r_ch and q₀ = 3.1 fm⁻¹ — Suelzle offers no ⁶Li HO pair to compare
them with; its rms 2.54 ± 0.05 F is the Su67 entry, one of the three de Vries
Table I analyses the `t3_li6_rms_devries` row already uses, now confirmed in
Suelzle's own text (abstract p. 992; p. 1001).

### 4.3 Moniz 1971 — the ⁶Li k_F row

Table I (p. 447): **₃Li⁶: k_F = 169 MeV/c, ε̄ = 17 MeV**, "The fitting
uncertainty in these numbers is approximately ±5 MeV/c" (k_F) and "±3 MeV"
(ε̄); 500 MeV, 60° (p. 445), a Fermi-gas fit with no overall normalisation
factor (p. 447). `RC_QE_KF_GEV = 0.169` (`rc.hpp`) **is that row**; the ±0.005 GeV
fitting uncertainty is not carried in the tree. Recorded only.

## 5. What this phase does NOT validate

The physics of any row: no verdict, tolerance, default or constant moved, and
no comparison was added. The COMPASS high-x printed values were not put into
any χ²; NMC's 2.5 % normalisation was not added to any error; NNPDFpol1.1's
leading-order evaluation was not tested; where 0.2859(3) came from (Ericson &
Rosa-Clot, not on disk) is not established; de Vries ADNDT 36 was not
re-read (its (1.77(2), 0.327) is reproduced by arithmetic from Suelzle, §4.1,
under an assumed Gaussian proton); NNPDFpol1.0's deuteron convention was not
read (not on disk).

## 6. Files touched

* `validation/benchmarks/t1_g1d_world.py` — the header (WHAT IS COMPARED per
  set, read; THE D-STATE CONVENTION … RESOLVED, the rule, the four papers'
  statements and pages, the phase-0 preview; rule 1 per nucleon from the
  papers; rule 3 SMC's target corrected; E143's 'sys_beam' named; SMC's
  recorded-only reason) and two `EXPERIMENTS` reference labels (E143: "E143
  PRD 58 (1998) 112003 (hep-ph/9802357) Table XIV, HEPData ins467140 t18";
  SMC: "SMC PRD 58 (1998) 112001 Table IX, HEPData ins471981 t7 (recorded
  only)"). No number, tolerance, cut or code path; `plot()` untouched.
* `validation/benchmarks/t1_nmc_f2d.py` — header only (Table 6 read; the
  normalisation resolved; the R NMC used; the nuclear-effects reading; the
  CT18NLO preview).
* `validation/benchmarks/t1_nmc_f2d_over_f2p.py` — header only (Table 3 read;
  the legacy components identified; which R; the normalisation; nuclear
  effects; the CT18NLO preview).
* `validation/benchmarks/t2_deuteron_static.py` — header (the papers read; the
  Q_d correction and its effect), `FILE_SHA256`/`BODY_SHA256`, the reference
  string of the REPORT row.
* `validation/benchmarks/data/deuteron_static_measured.json` — the dated Q_d
  correction (§2.3).
* `validation/benchmarks/t5_sum_rules.py` — the α_s provenance lines of the
  header only, line for line (§3); `REFERENCE` and its digest unchanged.
* `python/tests/test_bench_t1_nucleon.py` — the module docstring; the
  D-state-label test's docstring; five new tests:
  `test_g1d_papers_read_relabel_the_harness_not_the_files`,
  `test_g1d_pinned_ct18nlo`, `test_g1d_pinned_nnpdfpol_and_the_d_state_reading`,
  `test_f2d_pinned_ct18nlo`, `test_dp_pinned_ct18nlo` (the LHAPDF pins skip with
  the harness's blocked reason where a set is missing). No existing pin moved.
* `python/tests/test_bench_static.py` — the deuteron tests only:
  `test_deuteron_static_vendored_file_is_the_transcription` (the new values,
  provenance and correction, the flip target `FIRST-HAND`), the Q_d relative
  distances re-pinned to −0.05422 / −0.05710 / +0.01521 (old values in a
  comment), the verdict test's 100×-off fake Q_d derived from the reference,
  the exit-code flip on `"value": 0.286,`.
* `python/tests/test_bench_t5_closures.py` — `test_sr_reference_constants_are_the_recorded_ones`:
  a comment and `alpha_s_mz_err == 0.0009`.
* `docs/benchmarking/figures/t2_deuteron_static.png` — regenerated once by
  `t2_deuteron_static.py --plot` (the drawn Q_d reference is now
  0.2860 ± 0.0015).
* This record.

## 7. Commands

```sh
cd /home/cpeng/Projects/polli/PoliGen && source env.sh
python3 validation/benchmarks/t1_g1d_world.py            # exit 0, FAIL row
python3 validation/benchmarks/t1_nmc_f2d.py              # exit 0, FAIL row
python3 validation/benchmarks/t1_nmc_f2d_over_f2p.py     # exit 0, FAIL row
python3 validation/benchmarks/t2_deuteron_static.py --plot   # exit 0, FAIL row; the figure
python3 validation/benchmarks/t5_sum_rules.py            # exit 0, FAIL row
python3 -m pytest python/tests/test_bench_t1_nucleon.py -q          # 49 passed, 1 skipped
python3 -m pytest python/tests/test_bench_static.py -q              # 18 passed
python3 -m pytest python/tests/test_bench_t5_closures.py -q -k "sr_ or sum_rules"   # 14 passed
python3 validation/benchmarks/make_report.py \
    --only t1_g1d_world,t1_nmc_f2d,t1_nmc_f2d_over_f2p,t2_deuteron_static,t5_sum_rules \
    --out-dir <scratch>/report_D --compare docs/benchmarking/REPORT.json
pdftotext -layout <pdf> <txt>; pdftoppm -r 150 -f N -l N -png <pdf> <prefix>   # every page cited
```

## 8. For the integrator (shared files this phase may not edit)

The exact text is in the phase's returned `notes_for_integrator`; in short:
`validation/benchmarks/README.md`'s T2 §5.2 row (the reference cell names
0.2859(3) and "the papers not re-read"; the Q_d percentages are now −5.42 /
−5.71 / +1.52 %); `../../benchmarking/08_channel_priority.md` §4 rows H, A and
G (the D-state cleared; the normalisation reading; η_d confirmed);
`../../benchmarking/02_data_nucleon_deuteron.md` §5.2 (the "d + ⁴He" method is
wrong: ²⁰⁸Pb and ¹³⁶Xe); `phase_P_figures.md`'s `t2_deuteron_static` row
(written against 0.2859(3); the figure was regenerated here); the comments
that attribute 0.2859(3) to Bishop–Cheung (`include/lipolgen/b1_nuclear.hpp`)
and 0.0256(4) to CD-Bonn's Table XV (`tests/test_cluster.cpp`) — comments only,
no constant; the vendoring script's provenance strings for the T1 files
(E143 "CONFLICT", SMC "NOT IDENTIFIED", D-state "UNVERIFIED"), which only a
regeneration of the pinned files could change and which this phase therefore
leaves as they are.
