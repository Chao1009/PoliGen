<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# 03 — External data: the UNPOLARIZED nuclear inputs

**Survey date 2026-09-06.** Domain: measured data on ⁶Li, ⁷Li, ⁴He and ³He that
constrains LiPolGen's *nuclear* inputs. Companion to
`00_in_tree_checks.md`, which inventories what the tree already does; this file
lists only **external** resources and says, for each, exactly which LiPolGen
number it pins and which it does not.

Everything below is **unpolarized**, because §0 of `00_in_tree_checks.md` holds:
there is no polarized e + ⁶Li or e + ⁷Li DIS data anywhere, no other generator
does tensor-polarized lithium, and the world's only tensor-DIS datum is HERMES's
deuteron b₁. Nothing in this file changes any of that. What this file *does*
change is a claim that has been made loosely in this tree several times — that
"there is no lithium data". There is a great deal of lithium data. It is
unpolarized, it is old, most of it is not machine-readable, and **one dataset —
the single most important one — is on HEPData under a CC0 licence and was not
known to this repository.**

## 0. The headline

| | |
|---|---|
| **Found, and it is a real benchmark** | **NMC measured F₂(⁶Li)/F₂(D) directly.** Arneodo *et al.* (NMC), Nucl. Phys. **B441** (1995) 12, `hep-ex/9504002`, HEPData `ins394050` Table 1: **24 points**, x = 1.4×10⁻⁴ … 0.65, Q² = 0.034 … 39 GeV², stat+syst per point, CC0. The target was **⁶Li** — three 13 cm cylinders, 17.4 g/cm² total, 4.5 % ⁷Li contamination (paper §2, quoted below). This is a direct, isoscalar-on-isoscalar measurement of the EMC/shadowing ratio for the *exact nucleus the generator is built around*, and the tree's `EMC_VALENCE_DEPLETION_EPPS21` presently rests on an **nPDF fit** instead. |
| **Found, and it directly measures the tagged channel's own wave function** | **⁶Li(e,e′d)⁴He and ⁶Li(e,e′α)²H were both measured.** Ent *et al.*, Nucl. Phys. **A578** (1994) 93 (d knockout); Mitchell *et al.*, Phys. Rev. C **44** (1991) 2002 (α knockout). Plus ⁶Li(p,pd)⁴He at 670 MeV and ⁶Li(α,2α)²H at 700 MeV. These measure the **α–d relative momentum distribution** that `Wave::radial` / `data/vmc/li6_alpha_d/li6.ad` model and that the T1 tagged channel samples. |
| **Found, and it is a data-vs-VMC test of the exact wave-function family in `data/vmc/`** | Lapikás, Wesseling, **Wiringa**, Phys. Rev. Lett. **82** (1999) 4404, `nucl-th/9904008`: ⁷Li(e,e′p)⁶He momentum distributions vs **VMC**; summed spectroscopic factor **0.58 ± 0.05** measured against **0.60** VMC. |
| **Confirmed absent** | No Li in SLAC E139 (targets ᴰ, He, Be, C, Al, Ca, Fe, Ag, Au — HEPData `ins359103` qualifiers). No Li in JLab E03-103 (²H, ³He, ⁴He, ⁹Be, ¹²C). No Li in JLab E12-10-008 (⁹Be, ¹⁰B, ¹¹B, ¹²C). **No ⁷Li DIS data at any x.** No sum-of-Gaussians charge-density parameterisation for ⁶Li or ⁷Li in the de Vries compilation. |

---

## 1. The EMC/shadowing input — `EmcBaseline::Epps21`

### 1.1 What the tree does today

`EMC_VALENCE_DEPLETION_EPPS21 = 0.031052077003862335` (`include/lipolgen/sf.hpp`)
is ⟨1 − R_unpol⟩ over x ∈ [0.35, 0.65] at Q² = 5, with
R = F₂(`EPPS21nlo_CT18Anlo_Li6`)/F₂(CT18ANLO). It is the **scale** by which the
CBT (⁷Li) and TMT (nuclear matter) polarized-EMC curves are transferred onto a
common baseline, so every polarized-EMC number in the generator is proportional
to it. It comes from a **global nPDF fit**, and EPPS21's ⁶Li grid is constrained
by the very NMC data below plus A-dependence systematics — it is not a
measurement of ⁶Li, it is an interpolation *informed by* one.

### 1.2 The measurement, and a zero-free-parameter comparison run for this survey

HEPData `ins394050` Table 1 (`10.17182/hepdata.47955.v1/t1`), all 24 points, is
reproduced below beside EPPS21 evaluated **in this checkout** through
`lipolgen.Epps21Ratio` (`EPPS21nlo_CT18Anlo_Li6` over the shipped `CT18NLO`
proton, since CT18ANLO is not installed here), with the free-nucleon denominator
built **isoscalar** as ½(F₂ᵖ + F₂ⁿ) from `lipolgen.LhapdfSF`:

| x | Q² [GeV²] | NMC F₂(⁶Li)/F₂(D) | EPPS21 R_iso at NMC's Q² | at Q² = 5 |
|---|---|---|---|---|
| 0.00014 | 0.034 | 0.935 ± 0.069 | *below grid* | 0.9034 |
| 0.00028 | 0.066 | 0.922 ± 0.059 | *below grid* | 0.9045 |
| 0.00045 | 0.11 | 0.865 ± 0.044 | *below grid* | 0.9054 |
| 0.00067 | 0.15 | 0.897 ± 0.039 | *below grid* | 0.9065 |
| 0.00090 | 0.21 | 0.894 ± 0.034 | *below grid* | 0.9080 |
| 0.0015 | 0.34 | 0.888 ± 0.025 | *below grid* | 0.9116 |
| 0.0035 | 0.61 | 0.921 ± 0.013 | *below grid* | 0.9246 |
| 0.0055 | 1.0 | 0.947 ± 0.011 | *below grid* | 0.9357 |
| 0.0085 | 1.4 | 0.951 ± 0.011 | *below grid* (1.4 < 1.69 GeV²; the LHAPDF extrapolation reads 0.9469) | 0.9488 |
| 0.0125 | 1.8 | 0.980 ± 0.010 | 0.9545 | 0.9616 |
| 0.0175 | 2.3 | 0.992 ± 0.010 | 0.9679 | 0.9729 |
| 0.025 | 2.8 | 0.994 ± 0.009 | 0.9804 | 0.9841 |
| 0.035 | 3.6 | 0.971 ± 0.009 | 0.9921 | 0.9938 |
| 0.045 | 4.2 | 0.997 ± 0.011 | 0.9997 | 1.0003 |
| 0.055 | 4.8 | 0.986 ± 0.012 | 1.0041 | 1.0042 |
| 0.070 | 6.0 | 1.003 ± 0.011 | 1.0078 | 1.0076 |
| 0.090 | 7.2 | 1.004 ± 0.013 | 1.0097 | 1.0100 |
| 0.125 | 9.2 | 1.009 ± 0.011 | 1.0095 | 1.0102 |
| 0.175 | 12.0 | 1.002 ± 0.015 | 1.0054 | 1.0071 |
| 0.25 | 17.0 | 0.992 ± 0.016 | 0.9969 | 0.9992 |
| 0.35 | 23.0 | 0.990 ± 0.027 | 0.9844 | 0.9868 |
| 0.45 | 29.0 | 0.910 ± 0.043 | 0.9728 | 0.9745 |
| 0.55 | 33.0 | 0.988 ± 0.082 | 0.9651 | 0.9642 |
| 0.65 | 39.0 | 0.801 ± 0.109 | 0.9698 | 0.9596 |

(errors stat ⊕ syst; NMC additionally quotes a 0.004 overall normalisation
uncertainty not included. "Below grid" = Q² under the EPPS21 grid minimum.)

Numbers this survey computed from that comparison, all **new** and none of them
in the tree:

* Over the four points with x ≥ 0.30, **χ² = 4.63 for 4 points with zero free
  parameters.** EPPS21's ⁶Li valence depletion is consistent with the only
  measurement of it. Over all 15 on-grid points χ² = 26.94/15 (p = 0.029), carried by the shadowing/anti-shadowing crossover (pulls +2.55, +2.34, +1.59, −2.23 at x = 0.0125 … 0.035); wired as `validation/benchmarks/t3_nmc_li6_over_d.py` 2026-09-23 (pass at p ≥ 0.01 on both sets; the x ≥ 0.30 χ² reproduced as 4.6269/4).
* Weighted mean of those four NMC points: **0.9618 ± 0.0219**, i.e. a measured
  valence depletion **0.038 ± 0.022** against the library's **0.031052**. The
  agreement is 0.3 σ (against the library's window mean at Q² = 5 on CT18ANLO; like for like — the tree at the same four points, Q² and weights, 0.9795 — the distance is −0.81 σ) — and the **±0.022 is the honest size of the external
  constraint**: the data pin ⟨1 − R⟩ to about ±70 % of its own value, so this
  benchmark bounds a gross error and cannot refine the constant.
* Recomputing the library's own definition here with the *shipped* CT18NLO
  proton instead of CT18ANLO gives ⟨1 − R⟩ = **0.029803** against the recorded
  **0.031052** — a **4.0 % relative** move from the baseline swap alone, which
  is the concrete size of the "few tenths of a percent" the header claims for
  the two proton baselines (a few tenths of a percent *on R* is a few percent
  *on 1 − R*, and it is 1 − R that is used).
* **Shadowing cannot be compared point by point.** NMC's saturation value is
  **0.890 ± 0.010(stat) ± 0.021(syst) for x < 0.002** (paper §6, verbatim), but
  those points sit at Q² = 0.03 … 0.34 GeV², far below the EPPS21 grid. EPPS21
  at Q² = 5 gives 0.903–0.912 there. The numbers are close; the comparison is
  **not** apples to apples and must not be pinned as one.

### 1.3 A trap anyone wiring this gate will hit

`Epps21Ratio::f2_per_nucleon` returns F₂^A **per nucleon**, correctly. But
`EPPS21nlo_CT18Anlo_Li6` in LHAPDF is the **isoscalar-averaged nucleon** of ⁶Li,
not the bound proton: measured here, R_u = 0.578 and R_d = 2.829 at x = 0.65,
Q² = 5, which are exactly (1 + d/u)/2 and (1 + u/d)/2 — the n/p isospin average,
not a nuclear effect. Dividing `f2_per_nucleon` by the **free proton** F₂
therefore gives 0.6880 at x = 0.65, Q² = 5 (0.6924 at NMC's Q² = 39) and an apparent ⟨1 − R⟩ of **0.268** over the
valence window — nine times the real number, and entirely an artefact. The
denominator must be the isoscalar ½(F₂ᵖ + F₂ⁿ). This is written down because the
mistake is silent, reproduces a plausible-looking EMC curve, and would corrupt
every polarized-EMC transfer downstream.

### 1.4 What the NMC gate does NOT validate

* **The deuteron denominator.** NMC measures Li/D, not Li/(free isoscalar
  nucleon). The deuteron carries its own few-percent EMC/binding effect at
  x ≳ 0.5, so the true ⁶Li-to-free-nucleon depletion is **larger** than the NMC
  ratio by that amount, in the direction that would *increase* 0.031. Any gate
  must either state the omission or apply a deuteron correction.
* **The polarized EMC effect.** Nothing here touches it. The CBT-vs-TMT "×2 vs
  ×1" disagreement (`00_in_tree_checks.md` E-11) is untouched by a factor of 2
  in either direction by this data.
* **⁷Li.** There is no ⁷Li DIS measurement. CBT's calculation is *for* ⁷Li and
  has no ⁷Li data under it either.
* **Q² evolution at large x.** The x ≥ 0.30 points live at Q² = 23–39 GeV², not
  5; measured here, that costs ≤ 0.006 in R, so it is not the limiting
  systematic — the ±0.022 statistics are.

### 1.5 The other NMC record, and the A-dependence programme that has no Li

* **HEPData `ins393377`** (Amaudruz *et al.*, NPB **441** (1995) 3,
  `hep-ph/9503291`), Tables 4 and 5: **F₂(C)/F₂(⁶Li)** and **F₂(Ca)/F₂(⁶Li)**,
  25 points each, 90 GeV, 0.0085 < x < 0.6, 0.84 < Q² < 17 GeV². Same ⁶Li
  target. These are *A-dependence* ratios with ⁶Li as the light reference — a
  weaker but independent handle, and their Q² coverage in the valence region is
  lower (better matched to Q² = 5) than Table 1's.
* **SLAC E139** (Gomez *et al.*, PRD **49** (1994) 4348, HEPData `ins359103`):
  reaction qualifiers are `E- P/HE/BE/C/AL/CA/FE/AG/AU --> E- X`. **No lithium.**
  Listed here so the negative is on the record with its source.
* **JLab E03-103** (Seely *et al.*, PRL **103** (2009) 202301, `0904.4448`):
  ²H, ³He, ⁴He, ⁹Be, ¹²C at 0.3 < x < 0.9, Q² ≈ 3–6 GeV². **No Li**, and no
  HEPData record — the EMC slopes are in **Fig. 4 only** (checked: the PRL has
  no numerical table). ⁹Be is the cluster analogue: two α + one neutron, an EMC
  effect that tracks *local* density rather than average density, which is the
  same argument that would apply to ⁶Li's α + d. It is the nearest measured
  statement that cluster structure matters for the EMC effect, and it is
  **qualitative** for our purposes.
* **JLab E12-10-008** (Arrington *et al.*, PRC **104** (2021) 065203,
  `2110.08399`): ²H, ³He, ⁴He, ⁹Be, ¹²C, ⁶³Cu, ¹⁹⁷Au, 0.3 ≤ x ≤ 1, Q² up to
  8.3 GeV²; ¹⁰B and ¹¹B measured for the first time. **No Li.** No HEPData
  record; tables are in the paper.

**So: the light-nucleus EMC programme has systematically skipped lithium since
E139, and the one lithium measurement in existence is a 1995 muon experiment.**

---

## 2. Elastic form factors — the `HoSpin1FF` open item (Q1)

`include/lipolgen/rc.hpp` says of `LI6_FF_HO_A_FM = 1.9069` /
`LI6_FF_HO_ALPHA = 0.13822`: *"TO BE REFIT by the implementer against the
Suelzle-Yearian-Crannell … and Li-Sick-Whitney-Yearian … elastic data"*, and
that those data are *"STILL not in this repository in any machine-readable
form"*. Both papers are real and this survey verified them; below is what is
actually obtainable.

| resource | what it is | obtainable |
|---|---|---|
| **L. R. Suelzle, M. R. Yearian, H. Crannell, Phys. Rev. 162 (1967) 992–1005** (`10.1103/PhysRev.162.992`) | "Elastic Electron Scattering from Li⁶ and Li⁷", Stanford HEPL, beam **100–600 MeV**. The origin of the ⁶Li **and** ⁷Li charge form factors. | Journal paywalled; APS returns 403 to non-browser clients from here. **I could not open the article body, so I cannot confirm whether it prints cross-section tables.** Contemporaneous *Phys. Rev.* practice was to tabulate; treat as "probably table-in-paper, unverified". |
| **G. C. Li, I. Sick, R. R. Whitney, M. R. Yearian, Nucl. Phys. A162 (1971) 583** | "High-energy electron scattering from ⁶Li", **500 MeV**, q² up to 13 fm⁻², diffraction minimum at q² = 8 fm⁻² (**q ≈ 2.83 fm⁻¹**). | Elsevier paywalled. Same caveat. **A DTIC report version exists** (AD0714022, "High Energy Electron Scattering from ⁶Li") which may be openly downloadable; DTIC returned 403 to this environment, so **unverified**. |
| **H. de Vries, C. W. de Jager, C. de Vries, ADNDT 36 (1987) 495** | The standard compilation. **Fetched and read for this survey** (a scan is mirrored at `discovery.phys.virginia.edu/research/groups/ncd/dldata/1987_source.pdf`). | **Table-in-paper, and I read the rows.** See below — this is the concrete answer to Q1. |
| **G. J. C. van Niftrik, L. Lapikás, H. de Vries, G. Box, Nucl. Phys. A174 (1971) 173–192** | "Magnetization distribution of the ⁷Li nucleus … through 180°. The electric quadrupole moment of ⁷Li." 25–90 MeV, < 200 keV resolution, ground-state doublet. | Paywalled; title/authors/volume verified, body not read. |
| **R. B. Wiringa, R. Schiavilla, PRL 81 (1998) 4317, `nucl-th/9807037`** ("WS98") | VMC AV18+UIX ⁶Li longitudinal and transverse form factors and transition form factors. Verified. | The tree already records that it is **figures only**. This survey confirms nothing in the abstract or metadata contradicts that. It is **theory**, not data, and its Q(⁶Li) = −0.23(9) fm² is 3× the measured value — which is why `rc.hpp` takes only the *monopole shape* from it. |

### What de Vries actually contains for lithium (read, not inferred)

**Table I** (charge-distribution parameters), ⁶Li rows:

| model | ⟨r²⟩^½ [fm] | q-range [fm⁻¹] | source key |
|---|---|---|---|
| MI30 | 2.54(5) | 0.69 – 2.52 | Su67 (Suelzle 1967) |
| MI30 | **2.56(5)** | **0.56 – 3.66** | Li71a (Li–Sick–Whitney–Yearian) |
| MI | 2.57(10) | 0.09 – 0.90 | Bu72 |

**⁷Li row:** model **HO**, ⟨r²⟩^½ = **2.39(3) fm**, **a = 1.77(2) fm,
α = 0.327**, q-range 0.69 – 2.62 fm⁻¹, from Su67.

**Table II** (isotope differences), ⁷Li − ⁶Li Δ⟨r²⟩^½: −0.08(2) (Be65),
−0.13(2) (Su67), −0.003(20) (Ni71) fm — three analyses that do not agree.

**Table V (sum-of-Gaussians) covers ³H, ³He, ⁴He, ¹²C, ¹⁶O, ⁴⁸Ca, ⁵⁸Ni, ¹¹⁶Sn,
¹²⁴Sn, ²⁰⁸Tl/²⁰⁸Pb — and NO lithium isotope.** Verified by reading both Table V
pages.

### What that means for Q1, concretely

1. **There is no published closed-form ⁶Li charge form factor to drop in.** The
   compilation gives a *radius*, not a shape. `TabulatedSpin1FF` therefore has
   nothing to load unless Suelzle's or Li's printed cross sections are
   transcribed from the journal, or WS98 Fig. 1 is digitized.
2. **The de Vries ⁶Li rms is a real, immediate gate.** ⟨r²⟩^½_ch = 2.54–2.57 fm
   from three analyses against Angeli & Marinova's 2.589(39) fm, which is what
   `LI6_R2_POINT_FM2 = 6.0788 fm²` is built on. A gate asserting that the
   `HoSpin1FF` charge radius lies inside the *spread of the electron-scattering
   analyses* rather than on one compilation's central value is a one-line test
   and is more honest than the present single anchor.
3. **⁷Li's HO parameters exist and ⁶Li's do not.** `HoSpin1FF::for_ion` throws
   for ⁷Li today. Suelzle's (a, α) = (1.77(2), 0.327) is exactly the shipped
   functional form, so **⁷Li is the isotope for which the shape is fitted and
   published** — a fitted ⁷Li `HoSpin1FF` is a *cheaper and better-founded*
   object than the ⁶Li one now shipped, and it would give the ⁶Li starting
   values an external sanity check (⁶Li's α = 0.138 against ⁷Li's fitted
   0.327 is a factor 2.4; the shell-model value would be 1/3 for ⁶Li).
4. **The first C₀ zero can be bounded from data after all.** Li *et al.* report
   the ⁶Li diffraction minimum at q² = 8 fm⁻², i.e. **q ≈ 2.83 fm⁻¹**. The
   shipped model puts its zero at **3.0998 fm⁻¹** and T11 asserts only
   [2.9, 3.3] fm⁻¹ — **which excludes the measured minimum**. That is a
   substantive finding: the assertion window is on the wrong side of the one
   number the elastic data supply. It is *not* a contradiction yet — the
   measured minimum is in |F_L|² = F_C0² + F_C2² folded with the nucleon form
   factor, while q₀ is the zero of the *point* C0 — but the gap must be
   resolved before the [2.9, 3.3] window is defended, and the resolution is
   arithmetic the tree can do today. Done 2026-09-23: in `HoSpin1FF` the C2 shares the C0 monopole, so the model's own \|F_L\|² minimum IS q₀ (to 1.7e−11 fm⁻¹) — the C2 fill-in cannot reconcile it with 2.83; only Coulomb distortion, not computed, remains between the two (`../open_items/run_2026-09-23/phase_B2_nuclear.md` §3.2).

### What none of this validates

The **quadrupole shape**. No experiment has separated F_C2 for ⁶Li: elastic
longitudinal data measure F_C0² + F_C2², and C2 is sub-dominant exactly where
those data are good. `rc.hpp` already says this. Nothing found in this survey
changes it — the C2 has a measured q → 0 limit (the quadrupole moment) and a
model in between, and `fq_scale` remains the honest way to carry it.

---

## 3. The α–d momentum distribution — the tagged channel's own observable

This is the strongest finding after §1 and it appears to be new to the tree:
**the α–d relative momentum distribution of ⁶Li has been measured directly, in
four independent reactions.** The T1 tagged channel samples exactly this
distribution (`ClusterChannel::momentum_density`, `Wave::radial`,
`data/vmc/li6_alpha_d/li6.ad`), and `00_in_tree_checks.md` E-13 validates it
only against **ANL VMC**, i.e. against theory.

| # | measurement | reference | what it gives |
|---|---|---|---|
| 1 | **⁶Li(e,e′d)⁴He** | R. Ent, B. L. Berman, H. P. Blok, J. F. J. van den Brand, W. J. Briscoe, M. N. Harakeh, E. Jans, P. D. Kunz, L. Lapikás, **Nucl. Phys. A578 (1994) 93–133** ("The (e,e′d) reaction on ⁴He, ⁶Li and ¹²C") | The deuteron-cluster momentum distribution in ⁶Li, from an exclusive electron-scattering measurement, analysed both in a quasi-elastic cluster-knockout picture and microscopically. The abstract states the quasi-elastic description is valid *"only in limiting cases like the ⁶Li(e,e′d) reaction"* — i.e. ⁶Li is the case where the cluster picture works. Earlier letter: R. Ent *et al.*, **PRL 57 (1986) 2367**. |
| 2 | **⁶Li(e,e′α)²H** | J. H. Mitchell, H. P. Blok, B. L. Berman, W. J. Briscoe, M. A. Daman, R. Ent, E. Jans, L. Lapikás, J. J. M. Steijger, **Phys. Rev. C 44 (1991) 2002** ("Mechanism of the ⁶Li(e,e′α) reaction") | The **α**-knockout mirror: parallel kinematics, Q² dependence of two- and three-body breakup at fixed recoil momentum, and the **recoil-momentum dependence of the two-body α–d breakup**; the data indicate the mechanism is quasi-elastic. This is the *same* channel LiPolGen tags, with the roles of detected and undetected fragment exchanged — by momentum conservation it is the same relative-momentum distribution. |
| 3 | **⁶Li(e,e′α) and ⁶Li(e,e′d) at 520 MeV** | Phys. Lett. B **51** (1974), Elsevier PII `037026937490714X` | The first-generation version of 1 and 2. Superseded; listed for completeness. |
| 4 | **⁶Li(p,pd)⁴He at 670 MeV** | Nucl. Phys. A **(1980)**, Elsevier PII `0375947480900457` ("Large-angle quasi-free scattering in ⁶Li(p,pd)⁴He at 670 MeV") | Hadronic probe of the same distribution: energy-sharing and angular-correlation measurement of the recoil momentum distribution. The literature carries a **factor-two disagreement on its width** — an earlier 590 MeV analysis quoted FWHM ≈ 120 MeV/c, this one finds the data compatible with ≈ 70 MeV/c, *"in agreement with other experimental findings and with the prediction of the cluster model."* That disagreement is itself a useful band. |
| 5 | **⁶Li(α,2α)²H at 700 MeV** | Phys. Lett. B **(1975)**, Elsevier PII `0370269375905791` ("The ⁶Li(α,2α)²H reaction at 700 MeV and the α–d momentum distribution") | Third probe, α-induced. |
| 6 | **Trinucleon knockout ⁶Li(e,e′³H)³He and ⁶Li(e,e′³He)³H** | J. P. Connelly, B. L. Berman, W. J. Briscoe, K. S. Dhuga, A. Mokhtari, D. Zubanov, H. P. Blok, R. Ent, J. H. Mitchell, L. Lapikás, **Phys. Rev. C 57 (1998) 1569–1573** (full text read for this survey) | The *other* cluster decomposition of ⁶Li. Relevant as a **caution**, not a benchmark: the paper reports *"a significant deviation"* between mirror channels at low momentum transfer, i.e. two-step processes contaminate cluster knockout. Any use of 1–5 as a wave-function measurement inherits that reaction-mechanism systematic. |

**What these validate.** The *shape and width* of the α–d relative momentum
distribution — the object `ClusterConfigSampler` and the VMC overlap encode, and
the thing that sets the spectator-momentum spectrum of every tagged event and
therefore the FSI-sensitive region (p_spec ≳ 300 MeV).

**What they do NOT validate.** (a) The **D-wave**. These are unpolarized
cross sections; they constrain |φ₀|² + |φ₂|² summed, not the S–D interference
that b₁ is built from. The 1.93× distance of the undialled tables from George–Knutson's phase-shift-analysis η (C-1, a consistency band on the quadrupole dial) is untouched. (b) The
**absolute normalisation** — spectroscopic factors from cluster knockout carry
the reaction-mechanism systematic item 6 exhibits. (c) **⁷Li → α + t**: item 6 is
⁶Li's ³H/³He decomposition, not ⁷Li's, and the tree's ⁷Li triton remnant model
(C-11, E-14) has no measurement here.

**Obtainability.** All are paywalled journal articles from 1974–1998 with **no
HEPData records**; the momentum distributions are **figures**, digitizable, and
the papers are the only source. The Connelly PRC 57 full text is openly
available through the VU Research Portal (verified, read).

### The one clean data-vs-VMC test that already exists

**L. Lapikás, J. Wesseling, R. B. Wiringa, PRL 82 (1999) 4404,
`nucl-th/9904008`** — ⁷Li(e,e′p)⁶He, missing momentum −70 to 260 MeV/c,
transitions to the ⁶He ground state and first excited state, summed
spectroscopic factor **0.58 ± 0.05 measured vs 0.60 VMC**. The VMC is Wiringa's,
i.e. **the same AV18+UX/UIX family as `data/vmc/`**. This is the closest thing in
existence to an experimental validation of the wave functions E-13 reads, it is
on ⁷Li, and it is a *proton* removal rather than a *cluster* removal — so it
tests the one-body density and the normalisation, not the α–t or α–d overlap.
It costs one number and one sentence to cite, and it converts E-13 from
"transform convention checked" to "the underlying wave-function family has been
compared with data by its own authors, at the 8 % level, on the sister isotope."

---

## 4. Quasi-elastic (e,e′) — the RC tail's dominant knob

`rc.hpp` states the quasi-elastic term is **22 % / 73 % / 99.9 %** of the
radiative tail at x = 0.01 / 0.10 / 0.30, so `qe_kf_gev` and `qe_suppression` are
the tail's dominant knobs and are carried with a flat ±100 % band.

| resource | content | obtainable |
|---|---|---|
| **Quasielastic Electron Nucleus Scattering Archive** (Benhar, Day, Sick, `nucl-ex/0603032`), **⁶Li page** | **133 measured cross-section points**, plain ASCII, columns `Z A E[GeV] θ[deg] ω[GeV] σ[nb/sr/GeV] δσ`. Three settings: **2.5 GeV/12.0°** (40 pts, Q²|peak = 0.258 GeV²), **2.7 GeV/13.8°** (47 pts, 0.389), **2.7 GeV/15.0°** (46 pts, 0.452). Single source: `Heimlich:1973` = **F. H. Heimlich *et al.*, Nucl. Phys. A231 (1974) 509**, "High-energy electron scattering from ⁶Li and ¹²C". | **Downloaded and parsed for this survey**: `discovery.phys.virginia.edu/research/groups/qes-archive/data/6Li.dat`. Machine-readable, free. |
| **R. R. Whitney, I. Sick, J. R. Ficenec, R. D. Kephart, W. P. Trower, Phys. Rev. C 9 (1974) 2230** | The full paper behind the k_F fit: quasi-elastic (e,e′) at **500 MeV, 60°** on nine nuclei **including ⁶Li** — i.e. the large-angle, transverse-dominated kinematics the archive's ⁶Li entry does **not** cover. | Paywalled, no HEPData, **not in the QES archive under ⁶Li**. Digitizable from the paper's figures. |
| **E. J. Moniz, I. Sick, R. R. Whitney, J. R. Ficenec, R. D. Kephart, W. P. Trower, PRL 26 (1971) 445** | The letter; source of `RC_QE_KF_GEV = 0.169`. Verified: 500 MeV, 60°, nine nuclei Li → Pb, Fermi-gas fit. | Already in the tree as an input constant (C-10). |
| **QES archive ⁴He page** | **2729 points**, 8 independent experiments (`Zghiche:1993xg`, `Day:1993md`, `Mccarthy:1976re`, `Meziani:1992xr`, `O'Connell:1987ag`, `Rock:1981aa`, `Sealock:1989nx`, `vonReden:1990ah`), 425 (E, θ) settings. Same ASCII format. | Downloaded and parsed. Free. |

**What the ⁶Li archive data validate.** The Fermi-gas + Pauli-suppression
construction (`de Forest–Walecka`, Q9) against a real ⁶Li quasi-elastic peak, at
three kinematic points, absolutely normalised. That would replace an *input*
(k_F = 0.169 from a 1971 fit) with a **fit residual** — a genuine upgrade, and
the first time the ⁶Li QE tail would be compared with a ⁶Li QE measurement.

**What they do NOT validate.** (a) The kinematic reach: Q² ≈ 0.26–0.45 GeV², a
tiny corner of the eta_A quadrature that runs to q ~ 278 fm⁻¹. (b) The
**longitudinal/transverse separation** — the archive supplies cross sections
only, no separated responses, and explicitly no Coulomb corrections. (c) Anything
tensor: the polarized quasi-elastic tail remains uncomputed for any A = 6 spin-1
nucleus, and `qe_tensor_scale` remains a borrowed magnitude.

---

## 5. Moments, radii, and the α core

| quantity | LiPolGen symbol | best external source, verified | value found |
|---|---|---|---|
| μ(⁶Li) | `LI6_MU_N` = 0.8220473 | **N. J. Stone, Table of Nuclear Magnetic Dipole and Electric Quadrupole Moments** (NNDC mirror, read for this survey) | **+0.8220473(6)** μ_N (AB/D, 1974Be50, ZP 270 173) — exact agreement, digit for digit |
| Q(⁶Li) | `LI6_QUADRUPOLE_FM2` = −0.0818 | Stone (read) gives **−0.00083(8) b = −0.083(8) fm²** (molecular beam, ratio to ⁷Li, CPL 112 1 (1984)); TUNL A = 6 gives −0.0818(17) fm²; Pyykkö gives −0.0806(6) fm² | The three compilations span **−0.0806 … −0.083**, i.e. **± ~2 %** — negligible beside the model's 7.52× discrepancy (C-2). The choice of compilation is **not** a live systematic and the tree can stop worrying about it. |
| μ(⁷Li), Q(⁷Li) | ⁷Li moments block | Stone (read): **μ = +3.256427(2)** μ_N; **Q = −0.0400(3) b = −4.00(3) fm²** (CER, Voelk *et al.*, NPA 530 (1991) 475), with seven independent determinations listed spanning −0.037 to −0.059 b | The modern CER values (−0.0400(3), −0.0400(6)) are 6 σ from the old optical −0.059(8). Also independently determined from **electron scattering**: van Niftrik *et al.*, NPA 174 (1971) 173. **2026-09-06 (task B3): ⁷Li's Q is now IN the tree** as `LI7_QUADRUPOLE_FM2` = **−4.06 fm²**, and it is *not* the CER value of this row — it is TUNL's, `Q = −40.6 ± 0.8 mb` from the A = 7 half of the same Tilley *et al.* NPA 708 (2002) 3 evaluation that supplies `LI6_QUADRUPOLE_FM2`, chosen so both lithium quadrupoles share one document and one convention. The two differ by 1.5 % and move the A = 7 α–t gate ratio 0.858389 ↔ 0.871265; both are recorded in `rc.hpp` and pinned in `run_2026-09-06/phase_B_numbers.md` §B3. |
| r_ch(⁶Li) | behind `LI6_R2_POINT_FM2` | **I. Angeli, K. Marinova, ADNDT 99 (2013) 69–95** (verified); **de Vries ADNDT 36 (1987) 495** Table I (read) | 2.589(39) fm (Angeli) vs **2.54(5) / 2.56(5) / 2.57(10)** fm (three electron-scattering analyses). Angeli is ~1 % higher than the scattering analyses; the tree uses Angeli. |
| η(⁶Li → α+d) | `LI6_ETA_DS_GK` | **E. A. George, L. D. Knutson, Phys. Rev. C 59 (1999) 598** — verified; "Determination of the asymptotic D- to S-state ratio by a restricted phase shift analysis", **η = −0.025 ± 0.006 ± 0.010** | Already C-1 — a consistency band on the quadrupole dial, not a D-wave benchmark. Verified independently here; no better determination found. |
| ⁴He charge form factor (the α core) | `AlphaCoreSource::VmcHe4Density`, `he4.density` | **de Vries Table V has a full sum-of-Gaussians for ⁴He**: rms **1.676(8)** fm, 12 (R_i, Q_i) pairs (read for this survey, from Si82). Underlying data: Frosch *et al.*, PR 160 (1967) 874; **Ottermann *et al.*, NPA 436 (1985) 688**; high-Q²: Camsonne *et al.* (`1309.5297`) | **This is the cleanest, most immediately actionable item in the whole survey.** A published closed-form ⁴He charge form factor exists, in the standard SOG parameterisation, valid over the whole measured q range. It can be transcribed in an hour and used as an **external check on the α core** that today is validated only against ANL VMC. |

### Polarized ³He — the honest scope

**E97-110** (Sulkosky *et al.*, `1908.05709`) and **E06-014** (Flay *et al.*;
Parno *et al.*, PLB 744 (2015) 309) are the closest polarized light-nucleus
analogues, and their value here is **narrow and worth stating as such**: ³He is
spin-½, so it has **no tensor sector at all**. What they exercise is the
*effective-polarization* formalism — P_n ≈ 0.86, P_p ≈ −0.028 from the same
class of ab-initio wave functions — which is the vector-sector analogue of
E-17's VMC polarizations. They validate that a nuclear-structure correction of
that kind, computed from a realistic wave function, survives contact with data
in a *A = 3, J = ½* system. They validate **nothing** about b₁, A_zz, P_zz, or
any rank-2 quantity, and citing them as support for the tensor sector would be
exactly the kind of transfer this survey exists to prevent.

---

## 6. Peripheral, verified, low priority

* **T. Hotta *et al.*, Nucl. Phys. A645 (1999) 492–508**, "Measurement of the
  ⁶Li(e,e′p) reaction cross sections at low momentum transfer"
  (`nucl-ex/9810016`): triple-differential cross sections, E_x = 27–46 MeV, the
  giant-resonance region; DWIA direct-knockout reproduces the data. Relevant only
  to the low-ω end of the RC tail, which is dominated by the quasi-elastic and
  elastic pieces the tree already carries.
* **⁶Li photoabsorption**: B. L. Berman, R. L. Bramblett, J. T. Caldwell,
  R. R. Harvey, S. C. Fultz, PRL 15 (1965) 727; D. D. Faul, B. L. Berman,
  P. Meyer, D. L. Olson, PRC 24 (1981) 849. Cited in the Connelly paper as
  evidence that ⁶Li's photonucleon cross section *"shows no evidence for ²H or
  ⁴He substructures"* and looks like ³He/³H photodisintegration — i.e. it is
  **evidence against** the naive α–d picture at photon energies, and it is the
  clearest published caution about over-reading the cluster model. Bibliographic
  data verified via that reference list only; papers not opened.
* **NMC's own radiative correction used Moniz for lithium.** Verified in the
  ⁶Li/D paper's §3: quasi-elastic suppression from Bernabeu for D and C,
  *"whereas for lithium the result of Moniz"* — the same PRL 26 (1971) 445 that
  supplies `RC_QE_KF_GEV`. A pleasing consistency and worth one sentence if the
  NMC gate is wired.

---

## 7. Priority, and the two things worth doing first

| rank | action | effort | what it buys |
|---|---|---|---|
| 1 | *(Done 2026-09-23: `validation/benchmarks/t3_nmc_li6_over_d.py`, pass.)* Wire **HEPData `ins394050` Table 1** as a gate on `Epps21Ratio`, with the isoscalar denominator of §1.3 and the deuteron caveat of §1.4 stated in the test | ~half a day (24 points, CC0, one JSON download) | Converts the tree's **only** ⁶Li-specific unpolarized-DIS input from "an nPDF fit" to "an nPDF fit that reproduces the one measurement of this nucleus, χ² = 4.63/4". Also puts the ±0.022 external uncertainty on the record, which is currently missing everywhere. |
| 2 | Transcribe **de Vries Table V's ⁴He sum-of-Gaussians** and gate the ⁴He core against it | ~1 hour | The α core's only current validation is ANL VMC (i.e. theory checking theory). This makes it data. |
| 3 | *(Partly done 2026-09-23: the C2 fill-in is excluded inside the model, and the UVa FB zero 2.694 fm⁻¹ is a second number below the window — `validation/benchmarks/t3_li6_charge_ff_fb.py`, recorded FAIL; Coulomb distortion not computed; the window was not moved.)* Resolve **q₀ = 3.0998 vs the measured ⁶Li diffraction minimum at q ≈ 2.83 fm⁻¹** (§2.4) | ~1 day of arithmetic, no new data | Either the [2.9, 3.3] window is defended with the folding argument, or it moves. Right now the assertion window excludes the one measured number. |
| 4 | Cite **Lapikás–Wesseling–Wiringa PRL 82 4404** in E-13's block | ~15 minutes | The VMC family LiPolGen reads has been compared with data by its own author; the tree does not say so. |
| 5 | Add the **QES archive ⁶Li 133 points** as a quasi-elastic-tail benchmark | 1–2 days | Turns `RC_QE_KF_GEV` from an input into a fit residual, at three real ⁶Li kinematic points. |
| 6 | Fit a **⁷Li `HoSpin1FF`** on Suelzle's published HO parameters (a = 1.77(2), α = 0.327) | ~1 day | ⁷Li is the isotope whose shape is published in exactly the shipped functional form. It also gives the ⁶Li starting values their only external comparison. |
| 7 | Digitize **Ent NPA 578 / Mitchell PRC 44** α–d momentum distributions | 3–5 days | First data comparison for the tagged channel's own kinematic distribution. Unpolarized only — the D wave stays untested. |

## 8. What still does not exist, after all of the above

1. **Polarized ⁶Li or ⁷Li anything.** Unchanged.
2. **Any ⁷Li DIS measurement.** The polarized-EMC calculation the tree carries
   (CBT) is *for* ⁷Li and has no ⁷Li unpolarized data under it either.
3. **A separated ⁶Li C2 form factor.** No experiment has done it, and the
   longitudinal data cannot.
4. **A machine-readable ⁶Li elastic form factor.** No HEPData, no SOG in de
   Vries, no archive. Transcription from 1967/1971 journal pages is the only
   route, and this survey could not confirm from outside a paywall that the
   printed tables exist.
5. **Polarized quasi-elastic or elastic radiative tails for A = 6.** Unchanged.
6. **The D-wave of the α–d overlap, from data.** Every measurement in §3 is
   unpolarized. The 1.93× distance of the undialled tables from George–Knutson's phase-shift-analysis η (C-1, a consistency band on the quadrupole dial)
   has no better external anchor than that phase-shift analysis, and this survey found none.
