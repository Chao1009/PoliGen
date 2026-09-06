<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# 02 — External measured data at the nucleon and deuteron level

**Survey date 2026-09-06.** Companion to `00_in_tree_checks.md`, which inventories
what this repository already checks. This file inventories what the *world* has
measured at A = 1 and A = 2 that LiPolGen's kernels must reproduce, verifies for
each whether the numbers are actually obtainable, and states for each exactly
which LiPolGen symbol it would test — and which it would **not**.

Everything below was verified in this session against a primary index (arXiv
abstract page, INSPIRE literature record, or INSPIRE's HEPData `data` collection).
Nothing is quoted from memory. Items that could not be verified are marked
**UNVERIFIED** and say why.

---

## 0. Method, and the one place it is weaker than it looks

**How existence was verified.** arXiv abstracts through the arXiv API
(`export.arxiv.org/api/query`, id_list and search_query); journal reference,
collaboration and record identity through the INSPIRE REST API
(`inspirehep.net/api/literature`).

**How HEPData availability was verified — and the caveat.** `hepdata.net` itself
is behind a Cloudflare interactive challenge from this environment: both `curl`
(browser `User-Agent`) and the `WebFetch` tool get `HTTP 403` on
`/search/` and on `/record/insNNNNNN`. This is the *same* mitigation
`data/vmc/README.md` already records for `phy.anl.gov`. The availability column
below therefore comes from **INSPIRE's `data` collection**
(`inspirehep.net/api/data`), which mirrors every HEPData submission with its
`10.17182/hepdata.NNNNN` DOI, its per-table `…/tN` DOIs, and HEPData's own
`observables:` / `cmenergies:` keywords. For every record listed I resolved the
`literature` back-link to confirm the paper it belongs to, so record ↔ paper is
verified, not guessed. What I could **not** do is read the table *contents*
(bin centres, column names, correlations). Wherever a claim depends on contents
rather than existence, it says so.

**A negative from that index is a real negative, with one qualifier.** INSPIRE's
`data` collection returns *nineteen* HERMES submissions (§1.1 lists why that
matters). A HEPData record that INSPIRE has not ingested would be invisible to
this method; I know of no such case for the collaborations here, but the honest
statement is "not in INSPIRE's HEPData index", not "does not exist on HEPData".

**The frame, restated.** There is no polarized e + ⁶Li or e + ⁷Li DIS data of any
kind, and no other generator produces tensor-polarized lithium. Nothing in this
file changes that. Everything here is a check on a *kernel* — a structure
function, a wave function, an R, a spectral function — that LiPolGen uses on the
way to an observable nobody has measured.

---

## 1. The tensor sector at A = 2 — one datum, and no HEPData record

### 1.1 HERMES b₁ᵈ / A_zz — the only tensor-polarized DIS measurement in existence

| field | value |
|---|---|
| reference | A. Airapetian *et al.* (HERMES), *First measurement of the tensor structure function b₁ of the deuteron*, **PRL 95 (2005) 242001**, `hep-ex/0506018`, DESY-05-077 (INSPIRE 684394) — **VERIFIED** (arXiv abstract + INSPIRE record) |
| observable | A_zz and b₁ᵈ, six x bins, per nucleon |
| kinematics | paper range 0.01 < x < 0.45, 0.5 < Q² < 5 GeV²; bin means (as transcribed in-tree, §1.2) x = 0.012 / 0.032 / 0.063 / 0.128 / 0.248 / 0.452 at Q² = 0.51 / 1.06 / 1.65 / 2.33 / 3.11 / 4.69 GeV² |
| **HEPData** | **NONE.** INSPIRE's `data` collection holds **19** HERMES submissions; the b₁ paper is not among them. (The 19 include HERMES g₁ — `10.17182/hepdata.11211` — the unpolarized F₂ᵈ JHEP paper, GDH, ρ⁰/ω/φ production, SIDIS multiplicities, hadronization.) Obtainability is therefore **table-in-paper**: Table II of the PRL. |
| in-tree already? | The *numbers* are, in prose: `docs/open_items/run_2026-09-03/phase_A_miller_normalisation.md` §2.3 carries all six rows (x, Q², A_zz×10², b₁×10²) plus the (11.20 ± 5.51 ± 2.77)×10⁻² first bin and its five siblings, and `rc.hpp:272` uses the measured \|A_zz\| = 1.06×10⁻² at ⟨x⟩ = 0.012. The **data are not a file and not an assertion**: `tests/test_sf.cpp` "the digitized Miller b1 curve…" pins the *digitization* of Miller's curve and states the 0.97 σ comparison only in a comment (`00_in_tree_checks.md` row C-7). |
| **validates** | If wired as an assertion: `MillerB1`/`toy_b1` and `CdksB1`/`b1_convolution` evaluated per nucleon at the six HERMES (x, Q²) against a measured b₁ᵈ with its published errors; `azz()` (`asymmetries.hpp:120`) through Cosyn Eq. (27) A_zz = −(2/3) b₁/F₁; and, decisively, the **per-deuteron ÷ per-nucleon factor of 2** that is currently an open author decision (`constants.hpp:89`, `B1_MILLER_TABLE_TO_PER_NUCLEON`). |
| **does NOT validate** | Anything about ⁶Li or ⁷Li. b₁ᵈ is A = 2, spin 1 from one S–D pair; the ⁶Li tensor structure function is built by `Li6B1` / `b1_li6_from_deuteron` from an α + d cluster convolution whose *input* is b₁ᵈ. Reproducing the input does not test the convolution. It also does not validate the low-x rise as *shadowing*: five of six bins are consistent with zero at ≲ 1 σ and the two lowest carry the whole signal. |
| kinematic reach vs EIC | Poor. The EIC window this generator runs (x ≥ 0.001, Q² ≥ 1) contains **five** of the six bins (the x = 0.012 bin sits at Q² = 0.51, *below* the cut) and **no** bin below x = 0.012. The b₁ region the EIC exists to open — x ≲ 10⁻² at Q² ≳ 1 — is entirely unmeasured. |
| effort | **Low.** Six rows, already transcribed. Ship as `data/hermes_b1_2005.dat` (or a `DigitizedTable`) with stat and syst columns, then one TEST_CASE that reports χ² per configuration. |
| **priority** | **HIGH.** This is the single highest-value external item in the whole survey: it is the only measurement of the class of observable LiPolGen exists to produce, and it is presently *not asserted anywhere*. |

### 1.2 What the second measurement will be, and that it does not exist yet

* **JLab E12-13-011** ("The Deuteron Tensor Structure Function b₁", Hall C, PAC-40).
  **Approved; not run; no data.** Verified indirectly: the JLab proposal PDF
  (`jlab.org/exp_prog/proposals/13/PR12-13-011.pdf`) and the 2023 jeopardy
  document are live, and the target R&D review of Aug 2022 converted a conditional
  to a full approval. The strongest primary statement I verified is Cosyn, Roldán
  Tomei, Sosa, Zec, **EPJ A 61 (2025) 83**, `arXiv:2410.12764` (already in
  `refs/`, entry `cosyn2025-tensor-polarization-options`), whose abstract opens
  *"In the near future, the Jefferson Lab b₁ experiment will provide **the second**
  measurement of tensor polarized asymmetries in inclusive DIS on the deuteron"* —
  i.e. as of that paper there is exactly one. **Obtainable: does-not-exist (yet).**
* **CLAS12 RG-C Approved Analysis**, `arXiv:2502.20044` (v3, 2025) — proposes the
  *first* SIDIS extraction of tensor structure functions from existing ND₃ data,
  plus an inclusive b₁. A proposal, not a result. **does-not-exist (yet).**
* **JLab Hall D frozen-spin tensor target**, Dalton, Deur, Keith,
  `arXiv:2504.21177` (2025) — a proposal for coherent ρ photoproduction on a
  tensor-polarized deuteron. Photoproduction, not DIS. **does-not-exist (yet).**
* **Review of the programme**: `arXiv:2506.04506`, EPJ A 61, 81 (2025),
  *Experimental Study of Tensor Structure Function of Deuteron* — **VERIFIED**;
  useful as a citable status statement, carries no data.

### 1.3 The near-miss: "the Hall C A_zz on the deuteron" is not what it sounds like

The 2017 tensor-asymmetry paper is **A. DeGrush *et al.*, PRL 119, 182501 (2017)**,
`arXiv:1707.03028`, *Measurement of the Vector and Tensor Asymmetries at Large
Missing Momentum in Quasielastic (e⃗, e′p) Electron Scattering from Deuterium* —
**MIT-Bates / BLAST**, not JLab Hall C, and **quasielastic**, not DIS:
0.1 < Q² < 0.5 (GeV/c)², missing momentum to 500 MeV/c. **VERIFIED** (arXiv
abstract). It measures A_d^T sensitive to *np final-state interactions* and the
D-wave, not a partonic b₁.

**Where it is genuinely useful, and where it is not.** It is a real external
constraint on the **D-state at high k** and on **FSI** — i.e. on
`ClusterWaveSource::VmcAV18` / `Hulthen` at k ≳ 300 MeV/c and on
`GlauberFsiWeight` — in a regime where LiPolGen's tagged channel lives. It is
**not** a check on `toy_b1`, `b1_convolution`, `DeuteronConvolutionB1`, `azz()`,
or any structure function. Q² ≤ 0.5 GeV² is an order of magnitude below the
generator's DIS floor. **HEPData: not in INSPIRE's index.** Obtainability:
digitizable figures. **Priority: LOW** (wrong process, wrong Q²; the FSI content
is better attacked at §5.4).

A pre-HERMES tensor limit in *DIS* does not exist. Searching INSPIRE for tensor
asymmetries before 2005 returns only quasielastic / photodisintegration
electro-disintegration work (e.g. AIP Conf. Proc. 334 (1995) 768). HERMES's own
title — *First measurement* — is accurate.

---

## 2. Vector-polarized deuteron: g₁ᵈ and A₁ᵈ

These test `ToyG1`/`LhapdfG1` (`a1p`, `a1n`, `g1p`, `g1n`, `g1_nucleus`) and
`g2_ww`, i.e. the **vector** half of the polarized kernel — the half that is *not*
the reason this generator exists, but that carries the rates and the systematic
floor for every tensor asymmetry. All records below were verified by resolving the
INSPIRE data record's `literature` back-link to its paper.

| # | measurement | HEPData | tables | kinematics | LiPolGen symbol tested |
|---|---|---|---|---|---|
| D-1 | **HERMES**, *Precise determination of g₁ of the proton, deuteron and neutron*, **PRD 75 (2007) 012007**, `hep-ex/0609039` | `10.17182/hepdata.11211` | **23** (observables A1, G1, ASYM, CORR) | 0.0041 ≤ x ≤ 0.9, 0.18 ≤ Q² ≤ 20 GeV² | `ToyG1::g1_nucleus` on `DEUTERON()`, `a1p`/`a1n`, and the isoscalar combination |
| D-2 | **COMPASS final**, PLB 769 (2017) 34, `arXiv:1612.00620` | `10.17182/hepdata.78374` | 6 | 1 < Q² < 100 GeV², 0.004 < x < 0.7, W > 4 GeV | same; **the best overlap with the EIC window of any g₁ᵈ set** |
| D-3 | **COMPASS**, PLB 647 (2007) 8, `hep-ex/0609038` | `10.17182/hepdata.48555` | 1 | superseded by D-2 | — |
| D-4 | **COMPASS**, PLB 612 (2005) 154, `hep-ex/0501073` | `10.17182/hepdata.48552` | 13 | 1 < Q² < 100 GeV², 0.004 < x < 0.7 | superseded by D-2 |
| D-5 | **COMPASS low-x**, PLB 647 (2007) 330, `hep-ex/0701014` | `10.17182/hepdata.48534` | 1 | **Q² < 1 GeV², 4×10⁻⁵ < x < 2.5×10⁻²** | outside the DIS cut, but the **only** g₁ᵈ below x = 4×10⁻³; A₁ᵈ and g₁ᵈ consistent with zero throughout |
| D-6 | **SLAC E143**, PRD 58 (1998) 112003, `hep-ph/9802357` | `10.17182/hepdata.22265` | **33** (G1, G2, A1, A2) | beam 29.1 / 16.2 / 9.7 GeV; the g₂ and A₂ set is at 29.1 GeV | `g2_ww` — E143's own conclusion is that g₂ᵈ is described by Wandzura–Wilczek, which is *exactly* what `g2_ww` implements |
| D-7 | **SLAC E155**, PLB 463 (1999) 339, `hep-ex/9904002` | `10.17182/hepdata.41630` | 4 | **0.01 < x < 0.9, 1 < Q² < 40 GeV²** | `g1_nucleus`; the widest-Q² deuteron g₁ set |
| D-8 | **SLAC E155x**, PLB 553 (2003) 18 | `10.17182/hepdata.27033` | 7 | g₂ and A₂, p and d | `g2_ww` twist-3 residual |
| D-9 | **CLAS EG1b**, Guler *et al.*, **PRC 90 (2014) 025212**, `arXiv:1404.6231` | `10.17182/hepdata.64411` | 4 (G1, F1) | low Q², resonance-dominated | marginal for the EIC; useful only as a low-Q² bound |
| D-10 | **CLAS EG4**, PRL 120 (2018) 062501, `arXiv:1711.01974` | not in INSPIRE's index | — | very low Q² moments | **LOW**; outside the generator's domain |

**What §2 validates:** that the *vector* polarized structure function the generator
convolutes into ⁶Li/⁷Li (`g1_nucleus` = Z·P_p·g₁ᵖ + N·P_n·g₁ⁿ, times the medium
ratio) starts from a correct free-nucleon g₁ and a correct A = 2 combination; and
that `g2_ww` is the right twist-2 relation to be using.

**What §2 does NOT validate:** the tensor sector at all. g₁ᵈ is blind to b₁ by
construction — a tensor-polarized target with P_z = 0 has no g₁ signal, and the
vector-polarized measurements above have no A_zz sensitivity. It also does not
validate the **effective polarizations** `eff_pol_p` / `eff_pol_n`
(`beams.hpp:162`) that carry ⁶Li's α+d structure: those are VMC sums, and D-1…D-10
constrain only the nucleon-level g₁ they multiply.

**A pointed aside.** E155 (D-7) and COMPASS (D-2…D-5) both took their deuteron data
on **⁶LiD** targets, and both model ⁶Li as α + d to compute the dilution factor —
the same cluster picture as `li6_alpha_channel`. That is a *conceptual*
corroboration of the α+d picture from two experiments' target bookkeeping, not a
measurement of anything LiPolGen computes, and it must not be quoted as one.

---

## 3. Unpolarized F₂ᵈ, F₂ᵖ and F₂ⁿ/F₂ᵖ

These test `ToyF2::f2p`/`f2n`/`f2n_over_f2p` (an explicitly labelled TOY,
"anchored by eye", `sf.hpp:114`), and the two real backends `LhapdfSF` (CT18NLO)
and `MstwSF` (MSTW2008 LO, read from PYTHIA 8's `pdfdata/mstw2008lo.00.dat`).

| # | measurement | HEPData | tables | kinematics | notes |
|---|---|---|---|---|---|
| F-1 | **NMC**, Arneodo *et al.*, *F₂ᵖ, F₂ᵈ and σ_L/σ_T*, **NPB 483 (1997) 3**, `hep-ph/9610231` | `10.17182/hepdata.32752` | **33** (F2, R, SIG) | **0.002 < x < 0.60, 0.5 < Q² < 75 GeV²**; R for 0.002 < x < 0.12 | the best single F₂ᵈ set for the EIC window; carries **both** F₂ and R |
| F-2 | **NMC**, *Accurate measurement of F₂ᵈ/F₂ᵖ and R_d − R_p* | `10.17182/hepdata.32750` | 23 | µp/µd, 29.7 GeV c.m. | the ratio, with its own R difference |
| F-3 | **NMC**, *The ratio F₂ⁿ/F₂ᵖ in deep inelastic muon scattering* | `10.17182/hepdata.32955` | 16 | µ beam | **directly `ToyF2::f2n_over_f2p`**, which is the single-argument, Q²-frozen (Q² = 10) hook that `LhapdfSF`/`MstwSF` reproduce verbatim |
| F-4 | **BCDMS**, *High statistics F₂ᵈ and R*, **PLB 237 (1990) 592** | `10.17182/hepdata.6191` | **87** (F2, R) | high x, high Q² | the high-x anchor |
| F-5 | **SLAC**, Whitlow *et al.*, *Combined analysis of SLAC ep and ed* , PLB 282 (1992) 475 | `10.17182/hepdata.2721` | 44 | fixed target | |
| F-6 | **SLAC**, Whitlow, SLAC-357 (thesis) — *Deep inelastic structure functions from e on H, D, Fe* | `10.17182/hepdata.2722` | **140** | fixed target | the largest single F₂ table set |
| F-7 | **HERMES**, *Inclusive inelastic e±  on unpolarized H and D*, **JHEP 05 (2011) 126**, `arXiv:1103.5704` | `10.17182/hepdata.66147` | 3 | HERMES kinematics, √s = 7.2 GeV | F₂ᵈ **in exactly the kinematics of the b₁ measurement** — this is the F₁ᵈ that HERMES's own Eq. (5) divides by |
| F-8 | **MARATHON**, Abrams *et al.*, **PRL 128 (2022) 132003**, `arXiv:2104.05850` | **not in INSPIRE's HEPData index** | — | **0.19 < x < 0.83** | the modern F₂ⁿ/F₂ᵖ from ³H/³He mirror symmetry; obtainability **table-in-paper / digitizable-figure** |
| F-9 | **CLAS BONuS**, Baillie *et al.*, PRL 108 (2012) 142001, `arXiv:1110.2770`; Tkachenko *et al.*, PRC 89 (2014) 045206, `arXiv:1402.2477` | **neither in INSPIRE's HEPData index** | — | 0.65 < Q² < 4.52 GeV², 0.2 < x < 0.8, spectator p below 100 MeV/c at > 100° | see §5.5 — this is the **spectator-tagging** channel, not just an F₂ⁿ source |

**F-7 deserves special emphasis.** `constants.hpp:89` records that HERMES's
published b₁ᵈ is per nucleon because their Eq. (5) divides by an F₁ᵈ built from
F₂ᵈ = (F₂ᵖ + F₂ⁿ)/2, "confirmed by inverting their own Table II in all six bins".
F-7 is HERMES's *own* F₂ᵈ measurement on HEPData, in the same apparatus and
kinematics. Wiring F-7 would turn that inversion from an argument in a comment
into a check against a published table — and it is the cleanest available way to
settle the factor-of-2 in `B1_MILLER_TABLE_TO_PER_NUCLEON`.

**What §3 does NOT validate.** None of it constrains the *nuclear* F₂:
`NuclearF2::f2a`, `unpolarized_emc_ratio`, `cbt_polarized_emc_ratio`. Those need
A > 2 data (the next file in this series). And F-1…F-9 are all fixed-target: none
reaches x < 2×10⁻³, so the EIC's low-x decade is covered only by HERA (§4.4).

---

## 4. R = σ_L/σ_T — where the tree is weakest relative to what exists

`r_sigma_lt` (`sf.hpp:46`) is **the default R of the whole program** and is
declared in its own docstring to be "a placeholder of the right magnitude …
**not a fit**": 0.18/(1 + Q²/50). `r1998` (`src/core/sf.cpp:14–95`) is the real
thing — Abe *et al.* Table II, three six-parameter forms averaged, with
`r1998_spread` and `r1998_fit_error`. `b1_nuclear.hpp:487–514` records that the
CDKS A = 2 gate deliberately runs `r_sigma_lt` in one slot and `r1998` in the
other, and that (1 + r1998)/(1 + r_sigma_lt) = 1.088 at x = 0.1, Q² = 2.5.

**The gap.** `tests/test_sf.cpp:72–86` tests R1998 for *internal* consistency —
that it is the average of its three forms, that the clip works, and a regression
pin `r1998(0.1, 10) = 0.12654600717127532` at rtol 1e-12. **Nothing compares it
to a measured R.** The fit is transcribed correctly, and that is all that is
established.

| # | measurement | HEPData | tables | kinematics | what it would test |
|---|---|---|---|---|---|
| R-1 | Whitlow *et al.*, *Precise extraction of R from a global analysis of SLAC ep and ed*, **PLB 250 (1990) 193** — this is **R1990** | `10.17182/hepdata.29555` | **13** (R, F2) | SLAC fixed target | `r1998` against the *previous* world fit's data; also the closest thing to a reference for `sf.hpp:43`'s "simplified R1990-like magnitude" claim about `r_sigma_lt` |
| R-2 | Dasu *et al.* (SLAC **E140**), *Kinematic and nuclear dependence of R*, **PRD 49 (1994) 5641** | `10.17182/hepdata.22468` | **31** | 0.2 ≲ x ≲ 0.5, 1 ≲ Q² ≲ 10 GeV² | the largest single R data set feeding R1998; also carries R_A − R_D, the nuclear-dependence handle |
| R-3 | **NMC** R, inside F-1 | `10.17182/hepdata.32752` | (within 33) | **0.002 < x < 0.12** | the **only** R data at NMC's low x; R1998's support starts at x = 0.005 (`R1998_X_MIN`), so this is where the fit is thinnest and the tree clips |
| R-4 | Abe *et al.* (E143), *Measurements of R for 0.03 < x < 0.1 and fit to world data*, **PLB 452 (1999) 194**, `hep-ex/9808028` — the **source of `r1998`** | **not in INSPIRE's HEPData index** | — | fit support 0.005 < x < 0.86, 0.5 < Q² < 130 GeV² (as encoded in `sf.hpp:53–56`) | the paper (PDF already in `PolarizedLithiumSim/refs/hep-ex_9808028.pdf`, dict entry `abe1999-r1998-world-fit`) carries its own measured R points in tables — **table-in-paper** |
| R-5 | **H1**, *Measurement of the proton structure function F_L at low x*, **PLB 665 (2008) 139**, `arXiv:0805.2809` | `10.17182/hepdata.45340` | 9 (FL) | **12 < Q² < 90 GeV², 0.00024 < x < 0.0036** | **the R data at EIC x.** R = F_L/(F₂ − F_L). This is the only measurement in the whole survey that sits inside the EIC's low-x window |
| R-6 | **H1**, *Determination of F_L at low x* | `10.17182/hepdata.44694` | 3 (FL, F2, d²σ) | low x | corroborates R-5 |
| R-7 | **H1**, *Inclusive ep at √s = 225, 252 GeV and F_L* | `10.17182/hepdata.62536` | **52** (SIG, F2, FL) | high Q² | the widest H1 F_L set |

**Priority: HIGH for R-5 (and R-1/R-2 as the fixed-target anchor).** Two concrete
consequences follow from wiring them:

1. **`r1998` is being *extrapolated* over most of the EIC window.** Its declared
   support stops at x = 0.005 and Q² = 130 GeV²; the generator runs to x = 0.001
   and well past Q² = 130. `r1998(..., clip = true)` freezes at the boundary, so
   the code is honest about it, but nothing measures what that costs. R-5 gives a
   measured R at x ≈ 3×10⁻⁴ that the frozen fit can be compared against.
2. **`r_sigma_lt` is the *default*, and it is not a fit.** Every quantity that
   goes through `f1_from_f2` = F₂/(2x(1+R)) — every F₁, every b₁/F₁ ratio, the
   whole A_zz normalisation — carries it. A one-page comparison of `r_sigma_lt`
   against R-1/R-2/R-3/R-5 would put a number on the resulting bias where today
   there is only the 1.088 spot value at one (x, Q²).

**What §4 does NOT validate.** R on a *nucleus*. R_A − R_D is in R-2 for A > 2,
but there is no R for ⁶Li anywhere, and no tensor R at all: the L/T separation of
the *tensor* structure functions (Cosyn Eq. (16)–(17)'s F_UTLL,L vs F_UTLL,T) has
never been measured and is a pure model choice in `cosyn_tensor_sfs`.

---

## 5. Deuteron structure: the A = 2 machinery under the tagged channel

### 5.1 The wave function itself — CD-Bonn and AV18

`CdBonnWave` (`cluster.hpp:139–276`, `src/core/cluster.cpp:305–441`) is a
transcription of **R. Machleidt, *The high-precision, charge-dependent Bonn
nucleon–nucleon potential (CD-Bonn)*, PRC 63 (2001) 024001, `nucl-th/0006014`,
Table XX** — **VERIFIED** (INSPIRE 528472). `data/vmc/deuteron/fdeut.av18` is
Wiringa's AV18 deuteron (11 028 rows, header `deuteron for vij = Argonne v18`).

**Live numbers from this checkout** (`./build/lipolgen_tests -tc="cluster: the
CD-Bonn deuteron coefficients (Machleidt Table XX)"`, run for this survey, passed):

```
CD-Bonn P_D  = 4.85621 %      (Machleidt published 4.85 %)
CD-Bonn A_S  = 0.88473        (published 0.8846)
CD-Bonn eta  = 0.0255714      (published 0.0256)
norm = 1 (S 0.951438, D 0.0485621);  on the fdeut grid, int k^2(u^2+w^2)dk = 0.999978
```

This is already an external check — but only against **Machleidt's own table**,
i.e. a transcription gate (`00_in_tree_checks.md` row C-6 says exactly this).
Two of those four numbers are secretly *empirical*, and that is where an external
resource adds something the transcription cannot (§5.2).

The same run reported the A = 2 b₁ gate's configuration dependence, which is the
context for everything in §3 and §4 — `-tc="b1_nuclear T1 gate checklist item 5v:
CD-Bonn + MSTW2008 LO -- gate condition 3"`, passed:

```
G3b: max|x b1| over [0.10,0.80] = 0.00108558 ; ratio to the digitized CDKS peak
     = 1.00034   (AV18 on the same row is 0.843243)
G3c: int b1 dx (x >= 0.01) = 0.00044858 ; digitized CDKS 0.0004592
     (AV18 2.24896e-4)
```

i.e. the *wave function* moves the A = 2 gate by 19 %, and the *nucleon F₁*
(toy vs MSTW2008 LO) moves it by roughly a factor of two (`00_in_tree_checks.md`
§3: shipped ToyF2 gives 0.440 and fails the window). Both of those levers are
exactly what §3, §4 and §7 constrain — which is the argument for wiring them.

### 5.2 The asymptotic D/S ratio η_d — a measurement hiding inside a published constant

| field | value |
|---|---|
| reference | **N. L. Rodning, L. D. Knutson, *Asymptotic D-state to S-state ratio of the deuteron*, PRC 41 (1990) 898** (INSPIRE 312974) — **VERIFIED**; and the earlier PRL 57 (1986) 2248 (INSPIRE 1476759) |
| value | η_d = 0.0256(4), from d + ⁴He sub-Coulomb tensor analysing powers — **the same technique and the same group** as George & Knutson's η(⁶Li → α+d) = −0.025 ± 0.006 ± 0.010, which `00_in_tree_checks.md` row C-1 already gates |
| LiPolGen symbol | `CdBonnWave::eta()` (`cluster.hpp:267`), `CdBonnWave::a_s()`; indirectly `deuteron_av18_p_d()` and `P_D_DEUTERON` = 0.045 (`beams.hpp:98`) |
| **validates** | that the S–D mixture of the deuteron wave function the tagged channel and the A = 2 b₁ gate both ride on is **empirically right at the 0.2 % level** — CD-Bonn's 0.0255714 against the measured 0.0256(4) is 0.07 σ. Not a fit residual: η is an *asymptotic* observable, so it constrains the tail of w(r) directly |
| **does NOT validate** | P_D. The D-state probability is **not an observable** — it is representation-dependent — and CD-Bonn's 4.856 % vs AV18's ≈ 5.76 % vs the code's own scenario `P_D_DEUTERON` = 0.045 differ by more than any error bar. η pins the *asymptotic* D/S; it says nothing about w(k) at the k ≳ 300 MeV/c that dominates the tagged channel's high-k tail. It also does not validate the ⁶Li η, which the tree already knows is off by 1.93× (row C-1) |
| obtainability | **table-in-paper** (a single quoted number; the paper is not on arXiv — 1990 PRC). Effort: trivial — it is already quoted, via Machleidt, at `test_cluster.cpp:386` |
| priority | **MEDIUM.** The number is already in the test; what is missing is the *attribution* that 0.0256(4) is a measurement, not a model output. That is a two-line comment change, and it upgrades row C-6 from "transcription gate" to "transcription gate **plus** one genuine empirical contact" |

### 5.3 Deuteron elastic form factors — a real gap in the radiative-correction tail

| field | value |
|---|---|
| references | **D. Abbott *et al.* (JLab t₂₀), *Measurement of tensor polarization in elastic e–d scattering at large momentum transfer*, PRL 84 (2000) 5053, `nucl-ex/0001006`** — **VERIFIED**, HEPData **`10.17182/hepdata.40433`, 5 tables** (observables FORMFACTOR, POL), Q² = 0.66–1.7 (GeV/c)²; and **D. Abbott *et al.*, *Phenomenology of the deuteron electromagnetic form factors*, EPJ A 7 (2000) 421, `nucl-ex/0002003`** — **VERIFIED**, three world-data parametrizations of G_C, G_M, G_Q over 0–7 fm⁻¹ |
| also | `10.17182/hepdata.31372` (*Measurement of T₂₀ in elastic e–d scattering*, 1 table) |
| LiPolGen symbol | `Spin1ElasticFF` and its implementations `HoSpin1FF` / `TabulatedSpin1FF` (`rc.hpp:357–664`), reached through the deuteron control channel; `TEST_CASE("T10: the DEUTERON fixes the (F_c, F_m, F_q) normalisation")` and `T8(a)` in `tests/test_rc.cpp` |
| **the finding** | `tests/test_rc.cpp:374–397` says it outright: *"The SHAPE is a stand-in for POLRAD's `ffdeu` — every T8(a)/T10 gate is shape-blind"*. The deuteron form factor in the RC elastic tail is a harmonic-oscillator shape normalised on the measured moments (μ = 0.857406, Q_d = 0.2859 fm²), with `a_fm = 1.593, alpha = 0.02` chosen only to reproduce ⟨r²⟩_point. **The world data for G_C, G_M, G_Q exist, are parametrized in closed form by Abbott EPJ A 7, and are not used.** |
| **validates** (if wired) | the *magnitude and shape* of σ^el for the deuteron control channel, hence the elastic-tail fraction that `RcModel` subtracts — the one place where the deuteron channel's radiative correction can be checked against something other than POLRAD's own compiled numbers (T8(c')) |
| **does NOT validate** | the **⁶Li** form factors, which is where the RC tail's real uncertainty lives: `HoSpin1FF` for ⁶Li runs a ±100 % `fq_scale` band precisely because `eps_b0 = −0.08` corresponds to a quadrupole 11.4× the measured one (`00_in_tree_checks.md` row C-5). A correct deuteron F_c/F_m/F_q does not constrain that band. It also does not validate `nucleon_ff` (§8) |
| obtainability | **hepdata** for t₂₀ (`10.17182/hepdata.40433`); **table-in-paper** for the Abbott closed-form parametrizations (three parameter sets, EPJ A 7 (2000) 421) |
| priority | **MEDIUM.** Contained, well-defined, and it removes a stand-in. But it improves a *control channel*, not the physics channel |

### 5.4 Deuteron momentum distribution from d(e,e′p)n

| field | value |
|---|---|
| reference | **W. U. Boeglin *et al.*, *Probing the high momentum component of the deuteron at high Q²*, PRL 107 (2011) 262501, `arXiv:1106.0275`** — **VERIFIED** (JLab Hall A) |
| observable | reduced d(e,e′p) cross sections at Q² = 3.5 (GeV/c)², binned in neutron recoil angle θ_nq, giving **missing-momentum distributions up to 0.55 GeV/c**; in 35° ≤ θ_nq ≤ 45° FSI are predicted small and the reduced cross section is a direct read of the momentum distribution |
| LiPolGen symbol | `ClusterWaveSource::Hulthen` vs `VmcAV18` (`cluster.hpp`), `Wave::radial`, `TaggedModel`'s k-grid and `n_of_kc`, `momentum_density`; and — because the 35°–45° selection is exactly a *minimal-FSI* window — `GlauberFsiWeight` (`fsi.hpp`) |
| **validates** | that the S+D Hulthén at `P_D_DEUTERON` and the AV18 table give the right **shape and magnitude of n(k) out to 550 MeV/c**, which is the region the tagged spectator spectrum populates. Today `test_tagged.cpp` T25 compares "Hulthén against AV18" — *two models against each other*; Boeglin makes it model-against-data |
| **does NOT validate** | the ⁶Li α–d relative-motion wave function (a different object with a Coulomb tail and a node structure the deuteron has not got); nor the *polarized* momentum distribution — d(e,e′p) is unpolarized, so `eff_pol_p`/`eff_pol_n` and the D-state dilution `vector_dilution_of` are untouched; nor the tagged channel's **light-cone** variables, since PRL 107 is at fixed Q² in the target rest frame |
| **HEPData** | not in INSPIRE's index. Obtainability: **digitizable-figure** (the paper's reduced-cross-section vs p_miss panels) |
| priority | **MEDIUM–HIGH.** It is the only measurement that turns T25 from model-vs-model into model-vs-data, and its FSI-suppressed window is a rare clean handle |

### 5.5 Spectator tagging on the deuteron — the closest experimental analogue to the tagged channel

| field | value |
|---|---|
| references | **N. Baillie *et al.* (CLAS BONuS), PRL 108 (2012) 142001, `arXiv:1110.2770`** and **S. Tkachenko *et al.* (CLAS BONuS), PRC 89 (2014) 045206, `arXiv:1402.2477`** — both **VERIFIED** |
| what it is | DIS on deuterium with a **backward, low-momentum spectator proton** detected (p_s < 100 MeV/c, θ_pq > 100°), 0.65 < Q² < 4.52 GeV², 0.2 < x < 0.8 — i.e. **exactly LiPolGen's tagged topology at A = 2**, with the same physics argument (a slow backward spectator selects a nearly on-shell struck nucleon and suppresses FSI) |
| LiPolGen symbol | `TaggedChannel`, `TaggedModel` (the (α_s, p_T) grid, `struck_populations`, `population_integrated`, the dilutions), `TaggedSampler`, `boost_spectator`, and `GlauberFsiWeight`'s "never-rescattered" limit (`tagged.hpp:372` names Cosyn–Weiss and BeAGLE for this) |
| **validates** | the *shape of the tagged spectator spectrum* against a measured one, and the on-shell-extrapolation logic. Crucially it tests the **direction of the FSI correction**: BONuS's own systematic is dominated by exactly the rescattering `GlauberFsiWeight` models |
| **does NOT validate** | anything tensor, anything polarized, anything at ⁶Li. Nor the *forward* spectator kinematics the EIC's far-forward detectors see: BONuS tags **backward** in the target rest frame at p_s < 100 MeV/c; the EIC tags a boosted α or d in the B0/off-momentum/ZDC acceptance. The two are the same physics in different frames, but the acceptance-driven systematics do not transfer |
| **HEPData** | **neither paper is in INSPIRE's HEPData index.** Obtainability: **digitizable-figure**, or **on-request** from CLAS |
| priority | **MEDIUM.** High conceptual value, moderate wiring cost, and the frame mismatch limits what a numerical comparison would mean |

---

## 6. A = 3: the triton spectral function

| field | value |
|---|---|
| reference | **C. Ciofi degli Atti, S. Simula, *Realistic model of the nucleon spectral function in few- and many-nucleon systems*, PRC 53 (1996) 1689, `nucl-th/9507024`** — **VERIFIED** (INSPIRE 397215, arXiv abstract read) |
| how the tree uses it | `CiofiSimulaTriton` (`src/core/triton_sf.cpp`): n₀ from their **Table A.1** (A = 2, 3, 4 — A = 3 is the *"³He (proton)"* column), n₁ from **Eq. (76)** for A = 3 and **Eq. (75) + Table A.3** for A = 4. Reached by `--triton-sf ciofi-simula` for the ⁷Li α-tag's T1 breakup |
| **what it was fit to — the honest answer** | **Not to data.** Their abstract states the model is a *factorised ansatz* whose spectral function is validated by comparing the convolution result **against many-body calculations** — "The Spectral Functions of ³He and infinite nuclear matter resulting from the convolution formula and from many-body calculations are compared". The A = 3 parameters describe a **Faddeev/variational ³He momentum distribution**, and the paper's contact with experiment is through *applications* to inclusive and exclusive cross sections, not through a fit. So `cs_n0_terms(3)` / `cs_n1_terms(3)` are a **fit to a calculation**, and the LiPolGen use is theory-against-theory |
| **the A = 3 mismatch nobody should paper over** | the tree draws a **triton** (³H) spectral function from CS's **³He proton** column. That is isospin-mirror reasoning, and it is defensible, but it is an assumption, not a transcription. The paper offers no ³H set |
| external data that WOULD test it | **M. Rvachev *et al.*, *The quasielastic ³He(e,e′p)d reaction at Q² = 1.5 GeV² for recoil momenta up to 1 GeV/c*, PRL 94 (2005) 192302, `nucl-ex/0409005`** — **VERIFIED**. Measured cross sections to p_miss = 1000 MeV/c; the paper's own conclusion is that a **variational ³He ground-state wave function with three-body forces describes p_miss ≲ 150 MeV/c well, and that 150–750 MeV/c is dominated by strong FSI**. That is precisely the regime `CiofiSimulaTriton::n_of_k` and `p_two_body` claim to describe |
| **validates** (if wired) | the **low-k** ( ≲ 150 MeV/c) normalisation and shape of `n_of_k` for A = 3, and — from the FSI-dominated middle region — a *bound* on how far the plane-wave spectral function can be trusted before the T1 breakup needs an FSI weight of its own |
| **does NOT validate** | the removal-energy dependence `e_rel_density` / `e_max` (the (e,e′p) data at high p_miss are FSI-contaminated, so they cannot cleanly separate S(k, E)); nor the ³H ↔ ³He mirror step; nor anything about the *α + t* cluster overlap that gates whether a triton is drawn at all |
| **HEPData** | not in INSPIRE's index for either. Obtainability: **digitizable-figure**. MARATHON (F-8) is the modern ³H/³He pair but measures F₂ⁿ/F₂ᵖ, not a spectral function |
| priority | **LOW–MEDIUM.** The T1 triton channel is a small branch, and the honest documentation fix (§6's "fit to a calculation, and a mirror assumption") is worth more than the numerical comparison |

Also verified and available if the ³He side is ever needed on the *polarized* axis:
**HERMES ³He g₁ⁿ**, PLB 404 (1997) 383, `hep-ex/9703005`, HEPData
`10.17182/hepdata.44586`; and **JLab E06-014**, Posik *et al.*, *A precision
measurement of the neutron twist-3 matrix element d₂ⁿ*, `arXiv:1404.4003`,
PRL 113 (2014) 022002 — 0.25 ≤ x ≤ 0.90, ⟨Q²⟩ = 3.21 and 4.32 GeV², polarized ³He,
**not in INSPIRE's HEPData index**. E06-014 tests `g2_ww`'s twist-3 *residual*
(d₂ is by construction the piece `g2_ww` omits) — a genuine but narrow check on a
function that is currently only self-tested against its own analytic power-law
limit (`test_sf.cpp:147`).

---

## 7. The nucleon PDF backends — pedigree, not data

These are **generators/fits**, not measurements. Their value is that each carries
its own published validation against exactly the data of §2–§4, so attaching one
imports that validation wholesale — and *not* attaching one leaves the toy in
place.

| set | reference | role in LiPolGen | status |
|---|---|---|---|
| **CT18NLO** | T.-J. Hou *et al.*, *New CTEQ global analysis … with high-precision data from the LHC*, **PRD 103 (2021) 014013**, `arXiv:1912.10053` — **VERIFIED** (INSPIRE 1773096) | `LhapdfSF` (`src/lhapdf/lhapdf_sf.cpp`); `--unpol-sf ct18nlo`, `--b1-unpol ct18nlo` | **runnable-code**, optional tier: needs `-DLIPOLGEN_WITH_LHAPDF=ON` and the CT18NLO set in LHAPDF's store |
| **MSTW2008 LO** | A. D. Martin, W. J. Stirling, R. S. Thorne, G. Watt, *Parton distributions for the LHC*, **EPJC 63 (2009) 189**, `arXiv:0901.0002` — **VERIFIED** (INSPIRE 810127) | `MstwSF`, reading PYTHIA 8's own `pdfdata/mstw2008lo.00.dat`; **this is what CDKS computed their b₁ᵈ with** (`pipeline.cpp:1156`) | **in-tree via PYTHIA**, optional tier |
| **NNPDFpol1.1** | E. R. Nocera, R. D. Ball, S. Forte, G. Ridolfi, J. Rojo, *A first unbiased global determination of polarized PDFs and their uncertainties*, **NPB 887 (2014) 276**, `arXiv:1406.5539` — **VERIFIED** (INSPIRE 1302398) | `LhapdfG1`; `--pol-sf nnpdfpol` | **runnable-code**, optional tier |

**The honest caveat the tree already records.** `pipeline.cpp:1323–1347` states
that `sigma_per_category_pb` is **bit-identical** between `--pol-sf toy` and
`--pol-sf nnpdfpol` — the polarized set does not reach the cross-section
normalisation — and `pipeline.cpp:1062` records a run "built on CT18NLO" that in
fact ran at `F2A(0.3, 10) = 0.3691149345`, the **toy** value, with CT18NLO's
0.4639703442 appearing nowhere. So "we use CT18NLO" is a claim about
configuration, not automatically about the numbers in a given run. Any external
comparison in §3 must first assert *which* backend was live.

**What §7 validates:** that the free-nucleon F₂ᵖ, F₂ⁿ and g₁ feeding every
convolution are a globally fitted set with published χ² against §2–§4, rather than
`ToyF2`'s by-eye anchors. **What it does not validate:** anything nuclear,
anything tensor, and the *interface* — an LHAPDF set can be attached and still not
reach the number you are quoting (above).

---

## 8. Nucleon elastic form factors — the other stand-in in the RC tail

`nucleon_ff` (`rc.hpp:665–680`) is dipole G_E^p, G_M^p/μ_p, G_M^n/μ_n plus
**Galster** G_E^n, and its own docstring calls it "a 5 % dipole … adequate for
BOTH" consumers (the isoscalar folding inside `HoSpin1FF`, and POLRAD Eq. (44)'s
quasi-elastic tail carried with a ±100 % band).

The world fits that would replace it, both **VERIFIED**:

* **Z. Ye, J. Arrington, R. J. Hill, G. Lee, *Proton and neutron electromagnetic
  form factors and uncertainties*, PLB 777 (2018) 8, `arXiv:1707.09063`**
  (INSPIRE 1613527) — a world fit **with uncertainties**, which is what a band
  needs.
* **J. J. Kelly, *Simple parametrization of nucleon form factors*, PRC 70 (2004)
  068202** (INSPIRE 669316) — the standard closed form.

**Validates:** the 5 % dipole claim, quantitatively, and it replaces Galster with
a fit that has an error band. **Does NOT validate:** the *nuclear* form factors
(§5.3) or the ±100 % `qe_suppression` band, which is a policy choice about the
quasi-elastic ridge, not a form-factor question. **Obtainability:
table-in-paper / runnable-code** (both parametrizations are closed forms; Ye *et
al.* also ship code). **Priority: LOW** — the docstring's own argument that 5 % is
adequate is sound, and the band dominates.

---

## 9. Kinematic coverage against the EIC window (x ≥ 0.001, Q² ≥ 1 GeV²)

| region | who covers it | for which observable |
|---|---|---|
| x < 10⁻³ | **nobody**, for any deuteron observable | — |
| 10⁻³ ≲ x ≲ 4×10⁻³, Q² ≥ 1 | **nobody** for g₁ᵈ (COMPASS D-5 is there but at Q² < 1); **H1 F_L** (R-5) for R, on the *proton* | R only, proton only |
| 4×10⁻³ ≲ x ≲ 10⁻², Q² ≥ 1 | COMPASS D-2 (g₁ᵈ), NMC F-1 (F₂ᵈ, R) | vector + unpolarized only |
| 10⁻² ≲ x ≲ 0.5, Q² 1–40 | dense: D-1, D-2, D-6, D-7, F-1…F-7, R-1…R-4 | vector + unpolarized only |
| **the whole window, tensor** | **HERMES b₁, five usable bins, 0.032 ≤ x ≤ 0.452** | that is the entirety of it |
| x > 0.5 | BCDMS F-4, MARATHON F-8, E06-014 | unpolarized + twist-3 |

The one-sentence version: **the tensor sector has five measured points in the
entire EIC window, all at A = 2, none below x = 0.03; every other row of this
survey constrains a factor that multiplies b₁ rather than b₁ itself.**

---

## 10. Priority — what to wire, in order

1. **HERMES b₁ᵈ Table II as an assertion (§1.1).** Six rows, already transcribed,
   zero fetching required, and it converts the only measurement of the target
   observable from a comment into a gate. It also forces the
   `B1_MILLER_TABLE_TO_PER_NUCLEON` factor-of-2 decision to be made against data.
2. **HERMES F₂ᵈ, JHEP 05 (2011) 126, HEPData `10.17182/hepdata.66147` (F-7).**
   Same experiment, same kinematics; makes the "HERMES's b₁ is per nucleon"
   inversion checkable rather than argued.
3. **R at EIC x: H1 F_L, HEPData `10.17182/hepdata.45340` (R-5), plus Whitlow
   R1990 `…29555` and Dasu E140 `…22468` (R-1, R-2).** Puts a measured number on
   the `r_sigma_lt`-vs-`r1998` bias that today has one spot value, and measures
   what `r1998`'s clip costs below x = 0.005.
4. **COMPASS final g₁ᵈ, HEPData `10.17182/hepdata.78374` (D-2).** The best
   vector-sector overlap with the EIC window; one clean external gate on
   `g1_nucleus` at A = 2.
5. **Boeglin d(e,e′p)n (§5.4).** Turns `T25 Hulthén against AV18` from
   model-vs-model into model-vs-data out to 550 MeV/c.
6. **Abbott deuteron form factors (§5.3).** Removes the acknowledged HO stand-in
   in the RC elastic tail's deuteron control channel.
7. Everything else: documentation fixes worth more than the wiring — attribute
   η_d = 0.0256(4) to Rodning–Knutson (§5.2), and say in `triton_sf.cpp` that
   Ciofi–Simula A = 3 is a fit to a many-body **calculation** and that a ³He
   column is standing in for ³H (§6).

---

## 11. Verified absences — state these, do not soften them

* **No HEPData record for HERMES b₁** (19 HERMES submissions in INSPIRE's index;
  the b₁ paper is not one). Table-in-paper only.
* **No HEPData record found for**: Abe *et al.* R1998 (R-4), MARATHON (F-8),
  BONuS (F-9), Boeglin (§5.4), DeGrush/BLAST (§1.3), Rvachev ³He(e,e′p) (§6),
  E06-014 (§6), CLAS EG4 (D-10) — subject to the §0 qualifier.
* **No tensor-polarized DIS measurement other than HERMES 2005.** E143 and E155
  measured g₁ᵈ, g₂ᵈ, A₁ᵈ, A₂ᵈ, A_∥ and A_⊥ — **vector and transverse-vector
  asymmetries only**. Neither published a tensor asymmetry or a b₁ limit. Nothing
  before 2005 did either: the pre-HERMES tensor-polarized deuteron electron
  scattering literature is quasielastic and photodisintegration.
* **No second b₁ measurement exists.** JLab E12-13-011 is fully approved (PAC-40;
  target R&D review Aug 2022) and has **not run**; the CLAS12 RG-C tensor analysis
  (`2502.20044`) and the Hall D frozen-spin proposal (`2504.21177`) are proposals.
  Cosyn *et al.* EPJ A 61 (2025) 83 call the JLab experiment "the **second**
  measurement", in the future tense.
* **No R measurement on any nucleus at EIC kinematics, and no tensor L/T
  separation ever.** `cosyn_tensor_sfs`'s split of the tensor response into T and
  L pieces is unconstrained by data.
* **No polarized deuteron spectator-tagged DIS.** BONuS is unpolarized; the
  Cosyn–Weiss tagged-tensor framework (`arXiv:2603.23699`, `2603.23700`, already
  the tree's strongest external gate) is theory with no data behind it.
* **No ⁶Li or ⁷Li polarized DIS, tagged or inclusive, at any x or Q².**

---

*Sources for every claim above were fetched in this session: arXiv abstract pages
via the arXiv API; journal reference, collaboration and record identity via the
INSPIRE literature API; HEPData record existence, DOI, table count and observable
keywords via INSPIRE's `data` collection with the `literature` back-link resolved
to its paper. `hepdata.net` itself was unreachable from this environment
(Cloudflare challenge, HTTP 403 on both `curl` and `WebFetch`), so no table
contents were read. Nothing was downloaded into either repository. The two live
test invocations quoted in §5.1 and §1 were run read-only against the existing
`build/`: `-tc="R1998 is the average of its three published forms, …"` (3 cases,
1173 assertions) and `-tc="cluster: the CD-Bonn deuteron coefficients (Machleidt
Table XX), tagged: the Cosyn-Weiss deuteron tensor gate (CW TABLE II),
b1_nuclear T1 gate checklist item 5v…, T8(a): the DEUTERON gate…"` (4 cases,
54 633 assertions) — 7 cases, 0 failures.*
