<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-26 — phase C2: the one bound on the Δ sector (a moment-level gluonic Soffer-type bound, under stated assumptions)

Worker C2, 2026-10-26, in two passes. Plan: `PLAN.md` phase C2 (rewritten
after the phase-0r review); ranking: `../../benchmarking/08_channel_priority.md`
§4 ingredient C, under rank 1 (inclusive ⁶Li, cos 2φ — the letter's claim 1).
New row: `validation/benchmarks/t5_delta_moment_bound.py`, **PASS**. Nothing in
the tree moved: no physics default, constant, registry label, tolerance or
reference-JSON number changed; no decision was taken.

The **first pass** wrote the harness, its vendored sources file, its tests and
this record's first version, with every assumption and the tolerance in the
header before any saturation was computed. The **second pass** (the same day,
after every first-pass number had been computed and seen) re-read the sources
and found two things, both recorded in the harness header under SECOND PASS
and in §4 below: the first pass's variation V4 (a "second normalisation" of
the parton-model relation) is a misreading, and the first pass's assumption A7
(the kernel's Δ is the Δ the bound is written for) disagrees with Detmold and
Shanahan's own definition of Δ by **a factor 2**. The second pass changed no
central assumption, no central number, not the verdict rule and not its
tolerance; it added a source (KS20), the variation V6 and the fields that
carry it, and re-defined the two recorded band edges.

> **A bound says the scenario is ALLOWED, never that it is right**
> (`../../benchmarking/04_theory.md` §11 item 2). The PASS below says the
> shipped toy Δ does not ask for more linearly polarized glue than positivity
> permits under the stated assumptions. It says nothing about whether exotic
> glue in ⁶Li exists, or has the toy's size, shape or sign; and a scenario of
> the letter being *allowed* (or not) under a reading of the normalisation is
> not evidence for (or against) the physics it describes.

## 1. What was compared

The n = 2 Mellin moment of the generator's double-helicity-flip structure
function, M₂(Q²) = ∫₀¹ x Δ(x, Q²) dx (Detmold–Shanahan [DS16] Eq. (12) at
n = 2; the plain integral of the xΔ the letter extracts), at Q² = 2, 5,
10 GeV², on ⁶Li **per nucleon in the per-nucleon x**, against

  B(Q²) = (α_s/8π) Tr Q² · ½ ∫₀¹ x [f₁^g + ½ f₁LL^g + g₁^g] dx ,

with the saturation fraction S = |M₂|/B (how much of the allowed room the
scenario uses) and the margin 1 − S. B is the n = 2 moment of DS16's x-space
bound Eq. (32), |δG(x)| ≤ ½(f₁ + ½f₁LL + g₁), carried onto Δ by DS16's
restatement of the Jaffe–Manohar relation, Eq. (13), whose n = 2 moment is
M₂ = −(α_s/8π) Tr Q² ∫ x δG dx (the same as DS16 Eqs. (9) + (11) + (12) with
A₂ = ∫ x δG dx). DS16's own Eq. (36), |A₂| ≤ (5B₂,₁ − 6B₂,₂)/24, is that
statement written for the φ meson in lattice reduced matrix elements (the
map from B₂,₁, B₂,₂ to the moments of f₁, f₁LL is not printed), so it is
quoted, not evaluated.

**Which Δ, in which normalisation (read from the kernel).** The kernel's
cos 2φ term is a₂ += −(1−y)/y² · c_eff sin²θ_S · Δ/D_φ with the `qc.second`
moment (`src/core/xsec.cpp:281`), c_eff = 3Q_NN = 3m² − 2 at J = 1
(`src/core/xsec.cpp:178-179`, `tensor_moments`); Δ is evaluated on the
per-nucleon F₁ (`t.f1 = nf2_.f1a(x, q2) / static_cast<double>(ion_.A)`,
`src/core/xsec.cpp:95`; `t.delta = (*df) ? (*df)(x, q2, t.f1) : 0.0`,
`src/core/xsec.cpp:117`). The pipeline fills the slot with
`toy_delta_gluon(x, q2, f1, 1e-2)` (`src/core/pipeline.cpp:1264-1265`; its run
banner says "toy_delta_gluon at the 1e-2 discovery scale",
`src/core/pipeline.cpp:1292`), not the function's default 1e-3
(`include/lipolgen/sf.hpp:430`); `toy_delta_gluon` is
`scale * f1 * std::pow(x, 0.3) * std::pow(1.0 - x, 4.0)` (`src/core/sf.cpp:377-378`).
**The tree applies no dilution to Δ** (P_zz is the whole-nucleus fill; nothing
multiplies the slot). The row takes Δ per nucleon exactly as the kernel
carries it; the harness reads it from `default_inclusive_kernel(li6())`'s own
tables and measures the scale: 0.01 at every quadrature node (maximum relative
deviation 4.4e−16). The standalone identity `a_cos2phi` returns the same
coefficient (`src/core/asymmetries.cpp:90-92`), and the pytest checks that the
kernel's a₂ at m = +1, θ_S = 90° equals it to 1e−14 at three (x, Q², y).
Which *cross-section* normalisation that Δ carries is §4.2.

**The models.**
1. **Shipped (the verdict):** the toy above.
2. **The letter's** `moment_A` and `moment_B` (PolarizedLithiumSim `plans/07`
   §7.1, "moment_A, P_zz = 0.6, dilution 1/3"; the sibling's
   `polli_fastsim/delta_models.py`, read, never run) through the tree's port
   `make_moment_a` / `make_moment_b`: variant mid_x, c = `C_BAG` = −0.012
   (`include/lipolgen/constants.hpp:141`, Sather–Schmidt's bag value for the
   Δ⁺⁺), α_s the CT18NLO table, dilution 1/3 (the counting fraction the
   sibling folds into P_zz = 0.8/3 and the tree expresses on Δ with a
   whole-nucleus P_zz — `sf.hpp`'s `dilution` comment); `moment_A` solved on
   the shipped kernel's F₁ at q²_ref = 4.474 GeV² (the letter's toy-column
   ⟨Q²⟩). **Reproduced:** A = −0.292399 against the letter's printed −0.2924
   (`plans/07` addendum table, toy column).
3. **Recorded:** both at dilution 1 (the library default; the bag value per
   nucleon, undiluted).

## 2. The sources, read on disk (no PDF in the repository)

| key | document | file (sha256) | what this row takes from it |
|---|---|---|---|
| DS16 | Detmold & Shanahan, PRD 94 (2016) 014507, arXiv:1606.04505v3 | `PolarizedLithiumSim/refs/1606.04505.pdf` (`6cea5915…`) | the **definition** of Δ, Eqs. (3)–(5) (p. 2, page image: Δ = A₊₋,₋₊ with the projection operators printed above Eq. (3), and the double-helicity-flip tensor with its ½); Eqs. (9), (11)–(14) (pp. 2–3) the relation; Eq. (32) (p. 6) the bound, "δG(x) is the gluonic transversity distribution defined in Eq. (14)", named h₁TT in its Refs. [28, 29]; Eq. (36) (p. 7) its φ-meson moment form |
| SS90 | Sather & Schmidt, PRD 42 (1990) 1424 | session scratch, APS, fetched 2026-10-26 (`63378c24…`) | Eq. (3.5): ∫x Δ dx = −0.012 α_s for the Δ⁺⁺ with three flavours (C_BAG and Tr Q² = 2/3); Eq. (2.11), written (p. 1425) for a target "aligned (i.e., circularly polarized)" along x; Eq. (2.5) |
| LOI18 | Maxwell *et al.*, JLab PAC 44 LOI, arXiv:1803.11206 | `PolarizedLithiumSim/refs/1803.11206.pdf` (`30d2bd48…`) | Eq. (6) (p. 4 image): the Bjorken-limit cross section for a target with its spin along x̂, Δ term −x(1−y)/2 · Δ cos 2φ; Eq. (7) the relation with Q = diag(2/3, −1/3, −1/3); Eq. (9) (p. 5 image): the per-state cross section, the tree's a₂ term in form |
| CvDM17 | Cotogno, van Daal, Mulders, JHEP 11 (2017) 185, arXiv:1709.07827 | `PolarizedLithiumSim/refs/1709.07827.pdf` (`7ff97995…`) | Eqs. (34)–(36) (p. 6 image): \|g₁\| ≤ f₁ + f₁LL/2, f₁LL ≤ f₁, \|h₁TT\| ≤ ½(f₁ + f₁LL/2 + g₁) — DS16's (32) as a positivity bound |
| CRSZ25 | Cosyn, Roldan Tomei, Sosa, Zec, EPJ A 61 (2025) 83, arXiv:2410.12764 | `PolarizedLithiumSim/refs/2410.12764v1.pdf` (`f05fc42b…`) | Eq. (2) (p. 1): A^T = [dσ(+1) + dσ(−1) − 2dσ(0)]/[dσ(+1) + dσ(0) + dσ(−1)], printed with **[−2 ≤ A^T ≤ 1]** — the helicity-state positivity of the spin-1 cross sections |
| KS20 (second pass) | Kumano & Song, PRD 101 (2020) 054011, arXiv:1910.12523v3 | `PolarizedLithiumSim/refs/1910.12523.pdf` (`f13001ae…`) | Eq. (68) (p. 11 image): the dictionary of the gluon-transversity notations — DS16's δG [its ref. 28] and SS90 [its ref. 48] on the line of h₁TT [39, 41, 49], −g_x̂/x̂ + g_ŷ/x̂, the hadron polarization x̂ being the **linear** E_x of Eq. (22) (p. 6 image); Eq. (67) the gluon correlator |

- **The Jaffe–Manohar relation is sourced on disk**: DS16 Eq. (13) restates it
  (Jaffe–Manohar PLB 223 and HJM NPB 312 are not on disk); SS90 Eq. (2.11) and
  LOI18 Eq. (7) restate it too. The first pass read SS90's coefficient as twice
  DS16's; the second pass shows the two are one relation (§4.1).
- **The spin-1 helicity positivity is sourced** (CRSZ25 Eq. (2)): it went in
  as the second sub-row (§6). No b₁-specific positivity without that source
  was added.
- Every equation the harness uses was read off a page image; all are recorded
  verbatim (ASCII) in `validation/benchmarks/data/delta_moment_bound_sources.json`
  with the file sha256, the pages read and the role.

## 3. The assumptions and the tolerance (first pass, declared before any number)

A1 the relation and bound of DS16 (above). A2 Tr Q² = 2/3 (SS90 above
Eq. (3.5), LOI18 Eq. (7); DS16 gives none). A3 α_s(Q²) of CT18NLO member 0
(LHAPDF's `alphasQ2`), also handed to the letter's models. A4 f₁^g = the
CT18NLO gluon per nucleon, free nucleon, read from LHAPDF 6.5.5 itself through
PYTHIA 8.317's LHAPDF6 plugin (`_report.import_pythia8`, then `PDF::xfUpdate` +
`PDF::xfRaw` — PYTHIA's `PDF::xf` clips the gluon at zero); the tree's
`LhapdfSF`/`LhapdfG1` expose F₂ and g₁ only and the `lhapdf` module is not
built, so the brief's MSTW fallback was not needed; a same-grid guard holds the
plugin to `LhapdfSF('CT18NLO').f2p` and `LhapdfG1('NNPDFpol11_100').g1p` at six
points each (worst relative difference 0.0). A5 g₁^g = P_eff · Δg, Δg from
NNPDFpol11_100 member 0 per nucleon, P_eff = (Z P_p + N P_n)/A of `li6()` =
0.2704095 (the 6Li m = +1 state). A6 f₁LL^g = 0; band ±f₁^g, the largest band
symmetric about zero CvDM17 (34)–(35) allow. A7 the tree's Δ is DS16's (via
LOI18 Eq. (9), the kernel's a₂ with c_eff = 3m² − 2). A8 per nucleon on both
sides: with x = A x_A and per-nucleon densities D(x) = D_A(x_A)/A², Eq. (13)
and Eq. (32) are unchanged — asserted numerically by the pytest at A = 6.

**The tolerance** (rule (a) fixed by the brief, rule (b) the first pass's own
with the sourced sub-row the brief asked for, both before any number;
unchanged by the second pass): the row PASSES iff the shipped Δ has margin 1 − |M₂|/B ≥ −1e−12
at each of Q² = 2, 5, 10 GeV² under A1–A8, **and** every helicity-state rate
of the shipped ⁶Li kernel is ≥ −1e−12 relative to the unpolarized rate
(−2 ≤ A_T ≤ 1) on the vendored grid. The letter's models get the outcome of the
first rule, recorded beside the verdict (since the second pass under both A7
and V6). **Preview, as the rules require:** a back-of-envelope estimate
(⟨x⟩_g ≈ 0.42) preceded the header — toy ~0.1, bag models ~0.7 at dilution 1/3
and ~2 at dilution 1 — and is stated in the header. One numerics change was
made after the first run: the quadrature (72 + 36 panels) failed the harness's
own convergence check (1.7e−7 > 1e−7 on ⟨x⟩_g at Q² = 2) before any saturation
was computed, and the panels were doubled (144 + 72; plateau ~1e−9 relative).

## 4. The second pass: what was found after the numbers, and what it changed

### 4.1 S1 — V4 is a misreading, not a second normalisation

SS90 write Eq. (2.11), Δ = (α_s/2π) Tr Q² x² ∫(dy/y³) a(y) with
a = −2(g_x̂ − g_ŷ), for a spin-1 target "aligned (i.e., circularly polarized)"
along x (SS90 p. 1425). DS16's δG = g_x̂ − g_ŷ (Eq. (14), "a target ...
polarized in the x̂ direction") is read by KS20 as the linear-polarization
difference in the **linearly** polarized hadron: KS20 Eq. (68) puts DS16's δG
and SS90's function on the line of h₁TT, −g_x̂/x̂ + g_ŷ/x̂, the hadron
polarization x̂ being the linear E_x of its Eq. (22) — which is also the only
state in which DS16's own identification δG = h₁TT holds in magnitude. For a spin-1 target the transverse tensor polarization T_xx − T_yy
is **+½** in the circular (m_x = ±1) state and **−1** in the linear (E_x,
m_x = 0) state (the pytest computes both from the spin-1 matrices), so
(g_x − g_y)_circular = −½ (g_x − g_y)_linear, SS90's a = (g_x − g_y)_linear,
and SS90 (2.11) **is** DS16 (13). The first pass's V4 (SS90's coefficient
applied to DS16's δG, ×2 room) mixes the two states. It stays a recorded
number (`saturation_ss90`, labelled a misreading) and leaves the LOOSE edge.

### 4.2 S2 — the kernel's Δ is half of DS16's Δ at equal cos 2φ amplitude (V6)

DS16 *defines* Δ by Eqs. (3)–(4): W_μν,αβ = Σ P(hH, h′H′) A_hH,h′H′ with
P(hH, h′H′)_μν,αβ = ε*_μ(h′) E*_α(H′) E_β(H) ε_ν(h), and Δ ≡ A₊₋,₋₊ = A₋₊,₊₋.
Three checks, each in the pytest (numpy only):

1. Eqs. (3)–(4) give Eq. (5)'s tensor, with its printed ½: for E = x̂ the
   transverse block is W_xx = −W_yy = Δ/2 (and −Δ/4, +Δ/4 for m_x = +1); the
   same operators with A_hH,hH = 1 give, for an unpolarized target,
   W_ij = δ_ij · ⅓Σ_H A₊H,₊H, i.e. the standard −F₁ g_μν with
   F₁ = ⅓Σ_H A₊H,₊H. So DS16's Δ and F₁ are normalised together.
2. Eq. (5) and −F₁ g_μν contracted with the unpolarized lepton tensor give,
   for a pure m = ±1 state along an axis perpendicular to q, a cos 2φ term
   **−½ (1−y)/y² · Δ** relative to F₁ + (1−y)F₂/(xy²) (the Bjorken limit,
   approached as 1/Q²: this coefficient over the kernel's is 0.4961 at
   x = 0.1, y = 0.5, Q² = 5 GeV² and within 2e−5 of ½ at Q² = 1e5 GeV²; the
   m_x = 0 state carries −2 times the m_x = ±1 one, as c_eff does). That is
   **LOI18 Eq. (6)**, the Bjorken limit the LOI prints for a target with its
   spin along x̂: [xy²F₁ + (1−y)F₂ − (x(1−y)/2) Δ cos 2φ] × e⁴ME/(4π²Q⁴).
3. The kernel — and LOI18 Eq. (9), which A7 read — carries **−(1−y)/y² · Δ**
   for the same state (c_eff = 1, θ_S = 90°): the LOI's own Eqs. (6) and (9)
   differ by a factor 2 in the Δ term.

So, per unit cos 2φ amplitude, **Δ_DS16 = 2 Δ_tree**, and the bound — written
for DS16's Δ — reads on the tree's Δ as |2M₂| ≤ B. The harness carries this as
**V6** (`DS16_OVER_TREE = 2`, a derived ratio re-derived by the pytest, never
set from a result): every saturation doubles. The sibling's own code review
of 2026-08-25 had checked the y-structure and the m-dependence of the same term
and left its overall normalisation open for want of Jaffe–Manohar 1989
(PolarizedLithiumSim `docs/code_review_2026-08-25.md` §3, the master-formula
row); this pass settles it against DS16's printed definition, not against JM.
**A7 stays the pre-registered central reading** (the rules forbid moving a
central assumption after its numbers were seen); the on-disk definitions
support V6, and §7 reads the result under both.

### 4.3 S3 — the edges

TIGHT = V1 (Δg = 0) + V2(−) (f₁LL^g = −f₁^g) + V6, on the smaller of the free
and EPPS21 gluon moments; LOOSE = V2(+) (f₁LL^g = +f₁^g) + V3 (Tr Q² = 10/9) on
the larger, with the NNPDF Δg term only if it widens B. First-pass edges, for
the record: TIGHT without V6 (toy 0.136 / 0.168 / 0.197; moment_A, letter
1.508 / 1.466 / 1.469), LOOSE with V4 (toy 0.0130 / 0.0162 / 0.0190; moment_A,
letter 0.144 / 0.141 / 0.142; B_loose 9.521196e−3 / 8.097994e−3 /
7.280604e−3).

## 5. The numbers (one harness run, second pass; windows as stated)

### 5.1 The bound side (per nucleon, central assumptions)

| Q² [GeV²] | α_s (CT18NLO) | ⟨x⟩_g CT18NLO (EPPS21 ⁶Li) | ⟨x⟩_Δg NNPDFpol11 | P_eff ⟨x⟩_Δg | **B** | B, Δg = 0 | B, f₁LL = −f₁ / +f₁ | B, Tr Q² = 10/9 | B, P_eff = ⅓ | tight / loose B |
|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 0.354535 | 0.39760193 (0.38739898) | 0.04087300 | 0.01105245 | **1.921563e−3** | 1.869592e−3 | 9.867666e−4 / 2.856359e−3 | 3.202604e−3 | 1.933656e−3 | 9.108081e−4 / 4.760598e−3 |
| 5 | 0.283893 | 0.42296138 (0.41331810) | 0.03984385 | 0.01077415 | **1.633121e−3** | 1.592554e−3 | 8.368443e−4 / 2.429398e−3 | 2.721869e−3 | 1.642561e−3 | 7.781223e−4 / 4.048997e−3 |
| 10 | 0.247420 | 0.43670379 (0.42745569) | 0.03900526 | 0.01054739 | **1.467658e−3** | 1.433047e−3 | 7.511347e−4 / 2.184181e−3 | 2.446097e−3 | 1.475712e−3 | 7.013495e−4 / 3.640302e−3 |

(The TIGHT edge's V6 acts on the model side: S_tight = 2|M₂|/B_tight.)
Quadrature errors (fine − coarse): ⟨x⟩_g ≤ 1.2e−9, ⟨x⟩_Δg ≤ 8.4e−11, EPPS21
⟨x⟩_g ≤ 1.2e−8; every model M₂ ≤ 8.7e−19 absolute. Tails below x_lo:
x_lo|f(x_lo)| ≤ 1.1e−7 (gluon, x_lo = 1e−9), ≤ 1.3e−6 (Δg, below its grid's
x_min = 1e−5), ≤ 3.0e−13 (toy). **The mixed central inputs leave positivity at
large x:** |g₁^g| > f₁^g/2 from x = 0.9494 / 0.9431 / 0.9412 and |g₁^g| > f₁^g
(CvDM17 (34) itself) from x = 0.9713 / 0.9694 / 0.9681 at Q² = 2 / 5 / 10 — a
region carrying 2.1e−7 / 9.7e−8 / 5.2e−8 of ⟨x⟩_g; recorded, not repaired.

### 5.2 Saturation |M₂|/B per model

| model | M₂ at Q² = 2 / 5 / 10 | **S central (A7)** | **S under V6** | S tight (with V6) | S loose | S = 1 at f₁LL^g/f₁^g = (central; V6) | rule (a) |
|---|---|---|---|---|---|---|---|
| **toy @ 1e−2 (shipped)** | 1.2383e−4 / 1.3108e−4 / 1.3805e−4 | **0.0644 / 0.0803 / 0.0941** | **0.129 / 0.161 / 0.188** | 0.272 / 0.337 / 0.394 | 0.0260 / 0.0324 / 0.0379 | −1.92 / −1.89 / −1.86; −1.79 / −1.72 / −1.66 | **PASS (verdict)**; PASS under V6 |
| moment_A, letter (dilution ⅓) | −1.3733e−3 / −1.1411e−3 / −1.0303e−3 | 0.715 / 0.699 / 0.702 | **1.429 / 1.397 / 1.404** | 3.016 / 2.933 / 2.938 | 0.288 / 0.282 / 0.283 | −0.586 / −0.618 / −0.610; **+0.883 / +0.815 / +0.828** | recorded: PASS centrally, **FAIL under V6** |
| moment_B, letter (dilution ⅓) | −1.4181e−3 / −1.1356e−3 / −9.8968e−4 | 0.738 / 0.695 / 0.674 | **1.476 / 1.391 / 1.349** | 3.114 / 2.919 / 2.822 | 0.298 / 0.280 / 0.272 | −0.539 / −0.625 / −0.667; **+0.979 / +0.801 / +0.714** | recorded: PASS centrally, **FAIL under V6** |
| moment_A, dilution 1 | −4.1199e−3 / −3.4233e−3 / −3.0909e−3 | 2.144 / 2.096 / 2.106 | 4.288 / 4.192 / 4.212 | 9.05 / 8.80 / 8.81 | 0.865 / 0.845 / 0.849 | +2.35 / +2.25 / +2.27; +6.76 / +6.55 / +6.58 | recorded: FAIL / FAIL |
| moment_B, dilution 1 | −4.2544e−3 / −3.4067e−3 / −2.9690e−3 | 2.214 / 2.086 / 2.023 | 4.428 / 4.172 / 4.046 | 9.34 / 8.76 / 8.47 | 0.894 / 0.841 / 0.816 | +2.50 / +2.23 / +2.10; +7.05 / +6.51 / +6.24 | recorded: FAIL / FAIL |

Single variations on the central reading (S at Q² = 2 / 5 / 10), for the
letter's moment_A: Δg = 0 0.735 / 0.717 / 0.719; f₁LL = −f₁ 1.392 / 1.364 /
1.372; f₁LL = +f₁ 0.481 / 0.470 / 0.472; Tr Q² = 10/9 0.429 / 0.419 / 0.421;
EPPS21 0.733 / 0.715 / 0.717; P_eff = ⅓ 0.710 / 0.695 / 0.698; with the
library-default LO α_s on the model side 0.817 / 0.800 / 0.803; V4 (the
misreading) 0.357 / 0.349 / 0.351. Under V6 each doubles.

Cross-checks: `moment_B`'s M₂ equals dilution · C_BAG · α_s(Q²) to 4.4e−16
(closed form, asserted by the harness); the toy's M₂ recomputed by scipy's
adaptive `quad` in ln x on the kernel's own Δ agrees to ≤ 1e−8 relative (the
pytest `test_toy_moment_by_an_independent_quadrature`); the same-grid guard
0.0 on f₂ᵖ and g₁ᵖ.

## 6. The second, sourced sub-row: helicity-state positivity

The shipped ⁶Li kernel's per-state rates W_m(φ′) = 1 + w_avg + a₁ cos φ′ +
a₂ cos 2φ′ (unpolarized electrons, pure m = +1, 0, −1 along an axis at
θ_S = 0°, 45°, 90°; x ∈ {0.001, 0.01, 0.05, 0.1, 0.3, 0.6, 0.9}, Q² ∈ {2, 5,
10}, y ∈ {0.05, 0.5, 0.95}, φ′ ∈ {0°, 45°, 90°}: 567 points, 1701 rates):
**min W_m = 0.9494** (at x = 0.9, Q² = 2, y = 0.95, θ_S = 90°, φ′ = 90°,
m = 0), **A_T ∈ [−0.1012, 0.0506]** against [−2, 1] — **PASS**. The rates carry
b₁, b₂ and Δ; the margin is set by b₁/F₁ ≈ 0.15 at x = 0.9, not by Δ. Under V6
nothing here changes: the sub-row reads the kernel's own rates.

## 7. What the verdict means for the channel (rank 1, claim 1; ingredient C)

- **The shipped Δ is allowed under either normalisation**: the toy at the
  pipeline's 1e−2 uses 6.4–9.4 % of the room under A7 and 12.9–18.8 % under V6
  (27–39 % at the tight edge, which includes V6). The toy is an unnormalised
  ansatz, so V6 only re-labels what "1e−2" means in DS16's units (2e−2).
- **The letter's scenario is allowed under A7 and beyond the bound under V6.**
  moment_A at dilution ⅓ uses 0.70–0.71 of the room under A7 and **1.40–1.43**
  under V6 (moment_B 0.67–0.74 and **1.35–1.48**). Under V6 it is allowed only
  if the ⁶Li gluons are tensor-polarized at f₁LL^g ≥ +0.82 … +0.88 f₁^g
  (moment_B +0.71 … +0.98) — close to the CvDM17 (35) maximum f₁LL ≤ f₁ — or
  at the loose edge (f₁LL^g = +f₁^g with Tr Q² = 10/9: 0.56–0.58 under V6).
  Nothing measures or calculates f₁LL^g for any nucleus. The largest dilution
  the central bound admits, from these saturations (S is linear in the
  dilution): 0.466 (moment_A) / 0.452 (moment_B) under A7, **0.233 / 0.226**
  under V6.
- **What V6 means for claim 1, if it holds.** The kernel's cos 2φ amplitude is
  linear in Δ, and per unit of DS16's Δ — Jaffe–Manohar's ("We follow that
  reference in defining … Δ", DS16 Sec. II), the Δ SS90's bag value and DS16's
  bound are written for — it is twice what DS16's Eqs. (3)–(5) give. Then
  either (i) the letter's scenario is SS90's bag value in those units: the
  tree should carry C_BAG/2 in its own Δ, the scenario uses 0.67–0.74 of the
  room (the A7 column), and every cos 2φ amplitude the letter predicts for it
  is **twice** what that scenario gives; or (ii) the letter's amplitudes stand
  as computed: the scenario they encode is twice SS90's bag value at the stated
  dilution and uses 1.35–1.48 of the room under the central assumptions (the
  V6 column). Which of the two the letter means is the author's; this row
  records both and decides neither.
- **The dilution ⅓ is load-bearing for admissibility.** Undiluted (C_BAG · α_s,
  the library default `dilution = 1`) the bag value exceeds the bound ×2.0–2.2
  under A7 and ×4.0–4.4 under V6; it would need f₁LL^g ≥ +2.1 … +2.5 f₁^g
  (A7), which CvDM17 (35) forbids. The open question whether Δ should carry
  the rank-2 transfer as b₁ does (⅓ → 0.3073, `include/lipolgen/sf.hpp:541-543`,
  `LI6_B1_RANK2_TRANSFER`) would scale every letter saturation by 0.921947:
  0.659 / 0.644 / 0.647 (A7) and 1.318 / 1.288 / 1.294 (V6) for moment_A.
- Context only, not a comparison (04 §4 forbids comparing the lattice φ number
  to a nuclear Δ): DS16 find the φ meson's moment bound saturated at 80–100 %
  (bare lattice operators, m_π = 450 MeV).
- For `08_channel_priority.md` §4 row C: "no gate at all" becomes "one
  moment-level bound under stated assumptions: PASS for the shipped toy (6–9 %
  of the room; 13–19 % under DS16's normalisation of Δ); the letter's scenario
  0.70 under the first-pass normalisation and 1.35–1.48 — beyond the bound —
  under DS16's own definition of Δ, which makes the kernel's cos 2φ
  coefficient twice DS16's per unit Δ (a factor-2 question on claim 1's
  amplitudes, filed for the author)". It does not move rank 1's C3 score (Δ is
  still a scenario with no external calculation).

## 8. What this does NOT validate

The size, shape or sign of Δ(⁶Li); which of A7 and V6 the tree's Δ follows (a
property of the kernel's coefficient, recorded both ways); the x-space bound
Eq. (32) itself (the moment form is weaker: a model can pass the moment and
fail pointwise — any Δ falling more slowly than the gluon as x → 1 does); the
parton-model relation Eq. (13) itself (taken from DS16, SS90 and KS20 as
printed; not re-derived from the photon–gluon box); PDF uncertainties
(CT18NLO's 58 eigenvectors, NNPDFpol1.1's 100 replicas — not propagated; the
declared Δg = 0 edge and the f₁LL band are what the result is read against);
n = 4, 6, … moments; nuclear modification beyond the EPPS21 gluon ratio (V5);
RC on the Δ term (none exists); the finite-γ b₁ leakage into the cos 2φ
amplitude; ⁷Li (its rank-2 sector is identically zero; T-13's spin-3/2
positivity is not this row); b₃/b₄. The helicity-positivity sub-row is a
necessary condition on the generator's per-state rates, not a validation of
b₁(⁶Li).

## 9. Decisions priced (none taken)

- **Registry row 23** (`LI6_CLUSTER_POLARIZATION`, `include/lipolgen/beams.hpp:136-138`,
  against the naive 1/3, `LI6_NAIVE_ONE_THIRD` at `include/lipolgen/beams.hpp:155`):
  on this row it enters only through g₁^g; P_eff = ⅓ instead of 0.2704 widens B
  by +0.63 / +0.58 / +0.55 % at Q² = 2 / 5 / 10 (`B_peff_naive`) and lowers
  every saturation by the same factor. No verdict moves.
- **A registry candidate (new): the cross-section normalisation of Δ in the
  kernel's cos 2φ term.** Options and measured cost: (i) keep the kernel and
  read its Δ as half of DS16's — the letter's bag scenario in DS16's units then
  needs C_BAG/2 in the tree, every letter cos 2φ amplitude at that
  normalisation halves, its saturation is 0.67–0.74; (ii) halve the kernel's
  Δ coefficient to DS16's — the same physics as (i), and the toy's "1e−2"
  produces half its present amplitude; (iii) keep both as they are — the letter's
  amplitudes encode twice the bag value, which uses 1.35–1.48 of the gluonic
  Soffer room under the central assumptions (allowed only for
  f₁LL^g ≳ +0.7 … +1.0 f₁^g). A copy of Jaffe–Manohar PLB 223 (1989) 218
  (`../../references/REFERENCES.md` §1 item 4, "PLEASE DOWNLOAD") would show which
  of the LOI's Eqs. (6) and (9) is theirs.
- Not a registry row (open in `sf.hpp` and `plans/07` §7.1 row 3): the Δ
  dilution ⅓ vs the rank-2 transfer 0.3073 — priced in §7 (×0.921947).

## 10. Files touched, commands, runtimes

Owned by this task (created in the first pass, revised in the second):
- `validation/benchmarks/t5_delta_moment_bound.py` — the harness (`--plot`),
  EXPENSIVE = False, 2.8–3.1 s wall (`run()` 2.76 s).
- `validation/benchmarks/data/delta_moment_bound_sources.json` — the six
  sources (equations verbatim, each PDF's sha256, pages read), the two PDF sets
  (set index, data version, `.info` and member-file sha256, checked against the
  installed files on every run), the letter's model configuration with
  provenance, the grids; file sha256
  `794d1429662a41bf717a7269071d9e3bdafa6c8f4933a92557e38807912218e7`, body
  sha256 `d355287d1cbb0fd06269ed4ee5af58e5b1d8790b8aa314bce112cdf109364594`,
  both re-checked on every read (first pass: `fa3b171d…` / `9c05e166…`).
- `python/tests/test_bench_t5_delta_bound.py` — 21 tests: the file guards (Q²
  point count 3, first 2.0, last 10.0, the positivity grid's 7 x points, one
  flipped byte in the body and one in the provenance refused); the equations;
  the n = 2 moment algebra and the per-nucleon invariance; the second pass's
  S1 (the spin-1 tensor polarizations) and S2 (DS16 Eqs. (3)–(5) against the
  kernel's coefficient, and the kernel against `a_cos2phi`); the headline and
  every sub-row pinned, the letter's models pinned as recorded; the closed
  forms of B and of both edges; the verdict rule on synthetic models; the toy's
  M₂ by an independent quadrature; the BLOCKED path; the REPORT line; exit 0 on
  a run and 2 on an altered file; the PNG.
- `docs/benchmarking/figures/t5_delta_moment_bound.png` — written once by
  `--plot` (57.0 kB): |M₂|/B per Q² and model, filled = central, open = V6,
  whiskers loose → tight, the bound at 1.
- this record.

Commands (from `/home/cpeng/Projects/polli/PoliGen`, after `source env.sh`):

    python3 validation/benchmarks/t5_delta_moment_bound.py --plot          # exit 0, PASS
    python3 -m pytest -p no:cacheprovider python/tests/test_bench_t5_delta_bound.py -q
                                                                           # 21 passed (~10 s)
    python3 validation/benchmarks/make_report.py --only t5_delta_moment_bound \
        --out-dir <session scratch>/report_C2                              # 1 pass
    python3 validation/check_spdx_headers.py                              # 131 files, all carry the line
    python3 validation/check_physics_channels_links.py --records-only     # this record's citations

Outside the repository (session scratch, `c2v/`): the second pass's JSON
generator (`scripts/make_sources_json.py`, which rebuilds the vendored file
from the first pass's and prints both digests), the two derivation checks the
pytest repeats (`scripts/delta_norm_check.py`, `scripts/delta_projector_check.py`),
the page renders read, and a copy of every first-pass file with its full
`run()` output (`first_pass/`).
