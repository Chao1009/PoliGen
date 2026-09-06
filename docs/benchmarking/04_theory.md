<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# 04 — THEORY: published calculations LiPolGen implements, or could be gated against

**Survey date 2026-09-06.** Companion to `00_in_tree_checks.md`. That file
inventories what the tree *already* checks; this one asks a different question of
the same literature: **for each theoretical calculation, is there a NUMBER — a
table, a closed form, a digitizable curve, a runnable code — that LiPolGen could
be gated against, what is that number's own uncertainty, and which LiPolGen
symbol would it constrain?**

Every reference below was verified in this session against the arXiv abstract
page, the INSPIRE record, the LHAPDF set index, the source repository, or (where
stated) the paper's own PDF text. Items I could **not** verify are marked
UNVERIFIED rather than repeated from memory.

---

## 0. What a "theory benchmark" can and cannot be here

Three frame facts, inherited from `00_in_tree_checks.md` §0 and re-checked:

1. **There is no polarized e + ⁶Li/⁷Li DIS calculation, by anyone.** An INSPIRE
   search for a tagged-DIS calculation with a ⁶Li / lithium / α-spectator target
   returns **zero** hits, and the same is true of b₁ for any A > 2. So no
   theory paper computes the observable this generator produces, and none can.
2. Consequently every entry below is a benchmark of a **component**: an A = 2
   tensor calculation, an ab-initio ⁶Li wave-function observable, a Glauber
   profile, an unpolarized nuclear-PDF ratio, a kinematic identity. The
   decomposition is the deliverable, not a headline number.
3. A theory comparison is worth less than a data comparison **and more than
   nothing**, but only when its own pedigree is stated. Three of the strongest
   candidates below (Gakh–Shekhovtsova, Sather–Schmidt, the WS98 ⁶Li quadrupole)
   are single, unchecked, or openly-discrepant calculations, and are labelled so.

### Verdict vocabulary used in the tables

| tag | meaning |
|---|---|
| **TABLE** | numbers printed in the paper (or an online table), transcribable |
| **CLOSED FORM** | an equation one can evaluate independently — the cheapest and strongest kind of gate |
| **FIGURE** | curve only; needs digitization, with the digitizer's error added |
| **CODE** | a released program one can run and diff against |
| **NONE** | the calculation exists but publishes nothing gateable |
| **DOES NOT EXIST** | nobody has done the calculation |

---

## 1. The A = 2 tensor structure function b₁ — the only place a tensor DIS calculation exists

### 1.1 What is already gated (see `00` E-2, E-3, E-4) and what would strengthen it

| calculation | reference (verified) | gateable form | LiPolGen symbol | status |
|---|---|---|---|---|
| **CDKS** — Cosyn, Dong, Kumano, Sargsian, *Tensor-polarized structure function b₁ in the standard convolution description of the deuteron* | PRD **95** (2017) 074036, `arXiv:1702.05337` (abstract page fetched; "12 pages, 7 eps figures", **no supplementary data, no HEPData record**) | **FIGURE only.** Fig. 4 (x·b₁ᵈ, Q² = 2.5 GeV²) is digitized in-tree as `tables::kB1CdksQ2p5`; Eqs. (16), (17), (19)–(21) are **CLOSED FORM** and are the actual convolution kernel the tree implements | `b1_convolution`, `CdksB1`, `tables::kB1CdksQ2p5` | in-tree, E-2 |
| **Miller** — *Pionic and hidden-color, six-quark contributions to the deuteron b₁* | PRC **89** (2014) 045203, `arXiv:1311.4561` (abstract page fetched; "15 pages, 6 figures", **no tables of b₁ values**) | **FIGURE only.** Fig. 5 digitized as `tables::kB1Miller` | `toy_b1`, `MillerB1` | in-tree, E-3 |
| **Close–Kumano sum rule** ∫b₁ dx = 0 | Close & Kumano, PRD **42** (1990) 2377 (INSPIRE: *"A sum rule for the spin dependent structure function b₁(x) for spin one hadrons"*, 105 citations; no arXiv) | **CLOSED FORM** | `close_kumano_integral` | in-tree, E-4 (reported, not enforced) |

**The honest position on the two b₁ camps:** neither CDKS nor Miller publishes a
table. Both are figure-only, which is why the tree digitizes them, and the
digitization error rides on every downstream comparison. **No HEPData record
exists for either** (HEPData full-text search returns only g₁ measurements and
unrelated b₁(1235) meson records). Nothing in this section can be upgraded from
FIGURE to TABLE by more searching; it can only be upgraded by asking the authors.

### 1.2 What is NOT yet used, and would be a genuine second opinion

| # | calculation | reference (verified) | what it would gate | uncertainty | obtainable | effort | priority |
|---|---|---|---|---|---|---|---|
| T-1 | **Kumano & Kuroki 2026** — *Tensor-polarized parton distribution functions of the deuteron by a convolution model* | `arXiv:2607.09237` (10 Jul 2026; 7 pages, 5 figures; PDF text checked — **no tables**) | An **independent second convolution calculation of δ_T q at exactly Q² = 2.5 GeV²**, i.e. the same point as the digitized CDKS column. Would turn E-2's single-camp shape test into a two-calculation spread, and its stated conclusion — the convolution result is "very different" from the HERMES-fitted PDFs — is the same tension E-2 already reports | none quoted (a model calculation); the paper's own point is the **disagreement with the HERMES fit** | **digitizable figure** (vector PDF) | 0.5–1 d to digitize + a `tables::kB1KumanoKuroki26` column and a spread test beside G3b | **high** |
| T-2 | **Kumano 2010 tensor-polarized PDF fit** — δ_T q_v, δ_T q̄ fitted to HERMES b₁ | PRD **82** (2010) 017501, `arXiv:1005.4524` (abstract fetched) | A **parametrized b₁ᵈ(x, Q²)** — the third camp, and the only one anchored on the measurement rather than on a model. Would give `B1Model` a *data-driven* option to sit beside `Miller` and `Cdks` | carries HERMES's own ±50 % errors, propagated through a 2-parameter fit | **table-in-paper** (fit parameters; 4 pages) | 1–2 d | medium |
| T-3 | **Khan & Hoodbhoy** — *Convenient parametrization for deep inelastic structure functions of the deuteron* | PRC **44** (1991) 1219 (INSPIRE, 55 citations; no arXiv) | The original convolution b₁ as a **few-parameter closed form in wave-function moments** — the form one would refit for an α–d system | none published | **table-in-paper**, but the paper is paywalled (no arXiv, no preprint scan on INSPIRE) | 1 d if a library copy is obtained | low |
| T-4 | **Umnikov** — *Relativistic calculation of b₁,₂(x) of the deuteron* | PLB **391** (1997) 177, `arXiv:hep-ph/9605291` (abstract fetched) | The **validity floor**: shows non-relativistic convolution gets small x wrong and **violates the exact sum rules**. Would justify a documented x ≳ 0.1 validity band on `b1_convolution` rather than the current silent extrapolation | n/a — a statement about the method | **digitizable figure** (3 figures) | 0.5 d for a doc paragraph + an x-range guard | medium |
| T-5 | **Hirai, Kumano, Saito, Watanabe** — *Clustering aspects in nuclear structure functions* | PRC **83** (2011) 035202, `arXiv:1008.1313` (abstract fetched) | **The closest A > 2 analogue that exists.** Not b₁ — it convolves *unpolarized* nucleon SFs with AMD **cluster** momentum distributions for ⁹Be and other light nuclei. It is the methodological precedent for LiPolGen's α–d convolution and the only published check that a cluster light-cone density folded into a DIS convolution behaves sanely | model spread between AMD and shell model is the paper's own error bar | **digitizable figure** | 2–3 d for an unpolarized F₂(⁶Li) cross-check | medium |

**What none of T-1…T-5 validates:** b₁(⁶Li). Every one of them is A = 2 (or, for
T-5, unpolarized A = 9). LiPolGen's ⁶Li number is the same convolution machinery
applied to the α–d overlap through `LI6_B1_RANK2_TRANSFER` and
`LI6_B1_PER_NUCLEON`, and **there is nothing anywhere to compare that step
against**.

---

## 2. Tagged spin-1 DIS — Cosyn–Weiss, the strongest external anchor in the tree

| # | calculation | reference (verified) | gateable form | LiPolGen symbol | status |
|---|---|---|---|---|---|
| T-6 | **Cosyn & Weiss II** — *SIDIS on a polarized spin-1 target. II. Deuteron and spectator nucleon tagging* | `arXiv:2603.23700`, JLAB-THY-26-4661, submitted 24 Mar 2026; 45 pages, 17 figures (abstract page fetched; **Table II read directly out of the local PDF**) | **TABLE.** Table II has three rows × (A_T∥, f₂/f₀, k, θ_k, \|α_p−1\|, p_pT): A_T∥ = **−2** at f₂/f₀ = √2, k = 0.3 GeV, θ_k = 0, \|α_p−1\| = 0.3; A_T∥ = **+1** at f₂/f₀ = √2, k = 0.3 GeV, θ_k = π/2, p_pT = 0.3 GeV; A_T∥ = **+1** at f₂/f₀ = −1/√2, **k = 1 GeV**. AV18 radial wave functions. Eqs. (6.11)–(6.15) are **CLOSED FORM**. **Mapping settled 2026-09-06:** `A_T∥ = +1 · A_zz^wf` exactly (`07_cw_sign_investigation.md` §2, derived twice); CW's −2 is their angular factor at θ_k = 0, a node of the Λ = ±1 densities (Eq. 6.13), **not** a conversion factor. Reading it as one masked an inverted S–D interference sign in `build_amp2`, now fixed | `TaggedModel`, `azz_tensor_curve`, **`deuteron_channel(…, VmcAV18)`** — the gate moved onto the AV18 control on 2026-09-06, because TABLE II is quoted for AV18 and the Hulthén pair's f₂/f₀ never reaches √2 | **in-tree, E-1** |
| T-7 | Cosyn & Weiss I — *…I. Cross section and spin observables* | `arXiv:2603.23699`, JLAB-THY-26-4663, 28 pages, 2 figures (abstract page fetched) | **CLOSED FORM** — the covariant cross-section decomposition; Appendix D is the spin-j multipole counting | the finite-γ decomposition the tree transcribes | partially in-tree (§3 below) |
| T-8 | Cosyn & Weiss 2020 — *Polarized electron-deuteron DIS with spectator nucleon tagging* | PRC **102** (2020) 065204, `arXiv:2006.03033`, 52 pages, 19 figures (abstract page fetched) | **FIGURE + CLOSED FORM.** The LF **polarized spectral function** and nucleon momentum distributions, analytic and numerical | would gate `TaggedModel`'s *vector* sector (`f₀`, `f₂` momentum densities), which E-1 does not touch | **not used** |
| T-9 | Cosyn, Roldan Tomei, Sosa, Zec — *Polarization options in inclusive DIS off tensor polarized deuteron* | EPJ A **61** (2025) 83, `arXiv:2410.12764` (local copy `refs/2410.12764v1.pdf`) | **CLOSED FORM** (Table 1 rows 2, 3) | `tensor_gamma`, `InclusiveKernel` | **in-tree, E-5** |

**T-8 is the one clear, cheap, unexploited win in this section.** E-1 gates the
*tensor* asymmetry against Table II; nothing gates the **vector** tagged sector,
and PRC 102 065204 publishes exactly that, with 19 figures and a closed-form LF
spectral function. Effort: 2–4 d to digitize two or three panels and add a
`test_tagged.cpp` case beside the CW Table II one. **Priority: high.**

**Two honesty notes on T-6.** (i) Cosyn & Weiss themselves disclaim Table II's
third row — the caption says the k = 1 GeV setting "is not presumed to be a
realistic prediction and is listed only for completeness". The in-tree gate uses
rows 1 and 2 only, which is correct; a future gate must not reach for row 3.
(ii) CW II §VI C states the framework is **IA only** and lists the spin
dependence of FSI as an open question — so E-1 validates the tensor asymmetry of
a wave function, in a paper that does not claim FSI is small.

---

## 3. Arbitrary spin and the J = 3/2 (⁷Li) sector

| # | calculation | reference (verified) | gateable form | what it would settle | obtainable | priority |
|---|---|---|---|---|---|---|
| T-10 | **Hoodbhoy, Jaffe, Manohar** — *Novel effects in deep inelastic scattering from spin-1 hadrons* | NPB **312** (1989) 571 (INSPIRE, **233 citations**; no arXiv) | **CLOSED FORM** — b₁…b₄, the (1, −2, 1) helicity pattern, and Δ | already the tree's J = 1 geometry transcription | in-tree (E-7) | — |
| T-11 | **Jaffe & Manohar** — *Deep inelastic scattering from arbitrary spin targets* | NPB **321** (1989) 343 (INSPIRE, **87 citations**; no arXiv, **no preprint scan**) | **CLOSED FORM** — the complete basis for arbitrary J: 2J+1 quark distributions per flavour | the **J = 3/2 rank-2 normalisation**, `00` §9 item 4 / plans/04 #14 | **on request / library copy** — this is the one genuinely hard-to-get item in the survey | high |
| T-12 | **Fu, Sun, Dong** — *Generalized parton distributions in spin-3/2 particles* | PRD **106** (2022) 116012, `arXiv:2209.12161` (abstract fetched: "Eight unpolarized and eight polarized GPDs… In the forward limit… the structure functions and parton distribution functions are obtained") | **CLOSED FORM** — the explicit J = 3/2 forward-limit structure functions as helicity densities | a written-down J = 3/2 b₁ definition to normalise the ⁷Li rank-2 block against. Their "g₂" is a **rank-3** function and collides with the twist-3 nucleon g₂ — rename before use | **arXiv, free** | high |
| T-13 | **Fu, Dong, Kumano, Xie** — *Generalizing the Soffer bound: positivity constraints on PDFs of spin-3/2 particles* | `arXiv:2602.11587` (abstract fetched: *"derive the complete set of positivity bounds for the leading-twist PDFs of a spin-3/2 hadron for the first time"*) | **CLOSED FORM inequalities** | **A free, zero-model gate.** LiPolGen's ⁷Li rank-2 sector and every invented scenario shape (b₃ = 0.05·f₁, b₄ = −0.02·f₁, `toy_delta_gluon`) could be checked to satisfy published positivity bounds. Nothing in the tree currently does this | **arXiv, free** (the inequalities are in the body, not the abstract — must be read out of the PDF) | **high** — cheapest new external gate in the whole survey |
| T-14 | Zhao, Zhang, Liang, Liu, Zhou — covariant spin-3/2 polarization parametrization | PRD **106** (2022) 094006 (`arXiv:2206.11742`); PRD **109** (2024) 074017 (`arXiv:2401.10031`) — UNVERIFIED in this session (carried from `docs/open_items/physics_literature.md`, which reports them as read) | CLOSED FORM — the 3+5+7 = 15-parameter S/T/R decomposition | the geometry `spin` would need for a full J = 3/2 treatment | arXiv, free | medium |
| T-15 | **Cosyn & Weiss I, Appendix D** — extension to higher-spin targets | `arXiv:2603.23699` (verified; local copy in `PolarizedLithiumSim/refs/`) | **CLOSED FORM** counting: spin-j adds an l = 2j multipole; SIDIS structure-function counts 5, 18, 41, **72**, 113 for j = 0, ½, 1, **3/2**, 2 | the parity rule that makes the **rank ≤ 2 truncation exact** for an unpolarized beam — i.e. it turns a caveat into a theorem | arXiv, free | medium |

**What none of this validates:** the ⁷Li rank-2 *magnitude*. T-12/T-13 give a
basis and bounds; the value of the ⁷Li alignment and its structure functions is
still a model. `tests/test_pipeline.cpp` D1 pins the ⁷Li rank-2 sector as
identically zero today, so T-13's bounds would be satisfied trivially — they
become a real gate only once that sector is filled.

---

## 4. The double-helicity-flip gluon distribution Δ ("exotic glue")

This is the least-supported number in the generator and the section where the
gap between what is cited and what is checkable is widest.

| # | source | reference (verified) | what it gives | relation to LiPolGen |
|---|---|---|---|---|
| T-16 | **Jaffe & Manohar** — *Nuclear gluonometry* | PLB **223** (1989) 218 (INSPIRE, **107 citations**; no arXiv) | The definition of Δ as an observable with no nucleonic counterpart | the concept behind `toy_delta_gluon` and the coherent cos 2φ term |
| T-17 | **Sather & Schmidt** — *Size and scaling of the double helicity flip hadronic structure function* | PRD **42** (1990) 1424 (INSPIRE, **31 citations**; no arXiv) — **this is the actual source of `C_BAG`** (`constants.hpp` cites it) | The bag-model sum rule ∫₀¹ x Δ dx = C·α_s(Q²), with LiPolGen's C = −0.012 | `C_BAG`, `toy_delta_gluon` |
| T-18 | **Detmold & Shanahan** — *Gluonic transversity from lattice QCD* | PRD **94** (2016) 014507, Erratum PRD **95** (2017) 079902, `arXiv:1606.04505` (arXiv HTML read) | **TABLE/number.** Reduced matrix element **A₂ ≈ 0.23(2)(5)** for the **φ meson** at m_π = 450(5) MeV — statistical, then a 20 % renormalization estimate; the gluonic **Soffer bound is saturated at 80–100 %** | a magnitude scale for Δ in *any* spin-1 hadron, and a **bound** that any Δ ansatz must respect |
| T-19 | **Kumano & Song** — gluon transversity in polarized pd Drell–Yan | *Gluon transversity in polarized proton-deuteron Drell-Yan process*, PRD **101** (2020) 054011, `arXiv:1910.12523` (arXiv journal-ref verified); and *Deuteron polarizations in the proton-deuteron Drell-Yan process for finding the gluon transversity*, PRD **101** (2020), `arXiv:2003.06623` (INSPIRE-verified; page number not checked). 24 and 19 citations | Cross-section-level estimates for a *deuteron* Δ at Fermilab | the "Drell–Yan estimates" side of `coherent.hpp`'s "bounded by lattice + Drell-Yan estimates" comment |

**The honest reading.** `C_BAG = −0.012` is a **single bag-model number from
1990 with 31 citations**, and the lattice number T-18 is **not a benchmark for
it**: A₂ = 0.23(2)(5) is a *valence-gluon* matrix element in a **φ meson**,
whereas C_BAG parametrises *exotic glue in a nucleus* — gluons not associated
with individual nucleons. They are different objects and a numerical comparison
between them would be exactly the kind of wrong-observable benchmark this survey
exists to prevent. What T-18 **can** legitimately do is supply the **gluonic
Soffer bound**, a genuine inequality that `toy_delta_gluon` should be checked
against. That is the only defensible Δ gate available today.
**Effort: 1–2 d. Priority: medium-high.** Everything else about Δ remains a
scenario, and `RcScope::TensorAll`'s own comment already says nobody has computed
radiative corrections for it.

---

## 5. Ab-initio ⁶Li / ⁷Li structure (VMC, GFMC, NCSMC)

### 5.1 Already in the tree

`data/vmc/` (E-13, E-17) carries the ANL cluster overlaps and momentum
distributions. Their published home is verified: **Wiringa, Schiavilla, Pieper,
Carlson, PRC 89 (2014) 024305, `arXiv:1309.3794`**, whose abstract states
*"nucleon-cluster momentum distributions include … alpha-d in 6Li, alpha-t in
7Li, and alpha-alpha in 8Be. Detailed tables are provided on-line for
download."* The in-tree files **are** that online supplement (fetched through the
Wayback Machine; `data/vmc/README.md` records every URL). Wave-function pedigree:
Pudliner *et al.*, PRC **56** (1997) 1720, `arXiv:nucl-th/9705009`; cluster-overlap
method: Forest *et al.*, PRC **54** (1996) 646.

### 5.2 What exists and is not used

| # | calculation | reference (verified) | gateable form | LiPolGen symbol it would constrain | uncertainty | obtainable | effort | priority |
|---|---|---|---|---|---|---|---|---|
| T-20 | **Wiringa & Schiavilla** — *Microscopic calculation of ⁶Li elastic and transition form factors* | PRL **81** (1998) 4317, `arXiv:nucl-th/9807037` (abstract page fetched **and the PDF text read**) | **FIGURE.** Fig. 1: VMC AV18+UIX F_L(q) and F_T(q) vs data, **with the C0 and C2 multipoles shown separately** and MC error bars. From the text: C2 ≪ C0 below 3 fm⁻¹ and **C2 dominant for q ≥ 3 fm⁻¹** (the shoulder in the data "is entirely due to this component"); F_T's first peak at q = 0.5 fm⁻¹ is reproduced but "the zero comes a little too early" and the second peak at q = 2 fm⁻¹ is overpredicted; **Q(⁶Li) = −0.23(9) fm²** (vs measured −0.08) and **μ = 0.829 (IA) / 0.832 (with MEC) μ_N**, ~1 % above experiment | `HoSpin1FF`, `TabulatedSpin1FF`, `LI6_FF_HO_A_FM`, `LI6_FF_HO_ALPHA`, and the **unfitted first-C0-zero band [2.9, 3.3] fm⁻¹** | the paper's own **50 % statistical error on Q**; the C2 low-q normalisation is wrong by ~2.8× (its own quadrupole vs the measured one) while the high-q shoulder is right | **digitizable figure** — and I checked the mechanics: the arXiv PDF's figures are **vector** (Ghostscript-converted PostScript, only 1×1 stencil masks in `pdfimages -list`), so the sibling repository's `PolarizedLithiumSim/tools/digitize_figure.py` applies verbatim — its own docstring already names CDKS Figs. 4/5 and Miller Figs. 5/6 as targets drawn "as PDF path operators, not as a raster". It is a dev-time script (PyMuPDF), the CSV is what ships, and provenance goes in that repo's `SOURCES.md`; the table itself would land in LiPolGen's `data/ff/li6_elastic.csv` | **1–2 d** | **HIGHEST in this file** |
| T-21 | **Pastore, Pieper, Schiavilla, Wiringa** — QMC EM moments, A ≤ 9, with χEFT MEC | PRC **87** (2013) 035503, `arXiv:1212.3375` (abstract fetched) | **TABLE** — GFMC AV18+IL7 moments | `LI6_QUADRUPOLE_FM2`'s model counterpart; already quoted in `00` C-2 as Q = −0.20(6) fm² | 30 % | **table-in-paper** | 0.5 d | medium (mostly already used) |
| T-22 | **Piarulli, Pastore, Wiringa** *et al.* — *Densities and momentum distributions in A ≤ 12 nuclei from chiral EFT interactions* | `arXiv:2210.02421` (abstract fetched; PRC **107** (2023) 014314) | **TABLE/data files** — the same observables from **five Norfolk NV2+3 Hamiltonians** | turns E-13's single-Hamiltonian VMC input into a **Hamiltonian band**: the AV18+UX vs chiral spread is the honest systematic on P_D(⁶Li) = 1.94 %, ⟨k⟩ and the polarizations | inter-model spread ≲ 0.01 on polarizations | **table-in-paper**; data files may need a request | 2–3 d | **high** — it is the only way to price the wave-function systematic with something other than the MC error |
| T-23 | **Hebborn, Brune, Phillips** — *Connecting ground-state properties of ⁶Li to each other and to scattering data* | `arXiv:2510.19067` (abstract fetched; 28 pages, 7 figures, submitted 21 Oct 2025) | **FIGURE** (correlations). NCSMC ab-initio; the paper's stated result is a strong ANC²–separation-energy correlation over S_d = 1.3–1.8 MeV | would put an **ab-initio prior** on the α–d asymptotics that C-1's η gate currently reads only from the George–Knutson *measurement*. **The abstract does not promise a D/S ratio or Q correlation** — I could not confirm that η or Q appear, so this is a "read the paper first" item, not a ready gate | not stated in the abstract | **arXiv, free**; the specific numbers must be read out | 1 d to read, then 1–2 d if η is in there | medium |
| T-24 | **Brida, Pieper, Wiringa** — QMC spectroscopic overlaps, A ≤ 7 | PRC **84** (2011) 024319, `arXiv:1106.3121` (abstract fetched) | **TABLE (online)** — but **single-nucleon overlaps only** | nothing new: `data/vmc/README.md` already establishes that the modern ANL overlap pages carry no two-cluster (α+d, α+t) overlaps, which is why the tree uses the 2004 `overlap_old` files | — | in-tree conclusion, verified | — | — |

### 5.3 What does not exist

- **GFMC or VMC structure functions for ⁶Li.** None. QMC computes moments,
  densities, momentum distributions and form factors — not DIS structure
  functions.
- **GFMC Euclidean electromagnetic response functions for ⁶Li.** Verified absent:
  the GFMC response programme covers A = 3, 4 and ¹²C (Lovato *et al.*), and ⁶Li
  appears in the literature only as a *target of data-driven neural-network
  models* (Sobczyk, Rocco, Lovato, `arXiv:2406.06292`, which trains on the
  quasielastic e–A database and lists ⁶Li among its targets). So the
  quasi-elastic piece of the radiative tail — LiPolGen's `RC_QE_KF_GEV` Fermi-gas
  model (C-10) — has **no ab-initio ⁶Li counterpart**, and the only realistic
  upgrade path is that data-driven response, not a QMC calculation.
- **A tabulated coordinate-space ⟨⁴He+d\|⁶Li⟩ overlap** beyond the 2004
  `overlap_old` files already in the tree.
- **An α–d D/S ratio from an ab-initio calculation** that could replace C-1's
  "gate a documented 1.93× discrepancy" with an agreement test — pending T-23.

---

## 6. Glauber / FSI for tagged spectators

| # | calculation | reference (verified) | gateable form | LiPolGen symbol | what it does NOT give |
|---|---|---|---|---|---|
| T-25 | **Cosyn & Sargsian** — *Nuclear FSI in DIS off the lightest nuclei* | IJMPE **26** (2017) 1730004, `arXiv:1704.06117` (abstract fetched) | **CLOSED FORM + fitted parameters**: the eikonal amplitude f_XN = σ_tot(Q², M_X)(i + ε)e^{Bt/2}, ε = −0.5, B ≈ 6 GeV⁻², σ_tot fitted to JLab "Deeps" | `GlauberFsiWeight` (in-tree, E-15) | any ⁶Li number; the fit is deuteron, nucleon-spectator |
| T-26 | Cosyn & Sargsian 2011; Cosyn, Melnitchouk, Sargsian 2014 | PRC **84** (2011) 014601, `arXiv:1012.0293`; PRC **89** (2014) 014612, `arXiv:1311.3550` (both abstracts fetched) | **FIGURE** — GEA FSI vs JLab data, including the forward-angle rise | would gate the *shape* (θ_k, α_s dependence) of `GlauberFsiWeight`, which E-15 checks only for analytic limits | same |
| T-27 | **Ciofi degli Atti & Kaptari** | PRC **83** (2011) 044602, `arXiv:1011.5960` (abstract fetched) | **CLOSED FORM** — the per-nucleon survival product and Gaussian profile; **the only published *cluster*-spectator FSI**, ³He(e,e′d)X | the construction `GlauberFsiWeight` implements | it is A = 3 with a **d** spectator; no α or t spectator anywhere |
| T-28 | **Kaptari, Del Dotto, Pace, Salmè, Scopetta** — *Distorted spin-dependent spectral function of an A=3 nucleus* | PRC **89** (2014) 035206, `arXiv:1307.2848` (abstract fetched) | **FIGURE + CLOSED FORM** — the **polarized** cluster-spectator distorted spectral function | the closest published object to "FSI acting on a tagged polarized cluster", i.e. the missing spin dependence of `GlauberFsiWeight` | A = 3, deuteron spectator; and it is SIDIS, not tensor |
| T-29 | **Strikman & Weiss** | PRC **97** (2018) 035209, `arXiv:1706.02244` (abstract fetched) | **CLOSED FORM** — the factorized distorted spectral function S[IA] + S[FSI] + S[FSI²], the unitarity sum rule ∫dΓ S[FSI+FSI²] = 0, and **FSI vanishing at the pole** | the second Glauber variant already pinned in E-15, and the pole-extrapolation argument for the α-tagged channel | a ⁶Li number |
| T-30 | **Sargsian & Strikman** — the no-loop / pole-extrapolation theorem | PLB **639** (2006) 223, `arXiv:hep-ph/0511054` (abstract fetched) | **CLOSED FORM** | justifies the α-tagged pole-extrapolation diagnostic (S_αd = 1.474 MeV, comparable with the deuteron's 2.22 MeV) | — |
| T-31 | **Sargsian** — GEA foundations | IJMPE **10** (2001) 405, `arXiv:nucl-th/0110053` (INSPIRE, 117 citations) | CLOSED FORM | the eikonal formalism itself | — |
| T-32 | **`github.com/wcosyn/physics-code`** — Cosyn's released code | repository fetched, **HTTP 200** | **CODE** — `DIS/DeuteronCross`, `Glauber/FastParticle`, with a self-contained energy-dependent NN parameterization | a **runnable second implementation** of the deuteron tagged FSI: the only way to turn E-15 from "analytic limits + a Python-prototype self-pin (B-2)" into a genuine code-vs-code diff | it is deuteron/nucleon-spectator; it has no α spectator and no tensor sector |

**The one recommendation here.** T-32 is the highest-value FSI item: E-15's
thirteen test cases check **limits and symmetries**, and B-2 pins the C++ against
this project's own Python prototype. Running `physics-code` at a matched deuteron
kinematic point and comparing the distorted-to-IA ratio would be the tree's
**first external FSI number**. Effort: 3–5 d (build + a matched configuration +
one comparison table). **Priority: high.** It validates the deuteron control
channel only — the α-spectator weight remains unvalidated by construction,
because no calculation of it exists.

---

## 7. Radiative corrections

| # | source | reference (verified) | gateable form | LiPolGen symbol | status |
|---|---|---|---|---|---|
| T-33 | **POLRAD 2.0** | CPC **104** (1997) 201, `arXiv:hep-ph/9706516` (abstract fetched) | **CODE** — carries the full spin-1 tensor Born sector (Eqs. 9/10 with the −1/3 : +1/6 = P₂ ratio), the tensor elastic tail via F_q, and the b₁–b₄ kernels | `rc_delta`, `HoSpin1FF`, the σ^el gates | **in-tree, D-3** — compiled standalone with its own `ffdeu` |
| T-34 | **Gakh & Shekhovtsova** — *RC to DIS ed scattering: case of tensor polarized deuteron* | JETP **99** (2004) 898, `arXiv:hep-ph/0403262` (abstract fetched: "model-independent RC … unpolarized electron beam off the tensor polarized deuteron target … including the elastic radiative tail") | **FIGURE** (5 figures) | `RC_DELTA_LOW_X` — the band half-width | in-tree; **zero INSPIRE citations**, i.e. never independently checked. The tree already labels the shape as rejected and uses only the band width |
| T-35 | **Mo & Tsai** | RMP **41** (1969) 205 (INSPIRE, **1088 citations**) — **paywalled at APS, no arXiv, no free scan of the RMP article** | — | the gate `tests/test_rc.cpp` says was deliberately not written | **paywalled** |
| T-36 | **Tsai, SLAC-PUB-848** — *Radiative corrections to electron scatterings* (1971) | INSPIRE recid **67278**, 93 citations. **I verified the full text is freely downloadable**: `https://inspirehep.net/files/3ce239706be17b0eefec145433700c64` returns HTTP 200, `application/pdf`, 4.1 MB | **CLOSED FORM** — Tsai's own expanded treatment of the Mo–Tsai formalism, including the equivalent-radiator and exact tail formulas | the same gate T-35 blocks | **THIS CLOSES AN OPEN ITEM.** `00` §9 item 5 records Mo–Tsai as "named … as the gate that was deliberately not written and was not obtained for this run". It is obtainable — not as RMP 41 205, but as the freely-hosted SLAC report by the same author, which is what most of the community actually cites for the formulas |
| T-37 | **Afanasev, Bernauer, Blunden, Blümlein** *et al.* — *Radiative corrections: from medium to high energy experiments* | `arXiv:2306.14578` (abstract fetched; EPJ A topical review, 63 pages, 27 figures, **3 tables**) | **TABLE + review** — the modern state of QED corrections in DIS and the associated Monte-Carlo codes | a citable modern baseline for what "1.5 % on A_zz" means, and the place to find which codes are maintained | **arXiv, free** |
| T-38 | Akushevich & Ilyichev, next-to-leading-order RC codes | Verified via INSPIRE: PRD **90** (2014) 013006-class work on exclusive photon electroproduction with next-to-leading accuracy (`arXiv:1403.3421`), and NLO mass-effect papers (`arXiv:hep-ph/0101126`). **No spin-1/tensor NLO code found** | — | — | **the tensor sector has no modern successor to POLRAD**; the searches return only nucleon-target work |

**Net position on RC.** Two things change from `00` §9 item 5. (a) **Mo–Tsai is
obtainable** in the form the community uses (T-36) — the "not obtained" status
should be retired. (b) The *tensor* RC situation is unchanged and remains the
weakest link: one 2004 calculation with zero citations (T-34), no modern code
(T-38), and nothing at all for the φ-dependent (cos φ_TL, cos 2φ_TT) observables
or for J = 3/2. `RcOptions::qe_tensor_scale`'s "borrowed magnitude and not a
derived bound" label is the correct one and should stay.

---

## 8. Nuclear PDFs — the one place ⁶Li is genuinely over-supplied

I enumerated the **live LHAPDF set index** (`lhapdfsets.web.cern.ch/current/`,
1590 sets) rather than trusting memory. Grids that exist for A = 6, Z = 3 and
A = 7, Z = 3:

| family | ⁶Li set(s) | ⁷Li set(s) | reference (verified) |
|---|---|---|---|
| **EPPS21** | `EPPS21nlo_CT18Anlo_Li6` | **none** | Eskola, Paakkinen, Paukkunen, Salgado, EPJC (2022), `arXiv:2112.12462` |
| **EPPS16** | `EPPS16nlo_CT14nlo_Li6` | **none** | — |
| **nCTEQ15** | `nCTEQ15_6_3` **+ 6 variant families** (`FullNuc`, `HIX`, `HQ`, `np`, `WZ`, `WZSIH`, each ± `FullNuc`) | **`nCTEQ15_7_3` + the same 6 variants** | Kovarik *et al.*, PRD **93** (2016) 085037, `arXiv:1509.00792` |
| **nNNPDF** | `nNNPDF10_nlo/nnlo_as_0118_Li6`, `nNNPDF20_nlo_as_0118_Li6`, `nNNPDF30_nlo_as_0118_A6_Z3` | **none** | Abdul Khalek *et al.* (nNNPDF3.0), `arXiv:2201.12363` |
| **TUJU19** | none | **`TUJU19_nlo_7_3`, `TUJU19_nnlo_7_3`** | UNVERIFIED reference (set names verified in the index; I did not fetch the TUJU19 paper) |

| # | opportunity | what it gates | what it does NOT gate | effort | priority |
|---|---|---|---|---|---|
| T-39 | **A four-fit ⁶Li EMC cross-check.** EPPS21 (in-tree, D-5) vs nCTEQ15 vs nNNPDF3.0 vs EPPS16 on the same x, Q² grid | The **unpolarized nuclear modification** of ⁶Li as a *fit spread* rather than a single central curve. Today `Epps21Ratio` carries one fit and its Hessian errors; the fit-to-fit spread is a larger and unmeasured systematic on every unpolarized rate in the generator | anything polarized, anything tensor. And an important caveat: **almost no data constrains A ≈ 6** — the ⁶Li grids are evaluations of an A-dependent parametrization, anchored on one DIS ratio measurement (NMC ⁶Li/D), so four "independent" fits share one dataset and their agreement is weaker evidence than it looks | 2–3 d (all four sets are LHAPDF-installable; the tree already has the wrapper) | **high** |
| T-40 | **⁷Li nuclear PDFs exist after all** — `nCTEQ15_7_3` and `TUJU19_*_7_3` | This is a **correction to an assumption**: ⁷Li is not PDF-less. Anything in the tree that treats ⁷Li's unpolarized baseline as unavailable, or reuses the ⁶Li ratio for ⁷Li, can be replaced with a real grid | same caveat as T-39, more strongly: there is no ⁷Li DIS ratio measurement at all, so `_7_3` is pure A-interpolation | 1 d | medium |
| T-41 | **Polarized nuclear PDFs** | — | **They do not exist for any nucleus**, lithium included. INSPIRE returns zero fits under "polarized nuclear parton distributions". The nearest objects are (i) the ³He g₁ convolution analysis — Bissey, Guzey, Strikman, Thomas, PRC **65** (2002) 064317, `arXiv:hep-ph/0109069` (verified, 74 citations) — and (ii) the polarized-EMC model curves already digitized in-tree (E-11). LiPolGen's `P_p`/`P_n` route is the right one and has no PDF-level alternative | — | — |

---

## 9. Coherent diffraction on light nuclei

| # | calculation | reference (verified) | gateable form | LiPolGen symbol | what it does NOT give |
|---|---|---|---|---|---|
| T-42 | **Mäntysaari, Salazar, Schenke, Shen, Zhao** — *Spatial imaging of polarized deuterons at the EIC* | PLB **858** (2024) 139053, `arXiv:2408.13213` (local copy) | **FIGURE** — a₂(\|t\|) per m-state | `a2_from_quadrupole` (in-tree, E-12) | any A > 2. **Re-verified 2026-09-06: 8 forward citations on INSPIRE, and none of them extends the calculation to A > 2 or to another polarized nucleus.** It remains the only polarization-dependent coherent-diffraction calculation in existence |
| T-43 | **Mäntysaari, Roch, Schenke, Shen, Zhao** — *Nuclear structure and saturation effects from diffractive vector meson production* | PRD **114** (2026) 014068, `arXiv:2605.00454` (abstract fetched **and the arXiv HTML read for the nucleus list**) | **FIGURE** — coherent/incoherent J/ψ; the nucleus/model map is ³He and ⁴He with **GFMC** configurations, ¹²C with VMC, ¹⁶O with VMC/PGCM/NLEFT, ²⁰Ne with PGCM/NLEFT, ⁴⁰Ar with NLEFT, A > 40 Woods–Saxon | the **α-clustering** discussion is about ¹⁶O and ²⁰Ne | **NO ⁶Li and NO ⁷Li.** And the observable is **ultraperipheral AA photoproduction**, not e+A. A repo note that cites this paper for "A = 3, 4" is half right — GFMC ³He/⁴He configurations are used — but it is not a light-ion e+A paper and contains no lithium |
| T-44 | **Guzey, Rinaldi, Scopetta, Strikman, Viviani** — *Coherent J/ψ electroproduction on ⁴He and ³He at the EIC* | PRL **129** (2022) 242503, `arXiv:2202.12200` (abstract fetched; comments record **"Supplemental Material"**) | **FIGURE + possible TABLE in the supplement** — coherent J/ψ t-distributions on ³He/⁴He with **realistic wave functions** | an independent, non-eSTARlight, non-CGC baseline for coherent VM production on a *light* nucleus — the closest published cross-check for D-4's eSTARlight ⁶Li rates | unpolarized; A = 3, 4 only; no cluster structure of the ⁶Li kind |
| T-45 | **`github.com/hejajama/subnucleondiffraction`** — the code behind T-42 | repository fetched, **HTTP 200** | **CODE** — A-generic `Nucleons` amplitude, a VMC wave-function branch, Good–Walker averaging over an external configuration table | the only realistic route to a **tensor-polarized ⁶Li** coherent amplitude: replace the pn density with an α+d configuration sampler | it does not do this today; the author asks to be contacted before use |
| T-46 | **Tensor-polarized coherent diffraction for any A > 2** | — | **DOES NOT EXIST.** Checked twice: the 8 forward citations of T-42 (listed above, September 2026) and the light-nucleus diffraction literature | this is why `eps_b0 = −0.08` is a **deuteron** number applied to ⁶Li and why T10b prints that it corresponds to a ⁶Li quadrupole **11.4× the measured one** | — |

**Recommendation.** T-44 is the cheap one: 1–2 d to pull the supplemental
material and compare a t-slope against `estarlight_li6_coherent()`'s. It would
give D-4 its first non-eSTARlight point of contact. T-45 is the expensive one and
is correctly scoped in `docs/open_items/` as a collaboration, not an
implementation.

---

## 10. Priority table — what to do, in order

| rank | item | why it is first | effort | what it would change |
|---|---|---|---|---|
| 1 | **T-20** digitize Wiringa–Schiavilla Fig. 1 (⁶Li F_C0 / F_C2 / F_M1) | The slot is already built and empty: `TabulatedSpin1FF::from_data_dir("ff/li6_elastic.csv", …)` is implemented, `data/ff/` does not exist, and `tests/test_rc.cpp` has a TEST_CASE titled *"the API exists, and NO 6Li table is shipped"* whose message asserts "FIGURES ONLY". I verified the figures are **vector** (`pdfimages -list` finds only 1x1 stencil masks), so the sibling repo's `digitize_figure.py` applies. It would replace the **unfitted** C0-zero guess q₀ = 3.1 fm⁻¹ (band [2.9, 3.3]) with an ab-initio shape | 1–2 d | closes `00` §9 item 9 at the *theory* level (the *data* level still needs elastic ⁶Li measurements) |
| 2 | **T-13** apply the spin-3/2 positivity bounds to the ⁷Li sector and every scenario shape | Free, zero-model, closed-form inequalities. Nothing in the tree currently checks that b₃, b₄, Δ scenarios are even *allowed* | 1–2 d | first-ever external constraint on the invented shapes |
| 3 | **T-32** run `wcosyn/physics-code` at a matched deuteron point and diff `GlauberFsiWeight` | Turns the FSI sector from "analytic limits + a self-pin" into a code-vs-code comparison | 3–5 d | first external FSI number in the tree |
| 4 | **T-1** digitize Kumano–Kuroki 2026 at Q² = 2.5 GeV² | A second, brand-new convolution b₁ᵈ at exactly the digitized CDKS point | 0.5–1 d | E-2 becomes a two-calculation spread instead of one camp |
| 5 | **T-39** four-fit ⁶Li nPDF spread | Four ⁶Li grids are installable today; the fit-to-fit spread is currently unmeasured | 2–3 d | prices the unpolarized nuclear systematic honestly (with the shared-data caveat stated) |
| 6 | **T-22** Norfolk chiral wave functions as a Hamiltonian band on `data/vmc` | The only way to price the wave-function systematic with something other than the ANL MC error | 2–3 d | E-13's "Hamiltonian difference, not an error" becomes a number |
| 7 | **T-8** digitize the Cosyn–Weiss 2020 polarized LF spectral function | E-1 gates the tensor tagged sector; the **vector** tagged sector has no external gate at all | 2–4 d | closes the other half of the tagged channel |
| 8 | **T-36** obtain Tsai SLAC-PUB-848 and retire the "Mo–Tsai not obtained" status | Verified freely downloadable from INSPIRE | 0.5 d + whatever gate is then written | corrects a documented gap |
| 9 | **T-18** gluonic Soffer bound on `toy_delta_gluon` | The only defensible Δ gate that exists | 1–2 d | the Δ sector stops being entirely ungated |
| 10 | **T-44** Guzey *et al.* supplemental t-slopes vs `estarlight_li6_coherent()` | First non-eSTARlight coherent light-nucleus contact | 1–2 d | D-4 gains an external cross-check |
| 11 | **T-11** obtain Jaffe–Manohar NPB 321 (library copy) | The J = 3/2 normalisation open item (`00` §9 item 4) cannot close without it | on request | unblocks plans/04 #14 |

---

## 11. After all of it, what is still not validated

Doing every item above would leave these untouched, and they are the ones the
generator exists for:

1. **b₁(⁶Li), b₁(⁷Li).** No calculation exists for any A > 2. Verified by
   INSPIRE search in this session. LiPolGen would be the first, with nothing to
   check it against; the α–d rank-2 transfer step (`LI6_B1_RANK2_TRANSFER`) is
   validated by *nothing*.
2. **b₂, b₃, b₄ and Δ for lithium.** Unmeasured and uncalculated. T-13's
   positivity bounds would say the scenarios are *allowed*, never that they are
   *right*.
3. **Tagged tensor asymmetries for an α or t spectator.** T-27/T-28 do a
   *deuteron* spectator in A = 3. No α-spectator FSI calculation exists.
4. **Tensor-sector radiative corrections at A = 6**, and any RC at all for the
   φ-dependent tensor observables (T-34, T-38).
5. **Tensor-polarized coherent diffraction for A > 2** (T-46).
6. **The J = 3/2 rank-2 normalisation** until T-11 and T-12 are transcribed.

The pattern is the same one `00` §10 records: the external theory support is
real, it is dense at A = 2 and on ab-initio ⁶Li *structure*, and it stops at the
exact boundary where tensor polarization meets A > 2. Everything in this file is
a component check, and no accumulation of component checks turns into a
validation of the observable.
