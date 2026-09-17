<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# 01 — TENSOR-POLARIZED DIS AND b₁: the primary sources for the tensor sector

**Research pass 2026-09-16.** Domain file of the `docs/references/` series;
companion to `00_corpus.md` (what is already held) and to
`../benchmarking/04_theory.md` (which calculations could be gated, and how).

This file answers one question for the tensor sector: **for every number the
generator carries in its tensor block — or every number it is missing — which
paper is the primary source, is that paper obtainable, and which LiPolGen
symbol or gate does it serve?**

**Verification protocol.** Every entry below was checked in this session against
the arXiv API (`export.arxiv.org/api/query`, title + authors + journal-ref +
DOI + abstract), the INSPIRE literature API (`inspirehep.net/api/literature`,
title + authors + publication_info + DOI + citation count), or — where a
quantitative claim is made about a paper's *contents* — against the text of the
PDF itself, extracted with `pdftotext`. Nothing here is repeated from memory.
Claims sourced from a PDF body rather than an abstract are marked
**(PDF-read)**. Items that could not be verified, or could not be obtained
without a paywall, are marked so explicitly.

**What was downloaded.** **34** new open-access files landed in the sibling
corpus directory `../../../PolarizedLithiumSim/refs/` during this pass: **29**
arXiv e-prints, **3** IOP *J. Phys. Conf. Ser.* open-access articles, and the
**2** public JLab proposal PDFs. Paths are given per entry. One further file,
`hep-ex_0506018.pdf` (HERMES), appeared in that directory during this session
from a concurrent domain pass and is used here as found — it is still absent
from `refs_dict.json`. No paywall was approached; the closed items are listed
in §7 with the exact DOI page to buy or borrow from.

---

## 0. The headline finding, stated first

> **There is no calculation of b₁, of any tensor structure function, or of any
> tensor-polarized DIS observable, for any nucleus with A > 2 — anywhere in the
> literature.** Every tensor-DIS calculation that exists is for the deuteron
> (A = 2) or for a meson (the φ, on the lattice).

This was tested in this session with eight INSPIRE queries designed to fail
loudly if such a paper existed (`t "b_1" and t nuclei`; `t tensor and t
"structure function" and t "Li"`; `t "spin-1 nuclei" and t "deep inelastic"`;
`t "tensor" and t "deep inelastic" and t "nucleus"`; `ft "b_1" and ft "6Li" and
t tensor`; `ft "b_1" and t "nitrogen"`; `t "spin-1 nucleus" or t "spin-one
nucleus"`; `a Kumano and t "nuclear" and t "tensor"`). **All eight returned
zero hits.** A ninth, `ft "tensor polarized" and ft "lithium"`, returns 12
records, every one of which is a *target-material* mention (⁶LiD polarized
targets, EIC polarized-beam programmes, EDM storage rings) and not a structure
calculation.

This confirms, independently and by a different route, the statement already
standing in `../benchmarking/04_theory.md` §0 item 1. It is the single most
important fact about this domain and it is a **finding, not an absence of
effort**: LiPolGen's b₁(⁶Li) — `Li6B1`, built as
`LI6_B1_RANK2_TRANSFER × LI6_B1_PER_NUCLEON × b₁ᵈ` — has no external
counterpart and cannot acquire one by more searching. The honest framing for
any paper is that the A = 2 machinery is validated and the A = 6 step is a
documented, single-assumption extrapolation.

**Two partial exceptions worth knowing**, both new to this project's corpus:

1. **Lattice QCD has now computed a gluonic tensor quantity in a *nucleus*,**
   not just in a meson — Winter *et al.* (NPLQCD), `arXiv:1709.00395`, gives
   the deuteron gluon-transversity moment and a bound on the gluonic-b₁
   combination (§5.2). That closes the wrong-observable objection that
   `../benchmarking/04_theory.md` §4 correctly raises against using the φ-meson
   number for a nuclear Δ.
2. **Spin-3/2 now has a complete, published PDF basis and a complete set of
   positivity bounds** (§6). LiPolGen's ⁷Li rank-2 sector is identically zero
   today, so the bounds are satisfied trivially — but the moment that sector is
   filled they become a free, zero-model gate.

---

## 1. Table 1 — the tensor-b₁ domain, complete

`OA` = open access. `local` = path under `PolarizedLithiumSim/refs/`
(**new** = downloaded in this pass). `serves` = the LiPolGen symbol, gate or
open item the reference answers to.

### 1a. Foundations — the formalism the tree transcribes

| key | reference | id | OA | local | serves | tier |
|---|---|---|---|---|---|---|
| hoodbhoy-jaffe-manohar1989 | Hoodbhoy, Jaffe, Manohar, *Novel effects in DIS from spin-1 hadrons*, NPB **312** (1989) 571 | DOI 10.1016/0550-3213(89)90572-5 | ✗ | — (in corpus, `file: null`) | b₁–b₄ and Δ definitions; the (1, −2, 1) helicity pattern; `TensorSF` geometry, E-7 | theory-input |
| jaffe-manohar1989-arbitrary-spin | Jaffe, Manohar, *DIS from arbitrary spin targets*, NPB **321** (1989) 343 | DOI 10.1016/0550-3213(89)90347-7 | ✗ | — (in corpus, `file: null`) | the J = 3/2 basis for the ⁷Li rank-2 block | theory-input |
| **jaffe-manohar1989-gluonometry** | **Jaffe, Manohar, *Nuclear gluonometry*, PLB 223 (1989) 218** | DOI 10.1016/0370-2693(89)90242-6 | ✗ | — **GAP** | **the definition of Δ as "exotic glue"; the concept behind `toy_delta_gluon`** | theory-input |
| **close-kumano1990-sum-rule** | **Close, Kumano, *A sum rule for the spin-dependent structure function b₁(x) for spin-one hadrons*, PRD 42 (1990) 2377** | DOI 10.1103/PhysRevD.42.2377 | ✗ | — **GAP** | **`close_kumano_integral` — the ∫b₁ dx = 0 gate (E-4) cites a paper the corpus does not hold** | T5 |
| **khan-hoodbhoy1991** | **Khan, Hoodbhoy, *Convenient parametrization for DIS structure functions of the deuteron*, PRC 44 (1991) 1219** | DOI 10.1103/PhysRevC.44.1219 | ✗ | — **GAP** | the original convolution b₁ as a few-parameter closed form in wave-function moments — the form one would refit for α–d | theory-input |
| khan-hoodbhoy1993-shadowing | Khan, Hoodbhoy, *Shadowing of deuteron spin structure functions*, PLB **298** (1993) 181 | DOI 10.1016/0370-2693(93)91727-5 | ✗ | — | the same authors' small-x sequel; b₁ from double scattering in the LPS model | theory-input |
| kumano1993-b1-proposal | S. Kumano, *Tensor structure function b₁(x) for spin-one hadrons*, MKPH-T-93-03 (1993) | `arXiv:hep-ph/9302320` | ✓ | **new** `hep-ph_9302320.pdf` | the earliest statement of the b₁–b₄ measurement case; historical head of the Kumano series | theory-input |
| bacchetta-mulders2000 | Bacchetta, Mulders, *Deep inelastic leptoproduction of spin-one hadrons*, PRD **62** (2000) 114004 | `arXiv:hep-ph/0007120`, DOI 10.1103/PhysRevD.62.114004 | ✓ | **new** `hep-ph_0007120.pdf` | the spin-1 TMD basis (f₁LL, h⊥₁LL…) and the Trento-convention names Kumano–Kuroki convert into | theory-input |

### 1b. The measurement, and the experimental programme

| key | reference | id | OA | local | serves | tier |
|---|---|---|---|---|---|---|
| **hermes2005-b1** | **HERMES Collab. (A. Airapetian *et al.*), *First measurement of the tensor structure function b₁ of the deuteron*, PRL 95 (2005) 242001** | `arXiv:hep-ex/0506018`, DOI 10.1103/PhysRevLett.95.242001 | ✓ | `hep-ex_0506018.pdf` | **the only tensor-DIS datum in the world.** Table II (6 rows, PDF-read below) is BENCHMARK_PLAN §4 item 3 | **T2** |
| **jlab-pr12-13-011-proposal** | K. Slifer *et al.*, *The Deuteron Tensor Structure Function b₁*, PR12-13-011, JLab PAC-40 (2013) | jlab.org/exp_prog/proposals/13/PR12-13-011.pdf | ✓ | **new** `JLab_PR12-13-011.pdf` | the approved b₁ experiment; 0.16 < x < 0.49, 0.8 < Q² < 5.0 GeV², 30 days at 11 GeV, P_zz = 20 % (PDF-read). **Corpus entry `jlab-pr12-13-011` has `file: null` — now fillable** | T2 |
| **jlab-pr12-15-005-proposal** | E. Long, K. Slifer, P. Solvignon *et al.*, *Measurements of the Quasi-Elastic and Elastic Deuteron Tensor Asymmetries*, PR12-15-005, JLab PAC 43 (2015) | jlab.org/exp_prog/proposals/15/PR12-15-005.pdf | ✓ | **new** `JLab_PR12-15-005.pdf` | the A_zz proposal: quasi-elastic A_zz at x > 1 plus elastic T₂₀ over 0.2 < Q² < 1.8 GeV² (PDF-read) | T2 |
| slifer2014-b1-jpcs | K. Slifer, *The deuteron polarized tensor structure function b₁*, J. Phys. Conf. Ser. **543** (2014) 012003 | DOI 10.1088/1742-6596/543/1/012003 | ✓ | **new** `JPCS_543_012003.pdf` | the 4-page public summary of PR12-13-011 — the citable version | T2 |
| long2014-azz-jpcs | E. Long, *Potential for a tensor asymmetry A_zz measurement in the x > 1 region at JLab*, J. Phys. Conf. Ser. **543** (2014) 012010 | DOI 10.1088/1742-6596/543/1/012010 | ✓ | **new** `JPCS_543_012010.pdf` | the citable version of the A_zz case; names Frankfurt–Strikman 1988 as the first QE A_zz calculation | T2 |
| kalantarians2014-eic-jpcs | N. Kalantarians, *Tensor polarized deuteron at an electron-ion collider*, J. Phys. Conf. Ser. **543** (2014) 012008 | DOI 10.1088/1742-6596/543/1/012008 | ✓ | **new** `JPCS_543_012008.pdf` | the earliest written case for tensor polarization *at an EIC* — the direct ancestor of this project's premise | T2 |
| long-higinbotham-solvignon2014-proceedings | E. Long, D. Higinbotham, P. Solvignon (eds.), *Proceedings, Tensor Polarized Solid Target Workshop* (JLab, 10–12 Mar 2014), J. Phys. Conf. Ser. **543** (2014) | iopscience.iop.org/1742-6596/543/1 (INSPIRE recid 1324998) | ✓ | — (volume; individual papers held) | the volume the four entries above and Cosyn–Sargsian all sit in | T2 |
| slifer-long2013-pstp | K. Slifer, E. Long, *Novel physics with tensor polarized targets*, PoS **PSTP2013** (2013) 008 | `arXiv:1311.4835`, DOI 10.22323/1.182.0008 | ✓ | **new** `1311.4835.pdf` | the programme-level motivation; target-technology context for ⁶LiD/ND₃ | T2 |
| poudel2025-clas12-tensor-tmd | J. Poudel, A. Bacchetta, J.-P. Chen, D. Keller, I. Fernando, E. Long, D. Ruth, N. Santiesteban, K. Slifer, *Spin-1 transverse-momentum-dependent tensor structure functions in CLAS12* (2025) | `arXiv:2502.20044` | ✓ | **new** `2502.20044.pdf` | the *next* b₁ datum: an RG-C ND₃ analysis extracting b₁ inclusively and F_U(LL),T, F^{cos2φ}_U(LL) in SIDIS | T2 |
| **maxwell2018-exotic-glue-loi** | J. Maxwell, D. Crabb, D. Day, W. Detmold, R. Jaffe, M. Jones, C. Keith, D. Keller, D. Meekins, R. Milner *et al.*, *Search for exotic gluonic states in the nucleus*, Letter of Intent to JLab PAC 44, dated 6 Jun 2016, posted 2018 | `arXiv:1803.11206` | ✓ | **new** `1803.11206.pdf` | **the only experimental programme for Δ.** Unpolarized beam + transversely polarized spin-1 target, x < 0.3 | T2 |

### 1c. Calculations of b₁ᵈ — the camps

| key | reference | id | OA | local | serves | tier |
|---|---|---|---|---|---|---|
| cosyn2017-b1-convolution | Cosyn, Dong, Kumano, Sargsian, PRD **95** (2017) 074036 | `arXiv:1702.05337` | ✓ | `1702.05337.pdf` (in corpus) | `b1_convolution`, `CdksB1`, `tables::kB1CdksQ2p5`; E-2 | theory-input |
| miller2014-b1-pion-hidden-color | G. A. Miller, PRC **89** (2014) 045203 | `arXiv:1311.4561` | ✓ | `1311.4561.pdf` (in corpus) | `toy_b1`, `MillerB1`, `tables::kB1Miller`; E-3 | theory-input |
| **nikolaev-schafer1997** | **Nikolaev, Schäfer, *Nonvanishing tensor polarization of sea quarks in polarized deuterons*, PLB 398 (1997) 245; Erratum PLB 407 (1997) 453** | `arXiv:hep-ph/9611460`, DOI 10.1016/S0370-2693(97)00250-5 | ✓ | **new** `hep-ph_9611460.pdf` | the shadowing camp: b₂ rising at small x, A₂ ≈ 1 %, **two orders of magnitude above impulse approximation**, and an explicit claim that the Close–Kumano sum rule is *broken* | theory-input |
| **bora-jaffe1998** | **Bora, Jaffe, *The double-scattering contribution to b₁(x, Q²) in the deuteron*, PRD 57 (1998) 6906** | `arXiv:hep-ph/9711323`, DOI 10.1103/PhysRevD.57.6906 | ✓ | **new** `hep-ph_9711323.pdf` | the VMD/double-scattering camp, with the opposite small-x verdict to Nikolaev–Schäfer: b₁ → 0 as x → 0 by rotational symmetry | theory-input |
| **edelmann-piller-weise1997** | **Edelmann, Piller, Weise, *Polarized deuteron structure functions at small x*, Z. Phys. A 357 (1997) 129** | `arXiv:nucl-th/9701026`, DOI 10.1007/s002180050226 | ✓ | **new** `nucl-th_9701026.pdf` | third shadowing calculation: b₁ "surprisingly large at x < 0.1", dominated by coherent double scattering | theory-input |
| **edelmann-piller-weise1998** | **Edelmann, Piller, Weise, *Deuteron spin structure functions at small Bjorken x*, PRC 57 (1998) 3392** | `arXiv:hep-ph/9709455`, DOI 10.1103/PhysRevC.57.3392 | ✓ | **new** `hep-ph_9709455.pdf` | the long version: spin-extended Glauber–Gribov multiple scattering | theory-input |
| **umnikov1997** | **Umnikov, *Relativistic calculation of structure functions b₁,₂(x) of the deuteron*, PLB 391 (1997) 177** | `arXiv:hep-ph/9605291`, DOI 10.1016/S0370-2693(96)01440-2 | ✓ | **new** `hep-ph_9605291.pdf` | **the validity floor on `b1_convolution`**: non-relativistic convolution gets small x wrong and violates the exact sum rules | theory-input |
| **kumano-kuroki2026** | **Kumano, Kuroki, *Tensor-polarized parton distribution functions of the deuteron by a convolution model* (2026)** | `arXiv:2607.09237` | ✓ | **new** `2607.09237.pdf` | **an independent second convolution calculation at exactly Q² = 2.5 GeV²** — the same point as `tables::kB1CdksQ2p5` | theory-input |

### 1d. Tensor-polarized PDFs, sum rules, the modern spin-1 basis

| key | reference | id | OA | local | serves | tier |
|---|---|---|---|---|---|---|
| **kumano2010-tensor-pdf-fit** | **S. Kumano, *Tensor-polarized quark and antiquark distribution functions in a spin-one hadron*, PRD 82 (2010) 017501** | `arXiv:1005.4524`, DOI 10.1103/PhysRevD.82.017501 | ✓ | **new** `1005.4524.pdf` | **a TABLE, not a figure**: a closed-form δ_T w(x) with fitted parameters at Q² = 2.5 GeV². A data-driven third `TensorSF` camp | theory-input |
| kumano-song2016-pd-dy | Kumano, Song, *Theoretical estimate on tensor-polarization asymmetry in proton-deuteron Drell-Yan*, PRD **94** (2016) 054022 | `arXiv:1606.03149` | ✓ | **new** `1606.03149.pdf` | the Fermilab route to δ_T q̄, using the 1005.4524 PDFs | theory-input |
| kumano-song2021-tmd-twist4 | Kumano, Song, *TMD parton distribution functions up to twist 4 for spin-1 hadrons*, PRD **103** (2021) 014025 | `arXiv:2011.08583` | ✓ | **new** `2011.08583.pdf` | the complete spin-1 correlator decomposition: **40 TMDs** in a tensor-polarized spin-1 hadron over twists 2–4 | theory-input |
| **kumano-song2021-sum-rule** | **Kumano, Song, *Twist-2 relation and sum rule for tensor-polarized PDFs of spin-1 hadrons*, JHEP 09 (2021) 141** | `arXiv:2106.15849`, DOI 10.1007/JHEP09(2021)141 | ✓ | **new** `2106.15849.pdf` | **a second free, closed-form sum rule** beside Close–Kumano: ∫dx f₂LT = 0 with f₂LT = (2/3)f_LT − f₁LL, plus a WW-type twist-2 relation | T5 |
| kumano2014-tensor-2020s-jpcs | S. Kumano, *Tensor-polarized structure functions: tensor structure of deuteron in 2020's*, J. Phys. Conf. Ser. **543** (2014) 012001 | `arXiv:1407.3852`, DOI 10.1088/1742-6596/543/1/012001 | ✓ | **new** `1407.3852.pdf` | **closed-form projection operators extracting b₁–b₄ from the hadron tensor W_μν** — a free, independent check on `InclusiveKernel`'s transcription of the Cosyn *et al.* decomposition | theory-input |
| kumano2024-spin1-review | S. Kumano, *Parton distribution functions and fragmentation functions of spin-1 hadrons*, EPJ A **60** (2024) 205 | `arXiv:2406.01180`, DOI 10.1140/epja/s10050-024-01411-6 | ✓ | **new** `2406.01180.pdf` | the single best entry point to the whole spin-1 sector (20 figures); the review to cite once instead of six primaries | theory-input |
| zhao2025-sidis-tensor | Zhao, Bacchetta, Kumano, Liu, Zhou, *SIDIS off a tensor-polarized spin-1 target*, JHEP **12** (2025) 067 | `arXiv:2508.06134`, DOI 10.1007/JHEP12(2025)067 | ✓ | **new** `2508.06134.pdf` | 23 structure functions, 21 non-vanishing at tree level to twist 3 — the SIDIS counterpart of the inclusive decomposition the tree implements | theory-input |

### 1e. The double-helicity-flip gluon Δ ("exotic glue")

| key | reference | id | OA | local | serves | tier |
|---|---|---|---|---|---|---|
| sather-schmidt1990 | Sather, Schmidt, PRD **42** (1990) 1424 | no arXiv | ✗ | — (in corpus, `file: null`) | **the actual source of `C_BAG` = −0.012** | theory-input |
| detmold-shanahan2016-lattice-phi | Detmold, Shanahan, *Gluonic transversity from lattice QCD*, PRD **94** (2016) 014507; Erratum PRD **95** (2017) 079902 | `arXiv:1606.04505` | ✓ | **new** `1606.04505.pdf` | A₂ ≈ 0.23(2)(5) in the **φ meson**; gluonic Soffer bound saturated at ≈ 80 % (PDF-read) | theory-input |
| **winter2017-lattice-light-nuclei** | **Winter, Detmold, Gambhir, Orginos, Savage, Shanahan, Wagman (NPLQCD), *First lattice QCD study of the gluonic structure of light nuclei*, PRD 96 (2017) 094512** | `arXiv:1709.00395`, DOI 10.1103/PhysRevD.96.094512 | ✓ | **new** `1709.00395.pdf` | **the first exotic-glue number in a NUCLEUS**: a₂^(d) = −0.010(3) bare at m_π ≈ 806 MeV against ⟨x⟩_g^(d) = 0.51(5), and c₂^(d)/b₂^(d) ≲ 1/20 for the gluonic-b₁ combination (PDF-read) | theory-input |
| kumano-song2020-gluon-transversity-dy | Kumano, Song, *Gluon transversity in polarized pd Drell-Yan*, PRD **101** (2020) 054011 | `arXiv:1910.12523` | ✓ | **new** `1910.12523.pdf` | the "Drell–Yan estimates" half of `coherent.hpp`'s bound comment | theory-input |
| kumano-song2020-deuteron-polarizations-dy | Kumano, Song, *Deuteron polarizations in the pd Drell-Yan process for finding the gluon transversity*, PRD **101** (2020) 094013 | `arXiv:2003.06623` | ✓ | **new** `2003.06623.pdf` | the same, re-expressed in conventional deuteron polarizations (measurable ones) | theory-input |
| **cotogno-vandaal-mulders2017-positivity** | **Cotogno, van Daal, Mulders, *Positivity bounds on gluon TMDs for hadrons of spin ≤ 1*, JHEP 11 (2017) 185** | `arXiv:1709.07827`, DOI 10.1007/JHEP11(2017)185 | ✓ | **new** `1709.07827.pdf` | **a free positivity gate on `toy_delta_gluon`** — leading-twist gluon-distribution bounds including tensor polarization, with a small-x limit | T5 |

### 1f. Spin-3/2 — the ⁷Li rank-2 sector

| key | reference | id | OA | local | serves | tier |
|---|---|---|---|---|---|---|
| fu-sun-dong2022-spin32-gpd | Fu, Sun, Dong, *Generalized parton distributions in spin-3/2 particles*, PRD **106** (2022) 116012 | `arXiv:2209.12161` | ✓ | **new** `2209.12161.pdf` | the J = 3/2 forward-limit structure functions; the b₁ definition to normalise ⁷Li against. **Their "g₂" is rank-3 and collides with the twist-3 nucleon g₂** | theory-input |
| fu-dong-kumano2024-spin32-transversity-gpd | Fu, Dong, Kumano, *Transversity GPDs in spin-3/2 particles*, PRD **109** (2024) 096006 | `arXiv:2402.11561` | ✓ | **new** `2402.11561.pdf` | 16 transversity GPDs per parton; the gluon-transversity analogue of Δ for J = 3/2 | theory-input |
| **fu-dong-kumano-xie2026-positivity** | **Fu, Dong, Kumano, Xie, *Generalizing the Soffer bound: positivity constraints on PDFs of spin-3/2 particles*, PRD 113 (2026) L111901** | `arXiv:2602.11587`, DOI 10.1103/qq1h-snhk | ✓ | **new** `2602.11587.pdf` | **the complete set of spin-3/2 positivity bounds, first time.** The cheapest new external gate available to this project | T5 |

### 1g. Spin-1 GPDs, tagging and FSI

| key | reference | id | OA | local | serves | tier |
|---|---|---|---|---|---|---|
| **berger-cano-diehl-pire2001** | **Berger, Cano, Diehl, Pire, *Generalized parton distributions in the deuteron*, PRL 87 (2001) 142302** | `arXiv:hep-ph/0106192`, DOI 10.1103/PhysRevLett.87.142302 | ✓ | **new** `hep-ph_0106192.pdf` | the spin-1 GPD basis (5 quark + 9 gluon), incl. the gluon-transversity GPD; the exclusive counterpart of the coherent ⁶Li channel | theory-input |
| **cosyn-weiss2020-polarized-tagged** | **Cosyn, Weiss, *Polarized e–d DIS with spectator nucleon tagging*, PRC 102 (2020) 065204** | `arXiv:2006.03033`, DOI 10.1103/PhysRevC.102.065204 | ✓ | **new** `2006.03033.pdf` | 52 pp, 19 figs: the LF polarized spectral function. **The one cheap unexploited win — it gates `TaggedModel`'s VECTOR sector, which E-1 does not touch** | T2 |
| **cosyn-sargsian2014-tensor-fsi** | **Cosyn, Sargsian, *FSI in DIS from a tensor polarized deuteron target*, J. Phys. Conf. Ser. 543 (2014) 012006** | `arXiv:1407.1653`, DOI 10.1088/1742-6596/543/1/012006 | ✓ | **new** `1407.1653.pdf` | **the only calculation of A_zz with any non-Born effect.** Inclusive: sizeable FSI for x > 0.2 but A_zz stays small; tagged: FSI largest at p_spec ≈ 300 MeV, forward angles | theory-input |

### 1h. HERMES Table II, transcribed

**(PDF-read)** from `hep-ex_0506018.pdf`. This is the whole of the world's
tensor-DIS data. Both A_zz^d and b₁^d are in units of 10⁻², statistical then
systematic uncertainty. Kinematic selection for the sample: 0.1 < Q² < 20 GeV²,
W > 1.8 GeV, 0.002 < x < 0.85, 0.1 < y < 0.91.

| ⟨x⟩ | ⟨Q²⟩ [GeV²] | A_zz^d | δA_zz^stat | δA_zz^sys | b₁^d | δb₁^stat | δb₁^sys |
|---|---|---|---|---|---|---|---|
| 0.012 | 0.51 | −1.06 | 0.52 | 0.26 | 11.20 | 5.51 | 2.77 |
| 0.032 | 1.06 | −1.07 | 0.49 | 0.36 | 5.50 | 2.53 | 1.84 |
| 0.063 | 1.65 | −1.32 | 0.38 | 0.21 | 3.82 | 1.11 | 0.60 |
| 0.128 | 2.33 | −0.19 | 0.34 | 0.29 | 0.29 | 0.53 | 0.44 |
| 0.248 | 3.11 | −0.39 | 0.39 | 0.32 | 0.29 | 0.28 | 0.24 |
| 0.452 | 4.69 | +1.57 | 0.68 | 0.13 | −0.38 | 0.16 | 0.03 |

Three points the tree should carry with these numbers. (i) The **sign change**
between ⟨x⟩ = 0.248 and 0.452 in both columns is the feature any b₁ camp must
reproduce, and it is a 2-sigma-scale statement, not a measurement of the zero
crossing. (ii) The **definition** is A_zz^d = (2σ¹ − 2σ⁰)/(3 σ_U P_zz^eff)
(Eq. 2), with σ_U = (2σ¹ + σ⁰)/3 the polarization-averaged yield — so the sign
convention of `TENSOR_LL_SIGN` is fixed against *this* equation and not against
a generic one. (iii) Radiative background is ≈ 50 % of the statistics in the
lowest-x bin, and the paper's own residual RC systematic on A_zz is ≈ 2 × 10⁻³
at low x against |A_zz| ≤ 0.02 — so the first two rows carry a correction
comparable to their own central values. **The six rows above are BENCHMARK_PLAN
§4 item 3: transcribe them as an assertion and the only measurement of the
target observable stops being a comment.**

---

## 2. Foundations — what the tree transcribes, and the two gaps in it

### 2.1 Hoodbhoy–Jaffe–Manohar (NPB 312) and Jaffe–Manohar (NPB 321)

Both are in the corpus as one folded entry with `file: null`, and both are
genuinely closed: no arXiv, no preprint scan on INSPIRE (verified this session —
INSPIRE recids 262935 and 266865, 236 and 87 citations respectively, `arxiv:
None` on both). NPB 312 defines b₁…b₄ and Δ and gives the (1, −2, 1) helicity
pattern the tree's rank-2 geometry implements; NPB 321 generalises to arbitrary
J and is the formal home of the ⁷Li J = 3/2 basis. Nothing found in this pass
changes their status; §7 lists the DOI pages.

### 2.2 Jaffe–Manohar, *Nuclear gluonometry*, PLB 223 (1989) 218 — **a corpus gap**

Verified: INSPIRE recid 277858, DOI 10.1016/0370-2693(89)90242-6, 107 citations,
no arXiv. This is the paper that *names* Δ an "exotic glue" observable and
argues it has no nucleonic counterpart — the conceptual source of
`toy_delta_gluon` and of the coherent cos 2φ term. It is **not** in
`refs_dict.json` and it is **not** the same paper as either 1989 NPB entry, so
the corpus's folded `hoodbhoy-jaffe-manohar1989-and-jaffe-manohar1989` entry
does not cover it. The JLab Letter of Intent (§1b) quotes it directly —
"first identified by Jaffe and Manohar in 1989 as *a clear signature for exotic
gluonic components in the target*" — so `arXiv:1803.11206`, which is free and
now on disk, is a usable secondary citation while the primary is unobtainable.
**Priority: medium.** Cheap to fix at a library; nothing numeric depends on it.

### 2.3 Close–Kumano, PRD 42 (1990) 2377 — **a corpus gap, and the one that matters**

Verified: INSPIRE recid 296530, DOI 10.1103/PhysRevD.42.2377, 105 citations, no
arXiv. The abstract, read in full this session, states the result exactly:
∫dx b₁(x) = −(5/3) lim_{t→0} (t/4M²) F_Q(t) = 0 for isoscalar targets **if the
quark/antiquark sea is unpolarized**, and shows how a polarized sea modifies it.

This is the source of `close_kumano_integral` and of in-tree check E-4, and the
corpus does not hold it. That is the sharpest citation gap in this domain: a
gate whose defining paper is uncited-in-corpus. Two things follow:

- The sum rule is **conditional**, and the condition is precisely what
  Kumano's own 2010 fit finds to be violated (§4.1): a finite tensor-polarized
  antiquark distribution gives ∫b₁ = 0.0058 rather than 0. `close_kumano_integral`
  is documented as "reported, not enforced", which is the correct behaviour —
  but the docstring should say *why*, and the reason is in this paper plus
  1005.4524.
- Nikolaev–Schäfer (§3.1) claim their two shadowing mechanisms **break** the
  sum rule outright. So the rule is one of three positions in the literature,
  not a theorem. **Priority: high** — obtain it, or cite it through
  `arXiv:1005.4524` §I, which restates it with the antiquark correction.

### 2.4 Khan–Hoodbhoy, PRC 44 (1991) 1219 — paywalled

Verified: INSPIRE recid 324952, DOI 10.1103/PhysRevC.44.1219, 56 citations, no
arXiv, no preprint scan. Abstract (read in full): a convolution model with
relativistic and binding corrections, giving spin and spin-averaged twist-2
deuteron structure functions as "a simple parametrization … in terms of a few
deuteron wave-function parameters and the free nucleon structure functions".
That form — b₁ as an explicit function of wave-function moments — is exactly the
shape one would refit for an α–d system, which is why it is worth a library
copy even though `b1_convolution` implements CDKS instead. Their 1993 sequel
(PLB 298, 181) is also paywalled and adds the small-x double-scattering piece.
**Priority: low-medium.**

---

## 3. The b₁ᵈ camps — and the small-x disagreement the tree does not currently represent

`B1Mode`/`TensorSF` offers two camps today: `MillerB1` (pion + hidden colour,
tuned to reproduce HERMES) and `CdksB1` (standard convolution). Both are
figure-only and both are already digitized. **Five further calculations exist,
all now on disk, and three of them disagree with each other about the sign and
size of b₁ at small x** — which is the region where LiPolGen's own EIC reach
lives and where the two in-tree camps differ by two orders of magnitude.

### 3.1 Nikolaev & Schäfer, PLB 398 (1997) 245 (+ Erratum PLB 407, 453)

`arXiv:hep-ph/9611460`, verified against both the arXiv API and INSPIRE
(recid 426591, 36 citations, two DOIs for article + erratum). Two mechanisms —
diffractive nuclear shadowing sensitive to the alignment of the nucleons, and
the deuteron's nuclear pion excess — give a b₂ that **rises towards small x**
and a tensor asymmetry A₂ = b₂/F₂ᵈ of about **one percent**, which the abstract
says "by almost two orders of magnitude exceeds the effect evaluated earlier in
the impulse approximation". The abstract also states outright that **both
mechanisms break the Close–Kumano sum rule**.

*What it serves*: the "HERMES-like" camp's *physical* justification. Miller
reproduces HERMES with hidden colour; this paper reproduces a HERMES-sized
effect with shadowing instead, from 1997 — i.e. **before** the measurement.
It is the strongest available argument that a large small-x b₁ is not exotic,
and it belongs in any discussion of `toy_b1`'s magnitude. It also supplies the
honest caveat on `close_kumano_integral`.

### 3.2 Bora & Jaffe, PRD 57 (1998) 6906

`arXiv:hep-ph/9711323`, verified (INSPIRE recid 451055, 38 citations,
MIT-CTP-2692). A VMD model of shadowing by double scattering of vector mesons,
plus the D-state admixture. Its conclusion is the **opposite** of
Nikolaev–Schäfer's at the smallest x: "the restoration of rotational symmetry
for small x (≤ 10⁻³) requires that b₁⁽²⁾ approach zero as x → 0 in this model",
while remaining "significant at large Bjorken x". Note the v2 comment records a
correction in the other direction — "discussion of impact parameter dependence
of shadowing has been corrected, confirming a more significant enhancement of
b₁ at small-x".

*What it serves*: a **sign/limit constraint** on any b₁ extrapolation below the
digitized floor. `b1_convolution` freezes its table at x = 0.010 and `toy_b1`
freezes at 0.0647; neither behaviour is derived, and this paper plus §3.1 and
§3.3 bracket what the literature actually permits there. Cheap, real, and
currently absent from the tree: a documented small-x extrapolation band with
three citations behind it rather than a flat freeze.

### 3.3 Edelmann, Piller & Weise — Z. Phys. A 357 (1997) 129 and PRC 57 (1998) 3392

Both verified (`arXiv:nucl-th/9701026` / recid 439601, 29 citations; and
`arXiv:hep-ph/9709455` / recid 448857, 25 citations). The pair is one
calculation published twice, short then long: shadowing corrections to g₁ᵈ and
b₁ from coherent double scattering, with the Glauber–Gribov multiple-scattering
formalism **extended to include spin degrees of freedom** — which is the
technical piece Cosyn–Weiss II §VI C lists as an open question for the tagged
case. Their result: b₁ is "surprisingly large at x < 0.1" and coherent double
scattering dominates it. A third, independent shadowing camp agreeing with
Nikolaev–Schäfer on the direction.

*What it serves*: the third leg of the small-x band; and, methodologically, a
worked precedent for spin-dependent Glauber, which is what LiPolGen's Glauber
FSI module would need if it ever carried spin.

### 3.4 Umnikov, PLB 391 (1997) 177 — the validity floor

`arXiv:hep-ph/9605291`, verified. A covariant Bethe–Salpeter calculation of
b₁,₂ᵈ whose stated purpose is negative: "It is shown that usual nonrelativistic
convolution model result in incorrect behavior of this structure functions at
small X and **violates the exact sum rules**."

*What it serves*: `b1_convolution` is exactly a non-relativistic convolution,
and it is currently extrapolated below its table with no stated validity range.
This paper is the citation for putting one in — an x ≳ 0.1 band, or at minimum
a documented warning — and for explaining why `close_kumano_integral` is
reported rather than enforced. **Effort: 0.5 d for a docstring and an x-range
guard. Priority: medium**, and it is the single cheapest honesty upgrade in
this file.

### 3.5 Kumano & Kuroki, `arXiv:2607.09237` (2026) — the second convolution

Verified against the arXiv API: 7 pages, 5 figures, submitted 10 Jul 2026, no
journal reference yet. Tensor-polarized PDFs for the deuteron from a
convolution of nucleon unpolarized PDFs with the tensor-polarized nucleon
momentum distribution, computed **at Q² = 2.5 GeV² specifically in order to
compare with the HERMES-fitted PDFs**. Their headline is a disagreement: "The
obtained distributions are very different from the ones determined from the
HERMES data, which indicates further studies are needed … possibly by
considering a new mechanism beyond the simple bound system of a proton and a
neutron." They also convert δ_T q → the Trento-convention f₁LL^q and estimate
the twist-3 f_LT by a WW-like relation.

*What it serves*: **E-2's single-camp shape test becomes a two-calculation
spread.** Q² = 2.5 GeV² is the same point as `tables::kB1CdksQ2p5`, so the
comparison is direct with no evolution. And their stated conclusion is the same
tension E-2 already reports between the convolution camp and the HERMES-anchored
camp — an independent group reaching the project's own conclusion is worth more
than another curve. Figure-only (no tables, checked); digitizable with the
sibling repo's `tools/digitize_figure.py`. **Effort 0.5–1 d. Priority: high.**

---

## 4. Tensor-polarized PDFs and the modern spin-1 basis

### 4.1 Kumano, PRD 82 (2010) 017501 — **the only b₁ source in this domain that publishes numbers**

`arXiv:1005.4524`, verified against the arXiv API and read from the PDF. This
is the one entry in §1c–1d that is **TABLE**, not **FIGURE**, and it is the
reason it ranks high despite being four pages.

**(PDF-read)** The parametrization is a closed form:

> δ_T q_v^D(x) = δ_T w(x) · [u_v(x) + d_v(x)]/2,
> δ_T q̄^D(x) = α_q̄ · δ_T w(x) · [2ū + 2d̄ + s + s̄]/6, with
> **δ_T w(x) = a x^b (1 − x)^c (x₀ − x)** (Eq. 7),

and Table I gives the fitted parameters **at Q² = 2.5 GeV²** — the same Q² as
`tables::kB1CdksQ2p5` and as Kumano–Kuroki:

| analysis | χ²/d.o.f. | a | α_q̄ | b | c | x₀ |
|---|---|---|---|---|---|---|
| Set 1 (no tensor sea) | 2.83 | 0.378 ± 0.212 | 0.0 (fixed) | 0.706 ± 0.324 | 1.0 (fixed) | 0.229 |
| Set 2 (tensor sea free) | 1.57 | 0.221 ± 0.174 | 3.20 ± 2.75 | 0.648 ± 0.342 | 1.0 (fixed) | 0.221 |

Two facts from the body matter to this project. First, **Set 1 does not fit**
(χ²/dof = 2.83) and Set 2 does (1.57): HERMES requires a finite tensor-polarized
*antiquark* distribution. Second, the consequence for the sum rule is quantified
— with the Set-2 sea, ∫dx b₁(x) = **0.0058**, not zero (Eq. 9), against HERMES's
own ∫₀.₀₂^0.85 b₁ dx = [0.35 ± 0.10(stat) ± 0.18(sys)] × 10⁻².

*What it serves*: a third `TensorSF` implementation, `KumanoFitB1`, that is
**data-anchored rather than model-anchored** and needs no digitization — three
parameters, a published functional form, a stated Q², and a stated χ². Beside
`MillerB1` and `CdksB1` it turns the b₁ registry from "two model camps" into
"two models and the fit to the only measurement". It also gives
`close_kumano_integral` its documented non-zero expectation.
**Effort: 1–2 d (it needs an unpolarized PDF set for u_v + d_v, which the tree
already has via LHAPDF). Priority: high.**

### 4.2 Kumano & Song, JHEP 09 (2021) 141 — a second free sum rule

`arXiv:2106.15849`, verified. Abstract read in full: an analogue of the
Wandzura–Wilczek relation and of the Burkhardt–Cottingham sum rule for the
tensor sector. The twist-2 part of the twist-3 f_LT is an integral of the
twist-2 f₁LL (equivalently b₁), and **∫dx f₂LT = 0 with f₂LT ≡ (2/3) f_LT −
f₁LL**. If the parton-model sum rule for f₁LL is imposed (vanishing tensor sea),
a further sum rule for f_LT itself follows.

*What it serves*: **T5.** The tree currently checks exactly one tensor sum rule
(Close–Kumano, reported not enforced). This adds a second, in a different
function, that is free to evaluate from quantities the tree already computes.
Note the conditional structure mirrors §2.3 and §4.1: the stronger f_LT rule
holds only if the tensor sea vanishes, which the fit says it does not — so the
useful gate is the unconditional ∫f₂LT = 0.

### 4.3 Kumano & Song, PRD 103 (2021) 014025 — the full TMD basis

`arXiv:2011.08583`, verified. The complete Lorentz decomposition of the quark
correlator for spin-1 up to twist 4, with Hermiticity and parity imposed but
**not** time reversal (so T-odd functions survive). Headline counts from the
abstract: **40 TMDs in a tensor-polarized spin-1 hadron across twists 2, 3, 4**,
of which **30 structure functions in twists 3 and 4 are new in that work**; and
some twist-2 expressions are *corrected* relative to earlier derivations because
of previously-missing terms involving the light-cone vector n.

*What it serves*: the authority for any claim this project makes about how many
tensor structure functions there are, and a warning that pre-2021 twist-2 spin-1
expressions may be incomplete. Worth checking `InclusiveKernel`'s transcription
of the Cosyn *et al.* decomposition against it once.

### 4.4 Kumano, EPJ A 60 (2024) 205 — the review

`arXiv:2406.01180`, verified (12 pages, 20 figures; EPJ A topical issue on
tensor spin observables). Current theoretical status of spin-1 structure
functions, TMDs and fragmentation functions, written against the late-2020s
experimental programme (JLab, Fermilab, NICA, EICs). **This is the one citation
to use when a paper needs a single pointer to the tensor sector rather than six
primaries** — and per the project's own rule (prefer the paper a number comes
from), it is a *navigation* reference, not a source of numbers.

### 4.5 Zhao, Bacchetta, Kumano, Liu, Zhou, JHEP 12 (2025) 067 — SIDIS

`arXiv:2508.06134`, verified. The complete differential cross section for SIDIS
off a tensor-polarized spin-1 target: **23 structure functions, 21 non-vanishing
at tree level** through subleading twist, as convolutions of TMD PDFs and FFs.

*What it serves*: not the inclusive channel LiPolGen generates, but it is the
theory behind the CLAS12 proposal (§1b) that will produce the next tensor data,
and it is the natural cross-check on Cosyn–Weiss I Appendix D's SIDIS counting
(which gives 41 structure functions for j = 1 — a *different* count because it
includes beam polarization and all twists; the two are not in conflict but the
difference must be explained if both are cited).

### 4.6 Bacchetta & Mulders, PRD 62 (2000) 114004

`arXiv:hep-ph/0007120`, verified. The original leading-twist spin-1 TMD
analysis, including T-odd functions — the paper that named f₁LL and h⊥₁LL. It is
the convention Kumano–Kuroki convert *into* when they quote Trento-convention
PDFs, so holding it is what makes that conversion checkable rather than
asserted. *Serves*: the naming of the tensor TMDs; relevant if LiPolGen ever
emits a φ_h-dependent observable.

---

## 5. The double-helicity-flip gluon Δ — and the finding that changes the story

This is where the pass produced the most: `../benchmarking/04_theory.md` §4
concludes, correctly on the evidence it had, that the only defensible Δ gate is
the gluonic Soffer bound from a **φ meson**, and that comparing that number to
`C_BAG` would be a wrong-observable benchmark. **That conclusion can now be
strengthened rather than merely stated.**

### 5.1 Detmold & Shanahan, PRD 94 (2016) 014507 (Erratum PRD 95 (2017) 079902)

`arXiv:1606.04505`, verified against the arXiv API (journal-ref field carries
the erratum; INSPIRE recid carries both DOIs, 27 citations) and **(PDF-read)**:
"we find A₂ ∼ 0.23(2)(5), where the first un[certainty]…" at §III, and the
abstract's gluonic-Soffer statement, "finding it to be saturated at the level of
80 %". The target is the φ meson at m_π = 450 MeV.

### 5.2 Winter *et al.* (NPLQCD), PRD 96 (2017) 094512 — **the new one**

`arXiv:1709.00395`, verified against the arXiv API and read from the PDF. This
is the paper `04_theory.md` §4 did not have, and it changes the Δ discussion in
three ways.

**(PDF-read), §V and §IV:**

- The deuteron gluon-transversity reduced matrix element is
  **a₂^(d) = −0.010(3)** (unrenormalised, m_π ≈ 806 MeV), "statistically
  non-zero to three standard deviations … the first evidence for non-nucleonic
  gluon contributions to nuclear structure". No signal is resolved at
  m_π ≈ 450 MeV.
- The comparison quantity in the same calculation is the bare gluon momentum
  fraction **⟨x⟩_g^(d) = b₂^(d) = 0.51(5)**, so the transversity-to-unpolarised
  ratio in the deuteron is ≈ 2 % — which the paper states is "approximately an
  order of magnitude smaller than seen in previous studies of the φ meson",
  and attributes to large-N_c scaling plus the deuteron's loose binding.
- For the **gluonic analogue of b₁** — c₂^(d), the first moment of the
  difference between the unpolarised gluon PDF in the j_z = ±1 and j_z = 0
  deuteron states — no signal is resolved at either pion mass, and the paper
  quotes the resulting bound **c₂^(d)/b₂^(d) ≲ 1/20**.
- The unpolarised gluon EMC effect in A ≤ 3 is constrained to < ~10 %.

*What it serves, concretely.*

1. **It kills the transfer of the φ number to a nucleus, quantitatively.** The
   φ saturates the gluonic Soffer bound at 80–100 %; the deuteron's ratio is
   ~2 %, an order of magnitude down. Any ⁶Li Δ scenario scaled from
   Detmold–Shanahan is therefore high by roughly 10×, and this is now a
   citable statement rather than a caution. `toy_delta_gluon`'s discovery
   scenarios (`scale` = 1e-3, 3e-3, 1e-2 as peak Δ/F₁) should be read against
   it.
2. **It supplies a nuclear bound for the b₁-like gluon moment**, c₂/b₂ ≲ 1/20,
   which is the first external number of any kind for a gluonic tensor
   observable in a nucleus.
3. It is the calculation the JLab LoI (§5.3) cites as having "prompted renewed
   interest", so the two belong together in any writeup.

**Caveats the paper states about itself, and which must travel with the number:
** unrenormalised matrix elements; m_π ≈ 806 MeV, far from physical; a short fit
window with "contamination from excited states that is not fully quantified";
and the ratio argument assumes the two operators' multiplicative renormalisation
roughly cancels. It is an indication, not a measurement.
**Priority: high** — 1 d to write the Δ section of the docs against it;
it does not require any code change to be worth having.

### 5.3 Maxwell *et al.*, LoI to JLab PAC 44, `arXiv:1803.11206`

Verified; the document itself is dated 6 June 2016 (PDF-read) and was posted to
arXiv in March 2018. Renewed intent to search for a non-zero Δ(x, Q²) with an
**unpolarized electron beam and a transversely polarized spin-1 nuclear
target**, inclusive DIS below x = 0.3, via single-spin tensor asymmetries, at
12 GeV with the JLab/UVa solid polarized target. Authors include Detmold and
Jaffe — i.e. the lattice and the 1989 theory, on the same document.

*What it serves*: the **experimental** side of `RcScope::TensorAll` and
`toy_delta_gluon`. It is the only programme anywhere aiming at Δ, it states the
observable in the same geometry LiPolGen's coherent cos 2φ term uses, and it is
a free, quotable secondary source for the `Nuclear gluonometry` definition the
corpus lacks (§2.2).

### 5.4 Kumano & Song, PRD 101 (2020) 054011 and 094013

Both verified. 054011 (`arXiv:1910.12523`, 21 pp, 16 figures) proposes finding
the gluon transversity in polarized proton–deuteron Drell–Yan, noting it cannot
exist for spin-1/2 by helicity conservation; 094013 (`arXiv:2003.06623`) recasts
the same cross sections in terms of **conventional deuteron polarizations**, so
that the observable is the cross-section difference between deuteron spin
polarizations along two transverse axes — the experimentally realisable version.

*What it serves*: the "Drell–Yan estimates" half of the bound comment in
`coherent.hpp`. Holding both means that comment can now name a paper and an
equation instead of a genre.

### 5.5 Cotogno, van Daal & Mulders, JHEP 11 (2017) 185 — a free positivity gate

`arXiv:1709.07827`, verified. The light-front gluon–gluon correlator extended
from unpolarized and vector-polarized targets **to tensor-polarized targets**,
with process-dependent gauge links, and — the operative part — "positivity
bounds for combinations of leading-twist gluon distributions that may be used to
estimate their maximal contribution to observables", plus the small-x limit for
dipole-type gauge links via a single Wilson loop.

*What it serves*: **the only model-free constraint on `toy_delta_gluon` that
does not go through a lattice calculation at an unphysical pion mass.** Nothing
in the tree currently bounds Δ at all; the scenarios are asserted. Reading the
bounds out of the body (they are not in the abstract) and adding a
`test_sf.cpp` assertion that every shipped Δ scenario satisfies them is, with
§6.3, the cheapest new external gate in this domain. **Effort 1–2 d.
Priority: high.**

---

## 6. Spin-3/2 — the ⁷Li rank-2 sector, which is identically zero today

### 6.1 Fu, Sun & Dong, PRD 106 (2022) 116012

`arXiv:2209.12161`, verified. Eight unpolarized and eight polarized GPDs for
spin-3/2, defined for the first time; the forward limit gives the structure
functions and PDFs. Already read and summarised in
`../open_items/physics_literature.md`, which records the four leading-twist
J = 3/2 DIS structure functions and the **naming collision**: their "g₂" is the
rank-3 octupole partner of g₁, not the twist-3 nucleon g₂, and must be renamed
(`g1_rank3`) before any use in code.

### 6.2 Fu, Dong & Kumano, PRD 109 (2024) 096006

`arXiv:2402.11561`, verified. Quark and gluon **transversity** GPDs for
spin-3/2 in light-cone gauge: 16 independent components per parton, 16
helicity-flip amplitudes, and the result that the odd-skewness transversity GPDs
vanish in the forward limit. *Serves*: the J = 3/2 analogue of Δ — i.e. whether
a ⁷Li "exotic glue" term exists at all and in what forward-limit combination.
The tree's ⁷Li sector does not need this yet; the moment it does, this is the
definition.

### 6.3 Fu, Dong, Kumano & Xie, PRD 113 (2026) L111901 — **the gate**

`arXiv:2602.11587`, verified against INSPIRE (recid 3118921, DOI
10.1103/qq1h-snhk) and read from the PDF. "We derive the complete set of
positivity bounds for the leading-twist parton distribution functions of a
spin-3/2 hadron for the first time", generalizing the Soffer bound
|h₁| ≤ ½(f₁ + g₁) to quark *and gluon* distributions in higher-spin systems.
Method **(PDF-read)**: antiparton–hadron scattering amplitudes expressed via the
PDFs and the spin density matrix; positive-definiteness of the amplitude matrix
yields the inequalities; Cauchy–Schwarz relates PDFs to GPDs. The paper is
explicit that the spin-1 constraints are the intermediate case it generalizes.

*What it serves*: LiPolGen's ⁷Li rank-2 block and **every invented scenario
shape** — `toy_delta_gluon`, and the b₃ = 0.05·f₁ / b₄ = −0.02·f₁ placeholders.
Today `tests/test_pipeline.cpp` D1 pins the ⁷Li rank-2 sector at identically
zero, so the bounds pass trivially; they become a real gate the moment that
sector is filled, and they cost nothing to add now as a guard that will catch a
future mistake. Combined with §5.5 (spin ≤ 1) this project can cover both
isotopes with published inequalities and no model at all.
**Priority: high** — reaffirming `04_theory.md` T-13, now with the published
journal reference rather than a bare preprint id.

---

## 7. What could not be obtained, and exactly where to get it

No paywall was approached. Four items in this domain have no free copy, all
pre-arXiv, all verified to exist with the citation given:

| reference | why the project wants it | DOI page to use |
|---|---|---|
| Hoodbhoy, Jaffe, Manohar, NPB **312** (1989) 571 | the b₁…b₄ and Δ definitions the tree transcribes | https://doi.org/10.1016/0550-3213(89)90572-5 |
| Jaffe, Manohar, NPB **321** (1989) 343 | the arbitrary-J basis for the ⁷Li rank-2 block | https://doi.org/10.1016/0550-3213(89)90347-7 |
| Jaffe, Manohar, PLB **223** (1989) 218, *Nuclear gluonometry* | the definition of Δ as exotic glue — **not in the corpus at all** | https://doi.org/10.1016/0370-2693(89)90242-6 |
| Close, Kumano, PRD **42** (1990) 2377 | the ∫b₁ dx = 0 sum rule behind `close_kumano_integral` (E-4) | https://doi.org/10.1103/PhysRevD.42.2377 |
| Khan, Hoodbhoy, PRC **44** (1991) 1219 | convolution b₁ parametrized in wave-function moments | https://doi.org/10.1103/PhysRevC.44.1219 |
| Khan, Hoodbhoy, PLB **298** (1993) 181 | its small-x double-scattering sequel | https://doi.org/10.1016/0370-2693(93)91727-5 |
| Sather, Schmidt, PRD **42** (1990) 1424 | the source of `C_BAG` = −0.012 (already a known corpus gap) | https://doi.org/10.1103/PhysRevD.42.1424 |

For every one of these except Sather–Schmidt and the two NPB papers, a **free
secondary** that restates the needed content is now on disk: Close–Kumano's sum
rule is restated with its antiquark correction in `1005.4524.pdf` §I and Eq. (9);
*Nuclear gluonometry*'s definition of Δ is quoted in `1803.11206.pdf`;
Khan–Hoodbhoy's convolution structure is superseded in `1702.05337.pdf`
Eqs. (16)–(21).

**One corpus-hygiene item, for the author to decide.** Two files downloaded in
this pass fill gaps that `refs_dict.json` records as `file: null`:
`JLab_PR12-13-011.pdf` (entry `jlab-pr12-13-011-b1-and-2023-jeopardy`) and —
indirectly — nothing else. This pass is append-only on that dictionary and did
**not** edit the existing entry, so its `file` field still reads `null` while the
PDF now sits beside it. One-line fix, author's call.

---

## 8. Priority list — what this domain is worth doing next

Ranked by (value to a defensible tensor claim) ÷ (effort), using the four-camp
structure the corpus now supports.

| # | action | reference(s) | effort | why |
|---|---|---|---|---|
| 1 | **Add `KumanoFitB1` as a third `TensorSF`** — closed form δ_T w = a x^b (1−x)^c (x₀−x), Set 2 parameters, Q² = 2.5 GeV² | `1005.4524` | 1–2 d | the only b₁ source with published numbers; makes the registry "two models + the fit to the measurement", no digitization |
| 2 | **Write the Δ section against the nuclear lattice number**, retiring the φ-meson scaling | `1709.00395`, `1606.04505`, `1803.11206` | 1 d, docs only | turns `04_theory.md` §4's correct caution into a citable quantitative statement: the nuclear ratio is ~10× below the φ's |
| 3 | **Positivity guards on every invented tensor shape** | `1709.07827` (spin ≤ 1), `2602.11587` (spin-3/2) | 1–2 d | free, model-free, and the only thing that would catch a wrong-sign or over-large scenario today |
| 4 | **Digitize Kumano–Kuroki at Q² = 2.5 GeV²** beside `tables::kB1CdksQ2p5` | `2607.09237` | 0.5–1 d | second convolution at the identical Q²; converts E-2 into a two-calculation spread |
| 5 | **A documented small-x validity band on `b1_convolution`**, replacing the silent table freeze | `hep-ph/9605291`, `hep-ph/9611460`, `hep-ph/9711323`, `nucl-th/9701026` | 0.5–1 d | four papers, two of which disagree about the small-x limit: the band is the honest object, not any one curve |
| 6 | **Gate the tagged VECTOR sector** (E-1 gates only the tensor one) | `2006.03033` | 2–4 d | 19 figures and a closed-form LF polarized spectral function; the clearest unexploited external anchor in the tree |
| 7 | **Second sum rule: ∫dx f₂LT = 0** beside `close_kumano_integral` | `2106.15849` | 1 d | free closed form, different function, same machinery |
| 8 | **Quote FSI on A_zz with a number** instead of "IA only" | `1407.1653` | 0.5 d | the only A_zz calculation with a non-Born effect: sizeable for x > 0.2 inclusive, peaking at p_spec ≈ 300 MeV tagged — exactly the tagged region LiPolGen cares about |
| 9 | Obtain Close–Kumano and *Nuclear gluonometry* from a library | §7 | hours | two gates and one central concept currently cited to papers the project does not hold |

---

## 9. Provenance

Searches run this session: 14 INSPIRE literature-API queries (author, title,
fulltext and journal-reference forms) and 8 arXiv-API `id_list` batches
covering 33 distinct records. Every citation line above was written from the
API response, not from memory. Quantitative statements marked **(PDF-read)**
were extracted with `pdftotext -layout` from the PDF now on disk at the path
named in Table 1. Downloads were restricted to arXiv e-prints, IOP
open-access articles (*J. Phys. Conf. Ser.* is CC-BY) and the two public JLab
proposal PDFs on `jlab.org/exp_prog/proposals/`.
