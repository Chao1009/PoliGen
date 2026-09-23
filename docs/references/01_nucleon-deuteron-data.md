<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# 01 — Nucleon and deuteron data and theory: the references the T1/T2 benchmarks rest on

**Compiled 2026-09-16.** Companion to `00_corpus.md` (what the corpus already
holds) and to `../benchmarking/02_data_nucleon_deuteron.md` (which *surveyed*
this layer and verified that the measurements exist, without fetching a single
paper). This file closes that loop: every reference below was re-verified in
this session against a primary index, and every one that is open access has
been **downloaded into the sibling corpus** at
`/home/cpeng/Projects/polli/PolarizedLithiumSim/refs/`.

**Method.** Identity (title, authors, year, journal, DOI) from the arXiv API
(`export.arxiv.org/api/query`, `id_list`) and the INSPIRE literature API
(`inspirehep.net/api/literature`, by arXiv id or by `j <journal>,<vol>,<page>`).
Claims about *content* below are either quoted from the abstract that was
fetched, or — where the sentence says so — read out of the downloaded PDF's own
text layer in this session. Nothing is quoted from memory. Items with no free
copy carry the exact publisher link to use and are marked **paywalled**, never
silently dropped.

**What is new on disk (31 files fetched in this session).** Thirty arXiv PDFs
plus one freely-hosted SLAC preprint scan (`SLAC-PUB-5814.pdf`, Dasu E140, from
INSPIRE's own file store). A thirty-second, Whitlow R1990 (SLAC-PUB-5284), was
fetched from the same store and then removed as a duplicate: a concurrent
sibling workflow had saved the byte-identical file as `whitlow1990_R.pdf` in the
same minute (md5 `d4b4fcdb…`, matching INSPIRE's own file hash), and that is the
copy referenced below. Three further rows (I-1 … I-3, the PDF-set papers) were
on disk from a sibling workflow before this pass reached them and are listed
for completeness, not claimed as this pass's downloads.

---

## 0. Three findings that fell out of the verification

These are not bibliography. They are corrections and constraints that only
became visible once the papers were on disk.

**(0.1) HERMES's b₁ᵈ is per nucleon, and the paper says so in one line.**
`constants.hpp:89` (`B1_MILLER_TABLE_TO_PER_NUCLEON`) records the per-deuteron /
per-nucleon factor of 2 as an open author decision, argued from an inversion of
Table II. The primary source settles it directly. Read out of
`refs/hep-ex_0506018.pdf` in this session, the sentence under Eq. (5) is:

> "The structure function F2d is calculated as F2d = F2p (1 + F2n/F2p)/2 using
> the parameterizations of the precisely measured structure function F2p [16]
> and F2n/F2p ratio [17]."

with Eq. (5) itself `b₁ᵈ = −(3/2) A_zz^d F₁ᵈ`, `F₁ᵈ = (1 + Q²/ν²) F₂ᵈ / (2x(1+R))`.
F₂ᵈ = (F₂ᵖ + F₂ⁿ)/2 **is** the per-nucleon convention. The registry decision can
now be made against a quoted line from the measurement, not against an inversion.

**(0.2) The R inside HERMES's own b₁ is R1990, not R1998.** Reference [18] of
the same paper is **L. W. Whitlow et al., Phys. Lett. B 250, 193 (1990)** —
verified by reading the reference list out of the PDF. LiPolGen's default R is
`r_sigma_lt` (`sf.hpp:46`, "not a fit") and its real R is `r1998`
(`src/core/sf.cpp:14–95`, the E143 world fit). **Neither is the R the HERMES
b₁ᵈ numbers were extracted with.** Any comparison of a LiPolGen b₁ against
HERMES Table II must either rebuild F₁ᵈ with R1990 or state the R-induced shift
as a systematic; the benchmark harness header (`BENCHMARK_PLAN.md` §8) must
carry the R convention the same way it carries the per-nucleon convention.
Whitlow R1990 is now on disk (`refs/whitlow1990_R.pdf`), so this is a
one-afternoon fix rather than a caveat. The same reference list pins the other
two inputs: [16] is **ALLM97** (Abramowicz–Levy, `hep-ph/9712415`, fetched) for
F₂ᵖ and [17] is **NMC, Amaudruz et al., NPB 371 (1992) 3** for F₂ⁿ/F₂ᵖ (the one
input of the four that has no free copy).

**(0.3) The CD-Bonn coefficients in the tree are a faithful transcription.**
`include/lipolgen/cluster.hpp:155–161` types ten `CD_BONN_C` and eight
`CD_BONN_D`. Machleidt's **Table XX** (page 58 of `refs/nucl-th_0006014.pdf`,
extracted in this session) reads
`C₁ = 0.88472985D+00, C₂ = −0.26408759D+00, C₃ = −0.44114404D−01,
C₄ = −0.14397512D+02, …, C₁₀ = −0.25958894D+03` and
`D₁ = 0.22623762D−01, …, D₈ = −0.16114995D+04`, with C₁₁ given by his Eq. (D23)
and D₉–D₁₁ by Eq. (D24) — **exactly** the tree's arrays and exactly the tree's
"computed, never typed" boundary-condition rule. Eighteen of eighteen digits
match. `00_corpus.md` gap **G8** (the single most load-bearing missing
reference) is closed by `refs/nucl-th_0006014.pdf`, and the transcription gate
`test_cluster.cpp` asserts against a source that is now in the repository.

---

## 1. Table 1 — the references

Tier codes are `BENCHMARK_PLAN.md` §2's (T1 nucleon, T2 deuteron, T4
generator-vs-generator); `theory-input` marks a calculation the code
transcribes rather than a measurement it is compared against. "OA" = open
access. Local paths are relative to `PolarizedLithiumSim/`.

| # | key | reference | arXiv | DOI | OA | local file (new unless noted) | tier | serves | pri |
|---|---|---|---|---|---|---|---|---|---|
| **A. The tensor sector at A = 2** |
| A-1 | `hermes2005-b1-tensor` | A. Airapetian *et al.* (HERMES), *First measurement of the tensor structure function b₁ of the deuteron*, PRL **95** (2005) 242001 | hep-ex/0506018 | 10.1103/PhysRevLett.95.242001 | yes | `refs/hep-ex_0506018.pdf` | T2 | `MillerB1`, `CdksB1`, `azz()`, `B1_MILLER_TABLE_TO_PER_NUCLEON`, `close_kumano_integral` | **HIGH** |
| A-2 | `allm97-f2-parametrization` | H. Abramowicz, A. Levy, *The ALLM parameterization of σ_tot(γ\*p) — an update*, DESY 97-251 | hep-ph/9712415 | — | yes | `refs/hep-ph_9712415.pdf` | T1 | reproducing A-1's F₁ᵈ denominator | MEDIUM |
| A-3 | `poudel2025-tensor-sf-review` | J. Poudel, A. Bacchetta, J.-P. Chen, N. Santiesteban, *Experimental study of tensor structure function of deuteron*, EPJ A **61** (2025) 81 | 2506.04506 | 10.1140/epja/s10050-025-01558-w | yes | `refs/2506.04506.pdf` | T2 | citable status statement: one datum exists | LOW |
| **B. Vector-polarized deuteron (the half that carries the rates)** |
| B-1 | `hermes2007-g1` | A. Airapetian *et al.* (HERMES), *Precise determination of the spin structure function g₁ of the proton, deuteron and neutron*, PRD **75** (2007) 012007 | hep-ex/0609039 | 10.1103/PhysRevD.75.012007 | yes | `refs/hep-ex_0609039.pdf` | T1/T2 | `ToyG1::g1_nucleus`, `a1p`/`a1n` | MEDIUM |
| B-2 | `compass2017-g1d-final` | C. Adolph *et al.* (COMPASS), *Final COMPASS results on the deuteron spin-dependent structure function g₁ᵈ and the Bjorken sum rule*, PLB **769** (2017) 34 | 1612.00620 | 10.1016/j.physletb.2017.03.018 | yes | `refs/1612.00620.pdf` | T1/T2 | `g1_nucleus` at A = 2; best EIC-window overlap | **HIGH** |
| B-3 | `compass2007-g1d` | V. Yu. Alexakhin *et al.* (COMPASS), *The deuteron spin-dependent structure function g₁ᵈ and its first moment*, PLB **647** (2007) 8 | hep-ex/0609038 | 10.1016/j.physletb.2006.12.076 | yes | `refs/hep-ex_0609038.pdf` | T1 | superseded by B-2; kept for the ⁶LiD dilution bookkeeping | LOW |
| B-4 | `e143-1998-g1-g2` | K. Abe *et al.* (E143), *Measurements of the proton and deuteron spin structure functions g₁ and g₂*, PRD **58** (1998) 112003 | hep-ph/9802357 | 10.1103/PhysRevD.58.112003 | yes | `refs/hep-ph_9802357.pdf` | T1 | `g2_ww` — the Wandzura–Wilczek check | MEDIUM |
| B-5 | `e155-1999-deuteron-g1-and-lid-target` | P. L. Anthony *et al.* (E155), PLB **463** (1999) 339 | hep-ex/9904002 | 10.1016/S0370-2693(99)00940-5 | yes | `refs/hep-ex_9904002.pdf` **(already in corpus)** | T1 | widest-Q² g₁ᵈ | — |
| **C. Unpolarized F₂ and the ratios** |
| C-0 | `hermes2011-f2p-f2d` | A. Airapetian *et al.* (HERMES), *Inclusive measurements of inelastic electron and positron scattering from unpolarized hydrogen and deuterium targets*, JHEP **05** (2011) 126 | 1103.5704 | 10.1007/JHEP05(2011)126 | yes | `refs/1103.5704.pdf` | T1/T2 | **HERMES's own F₂ᵈ in the kinematics of A-1** — the F₁ᵈ Eq. (5) divides by | **HIGH** |
| C-1 | `nmc1997-f2p-f2d-r` | M. Arneodo *et al.* (NMC), *Measurement of the proton and deuteron structure functions F₂ᵖ and F₂ᵈ, and of the ratio σ_L/σ_T*, NPB **483** (1997) 3 | hep-ph/9610231 | 10.1016/S0550-3213(96)00538-X | yes | `refs/hep-ph_9610231.pdf` | T1 | `ToyF2::f2p/f2n`, `LhapdfSF`, `MstwSF`; **and R at low x** | **HIGH** |
| C-2 | `nmc1997-f2d-over-f2p` | M. Arneodo *et al.* (NMC), *Accurate measurement of F₂ᵈ/F₂ᵖ and R^d − R^p*, NPB **487** (1997) 3 | hep-ex/9611022 | 10.1016/S0550-3213(96)00673-6 | yes | `refs/hep-ex_9611022.pdf` | T1 | `f2n_over_f2p`; R^d − R^p ≈ 0 justifies one R for p and d | MEDIUM |
| C-3 | `nmc1995-f2-parametrization` | M. Arneodo *et al.* (NMC), *Measurement of the proton and the deuteron structure functions F₂ᵖ and F₂ᵈ*, PLB **364** (1995) 107 | hep-ph/9509406 | 10.1016/0370-2693(95)01318-9 | yes | `refs/hep-ph_9509406.pdf` | T1 | the closed-form F₂ parametrisation `ToyF2` is "anchored by eye" against | MEDIUM |
| C-4 | `nmc1992-f2n-over-f2p` | P. Amaudruz *et al.* (NMC), *The ratio F₂ⁿ/F₂ᵖ in deep inelastic muon scattering*, NPB **371** (1992) 3 | — | 10.1016/0550-3213(92)90227-3 | **no** | — (paywalled) | T1 | **the F₂ⁿ/F₂ᵖ inside HERMES's b₁ extraction** (A-1 ref. [17]) | MEDIUM |
| C-5 | `bcdms1990-f2d-r` | A. C. Benvenuti *et al.* (BCDMS), *A high statistics measurement of the deuteron structure functions F₂(x,Q²) and R from deep inelastic muon scattering at high Q²*, PLB **237** (1990) 592 | — | 10.1016/0370-2693(90)91231-Y | **no** | — (paywalled) | T1 | the high-x, high-Q² F₂ᵈ anchor | MEDIUM |
| C-6 | `marathon2022-f2n-f2p` | D. Abrams *et al.* (JLab Hall A Tritium / MARATHON), PRL **128** (2022) 132003 | 2104.05850 | 10.1103/PhysRevLett.128.132003 | yes | `refs/2104.05850.pdf` | T1 | `f2n_over_f2p` at 0.19 < x < 0.83, modern | MEDIUM |
| C-7 | `h1-2008-fl-low-x` | F. D. Aaron *et al.* (H1), *Measurement of the proton structure function F_L at low x*, PLB **665** (2008) 139 | 0805.2809 | 10.1016/j.physletb.2008.05.070 | yes | `refs/0805.2809.pdf` | T1 | **R at EIC x**: 0.00024 < x < 0.0036 | **HIGH** |
| C-8 | `h1-2014-fl` | V. Andreev *et al.* (H1), *Measurement of inclusive ep cross sections at high Q² at √s = 225 and 252 GeV and of the longitudinal proton structure function F_L at HERA*, EPJC **74** (2014) 2814 | 1312.4821 | 10.1140/epjc/s10052-014-2814-6 | yes | `refs/1312.4821.pdf` | T1 | the widest H1 F_L set, 1.5 ≤ Q² ≤ 800 GeV² | **HIGH** |
| **D. R = σ_L/σ_T — the factor that multiplies every F₁** |
| D-1 | `whitlow1990-r1990` | L. W. Whitlow, S. Rock, A. Bodek, E. M. Riordan, S. Dasu, *A precise extraction of R = σ_L/σ_T from a global analysis of the SLAC deep inelastic e-p and e-d scattering cross sections*, PLB **250** (1990) 193; SLAC-PUB-5284 | — | 10.1016/0370-2693(90)91176-C | **yes, preprint** | `refs/whitlow1990_R.pdf` | T1/T2 | **the R of the HERMES b₁ extraction**; the "R1990-like magnitude" `sf.hpp:43` claims for `r_sigma_lt` | **HIGH** |
| D-2 | `dasu1994-e140-r` | S. Dasu *et al.* (SLAC E140), *Measurement of kinematic and nuclear dependence of R = σ_L/σ_T in deep inelastic electron scattering*, PRD **49** (1994) 5641; SLAC-PUB-5814 | — | 10.1103/PhysRevD.49.5641 | **yes, preprint** | `refs/SLAC-PUB-5814.pdf` | T1/T3 | the largest single R set feeding `r1998`; also R_A − R_D | MEDIUM |
| **E. Tagged spectator DIS at A = 2 — the closest analogue of the tagged channel** |
| E-1 | `bonus2012-f2n-tagged` | N. Baillie *et al.* (CLAS BONuS), *Measurement of the neutron F₂ structure function via spectator tagging with CLAS*, PRL **108** (2012) 142001 | 1110.2770 | 10.1103/PhysRevLett.108.142001 | yes | `refs/1110.2770.pdf` | T2 | `TaggedChannel`/`TaggedModel` spectator spectrum shape | MEDIUM |
| E-2 | `bonus2014-f2n-tagged` | S. Tkachenko *et al.* (CLAS BONuS), PRC **89** (2014) 045206 | 1402.2477 | 10.1103/PhysRevC.89.045206 | yes | `refs/1402.2477.pdf` | T2 | same, with the full systematics and the on-shell extrapolation | MEDIUM |
| E-3 | `melnitchouk1997-tagged-sf` | W. Melnitchouk, M. Sargsian, M. Strikman, *Probing the origin of the EMC effect via tagged structure functions of the deuteron*, Z. Phys. A **359** (1997) 99 | nucl-th/9609048 | 10.1007/s002180050372 | yes | `refs/nucl-th_9609048.pdf` | theory-input | the pole-extrapolation logic under `TaggedModel`'s backward-spectator argument | MEDIUM |
| E-4 | `sargsian2010-vna` | M. M. Sargsian, *Large Q² electrodisintegration of the deuteron in virtual nucleon approximation*, PRC **82** (2010) 014612 | 0910.2016 | 10.1103/PhysRevC.82.014612 | yes | `refs/0910.2016.pdf` | theory-input | `GlauberFsiWeight`: the GEA the tree's FSI weight is an approximation of | **HIGH** |
| E-5 | `cosyn-sargsian2017-fsi-review` | W. Cosyn, M. Sargsian, *Nuclear final-state interactions in deep inelastic scattering off the lightest nuclei*, IJMPE **26** (2017) 1730004 | 1704.06117 | 10.1142/S0218301317300041 | yes | `refs/1704.06117.pdf` | theory-input | `[CS17]`: the σ_XN(W) Deeps fit `fsi.hpp` uses, Eq. (11) | **HIGH** |
| E-6 | `egiyan2007-deeps` | K. S. Egiyan *et al.* (CLAS), *Experimental study of exclusive ²H(e,e′p)n reaction mechanisms at high Q²*, PRL **98** (2007) 262502 | nucl-ex/0701013 | 10.1103/PhysRevLett.98.262502 | yes | `refs/nucl-ex_0701013.pdf` | T2 | the data E-4/E-5 fit: where PWIA ends and rescattering begins in p_n | MEDIUM |
| **F. The deuteron wave function and spectral function** |
| F-1 | `machleidt2001-cdbonn` | R. Machleidt, *The high-precision, charge-dependent Bonn nucleon-nucleon potential (CD-Bonn)*, PRC **63** (2001) 024001 | nucl-th/0006014 | 10.1103/PhysRevC.63.024001 | yes | `refs/nucl-th_0006014.pdf` | theory-input | `CdBonnWave`, `CD_BONN_C/D` (Table XX), `A_S`/`η`/`P_D` (Table XV) — **corpus gap G8** | **HIGH** |
| F-2 | `wiringa1995-av18` | R. B. Wiringa, V. G. J. Stoks, R. Schiavilla, *An accurate nucleon-nucleon potential with charge-independence breaking*, PRC **51** (1995) 38 | nucl-th/9408016 | 10.1103/PhysRevC.51.38 | yes | `refs/nucl-th_9408016.pdf` | theory-input | `[AV18]`: the Hamiltonian behind every ANL VMC table in `data/vmc/` and behind `fdeut.av18` | **HIGH** |
| F-3 | `ciofi-simula1996-spectral` | C. Ciofi degli Atti, S. Simula, *Realistic model of the nucleon spectral function in few- and many-nucleon systems*, PRC **53** (1996) 1689 | nucl-th/9507024 | 10.1103/PhysRevC.53.1689 | yes | `refs/nucl-th_9507024.pdf` | theory-input | `CiofiSimulaTriton`, `cs_n0_terms`/`cs_n1_terms` (Table A.1, Eqs. 75–76) | MEDIUM |
| F-4 | `rodning-knutson1990-eta-d` | N. L. Rodning, L. D. Knutson, *Asymptotic D-state to S-state ratio of the deuteron*, PRC **41** (1990) 898 | — | 10.1103/PhysRevC.41.898 | **no** | — (paywalled) | T2 | the *measured* η_d = 0.0256(4) that `CdBonnWave::eta()` = 0.0255714 is 0.07 σ from | MEDIUM |
| F-5 | `boeglin2011-deuteron-high-k` | W. U. Boeglin *et al.* (Hall A), *Probing the high momentum component of the deuteron at high Q²*, PRL **107** (2011) 262501 | 1106.0275 | 10.1103/PhysRevLett.107.262501 | yes | `refs/1106.0275.pdf` | T2 | `Hulthen` vs `VmcAV18` out to 550 MeV/c — turns T25 from model-vs-model into model-vs-data | MEDIUM–HIGH |
| **G. Deuteron elastic form factors (the RC elastic tail's control channel)** |
| G-1 | `abbott1999-deuteron-a` | D. Abbott *et al.* (JLab t₂₀), *A precise measurement of the deuteron elastic structure function A(Q²)*, PRL **82** (1999) 1379 | nucl-ex/9810017 | 10.1103/PhysRevLett.82.1379 | yes | `refs/nucl-ex_9810017.pdf` | T2 | A(Q²) at 0.66–1.8 (GeV/c)²; the magnitude of σ^el(d) | MEDIUM |
| G-2 | `abbott2000-t20` | D. Abbott *et al.* (JLab t₂₀), *Measurement of tensor polarization in elastic electron-deuteron scattering at large momentum transfer*, PRL **84** (2000) 5053 | nucl-ex/0001006 | 10.1103/PhysRevLett.84.5053 | yes | `refs/nucl-ex_0001006.pdf` | T2 | t₂₀, t₂₁, t₂₂ → G_C and G_Q separately | MEDIUM |
| G-3 | `abbott2000-ff-parametrization` | D. Abbott *et al.*, *Phenomenology of the deuteron electromagnetic form factors*, EPJ A **7** (2000) 421 | nucl-ex/0002003 | 10.1007/PL00013629 | yes | `refs/nucl-ex_0002003.pdf` | T2 | **three closed-form world-data parametrisations of G_C, G_M, G_Q (0–7 fm⁻¹)** to replace `HoSpin1FF`'s stand-in shape for the deuteron | MEDIUM |
| G-4 | `garcon-vanorden2001-deuteron` | M. Garçon, J. W. Van Orden, *The deuteron: structure and form factors*, Adv. Nucl. Phys. **26** (2001) 293 | nucl-th/0102049 | 10.1007/0-306-47915-X_4 | yes | `refs/nucl-th_0102049.pdf` | theory-input | the review that fixes conventions (G_C/G_M/G_Q vs A/B/t₂₀, spherical vs Cartesian) | MEDIUM |
| **H. Measured rank-2 observables on a spin-1 nucleus — the Dubna proxy** |
| H-1 | `ladygin2006-ayy-9gev` | V. P. Ladygin *et al.*, *Tensor A_yy and vector A_y analyzing powers in H(d,d′)X and ¹²C(d,d′)X at 9 GeV/c*, Phys. At. Nucl. **69** (2006) 852 | nucl-ex/0510050 | 10.1134/S1063778806050073 | yes | `refs/nucl-ex_0510050.pdf` | T2 (proxy `U-2`) | S–D interference at internal momenta 0.2–1.0 GeV/c; `A_zz^tag`'s only measured analogue | MEDIUM |
| H-2 | `azhgirey2004-ayy-breakup` | L. S. Azhgirey *et al.*, *New data on tensor analyzing power A_yy of the relativistic deuteron breakup…*, PLB **595** (2004) 151 | — | 10.1016/j.physletb.2004.05.057 | **no** | — (paywalled) | T2 (proxy) | same, at 4.5–9 GeV/c | LOW–MED |
| H-3 | `azhgirey1996-t20-breakup` | L. S. Azhgirey *et al.*, *Measurement of the tensor analyzing power T₂₀ in inclusive deuteron breakup at 9 GeV/c on hydrogen and carbon*, PLB **387** (1996) 37 | — | 10.1016/0370-2693(96)01007-6 | **no** | — (paywalled) | T2 (proxy) | the spherical-convention sibling of H-1/H-2 | LOW–MED |
| H-4 | `ladygin2002-ayy-4.5gev` | V. P. Ladygin *et al.*, *Measurement of the tensor-analyzing power A_yy in deuteron breakup at 4.5 GeV/c and 80 mr*, Few Body Syst. **32** (2002) 127 | — | 10.1007/s00601-002-0115-3 | **no** | — (paywalled) | T2 (proxy) | same | LOW |
| **I. The fits the nucleon inputs come from (pedigree, not data)** |
| I-1 | `ct18-2021` | T.-J. Hou *et al.*, *New CTEQ global analysis of QCD with high-precision data from the LHC*, PRD **103** (2021) 014013 | 1912.10053 | 10.1103/PhysRevD.103.014013 | yes | `refs/1912.10053.pdf` (fetched by a sibling workflow this session) | T1 | `LhapdfSF`, `--unpol-sf ct18nlo` | MEDIUM |
| I-2 | `mstw2008` | A. D. Martin, W. J. Stirling, R. S. Thorne, G. Watt, *Parton distributions for the LHC*, EPJC **63** (2009) 189 | 0901.0002 | 10.1140/epjc/s10052-009-1072-5 | yes | `refs/0901.0002.pdf` (sibling) | T1/T2 | `MstwSF` — **the PDF CDKS computed b₁ᵈ with**, hence the A = 2 gate's passing configuration | **HIGH** |
| I-3 | `nnpdfpol11` | E. R. Nocera, R. D. Ball, S. Forte, G. Ridolfi, J. Rojo, *A first unbiased global determination of polarized PDFs and their uncertainties*, NPB **887** (2014) 276 | 1406.5539 | 10.1016/j.nuclphysb.2014.08.008 | yes | `refs/1406.5539.pdf` (sibling) | T1 | `LhapdfG1`, `--pol-sf nnpdfpol` | LOW |
| I-4 | `ye2018-nucleon-ff` | Z. Ye, J. Arrington, R. J. Hill, G. Lee, *Proton and neutron electromagnetic form factors and uncertainties*, PLB **777** (2018) 8 | 1707.09063 | 10.1016/j.physletb.2017.11.023 | yes | `refs/1707.09063.pdf` | T1 | `nucleon_ff` (`rc.hpp:665`) — a world fit **with an error band**, replacing dipole + Galster | LOW |

---

## 2. What each one contains, why the project needs it, what it serves

### A-1 HERMES b₁ᵈ — the only tensor-polarized DIS measurement in existence

Twenty-seven-GeV positrons on a tensor-polarized deuterium gas target with
negligible residual vector polarization; the abstract states the range
0.01 < ⟨x⟩ < 0.45, 0.5 < ⟨Q²⟩ < 5 GeV², and that A_zz^d and b₁ᵈ are non-zero.
**Table II, read from the PDF in this session, is six rows** — ⟨x⟩ = 0.012,
0.032, 0.063, 0.128, 0.248, 0.452 at ⟨Q²⟩ = 0.51, 1.06, 1.65, 2.33, 3.11,
4.69 GeV², A_zz^d ×10² = −1.06 ± 0.52 ± 0.26, −1.07 ± 0.49 ± 0.36,
−1.32 ± 0.38 ± 0.21, −0.19 ± 0.34 ± 0.29, −0.39 ± 0.39 ± 0.32, +1.57 ± 0.68 ± 0.13
and b₁ᵈ ×10² = 11.20 ± 5.51 ± 2.77, 5.50 ± 2.53 ± 1.84, 3.82 ± 1.11 ± 0.60,
0.29 ± 0.53 ± 0.44, 0.29 ± 0.28 ± 0.24, −0.38 ± 0.16 ± 0.03. These are
digit-for-digit the rows already transcribed in
`docs/open_items/run_2026-09-03/phase_A_miller_normalisation.md` §2.3, so that
transcription is now verified against the source. The paper also publishes two
integrals, both after evolving to Q₀² = 5 GeV² assuming b₁ᵈ/F₁ᵈ is Q²-independent:
∫₀.₀₀₂^0.85 b₁ dx = (1.05 ± 0.34 ± 0.35)×10⁻² over the measured range, and
∫₀.₀₂^0.85 b₁ dx = (0.35 ± 0.10 ± 0.18)×10⁻², a 1.7 σ result, over the restricted
range where Q² > 1 GeV². **That second number is the one to compare
`close_kumano_integral` against** — LiPolGen reports ∫b₁ dx over [0.01, 1.2] and
never enforces it; HERMES gives a measured value with a stated x-window and a
stated evolution assumption, which is what a sum-rule comparison needs.
It serves `MillerB1`/`toy_b1`, `CdksB1`/`b1_convolution`, `azz()`
(`asymmetries.hpp:120`) through A_zz = −(2/3) b₁/F₁, the |A_zz| = 1.06×10⁻² spot
value at `rc.hpp:272`, and — per §0.1 — the `B1_MILLER_TABLE_TO_PER_NUCLEON`
decision. It does **not** validate anything at A = 6: b₁ᵈ is the *input* to
`Li6B1`/`b1_li6_from_deuteron`, so reproducing it does not test the convolution.

### A-2 ALLM97 — the F₂ᵖ parametrisation inside A-1

Abramowicz and Levy's update of the ALLM Regge-plus-hard-Pomeron form for
σ_tot(γ\*p), fitted to 1356 F₂ points with χ²/ndf = 0.97 over
3×10⁻⁶ < x < 0.85, 0 ≤ Q² < 5000 GeV² (abstract, verified). It is in this file
for exactly one reason: it is reference [16] of the HERMES b₁ paper, so
rebuilding HERMES's F₁ᵈ — the denominator of every b₁ᵈ row in Table II — needs
it, together with C-4's ratio and D-1's R. Nothing in LiPolGen currently uses
ALLM; it is a benchmark-harness input, not a kernel.

### A-3 Poudel *et al.* 2025 — the citable status statement

A review of the experimental tensor programme: inclusive and semi-inclusive DIS
on tensor-polarized deuterons, the polarization technology, and what is
approved. No data of its own. Its value is that release notes and the letter can
cite a 2025 peer-reviewed sentence for "one measurement exists, the second has
not run", rather than asserting it.

### B-1 … B-4 The vector-polarized deuteron

**B-1 (HERMES PRD 75)** gives g₁ᵖ, g₁ᵈ and an extracted g₁ⁿ over
0.0041 ≤ x ≤ 0.9, 0.18 ≤ Q² ≤ 20 GeV² (abstract) with 23 HEPData tables
(`10.17182/hepdata.11211`). **B-2 (COMPASS PLB 769)** is the final COMPASS g₁ᵈ,
1 < Q² < 100 (GeV/c)², 0.004 < x < 0.7, W > 4 GeV/c² — **the best overlap with
the EIC window of any g₁ᵈ set**, and the single cleanest external gate on
`ToyG1::g1_nucleus` evaluated on `DEUTERON()`. **B-3** is its superseded 2007
predecessor, kept because it documents the ⁶LiD target bookkeeping COMPASS used.
**B-4 (E143)** carries g₂ᵈ at 29.1 GeV and states in its own abstract that "the
g₂ data are well-described by the Wandzura–Wilczek twist-2 contribution" —
which is precisely what `g2_ww` implements, and today `g2_ww` is tested only
against its own analytic power-law limit (`test_sf.cpp:147`). None of B-1…B-4
constrains anything tensor: a vector-polarized target has no A_zz sensitivity,
and they do not touch `eff_pol_p`/`eff_pol_n`, which are VMC sums.

### C-0 HERMES F₂ᵈ — the same apparatus as the b₁ measurement

Inclusive inelastic e± scattering on unpolarized hydrogen and deuterium at
HERMES, 0.006 ≤ x ≤ 0.9, 0.1 ≤ Q² ≤ 20 GeV², with F₂ᵖ and F₂ᵈ determined using a
parametrisation of existing data for R (abstract), and the deuteron-to-proton
cross-section ratio given as well; three HEPData tables
(`10.17182/hepdata.66147`). Its value here is specific: `constants.hpp:89`
argues that HERMES's published b₁ᵈ is per nucleon by *inverting* their Table II
through Eq. (5). §0.1 now replaces that inversion with the paper's own sentence,
and C-0 is the cross-check — HERMES's own measured F₂ᵈ, in the same apparatus
and overlapping kinematics, so the inversion can be re-done against a published
table rather than a parametrisation. It also serves `ToyF2`/`LhapdfSF`/`MstwSF`
directly in the x range where HERMES and the EIC overlap. It does **not**
supersede C-1 for the low-x F₂ᵈ: HERMES's √s = 7.2 GeV reach stops well short of
NMC's.

### C-1 … C-3 NMC — the F₂ᵈ and the ratio the toy is anchored to

**C-1** is the full NMC F₂ᵖ and F₂ᵈ set, 0.002 < x < 0.60, 0.5 < Q² < 75 GeV²,
**plus R for 0.002 < x < 0.12** (abstract) — it is simultaneously the best
fixed-target F₂ᵈ for the EIC window and the only R data at NMC's low x, where
`r1998`'s support (`R1998_X_MIN` = 0.005) runs out and the tree clips. **C-2**
measures F₂ᵈ/F₂ᵖ over 0.001 < x < 0.8, 0.1 < Q² < 145 GeV² and finds R^d − R^p
compatible with zero over 0.002 < x < 0.4 — the empirical licence for using one
R for both p and d, which `f1_from_f2` does implicitly everywhere. **C-3** is
the NMC paper whose closed-form parametrisation of F₂ᵖ and F₂ᵈ (0.006 < x < 0.9)
is the natural replacement for `ToyF2`'s by-eye anchors when LHAPDF is not
built. Together they serve `ToyF2::f2p`/`f2n`/`f2n_over_f2p`, `LhapdfSF` and
`MstwSF`, and they bound nothing nuclear: `NuclearF2::f2a` and
`unpolarized_emc_ratio` need A > 2 data (file `02` of this series).

### C-4 NMC 1992 F₂ⁿ/F₂ᵖ — paywalled, and load-bearing anyway

Amaudruz *et al.*, NPB 371 (1992) 3, verified through INSPIRE (recid 321412,
69 authors, CERN-PPE-91-167). No arXiv (1991), no free copy located: INSPIRE
holds no fulltext, and CERN's document server was behind an interactive bot
challenge from this environment, so the preprint could not be checked there.
**Download link for the user: <https://doi.org/10.1016/0550-3213(92)90227-3>.**
It matters because it is the F₂ⁿ/F₂ᵖ that HERMES used to build the F₂ᵈ in their
Eq. (5) (§0.1): reproducing their b₁ᵈ from their A_zz^d exactly requires it.
For the tree's own `f2n_over_f2p`, C-2 and C-6 are adequate substitutes and are
both free.

### C-5 BCDMS — paywalled

Benvenuti *et al.*, PLB 237 (1990) 592 (INSPIRE 285497, 39 authors, BCDMS,
CERN-EP-89-170). High-x, high-Q² F₂ᵈ and R from muon scattering; the anchor at
the top end of the x range where `ToyF2`'s by-eye shape is least constrained.
Same availability story as C-4. **Download link:
<https://doi.org/10.1016/0370-2693(90)91231-Y>.**

### C-6 MARATHON — the modern F₂ⁿ/F₂ᵖ

³H/³He mirror symmetry at 0.19 < x < 0.83, which "essentially eliminates many
theoretical uncertainties in the extraction of the ratio" (abstract). It is the
cleanest available test of `f2n_over_f2p` at large x — the region where the
toy's single-argument, Q²-frozen hook is most exposed — and, being a tritium
measurement, it is also the only modern experimental contact with the A = 3
system that `CiofiSimulaTriton` models. It is *not* a spectral function: it
measures a ratio of structure functions, not S(k, E).

### C-7, C-8 H1 F_L — the only R inside the EIC's low-x window

**C-7** measures F_L over 12 < Q² < 90 GeV², 0.00024 < x < 0.0036 (abstract);
**C-8** extends the extraction to 1.5 ≤ Q² ≤ 800 GeV² using the 225 and 252 GeV
runs. With R = F_L/(F₂ − F_L), these are the only measured R at x below 10⁻³ —
where LiPolGen runs and where `r1998` is frozen at its boundary by
`clip = true`. Wiring them puts a number on the cost of that extrapolation, and
on the bias `r_sigma_lt` (the *default*, and by its own docstring "not a fit")
carries into every `f1_from_f2`, every b₁/F₁ ratio and the whole A_zz
normalisation. Caveat to state on the harness row: these are **proton** F_L; no
deuteron F_L exists at that x.

### D-1 Whitlow R1990 — free, and more important than it looks

The global analysis of eight SLAC e-p and e-d experiments (1970–1985); the
preprint's own abstract, read from `refs/whitlow1990_R.pdf`, gives R^p, R^d and
R^d − R^p "over the entire SLAC kinematic range: 0.1 ≤ x ≤ 0.9 and
0.6 ≤ Q² ≤ 20.0 (GeV/c)²" and finds R^p = R^d. Two distinct jobs: (i) it is the
R that the HERMES b₁ᵈ values were extracted with (§0.2), so it belongs in any
comparison against Table II; (ii) it is the actual reference behind
`sf.hpp:43`'s unsourced claim that `r_sigma_lt` has a "simplified R1990-like
magnitude" — a claim that can now be checked rather than repeated. The freely
hosted copy is the SLAC-PUB-5284 scan in INSPIRE's file store
(`https://inspirehep.net/files/d4b4fcdb5d834b34d41d3226090c37ea`); the journal
version is paywalled at <https://doi.org/10.1016/0370-2693(90)91176-C>.

### D-2 Dasu E140 — the largest R data set feeding r1998

Kinematic *and nuclear* dependence of R (INSPIRE 360765, SLAC-PUB-5814, 30
authors; the 0.2 ≲ x ≲ 0.5, 1 ≲ Q² ≲ 10 GeV² range is `02`'s row R-2, carried
forward, not re-read from the scan). It is the dominant input to the
E143 R1998 world fit that `src/core/sf.cpp` transcribes, so comparing `r1998`
to E140 is the closest thing available to testing the transcription against the
data rather than against itself — today `tests/test_sf.cpp:72–86` only checks
internal consistency and pins `r1998(0.1, 10) = 0.12654600717127532`. It also
carries R_A − R_D, which is the only handle anywhere on whether R on a *nucleus*
differs from R on the deuteron — the assumption `NuclearF2` and the ⁶Li R
default both make silently. The INSPIRE-hosted preprint scan is a 4.8 MB image
PDF; its OCR text layer is poor, so the tables must be read from the page
images.

### E-1, E-2 CLAS BONuS — the tagged topology at A = 2

DIS on deuterium with a backward spectator proton, p_s < 100 MeV/c,
θ_pq > 100°, 0.65 < Q² < 4.52 GeV², 0.2 < x < 0.8 (abstracts). This is
LiPolGen's tagged topology with a nucleon instead of an α, and with the same
physical argument — a slow backward spectator selects a nearly on-shell struck
nucleon and suppresses FSI. E-2 adds the full systematic treatment and the
radial-TPC acceptance (protons 70–100 MeV/c). They serve `TaggedChannel`,
`TaggedModel`'s (α_s, p_T) grid and dilutions, `TaggedSampler`,
`boost_spectator`, and the never-rescattered limit of `GlauberFsiWeight`.
The frame mismatch must be stated on any harness row: BONuS tags backward in
the target rest frame, the EIC tags a boosted α or d in B0/OMD/ZDC acceptance,
and the acceptance-driven systematics do not transfer.

### E-3 Melnitchouk–Sargsian–Strikman — why backward tagging works

The paper that established tagged structure functions of the deuteron as an EMC
discriminant: choosing extreme backward spectator kinematics minimises wave
function and FSI sensitivity so that bound-nucleon modification can be isolated
within the impulse approximation (abstract). It is the theory citation for the
argument `tagged.hpp` makes when it treats a slow backward spectator as a clean
tag, and it is the origin of the pole-extrapolation idea the tree references
through `[Sargsian–Strikman]`.

### E-4 Sargsian VNA — the FSI calculation `GlauberFsiWeight` approximates

Two-body deuteron breakup at high Q² in the virtual nucleon approximation, with
final-state interactions computed in the **generalised eikonal approximation**
and an explicit study of the calculation's uncertainties against the first
large-Q² electrodisintegration data (abstract). LiPolGen's `fsi.hpp` implements
a Glauber weight with a σ_XN(W) fit; this is the paper that says what that
weight is an approximation *of*, in which kinematic regions it is controlled,
and where the eikonal breaks. It is the single most important theory reference
for the tagged channel's systematic band.

### E-5 Cosyn–Sargsian review — the σ_XN(W) fit in the code

`PHYSICS_CHANNELS.md` cites `[CS17]` Eq. (11) —
f_XN = σ_tot(i + ε) e^{Bt/2}, ε = −0.5, B ≈ 6 GeV⁻², fitted to Deeps — as the
source of the effective rescattering amplitude. The review's abstract confirms
the framing: the FSI of DIS products with nuclear fragments described by an
effective eikonal amplitude "whose parameters can be extracted from the analysis
of semi-inclusive DIS off the deuteron target", with the Q² and W dependences as
a new observable. The paper was cited in the tree and absent from the corpus;
it is now on disk, so the Eq. (11) parameters `fsi.hpp` hard-codes can be
checked against their source.

### E-6 Egiyan Deeps — the data under E-4 and E-5

²H(e,e′p)n with full kinematic coverage, 1.75 < Q² < 5.5 GeV². Its abstract
gives the three-region structure the tagged channel lives in: p_n < 100 MeV/c
the neutron is a spectator and PWIA describes the reaction; 100 < p_n < 750
MeV/c pn rescattering dominates; above that Δ production followed by NΔ → NN
takes over. That is a directly usable, *measured* statement of where a
plane-wave tagged calculation stops being valid — i.e. an external bound on the
regime in which `GlauberFsiWeight` must carry its full band rather than a small
correction.

### F-1 CD-Bonn — corpus gap G8, closed

The potential paper whose Appendix D and **Table XX** are transcribed verbatim
into `cluster.hpp` (verified digit-by-digit in §0.3), and whose **Table XV**
publishes the deuteron observables the tree reproduces: B_d, A_S = 0.8846,
η = 0.0256, P_D = 4.85 %, Q_d. The live test run recorded in
`../benchmarking/02_data_nucleon_deuteron.md` §5.1 gives 4.85621 %, 0.88473 and
0.0255714 against those published values. Until this session that gate compared
the code to numbers typed from a paper nobody held; the paper is now in the
corpus.

### F-2 AV18 — the Hamiltonian under every VMC table

Argonne v₁₈: fourteen charge-independent operators plus three charge-dependent
and one charge-asymmetric, fitted to the Nijmegen pp and np database with
χ²/datum = 1.09 for 4301 data (abstract). `PHYSICS_CHANNELS.md` note 27 flags it
as "not cited with a specific reference in the repository" — yet
`data/vmc/deuteron/fdeut.av18`, every `momenta/*.momentum` file, the `overlap_old`
sign reference and `deuteron_av18_p_d()` are all AV18(+UX/UIX) outputs. It is
the provenance of the tree's largest single data directory, and it is now on
disk.

### F-3 Ciofi degli Atti–Simula — a fit to a calculation, honestly labelled

The factorised-ansatz spectral function whose Table A.1 (n₀ for A = 2, 3, 4),
Eq. (76) (n₁ for A = 3) and Eq. (75) + Table A.3 (A = 4) are transcribed into
`src/core/triton_sf.cpp`. Its abstract is explicit that validation is against
**many-body calculations** ("the Spectral Functions of ³He and infinite nuclear
matter resulting from the convolution formula and from many-body calculations
are compared"), not against data — so `cs_n0_terms(3)`/`cs_n1_terms(3)` are a
fit to a calculation, and the tree additionally uses the ³He *proton* column for
a triton by isospin-mirror reasoning. Both statements belong in
`triton_sf.cpp`'s docstring; the paper being on disk is what lets that be
written with a page reference rather than from memory.

### F-4 Rodning–Knutson η_d — paywalled, one number, high leverage

PRC 41 (1990) 898 (INSPIRE 312974). The measurement of the deuteron's asymptotic
D/S ratio, η_d = 0.0256(4), from sub-Coulomb d + ⁴He tensor analysing powers —
the same group as the George–Knutson η(⁶Li → α+d) (a restricted phase-shift
analysis, not this technique), which the tree carries as a consistency band on the quadrupole dial. η is an *asymptotic* observable, so it constrains the tail of
w(r) directly and is not representation-dependent the way P_D is. Its role here
is attribution: `test_cluster.cpp:386` already carries 0.0256(4) via Machleidt's
Table XV, but as a model output rather than as a measurement. **Download link:
<https://doi.org/10.1103/PhysRevC.41.898>** (no arXiv; 1990 PRC; no free copy
located).

### F-5 Boeglin — momentum distribution vs data out to 550 MeV/c

d(e,e′p) at Q² = 3.5 (GeV/c)², binned in the neutron recoil angle θ_nq, giving
missing-momentum distributions to 0.55 GeV/c; in 35° ≤ θ_nq ≤ 45° FSI are
predicted small, so the reduced cross sections "provide direct access to the
high momentum component of the deuteron momentum distribution" (abstract). This
is the reference that turns `test_tagged.cpp` T25 — today Hulthén against AV18,
two models — into model-against-data, over exactly the k range the tagged
spectator spectrum populates. It is unpolarized, so it says nothing about
`eff_pol_p`/`eff_pol_n` or `vector_dilution_of`, and it is a rest-frame
measurement, so the light-cone variables do not transfer directly.

### G-1, G-2, G-3, G-4 The deuteron form factors

`tests/test_rc.cpp:374–397` states that the deuteron form factor in the RC
elastic tail is a harmonic-oscillator **shape stand-in** for POLRAD's `ffdeu`,
normalised on the measured μ and Q_d, and that every T8(a)/T10 gate is
shape-blind. **G-3 removes the stand-in**: it publishes three closed-form
parametrisations of G_C, G_M and G_Q over 0–7 fm⁻¹ built from the world data
including the then-new polarization measurements (abstract) — a drop-in
replacement for `HoSpin1FF`'s deuteron branch, and the natural filling for
`TabulatedSpin1FF`. **G-1** and **G-2** are the two JLab measurements that made
that parametrisation possible: A(Q²) at six points between 0.66 and 1.80
(GeV/c)², and t₂₀/t₂₁/t₂₂ at six points over the same range, which separate G_C
from G_Q. **G-4** is the review to read first, because it is where the
conventions live — A and B versus G_C, G_M, G_Q; spherical t₂₀ versus Cartesian
A_zz — and the project has a documented history of convention traps (`06_critic.md`
trap T-8). None of these constrains the **⁶Li** form factors, where the RC
tail's real uncertainty is (`HoSpin1FF`'s ±100 % `fq_scale` band).

### H-1 … H-4 Dubna/JINR — measured rank-2 observables on a spin-1 nucleus

`06_critic.md` §4.3 identified this literature as the closest measured proxy for
`U-2` (the tagged tensor asymmetry with a cluster spectator) and noted that
files `00`–`05` did not contain it. **H-1 is the only one with a free copy**:
tensor A_yy and vector A_y in H(d,d′)X and ¹²C(d,d′)X at 9.0 GeV/c over the
baryon-resonance region, with an approximate t-scaling up to −1.5 (GeV/c)² and
a comparison to plane-wave impulse approximation predictions (abstract).
H-2, H-3 and H-4 are the A_yy and T₂₀ breakup measurements at 4.5–9 GeV/c;
all three are pre-arXiv or non-deposited and paywalled —
<https://doi.org/10.1016/j.physletb.2004.05.057>,
<https://doi.org/10.1016/0370-2693(96)01007-6>,
<https://doi.org/10.1007/s00601-002-0115-3>. What they measure is a rank-2
observable driven by S–D interference at internal momenta 0.2–1.0 GeV/c — the
same `f₀`, `f₂` structure as `A_zz^wf`, far beyond BLAST's ≤ 500 MeV/c. What
they are not: an electromagnetic probe, a tagged α, or the target-polarization
convention. **The analysing-power convention is not the target-polarization
convention** (Madison convention; t₂₀ versus P_zz differ by ~√2 and by sign
usage) — resolve it before any number is compared, exactly as `06_critic.md`
trap T-8 demands.

### I-1 … I-4 The fits

These are not measurements; each carries its own published validation against
the data of §§B–D, so attaching one imports that validation. **I-2 (MSTW2008
LO)** deserves the HIGH: it is the PDF CDKS used to compute the b₁ᵈ the A = 2
gate is defined against, and `OPEN_ITEMS_SOLUTIONS.md` records that the gate's
peak ratio is 0.843 with MSTW and 0.440 — outside the window — with the shipped
`ToyF2`. The gate is therefore a statement about a *configuration*, and the
reference for that configuration now sits in the corpus. **I-4** replaces
`nucleon_ff`'s dipole+Galster with a world fit that has an error band; its own
priority stays LOW because the `qe_suppression` ±100 % band dominates anyway.

---

## 3. Verified absences and unobtainables

* **Paywalled, no free copy located** (each verified to exist through INSPIRE;
  use the DOI link given in §2): C-4 NMC 1992 F₂ⁿ/F₂ᵖ; C-5 BCDMS; F-4
  Rodning–Knutson; H-2/H-3/H-4 the three Dubna breakup papers. For C-4 and C-5
  a CERN preprint (CERN-PPE-91-167, CERN-EP-89-170) plausibly exists on
  `cds.cern.ch`, but CDS answered with an interactive bot challenge from this
  environment and could not be checked — that is "not verified", not "does not
  exist".
* **No HEPData table for HERMES b₁** (A-1). The six rows are table-in-paper,
  now on disk.
* **No R on a nucleus at EIC kinematics, and no tensor L/T separation ever.**
  D-2's R_A − R_D is fixed-target and A > 2; the split of the tensor response
  into T and L pieces (`cosyn_tensor_sfs`) has never been measured.
* **No polarized tagged DIS on the deuteron.** E-1/E-2 are unpolarized.
* **No deuteron F_L at low x.** C-7/C-8 are proton.
* **No second b₁ measurement.** A-3 is the 2025 review that says so.

## 4. Suggested wiring order for this file's rows

1. **A-1 Table II as an assertion** (six rows, now verified against the PDF),
   built with **D-1's R** and, if exactness is wanted, A-2 + C-4 for F₂ᵈ,
   cross-checked against **C-0**. This
   simultaneously settles registry row `B1_MILLER_TABLE_TO_PER_NUCLEON` (§0.1)
   and gives `close_kumano_integral` a measured number to report against.
2. **C-7/C-8 (H1 F_L) + D-1/D-2** — a one-page R comparison that prices
   `r_sigma_lt` against measured R across the whole window, including below
   `R1998_X_MIN` where the fit is clipped.
3. **B-2 (COMPASS final g₁ᵈ)** — one clean external gate on `g1_nucleus`.
4. **F-5 (Boeglin)** — T25 from model-vs-model to model-vs-data.
5. **G-3 (Abbott parametrisations)** — remove the acknowledged HO stand-in in
   the deuteron control channel of the RC tail.
6. **E-5/E-4** — check the σ_XN(W) parameters `fsi.hpp` hard-codes against their
   published source, and state the eikonal's validity window from E-6's measured
   three-region result.
7. Documentation, now sourceable: attribute η_d = 0.0256(4) to F-4; say in
   `triton_sf.cpp` that F-3 is a fit to a calculation with a ³He column standing
   in for ³H; cite F-2 for `data/vmc/`; cite F-1 for `CD_BONN_C/D`.

---

*Every identity claim in this file (title, authors, year, journal, DOI, arXiv
id) was checked in this session against the arXiv API or the INSPIRE literature
API; every content claim is from an abstract fetched in this session or from the
text layer of a PDF downloaded in this session, and the ones read from PDFs say
so. The two HEPData DOIs quoted in §2 (B-1, C-0), and the table counts and
kinematic ranges attributed to `../benchmarking/02_data_nucleon_deuteron.md`,
are that file's INSPIRE-sourced numbers carried forward, not re-checked here:
`hepdata.net` and `cds.cern.ch` were both unreachable from this environment
(interactive bot challenges), so no HEPData table was read and no CERN preprint
was checked.*
