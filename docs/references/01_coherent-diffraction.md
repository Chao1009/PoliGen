# 01 — Coherent diffraction and exclusive production on light nuclei

**Domain.** Coherent (intact-ground-state) diffraction and exclusive vector-meson
production on light nuclei, and the two `cos 2φ` mechanisms the ⁶Li coherent channel
has to tell apart. This is the reference set for `include/lipolgen/coherent.hpp`,
`include/lipolgen/cluster_config.hpp` and open item **O5** (is the coherent tensor
a₂ measurable in one EIC year — `docs/OPEN_ITEMS_SOLUTIONS.md` §11.3).

**Date.** 2026-09-16. **Stage.** Domain search, after `00_corpus.md`'s inventory.

**Method, and what "verified" means here.** Every entry below was checked against
the arXiv API (`export.arxiv.org/api/query`, `id_list=`) or the INSPIRE REST API
(`inspirehep.net/api/literature`) in this pass: title, author list, year, journal
reference, DOI and the abstract text were fetched and read, and the claim each row
makes about the paper's content is either (a) in the fetched abstract, or (b) quoted
from the PDF text of the copy now on disk, with the locating grep stated in the
paragraph. Rows whose content claim rests on the *title alone* say so. Twenty-six
PDFs were downloaded to the sibling corpus `PolarizedLithiumSim/refs/` — all of them
arXiv e-prints (open access); nothing paywalled was fetched.

**Relation to the corpus.** Of the 54 corpus entries, exactly two are in this domain:
`mantysaari2024-polarized-deuteron-imaging` (arXiv:2408.13213, `refs/2408.13213v1.pdf`)
and `chang2026-ir8-light-nuclei` (arXiv:2511.05638, `refs/2511.05638.pdf`). Both are
listed below for completeness with `in_corpus_already = true` and were **not**
re-downloaded. Everything else here is new to the corpus. Several items the domain
brief asked about were already named in `docs/open_items/physics_literature.md` §7 —
named, but never fetched; those are now on disk.

**One correction to the brief.** "Jentsch–Tu–Zhang PRC 104 (2021) 065205" resolves,
on INSPIRE, to **Jentsch, Tu, Weiss**, *Deep-inelastic electron-deuteron scattering
with spectator nucleon tagging*, PRC 104 (2021) 065205, arXiv:2108.08314 — which is
already corpus entry `jentsch-tu-weiss2021-ed-tagging`. There is no separate
Jentsch–Tu–Zhang paper at that citation. Likewise, "Mäntysaari–Roy–Salazar–Schenke
2021 on ³He/⁴He" is arXiv:2011.02464, PRD 103 (2021) 094026 — a real and relevant
paper, but its targets are the **proton and gold**, not ³He/⁴He; the ³He/⁴He coherent
calculation the brief is reaching for is Guzey *et al.*, arXiv:2202.12200 (§B1), and
the ab-initio A = 3, 4 configurations are in Mäntysaari *et al.*, arXiv:2605.00454 (§B2).

---

## Table 1 — the domain's references

`OA` = open access. `local` = file now committed under `PolarizedLithiumSim/refs/`
(⭑ = fetched by this pass; ⊙ = already in the corpus; ✗ = no free copy, link given).
`Serves` names the `coherent.hpp` / `cluster_config.hpp` symbol, the O5 input, or the
open item the reference is for.

| # | Key | Citation | arXiv / DOI | OA | local | Tier | Serves |
|---|---|---|---|---|---|---|---|
| **A. The tensor / polarized question** |
| A1 | `mantysaari2024-polarized-deuteron-imaging` | H. Mäntysaari, F. Salazar, B. Schenke, C. Shen, W. Zhao, *Spatial imaging of polarized deuterons at the EIC*, Phys. Lett. B 858 (2024) 139053 | [2408.13213](https://arxiv.org/abs/2408.13213) · [10.1016/j.physletb.2024.139053](https://doi.org/10.1016/j.physletb.2024.139053) | yes | ⊙ `2408.13213v1.pdf` | theory-input | `mantysaari_a2_deuteron()`, `a2_from_quadrupole` Eq. (G8), the whole `eps_b0` anchor |
| A2 | `dalton2025-halld-tensor-deuteron` | M. M. Dalton, A. Deur, C. Keith, *Potential for Tensor Polarized Deuterons in Hall D at Jefferson Lab*, Eur. Phys. J. A 61 (2025) 111 | [2504.21177](https://arxiv.org/abs/2504.21177) · [10.1140/epja/s10050-025-01580-y](https://doi.org/10.1140/epja/s10050-025-01580-y) | yes | ⭑ `2504.21177.pdf` | theory-input | the **only** other proposal to measure a tensor observable in coherent VM production; the m = 0 / intermediate-\|t\| sensitivity argument behind `a2_m_state` |
| A3 | `good-walker1960` | M. L. Good, W. D. Walker, *Diffraction dissociation of beam particles*, Phys. Rev. 120 (1960) 1857 | no arXiv · [10.1103/PhysRev.120.1857](https://doi.org/10.1103/PhysRev.120.1857) | **no** | ✗ | theory-input | the eigenstate decomposition `write_snd_configs` / `ClusterConfigSampler` exists to average over |
| A4 | `miettinen-pumplin1978` | H. I. Miettinen, J. Pumplin, *Diffraction scattering and the parton structure of hadrons*, Phys. Rev. D 18 (1978) 1696 | no arXiv · [10.1103/PhysRevD.18.1696](https://doi.org/10.1103/PhysRevD.18.1696) | **no** (free scan on INSPIRE, see §A4) | ✗ | theory-input | coherent = ⟨A⟩², incoherent = variance — the statement `ClusterConfigSet` implements |
| A5 | `rachek2017-t20-coherent-pi0` | I. A. Rachek et al., *Measurement of Tensor Analyzing Power T₂₀ in Coherent π⁰ Photoproduction on Deuteron*, Few Body Syst. 58 (2017) 29 | no arXiv · [10.1007/s00601-016-1191-0](https://doi.org/10.1007/s00601-016-1191-0) | **no** | ✗ | theory-input | precedent only: the sole *measured* tensor observable in a coherent reaction on a nucleus |
| **B. Light nuclei with realistic structure — the ⁶Li analogues** |
| B1 | `guzey2022-coherent-jpsi-3he-4he` | V. Guzey, M. Rinaldi, S. Scopetta, M. Strikman, M. Viviani, *Coherent J/ψ electroproduction on ⁴He and ³He at the EIC: probing nuclear shadowing one nucleon at a time*, Phys. Rev. Lett. 129 (2022) 242503 | [2202.12200](https://arxiv.org/abs/2202.12200) · [10.1103/PhysRevLett.129.242503](https://doi.org/10.1103/PhysRevLett.129.242503) | yes | ⭑ `2202.12200.pdf` | theory-input | the closest published calculation to a ⁶Li coherent amplitude: AV18 few-body wave functions, `slope_b` / `sample_t` shape, O5's rate model |
| B2 | `mantysaari2026-nuclear-structure-saturation` | H. Mäntysaari, H. Roch, B. Schenke, C. Shen, W. Zhao, *Nuclear structure and saturation effects from diffractive vector meson production*, Phys. Rev. D 114 (2026) 014068 | [2605.00454](https://arxiv.org/abs/2605.00454) · [10.1103/2628-2spx](https://doi.org/10.1103/2628-2spx) | yes | ⭑ `2605.00454.pdf` | theory-input | kills the "lightest published is Ca" note: GFMC A = 3, 4 and VMC/PGCM/NLEFT configurations, incl. a "PGCM clustering" sampling — the direct precedent for `ClusterConfigSet` |
| B3 | `chang2026-ir8-light-nuclei` | W. Chang, E.-C. Aschenauer, A. Jentsch, A. Kumar, Z. Tu, Z. Yin, *Opportunities for imaging light nuclei at a second IR*, Phys. Rev. D 113 (2026) 032018 | [2511.05638](https://arxiv.org/abs/2511.05638) · [10.1103/y4yv-y9dn](https://doi.org/10.1103/y4yv-y9dn) | yes | ⊙ `2511.05638.pdf` | T4 | `COHERENT_JPSI_EFF_IR8_LI7`, `Chang26EffEnergyRow`, `Chang26SpeciesEffRow` — the whole O5 efficiency chain |
| B4 | `mondal2026-shell-structure-vm` | A. Mondal, A. Kumar, D. Sarkar, *Imprints of nuclear shell structure in exclusive vector meson production*, arXiv:2608.23445 (Aug 2026) | [2608.23445](https://arxiv.org/abs/2608.23445) | yes | ⭑ `2608.23445.pdf` | theory-input | how far a realistic light-nucleus density moves the coherent \|t\| shape — a systematic on `slope_b` beyond the Gaussian |
| B5 | `mantysaari2023-deformation-prl` | H. Mäntysaari, B. Schenke, C. Shen, W. Zhao, *Multi-scale imaging of nuclear deformation at the EIC*, Phys. Rev. Lett. 131 (2023) 062301 | [2303.04866](https://arxiv.org/abs/2303.04866) · [10.1103/PhysRevLett.131.062301](https://doi.org/10.1103/PhysRevLett.131.062301) | yes | ⭑ `2303.04866.pdf` | theory-input | deformation → coherent/incoherent ratio; the *unpolarized* half of the `delta_b_m` geometry |
| B6 | `mantysaari2022-nuclear-geometry` | H. Mäntysaari, F. Salazar, B. Schenke, *Nuclear geometry at high energy from exclusive vector meson production*, Phys. Rev. D 106 (2022) 074019 | [2207.03712](https://arxiv.org/abs/2207.03712) · [10.1103/PhysRevD.106.074019](https://doi.org/10.1103/PhysRevD.106.074019) | yes | ⭑ `2207.03712.pdf` | theory-input | strong-interaction vs charge radius — a named systematic on `gaussian_slope(r_rms_fm)` |
| B7 | `mantysaari2021-azimuthal-correlations` | H. Mäntysaari, K. Roy, F. Salazar, B. Schenke, *Gluon imaging using azimuthal correlations in diffractive scattering at the EIC*, Phys. Rev. D 103 (2021) 094026 | [2011.02464](https://arxiv.org/abs/2011.02464) · [10.1103/PhysRevD.103.094026](https://doi.org/10.1103/PhysRevD.103.094026) | yes | ⭑ `2011.02464.pdf` | theory-input | the *unpolarized* e–V azimuthal modulation: a third cos-type background to `cos2phi_coefficient`, and the size to beat |
| B8 | `mantysaari-schenke2016-subnucleon` | H. Mäntysaari, B. Schenke, *Revealing proton shape fluctuations with incoherent diffraction at high energy*, Phys. Rev. D 94 (2016) 034042 | [1607.01711](https://arxiv.org/abs/1607.01711) · [10.1103/PhysRevD.94.034042](https://doi.org/10.1103/PhysRevD.94.034042) | yes | ⭑ `1607.01711.pdf` | theory-input | the `subnucleondiffraction` code `write_snd_configs()` writes configurations for |
| B9 | `tu2020-src-incoherent-jpsi-tagging` | Z. Tu, A. Jentsch, M. Baker, L. Zheng, J.-H. Lee, R. Venugopalan, O. Hen, D. Higinbotham, E.-C. Aschenauer, T. Ullrich, *Probing short-range correlations in the deuteron via incoherent diffractive J/ψ production with spectator tagging at the EIC*, Phys. Lett. B 811 (2020) 135877 | [2005.14706](https://arxiv.org/abs/2005.14706) · [10.1016/j.physletb.2020.135877](https://doi.org/10.1016/j.physletb.2020.135877) | yes | ⭑ `2005.14706.pdf` | T4 | the incoherent partner of the coherent channel: what a *tagged* diffractive VM event looks like, and the BeAGLE precedent for `CoherentSampler` + spectator tagging |
| **C. The generators the channel is scaled against** |
| C1 | `lomnitz-klein2019-estarlight` | M. Lomnitz, S. Klein, *Exclusive vector meson production at an electron-ion collider*, Phys. Rev. C 99 (2019) 015203 | [1803.06420](https://arxiv.org/abs/1803.06420) · [10.1103/PhysRevC.99.015203](https://doi.org/10.1103/PhysRevC.99.015203) | yes | ⭑ `1803.06420.pdf` | T4 | `estarlight_li6_coherent()`, `estarlight_li6_q2_floors()` — the code every O5 rate comes from |
| C2 | `klein2017-starlight` | S. R. Klein, J. Nystrand, J. Seger, Y. Gorbunov, J. Butterworth, *STARlight: a Monte Carlo simulation program for ultra-peripheral collisions*, Comput. Phys. Commun. 212 (2017) 258 | [1607.03838](https://arxiv.org/abs/1607.03838) · [10.1016/j.cpc.2016.10.016](https://doi.org/10.1016/j.cpc.2016.10.016) | yes | ⭑ `1607.03838.pdf` | T4 | the **Z ≤ 6 Gaussian form factor** branch — the published basis for "⁶Li runs today", and for `slope_b`'s Gaussian model |
| C3 | `toll-ullrich2013-dipole-ea` | T. Toll, T. Ullrich, *Exclusive diffractive processes in electron-ion collisions*, Phys. Rev. C 87 (2013) 024913 | [1211.3048](https://arxiv.org/abs/1211.3048) · [10.1103/PhysRevC.87.024913](https://doi.org/10.1103/PhysRevC.87.024913) | yes | ⭑ `1211.3048.pdf` | T4 | the method behind Sartre; F(b) from coherent dσ/dt — the imaging claim `slope_b` stands in for |
| C4 | `toll-ullrich2014-sartre` | T. Toll, T. Ullrich, *The dipole model Monte Carlo generator Sartre 1*, Comput. Phys. Commun. 185 (2014) 1835 | [1307.8059](https://arxiv.org/abs/1307.8059) · [10.1016/j.cpc.2014.03.010](https://doi.org/10.1016/j.cpc.2014.03.010) | yes | ⭑ `1307.8059.pdf` | T4 | the manual whose Woods–Saxon-only nucleus model is *why* `physics_literature.md` Item 7 rules Sartre out for a tensor axis |
| C5 | `toll2026-sartre-fast-tables` | T. Toll, D. Ghosh, A. Srivastav, *Efficient calculation of exclusive diffractive cross sections at the EIC and LHeC with the Sartre event generator*, arXiv:2606.14633 (June 2026) | [2606.14633](https://arxiv.org/abs/2606.14633) | yes | ⭑ `2606.14633.pdf` | T4 | the 3–4 orders-of-magnitude table speedup — the cost half of the "Sartre route" decision |
| **D. Dipole-model inputs a ⁶Li amplitude would need** |
| D1 | `kowalski-teaney2003-ipsat` | H. Kowalski, D. Teaney, *An impact parameter dipole saturation model*, Phys. Rev. D 68 (2003) 114005 | [hep-ph/0304189](https://arxiv.org/abs/hep-ph/0304189) · [10.1103/PhysRevD.68.114005](https://doi.org/10.1103/PhysRevD.68.114005) | yes | ⭑ `hep-ph_0304189.pdf` | theory-input | IPsat/bSat itself — Sartre's `bSat`, and the b-dependence any ⁶Li amplitude inherits |
| D2 | `kowalski-motyka-watt2006-bsat-exclusive` | H. Kowalski, L. Motyka, G. Watt, *Exclusive diffractive processes at HERA within the dipole picture*, Phys. Rev. D 74 (2006) 074016 | [hep-ph/0606272](https://arxiv.org/abs/hep-ph/0606272) · [10.1103/PhysRevD.74.074016](https://doi.org/10.1103/PhysRevD.74.074016) | yes | ⭑ `hep-ph_0606272.pdf` | theory-input | the bSat fit Sartre cites (its ref. [18]); VM wave functions and t-slopes |
| D3 | `watt-kowalski2008-bcgc` | G. Watt, H. Kowalski, *Impact parameter dependent colour glass condensate dipole model*, Phys. Rev. D 78 (2008) 014016 | [0712.2670](https://arxiv.org/abs/0712.2670) · [10.1103/PhysRevD.78.014016](https://doi.org/10.1103/PhysRevD.78.014016) | yes | ⭑ `0712.2670.pdf` | theory-input | Sartre's `DipoleModel_bCGC` — the model-systematic partner of D1/D2 |
| D4 | `rezaeian2013-ipsat-refit` | A. H. Rezaeian, M. Siddikov, M. Van de Klundert, R. Venugopalan, *Analysis of combined HERA data in the impact-parameter dependent saturation model*, Phys. Rev. D 87 (2013) 034002 | [1212.2974](https://arxiv.org/abs/1212.2974) · [10.1103/PhysRevD.87.034002](https://doi.org/10.1103/PhysRevD.87.034002) | yes | ⭑ `1212.2974.pdf` | theory-input | the current IPsat parameters (combined H1+ZEUS); what a new ⁶Li table would be generated with |
| **E. `f0` — the coherent fraction, and gluon imaging** |
| E1 | `frankfurt-guzey-strikman2012-lt-shadowing` | L. Frankfurt, V. Guzey, M. Strikman, *Leading twist nuclear shadowing phenomena in hard processes with nuclei*, Phys. Rept. 512 (2012) 255 | [1106.2091](https://arxiv.org/abs/1106.2091) · [10.1016/j.physrep.2011.12.002](https://doi.org/10.1016/j.physrep.2011.12.002) | yes | ⭑ `1106.2091.pdf` | theory-input | **the missing `f0` calculation class**: nuclear diffractive PDFs, §6.1.1–6.1.3, Eq. (189) and Fig. 69 = probability of diffraction for a nucleus |
| E2 | `frankfurt-strikman-weiss2005-smallx` | L. Frankfurt, M. Strikman, C. Weiss, *Small-x physics: from HERA to LHC and beyond*, Ann. Rev. Nucl. Part. Sci. 55 (2005) 403 | [hep-ph/0507286](https://arxiv.org/abs/hep-ph/0507286) · [10.1146/annurev.nucl.53.041002.110615](https://doi.org/10.1146/annurev.nucl.53.041002.110615) | yes | ⭑ `hep-ph_0507286.pdf` | theory-input | the review the "coherent diffraction images the transverse gluon distribution" framing comes from (a review — cite E1/C3 for numbers) |
| E3 | `caldwell-kowalski2010-jpsi-imaging` | A. Caldwell, H. Kowalski, *Investigating the gluonic structure of nuclei via J/ψ scattering*, Phys. Rev. C 81 (2010) 025203; preprint *The J/ψ way to nuclear structure*, arXiv:0909.1254 | [0909.1254](https://arxiv.org/abs/0909.1254) (preprint) · [10.1103/PhysRevC.81.025203](https://doi.org/10.1103/PhysRevC.81.025203) (journal, **not** OA) | preprint only | ⭑ `0909.1254.pdf` | theory-input | the \|t\| → b Fourier-Bessel inversion that gives `slope_b` its meaning as an imaging observable |
| **F. The photon-polarization `cos 2φ` — the O5 background** |
| F1 | `star2023-polarized-photon-tomography` | STAR Collaboration (M. S. Abdallah et al.), *Tomography of ultra-relativistic nuclei with polarized photon-gluon collisions*, Sci. Adv. 9 (2023) eabq3903 | [2204.01625](https://arxiv.org/abs/2204.01625) · [10.1126/sciadv.abq3903](https://doi.org/10.1126/sciadv.abq3903) | yes (Sci. Adv. is OA) | ⭑ `2204.01625.pdf` | theory-input | the measured size of mechanism (ii); `tensor_flip_plan`'s reason for existing |
| F2 | `xing2020-cos2phi-rho0` | H. Xing, C. Zhang, J. Zhou, Y.-J. Zhou, *The cos 2φ azimuthal asymmetry in ρ⁰ meson production in ultraperipheral heavy ion collisions*, JHEP 10 (2020) 064 | [2006.06206](https://arxiv.org/abs/2006.06206) · [10.1007/JHEP10(2020)064](https://doi.org/10.1007/JHEP10(2020)064) | yes | ⭑ `2006.06206.pdf` | theory-input | the dipole-model calculation of the linear-polarization cos 2φ — the background amplitude `phase_C_numbers.md` §C2.5 bounds |
| F3 | `zha2021-double-slit` | W. Zha, J. D. Brandenburg, L. Ruan, Z. Tang, Z. Xu, *Exploring the double-slit interference with linearly polarized photons*, Phys. Rev. D 103 (2021) 033007 | [2006.12099](https://arxiv.org/abs/2006.12099) · [10.1103/PhysRevD.103.033007](https://doi.org/10.1103/PhysRevD.103.033007) | yes | ⭑ `2006.12099.pdf` | theory-input | the interference reading of the same modulation — the alternative that makes the two-nucleus part of F1 inapplicable to e+A |
| F4 | `hagiwara2021-dipion-tomography` | Y. Hagiwara, C. Zhang, J. Zhou, Y.-J. Zhou, *Probing the gluon tomography in photoproduction of di-pions*, Phys. Rev. D 104 (2021) 094021 | [2106.13466](https://arxiv.org/abs/2106.13466) · [10.1103/PhysRevD.104.094021](https://doi.org/10.1103/PhysRevD.104.094021) | yes | ⭑ `2106.13466.pdf` | theory-input | cos 4φ and elliptic-gluon vs QED-radiation splitting — the higher harmonic the ρ⁰ control has to survive |
| **G. `amp` — the flat gluon-transversity scenario** (shared with the tensor-SF domain) |
| G1 | `detmold-shanahan2016-gluon-transversity` | W. Detmold, P. E. Shanahan, *Gluonic transversity from lattice QCD*, Phys. Rev. D 94 (2016) 014507 [erratum PRD 95 (2017) 079902] | [1606.04505](https://arxiv.org/abs/1606.04505) · [10.1103/PhysRevD.94.014507](https://doi.org/10.1103/PhysRevD.94.014507) | yes | ⭑ `1606.04505.pdf` (fetched by a sibling domain) | theory-input | the lattice half of `CoherentScenario::amp`'s 3e-3 … 1e-2 band |
| G2 | `kumano-song2020-gluon-transversity-dy` | S. Kumano, Qin-Tao Song, *Gluon transversity in polarized proton-deuteron Drell-Yan process*, Phys. Rev. D 101 (2020) 054011 | [1910.12523](https://arxiv.org/abs/1910.12523) · [10.1103/PhysRevD.101.054011](https://doi.org/10.1103/PhysRevD.101.054011) | yes | ⭑ `1910.12523.pdf` (fetched by a sibling domain) | theory-input | the Drell-Yan half of the same band |

---

## 2. What each reference contains, and what it is for

### A. The tensor / polarized question

**A1 — Mäntysaari, Salazar, Schenke, Shen, Zhao, arXiv:2408.13213, PLB 858 (2024) 139053.**
*In the corpus already* (`refs/2408.13213v1.pdf`); repeated here because it is the
anchor every other row is measured against. It computes coherent J/ψ production off a
**polarized deuteron** in the CGC with nucleon-level (u, w) structure, writes the
azimuthal dependence as dσ/dΦ d\|t\| ∝ 1 + 2 Σ aₙ cos nΦ (its Eq. 9), and gives a₂, a₄ per
deuteron spin state at x_P = 1.7·10⁻³. `coherent.hpp`'s `mantysaari_a2_deuteron()` is four
digitized rows of its Fig. 4; `cluster_config.hpp`'s `a2_from_quadrupole` is its Eq. (9)
re-derived with a quadrupole in place of the numerical amplitude, and reproduces its
a₂(±1) to 8 % with no free parameter. **It is still the only polarization-dependent
coherent-diffraction calculation for any nucleus** — see §3 for how that was re-checked.

**A2 — Dalton, Deur, Keith, arXiv:2504.21177, EPJ A 61 (2025) 111.** New to this pass, and
the domain brief's "any paper on tensor-polarized-target coherent diffraction" has exactly
one more answer than the tree knew. Hall D proposes a frozen-spin **tensor-polarized
deuteron target** in the 12 GeV bremsstrahlung beam, with AFP spin manipulation to enhance
the m = 0 population, and its flagship measurement is **coherent ρ photoproduction off the
deuteron**: the abstract states the coherent ρ cross section is sensitive to double
scattering at high momentum transfer and that "in the m = 0 spin state an additional
sensitivity at intermediate momentum transfer opens up". That is the same structure as
`a2_m_state` / `delta_b_m` — an m-dependent modification of the coherent \|t\| distribution —
at Q² = 0 and A = 2, from an experimental group. Two uses: it is the nearest external
precedent for the O5 proposal's *shape* of argument (m-resolved coherent \|t\|), and it is a
live collaboration target for the ⁶Li ask recorded in
`run_2026-09-03/mantysaari_collaboration_draft.md`. Verified: arXiv abstract fetched;
journal reference and DOI from INSPIRE recid 2916899.

**A3 — Good & Walker, Phys. Rev. 120 (1960) 1857.** The origin of the eigenstate picture:
a projectile is a superposition of states with different absorption, the coherent (elastic)
amplitude is the average, and the difference produces diffractive dissociation. Verified on
INSPIRE (recid 9375, 750 citations; abstract fetched and read). **No arXiv, and INSPIRE
holds no fulltext for it** — the download link is the APS DOI page,
<https://doi.org/10.1103/PhysRev.120.1857>, through an institutional subscription. It is
cited by C4 (Sartre's §"According to the Good-Walker picture") and is the formal
justification for what `ClusterConfigSampler` / `write_snd_configs` do: produce an ensemble
of α+d configurations whose amplitude-average is the coherent cross section.

**A4 — Miettinen & Pumplin, Phys. Rev. D 18 (1978) 1696.** The quantitative companion to A3:
coherent = ⟨A⟩², total diffractive = ⟨A²⟩, so incoherent = variance over the eigenstate
ensemble. Every code in §B–C implements this sentence, and so does `ClusterConfigSet`.
Verified on INSPIRE (recid 129259, 271 citations). Paywalled at
<https://doi.org/10.1103/PhysRevD.18.1696>; INSPIRE serves a free scan at
<https://inspirehep.net/files/8c3342725e4cd03851b56f8078d4e1fa> (checked this pass:
HTTP 200, `application/pdf`, 1.37 MB). **Not downloaded** — an INSPIRE-hosted scan of an
APS article is not unambiguously an open-access copy, and this pass only fetched arXiv
e-prints. Use whichever of the two links your institution's rules allow.

**A5 — Rachek et al., Few Body Syst. 58 (2017) 29.** Listed as precedent, not as an input.
VEPP-3 measured the **tensor analyzing power T₂₀ in coherent π⁰ photoproduction on a
tensor-polarized deuteron**, 200 < E_γ < 500 MeV, 90° < θ_π^cm < 140° (INSPIRE recid
1508664; abstract fetched). It is the only *measured* tensor observable of a coherent
production reaction on a nucleus anywhere — but at 0.2–0.5 GeV it is a resonance-region
nuclear-structure measurement, not diffraction, and it does not constrain any number in
`coherent.hpp`. Its value is rhetorical: "a tensor observable in coherent production has
been measured before" is true, and the O5 write-up should say where. No arXiv; Springer
paywall at <https://doi.org/10.1007/s00601-016-1191-0>. (The same group's earlier JETP Lett.
89 (2009) 432, <https://doi.org/10.1134/S0021364009090021>, is the companion.)

### B. Light nuclei with realistic structure

**B1 — Guzey, Rinaldi, Scopetta, Strikman, Viviani, arXiv:2202.12200, PRL 129 (2022) 242503.**
The single most useful new reference in this domain. Coherent **J/ψ electroproduction on
⁴He and ³He** at EIC kinematics, computed in a Gribov–Glauber multiple-scattering expansion
in which the k-nucleon terms are separated explicitly, with the one-, two- and three-body
nuclear form factors Φ₁, Φ₂, Φ₃ built from **realistic AV18-based few-body wave functions**
(PDF text, p. 4: "the nuclear ffs Φ₁, Φ₂, and Φ₃ have been calculated using a realistic wave
function … using the AV18 nucleon-nucleon" interaction). It predicts a shift of the coherent
\|t\| distribution toward smaller \|t\|, and shows the ⁴He diffractive minimum is *not* at the
charge-form-factor minimum (p. 2: the cross section "does not present a minimum at −t ≃ 0.6
GeV², where the ⁴He charge ff has a minimum"). Three uses: (i) it is the existence proof
that a light-nucleus coherent amplitude can be built from the *same class of wave functions*
the tree already carries (ANL AV18 VMC), which is the technical core of the ⁶Li ask;
(ii) its inputs are all nucleon-level and published — the γ*p → J/ψ p slope
B₀(x) = 4.5 GeV⁻², the γ*p → Xp diffractive slope B ≃ 6 GeV⁻² that together with the
proton's gluon-channel diffractive probability fixes the two-nucleon effective cross
section σ₂(x = 10⁻³) = 25 mb ± 15 % (and σ₃ = 30–50 mb, to which light nuclei are shown
to be insensitive at small \|t\|) — none of which is `slope_b` = 50 GeV⁻², a *nuclear*
\|F(t)\|² slope B = ⟨r²⟩/3, and the two must not be conflated; (iii) the "charge-ff minimum ≠ diffractive minimum"
result is a direct warning against the tree's Gaussian `gaussian_slope(r_rms_fm)` shortcut
at \|t\| beyond the first lobe — which is one more reason `COHERENT_T_MAX_DEFAULT` = 0.2
is the right ceiling.

**B2 — Mäntysaari, Roch, Schenke, Shen, Zhao, arXiv:2605.00454, PRD 114 (2026) 014068.**
The paper `physics_literature.md` Item 7's 2026 annotation flagged as retiring the
"lightest published is Ca" note, now on disk and checked. Its headline is O+O and Ne+Ne UPCs
at the LHC, but the model section (PDF p. 4) reads: "We employ several different models to
generate initial-state nucleon configurations, including Woods-Saxon … ab-initio Variational
Monte Carlo (VMC) calculations, Projected Generator Coordinate Method (PGCM) calculations
with two different sampling methods ('clustering' and 'uniform'), as well as Nuclear Lattice
Effective Field Theory (NLEFT) configurations and Green's Function Monte Carlo (GFMC)
calculations using the AV18+UIX model interaction. … **The GFMC model is used to describe
He-3 and He-4 nuclei.**" Its suppression-factor figure runs over **A ≥ 3**. So: ab-initio
nucleon configurations of exactly the ANL kind, inside a Good–Walker CGC amplitude, down to
A = 3 — the machinery of the ⁶Li ask, minus the polarization. Two concrete uses for the
collaboration draft: it shows the group already ingests VMC/GFMC configuration files (so
"supply an α+d ⁶Li configuration file" is a small ask, not a new interface), and its
"PGCM clustering" vs "uniform" comparison is the closest published measurement of how much
a *clustered* sampling changes a coherent \|t\| distribution — the systematic
`cluster_config.hpp`'s diagonal truncation would inherit. Caveat, stated: it is unpolarized,
UPC-centric, and its A = 3, 4 results are GFMC-spherical, so it carries no tensor axis.

**B3 — Chang et al., arXiv:2511.05638, PRD 113 (2026) 032018.** *In the corpus already.*
The source of `COHERENT_JPSI_EFF_IR8_LI7` = 0.1775 and of the two efficiency tables
(`Chang26EffEnergyRow`, `Chang26SpeciesEffRow`) that carry O5's whole detection chain, and
of the 0.1 < Q² < 100 GeV² acceptance window that `estarlight_li6_q2_floors()` exists to
step outside of. Nothing new to fetch; listed because half of §11.3's open items are about
how far its ⁷Li-at-18×118 number can be transported to ⁶Li-at-10×99.5.

**B4 — Mondal, Kumar, Sarkar, arXiv:2608.23445.** Within a saturation framework, coherent
VM production is shown to be sensitive to nuclear **shell** structure through self-consistent
QMC single-particle densities, "enhancing the secondary diffractive lobes in light nuclei"
(abstract, fetched). Mean-field, unpolarized, no clusters. Use: it bounds how much of a
light-nucleus coherent \|t\| shape is *not* captured by a one-parameter Gaussian, independently
of B1's few-body route — i.e. it is the second leg of the `slope_b` model systematic. Not
yet journal-published (INSPIRE recid 3195389, no DOI as of this pass).

**B5 / B6 — Mäntysaari et al. arXiv:2303.04866 (PRL 131 062301) and arXiv:2207.03712
(PRD 106 074019).** The deformation machinery. B5 varies the deformation of a U target and
of ²⁰Ne vs ¹⁶O and finds incoherent production the deformation-sensitive channel, with
different multipoles populating different \|t\| ranges; B6 fits ALICE UPC J/ψ spectra and finds
a preference for a **strong-interaction radius larger than the charge radius**. Both are
unpolarized-geometry analogues of what `delta_b_m` does with a spin-aligned quadrupole. B6 is
the citation for a systematic the tree currently does not carry: `gaussian_slope()` is fed a
*charge* radius, and B6 says the gluonic radius is the larger one — a shift in `slope_b`, in
the direction that helps (larger B ⇒ smaller ⟨\|t\|⟩ ⇒ smaller a₂ at fixed eps_b0, so it must
be banded, not ignored).

**B7 — Mäntysaari, Roy, Salazar, Schenke, arXiv:2011.02464, PRD 103 (2021) 094026.**
Coherent diffractive photon and VM production in the CGC with MV + JIMWLK configurations,
computing the **electron–vector-meson azimuthal correlation** differentially in Q² and \|t\|:
v₁,₂ ≈ 2–10 % for e+p, and **< 0.1 % for e+Au** (abstract, fetched). This is a *third*
cos-modulation mechanism, distinct from both of the two the tree tracks: it is unpolarized,
comes from spatial gluon correlations, and is referred to the **lepton plane**. Its
usefulness is the number: the nuclear suppression from p to Au takes it below 0.1 %, i.e.
below the 0.25 % information-weighted tensor modulation O5 lives on, but a ⁶Li-sized target
is much nearer the proton than gold and this background has never been evaluated there.
That is an honest, citable open item for `cos2phi_coefficient`, and it is not in the tree's
current background list (which has only the photon-polarization term).

**B8 — Mäntysaari & Schenke, arXiv:1607.01711, PRD 94 (2016) 034042.** The paper behind
the `subnucleondiffraction` code that `write_snd_configs()` targets: coherent vs incoherent
J/ψ at HERA from a fluctuating-hot-spot proton in the IPsat/CGC framework. It is the
methods citation for the file format and for the Good–Walker averaging the exported ⁶Li
configurations are meant to enter.

**B9 — Tu, Jentsch, Baker, Zheng, Lee, Venugopalan, Hen, Higinbotham, Aschenauer, Ullrich,
arXiv:2005.14706, PLB 811 (2020) 135877.** The e+d companion to B3: BeAGLE events for
**incoherent** diffractive J/ψ production with **spectator-nucleon tagging**, using the tagged
nucleon four-momentum in the plane-wave impulse approximation to reach the high-momentum part
of the deuteron wave function (abstract, fetched). It is not a coherent calculation, and it is
listed for two reasons: it is the published precedent for combining a diffractive VM final state
with far-forward spectator tagging — the exact event topology `CoherentSampler` + `recoil_lab()`
produce, one channel over — and it is the reference for why the **incoherent** background to a
coherent ⁶Li sample is not negligible and has to be vetoed rather than assumed away (the tree's
`veto_table()` / `fragment_rigidity()` machinery). Its A = 2 target means none of its numbers
transfer to ⁶Li directly.

### C. The generators

**C1 — Lomnitz & Klein, arXiv:1803.06420, PRC 99 (2019) 015203 (eSTARlight).** The generator
that produced every number in `estarlight_li6_coherent()` and `estarlight_li6_q2_floors()`,
and hence the f₀ lower bound and the whole O5 rate. Exclusive ρ, φ, J/ψ, ψ′, Υ
electroproduction at an EIC with acceptance studies vs detector rapidity coverage.

**C2 — Klein, Nystrand, Seger, Gorbunov, Butterworth, arXiv:1607.03838, CPC 212 (2017) 258
(STARlight).** eSTARlight's parent, and the *published* home of the statement the ⁶Li plan
rests on: PDF text, §2 — "**For light nuclei, with Z ≤ 6, a Gaussian mass distribution is
used**", and again under form factors, "For light nuclei, with Z ≤ 6, a Gaussian form factor"
is used, with heavy ions instead a hard sphere convolved with a Yukawa. ⁶Li has Z = 3, so
the claim in `physics_literature.md` Item 7 — that adding ⁶Li is a table entry inside an
existing branch rather than new physics code — now has a citation, and `slope_b`'s Gaussian
\|F(t)\|² is the same functional form the generator itself would use.

**C3 / C4 / C5 — Sartre.** C3 (Toll–Ullrich, PRC 87 (2013) 024913) is the method: dipole
amplitudes for eA exclusive VM/DVCS, coherent dσ/dt as the route to the gluon source
distribution F(b). C4 (CPC 185 (2014) 1835) is the generator manual; its §2 names the
implemented models (`bSat`, `bNonSat`, `bCGC`, referencing D1 and D2), states the Good–Walker
coherent/incoherent split explicitly, and describes nuclei only through "the Woods-Saxon
function projected onto the transverse plane" with configurations Ω = {b₁…b_A} — which is
precisely why `physics_literature.md` Item 7's verdict (Sartre cannot carry a tensor axis;
nucleon sampling is spherically symmetric by construction) is right and should cite C4 §3
rather than the tables README alone. C5 (arXiv:2606.14633) removes the lookup-table
bottleneck by 3–4 orders of magnitude ("a few CPU-year" → "a few hours", abstract) but the
code is unreleased; it changes the cost of the Sartre route, not its feasibility for a
polarized target.

### D. Dipole-model inputs

**D1–D4.** If a ⁶Li amplitude is ever built — in Sartre, in `subnucleondiffraction`, or in
the Mäntysaari-group code — these four fix the dipole cross section it is built from.
**D1** (Kowalski–Teaney) introduces the impact-parameter dipole saturation model (IPsat/bSat)
and its linearized bNonSat; **D2** (Kowalski–Motyka–Watt) is the exclusive-process fit,
including the vector-meson wave functions and t-slopes, and is Sartre's own citation [18];
**D3** (Watt–Kowalski) is the bCGC alternative, i.e. the model systematic; **D4** (Rezaeian,
Siddikov, Van de Klundert, Venugopalan) is the IPsat refit to the *combined* H1+ZEUS reduced
cross sections and is the parameter set a table generated today would use. None of them
touches polarization; they are listed because the O5 collaboration draft asks for a
calculation, and a calculation has to say which dipole model and which fit.

### E. `f0` and gluon imaging

**E1 — Frankfurt, Guzey, Strikman, arXiv:1106.2091, Phys. Rept. 512 (2012) 255.** The
answer to a gap `coherent.hpp` states in prose and then leaves open: "Determining f0 needs a
coherent diffractive-DIS calculation (a coherent-A analogue of the H1/ZEUS diffractive PDFs),
which exists in neither eSTARlight nor Sartre." It exists — in this report. §6.1 builds
**nuclear diffractive structure functions and diffractive PDFs** from the leading-twist
Gribov–Glauber generalization of the H1/ZEUS diffractive PDFs; §6.1.1 is coherent diffraction,
§6.1.2 incoherent, §6.1.3 the numerical predictions, where **Eq. (189)** defines the
probability of diffraction P_diff^j = ∫_x^0.1 dx_IP β f_j^{D(3)}(β,Q²,x_IP) / x f_j(x,Q²) per
parton flavour and **Fig. 69** plots it vs x at Q² = 4 GeV² for ⁴⁰Ca and ²⁰⁸Pb, gluon and
quark channels, against the free proton (verified by reading the PDF text around §6.1.3;
the report also records that the gluon-channel probability is about twice the quark-channel
one, and that for nuclei it is at or below the proton value depending on the FGS10 L/H model).
That is `f0`'s definition, computed — modulo two honest caveats that must ride with any use:
the lightest nucleus treated there is ⁴⁰Ca (deuterium is handled separately in §4), so ⁶Li
needs an A-extrapolation, and P_diff is *total* diffraction, with §6.1.1–6.1.2 supplying the
coherent/incoherent split. This is a several-day task with a citable answer at the end of it,
and it replaces the scenario band {0.02, 0.08} that `CoherentScenario::f0` currently ships.

**E2 — Frankfurt, Strikman, Weiss, hep-ph/0507286, Ann. Rev. Nucl. Part. Sci. 55 (2005) 403.**
The review the brief asked for: exclusive and diffractive processes in DIS in a unified
parton/dipole treatment, colour transparency, nuclear shadowing, the black-disk limit, and
the three-dimensional gluon structure that coherent \|t\| distributions measure. Use it for
framing and for the physical argument that coherent diffraction images transverse gluon
structure; take *numbers* from E1 or C3, per the project rule that a primary source beats a
review.

**E3 — Caldwell & Kowalski.** The proposal to measure the nuclear gluon distribution by
elastic J/ψ scattering and invert the coherent \|t\| distribution to F(b). The refereed version
is PRC 81 (2010) 025203 (INSPIRE recid 850035) and has **no arXiv e-print**; the freely
available copy is the earlier preprint arXiv:0909.1254, *The J/ψ way to nuclear structure*
(abstract fetched and read: "We propose to investigate the properties of nuclear matter by
measuring the elastic scattering of J/ψ on nuclei with high precision … at the future ENC,
EIC or LHeC"), which is what was downloaded. If the journal version is wanted:
<https://doi.org/10.1103/PhysRevC.81.025203>. It gives `slope_b` its interpretation — B is
not a fit parameter but the second moment of the gluon transverse profile — which is exactly
the argument `cluster_config.hpp`'s `eps_b0_equivalent()` makes when it converts a
configuration ensemble into a slope.

### F. The photon-polarization `cos 2φ`

These four are the background side of O5, i.e. mechanism (ii) in
`physics_literature.md` §7's warning and the thing `tensor_flip_plan`'s two-fill difference
is designed to cancel.

**F1 — STAR, arXiv:2204.01625, Sci. Adv. 9 (2023) eabq3903.** The measurement. Quasi-real
photons from a relativistic nucleus are 100 % linearly polarized along their transverse
momentum; diffractive ρ⁰ → π⁺π⁻ photoproduction then shows a second-order azimuthal
modulation, fitted in the paper as f(φ) = 1 + A cos(2φ) and reported as 2⟨cos 2φ⟩ (both
strings verified in the PDF text), from which strong-interaction radii of 6.53 ± 0.06 fm
(¹⁹⁷Au) and 7.29 ± 0.08 fm (²³⁸U) are extracted — larger than the charge radii, the same
conclusion as B6. Published in Science Advances, which is fully open access. It is the
empirical anchor for "the background is real and this big"; note for honesty that its
geometry is A+A UPC with two possible photon emitters, so the interference part of its signal
has no e+A analogue (see F3).

**F2 — Xing, Zhang, Zhou, Zhou, arXiv:2006.06206, JHEP 10 (2020) 064.** The dipole-model
calculation of the same effect: joint impact-parameter and transverse-momentum dependent
cross sections for ρ⁰ photoproduction in UPCs, unpolarized cross section plus the cos 2φ
correlation from linearly polarized photons, compared to STAR. This is the reference to cite
for the *amplitude* of the background modulation when bounding it against the tensor term —
`phase_C_numbers.md` §C2.5's "background amplitude above 1.9 (J/ψ), 0.13 (φ), 0.017 (ρ⁰)"
thresholds are only meaningful against a published expectation, and this is it.

**F3 — Zha, Brandenburg, Ruan, Tang, Xu, arXiv:2006.12099, PRD 103 (2021) 033007.** Reads the
same STAR modulation as a **double-slit interference** between the two nuclei acting as
photon emitter and target, predicting a periodic oscillation with the VM transverse momentum.
Its importance for LiPolGen is negative and load-bearing: the part of F1's signal that comes
from two-emitter interference *cannot* occur in e+A, where there is exactly one photon source.
So F1's measured size is an upper bound, not an estimate, for the e+⁶Li background — and the
tree should not import F1's amplitude directly into `cos2phi_coefficient`'s systematic.

**F4 — Hagiwara, Zhang, Zhou, Zhou, arXiv:2106.13466, PRD 104 (2021) 094021.** The cos **4**φ
companion: STAR's di-pion cos 4φ asymmetry gets contributions from both an elliptic gluon
Wigner distribution and final-state soft-photon radiation, and QED alone underestimates the
data. Two uses: it is the reference for the higher harmonic a ρ⁰-channel control would have
to fit simultaneously (the tree's ρ⁰ path already needs \|t\| binning), and its elliptic-gluon
term is the same physics as B7's v₂ — a target-structure cos 2φ/cos 4φ that is not the tensor
one.

### G. `amp` — the flat gluon-transversity scenario

`CoherentScenario::amp` (band 3e-3 … 1e-2) is documented as "bounded by lattice + Drell-Yan
estimates" with no citation in the header. **G1** (Detmold & Shanahan, PRD 94 (2016) 014507,
erratum PRD 95 079902) is the lattice calculation of the gluon transversity matrix element in
a spin-1 system; **G2** (Kumano & Song, PRD 101 (2020) 054011) is the polarized
proton–deuteron Drell-Yan route to the same quantity. Both were fetched into the corpus by a
sibling domain pass during this run (`1606.04505.pdf`, `1910.12523.pdf`, both 2026-09-16) and
are listed here only to close the citation loop for `amp`; the tensor-structure-function
domain owns them. Related: `PHYSICS_CHANNELS.md` records that the repository carries **no**
literature reference for the double-helicity-flip Δ at all and that `toy_delta_gluon` is an
unsourced ansatz — the same gap, in the inclusive channel.

---

## 3. What does **not** exist — re-checked, not assumed

The tree asserts (`physics_literature.md` §7, and the 2026-09 annotation) that
arXiv:2408.13213 remains the only polarization-dependent coherent-diffraction calculation
for any nucleus, and that nothing exists for ⁶Li or for any α-cluster nucleus in coherent
diffraction. This pass re-ran the check independently:

1. **Forward citations of 2408.13213** (INSPIRE `refersto recid 2821172`, fetched this pass):
   the citing set is heavy-ion structure and deformation work — *Probing nuclear structure of
   heavy ions at the LHC* (2409.19064), *Energy dependence of the deformed nuclear structure
   at small-x* (2411.14934), the symmetry-restoration papers (2509.09549, 2606.02412),
   *Quantum stress and torsion distributions in the deuteron* (2602.18298), B2 (2605.00454),
   A2 (2504.21177) and a Pb+d↑ hydrodynamics proceedings. **None extends a polarized coherent
   amplitude to A > 2, and none is a ⁶Li calculation.**
2. **A = 3, 4 coherent diffraction now exists but is unpolarized** — B1 (few-body AV18 wave
   functions) and B2 (GFMC configurations). The "lightest published is Ca" note is dead at the
   *structure* level as well as the rate level; it should be retired wherever it still appears.
3. **α-clustering in coherent diffraction** appears only as a *sampling option* in B2
   ("PGCM clustering") and as shell/mean-field structure in B4. No coherent-diffraction
   calculation uses an explicit α+d or α-cluster wave function, for ⁶Li or anything else.
4. **Tensor-polarized target + coherent production** exists exactly twice outside A1: A2
   (proposed, A = 2, real photon, ρ) and A5 (measured, A = 2, real photon, π⁰, 0.2–0.5 GeV).
   Neither is diffractive and neither is A > 2.

So the project's first-mover claim stands, and is now sourced at four points rather than
asserted.

---

## 4. What this changes in the tree (suggestions, not edits)

Nothing in `coherent.hpp` or `cluster_config.hpp` was touched by this pass. Four items are
worth an author decision:

- **`f0` has a calculation class after all** (E1 §6.1). The header's "exists in neither
  eSTARlight nor Sartre" is true of *generators* and should be narrowed to that; the sentence
  as written reads as "does not exist".
- **`slope_b` needs a second systematic**: B6's gluonic-radius-larger-than-charge-radius
  result, and B1/B4's "the coherent minimum is not the charge-form-factor minimum". Both push
  in the direction of using the Gaussian only below the first lobe — which
  `COHERENT_T_MAX_DEFAULT` = 0.2 already does, for a different stated reason.
- **The background list for `cos2phi_coefficient` is one mechanism short**: B7's unpolarized
  e–V azimuthal correlation (v₁,₂ ≈ 2–10 % on a proton, < 0.1 % on gold) is referred to the
  lepton plane and has never been evaluated for a light nucleus. `tensor_flip_plan`'s P_zz-odd
  argument removes it as well — it is P_zz-even like the photon term — but it should be named.
- **F1's amplitude must not be imported directly** as the e+⁶Li background scale: F3 shows a
  two-emitter interference component with no e+A analogue.

---

## 5. Download pointers for the non-OA items

| Item | Where the user should download it |
|---|---|
| A3 Good–Walker PR 120 (1960) 1857 | <https://doi.org/10.1103/PhysRev.120.1857> (APS; no arXiv, no INSPIRE fulltext) |
| A4 Miettinen–Pumplin PRD 18 (1978) 1696 | <https://doi.org/10.1103/PhysRevD.18.1696> (APS); free scan served by INSPIRE at <https://inspirehep.net/files/8c3342725e4cd03851b56f8078d4e1fa> |
| A5 Rachek et al. Few Body Syst. 58 (2017) 29 | <https://doi.org/10.1007/s00601-016-1191-0> (Springer) |
| E3 Caldwell–Kowalski PRC 81 (2010) 025203 | <https://doi.org/10.1103/PhysRevC.81.025203> (APS) — or use the OA preprint `refs/0909.1254.pdf` |

All other rows in Table 1 are arXiv e-prints and are on disk.
