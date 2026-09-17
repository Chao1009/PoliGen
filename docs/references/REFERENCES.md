<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# REFERENCES — LiPolGen's annotated bibliography

**Assembled 2026-09-16** from the reference-research series in this directory:
`00_corpus.md` (inventory of what was already held) and the six domain files
`01_tensor-b1.md`, `01_li-nuclear-data.md`, `01_radiative.md`,
`01_coherent-diffraction.md`, `01_generators-chain.md`,
`01_nucleon-deuteron-data.md`. Every identity claim below (authors, title,
journal, volume, year, DOI, arXiv id) and every statement about what a paper
contains was verified in that series against the arXiv API
(`export.arxiv.org/api/query`), the INSPIRE literature API
(`inspirehep.net/api/literature`), Crossref, or the text of the PDF itself.
Items that were **not** verified are collected in §9 and are labelled
`unverified`; nothing here asserts a paper's content that was not read from its
abstract or its text.

**Organisation.** Grouped by the `BENCHMARK_PLAN.md` tier scheme, plus
**theory inputs** (formalism the code transcribes rather than data it is
compared against) and **documents / standards**:
§1 the purchase list · §2 **T1** nucleon · §3 **T2** deuteron and b₁ ·
§4 **T3** nuclear inputs (⁶Li / ⁷Li / ⁴He) · §5 **T4** generators ·
§6 theory inputs · §7 **T5/T6** chain and software ·
§8 documents, standards and data archives · §9 unverified · §10 counts.

**Status vocabulary.**

| status | meaning |
|---|---|
| `held` | a PDF is on disk under `/home/cpeng/Projects/polli/PolarizedLithiumSim/refs/` (file named in the entry) |
| `downloaded 2026-09-16` | fetched by this reference stage into that same directory |
| `PLEASE DOWNLOAD` | no copy obtainable from this environment — the exact publisher/DOI link to use is given |
| `free, not fetched` | open access and reachable, but no PDF was committed — a data archive, a free report, or a document too large to commit; the link is given |
| `unverified` | named in the LiPolGen tree, **not** checked against a primary index in this pass and not on disk |

Counts are in §10.

---

## 1. PLEASE DOWNLOAD — act on these first

These are the references the user has to obtain. Sections **1a–1c** are ordered
by how much they unblock: 1a closes a named gate or supplies a number the code
already ships, 1b bounds something, 1c is completeness. Every link is the
publisher/DOI page — no paywall was approached in assembling this list.

The list is numbered **1–58**; **56** of those still need obtaining. Item **43**
was **RESOLVED** in the final verification pass — a free preprint scan turned up
on INSPIRE and was fetched — and item **27** was found in the reconciliation
pass (§10.4) to have a free INSPIRE-hosted JINR preprint, verified but not
committed, so it is `free, not fetched`; every other item here was re-queried
against the INSPIRE literature API and the arXiv API and confirmed to have
**no arXiv e-print**; where INSPIRE serves a fulltext, the entry now says what
that file actually is. The method and its findings are in §10.2 and §10.4.

### 1a. Highest value — each closes a named gate, number or open item

1. **Close, Kumano**, *A sum rule for the spin-dependent structure function b₁(x) for spin-one hadrons*, Phys. Rev. D **42** (1990) 2377 — <https://doi.org/10.1103/PhysRevD.42.2377>. **PLEASE DOWNLOAD**. The defining paper of `close_kumano_integral` and in-tree check **E-4**: its abstract gives ∫dx b₁(x) = −(5/3) lim_{t→0}(t/4M²)F_Q(t) = 0 for an isoscalar target *provided the quark/antiquark sea is unpolarized*, and shows how a polarized sea modifies it. The gate exists; the paper behind it does not. (Free secondary with the antiquark correction: `refs/1005.4524.pdf` §I and Eq. 9.)
2. **Akushevich, Shumeiko**, *Radiative effects in deep inelastic scattering of polarized leptons by polarized light nuclei*, J. Phys. G **20** (1994) 513 — <https://doi.org/10.1088/0954-3899/20/4/001>. **PLEASE DOWNLOAD**. POLRAD's own theory paper **for nuclear targets** (H, D, ³He), in covariant variables, with the transverse-target case and a discussion of correct-the-data vs correct-the-model. It is the only published place the A-dependence of a polarized radiative correction is worked out, i.e. the source `RcOptions::a_transfer_frac` (the A = 2 → A = 6 transfer of the tensor band, default 0.0) currently has none for. Highest-value paywalled item in the radiative domain.
3. **Hoodbhoy, Jaffe, Manohar**, *Novel effects in deep inelastic scattering from spin-one hadrons*, Nucl. Phys. B **312** (1989) 571 — <https://doi.org/10.1016/0550-3213(89)90572-5>. **PLEASE DOWNLOAD** (no arXiv, no preprint scan; INSPIRE recid 262935). Defines b₁…b₄ and Δ and gives the (1, −2, 1) helicity pattern that `TensorSF`'s rank-2 geometry implements; cited across `CONVENTIONS.md`, `PHYSICS_CHANNELS.md` `[HJM89]` and `BENCHMARK_PLAN.md`.
4. **Jaffe, Manohar**, *Nuclear gluonometry*, Phys. Lett. B **223** (1989) 218 — <https://doi.org/10.1016/0370-2693(89)90242-6>. **PLEASE DOWNLOAD** (INSPIRE recid 277858, no arXiv). The paper that names Δ an "exotic glue" observable with no nucleonic counterpart — the conceptual source of `toy_delta_gluon` and of the coherent cos 2φ term. **Not** covered by the corpus's folded 1989 Jaffe–Manohar entry. Free secondary that quotes it: `refs/1803.11206.pdf`.
5. **Jaffe, Manohar**, *Deep inelastic scattering from arbitrary spin targets*, Nucl. Phys. B **321** (1989) 343 — <https://doi.org/10.1016/0550-3213(89)90347-7>. **PLEASE DOWNLOAD** (INSPIRE recid 266865, no arXiv). The arbitrary-J generalisation, and the formal home of the J = 3/2 basis LiPolGen's ⁷Li rank-2 block would need.
6. **Sather, Schmidt**, *Size and scaling of the double-helicity-flip hadronic structure function*, Phys. Rev. D **42** (1990) 1424 — <https://doi.org/10.1103/PhysRevD.42.1424>. **PLEASE DOWNLOAD** (no arXiv). The actual source of `C_BAG` = −0.012, a constant the generator ships and cannot currently re-derive.
7. **Suelzle, Yearian, Crannell**, *Elastic Electron Scattering from Li⁶ and Li⁷*, Phys. Rev. **162** (1967) 992 — <https://doi.org/10.1103/PhysRev.162.992>. **PLEASE DOWNLOAD**. The origin of both lithium charge form factors. Its abstract (read via OpenAlex in this pass) says the ⁷Li form factor *is* described by a simple harmonic-oscillator shell model (r_rms 2.39 ± 0.03 fm) and that the ⁶Li one **cannot** be fitted by either the simple or the modified harmonic-well model. That is a published statement about the functional form `HoSpin1FF` ships for ⁶Li, and it is the paper behind `rc.hpp`'s "TO BE REFIT". *Serves:* `LI6_FF_HO_A_FM`, `LI6_FF_HO_ALPHA`, `HoSpin1FF::for_ion(LI7)`, open item Q1.
8. **Rand, Frosch, Yearian**, *Elastic Electron Scattering from the Magnetic Multipole Distributions of Li⁶, Li⁷, Be⁹, B¹⁰, B¹¹ and N¹⁴*, Phys. Rev. **144** (1966) 859 — <https://doi.org/10.1103/PhysRev.144.859>. **PLEASE DOWNLOAD**. `BENCHMARK_PLAN.md` §4 item 6, and the only measured electromagnetic signature of α+d clustering in the magnetic sector. `rc.hpp` builds `F_m(0)` from μ and then has nothing measured to constrain its q dependence.
9. **de Jager, de Vries, de Vries**, *Nuclear charge- and magnetization-density-distribution parameters from elastic electron scattering*, At. Data Nucl. Data Tables **14** (1974) 479 — <https://doi.org/10.1016/S0092-640X(74)80002-1>. **PLEASE DOWNLOAD**. The compiled form of entry 8 — and the **only** free-standing magnetization compilation, since the UVa machine-readable archive is charge-only (verified: `FB_data.dat`, `SOG_data.dat`, `tabletr.txt` carry no magnetization rows).
10. **Ent** *et al.*, *The (e,e′d) reaction on ⁴He, ⁶Li, and ¹²C*, Nucl. Phys. **A578** (1994) 93 — <https://doi.org/10.1016/0375-9474(94)90971-7>. **PLEASE DOWNLOAD**. The (e,e′d) measurement of the α–d relative momentum distribution. LiPolGen validates `Wave::radial` / `data/vmc/li6_alpha_d/li6.ad` against ANL VMC only (E-13): this is theory-vs-data for the tagged channel's own observable.
11. **Mitchell** *et al.*, *Mechanism of the ⁶Li(e,e′α) reaction*, Phys. Rev. C **44** (1991) 2002 — <https://doi.org/10.1103/PhysRevC.44.2002>. **PLEASE DOWNLOAD**. The α-knockout mirror of entry 10. Its abstract states the two-body and three-body Q² dependence and the recoil-momentum dependence of the two-body α–d breakup were measured, and that "the data indicate that the reaction mechanism is quasielastic in these kinematics" — i.e. the same relative-momentum distribution the T1 channel samples, fragments exchanged.
12. **Heimlich** *et al.*, *High-energy electron scattering from ⁶Li and ¹²C*, Nucl. Phys. **A231** (1974) 509 — <https://doi.org/10.1016/0375-9474(74)90514-4>. **PLEASE DOWNLOAD**. The primary source of every one of the 133 points in the free QES-archive ⁶Li dataset (three settings, source key `Heimlich:1973`, re-fetched and re-counted in this pass). Cite it, not only the archive, when `RC_QE_KF_GEV` becomes a fit residual.
13. **Moniz, Sick, Whitney, Ficenec, Kephart, Trower**, *Nuclear Fermi Momenta from Quasielastic Electron Scattering*, Phys. Rev. Lett. **26** (1971) 445 — <https://doi.org/10.1103/PhysRevLett.26.445>. **PLEASE DOWNLOAD** (INSPIRE recid 67648, no fulltext). The paper `RC_QE_KF_GEV = 0.169` (⁶Li) is taken from; `rc.hpp` calls the QE Fermi momentum the tail's dominant knob. Read the ⁶Li row off its own table before the value is next quoted with a fit-implying digit count.
14. **Amaudruz** *et al.* (NMC), *The ratio F₂ⁿ/F₂ᵖ in deep inelastic muon scattering*, Nucl. Phys. B **371** (1992) 3 — <https://doi.org/10.1016/0550-3213(92)90227-3>. **PLEASE DOWNLOAD** (INSPIRE recid 321412, CERN-PPE-91-167; CDS unreachable from this environment). Reference [17] of the HERMES b₁ paper: rebuilding HERMES's F₂ᵈ = F₂ᵖ(1 + F₂ⁿ/F₂ᵖ)/2, and therefore their F₁ᵈ denominator, exactly, requires it.
15. **Ent, Filippone, Makins, Milner, O'Neill, Wasson**, *Radiative corrections for (e,e′p) reactions at GeV energies*, Phys. Rev. C **64** (2001) 054610 — <https://doi.org/10.1103/PhysRevC.64.054610>. **PLEASE DOWNLOAD**. The coincidence-RC practice against which the tagged channels' `rc_tail ≡ 1` — the one place the generator asserts a correction is exactly 1.0 — should be priced, alongside the held Afanasev *et al.* calculation.
16. **Bardin, Shumeiko**, *An Exact Calculation of the Lowest Order Electromagnetic Correction to the Elastic Scattering*, Nucl. Phys. B **127** (1977) 242 — <https://doi.org/10.1016/0550-3213(77)90213-9>. **PLEASE DOWNLOAD** for the **English** text. The covariant infrared separation POLRAD Eq. (18) and `polrad_tpeak_quadrature` are built on. The free JINR-P2-10113 preprint *is* on disk (`refs/JINR-P2-10113.pdf`) but is Russian-language.

### 1b. Bounds something the tree ships as a choice

17. **Khan, Hoodbhoy**, *Convenient parametrization for deep-inelastic structure functions of the deuteron*, Phys. Rev. C **44** (1991) 1219 — <https://doi.org/10.1103/PhysRevC.44.1219>. **PLEASE DOWNLOAD** (INSPIRE recid 324952, no arXiv). A convolution model with relativistic and binding corrections giving b₁ as "a simple parametrization … in terms of a few deuteron wave-function parameters and the free nucleon structure functions" — exactly the shape one would refit for an α–d system.
18. **Khan, Hoodbhoy**, *Shadowing of deuteron spin structure functions*, Phys. Lett. B **298** (1993) 181 — <https://doi.org/10.1016/0370-2693(93)91727-5>. **PLEASE DOWNLOAD**. The same authors' small-x sequel: b₁ from double scattering in the light-cone parton model.
19. **de Vries, de Jager, de Vries**, *Nuclear charge-density-distribution parameters from elastic electron scattering*, At. Data Nucl. Data Tables **36** (1987) 495 — <https://doi.org/10.1016/0092-640X(87)90013-1>. **PLEASE DOWNLOAD**. The standard charge-density compilation (⁶Li MI30 rows; the ⁷Li HO row). Its *data* are free and machine-readable in the UVa archive (§8), but the article text is not.
20. **Li, Sick, Whitney, Yearian**, *High-energy electron scattering from ⁶Li*, Nucl. Phys. **A162** (1971) 583 — <https://doi.org/10.1016/0375-9474(71)90257-0>. **PLEASE DOWNLOAD** (already a corpus entry with `file: null`; the DOI above is new from this pass, and no free copy exists — OpenAlex closed, no arXiv, DTIC not served here). 500 MeV, q² to 13 fm⁻², the ⁶Li diffraction minimum at q² = 8 fm⁻² (q ≈ 2.83 fm⁻¹). *Serves:* the T11 `q₀` window and `T_DIP` in the coherent channel.
21. **Bergstrom**, *⁶Li electromagnetic form factors and phenomenological cluster models*, Nucl. Phys. **A327** (1979) 458 — <https://doi.org/10.1016/0375-9474(79)90269-0>. **PLEASE DOWNLOAD**. A phenomenological **α–d cluster-model** fit to ⁶Li's electromagnetic form factors — the published attempt to make LiPolGen's two ⁶Li descriptions (α+d clustering and a harmonic-oscillator form factor) the same object. The natural external check on `ClusterConfigSampler`'s charge distribution; absent from the tree entirely.
22. **Bergstrom, Kowalski, Neuhausen**, *Elastic magnetic form factor of Li⁶*, Phys. Rev. C **25** (1982) 1156 — <https://doi.org/10.1103/PhysRevC.25.1156>. **PLEASE DOWNLOAD**. Supersedes entry 8's ⁶Li half at higher momentum transfer; the better single purchase if only one magnetization paper is bought.
23. **van Niftrik, Lapikás, de Vries, Box**, *Magnetization distribution of the ⁷Li nucleus as obtained from electron scattering through 180°. The electric quadrupole moment of ⁷Li*, Nucl. Phys. **A174** (1971) 173 — <https://doi.org/10.1016/0375-9474(71)91011-6>. **PLEASE DOWNLOAD**. The *electron-scattering* determination of Q(⁷Li) — an independent route to `LI7_QUADRUPOLE_FM2`, which today rests on TUNL alone.
24. **Lichtenstadt, Alster, Moinester, Dubach, Hicks, Peterson, Kowalski**, *High momentum transfer longitudinal and transverse form factors of the ⁷Li ground-state doublet*, Phys. Lett. B **219** (1989) 394 — <https://doi.org/10.1016/0370-2693(89)91083-6>. **PLEASE DOWNLOAD**. ⁷Li F_L and F_T above Suelzle's q range.
25. **Ent** *et al.*, *Reaction ⁶Li(e,e′d)⁴He and the α–d Momentum Distribution in the Ground State of ⁶Li*, Phys. Rev. Lett. **57** (1986) 2367 — <https://doi.org/10.1103/PhysRevLett.57.2367>. **PLEASE DOWNLOAD**. The letter version of entry 10; its title is literally the observable the tagged channel samples (verified through Crossref).
26. **Genin, Julien, Rambaut, Samour, Palmeri, Vinciguerra**, *The ⁶Li(e,e′α) and ⁶Li(e,e′d) reactions at 520 MeV*, Phys. Lett. B **52** (1974) 46 — <https://doi.org/10.1016/0370-2693(74)90714-X>. **PLEASE DOWNLOAD**. First-generation version of entries 10/11. **Note a citation fix:** `benchmarking/03_data_nuclear.md` lists this as Phys. Lett. B **51** (1974); the PII resolves to volume **52**, page 46. PLB **51** (1974) 217 is a different, theoretical paper (Toyama & Sakamoto).
27. **Albrecht** *et al.*, *Large-angle quasi-free scattering in ⁶Li(p,pd)⁴He at 670 MeV*, Nucl. Phys. **A338** (1980) 477 — <https://doi.org/10.1016/0375-9474(80)90045-7>. **FREE PREPRINT LOCATED 2026-09-16 — `free, not fetched`.** INSPIRE record 143652 serves the Dubna preprint **JINR-E1-12727** (1979) at <https://inspirehep.net/files/8e09e655183d322b5ea66f25679dc428> — checked in the reconciliation pass (§10.4): HTTP 200, `application/pdf`, 972,423 bytes, PDF 1.6, 10 pp., page 1 reads *LARGE-ANGLE QUASI-FREE SCATTERING IN ⁶Li(p,pd)⁴He AT 670 MeV*, D. Albrecht, M. Csatlós, J. Erő, Z. Fodor *et al.* That pass was not permitted to add files to the corpus, so it is **not on disk**: fetch it (preprint text, not the NPA pagination) before buying the journal article. The first draft listed this item as purchase-only. A hadronic probe of the same α–d distribution, and where the factor-two width disagreement (FWHM ≈ 70 vs 120 MeV/c) lives — the honest external band on the distribution's width.
28. **Dollhopf, Perdrisat, Kitching, Olsen**, *The ⁶Li(α,2α)²H reaction at 700 MeV and the α–d momentum distribution*, Phys. Lett. B **58** (1975) 425 — <https://doi.org/10.1016/0370-2693(75)90579-1>. **PLEASE DOWNLOAD**. The third, α-induced probe of the same distribution.
29. **Whitney, Sick, Ficenec, Kephart, Trower**, *Quasielastic electron scattering*, Phys. Rev. C **9** (1974) 2230 — <https://doi.org/10.1103/PhysRevC.9.2230>. **PLEASE DOWNLOAD**. The full-length companion to entry 13, containing the large-angle ⁶Li data the QES archive's ⁶Li page does not carry.
30. **Rodning, Knutson**, *Asymptotic D-state to S-state ratio of the deuteron*, Phys. Rev. C **41** (1990) 898 — <https://doi.org/10.1103/PhysRevC.41.898>. **PLEASE DOWNLOAD** (INSPIRE recid 312974, no arXiv). The **measured** η_d = 0.0256(4), from sub-Coulomb d + ⁴He tensor analysing powers — same technique and group as the George–Knutson η(⁶Li → α+d) the tree gates on. `test_cluster.cpp` already carries 0.0256(4), but via Machleidt's Table XV, i.e. as a model output rather than a measurement.
31. **George, Knutson**, *Determination of the ⁶Li → α + d asymptotic D- to S-state ratio by a restricted phase shift analysis*, Phys. Rev. C **59** (1999) 598 — <https://doi.org/10.1103/PhysRevC.59.598>. **PLEASE DOWNLOAD** (INSPIRE 522531, no arXiv). `LI6_ETA_DS_GK` and gate **C-1**. The title, verified through Crossref, confirms `06_critic.md`'s objection: it is a *phase-shift analysis*, not a direct tensor-analysing-power measurement, which is why `BENCHMARK_PLAN.md` §4 item 8 asks for it to be retitled in the tree. §10.4 of `01_li-nuclear-data.md` records that **no newer determination of η(⁶Li → α+d) has been published since**, searched three ways.
32. **Benvenuti** *et al.* (BCDMS), *A high statistics measurement of the deuteron structure functions F₂(x,Q²) and R from deep inelastic muon scattering at high Q²*, Phys. Lett. B **237** (1990) 592 — <https://doi.org/10.1016/0370-2693(90)91231-Y>. **PLEASE DOWNLOAD** (INSPIRE 285497, CERN-EP-89-170). The high-x, high-Q² F₂ᵈ anchor, in the region where `ToyF2`'s by-eye shape is least constrained.
33. **Kuraev, Fadin**, *On Radiative Corrections to e⁺e⁻ Single Photon Annihilation at High Energy*, Sov. J. Nucl. Phys. **41** (1985) 466 [Yad. Fiz. **41** (1985) 733] — no DOI, no e-print; INSPIRE record <https://inspirehep.net/literature/217313>. **PLEASE DOWNLOAD**. The electron structure-function (equivalent-radiator) method and the **exponentiated** soft-photon form. `ll_radiator` implements the single-z leading log and its own comment concedes the soft 1/(1−z) is uncancelled; this is the citation the code does not carry. **Convention warning carried forward from `01_radiative.md` §3.6:** no copy was read in this pass, sources differ over α/2π vs α/π and over the −1 in β, and `rc.hpp` uses α/π — whoever writes the citation must pin the convention against a copy they have actually read.
34. **de Forest Jr., Walecka**, *Electron scattering and nuclear structure*, Adv. Phys. **15** (1966) 1 — <https://doi.org/10.1080/00018736600101254>. **PLEASE DOWNLOAD** (INSPIRE recid 50061). The source of `pauli_suppression`'s S(q) = (3/4)u − u³/16.
35. **Kwiatkowski, Spiesberger, Möhring**, *HERACLES: an event generator for ep interactions at HERA energies including radiative processes*, Comput. Phys. Commun. **69** (1992) 155 (DESY-90-041) — <https://doi.org/10.1016/0010-4655(92)90136-M>. **PLEASE DOWNLOAD** (INSPIRE recid 296225, no arXiv, no fulltext). Single-photon emission from the lepton line plus the complete one-loop weak corrections at **collider** kinematics — the actual independent RC engine inside DJANGOH, i.e. the second opinion gates `D-3`/`C-8` want.
36. **Charchula, Schuler, Spiesberger**, *Combined QED and QCD radiative effects in deep inelastic lepton–proton scattering: the Monte Carlo generator DJANGO6*, Comput. Phys. Commun. **81** (1994) 381 (CERN-TH-7133-94) — <https://doi.org/10.1016/0010-4655(94)90086-8>. **PLEASE DOWNLOAD** (INSPIRE recid 372027; CDS served an anti-bot challenge here). Defines the methodology of the published Rad/noRad table that `BENCHMARK_PLAN.md` §4 item 7 rates a top-ten, zero-dependency item.
37. **Mankiewicz, Schäfer, Veltri**, *PEPSI: a Monte Carlo generator for polarized leptoproduction*, Comput. Phys. Commun. **71** (1992) 305 — <https://doi.org/10.1016/0010-4655(92)90016-R>. **PLEASE DOWNLOAD** (INSPIRE recid 321673, HD-THEP-91-47 / MPIH-V38-1991). This pass **lifts a documented UNVERIFIED**: `benchmarking/01` §2.5 had the CPC citation from a README only; it is now confirmed against INSPIRE, and the abstract confirms PEPSI is a modification of LEPTO 4.3 generating the hard γ*–parton scattering with the polarization-dependent QCD cross section. It is the only code with a live tensor branch (`pnrun = ±2`).

### 1c. Completeness, context and provenance

38. **Mo, Tsai**, *Radiative Corrections to Elastic and Inelastic ep and μp Scattering*, Rev. Mod. Phys. **41** (1969) 205 — <https://doi.org/10.1103/RevModPhys.41.205>. **PLEASE DOWNLOAD** *only for the published equation numbering* — the January 1968 preprint SLAC-PUB-380, which carries Section III's exact elastic radiative tail, Appendix B's exact bremsstrahlung expression and the paper's own exact-vs-peaking comparison, **is on disk** (`refs/SLAC-PUB-0380.pdf`). Cite a gate as "SLAC-PUB-380 Eq. (B.5)" unless someone checks the RMP text.
39. **Tsai**, *Radiative corrections to electron–proton scattering*, Phys. Rev. **122** (1961) 1898 — <https://doi.org/10.1103/PhysRev.122.1898>. **PLEASE DOWNLOAD**, historical only (title and journal reference confirmed in this pass against INSPIRE record 46623, which carries no e-print and no fulltext): the root of the elastic-tail formula, superseded by SLAC-PUB-380/848 which are both held.
40. **Angeli, Marinova**, *Table of experimental nuclear ground state charge radii: An update*, At. Data Nucl. Data Tables **99** (2013) 69 — <https://doi.org/10.1016/j.adt.2011.12.006>. **PLEASE DOWNLOAD**. The source of the r_ch(⁶Li) = 2.589(39) fm behind `LI6_R2_POINT_FM2` = 6.0788 fm².
41. **Sick**, *Model-independent nuclear charge densities from elastic electron scattering*, Nucl. Phys. **A218** (1974) 509 — <https://doi.org/10.1016/0375-9474(74)90039-6>. **PLEASE DOWNLOAD**. The method paper behind the Fourier–Bessel rows of the UVa archive; it defines the incompleteness error that any gate built on those coefficients must carry.
42. **Stone**, *Table of nuclear electric quadrupole moments*, At. Data Nucl. Data Tables **111–112** (2016) 1 — <https://doi.org/10.1016/j.adt.2015.12.002>. **PLEASE DOWNLOAD**. The journal twin of the free IAEA report already recorded in the corpus; listed for citation completeness behind `LI6_QUADRUPOLE_FM2` / `LI7_QUADRUPOLE_FM2`.
43. **Gomez** *et al.* (E139), *Measurement of the A-dependence of deep-inelastic electron scattering*, Phys. Rev. D **49** (1994) 4348 — <https://doi.org/10.1103/PhysRevD.49.4348>. **RESOLVED 2026-09-16 — no purchase needed.** This entry was wrong in the first pass. INSPIRE record 359103 serves the **SLAC-PUB-5813** preprint free (<https://inspirehep.net/files/5634062e31e8530af3a67a9b238f5bd0>; checked here: HTTP 200, `application/pdf`, 5,687,755 bytes, 86 pp., page 1 reads *SLAC-PUB-5813 / August 1993 / MEASUREMENT OF THE A-DEPENDENCE OF DEEP-INELASTIC ELECTRON SCATTERING*), and it is now on disk as `refs/SLAC-PUB-5813.pdf`. Cited for a **negative**: the canonical A-dependence programme has no lithium target, which is why the ⁶Li EMC baseline is an nPDF fit rather than a measurement — and that negative is readable in the preprint. Buy the RMP-style published text only if exact Phys. Rev. D pagination is needed. Entry in §4c.
44. **Frosch, McCarthy, Rand, Yearian**, *Structure of the He⁴ Nucleus from Elastic Electron Scattering*, Phys. Rev. **160** (1967) 874 — <https://doi.org/10.1103/PhysRev.160.874>. **PLEASE DOWNLOAD**. One of the two datasets under the free ⁴He sum-of-Gaussians (§8) that `AlphaCoreSource` should be gated on.
45. **Ottermann, Köbschall, Maurer, Röhrich, Schmitt, Walther**, *Elastic electron scattering from ³He and ⁴He*, Nucl. Phys. **A436** (1985) 688 — <https://doi.org/10.1016/0375-9474(85)90554-8>. **PLEASE DOWNLOAD**. The low-q half of the same ⁴He fit.
46. **Goertz, Meyer, Reicherz**, *Polarized H, D and ³He targets for particle physics experiments*, Prog. Part. Nucl. Phys. **49** (2002) 403 — <https://doi.org/10.1016/S0146-6410(02)00159-X>. **PLEASE DOWNLOAD**. The standard review containing the equal-spin-temperature closed form behind `spin_temperature_pzz`.
47. **Goertz** *et al.*, *Investigations in high temperature irradiated ⁶,⁷LiH and ⁶LiD, its dynamic nuclear polarization and radiation resistance*, Nucl. Instrum. Meth. A **356** (1995) 20 — <https://doi.org/10.1016/0168-9002(94)01437-X>. **PLEASE DOWNLOAD**. The ⁶LiD material paper itself.
48. **Goertz**, *The dynamic nuclear polarization process*, Nucl. Instrum. Meth. A **526** (2004) 28 — <https://doi.org/10.1016/j.nima.2004.03.147>. **PLEASE DOWNLOAD**. Where the EST relation is derived; `BENCHMARK_PLAN.md` §4 item 9. (Whose honest half remains that **P_zz of ⁶Li has never been measured by anyone** — nothing found in this pass changes that.)
49. **Good, Walker**, *Diffraction dissociation of beam particles*, Phys. Rev. **120** (1960) 1857 — <https://doi.org/10.1103/PhysRev.120.1857>. **PLEASE DOWNLOAD** (INSPIRE recid 9375, no arXiv, no INSPIRE fulltext). The eigenstate decomposition — coherent amplitude = the average over absorption eigenstates — that `write_snd_configs` / `ClusterConfigSampler` exist to perform.
50. **Miettinen, Pumplin**, *Diffraction scattering and the parton structure of hadrons*, Phys. Rev. D **18** (1978) 1696 — <https://doi.org/10.1103/PhysRevD.18.1696>. **PLEASE DOWNLOAD**. *Reconciliation note 2026-09-16:* the first draft said INSPIRE serves a free scan at <https://inspirehep.net/files/8c3342725e4cd03851b56f8078d4e1fa> (HTTP 200, `application/pdf`, 1.37 MB) — the headers are as stated, but the 1,371,932 bytes served to this environment (fetched twice, two user agents) carry no `%PDF` header and do not parse (`pdfinfo`: not a PDF), so that link is **not** a usable copy from here. INSPIRE record 129259 names the original as the Fermilab preprint **FERMILAB-PUB-78-021-T**, <https://lss.fnal.gov/archive/1978/pub/Pub-78-021-T.pdf>, which answered an automated fetch with HTTP 403 (bot filter); a browser will probably fetch it, and it is the free copy to try before buying the journal text. Coherent = ⟨A⟩², total diffractive = ⟨A²⟩, so incoherent = the variance — the sentence `ClusterConfigSet` implements.
51. **Rachek** *et al.*, *Measurement of Tensor Analyzing Power T₂₀ in Coherent π⁰ Photoproduction on Deuteron*, Few Body Syst. **58** (2017) 29 — <https://doi.org/10.1007/s00601-016-1191-0>. **PLEASE DOWNLOAD** (INSPIRE recid 1508664, no arXiv). Precedent, not input: the only *measured* tensor observable of a coherent production reaction on a nucleus anywhere. At 0.2–0.5 GeV it is resonance-region nuclear structure, not diffraction, and it constrains no number in `coherent.hpp` — but the O5 write-up should be able to say it exists.
52. **Azhgirey** *et al.*, *New data on tensor analyzing power A_yy of the relativistic deuteron breakup as an additional test of deuteron structure at small distances*, Phys. Lett. B **595** (2004) 151 — <https://doi.org/10.1016/j.physletb.2004.05.057>. **PLEASE DOWNLOAD** for a readable PDF. Part of the Dubna rank-2 proxy set for `U-2`. **Partly resolved in this pass:** INSPIRE record 1464537 serves the **full text as Elsevier XML**, not as a PDF, at <https://inspirehep.net/files/bd01c549fb5585bac195643e31241376> (checked here: HTTP 200, `text/plain`, 71,506 bytes, whose `<ce:title>`, `<prism:doi>` and author surnames match). It was **not** committed — an XML document does not belong in a PDF corpus — but the text is obtainable without the paywall if a number is urgently needed.
53. **Azhgirey** *et al.*, *Measurement of the tensor analyzing power T₂₀ in inclusive deuteron breakup at 9 GeV/c on hydrogen and carbon*, Phys. Lett. B **387** (1996) 37 — <https://doi.org/10.1016/0370-2693(96)01007-6>. **PLEASE DOWNLOAD**. The spherical-convention sibling of entry 52.
54. **Ladygin** *et al.*, *Measurement of the tensor-analyzing power A_yy in deuteron breakup at 4.5 GeV/c and 80 mr*, Few Body Syst. **32** (2002) 127 — <https://doi.org/10.1007/s00601-002-0115-3>. **PLEASE DOWNLOAD**. Same proxy set. *For all three:* the analysing-power convention is **not** the target-polarization convention (Madison vs P_zz differ by ≈ √2 and by sign usage) — resolve it before comparing any number, per `06_critic.md` trap T-8.

### 1d. Open access, but blocked here — a browser will fetch these

55. **Connelly** *et al.*, *Trinucleon cluster knockout from ⁶Li*, Phys. Rev. C **57** (1998) 1569 — OA copies at <http://link.aps.org/pdf/10.1103/PhysRevC.57.1569> (bronze) and <https://research.vu.nl/files/1797325/144572.pdf> (green). **PLEASE DOWNLOAD** — both refused an automated fetch here (HTTP 403 bot filter). It is the *caution* on every α–d knockout number: a significant mirror-channel deviation at low momentum transfer, i.e. two-step contamination that any use of the cluster-knockout data as a wave-function measurement inherits.
56. **Pyykkö**, *Year-2008 nuclear quadrupole moments*, Mol. Phys. **106** (2008) 1965 — green OA at <https://hal.science/hal-00513184>; doi:<https://doi.org/10.1080/00268970802018367>. **PLEASE DOWNLOAD** (bot-blocked here). The compilation whose Q(⁶Li) = −0.0806 fm² gate **T11 explicitly tests against** — the tree asserts TUNL's −0.0818, not this.
57. **Borremans** *et al.*, *New measurement and reevaluation of the nuclear magnetic and quadrupole moments of ⁸Li and ⁹Li*, Phys. Rev. C **72** (2005) 044309 — bronze OA at <http://link.aps.org/pdf/10.1103/PhysRevC.72.044309>. **PLEASE DOWNLOAD** (bot-blocked here). The primary measurement Stone's ⁶Li quadrupole entry traces to (`2005Bo45`); it closes the provenance chain `LI6_QUADRUPOLE_FM2` → TUNL/Stone → measurement.
58. **Fricke, Bernhardt, Heilig, Schaller, Schellenberg, Shera, de Jager**, *Nuclear Ground State Charge Radii from Electromagnetic Interactions*, At. Data Nucl. Data Tables **60** (1995) 177 — <https://doi.org/10.1006/adnd.1995.1007>; reported green OA through a repository. **PLEASE DOWNLOAD**. The third charge-radius compilation and a useful independent cross-check on `LI6_R2_POINT_FM2`.

> **Open access and simply not fetched** (no action needed unless wanted): Tsai, *Pair production and bremsstrahlung of charged leptons*, Rev. Mod. Phys. **46** (1974) 815 = SLAC-PUB-1365, free at <https://inspirehep.net/files/be2ffc45ed59b1b28208cc386c1e16b3> (6.7 MB — external bremsstrahlung does not exist at a collider, so only the radiator half applies); DD4hep and GEANT4 (§7); the PDG *Review of Particle Physics* (§8); and the EIC CDR / ePIC PDR / ePIC pTDR, deliberately uncommitted on size grounds (§8).

---

## 2. T1 — the nucleon layer

Structure functions, R = σ_L/σ_T, nucleon form factors and the proton PDF fits
the free-nucleon backends read.

- **Abe** *et al.* (E143), *Measurements of R = σ_L/σ_T for 0.03 < x < 0.1 and Fit to World Data*, Phys. Lett. B **452** (1999) 194 (SLAC-PUB-7927). [arXiv:hep-ex/9808028](https://arxiv.org/abs/hep-ex/9808028) — `held`, `refs/hep-ex_9808028.pdf`.
  The R1998 world fit. It is LiPolGen's default R, transcribed as `r1998()` (`src/core/sf.cpp`, forms `kR1998A/B/C` and their spread). *Serves:* `f1_from_f2` everywhere, and the A_zz normalisation through F₁.
- **Whitlow, Rock, Bodek, Riordan, Dasu**, *A precise extraction of R = σ_L/σ_T from a global analysis of the SLAC deep inelastic e-p and e-d scattering cross sections*, Phys. Lett. B **250** (1990) 193; preprint SLAC-PUB-5284. Journal doi:[10.1016/0370-2693(90)91176-C](https://doi.org/10.1016/0370-2693(90)91176-C) (paywalled); free preprint at <https://inspirehep.net/files/d4b4fcdb5d834b34d41d3226090c37ea> — `held`, `refs/whitlow1990_R.pdf`.
  The preprint's own abstract gives R^p, R^d and R^d − R^p "over the entire SLAC kinematic range: 0.1 ≤ x ≤ 0.9 and 0.6 ≤ Q² ≤ 20.0 (GeV/c)²" and finds R^p = R^d. Two jobs: it is **the R the HERMES b₁ᵈ values were extracted with** (their ref. [18], read from the PDF), so it belongs in any comparison against HERMES Table II; and it is the actual reference behind `sf.hpp:43`'s unsourced "simplified R1990-like magnitude" claim for `r_sigma_lt`.
- **Dasu** *et al.* (SLAC E140), *Measurement of kinematic and nuclear dependence of R = σ_L/σ_T in deep inelastic electron scattering*, Phys. Rev. D **49** (1994) 5641; preprint SLAC-PUB-5814. doi:[10.1103/PhysRevD.49.5641](https://doi.org/10.1103/PhysRevD.49.5641) — `held`, `refs/SLAC-PUB-5814.pdf` (INSPIRE-hosted scan).
  The dominant input to the E143 R1998 world fit, so comparing `r1998` to E140 is the closest available test of the transcription against data rather than against itself (today `tests/test_sf.cpp` only checks internal consistency). It also carries R_A − R_D, the only handle anywhere on whether R on a nucleus differs from R on the deuteron — an assumption `NuclearF2` and the ⁶Li R default both make silently. Caveat: the scan's OCR layer is poor; read the tables from the page images.
- **Aaron** *et al.* (H1), *Measurement of the proton structure function F_L at low x*, Phys. Lett. B **665** (2008) 139. [arXiv:0805.2809](https://arxiv.org/abs/0805.2809) · doi:[10.1016/j.physletb.2008.05.070](https://doi.org/10.1016/j.physletb.2008.05.070) — `held`, `refs/0805.2809.pdf`.
  F_L over 12 < Q² < 90 GeV², 0.00024 < x < 0.0036 — with R = F_L/(F₂ − F_L), among the only measured R inside the EIC's low-x window, where `r1998` is frozen at `R1998_X_MIN` = 0.005 by `clip = true`. Caveat for any harness row: proton, not deuteron.
- **Andreev** *et al.* (H1), *Measurement of inclusive ep cross sections at high Q² at √s = 225 and 252 GeV and of the longitudinal proton structure function F_L at HERA*, Eur. Phys. J. C **74** (2014) 2814. [arXiv:1312.4821](https://arxiv.org/abs/1312.4821) · doi:[10.1140/epjc/s10052-014-2814-6](https://doi.org/10.1140/epjc/s10052-014-2814-6) — `held`, `refs/1312.4821.pdf`.
  The widest H1 F_L extraction, 1.5 ≤ Q² ≤ 800 GeV². Together with the entry above it prices the cost of `r1998`'s low-x extrapolation and the bias `r_sigma_lt` — the default, and by its own docstring "not a fit" — carries into every `f1_from_f2`.
- **Airapetian** *et al.* (HERMES), *Inclusive measurements of inelastic electron and positron scattering from unpolarized hydrogen and deuterium targets*, JHEP **05** (2011) 126. [arXiv:1103.5704](https://arxiv.org/abs/1103.5704) · doi:[10.1007/JHEP05(2011)126](https://doi.org/10.1007/JHEP05(2011)126) — `held`, `refs/1103.5704.pdf`.
  HERMES's own F₂ᵖ and F₂ᵈ, 0.006 ≤ x ≤ 0.9, 0.1 ≤ Q² ≤ 20 GeV², in the same apparatus and overlapping kinematics as the b₁ measurement — so the `B1_MILLER_TABLE_TO_PER_NUCLEON` inversion can be redone against a published table rather than a parametrisation. Also serves `ToyF2` / `LhapdfSF` / `MstwSF` directly.
- **Arneodo** *et al.* (NMC), *Measurement of the proton and deuteron structure functions F₂ᵖ and F₂ᵈ, and of the ratio σ_L/σ_T*, Nucl. Phys. B **483** (1997) 3. [arXiv:hep-ph/9610231](https://arxiv.org/abs/hep-ph/9610231) · doi:[10.1016/S0550-3213(96)00538-X](https://doi.org/10.1016/S0550-3213(96)00538-X) — `held`, `refs/hep-ph_9610231.pdf`.
  0.002 < x < 0.60, 0.5 < Q² < 75 GeV², **plus R for 0.002 < x < 0.12** — simultaneously the best fixed-target F₂ᵈ for the EIC window and the only R data at NMC's low x, where `r1998`'s support runs out. *Serves:* `ToyF2::f2p/f2n`, `LhapdfSF`, `MstwSF`.
- **Arneodo** *et al.* (NMC), *Accurate measurement of F₂ᵈ/F₂ᵖ and R^d − R^p*, Nucl. Phys. B **487** (1997) 3. [arXiv:hep-ex/9611022](https://arxiv.org/abs/hep-ex/9611022) · doi:[10.1016/S0550-3213(96)00673-6](https://doi.org/10.1016/S0550-3213(96)00673-6) — `held`, `refs/hep-ex_9611022.pdf`.
  F₂ᵈ/F₂ᵖ over 0.001 < x < 0.8, 0.1 < Q² < 145 GeV², with R^d − R^p compatible with zero over 0.002 < x < 0.4 — the empirical licence for using one R for both p and d, which `f1_from_f2` does implicitly everywhere.
- **Arneodo** *et al.* (NMC), *Measurement of the proton and the deuteron structure functions F₂ᵖ and F₂ᵈ*, Phys. Lett. B **364** (1995) 107. [arXiv:hep-ph/9509406](https://arxiv.org/abs/hep-ph/9509406) · doi:[10.1016/0370-2693(95)01318-9](https://doi.org/10.1016/0370-2693(95)01318-9) — `held`, `refs/hep-ph_9509406.pdf`.
  The closed-form NMC F₂ parametrisation (0.006 < x < 0.9) — the natural replacement for `ToyF2`'s by-eye anchors when LHAPDF is not built.
- **Abrams** *et al.* (JLab Hall A Tritium / MARATHON), *Measurement of the nucleon F₂ⁿ/F₂ᵖ structure function ratio by the Jefferson Lab MARATHON tritium/helium-3 deep inelastic scattering experiment*, Phys. Rev. Lett. **128** (2022) 132003. [arXiv:2104.05850](https://arxiv.org/abs/2104.05850) · doi:[10.1103/PhysRevLett.128.132003](https://doi.org/10.1103/PhysRevLett.128.132003) — `held`, `refs/2104.05850.pdf`.
  ³H/³He mirror symmetry at 0.19 < x < 0.83, which "essentially eliminates many theoretical uncertainties in the extraction of the ratio". The cleanest modern test of `f2n_over_f2p` at large x, and the only modern experimental contact with the A = 3 system `CiofiSimulaTriton` models. It measures a structure-function ratio, **not** a spectral function.
- **Abramowicz, Levy**, *The ALLM parameterization of σ_tot(γ*p) — an update*, DESY 97-251. [arXiv:hep-ph/9712415](https://arxiv.org/abs/hep-ph/9712415) — `held`, `refs/hep-ph_9712415.pdf`.
  The Regge-plus-hard-Pomeron form fitted to 1356 F₂ points with χ²/ndf = 0.97 over 3×10⁻⁶ < x < 0.85, 0 ≤ Q² < 5000 GeV². It is in the corpus for one reason: it is reference [16] of the HERMES b₁ paper, so rebuilding HERMES's F₁ᵈ denominator needs it. Nothing in LiPolGen uses ALLM today; it is a benchmark-harness input.
- **Airapetian** *et al.* (HERMES), *Precise determination of the spin structure function g₁ of the proton, deuteron and neutron*, Phys. Rev. D **75** (2007) 012007. [arXiv:hep-ex/0609039](https://arxiv.org/abs/hep-ex/0609039) · doi:[10.1103/PhysRevD.75.012007](https://doi.org/10.1103/PhysRevD.75.012007) — `held`, `refs/hep-ex_0609039.pdf`.
  g₁ᵖ, g₁ᵈ and an extracted g₁ⁿ over 0.0041 ≤ x ≤ 0.9, 0.18 ≤ Q² ≤ 20 GeV². It is also the **radiative** reference for the b₁ measurement: its section *"Unfolding of Radiative and Instrumental Smearing"* and Appendix A (headings read from the PDF) give in full the migration/smearing algorithm the b₁ PRL compresses into one sentence. An INSPIRE title search confirms there is **no** long companion paper to the b₁ PRL; this is the document that carries its method. *Serves:* `ToyG1::g1_nucleus`, `a1p`/`a1n`, and the definition of what "Born" means in any comparison against HERMES Table II.
- **Abe** *et al.* (E143), *Measurements of the proton and deuteron spin structure functions g₁ and g₂*, Phys. Rev. D **58** (1998) 112003. [arXiv:hep-ph/9802357](https://arxiv.org/abs/hep-ph/9802357) · doi:[10.1103/PhysRevD.58.112003](https://doi.org/10.1103/PhysRevD.58.112003) — `held`, `refs/hep-ph_9802357.pdf`.
  Carries g₂ᵈ at 29.1 GeV and states in its abstract that "the g₂ data are well-described by the Wandzura–Wilczek twist-2 contribution" — precisely what `g2_ww` implements, and today `g2_ww` is tested only against its own analytic power-law limit. This is a **different** E143 paper from the R1998 fit above.
- **Anthony** *et al.* (E155), *Measurement of the deuteron spin structure function g₁ᵈ(x) for 1 (GeV/c)² < Q² < 40 (GeV/c)²*, Phys. Lett. B **463** (1999) 339. [arXiv:hep-ex/9904002](https://arxiv.org/abs/hep-ex/9904002) · doi:[10.1016/S0370-2693(99)00940-5](https://doi.org/10.1016/S0370-2693(99)00940-5) — `held`, `refs/hep-ex_9904002.pdf` (already a corpus entry).
  The widest-Q² g₁ᵈ set; in the tree it is cited for the ⁶LiD target dilution/polarization bookkeeping as much as for g₁.
- **Ye, Arrington, Hill, Lee**, *Proton and neutron electromagnetic form factors and uncertainties*, Phys. Lett. B **777** (2018) 8. [arXiv:1707.09063](https://arxiv.org/abs/1707.09063) · doi:[10.1016/j.physletb.2017.11.023](https://doi.org/10.1016/j.physletb.2017.11.023) — `held`, `refs/1707.09063.pdf`.
  A world fit of the nucleon electromagnetic form factors **with an error band** — a drop-in replacement for `nucleon_ff`'s dipole + Galster (`rc.hpp`). Low priority only because the `qe_suppression` ±100 % band dominates the quasi-elastic tail anyway.
- **Hou** *et al.* (CTEQ-TEA), *New CTEQ global analysis of quantum chromodynamics with high-precision data from the LHC*, Phys. Rev. D **103** (2021) 014013. [arXiv:1912.10053](https://arxiv.org/abs/1912.10053) · doi:[10.1103/PhysRevD.103.014013](https://doi.org/10.1103/PhysRevD.103.014013) — `held`, `refs/1912.10053.pdf`.
  `LhapdfSF`'s **default `CT18NLO`**, and `CT18Anlo` is EPPS21's proton baseline — so one paper stands behind both the free-nucleon F₂ and the nuclear-ratio denominator. The CT18/CT18A distinction (the latter includes the ATLAS 7 TeV W/Z data) is why the grid name reads `EPPS21nlo_CT18Anlo_Li6`.
- **Nocera, Ball, Forte, Ridolfi, Rojo** (NNPDF), *A first unbiased global determination of polarized PDFs and their uncertainties*, Nucl. Phys. B **887** (2014) 276. [arXiv:1406.5539](https://arxiv.org/abs/1406.5539) · doi:[10.1016/j.nuclphysb.2014.08.008](https://doi.org/10.1016/j.nuclphysb.2014.08.008) — `held`, `refs/1406.5539.pdf`.
  `LhapdfG1`'s default `NNPDFpol11_100`, including the replica-based uncertainties and the Δb = 0 convention the tree's header note attributes to this fit. Note the arXiv record carries the DOI but no `journal_ref`.

---

## 3. T2 — the deuteron, b₁ and the tagged A = 2 analogue

### 3a. The measurement, and the experimental programme

- **Airapetian** *et al.* (HERMES), *First Measurement of the Tensor Structure Function b₁ of the Deuteron*, Phys. Rev. Lett. **95** (2005) 242001. [arXiv:hep-ex/0506018](https://arxiv.org/abs/hep-ex/0506018) · doi:[10.1103/PhysRevLett.95.242001](https://doi.org/10.1103/PhysRevLett.95.242001) — `held`, `refs/hep-ex_0506018.pdf`. **Corpus gap G2, closed.**
  **The only tensor-DIS datum in the world.** 27 GeV positrons on a tensor-polarized deuterium gas target with negligible residual vector polarization; Table II is six rows at ⟨x⟩ = 0.012 … 0.452, ⟨Q²⟩ = 0.51 … 4.69 GeV², read digit-for-digit from the PDF in this pass, plus two published integrals (∫₀.₀₀₂^0.85 b₁ dx = (1.05 ± 0.34 ± 0.35)×10⁻², and ∫₀.₀₂^0.85 b₁ dx = (0.35 ± 0.10 ± 0.18)×10⁻² over the Q² > 1 GeV² range). Three further things were read out of it and matter: (i) the sentence under Eq. (5), "F₂ᵈ is calculated as F₂ᵈ = F₂ᵖ(1 + F₂ⁿ/F₂ᵖ)/2", settles `B1_MILLER_TABLE_TO_PER_NUCLEON` — the published b₁ᵈ **is** per nucleon — against a quoted line rather than an inversion; (ii) the A_zz definition in their Eq. (2) is what fixes the sign convention of `TENSOR_LL_SIGN`, not a generic one; (iii) the RC paragraph (quoted in full in `01_radiative.md` §7.1) says the radiative background reaches almost 50 % of the statistics in the lowest-x bin and that the residual RC systematic on A_zz is ≈ 2 × 10⁻³ at low x — which is where `RC_DELTA_LOW_X_OPTIMISTIC` = 0.19 comes from. *Serves:* `MillerB1`, `CdksB1`, `azz()`, `close_kumano_integral`, `RC_DELTA_LOW_X_OPTIMISTIC`, `qe_tensor_scale`'s default, and `BENCHMARK_PLAN.md` §4 item 3.
- **Poudel, Bacchetta, Chen, Santiesteban**, *Experimental Study of Tensor Structure Function of Deuteron*, Eur. Phys. J. A **61** (2025) 81. [arXiv:2506.04506](https://arxiv.org/abs/2506.04506) · doi:[10.1140/epja/s10050-025-01558-w](https://doi.org/10.1140/epja/s10050-025-01558-w) — `held`, `refs/2506.04506.pdf`.
  A review of the experimental tensor programme with no data of its own. Two specific uses: it is the citable peer-reviewed 2025 sentence for "one measurement exists, the second has not run", and `RC_X_HIGH` = 0.16 comes from its p. 8 (0.16 < x < 0.49, 0.8 < Q² < 5.0 GeV²). **It does *not* contain the 1.5 %**: `rc.hpp` records that `RC_DELTA_HIGH_X` = 0.015 appears only in the unpublished E12-13-011 proposal and that this paper has no radiative-correction discussion at all. Do not conflate the two provenances.
- **Slifer** *et al.*, *The Deuteron Tensor Structure Function b₁*, JLab proposal PR12-13-011 to PAC 40 (2013) — <https://www.jlab.org/exp_prog/proposals/13/PR12-13-011.pdf> — `held`, `refs/JLab_PR12-13-011.pdf`.
  The approved b₁ experiment: 0.16 < x < 0.49, 0.8 < Q² < 5.0 GeV², 30 days at 11 GeV, P_zz = 20 % (read from the PDF). The corpus entry `jlab-pr12-13-011-b1-and-2023-jeopardy` records `file: null`; the file now sits beside it.
- **Long, Slifer, Solvignon** *et al.*, *Measurements of the Quasi-Elastic and Elastic Deuteron Tensor Asymmetries*, JLab proposal PR12-15-005 to PAC 43 (2015) — <https://www.jlab.org/exp_prog/proposals/15/PR12-15-005.pdf> — `held`, `refs/JLab_PR12-15-005.pdf`.
  The A_zz proposal: quasi-elastic A_zz at x > 1 plus elastic T₂₀ over 0.2 < Q² < 1.8 GeV² (read from the PDF). Relevant to the quasi-elastic tensor tail the generator sets to zero.
- **Slifer**, *The deuteron polarized tensor structure function b₁*, J. Phys. Conf. Ser. **543** (2014) 012003. doi:[10.1088/1742-6596/543/1/012003](https://doi.org/10.1088/1742-6596/543/1/012003) — `held`, `refs/JPCS_543_012003.pdf` (IOP CC-BY).
  The four-page public summary of PR12-13-011 — the *citable* version of the proposal.
- **Long**, *Potential for a tensor asymmetry A_zz measurement in the x > 1 region at JLab*, J. Phys. Conf. Ser. **543** (2014) 012010. doi:[10.1088/1742-6596/543/1/012010](https://doi.org/10.1088/1742-6596/543/1/012010) — `held`, `refs/JPCS_543_012010.pdf`.
  The citable version of the A_zz case; it names Frankfurt–Strikman 1988 as the first quasi-elastic A_zz calculation.
- **Kalantarians**, *Tensor polarized deuteron at an electron-ion collider*, J. Phys. Conf. Ser. **543** (2014) 012008. doi:[10.1088/1742-6596/543/1/012008](https://doi.org/10.1088/1742-6596/543/1/012008) — `held`, `refs/JPCS_543_012008.pdf`.
  The earliest written case for tensor polarization *at an EIC* — the direct ancestor of this project's premise.
- **Slifer, Long**, *Novel physics with tensor polarized targets*, PoS **PSTP2013** (2013) 008. [arXiv:1311.4835](https://arxiv.org/abs/1311.4835) · doi:[10.22323/1.182.0008](https://doi.org/10.22323/1.182.0008) — `held`, `refs/1311.4835.pdf`.
  Programme-level motivation plus the target-technology context for ⁶LiD / ND₃.
- **Poudel, Bacchetta, Chen, Keller, Fernando, Long, Ruth, Santiesteban, Slifer**, *Spin-1 transverse-momentum-dependent tensor structure functions in CLAS12* (2025). [arXiv:2502.20044](https://arxiv.org/abs/2502.20044) — `held`, `refs/2502.20044.pdf`.
  The *next* b₁ datum: an RG-C ND₃ analysis extracting b₁ inclusively and F_U(LL),T and F^{cos2φ}_U(LL) in SIDIS.
- **Maxwell, Crabb, Day, Detmold, Jaffe, Jones, Keith, Keller, Meekins, Milner** *et al.*, *Search for exotic gluonic states in the nucleus*, Letter of Intent to JLab PAC 44 (dated 6 June 2016, posted 2018). [arXiv:1803.11206](https://arxiv.org/abs/1803.11206) — `held`, `refs/1803.11206.pdf`.
  **The only experimental programme anywhere aiming at Δ**: unpolarized electron beam on a transversely polarized spin-1 nuclear target, inclusive DIS below x = 0.3, via single-spin tensor asymmetries. Its author list carries both the lattice (Detmold) and the 1989 theory (Jaffe). It is also a free, quotable secondary source for *Nuclear gluonometry*'s definition of Δ, which the corpus lacks. *Serves:* `toy_delta_gluon`, `RcScope::TensorAll`, and the geometry of the coherent cos 2φ term.

### 3b. Tensor DIS on the deuteron — the observable side of the A = 2 gate

- **Cosyn, Roldan Tomei, Sosa, Zec**, *Polarization options in inclusive DIS off tensor polarized deuteron*, Eur. Phys. J. A **61** (2025) 83. [arXiv:2410.12764](https://arxiv.org/abs/2410.12764) — `held`, `refs/2410.12764v1.pdf` (corpus entry).
  The polarization-option analysis the tree's tensor observables are defined against; cited across `PHYSICS_CHANNELS.md`, `benchmarking/00,02,04` and `open_items/physics_literature.md`.
- **Cosyn, Weiss**, *Semi-inclusive deep-inelastic scattering on a polarized spin-1 target. I. Cross section and spin observables* and *II. Deuteron and spectator nucleon tagging*, Phys. Rev. C **114** (2026) 025201 and 025202. [arXiv:2603.23699](https://arxiv.org/abs/2603.23699), [arXiv:2603.23700](https://arxiv.org/abs/2603.23700) — `held`, `refs/2603.23699.pdf`, `refs/2603.23700.pdf` (corpus entries).
  The spin-1 SIDIS decomposition and the tagged formalism LiPolGen's `InclusiveKernel` and `TaggedModel` are built against; `[CW26]`/`[CWappD]` throughout the tree, and the source of gate **E-1**'s Table II.
- **Cosyn, Weiss**, *Polarized electron–deuteron deep-inelastic scattering with spectator nucleon tagging*, Phys. Rev. C **102** (2020) 065204. [arXiv:2006.03033](https://arxiv.org/abs/2006.03033) · doi:[10.1103/PhysRevC.102.065204](https://doi.org/10.1103/PhysRevC.102.065204) — `held`, `refs/2006.03033.pdf`.
  52 pages and 19 figures: the light-front polarized spectral function in closed form. `01_tensor-b1.md` rates it "the one cheap unexploited win" — it gates `TaggedModel`'s **VECTOR** sector, which E-1 does not touch at all.
- **Jentsch, Tu, Weiss**, *Deep-inelastic electron-deuteron scattering with spectator nucleon tagging at the future Electron-Ion Collider: Extracting free nucleon structure*, Phys. Rev. C **104** (2021) 065205. [arXiv:2108.08314](https://arxiv.org/abs/2108.08314) — `held`, `refs/2108.08314.pdf` (corpus entry).
  The tagged e+d design comparator for LiPolGen's tagged channels. **Citation correction from this pass:** "Jentsch–Tu–Zhang PRC 104 (2021) 065205" does not exist; the authors are Jentsch, Tu and **Weiss**.
- **Mäntysaari, Salazar, Schenke, Shen, Zhao**, *Spatial imaging of polarized deuterons at the Electron-Ion Collider*, Phys. Lett. B **858** (2024) 139053. [arXiv:2408.13213](https://arxiv.org/abs/2408.13213) · doi:[10.1016/j.physletb.2024.139053](https://doi.org/10.1016/j.physletb.2024.139053) — `held`, `refs/2408.13213v1.pdf` (corpus entry).
  Coherent J/ψ production off a **polarized deuteron** in the CGC with nucleon-level (u, w) structure; its Eq. (9) writes dσ/dΦ d|t| ∝ 1 + 2 Σ aₙ cos nΦ and it gives a₂, a₄ per deuteron spin state at x_P = 1.7×10⁻³. `coherent.hpp`'s `mantysaari_a2_deuteron()` is four digitized rows of its Fig. 4, and `cluster_config.hpp`'s `a2_from_quadrupole` is its Eq. (9) re-derived with a quadrupole, reproducing a₂(±1) to 8 % with no free parameter. **Re-checked in this pass through its forward citations on INSPIRE: it is still the only polarization-dependent coherent-diffraction calculation for any nucleus.**

### 3c. Vector-polarized deuteron, tagging and the deuteron control channel

- **Adolph** *et al.* (COMPASS), *Final COMPASS results on the deuteron spin-dependent structure function g₁ᵈ and the Bjorken sum rule*, Phys. Lett. B **769** (2017) 34. [arXiv:1612.00620](https://arxiv.org/abs/1612.00620) · doi:[10.1016/j.physletb.2017.03.018](https://doi.org/10.1016/j.physletb.2017.03.018) — `held`, `refs/1612.00620.pdf`.
  1 < Q² < 100 (GeV/c)², 0.004 < x < 0.7, W > 4 GeV/c² — **the best overlap with the EIC window of any g₁ᵈ set**, and the single cleanest external gate on `ToyG1::g1_nucleus` evaluated on `DEUTERON()`.
- **Alexakhin** *et al.* (COMPASS), *The deuteron spin-dependent structure function g₁ᵈ and its first moment*, Phys. Lett. B **647** (2007) 8. [arXiv:hep-ex/0609038](https://arxiv.org/abs/hep-ex/0609038) · doi:[10.1016/j.physletb.2006.12.076](https://doi.org/10.1016/j.physletb.2006.12.076) — `held`, `refs/hep-ex_0609038.pdf`.
  Superseded by the 2017 paper for g₁ᵈ; kept because it documents the ⁶LiD target bookkeeping COMPASS used.
- **Baillie** *et al.* (CLAS BONuS), *Measurement of the neutron F₂ structure function via spectator tagging with CLAS*, Phys. Rev. Lett. **108** (2012) 142001. [arXiv:1110.2770](https://arxiv.org/abs/1110.2770) · doi:[10.1103/PhysRevLett.108.142001](https://doi.org/10.1103/PhysRevLett.108.142001) — `held`, `refs/1110.2770.pdf`.
  DIS on deuterium with a backward spectator proton, p_s < 100 MeV/c, θ_pq > 100°, 0.65 < Q² < 4.52 GeV², 0.2 < x < 0.8 — LiPolGen's tagged topology with a nucleon in place of an α and the same physical argument. *Serves:* `TaggedChannel`, `TaggedModel`'s (α_s, p_T) grid, `TaggedSampler`, `boost_spectator`.
- **Tkachenko** *et al.* (CLAS BONuS), *Measurement of the structure function of the nearly free neutron using spectator tagging in inelastic ²H(e,e′p_s)X scattering with CLAS*, Phys. Rev. C **89** (2014) 045206. [arXiv:1402.2477](https://arxiv.org/abs/1402.2477) · doi:[10.1103/PhysRevC.89.045206](https://doi.org/10.1103/PhysRevC.89.045206) — `held`, `refs/1402.2477.pdf`.
  The same measurement with full systematics and the on-shell extrapolation. The frame mismatch must be stated on any harness row: BONuS tags backward in the target rest frame, the EIC tags a boosted α or d in B0/OMD/ZDC acceptance, and the acceptance-driven systematics do not transfer.
- **Egiyan** *et al.* (CLAS), *Experimental study of exclusive ²H(e,e′p)n reaction mechanisms at high Q²*, Phys. Rev. Lett. **98** (2007) 262502. [arXiv:nucl-ex/0701013](https://arxiv.org/abs/nucl-ex/0701013) · doi:[10.1103/PhysRevLett.98.262502](https://doi.org/10.1103/PhysRevLett.98.262502) — `held`, `refs/nucl-ex_0701013.pdf`.
  1.75 < Q² < 5.5 GeV², and its abstract gives the measured three-region structure the tagged channel lives in: p_n < 100 MeV/c the neutron is a spectator and PWIA describes the reaction; 100 < p_n < 750 MeV/c pn rescattering dominates; above that Δ production followed by NΔ → NN. That is a *measured* statement of where a plane-wave tagged calculation stops being valid — an external bound on where `GlauberFsiWeight` must carry its full band.
- **Boeglin** *et al.* (JLab Hall A), *Probing the high momentum component of the deuteron at high Q²*, Phys. Rev. Lett. **107** (2011) 262501. [arXiv:1106.0275](https://arxiv.org/abs/1106.0275) · doi:[10.1103/PhysRevLett.107.262501](https://doi.org/10.1103/PhysRevLett.107.262501) — `held`, `refs/1106.0275.pdf`.
  d(e,e′p) at Q² = 3.5 (GeV/c)² binned in recoil angle; in 35° ≤ θ_nq ≤ 45° FSI are predicted small, so the reduced cross sections give direct access to the high-momentum component out to 0.55 GeV/c. This turns `test_tagged.cpp` T25 from Hulthén-against-AV18 (two models) into model-against-data over exactly the k range the tagged spectator spectrum populates. Unpolarized and rest-frame, so it says nothing about `eff_pol_p`/`eff_pol_n` and the light-cone variables do not transfer directly.
- **Abbott** *et al.* (JLab t₂₀), *A precise measurement of the deuteron elastic structure function A(Q²)*, Phys. Rev. Lett. **82** (1999) 1379. [arXiv:nucl-ex/9810017](https://arxiv.org/abs/nucl-ex/9810017) · doi:[10.1103/PhysRevLett.82.1379](https://doi.org/10.1103/PhysRevLett.82.1379) — `held`, `refs/nucl-ex_9810017.pdf`.
  A(Q²) at six points between 0.66 and 1.80 (GeV/c)² — the magnitude of σ^el(d) in the RC elastic tail's control channel.
- **Abbott** *et al.* (JLab t₂₀), *Measurement of tensor polarization in elastic electron-deuteron scattering at large momentum transfer*, Phys. Rev. Lett. **84** (2000) 5053. [arXiv:nucl-ex/0001006](https://arxiv.org/abs/nucl-ex/0001006) · doi:[10.1103/PhysRevLett.84.5053](https://doi.org/10.1103/PhysRevLett.84.5053) — `held`, `refs/nucl-ex_0001006.pdf`.
  t₂₀, t₂₁, t₂₂ at six points over the same range — the measurements that separate G_C from G_Q.
- **Abbott** *et al.*, *Phenomenology of the deuteron electromagnetic form factors*, Eur. Phys. J. A **7** (2000) 421. [arXiv:nucl-ex/0002003](https://arxiv.org/abs/nucl-ex/0002003) · doi:[10.1007/PL00013629](https://doi.org/10.1007/PL00013629) — `held`, `refs/nucl-ex_0002003.pdf`.
  Three closed-form world-data parametrisations (I, II, III) of G_C, G_M and G_Q over 0–7 fm⁻¹, plus the extraction of the node of G_C. Two independent pulls: it is HERMES's ref. [11], the form factors behind the coherent and quasi-elastic radiative tails they subtracted; and it removes the acknowledged stand-in in LiPolGen's own deuteron control channel, where `tests/test_rc.cpp` states the deuteron form factor is a harmonic-oscillator **shape** normalised on μ and Q_d and that every T8(a)/T10 gate is therefore shape-blind. `benchmarking/02` put it bluntly: these parametrisations "exist … and are not used".
- **Zhou, Bouwhuis, Ferro-Luzzi, Passchier, Alarcon, Anghinolfi, Arenhövel** *et al.*, *Tensor Analyzing Powers for Quasi-Elastic Electron Scattering from Deuterium*, Phys. Rev. Lett. **82** (1999) 687. [arXiv:nucl-ex/9809002](https://arxiv.org/abs/nucl-ex/9809002) · doi:[10.1103/PhysRevLett.82.687](https://doi.org/10.1103/PhysRevLett.82.687) — `held`, `refs/nucl-ex_9809002.pdf`. **Corpus gap G7, closed.**
  The sole citation under `RcOptions::qe_tensor_scale = 0.0`. Its abstract: a first measurement of tensor analyzing powers in quasi-elastic ed scattering at ⟨q⟩ = 1.7 fm⁻¹, for missing momenta up to 150 MeV/c, well described by a calculation including FSI, meson-exchange and isobar currents and leading-order relativistic contributions. **Read carefully, it is a measurement of a non-zero effect, not a null result** — HERMES's sentence ("no net tensor effect by inclusive scattering on weakly-bound spin-1/2 objects [13]") is the claim, and this paper is its citation, not its proof. `phase_B_numbers.md` §B3 already withdrew the A = 2 → A = 6 transfer of the argument and records the quasi-elastic term as 22 % / 73 % / 99.9 % of the tail at x = 0.01 / 0.10 / 0.30.
- **Ladygin** *et al.*, *Tensor A_yy and vector A_y analyzing powers in H(d,d′)X and ¹²C(d,d′)X at 9 GeV/c*, Phys. At. Nucl. **69** (2006) 852. [arXiv:nucl-ex/0510050](https://arxiv.org/abs/nucl-ex/0510050) · doi:[10.1134/S1063778806050073](https://doi.org/10.1134/S1063778806050073) — `held`, `refs/nucl-ex_0510050.pdf`.
  The only free member of the Dubna/JINR set — measured rank-2 observables on a spin-1 nucleus, driven by S–D interference at internal momenta 0.2–1.0 GeV/c, i.e. the same f₀/f₂ structure as `A_zz^wf` and far beyond BLAST's ≤ 500 MeV/c. `06_critic.md` §4.3 identifies this literature as the closest measured proxy for `U-2`. It is not an electromagnetic probe, not a tagged α, and **not the target-polarization convention**.
- **Martin, Stirling, Thorne, Watt**, *Parton distributions for the LHC* (MSTW 2008), Eur. Phys. J. C **63** (2009) 189. [arXiv:0901.0002](https://arxiv.org/abs/0901.0002) · doi:[10.1140/epjc/s10052-009-1072-5](https://doi.org/10.1140/epjc/s10052-009-1072-5) — `held`, `refs/0901.0002.pdf`.
  A **gate-defining** reference, not a convenience: `b1_nuclear.hpp` states the CDKS-comparison gate passes for the MSTW2008 LO configuration specifically, because MSTW2008 LO is the unpolarized PDF CDKS computed their b₁ with. `OPEN_ITEMS_SOLUTIONS.md` records the gate's peak ratio as 0.843 with MSTW and 0.440 — outside the window — with the shipped `ToyF2`. The gate is therefore a statement about a *configuration*, and this is that configuration's reference.

---

## 4. T3 — the nuclear inputs (⁶Li, ⁷Li, ⁴He)

Everything the generator needs about the target itself: the elastic form
factors, the α–d cluster structure the tagged channels sample, the EMC
baseline, the quasi-elastic knobs of the radiative tail, and the ab-initio
chain under `data/vmc/`. Sourced from `01_li-nuclear-data.md`, whose master
table (A1–G3) carries the per-row verification.

### 4a. Elastic form factors and charge densities

- **Wiringa, Schiavilla**, *Microscopic calculation of ⁶Li elastic and transition form factors*, Phys. Rev. Lett. **81** (1998) 4317. [arXiv:nucl-th/9807037](https://arxiv.org/abs/nucl-th/9807037) — `held`, `refs/nucl-th_9807037.pdf`.
  VMC with AV18+UIX, the same Hamiltonian as the `data/vmc/` tables. It is the theory `rc.hpp` borrows its monopole **shape** from, and `01_li-nuclear-data.md` §2.4 records that its Q(⁶Li) = −0.23(9) fm² is a caution on the `C-2` quadrupole comparison, not an endorsement of it. The one ⁶Li form-factor reference in this group that is free.
- **Lapikás, Wesseling, Wiringa**, *Nuclear structure studies with the ⁷Li(e,e′p) reaction*, Phys. Rev. Lett. **82** (1999) 4404. [arXiv:nucl-th/9904008](https://arxiv.org/abs/nucl-th/9904008) — `held`, `refs/nucl-th_9904008.pdf`.
  ⁷Li(e,e′p)⁶He momentum distributions measured against VMC, spectroscopic factor 0.58 ± 0.05 vs 0.60 from VMC (`01_li-nuclear-data.md` §3d table). It is **gate E-13**: the only data-vs-VMC test anywhere of the `data/vmc/` family the tagged channel reads.
- **Hotta, Tamae, Miura, Miyase, Nakagawa, Suda, Sugawara, Tadokoro**, *Measurement of the ⁶Li(e,e′p) reaction cross sections at low momentum transfer*, Nucl. Phys. **A645** (1999) 492. [arXiv:nucl-ex/9810016](https://arxiv.org/abs/nucl-ex/9810016) — `held`, `refs/nucl-ex_9810016.pdf`.
  Proton knockout from ⁶Li at low q — the low-ω end of the region `RcTailModel`'s quasi-elastic term has to describe.

> **Also T3, PLEASE DOWNLOAD** (full entries in §1): **#7** Suelzle–Yearian–Crannell 1967 (the origin of both lithium charge form factors and of `HoSpin1FF`'s parameters), **#20** Li–Sick–Whitney–Yearian 1971 (the ⁶Li diffraction minimum behind `T_DIP`), **#21** Bergstrom 1979 (the α–d cluster-model fit), **#19** de Vries 1987 and **#9** de Jager 1974 (the ADNDT compilations), **#58** Fricke 1995 and **#40** Angeli–Marinova 2013 (charge radii), **#41** Sick 1974 (the Fourier–Bessel method), **#8** Rand–Frosch–Yearian 1966, **#22** Bergstrom–Kowalski–Neuhausen 1982, **#23** van Niftrik 1971 and **#24** Lichtenstadt 1989 (the magnetization/magnetic sector, which has no free source at all).

### 4b. The α–d cluster and the tagged channel's own observable

**This subsection has no free entry at all.** Every measurement of the ⁶Li
α–d relative-momentum distribution — the exact observable `TaggedChannel`
samples — is pre-arXiv and paywalled, which is the single sharpest gap in
the T3 tier. The two free cluster references (Lapikás–Wesseling–Wiringa and
Hotta *et al.*) are in §4a and are *nucleon* knockout, not α–d.

> **Also T3, PLEASE DOWNLOAD** (full entries in §1): **#10** Ent *et al.* 1994 NPA 578 and **#25** Ent *et al.* 1986 PRL 57 (the (e,e′d) measurement of the α–d relative-momentum distribution — the observable `Wave::radial` and `data/vmc/li6_alpha_d/li6.ad` model), **#11** Mitchell *et al.* 1991 (the α-knockout mirror), **#26** Genin *et al.* 1974, **#27** Albrecht *et al.* 1980 (a free JINR preprint was located in the reconciliation pass — see §1 #27) and **#28** Dollhopf *et al.* 1975 (the hadron- and α-induced probes, and the factor-two FWHM disagreement that is the honest band on the distribution's width), **#55** Connelly *et al.* 1998 (the reaction-mechanism systematic on all of them).

### 4c. EMC, DIS on lithium, and the nuclear-PDF baseline

- **Arneodo** *et al.* (NMC), *The structure function ratios F₂(Li)/F₂(D) and F₂(C)/F₂(D) at small x*, Nucl. Phys. B **441** (1995) 12. [arXiv:hep-ex/9504002](https://arxiv.org/abs/hep-ex/9504002) · HEPData [ins394050](https://www.hepdata.net/record/ins394050) — `held`, `refs/hep-ex_9504002.pdf`.
  **The only DIS measurement on ⁶Li in existence**, and the isoscalar one at that: 24 points of F₂(⁶Li)/F₂(D), free and machine-readable on HEPData. `EMC_VALENCE_DEPLETION_EPPS21` is today an nPDF-fit proxy for exactly this ratio; this is the measurement it should be gated against.
- **Amaudruz** *et al.* (NMC), *A re-evaluation of the nuclear structure function ratios for D, He, ⁶Li, C and Ca*, Nucl. Phys. B **441** (1995) 3. [arXiv:hep-ph/9503291](https://arxiv.org/abs/hep-ph/9503291) · HEPData [ins393377](https://www.hepdata.net/record/ins393377) — `held`, `refs/hep-ph_9503291.pdf`.
  The same ⁶Li target at lower Q², as F₂(C)/F₂(⁶Li) and F₂(Ca)/F₂(⁶Li) — the second lithium ratio, and the one that makes the first a two-point consistency check rather than a single number.
- **Gomez** *et al.* (SLAC E139), *Measurement of the A-dependence of deep-inelastic electron scattering*, SLAC-PUB-5813 (August 1993); published as Phys. Rev. D **49** (1994) 4348 · doi:[10.1103/PhysRevD.49.4348](https://doi.org/10.1103/PhysRevD.49.4348) — `downloaded 2026-09-16`, `refs/SLAC-PUB-5813.pdf` (86 pp., 5.7 MB, free from INSPIRE recid 359103).
  Cited for a **negative**: the canonical A-dependence programme has no lithium target, which is why the ⁶Li EMC baseline is an nPDF fit rather than a measurement. **Correction from this pass:** §1 #43 listed it as purchase-only; INSPIRE serves the SLAC-PUB-5813 preprint scan free, and it is now on disk (see §1 #43, amended).
- **Seely, Daniel, Gaskell, Arrington, Fomin, Solvignon** *et al.*, *New measurements of the EMC effect in very light nuclei*, Phys. Rev. Lett. **103** (2009) 202301. [arXiv:0904.4448](https://arxiv.org/abs/0904.4448) — `held`, `refs/0904.4448.pdf`.
  The ⁹Be result: the nearest measured statement that **clustering**, not average density, drives the EMC effect in a light nucleus — the empirical argument for treating ⁶Li as α+d rather than as six nucleons at mean density.
- **Arrington, Bane, Daniel, Fomin, Gaskell, Seely** *et al.*, *Measurement of the EMC effect in light and heavy nuclei*, Phys. Rev. C **104** (2021) 065203. [arXiv:2110.08399](https://arxiv.org/abs/2110.08399) — `held`, `refs/2110.08399.pdf`.
  The first ¹⁰B/¹¹B measurement, and the modern confirmation that **there is still no lithium point**.
- **Eskola, Paakkinen, Paukkunen, Salgado**, *EPPS21: a global QCD analysis of nuclear PDFs*, Eur. Phys. J. C **82** (2022) 413. [arXiv:2112.12462](https://arxiv.org/abs/2112.12462) · doi:[10.1140/epjc/s10052-022-10359-0](https://doi.org/10.1140/epjc/s10052-022-10359-0) — `held`, `refs/2112.12462.pdf`.
  `Epps21Ratio` / `EmcBaseline::Epps21` and the grid `EPPS21nlo_CT18Anlo_Li6` — the shipped unpolarized ⁶Li EMC baseline (`D-5`), whose proton baseline is the CT18A fit of §2.
- **Kovarik** *et al.*, *nCTEQ15 — global analysis of nuclear parton distributions with uncertainties in the CTEQ framework*, Phys. Rev. D **93** (2016) 085037. [arXiv:1509.00792](https://arxiv.org/abs/1509.00792) · doi:[10.1103/PhysRevD.93.085037](https://doi.org/10.1103/PhysRevD.93.085037) — `held`, `refs/1509.00792.pdf`.
  Second fit of the `T-39` four-fit EMC spread, and the source of a real **⁷Li** grid (`nCTEQ15_7_3`) — `01_generators-chain.md` `T-40` records that ⁷Li is therefore *not* PDF-less, contrary to the tree's current assumption.
- **Abdul Khalek, Ethier, Nocera, Rojo** *et al.*, *nNNPDF3.0: evidence for a modified partonic structure in heavy nuclei*, Eur. Phys. J. C **82** (2022) 507. [arXiv:2201.12363](https://arxiv.org/abs/2201.12363) · doi:[10.1140/epjc/s10052-022-10417-7](https://doi.org/10.1140/epjc/s10052-022-10417-7) — `held`, `refs/2201.12363.pdf`.
  The third fit of the same spread; `nNNPDF30_nlo_as_0118_A6_Z3` is its ⁶Li member.
- **Walt, Helenius, Vogelsang**, *Open-source QCD analysis of nuclear parton distribution functions at NLO and NNLO* (TUJU19), Phys. Rev. D **100** (2019) 096015. [arXiv:1908.03355](https://arxiv.org/abs/1908.03355) · doi:[10.1103/PhysRevD.100.096015](https://doi.org/10.1103/PhysRevD.100.096015) — `held`, `refs/1908.03355.pdf`.
  The source of the `TUJU19_*_7_3` ⁷Li grids named in `T-40`.

### 4d. Quasi-elastic (e,e′) — the radiative tail's dominant knob

- **Benhar, Day, Sick**, *Inclusive quasi-elastic electron–nucleus scattering*, Rev. Mod. Phys. **80** (2008) 189. [arXiv:nucl-ex/0603029](https://arxiv.org/abs/nucl-ex/0603029) — `held`, `refs/nucl-ex_0603029.pdf`.
  The review behind `qe_kf_gev` and `qe_suppression`, and the framing document for open item **Q9** (how much of the tail the quasi-elastic term really carries).
- **Benhar, Day, Sick**, *An archive for quasi-elastic electron–nucleus scattering data*. [arXiv:nucl-ex/0603032](https://arxiv.org/abs/nucl-ex/0603032) — `held`, `refs/nucl-ex_0603032.pdf`.
  The citation for the QES archive, whose free ⁶Li file is the dataset that would turn `RC_QE_KF_GEV` from an input into a fit residual.

> **Also T3, PLEASE DOWNLOAD** (full entries in §1): **#12** Heimlich *et al.* 1974 (the primary source of every point in the archive's ⁶Li set), **#13** Moniz *et al.* 1971 (the paper `RC_QE_KF_GEV` = 0.169 is taken from) and **#29** Whitney *et al.* 1974 (its full-length companion, with the large-angle ⁶Li data the archive does not carry).

### 4e. Moments, radii, the α core, and the ab-initio chain

- **Tilley, Cheves, Godwin, Hale, Hofmann, Kelley, Sheu, Weller** (TUNL), *Energy levels of light nuclei, A = 6*, Nucl. Phys. **A708** (2002) 3 — `held`, `refs/TUNL_A6_2002.pdf` (the A = 6 chapter manuscript, revised 5 Oct 2017; first page read in this pass).
  The evaluated ⁶Li ground-state properties the tree quotes — `LI6_QUADRUPOLE_FM2` rests on this alone, which is why §1 #56/#57 matter.
- **Pastore, Pieper, Schiavilla, Wiringa**, *Quantum Monte Carlo calculations of electromagnetic moments and transitions in A ≤ 9 nuclei with meson-exchange currents*, Phys. Rev. C **87** (2013) 035503. [arXiv:1212.3375](https://arxiv.org/abs/1212.3375) — `held`, `refs/1212.3375.pdf`.
  The GFMC electromagnetic moments; `01_li-nuclear-data.md` §7.1 quotes its Q(⁶Li) = −0.20(6) fm² against gate **C-2**. A calculation, not a measurement — the distinction the `C-2` row has to state.
- **Wiringa, Stoks, Schiavilla**, *An accurate nucleon–nucleon potential with charge-independence breaking* (AV18), Phys. Rev. C **51** (1995) 38. [arXiv:nucl-th/9408016](https://arxiv.org/abs/nucl-th/9408016) · doi:[10.1103/PhysRevC.51.38](https://doi.org/10.1103/PhysRevC.51.38) — `held`, `refs/nucl-th_9408016.pdf`.
  `[AV18]` — the interaction under **every** ANL VMC table in `data/vmc/` and under `fdeut.av18`. `PHYSICS_CHANNELS.md` itself flags that the tree used it without a specific citation; this is that citation.
- **Pudliner, Pandharipande, Carlson, Pieper, Wiringa**, *Quantum Monte Carlo calculations of nuclei with A ≤ 7*, Phys. Rev. C **56** (1997) 1720. [arXiv:nucl-th/9705009](https://arxiv.org/abs/nucl-th/9705009) — `held`, `refs/nucl-th_9705009.pdf`.
  The AV18+UIX Hamiltonian paper the `data/vmc/` files cite — part of gate **E-13**'s pedigree.
- **Forest, Pandharipande, Pieper, Wiringa, Schiavilla, Arriaga**, *Femtometer toroidal structures in nuclei*, Phys. Rev. C **54** (1996) 646. [arXiv:nucl-th/9603035](https://arxiv.org/abs/nucl-th/9603035) — `held`, `refs/nucl-th_9603035.pdf`.
  `[Forest96]` — the two-nucleon densities in ⁶Li/⁷Li and the cluster-overlap method, also part of **E-13**.
- **Wiringa, Schiavilla, Pieper, Carlson**, *Nucleon and nucleon-pair momentum distributions in A ≤ 12 nuclei*, Phys. Rev. C **89** (2014) 024305. [arXiv:1309.3794](https://arxiv.org/abs/1309.3794) · doi:[10.1103/PhysRevC.89.024305](https://doi.org/10.1103/PhysRevC.89.024305) — `held`, `refs/1309.3794.pdf` (corpus entry).
  `[Wiringa14]` — the published momentum distributions the tree's VMC tables are read against.
- **Carlson, Gandolfi, Pederiva, Pieper, Schiavilla, Schmidt, Wiringa**, *Quantum Monte Carlo methods for nuclear physics*, Rev. Mod. Phys. **87** (2015) 1067. [arXiv:1412.3081](https://arxiv.org/abs/1412.3081) — `held`, `refs/1412.3081.pdf`.
  The method review for everything under `data/vmc/`.
- **Pieper, Wiringa**, *Quantum Monte Carlo calculations of light nuclei*, Ann. Rev. Nucl. Part. Sci. **51** (2001) 53. [arXiv:nucl-th/0103005](https://arxiv.org/abs/nucl-th/0103005) — `held`, `refs/nucl-th_0103005.pdf`.
  The same at review length, and the older statement of the ⁶Li α+d picture inside VMC.
- **Camsonne, Katramatou, Olson, Sparveris** *et al.*, *JLab measurement of the ⁴He charge form factor at large momentum transfers*, Phys. Rev. Lett. **112** (2014) 132503. [arXiv:1309.5297](https://arxiv.org/abs/1309.5297) — `held`, `refs/1309.5297.pdf`.
  ⁴He above the sum-of-Gaussians fit's q range — the high-q half of any gate on `AlphaCoreSource`.

> **Also T3, PLEASE DOWNLOAD** (full entries in §1): **#42** Stone 2016 (the ADNDT quadrupole table; its free IAEA twin INDC(NDS)-0833 is a corpus entry with no file), **#56** Pyykkö 2008 (the Q(⁶Li) = −0.0806 fm² that gate **T11** explicitly tests against) and **#57** Borremans *et al.* 2005 (the measurement Stone's ⁶Li row traces to); **#31** George–Knutson 1999 (`LI6_ETA_DS_GK`, gate **C-1**) and **#30** Rodning–Knutson 1990 (the measured deuteron twin η_d = 0.0256(4), same technique and group, which `test_cluster.cpp` currently carries only through Machleidt's Table XV); **#44** Frosch *et al.* 1967 and **#45** Ottermann *et al.* 1985 (the two datasets under the free ⁴He sum-of-Gaussians parameterisation).

### 4f. The polarized ⁶Li target material

- **Abbon** *et al.* (COMPASS), *The COMPASS experiment at CERN*, Nucl. Instrum. Meth. A **577** (2007) 455. [arXiv:hep-ex/0703049](https://arxiv.org/abs/hep-ex/0703049) · doi:[10.1016/j.nima.2007.03.026](https://doi.org/10.1016/j.nima.2007.03.026) — `held`, `refs/hep-ex_0703049.pdf` (corpus entry).
  The spectrometer paper `PHYSICS_CHANNELS.md` takes the COMPASS ⁶LiD **dilution-factor** figures from — the only operating-experiment numbers the tree has for a ⁶Li target.

> **Also T3, PLEASE DOWNLOAD** (full entries in §1): **#46** Goertz–Meyer–Reicherz 2002 (the equal-spin-temperature closed form behind `spin_temperature_pzz`), **#47** Goertz *et al.* 1995 (the ⁶LiD material paper) and **#48** Goertz 2004 (where the EST relation is derived; `BENCHMARK_PLAN.md` §4 item 9). Their honest half is unchanged by anything in this pass: **P_zz of ⁶Li has never been measured by anyone.**
>
> The shared corpus additionally records five polarized-target sources with **no free copy and no PDF** — Ball *et al.* 2003 (first COMPASS ⁶LiD results, NIM A 498, 101), Bultmann *et al.* 1999 (the SLAC high-density ⁶LiD target, NIM A 425, 23), Chaumette *et al.* 1989 (Saclay irradiated ⁶LiD/⁷LiH, AIP Conf. Proc. 187, 1275), Koivuniemi *et al.* 2004 (COMPASS polarization build-up) and the JLab PR12-14-001 polarized-EMC proposal. None is cited anywhere in the LiPolGen tree; they are listed here so they are not re-searched.

---

## 5. T4 — the generators

The codes LiPolGen hadronizes with, bridges to, and is benchmarked against,
plus the radiative-correction Monte Carlos that define what the field means
by "the radiative correction". Sourced from `01_generators-chain.md`
(GC-1…GC-23), `01_radiative.md` Tier D and `01_coherent-diffraction.md` §C.

### 5a. The PYTHIA backend and the LHAup bridge

- **Bierlich, Chakraborty, Desai, Gellersen, Helenius, Ilten, Lönnblad, Mrenna, Prestel, Preuss, Sjöstrand, Skands, Utheim, Verheyen**, *A comprehensive guide to the physics and usage of PYTHIA 8.3*, SciPost Phys. Codebases **8** (2022). [arXiv:2203.11601](https://arxiv.org/abs/2203.11601) · doi:[10.21468/SciPostPhysCodeb.8](https://doi.org/10.21468/SciPostPhysCodeb.8) — `held`, `refs/2203.11601.pdf`.
  The hadronizer itself and the manual for everything in `src/pythia/` and `docs/PYTHIA_BRIDGE.md`; gate **D-1**.
- **Alwall** *et al.*, *A standard format for Les Houches Event Files*, Comput. Phys. Commun. **176** (2007) 300. [arXiv:hep-ph/0609017](https://arxiv.org/abs/hep-ph/0609017) · doi:[10.1016/j.cpc.2006.11.010](https://doi.org/10.1016/j.cpc.2006.11.010) — `held`, `refs/hep-ph_0609017.pdf`.
  The LHEF record `LhaupDis::setEvent` writes (`PYTHIA_BRIDGE.md` §6).
- **Boos** *et al.*, *Generic user process interface for event generators* (Les Houches Accord 1). [arXiv:hep-ph/0109068](https://arxiv.org/abs/hep-ph/0109068) — `held`, `refs/hep-ph_0109068.pdf`.
  Defines `HEPRUP`/`HEPEUP` and the **strategy-3** weighting convention the bridge uses — the normative document for the one interface LiPolGen must get exactly right.
- **Cabouat, Sjöstrand**, *Some dipole shower studies*, Eur. Phys. J. C **78** (2018) 226. [arXiv:1710.00391](https://arxiv.org/abs/1710.00391) · doi:[10.1140/epjc/s10052-018-5645-z](https://doi.org/10.1140/epjc/s10052-018-5645-z) — `held`, `refs/1710.00391.pdf`.
  The DIS-specific pedigree of the shower the bridge inherits (dipole recoil) — what `D-1` is actually testing when it tests "PYTHIA's DIS shower".
- **Helenius, Laulainen, Preuss**, *Multi-jet production in deep inelastic scattering with PYTHIA*, JHEP **05** (2025) 153. [arXiv:2410.20950](https://arxiv.org/abs/2410.20950) · doi:[10.1007/JHEP05(2025)153](https://doi.org/10.1007/JHEP05(2025)153) — `held`, `refs/2410.20950.pdf`.
  PYTHIA's DIS shower validated against **H1** data — the published pedigree gate `D-1` borrows rather than re-deriving.
- **Helenius, Laulainen, Preuss**, *Multiplicative matching of neutral current deep-inelastic scattering processes at next-to-leading order in PYTHIA 8* (2026). [arXiv:2605.00502](https://arxiv.org/abs/2605.00502) — `downloaded 2026-09-16`, `refs/2605.00502.pdf`.
  The NLO successor of the entry above: Born-level events are reweighted so the first shower emission follows the real matrix element, usable with both the default shower and Vincia, validated against HERA reduced cross sections. **No journal reference or DOI exists yet** — the arXiv record carries none as of 2026-09-16, so cite it as a preprint.
- **Aktas** *et al.* (H1), *Measurement and QCD analysis of the diffractive deep-inelastic scattering cross section at HERA* (the H1 2006 DPDF Fit A/B), Eur. Phys. J. C **48** (2006) 715. [arXiv:hep-ex/0606004](https://arxiv.org/abs/hep-ex/0606004) · doi:[10.1140/epjc/s10052-006-0035-3](https://doi.org/10.1140/epjc/s10052-006-0035-3) — `held`, `refs/hep-ex_0606004.pdf`.
  `[H1FitB]` — the diffractive PDF family behind `PDF:PomSet = 6`, i.e. the definition of the **T2 Pomeron tier** (`F-8`).
- **Bierlich, Gustafson, Lönnblad, Shah**, *The Angantyr model for heavy-ion collisions in PYTHIA8*, JHEP **10** (2018) 134. [arXiv:1806.10820](https://arxiv.org/abs/1806.10820) · doi:[10.1007/JHEP10(2018)134](https://doi.org/10.1007/JHEP10(2018)134) — `held`, `refs/1806.10820.pdf`.
  Documents *why* a nucleus id must never reach `Beams:idA/idB` — the published basis for the negative result recorded in `benchmarking/01` §2.2.

### 5b. Peer DIS / eA / polarized generators

- **Chang, Aschenauer, Baker, Jentsch, Lee, Tu, Yin, Zheng** *et al.*, *Benchmark eA generator for leptoproduction in high-energy lepton-nucleus collisions* (BeAGLE), Phys. Rev. D **106** (2022) 012007. [arXiv:2204.11998](https://arxiv.org/abs/2204.11998) · doi:[10.1103/PhysRevD.106.012007](https://doi.org/10.1103/PhysRevD.106.012007) — `held`, `refs/2204.11998.pdf`.
  The standard eA generator for the EIC (DPMJet-III + PyQM + FLUKA evaporation on a nuclear spectral function). LiPolGen needs it for three separable things: the A = 2 tagged-spectator path, the Ciofi–Simula n(k) and Strikman–Weiss light-front prescription it implements, and far-forward fragment spectra to set beside `E-16`. **Citation correction from `01_generators-chain.md` §2.2:** the arXiv record *does* carry the PRD DOI, so the sibling survey's "no journal-ref" note is stale.
- **Ingelman, Edin, Rathsman**, *LEPTO 6.5 — a Monte Carlo generator for deep inelastic lepton–nucleon scattering*, Comput. Phys. Commun. **101** (1997) 108. [arXiv:hep-ph/9605286](https://arxiv.org/abs/hep-ph/9605286) · doi:[10.1016/S0010-4655(96)00157-9](https://doi.org/10.1016/S0010-4655(96)00157-9) — `held`, `refs/hep-ph_9605286.pdf`.
  The DIS hard-process + Lund-fragmentation layer under **both** DJANGOH and PEPSI, and the closest published precedent for LiPolGen's own architecture (parton-level cross section handed to a string hadronizer).
- **Akushevich, Böttcher, Ryckbosch**, *RADGEN 1.0 — Monte Carlo generator for radiative events in DIS on polarized and unpolarized targets*. [arXiv:hep-ph/9906408](https://arxiv.org/abs/hep-ph/9906408) — `held`, `refs/hep-ph_9906408.pdf`.
  The POLRAD-2.0-derived radiative generator used by HERMES and COMPASS, and the code that physically ships inside **both** BeAGLE and PEPSI. It is the engine behind the radiative treatment of the only tensor-DIS datum in the world, so it is what the HERMES paper's ≈ 2 × 10⁻³ residual RC systematic actually means. Vector polarization only.
- **Jung**, *Hard diffractive scattering in high-energy ep collisions and the Monte Carlo generator RAPGAP*, Comput. Phys. Commun. **86** (1995) 147; preprint DESY 93-182 — `held`, `refs/jung1995_rapgap_desy93-182.pdf` (the free INSPIRE preprint, converted from gzipped PostScript; first page verified).
  The diffractive-DIS generator behind the official ePIC `DDIS/rapgap3.310` samples, and the closest external analogue of the T2 γ*–Pomeron tier's *hadronic final state* (`CH-9`).
- **Fucini, Scopetta, Viviani**, *Coherent deeply virtual Compton scattering off ⁴He*, Phys. Rev. C **98** (2018) 015203. [arXiv:1805.05877](https://arxiv.org/abs/1805.05877) · doi:[10.1103/PhysRevC.98.015203](https://doi.org/10.1103/PhysRevC.98.015203) — `held`, `refs/1805.05877.pdf`.
  The impulse-approximation convolution that **TOPEG** implements — the citable primary source for a generator that has no code paper.
- **Fucini, Scopetta, Viviani**, *Incoherent deeply virtual Compton scattering off ⁴He*, Phys. Rev. C **102** (2020) 065205. [arXiv:2008.11437](https://arxiv.org/abs/2008.11437) · doi:[10.1103/PhysRevC.102.065205](https://doi.org/10.1103/PhysRevC.102.065205) — `held`, `refs/2008.11437.pdf`.
  Its incoherent companion, and the paper that names the Orsay-Perugia generator in text. **Correction the pass carries:** TOPEG is *not* a tagged-deuteron generator — it is coherent DVCS off ⁴He, unpolarized, exclusive, untagged. The code the project brief was reaching for is STEG.
- **Lomnitz, Klein**, *Exclusive vector meson production at an electron-ion collider* (eSTARlight), Phys. Rev. C **99** (2019) 015203. [arXiv:1803.06420](https://arxiv.org/abs/1803.06420) · doi:[10.1103/PhysRevC.99.015203](https://doi.org/10.1103/PhysRevC.99.015203) — `held`, `refs/1803.06420.pdf`.
  `estarlight_li6_coherent()` and `estarlight_li6_q2_floors()` — the code every O5 coherent rate in the tree comes from (`D-4`).
- **Klein, Nystrand, Seger, Gorbunov, Butterworth**, *STARlight: a Monte Carlo simulation program for ultra-peripheral collisions of relativistic ions*, Comput. Phys. Commun. **212** (2017) 258. [arXiv:1607.03838](https://arxiv.org/abs/1607.03838) · doi:[10.1016/j.cpc.2016.10.016](https://doi.org/10.1016/j.cpc.2016.10.016) — `held`, `refs/1607.03838.pdf`.
  The Q² → 0 parent of eSTARlight, with a decade of RHIC/LHC UPC data behind it, and the documented **Z ≤ 6 Gaussian form-factor branch** — the published basis both for "⁶Li runs today" and for `slope_b`'s Gaussian model.
- **Toll, Ullrich**, *Exclusive diffractive processes in electron-ion collisions*, Phys. Rev. C **87** (2013) 024913. [arXiv:1211.3048](https://arxiv.org/abs/1211.3048) · doi:[10.1103/PhysRevC.87.024913](https://doi.org/10.1103/PhysRevC.87.024913) — `held`, `refs/1211.3048.pdf`.
  The method behind Sartre, and the statement that F(b) is recovered from the coherent dσ/dt — the imaging claim `slope_b` stands in for.
- **Toll, Ullrich**, *The dipole model Monte Carlo generator Sartre*, Comput. Phys. Commun. **185** (2014) 1835. [arXiv:1307.8059](https://arxiv.org/abs/1307.8059) · doi:[10.1016/j.cpc.2014.03.010](https://doi.org/10.1016/j.cpc.2014.03.010) — `held`, `refs/1307.8059.pdf`.
  Kept for its **negative**: §3.3.4 enumerates the supported nuclei and A = 6 is not among them, and its nucleus model is Woods–Saxon only — which is why `open_items/physics_literature.md` Item 7 rules Sartre out as the route to a tensor axis.
- **Toll, Ghosh, Srivastav**, *Efficient calculation of exclusive diffractive cross sections at the EIC and LHeC with the Sartre event generator* (2026). [arXiv:2606.14633](https://arxiv.org/abs/2606.14633) — `held`, `refs/2606.14633.pdf`.
  The table-speedup paper — the cost half of the "Sartre route" decision, and the paper whose abstract prompted (and whose source code then disproved) the Sartre-feasibility claim.
- **Tu, Jentsch, Baker, Zheng, Lee, Venugopalan, Hen, Higinbotham, Aschenauer, Ullrich**, *Probing short-range correlations in the deuteron via incoherent diffractive J/ψ production with spectator tagging at the EIC*, Phys. Lett. B **811** (2020) 135877. [arXiv:2005.14706](https://arxiv.org/abs/2005.14706) · doi:[10.1016/j.physletb.2020.135877](https://doi.org/10.1016/j.physletb.2020.135877) — `held`, `refs/2005.14706.pdf`.
  The incoherent, *tagged* partner of the coherent channel — what a tagged diffractive vector-meson event looks like at the EIC, and the BeAGLE precedent for `CoherentSampler` plus spectator tagging.

> **Code with no paper — do not cite literature that does not exist.** `01_generators-chain.md` §3.2 establishes by search that three of the tools the tree leans on have **no** literature record: **STEG** (`github.com/JeffersonLab/LightIonEIC`, public, no license, last push 2016-09-06 — cite the repository plus the Cosyn–Weiss theory it implements), **lAger** (`eicweb.phy.anl.gov/monte_carlo/lager`, public, active — cite the repository plus the Yellow Report), and **eic-smear** / **afterburner** (§7). **TOPEG** likewise has no code paper; cite the two Fucini–Scopetta–Viviani papers above.
>
> **Also T4, PLEASE DOWNLOAD** (full entries in §1): **#35** HERACLES (the collider-kinematics O(α) engine *inside* DJANGOH — the actual second opinion `D-3`/`C-8` want), **#36** DJANGO6 (which defines the published Rad/noRad table `BENCHMARK_PLAN.md` §4 item 7 rates a top-ten zero-dependency item) and **#37** PEPSI (the only code anywhere with a live tensor branch, `pnrun = ±2`). All three are Elsevier-era CPC articles with no arXiv entry and no INSPIRE fulltext.

### 5c. The radiative Monte Carlos

- **Akushevich, Ilyichev, Shumeiko, Soroko, Tolkachev**, *POLRAD 2.0 — FORTRAN code for the radiative corrections calculation in deep inelastic scattering of polarized particles*, Comput. Phys. Commun. **104** (1997) 201. [arXiv:hep-ph/9706516](https://arxiv.org/abs/hep-ph/9706516) — `held`, `refs/hep-ph_9706516.pdf`.
  **The single most load-bearing reference in the radiative sector**: `01_radiative.md` Table 1 lists Eqs. (9)/(10), (18), (37)–(39), (43), (44), (A.1)–(A.4) and Appendix B as the equations LiPolGen transcribes. Everything `RcModel` does in the tensor Born and elastic-tail sector comes from this paper. It was a corpus gap until this reference series; it is now on disk.
- **Afanasev, Akushevich, Ilyichev, Merenkov**, *ESFRAD* (2022). [arXiv:2212.04730](https://arxiv.org/abs/2212.04730) — `held`, `refs/2212.04730.pdf`.
  Higher-order QED against POLRAD and RADGEN — the published measure of the exponentiation `TPeakPlusLL` does not have, and the bound on how much that costs.
- **Byer, Khachatryan, Gao, Akushevich, Ilyichev** *et al.*, *SIDIS-RC EvGen: a Monte Carlo event generator for semi-inclusive DIS with radiative corrections*, Comput. Phys. Commun. **287** (2023) 108702. [arXiv:2210.03785](https://arxiv.org/abs/2210.03785) — `held`, `refs/2210.03785.pdf`.
  The maintained modern implementation of the same O(α) formalism — the "Byer" reference, and the living code to compare a transcription against.
- **Akushevich, Ilyichev, Shumeiko**, *Radiative effects in deep inelastic scattering of polarized leptons on polarized nucleons* (review, 2001). [arXiv:hep-ph/0106180](https://arxiv.org/abs/hep-ph/0106180) — `held`, `refs/hep-ph_0106180.pdf`.
  The code map — POLRAD, RADGEN, HAPRAD, DIFFRAD, MASCARAD — and the covariant formulae in one place; the document that says which of these codes does what.

> **Not fetched, deliberately:** the **DJANGOH** manual (v4.6.x) is free on the Mainz web page and `01_generators-chain.md` records that ≥ 4.6.10 does polarized NC/CC; the **POLRAD 2.0** FORTRAN source is the CPC Program Library entry `ADGH_v1_0` (a catalogue id, not a paper) and the in-tree source is what `polrad_transcription_check.md` checks against.

---

## 6. Theory inputs

Formalism the code **transcribes** rather than data it is compared against.
Sourced from `01_tensor-b1.md` (Tables 1a–1g), `01_nucleon-deuteron-data.md`
(E, F, G, I), `01_coherent-diffraction.md` (A, B, D, E, F, G) and
`01_radiative.md` (Tiers A–C).

### 6a. Spin-1 foundations

- **Kumano**, *Tensor structure function b₁(x) for spin-one hadrons*, MKPH-T-93-03 (1993). [arXiv:hep-ph/9302320](https://arxiv.org/abs/hep-ph/9302320) — `held`, `refs/hep-ph_9302320.pdf`.
  The earliest statement of the b₁–b₄ measurement case and the historical head of the Kumano series that runs through §6c — useful for dating the field's claims, not for numbers.
- **Bacchetta, Mulders**, *Deep inelastic leptoproduction of spin-one hadrons*, Phys. Rev. D **62** (2000) 114004. [arXiv:hep-ph/0007120](https://arxiv.org/abs/hep-ph/0007120) · doi:[10.1103/PhysRevD.62.114004](https://doi.org/10.1103/PhysRevD.62.114004) — `held`, `refs/hep-ph_0007120.pdf`.
  The spin-1 TMD basis (f₁LL, h⊥₁LL, …) and the Trento-convention names the modern spin-1 literature converts into — the dictionary between `TensorSF`'s labels and everyone else's.

> **The foundations themselves are all PLEASE DOWNLOAD** (full entries in §1), and this is the uncomfortable fact of the tensor sector: **#3** Hoodbhoy–Jaffe–Manohar 1989 (which *defines* b₁…b₄ and Δ and gives the (1, −2, 1) helicity pattern `TensorSF` implements), **#5** Jaffe–Manohar NPB 321 (the arbitrary-J generalisation, and the formal home of the J = 3/2 basis ⁷Li would need), **#4** Jaffe–Manohar *Nuclear gluonometry* (which names Δ an exotic-glue observable — the concept behind `toy_delta_gluon`), **#1** Close–Kumano 1990 (the ∫b₁ dx = 0 sum rule that *is* in-tree check **E-4**), and **#17**/**#18** Khan–Hoodbhoy 1991/1993. Free secondaries that quote them are on disk: `refs/1803.11206.pdf` for *Nuclear gluonometry*, and `refs/1005.4524.pdf` §I and Eq. (9) for the Close–Kumano sum rule with its antiquark correction.

### 6b. Calculations of b₁ᵈ — the camps

- **Cosyn, Dong, Kumano, Sargsian**, *Deuteron tensor structure function b₁*, Phys. Rev. D **95** (2017) 074036. [arXiv:1702.05337](https://arxiv.org/abs/1702.05337) — `held`, `refs/1702.05337.pdf` (corpus entry).
  `[CDKS17]` — the convolution calculation behind `b1_convolution`, `CdksB1` and `tables::kB1CdksQ2p5`, and therefore behind gate **E-2** and the whole A = 2 CDKS gate. Its unpolarized PDF is MSTW2008 LO (§2), which is why that PDF is part of the gate's definition.
- **Miller**, *Pionic and hidden-color, six-quark contributions to the deuteron b₁ structure function*, Phys. Rev. C **89** (2014) 045203. [arXiv:1311.4561](https://arxiv.org/abs/1311.4561) — `held`, `refs/1311.4561.pdf` (corpus entry).
  `[Miller14]` — the second camp, behind `toy_b1`, `MillerB1`, `tables::kB1Miller` and gate **E-3**.
- **Nikolaev, Schäfer**, *Nonvanishing tensor polarization of sea quarks in polarized deuterons*, Phys. Lett. B **398** (1997) 245; Erratum PLB **407** (1997) 453. [arXiv:hep-ph/9611460](https://arxiv.org/abs/hep-ph/9611460) · doi:[10.1016/S0370-2693(97)00250-5](https://doi.org/10.1016/S0370-2693(97)00250-5) — `held`, `refs/hep-ph_9611460.pdf`.
  The shadowing camp, and the reason the small-x end of the tree's b₁ is not a settled number: `01_tensor-b1.md` §3.1 records that it predicts b₂ rising at small x and A₂ ≈ 1 %, **two orders of magnitude above the impulse approximation**, with an explicit claim that the Close–Kumano sum rule is broken.
- **Bora, Jaffe**, *The double-scattering contribution to b₁(x, Q²) in the deuteron*, Phys. Rev. D **57** (1998) 6906. [arXiv:hep-ph/9711323](https://arxiv.org/abs/hep-ph/9711323) · doi:[10.1103/PhysRevD.57.6906](https://doi.org/10.1103/PhysRevD.57.6906) — `held`, `refs/hep-ph_9711323.pdf`.
  The VMD/double-scattering camp, with the **opposite** small-x verdict to Nikolaev–Schäfer: b₁ → 0 as x → 0 by rotational symmetry. The disagreement between these two is real and the tree does not currently represent it.
- **Edelmann, Piller, Weise**, *Polarized deuteron structure functions at small x*, Z. Phys. A **357** (1997) 129. [arXiv:nucl-th/9701026](https://arxiv.org/abs/nucl-th/9701026) · doi:[10.1007/s002180050226](https://doi.org/10.1007/s002180050226) — `held`, `refs/nucl-th_9701026.pdf`.
  A third shadowing calculation: b₁ "surprisingly large at x < 0.1", dominated by coherent double scattering.
- **Edelmann, Piller, Weise**, *Deuteron spin structure functions at small Bjorken x*, Phys. Rev. C **57** (1998) 3392. [arXiv:hep-ph/9709455](https://arxiv.org/abs/hep-ph/9709455) · doi:[10.1103/PhysRevC.57.3392](https://doi.org/10.1103/PhysRevC.57.3392) — `held`, `refs/hep-ph_9709455.pdf`.
  Its long version — spin-extended Glauber–Gribov multiple scattering worked out in full.
- **Umnikov**, *Relativistic calculation of the structure functions b₁,₂(x) of the deuteron*, Phys. Lett. B **391** (1997) 177. [arXiv:hep-ph/9605291](https://arxiv.org/abs/hep-ph/9605291) · doi:[10.1016/S0370-2693(96)01440-2](https://doi.org/10.1016/S0370-2693(96)01440-2) — `held`, `refs/hep-ph_9605291.pdf`.
  **The validity floor on `b1_convolution`**: a non-relativistic convolution gets small x wrong and violates the exact sum rules. Any statement about where the tree's b₁ may be trusted should cite this, not a camp.
- **Kumano, Kuroki**, *Tensor-polarized parton distribution functions of the deuteron by a convolution model* (2026). [arXiv:2607.09237](https://arxiv.org/abs/2607.09237) — `held`, `refs/2607.09237.pdf`.
  An independent second convolution calculation **at exactly Q² = 2.5 GeV²** — the same point as `tables::kB1CdksQ2p5`, which makes it a directly comparable external check rather than another camp.

### 6c. Tensor-polarized PDFs, sum rules and the modern spin-1 basis

- **Kumano**, *Tensor-polarized quark and antiquark distribution functions in a spin-one hadron*, Phys. Rev. D **82** (2010) 017501. [arXiv:1005.4524](https://arxiv.org/abs/1005.4524) · doi:[10.1103/PhysRevD.82.017501](https://doi.org/10.1103/PhysRevD.82.017501) — `held`, `refs/1005.4524.pdf`.
  `01_tensor-b1.md` §4.1 calls it "the only b₁ source in this domain that publishes numbers": a closed-form δ_T w(x) with fitted parameters at Q² = 2.5 GeV², fitted to HERMES. A **data-driven third camp** for `TensorSF`, and the free secondary that states the Close–Kumano sum rule with its antiquark correction.
- **Kumano, Song**, *Twist-2 relation and sum rule for tensor-polarized parton distribution functions of spin-1 hadrons*, JHEP **09** (2021) 141. [arXiv:2106.15849](https://arxiv.org/abs/2106.15849) · doi:[10.1007/JHEP09(2021)141](https://doi.org/10.1007/JHEP09(2021)141) — `held`, `refs/2106.15849.pdf`.
  **A second free, closed-form sum rule beside Close–Kumano** — ∫dx f₂LT = 0 with f₂LT = (2/3)f_LT − f₁LL, plus a Wandzura–Wilczek-type twist-2 relation. Since the Close–Kumano paper itself is paywalled, this is the cheapest new sum-rule gate the project can actually read.
- **Kumano, Song**, *TMD parton distribution functions up to twist 4 for spin-1 hadrons*, Phys. Rev. D **103** (2021) 014025. [arXiv:2011.08583](https://arxiv.org/abs/2011.08583) — `held`, `refs/2011.08583.pdf`.
  The complete spin-1 correlator decomposition — 40 TMDs over twists 2–4 — i.e. the map of everything a spin-1 target *can* carry, against which the tree implements a handful.
- **Kumano, Song**, *Theoretical estimate on tensor-polarization asymmetry in the proton–deuteron Drell–Yan process*, Phys. Rev. D **94** (2016) 054022. [arXiv:1606.03149](https://arxiv.org/abs/1606.03149) — `held`, `refs/1606.03149.pdf`.
  The Fermilab route to δ_T q̄ using the 2010 PDFs above — the non-DIS way the same distributions could be measured.
- **Kumano**, *Tensor-polarized structure functions: tensor structure of deuteron in 2020's*, J. Phys. Conf. Ser. **543** (2014) 012001. [arXiv:1407.3852](https://arxiv.org/abs/1407.3852) · doi:[10.1088/1742-6596/543/1/012001](https://doi.org/10.1088/1742-6596/543/1/012001) — `held`, `refs/1407.3852.pdf`.
  Carries **closed-form projection operators** extracting b₁–b₄ from the hadron tensor W_μν — a free, independent check on `InclusiveKernel`'s transcription of the Cosyn *et al.* decomposition, which is otherwise checked only against itself.
- **Kumano**, *Parton distribution functions and fragmentation functions of spin-1 hadrons*, Eur. Phys. J. A **60** (2024) 205. [arXiv:2406.01180](https://arxiv.org/abs/2406.01180) · doi:[10.1140/epja/s10050-024-01411-6](https://doi.org/10.1140/epja/s10050-024-01411-6) — `held`, `refs/2406.01180.pdf`.
  The single best entry point to the whole spin-1 sector — the review to cite once instead of six primaries, when a primary is not what is wanted.
- **Zhao, Bacchetta, Kumano, Liu, Zhou**, *Semi-inclusive deep inelastic scattering off a tensor-polarized spin-1 target*, JHEP **12** (2025) 067. [arXiv:2508.06134](https://arxiv.org/abs/2508.06134) · doi:[10.1007/JHEP12(2025)067](https://doi.org/10.1007/JHEP12(2025)067) — `held`, `refs/2508.06134.pdf`.
  23 structure functions, 21 non-vanishing at tree level through twist 3 — the SIDIS counterpart of the inclusive decomposition the tree implements, and the map of what a SIDIS extension would owe.

### 6d. The double-helicity-flip gluon Δ ("exotic glue")

- **Detmold, Shanahan**, *Gluonic transversity from lattice QCD*, Phys. Rev. D **94** (2016) 014507; Erratum PRD **95** (2017) 079902. [arXiv:1606.04505](https://arxiv.org/abs/1606.04505) · doi:[10.1103/PhysRevD.94.014507](https://doi.org/10.1103/PhysRevD.94.014507) — `held`, `refs/1606.04505.pdf`.
  The lattice half of `CoherentScenario::amp`'s 3 × 10⁻³ … 1 × 10⁻² band. `01_tensor-b1.md` §5.1 records (PDF-read) A₂ ≈ 0.23(2)(5) in the **φ meson**, with the gluonic Soffer bound saturated at ≈ 80 % — a meson, not a nucleus, which is exactly why the next entry matters.
- **Winter, Detmold, Gambhir, Orginos, Savage, Shanahan, Wagman** (NPLQCD), *First lattice QCD study of the gluonic structure of light nuclei*, Phys. Rev. D **96** (2017) 094512. [arXiv:1709.00395](https://arxiv.org/abs/1709.00395) · doi:[10.1103/PhysRevD.96.094512](https://doi.org/10.1103/PhysRevD.96.094512) — `held`, `refs/1709.00395.pdf`.
  **The first exotic-glue number in a nucleus.** `01_tensor-b1.md` §5.2 records (PDF-read) a₂^(d) = −0.010(3) bare at m_π ≈ 806 MeV against ⟨x⟩_g^(d) = 0.51(5), and c₂^(d)/b₂^(d) ≲ 1/20 for the gluonic-b₁ combination. This is the reference that turns `toy_delta_gluon`'s scenario band from a guess into a bracketed one.
- **Cotogno, van Daal, Mulders**, *Positivity bounds on gluon TMDs for hadrons of spin ≤ 1*, JHEP **11** (2017) 185. [arXiv:1709.07827](https://arxiv.org/abs/1709.07827) · doi:[10.1007/JHEP11(2017)185](https://doi.org/10.1007/JHEP11(2017)185) — `held`, `refs/1709.07827.pdf`.
  **A free positivity gate on `toy_delta_gluon`**: leading-twist gluon-distribution bounds including tensor polarization, with a small-x limit. Cheap to implement and entirely external.
- **Kumano, Song**, *Gluon transversity in polarized proton–deuteron Drell–Yan processes*, Phys. Rev. D **101** (2020) 054011. [arXiv:1910.12523](https://arxiv.org/abs/1910.12523) · doi:[10.1103/PhysRevD.101.054011](https://doi.org/10.1103/PhysRevD.101.054011) — `held`, `refs/1910.12523.pdf`.
  The Drell–Yan half of the same band — the "Drell–Yan estimates" `coherent.hpp`'s bound comment refers to.
- **Kumano, Song**, *Deuteron polarizations in the proton–deuteron Drell–Yan process for finding the gluon transversity*, Phys. Rev. D **101** (2020) 094013. [arXiv:2003.06623](https://arxiv.org/abs/2003.06623) — `held`, `refs/2003.06623.pdf`.
  The same, re-expressed in conventional (measurable) deuteron polarizations.

> **Also a theory input, PLEASE DOWNLOAD** (full entry in §1): **#6** Sather–Schmidt 1990 — the actual source of `C_BAG` = −0.012, a constant the generator ships and cannot currently re-derive.

### 6e. Spin-3/2 — the ⁷Li rank-2 sector

- **Fu, Sun, Dong**, *Generalized parton distributions in a spin-3/2 particle*, Phys. Rev. D **106** (2022) 116012. [arXiv:2209.12161](https://arxiv.org/abs/2209.12161) — `held`, `refs/2209.12161.pdf`.
  `[Fu22]` — the J = 3/2 forward-limit structure functions, i.e. the b₁ definition ⁷Li would be normalised against. **Naming trap recorded in `01_tensor-b1.md` §6.1:** their "g₂" is a rank-3 object and collides with the twist-3 nucleon g₂.
- **Fu, Dong, Kumano**, *Transversity generalized parton distributions in a spin-3/2 particle*, Phys. Rev. D **109** (2024) 096006. [arXiv:2402.11561](https://arxiv.org/abs/2402.11561) — `held`, `refs/2402.11561.pdf`.
  16 transversity GPDs per parton — the gluon-transversity analogue of Δ for J = 3/2.
- **Fu, Dong, Kumano, Xie**, *Generalizing the Soffer bound: positivity constraints on parton distribution functions of spin-3/2 particles*, Phys. Rev. D **113** (2026) L111901. [arXiv:2602.11587](https://arxiv.org/abs/2602.11587) · doi:[10.1103/qq1h-snhk](https://doi.org/10.1103/qq1h-snhk) — `held`, `refs/2602.11587.pdf`.
  The complete set of spin-3/2 positivity bounds, for the first time. `01_tensor-b1.md` §8 rates it **the cheapest new external gate available to this project**, and it is the only thing in the corpus that constrains the ⁷Li rank-2 block — which is identically zero today.

### 6f. Nuclear spin structure, tagging, FSI and the wave functions

- **Berger, Cano, Diehl, Pire**, *Generalized parton distributions in the deuteron*, Phys. Rev. Lett. **87** (2001) 142302. [arXiv:hep-ph/0106192](https://arxiv.org/abs/hep-ph/0106192) · doi:[10.1103/PhysRevLett.87.142302](https://doi.org/10.1103/PhysRevLett.87.142302) — `held`, `refs/hep-ph_0106192.pdf`.
  The spin-1 GPD basis — 5 quark and 9 gluon distributions, including the **gluon-transversity GPD** — i.e. the exclusive counterpart of the coherent ⁶Li channel and the formal home of the object `CoherentScenario::amp` parameterises.
- **Cosyn, Sargsian**, *Final-state interactions in DIS from a tensor-polarized deuteron target*, J. Phys. Conf. Ser. **543** (2014) 012006. [arXiv:1407.1653](https://arxiv.org/abs/1407.1653) · doi:[10.1088/1742-6596/543/1/012006](https://doi.org/10.1088/1742-6596/543/1/012006) — `held`, `refs/1407.1653.pdf`.
  **The only calculation anywhere of A_zz with any non-Born effect in it.** `01_tensor-b1.md` §1g records its two statements: inclusively, FSI are sizeable for x > 0.2 but A_zz stays small; in the tagged channel, FSI are largest at p_spec ≈ 300 MeV and forward angles. That is the external shape `GlauberFsiWeight`'s band should be argued against.
- **Melnitchouk, Sargsian, Strikman**, *Probing the origin of the EMC effect via tagged structure functions of the deuteron*, Z. Phys. A **359** (1997) 99. [arXiv:nucl-th/9609048](https://arxiv.org/abs/nucl-th/9609048) · doi:[10.1007/s002180050372](https://doi.org/10.1007/s002180050372) — `held`, `refs/nucl-th_9609048.pdf`.
  The pole-extrapolation logic under `TaggedModel`'s backward-spectator argument — why tagging a slow backward spectator selects a nearly-free nucleon at all.
- **Sargsian**, *Large Q² electrodisintegration of the deuteron in the virtual-nucleon approximation*, Phys. Rev. C **82** (2010) 014612. [arXiv:0910.2016](https://arxiv.org/abs/0910.2016) · doi:[10.1103/PhysRevC.82.014612](https://doi.org/10.1103/PhysRevC.82.014612) — `held`, `refs/0910.2016.pdf`.
  The generalized eikonal approximation that `GlauberFsiWeight` is a cheap approximation *of* — the calculation the tree's FSI weight owes its form to.
- **Cosyn, Sargsian**, *Nuclear final-state interactions in deep inelastic scattering off the lightest nuclei*, Int. J. Mod. Phys. E **26** (2017) 1730004. [arXiv:1704.06117](https://arxiv.org/abs/1704.06117) · doi:[10.1142/S0218301317300041](https://doi.org/10.1142/S0218301317300041) — `held`, `refs/1704.06117.pdf`.
  `[CS17]` — the review carrying the σ_XN(W) Deeps fit that `fsi.hpp` uses as its Eq. (11). The numbers in the code come from here.
- **Machleidt**, *The high-precision, charge-dependent Bonn nucleon–nucleon potential (CD-Bonn)*, Phys. Rev. C **63** (2001) 024001. [arXiv:nucl-th/0006014](https://arxiv.org/abs/nucl-th/0006014) · doi:[10.1103/PhysRevC.63.024001](https://doi.org/10.1103/PhysRevC.63.024001) — `held`, `refs/nucl-th_0006014.pdf`. **Corpus gap G8, closed.**
  `[CDBonn01]` — the direct source of `CdBonnWave`'s twenty typed coefficients (Table XX) and of the published B_d, A_S, η, P_D, Q_d (Table XV). `00_corpus.md` called it "the single most load-bearing reference missing from the corpus"; it is now on disk.
- **Ciofi degli Atti, Simula**, *Realistic model of the nucleon spectral function in few- and many-nucleon systems*, Phys. Rev. C **53** (1996) 1689. [arXiv:nucl-th/9507024](https://arxiv.org/abs/nucl-th/9507024) · doi:[10.1103/PhysRevC.53.1689](https://doi.org/10.1103/PhysRevC.53.1689) — `held`, `refs/nucl-th_9507024.pdf`.
  `CiofiSimulaTriton` and the `cs_n0_terms` / `cs_n1_terms` coefficients (Table A.1, Eqs. 75–76). `01_nucleon-deuteron-data.md` labels it honestly for what it is: a **fit to a calculation**, not to data.
- **Garçon, Van Orden**, *The deuteron: structure and form factors*, Adv. Nucl. Phys. **26** (2001) 293. [arXiv:nucl-th/0102049](https://arxiv.org/abs/nucl-th/0102049) · doi:[10.1007/0-306-47915-X_4](https://doi.org/10.1007/0-306-47915-X_4) — `held`, `refs/nucl-th_0102049.pdf`.
  The review that fixes conventions — G_C/G_M/G_Q against A/B/t₂₀, spherical against Cartesian. The document to settle a sign argument with before it starts.
- **Cloët, Bentz, Thomas**, *EMC and polarized EMC effects in nuclei*, Phys. Lett. B **642** (2006) 210. [arXiv:nucl-th/0605061](https://arxiv.org/abs/nucl-th/0605061) · doi:[10.1016/j.physletb.2006.08.076](https://doi.org/10.1016/j.physletb.2006.08.076) — `held`, `refs/nucl-th_0605061.pdf` (corpus entry).
  `[CBT06]` — nuclear structure functions and quark distributions computed for **⁷Li**, ¹¹B, ¹⁵N and ²⁷Al in a relativistic shell model with a confining NJL nucleon, and the prediction of a polarized EMC effect larger than the unpolarized one. **Citation correction from this pass** (verified on INSPIRE, recid for nucl-th/0605061): the corpus entry records this arXiv id under the title *Spin-dependent structure functions in nuclear matter and the polarized EMC effect*, which is a **different** Cloët–Bentz–Thomas paper (PRL 95 (2005) 052302). PLB 642 (2006) 210 is *EMC and polarized EMC effects in nuclei*.
- **Tronchin, Matevosyan, Thomas**, *Polarized EMC effect in the QMC model*, Phys. Lett. B **783** (2018) 247. [arXiv:1806.00481](https://arxiv.org/abs/1806.00481) · doi:[10.1016/j.physletb.2018.06.065](https://doi.org/10.1016/j.physletb.2018.06.065) — `held`, `refs/1806.00481.pdf` (corpus entry).
  `[TMT18]` — the same question in the quark–meson-coupling model, with the free-nucleon case from the MIT bag. The second of the two polarized-EMC predictions the tree's `00_in_tree_checks.md` weighs.
- **Wang, Bentz, Cloët, Thomas**, *Gluon EMC effects in nuclear matter*, J. Phys. G **49** (2022) 03LT01. [arXiv:2109.03591](https://arxiv.org/abs/2109.03591) · doi:[10.1088/1361-6471/ac4c90](https://doi.org/10.1088/1361-6471/ac4c90) — `held`, `refs/2109.03591.pdf` (corpus entry).
  The gluonic structure of nuclei in the same mean-field/NJL framework — the gluon-spin arm of `plans/02` step 1.2.2. **Citation correction:** the corpus records the title as *Polarized gluon EMC effect*; the published title (verified on INSPIRE) is *Gluon EMC effects in nuclear matter*.
- **Bacchetta, Diehl, Goeke, Metz, Mulders, Schlegel**, *Semi-inclusive deep inelastic scattering at small transverse momentum*, JHEP **02** (2007) 093. [arXiv:hep-ph/0611265](https://arxiv.org/abs/hep-ph/0611265) · doi:[10.1088/1126-6708/2007/02/093](https://doi.org/10.1088/1126-6708/2007/02/093) — `held`, `refs/hep-ph_0611265.pdf` (corpus entry).
  `[Bacchetta07]` — the tree-level structure-function decomposition of SIDIS at leading and first subleading twist. Cited in LiPolGen's code for the azimuthal-modulation conventions; `PHYSICS_CHANNELS.md` flags only that it is unregistered in the tree's own literature list, not that it is missing.
- **Bacchetta, D'Alesio, Diehl, Miller**, *Single-spin asymmetries: the Trento conventions*, Phys. Rev. D **70** (2004) 117504. [arXiv:hep-ph/0410050](https://arxiv.org/abs/hep-ph/0410050) · doi:[10.1103/PhysRevD.70.117504](https://doi.org/10.1103/PhysRevD.70.117504) — `held`, `refs/hep-ph_0410050.pdf` (corpus entry).
  The definitions and notations `CONVENTIONS.md` binds the tree to — the reason a φ in LiPolGen means the same thing as a φ in anyone else's paper.
- **Diehl, Sapeta**, *On the analysis of lepton scattering on longitudinally or transversely polarized protons*, Eur. Phys. J. C **41** (2005) 515. [arXiv:hep-ph/0503023](https://arxiv.org/abs/hep-ph/0503023) · doi:[10.1140/epjc/s2005-02242-9](https://doi.org/10.1140/epjc/s2005-02242-9) — `held`, `refs/hep-ph_0503023.pdf` (corpus entry).
  The companion that treats target polarization defined relative to the **lepton beam** against the **virtual photon** — the distinction that changes azimuthal distributions, and the one `CONVENTIONS.md` has to keep straight for `TENSOR_LL_SIGN`.

> **Also theory inputs, PLEASE DOWNLOAD** (full entries in §1): **#30** Rodning–Knutson 1990 (the measured η_d behind `test_cluster.cpp`) — see §4e.

### 6g. Coherent diffraction

- **Guzey, Rinaldi, Scopetta, Strikman, Viviani**, *Coherent J/ψ electroproduction on ⁴He and ³He at the Electron-Ion Collider: probing nuclear shadowing one nucleon at a time*, Phys. Rev. Lett. **129** (2022) 242503. [arXiv:2202.12200](https://arxiv.org/abs/2202.12200) · doi:[10.1103/PhysRevLett.129.242503](https://doi.org/10.1103/PhysRevLett.129.242503) — `held`, `refs/2202.12200.pdf`.
  The closest published thing to a ⁶Li coherent amplitude: a Gribov–Glauber expansion on ⁴He and ³He with one-, two- and three-body form factors built from **realistic AV18 few-body wave functions** — the same class the tree already carries. `01_coherent-diffraction.md` §B1 records two results that ride with any use of it: the coherent |t| distribution shifts to smaller |t|, and the ⁴He **diffractive minimum is not at the charge-form-factor minimum** — a direct warning against `gaussian_slope(r_rms_fm)` beyond the first lobe, and part of why `COHERENT_T_MAX_DEFAULT` = 0.2 is the right ceiling.
- **Mäntysaari, Roch, Schenke, Shen, Zhao**, *Nuclear structure and saturation effects from diffractive vector meson production*, Phys. Rev. D **114** (2026) 014068. [arXiv:2605.00454](https://arxiv.org/abs/2605.00454) · doi:[10.1103/2628-2spx](https://doi.org/10.1103/2628-2spx) — `held`, `refs/2605.00454.pdf`.
  Retires the tree's "lightest published nucleus is Ca" note: it generates initial-state nucleon configurations from VMC, PGCM (with *clustering* and *uniform* sampling), NLEFT and **GFMC (used for ³He and ⁴He)** inside a Good–Walker amplitude. Two uses: it shows the machinery already ingests ab-initio configuration files, so supplying an α+d ⁶Li file is a small ask; and its clustering-vs-uniform comparison is the closest published measure of how much a clustered sampling moves a coherent |t| distribution. Unpolarized, so no tensor axis.
- **Mondal, Kumar, Sarkar**, *Imprints of nuclear shell structure in exclusive vector meson production* (2026). [arXiv:2608.23445](https://arxiv.org/abs/2608.23445) — `held`, `refs/2608.23445.pdf`.
  How far a realistic light-nucleus density moves the coherent |t| shape — a named systematic on `slope_b` beyond the Gaussian.
- **Mäntysaari, Schenke, Shen, Zhao**, *Multiscale imaging of nuclear deformation at the Electron-Ion Collider*, Phys. Rev. Lett. **131** (2023) 062301. [arXiv:2303.04866](https://arxiv.org/abs/2303.04866) · doi:[10.1103/PhysRevLett.131.062301](https://doi.org/10.1103/PhysRevLett.131.062301) — `held`, `refs/2303.04866.pdf`.
  Deformation → coherent/incoherent ratio: the *unpolarized* half of the `delta_b_m` geometry, and the precedent for reading a shape off a diffractive pattern.
- **Mäntysaari, Salazar, Schenke**, *Nuclear geometry at high energy from exclusive vector meson production*, Phys. Rev. D **106** (2022) 074019. [arXiv:2207.03712](https://arxiv.org/abs/2207.03712) · doi:[10.1103/PhysRevD.106.074019](https://doi.org/10.1103/PhysRevD.106.074019) — `held`, `refs/2207.03712.pdf`.
  Strong-interaction radius against charge radius — a named systematic on `gaussian_slope(r_rms_fm)`, which takes the charge radius and uses it as a gluonic one.
- **Mäntysaari, Roy, Salazar, Schenke**, *Gluon imaging using azimuthal correlations in diffractive scattering at the Electron-Ion Collider*, Phys. Rev. D **103** (2021) 094026. [arXiv:2011.02464](https://arxiv.org/abs/2011.02464) · doi:[10.1103/PhysRevD.103.094026](https://doi.org/10.1103/PhysRevD.103.094026) — `held`, `refs/2011.02464.pdf`.
  The *unpolarized* e–V azimuthal modulation — a third cos-type background under `cos2phi_coefficient`, and the size any tensor signal has to beat.
- **Mäntysaari, Schenke**, *Revealing proton shape fluctuations with incoherent diffraction at high energy*, Phys. Rev. D **94** (2016) 034042. [arXiv:1607.01711](https://arxiv.org/abs/1607.01711) · doi:[10.1103/PhysRevD.94.034042](https://doi.org/10.1103/PhysRevD.94.034042) — `held`, `refs/1607.01711.pdf`.
  The paper behind `subnucleondiffraction`, the code `write_snd_configs()` writes configuration files for.
- **Kowalski, Teaney**, *An impact parameter dipole saturation model*, Phys. Rev. D **68** (2003) 114005. [arXiv:hep-ph/0304189](https://arxiv.org/abs/hep-ph/0304189) · doi:[10.1103/PhysRevD.68.114005](https://doi.org/10.1103/PhysRevD.68.114005) — `held`, `refs/hep-ph_0304189.pdf`.
  IPsat/bSat itself — Sartre's `bSat`, and the b-dependence any ⁶Li amplitude would inherit.
- **Kowalski, Motyka, Watt**, *Exclusive diffractive processes at HERA within the dipole picture*, Phys. Rev. D **74** (2006) 074016. [arXiv:hep-ph/0606272](https://arxiv.org/abs/hep-ph/0606272) · doi:[10.1103/PhysRevD.74.074016](https://doi.org/10.1103/PhysRevD.74.074016) — `held`, `refs/hep-ph_0606272.pdf`.
  The bSat fit Sartre cites, with the vector-meson wave functions and t-slopes.
- **Watt, Kowalski**, *Impact parameter dependent colour glass condensate dipole model*, Phys. Rev. D **78** (2008) 014016. [arXiv:0712.2670](https://arxiv.org/abs/0712.2670) · doi:[10.1103/PhysRevD.78.014016](https://doi.org/10.1103/PhysRevD.78.014016) — `held`, `refs/0712.2670.pdf`.
  Sartre's `DipoleModel_bCGC` — the model-systematic partner of the two above.
- **Rezaeian, Siddikov, Van de Klundert, Venugopalan**, *Analysis of combined HERA data in the impact-parameter dependent saturation model*, Phys. Rev. D **87** (2013) 034002. [arXiv:1212.2974](https://arxiv.org/abs/1212.2974) · doi:[10.1103/PhysRevD.87.034002](https://doi.org/10.1103/PhysRevD.87.034002) — `held`, `refs/1212.2974.pdf`.
  The current IPsat parameters from combined H1+ZEUS data — what a new ⁶Li dipole table would be generated with.
- **Frankfurt, Guzey, Strikman**, *Leading twist nuclear shadowing phenomena in hard processes with nuclei*, Phys. Rept. **512** (2012) 255. [arXiv:1106.2091](https://arxiv.org/abs/1106.2091) · doi:[10.1016/j.physrep.2011.12.002](https://doi.org/10.1016/j.physrep.2011.12.002) — `held`, `refs/1106.2091.pdf`.
  **The missing `f0` calculation class.** `coherent.hpp` states in prose that a coherent-A analogue of the H1/ZEUS diffractive PDFs "exists in neither eSTARlight nor Sartre"; it exists here. `01_coherent-diffraction.md` §E1 records §6.1's nuclear diffractive PDFs, **Eq. (189)**'s probability of diffraction per parton flavour, and **Fig. 69**'s x dependence at Q² = 4 GeV² for ⁴⁰Ca and ²⁰⁸Pb. Two caveats ride with it: the lightest nucleus treated there is ⁴⁰Ca, so ⁶Li needs an A-extrapolation, and P_diff is *total* diffraction, with §6.1.1–6.1.2 supplying the coherent/incoherent split. It would replace the scenario band {0.02, 0.08} that `CoherentScenario::f0` currently ships.
- **Frankfurt, Strikman, Weiss**, *Small-x physics: from HERA to LHC and beyond*, Ann. Rev. Nucl. Part. Sci. **55** (2005) 403. [arXiv:hep-ph/0507286](https://arxiv.org/abs/hep-ph/0507286) · doi:[10.1146/annurev.nucl.53.041002.110615](https://doi.org/10.1146/annurev.nucl.53.041002.110615) — `held`, `refs/hep-ph_0507286.pdf`.
  The review behind the framing that coherent diffraction images transverse gluon structure. Use it for framing; take numbers from the entries above, per the project's own rule that a primary beats a review.
- **Caldwell, Kowalski**, *The J/ψ way to nuclear structure* (preprint of *Investigating the gluonic structure of nuclei via J/ψ scattering*, Phys. Rev. C **81** (2010) 025203). [arXiv:0909.1254](https://arxiv.org/abs/0909.1254) — `held`, `refs/0909.1254.pdf`; the refereed version has **no** e-print, doi:[10.1103/PhysRevC.81.025203](https://doi.org/10.1103/PhysRevC.81.025203) (not open access).
  The |t| → b Fourier–Bessel inversion that gives `slope_b` its meaning as an imaging observable rather than a fit parameter — the argument `cluster_config.hpp`'s `eps_b0_equivalent()` makes when it turns a configuration ensemble into a slope.
- **Dalton, Deur, Keith**, *Potential for tensor polarized deuterons in Hall D at Jefferson Lab*, Eur. Phys. J. A **61** (2025) 111. [arXiv:2504.21177](https://arxiv.org/abs/2504.21177) · doi:[10.1140/epja/s10050-025-01580-y](https://doi.org/10.1140/epja/s10050-025-01580-y) — `held`, `refs/2504.21177.pdf`.
  The **only** other proposal anywhere to measure a tensor observable in coherent vector-meson production, and the m = 0 / intermediate-|t| sensitivity argument behind `a2_m_state`.
- **Abdallah** *et al.* (STAR), *Tomography of ultra-relativistic nuclei with polarized photon-gluon collisions*, Sci. Adv. **9** (2023) eabq3903. [arXiv:2204.01625](https://arxiv.org/abs/2204.01625) · doi:[10.1126/sciadv.abq3903](https://doi.org/10.1126/sciadv.abq3903) — `held`, `refs/2204.01625.pdf` (Sci. Adv. is open access).
  The **measured** size of the photon-polarization cos 2φ background — mechanism (ii) in `physics_literature.md` §7's warning, and the reason `tensor_flip_plan` exists at all.
- **Xing, Zhang, Zhou, Zhou**, *The cos 2φ azimuthal asymmetry in ρ⁰ meson production in ultraperipheral heavy ion collisions*, JHEP **10** (2020) 064. [arXiv:2006.06206](https://arxiv.org/abs/2006.06206) · doi:[10.1007/JHEP10(2020)064](https://doi.org/10.1007/JHEP10(2020)064) — `held`, `refs/2006.06206.pdf`.
  The dipole-model calculation of the same linear-polarization cos 2φ — the background amplitude `phase_C_numbers.md` §C2.5 bounds.
- **Zha, Brandenburg, Ruan, Tang, Xu**, *Exploring the double-slit interference with linearly polarized photons*, Phys. Rev. D **103** (2021) 033007. [arXiv:2006.12099](https://arxiv.org/abs/2006.12099) · doi:[10.1103/PhysRevD.103.033007](https://doi.org/10.1103/PhysRevD.103.033007) — `held`, `refs/2006.12099.pdf`.
  The interference reading of the same modulation — the alternative that makes the *two-nucleus* part of the STAR result inapplicable to e+A, which is what keeps the O5 background estimate honest.
- **Hagiwara, Zhang, Zhou, Zhou**, *Probing the gluon tomography in photoproduction of dipions*, Phys. Rev. D **104** (2021) 094021. [arXiv:2106.13466](https://arxiv.org/abs/2106.13466) · doi:[10.1103/PhysRevD.104.094021](https://doi.org/10.1103/PhysRevD.104.094021) — `held`, `refs/2106.13466.pdf`.
  cos 4φ and the elliptic-gluon vs QED-radiation splitting — the higher harmonic a ρ⁰ control measurement has to survive.

> **Also theory inputs, PLEASE DOWNLOAD** (full entries in §1): **#49** Good–Walker 1960 (the eigenstate decomposition `write_snd_configs` / `ClusterConfigSampler` exist to average over), **#50** Miettinen–Pumplin 1978 (coherent = ⟨A⟩², total diffractive = ⟨A²⟩, so incoherent = the variance — the sentence `ClusterConfigSet` implements; a free scan is served by INSPIRE, link in §1) and **#51** Rachek *et al.* 2017 (precedent only: the sole *measured* tensor observable in a coherent reaction on a nucleus).

### 6h. Radiative-correction formalism

- **Mo, Tsai**, *Radiative corrections to elastic and inelastic ep and μp scattering*, SLAC-PUB-380 (January 1968); published as Rev. Mod. Phys. **41** (1969) 205 — `held`, `refs/SLAC-PUB-0380.pdf` (the free preprint; INSPIRE recid 52657 serves it).
  `[MT69]` — Section III's exact elastic radiative tail, Appendix B's exact bremsstrahlung expression, and the paper's own exact-vs-peaking comparison. It is the unwritten `tests/test_rc.cpp` `[MT69]` gate and the external bound on `RcTailModel::TPeak`. Cite equations as "SLAC-PUB-380 Eq. (B.5)" unless someone has checked the RMP text (§1 #38).
- **Tsai**, *Radiative corrections to electron scatterings*, SLAC-PUB-848 (1971) — `held`, `refs/SLAC-PUB-0848.pdf` (INSPIRE recid 67278, free).
  The same formalism restated, and what most of the field actually cites for the exact formulas. Together with the entry above it closes `04_theory.md` T-36's action item without a purchase.
*(POLRAD 2.0 — Akushevich, Ilyichev, Shumeiko, Soroko, Tolkachev, CPC **104** (1997) 201, `refs/hep-ph_9706516.pdf` — is the radiative* formalism *the tree transcribes as much as it is a code; its entry is in §5c and is not repeated here.)*
- **Maximon, Tjon**, *Radiative corrections to electron–proton scattering*, Phys. Rev. C **62** (2000) 054320. [arXiv:nucl-th/0002058](https://arxiv.org/abs/nucl-th/0002058) — `held`, `refs/nucl-th_0002058.pdf`.
  How much Mo–Tsai's soft-photon/peaking treatment leaves out — i.e. a published bound on the error `RcTailModel::TPeak` carries.
- **Gakh, Shekhovtsova**, *Radiative corrections to deep-inelastic ed⁻ scattering. Case of tensor polarized deuteron*, JETP **99** (2004) 898. [arXiv:hep-ph/0403262](https://arxiv.org/abs/hep-ph/0403262) — `held`, `refs/hep-ph_0403262.pdf`.
  The source of `RC_DELTA_LOW_X` = 0.30 and of the *shape* the tree examined and rejected — one of only two tensor-target radiative calculations that exist.
- **Gakh, Shekhovtsova**, *Radiative events in DIS of unpolarized electrons by a tensor-polarized deuteron*, eConf **C030626** (2003) 351. [arXiv:hep-ph/0309123](https://arxiv.org/abs/hep-ph/0309123) — `downloaded 2026-09-16`, `refs/hep-ph_0309123.pdf` (3 pp.).
  The conference precursor of the entry above; useful only to date the calculation and to see which pieces were already present in 2003.
- **Gakh, Konchatnij, Merenkov**, *Model-independent QED corrections to elastic electron–deuteron scattering with a tensor-polarized target*, JETP **115** (2012) 212. [arXiv:1202.2225](https://arxiv.org/abs/1202.2225) — `held`, `refs/1202.2225.pdf`.
  **The only independent tensor-target radiative calculation besides Gakh–Shekhovtsova** — the one thing available to check POLRAD Eq. (A.4) against.
- **Afanasev, Akushevich, Gakh, Merenkov**, *Radiative corrections in coincidence (e,e′X) processes*, JETP **93** (2001) 449. [arXiv:hep-ph/0105032](https://arxiv.org/abs/hep-ph/0105032) — `held`, `refs/hep-ph_0105032.pdf`.
  Radiative corrections in **coincidence**, polarized and model-independent — the calculation against which the tagged channels' `rc_tail ≡ 1` (the one place the generator asserts a correction is exactly 1.0) should be priced.
- **Afanasev, Akushevich, Merenkov**, *Model-independent radiative corrections in processes of polarized electron–nucleon elastic scattering* (MASCARAD), Phys. Rev. D **64** (2001) 113009. [arXiv:hep-ph/0102086](https://arxiv.org/abs/hep-ph/0102086) — `held`, `refs/hep-ph_0102086.pdf`.
  Polarized **elastic** radiative corrections with realistic acceptance — the acceptance question buried inside POLRAD Eq. (44).
- **Liu, Melnitchouk, Qiu, Sato**, *Factorized approach to radiative corrections for inelastic lepton–hadron collisions*, Phys. Rev. D **104** (2021) 094033. [arXiv:2008.02895](https://arxiv.org/abs/2008.02895) — `held`, `refs/2008.02895.pdf`.
  Whether a fixed-target O(α) tail ports to collider kinematics at all — the question that decides whether the tree's POLRAD transcription is the right object for the EIC.
- **Cammarota, Qiu, Watanabe, Zhang**, *Factorized QED and QCD contribution to deeply inelastic scattering*, Phys. Rev. D **112** (2025) 056007. [arXiv:2505.23487](https://arxiv.org/abs/2505.23487) · doi:[10.1103/2d8y-ljwx](https://doi.org/10.1103/2d8y-ljwx) — `held`, `refs/2505.23487.pdf`.
  The successor framework to the entry above — NLO joint QED ⊗ QCD factorization. *(Reconciliation fix 2026-09-16: the first draft of this line carried the description "Joint QED and QCD factorization at NLO" in the title slot; the title above is the one on the PDF's title page, on Crossref and on INSPIRE, and the journal reference was added from INSPIRE/Crossref.)*
- **Afanasev, Bernauer, Blunden, Blümlein** *et al.*, *Radiative corrections: from medium to high energy experiments* (topical review), Eur. Phys. J. A **60** (2024) 91. [arXiv:2306.14578](https://arxiv.org/abs/2306.14578) · doi:[10.1140/epja/s10050-024-01281-y](https://doi.org/10.1140/epja/s10050-024-01281-y) — `held`, `refs/2306.14578.pdf`. *(Journal reference added in the reconciliation pass from Crossref; the draft said only "Eur. Phys. J. A (2023)".)*
  The modern baseline for "what 1.5 % on A_zz means", and the survey of which radiative codes are still alive.
- **Afanasev** *et al.*, *Radiative corrections for the Electron-Ion Collider* (CFNS ad-hoc workshop whitepaper, 2020). [arXiv:2012.09970](https://arxiv.org/abs/2012.09970) — `held`, `refs/2012.09970.pdf` (21.7 MB).
  The EIC-specific radiative problem statement, including nuclear targets — the document that says what the community thinks is unsolved.
- **Stein** *et al.*, *Electron scattering at 4° with energies of 4.5–20 GeV*, SLAC-PUB-1528 (1975); published as Phys. Rev. D **12** (1975) 1884 — `held`, `refs/SLAC-PUB-1528.pdf` (the free preprint).
  Reference [12] of the HERMES b₁ paper: the inclusive elastic + quasi-elastic radiative-tail recipe HERMES's correction is built on.

> **Also theory inputs, PLEASE DOWNLOAD** (full entries in §1): **#2** Akushevich–Shumeiko 1994 (POLRAD's own theory paper *for nuclear targets*, the only published place the A-dependence of a polarized radiative correction is worked out — the source `RcOptions::a_transfer_frac` does not have), **#16** Bardin–Shumeiko 1977 (the covariant infrared separation POLRAD Eq. (18) and `polrad_tpeak_quadrature` are built on; the free JINR-P2-10113 preprint is on disk but Russian-language), **#33** Kuraev–Fadin 1985 (the electron structure-function method and the exponentiated soft-photon form behind `ll_radiator` — **with the convention warning of §1**), **#34** de Forest–Walecka 1966 (`pauli_suppression`'s S(q)), **#15** Ent *et al.* 2001 (coincidence-RC practice, the price of `rc_tail ≡ 1`), **#38**/**#39** Mo–Tsai and Tsai 1961 for the published texts.

---

## 7. T5 / T6 — the chain and the software

The event record LiPolGen writes, the libraries it links, and the ePIC
simulation chain its output has to survive. Sourced from
`01_generators-chain.md` §1c–1e and the corpus inventory.

### 7a. The event record and the analysis chain

- **Buckley, Ilten, Konstantinov, Lönnblad, Monk, Pokorski, Przedzinski, Verbytskyi**, *The HepMC3 event record library for Monte Carlo event generators*, Comput. Phys. Commun. **260** (2021) 107310. [arXiv:1912.08005](https://arxiv.org/abs/1912.08005) · doi:[10.1016/j.cpc.2020.107310](https://doi.org/10.1016/j.cpc.2020.107310) — `held`, `refs/1912.08005.pdf`.
  The library LiPolGen writes with (3.3.0 in `deps/install`) and the format the whole ePIC chain consumes. It specifies the **arbitrary-attribute mechanism** `HEPMC3_CONVENTION.md` uses to carry ion spin, and the Asciiv3/ROOT/protobuf back-ends gate `F-3`'s round trips exercise. Gates `F-1`…`F-4`.
- **Gardiner, Isaacson, Pickering**, *NuHepMC: a standardized event record format for neutrino event generators*, SciPost Phys. Codebases **57** (2025). [arXiv:2310.13211](https://arxiv.org/abs/2310.13211) · doi:[10.21468/SciPostPhysCodeb.57](https://doi.org/10.21468/SciPostPhysCodeb.57) — `held`, `refs/2310.13211.pdf`.
  The **template** for a community attribute convention layered on HepMC3 — the published precedent for the ion-spin attribute convention (`F-4`, `CH-18`). LiPolGen is doing for polarized ions what NuHepMC did for neutrinos; this is the document that shows how such a convention is written down.
- **Bierlich** *et al.*, *Rivet version 4 release note*, SciPost Phys. Codebases **36** (2024). [arXiv:2404.15984](https://arxiv.org/abs/2404.15984) · doi:[10.21468/SciPostPhysCodeb.36](https://doi.org/10.21468/SciPostPhysCodeb.36) — `held`, `refs/2404.15984.pdf`.
  Rivet 4.1.2 ships in `eic_xl-nightly.sif`; it is the route to the 46 H1/ZEUS DIS analyses behind `CH-1`…`CH-5` — the cheapest external validation surface the chain has.
- **Buckley, Ferrando, Lloyd, Nordström, Page, Rüfenacht, Schönherr, Watt**, *LHAPDF6: parton density access in the LHC precision era*, Eur. Phys. J. C **75** (2015) 132. [arXiv:1412.7420](https://arxiv.org/abs/1412.7420) · doi:[10.1140/epjc/s10052-015-3318-8](https://doi.org/10.1140/epjc/s10052-015-3318-8) — `held`, `refs/1412.7420.pdf`.
  `[LHAPDF6]` — the grid library behind `src/lhapdf/lhapdf_sf.cpp`, i.e. the whole optional LHAPDF tier and every PDF-set name quoted in §2 and §4c.
- **Frank, Gaede, Grefe, Mato**, *DD4hep: a detector description toolkit for high energy physics experiments*, J. Phys. Conf. Ser. **513** (2014) 022010 · doi:[10.1088/1742-6596/513/2/022010](https://doi.org/10.1088/1742-6596/513/2/022010) — `free, not fetched` (IOP open access).
  The geometry layer under `npsim` / `eic/epic` — the first thing downstream of the HepMC3 file (`F-1`, `CH-13`).
- **Agostinelli** *et al.* (GEANT4 Collaboration), *GEANT4 — a simulation toolkit*, Nucl. Instrum. Meth. A **506** (2003) 250 · doi:[10.1016/S0168-9002(03)01368-8](https://doi.org/10.1016/S0168-9002(03)01368-8) — `free, not fetched` (Elsevier open access).
  The transport engine `npsim` runs (`F-1`). Both this and DD4hep are free and one click away; neither was committed because neither is read for a number.

> **Code with no paper (§7):** `01_generators-chain.md` §3.2 records that **eic-smear** (`github.com/eic/eic-smear`, GPL-3.0 — the generator-text ↔ ROOT ↔ HepMC3 bridge every EIC generator crosses, `CH-23`) and the **afterburner** / `abconv` (`github.com/eic/afterburner`, no declared license — crossing angle and beam effects, gate `F-2`, and the `CH-12` finding that it has no light-ion preset) have no literature record. Cite the repositories.

### 7b. The EIC machine, the detector and the far-forward region

- **Abdul Khalek** *et al.*, *Science requirements and detector concepts for the Electron-Ion Collider: EIC Yellow Report*, Nucl. Phys. A **1026** (2022) 122447. [arXiv:2103.05419](https://arxiv.org/abs/2103.05419) · doi:[10.1016/j.nuclphysa.2022.122447](https://doi.org/10.1016/j.nuclphysa.2022.122447) — `held`, `refs/2103.05419_part1.pdf` … `_part4.pdf` (corpus entry, split on size).
  `[YR]` — §8.1.6's generator verification, the TOPEG/Sartre/lAger provenance, and the far-forward acceptance tables. One of the two documents (with Chang *et al.*) the tree's far-forward efficiencies actually trace to.
- **Chang, Aschenauer, Jentsch, Kumar, Tu, Yin**, *Opportunities for imaging light nuclei with a second interaction region at the Electron-Ion Collider*, Phys. Rev. D **113** (2026) 032018. [arXiv:2511.05638](https://arxiv.org/abs/2511.05638) · doi:[10.1103/y4yv-y9dn](https://doi.org/10.1103/y4yv-y9dn) — `held`, `refs/2511.05638.pdf` (corpus entry).
  `[Chang26]` — the source of `COHERENT_JPSI_EFF_IR8_LI7`, of `Chang26EffEnergyRow` / `Chang26SpeciesEffRow`, and of the 0.1 < Q² < 100 GeV² acceptance window `estarlight_li6_q2_floors()` exists to step outside of. The whole O5 detection chain and row `E-16` rest on this one paper.
- **Gamage, Aschenauer, Berg, Burkert, Ent, Furletova, Higinbotham, Hutton, Morozov, Weiss** *et al.*, *Design concept for the second interaction region for the Electron-Ion Collider*, Proc. IPAC'21, TUPAB040, p. 1435. [arXiv:2105.13564](https://arxiv.org/abs/2105.13564) · doi:[10.18429/JACoW-IPAC2021-TUPAB040](https://doi.org/10.18429/JACoW-IPAC2021-TUPAB040) — `held`, `refs/2105.13564.pdf` (corpus entry).
  The IR-8 optics behind the second-IR acceptance the entry above assumes. Cited in the sibling repository's dispersive-tagging discussion rather than in the LiPolGen tree.
- **Jentsch** (for the ePIC Collaboration), *Far-forward detectors and physics with ePIC @ the EIC*, talk, DIS 2023, Michigan State University, 27–31 March 2023 — `held`, `refs/ePIC_far_forward_talk_DIS_2023_v2.pdf` (first page verified in this pass).
  The far-forward routing table `PHYSICS_CHANNELS.md` uses for α/d/t spectators — B0, off-momentum detectors, Roman pots, ZDC. A talk, not a paper: cite it as such, and prefer the Yellow Report or Chang *et al.* for anything quantitative.
- **Atoian, Buttimore, Ciullo, Cloët** *et al.*, *Realizing the scientific program with polarized ion beams at the future BNL Electron-Ion Collider* (EPIOS white paper), Phys. Rev. C **113** (2026) 060501. [arXiv:2510.10794](https://arxiv.org/abs/2510.10794) · doi:[10.1103/261w-8f38](https://doi.org/10.1103/261w-8f38) — `held`, `refs/2510.10794.pdf` (corpus entry).
  `[EPIOS]` — the community statement of which polarized ion species the EIC intends to deliver, and therefore the document that decides whether a polarized ⁶Li/⁷Li beam is a premise or a proposal. Cited in `PHYSICS_CHANNELS.md` and `benchmarking/00,06`.
- **Adam** *et al.* (ATHENA Collaboration), *ATHENA detector proposal — a totally hermetic electron nucleus apparatus proposed for IP6 at the Electron-Ion Collider*, JINST **17** (2022) P10019. [arXiv:2210.09048](https://arxiv.org/abs/2210.09048) · doi:[10.1088/1748-0221/17/10/P10019](https://doi.org/10.1088/1748-0221/17/10/P10019) — `held`, `refs/2210.09048.pdf` (corpus entry).
  A superseded detector proposal, kept because `PHYSICS_CHANNELS.md` cites it for acceptance and resolution figures that predate the ePIC documents. Check any number taken from it against the ePIC pTDR (§8) before it is quoted.
- **Chekanov** *et al.* (ZEUS), *Deep inelastic scattering with leading protons or large rapidity gaps at HERA*, Nucl. Phys. B **816** (2009) 1. [arXiv:0812.2003](https://arxiv.org/abs/0812.2003) · doi:[10.1016/j.nuclphysb.2009.03.003](https://doi.org/10.1016/j.nuclphysb.2009.03.003) — `held`, `refs/0812.2003v3.pdf` (corpus entry).
  Diffractive DIS measured both with a large rapidity gap and with the leading proton directly detected, Q² > 2 GeV², 40 < W < 240 GeV. `benchmarking/00_in_tree_checks.md` uses it for the diffractive azimuthal-asymmetry gate on the Pomeron tier.
- **Nikolaev, Pronyaev, Zakharov**, *Azimuthal asymmetry as a new handle on σ_L/σ_T in diffractive DIS*, Phys. Rev. D **59** (1999) 091501. [arXiv:hep-ph/9812212](https://arxiv.org/abs/hep-ph/9812212) · doi:[10.1103/PhysRevD.59.091501](https://doi.org/10.1103/PhysRevD.59.091501) — `held`, `refs/9812212v1.pdf` (corpus entry).
  The proposal to extract R^D from the azimuthal dependence of the diffractive cross section, resting on the model-independence of the LT-interference to transverse ratio. `06_critic.md` uses it as the theory statement behind the same diffractive azimuthal gate.

> **Shared-corpus entries the LiPolGen tree does not cite.** The corpus is shared with the sibling reconstruction/detector-R&D work, and about twenty of its rows are not referenced anywhere in this tree: the DIS kinematic-reconstruction family (Jacquet–Blondel / Bassler–Bernardi `hep-ex/9412004`, Arratia *et al.* `2110.05505`, Diefenthaler *et al.* `2108.11638`, Aggarwal–Caldwell `2206.04897`, Pecar–Vossen `2209.14489`, the ePIC inclusive-WG wiki resources and the Maple 2024 seminar slides), the near-beam superconducting-nanowire studies (`1907.13059`, `2312.13405`, `2410.00251`, `2510.11725`, `2601.03158`), the beam-polarization-transmission work (`2509.18558`), the tracking reference design (`li2023`), and two diffractive-DIS proceedings (`hep-ph/9808432`, `hep-ex/0206031`). They are on disk; their identities are as recorded in `refs_dict.json` and were **not** re-verified in this pass, and nothing in LiPolGen depends on them. Listed so they are neither re-searched nor mistaken for gaps.

---

## 8. Documents, standards and data archives

Not papers in the ordinary sense: programme documents, evaluated tables and
machine-readable archives. The three EIC/ePIC design documents and the PDG
review are open access and deliberately uncommitted on size grounds — the
numbers the tree actually uses come from §7b.

*(The two JLab PAC proposals — PR12-13-011 and PR12-15-005 — are programme documents too, but their entries are in §3a, beside the measurement they belong to, and are not repeated here.)*

- **Willeke** *et al.* (J. Beebe-Wang, ed.), *Electron-Ion Collider Conceptual Design Report*, BNL-221006-2021-FORE (2021), 972 pp. · doi:[10.2172/1765663](https://doi.org/10.2172/1765663) — `free, not fetched` (178 MB; corpus entry with `file: null`).
  Machine parameters and IR layout. The ion cap and ξ_p limit the sibling repository sources from it are the only numbers anything here takes from the CDR.
- **ePIC Collaboration**, *The ePIC Detector Preliminary Design Report* (September 2024 draft, chapters 2 and 8), Zenodo record 13866213 · doi:[10.5281/zenodo.13866213](https://doi.org/10.5281/zenodo.13866213) — `free, not fetched` (114 MB, 282 pp.; corpus entry with `file: null`).
- **ePIC Collaboration**, *The ePIC Detector Preliminary Technical Design Report, Version 3.1* (2026), Zenodo concept doi:[10.5281/zenodo.18271601](https://doi.org/10.5281/zenodo.18271601) (record 18271602, 2026-01-16, 553 pp., 260 MB; record 19496158 is the line-numbered build), CC-BY-4.0 — `free, not fetched` (corpus entry with `file: null`).
  The current detector document, and the right place to check any acceptance number the tree still takes from the ATHENA proposal.
- **Particle Data Group** (Navas *et al.*), *Review of Particle Physics*, Phys. Rev. D **110** (2024) 030001 · doi:[10.1103/PhysRevD.110.030001](https://doi.org/10.1103/PhysRevD.110.030001) — `free, not fetched`.
  Cited for one section only: *Monte Carlo particle numbering scheme*, which defines the `10LZZZAAAI` nuclear codes the tree emits for ⁶Li, ⁷Li, α, d and t (`CH-19`, `05_chain.md` §5.1). The normative source for the one convention a downstream reader cannot guess.
- **Nuclear Charge Density Archive** (Day *et al.*, University of Virginia) — `free, not fetched`. <https://discovery.phys.virginia.edu/research/groups/ncd/index.html>
  `FB_data.dat`, `SOG_data.dat`, `tabletr.txt` and C++ evaluators, collecting the fit parameters of ADNDT **14** (1974), **36** (1987) and **60** (1995). `01_li-nuclear-data.md` §1.3 read the ⁶Li Fourier–Bessel row (seven non-zero coefficients, cut-off R = 6.0 fm), the **⁷Li harmonic-oscillator row** `7Li 3 7 1 1.77 0.327 0 2.39 Su67` — the `HoSpin1FF::for_ion(LI7)` parameters, from data rather than from a paywalled paper — and the ⁴He sum-of-Gaussians row (rms 1.676 fm, twelve (R_i, Q_i) pairs, ΣQ_i = 1.000009). **There is no ⁶Li row in `tabletr.txt`**, and the archive carries **no magnetization rows at all**, which is why §1 #9 and #8 stay on the purchase list. The archive also mirrors scans of the three paywalled ADNDT volumes; those were deliberately not taken.
- **Quasi-elastic Electron–Nucleus Scattering Archive** (Benhar, Day, Sick), ⁶Li dataset `6Li.dat` — `free, not fetched`. <https://discovery.phys.virginia.edu/research/groups/qes-archive/data/6Li.dat>
  133 points in three settings, all from Heimlich *et al.* 1974 (§1 #12), re-fetched and re-counted in `01_li-nuclear-data.md` §6.2. Free ASCII; it is what turns `RC_QE_KF_GEV` from an asserted input into a fit residual. Cite the archive note (§4d) **and** the primary source.
- **HEPData** records `ins394050` (NMC F₂(⁶Li)/F₂(D)) and `ins393377` (NMC nuclear-ratio re-evaluation) — `free, not fetched`, CC0. <https://www.hepdata.net/record/ins394050>, <https://www.hepdata.net/record/ins393377>
  The machine-readable form of the only ⁶Li DIS measurements (§4c) — the actual numbers for a `D-5` comparison.
- **Stone**, *Table of nuclear electric quadrupole moments*, IAEA INDC(NDS)-0833 (2021) · doi:[10.61092/iaea.a6te-dg7q](https://doi.org/10.61092/iaea.a6te-dg7q); and *Table of recommended nuclear magnetic dipole moments: Part I, long-lived states*, IAEA INDC(NDS)-0794 (2019) — `free, not fetched` (both corpus entries with `file: null`).
  The free IAEA twins of the paywalled ADNDT table (§1 #42). `00_corpus.md` G10 notes that these are normally downloadable from `nds.iaea.org` and that nobody has attempted the fetch — the cheapest of all the outstanding items.
- **Tsai**, *Pair production and bremsstrahlung of charged leptons*, Rev. Mod. Phys. **46** (1974) 815 = SLAC-PUB-1365 — `free, not fetched`. Free at <https://inspirehep.net/files/be2ffc45ed59b1b28208cc386c1e16b3> (6.7 MB).
  The equivalent-radiator whose single-z form `ll_radiator` implements. Uncommitted because the external-bremsstrahlung half of it does not exist at a collider; take only the radiator.
- **Long, Higinbotham, Solvignon** (eds.), *Proceedings, Tensor Polarized Solid Target Workshop*, JLab, 10–12 March 2014, J. Phys. Conf. Ser. **543** (2014) (INSPIRE recid 1324998) — `free, not fetched` as a volume; the four individual papers this project uses have their own entries in §3a and §6f.
  Listed because it is the volume the tensor-target programme's written record lives in, and a reader who wants the context rather than four separate citations should go to it.

---

## 9. Unverified — named in the tree, not checked in this pass

Everything above was verified. These were **not**. Each is named somewhere in
the LiPolGen tree — in `PHYSICS_CHANNELS.md` §14's `[KEY]` list, or in
`benchmarking/04_theory.md`'s `T-1`…`T-46` survey, or in
`open_items/physics_literature.md` — and each was checked against the corpus
and against `PolarizedLithiumSim/refs/` in this pass and is **not on disk**;
none had its identity, its journal reference or its content confirmed against
a primary index here. The citation strings below are transcribed from the
tree's own text and should be treated as leads, not as citations, until
someone fetches the record. They are listed so the next pass starts here.

- **Wandzura, Wilczek**, *Sum rules for spin-dependent electroproduction — test of relativistic constituent quarks*, Phys. Lett. B **72** (1977) 195 — `unverified`. `[WW]`. The relation `g2_ww` implements; `00_corpus.md` §3a records that **no identifier for it exists anywhere in the repository**. Pre-arXiv; likely a free INSPIRE scan. The highest-value row in this section: the tree implements the relation and cites nothing.
- **Bissey**, Phys. Rev. C **65** (2002) 064317 — `unverified`. `[Bissey02]`, ³He per-nucleon polarizations.
- **Bissey, Guzey, Strikman** *et al.*, ³He g₁ convolution — `unverified`. `04_theory.md` T-41: the nearest thing to a polarized nuclear PDF for a light nucleus. No journal reference recorded in the tree.
- **Schellingerhout**, Phys. Rev. C **48** (1993) 2714 — `unverified`. `[Schell93]`, a ⁶Li cluster product named beside a scenario value.
- **Hulthén**, Ark. Mat. Astron. Fys. **28B** No. 5 (1942) — `unverified`. `[Hulthen42]`, the analytic radial form `Hulthen` is named after, used by name only.
- **Glauber**, *High-energy collision theory*, in *Lectures in Theoretical Physics* Vol. 1 (1959) — `unverified`. `[Glauber59]`, the eikonal formalism, used by name. A book chapter; there may be no citable digital copy.
- **Varshalovich, Moskalev, Khersonskii**, *Quantum Theory of Angular Momentum* — `unverified`. `[Varsh]`, the Wigner-d and Clebsch–Gordan conventions. `PHYSICS_CHANNELS.md` states plainly that no bibliographic identifier for it exists in the repository.
- **Ciofi degli Atti, Kaptari**, Phys. Rev. C **83** (2011) 044602. [arXiv:1011.5960](https://arxiv.org/abs/1011.5960) — `unverified`. `[CK11]`, the cluster-spectator Glauber survival product.
- **Ciofi degli Atti, Kaptari**. [arXiv:nucl-th/0407024](https://arxiv.org/abs/nucl-th/0407024) — `unverified`. `[CK04]`, generalized-eikonal vs plain-Glauber comparison.
- **Strikman, Weiss**, Phys. Rev. C **97** (2018) 035209. [arXiv:1706.02244](https://arxiv.org/abs/1706.02244) — `unverified`. `[SW18]`, the S[IA] + S[FSI] + S[FSI²] decomposition and its unitarity sum rule — the structure `fsi.hpp` is organised around.
- **Blinov** *et al.*, Phys. At. Nucl. **64** (2001) 907. [arXiv:nucl-ex/9910012](https://arxiv.org/abs/nucl-ex/9910012) — `unverified`. `[Blinov01]`, measured α–p σ_tot, σ_el and slope B — the *measured* inputs a ⁶Li Glauber weight would want.
- **Hirai, Kumano, Saito, Watanabe**, Phys. Rev. C **83** (2011) 035202. [arXiv:1008.1313](https://arxiv.org/abs/1008.1313) — `unverified`. `04_theory.md` T-5: the closest A > 2 unpolarized structure-function clustering analogue.
- **Zhao, Zhang, Liang, Liu, Zhou**. [arXiv:2206.11742](https://arxiv.org/abs/2206.11742), [arXiv:2401.10031](https://arxiv.org/abs/2401.10031) — `unverified`. T-14, covariant spin-3/2 polarization parametrization. **Already flagged UNVERIFIED in the tree itself**; this pass did not lift that.
- **Piarulli, Pastore, Wiringa** *et al.*, Phys. Rev. C **107** (2023) 014314. [arXiv:2210.02421](https://arxiv.org/abs/2210.02421) — `unverified`. T-22, newer chiral-EFT densities and momentum distributions for A ≤ 12.
- **Hebborn, Brune, Phillips** (2025). [arXiv:2510.19067](https://arxiv.org/abs/2510.19067) — `unverified`. T-23, NCSMC ab-initio ⁶Li correlations. The tree's own note says "read the paper first — the abstract does not confirm a D/S or Q result".
- **Brida, Pieper, Wiringa**, Phys. Rev. C **84** (2011) 024319. [arXiv:1106.3121](https://arxiv.org/abs/1106.3121) — `unverified`. T-24, QMC spectroscopic overlaps for A ≤ 7 (single-nucleon, not α–d).
- **Cosyn, Sargsian** (2011). [arXiv:1012.0293](https://arxiv.org/abs/1012.0293); **Cosyn, Melnitchouk, Sargsian** (2014). [arXiv:1311.3550](https://arxiv.org/abs/1311.3550) — `unverified`. T-26, generalized-eikonal FSI against JLab data.
- **Kaptari, Del Dotto, Pace, Salmè, Scopetta**, Phys. Rev. C **89** (2014) 035206. [arXiv:1307.2848](https://arxiv.org/abs/1307.2848) — `unverified`. T-28, the **polarized** cluster-spectator distorted spectral function at A = 3 — the nearest published object to what a polarized tagged ⁶Li calculation would need.
- **Sargsian, Strikman**, Phys. Lett. B **639** (2006) 223. [arXiv:hep-ph/0511054](https://arxiv.org/abs/hep-ph/0511054) — `unverified`. T-30, the no-loop / pole-extrapolation theorem.
- **Sargsian**, Int. J. Mod. Phys. E **10** (2001) 405. [arXiv:nucl-th/0110053](https://arxiv.org/abs/nucl-th/0110053) — `unverified`. T-31, generalized-eikonal foundations.
- **Akushevich, Ilyichev**. [arXiv:1403.3421](https://arxiv.org/abs/1403.3421) — `unverified`. T-38, NLO radiative corrections for exclusive photon electroproduction.
- **Magdy** *et al.*, Eur. Phys. J. A **60** (2024). [arXiv:2405.07844](https://arxiv.org/abs/2405.07844) — `unverified`. α-clustering work at the EIC; the tree's own note records it is BeAGLE-based and not diffraction.

> **Code, not literature:** `github.com/wcosyn/physics-code` (T-32, a second independent FSI implementation) and `github.com/hejajama/subnucleondiffraction` (T-45, the code behind Mäntysaari *et al.*) are repositories named in the survey. Neither was checked in this pass and neither has a paper of its own — cite the repository and the physics papers of §6f and §6g.

---

## 10. Counts, and how this file was checked

### 10.1 Counts

| status | count | where |
|---|---|---|
| `held` | **172** | §2–§9 — a PDF is on disk under `PolarizedLithiumSim/refs/`, named in the entry |
| `downloaded 2026-09-16` | **3** | §4c (Gomez E139, `SLAC-PUB-5813.pdf`), §5a (Helenius *et al.* 2026, `2605.00502.pdf`), §6h (Gakh–Shekhovtsova 2003, `hep-ph_0309123.pdf`) |
| `PLEASE DOWNLOAD` | **56** | §1, items 1–58 less the one resolved and the one found free below; cross-referenced from the tier sections |
| `free, not fetched` | **13** | §7a, §8 (12 bullets) and §1 item 27 (found free in the reconciliation pass, §10.4) — open access and reachable, no PDF committed |
| `unverified` | **22** | §9 |
| **total entries** | **266** | 209 status-bearing entries in §2–§9 plus §1's 56 to obtain plus §1 item 27 |

Two bookkeeping notes, so the arithmetic is not mistaken for precision it
does not have. **(i)** §1 numbers **58** items; item **43** (Gomez *et al.*,
E139) was **resolved in this pass** — a free preprint scan was found and
fetched — so it is a `downloaded 2026-09-16` entry in §4c and is *not* one of
the 56 references to obtain; item **27** (Albrecht *et al.*) likewise left the
purchase list in the reconciliation pass (§10.4), which located a free JINR
preprint on INSPIRE but was not allowed to commit it. **(ii)** One work is deliberately counted twice:
Mo & Tsai appears as a `held` entry in §6h (the free SLAC-PUB-380 preprint,
which is what the gates should cite) and as §1 item **38** (the Rev. Mod.
Phys. text, wanted only for published equation numbering). Netting that pair
gives **265 distinct works**.

These reproduce by grep over this file:

```sh
grep -cE '^- \*\*.*`held`'                        REFERENCES.md   # 172
grep -cE '^- \*\*.*`downloaded 2026-09-16`'       REFERENCES.md   #   3
grep -cE '^- \*\*.*`free, not fetched`'           REFERENCES.md   #  12
grep -cE '^- \*\*.*`unverified`'                  REFERENCES.md   #  22
grep -cE '^[0-9]+\. \*\*.*\*\*PLEASE DOWNLOAD\*\*' REFERENCES.md  #  56
grep -cE '^[0-9]+\. \*\*.*\*\*RESOLVED'            REFERENCES.md  #   1
grep -cE '^[0-9]+\. \*\*.*`free, not fetched`'    REFERENCES.md  #   1  (§1 item 27)
```

Every entry bullet carries **exactly one** status token, which is what makes
the counts exclusive.

### 10.2 What was checked in this final pass

1. **Every one of §1's 58 purchase items was re-queried against the INSPIRE
   literature API** (by title, and where a title search failed, by journal
   reference `j Journal,volume,page` or by the record id the entry already
   named). Result: **not one of the 58 has an arXiv e-print.** Four records
   were found in that pass to carry an INSPIRE-hosted fulltext, and each was
   opened and identified (the reconciliation pass, §10.4, found two more:
   **#27**, a genuine JINR preprint, and **#50**, a file that does not parse):
   - **#43 Gomez *et al.* (E139)** — a genuine free PDF (SLAC-PUB-5813, 86 pp.).
     The first pass was **wrong** to list it as purchase-only. Fetched, now on
     disk, entry rewritten, and a new key appended to the shared corpus index.
   - **#52 Azhgirey *et al.* 2004** — full text, but as **Elsevier XML**, not a
     PDF. Identity confirmed from the XML's own `<ce:title>`, `<prism:doi>` and
     author surnames. Not committed; the link is recorded in the entry.
   - **#38 Mo & Tsai** — the hosted file is SLAC-PUB-380, **already on disk**,
     which is exactly what §1 item 38 already said.
   - **#16 Bardin & Shumeiko** — the hosted file is the JINR-P2-10113 preprint,
     **already on disk** and Russian-language, exactly as §1 item 16 says.
   In addition, **#39** (Tsai 1961) had no title in the first pass; INSPIRE
   record 46623 supplies it (*Radiative corrections to electron–proton
   scattering*) and confirms no e-print and no fulltext.
2. **Every `held` and `downloaded` entry was checked against the filesystem.**
   The file names 177 distinct paths of the form `refs/<name>.pdf`; **all 177
   exist** under `/home/cpeng/Projects/polli/PolarizedLithiumSim/refs/`.
   Two files were read directly to confirm what they are, because neither has
   an arXiv record: `TUNL_A6_2002.pdf` (the TUNL *Energy Levels of Light
   Nuclei, A = 6* chapter) and `ePIC_far_forward_talk_DIS_2023_v2.pdf` (the
   Jentsch DIS 2023 talk).
3. **Three corpus citation strings were corrected** against INSPIRE, and the
   corrections are recorded in the entries rather than only here:
   `nucl-th/0605061` is *EMC and polarized EMC effects in nuclei* (PLB 642
   (2006) 210), **not** *Spin-dependent structure functions in nuclear matter
   and the polarized EMC effect*, which is a different Cloët–Bentz–Thomas
   paper; `2109.03591` is *Gluon EMC effects in nuclear matter* (J. Phys. G 49
   (2022) 03LT01), not *Polarized gluon EMC effect*; and §1 item 43's title is
   *…the A-dependence…*, hyphenated as INSPIRE and the preprint have it.
4. **Two changes were made outside this directory**, both of the kind the
   reference stage is allowed to make: the new file
   `PolarizedLithiumSim/refs/SLAC-PUB-5813.pdf`, and one appended entry
   (`gomez1993-e139-slac-pub-5813-preprint`) in
   `PolarizedLithiumSim/refs/refs_dict.json`, matching that file's existing
   schema and formatting. The append was confirmed append-only by `git diff`:
   the only deleted line anywhere in that file belongs to the concurrent
   workflow's `_updated` bump, not to this one.

### 10.3 What this file does not claim

No statement above about what a paper *contains* rests on anything but its
abstract, its INSPIRE/arXiv record, or the text of the PDF — either read in
the domain-file series this bibliography is assembled from, or read here.
Where a specific number or a figure/equation number is quoted, the entry says
which document it was read from. The **56** references of §1 still to obtain were not read at
all — their entries state only what the project needs them *for*, which is
established from the tree's own code and documents, plus whatever the title,
the abstract or a free secondary source supports; those are the claims to
re-check first when the PDFs arrive (item 27's free preprint was opened in the
reconciliation pass only far enough to identify it). §9's 22 rows are leads
transcribed from the tree, not citations.

### 10.4 Reconciliation pass (2026-09-16, second reviewer)

An adversarial re-check of this file against the corpus and the primary
indices, made after §10.1–10.3 were written. What it did, and the six
changes it made in place:

1. **Corpus.** All 194 `refs/*.pdf` files carry the `%PDF` magic and exceed
   20 kB; for each of the 93 file paths named by a `refs_dict.json` entry,
   `pdftotext` of the first pages contains the entry's title words (the only
   sub-threshold hits are the Yellow Report's split parts 1/2/4, whose title
   page is in part 3). No entry names a missing file. 101 PDFs have no entry;
   all 101 are untracked files fetched on 2026-09-16 by the concurrent
   reference workflow (the 43 PDFs tracked at `HEAD` are all indexed).
   `git diff refs/refs_dict.json refs/README.md` is append-only: `HEAD`'s 54
   entries are byte-identical and in order, 83 keys were appended (82 by the
   concurrent workflow, `gomez1993-e139-slac-pub-5813-preprint` by this
   series), the README gained 17 trailing lines, and the only changed
   pre-existing line is the top-level `_updated` date.
2. **Identity.** Every `held`/`downloaded` bullet's title was matched against
   the first pages of its own PDF (175 bullets, all ≥ 70 % of title words
   present). Every DOI in the file was resolved through Crossref (DataCite for
   the JACoW and Zenodo DOIs) and its volume/year compared with the citation;
   all agree. Every INSPIRE record id named in the file was fetched and is the
   paper claimed.
3. **Corrections made.** (a) §6h, Cammarota *et al.* 2505.23487: the title
   slot held a description, *Joint QED and QCD factorization at NLO*; the
   PDF's title page, Crossref and INSPIRE give *Factorized QED and QCD
   contribution to deeply inelastic scattering*, published Phys. Rev. D 112
   (2025) 056007 — corrected and the journal reference added. (b) §6h,
   Afanasev *et al.* 2306.14578: "Eur. Phys. J. A (2023)" replaced by the
   published Eur. Phys. J. A 60 (2024) 91 with its DOI. (c) §3c, Tkachenko
   *et al.*: the bullet had no title; added from the PDF. (d) §1 #50,
   Miettinen–Pumplin: the INSPIRE "fulltext" the first draft called a free
   scan is served as `application/pdf` but the bytes are not a PDF (fetched
   twice, no `%PDF` header, `pdfinfo` rejects it); the entry now says so and
   records the original Fermilab preprint URL (`lss.fnal.gov`, HTTP 403 to an
   automated fetch). Consequently §10.2 item 1's "four records carry an
   INSPIRE-hosted fulltext" is five: the fifth is #50 and it does not parse.
   (e) §1 #27, Albrecht *et al.* 1980: INSPIRE record 143652 carries the Dubna preprint JINR-E1-12727 (1979), fetched here and read — a genuine 10-page PDF of the paper — which the first draft missed exactly as it had missed #43; the item now carries `free, not fetched` (this pass could not add files to the corpus), and §1's preamble, §4b, §10.1 (PLEASE DOWNLOAD 57 → 56, free-not-fetched 12 → 13, total unchanged at 266) and §10.2 were updated to match. (f) This section.
4. **§1 links and arXiv.** Each of the 57 remaining §1 links was exercised: every DOI resolves (HTTP 302 to the publisher landing page; the two non-DOI items resolve on INSPIRE), and Crossref returns the cited title, volume and year for each. Each of the 58 titles was searched on the arXiv API (exact-phrase and keyword `ti:` queries) and each DOI was looked up on INSPIRE for an `arxiv_eprints` field: **no item has an arXiv copy** — every title hit is a different-author namesake, except hep-ph/9804361 under #2, which is the same authors' *1998 nucleon* paper (J. Phys. G 24 (1998) 1995), not the 1994 light-nuclei paper. INSPIRE has no record for #47, #48 and #56 (non-HEP journals), which is why those three rest on Crossref alone. The bare links in §1d and §8 were also exercised: INSPIRE's Tsai 1974 file is a real 131-page PDF (SLAC-PUB-1365 with erratum); the UVa and QES archives and both JLab proposals answer 200; HEPData, HAL, `link.aps.org` and `research.vu.nl` answer automated requests with a bot challenge, as the entries say; Zenodo answered 504 to every request during this pass, so its record ids rest on the DataCite lookups of the two concept DOIs, which return the ePIC PDR draft (2024) and pTDR v3.1 (2026).

- 2026-09-16 (corpus indexing pass): two titles in this file were corrected from the PDFs' own first pages — arXiv:2202.12200 (the earlier line carried the title of a different Guzey et al. paper; the paper is *Coherent J/ψ electroproduction on ⁴He and ³He at the Electron-Ion Collider: probing nuclear shadowing one nucleon at a time*) and arXiv:hep-ph/0403262 (the earlier line was a paraphrase; the paper's title is *Radiative corrections to deep-inelastic ed⁻ scattering. Case of tensor polarized deuteron*). Journal references unchanged.
