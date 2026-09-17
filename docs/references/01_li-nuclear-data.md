# 01 — External references: ⁶Li / ⁷Li nuclear structure and electron-scattering data

**Search date 2026-09-16.** Domain: the measured ⁶Li, ⁷Li and ⁴He nuclear
inputs LiPolGen's kernels are built on — elastic charge and magnetic form
factors, charge-density compilations, α–d cluster knockout, DIS on lithium,
quasi-elastic (e,e′), moments and radii, the ab-initio wave functions under
`data/vmc/`, and polarized ⁶Li target material. Companion to `00_corpus.md`
(what the project already holds) and to `../benchmarking/03_data_nuclear.md`
(which said what the tree *needs*; this file says where each item actually
is, whether it can be had for free, and what it pins).

**Method.** Every entry below was verified in this pass over the network:
Crossref (`api.crossref.org/works/{DOI}` and bibliographic search) for
title / authors / year / journal / volume / pages / DOI; the arXiv API
(`export.arxiv.org/api/query`) for preprint identity, abstract and
journal-ref; OpenAlex and Semantic Scholar for abstracts and open-access
status of the pre-arXiv literature; INSPIRE for record ids; direct HTTP
fetch for the machine-readable archives. Four resources were verified by
*reading their contents*, not just their metadata (§1.3, §5.2, §6.2, §7.3).
Nothing below is quoted from memory, and the three items that could not be
verified in full are marked **UNVERIFIED** with the reason.

**What is new in this pass**, in one list, because most of this domain is
1966–1999 journal literature the tree has only ever cited second-hand:

1. **The ⁷Li harmonic-oscillator parameters are machine-readable and free** —
   `tabletr.txt` of the UVa Nuclear Charge Density archive, row `7Li 3 7 1
   1.77 0.327 0 2.39 Su67`, in exactly the functional form `HoSpin1FF`
   ships (§1.3). Survey priority 6 costs a download, not a transcription.
2. **Suelzle's abstract says a harmonic-oscillator model does not fit ⁶Li** —
   verbatim, and it says the *modified* HO does not either (§2.1). That is a
   published statement about the functional form `HoSpin1FF` uses for ⁶Li,
   and it is stronger than `rc.hpp`'s "TO BE REFIT".
3. **The ⁶Li Fourier–Bessel charge density is free, machine-readable, and
   reproduces its own published radius and diffraction zero** — fetched and
   integrated here: ∫ρ dV = 2.991 (Z = 3), r_rms = 2.5206 fm, first zero of
   F_ch at **2.6944 fm⁻¹** (§1.3).
4. **The ⁴He sum-of-Gaussians is in the same archive** — 12 (R_i, Q_i) pairs,
   ΣQ_i = 1.000009, rms 1.676 fm — so survey priority 2 needs no transcription
   from the de Vries pages at all (§7.3).
5. **A citation error in the tree is corrected**: the "⁶Li(e,e′α) and (e,e′d)
   at 520 MeV, PLB **51** (1974)" of `03_data_nuclear.md` is Genin *et al.*,
   PLB **52** (1974) 46. PLB **51** (1974) 217 is a different paper —
   Toyama & Sakamoto, *theory* (§5.4).
6. **The QES-archive ⁶Li dataset is exactly as the survey described it** —
   re-fetched and re-counted here: 133 rows, 40/47/46 at three settings,
   single source `Heimlich:1973` (§6.2).
7. **No η(⁶Li → α+d) determination has been published since George–Knutson
   1999** — searched and found nothing but astrophysical S-wave ANC work
   (§10.4).

---

## 1. Master table

`OA` = open access. `disk` = a PDF is now committed under
`PolarizedLithiumSim/refs/` (filename given); "—" means no free copy exists
and the DOI link is the purchase route. `serves` names the LiPolGen symbol,
gate or open item. Tier is `BENCHMARK_PLAN.md`'s **T3** unless stated.

| # | key | citation | link | OA | disk | serves |
|---|---|---|---|---|---|---|
| **A. Elastic charge form factor and charge densities** |
| A1 | `suelzle1967-li6-li7-elastic` | L. R. Suelzle, M. R. Yearian, H. Crannell, *Elastic Electron Scattering from Li⁶ and Li⁷*, Phys. Rev. **162** (1967) 992–1005 | [doi:10.1103/PhysRev.162.992](https://doi.org/10.1103/PhysRev.162.992) | no | — | `LI6_FF_HO_A_FM`, `LI6_FF_HO_ALPHA`, `HoSpin1FF::for_ion(LI7)`, Q1 |
| A2 | `li-sick-whitney-yearian1971-6li-form-factor` | G. C. Li, I. Sick, R. R. Whitney, M. R. Yearian, *High-energy electron scattering from ⁶Li*, Nucl. Phys. **A162** (1971) 583–592 | [doi:10.1016/0375-9474(71)90257-0](https://doi.org/10.1016/0375-9474(71)90257-0) | no | — (in corpus, `file: null`) | T11 `q₀` window; `T_DIP` in the coherent channel |
| A3 | `devries1987-adndt36-charge-densities` | H. de Vries, C. W. de Jager, C. de Vries, *Nuclear charge-density-distribution parameters from elastic electron scattering*, At. Data Nucl. Data Tables **36** (1987) 495–536 | [doi:10.1016/0092-640X(87)90013-1](https://doi.org/10.1016/0092-640X(87)90013-1) | no (scan mirrored, §1.3) | — | `LI6_R2_POINT_FM2`; the ⁷Li HO row |
| A4 | `dejager1974-adndt14-charge-magnetization` | C. W. de Jager, H. de Vries, C. de Vries, *Nuclear charge- and magnetization-density-distribution parameters from elastic electron scattering*, At. Data Nucl. Data Tables **14** (1974) 479–508 | [doi:10.1016/S0092-640X(74)80002-1](https://doi.org/10.1016/S0092-640X(74)80002-1) | no (scan mirrored, §1.3) | the **magnetization** half — `F_m`, and the Rand–Frosch–Yearian numbers |
| A5 | `fricke1995-adndt60-charge-radii` | G. Fricke, C. Bernhardt, K. Heilig, L. A. Schaller, L. Schellenberg, E. B. Shera, C. W. de Jager, *Nuclear Ground State Charge Radii from Electromagnetic Interactions*, At. Data Nucl. Data Tables **60** (1995) 177–285 | [doi:10.1006/adnd.1995.1007](https://doi.org/10.1006/adnd.1995.1007) | green (repository) | — | cross-check on `LI6_R2_POINT_FM2` |
| A6 | `uva-ncd-archive` | Nuclear Charge Density Archive (D. Day *et al.*, U. Virginia) — `FB_data.dat`, `SOG_data.dat`, `tabletr.txt`, C++ evaluators, and scans of ADNDT 14/36/60 | [discovery.phys.virginia.edu/research/groups/ncd](https://discovery.phys.virginia.edu/research/groups/ncd/index.html) | **yes, free + machine-readable** | not committed (data, not a paper) | ⁶Li C0 from data; **⁷Li HO parameters**; ⁴He SOG |
| A7 | `angeli-marinova2013-charge-radii` | I. Angeli, K. P. Marinova, *Table of experimental nuclear ground state charge radii: An update*, At. Data Nucl. Data Tables **99** (2013) 69–95 | [doi:10.1016/j.adt.2011.12.006](https://doi.org/10.1016/j.adt.2011.12.006) | no | — | `LI6_R2_POINT_FM2` = 6.0788 fm² (r_ch = 2.589(39) fm) |
| A8 | `sick1974-model-independent-densities` | I. Sick, *Model-independent nuclear charge densities from elastic electron scattering*, Nucl. Phys. **A218** (1974) 509–541 | [doi:10.1016/0375-9474(74)90039-6](https://doi.org/10.1016/0375-9474(74)90039-6) | no | — | the method behind A6's FB rows (cut-off, incompleteness error) |
| A9 | `wiringa-schiavilla1998-6li-form-factors` | R. B. Wiringa, R. Schiavilla, *Microscopic calculation of ⁶Li elastic and transition form factors*, Phys. Rev. Lett. **81** (1998) 4317–4320 | [arXiv:nucl-th/9807037](https://arxiv.org/abs/nucl-th/9807037) | **yes** | `nucl-th_9807037.pdf` | the monopole **shape** `rc.hpp` borrows; Q(⁶Li) = −0.23(9) fm² caution |
| A10 | `bergstrom1979-6li-cluster-form-factors` | J. C. Bergstrom, *⁶Li electromagnetic form factors and phenomenological cluster models*, Nucl. Phys. **A327** (1979) 458–476 | [doi:10.1016/0375-9474(79)90269-0](https://doi.org/10.1016/0375-9474(79)90269-0) | no | — | the α–d cluster picture *fitted to the elastic data* — the closest published analogue of `HoSpin1FF` + `ClusterConfigSampler` |
| **B. Magnetization and magnetic form factors** |
| B1 | `rand-frosch-yearian1966-magnetization` | R. E. Rand, R. Frosch, M. R. Yearian, *Elastic Electron Scattering from the Magnetic Multipole Distributions of Li⁶, Li⁷, Be⁹, B¹⁰, B¹¹ and N¹⁴*, Phys. Rev. **144** (1966) 859–873 | [doi:10.1103/PhysRev.144.859](https://doi.org/10.1103/PhysRev.144.859) | no | — | `F_m(0)`/`LI6_MU_N` slope; **BENCHMARK_PLAN §4 item 6** |
| B2 | `bergstrom1982-6li-magnetic-form-factor` | J. C. Bergstrom, S. B. Kowalski, R. Neuhausen, *Elastic magnetic form factor of Li⁶*, Phys. Rev. C **25** (1982) 1156–1167 | [doi:10.1103/PhysRevC.25.1156](https://doi.org/10.1103/PhysRevC.25.1156) | no | — | the modern (higher-q) replacement for B1's ⁶Li half |
| B3 | `vanniftrik1971-7li-magnetization` | G. J. C. van Niftrik, L. Lapikás, H. de Vries, G. Box, *Magnetization distribution of the ⁷Li nucleus as obtained from electron scattering through 180°. The electric quadrupole moment of ⁷Li*, Nucl. Phys. **A174** (1971) 173–192 | [doi:10.1016/0375-9474(71)91011-6](https://doi.org/10.1016/0375-9474(71)91011-6) | no | — | `LI7_QUADRUPOLE_FM2` — the *electron-scattering* determination |
| B4 | `lichtenstadt1989-7li-doublet-form-factors` | J. Lichtenstadt, J. Alster, M. A. Moinester, J. Dubach, R. S. Hicks, G. A. Peterson, S. Kowalski, *High momentum transfer longitudinal and transverse form factors of the ⁷Li ground-state doublet*, Phys. Lett. B **219** (1989) 394–398 | [doi:10.1016/0370-2693(89)91083-6](https://doi.org/10.1016/0370-2693(89)91083-6) | no | — | ⁷Li F_L and F_T above Suelzle's q range |
| **C. α–d cluster knockout — the tagged channel's own observable** |
| C1 | `ent1994-eed-4he-6li-12c` | R. Ent, B. L. Berman, H. P. Blok, J. F. J. van den Brand, W. J. Briscoe, M. N. Harakeh, E. Jans, P. D. Kunz, L. Lapikás, *The (e,e′d) reaction on ⁴He, ⁶Li, and ¹²C*, Nucl. Phys. **A578** (1994) 93–133 | [doi:10.1016/0375-9474(94)90971-7](https://doi.org/10.1016/0375-9474(94)90971-7) | no | — | `Wave::radial`, `data/vmc/li6_alpha_d/li6.ad`, T1 spectator spectrum |
| C2 | `ent1986-alpha-d-momentum-distribution` | R. Ent, H. P. Blok, J. F. A. van Hienen, G. van der Steenhoven, J. F. J. van den Brand, J. W. A. den Herder, E. Jans, P. H. M. Keizer, L. Lapikás, E. N. M. Quint, P. K. A. de Witt Huberts, B. L. Berman, W. J. Briscoe, C. T. Christou, D. R. Lehman, B. E. Norum, A. Saha, *Reaction ⁶Li(e,e′d)⁴He and the α–d Momentum Distribution in the Ground State of ⁶Li*, Phys. Rev. Lett. **57** (1986) 2367–2370 | [doi:10.1103/PhysRevLett.57.2367](https://doi.org/10.1103/PhysRevLett.57.2367) | no | — | same; the letter, and the one whose **title** is the observable |
| C3 | `mitchell1991-6li-eealpha` | J. H. Mitchell, H. P. Blok, B. L. Berman, W. J. Briscoe, M. A. Daman, R. Ent, E. Jans, L. Lapikás, *Mechanism of the ⁶Li(e,e′α) reaction*, Phys. Rev. C **44** (1991) 2002–2005 (INSPIRE 331689) | [doi:10.1103/PhysRevC.44.2002](https://doi.org/10.1103/PhysRevC.44.2002) | no | — | the α-knockout mirror of C1/C2 |
| C4 | `genin1974-6li-eealpha-eed-520mev` | J. P. Genin, J. Julien, M. Rambaut, C. Samour, A. Palmeri, Vinciguerra, *The ⁶Li(e,e′α) and ⁶Li(e,e′d) reactions at 520 MeV*, Phys. Lett. B **52** (1974) 46–48 | [doi:10.1016/0370-2693(74)90714-X](https://doi.org/10.1016/0370-2693(74)90714-X) | no | — | first-generation C1/C3; **corrects a volume number in the tree** |
| C5 | `albrecht1980-6li-ppd-670mev` | D. Albrecht, M. Csatlós, J. Erö, Z. Fodor, I. Hernyes, Mu Hongsung, B. A. Khomenko, N. N. Khovanskij *et al.*, *Large-angle quasi-free scattering in ⁶Li(p,pd)⁴He at 670 MeV*, Nucl. Phys. **A338** (1980) 477–494 | [doi:10.1016/0375-9474(80)90045-7](https://doi.org/10.1016/0375-9474(80)90045-7) | no | — | hadronic probe of the same α–d distribution; the FWHM 70 vs 120 MeV/c band |
| C6 | `dollhopf1975-6li-alpha2alpha` | W. Dollhopf, C. F. Perdrisat, P. Kitching, W. C. Olsen, *The ⁶Li(α,2α)²H reaction at 700 MeV and the α–d momentum distribution*, Phys. Lett. B **58** (1975) 425–427 | [doi:10.1016/0370-2693(75)90579-1](https://doi.org/10.1016/0370-2693(75)90579-1) | no | — | third, α-induced probe of the same distribution |
| C7 | `connelly1998-trinucleon-knockout` | J. P. Connelly, B. L. Berman, W. J. Briscoe, K. S. Dhuga, A. Mokhtari, D. Zubanov, H. P. Blok, R. Ent, J. H. Mitchell, L. Lapikás, *Trinucleon cluster knockout from ⁶Li*, Phys. Rev. C **57** (1998) 1569–1573 | [doi:10.1103/PhysRevC.57.1569](https://doi.org/10.1103/PhysRevC.57.1569) · OA copies §5.5 | **bronze/green** | — (bot-walled, §11) | the reaction-mechanism **systematic** on C1–C6 |
| C8 | `lapikas1999-7li-eep-vmc` | L. Lapikás, J. Wesseling, R. B. Wiringa, *Nuclear Structure Studies with the ⁷Li(e,e′p) Reaction*, Phys. Rev. Lett. **82** (1999) 4404–4407 | [arXiv:nucl-th/9904008](https://arxiv.org/abs/nucl-th/9904008) | **yes** | `nucl-th_9904008.pdf` | **E-13** — the only data-vs-VMC test of the `data/vmc/` family |
| C9 | `hotta1999-6li-eep-low-q` | T. Hotta, T. Tamae, T. Miura, H. Miyase, I. Nakagawa, T. Suda, M. Sugawara, T. Tadokoro, *Measurement of the ⁶Li(e,e′p) reaction cross sections at low momentum transfer*, Nucl. Phys. **A645** (1999) 492–508 | [arXiv:nucl-ex/9810016](https://arxiv.org/abs/nucl-ex/9810016) | **yes** | `nucl-ex_9810016.pdf` | the low-ω end of the RC tail |
| **D. EMC / DIS on lithium** |
| D1 | `nmc1995-f2-li-d` | M. Arneodo *et al.* (New Muon Collaboration), *The structure function ratios F2(Li)/F2(D) and F2(C)/F2(D) at small x*, Nucl. Phys. **B441** (1995) 12–30 | [arXiv:hep-ex/9504002](https://arxiv.org/abs/hep-ex/9504002) · HEPData [ins394050](https://www.hepdata.net/record/ins394050) | **yes (CC0 data)** | `hep-ex_9504002.pdf` | **`EMC_VALENCE_DEPLETION_EPPS21`** — the only ⁶Li DIS measurement |
| D2 | `nmc1995-a-dependence-li-c-ca` | P. Amaudruz *et al.* (New Muon Collaboration), *A re-evaluation of the nuclear structure function ratios for D, He, ⁶Li, C and Ca*, Nucl. Phys. **B441** (1995) 3–11 | [arXiv:hep-ph/9503291](https://arxiv.org/abs/hep-ph/9503291) · HEPData [ins393377](https://www.hepdata.net/record/ins393377) | **yes** | `hep-ph_9503291.pdf` | F₂(C)/F₂(⁶Li), F₂(Ca)/F₂(⁶Li) — same ⁶Li target, lower Q² |
| D3 | `gomez1994-e139-a-dependence` | J. Gomez, R. G. Arnold, P. E. Bosted, C. C. Chang, A. T. Katramatou, G. G. Petratos, A. A. Rahbar, S. E. Rock *et al.*, *Measurement of the A dependence of deep-inelastic electron scattering*, Phys. Rev. D **49** (1994) 4348–4372 | [doi:10.1103/PhysRevD.49.4348](https://doi.org/10.1103/PhysRevD.49.4348) | no | — | the **negative**: E139 has no lithium |
| D4 | `seely2009-emc-light-nuclei` | J. Seely, A. Daniel, D. Gaskell, J. Arrington, N. Fomin, P. Solvignon *et al.*, *New measurements of the EMC effect in very light nuclei*, Phys. Rev. Lett. **103** (2009) 202301 | [arXiv:0904.4448](https://arxiv.org/abs/0904.4448) | **yes** | `0904.4448.pdf` | the ⁹Be cluster/EMC argument — the nearest measured statement that clustering matters |
| D5 | `arrington2021-emc-light-heavy` | J. Arrington, J. Bane, A. Daniel, N. Fomin, D. Gaskell, J. Seely *et al.*, *Measurement of the EMC effect in light and heavy nuclei*, Phys. Rev. C **104** (2021) 065203 | [arXiv:2110.08399](https://arxiv.org/abs/2110.08399) | **yes** | `2110.08399.pdf` | ¹⁰B/¹¹B first measurement; still no Li |
| **E. Quasi-elastic (e,e′) — the RC tail's dominant knob** |
| E1 | `benhar2008-quasielastic-review` | O. Benhar, D. Day, I. Sick, *Inclusive quasi-elastic electron–nucleus scattering*, Rev. Mod. Phys. **80** (2008) 189–224 | [arXiv:nucl-ex/0603029](https://arxiv.org/abs/nucl-ex/0603029) | **yes** | `nucl-ex_0603029.pdf` | `qe_kf_gev`, `qe_suppression`, Q9 |
| E2 | `benhar2006-qes-archive-note` | O. Benhar, D. Day, I. Sick, *An archive for quasi-elastic electron–nucleus scattering data* | [arXiv:nucl-ex/0603032](https://arxiv.org/abs/nucl-ex/0603032) | **yes** | `nucl-ex_0603032.pdf` | the citation for E3 |
| E3 | `qes-archive-6li` | QES archive ⁶Li dataset, `6Li.dat` — 133 points, three settings, source `Heimlich:1973` | [discovery.phys.virginia.edu/…/qes-archive/data/6Li.dat](https://discovery.phys.virginia.edu/research/groups/qes-archive/data/6Li.dat) | **yes, free ASCII** | not committed (data) | turns `RC_QE_KF_GEV` from an input into a fit residual |
| E4 | `heimlich1974-6li-12c-electron-scattering` | F. H. Heimlich, M. Köbberling, J. Moritz, K. H. Schmidt, D. Wegener, D. Zeller, J. K. Bienlein, J. Bleckwenn, *High-energy electron scattering from ⁶Li and ¹²C*, Nucl. Phys. **A231** (1974) 509–520 | [doi:10.1016/0375-9474(74)90514-4](https://doi.org/10.1016/0375-9474(74)90514-4) | no | — | the **primary source** of every point in E3 |
| E5 | `moniz1971-fermi-momenta` | E. J. Moniz, I. Sick, R. R. Whitney, J. R. Ficenec, R. D. Kephart, W. P. Trower, *Nuclear Fermi Momenta from Quasielastic Electron Scattering*, Phys. Rev. Lett. **26** (1971) 445–448 | [doi:10.1103/PhysRevLett.26.445](https://doi.org/10.1103/PhysRevLett.26.445) | no | — | **`RC_QE_KF_GEV` = 0.169** (C-10) — the number's own paper |
| E6 | `whitney1974-quasielastic` | R. R. Whitney, I. Sick, J. R. Ficenec, R. D. Kephart, W. P. Trower, *Quasielastic electron scattering*, Phys. Rev. C **9** (1974) 2230–2235 | [doi:10.1103/PhysRevC.9.2230](https://doi.org/10.1103/PhysRevC.9.2230) | no | — | the full paper behind E5's k_F fit; large-angle ⁶Li data E3 lacks |
| **F. Moments, radii, the α core, ab initio** |
| F1 | `stone2016-adndt-quadrupole` | N. J. Stone, *Table of nuclear electric quadrupole moments*, At. Data Nucl. Data Tables **111–112** (2016) 1–28 | [doi:10.1016/j.adt.2015.12.002](https://doi.org/10.1016/j.adt.2015.12.002) | no | — (free IAEA twin in corpus) | `LI6_QUADRUPOLE_FM2`, `LI7_QUADRUPOLE_FM2` |
| F2 | `pyykko2008-quadrupole-moments` | P. Pyykkö, *Year-2008 nuclear quadrupole moments*, Mol. Phys. **106** (2008) 1965–1974 | [doi:10.1080/00268970802018367](https://doi.org/10.1080/00268970802018367) · OA [hal-00513184](https://hal.science/hal-00513184) | **green** | — (bot-walled, §11) | the −0.0806 fm² that T11 explicitly is **not** |
| F3 | `borremans2005-li-moments` | D. Borremans, D. L. Balabanski, K. Blaum, W. Geithner, S. Gheysen, P. Himpe, M. Kowalska, J. Lassen, P. Lievens, S. Mallion, R. Neugart, G. Neyens, N. Vermeulen, D. Yordanov, *New measurement and reevaluation of the nuclear magnetic and quadrupole moments of ⁸Li and ⁹Li*, Phys. Rev. C **72** (2005) 044309 | [doi:10.1103/PhysRevC.72.044309](https://doi.org/10.1103/PhysRevC.72.044309) | **bronze** | — | the measurement Stone's Q(⁶Li) = −0.000806(6) b traces to |
| F4 | `george-knutson1999-eta` | E. A. George, L. D. Knutson, *Determination of the ⁶Li → α + d asymptotic D- to S-state ratio by a restricted phase shift analysis*, Phys. Rev. C **59** (1999) 598–606 | [doi:10.1103/PhysRevC.59.598](https://doi.org/10.1103/PhysRevC.59.598) | no | — | **C-1**, `LI6_ETA_DS_GK` — and see §10.4 |
| F5 | `pastore2013-em-moments-qmc` | S. Pastore, S. C. Pieper, R. Schiavilla, R. B. Wiringa, *Quantum Monte Carlo calculations of electromagnetic moments and transitions in A ≤ 9 nuclei…*, Phys. Rev. C **87** (2013) 035503 | [arXiv:1212.3375](https://arxiv.org/abs/1212.3375) | **yes** | `1212.3375.pdf` | the GFMC Q(⁶Li) = −0.20(6) fm² quoted against C-2 |
| F6 | `carlson2015-qmc-review` | J. Carlson, S. Gandolfi, F. Pederiva, S. C. Pieper, R. Schiavilla, K. E. Schmidt, R. B. Wiringa, *Quantum Monte Carlo methods for nuclear physics*, Rev. Mod. Phys. **87** (2015) 1067 | [arXiv:1412.3081](https://arxiv.org/abs/1412.3081) | **yes** | `1412.3081.pdf` | the method review for `data/vmc/` |
| F7 | `pieper-wiringa2001-qmc-light-nuclei` | S. C. Pieper, R. B. Wiringa, *Quantum Monte Carlo Calculations of Light Nuclei*, Ann. Rev. Nucl. Part. Sci. **51** (2001) 53–90 | [arXiv:nucl-th/0103005](https://arxiv.org/abs/nucl-th/0103005) | **yes** | `nucl-th_0103005.pdf` | same, at review length |
| F8 | `pudliner1997-qmc-a-le-7` | B. S. Pudliner, V. R. Pandharipande, J. Carlson, S. C. Pieper, R. B. Wiringa, *Quantum Monte Carlo calculations of nuclei with A ≤ 7*, Phys. Rev. C **56** (1997) 1720–1750 | [arXiv:nucl-th/9705009](https://arxiv.org/abs/nucl-th/9705009) | **yes** | `nucl-th_9705009.pdf` | **E-13** — the AV18+UIX Hamiltonian paper `data/vmc/` cites |
| F9 | `forest1996-toroidal-structures` | J. L. Forest, V. R. Pandharipande, S. C. Pieper, R. B. Wiringa, R. Schiavilla, A. Arriaga, *Femtometer toroidal structures in nuclei*, Phys. Rev. C **54** (1996) 646–667 | [arXiv:nucl-th/9603035](https://arxiv.org/abs/nucl-th/9603035) | **yes** | `nucl-th_9603035.pdf` | **E-13** — the two-nucleon densities in ⁶Li/⁷Li |
| F10 | `wiringa1995-av18` | R. B. Wiringa, V. G. J. Stoks, R. Schiavilla, *An accurate nucleon–nucleon potential with charge-independence breaking*, Phys. Rev. C **51** (1995) 38–51 | [arXiv:nucl-th/9408016](https://arxiv.org/abs/nucl-th/9408016) | **yes** | `nucl-th_9408016.pdf` | the AV18 interaction under every `data/vmc/` file |
| F11 | `frosch1967-4he-structure` | R. F. Frosch, J. S. McCarthy, R. E. Rand, M. R. Yearian, *Structure of the He⁴ Nucleus from Elastic Electron Scattering*, Phys. Rev. **160** (1967) 874–879 | [doi:10.1103/PhysRev.160.874](https://doi.org/10.1103/PhysRev.160.874) | no | — | `AlphaCoreSource`; the `Fr67` row of A6's `tabletr.txt` |
| F12 | `ottermann1985-3he-4he-elastic` | C. R. Ottermann, G. Köbschall, K. Maurer, K. Röhrich, Ch. Schmitt, V. H. Walther, *Elastic electron scattering from ³He and ⁴He*, Nucl. Phys. **A436** (1985) 688–698 | [doi:10.1016/0375-9474(85)90554-8](https://doi.org/10.1016/0375-9474(85)90554-8) | no | — | `AlphaCoreSource`, the low-q half of the SOG fit |
| F13 | `camsonne2014-4he-form-factor` | A. Camsonne, A. T. Katramatou, M. Olson, N. Sparveris *et al.*, *JLab measurement of the ⁴He charge form factor at large momentum transfers*, Phys. Rev. Lett. **112** (2014) 132503 | [arXiv:1309.5297](https://arxiv.org/abs/1309.5297) | **yes** | `1309.5297.pdf` | `AlphaCoreSource` above the SOG's q range |
| **G. Polarized ⁶Li target material** (T5 context) |
| G1 | `goertz2002-polarized-targets-review` | St. Goertz, W. Meyer, G. Reicherz, *Polarized H, D and ³He targets for particle physics experiments*, Prog. Part. Nucl. Phys. **49** (2002) 403–489 | [doi:10.1016/S0146-6410(02)00159-X](https://doi.org/10.1016/S0146-6410(02)00159-X) | no | — | the EST / spin-temperature closed form behind `spin_temperature_pzz` |
| G2 | `goertz1995-6lid-polarization` | St. Goertz, Ch. Bradtke, H. Dutz, R. Gehring, W. Meyer, M. Plückthun, G. Reicherz, K. Runkel, *Investigations in high temperature irradiated ⁶,⁷LiH and ⁶LiD, its dynamic nuclear polarization and radiation resistance*, Nucl. Instrum. Meth. A **356** (1995) 20–28 | [doi:10.1016/0168-9002(94)01437-X](https://doi.org/10.1016/0168-9002(94)01437-X) | no | — | the ⁶LiD material itself |
| G3 | `goertz2004-dnp-process` | S. T. Goertz, *The dynamic nuclear polarization process*, Nucl. Instrum. Meth. A **526** (2004) 28–42 | [doi:10.1016/j.nima.2004.03.147](https://doi.org/10.1016/j.nima.2004.03.147) | no | — | EST derivation; BENCHMARK_PLAN §4 item 9 |

### 1.3 The free machine-readable archive (A6), read for this survey

The **UVa Nuclear Charge Density archive** collects the fit parameters of
At. Data Nucl. Data Tables **14** (1974), **36** (1987) and **60** (1995)
and ships C++ evaluators for them. Three of its files matter here. All
three were fetched and their lithium/helium rows read in this pass.

* **`dldata/FB_data.dat`** (18 306 B, 81 rows). Fourier–Bessel coefficients.
  The **⁶Li** row is
  `6Li 6 3  1.6353e-02 2.9603e-02 2.0807e-02 7.2731e-03 4.7580e-04 −8.7510e-04 1.1413e-03 … 6.0`
  — **seven non-zero coefficients and a cut-off R = 6.0 fm**. (The first two
  integer columns are labelled `Z A` but are filled `A Z`; the ⁴He row is
  mis-keyed `4 4`. Parse by the atom string, not by Z.) Integrating it here:
  **∫ρ dV = 2.991** against Z = 3 (0.3 %, the FB truncation), **r_rms =
  2.5206 fm**, and the first two zeros of the point-charge F_ch at
  **q = 2.6944 and 3.3774 fm⁻¹**. The first of those is the number
  `BENCHMARK_PLAN.md` §4 item 5 quotes as 2.695 — confirmed independently
  here, and it sits **outside** T11's asserted [2.9, 3.3] fm⁻¹ window, on
  the same side as Li *et al.*'s measured diffraction minimum (2.83 fm⁻¹).
* **`dldata/tabletr.txt`** (274 rows). The analytic-model parameter table.
  Row 5 is the one this project has been looking for:
  `7Li 3 7 1 1.77 0.327 0 2.39 Su67` — model 1 = **harmonic oscillator**,
  a = **1.77 fm**, α = **0.327**, w = 0, ⟨r²⟩^½ = **2.39 fm**, from Suelzle
  1967. **There is no ⁶Li row anywhere in the file** (§10.1).
* **`dldata/SOG_data.dat`**. Sum-of-Gaussians. The **⁴He** row gives
  rms **1.676 fm**, γ = **1.0 fm**, and **12 (R_i, Q_i) pairs** whose
  **ΣQ_i = 1.000009** — i.e. the published parameterisation is complete and
  self-normalised as shipped.

The archive also mirrors scans of the three ADNDT source volumes at
`dldata/1974_source.pdf`, `1987_source.pdf`, `1995_source.pdf` (HTTP 200,
4.1 / 3.0 / 4.1 MB). Those are scans of paywalled Elsevier articles hosted
by a third party; they were **deliberately not downloaded into the corpus**
(§11). The *data* files are facts and the archive publishes them as a
research resource, so they are the route to use.

---

## 2. The elastic form factor and the `HoSpin1FF` open item (Q1)

### 2.1 Suelzle, Yearian & Crannell 1967 (A1) — and what its abstract says

100–600 MeV electrons on ⁶Li and ⁷Li at Stanford HEPL; maximum momentum
transfer squared **6.9 fm⁻²** (q ≤ 2.63 fm⁻¹). This is the origin of both
lithium charge form factors and of the ⁷Li parameters in every compilation
since. Its abstract — retrieved and read in this pass, which
`03_data_nuclear.md` could not do from behind the APS paywall — states
four things the tree should act on:

> "The charge form factor for Li⁷ could be interpreted in terms of a simple
> harmonic-oscillator shell model with a quadrupole contribution described
> by the undeformed p-shell model. For this model, the rms charge radius was
> 2.39 ± 0.03 F, and the electric-quadrupole moment required for the fit was
> in excellent agreement with the spectroscopic measurements. **The Li⁶
> charge form factor could not be interpreted either in terms of the simple
> harmonic-well shell model, or in terms of the modified harmonic-well model**
> in which the s and p nucleons are permitted to move in potential wells of
> different strengths. The rms radius obtained for a phenomenological fit to
> the data was 2.54 ± 0.05 F. The ratio of the Li⁶ and Li⁷ rms radii … was
> 1.055 ± 0.008."

`HoSpin1FF` is a harmonic-oscillator charge form factor, and LiPolGen ships
it **for ⁶Li** with `LI6_FF_HO_A_FM = 1.9069`, `LI6_FF_HO_ALPHA = 0.13822`,
while `HoSpin1FF::for_ion` *throws* for ⁷Li. The primary paper says the
model works for ⁷Li and does not work for ⁶Li. That is the cleanest
available statement of why `rc.hpp`'s "TO BE REFIT" note exists, it explains
why the compilations carry a ⁷Li HO row and no ⁶Li one (§1.3, §10.1), and it
argues that the refit `rc.hpp` asks for should be **⁷Li's** — where the
parameters are published, free and machine-readable — with ⁶Li kept as an
explicitly-labelled effective shape or replaced by the FB density of §1.3.
Note also that Suelzle's q range stops *below* the ⁶Li diffraction minimum,
so A1 alone cannot settle T11's `q₀` window; A2 can.

### 2.2 Li, Sick, Whitney & Yearian 1971 (A2) — already a corpus entry

500 MeV, q² to 13 fm⁻², the ⁶Li diffraction minimum at q² = 8 fm⁻²
(q ≈ 2.83 fm⁻¹). Already in `refs_dict.json` with `file: null`; this pass
adds the **DOI, 10.1016/0375-9474(71)90257-0**, which the corpus entry
lacks, and confirms there is no free copy anywhere (OpenAlex: closed; no
arXiv; the DTIC report version is not served to this environment).

### 2.3 The compilations (A3, A4, A5, A7) and the method paper (A8)

`de Vries 1987` (A3) is the standard charge-density compilation and the
source of the ⁶Li MI30 rows (2.54(5) fm from Su67, 2.56(5) fm from Li71a)
and the ⁷Li HO row. `de Jager 1974` (A4) is its predecessor and carries the
other half the tree needs — the **magnetization**-density parameters, i.e.
the compiled form of B1. **Neither is open access** (both verified closed),
and **neither's magnetization table is in the UVa machine-readable files** —
`FB_data.dat`, `SOG_data.dat` and `tabletr.txt` are all charge-only. So for
`F_m` there is no free machine-readable route and A4 or B1 must be bought or
read in a library. `Fricke 1995` (A5) is the third compilation, is reported
green open access through the UvA repository, and is a useful cross-check on
`LI6_R2_POINT_FM2`; `Angeli & Marinova 2013` (A7) is the source of the
2.589(39) fm the tree actually uses. `Sick 1974` (A8) is the method paper
for the FB analysis — it defines the incompleteness error that should be
carried on any gate built on §1.3's coefficients.

### 2.4 Theory that the tree already leans on (A9) and the cluster fit (A10)

`Wiringa & Schiavilla 1998` (A9, **downloaded**) is the VMC AV18+UIX
calculation of ⁶Li's longitudinal and transverse elastic form factors whose
monopole *shape* `rc.hpp` borrows; its abstract (verified) confirms one- and
two-body charge and current operators constructed consistently with the NN
interaction, and "good agreement with available experimental data" — the
figures-only caveat the tree already records stands. `Bergstrom 1979` (A10)
is the complementary object and is **not** in the tree anywhere: a
phenomenological **α–d cluster-model** fit to ⁶Li's electromagnetic form
factors. LiPolGen models ⁶Li as α + d *and* carries a separate
harmonic-oscillator form factor; A10 is the published attempt to make those
two the same object, and it is the natural external check on
`ClusterConfigSampler`'s charge distribution.

---

## 3. Magnetization — the half of the form factor with no data behind it

`rc.hpp` builds `F_m(0) = (M_A/m_p)·LI6_MU_N = 4.90765` and then has nothing
measured to constrain its q dependence. Four papers cover this, none free.

**B1, Rand, Frosch & Yearian 1966** is the source for both isotopes at once
— elastic scattering from the magnetic multipole distributions of ⁶Li, ⁷Li,
⁹Be, ¹⁰B, ¹¹B and ¹⁴N (title verified; Phys. Rev. 144, 859–873). This is
`BENCHMARK_PLAN.md` §4 item 6, and it is the *only* measured electromagnetic
signature of α+d clustering in the magnetic sector. **B2, Bergstrom,
Kowalski & Neuhausen 1982** supersedes its ⁶Li half at higher momentum
transfer (PRC 25, 1156–1167) and is the better number to gate against if
only one is bought. **B3, van Niftrik, Lapikás, de Vries & Box 1971** is the
⁷Li magnetization distribution from 180° scattering, and it determines
Q(⁷Li) *from electron scattering* — an independent route to
`LI7_QUADRUPOLE_FM2 = −4.06`, which today rests on TUNL alone. **B4,
Lichtenstadt *et al.* 1989** extends ⁷Li's longitudinal and transverse form
factors to high q, past where Suelzle's HO fit was made.

---

## 4. The α–d momentum distribution — the tagged channel's own observable

`ClusterChannel::momentum_density`, `Wave::radial` and
`data/vmc/li6_alpha_d/li6.ad` encode the α–d relative momentum distribution;
**E-13 validates them against ANL VMC only**, i.e. theory against theory.
It has been measured four independent ways, and this pass pins every one to
a DOI.

**C1/C2 (Ent 1994 / 1986)** are the (e,e′d) measurements. C2's *title* is
literally "the α–d Momentum Distribution in the Ground State of ⁶Li"
(verified through Crossref), and C1 is the full NPA paper covering ⁴He, ⁶Li
and ¹²C together. **C3 (Mitchell 1991)** is the α-knockout mirror; its
abstract, read in this pass, states verbatim: *"The ⁶Li(e,e′α) reaction has
been measured in parallel kinematics. The Q² dependence of both the
two-body and three-body breakup of ⁶Li has been studied at fixed recoil
momentum. In addition, the recoil momentum dependence of the two-body α-d
breakup has been investigated. The data indicate that the reaction mechanism
is quasielastic in these kinematics."* By momentum conservation that is the
same relative-momentum distribution LiPolGen's T1 channel samples, with the
detected and undetected fragments exchanged.

**C4 (Genin 1974)** is the first-generation version, and here the tree's
citation must be fixed: `03_data_nuclear.md` lists it as "Phys. Lett. B **51**
(1974), PII 037026937490714X". That PII resolves to
`10.1016/0370-2693(74)90714-X` = Genin *et al.*, Phys. Lett. B **52** (1974)
46–48. Phys. Lett. B **51** (1974) 217 is a *different* paper — Toyama &
Sakamoto, "⁶Li(e,e′α) and (e,e′d) reactions", which is **theory**. Both were
verified; the volume number in the survey is wrong and the two papers are
easy to conflate.

**C5 (Albrecht 1980, ⁶Li(p,pd)⁴He at 670 MeV)** and **C6 (Dollhopf 1975,
⁶Li(α,2α)²H at 700 MeV)** are the hadronic probes; C5 is where the
factor-two width disagreement (FWHM ≈ 70 vs 120 MeV/c) lives, and that
disagreement is the honest external band on the distribution's width.

**C7 (Connelly 1998)** is the caution, not the benchmark: the trinucleon
decomposition of the same nucleus, reporting a significant mirror-channel
deviation at low momentum transfer, i.e. two-step contamination in cluster
knockout. Any use of C1–C6 as a *wave-function* measurement inherits that
systematic. It is the one pre-arXiv paper in this domain with an open copy
(§5.5).

**C8 (Lapikás, Wesseling & Wiringa 1999, downloaded)** is the single
cleanest item in this file. ⁷Li(e,e′p)⁶He momentum distributions from −70
to 260 MeV/c compared with **VMC wave functions from the same author and the
same Hamiltonian family as `data/vmc/`**; the abstract (verified) says the
VMC calculations "provide a parameter-free prediction of the momentum
distribution that reproduces the measured data, including its
normalization", with a summed spectroscopic factor **0.58 ± 0.05** against
the VMC value. Citing it converts E-13 from "the transform convention is
right" to "the wave-function family has been compared with data by its own
author, on the sister isotope, at the 8 % level". Note the journal spells
the middle author **Wesseling**; the arXiv listing spells it *Wessling*.

**C9 (Hotta 1999, downloaded)** is ⁶Li(e,e′p) at low momentum transfer in
the giant-resonance region — relevant only to the low-ω end of the RC tail,
and listed so the negative ("it does not constrain the cluster overlap") is
on the record with its source.

---

## 5. EMC and DIS on lithium

### 5.1 The one measurement (D1)

**NMC, Nucl. Phys. B441 (1995) 12, `hep-ex/9504002`, downloaded.** The
abstract, verified: *"We present the structure function ratios F2(Li)/F2(D)
and F2(C)/F2(D) measured in deep inelastic muon-nucleus scattering at a
nominal incident muon energy of 200 GeV. The kinematic range 0.0001 < x <
0.7 and 0.01 < Q² < 70 GeV² is covered. For values of x less than 0.002 both
ratios indicate saturation of shadowing at values compatible with
photoabsorption results."* The target really is lithium and the paper says
so in its §3, read from the downloaded PDF in this pass: *"an upstream one
with Li and D targets and a downstream one with C and D targets … one
lithium and one deuterium target in the upstream set"*. HEPData record
`ins394050` is live (HTTP 200). This is the only ⁶Li DIS measurement in
existence and it is the external constraint on
`EMC_VALENCE_DEPLETION_EPPS21`.

### 5.2 The A-dependence twin (D2)

**NMC, Nucl. Phys. B441 (1995) 3, `hep-ph/9503291`, downloaded.** Abstract
verified: the re-evaluation presents "the ratios F2(C)/F2(Li), F2(Ca)/F2(Li)
and F2(Ca)/F2(C) measured at 90 GeV", kinematic range "0.0085 < x < 0.6,
0.84 < Q² < 17 GeV² for the Li/C/Ca ones" — the same ⁶Li target, and a Q²
range in the valence region better matched to the Q² = 5 GeV² at which the
tree's constant is defined. HEPData `ins393377` is live.

### 5.3 The negatives, with sources (D3, D4, D5)

D3 (E139) and D5 (E12-10-008) are the A-dependence programmes that skipped
lithium; both verified, both cited here so the negative has a reference
rather than an assertion. D4 (Seely 2009, downloaded) is the ⁹Be result —
two α plus a neutron, an EMC effect tracking local rather than average
density. It is the nearest *measured* statement that cluster structure
matters for the EMC effect, and the argument transfers to ⁶Li's α + d
qualitatively and only qualitatively.

### 5.4 Where a volume number was wrong

See §4 on Genin/Toyama. Listed here too because it is the sort of error that
propagates: `03_data_nuclear.md` §3 row 3 should read **Phys. Lett. B 52
(1974) 46**.

### 5.5 Connelly's open copies

`10.1103/PhysRevC.57.1569` is reported open access (bronze at APS,
`http://link.aps.org/pdf/10.1103/PhysRevC.57.1569`; green at the VU Research
Portal, record `4d7aadf3-5212-4ce0-9322-dae82436721d`, file
`https://research.vu.nl/files/1797325/144572.pdf`). Both refused an
automated fetch here (HTTP 403 from the portal's bot filter). A browser will
get it; see §11.

---

## 6. Quasi-elastic (e,e′) — the RC tail's dominant knob

### 6.1 The review and the archive note (E1, E2)

Both downloaded. E1 (RMP 80, 189) is the review, and its abstract states it
"includes a compilation of data available in numerical form" — that
compilation is E3. E2 is the two-page note that exists to be cited for the
archive itself.

### 6.2 The ⁶Li dataset (E3), re-counted here

`6Li.dat` was fetched in this pass and parsed: **133 data rows**, columns
`Z A E[GeV] θ[deg] ω[GeV] σ δσ source`, split **40 / 47 / 46** across the
three settings **(2.5 GeV, 12.0°)**, **(2.7 GeV, 13.8°)**, **(2.7 GeV,
15.0°)**, every row carrying the single source key `Heimlich:1973`. This
reproduces `03_data_nuclear.md`'s description exactly and confirms the file
is still served. It is free, ASCII and absolutely normalised, and it is the
route by which `RC_QE_KF_GEV` stops being an input and becomes a residual.

### 6.3 The primary sources (E4, E5, E6)

E4 (Heimlich 1974) is the paper every point in E3 comes from — cite it, not
just the archive. E5 (Moniz 1971) is the paper `RC_QE_KF_GEV = 0.169` is
taken from, and E6 (Whitney 1974) is its full-length companion, containing
the 500 MeV / 60° large-angle ⁶Li data that the archive's ⁶Li page does
**not** carry. All three verified, all three closed.

---

## 7. Moments, radii and the α core

### 7.1 Quadrupole moments (F1, F2, F3)

The corpus already holds Stone's free IAEA reports; F1 is the *journal*
version of the quadrupole table (ADNDT 111–112 (2016) 1, closed), listed for
citation completeness. F2 (Pyykkö 2008) is the compilation whose
Q(⁶Li) = −0.0806 fm² T11 explicitly tests *against* — the tree asserts
TUNL's −0.0818, not Pyykkö's — and it is green open access at HAL. F3
(Borremans 2005) is the primary measurement Stone's ⁶Li entry traces to
(`2005Bo45`), and it is bronze open access at APS. Together these close the
provenance chain `LI6_QUADRUPOLE_FM2` → TUNL/Stone → the measurement.

### 7.2 The asymptotic D/S ratio (F4)

`George & Knutson 1999` — title verified exactly as *"Determination of the
⁶Li → α + d asymptotic D- to S-state ratio by a restricted phase shift
analysis"*, PRC 59, 598–606. Note the title confirms
`06_critic.md`'s refutation: it is a **phase-shift analysis**, not a direct
measurement of d+α tensor analysing powers, which is exactly why
BENCHMARK_PLAN §4 item 8 asks for it to be retitled in the tree. §10.4
records that no newer determination exists.

### 7.3 The α core (F11, F12, F13, and A6's SOG row)

`AlphaCoreSource::VmcHe4Density` / `he4.density` is validated today against
ANL VMC alone. The published ⁴He charge form factor is available **for free
and in closed form**: the sum-of-Gaussians of §1.3, 12 (R_i, Q_i) pairs,
rms 1.676 fm, ΣQ_i = 1.000009, machine-readable in `SOG_data.dat`. Its
underlying data are F11 (Frosch 1967) and F12 (Ottermann 1985); F13
(Camsonne 2014, downloaded) extends the measurement to 29 ≤ Q² ≤ 77 fm⁻²,
beyond where the SOG fit is valid. `tabletr.txt` additionally carries two
analytic ⁴He rows (`Fr67`: c = 1.008, z = 0.327, w = 0.445, rms 1.71;
`MC74`: 0.964, 0.322, 0.517, 1.71) as an independent cross-check on the SOG.
This is the cheapest data-for-theory swap available anywhere in this file.

### 7.4 The ab-initio chain under `data/vmc/` (F5–F10)

E-13 names "Pudliner *et al.* PRC 56 (1997) 1720" and "Forest *et al.* PRC
54 (1996) 646" as the provenance of the wave functions LiPolGen reads, and
neither was in the corpus. Both are on arXiv and are now downloaded (F8,
F9), together with the AV18 paper they build on (F10), the QMC review
(F6, F7), and the GFMC electromagnetic-moments paper (F5) that supplies the
Q(⁶Li) = −0.20(6) fm² the tree quotes against C-2. With these the chain
`data/vmc/*` → AV18+UIX/UX → published Hamiltonian → published method is
citable end to end without leaving the corpus.

---

## 8. Polarized ⁶Li target material (G1–G3)

⁶LiD is the material every polarized-lithium experiment has used, and the
corpus already holds the COMPASS, SLAC and Saclay target papers. What it
lacks is the **physics of the polarization mechanism**: G1 (Goertz, Meyer &
Reicherz 2002) is the standard review, G3 (Goertz 2004) is the DNP-process
paper where the equal-spin-temperature relation the tree implements as
`spin_temperature_pzz` is derived, and G2 (Goertz 1995) is the ⁶LiD material
paper — irradiation, DNP and radiation resistance. All three are closed
access. They matter to BENCHMARK_PLAN §4 item 9, whose honest half is that
**P_zz of ⁶Li has never been measured by anyone**; nothing found in this
pass changes that.

---

## 9. What this domain does *not* constrain

Stated once, because every item above is unpolarized:

1. **Nothing here touches b₁, A_zz or any tensor observable of lithium.**
   C1–C6 are unpolarized cross sections: they constrain |φ₀|² + |φ₂|²
   summed, never the S–D interference b₁ is built from. The 1.93× η
   discrepancy (C-1) is untouched by all of it.
2. **No ⁷Li DIS measurement exists**, so the CBT polarized-EMC curve the
   tree carries has no ⁷Li unpolarized data under it either.
3. **No separated ⁶Li C2 form factor exists.** Longitudinal elastic data
   measure F_C0² + F_C2², and C2 is sub-dominant where those data are good.
   `fq_scale` remains the honest carrier.
4. **Absolute normalisations from cluster knockout carry C7's systematic.**

---

## 10. Verified negatives

**10.1 There is no ⁶Li harmonic-oscillator row in any compilation.**
`tabletr.txt` (274 rows, all eight analytic models of ADNDT 14/36/60)
contains **one** lithium row and it is ⁷Li's. Combined with §2.1's abstract
this is not an omission — it is the compilers recording that the model does
not fit ⁶Li.

**10.2 There is no ⁷Li Fourier–Bessel or sum-of-Gaussians density.**
`FB_data.dat` has exactly one lithium row (⁶Li); `SOG_data.dat` has no
lithium at all. So the data-derived C0 route of §1.3 exists for ⁶Li only,
and ⁷Li has only the HO parameters of §1.3.

**10.3 The magnetization compilations are not in the machine-readable
archive.** All three UVa data files are charge-only. B1/A4 must be read from
the journal.

**10.4 No determination of η(⁶Li → α+d) has been published since George &
Knutson 1999.** Searched via INSPIRE (`t "asymptotic normalization" and t
"6Li"`, and a date-restricted variant), Crossref bibliographic search and
the arXiv API. Everything returned is the *astrophysical* S-wave ANC
programme for α + d → ⁶Li extracted from d(α,γ)⁶Li capture (Blokhintsev,
Mukhamedzhanov *et al.*, e.g. EPJ Web Conf. **227** (2020) 02016,
[doi:10.1051/epjconf/202022702016](https://doi.org/10.1051/epjconf/202022702016);
Int. J. Mod. Phys. Conf. Ser. **49** (2019),
[doi:10.1142/S2010194519600176](https://doi.org/10.1142/S2010194519600176)),
which determines the **S-wave normalisation** for a Big-Bang S-factor
extrapolation and does not quote the D/S ratio. F4 stands alone, and C-1's
band is as good as the external world gets.

**10.5 Almost none of the pre-1990 nuclear literature is open access.**
Checked individually through OpenAlex: A1, A2, A3, A4, A7, A8, A10, B1–B4,
C1–C6, D3, E4, E5, E6, F1, F11, F12, G1–G3 are all `closed`. The exceptions
are C7 (bronze + green), F2 (green), F3 (bronze) and A5 (green). Everything
on arXiv is of course open, and that is the entire D/E1–E2/F5–F10 block.

---

## 11. What was downloaded, and what has to be bought

**Committed to `PolarizedLithiumSim/refs/` in this pass** (all arXiv, all
open access, filenames in the house style):

| file | entry |
|---|---|
| `hep-ex_9504002.pdf` | D1 NMC F₂(⁶Li)/F₂(D) |
| `hep-ph_9503291.pdf` | D2 NMC A-dependence with the ⁶Li reference |
| `nucl-th_9904008.pdf` | C8 Lapikás–Wesseling–Wiringa ⁷Li(e,e′p) vs VMC |
| `nucl-th_9807037.pdf` | A9 Wiringa–Schiavilla ⁶Li form factors |
| `nucl-ex_9810016.pdf` | C9 Hotta ⁶Li(e,e′p) |
| `nucl-ex_0603029.pdf` | E1 Benhar–Day–Sick RMP review |
| `nucl-ex_0603032.pdf` | E2 the QES archive note |
| `1212.3375.pdf` | F5 Pastore GFMC moments |
| `1412.3081.pdf` | F6 Carlson QMC review |
| `nucl-th_0103005.pdf` | F7 Pieper–Wiringa QMC review |
| `nucl-th_9705009.pdf` | F8 Pudliner AV18+UIX QMC A ≤ 7 |
| `nucl-th_9603035.pdf` | F9 Forest two-nucleon densities |
| `0904.4448.pdf` | D4 Seely EMC in light nuclei |
| `2110.08399.pdf` | D5 Arrington EMC light and heavy |
| `1309.5297.pdf` | F13 Camsonne ⁴He form factor |

(`nucl-th_9408016.pdf`, F10 AV18, was already present.)

**Deliberately not downloaded.** The ADNDT scans at
`discovery.phys.virginia.edu/research/groups/ncd/dldata/{1974,1987,1995}_source.pdf`
are third-party scans of paywalled Elsevier articles; they are reachable but
they are not publisher open access, so they are linked and not committed.
The UVa `.dat` files and the QES `6Li.dat` are data, not papers, and belong
under `validation/benchmarks/data/` with a provenance header when a gate is
wired — not in the paper corpus.

**Blocked by bot filters, not by a paywall** (a browser will fetch these):
C7 Connelly at `https://research.vu.nl/files/1797325/144572.pdf` or
`http://link.aps.org/pdf/10.1103/PhysRevC.57.1569`; F2 Pyykkö at
`https://hal.science/hal-00513184`; F3 Borremans at
`http://link.aps.org/pdf/10.1103/PhysRevC.72.044309`.

**Must be bought or read in a library** — the DOI page is the download link
in every case: A1, A2, A3, A4, A7, A8, A10, B1, B2, B3, B4, C1, C2, C3, C4,
C5, C6, D3, E4, E5, E6, F1, F11, F12, G1, G2, G3. Of these the six worth
paying for first are **A1** (the ⁷Li HO fit and the ⁶Li statement), **A4 or
B1** (the only magnetization data), **C1 + C3** (the α–d distribution), and
**E4** (the ⁶Li quasi-elastic cross sections behind the archive).

---

## 12. Priority

| rank | action | reference | effort | what it buys |
|---|---|---|---|---|
| 1 | Fit a **⁷Li `HoSpin1FF`** on the published HO parameters, straight from `tabletr.txt` | A6 / A1 | hours, not days | `HoSpin1FF::for_ion(LI7)` stops throwing, on *published, fitted* parameters in the shipped functional form — and ⁶Li's unfitted (a, α) finally get an external comparison |
| 2 | Load the **⁶Li Fourier–Bessel density** as a data-derived C0 and re-open T11's `q₀` window | A6 / A8 | hours | a form factor from data instead of a model; the zero at **2.6944 fm⁻¹** and the measured minimum at 2.83 fm⁻¹ both sit outside the asserted [2.9, 3.3] |
| 3 | Gate the **⁴He core** on the free sum-of-Gaussians | A6 / F11 / F12 | ~1 hour | `AlphaCoreSource` stops being theory-checked-by-theory; no transcription needed |
| 4 | Wire **HEPData `ins394050`** with the isoscalar denominator | D1 | ½ day | the only ⁶Li DIS measurement, against `EMC_VALENCE_DEPLETION_EPPS21` |
| 5 | Cite **C8** in E-13 and **E4/E5** in the `RC_QE_KF_GEV` block | C8, E4, E5 | minutes | the VMC family and the k_F constant both get their own paper |
| 6 | Add the **QES archive ⁶Li 133 points** as a quasi-elastic gate | E3 / E4 | 1–2 days | `RC_QE_KF_GEV` becomes a fit residual at three real ⁶Li settings |
| 7 | Fix the **Genin volume number** in `03_data_nuclear.md` | C4 | minutes | PLB 52 (1974) 46, not PLB 51 |
| 8 | Buy **A4 or B1** and put a q dependence under `F_m` | A4, B1, B2 | ½–1 day after purchase | the magnetization form factor has no data behind it at all today |
