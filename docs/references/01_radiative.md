# 01 — Radiative corrections: the references the tail rests on

**Domain:** the QED radiative corrections LiPolGen applies (`RcModel`, `src/core/rc.cpp`,
`include/lipolgen/rc.hpp`) — the tensor BAND `rc_delta`, the elastic/quasi-elastic
radiative TAIL on its three models (`RcTailModel::TPeak`, `TPeakPlusLL`, `PolradFull`),
and the pieces the tree has priced but not computed (the polarised quasi-elastic tail,
the A = 2 → A = 6 transfer of the band, the tagged `rc_tail ≡ 1` exclusion).

**Pass date:** 2026-09-16. **Written by:** the reference-research stage, radiative domain.

---

## 0. What this pass did, and what it did not

*Method.* Every entry below was verified over the network in this pass — arXiv API
(`export.arxiv.org/api/query`), the INSPIRE literature API
(`inspirehep.net/api/literature`), or a DOI resolution — and the title, author list,
year, journal reference and abstract were read before the claim in the "what it
contains" paragraph was written. Where a PDF was downloaded, its first page was
re-read with `pdftotext` to confirm the file is the paper it claims to be. Entries
that could **not** be verified, or for which no free copy exists, say so in the row
and in §8; nothing is listed silently.

*Two results worth putting first.*

1. **`[MT69]` is obtainable and is now on disk.** The tree records Mo–Tsai, Rev. Mod.
   Phys. **41** (1969) 205 as paywalled-with-no-arXiv and the `tests/test_rc.cpp` gate
   against it as "deliberately not written" (`benchmarking/00` §9 item 5, `04_theory.md`
   T-35, `references/00_corpus.md` G1). The INSPIRE record for that exact article
   (recid 52657) carries a **freely served fulltext**: `SLAC-PUB-0380` (January 1968),
   whose title page reads *"RADIATIVE CORRECTIONS TO ELASTIC AND INELASTIC ep AND μp
   SCATTERING … (To be submitted to Rev. Mod. Phys.)"* — i.e. the preprint of the RMP
   article itself, Section III's exact elastic radiative tail, Appendix B's Eq. (B.5)
   and the paper's own section *"Comparisons of Various Versions of Peaking
   Approximations with the Exact …"* all present. It is committed as
   `refs/SLAC-PUB-0380.pdf` (2.79 MB). **Caveat to carry into any gate:** the equation
   numbers in the preprint need not coincide with the published RMP ones; cite the gate
   as "SLAC-PUB-380 Eq. (B.5)" unless someone checks the RMP text, and keep the RMP DOI
   in the citation line.

2. **The polarised quasi-elastic neglect now has its paper.** `G7` — Z.-L. Zhou *et al.*,
   PRL **82** (1999) 687, the deuteron measurement both HERMES and this generator's
   `qe_tensor_scale = 0` default lean on — is `arXiv:nucl-ex/9809002`, open access, and
   is now `refs/nucl-ex_9809002.pdf`. Its actual title is *"Tensor Analyzing Powers for
   Quasi-Elastic Electron Scattering from Deuterium"*: a **measurement at ⟨q⟩ = 1.7 fm⁻¹
   with non-zero A_zz**, not a null result. Anyone re-opening `phase_B_numbers.md` §B3
   should read it before re-asserting "no net tensor effect"; HERMES's own sentence
   (quoted verbatim in §7.1 below) is the claim, and Zhou is its citation, not its proof.

*Downloaded this pass (19 new files, all open access, all under
`PolarizedLithiumSim/refs/`):* `SLAC-PUB-0380.pdf`, `SLAC-PUB-0848.pdf`,
`SLAC-PUB-1528.pdf`, `JINR-P2-10113.pdf`, `hep-ph_9706516.pdf`, `nucl-ex_9809002.pdf`,
`hep-ph_0403262.pdf`, `hep-ph_9906408.pdf`, `hep-ph_0106180.pdf`, `hep-ph_0105032.pdf`,
`hep-ph_0102086.pdf`, `nucl-th_0002058.pdf`, `1202.2225.pdf`, `2306.14578.pdf`,
`2008.02895.pdf`, `2210.03785.pdf`, `2012.09970.pdf`, `2212.04730.pdf`, `2505.23487.pdf`.

*Present on disk but fetched by a sibling domain in the same pass* (listed here because
this domain needs them, not claimed as this domain's downloads): `hep-ex_0506018.pdf`
(HERMES b₁), `hep-ex_0609039.pdf` (HERMES g₁ long paper), `nucl-ex_0002003.pdf`
(Abbott deuteron form factors), `2506.04506.pdf` (Poudel *et al.*). All four are already
in `refs_dict.json` under the sibling domain's own keys.

*Not appended to `refs_dict.json` by this domain.* Parallel domain agents were writing
into `refs/` while this ran — the dictionary grew from 54 to 92 entries at 11:40 during
this pass — and a concurrent read-modify-write would race. The **31** ready-to-merge
entries for everything below are in `docs/references/01_radiative_refs_dict_entries.json`,
schema-matched to the existing `entries` map (`file`, `arxiv`, `journal`, `title`,
`authors`, `year`, `keywords`, `key_content`, `used_in`, optional `url` /
`repo_usage_note`), for one atomic merge after a fresh read. Four of this domain's
references — R4, R29, R30, R34 — were appended by a sibling domain during the pass, under
the keys `hermes2005-b1-tensor`, `hermes2007-g1`, `abbott2000-ff-parametrization` and
`poudel2025-tensor-sf-review`; the merge file lists them separately so nothing is
duplicated.

---

## 1. Table 1 — the radiative corpus, in one table

Priority: **H** = closes a documented in-tree gap or gate; **M** = bounds something the
tree ships as a choice; **L** = context, cite-only.

| # | Key | Reference | Id | OA | File | Serves (RcModel piece) | Pri |
|---|---|---|---|---|---|---|---|
| R1 | `motsai1969-rc-elastic-inelastic` | L. W. Mo, Y.-S. Tsai, RMP **41** (1969) 205 | DOI 10.1103/RevModPhys.41.205 | article no; **preprint yes** | `refs/SLAC-PUB-0380.pdf` | the unwritten `tests/test_rc.cpp` `[MT69]` gate; the exact-vs-peaking bound on `RcTailModel::TPeak` | **H** |
| R2 | `tsai1971-slac-pub-848` | Y.-S. Tsai, SLAC-PUB-848 (1971) | INSPIRE 67278 | yes | `refs/SLAC-PUB-0848.pdf` | the same gate; `04_theory.md` T-36's action item | **H** |
| R3 | `polrad2-1997` | Akushevich, Ilyichev, Shumeiko, Soroko, Tolkachev, CPC **104** (1997) 201 | `hep-ph/9706516` | yes | `refs/hep-ph_9706516.pdf` | **everything**: Eqs. (9)/(10), (18), (37)–(39), (43), (44), (A.1)–(A.4), App. B — the whole transcribed sector | **H** |
| R4 | `hermes2005-b1-tensor` | HERMES, A. Airapetian *et al.*, PRL **95** (2005) 242001 | `hep-ex/0506018` | yes | `refs/hep-ex_0506018.pdf` | `RC_DELTA_LOW_X_OPTIMISTIC` (0.19, from Table II); the RC-treatment precedent | **H** |
| R5 | `zhou1999-tensor-analyzing-powers` | Z.-L. Zhou *et al.*, PRL **82** (1999) 687 | `nucl-ex/9809002` | yes | `refs/nucl-ex_9809002.pdf` | the sole citation under `qe_tensor_scale = 0.0` | **H** |
| R6 | `gakh-shekhovtsova2004-tensor-rc` | G. I. Gakh, O. Shekhovtsova, JETP **99** (2004) 898 | `hep-ph/0403262` | yes | `refs/hep-ph_0403262.pdf` | `RC_DELTA_LOW_X` = 0.30 and its rejected SHAPE | **H** |
| R7 | `bardin-shumeiko1977-exact-ir` | D. Yu. Bardin, N. M. Shumeiko, NPB **127** (1977) 242 | DOI 10.1016/0550-3213(77)90213-9 | article no; **JINR preprint yes (Russian)** | `refs/JINR-P2-10113.pdf` | the covariant IR separation POLRAD Eq. (18) and `polrad_tpeak_quadrature` are built on | M |
| R8 | `akushevich-shumeiko1994-light-nuclei` | I. V. Akushevich, N. M. Shumeiko, J. Phys. G **20** (1994) 513 | DOI 10.1088/0954-3899/20/4/001 | **no** | — | POLRAD's polarized-**light-nucleus** predecessor: the paper the A = 2 → A = 6 transfer should be argued from | **H** |
| R9 | `akushevich2001-polarized-rc-review` | I. Akushevich, A. Ilyichev, N. Shumeiko, review (2001) | `hep-ph/0106180` | yes | `refs/hep-ph_0106180.pdf` | the code map (POLRAD/RADGEN/HAPRAD/DIFFRAD/MASCARAD) and the covariant formulae in one place | M |
| R10 | `maximon-tjon2000-rc-ep` | L. C. Maximon, J. A. Tjon, PRC **62** (2000) 054320 | `nucl-th/0002058` | yes | `refs/nucl-th_0002058.pdf` | the size of what Mo–Tsai's soft-photon/peaking treatment misses → a bound on `TPeak` | M |
| R11 | `tsai1961-rc-ep` | Y.-S. Tsai, Phys. Rev. **122** (1961) 1898 | DOI 10.1103/PhysRev.122.1898 | **no** | — | historical root of the elastic-tail formula; cite-only | L |
| R12 | `tsai1974-pair-bremsstrahlung` | Y.-S. Tsai, RMP **46** (1974) 815 (+ erratum RMP **49** (1977) 421) | INSPIRE 81153 (`SLAC-PUB-1365`) | yes (not fetched) | — | the equivalent-radiator whose single-z form is `ll_radiator`; external brems is collider-irrelevant | L |
| R13 | `kuraev-fadin1985-structure-function` | E. A. Kuraev, V. S. Fadin, Sov. J. Nucl. Phys. **41** (1985) 466 | INSPIRE 217313 | **no** | — | the citation `ll_radiator` does not carry: D(z) and its **exponentiated** form | M |
| R14 | `gakh-konchatnij-merenkov2012-elastic-ed` | G. I. Gakh, M. I. Konchatnij, N. P. Merenkov, JETP **115** (2012) 212 | `1202.2225` | yes | `refs/1202.2225.pdf` | the **only independent** tensor-target RC calculation besides R6 — checks POLRAD Eq. (A.4) | **H** |
| R15 | `gakh-shekhovtsova2003-proceedings` | G. I. Gakh, O. Shekhovtsova, eConf **C030626** (2003) | `hep-ph/0309123` | yes (not fetched) | — | the conference precursor of R6; useful only to date the calculation | L |
| R16 | `afanasev2001-coincidence-rc` | A. V. Afanasev, I. Akushevich, G. I. Gakh, N. P. Merenkov, JETP **93** (2001) 449 | `hep-ph/0105032` | yes | `refs/hep-ph_0105032.pdf` | the **tagged** channels' `rc_tail ≡ 1`: RC in coincidence, polarized, model-independent | **H** |
| R17 | `afanasev2001-mascarad` | A. Afanasev, I. Akushevich, N. Merenkov, PRD **64** (2001) 113009 | `hep-ph/0102086` | yes | `refs/hep-ph_0102086.pdf` | polarized **elastic** RC with realistic acceptance — the acceptance question inside Eq. (44) | M |
| R18 | `radgen1999` | I. Akushevich, H. Böttcher, D. Ryckbosch, RADGEN 1.0 (1999) | `hep-ph/9906408` | yes | `refs/hep-ph_9906408.pdf` | the engine behind HERMES's b₁ RC → what R4's 2 × 10⁻³ residual actually means | **H** |
| R19 | `heracles1992` | A. Kwiatkowski, H. Spiesberger, H.-J. Möhring, CPC **69** (1992) 155 | DOI 10.1016/0010-4655(92)90136-M | **no** | — | the collider-kinematics O(α) ep generator under DJANGOH | M |
| R20 | `django6-1994` | K. Charchula, G. A. Schuler, H. Spiesberger, CPC **81** (1994) 381 | DOI 10.1016/0010-4655(94)90086-8 | **no** | — | `BENCHMARK_PLAN.md` §4 item 7's Rad/noRad comparison | **H** |
| R21 | `djangoh-manual` | H. Spiesberger, DJANGOH homepage + manual (v4.6.x) | uni-mainz.de | yes (web) | — | the maintained descendant; ≥ 4.6.10 does **polarized** NC/CC | M |
| R22 | `pepsi1992` | L. Mankiewicz, A. Schäfer, M. Veltri, CPC **71** (1992) 305 | DOI 10.1016/0010-4655(92)90016-R | **no** | — | the polarized-leptoproduction MC the plan pairs with RADGEN (T4) | M |
| R23 | `esfrad2022` | A. Afanasev, I. Akushevich, A. Ilyichev, N. Merenkov, ESFRAD (2022) | `2212.04730` | yes | `refs/2212.04730.pdf` | **higher-order** QED vs POLRAD/RADGEN — the missing exponentiation in `TPeakPlusLL` | **H** |
| R24 | `byer2022-sidis-rc-evgen` | D. Byer, V. Khachatryan, H. Gao, I. Akushevich, A. Ilyichev *et al.*, CPC **287** (2023) 108702 | `2210.03785` | yes | `refs/2210.03785.pdf` | the maintained modern implementation of the same O(α) formalism; the "Byer" reference | M |
| R25 | `afanasev2023-rc-review` | A. Afanasev, J. C. Bernauer, P. Blunden, J. Blümlein *et al.*, EPJ A topical review (2023) | `2306.14578` | yes | `refs/2306.14578.pdf` | the modern baseline for "what 1.5 % on A_zz means"; which codes are alive | M |
| R26 | `liu2021-factorized-rc` | T. Liu, W. Melnitchouk, J.-W. Qiu, N. Sato, PRD **104** (2021) 094033 | `2008.02895` | yes | `refs/2008.02895.pdf` | whether a fixed-target O(α) tail ports to collider kinematics at all | **H** |
| R27 | `cammarota2025-factorized-qed-qcd` | J. Cammarota, J.-W. Qiu, K. Watanabe, J.-Y. Zhang (2025) | `2505.23487` | yes | `refs/2505.23487.pdf` | NLO joint QED⊗QCD factorization — the successor framework to R26 | L |
| R28 | `cfns2020-rc-whitepaper` | A. Afanasev *et al.*, CFNS ad-hoc RC whitepaper (2020) | `2012.09970` | yes | `refs/2012.09970.pdf` (21.7 MB) | the EIC-specific RC problem statement, incl. nuclear targets | L |
| R29 | `hermes2007-g1` | HERMES, A. Airapetian *et al.*, PRD **75** (2007) 012007 | `hep-ex/0609039` | yes | `refs/hep-ex_0609039.pdf` | HERMES's RC method in full ("Unfolding of Radiative and Instrumental Smearing" + App. A) | **H** |
| R30 | `abbott2000-ff-parametrization` | D. Abbott *et al.* (JLab t₂₀), EPJ A **7** (2000) 421 | `nucl-ex/0002003` | yes | `refs/nucl-ex_0002003.pdf` | R4's ref. [11]: the deuteron FF parametrization for the coherent/QE tails; the tree's own "Abbott … not used" finding | **H** |
| R31 | `stein1975-slac-4deg` | S. Stein *et al.*, PRD **12** (1975) 1884 | DOI 10.1103/PhysRevD.12.1884 | article no; **preprint yes** (`SLAC-PUB-1528`) | `refs/SLAC-PUB-1528.pdf` | R4's ref. [12]: the inclusive elastic + quasi-elastic radiative-tail recipe | **H** |
| R32 | `moniz1971-fermi-momenta` | E. J. Moniz, I. Sick, R. R. Whitney, J. R. Ficenec, R. D. Kephart, W. P. Trower, PRL **26** (1971) 445 | DOI 10.1103/PhysRevLett.26.445 | **no** | — | `RC_QE_KF_GEV = 0.169` — the knob the QE tail is most sensitive to | M |
| R33 | `deforest-walecka1966` | T. de Forest Jr., J. D. Walecka, Adv. Phys. **15** (1966) 1 | DOI 10.1080/00018736600101254 | **no** | — | `pauli_suppression`'s S(q) = (3/4)u − u³/16 | M |
| R34 | `poudel2025-tensor-sf-review` | J. Poudel, A. Bacchetta, J.-P. Chen, N. Santiesteban, EPJ A **61** (2025) 81 | `2506.04506` | yes | `refs/2506.04506.pdf` | `RC_X_HIGH` = 0.16 (the published kinematic edge); does **not** contain the 1.5 % | M |
| R35 | `polrad-cpc-library-adgh` | CPC Program Library entry `ADGH_v1_0` (POLRAD 2.0 source) | catalogue id only | source in-tree | — | the FORTRAN `polrad_transcription_check.md` checks against | L |

---

## 2. Tier A — the two gates this closes

### 2.1 R1 · Mo & Tsai (1969), and its free preprint

> L. W. Mo, Y.-S. Tsai, *Radiative Corrections to Elastic and Inelastic ep and μp
> Scattering*, Rev. Mod. Phys. **41** (1969) 205. DOI
> [10.1103/RevModPhys.41.205](https://doi.org/10.1103/RevModPhys.41.205) (APS, paywalled;
> 1090 INSPIRE citations). **Free preprint:** L. W. Mo, Y. S. Tsai, SLAC-PUB-380,
> January 1968 — [INSPIRE recid 52657](https://inspirehep.net/literature/52657),
> PDF [inspirehep.net/files/1fcaa81f63f50d7bf56a22ce2c6b8b58](https://inspirehep.net/files/1fcaa81f63f50d7bf56a22ce2c6b8b58)
> (HTTP 200, `application/pdf`, 2 791 721 B). Committed as `refs/SLAC-PUB-0380.pdf`.

*What it contains.* The paper's own abstract: it "investigated and improved the
reliability of many formulae used in the radiative corrections to elastic and inelastic
electron scatterings when only the scattered electrons are detected", and gives "a
practical and reliable recipe for unfolding the entire inelastic spectra, including
effects due to virtual photons, internal and external bremsstrahlungs". Concretely, and
checked in the downloaded text: Section III computes the **elastic radiative tail with
the exact lowest-order bremsstrahlung formula**, Appendix B gives that exact expression
(Eq. (B.5) in the preprint's numbering, written in terms of the target's G_E/G_M), and a
later section is titled *"Comparisons of Various Versions of Peaking Approximations with
the Exact"* — i.e. the paper measures its own peaking approximation against its own exact
answer, on the elastic peak and on the 3-3 resonance.

*Why the project needs it.* `tests/test_rc.cpp` names `[MT69]` as the gate that was
deliberately left unwritten because the paper could not be obtained
(`benchmarking/00` §9 item 5, re-stated in `run_2026-09-03/SUMMARY.md` item 5 and
`phase_B_numbers.md`: *"Mo–Tsai is still not obtained, so the absolute normalisation of
`rc_tail` rests on POLRAD alone"*). It is the single most-cited RC reference in the
field and the only external normalisation the elastic tail can be checked against that
is not POLRAD's own compiled numbers.

*What it bounds or validates.* Two distinct things, and they should not be conflated.
(a) **Normalisation** — the σ^el of `polrad_sigma_el_u` at Q_N = 0 against Mo–Tsai's
exact elastic tail: an independent-source check on the number `RcModel::tail_ratio_at`
divides by. (b) **The peaking approximation itself** — the tree ships `TPeak` (one of
three peaks) as an explicit LOWER BOUND, with `TPeakPlusLL` and `PolradFull` disagreeing
by factors up to ×6444 in the low-y corner (`rc.hpp` lines ~60–180, `phase_B_numbers.md`
§B2). Mo–Tsai's own exact-vs-peaking comparison is the classical statement of exactly
that error, on a target where both sides are known; it is the right place to calibrate
how much of the ×6444 is physics and how much is `TPeakPlusLL` breaking down.

*Caveat.* SLAC-PUB-380 is the **preprint**; its equation numbers may differ from the
published RMP. Cite the gate against the preprint's numbering explicitly, and keep the
RMP DOI in the bibliography line so a reader with APS access can cross-check.

### 2.2 R2 · Tsai, SLAC-PUB-848 (1971)

> Y.-S. Tsai, *Radiative Corrections to Electron Scatterings*, SLAC-PUB-848 (January
> 1971). [INSPIRE recid 67278](https://inspirehep.net/literature/67278) (93 citations),
> PDF [inspirehep.net/files/3ce239706be17b0eefec145433700c64](https://inspirehep.net/files/3ce239706be17b0eefec145433700c64)
> (4 112 884 B). Committed as `refs/SLAC-PUB-0848.pdf`. No DOI, no journal version.

*What it contains.* From the record's own abstract: radiative corrections to electron
scattering from nucleons **and nuclei** above ~50 MeV; "many formulas in the Mo and
Tsai's article in Review of Modern Physics are improved and better derivations of them
are presented"; straggling from external bremsstrahlung and ionization is folded in; and
a method is proposed for high-Z targets.

*Why the project needs it.* `04_theory.md` T-36 already identified this as the free
substitute for R1 and recorded the action *"obtain Tsai SLAC-PUB-848 and retire the
'Mo–Tsai not obtained' status"* (§"Ten things to do", item 8, ½ day). It is done. With
R1's preprint also in hand, both halves of that action are closed and the gate can be
written against the better-derived forms.

*What it bounds or validates.* The same `polrad_sigma_el_u` normalisation as R1, plus —
because Tsai treats **nuclei**, not only the proton — the Z-dependence question that the
A = 2 → A = 6 transfer (`RcOptions::a_transfer_frac`, default 0.0, registry row B6)
prices without a source. Note that the straggling/external-bremsstrahlung half of the
paper is **not** applicable: the EIC is a collider, there is no target material, and
`RcModel` correctly implements internal radiation only.

---

## 3. Tier B — the O(α) chain POLRAD sits on

### 3.1 R3 · POLRAD 2.0 — the code the whole sector is transcribed from

> I. Akushevich, A. Ilyichev, N. Shumeiko, A. Soroko, A. Tolkachev, *POLRAD 2.0. FORTRAN
> code for the Radiative Corrections Calculation to Deep Inelastic Scattering of
> Polarized Particles*, Comput. Phys. Commun. **104** (1997) 201–244.
> [arXiv:hep-ph/9706516](https://arxiv.org/abs/hep-ph/9706516), DOI
> [10.1016/S0010-4655(97)00062-3](https://doi.org/10.1016/S0010-4655(97)00062-3).
> Open access. Committed as `refs/hep-ph_9706516.pdf`.

*What it contains.* Abstract, verbatim: "The FORTRAN code POLRAD 2.0 for radiative
correction calculation in inclusive and semi-inclusive deep inelastic scattering of
polarized leptons by polarized nucleons **and nuclei** is described. Its theoretical
basis, structure and algorithms are discussed in details." For this project specifically
(per `open_items/physics_literature.md` §3, `04_theory.md` T-33 and the in-tree
transcription check): Born Eqs. (9)/(10) with the −1/3 : +1/6 = P₂(cos Θ) tensor
geometry; Eqs. (A.1)/(A.3) the eight generalized structure functions ℑ₁…ℑ₈ with ℑ₅–ℑ₈
the quadrupole ones in b₁–b₄; **Eq. (A.4) the tensor elastic radiative tail via the
quadrupole form factor F_q**; Eq. (18) + Appendix B the exact τ_A-integrated tail;
Eqs. (37)–(39) + (43) the t-peak closed forms; Eq. (44) the quasi-elastic tail;
Eqs. (110)/(111) the iterative b₁ unfolding.

*Why the project needs it.* It is a **surprising gap**: `references/00_corpus.md`
records POLRAD as "in-tree, D-3" — the source is compiled and checked against — while
the *paper* had no PDF and no corpus entry. Every equation number quoted in
`rc.hpp`, `design_C_tensor_rc.md`, `polrad_transcription_check.md` and
`phase_B_numbers.md` points at a document nobody could open from the repository.

*What it bounds or validates.* It does not bound anything; it **is** the implementation
target for `polrad_sigma_el_u`, `polrad_sigma_el_t`, `polrad_sigma_qe_u`,
`polrad_im_el_spin1`, `polrad_full_sigma_el`, `polrad_full_sigma_qe_u`,
`polrad_tpeak_quadrature`, `polrad_eta_limits`/`polrad_tau_limits` and
`rosenbluth_spin1`. Note the tree's own finding that `polrad2t.tex`'s Appendix B is
**wrong in five places against POLRAD's own FORTRAN** (`polrad_transcription_check.md`
§7 preamble) — so the paper is necessary but not sufficient, and the FORTRAN (R35)
remains the arbiter.

### 3.2 R7 · Bardin & Shumeiko (1977) — the covariant infrared separation

> D. Yu. Bardin, N. M. Shumeiko, *An Exact Calculation of the Lowest Order
> Electromagnetic Correction to the Elastic Scattering*, Nucl. Phys. B **127** (1977)
> 242. DOI [10.1016/0550-3213(77)90213-9](https://doi.org/10.1016/0550-3213(77)90213-9)
> (Elsevier, paywalled). **Free preprint:** JINR-P2-10113,
> [INSPIRE recid 111694](https://inspirehep.net/literature/111694), PDF
> [inspirehep.net/files/367ff4ac8f411583c055fa795497645a](https://inspirehep.net/files/367ff4ac8f411583c055fa795497645a).
> Committed as `refs/JINR-P2-10113.pdf` — **the preprint is in Russian** (Cyrillic title
> page, "О точном вычислении электромагнитной поправки низшего порядка к упругому
> рассеянию"); the English text is the paywalled NPB article.

*What it contains.* The covariant procedure for separating the infrared-divergent part
of the lowest-order bremsstrahlung cross section, and the exact evaluation of that part
for experiments that measure the energy (or energy and angle) of one final particle and
cannot distinguish elastic from bremsstrahlung events.

*Why the project needs it.* This is the method POLRAD is built on — the reason POLRAD's
tail is stated as an exact integral with a covariantly-subtracted IR piece rather than a
peaking formula. `rc.hpp` describes `polrad_tpeak_quadrature` as "POLRAD's η_A
quadrature of Eqs. (37)–(39), (43): exact in the peaking sense"; the distinction between
"exact" and "exact in the peaking sense" is precisely Bardin–Shumeiko's subtraction.

*What it bounds or validates.* The structural claim behind `RcTailModel::PolradFull`'s
status as "ONE exact τ_A quadrature, with the lepton mass in". If anyone needs to argue
that `PolradFull` is not merely a third model but the controlled O(α) answer, this is
the paper the argument is made from.

### 3.3 R8 · Akushevich & Shumeiko (1994) — polarized DIS off polarized **light nuclei**

> I. V. Akushevich, N. M. Shumeiko, *Radiative effects in deep inelastic scattering of
> polarized leptons by polarized light nuclei*, J. Phys. G **20** (1994) 513–530. DOI
> [10.1088/0954-3899/20/4/001](https://doi.org/10.1088/0954-3899/20/4/001). **No arXiv,
> no free copy located** (INSPIRE recid 383295 carries no fulltext and no report number;
> checked this pass). IOP paywall — the DOI page above is the link to use.

*What it contains.* Per the abstract: the principal QED radiative effects in DIS of
polarized leptons off polarized H, **D and ³He**, at Born level and with RC; the
transversely-polarized-target case; everything in covariant variables; exact low-order
results; numerical analysis for the then-forthcoming polarization experiments; and a
discussion of "two points of view on calculation of RC to experimental data" — i.e. the
correct-the-data vs. correct-the-model choice `rc.hpp`'s opening comment also makes.

*Why the project needs it.* This is POLRAD's own theory paper **for nuclear targets**,
and it is the closest published statement to what LiPolGen does: a polarized lepton on a
polarized light nucleus, with the nuclear (not nucleon) invariants. The tree's most
expensive unpriced assumption — carrying a deuteron-derived δ(x) to A = 6
(`RcOptions::a_transfer_frac`, default 0.0, "a defensible, conservative EXTRAPOLATION"
per `rc.hpp`) — has no cited source; this paper is where the A-dependence of the
polarized RC is actually worked out, for D and ³He, and is the right place to start.

*What it bounds or validates.* `a_transfer_frac`'s default of 0 and the band's implicit
A-independence; secondarily, the transverse-target formulae behind the `Q_N → P_zz^eff`
substitution `rc.hpp` describes at every θ_S through POLRAD Eq. (43).

### 3.4 R9 · Akushevich, Ilyichev & Shumeiko (2001) — the review

> I. Akushevich, A. Ilyichev, N. Shumeiko, *Radiative effects in scattering of polarized
> leptons by polarized nucleons and light nuclei*,
> [arXiv:hep-ph/0106180](https://arxiv.org/abs/hep-ph/0106180) (2001). Open access; no
> journal reference on the arXiv record. Committed as `refs/hep-ph_0106180.pdf`.

*What it contains.* A review of the covariant-approach formulae for inclusive,
semi-inclusive, diffractive and elastic polarized scattering, with the FORTRAN codes
POLRAD, RADGEN, HAPRAD, DIFFRAD and MASCARAD described and applications to CERN, DESY,
SLAC and TJNAF data shown numerically.

*Why the project needs it.* One document that maps formula → code → experiment for the
entire Minsk-school RC family. When `rc.hpp` says a term is "POLRAD Eq. (44)" or that
"POLRAD supplies no σ_T counterpart at the s- or p-peak", this is the fastest way to
confirm the statement is about the *formalism* and not about one code's options.

*What it bounds or validates.* Nothing numerically; it is the orientation document for
anyone extending `RcModel` beyond the inclusive channel (the diffractive and
semi-inclusive sections are the coherent-⁶Li and tagged analogues).

### 3.5 R10 · Maximon & Tjon (2000) — how much Mo–Tsai's treatment leaves out

> L. C. Maximon, J. A. Tjon, *Radiative Corrections to Electron-Proton Scattering*,
> Phys. Rev. C **62** (2000) 054320. [arXiv:nucl-th/0002058](https://arxiv.org/abs/nucl-th/0002058),
> DOI [10.1103/PhysRevC.62.054320](https://doi.org/10.1103/PhysRevC.62.054320). Open
> access. Committed as `refs/nucl-th_0002058.pdf`.

*What it contains.* Elastic ep RC in a hadronic model with nucleon finite size; the soft-
photon contribution calculated **exactly**; an explicit comparison with "the generally
used expressions previously obtained by Mo and Tsai"; and the finding that above 8 GeV
at large angles the proton vertex correction raises the Rosenbluth factor by ≥ 2 %.

*Why the project needs it.* It is the standard modern quantification of the error in the
Mo–Tsai soft-photon/peaking package — the package `TPeak` and `TPeakPlusLL` inherit. The
tree currently has **no external number** for "how wrong is the peaking approximation";
it has only the internal spread between its own three models.

*What it bounds or validates.* The systematic floor under `RcTailModel::TPeak` and
`TPeakPlusLL`, and the ≈ 2 % scale at which a hadronic-structure vertex correction enters
— useful context for the statement that the whole tail is a small dilution on `A_zz`
(the band leads it by ×768 at x = 0.01).

### 3.6 R13 · Kuraev & Fadin (1985) — the citation `ll_radiator` does not carry

> E. A. Kuraev, V. S. Fadin, *On Radiative Corrections to e⁺e⁻ Single Photon
> Annihilation at High Energy*, Sov. J. Nucl. Phys. **41** (1985) 466 [Yad. Fiz. **41**
> (1985) 733]. [INSPIRE recid 217313](https://inspirehep.net/literature/217313). **No
> DOI, no arXiv, no free copy located** (verified this pass).

*What it contains.* The electron structure-function (equivalent-radiator) method: the
leading-log radiator D(z, Q²) ∝ ln(Q²/m_e²)(1+z²)/(1−z), its soft-photon
**exponentiated** form, and the higher orders that follow from iterating it. **Flagged,
not asserted:** no free copy of this paper was obtained, so its text was not read in
this pass; the coefficient convention differs between sources (α/2π vs α/π, and whether
the −1 in β = (2α/π)(ln(Q²/m²) − 1) is kept), and `ll_radiator` uses α/π. Whoever writes
the citation should pin the convention against a source they have actually read, and
satisfy themselves that the α/π in `rc.hpp` is the intended normalisation and not a
factor-2 convention slip — this pass did not check it.

*Why the project needs it.* `ll_radiator` in `rc.hpp` implements exactly
D(z) = (α/π) ln(Q²/m_e²)(1+z²)/(1−z) and its own comment states the limitation
honestly: *"It is a SINGLE-z LEADING LOG and nothing more: no soft-photon exponent, no
non-log O(α) piece, no second emission… AND ITS SOFT 1/(1−z) IS UNCANCELLED"*, reaching
13.5 at x = 0.744, y = 0.0071. The exponentiated form is the standard cure, and it has a
name and a source; the code cites neither.

*What it bounds or validates.* `RcTailModel::TPeakPlusLL`'s breakdown corner, and — with
R23 (ESFRAD), which is the same method applied to polarized ep — supplies the estimate
of what exponentiation would change.

---

## 4. Tier C — tensor-target radiative corrections (all of them)

The honest statement in `open_items/physics_literature.md` §3 is that **one** dedicated
analytic tensor-DIS RC calculation exists. That remains true for **deep-inelastic**
scattering. This pass found a **second** tensor-polarized-deuteron RC calculation for
**elastic** ed scattering, by an overlapping group, which is the only external check
POLRAD Eq. (A.4)'s tensor elastic tail can be given.

### 4.1 R6 · Gakh & Shekhovtsova (2004)

> G. I. Gakh, O. Shekhovtsova, *Radiative corrections to deep-inelastic ed scattering.
> Case of tensor polarized deuteron*, J. Exp. Theor. Phys. **99** (2004) 898–914
> (the Russian ZhETF volume is given inconsistently by the two records — arXiv's
> `journal_ref` says Zh. Eksp. Teor. Fiz. **99** (2004) 1034–1050, INSPIRE recid 647050
> says Zh. Eksp. Teor. Fiz. **123** (2004) 1034; cite the English JETP reference).
> [arXiv:hep-ph/0403262](https://arxiv.org/abs/hep-ph/0403262),
> DOI [10.1134/1.1842872](https://doi.org/10.1134/1.1842872). Open access. Committed as
> `refs/hep-ph_0403262.pdf`. **INSPIRE recid 647050, citation count 0** (re-verified
> 2026-09-15 per `run_2026-09-03/STATUS.md` row 7).

*What it contains.* Abstract, verbatim: "The model-independent radiative corrections to
deep-inelastic scattering of unpolarized electron beam off the tensor polarized deuteron
target have been considered. The contribution to the radiative corrections due to the
hard-photon emission from the elastic electron-deuteron scattering (the so-called
elastic radiative tail) is also investigated. The calculation is based on the covariant
parametrization of the deuteron quadrupole polarization tensor. The numerical estimates
of the radiative corrections to the polarization observables have been done for the
kinematical conditions of the current experiment at HERA." Appendix A tabulates the
hard-photon coefficients; four figure panels at Q² = 0.1, 1, 4, 10 GeV².

*Why the project needs it.* It is the sole source of `RC_DELTA_LOW_X = 0.30`, the low-x
anchor of the whole tensor band — and the tree already holds a detailed, corrected
reading of what that number is (`rc.hpp` lines ~326–393, `phase_B_numbers.md` §B6.3):
the 10 %/30 % pair is the two ends of **one panel's x window at Q² = 0.1**, |δ| = 0.1133
at x = 0.00966 and 0.2662 at x = 0.00226, and the shipped 0.30 is ×2.65 the panel's
reading at the x nearest the anchor. Having the PDF in the corpus is what lets the next
reader check that reading instead of inheriting it.

*What it bounds or validates.* `rc_delta`'s low-x anchor and `RcOptions::delta_low_x`;
and, on the record, **what it does not**: the SHAPE is a documented rejection (singular
across the b₁ zero-crossing, −1.691…+7.810 at Q² = 4; non-monotone; no panel covering
x = 0.01, Q² = 5; a ratio against *their* Born).

### 4.2 R14 · Gakh, Konchatnij & Merenkov (2012) — the second tensor calculation

> G. I. Gakh, M. I. Konchatnij, N. P. Merenkov, *Radiative Corrections to Polarization
> Observables in Elastic Electron–Deuteron Scattering in Leptonic Variables*, J. Exp.
> Theor. Phys. **115** (2012) 212. [arXiv:1202.2225](https://arxiv.org/abs/1202.2225),
> DOI [10.1134/S1063776112070060](https://doi.org/10.1134/S1063776112070060). Open
> access, 13 INSPIRE citations. Committed as `refs/1202.2225.pdf`.

*What it contains.* Model-independent QED radiative corrections to polarization
observables in elastic scattering of unpolarized and longitudinally-polarized electrons
off a deuteron target, in **leptonic variables**, with the target **arbitrarily
polarized** and an explicit procedure for applying the results to the vector *or tensor*
polarization; built on the Drell–Yan-like representation (the same electron
structure-function machinery as R16/R23).

*Why the project needs it.* `rc.hpp`'s tensor elastic tail is POLRAD Eq. (A.4) — a
single source, transcribed once, with the tree's own record of a spurious η_A found in
the design document (§B5) and five Appendix-B errors found in the paper against its own
FORTRAN. There is **no second implementation** of the tensor elastic tail anywhere in
the project. This paper is an independent analytic treatment of the same physical object
(hard-photon emission off elastic ed with a tensor-polarized target, in leptonic
variables — POLRAD's variables) by a different group and a different method.

*What it bounds or validates.* `polrad_im_el_spin1` / `polrad_full_sigma_el`'s tensor
half, i.e. the `t` component of `PolradFullPair`, and therefore `sp_tensor_scale`'s
bound on the two non-exact models. It is a deuteron calculation, so it validates the
**formalism and the elastic-vertex structure**, not the ⁶Li form factors; that is the
right division of labour, since the ⁶Li C0/C2 shapes are already banded
(`C0Shape::{Ho, VmcFt}`).

*Caveat.* Elastic kinematics in **leptonic** variables; extracting a number comparable
to Eq. (A.4)'s tail requires mapping their observables onto POLRAD's ℑ^el basis. This is
a real piece of work (estimate: 2–4 days), not a table lookup.

### 4.3 R16 · Afanasev, Akushevich, Gakh & Merenkov (2001) — RC in **coincidence**

> A. V. Afanasev, I. Akushevich, G. I. Gakh, N. P. Merenkov, *Radiative Corrections to
> Polarized Inelastic Scattering in Coincidence*, J. Exp. Theor. Phys. **93** (2001)
> 449–461. [arXiv:hep-ph/0105032](https://arxiv.org/abs/hep-ph/0105032), DOI
> [10.1134/1.1410589](https://doi.org/10.1134/1.1410589). Open access. Committed as
> `refs/hep-ph_0105032.pdf`.

*What it contains.* Abstract: "The complete analysis of the model-independent leading
radiative corrections to cross-section and polarization observables in semi-inclusive
deep-inelastic **electron-nucleus** scattering with detection of a proton and scattered
electron **in coincidence**", from the Drell–Yan-like representation for both the
spin-independent and spin-dependent parts, with applications to polarization transfer
and to a polarized target.

*Why the project needs it.* The tagged channels ship `rc_tail ≡ 1` exactly. `rc.hpp` and
`phase_B_numbers.md` §B4 are explicit that this is "half a kinematic fact — the elastic
half is a veto, the quasi-elastic half is an **omission** of the same order as the
inclusive quasi-elastic dilution", i.e. a documented exclusion with no calculation
behind it. This paper is the calculation for the coincidence topology: it is what says
how a detected-hadron coincidence changes the radiative kernel, on a **nucleus**, with
polarization.

*What it bounds or validates.* `rc_tail_applies(channel, with_tail) == false` for every
`Tagged*` channel — the one place the generator asserts a correction is exactly 1.0.
If any number in this project deserves a bound rather than a paragraph, it is that one.

### 4.4 R17 · Afanasev, Akushevich & Merenkov (2001) — MASCARAD

> A. Afanasev, I. Akushevich, N. Merenkov, *Model independent radiative corrections in
> processes of polarized electron–nucleon elastic scattering*, Phys. Rev. D **64** (2001)
> 113009. [arXiv:hep-ph/0102086](https://arxiv.org/abs/hep-ph/0102086), DOI
> [10.1103/PhysRevD.64.113009](https://doi.org/10.1103/PhysRevD.64.113009). Open access.
> Committed as `refs/hep-ph_0102086.pdf`.

*What it contains.* Explicit RC formulae for elastic ep with two typical polarization
measurements (beam-target asymmetry; recoil polarization), an explicit treatment of
**realistic experimental acceptances**, and the FORTRAN code MASCARAD.

*Why the project needs it.* Two uses. (a) The acceptance discussion is the published
treatment of the question `phase_B_numbers.md` raises and leaves open: *"The tag
acceptance inside Eq. (44). POLRAD Eq. (44) is an inclusive [formula]"* — a radiative
tail computed inclusively but then compared against a tagged or binned sample.
(b) It is the elastic, polarized, nucleon-level member of the same formula family that
`rosenbluth_nucleon` and the quasi-elastic peaks use.

*What it bounds or validates.* The acceptance assumption folded into
`polrad_sigma_qe_u` and `ll_peaks_qe`, and the nucleon-level Rosenbluth inputs of the
quasi-elastic tail.

### 4.5 R15 · the precursor proceedings

> G. I. Gakh, O. Shekhovtsova, *Radiative events in DIS of unpolarized electron by tensor
> polarized deuteron: Radiative corrections*, eConf **C030626** (2003) 351,
> [arXiv:hep-ph/0309123](https://arxiv.org/abs/hep-ph/0309123). Open access; not
> downloaded (superseded by R6). Listed so the dating of the calculation is on the
> record: the JETP paper is a 2004 write-up of a 2003 result, and the zero-citation
> status applies to both.

---

## 5. Tier D — the Monte Carlo codes (what the field does, and what the band must cover)

### 5.1 R18 · RADGEN 1.0 — the engine behind the only tensor-DIS datum

> I. Akushevich, H. Böttcher, D. Ryckbosch, *RADGEN 1.0. Monte Carlo Generator for
> Radiative Events in DIS on Polarized and Unpolarized Targets*,
> [arXiv:hep-ph/9906408](https://arxiv.org/abs/hep-ph/9906408) (1999). Open access; no
> journal reference on the arXiv record. Committed as `refs/hep-ph_9906408.pdf`.

*What it contains.* A MC generator including real radiated photons in DIS on polarized
and unpolarized targets, with analytical and numerical tests. It is POLRAD's formulae
turned into an event generator — the natural comparison object for LiPolGen, which is
also an event generator applying POLRAD's formulae as weights.

*Why the project needs it.* HERMES's b₁ result — the only tensor-DIS measurement in
existence, and the anchor of `RC_DELTA_LOW_X_OPTIMISTIC` — was radiatively corrected
with "a Monte Carlo generator based on RADGEN [9]" (R4, verbatim). Every statement the
tree makes about "HERMES's residual RC systematic" is a statement about **this code's**
output. Reading it is how one learns whether HERMES's 2 × 10⁻³ is comparable to what
`RcModel` reports.

*What it bounds or validates.* `RC_DELTA_LOW_X_OPTIMISTIC = 0.19` (the "as good as
HERMES" edge), and the event-level question of weights-vs-generated-photons that
separates `RcModel`'s weight approach from a radiative event generator.

### 5.2 R19–R22 · HERACLES, DJANGO6/DJANGOH, PEPSI

> A. Kwiatkowski, H. Spiesberger, H.-J. Möhring, *HERACLES: an event generator for ep
> interactions at HERA energies including radiative processes*, Comput. Phys. Commun.
> **69** (1992) 155 (DESY-90-041). DOI
> [10.1016/0010-4655(92)90136-M](https://doi.org/10.1016/0010-4655(92)90136-M). **No
> arXiv, no free copy located** (INSPIRE recid 296225, no fulltext) — use the DOI page.
>
> K. Charchula, G. A. Schuler, H. Spiesberger, *Combined QED and QCD radiative effects in
> deep inelastic lepton–proton scattering: the Monte Carlo generator DJANGO6*, Comput.
> Phys. Commun. **81** (1994) 381 (CERN-TH-7133-94). DOI
> [10.1016/0010-4655(94)90086-8](https://doi.org/10.1016/0010-4655(94)90086-8). **No
> arXiv, no free copy located** (INSPIRE recid 372027) — use the DOI page.
>
> H. Spiesberger, *DJANGOH* (maintained descendant),
> <http://wwwthep.physik.uni-mainz.de/~hspiesb/djangoh/djangoh.html> — free web
> manual; from version 4.6.10 it simulates **longitudinally polarized** NC and CC DIS
> with QED and QCD radiative effects. EIC wiki mirror: <https://wiki.bnl.gov/eic/index.php/DJANGOH>.
>
> L. Mankiewicz, A. Schäfer, M. Veltri, *PEPSI: a Monte Carlo generator for polarized
> leptoproduction*, Comput. Phys. Commun. **71** (1992) 305. DOI
> [10.1016/0010-4655(92)90016-R](https://doi.org/10.1016/0010-4655(92)90016-R). **No
> arXiv, no free copy located** (INSPIRE recid 321673) — use the DOI page.

*What they contain.* HERACLES: single-photon emission from the lepton line, self-energy
and the complete one-loop weak corrections, at **collider** (HERA) kinematics. DJANGO6:
HERACLES ⊗ LEPTO, i.e. QED radiation together with QCD radiation and string
fragmentation. PEPSI: polarized leptoproduction on top of LEPTO/JETSET — the polarized
counterpart the plan pairs with RADGEN.

*Why the project needs them.* `BENCHMARK_PLAN.md` T4 names exactly this trio ("DJANGOH
paired Rad/noRad tables; POLRAD ADGH compiled numbers; PEPSI+RADGEN") and §4 item 7
makes the **DJANGOH published Rad/noRad four-bin table** one of the top-ten items to
wire, at "hours" of effort and "zero dependencies", because it "turns the RC
'must-not-contradict' bound into a Q²-trend comparison". The DJANGO6 paper is where that
table's methodology is defined, and HERACLES is where its QED content is defined.

*What they bound or validate.* The whole of `rc_tail` at **collider** kinematics — the
single most important caveat in `rc.hpp` is that POLRAD is a fixed-target code and the
comment warns in capitals never to port `TPeak` to fixed-target normalisations without
re-deriving the 1/A factors. DJANGOH is the code that does this at HERA/EIC energies,
and it is the only cross-check in the plan that is not another transcription of POLRAD.

### 5.3 R23 · ESFRAD — higher-order QED against POLRAD and RADGEN

> A. Afanasev, I. Akushevich, A. Ilyichev, N. Merenkov, *ESFRAD. FORTRAN code for
> calculation of QED corrections to polarized ep-scattering by the electron structure
> function method*, [arXiv:2212.04730](https://arxiv.org/abs/2212.04730) (2022). Open
> access. Committed as `refs/2212.04730.pdf`.

*What it contains.* The electron structure-function method for **higher-order** QED
effects in polarized DIS, the code ESFRAD, and — the part this project wants — "a
detailed quantitative comparison between the results of ESFRAD and other methods
implemented in the codes POLRAD and RADGEN for calculation of the higher order radiative
corrections".

*Why the project needs it.* `RcModel` is O(α) with a single-emission leading log
(`ll_radiator`) and **no exponentiation**; `rc.hpp` says so explicitly and flags the
uncancelled 1/(1−z). The size of what that leaves out is exactly what this paper
measures, against the two codes the tree already treats as authoritative.

*What it bounds or validates.* `RcTailModel::TPeakPlusLL`'s missing higher orders, and
the implicit assumption that an O(α) band is wide enough to cover them.

### 5.4 R24 · SIDIS-RC EvGen — the maintained modern implementation

> D. Byer, V. Khachatryan, H. Gao, I. Akushevich, A. Ilyichev, C. Peng, A. Prokudin,
> S. Srednyak, Z. Zhao, *SIDIS-RC EvGen: a Monte-Carlo event generator of semi-inclusive
> deep inelastic scattering with the lowest-order QED radiative corrections*, Comput.
> Phys. Commun. **287** (2023) 108702. [arXiv:2210.03785](https://arxiv.org/abs/2210.03785),
> DOI [10.1016/j.cpc.2023.108702](https://doi.org/10.1016/j.cpc.2023.108702). Open
> access. Committed as `refs/2210.03785.pdf`.

*What it contains.* A C++ standalone generator for SIDIS with unpolarized/longitudinally
polarized beam and unpolarized/longitudinally/transversely polarized target, built on
"recent elaborate QED calculations of the lowest-order radiative effects applied to the
leading order Born cross section", with the multi-dimensional integration machinery
described.

*Why the project needs it.* It is the closest existing software analogue of what
LiPolGen's RC sector is: a modern C++ generator that applies the Minsk-school O(α)
formulae event by event. Its structure (how the radiative kernel is integrated, how the
weights are bounded) is a design reference, and it is the paper the task's
"Byer–Afanasev–Kalantarians" pointer resolves to — see §8 for the correction.

*What it bounds or validates.* Nothing numerically for a spin-1 target (it is
spin-1/2 only, and semi-inclusive); it is an implementation reference for
`build_tail_tables` / `tail_ratio_at` and for `RcOptions::tail_max`-style weight bounds.

---

## 6. Tier E — the modern / EIC-era picture

### 6.1 R25 · Afanasev *et al.* (2023) — the topical review

> A. Afanasev, J. C. Bernauer, P. Blunden, J. Blümlein, E. W. Cline, J. M. Friedrich,
> F. Hagelstein, T. Husek, M. Kohl, F. Myhrer, G. Paz, S. Schadmand *et al.*, *Radiative
> Corrections: From Medium to High Energy Experiments*,
> [arXiv:2306.14578](https://arxiv.org/abs/2306.14578) (2023); EPJ A topical review.
> Open access. Committed as `refs/2306.14578.pdf`.

*What it contains.* The state of the field across lepton–proton scattering, QED
corrections in DIS, and radiative light-hadron decays, with emphasis on two-photon
exchange and the associated Monte Carlo codes.

*Why the project needs it.* `04_theory.md` T-37 already named it as "a citable modern
baseline for what '1.5 % on A_zz' means, and the place to find which codes are
maintained". Both matter: the `RC_DELTA_HIGH_X = 0.015` anchor is taken from an
**unpublished** proposal (see R34), and the tree's own audit found no modern successor
to POLRAD in the tensor sector (T-38).

*What it bounds or validates.* `RC_DELTA_HIGH_X` — not by supplying a number, but by
supplying the context in which 1.5 % is or is not a credible modern RC uncertainty.

### 6.2 R26 · Liu, Melnitchouk, Qiu & Sato (2021) — does fixed-target RC port to a collider?

> T. Liu, W. Melnitchouk, J.-W. Qiu, N. Sato, *Factorized approach to radiative
> corrections for inelastic lepton–hadron collisions*, Phys. Rev. D **104** (2021)
> 094033. [arXiv:2008.02895](https://arxiv.org/abs/2008.02895), DOI
> [10.1103/PhysRevD.104.094033](https://doi.org/10.1103/PhysRevD.104.094033). Open
> access. Committed as `refs/2008.02895.pdf`.

*What it contains.* A factorization-based treatment of QED radiation for inclusive and
semi-inclusive DIS, resumming logarithmically enhanced QED radiation into universal
lepton distribution and fragmentation/jet functions, "to systematically account for QED
and QCD radiation contributions to both processes on equal footing". The headline
number: QED effects from the rotational distortion of the hadron transverse momentum —
the mismatch between the experimental and true photon–hadron frames — "can be as large
as 50 % for moderate Q".

*Why the project needs it.* LiPolGen applies a **fixed-target** O(α) formalism at
**collider** kinematics. `rc.hpp` already flags the direction-of-port hazard for the
normalisation; this paper flags a different and larger one — that at a collider the
radiative distortion shows up as a *frame* effect on reconstructed kinematics, not only
as a multiplicative weight. For a generator that hands events to the ePIC reconstruction
chain (where `jacquet-blondel1979-and-bassler-bernardi1995` and the DNN/kinematic-fitting
references in the corpus already live), that is directly on-topic.

*What it bounds or validates.* The **scope** of `RcMode::TensorBand` — specifically
`rc.hpp`'s statement that the band is a weight and "a shift would need a photon-energy
sampler", which this paper is the argument for or against.

### 6.3 R27 · Cammarota, Qiu, Watanabe & Zhang (2025)

> J. Cammarota, J.-W. Qiu, K. Watanabe, J.-Y. Zhang, *Factorized QED and QCD Contribution
> to Deeply Inelastic Scattering*, [arXiv:2505.23487](https://arxiv.org/abs/2505.23487)
> (2025). Open access; no journal reference on the arXiv record as of this pass.
> Committed as `refs/2505.23487.pdf`.

*What it contains.* The first NLO calculation of jointly factorized QED and QCD
short-distance coefficients for inclusive DIS, treating QED radiation from all charged
leptons and quarks on the same footing, with the factorized QED contribution shown to be
infrared safe.

*Why the project needs it.* It is the current state of R26's programme and the frame in
which a future EIC RC treatment will most likely be stated. Cite-only today; relevant if
the project ever moves from a band to a correction.

### 6.4 R28 · CFNS RC whitepaper (2020)

> A. Afanasev, J. Ahmed, I. Akushevich, J. C. Bernauer, P. G. Blunden, A. Bressan,
> D. Byer, E. Cline, M. Diefenthaler, J. M. Friedrich, H. Gao, A. Ilyichev *et al.*,
> *CFNS Ad-Hoc meeting on Radiative Corrections Whitepaper*,
> [arXiv:2012.09970](https://arxiv.org/abs/2012.09970) (2020). Open access. Committed as
> `refs/2012.09970.pdf` — **21.7 MB**, the largest file this domain added; drop it if the
> sibling repo's size policy bites.

*What it contains.* A collection of community contributions written as input to the EIC
Yellow Report process, on the current state of RC technique across nuclear physics.

*Why the project needs it.* It is the document that frames RC as a *systematics
limitation* for the EIC rather than a calculational exercise — the framing this
project's band adopts. Superseded in most respects by R25.

---

## 7. Tier F — HERMES's own treatment, and the two references inside it

### 7.1 R4 · HERMES, *First Measurement of the Tensor Structure Function b₁* (2005)

> HERMES Collaboration, A. Airapetian *et al.*, *First Measurement of the Tensor
> Structure Function b₁ of the Deuteron*, Phys. Rev. Lett. **95** (2005) 242001.
> [arXiv:hep-ex/0506018](https://arxiv.org/abs/hep-ex/0506018), DOI
> [10.1103/PhysRevLett.95.242001](https://doi.org/10.1103/PhysRevLett.95.242001). Open
> access. On disk as `refs/hep-ex_0506018.pdf`, dictionary key `hermes2005-b1-tensor` —
> both added by a sibling domain during this pass. Gap **G2** of
> `references/00_corpus.md`, now closed.

*What its RC paragraph actually says* — read from the PDF this pass, verbatim, because
several tree statements paraphrase it:

> "…due to radiative and detector smearing is treated using an unfolding algorithm,
> which is only sensitive to the detector model, the known unpolarized cross section, and
> the models for the background processes [10]. The radiative background is negligible at
> high x but increases as x → 0 and reaches almost 50 % of the statistics in the lowest-x
> bin. The radiative corrections are calculated using a Monte Carlo generator based on
> RADGEN [9]. The coherent and quasi-elastic radiative tails are estimated using
> parameterizations of the deuteron form factors [11, 12] and corrected for the tracking
> inefficiency due to showering of the radiated photons. The polarized part of the
> quasi-elastic radiative tail is neglected since there is no net tensor effect by
> inclusive scattering on weakly-bound spin-1/2 objects [13]."

and later:

> "The subtraction of the radiative background inflates the size of the statistical and
> the above mentioned systematic uncertainties by almost a factor of 2 at low x. The
> systematic uncertainty of the radiative corrections is ≈ 2 × 10⁻³ for the three bins at
> low x and negligible at high x."

Refs [9] = RADGEN (R18), [10] = HERMES PRD **71** (2005) 012003, [11] = Abbott *et al.*
(R30), [12] = Stein *et al.* (R31), [13] = Zhou *et al.* (R5).

*Why the project needs it.* Three separate loads. (a) It is the **only tensor-DIS datum
in the world** (`BENCHMARK_PLAN.md` T2) and `BENCHMARK_PLAN.md` §4 item 3 wants its
Table II wired as an assertion. (b) `RC_DELTA_LOW_X_OPTIMISTIC = 0.19` is derived in
`rc.hpp` from *this paper's Table II* — 2 × 10⁻³ against a measured |A_zz| = 1.06 × 10⁻²
at ⟨x⟩ = 0.012 — so the constant cannot be re-derived or defended without the PDF.
(c) It is the precedent for neglecting the polarized quasi-elastic tail.

*What it bounds or validates.* `RC_DELTA_LOW_X_OPTIMISTIC`; the `rc_delta(0.063) = 0.1108`
consistency claim in `rc.hpp` ("it does not contradict HERMES's measured 15 % residual
there"); and `qe_tensor_scale`'s default.

### 7.2 R29 · HERMES, the g₁ long paper — the unfolding method in full

> HERMES Collaboration, A. Airapetian *et al.*, *Precise determination of the spin
> structure function g₁ of the proton, deuteron and neutron*, Phys. Rev. D **75** (2007)
> 012007. [arXiv:hep-ex/0609039](https://arxiv.org/abs/hep-ex/0609039), DOI
> [10.1103/PhysRevD.75.012007](https://doi.org/10.1103/PhysRevD.75.012007). Open access.
> On disk as `refs/hep-ex_0609039.pdf`, dictionary key `hermes2007-g1` (sibling domain,
> this pass).

*What it contains* (section headings read from the PDF this pass): a section
**"Unfolding of Radiative and Instrumental Smearing"** plus **Appendix A**, which give
the algorithm the b₁ PRL compresses into one sentence — the migration/smearing matrix
built from a RADGEN-based simulation, the effective subtraction of the radiative
background, the statistical correlation matrix the unfolding produces, and Fig. 6
showing where events from one x-bin actually land. The paper is the "HERMES g₁ paper, in
preparation" cited as ref. [8] of the b₁ PRL.

*Why the project needs it.* This is the answer to the task's "the Airapetian long
paper": there is **no** long companion to the b₁ PRL — an INSPIRE search for HERMES papers with
"tensor" in the title returns exactly one record (recid 684394, the PRL itself) — and
this is the paper that documents the RC/unfolding machinery used for it.

*What it bounds or validates.* Nothing in `RcModel` directly — LiPolGen is a generator,
not an unfolding analysis. It matters for the **comparison**: when a b₁ prediction is
compared against HERMES Table II, this defines what "Born" means in that table, i.e. what
`RcMode::Off` has to correspond to for the comparison to be apples-to-apples.

### 7.3 R30 · Abbott *et al.* (2000) — the deuteron form factors HERMES's tails used

> D. Abbott *et al.* (JLab t₂₀ Collaboration), *Phenomenology of the deuteron
> electromagnetic form factors*, Eur. Phys. J. A **7** (2000) 421.
> [arXiv:nucl-ex/0002003](https://arxiv.org/abs/nucl-ex/0002003), DOI
> [10.1007/PL00013629](https://doi.org/10.1007/PL00013629). Open access. On disk as
> `refs/nucl-ex_0002003.pdf`, dictionary key `abbott2000-ff-parametrization` (sibling
> domain, this pass).

*What it contains.* Closed-form parameterizations (I, II and III, §3 of the paper) of
the world data for the deuteron's G_C, G_M and G_Q, plus the extraction of the node of
G_C. Read from the downloaded PDF this pass.

*Why the project needs it.* Two independent pulls. (a) It is HERMES's ref. [11] — the
form factors behind the coherent and quasi-elastic radiative tails they subtracted, i.e.
the A = 2 analogue of what `HoSpin1FF` supplies for ⁶Li. (b) The tree's own finding in
`benchmarking/02_data_nucleon_deuteron.md`: *"The deuteron form factor in the RC elastic
tail is a harmonic-oscillator shape normalised on the measured moments… The world data
for G_C, G_M, G_Q exist, are parametrized in closed form by Abbott EPJ A 7, and are not
used"*, with the consequence that the deuteron control channel's σ^el is validated only
against POLRAD's own compiled numbers.

*What it bounds or validates.* `HoSpin1FF` on the **deuteron control channel** — the one
place `polrad_sigma_el_u`'s magnitude and shape can be checked against measurement rather
than against POLRAD's `ffdeu`; and hence the T8(a)/T10 gates that `tests/test_rc.cpp`
itself calls "shape-blind".

### 7.4 R31 · Stein *et al.* (1975) — the inclusive radiative-tail recipe

> S. Stein, W. B. Atwood, E. D. Bloom, R. L. A. Cottrell, H. DeStaebler, C. L. Jordan,
> H. G. Piel, C. Y. Prescott, R. Siemann, R. E. Taylor, *Electron scattering at 4° with
> energies of 4.5–20 GeV*, Phys. Rev. D **12** (1975) 1884. DOI
> [10.1103/PhysRevD.12.1884](https://doi.org/10.1103/PhysRevD.12.1884) (APS, paywalled).
> **Free preprint:** SLAC-PUB-1528 (January 1975),
> [INSPIRE recid 100597](https://inspirehep.net/literature/100597), PDF
> [inspirehep.net/files/3b475998b3e06e92047bc4d758c8f0be](https://inspirehep.net/files/3b475998b3e06e92047bc4d758c8f0be).
> Committed as `refs/SLAC-PUB-1528.pdf` (5.6 MB; title page verified).

*What it contains.* The classic SLAC 4° inclusive electron-scattering measurement and,
in its appendices, the working recipe for the **elastic and quasi-elastic radiative
tails** of inclusive electron–nucleus scattering built on Mo–Tsai: how the elastic peak's
tail, the quasi-elastic tail and the inelastic continuum are unfolded in practice, with
the nuclear form-factor and Fermi-smearing inputs spelled out.

*Why the project needs it.* It is HERMES's ref. [12], paired with Abbott, as the source
of the tail parametrizations. It is also the reference implementation of the object
`polrad_sigma_qe_u` + `ll_peaks_qe` + `pauli_suppression` collectively are: a
quasi-elastic tail off a nucleus with Fermi suppression. The tree currently derives that
object from POLRAD Eq. (44) alone and prices its knobs (`qe_kf_gev`, `qe_suppression`) as
"the tail's dominant knobs" with no external recipe to compare against.

*What it bounds or validates.* `polrad_sigma_qe_u`'s magnitude on a nucleus, and the
Pauli-suppression treatment — the piece `rc.hpp` says is "NOT a small choice: the t-peak
reaches down to t_min ~ (x M_N)², so at x ≤ 0.1 most of the QRT integral sits below 2 k_F".

### 7.5 R5 · Zhou *et al.* (1999) — what the neglect actually rests on

> Z.-L. Zhou, M. Bouwhuis, M. Ferro-Luzzi, E. Passchier, R. Alarcon, M. Anghinolfi,
> H. Arenhövel *et al.*, *Tensor Analyzing Powers for Quasi-Elastic Electron Scattering
> from Deuterium*, Phys. Rev. Lett. **82** (1999) 687.
> [arXiv:nucl-ex/9809002](https://arxiv.org/abs/nucl-ex/9809002), DOI
> [10.1103/PhysRevLett.82.687](https://doi.org/10.1103/PhysRevLett.82.687). Open access,
> 19 INSPIRE citations. Committed as `refs/nucl-ex_9809002.pdf`. Gap **G7**.

*What it contains.* Abstract, verbatim: "We report on a first measurement of tensor
analyzing powers in quasi-elastic electron–deuteron scattering at an average
three-momentum transfer of 1.7 fm⁻¹. Data sensitive to the spin-dependent nucleon
density in the deuteron were obtained for missing momenta up to 150 MeV/c with a tensor
polarized ²H target internal to an electron storage ring. The data are well described by
a calculation that includes the effects of final-state interaction, meson-exchange and
isobar currents, and leading-order relativistic contributions."

*Why the project needs it.* `RcOptions::qe_tensor_scale` defaults to 0.0 — the polarized
quasi-elastic tail is set to zero — and the justification chain is HERMES's sentence
("no net tensor effect by inclusive scattering on weakly-bound spin-1/2 objects [13]")
→ this paper. `phase_B_numbers.md` §B3 already withdrew the A = 2 → A = 6 transfer of
that argument and recorded that the quasi-elastic term is 22 % / 73 % / 99.9 % of the
tail at x = 0.01 / 0.10 / 0.30, i.e. "the largest unpriced piece". The paper is now in
the corpus so the argument can be read rather than inherited.

*What it bounds or validates.* `qe_tensor_scale`'s default, and — read carefully — the
scope of HERMES's claim: this is an **exclusive-kinematics coincidence** measurement at
low momentum transfer with a **non-zero** measured tensor analyzing power, described by a
calculation full of FSI/MEC/IC. Whatever supports "no net tensor effect in *inclusive*
scattering" must be an integral statement over that measurement, not the measurement
itself, and the integral has not been done for A = 6.

### 7.6 R32, R33 · the two knobs' sources

> E. J. Moniz, I. Sick, R. R. Whitney, J. R. Ficenec, R. D. Kephart, W. P. Trower,
> *Nuclear Fermi momenta from quasielastic electron scattering*, Phys. Rev. Lett. **26**
> (1971) 445. DOI [10.1103/PhysRevLett.26.445](https://doi.org/10.1103/PhysRevLett.26.445)
> (APS, paywalled; no free copy located — INSPIRE recid 67648 carries no fulltext).
>
> T. de Forest Jr., J. D. Walecka, *Electron scattering and nuclear structure*, Adv. Phys.
> **15** (1966) 1. DOI
> [10.1080/00018736600101254](https://doi.org/10.1080/00018736600101254) (Taylor &
> Francis, paywalled; no free copy located — INSPIRE recid 50061).

Moniz is the source of `RC_QE_KF_GEV = 0.169` (⁶Li), which `rc.hpp` calls not a small
choice; de Forest–Walecka is the source of `pauli_suppression`'s
S(q) = (3/4)u − u³/16. Both are named in `rc.hpp` and in `PHYSICS_CHANNELS.md`'s
reference list (`[Moniz71]`, `[dFW66]`) with no copy in the repository. Neither is
obtainable open-access; the DOI pages above are the links to use. Priority is **M**, not
H: the constants are standard and quoted correctly, but the k_F value for ⁶Li
specifically should be read off Moniz's own table before the next time it is quoted with
a digit-count that implies a fit.

### 7.7 R34 · Poudel *et al.* (2025) — and what it does **not** contain

> J. Poudel, A. Bacchetta, J.-P. Chen, N. Santiesteban, *Experimental Study of Tensor
> Structure Function of Deuteron*, Eur. Phys. J. A **61** (2025) 81.
> [arXiv:2506.04506](https://arxiv.org/abs/2506.04506), DOI
> [10.1140/epja/s10050-025-01558-w](https://doi.org/10.1140/epja/s10050-025-01558-w).
> Open access. On disk as `refs/2506.04506.pdf`, dictionary key
> `poudel2025-tensor-sf-review` (sibling domain, this pass).

`RC_X_HIGH = 0.16` is taken from this paper (p. 8: 0.16 < x < 0.49, 0.8 < Q² < 5.0 GeV²).
`RC_DELTA_HIGH_X = 0.015` is **not**: `rc.hpp` records that the 1.5 % appears only in the
unpublished E12-13-011 proposal and that "the published companion arXiv:2506.04506
contains no radiative-correction discussion at all", leaving the action "cite it by page
for the 1.5 %, or DROP the anchor". This pass did not re-open that question; it is
recorded here so the two constants' provenances are not conflated in any future citation.

---

## 8. What could not be obtained, and one correction to the brief

**No free copy exists (paywalled; DOI link is the one to use):**

| Reference | Where to buy/borrow | Why it still matters |
|---|---|---|
| Mo–Tsai, RMP **41** (1969) 205 | <https://doi.org/10.1103/RevModPhys.41.205> | the published equation numbering; the preprint (R1) is otherwise equivalent |
| Tsai, PR **122** (1961) 1898 | <https://doi.org/10.1103/PhysRev.122.1898> | historical only; superseded by R1/R2 |
| Akushevich & Shumeiko, JPG **20** (1994) 513 | <https://doi.org/10.1088/0954-3899/20/4/001> | **the A-dependence of polarized RC** — the highest-value paywalled item in this domain |
| Bardin & Shumeiko, NPB **127** (1977) 242 | <https://doi.org/10.1016/0550-3213(77)90213-9> | the free JINR preprint (R7) is **Russian-language**; the English text is here |
| HERACLES, CPC **69** (1992) 155 | <https://doi.org/10.1016/0010-4655(92)90136-M> | the QED content of the DJANGOH comparison |
| DJANGO6, CPC **81** (1994) 381 | <https://doi.org/10.1016/0010-4655(94)90086-8> | the Rad/noRad methodology for T4 item 7 |
| PEPSI, CPC **71** (1992) 305 | <https://doi.org/10.1016/0010-4655(92)90016-R> | the polarized MC in the T4 pairing |
| Ent, Filippone, Makins, Milner, O'Neill, Wasson, *Radiative corrections for (e,e′p) reactions at GeV energies*, PRC **64** (2001) 054610 (JLAB-PHY-00-07) | <https://doi.org/10.1103/PhysRevC.64.054610> | the **coincidence** RC practice that the tagged `rc_tail ≡ 1` should be argued against, alongside R16. JLab record: <https://www1.jlab.org/Ul/publications/view_pub.cfm?pub_id=1090> |
| Kuraev & Fadin, Sov. J. Nucl. Phys. **41** (1985) 466 | INSPIRE recid 217313 — no DOI, no e-print | the exponentiated D(z) `ll_radiator` omits |
| Moniz *et al.*, PRL **26** (1971) 445 | <https://doi.org/10.1103/PhysRevLett.26.445> | `RC_QE_KF_GEV` |
| de Forest & Walecka, Adv. Phys. **15** (1966) 1 | <https://doi.org/10.1080/00018736600101254> | `pauli_suppression` |

**Free but not fetched (deliberate):** Tsai, RMP **46** (1974) 815 = SLAC-PUB-1365,
6.7 MB, at <https://inspirehep.net/files/be2ffc45ed59b1b28208cc386c1e16b3> — the
equivalent-radiator and external-bremsstrahlung reference. External bremsstrahlung does
not exist at a collider, so only the radiator half applies and R13/R23 cover it better.
Fetch it if the `ll_radiator` citation is ever written in full. Also
`hep-ph/0309123` (R15), superseded by R6.

**Dead link, recorded:** the old CPC Program Library summary host
`cpc.cs.qub.ac.uk/summaries/ADGH_v1_0.html` does not resolve (checked this pass); the
Elsevier CPC Program Library has moved to Mendeley Data. The catalogue identifier
`ADGH_v1_0` is still correct and `polrad_transcription_check.md` records that the tree
already holds the FORTRAN (decks `apptai`, `elu`, `elp`, `elq`; steering `polrad20.cra`),
which is the authoritative artifact anyway.

**One correction to the brief.** The task named *"Byer–Afanasev–Kalantarians"*. No such
paper exists: searches on INSPIRE (author `Byer, Duane`, all 4 records; author
`Kalantarians` + title `radiative`) and on the arXiv API (`au:Byer AND au:Afanasev`)
return no joint Byer–Kalantarians radiative-corrections paper. The two real referents are
**R24** (Byer *et al.*, SIDIS-RC EvGen, `2210.03785`) and **R28** (the CFNS whitepaper
`2012.09970`, which has both Afanasev and Byer among its authors). N. Kalantarians'
radiative-corrections work appears in that community but not in a paper with those three
names; if a specific result was meant, it needs a title.

**Also checked and reported as absent:** an INSPIRE title search (HERMES + "tensor")
returns only the b₁ PRL itself (recid 684394), i.e. no long companion paper was found; no NLO or modern successor
to POLRAD in the **tensor** sector was found, confirming `04_theory.md` T-38's finding;
and no radiative-correction calculation of any kind exists for the φ-dependent tensor
observables (cos φ_TL, cos 2φ_TT) or for spin-3/2, confirming
`open_items/physics_literature.md` §3's verdict, which this pass re-searched and did not
overturn.

---

## 9. What to do with these, in order

| # | Action | References | Effort | What it closes |
|---|---|---|---|---|
| 1 | Write the `[MT69]` gate in `tests/test_rc.cpp` against SLAC-PUB-380 Eq. (B.5) — σ^el normalisation for a point target, then the deuteron control channel | R1, R2 | 1–2 d | `benchmarking/00` §9 item 5; `SUMMARY.md` item 5; `04_theory.md` T-35/T-36 |
| 2 | Wire the DJANGOH published Rad/noRad four-bin table as a Q²-trend comparison | R20, R21, R19 | hours | `BENCHMARK_PLAN.md` §4 item 7 (top-ten, zero dependencies) |
| 3 | Wire HERMES Table II as an assertion **and** re-derive `RC_DELTA_LOW_X_OPTIMISTIC` from the PDF rather than from the comment | R4 | ½ d | `BENCHMARK_PLAN.md` §4 item 3 |
| 4 | Read R8 and put a **sourced** number on `a_transfer_frac` — the A = 2 → A = 6 transfer of δ(x) | R8 (paywalled: needs library access), R2 | 2–3 d + access | registry row B6 / `AUTHOR_DECISIONS.md` §B14 |
| 5 | Wire Abbott's G_C/G_M/G_Q into the deuteron control channel's `Spin1ElasticFF` and un-blind the T8(a)/T10 shape gates | R30 | 1–2 d | `02_data_nucleon_deuteron.md`'s "exist … and are not used" |
| 6 | Map R14's tensor elastic-tail observables onto POLRAD Eq. (A.4) and get one independent number for `polrad_im_el_spin1`'s `t` component | R14, R3 | 2–4 d | the only single-source-transcription risk left in the tail |
| 7 | Price the tagged `rc_tail ≡ 1` against the coincidence RC literature instead of a paragraph | R16, Ent *et al.* (paywalled) | 2–3 d | `phase_B_numbers.md` §B4 |
| 8 | Quantify the missing exponentiation in `TPeakPlusLL` from ESFRAD's own POLRAD/RADGEN comparison; cite Kuraev–Fadin in `ll_radiator` | R23, R13 | 1 d | `rc.hpp`'s uncancelled 1/(1−z) caveat |
| 9 | Decide whether a collider RC needs a photon-energy **shift** and not only a weight | R26, R27, R25 | reading | `rc.hpp`'s "a shift would need a photon-energy sampler" |
| 10 | Re-read R5 before the next time `qe_tensor_scale = 0` is defended in prose | R5, R4 | hours | `phase_B_numbers.md` §B3's withdrawn transfer clause |

---

*Provenance of this file: written 2026-09-16 by the reference-research stage (radiative
domain). Every arXiv id, DOI, INSPIRE recid, byte count and HTTP status quoted above was
obtained in this pass; every verbatim quotation was read from the downloaded PDF or from
the arXiv/INSPIRE abstract record named beside it. Claims about LiPolGen's own code and
documents are cited to the file that makes them and were read from the tree, not
inferred.*
