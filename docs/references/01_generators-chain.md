# 01 — Generators, the software chain, and the EIC/ePIC documents

**Domain.** The external *codes* LiPolGen leans on or is benchmarked
against, the *event-record and software chain* it writes into, and the
*parton-distribution and R fits* its structure-function backends read.
Companion to `00_corpus.md` (what the project already holds) and to the
benchmarking series `benchmarking/01_generators.md` (what each generator
could and could not check) and `benchmarking/05_chain.md` (the chain).

**Search date 2026-09-16.** Every entry below was verified in this pass by
fetching the arXiv API record (`export.arxiv.org/api/query`) or the
INSPIRE literature record (`inspirehep.net/api/literature`) and reading
the returned title, author list, year, journal reference, DOI and
abstract. Where a claim in the tree is *corrected* or *confirmed* by that
fetch, it is said so in place. Anything not confirmed is marked
**UNVERIFIED** with the reason; nothing is listed silently.

**Downloads.** 24 open-access files were added to the shared corpus at
`/home/cpeng/Projects/polli/PolarizedLithiumSim/refs/` in this pass
(§5). No existing file in either repository was modified, and
`refs_dict.json` was **not** touched — appending the dictionary entries is
left to the final stage, because several domain passes were running
concurrently against the same file.

---

## 1. Table

Columns: **tier** maps onto `BENCHMARK_PLAN.md` §2 (T1 nucleon, T2
deuteron/b₁, T3 nuclear inputs, T4 generator-vs-generator, T5 MC
closures, T6 chain); **serves** names the LiPolGen symbol, backend,
bridge or benchmark row it is needed for; **OA** = open access, and
"local" = a PDF is now on disk under `PolarizedLithiumSim/refs/`.

### 1a. The PYTHIA backend and the LHAup bridge

| # | reference | id | OA / local | tier | serves |
|---|---|---|---|---|---|
| GC-1 | Bierlich *et al.*, **PYTHIA 8.3** manual, SciPost Phys. Codebases **8** (2022) | [2203.11601](https://arxiv.org/abs/2203.11601), doi:[10.21468/SciPostPhysCodeb.8](https://doi.org/10.21468/SciPostPhysCodeb.8) | yes / `2203.11601.pdf` | T4 | the hadronizer itself; `src/pythia/`, `docs/PYTHIA_BRIDGE.md`, gate `D-1` |
| GC-2 | Alwall *et al.*, **A standard format for Les Houches Event Files**, CPC **176** (2007) 300 | [hep-ph/0609017](https://arxiv.org/abs/hep-ph/0609017), doi:[10.1016/j.cpc.2006.11.010](https://doi.org/10.1016/j.cpc.2006.11.010) | yes / `hep-ph_0609017.pdf` | T4 | the LHEF record `LhaupDis::setEvent` writes (`PYTHIA_BRIDGE.md` §6) |
| GC-3 | Boos *et al.*, **Generic User Process Interface for Event Generators** (Les Houches Accord 1), hep-ph/0109068 | [hep-ph/0109068](https://arxiv.org/abs/hep-ph/0109068) | yes / `hep-ph_0109068.pdf` | T4 | defines `HEPRUP`/`HEPEUP` and the **strategy-3** weighting the bridge uses |
| GC-4 | Cabouat & Sjöstrand, **Some Dipole Shower Studies**, EPJ C **78** (2018) 226 | [1710.00391](https://arxiv.org/abs/1710.00391), doi:[10.1140/epjc/s10052-018-5645-z](https://doi.org/10.1140/epjc/s10052-018-5645-z) | yes / `1710.00391.pdf` | T4 | the *DIS-specific* pedigree of the shower the bridge inherits (dipole recoil) |
| GC-5 | Helenius, Laulainen & Preuss, **Multi-Jet Production in DIS with Pythia**, JHEP **05** (2025) 153 | [2410.20950](https://arxiv.org/abs/2410.20950), doi:[10.1007/JHEP05(2025)153](https://doi.org/10.1007/JHEP05(2025)153) | yes / `2410.20950.pdf` | T4 | PYTHIA's DIS shower validated against **H1** — the published pedigree `D-1` borrows |
| GC-6 | Helenius, Laulainen & Preuss, **Multiplicative matching of NC DIS at NLO in PYTHIA 8** | [2605.00502](https://arxiv.org/abs/2605.00502) | yes / not fetched | T4 | the NLO successor of GC-5; **no journal ref yet** (arXiv record carries none) |
| GC-7 | H1 Collaboration, **Measurement and QCD Analysis of Diffractive DIS at HERA** (H1 2006 DPDF Fit A/B), EPJ C **48** (2006) 715 | [hep-ex/0606004](https://arxiv.org/abs/hep-ex/0606004), doi:[10.1140/epjc/s10052-006-0035-3](https://doi.org/10.1140/epjc/s10052-006-0035-3) | yes / `hep-ex_0606004.pdf` | T4 | the DPDF family behind `PDF:PomSet = 6`, i.e. the **T2 Pomeron tier** (`F-8`) |
| GC-8 | Bierlich, Gustafson, Lönnblad & Shah, **The Angantyr model for Heavy-Ion Collisions in PYTHIA8**, JHEP **10** (2018) 134 | [1806.10820](https://arxiv.org/abs/1806.10820), doi:[10.1007/JHEP10(2018)134](https://doi.org/10.1007/JHEP10(2018)134) | yes / `1806.10820.pdf` | T4 | documents *why* a nucleus id must never reach `Beams:idA/idB` — the negative result of `benchmarking/01` §2.2 |

### 1b. Peer DIS / eA / polarized generators

| # | reference | id | OA / local | tier | serves |
|---|---|---|---|---|---|
| GC-9 | Chang *et al.*, **BeAGLE**, PRD **106** (2022) 012007 | [2204.11998](https://arxiv.org/abs/2204.11998), doi:[10.1103/PhysRevD.106.012007](https://doi.org/10.1103/PhysRevD.106.012007) | yes / `2204.11998.pdf` | T4 | the eA reference generator; the eD tagged-spectator and far-forward comparison (`E-16`, `CH-8`) |
| GC-10 | Charchula, Schuler & Spiesberger, **DJANGO6**, CPC **81** (1994) 381 | doi:[10.1016/0010-4655(94)90086-8](https://doi.org/10.1016/0010-4655(94)90086-8) | **no** / no free copy | T4 | the code paper of DJANGOH, the independent one-loop RC engine of `D-3`/`C-8` |
| GC-11 | Kwiatkowski, Spiesberger & Möhring, **HERACLES**, CPC **69** (1992) 155 | doi:[10.1016/0010-4655(92)90136-M](https://doi.org/10.1016/0010-4655(92)90136-M) | **no** / no free copy | T4 | the RC engine *inside* DJANGOH — the actual second opinion on the radiative correction |
| GC-12 | Ingelman, Edin & Rathsman, **LEPTO 6.5**, CPC **101** (1997) 108 | [hep-ph/9605286](https://arxiv.org/abs/hep-ph/9605286), doi:[10.1016/S0010-4655(96)00157-9](https://doi.org/10.1016/S0010-4655(96)00157-9) | yes / `hep-ph_9605286.pdf` | T4 | the DIS hard-process + fragmentation layer under DJANGOH **and** PEPSI |
| GC-13 | Mankiewicz, Schäfer & Veltri, **PEPSI: a MC generator for polarized leptoproduction**, CPC **71** (1992) 305 | doi:[10.1016/0010-4655(92)90016-R](https://doi.org/10.1016/0010-4655(92)90016-R) | **no** / no free copy | T4 | the only code with a live tensor (`pnrun = ±2`) branch — the RC *kernel* cross-check |
| GC-14 | Akushevich, Böttcher & Ryckbosch, **RADGEN 1.0**, hep-ph/9906408 | [hep-ph/9906408](https://arxiv.org/abs/hep-ph/9906408) | yes / `hep-ph_9906408.pdf` | T4 | the POLRAD-derived RC generator shipped inside both BeAGLE and PEPSI (`CH-22`) |
| GC-15 | Jung, **Hard diffractive scattering … and the MC generator RAPGAP**, CPC **86** (1995) 147 | doi:[10.1016/0010-4655(94)00150-Z](https://doi.org/10.1016/0010-4655(94)00150-Z); preprint DESY 93-182 | yes (preprint) / `jung1995_rapgap_desy93-182.pdf` | T4 | the DDIS generator behind the ePIC `rapgap3.310` samples — the T2 tier's final-state analogue (`CH-9`) |
| GC-16 | **STEG** / `JeffersonLab/LightIonEIC` — the JLab tagged-deuteron generator | repository only | code / — | T4 | `E-1` (Cosyn–Weiss TABLE II) made executable. **No paper exists** — see §3.2 |
| GC-17 | **TOPEG** (Orsay-Perugia Event Generator) — coherent DVCS off ⁴He | no code paper; physics in GC-18/GC-19 | — | T4 | **corrects the brief**: TOPEG is *not* a tagged-deuteron generator (§3.3) |
| GC-18 | Fucini, Scopetta & Viviani, **Coherent DVCS off ⁴He**, PRC **98** (2018) 015203 | [1805.05877](https://arxiv.org/abs/1805.05877), doi:[10.1103/PhysRevC.98.015203](https://doi.org/10.1103/PhysRevC.98.015203) | yes / `1805.05877.pdf` | T4 | the IA convolution TOPEG implements — the citable primary for "TOPEG" |
| GC-19 | Fucini, Scopetta & Viviani, **Incoherent DVCS off ⁴He**, PRC **102** (2020) 065205 | [2008.11437](https://arxiv.org/abs/2008.11437), doi:[10.1103/PhysRevC.102.065205](https://doi.org/10.1103/PhysRevC.102.065205) | yes / `2008.11437.pdf` | T4 | the paper that names the Orsay-Perugia generator in text |
| GC-20 | **lAger** 3.7.0, `eicweb.phy.anl.gov/monte_carlo/lager`, GPL-3.0-or-later | repository only | code / — | T4 | `fermi87(P, A)` as an independent n(k) shape check; the ⁴He J/ψ grid. **No paper** (§3.2) |
| GC-21 | Lomnitz & Klein, **eSTARlight**, PRC **99** (2019) 015203 | [1803.06420](https://arxiv.org/abs/1803.06420), doi:[10.1103/PhysRevC.99.015203](https://doi.org/10.1103/PhysRevC.99.015203) | yes / `1803.06420.pdf` | T4 | the in-tree coherent ⁶Li baseline `D-4` |
| GC-22 | Klein *et al.*, **STARlight**, CPC **212** (2017) 258 | [1607.03838](https://arxiv.org/abs/1607.03838), doi:[10.1016/j.cpc.2016.10.016](https://doi.org/10.1016/j.cpc.2016.10.016) | yes / `1607.03838.pdf` | T4 | the Q²→0 limit of `D-4` |
| GC-23 | Toll & Ullrich, **Sartre 1**, CPC **185** (2014) 1835 | [1307.8059](https://arxiv.org/abs/1307.8059), doi:[10.1016/j.cpc.2014.03.010](https://doi.org/10.1016/j.cpc.2014.03.010) | yes / `1307.8059.pdf` | T4 | the *negative* result: §3.3.4 lists the supported nuclei and A = 6 is not among them |

### 1c. The event record and the ePIC chain

| # | reference | id | OA / local | tier | serves |
|---|---|---|---|---|---|
| GC-24 | Buckley *et al.*, **The HepMC3 Event Record Library**, CPC **260** (2021) 107310 | [1912.08005](https://arxiv.org/abs/1912.08005), doi:[10.1016/j.cpc.2020.107310](https://doi.org/10.1016/j.cpc.2020.107310) | yes / `1912.08005.pdf` | T6 | the output format itself — `docs/HEPMC3_CONVENTION.md`, gates `F-1`…`F-4` |
| GC-25 | Gardiner, Isaacson & Pickering, **NuHepMC**, SciPost Phys. Codebases **57** (2025) | [2310.13211](https://arxiv.org/abs/2310.13211), doi:[10.21468/SciPostPhysCodeb.57](https://doi.org/10.21468/SciPostPhysCodeb.57) | yes / `2310.13211.pdf` | T6 | the **template** for the ion-spin attribute convention (`F-4`, `CH-18`) |
| GC-26 | Bierlich *et al.*, **Rivet version 4 release note**, SciPost Phys. Codebases **36** (2024) | [2404.15984](https://arxiv.org/abs/2404.15984), doi:[10.21468/SciPostPhysCodeb.36](https://doi.org/10.21468/SciPostPhysCodeb.36) | yes / `2404.15984.pdf` | T6 | Rivet 4.1.2 in `eic_xl-nightly.sif`; `CH-1`…`CH-5`, the 46 H1/ZEUS DIS analyses |
| GC-27 | **eic-smear**, `github.com/eic/eic-smear`, GPL-3.0 | repository only | code / — | T6 | the generator-text ↔ ROOT ↔ HepMC3 bridge every other EIC generator crosses (`CH-23`). **No paper** (§3.2) |
| GC-28 | **afterburner / `abconv`**, `github.com/eic/afterburner`, no license declared | repository only | code / — | T6 | crossing angle + beam effects; gate `F-2`, and `CH-12` (no light-ion preset) |
| GC-29 | Frank, Gaede, Grefe & Mato, **DD4hep**, J. Phys. Conf. Ser. **513** (2014) 022010 | doi:[10.1088/1742-6596/513/2/022010](https://doi.org/10.1088/1742-6596/513/2/022010) | yes (IOP OA) / not fetched | T6 | the geometry layer under `npsim`/`eic/epic` (`F-1`, `CH-13`) |
| GC-30 | Agostinelli *et al.*, **GEANT4 — a simulation toolkit**, NIM A **506** (2003) 250 | doi:[10.1016/S0168-9002(03)01368-8](https://doi.org/10.1016/S0168-9002(03)01368-8) | yes (Elsevier OA) / not fetched | T6 | the transport engine `npsim` runs (`F-1`) |
| GC-31 | Particle Data Group (S. Navas *et al.*), **Review of Particle Physics**, PRD **110** (2024) 030001 — §"Monte Carlo particle numbering scheme" | doi:[10.1103/PhysRevD.110.030001](https://doi.org/10.1103/PhysRevD.110.030001) | yes | T6 | the `10LZZZAAAI` ion codes the tree emits for ⁶Li/⁷Li/α/d/t (`CH-19`, `05_chain.md` §5.1) |

### 1d. PDFs, nuclear PDFs and R — the structure-function backends

| # | reference | id | OA / local | tier | serves |
|---|---|---|---|---|---|
| GC-32 | Buckley *et al.*, **LHAPDF6**, EPJ C **75** (2015) 132 | [1412.7420](https://arxiv.org/abs/1412.7420), doi:[10.1140/epjc/s10052-015-3318-8](https://doi.org/10.1140/epjc/s10052-015-3318-8) | yes / `1412.7420.pdf` | T1 | the grid library behind `src/lhapdf/lhapdf_sf.cpp` (the whole optional LHAPDF tier) |
| GC-33 | Hou *et al.*, **CT18**, PRD **103** (2021) 014013 | [1912.10053](https://arxiv.org/abs/1912.10053), doi:[10.1103/PhysRevD.103.014013](https://doi.org/10.1103/PhysRevD.103.014013) | yes / `1912.10053.pdf` | T1 | `LhapdfSF`'s **default `CT18NLO`**, and the `CT18ANLO` proton baseline of EPPS21 |
| GC-34 | Martin, Stirling, Thorne & Watt, **MSTW 2008**, EPJ C **63** (2009) 189 | [0901.0002](https://arxiv.org/abs/0901.0002), doi:[10.1140/epjc/s10052-009-1072-5](https://doi.org/10.1140/epjc/s10052-009-1072-5) | yes / `0901.0002.pdf` | T2 | `MstwSF` — **CDKS's own** unpolarized PDF, and therefore part of the b₁/CDKS `A = 2` gate's definition |
| GC-35 | Nocera, Ball, Forte, Ridolfi & Rojo, **NNPDFpol1.1**, NPB **887** (2014) 276 | [1406.5539](https://arxiv.org/abs/1406.5539), doi:[10.1016/j.nuclphysb.2014.08.008](https://doi.org/10.1016/j.nuclphysb.2014.08.008) | yes / `1406.5539.pdf` | T1 | `LhapdfG1`'s default `NNPDFpol11_100`, including its Δb = 0 convention |
| GC-36 | Eskola, Paakkinen, Paukkunen & Salgado, **EPPS21**, EPJ C **82** (2022) 413 | [2112.12462](https://arxiv.org/abs/2112.12462), doi:[10.1140/epjc/s10052-022-10359-0](https://doi.org/10.1140/epjc/s10052-022-10359-0) | yes / `2112.12462.pdf` | T3 | `Epps21Ratio` / `EmcBaseline::Epps21` — the unpolarized ⁶Li EMC baseline (`D-5`) |
| GC-37 | Kovarik *et al.*, **nCTEQ15**, PRD **93** (2016) 085037 | [1509.00792](https://arxiv.org/abs/1509.00792), doi:[10.1103/PhysRevD.93.085037](https://doi.org/10.1103/PhysRevD.93.085037) | yes / `1509.00792.pdf` | T3 | `T-39` four-fit EMC spread; and `nCTEQ15_7_3` is a real **⁷Li** grid (`T-40`) |
| GC-38 | Abdul Khalek *et al.*, **nNNPDF3.0**, EPJ C **82** (2022) 507 | [2201.12363](https://arxiv.org/abs/2201.12363), doi:[10.1140/epjc/s10052-022-10417-7](https://doi.org/10.1140/epjc/s10052-022-10417-7) | yes / `2201.12363.pdf` | T3 | the second fit in `T-39`; `nNNPDF30_nlo_as_0118_A6_Z3` |
| GC-39 | Walt, Helenius & Vogelsang, **TUJU19**, PRD **100** (2019) 096015 | [1908.03355](https://arxiv.org/abs/1908.03355), doi:[10.1103/PhysRevD.100.096015](https://doi.org/10.1103/PhysRevD.100.096015) | yes / `1908.03355.pdf` | T3 | the source of the `TUJU19_*_7_3` ⁷Li grids named in `T-40` |
| GC-40 | Whitlow, Rock, Bodek, Riordan & Dasu, **R1990**, PLB **250** (1990) 193 | doi:[10.1016/0370-2693(90)91176-C](https://doi.org/10.1016/0370-2693(90)91176-C); SLAC-PUB-5284 | yes (preprint) / `whitlow1990_R.pdf` | T1 | **closes a named gap** — the source of `sf.hpp:43`'s "R1990-like" toy and of HERMES's ref [18] |
| GC-41 | Abe *et al.* (E143), **R1998 world fit**, PLB **452** (1999) 194 | [hep-ex/9808028](https://arxiv.org/abs/hep-ex/9808028) | yes / `hep-ex_9808028.pdf` | T1 | `r1998()` in `src/core/sf.cpp` — the default R. **ALREADY IN CORPUS** |

### 1e. EIC / ePIC programme documents — all already in the corpus

| # | reference | id | local | tier | serves |
|---|---|---|---|---|---|
| GC-42 | Abdul Khalek *et al.*, **EIC Yellow Report**, NPA **1026** (2022) 122447 | [2103.05419](https://arxiv.org/abs/2103.05419) | `2103.05419_part{1..4}.pdf` | T6 | generator verification §8.1.6; TOPEG/Sartre/lAger provenance; far-forward acceptance. **IN CORPUS** |
| GC-43 | Willeke *et al.*, **EIC Conceptual Design Report** (2021) | doi:[10.2172/1765663](https://doi.org/10.2172/1765663) | none (178 MB) | T6 | machine parameters, IR layout. **IN CORPUS** (`file: null`) |
| GC-44 | ePIC Collaboration, **Preliminary Design Report** (2024) | doi:[10.5281/zenodo.13866213](https://doi.org/10.5281/zenodo.13866213) | none (114 MB) | T6 | detector status. **IN CORPUS** (`file: null`) |
| GC-45 | ePIC Collaboration, **Preliminary TDR v3.1** (2026) | doi:[10.5281/zenodo.18271601](https://doi.org/10.5281/zenodo.18271601) | none (260 MB) | T6 | the current detector document. **IN CORPUS** (`file: null`) |
| GC-46 | Jentsch, **Far-Forward Detectors and Physics** (DIS 2023) | — | `ePIC_far_forward_talk_DIS_2023_v2.pdf` | T6 | far-forward acceptance for α/d/t spectators. **IN CORPUS** |
| GC-47 | Chang *et al.*, **Opportunities for imaging light nuclei at IR-8** | [2511.05638](https://arxiv.org/abs/2511.05638) | `2511.05638.pdf` | T6 | the IR-8 secondary focus; the `E-16` efficiency table. **IN CORPUS** |
| GC-48 | Gamage *et al.*, **Second interaction region** | — | on disk | T6 | the IR-8 optics. **IN CORPUS** |

---

## 2. What each reference contains, and why the project needs it

### 2.1 The PYTHIA backend and the LHAup bridge

**GC-1 — PYTHIA 8.3 (2203.11601).** The 600-page physics-and-usage manual
for the version LiPolGen links (8.317, built at
`/home/cpeng/Projects/polli/deps/install`). Everything downstream of
LiPolGen's hard process — the spacelike/timelike showers, MPI, the Lund
string, the decay tables — is specified here, so this is the reference
that stands behind gate `D-1` and every hadronic-final-state number the
generator prints. A verification note worth recording: the **arXiv API
record carries no `journal_ref` and no DOI**, and `benchmarking/01`
§2.1 therefore cites it as an arXiv item. That is now out of date —
**INSPIRE recid 2056998 gives SciPost Phys. Codebases 2022 (2022) 8,
doi:10.21468/SciPostPhysCodeb.8**, 1778 citations. Cite the SciPost
Codebases item.

**GC-2 / GC-3 — the Les Houches accords.** Two references, and the
project needs both, because the bridge uses both halves. GC-3
(hep-ph/0109068, 18 authors, Les Houches 2001) defines the `HEPRUP` /
`HEPEUP` Fortran common blocks and, with them, the **weighting
strategies** — including strategy 3 (unit weights, process choice made
by the parton-level generator), which is exactly what
`src/pythia/lhaup_dis.hpp` declares and why `info.sigmaGen()` is not
physical in the bridge (the "does not validate" column of `D-1`). GC-2
(CPC 176 (2007) 300) is the XML re-expression of the *same* information
content — its abstract says so explicitly ("The information content is
identical with what was already defined by the Les Houches Accord five
years ago, but then in terms of Fortran commonblocks") — and it is the
document that fixes the four-line-per-process `<event>` layout that
`LhaupDis::setEvent` writes, described in `PYTHIA_BRIDGE.md` §6. Neither
was in the corpus; both are short and open access.

**GC-4 — dipole recoil (1710.00391).** The tree's claim that PYTHIA can
do DIS at all rests on this paper: the initial-state dipole-recoil
scheme is what `SpacelikeShowers.xml` says "for the first time allows the
simulation of Deeply Inelastic Scattering processes in PYTHIA 8". The
abstract's own framing is a comparison of dipole and global-recoil
showers, with DIS as the test case. This is the pedigree LiPolGen
inherits **for the shower only** — it says nothing about a polarized or
nuclear beam.

**GC-5 / GC-6 — the modern DIS validation (2410.20950, 2605.00502).**
GC-5 is now published (JHEP 05 (2025) 153) and introduces multi-jet
merging for DIS in Vincia, compared against H1; GC-6 is the NLO
multiplicative-matching successor and, as of this pass, **carries no
journal reference and no DOI on its arXiv record**. Together they are the
only recent, explicit "PYTHIA 8 vs HERA DIS data" statements in the
literature, and they are what the `CH-1` Rivet plan would reproduce
locally. Value: they turn `D-1`'s one-observable ±20 % window into a
citable external pedigree for the shower half of the bridge.

**GC-7 — H1 2006 DPDF (hep-ex/0606004).** PYTHIA's `PDF:PomSet = 6`
*is* this fit. The T2 coherent tier's Pomeron flux and diffractive parton
densities therefore come from an H1 measurement of `ep → eXY` with
`x_pom < 0.05`, and `F-8` ("the DPDF tier is stock physics") is really a
statement about this paper. The tree cites `[H1FitB]` by key in
`PHYSICS_CHANNELS.md` §14 with no copy on disk; that is closed.

**GC-8 — Angantyr (1806.10820).** Included for the *negative* result.
`benchmarking/01` §2.2 establishes that a nucleus id on `Beams:idA/idB`
diverts into Angantyr, which has no lepton physics, and that Angantyr's
Glauber machinery counts wounded nucleons rather than computing an
eikonal survival factor for a composite cluster spectator
(`include/lipolgen/fsi.hpp`). The paper is the primary source for what
Angantyr actually models, so the next reader does not re-derive the
conclusion. Note the author list is **Bierlich, Gustafson, Lönnblad,
Shah** — verified against the arXiv record.

### 2.2 Peer generators

**GC-9 — BeAGLE (2204.11998).** The standard eA generator for the EIC:
DPMJet-III + PyQM + FLUKA evaporation with a nuclear spectral function.
LiPolGen needs it for three separable things: the **A = 2** tagged
spectator path, the Ciofi degli Atti–Simula n(k) and Strikman–Weiss
light-front prescription it implements, and far-forward fragment η /
rigidity spectra to place beside `E-16`'s published table. A
bibliographic correction, verified in this pass and already recorded in
`benchmarking/01` §2.3: the arXiv record **does** carry
doi:10.1103/PhysRevD.106.012007, so the sibling survey's "no journal-ref,
no DOI" note is stale. Cite PRD 106 (2022) 012007.

**GC-10 / GC-11 — DJANGO6 and HERACLES.** These are the two papers behind
"DJANGOH", the code `benchmarking/01` §2.4 calls "the strongest L1
candidate" and whose published Rad = 1 / Rad = 0 cross-section table is
recommendation #1 of that survey *(2026-09-23: wired and BLOCKED — the published table has the elastic radiative tail OFF, IEL2 = IEL31..33 = 0, so it shares no term with `rc_tail`; `../open_items/run_2026-09-23/phase_B3_chain_rc.md` §1.2.)*. Verified on INSPIRE: **DJANGO6** is
Charchula, Schuler & Spiesberger, CPC 81 (1994) 381 (recid 372027, 269
citations, CERN-TH-7133-94), an interface of HERACLES 4.4 to LEPTO 6.1;
**HERACLES** is Kwiatkowski, Spiesberger & Möhring, CPC 69 (1992) 155
(recid 296225, 451 citations, DESY-90-041), and its abstract is the one
that names "single-photon emission from the lepton line" and the complete
one-loop weak corrections — i.e. HERACLES, not DJANGO6, is the actual
independent RC engine that `D-3` and `C-8` would be compared against.
**Neither has a free copy**: both are Elsevier-era CPC articles with no
arXiv entry and no INSPIRE fulltext, and CERN's document server (which
would hold CERN-TH-7133-94) is behind an anti-bot challenge from this
machine. Download route: the DOI pages above, through an institutional
subscription.

**GC-12 — LEPTO 6.5 (hep-ph/9605286).** The DIS hard-process generator
that sits underneath *both* DJANGOH and PEPSI. It matters here because
LiPolGen's architecture is the same shape — a parton-level DIS cross
section handed to a Lund-string hadronizer — so LEPTO's treatment of the
electroweak DIS cross sections, the QCD matrix-element/parton-shower
options, and the target-remnant treatment is the closest published
precedent for what `LhaupDis` does. Free, on arXiv (DESY 96-057).

**GC-13 — PEPSI (CPC 71 (1992) 305).** This entry **closes a documented
UNVERIFIED**: `benchmarking/01` §2.5 says "the 1992 CPC citation itself
was taken from the README, not independently confirmed against the
journal". It is now confirmed — INSPIRE recid 321673, Mankiewicz, Schäfer
& Veltri, *PEPSI: A Monte Carlo generator for polarized leptoproduction*,
Comput. Phys. Commun. **71** (1992) 305,
doi:10.1016/0010-4655(92)90016-R, 63 citations, preprints
HD-THEP-91-47 / MPIH-V38-1991; the abstract confirms it is a modification
of LEPTO 4.3 in which "the hard virtual gamma-parton scattering is
generated according to the polarization-dependent QCD cross-section". No
free copy was located (no arXiv, no INSPIRE fulltext). PEPSI is the code
whose bundled `radgen.f` carries the only live tensor (`b₁`, `b₂`) branch
found anywhere, and §2.5's honest warning stands: that branch sets
`b₁ = −(3/2) F₁ A₁` from the *vector* asymmetry, so it is a placeholder,
and its value is the RC **kernel**, not its content.

**GC-14 — RADGEN 1.0 (hep-ph/9906408).** The POLRAD-2.0-derived radiative
event generator used by HERMES and COMPASS, and the code that physically
ships inside both BeAGLE (`radgen-6.4.28`) and PEPSI. It is the natural
external partner for the tree's POLRAD gate `D-3`, and the published
statement of how the polarized RC integrals are organized. Vector
polarization only; open access. (A copy is already on disk — a
concurrent pass fetched it; content verified.)

**GC-15 — RAPGAP (CPC 86 (1995) 147).** The diffractive-DIS generator
behind the official ePIC `DDIS/rapgap3.310` samples listed live in
`05_chain.md` §2.4, and the closest external analogue of the T2
γ*–Pomeron tier's *hadronic final state*. 598 citations. The journal
article is Elsevier, but **the DESY 93-182 preprint is free from INSPIRE**
as gzipped PostScript; it was fetched and converted to PDF here
(`jung1995_rapgap_desy93-182.pdf`, first page verified to be the RAPGAP
report). Jung's later v3.2 manual (INSPIRE recid 2653078) is also freely
downloadable if the newer parameter set is needed.

**GC-16 — STEG.** Searched and **not found as a paper.** INSPIRE returns
zero hits for a fulltext search on `SpectatorMC` or `LightIonEIC`, and no
literature record matches the generator. What exists is the repository
`github.com/JeffersonLab/LightIonEIC` (verified live via the GitHub API:
public, description "Light Ion EIC Code development", **no license**,
last push **2016-09-06**) and the JLab Spectator Tagging project page.
The physics it implements is C. Weiss's `TAG` package, so the citable
literature is the Cosyn–Weiss tagged-deuteron series **already in the
corpus** (`cosyn-weiss-tagged-tensor2026`) plus `[CW20]` PRC 102 (2020)
065204 — a theory-domain item, not fetched here. Treat STEG as *code with
no paper*: cite the repository and the theory it implements, never a
phantom reference.

**GC-17 / GC-18 / GC-19 — TOPEG, and the correction it carries.** The
task brief (and the earlier project brief) described TOPEG as the
Cosyn/Sargsian tagged-deuteron generator. **That is wrong**, and
`benchmarking/01` §5.2 already established it from the Yellow Report's
own text: TOPEG is the *Orsay-Perugia Event Generator* for **coherent
DVCS off ⁴He**, unpolarized, exclusive, no tagging, no deuteron. This
pass confirms there is **no TOPEG code paper**: an arXiv full-text search
for "TOPEG" returns nothing, and an INSPIRE full-text search for
"Orsay-Perugia event generator" returns only papers that *use* or *name*
it — the Yellow Report (2103.05419), the ECCE detector paper
(2208.14575), and Fucini–Scopetta–Viviani's incoherent-DVCS paper
(2008.11437). The citable primary sources are therefore the physics
papers TOPEG implements: GC-18 (coherent, PRC 98 (2018) 015203, 24
citations, an IA convolution for the nuclear GPD built on the ⁴He
one-body non-diagonal spectral function from AV18) and GC-19 (incoherent,
PRC 102 (2020) 065205). Both open access; both now on disk. The code the
brief was reaching for is STEG (GC-16).

**GC-20 — lAger.** Verified live via the GitLab API:
`eicweb.phy.anl.gov/monte_carlo/lager`, public, created 2020-02-07, last
activity 2025-11-21. **No INSPIRE literature record exists** — a search
for a lAger generator paper returns zero hits. Its published pedigree is
indirect: the Yellow Report cites it for the J/ψ and Υ projections, and
its process modules wrap published calculations (jpacPhoto, holographic
VM, the Lee ⁴He J/ψ grid). Cite the repository plus the YR; do not invent
a code paper.

**GC-21 / GC-22 / GC-23 — the exclusive/coherent trio.** These belong
primarily to the coherent-diffraction domain and are listed here because
they are *generators*; expect overlap with that domain's file.
eSTARlight (GC-21) is the in-tree `D-4` exclusive-VM rate scale (85–97 % of its rate at Q² < 0.1 GeV², below `coherent.hpp`'s q2_min = 0.7 — not a coherent-channel baseline); STARlight (GC-22) is
its Q²→0 parent with by far the strongest data pedigree (a decade of
RHIC/LHC UPC measurements); Sartre (GC-23) matters for its **negative**
statement — §3.3.4 of the code paper enumerates the supported nuclei (p,
O, Al, Ca, Cu, Cd, Au, Pb) and A = 6 is not among them, which is why
"use Sartre for ⁶Li" is closed as a line of enquiry. All three were on
disk already or arrived in this batch.

### 2.3 The record and the chain

**GC-24 — HepMC3 (1912.08005).** The library LiPolGen writes with
(3.3.0 in `deps/install`) and the format the whole ePIC chain consumes.
The paper specifies the simplified event-record structure, the
**arbitrary-attribute mechanism** that `HEPMC3_CONVENTION.md` uses to
carry ion spin, and the I/O back-ends (Asciiv3, ROOT, protobuf) that gate
`F-3`'s round trips exercise. `05_chain.md` §1.3 also points at the
paper's own validation tooling — a generation/analysis/plotting/**timing**
suite — as the template for the cross-version regression LiPolGen does
not yet have. A core dependency with no copy in the corpus; closed.

**GC-25 — NuHepMC (2310.13211).** The single most directly reusable
document for `F-4`. It is a *standard* layered on HepMC3: conventions
named `<Component>.<Category>.<Index>` over G/E/V/P, split into
**R**equirements, **C**onventions and **S**uggestions, with a reference
implementation and a `Validator` tool that checks a file against the
standard. LiPolGen's ion-spin attribute schema is today a proposal with
no machine-checkable conformance test; NuHepMC is what "propose it to the
ePIC MC group" should look like as a deliverable (`CH-18`). Published:
SciPost Phys. Codebases 57 (2025).

**GC-26 — Rivet 4 (2404.15984).** Rivet 4.1.2 is already installed on
this box inside `eic_xl-nightly.sif`, with 46 H1/ZEUS DIS analyses, and
`05_chain.md` §3 calls it "the largest finding in this survey". The
release note is the citable reference for the version actually present
and for the API changes (the projection/analysis generalization) that a
LiPolGen Rivet plugin would target. Also relevant to §3.2 of that file:
Rivet cannot read a ⁶Li beam — a measured limitation, so the Rivet route
runs on the ep surrogate, not on the lithium record. Published: SciPost
Phys. Codebases 36 (2024).

**GC-27 / GC-28 — eic-smear and abconv.** Both verified live via the
GitHub API in this pass: `eic/eic-smear` **GPL-3.0**, last push
2026-07-07; `eic/afterburner` **no license declared**, last push
2026-07-05. **Neither has a paper** — INSPIRE returns no literature
record for either, and a Zenodo search for an eic-smear software DOI
returns nothing relevant. This matters for how the benchmark report is
written: `CH-12` and `CH-23` must cite repositories and commit hashes,
not literature. eic-smear is the leg LiPolGen deliberately skips (it
writes HepMC3 directly), which is why `F-1`/`F-2` are the *only* evidence
that its record is shaped like everyone else's; abconv is the crossing-
angle afterburner with **no light-ion preset**, the open item behind
`CH-12`.

**GC-29 / GC-30 — DD4hep and Geant4.** The two layers `npsim` is built
from, and therefore what gate `F-1` actually exercises once the record
leaves LiPolGen. Neither is cited anywhere in the tree today. Both are
open access (IOP CC-BY and an Elsevier open article respectively); low
priority, listed so the chain section of the benchmark report can cite
the software it ran rather than naming binaries.

**GC-31 — the PDG numbering scheme.** `05_chain.md` §5.1 records that the
tree's ion PDG codes are "settled, and the tree is right"; the authority
for `10LZZZAAAI` is the PDG's *Monte Carlo particle numbering scheme*
review. Cite the 2024 edition (Navas *et al.*, PRD 110 (2024) 030001); a
2026 edition exists (Int. J. Mod. Phys. A **41** (2026) 2630011, INSPIRE
recid 3193100) if a current-edition citation is wanted.

### 2.4 PDFs, nuclear PDFs and R

These are not "background reading" — each one is a *named set string* in
the code, so the citation and the run configuration have to agree.

**GC-32 — LHAPDF6 (1412.7420).** The library `src/lhapdf/lhapdf_sf.cpp`
links (6.5.5 here). Needed for the interpolation and error-set semantics
that `Epps21Ratio`'s Hessian band relies on.

**GC-33 — CT18 (1912.10053).** `LhapdfSF`'s default is **`CT18NLO`**
(`lhapdf_sf.cpp` header comment), and EPPS21's proton baseline is
`CT18Anlo` — so this one paper stands behind both the free-nucleon F₂ and
the nuclear-ratio denominator. The CT18/CT18A distinction (the latter
includes the ATLAS 7 TeV W/Z data) is the reason the grid name in
`EPPS21nlo_CT18Anlo_Li6` reads as it does.

**GC-34 — MSTW 2008 (0901.0002).** This is a *gate-defining* reference,
not a convenience: `include/lipolgen/b1_nuclear.hpp` states that the
CDKS-comparison gate passes for the **MSTW2008 LO** configuration
specifically ("specific to CD-Bonn AND MSTW together"), because MSTW2008
LO is the unpolarized PDF CDKS computed their b₁ with. The `--b1-unpol
mstw` path reads PYTHIA's own bundled `pdfdata/mstw2008lo.00.dat`, so the
gate's numbers are tied to this fit and a substitution changes them.

**GC-35 — NNPDFpol1.1 (1406.5539).** `LhapdfG1`'s default
`NNPDFpol11_100`. The tree's header note that "Δb = 0 is NNPDFpol1.1's
own convention" is a statement about this fit, and the paper is where the
determination (and its replica-based uncertainties) is defined. Note the
arXiv record carries the DOI but **no journal_ref**; the published
reference is Nucl. Phys. B 887 (2014) 276.

**GC-36 — EPPS21 (2112.12462).** The nuclear PDF fit behind
`Epps21Ratio` and `EmcBaseline::Epps21` — the generator's single
unpolarized ⁶Li EMC baseline today, with the caveat `04_theory.md` T-39
insists on: almost no data constrains A ≈ 6, so `EPPS21nlo_CT18Anlo_Li6`
is an A-dependent parametrization anchored on essentially one NMC ⁶Li/D
ratio measurement.

**GC-37 / GC-38 / GC-39 — nCTEQ15, nNNPDF3.0, TUJU19.** The other three
fits of the `T-39` four-fit cross-check, and the two that carry genuine
**⁷Li** grids (`nCTEQ15_7_3`, `TUJU19_*_7_3`) for `T-40`'s correction
that ⁷Li is not PDF-less. nNNPDF3.0's journal reference had to come from
INSPIRE (the arXiv record has neither journal_ref nor DOI, and its title
there carries a typo, "NNNPDF3.0"); the published item is EPJ C 82 (2022)
507, doi:10.1140/epjc/s10052-022-10417-7. The honest framing from T-39
belongs in any use of these: four fits sharing one dataset agree more
than the evidence warrants.

**GC-40 — R1990 (Whitlow *et al.*, PLB 250 (1990) 193).** A real gap
closed. `include/lipolgen/sf.hpp:43` describes its toy `r_sigma_lt` as a
"simplified R1990-like magnitude", `02_data_nucleon_deuteron.md` row R-1
calls this "the closest thing to a reference" for that claim, and
`run_2026-09-03/phase_A_miller_normalisation.md:527` says outright
"Whitlow's paper, which I did not fetch" while nevertheless
re-implementing `R1990(x,Q2)` at line 171 as HERMES's ref [18]. The paper
is Elsevier, but **INSPIRE hosts the SLAC-PUB-5284 preprint free**; it
was fetched here and its first page verified (title, five authors,
abstract: "R = σ_L/σ_T from a global analysis of eight SLAC deep
inelastic experiments … 0.1 ≤ x ≤ 0.9 and 0.6 ≤ Q² ≤ 20.0"). The
transcribed `R1990` function can now be checked against the published
parametrization rather than against a memory of it.

**GC-41 — R1998.** Already in the corpus with a PDF; listed only so the
R pair is complete in one place. It is the default R in
`src/core/sf.cpp` (`kR1998A/B/C`, the three forms and their spread).

### 2.5 The EIC/ePIC documents

All seven (GC-42…GC-48) are already corpus entries; see `00_corpus.md`
Table 1 for their local-file status. Three points are worth carrying
forward rather than re-deriving:

1. The **Yellow Report** (GC-42) is the only one of them the generator
   domain actually reads: §8.1.6's DJANGOH/PYTHIA6-vs-HERA comparison is
   the published "generator verification" for the EIC, and the YR text is
   also where TOPEG, Sartre 1.34 and lAger are described.
2. The **CDR, PDR and TDR** (GC-43…GC-45) are deliberately not committed
   on size grounds (178 / 114 / 260 MB) and are correctly recorded with
   `file: null` and a DOI. The far-forward *numbers* the tree uses come
   from Chang *et al.* (GC-47) and the YR, not from these, so nothing is
   blocked by their absence.
3. There is **no single "ePIC detector paper"** to cite: an INSPIRE
   search returns only sub-detector and computing proceedings. The
   detector document is the TDR (GC-45).

---

## 3. Negative and unverified results

### 3.1 Paywalled, no free copy located

| reference | why no copy | download route for the user |
|---|---|---|
| **DJANGO6**, CPC 81 (1994) 381 | no arXiv; no INSPIRE fulltext. Report number CERN-TH-7133-94 would be on CDS, but `cds.cern.ch` served an Anubis anti-bot challenge from this machine and was not bypassed | doi:[10.1016/0010-4655(94)90086-8](https://doi.org/10.1016/0010-4655(94)90086-8); or search CDS for `CERN-TH-7133-94` in a browser |
| **HERACLES**, CPC 69 (1992) 155 | no arXiv; no INSPIRE fulltext. Report number DESY-90-041 | doi:[10.1016/0010-4655(92)90136-M](https://doi.org/10.1016/0010-4655(92)90136-M); or the DESY publication database for `DESY-90-041` |
| **PEPSI**, CPC 71 (1992) 305 | no arXiv; no INSPIRE fulltext. Report numbers HD-THEP-91-47, MPIH-V38-1991 | doi:[10.1016/0010-4655(92)90016-R](https://doi.org/10.1016/0010-4655(92)90016-R) |

For all three the *bibliographic record* is now verified against INSPIRE
even though the PDF is not obtainable — which is what `benchmarking/01`
§2.5 asked for in the PEPSI case.

### 3.2 Code with no paper — do not cite literature that does not exist

| code | verified state (2026-09-16) | what to cite instead |
|---|---|---|
| **STEG** / `JeffersonLab/LightIonEIC` | GitHub API: public, no license, last push 2016-09-06. INSPIRE fulltext search for `SpectatorMC` / `LightIonEIC`: **0 hits** | the repository + the Cosyn–Weiss tagged-deuteron theory papers it implements |
| **lAger** | GitLab API: `monte_carlo/lAger`, public, created 2020-02-07, last activity 2025-11-21. INSPIRE: **0 hits** for a code paper | the repository (GPL-3.0-or-later) + the Yellow Report's J/ψ/Υ projections |
| **eic-smear** | GitHub API: GPL-3.0, last push 2026-07-07. INSPIRE: no record; Zenodo: no software DOI found | the repository + commit hash |
| **afterburner / abconv** | GitHub API: no license, last push 2026-07-05 | the repository + the version string (`abconv 0.1.3` in `eic_xl-nightly.sif`) |
| **TOPEG** | arXiv full-text "TOPEG": **0 hits**. INSPIRE full-text "Orsay-Perugia event generator": 5 hits, all *users* of it (YR, ECCE, CORE, the 22 GeV workshop, Fucini *et al.*) | GC-18 / GC-19 (Fucini, Scopetta & Viviani) — the physics TOPEG implements — plus the YR §8 description |

### 3.3 Corrections this pass makes to statements in the tree

1. **PYTHIA 8.3 is published.** `benchmarking/01` §2.1 cites it as an
   arXiv-only item because the arXiv record has no `journal_ref`. INSPIRE
   gives SciPost Phys. Codebases 8 (2022),
   doi:10.21468/SciPostPhysCodeb.8.
2. **PEPSI's CPC citation is confirmed.** §2.5's UNVERIFIED flag can be
   lifted: CPC 71 (1992) 305, doi:10.1016/0010-4655(92)90016-R, INSPIRE
   recid 321673.
3. **TOPEG is not a tagged-deuteron generator**, and has no code paper —
   confirming §5.2 independently and adding the two Fucini–Scopetta–
   Viviani papers as the citable primaries.
4. **nNNPDF3.0's arXiv title contains a typo** ("NNNPDF3.0") and the
   record has no journal_ref; the published item is EPJ C 82 (2022) 507.
5. **Angantyr's fourth author is Shah**, not Tarasov (arXiv record).

### 3.4 Not fetched, deliberately

GC-6 (2605.00502, no journal ref yet), GC-29 and GC-30 (DD4hep, Geant4 —
open access, but low priority and not cited anywhere in the tree today),
GC-31 (the PDG review, a 2000-page document better used online), and
GC-43…GC-45 (the CDR/PDR/TDR, 178–260 MB, already recorded with
`file: null` in the corpus by deliberate policy).

---

## 4. Priority, if only a few of these are acted on

1. **GC-40 R1990** — closes a gap the tree names twice in its own words
   ("Whitlow's paper, which I did not fetch"), and makes `sf.hpp:43`'s
   "R1990-like" claim checkable. Free, on disk, hours of work.
2. **GC-24 HepMC3 + GC-25 NuHepMC** — the record LiPolGen writes, and the
   only existing template for turning `HEPMC3_CONVENTION.md`'s spin
   attributes into a machine-checkable standard (`F-4`, `CH-18`).
3. **GC-1 PYTHIA 8.3 + GC-2/GC-3 Les Houches + GC-4 dipole recoil** — the
   bridge's full documentary basis, all previously uncited.
4. **GC-34 MSTW2008 + GC-33 CT18 + GC-35 NNPDFpol1.1** — the three PDF
   fits that *define* what the b₁/CDKS gate and the two LHAPDF backends
   compute. Without them a gate number cannot be reproduced.
5. **GC-10/GC-11 DJANGOH+HERACLES** — the only independent one-loop RC
   engine; needs a subscription, but the bibliographic record is now
   exact enough to request it.
6. **GC-7 H1 2006 DPDF** — the actual content of `PDF:PomSet = 6`, i.e.
   of the whole T2 Pomeron tier.

---

## 5. Files added to the shared corpus in this pass

All under `/home/cpeng/Projects/polli/PolarizedLithiumSim/refs/`, all
open access, all fetched from arXiv or INSPIRE and content-checked
(first page read) after download. **`refs_dict.json` was not modified.**

```
2203.11601.pdf   PYTHIA 8.3 manual                      3.4 MB
hep-ph_0609017.pdf  Les Houches Event Files (LHEF)      103 kB
hep-ph_0109068.pdf  Les Houches Accord 1 (common blocks) 172 kB
2204.11998.pdf   BeAGLE                                  1.3 MB
hep-ph_9605286.pdf  LEPTO 6.5                            296 kB
hep-ex_0606004.pdf  H1 2006 diffractive PDFs             1.1 MB
1710.00391.pdf   Cabouat & Sjostrand, dipole showers     1.2 MB
1806.10820.pdf   Angantyr                                1.0 MB
2410.20950.pdf   PYTHIA multi-jet merging in DIS         1.9 MB
1805.05877.pdf   Fucini et al., coherent DVCS off 4He    211 kB
2008.11437.pdf   Fucini et al., incoherent DVCS off 4He  435 kB
1912.08005.pdf   HepMC3                                  533 kB
2310.13211.pdf   NuHepMC                                 2.0 MB
2404.15984.pdf   Rivet 4 release note                    172 kB
1412.7420.pdf    LHAPDF6                                 697 kB
1912.10053.pdf   CT18                                  12.5 MB
0901.0002.pdf    MSTW 2008                               2.2 MB
1406.5539.pdf    NNPDFpol1.1                             965 kB
2112.12462.pdf   EPPS21                                  5.4 MB
1509.00792.pdf   nCTEQ15                                 3.7 MB
2201.12363.pdf   nNNPDF3.0                               3.3 MB
1908.03355.pdf   TUJU19                                  2.8 MB
whitlow1990_R.pdf            R1990 (SLAC-PUB-5284)       580 kB
jung1995_rapgap_desy93-182.pdf  RAPGAP (DESY 93-182)     248 kB
```

Already on disk before or during this pass (verified, not re-fetched):
`hep-ph_9906408.pdf` (RADGEN), `1803.06420.pdf` (eSTARlight),
`1607.03838.pdf` (STARlight), `1307.8059.pdf` (Sartre),
`hep-ex_9808028.pdf` (R1998), `2103.05419_part{1..4}.pdf` (Yellow
Report), `2511.05638.pdf` (IR-8 light nuclei).

---

## 6. What this domain still cannot supply

Nothing in this file is a *polarized*-nucleus generator reference,
because no such generator exists (`benchmarking/01` §5.1: no MCEG
anywhere implements a tensor-polarized target of any species). The
strongest statements available from the generator/chain domain remain:

- the **shower and hadronization** half of the bridge has a HERA
  pedigree (GC-4, GC-5);
- the **unpolarized RC magnitude** can be bounded against an independent
  one-loop engine (GC-10/GC-11), and the RC *kernel* against a second
  copy of POLRAD (GC-13/GC-14);
- the **record** is well formed and conforms to a published standard
  (GC-24), with a template for making the spin convention checkable
  (GC-25);
- the **PDF inputs** are pinned to named, citable fits (GC-32…GC-39).

The tensor-polarized ⁶Li/⁷Li observable the generator exists to produce
has no external code and no external data, and no reference in this file
changes that.
