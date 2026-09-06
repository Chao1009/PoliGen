<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# 01 — External generators: what each one could and could not benchmark

**Survey date 2026-09-06.** Second file in the benchmarking series; the baseline
is `docs/benchmarking/00_in_tree_checks.md`, whose row labels (`D-1` … `D-6`,
`E-1` … `E-18`, `F-1` … `F-8`) are used here without repeating their content.
Scope: **event generators** — external codes that produce events or cross
sections overlapping a LiPolGen channel. Measured data, theory tables and
software-chain artefacts get their own files.

Everything below was verified live on the survey date by fetching the
repository, the API record, the manual, or the source file named in the entry.
Anything that could not be verified is marked **UNVERIFIED** in place, with the
reason. No generator was built or downloaded for this survey; the only code run
was LiPolGen's own test binary (§3.1).

---

## 0. The frame this file has to be read in

`00_in_tree_checks.md` §0 established three facts. This file adds a fourth,
which is the one that governs the whole generator domain.

1. There is **no polarized e + ⁶Li or e + ⁷Li DIS data anywhere.**
2. **No other generator does tensor-polarized lithium.**
3. The only tensor-polarized DIS datum in the world is HERMES's deuteron b₁.
4. **NEW, and stronger than (2): no Monte-Carlo event generator anywhere
   implements a tensor-polarized target of any species.** Not deuteron, not
   ³He, not ⁶Li. The search behind this claim is recorded in §5.1; the single
   near-miss — a dormant `b₁`/`b₂` branch inside the POLRAD-derived RADGEN
   found in both BeAGLE and PEPSI — is §2.5, and it computes b₁ from the
   *vector* asymmetry, i.e. it is a placeholder, not a tensor model.

So **nothing in this file benchmarks the observable LiPolGen exists to
produce.** What follows is the decomposition: the four levels at which an
external generator *can* say something, and exactly where each stops.

| level | what an external generator can check | best available code |
|---|---|---|
| **L1 — nucleon** | polarized inclusive DIS: g₁/A₁ from polarized PDFs, the depolarization factor, and the **size of the radiative correction** | DJANGOH 4.6.10+ (§2.4), PEPSI (§2.5) |
| **L2 — deuteron** | tagged-spectator kinematics from a deuteron wave function, light-front spectator prescription, **vector**-polarized tagged asymmetry | STEG / `JeffersonLab/LightIonEIC` (§2.6), BeAGLE eD (§2.3) |
| **L3 — unpolarized nucleus** | coherent VM production and \|t\| slopes, Fermi smearing, nuclear breakup / far-forward fragment kinematics | eSTARlight (§2.7), STARlight (§2.8), Sartre (§2.9), lAger (§2.10), BeAGLE (§2.3) |
| **L4 — chain** | the record is well formed and every downstream EIC tool accepts it | eic-smear (§2.12), abconv + npsim + EICrecon (§2.13) |

Nothing sits at an "L5 — polarized nucleus" level, because no such code exists.

---

## 1. The map — which LiPolGen channel has an external code at all

| LiPolGen channel | closest external generator | the most it could establish | the gap that stays open |
|---|---|---|---|
| Inclusive DIS, **vector** (g₁, A₁) on ⁶Li/⁷Li | DJANGOH (free-nucleon polarized PDFs), PEPSI | that the **nucleon-level** g₁/A₁/D construction and its RC are standard | the nuclear convolution, the effective polarizations, the cluster decomposition |
| Inclusive DIS, **tensor** (b₁, b₂, A_zz) | **none** | — | everything. §5.1 |
| Tagged α/d/t spectator with cluster wave functions | STEG (deuteron, S-wave Hulthén), BeAGLE eD | that the tagged **kinematics and spectral-function machinery** is the standard one at A = 2 | the α+d / α+t cluster wave function; the D-state; the tensor spectator asymmetry |
| T1 breakup | BeAGLE (FLUKA evaporation), Sartre (GEMINI++) | that ⁶Li* → α+d energetics and far-forward fragment η/rigidity are physical | pre-formed clusters — both codes only make clusters by *statistical evaporation* |
| Glauber FSI on the spectator | BeAGLE (INC, A > 12 only), Sartre (Glauber for the dipole) | nothing directly: neither implements eikonal *cluster*-spectator survival | the whole FSI weight. §5.4 |
| Coherent ⁶Li diffraction | eSTARlight (already in tree, `D-4`), STARlight (γ limit), Sartre (saturation) | the unpolarized coherent rate, the \|t\| slope and the VM ratios; and a **second dipole model** at the deuteron | ⁶Li in Sartre (§5.3); the tensor a₂ modulation everywhere |
| PYTHIA Pomeron tier (T2 coherent) | stock PYTHIA 8 diffraction, RAPGAP | that the DPDF tier is stock physics (already `F-8`) | — |
| Radiative corrections, **vector** | DJANGOH/HERACLES (independent engine), PEPSI/RADGEN, POLRAD 2.0 itself | the **magnitude and Q² trend** of the RC against a second, independently written one-loop engine | — |
| Radiative corrections, **tensor band** | RADGEN's `pnrun = ±2` branch (§2.5) | the *kernel structure* — how b₁, b₂ enter the RC integrals | the tensor structure functions themselves; RADGEN's b₁ is a placeholder |
| HepMC3 → ePIC | eic-smear `TreeToHepMC`, abconv, npsim | that LiPolGen's record has the same shape every other EIC generator's does | the spin attributes, which no downstream tool carries (already `F-4`) |

---

## 2. Verified inventory

Each entry: **repository / obtainability → activity → license → validation
pedigree (cited) → exact observable overlap → runnable here? → what it does NOT
validate.**

### 2.1 PYTHIA 8.317 — already in tree

* **Repository / obtainability.** `IN-TREE`: built at
  `/home/cpeng/Projects/polli/deps/install` (`libpythia8.so`, `pythia8.so`,
  `libpythia8hepmc3.so`, `libpythia8lhapdf6.so`), sources at
  `/home/cpeng/Projects/polli/pythia/pythia8317`. Upstream `pythia.org`,
  GPL-2-or-later plus the MCnet Guidelines (`COPYING`, `GUIDELINES` in the
  source tree; recorded in `docs/surveys/pythia8_survey.md` §6).
* **Activity.** 8.317, copyright 2026; version constant
  `PYTHIA_VERSION 8.317` in `include/Pythia8/Pythia.h`.
* **Pedigree.** The 8.3 physics-and-usage manual, arXiv:2203.11601 (Bierlich
  *et al.*; verified via the arXiv API — the record carries **no journal
  reference**, so it is cited as the arXiv item). The DIS-specific pedigree is
  the dipole-recoil initial-state shower of Cabouat & Sjöstrand, *Some Dipole
  Shower Studies*, **EPJ C 78 (2018) 226**, arXiv:1710.00391 (verified), whose
  abstract states the DIS comparison to data explicitly. `SpacelikeShowers.xml`
  in the shipped manual says dipole recoil "for the first time allows the
  simulation of Deeply Inelastic Scattering processes in PYTHIA 8".
* **Overlap.** Shower + hadronization of the LiPolGen hard process (the LHAup
  bridge, `docs/PYTHIA_BRIDGE.md`), and the coherent T2 Pomeron tier.
* **Runnable here.** Yes — `D-1` is a live gate. Re-run on the survey date, §3.1.
* **Does NOT validate.** Anything spin-dependent (`SigmaProcess` is
  helicity-summed; there are no polarized PDFs and no `Beams:polA/polB` — see
  `pythia8_survey.md` §5); any cross section under LHAup strategy 3
  (`info.sigmaGen()` is not physical there — `D-1`'s "does not validate"
  column); e + A beams (a nucleus id diverts into Angantyr, which then fails
  `BeamSetup::checkBeams`); and α+d geometry (`HINucleusModel.cc`'s
  `ClusterModel` is ⁴He-only, unreachable from `create()`, and would
  null-deref).

### 2.2 Angantyr (inside PYTHIA 8) — for completeness, and it does not help

* **What it is.** PYTHIA 8's Glauber heavy-ion module, `src/HeavyIons.cc`.
* **Overlap with LiPolGen.** Superficially the Glauber FSI module. In fact
  none: Angantyr counts *wounded nucleons* in a hadron-nucleus geometry, while
  `include/lipolgen/fsi.hpp` computes an **eikonal survival product for a
  composite cluster spectator** (`E-15`). Different object, different
  observable.
* **Blocking facts.** `HeavyIons.xml` restricts `idA` to a hadron; there are
  zero occurrences of lepton/DIS/photon physics in `HeavyIons.cc`. Light-nucleus
  geometry `HOShellModel` has a built-in ⁶Li charge radius but **none for ⁷Li**.
* **Priority.** Low. Listed so the next reader does not re-derive this.

### 2.3 BeAGLE

* **Repository.** Upstream is **`https://gitlab.in2p3.fr/BeAGLE/BeAGLE`**
  (verified via the GitLab API: public, created 2018-09-14, **last activity
  2023-08-17**, latest commits `dbd4750e` 2023-08-16 "Improve (reduce) default
  printout", `522d47e4` 2023-08-10 "Fix E* for special cases", `120fee86`
  2023-07-10 "Add ESTARFIX card"). `https://github.com/eic/BeAGLE` is a
  **snapshot mirror**: created *and* last pushed 2024-01-17, 95 commits.
  `docs/surveys/beagle_survey.md` inventories the GitHub clone; the GitLab
  upstream is the newer of the two homes and is worth noting when that survey
  is next revised.
* **License.** **None**, on either host (GitHub API `license: null`; GitLab API
  `license: None`). Redistribution of a fork is legally murky; FLUKA is
  non-redistributable outright.
* **Pedigree.** Chang *et al.*, arXiv:2204.11998 — and it **is published**:
  the arXiv record carries **DOI 10.1103/PhysRevD.106.012007**, i.e. *Phys.
  Rev. D* **106** (2022) 012007. (`beagle_survey.md` §6, written against v1,
  says "no journal-ref, no DOI"; that is now stale.) The paper's own validation
  is against HERA/EIC-projection distributions and eD spectator kinematics; it
  is the standard eA reference generator for the EIC.
* **Overlap.** The genuinely reusable pieces are all at **A = 2**: the tagged
  deuteron spectator path (`src/dpm_pythia.f:1180-1305`, including the
  design rule "leave the light-nucleus spectator unmodified — do not ptkick or
  recoil-correct it"), the Ciofi degli Atti–Simula n(k) for A = 2, 3, 4
  (`DT_KFERMI`), the Strikman–Weiss light-front spectator prescription
  (`DT_SPECTRALFUNC`), and the far-forward fragment η / rigidity distributions
  that `E-16` currently gets from a paper table.
* **Runnable here.** **Yes, with one registration.** `eic/beagle-container`
  (pushed 2026-06-12) is a complete `Dockerfile` — verified line by line — that
  builds CERNLIB (the free 2024.06.12.0 source), LHAPDF **5.9.1** with seven
  `sed` fixes, PYTHIA 6.4.28 from BeAGLE's bundle, RAPGAP 3.302, and BeAGLE
  itself from the GitLab upstream, on AlmaLinux 9. The **only** thing it cannot
  fetch is FLUKA: `COPY fluka2011.2x-linux-gfor64bit-8.3-AA.tar.gz /tmp/`
  requires the user to place that licensed 2011-era tarball in the build
  context. `singularity-ce 3.11.2` is installed here, and the README gives the
  `docker save` → `singularity build` route. So the FLUKA registration is the
  entire cost; every other blocker listed in `beagle_survey.md` §1 is solved
  upstream.
* **Does NOT validate.** Anything polarized (grepping BeAGLE's own sources for
  `polariz|helicit|spin|g1|b1` returns **zero** hits outside RADGEN); ⁶Li/⁷Li
  structure of any kind (A ≥ 5 takes a plain Woods–Saxon branch, lithium has no
  entry in PyQM's density table and falls to a generic Bialas parameterisation,
  and the A > 4 Fermi dispatcher has an `.OR.`/`.AND.` bug that hands A = 6, 7
  the **Fe–Pb** momentum distribution — `beagle_survey.md` §2c); pre-formed
  clusters (clusters appear only as FLUKA evaporation products); and coherent
  or intact-nucleus recoil, which the paper says is impossible by construction.

### 2.4 DJANGOH 4.6.10 / 4.6.21 — **the strongest L1 candidate**

* **What it is.** HERACLES (complete one-loop electroweak + QED radiative
  corrections and radiative scattering) interfaced to LEPTO 6.5.1 and JETSET
  fragmentation, with SOPHIA at low hadronic mass.
* **Obtainability.** **ON REQUEST / partially blocked.** The author's page
  `www.desy.de/~hspiesb/mcp.html` returns **403** from this machine (verified
  twice, with a browser user agent), and there is no public source mirror:
  a GitHub search for a DJANGOH/HERACLES repository returns nothing, there is
  **no `djangoh` package in `eic-spack`**, and `eic-shell`'s environment
  (`eic-spack/environments/eic/spack.yaml`, verified) ships only `eic-smear`,
  `afterburner`, `root +pythia8` and `hepmc3` — not DJANGOH. What **is** public
  and was verified: the three manuals in `github.com/eic/documents`
  (`software/general/Djangoh_m.pdf` 196 kB, `Djangoh_Updateh-4.6.10.pdf`
  60 kB, `Lepto-6.5.1.pdf` 308 kB), the ePIC steering cards and production
  notes (`eic/djangoh_production`, pushed 2026-09-04; `eic/InclusiveDjangohSamples`,
  pushed 2025-12-07), and the software page `eic.github.io/software/djangoh.html`.
* **Activity.** ePIC's 2026 production release tag is `DJANGOH4.6.10-2.0`
  (charged-current 9×275); the inclusive NC samples were produced with
  **4.6.21** (`InclusiveDjangohSamples/README.md`). So the code is alive and in
  routine EIC production, even though the source is not publicly posted.
* **License.** Not stated on any page reachable from here. **UNVERIFIED.**
* **Pedigree.** EIC Yellow Report §8.1.6 "Generator verification"
  (arXiv:2103.05419, read locally from `PolarizedLithiumSim/refs/`): PYTHIA6 and
  DJANGOH reduced inclusive e⁺p NC cross sections were compared **to HERA data**,
  Fig. 8.23, ≈10 pb⁻¹ of pseudo-data each, `cteq61` grid. The YR also records
  that DJANGOH's forward particle yields are "nearly identical" to PYTHIA6's.
* **Overlap — and this is the specific, valuable part.** Two things:
  1. **Polarized inclusive DIS with radiative corrections.** From 4.6.10 the
     *proton* can be polarized (`Djangoh_Updateh-4.6.10.pdf` §2, verified):
     `PR-BEAM` takes a second real `HPOLAR ∈ [-1, +1]`, `POLPDF` selects one of
     nine polarized PDF sets (DSSV, DNS, DS, GSLO, GSNLO, BB, AAC, LSS, GRSV).
     The per-event output carries **`G1NC`, `G3NC`, `A1NC`, `G1CC`, `G5CC` and
     the depolarization factor `D`** as named columns
     (`eic.github.io/software/djangoh.html`, verified). That is a
     column-for-column comparison target for LiPolGen's `g1`, `A1` and the
     finite-γ depolarization of `E-8`, computed by an engine that shares no
     code with this tree.
  2. **A published, reusable RC magnitude.** `eic/InclusiveDjangohSamples`
     publishes the four-bin cross-section table for e 18 GeV × p 275 GeV NC,
     CTEQ61M (`10150`), x ∈ [1e-5, 1], y ∈ [1e-4, 1], W > 3.0 GeV (kinematics
     read from `runcards/ep.Rad=1.NC.Q2_1_10.in`), with radiation on and off:

     | Q² [GeV²] | Rad = 1 [µb] | Rad = 0 [µb] | ratio |
     |---|---|---|---|
     | 1 – 10 | 0.6733778912 | 0.6278853416 | 1.0725 |
     | 10 – 100 | 0.07405000239 | 0.06797276736 | 1.0894 |
     | 100 – 1000 | 0.003760005602 | 0.003239636940 | 1.1606 |
     | 1000 – 10000 | 0.00008615325533 | 0.00007015936779 | 1.2280 |

     (ratios computed here from the published numbers.) This is an
     **obtainable-today table** against which the magnitude and Q² trend of
     LiPolGen's unpolarized RC can be sanity-checked without running anything.
* **Runnable here.** Not without the source. If obtained: gfortran 11.4.0 is
  present, but it needs **LHAPDF 5** (tested against 5.8.6; only LHAPDF 6.5.5 is
  installed here) and CERNLIB `ranlux` — the same two dependencies the BeAGLE
  container builds from source, so that container is the natural host.
* **Does NOT validate.** Any tensor quantity — DJANGOH has none. Any *nuclear*
  polarization: the 4.6.10 note is explicit that polarization enters through
  **free-nucleon** polarized PDFs, while the nuclear machinery (`NUCL-MOD`,
  EKS98/EPS08/EPS09/HKN) modifies only the **unpolarized** PDFs. Running
  `POLPDF` with `NUCLEUS A>1` therefore gives free-nucleon Δq inside a nucleus
  with no nuclear correction — a documented limitation, not a model of a
  polarized nucleus. Also: no clusters, no spectator tagging, no coherent
  diffraction.

### 2.5 PEPSI + its RADGEN — **the only code with a live tensor switch**

* **Repository.** `https://github.com/eic/PEPSI` (verified: description "MCEG
  for polarized DIS", **no license**, last push **2020-04-10**, default branch
  `master`). Packaged in `eic-spack` as `pepsi` (`depends_on("cernlib")`,
  builds with `make`, patches `-g -m64` → `-fallow-argument-mismatch` for
  gcc ≥ 10).
* **Pedigree.** The README cites Mankiewicz, Schäfer & Veltri, *Comp. Phys.
  Comm.* **71** (1992) 305, and links `PEPSI.paper.pdf` on the BNL EIC wiki.
  PEPSI is LEPTO 6.5 with polarized structure functions and is the generator
  behind CLAS's `clasdis`; the EIC Yellow Report records "detailed impact
  studies that use PEPSI as polarized" MC (arXiv:2103.05419, read locally).
  **UNVERIFIED:** the 1992 CPC citation itself was taken from the README, not
  independently confirmed against the journal.
* **Overlap — three distinct things, of very different value.**
  1. **Polarized DIS event generation.** `pepsiMaineRHIC_noradcorr.f` with
     polarized PDFs; the same L1 target as DJANGOH, and the two are genuinely
     independent implementations.
  2. **A *working* polarized-asymmetry hook.** `pepsi_radgen_extras.f:250-303`
     implements `mkasym(dQ2, dX, A, Z, dA1, dA2)` for real: it calls `mkf2`,
     then builds g₁ from the polarized parton array `xdpq` with an explicit
     Z·(proton) + (A−Z)·(neutron) charge weighting, returns `dA1 = g1/F1`, and
     supplies a twist-3-flavoured `dA2 = c·M x/√Q²` with **c = 0.53 for the
     proton, 0.22 for the deuteron**. This is exactly the routine BeAGLE ships
     as an **empty stub** ("Routine is empty because Pythia is unpolarised").
     Two copies of the same POLRAD-derived engine, one live and one dead, is a
     usable transcription cross-check.
  3. **A live tensor branch.** `radgen_init.f` (verified) selects, for one run
     mode, `plrun = 0.` and **`pnrun = mcSet_PTValue * 2.`** — and `radgen.f`
     keys the tensor path on exactly that: line 1423
     `elseif((plrun.eq.0).and.(abs(pnrun).eq.2))`, line 1425
     `b1 = -3./2.*f1*asym1`, line 1426 `b2 = 2.*aks*b1`. The tensor
     polarization enters the structure-function sum with the POLRAD weight at
     line 793 / 1265 (`if(isf.ge.5) ppol = qn/6`) and line 1574
     (`sfm(1) = f1 + qn/6.*b1`), and the **elastic tensor tail** appears at line
     1514 (`b1 = 2.*tau**2*fm**2`). The target table
     (`radgen_init.f:45-59`) covers n, p, **D**, ³He, ⁴He, ¹⁴N, ²⁰Ne, ⁸⁴Kr —
     **no A = 6 or 7**.
* **What that tensor branch is honestly worth.** It sets
  `b₁ = −(3/2) F₁ A₁ = −(3/2) g₁`, i.e. it computes the *tensor* structure
  function from the *vector* asymmetry. That is not a tensor model; it is a
  placeholder. Its value to LiPolGen is **the kernel, not the content**: it is
  an independent statement of *how* b₁ and b₂ enter the POLRAD radiative
  integrals (the `qn/6` weight, the `sfm` assembly, the elastic
  `2τ²F_m²` tensor tail), which is precisely the structure `include/lipolgen/rc.hpp`
  transcribed from `polrad2t.tex` and the `adgh` FORTRAN. Feed it LiPolGen's
  own b₁ and the RC *kernel* becomes checkable; take its b₁ at face value and
  the check is meaningless.
* **Runnable here.** **Not today.** Needs CERNLIB (`packlib`, `mathlib`,
  `kernlib`), which is absent (`ldconfig` has no `cernlib`); the BeAGLE
  container builds a free CERNLIB from source and is the obvious host. And the
  README states plainly: **"DO NOT TRY radiative corrections. Currently does
  not work"** — the RC path aborts with a Fortran formatted-transfer error at
  `pepsiMaineRHIC_radcorr.v2.F:514`. So the tensor branch is reachable only by
  linking `radgen.f` standalone, not by running PEPSI end to end.
* **Does NOT validate.** Any lithium quantity (no A = 6/7 target); any real
  tensor structure function; the tagged, coherent or breakup channels.

### 2.6 STEG — `JeffersonLab/LightIonEIC` — **the actual tagged-deuteron generator**

> The task brief named "TOPEG" as "the deuteron tagged DIS generator by the
> Cosyn/Sargsian community". **That is wrong** — see §5.2. This is the code
> that brief was looking for.

* **Repository.** `https://github.com/JeffersonLab/LightIonEIC` (verified:
  "Light Ion EIC Code development", **no license**, **last push 2016-09-06**,
  ~15 MB, default branch `master`). Directories: `Event-Generator`,
  `Onshell-Extrapolation`, `Theory-Codes`, `Toolbox`.
* **Provenance.** The JLab Spectator Tagging project page
  (`jlab.org/theory/tag/event`, fetched with a browser user agent — WebFetch
  gets 403, curl succeeds) describes exactly these four categories and says
  *"The codes and simulation results are maintained on github. Scientists
  interested in using these resources … are requested to contact Doug
  Higinbotham or Kijun Park."* The project is the JLab LDRD FY14/15
  "Physics potential of polarized light ions with EIC at JLab".
* **What is in it (verified file listing and sources).**
  * `Event-Generator/SpectatorMC_collMC.cpp` — collider-frame e + D → e′ + p + X
    with crossing angle and beam momentum spread; `SpectatorMC_fixMC.cpp` — the
    fixed-target twin, used as the zero-ion-momentum cross-check (the README
    states the two agree in that limit).
  * `Event-Generator/cweiss_pol/tag_core_pol.f` — **C. Weiss's `TAG` package**,
    header verified: `TAGXP` (tagged cross section, polarized), `TAGGD`
    (tagged structure functions, polarized), `TAGX`/`TAGFD` (unpolarized),
    `TAGSP` (deuteron spectral function); `V2 30JUN14 ADDED POLARIZED E-D`.
    The polarized cross section is `FSIG = FACD*(F2/XD - IPOL*2*D1*G1)` with
    `IPOL = (2·helicity_e)·helicity_D = ±1` and `D1` the depolarization factor.
  * `cweiss_unpol/`, `msargsian_unpol/` — the two unpolarized model families
    (Weiss and Sargsian), selected by `imodel` in `crs_weight_pol.f`.
  * Grids: `std2000_lo_g1.grid` (GRSV2000 g₁), `grv98lo.grid`,
    `JR14NLO08SF.grd`, plus `Deuteron_shadowing_tagged_F2ratio_2014_grid.dat`.
  * Documentation: `Event-Generator/Documenation/STEGs.pdf` (sic).
* **Overlap.** This is the closest external code to LiPolGen's **tagged**
  channel, and it is written by the same author as the theory the tree already
  gates against (`E-1`, Cosyn–Weiss II TABLE II). It could check: the tagged
  cross-section normalisation and its `1/GeV⁴` convention, the spectral-function
  → light-cone (α_R, p_T) mapping, the depolarization factor multiplying g₁,
  and the collider-frame boost with a crossing angle.
* **The limitation that decides its value.** The shipped wave function
  `cweiss_pol/tag_user_wf.f` is a **Hulthén S-wave**, verified in full:
  `A = 0.045647`, `B = 0.2719`, `PSI = (1/(p²+A²) − 1/(p²+B²))/√C`. **There is
  no D-wave.** A deuteron with no D state has no tensor structure at all, so
  this code cannot reach b₁, A_zz, or the D/S ratio η that `C-1` and `E-1`
  live on. `TAGWF` is declared a user-defined routine and is meant to be
  replaced — replacing it with CD-Bonn (which this tree already has,
  `C-6`) is a genuinely small piece of work and would turn STEG into a real
  cross-check of the *vector* tagged sector.
* **Runnable here.** Plausibly. Fortran + C++ with ROOT-era conventions;
  gfortran 11.4.0, g++ 11.4.0 and ROOT 6.28/02 are installed. Unverified
  whether it builds unmodified after ten years; the `.kumac` files imply PAW/HBOOK
  for the histogramming layer, which is a CERNLIB dependency that is **not**
  present here (the cross-section weight program itself writes plain
  `ed_semi_eic**.dat`, so the generator may be usable without PAW).
* **Does NOT validate.** Tensor anything; ⁶Li/⁷Li; α+d clusters; FSI (the
  package is impulse-approximation, `TAGFD` says so in its header); coherent
  diffraction.

### 2.7 eSTARlight — already in tree (`D-4`)

* **Repository.** `https://github.com/eic/estarlight` (verified: **no license**
  reported by the API, last push **2025-05-20**, created 2020-10-09). The
  in-tree run used commit `939b11a…` (2025-05-02).
* **Pedigree.** Lomnitz & Klein, **Phys. Rev. C 99, 015203 (2019)**,
  arXiv:1803.06420 (verified). The EIC Yellow Report says eSTARlight
  "accurately reproduces the essential features of the vector meson production
  and decay" and is "based on parameterized HERA data" — i.e. its pedigree is a
  **fit to HERA γ*p exclusive data**, extended to nuclei by a Glauber
  calculation, not a first-principles nuclear calculation.
* **Overlap.** Already exhausted for the coherent ⁶Li rate, the \|t\| slopes and
  a nine-row Q²-floor scan (`D-4`, `docs/open_items/run_2026-09-02/estarlight_li6.md`).
* **What would *strengthen* `D-4`.** Not more eSTARlight. Two things: (a) a
  **second, independent** dipole model at a species both codes support — Sartre
  cannot do ⁶Li but both do Au/Ca, so a Au cross-check would at least
  test the LiPolGen-side plumbing of the coherent tier; (b) real UPC/HERA data,
  which belongs to the data file, not this one.
* **Does NOT validate.** ⁶Li structure: eSTARlight's Z < 7 branch uses a plain
  **Gaussian** form factor with `R = 1.2·A^{1/3}` (2.18 fm for ⁶Li, well below
  the measured 2.589 fm rms), which is why the in-tree run carried a radius
  patch as a *sensitivity variant*. No clusters, no polarization, no tensor a₂.

### 2.8 STARlight

* **Repository.** `https://github.com/STARlightsim/STARlight` (verified:
  **AGPL-3.0**, last push **2026-04-22**, created 2023-01-06 — the most
  actively maintained code in this whole list).
* **Pedigree.** Klein, Nystrand, Seger, Gorbunov & Butterworth, *Comput. Phys.
  Commun.* **212** (2017) 258, arXiv:1607.03838 (verified). Validated against a
  decade of RHIC and LHC ultra-peripheral-collision measurements — by far the
  best-tested pedigree of any generator here.
* **Overlap.** The **Q² → 0 limit** of the coherent channel: eSTARlight is
  STARlight extended to finite photon virtuality, so STARlight is the natural
  check that eSTARlight's photoproduction limit is unchanged. Also the ρ/φ/J/ψ
  cross-section ratios and the \|t\| slopes at the same nucleus.
* **Runnable here.** Yes in principle — C++ / CMake, ROOT optional.
* **Does NOT validate.** Electroproduction at the EIC's Q² (that is precisely
  what eSTARlight added); lithium beyond the same Gaussian-form-factor default;
  anything polarized.

### 2.9 Sartre 1.39 — **cannot do A = 6, and does not ship A = 2 or 4 either**

* **Repository / obtainability.** `https://sartre.hepforge.org` — **the site is
  behind an Anubis proof-of-work anti-bot challenge and could not be fetched
  from here** (verified: the page returns the challenge shell, not content).
  The version and dependency set are nevertheless verified from `eic-spack`'s
  package recipe: `version("1.39", sha256="82ed7724…")`, url
  `sartre.hepforge.org/downloads/?f=sartre-1.39-src.tgz`, `depends_on` gsl,
  root, `boost@1.39: +thread`, `cuba@4:`, plus two shipped patches
  (`add-cmath-headers.diff`, and a `FindROOT` fix needed only for root ≥ 6.30).
* **Pedigree.** Toll & Ullrich, *Comput. Phys. Commun.* **185** (2014) 1835,
  arXiv:1307.8059 (verified, DOI 10.1016/j.cpc.2014.03.010); bSat/bNonSat
  dipole models. The EIC Yellow Report used **version 1.34** for its coherent
  diffraction projections (arXiv:2103.05419, read locally).
* **The blocking fact.** The code paper, §3.3.4 (read from the arXiv PDF):
  *"At the moment this class is able to describe the following nuclei: proton
  (1), oxygen (16), aluminum (27), calcium (40), copper (63), cadmium (110),
  gold (197), and lead (208)."* **No deuteron, no ⁴He, no A = 6 or 7.** ePIC's
  own Sartre productions are Au and Ca only (`eic/SARTREdataset/runcards`,
  verified: eight run cards, all Au or Ca). Amplitudes live in 3-D lookup
  tables in (Q², W², t) that must be generated **per nuclear species** by a
  separate table generator — adding ⁶Li is a table-generation campaign, not a
  configuration change, and the density it would use is a Woods–Saxon with no
  cluster substructure. **UNVERIFIED:** whether 1.39 added species beyond the
  1.x paper's list, because the download page is unreachable from here.
* **Overlap, if it were run.** An **independent saturation-model** coherent /
  incoherent cross section and \|t\| distribution, plus GEMINI++ statistical
  de-excitation of the struck nucleus — a second opinion on the breakup
  spectrum next to BeAGLE's FLUKA.
* **Runnable here.** Everything but `cuba` is present (gsl and boost 1.74 dev
  headers verified, ROOT 6.28/02 is below the 6.30 patch threshold). `cuba` is
  a small self-contained build. The obstacle is fetching the tarball past the
  anti-bot page, not the build.
* **Does NOT validate.** ⁶Li in any form today; polarization; tagging.

### 2.10 lAger 3.7.0

* **Repository.** `https://eicweb.phy.anl.gov/monte_carlo/lager` (verified via
  the GitLab API: public, created 2020-02-07, **last activity 2025-11-21**;
  latest tag **3.7.0**, 2024-03-27). **GPL-3.0-or-later** (declared in
  `eic-spack`'s `license()` and in the README's copyright block, © 2016-2021
  S. Joosten).
* **Pedigree.** Cited in the EIC Yellow Report for the J/ψ and Υ projections
  (arXiv:2103.05419 §8, read locally). Its process modules are wrappers on
  published calculations: `jpacPhoto_pomeron`, `jpacPhoto_pentaquark`,
  `holographic_vm`, `brodsky_2vmX`, `phi_hatta`, `phi_clas12`,
  `lee_4He_jpsi_grid`.
* **Overlap.** Two narrow but real ones. (a) `src/lager/gen/lA/lee_4He_jpsi_grid.cc`
  is a **⁴He** J/ψ grid — the only light-nucleus exclusive module found in any
  of these codes. (b) `src/lager/physics/fermi.{hh,cc,f}` exposes
  `fermi87(P, A)`, the O'Connell & Lightbody (NBS, 1987) single-nucleon Fermi
  momentum parameterisation for arbitrary A, i.e. an independent n(k) to place
  beside the cluster relative-momentum distributions of `E-13`/`E-14`.
* **Honesty on (b).** `fermi87` is a *single-nucleon* momentum distribution
  from a 1987 fit; LiPolGen's tagged sector needs the **α–d relative** momentum
  from a cluster wave function. Comparing them is a shape sanity check, not a
  benchmark, and must be labelled as such or it is worse than nothing.
* **Runnable here.** Two routes: a prebuilt Singularity container via the
  README's `install.sh` (singularity-ce 3.11.2 is installed), or a local build
  needing ROOT + boost + gsl + HepMC2 **and** HepMC3 + PHOTOS++ + jpacPhoto.
* **Does NOT validate.** Anything polarized or tensor; DIS (it is an
  exclusive/photoproduction generator); ⁶Li.

### 2.11 Exclusive-process generators — surveyed, low relevance

Kept short because none of them touches a LiPolGen observable.

| code | verified location / state | why it is here | why it does not help |
|---|---|---|---|
| **DEMPgen** | `eic/DEMPgen` GPL-3.0, pushed 2025-07-23; `JeffersonLab/DEMPgen` GPL-3.0, pushed 2026-06-19 | deep exclusive meson production for SoLID + EIC; ships `reference_output/` | exclusive π/K on the nucleon; no nucleus, no polarization |
| **MILOU / Milou3d** | `eic/Milou3d` (no license, pushed 2021-05-23); spack `milou` → `gitlab.com/eic/mceg/milou`, needs cernlib + pythia6 | DVCS with GPD models; YR §8 used MILOU 3D | DVCS on the proton |
| **HEPGEN** | arXiv:1207.0333 (Sandacz & Sznajder, verified; **no journal ref in the arXiv record**). HEPGen++ described on the COMPASS TGEANT pages. **No public repository located.** | hard exclusive photon/meson leptoproduction, COMPASS kinematics | not distributed publicly; nucleon target |
| **TOPEG** | see §5.2 — coherent DVCS off ⁴He, no public repository located | the only nuclear exclusive generator in the YR | ⁴He, DVCS, unpolarized |
| **ELRADGEN 2.0** | arXiv:1104.0039 (Akushevich, Filoti, Ilyichev, Shumeiko; verified, **no journal ref in the arXiv record**) | same authors and same formalism family as POLRAD, i.e. a sibling of this tree's RC engine | **elastic** ep only, spin-½ target; no DIS tail, no spin-1 form factors |

### 2.12 eic-smear — the format bridge

* **Repository.** `https://github.com/eic/eic-smear` (verified: **GPL-3.0**,
  last push **2026-07-07**; spack versions up to 1.2.2). In `eic-shell`'s
  default environment.
* **Overlap.** `src/erhic/` (verified listing) contains readers for
  **Djangoh, Pepsi, Dpmjet (= BeAGLE), Sartre, Milou, DEMP, Pythia6, Rapgap,
  GmcTrans, HepMC, Simple**, and `TreeToHepMC.cxx` writes HepMC3. That means
  every generator in this file lands in the same HepMC3 shape LiPolGen writes,
  and the ePIC production chain (`eic/djangoh_production/versions`, verified)
  is literally *generator → eic-smear `BuildTree` → `TreeToHepMC` → HepMC3
  Asciiv3 → `abconv` → `hepmc3.tree.root` → `npsim`*.
* **Value to LiPolGen.** LiPolGen writes HepMC3 **directly**, skipping the
  eic-smear leg. That is a design advantage, but it also means the tree's `F-1`
  / `F-2` chain evidence is the *only* evidence that its record is shaped like
  everyone else's. Running one external generator through the standard leg and
  diffing the resulting HepMC3 against LiPolGen's convention
  (`docs/HEPMC3_CONVENTION.md`) is a cheap, decisive chain check.
* **Does NOT validate.** Any physics whatsoever.

### 2.13 The ePIC full chain

* **Verified state (2026-09-06).** `eic/epic` LGPL-3.0, pushed 2026-09-05;
  `eic/npsim` pushed 2026-09-05; `eic/EICrecon` LGPL-3.0, pushed 2026-09-05;
  `eic/afterburner` (the `abconv` crossing-angle/beam-effects tool) **no
  license**, pushed 2026-07-05; `eic/eic-shell` is the container entry point;
  `eic/genpythia` (GPL-3.0, pushed 2026-09-02) is the reference "PYTHIA 8 →
  `hepmc3.tree.root` → npsim" app.
* **Relation to the tree.** This *is* `F-1` / `F-2`, which passed once on
  2026-09-02 and is not re-run by any test invocation. Nothing new is needed
  here except the container: **there is no `/cvmfs` and no `eic-shell` on this
  machine**, so re-running the chain means pulling `eic_xl` with the installed
  `singularity-ce 3.11.2`.
* **What genpythia would add.** A same-machine reference `hepmc3.tree.root`
  from a *known-good* generator, so that a future `npsim` failure can be
  attributed to LiPolGen's record rather than to the container.

---

## 3. Runnable on this machine, today

| dependency | state on this machine |
|---|---|
| gfortran / g++ | 11.4.0 / 11.4.0 (Ubuntu 22.04) |
| cmake | 3.22.1 |
| ROOT | 6.28/02, cxx17 (below the 6.30 threshold of Sartre's spack patch) |
| GSL, Boost | headers present (`/usr/include/gsl/`, Boost 1.74) |
| LHAPDF, HepMC3, PYTHIA 8 | 6.5.5 / 3.3.0 / 8.317 in `/home/cpeng/Projects/polli/deps/install` |
| CERNLIB | **absent** — blocks PEPSI, MILOU, DJANGOH's `ranlux`, and STEG's PAW layer |
| Cuba | **absent** — blocks Sartre |
| LHAPDF 5 | **absent** — blocks DJANGOH and BeAGLE |
| FLUKA | **absent, and licensed** — blocks BeAGLE |
| Singularity | `singularity-ce 3.11.2-jammy` present |
| `/cvmfs`, `eic-shell` | **absent** |

**Verdict by code.** STARlight, eSTARlight and lAger (container route) are
buildable/runnable today with no new licensed dependency. Sartre needs `cuba`
plus a way past hepforge's anti-bot page. BeAGLE needs one licensed FLUKA
tarball and then the `eic/beagle-container` recipe does everything else.
PEPSI needs CERNLIB (free, buildable — the BeAGLE container already does it).
DJANGOH needs the source, which is not publicly posted. STEG needs nothing
obvious but is ten years stale and unverified against a modern toolchain.

### 3.1 Live check run for this survey

The one in-tree generator-comparison gate (`D-1`) was re-run read-only on the
survey date:

```
$ source env.sh && ./build/lipolgen_tests -tc="pythia: the charged multiplicity*"
[pythia] charged multiplicity, e + p at 10 x 99.5 GeV/u, Q2 in [4.0, 30.0], x in [0.005, 0.10]
         stock WeakBosonExchange : 8.091   (1153 events, <x> = 0.0312, <Q2> = 8.88)
         LiPolGen PythiaBridge   : 7.851   (2000 events at that exact (x, Q2))
         ratio bridge/stock      : 0.970
[doctest] test cases:  1 |  1 passed | 0 failed | 402 skipped
[doctest] assertions: 25 | 25 passed | 0 failed |
```

So `D-1` is measurable and passing in this checkout; the ratio quoted in
`docs/PYTHIA_BRIDGE.md` §10 reproduces.

---

## 4. Where an external generator would strengthen or replace an in-tree check

| in-tree row | its present limit | external generator that would move it | what would actually change |
|---|---|---|---|
| `D-1` PYTHIA closure | one observable (charged multiplicity), one (x, Q²) point, ±20 % window | the same PYTHIA 8 | broaden, don't replace: add HFS energy/pT flow and a second kinematic point. No new code needed. |
| `D-3` POLRAD elastic tail | magnitudes gated 0.25×–4× because the tree uses a different deuteron form factor than POLRAD's `ffdeu` | **DJANGOH/HERACLES** | a *second, independently written* one-loop RC engine. Even the published Rad=1/Rad=0 table (§2.4) constrains the RC magnitude without running anything. |
| `D-3` / `C-8` RC band δ(x) | anchored on one HERMES residual and a must-not-contradict bound | DJANGOH `Rad=1` vs `Rad=0` ratios (§2.4 table); PEPSI/RADGEN for the polarized case | turns a one-point bound into a Q²-trend comparison |
| the tensor RC band | pinned only against POLRAD's own paper and FORTRAN | **PEPSI's `radgen.f`** `pnrun = ±2` branch | checks the *kernel* — `qn/6` weighting, `sfm` assembly, elastic `2τ²F_m²` tail — against a second copy of the same engine, with LiPolGen's own b₁ substituted for RADGEN's placeholder |
| `E-1` CW TABLE II tagged tensor | a **paper table**, not running code | **STEG `cweiss_pol/tag_core_pol.f`** — same author | the tagged cross-section normalisation and the spectral-function → (α_R, p_T) map become executable. Requires replacing the Hulthén S-wave `TAGWF` with CD-Bonn (already in the tree, `C-6`). Reaches the **vector** sector only. |
| `D-4` eSTARlight ⁶Li baseline | a single code, with a Gaussian ⁶Li form factor and `R = 1.2A^{1/3}` | STARlight (Q²→0 limit); Sartre **at Au/Ca only** | a photoproduction-limit consistency check and a dipole-model second opinion at a species both support. **Not** a ⁶Li cross-check — no second code does ⁶Li. |
| `E-16` far-forward efficiencies | a published table (Chang *et al.*), with ⁷Li standing in for ⁶Li | **BeAGLE** via the container | fragment η and rigidity distributions generated rather than read off a table — though for A ≥ 5 BeAGLE's own nuclear input is a generic Woods–Saxon with a known Fermi-distribution bug |
| `E-15` Glauber FSI | analytic limits and a Python prototype (`B-2`) | **none** — see §5.4 | nothing today |
| `F-1` / `F-2` chain | passed once, not re-run, needs external containers | `eic/genpythia` + `eic-shell`, or one external generator through `eic-smear TreeToHepMC` | a same-machine known-good reference file, so a future `npsim` failure is attributable |
| `F-4` HepMC3 spin attributes | a proposed convention with no ePIC standard | none | no generator writes target-spin attributes; this stays a proposal |

---

## 5. Negative results — checked, and they do not exist

These are the most useful part of this file, because each one closes a line of
enquiry that looks promising from the outside.

### 5.1 No generator implements a tensor-polarized target

Searched: the PYTHIA 8.317 source tree (already inventoried —
`SigmaProcess` has zero helicity/spin/polarization hits and there is no
polarized PDF machinery), BeAGLE's own sources (zero hits for
`polariz|helicit|spin|g1|b1` outside RADGEN), DJANGOH's 4.6.10 update note
(polarization is **longitudinal nucleon** only, via free-proton polarized PDFs),
Sartre (dipole models, unpolarized), eSTARlight / STARlight (unpolarized),
lAger (unpolarized), and a literature search for a tensor-polarized MCEG. The
only tensor code found anywhere is the RADGEN `b₁`/`b₂` branch of §2.5, present
in **two** copies (BeAGLE's `radgen-6.4.28`, PEPSI's `radgen.f`), which sets
`b₁ = −(3/2) F₁ A₁` from the vector asymmetry and is therefore a placeholder.
**LiPolGen has no peer implementation to diff against, and this file cannot
manufacture one.**

### 5.2 TOPEG is not a tagged deuteron generator

The task brief described TOPEG as "the deuteron tagged DIS generator by the
Cosyn/Sargsian community — the most relevant cross-check for the tagged
channel". It is not. Verified from the EIC Yellow Report's own text
(arXiv:2103.05419, read locally from `PolarizedLithiumSim/refs/`, part 2):

> *"The Orsay-Perugia Event Generator (TOPEG). For this study, a new Monte
> Carlo event generator for the coherent DVCS off the ⁴He nucleus has been
> developed. This tool, called the Orsay-Perugia event generator (TOPEG), is
> based on the Foam ROOT library … Checks at the JLab kinematics with an
> electron beam energy of 6 GeV have been successfully performed."*

So: **coherent DVCS off ⁴He**, an unpolarized exclusive generator, no tagging,
no deuteron, no spin. No public repository was located. The code the brief was
reaching for is STEG / `JeffersonLab/LightIonEIC` (§2.6).

### 5.3 Sartre has no light nucleus

Sartre's `Nucleus` class supports p, O, Al, Ca, Cu, Cd, Au, Pb (§2.9, quoted
from the code paper). There is no deuteron, no ⁴He, and no A = 6 or 7, and
adding one is a lookup-table generation campaign over (Q², W², t) using a
Woods–Saxon density that could not represent α+d anyway. Sartre is **not** a
route to a ⁶Li coherent cross-check.

### 5.4 No generator computes cluster-spectator Glauber FSI

BeAGLE's INC is switched on only for A > 12 for SRC pairs and is a
formation-zone cascade, not an eikonal survival factor; Sartre's Glauber is for
the dipole–nucleus amplitude; Angantyr's is wounded-nucleon counting. LiPolGen's
`fsi.hpp` computes something none of them computes. `E-15` (analytic limits +
published formulae) remains the only available check.

### 5.5 No public GCF event generator was located

BeAGLE has a `MODEL GCF-FT` mode (`src/dpm_gcf.f`) that **reads** an external
Generalized-Contact-Formalism event file in BeAGLE's own ASCII format and runs
INC + FLUKA evaporation on it — but the GCF generator that writes those files
was not found in any public repository. Obtainability: **on request** from the
SRC/GCF community. Relevance to LiPolGen is secondary anyway (SRC pairs, not
α+d clusters).

### 5.6 DJANGOH's source is not publicly posted

The author's distribution page is 403 from here, there is no mirror on GitHub,
and it is not an `eic-spack` package nor in the `eic-shell` environment — even
though ePIC produces with it routinely (4.6.10 and 4.6.21). The **manuals** and
the **steering cards** and one **cross-section table** are public (§2.4); the
code is obtained by contacting the author or the ePIC production team.

---

## 6. Recommended order, with the honest payoff of each

1. **DJANGOH's published Rad=1 / Rad=0 table** (§2.4). Zero code, zero
   dependencies, obtainable this minute. Constrains the *magnitude and Q²
   trend* of the unpolarized RC against an independent one-loop engine.
   Payoff: turns `C-8`'s single must-not-contradict bound into a trend check.
2. **PEPSI's `radgen.f`, linked standalone** (§2.5). Needs CERNLIB only for the
   full PEPSI; the RC kernel can be exercised on its own. Payoff: the only
   external statement anywhere of *how* b₁ and b₂ enter the POLRAD radiative
   integrals — a genuine check of the tree's transcription, with LiPolGen's own
   b₁ substituted for RADGEN's placeholder. Cost: must be written up so that
   nobody mistakes RADGEN's `b₁ = −(3/2) g₁` for a tensor model.
3. **STEG `cweiss_pol` with `TAGWF` replaced by CD-Bonn** (§2.6). Payoff: makes
   `E-1` executable rather than a paper table, from the same author. Cost: a
   ten-year-old unlicensed tree, and it reaches the **vector** sector only.
4. **BeAGLE in the container, eD mode** (§2.3). Payoff: independent tagged
   spectator momentum spectra and far-forward fragment kinematics at A = 2, and
   an evaporation model to place beside the T1 breakup. Cost: one licensed
   FLUKA tarball; and every A ≥ 5 number it produces carries the Woods–Saxon
   and Fermi-distribution caveats of §2.3.
5. **`eic/genpythia` + `eic-shell` re-run of `F-1`/`F-2`** (§2.13). Payoff: a
   same-machine known-good chain reference. Cost: a container pull.
6. **DJANGOH proper, if the source can be obtained** (§2.4). Payoff: the only
   maintained code that generates *polarized* DIS with full radiative
   corrections and prints `G1NC`, `A1NC` and `D` as columns. Cost: contact the
   author or ePIC production; LHAPDF 5 + CERNLIB.
7. **STARlight for the Q² → 0 limit of `D-4`** (§2.8). Cheap, AGPL, actively
   maintained, but a narrow check.
8. Everything else in §2.11 and §5: record and move on.

**Nothing on this list validates a tensor observable, a lithium observable, or
a cluster wave function.** The list is worth doing anyway, because it converts
several in-tree rows from "our transcription of a paper" into "our code against
someone else's running code" — but it must never be summarised as validating
the generator's purpose.
