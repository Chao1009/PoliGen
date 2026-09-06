<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# 00 — Cross-checks that already exist in this tree

**Survey date 2026-09-06.** Baseline inventory for the benchmarking series: every
comparison the repository *already* performs, and precisely what each one does and
does not establish. Nothing here is a proposal; everything is code that runs today.
Later files in `docs/benchmarking/` are measured against this one.

## 0. The frame this survey has to be read in

Three facts bound every claim below.

1. **There is no polarized e + ⁶Li or e + ⁷Li DIS data, anywhere.** Not inclusive,
   not tagged, not coherent. There never has been. No number in this repository
   can be validated against a measurement of the process it generates.
2. **No other generator does tensor-polarized lithium.** BeAGLE, eSTARlight,
   PYTHIA, DJANGOH, Sartre and the ePIC chain all carry unpolarized (or at best
   vector-polarized) nuclei. There is no second implementation to diff against.
3. **The only tensor-polarized DIS datum in the world is HERMES's deuteron b₁**
   (Airapetian *et al.*, PRL **95** (2005) 242001, `hep-ex/0506018`) — six bins,
   one nucleus, A = 2, and it is **not** pinned as an agreement test anywhere in
   this tree (§3, row C-7).

So the value of this survey is the *decomposition*: what can be checked at the
**nucleon**, **deuteron**, **unpolarized-nucleus**, **numerical-identity** and
**software-chain** levels, and where each of those checks stops.

### The five verdict classes used below

| class | meaning |
|---|---|
| **DATA** | compared against a measured number (with its published uncertainty) |
| **GENERATOR** | compared against another, independently maintained code |
| **THEORY** | compared against someone else's published calculation (table, closed form, or a figure digitized from their PDF) |
| **CHAIN** | a software-interface check: the record is well formed and downstream tools accept it |
| **INTERNAL** | closure against this project's own earlier implementation or its own frozen output — a **regression lock**, never a benchmark |

`polligen` (`PolarizedLithiumSim/evgen/polligen`, `fastsim/polli_fastsim`) is
**this project's own NumPy predecessor**. Every `validation/reference/*.json`
comparison is therefore **INTERNAL**: it proves the C++ port did not change the
numbers, and it proves nothing whatsoever about whether those numbers are right.
It is listed first, and at length, precisely so that it is never mistaken for
external validation.

---

## 1. Table A — INTERNAL: the `polligen` reference tables (rtol 1e-12)

Generator: `validation/dump_polligen_reference.py` (imports `polligen` /
`polli_fastsim` read-only; never writes into the Python repository).
Provenance and per-table caveats: `validation/README.md`; inventory and
"functions this script could not call" list: `validation/reference/_manifest.json`.
Every file is full-double-precision `repr()` JSON, so the C++ side (`jsonmin`,
`strtod`) reads back identical bits.

| # | reference file | what it freezes | consumer TEST_CASE(s) | tolerance | class |
|---|---|---|---|---|---|
| A-1 | `spin.json` | `spin.wigner_d`, `clebsch_gordan`, ρ(lab) for J = 1 and 3/2 on five axes, vector/tensor/octupole moments, `SpinDensity.lab_moments` | `tests/test_reference.cpp` "reference tables: spin.json"; `python/tests/test_kernel_reference.py` | rtol **1e-12**; max-entropy fills **1e-13** (bisection tol 1e-13 upstream) | INTERNAL |
| A-2 | `xsec.json` | `InclusiveKernel.tables/amplitudes/dsigma`, `asymmetries.*`, the exact finite-γ tensor block (`tensor_gamma`), the `target_mass` block, the J = 3/2 rank-2 scenario | `tests/test_reference.cpp` "reference tables: xsec.json" | rtol **1e-12** (incl. `g2_ww`, which must use the *same* 96-point trapezoid) | INTERNAL |
| A-3 | `beams.json` | `beams.IONS`, `Ion.mass_per_nucleon`, `momentum_per_nucleon_max`, `default_configs`, `√s/nucleon`, `gamma_of` | `tests/test_reference.cpp` "reference tables: beams.json" | rtol **1e-12** | INTERNAL |
| A-4 | `tagged.json` | `tagged.Wave.radial` on a 40-point k grid, `TaggedModel` grid cells (`n_of_kc`, `struck_populations`, `population_integrated`, dilutions, `p2_moment`), `boost_spectator` | `tests/test_tagged.cpp` "tagged: channel construction against polligen", "…the model grid and its tables against polligen", "…boost_spectator against polligen"; `tests/test_cluster.cpp` "cluster: `Wave::radial` reproduces the polligen radial tables" | closed forms **1e-12**; `norm`, `p2_moment_mixture` are grid quadratures → **1e-9 / 1e-6** by design | INTERNAL |
| A-5 | `spectator.json` | `M_U`, `MASSES`, `NUCLEUS_MASS` (AME2020-derived), every `ClusterChannel`'s masses/κ/R(k=0), unnormalized `momentum_density`, `_boost_fragment` at fixed (kₓ,k_y,k_z) | `tests/test_spectator.cpp` (5 cases) | rtol **1e-12** | INTERNAL |
| A-6 | `coherent.json` | `CoherentScenario` defaults and formulas (`coherent_fraction`, `tag_acceptance`, `mean_t_tagged`, `a2_deformation`, `cos2phi_*`), `recoil_lab`, `fragment_rigidity`, the ⁶Li/⁷Li breakup tables, `MANTYSAARI_A2_DEUTERON` | `tests/test_coherent.cpp` "coherent: constants and scenario defaults against polligen", "…the scenario formulas against polligen", "…fragment rigidities are mass-to-charge ratios" | rtol **1e-12** | INTERNAL |
| A-7 | `bookkeeping.json` | the four run-plan constructors, `SpinCategory.moments()`, `lumi_shares`, the t = 3 spin-temperature ladder, `azz_rel_lumi_bias`, `apar_rel_lumi_bias` | `tests/test_bookkeeping.cpp` "reference tables: bookkeeping.json plans", "…spin-temperature ladder and rel-lumi biases" | **1e-12**; maxent anchors **1e-13** | INTERNAL |

**What Table A validates:** that the C++ reimplementation of ~40 numerical
routines is bit-compatible with the NumPy original, including NumPy's pairwise
summation order and `np.interp`'s clamped linear interpolation
(`tests/test_sf.cpp`: "np_interp reproduces numpy's clamped linear interpolation",
"the NumPy pairwise summation order is reproduced" — pinned to NumPy's own
`6.2826638802995038`).

**What Table A does NOT validate:** any physics. A shared bug in
polligen and LiPolGen passes every one of these at 1e-12. Three quantities are
explicitly *not* at 1e-12 and say so — `TaggedModel.norm`,
`p2_moment_mixture`, `spin.populations_maxent` — because upstream they are
themselves iterative or grid quadratures.

**Skip behaviour:** A-1/2/3 and the `tagged.json`/`spectator.json` consumers
print "…not found — skipping" and pass with zero assertions if
`LIPOLGEN_REFERENCE_DIR` is empty; the `coherent.json` and `cluster.cpp`
consumers use `doctest::skip(...)` and are reported as SKIPPED instead.

---

## 2. Table B — INTERNAL: LiPolGen's own frozen output

| # | artefact | what it freezes | consumer | tolerance | class |
|---|---|---|---|---|---|
| B-1 | `validation/reference/b1_default_li6.json` | `default_inclusive_kernel(LI6())`'s `f1/b1/b2/delta` on a 200-point log-x grid at Q² = 1, 2.5, 10 — the **library default** `B1Model::Miller` path | `tests/test_b1_nuclear.cpp` T9 "default_inclusive_kernel(LI6()) is unchanged" | rtol **1e-12** | INTERNAL (self-pin) |
| B-2 | `docs/open_items/prototypes/fsi_alpha.py` table | Glauber X–α profile (σ_cluster = 131.0 mb at σ_XN = 40 mb, σ_el = 35.2 mb, B = 27.2 GeV⁻²) and 7 FSI/IA rows | `tests/test_fsi.cpp` "fsi: the prototype FSI/IA table and profile are reproduced" | **5e-3** absolute on the ratio rows; ≤ 0.1 mb on the profile | INTERNAL (Python prototype → C++ port) |
| B-3 | measured C++ literals throughout `tests/` | several hundred `CHECK_CLOSE` pins on quantities this build itself computed (e.g. `T22b`'s η = −0.0482160903696 at 1e-9, `T24`'s VMC MC band 0.4705 %/1.6176 %, `T10b`'s q_charge/Q_TUNL = 11.42) | the owning TEST_CASE | typically 1e-9 … 1e-3 | INTERNAL (regression anchors, explicitly labelled as such in the comments) |

`dump_b1_default_li6.py`'s own docstring is the clearest statement of why B-1 is
not a benchmark: it "runs through the INSTALLED pybind11 module (`import
lipolgen`), i.e. the same C++ `default_inclusive_kernel` the test calls, so the
file records what the library does TODAY."

---

## 3. Table C — DATA: comparisons against measured numbers

Every row here is A ≤ 6 and none of them is the polarized-lithium DIS observable
the generator exists to produce.

| # | measured quantity + source | LiPolGen symbol | TEST_CASE | tolerance / verdict | validates | does **not** validate |
|---|---|---|---|---|---|---|
| C-1 | **η(⁶Li → α+d) = −0.025 ± 0.006(stat) ± 0.010(syst)** — E. A. George, L. D. Knutson, PRC **59** (1999) 598, from d+α elastic tensor analysing powers | `ClusterConfigSampler::asymptotic_ds_ratio()`, `LI6_ETA_DS_GK{,_STAT,_SYST}` | `tests/test_cluster_config.cpp` T22 "the asymptotic D/S ratio, Whittaker-divided"; T22b "eta rides the quadrupole dial" | model gives **−0.0482**, i.e. **1.93×** the measurement. The test asserts `1.2 < |η/η_GK| < 3.0` — it **gates a documented discrepancy, not agreement** | that the α–d D/S asymptotic ratio is computed correctly (Whittaker-divided, not the naive R₂/R₀ = −0.149) and that the model's D wave is a *moderate* — not order-of-magnitude — excess | anything about the ⁶Li tensor observables themselves. T22b shows η is **exactly linear** in the quadrupole dial, so it is a relabelling of that dial, and GK's ±1σ band maps onto a Q(⁶Li) interval that **crosses zero** |
| C-2 | **Q(⁶Li) = −0.0818(17) fm²** — TUNL A = 6 evaluation, Tilley *et al.*, NPA **708** (2002) 3 (1998CE04) | `LI6_QUADRUPOLE_FM2`; `HoSpin1FF::for_ion(LI6()).fq(0)` | `tests/test_rc.cpp` T11 "the 6Li form-factor anchors"; `tests/test_cluster_config.cpp` T22b | F_q(0) = −65.914 at 1e-4, **derived** from the constant (never retyped); explicitly checked to be TUNL's value and not Pyykkö's −0.0806 (which would give −64.947) | that the spin-1 form-factor normalisation is anchored on the measured quadrupole | the ⁶Li wave function: the same T22b prints the cluster model at **7.52× the measured Q**, and the GFMC value (Pastore *et al.* −0.20(6) fm²) at 3.08× — which is why §10's quadrupole is a **dial with a band**, not a model output |
| C-3 | **μ(⁶Li) = +0.8220473 μ_N** — TUNL A = 6 | `LI6_MU_N`; `fm(0)` = (m_A/m_p)·μ | `tests/test_rc.cpp` T11 | F_m(0) = 4.90765 at 1e-4, derived | the magnetic form-factor normalisation | the q-dependence: the C0 shape and its first zero are an unfitted **model band** ([2.9, 3.3] fm⁻¹), with **no elastic ⁶Li data in this repository** to close it (the test says so in full) |
| C-4 | **AME2020 masses and separation energies**: S_d = 2.2246 MeV, S_n(³H) = 6.2572, S(p+n+n) = 8.4818 | `nuclear_mass(Z,A)`, `cluster_mass`, `ClusterChannel::separation_energy` | `tests/test_breakup.cpp` "species table and the AME2020 separation energies"; `tests/test_spectator.cpp` "the AME2020 mass tables are the Python's"; `tests/test_beams.cpp` "mass per nucleon uses the physical nuclear mass" | 1e-3 MeV on the S values; 1e-12 (rtol) on the table itself | that separation energies are *derived* from one mass table rather than hard-coded twice, and that binding is carried at the 0.5 % level that separates a 41 GeV proton from a γ-matched 40.8 GeV/u ⁶Li | nothing dynamical |
| C-5 | **TUNL A = 6 breakup thresholds / levels** (2.186 MeV 3⁺;0, 5.366 MeV 2⁺;1, Q_m = 1.4743 MeV, the 3.5629 MeV 0⁺ α+d channel) | `coherent::LI6_BREAKUP`, `LI7_BREAKUP`, `fragment_rigidity` | `tests/test_coherent.cpp` "coherent: the veto table routing", "…fragment rigidities are mass-to-charge ratios" (against `coherent.json`) | exact table membership; rigidities at rtol 1e-12 | that the incoherent-breakup veto table's species list and rigidity ratios are the evaluated ones | branching fractions or any breakup dynamics — the table is a *routing* table |
| C-6 | **CD-Bonn deuteron properties**: P_D = 4.85 %, A_S = 0.8846(9) fm^{−1/2}, η = 0.0256(4), Q_d = 0.270 fm² — Machleidt, PRC **63** (2001) 024001, Tables XV and XX | `CdBonnWave` (20 typed coefficients + 4 derived), `cdbonn_fdeut_table()`, `alpha_d_quadrupole_fm2` | `tests/test_cluster.cpp` "the CD-Bonn deuteron coefficients (Machleidt Table XX)"; `tests/test_b1_nuclear.cpp` checklist item 5 | published-value clauses: P_D 1.5e-4, A_S 2e-4, η 5e-5, Q_d 5e-4 abs. **Regression pins at 1e-9** are what actually catch a mistyped digit; the file *measures* the sharpness (one unit in the last printed digit of each coefficient moves the norm by 0.2–243 ×10⁻⁹) and names the **one** digit (C₃'s 8th figure) that is unresolvable | a **transcription gate** on twenty published numbers, plus partial contact with the *empirical* A_S and η (the parenthesised errors on those two are the experimental ones); and the S–D relative-sign convention, settled against CD-Bonn's own Q_d | the deuteron wave function as physics: the parameterisation is a **fit** to Machleidt's numerical solution (his quoted L2 quality 2.2e-4 / 1.1e-4), so P_D = 4.8562 % does **not** round to the published 4.85 and the test says so rather than loosening a tolerance |
| C-7 | **HERMES b₁ᵈ(x = 0.012, Q² = 2.5) = (11.20 ± 5.51 ± 2.77)×10⁻²** — Airapetian *et al.*, PRL **95** (2005) 242001 | `toy_b1` on `tables::kB1Miller` | `tests/test_sf.cpp` "the digitized Miller b1 curve reproduces the published figure" | **NOT ASSERTED.** The two `CHECK`s are digitization checks (0.11429317074113018 at rtol 1e-12; 0.114 at 1 % absolute). The comparison to the datum — 0.1143/2 = 0.057 per nucleon vs 0.112 ± 0.055 ± 0.028, **0.97 σ low** — lives only in the comment | that reading the digitized Miller column at x = 0.012 returns the published curve value | **agreement with HERMES.** The comment says so explicitly. The world's only tensor-DIS measurement is present in this tree as prose, not as a gate |
| C-8 | **HERMES's achieved radiative residual, 15 % at x = 0.063** (same reference) | `rc_delta(x)` interpolation | `tests/test_rc.cpp` T3 "the RC band delta(x) is monotone, clamped and anchored" | `CHECK(rc_delta(0.063) > 0.10)` — a **must-not-contradict** bound, not a fit | that the RC band's x-interpolation cannot claim to beat a published residual | the RC calculation itself (see D-3) |
| C-9 | **α–p σ_tot = 121.5 mb, σ_el = 31.4 mb → σ_el/σ_tot = 0.258** — Blinov *et al.*, Phys. At. Nucl. **64** (2001) 907, `nucl-ex/9910012` | `GlauberFsiWeight::sigma_cluster_{,el}_mb` | `tests/test_fsi.cpp` "fsi: the prototype FSI/IA table and profile are reproduced" | Glauber gives **0.269**; the assertion pins 0.269 at 2e-3 and the 0.258 appears in the comment — **reported, not enforced** | that the elastic/total split of the composite X–α profile is an automatic consequence of the Glauber construction and lands within 4 % of the measured α–p ratio | the FSI weight itself (there is no tagged-spectator FSI datum for ⁶Li) |
| C-10 | **k_F(⁶Li) = 0.169 GeV** — Moniz *et al.*, PRL **26** (1971) 445 (quasi-elastic e–A scaling fit) | `RC_QE_KF_GEV` | `tests/test_rc.cpp` Q9 "the QUASI-ELASTIC Pauli suppression S(q)" | used as an **input**; the test scans 0.169 vs 0.221 and reports the tail ratio (0.472 / 0.868 / 1.000 at x = 0.01 / 0.1 / 0.3) at 2e-3 | that the de Forest–Walecka Fermi-gas suppression is implemented with a measured k_F and that turning it off is a 15–60 % effect on the dominant tail piece | the suppression factor against any quasi-elastic ⁶Li measurement |
| C-11 | **³He(e,e′p)d spectroscopic factor S₀ ≈ 2/3** | `cs_norm(cs_n0_terms(3))` | `tests/test_triton_sf.cpp` "the CS spectroscopic factors" | 0.6525 at 1e-3; S₀+S₁ = 1 to 5e-4 **untuned** | that the Ciofi degli Atti–Simula A = 3 parameterisation is transcribed correctly and closes its own sum rule without fitting | the triton breakup dynamics; the ⁷Li → α + t remnant model remains "crude, flagged" |
| C-12 | **EIC machine menu**: 41/100/275 GeV proton ring, γ-matched light-ion energies, rigidity caps (d 137.5, ³He 183.3, ⁶Li 137.5, ⁷Li 117.9 GeV/u), EPIOS γ-synchronisation windows (`arXiv:2510.10794`) | `default_configs`, `Ion::momentum_per_nucleon_max`, `gamma_of`, `epios_window_of` | `tests/test_beams.cpp` (6 cases) | exact / 0.15 GeV; the 100 GeV point (γ = 106.6) is checked to sit in **neither** EPIOS window — a **known conflict flagged, not resolved** | that the beam menu is the published one and that the γ-matching prescription is applied consistently | nothing physics-side |

**Nothing in Table C is polarized DIS off lithium.** The closest measured contact
with the tensor sector anywhere in the tree is C-7, and it is a comment.

---

## 4. Table D — GENERATOR / external-code comparisons

| # | external code | what is compared | TEST_CASE / artefact | tolerance / verdict | validates | does **not** validate |
|---|---|---|---|---|---|---|
| D-1 | **PYTHIA 8.317, stock `WeakBosonExchange` run** at `gen_dis_hfs.py`'s own settings | mean charged multiplicity of the hadronic final state in 4 < Q² < 30 GeV², 0.005 < x < 0.10 | `tests/test_pythia.cpp` "pythia: the charged multiplicity agrees with a stock WeakBosonExchange run in the same (x, Q2) window" | ratio bridge/stock ∈ **[0.80, 1.20]**, measured ≈ 0.97 (`docs/PYTHIA_BRIDGE.md` §10) | that the `LHAup` bridge's **shower + hadronization** is the same physics as a stock PYTHIA DIS run at the same point | the hard process or any cross section. `DEVELOPMENT_PLAN.md` §4 row 7 was **rewritten 2026-09-05** to say so: under LHAup strategy 3, `info.sigmaGen()` is not a physical cross section, so the original "PYTHIA cross section matches `gen_dis_hfs.py`" clause is **undemonstrable in this architecture** and was withdrawn |
| D-2 | **PYTHIA's own MSTW2008 LO grid** (`mstw2008lo.00.dat`, via `Pythia8::MSTWpdf`) | `MstwSF::f2p/f2n` against a raw `MSTWpdf` probe; and MSTW/CT18NLO F2p ratios | `tests/test_mstw_sf.cpp` (4 cases) | **1e-14** against the raw probe (same grid, same arithmetic order); the MSTW/CT18 ratio 0.932/0.885/0.984/1.168/1.418/1.556 at x = 0.01…0.80 | that the wrapper is bit-for-bit the library it wraps, so the b₁ gate's PDF dependence is the *grid* and nothing else | the PDF itself |
| D-3 | **POLRAD 2.0** (Akushevich *et al.*, CPC **104** (1997) 201), compiled standalone with **its own** `ffdeu` | σ^el_U at three (E, x, y) points; the sign of σ^el_T/σ^el_U; G_C(0)=1, G_M(0)=1.7139610634, G_Q(0)=25.84 | `tests/test_rc.cpp` T8(c′) "POLRAD's OWN compiled deuteron numbers"; T10 "the DEUTERON fixes the (F_c, F_m, F_q) normalisation" | magnitudes gated **loosely** (0.25×–4×) because this repository uses a *different* deuteron form factor (an HO stand-in, not `ffdeu`); measured 0.883 / 0.983 / 0.390 of POLRAD's. The two low-x points are tightened to **20 %**. The **sign** of σ_q/σ_u — which changes between x = 0.05 and 0.20 — is gated **tightly** | that Eq. (38) is transcribed with the right power of A, the right prefactor, no missing Y₊ and no dropped leading minus (any of those moves the magnitude ≥ 6×); and that G_M(0) matches POLRAD's compiled value to 1e-7 | the radiative correction as physics for ⁶Li. There is **no** tensor-sector RC calculation to compare against except Gakh–Shekhovtsova (`hep-ph/0403262`), which has **zero INSPIRE citations** and supplies only the band half-width `RC_DELTA_LOW_X` (see §7) |
| D-4 | **eSTARlight** (Lomnitz & Klein, PRC **99** (2019) 015203) — run in-house 2026-09-02 at commit `939b11a…`, 2×10⁵ events/channel, e 10 GeV × ⁶Li 99.5 GeV/u | coherent J/ψ, φ, ρ cross sections (3 Q² windows × 2 densities), dN/d|t| slopes, and a Q²-floor scan (9 rows) | `include/lipolgen/coherent.hpp` `estarlight_li6_coherent()` / `estarlight_li6_q2_floors()`; `tests/test_coherent.cpp` "the eSTARlight 6Li baseline, and O5's reach arithmetic", "the eSTARlight 6Li Q²-floor scan"; `python/tests/test_o5_reach.py` | table values pinned at 1e-12 (they are the single code home of those numbers); the *derived* reach arithmetic at 1e-3…2e-3 | an **unpolarized rate and slope baseline** for the coherent channel, and the O5 feasibility arithmetic built on it (verdict: J/ψ tensor a₂ is **marginal**, 2.62 σ at 10 fb⁻¹/u, 3 σ at 13.1 fb⁻¹/u) | anything polarized. The caveats ride with every row and are stated in the header: **one spherically symmetric Gaussian density → no α+d clustering, no polarization axis, no diffractive minimum at any A; unpolarized throughout; no saturation.** A rate baseline, never an imaging one |
| D-5 | **LHAPDF 6.5.5 + the Python `parton` package** (CT18NLO, NNPDFpol11_100, EPPS21nlo_CT18Anlo_Li6) | `LhapdfSF::f2p/f2n`, `LhapdfG1::g1p/g1n`, `Epps21Ratio::ratio` against `PartonF2`/`PartonG1` scalar workers | `tests/test_lhapdf.cpp` (7 cases) | **1e-3** (grid vs analytic bicubic); measured ≤ 5e-4, mostly ~1e-5. The one row at CT18NLO's Q_min boundary disagrees at 2.7e-3 and is **reported rather than masked** | that the C++ PDF wrappers evaluate the same sets the Python did, including the finite-flavour convention difference (NNPDFpol11 carries c/b, `PartonG1` drops them) | the PDFs |
| D-6 | **BeAGLE / eSTARlight / PYTHIA8 feasibility surveys** | what each code can and cannot do for polarized light ions | `docs/surveys/beagle_survey.md`, `pythia8_survey.md`, `needs_survey.md` | prose | that the "no other generator does this" claim was **checked**, not assumed | — |

---

## 5. Table E — THEORY: published calculations, tables, closed forms and digitized figures

| # | calculation | what is compared | TEST_CASE | tolerance / verdict | validates | does **not** validate |
|---|---|---|---|---|---|---|
| E-1 | **Cosyn–Weiss II, `arXiv:2603.23700` p. 35, Eqs. (6.12)–(6.14) and TABLE II** — tagged deuteron tensor asymmetry A_T∥ | (a) exactness of the P₂(cos θ_k) factorization; (b) the k envelope peaks at 1 where f₂/f₀ = √2, at k = 0.30 GeV for AV18; (c) A_T∥ = −2 A_zz^wf reaches CW's **+1** and **−2** and stays in [−2, 1] everywhere | `tests/test_tagged.cpp` "tagged: the Cosyn-Weiss deuteron tensor gate (CW TABLE II)" | (a) spread across angle bins < **1e-5**, mean ratio 0.99940 at 1e-4; (b) peak at k = 0.3098 GeV, within **0.02** of CW's 0.30; (c) extremes 0.9997 and −2 at **3e-3** | **the single strongest external physics check in the tree.** It validates the tagged spin-1 wave-function tensor asymmetry — the D-wave/spin correlation that drives A_zz^tag — against an independent analytic calculation, at A = 2 | A = 6 or 7. It runs on the **deuteron control channel** (`deuteron_channel()`, which is guaranteed to carry no VMC table so the gate is always on the analytic pair). It also does **not** compare against CW's Fig. 13, which is in light-front (α_p, p_pT) variables the generator does not carry |
| E-2 | **CDKS Fig. 4** — Cosyn, Dong, Kumano, Sargsian, PRD **95** (2017) 074036, x·b₁(deuteron) at Q² = 2.5 GeV², theory-1 sum. Digitized from the PDF's own path operators (`tools/digitize_figure.py`), embedded as `tables::kB1CdksQ2p5` | **G3a (hard)**: exactly two sign changes in (0, 1], first within 0.08 of 0.06564, second within 0.10 of 0.45718, peak within 0.10 of 0.7657. **G3b**: max\|x·b₁\| over [0.10, 0.80] vs the digitized peak, window [0.5, 2]. **G3c (reported)**: ∫b₁ dx vs the digitized 4.592e-4 | `tests/test_b1_nuclear.cpp` "T1 gate layer 3: CDKS Fig. 4"; **T1v** "checklist item 4, THE VERDICT ROW: MSTW2008 LO at CDKS Eq. (21)"; checklist item 5 / **5v** (CD-Bonn) | **The verdict is configuration-dependent and the tree says so.** With the shipped default `ToyF2`: ratio **0.440**, i.e. **FAILS** G3b's [0.5, 2] window (the test asserts `ratio < 0.5`). With **CDKS's own MSTW2008 LO** at CDKS Eq. (21)'s finite-\|q\| δ-function: **0.843 — PASSES**. With MSTW + **CD-Bonn**: **1.000**, and every G3a landmark within 0.005 of the digitized column | that the α–d/deuteron b₁ convolution reproduces an independent published b₁ᵈ calculation in **shape** (two sign changes, positions, peak) and, in the like-for-like configuration, in **magnitude** to 16 % | that b₁(⁶Li) is right. This is an **A = 2** gate; the ⁶Li number is the same convolution machinery applied to the α–d overlap, with nothing external to check it against. And the pass is a statement about **one configuration**: `T1r` is a dedicated TEST_CASE that prints, in every build, whether the verdict row was actually measurable (it is in this build — MSTW2008 LO is on disk under `deps/install/share/Pythia8/pdfdata`) |
| E-3 | **Miller, PRC 89 (2014) 045203, Fig. 5** — deuteron b₁ = b₁(π) + b₁(6q); digitized as `tables::kB1Miller` | the digitized column reads back; the taper above x = 0.9; the constant extrapolation below x = 0.01 | `tests/test_sf.cpp` "the digitized Miller b1 curve reproduces the published figure" | rtol **1e-12** on the read-back | the digitization and the two extrapolation conventions | Miller's calculation, and — see C-7 — not agreement with HERMES. The **per-deuteron vs per-nucleon normalisation** of this curve is flagged as **LIKELY, NOT CERTAIN** in `constants.hpp` and recorded as an open **author decision** (`OPEN_ITEMS_SOLUTIONS.md` §10): Miller's own paper is self-inconsistent by exactly this factor of 2 |
| E-4 | **Close–Kumano sum rule**, PRD **42** (1990) 2377: ∫b₁ dx = 0 | the digitized CDKS column integrates to 4.592e-4; Miller's to 5.915e-3 | `tests/test_sf.cpp` "the digitized CDKS convolution is the other b1 camp"; T1 gate G3c | 1e-14 on the integrals; the comparison `|I_CDKS| < 1e-3` and `|I_Miller| > 5|I_CDKS|` is **reported, not enforced** | that the two b₁ camps differ by an order of magnitude on a sum rule one of them (Miller, his §V) openly violates | the sum rule itself — the digitized table stops at x = 1.59 while the deuteron's x runs to 2, so `close_kumano_integral` **reports and does not enforce** |
| E-5 | **Cosyn *et al.*, EPJ A 61 (2025) 83 (`arXiv:2410.12764`) Table 1, rows 2 and 3** — the exact finite-γ tensor decomposition | b₁/(F₁ A_T) against the paper's **derived** closed forms, for a target polarized along **q** and along the **beam** | `tests/test_tensor_gamma.cpp` "Table 1 row 2, a target polarized along q", "…row 3, a target polarized along the beam" | rtol **1e-10** | the **transcription of Eqs. (9)/(10)/(14)/(16)/(17)/(24)** against a result the module has never seen. The file states that these two rows — and not the re-typing tests beside them — are what caught the two transcription errors in its history (the bracketing of Eq. 17a; the sign of Eq. 22b / the frame orientation) | the physics of b₃, b₄ — those slots are filled with **scenario shapes** (0.05·f₁, −0.02·f₁) because b₃ and b₄ are unmeasured |
| E-6 | **Cosyn Eq. (27) / HERMES Eq. (6) / HJM / POLRAD Eqs. 9–10**: A_zz = −(2/3)b₁/F₁ | `azz(...)·(1 + εR)` against the literature relation, written **without reference to `TENSOR_LL_SIGN`** so flipping the constant fails | `tests/test_xsec.cpp` "Cosyn Eq. 27: A_zz(theta_S=0)(1 + eps R) = -(2/3) b1/F1", "the program sign IS the literature sign, deliberately" | rtol **1e-12** at every y | the **sign and normalisation convention** of the whole tensor sector, agreed by four independent sources. This was a real change: the program's private +1 was flipped to the literature −1 on 2026-08-29 | — |
| E-7 | **Hoodbhoy–Jaffe–Manohar, NPB 312 (1989) 571** — the massless J = 1 tensor geometry | the (1, −2, 1) c_m pattern; `tensor_moments(m) = ((3m²−J(J+1))/3, 3×that)` for **both** J = 1 and J = 3/2 | `tests/test_xsec.cpp` "the J = 1 geometry is the HJM transcription digit for digit", "the rank-2 geometry is one formula for both spins" | 1e-13 / rtol 1e-12 | that one Cosyn Eq. (9) formula covers both spins and reduces to HJM at J = 1 | the J = 3/2 rank-2 **normalisation**, which is an open item (plans/04 #14). `tests/test_xsec.cpp` "the spin-3/2 rate and cos-2phi channels are mutually consistent" is labelled **CHARACTERIZATION, not a physics assertion** |
| E-8 | **E143 (Abe *et al.*, PRD 58 (1998) 112003)** lab-frame definitions of ε, D, η | the finite-γ factors against ε = 1/[1+2(1+ν²/Q²)tan²(θ/2)], D = (1−E′ε/E)/(1+εR), η = ε√Q²/(E−E′ε) on 64 (E, E′, θ) points | `tests/test_xsec.cpp` "the finite-gamma factors match E143's lab-frame definitions" | rtol **1e-12** | that the target-mass (finite-γ) vector sector is the standard one written in different variables | — |
| E-9 | **R1998 world fit** (Abe *et al.*, PLB **452** (1999) 194, `hep-ex/9808028`) | `r1998` = mean of its three published forms; `r1998_spread`, `r1998_fit_error`; support clipping | `tests/test_sf.cpp` "R1998 is the average of its three published forms" | rtol **1e-12** on all three forms and the average | the transcription of a published fit, and its stated support | — |
| E-10 | **Wandzura–Wilczek relation** | `g2_ww` on a power law against the analytic −xᵃ + (1−xᵃ)/a | `tests/test_sf.cpp` "g2_ww on a power law matches its analytic Wandzura-Wilczek form" | **2e-3** (96-point trapezoid), plus a 1e-12 pin against the Python's own quadrature | that the quadrature reproduces the closed form where one exists | — |
| E-11 | **Cloët–Bentz–Thomas, PLB 642 (2006) 210, Fig. 6** (⁷Li polarized EMC, Q² = 5) and **Tronchin–Matevosyan–Thomas, PLB 783 (2018) 247, Fig. 4** (nuclear matter, Q² = 10) — digitized as `kCbtPolemc7liQ5`, `kTmtPolemcNmQ10` | curve read-back; the ratio-of-effects (CBT 2.25/1.69/1.41/1.14 vs TMT 1.01/0.98/1.00/1.08 at x = 0.40/0.45/0.50/0.60); the valence transfer onto a common EPPS21 baseline | `tests/test_sf.cpp` "the digitized polarized-EMC curves and the valence transfer", "…names its unpolarized baseline" | rtol **1e-12** on the read-back; **0.05** on the ratio-of-effects | that the two polarized-EMC camps are carried as digitized curves rather than constants, and that the "≈ 2×" vs "≈ 1×" disagreement is reproduced pointwise | either calculation. Also: the digitization's own self-validation (CDKS Fig. 5 vs Fig. 4, two pages, two calibrations, agreeing to 1e-7) lives in the **Python** repo's `data/SOURCES.md`, not in a LiPolGen test |
| E-12 | **Mäntysaari *et al.*, PLB 858 (2024) 139053, Fig. 4** — polarized-deuteron coherent a₂(\|t\|) at 4 \|t\| points | `a2_from_quadrupole(2Q_d, 2, |t|, m)` against the four digitized rows, **zero free parameters** | `tests/test_cluster_config.cpp` T9 "the coherent.hpp bridge: a2_from_quadrupole IS the published a2" | **10 %** on m = ±1, **25 %** on m = 0; the exact m-relation a₀ = −2a₁ at 1e-14 | that the closed-form quadrupole → a₂ map reproduces an independent dipole-model calculation **for the deuteron** | ⁶Li. The map is then applied to ⁶Li with **no published calculation to check it against** (`coherent.hpp`: "No published calculation exists for any polarized A > 2 nucleus"). T10b shows the shipped default `eps_b0 = −0.08` corresponds to a ⁶Li quadrupole **11.4× the measured one** — i.e. it is a *deuteron* number and the test says so in its title |
| E-13 | **ANL VMC AV18+UIX / AV18+UX cluster overlaps** (R. B. Wiringa; Pudliner *et al.* PRC 56 (1997) 1720; Forest *et al.* PRC 54 (1996) 646) — `data/vmc/`, raw text tables fetched via the Wayback Machine (provenance in `data/vmc/README.md`) | (a) the r-space block Fourier–Bessel-transformed reproduces the **independently estimated** k-space block; (b) the S/D nodes at 0.678 and 2.250 fm⁻¹; (c) the S–D relative **sign**; (d) P_D(⁶Li) = 1.94 % against the file's own printed integrals; (e) the k-moments ⟨k⟩ = 0.1225 (⁶Li), 0.1864 (⁷Li) GeV, tail fractions; (f) the file's own MC errors carried as a band | `tests/test_cluster.cpp` (5 cases); `tests/test_tagged.cpp` "the VMC channels…", "VMC P_D(6Li) = 1.94%…", "the VMC moments reproduce the reconciled table", "T24 the VMC tables carry their Monte Carlo band", "T25 the deuteron control channel, Hulthen against AV18" | (a) **5e-3**; (d) 2e-4 abs, against the file's own ~1 % MC error on the D block; (e) 5e-4 against the table, **2 %** against an independent 4000-point Python integration on a different grid/cap | that an **ab-initio** wave function is read, transformed and normalised correctly; that the transform convention is right (validated first on the clean, non-MC `fdeut.av18`); that the Hulthen scenario P_D is **4.48× too large**; and that the ANL statistics are **not** the dominant systematic (1.1 %/σ on P_D vs 6.6 % for the wave-function choice, T26) | the ⁶Li wave function as truth: the independent 2004 AV18+UIX overlap gives P_D = 2.01 % vs the 2014/2024 momentum files' 1.94 % — a **Hamiltonian difference**, not an error. And T6 gates a **documented 4 % discrepancy** between the assembled configuration sampler's r_rms and `li6.density`'s |
| E-14 | **Ciofi degli Atti–Simula, PRC 53 (1996) 1689**, Tables A.1/A.3 and Eq. (28) | S₀/S₁ for A = 2, 3, 4; the n₀ high-k tail P(k > 0.1/0.2/0.3/0.45); ⟨k⟩ | `tests/test_triton_sf.cpp` (4 cases) | 1e-3 … 2e-3 | the transcription of a published spectral-function parameterisation, and that its own normalisation closes untuned | that the ⁷Li triton remnant is realistic — a Faddeev/AV18 three-body correlated (p, n, n) distribution remains the full replacement |
| E-15 | **Ciofi degli Atti–Kaptari, PRC 83 (2011) 044602, Eqs. (7)–(8)**; **Cosyn–Sargsian, `arXiv:1704.06117` Eq. (11)**; **Strikman–Weiss, PRC 97 (2018) 035209** | the Glauber cluster-spectator survival product and the f_XN profile (σ_tot(i+ε)e^{Bt/2}, ε = −0.5, B ≈ 6 GeV⁻²) | `tests/test_fsi.cpp` (13 cases): the FSI/IA → 1 as σ_XN → 0 limit, exact θ_k → π−θ_k symmetry, `weight/survival` averaging to 1, the pointwise vs integrated bracket | analytic limits at machine precision; the table rows at 5e-3 (B-2) | that the eikonal weight is built and applied correctly, and never moves a four-vector | the FSI physics: **no measurement and no independent code**. The two Glauber variants are pinned as *differing on essentially every event*, i.e. the model spread is exposed rather than hidden |
| E-16 | **Chang *et al.*, PRD 113 (2026) 032018 (`arXiv:2511.05638`)** — IR-8 far-forward intact-recoil efficiencies (d 47.12 %, ³He 32.23, ⁴He 29.42, **⁷Li 17.75**, ⁹Be 12.37, ¹²C 6.36, ¹⁶O 1.59) and the ³He beam-energy scan (0.3223 / 0.5438 / 0.9977 at 183 / 100 / 41 GeV) | table values, monotonicity in A and in E, the two local power laws (−0.866, −0.681), the fixed-rigidity structure | `tests/test_coherent.cpp` "the far-forward efficiency is not flat in beam energy", "the species list is fixed RIGIDITY, NOT fixed E/u" | 1e-12 on the table; 1e-3 on the derived slopes | that another group's published detector-simulation result is carried verbatim and used with the right lever arm (1.12–1.16 for the 117.9 → 99.5 GeV/u step, **not** the paper's own 1.687) | ⁶Li — **there is no ⁶Li row in that paper**, so ⁷Li is used as a stand-in. And the efficiency is a **recoil-nucleus number only**: it carries no central-detector acceptance and no decay-lepton reconstruction efficiency, so the O5 chain is missing a factor ≤ 1, direction **down**. The header says all of this; the test block titled "THE RETRACTED SENTENCE" is the arithmetic retraction of an earlier over-claim |
| E-17 | **Wiringa *et al.*, PRC 89 (2014) 024305 Table I** (VMC AV18+UX polarizations) | ⁷Li whole-nucleus P_p = +0.866, P_n = −0.037; ⁶Li cluster product 0.81123 = (1−1.5P_D^{αd})(1−1.5P_D^{d}) | `tests/test_beams.cpp` "the effective polarizations are per nucleon"; `tests/test_sf.cpp` "ToyG1 and the nuclear g1A" | 1e-14 / 1e-13 | that the effective polarizations are stored per nucleon and reassemble to the published whole-nucleus values, and that ⁶Li's comes from the **same two D-state probabilities the tagged sector uses** rather than a hard-coded literal | the neutron half: `DEVELOPMENT_PLAN.md` §4 row 4 calls P_n ≈ −0.037 an **open** gate, and `tests/test_tagged.cpp` "7Li triton polarization and the forward-limit gate" is where that open status is printed |
| E-18 | **George–Knutson band ↔ quadrupole dial ↔ a₂ map** (composite) | the two-factor Q budget (3.317 × 2.269 = 7.524), the η ↔ dial exact linearity, the GK ±1σ mapping | `tests/test_cluster_config.cpp` T22b, T23, "the design's sec. 8 table, all three alpha-d sources" | 1e-9 … 1e-13 | that the ⁶Li quadrupole is an **explicitly dialled, banded** input and that every downstream number (a₂, ε_b0, P_D^{αd}) moves with it coherently | it is a *bookkeeping* validation of a model band, not a physics check |

> **E-1 CORRECTED AND STRENGTHENED, 2026-09-06.** **Row E-1 of the table above** records the gate as it
> stood on the survey date. Its criterion (c) — *"A_T∥ = −2 A_zz^wf reaches
> CW's +1 and −2"* — was **wrong twice over**, and `07_cw_sign_investigation.md`
> settled both: the mapping is **`A_T∥ = +1 · A_zz^wf` exactly**, and the −2 the
> gate applied was masking an **inverted S–D interference sign** in
> `TaggedModel::build_amp2` (the partial-wave sum consumed ψ_L where it needs
> φ_L = i^L ψ_L). The three quoted pass values inherit that: the mean ratio
> **0.99940**, the envelope peak at **k = 0.3098 GeV**, and the extremes
> **0.9997 / −2** are all pre-fix, and the peak was matching CW's Eq. (6.14)
> *minimum* at f₂/f₀ = 1/√2 while calling it Eq. (6.13)'s *maximum* at √2 — on
> the **Hulthén** pair, whose f₂/f₀ never reaches √2 anywhere on the grid.
> **The gate was rewritten (`07` §8) and the code fixed** with one
> `(-1)^floor(L/2)` in `build_amp2`. It now runs on the **AV18** control CW
> quote TABLE II for, and criterion (a) is Eq. (6.12) as an *identity*:
> measured `max|A_zz^wf − CW| = 8.881784e−16` over all 280 × 96 = 26 880 cells
> on both deuteron controls, against **2.740499e+00** before. CW's own k
> landmarks reproduce at **0.298121** and **1.034872 GeV** (CW: "0.30",
> "≈ 1"), TABLE II's three rows at **−1.937124 / +0.999313 / +0.967340**
> against −2 / +1 / +1, and the curve stays inside CW's [−2, 1]. **This row's
> last column — "the single strongest external physics check in the tree" —
> is now true in a way it was not**: nothing in the rewritten gate is a
> self-consistency restatement.


---

## 6. Table F — CHAIN: software-interface checks

| # | check | TEST_CASE / record | result | validates | does **not** validate |
|---|---|---|---|---|---|
| F-1 | **HepMC3 → `npsim` (DD4hep, `epic_craterlake_10x100.xml`)** | `docs/T2_CHAIN.md` §3 item 1; `docs/OPEN_ITEMS_SOLUTIONS.md` §2; `DEVELOPMENT_PLAN.md` validation matrix | **PASSED 2026-09-02, 10/10 events.** Accepts the file as written — 10-digit ion codes for beam ⁶Li at status 4 and the α spectator at status 1 — with EDM4hep `MCParticles` carrying the full role chain | that the generator's HepMC3 Asciiv3 output is a legal ePIC-chain input | anything physics-side. **Not re-run by any `pytest`/`lipolgen_tests` invocation** — it needs the external `eic_xl`/`jug_xl` containers, which are outside this repository's dependency tree. The doc paragraph *is* the record of the one time it ran |
| F-2 | **HepMC3 → `abconv` → `npsim`** | same | **PASSED 2026-09-02, 10/10.** `abconv -p 1`'s automatic energy detection **cannot decode a ⁶Li ion** (per-nucleon energy comes out 0); the manual preset `abconv -p ip6_hiacc_100x10` works and its output also runs through `npsim` | that the Geant4-primary-generator bridge accepts the file with an explicit preset | same caveat as F-1; and the automatic-detection failure is a **known, unfixed** limitation of the downstream tool |
| F-3 | whole-record 4-momentum and charge conservation, every channel, every tier | `tests/test_t2.cpp` (9 cases), `tests/test_pipeline.cpp` "every channel conserves four-momentum and charge", `tests/test_breakup.cpp` T1 pipeline cases | worst relative 4p residual **1.4e-13** (⁶Li α), 5.9e-14 (⁷Li α), 2.2e-13 (d), 4.8e-14 (coherent); charge exactly 0 | that T1 (struck-cluster resolution) and T2 (PYTHIA/Pomeron tiers) close the record exactly | — |
| F-4 | HepMC3 round trip (write → read → compare) | `tests/test_hepmc.cpp` (6 cases), `tests/test_t2.cpp` round-trip cases, `tests/test_pipeline.cpp` "100 tagged events survive a HepMC3 round trip" | exact; HepMC2 ascii **throws at construction** by design (`npsim` would hit immediate EOF) | the writer/reader convention of `docs/HEPMC3_CONVENTION.md` | the **spin labels**: `spin_J`/`spin_M`/`P_z`/`P_zz`/`spin_axis_*` are a **proposed** convention with no ePIC-MC-group standard, and `abconv`/`npsim` do not carry event attributes through to EDM4hep |
| F-5 | determinism and thread-invariance | `tests/test_rng.cpp` (5), `tests/test_pipeline.cpp` "threaded generation is bit-identical to single-threaded", `tests/test_rc_pipeline.cpp` T18b, `tests/test_t2.cpp` "determinism — same seed gives the identical HepMC3" | bit-identical at 1 and 8 threads, including the RC weights and the byte-level HepMC3 output | reproducibility | — |
| F-6 | `--rc tensor-band` is `--rc off` bit for bit, RNG included | `tests/test_rc_pipeline.cpp` T6 | exact | that the RC band is a **weight family and never a shift** — no four-vector moves | — |
| F-7 | estimator closure (Monte-Carlo self-consistency) | `tests/test_sampler.cpp` (6 closure cases) | pulls unbiased; spreads within **15 %** of 1/(P_e P_z √N), √(2/N)/P_zz; the relative-luminosity biases −(2/3)δ/P_zz and δ/(2P_eP_z) reproduced **and removed** | that the analysis chain recovers an injected asymmetry from generated events | INTERNAL closure — the generator checking its own estimator |
| F-8 | PYTHIA `PDF:PomSet` systematic scan | `python/tests/test_pom_set_band.py`; full scan in `docs/open_items/run_2026-09-03/phase_D_numbers.md` §D4 | 15 sets × 20 000 events; T0 columns bit-identical across the fourteen that run (set 11 refused); the band is over the **twelve DPDF fits** (3–10, 12–15) about set 6 | that the coherent T2 tier's largest model systematic is measured, not asserted | — |

---

## 7. Table G — documentation / provenance integrity gates

These are not physics checks, but they are what keeps every number in Tables C–F
traceable, so they belong in the baseline.

| # | gate | script / test | what it enforces |
|---|---|---|---|
| G-1 | citation gate | `validation/check_physics_channels_links.py`; `python/tests/test_doc_link_gate.py` | every `` `path:line` `name` `` citation in `docs/PHYSICS_CHANNELS.md` (1075 of them), `docs/theory/SPIN32_FINITE_GAMMA.md` and `docs/PYTHIA_BRIDGE.md` still lands on a **declaration line of that symbol, as code** — not on a comment that merely names it |
| G-2 | line-range gate | `validation/physics_channels_ranges.json` (sha256 + first/last line of every cited range) | a quoted block cannot drift without the gate noticing. **Do not run with `--fix` or `--record-ranges`** |
| G-3 | SPDX gate | `validation/check_spdx_headers.py`; `python/tests/test_spdx_headers.py` | every first-party source file declares GPL-3.0-or-later |
| G-4 | knob provenance | `python/tests/test_knob_provenance.py`, `test_meta_provenance.py` | **a knob that did not run may not be recorded in `meta` or printed in the banner as if it had.** The matrix is every knob × every channel × every plan, checked against the **output hash** of two runs. It caught five distinct failures of that rule in one run |
| G-5 | release metadata | `python/tests/test_release_metadata.py` | `CITATION.cff` / `pyproject.toml` / version consistency |

---

## 8. Which gates are build-dependent (and their status in this checkout)

| gate | condition | status here (verified 2026-09-06) |
|---|---|---|
| E-2 verdict row (T1v, 5v) | `LIPOLGEN_HAVE_PYTHIA8` **and** `mstw2008lo.00.dat` on disk | **MEASURED.** `T1r` prints: "VERDICT ROW MEASURED: MSTW2008 LO is on disk under `…/deps/install/share/Pythia8/pdfdata`, so T1v ran and G3b's 0.843243 was checked in THIS build." |
| E-2 CT18NLO row, D-5 | `LIPOLGEN_HAVE_LHAPDF` | compiled in (`libLiPolGenLHAPDF.so` present) |
| E-13, B-2 VMC cases | `data/vmc/*` present | present and committed (1.1 MB) |

| Table A | `LIPOLGEN_REFERENCE_DIR` populated | populated (8 JSON files, `_manifest.json` dated 2026-09-04) |
| F-1, F-2 | external `eic_xl`/`jug_xl` containers | **not runnable from this tree**; the 2026-09-02 pass is recorded in prose |

Spot-checked live for this survey (`./build/lipolgen_tests -tc=…`, 6 cases,
54 522 assertions, all passing): E-1 (CW Table II), E-2 (CDKS Fig. 4 layer 3),
C-6 (Machleidt Table XX), C-1 (T22 η), D-4 (eSTARlight baseline), and T1r.

---

## 9. The honest gaps — what nothing in this tree checks

1. **Polarized e + ⁶Li/⁷Li DIS, at any level.** No data, no second generator, no
   published calculation of the tensor observables for A > 2. Every ⁶Li/⁷Li
   tensor number is a model extrapolation from A = 2, and the only external
   contact points are E-1 (deuteron, tagged) and E-2 (deuteron, inclusive b₁).
2. **b₁(⁶Li).** The α–d convolution is validated **only** at A = 2 (E-2), and
   even there the pass is configuration-specific (MSTW2008 LO + CD-Bonn, 1.000;
   shipped `ToyF2` default, 0.440 — outside the window).
3. **b₂, b₃, b₄, Δ.** Unmeasured everywhere. `b₂ = 2x·b₁` is a default; b₃/b₄
   carry scenario shapes; Δ is a scenario ansatz (`toy_delta_gluon`, `C_BAG`).
4. **The J = 3/2 rank-2 normalisation** (⁷Li). Open item plans/04 #14; the
   corresponding test is labelled CHARACTERIZATION. `tests/test_pipeline.cpp`
   D1 pins that the ⁷Li rank-2 sector is **identically zero and says so loudly**.
5. **Tensor-sector radiative corrections.** One dedicated calculation exists
   (Gakh–Shekhovtsova, `hep-ph/0403262`), with **zero INSPIRE citations**, and it
   supplies only a band half-width read off *one panel's x window at a single
   Q² = 0.1*. Mo–Tsai (RMP 41 (1969) 205) is named in `tests/test_rc.cpp` as the
   gate that was **deliberately not written** and was not obtained for this run.
6. **Spectator FSI.** No measurement, no independent code (E-15). Cosyn–Weiss II
   §VI C lists the spin dependence of FSI as an open question.
7. **Coherent ⁶Li tensor amplitude.** No published polarized A > 2 diffraction
   calculation. The a₂ map is anchored on the **deuteron** (E-12) and the shipped
   `eps_b0` default is a deuteron number 11.4× the measured ⁶Li quadrupole.
8. **⁶Li far-forward efficiency.** Chang *et al.* publish no ⁶Li sample and no
   absolute cross section (E-16); ⁷Li is a stand-in, and the number carries no
   decay-lepton reconstruction efficiency.
9. **Elastic ⁶Li form-factor data.** None in this repository; the C0 zero is a
   model band ([2.9, 3.3] fm⁻¹) around an admitted starting guess.
10. **The linearly-polarized-photon cos 2φ background.** ZEUS's measured
    diffractive azimuthal asymmetries (A_LT, A_TT, `arXiv:0812.2003`), which
    bound the `u₁ = 0.05` / `u₂ = 0.02` background in the **Python** `fastsim`,
    have **no counterpart in LiPolGen**: grep finds no `A_LT`, `u1`/`u2`
    anywhere in `src/`, `include/` or `tests/`. The only H1/ZEUS mention in the
    C++ is `coherent.hpp:394`, and it is a statement that the calculation the
    coherent fraction `f0` would need — a coherent-A analogue of the H1/ZEUS
    diffractive PDFs — **exists in neither eSTARlight nor Sartre**. The P_zz
    flip is the separation argument (D-4), and it is arithmetic, not a
    measured background subtraction.
11. **Spin labels in HepMC3** have no community convention (F-4), so nothing
    downstream of `npsim` can be checked against them.

## 10. One-line summary

The tree's external validation is **real but narrow**: it is strongest at
A = 2 (Cosyn–Weiss Table II, CDKS Fig. 4, CD-Bonn Table XX, Mäntysaari a₂,
POLRAD), solid on conventions and transcriptions (Cosyn Eq. 27, Table 1 rows 2/3,
E143, R1998, HJM), credible on unpolarized nuclear inputs (AME2020, VMC overlaps,
LHAPDF/MSTW, eSTARlight, Chang *et al.*), and complete at the software-chain level
(npsim/abconv, HepMC3, conservation, determinism). Everything specific to
**tensor-polarized lithium** rests on model extrapolation from those A ≤ 2 anchors
plus internal `polligen` regression locks, and the tree is unusually explicit
about saying so — several of its "gates" (C-1, C-2, E-13/T6, E-16) exist to pin a
**documented discrepancy** rather than an agreement.
