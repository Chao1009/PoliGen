# Benchmarking LiPolGen — the plan

*Written 2026-09-06 from the survey in `00_in_tree_checks.md` … `06_critic.md`
(201 resources, 16 high-priority ones adversarially verified, 4 refuted). The
survey files are the evidence; this file is the plan. Re-read `06_critic.md`
§5 before wiring any data row — the per-nucleon / per-nucleus / per-deuteron
trap is where this project has failed repeatedly.*

## 1. Why this is structured by decomposition

There is no polarized e + ⁶Li or ⁷Li DIS data, and no other Monte Carlo does a
tensor-polarized target of any nucleus (`01_generators.md` §1 — verified as a
negative by INSPIRE full-text search, `06_critic.md` §4.2). **Nothing in this
tree validates, end to end, the observable the generator exists to produce.**
That sentence goes into every release note until data exists.

What CAN be validated is every layer the observable is built from:

- the nucleon inputs (F₂, R, g₁) against world data and the PDF fits;
- the deuteron machinery (b₁ᵈ, A_zzᵈ, tagged spin-1 structure, the wave
  function) against HERMES, JLab, and the CDKS / Miller / Cosyn–Weiss
  calculations — this is where the tensor physics is actually anchored;
- the nuclear inputs (radii, moments, form factors, EMC ratio, cluster
  momentum distributions) against their measurements;
- overlapping channels against well-tested generators (PYTHIA 8 ep, DJANGOH
  RC, STEG/BeAGLE tagged e+d, eSTARlight exclusive VM);
- MC-internal closures; and the ePIC chain.

And what cannot, listed in §6 with the closest proxy and its distance.

## 2. The tiers

| tier | what is compared | reference kind | owns the reference's own validation |
|---|---|---|---|
| **T1 nucleon** | F₂ᵖ, F₂ᵈ, F₂ⁿ/F₂ᵖ, R, g₁ᵖ/ᵈ at EIC (x, Q²) | HEPData: NMC, BCDMS, HERMES, SLAC E143/E155, COMPASS, H1 F_L, R1990/R1998 | the experiments; the PDF fits |
| **T2 deuteron** | b₁ᵈ/A_zzᵈ (HERMES Table II), tagged spin-1 structure (Cosyn–Weiss Table II), b₁ᵈ shape and magnitude (CDKS Fig. 4, Miller Fig. 5), the wave function (CD-Bonn/AV18 A(Q²), t₂₀, (e,e′p) momentum distributions) | data + published calculations | the A = 2 gate in-tree; HERMES is the ONLY tensor-DIS datum in the world |
| **T3 nuclear inputs** | ⁶Li/⁷Li ⟨r²⟩, μ, Q; elastic F_C, F_M (Suelzle–Yearian–Crannell, Li–Sick–Whitney–Yearian, the UVa Fourier–Bessel density, Rand–Frosch–Yearian magnetization); NMC F₂(⁶Li)/F₂(D) (HEPData ins394050); ⁶Li(e,e′d)α and (e,e′α)d cluster-knockout spectra (Ent, Mitchell); ⁶Li QE (e,e′) (QES archive, 133 points); ANL VMC/GFMC overlaps, momentum distributions, densities | data + ab initio | TUNL/Stone; the (e,e′) groups; ANL |
| **T4 generator-vs-generator** | unpolarized inclusive closure (PYTHIA 8 ep, the official ePIC EVGEN samples); RC (DJANGOH paired Rad/noRad tables; POLRAD ADGH compiled numbers; PEPSI+RADGEN); tagged e+d spectra (STEG; BeAGLE eH2 official samples); exclusive VM rates/slopes (eSTARlight) — **as a rate scale, NOT the coherent continuum channel** | each generator's published validation vs HERA/E665/NMC | MCnet-style: cite their validation pages |
| **T5 MC closures** | estimator errors, weighted/unweighted, conservation, spin-density trace/positivity, sum rules, determinism, the knob-provenance matrix | analytic | this repository |
| **T6 chain** | HepMC3 → abconv → npsim → **EICrecon (never run)**; Rivet with the 46 HERA DIS analyses on the PYTHIA-tier output; NuHepMC validator conventions; ion-spin attribute convention | ePIC software | ePIC MC/software groups |

## 3. Milestones, and what is re-run at each

The project's hard date is the PLB letter / INT window, **22 March – 2 April
2027** (`../../PolarizedLithiumSim/plans/07`). Phase 1 (this generator) is
in progress; Phase 2 (full simulation, `plans/03`) has not started.

| milestone | when | trigger | tiers | pass rule |
|---|---|---|---|---|
| **M1 kernel closure** | now — the 2026-09-06 run's close-out | this plan adopted | T1, T2, T5 in full; T3 for what is in-tree; the §4 top-ten wired | every gate green; every new data comparison recorded with tolerance and **no tuning**; §5's owed corrections applied |
| **M2 numbers frozen** | before the letter's numbers are computed (target Jan 2027) | plans/07 WP freeze | T1–T5 in full; T6 smoke | no benchmark moved beyond its stated tolerance since M1, or the move traces to a registry decision |
| **M3 full-simulation chain** | Phase 2's first end-to-end Li run | fullsim available for Li species | T6 in full incl. EICrecon; T4 tagged/coherent vs fullsim acceptance | chain identities hold; acceptance-folded rates agree with the fastsim to a stated tolerance |
| **M4 letter** | 22 Mar – 2 Apr 2027 | plans/07 | everything, as the release's benchmark report | the report ships with the release; §6's list is in it verbatim |
| **M5 recurring** | every tagged release; any change touching a kernel, wave function, PDF set, RC model, or the chain | CI release job + checklist | the tiers the change touches (§7 map) | as M2 |

**A benchmark that is not automated is not a benchmark.** Each row of the
resource matrix gets `validation/benchmarks/<tier>_<name>.py`, which (i) reads
its reference from a vendored table under `validation/benchmarks/data/` with a
provenance header stating the source, the fetch date, the license, and the
**publishing convention** (per nucleon / per nucleus / per deuteron / ratio;
Cartesian vs spherical for tensor observables), (ii) runs the generator or the
kernel at the reference's kinematics, (iii) writes one row of
`docs/benchmarking/REPORT.md` — value, reference, tolerance, pass/fail, and
the configuration that made it. Cheap rows run in pytest always; expensive
rows (a full generation, a chain leg, a generator-vs-generator sample) behind
a `benchmark` marker that the milestone checklist and the CI release job
enable. The report is regenerated, never hand-edited.

## 4. The first ten to wire (most validation per unit effort, `06_critic.md`)

| # | benchmark | tier | effort | what it buys |
|---|---|---|---|---|
| 1 | **EPIOS Table II source modes vs `spin1_populations`** | T5/T1 | hours | the FIRST external anchor on the spin bookkeeping (today only a polligen self-pin); also shows the `ladder` fill is the wrong model for an EIC ion beam |
| 2 | **Cosyn–Weiss II Table II — re-derive the A_T∥ ↔ A_zz^wf mapping** | T2 | hours to read, unknown to fix | **blocker.** The verification says the tree's −2 factor has no basis, the correct mapping is +1, and on the AV18 control the gate would FAIL by the sign of the S–D interference — the correlation that drives A_zz^tag. See §5 |
| 3 | **HERMES b₁ᵈ Table II as an assertion** (six rows, transcribed) | T2 | low | converts the only measurement of the target observable from a comment into a gate; forces registry row 2 against data |
| 4 | **NMC F₂(⁶Li)/F₂(D), HEPData ins394050** (24 points, CC0) | T3 | ½ day | the only ⁶Li DIS measurement; tests the EPPS21 baseline (χ² = 4.63/4 measured in the survey). Divide by the isoscalar nucleon — the wrong denominator gives a plausible EMC curve 9× too large |
| 5 | **UVa ⁶Li Fourier–Bessel charge density** (7 coefficients) | T3 | hours | a data-derived C0 form factor; its zero at 2.695 fm⁻¹ and the diffraction minimum at 2.828 both sit OUTSIDE the tree's asserted [2.9, 3.3] window |
| 6 | **Rand–Frosch–Yearian ⁶Li/⁷Li magnetization densities** | T3 | ½–1 d | the q → 0 slope of F_M that has no data behind it today; delivers ⁷Li, for which `HoSpin1FF` throws; the only measured electromagnetic signature of α+d clustering |
| 7 | **DJANGOH published Rad/noRad four-bin table** | T4 | hours | turns the RC "must-not-contradict" bound into a Q²-trend comparison; zero dependencies |
| 8 | **George–Knutson η — retitle** | T3 | hours | the verification refuted it AS A BENCHMARK (it is a ⁶Li+⁴He phase-shift analysis, not d+α tensor analysing powers, and η is exactly linear in the quadrupole dial): relabel "consistency band on a dial" everywhere |
| 9 | **EST identity + COMPASS common-spin-temperature citation** | T5 | hours | `spin_temperature_pzz` IS the polarized-target closed form (8e-15); also the honest limitation: P_zz of ⁶Li has never been measured by anyone |
| 10 | **CLAS BONuS / EG1b database query** | T2 | hours | if machine-readable, the closest experimental analogue of the tagged topology at A = 2 |

Then, by tier, everything in `01`–`05` marked high that is `hepdata`,
`table-in-paper` or `runnable-code`; `digitizable-figure` rows only where no
table exists and the digitization error is stated on the row.

## 5. Corrections the survey found owed in the tree (apply at M1)

1. **Cosyn–Weiss gate (`tests/test_tagged.cpp:492-548`)** — the mapping
   A_T∥ = −2·A_zz^wf must be re-derived (CW Eq. 6.12 is the ordinary
   (n₊₁ + n₋₁ − 2n₀)/Σ ratio, i.e. +1·A_zz^wf). If the independent
   re-derivation confirms it, `build_amp2` (`src/core/tagged.cpp:388-428`)
   lacks the i^L phase CDKS state explicitly ("φ₂(p) < 0 due to the i^L
   factor"), the S–D interference sign in the tagged sector is inverted, and
   `b1_nuclear.hpp:335` (which writes +U W/√2, the CDKS sign) may carry the
   opposite sign from the tagged sector. **This is being re-derived
   independently before anything is changed.**
2. **George–Knutson** — technique misattributed at every site (see #8).
3. **eSTARlight** — the claim line "benchmarks `coherent.hpp`'s channel rate"
   is false: eSTARlight is exclusive ρ/φ/J/ψ with 85–97 % of its rate at
   Q² < 0.1 where `coherent.hpp` generates nothing; it is a rate SCALE for
   exclusive VM, not a check on the M_X ≥ 1.2 GeV continuum channel.
4. **POLRAD ADGH entry** — scope overstated in five places, one stale
   number; the three numbers themselves reproduce (re-downloaded, re-driven).
5. **`b1_default_li6.json`** is a self-pin dumped from LiPolGen's own C++ —
   label it so; it is not a reference of any kind.

## 6. What cannot be benchmarked until data exists (`06_critic.md` §2)

| # | observable | closest proxy | distance |
|---|---|---|---|
| U-1 | b₁(⁶Li), b₂(⁶Li) | b₁ᵈ (HERMES; CDKS/Miller) | A = 2 → 6 through `LI6_B1_PER_NUCLEON` = 2/6 — a modelling ASSUMPTION (inert α), not a convention |
| U-2 | tagged tensor asymmetry with an α/t spectator | Cosyn–Weiss A_T∥ with a NUCLEON spectator | different spectator, different wave function; and see §5.1 |
| U-3 | coherent ⁶Li tensor cos 2φ / a₂ | Mäntysaari *et al.* polarized-DEUTERON a₂(|t|) | the deuteron's quadrupole is 3.5× ⁶Li's and of the opposite sign; O5 is a closed-form estimate, not a Good–Walker amplitude |
| U-4 | ⁷Li rank-2 sector | Jaffe–Manohar (formalism only); Fu–Sun–Dong GPD forward limits | no number for ⁷Li anywhere; the sector is identically zero in-tree |
| U-5 | b₃, b₄, Δ (double-helicity-flip gluon) | Detmold–Shanahan lattice A₂ for the φ meson | a different hadron |
| U-6 | tensor-sector radiative corrections | Gakh–Shekhovtsova (zero citations, one Q² panel, deuteron) | no calculation for any A > 2 |
| U-7 | cluster-spectator Glauber FSI | Ciofi degli Atti–Kaptari ³He(e,e′d)X | the only cluster-spectator FSI calculation; a d spectator on A = 3 |
| U-8 | the α–d D-wave from data | George–Knutson η | refuted as a benchmark; η is a relabelling of the quadrupole dial |
| U-9 | A_zz(⁶Li) inclusive | the relation A_zz = −(2/3) b₁/F₁ (four sources agree) | a convention gate on the RELATION, not on either side |
| U-10 | ⁶Li elastic FF shape in the RC tail | F_q(0), F_m(0) on measured moments | the C0 SHAPE is a model band (now: UVa density is a data-derived shape — #5 above) |
| U-11 | the spin bookkeeping | nothing external today | #1 above fixes this |
| U-12 | polarized ⁶Li QE/elastic tails | 133 unpolarized QE points | the tensor QE tail is priced, not computed |

## 7. The change → tier map (what a diff must re-run)

- `sf.hpp`, `lhapdf_sf.cpp`, `mstw_sf.cpp`, any PDF-set version → T1, T2, T3-EMC
- `b1_nuclear.*`, `cluster.*`, `cluster_config.*`, `tagged.*`, `data/vmc/**` → T2, T3
- `breakup.*`, `spectator.*`, `fsi.*`, `triton_sf.*` → T3-cluster, T4-tagged
- `coherent.*`, the Pomeron tier, `pythia_bridge.*` → T4-inclusive, T4-VM, T6
- `rc.*` → T4-RC, T2 (HERMES b₁ was RC-corrected — the comparison states what RC it assumes)
- `spin.*`, `bookkeeping.*`, `sampler.*`, `rng.*` → T5 (and T1 #1)
- `hepmc_writer.*`, `cli.py`, `export.py` → T6

## 8. The three normalisation rules (`06_critic.md` §5.2), on every harness file

1. Never compare two b₁ values without naming both normalisations AND the
   A-scaling assumption (0.5 per-deuteron→per-nucleon; 2/6 the inert-α
   assumption; CDKS's x-to-2 convention).
2. Nuclear PDF grids are the average nucleon: divide by ½(F₂ᵖ + F₂ⁿ), always,
   and assert it in the test.
3. Cross sections on a nucleus are per nucleus; structure functions on a
   nucleus are per nucleon. Nothing in either convention is written on the
   file — the harness header writes it.
