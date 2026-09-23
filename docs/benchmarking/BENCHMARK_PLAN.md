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
| **M1 kernel closure** | now — the 2026-09-06 run's close-out. **EXECUTED 2026-09-23** (run `../open_items/run_2026-09-23/`, phases A, B1–B3): the §4 top ten wired or recorded (§4 status column), §5 items 1–5 applied (§5), every in-tree gate green, no row tuned. **NOT re-run at M1:** no T1 harness exists (the HEPData F₂/R/g₁ comparisons of `02_data_nucleon_deuteron.md` are unwired — T1 is covered only by the in-tree suites); T2's wave-function data (A(Q²), t₂₀, (e,e′p)) unwired; rows 6, 7, 10 BLOCKED on a document or a run; T4 beyond row 7 and T6 not run; no `REPORT.md` generator exists, so the rows live in the harness output and in the run record | this plan adopted | T1, T2, T5 in full; T3 for what is in-tree; the §4 top-ten wired | every gate green; every new data comparison recorded with tolerance and **no tuning**; §5's owed corrections applied |
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
enable. The report is regenerated, never hand-edited. (As built 2026-09-23 the harnesses print and return their row — `run()` returns a dict, a `REPORT |` line is printed — and no `REPORT.md` generator exists yet. The expensive-row switch is the environment variable `LIPOLGEN_BENCH=1`, not a pytest marker; no row needs it yet — the slowest, row 7, takes 17.3–19.6 s — so no test reads it today. Conventions and the harness list: `validation/benchmarks/README.md`.)

## 4. The first ten to wire (most validation per unit effort, `06_critic.md`)

*Status 2026-09-23: five wired (1, 3, 4, 5, 9 — three PASS, two FAIL recorded), one wired but BLOCKED because it cannot be a benchmark as planned (7), two BLOCKED on a document (6, 10), one resolved in code (2), one a relabel (8). Nothing was tuned. Every harness lives in `validation/benchmarks/` (`README.md` there); runtimes measured 2026-09-23.*

| # | benchmark | tier | effort | what it buys | status (2026-09-23) |
|---|---|---|---|---|---|
| 1 | **EPIOS Table II source modes vs `spin1_populations`** | T5/T1 | hours | the FIRST external anchor on the spin bookkeeping (today only a polligen self-pin). **WIRED 2026-09-23** (`validation/benchmarks/t5_epios_source_modes.py`, pass): all eight modes round-trip with residual 0 at tolerance 1e-12. Measured: the `ladder` (EST) fill reproduces the source's ideal P_zz at 3 of 8 modes (0, and the two \|P_z\| = 1 pure states where the tree throws) and misses the five partially polarized ones by 0.37 to 1.09 in P_zz (EST − ideal = +0.367, −0.915, +1.085, +0.697, +0.697 at modes 1, 2, 3, 4, 7); EST cannot produce P_zz < 0 at all. An atomic-beam source is not an EST system; which fill an EIC lithium beam has is not measured (no Li source table exists) | **wired — PASS.** Harness `t5_epios_source_modes.py` (0.3 s); pytest `test_bench_spin_deuteron.py::test_epios_report_row_passes` (+3 `test_epios_*`), doctest *"bench row 1: …"* (`tests/test_bench_spin.cpp`). Measured: max round-trip residual **0** over 8 modes, tolerance 1e-12 |
| 2 | **Cosyn–Weiss II Table II — re-derive the A_T∥ ↔ A_zz^wf mapping** | T2 | hours to read, unknown to fix | **blocker — RESOLVED AND FIXED 2026-09-06.** The verification says the tree's −2 factor has no basis, the correct mapping is +1, and on the AV18 control the gate would FAIL by the sign of the S–D interference — the correlation that drives A_zz^tag. See §5. All of it confirmed by `07_cw_sign_investigation.md` and fixed in one line of `src/core/tagged.cpp` `build_amp2`; Eq. (6.12) now holds as an identity (8.88e−16 over 26 880 cells, was 2.74) and ⁶Li A_zz^tag(0.20 GeV) moved +0.8450 → −1.2069 / +0.4518 → −0.5191 | **resolved and fixed 2026-09-06** (not a harness; §5.1). The CW TABLE II gate is `tests/test_tagged.cpp` *"the Cosyn-Weiss deuteron tensor gate"* |
| 3 | **HERMES b₁ᵈ Table II as an assertion** (six rows, transcribed) | T2 | low | converts the only measurement of the target observable from a comment into a gate; forces registry row 2 against data. **WIRED 2026-09-23** (`validation/benchmarks/t2_hermes_b1_table2.py`, fail — a measurement, no tuning): χ²/6 on b₁ᵈ = 22.26 (A = 2 convolution, MSTW; indistinguishable from b₁ = 0 at 21.80), 5.26 (shipped Miller ×0.5), 5.30 (Miller ×1). HERMES does not decide registry row 2 (`../open_items/run_2026-09-23/phase_B1_spin_deuteron.md` §3) | **wired — FAIL, recorded.** Harness `t2_hermes_b1_table2.py` (4.0–4.4 s); pytest `test_bench_spin_deuteron.py::test_hermes_report_row` (+5 `test_hermes_*`, `test_r1990_typo_guard`). Tolerance: every bin inside its stat ⊕ syst error — no configuration meets it (2/6 … 5/6 bins within). Measured χ²(b₁ᵈ)/6 = 22.264 / 22.669 / **5.262** / **5.297** / 21.782 / 21.799 (conv. MSTW / conv. ToyF2 / Miller ×0.5 shipped / Miller ×1 / CDKS digitized / b₁ = 0) |
| 4 | **NMC F₂(⁶Li)/F₂(D), HEPData ins394050** (24 points, CC0) | T3 | ½ day | the only ⁶Li DIS measurement; tests the EPPS21 baseline (χ² = 4.6269/4 on x ≥ 0.30, p = 0.33, and 26.94/15 on all 15 on-grid points, p = 0.029 — wired 2026-09-23, `validation/benchmarks/t3_nmc_li6_over_d.py`, pass). Divide by the isoscalar nucleon — the wrong denominator gives a plausible EMC curve 9× too large | **wired — PASS.** Harness `t3_nmc_li6_over_d.py` (1.3–1.4 s); pytest `test_bench_nuclear.py::test_nmc_chi2_reproduces_the_survey_and_passes` (+3 `test_nmc_*`). Tolerance p(χ², ndf) ≥ 0.01 on both sets, chosen after the survey's number was known and stated in the header. Measured **4.6269/4** (p = 0.328) and **26.941/15** (p = 0.029); 9 of 24 points below the EPPS21 grid, reported not compared |
| 5 | **UVa ⁶Li Fourier–Bessel charge density** (7 coefficients) | T3 | hours | a data-derived C0 form factor; its zero at 2.694 fm⁻¹ and the diffraction minimum at 2.828 both sit below the tree's asserted [2.9, 3.3] window — wired 2026-09-23, `validation/benchmarks/t3_li6_charge_ff_fb.py`: the shipped q₀ = 3.0998 is +0.2713 fm⁻¹ above the band [2.694, 2.828] (recorded FAIL, nothing moved; Q1 priced in `../open_items/run_2026-09-23/phase_B2_nuclear.md` §3.3) | **wired — FAIL, recorded.** Harness `t3_li6_charge_ff_fb.py` (0.3 s); pytest `test_bench_nuclear.py::test_fb_tree_side_and_verdict` (+2 `test_fb_*`). Tolerance: q₀ inside the band [2.6944, 2.8284] fm⁻¹ (must not contradict). Measured q₀ = **3.0998** fm⁻¹, **+0.2713 fm⁻¹ (+9.6 %)** above the band; T11's window not moved; the FB row's provenance is still UNVERIFIED |
| 6 | **Rand–Frosch–Yearian ⁶Li/⁷Li magnetization densities** | T3 | ½–1 d | the q → 0 slope of F_M that has no data behind it today; delivers ⁷Li, for which `HoSpin1FF` throws; the only measured electromagnetic signature of α+d clustering | **blocked.** Harness `t3_li_magnetization_rfy.py` exists (0.3 s); pytest `test_bench_nuclear.py::test_rfy_magnetization_benchmark` SKIPS naming both documents: Rand–Frosch–Yearian, Phys. Rev. 144 (1966) 859 and de Jager *et al.*, ADNDT 14 (1974) 479 — `../references/REFERENCES.md` §1 **#8, #9**. Tree side measured, uncompared: r_mag(⁶Li) = 2.9469 fm |
| 7 | **DJANGOH published Rad/noRad four-bin table** | T4 | hours | **wired 2026-09-23 as `t4_djangoh_rad_noRad` and BLOCKED**: the ePIC samples ran with the elastic radiative tail OFF (IEL2=IEL31=IEL32=IEL33=0), the only O(α) term `rc_tail` computes. The two share no term, so no tolerance exists (`../open_items/run_2026-09-23/phase_B3_chain_rc.md` §1). It unblocks with a DJANGOH run at IEL31..33 ≠ 0. | **wired — BLOCKED (not a benchmark as planned).** Harness `t4_djangoh_rad_noRad.py` (17.3–19.6 s, under the 30 s opt-in line); pytest `test_bench_chain_rc.py::test_djangoh_row_is_blocked_and_prints_both_sides` (+2). Measured: all eight ePIC DJANGOH logs have IEL2 = IEL31 = IEL32 = IEL33 = 0, so the published Rad/noRad ratios (1.0725 / 1.0894 / 1.1606 / 1.2280) are the inelastic correction with the elastic tail off, and LiPolGen's proton elastic tail (t-peak 1.0162 / 1.0132 / 1.0150 / 1.0207) shares no term with them — no tolerance can be defined. Unblocks with a DJANGOH run at IEL31..33 ≠ 0 |
| 8 | **George–Knutson η — retitle** | T3 | hours | the verification refuted it AS A BENCHMARK (it is a restricted phase-shift analysis — PRC 59 (1999) 598's own title; the scattering system it analyses is not verified in this tree, the paper is not held — not a direct d+α tensor-analysing-power measurement, and η is exactly linear in the quadrupole dial): relabel "consistency band on a dial" everywhere — **done 2026-09-23 for the code comments** (`../open_items/run_2026-09-23/phase_B2_nuclear.md` §5.2); the doc sites are G-1 … G-17 there. | **recorded, not a harness.** Code comments relabelled 2026-09-23 (`include/lipolgen/cluster_config.hpp`, three `///` blocks); doc sites G-1 … G-16 applied 2026-09-23; the T22 test name now says what it checks; G-17 (the Mäntysaari draft) left to the author, registry row 21 |
| 9 | **EST identity + COMPASS common-spin-temperature citation** | T5 | hours | `spin_temperature_pzz` IS the polarized-target closed form: max deviation 1.14e-14 over 19 999 points (8e-15 was five points; the bisection bound is 2.13e-14), **WIRED 2026-09-23** (`validation/benchmarks/t5_est_identity.py`, pass at 2.31e-14), with the COMPASS ⁶LiD populations, T_S and P(⁶Li), P(⁷Li) reproduced at the printed digit. The honest limitation: no measurement of the tensor polarization of ⁶Li was found by any search this project ran (see 06_critic.md §4.1 N-2) | **wired — PASS.** Harness `t5_est_identity.py` (0.35–0.37 s); pytest `test_bench_spin_deuteron.py::test_est_identity_to_the_bisection_bound` (+4 `test_est_*`), doctest *"bench row 9: …"*. Measured max \|diff\| = **1.138e-14** over 19 999 interior points, tolerance **2.31e-14** (the bisection bound 2.13e-14 plus margin); **the planned 1e-14 does NOT hold** (496 points exceed it) |
| 10 | **CLAS BONuS / EG1b database query** | T2 | hours | if machine-readable, the closest experimental analogue of the tagged topology at A = 2 | **blocked.** Harness `t2_bonus_spectator_shape.py` exists (0.06–0.08 s); pytest `test_bench_spin_deuteron.py::test_bonus_row_is_blocked_by_name` SKIPS naming the document: Tkachenko *et al.*, PRC 89 (2014) 045206 **Supplemental Material** (HTTP 401 without an APS login) — `../references/REFERENCES.md` §1 **#59**. Recorded instead: the CLAS Physics Database holds BONuS only as F₂ⁿ/F₂ᵖ (one measurement); Deeps (Klimenko 2006, F₂ᴺ × P(p_s), 115 tables) vendored, unwired |

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

   > **RESOLVED AND APPLIED, 2026-09-06** (`07_cw_sign_investigation.md`,
   > `../open_items/run_2026-09-06/phase_CW_numbers.md`). Every clause of this
   > item was confirmed. The mapping is **+1**, exactly and for every θ_k, f₂/f₀
   > and α_p; CW's −2 is the value of their angular factor at θ_k = 0, which
   > `A_zz^wf` already carries. `build_amp2` did lack the i^L phase, the tagged
   > S–D interference sign was inverted, and `b1_nuclear` was the sector in the
   > right. **Fixed with one line** in `src/core/tagged.cpp` `build_amp2` — the
   > phase `(-1)^floor(L/2)`, the observable relative part of φ_L = i^L ψ_L —
   > and `rad_[2]` / `Wave::psi()` were **not** negated, so ψ₂ = +W stays the
   > stored convention and the two tests asserting it still pass. Measured on
   > the fixed library: `max|A_zz^wf − CW Eq. (6.12)| = 8.881784e−16` over all
   > 280 × 96 = 26 880 grid cells on **both** deuteron controls, against
   > **2.740499e+00** before; TABLE II's three rows −1.937124 / +0.999313 /
   > +0.967340 against CW's −2 / +1 / +1 (they read +0.617024 / −0.318307 /
   > −1.656113 before). The shipped observable moved: ⁶Li A_zz^tag at
   > k = 0.20 GeV, acceptance-weighted at the YR high-acceptance optics,
   > **+0.8450 → −1.2069** (Hulthén β = 0.30) and **+0.4518 → −0.5191**
   > (VMC AV18). ⁷Li does not move at all (one L = 1 wave), and neither does
   > the spin-blind rate (≤ 4.0e−16 relative). The gate itself is rewritten on
   > the **AV18** control CW quote TABLE II for, per `07` §8, and carries a
   > regression guard. *(2026-09-23: `tagged.json`'s port gate against
   > `polligen` is restored in every block too — the sibling took the same
   > phase; `../open_items/run_2026-09-23/phase_A_port_gate.md`.)*
2. **George–Knutson** — technique misattributed at every site (see #8).

   > **APPLIED 2026-09-23** (`../open_items/run_2026-09-23/phase_B2_nuclear.md` §5, §8). Code:
   > `include/lipolgen/cluster_config.hpp` (three `///` blocks, line count kept);
   > `tests/test_cluster_config.cpp` T22's name, `MESSAGE` and comment. Docs
   > (G-1 … G-16): `PHYSICS_CHANNELS.md` [GK99] and two further η sites,
   > `00_in_tree_checks.md` (C-1 moved to Table E; E-18; §8; §10),
   > `03_data_nuclear.md`, `04_theory.md`, `02_data_nucleon_deuteron.md`,
   > `01_generators.md`, `../references/{REFERENCES,00_corpus,01_nucleon-deuteron-data,01_li-nuclear-data}.md`,
   > `../OPEN_ITEMS_SOLUTIONS.md` §11.2, `../USAGE.md` §9, and dated notes in three
   > run records. **Not applied:** G-17, the outward-facing Mäntysaari draft —
   > the author's, registry row 21.
3. **eSTARlight** — the claim line "benchmarks `coherent.hpp`'s channel rate"
   is false: eSTARlight is exclusive ρ/φ/J/ψ with 85–97 % of its rate at
   Q² < 0.1 where `coherent.hpp` generates nothing; it is a rate SCALE for
   exclusive VM, not a check on the M_X ≥ 1.2 GeV continuum channel.

   > **APPLIED 2026-09-23** (`../open_items/run_2026-09-23/phase_B3_chain_rc.md` §2, §6.2). The exact
   > sentence exists nowhere; its equivalents were corrected. Measured from
   > `estarlight_li6_q2_floors()`: J/ψ 85.2 %, φ 95.4 %, ρ 97.2 % of the rate is at
   > Q² < 0.1 GeV², and `make_config(channel="coherent").scenario.q2_min` = 0.7.
   > Sites: `include/lipolgen/coherent.hpp` (`EstarlightLi6Row` comment and the
   > ESTARLIGHT SCOPE block), `tests/test_coherent.cpp` (comment),
   > `validation/o5_a2_reach.py` (docstring), `00_in_tree_checks.md` D-4,
   > `01_generators.md` §1 and §2.7, `PHYSICS_CHANNELS.md` (the eSTARlight row and
   > [53]), `OPEN_ITEMS_SOLUTIONS.md` §11.1, `README.md`, `DEVELOPMENT_PLAN.md` (two).
   > **Not applied:** the Mäntysaari draft's line 67 (outward-facing; registry row 21).
4. **POLRAD ADGH entry** — scope overstated in five places, one stale
   number; the three numbers themselves reproduce (re-downloaded, re-driven).

   > **APPLIED 2026-09-23 — with a correction to this item** (`../open_items/run_2026-09-23/phase_B3_chain_rc.md`
   > §3, §6.3). Re-downloaded (sha256 `83703668…1ea7b`; the recorded `curl` now
   > needs a browser user agent) and re-driven with POLRAD's own `ffdeu`/`qunc8`:
   > rows 1 and 3 reproduce to every printed digit, **row 2 does NOT** —
   > σ_u = −1.0303e−03, not −1.001e−03 (+2.93 %), σ_q/σ_u = +0.0623, not +0.0642.
   > So "the three numbers reproduce" held for two. The D-3 row of
   > `00_in_tree_checks.md` now states what is gated (G_M(0) at 1e−7, the sign;
   > magnitudes 20 % only at the two x < 0.1 rows), and row 2 is corrected at
   > every site found by a tree-wide grep: `tests/test_rc.cpp` T8(c′) (value, ratio, comments — su/POLRAD
   > 0.9548, still inside the 20 % gate), `include/lipolgen/rc.hpp` (two comments),
   > `src/core/rc.cpp` (one), `polrad_transcription_check.md` (two), the 2026-09-02 and
   > 2026-09-03 phase records (one each), `run_2026-09-03/SUMMARY.md` (two), `run_2026-09-03/STATUS.md` (one),
   > `OPEN_ITEMS_SOLUTIONS.md` (two), `README.md` (one), the D-3 row of `00_in_tree_checks.md`, and `docs/PHYSICS_CHANNELS.md` (one, the polarised-QE row — plus the second half of `run_2026-09-03/STATUS.md`'s sentence — both missed by the Documents stage and corrected by its fix stage the same day) (+0.062; factors 67 / ≈ 255).
5. **`b1_default_li6.json`** is a self-pin dumped from LiPolGen's own C++ —
   label it so; it is not a reference of any kind.

   > **APPLIED 2026-09-23** (provenance field + README section; numeric blocks
   > byte-identical; sha256 d7bd8ce6… → 8373db8e…; `../open_items/run_2026-09-23/phase_B3_chain_rc.md` §4).
   > Both dump scripts now keep the label on regeneration
   > (`validation/dump_b1_default_li6.py` writes the `provenance` key,
   > `validation/dump_polligen_reference.py` writes the README section) — both
   > re-run into a scratch directory on 2026-09-23 reproduced every reference
   > file and `validation/README.md` byte for byte. Sites that called it a
   > reference now say self-pin: `PHYSICS_CHANNELS.md`, `USAGE.md` (`--b1-model`
   > table), `README.md`, and dated notes in the 2026-09-03 registry prose and
   > phase records.

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
| U-10 | ⁶Li elastic FF shape in the RC tail | F_q(0), F_m(0) on measured moments | the C0 SHAPE is a model band (now: the UVa density is a data-derived shape, and it lies OUTSIDE the [ho, vmc-ft] band on every priced tail number — `../open_items/run_2026-09-23/phase_B2_nuclear.md` §3.3) |
| U-11 | the spin bookkeeping | nothing external today | #1 above (wired 2026-09-23): the bookkeeping now has an external anchor, EPIOS Table II, passing at residual 0; the anchor is for an atomic-beam deuteron source, not for lithium |
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
