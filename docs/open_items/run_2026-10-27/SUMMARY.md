<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-27 — one page

**What this run is.** The author took the four options run 2026-10-26 had
priced (`../run_2026-10-26/DECISIONS_PRICED.md`), in the words quoted in
`DECISIONS.md`: D4 the Δ normalisation (first "keep the factor 2, is this
more common?", then, on the literature reading, the Jaffe–Manohar ½), D5
`PolradFull` as the `rc_tail` default, D6 the VMC AV18+UX cluster waves as
the default, D7 the grid structure-function backends (CT18NLO, NNPDFpol1.1)
as the defaults with the toys opt-in. Each decision was researched where it
asked a question (R1), designed with an adversarial check (D4–D7), implemented
by one Opus worker, verified by a second and fixed by a third (I4–I7), and
committed one per decision; then four adversarial lenses read the whole run
and one fixer applied every finding (V, F). Branch `claude/run-2026-10-26`,
on top of `ca9b23a`: `fbc24f2` (R, D), `9779edd` (I4), `cb75ec7` (I5),
`c52c164` (I6), `7fd441c` (I7), `97c521c` (V, F), this page after it; all
pushed, no pull request opened. **The author's four defaults are the only
physics that moved. No constant, tolerance, registry label, reference JSON,
vendored file or hashed design block was touched; every number that depends on
a moved default was re-taken with a dated comment naming the decision, or the
test passes the old option explicitly where its purpose is to pin the old
physics; every benchmark that disagrees is recorded, none tuned. Nothing was
published or sent.** The board is `STATUS.md`; the whole-run verification is
`phase_V_verification.md`.

## What moved, and what it measured

- **D4 (`9779edd`).** `DELTA_COS2PHI_NORM` = ½ at the two evaluation sites
  (`src/core/xsec.cpp`, `src/core/asymmetries.cpp`), its label tied by a
  `static_assert`. Every Δ-induced cos 2φ amplitude the generator produces
  exactly halves; everything Δ-free, the θ_S = 0 output among it, is bit
  identical. External Δ numbers (the bag −0.012 α_s, the lattice bound) now
  enter unconverted; row 13 (`t5_delta_moment_bound`) keeps PASS under its
  now-consistent identification. The port gate reads polligen's coefficient
  off `validation/reference/xsec.json` against a literal ½, so a sibling
  re-dump in either convention is caught. The sibling's `asymmetries.py`,
  `evgen/polligen/xsec.py` and the letter's Eqs. (1), (3) still carry the old
  coefficient: until they change, the letter's printed amplitudes encode
  2 × C_BAG. R1's finding: eight sources carry the ½ (Jaffe–Manohar 1989 read
  as CTP#1712 Eq. (16), Sather–Schmidt 1990, Ma–Wang–Zhang 2013,
  Detmold–Shanahan 2016, the LOI's Eq. (6) …); the doubled coefficient appears
  only in the LOI's underived Eq. (9) and one reprint.
- **D5 (`cb75ec7`).** `RcOptions::tail_model` = `PolradFull` in one place,
  the CLI default read from it. `--rc off` and every explicit-model run bit
  for bit; row 7 (DJANGOH elastic tail) PASS, relabelled "polrad-full
  (shipped), t-peak (option)"; T6 16.3 s on the default. Every t-peak number
  in the tree labelled with the default's beside it (the b₁-band lead at x = 0.01
  ×768 → ×1152; the QE shares 22 / 73 / 99.9 → 39 / 99.8 / 100 %).
- **D6 (`c52c164`).** `CLUSTER_WAVE_DEFAULT` = `VmcAV18` on the run surface;
  the channel factories, `BreakupOptions` and the polligen port gate keep the
  Hulthén forms, so every reference JSON stands. ⁶Li tagged A_zz −1.2069 /
  −0.9533 → −0.5191 / −0.1993, the recoil-proton tag fraction 0.0247 → 0.0342, ⁷Li ⟨k⟩
  0.2885 → 0.1861 GeV. The α–d asymptotic D/S ratio decided the physics:
  Hulthén +0.070 (the wrong sign), VMC −0.048, measured −0.025 ± 0.012 (VMC
  itself 1.9 σ off; CD-Bonn exists in the tree for the deuteron only). Deeps
  keeps its declared Hulthén verdict, AV18 recorded. `LI6_CLUSTER_POLARIZATION`
  0.811228 and `LI6_B1_RANK2_TRANSFER` 0.921947 were KEPT (the VMC wave
  gives 0.887 / 0.983); the inconsistency is the author's first item.
- **D7 (`7fd441c`).** `UNPOL_SF_DEFAULT` = `Ct18Nlo`, `POL_SF_DEFAULT` =
  `NnpdfPol` in one place; the grids attached lazily when a `Pipeline` is
  built and FROZEN at their knot edges in code (LHAPDF `nearest` set on the
  object plus a clamp), so no `lhapdf.conf` moves a default number (the
  freeze was the supervisor's choice: the design had proposed continuation,
  its check the freeze, after a conf difference moved the pre-freeze ⁶Li rate
  by +5.03 %). ⁶Li inclusive σ ×0.83862 (591846.17 → 496334.16 pb), coherent
  ×0.75541, ⁷Li ×0.83604, A_zz ×1.19243; 39.23 % of the default ⁶Li rate sits
  on the frozen edge below CT18NLO's Q² floor of 1.68 GeV² (48.37 % coherent).
  A missing set is refused by name before any LHAPDF call and blocks the four
  T1 rows "blocked (environment)". The development item carries its measured
  cost: at LO the shipped pair breaks |A₁| ≤ 1 from x ≳ 0.83 on the free
  neutron and x ≳ 0.87 on the deuteron and proton; `--x-max 0.95` prevents the
  refusal, not the violation.
- **The report** (`docs/benchmarking/REPORT.md`, regenerated after each
  phase): **23 harnesses — 14 pass, 9 fail (recorded), 0 blocked** (run
  2026-10-26: 13 / 10 / 0). The one verdict that moved, `t1_nmc_f2d_over_f2p`
  FAIL → PASS, did so because its declared rule names the shipped default,
  which is now CT18NLO (215.5/189); rows 7 and 13 keep PASS relabelled. F's
  regeneration moved three printed cut strings only (Q² ≥ 1 → the 1.677 the
  numbers were computed on); every χ², ndf and p is unchanged.
- **Gates on the final tree** (`phase_V_verification.md` §7): C++ 428 / 428
  (17 244 386 assertions), pytest 1515 passed / 155 skipped / 0 failed, the
  bench subset 368 / 1 / 0, the strict citation gate 1241 / 95 / 7 / 0 broken
  / 24, R5 209 ranges on 179 blocks 0 refused, R6 1432 citations in 79 records
  0 broken, SPDX 136 files, `make_report.py --check` OK.

## What it did NOT establish

- Whether the LO evaluation of the NLO sets is acceptable for production: the
  |A₁| violation above x ≈ 0.83 is measured and filed, not cured.
- Whether the T1 "shipped default" rule is to be read as I7 read it (follow
  the default: 14 / 9) or as I6 read it for Deeps (keep the declared option:
  13 / 10); both readings are recorded, the choice is the author's.
- The ⁶Li cluster constants on the VMC wave: kept, not re-derived.
- The EPPS21 ⁶Li ratio and R1998 as defaults: not moved; the new default is
  free-nucleon CT18NLO with the toy R, neither the letter's grid chain nor its
  published toy footing.
- Three measured numbers stand unreproduced or unlogged beside their re-takes
  (I5's 5.18398e−04, I7's 5000-event rc_tail pair, USAGE §8's throughput).

## For the author

`STATUS.md` "For the author" lists all 17 items with their record sections.
The ones that change physics output, with the supervisor's recommendation
where it has one:

1. Re-derive `LI6_CLUSTER_POLARIZATION` / `LI6_B1_RANK2_TRANSFER` on the VMC
   wave (0.887 / 0.983), so the default inclusive ⁶Li run and the default
   tagged run use one wave. Recommended yes; the Hulthén values were kept
   only because they are not among the four decisions.
2. Decide the LO development item's fate for production: NLO coefficient
   functions, a positivity clamp, or the toys for x > 0.8.
3. EPPS21 ⁶Li ratio by default for lithium beams, and R1998 in place of the
   toy `r_sigma_lt`.
4. The T1 "shipped default" rule reading (item 5 of the list).
5. The sibling sites (`asymmetries.py`, `evgen/polligen/xsec.py`, the letter's
   Eqs. (1), (3)): whether the letter means the amplitude or Δ to halve.

The rest of the list: the t5 TIGHT edge and the grid row's missing bound,
the HepMC3 provenance gap (changes the sha256-pinned `docs/HEPMC3_CONVENTION.md`),
the two CLI plans refused at the default, the harness exit code 2, the CMake
`FATAL_ERROR`, a finite `cluster_vmc_mc_sigma` bound, the ⁷Li FSI ceiling,
a VMC port-gate block, I5's two changed readings, the knob asymmetry across
D5 / D6 / D7 (stated, not unified), and the unlogged throughput. Outside the
tree: the papers fetched for R1 sit in the session scratchpad (`refs/`) for
filing into `PolarizedLithiumSim/refs/`; the licence questions of run
2026-10-26 stand.
