# Run 2026-09-06 — what is still open and does not wait on the author

Follow-on to the 2026-09-03 run (`../run_2026-09-03/`, 12 commits
`a94fd6e..5af0427`, nothing pushed). Scope is `../run_2026-09-03/SUMMARY.md`'s
"What the run did NOT establish", restricted to what is (a) offline-doable and
(b) independent of the 25 registry decisions — i.e. work that makes those
decisions *decidable* or closes an item outright, never work that takes one.
Progress in `STATUS.md` next to this file.

## Baseline (measured at the start, commit 5af0427, tree clean)

401 doctest cases / 17 240 286 assertions / 1 skipped (T9, `doctest::skip(true)`,
reported by count only); pytest 926 passed / 112 refused-base skips; docs gate
1206 references (strict) / 96 ranges / 7 external / 0 broken, +19/115
(SPIN32_FINITE_GAMMA.md), +8 external (PYTHIA_BRIDGE.md); SPDX 101/101.

## Ground rules

Those of `../run_2026-09-03/PLAN.md`, unchanged, plus the run's lesson stated
as a rule: **write the claim you measured, with its window and exceptions, at
every site that carries the old one; a knob that did not run is refused or
labelled, never recorded as run; a count that does not reproduce from committed
state is not published.** No default moves; every decision stays in the
registry (`../run_2026-09-03/AUTHOR_DECISIONS.md`); nothing is pushed or sent.

## Phase A — make the batch decidable: the three unpriced costs, and the tables

- **A1** Registry row 3, option (iii): ONE R hook threaded through both the
  numerator (`Li6ConvolutionOptions`) and the denominator (`InclusiveKernel`)
  of the ⁶Li tensor weight. Implement opt-in, measure against (i)
  `r_sigma_lt` and (ii) `r1998`, price it. Default unchanged.
- **A2** Row 17: the RC band's low-x anchor alternatives 0.266 / 0.113 (the
  panel's own values at x = 0.00226 / 0.00966) against the shipped 0.30 —
  re-evaluate the band half-width on A_zz at the standard points, publish the
  three-row table, pin it. The knob exists (`--rc-delta-low-x`); this is
  measurement.
- **A3** Row 20 / §15.5 D13: honour `--pzz` on `helicity-flip`. Implement as
  an opt-in fill mode, measure the effect on every observable that reads the
  fill, price it. Default (the max-entropy ladder) unchanged.
- **A4** The ⁶Li b₁ tables in `../run_2026-09-02/phase_D_numbers.md` were made
  with `ToyF2` and are marked "NOT covered by the lift until rerun". Regenerate
  them under the configuration the lift was granted on and publish them as the
  quotable **band**: `--b1-unpol {toy, ct18nlo, mstw}` × `--unpol-sf {toy,
  ct18nlo, mstw}`. Registry row 6 (B11): recompute the `A_zz(Born)` column with
  the shipped kernel and publish it BESIDE the original with a dated correction
  — the original stays, so no decision is taken.

## Phase B — the physics that is bounded, not computed

- **B1** The **tensor** fraction of the s-/p-peaks is "unknown, not zero".
  Bound it the way `qe_tensor_scale` bounds the quasi-elastic tail: a stated
  scale whose value 1 is the elastic tensor fraction, priced at the standard
  points, opt-in, default 0.
- **B2** `RcTailModel::PolradFull` — Eq. (18) + Appendix B + Eq. (A.4) — is
  `NOT IMPLEMENTED` and T9's half is `doctest::skip(true)`. The transcription
  was checked line by line against the POLRAD LaTeX and the ADGH FORTRAN in
  the 2026-09-02 run (`polrad_transcription_check.md`), so the formulas are in
  this tree. Implement it opt-in, gate it against `TPeak` in the limit where
  they must agree, un-skip T9, and measure where the two differ.
- **B3** Q(⁷Li) is "quoted from memory and absent from this tree". Source it
  (network is available: Stone's moment tables / the original measurement),
  put the constant in `constants.hpp` with its provenance, and make the A = 7
  quadrupole gate a real, cited test. This does NOT implement b₁(⁷Li) — that
  waits on row 3.

## Phase C — the efficiency the O5 verdict is missing

- **C1** No decay-lepton reconstruction efficiency exists in this tree. The
  sibling `../PolarizedLithiumSim/fastsim/` and `tools/fullsim` carry the
  project's own detector assumptions. Derive a BOUND from them — stated as a
  bound, with its source — and re-run `validation/o5_a2_reach.py` with it as
  a fourth ladder row. If the sibling has nothing usable, say so and stop.

## Phase D — the gate residue class

- **D1** `--record-ranges` blesses content at recording time; a range recorded
  onto the wrong block stays "0 broken". Require, at recording, that the
  citing sentence's quoted text or named symbol occurs inside the block, and
  refuse to record otherwise.
- **D2** The dated run records under `docs/open_items/run_*/` are un-gated and
  drift with every line shift. Add a relaxed gate mode for them: each
  `path:line` must land on a non-blank line sharing a token with the citing
  sentence; fix what it finds or annotate as historical.

## Phase E — close-out

Whole-diff review, `STATUS.md` final, registry rows added ONLY where a phase
surfaced a genuinely new decision, and the batch's numbers re-verified.
