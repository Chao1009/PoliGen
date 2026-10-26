<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-26 — status board

The plan is `PLAN.md`; the one page is `SUMMARY.md`; the ranking this run
was asked for is `../../benchmarking/08_channel_priority.md`. The decision
registry is still `../run_2026-09-03/STATUS.md` (27 rows); cite a decision
by its row number there; this run's priced options are in
`DECISIONS_PRICED.md`. **No physics default, registry status label,
constant, tolerance or reference-JSON number moves in this run; no decision
is taken.** Work is on branch `claude/run-2026-10-26`, committed and pushed
phase by phase.

## Phase board

| phase | state | commit | what it did, with its window | record |
|---|---|---|---|---|
| 0 — environment, ranking, plan | DONE | `82b3237` | the dependency stack built with the project's script (HepMC3 3.3.0, LHAPDF 6.5.5 + the three sets, PYTHIA 8.317; 439 MB, 46 min on 4 cores); LiPolGen built with every tier on; **C++ 416/416 cases, 17 242 712 assertions, 0 failed; pytest 1269 passed / 158 skipped / 0 failed**; the first report with the LHAPDF sets (5 pass / 8 fail / 5 blocked before the PYTHIA fix); PYTHIA 8.317's cwd-relative python loader found and worked around in `_report.py` (`import_pythia8`), T4 PASS at worst 0.076 on 8.317, T6 23/23 | `phase_0_environment.md`, `../../benchmarking/08_channel_priority.md`, `PLAN.md` |
| 0r — the ranking and plan reviewed by two adversarial Opus readers | DONE | `87787dd` | physics lens 24 findings (8 major), honesty/feasibility lens 25 (5 major); 49 in all, 45 applied in the rewrite of `08_channel_priority.md` and `PLAN.md`, 4 applied in part (below); the DJANGOH worker, started on the pre-review brief, was stopped before it built anything and relaunched on the reviewed one | this file (below) |
| A — baseline and price table | TODO | — | | `phase_A_baseline.md` |
| B1 — row 10 BONuS | TODO | — | | `phase_B1_bonus.md` |
| B2 — row 6 RFY | TODO | — | | `phase_B2_rfy.md` |
| B3 — row 7 DJANGOH with the elastic tail on | DONE (an adversarial reader of it is still running) | (this commit) | DJANGOH 4.6.21 built from `e356b84` (gfortran 11.4, LHAPDF 6.5.5 + cteq61 = 10150); the eight ePIC numbers reproduced (\|pull\| ≤ 1.74); tail isolated (IEL31+32+33, IEL2 = 0, WMIN 0 → M_p, x_l ≤ 0.1096 / 0.5518 / 0.9249 / 0.95, y ≤ 0.985, Q² bins by subtraction since the IEL channels ignore Q2MAX); row 7 BLOCKED → PASS under the tolerance fixed in the header before the runs: `polrad-full` R 1.0199 / 1.0261 / 1.0219 / 1.0098, `t-peak` 0.6118 / 0.8968 / 1.0112 / 1.0067; straight from the quadrature 1.0004 / 1.0010 / 0.9995 / 0.9866; the tail table adds +1.9 … +2.5 %; `test_bench_chain_rc.py` re-pointed to the history block; plan row 7, the README rows, `08_channel_priority.md` row F and REFERENCES updated | `phase_B3_djangoh.md` |
| C1 — R data; HERMES F₂ᵈ; COMPASS g₁ᵈ | TODO | — | | `phase_C1_r_world.md` |
| C2 — positivity and the Soffer bound | TODO | — | | `phase_C2_positivity.md` |
| C3 — Abbott t₂₀ in the RC control | DEFERRED WITH DESIGN (review) | — | no generation path reads a deuteron elastic form factor; a deuteron `TabulatedSpin1FF` is a registry candidate | `DECISIONS_PRICED.md` |
| C4/C5 — Wiringa–Schiavilla; Boeglin | TODO | — | | `phase_C4_wiringa_schiavilla.md`, `phase_C5_boeglin.md` |
| D — the least-verified inputs read | TODO | — | | `phase_D_inputs.md` |
| P — figures for the 18 existing rows; `make_report.py --plots` | TODO | — | | `phase_P_figures.md` |
| E — documents, report, priced decisions | TODO | — | | this file, `SUMMARY.md`, `DECISIONS_PRICED.md` |
| F — verification and fix stage | TODO | — | | this file (Verification) |

States: TODO / IN PROGRESS / DONE (with commit hash).

## The review of the ranking and the plan (phase 0r)

Two Opus readers, read-only, 35 minutes, 223 tool uses, each told to refute
the two documents against the tree (the workflow `critique-channel-ranking`,
its transcript in the session's workflow directory). **Physics lens: 24
findings (8 major, 12 minor, 4 nit); honesty/feasibility lens: 25 (5 major,
18 minor, 2 nit).** The majors, and what was done:

| # | finding (abridged) | disposition |
|---|---|---|
| P1 | the fill model (ingredient B) ranked second on a misreading: the tensor plans build the (+0.6, −1.2) pair by default; `--pzz-mode` is read by helicity-flip alone; EST's P_zz ≥ 0 concerns that fill only | applied: B demoted to a beam-assumption line; §0(c) rewritten |
| P2 | F₂ cancels exactly in the cos 2φ amplitude; only R moves it (≤ 18 % sibling, −3 … −8 % generator-side); ×0.80 / ×1.25 are σ and A_zz | applied: rank 1 split into cos 2φ (claim 1) and A_zz/b₁ (paper 2); ingredient A restated per observable |
| P3, H2 | the rubric not applied literally (seven cells); the d–p control outranks tagged ⁶Li α when it is | applied: levels 0–3 defined per criterion; the table rescored literally (28/25/25/23/20/18/19); the control placed last by a stated rule R1; ties by R2 |
| P4, H6 | the deuteron wave function priced at 6.6 % where A_zz^tag moves ×2.3–4.8; −2.03 % belongs to the d–p control | applied |
| P5, H16 | the positivity/Soffer gate not well-defined for the generator's objects (gluon distributions needed; T-13 is spin-3/2 on a zero sector; b₃/b₄ higher twist) | applied: C2 re-specified as the Detmold–Shanahan n = 2 moment bound under declared assumptions; T-13 and b₃/b₄ dropped; spin-1 positivity only if sourced |
| P6, H21 | Abbott t₂₀ against a deuteron form factor no generation path uses (the d–p control has `rc_tail` ≡ 1; the inclusive deuteron run throws) | applied: C3 deferred with design, filed as a registry candidate |
| P7 | the R data chosen miss all four sweet spots; NMC's R covers them; E140/R1990 are `r1998`'s own inputs | applied: NMC R primary, H1 = low-x extrapolation, E140/R1990 = closure, the price at the sweet-spot cells |
| P8 | "no RC statement for a φ-dependent tensor observable" — the letter carries a leading-log ISR bound ≤ 2.9 % (WP4) and `RcScope::TensorAll` prices the omission | applied |
| H1 | B3's design: IEL2 is the elastic peak, not a tail; the inelastic pieces would dominate; the tolerance would be set after the number | applied: the elastic tail isolated (IEL31–33 on / INC off, or an on/off pair), IEL2 = 0, W_l² matched, tolerance from the Eq. (44) accounting written before the run; the worker relaunched |
| H3 | the §0 PASS list wrong (NMC ⁶Li/D blocked on 09-27; r_ch is a passing data comparison; HepMC3 passed by re-scope; Deeps p-values per wave) | applied |
| H4, P14 | A3 priced σ, not the ranked observables; most knob cells refused or not read | applied: A3 re-specified per channel and observable with the knob-provenance labels |
| H5 | RFY prints a₀, not r_mag — the radius rule would be a derivation with the tree value previewed | applied: the measured F_M(q) points are the primary comparison; the radius rule recorded only with the formula declared first, else retired |

Applied in part: H8 (C3 = C1 weighting stated as this document's choice and the order shown robust to it — the weighting kept); P9/H12 (the O5 J/ψ reach labelled MARGINAL and as another process — kept in the dossier for context); H19 (the M2 clause made a proposal for the registry — kept in §5.1 as a proposal); H20 (the author's 2026-10-26 request cited as the release of D3's "left to the author"; the refs/ filing and licence question listed for the author — the fetches kept). Every numeric correction (triton polarizations 0.86/−0.028 and the OPEN P_n gate; the Δ scale 10⁻²; δA single-fill vs two-fill; c₂ ≈ ×6; Deeps p-values; the POLRAD wording; COMPASS 2017 already vendored; the BONuS binning; the `data/ff/` slot) was applied.

## What the run did not apply, and where it went

Filled in by the integration stage.

## Verification

Filled in by the fix stage.
