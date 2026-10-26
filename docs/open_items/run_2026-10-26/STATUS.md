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
| 0r — the ranking and plan reviewed by two adversarial Opus readers | IN PROGRESS | — | physics lens; honesty and feasibility lens | this file (below) |
| A — baseline and price table | TODO | — | | `phase_A_baseline.md` |
| B1 — row 10 BONuS | TODO | — | | `phase_B1_bonus.md` |
| B2 — row 6 RFY | TODO | — | | `phase_B2_rfy.md` |
| B3 — row 7 DJANGOH with the elastic tail on | TODO | — | | `phase_B3_djangoh.md` |
| C1 — R data; HERMES F₂ᵈ; COMPASS g₁ᵈ | TODO | — | | `phase_C1_r_world.md` |
| C2 — positivity and the Soffer bound | TODO | — | | `phase_C2_positivity.md` |
| C3 — Abbott t₂₀ in the RC control | TODO | — | | `phase_C3_abbott.md` |
| C4/C5 — Wiringa–Schiavilla; Boeglin | TODO | — | | `phase_C4_wiringa_schiavilla.md`, `phase_C5_boeglin.md` |
| D — the least-verified inputs read | TODO | — | | `phase_D_inputs.md` |
| P — figures for the 18 existing rows; `make_report.py --plots` | TODO | — | | `phase_P_figures.md` |
| E — documents, report, priced decisions | TODO | — | | this file, `SUMMARY.md`, `DECISIONS_PRICED.md` |
| F — verification and fix stage | TODO | — | | this file (Verification) |

States: TODO / IN PROGRESS / DONE (with commit hash).

## The review of the ranking and the plan (phase 0r)

Filled in when the two readers return.

## What the run did not apply, and where it went

Filled in by the integration stage.

## Verification

Filled in by the fix stage.
