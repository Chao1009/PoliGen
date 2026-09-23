<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-09-23 — status board

This run, as the record counts it, is four pieces of work since `b647c16`:
the thorough review of 2026-09-16 (committed `239ac88`), the reference search
(committed `d735abc`), and — uncommitted in the working tree at the time of
writing, on base `d735abc` — the polligen port gate (phase A), the M1
benchmark harnesses (phases B1–B3) and this Documents stage with its fix stage. The one-page
account is `SUMMARY.md`. **No physics default, registry status label,
constant or reference-JSON number moved in any of them.** The decision
registry is still `../run_2026-09-03/STATUS.md` (27 rows); cite a decision by
its row number there.

## Phase board

| phase | state | commit | what it did, with its window | record |
|---|---|---|---|---|
| review 2026-09-16 | DONE | `239ac88` | six lenses over `b647c16`: **78** findings, **28** serious, **16** confirmed, **15** applied; **34** of 40 minors applied; the eight residuals closed. Every CLI default, plan and reference gate bit-identical on **207** configurations hashed against a from-scratch `b647c16` build (threaded = serial at 1/4/8 threads). One behaviour change, on a non-default configuration: a TILTED fill on a tagged channel with a spin-dependent struck-cluster term is refused | the commit message; `../run_2026-09-06/STATUS.md` (the two costs priced 2026-09-16) |
| references 2026-09-16 | DONE | `d735abc` | `docs/references/REFERENCES.md`: **266** references, **56** PLEASE DOWNLOAD (§1, numbered 1–58; #43 resolved, #27 free-not-fetched); six domain searches `01_*.md`; the sibling corpus `PolarizedLithiumSim/refs/` at **200** PDFs / **244** `refs_dict.json` entries (committed there since: the sibling's HEAD is `4726812` "add references", 2026-09-23 08:04, which touches `refs/` only) | `../../references/REFERENCES.md` |
| A — port gate | DONE (uncommitted) | — | `validation/reference/tagged.json` dumped from `polligen` again in every block (polligen code as of the sibling's `c0f86a8`, identical at its HEAD `4726812`, which differs only under `refs/`); the two spin-1 model blocks moved by ≤ **1.33e−13** (⁶Li) / **2.28e−13** (deuteron) relative, every rtol-1e−12-gated entry by ≤ 6.1e−15; the dump script refuses a pre-fix `polligen` | `phase_A_port_gate.md` |
| B1 — spin and deuteron | DONE (uncommitted) | — | plan rows 1, 3, 9, 10: PASS (residual 0 at 1e-12), FAIL recorded (HERMES χ²/6 = 5.26 … 22.67), PASS (1.138e-14 at 2.31e-14), BLOCKED (BONuS supplemental) | `phase_B1_spin_deuteron.md` |
| B2 — nuclear inputs | DONE (uncommitted) | — | plan rows 4, 5, 6, 8: PASS (NMC 4.6269/4 and 26.941/15), FAIL recorded (C0 zero +0.2713 fm⁻¹ above the data band), BLOCKED (RFY / ADNDT 14), relabel (three comment blocks) | `phase_B2_nuclear.md` |
| B3 — chain and RC | DONE (uncommitted) | — | plan row 7 wired and BLOCKED (DJANGOH ran with the elastic tail off); §5.3 eSTARlight code site; §5.4 POLRAD row 2 re-driven — stale (−1.0303e−03 not −1.001e−03); §5.5 `b1_default_li6.json` labelled a self-pin | `phase_B3_chain_rc.md` |
| Documents | DONE (uncommitted) | — | every requested edit of the four records applied or, where outward-facing, put in the registry (§ below); `BENCHMARK_PLAN.md` status column; registry prices; state board; one `validation/benchmarks/README.md`; SPDX gate extended to `validation/benchmarks/*.py` | this file, `SUMMARY.md` |
| Documents — verification and fix stage | DONE (uncommitted) | — | four adversarial lenses: numbers, gates and defaults held; 16 wording/bookkeeping findings (F1–F16) fixed the same day — stale POLRAD +0.064 at two more sites, the ⁶Li P_zz absence restated with its searches, George–Knutson / eSTARlight residues, tallies re-based on `d735abc`, sibling HEAD `4726812`, uniform harness exit convention (0 ran / 2 broken), runtime ranges, the dump's sign guard moved to a pre-flight before any write. Suites 416 / 17 242 712 / 0 and 1134 / 153 / 0; strict 0 broken, `--audit-ranges` 0 REFUSED, `--records` 799 / 49 / 0 broken, SPDX 112 | `SUMMARY.md` §Verification |

## What the Documents stage did not apply (or applied differently), and where it went

| requested edit | why not | where it is instead |
|---|---|---|
| B2 **G-17** — the Mäntysaari draft's "anchored on the *measured* … η" | the draft is outward-facing; the author decides (registry row 21) | `../run_2026-09-03/AUTHOR_DECISIONS.md` §B25 box (with the letter's second site, its line 531, which the record did not list) |
| B3 §6.2 — the same draft's line 67, "coherent-rate baseline" | outward-facing, same reason | the same §B25 box |
| B3 §6.1 — `eic/InclusiveDjangohSamples`' README sentence | external repository; the record asks for a note only | `phase_B3_chain_rc.md` §1.1 |
| B2 **F-6** — the `LI6_FF_HO_A_FM` block's "Fit q_0 over [2.9, 3.3] …" sentence | **applied, with different wording** (a deviation, not an omission; recorded by the 2026-09-23 fix stage): the requested "Two data-derived numbers sit BELOW [2.9, 3.3]: …" became "T11 gates q_0 over [2.9, 3.3] fm^-1; two data-derived numbers sit BELOW it: the UVa FB zero 2.694 and Li71's minimum 2.828 fm^-1 (t3_li6_charge_ff_fb.py, 2026-09-23)" — content equivalent (same two numbers, harness and date), line count kept | `include/lipolgen/rc.hpp:551-552` |

No requested edit would have moved a default or taken a decision, so none
had to be converted into a priced option; the prices the build stage
produced were filed as dated boxes (row 2 in §B2(c); the C0-zero price under
§B15, no row existing for `LI6_FF_HO_*`).

## Gates at the end of the fix stage (2026-09-23)

See `SUMMARY.md` §"Tallies" for the numbers and their baselines.
