<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-10 — summary (one page)

Supervisor, 2026-11-17, on `claude/run-2026-10-26` from `b3f8858` (the end of
run 2026-11-03) to `03f7c7d` plus this file's own commit, the last of the
run. The author's instruction on the nine items run 2026-11-03 left,
verbatim: "1. go with your recommendation. 2. re-dump, for the others go
with your recommendation". So **J1 and J2 are the author's decisions**; **J3–J9
and every dated note and ruling under J1, J2 and J6 are the SUPERVISOR's
choices** under that delegation, recorded with their reasons in
`DECISIONS.md` and never to be quoted as the author's words. Commits: A
`dd7f4d4` · P1 `1c98f62` · P2 `d86fbf2` · P2b `8b72cd5` · P3 `8b45e3c` · P4
`2bd78e2` · P5 `c87b5ea` (board `931e0fb`) · V/F `03f7c7d` · S this file's
commit. Phase records `phase_P*.md`, `phase_V_verification.md`; the board
`STATUS.md`.

## What moved, and what it measured

| decision | what moved | measured consequence |
|---|---|---|
| J2 (P1), the author's | the polligen reference tables re-dumped in place, committed as one run wrote them | 3154 of 19928 numbers at the rounding level (`xsec.json` up to 9.36e−13 relative among values above 1e−6) and one manifest string; the cause is ONE dump input, `arccos(1/√3)`, one ULP higher than on 2026-09-23 (a libm change the likely, unproven cause); every port gate passes at its present rtol and uses less of it |
| J1 (P2), the author's | on an inclusive deuteron the kernel's rank-2 slot is `MillerB1` itself, not the ⁶Li mapping | the deuteron's b₁, b₂ × 3.0532, nothing else; HERMES on the kernel's own b₁ χ²/6 14.56 → 5.26 (= the headline); and every default inclusive deuteron run (24 of 24 plan × configuration; 4 before) refused at `x_max` 1.0 on the x = 0.955 row |
| J1's note of 2026-11-12 (P2b), the supervisor's | the deuteron's b₁ confined to the leading-order positivity domain of a spin-1 target, −(3/2)F₁ ≤ b₁ ≤ 3(F₁ − \|g₁\|), visibly (`b1_clipped_rate_fraction`, a banner line); the refusal stays as backstop | all 24 runs build, also at P_z = P_zz = P_e = 1; the confinement binds on 3.5–4.8 × 10⁻⁷ of the rate (rows x = 0.871, 0.955); HERMES unchanged; `--pol-sf` is read on every default inclusive deuteron run |
| J3, J4, J7 (P3) | a bench-gated test of the BONuS headline's robustness; a recorded in-grid sub-row of the HERMES b₁ row; exit 2 and `blocked (environment)` for a verdict set installed but not as vendored | 0 of 28 BONuS verdict blocks flip under the 321 × 321 table or bin averaging (closest block at 1.027 of its threshold); in-grid χ²(b₁) 1.89/3, χ²(A_zz) 1.18/3, the headline unchanged (5.26/6, FAIL) |
| J5, J6 (P4) | the two clip shares said to be cell-centre estimates, with a measured bracket; a HepMC3 file written without a pipeline's attributes says so in one line; the two extrapolation keys come from `run_provenance`; `HEPMC3_CONVENTION.md` edited once, sha re-pinned, T6 clause C25 | g₁ share (⁶Li, config 1) 3.864e−7 at the centre, 3.847e−7 / 5.10e−6 all / any corner; 69 npz and 69 CLI HepMC3 files byte-identical, the six example files gain one or two lines; T6 25/25 |
| J8, J9 (P5) | the LO development item stays filed; USAGE §8's throughput recipe is 2 M events per run | inclusive 0.605, tagged ⁶Li α 0.321, coherent 1.01 M ev/s (2026-11-17, a quiet but shared machine; run-to-run spread up to 6.5 %) |

⁶Li, ⁷Li, the tagged and coherent channels: byte-identical in arrays, `meta`
and HepMC3 through the whole run, except the example files' provenance
lines (J6). The report stays **15 pass / 8 fail / 0 blocked**; its only
changes are the HERMES sub-row (J4) and the T6 clause count (J6).

Gates at the end: C++ 441 / 441 cases (17,244,832
assertions), pytest 1535 passed / 156 skipped (six parts), bench subset
372 / 1 / 0, citations 1243 / 0 broken (records 1432 in 100 / 0), SPDX 136,
doc-link 112; `make_report --check` OK at 15 / 8 / 0; the out-of-tree
LHAPDF-off build, re-run after F: 0 warnings, 422 / 422 cases.

## Verification (V) and fixes (F)
Two read-only Opus lenses over `b3f8858..931e0fb` (physics and honesty; code
and regression) filed 6 + 7 findings; one skeptic per lens, told to refute
them, left all 13 standing at least in part (7 confirmed, 4 partly, 2 true
of the code but never reached or of no consequence). None is a wrong
behaviour on a default run. Phase F (`03f7c7d`) applied 10 and recorded 3 as
not applied, each as the supervisor disposed (`phase_V_verification.md`
carries the table):
* **Stale scope statements (rated major).** About ten user-facing passages
  (USAGE, PHYSICS_CHANNELS, header comments, three Python help strings) still
  said `--pol-sf` is read only under a beam helicity — false for the default
  inclusive deuteron since P2b. Each gained a dated clause on its own line;
  the code and `meta` were already right.
* **A wrong tally.** The re-dump moved **3154** numbers, not 3170: the count
  had taken 16 leaves that are NaN before and after as changed. No value or
  gate is affected; corrected in the records of both runs and by dated notes
  in both `DECISIONS.md`. The commit messages of `1c98f62` and of run
  2026-11-03's P5 keep 3170 and cannot be changed.
* **The one behaviour change** (the ruling of 2026-11-17 under J6): a
  caller's LHAPDF grid wrapped in `PositivityClippedG1` was reported as "not
  an LHAPDF grid" with NaN grid keys; the wrapper now answers for the
  backend it wraps (three `meta` keys, on a caller's wrapped object only;
  no shipped path builds one). A/B on the run list: 69 npz and 75 HepMC3
  files byte-identical.
* Four record / USAGE wordings corrected, one test renamed to what it pins,
  two comment clauses.
* **Not applied, the reason recorded:** the `pol_sf` reason string (literally
  true; the ruling of 2026-11-12 keeps the configuration-level rule); the
  empty attribute map (P4's choice stands — item 3 below); a manifest string
  in a record (a quotation of data, not a citation).

## What went differently from the plan
A reboot on 2026-11-12 wiped the session scratch (briefs, the P1/P2 raw
outputs and logs); the briefs were rebuilt from the transcript and the P2
record says which gate ran when. Workers were stopped by the usage limit in
P2 and P4 (twice); the supervisor ran the gates they had not reached, on
the same tree, and the records' tables say who ran what. The machine is
shared with other sessions since 2026-11-12: gates run 1.5–2 × slower, one
bench run was killed by its timeout and repeated, and in F the T6 row's
30 s runtime ceiling failed once under load (35.08 s) and passed on
re-run — nothing loosened. In P4 one bench run preceded the report's
regeneration and failed the committed-report test, which passes after it.

## For the author — exposed by the run, not decided
1. **The toy b₁ has no scale dependence.** `MillerB1` is a Q²-independent
   table used against an evolving F₁, so b₁/F₁ grows with Q² at fixed x
   (0.54 → 6.99 at x = 0.955 from Q² = 2.5 to 1897 GeV² on the deuteron) —
   the cause of every high-x rank-2 violation. *Recommendation:* read the
   toy as the ratio b₁/F₁ at a reference scale (a scale-independent A_zz) —
   as its own decision: it moves every rank-2 number of ⁶Li at every Q².
   Until then the deuteron's b₁ above x ≈ 0.85 is the positivity domain's
   edge, not a model (b₁/F₁ 0.56 at x = 0.85, 0.26 at 0.955; not monotone).
2. **⁶Li was left exactly as it was** (J1's list): three cells of the
   default window sit outside the joint half of the domain (x = 0.955,
   Q² 1536–1897; \|g₁\|/F₁ + b₁/(3F₁) up to 1.046; 2.4–2.7 × 10⁻¹¹ of the
   rate; no run refused), and ⁶Li's b₁ at x ≥ 0.87 is built from
   deuteron-toy values the deuteron's own domain excludes. Item 1 would
   settle both; otherwise extend the confinement to ⁶Li.
3. **Smaller, each recorded:** `--b1-model cdks` on the deuteron would now be
   well defined and stays refused; `pol_sf` is reported read on every default
   inclusive deuteron run, also where no cell binds (a per-run rule is
   possible); `set_run_attributes` with an explicitly empty map writes
   neither attributes nor the absence line; the g₁ clip share's any-corner
   bracket is 13 × its centre value (one row, x = 0.794, binds inside its
   cells); the shipped F₁ᵈ is 10 % below to 28 % above what HERMES's own
   columns imply even inside the grid; the BONuS headline's closest block
   sits 2.7 % from its threshold; the dump's `arccos(1/√3)` input is
   libm-sensitive (a literal would make the dump byte-stable); the T6
   measured row's wall-clock ceiling (30 s) is load-sensitive on a shared
   machine, and throughput there is "quiet", not idle.
