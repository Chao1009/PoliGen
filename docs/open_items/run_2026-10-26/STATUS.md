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
| A — baseline and price table | DONE | uncommitted (phase E tree) | tallies C++ 416/416 (17 242 712 assertions), pytest 1269 passed / 158 skipped / 0 failed, against 2026-09-27's 13 and 39 LHAPDF failures; the first `--compare` with the LHAPDF sets, 38 changes, one line each; T4/T6 re-run after the PYTHIA workaround (PASS, 0.076; 23/23), the 7 / 8 / 3 tally inferred for the 18 harnesses (with B1–B3 the same 18 read 9 / 9 / 0, measured in phase E; 13 / 10 / 0 with the five new rows); the default-input price table of the ranked observables at config 1, N = 200 000, seed 20261026 (`validation/price_defaults.py`, `price_defaults.png`): R −15.2 / −14.7 / −7.8 / +4.2 % on the cos 2φ amplitude over the four super-bins x = 0.011 / 0.028 / 0.071 / 0.141, F₂ ≤ 0.23 %; A_zz ×1.2524 / ×1.2604; A_∥ ×0.532 / ×0.514 under NNPDFpol1.1; A_zz^tag ×0.430 / ×0.209 under VMC; c₂ ×0.270 / ×0.158 at the measured quadrupole; `test_price_defaults.py` (1) | `phase_A_baseline.md` |
| B1 — row 10 BONuS | DONE | uncommitted (phase E tree) | the BONuS Supplemental Material vendored as numbers (100 R_D/S blocks + the DIS F₂ⁿ/F₂ᵖ table, the original's sha256 recorded); row 10 BLOCKED → **FAIL, recorded** under a rule declared before any number: the Hulthén 40 mb `GlauberFsiWeight` over its backward-window mean inside the data's 1σ flat band in 4/12 (4 GeV) and 0/16 (5 GeV) DIS blocks; the backward data are not flat (χ² 739.4/84) and fall toward 90° as the weight does (χ² 226.3/84 against its shape; 20 mb 403.3); `test_bench_bonus.py` 26 passed | `phase_B1_bonus.md` |
| B2 — row 6 RFY | DONE | uncommitted (phase E tree) | row 6 BLOCKED → **PASS**: RFY Table I's seven ⁶Li F_M² points against `HoSpin1FF`'s unfitted F_mag², χ² = 13.35/6 with one free scale s = 1.419 ± 0.206 (p = 0.0378, against p ≥ 0.01 declared in a hashed design block); s = 1 recorded 17.49/7 (p = 0.0145); ⁷Li not compared; the radius rule retired (ADNDT 14's lithium 'magnetization radii' are √(5/2)·a₀); `test_bench_rfy.py` 16 | `phase_B2_rfy.md` |
| B3 — row 7 DJANGOH with the elastic tail on | DONE, verified | `372b915`, `2d89a64` | DJANGOH 4.6.21 built from `e356b84` (gfortran 11.4, LHAPDF 6.5.5 + cteq61 = 10150); the eight ePIC numbers reproduced (\|pull\| ≤ 1.74); tail isolated (IEL31+32+33, IEL2 = 0, WMIN 0 → M_p, x_l ≤ 0.1096 / 0.5518 / 0.9249 / 0.95, y ≤ 0.985, Q² bins by subtraction since the IEL channels ignore Q2MAX); row 7 BLOCKED → PASS under the tolerance fixed in the header before the runs: `polrad-full` R 1.0199 / 1.0261 / 1.0219 / 1.0098, `t-peak` 0.6118 / 0.8968 / 1.0112 / 1.0067 (a one-sided ceiling only); straight from the quadrature 1.0004 / 1.0010 / 0.9995 / 0.9866 (bin 4's −1.3 % real, cause unknown); the tail table adds +1.9 … +2.5 %; `test_bench_chain_rc.py` re-pointed to the history block; plan row 7, the README rows, `08_channel_priority.md` row F and REFERENCES updated; an adversarial reader confirmed the PASS as a measurement and found 9 items (`phase_B3_djangoh.md` §11): the tolerance constants now read back from the hashed block, the IEL32 cause corrected, overstatements removed, one stale test fixed | `phase_B3_djangoh.md` |
| C1 — R data; HERMES F₂ᵈ; COMPASS g₁ᵈ | DONE | uncommitted (phase E tree) | two new rows: `t1_r_world` (row 12) **PASS** — `r_sigma_lt` on NMC R 12.54/10 (p = 0.25, diagonal), `r1998` 6.20/10; NMC sys correlated 60.22/10 vs 21.28/10; `r_sigma_lt` fails Whitlow 209.76/68 and E140 29.66/9 where `r1998` (their fit) passes; H1 F_L consistent with both; the (1 + R) bias at the four sweet-spot cells 1.1848 / 1.1780 / 1.0762 / 0.9563 (the kernel's cos 2φ factor ×0.844 / 0.849 / 0.929 / 1.046); `t1_hermes_f2d` (row 14) **FAIL, recorded** — ToyF2 5303.5/68, MSTW 318.9/68, CT18NLO-LO 1788.2/45 on HERMES F₂ᵈ; the b₁ inversion reads per nucleon (0.929–1.039); the final COMPASS g₁ᵈ = HEPData ins1501480 Table 1 = the vendored file, 15/15 identical; `test_bench_t1_r_world.py` 51 | `phase_C1_r_world.md` |
| C2 — the Δ moment bound | DONE (two passes) | uncommitted (phase E tree) | the gluonic Soffer-type bound on Δ's n = 2 moment wired as `t5_delta_moment_bound` (row 13, **PASS**): the shipped toy uses 0.0644 / 0.0803 / 0.0941 of the room at Q² = 2 / 5 / 10 GeV² (0.129 / 0.161 / 0.188 under V6); helicity-state positivity min W 0.9494; the letter's `moment_A` / `moment_B` (dilution 1/3) 0.67–0.74 under A7 and 1.35–1.48 under V6; second pass: V4 shown a misreading (KS20 Eq. (68)), V6 = DS16 Eqs. (3)–(5) against the kernel's cos 2φ coefficient, a factor 2 (a registry candidate); `test_bench_t5_delta_bound.py` 21 | `phase_C2_delta_bound.md` (the board named `phase_C2_positivity.md` until phase E) |
| C3 — Abbott t₂₀ in the RC control | DEFERRED WITH DESIGN (review) | — | no generation path reads a deuteron elastic form factor; a deuteron `TabulatedSpin1FF` is a registry candidate | `DECISIONS_PRICED.md` |
| C4/C5 — Wiringa–Schiavilla; Boeglin | DONE | uncommitted (phase E tree) | two digitized vector figures: `t3_li6_ff_wiringa_schiavilla` (row 15) **PASS, one-sided** — the shipped C0 zero 3.0998 fm⁻¹ ≥ q_last − (δ + 0.1) = 2.8143 (WS98 draws IA+MEC C0 non-zero to 2.9190 and never draws the zero); row 5's FB zero would fail it; `t2_boeglin_nk` (row 16) **PASS** under the paper's 20 % FSI/off-shell allowance — Hulthén 21.98/21 (p = 0.40), AV18 11.54, CD-Bonn 8.89; all fail without it; only CD-Bonn passes with the MS FSI factor; `test_bench_li6_ff_theory.py` 18, `test_bench_boeglin.py` 13 | `phase_C4_wiringa_schiavilla.md`, `phase_C5_boeglin.md` |
| D — the least-verified inputs read | DONE | uncommitted (phase E tree) | D-state convention held (COMPASS PLB p. 39, HERMES Eq. (23), E143 Eq. (56), SMC Eq. (13)): verdict column stands, nothing swapped; E143 = hep-ph/9802357, SMC = PRD 58 112001 Table IX; every vendored T1 point equals its paper's table except 4 high-x COMPASS points (HEPData self-consistent, kept); NMC 2.5 % normalisation outside 'sys'; η_d confirmed; Q_d corrected 0.2859(3) → 0.2860(15) (recorded only: −5.42 / −5.71 / +1.52 %); α_s 0.1180(9) = PDG 2026; Suelzle ⁷Li a₀ = 1.755(20) F, ⁶Li no HO fit; Moniz k_F = 169 ± 5 MeV/c | `phase_D_inputs.md` |
| P — figures for the 18 existing rows; `make_report.py --plots` | DONE | uncommitted (phase E tree) | 11 figures (t1_g1d_world, t1_nmc_f2d, t1_nmc_f2d_over_f2p, t2_deeps_spectator_tail, t2_deuteron_static, t2_hermes_b1_table2, t3_li6_charge_ff_fb, t3_li6_rms_devries, t3_nmc_li6_over_d, t4_pythia_ep_closure, t5_sum_rules); none by design for EPIOS, EST, weighted/unweighted, T6; B1/B2/B3's and the new rows' drawn by their owners; `make_report.py --plots` (relative image links under each row's Detail heading; `--check` ignores every figure field; without `--plots` byte-identical to before); no `run()`, printed line or number moved (A/B 11/11); 19 of 23 harnesses define `plot()` after phase E | `phase_P_figures.md` |
| E — documents, report, priced decisions | DONE | uncommitted (phase E tree) | the report regenerated with figures (23 harnesses: **13 pass, 10 fail, 0 blocked**; 41 changes against 2026-09-27; `--check` OK); the shared files integrated from the ten phase records (below); the gates run once each (below) | this file, `SUMMARY.md`, `DECISIONS_PRICED.md`, `phase_E_report_compare.md` |
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

## Phase E — the integration stage (2026-10-26; times UTC, 2026-10-27 01:00–02:00)

**The report.** `LIPOLGEN_BENCH=1 python3 validation/benchmarks/make_report.py
--all --plots --compare docs/benchmarking/REPORT.json --out-dir
docs/benchmarking` (77 s; `--out-dir` spelled out because `make_report.py`
refuses `--plots` beside `--compare` without it — exit 2 by phase P's design
— and `--compare` alone writes nothing): **23 harnesses, 13 pass, 10 fail,
0 blocked, 0 blocked (environment), 0 error, 0 skipped (opt-in)**; 41 changes
(4 STATUS, 8 MOVED, 5 SUBROW, 24 SUBMOVED; 5 new rows and 52 new sub-rows
informational), every one traced to a phase; 19 of 23 rows carry a figure,
none failed to draw. Output verbatim: `phase_E_report_compare.md`. Then
`make_report.py --check`: OK (64 s).

**Shared files integrated** (each from a phase record's notes, the numbers as
the record states them; where two phases measured one quantity in different
windows both are given and the window named):

- `docs/benchmarking/BENCHMARK_PLAN.md` — §3: a dated status note in the M1
  row's "when" cell (the M2 pass rule untouched; the proposed clause is in
  `DECISIONS_PRICED.md`); §4: a dated 2026-10-26 status line, rows 6 and 10
  re-written from B2 and B1 (title, buys and status cells — both re-designed
  rows; the 2026-09-23 titles kept as "planned …"), row 7's test count
  (17 → 19 since the B3 verification), rows **12 `t1_r_world`, 13
  `t5_delta_moment_bound`, 14 `t1_hermes_f2d`, 15
  `t3_li6_ff_wiringa_schiavilla`, 16 `t2_boeglin_nk`** added — each names
  its own harness and no other, so `make_report.py`'s `plan_row_map` maps
  them one to one (checked; the report's plan-row column reads 12–16).
- `python/tests/test_bench_report.py` — `test_the_real_plan_maps_the_section_4_rows`
  now expects rows 1–16 and the five new harnesses, and that no row names
  two; **new** `test_import_pythia8_works_from_the_repository_root` (a fresh
  interpreter with cwd = the root, where a plain `import pythia8` fails with
  "ImportError: initialization failed" on this machine: the helper returns
  the module, `Pythia` is there, the cwd is unchanged, a second call returns
  the same object; then again in process; skips, naming the tier, without
  PYTHIA 8).
- `python/tests/test_bench_t4_t6.py` — `T4_PINNED` re-taken from this
  machine's run on PYTHIA 8.317 (LIPOLGEN_BENCH=1, 14 s, PASS, worst 0.076)
  with a dated comment keeping 8.312's values; tolerance (rel 3 %) and the
  verdict logic untouched. **The premise "the pin fails on 8.317" did not
  hold:** the 8.312 means sit within rel 3 % of 8.317's (largest move 1.97 %,
  W2 n_ch stock) and the test passed before the re-pin; the worst-distance pin
  was re-taken in phase 0 (`82b3237`). Re-taken anyway so that the pin is the
  pinned version's.
- `python/tests/test_bench_spin_deuteron.py` — `test_bonus_row_is_blocked_by_name`
  retired (B1; it asserted BLOCKED, row 10 is a FAIL measurement), the three
  history guards kept, the docstring's row-10 line re-pointed.
- `python/tests/test_bench_nuclear.py` — `test_rfy_tree_side_closed_form` and
  `test_rfy_magnetization_benchmark` retired (B2; superseded by
  `test_bench_rfy.py::test_tree_side_and_rule3` and `::test_verdict_is_pinned`),
  the docstring re-pointed. `test_bench_chain_rc.py` was re-pointed in the B3
  commit already (4 tests on the history block).
- `validation/benchmarks/README.md` — the harness count (23; 19 draw a
  figure); convention 4's `--plots` line and **convention 7** (`plot()`, phase
  P's note); the tools paragraph (`_plots.py`, `../price_defaults.py`);
  rows 4, 6, 7 (test count), 10 re-written; five new rows in the by-tier
  table (12, 14 under T1; 16 under T2; 15 under T3; 13 under T5); the LHAPDF
  sub-rows now measured on the T1 rows and T5 (dated); phase D's readings on
  the T1 and T2 rows; T4 on 8.317 (0.076) and its re-pinned test; the T6 test
  count (8 → 9; the ninth test predates this run); the PYTHIA paragraph; this run's previewed numbers and design blocks named in the tolerance paragraph.
- `README.md` — the Benchmarks paragraph: "As of 2026-10-26 (23 harnesses,
  every row run, the three LHAPDF sets installed): 13 pass, 10 fail
  (recorded), 0 blocked" (it read "As of 2026-09-23: 3 pass, 2 fail, 3
  blocked", and its "no test reads the variable" had been untrue since
  2026-09-26); the `--cluster-wave` row's −2.027 % sentence re-worded (phase A:
  it is the 2026-09-04 embedded-deuteron fix, not the hulthen → vmc price,
  which is ×0.430 / ×0.209 on A_zz^tag).
- `docs/references/REFERENCES.md` — items 1, 6, 7, 8, 13, 30, 59 **OBTAINED
  2026-10-26** (kept outside the repository, sha256 named; item 9 not
  obtained, no longer needed for row 6), with what each phase read from them;
  item 4 (the factor-2 Δ question), item 22; the R1998 entry ("LiPolGen's
  default R" was wrong: the default is `r_sigma_lt`) and E140 ("the dominant
  input" is not what the papers show); the vendored HEPData records named on
  the NMC, H1, Whitlow, E140 and HERMES entries; Boeglin, Wiringa–Schiavilla,
  DS16, CvDM, KS20; three new entries (SMC 1998, Bishop–Cheung 1979, PDG 2026)
  with a new status token `obtained 2026-10-26` in the vocabulary; a dated
  bookkeeping note under §10.1 (the greps now give 7 OBTAINED, 50 PLEASE
  DOWNLOAD, 3 `obtained`; nothing renumbered).
- `docs/benchmarking/08_channel_priority.md` — the R price corrected in §2
  rank 1 and §4 A (the function-level −8.3 / −6.7 / −3.1 % → −15.6 / −15.1 /
  −7.1 / +4.6 % at the cells, −15.2 / −14.7 / −7.8 / +4.2 % over the
  generated super-bins); dated status additions in §3.1, §3.4, §3.6, §3.7, §4
  rows A, C, F, G, H, I and §5 item 2 from phases A, B1, B2, C1, C2, C4, C5,
  D; one §7 log line. **No rank moved; no criterion score changed** (rank 1
  keeps C3 = 3 and C4 = 1; ranks 3, 5, 6 keep C3 = 3) — the log line says
  why.
- `docs/PHYSICS_CHANNELS.md` — one sentence the B3 record said went stale
  ("never against an external exact tail"): a dated italic sentence appended
  **on the same line** (no line shifted), giving row 7's external check; the
  strict gate re-run (below).
- Survey documents, dated corrections the workers asked for:
  `02_data_nucleon_deuteron.md` (§5.2 η_d's method — sub-Coulomb
  ²⁰⁸Pb(d,p), ¹³⁶Xe(d,p), not d + ⁴He; R-2 "the largest single" input; R-3,
  F-7 and §5.4 statuses), `04_theory.md` (T-20, §10 item 1),
  `00_in_tree_checks.md` (C-3), `06_critic.md` (§4.6: ADNDT 14's lithium
  "magnetization radii" are √(5/2)·a₀).
- Run records corrected minimally, each marked in place:
  `phase_0_environment.md` (phase A: the NNPDFpol11 Bjorken integral is
  **not** converged — the last decade adds 2.25 / 2.20 / 2.16 % at Q² =
  10 / 5 / 2 against ≤ 1 %, the first wording said "converged (last decade
  ≤ 0.4 %)"; and phase A/D's g₁ᵈ readings appended); `phase_P_figures.md`
  (the `t2_deuteron_static` row's Q_d reference and distances after phase D's
  correction).
- `tiers.py` **not edited**: `tiers.py --since master` prints no UNMAPPED
  file (every new harness matches `t[1-6]_*.py` and every new data file is
  named by its harness's constants).

**The gates, one at a time** (this machine, `source env.sh`; the full pytest
and the bench subset run sequentially, nothing else running):

| gate | result | when (UTC, 2026-10-27) |
|---|---|---|
| `./build/lipolgen_tests` (C++ unchanged) | **416 / 416 cases passed, 0 failed, 0 skipped; 17 242 712 / 17 242 712 assertions** | 01:01:14–01:08:18 (424 s) |
| `python3 -m pytest python/tests -q` | **1453 passed, 155 skipped, 0 failed**, 3 warnings (matplotlib's Axes3D import; two `IntegrationWarning`s from `t5_sum_rules`'s quad, as in phase 0); phase 0 had 1269 / 158 / 0 | 01:25:10–01:32:58 (467 s) |
| `LIPOLGEN_BENCH=1 python3 -m pytest python/tests -q -k bench` | **348 passed, 1 skipped, 0 failed** (1259 deselected; the skip is the nnpdf-wheel reproduction, which needs `$LIPOLGEN_NNPDF_WHEEL`), 3 warnings as above; the opt-in tests ran: T4's full row (re-pinned), the committed report's `--check`, the DJANGOH table-vs-quadrature diagnostic | 01:32:58–01:38:03 (304 s) |
| `python3 validation/check_physics_channels_links.py` (strict) | **exit 0** — `PHYSICS_CHANNELS.md` 1240 references checked (strict), 95 ranges, 7 external, **0 broken**, 24 allow-listed; `SPIN32_FINITE_GAMMA.md` 19 / 115 / 0 / **0 broken** / 10; `PYTHIA_BRIDGE.md` 0 / 0 / 8 / **0 broken** / 0 (run before and after the one same-line `PHYSICS_CHANNELS.md` edit) | 01:38:29 |
| `python3 validation/check_physics_channels_links.py --records-only` | **exit 0** — 830 citations in 67 dated run records: 335 checked against a live file, 57 into another record, 110 annotated historical, 129 read by hand, 199 skipped, **0 broken** | 01:38:30 |
| `python3 validation/check_physics_channels_links.py --audit-ranges` | **exit 0** — 209 range citations checked against their own citing sentences on 179 blocks, **0 REFUSED** | 01:38:32 |
| `python3 validation/check_spdx_headers.py` | **exit 0** — 133 files checked, all carry the SPDX line | 01:38:34 |
| `python3 validation/benchmarks/make_report.py --check` | **OK** — the committed-to-be `REPORT.md` / `REPORT.json` are what regeneration gives | 01:12:21–01:13:25 (64 s) |

(The task's `python -m pytest` is `python3 -m pytest` here: this machine has
no `python`.)

## Files, by phase (the working tree the supervising session commits)

Paths from the repository root; "new" = untracked before phase E. Every
`docs/benchmarking/figures/<harness>.png` was rewritten once more by phase E's
`--plots` run (19 files; `price_defaults.png` is not a harness's and was not).

| phase | new | modified |
|---|---|---|
| A | `validation/price_defaults.py`, `python/tests/test_price_defaults.py`, `docs/benchmarking/figures/price_defaults.png`, `docs/open_items/run_2026-10-26/phase_A_baseline.md` | — |
| B1 | `validation/benchmarks/data/bonus2014_tkachenko_supplemental.json`, `python/tests/test_bench_bonus.py`, `docs/benchmarking/figures/t2_bonus_spectator_shape.png`, `docs/open_items/run_2026-10-26/phase_B1_bonus.md` | `validation/benchmarks/t2_bonus_spectator_shape.py` |
| B2 | `validation/benchmarks/data/rfy1966_li_magnetic_ff.csv`, `validation/benchmarks/data/rfy1966_li_magnetization.csv`, `python/tests/test_bench_rfy.py`, `docs/benchmarking/figures/t3_li_magnetization_rfy.png`, `docs/open_items/run_2026-10-26/phase_B2_rfy.md` | `validation/benchmarks/t3_li_magnetization_rfy.py` |
| B3 | committed in `372b915`, `2d89a64` | `docs/benchmarking/figures/t4_djangoh_rad_noRad.png` (redrawn by phase E's `--plots`) |
| C1 | `validation/benchmarks/t1_r_world.py`, `validation/benchmarks/t1_hermes_f2d.py`, `validation/benchmarks/data/hepdata_ins424154_nmc_r.json`, `…/hepdata_ins786161_h1_fl.json`, `…/hepdata_ins296980_whitlow_r.json`, `…/hepdata_ins360765_e140_r.json`, `…/hepdata_ins894309_hermes_f2.json`, `python/tests/test_bench_t1_r_world.py`, `docs/benchmarking/figures/t1_r_world.png`, `…/t1_hermes_f2d.png`, `docs/open_items/run_2026-10-26/phase_C1_r_world.md` | — |
| C2 | `validation/benchmarks/t5_delta_moment_bound.py`, `validation/benchmarks/data/delta_moment_bound_sources.json`, `python/tests/test_bench_t5_delta_bound.py`, `docs/benchmarking/figures/t5_delta_moment_bound.png`, `docs/open_items/run_2026-10-26/phase_C2_delta_bound.md` | — |
| C4/C5 | `validation/benchmarks/t3_li6_ff_wiringa_schiavilla.py`, `validation/benchmarks/data/ws1998_li6_formfactors.csv`, `python/tests/test_bench_li6_ff_theory.py`, `validation/benchmarks/t2_boeglin_nk.py`, `validation/benchmarks/data/boeglin2011_deuteron_pmiss.csv`, `python/tests/test_bench_boeglin.py`, `docs/benchmarking/figures/t3_li6_ff_wiringa_schiavilla.png`, `…/t2_boeglin_nk.png`, `docs/open_items/run_2026-10-26/phase_C4_wiringa_schiavilla.md`, `…/phase_C5_boeglin.md` | — |
| D | `docs/open_items/run_2026-10-26/phase_D_inputs.md` | `validation/benchmarks/t1_g1d_world.py`, `t1_nmc_f2d.py`, `t1_nmc_f2d_over_f2p.py`, `t2_deuteron_static.py`, `t5_sum_rules.py` (headers; `t2_deuteron_static.py`'s digests), `validation/benchmarks/data/deuteron_static_measured.json` (the dated Q_d correction), `python/tests/test_bench_t1_nucleon.py`, `test_bench_static.py`, `test_bench_t5_closures.py` |
| P | `docs/open_items/run_2026-10-26/phase_P_figures.md`, eleven `docs/benchmarking/figures/*.png` | `validation/benchmarks/make_report.py` (`--plots`), the `plot()` of eleven harnesses (`t1_g1d_world`, `t1_nmc_f2d`, `t1_nmc_f2d_over_f2p`, `t2_deeps_spectator_tail`, `t2_deuteron_static`, `t2_hermes_b1_table2`, `t3_li6_charge_ff_fb`, `t3_li6_rms_devries`, `t3_nmc_li6_over_d`, `t4_pythia_ep_closure`, `t5_sum_rules`), their tests in `python/tests/test_bench_{t1_nucleon,static,t5_closures,spin_deuteron,deeps,nuclear,t4_t6,report}.py` |
| E | `docs/open_items/run_2026-10-26/SUMMARY.md`, `DECISIONS_PRICED.md`, `phase_E_report_compare.md` | `docs/benchmarking/REPORT.md`, `REPORT.json`, `BENCHMARK_PLAN.md`, `08_channel_priority.md`, `00_in_tree_checks.md`, `02_data_nucleon_deuteron.md`, `04_theory.md`, `06_critic.md`; `validation/benchmarks/README.md`; `README.md`; `docs/references/REFERENCES.md`; `docs/PHYSICS_CHANNELS.md` (one line); `python/tests/test_bench_report.py`, `test_bench_t4_t6.py`, `test_bench_spin_deuteron.py`, `test_bench_nuclear.py`; this file; `phase_0_environment.md` and `phase_P_figures.md` (marked corrections) |

## What the run did not apply, and where it went

| item | from | why not applied | where it went |
|---|---|---|---|
| `BENCHMARK_PLAN.md` §6 U-10, a dated note that F_m's shape now has data over q = 0.85–1.39 fm⁻¹ (row 6) and the second lobe none | B2 | phase E's scope in the plan is §3's status note and §4; §6 is also copied verbatim into `REPORT.json`, so an edit needs a regeneration | the exact text below (E-1); the author or the next run |
| `BENCHMARK_PLAN.md` §7 coverage: `fsi.*` → T2 (rows 10, 11 evaluate `GlauberFsiWeight`), with the matching `tiers.py` SECTION7 entry (the two change together, `test_section7_as_data_matches_the_plan`) | B1 | §7 outside phase E's plan scope; `tiers.py` is edited only for an unmapped file and none was | `phase_B1_bonus.md` §10 (the gap); the exact text below (E-2) |
| `BENCHMARK_PLAN.md` §7: T5 for `sf.hpp` / `sf.cpp` / `xsec.*` / `pipeline.*` / a PDF-set upgrade (the Δ bound reads them) | C2 | as above | the exact text below (E-3) |
| `BENCHMARK_PLAN.md` §7: T3 for `rc.*` (`HoSpin1FF` lives in `rc.hpp`; rows 5, 6, 15 and `t3_li6_rms_devries` read it) | C4 | as above | below (E-4; C4 asked it as a question, no text proposed) |
| `include/lipolgen/rc.hpp`: `RcTailModel::PolradFull` point (3) and the first bullet of "FIVE THINGS THIS IS NOT" are partly stale since row 7 | B3 | no C++ file changes in this run (the C++ suite was run on the unchanged build); the PHYSICS_CHANNELS gate cites lines of `rc.hpp` | `phase_B3_djangoh.md` §10; the author |
| `include/lipolgen/b1_nuclear.hpp`, the "THE 100 % BAND IS MANDATORY" comment: Q_d +0.2859(3) fm² is the AV18 / CD-Bonn attribution, Bishop–Cheung conclude 0.2860 ± 0.0015 | D | no C++ change | `phase_D_inputs.md` §8 names it; the exact text below (E-5) |
| `tests/test_cluster.cpp`, the CD-Bonn asymptotics comment (0.0256(4) is Rodning–Knutson's measurement, sub-Coulomb ²⁰⁸Pb(d,p), ¹³⁶Xe(d,p)) | D | no C++ change | the exact text below (E-6) |
| `tests/test_rc.cpp`'s "NO 6Li table is shipped" comment and `rc.cpp`'s error message (WS98 now digitized, under `validation/benchmarks/data/`) | C4 | no C++ change; the test still holds | C4 asked "reword?", no text proposed |
| `docs/PHYSICS_CHANNELS.md` §2 Δ row: "twice what DS16 Eqs. (3)–(5) give per unit Δ" | C2 | offered "only if judged stale; nothing in it is false" — not stale | the exact text below (E-7); `DECISIONS_PRICED.md` §3 item 4 |
| the vendored nnpdf T1 files' provenance strings (E143 "CONFLICT", SMC "NOT IDENTIFIED", d_state "UNVERIFIED") | D | a reference-file change (vendor script, regeneration, new sha256 pins); the harness header carries the resolution | `phase_D_inputs.md` §1 (the readings); the author |
| SMC into `t1_g1d_world`'s verdict | D | would change a declared rule after its χ² (toy 16.7/12) was seen | `SUMMARY.md` item 6 |
| NMC R's correlated-systematic reading as row 12's verdict (60.22/10 would flip it to FAIL) | C1 | the declared tolerance is the diagonal one; no tolerance moves | `SUMMARY.md` item 6 |
| a normalisation-profiled χ² for MSTW on NMC F₂ᵈ (1.0 normalisation σ) | D | not computed (would need a tolerance declared first) | `phase_D_inputs.md` |
| KS20's DOI in `delta_moment_bound_sources.json` | C2 | changes the vendored file's digests | C2's return: REFERENCES.md §6d gives 10.1103/PhysRevD.101.054011 |
| `tiers.py` printing "regenerate with `--plots`" when `_plots.py` changes | P | `tiers.py` edited only for an unmapped file | P's return (an open question); this table |
| `t2_hermes_b1_table2`'s figure in b₁ᵈ rather than x·b₁ᵈ; error bars on the Deeps figure (needs a `run()` change) | P | presentation choices; a `run()` change is a later task's | `phase_P_figures.md` §2 (the row table) |
| row 10's tolerance re-declared (a shape χ² of the weight, or stat ⊕ smooth syst) | B1 | the declared rule is recorded, not re-written after its number | `SUMMARY.md` item 6 |
| HEPData ins894309's Table 2 defect; COMPASS PLB 769 Table A.1 vs HEPData; PYTHIA 8.317's loader | C1, D, phase 0 | nothing is sent from this run | `SUMMARY.md` item 5 |
| the fetched documents filed into `PolarizedLithiumSim/refs/` | B1, B2, D, C2 | the sibling repository is read-only | `DECISIONS_PRICED.md` §3 item 7; `SUMMARY.md` item 3 |
| C4's cross-check sentence for row 6 (WS98 drew RFY Table I at 0.7116 × RFY's normalisation; row 6's free scale is that factor) | C4 | applied to `BENCHMARK_PLAN.md` row 6 and the README row instead of B2's record (a worker's record is not re-written for a cross-check) | plan §4 row 6, README row 6 |
| the ePIC chain legs; Abbott's t₂₀ / A(Q²) | PLAN | not in this run (no container stack; no generation path reads a deuteron elastic form factor) | `PLAN.md` "Not in this run"; `DECISIONS_PRICED.md` §3 item 2 |

### The exact texts not applied (from the phases' returns, which are not in the repository)

- **E-1** (B2) — `BENCHMARK_PLAN.md` §6 row U-10, 'distance' cell, append:
  "; *2026-10-26:* F_m's SHAPE now has data over q = 0.85–1.39 fm⁻¹ (RFY Table
  I; row 6 PASS at p = 0.0378 with one free scale, 0.0145 at s = 1); the second
  lobe near 2 fm⁻¹ still has none".
- **E-2** (B1) — `BENCHMARK_PLAN.md` §7, after the bullet "`breakup.*`,
  `spectator.*`, `fsi.*`, `triton_sf.*` → T3-cluster, T4-tagged", the bullet
  "- `fsi.*` → T2 (rows 10 and 11 evaluate `GlauberFsiWeight` on the deuteron
  control; added 2026-10-26)"; and in `validation/benchmarks/tiers.py`
  SECTION7 at the same position the entry `(('fsi.*',), (), ('T2',), 'rows 10
  and 11 evaluate GlauberFsiWeight on the deuteron control (2026-10-26)'),`.
- **E-3** (C2) — `BENCHMARK_PLAN.md` §7, a proposal: "`t5_delta_moment_bound`
  reads `toy_delta_gluon` and `make_moment_a` / `make_moment_b` (sf.hpp /
  sf.cpp), the kernel's cos 2φ term (xsec.cpp) and the pipeline's Δ slot
  (pipeline.cpp), and the CT18NLO / NNPDFpol11_100 sets; the sf.hpp bullet and
  the 'any PDF-set version' note re-run T1, T2, T3-EMC only and xsec.* /
  pipeline.* / sf.cpp map to no tier, so add T5 to the sf.hpp bullet and the
  PDF-set note (or map xsec.* and pipeline.* to T5)."
- **E-4** (C4, a question) — §7 maps `rc.*` to T4-RC and T2 only, so a change
  to `HoSpin1FF` in `rc.hpp` / `rc.cpp` re-runs none of the T3 form-factor rows
  (5, 6, 15 and `t3_li6_rms_devries`): add T3 to that entry?
- **E-5** (D) — `include/lipolgen/b1_nuclear.hpp`, the "THE 100 % BAND IS
  MANDATORY" block (comment only, not fingerprinted by
  `physics_channels_ranges.json`): replace "Q_d: Bishop and Cheung, Phys. Rev.
  A 20 (1979) 381.)" with "Q_d: +0.2859(3) fm^2 as AV18 (PRC 51 (1995) 38,
  Table X) and CD-Bonn (PRC 63 (2001) 024001, Table XV) quote it, citing Bishop
  and Cheung, Phys. Rev. A 20 (1979) 381, whose own conclusion is 0.2860 +-
  0.0015 fm^2 (read 2026-10-26) -- a 0.035 % difference that changes nothing
  here.)". Leave `include/lipolgen/pipeline.hpp`'s comment alone (its line is
  fingerprinted).
- **E-6** (D) — `tests/test_cluster.cpp`, the CD-Bonn asymptotics comment:
  replace "Table XV prints A_S = 0.8846(9) fm^-1/2 and eta = 0.0256(4) for
  CD-Bonn (the parentheses are the EMPIRICAL errors it is compared against)."
  with "Table XV prints A_S = 0.8846 fm^-1/2 and eta = 0.0256 for CD-Bonn,
  beside the EMPIRICAL 0.8846(9) fm^-1/2 (its refs. [29, 73]) and 0.0256(4) --
  the latter Rodning & Knutson's MEASUREMENT, PRC 41 (1990) 898, Eq. (10)
  (tensor analysing powers of sub-Coulomb 208Pb(d,p) and 136Xe(d,p) analysed
  in DWBA; read 2026-10-26)."
- **E-7** (C2, optional) — `docs/PHYSICS_CHANNELS.md` §2 Δ row, after
  "Standalone identity: a_cos2phi(Delta,F1,F2,x,y) = -(1-y)/y^2 Delta/D_phi at
  lambda_m = +-1, theta_m = 90 deg" add "(twice what DS16 Eqs. (3)-(5) give
  per unit Delta: run 2026-10-26 phase C2 section 4.2, recorded, not changed)".

## Verification

Filled in by the fix stage.
