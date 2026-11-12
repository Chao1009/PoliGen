<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-26 — status board

The plan is `PLAN.md`; the one page is `SUMMARY.md`; the ranking this run
was asked for is `../../benchmarking/08_channel_priority.md`. The decision
registry is still `../run_2026-09-03/STATUS.md` (27 rows); cite a decision
by its row number there; this run's priced options are in
`DECISIONS_PRICED.md`. **No physics default, registry status label,
constant or tolerance moves in this run, and no `validation/reference/*.json`
number; one vendored benchmark reference moves, the dated transcription
correction PLAN.md D2 allows (Q_d in
`validation/benchmarks/data/deuteron_static_measured.json`, 0.2859(3) →
0.2860(15) fm², recorded only, no verdict moves); no decision is taken.**
Work is on branch `claude/run-2026-10-26`, committed and pushed: phases 0,
0r and B3 one by one, A–E together as `6b1b63e` (worded so in phase F; the
board said "phase by phase" and "uncommitted (phase E tree)" until then);
phase F's fixes, with the §7 map extension, as `1c5cf68` (pushed).

## Phase board

| phase | state | commit | what it did, with its window | record |
|---|---|---|---|---|
| 0 — environment, ranking, plan | DONE | `82b3237` | the dependency stack built with the project's script (HepMC3 3.3.0, LHAPDF 6.5.5 + the three sets, PYTHIA 8.317; 439 MB, 46 min on 4 cores); LiPolGen built with every tier on; **C++ 416/416 cases, 17 242 712 assertions, 0 failed; pytest 1269 passed / 158 skipped / 0 failed**; the first report with the LHAPDF sets (5 pass / 8 fail / 5 blocked before the PYTHIA fix); PYTHIA 8.317's cwd-relative python loader found and worked around in `_report.py` (`import_pythia8`), T4 PASS at worst 0.076 on 8.317, T6 23/23 | `phase_0_environment.md`, `../../benchmarking/08_channel_priority.md`, `PLAN.md` |
| 0r — the ranking and plan reviewed by two adversarial Opus readers | DONE | `87787dd` | physics lens 24 findings (8 major), honesty/feasibility lens 25 (5 major); 49 in all, 45 applied in the rewrite of `08_channel_priority.md` and `PLAN.md`, 4 applied in part (below); the DJANGOH worker, started on the pre-review brief, was stopped before it built anything and relaunched on the reviewed one | this file (below) |
| A — baseline and price table | DONE | `6b1b63e` | tallies C++ 416/416 (17 242 712 assertions), pytest 1269 passed / 158 skipped / 0 failed, against 2026-09-27's 13 and 39 LHAPDF failures; the first `--compare` with the LHAPDF sets, 38 changes, one line each; T4/T6 re-run after the PYTHIA workaround (PASS, 0.076; 23/23), the 7 / 8 / 3 tally inferred for the 18 harnesses (with B1–B3 the same 18 read 9 / 9 / 0, measured in phase E; 13 / 10 / 0 with the five new rows); the default-input price table of the ranked observables at config 1, N = 200 000, seed 20261026 (`validation/price_defaults.py`, `price_defaults.png`): R −15.2 / −14.7 / −7.8 / +4.2 % on the cos 2φ amplitude over the four super-bins x = 0.011 / 0.028 / 0.071 / 0.141, F₂ ≤ 0.23 %; A_zz ×1.2524 / ×1.2604; A_∥ ×0.532 / ×0.514 under NNPDFpol1.1; A_zz^tag ×0.430 / ×0.209 under VMC; c₂ ×0.270 / ×0.158 at the measured quadrupole; `test_price_defaults.py` (1) | `phase_A_baseline.md` |
| B1 — row 10 BONuS | DONE | `6b1b63e` | the BONuS Supplemental Material vendored as numbers (100 R_D/S blocks + the DIS F₂ⁿ/F₂ᵖ table, the original's sha256 recorded); row 10 BLOCKED → **FAIL, recorded** under a rule declared before any number: the Hulthén 40 mb `GlauberFsiWeight` over its backward-window mean inside the data's 1σ flat band in 4/12 (4 GeV) and 0/16 (5 GeV) DIS blocks; the backward data are not flat (χ² 739.4/84) and fall toward 90° as the weight does (χ² 226.3/84 against its shape; 20 mb 403.3); `test_bench_bonus.py` 26 passed | `phase_B1_bonus.md` |
| B2 — row 6 RFY | DONE | `6b1b63e` | row 6 BLOCKED → **PASS**: RFY Table I's seven ⁶Li F_M² points against `HoSpin1FF`'s unfitted F_mag², χ² = 13.35/6 with one free scale s = 1.419 ± 0.206 (p = 0.0378, against p ≥ 0.01 declared in a hashed design block); s = 1 recorded 17.49/7 (p = 0.0145); ⁷Li not compared; the radius rule retired (ADNDT 14's lithium 'magnetization radii' are √(5/2)·a₀); `test_bench_rfy.py` 16 | `phase_B2_rfy.md` |
| B3 — row 7 DJANGOH with the elastic tail on | DONE, verified | `372b915`, `2d89a64` | DJANGOH 4.6.21 built from `e356b84` (gfortran 11.4, LHAPDF 6.5.5 + cteq61 = 10150); the eight ePIC numbers reproduced (\|pull\| ≤ 1.74); tail isolated (IEL31+32+33, IEL2 = 0, WMIN 0 → M_p, x_l ≤ 0.1096 / 0.5518 / 0.9249 / 0.95, y ≤ 0.985, Q² bins by subtraction since the IEL channels ignore Q2MAX); row 7 BLOCKED → PASS under the tolerance fixed in the header before the runs: `polrad-full` R 1.0199 / 1.0261 / 1.0219 / 1.0098, `t-peak` 0.6118 / 0.8968 / 1.0112 / 1.0067 (a one-sided ceiling only); straight from the quadrature 1.0004 / 1.0010 / 0.9995 / 0.9866 (bin 4's −1.3 % real, cause unknown); the tail table adds +1.9 … +2.5 %; `test_bench_chain_rc.py` re-pointed to the history block; plan row 7, the README rows, `08_channel_priority.md` row F and REFERENCES updated; an adversarial reader confirmed the PASS as a measurement and found 9 items (`phase_B3_djangoh.md` §11): the tolerance constants now read back from the hashed block, the IEL32 cause corrected, overstatements removed, one stale test fixed | `phase_B3_djangoh.md` |
| C1 — R data; HERMES F₂ᵈ; COMPASS g₁ᵈ | DONE | `6b1b63e` | two new rows: `t1_r_world` (row 12) **PASS** — `r_sigma_lt` on NMC R 12.54/10 (p = 0.25, diagonal), `r1998` 6.20/10; NMC sys correlated 60.22/10 vs 21.28/10; `r_sigma_lt` fails Whitlow 209.76/68 and E140 29.66/9 where `r1998` (their fit) passes; H1 F_L consistent with both; the (1 + R) bias at the four sweet-spot cells 1.1848 / 1.1780 / 1.0762 / 0.9563 (the kernel's cos 2φ factor ×0.844 / 0.849 / 0.929 / 1.046); `t1_hermes_f2d` (row 14) **FAIL, recorded** — ToyF2 5303.5/68, MSTW 318.9/68, CT18NLO-LO 1788.2/45 on HERMES F₂ᵈ; the b₁ inversion reads per nucleon (0.929–1.039); the final COMPASS g₁ᵈ = HEPData ins1501480 Table 1 = the vendored file, 15/15 identical; `test_bench_t1_r_world.py` 51 | `phase_C1_r_world.md` |
| C2 — the Δ moment bound | DONE (two passes) | `6b1b63e` | the gluonic Soffer-type bound on Δ's n = 2 moment wired as `t5_delta_moment_bound` (row 13, **PASS**): the shipped toy uses 0.0644 / 0.0803 / 0.0941 of the room at Q² = 2 / 5 / 10 GeV² (0.129 / 0.161 / 0.188 under V6); helicity-state positivity min W 0.9494; the letter's `moment_A` / `moment_B` (dilution 1/3) 0.67–0.74 under A7 and 1.35–1.48 under V6; second pass: V4 shown a misreading (KS20 Eq. (68)), V6 = DS16 Eqs. (3)–(5) against the kernel's cos 2φ coefficient, a factor 2 (a registry candidate); `test_bench_t5_delta_bound.py` 21 | `phase_C2_delta_bound.md` (the board named `phase_C2_positivity.md` until phase E) |
| C3 — Abbott t₂₀ in the RC control | DEFERRED WITH DESIGN (review) | — | no generation path reads a deuteron elastic form factor; a deuteron `TabulatedSpin1FF` is a registry candidate | `DECISIONS_PRICED.md` |
| C4/C5 — Wiringa–Schiavilla; Boeglin | DONE | `6b1b63e` | two digitized vector figures: `t3_li6_ff_wiringa_schiavilla` (row 15) **PASS, one-sided** — the shipped C0 zero 3.0998 fm⁻¹ ≥ q_last − (δ + 0.1) = 2.8143 (WS98 draws IA+MEC C0 non-zero to 2.9190 and never draws the zero); row 5's FB zero would fail it; `t2_boeglin_nk` (row 16) **PASS** under A = 0.20, the upper end of the paper's 10–20 % FSI estimate (FSI and the CC1 off-shell prescription together; the Hulthén fails at 0.10, phase F) — Hulthén 21.98/21 (p = 0.40), AV18 11.54, CD-Bonn 8.89; all fail without it; only CD-Bonn passes with the MS FSI factor; `test_bench_li6_ff_theory.py` 18, `test_bench_boeglin.py` 13 | `phase_C4_wiringa_schiavilla.md`, `phase_C5_boeglin.md` |
| D — the least-verified inputs read | DONE | `6b1b63e` | D-state convention held (COMPASS PLB p. 39, HERMES Eq. (23), E143 Eq. (56), SMC Eq. (13)): verdict column stands, nothing swapped; E143 = hep-ph/9802357, SMC = PRD 58 112001 Table IX; every vendored T1 point equals its paper's table except 4 high-x COMPASS points (HEPData self-consistent, kept); NMC 2.5 % normalisation outside 'sys'; η_d confirmed; Q_d corrected 0.2859(3) → 0.2860(15) (recorded only: −5.42 / −5.71 / +1.52 %); α_s 0.1180(9) = PDG 2026; Suelzle ⁷Li a₀ = 1.755(20) F, ⁶Li no HO fit; Moniz k_F = 169 ± 5 MeV/c | `phase_D_inputs.md` |
| P — figures for the 18 existing rows; `make_report.py --plots` | DONE | `6b1b63e` | 11 figures (t1_g1d_world, t1_nmc_f2d, t1_nmc_f2d_over_f2p, t2_deeps_spectator_tail, t2_deuteron_static, t2_hermes_b1_table2, t3_li6_charge_ff_fb, t3_li6_rms_devries, t3_nmc_li6_over_d, t4_pythia_ep_closure, t5_sum_rules); none by design for EPIOS, EST, weighted/unweighted, T6; B1/B2/B3's and the new rows' drawn by their owners; `make_report.py --plots` (relative image links under each row's Detail heading; `--check` ignores every figure field; without `--plots` byte-identical to before); no `run()`, printed line or number moved (A/B 11/11); 19 of 23 harnesses define `plot()` after phase E | `phase_P_figures.md` |
| E — documents, report, priced decisions | DONE | `6b1b63e` | the report regenerated with figures (23 harnesses: **13 pass, 10 fail, 0 blocked**; 41 changes against 2026-09-27; `--check` OK); the shared files integrated from the ten phase records (below); the gates run once each (below) | this file, `SUMMARY.md`, `DECISIONS_PRICED.md`, `phase_E_report_compare.md` |
| F — verification and fix stage | DONE | `1c5cf68` | four adversarial lenses over `6b1b63e`: 38 findings (physics 8, one major; honesty 14; code 9; regression 7), each re-checked against the tree before anything was touched — none refuted, 33 applied, 4 in part, 1 not (outside this run's files); no default, constant, tolerance, registry label, verdict or reference number moved, no vendored data file edited; the report regenerated (23 harnesses, 13 / 10 / 0, **0 changes** against `6b1b63e`'s, `--check` OK); pytest 1456 / 155 / 0, the bench subset 351 / 1 / 0, the three citation gates and SPDX exit 0 | this file (Verification) |

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
| `python3 validation/check_physics_channels_links.py --records-only` | **exit 0** — 830 citations in 67 dated run records: 335 checked against a live file, 57 into another record, 110 annotated historical, 129 read by hand, 199 skipped, **0 broken** *(phase F: the gate reads only backticked citations; `phase_C1_r_world.md`'s 18 were plain text, and the ones into `08_channel_priority.md` and `REFERENCES.md`, which phase E edited, had drifted; all 18 re-read in phase F, the seven the gate checks unannotated backticked, the drifted ones re-pointed in plain text with the commit they held at — "Verification", H-14)* | 01:38:30 |
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
| F | — (no new file in the repository; scratch scripts only) | `validation/benchmarks/t5_delta_moment_bound.py` (a missing set raises, LHAPDF's search order), `python/tests/test_bench_t5_delta_bound.py` (21 → 23), `validation/benchmarks/make_report.py` (one legend sentence), `python/tests/test_bench_report.py` (+1 test), the headers of `t4_pythia_ep_closure.py`, `t5_sum_rules.py`, `t1_nmc_f2d.py`, `t2_deuteron_static.py` (dated notes, no number); `docs/benchmarking/REPORT.md`, `REPORT.json` and the 19 figures (regenerated: date, revision, runtimes, the legend sentence, the figure footers); `docs/benchmarking/08_channel_priority.md`, `BENCHMARK_PLAN.md`, `02_data_nucleon_deuteron.md`; `validation/benchmarks/README.md`; this file, `SUMMARY.md`, `DECISIONS_PRICED.md`, `PLAN.md`, `phase_0_environment.md`, `phase_A_baseline.md`, `phase_C1_r_world.md`, `phase_C4_wiringa_schiavilla.md`, `phase_C5_boeglin.md`, `phase_D_inputs.md` (marked corrections) |
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
| KS20's and CRSZ25's DOIs in `delta_moment_bound_sources.json` (both `null` for published articles, PRD 101 (2020) 054011 and EPJ A 61 (2025) 83; CRSZ25 added in phase F) | C2; phase F | changes the vendored file's digests | C2's return: REFERENCES.md §6d gives 10.1103/PhysRevD.101.054011; fill both at the file's next revision |
| `tiers.py` printing "regenerate with `--plots`" when `_plots.py` changes | P | `tiers.py` edited only for an unmapped file | P's return (an open question); this table |
| `t2_hermes_b1_table2`'s figure in b₁ᵈ rather than x·b₁ᵈ; error bars on the Deeps figure (needs a `run()` change) | P | presentation choices; a `run()` change is a later task's | `phase_P_figures.md` §2 (the row table) |
| row 10's tolerance re-declared (a shape χ² of the weight, or stat ⊕ smooth syst) | B1 | the declared rule is recorded, not re-written after its number | `SUMMARY.md` item 6 |
| HEPData ins894309's Table 2 defect; COMPASS PLB 769 Table A.1 vs HEPData; PYTHIA 8.317's loader | C1, D, phase 0 | nothing is sent from this run | `SUMMARY.md` item 5 |
| the fetched documents filed into `PolarizedLithiumSim/refs/` | B1, B2, D, C2 | the sibling repository is read-only | `DECISIONS_PRICED.md` §3 item 7; `SUMMARY.md` item 3 |
| C4's cross-check sentence for row 6 (WS98 drew RFY Table I at 0.7116 × RFY's normalisation; row 6's free scale is that factor) | C4 | applied to `BENCHMARK_PLAN.md` row 6 and the README row instead of B2's record (a worker's record is not re-written for a cross-check) | plan §4 row 6, README row 6 |
| the ePIC chain legs; Abbott's t₂₀ / A(Q²) | PLAN | not in this run (no container stack; no generation path reads a deuteron elastic form factor) | `PLAN.md` "Not in this run"; `DECISIONS_PRICED.md` §3 item 2 |
| a dated note in `deuteron_static_measured.json`'s `_provenance.corrections` that Q_d = 0.285699(15)(18) fm² (Puchalski–Komasa–Pachucki 2020) supersedes both values | F (P-7) | the fix stage edits no vendored reference file (and the note would re-pin FILE_SHA256); the note is in `t2_deuteron_static.py`'s header and `SUMMARY.md` | the file's next revision; re-vendoring is the author's call |
| the Licence lines of `ws1998_li6_formfactors.csv` and `boeglin2011_deuteron_pmiss.csv` ("no figure, text or table is reproduced", while their headers quote captions and sentences) | F (H-20) | re-pins DATA_SHA256 in two harnesses and their tests | noted in `phase_C4_wiringa_schiavilla.md` §8, `phase_C5_boeglin.md` §8, `DECISIONS_PRICED.md` §3 item 8; the files' next revision |
| E-2 and E-4 (§7 `fsi.*` → T2, `rc.*` → +T3) in `BENCHMARK_PLAN.md` §7 and `tiers.py` SECTION7 together | F (C-25), B1, C4 | the fix stage left them (`tiers.py` outside its files); **applied by the supervising session in the phase-F commit**, both files together, with `test_a_diff_maps_to_tiers_and_harnesses` extended (`rc.cpp` → T2, T3, T4-RC; `fsi.hpp` → T2, T3-cluster, T4-tagged) | done |
| matplotlib in CI's two pip lines (`.github/workflows/ci.yml`, lines 206 and 308) | F (R-37) | `ci.yml` is not a file this run touched | `validation/benchmarks/README.md` convention 7 now says every `plot()` test skips in CI |
| `sys.dont_write_bytecode` around the `polligen` fixture's import of the sibling packages | F (H-19, optional) | `python/tests/conftest.py` is not a file this run touched | `DECISIONS_PRICED.md` §3 item 7 names the gitignored `.pyc` it leaves |
| the full annotation of `phase_C1_r_world.md`'s citations (four `(anchor read …)` pins, seven `(as of 2d89a64)`) | F (R-32) | applied first, then reverted: it moves `python/tests/test_doc_link_gate.py`'s census (`pins == 62`, `checked == 110`, and the `ASOF_READ_BY_HAND` list — the `08_channel_priority.md:91` anchor shares only "sweet" with its line at `2d89a64`), a file this run did not touch | "Verification" below; whoever next updates that census |

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

**Phase F** (2026-10-26; times UTC, 2026-10-27 02:45–04:00). Four adversarial
lenses read `6b1b63e`'s tree — physics (8 findings, one major), honesty (14),
code (9), regression (7): **38 findings, 29 distinct issues once the
duplicates are merged**. The fix stage re-checked every one against the tree
first, re-running the verifiers' scripts (copied into the session scratch
directory, `fixG/`) or the harnesses wherever a number was involved, and only
then edited. **None was refuted; 33 were applied as proposed (or by the
alternative the finding itself offered), 4 in part, 1 not at all** (it needs
a file this run did not touch). No physics default, constant, tolerance,
registry label, verdict or reference number moved; no file under
`validation/benchmarks/data/` or `validation/reference/` was edited; no C++
file changed. Numbering: P-1 … P-8, H-9 … H-22, C-23 … C-31, R-32 … R-38, in
the order the findings were handed over.

| # | lens, severity | finding (abridged) | re-checked | what was done |
|---|---|---|---|---|
| P-1 | physics, **major** | "F₂ cancels exactly in the cos 2φ amplitude and R is the lever" (`08` §0 (a), §2 rank 1, §3.1, §4 A; SUMMARY, DECISIONS_PRICED, plan row 12, phases A, D) holds only for the shipped Δ = 10⁻² F₁ x^0.3 (1−x)^4; per unit Δ F₂ is the 1:1 lever and R ≤ 5 × 10⁻⁵ | **confirmed by running** (the verifier's `rprice.py`, `f2price.py`): the R swap moves a₂/Δ by +4.4 / +0.7 / +0.4 / −1.3 × 10⁻⁵ at (0.0112, 1.14), (0.0282, 1.14), (0.0708, 3.13), (0.141, 14.3), while the toy's Δ itself moves ×0.84405 / 0.84893 / 0.92916 / 1.04568 — the F₁ ratio, i.e. the quoted R price; at fixed Δ toy → CT18NLO moves a₂ ×1.3595 / 1.1878 / 0.9612 / 0.9029 = F₂(toy)/F₂(CT18NLO) (Q² = 1.7 at the two low-x cells, CT18NLO's floor), the toy-Δ kernel 1.000000; the sibling's `moment_A` moves ×0.829 / 0.831 / 0.829 / 0.829 under F₂ alone (`plans/07` WP1 addendum, read). Phase C1 had said it (its §8 item 1: "right only with Δ/F₁ held fixed") | applied: `08` §0 (a), §2 rank 1 (C4 1 → 2 — F₂(⁶Li)/F₂(D), row 4, is a same-A anchor of an input the extraction divides by — score 28 → 29, **no rank moved**), §3.1 twice, §4 A (moves, development item), a §7 log line; SUMMARY; DECISIONS_PRICED row 3 and the unpolarized-backend bullet; plan row 12 "buys"; marked corrections in `phase_A_baseline.md` (§4.3 and two sentences) and `phase_D_inputs.md`. Not edited, in files this run did not touch: `docs/DEVELOPMENT_PLAN.md`'s "F₂ cancels in a₂" and `docs/surveys/needs_survey.md`'s reading of `test_xsec_identity.py` (a closed form in Δ/F₁, right as written) |
| P-2 | physics, minor | CT18NLO's 8.6 % NMC F₂ᵈ excess "of the size of the NLO coefficient functions" was never tested; `LhapdfSF`'s LO F₂ carries massless c and b | **confirmed by running** (`charm.py`, CT18NLO member 0 through PYTHIA's LHAPDF6 plugin, the 146 points): c + b = 3.56 % of the model F₂ᵈ on average (7.0 % at x < 0.02, 7.6 % at 0.02–0.05, 12.7 % at most); without them the mean model/data 1.0862 → 1.0459, per x bin 1.064 / 1.084 / 1.084 / 1.072, 0.977 above x = 0.2; the five weights include c (4/9) and b (1/9) in `src/lhapdf/lhapdf_sf.cpp`. The verifier's "about 40 %" of the excess is 47 % (0.040 of 0.086): written "about half" | applied: SUMMARY, `08` §4 A status, a dated note in `t1_nmc_f2d.py`'s header (no number or tolerance moved), DECISIONS_PRICED's backend bullet |
| P-3 | physics, minor | A_zz^tag ×0.43 / ×0.21 attributed to "the AV18 deuteron"; `--cluster-wave vmc` on tagged ⁶Li α swaps the α–d overlap (P_D 0.0867 → 0.019355), which no p–n anchor tests | confirmed by reading `src/core/tagged.cpp` (the VmcAV18 branch, `li6_vmc_waves`, and `DEUTERON_AV18()` for the embedded deuteron) and `07_cw_sign_investigation.md` §6.4 (Hulthén at the VMC P_D: −0.6001 of the −1.2069 → −0.5191 move at k = 0.198 GeV) | applied: SUMMARY (twice); `08` §4 G renamed "the cluster relative waves" (the old name kept in the cell) with the U-2 status sentence; DECISIONS_PRICED row 9 |
| P-4 | physics, minor | row 13's headline "6–9 % of the room" is the first-pass identification A7, which the harness's own S2 shows inconsistent with DS16 | confirmed by reading the header's S2 and running its tests (`test_s2_ds16_projection_operators_give_eq5_and_the_f1_normalisation`, `test_s2_v6_factor_is_ds16_eq5_against_the_kernel_coefficient` pass); the verdict, computed under A7, holds under V6 | applied: `08` §4 C, SUMMARY and plan row 13 headline V6's 13–19 % (0.1289 / 0.1605 / 0.1881) and name A7's 6–9 % as the first-pass reading; no assumption, verdict or tolerance changed |
| P-5 | physics, minor | phase C5 §4.1: "the scales are not comparable across waves" — false for the smooth reading | **confirmed by running** (`boeglin_check.py`, the harness's `Tree` and `scaled_chi2`): ∫4πk² n dk over 0–1.19 GeV = 0.99992 / 0.99983 / 0.99998 read smoothly (1.0828 / 1.0856 cell-wise); smooth scales 0.5015 / 0.5866 / 0.6714; Hulthén / CD-Bonn = 1.243 / 1.495 / 1.403 at k = 0.2 / 0.3 / 0.4 GeV; (r1) × (r5) = 0.902 × 0.855 = 0.77 | applied: phase C5 §4.1 corrected in place (marked) and §4.5 added (post hoc, not a verdict); `08` §3.6, §4 G and SUMMARY call row 16 a shape test whose free scale absorbs the Hulthén's ×1.34 |
| P-6 | physics, minor | `08` §0 still at the 2026-09-27 state ("every other data comparison of a shipped default fails"; the C0-zero item) | confirmed by reading (rows 6, 12, 16 pass; phase C4 §5) | applied: a dated note on §0's last line (no line shifted) |
| P-7 | physics, nit | the corrected Q_d reference 0.2860(15) is itself superseded: Puchalski, Komasa, Pachucki, PRL 125 (2020) 253001, 0.285699(15)(18) fm² | confirmed on INSPIRE (record 1822750, its abstract; the paper not read): the IA distances become −5.32 / −5.61 / +1.63 % | **applied in part**: the dated note in `t2_deuteron_static.py`'s TOLERANCE paragraph and in SUMMARY; the note in the JSON's `_provenance.corrections` **not applied** — the fix stage edits no vendored reference file (and it would re-pin FILE_SHA256); listed under "What the run did not apply" |
| P-8 | physics, nit | `08` §4 I reads row 10's FAIL as the weight's; the band is the fitted mean's error, which a shaped model fails even where the data prefer it | confirmed by reading the harness's rule (b) and (d) (δ = S·σ(R̄)/R̄) and its own Δχ² = 513 / 177 | applied: `08` §4 I |
| H-9 | honesty, minor | "the paper's 20 % allowance": the paper gives a 10–20 % FSI estimate and the shipped Hulthén passes only at its upper end | **confirmed by running** and reading (Boeglin p. 1: "all predict small FSI contributions (10 - 20%) for 35° ≤ θnq ≤ 45°"): A = 0.10 gives 51.57 / 27.81 / 19.95 on 21 (p = 2.2 × 10⁻⁴ / 0.146 / 0.525), A = 0.15 32.31 (p = 0.054) / 17.32 / 12.86, A = 0.20 the recorded 21.98 / 11.54 / 8.89 | applied in the eight places named (SUMMARY, this file's C4/C5 row, DECISIONS_PRICED, `BENCHMARK_PLAN.md` ×2, `08` §3.6, `02_data_nucleon_deuteron.md`, the README) and `08` §3.7; the dated post-hoc sentence in phase C5 §4.5 and DECISIONS_PRICED row 9; the hashed design block untouched |
| H-10 | honesty, minor | T4's header still on 8.312 without 0.076; T5's "never run here" / BLOCKED for the NNPDFpol11 sub-row that now runs and fails; the README lists T4's 0.076 as seen before its tolerance (0.20 dates from 2026-09-26) | confirmed by reading (0.076 = W2 n_ch_eta[−2,0), 0.924 ± 0.062, 776 stock particles; ⟨n_ch⟩ 1.007 / 0.997 — REPORT.json and `phase_A_baseline.md` §2.3; NNPDFpol11 0.2087 / 0.2084 / 0.2083 — REPORT.json) | applied: `t4_pythia_ep_closure.py` (the version sentence; a MEASURED 2026-10-26 block), `t5_sum_rules.py` (the LhapdfG1 NOTE's measurement; "blocked here until 2026-10-26"; the 2026-09-26 block's BLOCKED "on that machine"), the README's tolerance paragraph; `phase_A_baseline.md`'s two citations re-pointed to the shifted lines (`BJ_CONV_REL` 207 → 214, `bjorken_verdict` 462-472 → 469-479; content identical) |
| H-11 | honesty, minor | SUMMARY, DECISIONS_PRICED and this file say no reference number moved; Q_d's did (PLAN D2) | confirmed (`git diff master` on `deuteron_static_measured.json`) | applied in all three (= C-26, R-34) |
| H-12 | honesty, minor | DECISIONS_PRICED's five lookup-bias numbers at the Deeps momenta are in no record | **confirmed by running**: +5.71 / +4.78 / +5.58 / +2.39 / +2.97 % (Hulthén), +5.22 / +3.71 / +3.30 / +1.62 / +3.76 % (AV18) | applied: provenance and the one-line command there |
| H-13 | honesty, minor | the report's legend "one of the first ten" against plan rows 11–16 | confirmed | applied in `make_report.py` ("when the harness wires a §4 row"), the report regenerated (= C-27, R-38) |
| H-14 | honesty, minor | phase C1's plain-text path:line citations into `08_channel_priority.md` and `REFERENCES.md` drifted under phase E's edits | confirmed (`git show 2d89a64:` of `08` lines 91, 117, 308 and `REFERENCES.md` 136, 140 against the tree) | **applied in part**, with R-32: the drifted ones re-pointed in plain text with the commit they held at (`08` :91 and :308 as of `2d89a64`, row A at line 327 since; :117 reworded in phase F; `REFERENCES.md` 136 → 137, 140 → 141); the phase-E records-gate row of this file notes the re-read |
| H-15 | honesty, nit | `PLAN.md` says the T4 pin re-take is in `phase_0_environment.md` with the sub-rows; that record says `82b3237` changed no number | confirmed (`git show 82b3237 -- python/tests/test_bench_t4_t6.py`: 0.023 → 0.076) | applied: `PLAN.md` A0 and `phase_0_environment.md`, both marked |
| H-16 | honesty, nit | SUMMARY's price lines drop the `li6-convolution` label, merge registry row 3 with the claim-1 R choice, and say "AV18 deuteron" | confirmed | applied (with P-1, P-3) |
| H-17 | honesty, nit | the board and SUMMARY say phases A–E are uncommitted | confirmed (HEAD = origin = `6b1b63e`) | applied: the nine commit cells, the header sentence, SUMMARY (= C-30, R-35) |
| H-18 | honesty, nit | "three documents" names two; the deuteron JSON's two APS values left out | confirmed (its `source_url`) | applied: DECISIONS_PRICED §3 item 8, SUMMARY item 4 |
| H-19 | honesty, nit | "no file was written into the sibling repository": the `polligen` fixture leaves 22 gitignored `.pyc` | confirmed (all 22 stamped 2026-10-26 13:16:14 CDT; the sibling's `git status` clean) | applied: DECISIONS_PRICED §3 item 7 reworded; the optional `sys.dont_write_bytecode` not applied (`conftest.py` is not a file this run touched) |
| H-20 | honesty, nit | the WS98 and Boeglin CSVs' Licence lines ("no … text … is reproduced") against their quoted captions and sentences | confirmed by reading both headers | applied as the finding proposes until the files' next revision: noted in phases C4 and C5 §8, DECISIONS_PRICED §3 item 8, SUMMARY item 4; the files and their digests untouched |
| H-21 | honesty, nit | README: `test_sr_*` "(14 …)" | confirmed (`--collect-only`: 13, plus `test_t5_sum_rules_plot_writes_png`) | applied (= part of R-36) |
| H-22 | honesty, nit | CRSZ25's DOI is null too, not only KS20's | confirmed (`body.sources`) | applied: the not-applied row names both |
| C-23 | code, minor | row 13 reports a missing LHAPDF set as a plain BLOCKED row with exit 0, against README convention 6 | **confirmed by running** the verifier's scenario (`LHAPDF_DATA_PATH='<set-less dir>::'`, `make_report.py --only t5_delta_moment_bound,t3_nmc_li6_over_d`: 1 blocked, 1 blocked (environment)) | **applied in part**: `check_set`'s not-installed raise is now LHAPDF's own `RuntimeError("Info file not found for PDF set '…'")`, so `main()` exits 2 and the report reads blocked (environment) — the same scenario now gives 0 blocked, 2 blocked (environment); `PluginPdf`'s "could not load" **stays** a BLOCKED row — with C-31's search order a missing set is caught by `check_set` first, and a plugin that cannot load a set LHAPDF finds is the PYTHIA tier's failure, BLOCKED as in T4/T6; the header's exit paragraph, README convention 6 and the tests (`test_a_missing_set_is_blocked_by_the_environment`, `test_a_set_of_another_version_is_a_blocked_row`; the `rep` fixture and `test_exit_codes` skip on, or accept, exit 2 for a missing set) |
| C-24 | code, minor | no gate checks that the figures the committed report links exist | confirmed (`--check` copies the figure fields, never opens a PNG) | applied: `test_bench_report.py::test_the_committed_report_links_only_figures_that_exist` (not opt-in; it fails on a copy of the report with one PNG removed) |
| C-25 | code, minor | §7's change → tier map: a change to `fsi.*` re-runs no T2 row (10, 11), to `rc.*` no T3 row (5, 6, 15, `t3_li6_rms_devries`) | **confirmed by running** `tiers.py --files` on each file alone | **not applied**: `tiers.py` is not a file this run touched, and it must change with `BENCHMARK_PLAN.md` §7 together (`test_section7_as_data_matches_the_plan`); E-2 and E-4 stand as the exact text |
| C-26 | code, minor | = H-11 | — | — |
| C-27 | code, nit | = H-13 | — | — |
| C-28 | code, nit | README convention 4: one line of the code block indented 4, not 7 (CommonMark ends the block there) | confirmed (`cat -A`) | applied |
| C-29 | code, nit | README row 11: `test_bench_deeps.py` "(15)" | confirmed (16 collected) | applied (with R-36) |
| C-30 | code, nit | = H-17 | — | — |
| C-31 | code, nit | `lhapdf_dirs()` always appends the install datadir, unlike LHAPDF (`$LHAPATH` fallback, a trailing `::`) | confirmed (the C-23 scenario: `check_set` certified a copy the plugin never loaded) | applied, with `test_lhapdf_dirs_follow_lhapdf_search_order` |
| R-32 | regression, minor | phase C1's 18 citations are plain text, so R6 never reads them; backtick all and annotate the drifted and the tokenless | confirmed | **applied in part**: the seven that pass R6 unannotated are backticked (the records gate: 838 citations, 0 broken; 343 against a live file); the drifted ones fixed in plain text (H-14). The full form — four `(anchor read 2026-10-26 "…")` pins and seven `(as of 2d89a64)` — was applied first and **reverted**: the first full pytest of this phase failed `test_doc_link_gate.py::test_every_anchor_read_pin_in_the_records_is_still_on_its_line` (`pins == 62`, now 66) and `::test_every_as_of_anchor_passes_the_bar_at_its_commit` (the `08` :91 anchor shares only "sweet" with its line at `2d89a64`, so it would need an `ASOF_READ_BY_HAND` entry) — a census in a file this run did not touch |
| R-33 | regression, nit | phase 0's pytest row: "the T4/T6 tests may have seen either state — both pass" | **confirmed by running** master's tree (`git archive master` into the scratch directory, `test_bench_t4_t6.py` from its root: 15 passed, 3 skipped — the two pythia8 tests skip) | applied: a marked correction in that row |
| R-34 | regression, nit | = H-11 | — | — |
| R-35 | regression, nit | = H-17 | — | — |
| R-36 | regression, nit | stale test counts: README rows T1 (`test_f2d_*` 6, `test_dp_*` 5, `test_g1d_*` 8), 11 (15), T5 (14); plan row 11 (15) | confirmed (`--collect-only`: 7 / 6 / 11; 16; 13 + the figure test) | applied; also row 13's 21 → 23 (this stage's two new tests) in the plan and the README |
| R-37 | regression, nit | CI installs no matplotlib, so every `plot()` test skips there | confirmed (`.github/workflows/ci.yml`'s two pip lines) | applied as the finding's alternative: README convention 7 says so; `ci.yml` not edited (not a file this run touched) |
| R-38 | regression, nit | = H-13 | — | — |

**Gates, on the final tree** (this machine, `source env.sh`; one at a time,
nothing else running):

| gate | result | when (UTC, 2026-10-27) |
|---|---|---|
| `LIPOLGEN_BENCH=1 python3 validation/benchmarks/make_report.py --all --plots --compare <6b1b63e's REPORT.json> --out-dir docs/benchmarking` | **23 harnesses: 13 pass, 10 fail, 0 blocked, 0 blocked (environment), 0 error, 0 skipped (opt-in); compare: 0 changes**; the 19 figures redrawn — pixel-identical to `6b1b63e`'s except the footer (date, revision); `REPORT.md` differs in the header, the runtimes and the legend sentence | 03:16:27–03:17:41 (74 s) |
| `python3 validation/benchmarks/make_report.py --check` | **OK** | 03:17:47–03:18:50 (63 s) |
| `python3 -m pytest python/tests -q` | **1456 passed, 155 skipped, 0 failed**, 3 warnings (Axes3D; two `IntegrationWarning`s, as in phases 0 and E). The first run (03:24–03:32) failed the two census tests of R-32, before the annotations were reverted | 03:37:16–03:45:06 (470 s) |
| `LIPOLGEN_BENCH=1 python3 -m pytest python/tests -q -k bench` | **351 passed, 1 skipped, 0 failed** (1259 deselected; the skip is the nnpdf-wheel reproduction); phase E's 348 + this stage's three tests | 03:45:19–03:50:22 (303 s) |
| `python3 validation/check_physics_channels_links.py` (strict) | **exit 0** — `PHYSICS_CHANNELS.md` 1240 / 95 ranges / 7 external / **0 broken** / 24 allow-listed; `SPIN32_FINITE_GAMMA.md` 19 / 115 / 0 / **0** / 10; `PYTHIA_BRIDGE.md` 0 / 0 / 8 / **0** / 0 | 03:37:04 |
| `python3 validation/check_physics_channels_links.py --records-only` | **exit 0** — 838 citations in 67 records: 343 against a live file, 57 into another record, 110 historical, 129 read by hand, 199 skipped, **0 broken** (830 before: + the seven of R-32, + one range in DECISIONS_PRICED §3 item 7) | 03:37:06; again after this section was written |
| `python3 validation/check_physics_channels_links.py --audit-ranges` | **exit 0** — 209 range citations on 179 blocks, **0 REFUSED** | 03:37:07 |
| `python3 validation/check_spdx_headers.py` | **exit 0** — 133 files | 03:37:09 |
| the C++ suite | not re-run: no C++ file changed since phase E's 416 / 416 | — |

The scripts the fix stage ran (the verifiers' `rprice.py`, `f2price.py`,
`charm.py`; its own `boeglin_check.py`; the exact-replacement edit scripts)
are in the session scratch directory, not in the repository.
