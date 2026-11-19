<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-10 — phase P3: harness rows (decisions J3, J4, J7)

Written by the P3 implementer (an Opus worker under the supervising session)
on 2026-11-12, on base commit `8b72cd5` with a clean tree (phases P1, P2, P2b
of this run committed; this phase uncommitted, the supervising session
commits). It implements **J3, J4 and J7 — the SUPERVISOR's choices under the
author's delegation** ("for the others go with your recommendation";
`DECISIONS.md` J3, J4, J7), never the author's words. Scripts and raw logs
are in the session scratch `run_1110/p3/` (volatile); every number a reader
needs is below. No C++ file changed, so no build ran and the C++ binary was
not run.

## 0. What changed, in one paragraph per entry

**J3 (BONuS, row 10).** The 2026-11-06 hand check — no block verdict of the
AV18 headline flips under the 321 × 321 weight table or under bin averaging —
is now a bench-gated test,
`python/tests/test_bench_bonus.py::test_no_headline_block_verdict_flips_under_the_finer_table_or_bin_averaging`.
It calls the harness's own `judge_block` on the harness's own block
selection and varies only the weight function. Two small hooks were added to
`validation/benchmarks/t2_bonus_spectator_shape.py`, both parameters and no
copy of the logic: `verdict_blocks(doc)`, which is rule (a)'s selection that
`measure()` used inline and now calls; and `Tree.weight(..., fsi=None)`,
which takes a variant weight object. The harness header carries a dated
`J3 NOTE`, and README row 10 a dated clause. 0 of 28 blocks flip. The
thresholds (5e−4, 2 %), the 161 × 161 table, the headline and the status are
unchanged.

**J4 (HERMES b₁, row 3).** One RECORDED sub-row sits beside the six-bin
headline:
`miller_shipped[in-grid: <Q2> >= the shipped grid's Q2 floor] (J4, recorded)`.
It holds the headline model's χ²(b₁) and χ²(A_zz) (shipped F₁ᵈ) on the bins
at or above the shipped CT18NLO grid's own Q² floor. The floor is read from
`make_lhapdf_unpol_sf(UNPOL_SF_DEFAULT).q2_min`, the process's one shared
object, which is the one `shipped_inclusive_kernel` hands the kernel. It was
checked to be the same Python object on a second call, and it is never
typed. The label says the sub-row was chosen after the values were seen and
is a diagnostic with no verdict. It also names the three bins below the
floor and states what the frozen edge does there. In the figure, the
below-floor HERMES points are ringed by an open grey square with a legend
entry (6 lines in `plot()`). The harness has a dated `J4 NOTE` after the J1
NOTE, outside the first 160 lines. README row 3 has a dated clause and
`test_bench_spin_deuteron.py` has a new pin. The headline (all six bins, its
rule, FAIL) is unchanged. `docs/benchmarking/02_data_nucleon_deuteron.md`
states none of the row's numbers, so it was not edited.

**J7 (`t5_delta_moment_bound`, row 13).** A verdict set can be installed but
unusable: its `.info` or member file differs from the vendored sha256, or its
member file is absent. Such a set is now treated as a set that is not
installed. `check_set` raises a `RuntimeError` naming the set, the directory
and what differs, `main()` exits 2, and `make_report` records
`blocked (environment)` with that sentence as the reason. Until now these
cases gave a `BLOCKED` row with exit 0. The tier-off, module-missing and
PYTHIA-refused branches keep their `BLOCKED` row with exit 0.

For the report to say this, the shared classifier needed one branch.
`validation/benchmarks/_report.py` gained the constant
`LHAPDF_NOT_AS_VENDORED = "is installed but not usable as vendored"`, and
`classify_exception` returns `(blocked (environment), <the message's first
line>)` for an exception that starts with `environment: LHAPDF set ` and
carries that phrase. Without that branch the report would have said `error`
or a false "not installed". J7's "what changes" list does not name
`_report.py`, so this is listed in §5. `make_report.py`'s docstring gained one
clause on the same line. README convention 6 has a line-neutral dated clause.

No other harness has the same two-outcome split. No `t1_*`, `t3_*`,
`t4_*` or `t5_sum_rules` harness checks an LHAPDF set's sha256 or version. They
only distinguish a missing set ("Info file not found"). `tiers.py`'s
`STANDING_NOTE` ("any PDF-set version ... re-runs T1, T2, T3-EMC by hand") has
a different, documented meaning — re-run guidance, not a check — and was left
alone.

## 1. Every number or label that moved, before → after

### J3 — the robustness pin (recipe: the test itself; exploration `scripts/j3_flip.py`)

| quantity | value |
|---|---|
| verdict blocks (`verdict_blocks`) | 28 (12 at 4 GeV, 16 at 5 GeV), the report's order and verdicts |
| headline passes | 3 / 28 (4 GeV Q² = 3.38, W* = 2.03, p_s = 0.078 / 0.093 / 0.11) |
| (i) 321 × 321 on [0, 0.6] GeV, AV18, 40 mb: blocks whose verdict differs | **0 of 28** (3 pass, the same three) |
| (ii) weight averaged over each p_s bin × cos bin (10 × 10 midpoints): blocks whose verdict differs | **0 of 28** (3 pass, the same three) |
| closest margins, \|v/δ − 1\| | headline 0.0272; (i) 0.0328; (ii) 0.0241 |
| the block nearest the threshold | 4 GeV Q² = 3.38, W* = 2.03, p_s = 0.135: v/δ 1.0272 (fail) → 1.0328 / 1.0241 (fail); next, p_s = 0.11: 0.9583 (pass) → 0.9643 / 0.9662 (pass) |
| test wall time | 2.91 s (call); bench-gated (`@pytest.mark.bench_expensive`), so it runs in the bench subset / `LIPOLGEN_BENCH=1` only and skips in the plain suite |

The definition of "a block's verdict" (the harness's `judge_block(b,
wfun)["passed"]`, rules (b)–(d), stat errors, headline weight, only wfun
varied) is written in the test's docstring. In the order of the work, an
exploration script (`scripts/j3_flip.py`) was run BEFORE the test was written.
It computed the same quantities, so the numbers above were seen before the
docstring existed. The definition is not a choice: J3 and the brief fix it as
the harness's own rule. Nothing was tuned. The run-2026-11-03 hand check is
described in `../run_2026-11-03/phase_P3_report_rules.md` §1.2.1 (BONuS), not
§1.2.2 (deuteron statics) as the brief says.

The harness's printed output is **byte-identical** before and after (the
`8b72cd5` harness run from a scratch copy against the same data, `cmp`), and
the report row did not move (§2).

### J4 — the in-grid sub-row (recipe: `t2_hermes_b1_table2.run()`; `scripts/j4_explore.py`)

Floor: `make_lhapdf_unpol_sf(UNPOL_SF_DEFAULT).q2_min` = 1.677025 GeV²
(CT18NLO). Bins with ⟨Q²⟩ below it: **0.51, 1.06, 1.65 GeV²** (x = 0.012,
0.032, 0.063), as expected. For each of them the shipped kernel's F₂ᵈ at
⟨Q²⟩ equals its value at the floor bit for bit (`f2d_frozen_at_floor` true).

| x | ⟨Q²⟩ GeV² | below floor | shipped F₁ᵈ (HERMES's construction, CT18NLO F₂, `r_sigma_lt`) | implied F₁ᵈ = −(2/3) b₁/A_zz (HERMES's columns) | ratio |
|---|---|---|---|---|---|
| 0.012 | 0.51 | yes | 12.402178 | 7.044025 | 1.760666 |
| 0.032 | 1.06 | yes | 4.843366 | 3.426791 | 1.413382 |
| 0.063 | 1.65 | yes | 2.582180 | 1.929293 | 1.338407 |
| 0.128 | 2.33 | no | 1.298265 | 1.017544 | 1.275881 |
| 0.248 | 3.11 | no | 0.563796 | 0.495726 | 1.137313 |
| 0.452 | 4.69 | no | 0.145109 | 0.161359 | 0.899293 |

The sub-row (new; nothing existed before):

| quantity | headline, 6 bins (unchanged) | in-grid sub-row, 3 bins (J4, recorded) |
|---|---|---|
| χ²(b₁) / ndf, bins within | 5.261863562 / 6, 4/6 | **1.894165012 / 3, 2/3** |
| χ²(A_zz) / ndf on the shipped F₁ᵈ, bins within | 7.307981761 / 6, 3/6 | **1.177692945 / 3, 3/3** |
| status | FAIL | recorded (no verdict) |

The in-grid bins also differ from the implied F₁ᵈ: −10 % to +28 %. The
frozen edge is not the only source of the mismatch. It is the largest
(+34 % to +76 %).

### J7 — exit codes (recipe: the harness with `LHAPDF_DATA_PATH=<scratch copy>::`)

| situation | before (`8b72cd5`) | after |
|---|---|---|
| right sets installed | row as printed, exit 0 | **byte-identical stdout**, exit 0 |
| CT18NLO `.info` one byte changed (scratch copy) | `BLOCKED` row, exit 0 (`... at <dir> is not the version this row was measured with (sha256 info 6c719671301f / member 375db856d2f8, vendored be60232d8e6c / 375db856d2f8)`) | **exit 2**, no row: `environment: LHAPDF set CT18NLO at <dir> is installed but not usable as vendored: .info not the version this row was measured with (sha256 info 6c719671301f / member 375db856d2f8, vendored be60232d8e6c / 375db856d2f8)`; report `blocked (environment)`, that sentence the reason |
| NNPDFpol11_100 member file one byte changed (pytest) | `BLOCKED`, exit 0 | **exit 2**: `... NNPDFpol11_100 at <dir> is installed but not usable as vendored: member file NNPDFpol11_100_0000.dat not the version ...` |
| CT18NLO `.info` present, member file absent | `BLOCKED` row, exit 0 (`... member file <path> not installed`) | **exit 2**: `... CT18NLO at <dir> is installed but not usable as vendored: member file CT18NLO_0000.dat absent`; `blocked (environment)` |
| a set not installed | exit 2, `Info file not found for PDF set 'CT18NLO' (searched <dir>)`, `blocked (environment)`, reason `environment: LHAPDF set CT18NLO not installed` | unchanged |
| tier off / `pythia8` missing / PYTHIA refuses a set | `BLOCKED` row, exit 0 | unchanged |

### Labels, documents and pins (each with its dated comment)

- `t2_bonus_spectator_shape.py`: the `J3 NOTE` after the A5 NOTE, plus the
  `verdict_blocks` and `fsi=` hooks (`J3 (2026-11-10)`).
- `t2_hermes_b1_table2.py`: the `J4 NOTE` after the J1 NOTE, the constants
  `INGRID_SUBROW` / `INGRID_LABEL`, and `ingrid_subrow()`. `per_bin` gained
  `q2_grid_floor`, `below_grid_floor`, `f2d_frozen_at_floor` and
  `f1d_shipped_over_implied`. The report gained `subrows`. There is one
  printed plain line (not a `REPORT |` or `SUBROW |` line, so convention 4's
  "row 3 prints only its six configuration lines" still holds), and the
  figure's open squares.
- `t5_delta_moment_bound.py`: the header A4 sentence, the exit-status
  paragraph, the `EnvironmentBlocked` and `check_set` docstrings and the two
  `check_set` branches were rewritten line-neutrally. A `J7` comment block
  follows the A7 note.
- `_report.py`: `LHAPDF_NOT_AS_VENDORED` and the `classify_exception` branch.
  The constant shares an existing line, so the line-cited part of the file did
  not shift. `make_report.py` gained a docstring clause on one existing line.
- `validation/benchmarks/README.md`: convention 6 has the J7 clause, row 3
  the J4 clause (test count 7 → 8 `test_hermes_*`), and row 10 the J3 clause
  (`test_bench_bonus.py` 26, + 1 bench-gated). All edits are line-neutral.
- Tests:
  - `test_bench_bonus.py` gained the J3 test, appended at the end of the file
    with module constants `J3_*`.
  - `test_bench_spin_deuteron.py` gained `test_hermes_in_grid_subrow_recorded`
    (J4).
  - `test_bench_t5_delta_bound.py`: `test_a_set_of_another_version_is_a_blocked_row`
    pinned the old two-outcome behaviour, so it was replaced. Its place is
    taken by three tests:
    - `test_a_tier_or_module_that_is_missing_is_still_a_blocked_row`: the
      unchanged `EnvironmentBlocked` branch, exit 0.
    - `test_a_set_of_another_version_exits_2_blocked_by_the_environment`: a
      `.info` copy and a member copy, each with one byte changed.
    - `test_a_set_without_its_member_file_exits_2_blocked_by_the_environment`.

    Each of the two new exit-2 tests checks `run()`'s sentence, the
    classifier, `main() == 2`, the script's exit 2 with no `REPORT |` line,
    and `make_report.run_harness`'s row. The installed sets are only read:
    the copies sit in a pytest `tmp_path` with `LHAPDF_DATA_PATH=<copy>::`.
    The module docstring was edited line-neutrally.

## 2. What was verified and how

- **BONuS headline unchanged.** `t2_bonus_spectator_shape.py` stdout is
  byte-identical, comparing `8b72cd5`'s harness (a scratch copy against the
  same `data/`) with the edited one.
- **The J3 test detects what it claims.** Its recomputed headline verdicts
  equal `rep["measured"]["judged"]["av18-40mb"]` in block order, and the
  variant verdicts are compared block by block. Its result reproduces the
  exploration script's 0 / 0.
- **HERMES.** The printed `REPORT |` lines and per-bin tables are identical
  to `8b72cd5`'s run. Only the J4 line was added (diff of the two outputs).
- **`t5_delta_moment_bound` with the right sets.** Stdout is byte-identical
  (`cmp`). The three broken environments were run on scratch copies with both
  the `8b72cd5` harness and the edited one (the table in §1).
- **Report.** `make_report.py --all --plots --out-dir docs/benchmarking
  --compare <pre copy of REPORT.json>` printed one line,
  `SUBROW t2_hermes_b1_table2[miller_shipped[in-grid: ...] (J4, recorded)]: new (informational)`,
  and `compare: 0 changes`, exit 0. Status stays 15 pass / 8 fail / 0 blocked.
  A JSON diff of the old and new `REPORT.json` shows only `generated.date`,
  `generated.git_rev`, the runtimes and that one new sub-row. REPORT.json and
  REPORT.md are kept, because they carry the new sub-row. Then
  `make_report.py --check` printed "OK".
- **Figures.** 19 PNGs were rewritten by `--plots`. 18 differ from `HEAD`
  only in the footer stamp: every differing pixel lies in the bottom 2.0–2.8 %
  of the height (`pixdiff.txt`), so they were restored with
  `git checkout -- docs/benchmarking/figures/<file>`.
  `t2_hermes_b1_table2.png` changed in content (the J4 squares and legend
  entry) and is kept.
- **Citations.** All three modes give the same counts as before the edits
  (strict 1243 / 95 / 7 / 0 broken / 24; records-only 0 broken; audit-ranges
  209, 0 REFUSED). Edits to cited files were kept line-neutral, or appended
  after every cited line.

## 3. Files

`validation/benchmarks/t2_bonus_spectator_shape.py`,
`validation/benchmarks/t2_hermes_b1_table2.py`,
`validation/benchmarks/t5_delta_moment_bound.py`,
`validation/benchmarks/_report.py`, `validation/benchmarks/make_report.py`,
`validation/benchmarks/README.md`, `python/tests/test_bench_bonus.py`,
`python/tests/test_bench_spin_deuteron.py`,
`python/tests/test_bench_t5_delta_bound.py`, `docs/benchmarking/REPORT.json`,
`docs/benchmarking/REPORT.md`,
`docs/benchmarking/figures/t2_hermes_b1_table2.png`, this record,
`STATUS.md` (the P3 row). No C++ source, header, test or binding was
changed, and no vendored data or reference JSON either.

## 4. Gates (machine shared, load ~1–3)

| gate | result | wall |
|---|---|---|
| build / `lipolgen_tests` | not run: no C++ change | — |
| pytest, six parts | A 359 passed / 4 skipped; B 138 / 1; C 186 / 0; D 249 / 90; E 277 / 61; F 324 / 0; **sum 1533 passed / 156 skipped / 0 failed** (1530 / 155 at `8b72cd5`: +1 J4 pin, +2 J7 tests net, +1 J3 skipped without the opt-in) | 171 / 96 / 68 / 341 / 248 / 162 s |
| bench subset (`LIPOLGEN_BENCH=1 -k bench`) | **372 passed / 1 skipped** (368 / 1 at `8b72cd5`: +1 J3, +1 J4, +2 J7 net; the J3 test 2.96 s there) | 383 s |
| report `--all --plots --compare` | 0 changes, 1 new informational sub-row; 15 / 8 / 0 | 99 s |
| report `--check` | OK | 87 s |
| citations strict / `--records-only` / `--audit-ranges` | 1243 / 95 / 7 / 0 broken / 24; 1432 citations in 96 records, 0 broken; 209, 0 REFUSED | ~3 s each |
| `check_spdx_headers.py` | 136 files, all carry the line | — |
| `test_doc_link_gate.py` | 112 passed (after this record and the board row were written) | — |

## 5. Open points

1. **`_report.py` is outside J7's list.** J7's "what changes" list names the
   harness, its pytest and README convention 6. It does not name
   `_report.py` / `make_report.py`. The required outcome ("`make_report` shows
   the row as `blocked (environment)` with the reason naming the set, the
   directory and what differs") could not be had without one classifier
   branch: before it, the new exception would have read `error`. The branch is
   keyed on the harness's phrase, so no other harness's classification
   changes.
2. **J3's numbers were seen first.** They were seen in an exploration run
   before the test's docstring was written. The definition is the harness's
   own rule, as J3 prescribes, so nothing was chosen after seeing.
3. **J3's margin is thin.** The block nearest the threshold (4 GeV, Q² =
   3.38, W* = 2.03, p_s = 0.135) sits at v/δ 1.024–1.033 under the three
   weights. A future change of the weight of a few percent could flip it, and
   the test will then fail. That is its purpose.
4. **J4: the shipped F₁ᵈ is off in-grid too.** The shipped F₁ᵈ differs from
   HERMES's implied F₁ᵈ by up to +28 % even on the in-grid bins (1.276 at
   x = 0.128). The sub-row isolates the frozen edge but not the LO CT18NLO /
   `r_sigma_lt` construction. This is recorded, not acted on.
5. **`BENCHMARK_PLAN.md` row 3 was not edited.** It quotes the headline's
   numbers but not a sub-row. The brief did not name it.
