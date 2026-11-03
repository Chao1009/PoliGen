<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-03 — phase P1: the inclusive polarization factors (decision A1)

Worker P1, 2026-11-03, on `claude/run-2026-10-26` at `8eea2c7` (uncommitted;
the supervising session commits). Decision: `DECISIONS.md` A1 — the
SUPERVISOR's choice under the author's delegation ("Go with your
recommendations for the 17 items, just record the choice and why"). Every
number below was measured in this phase with `source env.sh` on this machine
(4 cores, one long gate at a time); scripts and logs are in the session
scratch `run_1103/p1/` (`ab_runs.sh`, `ab_compare.py`, `ab_ratios.py`,
`model_dump.py`, `pre/`, `post/`, `flip_*.log`, `final_*.log`, `report.log`).
Conventions: **pre** = `8eea2c7`'s build, taken before any edit; **post** =
this tree.

## 0. What changed, in one paragraph

The three factors an INCLUSIVE run reads moved exactly as A1 names, each with
one home and its pre-A1 value under a name. `LI6_CLUSTER_POLARIZATION`
(include/lipolgen/beams.hpp) is now DEFINED as `LI6_POLARIZATION_VMC_SIX_BODY`
= 0.848 (was the cluster product 0.811228375, now
`LI6_CLUSTER_POLARIZATION_HULTHEN`, bit for bit the same expression);
`LI6_B1_RANK2_TRANSFER` (include/lipolgen/constants.hpp) is
`LI6_B1_RANK2_TRANSFER_VMC` = 0.982576 (was 0.921947, now
`LI6_B1_RANK2_TRANSFER_HULTHEN`), with the unquantified systematic stated
beside it (the 4.4 % non-cluster depletion of the vector sector, cluster 0.887
→ six-body 0.848, has no rank-2 counterpart and none is applied);
`P_D_DEUTERON` = 0.0576, AV18's (was 0.045, now `P_D_DEUTERON_SCENARIO`), so
`DEUTERON_VECTOR_POLARIZATION` 0.9325 → 0.9136
(`DEUTERON_VECTOR_POLARIZATION_SCENARIO` keeps 0.9325). The kernels take the
factors as PARAMETERS with the new values as defaults — the vector factors
through the `Ion` (`li6_ion(P)`, `deuteron_ion(P_D)`, `DEUTERON_SCENARIO()`,
beams.hpp / src/core/beams.cpp), the transfer through a new last argument
`b1_transfer` of `default_inclusive_kernel` (pipeline.hpp; Python
`default_inclusive_kernel(..., b1_transfer=)` and
`lipolgen.shipped_inclusive_kernel(ion, b1_transfer=None)`) — and every polligen
port gate passes polligen's values explicitly with a dated comment. The Hulthén
TAGGED physics is wired to the scenario names (`li6_alpha_channel`'s Hulthén
`dis_target` = `DEUTERON_SCENARIO()`, `deuteron_channel`'s and
`BreakupOptions::p_d`'s default = `P_D_DEUTERON_SCENARIO`, the `Pipeline`'s
d–p factory call, the Python `deuteron_channel` default), so
`--cluster-wave hulthen` reaches the old physics bit for bit and the Hulthén
cluster reading is exactly 0.811228375. No CLI knob was added. Provenance: an
inclusive npz carries `meta["li6_cluster_polarization"]`, `["deuteron_p_d"]`
and `["li6_b1_rank2_transfer"]` where each is read and known (one definition,
`inclusive_polarization_factors`, pipeline.hpp), and the CLI prints a line
"inclusive polarization factors: …". No reference JSON was edited (`git
status validation/reference` clean) and `dump_polligen_reference.py` was not
run. Registry row 10 (run_2026-09-03/STATUS.md) is closed with a pointer to A1.
**One consequence A1 does not name (§5 item 1):** the default inclusive ⁶Li
`--plan helicity-flip` (and `--plan apar`) at P_z = P_e = 0.7 is now REFUSED at
the window's top cell by the sampler's positivity check.

## 1. Every number that moved, before → after

| quantity (where) | pre-A1 | post-A1 | recipe |
|---|---|---|---|
| `LI6_CLUSTER_POLARIZATION` (beams.hpp) | 0.811228375 | 0.848 (= `LI6_POLARIZATION_VMC_SIX_BODY`) | constant |
| `LI6()` eff_pol_p = eff_pol_n | 0.270409458 | 0.282666667 | `model_dump.py` |
| `LI6_B1_RANK2_TRANSFER` (constants.hpp) | 0.921947 | 0.982576 (= `_VMC`) | constant |
| `P_D_DEUTERON` / `DEUTERON_VECTOR_POLARIZATION` | 0.045 / 0.9325 | 0.0576 / 0.9136 | constant |
| `DEUTERON()` eff_pol | 0.9325 | 0.9136 | `model_dump.py` |
| default ions' per-nucleon g₁(⁶Li)/g₁(d) | 0.290 | 0.309399 | (0.848/3)/0.9136 |
| CLI `B1_TOP_CELL_B1_OVER_F1` (cli.py) | cdks 48.67, miller 1.085 | cdks 51.87, miller 1.156 (li6-convolution 42.85 unmoved) | shipped kernel at x = 0.954993, Q² = 167.344 |
| test pins | — | every one re-taken with `A1 (2026-11-03)` | §3 |

Unchanged and verified: `LI6_CLUSTER_POLARIZATION_HULTHEN` = 0.811228375 (==
the pre-A1 constant), `LI6_CLUSTER_POLARIZATION_VMC` = 0.905427287 (it now
reads `P_D_DEUTERON_SCENARIO`), the Hulthén tag's dis_target 0.9325, the
Hulthén d–p control (vector 0.932496109, tensor 0.959488878), the VMC α–d
(vector 0.970966086, tensor 0.982575872), `deuteron_av18_p_d()` = 0.0575998920.

**The literal pins (A1 (b), (c)).** Declared tolerance 5e-7 absolute on each
— half a unit of the sixth decimal the literal is printed to — set by the
literal's rounding. Said plainly: the computed values had been seen before
the tolerance was written (the pre dump `pre/model.json` printed 0.98257587
and 0.05759989, and beams.hpp already quoted 0.9826 / 0.057600). Measured
(tests/test_tagged.cpp "A1: the shipped inclusive literals are the tree's own
numbers"): VMC `tensor_dilution()` − 0.982576 = −1.277e-7; `deuteron_av18_p_d()`
− 0.0576 = −1.080e-7. Both inside; neither literal is wrong. (The Hulthén
analogue, for scale: 7.84e-6 in beams.hpp's note.)

**Benchmark numbers that moved** (`make_report.py --compare`, §4 verbatim
there): t1_g1d_world (every with-D-state χ², all four backends; the theory is
the inclusive deuteron's `g1_nucleus(d)/2`, factor 0.9325 → 0.9136; the
no-D-state column unmoved): toy world 232.16 → 228.32, mstw 150.61 → 141.72,
ct18nlo 143.43 → 133.36, nnpdfpol (the headline) 148.53 → 140.65, compass
45.0 → 42.4 / hermes 52.6 → 49.6 / e143 51.0 → 48.6. t5_delta_moment_bound:
the bound's P_eff 0.270409 → 0.282667 (|M₂|/B 0.0760/0.0907/0.1023 →
0.0759/0.0906/0.1022) and the helicity-state rates' min W 0.913494 → 0.907806
(A_T [−0.173, 0.0865] → [−0.184, 0.0922]; b₁ ×1.06576). t5_weighted_unweighted:
N_U 2712940 → 2714846, σ_u 79.7924 → 79.8484, σ_w 79.7677 → 79.8237, exact
79.767466 → 79.823531, rate shift 1.5850e-2 → 1.6564e-2, no-op pull 18.61 →
19.42, pull +0.359 unmoved. The pre values of rc_low_x_anchor's HALF_WIDTH
(anchor 0.30: 1.358472494654602e-04, 1.4590673286701302e-04,
9.680015652051017e-05, 1.920690621629717e-05; 0.266: 1.2045122785937472e-04,
1.3085664856817475e-04, 8.798803601270058e-05, same; 0.113:
5.116913063199001e-05, 6.313126922340236e-05, 4.83334937275573e-05, same) all
moved ×1.0657619 with A_zz.

## 2. What was verified and how (A/B)

**Commands** (`ab_runs.sh`, 4000 events, seed 1234, CLI, pre and post): ⁶Li
inclusive tensor-thirds at configs 0/1/2, `helicity-flip`, `apar`, `azz`,
`transverse-tensor`, `tensor-flip`, `cos2phi`; deuteron inclusive
`helicity-flip` (REFUSED at the default x_max on pre already — vector term at
the top cell, A4's case — so also at `--x-max 0.95`) and tensor-thirds; ⁷Li
inclusive helicity-flip; coherent ⁶Li; tagged ⁶Li α (default, `--cluster-wave
hulthen`, helicity-flip), tagged ⁷Li α, tagged d–p (default, hulthen).
Columns compared byte for byte over every array (`ab_compare.py`).

* **Column-identical pre → post (9 / 9):** coherent ⁶Li, inclusive ⁷Li, tagged
  ⁶Li α (default, hulthen, helicity-flip), tagged ⁷Li α, tagged d–p (default,
  hulthen). No new meta key on any of them.
* **Explicit old values through the API reproduce PRE bit for bit:**
  `cfg.kernel = shipped_inclusive_kernel(li6_ion(LI6_CLUSTER_POLARIZATION_HULTHEN),
  b1_transfer=LI6_B1_RANK2_TRANSFER_HULTHEN)` against pre `inc6_tt_c1`, and the
  same with `deuteron_scenario()` against pre `incd_tt`: every column identical.
  At the default window the same API run of ⁶Li helicity-flip gives the pre
  category σ 496234.389 / 496308.089 pb (the pre CLI printed 496234 / 496308).
* **Cell-level ratios** (difference of the per-category σ, which isolates the
  sector; expected from the constants: vector ⁶Li 0.848/0.811228375 =
  1.04532833, rank-2 0.982576/0.921947 = 1.06576191, deuteron vector
  0.9136/0.9325 = 0.97973190):

| run | sector (categories) | measured ratio |
|---|---|---|
| ⁶Li tensor-thirds c0 / c1 / c2 | rank-2, σ(azz0) − σ(azz±) | 1.06576 (1.0657619147328 / …291 / …323) |
| ⁶Li `azz`, `tensor-flip` | rank-2 | 1.06576 / 1.06576 |
| deuteron tensor-thirds | rank-2 (the deuteron's class toy b₁ IS the ⁶Li `Li6B1`) | 1.06576 (1.0657619147338) |
| ⁶Li helicity-flip, x_max 0.95, default vs old API | vector, σ(apar+) − σ(apar−) | 1.04533 (1.0453283270) |
| deuteron helicity-flip, x_max 0.95, pre vs post CLI | vector | 0.979732 (0.9797319035) |
| sweet-spot cells (test_d7_shipped_kernel) | a₂ (cos 2φ, Δ) / w_avg at θ_S = 0 | 1 exactly / 1.06576 |
| kernel tables at six points | g₁ / b₁ / f₁, f₂, Δ | 1.04533 / 1.06576 / 1 |
| total σ, tensor-thirds | spin-averaged | 1 (to 1e-15) |
| total σ, transverse-tensor, cos2phi plans | b₁ at θ_S = 90° | 1.0000061 |

  Nothing moved that is not on A1's list, with ONE item to name: the
  inclusive DEUTERON's b₁ / A_zz moved ×1.06576 too. That is not a new choice:
  the class-level toy kernel fills b₁ = `Li6B1(MillerB1)` for every spin-1 ion
  (PHYSICS_CHANNELS.md row 205: "the bare deuteron included — a nucleus with
  no α core receives the ⁶Li transfer"), so the deuteron reads
  `LI6_B1_RANK2_TRANSFER` as ⁶Li does. Recorded, not changed (§5 item 2).
* **Banner / meta:** every post inclusive ⁶Li log carries "inclusive
  polarization factors: 6Li whole-nucleus vector polarization 0.848 (…); 6Li
  rank-2 transfer 0.982576 (…)", the deuteron's "deuteron P_D 0.0576, eff_pol
  0.9136 (…)"; the explicit-old API runs report 0.811228 / 0.045 and "rank-2
  transfer: caller-supplied kernel (not knowable from the run)".

## 3. Flip-only gate, the re-takes, and the knob matrix

The flipped tree before any pin was touched: `./build/lipolgen_tests` 413/428
(15 failed, 2443 assertions), pytest 60 failed / 1446 passed / 164 skipped
(21 min 35 s). One bug of mine was found BEFORE that pytest run, in the model
dump, and fixed: the Python `deuteron_channel` default `p_d` still read
`P_D_DEUTERON` (it moved the Hulthén d–p control to 0.0576); it reads
`P_D_DEUTERON_SCENARIO` now. During the flip pytest run I was already editing
C++ test files and documents, so its 13 `test_doc_link_gate` failures mix
the flip's line shifts with those edits; they are resolved by §4's gates.
Each failure, classified:

* **Pins on the old constants (re-taken, dated `A1 (2026-11-03)`):** C++ —
  test_beams "the effective polarizations…", "the 6Li slot…"; test_sf "the 6Li
  b1 transfer constants"; test_b1_nuclear T6; test_breakup "…D-state-diluted
  polarization"; test_tagged "6Li embedded-deuteron dilutions", T25, T26, T27.
  Python — test_module (4), test_bench_static (1), test_bench_t1_nucleon (6),
  test_bench_t5_closures (4), test_bench_t5_delta_bound (5), test_d7_shipped_kernel
  (2), test_b1_model top-cell note, test_price_defaults, test_rc_low_x_anchor,
  test_pzz_mode "measured price" (2), test_sf_backend's two σ pins (the inclusive
  tensor-thirds σ moved by one ulp: 591846.1688784091 → 591846.168878409 and
  496334.1559540481 → 496334.15595404804).
* **Port gates on polligen's dump (now pass the old values explicitly):**
  test_reference beams.json / xsec.json (a `polligen_ion` helper:
  `li6_ion(LI6_CLUSTER_POLARIZATION_HULTHEN)`, `DEUTERON_SCENARIO()`), test_tagged
  "channel construction / model grid against polligen" (`deut_model()` at
  `P_D_DEUTERON_SCENARIO`), T9 (`b1_default_li6.json`, transfer 0.921947 passed),
  test_b1_model's twin of T9, test_kernel_reference (5). Those twelve
  pre-identified gate ids pass on the final tree.
* **Neither — a FINDING:** 13 tests (test_pzz_mode 6, test_sf_backend 1, the
  knob matrix's `inclusive-6Li/helicity-flip` spec 2 + its 12 cells skipped)
  failed because the default inclusive ⁶Li
  helicity-flip run is refused (§5 item 1). Their subject is not that cell, so
  each now runs that one configuration at x_max = 0.95, through a helper
  appended to the file with the reason (`_x95` in test_pzz_mode.py and
  test_sf_backend.py; `Spec(..., x_max=0.95)` + `_spec_x_max` in
  test_knob_provenance.py). No gate threshold was moved; the refusal itself
  is untouched.
* **Load:** test_bench_t4_t6 `runtime_s < 30` (35.0 s while I was building
  and editing alongside); passes on the final, idle run.
* **Knob matrix, one to one** (ids `test_the_table_matches_the_output_hash[inclusive-6Li/helicity-flip|…]`):
  `pz:0.7 -> 0.5`, `pzz:0.6 -> 0.4`, `pzz:0.6 -> 0.4 at --pzz-mode typed`,
  `pzz_mode:ladder -> typed`, `pe:-> 0.35`, `rel_lumi_offset:0 -> 0.02`,
  `unpol_sf:ct18nlo -> toy`, `pol_sf:nnpdfpol -> toy`, `rc:off -> tensor-band`
  — pass (flip: skipped, base refused by positivity) → pass at x_max 0.95; the
  three `rc` cells of the spec skip on the final tree as their base is refused
  by `RcModel` (polarised beam), the reason they skipped before A1; and
  `test_every_knob_row_is_covered[…]`, `test_every_row_is_well_formed[…]`:
  failed on the flip → pass. Every other id of the matrix: unchanged outcome
  (the file: 515 passed / 2 failed / 160 skipped on the flip; the full-suite
  tallies are back to the pre-A1 1515 / 155 exactly).

## 4. Gates (final tree), tallies and wall times

| gate | result | wall |
|---|---|---|
| `cmake --build build -j2` | 0 warnings, 0 errors | 2 m 9 s |
| `./build/lipolgen_tests` | 430 / 430 test cases, 17 244 416 assertions, 0 failed | 7 m 19 s |
| pytest, full | 1515 passed, 155 skipped, 0 failed | 14 m 42 s |
| bench subset (`LIPOLGEN_BENCH=1 -k bench`) | 368 passed, 1 skipped (after the report regeneration; before it, 367 / 1 / 1 failed — `test_the_committed_report_is_what_regeneration_gives`, as expected) | 8 m 47 s |
| `check_spdx_headers.py` | 136 files, all carry the line | <1 s |
| `pytest python/tests/test_doc_link_gate.py` | 112 passed | 50 s |
| strict citation gate | 1241 / 95 / 7 / **0 broken** / 24 (SPIN32: 19 / 115 / 0 / 0 / 10; PYTHIA_BRIDGE 0 / 0 / 8 / 0) | 2 s |
| `--audit-ranges` | 209 range citations, 0 REFUSED | 2 s |
| `--records-only` | 1432 citations in 82 records: 922 live, 61 into another record, 121 annotated, 129 read by hand, 199 skipped, **0 broken** | 2 s |
| `make_report.py --all --plots --out-dir docs/benchmarking --compare REPORT_70982c9.json` | 14 pass / 9 fail / 0 blocked; 26 changes, all in three rows | 1 m 46 s |
| `make_report.py --check` | OK | 1 m 34 s |

Citations: the code edits of beams.hpp, beams.cpp, constants.hpp, pipeline.hpp,
pipeline.cpp, bindings.cpp (above line ~5190), the test files and the cited
documents are LINE-NEUTRAL wherever a citation reaches (`git diff --numstat`:
equal insertions and deletions in each), the new code sits at the ends of
files; `--fix` re-pointed nothing that stayed broken; one cited range was
edited and re-recorded after re-reading it — include/lipolgen/beams.hpp 99-100
(PHYSICS_CHANNELS.md row 302's "AV18-like P_D = 0.045 carries no primary
reference", still true of `P_D_DEUTERON_SCENARIO`, which the citation now
names) — and one cited block, docs/USAGE.md 163-533 (row 208's §2a band rules;
the edit is the miller row's transfer), re-recorded. The `(as of)` census is
untouched (121), `ASOF_READ_BY_HAND` untouched.

`--compare`: every status change is attributed to A1. Verdicts: none moved
(14 / 9). Sub-row status: `t1_g1d_world[mstw/e143]` FAIL → PASS — its declared
rule is p(χ², ndf) ≥ 0.01 on the tree's own deuteron g₁ (`g1_nucleus(d)/2`,
the shipped deuteron factor), and χ² 50.24 → 47.86 (p 0.0061 → 0.0111) with
the factor 0.9325 → 0.9136; it is a sub-row of an opt-in backend, the row
stays FAIL. MOVED rows: t1_g1d_world, t5_delta_moment_bound,
t5_weighted_unweighted (numbers in §1). Figures: 19 regenerated, 17 differ
from the committed ones only in the footer stamp (pixel rows within the bottom
2.7 %) and were restored; t1_g1d_world.png and t5_delta_moment_bound.png carry
real changes and are kept. `validation/price_defaults.py --plot` (all
families) re-run: the figure differs by its footer only (restored); its
li6-incl table moved A_zz ×1.06576 on every backend row and no ratio
(test_price_defaults re-taken).

## 5. Open points

1. **The default inclusive ⁶Li `helicity-flip` / `apar` runs are refused
   since A1.** At config 1, P_z = P_e = 0.7, the shipped window's top cell
   x = 0.955, Q² = 1897 gives 1 + w_avg = −0.0185 (vector −0.2578, rank-2
   −0.7607; with the pre-A1 factors +0.0396 by the measured ratios, and the
   pre CLI run was accepted). A₁ = g₁/F₁ = 0.616 there, so A3's |g₁| ≤ F₁ clip
   will not lift it; by A4's text ("any refusal that survives A3 (a rank-2
   one) is the right answer and stays") it stays. The CLI and USAGE say so;
   `--x-max 0.95` runs. P2 should re-measure it with the deuteron cases.
2. **The inclusive deuteron reads the ⁶Li rank-2 transfer** (pre-existing:
   its class-level toy b₁ is `Li6B1(MillerB1)`), so A1 moved the deuteron's
   b₁ / A_zz ×1.06576 as well as its vector sector ×0.97973. Not a choice A1
   makes; filed.
3. A caller-supplied kernel cannot report its transfer (a closure); the meta
   key is then absent and the banner says "not knowable".
4. Not relabelled: historical run records quoting 0.811228 / 0.921947 / 0.9325
   (they record their own date), `docs/surveys/needs_survey.md` (a survey of
   the sibling), BENCHMARK_PLAN.md's MSTW g₁ᵈ sub-row figures and
   validation/benchmarks/README.md row 13's letter sub-rows (third-decimal
   moves, 0.715 → 0.714 etc.).
5. A git note: one `git checkout -p` was typed by mistake with stdin closed;
   it discarded nothing (`git diff --stat` unchanged immediately after).
