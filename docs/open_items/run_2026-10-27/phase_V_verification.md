<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-27 — phases V and F: the adversarial verification of the whole run, and its fixes

Branch `claude/run-2026-10-26`, over `7fd441c` (D4 `9779edd`, D5 `cb75ec7`,
D6 `c52c164`, D7 `7fd441c`, on top of `ca9b23a`). 2026-11-03.

**V.** Four adversarial verifiers read the whole run, each through one lens —
honesty, physics, code, regression. Their reports live in the session
scratchpad (`final_verify/{honesty,physics,code,regression}.md`), outside the
repository. Verdicts: honesty "the run holds up on substance; fixes needed in
documents"; physics "the run holds on physics; no decision, combination or
re-derived number refuted; all findings are wording or scope"; code "mostly
holds; no gate refuted; one medium regression"; regression "could not refute
the run; every gate reproduces the records' final tallies exactly". In all:
0 high; 5 medium (honesty H1–H3, physics P1, code C1) and one
low-medium each in honesty (H4) and code (C2); the rest low or info.

**F.** One fixer re-checked every finding against the tree before touching
anything, applied the minimal fix, and ran every gate on the final tree. No
finding was refuted. Nothing was committed by this stage (the supervisor
commits), and nothing was published or sent. The sibling repository was read
through shell tools only, and not edited.

Rules kept: no physics constant, tolerance, registry label, reference-JSON
number, vendored file or hashed design block moved; no census test loosened;
the two re-pins (`tests/test_pipeline.cpp`'s D2 case, see code 1) carry a
dated comment naming the decision; edits to cited documents are line-neutral,
or are appended where nothing below is cited. Neither `--fix` nor
`--record-ranges` was needed. In `DECISIONS.md` the author's quoted words are
verbatim and the bold decision texts keep their meaning. The D7 bold text was
reworded as the brief directed ("grid sets … (free nucleon, without its
EPPS21 ratio)").

## 1. Honesty lens

| # | finding | severity | re-check (2026-11-03) | done |
|---|---|---|---|---|
| H1 | `OPEN_ITEMS_SOLUTIONS.md` §14: ct18nlo rows 472571.9 / ×0.79847 / A_zz ×1.25239 / A_∥ ×0.82519 are pre-D7 continued numbers, unlabelled; the grid clause says, in the present tense, that LHAPDF continues below the grid (36.18 %) | medium | confirmed. Re-measured on the frozen tree: inclusive σ 496334.156 (×0.83862), A_zz −6.192590e−4 (×1.19243), A_∥ −9.848331e−4 (×0.84050); tagged-6Li-α, VMC wave: A_∥ ×0.84049; F₂ᵖ(3e−4, 0.7) toy / frozen CT18NLO = 1.2516 | both rows labelled "(continued, pre-D7; frozen since D7: …)"; the at-0.95 number labelled; the grid clause put in the past tense, with the frozen 39.23 % and ×1.2516 beside it (line-neutral) |
| H2 | I7 §1 quotes the `fbc24f2` D7 entry as "verbatim"; it says the author chose the freeze | medium | confirmed: the author's only D7 words are "Your recommendation sounds the best approach, go with it"; `plan_D7` §1.3 proposed continuation; its check proposed `nearest` | §1 re-quotes the current entry and says so; the freeze is attributed to the supervisor (design continuation, check nearest); Q² only for the letter's Q₀ freeze. `include/lipolgen/lhapdf_sf.hpp`'s "the author's choice (D7)" → "the supervisor's choice under D7" (same lines) |
| H3 | no consolidated "for the author" list; STATUS's V row points to a "Verification" section that does not exist | medium | confirmed | `STATUS.md` gains `## For the author`, every open item with its record section; the V row points here |
| H4 | `DECISIONS.md` D4 overstates: "every cos 2φ amplitude halves"; `moment_A/B` "only after the sibling changes"; I4 §8 item 8's (1)–(6) never applied | low-medium | confirmed against I4 §8 item 8 and the commit (`tensor_gamma`'s h₂ unchanged) | (1)–(6) applied in the italic annotation: the general c_eff(m) formula; CTP#1712 Eq. (16); "every **Δ-induced** amplitude", θ_S = 0 bit-identical, leakage's Δ-equivalent doubles; `moment_A/B` the bag model at face value in this tree from D4 on, the letter encoding 2 × C_BAG until the sibling changes; the HJM attribution and the b-sector reading; the port gate reading `xsec.json` |
| H5 | `DECISIONS.md`'s "5 %" and STATUS / I7's "opposite lhapdf.conf … 4.5e−4" are different comparisons, neither says so; the test's "opposite" conf keeps `continuation` | low | confirmed: this store's conf is `logcubic` / `continuation` / ForcePositive 0; `_OPPOSITE_CONF` differs in interpolator and ForcePositive only; +5.03 % is `plan_D7`'s continuation → nearest | each comparison named where quoted (DECISIONS D7, I7 §2.3, the STATUS I7 row); a conf variant with a DIFFERENT extrapolator added as a test (`test_the_default_numbers_do_not_depend_on_the_conf_extrapolator`, `Extrapolator: error`, `linear`, ForcePositive 2: identical, 496334.155954048) |
| H6 | USAGE §8's D7 throughput "idle machine, median of 3, 2026-10-30" was taken inside F7's gate runs, raw log lost | low | confirmed by I7 §16 V7-6 | labelled "taken while F7's gate runs shared the machine, the raw log not kept and not re-measured"; listed for the author |
| H7 | D6 entry omits that VMC is itself ~1.9 σ from the measurement, that "well accepted" was not surveyed, that CD-Bonn is the deuteron-only alternative | low | confirmed: (−0.048 + 0.025)/0.012 = −1.92; `plan_D6` says usage was not surveyed and CD-Bonn exists as a deuteron only | one clause added to the D6 annotation |
| H8 | the corrected STATUS row is an uncommitted edit | low | true | not a tree fix; the supervisor commits this stage's tree with the V / F rows |
| H9 | D7 bold "grid backends the letter's production uses" vs "neither the letter's grid chain"; the letter freezes Q² only, on a nuclear set | info | confirmed (sibling `structure.py`, read with grep) | bold text reworded as directed; the Q₀ freeze marked Q²-only in DECISIONS, I7 §1, USAGE (two sites), `pipeline.hpp`, `lhapdf_sf.hpp` (line-neutral) |

## 2. Physics lens

| # | finding | severity | re-check | done |
|---|---|---|---|---|
| P1 | the positivity violation is described too narrowly: |A₁| ≤ 1 is broken from x ≳ 0.83 (n) and ≳ 0.87 (d, p); `--x-max 0.95` cures the refusal, not the violation | medium | re-taken on `LhapdfSF("CT18NLO")` / `LhapdfG1("NNPDFpol11_100")`: A₁ⁿ 1.046 / 1.222 / 1.455 / 1.914 at x = 0.83 / 0.85 / 0.871 / 0.90 (Q² = 10); A₁ᵖ 1.006 at 0.871 (Q² = 10), 0.958 at 0.85; matches | stated in I7 §1 (the development item) and §3.5, in USAGE's refusal paragraph and §10, and filed under the LO development item in `DECISIONS.md` D7. `plan_D7` (a design record) is not edited |
| P2 | the T1 headlines say "on Q2 >= 1" where the cut is Q² ≥ 1.677 | low | confirmed | `t1_nmc_f2d`, `t1_nmc_f2d_over_f2p`, `t1_hermes_f2d` format the cut from `h["q2_floor"]` (`%.4g`); the report regenerated (§6) |
| P3 | I7 §6 quotes the robustness number of the continued backend | low | the verifier's frozen 244.41/211, p = 0.0571 taken as measured | I7 §6 quotes 244.41/211, p = 0.057 (frozen) and adds that the declared verdict does not cover the frozen region (39 % of the default rate) |
| P4 | "the declared rule names the shipped default" read two opposite ways in I6 and I7 | low | confirmed | one item in I6 §9 (item 12) and one sentence in I7 §6 naming both readings, the choice left to the author |
| P5 | `PHYSICS_CHANNELS.md`'s `--unpol-sf` row: continued numbers unlabelled; "the kernel's default ToyG1" stale | low | confirmed | labelled "2026-09-04, continued, pre-D7", frozen values beside (496334.2, ×0.83862; A_zz ×1.19243; A_∥ ×0.84050); the ToyG1 sentence qualified as `--pol-sf toy` (in-line) |
| P6 | `--isotope d --plan helicity-flip` and `--channel tagged-d-p --plan helicity-flip` exit 1 at the default | info | documented in I7 §3.5 | listed for the author |
| P7 | I7 §5's "5000-event CLI files" row has no recipe | info | re-taken with the CLI default seed: t-peak 1.027943 (max 3.7617), polrad-full 1.037955 (max 4.6656), the verifier's numbers | I7 §5's row gives the recipe and the re-take beside the unreproduced 1.027211 / 1.037150 |

## 3. Code lens

| # | finding | severity | re-check | done |
|---|---|---|---|---|
| C1 | `validate()` accepts any object in `unpol_sf_obj` / `pol_sf_obj` under a grid selector: `ToyF2` under `ct18nlo` runs 591846.17 pb labelled ct18nlo | medium | reproduced in Python before the fix (validates, the toy's σ, knob row "ct18nlo") | the LHAPDF tier registers the ONE cached object it builds (`register_shipped_grid_sf`, called in `make_lhapdf_unpol_sf` / `make_lhapdf_pol_sf`); `validate()` calls `refuse_foreign_grid_sf_objects` (pipeline.hpp D7 block, defined at the end of pipeline.cpp; the call takes a line freed from clause (5)'s comment, so validate() is line-neutral). A grid selector holding any object but the shipped one is refused and names `Custom`; the tier's object (also the one `attach_shipped_sf_backends` puts there) is accepted. Tests: C++ "D7 verification (2026-11-03, code finding 1)" in `tests/test_lhapdf.cpp` (a fresh `LhapdfSF("CT18NLO", 0)` copy is refused as well); `tests/test_pipeline.cpp`'s D2 case RE-PINNED (dated comment): ScaledF2 under Ct18Nlo and ScaledG1 under NnpdfPol now throw. Python `test_a_grid_selector_refuses_any_object_but_the_shipped_one` (3 channels: refused; as `Custom` it runs the toy's σ and the knob row says custom; the shipped objects run 496334.155954048). USAGE §10 says so |
| C2 | the `pol_sf_reach` reason says ×0.705841 (continued) | low-med | confirmed; frozen 9424.476321 / 12476.025113 = 0.755407 | `pipeline.cpp`'s reason and `PHYSICS_CHANNELS.md`'s coherent row say ×0.755407 (×0.705841 labelled continued); the bindings' below-grid docstring, which said LHAPDF "extrapolates downward rather than freezing" with the continued 0.447460 / 0.361791, now says frozen, 0.4837146 / 0.3923453 (continued beside) |
| C3 | toy-era numbers unlabelled: the rc_tail triples (USAGE, `rc.hpp`, `--rc-tail-model` help), `--r-source`'s −0.6847 %, "r_func reaches … ToyG1" | low | confirmed; also two more copies of the triple in USAGE's `--rc-tail-model` table row; −0.6847 % (toy) / **−0.59899 %** (shipped, `--b1-unpol ct18nlo`, x_max 0.95) re-measured | labelled "toy SFs" with D7's grids beside wherever measured (event-weighted triples at P_z 0.7 and 0); the cell-weighted triple and K/D_φ / cos 2φ percentages labelled toy only (not re-measured); the ToyG1 clause qualified `--pol-sf toy`; I7 §5's row corrected (it named those sites as labelled) |
| C4 | an explicit default on an unread knob: D6 refuses, D5 accepts and records `refused`, D7 accepts and labels | low | reproduced all three | STATED, not unified, in `docs/CONVENTIONS.md` ("Refuse or label", line-neutral); pinned by `test_an_explicit_default_on_an_unread_knob_is_handled_three_ways`; for the author |
| C5 | the `b1_unpol ct18nlo -> mstw` cell passes vacuously in the PYTHIA-off job | low | confirmed (`make_config` refuses mstw without the tier → "refused" branch) | `V(…, pythia=True)`; the cell skips with `NO_PYTHIA_CELL`'s reason when `HAVE_PYTHIA8` is false |
| C6 | no tail sub-knob cell runs on the shipped `polrad-full` | low | confirmed (RC_TPK) | two cells per full spec on the default tail: `rc_qe_suppression 1 -> 0.5` and `rc_fq_scale 1 -> 2` (14 cells: 8 pass, 6 skip because the base is refused on 7Li helicity-flip, d inclusive and 7Li tagged, as for every RC cell); the matrix is 641 cells (CONVENTIONS and USAGE counts updated in-line) |
| C7 | retyped defaults: D5 has no static_assert; `price_defaults`' "toy" typed literally; `escape_for` keyed on the set name | low | `RcOptions` holds a `shared_ptr`, so no `static_assert` | `test_the_texts_that_call_polrad_full_the_default_are_tied_to_it` (`test_rc.py`): fails when `RcOptions::tail_model` moves and names every retyped site; `price_defaults.py` asserts the defaults are not `toy` on the line that reads them; `escape_for` derives "shipped" from `UNPOL_SF_DEFAULT` / `POL_SF_DEFAULT`. The examples' default-enum factories throw loudly on a move, so they are left as they are |
| C8 | a store with sets but no `lhapdf.conf` gives `error`, not `blocked (environment)` | low | confirmed | `_report.py` also matches the pre-flight's "lhapdf.conf was not found" (reason "… present, lhapdf.conf was not found"); the five T1 harnesses classify it likewise (line-neutral); `test_a_store_without_lhapdf_conf_blocks_the_row_by_the_environment` (4 harnesses, exit 2, classified) |
| C9 | a kernel beside an explicit `unpol_sf='ct18nlo'` is accepted unread | info | matches rule (3) | none (consistent with the documented rule) |

## 4. Regression lens

| # | finding | severity | re-check | done |
|---|---|---|---|---|
| R1 | the missing-set message offers `--unpol-sf toy` to the `generate_*` examples, which have no such option | low | confirmed | the escape text adds "(the generate_* example programs have no toy switch: use lipolgen-run or lipolgen::toy_config())" |
| R2 | I7 §16 "Nothing is committed." is stale | low | true | "Committed as `7fd441c`" |
| R3 | CI `DEPS_CACHE_KEY` still `-v1` after `build_deps.sh` changed | info | true | bumped to `-v2` |
| R4 | `ci.yml`'s "9-11 min" extrapolation is stale | info | true (C++ ~8 min, pytest ~13 min here) | labelled STALE with the 2026-11-03 timings, not re-extrapolated |
| R5 | the LHAPDF-off build after D7 was never recorded | info | the verifier's measurement (0 warnings, examples not built, 411 / 411, 17 242 996 assertions), not repeated | noted in I7's closing paragraph |

## 5. Files changed by F

Code: `include/lipolgen/pipeline.hpp` (the grid-object API, appended; one
comment line-neutral), `src/core/pipeline.cpp` (`validate()` line-neutral;
the registry and the refusal appended; the reach string), `src/lhapdf/shipped_sf.cpp`
(registration), `src/lhapdf/lhapdf_sf.cpp` (`escape_for`), `include/lipolgen/rc.hpp`,
`include/lipolgen/lhapdf_sf.hpp` (comments), `python/bindings.cpp` (a
docstring), `python/lipolgen/{cli.py, __init__.py}` (help texts, comments).
Tests: `tests/test_lhapdf.cpp`, `tests/test_pipeline.cpp`,
`python/tests/{test_sf_backend, test_knob_provenance, test_rc,
test_bench_t1_nucleon, test_d7_shipped_kernel}.py`. Validation:
`validation/benchmarks/{_report, t1_nmc_f2d, t1_nmc_f2d_over_f2p,
t1_hermes_f2d, t1_g1d_world, t1_r_world}.py`, `validation/price_defaults.py`;
the regenerated `docs/benchmarking/REPORT.{json,md}`. Documents:
`docs/{USAGE, PHYSICS_CHANNELS, OPEN_ITEMS_SOLUTIONS, CONVENTIONS}.md`; this
run's `DECISIONS.md`, `STATUS.md`, `phase_I6_cluster_wave.md`,
`phase_I7_backends.md`, and this record. CI: `.github/workflows/ci.yml`.

## 6. The report

`make_report.py --all --plots --compare <REPORT.json of 7fd441c>`, run once on
the fixed tree (90 s wall): **14 pass, 9 fail, 0 blocked, 0 blocked
(environment), 0 error, 0 skipped (opt-in)** — 7fd441c's summary. `--compare`:

```
MOVED    t1_hermes_f2d:   … p = 0 on Q2 >= 1, W2 >= 4 …   ->  … p = 0 on Q2 >= 1.677, W2 >= 4 …
MOVED    t1_nmc_f2d:      … p = 0 on Q2 >= 1, W2 >= 4 …   ->  … p = 0 on Q2 >= 1.677, W2 >= 4 …
MOVED    t1_nmc_f2d_over_f2p: … p = 0.0904 on Q2 >= 1, …  ->  … p = 0.0904 on Q2 >= 1.677, …
compare: 3 changes (the M2 rule: any change exits 1)
```

Every one of the three is physics finding P2's label: the χ², ndf, p and
status of each row are unchanged, and only the printed cut moved to the one
the numbers were computed on. Nothing else moved; the runtimes are volatile.
All 19 regenerated figures differed from the committed ones only in the
footer stamp (pixel rows 417–632 by figure, the bottom strip only) and were
restored, as I7 did. `REPORT.json` records `git_rev` `7fd441c` with
`git_dirty` true: it was regenerated before this stage's commit.

## 7. Gates on the final tree (2026-11-03, one at a time)

| gate | before (7fd441c, the records' tallies) | after F (this tree) | wall |
|---|---|---|---|
| `cmake --build build -j2` | 0 warnings | nothing to rebuild after the last build (two rebuilds during F: 0 warnings) | 0.3 s |
| `./build/lipolgen_tests` | 427 / 427, 17 244 377 assertions | **428 / 428, 17 244 386 assertions**, 0 failed, 0 skipped (+1 case, the code-1 refusal, 9 assertions; the D2 re-pin is count-neutral) | 417.1 s |
| full pytest | 1497 passed / 149 skipped / 0 failed | **1515 passed / 155 skipped / 0 failed**, 3 warnings (matplotlib Axes3D, pre-existing) — +18 passed: `test_sf_backend.py` 3, the knob matrix's 14 shipped-tail cells (8 pass, 6 skip: base refused) and the three-ways test, `test_rc.py` 1, `test_bench_t1_nucleon.py` 4, `test_d7_shipped_kernel.py` 1 | 784.9 s |
| bench subset (`LIPOLGEN_BENCH=1 -k bench`) | 364 / 1 / 0 | **368 passed / 1 skipped / 0 failed** (+4: the lhapdf.conf-missing harness tests) | 345.5 s |
| strict citation gate | 1241 / 95 / 7 / 0 / 24; SPIN32 19 / 115 / 0 / 0 / 10; PYTHIA_BRIDGE 0 / 0 / 8 / 0 / 0 | **the same** | 2.4 s |
| `--audit-ranges` (R5) | 209 on 179 blocks, 0 REFUSED | **209 on 179 blocks, 0 REFUSED** | 2.5 s |
| `--records-only` (R6) | 1432 citations in 78 records, 0 broken | **1432 citations in 79 records** (this one cites none by line): 922 / 61 / 121 / 129 / 199, **0 broken** | 0.3 s |
| SPDX | 136 files | **136 files** (Markdown records are not in its scope) | 0.04 s |
| `make_report.py --check` | OK | **OK** | 77.5 s |
| `make_report.py --all --plots --compare` | 14 / 9 | **14 pass / 9 fail**, 3 MOVED strings (§6) | 90 s |

The citation gates, SPDX and `--check` were re-run after this record's
tallies were written in (the only edit after the long gates): unchanged.

## 8. Left open

Everything in `STATUS.md` "For the author". From this stage in particular:
the knob asymmetry (C4) is stated, not unified; the K/D_φ / cos 2φ
percentages of `--r-source` and the cell-weighted rc_tail triple are labelled
toy and not re-measured on D7; the 1.027211 / 1.037150 5000-event pair
(P7) remains unreproduced beside its re-take; USAGE §8's D7 throughput has no
log (H6); `plan_D7`'s statements of the top-cell positivity (P1) are left as
the design wrote them.
