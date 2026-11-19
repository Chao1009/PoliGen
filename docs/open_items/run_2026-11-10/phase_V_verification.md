<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-10 — phase V (verification) and phase F (fixes)

Written by the F implementer (an Opus worker under the supervising session)
on 2026-11-17, on base commit `931e0fb` (the tree clean but for the
supervisor's dated notes of 2026-11-17 in `DECISIONS.md` and
`../run_2026-11-03/DECISIONS.md`, which F did not touch). V reviewed the
run's implementation of J1–J9; F applied the SUPERVISOR's disposition of
each finding (J1 and J2 are the author's decisions, J3–J9 and the dated
notes and rulings under them the supervisor's under the author's
delegation; nothing here re-opens one). The supervisor's ruling of
2026-11-17 under J6 (findings CODE-1, CODE-4) and correction under J2
(PHYSICS-2) are in `DECISIONS.md`; the findings and verdict files are in
the session scratch only, so this record carries the table.

## 0. What V was, and what F did, in one paragraph

Two read-only Opus lenses reviewed `b3f8858..931e0fb`: lens 1 (physics and
honesty) filed PHYSICS-1 … 6, lens 2 (code and regression) CODE-1 … 7; one
skeptic per lens was told to refute each finding and re-ran the evidence.
Of the 13: 7 confirmed, 4 partly, 2 confirmed in code but never reached or
of no consequence (CODE-5, CODE-7); no BLOCKER survived (PHYSICS-1 was filed
as one and rated MAJOR). F applied 10 (PHYSICS-1 … 6, CODE-1, 3, 5, 6) and
recorded 3 as not applied (CODE-2, CODE-4, CODE-7), each as the supervisor
disposed. The one behaviour change is CODE-1: `PositivityClippedG1`
answers `extrapolation()` for the backend it wraps and the polarised grid
keys are read through it — it moves three `meta` keys on a CALLER's wrapped
object only, and every npz and HepMC3 file of the run list is byte-identical
before and after (§2). Everything else is a same-line dated clause in a
document, a comment or a Python help string, one test rename, and one new
C++ case and one new Python test. No physics default, constant, tolerance,
registry label, verdict rule or reference-JSON number moved; no check was
loosened; `docs/HEPMC3_CONVENTION.md`, the hashed design blocks, the first
160 lines of `t2_hermes_b1_table2.py` and `validation/reference/*.json`
were not touched; every edit of an existing file is inside an existing line
(each changed file is +n/−n, except the two test files, which gained the
new cases at their ends); no runtime string on a default path moved.

## 1. Every finding

| id | lens severity | skeptic: verdict, severity | disposition | what was done |
|---|---|---|---|---|
| PHYSICS-1 | BLOCKER | CONFIRMED (passage (d) partly), MAJOR | APPLIED, widened by the skeptic's list | every present-tense "`--pol-sf` / g₁ is read only under a beam helicity" statement got a same-line `V (2026-11-17), PHYSICS-1; J1 note (2026-11-12)` clause: on the default inclusive deuteron g₁ also bounds the rank-2 slot, so `--pol-sf` is read there under every fill — `docs/USAGE.md` (§2b's two sentences, the "`--pol-sf` does not run on the coherent channel …" section's lead and its "what never happens" sentence, the "What is not refused" paragraph; the 2026-09-05 inclusive-d measurement marked as of its date, superseded for the deuteron), `docs/PHYSICS_CHANNELS.md` (the `--pol-sf` row), `include/lipolgen/pipeline.hpp` (`PipelineConfig::pol_sf`'s doc; the comment above `pol_sf_is_read`: the "FALSE ON `CoherentLi6` AND NOWHERE ELSE" line, the "exactly one product" sentence, the 2026-09-05 measurement), `python/bindings.cpp` (the `PolSfSource`, `pol_sf_is_read` and `PipelineConfig.pol_sf` help strings), `python/lipolgen/__init__.py` (the `POL_SF` comment), `src/core/pipeline.cpp` (the "`plan_has_beam_helicity` is the whole of it" comment). The runtime reason string of `pol_sf_reach_report` ("… nothing else reads a PolSF") is NOT changed: a same-line C++ comment says it is emitted only off the default inclusive deuteron, where it is true. `docs/CONVENTIONS.md`'s label example (true off the deuteron) left as it is (optional) |
| PHYSICS-2 | MAJOR | CONFIRMED (the supervisor's check reproduced), MINOR | APPLIED | re-counted here (§2): 3154 numbers, not 3170 — the 3170 counted 16 NaN-to-NaN leaves (13 `spectator.json`, 3 `coherent.json`) as changed. Dated same-line clauses keeping the old numbers visible: `phase_P1_redump.md` (§0, the table's spectator 28 → 15, coherent 9 → 6 and total rows, the `diff_stats.py` description, "Against the picture"), the P1 row of `STATUS.md`; `../run_2026-11-03/phase_P5_port_gate.md` (the three "3170" places and its two table rows), `../run_2026-11-03/SUMMARY.md` item 2, `../run_2026-11-03/STATUS.md`. The check was the supervisor's (recorded in the V findings); the corrections under J2 and A13 are the supervisor's notes. The commit messages of `1c98f62` and of run 2026-11-03's phase P5 keep "3170": they cannot be changed |
| PHYSICS-3 | MINOR | CONFIRMED, MINOR | APPLIED | `phase_P1_redump.md` ("Where the 9.4e−13 comes from" and §5 item 1): today's value is the correctly rounded arccos of the DOUBLE `1/np.sqrt(3)`; the magic angle itself is 0.82 ULP from it (the old value 1.82) — dated correction (numbers from the lens's and skeptic's mpmath runs, 50 digits; not re-run in F) |
| PHYSICS-4 | MINOR | CONFIRMED, MINOR | APPLIED | `phase_P2b_spin1_domain.md` N3: dated correction in the skeptic's wording — the LO density is exactly 0 there by construction; the kernel's corrections (vector part 0.9091 against 0.9136 at LO, tensor part −0.08628 against −0.0864) leave +0.0046; re-measured in F (§2) |
| PHYSICS-5 | MINOR | PARTLY, MINOR | APPLIED (USAGE only) | `docs/USAGE.md`, the `b1_clipped_rate_fraction` sentence, same line: where it binds b₁ IS the bound 3 (F₁ − \|g₁\|), not the toy model, and the confined b₁/F₁ is not monotone in x (re-measured, §2). The runtime `rank2_input` string and `PHYSICS_CHANNELS.md` (optional) not changed |
| PHYSICS-6 | MINOR | PARTLY + one missed sentence, MINOR | APPLIED as the skeptic said | USAGE refusal paragraph: the bold A3 lead kept, followed on its line by "(J1 note (2026-11-12): today none is refused — see the end of this parenthesis; V (2026-11-17), PHYSICS-6)"; the "with the toy g₁ … configs 0 and 1 are still refused" sentence gained a dated clause quoting today's exits (re-run, §2) |
| CODE-1 | MAJOR | CONFIRMED, MINOR | APPLIED, completed, under the supervisor's ruling of 2026-11-17 under J6 | `PositivityClippedG1::extrapolation()` forwards to the inner backend (appended to the class's `g1n` line in `include/lipolgen/sf.hpp`; the cited `class` line not moved); `pol_sf_grid_report` (`python/bindings.cpp`) casts to `LhapdfG1` through the wrapper. Pinned by the C++ case "V CODE-1: PositivityClippedG1 answers extrapolation() for the backend it wraps" (`tests/test_lhapdf.cpp`) and `test_v_code1_a_wrapped_lhapdf_grid_reports_the_grid_it_wraps` (`python/tests/test_hepmc.py`) |
| CODE-2 | MINOR | PARTLY, NONE (MINOR at most) | NOT APPLIED (the supervisor's disposition) | the shipped reason ("g1 also bounds the rank-2 slot … so it is read under every fill") is literally true — g₁ is evaluated under every fill — and the supervisor's ruling of 2026-11-12 under J1 keeps the configuration-level rule; the proposed clause would move `meta["pol_sf_reach"]` and the knob-row reason on every inclusive deuteron run. Left for the summary (§5) |
| CODE-3 | MINOR | CONFIRMED, MINOR | APPLIED | `test_every_default_deuteron_plan_is_refused_at_x_max_1_by_the_rank2_term` → `test_the_default_deuteron_builds_and_the_unconfined_kernel_is_refused` (`python/tests/test_j1_deuteron_b1.py`, the old name in a dated comment on the same line); no document, record or script cites the old name (grep) |
| CODE-4 | MINOR | PARTLY, NONE | NOT APPLIED, by the supervisor's ruling of 2026-11-17 under J6 | `set_run_attributes({})` writes neither attributes nor the absence line — phase P4's choice, pinned by its test; the convention document's two sentences hold as written and it is not edited again; no shipped path passes an empty map. Listed for the author (§5) |
| CODE-5 | MINOR | CONFIRMED in code, never reached, NONE | APPLIED (comment only) | appended to the existing line of the clamp in `InclusiveKernel::tables` (`src/core/xsec.cpp`): the b₁ == 0 / set b₂ / upper bound below zero corner leaves "b₂ in step" undefined and is not reachable on a shipped kernel (\|g₁\| ≤ F₁ since A3; the deuteron's b₂ is unset) |
| CODE-6 | MINOR | CONFIRMED, NONE | APPLIED | one clause in `phase_P4_provenance.md` §4's OFF-build row: the helper's own exit 1 came from its last line, `tail -1 f1 f2 f3`; the three components are as tabulated (the scratch helper already uses `tail -n 1`) |
| CODE-7 | MINOR | CONFIRMED, NONE | NOT APPLIED (the supervisor's disposition) | `phase_P1_redump.md` quotes the manifest's changed DATA string (the line reference into the sibling's `tagged.py` that the manifest carries, old and new); a quotation, not a citation, and the citation gate passes it |

## 2. What moved in F, and the evidence

**Strings, labels, numbers on a run.** Only CODE-1's, on a caller's
`PositivityClippedG1` wrapping an LHAPDF grid passed as `pol_sf_obj` (6Li,
`helicity_flip_plan(1.0, 0.7, 0.7)`, 20 events):

| key | before (the skeptic's reproduction at `931e0fb`) | after (F, scratch `code1_probe.py`) |
|---|---|---|
| `pol_sf_extrapolation` (npz meta; `lipolgen.pol_sf_extrapolation` through `run_provenance`) | "the backend's own (not an LHAPDF grid)" | "nearest" |
| `pol_sf_grid_q2_min` | NaN | 1.0 |
| `pol_sf_below_grid_frac` | NaN | 0.0 |

After F the wrapped grid reports exactly what the bare grid object and the
default selector report (`nearest`, 1.0, 0.0; `pol_sf` stays `custom` for a
caller's object); a wrapped `ToyG1` reports what a bare `ToyG1` object does
("the backend's own (not an LHAPDF grid)", NaN, NaN); the default
(`nnpdfpol`: `nearest`, 1.0, 0.0) and the toy selector (`none (closed-form
toy)`, 0, 0) are unchanged — all four pinned by the new Python test.

**A/B, the default paths.** The run list `run_1110/p2/scripts/ab_runs.sh`
(22 CLI runs over every channel and the RC / cluster-wave / FSI knobs, the
six C++ examples, the 24 inclusive-d plan × config runs at x_max 0.85 and
the 24 inclusive-⁶Li ones; 300 events, seed 11): pre = `931e0fb`'s tree and
build before the first edit, post = the final tree after the rebuild.
Exits identical (all 0 but `inc6_hul`, 1, the A17 refusal, pre and post);
`ab_compare.py`: **all 69 npz identical** (every array byte for byte, meta
and knob rows) and **all 75 HepMC3 files identical**; `cmp`: 144 of 144
files byte-identical. The wrapper path is not in the list: no selector
attaches a wrapper (the kernel's own A3 clip lives inside the kernel).

**Re-measured for the document clauses (F, this tree; scripts in the
session scratch):**
* PHYSICS-2 (`p2_count.py`, leaf by leaf over `git show b3f8858` /
  `1c98f62`, a NaN pair counted unchanged): xsec 2878, spin 202, tagged 51,
  spectator 15 (13 NaN pairs), coherent 6 (3 NaN pairs), bookkeeping 2,
  beams 0, b1_default_li6 0 — **3154**; the naive `!=` count gives 3170;
  `git diff --stat` 3155 insertions = 3154 one-number lines + the manifest
  line.
* PHYSICS-4 (`n3row.py`, the skeptic's script): smallest density 0.00462 at
  x = 0.9549926, Q² = 733.63 (λ_e +1, m −1 and λ_e −1, m +1); b₁/F₁ =
  0.2592000, g₁/F₁ = 0.9136, LO value 1 − b₁/(3F₁) − \|g₁\|/F₁ = 0.0; kernel
  tensor part −0.0862753 (LO −0.0864), vector part 0.909106 (LO 0.9136).
* PHYSICS-5 (`p5_scan.py`, the default deuteron kernel, config 0, Q² =
  733.6): confined b₁/F₁ = 0.0129 / 0.1165 / 0.3840 / 0.5582 / 0.4730 /
  0.3076 / 0.2592 at x = 0.6 / 0.7 / 0.794 / 0.85 / 0.871 / 0.9 / 0.955,
  equal to 3(1 − \|g₁\|/F₁) from x = 0.85 up; the run banner's "max
  unconfined b1/F1 = 5.422 at x = 0.955, Q2 = 733.6".
* PHYSICS-6: `lipolgen-run --isotope d --channel inclusive --config 0|1
  --plan helicity-flip --events 200 --seed 1234 --pol-sf toy`: exit 0 and
  exit 0; the banner's spin-1 domain line binds on 4.792e−5 % / 3.871e−5 %
  of the rate, 40 / 59 cells, all upper.

## 3. Files

Documents: `docs/USAGE.md`, `docs/PHYSICS_CHANNELS.md`. Records and
boards: this record (new); `phase_P1_redump.md`,
`phase_P2b_spin1_domain.md`, `phase_P4_provenance.md`, `STATUS.md` (the P1
row's clause, the V and F rows); `../run_2026-11-03/phase_P5_port_gate.md`,
`../run_2026-11-03/SUMMARY.md`, `../run_2026-11-03/STATUS.md`. Sources:
`include/lipolgen/sf.hpp` (CODE-1), `include/lipolgen/pipeline.hpp`
(comments), `src/core/pipeline.cpp` (comments), `src/core/xsec.cpp`
(comment), `python/bindings.cpp` (CODE-1's cast; three help strings),
`python/lipolgen/__init__.py` (comment). Tests: `tests/test_lhapdf.cpp`
(+1 case), `python/tests/test_hepmc.py` (+1 test),
`python/tests/test_j1_deuteron_b1.py` (rename). Not touched by F:
`DECISIONS.md` of either run (the supervisor's), `docs/HEPMC3_CONVENTION.md`,
`validation/reference/*.json`, any harness, the report.

## 4. Gates (final tree)

All run by F on the final tree, in the foreground, on the shared machine
(load 1.6 at the start, 5–12 during the pytest parts).

| gate | result | wall |
|---|---|---|
| build (`cmake --build build -j2`, after the header change; one compile error in the new test case fixed, then the incremental rebuild) | 0 warnings in both logs | 283 s + 78 s |
| `lipolgen_tests` (full) | **441/441 cases, 17 244 832 assertions** (440 / 17 244 828 before: + the CODE-1 case, 4 assertions) | ≈ 15 min (13:01:47 to ≈ 13:17) |
| pytest, part A | 359 passed, 4 skipped | 325 s |
| pytest, part B | 138 passed, 1 skipped (first run: 137 passed, 1 skipped, **1 failed** — `test_bench_t4_t6.py::test_t6_measured_row_is_pinned`, its `runtime_s < 30.0` ceiling at 35.08 s under load ~5–10; alone at load 6.4: passed, 26.97 s; part B re-run at load 5.0: as tabulated) | 180 s; 123 s |
| pytest, part C | 188 passed | 159 s |
| pytest, part D | 249 passed, 90 skipped | 814 s |
| pytest, part E | 277 passed, 61 skipped | 857 s |
| pytest, part F | 324 passed | 308 s |
| pytest, sum | **1535 passed, 156 skipped, 0 failed** (1534 / 156 before: + `test_v_code1_a_wrapped_lhapdf_grid_reports_the_grid_it_wraps`) | |
| bench subset (`LIPOLGEN_BENCH=1 … -k bench`) | **372 passed, 1 skipped** | 460 s |
| report `--all --plots --compare` (baseline: `931e0fb`'s `REPORT.json`, copied before the first edit) | 15 pass / 8 fail / 0 blocked; **0 changes**; `--check` OK | 185 s |
| report files and figures | `REPORT.json` / `REPORT.md` differed by date, revision and runtimes only, 19 PNGs by the footer rows only (pixel diff: first differing row at ≥ 0.972 of the height) — all restored | |
| A/B (§2) | 69 npz + 75 HepMC3 byte-identical; exits identical | 85 s + 120 s |
| citations | strict 1243 references / 95 ranges / 7 external / 0 broken / 24 allow-listed; `--audit-ranges` 209, 0 REFUSED — no cited range block re-recorded; `--records-only` 1432 citations in 99 records (this record added), 0 broken | |
| SPDX; doc-link gate | 136 files, all carry the line; `test_doc_link_gate.py` 112 passed (31 s); `git diff --check` clean — these three re-run after this record and the board rows were written | |

## 5. Open points (for the summary)

1. **The empty-map corner (CODE-4).** `set_run_attributes({})` gives a file
   with no `lipolgen.*` attribute and no absence line — indistinguishable
   from a pre-A10 file. It stands by the supervisor's ruling of 2026-11-17
   under J6; no shipped path passes an empty map. For the author.
2. **A run-level `pol_sf` rule (CODE-2).** On a default inclusive deuteron
   run below the binding rows (x_max ≲ 0.87) g₁ is evaluated and moves
   nothing, yet `pol_sf` is "read" (the configuration-level rule the
   ruling of 2026-11-12 keeps); a "read where it binds" rule would move
   `meta["pol_sf_reach"]` and the knob row on those runs. A possible
   refinement, not taken.
3. **Found on the way.** The C++ case run as
   `./build/lipolgen_tests "<name>"` runs the WHOLE suite (doctest ignores a
   positional name; `-tc=` selects one) — RULES.md's "a quoted name runs one
   case" is not so on this build; F's full C++ gate is that run (§4).
   `pol_sf_grid_report` and `PositivityClippedG1::extrapolation` unwrap one
   level; a wrapper of a wrapper is not a shipped object.
