<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-03 — phase V (verification) and phase F (fixes)

Two adversarial lenses read the implemented tree at `5372f2d`; every finding
then had one or two adversarial verdicts (reproduction; rules); phase F
applied or refuted each under the run's rules and the supervisor's addendum
of 2026-11-10 to the F brief. The decisions are the supervisor's under the
author's delegation (`DECISIONS.md`); nothing here re-opens one.

## 0. The two lenses and the verification, in one paragraph

Lens 1 (physics and honesty) filed 11 findings — 4 MAJOR (PHYSICS-1 … 4),
7 MINOR (PHYSICS-5 … 11); lens 2 (code and regression) filed 8 — 4 MAJOR
(CODE-1 … 4), 4 MINOR (CODE-5 … 8); no BLOCKER in either. The verifiers gave
27 verdicts on the 19 (21 confirmed, 5 partly, 1 refuted): 18 findings
confirmed or partly confirmed by every verifier that read them, PHYSICS-2
disputed (its reproduction verdict refuted, its rules verdict partly). Every verifier
that rated severity rated the defects MINOR, except PHYSICS-3's
reproduction verifier (MAJOR). Phase F: 15 APPLIED by F (PHYSICS-1, 3, 4, 6,
8, 9, 10, 11; CODE-1 … 7, CODE-2 together with PHYSICS-1), PHYSICS-2 APPLIED
IN PART (the code-comment sentence) and REFUTED IN PART (the CONVENTIONS
paragraph, already marked as history), PHYSICS-5 and PHYSICS-7 APPLIED BY
THE SUPERVISOR (the dated notes of 2026-11-10 in `DECISIONS.md`; F edited
nothing for them), CODE-8 RECORDED AS ADEQUATE. No physics default,
constant, tolerance, registry label, verdict rule or reference-JSON number
moved; no check loosened (one check added, PHYSICS-8); no hashed block,
`docs/HEPMC3_CONVENTION.md`, `docs/theory/SPIN32_FINITE_GAMMA.md` or
`validation/reference/*.json` touched; every edit inside an existing line
(each changed file is +n/−n); two cited range blocks edited in place and
re-recorded (§2). No harness, kernel or report input changed, so the
report was not regenerated.

## 1. Every finding

**PHYSICS-1 (MAJOR) / CODE-2 (MAJOR) — APPLIED.** Claim: `docs/CONVENTIONS.md`
(the C5.5b sentence of the cluster-polarization bullet) and
`docs/PHYSICS_CHANNELS.md` (the D6/C5.5 dilution bullet of §3) still say an
explicit `cluster_wave` is refused off the tags (and PHYSICS_CHANNELS that
`cluster_vmc_mc_sigma` "was already refused the same way"), which A17
reversed. Verdicts: PHYSICS-1 reproduction confirmed (reproduced the
exit-0 run and the not-read meta; would rate MINOR), rules confirmed (stale,
no record chose it); CODE-2 reproduction and rules confirmed. Re-measured
here: an inclusive `make_config` Pipeline's knob rows read `cluster_wave`
NotRead "not read on channel inclusive" and `cluster_vmc_mc_sigma` NotRead
at 0 (scratch `probe1.py`). Changed: one dated clause on the same line at
each site, `V-PHYSICS-1 (2026-11-10), A17 (2026-11-03)` — only a moved value
(`hulthen`) is refused, an explicit `vmc` is accepted and labelled
`not-read`; at the PHYSICS_CHANNELS site also `cluster_vmc_mc_sigma` at 0
reads `not-read` since the supervisor's note of 2026-11-06. CODE-2's second
half (the "+11.61 %" present tense in CONVENTIONS) is PHYSICS-2's and is
judged there. Note: the CONVENTIONS sentence sits inside the bullet whose
head says "the rest of this bullet is the PRE-A1 record"; that label covers
A1's numbers, not the A17 rule, so the clause was added as ruled.

**PHYSICS-2 (MAJOR) — APPLIED IN PART / REFUTED IN PART.** Claim: the
CONVENTIONS C5.5 paragraph states the pre-A1 constants and band (+11.61 %,
0.921947, "the shipped 0.811228", 0.81–0.91) in the present tense; the same
in `include/lipolgen/beams.hpp`'s C5.5 block and `src/core/beams.cpp`
`LI6()`. Verdicts: reproduction refuted (the numbers are right but the
bullet's head, "since A1 (2026-11-03) … the single copies are … The rest of
this bullet is the PRE-A1 record", already marks them history; the code
comments are followed by A1 sentences); rules partly (the CONVENTIONS part
holds as an unaccepted inconsistency, the code-comment part is weak). Judged
sentence by sentence, as the addendum asks: (a) CONVENTIONS — the sentences
sit under the bullet's explicit PRE-A1-record label, with the current values
at the head of the same bullet: already history, left as is (REFUTED for
this part). (b) `src/core/beams.cpp` `LI6()` — the paragraph is followed by
"A1 (2026-11-03): the paragraph above is the PRE-A1 reading": already
history, left (REFUTED). (c) `include/lipolgen/beams.hpp` — "while the
shipped 0.811228 is 4.3 % BELOW it" was a present-tense sentence with no
label of its own: APPLIED, now "the pre-A1 shipped 0.811228 (V-PHYSICS-2,
2026-11-10: the inclusive default until A1, now
`LI6_CLUSTER_POLARIZATION_HULTHEN`)".

**PHYSICS-3 (MAJOR) — APPLIED, as the supervisor's clarification of
2026-11-10 under A5 rules.** Claim: README row 3 calls the HERMES b₁ headline
(χ²(b₁)/6 = 5.26) "the shipped Miller", but a default inclusive deuteron run
carries 0.32753 × MillerB1, on which χ²/6 = 14.56. Verdicts: reproduction
confirmed (0.327525 × MillerB1 at all three points; rated MAJOR); rules
partly (the number is right, the records cover it in P3 §5 item 1 and the
harness's A5 NOTE; the README does not). Changed in
`validation/benchmarks/README.md`: row 3 and convention 8 each gain a dated
clause (`V-PHYSICS-3 (2026-11-10)`, the clarification named): row 3 is
convention 8's opt-in clause — the b₁ model headlined as `MillerB1` in the
form `--b1-model miller` names; a default deuteron run carries
0.32753 × MillerB1 (the ⁶Li α–d transfer, P3 §5 item 1), χ²(b₁)/6 = 14.56
(2/6) recorded; FAIL either way. No status, no harness number moved; the
mapping stays the author's item.

**PHYSICS-4 (MAJOR) — APPLIED, as the same clarification rules.** Claim: row
4 (NMC ⁶Li/d) and the T3 rows keep an EPPS21 headline no shipped run applies,
unnamed in convention 8; P3 §1.2 calls `Epps21Ratio` "the shipped EMC
baseline". Verdicts: reproduction partly (the evidence reproduces;
`EMC_BASELINE_DEFAULT` feeds the polarized-EMC transfer only); rules partly
(the facts hold, "unnamed exception" overreaches). Changed: convention 8's
clause names row 4 and every T3 row on `Epps21Ratio` as the opt-in clause's
case, not an exception (the inclusive rate applies no F₂ ratio; NMC on the
ratio 1 gives 7.815/4 and 24.858/15, recorded — numbers from the P3 scratch
`a5_survey.json`, read here); P3 §1.2's survey row gains a dated
`V-PHYSICS-4 (2026-11-10)` correction on the same table line. No headline,
no status moved.

**PHYSICS-5 (MINOR) — APPLIED BY THE SUPERVISOR.** Claim: the A13 note
quotes 9.4e−13 in `xsec.json` without the record's "among values > 1e−6"
qualifier. Verdict: reproduction confirmed. The supervisor's "Qualifier,
2026-11-10" under A13 is in `DECISIONS.md`; F edited nothing.

**PHYSICS-6 (MINOR) — APPLIED.** Claim: P5 §5 item 3 says the A13 note reads
"at 21:04:39 UTC" and was left unedited; the committed note never carried
that wording. Verdict: reproduction confirmed. Changed: a dated
`V-PHYSICS-6 (2026-11-10)` sentence appended on the item's last line of
`phase_P5_port_gate.md` (resolved before commit; the note carries the local
time and its UTC equivalent; item closed).

**PHYSICS-7 (MINOR) — APPLIED BY THE SUPERVISOR.** Claim: the A1 correction
note says pins were re-taken when the x_max 0.95 helpers were removed; P2
says none moved. Verdict: reproduction confirmed. The supervisor's
"Clarification, 2026-11-10" under A1 is in `DECISIONS.md`; F edited nothing.

**PHYSICS-8 (MINOR) — APPLIED.** Claim: the A11 refusal prints 3.0000001 as
"3" (`%g`). Verdict: reproduction confirmed. Changed: in
`PipelineConfig::validate()` (`src/core/pipeline.cpp`) the value is printed
with the shortest `%.*g` (6 … 17 digits) that round-trips, on the same line,
dated `V-PHYSICS-8 (2026-11-10)`; the text is otherwise unchanged and the
bound `fabs(sigma) > 3.0` untouched. Measured after the build: 3.0000001 →
"cluster_vmc_mc_sigma = 3.0000001 is outside …", −3.1 → "−3.1", −4 → "−4",
1000 → "1000", 3.0 builds. Pin added (a check, nothing loosened): the A11
case of `tests/test_pipeline.cpp` now also asserts the message contains
"cluster_vmc_mc_sigma = 3.0000001 is outside" (same line as the loop's close;
the case reads 12 assertions, passing).

**PHYSICS-9 (MINOR) — APPLIED.** Claim: README row 10 quotes the BONuS AV18
headline without saying its weight table misses two declared CHECKS.
Verdict: reproduction confirmed. Changed: row 10 gains a dated
`V-PHYSICS-9 (2026-11-10)` clause — convergence |dv| 5.5e−4 / 8.3e−4 > 5e−4
at the last two p_s bins and bin averaging −3.2 % > 2 % at p_s = 0.078,
recorded, no threshold moved, no block's verdict flips (numbers from the
harness docstring and P3 §5 item 2).

**PHYSICS-10 (MINOR) — APPLIED.** Claim: about 54 re-pin comments added by
the run name the entry without the date. Verdict: reproduction confirmed (54).
F's scan (every added comment line of `git diff 70982c9 5372f2d -- python/tests
tests` that names an entry A1–A17 and lacks "2026-11-03", excluding the
asymmetry A₁ and lines whose contiguous comment block already carries the
dated form of the same entry) found 75; 71 dated, 4 left because their block
is dated on an adjacent line. Each edit inserts " (2026-11-03)" after the
entry token on its own line; checked mechanically that every changed line
differs from its original only by that insertion. Files (lines): Python —
`test_bench_t5_delta_bound.py` 19, `test_bench_t5_closures.py` 8,
`test_bench_bonus.py` 7, `test_bench_t1_nucleon.py` 6, `test_module.py` 6,
`test_bench_static.py` 5, `test_bench_spin_deuteron.py` 2,
`test_kernel_reference.py` 2, `test_b1_model.py` 1, `test_bench_deeps.py` 1,
`test_price_defaults.py` 1, `test_rc_low_x_anchor.py` 1, `test_sf_backend.py`
1; C++ — `test_tagged.cpp` 6, `test_beams.cpp` 3, `test_b1_nuclear.cpp` 1,
`test_sf.cpp` 1. No tolerance on those lines changed.

**PHYSICS-11 (MINOR) — APPLIED.** Claim: the P4b board row still reads
"(next commit)" (it is `b591880`), and registry row 10's title still says the
constant "stays 0.811228". Verdict: reproduction confirmed. Changed: the
board's P4b row → `b591880` and P6's → `5372f2d` (the run's back-fill of the
previous phase, which the finding's fix names too); registry row 10 of
`docs/open_items/run_2026-09-03/STATUS.md` gains a dated prefix on its title
cell (closed 2026-11-03 by A1, 0.848 / 0.982576; the title is the pre-A1 one).

**CODE-1 (MAJOR) — APPLIED.** Claim: USAGE (the b₁ provenance paragraph of
§2a, the refusal-table row on a kernel beside an off-default selector, and
the §10 configuration rules) still gives the label "caller-supplied kernel";
since A17 it reads "not read with a caller-supplied kernel". Verdicts:
reproduction and rules confirmed. Re-measured here: a caller kernel at the
Hulthén constants gives `run_provenance` `unpol_sf`, `pol_sf`, `b1_model` =
"not read with a caller-supplied kernel". Changed: one dated
`V-CODE-1 (2026-11-10)`, `A17 (2026-11-03)` clause at each of the three
sentences, same line.

**CODE-3 (MAJOR) — APPLIED.** Claim: USAGE §10 quotes the b₁ band ×2 refusal
at D7's 1 + w_avg = −0.008234. Verdicts: reproduction confirmed, rules partly
(stale number). Measured here, `lipolgen-run --b1-model cdks --b1-band-scale 2
--x-max 0.95|0.9|0.85 --events 100`: exit 1 / 1 / 0, the refused cell
x = 0.871, Q² = 150.6, 1 + w_avg = −0.01927 at both. Changed: the sentence
labels −0.008234 as D7's (Q² = 186) and gives the measured value beside it,
with the recipe, dated `V-CODE-3 (2026-11-10)`.

**CODE-4 (MAJOR) — APPLIED.** Claim: PHYSICS_CHANNELS' `TaggedDeuteronP` row
says the ct18nlo|mstw + nnpdfpol pair needs `--x-max 0.95` under helicity-flip;
A3 lifted it. Verdicts: reproduction and rules confirmed. Measured here
(`--channel tagged-d-p --plan helicity-flip --events 200 --seed 1234`):
configs 0 / 1 exit 0, clip 2.635e−5 % / 2.092e−5 % of the rate (max unclipped
A₁ⁿ 4.385 / 4.436); config 1 on `--unpol-sf mstw` exit 0, 2.237e−5 % (4.912).
Changed: a dated `V-CODE-4 (2026-11-10)` clause in the same cell (until A3;
builds at 1.0 since, with the measured clip shares).

**CODE-5 (MINOR) — APPLIED.** Claim: USAGE §10's six-row refusal table is
D7's (2026-10-29) with no label, and the comment in `python/lipolgen/cli.py`
`main` still speaks of six refusals and "NOT clamped". Verdict: reproduction
confirmed. Measured here (`--plan helicity-flip --events 200 --seed 1234`): d
inclusive configs 0 / 1 exit 1 at x = 0.955, Q² = 534.5 / 1536,
1 + w_avg = −0.05837 / −0.02078; d config 2, `tagged-d-p` 0 / 1 and
`tagged-6Li-alpha` 0 exit 0 — P2 §1's POST column. Changed: the sentence that
introduces the table gains a dated caption (the table is D7's, before A1 and
A3, with today's two refusals and the four builds beside it); the cli.py
comment gains a dated A3 clause on its last line (g₁, not W, is clipped; two
of the six still refuse).

**CODE-6 (MINOR) — APPLIED.** Claim: `include/lipolgen/tagged.hpp` (the block
above `FsiWeight`, and the `DEUTERON_AV18` doc) and `include/lipolgen/beams.hpp`
(the C5.5b paragraph) state the pre-A1 identity of the inclusive constants
in the present tense. Verdict: reproduction confirmed. Changed (dated
`V-CODE-6 (2026-11-10)`, same lines): the product named
`LI6_CLUSTER_POLARIZATION_HULTHEN` since A1; "the names, values and meanings
are unchanged" scoped to until A1 (since then 0.848, `P_D_DEUTERON` = AV18's
0.0576, the Hulthén 0.045 = `P_D_DEUTERON_SCENARIO`); "THEY DRIFT APART"
scoped to until A1 (since: +4.6 % vector, 0 rank-2); the `DEUTERON_AV18` doc
names `P_D_DEUTERON_SCENARIO`; beams.hpp's "stayed on `P_D_DEUTERON` = 0.045"
names the constant's A1 name.

**CODE-7 (MINOR) — APPLIED (a), (b), (c), (d).** Claim: refactor residue.
Verdict: reproduction confirmed (all four). (a)
`validation/dump_polligen_reference.py`: the original `build_tagged()` is
renamed `_build_tagged_hulthen()` on its own line; the alias line and the
`# noqa: F811` became dated comments; the A13 wrapper `build_tagged()` is the
only definition of that name (`grep` of duplicate `def`s: none). Checked: a
dry import (module loaded by spec, not run) resolves `build_tagged` to the
wrapper, which calls `_build_tagged_hulthen()`; `repin_tagged_from_lipolgen.py`
agrees (P_D rel 0, dilutions 5.389e−9 / 3.196e−9 inside 1e−8, "nothing
written"); the dump itself was NOT run; `validation/reference/` unchanged.
(b) `include/lipolgen/pipeline.hpp`: the namespace close, `#include <map>`
and re-open before `run_provenance` replaced line for line by a dated comment
(`<map>` is in the include list at the top); the same close/re-open in
`src/core/pipeline.cpp` before `run_provenance` replaced line for line —
line-neutral, so applied rather than refused; no line citation falls in
either region. (c) `default_inclusive_kernel`'s Miller branch builds the
`Li6B1` only when `b1_transfer` differs from `LI6_B1_RANK2_TRANSFER`
(conditional on the same line; the default branch reads the `static`, as
before). (d) `tests/test_b1_nuclear.cpp` T9 renamed "default_inclusive_kernel
at the pre-A1 constants, passed explicitly, is unchanged" (same line, dated),
and `docs/benchmarking/00_in_tree_checks.md` B-1 quotes the new name with the
old one. All four behaviour-free; the C++ suite passes.

**CODE-8 (MINOR) — RECORDED AS ADEQUATE.** Claim: the knob matrix has only
0 → 1 cells for `cluster_vmc_mc_sigma`, neither touching |σ| ≤ 3. Verdict:
reproduction confirmed (both sides pinned elsewhere). Not added: the matrix
asserts read / not-read against output hashes of two runs, not
`validate()`'s bounds; both sides of A11 are pinned in `tests/test_pipeline.cpp`
(±3 build and construct a `Pipeline`; 3.0000001, −4, 1000 refused naming the
band and A11; NaN refused as non-finite; and now the printed value) and in
`python/tests/test_module.py` (−3 accepted, 3.5 refused); the knob has no CLI
flag; a refused cell would test `validate()` a third time inside the suite's
slowest file (≈ 11 min, split by node id).

## 2. What was verified and how

- Re-measured before editing each statement (scratch `f/cli/`, `f/probe1.py`):
  CODE-1's labels, PHYSICS-1's not-read rows, CODE-3's three x_max runs,
  CODE-4's three tagged-d-p runs, CODE-5's four helicity-flip runs; after the
  build, PHYSICS-8's five messages (`f/physics8_probe.log`).
- Line-neutrality: `git diff --numstat` shows +n/−n for every file the phase
  edited. The citation gate's R4 reported two cited blocks EDITED in place —
  `docs/USAGE.md` 163–533 (CODE-1's first site, inside §2a, cited by the
  PHYSICS_CHANNELS `--b1-model` wiring row) and `include/lipolgen/tagged.hpp`
  54–62 (CODE-6, cited by the "Integrated dilutions" row as the "same wave
  function seen in two experiments" rule). Both rows were re-read and still
  say something true (§2a still documents the wiring; the block still states
  the Hulthén identity, now dated); re-recorded with `--record-ranges`:
  218 fingerprints, 0 new, 2 changed, 0 dropped, 0 REFUSED
  (`validation/physics_channels_ranges.json`, +3/−3). The first part-C pytest
  run before the re-record failed 6 doc-link-gate tests on exactly those two
  blocks; after it, 186 passed.
- PHYSICS-10: the changed test lines compared pairwise (70 Python/C++ lines of
  the date-only edits, plus the one in `test_b1_nuclear.cpp`): each equals its
  original with " (2026-11-03)" removed.
- Reference JSONs: `git status` shows nothing under `validation/reference/`.
- The report: no harness, kernel or report input changed (README, records and
  comments only; the harnesses' files untouched), so `make_report.py` was not
  run and `docs/benchmarking/` (REPORT.json, figures) is unchanged.

## 3. Files

Docs: `docs/CONVENTIONS.md`, `docs/PHYSICS_CHANNELS.md`, `docs/USAGE.md`,
`docs/benchmarking/00_in_tree_checks.md`, `validation/benchmarks/README.md`,
`docs/open_items/run_2026-09-03/STATUS.md` (registry row 10),
`docs/open_items/run_2026-11-03/STATUS.md` (P4b, P6, V, F rows),
`docs/open_items/run_2026-11-03/phase_P3_report_rules.md` (§1.2),
`docs/open_items/run_2026-11-03/phase_P5_port_gate.md` (§5 item 3), this
record (new). Code (comments and one message): `include/lipolgen/beams.hpp`,
`include/lipolgen/tagged.hpp`, `include/lipolgen/pipeline.hpp`,
`src/core/pipeline.cpp`, `python/lipolgen/cli.py`,
`validation/dump_polligen_reference.py`. Tests (comments, one added check,
one renamed case): `tests/test_pipeline.cpp`, `tests/test_b1_nuclear.cpp`,
`tests/test_beams.cpp`, `tests/test_sf.cpp`, `tests/test_tagged.cpp`, and the
13 Python test files of PHYSICS-10. Gate sidecar:
`validation/physics_channels_ranges.json` (2 fingerprints re-recorded).
`docs/open_items/run_2026-11-03/DECISIONS.md` carries the supervisor's three
notes of 2026-11-10, not F's.

## 4. Gates (final tree; foreground, one at a time, logs in the phase's scratch)

| gate | result | wall |
|---|---|---|
| `cmake --build build -j2` | exit 0, 0 warning lines | 143 s |
| `./build/lipolgen_tests` (full) | 435 / 435 cases, 17 244 595 assertions, 0 failed | 451 s |
| pytest A (14 files, `test_b1_model` … `test_bench_t1_r_world`) | 358 passed, 3 skipped | 184 s |
| pytest B (7 files, `test_bench_t4_t6` … `test_delta_convention`) | 136 passed, 1 skipped | 103 s |
| pytest C (7 files, `test_doc_link_gate` … `test_meta_provenance`, without `test_knob_provenance`) | 186 passed (first run, before the range re-record: 180 passed, 6 failed in `test_doc_link_gate.py`, §2) | 78 s (72 s) |
| pytest D (`test_knob_provenance.py`, 339 node ids) | 249 passed, 90 skipped | 395 s |
| pytest E (`test_knob_provenance.py`, the other 338 node ids) | 277 passed, 61 skipped | 291 s |
| pytest F (15 files, `test_module` … `test_throughput`) | 317 passed | 201 s |
| **pytest sum** (44 files; the 677 knob ids checked complete and disjoint against `--collect-only`) | **1523 passed / 155 skipped / 0 failed** | 1252 s |
| bench subset (`LIPOLGEN_BENCH=1 … -k bench python/tests`) | 368 passed, 1 skipped, 1309 deselected | 492 s |
| `repin_tagged_from_lipolgen.py` (check-only) | agrees, nothing written | < 1 s |
| dump script dry import (not run) | `build_tagged` → the A13 wrapper; no duplicate `def` | < 5 s |
| `check_physics_channels_links.py` (strict) | 1243 references, 95 ranges, 7 external, 0 broken, 24 allow-listed (+ SPIN32 19 / 115 / 0 broken; PYTHIA_BRIDGE 0 broken) | 2.9 s |
| `--records-only` | 1432 citations in 89 dated run records, 0 broken | 0.35 s |
| `--audit-ranges` | 209 range citations on 179 blocks, 0 REFUSED | 3.2 s |
| `--record-ranges` (once, §2) | 218 fingerprints: 0 new, 2 changed, 0 dropped | ≈ 3 s |
| `check_spdx_headers.py` | 136 files, all carry the SPDX line | 0.04 s |
| `test_doc_link_gate.py` (final tree) | 112 passed | 32 s |

The report was not regenerated (no harness, kernel or report input
changed); `docs/benchmarking/` is unchanged.

## 5. Open points

1. **For the supervisor — PHYSICS-2's CONVENTIONS paragraph.** Left as
   history under the bullet's own PRE-A1-record label (the reproduction
   verifier's reading). The rules verifier would turn each number into
   "until A1 … since A1 …" as USAGE does. Recommendation: leave it; if the
   supervisor wants the paragraph rewritten, it is a line-neutral edit of
   one bullet, but the bullet's head already states the current values.
2. **For the supervisor — the PHYSICS-1 site in CONVENTIONS** sits in that
   same PRE-A1-record bullet; the A17 clause was added as ruled, since the
   label covers A1's numbers and not the A17 rule. Recommendation: keep.
3. **For the supervisor — the range re-record.** Two cited blocks
   (`docs/USAGE.md` §2a, the `tagged.hpp` block above `FsiWeight`) were
   edited in place and re-recorded after the citing rows were re-read
   (§2). Recommendation: accept; the alternative, moving the clauses outside
   the cited lines, would put them away from the sentences they correct.
4. **For the supervisor — the P6 back-fill.** The brief named the P4b hash;
   F also back-filled P6's `5372f2d` (the finding's fix names it, and every
   phase back-filled the previous row). Recommendation: keep.
5. **PHYSICS-10's count** is 71 dated lines against the lens's ~54: F's scan
   also covered the C++ test files and lines whose entry token sits inside
   a parenthesis ("(A1; …"). Four lines whose block carries the dated form
   on an adjacent line were left.
6. **Unchanged, for the author (already filed):** the deuteron's inclusive
   b₁ reading the ⁶Li transfer (0.32753 × MillerB1; PHYSICS-3 labels it, does
   not decide it), and the dump's non-reproducibility (P5 §5 item 1).

## 6. The correct-lists of both lenses, verbatim

### Lens 1 (physics and honesty)


- **A1 constants and their homes:** `LI6_CLUSTER_POLARIZATION` = `LI6_POLARIZATION_VMC_SIX_BODY` = 0.848 (beams.hpp). `LI6_B1_RANK2_TRANSFER` = `_VMC` = 0.982576 (constants.hpp). `P_D_DEUTERON` = 0.0576 and `DEUTERON_VECTOR_POLARIZATION` = 0.9136. The pre-A1 values are named `_HULTHEN`/`_SCENARIO`. `deuteron_av18_p_d()` = 0.0575998920. `deuteron_scenario().eff_pol_p` = 0.9325.
- **A1 ratios reproduced from the kernel tables** at (0.3, 5): ⁶Li g₁ and g₂ ×1.0453283 (= 0.848/0.811228375), b₁ and b₂ ×1.0657619 (= 0.982576/0.921947), f₁, f₂ and Δ ×1. Deuteron g₁ ×0.9797319 (= 0.9136/0.9325). The deuteron's b₁ moves ×1.0657619 with `b1_transfer`, as the A1 note (2) says.
- **The Hulthén tag reads the scenario values:** the pipeline's d–p factory uses `P_D_DEUTERON_SCENARIO`, and beams.cpp's `DEUTERON_SCENARIO()` is built on `DEUTERON_VECTOR_POLARIZATION_SCENARIO`.
- **A1 meta and banner:** the d run's meta carries `deuteron_p_d` 0.0576 and `li6_b1_rank2_transfer` 0.982576, and its banner says the deuteron's b₁ reads the ⁶Li transfer. The ⁶Li run carries 0.848 and 0.982576.
- **DECISIONS.md is append-only:** no line was removed in any of the 8 run commits. The A1 correction of 2026-11-06 is appended under the 2026-11-03 note, and it (not the prediction) describes the tree. The default ⁶Li `--plan helicity-flip` and `--plan apar` build at x_max 1.0, and the banner reports max unclipped A₁ᵖ 1.439 and A₁ⁿ 4.436. USAGE "The top x cell" says "it builds at 1.0 (… +0.1210)".
- **A3:** `PositivityClippedG1` wraps g₁ᵖ and g₁ⁿ at ±F₁ on the kernel's own base and R. It is installed once per `InclusiveKernel` (xsec.cpp, the constructor), so `g1a`/`g1_nucleus` and g₂^WW read it. The sampler's "NO CLAMP" comment and refusal text say what the code does. The banner's clip share (3.864e−5 % = 3.864e−7) matches P2 §1. The definition section precedes the measurements in the record. The `medium_ratio` analysis is recorded as code reading: `g1_nucleus` is called from xsec.cpp only.
- **A4 outcomes re-run:** `tagged-d-p` helicity-flip builds. `--isotope d --plan helicity-flip` config 1 is refused (1 + w_avg −0.02078; vector −0.3019 at A₁ 0.9136, rank-2 −0.7189). Config 0 is refused (−0.05837). Config 2 builds, as does `tagged-6Li-alpha` config 0. Config 1 builds at `--x-max 0.95`. All of these match P2 §1.
- **A6:** the TIGHT saturations in REPORT.json, df56035 → HEAD, are exactly the amendment's numbers: 0.3207/0.3808/0.4281 → 0.1603/0.1904/0.2141; 3.016/2.933/2.938 → 1.508/1.466/1.469; 3.114/2.919/2.822 → 1.557/1.459/1.411; 9.05/8.80/8.81 → 4.52/4.40/4.41; 9.34/8.76/8.47 → 4.67/4.38/4.23. `saturation_tight_v6` equals the pre-A6 values.
- **A7:** `m2_error_bound` = `QUAD_REL_TOL` × |M₂|. Errors 3.26e−14 / 4.90e−14 / 5.53e−14 sit 448× / 302× / 271× inside, and the record says the bound was chosen after the values were seen.
- **Report tallies:** 70982c9 14/9, df56035 14/9, 5cd19d2 15/8, HEAD 15/8. Deuteron statics: headline av18 PASS, the declared all-waves sub-row FAIL. Deeps: av18/lc is in the verdict. BONuS: av18-40mb is in the verdict.
- **A9:** CMakeLists.txt raises `message(FATAL_ERROR …)` when `lhapdf-config` is missing and `LIPOLGEN_WITH_LHAPDF` is ON, with a dated A9 comment.
- **A10:** HEPMC3_CONVENTION.md's current sha256 20de1d3b… equals `DOC_SHA256_AT_AUTHORING`. The pre-A10 digest 56d3a91a… (= `git show 70982c9:` of the file) is in the comment. The diff is +32/−1, the one changed line being the stated "byte-identical" sentence. On a `tagged-6Li-alpha` run, all 13 `lipolgen.*` attributes equal the npz meta values. On an inclusive run the A1 factors appear as attributes.
- **A11:** ±3 build; 3.0000001 and −3.1 are refused; NaN is refused as non-finite.
- **A12:** the banner fractions re-measured at 2000 events: ⁷Li 4.41 / 9.52 %, ⁶Li 0.70 / 5.25 %, d–p 1.35 %, each with "meta fsi_clipped_grid_fraction; an unpriced systematic, A12". Registry row 28 is present.
- **A17 / P4b:** coherent, d inclusive, tagged-d-p and ⁶Li inclusive knob tables have only Read or NotRead rows. A moved value is refused (`b1_model=cdks` on coherent; `cluster_wave=hulthen` on inclusive). The three-ways test was replaced by `test_an_explicit_default_on_an_unread_knob_is_handled_one_way`.
- **A13:** `git diff 70982c9 HEAD -- validation/reference/` is tagged.json +15/−0 only. The block's P_D 0.019354933183319687 equals the sibling's `li6_vmc_tables()` printed split 0.015861/(0.80362 + 0.015861) (sed/grep of evgen/polligen/tagged.py at 4726812; computed here). The quantum numbers (1, 1, 0, 1) match the sibling's `TaggedChannel(…, 1.0, 1.0, 0.0, 1.0, li6_vmc_waves(), …)`. I re-measured the residuals against `TaggedModel(li6_alpha_channel(source=VmcAV18))`: 5.3888e−9 / 3.1958e−9. DECLARED.txt predates dump.log. The drift is filed for the author (P5 §5 item 1).
- **A16:** USAGE §8 carries the date, the idle-machine condition, the record pointer and the 2026-10-30 history. The medians recomputed from the raw logs: 0.540 / 0.3305 / 1.071. The A/B ratios ×0.98 / ×0.96 / ×0.98 check out. The record's raw log is identical to run_1103/p6/throughput.log except the summary header line.
- **Labelling:** every DECISIONS entry and note, and every phase record header, labels the choice as the supervisor's under the author's delegation and quotes the author verbatim.
- **Pre-existing physics items are filed:** the deuteron reading the ⁶Li transfer (A1 note (2) and correction; P1 §5 item 2, P2 §5 item 2, P3 §5 item 1), the toy b₁ behind the surviving top-cell refusals (P2 §5 item 2), the NLO/LO development item (A3; P2 §5 item 5), and the dump's non-reproducibility (P5 §5 item 1). STATUS.md has no "For the author" list yet; phase S is TODO.

### Lens 2 (code and regression)
- A1 parameterisation: `li6_ion` / `deuteron_ion` / `DEUTERON_SCENARIO()` defined in beams.cpp, header defaults (`LI6_CLUSTER_POLARIZATION` 0.848, `P_D_DEUTERON` 0.0576) equal the pybind `py::arg` defaults; `default_inclusive_kernel(..., b1_transfer = LI6_B1_RANK2_TRANSFER)` header = binding = Python `shipped_inclusive_kernel(b1_transfer=None -> LI6_B1_RANK2_TRANSFER)`; `Li6B1` / `b1_li6_from_deuteron` bindings default to the same constant. Hulthen chain: `li6_alpha_channel` Hulthen -> `dis_target = DEUTERON_SCENARIO()`; `deuteron_channel` default and binding `P_D_DEUTERON_SCENARIO`; `BreakupOptions::p_d = P_D_DEUTERON_SCENARIO`; the Pipeline d-control factory passes `P_D_DEUTERON_SCENARIO`; `DEUTERON_AV18()` only overwrites eff_pol. Port gates pass the old values explicitly (test_reference.cpp `polligen_ion`, test_kernel_reference.py `_polligen_ion`, T9 with both Hulthen constants). `li6_ion(LI6_CLUSTER_POLARIZATION_HULTHEN)` reads back 0.811228375 exactly (Z*eff_pol_p == constant, measured), `deuteron_ion(0.045).eff_pol_p == DEUTERON_VECTOR_POLARIZATION_SCENARIO`.
- A1 meta/banner: `inclusive_polarization_factors` reads the run's own kernel's Ion; a caller kernel at the Hulthen value reports 0.811228375 "(the pre-A1 Hulthen cluster product, passed explicitly)" and the transfer as not knowable (no meta key); default 6Li run meta li6_cluster_polarization 0.848, li6_b1_rank2_transfer 0.982576; default d run deuteron_p_d 0.0576 (all measured).
- A3 clip: `PositivityClippedG1` holds shared_ptrs to the inner PolSF and to `nf2_.base()` (no dangling reference), is stateless/const (thread-safe), uses `options.r_func` and the kernel's own UnpolSF (the same object `NuclearF2` uses, same F2/(2x(1+R)) expression), returns the inner value untouched where it does not bind (plain comparisons, no min/max arithmetic; NaN passes); `g1_nucleus` and WW `g2p/g2n` are PolSF base methods calling the virtual clipped `g1p/g1n` (no backend overrides them); the tagged struck kernel is built by the same `InclusiveKernel` ctor (pipeline.cpp line 289) so it is clipped; double wrapping is prevented by the unwrap in the ctor; `pol_sf_grid_report` casts `cfg.pol_sf_obj`, not the wrapped kernel object, so D7's floor report is unaffected. Clipped-share sums over `cell_rate_weights_pb()` and `dis_sampler().x_cells()/q2_cells()` of the same sampler; NaN where g1 is not read.
- A3/A4 CLI outcomes of phase_P2 sec. 1 POST, re-run here (200 events, seed 1234): d inclusive config 0 / 1 refused (Q2 534.5 / 1536, 1 + w_avg -0.05837 / -0.02078); d config 2, tagged-d-p 0 / 1, tagged-6Li-alpha 0, 6Li inclusive config 1 helicity-flip and apar build, clip shares 3.536e-7 / 2.635e-7 / 2.092e-7 / 4.793e-7 / 3.864e-7 / 3.864e-7 as recorded. P1 helpers `_x95`, `_spec_x_max` and a Spec `x_max`: no residue (grep).
- A9: CMakeLists FATAL_ERROR only inside `if(LIPOLGEN_WITH_LHAPDF)` when `lhapdf-config` is not found; p4/build_nolhapdf.log ends "Configuring incomplete, errors occurred! exit 1 wall 0.21 s"; build_off_configure/build/tests logs: configure exit 0, build exit 0 with 0 warning lines, 417/417 test cases, 17 243 193 assertions -- as phase_P4 sec. 4 says. ci.yml's CI and wheel builds pass the deps prefix (no OFF-less build without LHAPDF).
- A10: `set_run_attributes` stores the map and the first `write` attaches it right after `ensure_run_info`, before the GenRunInfo is serialized; a later call throws; keys come from `run_provenance`, from which bindings fill the npz meta (`pv(...)`). `sha256sum docs/HEPMC3_CONVENTION.md` = 20de1d3b...ad0f = `DOC_SHA256_AT_AUTHORING`; C24 anchor (262, "`lipolgen.<key>` holds exactly ...") is line 262, the section starts at 255 and the file has 284 lines; C24 compares numeric keys by float, matching the doc's "%.17g so they round-trip" (measured: the inc6 HepMC file's 15 lipolgen.* lines equal the npz meta, the one textual difference 0.84799999999999998 vs 0.848 being float-equal). The rewritten line 232 says the weight block (what T15/T15b assert) is unchanged.
- A11: `validate()` refuses only `fabs(sigma) > 3.0` after the non-finite check; USAGE line 1399 and PHYSICS_CHANNELS line 708 state |sigma| <= 3, refused beyond.
- A12: banner `fsi_clip_line` prints `fw.clipped_grid_fraction()`, the same accessor bindings write to `meta["fsi_clipped_grid_fraction"]`; tagged-7Li-alpha glauber-cluster 2000 events: banner 4.41 %, meta 0.04405694; registry row 28 appended after 27, no row renumbered; the two OPEN_ITEMS_SOLUTIONS "27 rows" sentences carry "(28 since A12)" (docs/benchmarking/08_channel_priority.md line 9's "27 rows" sits in its dated "Written 2026-10-26 ... as it stood" paragraph).
- A17: every `KnobStatus::Refused` passed to `add()` at `at_default == true` is converted to NotRead inside `add()` (the only `rows.push_back`), and the rc rows pass Refused only off-default; inclusive 6Li and coherent runs' meta contain no "refused" (measured); `test_an_explicit_default_on_an_unread_knob_is_handled_one_way` exists and passes (1 passed, 5.3 s); p4/ab_explicit_defaults.txt reports the five pairs arrays/meta/knob rows/HepMC identical; re-run here: `--cluster-wave vmc` vs unnamed inclusive 6Li, `cmp` identical for npz and HepMC; USAGE's knob-provenance example re-run: 69 knobs, read (16), not read at their defaults (51), no refused line; price_defaults.py moves only the three cluster_wave=vmc probes to not-read, r_source=r1998 (non-default) stays refused; CONVENTIONS "Refuse or label" states the one rule.
- A13: `git diff --numstat 70982c9..HEAD -- validation/reference/` = tagged.json 15/0 only; validation/README.md +24/-0; `_manifest.json` unchanged and still lists tagged.json; C++ case "tagged: the VMC alpha-d dilutions against polligen (A13, rtol 1e-8)" passes (11 assertions), `kVmcPortRtol = 1e-8` and the file's `port_gate_rtol` checked equal; Hulthen loops in test_tagged.cpp (3) and test_cluster.cpp skip the key, `model_by_key` throws on unknown, those cases pass (12 cases, 516 assertions); test_port_gate_li6_vmc.py 4 passed; repin script is check-only (prints "nothing written", SystemExit on disagreement); dump script carries `build_tagged_li6_alpha_vmc`; no other test reads tagged.json's channels.
- A5/A6/A7: `DEFAULT_WAVE` is per harness (deeps "av18", boeglin "hulthen" unchanged); `saturation_tight` has no `PRE_D4_OVER_JM`, `saturation_tight_v6` keeps it; `m2_error_bound = QUAD_REL_TOL * |M2|` and the test asserts err <= bound, the 1e-15 pin kept on the toy integrand; pins carry dated A-entry comments.
- Regression sweep: build logs of P1, P2, P4, P4b, P5 and the OFF build have 0 warning lines; REPORT.json summary 15 pass / 8 fail and every row status equals REPORT.md's; the six figures in `git diff --stat -- docs/benchmarking/figures` are exactly the six rows with a figure whose subrows changed (t1_g1d_world, t2_bonus, t2_deeps, t2_deuteron_static, t2_hermes_b1, t5_delta_moment_bound); t2_deuteron_static is the one status move (FAIL -> PASS, A5).
