<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-03 — phase P2: the positivity bound |g₁| ≤ F₁ (decisions A3, A4)

Worker P2, 2026-11-06, on `claude/run-2026-10-26` at `3d059d3` (phase P1
committed; this phase uncommitted, the supervising session commits).
Decisions: `DECISIONS.md` A3 and A4 — the SUPERVISOR's choices under the
author's delegation ("Go with your recommendations for the 17 items, just
record the choice and why"), not the author's words. Every number below was
measured in this phase with `source env.sh` on this machine (4 cores, one long
gate at a time); scripts and logs are in the session scratch `run_1103/p2/`
(`ab_runs.sh`, `ab_compare.py`, `a1_dump.py`, `split.py`, `pre/`, `post/`,
`a1_pre.json`, `a1_post.json`, `split_post.json`, `flip_*.log`, `final_*.log`,
`report.log`, `pixdiff.txt`, `result.json`). **pre** = `3d059d3`'s build,
taken before any edit; **post** = this tree.

## 0. What changed, in one paragraph

A3 is implemented as ONE clip at the nucleon level: a new `PolSF`,
`PositivityClippedG1` (include/lipolgen/sf.hpp, src/core/sf.cpp), wraps the
backend a kernel reads and returns g₁ᵖ, g₁ⁿ clipped to ±F₁ of that nucleon,
F₁ = F₂/(2x(1+R)) on the kernel's OWN base `UnpolSF` (CT18NLO frozen at its
grid edges as D7 left it, or the toy) and the kernel's own R (`r_func`); where
the bound holds the inner value is returned as is. `InclusiveKernel`'s
constructor (src/core/xsec.cpp) installs it on EVERY kernel — inclusive and
tagged struck-cluster alike, toy, grid and custom backends alike, never twice
— so `PolSF::g1_nucleus` and the Wandzura–Wilczek g₂ (the kernel's
`g2_ww` over `g1a`) see the clipped g₁; `InclusiveKernel::g1_model()` now
returns the wrapper. The clipped share is defined in §1 before it was
measured, computed by one function (`g1_positivity_clip`, `G1PositivityClip`,
include/lipolgen/pipeline.hpp / src/core/pipeline.cpp), stored under ONE meta
key `meta["g1_clipped_rate_fraction"]` (python/bindings.cpp; NaN where g₁ is
not read), exposed as `lipolgen.g1_positivity_clip(pipeline)` (fraction, the
largest unclipped |A₁| per species, the sentence) and printed by the CLI in
one line, "g1 positivity clip (A3, |g1| <= F1 per nucleon): binds on … % of
the rate (max unclipped A1^p = …, A1^n = …)". The sampler's top-cell refusal
stays; its comment and message (src/core/sampler.cpp `build_state`) now say
that g₁ is clipped per nucleon since A3 and that what the check still catches
is the rank-2 term. A4: `x_max` stays 1.0, no CLI knob was added, no reference
JSON edited (`git status validation/reference` clean). **Measured outcome:**
the toys are bit-identical (columns and meta, bar the new key); of the seven
default top-cell refusals (the six of I7 §3.5 and P1's ⁶Li one) FIVE now build
and two stay — the inclusive deuteron at configs 0 and 1, where the rank-2
term is what drives 1 + w_avg negative at the clipped A₁; of A4's two CLI
plans `--channel tagged-d-p --plan helicity-flip` builds and `--isotope d
--plan helicity-flip` stays refused. **Contrary to the expectation written
under A1, the default inclusive ⁶Li `helicity-flip` / `apar` refusal is
LIFTED** (§5 item 1): its A₁ = 0.616 is the nucleus' (eff_pol 0.2827 per
nucleon), the nucleon A₁ⁿ there is 4.4, and the clip binds.

## 1. The definition (written before the measurement) and every number that moved

**The clipped share** (pipeline.hpp `G1PositivityClip`): over the run's
accepted cells c with w_c = `Pipeline::cell_rate_weights_pb()` (the
spin-averaged cell cross section, the channel's own rate),
fraction = Σ_{c: the clip binds} w_c / Σ_c w_c, "binds" being |g₁ᴺ| > F₁ᴺ of
the backend at the cell centre for a species N ∈ {p, n} that the kernel's DIS
target carries with a non-zero polarization — the test the clip itself makes;
plus the largest unclipped |A₁ᴺ| = |g₁ᴺ|/F₁ᴺ over those cells per species.
NaN where g₁ is not read (`pol_sf_is_read` false: the coherent channel, every
plan at λ_e P_e = 0). The definition counts cell CENTRES: g₂^WW at a cell
integrates g₁ over u ∈ [x, 1], so the clip also reaches cells below the ones
it binds at (sizes in §2).

**Nucleon A₁ as the kernel sees it** (`a1_dump.py`; per-nucleon A₁ = (g₁/F₁)
of the kernel divided by the ion's mean polarization; Q² = 10 / 1000):

| species | x = 0.83 | 0.85 | 0.871 | 0.90 | 0.93 | 0.955 |
|---|---|---|---|---|---|---|
| p, pre | 0.922 / 0.818 | 0.958 / 0.850 | 1.01 / 0.892 | 1.10 / 0.974 | 1.25 / 1.11 | 1.55 / 1.41 |
| p, post | unchanged | unchanged | **1** / 0.892 | **1** / 0.974 | **1** / **1** | **1** / **1** |
| n, pre | 1.05 / 1.02 | 1.22 / 1.18 | 1.45 / 1.40 | 1.91 / 1.83 | 2.76 / 2.65 | 4.43 / 4.40 |
| n, post | **1** / **1** | **1** / **1** | **1** / **1** | **1** / **1** | **1** / **1** | **1** / **1** |
| d, pre | 0.955 / 0.871 | 1.03 / 0.938 | 1.12 / 1.02 | 1.31 / 1.20 | 1.64 / 1.51 | 2.28 / 2.17 |
| d, post | 0.943 / 0.866 | 0.969 / 0.890 | **1** / 0.921 | **1** / 0.981 | **1** / **1** | **1** / **1** |

The toy pair: all 165 (x, Q², species) points bit-identical pre → post
(max A₁ᵖ 0.968, A₁ⁿ 0.723 at x = 0.955). At x ≤ 0.80 nothing moved on either pair.

**The seven top-cell refusals, the default pair at x_max 1.0, helicity-flip at
P_z = P_e = 0.7** (pre: the sampler's message, i.e. the FIRST refused cell;
post: the message where still refused, and `split.py`'s WORST cell over the
pure m = ±1 states at λ_e = −1 with vector = w(P_e) − w(0), rank-2 = w(0);
"cells" = refused cells of the state):

| run | pre: Q², 1 + w_avg, vector (A₁), rank-2, cells | post |
|---|---|---|
| d inclusive, config 0 | 283.7, −0.02003, −0.5794 (1.985), −0.4406, 10 | **refused**: 534.5, −0.05837, −0.5213 (0.9136), −0.537, 4; worst 733.6, −0.2275, −0.6364, −0.5911 |
| d inclusive, config 1 (= A4's `--isotope d --plan helicity-flip`) | 1007, −0.05055, −0.4101 (1.983), −0.6404, 7 | **refused**: 1536, −0.02078, −0.3019 (0.9136), −0.7189, 3; worst 1897, −0.1410, −0.3803, −0.7607 |
| d inclusive, config 2 | 1707, −0.01199, −0.2729 (1.99), −0.7391, 2 | builds; worst +0.0994 (−0.1405, −0.7601) |
| `tagged-d-p`, config 0 | 229.7, −0.005726, −1.006 (4.366), 0, 12 | builds; worst +0.3034 (−0.6966 at A₁ 1) |
| `tagged-d-p`, config 1 (= A4's `--channel tagged-d-p --plan helicity-flip`) | 1119, −0.02414, −1.024 (4.406), 0, 6 | builds; worst +0.5414 (x 0.871, −0.4586) |
| `tagged-6Li-alpha`, config 0 | 480.9, −0.02789, −1.028 (1.98), 0, 5 | builds; worst +0.3629 (−0.6371 at A₁ 0.9136) |
| ⁶Li inclusive, config 1, `helicity-flip` and `apar` (P1's) | 1897, −0.0185, −0.2578 (0.6161), −0.7607, 1 | builds; worst +0.1210 (−0.1183 at A₁ 0.2827, −0.7607) |

The two that stay are A4's "refusal that survives A3": at the clipped A₁ (the
deuteron's eff_pol 0.9136, nucleon A₁ = 1) the vector term alone leaves
1 + w ≥ +0.36 at the worst cell; the rank-2 term (the class toy b₁
`Li6B1(MillerB1)` on CT18NLO's F₁, P1 §5 item 2) is what takes it negative.
Neither term alone is below −1; the sum is. They stay (A4); `--x-max 0.95`
builds both. The toy-g₁ control `--pol-sf toy` on d config 0 is refused at the
same cell, identically pre and post (534.5, −0.009225, vector −0.4722 at A₁
0.8275, rank-2 −0.537, 4 cells): the clip does not bind there.

**The rank-2 refusals, unmoved** (same cell, 1 + w_avg, rank-2 term and cell
count pre and post; only the printed cell A₁ moved, the vector term being 0 at
λ_e = 0): the b₁ band's ×2 row, `cdks` at `--x-max` 0.95 and 0.9 (x = 0.871,
Q² 150.6, −0.01927, rank-2 −1.019, 25 cells; A₁ 0.2961 → 0.2662),
`li6-convolution --b1-unpol ct18nlo` ×2 at 0.95 (Q² 315.3, −0.005282, −1.005,
18 cells; 0.292 → 0.2627), `--inclusive-b1` on `tagged-6Li-alpha` config 1
(x 0.955, Q² 167.3, −0.1302, −1.13, 24 cells; 1.997 → 0.9136).

**The clipped share on every default run that builds** (`meta`, post):
⁶Li inclusive c1 helicity-flip / apar 3.8640e−7 (A₁ᵖ 1.439, A₁ⁿ 4.436 max
unclipped); at x_max 0.95 3.8474e−7; ⁷Li inclusive 3.6174e−7; d inclusive c2
3.5361e−7, at 0.95 (c1) 3.8546e−7; `tagged-d-p` c0 2.6348e−7, c1 2.0925e−7, c2
1.8936e−7, c1 at 0.95 2.0838e−7; `tagged-6Li-alpha` c0 4.7926e−7, at 0.95
4.7758e−7, c1 3.8640e−7; every tensor plan and the coherent channel NaN (g₁
not read); the toy runs 0. So the clip binds on ~2–5 × 10⁻⁷ of the rate.

**Rates and shapes** (4000 events, seed 1234): the A4 plan `tagged-d-p`
helicity-flip, refused pre, builds post at σ 479692.8571 pb, ⟨x⟩ 0.0296589,
⟨Q²⟩ 4.41799; ⁶Li inclusive config 1 helicity-flip at x_max 1.0 builds at
σ 496267.1016 pb, ⟨x⟩ 0.0335338, ⟨Q²⟩ 4.58278 (pre refused; at x_max 0.95
pre and post give the same columns, σ 496267.10088168003 → 496267.1008816811,
the per-category σ 496228.5804537803 / 496305.62130957976 →
496228.58033517894 / 496305.6214281833, ±2.4e−10); the ⁶Li tensor-thirds run
is bit-identical including σ 496334.156.

**Test pins re-taken** (each with an `A3 (2026-11-03)` comment):
`test_bench_t5_closures.py::test_wu_row_passes_and_verdict_follows` σ_w
79.823733 → 79.823731 and σ_cat 79.823531 → 79.823529 (measured
79.8237306934377 / 79.82352911238141; rel 1e−8 kept; the Mode-W weights carry
g₂^WW, whose integral reads the clipped g₁ above x ≈ 0.83; N, σ_u, err_w, the
pull, the rate shift and the no-op pull did not move); tests/test_pipeline.cpp
"D2: the tagged struck-cluster kernel takes a structure-function backend …"
`ScaledG1(3.0)` → `ScaledG1(1.5)` (×3 breaks the bound at x = 0.5, toy A₁ᵖ
0.616, and is clipped; the test's subject is that the slot reaches g₁);
test_sf_backend.py's three `tagged-d-p` refusal tests re-taken: the pair now
builds at 1.0 with its clip share pinned (2.0924992e−7, rel 1e−6, a
reproducibility tolerance chosen after the value was seen, said here), the
message checks moved to the inclusive deuteron config 1, which A3 still
refuses (A₁ 0.9136, 1 + w −0.02078, rank-2 −0.7189), and the CLI's
configuration-error test to the same run. No tolerance was loosened. **P1's x_max 0.95 helpers removed** (the supervisor's correction of 2026-11-06 under A1; §5 item 1): the tests they served now run the shipped window (x_max 1.0) and NO pin moved — none of them pins a number that depended on the 0.95 window (test_pzz_mode's bit-identity of the default vs pre-flag run, the `pzz_mode` / `pzz` rows and the meta record; test_sf_backend's `pol_sf` grid keys, `pol_sf_below_grid_frac` 0.0 and floor 1.0 unchanged; the knob matrix's `inclusive-6Li/helicity-flip` spec, 25 passed / 3 skipped as on P1's final tree, the three `rc` cells skipping on `RcModel`'s polarised-beam refusal as before). New:
tests/test_sf.cpp and tests/test_xsec.cpp "A3: …" (identity on the toy at nine x
× three Q², the ±F₁ values, g1_nucleus and g₂^WW on the clipped g₁, one
clip per kernel), test_meta_provenance.py `test_a3_*` (the key, NaN where not
read, the default ⁶Li helicity-flip building with its share and the banner
line).

**Unchanged and verified:** the report (`--compare`: 0 changes; 14 pass /
9 fail / 0 blocked; t5_weighted_unweighted prints σ_w to 4 decimals, 79.8237);
`B1_TOP_CELL_B1_OVER_F1` (b₁/F₁ reads no g₁); the polligen port gates (toy:
bit for bit).

## 2. What was verified and how (A/B)

`ab_runs.sh pre|post` (CLI, 4000 events, seed 1234, 42 runs) and
`ab_compare.py` (every npz array byte for byte, σ, ⟨x⟩, ⟨Q²⟩, meta key by key):

* **Toys bit-identical (13 runs):** `--unpol-sf toy --pol-sf toy` on inclusive
  ⁶Li config 1 tensor-thirds / helicity-flip / apar, inclusive d helicity-flip
  configs 0 and 1, `tagged-6Li-alpha` tensor-thirds / helicity-flip c0 /
  `--inclusive-b1`, `tagged-d-p` tensor-thirds / helicity-flip, ⁷Li
  helicity-flip, coherent, the `cdks` ×2 row at 0.95: every column identical,
  every meta key identical except the new one (0.0, or NaN where g₁ is not
  read). The mixed `ct18nlo` + toy g₁ run likewise (0.0: `ToyG1` sits on the
  kernel's own base, |A₁| ≤ 1).
* **The default pair where it built pre (12 runs):** columns identical in
  every one; the tensor plans and the coherent run identical including σ; on
  the helicity plans the per-category σ moved by ≤ 4e−10 relative (g₂^WW, §1),
  e.g. `tagged-6Li-alpha` c1 496213.2698292774 → 496213.2694151783. The `toy F₂
  + nnpdfpol` d run records 0.0 (no cell centre binds) yet moved by 1e−11: the
  WW integral above the top cell, where the toy F₁ is below NNPDFpol's g₁.
* **Seven refused → five build** (§1); the post refusal texts carry "g1
  clipped to |g1| <= F1 per nucleon since A3".
* **The polarised-EMC `medium_ratio` path:** applied after the clip in
  `PolSF::g1_nucleus`, but no shipped route passes one (`InclusiveKernel::g1a`
  calls it without, and no other caller exists in src/ or include/), so it
  cannot re-break the bound in a shipped configuration; a caller-supplied
  ratio > 1 could. Not changed.

## 3. Files

C++: include/lipolgen/sf.hpp (`PositivityClippedG1`, appended),
src/core/sf.cpp (its definition, appended), src/core/xsec.cpp (the constructor
installs it; line-neutral), include/lipolgen/xsec.hpp (`g1_model()` comment),
include/lipolgen/pipeline.hpp (`G1PositivityClip`, appended),
src/core/pipeline.cpp (`g1_positivity_clip`, appended), src/core/sampler.cpp
(`build_state`'s comment and refusal text, line-neutral), python/bindings.cpp
(the meta key on the existing provenance line; `PositivityClippedG1` and
`g1_positivity_clip` bound beside `inclusive_polarization_factors`),
python/lipolgen/cli.py (the banner line; `--x-max`, `--unpol-sf`, `--pol-sf`
help re-worded, line-neutral). Tests: tests/test_sf.cpp, tests/test_xsec.cpp,
tests/test_pipeline.cpp, python/tests/test_sf_backend.py,
python/tests/test_bench_t5_closures.py, python/tests/test_meta_provenance.py.
Documents: docs/PHYSICS_CHANNELS.md (the "g₁ and Wandzura–Wilczek g₂" row),
docs/USAGE.md (§1's meta paragraph, "The top x cell", the `tagged-d-p` pair
section, the refusal table; line-neutral), docs/benchmarking/REPORT.json /
REPORT.md (regenerated: date, revision, runtimes only), this record, STATUS.md
(P2 row). P1's x_max 0.95 helpers removed (§5 item 1): python/tests/test_pzz_mode.py and
python/tests/test_sf_backend.py (`_x95` and its five / one calls), python/tests/test_knob_provenance.py
(the `Spec` field `x_max`, which P1 added for this alone, and `_spec_x_max`).

## 4. Gates (final tree), tallies and wall times

| gate | result | wall |
|---|---|---|
| `cmake --build build -j2` | 0 warnings, 0 errors | 2 m 13 s (first), 28 s (last) |
| flip-only `./build/lipolgen_tests` | 429 / 430, 1 failed (the `ScaledG1(3.0)` pin, re-taken) | 7 m 51 s |
| flip-only pytest | 7 failed / 1508 passed / 155 skipped: t5 σ_w pin (re-taken), 3 `tagged-d-p` refusal tests (re-taken, A4), 3 `test_doc_link_gate` (two record citations of src/core/xsec.cpp line 81 broken by my first edit order; fixed by keeping the `ToyG1` line on 81) | 18 m 31 s |
| `./build/lipolgen_tests` | 432 / 432 test cases, 17 244 559 assertions, 0 failed | 9 m 19 s |
| pytest, full | 1517 passed, 155 skipped, 0 failed | 22 m 25 s |
| bench subset (`LIPOLGEN_BENCH=1 -k bench`) | 368 passed, 1 skipped, 0 failed | 8 m 40 s |
| strict citation gate | 1243 / 95 / 7 / **0 broken** / 24 (+2: the row's two new citations; SPIN32 19 / 115 / 0 / 0 / 10; PYTHIA_BRIDGE 0 / 0 / 8 / 0) | 4 s |
| `--audit-ranges` | 209 range citations, 0 REFUSED | 2 s |
| `--records-only` | 1432 citations in 84 records: 922 live, 61 into another record, 121 annotated, 129 read by hand, 199 skipped, **0 broken** | 2 s |
| `check_spdx_headers.py` | 136 files, all carry the line | <1 s |
| `pytest python/tests/test_doc_link_gate.py` | 112 passed | 38 s |
| after the helpers' removal: the three test files (test_pzz_mode, test_sf_backend, test_knob_provenance) | 611 passed, 151 skipped, 0 failed | 13 m 29 s |
| after the helpers' removal: pytest, full | 1517 passed, 155 skipped, 0 failed | 22 m 1 s |
| after the helpers' removal: `test_doc_link_gate.py`; strict / `--records-only` / `--audit-ranges`; SPDX | 112 passed (34 s); 1243 / 95 / 7 / 0 broken / 24; 0 broken; 0 REFUSED; 136 files | |
| `make_report.py --all --plots --out-dir docs/benchmarking --compare <pre REPORT.json>` | 14 / 9 / 0; **0 changes** | 2 m 1 s |
| `make_report.py --check` | OK | 1 m 46 s |

Figures: 19 regenerated, all 19 differ from the committed ones only in the
footer (the differing pixel rows start below 97.2 % of the height,
`pixdiff.txt`) and were restored.

## 5. Open points

1. **The default inclusive ⁶Li `helicity-flip` / `apar` builds again.** The
   supervisor's note under A1 expected A3 not to lift it ("A₁ = 0.616 at that
   cell"); 0.616 is the nucleus' A₁, diluted by the per-nucleon 0.2827, and the
   clip binds on the nucleons (A₁ⁿ 4.4, A₁ᵖ 1.4 unclipped): worst cell
   −0.0185 → +0.1210. P1's helpers (`_x95` in test_pzz_mode.py and
   test_sf_backend.py, `Spec(..., x_max=0.95)` / `_spec_x_max` in
   test_knob_provenance.py) were REMOVED under the supervisor's correction of
   2026-11-06 under A1: every test they served runs the shipped window again,
   and no pin moved (§1). DECISIONS.md was not edited by this phase.
2. **The two surviving refusals ride on the class toy b₁** that the
   deuteron reads (`Li6B1(MillerB1)`, P1 §5 item 2): with a deuteron b₁ of its
   own they might build or not; not measured.
3. **The share is tiny (≤ 5e−7 of the rate)** while the bound binds over
   x ≳ 0.83: the shipped window's rate sits at ⟨x⟩ ≈ 0.03. A per-observable
   share (e.g. of a high-x bin) is not computed.
4. **g₂^WW below the bound's cells moves** (≤ 4e−10 relative in σ per
   category); the meta share does not count it (§1's definition).
5. The NLO coefficient functions stay the filed development item (A3).
