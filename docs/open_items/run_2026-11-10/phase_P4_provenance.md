<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-10 — phase P4: provenance (decisions J5, J6)

Written by the P4 implementer (an Opus worker under the supervising session),
started 2026-11-12 (the "pre" A/B only; the session then ended) and carried
out on 2026-11-16, on base commit `8b45e3c` with a clean tree (phases P1, P2,
P2b, P3 of this run committed; this phase uncommitted, the supervising
session commits). It implements **J5 and J6 — the SUPERVISOR's choices under
the author's delegation** ("for the others go with your recommendation";
`DECISIONS.md` J5, J6), never the author's words; J5 reaches
`b1_clipped_rate_fraction` through the supervisor's note of 2026-11-12 under
J1, which defines that key by "the definition `g1_clipped_rate_fraction`
has". Scripts and raw logs are in the session scratch `run_1110/p4/`
(volatile); every number a reader needs is below.

## 0. What changed, in one paragraph

**J5** changes comments and documents only: the definition of
`g1_clipped_rate_fraction` (include/lipolgen/pipeline.hpp, the comment above
`G1PositivityClip`, kept line-neutral; docs/USAGE.md, the sentence that
defines the key, same line; docs/PHYSICS_CHANNELS.md, the g₁ row, same line)
and of `b1_clipped_rate_fraction` (the comment above `B1RankTwoConfinement`;
the USAGE sentence) now say that the share is a **cell-centre estimate** (a
cell is counted whole when the bound binds at its centre) and carry the
bracket measured in §1. No key, value, banner number or code path of a run
moved. **J6**: (i) `HepMC3Writer` records ONE `GenRunInfo` string attribute,
`lipolgen.run_provenance`, with the sentence `HepMC3Writer::kNoRunProvenance`
(Python `HepMC3Writer.NO_RUN_PROVENANCE`), when `set_run_attributes` was never
called before the first `write`; a call — even with an empty map — suppresses
it. (ii) The classification of the two `*_sf_extrapolation` keys moved into
the core: a virtual `extrapolation()` on `UnpolSF` and `PolSF` (default
`SF_BACKEND_OWN_EXTRAPOLATION`, "the backend's own (not an LHAPDF grid)"),
overridden by `LhapdfSF` / `LhapdfG1` with `LHAPDF_EXTRAPOLATOR` in the LHAPDF
tier (src/lhapdf/lhapdf_sf.cpp), and one core function `sf_extrapolation(p,
pol_row)` (pipeline.hpp / pipeline.cpp) that keeps today's order —
caller-supplied kernel, toy selector, the object's own answer, no object, and
the not-read label of the `pol_sf` row — so `run_provenance` returns the two
keys, `sf_grid_extras` (python/bindings.cpp) reads the same function, and
`hepmc_run_attributes` is now `run_provenance` itself. Every string is
byte-identical on every path measured (§2); the core does not link LHAPDF
(checked with `ldd` and the OFF build, §2). (iii) The one edit of
docs/HEPMC3_CONVENTION.md, inside "Run provenance attributes"; T6's
`DOC_SHA256_AT_AUTHORING` re-pinned with the A10 digest kept in a
`J6 (2026-11-10)` comment; C24's range 255-284 → 255-306 (anchor line 262
unmoved); new clause **C25** (the absence attribute, read back). (iv) The
`-DLIPOLGEN_WITH_LHAPDF=OFF` build out of tree builds and passes its tests.

## 1. Every number or label that moved, before → after

### J5 — the bracket (measured, scratch `j5/`)

Recipe: a scratch pybind11 module (`j5ext.cpp`, linked against
`build/libLiPolGenCore.so`, never added to the tree) takes the `Pipeline` the
CLI builds (`lipolgen.cli.main` with `--events 20 --seed 11 --npz`, the
`Pipeline` captured) and, for every accepted cell, evaluates the shipped
test at the cell centre and at its four corners, the edges taken from
`dis_sampler().grid()` (`LogGrid::x_edges` / `q2_edges`): for g₁ the
expressions of `g1_positivity_clip` (|g₁ᴺ| of the clip's inner backend >
F₁ᴺ of the clip, for each species the kernel's ion carries with non-zero
polarization); for b₁ those of `b1_rank2_confinement` (F₁ > 0, window
[−(3/2)F₁, 3(F₁ − |g₁|)] not empty, the unconfined b₁ outside it, with F₁,
g₁ from `kernel().tables(x, Q²)` at the corner and from the sampler's cell
table at the centre). Weights are `cell_rate_weights_pb()`, as the shipped
share. **centre** = the shipped definition; **lower** = bound at the centre
AND all four corners; **upper** = at the centre OR any corner. The centre
column reproduces `meta[...]` of the same CLI run to 0.4–5.2e−15 relative
(summation order; equal at the printed `%.4e`) on all 15 + 6 runs; for b₁ the
centre recomputed from `kernel().tables` at the centre gives the same cells
as the sampler's table on all six runs.

**g₁ share** (`g1_clipped_rate_fraction`, all `--plan helicity-flip`, x_max 1.0):

| run | config | centre (= meta) | lower | upper | cells c / lower / upper | rows (species binding) |
|---|---|---|---|---|---|---|
| inclusive ⁶Li | 0 | 4.7926e−7 | 4.7758e−7 | 6.6455e−6 | 40 / 25 / 70 | 0.871 (n, p+n), 0.955 (p+n); upper adds 0.794 (n) |
| inclusive ⁶Li | 1 | 3.8640e−7 | 3.8474e−7 | 5.1014e−6 | 59 / 35 / 100 | same |
| inclusive ⁶Li | 2 | 3.5361e−7 | 3.5205e−7 | 3.6343e−6 | 59 / 35 / 98 | same |
| inclusive ⁷Li | 0 | 4.4972e−7 | 4.4815e−7 | 6.2551e−6 | 40 / 25 / 70 | same |
| inclusive ⁷Li | 1 | 3.6174e−7 | 3.6019e−7 | 4.7905e−6 | 59 / 35 / 100 | same |
| inclusive ⁷Li | 2 | 3.3747e−7 | 3.3599e−7 | 4.4372e−6 | 59 / 35 / 100 | same |
| inclusive d | 0 | 4.7920e−7 | 4.7752e−7 | 6.6431e−6 | 40 / 25 / 70 | same |
| inclusive d | 1 | 3.8711e−7 | 3.8546e−7 | 5.1105e−6 | 59 / 35 / 100 | same |
| inclusive d | 2 | 3.5361e−7 | 3.5205e−7 | 3.6343e−6 | 59 / 35 / 98 | same |
| tagged-d-p | 0 | 2.6348e−7 | 2.6258e−7 | 3.7931e−6 | 40 / 25 / 70 | n only (struck neutron) |
| tagged-d-p | 1 | 2.0925e−7 | 2.0838e−7 | 2.8682e−6 | 59 / 35 / 100 | n only |
| tagged-d-p | 2 | 1.8936e−7 | 1.8855e−7 | 2.0182e−6 | 59 / 35 / 98 | n only |
| tagged-6Li-alpha | 0 / 1 / 2 | = inclusive ⁶Li to the printed digits | | | | |

Where it binds (⁶Li config 1; the others alike): the x = 0.871 row binds at
the centre and at all four corners (35 cells, 3.8474e−7: it IS the lower
bracket); the x = 0.955 row binds at the centre and at its two x = 0.912
corners but not at its two x = 1 corners, where the shipped F₁ = g₁ = 0 (the
bound holds as 0 ≤ 0; measured F₁ = g₁ = 0 at x = 1, Q² = 100 and 1000), so
it drops out of the lower bracket (24 cells, 1.65e−9); the x = 0.794 row
binds at no centre but at both upper corners (x = 0.832) on every accepted
cell (41 cells, 4.715e−6 of the rate) — the neutron's bound, which A3 recorded
binding "from x ≈ 0.83". So the centre estimate sits 0.4 % above the lower
bracket and the upper one is ×13 above it, entirely from that one row.

**b₁ share** (`b1_clipped_rate_fraction`, default inclusive deuteron; identical
on the default plan and on `--plan helicity-flip`):

| config | centre (= meta) | lower | upper | cells c / lower / upper |
|---|---|---|---|---|
| 0 | 4.7920e−7 | 5.4132e−9 | 5.1532e−7 | 40 / 3 / 43 |
| 1 | 3.8711e−7 | 1.4567e−8 | 4.8745e−7 | 59 / 13 / 73 |
| 2 | 3.5361e−7 | 1.5014e−8 | 4.5729e−7 | 59 / 13 / 73 |

Where (config 1): x = 0.871, 35 cells at the centre and the two upper
corners, 13 also at both lower (x = 0.832) corners (1.4567e−8); x = 0.955, 24
cells at the centre and the x = 0.912 corners, never at x = 1 (F₁ = 0 there,
excluded by the definition's F₁ > 0); x = 0.794, 14 cells at an upper corner
only (1.003e−7). Config 0: 25 / 3 / 15 / 3 cells. Every binding corner is on
the upper half, as at the centres (P2b N4).

These numbers are written, rounded to four digits, beside the two
definitions (pipeline.hpp, USAGE line of the meta keys; PHYSICS_CHANNELS
quotes the ⁶Li config-1 bracket). Nothing else moved under J5.

### J6

* HepMC3 files: the three C++ examples that hold a `Pipeline`
  (`generate_tagged`, `generate_coherent`, `generate_full`) gain exactly two
  lines each, `A lipolgen.pol_sf_extrapolation <label>` and
  `A lipolgen.unpol_sf_extrapolation nearest` (tagged and full: `not read
  under this run's fill`; coherent: `not read on channel coherent-6Li` — each
  equal to the file's own `lipolgen.pol_sf` label, as the npz rule F7 has
  it); the three that do not (`generate_inclusive`, `write_hepmc_example`,
  `hadronize_example`) gain exactly one line, `A lipolgen.run_provenance not
  recorded -- this file was written through the LiPolGen HepMC3 writer without
  a Pipeline's run attributes, so it carries no record of the physics settings
  (structure functions, conventions, models) that produced its events`. Every
  CLI-written file: unchanged, byte for byte.
* T6: `DOC_SHA256_AT_AUTHORING` 20de1d3b3cbc03a7… → **2727952ae5a24ec5…**
  (full: `2727952ae5a24ec59ca5f6eaf81a159d65517ab2ab0fd3a8c230f36634964eb8`;
  the old one kept in the comment); clauses 24 → 25; C24's range "255-284" →
  "255-306"; the row: "24/24 clauses … 24 clauses" → "25/25 clauses hold on 8
  files, 800 events … 25 clauses" (the C25 file is reported under its own key
  `c25_file`, so the file / event counts stay the CLI's). **C25, the rule
  declared before it was run** (its `CLAUSES` comment): every `lipolgen-run`
  file carries no `lipolgen.run_provenance`; one file written through
  `lipolgen.HepMC3Writer` without `set_run_attributes` (20 T0
  tagged-6Li-alpha events, `C25_N_EVENTS`) carries exactly ONE `lipolgen.*`
  run-level line before its first event, `lipolgen.run_provenance`, whose
  value is the sentence the document quotes on its anchor line (290), and
  every event reads that attribute through its `GenRunInfo`. Measured: PASS,
  n = 9 (8 + 1).
* docs/HEPMC3_CONVENTION.md, "Run provenance attributes" (lines 270–298 now;
  nothing above line 270 moved): the extrapolation keys no longer "from the
  Python module and `lipolgen-run`" only — every path that holds a `Pipeline`
  writes the same keys; the old last sentence "A writer that is not given the
  map … writes no `lipolgen.*` attribute at all" replaced by a dated J6
  paragraph: the absence attribute and its value (quoted), told apart from a
  pre-A10 file; a set file never carries it; on an inclusive deuteron run
  `lipolgen.li6_b1_rank2_transfer` is a `not read on isotope d -- …` label
  (J1) and `lipolgen.pol_sf` / `lipolgen.pol_sf_extrapolation` name the
  backend under every fill, the tensor plans included (the J1 note of
  2026-11-12); `g1_clipped_rate_fraction` and `b1_clipped_rate_fraction` are
  npz `meta` keys only. Each statement checked against a file written in this
  phase (§2).
* Test pin re-taken: tests/test_hepmc.cpp "A10: run provenance attributes
  round-trip through ReaderAscii" asserted "no map, no attribute"; it now
  asserts exactly one `lipolgen.*` attribute, `lipolgen.run_provenance` =
  `kNoRunProvenance`, on both events of a bare file, and none on the set file
  (`J6 (2026-11-10)` comment). test_bench_t4_t6.py pins `C25.files == 9`.

## 2. What was verified and how

* **A/B, the run list** (`run_1110/p2/scripts/ab_runs.sh`, the corrected list
  phase P2b used: 22 CLI runs + 6 examples + 24 inclusive d plan × config at
  x_max 0.85 + 24 inclusive ⁶Li plan × config; 300 events, seed 11). "Pre"
  taken on the clean tree at `8b45e3c` before the first edit (2026-11-16; an
  earlier pre set of 2026-11-12 was taken with a stale list — no ⁷Li output —
  and was discarded); exits: all 0 but `inc6_hul` (1, the A17 refusal), pre
  and post. `ab_compare.py`: **all 69 npz identical (every array byte for
  byte, every `meta` key and value, every knob row); all 69 CLI HepMC3 files
  byte-identical**; the six example files differ only by the added lines of
  §1 (`diff`: pure insertions in the run block, the attribute map's sorted
  order). Run twice post (after the first build and on the final build):
  the same result, the post files identical between the two.
* **The C++ keys equal the npz's**: python/tests/test_hepmc.py
  `test_j6_extrapolation_keys_and_the_absence_attribute` — on default
  inclusive d (tensor-thirds), inclusive ⁶Li helicity-flip, coherent and
  tagged-d-p helicity-flip with `--unpol-sf toy --pol-sf toy`, the file's two
  keys equal the npz `meta`'s (nearest / nearest; label / nearest; the
  coherent label; "none (closed-form toy)" ×2); no `run_provenance` and no
  clip share among the attributes; `hepmc_run_attributes(p) ==
  run_provenance(p)`; a file written by `lipolgen.HepMC3Writer` with no map
  has exactly the one line, on every event. The document's deuteron
  sentences were also read off `ab/post/incd.hepmc` and
  `d_{tensor-thirds_c0, azz_c1, helicity-flip_c2}.hepmc`
  (`lipolgen.li6_b1_rank2_transfer not read on isotope d -- the deuteron's
  rank-2 slot is MillerB1 itself, per nucleon, unmapped (J1, 2026-11-10)`,
  `lipolgen.pol_sf nnpdfpol`, `lipolgen.pol_sf_extrapolation nearest`), the
  npz of `incd` carrying `b1_clipped_rate_fraction` 3.871134487837413e−7 and
  no file carrying either share.
* **The strings on the paths the A/B does not reach**: tests/test_pipeline.cpp
  "J6: run_provenance carries the two extrapolation keys, strings as before"
  (toy selector read / not read, and a caller's `ToyF2` / `ToyG1` under
  `Custom`, which keeps "the backend's own (not an LHAPDF grid)" — the toys do
  NOT override the virtual, see §5 item 1); tests/test_lhapdf.cpp "J6: LhapdfSF
  / LhapdfG1 answer LHAPDF_EXTRAPOLATOR for themselves"; "J6: set_run_attributes
  with an empty map writes no lipolgen attribute". The caller-supplied-kernel
  branch is the pre-J6 code moved verbatim and is not separately exercised.
* **No LHAPDF in the core**: `ldd build/libLiPolGenCore.so` lists no LHAPDF and
  `nm -D` no LHAPDF symbol; the overrides are defined in
  src/lhapdf/lhapdf_sf.cpp (the tier library).
* **(iv) OFF build out of tree** (`run_1103/p4/build_off.sh <scratch>`;
  `./build` untouched): configure with `-DLIPOLGEN_WITH_LHAPDF=OFF`,
  build, its `lipolgen_tests`: configure exit 0; build exit 0, 0 warnings
  (156 s); 422 cases / 17 243 435 assertions, all passed (354 s) — run by the
  supervisor, 2026-11-17 (§4).
* **Line neutrality**: pipeline.hpp (cited at the `G1PositivityClip` line,
  which stays where it was), sf.hpp (every addition appended to an existing
  line), docs/USAGE.md (one line edited in place), docs/PHYSICS_CHANNELS.md
  (one row in place), hepmc_writer.hpp (the comment block re-flowed in its
  own lines); new code only at the ends of pipeline.hpp, pipeline.cpp,
  lhapdf_sf.cpp and the test files. The three citation modes after the edits:
  0 broken, 0 REFUSED (§4).

## 3. Files

include/lipolgen/sf.hpp, src/lhapdf/lhapdf_sf.cpp (the virtual and its LHAPDF
answers); include/lipolgen/pipeline.hpp, src/core/pipeline.cpp
(`sf_extrapolation`, `run_provenance`'s two keys, the A10 comment, the J5
comments); python/bindings.cpp (`sf_grid_extras`, `hepmc_run_attributes`,
`HepMC3Writer.NO_RUN_PROVENANCE`); include/lipolgen/hepmc_writer.hpp,
src/hepmc/hepmc_writer.cpp (the absence attribute); docs/HEPMC3_CONVENTION.md,
validation/benchmarks/t6_hepmc3_convention.py (re-pin, C24 range, C25);
docs/USAGE.md, docs/PHYSICS_CHANNELS.md (J5); tests/test_hepmc.cpp,
tests/test_pipeline.cpp, tests/test_lhapdf.cpp, python/tests/test_hepmc.py,
python/tests/test_bench_t4_t6.py; this record; STATUS.md (the P4 row).

## 4. Gates (final tree)

*Supervisor's addition, 2026-11-17:* the worker was stopped twice by the
usage limit (2026-11-16, after the C++ tests and again after pytest part E);
the supervisor ran the gates below marked so, on the same tree (no file
edited after the worker's last build of 12:22 on 2026-11-16 except this
record), and filled this table.

| gate | result | wall | run by |
|---|---|---|---|
| build | the worker's two builds (its logs `build1.log`, `build2.log`); the final one is what every gate below ran on | — | worker |
| `./build/lipolgen_tests` | 440 cases, 17 244 828 assertions, all passed (437 / 17 244 808 before: + three J6 cases) | — | worker, 2026-11-16 |
| pytest parts A / B / C / D / E | 359 passed, 4 skipped / 138, 1 / 187, 0 / 249, 90 / 277, 61 | 196 / 109 / 79 / 410 / 263 s | worker, 2026-11-16 |
| pytest part F | 324 passed | 175 s | supervisor, 2026-11-17 |
| pytest, sum | **1534 passed, 156 skipped, 0 failed** (1533 / 156 before: + `test_j6_extrapolation_keys_and_the_absence_attribute`) | | |
| bench subset | 371 passed, 1 skipped, 1 failed in the supervisor's run — the failure `test_bench_report.py::test_the_committed_report_is_what_regeneration_gives`, because that run preceded the report's regeneration (the committed T6 row still said 24 clauses); after the regeneration `test_bench_report.py` alone: 26 passed. So **372 passed, 1 skipped** on the tree as committed (372 / 1 before) | 425 s + 98 s | supervisor, 2026-11-17 |
| report `--all --plots --compare` (baseline `8b45e3c`'s `REPORT.json`) | 15 pass / 8 fail / 0 blocked; **1 change, the T6 row** ("24/24 clauses … 24 clauses" → "25/25 clauses … 25 clauses", the sha256 prefix) and one new informational sub-row, `t6_hepmc3_convention[C25]` — what the brief expected; `--check` OK | 110 s + 94 s | supervisor, 2026-11-17 |
| figures | 19 regenerated PNGs differed from the committed ones in the footer rows only (pixel diff: every differing row within the bottom 4 % of the image) and were restored; none kept | | supervisor |
| OFF build (out of tree) | configure 0; build 0, 0 warnings; 422 cases / 17 243 435 assertions passed *(V (2026-11-17), CODE-6: the helper script itself exited 1, from its last line `tail -1 f1 f2 f3` — an invalid `tail` option, fixed to `tail -n 1` in the scratch helper; the three components are as tabulated)* | 1 / 156 / 354 s | supervisor, 2026-11-17 |
| citations | strict 1243 references / 95 ranges / 7 external / 0 broken / 24 allow-listed; `--records-only` 1432 citations in 97 records, 0 broken; `--audit-ranges` 209, 0 REFUSED | | supervisor, 2026-11-17 |
| SPDX; doc-link gate; `git diff --check` | 136 files; 112 passed; clean | 28 s | supervisor, 2026-11-17 |

## 5. Open points

1. **The toys do not answer "none (closed-form toy)" through the virtual.**
   The brief listed it among the virtual's answers, but the toy SELECTOR holds
   no object (`validate()` refuses `toy` with an object), and a caller's
   `ToyF2` / `ToyG1` under `Custom` was classified "the backend's own (not an
   LHAPDF grid)" before J6; an override in the toys would have changed that
   string on that path. So "none (closed-form toy)" stays decided by the
   selector inside `sf_extrapolation` (as the brief allows for the other
   selector-level cases), and the toy classes inherit the base answer. Whether
   a caller's toy object should say "none (closed-form toy)" is a string
   change, for the supervisor.
2. The b₁ lower bracket is small (1.1–4.2 % of the centre value) for a
   reason of the window, not of the physics: the x = 0.955 row's upper corners
   sit at x = 1, where every structure function of the shipped pair vanishes
   and the definition's F₁ > 0 excludes them; the same holds for g₁. A
   bracket with the top row's x-edge pulled inside the grid would read
   differently; not measured.
3. docs/benchmarking/00_in_tree_checks.md F-4 says `tests/test_hepmc.cpp` has
   "6 cases"; it had 7 before this phase and has 8 now. Not edited (not this
   entry's number; a census row for whoever re-takes that table).
4. HepMC3_CONVENTION's statement that the deuteron's `pol_sf` names the backend
   "under every fill" rests on P2b's configuration-level rule (the
   supervisor's ruling of 2026-11-12 under J1); if a run-level "read where it
   binds" rule is ever chosen, that sentence moves with it.
