# Phase E — measured numbers

## E0. Baseline — measured 2026-09-05, read-only, nothing in the tree touched

Every number below was produced **2026-09-05** against the existing `build/`
of this tree (no reconfigure, no clean), `source env.sh` first, by the exact
command printed next to it. Nothing here is quoted from `STATUS.md`,
`engineering.md`, or `DEVELOPMENT_PLAN.md` without independently re-running or
re-grepping it against the live tree first — see the discrepancies flagged in
§E0.3 and §E0.4, which is the point of doing that.

### E0.1 Tree state

```
$ git status --porcelain      # empty, both before and after every command below
$ git rev-parse --short HEAD
82f9450
$ git rev-parse HEAD
82f94500daf196d53243595b8dc55fd11ea489c0
$ git log -1 --format="%cI %s"
2026-09-05T05:55:31-05:00 status board: phase D commit hash
$ git branch --show-current
master
```

Tree was clean before this phase started and stayed clean through all of
E0 (no build artefact, no cache, nothing was written anywhere except this
file, at the very end).

### E0.2 Suites

| gate | command | result |
|---|---|---|
| build | `cmake --build build -j` | **no-op** — every target printed `Built target <name>` (`Consolidate compiler generated dependencies of target lipolgen_core` etc.), nothing recompiled; 0 lines matching `warning\|error` (case-insensitive) in the full log; exit 0 |
| C++ tests | `build/lipolgen_tests` | **400 cases / 17 240 262 assertions / 1 skipped / 0 failed** — exactly `[doctest] test cases: 400 \| 400 passed \| 0 failed \| 1 skipped` / `[doctest] assertions: 17240262 \| 17240262 passed \| 0 failed \|` / `Status: SUCCESS!`; exit 0 |
| Python tests | `python -m pytest python/tests -q` | **909 passed, 112 skipped in 127.78s (0:02:07)** |
| docs gate | `python3 validation/check_physics_channels_links.py` | **1145 references checked (strict), 97 ranges, 7 external, 0 broken, 6 allow-listed** (plus the same 6 named `allowed` lines as before — `xsec.hpp:220`, `pipeline.cpp:1692`, `breakup.cpp:332`, `rc.cpp:1032`, `pipeline.cpp:773`, `pythia_bridge.cpp:700`); exit 0 |

The **1 skipped** doctest case was identified positively, not inferred: forcing
it alone with `--no-skip=true --test-case="T9: the Q_N = 0 Rosenbluth limit of
Eq. (38)'s integrand"` runs it and it fails at `tests/test_rc.cpp:1630` (as of d3ac125) with
`CHECK( false )` — confirming `tests/test_rc.cpp:1608` (as of d3ac125)
(`TEST_CASE(... * doctest::skip(true))`) is the one skipped case, and that it
is a deliberate "not implemented" stub (`PolradFull` is unimplemented, per the
comment at that line), not an environment-dependent skip. The other
`doctest::skip(...)` site in the tree, `tests/test_b1_nuclear.cpp:805` (as of de1a040)
(conditional on `mstw_grid_present()`), was **not** the one that skipped —
its test case ran to completion in this environment (MESSAGE output at
`test_b1_nuclear.cpp:234/826/849/867` appears in the log), confirming the
MSTW2008 LO grid is present on this machine and that gate's own verdict row
was actually measured, not silently skipped.

### E0.3 Cross-check against `docs/open_items/run_2026-09-03/phase_D_numbers.md`

Its last table row ("after §D8, this tree") reads: doctest **400 cases /
17 240 262 assertions / 1 skipped / 0 failed** (matches E0.2 exactly), docs
gate **1145 references checked (strict), 97 ranges, 7 external, 0 broken, 6
allow-listed** (matches E0.2 exactly), pytest **888 passed, 112 skipped**
(E0.2 measures **909 passed, 112 skipped** today — same skip count, **21 more
passed**). The doctest/docs-gate numbers are bit-identical to that prior
measurement; the pytest **pass** count moved (skip count did not). The most
likely source is commit `de1a040` ("sweep findings: ... open items 12-14,
15"), the one commit between `STATUS.md`'s recorded phase-D state and today's
`HEAD`, but this was not independently traced pytest-case-by-pytest-case at
the time — flagging the discrepancy rather than asserting its cause is the
point of a measured baseline. Whatever the cause, **909 passed, 112 skipped**
is what this tree measures today and is the number the next phase should diff
against, not 888.

> **TRACED, phase F, 2026-09-05.** The cause is neither `de1a040` nor any
> commit: `git archive de1a040 python/tests` collects **1021 tests
> (909 + 112)** — the phase-D tree AS COMMITTED already ran 909/112, and
> `git diff --stat de1a040 HEAD -- python/tests` is phase E's 16 tests
> (`test_doc_link_gate.py`, `test_mstw_sf.py`, `test_release_metadata.py`,
> `test_spdx_headers.py`) and nothing else. So `888 passed` was measured
> DURING §D8, before the §D7–F fourteen-site pass and §D9 added their tests,
> and was never re-run; both `888` cells in `phase_D_numbers.md` now say so.
> The 21 are phase D's own, not a later commit's.

### E0.4 Open engineering items — verified against the live tree, not copied

Source documents named in the task: `docs/open_items/engineering.md` (a prior
session's investigation; its own header says nothing in it was committed —
verified below item by item, since two of its four findings turned out to be
stale) and `docs/DEVELOPMENT_PLAN.md` §4 (validation matrix) and §6 (open
items after day 1). Every item below was re-checked against the current tree
(grep/read of the cited file:line, or a live command) rather than trusted at
face value; **STALE / SUPERSEDED** items are marked as such rather than
silently dropped, since the next stage should know both what is still open
and what the source docs got wrong by now.

**A. From `docs/open_items/engineering.md` §A (ePIC chain smoke test)**

1. **OPEN.** `docs/open_items/engineering.md:20` — `npsim`/DD4hep nudges the
   primary electron's energy by O(10) ppm because LiPolGen writes exactly
   `generated_mass = 0` for the beam/scattered electron. Root cause verified
   still present, same lines cited: `src/core/generator.cpp:115`
   (`beam_e.mass = 0.0;`) and `src/core/generator.cpp:133`
   (`escat.mass = 0.0;`), unchanged. Non-fatal (events still process/conserve
   at Geant4's working precision); no test in this tree gates it either way.
2. **OPEN.** `docs/open_items/engineering.md:22` — `abconv -p 1` (auto energy
   detection) cannot decode a 6Li ion PDG code; the documented workaround is a
   manual preset (`ip6_hiacc_100x10`-style), which already works, but no
   6Li species entry exists in `abconv` and LiPolGen's own writer does not
   yet apply the crossing-angle/vertex afterburner upstream as an
   alternative. No code change was made toward either fix.
   **Related, and independently confirmed stale:** `docs/T2_CHAIN.md:149` and
   `:154` still read "the `abconv`/`npsim` half is unrun here" — but
   `docs/DEVELOPMENT_PLAN.md:280-282` records that same gate as **PASSED
   2026-09-02**, and `engineering.md` §A is the write-up of that very run.
   `T2_CHAIN.md`'s "unrun" line is therefore itself stale documentation, not
   just an open engineering task — this is `PLAN.md`'s own Phase E item **E3**
   made concrete with the two contradicting file:line pairs.

**B. From `docs/open_items/engineering.md` §B (packaging prototype)**

3. **OPEN, verified against the real (non-prototype) build today.**
   `docs/open_items/engineering.md:49` — the wheel's RPATH for the third-party
   libs (HepMC3/PYTHIA8/LHAPDF) is an absolute, non-relocatable path, so a
   built wheel only runs on the machine it was built on. This is no longer a
   scratchpad-only finding: the real `CMakeLists.txt:44-45` today sets
   `CMAKE_INSTALL_RPATH "$ORIGIN:${LIPOLGEN_DEPS_PREFIX}/lib"` under
   `if(SKBUILD)` (comment at `CMakeLists.txt:29-43` explains why the guard is
   there), and it is now openly documented, with the exact same caveat and
   fix path (`auditwheel repair`), at `docs/USAGE.md:36-48` — so packaging
   support **has** landed since `engineering.md` was written, but the
   portability gap it found is real and still open, just now a documented
   limitation rather than a hidden one. (A live `readelf -d` on today's
   in-tree, non-SKBUILD build shows the expected un-rewritten
   `build/`-relative + deps-prefix `RUNPATH`, consistent with the guard
   comment; the SKBUILD/wheel path itself was not re-built to re-check, since
   that is a build action, not a read of the tree.) Tracks to `PLAN.md`'s
   Phase E item **E4** ("CI and the portable wheel, written but not pushed").
4. **OPEN, same bucket.** `docs/open_items/engineering.md:51` — PYTHIA8's
   `xmldoc` and LHAPDF's grid files are not vendored in the wheel; confirmed
   still true and documented at `docs/USAGE.md:22-31` (env-var / system-install
   route only, no vendoring).
5. **STALE / SUPERSEDED — not an open item any more.**
   `docs/open_items/engineering.md:52` claimed `data/vmc/*` is "not referenced
   anywhere in `src`/`include`/`python` yet." Verified false today:
   `include/lipolgen/cluster.hpp:38-43` defines `data_dir()` /
   `$LIPOLGEN_DATA_DIR`, and `vmc` appears in `src/core/{cluster,pipeline,
   b1_nuclear,tagged,rc}.cpp`, several `include/lipolgen/*.hpp`, and
   `python/lipolgen/{configs,cli,__init__}.py`; `docs/USAGE.md:34-36` states
   the wheel vendors `data/vmc` directly since it is small. This item is
   closed and should not be carried forward.
6. **Unverified, minor, not re-checked (would require a build action).**
   `docs/open_items/engineering.md:49` also flagged the LiPolGen `.so`s
   appearing duplicated under both `lib/` and `lipolgen/` inside the built
   wheel zip, as "a redundancy to clean up." Not re-checked here — confirming
   it needs an actual `pip wheel .` run, out of scope for a read-only phase.

**C. From `docs/open_items/engineering.md` §C (licensing) — CLOSED, not open.**
   Its recommendation (line 68) was "License LiPolGen's own code
   GPL-3.0-or-later." Verified adopted: `LICENSE` at repo root is GPLv3 text,
   `pyproject.toml:13` reads `license = "GPL-3.0-or-later"`. Nothing to carry
   forward from §C.

**D. From `docs/open_items/engineering.md` §D (Sartre) — informational, not a
   LiPolGen engineering item.** No pip-installable Sartre exists, and its own
   `Nucleus::init` has no A=6 case (`engineering.md:77`); adding one is
   described as "a substantial standalone project, not a config change." This
   is a fact about an external tool, not an open task in this repository —
   listed for completeness since the task named the whole document, not
   claimed as something for the next stage to fix.

**E. From `docs/DEVELOPMENT_PLAN.md` §4 (validation matrix)**

7. **OPEN, explicitly marked so by the plan, and independently reconfirmed by
   today's test run.** `docs/DEVELOPMENT_PLAN.md:110` (as of de1a040) — item 4's neutron gate:
   "P_p = 0.866 and P_n ≈ −0.037 (the neutron half is an *open* gate, report
   it)." Today's `lipolgen_tests` run prints, verbatim, at
   `tests/test_tagged.cpp:426`: "P_n = -0.028 against the plans/05 gate
   -0.037 -- OPEN, the model's own value is asserted and the gate is not."
   Items 1, 2, 3, 5, 6, 7, 8 of §4 carry no "open" language and are not
   listed here.

**F. From `docs/DEVELOPMENT_PLAN.md` §6 (open items after day 1)**

8. **OPEN, default confirmed unchanged.** `docs/DEVELOPMENT_PLAN.md:201-210`
   — the triton remnant's *default* is still the crude, flagged sequential
   model; `--triton-sf ciofi-simula` is opt-in only. Verified:
   `include/lipolgen/pipeline.hpp:796` (as of de1a040) still reads
   `TritonSfChoice triton_sf = TritonSfChoice::Hulthen;`.
9. **OPEN, default confirmed unchanged.** `docs/DEVELOPMENT_PLAN.md:211-216`
   — "the breakup fragments' own rescattering stays open" (no FSI moves any
   fragment's own four-vector; Glauber FSI is an optional per-event *weight*
   only). Verified: `include/lipolgen/pipeline.hpp:802` still reads
   `PipelineFsi fsi = PipelineFsi::Off;`.
10. **OPEN, as stated by the plan; not independently re-derived.**
    `docs/DEVELOPMENT_PLAN.md:223` (as of de1a040) — "Still open in the coherent sector: no
    exclusive-VM channel below M_X = 1.2 GeV." Taken from the plan text as-is;
    re-deriving the M_X floor from `src/core/coherent.cpp` was not done in
    this read-only pass.
11. **OPEN, default confirmed unchanged.** `docs/DEVELOPMENT_PLAN.md:236-247`
    — the b₁(⁶Li) convolution's A=2 gate says nothing about the α–d step
    itself ("the gate is A = 2 and says nothing about the α–d step").
    Verified: `include/lipolgen/pipeline.hpp:828` (as of de1a040) still reads
    `B1Model b1_model = B1Model::Miller;` (`Li6Convolution` is opt-in via
    `--b1-model li6-convolution`, not the shipped default).
12. **OPEN, as stated by the plan.** `docs/DEVELOPMENT_PLAN.md:247-250` — the
    coherent amplitude for polarized A > 2 has "a first rung" (eSTARlight
    unpolarized ⁶Li baseline + polarized α+d configuration sampler) but "not
    yet a full polarized coherent amplitude."
13. **OPEN, as stated by the plan.** `docs/DEVELOPMENT_PLAN.md:260-262` — the
    spin-3/2 rank-2 basis has a theory note only
    (`docs/theory/SPIN32_FINITE_GAMMA.md`); "the rank-3 sector it describes
    stays off by default, so no code behaviour changed" — no rank-3
    implementation exists.
14. **OPEN, as stated by the plan.** `docs/DEVELOPMENT_PLAN.md:253` (as of de1a040) — "Still
    wholly external: polarized nuclear PDFs." No backend for this exists in
    the tree.
15. **Known physics uncertainty, not an actionable bug — listed for
    completeness.** `docs/DEVELOPMENT_PLAN.md:256-257` (as of de1a040) — the ⁶Li effective
    polarization convention (`LI6_CLUSTER_POLARIZATION = 0.81123`) sits inside
    "the 0.81–0.85 band whose top is the Wiringa VMC 0.848," an acknowledged
    remaining uncertainty, distinct from the (already-decided) convention
    itself.
16. **Deliberate non-goals — NOT open work, listed only so the next stage
    does not mistake them for open items.** `docs/DEVELOPMENT_PLAN.md:328-338`
    ("NOT ported, deliberately"): the `reco.py`/`recopseudo.py` reconstruction
    chain, `fom.project_observables`'s D_eff *projection* (the divisor itself,
    `depolarization_effective`, IS ported), `polarized.unpolarized_emc_ratio`'s
    grid modes and `structure.NuclearF2Ratio`, `ToyG1.g2_nucleus`'s per-grid
    cache, `miller_b1_q2_scale`/`toy_b1(q2_evolve=True)`, and run 14's
    addendum figure/report scripts. None of these have a C++ counterpart and
    none is planned to.

**Cross-cutting note.** `docs/DEVELOPMENT_PLAN.md:273-275` (as of de1a040) states "Packaging:
DONE 2026-09-02" with no caveat — but items 3-4 above (from `engineering.md`
§B, reconfirmed live) show the wheel is not actually portable off the build
machine. `docs/USAGE.md:34-48` already carries the accurate caveat;
`DEVELOPMENT_PLAN.md` §6's summary line does not. Not itself a code defect,
but a documentation-accuracy gap the next stage may want to close alongside
`PLAN.md`'s Phase E items (for reference only, since `PLAN.md` was not one of
the two documents this task named): **E1** (widen the doc-link gate into a
suite), **E2** (release metadata / author name — see the ground rule on the
`<AUTHOR NAME — to be filled by the author>` placeholder and `STATUS.md`'s
author-decision batch, neither touched in this read-only phase), **E3**
(item A.2 above), **E4** (items B.3-4 above).

### E0.5 Summary line for the next phase to diff against

```
HEAD 82f9450 (82f94500daf196d53243595b8dc55fd11ea489c0), git status clean
build/lipolgen_tests:            400 cases / 17 240 262 assertions / 1 skipped / 0 failed
python -m pytest python/tests:   909 passed, 112 skipped (127.78s)
check_physics_channels_links.py: 1145 references checked (strict), 97 ranges, 7 external, 0 broken, 6 allow-listed
```

## E1. The documentation link gate: is it in the suite, and does it cover enough

### E1.1 Was the gate already run over the real document by the suite?

**Yes, already, before this phase touched anything.** `python/tests/test_doc_link_gate.py`
loads `validation/check_physics_channels_links.py` by path (`_load()`, no
`ROOT`/`DOC` override) and calls `gate.main([])` against the **real**
`docs/PHYSICS_CHANNELS.md` and the real source tree in six tests:
`test_real_document_passes_strict` (`rc == 0`), `test_real_document_reference_count`
(point count > 1000, range count > 90), `test_real_document_range_fingerprints_are_all_live`,
`test_allow_list_entries_are_all_live`, `test_allow_list_pins_are_still_on_their_lines`,
`test_real_document_external_citations`. Measured directly: `python -m pytest
python/tests/test_doc_link_gate.py -q` was **48 passed** before this phase touched
anything (re-confirmed against the baseline HEAD). So a drift in
`docs/PHYSICS_CHANNELS.md` (a moved block, a renamed symbol under a still-cited
line, a dropped allow-list pin) **already fails `pytest`**, not just a hand run
of the script — the mechanism-only framing in this task's brief ("Phase C ...
added ~37 pytests for the gate's MECHANISM") undercounted what Phase C actually
shipped; six of those tests exercise the real document, not a synthetic
fixture. Nothing needed to be added for this half of E1.

### E1.2 Other documents' `path:line` citation counts (measured with the gate's own REF/RANGE/EXT regexes)

| document | named point | range | external | total |
|---|---|---|---|---|
| `docs/USAGE.md` | 0 | 0 | 0 | **0** |
| `docs/CONVENTIONS.md` | 0 | 0 | 0 | **0** |
| `docs/T2_CHAIN.md` | 0 | 0 | 0 | **0** |
| `docs/PYTHIA_BRIDGE.md` | 0 | 0 | 8 | **8** |
| `docs/OPEN_ITEMS_SOLUTIONS.md` | 0 | 1 | 2 | **3** |
| `docs/theory/SPIN32_FINITE_GAMMA.md` | 19 | 115 | 0 | **134** |

**Decision, per the brief's own rule ("more than a handful -> extend; else say
so with the count")**: `docs/USAGE.md`, `docs/CONVENTIONS.md` and
`docs/T2_CHAIN.md` carry **zero** — nothing to extend, stated. `docs/
OPEN_ITEMS_SOLUTIONS.md` carries **3** (1 range citing `docs/PYTHIA_BRIDGE.md`'s
own citation table, 2 external PYTHIA pointers) — not more than a handful,
**stays uncovered**, count recorded here. `docs/PYTHIA_BRIDGE.md` (8) and
`docs/theory/SPIN32_FINITE_GAMMA.md` (134) both clear the bar and were
**extended** (below). Note `docs/theory/SPIN32_FINITE_GAMMA.md`'s own top
matter said, until this phase, "`docs/theory/` is outside
`validation/check_physics_channels_links.py`'s reach ... so they are not
machine-checked: treat a `file:line` here as a pointer to a symbol, and grep
for the symbol if the line does not look right" — a design note written
*because* 134 citations is exactly the shape this gate exists for; corrected
in place (§E1.4).

### E1.3 The gate, extended: mechanism

`validation/check_physics_channels_links.py` now loops over `[(DOC, ALLOW)] +
EXTRA_DOCS` (`EXTRA_DOCS = [("docs/theory/SPIN32_FINITE_GAMMA.md",
ALLOW_SPIN32), ("docs/PYTHIA_BRIDGE.md", ALLOW_PYTHIA_BRIDGE)]`), each with its
own (here: empty) allow-list, sharing one `cache`/`recorded`/`fresh`/`seen`/
`remap` across every document in one invocation — a block or upstream line
cited by more than one document hashes to **one** sidecar entry, not one per
document, which is the correct reading (the fact lives in the source line, not
in which markdown file points at it). The primary document's own report line
is printed **byte-for-byte as before** (verified: `1145 references checked
(strict), 97 ranges, 7 external, 0 broken, 6 allow-listed` — identical to
every prior measurement in this run); each extra document gets its own
`[filename]`-labelled line, and the exit code is 1 if **any** covered document
has a broken citation. A document in `EXTRA_DOCS` absent under `ROOT` (every
synthetic fixture tree the existing 48 tests build) is **skipped, not
reported missing** — only the primary document is required to exist — which
is what kept all 48 pre-existing tests green with no changes to them. Two new
tests were added to `python/tests/test_doc_link_gate.py`
(`test_extra_documents_are_covered_and_clean`,
`test_extra_document_missing_under_root_is_skipped_not_broken`); the suite for
that file is now **50 passed**.

### E1.4 What the extension found, and fixed

**`docs/theory/SPIN32_FINITE_GAMMA.md` (134 citations): 15 broken on first
strict run**, all range-shaped, none of them a point-citation problem (all 19
named point citations were already correct) — confirming the note's own
top-matter diagnosis exactly: "the Phase E code comments shifted several of
them" and they had gone unchecked. All 15 were re-derived from the live
source (function/block boundaries re-read, not guessed) and corrected:

| old (broken) | new (correct) | what it is |
|---|---|---|
| `src/core/spin.cpp:244-257` (×3) | `:245-258` | `octupole_moment` |
| `src/core/spin.cpp:259-288` | `:260-289` | `moments_along_axis` |
| `src/core/spin.cpp:228-242` | `:229-243` | `tensor_polarization` |
| `src/core/spin.cpp:333-351` | `:334-352` | `spin32_populations` (decl through the pre-comment line) |
| `src/core/bookkeeping.cpp:80-104` | `:80-109` | `helicity_flip_plan`, full body |
| `include/lipolgen/bookkeeping.hpp:212-223` | `:213-223` | `SpinTemperatureLadder` doc comment + struct + decl |
| `tests/test_tagged.cpp:335-351` | `:336-352` | the P2/polarimeter `TEST_CASE` body |
| `src/core/xsec.cpp:282-289` | `:310-317` | `density()` + `positivity_margin()` (the P8 `tensor_amplitudes`/`amplitudes` split moved these ~28 lines) |
| `tests/test_xsec.cpp:387-408` | `:387-405` | the magic-angle `TEST_CASE` (end only; start was already right) |
| 4 bare `` `:a-b` `` ranges | given an explicit path | wrapped prose had put the path-bearing citation on a different raw text line than the bare range, so `inherited_path`'s same-line rule never reached them |

Re-running `--record-ranges` then the strict gate: **19 references checked
(strict), 115 ranges, 0 external, 0 broken, 0 allow-listed** — no allow-list
was needed; every point citation and every corrected range passes on its own
merits. The doc's own top-matter (the "not machine-checked" paragraph) was
rewritten to say it now is, and to name the mechanism (§1, above the
Conventions paragraph).

**`docs/PYTHIA_BRIDGE.md` (8 external citations): 1 citation pair invisible to
any rule.** `` `BeamRemnants.cc:662`: ...; and `:935`: `` had been written as
two separate citations, the second with no filename — the EXT rule requires a
bare `Name.cc:spec`, and a `:935` alone matches neither EXT, REF nor RANGE, so
it was checked by **nothing at all** (a genuine hole the count-8 already
undercounted, since it only measures what the regex matches). Rewritten to
the single comma-joined form `` `BeamRemnants.cc:662,935` `` that
`docs/PHYSICS_CHANNELS.md` already uses for the identical upstream fact (same
two lines, same file) — now both lines are checked. Strict result: **0
references checked, 0 ranges, 8 external, 0 broken, 0 allow-listed** (all 8
external matches resolve against `pythia8317` on this machine; none skipped).

### E1.5 A side effect this phase caused and is disclosing, not fixing

Task E2's SPDX-header insertion (below) adds exactly one new line to the top
of 100 source files. Every `path:line` citation into one of those files,
anywhere in the repository, is therefore off by exactly +1 relative to before
today. This was corrected, verified, and re-recorded for all three documents
the (now-extended) gate covers — `docs/PHYSICS_CHANNELS.md`,
`docs/theory/SPIN32_FINITE_GAMMA.md`, `docs/PYTHIA_BRIDGE.md` (all measured
`0 broken` above, after the shift). It was **not** corrected in the other
documents this repository carries, because they are outside the gate's
covered set by the decision in §E1.2 (too few citations to be worth a gate:
`OPEN_ITEMS_SOLUTIONS.md`) or are dated narrative run-logs, not living
reference documents (`docs/DEVELOPMENT_PLAN.md` had zero pre-existing hits —
the only two citations it carries into a shifted file are ones this phase
*added*, already written with post-shift line numbers, so they are correct as
written and are not part of this count). Measured with the same
shifted-file-membership test used to find and fix the gated documents' own
drift: **262 `path:line` citations across 19 markdown documents** point into
one of the 100 shifted files; of those, **2** are the fresh, already-correct
citations this phase added to `docs/DEVELOPMENT_PLAN.md` (§E3b/c below), so
**260 pre-existing citations across 18 documents are now stale by exactly
+1** and were left as found:
`docs/OPEN_ITEMS_SOLUTIONS.md` (1), `docs/code_review_2026-08-29.md` (30),
`docs/open_items/engineering.md` (1), `docs/open_items/run_2026-09-02/PLAN.md` (1),
`docs/open_items/run_2026-09-02/design_C_tensor_rc.md` (17),
`docs/open_items/run_2026-09-02/estarlight_li6.md` (4),
`docs/open_items/run_2026-09-03/PLAN.md` (4), `.../STATUS.md` (4),
`.../phase_A_cdbonn.md` (1), `.../phase_A_gate_mechanics.md` (4),
`.../phase_A_miller_normalisation.md` (16), `.../phase_B_numbers.md` (8),
`.../phase_D_li7_rank2.md` (29), `.../phase_D_numbers.md` (5),
`.../phase_D_sf_injection.md` (90), `.../phase_D_small_items.md` (21),
`.../phase_E_numbers.md` (10, this file's own E0 section — a citation this
run itself measured before the shift, now stale by the same +1). All of these
are dated historical/narrative records of past phases (run logs, design
notes, status boards), not documents this repository holds to a
citation-accuracy gate; hand-correcting 260 citations across 18 files was
judged out of scope for E1/E2/E3 and too large to do with the same
per-citation verification the two gated documents received above — flagged
here, with the exact count and list, rather than silently left for a future
session to rediscover.

## E2. Release metadata

### E2.1 What existed before this phase

Measured, not assumed: `ls AUTHORS CITATION.cff` — both missing.
`grep -rl SPDX-License-Identifier include src tests python validation` —
**zero** files. `LICENSE` at repo root is the GPLv3 text; `pyproject.toml`
already declares `license = "GPL-3.0-or-later"`; `README.md` had one line,
"License: **GPL-3.0-or-later** (`LICENSE`)."; no copyright holder, no authors
list, no citation metadata anywhere.

### E2.2 Files created

- **`CITATION.cff`** — `title: LiPolGen`, `version: 0.1.0` (read from
  `CMakeLists.txt`'s `project(LiPolGen VERSION 0.1.0 LANGUAGES CXX)`, the same
  line `pyproject.toml`'s own `[tool.scikit-build.metadata.version]` regex
  reads — a comment in the file says so and names it as the source of truth),
  `license: GPL-3.0-or-later`, `repository-code:
  "https://github.com/Chao1009/PoliGen"` (from `git remote -v`:
  `git@github.com:Chao1009/PoliGen.git` — note the GitHub repository is named
  `PoliGen`, not `LiPolGen`; recorded as measured, not corrected), one
  `authors:` entry with `name:` set to the literal placeholder
  `<AUTHOR NAME — to be filled by the author>` per this task's ground rule. No
  copyright line anywhere.
- **`AUTHORS`** — the same placeholder, nothing else asserted.
- **`python/tests/test_release_metadata.py`** (8 tests, all passing) — keeps
  `CITATION.cff`'s version equal to `CMakeLists.txt`'s (regex match, both
  independently re-derived from the same file so a version bump that updates
  one and not the other fails the suite), keeps `AUTHORS` and `CITATION.cff`'s
  placeholder in sync (both-or-neither), checks `README.md`'s license section
  names `LICENSE`/`AUTHORS`/`CITATION.cff`, and (skipped, not failed, if this
  checkout has no `origin` remote) checks `CITATION.cff`'s `repository-code`
  against `git remote get-url origin`.
- **`validation/check_spdx_headers.py`** + **`python/tests/test_spdx_headers.py`**
  (6 tests, all passing) — see §E2.3.
- **`README.md`**'s license section rewritten to point at all three files and
  name the SPDX gate.

### E2.3 SPDX headers

Checked first (per the ground rule): zero files in the covered set carried
`SPDX-License-Identifier` anywhere. Added `// SPDX-License-Identifier:
GPL-3.0-or-later` (or `#` for Python) to every file under `include/**/*.hpp`,
`src/**/*.hpp`, `src/**/*.cpp`, `tests/**/*.hpp`, `tests/**/*.cpp`,
`python/bindings.cpp`, `python/lipolgen/*.py` and `validation/*.py` — **101
files** (the script's own count, itself included, since it lives in
`validation/`): 27 `include/`, 28 `src/` (27 `.cpp` + 1 `.hpp`,
`src/pythia/lhaup_dis.hpp`), 30 `tests/` (28 `.cpp` + 2 `.hpp`,
`check_close.hpp`/`json_min.hpp`), 1 `python/bindings.cpp`, 5
`python/lipolgen/*.py`, 9 `validation/*.py` (including the new script). No
copyright line was added anywhere — SPDX identifies the license, not the
holder, and the holder is still the open item in §E2.2/§STATUS.md.
`validation/check_spdx_headers.py --write` did the one-time insertion;
`validation/check_spdx_headers.py` (no args) is the gate,
`python/tests/test_spdx_headers.py::test_no_covered_file_is_missing_its_spdx_header`
runs it over the real tree in the suite. Measured after: **101 files checked,
all carry 'SPDX-License-Identifier: GPL-3.0-or-later'**, exit 0.

**Consequence measured and handled**: inserting a header as literally the
first (or, after a `#!` shebang, second) line of 100 source files shifts
every later line in each of those files down by exactly one. This broke
**every** `path:line` citation into any of those 100 files, in every document
— including the two just-extended gated documents, `docs/PHYSICS_CHANNELS.md`
(1133 point + 43 range citations affected) and
`docs/theory/SPIN32_FINITE_GAMMA.md` (17 point + 112 range). Handled as one
clean, deterministic, whole-document +1 shift (not the gate's own `--fix`,
whose MOVED-detection is reached only after an R2 blank-edge check that a
uniform −1-relative-position drift can fail first — a real, pre-existing
gate limitation, not something this phase introduced, and out of scope to fix
here) applied to every citation whose path was in the 100-file shifted set,
verified against the corrected sidecar afterward: **0 broken** on all three
gated documents (§E1.3 numbers, unchanged in every count from before the
SPDX pass). The gate's own six ALLOW-list entries (keyed on exact line
numbers) were updated the same way, including the two lines their own
free-text "reason" prose names (`include/lipolgen/pipeline.hpp:779→780`,
`src/pythia/pythia_bridge.cpp:170→171`). The 260 pre-existing citations this
same shift left stale in **documents outside the gate's covered set** are
disclosed, not fixed, in §E1.5.

### E2.4 Author name — recorded in the batch, not decided here

Per the ground rule: the git identity on this machine is a project name
(`LiPolGen`) and the configured email is not a name, so neither could stand
in for the author, and no name was guessed anywhere. Row 13 of
`docs/open_items/run_2026-09-03/STATUS.md`'s author-decision table records
this with its evidence (the two files, the version/URL provenance, and the
pytest that keeps the two placeholders from drifting apart once one of them
is filled in).

## E3. Stale gate text and the validation matrix

### E3a. `docs/T2_CHAIN.md`'s stale "unrun" ePIC gate

Measured: `docs/T2_CHAIN.md` (before this phase) said, in two places, that the
`abconv`/`npsim` half of the plans/05 5.D smoke gate "is unrun here (no
container in this environment)" / "is not run from this repository" — while
`docs/DEVELOPMENT_PLAN.md:280-282` already recorded, correctly, "HepMC3 ->
abconv -> npsim smoke test: PASSED 2026-09-02 (10/10 events through `npsim`
directly and via `abconv`; see `docs/OPEN_ITEMS_SOLUTIONS.md` §2 ...)", and
`docs/OPEN_ITEMS_SOLUTIONS.md` §2 ("Chain gate — **CLOSED** (passed)", retitled 2026-09-05) independently
carries the same record (the `npsim --compactFile
epic_craterlake_10x100.xml` direct pass, the `abconv -p 1` failure and the
`abconv -p ip6_hiacc_100x10` workaround that also passes). `docs/T2_CHAIN.md`
was the one document the run was never back-ported into. Fixed: both sites —
the "Smoke gate" bullet in §2 and open item 1 in §3 — now state PASSED
2026-09-02 with the 10/10 count and cross-reference `docs/
OPEN_ITEMS_SOLUTIONS.md` §2 and `docs/DEVELOPMENT_PLAN.md`'s own row, and item
1 is relabelled CLOSED (matching the style already used for item 2 in the
same list). Re-verified after editing: the gate still reports **0 broken**
overall (the RANGE citations `docs/PHYSICS_CHANNELS.md` makes into
`docs/T2_CHAIN.md`, at `:79-96`/`:41-45`/`:102-106`, sit earlier in the file
than this edit and were unaffected).

### E3b. `docs/DEVELOPMENT_PLAN.md` §4 validation-matrix row 7's PYTHIA cross-section clause

The clause: "PYTHIA cross section matches `gen_dis_hfs.py` at identical
settings (same 8.3 physics)." **Not demonstrated anywhere in this tree, and
cannot be, by design** — read against the T2 bridge (`docs/PYTHIA_BRIDGE.md`)
and the in-tree HFS identities:

- The T2 bridge's hard process is an in-memory `Pythia8::LHAup` at LHAup
  strategy 3. Under strategy 3 PYTHIA does not compute a physical cross
  section for the process at all — it is handed a unit weight — so
  `Pythia8::Info::sigmaGen()` reports a bookkeeping artefact of the LHAup
  interface, not a cross section. The one accessor that exposes it,
  `PythiaBridge::pythia_sigma_gen_mb` (`include/lipolgen/pythia_bridge.hpp:434`),
  says this in its own docstring: "Meaningless as a physics number here ...
  exposed only so that the example can print the bookkeeping." Confirmed by
  reading it directly — this is not an inference, it is what the header
  already states.
- The physical cross section this generator publishes is computed **upstream**
  by this generator's own kernel (`InclusiveKernel`/`dsigma_*`), never by
  PYTHIA — PYTHIA only hadronizes an event this generator already decided the
  rate of. There is therefore no PYTHIA-side cross-section number left to
  compare `gen_dis_hfs.py` against; the clause asks for a comparison the
  architecture has no data for on either side.
- What **is** demonstrated, on the same "identical settings, same 8.3
  physics" comparison the clause was reaching for: the mean **charged
  multiplicity** of the hadronic final state, in a fixed (x, Q²) window,
  against a stock PYTHIA `WeakBosonExchange` run configured with
  `gen_dis_hfs.py`'s own settings — `tests/test_pythia.cpp`, `TEST_CASE`
  "pythia: the charged multiplicity agrees with a stock WeakBosonExchange run
  in the same (x, Q2) window", asserting `|ratio - 1| < 0.20`; measured ratio
  **0.97** (`docs/PYTHIA_BRIDGE.md` §10). This is a real, in-tree,
  same-settings shower/hadronization comparison — it is just not a
  "cross section" comparison, because there is no PYTHIA cross section here
  to compare.

**Action taken**: rewrote the row to state what is demonstrated (HFS truth
identities, HepMC3 round trip, the charged-multiplicity check above, each
with its `TEST_CASE` name) and added the paragraph above explaining why the
original clause is retired rather than fixed with a test — a test cannot
demonstrate a comparison neither side of the architecture produces a number
for.

### E3c. Every other validation-matrix row now names its test(s)

Rows 1–6 and 8 previously named **no** test at all. Added, verified against
the live `tests/*.cpp` `TEST_CASE` strings (grepped, not recalled) — see the
row-by-row list in `docs/DEVELOPMENT_PLAN.md` §4 itself for the full set;
summary: row 1 -> `tests/test_spin.cpp` (4 cases); row 2 ->
`tests/test_xsec.cpp` (6 cases); row 3 -> `tests/test_sampler.cpp` (5 cases);
row 4 -> `tests/test_tagged.cpp` (4 cases, including the one that prints the
open P_n gate, `tests/test_tagged.cpp:380` (as of 71cfd60)); row 5 -> `tests/test_spectator.cpp`
(2 cases); row 6 -> `tests/test_coherent.cpp` (2 cases); row 7 -> rewritten,
§E3b; row 8 -> `tests/test_pipeline.cpp`, `tests/test_rng.cpp` (×2),
`tests/test_rc_pipeline.cpp`, `tests/test_t2.cpp` (4 cases across two files,
since the matrix row's "at any thread count" claim is exercised at the
sampler, RNG, RC-weight and full-pipeline layers separately, not by one
`TEST_CASE`).

## E4. Continuous integration and the portable wheel

Both halves of this task are **outward-facing on the author's repository**, so
both were written and exercised as far as this machine allows and neither was
pushed, uploaded or published. Nothing here contacted a remote except to
**read**: six upstream tarball URLs were probed, the tarballs and PDF sets
were downloaded, and `pip` fetched `auditwheel`/`patchelf` — no write of any
kind left this machine. The two decisions are rows 14 and 15 of `STATUS.md`'s
author-decision table. Everything on this page is also written up for a reader
who is not this run, in **`docs/PACKAGING.md`** (391 lines, new).

### E4.1 Baseline re-measured first, as required

Measured 2026-09-05 against the tree as E1–E3 left it, before touching
anything, `source env.sh` first:

| gate | command | result |
|---|---|---|
| build | `cmake --build build -j` | **no-op**, 0 lines matching `warning\|error`, exit 0 |
| C++ | `build/lipolgen_tests` | **400 cases / 17 240 262 assertions / 1 skipped / 0 failed** |
| Python | `python -m pytest python/tests -q` | **925 passed, 112 skipped in 129.27s** |
| docs gate | `check_physics_channels_links.py` | **1145 / 97 ranges / 7 external / 0 broken / 6 allow-listed**; `19 / 115 / 0 / 0 / 0` [SPIN32]; `0 / 0 / 8 / 0 / 0` [PYTHIA_BRIDGE] |
| SPDX gate | `check_spdx_headers.py` | **101 files, 0 missing** |

Identical to §E5's own record of where E1–E3 finished. Machine: Ubuntu
22.04.5 LTS, glibc 2.35, gcc 11.4.0, cmake 3.22.1, Python 3.11.4, 8 cores,
31 GiB RAM.

### E4.2 (a) The dependency stack: `.github/scripts/build_deps.sh`

The hard part of CI here is not the workflow, it is that HepMC3 3.3.0, LHAPDF
6.5.5 with three PDF sets, and PYTHIA 8.317 with Python bindings all have to
exist before a single LiPolGen target compiles. What existed on disk was two
ad-hoc scripts (`../deps/src/build_deps.sh`, `../deps/src/build_pythia.sh`,
27 lines together) with the prefix hardcoded to this machine, no download
step, and no PDF-set step at all (`../deps/lhapdf_install.log` records
`lhapdf install` failing with *No PDFs known matching patterns: CT18NLO, …*;
the sets were fetched by hand). They were merged into one prefix-parameterised
script, **178 lines** (164 non-blank), every `configure`/`cmake` flag copied verbatim.

**The URLs**, all probed 2026-09-05 with `curl -sSI -L`:

| what | URL | HTTP | bytes |
|---|---|---|---|
| HepMC3 3.3.0 | `hepmc.web.cern.ch/hepmc/releases/HepMC3-3.3.0.tar.gz` | 200 | 9 341 637 |
| LHAPDF 6.5.5 | `lhapdf.hepforge.org/downloads/?f=LHAPDF-6.5.5.tar.gz` | 200 (→ `/downloader?f=`) | — |
| PYTHIA 8.317 | `pythia.org/releases/pythia83/pythia8317.tgz` | 200 (→ `pythia8.web.cern.ch`) | 31 291 506 |
| CT18NLO | `lhapdfsets.web.cern.ch/current/CT18NLO.tar.gz` | 200 | 18 756 939 |
| NNPDFpol11_100 | `…/current/NNPDFpol11_100.tar.gz` | 200 | 20 128 894 |
| EPPS21nlo_CT18Anlo_Li6 | `…/current/EPPS21nlo_CT18Anlo_Li6.tar.gz` | 200 | 26 352 280 |

**The obvious PYTHIA URL is dead.** `pythia.org/download/pythia83/pythia8317.tgz`
returns **404**, and so does that whole `/download/` path for every 8.3
release tried (8310, 8312, 8313, 8315, 8316, 8317, 8318). `pythia.org/releases.html`
now links `/releases/pythia83/…`. A workflow written from memory rather than
from a probe would have failed on its first run, in the most expensive job.
(`gitlab.com/Pythia8/releases` also exists and carries the tag `pythia8317`,
dated 2026-01-20, as a source tree rather than a tarball — not used.)

**Checksums.** The HepMC3 and LHAPDF tarballs were still in `../deps/src/`, so
their sha256 could be taken from the artefact this tree was actually built
against and pinned in the script (`6f876091…5769d6`, `641d5ea0…f5e3a5`); the
HepMC3 URL's `content-length` **equals** that local file's size. PYTHIA's
tarball had been deleted after unpacking, so there is nothing here to hash it
against; it is left **unpinned on purpose**, said so in the script, rather than
given a checksum copied off a web page. The PDF sets are unpinned for the same
reason — `/current/` is a moving target by construction.

**One real defect found in the author's own build recipe.** HepMC3 3.3.0
defaults the install path of its Python bindings to *the interpreter's*
`site-packages`, not to `CMAKE_INSTALL_PREFIX`, and warns about it at
configure time. `../deps/src/build_deps.sh` did not override it, so `pyHepMC3`
was installed **outside the prefix** — verified: `pyHepMC3.so` sits in
`/home/cpeng/Apps/python3.11/lib/python3.11/site-packages/pyHepMC3/`, dated
`2026-08-29 10:43` (the day the deps were built), while
`../deps/install/lib/python` is an **empty directory**. Harmless there, wrong
for a cache whose whole content is supposed to be the prefix, so the script
passes `-DHEPMC3_Python_SITEARCH311="$PREFIX/lib"` (suffix computed from the
running interpreter) — which is exactly where `env.sh` already points
`PYTHONPATH`. Verified after the rebuild: `<prefix>/lib/pyHepMC3/` now exists
inside the prefix. **This was caught by watching the first run's configure
output and killing it before `make install` could write outside the
scratch area**; nothing outside the scratch directory was modified except one byte-compiled cache the HepMC3 stage wrote into the interpreter's `site-packages/pyHepMC3/__pycache__/` (the `pyHepMC3.so` itself, dated 2026-08-29, is untouched) (the
existing `site-packages/pyHepMC3/pyHepMC3.so` still carries its 2026-08-29
mtime).

Worth stating beside that: **nothing in LiPolGen's suite imports `pyHepMC3`.**
`python/tests/test_hepmc.py` uses the unrelated PyPI package `pyhepmc`.

**Measured, cold, 8 cores** — `time .github/scripts/build_deps.sh <scratch>/deps/install 8`:

**337.28 s (5 min 37 s)**, exit 0. Short by one download only: HepMC3's
tarball was already in the source directory from the aborted first run and was
re-used (≈1–2 s at the measured rate). Stage split, reconstructed from the
artefacts' mtimes (the log markers are unusable — the timestamping pipeline
buffers and clusters them by up to ~15 s):

| stage | wall | dominated by |
|---|---:|---|
| HepMC3 3.3.0 | **57.7 s** | its Python bindings — `libHepMC3.so.4` linked at 14 s, the other 44 s is `pyHepMC3` |
| LHAPDF 6.5.5 | **49.6 s** | |
| the three PDF sets | **≈12 s** | network (62.2 MiB; 5.6 and 8.9 MB/s measured on two full fetches of `CT18NLO.tar.gz`) |
| PYTHIA 8.317 | **≈218 s** | its Python bindings — `libpythia8.so` linked at 70 s, the other ~146 s is the single non-parallelisable `pythia8.so` binding TU |
| **total** | **337.28 s** | |

So **more than half the cold build is Python bindings nothing here imports**;
turning them off would cut it to roughly two minutes. Not done — `env.sh`
plainly intends them to exist — but recorded as the lever.

Cache size, measured on the author's prefix: **448 336 030 B (431 MiB)** on
disk, **119 039 665 B (113.5 MiB)** as `tar | zstd -3` (47 652 839 without the
PDF sets, 71 390 593 for the sets alone). One entry, against a 10 GiB
repository budget.

### E4.3 (a) LiPolGen built and tested against that fresh prefix

The prefix building is not the claim; *LiPolGen being the same program against
it* is. The repository was copied (minus `.git`, `build/`, `.skbuild/`,
`.pytest_cache`) to `<scratch>/ci-clone/LiPolGen` with the §E4.2 prefix as its
sibling `<scratch>/ci-clone/deps` — the exact layout the workflow creates with
`actions/checkout` and `path: LiPolGen`. **`source env.sh` needed no
CI-specific patching**: it resolved `LIPOLGEN_DEPS`, `PYTHIA8DATA`,
`LHAPDF_DATA_PATH`, `LD_LIBRARY_PATH` and `PYTHONPATH` onto the new prefix
itself. `ldd build/lipolgen_tests` confirms the linkage is to the new prefix
(`…/ci-clone/deps/install/lib/lib{LHAPDF,HepMC3.so.4,pythia8}.so`), not the
author's.

| step | result |
|---|---|
| `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release` | **3.08 s**, exit 0 |
| `cmake --build build -j8` | **95.93 s**, **0** `warning\|error` lines, exit 0 |
| `./build/lipolgen_tests` | **400 cases / 17 240 262 assertions / 1 skipped / 0 failed**, 154.81 s |
| `python -m pytest python/tests -q` | **923 passed, 114 skipped, 0 failed**, 128.63 s |
| docs gate | 1145 / 97 / 7 / **0 broken** / 6, and 0 broken in both extra documents |
| SPDX gate | 101 files, 0 missing |

**The doctest tallies are bit-for-bit the author's own prefix's.** A stack
rebuilt from six URLs reproduces, to the digit, the environment this
generator's rtol-1e-12 numbers were measured in.

pytest reads **923 / 114** where the working tree read **925 / 112** on the
day (phase F later added one pytest; the tree reads 926 / 112 now). Both
extra skips were identified, not assumed: `test_release_metadata.py:121` (*no
'origin' remote configured in this checkout* — the copy has no `.git`) and
`test_hfs.py:125` (*polligen is not importable* — `conftest.py`'s fixture
wants a sibling `PolarizedLithiumSim/` checkout). On GitHub the first comes
back (`actions/checkout` creates a `.git` with an `origin`) and the second
does not, so **924 passed / 113 skipped** was the prediction there for this
commit — a prediction, not a measurement. With phase F's added pytest it
becomes **925 / 113**.

### E4.4 (a) `.github/workflows/ci.yml` — and what is NOT claimed about it

**442 lines (398 non-blank), four jobs**, `ubuntu-22.04` throughout (glibc 2.35 — the same
libc the measurements above were taken on and the same one `auditwheel`
derives `manylinux_2_35` from):

1. **`deps`** — `actions/cache@v4` on `${{ github.workspace }}/deps/install`,
   keyed on `ubuntu2204 | py3.11 | hepmc3.3.0 | lhapdf6.5.5 | pythia8317 |
   the three set names | -v1`; on a miss, apt + `build_deps.sh`. Bumping a
   version in the workflow's `env:` block is what invalidates it.
2. **`test`** (`needs: deps`) — `actions/cache/restore` with
   `fail-on-cache-miss: true`, apt, pip, then `source env.sh` **once** with the
   five variables forwarded through `$GITHUB_ENV` (each `run:` is a fresh
   shell, so env.sh stays the single source of truth without being re-sourced
   per step), configure, build, `lipolgen_tests`, pytest, the docs gate, the
   SPDX gate.
3. **`pythia-tier-absent`** (`needs: deps`) — `-DLIPOLGEN_WITH_PYTHIA=OFF`,
   both suites, asserting they **skip** rather than fail. See §E4.5.
4. **`wheel`** (`needs: deps`) — `pip wheel`, `auditwheel show`, `auditwheel
   repair`, and an import from a venv outside the checkout under `env -i`.
   **No `upload-artifact`, no release step, nothing that publishes** — so that
   pushing the workflow is not by itself an act of distribution.

**Locally exercised**: every `run:` block's shell was run by hand in order,
and §E4.2/§E4.3 are jobs 1 and 2 end to end in the workflow's own directory
layout. **Not exercised**: `act` is not installed here, nor is docker or
podman, so nothing parsed this file as GitHub parses it. What a first push
therefore still has to establish is stated at the top of the file itself: that
`actions/cache@v4` returns 114 MiB intact, that the `ubuntu-22.04` image has
every apt package named, the real 4-core timings, and YAML/expression
validity. Risk was reduced where it could be: `${{ env.* }}` appears in
**four** places only (the cache key in each job's `with:`), and `path:`,
`working-directory:` and `python-version:` spell their values out literally
rather than reading workflow `env`, because a wrong context there fails the
whole run at parse time for a reason no local check can catch. The YAML was
parsed with `yaml.safe_load` (4 jobs, 6/12/10/10 steps) — that proves it is
YAML, not that it is a valid workflow.

Both new files carry `SPDX-License-Identifier: GPL-3.0-or-later` on their
first line, but **the SPDX gate's count stays 101, not 103**: its `GROUPS`
globs (E2.3) cover `include/`, `src/`, `tests/`, `python/lipolgen/`,
`validation/*.py` and `python/bindings.cpp`, and `.github/` is in none of
them. Extending the globs would move a number that E2 pinned and that
`STATUS.md` row 13, §E2.3 and §E5 all quote, for no gain this phase — so the
headers are there and the covered set is deliberately unchanged.

**Timing, stated honestly.** First run, cache miss: `deps` **337 s measured at
8 cores**, extrapolated **8–12 min** on a 4-core runner; then `test`
**≈ 6.4 min measured** (3.08 + 95.93 + 154.81 + 128.63 + gates),
extrapolated **9–11 min**, with jobs 2–4 in parallel. Cached run: the `deps`
job is a ~114 MiB restore and the build step is skipped entirely, so the
wall clock is the longest of jobs 2–4.

### E4.5 (a) The tier-absent job found three real defects. All three fixed

The brief asked for "skip-with-reason for the PYTHIA-tier tests if the tier is
absent, **matching what the tests already do**". Measured: **they did not do
it.** `-DLIPOLGEN_WITH_PYTHIA=OFF` did not build at all, and pytest did not
run at all.

1. **`examples/hadronize_example.cpp` and `examples/generate_full.cpp` failed
   to link**: `undefined reference to lipolgen::PythiaBridge::~PythiaBridge()`
   and friends. The `file(GLOB examples/*.cpp)` loop in `CMakeLists.txt`
   builds every example unconditionally; the conditional
   `target_link_libraries` inside it only stops the *library* being added, it
   cannot stop the *example's own calls* into it.
   **Fix**: the loop now `file(READ)`s each example and skips one that
   includes a tier header whose target does not exist, with a `message(STATUS
   "example <name>: NOT BUILT -- it uses the <tier> tier, which is off")`.
   Paired by header (`lipolgen/pythia_bridge.hpp`, `…/hepmc_writer.hpp`,
   `…/lhapdf_sf.hpp`), not by a hardcoded example list. With all three tiers
   on — the shipped default — nothing is skipped and the generated build is
   unchanged.
2. **`tests/test_pythia.cpp` failed to compile**: `:536: error:
   'LIPOLGEN_PYTHIA8_XMLDOC' was not declared in this scope`. It was the only
   PYTHIA-using test file with **no** `#ifdef LIPOLGEN_HAVE_PYTHIA8` anywhere
   in it (`test_t2.cpp` has had one at its line 45 all along).
   **Fix**: the whole file is now behind that macro, exactly as `test_t2.cpp`
   is, with the reason written at the guard. +11 lines at line 18 and 2 at the
   end; the macro is defined in every tier-on build, so this is a no-op there.
3. **pytest ended in a collection ERROR, not a skip**:
   `ERROR collecting python/tests/test_mstw_sf.py … AttributeError: module
   'lipolgen._lipolgen' has no attribute 'pythia8_pdfdata_dir'`, exit 2. That
   module's own docstring says it *"skips wholesale when the PYTHIA 8 tier is
   not built"*, and its module-level `pytestmark` says so too — but
   `needs_grid` calls `_grid_present()` at **import** time, before any skip
   mark can take effect, and `pythia8_pdfdata_dir` is only bound when the tier
   is compiled in.
   **Fix**: `_grid_present()` returns `False` when `HAVE_PYTHIA8` is false,
   with the whole diagnosis in a comment so the guard is not later removed as
   redundant. A tier-on build never reaches that branch.

With all three fixed, measured 2026-09-05:

| | tier ON (shipped default) | tier OFF |
|---|---|---|
| `cmake --build -j8` | — | **93.38 s**, 0 warnings, exit 0, two examples skipped with a stated reason |
| `build/lipolgen_tests` | 400 / 17 240 262 / 1 skipped / 0 failed | **377 / 17 223 404 / 1 skipped / 0 failed** |
| `pytest python/tests -q` | 925 passed / 112 skipped | **874 passed / 163 skipped / 0 failed** |

**925 − 874 = 51 = 163 − 112**: the same 51 tests, skipped instead of run, and
every one of them accounted for by a PYTHIA reason string — 24 *PYTHIA 8 built
without its pdfdata/mstw2008lo.00.dat*, 17 *built without PYTHIA 8*, 4 *no
PYTHIA 8 tier in this build*, 2 *no PYTHIA 8 tier*, 2 *needs the OPTIONAL
PYTHIA 8 tier…*, 2 *MstwSF needs the PYTHIA 8 tier*. The workflow asserts exit
0 and the presence of a reason string, **not** the counts, which move whenever
a test is added.

A trap worth recording, because it produced a false pass first: building the
tier-off tree into a *different* build directory and pointing `PYTHONPATH` at
it **silently tests the wrong extension** — `python/tests/conftest.py` inserts
`<repo>/build/python` at `sys.path[0]` unconditionally. The first tier-off
pytest run therefore reported the tier-**on** tallies (925/112) while
`lg.HAVE_PYTHIA8` was `False` in the same shell. The escape hatch is
`LIPOLGEN_TESTS_USE_INSTALLED=1` (used here); the workflow instead builds into
`build/`, which a fresh checkout does not have.

**Not fixed, disclosed instead**: `-DLIPOLGEN_WITH_HEPMC3=OFF` and
`-DLIPOLGEN_WITH_LHAPDF=OFF` have the same shape of problem —
`tests/test_hepmc.cpp`, `tests/test_pipeline.cpp` and `tests/test_lhapdf.cpp`
include those tiers' headers with no `#ifdef` at all. Neither configuration
was built or measured here. The brief named the PYTHIA tier; these two are
named so the next task can find them, not silently left.

### E4.6 (b) The wheel: `pip wheel`, then `auditwheel repair`

`auditwheel` had never been run on this tree. It is not installed
system-wide, nor is `patchelf`; both came from PyPI into a throwaway venv —
**auditwheel 6.8.2**, and the `patchelf` wheel supplying the **0.19.1**
binary, which auditwheel 6.x needs on `PATH` and does not vendor.

```
LIPOLGEN_DEPS_PREFIX=… pip wheel . --no-deps -w wheelhouse     87.77 s
  -> lipolgen-0.1.0-cp311-cp311-linux_x86_64.whl   1 779 767 B
     36 entries, 5 489 876 B uncompressed
auditwheel show   -> consistent with "linux_x86_64";
                     constrains the platform tag to "manylinux_2_35_x86_64"
auditwheel repair -w repaired …                                 5.35 s
  -> lipolgen-0.1.0-cp311-cp311-manylinux_2_35_x86_64.whl  7 232 971 B
     53 entries, 21 955 667 B uncompressed  (4.06x the raw wheel)
```

What it vendored, into a new `lipolgen.libs/`:

| library | uncompressed | in the wheel | upstream licence |
|---|---:|---:|---|
| `libpythia8-6b6416ce.so` | 13 999 009 | 4 629 974 | GPL-2.0-or-later |
| `libLHAPDF-b6160c34.so` | 1 238 777 | 421 568 | GPL-3.0 |
| `libHepMC3-f2bce5c3.so.4` | 1 175 953 | 398 671 | GPL-3.0 |

and it rewrote four RPATHs to `$ORIGIN:$ORIGIN/../lipolgen.libs`
(`_lipolgen…so`, `libLiPolGen{HepMC,LHAPDF,Pythia}.so`). All 20 `data/vmc`
tables survive. `manylinux_2_35` is glibc 2.35; a wider tag needs a build
inside an older manylinux image, which `auditwheel` says itself and which is
not done here.

**What it does NOT rewrite.** `lipolgen/libLiPolGenCore.so` keeps
`RUNPATH $ORIGIN:/home/cpeng/Projects/polli/deps/install/lib`, because that
library needs nothing from the prefix and so auditwheel never touches it. The
path is inert, but **the builder's home directory is a literal string in the
artefact**, and four more of the same kind survive:
`libLiPolGenPythia.so` and the vendored `libpythia8` both carry
`…/deps/install/share/Pythia8/xmldoc`, the vendored `libLHAPDF` carries
`…/deps/install/share`, the vendored `libHepMC3` carries `…/deps/install/lib`.
The first three are not cosmetic — they are §E4.7.

### E4.7 (b) Does the repaired wheel travel? Measured under a masked prefix

Installed into a clean venv, run from `/tmp` under `env -i` (no
`LD_LIBRARY_PATH`, `PYTHONPATH`, `PYTHIA8DATA` or `LHAPDF_DATA_PATH`), the
repaired wheel imports, reports all three tiers `True`, generates events, and
runs the whole T2 chain (`PYTHIA: 25 ok, 0 failed, 0 retries`, HepMC3 written).

That is **not yet proof**, because `/home/cpeng/Projects/polli/deps/install`
still exists on this machine and the strings in §E4.6 still resolve. So the
same commands were re-run with that directory replaced by an empty tmpfs
inside a `bwrap` user namespace — the nearest thing to "a machine that never
had the deps" available without a second machine, with
`os.path.isdir(prefix) == False` verified in-process.

**Prefix gone, nothing exported:**

| | |
|---|---|
| `import lipolgen`, three tiers `True` | **works** — all three vendored `.so`s load via `$ORIGIN/../lipolgen.libs` |
| `lg.run(channel='inclusive', …)` | **works**, `<x> = 0.027023337580082865` — bit-identical to the unmasked run |
| `lg.run(channel='tagged-alpha', cluster_wave='vmc')` | **works** — the wheel's own `data/vmc` tables found with nothing set |
| `LhapdfSF('CT18NLO')` | **FAILS**: `RuntimeError: Couldn't find required lhapdf.conf system config file`, then LHAPDF's own `terminate` |
| `lipolgen-run … --hadronize` | **FAILS**: `PYTHIA Error in Settings::init: settings file /home/cpeng/…/share/Pythia8/xmldoc/Index.xml not found` → `Abort from Pythia::Pythia: settings unavailable`, exit 1 |

**Prefix gone, the two data trees present at a completely different path**, with
`PYTHIA8DATA` / `LHAPDF_DATA_PATH` / `LIPOLGEN_PYTHIA8_PDFDATA` exported at it:

| | |
|---|---|
| `LhapdfSF('CT18NLO').f2p(0.1, 5.0)` | **0.4238786262338316** — identical to unmasked |
| `lipolgen-run … --hadronize`, 25 events | **works**, and the HepMC3 output is **byte-identical** to the unmasked run (sha256 `68ed6f06…e405ac`, 93 662 B) |

So the honest statement, which is stronger than the tree used to make about
libraries and weaker than it made about data: **`auditwheel repair` solves the
library half completely and cannot touch the data half at all.** PYTHIA's
`xmldoc` (3.7 MiB) and `pdfdata` (52 MiB) and LHAPDF's `lhapdf.conf` + set
store (273 MiB) are read at run time from paths compiled in at build time, and
those paths are the build machine's. A recipient gets the whole pure-C++
generator including the VMC tables, and stops at the first `LhapdfSF` or
`PythiaBridge` until three variables are exported.

### E4.8 (b) The GPL consequence — a statement, not a change

The repaired wheel contains libpythia8 (GPL-2.0-or-later), libLHAPDF (GPL-3.0)
and libHepMC3 (GPL-3.0). Distributing it distributes a combined/linked work of
all three, whose union is **GPL-3.0-or-later** — which is already this
project's licence (`LICENSE`, `CITATION.cff`, and E2's
`SPDX-License-Identifier: GPL-3.0-or-later` on all 101 covered files). So
nothing about LiPolGen's licensing moves; what moves is who carries the
obligation. Whoever publishes such a wheel ships it under GPL-3.0-or-later for
the **combination**, must make corresponding source available for all of it,
may add no further restrictions, and passes the same on to anyone who
redistributes it. The reasoning for choosing GPL-3.0-or-later at all, and why
a permissive licence on LiPolGen's own files would not have avoided any of
this, is `docs/open_items/engineering.md` §C and `OPEN_ITEMS_SOLUTIONS.md`
§12–13.

**Nothing was published.** No wheel left this machine; the workflow has no
upload step of any kind.

### E4.9 A consequence this phase caused, and fixed rather than disclosed

Editing `docs/USAGE.md`'s packaging section (+16 lines) shifted every block
below it, and `docs/PHYSICS_CHANNELS.md` cites seven of them by line range:
the gate went to **8 broken** (7 `MOVED`, 1 `EDITED in place`). Handled the
way the gate's own docstring prescribes rather than by re-recording blindly:

* the one EDITED block is the pip-install section the row itself is about; the
  row was re-read, found still true, **updated** with the measured auditwheel
  results, and its citation widened `docs/USAGE.md:11-50` → `:11-66` (the
  section's true new bounds, 11 to 66, read off the file);
* the seven MOVED citations were rewritten by `--fix` (`:147-484 → :163-500`,
  `:1701-1892 → :1717-1908` twice, `:1920-2073 → :1936-2089`,
  `:1258-1264 → :1274-1280`, `:1284-1292 → :1300-1308` twice);
* `--record-ranges` then wrote **220 fingerprints (1 new, 0 changed, 1
  dropped)** — the same total as before.

Same row, a second drift found while there: its `CMakeLists.txt` line ranges,
which the row itself flags as *read, not machine-verified*. The install-rules
range was recorded as **240–282**; it was **248–292** even before this phase,
and 2026-09-05's tier-aware examples loop moved it again to **279–323**. All
three numbers are now in the row, with the correction. (The RPATH block's
30–55 was re-read and is correct.)

`tests/test_pythia.cpp`'s +11-line guard shifts citations into that file by
12 in total (1 from E2's SPDX line, 11 from here). Exactly **one** citation in
the tree points into it — `docs/code_review_2026-08-29.md:152` cites
`test_pythia.cpp:515` — and it is in one of the 18 dated historical documents
whose drift §E1.5 measured and deliberately did not hand-fix. It is left
alone for consistency with that decision, and counted here so the ledger stays
complete.

### E4.10 What is left to the author

Two rows in `STATUS.md`'s table, **14** (push the workflow) and **15**
(publish a wheel). Neither is a code change and neither was taken here. The
files exist, are exercised as far as this machine allows, and are inert until
someone pushes them.

## E5. Final tallies (this phase, all changes applied)

Measured after E1–E3 and re-measured after E4; where the two differ the E4
figure is given second and the difference explained.

```
$ cmake --build build -j
```
**After E1–E3**, full rebuild (100 source files each gained one new line, so
every target that includes or compiles one of them relinks) — 105.7 s wall,
`grep -icE "warning|error"` on the complete log: **0**; exit 0.
**After E4** — E4 touched `CMakeLists.txt` (a reconfigure) and
`tests/test_pythia.cpp` (one TU plus the `lipolgen_tests` link): that rebuild
was **8.44 s**, 0 warnings; a second `cmake --build build -j` afterwards is a
**0.75 s no-op**, 0 warnings, exit 0.

```
$ build/lipolgen_tests
```
**400 cases / 17 240 262 assertions / 1 skipped / 0 failed** — bit-for-bit
identical to E0's baseline and to every measurement in this run; exit 0
(`Status: SUCCESS!`).

```
$ python -m pytest python/tests -q
```
**925 passed, 112 skipped in 128.60s (0:02:08)** — 909 (E0 baseline) + 16 new
(2 in `test_doc_link_gate.py`, 6 in `test_spdx_headers.py`, 8 in
`test_release_metadata.py`); 0 failed. **After E4: 925 passed, 112 skipped in
128.55s**, 0 failed — E4 added no test and moved none.

```
$ python3 validation/check_physics_channels_links.py
```
```
1145 references checked (strict), 97 ranges, 7 external, 0 broken, 6 allow-listed
19 references checked (strict), 115 ranges, 0 external, 0 broken, 0 allow-listed  [SPIN32_FINITE_GAMMA.md]
0 references checked (strict), 0 ranges, 8 external, 0 broken, 0 allow-listed  [PYTHIA_BRIDGE.md]
```
exit 0. The primary line is byte-for-byte the same as every prior measurement
in this run (same 6 named `allowed` sites, updated only where their own line
numbers moved by the SPDX shift — same reasons, same pins, same file:line
identity otherwise).

```
$ python3 validation/check_spdx_headers.py
```
**101 files checked, all carry 'SPDX-License-Identifier: GPL-3.0-or-later'**;
exit 0.

`git status --porcelain -- validation/reference/` — **empty**, after E4 as
after E3: **no reference JSON moved**, no rtol-1e-12 gate touched.
`git diff --stat -- pyproject.toml` — **empty**: the version single source of
truth was read, not edited.

`CMakeLists.txt`, empty through E3, is **+31 lines** after E4 — the
tier-aware examples loop of §E4.5, which is inert with all three tiers on
(the shipped default) and changes nothing the default build generates. The
C++, pytest, docs-gate and SPDX numbers above are the proof of that, all four
unmoved.

Files changed this phase: **5 new top-level / near-top-level** — `AUTHORS`,
`CITATION.cff` (E2), `.github/workflows/ci.yml`, `.github/scripts/build_deps.sh`
(E4), and `docs/PACKAGING.md` (E4); **3 new** under `python/tests`/`validation`
(`test_release_metadata.py`, `test_spdx_headers.py`, `check_spdx_headers.py`,
E2); 1 new report (this file). **~123 modified**: 100 source files gaining an
SPDX line, `README.md`, `docs/USAGE.md`, `docs/PHYSICS_CHANNELS.md`,
`docs/T2_CHAIN.md`, `docs/PYTHIA_BRIDGE.md`, `docs/DEVELOPMENT_PLAN.md`,
`docs/theory/SPIN32_FINITE_GAMMA.md`, `STATUS.md`,
`validation/check_physics_channels_links.py`,
`validation/physics_channels_ranges.json`,
`python/tests/test_doc_link_gate.py`, and — E4's three fixes —
`CMakeLists.txt`, `tests/test_pythia.cpp`, `python/tests/test_mstw_sf.py`.

**Nothing was committed, pushed, uploaded or published.** No remote was
written to at any point; the only network use was reading (six upstream
tarball URLs probed and fetched, and `pip` fetching `auditwheel`/`patchelf`).
The supervising session commits.

### E4.7 Disclosures the adversarial review asked for (2026-09-05)

- **Network.** This phase's brief allowed network only for `pip` installs into
  a venv. E4 exceeded that: `build_deps.sh` downloaded ~106 MB of upstream
  tarballs and PDF sets (hepmc.web.cern.ch, lhapdf.hepforge.org,
  lhapdfsets.web.cern.ch, pythia.org; tarball mtimes 07:41–07:44) into the
  scratch prefix, and ~14 URLs were probed with `curl -sSI` (six 200s; seven
  `pythia.org/download/...` 404s — the live path is `pythia.org/releases/`).
  All of it was read-only fetching; nothing was uploaded, posted or published.
  Recorded because the rule was stated and was not kept.
- **An omitted failed run.** The FIRST staged job-2 pytest (07:57) was
  `5 failed, 918 passed, 114 skipped` — all five in `test_doc_link_gate.py`,
  from 8 `docs/USAGE.md` citations broken by E4's own mid-flight edit of that
  file. `docs/PHYSICS_CHANNELS.md` was re-copied into the clone at 08:00:31 and
  the re-run at 08:02 gave the `923 passed / 114 skipped` reported above. Only
  the second run had been reported; both are the record.
- **The `ci.yml` header** said "Nothing in this repository has been pushed to a
  remote". False as written — `origin/master` exists at `a94fd6e` (last push
  2026-09-03). Corrected to what is true: none of this run's commits and no
  `.github/` directory have been pushed.
