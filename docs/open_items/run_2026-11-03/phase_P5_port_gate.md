<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-03 — phase P5: the VMC port-gate block (decision A13) — first pass STOPPED at the additivity check, second pass landed by splice

Date of this record: 2026-11-06. Decision A13 (`DECISIONS.md`) is the
supervisor's choice under the author's delegation. The phase STOPPED at the
additivity check of the brief's design point (c). Nothing of A13 is in the
tree. The work and every measurement are kept in the phase's scratch
directory (`run_1103/p5/`, named below by file).

## 0. What happened, in one paragraph

The dump script `validation/dump_polligen_reference.py` was extended with
A13's reduced `channels.li6_alpha_vmc` block. The code was appended below
`main()` so no cited line moved: `build_tagged_li6_alpha_vmc()`, a
`build_tagged()` that wraps the Hulthén one, and a `README_ADDENDA["tagged.json"]`
paragraph. The tolerance was declared before any comparison (§2). The
script was run in place with the tree's `env.sh`. Its three preflights
passed and the new block came out as designed. The run was NOT additive:
it also rewrote 3170 numbers of six reference JSONs, `tagged.json` among
them, at the rounding level, and one manifest string. Under the brief's
rule (c) every file was restored from `git show HEAD:<path>`, what moved is
recorded in §1, and the phase stopped. Its residual against LiPolGen,
measured read-only from the scratch copy, is 5.389e−9 (vector) and
3.196e−9 (tensor). Both are inside the declared 1e−8. So the gate itself is
sound, and only the landing failed. The dump-script change was also taken
back out of the tree and saved as a patch. No C++ or Python gate was
written. `STATUS.md`'s P5 row reads STOPPED.

## 1. What moved in the dump (not landed; all restored)

Run: `source env.sh && PYTHONDONTWRITEBYTECODE=1 python3 -B validation/dump_polligen_reference.py`
(exit 0; log `dump.log`). The sibling stayed at `4726812` with a clean
status and no new file. Compared entry by entry against `HEAD` (`cmp.py`,
`cmp2.py`; full diffs `dump_reference.diff`, `dump_readme.diff`; the dumped
files in `dumped/`):

| file | numbers changed | max abs diff | worst relative | where |
|---|---|---|---|---|
| `tagged.json` | 51 (+ the new block) | 4.4e−16 | 2.2e−16 (≤ 1 ULP) | `boost_spectator` of all three channels |
| `bookkeeping.json` | 2 | 6.9e−18 | 2.8e−16 (2 ULP) | `helicity_flip_j32_maxent_anchor` populations |
| `coherent.json` | 9 | 2.8e−17 | 1 ULP | `tag_acceptance*`, `fragment_rigidity` |
| `spectator.json` | 28 | 4.3e−19 | 1 ULP | `boost_fragment` of six channels |
| `spin.json` | 202 | 1.9e−16 | (values near 0) | `population_sets` of j1 and j3_2 |
| `xsec.json` | 2878 | 1.4e−9 | 9.4e−13 among values > 1e−6; 3.8 on 1e−19 cancellation residues | every kernel block |
| `_manifest.json` | 0 numbers | – | – | `could_not_call[0].reason`: `tagged.py:299-302` → `tagged.py:849-852` |
| `beams.json`, `b1_default_li6.json` | 0 | – | – | – |
| `validation/README.md` | +24 lines, 0 removed | – | – | the A13 paragraph only (additive) |

The manifest string is an older mismatch. The committed manifest was not
regenerated after the script's citation text changed, so it is not today's
drift. The numeric moves are rounding-level and come from outside both
repositories. The sibling is unchanged, and numpy 2.2.6 has been installed
since 2026-06-12. `dpkg.log` shows `libc6` 2.35-0ubuntu3.13 → .15 and
`libpython3.10` ~22.04.15 → .18 upgraded on 2026-10-26, after the last
in-tree dump (2026-09-23, whose run reproduced `xsec.json` bit for bit).
That is the likely cause, but it is NOT proven. In `xsec.json` the moves
reach 9.4e−13 relative on values above 1e−6, close to the rtol 1e−12 the
port gates apply. No reference number moved in the tree.

## 2. What was verified and how

- **The declared tolerance, before the comparison** (`DECLARED.txt`,
  2026-11-06T21:04:39−06:00, written before the dump ran):
  `vector_dilution` and `tensor_dilution` at relative 1e−8, never widened;
  `P_D` and `j_ion`, `s_struck`, `s_spec`, `s_channel` at rtol 1e−12; `label`
  not compared.
- **The dumped block** (`dumped/tagged.json`): `P_D` 0.019354933183319687,
  `vector_dilution` 0.970966091698428, `tensor_dilution` 0.9825758754851597,
  `j_ion` 1, `s_struck` 1, `s_spec` 0, `s_channel` 1, label "6Li alpha-tag
  (embedded d, VMC alpha-d)", `port_gate_rtol` 1e−08 with its reason and
  `scope`.
- **The residual** (`residual.log`), against the installed module's
  `TaggedModel(li6_alpha_channel(source=ClusterWaveSource.VmcAV18))`:
  vector 0.970966086466056 against 0.970966091698428, **5.389e−9**; tensor
  0.982575872345055 against 0.98257587548516, **3.196e−9**; `P_D` and
  the four quantum numbers 0 difference. The labels differ: LiPolGen's is
  "… VMC AV18+UX)". This reproduces the D6 design's 5.4e−9 / 3.2e−9
  (`../run_2026-10-27/plan_D6_cluster_wave_default.md` §4.5).
- **An additive candidate, NOT in the tree** (`tagged.candidate.json`):
  `HEAD`'s `tagged.json` round-trips byte-identically through
  `json.dumps(indent=1, sort_keys=True)`. With only the dumped
  `li6_alpha_vmc` block inserted, the diff against `HEAD` is 15 added lines
  and 0 removed (`tagged.candidate.diff`).
- **Restoration**: `git status` after the restore shows no modified file
  under `validation/`. File modes are unchanged (664).

## 3. Files

In the tree: this record (new) and `STATUS.md` (the P5 row only).
`validation/dump_polligen_reference.py`, `validation/reference/*.json` and
`validation/README.md` are byte-identical to `HEAD`. In scratch:
`dump_script_A13.patch`, which `git apply --check` accepts on `HEAD`, the
full script `dump_polligen_reference.A13.py`, `tagged.candidate.json` and
the logs above.

## 4. Gates

| gate | result | wall |
|---|---|---|
| the three dump preflights (`_preflight_sd_phase`, `_preflight_delta_convention`, `_preflight_cluster_wave`) | passed (printed: coefficient 1, the LOI's; defaults Hulthén) | the dump, < 1 min |
| additivity of the dump (design point c) | **FAILED**: 3170 numbers in 6 files + 1 manifest string; restored | – |
| residual against the declared 1e−8 | 5.389e−9 / 3.196e−9: within | < 1 min |
| C++ build, `lipolgen_tests`, pytest, bench, report | NOT RUN: no code or reference change in the tree | – |
| `check_spdx_headers.py` | 136 files, all carry the line | 0.04 s |
| `check_physics_channels_links.py` (strict) | 1243 / 95 / 7 / 0 broken / 24 | 2.4 s |
| `check_physics_channels_links.py --records-only` | 1432 citations in 87 records, 0 broken | 0.3 s |
| `pytest python/tests/test_doc_link_gate.py` | 112 passed | 24.5 s |

## 5. Open points, for the supervisor

1. **How A13 lands, given an environment that no longer reproduces the
   committed tables.** (i) Splice: commit `tagged.candidate.json`, i.e.
   `HEAD` plus the block the dump produced, purely additive. Then apply
   `dump_script_A13.patch`, and write the three gates the brief names at
   `kVmcPortRtol = 1e-8`. The Hulthén loops in `tests/test_tagged.cpp`
   (three cases) and `tests/test_cluster.cpp` (`Wave::radial`) iterate over
   EVERY channel key and `model_by_key` falls through to the deuteron, so
   they need a line-neutral skip of `li6_alpha_vmc`. The file would then no
   longer be exactly what one dump run writes in today's environment, and
   the record would have to say so. (ii) Re-dump everything and accept the
   rounding moves: this is a reference-number move, which A13 does not
   authorize. (iii) Find and fix the environment drift first.
   Recommendation: (i). The block's numbers come from the sibling through
   the dump's own code path. The Hulthén blocks stay the committed bytes,
   and the drift is filed as its own item.
2. **The environment drift itself** is a standing fact for every future
   dump: rounding-level, up to 9.4e−13 relative in `xsec.json` against
   rtol 1e−12 gates. The likely cause is the 2026-10-26 libc6 / libpython
   upgrade, not proven.
3. **The stale manifest string** (`tagged.py:299-302` against the script's
   `849-852`) predates this run. Any future dump rewrites it.

---

# Second pass, 2026-11-06: A13 landed by splice

The supervisor took open point 1's recommendation (i). It is recorded as
the supervisor's dated note under A13 in `DECISIONS.md`. Sections 0–5
above are the first pass, kept as written. What follows is the landing.

## 0. What changed, in one paragraph

`validation/reference/tagged.json` is now `HEAD`'s file with the dumped
`channels.li6_alpha_vmc` block inserted byte for byte (+15/−0). No other
byte of any reference JSON changed. `_manifest.json` is untouched, because
its test, `python/tests/test_kernel_reference.py`
`test_reference_manifest_lists_exactly_the_tables_on_disk`, lists files,
not blocks. `validation/README.md` gains exactly the A13 paragraph the
dump wrote (+24/−0, at the end of the tagged.json section), and the file
is now byte-identical to the README the dump produced. The dump script
carries the block (the first pass's patch, applied), so a future dump
writes it. Three gates read the block:
- `tests/test_tagged.cpp`: a new case at the end of the file, with a named
  `kVmcPortRtol = 1e-8` and `A13 (2026-11-03)` comments;
- `validation/repin_tagged_from_lipolgen.py`: `check_vmc_block`, check-only;
- `python/tests/test_port_gate_li6_vmc.py`: new, with the SPDX line.

The three Hulthén loops of `tests/test_tagged.cpp` and the `Wave::radial`
loop of `tests/test_cluster.cpp` skip the key. `model_by_key` now throws on
an unknown key instead of falling through to the deuteron. All five of
those edits are line-neutral. No physics default, constant or existing
tolerance moved. The report does not change, because no harness reads this
block, and it was not regenerated.

## 1. The block's numbers and their origin

Every value was produced by `build_tagged_li6_alpha_vmc()` in the in-place
dump of 2026-11-06, run against PolarizedLithiumSim `4726812` (first pass
§1; the dump's own file is `dumped/tagged.json` in the scratch directory).
Measured against `TaggedModel(li6_alpha_channel(source=VmcAV18))`:

| entry | polligen (the file) | LiPolGen | relative | declared |
|---|---|---|---|---|
| `model.vector_dilution` | 0.970966091698428 | 0.97096608646605587 | 5.389e−9 | 1e−8 |
| `model.tensor_dilution` | 0.9825758754851597 | 0.98257587234505528 | 3.196e−9 | 1e−8 |
| `P_D` | 0.019354933183319687 | `VMC_P_D_LI6`, same | 0 | 1e−12 |
| `j_ion`, `s_struck`, `s_spec`, `s_channel` | 1, 1, 0, 1 | same | 0 | 1e−12 |
| `label` | "6Li alpha-tag (embedded d, VMC alpha-d)" | "… VMC AV18+UX)" | not compared | – |
| `port_gate_rtol` | 1e−08 | `kVmcPortRtol` | asserted equal | – |

The C++ case prints `vector rel 5.38883e-09, tensor rel 3.19579e-09`.
Nothing that existed before moved.

## 2. What was verified and how

- **The splice is additive.** `git diff --numstat validation/reference/`
  shows `tagged.json` alone at 15/0. A per-file `git diff --quiet HEAD` over
  the eight other JSONs is silent. The inserted block equals the dumped one
  both as parsed JSON and as text.
- **The README matches the dump.** `cmp` of `validation/README.md` against
  the README the dump wrote shows them identical. The one line citation of
  `validation/README.md` in `docs/` (a D7 plan) points above the inserted
  paragraph, which goes at the end of the tagged.json section.
- **The new gate and the old loops.** The new C++ case passes with 11
  assertions. Every Hulthén case still passes: the tagged and cluster
  subset is 65 cases, 1 169 384 assertions.
- **The check-only script.** `repin_tagged_from_lipolgen.py` passes the
  three Hulthén blocks with the same worst ratios as before, plus the VMC
  block.
- **Citation ranges.** The `tests/test_tagged.cpp` edits sit above the
  fingerprinted ⁷Li ⟨P₂⟩ polarimeter case and are line-neutral, as are the
  edits to the cluster radial-forms bullet of `docs/CONVENTIONS.md` and the
  port-gate bullet of `README.md`. No range moved: the three citation modes
  are green.

## 3. Files

- `validation/reference/tagged.json` (+15)
- `validation/README.md` (+24)
- `validation/dump_polligen_reference.py` (+91, appended below `main()`)
- `validation/repin_tagged_from_lipolgen.py` (`check_vmc_block`,
  `main_with_vmc`, one docstring paragraph)
- `tests/test_tagged.cpp` (the new case and `kVmcPortRtol` appended; four
  line-neutral edits)
- `tests/test_cluster.cpp` (one line-neutral skip)
- `python/tests/test_port_gate_li6_vmc.py` (new: 4 tests)
- `docs/CONVENTIONS.md` (the cluster radial-forms bullet: one clause, A13)
- `README.md` (the port-gate bullet: one clause, A13)
- this record, and `STATUS.md` (the P5 row)

## 4. Gates (final tree)

| gate | result | wall |
|---|---|---|
| `cmake --build build -j2` | 0 warnings | 8.7 s |
| `lipolgen_tests` new case / tagged and cluster subset | 1/1 (11 assertions) / 65/65 (1 169 384) | 6.7 s (subset) |
| `lipolgen_tests` full | **435/435**, 17 244 594 assertions | 7 m 01 s |
| pytest run 1 (14 files, `test_b1_model` … `test_bench_t1_r_world`) | 358 passed, 3 skipped | 2 m 49 s |
| pytest run 2 (15 files, `test_bench_t4_t6` … `test_meta_provenance`) | 136 + 1 skip; 137; `test_knob_provenance` 249 + 90 skip and 277 + 61 skip; 12; 37 | 1 m 36 s + 28.8 s + 5 m 39 s + 4 m 07 s + 18.6 s + 20.8 s |
| pytest run 3 (15 files, `test_module` … `test_throughput`) | 317 passed | 2 m 41 s |
| pytest total, 44 files, disjoint | **1523 passed, 155 skipped, 0 failed** (P4b's 1519 + the 4 new) | ≈ 18 m 40 s |
| `repin_tagged_from_lipolgen.py` | 3 Hulthén blocks + the VMC block agree; nothing written | 0.2 s |
| citation strict | 1243 / 95 / 7 / 0 broken / 24 | 3.1 s |
| citation `--records-only` | 1432 in 87 records, 0 broken | 0.3 s |
| citation `--audit-ranges` | 209 range citations on 179 blocks, 0 REFUSED | 3.0 s |
| `check_spdx_headers.py` | 136 files, all carry the line | < 0.1 s |
| `test_doc_link_gate.py` | 112 passed | 24.9 s |
| bench subset, report | not run: no harness reads the block; the report was not regenerated | – |
| bench subset (`LIPOLGEN_BENCH=1 … -k bench`), run by the supervisor after the hand-back | 368 passed, 1 skipped, 0 failed (unchanged from P4b) | 6 m 26 s |

Run 2 was split further because its first form hit `timeout 590`. A
590 s attempt on 15 files ran out, as did a second on 8 files and a third
on `test_knob_provenance` plus `test_li6_unpol_band`. All three were
discarded with no result. The final sub-runs are disjoint and cover
`test_knob_provenance`'s 677 collected tests exactly, as two node-id halves
of 339 and 338.

## 5. Open points

1. **For the author: the reference dump does not reproduce on today's
   machine.** The in-place dump moves 3170 numbers at the rounding level,
   up to 9.4e−13 relative in `xsec.json` against the rtol 1e−12 port gates
   (first pass §1). The likely cause is the 2026-10-26 libc6 / libpython
   upgrade, not proven. The committed tables still pass every port gate.
2. **The stale manifest string** (`tagged.py:299-302` against the script's
   `849-852`) is left as it is. The next real dump rewrites it.
3. **The time zone in `DECISIONS.md`.** The supervisor's note under A13
   says the tolerance was written "at 21:04:39 UTC". `DECLARED.txt`
   records 2026-11-06T21:04:39−06:00, which is local time (CST). It was
   left unedited, as instructed.
4. **Run 2's pytest list** does not fit in one 590 s call on this machine:
   `test_knob_provenance.py` alone takes about 9 min 45 s. Later phases
   should split it by node id, which `run_ids.py` in this phase's scratch
   directory does.
5. **The splice is not reproducible from a single run.** The block lives
   in a file that one dump run in today's environment would not reproduce
   byte for byte outside it (item 1). The next full re-dump in a
   reproducing environment should leave `tagged.json` unchanged.
