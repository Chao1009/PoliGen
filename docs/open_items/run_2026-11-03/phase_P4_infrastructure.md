<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-03 — phase P4: infrastructure and provenance (decisions A17, A11, A12, A10, A9; registry rows for A2, A8, A14)

Worker P4, 2026-11-06, on `claude/run-2026-10-26` at `5cd19d2` (phases P1–P3
committed; this phase uncommitted, the supervising session commits).
Decisions: `DECISIONS.md` A2, A8, A9, A10, A11, A12, A14, A17 — the
SUPERVISOR's choices under the author's delegation ("Go with your
recommendations for the 17 items, just record the choice and why"), not the
author's words. Every number below was measured in this phase with
`source env.sh` on this machine (4 cores, one long gate at a time); scripts
and logs are in the session scratch `run_1103/p4/` (`ab_runs.sh`,
`ab_compare.py`, `ab_compare.txt`, `ab_explicit_defaults.txt`,
`a12_measure.txt`, `usage_knob_example.log`, `build_nolhapdf.log`,
`build_off_*.log`, `flip1_*.log`, `final_*.log`, `report*.log`,
`result.json`). **pre** = `5cd19d2`'s build, taken before any edit;
**post** = this tree.

## 0. What changed, in one paragraph

**A17**: one rule for an explicit DEFAULT on a knob the run never reads —
accepted, its knob row `not-read` with a label naming why, the `meta` key
carrying the label; a non-default value keeps its refusal. Applied to
`cluster_wave` off the tags (`make_config` / the CLI now accept `vmc` there
and refuse `hulthen`; `_cluster_wave_on_tagged_only` in
python/lipolgen/__init__.py; `validate()` unchanged), to `rc_tail_model` and
every other band / tail knob of `Pipeline::knob_provenance` at its default
where the band or the tail does not apply (`band_row`, `tail_row`, the
`!band` / `!tail` branches of `rc_scope`, the quasi-elastic rows,
`rc_sp_tensor_scale`, `rc_m_lepton`; one reason suffix `kA17`), and
`meta["rc_tail_model"]` / `meta["rc_c0_shape"]` now go through the table.
CONVENTIONS "Refuse or label" states the rule once (line-neutral) and names
the rows it does not yet reach (§5). **A11**: `validate()` refuses
|`cluster_vmc_mc_sigma`| > 3, naming the band and the entry. **A12**: the
ceiling stays; the CLI banner prints one line from the same call as
`meta["fsi_clipped_grid_fraction"]` (`fsi_clip_line`; `GlauberFsiOptions.w_max`
bound); registry row 28 appended. **A10**: `lipolgen::run_provenance(p, ctx)`
(pipeline.hpp) is the ONE list of the run's provenance strings; the npz
`meta` fills those keys from it (python/bindings.cpp, line-neutral), and
`HepMC3Writer::set_run_attributes` writes it as `lipolgen.<key>` `GenRunInfo`
string attributes; `lipolgen-run --hepmc`, `Pipeline.write_hepmc(context=)`
and the generate_tagged / _coherent / _full examples pass it (the Python path
adds the two `*_sf_extrapolation` keys from `sf_grid_extras`, the function
`meta` uses); docs/HEPMC3_CONVENTION.md gets an appended section and its
line 232 re-worded line-neutrally; T6 gets clause C24 and a re-pinned sha256.
**A9**: a missing `lhapdf-config` is a configure-time FATAL_ERROR unless
`-DLIPOLGEN_WITH_LHAPDF=OFF`. **A2, A8, A14**: records only. No physics
default, constant, tolerance or reference number moved; every event stream is
bit-identical pre → post (§2).

## 1. Every number or label that moved, before → after

**A17 — knob rows (status at the default; values unchanged).**
- `cluster_wave` on inclusive (6Li, 7Li, d) and coherent: `refused` →
  `not-read`, label "not read on channel <channel>" (the value stays the
  label, never a bare `vmc`); `meta["cluster_wave"]` unchanged (it was the
  label already).
- On `coherent --rc tensor-band` (band and tail absent): 19 rc rows
  `refused` → `not-read` with label "not read on channel coherent-6Li"
  (`rc_delta_low_x`, `rc_delta_high_x`, `rc_x_low`, `rc_x_high`,
  `rc_a_transfer_frac`, `rc_band_tau_max`, `rc_scope`, `rc_fq_scale`,
  `rc_tail_tensor_scale`, `rc_c0_shape`, `rc_qe_suppression`,
  `rc_qe_tensor_scale`, `rc_qe_kf_gev`, `rc_tail_model`,
  `rc_sp_tensor_scale`, `rc_with_qe_tail`, `rc_n_eta`, `rc_tail_max`,
  `rc_m_lepton`); on `tagged-6Li-alpha --rc tensor-band` the 12 tail rows of
  that list, label "not read on channel tagged-6Li-alpha".
- `meta` keys: `rc_tail_model` "polrad-full" → "not read on channel …",
  `rc_c0_shape` "ho" → the label, `rc_scope` "tensor-rate" → the label, on
  those runs (coherent / tagged with `--rc tensor-band`); everywhere the
  tail runs they are unchanged.
- `make_config(channel="inclusive"|"coherent", cluster_wave="vmc")` and
  `lipolgen-run --channel inclusive|coherent --cluster-wave vmc`: refused →
  accepted (exit 1 → 0); `hulthen` there stays refused, message now ending
  "(… D6; the default is accepted and labelled, A17)".
- USAGE's knob-provenance example banner (re-run 2026-11-06, `lipolgen-run --events 200
  --seed 5 --pol-sf toy`): not-read defaults 44 → 45, refused axes 7 → 6
  (`cluster_wave` moved).
- Pins re-taken with `A17 (2026-11-03)` comments:
  `test_an_explicit_default_on_an_unread_knob_is_handled_three_ways` →
  `…_handled_one_way` (three cases plus a band and a tail sub-knob, each
  hash-identical to the unnamed run); the matrix gate's refused branch
  accepts a `not-read` base row only when it is AT ITS DEFAULT and its reason
  carries the A17 sentence (`A17_LABELLED`) — the other value is still
  required to be refused, so nothing is loosened;
  test_cli.py (`test_d6_…_is_one_line`: `vmc` exit 0; the inclusive row's
  status "refused" → "not-read"); test_module.py
  (`test_c5_5b_…`: only `hulthen` raises off the tags); tests/test_pipeline.cpp
  (the off-tag `cluster_wave` row `Refused` → `NotRead`);
  validation/price_defaults.py (three `cluster_wave = vmc` probes "refused" →
  "not-read", which `_probe` checks bit-identical to the default).

**A11.** `cluster_vmc_mc_sigma` = 1000 (and any |σ| > 3): built (α–d P_D
0.133, I6 §12 V-8) → refused: "PipelineConfig: cluster_vmc_mc_sigma = 1000 is
outside the band it names -- a number of the tables' printed 1-sigma;
|sigma| <= 3, … (A11, 2026-11-03 …)". ±3 build; NaN / ±inf refused as before.

**A12.** New banner line on every FSI run (no number moved). Measured
2026-11-06 on the shipped default (CT18NLO / NNPDFpol1.1, 40 mb, 2000 events,
seed 1234; `a12_measure.txt`), banner = meta in each case:

| channel | glauber-cluster | glauber-nucleon |
|---|---|---|
| tagged-7Li-alpha (`helicity-flip`) | 4.41 % (0.044057) | 9.52 % (0.095212) |
| tagged-6Li-alpha | 0.70 % (0.0070213) | 5.25 % (0.052506) |
| tagged-d-p | 1.35 % (0.013464) | 1.35 % (0.013464) |

The ⁷Li and ⁶Li values equal I6 §9 item 3's (toy SFs, 2026-10-29); the d
control equals PHYSICS_CHANNELS row 347's 1.35 % (I6 §9 item 3 quotes the
Hulthén 1.45 %). Registry: 27 → 28 rows (OPEN_ITEMS_SOLUTIONS.md's two "27
rows" sentences carry "(28 since A12)").

**A10.** HepMC3 files written by `lipolgen-run`, `Pipeline.write_hepmc` and
three examples gain 11–19 `A lipolgen.<key> <value>` lines in the run block
(§2); nothing else in them moved. T6: `DOC_SHA256_AT_AUTHORING`
56d3a91a9bdf789b… → 20de1d3b3cbc03a7… (the old digest kept in the
`A10 (2026-11-03)` comment); clauses 23 → 24 (C24, doc lines 255–284, anchor
line 262 "`lipolgen.<key>` holds exactly the string the run's npz
`meta[<key>]` holds"); C23's anchor (line 230) unchanged; the per-run argv
the row prints gains `--npz <file>.npz` (C24's reference, written by the same
CLI run). test_bench_t4_t6.py: `C24.files == 8` added; `_violations` passes
the npz meta; two tamper cases for C24.

**A9.** Configure without `lhapdf-config`: WARNING + tier off → FATAL_ERROR
(exit 1), message quoted in §2.

## 2. What was verified and how

**Event streams, pre → post** (`ab_runs.sh`, 300 events seed 11; 22 CLI runs
over the seven channels, `--rc tensor-band`, the explicit defaults and both
FSI variants on the three tags, plus the six example programs;
`ab_compare.txt`): every npz array of every run that built pre is
byte-identical post (dtype and bytes); `meta` moved only on the coherent /
tagged `--rc tensor-band` runs (`rc_tail_model`, `rc_c0_shape`, `rc_scope`
to the label) and in knob-row statuses as in §1. Every HepMC3 file differs
from pre ONLY by added `A lipolgen.*` lines (CLI runs 13–19 lines, examples
generate_tagged / _coherent 11, _full 13); generate_inclusive,
write_hepmc_example and hadronize_example are byte-identical (they hold no
`Pipeline` and pass no map).

**A17 is label-only** (`ab_explicit_defaults.txt`): the accepted explicit
defaults `inclusive --cluster-wave vmc`, `coherent --cluster-wave vmc`,
`tagged-6Li-alpha --rc tensor-band --rc-tail-model polrad-full`,
`coherent --rc tensor-band --rc-fq-scale 1` and `--rc-delta-low-x 0.30`
against the same runs without the flag: arrays identical, `meta` identical,
knob tables identical, HepMC3 files identical. The one-way pytest asserts the
same through the matrix's output hash.

**A10 round trips**: C++ `A10: run provenance attributes round-trip through
ReaderAscii` (values byte for byte incl. labels with spaces, both events share
the one `GenRunInfo`, the weight block 1 + N, no attribute without a map,
`set_run_attributes` after the first write throws); pytest
`test_a10_a_cli_file_carries_the_npz_meta_as_run_attributes` (pyhepmc; every
attribute equals the npz meta of the same run, on the default tag, an RC tag
and an inclusive run with `--cluster-wave vmc`); T6's C24 on all eight files.
T15/T15b assert the weight names / values (1 + N, no `rc_*`), which still
holds; the document's "byte-identical to one written before `rc.hpp`" was no
longer true of a file carrying the attributes and now reads "its weight block
is exactly the pre-`rc.hpp` one (since A10 the file also carries the
`lipolgen.*` run attributes, last section)" — same line, C23's anchor on 230
standing; USAGE §7b and OPEN_ITEMS_SOLUTIONS §9 say the same.

**A9, out of tree** (scratch, `./build` untouched): (i) `env -i
PATH=/usr/bin:/bin cmake -S . -B build_nolhapdf -DLIPOLGEN_DEPS_PREFIX=/nonexistent`
(no `lhapdf-config` on that PATH): exit 1 in 0.2 s, "CMake Error at
CMakeLists.txt (the LHAPDF block's message): lhapdf-config not found (searched PATH and
/nonexistent/bin). The SHIPPED structure-function defaults (D7: CT18NLO /
NNPDFpol1.1) live in the LHAPDF tier, so a build without it refuses every
DEFAULT run. Either point -DLIPOLGEN_DEPS_PREFIX=<prefix> at an LHAPDF 6
install (or put its lhapdf-config on PATH), or pass -DLIPOLGEN_WITH_LHAPDF=OFF
explicitly for the core-only build, in which only runs naming --unpol-sf toy
--pol-sf toy (C++: lipolgen::toy_config()) work. (A9, 2026-11-03)". (ii)
`cmake -S . -B build_off -DLIPOLGEN_WITH_LHAPDF=OFF`: configures (exit 0, 1 s; the four `generate_*` examples "NOT BUILT -- it uses the LHAPDF tier, which is off"), `cmake --build build_off -j2` builds (exit 0, 137 s, 0 warning lines), and its `lipolgen_tests` pass 417 / 417 test cases, 17 243 193 assertions (320 s; the regression lens had 411 / 411 after D7).
.github/workflows/ci.yml read, not changed: the `test` and
`pythia-tier-absent` jobs configure against the deps prefix the `deps` job
built (its "what the cache now holds" step runs `$P/bin/lhapdf-config
--version`), the wheel job sets `LIPOLGEN_DEPS_PREFIX` to it — CI always has
the tier.

**A8** — looked for an "exit 0" statement about a missing verdict set in
docs/*.md, docs/benchmarking/*.md, validation/benchmarks/README.md and every
validation/benchmarks/*.py (run records excluded): none says a harness whose
VERDICT configuration needs a missing set exits 0; README convention 6
states the exit-2 rule. Nothing edited. Seen and left: t5_delta_moment_bound's
`EnvironmentBlocked` (a set present but with a missing member file or not the
vendored version → BLOCKED, exit 0; a set not installed → exit 2) — §5.

**A2, A14** — dated notes, the supervisor's choice quoted as such: registry
row 3 and PHYSICS_CHANNELS rows 200 (`--r-source`) and 210 (EPPS21) for A2;
PHYSICS_CHANNELS row 387 (the clause "which is why it is the default rather
than a band edge" now says "on `t-peak` only" with the A14 facts),
OPEN_ITEMS_SOLUTIONS §9 and registry row 17 for A14; A11 on registry row 11
and I6 §9 item 11; A17 on registry row 11, PHYSICS_CHANNELS row 270 and USAGE.

## 3. Files

C++: src/core/pipeline.cpp (knob rows, `kA17`, `tail_scope`, the A11 bound —
line-neutral; `run_provenance` appended), include/lipolgen/pipeline.hpp
(`KnobStatus::NotRead` doc, `cluster_vmc_mc_sigma` doc — line-neutral;
`run_provenance` declared at the end), include/lipolgen/hepmc_writer.hpp
(`set_run_attributes`, after every cited line), src/hepmc/hepmc_writer.cpp
(two in-place lines, one comment, the definition appended),
python/bindings.cpp (`meta` provenance keys from `run_provenance`,
`hepmc_run_attributes`, `run_provenance` / `hepmc_run_attributes` /
`HepMC3Writer.set_run_attributes` / `GlauberFsiOptions.w_max` bound,
`Pipeline.write_hepmc(context=)` — line-neutral above the appended helper),
examples/generate_tagged.cpp, generate_coherent.cpp, generate_full.cpp (one
line each), CMakeLists.txt (A9, line-neutral). Python:
python/lipolgen/__init__.py, python/lipolgen/cli.py (`--cluster-wave` help,
the FSI banner call, the `--hepmc` paths — line-neutral; `fsi_clip_line`
appended). Validation: validation/price_defaults.py,
validation/benchmarks/t6_hepmc3_convention.py, validation/physics_channels_ranges.json
(USAGE 1978–2169 re-recorded after the §7b sentence was re-read).
Tests: tests/test_pipeline.cpp, tests/test_hepmc.cpp,
python/tests/test_knob_provenance.py, test_cli.py, test_module.py,
test_hepmc.py, test_bench_t4_t6.py. Documents: docs/CONVENTIONS.md,
docs/USAGE.md, docs/PHYSICS_CHANNELS.md, docs/OPEN_ITEMS_SOLUTIONS.md,
docs/PACKAGING.md, README.md, docs/HEPMC3_CONVENTION.md (A10's append and
the one line), docs/benchmarking/08_channel_priority.md (row K's clause
count), docs/open_items/run_2026-09-03/STATUS.md (rows 3, 11, 17 notes; row
28 appended), docs/open_items/run_2026-10-27/phase_I6_cluster_wave.md (§9
item 11, one line), this record, STATUS.md (P4 row); the report
(docs/benchmarking/REPORT.json / REPORT.md, figures as §4 says).

## 4. Gates (final tree), tallies and wall times

| gate | result | wall |
|---|---|---|
| `cmake --build build -j2` | 0 warnings, 0 errors (after A17: 1 m 38 s; after A10–A12, one compile error of my own fixed, then 1 m 16 s; final no-op 0 s) | — |
| flip-only `./build/lipolgen_tests` (after A17) | 432 / 432 test cases, 17 244 559 assertions | 7 m 10 s |
| flip-only pytest (after A17) | 1 failed / 1516 passed / 155 skipped — the failure is the A11 assertion I had added to test_module.py before the A11 C++ was built (the Python tests run from source), not a regression; passes on the final tree | 17 m 36 s |
| `./build/lipolgen_tests` | 434 / 434 test cases (+2: A11, A10), 17 244 583 assertions, 0 failed | 6 m 54 s |
| pytest, full | 1519 passed, 155 skipped, 0 failed | 17 m 31 s |
| bench subset (`LIPOLGEN_BENCH=1 -k bench`) | 367 passed, 1 skipped, 1 failed — `test_the_committed_report_is_what_regeneration_gives`, run before the report was regenerated; after it, `test_bench_report.py` 26 passed | 6 m 21 s; 1 m 30 s |
| strict citation gate | 1243 / 95 / 7 / **0 broken** / 24 (SPIN32 19 / 115 / 0 / 0 / 10; PYTHIA_BRIDGE 0 / 0 / 8 / 0); one cited block (USAGE 1978–2169, §7b) re-recorded after re-reading rows 390 / 393 against the re-worded sentence | 2 s |
| `--audit-ranges` | 209 range citations, 0 REFUSED | 2 s |
| `--records-only` | 1432 citations in 86 records: 922 live, 61 into another record, 121 annotated, 129 read by hand, 199 skipped, **0 broken** | 2 s |
| `check_spdx_headers.py` | 136 files, all carry the line | <1 s |
| `pytest python/tests/test_doc_link_gate.py` | 112 passed | 24 s |
| `make_report.py --all --plots --out-dir docs/benchmarking --compare REPORT_5cd19d2.json` | 15 pass / 8 fail / 0 blocked; 1 change: `t6_hepmc3_convention` 23/23 → 24/24 clauses (C24 new, informational sub-row), status PASS → PASS; exit 1 by the "any change exits 1" rule; 19 figures differed by the footer stamp alone (pixel bbox inside the bottom 2 %), restored | 1 m 39 s |
| `make_report.py --check` | OK | 1 m 26 s |
| A9 (i) `build_nolhapdf` configure | exit 1 with the message (§2) | 0.2 s |
| A9 (ii) `build_off` configure / build / `lipolgen_tests` | exit 0 / exit 0, 0 warnings / 417 / 417 test cases, 17 243 193 assertions | 1 s / 2 m 17 s / 5 m 20 s |

No other row of the report moved.

## 5. Open points

1. **A17's reach.** The rule is applied to the rows the brief named and to
   the rc rows with the same cause; these rows still read `refused` AT THEIR
   DEFAULT on a run that does not read them: `kernel`, `b1_model`,
   `b1_band_scale`, `b1_alpha_d_dwave_weight`, `b1_unpol`, `r_source`,
   `cluster_vmc_mc_sigma`, `fsi` off the tags, and the rc sub-pieces closed by
   another knob (`rc_qe_*` at `rc_with_qe_tail = false`, `rc_sp_tensor_scale`
   on `t-peak` and — the opposite cause — on `polrad-full`, `rc_m_lepton` off
   `polrad-full`). Carrying the rule to them would leave the table no
   `refused` row at all (a refused non-default never builds), which is a
   decision about the status itself; CONVENTIONS names them. For the
   supervisor.
2. **A10's C++ examples without a Pipeline.** generate_inclusive,
   write_hepmc_example and hadronize_example write no `lipolgen.*` attribute;
   the C++ examples that do write them lack the two `*_sf_extrapolation` keys
   (the LHAPDF classification lives in the Python module's `sf_grid_extras`).
3. **A8.** t5_delta_moment_bound's exit 0 for a present-but-wrong set
   (`EnvironmentBlocked`) is outside A8's "missing set" and was left; whether
   the verdict configuration on a non-vendored set version should exit 2 is
   open.
4. The A12 band on the ceiling stays filed with the next FSI design
   (registry row 28).
