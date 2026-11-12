<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-27 — D5 design: `PolradFull` becomes the shipped `rc_tail` default

**Design only.** Written by an Opus worker on 2026-10-27 against `ca9b23a`
(clean tree, branch `claude/run-2026-10-26`); this file is the only repository
file written. Scratch measurements and scripts are in the session scratchpad
(`$SCR` = `/tmp/claude-1000/-home-cpeng-Projects-polli-PoliGen/7790c59b-ea07-4677-a36f-0183e2410ef4/scratchpad/d5/`),
outside the repository. No physics constant, tolerance, registry label or
reference number was moved.

**Authority.** The author, 2026-10-27, quoted in `DECISIONS.md` beside this
file: *"2. PolradFull is wanted."* — on `../run_2026-10-26/DECISIONS_PRICED.md`
§2, the bullet *"The `rc_tail` default model, `TPeak` (shipped) vs
`PolradFull`"* (its lines 181–188). `TPeak` stays available as an option.

## 0. The design in one paragraph

Flip one line, `include/lipolgen/rc.hpp:1503`
`RcTailModel tail_model = RcTailModel::TPeak;` → `RcTailModel::PolradFull`,
and the CLI's second definition of the same default,
`python/lipolgen/cli.py:817` `rc_tail_model="t-peak"` (derive it from the C++
default instead of retyping it). Everything else is labels, pins and
documents. **Nothing moves outside `--rc tensor-band` on the inclusive
channel**: `--rc off` (the generator's default) keeps every numeric column
bit for bit — one `meta` string changes, the not-read `knob_provenance` row
`rc_tail_model` (value `t-peak` → `polrad-full`); no reference JSON carries an
rc field; the polligen port gate reads no RC; the tail table is never built on
a tagged or coherent channel. **What moves** on a `--rc tensor-band` ⁶Li run
(measured, §2): the table build 0.23 s → 14.7 s, the mean dilution
⟨rc_tail⟩ − 1 ×1.36, the run-level A_zz with the tail −0.05 %, tail nodes
clipped 0 → 0.386 % (events: 0). Tests that pin **t-peak physics** pass
`TPeak` / `--rc-tail-model t-peak` explicitly; tests that pin **the default's
identity** are re-taken with a dated D5 comment; tests whose subject does not
depend on the tail model but which build a full-grid tail only incidentally
pass `t-peak` for runtime (each polrad-full build costs ~15 s). **Every edit to
a file the citation gates or the records cite by line is line-count-neutral**
(§7) — with that, the strict gate, the records gate and R5 stay green with no
`--fix` and no `--record-ranges`. The hashed design block of
`t4_djangoh_rad_noRad.py` is not touched; its rules name the models, not the
default, so row 7's verdict is unchanged and only its labels move.

## 1. What was read, and what it says

| source | what it says (quoted) |
|---|---|
| `include/lipolgen/rc.hpp:1503` | `RcTailModel tail_model = RcTailModel::TPeak;` — the default. |
| `rc.hpp:1229-1230, 1250, 1342` | `enum class RcTailModel : int {` `TPeak = 0, ///< DEFAULT: POLRAD Eqs. (37)-(39), (43)`; `PolradFull = 1, ///< OPT-IN: POLRAD Eq. (18) + Appendix B + Eq. (A.4)`; `TPeakPlusLL = 2, ///< OPT-IN`. The integer values are what the bindings and the npz carry; **they do not change**. |
| `rc.hpp:1226` | "`TPeak` remains the DEFAULT and is bit for bit every published number." |
| `rc.hpp:1331-1341` (PolradFull points (5), (6)) | "(5) IT COSTS, AND IT REFUSES RATHER THAN DEGRADES. ~10 s of tail-table build against ~0.15 s. It needs `n_eta` >= 64 … and a `long double` wider than `double` …; both are REFUSED by `RcModel`'s constructor, never silently degraded." "(6) NOT THE DEFAULT, for the same reason `TPeakPlusLL` is not: `TPeak` keeps every published number and every reference JSON bit for bit." |
| `rc.hpp:1269-1279` (point (3)) and `include/lipolgen/rc.hpp:1023-1029` (first bullet of "FIVE THINGS THIS IS NOT") | "ITS NORMALISATION IS STILL NOT CHECKED AGAINST MO-TSAI, or against any external exact tail"; "NOT validated against Mo-Tsai, or against any external exact tail … Nothing more." Both partly stale since row 7 (`../run_2026-10-26/SUMMARY.md` item 7). |
| `rc.hpp:60-61, 260-261` | "AND THE DEFAULT IS A LOWER BOUND, WITH TWO OTHER EDGES NOW SHIPPED. `rc_tail` at the shipped `RcTailModel::TPeak` is the t-PEAK ALONE." / "`rc_tail` at the DEFAULT `RcTailModel::TPeak` is the t-PEAK ONLY and is therefore a LOWER BOUND on the dilution." |
| `src/core/rc.cpp:1381-1424, 1813-1945, 1608-1619` | the constructor refuses `PolradFull` at `n_eta < 64` and where `std::numeric_limits<long double>::digits <= …<double>::digits` ("run --rc-tail-model t-peak or t-peak+ll"); `build_tail_tables()` fills `sigma_u_`, `sigma_t_`, `sigma_qe_`, `born_pb_` in two nested serial loops over `node_x_ × node_y_`; it is called only `if (tail_applies_)`. |
| `src/core/pipeline.cpp:1823-1831` | `if (cfg_.rc != PipelineRc::Off) { … rc_ = std::make_shared<RcModel>(…) }` — one model per `Pipeline`, built in its constructor. |
| `src/core/pipeline.cpp:731-777` | the tail sub-knobs, `{"--rc-tail-model", rc_options.tail_model != od.tail_model}` among them, are refused on a channel where `rc_tail_applies` is false whenever they differ from `RcOptions()`. |
| `../run_2026-10-26/phase_B3_djangoh.md` §0, §5, §6, §10.1 | Row 7 PASS. Through the shipped table R = 1.0199 / 1.0261 / 1.0219 / 1.0098 (`polrad-full`) and 0.6118 / 0.8968 / 1.0112 / 1.0067 (`t-peak`); straight from the quadrature "LiPolGen's exact tail is **1.0004 / 1.0010 / 0.9995 / 0.9866** of DJANGOH's" and t-peak 0.59926 / 0.87775 / 0.98948 / 0.98465 (§5 table). §6: "The shipped default `TPeak` under-weights the proton tail at low Q² … 39 % low at Q² = 1–10 and 10 % low at 10–100." §10.1: "Cost of switching: the tail table build goes from ~2 s to 10–21 s, and every published number moves (`rc.hpp` `RcTailModel::PolradFull` point (6))." Scope (§7): one beam (e 18 × p 275), one shared dipole, the proton (quasi-elastic) slot only; no coherent ⁶Li, no Pauli factor, no tensor term; bin 4's −1.34 % "cause unknown". |
| `../run_2026-10-26/DECISIONS_PRICED.md` §2 (lines 181–194) | the TPeak-vs-PolradFull bullet (numbers as above) "Row 17 untouched. None taken."; the next bullet: the (x, y) table's interpolation "raises the window-integrated proton tail by +2.0 … +3.7 % (+1.9 … +2.5 % net of the `tail_max` clip), both models" — a separate open candidate. |
| `validation/benchmarks/t4_djangoh_rad_noRad.py:41-187` (the hashed block) | line 148 "LiPolGen `polrad-full` (RcTailModel.PolradFull, opt-in)"; line 171 "PASS(polrad-full) iff \|R_b - 1\| <= 0.03 + 2 delta_b in EACH bin."; line 172 "LiPolGen `t-peak` (RcTailModel.TPeak, the SHIPPED default)"; line 178 "PASS(t-peak) iff R_b <= 1 + 0.03 + 2 delta_b in EACH bin"; lines 180–181 "its shortfall 1 - R_b is RECORDED as the price of the shipped default on this window, never asserted."; line 182 "The row PASSES iff both hold in all four bins". |
| `../run_2026-09-06/phase_B_numbers.md` §B2.5 (lines 900–913) | the price of the switch on ΔA_zz at Q² = 5: `polrad-full` − `t-peak` = 5.908189e−08 / 2.974583e−07 / 2.377717e−08 at x = 0.01 / 0.10 / 0.30 = **0.0435 % / 0.3073 % / 3.254 %** of the band half-width; "the band remains ×768–×1086 wider than either". |
| `../run_2026-09-03/STATUS.md` (the registry, lines 37–94) | no row names the tail model (grep of rows 1–27 for `tpeak\|t-peak\|tail_model\|rc_tail\|RcTailModel\|polrad-full\|PolradFull\|tail model\|rc-tail`: no hit). Row 17 is the RC **band's** low-x anchor; its text carries "the band still leads the whole radiative tail at x = 0.01 by ×768 / ×681 / ×289", a number measured against the t-peak tail. |
| `../run_2026-09-27/DECISIONS.md` | the pattern for a decision entry: "**D2 — the Hulthén deuteron stays the default; physics numbers use AV18.** <evidence> … The default is kept because every reference pin, the polligen port gate and the knob-provenance tables are defined on it … Revisit at M2". |

The POLRAD equation numbers above are the tree's (`rc.hpp`); the POLRAD 2.0
paper (Akushevich *et al.*, CPC 104 (1997) 201) was not re-read for this design.

## 2. What was measured on this machine (2026-10-27)

4 cores shared with another worker; load average 1.2–5.8 during the
measurements (each line says which). `source env.sh`, `build/` of 2026-10-26
13:04 (no C++ source changed since 2026-09-27); scripts in `$SCR`.

### 2.1 The tail-table build time — `Pipeline` construction, `--rc tensor-band` minus `--rc off`

| target, configuration | nodes | `t-peak` | `polrad-full` | script |
|---|---|---|---|---|
| ⁶Li config 1, CLI default grid (load 1.2–1.4) | 101 × 77 = 7777 | **0.226 s** | **14.71 s** | `m1_build_time.py` |
| same, repeated (load 3.6–5.1) | 7777 | 0.30 / 0.31 s | 15.29 / 14.78 s | `m3_azz_dilution.py`, `m5_ff_cost.py` |
| same, the CLI's own "setup" line (5000 events) | 7777 | 0.31 s | 15.34 s | `run5k_*.log` |
| ⁶Li, `--rc-c0-shape vmc-ft` (load ~5) | 7777 | 2.35 s | **21.01 s** | `m5_ff_cost.py` |
| ⁶Li, a Python `Spin1ElasticFF` subclass (load ~5) | 7777 | 2.36 s | **24.44 s** | `m5_ff_cost.py` |
| ⁶Li, grid 40 × 24 (C++ `rc_sampler()`'s shape) | 41 × 29 = 1189 | 0.039 s | 2.34 s | `m1_build_time.py grids` |
| ⁶Li, grid 40 × 28 (`test_rc_pipeline.cpp` `base_cfg`) | 41 × 33 = 1353 | 0.038 s | 2.69 s | same |
| ⁶Li, grid 20 × 12 | 21 × 29 = 609 | 0.019 s | 1.19 s | same |
| **proton**, row 7's harness (`t4_djangoh_rad_noRad.build_pipeline`, config 2, zero coherent FF as a Python subclass) | 7777 | 2.46 s | **23.08 s** | `m2_proton.py` |
| proton, CLI default grid, config 1 / config 2 (same zero FF) | 7777 | 2.67 / 2.76 s | 25.10 / 23.34 s | `m2_proton.py` |

≈ 1.9 ms per node for `polrad-full` against ≈ 0.03 ms for `t-peak` with the
C++ form factor; a Python form factor adds its trampoline cost per τ node.
The proton's "~2 s" of `DECISIONS_PRICED.md` is the t-peak with the harness's
Python zero form factor; its "10–21 s" is the polrad-full build (23 s here
under load). **Memory:** maxrss 41.7 MB either way (`/usr/bin/time`).

### 2.2 Caching — none

Read in `src/core/rc.cpp:1608-1619` (`if (tail_applies_) { … build_tail_tables(); }`)
and `src/core/rc.cpp:1813-1945` (`RcModel::build_tail_tables`): the tables are member
vectors of `RcModel`, filled once in its constructor by a **serial** loop and
held in memory for that object's life; `Pipeline` builds one `RcModel` per
construction (`pipeline.cpp:1823-1831`). There is no disk cache, no
cross-`Pipeline` cache and no thread pool; the only `static` caches in
`rc.cpp` are the per-thread `m_values` (`src/core/rc.cpp:1342-1345`) and the VMC C0 table
(`src/core/rc.cpp:947-951`). Every `lipolgen-run --rc tensor-band`, `lg.run(rc=…)` and
`lg.Pipeline(make_config(rc=…), plan)` on an inclusive ⁶Li run pays the build.
Nothing is built on a tagged or coherent channel (`if (tail_applies_)`), nor
on ⁷Li, whose inclusive tail throws in `HoSpin1FF::for_ion` on both models
("ion '7Li' has no measured-moment block in rc.hpp", measured).

### 2.3 Who builds the table — the census

**pytest**, counted by a scratch plugin (`$SCR/d5count_plugin.py`, wrapping
`_lipolgen.Pipeline`; runs `pytest -p no:cacheprovider`, bytecode off) over
`test_rc.py` + `test_rc_low_x_anchor.py` (37 passed, 72.6 s) and
`test_knob_provenance.py` + `test_cli.py::test_the_seven_untyped_flags_reach_the_run`
+ `test_bench_t4_t6.py` (530 passed / 152 skipped, 182 s). A "default-dependent
build" is a tail table built with `tail_model == TPeak` that nobody asked for
by name — it becomes a `polrad-full` build under D5:

| file | default-dependent builds | of which | explicit `polrad-full` builds today |
|---|---|---|---|
| `python/tests/test_rc.py` | 20 | rc_run fixture 1, a_transfer 4, meta_distinguishes 3, qe_tensor_scale_prices 3, c0_shape 2 (one `vmc-ft`), optimistic_delta 2, qe_refused 1, qe_sp_column 1, tail_model_knob 1, polrad_full 1, python_subclass 1 (Python FF) | 4 (57.2 s) |
| `python/tests/test_rc_low_x_anchor.py` | 10 | 3 + 3 + 3 + 1 | 0 |
| `python/tests/test_knob_provenance.py` | 34 | the RC cells of the inclusive-6Li specs: tensor-thirds 19, transverse-tensor 5, tensor-flip 5, helicity-flip pe=0 5 (one `vmc-ft`); each spec's first RC cell also builds that spec's cached `RC` base | 3 (41.2 s) |
| `python/tests/test_cli.py` | 2 | `test_the_seven_untyped_flags_reach_the_run` (one `vmc-ft`) | 0 |
| `python/tests/test_bench_t4_t6.py` | 1 | T6's `--rc tensor-band` file | 0 |
| `test_bench_djangoh_eltail.py`, `test_bench_chain_rc.py` | 0 | the t4 harness names both models explicitly | (the harness's own proton builds, unchanged) |
| **total** | **67** | ≈ +1000 s (≈ +17 min) if every one became `polrad-full` (63 × 14.8 s + 3 × 21.0 s + 1 × 24.4 s ≈ 1020 s, less the 11.7 + 12.6 s they cost today) — the phase-E full suite took 467 s | |

**C++** (static count from `tests/test_rc.cpp` / `tests/test_rc_pipeline.cpp`,
not run): ≈ 36 inclusive `RcModel` constructions with a default or local
`RcOptions` on the 40 × 24 grid (≈ 2.3 s each under `polrad-full`) and 2 on the
40 × 28 grid (T6, T18b, ≈ 2.7 s): ≈ +85 s if all flipped.

**CLI paths and harnesses that build it:** `lipolgen-run --rc tensor-band` on
`--channel inclusive` (⁶Li; any `--plan` with λ_e·P_e = 0); `lg.run(…, rc="tensor-band")`;
`make_config(rc=…)` + `Pipeline`; `validation/benchmarks/t6_hepmc3_convention.py`
(one file of eight, `RUNS[-1]`); `t4_djangoh_rad_noRad.py` (both models by
name — unaffected); `make_report.py --all` (runs both). No other harness and no
`validation/*.py` script builds an `RcModel` (grep for `RcModel|tail_ratio|rc_model|rc=|--rc|tensor-band|RcOptions|rc_tail`:
only `t3_li6_ff_wiringa_schiavilla.py:212`, a docstring about the t-peak's q reach).
`examples/*.cpp` and `.github/workflows/ci.yml` use no RC.

### 2.4 The size of the change on the RC dilution and on A_zz

**Event level** — `python3 -m lipolgen --isotope 6Li --channel inclusive --events 5000 --seed 20261027 --rc tensor-band --rc-tail-model {t-peak|polrad-full} --npz …`
(CLI defaults otherwise: config 1, `tensor-thirds`, P_z 0.7, P_zz 0.6; toy SF
backends). The two files carry **bit-identical events** (category, x, Q² —
`RcModel::fill` takes no `Rng&`):

| quantity | `t-peak` | `polrad-full` |
|---|---|---|
| wall / CLI "setup" | 0.49 s / 0.31 s | 15.56 s / 15.34 s |
| mean `rc_tail` (max) | 1.020341 (3.2380) | **1.028010** (4.0053) — ⟨rc_tail⟩ − 1 ×1.377 |
| per category azz+ / azz− / azz0 | 1.020453 / 1.023977 / 1.016594 | 1.028190 / 1.032881 / 1.022962 |
| ⟨rc_tail − 1⟩ by x: < 0.01 (n 2631) / 0.01–0.03 (1089) / 0.03–0.1 (997) / 0.1–0.3 (261) / 0.3–1 (22) | 3.83e−02 / 8.00e−04 / 1.75e−04 / 5.49e−05 / 7.10e−06 | 5.08e−02 / 3.18e−03 / 2.72e−03 / 5.61e−04 / 1.07e−04 (×1.33 / ×3.98 / ×15.6 / ×10.2 / ×15.0) |
| clipped tail nodes, global (y < 0.5 / 0.5–0.9 / > 0.9) | 0 % (0 / 0 / 0) | **0.386 %** (0.456 / 0.124 / 0 %) |
| clipped events, tail / band | 0 / 0 | 0 / 0 |
| banner knob block | read 33, refused axes 9 (`rc_m_lepton` among them) | read 34, refused axes 8 (`rc_m_lepton` READ) |

**Cell level, N-independent** (`$SCR/m3_azz_dilution.py`): per category k,
B_k = Σ_c Σ_m p_km σ_c (1 + w_avg,m(c)) and R_k = B_k + Σ_c Σ_m p_km σ_c
`tail_ratio(c, q_n(m))`, q_n(m) = 3m² − 2 at θ_S = 0; A_zz =
`estimators.azz_thirds(·, ·, ·, 0.6)`, the run-level estimator of
`validation/price_defaults.py:_thirds`. The recipe reproduces the documented
cell-σ-weighted means of `rc.hpp:1322-1325` exactly (1.02217080 t-peak,
1.03020526 polrad-full at `tail_max` = 10):

| quantity | `t-peak` | `polrad-full` | change |
|---|---|---|---|
| ⟨rc_tail⟩ − 1, cell-σ-weighted | 2.217080e−02 | 3.020526e−02 | ×1.3624 |
| run-level A_zz, Born | −5.193231e−04 | −5.193231e−04 | — |
| run-level A_zz with the tail | −5.173259e−04 (dilution 0.3846 %) | −5.170596e−04 (dilution 0.4359 %) | **×0.999485 (−0.0515 %)** |

The run-level A_zz with only 5000 events is dominated by sampling noise
(the category contrast is 4.7 × 10⁻⁴ of σ), so no event-level A_zz is quoted.

### 2.5 The per-point numbers the documents publish (`$SCR/m7_published_points.py`, `m8_qe_share.py`)

At Q² = 5 GeV², x = 0.01 / 0.03 / 0.10 / 0.30 (σ pieces straight from the
quadrature, `tail_sigma_at`; w_tail − 1 from the table, q_n = 0):

| quantity | `t-peak` (= the published numbers) | `polrad-full` |
|---|---|---|
| (1/6)σ^el_T/σ^el_U, `ho` | −5.0943e−04 / −7.3583e−04 / **+1.5595e−04** / +1.3287e−02 | −5.1540e−04 / −7.4709e−04 / **+4.4562e−05** / +1.0947e−02 |
| (1/6)σ^el_T/σ^el_U, `vmc-ft` | −5.4705e−04 / −8.0426e−04 / **−4.0702e−05** / +1.5556e−02 | −5.5461e−04 / −8.1840e−04 / **−1.6369e−04** / +1.1293e−02 |
| w_tail − 1, `ho` | 5.18222e−04 / 3.54397e−05 / 1.31963e−06 / 6.86176e−08 | 6.62870e−04 / 1.75050e−04 / 1.95963e−04 / 4.88356e−04 (×1.28 / ×4.94 / ×148 / ×7117) |
| σ^q_U/σ^el_U at x = 0.01 / 0.10 / 0.30 | 0.28237 / 2.7475 / 1000.2 | 0.63146 / 533 / 3.6846e+06 |
| quasi-elastic share of the tail, x = 0.01 / 0.10 / 0.30 | **22.0 % / 73.3 % / 99.90 %** | **38.7 % / 99.81 % / 100.00 %** |

The t-peak rows reproduce the C0 box rows of `docs/USAGE.md:2258-2259` (`ho` −5.0943e−04 …,
`vmc-ft` −5.4705e−04 …) and `docs/PHYSICS_CHANNELS.md:386-387`
digit for digit — **those published numbers are t-peak numbers**. The C0
band's sign flip at x = 0.10 survives on `polrad-full` (+4.46e−05 / −1.64e−04).
The published "w_tail − 1 = 5.18398e−04 (x = 0.01), 1.32298e−06 (0.10),
6.8624e−08 (0.30)" (`PHYSICS_CHANNELS.md:386`) is a t-peak number; this table's
table-route reading differs from it by −0.03 % (its recipe was not
reproduced exactly; recorded, not resolved).

### 2.6 T6 runtime (`$SCR/m4_t6_runtime.py`, the harness in process with its rc file's argv extended by `--rc-tail-model …`)

`t-peak` **6.06 s**, `polrad-full` **21.07 s** (load 4.5–5.2), 23 / 23 clauses
PASS on both. `python/tests/test_bench_t4_t6.py:test_t6_measured_row_is_pinned`
asserts `rep["runtime_s"] < 30.0`.

### 2.7 A/B column hashes for the implementation's verification (`$SCR/npz_hash.py`)

sha256 over every array column of the npz (sorted keys, dtype + bytes, `meta`
excluded), the commands of §2.4 on `ca9b23a`:

| run | column sha256 |
|---|---|
| `--rc tensor-band --rc-tail-model t-peak` (today's default) | `2a312c26af9f18d140f6f40cccc720b96ac76d3f4942e99769b8f53476ef480a` |
| `--rc tensor-band --rc-tail-model polrad-full` | `00c012b459e25cfe0d1259abc8c39c4d8d1c606da90dfa7efcde114fe793892e` |
| `--rc off` (same seed, no rc flag) | `2c48b4e516fea45503b09627451b72414b9e921309124dcb962c55d2b7759660` |

`meta` keys that differ between the first two: `knob_provenance` (rows
`rc_m_lepton`, `rc_n_eta`, `rc_sp_tensor_scale`, `rc_tail_model`),
`rc_clipped_cell_fraction`, `rc_clipped_fraction_by_y`, `rc_tail_model`. At
`--rc off` the row reads `{'value': 't-peak', 'status': 'not-read', 'reason':
'not read at rc = off: no RcModel is built', 'at_default': True, …}` (measured
with `lg.run(events=50, seed=3)`). **These hashes are valid only on a tree
where D4 / D6 / D7 have not changed the sampled events**: re-take all three on
the tree I5 starts from (§8 step 0).

### 2.8 One finding: `tests/test_rc.cpp` T8(0)'s x = 1 clause is t-peak-only

T8(0) (`test_rc.cpp:775-821`) asserts, at the table's last x edge, "the QRT is
in NUCLEON invariants, so at the grid's last edge x = 1 it vanishes
identically" (`CHECK(sq[i * ny + j] < 1e-30)`). Emulated from Python on the
same 40 × 24 generator-scenario grid (`$SCR/m6_t8zero.py`): the sign gate
holds on both models (min σ^el_U 9.07e−31 / 9.08e−31 > 0, σ^q_U > 0 off the
edge, `tail_ratio_at > 0` at every T8(0) point and q_n ∈ {0, 1, −2}); **σ^q_U at
the x = 1 edge is ≤ 2.3e−33 on `t-peak` and 3.124e+08 on `polrad-full`.** Not a
defect (Eq. (18)'s s-/p-peaks reach the nucleon elastic point; the Born is 0
there and `tail_ratio_at` clips at `tail_max`), but T8(0) must stay on `TPeak`.

## 3. Inventory — every site that encodes TPeak as the default or as "shipped"

Grepped over the whole tree for `TPeak|t-peak|t_peak|PolradFull|polrad-full|tail_model|tail-model|RcTailModel|rc_tail|rc-tail|shipped|default`.
**LN** = the edit must be line-count-neutral (the file is cited by line, §7).
Dated run records under `docs/open_items/run_*/` and the dated surveys are
history and are not rewritten (§3.6).

### 3.1 Library code

| file | line(s) | what | action |
|---|---|---|---|
| `include/lipolgen/rc.hpp` | 1503 | `RcTailModel tail_model = RcTailModel::TPeak;` | **→ `RcTailModel::PolradFull`** (LN). The switch. |
| `rc.hpp` | 1224-1227 | "All THREE are implemented since 2026-09-06; `TPeak` remains the DEFAULT and is bit for bit every published number." | LN rewrite: "`PolradFull` is the DEFAULT since 2026-10-27 (D5); `TPeak` was until then, and every RC-tail number published before that date without a model label is a `TPeak` number." |
| `rc.hpp` | 1230-1249 | `TPeak = 0, ///< DEFAULT: …` "IT IS ONE PEAK OF THE ELASTIC TAIL AND ITS ABSOLUTE NORMALISATION IS NOT VALIDATED AGAINST ANY EXTERNAL EXACT TAIL." | LN: "OPTION since D5 (the DEFAULT until 2026-10-27)"; replace the "NOT VALIDATED" sentence by row 7's measurement: 0.599 / 0.878 / 0.989 / 0.985 of DJANGOH's exact proton tail at Q² = 1–10 / 10–100 / 100–10³ / 10³–10⁴ GeV² (e 18 × p 275), i.e. low, as a lower bound should be. Enum value 0 unchanged. |
| `rc.hpp` | 1250 | `PolradFull = 1, ///< OPT-IN:` | LN: "DEFAULT since 2026-10-27 (D5; opt-in 2026-09-06 → 10-26):" |
| `rc.hpp` | 1269-1279 | point (3) "ITS NORMALISATION IS STILL NOT CHECKED AGAINST MO-TSAI, or against any external exact tail" | LN: "Mo-Tsai is not in this tree. Since 2026-10-26 the PROTON (quasi-elastic) slot is checked against DJANGOH 4.6.21's exact O(alpha) tail (BENCHMARK_PLAN.md row 7): 1.0004 / 1.0010 / 0.9995 / 0.9866 straight from the quadrature, +1.9 … +2.5 % through the (x, y) table; the coherent 6Li slot and the tensor term are not checked externally." Keep the POLRAD-internal sentence. |
| `rc.hpp` | 1331-1338 | point (5) "~10 s of tail-table build against ~0.15 s" | LN, optional: "(14.7 s against 0.23 s on 6Li config 1, 101 x 77 nodes, a loaded 4-core machine, 2026-10-27)". |
| `rc.hpp` | 1339-1341 | point (6) "NOT THE DEFAULT … `TPeak` keeps every published number and every reference JSON bit for bit." | LN: "(6) THE DEFAULT since 2026-10-27 (D5, docs/open_items/run_2026-10-27/DECISIONS.md). No reference JSON carries an rc field; RC-tail numbers published before D5 reproduce with tail_model = TPeak." |
| `rc.hpp` | 1416-1418 | `TPeakPlusLL` (3) "NOT THE DEFAULT. `TPeak` stays the shipped default so that every published number…" | LN: "NOT THE DEFAULT (`PolradFull` is, since D5; `TPeak` was until then)." |
| `rc.hpp` | 60-61, 68, 73-79 | header: "AND THE DEFAULT IS A LOWER BOUND, WITH TWO OTHER EDGES NOW SHIPPED. `rc_tail` at the shipped `RcTailModel::TPeak` is the t-PEAK ALONE." and the PolradFull bullet "ABSOLUTE NORMALISATION IS STILL NOT CHECKED AGAINST MO-TSAI OR ANY EXTERNAL EXACT TAIL (no such number is in this tree)" | LN rewrite: the default is Eq. (18) since D5; the t-peak pair stays as options; the PolradFull bullet gets row 7 as in point (3). |
| `rc.hpp` | 260-262 | HONEST FLAGS: "`rc_tail` at the DEFAULT `RcTailModel::TPeak` is the t-PEAK ONLY and is therefore a LOWER BOUND" | LN: "`rc_tail` under `RcTailModel::TPeak` (the default until D5) is the t-PEAK ONLY …". Keep line 268 as it is (a record cites it). |
| `rc.hpp` | 1021-1029 | "FIVE THINGS THIS IS NOT", first bullet "NOT validated against Mo-Tsai, or against any external exact tail … Nothing more." | LN: as point (3). (PHYSICS_CHANNELS.md does not cite these comment lines; it cites the declarations below them, from `PolradIm` at line 1070 onward — hence LN.) |
| `rc.hpp` | 47-51, 1544-1546, 1956 | "the quasi-elastic term is 22 % / 73 % / 99.9 % of the tail at x = 0.01 / 0.10 / 0.30" | LN: name the model: "(t-peak; 39 % / 99.8 % / 100.0 % on PolradFull, the default since D5; Q^2 = 5)". |
| `rc.hpp` | 1624-1627 | `sp_tensor_scale`: "unreachable under the shipped `TPeak` anyway" | LN: "unreachable under `TPeak` and refused on the default `PolradFull` (the term ran)". |
| `rc.hpp` | 1441-1443, 1503 | `rc_tail_model_name` — unchanged; the declaration line `:1503` keeps declaring `tail_model` (PHYSICS cites `:1503` `tail_model`, `:1229` `RcTailModel`). | — |
| `src/core/rc.cpp` | 1558 | the tagged-channel `exclusion_reason` string (written to `meta["rc_exclusion_reason"]`): "inclusive quasi-elastic one (22 % / 73 % / 99.9 % of the inclusive …" | LN: name the model as above. Runtime text; `test_rc.py::test_tagged_tail_exclusion_names_the_quasi_elastic_omission` checks tokens only (`x_L`, `FACT`, `OMISSION`, `QUASI-elastic`, `2/3`, `SPECTATOR`). |
| `rc.cpp` | 1986-1989 | comment "not one floating-point operation of the shipped default moved" (about `TPeak`) | LN: "of the t-peak model (the default until D5)". |
| `rc.cpp` | 1381-1424 | the `PolradFull` refusals (long double; `n_eta < 64`) advising "run --rc-tail-model t-peak or t-peak+ll" | no change — that advice becomes the default path's fallback on a narrow-long-double platform (§10 R2). |
| `src/core/pipeline.cpp` | 679-684 | `validate()` message "the shipped TPeak tail has no lepton-mass dependence" | LN: "the TPeak tail (an option since D5)". Tests match `"m_lepton"` only. |
| `pipeline.cpp` | 3126-3131 | comment "under the shipped `TPeak` the `u_sp` table is identically zero" | LN reword. |
| `pipeline.cpp` | 743, 3122-3125, 3133-3166, 3193-3197, 3203-3229 | the tail-knob refusal list and the knob-provenance rows (`rc_tail_model`, `rc_sp_tensor_scale`, `rc_n_eta`, `rc_m_lepton`) | no code change. Behaviour at the new default (measured, §2.4): `rc_m_lepton` READ instead of refused, `rc_n_eta`'s description is the tanh-sinh one, `rc_sp_tensor_scale`'s refusal reason is the "opposite reason" branch; on tagged channels `--rc-tail-model t-peak` becomes the refused value and `polrad-full` the accepted one. |
| `include/lipolgen/pipeline.hpp`, `include/lipolgen/constants.hpp:66-75`, `include/lipolgen/event.hpp` | — | no TPeak-default statement (`constants.hpp`'s "refuses a non-default m_lepton ON THOSE TWO ONLY" stays true) | none |
| `python/bindings.cpp` | 847, 857 | comments: "under the shipped `t-peak` tail they do not exist at all"; "the QRT is 73-99.9 % of the tail" | LN reword / name the model. |
| `bindings.cpp` | 3425-3441 | `RcTailModel` enum docstring "TPeak is the DEFAULT and is bit for bit every published number … PolradFull … still is NOT checked against Mo-Tsai or any external exact tail" | LN rewrite (string literal lines). |
| `bindings.cpp` | 3533-3542 | `tail_model` docstring "RcTailModel.TPeak (default, bit for bit every published number)" | LN: "RcTailModel.PolradFull (default since D5) … RcTailModel.TPeak (the default until 2026-10-27)". |
| `bindings.cpp` | 3602 | `sp_tensor_scale` docstring "unreachable under the shipped t-peak tail anyway" | LN reword. |
| `bindings.cpp` | 871-872 | `meta["rc_tail_model"] = rc_tail_model_name(...)` | none (reads "polrad-full" automatically). |
| `python/lipolgen/cli.py` | 817 | `DEFAULTS[...]`: `rc_tail_model="t-peak",` | **→ `rc_tail_model=_l.rc_tail_model_name(_l.RcOptions().tail_model),`** (LN; one definition, the C++ one — CONVENTIONS.md:20 "No physics number is hard-coded in two places", and `test_rc.py:999-1004`'s own rule "assert against the MODULE CONSTANTS, not against re-typed literals"). A plain `"polrad-full"` is the fallback if the author prefers a literal. |
| `cli.py` | 290-300 | `--rc` help: "At the DEFAULT --rc-tail-model t-peak, rc_tail is the t-PEAK ONLY: a LOWER BOUND …" | LN rewrite: the default is polrad-full (D5); t-peak is the lower-bound option. |
| `cli.py` | 356-423 | `--rc-tail-model` help: "\"t-peak\" (default, bit for bit every published number …" and "ITS ABSOLUTE NORMALISATION IS STILL NOT CHECKED AGAINST MO-TSAI OR ANY EXTERNAL EXACT TAIL" | LN rewrite: polrad-full is the default (D5), checked on the proton slot against DJANGOH (row 7); t-peak (the default until 2026-10-27) and t-peak+ll are options. |
| `cli.py` | 440-441 | `--rc-qe-tensor-scale` help "the piece that is 73 %% of rc_tail at x = 0.1 and 99.9 %% at x = 0.3" | LN: the default model's 99.8 % / 100 % (t-peak's 73 / 99.9 %). |
| `cli.py` | 1658-1662 | banner "tensor-blind the piece that is 73 % of it at x = 0.1 and 99.9 % at x = 0.3" | LN: print the running model's shares (t-peak 73 / 99.9 %, polrad-full 99.8 / 100 %), or name the model. |
| `cli.py` | 1673-1684 | banner, polrad-full paragraph: "ITS ABSOLUTE NORMALISATION IS STILL NOT CHECKED AGAINST MO-TSAI OR ANY EXTERNAL EXACT TAIL: no such number is in this tree." | LN: row 7 (proton slot checked; coherent and tensor not). Printed on every default `--rc tensor-band` run from D5 on. |
| `cli.py` | 1686-1700 | banner, t-peak paragraph ("Read rc_tail as a LOWER BOUND …") | keep; it is now the option's banner. |
| `python/lipolgen/__init__.py` | 347-378 | `RC_TAIL_MODELS` comment: '"t-peak" is the DEFAULT and is bit for bit every published number …'; '"polrad-full" … It is NOT the default and NOT validated against any external exact tail' | LN rewrite (the use-site pins at lines 422, 634 and 694 sit below). |
| `python/lipolgen/export.py`, `configs.py` | — | rc weight names only | none |

### 3.2 C++ tests

| file | line(s) | what | action |
|---|---|---|---|
| `tests/test_rc.cpp` | 346-351 | `rc_model()`: "/// The default 6Li inclusive model at the run's defaults." `static const RcModel m(RcMode::TensorBand, RcOptions(), …)` — the base of T8(0), T8(a'), T5, T18, T1, T4, B3/B1/B6 and T8(d)/T19's `mt` | **pass `TPeak` explicitly** (LN: an inline `[] { RcOptions o; o.tail_model = RcTailModel::TPeak; return o; }()` in the argument list); comment → "the 6Li inclusive model on the t-PEAK tail (the default until D5, 2026-10-27)". T8(a') compares `tail_sigma_at` with `polrad_sigma_el_u`/A² at 1e-12; T8(0)'s x = 1 clause is t-peak-only (§2.8). |
| `test_rc.cpp` | 1314-1331 (T8'), 1357-1511 (B3), 1563 and 1569 (B1), 1824-1830 (Q9), 2444-2449 (`band_model`, T12/T12') | local `RcOptions o;` / `RcOptions()` on the inclusive tail | **pass `TPeak`** where the case reads tail numbers (T8', B3, Q9, T12/T12', B1 `bad` — whose comment "// tail_model = TPeak" is what it means to test); LN above line 1877 (append `o.tail_model = RcTailModel::TPeak;` to the declaration line). **The authority is the flipped build**: flip, build, run, and every case that fails must be shown to pin t-peak physics before it is pinned (§8 step 1). |
| `test_rc.cpp` | 1544, 2724 | comments "// the shipped TPeak default", "// the shipped default, TPeak" | LN reword ("the t-peak model, the default until D5"). |
| `test_rc.cpp` | 2154 (T17), 2562 / 2628 (T14), 1991-2009 (B6), 3478-3480 (TensorAll), 1859 / 1867 (y_max refusal, band-only) | model-independent subjects | leave on the default (≈ 8 polrad-full builds × 2.3 s ≈ +20 s); B6 and TensorAll may take `TPeak` if the C++ runtime matters. |
| `test_rc.cpp` | end of file (after 3515) | — | **new** `TEST_CASE("D5: PolradFull is the shipped rc_tail default; TPeak stays an option")`: `CHECK(RcOptions().tail_model == RcTailModel::PolradFull)`; the name is `"polrad-full"`; one default-built model's node equals an explicit `PolradFull` model's (or T19(i)'s pinned value at one point); `rc_model().options().tail_model == RcTailModel::TPeak`. Appended at the end so no cited line moves. |
| `tests/test_rc_pipeline.cpp` | 247 | `cfg.rc_options.tail_model = RcTailModel::TPeak;` (the reset at the end of the m_lepton SUBCASE, "back to the default") | → `RcOptions().tail_model`. |
| `test_rc_pipeline.cpp` | T6 (`:107-157`), T6c (`:183-260`), T18b (`:262-275`) | bit-for-bit RNG gate, threading, validate | keep on the default (+2 × 2.7 s); T6c's `n_eta = 2` still throws (`validate()`'s general `n_eta < 8` check, `pipeline.cpp:544`, comes first); note that from D5 any `n_eta` in 8–63 is refused at the default (`src/core/pipeline.cpp:625`, `rc_options.n_eta < 64`). |

### 3.3 Python tests

| file | line(s) | what | action |
|---|---|---|---|
| `python/tests/test_rc.py` | 31-35 | `rc_run` fixture (`_cfg("tensor-band")`) | **keep on the default** — it is the suite's default-RC run (now polrad-full, +14.8 s once). |
| `test_rc.py` | 263 | `assert meta["rc_tail_model"] == "t-peak"` | **re-take** → `"polrad-full"`, dated D5 comment. |
| `test_rc.py` | 593-657 | `test_tail_model_knob_is_a_seven_file_chain`: 612 `base["meta"]["rc_tail_model"] == "t-peak"`; 650-651 `cli.DEFAULTS["rc_tail_model"] == "t-peak"`, `cli.resolve(["--rc", "tensor-band"])["rc_tail_model"] == "t-peak"`; the t-peak+ll comparisons against `base` (`rc_tail >= base`, "the three t-peak tables are untouched" against `m0` = default) | **re-take** the three default asserts (→ `"polrad-full"`; assert `cli.DEFAULTS[...] == _l.rc_tail_model_name(_l.RcOptions().tail_model)`); **pass t-peak** for the comparison base (an explicit `rc_tail_model="t-peak"` run and `m0`), +0.3 s. Docstring "The DEFAULT is asserted unchanged" → dated D5. |
| `test_rc.py` | 660-760 | `test_polrad_full_is_the_exact_tail_and_does_not_bracket_the_t_peak_band`: `base` = rc_run ("It ADDS the s- and p-peaks, so the tail never falls below the t-peak", `mean() > base.mean()`), `mt = lg.Pipeline(_cfg("tensor-band"), plan)` | **pass t-peak** for `base` and `mt` (the mean comparison fails trivially otherwise: equal means); `p` may reuse `rc_run` (saves one build). The T19(i) value pins are explicit polrad-full and stay. Docstring "opt-in since 2026-09-06", "(ii) that the DEFAULT is untouched" → dated D5. |
| `test_rc.py` | 190-206 | `test_c0_shape_reaches_the_tail_and_the_npz_meta`: `out["vmc-ft"][0] / out["ho"][0] == approx(1.1043, rel=1e-3)` and the x = 0.10 sign | **pass t-peak** (Eq. (38) numbers; the polrad-full sign also flips, §2.5, but its ratio was not measured). |
| `test_rc.py` | 499-529 | `test_the_quasi_elastic_sp_column_is_bounded_by_NEITHER_scale`: `tpk = lg.Pipeline(_cfg("tensor-band"), plan)` | **pass t-peak** (it is the t-peak side of a t-peak vs t-peak+ll comparison). |
| `test_rc.py` | 268-308, 341-377, 380-395, 532-571, 784-803, 935-974 | `meta_distinguishes` (3 builds), `qe_tensor_scale_prices` (3), `qe_tensor_scale_is_refused` (1), `a_transfer_frac` (4), `optimistic_delta` (2), `a_python_subclass_form_factor` (1, Python FF: 24 s on polrad-full) | **pass t-peak for runtime**: their subjects (scope label, the QE stand-in's algebra, the band, the FF's Python lifetime) do not read the tail model; dated comment "t-peak for runtime (D5): ~15 s per polrad-full build, subject model-independent". |
| `test_rc.py` | 105, 475-496 | comment / docstring "the shipped `t-peak` computes no s-/p-peaks" (`test_sp_tensor_scale_is_refused_when_the_sp_peaks_did_not_run` still raises on the default — the opposite-reason branch, `match="sp_tensor_scale"` holds) | reword. |
| `python/tests/test_rc_low_x_anchor.py` | 77-79 | `_pipeline()` helper, 10 builds; A1–A5 are band numbers | **pass t-peak for runtime** (`rc_tail_model="t-peak"` in the helper), dated comment. |
| `python/tests/test_knob_provenance.py` | 235 `RC = dict(rc="tensor-band")`, 351-441 | the RC cells | split: **`RC` (the shipped default)** for `rc: off -> tensor-band`, `rc_n_eta: 128 -> 64` (now a polrad-full variant, legal at 64) and the two `rc_tail_model` cells; **`RC_TPK = dict(RC, rc_tail_model="t-peak")`** for the model-independent band and tail sub-knob cells (`rc_delta_low_x … rc_tail_max`, the NO_QE cells); re-base and relabel: `"t-peak -> t-peak+ll"` → `"polrad-full -> t-peak+ll"`; `"t-peak -> polrad-full"` → **`"polrad-full -> t-peak"`** (`var=dict(RC, rc_tail_model="t-peak")`; as written it would be a no-op variant whose `read` row fails the gate); `"0 -> 1 at t-peak (the shipped tail)"` → base `RC_TPK`, label "at t-peak (an option since D5)"; `rc_m_lepton "M_ELECTRON -> the muon mass"` → base `RC_TPK` (keeps the t-peak refusal covered; on the default it would duplicate the "at polrad-full" cell). Comments 404-411 reworded. |
| `python/tests/test_cli.py` | 279-283 | `test_the_seven_untyped_flags_reach_the_run`: `common = [... "--rc", "tensor-band", "--quiet"]` (2 builds, one vmc-ft: 36 s on polrad-full) | **pass `"--rc-tail-model", "t-peak"`** in `common` for runtime (typed, identical in both runs; the untyped flags are the subject). |
| `python/tests/test_bench_t4_t6.py` | 377-403 | `test_t6_measured_row_is_pinned`: `rep["runtime_s"] < 30.0` | no edit: T6 stays on the user's own default path (21 s measured); fallback in §10 R1. |
| `python/tests/test_bench_djangoh_eltail.py` | 13, 277-278 | "the shipped t-peak tail's 39 % / 10 % shortfall", "# the shipped t-peak's shortfall" | reword ("the t-peak option's (the shipped default until D5)"); optional new assert: `"polrad-full (shipped)" in row["generator"] and "t-peak (option)" in row["generator"]`. R values (rel 2e−7) unchanged — both models are named. |
| `python/tests/test_bench_chain_rc.py` | — | both models by name | none |
| `python/tests/test_bench_report.py` | 162, 275, 288-291 | synthetic rows keyed `"t-peak"` | none |
| `python/tests/conftest.py:89-103` (the `polligen` fixture) and its users `test_kernel_reference.py`, `test_pipeline.py`, `test_sampler.py`, `test_hfs.py`, `test_cli.py` | — | the port gate | **none: no RC anywhere in it** (grep); the sibling repository has no `rc_tail` consumer either (grep of `PolarizedLithiumSim` outside `refs/`). |

### 3.4 Benchmark harnesses

| file | line(s) | what | action |
|---|---|---|---|
| `validation/benchmarks/t4_djangoh_rad_noRad.py` | 41-187 | the hashed design block ("the SHIPPED default", "opt-in") | **untouched** (`DESIGN_BLOCK_SHA256`, `:367`). |
| `t4_djangoh_rad_noRad.py` | after the ADDENDA (`:189-212`) | — | **new dated addendum** (text in §5). |
| `t4_djangoh_rad_noRad.py` | 955-958 | generator `"polrad-full %s; t-peak (shipped) %s"` | → `"polrad-full (shipped) %s; t-peak (option) %s"`. |
| `t4_djangoh_rad_noRad.py` | 972-981 | reason "… the shipped t-peak passes its one-sided CEILING only …" | → "polrad-full (the shipped default since D5) … the t-peak option passes its one-sided CEILING only …" — keep the words `CEILING` (asserted by `test_bench_djangoh_eltail.py:267`) and `one-sided` (in the tolerance). |
| `t4_djangoh_rad_noRad.py` | 1088-1093 | plot: `("polrad-full", "alt2", …)`, `("t-peak", "default", …)`; labels "LiPolGen polrad-full (opt-in)", "LiPolGen t-peak (shipped default)" | swap the style slots (the palette's `"default"` slot is "the shipped default", `_plots.py:50`) and relabel "polrad-full (shipped default)", "t-peak (option; a lower bound)". |
| `validation/benchmarks/_plots.py` | 56 | `"alt2": … # … (CD-Bonn, polrad-full, ...)` | → "(CD-Bonn, t-peak, ...)" with the swap. |
| `t6_hepmc3_convention.py` | 191-196, 114-115 | the `--rc tensor-band` file; "4.3 s wall for the whole row" | no edit to `RUNS`; optional dated note on the runtime (21 s on polrad-full). |
| `t3_li6_ff_wiringa_schiavilla.py` | 212 | "t-peak integrand reaches (q ~ 88-278 fm^-1 …)" | none (a statement about the t-peak region, true on both models). |
| `t5_sum_rules.py`, all `t1_*`, `t2_*`, other `t3_*`, `t4_pythia_ep_closure.py`, `validation/price_defaults.py` | — | no RC | none (`tiers.py --files` maps `rc.*` to T2, T3, T4-RC and `cli.py` to T6 — re-run, expect no change). |

### 3.5 Living documents (all LN; §7)

| file | line(s) | what | action |
|---|---|---|---|
| `README.md` | 51 | `--rc-tail-model` row: default column "`t-peak`, bit for bit" | → "`polrad-full` since 2026-10-27 (D5); `t-peak` until then — `--rc-tail-model t-peak` reproduces every earlier number"; the measured column (explicit-model numbers) stays. |
| `README.md` | 54 | `--rc-c0-shape` row: σ^el_T/σ^el_U at x = 0.10, Q² = 5 "+9.357e−04 on `ho`, −2.442e−04 on `vmc-ft`", "a factor 67 resp. ≈ 255 below … +0.062" | label "(t-peak)" and add the default's +2.674e−04 / −9.821e−04 (×232 / ×63 below +0.062) — 6 × §2.5. |
| `README.md` | 184-185 | "Still open: the exact elastic radiative tail (Mo–Tsai, beyond POLRAD's transcribed t-peak)" | → the shipped tail is POLRAD's exact Eq. (18) since D5, its proton slot checked against DJANGOH (row 7); Mo–Tsai not in the tree; the coherent ⁶Li tail unchecked externally. |
| `README.md` | 39-41 | "each selector below is off unless asked for" | D5 leaves `--rc` off by default, so the sentence holds for `--rc`; D6/D7 change this sentence anyway — coordinate. |
| `docs/USAGE.md` | 2247 | `--rc-tail-model` row ("`t-peak` is the default and the **lower** edge") | rewrite the default statement; keep the measured explicit-model numbers. |
| `USAGE.md` | 2244, 2246, 2249, 2253-2269 | `--rc-fq-scale` (the sign-change claim at x = 0.01), `--rc-c0-shape` ("Cost: `RcModel` construction 0.11 s → 1.29 s"), `--rc-qe-tensor-scale` (ΔA_zz −1.162e−07 / +3.018e−10 / +1.796e−09, 0.086 / 0.0003 / 0.25 % of the band), the C0 box table | label each "(t-peak)"; add the default's numbers where measured (§2.5: the box rows, in the same cells, LN); the c0-shape cost → 14.7 s → 21.0 s on the default; the fq-scale sign and the qe-tensor-scale ΔA_zz are **not re-measured on polrad-full** — label only. |
| `USAGE.md` | 2212, 2341 | "22 % / 73 % / 99.9 % of the inclusive tail" | name the model, add 39 / 99.8 / 100 %. |
| `USAGE.md` | 2352-2400 | "`--rc-tail-model t-peak` (the **default**, bit for bit every published number)"; "What `polrad-full` is still NOT checked against: Mo–Tsai, or any external exact tail" | rewrite the default sentence; add row 7 to the "NOT checked" paragraph (proton slot checked). |
| `USAGE.md` | 2540 | "Do not port the `t-peak` default to fixed-target kinematics" | → "the `t-peak` option". |
| `USAGE.md` | 2548-2579 | banner description ("the t-peak caveat above"); "At the current defaults **nothing clips on the tail**" | → the banner prints the running model's paragraph; at the default the nodes clip 0.386 % (0.456 / 0.124 / 0 % by y band), events 0 / 200 000 (the existing polrad-full sentence at 2564-2575 becomes the default's). |
| `USAGE.md` | 2581-2590 | "Cost: **0.106 s** at setup …; **9.84 s** on `--rc-tail-model polrad-full`" | → "~10–15 s at setup on the default (14.7 s measured 2026-10-27, a loaded 4-core machine); 0.106 s on `--rc-tail-model t-peak`"; the ≈ 14 % event-rate cost is the table lookup and does not depend on the model (not re-measured). |
| `USAGE.md` | 1978-2169 (fingerprinted, cited by PHYSICS:393 as §7b) | the §7b opening and the code example | **do not edit** (an edit inside a fingerprinted range is an EDITED block, §7). |
| `docs/PHYSICS_CHANNELS.md` | 386 | row heading "POLRAD's **t-peak** by default, and since 2026-09-06 its **exact Eq. (18)** as a third opt-in model"; "`RcTailModel::TPeak` is the **default** and is bit for bit every published number"; "its absolute normalisation is still **not** checked against [MT69]"; "Composite Gauss–Legendre quadrature in ln η_A (`n_eta` = 128 nodes in panels of 8)"; the (1/6)σ^el_T/σ^el_U and "w_tail − 1" numbers; "`rc_tail` at the DEFAULT `TPeak` is the t-peak alone" | LN (one table line): rewrite the default statements; label the numbers "(t-peak)"; add the default's (§2.5); keep the citations column untouched (`include/lipolgen/rc.hpp:1229` `RcTailModel` and `include/lipolgen/rc.hpp:1503` `tail_model` keep declaring their names). |
| `PHYSICS_CHANNELS.md` | 387, 841 | σ^q_U/σ^el_U = 0.28237 / 2.7475 / 1000.2 and "22 % … 73 % … **99.9 %**"; "at x = 0.30 the tail is 99.9 % quasi-elastic" | label (t-peak) and add the default's (0.63146 / 533 / 3.68e+06; 39 / 99.8 / 100 %). |
| `PHYSICS_CHANNELS.md` | 391 | "Measured … at the **default** `tail_model = TPeak`: **0 %** clipped nodes …" | → the default's 0.386 % (events 0); keep the t-peak 0 % labelled. |
| `PHYSICS_CHANNELS.md` | 393 | cost row: "construction 0.106 s … −13.9 % end-to-end throughput"; "the t-peak lower-bound warning" | → construction 14.7 s on the default (0.106 s on t-peak); the banner prints the running model's paragraph. |
| `PHYSICS_CHANNELS.md` | 747-812, 813-840 | "What ships instead is a **BAND OF TWO TAIL MODELS** … `RcTailModel::TPeak` (the default, bit for bit every published number)"; "`RcTailModel::PolradFull` … opt-in as `--rc-tail-model polrad-full`" | LN rewrite within the same lines (append to an existing line where more words are needed): the default is Eq. (18) since D5; the t-peak pair is a price range of options. |
| `docs/OPEN_ITEMS_SOLUTIONS.md` | 63, 87, 503, 562, 632-636, 749-755, 799, 830-832, 966, 986 | "the tail ships as three models"; "**And the DEFAULT `rc_tail` is a LOWER BOUND — it is the t-peak ALONE**"; the 749-755 table (t-peak numbers); "0 / 0 / 0 / 0 at the default `tail_model = TPeak`"; "**`rc_tail` at the DEFAULT is the t-peak only**" | dated D5 notes appended to the existing lines (LN); the table header labelled "(t-peak)". Records cite this file at lines 29–31, 414, 428, 1611 and 1991–2014. |
| `docs/benchmarking/BENCHMARK_PLAN.md` | 83 | row 7 status cell: "shipped `t-peak` 0.6118 / …" | → "`t-peak` (the shipped default until D5; an option since) …"; "`polrad-full` (shipped since 2026-10-27, D5)". |
| `docs/benchmarking/08_channel_priority.md` | 332, §7 revision log | row F: "the shipped `t-peak` R = …"; status "the `TPeak` default (−39 % / −10 %) … candidates for `DECISIONS_PRICED.md`" | row F: "taken: D5 (2026-10-27) — `polrad-full` is the shipped tail"; one revision-log row appended at the end of §7. |
| `validation/benchmarks/README.md` | 204, 238 | row 7: "shipped `t-peak`"; T6 runtime "3.7–4.2 s" | reword row 7; T6 runtime re-measured after step 6 (21 s here, dated). |
| `docs/benchmarking/REPORT.md`, `REPORT.json`, `figures/t4_djangoh_rad_noRad.png` | row 7 `generator` / `reason`, T6 `runtime_s` | **regenerated** by `make_report.py --all --plots` — never by hand. |
| `docs/CONVENTIONS.md` | 20, 55 | "No physics number is hard-coded in two places"; the refusal examples (`rc_sp_tensor_scale` refused on polrad-full "for the OPPOSITE cause") | none — both stay true under D5. |
| `docs/HEPMC3_CONVENTION.md` | 217-218, 239 | weight names only | none (T6's sha256-pinned reference is not touched). |
| `docs/open_items/run_2026-09-03/STATUS.md` (the registry) | after 94 (end of file) | no row names the tail model; row 17 untouched | optional dated blockquote line appended after the 2026-09-23 notes (append = no line shift): "2026-10-27 — D5 taken (`../run_2026-10-27/DECISIONS.md`): `PolradFull` is the shipped tail. Not a registry row. Row 17's ×768 / ×681 / ×289 was measured on the t-peak tail; `run_2026-09-06/phase_B_numbers.md` §B2.5 prices the switch on ΔA_zz at 0.0435 / 0.3073 / 3.254 % of the band half-width." |
| `docs/open_items/run_2026-10-27/DECISIONS.md` | D5 | the supervisor's file | text in §12. |

### 3.6 Dated records and surveys — not rewritten

`docs/open_items/run_*/…` (history; a drift the records gate flags is answered
by a `(as of <sha>)` annotation, never by re-pointing), `docs/references/01_radiative.md`
(pass date 2026-09-16; optional dated note at `docs/references/01_radiative.md:157-158`, "the tree ships `TPeak`
… as an explicit LOWER BOUND"), `docs/references/REFERENCES.md:565, 570`
(about TPeak itself — true), `docs/references/01_radiative_refs_dict_entries.json`,
`docs/benchmarking/01_generators.md`, `06_critic.md`, `docs/references/01_generators-chain.md`.

### 3.7 Checked and found to need nothing

`validation/reference/*.json` (no rc field: `grep` of the nine files); the
polligen port gate; `examples/`; `.github/workflows/ci.yml` (ubuntu-22.04 only,
no `timeout-minutes`); `python/README.md` (its rc weight-name line); `include/lipolgen/event.hpp` (the `rc_weights` comment);
the HepMC3 writer (writes no rc provenance — §10 R7).

## 4. Pin decisions — the rule, and the totals

1. **Re-take on the new default** (dated comment "re-taken 2026-10-27, D5:
   `PolradFull` is the shipped tail") where the assertion IS the default's
   identity: `test_rc.py:263, 612, 650-651`; the knob-provenance `rc` and
   `rc_tail_model` cells; the new C++ D5 case; `test_rc_pipeline.cpp:247`; the
   documents' default columns; T6 (it tests the user's own path).
2. **Pass `t-peak` explicitly** where the purpose is t-peak physics or a
   comparison against a t-peak-era number: `rc_model()` and the t-peak cases of
   `test_rc.cpp`; `test_c0_shape…`, `test_the_quasi_elastic_sp_column…`, the
   bases of `test_polrad_full…` and `test_tail_model_knob…`; the knob cells
   that gate the t-peak refusals (`rc_sp_tensor_scale`, `rc_m_lepton`). The
   t4 harness already names both models.
3. **Pass `t-peak` for runtime** where the subject is model-independent and
   the full-grid tail is built only incidentally: the rest of `test_rc.py`'s
   20, `test_rc_low_x_anchor.py`, the knob cells of the band and tail
   sub-knobs, `test_cli.py`'s untyped-flags run. Dated comment naming D5 and
   the measured ~15 s.
4. **Leave on the default** the model-independent C++ cases (T14, T17, B6,
   TensorAll, `test_rc_pipeline.cpp` T6 / T18b) — 2.3–2.7 s each on their
   small grids.

**Cost** (from §2.1 and §2.3; the phase-E full pytest was 467 s): this rule
adds ≈ 1.5–3 min to the pytest gate (rc_run +15 s; knob provenance ≈ +81 s:
four `rc`-cell variants, one cached `RC` base for the tensor-thirds spec, the
`n_eta` 64 variant at about half cost; T6 +15 s; one build saved in
`test_polrad_full…` −15 s; ≈ +96 s unloaded),
≈ +25 s to `lipolgen_tests`, ≈ +15 s to `LIPOLGEN_BENCH=1 -k bench` and to
`make_report.py --all`. Re-taking every default-dependent build on the new
default instead would add ≈ +17 min to pytest (67 builds) and ≈ +85 s to C++.

## 5. Benchmark rows, and how their declared rules read

- **Row 7, `t4_djangoh_rad_noRad` — verdict unchanged (PASS), numbers
  bit-identical, labels move.** The declared rules name the **models**, not the
  default (quoted as the block writes them):
  "PASS(polrad-full) iff |R_b - 1| <= 0.03 + 2 delta_b" (`validation/benchmarks/t4_djangoh_rad_noRad.py:171`),
  "PASS(t-peak) iff R_b <= 1 + 0.03 + 2 delta_b" (`validation/benchmarks/t4_djangoh_rad_noRad.py:178`),
  "The row PASSES iff both hold in all four bins" (`validation/benchmarks/t4_djangoh_rad_noRad.py:182`).
  Only the block's descriptions ("the SHIPPED default" at `validation/benchmarks/t4_djangoh_rad_noRad.py:172`,
  "opt-in" at `validation/benchmarks/t4_djangoh_rad_noRad.py:148` (`RcTailModel.PolradFull`), "the price of the shipped default" at
  `validation/benchmarks/t4_djangoh_rad_noRad.py:180-181`) name the default, and they are history. So **no rule has to follow
  the new default and nothing becomes a recorded sub-row**: both models stay
  in the verdict as declared; the t-peak sub-rows' `shortfall_1_minus_R` now
  reads as the price of the t-peak option. The row reads "polrad-full
  (shipped), t-peak (option)" through the generator string, the reason and the
  figure — all outside the block. **Addendum text** (append after the existing
  ADDENDA, whose last item is `validation/benchmarks/t4_djangoh_rad_noRad.py:207-212`): *"2026-10-27 (decision D5, docs/open_items/run_2026-10-27/
  DECISIONS.md): `PolradFull` is the shipped `rc_tail` default and `TPeak` an
  option. The block above is hashed and stays as written on 2026-10-26, when
  t-peak was 'the SHIPPED default' and polrad-full 'opt-in'. Its two PASS rules
  name the models, not the default, so the verdict is computed exactly as
  before; 'the price of the shipped default' (block (e)) now reads as the price
  of the t-peak option. The REPORT generator string, the reason and the figure
  legend say 'polrad-full (shipped)' and 't-peak (option)' since that date."*
- **T6, `t6_hepmc3_convention` — verdict unchanged** (23 / 23 measured on the
  polrad-full rc file); runtime 6.06 → 21.07 s; no rule names the tail model.
- **The rows `tiers.py --files` maps to `rc.*` / `cli.py`** (t2_*, t3_*,
  t4_pythia_ep_closure, t6): no harness other than t4_djangoh and t6 builds an
  `RcModel`; re-run by `make_report.py --all`, expected: no change.
- **No other row** names the tail model or runs `--rc` (`t5_sum_rules`
  included).

## 6. Published numbers that move

| number (site) | before (t-peak) | after (default) | regenerated by |
|---|---|---|---|
| row 7 REPORT generator / reason, figure legend and colours (`REPORT.md:39, 345`, `REPORT.json`, `figures/t4_djangoh_rad_noRad.png`) | "t-peak (shipped)" | "polrad-full (shipped)", "t-peak (option)" | script (`make_report.py --all --plots`), after the harness edit |
| T6 runtime (`REPORT.*`; the `t6_hepmc3_convention.py` row, `validation/benchmarks/README.md:238`) | 3.7–6 s | ≈ 21 s (loaded) | script / hand (README) |
| default-path setup cost (`USAGE.md:2581-2590`, `PHYSICS_CHANNELS.md:393`, `rc.hpp:1331-1338`, `OPEN_ITEMS_SOLUTIONS.md:828`) | 0.106 s | 14.7 s (0.106 s on `t-peak`) | hand |
| clipped tail nodes at the default (`USAGE.md:2556-2559`, `PHYSICS_CHANNELS.md:391`, `OPEN_ITEMS_SOLUTIONS.md:830`) | 0 % | 0.386 % (0.456 / 0.124 / 0 %); events 0 | hand |
| ⟨rc_tail⟩ at the default (USAGE, PHYSICS, rc.hpp — already published per model) | 1.021778527 (event, P_z 0.7) / 1.02217080 (cell) | 1.029702912 / 1.03020526 | hand: relabel which is the default |
| quasi-elastic share 22 / 73 / 99.9 % (rc.hpp:47-51, 1544-1546, 1956; rc.cpp:1558; cli.py:440-441, 1658-1662; bindings.cpp:857; USAGE:2212, 2341, 2248-2249; PHYSICS:387, 841; OPEN_ITEMS:562, 799, 966) | 22.0 / 73.3 / 99.90 % | 38.7 / 99.81 / 100.00 % | hand (code strings rebuilt) |
| (1/6)σ^el_T/σ^el_U box (`USAGE.md:2253-2259`, `PHYSICS:386`), README `--rc-c0-shape` row (`README.md:54`) | §2.5 t-peak rows | §2.5 polrad-full rows | hand (label + add) |
| w_tail − 1 at Q² = 5 (`PHYSICS:386`) | 5.18398e−04 / 1.32298e−06 / 6.8624e−08 | 6.6287e−04 / 1.9596e−04 / 4.8836e−04 (table route) | hand |
| σ^el_U shape ratios ×1.011 / ×1.022 / ×1.104 / ×3.03 and ΔA_zz −8 / +0.6 / +37 % (`USAGE.md:2264-2266`), `--rc-qe-tensor-scale` ΔA_zz price (`USAGE.md:2249`), the fq-scale sign change (`USAGE.md:2244`), the ΔA_zz spreads 1.680e−07 / 1.878e−07 (`PHYSICS:387`), OPEN_ITEMS 749-755 table | t-peak | **not re-measured** — label "(t-peak)" | hand |
| registry row 17's ×768 / ×681 / ×289 | t-peak | not re-measured (§B2.5's price of the switch is 0.04–3.3 % of the band) | optional dated note |

No reference JSON, no `validation/benchmarks/data/*` file and no tolerance moves.

## 7. The citation gates: how a line shift is handled, and this design's answer

**What the gate does** (read in `validation/check_physics_channels_links.py`,
its docstring `validation/check_physics_channels_links.py:1-247`
and code `validation/check_physics_channels_links.py:1485-1506` (`nearest_decl`), `:2025-2135`, `:2800-2813` (where `--fix` writes the `remap`); run
from the repository root, read-only unless `--fix` or `--record-ranges`):

- **Named point citation** (`` `path:line` `name` ``, strict): the cited line
  must declare the name (S1/S2) or be a pinned use site (S3). A shifted
  declaration is reported broken; **`--fix` rewrites the line number to the
  nearest declaration** (`nearest_decl`; refuses an exact tie; never touches an
  `ALLOW` entry). Hazard (inferred from the code): under a *negative* shift the
  nearest declaration of an overloaded name can be the other overload — e.g.
  `rc_weight_name` (`include/lipolgen/rc.hpp:1858`, `:1862`), `tensor_fraction`
  (`include/lipolgen/rc.hpp:1990`, `:1996`), `weights`
  (`include/lipolgen/rc.hpp:2092-2093`) — and the gate then accepts it.
- **Pinned use site / range** (sidecar `validation/physics_channels_ranges.json`,
  219 entries: sha256 of the line or block, its first/last line in clear): a
  block found byte-identical elsewhere is **MOVED** — `--fix` rewrites the
  document's citation and remaps the sidecar key, carrying the sha256 over; a
  block found nowhere is **EDITED in place** — a human re-reads the row, then
  `--record-ranges` (which re-applies R5 and refuses a block the citing
  sentence does not evidence). `RANGE_ALLOW` keys are not rewritten; `--fix`
  prints a note.
- **Unnamed point citations** (rule B, `UNNAMED_ALLOW` pins keyed on
  (document, citing line, `path:line`)) and **`ALLOW` entries** (keyed on
  (path, line, name)): **no `--fix`**; a moved pin is a hand edit of the gate
  script.
- **External PYTHIA citations**: no `--fix`.
- **Records** (`--records` / `--records-only`, rule R6, relaxed): **no `--fix`**;
  the answer to a drift is a dated `(as of <sha>)` / `(anchor read …)`
  annotation in the record. On `ca9b23a`: 838 citations in 69 records,
  0 broken; with this run's records added (this file included, checked
  2026-10-27): 1230 citations in 71 records, 0 broken.
- `--loose` (±2 lines; skips R4, S3, D) is "NOT the gate".

**Exposure of the D5 edits** (counted today): PHYSICS_CHANNELS.md cites
`rc.hpp` 89×, `rc.cpp` 43×, `pipeline.cpp` 43×, `pipeline.hpp` 98×,
`bindings.cpp` 17× (+4 ranges from SPIN32_FINITE_GAMMA.md), `cli.py` 71×,
`__init__.py` 13×, `USAGE.md` 8×, `test_rc.cpp` 2×, `README.md` 1×. Sidecar
entries in those files: `rc.hpp:305-325, 383-465, 414-421`; `rc.cpp:1435, 1481,
1506, 1525`; `pipeline.cpp:557, 811, 1149, 1537, 1680, 1710, 1781, 3749`;
`cli.py:1636`; `__init__.py:422, 634, 694`; `USAGE.md:11-66, 163-533,
1978-2169, 2059-2065, 2768-2775, 2795-2992`; `bindings.cpp:2187-2199, 2227-2229,
2341-2342, 2488-2501`. Hand-keyed exemptions: `ALLOW` `pipeline.cpp:847`
`coherent`, `pipeline.cpp:1877` `optics_lumi_factor`, `rc.cpp:1532`
`is_tagged_channel`; `UNNAMED_ALLOW` (PHYSICS:176 → `pipeline.cpp:1116`),
(PHYSICS:420 → tests/test_rc.cpp line 518), (PHYSICS:90 → README.md line 128);
`RANGE_ALLOW` (PHYSICS:451 → `USAGE.md:2795-2992`). 64 exemption keys carry a
PHYSICS_CHANNELS.md citing line, 26 of them below line 386. Ungated but
cited: README.md line 96 → tests/test_rc.cpp line 1876; docs/benchmarking/02_data_nucleon_deuteron.md
line 299 and docs/references/01_nucleon-deuteron-data.md line 429 → tests/test_rc.cpp line 374,
which is a blank line today (those two citations are already stale; not D5's to fix).
Records cite `test_rc.cpp:1257, 1267, 1608, 1630` and `OPEN_ITEMS_SOLUTIONS.md:1611, 1991-2014`.

**This design's answer: line-count-neutral edits** in every file above —
`rc.hpp`, `rc.cpp`, `pipeline.cpp`, `bindings.cpp`, `cli.py`, `__init__.py`,
`tests/test_rc.cpp` **above line 1877** (append new test cases at the end of the
file), `README.md`, `USAGE.md`, `PHYSICS_CHANNELS.md`, `OPEN_ITEMS_SOLUTIONS.md`
— with no edit inside a fingerprinted block (none of the edits in §3 lands in
one: checked against the sidecar list), and no edit of a cited declaration line
that removes its name (`rc.hpp:1503` keeps declaring `tail_model`). Prose that
needs more words is appended to an existing markdown line; code comments are
rewritten within their own line count. Then the three gate modes stay at
today's counts with nothing to fix: strict **1240 references / 95 ranges / 7
external / 0 broken / 24 allow-listed** (PHYSICS), **19 / 115 / 0 / 0 / 10**
(SPIN32), **0 / 0 / 8** (PYTHIA_BRIDGE); records **0 broken** (the citation
count grows with this run's own records). Check
line-neutrality before running them: `git diff --numstat -- <files> | awk '$1 != $2'`
must print nothing for the LN files. **If an edit cannot be neutral:** run the
strict gate, then `--fix` (named points, MOVED ranges and use sites), review the
rewritten citations of every overloaded name, hand-edit the `ALLOW` /
`UNNAMED_ALLOW` / `RANGE_ALLOW` keys the gate names, re-read and
`--record-ranges` any EDITED block, then `--records-only` and annotate each
newly broken record citation `(as of <base sha>)`.

## 8. Ordered implementation steps (I5), each with its gate

0. **Baseline on the tree I5 starts from** (after whichever of I4 / I6 / I7
   landed): the three A/B column hashes with the commands of §2.4 / §2.7
   (`$SCR/npz_hash.py`, copy it into the I5 scratch); `make_report.py --check`;
   the three gate modes (counts); save the base report with
   `git show HEAD:docs/benchmarking/REPORT.json > $SCR/REPORT_base.json`.
1. **The switch and the C++ tests.** `rc.hpp:1503` → `PolradFull`.
   `cmake --build build -j2 && ./build/lipolgen_tests -sf="*test_rc*"` —
   record every failing case; show each is a t-peak pin (§3.2), then pin
   `TPeak` (LN above 1877); append the D5 case; `test_rc_pipeline.cpp:247`.
   **Gate:** `./build/lipolgen_tests` whole suite, 0 failed; runtime against
   the baseline (+≈25 s expected).
2. **Library text and the CLI default.** `cli.py:817` (derived default), the
   help, the banner; `__init__.py`; `bindings.cpp` docstrings; `pipeline.cpp` /
   `rc.cpp` strings and comments; `rc.hpp` comments — all LN. **Gate:**
   rebuild; `./build/lipolgen_tests` (unchanged); the **A/B check** of §9 V2;
   `python3 -m pytest python/tests/test_rc.py python/tests/test_rc_low_x_anchor.py
   python/tests/test_knob_provenance.py python/tests/test_cli.py -q` — the
   failures must be exactly the re-take / pin list of §3.3.
3. **Python tests** per §3.3 / §4. **Gate:** those four files pass; then the
   full `python3 -m pytest python/tests -q` (0 failed; wall time against the
   pre-I5 baseline, +1.5–3 min expected).
4. **Harness labels** (t4 outside the block + the addendum; `_plots.py:56`;
   `test_bench_djangoh_eltail.py` wording and the optional label assert).
   **Gate:** `python3 validation/benchmarks/t4_djangoh_rad_noRad.py` exits 0
   (the block hash check passes) and prints R unchanged; `LIPOLGEN_BENCH=1
   python3 -m pytest python/tests -q -k bench` (0 failed).
5. **Documents** per §3.5 (LN). **Gate:** `git diff --numstat` neutrality;
   `python3 validation/check_physics_channels_links.py` (strict),
   `--records-only`, `--audit-ranges` — all exit 0 with today's counts (§7);
   the grep audit of §9 V8.
6. **The report.** `python3 validation/benchmarks/make_report.py --all --plots
   --compare $SCR/REPORT_base.json --out-dir docs/benchmarking`; expected
   changes: row 7's generator / reason strings and figure, T6's runtime, nothing
   else (no status, no number); then `make_report.py --check`; update
   the runtime of the `t6_hepmc3_convention.py` row, `validation/benchmarks/README.md:238` (dated).
7. **SPDX and the record.** `python3 validation/check_spdx_headers.py`; write
   `phase_I5_polradfull.md` (every measured number with its command); the
   supervisor commits and pushes.

## 9. Verification plan

- **V1 — the default.** `python3 -c "import lipolgen as lg; L = lg._lipolgen;
  print(L.RcOptions().tail_model, L.rc_tail_model_name(L.RcOptions().tail_model))"`
  → `RcTailModel.PolradFull polrad-full`; `cli.DEFAULTS["rc_tail_model"] ==
  "polrad-full"`; an `--rc tensor-band` npz has `meta["rc_tail_model"] ==
  "polrad-full"`.
- **V2 — A/B bit identity (the decisive check).** On the I5 tree, with the
  step-0 commands: the default `--rc tensor-band` run's column hash equals the
  step-0 explicit `polrad-full` hash; `--rc-tail-model t-peak` equals the
  step-0 default hash; `--rc off` equals the step-0 `--rc off` hash. `meta`
  differs only in the keys of §2.7 (rc on) and in
  `knob_provenance["rc_tail_model"]["value"]` (rc off). (On `ca9b23a` itself
  the hashes are those of §2.7.)
- **V3 — t-peak regression.** Every pin moved to explicit `TPeak` passes with
  its number unchanged (no tolerance touched).
- **V4 — row 7.** R pins at rel 2e−7 pass; `design_block_sha256()` unchanged;
  PASS; the new labels in `REPORT.md`.
- **V5 — T6.** 23 / 23; `runtime_s < 30`.
- **V6 — gates.** The counts of §7, 0 broken in all three modes; SPDX OK.
- **V7 — report.** `--compare` shows only the expected strings and runtimes.
- **V8 — no stale default claim.** `git grep -n -i -E "(t-peak|TPeak)[^.]{0,40}(DEFAULT|shipped)|(DEFAULT|shipped)[^.]{0,20}(t-peak|TPeak)|polrad-full[^.]{0,30}(opt-in|not the default)"
  -- include src python README.md docs/USAGE.md docs/PHYSICS_CHANNELS.md
  docs/OPEN_ITEMS_SOLUTIONS.md docs/benchmarking/BENCHMARK_PLAN.md
  docs/benchmarking/08_channel_priority.md validation/benchmarks` — every hit
  carries the D5 dating ("until 2026-10-27", "since D5") or is inside the
  hashed block.
- **V9 — the re-taken document numbers** reproduce from the recipes of §2.4–2.5
  (`m3`, `m7`, `m8`), copied into the I5 record with their commands.
- **V10 — the run's phase V** (four adversarial lenses) reads the diff against
  this plan.

## 10. Risks

- **R1 — runtime.** Every inclusive ⁶Li `--rc tensor-band` run pays ≈ 15 s of
  setup (21–24 s with `vmc-ft` or a Python form factor; 23–25 s for the proton
  harness), with no cache — batch jobs pay it per job. T6 rises to ≈ 21 s
  against its 30-s pin on a loaded machine. *Fallback, if the pin fails on the
  integration machine:* pass `--rc-tail-model t-peak` in T6's rc file with a
  dated comment (T6 tests HepMC3 weight names, which no tail model changes).
  A table cache or a threaded build is a separate design (not D5).
- **R2 — platform.** `PolradFull` is refused where `long double` is no wider
  than `double` (MSVC; Apple arm64): there `--rc tensor-band` fails by default,
  with a message naming `--rc-tail-model t-peak` (`rc.cpp:1401-1409`). On
  aarch64 Linux `long double` is a software quad — the build time there is
  **unmeasured**. CI (ubuntu-22.04 x86_64) is not affected.
- **R3 — the evidence's reach.** Row 7 checks the proton (quasi-elastic) slot
  only: one beam, one shared dipole, no Pauli factor, no coherent ⁶Li term, no
  tensor term; bin 4's −1.34 % is unexplained. The (x, y) table's +1.9 … +2.5 %
  (DECISIONS_PRICED §2, still open) now sits on the default path.
- **R4 — misattributed numbers.** Every RC-tail number published before D5
  without a model label is a t-peak number; a document not updated attributes
  it to the default. The shares and w_tail − 1 move by large factors at high x
  (×148 at x = 0.1, ×7100 at x = 0.3; small in absolute terms). Mitigated by
  V8 and the labels of §6.
- **R5 — clipping.** The default now clips 0.386 % of tail nodes (six accepted
  cells at x = 0.955, 4.96e−07 of σ; 0 / 200 000 events) — every "nothing clips
  at the default" sentence becomes false.
- **R6 — API behaviour.** At the default, `n_eta < 64` is refused (was legal
  from 8), `m_lepton` is read (was refused), `meta["rc_n_eta"]` means tanh-sinh
  nodes per panel; on tagged channels `--rc-tail-model t-peak` becomes the
  refused value. User scripts that set `n_eta` low for speed break loudly.
- **R7 — HepMC3 provenance.** The HepMC3 files carry `rc_tail` without the rc
  mode or tail model (`hepmc_writer.cpp` writes neither; the npz `meta` does):
  the same command line gives a different `rc_tail` before and after D5 with
  nothing in the file saying why. Out of D5's scope (it would change
  `docs/HEPMC3_CONVENTION.md`, T6's sha256-pinned reference); filed as an open
  item.
- **R8 — citation drift.** A non-neutral edit breaks up to 26 exemption keys,
  the `pipeline.cpp` / `rc.cpp` / `test_rc.cpp` / `README.md` pins and the
  records; `--fix` can re-point an overloaded name onto its neighbour (§7).
- **R9 — interactions** (§11): D7's Born moves every Born-dependent RC number
  again; D4 changes the sampled φ, which voids the §2.7 hashes.
- **R10 — coverage of the new default.** With rule 3, fewer Python tests run
  the full-grid polrad-full tail; the default stays covered by the `rc_run`
  chain, the knob `rc` / `rc_tail_model` / `n_eta` cells, T6, the C++ D5 case,
  and T14 / T17 / T6 / T18b in C++.
- **R11 — T8(0)** must stay on `TPeak` (§2.8); a polrad-full sign gate without
  the x = 1 clause can be added if wanted.

## 11. Interactions with D4, D6, D7 (`DECISIONS.md`)

- **D7 (grid SF backends default)** changes the Born that `rc_tail` = 1 +
  σ_tail/σ_Born divides by: w_tail − 1, ⟨rc_tail⟩, the A_zz dilution, the node
  clip fraction move again; (1/6)σ^el_T/σ^el_U and the QE shares do not (ratios
  of tail pieces). Row 7's R carries the Born only through the table's separate
  ln-interpolation of tail and Born — D7's plan must check whether row 7's
  rel 2e−7 pins move (the harness reads the default backend). Whichever of I5 /
  I7 lands second re-takes the Born-dependent RC numbers; every RC number in the
  documents should name its backend.
- **D4 (Jaffe–Manohar Δ)** changes the generated cos 2φ term, hence the
  sampled φ: the §2.7 hashes are void after I4; re-baseline at step 0. The RC
  band's `RcScope::TensorAll` prices the Δ term; the tail is untouched.
- **D6 (VMC cluster waves)** touches only the tagged channels, where
  `rc_tail ≡ 1`: no interaction with the tail model.

## 12. Decision text for `DECISIONS.md` D5

> **D5 — `RcTailModel::PolradFull` becomes the shipped `rc_tail` default;
> `TPeak` stays an option.** *(decision 2: "PolradFull is wanted.")*
> `RcOptions::tail_model` and `--rc-tail-model`'s default move from `TPeak` —
> POLRAD Eqs. (37)–(39), (43), the t-peak alone, a lower bound — to
> `PolradFull` — POLRAD Eq. (18) + Appendix B + Eq. (A.4), one exact τ_A
> quadrature carrying all three peaks with the lepton mass. The evidence is
> row 7 (`../run_2026-10-26/phase_B3_djangoh.md`): on a proton at e 18 × p 275,
> straight from the quadrature, `polrad-full` is 1.0004 / 1.0010 / 0.9995 /
> 0.9866 of DJANGOH 4.6.21's exact O(α) elastic tail at Q² = 1–10 / 10–100 /
> 100–10³ / 10³–10⁴ GeV², where `t-peak` is 0.599 / 0.878 / 0.989 / 0.985 —
> one beam, one shared dipole, the proton (quasi-elastic) slot only; the
> coherent ⁶Li and tensor terms remain unchecked externally, bin 4's −1.3 % is
> unexplained, and the (x, y) table's +2 % is a separate open candidate.
> **What does not move:** `--rc off`, the generator's default, keeps every
> column bit for bit (one `meta` string changes: the not-read
> `knob_provenance` row `rc_tail_model` names `polrad-full`); no reference JSON
> carries an rc field; the polligen port gate reads no RC; tagged and coherent
> channels build no tail. **What moves**, measured 2026-10-27 on ⁶Li config 1
> with the toy backends: the tail-table build 0.23 s → 14.7 s per `--rc
> tensor-band` run (21 s with `--rc-c0-shape vmc-ft`), with no cache; the mean
> dilution ⟨rc_tail⟩ − 1 ×1.36 (cell-σ-weighted 2.217 → 3.021 %); the run-level
> A_zz with the tail −0.05 % (its dilution 0.385 → 0.436 %); tail nodes
> clipped 0 → 0.386 % (events 0); the quasi-elastic share of the tail at
> Q² = 5 GeV² 22 / 73 / 99.9 % → 39 / 99.8 / 100 % at x = 0.01 / 0.10 / 0.30.
> Every RC-tail number published before 2026-10-27 without a model label is a
> t-peak number and reproduces with `--rc-tail-model t-peak`; the documents
> now say which model each number is. Tests that pin t-peak physics pass
> `t-peak` explicitly; tests of the default's identity are re-taken with a
> dated comment; row 7's verdict is unchanged (its declared rules name the
> models) and it now reads "polrad-full (shipped), t-peak (option)". Platforms
> whose `long double` is no wider than `double` refuse the new default and must
> run `--rc-tail-model t-peak`. Reversible: one line (`rc.hpp`
> `RcOptions::tail_model`), the CLI default that reads it, and the re-taken
> pins. Registry row 17 untouched; no registry row is added.

## 13. What this design did not establish

- The C++ suite under the flipped default was **not run** (no rebuild): which
  `test_rc.cpp` cases fail is predicted from reading (§3.2); step 1 decides.
- The build time on aarch64 Linux and the refusal on Apple arm64 / MSVC were
  not measured.
- Not re-measured on `polrad-full`: the σ^el_U shape ratios of the C0 band
  (incl. `test_c0_shape`'s 1.1043), the fq-scale sign change, the qe-tensor-scale
  ΔA_zz price, the ΔA_zz spreads of `PHYSICS_CHANNELS.md:387`, row 17's
  ×768 / ×681 / ×289 ordering — they stay labelled t-peak numbers.
- The published "w_tail − 1 = 5.18398e−04" recipe was not reproduced exactly
  (−0.03 %, §2.5).
- No event-level A_zz (5000 events cannot resolve it); the A_zz change is the
  cell-level computation of §2.4.
- The measurements ran on a shared, loaded machine: times are ±30 % (the same
  build measured 14.7–15.3 s).

## Adversarial check (2026-10-27)

The adversarial check of this design, run on 2026-10-27 against `ca9b23a`
with nothing in the repository edited, reproduced in the checker's own words
(nothing trimmed, no number changed). Its file and line references are the
checker's, at `ca9b23a`, and are written as plain text, not as citations the
records gate checks.

### Verdict

SOUND IN ITS CORE, NOT READY TO IMPLEMENT AS WRITTEN. These parts check out: the one-line switch and the CLI's second definition of the default (cli.py:817), the line-neutral strategy and the gate counts (1240/95/7/0/24, 19/115/0/0/10, 0/0/8 reproduced), the hashed t4 block (its PASS rules name the models, so row 7 cannot flip), the polligen and reference-JSON claims, the T8(0) finding, and the costs. I re-measured the costs unloaded: the build goes 0.157 s → 12.35 s, and T6 goes 6.07 → 19.11 s with 23/23 clauses passing. Four things must be fixed before I5: (1) the test_knob_provenance.py split is wrong. The proposed RC_TPK base makes 57 matrix cells that run today skip silently on the tagged and coherent specs, while pytest stays green. (2) The test_rc_low_x_anchor.py helper edit breaks A5, which builds tagged-alpha pipelines. (3) About 15 code and document sites that quote t-peak numbers as 'the whole tail' or as the default are missed: the ×768 ordering, −1.769797e−07, 'Nothing clips on the tail any more', the QE shares in PHYSICS:388/390, the 0.47/0.87 QRT cut, and the C0 numbers in __init__.py and cli.py. (4) The proposed claim that '--rc-tail-model t-peak reproduces every earlier number' becomes false once I7 lands. The platform risk is also understated: on a platform whose long double is no wider than double, the tagged and coherent band becomes unreachable.

### Sites the design missed (15)

| file | line | what |
|---|---|---|
| `include/lipolgen/rc.hpp` | 242-244 | HONEST FLAGS: 'the band still leads the WHOLE radiative tail by x768 / x681 / x289 there'. These ratios were measured on the t-peak tail; per run_2026-09-06 B2.5 the polrad-full tail's ΔA_zz differs by 5.9e-08 at x = 0.01. Label them or re-measure. |
| `docs/PHYSICS_CHANNELS.md` | 384, 385 | Row 384: 'The band still leads the whole radiative tail by ×768 / ×681 / ×289'. Row 385: 'against the whole tail's ΔA_zz = −1.769797e−07, i.e. the band leads by ×768 / ×858 / ×1086'. Both are t-peak numbers presented as 'the whole tail'. |
| `docs/USAGE.md` | 2273-2289 | '1.769797e−07 from the entire tail — a factor 768 over the tail', plus the correction box that restates −1.769797e−07 and 768. |
| `docs/OPEN_ITEMS_SOLUTIONS.md` | 806-822, 836-844, 907 | 'ΔA_zz from the whole tail … −1.769797e−07 … ×767.6'. The run-level ΔA_zz spread rows (qe_suppression, qe_kf_gev, c0_shape, qe_tensor_scale) all start from the t-peak tail's −1.769797e−07. Line 907: 'leads the whole radiative tail … by ×768 / ×681 / ×289'. The design lists 830-832, not these lines. |
| `docs/OPEN_ITEMS_SOLUTIONS.md` | 853-855 | '**Nothing clips on the tail any more**' becomes false at the D5 default. I measured 0.399 / 0.386 / 0.360 % of nodes clipped on 6Li configs 0/1/2, with 0 of 50 000 events clipped on each. The design's R5 lists USAGE, PHYSICS:391 and OPEN_ITEMS:830, not this sentence. |
| `docs/PHYSICS_CHANNELS.md` | 388, 390 | Row 388: 'the quasi-elastic term — **22 % / 73 % / 99.9 %** of `rc_tail`'. Row 390: '(99.9 % of the inclusive tail at x = 0.30)'. These are t-peak shares; the design lists only 387 and 841. |
| `include/lipolgen/rc.hpp` | 30-35, 52-57 | The header says rc_tail is 'Implemented as the additive weight `rc_tail` from POLRAD's t-peak closed forms, Eqs. (37)-(39) + (43)', which is no longer the default. It also quotes 'sigma^q_U/sigma^el_U = 0.28 / 2.7 / 1000' (t-peak). I reproduced polrad-full as 0.631 / 533 / 3.68e6. |
| `python/lipolgen/__init__.py` | 337-341 | RC_C0_SHAPES comment: 'ho/vmc-ft = -5.0943e-4 / -5.4705e-4 at x = 0.01, +1.5595e-4 / -4.0702e-5 at x = 0.10 …'. These are unlabelled t-peak numbers; the design's own §2.5 calls them t-peak, but it lists only lines 347-378 of this file. |
| `python/lipolgen/cli.py` | 353-355, 301-302, 458-460 | --rc-c0-shape help: 'The tensor fraction of the elastic tail at x = 0.1 is +1.6e-4 on one edge and -4.1e-5 on the other' (t-peak; polrad-full is +4.5e-5 / -1.6e-4). The --rc help's t-peak paragraph runs past 290-300 to line 302. The --rc-sp-tensor-scale help gives only the t-peak refusal reason ('the t-peak-only tail computes no s/p peaks'). |
| `python/lipolgen/cli.py` | 1651-1654 | Banner: 'clipped … of tail NODES (by y-band …; the y > 0.9 one is where Y+ ~ 1/(1-y) bites)'. At the new default the clipped nodes sit at y < 0.5 and 0.5–0.9 and none at y > 0.9, yet every default --rc tensor-band run prints this line. |
| python/bindings.cpp; docs/USAGE.md | bindings 3637-3640; USAGE 2251 | qe_kf_gev docstring and table row: 'It cuts the QRT to 0.47 at x = 0.01 and 0.87 at x = 0.1'. This is a t-peak number; the design pins C++ Q9 to TPeak for that reason but leaves the published text unlabelled. rc.hpp 1616-1621 gives the same t-peak reasoning for the default. |
| `docs/USAGE.md` | 2250 | --rc-sp-tensor-scale row: '**Refused unless `--rc-tail-model t-peak+ll`** — the t-peak-only tail computes no s/p peaks'. At the default the reason is now the opposite one (polrad-full computed the term). |
| `tests/test_rc.cpp` | 1330, 1355, 1550, 2716, 2727, 2878 | Comments and SUBCASE names: 'the QRT really is a large part of the default tail', 'the default is 0 and the SHIPPED tail cannot see it', '(a) the DEFAULT is untouched: the s+p columns are identically 0'. Once rc_model() is pinned to TPeak these test TPeak, not the default; the design rewords only 1544 and 2724. |
| `python/tests/test_rc.py` | 233-234, 342-343 | Docstrings 'the QRT is 73-99.9 % of the tail' and '73 % of rc_tail at x = 0.1 and 99.9 % at x = 0.3' (t-peak shares). |
| `validation/benchmarks/t4_djangoh_rad_noRad.py` | 318-320 | 'rc.hpp says ~10 s unloaded'. This cross-references the rc.hpp cost sentence that the design rewrites to 14.7 s. |

### Statements found wrong (8)

1. test_knob_provenance.py split, design §3.3: 'RC_TPK = dict(RC, rc_tail_model="t-peak") for the model-independent band and tail sub-knob cells (rc_delta_low_x … rc_tail_max, the NO_QE cells)'. VARIANTS is swept over every spec (`CELLS = [(s, v) for s in SPECS for v in VARIANTS if s.full or v.knob in PLAN_AXIS_KNOBS]`, :535-536). After D5, an explicit t-peak is a non-default tail knob, and validate() refuses it on tagged and coherent channels (pipeline.cpp:743, 767-783). The cell base is therefore refused and the cell hits `pytest.skip("the BASE of this cell is refused …")` (:547-549). I measured this on ca9b23a with the mirror case (a non-default polrad-full base): 57 cells that run today would skip, 18 each on coherent, tagged-6Li-alpha and tagged-d-p, and 3 on tagged-7Li-alpha pe=0. Only 6 of them have a polrad-full twin cell that takes over. These are exactly the class-4 'RC sub-knobs on the tagged channels' checks the file exists for. In addition, the 'rc_with_tail: true -> false' cell on inclusive-6Li would pass through the refused branch for the wrong reason: with with_tail = false, the explicit t-peak is a refused tail knob. The split must depend on the channel.
2. test_rc_low_x_anchor.py, design: '_pipeline() helper, 10 builds; A1–A5 are band numbers → pass t-peak for runtime (`rc_tail_model="t-peak"` in the helper)'. A5 calls the same helper with `channel="tagged-alpha"` (lines 174-175), and make_config validates. Mirror reproduced: 'RuntimeError PipelineConfig: --rc-tail-model is a knob on the RADIATIVE TAIL, and rc_tail == 1 exactly on this run (channel tagged-6Li-alpha)'. A5 would error.
3. Proposed published text: README:51 '`--rc-tail-model t-peak` reproduces every earlier number'; rc.hpp point (6) 'RC-tail numbers published before D5 reproduce with tail_model = TPeak'; §12 'Every RC-tail number published before 2026-10-27 without a model label is a t-peak number and reproduces with `--rc-tail-model t-peak`'. These become false in this same run once I7 lands. w_tail−1, ⟨rc_tail⟩, the ΔA_zz dilution and the node clip all divide by the Born, which D7 moves; the D7 plan §3.8 measures row 7's R moving, polrad-full 1.0199→1.0213. Event-level numbers would also need D4 undone, and there is no flag for that. The recipe must name `--unpol-sf toy --pol-sf toy` as well, or the claim must be limited to Born-free ratios.
4. §0: 'Nothing moves outside `--rc tensor-band` on the inclusive channel'. Not quite. On tagged and coherent `--rc tensor-band` runs, meta['rc_tail_model'] changes 't-peak' → 'polrad-full' (bindings.cpp:871-872 writes it on every rc-on run), even though no tail runs there, and the knob rows change with it. On narrow-long-double platforms those runs are refused outright (next item).
5. R2: 'there `--rc tensor-band` fails by default, with a message naming `--rc-tail-model t-peak` (rc.cpp:1401-1409)'. Both PolradFull preconditions fire whatever the channel (validate() pipeline.cpp:616-631 sits before the channel rule; RcModel rc.cpp:1390-1425 runs before tail_applies_ is used). On tagged and coherent channels, following the message is then refused by the tail sub-knob rule (pipeline.cpp:743). So on MSVC and Apple arm64 the tagged band becomes unreachable, and C++ tagged cases such as T2 (test_rc.cpp:2212) and the clamp case (:2307) would throw. The design does not see this deadlock.
6. §7: line-count-neutral edits mean 'the strict gate, the records gate and R5 stay green with no --fix'. That is incomplete. B2 (strict), R5 (`--audit-ranges`, part of the design's own step-5 gate) and R6 are content rules on the citing sentence or the cited line. The PHYSICS rows the design rewrites hold range citations whose evidence is a single quoted phrase within 200 characters, and that phrase must survive the rewrite: row 393 → USAGE.md:1978-2169 via `qe_suppression`; row 390 via `applies`; row 389 → rc.hpp:383-465 via `tests/test_rc.cpp`; row 384 → rc.hpp:305-325 via `RC_DELTA_HIGH_X`. The measured exposure is small, and --audit-ranges would catch a failure.
7. test_rc.py `test_sp_tensor_scale_is_refused_when_the_sp_peaks_did_not_run`: the design leaves its first clause on the default ('still raises on the default — the opposite-reason branch, `match="sp_tensor_scale"` holds') and only rewords it. That clause then passes through polrad-full's double-count branch, so it no longer tests what its name says. Pin t-peak in that clause.
8. Small quotation and line errors. rc.hpp '1224-1227' is lines 1225-1228 (1224 is blank). test_rc_pipeline.cpp:247 contains no 'back to the default' text. The 'stale test_rc.cpp:374' citations are the range `374–397`, not a point citation.

### Underrated risks (6)

1. Silent test-coverage loss. The 57 skipped knob-matrix cells leave pytest green, and the design's gates check only '0 failed' (steps 2–3) and the failure list. Add a gate that compares pass and skip counts per file against the pre-I5 baseline; otherwise this regression cannot be seen.
2. Narrow long double platforms. Besides the tagged/coherent deadlock above: on aarch64 Linux, long double is software binary128, and the Eq. (18) θ evaluation is long double throughout (67 `ld` uses in rc.cpp). The default build could take minutes per Pipeline there. This is inferred and unmeasured. CI is x86_64 only, so CI will not show it.
3. Reproduction recipe across D4/D5/D7. Every Born-dependent RC number published before 2026-10-27 needs `--rc-tail-model t-peak --unpol-sf toy --pol-sf toy` to reproduce. Event-level numbers cannot be reproduced after D4 at all, and HepMC3 files carry no rc mode, tail model or backend (R7). Documents that say 'reproduce with t-peak' will mislead.
4. Published 'whole tail' orderings: ×768/×681/×289, −1.769797e−07, ×858/×1086, the qe_suppression/qe_kf_gev/c0_shape ΔA_zz spreads. These stay as unlabelled t-peak numbers in rc.hpp, PHYSICS, USAGE and OPEN_ITEMS, and the registry note proposed for row 17 is only optional. Nobody regenerates them.
5. T6 runtime pin. I measured 19.1 s unloaded (the design: 21.1 s loaded) against `< 30 s`, and T6 runs in the default pytest of the main CI job on 4-core GitHub runners whose speed was never measured (ci.yml: 'EXTRAPOLATIONS'). That is about 30 % margin with no measurement of the runner.
6. The default path now carries the (x, y) table's +1.9…+2.5 % bias, and a 1+tail_max = 11 cap on six x = 0.955 cells where the true w_tail is about 4500 (USAGE:2564-2572). I measured 0 clipped events out of 50 000 on all three 6Li configs, but only in the default acceptance window. Any wider-x or other-acceptance study will meet the cap on the default path.

### The checker's notes: what it ran and confirmed, and what it did not

I edited no repository file. All checks were run on ca9b23a with `source env.sh`. The machine was nearly idle (load 0.0–0.6), so my timings are unloaded and run about 15–20 % below the design's loaded ones.

Measurements:
- Pipeline construction minus rc off, 6Li config 1, CLI grid: t-peak 0.157 s, polrad-full 12.35 s.
- On the 40×24 grid: 0.03 s and 1.81 s.
- On configs 0/1/2 (polrad-full): 14.06 / 12.56 / 12.31 s.
- T6 in process (rc file argv + --rc-tail-model): today 5.23 s; t-peak 6.07 s; polrad-full 19.11 s, of which the rc file is 14.08 s. 23/23 clauses pass on both.
- The design's §2.5 per-point numbers reproduce digit for digit on both models.
- T8(0) emulation: σ^q_U at the x = 1 edge is 2.287e-33 (t-peak) and 3.124e+08 (polrad-full). The §2.8 finding stands.
- Default node clip: 0.3986 / 0.3858 / 0.3600 % on configs 0/1/2, 0 of 50 000 tail events clipped on each.
- Gates: strict counts reproduce. `--audit-ranges`: 209 citations, 0 refused. Records: now 1389 citations in 74 records, 0 broken, because other workers have added records since the design ran.
- R6 exposure of the D5-rewritten lines: 16 checked record citations, mostly in the D5 plan itself, all with robust tokens. Low risk.

Mirror method: today's default is TPeak, so a non-default 'polrad-full' base on a tagged or coherent channel behaves exactly as an explicit t-peak will after D5. I used this for the knob-matrix and A5 findings.

Not done: no full pytest or C++ suite run, no rebuild, no D5 flip.

Scratch scripts are in /tmp/claude-1000/-home-cpeng-Projects-polli-PoliGen/7790c59b-ea07-4677-a36f-0183e2410ef4/scratchpad/d5chk/:
- mirror.py
- build_time.py
- t6_time.py (its T6 HepMC3 output is in t6out/)
- t8zero.py
- points.py
- clip_configs.py
- knob_split.py, knob_split2.py, knob_base.py
- cites.py, r6_exposure.py, allowkeys.py

Suggested fixes:
- Make RC_TPK and the low-x `_pipeline` change apply only when channel == 'inclusive'.
- Add a pass/skip-count gate to steps 2–3.
- Extend the §3 inventory with the missing sites listed.
- Change the reproduction sentences to name the backend flags as well as the tail model.
- Add the tagged/coherent deadlock to R2. A real fix would make the PolradFull preconditions apply only when the tail applies, but that is a code change for the author to decide, outside D5.
